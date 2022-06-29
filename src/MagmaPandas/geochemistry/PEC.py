import pandas as pd
import numpy as np
import elements as e
from typing import Union
from MagmaPandas.geochemistry.Kd_ol_melt import Kd_FeMg, Kd_FeMg_vectorised
from MagmaPandas.geochemistry.Fe_redox import Fe_redox
from MagmaPandas.geochemistry.fO2 import fO2_QFM
from MagmaPandas import MagmaSeries, Olivine, Melt_inclusion, Melt
from MagmaPandas.configuration import configuration


class PEC_configuration:

    # Initial stepsize
    stepsize_equilibration = 0.002
    stepsize_crystallisation = 0.005
    # Reduction factor for stepsize after overstepping
    decrease_factor = 5
    # Convergence values
    FeO_converge = 0.05
    Kd_converge = 0.001
    temperature_converge = 0.1

    @staticmethod
    def reset():
        PEC_configuration.stepsize_equilibration = 0.002
        PEC_configuration.stepsize_crystallisation = 0.005
        PEC_configuration.decrease_factor = 5
        PEC_configuration.FeO_converge = 0.05
        PEC_configuration.Kd_converge = 0.005
        PEC_configuration.temperature_converge = 0.1

    @classmethod
    def print(cls):
        """ 
        
        """

        variables = {
            "Stepsize equilibration (moles)": "stepsize_equilibration",
            "Stepsize crystallisation (moles)": "stepsize_crystallisation",
            "Decrease factor": "decrease_factor",
            "FeO convergence (wt. %)": "FeO_converge",
            "Kd convergence": "Kd_converge",
            "Temperature convergence (\N{DEGREE SIGN})": "temperature_converge",
        }

        names_length = max([len(i) for i in variables.keys()]) + 5
        values_length = max([len(str(getattr(cls, i))) for i in variables.values()])

        pad_right = 20
        pad_total = names_length + pad_right
        print(" Post-entrapment modification ".center(pad_total, "#"))
        print(" correction model ".center(pad_total, "#"))
        print("\nSettings".ljust(pad_total, "_"))
        for param, value in variables.items():
            value_str = f"{getattr(cls, value): <{values_length}}"
            print(f"{param:.<{names_length}}{value_str:.>{pad_right}}")


def calculate_olivine(forsterite):
    """
    Calculate olivine wt. % composition from forsterite content
    """
    if isinstance(forsterite, pd.Series):
        idx = forsterite.index
    else:
        try:
            idx = np.arange(0, len(forsterite))
        except:
            idx = [0]

    MgO = forsterite * 2
    FeO = (1 - forsterite) * 2
    SiO2 = 1
    total = MgO + FeO + SiO2

    cations = Olivine(
        {"MgO": MgO, "FeO": FeO, "SiO2": SiO2, "total": total},
        index=idx,
        units="mol fraction",
        datatype="oxide",
    )
    cations = cations.convert_moles_wtPercent

    return cations


def Fe_restore(
    inclusion: Union[Melt_inclusion, Melt],
    FeO_target: Union[int, float],
    P_bar: float,
    **kwargs,
):
    """
    Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
    Isothermal and isobaric.
    """
    # Grab model parameters
    converge = kwargs.get(
        "converge", getattr(PEC_configuration, "FeO_converge")
    )  # In wt. %
    temperature_converge = kwargs.get(
        "temperature_converge", getattr(PEC_configuration, "temperature_converge")
    )  # In degrees
    stepsize = kwargs.get(
        "stepsize", getattr(PEC_configuration, "stepsize")
    )  # In molar fraction
    # Select Fe loss or gain
    if inclusion["FeO"] > FeO_target:
        stepsize = -stepsize
    forsterite = kwargs.get("forsterite", 0.8)
    QFM_logshift = kwargs.get("QFM_logshift", configuration().QFMlogshift)
    # Parameters for the while loop
    olivine_corrected = 0
    olivine_stepsize_reduction = 4
    decrease_factor = getattr(PEC_configuration, "decrease_factor")
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements]
    inclusion = inclusion.normalise()
    # Calculate temperature and fO2
    temperature = inclusion.melt_temperature(P_bar=P_bar)
    fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)

    # Set up initial data
    moles = Melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    moles.loc[0] = inclusion.moles[inclusion.elements].values
    # Fe-Mg exchange vector
    FeMg_exchange = pd.Series(0, index=moles.columns)
    FeMg_exchange.loc[["FeO", "MgO"]] = stepsize, -stepsize
    # Inclusion starting FeO
    FeO = inclusion["FeO"]

    while not np.isclose(FeO, FeO_target, atol=converge):
        # Exchange Fe and Mg
        idx = moles.index[-1] + stepsize
        moles.loc[idx] = (moles.iloc[-1] + FeMg_exchange).values
        # New model temperature after Fe-Mg exchange
        temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )

        # Calculate equilibrium olivine Fo#
        # melt Fe3+/Fe2+
        Fe3Fe2 = Fe3Fe2_model(moles.iloc[-1], temperature, fO2)
        # FeMg Kd
        Kd = Kd_model(moles.iloc[-1], forsterite, temperature, P_bar, Fe3Fe2)
        Fe2_FeTotal = 1 / (1 + Fe3Fe2)
        Fe2Mg = moles.loc[idx, "FeO"] * Fe2_FeTotal / moles.loc[idx, "MgO"]
        # Equilibrium forsterite content
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
        # Equilibrium olivine composition in oxide mol fraction
        olivine = MagmaSeries(
            {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1}
        )
        olivine = olivine.normalise().reindex(moles.columns, fill_value=0.0)

        # Crystallise or melt olivine to remain isothermal
        add_olivine = moles.iloc[-1]
        # Set stepsize for olivine addition/removal
        olivine_stepsize = stepsize / olivine_stepsize_reduction
        while not np.isclose(temperature_new, temperature, atol=temperature_converge):
            # crystallise olivine until calculated temperature is back at initial
            add_olivine = add_olivine + olivine * olivine_stepsize
            temperature_new = add_olivine.convert_moles_wtPercent.melt_temperature(
                P_bar=P_bar
            )
            # Record added/removed olivine amount. Negative values for crystallisation
            olivine_corrected += olivine_stepsize
            T_overstepped = np.sign(temperature - temperature_new) != np.sign(stepsize)
            # Reverse one iteration and reduce stepsize if temperature was
            # overstepped by more than the convergence value
            if T_overstepped and not np.isclose(
                temperature_new, temperature, atol=temperature_converge
            ):
                add_olivine = add_olivine - olivine * olivine_stepsize
                olivine_corrected -= olivine_stepsize
                olivine_stepsize = olivine_stepsize / decrease_factor
                continue

        # Copy olivine corrected composition
        moles.iloc[-1] = add_olivine.values

        # New inclusion FeO
        FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

        overstepped = np.sign(FeO_target - FeO) != np.sign(stepsize)
        # Reverse one iteration and reduce stepsize if FeO content
        # gets oversteppend by more than the convergence value
        if overstepped and not np.isclose(FeO, FeO_target, atol=converge):
            moles.drop([idx], inplace=True)
            FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]
            stepsize = stepsize / decrease_factor
            FeMg_exchange = FeMg_exchange.div(decrease_factor)

    # Recalculate compositions to oxide wt. %
    wtPercent = moles.convert_moles_wtPercent

    return wtPercent, temperature, temperature_new, olivine_corrected


def Fe_equilibrate(
    inclusion: Union[Melt_inclusion, Melt],
    olivine: Union[float, MagmaSeries],
    P_bar: float,
    **kwargs,
):
    """
    Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
    Isothermal and isobaric.
    """
    # Grab model parameters
    stepsize = kwargs.get(
        "stepsize", getattr(PEC_configuration, "stepsize_equilibration")
    )  # In molar fraction, this is the maximum recommende stepsize.
    converge = kwargs.get(
        "converge", getattr(PEC_configuration, "Kd_converge")
    )  # Kd converge
    temperature_converge = kwargs.get(
        "temperature_converge", getattr(PEC_configuration, "temperature_converge")
    )  # In degrees
    QFM_logshift = kwargs.get("QFM_logshift", configuration().QFMlogshift)
    # Parameters for the while loop
    olivine_crystallised = 0
    olivine_stepsize_reduction = 4
    decrease_factor = getattr(PEC_configuration, "decrease_factor")
    # Normalise inclusion composition
    inclusion = inclusion.fillna(0.0)
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.normalise()
    # Calculate temperature and fO2
    temperature = inclusion.melt_temperature(P_bar=P_bar)
    fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

    # Get olivine molar oxide fractions
    if isinstance(olivine, float):
        if olivine < 0 or olivine > 1:
            raise ValueError(
                f"olivine host forsterite: {olivine:.3f} number is not between 0 and 1"
            )
        forsterite = olivine
    elif not isinstance(olivine, MagmaSeries):
        raise TypeError(
            f"Olivine host should be forsterite number as float, or full composition as MagmaSeries, not {type(olivine)}"
        )
    else:
        olivine = olivine.moles
        forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])

    # Function for calculating equilibrium forsterite content
    def calculate_Kd(
        mol_fractions,
        T_K=temperature,
        P_bar=P_bar,
        oxygen_fugacity=fO2,
        Fo_initial=forsterite,
    ):

        # melt Fe3+/Fe2+
        Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
        # FeMg Kd
        Kd = Kd_model(mol_fractions, Fo_initial, T_K, P_bar, Fe3Fe2)
        Fe2_FeTotal = 1 / (1 + Fe3Fe2)

        return Kd, Fe2_FeTotal

    # Set up initial data
    mi_moles = Melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    mi_moles.loc[0] = inclusion.moles[inclusion.elements].values
    mi_moles = mi_moles.normalise()
    # Equilibrium Kd
    Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.loc[0])
    # Real Kd
    olivine_MgFe = forsterite / (1 - forsterite)
    melt_MgFe = mi_moles.loc[0, "MgO"] / (mi_moles.loc[0, "FeO"] * Fe2_FeTotal)
    Kd_real = melt_MgFe / olivine_MgFe
    # Fe-Mg exchange vector
    FeMg_vector = pd.Series(0, index=mi_moles.columns)
    FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
    # Select Fe removal or addition
    if Kd_real < Kd_equilibrium:
        stepsize = -stepsize

    ##### MAIN LOOP #####
    #####################
    while not np.isclose(Kd_real, Kd_equilibrium, atol=converge):

        # Exchange Fe-Mg
        idx = mi_moles.index[-1] + stepsize
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + FeMg_vector.mul(stepsize)).values
        mi_moles = mi_moles.normalise()
        # New model temperature after Fe-Mg exchange
        temperature_new = mi_moles.loc[idx].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )

        # Equilibrium Kd and Fe speciation for new composition
        Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.loc[idx])
        melt_FeMg = (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal) / mi_moles.loc[idx, "MgO"]
        # Equilibrium olivine composition in oxide mol fractions
        Fo_equilibrium = 1 / (1 + Kd_equilibrium * melt_FeMg)
        olivine = MagmaSeries(
            {"MgO": Fo_equilibrium * 2, "FeO": (1 - Fo_equilibrium) * 2, "SiO2": 1}
        )
        olivine = olivine.normalise().reindex(mi_moles.columns, fill_value=0.0)

        ##### OLIVINE CRYSTALLISATION LOOP #####
        ########################################
        remove_olivine = mi_moles.loc[idx]
        # Set stepsize
        olivine_stepsize = stepsize / olivine_stepsize_reduction
        while not np.isclose(temperature_new, temperature, atol=temperature_converge):
            # Melt or crystallise olivine until the temperature is back to original.
            remove_olivine = remove_olivine + olivine * olivine_stepsize
            temperature_new = remove_olivine.convert_moles_wtPercent.melt_temperature(
                P_bar=P_bar
            )
            # Record removed/added olivine amount. Negative values for crystallisation
            olivine_crystallised += olivine_stepsize

            T_overstepped = np.sign(temperature - temperature_new) != np.sign(stepsize)
            # Reverse one iteration and reduce stepsize if temperature is
            # overstepped by more than the convergence value
            Temperature_mismatch = ~np.isclose(
                temperature_new, temperature, atol=temperature_converge
            )
            decrease_stepsize_T = np.logical_and(T_overstepped, Temperature_mismatch)
            if decrease_stepsize_T:
                remove_olivine = remove_olivine - olivine * olivine_stepsize
                olivine_crystallised -= olivine_stepsize
                olivine_stepsize = olivine_stepsize / decrease_factor

        # Copy olivine corrected composition
        mi_moles.loc[idx] = remove_olivine.values
        mi_moles = mi_moles.normalise()

        # New equilibrium Kd and Fe speciation
        Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.loc[idx])
        # Real Kd
        melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe

        disequilibrium = ~np.isclose(Kd_equilibrium, Kd_real, atol=converge)
        overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
        decrease_stepsize = np.logical_and(disequilibrium, overstepped)
        # Reverse one iteration and reduce stepsize if Kd
        # gets oversteppend by more than the convergence value
        if decrease_stepsize:
            mi_moles.drop(index=idx, inplace=True)
            # Reset equilibrium and real Kd
            Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.iloc[-1])
            idx = mi_moles.index[-1]
            melt_MgFe = mi_moles.loc[idx, "MgO"] / (
                mi_moles.loc[idx, "FeO"] * Fe2_FeTotal
            )
            Kd_real = melt_MgFe / olivine_MgFe
            stepsize = stepsize / decrease_factor

    # Recalculate compositions to oxide wt. %
    equilibrated_composition = mi_moles.convert_moles_wtPercent

    return (
        equilibrated_composition,
        olivine_crystallised,
        {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
    )


def crystallisation_correction(
    inclusion: Union[Melt_inclusion, Melt],
    olivine_host: Union[float, MagmaSeries],
    FeO_target: Union[int, float],
    P_bar: float,
    **kwargs,
):

    """
    Correct an olivine hosted melt inclusion for post entrapment crystallisation or melting by
    respectively melting or crystallising host olivine.
    Expects the melt inclusion is completely equilibrated with the host crystal.
    The models exits when the user input original melt inclusion FeO content is reached.
    Loosely based on the postentrapment reequilibration procedure in Petrolog:

    L. V. Danyushesky and P. Plechov (2011)
    Petrolog3: Integrated software for modeling crystallization processes
    Geochemistry, Geophysics, Geosystems, vol 12
    """
    # Grab model parameters
    stepsize = kwargs.get(
        "stepsize", getattr(PEC_configuration, "stepsize_crystallisation")
    )  # In molar fraction, 0.005 is the maximum recommended value
    converge = kwargs.get(
        "converge", getattr(PEC_configuration, "FeO_converge")
    )  # FeO convergence
    Kd_converge = kwargs.get(
        "Kd_converge", getattr(PEC_configuration, "Kd_converge")
    )  # Kd converge
    QFM_logshift = kwargs.get("QFM_logshift", configuration().QFMlogshift)
    # Parameters for the while loop
    olivine_corrected = 0
    FeMg_exchange_reduction = 4
    decrease_factor = 5
    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.normalise()
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

    # SET UP INITIAL DATA
    # Dataframe with new compositions ()
    mi_moles = Melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    mi_moles.loc[0] = inclusion.moles[inclusion.elements].values
    mi_moles = mi_moles.normalise()
    # Get olivine molar oxide fractions
    if isinstance(olivine_host, float):
        if olivine_host < 0 or olivine_host > 1:
            raise ValueError(
                f"olivine host forsterite: {olivine_host:.3f} number is not between 0 and 1"
            )
        olivine = MagmaSeries(
            {"MgO": olivine_host * 2, "FeO": (1 - olivine_host) * 2, "SiO2": 1},
            units="mol fraction",
            datatype="oxide",
        )
        olivine = olivine.normalise().reindex(mi_moles.columns, fill_value=0.0)
        forsterite = olivine_host
    elif not isinstance(olivine_host, MagmaSeries):
        raise TypeError(
            f"Olivine host should be forsterite number as float, or full composition as MagmaSeries, not {type(olivine_host)}"
        )
    else:
        olivine = olivine_host.moles
        olivine = olivine.reindex(mi_moles.columns, fill_value=0.0)
        olivine.recalculate()
        forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])

    # Fe-Mg exchange vector
    FeMg_vector = pd.Series(0, index=mi_moles.columns)
    FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
    # Inclusion starting FeO
    FeO = inclusion["FeO"]
    temperature_old = mi_moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
        P_bar=P_bar
    )

    def calculate_Kd(melt, pressure=P_bar, fObuffer_shift=QFM_logshift, Fo=forsterite):
        T_K = melt.convert_moles_wtPercent.melt_temperature(pressure)
        fO2 = fO2_QFM(fObuffer_shift, T_K, pressure)
        Fe3Fe2 = Fe3Fe2_model(melt, T_K, fO2)

        Fe2_FeTotal = 1 / (1 + Fe3Fe2)
        melt_MgFe = melt["MgO"] / (melt["FeO"] * Fe2_FeTotal)
        olivine_MgFe = Fo / (1 - Fo)
        Kd_observed = melt_MgFe / olivine_MgFe

        Kd_eq =  Kd_model(melt, Fo, T_K, P_bar, Fe3Fe2)

        return Kd_eq, Kd_observed        


    # Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.iloc[-1])

    if FeO > FeO_target:
        stepsize = -stepsize

    ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
    while not np.isclose(FeO, FeO_target, atol=converge):

        idx = round(mi_moles.index[-1] + stepsize, 4)
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + olivine.mul(stepsize)).values
        mi_moles = mi_moles.normalise()

        olivine_corrected += stepsize

        Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.loc[idx])

        ###### FE-MG EXCHANGE LOOP #####
        stepsize_FeMg = stepsize / FeMg_exchange_reduction
        FeMg_exchange = mi_moles.loc[idx]
        while not np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge):

            FeMg_exchange += FeMg_vector.mul(stepsize_FeMg)
            FeMg_exchange = FeMg_exchange.normalise()
            Kd_equilibrium, Kd_real = calculate_Kd(FeMg_exchange)

            Kd_overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
            Kd_mismatch = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge)
            decrease_stepsize_FeMg = np.logical_and(Kd_overstepped, Kd_mismatch)
            # Reverse one iteration and reduce stepsize if Kd
            # gets oversteppend by more than the convergence value
            if decrease_stepsize_FeMg:
                FeMg_exchange -= FeMg_vector.mul(stepsize_FeMg)
                stepsize_FeMg = stepsize_FeMg / decrease_factor

        mi_moles.loc[idx] = FeMg_exchange
        mi_moles = mi_moles.normalise()

        mi_wtPercent = mi_moles.convert_moles_wtPercent
        FeO = mi_wtPercent.loc[idx, "FeO"]

        disequilibrium = ~np.isclose(FeO, FeO_target, atol=converge)
        overstepped = np.sign(FeO_target - FeO) != np.sign(stepsize)
        decrease_stepsize = np.logical_and(disequilibrium, overstepped)
        # Reverse one iteration and reduce stepsize if FeO
        # gets oversteppend by more than the convergence value
        if decrease_stepsize:
            mi_moles.drop(index=idx, inplace=True)
            mi_wtPercent = mi_moles.convert_moles_wtPercent
            idx = mi_wtPercent.index[-1]
            FeO = mi_wtPercent.loc[idx, "FeO"]
            stepsize = stepsize / decrease_factor

    temperature_new = mi_wtPercent.loc[idx].melt_temperature(P_bar)

    return (
        mi_wtPercent,
        olivine_corrected,
        {"equilibrium": Kd_equilibrium, "real": Kd_real},
        {"old": temperature_old, "new": temperature_new},
    )


class PEC_olivine:
    def __init__(
        self,
        inclusions: Melt_inclusion,
        olivines: Union[Olivine, pd.Series, float],
        P_bar: Union[float, int, pd.Series],
        FeO_target: Union[float, pd.Series, callable],
        **kwargs,
    ):

        # Initialise variables
        self.olivine_corrected = pd.Series(0, index=inclusions.index)
        self._FeO_as_function = False

        # Process attributes
        ######################
        # For inclusions
        if not isinstance(inclusions, Melt):
            raise TypeError("Inclusions is not a Melt MagmaFrame")
        else:
            inclusions = inclusions.fillna(0.0)
            self.inclusions = inclusions.normalise()
        # For olivine
        if hasattr(olivines, "index"):
            if not olivines.index.equals(self.inclusions.index):
                raise ValueError("Inclusions and olivines indeces don't match")
        if not isinstance(olivines, Olivine):
            try:
                if len(olivines) != self.shape[0]:
                    raise ValueError("Number of olivines and inclusions does not match")
            except TypeError:
                pass
            forsterite = pd.Series(olivines, index=self.inclusions.index)
            if (~forsterite.between(0, 1)).any():
                raise ValueError(
                    "olivine host forsterite contents are not all between 0 and 1"
                )
            olivine = Olivine(
                {"MgO": forsterite * 2, "FeO": (1 - forsterite) * 2, "SiO2": 1},
                index=self.inclusions.index,
                units="mol fraction",
                datatype="oxide",
            )
            self.olivine = olivine.normalise()
        else:
            self.olivine = olivines.moles
        self.olivine = self.olivine.reindex(columns=self.inclusions.columns, fill_value=0.0)
        # For pressure
        try:
            if len(P_bar) != self.inclusions.shape[0]:
                raise ValueError(
                    "Number of pressure inputs and melt inclusions does not match"
                )
        except TypeError:
            pass
        if hasattr(P_bar, "index"):
            if not P_bar.index.equals(self.inclusions.index):
                raise ValueError("Pressure inputs and inclusion indeces don't match")
            self.P_bar = P_bar
        else:
            self.P_bar = pd.Series(P_bar, index=self.inclusions.index)
        # For FeO
        try:
            if len(FeO_target) != len(self.inclusions):
                raise ValueError(
                    "Number of initial FeO inputs and inclusions does not match"
                )
            if hasattr(FeO_target, "index"):
                if not FeO_target.equals(self.inclusions.index):
                    raise ValueError(
                        "FeO target inputs and inclusion indeces don't match"
                    )
                self.FeO_target = FeO_target
            else:
                self.FeO_target = pd.Series(FeO_target, index=self.inclusions.index)
        except TypeError:
            if isinstance(FeO_target, (int, float)):
                self.FeO_target = pd.Series(FeO_target, index=self.inclusions.index)
            elif hasattr(FeO_target, "__call__"):
                self._FeO_as_function = True
                self.FeO_function = FeO_target
                self.FeO_target = self.FeO_function(self.inclusions)

    @property
    def equilibrated(self):
        """
        Check which inclusions are in equilibrium with their host olivine
        """
        Kd_converge = getattr(PEC_configuration, "Kd_converge")
        Kd_equilibrium, Kd_real = self.calculate_Kds()

        return pd.Series(np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge), index=Kd_equilibrium.index, name="equilibrated")


    @property
    def Fe_loss(self):
        """
        Check which inclusions have experienced Fe loss
        """
        FeO_converge = getattr(PEC_configuration, "FeO_converge")
        return pd.Series(~np.isclose(
                self.FeO_target,
                self.inclusions.convert_moles_wtPercent["FeO"],
                FeO_converge,
            ), index = self.inclusions.index)

    def calculate_Kds(self):

            Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
            Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)
            dQFM = configuration().QFMlogshift
            moles = self.inclusions.moles
            # Calculate temperature and fO2
            temperatures = self.inclusions.temperature(P_bar=self.P_bar)
            fO2 = fO2_QFM(dQFM, temperatures, self.P_bar)
            forsterite = self.olivine.forsterite
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(moles, temperatures, fO2)
            # Real Kd
            olivine_MgFe = forsterite / (1 - forsterite)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            melt_MgFe = moles["MgO"] / (moles["FeO"] * Fe2_FeTotal)
            Kd_real = melt_MgFe / olivine_MgFe
            Kd_real.rename("real", inplace=True)
            # Equilibrium Kd
            Kd_equilibrium = Kd_model(moles, forsterite, temperatures, self.P_bar, Fe3Fe2)
            Kd_equilibrium.rename("equilibrium", inplace=True)
            return Kd_equilibrium, Kd_real

    def Fe_equilibrate(self, inplace=False):

        # Get settings
        stepsize = getattr(PEC_configuration, "stepsize_equilibration")
        stepsize = pd.Series(stepsize, index=self.inclusions.index)
        decrease_factor = getattr(PEC_configuration, "decrease_factor")
        olivine_stepsize_reduction = 4
        Kd_converge = getattr(PEC_configuration, "Kd_converge")
        temperature_converge = getattr(PEC_configuration, "temperature_converge")
        dQFM = configuration().QFMlogshift
        P_bar = self.P_bar
        # Get configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)
        # Calculate temperature and fO2
        temperatures = self.inclusions.temperature(P_bar=P_bar)
        fO2 = fO2_QFM(dQFM, temperatures, P_bar)
        # Get olivine forsterite contents
        forsterite_host = self.olivine.forsterite

        # Function for calculating equilibrium forsterite content and Fe2+/Fe(total)
        def calculate_Kd(
            mol_fractions,
            T_K,
            P_bar,
            oxygen_fugacity,
            forsterite_host,
        ):
            forsterite = forsterite_host.copy()
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
            # FeMg Kd
            Kd = Kd_model(mol_fractions, forsterite, T_K, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)

            return Kd, Fe2_FeTotal

        # Set up initial data
        mi_moles = self.inclusions.moles
        # Calculate equilibrium Kd
        Kd_equilibrium, Fe2_FeTotal = calculate_Kd(
            mi_moles, temperatures, P_bar, fO2, forsterite_host
        )
        # Real Kd
        olivine_MgFe = forsterite_host / (1 - forsterite_host)
        melt_MgFe = mi_moles["MgO"] / (mi_moles["FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe
        # Fe-Mg exchange vectors
        FeMg_exchange = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
        FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]
        # Find disequilibrium inclusions
        disequilibrium = ~np.isclose(Kd_real, Kd_equilibrium, atol=Kd_converge)
        # Set stepsizes acoording to Kd disequilibrium
        stepsize.loc[Kd_real < Kd_equilibrium] = -stepsize.loc[Kd_real < Kd_equilibrium].copy()

        ##### Main Fe-Mg exhange loop #####
        ###################################
        # Only reslice al the dataframes and series if the amount of disequilibrium inclusions has changed
        reslice = True
        while sum(disequilibrium) > 0:
            if reslice:
                stepsize_loop = stepsize.loc[disequilibrium]
                temperatures_loop = temperatures.loc[disequilibrium]
                P_loop = P_bar.loc[disequilibrium]
                FeMg_loop = FeMg_exchange.loc[disequilibrium, :]
                mi_moles_loop = mi_moles.loc[disequilibrium, :].copy()
                fO2_loop = fO2.loc[disequilibrium]
                Fo_host_loop = forsterite_host.loc[disequilibrium]
                Kd_eq_loop = Kd_equilibrium.loc[disequilibrium].copy()
                Kd_real_loop = Kd_real.loc[disequilibrium].copy()

            # Exchange Fe and Mg
            mi_moles_loop = mi_moles_loop + FeMg_loop.mul(stepsize_loop, axis=0)
            mi_moles_loop = mi_moles_loop.normalise()
            # Calculate new liquidus temperature
            temperatures_new = mi_moles_loop.convert_moles_wtPercent.temperature(
                P_bar=P_bar
            )

            # Calculate new equilibrium Kd and Fe speciation
            Kd_eq_loop, Fe2_FeTotal = calculate_Kd(
                mi_moles_loop, temperatures_loop, P_loop, fO2_loop, Fo_host_loop
            )
            melt_FeMg = (mi_moles_loop["FeO"] * Fe2_FeTotal) / mi_moles_loop["MgO"]
            # equilibrium olivine composition in oxide mol fractions
            Fo_equilibrium = 1 / (1 + Kd_eq_loop * melt_FeMg)
            olivine = Olivine(
                {
                    "MgO": Fo_equilibrium * 2,
                    "FeO": (1 - Fo_equilibrium) * 2,
                    "SiO2": 1,
                },
                index=mi_moles_loop.index,
                units="mol fraction",
                datatype="oxide",
            )
            olivine = olivine.normalise().reindex(
                columns=mi_moles_loop.columns, fill_value=0.0
            )

            # Find inclusions outside the equilibrium temperature range
            temperature_mismatch = ~np.isclose(
                temperatures_new, temperatures_loop, atol=temperature_converge
            )

            ############################# Olivine melting/crystallisation loop ############################
            # Melt or crystallise olivine until the liquidus temperature is back in the equilibrium range #
            # Initialise data for the olivine melting/crystallisation loop
            olivine_stepsize = stepsize_loop.div(olivine_stepsize_reduction)
            # olivine_correction = mi_moles_loop.loc[temperature_mismatch, :].copy()
            reslice_olivine = True
            T = 1
            while sum(temperature_mismatch) > 0:
                print(f"\rT correction {T:02}", end="")
                # Gather all loop data
                if reslice_olivine:
                    olivine_correction = mi_moles_loop.loc[temperature_mismatch, :]
                    olivine_loop = olivine.loc[temperature_mismatch, :]
                    ol_stepsize_loop = olivine_stepsize.loc[temperature_mismatch]
                    idx_olivine = olivine_correction.index

                # Crystallise or melt olivine
                mi_moles_loop.loc[
                    temperature_mismatch, :
                ] = olivine_correction + olivine_loop.mul(ol_stepsize_loop, axis=0)
                # Keep track of crystllised or melted olivine; negative values for crystallisation
                self.olivine_corrected.loc[idx_olivine] += olivine_stepsize
                # Calculate the new liquidus temperature
                temperatures_new = mi_moles_loop.convert_moles_wtPercent.temperature(
                    P_bar=P_loop
                )
                # Find inclusions outside the equilibrium temperature range
                temperature_mismatch_new = ~np.isclose(
                    temperatures_new,
                    temperatures_loop,
                    atol=temperature_converge,
                )
                # Find overcorrected inclusions
                T_overstepped = ~np.equal(
                    np.sign(temperatures_loop - temperatures_new),
                    np.sign(olivine_stepsize),
                )
                # Reverse one iteration and reduce stepsize for overcorrected inclusions
                decrease_stepsize_T = np.logical_and(
                    T_overstepped, temperature_mismatch_new
                )
                if sum(decrease_stepsize_T) > 0:
                    idx_reverse = mi_moles_loop.index[decrease_stepsize_T]
                    # Reverse meltind/crystallisation
                    reverse_olivine = olivine_loop.loc[idx_reverse, :].mul(
                        ol_stepsize_loop.loc[idx_reverse], axis=0
                    )
                    mi_moles_loop.loc[idx_reverse, :] = (
                        mi_moles_loop.loc[idx_reverse, :] - reverse_olivine
                    )
                    self.olivine_corrected.loc[idx_reverse] -= olivine_stepsize.loc[
                        idx_reverse
                    ]
                    # Decrease stepsize
                    olivine_stepsize.loc[idx_reverse] = olivine_stepsize.loc[
                        idx_reverse
                    ].div(decrease_factor)

                reslice_olivine = ~np.array_equal(
                    temperature_mismatch_new, temperature_mismatch
                )
                temperature_mismatch = temperature_mismatch_new
                T += 1

            # Renormalise after olivine correction
            mi_moles_loop = mi_moles_loop.normalise()

            # Recalculate equilibrium Kd and Fe speciation
            Kd_eq_loop, Fe2_FeTotal = calculate_Kd(
                mi_moles_loop, temperatures_loop, P_loop, fO2_loop, Fo_host_loop
            )
            # Real Kd
            melt_MgFe = mi_moles_loop["MgO"] / (mi_moles_loop["FeO"] * (Fe2_FeTotal))
            Kd_real_loop = melt_MgFe / olivine_MgFe.loc[disequilibrium]

            # Find inclusions outside the equilibrium forsterite range
            disequilibrium_loop = ~np.isclose(
                Kd_eq_loop, Kd_real_loop, atol=Kd_converge
            )
            # Find overcorrected inclusions
            overstepped = ~np.equal(
                np.sign(Kd_real_loop - Kd_eq_loop), np.sign(stepsize_loop)
            )
            # Reverse one iteration and decrease stepsize for overcorrected inclusions
            decrease_stepsize = np.logical_and(overstepped, disequilibrium_loop)
            if sum(decrease_stepsize) > 0:
                idx_FeMg = mi_moles_loop.index[decrease_stepsize]
                mi_moles_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                Kd_eq_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                Kd_real_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                stepsize.loc[idx_FeMg] = stepsize_loop.loc[idx_FeMg].div(
                    decrease_factor
                )

            mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
            Kd_equilibrium[Kd_eq_loop.index] = Kd_eq_loop.copy()
            Kd_real[Kd_real_loop.index] = Kd_real_loop.copy()
            disequilibrium_new = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge)
            reslice = ~np.array_equal(disequilibrium_new, disequilibrium)
            disequilibrium = disequilibrium_new
            # Make sure that a dataframe sliced with disequilibrium still returns a dataframe
            if not isinstance(disequilibrium, np.ndarray):
                disequilibrium = np.array([disequilibrium])

        corrected_compositions = mi_moles.convert_moles_wtPercent

        if inplace:
            # Convert corrected compositions back to oxide wt. %
            self.inclusions = corrected_compositions
            self._equilibrated = True
        else:
            return (
                corrected_compositions,
                self.olivine_corrected,
                {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
            )

    def correct_olivine(self, inplace=False, *args, **kwargs):
        """
        Correct an olivine hosted melt inclusion for post entrapment crystallisation or melting by
        respectively melting or crystallising host olivine.
        Expects the melt inclusion is completely equilibrated with the host crystal.
        The models exits when the user input original melt inclusion FeO content is reached.
        Loosely based on the postentrapment reequilibration procedure in Petrolog:

        L. V. Danyushesky and P. Plechov (2011)
        Petrolog3: Integrated software for modeling crystallization processes
        Geochemistry, Geophysics, Geosystems, vol 12
        """

        if not self._equilibrated:
            raise RuntimeError("Inclusion compositions have not yet been equilibrated")

        # Get settings
        stepsize = getattr(PEC_configuration, "stepsize_crystallisation")
        stepsize = pd.Series(stepsize, index=self.inclusions.index)
        decrease_factor = getattr(PEC_configuration, "decrease_factor")
        FeMg_exchange_reduction = 4
        Kd_converge = getattr(PEC_configuration, "Kd_converge")
        FeO_converge = kwargs.get("FeO_converge", getattr(PEC_configuration, "FeO_converge"))
        dQFM = configuration().QFMlogshift
        P_bar = self.P_bar
        # Get configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)
        # Inclusion compositions in oxide mol fractions
        mi_moles = self.inclusions.moles
        mi_moles = mi_moles.normalise()
        # Olivine forsterite and Mg/Fe
        forsterite = self.olivine.forsterite
        # olivine_MgFe = self.olivine["MgO"] / self.olivine["FeO"]
        # Fe-Mg exchange vectors
        FeMg_vector = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
        FeMg_vector.loc[:, ["FeO", "MgO"]] = [1, -1]
        # Starting FeO and termperature
        FeO = self.inclusions["FeO"]
        FeO_target = self.FeO_target
        temperature_old = self.inclusions.temperature(P_bar=P_bar)


        def calculate_Kd(melt, pressure=P_bar, fObuffer_shift=dQFM, Fo=forsterite):
            """
            Calculate observed and modelled Kds
            """

            T_K = melt.convert_moles_wtPercent.temperature(pressure)
            fO2 = fO2_QFM(fObuffer_shift, T_K, pressure)
            Fe3Fe2 = Fe3Fe2_model(melt, T_K, fO2)

            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            melt_MgFe = melt["MgO"] / (melt["FeO"] * Fe2_FeTotal)
            olivine_MgFe = Fo / (1 - Fo)
            Kd_observed = melt_MgFe / olivine_MgFe

            Kd_eq = Kd_model(melt, Fo.copy(), T_K, pressure, Fe3Fe2)

            return Kd_eq, Kd_observed
        
        # Kd_equilibrium, Kd_real = calculate_Kd(mi_moles)

        stepsize.loc[FeO > FeO_target] = -stepsize.loc[FeO > FeO_target].copy()

        FeO_mismatch = ~np.isclose(FeO, FeO_target, atol=FeO_converge)

        ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
        reslice = True
        while sum(FeO_mismatch) > 0:

            if reslice:
                idx_olivine = mi_moles_loop.index
                stepsize_loop = stepsize.loc[FeO_mismatch]
                
                mi_moles_loop = mi_moles.loc[FeO_mismatch, :].copy()

                olivine_loop = self.olivine.loc[FeO_mismatch,:]
                Fo_loop = forsterite.loc[FeO_mismatch]

                P_loop = P_bar.loc[FeO_mismatch]
                

            mi_moles_loop = mi_moles_loop + olivine_loop.mul(stepsize_loop, axis=0)
            mi_moles_loop = mi_moles_loop.normalise()

            self.olivine_corrected.loc[idx_olivine] += stepsize_loop

            Kd_eq_loop, Kd_real_loop = calculate_Kd(mi_moles_loop, pressure=P_loop, Fo=Fo_loop)
            Kd_mismatch = ~np.isclose(Kd_eq_loop, Kd_real_loop, atol=Kd_converge)

            ##### FE-MG EXCHANGE LOOP #####
            stepsize_FeMg = stepsize_loop.div(FeMg_exchange_reduction)
            reslice_FeMg = True
            while sum(Kd_mismatch) > 0:
                # collect loop variables
                if reslice_FeMg:
                    FeMg_idx = mi_moles_loop.index[Kd_mismatch]
                    FeMg_exchange = mi_moles_loop.loc[Kd_mismatch, :]
                    FeMg_vector_loop = FeMg_vector.loc[FeMg_idx, :]
                    stepsize_FeMg_loop = stepsize_FeMg.loc[Kd_mismatch]
                
                mi_moles_loop.loc[Kd_mismatch, :] = FeMg_exchange + FeMg_vector_loop.mul(stepsize_FeMg_loop, axis=0)
                mi_moles_loop = mi_moles_loop.normalise()

                Kd_eq_loop, Kd_real_loop = calculate_Kd(mi_moles_loop, pressure=P_loop, Fo=Fo_loop)

                Kd_mismatch_new = ~np.isclose(Kd_eq_loop, Kd_real_loop, atol=Kd_converge)
                Kd_overstepped = ~np.equal(np.sign(Kd_real_loop - Kd_eq_loop), np.sign(stepsize_loop))
                decrease_stepsize_Kd = np.logical_and(Kd_overstepped, Kd_mismatch_new)

                if sum(decrease_stepsize_Kd) > 0:
                    reverse_idx = mi_moles_loop.index[decrease_stepsize_Kd]
                    # Reverse Fe-Mg exchange                     
                    reverse_FeMg = FeMg_vector_loop.loc[reverse_idx, :].mul(stepsize_FeMg_loop.loc[reverse_idx], axis=0)
                    mi_moles_loop.loc[reverse_idx, :] = FeMg_exchange.loc[reverse_idx, :] - reverse_FeMg
                    stepsize_FeMg.loc[reverse_idx] = stepsize_FeMg.loc[reverse_idx].div(decrease_factor)
                
                reslice_FeMg = ~np.array_equal(Kd_mismatch, Kd_mismatch_new)
                Kd_mismatch = Kd_mismatch_new


            # Recalculate FeO
            mi_wtPercent = mi_moles_loop.convert_moles_wtPercent
            FeO = mi_wtPercent["FeO"]
            if self._FeO_as_function:
                FeO_target = self.FeO_function(mi_wtPercent)
                self.FeO_target.loc[FeO_mismatch] = FeO_target
            else:
                FeO_target = self.FeO_target.loc[FeO_mismatch]
            FeO_mismatch_loop = ~np.isclose(FeO, FeO_target, atol=FeO_converge)
            FeO_overstepped = ~np.equal(np.sign(FeO_target - FeO), np.sign(stepsize_loop))
            decrease_stepsize_FeO = np.logical_and(FeO_overstepped, FeO_mismatch_loop)

            if sum(decrease_stepsize_FeO) > 0:
                reverse_FeO = mi_moles_loop.index[decrease_stepsize_FeO]
                mi_moles_loop.drop(labels=reverse_FeO, axis=0, inplace=True)
                stepsize.loc[decrease_stepsize_FeO] = stepsize.loc[decrease_stepsize_FeO].div(decrease_factor)

            mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
            mi_wtPercent = mi_moles.convert_moles_wtPercent
            FeO = mi_wtPercent["FeO"]
            if self._FeO_as_function:
                FeO_target = self.FeO_function(mi_wtPercent)
                self.FeO_target = FeO_target
            else:
                FeO_target = self.FeO_target
            FeO_mismatch_new = ~np.isclose(FeO, FeO_target, atol=FeO_converge)
            reslice = ~np.array_equal(FeO_mismatch_new, FeO_mismatch)
            FeO_mismatch = FeO_mismatch_new
            # Make sure that a dataframe sliced with FeO_mismatch always returns a dataframe
            if not isinstance(FeO_mismatch, np.ndarray):
                FeO_mismatch = np.array([FeO_mismatch])

        corrected_compositions = mi_moles.convert_moles_wtPercent            

        if inplace:
            self._FeO_restored = True
            self.inclusions = corrected_compositions
        else:
            return (
                corrected_compositions,
                self.olivine_corrected,
                {"old": temperature_old, "new": corrected_compositions.temperature(P_bar)},
            )

    def correct(self):

        i = self.Fe_equilibrate(inplace=True)

        return self.crystallisation_correction(i)
