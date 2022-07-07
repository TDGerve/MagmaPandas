from re import A
import pandas as pd
import numpy as np
import elements as e
from typing import Union
from scipy.optimize import root_scalar

from MagmaPandas.geochemistry.Kd_ol_melt import Kd_FeMg, Kd_FeMg_vectorised
from MagmaPandas.geochemistry.Fe_redox import Fe_redox
from MagmaPandas.geochemistry.fO2 import fO2_QFM
from MagmaPandas import MagmaSeries, Olivine, Melt_inclusion, Melt
from MagmaPandas.configuration import configuration


Fe3_options = ["buffered", "incompatible"]

"""
ALL CONFIGURATION OBJECTS NEED TO BE REWRITTEN BECAUSE CLASS PROPERTIES DON'T WORK LIKE THIS.
IT EITHER NEEDS:
    - A METACLASS WHERE ALL CLASS VARIABLES AND PROPERTIES ARE DEFINED
    - A CUSTOM CLASSPROPERTY DECORATOR
    - ALL CONFIGURATION OBJECTS INITIALISED INSIDE THE MODULE __INIT__ WITH ALL ATTRIBUTES AS INSTANCE ATTRIBUTES
"""


class PEC_configuration:

    Fe2_behaviour = "buffered"
    # Initial stepsize
    stepsize_equilibration = 0.002
    stepsize_crystallisation = 0.05
    # Reduction factor for stepsize after overstepping
    decrease_factor = 5
    # Convergence values
    FeO_converge = 0.05
    Kd_converge = 1e-3
    temperature_converge = 0.1

    @classmethod
    def reset(cls):
        cls.Fe2_behaviour = "buffered"
        cls.stepsize_equilibration = 0.002
        cls.stepsize_crystallisation = 0.05
        cls.decrease_factor = 5
        cls.FeO_converge = 0.05
        cls.Kd_converge = 1e-4
        cls.temperature_converge = 0.1

    @classmethod
    def print(cls):
        """ """

        variables = {
            "Fe2+ behaviour": "Fe2_behaviour",
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
    Kd_converge = kwargs.get(
        "converge", getattr(PEC_configuration, "Kd_converge")
    )  # Kd converge
    temperature_converge = kwargs.get(
        "temperature_converge", getattr(PEC_configuration, "temperature_converge")
    )  # In degrees
    Fe2_behaviour = kwargs.get("Fe2_behaviour", PEC_configuration().Fe2_behaviour)
    QFM_logshift = kwargs.get("QFM_logshift", configuration().QFMlogshift)

    # Parameters for the while loop
    olivine_crystallised = np.array([0.0])
    olivine_stepsize_reduction = 0.25
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

    def calculate_Fe2(
        mol_fractions, Fe2_behaviour=Fe2_behaviour, temperature=temperature, fO2=fO2
    ):

        if Fe2_behaviour == "closed system":
            Fe2_FeTotal = 1
            Fe3Fe2 = mol_fractions.loc["Fe2O3"] * 2 / mol_fractions.loc["FeO"]
        elif Fe2_behaviour == "buffered":
            Fe3Fe2 = Fe3Fe2_model(mol_fractions, temperature, fO2)
            Fe2_FeTotal = 1 - Fe3Fe2

        return Fe3Fe2, Fe2_FeTotal

    # Function for calculating equilibrium forsterite content
    def calculate_Kd(
        mol_fractions,
        Fe3Fe2,
        T_K=temperature,
        P_bar=P_bar,
        # oxygen_fugacity=fO2,
        Fo_initial=forsterite,
    ):
        return Kd_model(mol_fractions, Fo_initial, T_K, P_bar, Fe3Fe2)

    # Calculate initial Fe speciation
    Fe3Fe2, Fe2_FeTotal = calculate_Fe2(inclusion.moles, Fe2_behaviour="buffered")
    # For non-buffered Fe2+
    if Fe2_behaviour == "closed system":
        # Calculate Fe speciation
        Fe3_FeTotal = 1 - Fe2_FeTotal
        # Calculate Fe2O3 wt. %
        Fe2O3 = inclusion["FeO"] * Fe3_FeTotal
        FeO_mass, Fe2O3_mass = e.compound_weights(["FeO", "Fe2O3"])
        inclusion["Fe2O3"] = Fe2O3 * (FeO_mass * 2 / Fe2O3_mass)
        # Recalculate FeO
        inclusion["FeO"] = inclusion["FeO"] * (1 - Fe3_FeTotal)
        inclusion.recalculate()

    # Calculate moles
    mi_moles = Melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    mi_moles.loc[0] = inclusion.moles[inclusion.elements].values
    mi_moles = mi_moles.normalise()
    # Equilibrium Kd
    Kd_equilibrium = calculate_Kd(mi_moles.iloc[-1], Fe3Fe2)
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
    while not np.isclose(Kd_real, Kd_equilibrium, atol=Kd_converge, rtol=0):

        # Exchange Fe-Mg
        idx = mi_moles.index[-1] + stepsize
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + FeMg_vector.mul(stepsize)).values

        # Equilibrium Kd and Fe speciation for new composition
        Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
        Kd_equilibrium = calculate_Kd(mi_moles.loc[idx], Fe3Fe2)
        melt_FeMg = (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal) / mi_moles.loc[idx, "MgO"]
        # Equilibrium olivine composition in oxide mol fractions
        Fo_equilibrium = 1 / (1 + Kd_equilibrium * melt_FeMg)
        olivine = MagmaSeries(
            {"MgO": Fo_equilibrium * 2, "FeO": (1 - Fo_equilibrium) * 2, "SiO2": 1},
            index=mi_moles.columns,
        )
        olivine = olivine.fillna(0.0).normalise()

        ###################################
        # Add or remove olivine to keep temperature constant
        olivine_amount = root_scalar(
            _root_temperature,
            args=(mi_moles.loc[idx], olivine, temperature, P_bar),
            x0=0,
            x1=1,
        ).root
        mi_moles.loc[idx] = mi_moles.loc[idx] + olivine.mul(olivine_amount)
        olivine_crystallised = np.append(
            olivine_crystallised, olivine_crystallised[-1] + olivine_amount
        )
        mi_moles = mi_moles.normalise()
        temperature_new = mi_moles.loc[idx].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )
        ###################################

        # ##### OLIVINE CRYSTALLISATION LOOP #####
        # ########################################
        # remove_olivine = mi_moles.loc[idx]
        # # Set stepsize
        # olivine_stepsize = stepsize / olivine_stepsize_reduction
        # while not np.isclose(temperature_new, temperature, atol=temperature_converge, rtol=0):
        #     # Melt or crystallise olivine until the temperature is back to original.

        #     remove_olivine = remove_olivine + olivine.mul(olivine_stepsize)
        #     temperature_new = remove_olivine.convert_moles_wtPercent.melt_temperature(
        #         P_bar=P_bar
        #     )
        #     # Record removed/added olivine amount. Negative values for crystallisation
        #     olivine_crystallised[-1] += olivine_stepsize

        #     T_overstepped = np.sign(temperature - temperature_new) != np.sign(stepsize)
        #     # Reverse one iteration and reduce stepsize if temperature is
        #     # overstepped by more than the convergence value
        #     Temperature_mismatch = ~np.isclose(
        #         temperature_new, temperature, atol=temperature_converge, rtol=0
        #     )
        #     decrease_stepsize_T = np.logical_and(T_overstepped, Temperature_mismatch)
        #     if decrease_stepsize_T:
        #         remove_olivine = remove_olivine - olivine * olivine_stepsize
        #         # if len(olivine_crystallised) > 1:
        #         #     olivine_crystallised = olivine_crystallised[:-1]
        #         # else:
        #         olivine_crystallised[-1] -= olivine_stepsize
        #         olivine_stepsize = olivine_stepsize / decrease_factor

        # # Copy olivine corrected composition
        # mi_moles.loc[idx] = remove_olivine.values
        # mi_moles = mi_moles.normalise()
        # ##################################################

        # New equilibrium Kd and Fe speciation
        Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
        Kd_equilibrium = calculate_Kd(mi_moles.loc[idx], Fe3Fe2)
        # Real Kd
        melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe

        disequilibrium = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0)
        overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
        decrease_stepsize = np.logical_and(disequilibrium, overstepped)
        # Reverse one iteration and reduce stepsize if Kd
        # gets oversteppend by more than the convergence value
        if decrease_stepsize:
            mi_moles.drop(index=idx, inplace=True)
            olivine_crystallised = olivine_crystallised[:-1]
            # Reset equilibrium and real Kd
            Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
            Kd_equilibrium = calculate_Kd(mi_moles.iloc[-1], Fe3Fe2)
            idx = mi_moles.index[-1]
            melt_MgFe = mi_moles.loc[idx, "MgO"] / (
                mi_moles.loc[idx, "FeO"] * Fe2_FeTotal
            )
            Kd_real = melt_MgFe / olivine_MgFe
            stepsize = stepsize / decrease_factor

    # Recalculate compositions to oxide wt. %
    equilibrated_composition = mi_moles.convert_moles_wtPercent
    if len(olivine_crystallised) == 1:
        olivine_crystallised = np.array([0])
        temperature_new = temperature

    return (
        equilibrated_composition,
        olivine_crystallised,
        {"old": temperature, "new": temperature_new},
        {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
    )


def _root_temperature(olivine_amount, melt_x_moles, olivine_x_moles, T_K, P_bar):

    melt_x_new = melt_x_moles + olivine_x_moles.mul(olivine_amount)
    temperature_new = melt_x_new.convert_moles_wtPercent.melt_temperature(P_bar=P_bar)

    return T_K - temperature_new


def _root_Kd(exchange_amount, melt_x_moles, exchange_vector, forsterite, P_bar, dQFM):

    melt_x_new = melt_x_moles + exchange_vector.mul(exchange_amount)
    Kd_equilibrium, Kd_real = calculate_Kds(melt_x_new, P_bar, dQFM, forsterite)

    return Kd_equilibrium - Kd_real


def calculate_Kds(melt_x_moles, P_bar, dQFM, forsterite):

    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Fe_redox, configuration().Kd_model)

    T_K = melt_x_moles.convert_moles_wtPercent.melt_temperature(P_bar)
    fO2 = fO2_QFM(dQFM, T_K, P_bar)
    Fe3Fe2 = Fe3Fe2_model(melt_x_moles, T_K, fO2)

    Fe2_FeTotal = 1 / (1 + Fe3Fe2)
    melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
    olivine_MgFe = forsterite / (1 - forsterite)
    Kd_observed = melt_MgFe / olivine_MgFe

    Kd_eq = Kd_model(melt_x_moles, forsterite, T_K, P_bar, Fe3Fe2)

    return Kd_eq, Kd_observed


def crystallisation_correction(
    inclusion: Union[Melt_inclusion, Melt],
    olivine_host: Union[float, MagmaSeries],
    FeO_target: Union[int, float, callable],
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
    calculate_FeO_target = False
    # Parameters for the while loop
    # olivine_corrected = pd.Series(0, index=[0])
    FeMg_exchange_reduction = 4
    decrease_factor = 5
    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.fillna(0.0)
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
        olivine = olivine.reindex(mi_moles.columns)
        olivine = olivine.fillna(0.0)
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

        Kd_eq = Kd_model(melt, Fo, T_K, P_bar, Fe3Fe2)

        return Kd_eq, Kd_observed

    if hasattr(FeO_target, "__call__"):
        calculate_FeO_target = FeO_target
        FeO_target = calculate_FeO_target(inclusion)

    if FeO > FeO_target:
        stepsize = -stepsize

    ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
    while not np.isclose(FeO, FeO_target, atol=converge, rtol=0):

        idx = round(mi_moles.index[-1] + stepsize, 4)
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + olivine.mul(stepsize)).values
        mi_moles = mi_moles.normalise()

        # olivine_corrected.loc[idx] = olivine_corrected.iloc[-1] + stepsize

        Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.loc[idx])

        ###### FE-MG EXCHANGE LOOP #####
        stepsize_FeMg = stepsize / FeMg_exchange_reduction
        FeMg_exchange = mi_moles.loc[idx]
        while not np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0):

            FeMg_exchange += FeMg_vector.mul(stepsize_FeMg)
            FeMg_exchange = FeMg_exchange.normalise()
            Kd_equilibrium, Kd_real = calculate_Kd(FeMg_exchange)

            Kd_overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
            Kd_mismatch = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0)
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

        if calculate_FeO_target:
            FeO_target = calculate_FeO_target(mi_wtPercent.loc[idx])

        disequilibrium = ~np.isclose(FeO, FeO_target, atol=converge, rtol=0)
        overstepped = np.sign(FeO_target - FeO) != np.sign(stepsize)
        decrease_stepsize = np.logical_and(disequilibrium, overstepped)
        # Reverse one iteration and reduce stepsize if FeO
        # gets oversteppend by more than the convergence value
        if decrease_stepsize:
            mi_moles.drop(index=idx, inplace=True)
            # olivine_corrected.drop(index=idx, inplace=True)
            mi_wtPercent = mi_moles.convert_moles_wtPercent
            idx = mi_wtPercent.index[-1]
            FeO = mi_wtPercent.loc[idx, "FeO"]
            stepsize = stepsize / decrease_factor

    olivine_corrected = mi_moles.index.values
    if olivine_corrected.max() == 0:
        mi_wtPercent = mi_moles.convert_moles_wtPercent
        Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.iloc[-1])

    temperature_new = mi_wtPercent.iloc[-1].melt_temperature(P_bar)

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
            self.inclusions_uncorrected = self.inclusions.copy()
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
            self._olivine = olivine.normalise()
        else:
            olivines = olivines.fillna(0.0)
            self._olivine = olivines.moles
        self._olivine = self._olivine.reindex(
            columns=self.inclusions.columns, fill_value=0.0
        )
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

    def reset(self):
        self.inclusions = self.inclusions_uncorrected.copy()
        self.olivine_corrected = pd.Series(0, index=self.inclusions.index)
        if self._FeO_as_function:
            self.FeO_target = self.FeO_function(self.inclusions)

    @property
    def olivine(self):
        return self._olivine.convert_moles_wtPercent

    @olivine.setter
    def olivine(self, value):
        print("Olivine is read only")

    @property
    def equilibrated(self):
        """
        Check which inclusions are in equilibrium with their host olivine
        based on modelled and observed Fe-Mg exchange Kd's.
        """
        Kd_converge = getattr(PEC_configuration, "Kd_converge") * 1.5
        Kd_equilibrium, Kd_real = self.calculate_Kds()

        return pd.Series(
            np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0),
            index=Kd_equilibrium.index,
            name="equilibrated",
        )

    @property
    def Fe_loss(self):
        """
        Check which inclusions have experienced Fe loss
        """
        FeO_converge = getattr(PEC_configuration, "FeO_converge")
        if self._FeO_as_function:
            FeO_target = self.FeO_function(self.inclusions)
        else:
            FeO_target = self.FeO_target
        return pd.Series(
            ~np.isclose(
                FeO_target,
                self.inclusions["FeO"],
                atol=FeO_converge,
                rtol=0,
            ),
            index=self.inclusions.index,
        )

    def calculate_Kds(self, **kwargs):
        """
        Calculate observed and modelled Kds
        """
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)
        dQFM = configuration().QFMlogshift

        melt = kwargs.get("melt", self.inclusions.moles)
        pressure = kwargs.get("P_bar", self.P_bar)
        forsterite = kwargs.get("forsterite", self._olivine.forsterite)

        T_K = kwargs.get("T_K", melt.convert_moles_wtPercent.temperature(pressure))
        fO2 = kwargs.get("fO2", fO2_QFM(dQFM, T_K, pressure))
        Fe3Fe2 = Fe3Fe2_model(melt, T_K, fO2)

        Fe2_FeTotal = 1 / (1 + Fe3Fe2)
        melt_MgFe = melt["MgO"] / (melt["FeO"] * Fe2_FeTotal)
        olivine_MgFe = forsterite / (1 - forsterite)
        Kd_observed = melt_MgFe / olivine_MgFe
        Kd_observed.rename("real", inplace=True)
        Kd_observed.index.name = "sample"

        Kd_equilibrium = Kd_model(melt, forsterite, T_K, pressure, Fe3Fe2)
        Kd_equilibrium.rename("equilibrium", inplace=True)
        Kd_equilibrium.index.name = "sample"

        return Kd_equilibrium, Kd_observed

    def Fe_equilibrate(self, inplace=False):
        """
        Docstring
        """

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
        # Calculate temperature and fO2
        temperatures = self.inclusions.temperature(P_bar=P_bar)
        fO2 = fO2_QFM(dQFM, temperatures, P_bar)
        # Get olivine forsterite contents
        forsterite_host = self._olivine.forsterite

        def calculate_Fe2FeTotal(mol_fractions, T_K, oxygen_fugacity):

            Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
            return 1 / (1 + Fe3Fe2)

        # Set up initial data
        mi_moles = self.inclusions.moles
        Kd_equilibrium, Kd_real = self.calculate_Kds(
            melt=mi_moles,
            T_K=temperatures,
            P_bar=P_bar,
            fO2=fO2,
            forsterite=forsterite_host,
        )
        # Fe-Mg exchange vectors
        FeMg_exchange = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
        FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]
        # Find disequilibrium inclusions
        disequilibrium = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0)
        # Set stepsizes acoording to Kd disequilibrium
        stepsize.loc[Kd_real < Kd_equilibrium] = -stepsize.loc[
            Kd_real < Kd_equilibrium
        ].copy()

        ##### Main Fe-Mg exhange loop #####
        ###################################
        total_inclusions = mi_moles.shape[0]
        reslice = True
        while sum(disequilibrium) > 0:
            print(
                f"{sum(~disequilibrium):0>3}/{total_inclusions:0>3} inclusions equilibrated",
                end="\r",
            )
            if reslice:
                # Only reslice al the dataframes and series if the amount of disequilibrium inclusions has changed
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
            # mi_moles_loop = mi_moles_loop.normalise()
            # Calculate new liquidus temperature
            temperatures_new = mi_moles_loop.convert_moles_wtPercent.temperature(
                P_bar=P_loop
            )

            # Calculate new equilibrium Kd and Fe speciation
            Kd_eq_loop, Kd_real_loop = self.calculate_Kds(
                melt=mi_moles_loop,
                T_K=temperatures_loop,
                fO2=fO2_loop,
                P_bar=P_loop,
                forsterite=Fo_host_loop,
            )
            Fe2_FeTotal = calculate_Fe2FeTotal(
                mi_moles_loop, temperatures_loop, fO2_loop
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
                columns=mi_moles_loop.columns,
                units="mol fraction",
                datatype="oxide",
            )
            olivine = olivine.fillna(0.0).normalise()
            # Find inclusions outside the equilibrium temperature range
            temperature_mismatch = ~np.isclose(
                temperatures_new, temperatures_loop, atol=temperature_converge, rtol=0
            )
            ############################# Olivine melting/crystallisation loop ############################
            # Melt or crystallise olivine until the liquidus temperature is back in the equilibrium range #
            # Initialise data for the olivine melting/crystallisation loop
            olivine_stepsize = stepsize_loop.div(olivine_stepsize_reduction)
            # olivine_correction = mi_moles_loop.loc[temperature_mismatch, :].copy()
            reslice_olivine = True
            while sum(temperature_mismatch) > 0:
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
                    rtol=0,
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
                # Make sure that temperature mismatch stays list-like
                # so that any sliced dataframe remains a dataframe
                try:
                    len(temperature_mismatch)
                except TypeError:
                    temperature_mismatch = [temperature_mismatch]

            # Renormalise after olivine correction
            mi_moles_loop = mi_moles_loop.normalise()
            Kd_eq_loop, Kd_real_loop = self.calculate_Kds(
                melt=mi_moles_loop,
                T_K=temperatures_loop,
                P_bar=P_loop,
                fO2=fO2_loop,
                forsterite=Fo_host_loop,
            )
            # Find inclusions outside the equilibrium forsterite range
            disequilibrium_loop = ~np.isclose(
                Kd_eq_loop, Kd_real_loop, atol=Kd_converge, rtol=0
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
            disequilibrium_new = ~np.isclose(
                Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0
            )
            reslice = ~np.array_equal(disequilibrium_new, disequilibrium)
            disequilibrium = disequilibrium_new
            # Make sure that disequilibrium stays list-like
            # so that any sliced dataframe remains a dataframe
            try:
                len(disequilibrium)
            except TypeError:
                disequilibrium = [disequilibrium]

        print(
            f"{sum(~disequilibrium):0>3}/{total_inclusions:0>3} inclusions equilibrated",
            end="\n",
        )
        corrected_compositions = mi_moles.convert_moles_wtPercent

        self.inclusions = corrected_compositions

        if not inplace:
            return (
                corrected_compositions,
                self.olivine_corrected,
                {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
            )

    def correct_olivine(self, inplace=False, **kwargs):
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

        if not self.equilibrated.all():
            raise RuntimeError("Inclusion compositions have not all been equilibrated")

        # Get settings
        stepsize = getattr(PEC_configuration, "stepsize_crystallisation")
        stepsize = pd.Series(stepsize, index=self.inclusions.index)
        decrease_factor = getattr(PEC_configuration, "decrease_factor")
        FeMg_exchange_reduction = 4
        Kd_converge = getattr(PEC_configuration, "Kd_converge")
        FeO_converge = kwargs.get(
            "FeO_converge", getattr(PEC_configuration, "FeO_converge")
        )
        P_bar = self.P_bar
        # Inclusion compositions in oxide mol fractions
        mi_moles = self.inclusions.moles
        mi_moles = mi_moles.normalise()
        # Olivine forsterite and Mg/Fe
        forsterite = self._olivine.forsterite
        # Fe-Mg exchange vectors
        FeMg_vector = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
        FeMg_vector.loc[:, ["FeO", "MgO"]] = [1, -1]
        # Starting FeO and termperature
        FeO = self.inclusions["FeO"]
        FeO_target = self.FeO_target
        temperature_old = self.inclusions.temperature(P_bar=P_bar)

        stepsize.loc[FeO > FeO_target] = -stepsize.loc[FeO > FeO_target]
        FeO_mismatch = ~np.isclose(FeO, FeO_target, atol=FeO_converge, rtol=0)

        ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
        reslice = True
        total_inclusions = mi_moles.shape[0]
        while sum(FeO_mismatch) > 0:
            print(
                f"{sum(~FeO_mismatch):0>3}/{total_inclusions:0>3} inclusions corrected",
                end="\r",
            )

            if reslice:
                stepsize_loop = stepsize.loc[FeO_mismatch]
                mi_moles_loop = mi_moles.loc[FeO_mismatch, :].copy()
                idx_olivine = mi_moles_loop.index
                olivine_loop = self._olivine.loc[FeO_mismatch, :]
                Fo_loop = forsterite.loc[FeO_mismatch]
                P_loop = P_bar.loc[FeO_mismatch]

            mi_moles_loop = mi_moles_loop + olivine_loop.mul(stepsize_loop, axis=0)
            mi_moles_loop = mi_moles_loop.normalise()
            self.olivine_corrected.loc[idx_olivine] += stepsize_loop

            Kd_eq_loop, Kd_real_loop = self.calculate_Kds(
                melt=mi_moles_loop, P_bar=P_loop, forsterite=Fo_loop
            )
            Kd_mismatch = ~np.isclose(
                Kd_eq_loop, Kd_real_loop, atol=Kd_converge, rtol=0
            )

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

                mi_moles_loop.loc[
                    Kd_mismatch, :
                ] = FeMg_exchange + FeMg_vector_loop.mul(stepsize_FeMg_loop, axis=0)
                mi_moles_loop = mi_moles_loop.normalise()

                Kd_eq_loop, Kd_real_loop = self.calculate_Kds(
                    melt=mi_moles_loop, P_bar=P_loop, forsterite=Fo_loop
                )
                Kd_mismatch_new = ~np.isclose(
                    Kd_eq_loop, Kd_real_loop, atol=Kd_converge, rtol=0
                )
                Kd_overstepped = ~np.equal(
                    np.sign(Kd_real_loop - Kd_eq_loop), np.sign(stepsize_loop)
                )
                decrease_stepsize_Kd = np.logical_and(Kd_overstepped, Kd_mismatch_new)

                if sum(decrease_stepsize_Kd) > 0:
                    reverse_idx = mi_moles_loop.index[decrease_stepsize_Kd]
                    # Reverse Fe-Mg exchange
                    reverse_FeMg = FeMg_vector_loop.loc[reverse_idx, :].mul(
                        stepsize_FeMg_loop.loc[reverse_idx], axis=0
                    )
                    mi_moles_loop.loc[reverse_idx, :] = (
                        FeMg_exchange.loc[reverse_idx, :] - reverse_FeMg
                    )
                    stepsize_FeMg.loc[reverse_idx] = stepsize_FeMg.loc[reverse_idx].div(
                        decrease_factor
                    )

                reslice_FeMg = ~np.array_equal(Kd_mismatch, Kd_mismatch_new)
                Kd_mismatch = Kd_mismatch_new

            # Recalculate FeO
            mi_wtPercent = mi_moles_loop.convert_moles_wtPercent
            FeO = mi_wtPercent["FeO"]
            if self._FeO_as_function:
                FeO_target_loop = self.FeO_function(mi_wtPercent)
                self.FeO_target.loc[FeO_mismatch] = FeO_target_loop
            else:
                FeO_target_loop = self.FeO_target.loc[FeO_mismatch]
            # Find mismatched and overcorrected inclusions
            FeO_mismatch_loop = ~np.isclose(
                FeO, FeO_target_loop, atol=FeO_converge, rtol=0
            )
            FeO_overstepped = ~np.equal(
                np.sign(FeO_target_loop - FeO), np.sign(stepsize_loop)
            )
            decrease_stepsize_FeO = np.logical_and(FeO_overstepped, FeO_mismatch_loop)

            if sum(decrease_stepsize_FeO) > 0:
                # Reverse one step and decrease stepsize
                reverse_FeO = mi_moles_loop.index[decrease_stepsize_FeO]
                mi_moles_loop.drop(labels=reverse_FeO, axis=0, inplace=True)
                self.olivine_corrected.loc[reverse_FeO] -= stepsize_loop.loc[
                    reverse_FeO
                ]
                stepsize.loc[reverse_FeO] = stepsize.loc[reverse_FeO].div(
                    decrease_factor
                )
            # Copy loop data to the main variables
            mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
            mi_wtPercent = mi_moles.convert_moles_wtPercent
            FeO = mi_wtPercent["FeO"]
            if self._FeO_as_function:
                FeO_target = self.FeO_function(mi_wtPercent)
                self.FeO_target = FeO_target
            else:
                FeO_target = self.FeO_target
            FeO_mismatch_new = ~np.isclose(FeO, FeO_target, atol=FeO_converge, rtol=0)
            reslice = ~np.array_equal(FeO_mismatch_new, FeO_mismatch)
            FeO_mismatch = FeO_mismatch_new

        print(
            f"{sum(~FeO_mismatch):0>3}/{total_inclusions:0>3} inclusions corrected",
            end="\n",
        )
        corrected_compositions = mi_moles.convert_moles_wtPercent

        self.inclusions = corrected_compositions
        if not inplace:
            return (
                corrected_compositions,
                self.olivine_corrected,
                {
                    "old": temperature_old,
                    "new": corrected_compositions.temperature(P_bar),
                },
            )

    def correct(self):

        self.Fe_equilibrate(inplace=True)
        corrected_compositions, olivine_corrected, temperatures = self.correct_olivine(
            inplace=False
        )
        return corrected_compositions, olivine_corrected, temperatures

    def correct_inclusion(self, index, plot=True, **kwargs):

        if type(index) == int:
            index = self.inclusions_uncorrected.index[index]

        inclusion = self.inclusions_uncorrected.loc[index].copy()
        olivine = self._olivine.loc[index].copy()
        FeO_target = self.FeO_target.loc[index]
        P_bar = self.P_bar.loc[index]

        if self._FeO_as_function:
            FeO_target = self.FeO_function

        equilibrated, olivine_equilibrated, *_ = Fe_equilibrate(
            inclusion, olivine, P_bar, **kwargs
        )
        corrected, olivine_corrected, *_ = crystallisation_correction(
            equilibrated.iloc[-1].copy(), olivine, FeO_target, P_bar, **kwargs
        )
        total_corrected = olivine_corrected[-1] + olivine_equilibrated[-1]

        equilibrated["correction"] = "equilibration"
        corrected["correction"] = "correction"

        total_inclusion = pd.concat([equilibrated, corrected], axis=0)

        if self._FeO_as_function:
            FeO_target = self.FeO_function(total_inclusion.iloc[-1])

        if plot:
            import matplotlib.pyplot as plt

            fontsize = 14

            fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=False)

            plt.title(index)

            plt.plot(
                equilibrated["MgO"], equilibrated["FeO"], "-o", label="equilibration"
            )
            plt.plot(corrected["MgO"], corrected["FeO"], "-D", label="correction")

            middle = sum(ax.get_xlim()) / 2

            plt.axhline(FeO_target, linestyle="--", color="k", linewidth=1.5)
            plt.text(
                middle,
                FeO_target,
                s="FeO target",
                size=fontsize,
                backgroundcolor=ax.get_facecolor(),
            )
            plt.text(
                0.6,
                0.15,
                s=f"{total_corrected * 100:.1f} mol %\ncrystallisation correction",
                size=fontsize,
                transform=ax.transAxes,
            )

            ax.set_ylim(ax.get_ylim()[0], FeO_target * 1.03)

            ax.set_xlabel("wt. % MgO")
            ax.set_ylabel("wt. % FeO")

            plt.legend()

        return total_inclusion
