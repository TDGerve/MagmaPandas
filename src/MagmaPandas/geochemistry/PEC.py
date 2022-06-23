from numpy import isin
import pandas as pd
import numpy as np
from functools import partial
import elements as e
from typing import Union
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg
from ..geochemistry.Fe_redox import Fe_redox
from ..geochemistry.fO2 import fO2_QFM
from MagmaPandas import MagmaSeries, Olivine, Melt_inclusion, Melt
from MagmaPandas.configuration import configuration


def calculate_olivine(forsterite):
    """
    Calculate olivine wt. % composition from forsterite contents
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
    FeO_initial: Union[int, float],
    P_bar: float,
    **kwargs,
):
    """
    Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
    Isothermal and isobaric.
    """
    # Grab model parameters
    converge = kwargs.get("converge", 0.02)  # In wt. %
    temperature_converge = kwargs.get("temperature_converge", 0.2)  # In degrees
    stepsize = kwargs.get("stepsize", 0.005)  # In molar fraction
    # Select Fe loss or gain
    if inclusion["FeO"] > FeO_initial:
        stepsize = -stepsize
    forsterite = kwargs.get("forsterite", 0.8)
    QFM_logshift = kwargs.get("QFM_logshift", 1)
    # Parameters for the while loop
    olivine_melted = 0
    olivine_stepsize_reduction = 4
    decrease_factor = 10
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

    while not np.isclose(FeO, FeO_initial, atol=converge):
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
            olivine_melted += olivine_stepsize
            T_overstepped = np.sign(temperature - temperature_new) != np.sign(stepsize)
            # Reverse one iteration and reduce stepsize if temperature was
            # overstepped by more than the convergence value
            if T_overstepped and not np.isclose(
                temperature_new, temperature, atol=temperature_converge
            ):
                add_olivine = add_olivine - olivine * olivine_stepsize
                olivine_melted -= olivine_stepsize
                olivine_stepsize = olivine_stepsize / decrease_factor
                continue

        # Copy olivine corrected composition
        moles.iloc[-1] = add_olivine.values

        # New inclusion FeO
        FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

        overstepped = np.sign(FeO_initial - FeO) != np.sign(stepsize)
        # Reverse one iteration and reduce stepsize if FeO content
        # gets oversteppend by more than the convergence value
        if overstepped and not np.isclose(FeO, FeO_initial, atol=converge):
            moles.drop([idx], inplace=True)
            FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]
            stepsize = stepsize / decrease_factor
            FeMg_exchange = FeMg_exchange.div(decrease_factor)

    # Recalculate compositions to oxide wt. %
    wtPercent = moles.convert_moles_wtPercent

    return wtPercent, temperature, temperature_new, olivine_melted

def Fe_equilibrate(
    inclusion: Union[Melt_inclusion, Melt],
    forsterite_host: float,
    P_bar: float,
    **kwargs,
):
    """
    Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
    Isothermal and isobaric.
    """
    # Grab model parameters
    stepsize = kwargs.get("stepsize", 0.001)  # In molar fraction, this is the maximum recommende stepsize.
    converge = kwargs.get("converge", 0.002)  # Kd converge
    temperature_converge = kwargs.get("temperature_converge", 0.1)  # In degrees
    QFM_logshift = kwargs.get("QFM_logshift", 1)
    # Parameters for the while loop
    olivine_crystallised = 0
    olivine_stepsize_reduction = 4
    decrease_factor = 5
    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.normalise()
    # Calculate temperature and fO2
    temperature = inclusion.melt_temperature(P_bar=P_bar)
    fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

    # Function for calculating equilibrium forsterite content
    def calculate_Kd(
        mol_fractions,
        T_K=temperature,
        P_bar=P_bar,
        oxygen_fugacity=fO2,
        Fo_initial=forsterite_host,
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
    olivine_MgFe = forsterite_host / (1 - forsterite_host)
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
            melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
            Kd_real = melt_MgFe / olivine_MgFe
            stepsize = stepsize / decrease_factor

    # Recalculate compositions to oxide wt. %
    equilibrated_composition = mi_moles.convert_moles_wtPercent

    return equilibrated_composition, olivine_crystallised, {"Equilibrium": Kd_equilibrium, "Real": Kd_real}


def crystallisation_correction(
    inclusion: Union[Melt_inclusion, Melt],
    olivine_host: Union[float, MagmaSeries],
    FeO_initial: Union[int, float],
    P_bar: float,
    **kwargs,
):

    """
    Correct a melt inclusion for post entrapment crystallisation or melting by 
    respectively melting or crystallising host olivine.
    Expects the melt inclusion is completely equilibrated with the host crystal.
    The models exits when the user input original melt inclusion FeO content is reached.
    Loosely based on the postentrapment reequilibration procedure in Petrolog:

    L. V. Danyushesky and P. Plechov (2011)
    Petrolog3: Integrated software for modeling crystallization processes
    Geochemistry, Geophysics, Geosystems, vol 12
    """
    # Grab model parameters
    stepsize = kwargs.get("stepsize", 0.01)  # In molar fraction, this is the maximum recommended value
    converge = kwargs.get("converge", 0.05)  # FeO convergence
    Kd_converge = kwargs.get("Kd_converge", 0.001)  # Kd converge
    QFM_logshift = kwargs.get("QFM_logshift", 1)
    # Parameters for the while loop
    olivine_melted = 0
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
            {"MgO": olivine_host * 2, "FeO": (1 - olivine_host) * 2, "SiO2": 1}
        )
        olivine = olivine.normalise().reindex(mi_moles.columns, fill_value=0.0)
    elif not isinstance(olivine_host, MagmaSeries):
        raise TypeError(
            f"Olivine host should be forsterite number as float, or full composition as MagmaSeries, not {type(olivine_host)}"
        )
    else:
        olivine = olivine_host.moles
        olivine = olivine.reindex(mi_moles.columns, fill_value=0.0)
        olivine.recalculate()        
    forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])
    olivine_MgFe = olivine["MgO"] / olivine["FeO"]
    # Fe-Mg exchange vector
    FeMg_vector = pd.Series(0, index=mi_moles.columns)
    FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
    # Inclusion starting FeO
    FeO = inclusion["FeO"]
    temperature_old = mi_moles.iloc[-1].convert_moles_wtPercent.melt_temperature(P_bar=P_bar)

    def calculate_Kd(
        melt, pressure=P_bar, fObuffer_shift=QFM_logshift, Fo=forsterite
    ):
        T_K = melt.convert_moles_wtPercent.melt_temperature(pressure)
        fO2 = fO2_QFM(fObuffer_shift, T_K, pressure)
        Fe3Fe2 = Fe3Fe2_model(melt, T_K, fO2)
        Fe2_FeTotal = 1 / (1 + Fe3Fe2)     

        return Kd_model(melt, Fo, T_K, P_bar, Fe3Fe2), Fe2_FeTotal

    Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.iloc[-1])
    Kd_real = mi_moles.loc[0, "MgO"] / (mi_moles.loc[0, "FeO"] * Fe2_FeTotal) / olivine_MgFe

    if FeO > FeO_initial:
        stepsize = -stepsize

    ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
    while not np.isclose(FeO, FeO_initial, atol=converge):

        idx = round(mi_moles.index[-1] + stepsize, 4)
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + olivine.mul(stepsize)).values
        mi_moles = mi_moles.normalise()

        olivine_melted += stepsize

        Kd_equilibrium, Fe2_FeTotal = calculate_Kd(mi_moles.loc[idx])
        melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe

        ###### FE-MG EXCHANGE LOOP #####
        stepsize_FeMg = stepsize / FeMg_exchange_reduction
        FeMg_exchange = mi_moles.loc[idx]
        while not np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge):

            FeMg_exchange += FeMg_vector.mul(stepsize_FeMg)
            FeMg_exchange = FeMg_exchange.normalise()
            Kd_equilibrium, Fe2_FeTotal = calculate_Kd(FeMg_exchange)
            melt_MgFe =  FeMg_exchange["MgO"] / (FeMg_exchange["FeO"] * Fe2_FeTotal)
            Kd_real = melt_MgFe / olivine_MgFe        

            FeMg_overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
            FeMg_mismatch = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge)
            decrease_stepsize_FeMg = np.logical_and(FeMg_overstepped, FeMg_mismatch)
            # Reverse one iteration and reduce stepsize if Kd
            # gets oversteppend by more than the convergence value 
            if decrease_stepsize_FeMg:
                FeMg_exchange -= FeMg_vector.mul(stepsize_FeMg)
                stepsize_FeMg = stepsize_FeMg / decrease_factor

        mi_moles.loc[idx] = FeMg_exchange
        mi_moles = mi_moles.normalise()

        mi_wtPercent = mi_moles.convert_moles_wtPercent
        FeO = mi_wtPercent.loc[idx, "FeO"]        

        disequilibrium = ~np.isclose(FeO, FeO_initial, atol=converge)
        overstepped = np.sign(FeO_initial - FeO) != np.sign(stepsize)
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

    return mi_wtPercent, olivine_melted, (Kd_real, Kd_equilibrium), (temperature_new, temperature_old)


