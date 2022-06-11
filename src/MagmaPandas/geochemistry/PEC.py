from numpy import isin
import pandas as pd
import numpy as np
from functools import partial
from ..MagmaFrames import olivine, melt_inclusion, melt
import elements as e
from typing import Union
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg
from ..geochemistry.Fe_redox import Fe_redox
from ..geochemistry.fO2 import fO2_QFM


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

    cations = olivine(
        {"MgO": MgO, "FeO": FeO, "SiO2": SiO2, "total": total},
        index=idx,
        units="mol fraction",
        datatype="oxide",
    )
    cations = cations.convert_moles_wtPercent

    return cations


def Fe_restore(
    inclusion: Union[melt_inclusion, melt],
    FeO_initial: Union[int, float],
    P_bar: float,
    **kwargs
):
    """
    Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
    Isothermal and isobaric.
    """
    # Grab model parameters
    converge = kwargs.get("converge", 0.02)  # In wt. %
    temperature_converge = kwargs.get("temperature_converge", 0.2)  # In degrees
    stepsize = kwargs.get("stepsize", 0.001)  # In molar fraction
    # Select Fe loss or gain
    if inclusion["FeO"] > FeO_initial:
        stepsize = -stepsize
    forsterite = kwargs.get("forsterite", 0.8)
    QFM_logshift = kwargs.get("QFM_logshift", 1)
    # Parameters for the while loop
    olivine_melted = 0
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
    moles = melt_inclusion(
        columns=inclusion.index, units="mol fraction", datatype="oxide"
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
        olivine = pd.Series(
            {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1},
            index=moles.columns,
        ).fillna(0.0)

        # Crystallise or melt olivine to remain isothermal
        add_olivine = moles.iloc[-1]
        # Set stepsize for olivine addition/removal
        olivine_stepsize = stepsize / 4
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
    inclusion: Union[melt_inclusion, melt],
    forsterite_host: float,
    P_bar: float,
    **kwargs
):
    """
    Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
    Isothermal and isobaric.
    """
    # Grab model parameters
    stepsize = kwargs.get("stepsize", -0.005)  # In molar fraction
    converge = kwargs.get("converge", 0.002)  # In molar fraction
    temperature_converge = kwargs.get("temperature_converge", 0.1)  # In degrees
    QFM_logshift = kwargs.get("QFM_logshift", 1)
    # Parameters for the while loop
    olivine_crystallised = 0
    olivine_stepsize_reduction = 4
    decrease_factor = 10
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
    def calc_forsterite_EQ(
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
        Fe2Mg = mol_fractions["FeO"] * Fe2_FeTotal / mol_fractions["MgO"]
        # Equilibrium forsterite content
        return 1 / (1 + Kd * Fe2Mg)

    # Equilibrium forsterite content
    forsterite_EQ = calc_forsterite_EQ(inclusion.moles)
    # Select Fe removal or addition
    if forsterite_EQ > forsterite_host:
        stepsize = -stepsize

    # Set up initial data
    moles = melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    moles.loc[0] = inclusion.moles[inclusion.elements].values
    # Fe-Mg exchange vector
    FeMg_exchange = pd.Series(0, index=moles.columns)
    FeMg_exchange.loc[["FeO", "MgO"]] = 1, -1

    while not np.isclose(forsterite_EQ, forsterite_host, atol=converge):
        # Exchange Fe-Mg
        idx = moles.index[-1] + stepsize
        moles.loc[idx] = (moles.iloc[-1] + FeMg_exchange.mul(stepsize)).values
        # New model temperature after Fe-Mg exchange
        temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )

        # Equilibrium forsterite content for new composition
        forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])
        # Equilibrium olivine composition in oxide mol fractions
        olivine = pd.Series(
            {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1},
            index=moles.columns,
        ).fillna(0.0)

        # Crystallise or melt olivine to remain isothermal
        remove_olivine = moles.iloc[-1]
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
            # Reverse one iteration and reduce stepsize if temperature was
            # overstepped by more than the convergence value
            Temperature_mismatch = ~ np.isclose(
                temperature_new, temperature, atol=temperature_converge
            )
            decrease_stepsize_T = np.logical_and(T_overstepped, Temperature_mismatch)
            if decrease_stepsize_T :
                remove_olivine = remove_olivine - olivine * olivine_stepsize
                olivine_crystallised -= olivine_stepsize
                olivine_stepsize = olivine_stepsize / decrease_factor

        # Copy olivine corrected composition
        moles.iloc[-1] = remove_olivine.values

        # New equilibrium forsterite content
        forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])

        disequilibrium = ~ np.isclose(forsterite_EQ, forsterite_host, atol=converge)
        overstepped = np.sign(forsterite_EQ - forsterite_host) != np.sign(stepsize)
        decrease_stepsize = np.logical_and(disequilibrium, overstepped)
        # Reverse one iteration and reduce stepsize if forsterite content
        # gets oversteppend by more than the convergence value
        if decrease_stepsize:
            moles.drop([idx], inplace=True)
            forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])
            stepsize = stepsize / decrease_factor

    # Recalculate compositions to oxide wt. %
    wtPercent = moles.convert_moles_wtPercent

    return wtPercent, forsterite_host, forsterite_EQ, olivine_crystallised
