from functools import partial
from typing import Union

import elementMass as e
import numpy as np
import pandas as pd
from scipy.optimize import root_scalar

from MagmaPandas.configuration import configuration
from MagmaPandas.Fe_redox import Fe_redox
from MagmaPandas.fO2 import fO2_QFM
from MagmaPandas.Kd.Ol_melt import calculate_FeMg_Kd
from MagmaPandas.MagmaFrames import Melt
from MagmaPandas.MagmaSeries import MagmaSeries
from MagmaPandas.PEC.PEC_configuration import PEC_configuration
from MagmaPandas.PEC.PEC_functions import _root_temperature


def Fe_equilibrate(
    inclusion: Melt,
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
        "Kd_converge", getattr(PEC_configuration, "Kd_converge")
    )  # Kd converge
    Fe2_behaviour = getattr(PEC_configuration, "_Fe2_behaviour")
    dQFM = kwargs.get("dQFM", configuration.dQFM)
    # Parameters for the while loop
    olivine_crystallised = np.array([0.0])
    decrease_factor = getattr(PEC_configuration, "decrease_factor")
    # Normalise inclusion composition
    inclusion = inclusion.fillna(0.0)
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.normalise()
    # Calculate temperature and fO2
    temperature = inclusion.melt_temperature(P_bar=P_bar)
    fO2 = fO2_QFM(dQFM, temperature, P_bar)
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
    Kd_model = calculate_FeMg_Kd
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
        m_fractions = mol_fractions.copy()
        m_fractions = m_fractions.normalise()
        if Fe2_behaviour == "closed system":
            Fe2_FeTotal = 1
            Fe3Fe2 = mol_fractions.loc["Fe2O3"] * 2 / mol_fractions.loc["FeO"]
        elif Fe2_behaviour == "buffered":
            Fe3Fe2 = Fe3Fe2_model(m_fractions, temperature, fO2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)

        return Fe3Fe2, Fe2_FeTotal

    # Fix some parameters for Kd calculation
    calculate_Kd = partial(
        Kd_model, forsterite_initial=forsterite, T_K=temperature, P_bar=P_bar
    )
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
        inclusion.recalculate(inplace=True)
    # Calculate moles
    mi_moles = Melt(columns=inclusion.elements, units="mol fraction", datatype="oxide")
    mi_moles.loc[0] = inclusion.moles[inclusion.elements].values
    mi_moles = mi_moles.normalise()
    # Equilibrium Kd
    Kd_equilibrium = calculate_Kd(Melt_mol_fractions=mi_moles.iloc[-1], Fe3Fe2=Fe3Fe2)
    # Real Kd
    olivine_MgFe = forsterite / (1 - forsterite)
    melt_MgFe = mi_moles.loc[0, "MgO"] / (mi_moles.loc[0, "FeO"] * Fe2_FeTotal)
    Kd_real = melt_MgFe / olivine_MgFe
    # Fe-Mg exchange vector
    FeMg_vector = pd.Series(0, index=mi_moles.columns, name="FeMg_exchange")
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
        Kd_equilibrium = calculate_Kd(
            Melt_mol_fractions=mi_moles.normalise().loc[idx], Fe3Fe2=Fe3Fe2
        )
        melt_FeMg = (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal) / mi_moles.loc[idx, "MgO"]
        # Equilibrium olivine composition in oxide mol fractions
        Fo_equilibrium = 1 / (1 + Kd_equilibrium * melt_FeMg)
        olivine = MagmaSeries(
            {"MgO": Fo_equilibrium * 2, "FeO": (1 - Fo_equilibrium) * 2, "SiO2": 1},
            index=mi_moles.columns,
        )
        olivine = olivine.fillna(0.0).normalise()

        ######################################################
        # Add or remove olivine to keep temperature constant #
        olivine_amount = root_scalar(
            _root_temperature,
            args=(mi_moles.loc[idx], olivine, temperature, P_bar),
            x0=0,
            x1=0.05,
        ).root
        mi_moles.loc[idx] = mi_moles.loc[idx] + olivine.mul(olivine_amount)
        olivine_crystallised = np.append(
            olivine_crystallised, olivine_crystallised[-1] + olivine_amount
        )
        # mi_moles = mi_moles.normalise()
        temperature_new = mi_moles.loc[idx].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )
        ######################################################
        # New equilibrium Kd and Fe speciation
        Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
        Kd_equilibrium = calculate_Kd(
            Melt_mol_fractions=mi_moles.normalise().loc[idx], Fe3Fe2=Fe3Fe2
        )
        # Real Kd
        melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe
        # Assess equilibrium
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
            Kd_equilibrium = calculate_Kd(
                Melt_mol_fractions=mi_moles.normalise().iloc[-1], Fe3Fe2=Fe3Fe2
            )
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

    equilibrated_composition.index = olivine_crystallised

    return (
        equilibrated_composition,
        olivine_crystallised,
        {"old": temperature, "new": temperature_new},
        {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
    )
