from functools import partial
from typing import Union

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar

from MagmaPandas.configuration import configuration
from MagmaPandas.MagmaFrames import Melt
from MagmaPandas.MagmaSeries import MagmaSeries
from MagmaPandas.PEC.PEC_configuration import PEC_configuration
from MagmaPandas.PEC.PEC_functions import _root_Kd, calculate_Kds


def crystallisation_correction(
    inclusion: Melt,
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
    equilibration_crystallisation = kwargs.get("eq_crystallised", 0.0)
    # Grab model parameters
    stepsize = kwargs.get(
        "stepsize", getattr(PEC_configuration, "stepsize_crystallisation")
    )  # In molar fraction, 0.005 is the maximum recommended value
    converge = kwargs.get(
        "converge", getattr(PEC_configuration, "FeO_converge")
    )  # FeO convergence
    dQFM = kwargs.get("dQFM", configuration.dQFM)
    calculate_FeO_target = False
    # Parameters for the while loop
    decrease_factor = getattr(PEC_configuration, "decrease_factor")
    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements].copy()
    inclusion = inclusion.fillna(0.0)
    inclusion = inclusion.normalise()
    # SET UP INITIAL DATA
    # Dataframe with new compositions
    mi_moles = Melt(columns=inclusion.elements, units="mol fraction", datatype="oxide")
    mi_moles.loc[equilibration_crystallisation] = inclusion.moles[
        inclusion.elements
    ].values
    # Covert to moles left after equilibration
    mi_moles = mi_moles.mul((1 + equilibration_crystallisation), axis=0)
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
        olivine.recalculate(inplace=True)
        forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])
    # Fe-Mg exchange vector
    FeMg_vector = pd.Series(0, index=mi_moles.columns, name="FeMg_exchange")
    FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
    # Inclusion starting FeO
    FeO = inclusion["FeO"]
    temperature_old = mi_moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
        P_bar=P_bar
    )
    # Function to calculate Kds
    calculate_Kd = partial(calculate_Kds, P_bar=P_bar, forsterite=forsterite)

    if hasattr(FeO_target, "__call__"):
        calculate_FeO_target = FeO_target
        FeO_target = calculate_FeO_target(inclusion)

    if FeO > FeO_target:
        stepsize = -stepsize

    ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
    while not np.isclose(FeO, FeO_target, atol=converge, rtol=0):

        idx = mi_moles.index[-1] + stepsize
        mi_moles.loc[idx] = (mi_moles.iloc[-1] + olivine.mul(stepsize)).values
        # mi_moles = mi_moles.normalise()
        ###### FE-MG EXCHANGE UNTIL KD EQUILIBRIIUM #####
        exchange_amount = root_scalar(
            _root_Kd,
            args=(mi_moles.loc[idx], FeMg_vector, forsterite, P_bar, {"dQFM": dQFM}),
            x0=0,
            x1=0.05,
            # bracket=[-0.5,0.5],
        ).root
        mi_moles.loc[idx] = mi_moles.loc[idx] + FeMg_vector.mul(exchange_amount)
        # mi_moles = mi_moles.normalise()
        #################################################
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

    Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.iloc[-1])
    temperature_new = mi_wtPercent.iloc[-1].melt_temperature(P_bar)
    olivine_corrected = mi_moles.index.values
    if olivine_corrected.max() == 0:
        mi_wtPercent = mi_moles.convert_moles_wtPercent

    return (
        mi_wtPercent,
        olivine_corrected,
        {"equilibrium": Kd_equilibrium, "real": Kd_real},
        {"old": temperature_old, "new": temperature_new},
    )
