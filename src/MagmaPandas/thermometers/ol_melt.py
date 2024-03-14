"""
Sub-module with olivine-melt thermometers
"""

import numpy as np
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma


def putirka2007_4(
    liquid: Magma, olivine: Magma, P_bar: float | pd.Series, **kwargs
) -> pd.Series:
    """
    Olivine-liquid thermometer

    Equation 4 from Putirka (2007)\ [16]_ calculates liquidus temperatures based on olivine-melt pairs.
    Identical to equation 22 from Putirka (2008)\ [15]_.

    SEE = hydrous: 29, anhydrous: 45, total: 43

    Total Fe is expressed as FeO

    Parameters
    ----------
    liquid  : Magma
        melt compositions in oxide wt. %
    olivine : pandas Series, pandas Dataframe
        olivine compositions in oxide wt. %
    P_bar : float, pandas Series
        pressures in bar

    Returns
    -------
    temperatures : pd.Series
        olivine liquidus temperatures in Kelvin.
    """
    from ..MagmaFrames import MagmaFrame
    from ..MagmaSeries import MagmaSeries

    composition = liquid.copy(deep=True)
    composition = composition.fillna(0.0)

    if isinstance(composition, MagmaFrame):
        elements = composition.columns
    elif isinstance(composition, MagmaSeries):
        elements = composition.index

    oxides = set(["MgO", "FeO", "Na2O", "K2O", "CaO", "SiO2", "K2O", "TiO2"])
    optional_oxides = {"MnO", "CoO", "NiO"}
    absentOxides = oxides.difference(elements)
    fill_oxides = optional_oxides.difference(elements)

    if len(absentOxides) > 0:
        raise KeyError(f"{absentOxides} not found in melt")

    # All oxides as mole/cation fractions on an anhydrous basis
    if len(fill_oxides) > 0:
        composition[list(fill_oxides)] = 0.0

    if "H2O" not in elements:
        H2O = 0.0
    else:
        H2O = composition["H2O"].copy()
        # moles are calculated on an anhydrous basis
        try:
            composition = composition.drop("H2O")
        except KeyError:
            composition = composition.drop(columns=["H2O"])

    # weight % oxides not renormalised to 100%, even when hydrous
    composition = composition.recalculate()

    # Calculate molar oxide fractions
    liquid_moles = composition.moles
    liquid_cations = composition.cations
    olivine_cations = olivine.cations

    C_NM = liquid_cations[["Fe", "Mn", "Mg", "Ca", "Co", "Ni"]].sum(
        axis=1
    )  # Fe = Fe total
    NF = 7 / 2 * np.log(1 - liquid_cations["Al"]) + 7 * np.log(1 - liquid_cations["Ti"])
    D_Mg = olivine_cations.loc[liquid_cations.index, "Mg"] / liquid_cations["Mg"]

    P_GPa = P_bar / 1e4

    numerator = 15294.6 + 1318.8 * P_GPa + 2.4834 * (P_GPa**2)
    denominator_a = 8.048 + 2.8532 * np.log(D_Mg) + 2.097 * np.log(1.5 * C_NM)
    denominator_b = (
        2.575 * np.log(3 * liquid_moles["SiO2"]) - 1.41 * NF + 0.222 * H2O + 0.5 * P_GPa
    )

    denominator = denominator_a + denominator_b

    return (numerator / denominator) + 273.15
