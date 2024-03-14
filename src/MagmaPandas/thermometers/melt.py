"""
Sub-module with melt-only thermometers
"""

import warnings as w

import elementMass as e
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma
from MagmaPandas.thermometers.data_parsing import (
    _anhydrous_composition,
    _get_oxides,
    parse_data,
)


def _check_calibration_range_Series(melt: pd.Series) -> None:

    SiO2_range = not (31 < melt["SiO2"] < 73.64)
    alkali_range = melt[["Na2O", "K2O"]].sum() > 14.3
    H2O_range = melt["H2O"] > 18.6

    outside_range = SiO2_range | alkali_range | H2O_range

    if not outside_range:
        return
    w.warn(f"sample has composition outside the thermometer calibration range")


def _check_calibration_range_dataframe(melt: Magma) -> None:

    SiO2_range = ~melt["SiO2"].between(31, 73.64)
    alkali_range = melt[["Na2O", "K2O"]].sum(axis=1) > 14.3
    H2O_range = melt["H2O"] > 18.6

    outside_range = SiO2_range | alkali_range | H2O_range

    if outside_range.sum() < 1:
        return

    w.warn(
        f"samples {*melt.index[outside_range],} have compositions outside the thermometer calibration range"
    )


def _check_calibration_range(melt: pd.Series | Magma) -> None:

    if isinstance(melt, pd.Series):
        _check_calibration_range_Series(melt)
    elif isinstance(melt, pd.DataFrame):
        _check_calibration_range_dataframe(melt)


def putirka2008_14(melt: Magma, *args, **kwargs) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 14 from Putirka (2008)\ [15]_ calculates liquidus temperatures based on melt compositions.
    Requires equilibrium with olivine.

    SEE = 58 degrees

    Applicable between:

    -    0 - 14.4 GPa

    -    729 - 2000 degrees C

    -    31 - 73.64 wt. % SiO\ :sub:`2`

    -    0 - 14.3 wt. % Na\ :sub:`2`\ O + K\ :sub:`2`\ O

    -    0 - 18.6 wt. % H\ :sub:`2`\ O

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """

    _check_calibration_range(melt)

    oxides = _get_oxides(melt)
    oxides_needed = set(["MgO", "FeO", "Na2O", "K2O"])

    composition = parse_data(melt, oxides_needed)

    if "H2O" not in oxides:
        H2O = 0.0
    else:
        H2O = composition["H2O"].copy()
        composition = _anhydrous_composition(composition)

    # Calculate molar oxide fractions
    mol_fractions = composition.moles
    # Melt Mg#
    Mg_no = mol_fractions["MgO"] / (
        mol_fractions["MgO"] + mol_fractions["FeO"]
    )  # Putirka doesn't specify whether this should be Fe2+ or Fe(total)

    if "Fe2O3" in oxides:
        FeO_w, Fe2O3_w = e.compound_weights(["FeO", "Fe2O3"])
        composition["FeO"] = composition["FeO"] + (
            2 * composition["Fe2O3"] * (FeO_w / Fe2O3_w)
        )

    T_K = (
        754
        + 190.6 * Mg_no
        + 25.52 * composition["MgO"]
        + 9.585 * composition["FeO"]
        + 14.87 * (composition["Na2O"] + composition["K2O"])
        - 9.176 * H2O
    ) + 273.15

    return pd.Series(T_K, name="T_K").squeeze()


def putirka2008_15(
    melt: Magma, P_bar: float | pd.Series = None, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 15 from Putirka (2008)\ [15]_ calculates liquidus temperatures based on melt compositions.
    Requires equilibrium with olivine.

    SEE = 46 degrees

    Applicable between:

    -    0 - 14.4 GPa

    -    729 - 2000 degrees C

    -    31 - 73.64 wt. % SiO\ :sub:`2`

    -    0 - 14.3 wt. % Na\ :sub:`2`\ O + K\ :sub:`2`\ O

    -    0 - 18.6 wt. % H\ :sub:`2`\ O

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %

    P_bar : float, pandas Series
        pressures in bar.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """
    _check_calibration_range(melt)

    if P_bar is None:
        P_bar = melt["P_bar"]

    oxides = _get_oxides(melt)
    oxides_needed = set(["MgO", "FeO", "Na2O", "K2O"])
    composition = parse_data(melt, oxides_needed)

    if "H2O" not in oxides:
        H2O = 0.0
    else:
        H2O = composition["H2O"].copy()
        composition = _anhydrous_composition(composition)

    P_GPa = P_bar / 1e4

    # Calculate molar oxide fractions
    mol_fractions = composition.moles
    # Melt Mg#
    Mg_no = mol_fractions["MgO"] / (
        mol_fractions["MgO"] + mol_fractions["FeO"]
    )  # Putirka doesn't specify whether this should be Fe2+ or Fe(total)

    if "Fe2O3" in oxides:
        FeO_w, Fe2O3_w = e.compound_weights(["FeO", "Fe2O3"])
        composition["FeO"] = composition["FeO"] + (
            2 * composition["Fe2O3"] * (FeO_w / Fe2O3_w)
        )

    T_K = (
        815.3
        + 265.5 * Mg_no
        + 15.37 * composition["MgO"]
        + 8.61 * composition["FeO"]
        + 6.646 * (composition["Na2O"] + composition["K2O"])
        + 39.16 * P_GPa
        - 12.83 * H2O
    ) + 273.15

    return pd.Series(T_K, name="T_K").squeeze()


def putirka2008_16(
    melt: Magma, P_bar: float | pd.Series = None, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 16 from Putirka (2008)\ [15]_ calculates liquiqdus temperatures based on melt compositions.
    Requires equilibrium with olivine + plagioclase + clinopyroxene.

    SEE = 19 degrees

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %

    P_bar : float, pandas Series
        pressures in bar

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """
    oxides = _get_oxides(melt)
    oxides_needed = set(["SiO2", "Al2O3", "MgO"])

    composition = parse_data(melt, oxides_needed)
    # import MagmaPandas as mp

    if P_bar is None:
        P_bar = melt["P_bar"]

    if isinstance(P_bar, pd.Series):
        if not melt.index.equals(P_bar.index):
            raise RuntimeError("Melt and P_bar indices don't match")

    if "H2O" in oxides:
        composition = _anhydrous_composition(composition)

    # Convert pressure from bars to GPa
    P_GPa = P_bar / 1e4

    # Calculate molar oxide fractions
    mol_fractions = composition.moles

    part_1 = (
        -583
        + 3141 * mol_fractions["SiO2"]
        + 15779 * mol_fractions["Al2O3"]
        + 1338.6 * mol_fractions["MgO"]
    )
    part_2 = -31440 * mol_fractions["SiO2"] * mol_fractions["Al2O3"] + 77.67 * P_GPa

    T_K = part_1 + part_2 + 273.15

    return pd.Series(T_K, name="T_K").squeeze()
