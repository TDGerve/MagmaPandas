"""
Sub-module with melt-only thermometers
"""

import warnings as w

import elementMass as e
import numpy as np
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma
from MagmaPandas.parse_io import check_components
from MagmaPandas.thermometers.data_parsing import (
    _anhydrous_composition,
    _check_calibration_range,
    _get_elements,
    _remove_elements,
    moles_per_oxygen,
)

calibration_range = {
    "putirka2008_14": [
        ["SiO2", 31, 73.64],
        [["Na2O", "K2O"], 0, 14.3],
        ["H2O", 0, 18.6],
    ],
    "putirka2008_15": [
        ["SiO2", 31, 73.64],
        [["Na2O", "K2O"], 0, 14.3],
        ["H2O", 0, 18.6],
    ],
}

errors = pd.Series(
    {"putirka2008_14": 58, "putirka2008_15": 46, "putirka2008_16": 26, "sun2020": 49}
)

components = {
    "putirka2008_14": ["MgO", "FeO", "Na2O", "K2O"],
    "putirka2008_15": ["MgO", "FeO", "Na2O", "K2O"],
    "putirka2008_16": ["SiO2", "Al2O3", "MgO"],
    "sun2020": ["MgO", "CaO", "K2O", "TiO2", "FeO", "CO2", "H2O"],
}


def _check_temperature(T_K):
    """
    Find negative or NaN temperatures
    """
    try:
        if (neg := any(T_K < 0)) or (nan := any(np.isnan(T_K))):
            str_arr = np.array(["negative", "NaN"])[np.array([neg, nan])]
            w.warn(f"{', '.join(str_arr)} temperatures found!")
    except TypeError:
        if (neg := (T_K < 0)) or (nan := (T_K != T_K)):
            str_arr = np.array(["negative", "NaN"])[np.array([neg, nan])]
            raise ValueError(f"{', '.join(str_arr)} temperature found!")


def putirka2008_14(
    melt: Magma, offset: float = 0.0, warn=True, *args, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 14 from Putirka (2008)\ [15]_ calculates liquidus temperatures based on melt compositions.
    Requires saturation in olivine.

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

    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """
    if warn:
        _check_calibration_range(
            melt=melt, calibration_range=calibration_range["putirka2008_14"]
        )

    elements = _get_elements(melt)
    composition = check_components(
        composition=melt, components=components["putirka2008_14"]
    )
    # oxides_needed = set(["MgO", "FeO", "Na2O", "K2O"])

    # composition = parse_data(melt, oxides_needed)

    if "H2O" not in elements:
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

    if "Fe2O3" in elements:
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

    _check_temperature(T_K)

    T_K = T_K + errors["putirka2008_14"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def putirka2008_15(
    melt: Magma,
    P_bar: float | pd.Series = None,
    offset: float = 0.0,
    warn=True,
    **kwargs,
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 15 from Putirka (2008)\ [15]_ calculates liquidus temperatures based on melt compositions.
    Requires saturation in olivine.

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

    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """
    if warn:
        _check_calibration_range(
            melt=melt, calibration_range=calibration_range["putirka2008_15"]
        )

    if P_bar is None:
        P_bar = melt["P_bar"]

    elements = _get_elements(melt)
    composition = check_components(
        composition=melt, components=components["putirka2008_15"]
    )
    # oxides_needed = set(["MgO", "FeO", "Na2O", "K2O"])
    # composition = parse_data(melt, oxides_needed)

    if "H2O" not in elements:
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

    if "Fe2O3" in elements:
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

    _check_temperature(T_K)

    T_K = T_K + errors["putirka2008_15"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def putirka2008_16(
    melt: Magma, P_bar: float | pd.Series = None, offset: float = 0.0, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 16 from Putirka (2008)\ [15]_ calculates liquiqdus temperatures based on melt compositions.
    Requires equilibrium with olivine + plagioclase + clinopyroxene and saturation with additional phases drastically increases the standard error of estimate.

    SEE = 26 degrees (saturation in olivine + plagioclase + clinopyroxene)
        = 60 degrees (saturation with additional phases)

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %

    P_bar : float, pandas Series
        pressures in bar

    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """
    elements = _get_elements(melt)
    composition = check_components(
        composition=melt, components=components["putirka2008_16"]
    )
    # oxides_needed = set(["SiO2", "Al2O3", "MgO"])

    # composition = parse_data(melt, oxides_needed)
    # import MagmaPandas as mp

    if P_bar is None:
        P_bar = melt["P_bar"]

    if isinstance(P_bar, pd.Series):
        if not melt.index.equals(P_bar.index):
            raise RuntimeError("Melt and P_bar indices don't match")

    if "H2O" in elements:
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

    _check_temperature(T_K)

    T_K = T_K + errors["putirka2008_16"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def sun2020(melt, P_bar, offset: float = 0.0, **kwargs):
    """
    Equation 4 from:

    Sun, C., Dasgupta, R. (2020) Thermobarometry of CO2-rich, silica-undersaturated melts constrains cratonic lithosphere thinning through time in areas of kimberlitic magmatism. Earth and Planetary Sience Letters. 550

    Calibrated at:
    ~ 2 - 10 GPa
    ~ 950 - 1600 degrees C
    """

    P_GPa = P_bar / 1e4

    composition = check_components(composition=melt, components=components["sun2020"])
    moles = composition.moles

    composition_volatile_free = _remove_elements(
        composition=moles, drop=["H2O", "CO2", "F", "S", "Cl"]
    )
    moles_unit_oxygen = moles_per_oxygen(moles=composition_volatile_free)

    Omega = (
        2.59
        + 3.5 * (moles_unit_oxygen["Ca1O"] - 2 * moles_unit_oxygen["K2O"])
        + 4.85 * moles_unit_oxygen["Ti1/2O"]
        + 1.4
        * (
            moles_unit_oxygen["Mg1O"]
            / (moles_unit_oxygen["Mg1O"] + moles_unit_oxygen["Fe1O"])
        )
        + 0.5 * moles_unit_oxygen["Mg1O"] * np.sqrt(composition["CO2"])
        + 5.7e-2 * composition["H2O"]
    )

    T_K = 1e4 / (
        Omega - 0.34 * np.sqrt(P_GPa) - 1.26 * np.log(moles_unit_oxygen["Mg1O"])
    )

    _check_temperature(T_K)

    T_K = T_K + errors["sun2020"] * offset

    return pd.Series(T_K, name="T_K").squeeze()
