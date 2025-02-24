"""
Sub-module with melt-only thermometers
"""

import inspect
import sys

import elementMass as e
import numpy as np
import pandas as pd

from MagmaPandas.magma_protocol import Magma
from MagmaPandas.parse_io import check_components
from MagmaPandas.thermometers.data_parsing import (
    _anhydrous_composition,
    _check_calibration_range,
    _check_temperature,
    _get_elements,
)
from MagmaPandas.tools.modify_compositions import _remove_elements, moles_per_oxygen

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
    {
        "putirka2008_13": 71,
        "putirka2008_14": 58,
        "putirka2008_15": 46,
        "putirka2008_16": 26,
        "sun2020": 49,
        "shea2022": 13,  # on calibration dataset, not validated
        "sugawara2000_3": 33,  # on calibration dataset
        "sugawara2000_6a": 30,  # on calibration dataset
    }
)

components = {
    "putirka2008_13": ["MgO"],
    "putirka2008_14": ["MgO", "FeO", "Na2O", "K2O", "H2O"],
    "putirka2008_15": ["MgO", "FeO", "Na2O", "K2O", "H2O"],
    "putirka2008_16": ["SiO2", "Al2O3", "MgO"],
    "sun2020": ["MgO", "CaO", "K2O", "TiO2", "FeO", "CO2", "H2O"],
    "shea2022": ["MgO"],
    "sugawara2000_3": ["MgO"],
    "sugawara2000_6a": ["MgO", "FeO", "CaO", "SiO2"],
}


def putirka2008_13(
    melt: Magma, offset: float = 0.0, *args, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 13 from Putirka (2008)\ [28]_ calculates liquidus temperatures based on melt compositions.
    Requires saturation in olivine.

    SEE = 71 degrees

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

    composition = check_components(
        composition=melt, components=components["putirka2008_13"]
    )

    T_K = 26.3 * composition["MgO"] + 994.4 + 273.15

    _check_temperature(T_K)

    T_K = T_K + errors["putirka2008_13"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def putirka2008_14(
    melt: Magma, warn=True, offset: float = 0.0, *args, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 14 from Putirka (2008)\ [28]_ calculates liquidus temperatures based on melt compositions.
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

    elements = _get_elements(melt)
    composition = check_components(
        composition=melt, components=components["putirka2008_14"]
    )

    if warn:
        _check_calibration_range(
            melt=composition, calibration_range=calibration_range["putirka2008_14"]
        )
    # oxides_needed = set(["MgO", "FeO", "Na2O", "K2O"])

    # composition = parse_data(melt, oxides_needed)

    if "H2O" not in elements:
        H2O = 0.0
    else:
        H2O = composition["H2O"].copy()
        composition = _anhydrous_composition(composition)

    # Calculate molar oxide fractions
    mol_fractions = composition.moles()
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
    P_bar: float | pd.Series,
    offset: float = 0.0,
    warn=True,
    **kwargs,
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 15 from Putirka (2008)\ [28]_ calculates liquidus temperatures based on melt compositions.
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

    elements = _get_elements(melt)
    composition = check_components(
        composition=melt, components=components["putirka2008_15"]
    )

    if warn:
        _check_calibration_range(
            melt=composition, calibration_range=calibration_range["putirka2008_15"]
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
    mol_fractions = composition.moles()
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
    melt: Magma, P_bar: float | pd.Series, offset: float = 0.0, **kwargs
) -> float | pd.Series:
    """
    melt-only thermometer

    Equation 16 from Putirka (2008)\ [28]_ calculates liquiqdus temperatures based on melt compositions.
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

    if isinstance(P_bar, pd.Series):
        if not melt.index.equals(P_bar.index):
            raise RuntimeError("Melt and P_bar indices don't match")

    if "H2O" in elements:
        composition = _anhydrous_composition(composition)

    # Convert pressure from bars to GPa
    P_GPa = P_bar / 1e4

    # Calculate molar oxide fractions
    mol_fractions = composition.moles()

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


def sun2020(melt, P_bar: float | pd.Series, offset: float = 0.0, **kwargs):
    """
    Equation 6 from:

    Sun and Dasgupta (2020)\ [22]_

    Calibrated at:
    ~ 2 - 10 GPa
    ~ 950 - 1600 degrees C

    SEE: 49 degrees C

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %
    P_bar   :
        pressures in bar
    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """

    P_GPa = P_bar / 1e4

    composition = check_components(composition=melt, components=components["sun2020"])
    moles = composition.moles()

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


def shea2022(melt, offset: float = 0.0, **kwargs):
    """
    Equation 1 from Shea et al. (2022)\ [30]_

    Calibrated at:
    1 bar
    1060 - 1500 degrees C

    SEE: 13 degrees C (on calibration dataset, not validated)

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

    composition = check_components(composition=melt, components=components["shea2022"])

    T_K = 21.2 * composition["MgO"] + 1017 + 273.15

    _check_temperature(T_K)

    T_K = T_K + errors["shea2022"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def sugawara2000_3(melt, P_bar: float | pd.Series, offset: float = 0.0, **kwargs):
    """
    Equation 3 with olivine-liquid parameters and corrected for H2O according to equation 7a from Sugawara (2000)\ [31]_

    Calibrated at:
    <= 3.5 GPa
    1266 - 1873 C

    SEE: 33 degrees C

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %
    P_bar   :
        pressures in bar
    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """

    composition = check_components(
        composition=melt, components=components["sugawara2000_3"]
    )
    moles_percent = composition.moles() * 100

    A = 1293
    B = 14.60
    C = 5.5e-3
    T_K = A + B * moles_percent["MgO"] + C * P_bar

    if "H2O" in moles_percent.elements:
        T_K = T_K - 5.403 * moles_percent["H2O"]

    _check_temperature(T_K)

    T_K = T_K + errors["sugawara2000_3"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


def sugawara2000_6a(melt, P_bar: float | pd.Series, offset: float = 0.0, **kwargs):
    """
    Equation 6a corrected for H2O according to equation 7a from Sugawara (2000)\ [31]_
    Calibrated at:
    <= 3.5 GPa
    1266 - 1873 C

    SEE: 30 degrees C

    Parameters
    ----------
    melt : Magma
        melt compositions in oxide wt. %
    P_bar   :
        pressures in bar
    offset : float
        offset value in standard deviations. Temperatures are calculated as ``temnperature + offset * thermometer error (SEE)``.

    Returns
    -------
    temperatures : float, pandas Series
        liquidus temperatures in Kelvin.
    """

    composition = check_components(
        composition=melt, components=components["sugawara2000_6a"]
    )
    moles_percent = composition.moles() * 100

    T_K = (
        1466
        - 1.44 * moles_percent["SiO2"]
        - 0.5 * moles_percent["FeO"]
        + 12.32 * moles_percent["MgO"]
        - 3.899 * moles_percent["CaO"]
        + 4.3e-3 * P_bar
    )

    if "H2O" in moles_percent.elements:
        T_K = T_K - 5.403 * moles_percent["H2O"]

    _check_temperature(T_K)

    T_K = T_K + errors["sugawara2000_6a"] * offset

    return pd.Series(T_K, name="T_K").squeeze()


_module = sys.modules[__name__]
_funcs = inspect.getmembers(_module, inspect.isfunction)
# collect all melt thermometers in a dictionary
melt_thermometers_dict = {
    f[0]: f[1] for f in _funcs if f[1].__module__ == _module.__name__
}
