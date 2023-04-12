import elementMass as e
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma

from .data_parsing import _anhydrous_composition, _get_oxides, parse_data


def putirka2008_14(melt: Magma, **kwargs):
    """
    Liquid thermometer

    Equation 14 from Putirka (2008) calculates liquiqdus temperature for bulk compositions.
    Requires equilibrium with olivine.

    SEE = 58 degrees

    Parameters
    ----------

    P_bar : int or list-like
        crystallisation pressures. Indices need to match melt if using pd.Series.


    Returns
    -------
    pd.Series
        liquidus temperatures in degrees Kelvin.
    """
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


def putirka2008_15(melt, P_bar, **kwargs):
    """
    Liquid thermometer

    Equation 15 from Putirka (2008) calculates liquiqdus temperature for bulk compositions.
    Requires equilibrium with olivine.

    Parameters
    ----------

    P_bar : int or list-like
        crystallisation pressures. Indices need to match melt if using pd.Series.


    Returns
    -------
    pd.Series
        liquidus temperatures in degrees Kelvin.
    """
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


def putirka2008_16(melt, P_bar=None, **kwargs):
    """Liquid thermometer

    Equation 16 from Putirka (2008) calculates liquiqdus temperature for bulk compositions.
    Requires equilibrium with olivine + plagioclase + clinopyroxene.

    Parameters
    ----------

    P_bar : int or list-like
        crystallisation pressures. Indices need to match melt if using pd.Series.


    Returns
    -------
    pd.Series
        liquidus temperatures in degrees Kelvin.
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
