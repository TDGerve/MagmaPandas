from typing import Union

import numpy as np
import pandas as pd

from MagmaPandas.parse_io import check_components

"""
Parameterisation by:

Giordano, D., Russel, J. & Dingwell, D. (2008) Viscosity of magmatic liquids a model. Earth and planetary science letters. Vol 271, issue 1-4, pp. 123-135.
"""

parameters = pd.Series(
    {
        "A": -4.55,
        "b1": 159.6,
        "b2": -173.3,
        "b3": 72.1,
        "b4": 75.7,
        "b5": -39.0,
        "b6": -84.1,
        "b7": 141.5,
        "b11": -2.43,
        "b12": -0.91,
        "b13": 17.6,
        "c1": 2.75,
        "c2": 15.7,
        "c3": 8.3,
        "c4": 10.2,
        "c5": -12.3,
        "c6": -99.5,
        "c11": 0.3,
    }
)

components = [
    "SiO2",
    "TiO2",
    "Al2O3",
    "FeO",
    "MnO",
    "P2O5",
    "MgO",
    "CaO",
    "Na2O",
    "H2O",
    "F2",
    "K2O",
]


def _calculate_B(melt_mol_fractions: Union[pd.DataFrame, pd.Series]) -> pd.Series:

    axis = [0, 1][isinstance(melt_mol_fractions, pd.DataFrame)]

    mol_percent = check_components(melt_mol_fractions, components=components) * 100

    part_1 = (
        parameters["b1"] * mol_percent[["SiO2", "TiO2"]].sum(axis=axis)
        + parameters["b2"] * mol_percent["Al2O3"]
        + parameters["b3"] * mol_percent[["FeO", "MnO", "P2O5"]].sum(axis=axis)
        + parameters["b4"] * mol_percent["MgO"]
        + parameters["b5"] * mol_percent["CaO"]
        + parameters["b6"] * mol_percent[["Na2O", "H2O", "F2"]].sum(axis=axis)
        + parameters["b7"]
        * (mol_percent[["H2O", "F2"]].sum(axis=axis) + np.log(1 + mol_percent["H2O"]))
    )

    part_2 = (
        parameters["b11"]
        * (
            mol_percent[["SiO2", "TiO2"]].sum(axis=axis)
            * mol_percent[["FeO", "MnO", "MgO"]].sum(axis=axis)
        )
        + parameters["b12"]
        * (
            mol_percent[["SiO2", "TiO2", "Al2O3", "P2O5"]].sum(axis=axis)
            * mol_percent[["Na2O", "K2O", "H2O"]].sum(axis=axis)
        )
        + parameters["b13"]
        * (mol_percent["Al2O3"] * mol_percent[["Na2O", "K2O"]].sum(axis=axis))
    )

    return part_1 + part_2


def _calculate_C(melt_mol_fractions: Union[pd.DataFrame, pd.Series]) -> pd.Series:

    axis = [0, 1][isinstance(melt_mol_fractions, pd.DataFrame)]

    mol_percent = check_components(melt_mol_fractions, components=components) * 100

    part_1 = (
        parameters["c1"] * mol_percent["SiO2"]
        + parameters["c2"] * mol_percent[["TiO2", "Al2O3"]].sum(axis=axis)
        + parameters["c3"] * mol_percent[["FeO", "MgO", "MnO"]].sum(axis=axis)
        + parameters["c4"] * mol_percent["CaO"]
        + parameters["c5"] * mol_percent[["Na2O", "K2O"]].sum(axis=axis)
        + parameters["c6"] * np.log(1 + mol_percent[["H2O", "F2"]].sum(axis=axis))
    )

    part_2 = parameters["c11"] * (
        (
            mol_percent[["Al2O3", "FeO", "MgO", "MnO", "CaO"]].sum(axis=axis)
            - mol_percent["P2O5"]
        )
        * mol_percent[["Na2O", "K2O", "H2O", "F2"]].sum(axis=axis)
    )

    return part_1 + part_2


def calculate_viscosity(
    melt_mol_fractions: Union[pd.DataFrame, pd.Series], T_K: Union[pd.Series, float]
) -> pd.Series:
    """
    return log(viscosity) in log(Pa.s)
    """

    B = _calculate_B(melt_mol_fractions=melt_mol_fractions)
    C = _calculate_C(melt_mol_fractions=melt_mol_fractions)

    return parameters["A"] + B / (T_K - C)
