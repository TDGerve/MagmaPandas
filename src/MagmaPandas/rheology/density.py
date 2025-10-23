import elementMass as e
import pandas as pd

from MagmaPandas.core.magma_protocol import Magma
from MagmaPandas.parse_io import check_components

"""
based on :cite:t:`Iacovino2019`
"""

T_reference = pd.Series(
    {
        "SiO2": 1773,
        "TiO2": 1773,
        "Al2O3": 1773,
        "Fe2O3": 1723,
        "FeO": 1723,
        "MgO": 1773,
        "CaO": 1773,
        "Na2O": 1773,
        "K2O": 1773,
        "H2O": 1273,
    }
)

molar_volumes = pd.Series(
    {
        "SiO2": 26.86,
        "TiO2": 28.32,
        "Al2O3": 37.42,
        "Fe2O3": 41.50,
        "FeO": 12.68,
        "MgO": 12.02,
        "CaO": 16.90,
        "Na2O": 29.65,
        "K2O": 47.28,
        "H2O": 22.9,
    }
)

dVdT = pd.Series(
    {
        "SiO2": 0.0,
        "TiO2": 0.00724,
        "Al2O3": 0.00262,
        "Fe2O3": 0.0,
        "FeO": 0.00369,
        "MgO": 0.00327,
        "CaO": 0.00374,
        "Na2O": 0.00768,
        "K2O": 0.01208,
        "H2O": 0.0095,
    }
)

dVdP = pd.Series(
    {
        "SiO2": -1.89e-4,
        "TiO2": -2.31e-4,
        "Al2O3": -2.26e-4,
        "Fe2O3": -2.5e-4,
        "FeO": -4.5e-5,
        "MgO": 2.7e-5,
        "CaO": 3.4e-5,
        "Na2O": -2.4e-4,
        "K2O": -6.75e-4,
        "H2O": -3.2e-4,
    }
)


def calculate_density(
    melt_wt_percent: Magma, T_K: float | pd.Series, P_bar: float | pd.Series
) -> pd.Series:
    """
    Calculate silicate liquid densities according to the model from :cite:t:`Iacovino2019`

    The model uses thermodynamic data from :cite:t:`Lange1987`, :cite:t:`Lange1997`, :cite:t:`kress_compressibility_1991` and :cite:t:`Ochs1999`

    Parameters
    ----------
    composition : Magma
        melt composition in oxide wt. %
    T_K : float, pandas Series
        temperatures in Kelvin
    P_bar : float, pandas Series
        pressures in bar

    Returns
    -------
    densities : pandas Series
        densities in kg/m\ :sup:`3`
    """
    axis = [0, 1][isinstance(melt_wt_percent, pd.DataFrame)]

    mol_fractions = check_components(melt_wt_percent, molar_volumes.index).moles()[
        molar_volumes.index
    ]

    oxide_masses = e.compound_weights(mol_fractions.columns)

    # equation 1 numerator
    mass_1_mol = mol_fractions.mul(oxide_masses, axis=1).sum(axis=axis)

    T_contribution = _calculate_T_contribution(T_K, index=mol_fractions.index)
    P_contribution = _calculate_P_contribution(P_bar, index=mol_fractions.index)

    V_liquid = (T_contribution + P_contribution).add(
        molar_volumes, axis=1
    ) * mol_fractions

    density = mass_1_mol / V_liquid.sum(axis=axis) * 1e3  # kg.m-3

    return density


def _calculate_temperature(
    melt_wt_percent: Magma, density: pd.Series, P_bar: pd.Series
) -> pd.Series:

    mol_fractions = check_components(
        melt_wt_percent, components=molar_volumes.index
    ).moles()[molar_volumes.index]

    oxide_masses = e.compound_weights(mol_fractions.columns)

    # eauation 1 numerator
    mass_1_mole = mol_fractions.mul(oxide_masses, axis=1).sum(axis=1)

    P_contribution = _calculate_P_contribution(P_bar, index=mol_fractions.index)

    V_liquid = mass_1_mole * 1e3 / density

    V_P_contribution = (P_contribution * mol_fractions).sum(axis=1)
    V_reference = mol_fractions.mul(molar_volumes, axis=1).sum(axis=1)
    temperature_contribution = V_liquid - V_P_contribution - V_reference

    T_contribution_1500K = mol_fractions.mul(
        _calculate_T_contribution(
            pd.Series(1500, index=mol_fractions.index), mol_fractions.index
        ),
        axis=1,
    ).sum(axis=1)
    V_missing = temperature_contribution - T_contribution_1500K
    V_per_degree = mol_fractions.mul(dVdT, axis=1).sum(axis=1)

    return 1500 + (V_missing / V_per_degree)


def _calculate_P_contribution(P_bar: pd.Series, index):
    P_contribution = pd.DataFrame(columns=dVdP.index, index=index)
    P_contribution.loc[:, :] = dVdP.values[None, :]
    P_contribution = P_contribution.mul(P_bar - 1, axis=0)

    return P_contribution


def _calculate_T_contribution(T_K: pd.Series, index):
    if isinstance(T_K, (int, float)):
        T_K = pd.Series(T_K, index=index)

    T_contribution = pd.DataFrame(columns=(T_reference.index), index=index)
    T_contribution.loc[:, :] = T_K.values[:, None]
    T_contribution = T_contribution.sub(T_reference, axis=1)
    T_contribution = T_contribution.mul(dVdT, axis=1)

    return T_contribution
