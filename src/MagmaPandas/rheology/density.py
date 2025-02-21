import elementMass as e
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma

"""
based on Iacovino and Till (2019)\ [23]_
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
        "SiO2": -0.000189,
        "TiO2": -0.000231,
        "Al2O3": -0.000226,
        "Fe2O3": -0.000253,
        "FeO": -0.000045,
        "MgO": 0.000027,
        "CaO": 0.000034,
        "Na2O": -0.00024,
        "K2O": -0.000675,
        "H2O": -0.00032,
    }
)


def calculate_density(
    composition: Magma, T_K: float | pd.Series, P_bar: float | pd.Series
) -> pd.Series:
    """
    Calculate silicate liquid densities according to the model from Iacovino and Till (2019)\ [23]_

    The model uses thermodynamic data from Lange and Carmichael (1987)\ [24]_, Lange (1997)\ [25]_, Kress and Carmichael (1991)\ [2]_ and Ochs III and Lange (1999)\ [26]_

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
    mole_fractions = composition.moles()[molar_volumes.index]

    oxide_masses = e.compound_weights(mole_fractions.columns)

    # equation 1 numerator
    mass_1_mole = mole_fractions.mul(oxide_masses, axis=1).sum(axis=1)

    T_contribution = _calculate_T_contribution(T_K, index=mole_fractions.index)
    P_contribution = _calculate_P_contribution(P_bar, index=mole_fractions.index)

    V_liquid = (T_contribution + P_contribution).add(
        molar_volumes, axis=1
    ) * mole_fractions

    density = mass_1_mole / V_liquid.sum(axis=1) * 1e3  # kg.m-3

    return density


def _calculate_temperature(
    composition: Magma, density: pd.Series, P_bar: pd.Series
) -> pd.Series:
    mole_fractions = composition.moles()[molar_volumes.index]

    oxide_masses = e.compound_weights(mole_fractions.columns)

    # eauation 1 numerator
    mass_1_mole = mole_fractions.mul(oxide_masses, axis=1).sum(axis=1)

    P_contribution = _calculate_P_contribution(P_bar, index=mole_fractions.index)

    V_liquid = mass_1_mole * 1e3 / density

    V_P_contribution = (P_contribution * mole_fractions).sum(axis=1)
    V_reference = mole_fractions.mul(molar_volumes, axis=1).sum(axis=1)
    temperature_contribution = V_liquid - V_P_contribution - V_reference

    T_contribution_1500K = mole_fractions.mul(
        _calculate_T_contribution(
            pd.Series(1500, index=mole_fractions.index), mole_fractions.index
        ),
        axis=1,
    ).sum(axis=1)
    V_missing = temperature_contribution - T_contribution_1500K
    V_per_degree = mole_fractions.mul(dVdT, axis=1).sum(axis=1)

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
