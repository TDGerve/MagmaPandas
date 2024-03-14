from functools import partial

import pandas as pd

from MagmaPandas import Fe_redox
from MagmaPandas.configuration import configuration
from MagmaPandas.fO2 import fO2_QFM
from MagmaPandas.Kd.Ol_melt import Kds_iterate as it
from MagmaPandas.Kd.Ol_melt import Kds_vectorised as vec
from MagmaPandas.MagmaFrames import MagmaFrame


def observed_FeMg_Kd(
    melt: MagmaFrame,
    forsterite: pd.Series,
    P_bar,
    T_K: None | float | pd.Series = None,
    **kwargs,
) -> pd.Series:
    """
    Calulate observed Kds based on given melt and olivine compositions.

    Parameters
    ----------
    melt : MagmaFrame
        melt composition in oxide wt. %
    forsterite : float, pandas Series
        olivine forsterite content
    P_bar   : float, pandas Series
        pressures in bar
    T_K : None, float, pandas Series
        temperatures in Kelvin. If set to None, temperatures are calculated according to the melt thermometer set in the global configuration.

    Returns
    -------
    Kds : pandas Series
    """
    if isinstance(P_bar, (int, float)):
        P_bar = pd.Series(P_bar, index=melt.index)
    if isinstance(forsterite, (int, float)):
        forsterite = pd.Series(forsterite, index=melt.index)

    forsterite = forsterite[melt.index]

    melt_x_moles = melt.moles

    Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
    dQFM = kwargs.get("dQFM", configuration.dQFM)

    if T_K is None:
        T_K = melt_x_moles.convert_moles_wtPercent().temperature(P_bar)
    fO2 = fO2_QFM(dQFM, T_K, P_bar)
    Fe3Fe2 = Fe3Fe2_model(melt_x_moles, T_K, fO2)

    melt_x_moles = melt_x_moles.normalise()
    Fe2_FeTotal = 1 / (1 + Fe3Fe2)
    melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
    olivine_MgFe = forsterite / (1 - forsterite)
    Kd_observed = melt_MgFe / olivine_MgFe

    return Kd_observed


def calculate_FeMg_Kd(
    Melt_mol_fractions: pd.Series | pd.DataFrame,
    forsterite_initial: float | pd.Series,
    T_K: float | pd.Series,
    Fe3Fe2: float | pd.Series,
    *args,
    **kwargs,
) -> pd.Series:
    """
    Calculate equilibrium Kds for given melt compositions.

    Kd models are set in the global configuration

    Parameters
    ----------
    Melt_mol_fractions : pandas Series, pandas Dataframe
        melt composition in oxide mol fractions
    forsterite_initial  : float, pandas Series
        initial forsterite contents. Forsterite values are iteratively adjusted and initial values are not necessarily in Fe-Mg equilibrium with melts.
    T_K : float, pandas Series
        temperatures in Kelvin
    Fe3Fe2  : float, pandas Series
        melt |Fe3Fe2| ratios

    Returns
    -------
    Kds : pandas Series
    """

    if isinstance(Melt_mol_fractions, pd.Series):
        Kd_func = it.calculate_Kd
    elif isinstance(Melt_mol_fractions, pd.DataFrame):
        Kd_func = vec.calculate_Kd

    return Kd_func(
        Melt_mol_fractions=Melt_mol_fractions,
        forsterite_initial=forsterite_initial,
        T_K=T_K,
        Fe3Fe2=Fe3Fe2,
        *args,
        **kwargs,
    )
