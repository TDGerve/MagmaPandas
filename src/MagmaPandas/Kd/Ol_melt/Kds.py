from functools import partial

import pandas as pd

from MagmaPandas import Fe_redox
from MagmaPandas.configuration import configuration
from MagmaPandas.fO2 import fO2_QFM

from . import Kds_iterate as it
from . import Kds_vectorised as vec


def observed_FeMg_Kd(melt: pd.DataFrame, forsterite: pd.Series, P_bar, **kwargs):
    """
    Calulate Fe-Mg exchange coefficients between olivine (ol) and melt (m) as:

    [Fe(ol) / Fe(m)] * [Mg(m) / Mg(ol)]
    """
    if isinstance(P_bar, (int, float)):
        P_bar = pd.Series(P_bar, index=melt.index)
    if isinstance(forsterite, (int, float)):
        forsterite = pd.Series(forsterite, index=melt.index)

    forsterite = forsterite[melt.index]

    melt_x_moles = melt.moles

    Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
    dQFM = kwargs.get("dQFM", configuration.dQFM)

    T_K = melt_x_moles.convert_moles_wtPercent.temperature(P_bar)
    fO2 = fO2_QFM(dQFM, T_K, P_bar)
    Fe3Fe2 = Fe3Fe2_model(melt_x_moles, T_K, fO2)

    melt_x_moles = melt_x_moles.normalise()
    Fe2_FeTotal = 1 / (1 + Fe3Fe2)
    melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
    olivine_MgFe = forsterite / (1 - forsterite)
    Kd_observed = melt_MgFe / olivine_MgFe

    return Kd_observed


def calculate_FeMg_Kd(
    Melt_mol_fractions, forsterite_initial, T_K, Fe3Fe2, *args, **kwargs
):

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
