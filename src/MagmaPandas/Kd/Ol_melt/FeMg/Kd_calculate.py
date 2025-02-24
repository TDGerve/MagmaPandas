import pandas as pd

from MagmaPandas.Fe_redox.Fe3Fe2_models import Fe3Fe2_models_dict
from MagmaPandas.fO2.fO2_calculate import calculate_fO2
from MagmaPandas.Kd.Ol_melt.FeMg.Kd_models import Kd_olmelt_FeMg_models_dict
from MagmaPandas.magma_protocol import Magma

# configuration needs to be imported last to avoid circular imports
from MagmaPandas.configuration import configuration  # isort: skip


def observed_FeMg_Kd(
    melt: Magma,
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

    melt_x_moles = melt.moles()

    Fe3Fe2_model = Fe3Fe2_models_dict[configuration.Fe3Fe2_model]
    dfO2 = kwargs.get("dfO2", configuration.dfO2)

    if T_K is None:
        T_K = melt_x_moles.temperature(P_bar)
    fO2 = calculate_fO2(T_K=T_K, P_bar=P_bar, dfO2=dfO2)
    Fe3Fe2 = Fe3Fe2_model.calculate_Fe3Fe2(
        melt_mol_fractions=melt_x_moles, T_K=T_K, fO2=fO2, P_bar=P_bar
    )

    melt_x_moles = melt_x_moles.normalise()
    Fe2_FeTotal = 1 / (1 + Fe3Fe2)
    melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
    olivine_MgFe = forsterite / (1 - forsterite)
    Kd_observed = melt_MgFe / olivine_MgFe

    return Kd_observed


def calculate_FeMg_Kd(
    melt_mol_fractions: pd.Series | pd.DataFrame,
    T_K: float | pd.Series,
    *args,
    **kwargs,
) -> pd.Series:
    """
    Calculate equilibrium Kds for given melt compositions.

    Kd models are set in the global configuration

    Parameters
    ----------
    melt_mol_fractions : pandas Series, pandas Dataframe
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

    Kd_model_name = kwargs.get("Kd_model", configuration.Kd_model)
    Kd_model = Kd_olmelt_FeMg_models_dict[Kd_model_name]

    if Kd_model_name == "blundy":
        dfO2 = kwargs.get("dfO2", configuration.dfO2)
        P_bar = kwargs["P_bar"]
        kwargs["fO2"] = calculate_fO2(T_K=T_K, P_bar=P_bar, dfO2=dfO2)

    return Kd_model.calculate_Kd(
        melt_mol_fractions=melt_mol_fractions, T_K=T_K, *args, **kwargs
    )

    # if isinstance(melt_mol_fractions, pd.Series):
    #     Kd_func = it.calculate_Kd
    # elif isinstance(melt_mol_fractions, pd.DataFrame):
    #     Kd_func = vec.calculate_Kd

    # return Kd_func(
    #     melt_mol_fractions=melt_mol_fractions,
    #     forsterite_initial=forsterite_initial,
    #     T_K=T_K,
    #     Fe3Fe2=Fe3Fe2,
    #     *args,
    #     **kwargs,
    # )
