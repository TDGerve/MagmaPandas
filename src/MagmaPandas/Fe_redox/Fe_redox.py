from MagmaPandas import fO2
from MagmaPandas.configuration import configuration
from MagmaPandas.Fe_redox.models import Fe3Fe2_models


def FeRedox_QFM(mol_fractions, T_K, P_bar, logshift, **kwargs):
    """
    Calculate melt |Fe3Fe2| ratios at the QFM *f*\ O2 buffer, with
    models from Borisov et al. (2018)\ [1]_ or Kress and Carmichael (1991)\ [2]_.

    Model choice is set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

    Parameters
    ----------
    mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
        Melt composition in oxide mol_fractions
    T_K  :       float, pd.Series-like
        temperature in Kelvin
    P_bar   :     float, pd.Series-like
        Pressure in bars
    logshift  :  int, pd.Series-like
        log units shift of QFM buffer
    model     :  string
        'kressCarmichael' or 'borisov'

    Returns
    -------
    melt Fe\ :sup:`3+`\ /Fe\ :sup:`2+` ratio
    """

    Fe_model_name = kwargs.get("Fe_model", configuration.Fe3Fe2_model)
    Fe_model = Fe3Fe2_models[Fe_model_name]  # getattr(Fe_redox, Fe_model_name)

    fO2_bar = kwargs.get("fO2", fO2.fO2_QFM(logshift, T_K, P_bar))

    return Fe_model.calculate_Fe3Fe2(mol_fractions, T_K, fO2_bar, P_bar)
