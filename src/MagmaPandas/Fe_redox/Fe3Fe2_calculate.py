from MagmaPandas import fO2 as f
from MagmaPandas.configuration import configuration
from MagmaPandas.Fe_redox.Fe3Fe2_models import Fe3Fe2_models_dict


def calculate_Fe3Fe2(mol_fractions, T_K, P_bar, fO2=None, **kwargs):
    """
    Calculate melt |Fe3Fe2| with the configured *f*\ O2 buffer and |Fe3Fe2| model.

    Model choices are set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

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
    Fe_model = Fe3Fe2_models_dict[Fe_model_name]  # getattr(Fe_redox, Fe_model_name)

    if fO2 is None:
        fO2 = f.calculate_fO2(T_K=T_K, P_bar=P_bar, **kwargs)

    return Fe_model.calculate_Fe3Fe2(
        melt_mol_fractions=mol_fractions, T_K=T_K, P_bar=P_bar, fO2=fO2
    )
