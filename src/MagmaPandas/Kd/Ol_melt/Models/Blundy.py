import numpy as np
import pandas as pd
from MagmaPandas import Fe_redox, configuration
from MagmaPandas.fO2 import fO2_QFM
from MagmaPandas.Kd.Kd_baseclass import Kd_model


class FeMg_blundy(Kd_model):
    """
    Calculate equilibrium Fe-Mg partition coefficients between olivine and melt according to equation 8 from Blundy (2020)\ [11]_.
    """

    errors = pd.Series({6: 0.019, 9: 0.04, 100: 0.063})

    @classmethod
    def _get_Fe3FeTotal(cls, Melt_mol_fractions, T_K, P_bar=1, **kwargs):
        Fe3Fe2_model = Fe_redox.Fe3Fe2_borisov
        dQFM = kwargs.get("dQFM", configuration.dQFM)

        fO2 = fO2_QFM(dQFM, T_K, P_bar)
        Fe3Fe2 = Fe3Fe2_model.calculate_Fe3Fe2(Melt_mol_fractions, T_K, fO2)

        return Fe3Fe2 / (1 + Fe3Fe2)

    @classmethod
    def calculate_Kd(
        cls, Melt_mol_fractions: pd.DataFrame, forsterite, T_K, P_bar=1, *args, **kwargs
    ) -> float | pd.Series:
        """
        calculate equilibrium Kds for given melt compositions.

        Parameters
        ----------
        Melt_mol_fractions  : pandas Dataframe
            melt composition in oxide mol fractions
        forsterite  : float, array-like
            initial olivine forsterite content. Forsterite values are iteratively adjusted and initial values are not necessarily in Fe-Mg equilibrium with melts.
        T_K : float, array-like
            temperatures in Kelvin
        P_bar   : float array-like
            pressures in bar

        Returns
        -------
        Kds : float, array-like
        """

        Fe3FeTotal = cls._get_Fe3FeTotal(Melt_mol_fractions, T_K, P_bar)

        return 0.3642 * (1 - Fe3FeTotal) * np.exp(312.7 * (1 - 2 * forsterite) / T_K)

    @classmethod
    def get_error(cls, melt_composition: pd.DataFrame) -> float | pd.Series:
        """
        Calculate one standard deviation errors on Kds

        Parameters
        ----------
        melt_composition    : pandas Dataframe
            melt composition in oxide wt. %

        Returns
        -------
        error:  float, array-like
            one standard deviation error
        """

        axis = [0, 1][isinstance(melt_composition, pd.DataFrame)]
        alkalis = melt_composition[["Na2O", "K2O"]].sum(axis=axis)

        if isinstance(alkalis, (int, float)):

            composition_filter = np.less(alkalis, cls.errors.index.values)
            idx = cls.errors.index.values[composition_filter][0]

            return cls.errors[idx]

        errors = pd.Series(index=melt_composition.index)

        for i, a in alkalis.items():
            composition_filter = np.less(a, cls.errors.index.values)
            idx = cls.errors.index.values[composition_filter][0]
            errors.loc[i] = cls.errors[idx]

        return errors

    @classmethod
    def get_offset_parameters(cls, n=1):
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, melt_composition, offset_parameters, *args, **kwargs):

        error = cls.get_error(melt_composition=melt_composition)

        return offset_parameters * error
