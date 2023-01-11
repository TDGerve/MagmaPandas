import numpy as np

from MagmaPandas import configuration
from MagmaPandas import Fe_redox
from MagmaPandas.fO2 import fO2_QFM

from ...Kd_baseclass import Kd_model


class FeMg_blundy(Kd_model):
    @classmethod
    def _get_Fe3Fe2Total(Melt_mol_fractions, T_K, P_bar=1, **kwargs):
        Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
        dQFM = kwargs.get("dQFM", configuration.dQFM)

        fO2 = fO2_QFM(dQFM, T_K, P_bar)
        Fe3Fe2 = Fe3Fe2_model(Melt_mol_fractions, T_K, fO2)

        return Fe3Fe2 / (1 + Fe3Fe2)

    @classmethod
    def calculate_Kd(
        cls, Melt_mol_fractions, forsterite, T_K, P_bar=1, *args, **kwargs
    ):
        """
        Blundy et al., 2020, equation 8
        """

        Fe3FeTotal = cls._get_Fe3Fe2_liquid(Melt_mol_fractions, T_K, P_bar)

        return 0.3642 * (1 - Fe3FeTotal) * np.exp(312.7 * (1 - 2 * forsterite) / T_K)
