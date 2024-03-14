import warnings as w

import numpy as np
import pandas as pd

from MagmaPandas.Fe_redox.Fe_redox_baseclass import Fe3Fe2_model


class Fe3Fe2_borisov(Fe3Fe2_model):
    """
    Calculate melt |Fe3Fe2| ratios according to equation 4 from Borisov et al. (2018)\ [1]_.
    """

    log_Fe3Fe2_error = 0.079

    @classmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions: pd.DataFrame, T_K, fO2, *args, **kwargs
    ):
        """
        Calculate melt |Fe3Fe2| ratios.

        Parameters
        ----------
        mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """

        moles = Melt_mol_fractions.copy()

        oxides = ["SiO2", "TiO2", "MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5"]

        if isinstance(moles, pd.DataFrame):
            missing_oxides = set(oxides).difference(moles.columns)
        elif isinstance(moles, pd.Series):
            missing_oxides = set(oxides).difference(moles.index)

        if len(missing_oxides) > 0:

            for oxide in missing_oxides:
                moles[oxide] = 0.0

            w.warn(
                f"{', '.join(str(i) for i in missing_oxides)} missing in composition and set to 0."
            )

        part1 = (
            0.207 * np.log10(fO2)
            + 4633.3 / T_K
            - 0.445 * moles["SiO2"]
            - 0.900 * moles["TiO2"]
            + 1.532 * moles["MgO"]
        )
        part2 = (
            0.314 * moles["CaO"]
            + 2.030 * moles["Na2O"]
            + 3.355 * moles["K2O"]
            - 4.851 * moles["P2O5"]
        )
        part3 = (
            -3.081 * moles["SiO2"] * moles["Al2O3"]
            - 4.370 * moles["SiO2"] * moles["MgO"]
            - 1.852
        )

        return 10 ** (part1 + part2 + part3)

    @classmethod
    def get_error(cls, Fe3Fe2, *args, **kwargs):
        """
        Calculate one standard deviation errors on |Fe3Fe2| ratios, based on a reported 0.079 calibration error on log(|Fe3Fe2|)

        Parameters
        ----------
        Fe3Fe2  :   float, array-like
            melt |Fe3Fe2| ratios

        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """

        return cls.log_Fe3Fe2_error * Fe3Fe2

    @classmethod
    def get_offset_parameters(cls, n: int = 1):
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, Fe3Fe2, offset_parameters, *args, **kwargs):
        return offset_parameters * cls.log_Fe3Fe2_error * Fe3Fe2
