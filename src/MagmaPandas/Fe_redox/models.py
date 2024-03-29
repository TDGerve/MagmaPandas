import inspect
import sys
import warnings as w

import numpy as np
import pandas as pd

from MagmaPandas.Fe_redox.Fe_redox_baseclass import Fe3Fe2_model


def _is_Fe3Fe2_model(cls):
    """
    check if class inherits from Kd_model
    """
    try:
        return isinstance(cls(), Fe3Fe2_model)
    except TypeError:
        return False


class borisov(Fe3Fe2_model):
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


class kressCarmichael(Fe3Fe2_model):
    """
    Calculate melt |Fe3Fe2| ratios according to equation 7 from Kress and Carmichael (1991)\ [2]_.
    """

    components = ["Al2O3", "FeO", "CaO", "Na2O", "K2O"]

    # Parameters from table 7
    a = 0.196
    b = 1.1492e4
    c = -6.675

    dCoefficients = pd.Series(
        {"Al2O3": -2.243, "FeO": -1.828, "CaO": 3.201, "Na2O": 5.854, "K2O": 6.215},
        name="dCoeff",
    )

    e = -3.36
    f = -7.01e-7
    g = -1.54e-10
    h = 3.85e-17
    T0 = 1673

    FeO_error = 0.21
    Fe2O3_error = 0.42

    @classmethod
    def calculate_Fe3Fe2(cls, Melt_mol_fractions, T_K, fO2, P_bar):
        """
        Calculate melt |Fe3Fe2| ratios\.

        Parameters
        ----------
        mol_fractions   :  :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            melt composition in oxide mol fractions
        T_K        :     float, array-like
            temperature in Kelvin
        fO2     :      float, array-like
            Oxygen fugacity
        P_Pa      :      float, array-like
            Pressure in Pascals


        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        P_Pa = P_bar * 1e-5

        LNfO2 = np.log(fO2)

        axis = [0, 1][isinstance(Melt_mol_fractions, pd.DataFrame)]

        sumComponents = (
            Melt_mol_fractions[cls.components].mul(cls.dCoefficients).sum(axis=axis)
        )

        part1 = cls.a * LNfO2 + cls.b / T_K + cls.c + sumComponents
        part2 = cls.e * (1 - cls.T0 / T_K - np.log(T_K / cls.T0))
        part3 = (
            cls.f * P_Pa / T_K
            + cls.g * ((T_K - cls.T0) * P_Pa) / T_K
            + cls.h * P_Pa**2 / T_K
        )

        return 2 * np.exp(part1 + part2 + part3)

    @classmethod
    def get_error(cls, melt_composition, Fe3Fe2):
        """
        Calculate one standard deviation errors on |Fe3Fe2| ratios, based on reported calibration errors of 0.21 and 0.42 wt.% on FeO and |Fe2O3|.

        Parameters
        ----------
        melt_composition : :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame`
            melt composition in wt. % oxides
        Fe3Fe2  :   float, array-like
            melt |Fe3Fe2| ratios

        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """

        melt = melt_composition.FeO_Fe2O3_calc(Fe3Fe2, total_Fe="FeO", inplace=False)

        Fe2O3_contribution = cls.FeO_error / 2 / melt["Fe2O3"]
        FeO_contribution = cls.FeO_error / melt["FeO"]

        return (Fe2O3_contribution + FeO_contribution) * Fe3Fe2

    @classmethod
    def get_offset_parameters(cls, n: int = 1):
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, melt_composition, Fe3Fe2, offset_parameters, *args, **kwargs):
        FeO_offset = offset_parameters * cls.FeO_error
        Fe2O3_offset = (
            FeO_offset * 2
        )  # Are FeO and Fe2O3 errors really correlated this way??

        melt = melt_composition.FeO_Fe2O3_calc(Fe3Fe2, total_Fe="FeO", inplace=False)

        Fe2O3_contribution = Fe2O3_offset / 2 / melt["Fe2O3"]
        FeO_contribution = FeO_offset / melt["FeO"]

        return (Fe2O3_contribution + FeO_contribution) * Fe3Fe2


_clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
# Collect all Fe3Fe2_models in a dictionary.
Fe3Fe2_models = {cls[0]: cls[1] for cls in _clsmembers if _is_Fe3Fe2_model(cls[1])}
