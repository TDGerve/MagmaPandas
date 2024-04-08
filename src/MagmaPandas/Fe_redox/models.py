import inspect
import math
import sys
import warnings as w

import numpy as np
import pandas as pd

from MagmaPandas.Fe_redox.Fe_redox_baseclass import Fe3Fe2_model
from MagmaPandas.parse_io import check_components


def _is_Fe3Fe2_model(cls):
    """
    check if class inherits from Kd_model
    """
    try:
        return isinstance(cls(), Fe3Fe2_model)
    except TypeError:
        return False


# def check_components(melt_mol_fractions, components):
#     moles = melt_mol_fractions.copy()

#     if isinstance(moles, pd.DataFrame):
#         missing_oxides = set(components).difference(moles.columns)
#     elif isinstance(moles, pd.Series):
#         missing_oxides = set(components).difference(moles.index)

#     if len(missing_oxides) == 0:
#         return moles

#     for oxide in missing_oxides:
#         moles[oxide] = 0.0

#     w.warn(
#         f"{', '.join(str(i) for i in missing_oxides)} missing in composition and set to 0."
#     )

#     return moles.recalculate()


class borisov(Fe3Fe2_model):
    """
    Calculate melt |Fe3Fe2| ratios according to equation 4 from Borisov et al. (2018)\ [1]_.
    """

    components = ["SiO2", "TiO2", "MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5"]

    log_Fe3Fe2_error = 0.079  # standard deviation on residuals

    Fe3Fe2_error = 0.07 * 2  # From Putirka (2016)

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

        moles = check_components(Melt_mol_fractions, cls.components)

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
    def get_error(cls, *args, **kwargs):
        """
        Returns one standard deviation error on |Fe3Fe2| ratios, based on Putirka (2016).


        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """

        return cls.log_Fe3Fe2_error  # cls.log_Fe3Fe2_error * Fe3Fe2

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs):
        return offset_parameters * cls.get_error()


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

    Fe3Fe2_error = 0.07 * 2  # From Putirka (2016)

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

        moles = check_components(Melt_mol_fractions, cls.components)

        P_Pa = P_bar * 1e-5

        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        sumComponents = moles[cls.components].mul(cls.dCoefficients).sum(axis=axis)

        part1 = cls.a * np.log(fO2) + cls.b / T_K + cls.c + sumComponents
        part2 = cls.e * (1 - cls.T0 / T_K - np.log(T_K / cls.T0))
        part3 = (
            cls.f * P_Pa / T_K
            + cls.g * ((T_K - cls.T0) * P_Pa) / T_K
            + cls.h * P_Pa**2 / T_K
        )

        return 2 * np.exp(part1 + part2 + part3)

    @classmethod
    def get_error(cls, *args, **kwargs):
        """
        Returns one standard deviation error on |Fe3Fe2| ratios, based on Putirka (2016).

        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """

        # melt = melt_composition.FeO_Fe2O3_calc(Fe3Fe2, total_Fe="FeO", inplace=False)

        # Fe2O3_contribution = cls.FeO_error / 2 / melt["Fe2O3"]
        # FeO_contribution = cls.FeO_error / melt["FeO"]

        return cls.Fe3Fe2_error  # (Fe2O3_contribution + FeO_contribution) * Fe3Fe2

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, melt_composition, Fe3Fe2, offset_parameters, *args, **kwargs):
        # FeO_offset = offset_parameters * cls.FeO_error
        # Fe2O3_offset = (
        #     FeO_offset * 2
        # )  # Are FeO and Fe2O3 errors really correlated this way??

        # melt = melt_composition.FeO_Fe2O3_calc(Fe3Fe2, total_Fe="FeO", inplace=False)

        # Fe2O3_contribution = Fe2O3_offset / 2 / melt["Fe2O3"]
        # FeO_contribution = FeO_offset / melt["FeO"]

        # return (Fe2O3_contribution + FeO_contribution) * Fe3Fe2

        return cls.get_error() * offset_parameters


class jayasuriya(Fe3Fe2_model):
    """
    Jayasuriya et al. (2004) equation 12.
    """

    components = ["MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5", "FeO"]

    # Parameters from table 7
    a = 0.1967
    b = 12420
    c = 7.054

    dCoefficients = pd.Series(
        {
            "MgO": -0.487,
            "CaO": 2.201,
            "Na2O": 6.610,
            "K2O": 8.214,
            "Al2O3": -3.781,
            "P2O5": -62.79,
            "FeO": 1.377,
        },
        name="dCoeff",
    )

    log_Fe3O2FeO_error = math.sqrt(
        3.1 * (289 - 10) / 289
    )  # reduced Chi-squared to RMSE

    Fe3Fe2_error = 0.07 * 2  # from Putirka (2016)

    @classmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(Melt_mol_fractions, cls.components)

        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        sumComponents = moles[cls.components].mul(cls.dCoefficients).sum(axis=axis)

        return 2 * np.exp(cls.a * np.log(fO2) + 12420 / T_K - cls.c + sumComponents)

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:

        return cls.Fe3Fe2_error

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        return offset_parameters * cls.Fe3Fe2_error

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)


class putirka2016_6b(Fe3Fe2_model):
    """
    Putirka (2016) equation 6b.
    """

    components = ["Na2O", "K2O", "Al2O3", "SiO2", "CaO"]

    a = 6.53
    b = 10813.8
    c = 0.19
    d = 12.4
    e = 3.44
    f = 4.15

    Fe3Fe2_error = 0.07 * 2

    @classmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(Melt_mol_fractions, cls.components)

        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        part1 = -cls.a + cls.b / T_K
        part2 = cls.c * np.log(fO2) + 12.4 * moles[["Na2O", "K2O"]].sum(axis=axis)
        part3 = (
            -cls.e * (moles["Al2O3"] / moles[["Al2O3", "SiO2"]].sum(axis=axis))
            + cls.f * moles["CaO"]
        )

        return 2 * np.exp(part1 + part2 + part3)

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:

        return cls.Fe3Fe2_error

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:

        return cls.get_error() * offset_parameters

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)


class putirka2016_6c(Fe3Fe2_model):
    """
    Putirka (2016) equation 6c.
    """

    components = [
        "Al2O3",
        "Na2O",
        "K2O",
        "CaO",
        "MgO",
        "SiO2",
        "TiO2",
        "Cr2O3",
        "FeO",
        "MnO",
        "P2O5",
    ]

    a = 6.75
    b = 10634.9
    c = 0.195
    d = 7.9
    e = 4.6
    f = 0.54
    g = 53.4
    h = 1.07

    Fe3Fe2_error = (
        0.07 * 2
    )  # On page 823 Putirka implies eqs. 6b and 6c have similar errors: 'Two new calibrations, though (using DS1 and DS3) reduce error by 20–30%'

    @classmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(
            melt_mol_fractions=Melt_mol_fractions, components=cls.components
        )
        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        NBO_T = cls._NBO_T(cations=moles.cations)

        part_1 = -cls.a + cls.b / T_K + cls.c * np.log(fO2)
        part_2 = cls.d * moles[["Na2O", "K2O"]].sum(axis=axis) - cls.e * moles["MgO"]
        part_3 = (
            cls.f * (moles["MgO"] / moles[["MgO", "FeO"]].sum(axis=axis))
            - cls.g * moles["P2O5"]
            + cls.h * NBO_T
        )

        return 2 * np.exp(part_1 + part_2 + part_3)

    @staticmethod
    def _NBO_T(cations):

        axis = [0, 1][isinstance(cations, pd.DataFrame)]

        if axis == 1:
            Al_IV = min(
                cations.apply(
                    lambda row: min(
                        row["Al"],
                        row[["Na", "K"]].sum() + 2 * row[["Ca", "Mg"]].sum(),
                    ),
                    axis=1,
                )
            )
        else:
            Al_IV = min(
                cations["Al"], cations[["Na", "K"]].sum() + 2 * cations[["Ca", "Mg"]]
            )

        tetrahedral = cations[["Si", "Ti"]].sum(axis=axis) + Al_IV
        O = (
            2 * cations[["Si", "Ti"]].sum(axis=axis)
            + 1.5 * cations[["Al", "Cr"]].sum(axis=axis)
            + cations[["Fe", "Mn", "Mg", "Ca"]].sum(axis=axis)
            + 0.5 * cations[["Na", "K"]].sum(axis=axis)
            + 2.5 * cations["P"]
        )
        NBO = 2 * O - 4 * tetrahedral
        return NBO / tetrahedral

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:

        return cls.Fe3Fe2_error

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:

        return cls.get_error() * offset_parameters

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)


_clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
# Collect all Fe3Fe2_models in a dictionary.
Fe3Fe2_models = {cls[0]: cls[1] for cls in _clsmembers if _is_Fe3Fe2_model(cls[1])}
