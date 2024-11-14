import inspect
import sys
from functools import partial

import numpy as np
import pandas as pd
from scipy.constants import R  # J*K-1*mol-1

from MagmaPandas import Fe_redox, configuration
from MagmaPandas.fO2 import calculate_fO2
from MagmaPandas.Kd.Kd_baseclass import Kd_model
from MagmaPandas.Kd.Ol_melt.iterative import iterate_Kd_scalar, iterate_Kd_vectorized
from MagmaPandas.parse_io import check_components
from MagmaPandas.thermometers.data_parsing import _remove_elements, moles_per_oxygen


def _is_Kd_model(cls):
    """
    check if class inherits from Kd_model
    """
    try:
        return isinstance(cls(), Kd_model)
    except TypeError:
        return False


class toplis(Kd_model):
    """
    Calculate equilibrium Fe-Mg partition coefficients between olivine and melt according to equation 10 from Toplis (2005)\ [10]_.
    """

    error = 0.02

    @classmethod
    def _Phi(cls, molar_SiO2, molar_Na2O, molar_K2O):
        """
        Equation 12 from Toplis (2005) calculates a Phi parameter to correct SiO2 for alkali-bearing liquids.

        Parameters
        ----------
        molar_SiO2 : int or list-like

        molar_Na2O : int or list-like

        molar_K2O : int or list-like


        Returns
        -------
            int or list-like
        """

        try:
            if sum(high_SiO2 := np.array(molar_SiO2) > 60) > 1:
                # raise RuntimeError("SiO2 >60 mol% present")
                results = cls._Phi_lowSiO2(molar_SiO2, molar_Na2O, molar_K2O)
                results[high_SiO2] = cls._Phi_highSiO2(
                    molar_SiO2[high_SiO2], molar_Na2O[high_SiO2], molar_K2O[high_SiO2]
                )
                return results

        except TypeError:
            if molar_SiO2 > 60:
                # raise RuntimeError("SiO2 >60 mol%")
                return cls._Phi_highSiO2(molar_SiO2, molar_Na2O, molar_K2O)

        return cls._Phi_lowSiO2(molar_SiO2, molar_Na2O, molar_K2O)

    @classmethod
    def _Phi_highSiO2(cls, molar_SiO2, molar_Na2O, molar_K2O):
        """
        Phi for SiO2 > 60 mol%
        """

        return (11 - 5.5 * 100 / (100 - molar_SiO2)) * np.exp(
            -0.31 * (molar_Na2O + molar_K2O)
        )

    @classmethod
    def _Phi_lowSiO2(cls, molar_SiO2, molar_Na2O, molar_K2O):
        """
        Phi for SiO2 < 60 mol%
        """

        return (0.46 * (100 / (100 - molar_SiO2)) - 0.93) * (molar_Na2O + molar_K2O) + (
            -5.33 * (100 / (100 - molar_SiO2)) + 9.69
        )

    @classmethod
    def _SiO2_A(cls, melt_mol_fractions):
        """returns adjusted SiO2 for Toplis (2005) Fe-Mg Kd calculations

        Equations 11 and 14 calculate adjusted molar SiO2 by correcting for akalis and water


        Parameters
        ----------
        molar_SiO2 : int or list-like

        molar_Na2O : int or list-like

        molar_K2O : int or list-like

        Phi : int or list-like
            coefficient for alkali correction, needs to be calculated according to Toplis (eq 12, 2005)

        H2O : int or list-like, optional
            wt. %


        Returns
        -------
        int or list-like
        """

        # Calculate melt molar concentrations
        # Molar fractions normalised to 1
        melt_mol_fractions = melt_mol_fractions.fillna(0.0)
        molar_concentrations = melt_mol_fractions * 100

        molar_SiO2 = molar_concentrations["SiO2"]
        molar_Na2O = molar_concentrations["Na2O"]
        molar_K2O = molar_concentrations["K2O"]

        Phi = cls._Phi(molar_SiO2, molar_Na2O, molar_K2O)
        # Equation 11
        SiO2_A = molar_SiO2 + Phi * (molar_Na2O + molar_K2O)

        try:
            # For dataframes
            if "H2O" in molar_concentrations.columns:
                SiO2_A = SiO2_A + 0.8 * molar_concentrations["H2O"]  # equation 14
        except:
            # For series
            if "H2O" in molar_concentrations.index:
                SiO2_A = SiO2_A + 0.8 * molar_concentrations["H2O"]  # equation 14

        return SiO2_A

    @classmethod
    @np.errstate(invalid="raise")
    def _calculate_Kd(
        cls,
        melt_mol_fractions: pd.DataFrame,
        forsterite: float | pd.Series,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        *args,
        **kwargs,
    ) -> float | pd.Series:
        """
        Calculate Kds for given melt compositions and fixed forsterite content.

        Parameters
        ----------
        melt_mol_fractions : pandas Dataframe
            melt compositions in oxide mol fractions
        forsterite : float, array-like
            olivine forsterite contents.
        T_K : float, array-like
            temperatures in Kelvin
        P_bar: float, array-like
            pressures in bar

        Returns
        -------
        Kds : array-like
            Fe-Mg partition coefficients
        """

        SiO2_A = cls._SiO2_A(melt_mol_fractions)

        return np.exp(
            (-6766 / (R * T_K) - 7.34 / R)
            + np.log(0.036 * SiO2_A - 0.22)
            + (3000 * (1 - 2 * forsterite) / (R * T_K))
            + (0.035 * (P_bar - 1) / (R * T_K))
        )

    @classmethod
    def calculate_Kd(
        cls,
        melt_mol_fractions: pd.Series | pd.DataFrame,
        Fe3Fe2: float | pd.Series,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        forsterite_initial: float | pd.Series = 0.85,
        *args,
        **kwargs,
    ) -> float | pd.Series:
        """
        Calculate Kds for given melt compositions and equilibriutm forsterite content.

        Parameters
        ----------
        melt_mol_fractions : pandas Dataframe
            melt compositions in oxide mol fractions
        forsterite_initial : float, array-like
            initial olivine forsterite contents. Forsterite values are iteratively adjusted until equilibrium with the melt is reached.
        Fe3Fe2 : float, array-like
            melt Fe3+/Fe2+ ratios
        T_K : float, array-like
            temperatures in Kelvin
        P_bar: float, array-like
            pressures in bar

        Returns
        -------
        Kds : array-like
            Fe-Mg partition coefficients
        """

        if isinstance(melt_mol_fractions, pd.Series):
            Kd_func = iterate_Kd_scalar
        elif isinstance(melt_mol_fractions, pd.DataFrame):
            Kd_func = iterate_Kd_vectorized

        return Kd_func(
            melt_mol_fractions=melt_mol_fractions,
            forsterite_initial=forsterite_initial,
            Fe3Fe2=Fe3Fe2,
            T_K=T_K,
            P_bar=P_bar,
            Kd_model=cls._calculate_Kd,
            *args,
            **kwargs,
        )

    @classmethod
    def get_error(cls) -> float:
        """
        Calculate one standard deviation errors on Kds

        Returns
        -------
        error:  float
            one standard deviation error
        """

        return cls.error

    @classmethod
    def get_offset_parameters(cls, n=1):
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs):
        error = cls.get_error()
        return offset_parameters * error


class blundy(Kd_model):
    """
    Calculate equilibrium Fe-Mg partition coefficients between olivine and melt according to equation 8 from Blundy (2020)\ [11]_.
    """

    errors = pd.Series({6: 0.019, 9: 0.04, 100: 0.063})

    @classmethod
    def _get_Fe3FeTotal(cls, Fe3Fe2, **kwargs):

        # Fe3Fe2_model = Fe_redox.borisov
        # dfO2 = kwargs.get("dfO2", configuration.dfO2)

        # fO2 = calculate_fO2(T_K=T_K, P_bar=P_bar, dfO2=dfO2)
        # Fe3Fe2 = Fe3Fe2_model.calculate_Fe3Fe2(
        #     melt_mol_fractions=melt_mol_fractions, T_K=T_K, fO2=fO2, P_bar=P_bar
        # )

        return Fe3Fe2 / (1 + Fe3Fe2)

    @classmethod
    def _calculate_Kd(
        cls, forsterite, T_K, Fe3Fe2, *args, **kwargs
    ) -> float | pd.Series:
        """
        calculate equilibrium Kds for given melt compositions.

        Parameters
        ----------
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

        Fe3FeTotal = cls._get_Fe3FeTotal(Fe3Fe2=Fe3Fe2)

        Kd_Fe_total = (
            0.3642 * (1 - Fe3FeTotal) * np.exp((312.7 * (1 - 2 * forsterite)) / T_K)
        )

        Kd_Fe2 = Kd_Fe_total / (1 - Fe3FeTotal)

        return Kd_Fe2

    @classmethod
    def calculate_Kd(
        cls,
        melt_mol_fractions: pd.Series | pd.DataFrame,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        forsterite_initial: float | pd.Series = 0.85,
        *args,
        **kwargs,
    ) -> float | pd.Series:
        """
        Calculate Kds for given melt compositions and equilibriutm forsterite content.

        Parameters
        ----------
        melt_mol_fractions : pandas Dataframe
            melt compositions in oxide mol fractions
        forsterite_initial : float, array-like
            initial olivine forsterite contents. Forsterite values are iteratively adjusted until equilibrium with the melt is reached.
        T_K : float, array-like
            temperatures in Kelvin
        P_bar: float, array-like
            pressures in bar

        Returns
        -------
        Kds : array-like
            Fe-Mg partition coefficients
        """
        # Fe3Fe2 needs to be calculated with 'borisov' and should be removed from kwargs if present.
        _ = kwargs.pop("Fe3Fe2", None)

        if isinstance(melt_mol_fractions, pd.Series):
            Kd_func = iterate_Kd_scalar
        elif isinstance(melt_mol_fractions, pd.DataFrame):
            Kd_func = iterate_Kd_vectorized

        dfO2 = kwargs.get("dfO2", configuration.dfO2)

        fO2 = calculate_fO2(T_K=T_K, P_bar=P_bar, dfO2=dfO2)
        Fe3Fe2 = Fe_redox.borisov.calculate_Fe3Fe2(
            melt_mol_fractions=melt_mol_fractions, T_K=T_K, fO2=fO2, P_bar=P_bar
        )

        model = partial(cls._calculate_Kd, Fe3Fe2=Fe3Fe2)

        return Kd_func(
            melt_mol_fractions=melt_mol_fractions,
            forsterite_initial=forsterite_initial,
            T_K=T_K,
            P_bar=P_bar,
            Fe3Fe2=Fe3Fe2,
            Kd_model=model,
            *args,
            **kwargs,
        )

    @classmethod
    def get_error(
        cls, melt_composition: pd.DataFrame, *args, **kwargs
    ) -> float | pd.Series:
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
        alkalis = melt_composition.wt_pc()[["Na2O", "K2O"]].sum(axis=axis)

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


class putirka2016_8a(Kd_model):

    Kd_error = 4.4e-2

    @staticmethod
    def calculate_Kd(melt_mol_fractions, *args, **kwargs) -> float:
        """
        Calculate mineral-melt partition coefficients

        Returns
        -------
        float
            mineral-melt partition coefficients
        """

        if isinstance(melt_mol_fractions, (pd.Series, np.ndarray)):
            return 0.33
        elif isinstance(melt_mol_fractions, pd.DataFrame):
            return pd.Series(0.33, index=melt_mol_fractions.index)

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:
        """
        Return one standard deviation errors on partition coefficients.

        Returns
        -------
        float, array-like
            partition coefficient errors
        """
        return cls.Kd_error

    @staticmethod
    def get_offset_parameters(n: int, *args, **kwargs) -> float | np.ndarray:
        """
        Randomly sample a standard normal distribution *n* times.

        n   : int
            sample amount.
        """
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate random samples of partition coefficient errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """
        return cls.get_error() * offset_parameters


class putirka2016_8b(Kd_model):
    """
    For P > 1 GPa
    """

    a = 0.21
    b = 8e-3
    c = 2.5e-3
    d = -3.63e-4

    components = ["SiO2", "Na2O", "K2O"]

    Kd_error = 4.4e-2

    @classmethod
    def calculate_Kd(
        cls, melt_mol_fractions, P_bar, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate mineral-melt partition coefficients

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        P_bar :   float, array-like
            pressure in bar

        Returns
        -------
        float, array-like
            mineral-melt partition coefficients
        """
        wt_pc = melt_mol_fractions.wt_pc()
        wt_pc = check_components(wt_pc, cls.components)
        P_GPa = P_bar / 1e4
        axis = [0, 1][isinstance(wt_pc, pd.DataFrame)]

        return (
            cls.a
            + cls.b * P_GPa
            + cls.c * wt_pc["SiO2"]
            + cls.d * (wt_pc[["Na2O", "K2O"]].sum(axis=axis) ** 2)
        )

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:
        """
        Return one standard deviation errors on partition coefficients.

        Returns
        -------
        float, array-like
            partition coefficient errors
        """
        return cls.Kd_error

    @staticmethod
    def get_offset_parameters(n: int, *args, **kwargs) -> float | np.ndarray:
        """
        Randomly sample a standard normal distribution *n* times.

        n   : int
            sample amount.
        """
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate random samples of partition coefficient errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """
        return cls.get_error() * offset_parameters


class putirka2016_8c(Kd_model):
    """
    for P < 1 GPa
    """

    a = 0.25
    b = 1.8e-3
    c = -3.25e-4

    components = ["SiO2", "Na2O", "K2O"]

    Kd_error = 4e-2

    @classmethod
    def calculate_Kd(cls, melt_mol_fractions, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate mineral-melt partition coefficients

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions


        Returns
        -------
        float, array-like
            mineral-melt partition coefficients
        """
        wt_pc = melt_mol_fractions.wt_pc()
        wt_pc = check_components(wt_pc, cls.components)
        axis = [0, 1][isinstance(wt_pc, pd.DataFrame)]

        return (
            cls.a
            + cls.b * wt_pc["SiO2"]
            + cls.c * (wt_pc[["Na2O", "K2O"]].sum(axis=axis) ** 2.0)
        )

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:
        """
        Return one standard deviation errors on partition coefficients.

        Returns
        -------
        float, array-like
            partition coefficient errors
        """
        return cls.Kd_error

    @staticmethod
    def get_offset_parameters(n: int, *args, **kwargs) -> float | np.ndarray:
        """
        Randomly sample a standard normal distribution *n* times.

        n   : int
            sample amount.
        """
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate random samples of partition coefficient errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """
        return cls.get_error() * offset_parameters


class putirka2016_8d(Kd_model):
    """
    For liquid compositions with <45 wt.% SiO2 and > 8 wt.% Na2O + K2O
    """

    a = 0.6
    b = 1.3e-2
    c = 1.6e-2
    d = -1.73e-4
    e = 1.79e-2
    f = -2.6
    g = 2.11e-1
    h = 3.19e-5

    components = ["SiO2", "Al2O3", "Na2O", "K2O"]

    components = ["SiO2", "Na2O", "K2O"]

    Kd_error = 4.2e-2

    @classmethod
    def calculate_Kd(
        cls, melt_mol_fractions, T_K, P_bar, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate mineral-melt partition coefficients

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        P_bar :   float, array-like
            pressure in bar

        Returns
        -------
        float, array-like
            mineral-melt partition coefficients
        """
        wt_pc = melt_mol_fractions.wt_pc()
        wt_pc = check_components(wt_pc, components=cls.components)

        P_GPa = P_bar / 1e4

        axis = [0, 1][isinstance(wt_pc, pd.DataFrame)]

        Al_number = wt_pc["Al2O3"] / wt_pc[["Al2O3", "SiO2"]].sum(axis=axis)

        return (
            cls.a
            + cls.b * P_GPa
            + cls.c * wt_pc["SiO2"]
            + cls.d * (wt_pc["SiO2"] ** 2.0)
            + cls.e * wt_pc["Al2O3"]
            + cls.f * Al_number
            + cls.g * np.log(Al_number)
            + cls.h * (wt_pc[["Na2O", "K2O"]].sum(axis=axis) ** 3.0)
        )

    @classmethod
    def get_error(cls, *args, **kwargs) -> float:
        """
        Return one standard deviation errors on partition coefficients.

        Returns
        -------
        float, array-like
            partition coefficient errors
        """
        return cls.Kd_error

    @staticmethod
    def get_offset_parameters(n: int, *args, **kwargs) -> float | np.ndarray:
        """
        Randomly sample a standard normal distribution *n* times.

        n   : int
            sample amount.
        """
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate random samples of partition coefficient errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """
        return cls.get_error() * offset_parameters


class sun2020(Kd_model):
    """
    Calculate equilibrium Fe-Mg partition coefficients between olivine and melt according to equation 7 from Sun & Dasgupta (2020)\ [10]_.

    Sun, C., Dasgupta, R. (2020) Thermobarometry of CO2-rich, silica-undersaturated melts constrains cratonic lithosphere thinning through time in areas of kimberlitic magmatism. Earth and Planetary Sience Letters. 550
    """

    volatiles = ["H2O", "CO2", "F", "S", "Cl"]
    components = ["MgO", "Na2O", "H2O"]
    error = 0.03

    @classmethod
    @np.errstate(invalid="raise")
    def calculate_Kd(
        cls, melt_mol_fractions, Fe3Fe2, *args, **kwargs
    ) -> float | pd.Series:
        """
        Calculate Kds for given melt compositions and fixed forsterite content.

        Parameters
        ----------
        melt_mol_fractions : pandas Dataframe
            melt compositions in oxide mol fractions
        forsterite : float, array-like
            olivine forsterite contents.
        T_K : float, array-like
            temperatures in Kelvin
        P_bar: float, array-like
            pressures in bar

        Returns
        -------
        Kds : array-like
            Fe-Mg partition coefficients
        """
        composition = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        melt_volatile_free = _remove_elements(composition, drop=cls.volatiles)
        melt_per_oxygen = moles_per_oxygen(moles=melt_volatile_free)
        melt_wtpc = composition.wt_pc()

        Kd_Fe_total = np.exp(
            -1.65
            + 1.22 * np.sqrt(melt_per_oxygen["Mg1O"])
            + 2.45 * melt_per_oxygen["Na2O"]
            + 0.54 * melt_wtpc["H2O"] / 100
        )
        Fe3FeTotal = Fe3Fe2 / (1 + Fe3Fe2)

        Kd_Fe2 = Kd_Fe_total / (1 - Fe3FeTotal)

        return Kd_Fe2

    @classmethod
    def get_error(cls) -> float:
        """
        Calculate one standard deviation errors on Kds

        Returns
        -------
        error:  float
            one standard deviation error
        """

        return cls.error

    @classmethod
    def get_offset_parameters(cls, n=1):
        return np.random.normal(loc=0, scale=1, size=n)

    @classmethod
    def get_offset(cls, offset_parameters, *args, **kwargs):
        error = cls.get_error()
        return offset_parameters * error


_clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
# Collect all Kd_models in a dictionary.
Kd_models = {cls[0]: cls[1] for cls in _clsmembers if _is_Kd_model(cls[1])}
