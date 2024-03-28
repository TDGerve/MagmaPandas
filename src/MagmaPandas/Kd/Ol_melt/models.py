import inspect
import sys

import numpy as np
import pandas as pd
from scipy.constants import R  # J*K-1*mol-1

from MagmaPandas import Fe_redox, configuration
from MagmaPandas.fO2 import fO2_QFM
from MagmaPandas.Kd.Kd_baseclass import Kd_model


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
    def calculate_Kd(
        cls,
        melt_mol_fractions: pd.DataFrame,
        forsterite: float | pd.Series,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        *args,
        **kwargs,
    ) -> float | pd.Series:
        """
        Calculate equilibrium Kds for given melt compositions

        Parameters
        ----------
        melt_mol_fractions : pandas Dataframe
            melt compositions in oxide mol fractions
        forsterite : float, array-like
            initial olivine forsterite contents. Forsterite values are iteratively adjusted and initial values are not necessarily in Fe-Mg equilibrium with melts.
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
    def _get_Fe3FeTotal(Melt_mol_fractions, T_K, P_bar=1, **kwargs):

        Fe3Fe2_model = Fe_redox.borisov
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


_clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

# Collect all Kd_models in a dictionary.
Kd_models = {cls[0]: cls[1] for cls in _clsmembers if _is_Kd_model(cls[1])}
