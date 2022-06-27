import numpy as np
import pandas as pd
from scipy.constants import R  # J*K-1*mol-1
from .Fe_redox import FeRedox_QFM
from ..parse.validate import _check_argument





##### Toplis (2005) Fe-Mg olivine - melt exchange coefficient #####
###################################################################


def FeMg_toplis(T_K, P_bar, forsterite, SiO2_A):
    """
    Toplis (2005) Equation 10
    """
    return np.exp(
        (-6766 / (R * T_K) - 7.34 / R)
        + np.log(0.036 * SiO2_A - 0.22)
        + (3000 * (1 - 2 * forsterite) / (R * T_K))
        + (0.035 * (P_bar - 1) / (R * T_K))
    )


def toplis_Phi(molar_SiO2, molar_Na2O, molar_K2O):
    """

    Equation 12 from Toplis (2005) calculates a Phi parameter to correct SiO2 for alkali-bearing liquids.
    This expression is only valid for SiO2 <= 60 mol %


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
        if sum(np.array(molar_SiO2) > 60) > 1:
            raise RuntimeError("SiO2 >60 mol% present")
    except:
        if molar_SiO2 > 60:
            raise RuntimeError("SiO2 >60 mol%")

    return (0.46 * (100 / (100 - molar_SiO2)) - 0.93) * (molar_Na2O + molar_K2O) + (
        -5.33 * (100 / (100 - molar_SiO2)) + 9.69
    )


def toplis_SiO2_A(melt_mol_fractions):
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
    molar_concentrations = melt_mol_fractions * 100

    molar_SiO2 = molar_concentrations["SiO2"]
    molar_Na2O = molar_concentrations["Na2O"]
    molar_K2O = molar_concentrations["K2O"]

    Phi = toplis_Phi(molar_SiO2, molar_Na2O, molar_K2O)
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


##### Blundy (2020) Fe-Mg olivine - melt exchange coefficient #####
###################################################################


def FeMg_blundy(forsterite, Fe3Fe2_liquid, T_K):
    """
    Blundy et al., 2020, equation 8
    """

    Fe3FeTotal = Fe3Fe2_liquid / (1 + Fe3Fe2_liquid)

    return 0.3642 * (1 - Fe3FeTotal) * np.exp(312.7 * (1 - 2 * forsterite) / T_K)


def equilibrium_forsterite(Kd, Fe2Mg):
    """
    Parameters
    ----------
    Kd :
        (Fe_olivine / Fe_liquid) * (Mg_liquid / Mg_olivine) partitioning coefficient
    Fe2Mg
        melt Fe2+/Mg ratio

    Returns
    -------
    Equilibrium forsterite fraction as Mg/(Mg + Fe)
    """

    return 1 / (1 + Kd * Fe2Mg)


class Kd_FeMg_vectorised:
    def blundy(
        melt_mol_fractions: pd.DataFrame,
        forsterite_initial,
        T_K,
        Fe3Fe2,
        **kwargs,
    ):

        """
        Equation 8 by Blundy (2020) iteratively solved for forsterite content

        Parameters
        ----------
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        olivine_forsterite
            forsterite fraction in olivine as Mg / (Mg + Fe)
        T_K :
            Temperature in Kelvin
        Fe3Fe2 :
            melt Fe2+/Fe3+ ratio
        """

        for name in ["T_K", "forsterite_initial"]:
            param = locals()[name]
            if isinstance(param, pd.Series):
                if not melt_mol_fractions.index.equals(param.index):
                    raise RuntimeError(f"Melt and {name} indices don't match")

        # Convert everything to Series for easier looping
        if isinstance(forsterite_initial, (int, float)):
            forsterite = pd.Series(forsterite_initial, index=melt_mol_fractions.index)
        else:
            forsterite = forsterite_initial.copy()
        if isinstance(T_K, (int, float)):
            T_K = pd.Series(T_K, index=melt_mol_fractions.index)
        if isinstance(Fe3Fe2, (int, float)):
            Fe3Fe2 = pd.Series(Fe3Fe2, index=melt_mol_fractions.index)

        fo_converge_default = 0.001
        fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

        # initialise Kds
        Kd = FeMg_blundy(forsterite, Fe3Fe2, T_K)
        # Liquid Fe2+/Fe(total)
        Fe2Fe_total = 1 / (1 + Fe3Fe2)
        Fe2Mg = (melt_mol_fractions["FeO"] / melt_mol_fractions["MgO"]) * Fe2Fe_total
        # Equilibrium forsterite content according to Kd
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
        # Difference between observed Fo and equilibrium Fo
        forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

        iterate = forsterite_delta > fo_converge
        # iterate until equilibrium forsterite content doesn't change any more
        while sum(iterate) > 1:
            Kd.loc[iterate] = FeMg_blundy(
                forsterite_EQ[iterate], Fe3Fe2[iterate], T_K[iterate]
            )
            forsterite[iterate] = forsterite_EQ[iterate].copy()
            forsterite_EQ.loc[iterate] = 1 / (1 + Kd[iterate] * Fe2Mg[iterate])
            forsterite_delta.loc[iterate] = (
                abs(forsterite[iterate] - forsterite_EQ[iterate]) / forsterite[iterate]
            )
            iterate = forsterite_delta > fo_converge

        return Kd

    def toplis(
        melt_mol_fractions: pd.DataFrame,
        forsterite_initial,
        T_K,
        P_bar,
        Fe3Fe2,
        **kwargs,
    ):
        """
        Equation 10 of Toplis (2005) iteratively solved for forsterite content

        Parameters
        ----------
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        olivine_forsterite
            forsterite fraction in olivine as Mg * 100 / (Mg + Fe)
        T_K :
            Temperature in Kelvin
        P_bar :
            Pressure in bar
        Fe3Fe2 :
            melt Fe3+/Fe2+ ratio
        """

        for name in ["T_K", "P_bar", "forsterite_initial"]:
            param = locals()[name]
            if isinstance(param, pd.Series):
                if not melt_mol_fractions.index.equals(param.index):
                    raise RuntimeError(f"Melt and {name} indices don't match")

        # Convert everything to Series for easier looping
        if isinstance(forsterite_initial, (int, float)):
            forsterite = pd.Series(forsterite_initial, index=melt_mol_fractions.index)
        else:
            forsterite = forsterite_initial.copy()
        if isinstance(T_K, (int, float)):
            T_K = pd.Series(T_K, index=melt_mol_fractions.index)
        if isinstance(Fe3Fe2, (int, float)):
            Fe3Fe2 = pd.Series(Fe3Fe2, index=melt_mol_fractions.index)
        if isinstance(P_bar, (int, float)):
            P_bar = pd.Series(P_bar, index=melt_mol_fractions.index)

        fo_converge_default = 0.001
        fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

        SiO2mol_A = toplis_SiO2_A(melt_mol_fractions)
        # initialise Kds
        Kd = FeMg_toplis(T_K, P_bar, forsterite, SiO2mol_A)
        # Liquid Fe2+/Fe(total)
        Fe2Fe_total = 1 / (1 + Fe3Fe2)
        # liquid Fe2+/Mg
        Fe2Mg = (melt_mol_fractions["FeO"] * Fe2Fe_total)/ melt_mol_fractions["MgO"]
        # Equilibrium forsterite content according to Kd
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
        # Difference between observed Fo and equilibrium Fo
        forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

        iterate = forsterite_delta > fo_converge
        # iterate until equilibrium forsterite content doesn't change any more
        while sum(iterate) > 1:
            Kd[iterate] = FeMg_toplis(
                T_K[iterate],
                P_bar[iterate],
                forsterite_EQ.loc[iterate],
                SiO2mol_A.loc[iterate],
            )
            forsterite[iterate] = forsterite_EQ[iterate].copy()
            forsterite_EQ[iterate] = 1 / (1 + Kd[iterate] * Fe2Mg[iterate])
            forsterite_delta[iterate] = (
                abs(forsterite[iterate] - forsterite_EQ[iterate]) / forsterite[iterate]
            )
            iterate = forsterite_delta > fo_converge

        return Kd


class Kd_FeMg:
    def blundy(
        melt_mol_fractions: pd.Series,
        forsterite: float,
        T_K: float,
        Fe3Fe2: float,
        **kwargs,
    ):

        """
        Equation 8 by Blundy (2020) iteratively solved for forsterite content

        Parameters
        ----------
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        olivine_forsterite
            forsterite fraction in olivine as Mg / (Mg + Fe)
        T_K :
            Temperature in Kelvin
        Fe3Fe2 :
            melt Fe2+/Fe3+ ratio
        """

        fo_converge_default = 0.001
        fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

        # initialise Kds
        Kd = FeMg_blundy(forsterite, Fe3Fe2, T_K)
        # Liquid Fe2+/Fe(total)
        Fe2Fe_total = 1 / (1 + Fe3Fe2)
        Fe2Mg = (melt_mol_fractions["FeO"] * Fe2Fe_total)/ melt_mol_fractions["MgO"]
        # Equilibrium forsterite content according to Kd
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
        # Difference between observed Fo and equilibrium Fo
        forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

        iterate = forsterite_delta > fo_converge
        # iterate until equilibrium forsterite content doesn't change any more
        while iterate:
            Kd = FeMg_blundy(forsterite_EQ, Fe3Fe2, T_K)
            forsterite = forsterite_EQ.copy()
            forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
            forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite
            iterate = forsterite_delta > fo_converge

        return Kd

    def toplis(
        melt_mol_fractions: pd.Series,
        forsterite: float,
        T_K: float,
        P_bar: float,
        Fe3Fe2: float,
        **kwargs,
    ):
        """
        Equation 10 of Toplis (2005) iteratively solved for forsterite content

        Parameters
        ----------
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        olivine_forsterite
            forsterite fraction in olivine as Mg * 100 / (Mg + Fe)
        T_K :
            Temperature in Kelvin
        P_bar :
            Pressure in bar
        Fe3Fe2 :
            melt Fe3+/Fe2+ ratio
        """

        fo_converge_default = 0.001
        fo_converge = kwargs.setdefault("fo_converge", fo_converge_default)

        SiO2mol_A = toplis_SiO2_A(melt_mol_fractions)

        # initialise Kd
        Kd = FeMg_toplis(T_K, P_bar, forsterite, SiO2mol_A)

        # Liquid Fe2+/Fe(total)
        Fe2Fe_total = 1 / (1 + Fe3Fe2)
        # liquid Fe2+/Mg
        Fe2Mg = (melt_mol_fractions["FeO"] * Fe2Fe_total)/ melt_mol_fractions["MgO"]
        # Equilibrium forsterite content according to Kd
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)

        # Difference between observed Fo and equilibrium Fo
        forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

        iterate = forsterite_delta > fo_converge
        # iterate until equilibrium forsterite content doesn't change any more
        while iterate:

            Kd = FeMg_toplis(
                T_K,
                P_bar,
                forsterite_EQ,
                SiO2mol_A,
            )

            forsterite = forsterite_EQ.copy()
            forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
            forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite
            iterate = forsterite_delta > fo_converge

        return Kd
