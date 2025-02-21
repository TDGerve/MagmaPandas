import inspect
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy.constants import Avogadro, R
from scipy.optimize import fsolve

from MagmaPandas.EOSs.birch_murnaghan import birch_murnaghan_4th_order
from MagmaPandas.Fe_redox.Fe3Fe2_baseclass import Fe3Fe2_model
from MagmaPandas.Fe_redox.Fe3Fe2_errors import (
    error_params_1bar,
    error_params_high_pressure,
)
from MagmaPandas.parse_io import check_components, make_equal_length, make_iterable


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

    components = ["SiO2", "TiO2", "MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5"]

    error_params = [0.03152663, 0.02589086, 0.6874305, 7.73489745]

    @classmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions: pd.DataFrame, T_K, fO2, *args, **kwargs
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

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
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
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


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

    error_params = [7.34976440e-02, 1.81801171e-02, 9.87644316e-01, 2.27525692e02]

    @classmethod
    def calculate_Fe3Fe2(cls, melt_mol_fractions, T_K, fO2, P_bar):
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
        P_bar      :      float, array-like
            Pressure in bar


        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )

        P_Pa = P_bar / 1e5

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
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class jayasuriya(Fe3Fe2_model):
    """
    Jayasuriya et al. (2004)\ [3]_ equation 12.
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

    error_params = [2.10272743e-01, -1.27359105e-02, 9.82853119e-01, 1.75114023e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )

        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        sumComponents = moles[cls.components].mul(cls.dCoefficients).sum(axis=axis)

        return 2 * np.exp(cls.a * np.log(fO2) + 12420 / T_K - cls.c + sumComponents)

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class putirka2016_6b(Fe3Fe2_model):
    """
    Putirka (2016)\ [4]_ equation 6b.
    """

    components = ["Na2O", "K2O", "Al2O3", "SiO2", "CaO"]

    a = 6.53
    b = 10813.8
    c = 0.19
    d = 12.4
    e = 3.44
    f = 4.15

    error_params = [1.87571331e-01, 2.58857703e-03, 8.65663335e-01, 2.43304009e01]

    @classmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )

        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        part1 = -cls.a + cls.b / T_K
        part2 = cls.c * np.log(fO2) + 12.4 * moles[["Na2O", "K2O"]].sum(axis=axis)
        part3 = (
            -cls.e * (moles["Al2O3"] / moles[["Al2O3", "SiO2"]].sum(axis=axis))
            + cls.f * moles["CaO"]
        )

        return 2 * np.exp(part1 + part2 + part3)

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class putirka2016_6c(Fe3Fe2_model):
    """
    Putirka (2016)\ [4]_ equation 6c.
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

    a = -6.75
    b = 10634.9
    c = 0.195
    d = 7.9
    e = -4.6
    f = 0.54
    g = -53.4
    h = 1.07

    error_params = [-3.12797578e-02, 6.71834626e-02, 9.81552910e-01, 1.30748389e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        axis = [0, 1][isinstance(moles, pd.DataFrame)]

        NBO_T = cls._NBO_T(cations=moles.cations())

        part_1 = cls.a + cls.b / T_K + cls.c * np.log(fO2)
        part_2 = cls.d * moles[["Na2O", "K2O"]].sum(axis=axis) + cls.e * moles["MgO"]
        part_3 = (
            cls.f * (moles["MgO"] / moles[["MgO", "FeO"]].sum(axis=axis))
            + cls.g * moles["P2O5"]
            + cls.h * NBO_T
        )

        return 2 * np.exp(part_1 + part_2 + part_3)

    @staticmethod
    def _NBO_T(cations):

        axis = [0, 1][isinstance(cations, pd.DataFrame)]

        if axis == 1:
            Al_IV = cations.apply(
                lambda row: min(
                    row["Al"],
                    row[["Na", "K"]].sum() + 2 * row[["Ca", "Mg"]].sum(),
                ),
                axis=1,
            )
        else:
            Al_IV = min(
                cations["Al"],
                cations[["Na", "K"]].sum() + 2 * cations[["Ca", "Mg"]].sum(),
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
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Deng2020(Fe3Fe2_model):
    """
    Deng et al. (2020)\ [5]_
    """

    components = [
        "SiO2",
        "Al2O3",
        "FeO",
        "MgO",
        "CaO",
        "K2O",
        "Na2O",
        "TiO2",
        "P2O5",
    ]
    gibbs_parameters = pd.Series(
        {
            "a": -3.310e5,
            "b": -190.379,
            "c": 14.785,
            "d": -1.649e-3,
            "e": 9.348e6,
            "f": 1.077e4,
        }
    )

    margules = pd.Series(  # fit 3
        {
            "Mg": 68629,
            "Si": 4601,
            "Al": 40923,
            "Ca": -58109,
            "Na": 0,
            "K": -59584,
            "P": 0,
            "Ti": 0,
        }
    )
    Fe_margules = -14210  # fit 3

    eos_params = pd.DataFrame(
        [
            [1180.114014, 1204.763652, 1192.011066, 1256.727179],
            [26.94713861, 23.19530062, 23.95435759, 16.12613905],
            [2.802531871, 3.216089358, 3.32104996, 4.584011905],
            [0.012313472, 0.009340183, -0.008912497, -0.177152954],
        ],
        columns=list(zip(*[["12.5molpc"] * 2 + ["25molpc"] * 2, ["Fe2", "Fe3"] * 2])),
        index=["V_0", "K_0", "Kprime_0", "Kprime_prime_0"],
    )

    thermal_pressure_params = pd.DataFrame(
        [
            [35.79397483, 34.52616394, 31.34712676, 30.38414264],
            [71.10313668, 68.64429623, 62.48520005, 59.10950152],
            [36.59545225, 35.27069116, 32.4675829, 29.64971394],
        ],
        columns=list(zip(*[["12.5molpc"] * 2 + ["25molpc"] * 2, ["Fe2", "Fe3"] * 2])),
        index=["a", "b", "c"],
    )
    formula_units = {"12.5molpc": 2, "25molpc": 4}
    T_K_ref = 3000
    A3_to_cm3 = 1e-24

    error_params = [2.17554822e-01, -4.67257455e-03, 9.83474292e-01, 2.27563235e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls,
        melt_mol_fractions,
        T_K: float | np.ndarray,
        P_bar: float | np.ndarray,
        fO2: float | np.ndarray,
        melt_Fe: str = "12.5molpc",
        Fe3Fe2_init=0.3,
        total_Fe="FeO",
        **kwargs,
    ):

        moles = (
            melt_mol_fractions.to_frame().T
            if isinstance(melt_mol_fractions, pd.Series)
            else melt_mol_fractions.copy()
        )  # force to DataFrame

        # TODO  CHECK UNITS: BARS, PASCAL OR GPA??

        dVdP = cls._dVdP(T_K=T_K, P_bar=P_bar, melt_Fe=melt_Fe, **kwargs)
        gibbs0 = cls._gibbs0(T_K)

        try:
            if len(moles) != (l := len(T_K)):
                moles = moles.loc[moles.index.repeat(l)].reset_index(drop=True)
        except TypeError:
            pass
        # force everything to an iteratble
        T_K, fO2, gibbs0, dVdP = make_equal_length(T_K, fO2, gibbs0, dVdP)

        Fe3Fe2_func = (
            lambda Fe3Fe2_guess, moles, T, fO2, G, dVdP: cls._Fe3Fe2_solver(
                melt_mol_fractions=moles,
                T_K=T,
                fO2=fO2,
                gibbs0=G,
                dVdP=dVdP,
                Fe3Fe2=Fe3Fe2_guess[0],
                total_Fe=total_Fe,
            )
            - Fe3Fe2_guess
        )

        Fe3Fe2 = pd.Series(dtype=float)

        # for (n, X), T, fO2, G, VP in zip(moles.iterrows(), T_K, fO2_GPa, gibbs0, dVdP):
        #     Fe3Fe2.loc[n] = cls._Fe3Fe2_solver(
        #         melt_mol_fractions=X, T_K=T, fO2_GPa=fO2, gibbs0=G, dVdP=VP
        #     )

        for (n, X), T, f, G, VP in zip(moles.iterrows(), T_K, fO2, gibbs0, dVdP):
            Fe3Fe2.loc[n] = fsolve(Fe3Fe2_func, x0=Fe3Fe2_init, args=(X, T, f, G, VP))[
                0
            ]

        return Fe3Fe2.squeeze()

    @classmethod
    def _Fe3Fe2_solver(
        cls,
        melt_mol_fractions: pd.Series,
        T_K: float,
        fO2: float,
        gibbs0: float,
        dVdP: float,
        Fe3Fe2=0.3,
        total_Fe="FeO",
    ) -> float | np.ndarray:
        """
        Equation 3 solved for Fe3Fe2
        """

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        cations = moles.FeO_Fe2O3_calc(Fe3Fe2, wtpc=False, total_Fe=total_Fe).cations()
        Fe_activities = cls._Fe_activities(cations, T_K)

        Fe3Fe2 = np.exp(-(gibbs0 + dVdP) / (R * T_K) - Fe_activities + np.log(fO2) / 4)

        return Fe3Fe2

    @classmethod
    def _gibbs0(cls, T_K):
        """
        Calculates the Gibbs free energy of equation 1 at reference P and T. Supplementary Note 1
        """

        return (
            cls.gibbs_parameters["a"]
            + cls.gibbs_parameters["b"] * T_K
            + cls.gibbs_parameters["c"] * T_K * np.log(T_K)
            + cls.gibbs_parameters["d"] * (T_K**2)
            + cls.gibbs_parameters["e"] / T_K
            + cls.gibbs_parameters["f"] * np.sqrt(T_K)
        )

    @classmethod
    def _thermal_pressure_coeff(cls, V, V0, params):
        """
        Calculates the thermal pressure coefficient in equation 2.
        """

        return (
            params["a"] - params["b"] * (V / V0) + params["c"] * np.power(V / V0, 2)
        ) / 1000

    @classmethod
    @np.vectorize(excluded=["cls", "phase", "melt_Fe"])
    def _calculate_volume(
        cls, T_K: float, P_bar: float, phase: str, melt_Fe: str
    ) -> float:
        """
        Equation 2 solved for volume.

        Returns
        -------
        volume  :   float
            volume in Angstrom^3
        """
        P_GPa = P_bar / 1e4
        eos_params = cls.eos_params[(f"{melt_Fe}", f"{phase}")]
        thermal_pressure_params = cls.thermal_pressure_params[
            (f"{melt_Fe}", f"{phase}")
        ]

        v_init = eos_params["V_0"] - 6 * P_GPa

        func = (
            lambda v: (
                birch_murnaghan_4th_order(V=v, params=eos_params)
                + cls._thermal_pressure_coeff(
                    V=v, V0=eos_params["V_0"], params=thermal_pressure_params
                )
                * (T_K - cls.T_K_ref)
            )
            - P_GPa
        )

        V_sol = fsolve(func, x0=v_init)

        return V_sol

    @classmethod
    def _deltaV(cls, T_K, P_bar, melt_Fe):
        """
        Volume change across equation 1.

        Returns
        -------
        dV : float
            delta Volume in cm^3/mol (Fe3 - Fe2).
        """
        if isinstance(T_K, (float, int)):
            T_K = np.ones_like(P_bar) * T_K

        Fe2_volume = (
            cls._calculate_volume(T_K, P_bar, "Fe2", melt_Fe)
            / cls.formula_units[melt_Fe]
            * Avogadro
            * cls.A3_to_cm3
        )  # cm3/mol
        Fe3_volume = (
            cls._calculate_volume(T_K, P_bar, "Fe3", melt_Fe)
            / cls.formula_units[melt_Fe]
            * Avogadro
            * cls.A3_to_cm3
        )  # cm3/mol

        return Fe3_volume - Fe2_volume

    @classmethod
    @np.vectorize(excluded=["cls", "melt_Fe", "P_bar_min", "P_bar_step"])
    def _dVdP(
        cls,
        T_K: float | np.ndarray,
        P_bar: float | np.ndarray,
        melt_Fe: str,
        Pbar_min=1.0,
        Pbar_step=5e2,
        **kwargs,
    ):
        """
        delta V integrated for pressure

        Returns
        -------
        dVdP : float, array-like
            dVdP in m3.Pascal
        """

        P_array = np.arange(Pbar_min, P_bar + Pbar_step, Pbar_step)

        dV = cls._deltaV(T_K=T_K, P_bar=P_array, melt_Fe=melt_Fe) * 1e-6  # cm3 to m3

        # TODO  CHECK UNITS: CM3, M3, A3, BAR, PASCAL OR GPA??
        return np.trapz(dV, P_array * 1e5)  # bar to Pascal

    @classmethod
    def _Fe_activities(cls, melt_cation_fractions, T_K):
        """
        Calculate FeO1.5 and FeO activities in silicate melt. Supplementary Note 2
        """

        axis = [0, 1][isinstance(melt_cation_fractions, pd.DataFrame)]

        sum_margules = melt_cation_fractions.mul(cls.margules).sum(axis=axis)

        LN_aFe3_aFe2 = sum_margules / (R * T_K) + (
            melt_cation_fractions["Fe"] - melt_cation_fractions["Fe3"]
        ) * cls.Fe_margules / (R * T_K)

        return LN_aFe3_aFe2

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Oneill2006(Fe3Fe2_model):
    """
    O'neill et al. (2006)\ [6]_
    """

    components = ["MgO", "CaO", "Na2O", "K2O", "Al2O3", "P2O5", "FeO"]

    error_params = [2.91079515e-01, -1.37913265e-02, 9.82946462e-01, 1.90242796e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls,
        melt_mol_fractions,
        P_bar,
        T_K,
        fO2,
        Fe3Fe2_init=0.3,
        total_Fe="FeO",
        *args,
        **kwargs,
    ):
        """
        Calculate melt |Fe3Fe2| ratios with equation 10.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        P_bar : float, array-like
            pressure in bar
        T_K :   float, array-like
            temperature in Kelvin
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """

        moles = (
            melt_mol_fractions.to_frame().T
            if isinstance(melt_mol_fractions, pd.Series)
            else melt_mol_fractions.copy()
        )  # force to DataFrame

        try:
            # extend the moles dataframe is multiple P,T, fO2 conditions are given for a single composition.
            if len(moles) != (l := len(T_K)):
                moles = moles.loc[moles.index.repeat(l)].reset_index(drop=True)
        except TypeError:
            pass
        # force everything into an iterable
        P_bar, T_K, fO2 = make_equal_length(P_bar, T_K, fO2)

        Fe3Fe2_func = (
            lambda Fe3Fe2_guess, moles, P, T, fO2: cls._calculate_Fe3Fe2(
                melt_mol_fractions=moles,
                P_bar=P,
                T_K=T,
                fO2=fO2,
                Fe3Fe2=Fe3Fe2_guess[0],
                total_Fe=total_Fe,
            )
            - Fe3Fe2_guess[0]
        )

        Fe3Fe2 = pd.Series(dtype=float)

        for (n, X), P, T, f in zip(moles.iterrows(), P_bar, T_K, fO2):
            Fe3Fe2.loc[n] = fsolve(Fe3Fe2_func, x0=Fe3Fe2_init, args=(X, P, T, f))[0]

        return Fe3Fe2.squeeze()

    @classmethod
    def _calculate_Fe3Fe2(
        cls, melt_mol_fractions, P_bar, T_K, fO2, Fe3Fe2=0.3, total_Fe: str = "FeO"
    ):

        P_GPa = P_bar / 1e4
        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        cations = moles.FeO_Fe2O3_calc(Fe3Fe2, wtpc=False, total_Fe=total_Fe).cations()

        part_1 = (
            -28144
            + 3905 * cations["Mg"]
            - 13359 * cations["Ca"]
            - 14858 * cations["Na"]
            - 9805 * cations["K"]
            + 10906 * cations["Al"]
            + 110971 * cations["P"]
            - 11952 * (cations["Fe"] - cations["Fe3"])
        ) / T_K

        part_2 = (
            13.95
            + (33122 / T_K - 5.24) * ((1 + 0.241 * P_GPa) ** (3 / 4) - 1)
            - (39156 / T_K - 6.17) * ((1 + 0.132 * P_GPa) ** (3 / 4) - 1)
        )

        return 10 ** ((np.log10(fO2) - part_1 - part_2) / 4)

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Oneill2018(Fe3Fe2_model):
    """
    O'Neill et al. (2018)\ [7]_
    """

    components = ["CaO", "Na2O", "K2O", "P2O5"]

    error_params = [-1.99768682e-02, 7.41622563e-02, 9.80282471e-01, 1.28905955e02]

    @classmethod
    def calculate_Fe3Fe2(cls, melt_mol_fractions, T_K, fO2, *args, **kwargs):
        """
        Calculate melt |Fe3Fe2| ratios with equation 9a.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        fO2 :   float, array-like
            Oxygen fugacity in bar

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        cations = check_components(
            composition=melt_mol_fractions, components=cls.components
        ).cations()
        deltaQFM = cls._deltaQFM(fO2, T_K)

        return 10 ** (
            0.25 * deltaQFM
            - 1.36
            + 2.4 * cations["Ca"]
            + 2.0 * cations["Na"]
            + 3.7 * cations["K"]
            - 2.4 * cations["P"]
        )

    @staticmethod
    def _deltaQFM(fO2, T_K):
        return np.log10(fO2) - (8.58 - 25050 / T_K)

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Armstrong2019(Fe3Fe2_model):
    """
    Armstrong et al. (2019)\ [8]_

    Calibrated with one andesitic and one MORB compositions + data from O'neill et al. (2006)\ [6]_ and Zhang et al. (2017)\ [9]_
    """

    components = ["MgO", "CaO", "Na2O", "K2O", "Al2O3"]

    eos_parameters = pd.DataFrame(
        {"K_0": [37, 12.6], "Kprime_0": [8, 1.3]}, index=["Fe2", "Fe3"]
    )
    # K_0: GPa, Kprime_0: GPa-1

    margules = pd.Series(
        {
            "Mg": -2248,
            "Ca": 7690,
            "Na": 8553,
            "K": 5644,
            "Al": -6278,
        }
    )  # Jayasuriya et al., (2004)
    Fe_margules = 6880

    error_params = [2.20466771e-01, 9.51885060e-03, 9.85463396e-01, 2.05555278e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls,
        melt_mol_fractions,
        T_K,
        P_bar,
        fO2,
        Fe3Fe2_init=0.3,
        total_Fe="FeO",
        *args,
        **kwargs,
    ):
        """
        Calculate melt |Fe3Fe2| ratios.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar : float, array-like
            pressure in bar
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """

        moles = (
            melt_mol_fractions.to_frame().T
            if isinstance(melt_mol_fractions, pd.Series)
            else melt_mol_fractions.copy()
        )  # force to DataFrame

        try:
            # extend the moles dataframe is multiple P,T, fO2 conditions are given for a single composition.
            if len(moles) != (l := len(T_K)):
                moles = moles.loc[moles.index.repeat(l)].reset_index(drop=True)
        except TypeError:
            pass
        # force everything into iterables
        P_bar, T_K, fO2 = make_equal_length(P_bar, T_K, fO2)

        Fe3Fe2_func = (
            lambda Fe3Fe2_guess, moles, T, P, fO2: cls._calculate_Fe3Fe2(
                melt_mol_fractions=moles,
                T_K=T,
                P_bar=P,
                fO2=fO2,
                Fe3Fe2=Fe3Fe2_guess[0],
                total_Fe=total_Fe,
            )
            - Fe3Fe2_guess[0]
        )

        Fe3Fe2 = pd.Series(dtype=float)

        for (n, X), P, T, f in zip(moles.iterrows(), P_bar, T_K, fO2):
            Fe3Fe2.loc[n] = fsolve(Fe3Fe2_func, x0=Fe3Fe2_init, args=(X, T, P, f))[0]

        return Fe3Fe2.squeeze()

    @classmethod
    def _calculate_Fe3Fe2(
        cls,
        melt_mol_fractions,
        T_K,
        P_bar,
        fO2,
        Fe3Fe2=0.3,
        total_Fe: str = "FeO",
        *args,
        **kwargs,
    ):
        """
        Calculate melt |Fe3Fe2| ratios with Supplementary Materials eq. S12

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar : float, array-like
            pressure in bar
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """

        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        cations = moles.FeO_Fe2O3_calc(Fe3Fe2, wtpc=False, total_Fe=total_Fe).cations()

        gibbs0 = cls._Gibbs0(T_K=T_K)
        dVdP = cls._dVdP(P_bar=P_bar, T_K=T_K)
        Fe_activities = cls._Fe_activities(melt_cation_fractions=cations, T_K=T_K)

        return np.exp(np.log(fO2) / 4 - (gibbs0 + dVdP) / (R * T_K) + Fe_activities)

    @staticmethod
    def _Gibbs0(T_K):
        """
        Supplementary Materials equation S4
        """
        return -(16201 / T_K - 8.031) * (R * T_K)

    @classmethod
    def _dVdP(cls, P_bar, T_K, **kwargs):

        VdP_Fe3 = cls._VdP(P_bar=P_bar, T_K=T_K, phase="Fe3")
        VdP_Fe2 = cls._VdP(P_bar=P_bar, T_K=T_K, phase="Fe2")

        return VdP_Fe3 - VdP_Fe2

    @classmethod
    def _VdP(cls, P_bar, T_K, phase: str):
        """
        Supplementary Materials equation S10
        """
        params = cls.eos_parameters.loc[phase]
        P_GPa = P_bar * 1e5 / 1e9
        Kprime_prime_0 = (
            -params["Kprime_0"] / params["K_0"]
        )  # Supplementary Materials ref. 32

        V_0 = cls._V_0(T_K, phase=phase)  # mm3/mol

        # Supplementary Materials eqs. S7-S9
        a = (1 + params["Kprime_0"]) / (
            1 + params["Kprime_0"] + params["K_0"] * Kprime_prime_0
        )
        b = params["Kprime_0"] / params["K_0"] - Kprime_prime_0 / (
            1 + params["Kprime_0"]
        )
        c = (1 + params["Kprime_0"] + params["K_0"] * Kprime_prime_0) / (
            params["Kprime_0"] ** 2
            + params["Kprime_0"]
            - params["K_0"] * Kprime_prime_0
        )

        part_1 = a * (1 - (1 + b * P_GPa) ** (1 - c))
        part_2 = b * (c - 1) * P_GPa

        # TODO CHECK UNITS M3, CM3, BAR, PASCAL? -> M3 * (pressure, GPa?)
        return P_GPa * V_0 * (1 - a + part_1 / part_2)  # mm3.GPa = m3.Pa

    @staticmethod
    def _V_0(T_K, phase: str):
        """
        Partial molar volume of FeO or FeO1.5 in mm3/mol.

        Supplementary Materials ref 60, 61
        Lange & Carmichael (1987) values from table 8
        Lange & Carmichael (1990)
        """
        # TODO CHECK UNITS
        return {
            "Fe2": (13650 + 2.92 * (T_K - 1673)),
            "Fe3": (21070 + 4.54 * (T_K - 1673)),
        }[phase]

    @classmethod
    def _Fe_activities(cls, melt_cation_fractions, T_K):
        """
        Calculate FeO1.5 and FeO activities in silicate melt. Supplementary Materials eq. S11
        """

        axis = [0, 1][isinstance(melt_cation_fractions, pd.DataFrame)]

        sum_margules = melt_cation_fractions.mul(cls.margules).sum(axis=axis)

        LN_aFe3_aFe2 = (
            sum_margules / T_K
            + (melt_cation_fractions["Fe"] - melt_cation_fractions["Fe3"])
            * cls.Fe_margules
            / T_K
        )

        return LN_aFe3_aFe2

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Zhang2017(Fe3Fe2_model):
    """
    Zhang et al. (2017)\ [9]_

    Only calibrated with an andesitic melt composition
    """

    parameters = pd.DataFrame(
        {
            "dVdT": [2.92, 3.69],
            "a": [-6.376, -6.627],
            "b": [107257, 110729],
            "c": [15095, 15243],
            "d": [8.27e-2, 1.137e-1],
            "K_0": [36.61, 27.11],
        },
        index=["LC", "Guo"],
    )  # dVdT: J/GPa/K

    error_params = [1.55144202e-01, 3.39642134e-03, 9.87246138e-01, 2.68123961e02]

    @classmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions, T_K, P_bar, fO2, parameters="LC", *args, **kwargs
    ):
        """
        Calculate melt |Fe3Fe2| ratios with equation 11.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar : float, array-like
            pressure in bar
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        params = cls.parameters.loc[parameters]
        d = params["d"]  # 4 / params["K_0"]

        # P_bar, T_K = make_iterable(P_bar, T_K)
        # array_length_one = np.array([len(T_K), len(P_bar)]) == 1

        # if array_length_one.any():
        #     array_length = np.array([len(T_K), len(P_bar)]).max()
        #     P_bar, T_K = [
        #         np.ones(array_length) * i[0] if len(i) == 1 else i for i in (P_bar, T_K)
        #     ]

        P_GPa = P_bar * 1e5 / 1e9

        part_1 = np.log(fO2) / 4 + params["a"] + params["b"] / (R * T_K)
        part_2 = (
            -(20170 + 4.54 * (T_K - 1673))
            * (16.6 / 3)
            * ((1 + 0.241 * P_GPa) ** (0.75) - 1)
            / (R * T_K)
        )
        part_3 = (params["c"] + params["dVdT"] * (T_K - 1673)) * (4 / (3 * d))
        part_4 = ((1 + d * P_GPa) ** (0.75) - 1) / (R * T_K)

        LN_Fe3Fe2 = part_1 + part_2 + part_3 * part_4

        return np.exp(LN_Fe3Fe2)

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Hirschmann2022(Fe3Fe2_model):
    """
    Hirschmann (2022)\ [10]_
    """

    components = ["SiO2", "TiO2", "MgO", "CaO", "Na2O", "K2O", "P2O5", "Al2O3"]

    params = pd.Series(
        {
            "a": 0.1917,
            "b": -1.961,
            "c": 4158.1,
            "dCp": 33.25,
            "T0": 1673.15,
            "y1": -520.46,
            "y2": -185.37,
            "y3": 494.39,
            "y4": 1838.34,
            "y5": 2888.48,
            "y6": 3473.68,
            "y7": -4473.6,
            "y8": -1245.09,
            "y9": -1156.86,
        }
    )
    # TODO calculate parameters
    error_params = [1, 1, 1, 1]

    @classmethod
    def calculate_Fe3Fe2(
        cls,
        melt_mol_fractions,
        T_K,
        P_bar,
        fO2,
        dVdP_method="Armstrong2019",
        *args,
        **kwargs,
    ):
        """
        Calculate melt |Fe3Fe2| ratios with equation 21.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar : float, array-like
            pressure in bar
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        cations = moles.cations()

        compositional_term = cls._compositional_term(cations=cations)
        dVdP = cls._dVdP(T_K=T_K, P_bar=P_bar, method=dVdP_method)

        part_1 = (
            cls.params["a"] * np.log10(fO2) + cls.params["b"] + cls.params["c"] / T_K
        )
        part_2 = (
            -cls.params["dCp"]
            / (R * np.log(10))
            * (1 - cls.params["T0"] / T_K - np.log(T_K / cls.params["T0"]))
        )
        part_3 = -dVdP / (R * T_K * np.log(10))
        part_4 = 1 / T_K * compositional_term

        return 10 ** (part_1 + part_2 + part_3 + part_4)

    @staticmethod
    def _dVdP(T_K, P_bar, method="Armstrong2019"):

        model = getattr(sys.modules[__name__], method)

        return model._dVdP(T_K=T_K, P_bar=P_bar, melt_Fe="12.5molpc")

    @classmethod
    def _compositional_term(cls, cations):

        axis = [0, 1][isinstance(cations, pd.DataFrame)]

        part_1 = (
            cations[["Si", "Ti", "Mg", "Ca", "Na", "K", "P"]]
            .mul(cls.params[["y1", "y2", "y3", "y4", "y5", "y6", "y7"]].values)
            .sum(axis=axis)
        )
        part_2 = (
            cls.params["y8"] * cations["Si"] * cations["Al"]
            + cls.params["y9"] * cations["Si"] * cations["Mg"]
        )

        return part_1 + part_2

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


class Sun2024(Fe3Fe2_model):
    """
    Sun and Yao (2024)\ [11]_
    """

    components = ["FeO", "SiO2", "Al2O3", "TiO2", "CaO", "MgO", "Na2O", "K2O"]

    params = pd.Series(
        {
            "a0": 2.1479,
            "a1": -230.2593,
            "a2": -1.8557e-4,
            "a3": 34.3293,
            "a4": 1.4138,
            "a5": -17.3040,
            "a6": -10.1820,
            "a7": -6.7463,
            "a8": -7.3886,
            "a9": -14.5430,
            "a10": -9.9776,
            "a11": -16.1506,
            "a12": -37.5572,
            "h": 2.1410,
        }
    )
    # TODO calculate parameters
    error_params = [1, 1, 1, 1]

    @classmethod
    def _Gamma(
        cls,
        T_K: float | np.ndarray,
        P_bar: float | np.ndarray,
        **kwargs,
    ):
        """
        Calculate dV according to Deng (2020)
        """
        return Deng2020._dVdP(T_K=T_K, P_bar=P_bar, melt_Fe="12.5molpc", **kwargs) / (
            R * T_K
        )

    @classmethod
    def _Phi(cls, cations):

        axis = [0, 1][isinstance(cations, pd.DataFrame)]

        return (
            cls.params["a4"] * np.log(cations["Fe"])
            + cls.params["a5"] * cations["Fe"] ** 0.5
            + cls.params["a6"] * cations["Si"] ** 3.0
            + cations[["Al", "Ti", "Ca", "Mg"]]
            .mul(cls.params[["a7", "a8", "a9", "a10"]].values)
            .sum(axis=axis)
            + (cls.params["a11"] + cls.params["a12"] * cations["Fe"])
            * (cations["Na"] + cations["K"])
        )

    @classmethod
    def _Omega(cls, T_K):

        return (
            cls.params["a1"]
            + cls.params["a2"] * T_K**1.5
            + cls.params["a3"] * np.log(T_K)
        )

    @classmethod
    def calculate_Fe3Fe2(cls, melt_mol_fractions, T_K, P_bar, fO2, *args, **kwargs):
        """
        Calculate melt |Fe3Fe2| ratios with equation 9.

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar : float, array-like
            pressure in bar
        fO2 :   float, array-like
            Oxygen fugacity
        Fe3Fe2 : float, array-like
            melt Fe3+/Fe2+ ratios

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        moles = check_components(
            composition=melt_mol_fractions, components=cls.components
        )
        cations = moles.cations()

        omega = cls._Omega(T_K=T_K)
        phi = cls._Phi(cations=cations)
        gamma = cls._Gamma(T_K=T_K, P_bar=P_bar)

        return 10 ** (
            (np.log10(fO2) - omega - phi - cls.params["h"] * gamma)
            / (4 + cls.params["a0"] * cations["Fe"] ** 0.5)
        )

    @classmethod
    def get_error(cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs):

        return super().get_error(
            Fe3Fe2=Fe3Fe2,
            error_params_1bar=error_params_1bar,
            error_params_high_pressure=error_params_high_pressure,
            pressure=pressure,
        )


_clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
# Collect all Fe3Fe2_models in a dictionary.
Fe3Fe2_models_dict = {cls[0]: cls[1] for cls in _clsmembers if _is_Fe3Fe2_model(cls[1])}
