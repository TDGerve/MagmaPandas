from typing import List
import numpy as np
from scipy.optimize import root_scalar, root
from scipy.constants import R
from MagmaPandas.parse.validate import _check_argument, _check_setter
from MagmaPandas.geochemistry.eos_volatiles import hollowayBlank
from elements.elements import compound_weights

"""
Equations from:

Allison, C. M., Roggensack, K., Clarke, A. B. (2022)
MafiCH: a general model for H2O–CO2 solubility in mafic magmas
Contributions to mineralogy and petrology 40

"""


def calculate_saturation(*args, **kwargs):
    """
    Docstring
    """
    model = Allison_configuration.model
    equation = globals()[model].calculate_saturation

    return equation(*args, **kwargs)


def calculate_solubility(*args, **kwargs):
    """
    Docstring
    """
    model = Allison_configuration.model
    equation = globals()[model].calculate_solubility

    return equation(*args, **kwargs)


fugacity_options = ["hollowayBlank"]
model_options = ["mixed", "h2o", "co2"]


class _meta_Allison_configuration(type):
    def __init__(cls, *args, **kwargs):
        cls._fugacity = "hollowayBlank"
        cls._model = "mixed"

    @property
    def fugacity(cls):
        return cls._fugacity

    @fugacity.setter
    @_check_setter(fugacity_options)
    def fugacity(cls, model: str):
        cls._fugacity = model

    @property
    def model(cls):
        return cls._model

    @model.setter
    @_check_setter(model_options)
    def model(cls, model: str):
        cls._model = model


class Allison_configuration(metaclass=_meta_Allison_configuration):
    @classmethod
    def reset(cls):
        cls._fugacity = "hollowayBlank"
        cls._model = "mixed"

    @classmethod
    def print(cls):
        """ """

        variables = {"Fugacity model": "_fugacity", "Species model": "_model"}
        pad_left = 20
        pad_right = 20
        pad_total = pad_left + pad_right
        print(" Allison (2022) volatile solubility ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("Settings".ljust(pad_total, "_"))
        for param, model in variables.items():
            print(f"{param:.<{pad_left}}{getattr(cls, model):.>{pad_right}}")
        print("\nCalibration range".ljust(pad_total, "_"))
        T_string = f"{1000+273.15:.0f}-{1400+273.14:.0f}\N{DEGREE SIGN}K"
        print(f"{'Temperature':.<{pad_left}}{T_string:.>{pad_right}}")
        print(f"{'Pressure':.<{pad_left}}{'< 7 kbar':.>{pad_right}}")


FeO_mass, Fe2O3_mass = compound_weights(["FeO", "Fe2O3"])
fugacity_model = hollowayBlank.fugacity


class h2o:
    @staticmethod
    def calculate_saturation(oxide_wtPercents, T_K, x_fluid=1.0):
        """
        Equation 8 from Allison 2022
        """

        if oxide_wtPercents["H2O"] <= 0:
            raise ValueError(f"H2O lower than 0: {oxide_wtPercents['H2O']}")
        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()

        H2O = composition["H2O"]
        fH2O = 104.98 * H2O**1.83
        fH2O_pure = fH2O / x_fluid

        P_saturation = root_scalar(
            _root_fugacity_pressure,
            args=(T_K, fH2O_pure, "H2O"),
            bracket=[1e-15, 1.5e4],
        ).root

        return P_saturation

    @staticmethod
    def calculate_solubility(P_bar, T_K, x_fluid=1.0):
        """
        Equation 8 from Allison 2022
        """
        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")
        if P_bar < 0:
            raise ValueError(f"Pressure is negative: '{P_bar}'")

        fH2O_pure = fugacity_model(T_K, P_bar, "H2O")
        fH2O = fH2O_pure * x_fluid

        return (fH2O / 104.98) ** (1 / 1.83)


class co2:
    @staticmethod
    def calculate_saturation(oxide_wtPercents, T_K, x_fluid=0.0):

        if oxide_wtPercents["CO2"] <= 0:
            raise ValueError(f"CO2 lower than 0: {oxide_wtPercents['CO2']}")
        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()
        CO2 = composition["CO2"]
        cations = co2._cation_fractions_Allison(composition)
        deltaV = co2._deltaV(cations)
        lnK0 = co2._lnK0(cations)
        FW = 36.594  # alkali basalt formula weight per 1 oxygen

        XCO3 = CO2 * (1 / 44.01) / ((100 / FW) - (CO2 / FW))
        Kf = XCO3 / (1 + XCO3)

        # Partial pressure of CO2
        P_CO2 = root_scalar(
            co2._root_partial_pressure,
            args=(T_K, deltaV, lnK0, Kf),
            bracket=[1e-15, 1.5e4],
        ).root
        if x_fluid <= 0:
            return P_CO2
        else:
            fCO2 = fugacity_model(T_K, P_CO2, "CO2")
            fCO2pure = fCO2 / (1 - x_fluid)

            # Saturation pressure
            P_saturation = root_scalar(
                _root_fugacity_pressure,
                args=(T_K, fCO2pure, "CO2"),
                bracket=[1e-15, 1.5e4],
            ).root

            return P_saturation

    @staticmethod
    def calculate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid=0.0):
        """ """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")
        if P_bar < 0:
            raise ValueError(f"Pressure is negative: '{P_bar}'")

        composition = oxide_wtPercents.copy()

        Ra = R * 10  # cm3.bar.K-1.mol-1
        FW = 36.594  # alkali basalt formula weight per 1 oxygen
        P0 = 1e3  # reference pressure in bars

        # CO2 fugacity of a pure fluid
        fCO2_pure = fugacity_model(T_K, P_bar, "CO2")
        # CO2 fugacity of the mixed fluid
        fCO2 = fCO2_pure * (1 - x_fluid)
        # Partial pressure of CO2
        if x_fluid > 0:
            P_CO2 = root_scalar(
                _root_fugacity_pressure, args=(T_K, fCO2, "CO2"), bracket=[1e-15, 1.5e4]
            ).root
        else:
            P_CO2 = P_bar

        cations = co2._cation_fractions_Allison(composition)
        deltaV = co2._deltaV(cations)
        lnK0 = co2._lnK0(cations)

        K = np.exp(lnK0) * np.exp(-deltaV * (P_CO2 - P0) / (Ra * T_K))
        Kf = K * fCO2
        XCO3 = Kf / (1 - Kf)
        CO2 = 44.01 * XCO3 / (44.01 * XCO3 + (1 - XCO3) * FW) * 100

        return CO2

    @staticmethod
    def _root_partial_pressure(P_bar, T_K, deltaV, lnK0, Kf):

        P0 = 1e3  # reference pressure
        Ra = R * 10  # Gas constant in cm3.bar.K-1.mol-1

        K_fugacity = Kf / fugacity_model(T_K, P_bar, "CO2")
        K_solubility = np.exp(lnK0) * np.exp(-deltaV * (P_bar - P0) / (Ra * T_K))
        return K_fugacity - K_solubility

    @staticmethod
    def _deltaV(cations):
        """ """
        NaK = cations["Na"] / (cations["Na"] + cations["K"])

        deltaV = (
            -3350.65
            + 2625.385 * cations["Ti"]
            + 3105.426 * cations["Al"]
            + 47.0037 * NaK
            + 3375.552 * np.sum([cations[i] for i in ["Si", "Na"]], axis=0)
            + 3795.115 * cations["K"]
            + 3628.018 * cations["Fe"]
            + 3323.32 * np.sum([cations[i] for i in ["Mg", "Ca"]], axis=0)
        )
        return deltaV

    @staticmethod
    def _lnK0(cations):
        """ """
        NaK = cations["Na"] / (cations["Na"] + cations["K"])

        lnK0 = (
            -128.365
            + 122.644 * np.sum([cations[i] for i in ["Fe", "Na", "Ca"]], axis=0)
            + 92.263 * np.sum([cations[i] for i in ["Ti", "Al"]], axis=0)
            + 114.098 * cations["Si"]
            + 111.549 * cations["Mg"]
            + 138.855 * cations["K"]
            + 2.239 * NaK
        )
        return lnK0

    @staticmethod
    def _cation_fractions_Allison(oxide_wtPercents):
        """ """
        elements = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O"]

        composition = oxide_wtPercents.copy()

        # All Fe2O3 as Fe2+
        if "Fe2O3" in composition.index:
            mass_ratio = Fe2O3_mass / FeO_mass
            composition["FeO"] = composition["FeO"] + composition["Fe2O3"] / mass_ratio
        try:
            # For Series
            composition = composition.loc[elements]
        except KeyError:
            # For Dataframes
            composition = composition.loc[:, elements]
        composition.recalculate(inplace=True)

        cations = composition.cations
        # Rounding to 3 decimals because Allison did the same
        # Results will be different if you don't
        cations = cations.round(3)

        return cations


class mixed:
    @staticmethod
    @_check_argument("output", [None, "both", "P", "x_fluid"])
    def calculate_saturation(oxide_wtPercents, T_K, output="P"):

        composition = oxide_wtPercents.copy()

        P_H2O_saturation = h2o.calculate_saturation(composition, T_K=T_K, x_fluid=1.0)
        P_CO2_saturation = co2.calculate_saturation(composition, T_K=T_K, x_fluid=0.0)

        if oxide_wtPercents["H2O"] <= 0:
            return P_CO2_saturation
        if oxide_wtPercents["CO2"] <= 0:
            return P_H2O_saturation

        P_guess = 0

        for species in (P_H2O_saturation, P_CO2_saturation):
            if np.isfinite(species):
                P_guess += species

        saturation = root(
            mixed._saturation_rootFunction,
            x0=[P_guess, 0.01],
            args=(composition, T_K),
        ).x

        if saturation[1] <= 0.0:
            saturation[0] = P_CO2_saturation
        elif saturation[1] >= 1.0:
            saturation[0] = P_H2O_saturation
        saturation[1] = np.clip(saturation[1], 0.0, 1.0)

        return_dict = {"P": saturation[0], "x_fluid": saturation[1], "both": saturation}

        return return_dict[output]

    @staticmethod
    @_check_argument("output", [None, "both", "CO2", "H2O"])
    def calculate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid, output="both"):
        """ """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()

        H2O = h2o.calculate_solubility(P_bar, T_K, x_fluid)
        CO2 = co2.calculate_solubility(composition, P_bar, T_K, x_fluid)

        return_dict = {"both": np.array([H2O, CO2]), "H2O": H2O, "CO2": CO2}
        return return_dict[output]

    @staticmethod
    def _saturation_rootFunction(P_x_fluid: List, oxide_wtPercents, T_K):

        P_bar, x_fluid = P_x_fluid
        # Keep x_fluid and P_bar within bounds
        x_fluid = np.clip(x_fluid, 0.0, 1.0)
        P_bar = np.clip(P_bar, a_min=1e-15, a_max=None)

        composition = oxide_wtPercents.copy()

        H2O = composition["H2O"]
        CO2 = composition["CO2"]

        sample_concentrations = np.array([H2O, CO2])
        calculated_concentrations = np.array(
            mixed.calculate_solubility(
                oxide_wtPercents=composition,
                P_bar=P_bar,
                T_K=T_K,
                x_fluid=x_fluid,
            )
        )

        return abs(calculated_concentrations - sample_concentrations)


def _root_fugacity_pressure(P_bar, T_K, fCO2, species):
    return fCO2 - fugacity_model(T_K, P_bar, species)
