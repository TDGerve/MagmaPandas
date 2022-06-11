import numpy as np
from scipy.optimize import root_scalar, root
from scipy.constants import R
from MagmaPandas.parse.validate import _check_argument, _check_setter
from MagmaPandas.geochemistry.eos_volatiles import hollowayBlank
from elements.elements import compound_weights


def calculate_saturation(oxide_wtPercents, T_K, **kwargs):
    """
    Docstring
    """

    return 0


def calculate_solubility(oxide_wtPercents, P_bar, T_K, **kwargs):
    """
    Docstring
    """

    return 0


fugacity_options = ["hollowayBlank"]
model_options = ["mixed", "h20", "co2"]


class Allison_configuration:
    """
    Configure settings for the Allison (2022) volatile solubility model
    """

    __fugacity = "hollowayBlank"
    __model = "mixed"

    @property
    def fugacity(cls):
        return Allison_configuration.__fugacity

    @fugacity.setter
    @_check_setter(fugacity_options)
    def fugacity(cls, model: str):
        Allison_configuration.__fugacity = model

    @property
    def model(cls):
        return Allison_configuration.__model

    @model.setter
    @_check_setter(model_options)
    def model(cls, model: str):
        Allison_configuration.__model = model

    @staticmethod
    def reset():
        Allison_configuration.__fugacity = "hollowayBlank"
        Allison_configuration.__model = "mixed"

    @staticmethod
    def print():
        """ """
        names = [
            "Fugacity model",
            "Species model",
        ]
        attributes = [
            f"_Allison_configuration{i}"
            for i in ["__fugacity", "__model"]
        ]
        pad_left = 20
        pad_right = 20
        pad_total = pad_left + pad_right
        print("Settings".ljust(pad_total, "_"))
        for param, model in zip(names, attributes):
            print(
                f"{param:.<{pad_left}}{getattr(Allison_configuration, model):.>{pad_right}}"
            )
        print("\nCalibration range".ljust(pad_total, "_"))
        T_string = f"{1200+273.15}\N{DEGREE SIGN}K"
        print(f"{'Temperature':.<{pad_left}}{T_string:.>{pad_right}}")
        print(f"{'Pressure':.<{pad_left}}{'<7 kbar':.>{pad_right}}")


FeO_mass, Fe2O3_mass = compound_weights(["FeO", "Fe2O3"])
fugacity_model = hollowayBlank.fugacity


class h2o:
    @staticmethod
    def calculuate_saturation(oxide_wtPercents, T_K, x_fluid=1.0):
        """
        Equation 8 from Allison 2022
        """

        if oxide_wtPercents["H2O"] <= 0:
            raise ValueError(f"H2O lower than 0: {oxide_wtPercents['H2O']}")

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        H2O = oxide_wtPercents["H2O"]
        fH2O = 104.98 * H2O ** (1 / 1.83)
        fH2O_pure = fH2O / x_fluid

        P_saturation = root_scalar(
            _root_fugacity_pressure,
            args=(T_K, fH2O_pure, "H2O"),
            bracket=[1e-15, 1.5e4],
        )

        return P_saturation

    @staticmethod
    def calculuate_solubility(P_bar, T_K, x_fluid=1.0):
        """
        Equation 8 from Allison 2022
        """
        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        if any(i <= 0 for i in [P_bar, x_fluid]):
            return 0

        fH2O_pure = fugacity_model(T_K, P_bar, "H2O")
        fH2O = fH2O_pure * x_fluid

        return (fH2O / 104.98) ** (1 / 1.83)


class co2:
    @staticmethod
    def calculuate_saturation(oxide_wtPercents, T_K, x_fluid=0.0):

        if oxide_wtPercents["CO2"] <= 0:
            raise ValueError(f"CO2 lower than 0: {oxide_wtPercents['CO2']}")

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()
        CO2 = composition["CO2"]
        deltaV = co2._parameter(composition, "DeltaV")
        lnK0 = co2._parameter(composition, "lnK0")
        FW = 36.594  # alkali basalt formula weight per 1 oxygen

        XCO3 = CO2 * (1 / 44.01) / (100 / FW) - (CO2 / FW)
        Kf = XCO3 / (1 + XCO3)

        # Partial pressure of CO2
        P_CO2 = root_scalar(
            co2._root_partial_pressure,
            args=(T_K, deltaV, lnK0, Kf),
            bracket=[1e-15, 1.5e4],
        )
        if x_fluid == 0:

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
    def calculuate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid=0.0):
        """ """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        if any(i <= 0 for i in [P_bar, (1 - x_fluid)]):
            return 0

        Ra = R * 10  # cm3.bar.K-1.mol-1
        FW = 36.594  # alkali basalt formula weight per 1 oxygen
        P0 = 1e3  # reference pressure

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

        deltaV = co2._parameter(oxide_wtPercents, "DeltaV")
        lnK0 = co2._parameter(oxide_wtPercents, "lnK0")

        K = np.exp(lnK0) * np.exp(-deltaV * (P_CO2 - P0) / (Ra * T_K))
        Kf = K * fCO2
        XCO3 = Kf / (1 - Kf)
        CO2 = 44.01 * XCO3 * (44.01 * XCO3 + (1 - XCO3) * FW)

        return CO2

    @staticmethod
    def _root_partial_pressure(P_bar, T_K, deltaV, lnK0, Kf):

        P0 = 1e3  # reference pressure
        Ra = R * 10  # Gas constant in cm3.bar.K-1.mol-1

        K_fugacity = Kf / fugacity_model(T_K, P_bar, "CO2")
        K_solubility = np.exp(lnK0) * np.exp(-deltaV * (P_bar - P0) / (Ra * T_K))
        return K_fugacity - K_solubility

    @staticmethod
    @_check_argument("parameter", [None, "deltaV", "lnK0"])
    def _parameter(oxide_wtPercents, parameter):
        cations = co2._cation_fractions_Allison(oxide_wtPercents)
        NaK = cations["Na"] / (cations["Na"] + cations["K"])

        if parameter == "deltaV":
            deltaV = (
                -3350.65
                + 2625.385 * cations["Ti"]
                + 3105.426 * cations["Al"]
                + 47.0037 * NaK
                + 3375.552 * cations[["Si", "Na"]].sum()
                + 3795.115 * cations["K"]
                + 3628.018 * cations["Fe"]
                + 3323.32 * cations[["Mg", "Ca"]].sum()
            )
            return deltaV
        elif parameter == "lnK0":
            lnK0 = (
                -128.365
                + 122.644 * cations[["Fe", "Na", "Ca"]].sum()
                + 92.263 * cations["Ti", "Al"].sum()
                + 114.098 * cations["Si"]
                + 111.549 * cations["Mg"]
                + 138.855 * cations["K"]
                + 2.239 * NaK
            )
            return lnK0
        else:
            raise ValueError(
                f"parameter: {parameter} not recognised, please choose 'deltaV' or lnK0'"
            )

    @staticmethod
    def _cation_fractions_Allison(oxide_wtPercents):
        elements = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O"]

        if "Fe2O3" in oxide_wtPercents.index:
            mass_ratio = Fe2O3_mass / FeO_mass
            oxide_wtPercents["FeO"] = (
                oxide_wtPercents["FeO"] + oxide_wtPercents["Fe2O3"] / mass_ratio
            )

        cations = oxide_wtPercents.loc[elements].cations

        return cations


class mixed:
    @staticmethod
    def calculate_saturation(oxide_wtPercents, T_K, **kwargs):

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

        P_saturation = root(
            mixed._saturation_rootFunction,
            x0=[P_guess, 0.5],
            args=(composition, T_K, kwargs),
        ).x[0]

        return P_saturation

    @staticmethod
    def calculuate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid):
        """ """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()

        H2O = h2o.calculuate_solubility(P_bar, T_K, x_fluid)
        CO2 = co2.calculuate_solubility(composition, P_bar, T_K, x_fluid)

        return (H2O, CO2)

    @staticmethod
    def _saturation_rootFunction(P_x_fluid, oxide_wtPercents, T_K, kwargs):

        P_bar, x_fluid = P_x_fluid
        # Keep x_fluid and P_bar within bounds
        x_fluid = np.clip(x_fluid, 0.0, 1.0)
        P_bar = np.clip(P_bar, a_min=1e-15, a_max=None)

        H2O = oxide_wtPercents["H2O"]
        CO2 = oxide_wtPercents["CO2"]

        sample_concentrations = np.array([H2O, CO2])
        calculated_concentrations = np.array(
            mixed.calculuate_solubility(
                oxide_wtPercents=oxide_wtPercents,
                P_bar=P_bar,
                T_K=T_K,
                x_fluid=x_fluid,
                **kwargs,
            )
        )

        return calculated_concentrations - sample_concentrations


def _root_fugacity_pressure(P_bar, T_K, fCO2, species):
    return fCO2 - fugacity_model(T_K, P_bar, species)
