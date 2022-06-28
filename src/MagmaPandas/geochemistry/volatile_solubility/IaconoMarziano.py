from ...parse.validate import _check_setter, _check_argument
import pandas as pd
import numpy as np
from scipy.optimize import root_scalar, root

"""
Equations from:

G. Iacono-Marziano, Y. Morizet, E. Le Trong & F. Gaillard (2012)
New experimental data and semi-empirical parameterization of H2O–CO2 solubility in mafic melts.
Geochim. Gosmochim. Actca 97: 1-23

Python code modified from VESIcal:
https://github.com/kaylai/VESIcal
"""


def calculate_saturation(oxide_wtPercents, T_K, **kwargs):
    """
    Docstring
    """

    model = IaconoMarziano_configuration().model
    equation = globals()[model].calculate_saturation
    

    return equation(oxide_wtPercents, T_K, **kwargs)


def calculate_solubility(oxide_wtPercents, P_bar, T_K, **kwargs):
    """
    Docstring
    """

    model = IaconoMarziano_configuration().model
    equation = globals()[model].calculate_solubility

    return equation(oxide_wtPercents, P_bar, T_K, **kwargs)


parameter_options = ["hydrous_webapp", "hydrous_manuscript", "anhydrous"]
fugacity_options = ["ideal"]
activity_options = ["ideal"]
model_options = ["mixed", "h20", "co2"]


class IaconoMarziano_configuration:
    """
    Configure settings for the Iacono-Marziano volatile solubility model
    """

    __parameters = "hydrous_webapp"
    __fugacity = "ideal"
    __activity = "ideal"
    __model = "mixed"

    @property
    def parameters(cls):
        return cls.__parameters

    @parameters.setter
    @_check_setter(parameter_options)
    def parameters(cls, model: str):
        cls.__parameters = model

    @property
    def fugacity(cls):
        return cls.__fugacity

    @fugacity.setter
    @_check_setter(fugacity_options)
    def fugacity(cls, model: str):
        cls.__fugacity = model

    @property
    def activity(cls):
        return cls.__activity

    @activity.setter
    @_check_setter(activity_options)
    def activity(cls, model: str):
        cls.__activity = model

    @property
    def model(cls):
        return cls.__model

    @model.setter
    @_check_setter(model_options)
    def model(cls, model: str):
        cls.__model = model

    @classmethod
    def reset(cls):
        cls.__parameters = "hydrous_webapp"
        cls.__fugacity = "ideal"
        cls.__activity = "ideal"
        cls.__model = "mixed"

    @classmethod
    def print(cls):
        """ 
        
        """

        variables = {
            "Parameterisation": "__parameters",
            "Fugacity model": "__fugacity",
            "Activity model": "__activity",
            "Species model": "__model",
        }

        pad_left = 20
        pad_right = 20
        pad_total = pad_left + pad_right

        print(" Iacono-Marziano volatile solubility ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("Settings".ljust(pad_total, "_"))
        for param, model in variables.items():
            model_attr = f"_IaconoMarziano_configuration{model}"
            print(
                f"{param:.<{pad_left}}{getattr(cls, model_attr):.>{pad_right}}"
            )
        print("\nCalibration range".ljust(pad_total, "_"))
        T_string = f"1373-1673\N{DEGREE SIGN}K"
        print(f"{'Temperature':.<{pad_left}}{T_string:.>{pad_right}}")
        print(f"{'Pressure':.<{pad_left}}{'0.1-5 kbar':.>{pad_right}}")


H2O_coefficients = {
    "hydrous_webapp": {
        "a": 0.52096846,
        "b": 2.11575907,
        "B": -3.24443335,
        "C": -0.02238884,
    },
    "hydrous_manuscript": {
        "a": 0.53,
        "b": 2.35,
        "B": -3.37,
        "C": -0.02,
    },
    "anhydrous": {
        "a": 0.54,
        "b": 1.24,
        "B": -2.95,
        "C": 0.02,
    },
}

CO2_coefficients = {
    "hydrous": {
        "d_H2O": -16.4,
        "d_AI": 4.4,
        "d_FM": -17.1,
        "d_NK": 22.8,
        "a": 1.0,
        "b": 17.3,
        "B": -6.0,
        "C": 0.12,
    },
    "anhydrous": {
        "d_H2O": 2.3,
        "d_AI": 3.8,
        "d_FM": -16.3,
        "d_NK": 20.1,
        "a": 1.0,
        "b": 15.8,
        "B": -5.3,
        "C": 0.14,
    },
}

# Lower case names for classes, as to not mix up with variable names
class h2o:
    @staticmethod
    def calculate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid=1):
        """
        Solve quation 13

        x_fluid : float, int
            fluid mol fraction H2O
        """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        if any(i <= 0 for i in [P_bar, x_fluid]):
            return 0

        composition = oxide_wtPercents.copy()

        if IaconoMarziano_configuration().parameters == "anhydrous":
            # Solve equation 13
            return h2o._solubility(composition, x_fluid, P_bar, T_K)
        else:
            # Find equilibrium H2O content
            return root_scalar(
                h2o._solubility_rootFunction,
                args=(composition, x_fluid, P_bar, T_K),
                x0=1.0,
                x1=2.0,
            ).root

    @staticmethod
    def calculate_saturation(oxide_wtPercents, T_K, **kwargs):
        """ """
        if "H2O" not in oxide_wtPercents.index:
            raise ValueError("H2O not found in sample")
        if oxide_wtPercents["H2O"] <= 0:
            raise ValueError(f"H2O lower than 0: {oxide_wtPercents['H2O']}")

        composition = oxide_wtPercents.copy()

        # Upper limit of 15kbar
        upper_limit = 1.5e4
        try:
            P_saturation = root_scalar(
                h2o._saturation_rootFunction,
                args=(composition, T_K, kwargs),
                bracket=[1e-15, upper_limit],
            ).root
        except:
            P_saturation = np.nan
        return P_saturation

    @staticmethod
    def _solubility(oxide_wtPercents, x_fluid, P_bar, T_K):
        """
        Equation 13 from Iacono-Marziano et al. (2012)
        """

        coefficients = H2O_coefficients[IaconoMarziano_configuration().parameters]
        a, b, B, C = [coefficients[key] for key in ["a", "b", "B", "C"]]

        mol_fractions = oxide_wtPercents.moles
        NBO_O = NBO_O_calculate(mol_fractions)
        

        if IaconoMarziano_configuration().fugacity == "ideal":
            P_H2O = x_fluid * P_bar

        H2O = np.exp(a * np.log(P_H2O) + b * NBO_O + B + C * P_bar / T_K)

        return H2O

    @staticmethod
    def _solubility_rootFunction(H2O, oxide_wtPercents, x_fluid, P_bar, T_K):
        """
        Compare input and output dissolved H2O
        """
        # Copy so pandas doesn't raise SettingWithCopyWarning
        composition = oxide_wtPercents.copy()
        composition["H2O"] = H2O

        return H2O - h2o._solubility(composition, x_fluid, P_bar, T_K)

    @staticmethod
    def _saturation_rootFunction(P_bar, oxide_wtPercents, T_K, kwargs):
        """ 
        
        """
        composition = oxide_wtPercents.copy()
        #
        return composition["H2O"] - h2o.calculate_solubility(
            oxide_wtPercents=composition, P_bar=P_bar, T_K=T_K, **kwargs
        )


# Lower case names for classes, as to not mix up with variable names
class co2:
    @staticmethod
    def calculate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid=0.0):
        """
        Equation 12 from Iacono-Marziano (2012)
        """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        if any(i <= 0 for i in [P_bar, (1 - x_fluid)]):
            return 0

        if "anhydrous" in IaconoMarziano_configuration().parameters:
            parameterisation = "anhydrous"
        else:
            parameterisation = "hydrous"

        coefficients = CO2_coefficients[parameterisation]
        d_H2O, d_AI, d_FM, d_NK, a, b, B, C = [
            coefficients[key]
            for key in ["d_H2O", "d_AI", "d_FM", "d_NK", "a", "b", "B", "C"]
        ]

        composition = oxide_wtPercents.copy()
        composition["H2O"] = h2o.calculate_solubility(composition, P_bar, T_K, x_fluid)
        mol_fractions = composition.moles
        NBO_O = NBO_O_calculate(mol_fractions)

        if "Fe2O3" in mol_fractions.index:
            Fe2O3 = mol_fractions["Fe2O3"]
        else:
            Fe2O3 = 0.0

        if IaconoMarziano_configuration().fugacity == "ideal":
            P_CO2 = (1 - x_fluid) * P_bar

        x_AI = mol_fractions["Al2O3"] / (
            mol_fractions["CaO"] + mol_fractions["K2O"] + mol_fractions["Na2O"]
        )
        x_FM = mol_fractions["FeO"] + mol_fractions["MgO"] + 2 * Fe2O3
        x_NK = mol_fractions["Na2O"] + mol_fractions["K2O"]

        CO3_ppm = np.exp(
            mol_fractions["H2O"] * d_H2O
            + x_AI * d_AI
            + x_FM * d_FM
            + x_NK * d_NK
            + a * np.log(P_CO2)
            + b * NBO_O
            + B
            + C * P_bar / T_K
        )

        return CO3_ppm / 1e4

    @staticmethod
    def calculate_saturation(oxide_wtPercents, T_K, **kwargs):
        """ """
        if "CO2" not in oxide_wtPercents.index:
            raise ValueError("CO2 not found in sample")
        if oxide_wtPercents["CO2"] <= 0:
            raise ValueError(f"CO2 lower than 0: {oxide_wtPercents['CO2']}")
        # Upper limit of 100kbar

        composition = oxide_wtPercents.copy()
        upper_limit = 1e5
        try:
            P_saturation = root_scalar(
                co2._saturation_rootFunction,
                args=(composition, T_K, kwargs),
                bracket=[1e-10, upper_limit],
            ).root
        except:
            P_saturation = np.nan
        return P_saturation

    @staticmethod
    def _saturation_rootFunction(P_bar, oxide_wtPercents, T_K, kwargs):
        """
        Compare calculated and sample CO2
        """
        composition = oxide_wtPercents.copy()
        return composition["CO2"] - co2.calculate_solubility(
            oxide_wtPercents=composition, P_bar=P_bar, T_K=T_K, **kwargs
        )


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
            x0=[P_guess, 0.],
            args=(composition, T_K),
        ).x

        return_dict = {"P": saturation[0], "x_fluid": saturation[1], "both": saturation}
        return return_dict[output]

    @staticmethod
    @_check_argument("output", [None, "both", "CO2", "H2O"])
    def calculate_solubility(oxide_wtPercents, P_bar, T_K, x_fluid, output="both"):

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")

        composition = oxide_wtPercents.copy()
        H2O = h2o.calculate_solubility(composition, P_bar, T_K, x_fluid)
        composition["H2O"] = H2O
        CO2 = co2.calculate_solubility(composition, P_bar, T_K, x_fluid)

        return_dict = {"both": (H2O, CO2), "H2O": H2O, "CO2": CO2}
        return return_dict[output]

    @staticmethod
    def _saturation_rootFunction(P_x_fluid, oxide_wtPercents, T_K):
        """
        compare calculated with sample concentrations
        """

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


def NBO_O_calculate(mol_fractions):
    """
    Non-bridging oxygen / oxygen ratio, following Marrochhi & Toplis (2005).
    See also Iacono-Marziano (2012) appendix.
    """
    if isinstance(mol_fractions, pd.Series):
        elements = mol_fractions.index
    elif isinstance(mol_fractions, pd.DataFrame):
        elements = mol_fractions.columns

    # All Fe is considered as FeO
    if "Fe2O3" in elements:
        Fe2O3 = mol_fractions["Fe2O3"]
    else:
        Fe2O3 = 0.0

    NBO = 2 * (
        mol_fractions["K2O"]
        + mol_fractions["Na2O"]
        + mol_fractions["CaO"]
        + mol_fractions["MgO"]
        + mol_fractions["FeO"]
        + 2 * Fe2O3
        - mol_fractions["Al2O3"]
    )
    O = (
        2 * mol_fractions["SiO2"]
        + 2 * mol_fractions["TiO2"]
        + 3 * mol_fractions["Al2O3"]
        + mol_fractions["MgO"]
        + mol_fractions["FeO"]
        + 2 * Fe2O3
        + mol_fractions["CaO"]
        + mol_fractions["Na2O"]
        + mol_fractions["K2O"]
    )

    if not IaconoMarziano_configuration().parameters == "anhydrous":
        NBO = NBO + 2 * mol_fractions["H2O"]
        O = O + mol_fractions["H2O"]

    return NBO / O
