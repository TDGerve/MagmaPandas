"""
|CO2|-|H2O| solubility models from Iacono-Marziano et al. (2012)\ [33]_
"""

from typing import Tuple

import numpy as np
import pandas as pd
from scipy.optimize import root, root_scalar

from MagmaPandas.magma_protocol import Magma
from MagmaPandas.parse_io.validate import _check_argument, _check_setter
from MagmaPandas.volatile_solubility.solubility_baseclass import Solubility_model

"""
Equations from:

Iacono-Marziano G., Morizet Y., Le Trong E., Gaillard F. (2012) New experimental data and semi-empirical parameterization of |H2O|–|CO2| solubility in mafic melts. Geochimica et cosmochimica Acta. 97.
"""

parameter_options = ["hydrous_webapp", "hydrous_manuscript", "anhydrous"]
fugacity_options = ["ideal"]
activity_options = ["ideal"]
model_options = ["mixed", "h2o", "co2"]


class _meta_IaconoMarziano_configuration(type):
    def __init__(cls, *args, **kwargs):
        cls._parameters = "hydrous_webapp"
        cls._fugacity = "ideal"
        cls._activity = "ideal"
        cls._model = "mixed"

    @property
    def parameters(cls):
        return cls._parameters

    @parameters.setter
    @_check_setter(parameter_options)
    def parameters(cls, model: str):
        cls._parameters = model

    @property
    def fugacity(cls):
        return cls._fugacity

    @fugacity.setter
    @_check_setter(fugacity_options)
    def fugacity(cls, model: str):
        cls._fugacity = model

    @property
    def activity(cls):
        return cls._activity

    @activity.setter
    @_check_setter(activity_options)
    def activity(cls, model: str):
        cls._activity = model

    @property
    def model(cls):
        return cls._model

    @model.setter
    @_check_setter(model_options)
    def model(cls, model: str):
        cls._model = model

    @classmethod
    def reset(cls):
        cls._parameters = "hydrous_webapp"
        cls._fugacity = "ideal"
        cls._activity = "ideal"
        cls._model = "mixed"


class IaconoMarziano_configuration(metaclass=_meta_IaconoMarziano_configuration):
    @classmethod
    def reset(cls):
        cls._parameters = "hydrous_webapp"
        cls._fugacity = "ideal"
        cls._activity = "ideal"
        cls._model = "mixed"

    @classmethod
    def print(cls):
        """ """

        variables = {
            "Parameterisation": "_parameters",
            "Fugacity model": "_fugacity",
            "Activity model": "_activity",
            "Species model": "_model",
        }

        pad_left = 20
        pad_right = 20
        pad_total = pad_left + pad_right

        print(" Iacono-Marziano volatile solubility ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("Settings".ljust(pad_total, "_"))
        for param, model in variables.items():
            # model_attr = f"_IaconoMarziano_configuration{model}"
            print(f"{param:.<{pad_left}}{getattr(cls, model):.>{pad_right}}")
        print("\nCalibration range".ljust(pad_total, "_"))
        T_string = f"1373-1673\N{DEGREE SIGN}K"
        print(f"{'Temperature':.<{pad_left}}{T_string:.>{pad_right}}")
        print(f"{'Pressure':.<{pad_left}}{'0.1-5 kbar':.>{pad_right}}")
        print("\n")


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
class h2o(Solubility_model):
    @staticmethod
    def calculate_solubility(
        oxide_wtPercents: Magma,
        P_bar: float,
        T_K: float,
        x_fluid: float = 1.0,
        **kwargs,
    ) -> float:
        """
        Calculate melt |H2O| solubility according to equation 13.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %.
        P_bar : float
            Pressure in bar.
        T_K : float
            temperature in Kelvin
        x_fluid : float
            fraction of |H2O| in the fluid. Default value is 1.0

        Returns
        -------
        H2O : float
            melt |H2O| solubility in wt.%.
        """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")
        if P_bar < 0:
            raise ValueError(f"Pressure is negative: '{P_bar}'")

        if any(i <= 0 for i in [P_bar, x_fluid]):
            return 0

        composition = oxide_wtPercents.copy()

        if IaconoMarziano_configuration.parameters == "anhydrous":
            # Solve equation 13
            return h2o._solubility(composition, x_fluid, P_bar, T_K)
        else:
            # Find equilibrium H2O content
            return np.float32(
                root_scalar(
                    h2o._solubility_rootFunction,
                    args=(composition, x_fluid, P_bar, T_K),
                    x0=1.0,
                    x1=2.0,
                ).root
            )

    @staticmethod
    def calculate_saturation(oxide_wtPercents: Magma, T_K: float, **kwargs) -> float:
        """
        Calculate melt |H2O| saturation pressure according to equation 13.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %. Needs to have an 'H2O' column.
        T_K : float
            temperature in Kelvin
        kwargs
            keyword arguments passed to :py:meth:`~MagmaPandas.volatile_solubility.models.IaconoMarziano.h2o.calculate_solubility`

        Returns
        -------
        P_saturation : float
            Saturation pressure in bar
        """

        if "H2O" not in oxide_wtPercents.index:
            raise ValueError("H2O not found in sample")
        if oxide_wtPercents["H2O"] < 0:
            raise ValueError(f"H2O lower than 0: {oxide_wtPercents['H2O']}")
        if oxide_wtPercents["H2O"] == 0.0:
            return 0.0

        composition = oxide_wtPercents.copy()

        # Upper limit of 15kbar
        upper_limit = 1.5e4
        # try:
        P_saturation = root_scalar(
            h2o._saturation_rootFunction,
            args=(composition, T_K, kwargs),
            bracket=[1e-15, upper_limit],
        ).root
        # except:
        #     P_saturation = np.nan
        return P_saturation

    @staticmethod
    def _solubility(oxide_wtPercents: Magma, x_fluid, P_bar, T_K):
        """
        Equation 13 from Iacono-Marziano et al. (2012)
        """

        coefficients = H2O_coefficients[IaconoMarziano_configuration.parameters]
        a, b, B, C = [coefficients[key] for key in ["a", "b", "B", "C"]]

        mol_fractions = oxide_wtPercents.moles()
        NBO_O = NBO_O_calculate(mol_fractions)

        if IaconoMarziano_configuration.fugacity == "ideal":
            P_H2O = x_fluid * P_bar

        H2O = np.exp(a * np.log(P_H2O) + b * NBO_O + B + C * P_bar / T_K)

        return H2O

    @staticmethod
    def _solubility_rootFunction(H2O, oxide_wtPercents: Magma, x_fluid, P_bar, T_K):
        """
        Compare input and output dissolved H2O
        """
        # Copy so pandas doesn't raise SettingWithCopyWarning
        composition = oxide_wtPercents.copy()
        composition["H2O"] = np.float32(H2O)

        return H2O - h2o._solubility(composition, x_fluid, P_bar, T_K)

    @staticmethod
    def _saturation_rootFunction(P_bar, oxide_wtPercents: Magma, T_K, kwargs):
        """ """
        composition = oxide_wtPercents.copy()
        #
        return composition["H2O"] - h2o.calculate_solubility(
            oxide_wtPercents=composition,
            P_bar=P_bar,
            T_K=T_K,
            **kwargs,
        )


# Lower case names for classes, as to not mix up with variable names
class co2(Solubility_model):
    @staticmethod
    def calculate_solubility(
        oxide_wtPercents: Magma,
        P_bar: float | np.ndarray,
        T_K: float | np.ndarray,
        x_fluid: float = 0.0,
        **kwargs,
    ) -> float:
        """
        Calculate melt |CO2| solubility according to equation 12.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %.
        P_bar : float, array-like
            Pressure in bar
        T_K : float, array-like
            temperature in Kelvin
        x_fluid : float
            fraction of |H2O| in the fluid. Default value is 0.0

        Returns
        -------
        solublity : float
            melt |CO2| solubility in wt. %.
        """

        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")
        if P_bar < 0:
            raise ValueError(f"Pressure is negative: '{P_bar}'")

        if any(i <= 0 for i in [P_bar, (1 - x_fluid)]):
            return 0

        if "anhydrous" in IaconoMarziano_configuration.parameters:
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
        mol_fractions = composition.moles()
        NBO_O = NBO_O_calculate(mol_fractions)

        if "Fe2O3" in mol_fractions.index:
            Fe2O3 = mol_fractions["Fe2O3"]
        else:
            Fe2O3 = 0.0

        if IaconoMarziano_configuration.fugacity == "ideal":
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
    def calculate_saturation(
        oxide_wtPercents: Magma,
        T_K: float,
        **kwargs,
    ) -> float:
        """
        Calculate melt |CO2| saturation pressure according to equation 12.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %. Needs to have a 'CO2' column.
        T_K : float
            temperature in Kelvin
        kwargs
            keyword arguments passed to :py:meth:`~MagmaPandas.volatile_solubility.models.IaconoMarziano.co2.calculate_solubility`

        Returns
        -------
        P_saturation : float
            Saturation pressure in bar
        """

        if "CO2" not in oxide_wtPercents.index:
            raise ValueError("CO2 not found in sample")
        if oxide_wtPercents["CO2"] < 0:
            raise ValueError(f"CO2 lower than 0: {oxide_wtPercents['CO2']}")
        if oxide_wtPercents["CO2"] == 0.0:
            return 0.0

        composition = oxide_wtPercents.copy()
        # Upper limit of 100kbar
        upper_limit = 1e5
        # try:
        P_saturation = root_scalar(
            co2._saturation_rootFunction,
            args=(composition, T_K, kwargs),
            bracket=[1e-10, upper_limit],
        ).root
        # except:
        #     P_saturation = np.nan
        return P_saturation

    @staticmethod
    def _saturation_rootFunction(P_bar, oxide_wtPercents: Magma, T_K, kwargs):
        """
        Compare calculated and sample CO2
        """
        composition = oxide_wtPercents.copy()
        return composition["CO2"] - co2.calculate_solubility(
            oxide_wtPercents=composition,
            P_bar=P_bar,
            T_K=T_K,
            **kwargs,
        )


class mixed(Solubility_model):
    @staticmethod
    @_check_argument("output", [None, "PXfl", "P", "Xfl"])
    def calculate_saturation(
        oxide_wtPercents: Magma, T_K: float, output: str = "P", **kwargs
    ) -> float | Tuple[float, float]:
        """
        Calculate volatile saturation pressure for systems with mixed |CO2|-|H2O| fluids.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composision in oxide wt. %. Needs to have 'H2O' and 'CO2' columns.
        T_K : float
            Temperature in kelvin
        output : str
            Output format. 'P' for pressure only, 'Xfl' for |H2O| fluid fraction only and 'PXfl' for both.

        Returns
        -------
        saturation : float, (float, float)
            Depending on the value of ``output``: saturation pressure in bar, |H2O| fluid fraction or (saturation pressure, fluid fraction)
        """
        composition = oxide_wtPercents.copy()

        P_H2O_saturation = h2o.calculate_saturation(composition, T_K=T_K, x_fluid=1.0)
        P_CO2_saturation = co2.calculate_saturation(composition, T_K=T_K, x_fluid=0.0)

        if oxide_wtPercents["H2O"] <= 0.0:
            return P_CO2_saturation
        if oxide_wtPercents["CO2"] <= 0.0:
            return P_H2O_saturation

        P_guess = 0

        for species in (P_H2O_saturation, P_CO2_saturation):
            if np.isfinite(species):
                P_guess += species

        saturation = root(
            mixed._saturation_rootFunction,
            x0=[P_guess, 0.0],
            args=(composition, T_K),
        ).x

        if saturation[1] <= 0.0:
            saturation[0] = P_CO2_saturation
        elif saturation[1] >= 1.0:
            saturation[0] = P_H2O_saturation
        saturation[1] = np.clip(saturation[1], 0.0, 1.0)

        return_dict = {"P": saturation[0], "Xfl": saturation[1], "PXfl": saturation}
        return return_dict[output]

    @staticmethod
    @_check_argument("output", [None, "both", "CO2", "H2O"])
    def calculate_solubility(
        oxide_wtPercents: Magma,
        P_bar: float,
        T_K: float,
        x_fluid: float,
        output: str = "both",
        **kwargs,
    ) -> float | Tuple[float, float]:
        """
        Calculate volatile solubilities for systems with mixed |CO2|-|H2O| fluids.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composision in oxide wt. %.
        P_bar : float
            pressure in bar
        T_K : float
            Temperature in kelvin
        x_fluid: float
            fraction of |H2O| in the fluid.
        output : str
            Output format. 'CO2' for |CO2| only, 'H2O' for |H2O| only and 'both' for both.

        Returns
        -------
        saturation : float, (float, float)
            Solubility in wt. %. Depending on the value of ``output``: |CO2|, |H2O| or (|CO2|, |H2O|).
        """
        if not 1 >= x_fluid >= 0:
            raise ValueError(f"x_fluid: {x_fluid} is not between 0 and 1")
        if P_bar < 0:
            raise ValueError(f"Pressure is negative: '{P_bar}'")

        composition = oxide_wtPercents.copy()
        H2O = h2o.calculate_solubility(composition, P_bar, T_K, x_fluid)
        composition["H2O"] = np.float32(H2O)
        CO2 = co2.calculate_solubility(composition, P_bar, T_K, x_fluid)

        return_dict = {"both": (H2O, CO2), "H2O": H2O, "CO2": CO2}
        return return_dict[output]

    @staticmethod
    def _saturation_rootFunction(P_x_fluid, oxide_wtPercents: Magma, T_K):
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

    if not IaconoMarziano_configuration.parameters == "anhydrous":
        NBO = NBO + 2 * mol_fractions["H2O"]
        O = O + mol_fractions["H2O"]

    return NBO / O
