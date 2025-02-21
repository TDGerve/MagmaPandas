"""
|CO2|-|H2O| solubility models from Shiskina et al. (2010)\ [34]_
"""

from typing import Tuple

import numpy as np
from scipy.optimize import root, root_scalar

from MagmaPandas.MagmaSeries.MagmaSeries_baseclass import MagmaSeries
from MagmaPandas.parse_io.validate import _check_argument, _check_setter
from MagmaPandas.volatile_solubility.solubility_baseclass import Solubility_model

"""
Equations from:

Shiskina, T. A., Botcharnikov, R. E., Holtz, F., Almeev, R. R., Jazwa, A. M., Jakubiak, A. A. (2014)
Compositional and pressure effects on the solubility of H2Oand CO2 in maficmelt
Chemical geology 388 (112 - 129)

"""
fugacity_options = ["ideal"]


class _meta_shiskina_configuration(type):
    def __init__(cls, *args, **kwargs):
        cls._fugacity = "ideal"

    @property
    def fugacity(cls):
        return cls._fugacity

    @fugacity.setter
    @_check_setter(fugacity_options)
    def fugacity(cls, model: str):
        cls._fugacity = model

    @classmethod
    def reset(cls):
        cls._fugacity = "ideal"
        cls._model = "mixed"


class shiskina_configuration(metaclass=_meta_shiskina_configuration):
    @classmethod
    def reset(cls):
        cls._fugacity = "ideal"

    @classmethod
    def print(cls):
        """ """

        variables = {
            "Fugacity model": "_fugacity",
        }

        pad_left = 20
        pad_right = 20
        pad_total = pad_left + pad_right

        print(" Shiskina volatile solubility ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("Settings".ljust(pad_total, "_"))
        for param, model in variables.items():
            print(f"{param:.<{pad_left}}{getattr(cls, model):.>{pad_right}}")
        print("\nCalibration range".ljust(pad_total, "_"))
        T_string = f"1423-1523\N{DEGREE SIGN}K"
        print(f"{'Temperature':.<{pad_left}}{T_string:.>{pad_right}}")
        print(f"{'Pressure':.<{pad_left}}{'0.2-5 kbar':.>{pad_right}}")
        print("\n")


# Model parameters
co2_parameters = {
    "pi": {"A": 1.167, "B": 0.671, "C": 0.65},
    "pi_star": {"A": 1.150, "B": 6.71, "C": -1.345},
}


class h2o(Solubility_model):
    @staticmethod
    def calculate_saturation(oxide_wtPercents: MagmaSeries, **kwargs) -> float:
        """
        Calculate melt |H2O| solubility according to equation 9.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %.
        kwargs
            keyword arguments passed to :py:meth:`~MagmaPandas.volatile_solubility.models.Shiskina.h2o.calculate_solubility`

        Returns
        -------
        H2O : float
            melt |H2O| solubility in wt.%.
        """

        if "H2O" not in oxide_wtPercents.index:
            raise ValueError("H2O not found in sample")
        if oxide_wtPercents["H2O"] < 0:
            raise ValueError(f"H2O lower than 0: {oxide_wtPercents['H2O']}")
        if oxide_wtPercents["H2O"] == 0.0:
            return 0.0

        if oxide_wtPercents["H2O"] < h2o.calculate_solubility(
            oxide_wtPercents, P_bar=0
        ):
            return np.nan

        composition = oxide_wtPercents.copy()

        # Upper limit of 15kbar
        upper_limit = 1.5e4
        # try:
        P_saturation = root_scalar(
            h2o._root_saturation,
            args=(composition, kwargs),
            bracket=[1e-15, upper_limit],
        ).root
        # except:
        #     P_saturation = np.nan
        return P_saturation

    @staticmethod
    def calculate_solubility(
        oxide_wtPercents: MagmaSeries, P_bar: float, x_fluid: float = 1.0, **kwargs
    ) -> float:
        """
        Calculate melt |H2O| solubility according to equation 9.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %. Needs to have an 'H2O' column.
        P_bar : float
            Pressure in bar.
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

        oxides = {"K2O", "Na2O"}
        if oxides.difference(oxide_wtPercents.index):
            raise KeyError(
                f"missing oxides: {oxides.difference(oxide_wtPercents.index)}"
            )

        # Calculate cation mol fractions on an anhydrous basis
        volatiles = {"H2O", "CO2"}.intersection(oxide_wtPercents.index)
        if volatiles:
            oxide_wtPercents = oxide_wtPercents.drop(volatiles)
            oxide_wtPercents.recalculate(inplace=True)

        mol_fractions = oxide_wtPercents.cations()

        P_MPa = P_bar / 10

        # Partial pressure (ideal gas fugacity)
        fH2O = x_fluid * P_MPa

        a = 3.36e-7 * fH2O**3 - 2.33e-4 * fH2O**2 + 0.0711 * fH2O - 1.1309
        b = mol_fractions["Na"] + mol_fractions["K"]
        c = -1.2e-5 * fH2O**2 + 0.0196 * fH2O + 1.1297

        return a * b + c

    def _root_saturation(P_bar, oxide_wtPercents, kwargs):

        composition = oxide_wtPercents.copy()

        return composition["H2O"] - h2o.calculate_solubility(
            oxide_wtPercents=composition, P_bar=P_bar, **kwargs
        )


class co2(Solubility_model):
    @staticmethod
    def calculate_saturation(oxide_wtPercents: MagmaSeries, **kwargs) -> float:
        """
        Calculate melt |CO2| saturation pressure according to equation 13.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %. Needs to have a 'CO2' column.
        kwargs
            keyword arguments passed to :py:meth:`~MagmaPandas.volatile_solubility.models.Shiskina.co2.calculate_solubility`

        Returns
        -------
        P_saturation : float
            Saturation pressure in bar
        """

        if "CO2" not in oxide_wtPercents.index:
            raise ValueError("CO2 not found in sample")
        if oxide_wtPercents["CO2"] < 0:
            raise ValueError(f"CO2 lower than 0: {oxide_wtPercents['H2O']}")
        if oxide_wtPercents["CO2"] == 0.0:
            return 0.0

        composition = oxide_wtPercents.copy()

        # Upper limit of 15kbar
        upper_limit = 1.5e4
        # try:
        P_saturation = root_scalar(
            co2._root_saturation,
            args=(composition, kwargs),
            bracket=[1e-15, upper_limit],
        ).root
        # except:
        #     P_saturation = np.nan
        return P_saturation

    @staticmethod
    def calculate_solubility(
        oxide_wtPercents: MagmaSeries, P_bar: float, x_fluid: float = 0.0, **kwargs
    ) -> float:
        """
        Calculate melt |CO2| solubility according to equation 13.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composition in oxide wt. %.
        P_bar : float, array-like
            Pressure in bar
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
        if x_fluid == 1.0:
            return 0.0

        P_MPa = P_bar / 10
        # Partial pressure (ideal gas fugacity)
        fCO2 = (1 - x_fluid) * P_MPa

        model = kwargs.get("model", "pi_star")
        A = co2_parameters[model]["A"]
        B = co2_parameters[model]["B"]
        C = co2_parameters[model]["C"]

        mol_fractions = oxide_wtPercents.cations()
        pi_star = co2._pi_star(mol_fractions)

        return np.exp(A * np.log(fCO2) + B * pi_star + C) / 1e4

    @staticmethod
    def _root_saturation(P_bar, oxide_wtPercents, kwargs):
        """ """
        composition = oxide_wtPercents.copy()

        return composition["CO2"] - co2.calculate_solubility(
            oxide_wtPercents=composition, P_bar=P_bar, **kwargs
        )

    @staticmethod
    def _pi_star(mol_fractions):
        """
        equation 11
        """
        cations = {"Ca", "K", "Na", "Mg", "Fe", "Si", "Al"}
        if cations.difference(mol_fractions.index):
            raise KeyError(
                f"missing cations: {cations.difference(mol_fractions.index)}"
            )

        a = (
            mol_fractions["Ca"]
            + 0.8 * mol_fractions["K"]
            + 0.7 * mol_fractions["Na"]
            + 0.4 * mol_fractions["Mg"]
            + 0.4 * mol_fractions["Fe"]
        )
        b = mol_fractions["Si"] + mol_fractions["Al"]

        return a / b


class mixed(Solubility_model):
    @staticmethod
    @_check_argument("output", [None, "PXfl", "P", "Xfl"])
    def calculate_saturation(
        oxide_wtPercents: MagmaSeries, output: str = "P", **kwargs
    ) -> float | Tuple[float, float]:
        """
        Calculate volatile saturation pressure for systems with mixed |CO2|-|H2O| fluids.

        Parameters
        ----------
        oxide_wtPercents : MagmaSeries
            melt composision in oxide wt. %. Needs to have 'H2O' and 'CO2' columns.
        output : str
            Output format. 'P' for pressure only, 'Xfl' for |H2O| fluid fraction only and 'PXfl' for both.

        Returns
        -------
        saturation : float, (float, float)
            Depending on the value of ``output``: saturation pressure in bar, |H2O| fluid fraction or (saturation pressure, fluid fraction)
        """

        composition = oxide_wtPercents.copy()

        P_H2O_saturation = h2o.calculate_saturation(composition, x_fluid=1.0)
        P_CO2_saturation = co2.calculate_saturation(composition, x_fluid=0.0)

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
            args=(composition),
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
        oxide_wtPercents: MagmaSeries,
        P_bar: float,
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
        H2O = h2o.calculate_solubility(composition, P_bar, x_fluid)
        composition["H2O"] = H2O
        CO2 = co2.calculate_solubility(composition, P_bar, x_fluid)

        return_dict = {"both": (H2O, CO2), "H2O": H2O, "CO2": CO2}
        return return_dict[output]

    @staticmethod
    def _saturation_rootFunction(P_x_fluid, oxide_wtPercents):
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
                x_fluid=x_fluid,
                output="both",
            )
        )

        return abs(calculated_concentrations - sample_concentrations)
