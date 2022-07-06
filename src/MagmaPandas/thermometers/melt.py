import pandas as pd
import numpy as np


class melt_thermometers:
    def putirka2008_14(self, inplace=False, **kwargs):

        """Liquid thermometer

        Equation 14 from Putirka (2008) calculates liquiqdus temperature for liquid compositions. 
        Requires equilibrium with olivine.

        Parameters
        ----------

        P_bar : int or list-like
            crystallisation pressures. Indices need to match melt if using pd.Series.


        Returns
        -------
        pd.Series
            liquidus temperatures in degrees Kelvin.
        """
        import MagmaPandas as mp

        composition = self.copy(deep=True)
        composition = composition.fillna(0.0)

        if isinstance(composition, mp.MagmaFrame):
            elements = composition.columns
        elif isinstance(composition, mp.MagmaSeries):
            elements = composition.index

        oxides = set(["MgO", "FeO", "Na2O", "K2O"])
        absentOxides = oxides.difference(elements)

        if "H2O" not in elements:
            H2O = 0.0
        else:
            H2O = composition["H2O"]
            # moles are calculated on an anhydrous basis
            try:
                composition = composition.drop("H2O")
            except KeyError:
                composition = composition.drop(columns=["H2O"])
        composition.recalculate()
        composition = composition.normalise()
  

        if len(absentOxides) > 0:
            raise KeyError(f"{absentOxides} not found in melt")

        # Calculate molar oxide fractions
        mol_fractions = composition.moles
        # Melt Mg#
        Mg_no = mol_fractions["MgO"] / (mol_fractions["MgO"] + mol_fractions["FeO"]) # SHOULD PROBABLY BE STRICTLY Fe2+?

        T_K = (
            754
            + 190.6 * Mg_no
            + 25.52 * composition["MgO"]
            + 9.585 * composition["FeO"]
            + 14.87 * (composition["Na2O"] + composition["K2O"])
            - 9.176 * H2O
        ) + 273.15

        return pd.Series(T_K, name="T_K").squeeze()

    def putirka2008_15(self, P_bar, inplace=False, **kwargs):

        """Liquid thermometer

        Equation 15 from Putirka (2008) calculates liquiqdus temperature for liquid compositions. 
        Requires equilibrium with olivine.

        Parameters
        ----------

        P_bar : int or list-like
            crystallisation pressures. Indices need to match melt if using pd.Series.


        Returns
        -------
        pd.Series
            liquidus temperatures in degrees Kelvin.
        """
        import MagmaPandas as mp

        composition = self.copy()
        composition = composition.fillna(0.0)

        if isinstance(self, mp.MagmaFrame):
            elements = self.columns
        elif isinstance(self, mp.MagmaSeries):
            elements = self.index

        oxides = set(["MgO", "FeO", "Na2O", "K2O"])
        absentOxides = oxides.difference(elements)

        if "H2O" not in elements:
            H2O = 0.0
        else:
            H2O = composition["H2O"]
            # moles are calculated on an anhydrous basis
            try:
                composition = composition.drop("H2O")
            except KeyError:
                composition = composition.drop(columns=["H2O"])
        composition.recalculate()
        composition = composition.normalise()

        P_GPa = P_bar / 1e4  

        if len(absentOxides) > 0:
            raise KeyError(f"{absentOxides} not found in melt")

        # Calculate molar oxide fractions
        mol_fractions = composition.moles
        # Melt Mg#
        Mg_no = mol_fractions["MgO"] / (mol_fractions["MgO"] + mol_fractions["FeO"]) # SHOULD PROBABLY BE STRICTLY Fe2+

        T_K = (
            815.3
            + 265.5 * Mg_no
            + 15.37 * composition["MgO"]
            + 8.61 * composition["FeO"]
            + 6.646 * (composition["Na2O"] + composition["K2O"])
            + 39.16 * P_GPa
            - 12.83 * H2O
        ) + 273.15

        return pd.Series(T_K, name="T_K").squeeze()

    def putirka2008_16(self, P_bar=None, inplace=False, **kwargs):

        """Liquid thermometer

        Equation 16 from Putirka (2008) calculates liquiqdus temperature for liquid compositions. 
        Requires equilibrium with olivine + plagioclase + clinopyroxene.

        Parameters
        ----------

        P_bar : int or list-like
            crystallisation pressures. Indices need to match melt if using pd.Series.


        Returns
        -------
        pd.Series
            liquidus temperatures in degrees Kelvin.
        """
        import MagmaPandas as mp

        if P_bar is None:
            P_bar = self["P_bar"]

        if isinstance(self, mp.MagmaFrame):
            elements = self.columns
        elif isinstance(self, mp.MagmaSeries):
            elements = self.index

        if isinstance(P_bar, pd.Series):
            if not self.index.equals(P_bar.index):
                raise RuntimeError("Melt and P_bar indices don't match")

        oxides = set(["SiO2", "Al2O3", "MgO"])
        absentOxides = oxides.difference(elements)

        if len(absentOxides) > 0:
            raise KeyError(f"{absentOxides} not found in melt")

        composition = self.copy()

        if "H2O" in elements:
            # moles are calculated on an anhydrous basis
            try:
                composition = composition.drop("H2O")
            except KeyError:
                composition = composition.drop(columns=["H2O"])
        composition.recalculate()
        composition = composition.normalise()

        # Convert pressure from bars to GPa
        P_GPa = P_bar / 1e4

        # Calculate molar oxide fractions
        mol_fractions = composition.moles
        mol_fractions = mol_fractions.fillna(0.0)

        part_1 = (
            -583
            + 3141 * mol_fractions["SiO2"]
            + 15779 * mol_fractions["Al2O3"]
            + 1338.6 * mol_fractions["MgO"]
        )
        part_2 = -31440 * mol_fractions["SiO2"] * mol_fractions["Al2O3"] + 77.67 * P_GPa

        T_K = part_1 + part_2 + 273.15

        if inplace:
            self["T_K"] = T_K
            if "P_bar" not in elements:
                self["P_bar"] = P_bar
        else:
            return pd.Series(T_K, name="T_K").squeeze()


