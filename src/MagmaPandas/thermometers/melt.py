import pandas as pd


class melt_thermometers:
    def temperature(self, P_bar=None, inplace=False):

        """Liquid thermometer

        Equation 16 of Putirka (2008) calculates liquiqdus temperature for liquid compositions. Requires equilibrium with olivine + plagioclase + clinopyroxene.

        Parameters
        ----------

        P_bar : int or list-like
            crystallisation pressures. Indices need to match melt if using pd.Series.


        Returns
        -------
        pd.Series
            liquidus temperatures in degrees Kelvin.
        """

        if P_bar is None:
            P_bar = self["P_bar"]

        if isinstance(P_bar, pd.Series):
            if not self.index.equals(P_bar.index):
                raise RuntimeError("Melt and P_bar indices don't match")

        oxides = set(["SiO2", "Al2O3", "MgO"])
        absentOxides = oxides.difference(self.columns)

        if len(absentOxides) > 0:
            raise KeyError(f"{absentOxides} not found in melt")

        # Convert pressure from bars to GPa
        P_GPa = P_bar / 1e4

        # Calculate molar oxide fractions
        mol_fractions = self.moles

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
            if "P_bar" not in self.columns:
                self["P_bar"] = P_bar
        else:
            return pd.Series(T_K, name="T_K").squeeze()
