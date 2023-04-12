import warnings as w
from multiprocessing import Pool

import numpy as np
import pandas as pd
from alive_progress import alive_bar

from MagmaPandas import volatile_solubility as vs
from MagmaPandas.configuration import configuration
from MagmaPandas.Fe_redox import FeRedox_QFM
from MagmaPandas.Kd.Ol_melt import calculate_FeMg_Kd
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.parse_io.validate import _check_argument
from MagmaPandas.thermometers import melt_thermometers


class Melt(MagmaFrame):
    # @property
    # def _constructor(self):
    #     """This is the key to letting Pandas know how to keep
    #     derivatives of `MagmaBase` the same type as yours.  It should
    #     be enough to return the name of the Class.  However, in
    #     some cases, `__finalize__` is not called and `new attributes` are
    #     not carried over.  We can fix that by constructing a callable
    #     that makes sure to call `__finalize__` every time."""

    #     def _c(*args, weights=None, **kwargs):
    #         if weights is None:
    #             weights = self._weights.copy(deep=True)
    #         return Melt(*args, weights=weights, **kwargs).__finalize__(self)

    #     return _c

    def temperature(self, *args, **kwargs):
        thermometer = melt_thermometers[configuration.melt_thermometer]

        return thermometer(self, *args, **kwargs)

    def tetrahedral_cations(self):
        """
        Tetrahedral cations
        Si, Ti, Al and P are assumed to be in tetrahedral coordination. Fe3+ is not taken into account. See Mysen (1983)\ [1]_ for additional information

        Returns
        -------
        pd.Series
            summed tertrahedral cations per 1 mole total cations

        References
        ----------
        .. [1] Mysen, B. O. (1983) The structure of silicate melts. Ann. Rev. Earth Planet. Sci. 11. 75-97
        """
        cations = self.cations

        tetrahedral_cations = {"Si", "Ti", "Al", "P"}

        available_elements = tetrahedral_cations.intersection(cations.elements)

        return cations[list(available_elements)].sum(axis=1)

    def NBO(self):
        """
        Non-bridging oxygen in the melt
        Formulation according to Mysen (1983)\ [1]_

        Returns
        -------
        pd.Series
            NBO per 1 mole cations

        References
        ----------
        .. [1] Mysen, B. O. (1983) The structure of silicate melts. Ann. Rev. Earth Planet. Sci. 11. 75-97
        """

        oxygen = self.oxygen

        return (2 * oxygen) - (4 * self.tetrahedral_cations())

    def NBO_T(self):
        """
        NBO/T
        The ratio of non-bridging oxygen and tetrahedral cations. Formulation according to Mysen (1983)\ [1]_

        Returns
        -------
        pd.Series
            NBO/T

        References
        ----------
        .. [1] Mysen, B. O. (1983) The structure of silicate melts. Ann. Rev. Earth Planet. Sci. 11. 75-97

        """

        return self.NBO() / self.tetrahedral_cations()

    def Fe3Fe2_QFM(self, T_K=None, P_bar=None, inplace=False):
        """
        Calculate Fe-redox equilibrium at QFM oxygen buffer for silicate liquids.
        Uses either equation 7 from Kress and Carmichael (1991) or equation 4 from Borisov et al. (2018).

        Parameters
        ----------

        T_K :   float, pd.Series-like
            temperature in Kelvin
        Pbar    :   float, pd.Series-like
            Pressure in bars
        logshift    :   int, pd.Series-like
            log units shift of QFM
        inplace :   bool
            return a new dataframe of add columns to the existing one
        Returns
        -------
        melt Fe3+/Fe2+ ratio

        """
        if T_K is None:
            T_K = self["T_K"]
        if P_bar is None:
            P_bar = self["P_bar"]

        for name in ["T_K", "P_bar"]:
            param = locals()[name]
            if isinstance(param, pd.Series):
                if not self.index.equals(param.index):
                    raise RuntimeError(f"Melt and {name} indices don't match")

        logshift = configuration.dQFM

        mol_fractions = self.moles

        Fe3Fe2 = FeRedox_QFM(
            mol_fractions=mol_fractions, T_K=T_K, P_bar=P_bar, logshift=logshift
        )

        if inplace:
            self["Fe3Fe2"] = Fe3Fe2
            if "T_K" not in self.columns:
                self["T_K"] = T_K
            if "Pbar" not in self.columns:
                self["Pbar"] = P_bar
        else:
            return Fe3Fe2

    @_check_argument("total_Fe", ["FeO", "Fe2O3"])
    def FeO_Fe2O3_calc(self, Fe3Fe2, total_Fe="FeO", inplace=False):
        """
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        Fe3Fe2 :
            melt Fe3+/Fe2+ ratio
        total_Fe    :
            column in melt_mol_fractions with total Fe
        """

        Fe2Fe_total = 1 / (1 + Fe3Fe2)
        melt_mol_fractions = self.moles

        if total_Fe == "FeO":
            Fe2 = melt_mol_fractions["FeO"] * Fe2Fe_total
            Fe3 = melt_mol_fractions["FeO"] * (1 - Fe2Fe_total) / 2
        if total_Fe == "Fe2O3":
            Fe2 = melt_mol_fractions["Fe2O3"] * Fe2Fe_total * 2
            Fe3 = melt_mol_fractions["Fe2O3"] * (1 - Fe2Fe_total)

        melt_mol_fractions["FeO"] = Fe2
        melt_mol_fractions["Fe2O3"] = Fe3
        melt_mol_fractions.recalculate(inplace=True)

        # Recalculate to wt. % (normalised)
        melt = melt_mol_fractions.convert_moles_wtPercent

        if inplace:
            self["FeO"] = melt["FeO"]
            self["Fe2O3"] = melt["Fe2O3"]
            self.recalculate(inplace=True)

        else:
            return melt

    def Kd_olivine_FeMg(self, forsterite, T_K, Fe3Fe2, **kwargs):
        """
        Calulate Fe-Mg exchange coefficients between olivine (ol) and melt (m) as:

        [Fe(ol) / Fe(m)] * [Mg(m) / Mg(ol)]
        """
        mol_fractions = self.moles

        return calculate_FeMg_Kd(
            melt_mol_fractions=mol_fractions,
            forsterite_initial=forsterite,
            T_K=T_K,
            Fe3Fe2=Fe3Fe2,
            **kwargs,
        )

    def volatile_saturation_pressure(self, T_K, inplace=False, **kwargs):
        """
        Calculate volatile (H2O and/or CO2) saturation pressures for given liquid compositions.

        Returns
        -------
        P_bar   : pd.Series
            Saturation pressures in bar
        """

        if isinstance(T_K, pd.Series):
            T_K = T_K[self.index]
        elif isinstance(T_K, (float, int)):
            T_K = pd.Series(T_K, index=self.index)

        P_bar = pd.Series(index=self.index, dtype=float)
        Xfl = pd.Series(index=self.index, dtype=float)
        names = self.index.values

        model = configuration.volatile_solubility
        species = configuration.volatile_species

        samples = [
            (name, temperature, (model, species), kwargs)
            for name, temperature in zip(names, T_K)
        ]
        total = self.shape[0]

        # Run calculations in a pool across multiple cpu cores
        # Progress bar from alive_bar
        with Pool() as pool, alive_bar(
            total,
            spinner=None,
            length=30,
            manual=True,
            theme="smooth",
            force_tty=True,
        ) as bar:
            results = pool.imap_unordered(self._saturation_multicore, samples)
            pool.close()

            finished = 0
            bar(0.0 / total)

            for name, result in results:
                finished += 1
                bar(finished / total)
                try:
                    pressure, fl = result
                    P_bar[name] = pressure
                    Xfl[name] = fl
                except TypeError:
                    P_bar[name] = result

        if inplace:
            self["P_bar"] = P_bar
            return
        else:
            output = P_bar
            if sum(Xfl.isna()) != Xfl.shape[0]:
                output = [P_bar, Xfl]
            return output

    def _saturation_multicore(self, sample):
        """
        Refactor of calculate_saturation for multiprocess calling
        """

        name, temperature, model, kwargs = sample
        solubility_model = vs.get_solubility_model(*model)
        try:
            result = vs.calculate_saturation(
                self.loc[name],
                T_K=temperature,
                solubility_model=solubility_model,
                **kwargs,
            )
        except Exception:
            result = np.nan
            w.warn(f"Saturation pressure not found for sample {name}")

        return name, result
