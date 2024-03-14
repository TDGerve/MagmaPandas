import warnings as w
from multiprocessing import Pool

import numpy as np
import pandas as pd
from alive_progress import alive_bar
from typing_extensions import Self

from MagmaPandas import volatile_solubility as vs
from MagmaPandas.configuration import configuration
from MagmaPandas.Fe_redox import FeRedox_QFM
from MagmaPandas.Kd.Ol_melt import calculate_FeMg_Kd
from MagmaPandas.liquid_density import calculate_density
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.parse_io.validate import _check_argument
from MagmaPandas.thermometers import melt_thermometers


class Melt(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with melt specific methods.
    """

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

    def temperature(self, P_bar: float | pd.Series = None, **kwargs) -> pd.Series:
        """
        Calculate melt liquidus temperatures.
        Model choice is set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

        Parameters
        ----------
        P_bar   : float, pandas Series
            pressure in bar

        Returns
        -------
        temperatures : pd.Series
            Liquidus temperatures in Kelvin
        """
        thermometer = melt_thermometers[configuration.melt_thermometer]

        return thermometer(self, P_bar=P_bar, **kwargs)

    def density(
        self,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        fO2_logshift: None | int = None,
    ) -> pd.Series:
        """
        Calculate silicate melts densities with the Iacovino and Till (2019)\ [3]_ model

        Parameters
        ----------
        T_K : float, pandas Series
            temperatures in Kelvin
        P_bar : float, pandas Series
            pressures in bar
        fO2_logshift : None, int
            |fO2| buffer shift in log units of QFM. If set to None, the value set in the global configuration is used.

        Returns
        -------
        densities : pd.Series
            densities in kg/m\ :sup:`3`
        """

        self._match_index(arg_names=("T_K", "P_bar"), kwargs=locals())

        if fO2_logshift is None:
            fO2_logshift = configuration.dQFM

        Fe3Fe2 = self.Fe3Fe2_QFM(T_K=T_K, P_bar=P_bar, fO2_logshift=fO2_logshift)
        # calculate melt composition with Fe2O3 and FeO
        melt = self.FeO_Fe2O3_calc(Fe3Fe2)

        return calculate_density(melt, T_K=T_K, P_bar=P_bar)

    def tetrahedral_cations(self) -> pd.Series:
        """
        Calculate tetrahedral cations based on Mysen (1983)\ [4]_

        Si, Ti, Al and P are assumed to be in tetrahedral coordination and Fe\ :sup:`3+` is not taken into account

        Returns
        -------
        tetrahedral cations : pd.Series
            summed tertrahedral cations per 1 mole cations

        """
        cations = self.cations

        tetrahedral_cations = {"Si", "Ti", "Al", "P"}

        available_elements = tetrahedral_cations.intersection(cations.elements)

        return cations[list(available_elements)].sum(axis=1)

    def NBO(self):
        """
        Non-bridging oxygen in the melt
        Formulation according to Mysen (1983)\ [4]_

        Returns
        -------
        pd.Series
            NBO per 1 mole cations

        """

        oxygen = self.oxygen

        return (2 * oxygen) - (4 * self.tetrahedral_cations())

    def NBO_T(self):
        """
        NBO/T
        The ratio of non-bridging oxygen and tetrahedral cations according to Mysen (1983)\ [4]_

        Returns
        -------
        pd.Series
            NBO/T
        """

        return self.NBO() / self.tetrahedral_cations()

    def Fe3Fe2_QFM(
        self,
        T_K: float | pd.Series,
        P_bar: float | pd.Series,
        fO2_logshift: None | int = None,
        inplace=False,
        **kwargs,
    ) -> pd.Series:
        """
        Calculate melt |Fe3Fe2| ratios at the QFM |fO2| buffer.
        Model choice is set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

        Parameters
        ----------

        T_K :   float, pd.Series-like
            temperatures in Kelvin
        Pbar    :   float, pd.Series-like
            Pressure in bars
        logshift    :   int, pd.Series-like
            |fO2| buffer shift in log units of QFM.
        inplace :   bool

        Returns
        -------
        melt Fe3+/Fe2+ ratios : pandas Series

        """
        # if T_K is None:
        #     T_K = self["T_K"]
        # if P_bar is None:
        #     P_bar = self["P_bar"]

        self._match_index(arg_names=["T_K", "P_bar"], kwargs=locals())

        # for name in ["T_K", "P_bar"]:
        #     param = locals()[name]
        #     if isinstance(param, pd.Series):
        #         if not self.index.equals(param.index):
        #             raise RuntimeError(f"Melt and {name} indices don't match")

        if fO2_logshift is None:
            fO2_logshift = configuration.dQFM

        mol_fractions = self.moles

        Fe3Fe2 = FeRedox_QFM(
            mol_fractions=mol_fractions,
            T_K=T_K,
            P_bar=P_bar,
            logshift=fO2_logshift,
            **kwargs,
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
    def FeO_Fe2O3_calc(
        self, Fe3Fe2: float | pd.Series, total_Fe: str = "FeO", inplace: bool = False
    ) -> Self:
        """
        Calculate melt FeO and |Fe2O3| based on total Fe.

        Parameters
        ----------
        Fe3Fe2 : pandas Series
            melt |Fe3Fe2| ratios
        total_Fe    : str
            columname in Melt frame with total Fe
        inplace : bool

        Returns
        -------
        Melt    : Self
            melt compositions inclusding FeO and |Fe2O3|
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
        melt = melt_mol_fractions.convert_moles_wtPercent()

        if inplace:
            self["FeO"] = melt["FeO"]
            self["Fe2O3"] = melt["Fe2O3"]
            self.recalculate(inplace=True)

        else:
            return melt

    def Kd_olivine_FeMg(
        self,
        forsterite: float | pd.Series,
        T_K: float | pd.Series,
        Fe3Fe2: float | pd.Series,
        **kwargs,
    ):
        """
        Calulate equilibrium Fe-Mg partitioning coefficients between olivine and melt as:

        (Fe\ :sup:`2+` / Mg)\ :sub:`ol` / (Fe\ :sup:`2+` / Mg)\ :sub:`melt`

        Model choice is set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

        Parameters
        ----------
        forsterite : pandas Series
            initial olivine forsterite contents
        T_K : pandas Series
            temperatures in Kelvin
        Fe3Fe2 : pandas Series
            melt |Fe3Fe2| ratios

        Returns
        -------
        Kds :   pandas Series
            Fe-Mg partitioning coefficients
        """
        mol_fractions = self.moles

        return calculate_FeMg_Kd(
            Melt_mol_fractions=mol_fractions,
            forsterite_initial=forsterite,
            T_K=T_K,
            Fe3Fe2=Fe3Fe2,
            **kwargs,
        )

    def volatile_saturation_pressure(
        self, T_K: float | pd.Series, inplace: bool = False, **kwargs
    ):
        """
        Calculate melt volatile (|CO2| and/or |H2O|) saturation pressures.

        Model choice is set in the global :py:class:`~MagmaPandas.configuration.configuration` class.

        Parameters
        ----------
        T_K : float, pandas Series
            temperatures in Kelvin
        inplace : bool

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
        except Exception as e:
            result = np.nan
            w.warn(f"Saturation pressure not found for sample {name}: {repr(e)}")

        return name, result

    def _match_index(self, arg_names: list[str], kwargs):

        for name in arg_names:
            param = kwargs[name]
            if isinstance(param, pd.Series):
                if not self.index.equals(param.index):
                    raise RuntimeError(f"Melt and {name} indices don't match")
