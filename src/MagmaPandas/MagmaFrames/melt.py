from typing import List
import pandas as pd
from alive_progress import alive_bar
import warnings as w
import numpy as np

from .magmaFrame_baseclass import MagmaFrame

from ..configuration import configuration

from ..geochemistry.Fe_redox import FeRedox_QFM
from ..geochemistry.Kd_ol_melt import Kd_FeMg_vectorised
from ..geochemistry.volatiles import calculate_saturation

from ..thermometers.melt import melt_thermometers

from ..parse.validate import _check_argument
from ..parse.readers import _read_file


def read_melt(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    units="wt. %",
    datatype="oxide",
    **kwargs,
) -> "Melt":
    """
    Read melt compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="Melt",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units=units,
        datatype=datatype,
        **kwargs,
    )


class Melt(MagmaFrame):
    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, weights=None, **kwargs):
            if weights is None:
                weights = self._weights.copy(deep=True)
            return Melt(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    def temperature(self, *args, **kwargs):

        thermometer = getattr(melt_thermometers, configuration.melt_thermometer)

        return thermometer(self, *args, **kwargs)

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
        Kd_model = getattr(Kd_FeMg_vectorised, configuration.Kd_model)

        mol_fractions = self.moles

        return Kd_model(
            melt_mol_fractions=mol_fractions,
            forsterite=forsterite,
            T_K=T_K,
            Fe3Fe2=Fe3Fe2,
            **kwargs,
        )

    def volatile_saturation_pressure(self, T_K, inplace=False):
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

        total = self.shape[0]

        with alive_bar(total, spinner=None, force_tty=True) as bar:
            for (name, composition), temperature in zip(self.iterrows(), T_K):
                bar.text = f"-> Processing sample '{name}'..."
                try:
                    P_bar[name] = calculate_saturation(composition, T_K=temperature)
                except Exception:
                    P_bar[name] = np.nan
                    w.warn(f"Saturation pressure not found for sample {name}")
                bar()

        if inplace:
            self["P_bar"] = P_bar
            return
        else:
            return P_bar
