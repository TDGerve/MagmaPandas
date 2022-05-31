from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..parse.readers import _read_file
from ..geochemistry.fO2 import fO2_QFM
from ..geochemistry.Fe_redox import FeRedox_KC, FeRedox_Boris
from ..thermometers.melt import melt_thermometers
from ..parse.validate import _check_argument
from ..geochemistry.Kd import Kd_blundy_iterator, Kd_toplis_iterator


def read_melt(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> "melt":
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="melt",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs,
    )


class melt(melt_thermometers, MagmaFrame):
    def Fe3Fe2_QFM(
        self, T_K=None, P_bar=None, logshift=0, model="Borisov", inplace=False
    ):
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
        model   :   string
            'KressCarmichael' or 'Borisov'
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

        mol_fractions = self.moles

        model_dict = {
            "KressCarmichael": FeRedox_KC,
            "Borisov": FeRedox_Boris,
        }
        equation = model_dict[model]

        fO2_bar = fO2_QFM(logshift, T_K, P_bar)

        Fe3Fe2 = equation(mol_fractions, T_K, fO2_bar, P_bar)

        if inplace:
            self["Fe3Fe2"] = Fe3Fe2
            if "T_K" not in self.columns:
                self["T_K"] = T_K
            if "Pbar" not in self.columns:
                self["Pbar"] = P_bar
        else:
            return Fe3Fe2

    @_check_argument("total_Fe", ["FeO", "Fe2O3"])
    def FeO_Fe2O3_calc(self, Fe3Fe2, total_Fe, inplace=False):
        """
        melt_mol_fractions : pd.DataFrame
            melt composition in oxide mol fraction.
        Fe3Fe2 :
            melt Fe3+/Fe2+ ratio
        total_Fe    :
            column in melt_mol_fractions with total Fe
        """

        if len(set(("FeO", "Fe2O3")).intersection(self.columns)) == 2:
            raise RuntimeError("FeO and Fe2O3 columns already exist")

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
        melt_mol_fractions.recalculate()

        # Recalculate to wt. % (normalised)
        melt_mol_fractions = melt_mol_fractions.convert_moles_wtPercent

        if inplace:
            self["FeO"] = melt_mol_fractions["FeO"]
            self["Fe2O3"] = melt_mol_fractions["Fe2O3"]
            self.recalculate()

        else:
            return melt_mol_fractions

    @_check_argument("model", [None, "toplis", "blundy"])
    def Kd_olivine_FeMg(self, forsterite, T_K, Fe3Fe2, P_bar=None, model=None):
        """
        Calulate Fe-Mg exchange coefficients between olivine (ol) and melt (m) as:

        [Fe(ol) / Fe(m)] * [Mg(m) / Mg(ol)]
        """
        if model is None:
            model = "toplis"

        if (model == "toplis") & (P_bar is None):
            raise ValueError("P_bar argument missing")

        model_dict = {"toplis": Kd_toplis_iterator, "blundy": Kd_blundy_iterator}
        Kd_model = model_dict[model]

        mol_fractions = self.moles

        return Kd_model(
            melt_mol_fractions=mol_fractions,
            forsterite=forsterite,
            T_K=T_K,
            Fe3Fe2=Fe3Fe2,
            P_bar=P_bar,
        )
