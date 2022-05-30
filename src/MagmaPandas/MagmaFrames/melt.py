from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..parse.readers import _read_file
from .geochemistry import fO2, Fe_redox
from .thermometers.melt import melt_thermometers


def read_melt(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = [],
    **kwargs
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
        **kwargs
    )


class melt(melt_thermometers, MagmaFrame):
    def FeRedox_QFM(
        self, T_K=None, P_bar=None, logshift=0, model="borisov", inplace=False, **kwargs
    ):
        """
        Calculate Fe-redox equilibrium at QFM oxygen buffer for silicate liquids.
        Uses either equation 7 from Kress and Carmichael (1991) or equation 4 from Borisov et al. (2018).

        Parameters
        ----------

        T_K         float, pd.Series-like
            temperature in Kelvin
        Pbar        float, pd.Series-like
            Pressure in bars
        logshift    int, pd.Series-like
            log units shift of QFM
        model       string
            'KressCarmichael' or 'Borisov'

        Returns
        -------
        liquid Fe3+/Fe2+ ratio

        """
        if T_K is None:
            T_K = self["T_K"]
        if P_bar is None:
            P_bar = self["P_bar"]        
        
        for param in ['T_K', 'P_bar']:
            if isinstance(param, pd.Series):
                if not self.index.equals(P_bar.index):
                    raise RuntimeError(f'Melt and {param} indices don\'t match')

        mol_fractions = self.moles

        model_dict = {
            "KressCarmichael": Fe_redox.FeRedox_KC,
            "Borisov": Fe_redox.FeRedox_Boris,
        }
        equation = model_dict[model]

        fO2_bar = fO2.fO2_QFM(logshift, T_K, P_bar)

        Fe3Fe2 = equation(mol_fractions, T_K, fO2_bar, P_bar)

        if inplace:
            self["Fe3Fe2"] = Fe3Fe2
            if "T_K" not in self.columns:
                self["T_K"] = T_K
            if "Pbar" not in self.columns:
                self["Pbar"] = P_bar
        else:
            return pd.DataFrame(
                {"Fe3Fe2": Fe3Fe2, "T_K": T_K, "P_bar": P_bar, "QFM": logshift}
            )
