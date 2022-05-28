from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..file_readers.readers import _read_file


def read_olivine(filepath: str, *args, index_col: List[str], total_col=None, keep_columns: List[str]=[], **kwargs):

    return _read_file(
        filepath=filepath,
        *args,
        phase="olivine",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs
    )


class olivine(MagmaFrame):
    @property
    def forsterite(self):
        """
        Docstrings
        """
        cations = self.cations
        return pd.Series(
            cations["Mg"] * 100 / (cations["Fe"] + cations["Mg"]), name="Fo#"
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=4)
