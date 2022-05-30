from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..parse.readers import _read_file


def read_plagioclase(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs
) -> "plagioclase":
    """
    Read plagioclase compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="plagioclase",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs
    )


class plagioclase(MagmaFrame):

    @property
    def anorthite(self):
        """
        Docstrings
        """
        cations = self.cations
        return pd.Series(cations['Ca'] * 100 / (cations['Ca'] + cations['Na']), name='An')

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=8)
