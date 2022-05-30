from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..parse.readers import _read_file


def read_clinopyroxene(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs
) -> "clinopyroxene":
    """
    Read cpx compositions in wt. % oxide from a .csv file

    """

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


class clinopyroxene(MagmaFrame):

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
