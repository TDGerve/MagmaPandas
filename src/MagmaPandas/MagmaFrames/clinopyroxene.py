from typing import List
import pandas as pd
from .magmaFrame import MagmaFrame
from ..parse_io.readers import _read_file


def read_clinopyroxene(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs
) -> "Clinopyroxene":
    """
    Read cpx compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="Clinopyroxene",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs
    )


class Clinopyroxene(MagmaFrame):
    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
