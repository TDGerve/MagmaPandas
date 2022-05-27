from typing import List
import pandas as pd
from .magma_baseclass import MagmaBase


def read_clinopyroxene(file: str, *args, index_col: List[str], keep_columns: List[str] = [], **kwargs):
    """
    Docstring
    """

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return clinopyroxene(df, *args, keep_columns=keep_columns, calculate_total=True, **kwargs)


class clinopyroxene(MagmaBase):

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
