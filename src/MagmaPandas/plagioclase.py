from typing import List
import pandas as pd
from .magma_baseclass import MagmaBase


def read_plagioclase(file: str, *args, index_col: List[str], keep_columns: List[str] = [], **kwargs):
    """
    Docstring
    """

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return plagioclase(df, *args, keep_columns=keep_columns, calculate_total=True, **kwargs)


class plagioclase(MagmaBase):

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
