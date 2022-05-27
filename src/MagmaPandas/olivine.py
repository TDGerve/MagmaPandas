from typing import List
import pandas as pd
from .magma_baseclass import MagmaBase, read_file


def read_olivine(file: str, *args, index_col: List[str], keep_columns: List[str] = [], units='wt. %', type='oxide', **kwargs):
    """
    Docstring
    """

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return olivine(df, *args, keep_columns=keep_columns, calculate_total=True, units=units, datatype=type,**kwargs)


class olivine(MagmaBase):

    @property
    def forsterite(self):
        """
        Docstrings
        """
        cations = self.cations
        return pd.Series(cations['Mg'] * 100 / (cations['Fe'] + cations['Mg']), name='Fo#')

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=4)
