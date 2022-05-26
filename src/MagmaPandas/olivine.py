from typing import List
import pandas as pd
from .magma_baseclass import MagmaBase
import elements as e

def read_olivine(file: str, index_col: List[str], keep_columns: List[str] = [], **kwargs):
    """
    Docstring
    """

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return olivine(df, keep_columns=keep_columns, calculate_total=True, **kwargs)


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
        # Get oxygens per cation
        oxygen_numbers = e.oxygen_numbers(self.elements) / e.cation_numbers(self.elements)
        oxygen_numbers.index = e.cation_names(self.elements)
        # Calculate cation fractions
        cations = self.cations[oxygen_numbers.index]
        # Normalise to 4 oxygen
        oxygen_total = cations.mul(oxygen_numbers).sum(axis=1)
        oxygen_factor = 4. / oxygen_total
        cations = cations.mul(oxygen_factor, axis=0)
        cations['O'] = 4.
        return cations
