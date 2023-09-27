import pandas as pd

from ..enums import Unit
from .magmaFrame import MagmaFrame


class Clinopyroxene(MagmaFrame):
    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)

    @property
    def mg_no(self):
        """
        Docstrings
        """
        if self._units == Unit.WT_PERCENT:
            moles = self.moles
        else:
            moles = self

        type = self._datatype.value

        Mg = {"oxide": "MgO", "cation": "Mg"}
        Fe = {"oxide": "FeO", "cation": "Fe"}
        self.recalculate(inplace=True)
        return pd.Series(
            moles[Mg[type]] / (moles[Fe[type]] + moles[Mg[type]]),
            name="Mg_no",
        )
