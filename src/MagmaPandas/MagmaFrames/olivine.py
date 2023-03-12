from typing import List

import pandas as pd

from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame

from ..enums import Unit


class Olivine(MagmaFrame):

    # @property
    # def _constructor(self):
    #     """This is the key to letting Pandas know how to keep
    #     derivatives of `MagmaBase` the same type as yours.  It should
    #     be enough to return the name of the Class.  However, in
    #     some cases, `__finalize__` is not called and `new attributes` are
    #     not carried over.  We can fix that by constructing a callable
    #     that makes sure to call `__finalize__` every time."""

    #     def _c(*args, weights=None, **kwargs):
    #         if weights is None:
    #             weights = self._weights.copy(deep=True)
    #         return Olivine(*args, weights=weights, **kwargs).__finalize__(self)

    #     return _c

    @property
    def forsterite(self):
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
            name="Fo#",
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=4)
