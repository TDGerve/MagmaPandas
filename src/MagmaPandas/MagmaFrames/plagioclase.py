from typing import List

import pandas as pd

from ..parse_io.readers import _read_file
from .magmaFrame import MagmaFrame


class Plagioclase(MagmaFrame):

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
    #         return Plagioclase(*args, weights=weights, **kwargs).__finalize__(self)

    #     return _c

    @property
    def anorthite(self):
        """
        Docstrings
        """
        cations = self.cations
        return pd.Series(
            cations["Ca"] * 100 / (cations["Ca"] + cations["Na"]), name="An"
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=8)
