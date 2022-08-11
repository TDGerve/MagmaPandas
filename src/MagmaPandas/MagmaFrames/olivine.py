from typing import List
import pandas as pd
from .magmaFrame_baseclass import MagmaFrame
from ..parse.readers import _read_file


def read_olivine(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs
) -> "Olivine":
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="Olivine",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs
    )


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
        if self._units == "wt. %":
            moles = self.moles
        else:
            moles = self
        Mg = {"oxide": "MgO", "cation": "Mg"}
        Fe = {"oxide": "FeO", "cation": "Fe"}
        self.recalculate(inplace=True)
        return pd.Series(
            moles[Mg[self._datatype]] / (moles[Fe[self._datatype]] + moles[Mg[self._datatype]]), name="Fo#"
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=4)
