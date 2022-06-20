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
        Type="oxide",
        **kwargs
    )


class Clinopyroxene(MagmaFrame):

    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, weights=self._weights, **kwargs):
            return Clinopyroxene(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
