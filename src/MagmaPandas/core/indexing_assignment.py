import pandas as pd
from typing_extensions import Self

from MagmaPandas.core.magma_protocol import Magma


def _recalculate_check(data: Magma) -> bool:
    """checks if data needs to be recalculated"""
    return hasattr(data, "recalculate") and getattr(data, "_recalc", True)


class _MagmaLocIndexer(pd.core.indexing._LocIndexer):
    """Indexer modified to update all metadata"""

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if _recalculate_check(self.obj):
            self.obj.recalculate(inplace=True)

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if _recalculate_check(result):
            return result.recalculate()
        return result


class _MagmaILocIndexer(pd.core.indexing._iLocIndexer):
    """Integer-location indexer that triggers recalc on read and write."""

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if _recalculate_check(self.obj):
            self.obj.recalculate(inplace=True)

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if _recalculate_check(result):
            return result.recalculate()
        return result


class indexing_assignment_mixin:
    """
    Mixing for adding updated indexing and assignment methods to MagmaFrames and MagmaSeries
    """

    @property
    def loc(self) -> _MagmaLocIndexer:
        """Extended version of pandas._LocIndexer. Ensures that metadata are updated"""
        return _MagmaLocIndexer("loc", self)

    @property
    def iloc(self) -> _MagmaILocIndexer:
        return _MagmaILocIndexer("iloc", self)

    def drop(self, *args, **kwargs) -> Self:
        """Extended version of pandas.DataFrame.drop. Ensures that metadata are updated"""
        inplace = kwargs.get("inplace", False)
        dropped = super().drop(*args, **kwargs)
        if inplace and _recalculate_check(self):
            self.recalculate(inplace=True)
            return
        return dropped.recalculate()

    def pop(self, item):
        """Extended version of pandas.DataFrame.pop. Ensures that metadata are updated"""
        result = super().pop(item)
        if _recalculate_check(self):
            self.recalculate(inplace=True)
        return result

    def __setitem__(self, key, value):
        """Extended version of pandas.DataFrame.__setitem__. Ensures that metadata are updated"""
        super().__setitem__(key, value)
        if _recalculate_check(self):
            self.recalculate(inplace=True)

    def __getitem__(self, key):
        """Extended version of pandas.DataFrame.__getitem__. Ensures that metadata are updated"""
        result = super().__getitem__(key)
        if _recalculate_check(result):
            return result.recalculate()
        return result
