import pandas as pd
from typing_extensions import Self


class _MagmaLocIndexer(pd.core.indexing._LocIndexer):
    """Indexer modified to update all metadata"""

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if hasattr(self.obj, "recalculate"):
            self.obj.recalculate(inplace=True)

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if hasattr(result, "recalculate"):
            return result.recalculate()
        return result


class _MagmaILocIndexer(pd.core.indexing._iLocIndexer):
    """Integer-location indexer that triggers recalc on read and write."""

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if hasattr(self.obj, "recalculate"):
            self.obj.recalculate(inplace=True)

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if hasattr(result, "recalculate"):
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
        if inplace:
            self.recalculate(inplace=True)
            return
        return dropped.recalculate()

    def pop(self, item):
        """Extended version of pandas.DataFrame.pop. Ensures that metadata are updated"""
        result = super().pop(item)
        if hasattr(self, "recalculate") and getattr(self, "_recalc", True):
            self.recalculate(inplace=True)
        return result

    def __setitem__(self, key, value):
        """Extended version of pandas.DataFrame.__setitem__. Ensures that metadata are updated"""
        super().__setitem__(key, value)
        if hasattr(self, "recalculate") and getattr(self, "_recalc", True):
            self.recalculate(inplace=True)

    def __getitem__(self, key):
        """Extended version of pandas.DataFrame.__getitem__. Ensures that metadata are updated"""
        result = super().__getitem__(key)
        if hasattr(result, "recalculate") and getattr(result, "_recalc", True):
            return result.recalculate()
        return result
