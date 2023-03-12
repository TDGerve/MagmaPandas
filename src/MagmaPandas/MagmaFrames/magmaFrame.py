"""
===========
MagmaFrames
===========
Module with the generic MagmaFrame class.
"""

from typing import List

import numpy as np
import pandas as pd

import elements as e
from MagmaPandas.Elements import element_weights, oxide_compositions

from ..Magma_baseclass import Datatype, Unit
from ..parse_io.validate import _check_argument, _check_attribute


class MagmaFrame(pd.DataFrame):
    """
    Generic MagmaPandas DataFrame class for processing geochemical data.
    """

    # New attributes
    _metadata = ["_weights", "_units", "_datatype"]

    @_check_argument("units", [None, "mol fraction", "wt. %", "ppm"])
    @_check_argument("datatype", [None, "cation", "oxide"])
    def __init__(
        self,
        data=None,
        *args,
        units: str = None,
        datatype: str = None,
        total_col: str = None,
        weights: pd.Series = None,
        **kwargs,
    ) -> None:

        self._units: Unit = Unit(units)
        self._datatype: Datatype = Datatype(datatype)

        super().__init__(data, **kwargs)

        if weights is not None:
            self._weights = weights.copy()
        elif not hasattr(self, "_weights"):

            self._weights = element_weights.weights_as_series(self.columns)

    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, weights=None, **kwargs):
            if weights is None:
                weights = getattr(self, "_weights", None).copy(deep=True)

            current_class = type(self)

            return current_class(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def _constructor_sliced(self):
        """
        Docstrings
        """

        from MagmaPandas.MagmaSeries import MagmaSeries

        def _c(*args, weights=None, **kwargs):
            if weights is None:
                weights = getattr(self, "_weights", None).copy(deep=True)

            return MagmaSeries(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def _no_data(self) -> List:
        """
        Names of all columns without chemical data
        """
        no_data = list(self.columns.difference(self.elements))
        if "total" in no_data:
            no_data.remove("total")
        return no_data

    @property
    def _total(self) -> bool:
        """
        Dataframe contains column with totals
        """

        return "total" in self.columns

    @property
    def units(self) -> str:
        """
        Data units
        """
        return f"{self._datatype.value} {self._units.value}"

    @units.setter
    def units(self, value):
        print("units are read only")

    @property
    def weights(self) -> pd.Series:
        """
        Molar mass of the elements in the dataframe
        """
        return self._weights.copy()

    @property
    def elements(self) -> List:
        """
        Names of the elements in the dataframe
        """
        return list(self._weights.index).copy()

    @property
    @_check_attribute("_units", ["wt. %", "ppm"])
    def moles(self):
        """
        Calculate molar fractions from oxide concentrations
        """
        if self._units != Unit.MOL_FRACTIONS:
            return self.convert_moles_wtPercent
        else:
            return self.copy()

    @property
    def cations(self):
        """
        Calculate cation fractions from oxide concentrations
        """
        # Calculate oxide moles
        if self._units != Unit.MOL_FRACTIONS:
            moles = self.moles[self.elements]
        else:
            moles = self[self.elements].copy()
        # Calculate cation moles
        cations = moles.mul(oxide_compositions.amount(moles.elements))
        # Rename columns to cations
        cations.columns = oxide_compositions.names(moles.elements)
        # Normalise to 1
        total = cations.sum(axis=1)
        cations = cations.div(total, axis=0)
        cations["total"] = 1.0
        # Set the right datatype and elements
        cations._datatype = Datatype.CATION
        cations.recalculate(inplace=True)

        return cations

    def convert_ppm_wtPercent(self):
        """
        Convert ppm to wt. % and vice versa
        """
        convert_dict = {
            Unit.WT_PERCENT: [1e4, Unit.PPM],
            Unit.PPM: [1e-4, Unit.WT_PERCENT],
        }

        converted = self.mul(convert_dict[self._units][0])
        converted._units = convert_dict[self._units][1]

        return converted

    @property
    def convert_moles_wtPercent(self):
        """
        Convert moles to wt. % or vice versa
        """

        converted = self[self.elements].copy()
        if self._units == Unit.WT_PERCENT:
            converted = converted.div(converted.weights)
        elif self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(converted.weights)
        # Normalise
        total = converted.sum(axis=1)
        converted = converted.div(total, axis=0)
        converted["total"] = converted.sum(axis=1)
        # Set the right units
        if self._units == Unit.WT_PERCENT:
            converted._units = Unit.MOL_FRACTIONS
        elif self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(100)
            converted._units = Unit.WT_PERCENT

        return converted

    def mineral_formula(self, O: int = None):
        """
        Calculate mineral formulas by normalising to oxygen
        """
        # Calculate cation fractions
        cations = self.cations
        cations = cations[cations.elements]
        # Calculate oxygens per cation
        oxygen_numbers = e.oxygen_numbers(self.elements) / e.cation_numbers(
            self.elements
        )
        oxygen_numbers.index = cations.elements
        # Normalise to oxygen
        oxygen_total = cations.mul(oxygen_numbers).sum(axis=1)
        oxygen_factor = O / oxygen_total
        cations = cations.mul(oxygen_factor, axis=0)
        cations["O"] = O

        return cations

    def recalculate(self, inplace=False):
        """
        Recalculate element masses and total weight.
        """
        df = self if inplace else self.copy()

        df._weights = element_weights.weights_as_series(self.columns)

        if df._total:
            totals = df.loc[:, df.elements].sum(axis=1)
            df["total"] = totals.values

        if not inplace:
            return df

    def normalise(self, to=None):
        """
        Normalise composition.
        """
        if to is not None:
            norm = to
        elif self._units == Unit.WT_PERCENT:
            norm = 100
        else:
            norm = 1

        # self = self.recalculate()
        normalised = self[self.elements].copy()

        total = normalised.sum(axis=1)

        normalised = normalised.div(total, axis=0)
        normalised = normalised.mul(norm, axis=0)
        normalised.loc[:, "total"] = normalised.sum(axis=1)

        return normalised
