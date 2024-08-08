"""
===========
MagmaFrames
===========
Module with the generic MagmaFrame class.
"""

import re
from typing import List

import elementMass as e
import numpy as np
import pandas as pd
from typing_extensions import Self

from MagmaPandas.Elements import element_weights, oxide_compositions
from MagmaPandas.enums import Datatype, Unit
from MagmaPandas.MagmaFrames.protocols import Magma
from MagmaPandas.parse_io.validate import _check_argument, _check_attribute


class MagmaFrame(pd.DataFrame):
    """
    Generic MagmaPandas DataFrame class for geochemical data.

    Parameters
    ----------
    data : ndarray (structured or homogeneous), Iterable, dict, or DataFrame
        geochemical data with elements or oxides in columns
    units : None, str
        data units, either "mol fraction", "wt. %" or "ppm"
    datatype : None, str
        datatype either "cation" or "oxide"
    weights : None, pandas Series
        atomic weights of elements or oxides in the MagmaFrame
    """

    # New attributes
    _metadata = ["_weights", "_units", "_datatype"]

    @_check_argument("units", [None, "mol fraction", "wt. %", "ppm"])
    @_check_argument("datatype", [None, "cation", "oxide"])
    def __init__(
        self,
        data=None,
        *args,
        units: None | str = None,
        datatype: None | str = None,
        weights: pd.Series = None,
        **kwargs,
    ) -> None:
        self._units: Unit = Unit(units)
        self._datatype: Datatype = Datatype(datatype)

        super().__init__(data, **kwargs)

        if not self._total:
            total_regex = re.search(
                "total", "".join(map(str, self.columns.to_list())), re.IGNORECASE
            )
            if total_regex:
                self.rename(columns={total_regex[0]: "total"}, inplace=True)

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

        def _c(*args, **kwargs):
            if (weights := getattr(self, "_weights", None)) is not None:
                weights = weights.copy(deep=True)

            current_class = type(self)

            return current_class(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def _constructor_sliced(self):
        from MagmaPandas.MagmaSeries import MagmaSeries

        def _c(*args, **kwargs):
            if (weights := getattr(self, "_weights", None)) is not None:
                weights = weights.copy(deep=True)

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
        Datatype and units.
        """
        return f"{self._datatype.value} {self._units.value}"

    @units.setter
    def units(self, value):
        print("units are read only")

    @property
    def weights(self) -> pd.Series:
        """
        Atomic weights of all elements in the MagmaFrame.
        """
        return self._weights.copy()

    @property
    def elements(self) -> List[str]:
        """
        Names of all elements in the MagmaFrame.
        """
        return list(self._weights.index).copy()

    # @_check_attribute("_units", ["wt. %", "ppm"])
    @property
    def moles(self) -> Self:
        """
        Data converted to mol fraction.
        """

        if self._units == Unit.MOL_FRACTIONS:
            return self.copy()

        if self._units == Unit.WT_PERCENT:
            return self._convert_moles_wtPercent()
        elif self._units == Unit.PPM:
            return self.convert_ppm_wtPercent()._convert_moles_wtPercent()

        return self.copy()

    @property
    def wt_pc(self) -> Self:
        """
        Data converted to wt. %.
        """

        if self._units == Unit.WT_PERCENT:
            return self.copy()

        if self._units == Unit.MOL_FRACTIONS:
            return self._convert_moles_wtPercent()
        elif self._units == Unit.PPM:
            return self.convert_ppm_wtPercent()

        return self.copy()

    @property
    def ppm(self) -> Self:
        """
        Data converted to ppm.
        """

        if self._units == Unit.PPM:
            return self.copy()

        if self._units == Unit.WT_PERCENT:
            return self.convert_ppm_wtPercent()
        elif self._units == Unit.MOL_FRACTIONS:
            return self._convert_moles_wtPercent().convert_ppm_wtPercent()

        return self.copy()

    @property
    def cations(self) -> Self:
        """
        Data converted to cation mol fraction
        """
        # Calculate oxide moles
        if self._datatype == Datatype.CATION:
            return self.copy()

        moles = self.moles[self.elements]

        # Calculate cation moles
        cations = moles.mul(oxide_compositions.cation_amount(moles.elements))
        # Rename columns to cations
        cations.columns = oxide_compositions.cation_names(moles.elements)
        # Normalise to 1
        total = cations.sum(axis=1)
        cations = cations.div(total, axis=0)
        cations["total"] = 1.0

        # Set the right datatype and elements
        cations._datatype = Datatype.CATION
        cations.recalculate(inplace=True)

        return cations

    @property
    def oxygen(self) -> pd.Series:
        """
        oxygen per 1 mole of cations
        """
        # Calculate oxide moles
        if self._datatype != Datatype.CATION:
            cations = self.cations
            cations = cations[cations.elements]
        else:
            cations = self[self.elements].copy()

        oxygen_per_mole = oxide_compositions.oxygen_amount(
            cations.elements, type="cation"
        )
        cation_per_mole = oxide_compositions.cation_amount(
            cations.elements, type="cation"
        )

        oxygen_per_cation = oxygen_per_mole / cation_per_mole

        oxygen = cations.mul(oxygen_per_cation).sum(axis=1)

        return oxygen

    @_check_attribute("_units", ["wt. %", "ppm"])
    def convert_ppm_wtPercent(self) -> Self:
        """
        ppm converted to wt. % and vice versa
        """
        convert_dict = {
            Unit.WT_PERCENT: [1e4, Unit.PPM],
            Unit.PPM: [1e-4, Unit.WT_PERCENT],
        }

        converted = self.mul(convert_dict[self._units][0])
        converted._units = convert_dict[self._units][1]

        return converted

    @_check_attribute("_units", ["wt. %", "mol fraction"])
    def _convert_moles_wtPercent(self) -> Self:
        """
        moles converted to wt. % and vice versa
        """

        converted = self[self.elements].copy()
        if self._units == Unit.WT_PERCENT:
            converted = converted.div(converted.weights)
            units = Unit.MOL_FRACTIONS
        elif self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(converted.weights)
            units = Unit.WT_PERCENT
        # Normalise
        total = converted.sum(axis=1)
        converted = converted.div(total, axis=0)
        converted["total"] = converted.sum(axis=1)
        # Set the right units

        if self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(100)

        converted._units = units

        return converted

    def mineral_formula(self, O: int = None) -> Self:
        """
        Calculate mineral formulas by normalising to oxygen per formula unit

        Parameters
        ----------
        O : int
            Amount of oxygen to normalise to.

        Returns
        -------
        mineral formulas : MagmaFrame
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

    def recalculate(self, inplace=False) -> Self:
        """
        Recalculate element masses and total weight.
        """
        df = self if inplace else self.copy()

        df._weights = element_weights.weights_as_series(self.columns)

        if df._total:
            totals = df.loc[:, df.elements].sum(axis=1)
            df.loc[:, "total"] = totals.values

        if not inplace:
            return df

    def normalise(self, to=None) -> Self:
        """
        Normalise compositions.

        Parameters
        ----------
        to :    float, int
            normalisation value

        Returns
        -------
        normalised data : MagmaFrame
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

    def random_sample(self, errors) -> Self:
        """
        Randomly resample compositions within errors.

        Sampling distribution is assumed normal with measured values as means and errors as standard deviations.

        Parameters
        ----------
        errors  : float, array-like
            standard deviation of the normal distributions. Use int for a fixed value for all elements or an array for specific values for all elements in :py:attr:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame.elements`

        Returns
        -------
        resampled data : MagmaFrame
            Randomly resampled compositions
        """

        random_sample = np.random.normal(self[self.elements], errors)
        random_sample[random_sample < 0] = 0

        df = self.copy()
        df[df.elements] = random_sample

        return df
