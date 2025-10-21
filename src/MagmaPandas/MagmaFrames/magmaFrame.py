"""
===========
MagmaFrames
===========
Module with the generic MagmaFrame class.
"""

import re
from typing import Dict, List

import elementMass as e
import numpy as np
import pandas as pd
from typing_extensions import Self

from MagmaPandas.Elements import element_weights, oxide_compositions
from MagmaPandas.enums import Datatype, Unit
from MagmaPandas.parse_io.validate import _check_argument, _check_attribute


class _MagmaLocIndexer(pd.core.indexing._LocIndexer):
    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.obj.recalculate(inplace=True)


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
    _metadata = ["_weights", "_units", "_datatype", "_recalc"]

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

        self._recalc = True

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
    def loc(self):
        """Extended version of pandas._LocIndexer. Ensures that metadata are updated"""
        return _MagmaLocIndexer(self)

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

    def drop(self, *args, **kwargs) -> Self:
        """Extended version of pandas.DataFrame.drop. Ensures that metadata are updated"""
        inplace = kwargs.get("inplace", False)
        dropped = super().drop(*args, **kwargs)
        if inplace:
            self.recalculate(inplace=True)
            return
        return dropped.recalculate()

    def __setitem__(self, key, value):
        """Extended version of pandas.DataFrame.__setitem__. Ensures that metadata are updated"""
        super().__setitem__(key, value)
        if getattr(self, "_recalc", True):
            self.recalculate(inplace=True)

    def recalculate(self, inplace=False) -> Self:
        """
        Recalculate element masses and total weight and updates metadata.
        """
        df = self if inplace else self.copy()
        # avoid __setitem__ recursion when setting df.loc[:, "total"]
        df._recalc = False

        try:
            df._weights = element_weights.weights_as_series(self.columns)

            if df._total:
                totals = df.loc[:, df.elements].sum(axis=1)
                df.loc[:, "total"] = totals.astype(df["total"].dtype).values
        finally:
            df._recalc = True

        if not inplace:
            return df

    # @_check_attribute("_units", ["wt. %", "ppm"])
    def moles(self, normalise=True) -> Self:
        """
        Data converted to mol fraction.
        """

        if self._units == Unit.MOL_FRACTIONS:
            return self.copy()

        if self._units == Unit.WT_PERCENT:
            return self._convert_moles_wtPercent(normalise=normalise)
        elif self._units == Unit.PPM:
            return self.convert_ppm_wtPercent()._convert_moles_wtPercent(
                normalise=normalise
            )

        return self.copy()

    def wt_pc(self, normalise=True) -> Self:
        """
        Data converted to wt. %.
        """

        if self._units == Unit.WT_PERCENT:
            return self.copy()

        if self._units == Unit.MOL_FRACTIONS:
            return self._convert_moles_wtPercent(normalise=normalise)
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

    def cations(self, normalise=True, norm_to=1, mol_fractions=True) -> Self:
        """
        Data converted to cation mol fraction
        """
        # Calculate oxide moles
        if (self._datatype == Datatype.CATION) & (
            mol_fractions & (self._units == Unit.MOL_FRACTIONS)
        ):
            return self.copy()

        moles = self[self.elements].moles(normalise=False)

        # Calculate cation moles
        cations_per_oxide = oxide_compositions.cation_amount(moles.elements)
        cations = moles[moles.elements].mul(cations_per_oxide)
        # Rename columns to cations
        cations.columns = oxide_compositions.cation_names(moles.elements)

        cations._datatype = Datatype.CATION
        cations = cations.recalculate()

        if not mol_fractions:
            cations = cations[cations.elements].mul(cations.weights)
            cations._units = Unit.WT_PERCENT
            norm_to = 100

        if not normalise:
            cations["total"] = cations.sum(axis=1)
            return cations

        # Normalise to 1
        total = cations.sum(axis=1)
        cations = cations.div(total, axis=0) * norm_to
        cations["total"] = norm_to

        # Set the right datatype and elements
        cations._datatype = Datatype.CATION
        cations.recalculate(inplace=True)

        return cations

    def oxides(self, normalise=True, oxidation_state: Dict[str, int] = {}) -> Self:
        """
        Data converted to oxides
        """

        if (self._datatype == Datatype.OXIDE) & (not bool(oxidation_state)):
            return self.copy()

        units = self._units

        cations = self[self.elements].cations(
            normalise=False
        )  # convert to cation mol fractions
        cation_names = cations.elements
        cation_element_names = [
            re.sub(r"\d+", "", e) for e in cation_names
        ]  # strip numbers/charges from the names

        cation_names_new = [
            (
                i
                if oxidation_state.get(j, None) is None
                else f"{j}{int(oxidation_state[j])}"
            )
            for i, j in zip(cation_names, cation_element_names)
        ]  # new names include oxidation state for elements with non-default values

        oxide_names = e.get_oxide_names(cation_names_new)
        cations_per_oxide = e.cation_numbers(oxide_names)

        oxides = cations.rename(
            columns={cation: oxide for cation, oxide in zip(cation_names, oxide_names)}
        ).recalculate()  # rename to oxides

        oxides = oxides.div(cations_per_oxide, axis=1)  # recalculate to oxides
        oxides["total"] = oxides[oxides.elements].sum(axis=1)
        oxides._datatype = Datatype.OXIDE

        if units == Unit.MOL_FRACTIONS:
            if not normalise:
                return oxides
            return oxides.normalise()

        oxides_wt_pc = oxides.wt_pc(normalise=False)
        # oxides_wt_pc["total"] = oxides_wt_pc[oxides_wt_pc.elements].sum(axis=1)
        if not normalise:
            return oxides_wt_pc
        return oxides_wt_pc.normalise()

    @property
    def oxygen(self) -> pd.Series:
        """
        oxygen per 1 mole of cations
        """
        # Calculate oxide moles
        if self._datatype != Datatype.CATION:
            cations = self.cations()
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
    def _convert_moles_wtPercent(self, normalise=True) -> Self:
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

        if not normalise:
            converted["total"] = converted[converted.elements].sum(axis=1)
            converted._units = units
            return converted.recalculate()

        # Normalise
        total = converted[converted.elements].sum(axis=1)
        converted = converted.div(total, axis=0)
        converted["total"] = converted[converted.elements].sum(axis=1)
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
        O = float(O)

        cations = self.cations()
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
            norm = float(to)
        elif self._units == Unit.WT_PERCENT:
            norm = 100.0
        else:
            norm = 1.0

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
        random_sample[random_sample < 0] = 0.0

        df = self.copy()
        df[df.elements] = random_sample

        return df
