from typing import List
from ..parse.validate import _check_attribute, _check_argument
import pandas as pd
import elements as e
from ..configuration import configuration
from ..thermometers.melt import melt_thermometers


def _MagmaSeries_expanddim(data=None, *args, **kwargs):
    from MagmaPandas.MagmaFrames import MagmaFrame

    df = MagmaFrame(*args, **kwargs)

    if isinstance(data, MagmaSeries):
        df._units = data._units.copy(deep=True)
        df._datatype = data._datatype.copy(deep=True)

    return df


# pd.concat (pandas/core/reshape/concat.py) requires this for the
# concatenation of series since pandas 1.1
# (https://github.com/pandas-dev/pandas/commit/f9e4c8c84bcef987973f2624cc2932394c171c8c)
_MagmaSeries_expanddim._get_axis_number = pd.DataFrame._get_axis_number


class MagmaSeries(pd.Series):
    """
    Docstrings
    """

    # New attributes
    _metadata = ["_weights", "_units", "_datatype"]

    @_check_argument("units", [None, "mol fraction", "wt. %", "ppm"])
    @_check_argument("datatype", [None, "cation", "oxide"])
    def __init__(
        self, data=None, *args, units: str = None, datatype: str = None, **kwargs
    ) -> None:

        self._units = units
        self._datatype = datatype
        if "weights" in kwargs.keys():
            self._weights = kwargs.pop("weights").copy(deep=True)

        super().__init__(data, *args, **kwargs)

        if not hasattr(self, "_weights"):
            self._weights = pd.Series(name="weight", dtype=float)
            for idx in self.index:
                try:
                    # Calculate element/oxide weight
                    self._weights[idx] = e.calculate_weight(idx)
                except:
                    pass
        # self.recalculate()


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
                weights = self._weights.copy(deep=True)
            return MagmaSeries(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def _constructor_expanddim(self):
        def _c(*args, weights=None, **kwargs):
            if weights is None:
                weights = self._weights.copy(deep=True)
            return _MagmaSeries_expanddim(
                *args, weights=weights, **kwargs
            ).__finalize__(self)

        return _c

    @property
    def _no_data(self) -> List:
        """
        Names of all index without chemical data
        """
        no_data = list(self.index.difference(self.elements))
        if "total" in no_data:
            no_data.remove("total")
        return no_data

    @property
    def _total(self) -> bool:
        """
        Dataframe contains totals in index
        """
        if "total" in self.index:
            return True
        else:
            return False

    @property
    def units(self) -> str:
        """
        Data units
        """
        return f"{self._datatype} {self._units}"

    @units.setter
    def units(self, value):
        print("units are read only")

    @property
    def weights(self):
        return self._weights

    @property
    def elements(self) -> List:
        """
        Names of the elements in the series
        """
        return list(self._weights.index)

    @property
    @_check_attribute("_datatype", ["oxide"])
    def moles(self):
        """
        Calculate molar fractions from oxide concentrations
        """
        if self._units != "mol fraction":
            return self.convert_moles_wtPercent
        else:
            return self

    @property
    @_check_attribute("_datatype", ["oxide"])
    def cations(self):
        """
        Calculate cation fractions from oxide concentrations
        """
        # Calculate oxide moles
        moles = self.moles[self.elements]
        # Calculate cation moles
        cation_numbers = e.cation_numbers(moles.elements)
        cations = moles.mul(cation_numbers)
        # Rename index to cations
        cations.index = e.cation_names(cations.elements)
        # Normalise to 1
        total = cations.sum()
        cations = cations.div(total)
        cations["total"] = cations.sum()
        # Set the right datatype and elements
        cations._datatype = "cation"
        cations.recalculate()

        return cations

    @property
    def convert_moles_wtPercent(self):
        """
        Convert moles to wt. % or vice versa
        """

        converted = self.copy()[self.elements]
        if self._units == "wt. %":
            converted = converted.div(converted.weights)
        elif self._units == "mol fraction":
            converted = converted.mul(converted.weights)
        # Normalise
        total = converted.sum()
        converted = converted.div(total)
        converted["total"] = converted.sum()
        # Set the right units
        if self._units == "wt. %":
            converted._units = "mol fraction"
        elif self._units == "mol fraction":
            converted = converted.mul(100)
            converted._units = "wt. %"

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
        oxygen_total = cations.mul(oxygen_numbers).sum()
        oxygen_factor = O / oxygen_total
        cations = cations.mul(oxygen_factor)
        cations["O"] = O

        return cations

    def recalculate(self):
        """
        Recalculate element masses and total.
        """
        missing_elements = self.index.difference(self._weights.index)
        extra_elements = self._weights.index.difference(self.index)

        if all(i.size == 0 for i in[missing_elements, extra_elements]):
            return

        if extra_elements.size > 0:
            self._weights = self._weights.drop(extra_elements)        

        if missing_elements.size > 0:
            new_weights = pd.Series(name="weight", dtype="float32")
            for element in missing_elements:
                try:
                    new_weights[element] = e.calculate_weight(element)
                except:
                    pass
            self._weights = pd.concat([self._weights, new_weights])        
        
        if self._total:
            self["total"] = self[self.elements].sum()

    def normalise(self):
        """
        Normalise composition.
        """
        self.recalculate()
        normalised = self.copy()[self.elements]
        total = normalised.sum()

        normalised = normalised.div(total)
        if self._units == "wt. %":
            normalised = normalised.mul(100)
        normalised["total"] = normalised.sum()

        return normalised

    def melt_temperature(self, *args, **kwargs):
        """
        calculate liquidus temperature for melts
        """

        thermometer = getattr(melt_thermometers, configuration().melt_thermometer)

        return thermometer(self, *args, **kwargs)
