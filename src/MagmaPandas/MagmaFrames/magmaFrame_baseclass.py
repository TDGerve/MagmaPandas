from typing import List
from ..parse.validate import _check_attribute, _check_argument
from ..parse.readers import _read_file
import pandas as pd
import elements as e


def read_file(*args, **kwargs):

    return _read_file(*args, **kwargs)


class MagmaFrame(pd.DataFrame):
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
            self._weights = kwargs.pop("weights")

        super().__init__(data, *args, **kwargs)

        if not hasattr(self, "_weights"):
            self._weights = pd.Series(name="weight", dtype=float)
            for col in self.columns:
                try:
                    # Calculate element/oxide weight
                    self._weights[col] = e.calculate_weight(col)
                except:
                    pass



    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, weights=self._weights, **kwargs):
            return MagmaFrame(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    @property
    def _constructor_sliced(self):
        """
        Docstrings
        """

        from MagmaPandas.MagmaSeries import MagmaSeries

        def _c(*args, weights=self._weights, **kwargs):
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
        if "total" in self.columns:
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
    def weights(self) -> pd.Series:
        """
        Molar mass of the elements in the dataframe
        """
        return self._weights

    @property
    def elements(self) -> List:
        """
        Names of the elements in the dataframe
        """
        return list(self._weights.index)

    @property
    @_check_attribute("_datatype", ["oxide"])
    @_check_attribute("_units", ["wt. %"])
    def moles(self):
        """
        Calculate molar fractions from oxide concentrations
        """
        return self.convert_moles_wtPercent

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
        # Rename columns to cations
        cations.columns = e.cation_names(cations.elements)
        # Normalise to 1
        total = cations.sum(axis=1)
        cations = cations.div(total, axis=0)
        cations["total"] = cations.sum(axis=1)
        # Set the right datatype and elements
        cations._datatype = "cation"
        cations.recalculate()

        return cations

    @property
    def convert_moles_wtPercent(self):
        """
        Convert moles to wt. % or vice versa
        """

        converted = self.loc[:, self.elements]
        if self._units == "wt. %":
            converted = converted.div(converted.weights)
        elif self._units == "mol fraction":
            converted = converted.mul(converted.weights)
        # Normalise
        total = converted.sum(axis=1)
        converted = converted.div(total, axis=0)
        converted["total"] = converted.sum(axis=1)
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
        oxygen_total = cations.mul(oxygen_numbers).sum(axis=1)
        oxygen_factor = O / oxygen_total
        cations = cations.mul(oxygen_factor, axis=0)
        cations["O"] = O

        return cations

    def recalculate(self):
        """
        Recalculate element masses and total weight.
        """
        weights = pd.Series(name="weight", dtype="float32")

        for col in self.columns:
            try:
                weights[col] = e.calculate_weight(col)
            except:
                pass

        self._weights = weights
        if self._total:
            self["total"] = self[self.elements].sum(axis=1)

    def normalise(self):
        """
        Normalise composition.
        """
        self.recalculate()
        normalised = self.copy()[self.elements]
        total = normalised.sum(axis=1)

        normalised = normalised.div(total, axis=0)
        if self._units == "wt. %":
            normalised = normalised.mul(100)
        normalised["total"] = normalised.sum(axis=1)

        return normalised
