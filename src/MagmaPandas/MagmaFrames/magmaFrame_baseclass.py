from typing import List
from ..parse.validate import _check_attribute
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

    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, **kwargs):
            return MagmaFrame(*args, **kwargs).__finalize__(self)

        return _c

    def __init__(
        self, df=None, *args, units: str = None, datatype: str = None, **kwargs
    ) -> None:

        super().__init__(df, *args, **kwargs)

        # A pandas series with the masses of all oxides and elements in the dataframe
        self._weights = pd.Series(name="weight", dtype="float32")
        # A list with the names of all columns that do not contain chemical data
        self._units = units
        self._datatype = datatype

        for col in self.columns:
            try:
                self._weights[col] = e.calculate_weight(col)
            except:
                self._no_data.append(col)

        # Recalculate total concentrattions
        if "total" in self.columns:
            self["total"] = self[self.elements].sum(axis=1)

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
        moles = self.copy()[self.elements]
        # Calculate moles
        moles = moles.div(moles.weights)
        # Normalise to 1
        total = moles.sum(axis=1)
        moles = moles.div(total, axis=0)
        moles["total"] = moles.sum(axis=1)
        # Set the right units
        moles._units = "mol fraction"

        return moles

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
