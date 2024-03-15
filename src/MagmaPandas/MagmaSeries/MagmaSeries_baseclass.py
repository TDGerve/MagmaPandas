from typing import List

import elementMass as e
import numpy as np
import pandas as pd
from typing_extensions import Self

from MagmaPandas.configuration import configuration
from MagmaPandas.Elements import element_weights, oxide_compositions
from MagmaPandas.enums import Datatype, Unit
from MagmaPandas.parse_io.validate import _check_argument, _check_attribute
from MagmaPandas.thermometers import melt_thermometers


def _MagmaSeries_expanddim(data=None, *args, **kwargs):
    from MagmaPandas.MagmaFrames import MagmaFrame

    if isinstance(data, MagmaSeries):
        kwargs["units"] = data._units
        kwargs["datatype"] = data._datatype
        kwargs["weights"] = data._weights.copy(deep=True)

    df = MagmaFrame(data, *args, **kwargs)

    # df._units = data._units.copy(deep=True)
    # df._datatype = data._datatype.copy(deep=True)

    return df


# # pd.concat (pandas/core/reshape/concat.py) requires this for the
# # concatenation of series since pandas 1.1
# # (https://github.com/pandas-dev/pandas/commit/f9e4c8c84bcef987973f2624cc2932394c171c8c)
# _MagmaSeries_expanddim._get_axis_number = pd.DataFrame._get_axis_number


class MagmaSeries(pd.Series):
    """
    Generic MagmaSeries class for chemical data.

    Parameters
    ----------
    data : array-like, Iterable, dict, or scalar value
        geochemical data with elements or oxides as index
    units : None, str
        data units, either "mol fraction", "wt. %" or "ppm"
    datatype : None, str
        datatype either "cation" or "oxide"
    weights : None, pandas Series
        atomic weights of elements or oxides in the MagmaSeries
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
        weights: pd.Series = None,
        **kwargs,
    ) -> None:
        self._units: Unit = Unit(units)
        self._datatype: Datatype = Datatype(datatype)

        super().__init__(data, *args, **kwargs)

        if weights is not None:
            self._weights = weights.copy(deep=True)
        elif not hasattr(self, "_weights"):
            self._weights = element_weights.weights_as_series(self.index)

    @property
    def _constructor(self):
        """
        This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over (see: https://github.com/pandas-dev/pandas/issues/13208).
        We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time.
        """

        def _c(*args, **kwargs):
            if (weights := getattr(self, "_weights", None)) is not None:
                weights = weights.copy(deep=True)

            # The finalize call crashes the Python kernel when using a dictionary as data:
            # 'Cannot execute code, session has been disposed. Please try restarting the Kernel.'
            # WHY?
            return MagmaSeries(*args, weights=weights, **kwargs)  # .__finalize__(self)

        return _c

    @property
    def _constructor_expanddim(self):
        def _c(*args, **kwargs):
            if (weights := getattr(self, "_weights", None)) is not None:
                weights = weights.copy(deep=True)

            return _MagmaSeries_expanddim(
                *args, weights=weights, **kwargs
            ).__finalize__(self)

        # pd.concat (pandas/core/reshape/concat.py) requires this for the
        # concatenation of series since pandas 1.1
        # (https://github.com/pandas-dev/pandas/commit/f9e4c8c84bcef987973f2624cc2932394c171c8c)
        _c._get_axis_number = pd.DataFrame._get_axis_number

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
        Datatype and units.
        """
        return f"{self._datatype.value} {self._units.value}"

    @units.setter
    def units(self, value):
        print("units are read only")

    @property
    def weights(self) -> pd.Series:
        """
        Atomic weights of all elements in the MagmaSeries.
        """
        return self._weights

    @property
    def elements(self) -> List[str]:
        """
        Names of all elements in the MagmaSeries
        """
        return list(self._weights.index)

    @property
    @_check_attribute("_datatype", ["oxide"])
    def moles(self) -> Self:
        """
        Data converted to mol fractions.
        """
        if self._units != Unit.MOL_FRACTIONS:
            return self.convert_moles_wtPercent()
        else:
            return self

    @property
    @_check_attribute("_datatype", ["oxide"])
    def cations(self) -> Self:
        """
        Data converted to cation mol fractions
        """
        # Calculate oxide moles
        if self._units != Unit.MOL_FRACTIONS:
            moles = self.moles[self.elements]
        else:
            moles = self[self.elements].copy()

        cations = moles.mul(
            pd.Series(
                oxide_compositions.cation_amount(moles.elements), index=moles.elements
            )
        )
        # Rename index to cations
        cations.index = oxide_compositions.cation_names(cations.elements)

        # Normalise to 1
        total = cations.sum()
        cations = cations.div(total)
        cations["total"] = 1
        # Set the right datatype and elements
        cations._datatype = Datatype.CATION
        cations.recalculate(inplace=True)

        return cations

    def convert_moles_wtPercent(self) -> Self:
        """
        moles converted to wt. % and vice versa
        """

        converted = self.copy()[self.elements]
        if self._units == Unit.WT_PERCENT:
            converted = converted.div(converted.weights)
        elif self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(converted.weights)
        # Normalise
        total = converted.sum()
        converted = converted.div(total)
        converted["total"] = converted.sum()
        # Set the right units
        if self._units == Unit.WT_PERCENT:
            converted._units = Unit.MOL_FRACTIONS
        elif self._units == Unit.MOL_FRACTIONS:
            converted = converted.mul(100)
            converted._units = Unit.WT_PERCENT

        return converted

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

    def mineral_formula(self, O: int = None) -> Self:
        """
        Calculate mineral formulas by normalising to oxygen per formula unit

        Parameters
        ----------
        O : int
            Amount of oxygen to normalise to.

        Returns
        -------
        mineral formulas : MagmaSeries
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

    def recalculate(self, inplace=False) -> Self:
        """
        Recalculate element masses and total weight.
        """

        series = self if inplace else self.copy()

        series._weights = element_weights.weights_as_series(series.index)

        if series._total:
            series["total"] = series[series.elements].sum()

        if not inplace:
            return series

    def normalise(self, to=None) -> Self:
        """
        Normalise compositions.

        Parameters
        ----------
        to :    float, int
            normalisation value

        Returns
        -------
        normalised data : MagmaSeries
        """
        if to is not None:
            norm = to
        elif self._units == Unit.WT_PERCENT:
            norm = 100
        else:
            norm = 1

        self.recalculate(inplace=True)
        normalised = self.copy()[self.elements]

        total = normalised.sum()
        normalised = normalised.div(total)
        normalised = normalised.mul(norm)
        normalised["total"] = normalised.sum()

        return normalised

    def temperature(self, *args, **kwargs) -> float:
        """
        Calculate melt liquidus temperature.
        Thermometer models are selected in the global configuration.

        #TODO add link to the config docs.

        Parameters
        ----------
        P_bar   : float, pandas Series
            pressure in bar

        Returns
        -------
        temperatures : float
            Liquidus temperature in Kelvin
        """

        thermometer = melt_thermometers[configuration.melt_thermometer]

        return thermometer(self, *args, **kwargs)

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
        resampled data : MagmaSeries
            Randomly resampled compositions
        """

        random_sample = np.random.normal(self[self.elements], errors)
        random_sample[random_sample < 0] = 0

        sr = self.copy()
        sr[sr.elements] = random_sample

        return sr
