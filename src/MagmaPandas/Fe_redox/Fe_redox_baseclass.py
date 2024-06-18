from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
import pandas as pd
from MagmaPandas.Fe_redox.Fe3Fe2_errors import (
    error_params_1bar,
    error_params_high_pressure,
)
from MagmaPandas.model_errors import _error_func
from MagmaPandas.parse_io.validate import _check_value
from scipy import interpolate

# Fe3+/Fe2+ limits of the moving standard deviation in a 30 point window of the validation dataset (provided at ./data/Fe3Fe2_validation_data.csv).
validation_limits_1bar = (0.0351966873706004, 5.948890681577911)
validation_limits_high_pressure = (0.052631579, 2.160641174)


@_check_value(var_name="Fe3Fe2", allowed_range=validation_limits_1bar, error=False)
def _Fe3Fe2_error_func(a, b, c, d, Fe3Fe2):
    return _error_func(data=Fe3Fe2, a=a, b=b, c=c, d=d)


@_check_value(
    var_name="Fe3Fe2", allowed_range=validation_limits_high_pressure, error=False
)
def _Fe3Fe2_spline(Fe3Fe2, parameters):
    return interpolate.splev(Fe3Fe2, parameters)


class Fe3Fe2_model(ABC):
    @abstractmethod
    def calculate_Fe3Fe2(
        cls, melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate melt |Fe3Fe2| ratios

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        pass

    @classmethod
    def get_error(
        cls, Fe3Fe2, pressure: Optional[pd.Series] = None, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Returns one standard deviation error on |Fe3Fe2| ratios, calculated from a compiled validation dataset.

        Parameters
        ----------
        Fe3Fe2  : array-like
            melt |Fe3Fe2| ratios
        pressure : array-like, optional
            pressures of each element in ``|Fe3Fe2|``. If this term is not included, errors will be calculated strictly at 1 bar.

        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """

        name = cls.__name__
        error_1bar = _Fe3Fe2_error_func(*error_params_1bar[name], Fe3Fe2=Fe3Fe2)

        if pressure is None:
            return error_1bar

        error_high_pressure = _Fe3Fe2_spline(
            Fe3Fe2=Fe3Fe2, parameters=error_params_high_pressure[name]
        )

        errors = error_1bar.copy()

        if isinstance(pressure, (float, int)):
            pressure = np.array([pressure])
            errors = np.array([errors])
            error_high_pressure = np.array([error_high_pressure])
        hp = pressure > 1
        errors[hp] = error_high_pressure[hp]

        return errors

    @classmethod
    def get_offset(
        cls, Fe3Fe2, offset_parameters, pressure=None, *args, **kwargs
    ) -> float | np.ndarray:

        return cls.get_error(Fe3Fe2=Fe3Fe2, pressure=pressure) * offset_parameters

    @staticmethod
    def get_offset_parameters(n: int = 1) -> float | np.ndarray:
        return np.random.normal(loc=0, scale=1, size=n)
