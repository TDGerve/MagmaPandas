import inspect
import warnings as w
from functools import wraps
from typing import List

import numpy as np
import pandas as pd


def _check_argument(var_name: str, allowed_values: List[str]):
    def decorator(func):
        """
        Check if var_name has a valid value
        """
        # Get values for default arguments
        sig = inspect.signature(func)
        defaults = {
            n: val
            for n, val in zip(
                sig.parameters.keys(), [val.default for val in sig.parameters.values()]
            )
        }

        @wraps(func)  # comment
        def wrapper(*args, **kwargs):

            for key in kwargs.keys():
                # Remove user kwargs from the default dict
                defaults.pop(key, None)

            var = {**kwargs, **defaults}.get(var_name, None)
            if var not in allowed_values:
                raise ValueError(
                    f"{var_name}: {var}, not recognised, please choose from: {allowed_values}"
                )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def _check_setter(allowed_values):
    def decorator(func):
        """
        check if set value is valid
        """

        @wraps(func)  # comment
        def wrapper(*args):
            var = args[1]
            if isinstance(var, (list, tuple)):
                # tuples or list are used when setting fixed values for Kd or Fe3Fe2
                var = var[0]
            if isinstance(var, (float, int)):
                min, max = allowed_values
                if not (min < var < max):
                    raise ValueError(
                        f"value: {var}, outside allowed range: {*allowed_values,}"
                    )
            elif isinstance(var, str):
                if var not in allowed_values:
                    raise ValueError(
                        f"'{var}' is not recognised, please choose from: {*allowed_values,}"
                    )
            return func(*args)

        return wrapper

    return decorator


def _check_attribute(attr_name: str, allowed_values: List[str]):
    def decorator(func):
        """
        Check if attr_name has a valid value
        """

        @wraps(func)  # comment
        def wrapper(self, *args, **kwargs):
            attr = getattr(self, attr_name)
            if attr.value not in allowed_values:
                raise ValueError(
                    f"Calculation is not valid with: {attr.value}, please use: {*allowed_values,}"
                )
            return func(self, *args, **kwargs)

        # wrapper.__doc__ = func.__doc__
        return wrapper

    return decorator


def _check_value(var_name: str, allowed_range: List[float], error=True):
    def decorator(func):
        """
        Check if var_name has a valid value
        """

        @wraps(func)  # comment
        def wrapper(*args, **kwargs):
            var = kwargs.get(var_name, None)
            xmin, xmax = allowed_range
            try:
                # ints and floats
                within_range = xmin < var < xmax
            except ValueError:
                # array-likes
                within_range = np.any((var > xmin) & (var < xmax))
            if not within_range:
                if error:
                    raise ValueError(
                        f"(some) values of {var_name} are outside allowed range: {xmin:.3f} - {xmax:.3f}"
                    )
                else:
                    w.warn(
                        f"(some) values of {var_name} are outside allowed range: {xmin:.3f} - {xmax:.3f})"
                    )
            # if not (min < var < max):
            #     raise ValueError(
            #         f"{var_name}: {var}, outside allowed range: {min} - {max}"
            #     )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def _match_index(x, arg_names: list[str], kwargs):

    for name in arg_names:
        param = kwargs[name]
        if not isinstance(param, (pd.Series, pd.DataFrame)):
            continue
        if not x.index.equals(param.index):
            raise RuntimeError(f"index doesn't match with {name}")
