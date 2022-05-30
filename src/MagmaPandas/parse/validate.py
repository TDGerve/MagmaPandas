import pandas as pd
import elements as e
from typing import List
from functools import wraps


def _check_argument(var_name: str, allowed_values: List[str]):
    def decorator(func):
        """
        Check if var_name has a valid value
        """
        @wraps(func)
        def wrapper(*args, **kwargs):
            var = kwargs.get(var_name, None)
            if var not in allowed_values:
                raise ValueError(
                    f"{var_name}: {var}, not recognised, please choose from: {allowed_values}"
                )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def _check_attribute(attr_name: str, allowed_values: List[str]):
    def decorator(func):
        """
        Check if attr_name has a valid value
        """

        def wrapper(self, *args, **kwargs):
            attr = getattr(self, attr_name)
            if attr not in allowed_values:
                raise ValueError(
                    f"Calculation is not valid with: {self.units}, please use: {*allowed_values,}"
                )
            return func(self, *args, **kwargs)

        return wrapper

    return decorator