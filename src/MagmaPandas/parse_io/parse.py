import itertools as it
import warnings as w
from typing import List, Union

import numpy as np
import pandas as pd


def match_indeces(data: Union[pd.DataFrame, pd.Series], variables: List[pd.Series]):

    for var in variables:

        try:
            if not data.index.equals(var.index):
                raise RuntimeError(f"Data and {var.name} indices don't match")
        except AttributeError:
            pass


def make_iterable(*args):
    return [np.array([a]).flatten() for a in args]


def convert_to_series(variables: List, index):

    new_vars = []
    for var in variables:
        if isinstance(var, (int, float)):
            new_var = pd.Series(var, index=index)
        else:
            new_var = var.copy()
        new_vars.append(new_var)

    return new_vars


def check_components(composition, components):
    comp = composition.copy()

    if isinstance(comp, pd.DataFrame):
        missing_oxides = set(components).difference(comp.columns)
    elif isinstance(comp, pd.Series):
        missing_oxides = set(components).difference(comp.index)

    if len(missing_oxides) == 0:
        return comp.fillna(0.0)

    for oxide in missing_oxides:
        comp[oxide] = 0.0

    w.warn(
        f"{', '.join(str(i) for i in missing_oxides)} missing in composition and set to 0."
    )

    return comp.fillna(0.0).recalculate()


def repeat_vars(var1, var2):

    var1_is_int, var2_is_int = (
        True if isinstance(var, (int, float)) else False for var in (var1, var2)
    )
    if bool(var1_is_int) & bool(var2_is_int):
        return make_iterable(var1, var2)
    # If only one variable, P or T, is a single value
    if not (bool(var1_is_int) ^ bool(var2_is_int)):
        return var1, var2

    var1, var2 = (np.array([var]).flatten() for var in (var1, var2))
    # Cycle the short variable
    var1 = np.repeat(var1, len(var2)) if var1_is_int else var1
    var2 = np.repeat(var2, len(var1)) if var2_is_int else var2

    return var1, var2
