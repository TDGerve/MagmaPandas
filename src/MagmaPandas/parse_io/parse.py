import warnings as w
from typing import List, Union

import pandas as pd


def match_indeces(data: Union[pd.DataFrame, pd.Series], variables: List[pd.Series]):

    for var in variables:

        try:
            if not data.index.equals(var.index):
                raise RuntimeError(f"Data and {var.name} indices don't match")
        except AttributeError:
            pass


def convert_to_series(variables: List, index):

    new_vars = []
    for var in variables:
        if isinstance(var, (int, float)):
            new_var = pd.Series(var, index=index)
        else:
            new_var = var.copy()
        new_vars.append(new_var)

    return new_vars


def check_components(melt_mol_fractions, components):
    moles = melt_mol_fractions.copy()

    if isinstance(moles, pd.DataFrame):
        missing_oxides = set(components).difference(moles.columns)
    elif isinstance(moles, pd.Series):
        missing_oxides = set(components).difference(moles.index)

    if len(missing_oxides) == 0:
        return moles

    for oxide in missing_oxides:
        moles[oxide] = 0.0

    w.warn(
        f"{', '.join(str(i) for i in missing_oxides)} missing in composition and set to 0."
    )

    return moles.recalculate()
