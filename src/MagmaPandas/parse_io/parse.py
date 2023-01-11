import pandas as pd
from typing import List, Union


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
