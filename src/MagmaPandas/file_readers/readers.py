import pandas as pd
import elements as e
from typing import List
from functools import wraps
from .. import MagmaFrames as mf


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


@_check_argument("units", ["wt. %", "mol fraction", "ppm"])
@_check_argument("Type", ["oxide", "cation"])
@_check_argument("phase", [None, "melt", "olivine", "clinopyroxene", "plagioclase"])
def _read_file(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str,
    keep_columns: List[str] = [],
    phase: str = None,    
    units: str = None,
    Type: str = None,
    **kwargs,
):
    """
    Read compositions from a .csv file into a pd.dataframe, clean the data and pass it to a magmaframe.

    Parameters
    ----------
    filepath    :   str
        path to a .csv file with sample compositions
    phase   :       str


    """
    if phase is None:
        phase = "MagmaFrame"

    df = pd.read_csv(filepath, *args, index_col=index_col, **kwargs)

    delete_columns = set()
    keep_columns = set(keep_columns)
    elements = []

    # Check which column names have valid names for elements or oxides
    for col in df.columns:
        try:
            _ = e.calculate_weight(col)
            elements.append(col)
        except:
            delete_columns.add(col)

    # Drop all columns without chemical data, unless explicitely specified otherwise in 'keep_columns'
    df = df.drop(delete_columns.difference(keep_columns), axis=1)
    df[elements] = df[elements].astype("float32")
    # Recalculate total concentrattions
    if total_col is not None:
        df.rename(columns={total_col: "total"})
        df["total"] = df[elements].sum(axis=1)

    create_class_instance = getattr(mf, phase)

    return create_class_instance(
        df,
        units=units,
        datatype=Type,
    )
