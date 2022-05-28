import pandas as pd
import elements as e
from typing import List
from .. import MagmaFrames as mf


def _check_units(func):
    def wrapper(*args, **kwargs):
        units = kwargs["units"]
        allowed_units = ["wt. %", "mol fraction"]
        if not units in allowed_units:
            raise ValueError(
                f'units: "{units}" not recognised, please choose from: {allowed_units}'
            )
        return func(*args, **kwargs)

    return wrapper


def _check_type(func):
    def wrapper(*args, **kwargs):
        type = kwargs["Type"]
        allowed_types = ["oxide", "cation"]
        if type not in allowed_types:
            raise ValueError(
                f'type: "{type}" not recognised, please choose from: {allowed_types}'
            )
        return func(*args, **kwargs)

    return wrapper


@_check_units
@_check_type
def _read_file(
    filepath: str,
    *args,
    phase: str = None,
    index_col: List[str],
    keep_columns: List[str] = [],
    total_col: str = None,
    units: str = None,
    Type: str = None,
    **kwargs,
):
    df = pd.read_csv(filepath, *args, index_col=index_col, **kwargs)

    delete_columns = set()
    keep_columns = set(keep_columns)
    elements = []

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
