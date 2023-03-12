from typing import List

import pandas as pd

import elements as e
from MagmaPandas import MagmaFrames

from .validate import _check_argument


@_check_argument("phase", [None, "Melt", "Olivine", "Clinopyroxene", "Plagioclase"])
def _read_file(
    filepath: str,
    *args,
    total_col: str,
    index_col: List[str] = None,
    keep_columns: List[str] = None,
    phase: str = None,
    units: str = None,
    datatype: str = None,
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
    if keep_columns is None:
        keep_columns = []

    df = pd.read_csv(filepath, *args, index_col=index_col, **kwargs)

    delete_columns = set()
    elements = []

    # Check which column names have valid names for elements or oxides
    for col in df.columns:
        try:
            _ = e.calculate_weight(col)
            elements.append(col)
        except:
            delete_columns.add(col)

    if total_col is not None:
        df.rename(columns={total_col: "total"})
        df["total"] = df[elements].sum(axis=1)

    if "total" in df.columns:
        keep_columns.append("total")

    # Drop all columns without chemical data, unless explicitely specified otherwise in 'keep_columns'
    df = df.drop(delete_columns.difference(set(keep_columns)), axis=1)
    df[elements] = df[elements].astype("float32")

    create_class_instance = getattr(MagmaFrames, phase)

    return create_class_instance(
        df, units=units, datatype=datatype, total_col=total_col
    )
