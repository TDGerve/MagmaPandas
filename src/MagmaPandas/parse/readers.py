from .validate import _check_argument
from typing import List
import elements as e
import pandas as pd
from .. import MagmaFrames as mf

@_check_argument("phase", [None, "melt", "melt_inclusion", "olivine", "clinopyroxene", "plagioclase"])
def _read_file(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str,
    keep_columns: List[str] = None,
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
    if keep_columns is None:
        keep_columns = []

    df = pd.read_csv(filepath, *args, index_col=index_col, **kwargs)

    delete_columns = set()
    # keep_columns = set().add(keep_columns)
    elements = []

    # Check which column names have valid names for elements or oxides
    for col in df.columns:
        try:
            _ = e.calculate_weight(col)
            elements.append(col)
        except:
            delete_columns.add(col)

    # Drop all columns without chemical data, unless explicitely specified otherwise in 'keep_columns'
    df = df.drop(delete_columns.difference(set(keep_columns)), axis=1)
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
