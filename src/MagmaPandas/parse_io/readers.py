from typing import List

import elementMass as e
import pandas as pd

from MagmaPandas import MagmaFrames
from MagmaPandas.MagmaFrames import (
    Clinopyroxene,
    MagmaFrame,
    Melt,
    Olivine,
    Plagioclase,
)
from MagmaPandas.parse_io.validate import _check_argument


@_check_argument("phase", [None, "Melt", "Olivine", "Clinopyroxene", "Plagioclase"])
def read_file(
    filepath: str,
    *args,
    total_col: str,
    index_col: List[str] = None,
    keep_columns: List[str] = None,
    phase: str = None,
    units: str = None,
    datatype: str = None,
    **kwargs,
) -> MagmaFrame:
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

    return create_class_instance(df, units=units, datatype=datatype)


def read_clinopyroxene(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> Clinopyroxene:
    """
    Read cpx compositions in wt. % oxide from a .csv file

    """

    return read_file(
        filepath=filepath,
        *args,
        phase="Clinopyroxene",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs,
    )


def read_melt(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    units="wt. %",
    datatype="oxide",
    **kwargs,
) -> Melt:
    """
    Read melt compositions in wt. % oxide from a .csv file

    """

    return read_file(
        filepath=filepath,
        *args,
        phase="Melt",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units=units,
        datatype=datatype,
        **kwargs,
    )


def read_olivine(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> Olivine:
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return read_file(
        filepath=filepath,
        *args,
        phase="Olivine",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs,
    )


def read_plagioclase(
    filepath: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> Plagioclase:
    """
    Read plagioclase compositions in wt. % oxide from a .csv file

    """

    return read_file(
        filepath=filepath,
        *args,
        phase="Plagioclase",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs,
    )
