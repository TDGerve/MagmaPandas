from importlib import resources

import pandas as pd


def get_periodic_table(masses_only: bool = True):
    """
    Docstring
    """

    if not masses_only:
        with resources.open_text("elements.data", "PeriodicTable.csv") as df:
            return pd.read_csv(df, index_col=["Symbol"])

    with resources.open_text("elements.data", "PeriodicTable.csv") as df:
        return pd.read_csv(
            df, usecols=["Symbol", "AtomicMass"], index_col=["Symbol"]
        ).squeeze("columns")
