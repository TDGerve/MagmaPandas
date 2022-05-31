from typing import List
from ..parse.readers import _read_file
from .melt import melt


def read_melt(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> "melt_inclusion":
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="melt_inclusion",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs,
    )


class melt_inclusion(melt):
    def Fe_loss_correction(self):
        """
        Docstrings

        """
        # Do some corrections
