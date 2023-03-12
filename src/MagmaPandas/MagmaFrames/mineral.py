from abc import ABC, abstractmethod
from typing import List

import pandas as pd

from MagmaPandas.Magma_baseclass import Unit
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.parse_io.readers import _read_file


class Mineral(ABC):
    @property
    @abstractmethod
    def formula(self) -> MagmaFrame:
        pass


def read_mineral(
    filepath: str,
    phase: str,
    *args,
    index_col: List[str] = None,
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs
) -> Mineral:
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase=phase,
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        datatype="oxide",
        **kwargs
    )


class Olivine(Mineral, MagmaFrame):
    @property
    def forsterite(self):
        """
        Docstrings
        """
        if self._units == Unit.WT_PERCENT:
            moles = self.moles
        else:
            moles = self

        type = self._datatype.value

        Mg = {"oxide": "MgO", "cation": "Mg"}
        Fe = {"oxide": "FeO", "cation": "Fe"}
        self.recalculate(inplace=True)
        return pd.Series(
            moles[Mg[type]] / (moles[Fe[type]] + moles[Mg[type]]),
            name="Fo#",
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=4)


class Clinopyroxene(Mineral, MagmaFrame):
    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)


class Plagioclase(Mineral, MagmaFrame):
    @property
    def anorthite(self):
        """
        Docstrings
        """
        cations = self.cations
        return pd.Series(
            cations["Ca"] * 100 / (cations["Ca"] + cations["Na"]), name="An"
        )

    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=8)
