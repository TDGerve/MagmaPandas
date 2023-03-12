from abc import abstractmethod
from enum import Enum, Flag, auto
from typing import Protocol

import pandas as pd


class Unit(Enum):
    MOL_FRACTIONS = "mol fraction"
    WT_PERCENT = "wt. %"
    PPM = "ppm"
    UNKNOWN = None


class Datatype(Enum):
    CATION = "cation"
    OXIDE = "oxide"
    UNKNOWN = None


class Magma(Protocol):
    _units: Unit
    _datatype: Datatype
    _weights: pd.Series

    @property
    def elements(self):
        ...

    @property
    def moles(self):
        ...

    @property
    def cations(self):
        ...

    def convert_ppm_wtPercent(self):
        ...

    def convert_moles_wtPercent(self):
        ...

    def normalise(self):
        ...

    def recalculate(self):
        ...
