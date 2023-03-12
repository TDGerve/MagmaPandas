from typing import Protocol

import pandas as pd

from MagmaPandas.enums import Datatype, Unit


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
