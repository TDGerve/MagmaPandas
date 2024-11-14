from typing import Protocol

import pandas as pd

from MagmaPandas.enums import Datatype, Unit


class Magma(Protocol):
    _units: Unit
    _datatype: Datatype
    _weights: pd.Series

    @property
    def _no_data(self): ...

    @property
    def _total(self): ...

    @property
    def units(self): ...

    @property
    def weights(self): ...

    @property
    def elements(self): ...

    def moles(self): ...

    @property
    def wt_pc(self): ...

    def cations(self): ...

    @property
    def oxygen(self): ...

    def convert_ppm_wtPercent(self): ...

    def _convert_moles_wtPercent(self): ...

    def recalculate(self): ...

    def normalise(self): ...

    def random_sample(self): ...
