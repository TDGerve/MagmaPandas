from typing import List

import pandas as pd

from ..parse_io.readers import _read_file
from .magmaFrame import MagmaFrame


class Clinopyroxene(MagmaFrame):
    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
