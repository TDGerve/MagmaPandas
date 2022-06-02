from numpy import isin
import pandas as pd
import numpy as np
from ..MagmaFrames import olivine
import elements as e
from typing import Union



def calculate_olivine(forsterite):
    """
    Calculate olivine wt. % composition from forsterite contents
    """
    if isinstance(forsterite, pd.Series):
        idx = forsterite.index
    else:
        try:
            idx = np.arange(0, len(forsterite))
        except:
            idx = [0]

    MgO = forsterite * 2
    FeO = (1 - forsterite) * 2
    SiO2 = 1
    total = MgO + FeO + SiO2

    cations = olivine({'MgO': MgO, 'FeO': FeO, 'SiO2': SiO2, 'total': total}, index=idx, units='mol fraction', datatype='oxide')
    cations = cations.convert_moles_wtPercent

    return cations

