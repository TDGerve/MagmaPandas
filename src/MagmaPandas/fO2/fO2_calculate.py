"""
Module with models for calculating melt |fO2| at various buffers.
"""

import numpy as np

from MagmaPandas.configuration import configuration
from MagmaPandas.fO2 import IW, QFM


def calculate_fO2(T_K, P_bar, **kwargs):

    fO2_buffer = kwargs.get("fO2_buffer", configuration.fO2buffer)
    dfO2 = kwargs.get("dfO2", configuration.dfO2)

    fO2_model = globals()[fO2_buffer]

    fO2 = fO2_model.calculate_fO2(logshift=dfO2, T_K=T_K, P_bar=P_bar)

    try:
        fO2 = fO2.astype(np.float32)
    except AttributeError:
        fO2 = np.float32(fO2)

    return fO2


def change_fO2_buffer(to: str, T_K: float | int, P_bar: float | int):

    current_fO2 = calculate_fO2(T_K=T_K, P_bar=P_bar)

    new_fO2_model = globals()[f"fO2_{to}"]

    new_fO2 = new_fO2_model(logshift=0, T_K=T_K, P_bar=P_bar)

    new_dfO2 = np.log10(current_fO2 / new_fO2)
    new_fO2 = np.float32(new_fO2 * 10**new_dfO2)

    if not np.isclose(new_fO2, current_fO2, rtol=0, atol=1e-5):
        raise ValueError("new fO2 does not match old fO2")

    configuration.fO2buffer = to
    configuration.dfO2 = round(new_dfO2, 2)

    print(configuration)
