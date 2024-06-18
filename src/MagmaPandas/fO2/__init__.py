"""
Module with models for calculating melt |fO2| at the QFM buffer.
"""

from MagmaPandas.configuration import configuration
from MagmaPandas.fO2.IW import fO2_IW
from MagmaPandas.fO2.QFM import fO2_QFM


def calculate_fO2(T_K, P_bar, **kwargs):

    fO2_buffer = kwargs.get("fO2_buffer", configuration.fO2buffer)
    dfO2 = kwargs.get("dfO2", configuration.dfO2)

    fO2_model = globals()[f"fO2_{fO2_buffer}"]

    return fO2_model(logshift=dfO2, T_K=T_K, P_bar=P_bar)
