import warnings as w

import numpy as np


def _check_temperature(T_K):
    """
    Find negative or NaN temperatures
    """
    try:
        negative = any(T_K < 0)
        nan = any(np.isnan(T_K))
        if negative or nan:
            str_arr = np.array(["negative", "NaN"])[np.array([negative, nan])]
            w.warn(f"{', '.join(str_arr)} temperatures found!")
    except TypeError:
        negative = T_K < 0
        nan = T_K != T_K
        if negative or nan:
            str_arr = np.array(["negative", "NaN"])[np.array([negative, nan])]
