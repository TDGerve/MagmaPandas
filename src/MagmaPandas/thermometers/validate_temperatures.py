import warnings as w

import numpy as np


def _check_temperature(T_K):
    """
    Find negative or NaN temperatures
    """
    try:
        if (neg := any(T_K < 0)) or (nan := any(np.isnan(T_K))):
            str_arr = np.array(["negative", "NaN"])[np.array([neg, nan])]
            w.warn(f"{', '.join(str_arr)} temperatures found!")
    except TypeError:
        if (neg := (T_K < 0)) or (nan := (T_K != T_K)):
            str_arr = np.array(["negative", "NaN"])[np.array([neg, nan])]
