import numpy as np
import pandas as pd


def Kd_isopleths(
    kd_min: float = 0.1, kd_max: float = 0.4, stepsize: float = 0.1, **kwargs
):
    Kds = np.round(np.arange(kd_min, kd_max + 0.01, stepsize), 2)
    mg_no_melt_min, mg_no_melt_max, mg_no_melt_stepsize = kwargs.pop(
        "mg_no_melt", (0.2, 1, 0.01)
    )
    mg_no_melt = np.arange(mg_no_melt_min, mg_no_melt_max, mg_no_melt_stepsize)
    df = pd.DataFrame(index=mg_no_melt)

    Fe2Mg_melt = (1 - mg_no_melt) / mg_no_melt
    for kd in Kds:
        df[round(kd, 2)] = np.array(1 / (1 + (kd * Fe2Mg_melt))) * 100

    return df
