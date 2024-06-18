import warnings as w

import numpy as np
from scipy.constants import R

from MagmaPandas.parse_io.validate import _check_value


def fO2_NNO(logshift, T_K, P_bar):
    """
    from:
    Campbell et al. (2009) High pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers. EPSL

    Supplementary table S5

    calibrated between 2.34-65.9 GPa and 293-2467 K
    """
    P_GPa = P_bar / 1e4
    offset = 10**logshift

    part_1 = (
        8.699 + 1.642e-2 * P_GPa - 3e-4 * P_GPa**2 + 2.7e-6 * P_GPa**3 - 1e-8 * P_GPa**4
    )
    part_2 = (-24205 + 444.73 * P_GPa - 5.929e-1 * P_GPa**2 + 1.53e-3 * P_GPa) / T_K

    log10fO2 = part_1 + part_2

    return 10**log10fO2 * offset


@_check_value(var_name="T_K", allowed_range=[700, 1700], error=False)
def fO2_NNO_1bar(
    logshift: float | np.ndarray, T_K: float | np.ndarray
) -> float | np.ndarray:
    """
    Equation 6 from:

    O'Neill and Pownceby (1993) Thermodynamic data from redox reactions at high temperatures. I. An experimental and theoretical assessment of the electrochemical method using stabilized zirconia electrolytes, with revised values for the Fe-“FeO”, Co-CoO, Ni-NiO and Cu-Cu2O oxygen buffers, and new data for the W-WO2 buffer. Contributions to mineralogy and petrology. 114
    """

    try:
        logshift = float(logshift)
    except TypeError:
        pass

    offset = 10**logshift

    muO2 = -478967 + 248.514 * T_K - 9.7961 * np.log(T_K)  # J/mol
    fO2 = np.exp(muO2 / (R * T_K)) * offset

    return fO2
