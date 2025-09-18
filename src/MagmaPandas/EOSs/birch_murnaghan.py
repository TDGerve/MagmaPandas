from typing import Dict

import numpy as np


def birch_murnaghan_4th_order(
    V: float | np.ndarray, params: Dict = None
) -> float | np.ndarray:
    """
    Formulation from Katsura & Tange (2019)\ [26]_ (eq. 27). Needs parameters for:

    - pre-compression volume (V_0)
    - isothermal bulk modulus at standard temperature (K_0)
    - pressure derivative of the isothermal bulk modulus (Kprime_0)
    - second pressure derivative of the bulk modulus (Kprime_prime_0)


    Parameters
    ----------
    V : float, array-like
        volume
    params : dictionary, pandas Series
        fitted parameters 'V_0', 'K_0', 'Kprime_0', 'Kprime_prime_0'

    Returns
    -------
    pressure in same units as K_0
    """

    x = params["V_0"] / V

    part_1 = (3 / 2) * params["K_0"] * (np.power(x, 7 / 3) - np.power(x, 5 / 3))
    part_2 = 1 + (3 / 4) * (params["Kprime_0"] - 4) * (np.power(x, 2 / 3) - 1)
    part_3 = (1 / 24) * (
        9 * np.power(params["Kprime_0"], 2)
        - 63 * params["Kprime_0"]
        + 9 * params["K_0"] * params["Kprime_prime_0"]
        + 143
    )
    part_4 = np.power((np.power(x, 2 / 3) - 1), 2)

    return part_1 * (part_2 + part_3 * part_4)


def birch_murnaghan_4th_order_stacey(
    V: float | np.ndarray, params: Dict = None
) -> float | np.ndarray:
    """
    Formulation from eq. 18 Stacey et al. (1981, Geophysical Surveys 4, pp. 189-232), after the matlab script from :cite:t:`Hirschmann2022`. Needs parameters for:

    - pre-compression volume (V_0)
    - isothermal bulk modulus at standard temperature (K_0)
    - pressure derivative of the isothermal bulk modulus (Kprime_0)
    - second pressure derivative of the bulk modulus (Kprime_prime_0)


    Parameters
    ----------
    V : float, array-like
        volume
    params : dictionary, pandas Series
        fitted parameters 'V_0', 'K_0', 'Kprime_0', 'Kprime_prime_0'

    Returns
    -------
    pressure in same units as K_0
    """

    P0 = 1 / 1e4  #

    a1 = 3 * params["V_0"] * P0
    a2 = 3 * params["V_0"] * (3 * params["K_0"] - 5 * P0) / 2
    a3 = (
        params["V_0"]
        * (9 * params["K_0"] * params["K_0"] - 36 * params["K_0"] + 35 * P0)
        / 2
    )
    a4 = (
        3
        * params["V_0"]
        * (
            9 * (params["K_0"] ** 2) * params["Kprime_prime_0"]
            + 9 * params["K_0"] * (params["K_0"] ** 2)
            - 63 * params["K_0"] * params["K_0"]
            + 143 * params["K_0"]
            - 105 * P0
        )
        / 8.0
    )
    f = 1 / 2 * ((params["V_0"] / V) ** (2 / 3) - 1)
    dfdV = -(params["V_0"] ** (2 / 3)) / 3 / (V ** (5 / 3))
    P = (a1 + 2 * a2 * f + 3 * a3 * (f**2) + 4 * a4 * (f**3)) * (-dfdV)

    return P
