import numpy as np
import scipy.optimize as opt

from MagmaPandas.EOSs.tools import thermal_expansivity


def Vinet_P(V, V_0, K_0, Kprime_0):
    """
    Calculate pressure at room temperature with Vinet EOS.

    Parameters
    ----------
    V : float, array-like
        volume in cm3/mol
    V_0 : float, array-like
        volume at reference conditions
    K : float, array-like
        bulk modulus in GPa
    Kprime_0
        pressure derivative of K in GPa-1

    Returns
    -------
    P_GPa : float, array-like
        pressure in GPa
    """

    x = (V / V_0) ** (1 / 3)
    P_GPa = 3 * K_0 * x**-2 * (1 - x) * np.exp(1.5 * (Kprime_0 - 1) * (1 - x))
    return P_GPa


def Vinet_thermal(V, V_0, alpha0, delta0, kappa, T_K, Tref):
    """
    Calculate thermal contribution to volume with Vinet EOS.


    Parameters
    ----------
    V : float, array-like
        volume in cm3/mol at reference temperature Tref
    V_0 : float
        volume in cm3/mol at reference conditions
    alpha0 : float
        thermal expansian at 1 bar
    delta0 : float
        Value of Anderson-Grüneisen parameter (delta_T) at 1 bar
    kappa : float
        dimensionless Anderson-Grüneisen parameter.
    T_K  : float, array-like
        temperature in Kelvin
    Tref : float
        reference temperature in Kelvin

    Returns
    -------
    V   : float, array-like
        volume in cm3/mol
    """

    # calculate thermal expansion
    alpha = thermal_expansivity(V=V, V_0=V_0, alpha0=alpha0, delta0=delta0, kappa=kappa)
    # calculate volume at temperature T
    V = V * np.exp(alpha * (T_K - Tref))

    return V


@np.vectorize(excluded=["bracket"])
def Vinet_V_roomtemperature(P_GPa, V_0, K_0, Kprime_0):
    """
    Calculate volume at pressure P_GPa and room temperature with Vinet EOS.

    Parameters
    ----------
    P_GPa   : float, array-like
        pressure in GPa
    V_0 : float, array-like
        volume in cm3/mol at reference conditions
    K_0 : float
        bulk modulus at reference conditions
    Kprime_0
        pressure derivative of K in GPa-1

    Returns
    -------
    V   : float, array-like
        volume in cm3/mol
    """

    Vinet_solve = lambda V: Vinet_P(V=V, V_0=V_0, K_0=K_0, Kprime_0=Kprime_0) - P_GPa

    return opt.fsolve(Vinet_solve, x0=V_0)


def Vinet_V(P_GPa, T_K, V_0, K_0, Kprime_0, delta0, alpha0, kappa, Tref=298.15):
    """
    Calculate volume at pressure P_GPa and temperature T_K with Vinet EOS.

    Parameters
    ----------
    P_GPa : float, array-like
        pressure in GPa
    T_K  : float, array-like
        temperature in Kelvin
    V_0 : float
        volume in cm3/mol at reference conditions
    K_0 : float
        bulk modulus at reference conditions
    Kprime_0 : float
        pressure derivative of K_0
    delta0 : float
        Value of Anderson-Grüneisen parameter (delta_T) at 1 bar
    alpha0 : float
        thermal expansian at 1 bar
    kappa : float
        dimensionless Anderson-Grüneisen parameter.
    Tref : float
        reference temperature in Kelvin

    Returns
    -------
    V   : float, array-like
        volume in cm3/mol
    """

    # calculate volume at pressure P_GPa and room temperature
    V_T0_P = Vinet_V_roomtemperature(P_GPa=P_GPa, V_0=V_0, K_0=K_0, Kprime_0=Kprime_0)

    # calculate thermal part
    V = Vinet_thermal(
        V=V_T0_P, V_0=V_0, alpha0=alpha0, delta0=delta0, kappa=kappa, T_K=T_K, Tref=Tref
    )

    return V


def Vinet_VdP(
    P_GPa, T_K, V_0, K_0, Kprime_0, alpha0, delta0, kappa, n_step=100, **kwargs
):
    """
    Parameters
    ----------
    P_GPa : float, array-like
        pressure in GPa
    T_K  : float, array-like
        temperature in Kelvin
    V_0 : float
        volume at reference conditions
    K_0 : float
        bulk modulus at reference conditions
    Kprime_0 : float
        pressure derivative of K_0
    alpha0 : float
        thermal expansian at 1 bar
    delta0 : float
        Value of Anderson-Grüneisen parameter (delta_T) at 1 bar
    kappa : float
        dimensionless Anderson-Grüneisen parameter.
    n_step : float
        number of steps between 1 bar and P_GPa in the numerical integration of V

    Returns
    -------
    VdP
        volume integrated to from 1 bar to P_GPa in J/mol
    """
    if P_GPa <= 1e-4:
        return 0

    P = np.linspace(start=1e-4, stop=P_GPa, num=n_step)
    V = Vinet_V(
        P_GPa=P,
        T_K=T_K,
        V_0=V_0,
        K_0=K_0,
        Kprime_0=Kprime_0,
        delta0=delta0,
        alpha0=alpha0,
        kappa=kappa,
    )

    VdP = np.trapz(V, P)  # (cm3/mol)*GPa = (1000J/GPa/mol)*GPa

    return VdP * 1000  # convert to J/mol
