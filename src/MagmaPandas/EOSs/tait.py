import numpy as np

from MagmaPandas.EOSs.parameters import EOSparams


def tait_VdP(phase, pkbar, T_K, Tref=298.15, **kwargs):
    """
    Pressure contribution to Gibb's free energy.
    Tait equation of state from Holland and Powell (2011)

    Parameters
    ----------
    phase       str
        Mineral phase
    pkbar       int, float
        Pressure in kilobar
    T_K           int, float
        Temperature in Kelvin
    tref        int, float
        Reference temperature in Kelvin

    Returns
    Pth         float
        Thermal pressure term
    a, b, c     float
        Equation of state parameters
    VdP         float
        Contribution of pressure to the Gibbs free energy
    """

    params = ["s", "v0", "n", "a0", "K0", "dKdP", "dKdP2"]
    s, v0, n, a0, K0, dKdP, dKdP2 = [getattr(EOSparams, phase)[i] for i in params]

    # Einstein temperature
    theta = 10636.0 / (s / n + 6.44)

    u0 = theta / Tref
    u = theta / T_K

    # Einstein function page 345
    xi0 = u0**2 * np.exp(u0) / (np.exp(u0) - 1) ** 2.0

    # Equation 3
    a = (1.0 + dKdP) / (1.0 + dKdP + K0 * dKdP2)
    b = dKdP / K0 - dKdP2 / (1.0 + dKdP)
    c = (1.0 + dKdP + K0 * dKdP2) / (dKdP**2.0 + dKdP - K0 * dKdP2)

    # thermal pressure term, equation 11
    Pth = a0 * K0 * theta / xi0 * (1 / (np.exp(u) - 1.0) - 1 / (np.exp(u0) - 1.0))

    # Intergral of volume, equation 13
    PV0 = pkbar * v0
    part1 = np.sign(1 - b * Pth) * abs(1 - b * Pth) ** (1 - c)
    part2 = np.sign(1 + b * (pkbar - Pth)) * abs(1 + b * (pkbar - Pth)) ** (1 - c)
    part3 = b * (c - 1) * pkbar

    VdP = PV0 * (1 - a + a * (part1 - part2) / part3)

    return VdP
