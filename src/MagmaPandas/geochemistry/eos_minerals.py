import numpy as np
import scipy.optimize as opt


class EOSparams:
    """
    h                       enthalpy of formation
    s                       entropy
    v0                      volume at 1bar, 298K
    n                       atoms per formula unit
    a0                      coefficient of thermal expansion
    K0                       bulk modulus at 1bar, 298K
    dKdP                    first derivative of K
    dKdP2                   second derivative of K
    cp_a, cp_b, cp_c, cp_d  coefficients of heat capacity polynomial: a + bT + cT**-2 + dT**-(1/2)
    Tc0                     Landau critical temperature at 1 bar
    vmax                    maximum volume of disorder
    smax                    maximum entropy of disorder

    data from Holland & Powell (2011)
    """

    fayalite = {
        "h": -1477.510,
        "s": 151.0,
        "v0": 4.631,
        "n": 7,
        "a0": 2.82e-5,
        "K0": 1256,
        "dKdP": 4.68,
        "dKdP2": -3.7e-3,
        "cp_a": 2.011e-1,
        "cp_b": 1.733e-5,
        "cp_c": -1960.6,
        "cp_d": -9.009e-1,
    }

    ringwoodite = {
        "h": -1477.510,
        "s": 140.0,
        "v0": 4.203,
        "n": 7,
        "a0": 2.22e-5,
        "K0": 1977,
        "dKdP": 4.92,
        "dKdP2": -2.5e-3,
        "cp_a": 1.668e-1,
        "cp_b": 4.2610e-5,
        "cp_c": -1705.4,
        "cp_d": -5.414e-1,
    }

    quartz = {
        "h": -910.710,
        "s": 41.43,
        "v0": 2.269,
        "n": 3,
        "a0": 0,
        "K0": 730,
        "dKdP": 6,
        "dKdP2": -8.2e-3,
        "smax": 4.95 / 1e3,
        "vmax": 1.188e-1,
        "Tc0": 847,
        "cp_a": 9.29e-2,
        "cp_b": -6.42e-7,
        "cp_c": -714.9,
        "cp_d": -0.7161,
    }

    coesite = {
        "h": -906.990,
        "s": 39.60,
        "v0": 2.064,
        "n": 3,
        "a0": 1.23e-5,
        "K0": 979,
        "dKdP": 4.19,
        "dKdP2": -4.3e-3,
        "cp_a": 1.078e-1,
        "cp_b": -3.279e-6,
        "cp_c": -190.3,
        "cp_d": -1.0416,
    }

    stishovite = {
        "h": -876.720,
        "s": 24.0,
        "v0": 1.401,
        "n": 3,
        "a0": 1.58e-5,
        "K0": 3090,
        "dKdP": 4.6,
        "dKdP2": -1.50e-3,
        "cp_a": 6.81e-2,
        "cp_b": 6.010e-6,
        "cp_c": -1978.2,
        "cp_d": -8.21e-2,
    }

    magnetite = {
        "h": -1114.510,
        "s": 146.9,
        "v0": 4.452,
        "n": 7,
        "a0": 3.71e-5,
        "K0": 1857,
        "dKdP": 4.05,
        "dKdP2": -2.2e-3,
        "smax": 35.0,
        "vmax": 0.0,
        "Tc0": 848,
    }


def tait_eos_pressure(phase, pkbar, t, tref=298.15, **kwargs):
    """    
    Pressure contribution to Gibb's free energy.
    Tait equation of state from Holland and Powell (2011)

    Parameters
    ----------
    phase       str
        Mineral phase
    pkbar       int, float
        Pressure in kilobar
    t           int, float
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

    u0 = theta / tref
    u = theta / t

    # Einstein function page 345
    xi0 = u0 ** 2 * np.exp(u0) / (np.exp(u0) - 1) ** 2.0

    # Equation 3
    a = (1.0 + dKdP) / (1.0 + dKdP + K0 * dKdP2)
    b = dKdP / K0 - dKdP2 / (1.0 + dKdP)
    c = (1.0 + dKdP + K0 * dKdP2) / (dKdP ** 2.0 + dKdP - K0 * dKdP2)

    # thermal pressure term, equation 11
    Pth = a0 * K0 * theta / xi0 * (1 / (np.exp(u) - 1.0) - 1 / (np.exp(u0) - 1.0))

    # Intergral of volume, equation 13
    PV0 = pkbar * v0
    part1 = np.sign(1 - b * Pth) * abs(1 - b * Pth) ** (1 - c)
    part2 = np.sign(1 + b * (pkbar - Pth)) * abs(1 + b * (pkbar - Pth)) ** (1 - c)
    part3 = b * (c - 1) * pkbar

    VdP = PV0 * (1 - a + a * (part1 - part2) / part3)

    return VdP


def enthalpy(phase, t, tref=298.15):
    """
    Enthalpy as Cp.dT, integrated from T to T(reference), with heat capacity as:

    Cp = a + bT + cT**-2 + dT**(-1/2)

    Cp, a, b, c, & d from Holland and Powell (2011)

    Parameters
    ----------
    phase   str

    t       int, float
        temperature in Kelvin

    tref    int, float
        reference temperature in Kelvin

    Returns
    -------
    enthalpy   float
    """

    a, b, c, d = [
        getattr(EOSparams, phase)[i] for i in ["cp_a", "cp_b", "cp_c", "cp_d"]
    ]

    integral = lambda T: a * T + 0.5 * b * T**2 - c * T**-1 + 2 * d * T**(1/2)
    enthalpy = integral(t) - integral(tref)

    return enthalpy


def entropy(phase, t, tref=298.15):
    """
    Entropy as (Cp/T)dT, integrated from T to T(reference), with:

    Cp/T = a/T + b + cT**-3 + dT**(-3/2)

    Cp, a, b, c & d from Holland and Powell (2011)

    Parameters
    ----------
    phase   str

    t       int, float
        temperature in Kelvin

    tref    int, float
        reference temperature in Kelvin

    Returns
    -------
    entropy   float
    """

    a, b, c, d = [
        getattr(EOSparams, phase)[i] for i in ["cp_a", "cp_b", "cp_c", "cp_d"]
    ]

    integral = lambda T: a * np.log(T) + b * T - c / 2 * T**-2 - 2 * d * T**(-1/2)
    entropy = integral(t) - integral(tref)

    return entropy


def landau(phase, pkbar, T_K, **kwargs):
    """
    Excess Gibbs free energy from Landau theory.

    Based on:
    Holland and Powell, 1998, p.312
    and
    Holland and Powell, 1990

    Rewritten according to Michael Anenberg:
    https://fo2.rses.anu.edu.au/fo2app/
    (see FMQ buffer details and references)

    Parameters
    ----------
    phase   str
        Mineral phase
    pkbar   int, float
        Pressure in kilobar
    t       int, float
        Temperature in Kelvin
    vmax    int, float
        Optional, maximum volume of disorder. Set to 0 to calculate the pressure independent contribution.

    Returns
    -------
    float
        Pressure dependent contribution to the excess Gibbs free enery from Landau theory
    """
    t = np.array([])
    t = np.append(t, T_K)

    vmax_default = getattr(EOSparams, phase)["vmax"] 
    vmax = kwargs.get("vmax", vmax_default)

    smax, tc0 = [getattr(EOSparams, phase)[i] for i in ["smax", "Tc0"]]

    Q2_0 = np.sqrt(1 - 298.15 / tc0)

    # Landau critical temperature at pkbar
    tc = tc0 + pkbar * vmax / smax

    Q2 = np.zeros(shape=[len(t),])
    if any(t < tc):
        Q2 = np.where(t > tc, 0, np.sqrt((tc - t) / tc0))

    G_Landau = (
        smax * (tc0 * (Q2_0 + (Q2 ** 3 - Q2_0 ** 3) / 3) - tc * Q2 - t * (Q2_0 - Q2))
        + pkbar * vmax * Q2_0
    )

    # Convert to float if there is only one item
    if len(G_Landau) == 1:
        G_Landau = G_Landau.item()

    return G_Landau


def landau_P_dependent(phase, pkbar, t, formulation='anenberg'):
    """
    Pressure dependent excess Gibbs free energy from Landau theory.
    Calculated by subtracting the pressure independent contribution (with vmax = 0) from the total.

    Parameters
    ----------
    phase       :   str
        Mineral phase
    t           :   int, float
        Temperature in Kelvin
    pkbar       :   int, float
        Pressure in kilobars
    formulation :   str
        'holland' for the formulation by Holland and Powell (1998), otherwise Michael Anenberg's
        formulation will be used: 
        https://fo2.rses.anu.edu.au/fo2app/ 
        (see FMQ buffer details and references)
    Returns
    -------
    float
        Pressure contribution to excess Gibbs free enery from Landau theory
    """
    if formulation == 'holland':
        landau_total = landau_Holland(phase, pkbar, t)
        landau_1bar = landau_Holland(phase, 0, t, vmax=0)
    else:
        landau_total = landau(phase, pkbar, t)
        landau_1bar = landau(phase, 0, t, vmax=0)

    return landau_total - landau_1bar


def landau_Holland(phase, pkbar, T_K, **kwargs):
    """
    Excess Gibbs free energy from Landau theory

    Equations from Holland and Powell (1998), p. 312

    Parameters
    ----------
    phase   str
        Mineral phase
    pkbar   int, float
        Pressure in kilobars
    t       int, float
        Temperature in Kelvin

    Returns
    -------
    float
        Excess Gibbs free enery from Landau theory
    """
    t = np.array([])
    t = np.append(t, T_K)

    vmax_default = getattr(EOSparams, phase)["vmax"] 
    vmax = kwargs.get("vmax", vmax_default)
    smax, tc0, a0, K0 = [
        getattr(EOSparams, phase)[i] for i in ["smax", "Tc0", "a0", "K0"]
    ]

    # Landau critical temperature
    tc = tc0 + vmax * pkbar / smax
    # Q: oder paramter in the landau model
    Q2_0 = np.sqrt(1 - 298.15 / tc0)
    
    # if t > tc:
    #     Q2 = 0
    # else:
    #     Q2 = np.sqrt((tc - t) / tc0)

    Q2 = np.where(t > tc, 0, np.sqrt((tc - t) / tc0))

    # Bulk modulus
    K = K0 * (1 - 1.5e-4 * (t - 298))

    # Excess enthalpy at 298K from Landau model disordering
    h = smax * tc0 * (Q2_0 - (Q2_0 ** 3) / 3)
    # Excess entropy at 298K from Landau model disordering
    s = smax * Q2_0

    # Excess volume from Landau model disordering
    vt = vmax * Q2_0 * (1 + a0 * (t - 298)) - 20 * a0 * (np.sqrt(t) - np.sqrt(298))

    vtdP = vt * K / 3 * ((1 + 4 * pkbar / K)**(3 / 4) - 1)

    delta_G_landau = smax * ((t - tc0) * Q2 + (tc * Q2**3) / 3)

    G_excess = h - t * s + vtdP + delta_G_landau

    # Convert to float if there is only one item
    if len(G_excess) == 1:
        G_excess = G_excess.item()

    return G_excess


def phaseTransition(pkbar, t, phase_1, phase_2):
    """
    Gibbs free energy of transition from phase_1 to phase_2

    Parameters
    ----------
    pkbar               int, float
        Pressure in kilobar
    t                   int, float
        Temperature in Kelvin
    phase_1, phase_2    str
        Mineral phases

    Returns
    -------
    float
        Gibbs free energy of phase transition
    """

    results = []

    for phase in [phase_1, phase_2]:

        h = getattr(EOSparams, phase)["h"]
        s = getattr(EOSparams, phase)["s"] / 1e3

        Gibbs = h + enthalpy(phase=phase, t=t) - t * (s + entropy(phase=phase, t=t))
        VdP = tait_eos_pressure(phase=phase, pkbar=pkbar, t=t)

        Gibbs = Gibbs + VdP

        if phase in ["quartz", "magnetite"]:

            Gibbs = Gibbs + landau(phase=phase, pkbar=pkbar, T_K=t)

        results.append(Gibbs)

    return results[0] - results[1]
