import numpy as np

from MagmaPandas.EOSs.parameters import EOSparams
from MagmaPandas.EOSs.tait import tait_VdP


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

    Q2 = np.zeros(
        shape=[
            len(t),
        ]
    )
    if any(t < tc):
        Q2 = np.where(t > tc, 0, np.sqrt((tc - t) / tc0))

    G_Landau = (
        smax * (tc0 * (Q2_0 + (Q2**3 - Q2_0**3) / 3) - tc * Q2 - t * (Q2_0 - Q2))
        + pkbar * vmax * Q2_0
    )

    # Convert to float if there is only one item
    if len(G_Landau) == 1:
        G_Landau = G_Landau.item()

    return G_Landau


def landau_P_dependent(phase, pkbar, T_K, formulation="anenberg"):
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
    if formulation == "holland":
        landau_total = landau_Holland(phase, pkbar, T_K)
        landau_1bar = landau_Holland(phase, 0, T_K, vmax=0)
    else:
        landau_total = landau(phase, pkbar, T_K)
        landau_1bar = landau(phase, 0, T_K, vmax=0)

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
    h = smax * tc0 * (Q2_0 - (Q2_0**3) / 3)
    # Excess entropy at 298K from Landau model disordering
    s = smax * Q2_0

    # Excess volume from Landau model disordering
    vt = vmax * Q2_0 * (1 + a0 * (t - 298)) - 20 * a0 * (np.sqrt(t) - np.sqrt(298))

    vtdP = vt * K / 3 * ((1 + 4 * pkbar / K) ** (3 / 4) - 1)

    delta_G_landau = smax * ((t - tc0) * Q2 + (tc * Q2**3) / 3)

    G_excess = h - t * s + vtdP + delta_G_landau

    # Convert to float if there is only one item
    if len(G_excess) == 1:
        G_excess = G_excess.item()

    return G_excess


def phase_transition(pkbar, T_K, phase_1, phase_2):
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

        Gibbs = (
            h
            + enthalpy(phase=phase, T_K=T_K)
            - T_K * (s + entropy(phase=phase, T_K=T_K))
        )
        VdP = tait_VdP(phase=phase, pkbar=pkbar, T_K=T_K)

        Gibbs = Gibbs + VdP

        if phase in ["quartz", "magnetite"]:

            Gibbs = Gibbs + landau(phase=phase, pkbar=pkbar, T_K=T_K)

        results.append(Gibbs)

    return results[0] - results[1]


def enthalpy(phase: str, T_K: int | float, Tref: int | float = 298.15):
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

    integral = lambda T: a * T + 0.5 * b * T**2.0 - c * T**-1.0 + 2 * d * T ** (1 / 2)
    enthalpy = integral(T_K) - integral(Tref)

    return enthalpy


def entropy(phase: str, T_K: int | float, Tref: int | float = 298.15):
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

    integral = lambda T: a * np.log(T) + b * T - c / 2 * T**-2.0 - 2 * d * T ** (-1 / 2)
    entropy = integral(T_K) - integral(Tref)

    return entropy


def thermal_expansivity(V, V_0, alpha0, delta0, kappa):
    """
    Equation 5 from Komabayashi (2014)

    Parameters
    ----------
    V : float, array-like
        volume
    V_0 : float
        volume at reference conditions
    alpha0 : float
        thermal expansion at 1 bar
    delta0 : float
        Value of Anderson-Grüneisen parameter (delta_T) at 1 bar
    kappa : float
        dimensionless Anderson-Grüneisen parameter.


    Returns
    -------
    alpha
        thermal expansivity
    """

    return alpha0 * np.exp((-delta0 / kappa) * (1 - (V / V_0) ** kappa))
