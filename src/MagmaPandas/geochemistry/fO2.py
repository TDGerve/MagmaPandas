from . import eos_minerals
import scipy.optimize as opt
import warnings as w
import numpy as np
import pandas as pd
from scipy.constants import R
import itertools as it


def VdP_QFM(T_K, P_bar):
    """
    Solve Tait equations of state for VdP of quartz, magnetite and fayalite
    at temperature, T_K and P_bar, ignoring phase transitions.
    """

    p_kbar = P_bar / 1e3

    # VdP SiO2
    VdP_qtz = eos_minerals.tait_eos_pressure(phase="quartz", pkbar=p_kbar, t=T_K)
    Gibbs_landau = eos_minerals.landau_P_dependent(phase="quartz", pkbar=p_kbar, t=T_K)
    VdP_qtz = VdP_qtz + Gibbs_landau

    VdP_fay = eos_minerals.tait_eos_pressure(phase="fayalite", pkbar=p_kbar, t=T_K)

    VdP_mt = eos_minerals.tait_eos_pressure(phase="magnetite", pkbar=p_kbar, t=T_K)

    return VdP_qtz, VdP_mt, VdP_fay


def VdP_QFM_phaseTransitions(T_K, P_bar):
    """
    Solve Tait equations of state for VdP for quartz, magnetite and fayalite
    at temperature, T_K and P_bar, taking into account phase transitions
    """

    p_kbar = P_bar / 1e3

    # Calculate pressure of transition for SiO2 polymorphs
    # Quartz --> coesite
    qtz_coe = lambda P: eos_minerals.phaseTransition(
        P, t=T_K, phase_1="quartz", phase_2="coesite"
    )
    P_qtz_coe = opt.fsolve(qtz_coe, 8)
    # Coesite --> stishovite
    coe_stish = lambda P: eos_minerals.phaseTransition(
        P, t=T_K, phase_1="coesite", phase_2="stishovite"
    )
    P_coe_stish = opt.fsolve(coe_stish, 8)

    # Calculate pressure of transition for polymorphs
    # Fayalite --> ringwoodite
    fay_ring = lambda P: eos_minerals.phaseTransition(
        P, t=T_K, phase_1="fayalite", phase_2="ringwoodite"
    )
    P_fay_ring = opt.fsolve(fay_ring, 8)

    # Integrate VdP for all polymorphs
    # SiO2 polymorphs
    # Quartz
    VdP_SiO2 = eos_minerals.tait_eos_pressure(
        phase="quartz", pkbar=min(p_kbar, P_qtz_coe), t=T_K
    )
    SiO2_landau = eos_minerals.landau_P_dependent(
        phase="quartz", pkbar=min(p_kbar, P_qtz_coe), t=T_K
    )
    VdP_SiO2 = VdP_SiO2 + SiO2_landau
    if p_kbar > P_qtz_coe:
        # Coesite
        VdP_coe = eos_minerals.tait_eos_pressure(
            phase="coesite", pkbar=min(p_kbar, P_coe_stish), t=T_K
        ) - eos_minerals.tait_eos_pressure(phase="coesite", pkbar=P_qtz_coe, t=T_K)
        VdP_SiO2 = VdP_SiO2 + VdP_coe
        if p_kbar > P_coe_stish:
            # Stishovite
            VdP_stish = eos_minerals.tait_eos_pressure(
                phase="stishovite", pkbar=p_kbar, t=T_K
            ) - eos_minerals.tait_eos_pressure(
                phase="stishovite", pkbar=P_coe_stish, t=T_K
            )
            VdP_SiO2 = VdP_SiO2 + VdP_stish

    # Fe2SiO4 polymorphs
    # Fayalite
    VdP_Fe2SiO4 = eos_minerals.tait_eos_pressure(
        phase="fayalite", pkbar=min(p_kbar, P_fay_ring), t=T_K
    )
    if p_kbar > P_fay_ring:
        # Ringwoodite
        VdP_ring = eos_minerals.tait_eos_pressure(
            phase="ringwoodite", pkbar=p_kbar, t=T_K
        ) - eos_minerals.tait_eos_pressure(phase="ringwoodite", pkbar=P_fay_ring, t=T_K)
        VdP_Fe2SiO4 = VdP_Fe2SiO4 + VdP_ring

    # Magnetite
    VdP_mt = eos_minerals.tait_eos_pressure(phase="magnetite", pkbar=p_kbar, t=T_K)

    return VdP_SiO2, VdP_mt, VdP_Fe2SiO4


def muO2_QFM_P(T_K, P_bar):
    """
    calculate chemical potential of oxygen at QFM and pressure P with equations of state
    """

    P_bar_is_int = isinstance(P_bar, (int, float))
    T_K_is_int = isinstance(T_K, (int, float))

    # If P and T are not both numbers
    if not (P_bar_is_int and T_K_is_int):

        # If only one variable, P or T, is a number
        if bool(P_bar_is_int) ^ bool(T_K_is_int):

            # Cycle the short variable
            T_K = [np.array(T_K), it.cycle(np.array([T_K]))][T_K_is_int]
            P_bar = [np.array(P_bar), it.cycle(np.array([P_bar]))][P_bar_is_int]

        try:
            length = len(T_K)
        except TypeError:
            length = len(P_bar)

        muO2 = np.zeros(
            shape=[
                length,
            ]
        )

        for i, (temperature, pressure) in enumerate(zip(T_K, P_bar)):
            VdP_quartz, VdP_magnetite, VdP_fayalite = VdP_QFM_phaseTransitions(
                temperature, pressure
            )
            # kiloJoule to Joule
            muO2[i] = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)

    else:
        VdP_quartz, VdP_magnetite, VdP_fayalite = VdP_QFM_phaseTransitions(T_K, P_bar)
        muO2 = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)

    return muO2


def muO2_QFM_1bar(T_K):
    """
    calculate chemical potential of oxygen at QFM a 1 bar. Equation from O'Neill 1987

    Parameters
    ----------
    T_K     list-like, float
        Temperature in Kelvin

    Returns
    -------
    Chemical potential
    """

    T = np.array([])
    T = np.append(T, T_K)

    if (np.array(T_K) < 900).any():
        w.warn("O'Neill fO2: temperatures below 900K present")
    if (np.array(T_K) > 1420).any():
        w.warn("O'Neill fO2: temperatures above 1420K present")

    muO2 = -587474 + 1584.427 * T - 203.3164 * T * np.log(T) + 0.092710 * T ** 2

    # Convert to float if there is only one item
    if len(muO2) == 1:
        muO2 = muO2.item()

    return muO2


def fO2_QFM_1bar(T_K, logshift=0):
    """
    calculate fO2 at QFM + logshift a 1 bar. Equation from O'Neill 1987

    Parameters
    ----------
    T_K     list-like, float
        Temperature in Kelvin
    logshift   int, float
        Log units by which QFM is shifted

    Returns
    -------
    fO2
    """
    mu_O2 = muO2_QFM_1bar(T_K)

    offset = 10 ** logshift

    return np.exp(mu_O2 / (R * T_K)) * offset


def fO2_QFM(logshift, T_K, P_bar):
    """ """
    offset = 10 ** logshift

    muO2_pressure = muO2_QFM_P(T_K, P_bar) - muO2_QFM_P(T_K, 1)

    muO2_1bar = muO2_QFM_1bar(T_K)

    return np.exp((muO2_1bar + muO2_pressure) / (R * T_K)) * offset
