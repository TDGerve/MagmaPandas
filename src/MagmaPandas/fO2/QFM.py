import warnings as w

import numpy as np
import pandas as pd
import scipy.optimize as opt
from scipy.constants import R

from MagmaPandas.EOSs.tait import tait_VdP
from MagmaPandas.EOSs.tools import landau_P_dependent, phase_transition
from MagmaPandas.parse_io import repeat_vars


def _VdP_QFM(T_K, P_bar):
    """
    Solve Tait equations of state for VdP of quartz, magnetite and fayalite
    at temperature, T_K and P_bar, ignoring phase transitions.
    """

    p_kbar = P_bar / 1e3

    # VdP SiO2
    VdP_qtz = tait_VdP(phase="quartz", pkbar=p_kbar, T_K=T_K)
    Gibbs_landau = landau_P_dependent(phase="quartz", pkbar=p_kbar, T_K=T_K)
    VdP_qtz = VdP_qtz + Gibbs_landau

    VdP_fay = tait_VdP(phase="fayalite", pkbar=p_kbar, T_K=T_K)

    VdP_mt = tait_VdP(phase="magnetite", pkbar=p_kbar, T_K=T_K)

    return VdP_qtz, VdP_mt, VdP_fay


def _VdP_QFM_phaseTransitions(T_K, P_bar):
    """
    Solve Tait equations of state for VdP for quartz, magnetite and fayalite
    at temperature, T_K and P_bar, taking into account phase transitions
    """

    p_kbar = P_bar / 1e3

    if isinstance(p_kbar, pd.Series):
        p_kbar = p_kbar.squeeze()
    if isinstance(T_K, pd.Series):
        T_K = T_K.squeeze()

    # Calculate pressure of transition for SiO2 polymorphs
    # Quartz --> coesite
    qtz_coe = lambda P: phase_transition(
        P, T_K=T_K, phase_1="quartz", phase_2="coesite"
    )
    P_qtz_coe = opt.fsolve(qtz_coe, 8)
    # Coesite --> stishovite
    coe_stish = lambda P: phase_transition(
        P, T_K=T_K, phase_1="coesite", phase_2="stishovite"
    )
    P_coe_stish = opt.fsolve(coe_stish, 8)

    # Calculate pressure of transition for polymorphs
    # Fayalite --> ringwoodite
    fay_ring = lambda P: phase_transition(
        P, T_K=T_K, phase_1="fayalite", phase_2="ringwoodite"
    )
    P_fay_ring = opt.fsolve(fay_ring, 8)

    # Integrate VdP for all polymorphs
    # SiO2 polymorphs
    # Quartz
    VdP_SiO2 = tait_VdP(phase="quartz", pkbar=min(p_kbar, P_qtz_coe), T_K=T_K)
    SiO2_landau = landau_P_dependent(
        phase="quartz", pkbar=min(p_kbar, P_qtz_coe), T_K=T_K
    )
    VdP_SiO2 = VdP_SiO2 + SiO2_landau
    if p_kbar > P_qtz_coe:
        # Coesite
        VdP_coe = tait_VdP(
            phase="coesite", pkbar=min(p_kbar, P_coe_stish), T_K=T_K
        ) - tait_VdP(phase="coesite", pkbar=P_qtz_coe, T_K=T_K)
        VdP_SiO2 = VdP_SiO2 + VdP_coe
        if p_kbar > P_coe_stish:
            # Stishovite
            VdP_stish = tait_VdP(phase="stishovite", pkbar=p_kbar, T_K=T_K) - tait_VdP(
                phase="stishovite", pkbar=P_coe_stish, T_K=T_K
            )
            VdP_SiO2 = VdP_SiO2 + VdP_stish

    # Fe2SiO4 polymorphs
    # Fayalite
    VdP_Fe2SiO4 = tait_VdP(phase="fayalite", pkbar=min(p_kbar, P_fay_ring), T_K=T_K)
    if p_kbar > P_fay_ring:
        # Ringwoodite
        VdP_ring = tait_VdP(phase="ringwoodite", pkbar=p_kbar, T_K=T_K) - tait_VdP(
            phase="ringwoodite", pkbar=P_fay_ring, T_K=T_K
        )
        VdP_Fe2SiO4 = VdP_Fe2SiO4 + VdP_ring

    # Magnetite
    VdP_mt = tait_VdP(phase="magnetite", pkbar=p_kbar, T_K=T_K)

    return VdP_SiO2, VdP_mt, VdP_Fe2SiO4


def _muO2_QFM_P(T_K, P_bar):
    """
    calculate chemical potential of oxygen at QFM and pressure P with equations of state
    """

    P_bar_is_int = True if isinstance(P_bar, (int, float)) else False
    T_K_is_int = True if isinstance(T_K, (int, float)) else False

    # If P and T are both single values
    if P_bar_is_int and T_K_is_int:
        VdP_quartz, VdP_magnetite, VdP_fayalite = _VdP_QFM_phaseTransitions(T_K, P_bar)
        muO2 = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)
        return muO2

    # repeat the short variable
    T_K, P_bar = repeat_vars(T_K, P_bar)

    muO2 = np.array([])

    for temperature, pressure in zip(T_K, P_bar):
        VdP_quartz, VdP_magnetite, VdP_fayalite = _VdP_QFM_phaseTransitions(
            temperature, pressure
        )
        # kiloJoule to Joule
        muO2 = np.append(
            muO2, 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)
        )

    return muO2


def _muO2_QFM_1bar(T_K, warning=False):
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

    if warning:
        if (np.array(T_K) < 900).any():
            w.warn("O'Neill fO2: temperatures below 900K present")
        if (np.array(T_K) > 1420).any():
            w.warn("O'Neill fO2: temperatures above 1420K present")

    muO2 = -587474 + 1584.427 * T - 203.3164 * T * np.log(T) + 0.092710 * T**2

    # Convert to float if there is only one item
    if len(muO2) == 1:
        muO2 = muO2.item()

    return muO2


def _fO2_QFM_1bar(T_K, logshift=0):
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
    mu_O2 = _muO2_QFM_1bar(T_K)

    offset = 10**logshift

    return np.exp(mu_O2 / (R * T_K)) * offset


def calculate_fO2(
    logshift: int | float, T_K: float | np.ndarray, P_bar: float | np.ndarray
) -> float | np.ndarray:
    """
    Calculate |fO2| at the QFM buffer.

    1 bar components is calculated according to O'Neill (1987)\ [13]_ and pressure contributions according to Holland and Powell (2011)\ [14]_, with Landau theory from Holland and Powell (1990, 1998)\ [15]_:sup:`,`\ [16]_ and thermal Tait equations of state parameters from Holland and Powell (2011)\ [14]_, updated by Jennings and Holland (2015)\ [17]_.

    Parameters
    ----------
    logshift    : int, float
        |fO2| buffer shift in log units of QFM.
    T_K : float, array-like
        temperatures in Kelvin
    P_bar : float, array-like
        pressure in bar

    Returns
    -------
    fO2 : float, array-like
        |fO2| in bar
    """
    try:
        logshift = float(logshift)
    except TypeError:
        pass

    offset = 10**logshift

    # Chemical potential of oxygen
    # 1 bar contribution from O'Neill
    muO2_1bar_Oneill = _muO2_QFM_1bar(T_K)
    # Pressue contribution from equations of state
    muO2_pressure_eos = _muO2_QFM_P(T_K, P_bar)
    # Remove any 1 bar contribution from the equation of state formulation,
    # since the O'Neill emperical formulation is used for this
    VdP_quartz, VdP_magnetite, VdP_fayalite = _VdP_QFM(T_K, 1)
    muO2_1bar_eos = 1e3 * (3 * VdP_quartz + 2 * VdP_magnetite - 3 * VdP_fayalite)
    muO2_pressure = muO2_pressure_eos - muO2_1bar_eos
    # Total chemical potential
    muO2 = muO2_1bar_Oneill + muO2_pressure

    fO2 = np.exp(muO2 / (R * T_K)) * offset

    try:
        fO2 = np.float32(fO2.item())
    except ValueError:
        fO2 = fO2.astype(np.float32)

    # if isinstance(fO2, pd.Series):
    #     return fO2.squeeze()

    return fO2
