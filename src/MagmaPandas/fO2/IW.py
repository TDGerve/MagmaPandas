from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.constants import R
from scipy.optimize import brentq

from MagmaPandas.EOSs.vinet import Vinet_VdP
from MagmaPandas.parse_io import make_iterable, repeat_vars

Fe_polymorphs = [
    "Fe_fcc",
    "Fe_bcc-alpha",
    "Fe_HCP",
    "Fe_bcc-delta",
    "Fe_liquid",
]

"""
Polynomial parameters for calculating Gibbs free energy at 1 bar and temperatures above 1000 Kelvin.
Values from Hidayat et al. 2015 for FeO and FeO1.5 and Dinsdale (1991) for O2
"""
# TODO check if data are copied via PyPI
_IW_G0_params_file = (
    files("MagmaPandas.fO2.data").joinpath("IW_G0_params.csv").open("r")
)
gibbs0_parameters = pd.read_csv(_IW_G0_params_file, index_col=["phase", "temperature"])

# Below 1000 Kelvin for O2:
O2_low_temperature = pd.Series(
    {
        "a": -6961.7445,  # constant
        "b": -51.0057,  # bT_K
        "c": -22.271,  # cT_Kln(T_K)
        "d": 0,  # d ln(T_K)
        "e": -1.01977e-2,  # dT_K**2
        "f": 1.32369e-8,  # fT_K**3
        "g": -7629.7484,  # g/T_K
        "h": 0.0,  # hT**7
        "i": 0.0,  # iT_K**9
        # "func": Gibbs0_polynomial_1,
    }
)


"""
Values from table Hirschmann et al. (2018) table S2

EOS is from Komabayashi, 2014, with V0 of FeO1.5 taken by extrapolation
of wustite stoichiometry-volume data compiled by Hirschmann et al. (2018, 16.372±0.070)
This is very similar to value from McCammon and Liu 1984 trend (16.425).
Fe parameters from Komabayashi (2014), with bcc the same as fcc, except 
that V0 of bcc is taken from Dorogokupets et al. (2017)

V_0     
    volume at reference conditions
K_0
    bulk modulus at reference conditions
Kprime_0
    first pressure derivative of K_0
alpha0
    thermal expansivity at 1 bar
delta0
    Value of Anderson-Grüneisen parameter (delta_T) at 1 bar
kappa
    dimensionless Anderson-Grüneisen parameter. Set to 1.4 following Komabayashi (2014) and Wood (1993)
"""
# TODO Move these parameters to EOSs.parameters? Make the parameter names uniform, e.g. Kprime_0 vs dKdP
eos_parameters = pd.DataFrame(
    {
        "V_0": [12.256, 16.372, 6.82, 7.092, 6.753, 7.092, 6.88],  # cm3/mol
        "K_0": [149, 149, 163.4, 163.4, 163.4, 163.4, 148],  # GPa
        "Kprime_0": [3.83, 3.83, 5.38, 5.38, 5.38, 5.38, 5.8],  # GPa-1
        "alpha0": [4.5e-5, 4.5e-5, 7e-05, 7e-05, 5.8e-05, 7e-05, 9e-5],
        "delta0": [4.25, 4.25, 5.5, 5.5, 5.1, 5.5, 5.1],
        "kappa": [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4],
    },
    index=[
        "FeO",
        "FeO1.5",
        "Fe_fcc",
        "Fe_bcc-alpha",
        "Fe_HCP",
        "Fe_bcc-delta",
        "Fe_liquid",
    ],
)


"""
Mixing parameters in J/mol for FeO-FeO1.5 solid solution from Hidayat et al. (2015)

RT ln gamma(X1) =((q00+2*q10*(1-X2))*X2**2)         eq. S5a
RT ln gamma(X2) =((1-X2)**2*(q00+q10-2*q10*X2))     eq. S5b

XFeO=X1; XFeO1.5=X2
"""
mixing_parameters = pd.Series({"q00": -5.94e4, "q10": 4.27e4})


def Gibbs0_polynomial(T_K: float | np.ndarray, a, b, c, d, e, f, g, h, i, **kwargs):

    return (
        a
        + b * T_K
        + c * T_K * np.log(T_K)
        + d * np.log(T_K)
        + e * T_K**2
        + f * T_K**3
        + g / T_K
        + h * T_K**7
        + i * T_K**-9.0
    )


def Gibbs0(T_K: float | np.ndarray, params: pd.Series | pd.DataFrame):
    """
    calculate Gibbs free energy at 1 bar
    """
    if isinstance(params, pd.DataFrame):
        func = lambda params: Gibbs0_polynomial(T_K=T_K, **params)
        return params.apply(func, axis=1)

    return Gibbs0_polynomial(T_K=T_K, **params)


def Gibbs_Fe_magnetic(T_K):
    """
    Adjust bcc-alpha for magnetic contribution.  Magnetic contribution for fcc
    is << 1 J at temperatures of interest and is ignored

    translated from Hirschmann et al. (2021) matlab script
    """
    Tc = 1043  # curie temperature
    P_factor = 0.4  # structure dependent parameter
    beta = 2.22  # magnetic moment
    A = 1.55828482

    if T_K < Tc:
        gibbs_magnetic = 1 - (1 / A) * (
            (79 * (T_K / Tc) ** (-1) / (140 * P_factor))
            + (474 / 497)
            * ((1 / P_factor) - 1)
            * ((T_K / Tc) ** 3 / 6 + (T_K / Tc) ** 9 / 135 + (T_K / Tc) ** 15 / 600)
        )
    else:
        # why divide T_K by T_K, should this be T_K / Tc?
        gibbs_magnetic = (-1 / A) * (
            (T_K / Tc) ** -5 / 10 + (T_K / Tc) ** -15 / 315 + (T_K / T_K) ** -25 / 1500
        )

    return gibbs_magnetic * (R * T_K * np.log(beta + 1))


def Gibbs_Fe_polymorphs(P_bar, T_K):
    """
    Calculate Gibbs free energy for Fe polymorphs

    translated from Hirschmann et al. (2021) matlab script
    """
    # TODO there is a tiny difference in pressure adjusted G, check where this comes from.

    P_GPa = P_bar / 1e4
    temperature_range = "high" if T_K > 1811 else "low"
    gibbs0_params = gibbs0_parameters.loc[(Fe_polymorphs, temperature_range), :]
    gibbs0 = pd.Series(index=Fe_polymorphs)
    gibbs_pressure = pd.Series(index=Fe_polymorphs)

    for ((name, _), G0), (_, eos_params) in zip(
        gibbs0_params.iterrows(), eos_parameters.loc[Fe_polymorphs, :].iterrows()
    ):
        gibbs0.loc[name] = Gibbs0(T_K=T_K, params=G0)
        gibbs_pressure.loc[name] = Vinet_VdP(P_GPa=P_GPa, T_K=T_K, **eos_params)

    gibbs_total = gibbs0.add(gibbs_pressure)
    gibbs_total["Fe_bcc-alpha"] += Gibbs_Fe_magnetic(T_K=T_K)

    return gibbs_total


def Gibbs_wustite_O2(P_bar, T_K):
    """
    Calculate gibbs free enery for FeO, FeO1.5 and O2

    translated from Hirschmann et al. (2021) matlab script
    """
    P_GPa = P_bar / 1e4
    temperature_range = "high" if T_K > 1811 else "low"

    gibbs0 = pd.Series(index=["FeO", "FeO1.5", "O2"])
    gibbs_pressure = pd.Series(0.0, index=["FeO", "FeO1.5"])

    gibbs0_params = gibbs0_parameters.loc[(gibbs0.index, temperature_range), :].copy()
    gibbs0_params.loc["O2"] = (
        O2_low_temperature.values
        if (T_K < 1000)
        else gibbs0_params.loc[("O2", temperature_range), :].values
    )

    for (name, _), G0 in gibbs0_params.iterrows():
        gibbs0.loc[name] = Gibbs0(T_K=T_K, params=G0)
    for name, eos_params in eos_parameters.loc[gibbs_pressure.index, :].iterrows():
        gibbs_pressure.loc[name] = Vinet_VdP(P_GPa=P_GPa, T_K=T_K, **eos_params)

    gibbs_total = gibbs0.add(gibbs_pressure, fill_value=0.0)

    return gibbs_total


def Gibbs_IW(P_bar, T_K, Fe_phase=False, suppress_Fe_liquid=True):
    """
    Calculate Gibbs free enery for FeO, FeO1.5, O2 and Fe. For Fe, the polymorph with lowest Gibbs free energy is used.

    translated from Hirschmann et al. (2021) matlab script
    """

    gibbs_total = Gibbs_wustite_O2(P_bar=P_bar, T_K=T_K)

    gibbs_Fe_all = Gibbs_Fe_polymorphs(P_bar=P_bar, T_K=T_K)
    addendum = ""
    if suppress_Fe_liquid:
        addendum = (
            "*" if gibbs_Fe_all.index[gibbs_Fe_all.argmin()] == "Fe_liquid" else ""
        )
        gibbs_Fe_all = gibbs_Fe_all.drop("Fe_liquid")
    Fe_phase, Fe_gibbs = (
        gibbs_Fe_all.index[gibbs_Fe_all.argmin()],
        gibbs_Fe_all.min(),
    )
    gibbs_total["Fe"] = Fe_gibbs

    output = (gibbs_total, Fe_phase + addendum) if Fe_phase else gibbs_total

    return output


def deltaGibbs_Fe_wustite(Gibbs: pd.Series):
    """
    Calculate delta Gibbs of the reaction 2 FeO1.5 +  Fe = 3 FeO

    translated from Hirschmann et al. (2021) matlab script

    deltaG : float
        Change of Gibbs free energy of the reaction
    """

    return (3 * Gibbs["FeO"]) - (2 * Gibbs["FeO1.5"]) - Gibbs["Fe"]


def deltaGibbs_FeOFeO1p5(Gibbs: pd.Series):
    """
    Change in Gibbs free energy of the reaction

    FeO + 1/4 O2 = FeO1.5
    """

    return Gibbs["FeO1.5"] - Gibbs["FeO"] - (Gibbs["O2"] / 4)


def calculate_XFeO1p5(T_K, deltaGibbs_Fe_wustite, q00, q10, maxiter=300):
    """
    Calculate equilibrium mol fraction FeO1.5 in the reaction

    FeO1.5 + 1/2 Fe = 1.5 FeO

    translated from Hirschmann et al. (2021) matlab script
    """

    func = lambda X_FeO1p5: _Fe_wustite_equilibrium(
        T_K=T_K,
        X_FeO1p5=X_FeO1p5,
        deltaGibbs_Fe_wustite=deltaGibbs_Fe_wustite,
        q00=q00,
        q10=q10,
    )  # Why does Hirschmann square the output? To avoid nagatives?
    try:
        # find solution bracketed by ~0.0 and 1.0
        X_FeO1p5_solution = brentq(
            func, a=1e-6, b=1.0 - 1e-6, xtol=1e-9, maxiter=maxiter
        )
    except RuntimeError as e:
        print(f"{e}\n\nXFeO1.5 did not converge, value set to 1e-6")
        X_FeO1p5_solution = 1e-6
    except ValueError as e:
        print(f"{e}\n\nXFeO1.5 solution outside 0.0-1.0, value set to 1e-6")
        X_FeO1p5_solution = 1e-6

    return X_FeO1p5_solution


def _Fe_wustite_equilibrium(T_K, X_FeO1p5, deltaGibbs_Fe_wustite, q00, q10):
    """
    equilibrium of the reaction

    FeO1.5 + 1/2 Fe = 1.5 FeO (eq. 14)

    as

    deltaG + RT*LN(a_FeO**1.5/aFeO1.5**1) = 0

    including solid solution mixing terms rewritten as

    deltaG + RT*LN(X_FeO**1.5/X_FeO1.5**1) + 1.5 * RT*LN*gammaFeO - RT*LN*gammaFeO1.5

    with:

    RT*LN*gammaFeO          = ((q00+2*q10*(1-X_FeO1.5))*X_FeO1.5**2),
    RT*LN*gammaFeO1.5       = ((1-X_FeO1.5)**2*(q00+q10-2*q10*X_FeO1.5)),
    q00, q10 = solid solution mixing parameters
    a_FeO = gamma * X_FeO = activity of FeO,
    X_FeO1.5 = mol fraction FeO1.5
    X_FeO = 1 - X_FeO1.5

    translated from Hirschmann et al. (2021) matlab script
    """

    RT_LN_gammaFeO = _gamma_FeO(X_FeO1p5=X_FeO1p5, q00=q00, q10=q10)
    RT_LN_gammaFeO1p5 = _gamma_FeO1p5(X_FeO1p5=X_FeO1p5, q00=q00, q10=q10)

    # Half stoichiometry according to Hirschmann et al. (2018) Matlab code
    part_1 = 0.5 * deltaGibbs_Fe_wustite + R * T_K * np.log(
        (1 - X_FeO1p5) ** 1.5 / X_FeO1p5
    )
    part_2 = 1.5 * RT_LN_gammaFeO - RT_LN_gammaFeO1p5

    return part_1 + part_2


def _gamma_FeO(X_FeO1p5, q00, q10):
    """
    equation S5a from Hirschmann et al. (2021)

    RT*LN(a_FeO) =      ((q00+2*q10*(1-X_FeO1.5))*X_FeO1.5**2)

    translated from Hirschmann et al. (2021) matlab script

    Parameters
    ----------
    X_FeO1p5 :  float, array-like
        mol fraction FeO1.5
    q00, q10 : float
        mixing parameters
    """
    return (q00 + 2 * q10 * (1 - X_FeO1p5)) * (X_FeO1p5**2)


def _gamma_FeO1p5(X_FeO1p5, q00, q10):
    """
    equation S5b from Hirschmann et al. (2021)

    RT*LN(a_FeO1.5) =   ((1-X_FeO1.5)**2*(q00+q10-2*q10*X_FeO1.5))

    translated from Hirschmann et al. (2021) matlab script

    Parameters
    ----------
    X_FeO1p5 :  float, array-like
        mol fraction FeO1.5
    q00, q10 : float
        mixing parameters

    """
    return ((1 - X_FeO1p5) ** 2) * (q00 + q10 - 2 * q10 * X_FeO1p5)


@np.vectorize(excluded=["full_output", "suppress_Fe_liquid"])
def muO2_IW(T_K, P_bar, full_output=False, suppress_Fe_liquid=False):
    """
    calculate chemical potential of oxygen at IW and pressure P with equations of state

    translated from Hirschmann et al. (2021) matlab script
    """

    q00 = mixing_parameters["q00"]
    q10 = mixing_parameters["q10"]

    Gibbs, Fe_phase = Gibbs_IW(
        P_bar=P_bar, T_K=T_K, Fe_phase=True, suppress_Fe_liquid=suppress_Fe_liquid
    )
    deltaG_FeO_FeO1p5 = deltaGibbs_FeOFeO1p5(Gibbs=Gibbs)
    deltaG_Fe_wustite = deltaGibbs_Fe_wustite(Gibbs=Gibbs)

    X_FeO1p5 = calculate_XFeO1p5(
        T_K=T_K,
        deltaGibbs_Fe_wustite=deltaG_Fe_wustite,
        q00=q00,
        q10=q10,
    )

    X_FeO = 1 - X_FeO1p5
    y = 1 - 2 / (X_FeO1p5 + 2)  # Fe(1-y)O

    gamma_1 = _gamma_FeO(X_FeO1p5=X_FeO1p5, q00=q00, q10=q10)
    gamma_2 = _gamma_FeO1p5(X_FeO1p5=X_FeO1p5, q00=q00, q10=q10)

    mu_O2 = 4 * (
        deltaG_FeO_FeO1p5 + R * T_K * np.log(X_FeO1p5 / X_FeO) + gamma_2 - gamma_1
    )  # chemical potential

    output = mu_O2, Fe_phase, y if full_output else mu_O2

    return output


def fO2_IW(logshift: float, T_K, P_bar, full_output=False, suppress_Fe_liquid=False):
    """
    Calculate oxygen fugacity at the Iron-Wustite buffer according to Hirschmann (2021)\ [18]_

    Parameters
    ----------
    logshift    : int, float
        log units shift of fO2
    T_K         : float, array-like
        temperature in Kelvin
    P_bar       : float, aray-like
        pressure in bar
    full_output : boolean
        ouputs fO2 only if False. fO2, stable Fe phase and Fe(1-y)O if True
    suppress_Fe_liquid  : boolean
        ignore Fe liquid if True


    Returns
    -------
    fO2 : float
        oxygen fugacity
    """
    try:
        logshift = logshift
    except TypeError:
        pass

    offset = 10**logshift

    T_K, P_bar = repeat_vars(T_K, P_bar)

    mu_O2, Fe_phase, y = muO2_IW(
        T_K=T_K, P_bar=P_bar, full_output=True, suppress_Fe_liquid=suppress_Fe_liquid
    )

    fO2 = np.exp(mu_O2 / (R * T_K)) * offset
    try:
        fO2 = np.float32(fO2.item())
    except ValueError:
        fO2 = fO2.astype(np.float32)

    idx = T_K.index if isinstance(T_K, (pd.Series)) else range(len(mu_O2))
    output = (
        (pd.DataFrame({"fO2": fO2, "Fe_phase": Fe_phase, "Fe(1-y)O": y}, index=idx))
        if full_output
        else fO2
    )

    return output


def _fO2_IW_Campbell(logshift, T_K, P_bar):
    """
    from:
    Campbell et al. (2009) High pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers. EPSL

    Supplementary table S4

    calibrated between 4-57.6 GPa and 873-2096 K
    """
    P_GPa = P_bar * 1e5 / 1e9
    offset = 10**logshift

    # TODO Armstrong et al. (2018) Supplementary Materials eq. S15 has a positive P**2 coefficient, check which one is correct.
    part_1 = 6.54106 + 1.23e-3 * P_GPa
    part_2 = (-28164 + 546.32 * P_GPa - 1.1341 * P_GPa**2 + 1.93e-3 * P_GPa**3) / T_K

    log10fO2 = part_1 + part_2

    return 10**log10fO2 * offset


def _fO2_FeFeO94_Oneill_Huebner(logshift, T_K, P_bar):
    """
    Fe-FeO(0.94) buffer with the 100 kPa relation from O'neill (1988) and pressure term derived from Huebner (1971).
    Formulation is from Zhang et al. (2017) Supplementary Materials eq. S1
    """

    P_GPa = P_bar * 1e5 / 1e9
    offset = 10**logshift

    part_1 = -28777.89 / T_K + 14.0572
    part_2 = -2.039 * np.log10(T_K) + 550 * (P_GPa - 1e-4) / T_K

    log10fO2 = part_1 + part_2

    return 10 ** (log10fO2) * offset


def _fO2_IW_Zhang(logshift, T_K, P_bar):
    """
    Following Zhang et al. (2017):

    Low pressure part (< 5 GPa) is calculated as an interpolation between Oneill+Huebner and Campbell, while the high pressure part (>5 GPa) is pure Campbell
    """

    P_bar, T_K = make_iterable(P_bar, T_K)

    array_length_one = np.array([len(T_K), len(P_bar)]) == 1
    if array_length_one.any():
        array_length = np.array([len(T_K), len(P_bar)]).max()
        P_bar, T_K = [
            np.ones(array_length) * i[0] if len(i) == 1 else i for i in (P_bar, T_K)
        ]

    low_pressure = P_bar < 5e4
    high_pressure = P_bar > 5e4

    ONeill_Huebner_low_pressure = np.log10(
        _fO2_FeFeO94_Oneill_Huebner(
            logshift=logshift, T_K=T_K[low_pressure], P_bar=P_bar[low_pressure]
        )
    )
    Campbell_low_pressure = np.log10(
        _fO2_IW_Campbell(
            logshift=logshift, T_K=T_K[low_pressure], P_bar=P_bar[low_pressure]
        )
    )

    fO2_low_pressure = 10 ** (
        ONeill_Huebner_low_pressure * (1 - 0.2 * P_bar[low_pressure] / 1e4)
        + (0.2 * P_bar[low_pressure] / 1e4) * Campbell_low_pressure
    )

    fO2_high_pressure = _fO2_IW_Campbell(
        logshift=logshift, T_K=T_K[high_pressure], P_bar=P_bar[high_pressure]
    )

    fO2 = np.concatenate([fO2_low_pressure, fO2_high_pressure])

    return fO2
