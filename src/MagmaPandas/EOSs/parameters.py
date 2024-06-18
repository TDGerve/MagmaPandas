class EOSparams:
    """
    h                       enthalpy of formation
    s                       entropy
    v0                      volume at 1bar, 298K
    n                       atoms per formula unit
    a0                      coefficient of thermal expansion
    K0                      bulk modulus at 1bar, 298K
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
