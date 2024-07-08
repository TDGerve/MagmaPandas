from typing import Callable

import pandas as pd

from MagmaPandas.parse_io import convert_to_series, match_indeces


def equilibrium_forsterite(melt_mol_fractions, Kd, Fe3Fe2):
    """
    Parameters
    ----------
    melt_mol_fractions  :
        melt composition on mol fractions.
    Kd  :
        (Fe2+/Mg)olivine / (Fe2+/Mg)liquid partitioning coefficient
    Fe3Fe2  :
        melt Fe3+/Fe2+ ratio

    Returns
    -------
    Equilibrium forsterite fraction as Mg/(Mg + Fe)
    """

    # Mg_no = 1 + (1 + Fe2Mg_liquid)
    Fe2FeTotal = 1 / (1 + Fe3Fe2)
    Fe2Mg_liquid = melt_mol_fractions["FeO"] * Fe2FeTotal / melt_mol_fractions["MgO"]

    return 1 / (1 + Kd * Fe2Mg_liquid)


def iterate_Kd_scalar(
    melt_mol_fractions: pd.Series,
    forsterite_initial,
    T_K,
    Fe3Fe2,
    Kd_model: Callable,
    *args,
    **kwargs,
):
    """

    Parameters
    ----------
    melt_mol_fractions : pd.DataFrame
        melt composition in oxide mol fraction.
    olivine_forsterite
        forsterite fraction in olivine as Mg * 100 / (Mg + Fe)
    T_K :
        Temperature in Kelvin
    P_bar :
        Pressure in bar
    Fe3Fe2 :
        melt Fe3+/Fe2+ ratio
    """
    # model_name = kwargs.pop("model", configuration.Kd_model)
    # Kd_model = Kd_models[model_name]

    fo_converge_default = 0.001
    fo_converge = kwargs.pop("fo_converge", fo_converge_default)

    for var, name in zip(
        (forsterite_initial, T_K, Fe3Fe2), ("forsterite_initial", "T_K", "Fe3Fe2")
    ):
        if not isinstance(var, (float, int)):
            raise ValueError(f"{name} should be float or int")

    # initialise Kd
    Kd = Kd_model(
        melt_mol_fractions=melt_mol_fractions,
        forsterite=forsterite_initial,
        T_K=T_K,
        **kwargs,
    )

    forsterite_EQ = equilibrium_forsterite(
        melt_mol_fractions=melt_mol_fractions, Kd=Kd, Fe3Fe2=Fe3Fe2
    )

    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = abs(forsterite_initial - forsterite_EQ) / forsterite_initial

    iterate = forsterite_delta > fo_converge
    # iterate until equilibrium forsterite content doesn't change any more
    while iterate:

        Kd = Kd_model(
            melt_mol_fractions=melt_mol_fractions,
            forsterite=forsterite_initial,
            T_K=T_K,
            **kwargs,
        )

        forsterite_initial = forsterite_EQ.copy()
        forsterite_EQ = equilibrium_forsterite(
            melt_mol_fractions=melt_mol_fractions, Kd=Kd, Fe3Fe2=Fe3Fe2
        )
        forsterite_delta = abs(forsterite_initial - forsterite_EQ) / forsterite_initial
        iterate = forsterite_delta > fo_converge

    return Kd


def iterate_Kd_vectorized(
    melt_mol_fractions: pd.DataFrame,
    forsterite_initial,
    T_K,
    Fe3Fe2,
    Kd_model: Callable,
    **kwargs,
):
    """

    Parameters
    ----------
    melt_mol_fractions : pd.DataFrame
        melt composition in oxide mol fraction.
    olivine_forsterite
        forsterite fraction in olivine as Mg * 100 / (Mg + Fe)
    T_K :
        Temperature in Kelvin
    P_bar :
        Pressure in bar
    Fe3Fe2 :
        melt Fe3+/Fe2+ ratio
    """

    # model_name = kwargs.pop("model", configuration.Kd_model)
    # Kd_model = Kd_models[model_name]

    match_indeces(melt_mol_fractions, [T_K, Fe3Fe2, *kwargs.values()])

    melt_mol_fractions = melt_mol_fractions.copy()

    forsterite, T_K, Fe3Fe2, *args = convert_to_series(
        [forsterite_initial, T_K, Fe3Fe2, *kwargs.values()], melt_mol_fractions.index
    )
    if len(kwargs.keys()) < 2:
        args = [args]

    kwargs = {key: value for key, value in zip(kwargs.keys(), args)}

    fo_converge_default = 0.001
    fo_converge = kwargs.pop("fo_converge", fo_converge_default)

    # initialise Kds
    Kd = Kd_model(
        melt_mol_fractions=melt_mol_fractions, forsterite=forsterite, T_K=T_K, **kwargs
    )

    forsterite_EQ = equilibrium_forsterite(
        melt_mol_fractions=melt_mol_fractions, Kd=Kd, Fe3Fe2=Fe3Fe2
    )
    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

    iterate = forsterite_delta > fo_converge
    # iterate until equilibrium forsterite content doesn't change any more
    while sum(iterate) > 0:

        Kd = Kd_model(
            melt_mol_fractions=melt_mol_fractions,
            forsterite=forsterite_EQ,
            T_K=T_K,
            **kwargs,
        )
        forsterite[iterate] = forsterite_EQ[iterate].copy()
        forsterite_EQ = equilibrium_forsterite(
            melt_mol_fractions=melt_mol_fractions, Kd=Kd, Fe3Fe2=Fe3Fe2
        )

        forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite
        iterate = forsterite_delta > fo_converge

    return Kd
