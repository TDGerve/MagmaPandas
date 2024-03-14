import pandas as pd

from MagmaPandas.configuration import configuration
from MagmaPandas.Kd.Ol_melt.Models import Kd_models, equilibrium_forsterite
from MagmaPandas.parse_io import convert_to_series, match_indeces


def calculate_Kd(
    Melt_mol_fractions: pd.DataFrame,
    forsterite_initial,
    T_K,
    Fe3Fe2,
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

    model_name = kwargs.pop("model", configuration.Kd_model)
    Kd_model = Kd_models[model_name]

    match_indeces(Melt_mol_fractions, [T_K, Fe3Fe2, *kwargs.values()])

    melt_mol_fractions = Melt_mol_fractions.copy()

    forsterite, T_K, Fe3Fe2, args = convert_to_series(
        [forsterite_initial, T_K, Fe3Fe2, *kwargs.values()], melt_mol_fractions.index
    )
    if len(kwargs.keys()) < 2:
        args = [args]

    kwargs = {key: value for key, value in zip(kwargs.keys(), args)}

    fo_converge_default = 0.001
    fo_converge = kwargs.pop("fo_converge", fo_converge_default)

    # initialise Kds
    Kd = Kd_model.calculate_Kd(melt_mol_fractions, forsterite, T_K, **kwargs)

    Fe2FeTotal = 1 / (1 + Fe3Fe2)
    Fe2Mg = melt_mol_fractions["FeO"] * Fe2FeTotal / melt_mol_fractions["MgO"]

    forsterite_EQ = equilibrium_forsterite(Kd, Fe2Mg)
    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = abs(forsterite - forsterite_EQ) / forsterite

    iterate = forsterite_delta > fo_converge
    # iterate until equilibrium forsterite content doesn't change any more
    while sum(iterate) > 0:

        Kd[iterate] = Kd_model.calculate_Kd(
            melt_mol_fractions.loc[iterate],
            forsterite_EQ.loc[iterate],
            T_K[iterate],
            **{key: value[iterate] for key, value in kwargs.items()},
        )
        forsterite[iterate] = forsterite_EQ[iterate].copy()
        forsterite_EQ[iterate] = forsterite_EQ = equilibrium_forsterite(
            Kd[iterate], Fe2Mg[iterate]
        )

        forsterite_delta[iterate] = (
            abs(forsterite[iterate] - forsterite_EQ[iterate]) / forsterite[iterate]
        )
        iterate = forsterite_delta > fo_converge

    return Kd
