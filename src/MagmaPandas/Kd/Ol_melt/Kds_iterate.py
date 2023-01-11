import pandas as pd

from MagmaPandas.configuration import configuration

from .Models import Kd_models, equilibrium_forsterite


def calculate_Kd(
    Melt_mol_fractions: pd.Series,
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

    fo_converge_default = 0.001
    fo_converge = kwargs.pop("fo_converge", fo_converge_default)

    for var, name in zip(
        (forsterite_initial, T_K, Fe3Fe2), ("forsterite_initial", "T_K", "Fe3Fe2")
    ):
        if not isinstance(var, (float, int)):
            raise ValueError(f"{name} should be float or int")

    # initialise Kd
    Kd = Kd_model.calculate_Kd(Melt_mol_fractions, forsterite_initial, T_K, **kwargs)

    Fe2Fe_total = 1 / (1 + Fe3Fe2)
    Fe2Mg = (Melt_mol_fractions["FeO"] * Fe2Fe_total) / Melt_mol_fractions["MgO"]
    forsterite_EQ = equilibrium_forsterite(Kd, Fe2Mg)

    # Difference between observed Fo and equilibrium Fo
    forsterite_delta = abs(forsterite_initial - forsterite_EQ) / forsterite_initial

    iterate = forsterite_delta > fo_converge
    # iterate until equilibrium forsterite content doesn't change any more
    while iterate:

        Kd = Kd_model.calculate_Kd(Melt_mol_fractions, forsterite_EQ, T_K, **kwargs)

        forsterite_initial = forsterite_EQ.copy()
        forsterite_EQ = equilibrium_forsterite(Kd, Fe2Mg)
        forsterite_delta = abs(forsterite_initial - forsterite_EQ) / forsterite_initial
        iterate = forsterite_delta > fo_converge

    return Kd
