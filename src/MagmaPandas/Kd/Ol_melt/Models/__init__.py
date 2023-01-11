from .Blundy import FeMg_blundy
from .Toplis import FeMg_Toplis

Kd_models = {"toplis": FeMg_Toplis, "blundy": FeMg_blundy}


def equilibrium_forsterite(Kd, Fe2Mg_liquid):
    """
    Parameters
    ----------
    Kd :
        (Fe_olivine / Fe_liquid) * (Mg_liquid / Mg_olivine) partitioning coefficient
    Fe2Mg
        melt Fe2+/Mg ratio

    Returns
    -------
    Equilibrium forsterite fraction as Mg/(Mg + Fe)
    """

    return 1 / (1 + Kd * Fe2Mg_liquid)
