from .volatile_solubility_models import Allison2022
from .volatile_solubility_models import IaconoMarziano
from .volatile_solubility_models import Shiskina
from ..configuration import configuration


def calculate_saturation(oxide_wtPercents, **kwargs):
    """
    Docstring
    """

    model = kwargs.get("model", configuration.volatile_solubility)
    species = kwargs.get("species", configuration.volatile_species)

    equation = getattr(globals()[model], species).calculate_saturation

    return equation(oxide_wtPercents, **kwargs)


def calculate_solubility(oxide_wtPercents, P_bar, **kwargs):
    """
    Docstring
    """

    model = configuration.volatile_solubility
    species = configuration.volatile_species
    equation = getattr(globals()[model], species).calculate_solubility

    return equation(oxide_wtPercents, P_bar, **kwargs)
