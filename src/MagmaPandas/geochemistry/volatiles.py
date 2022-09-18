from .volatile_solubility_models import Allison2022
from .volatile_solubility_models import IaconoMarziano
from .volatile_solubility_models import Shiskina
from .volatile_solubility_models.solubility_abc import Solubility_model
from ..configuration import configuration



def _get_solubility_model():

    model = configuration.volatile_solubility
    species = configuration.volatile_species

    return getattr(globals()[model], species)


def calculate_saturation(oxide_wtPercents, solubility_model=_get_solubility_model(), **kwargs):
    """
    Docstring
    """

    model = configuration.volatile_solubility
    species = configuration.volatile_species
    equation = solubility_model.calculate_saturation


    return equation(oxide_wtPercents, **kwargs)


def calculate_solubility(oxide_wtPercents, P_bar, solubility_model=_get_solubility_model(), **kwargs):
    """
    Docstring
    """

    model = configuration.volatile_solubility
    species = configuration.volatile_species
    equation = solubility_model.calculate_solubility

    return equation(oxide_wtPercents, P_bar, **kwargs)

class Volatiles:

    def __init__(self, solubility_model: Solubility_model=_get_solubility_model()):
        self.model = solubility_model

    def calculate_saturation(self, oxide_wtPercents, **kwargs):
        return self.model.calculate_saturation(oxide_wtPercents, **kwargs)

    def calculate_solubility(self, oxide_wtPercents, P_bar, **kwargs):
        return self.model.calculate_solubility(oxide_wtPercents, P_bar, **kwargs)

