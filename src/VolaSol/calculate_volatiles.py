import VolaSol.solubility_models as vsm

from MagmaPandas.configuration import configuration

def get_solubility_model(model_name: str=None, species: str=None) -> vsm.Solubility_model:

    if model_name is None:
        model_name = configuration.volatile_solubility
    model = getattr(vsm, model_name)
    if species is None:
        species = configuration.volatile_species

    return getattr(model, species)


def calculate_saturation(oxide_wtPercents, *args, solubility_model: vsm.Solubility_model=None, **kwargs):
    """
    Docstring
    """
    if solubility_model is None:
        solubility_model = get_solubility_model()

    equation = solubility_model.calculate_saturation

    return equation(oxide_wtPercents, *args, **kwargs)


def calculate_solubility(oxide_wtPercents, P_bar, *args, solubility_model: vsm.Solubility_model=None, **kwargs):
    """
    Docstring
    """
    if solubility_model is None:
        solubility_model = get_solubility_model()

    equation = solubility_model.calculate_solubility

    return equation(oxide_wtPercents, P_bar, *args, **kwargs)



