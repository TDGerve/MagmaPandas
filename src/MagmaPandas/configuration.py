"""
Global configuration of MagmaPandas settings.
"""

from MagmaPandas.parse_io.validate import _check_setter

Fe3Fe2_models = [
    "borisov",
    "kressCarmichael",
    "jayasuriya",
    "putirka2016_6b",
    "putirka2016_6c",
    "Deng2020",
    "Oneill2006",
    "Oneill2018",
    "Armstrong2019",
    "Zhang2017",
    "Hirschmann2022",
    "Sun2024",
]
Kd_ol_FeMg_models = [
    "toplis",
    "blundy",
    "putirka2016_8a",
    "putirka2016_8b",
    "putirka2016_8c",
    "putirka2016_8d",
    "sun2020",
]
fO2_buffers = ["QFM", "IW"]
melt_thermometers = ["putirka2008_14", "putirka2008_15", "putirka2008_16", "sun2020"]
volatile_solubility_models = ["IaconoMarziano", "Allison2022", "Shiskina"]
volatile_species_options = ["co2", "h2o", "mixed"]


_variables = {
    "fO2 buffers": fO2_buffers,
    "Melt Fe3+/Fe2+ models": Fe3Fe2_models,
    "Ol-melt Fe-Mg Kd models": Kd_ol_FeMg_models,
    "Melt thermometers": melt_thermometers,
    "Volatile solubility models": volatile_solubility_models,
    "Volatile species": volatile_species_options,
}

_names_length = max([len(i) for i in _variables.keys()]) + 5
_pad_right = max([len(", ".join(i)) for i in _variables.values()])
# _pad_total = _names_length + _pad_right

[len(", ".join(i)) for i in _variables.values()]
_pad_total = 70
_pad_right = _pad_total - _names_length
_new_line = "\n"

configuration_options = (
    f"{_new_line}{' MagmaPandas ':#^{_pad_total}}"
    f"{_new_line}{'':#^{_pad_total}}"
    f"{_new_line}{'Configuration options':_<{_pad_total}}"
)

for param, values in _variables.items():

    val_string_2 = ""
    idxs = [0, 0]
    start_idx = 0
    while idxs[-1] < len(values):

        for i, _ in enumerate(values):
            if len(", ".join(values[start_idx : i + 1])) > (_pad_total - _names_length):
                idxs.insert(-1, i)
                start_idx = i
                continue
            idxs[-1] = i + 1

    if len(idxs) > 2:
        val_string_1 = ", ".join(values[: idxs[1]])
        val_string_2 = f""
        for i, idx in enumerate(idxs[2:]):
            val_string_2 += (
                f"{_new_line}{', '.join(values[idxs[i+1]:idx]):>{_pad_total}}"
            )
        val_string_2 += f"{_new_line}"
    else:
        val_string_1 = ", ".join(values)

    line = f"{_new_line}{param:.<{_names_length}}{val_string_1:.>{_pad_right}}{val_string_2}"
    configuration_options += line


class _meta_configuration(type):
    def __init__(cls, *args, **kwargs):
        cls._fO2buffer = "QFM"
        cls._dfO2 = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "Sun2024"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"
        cls._volatile_species = "mixed"

    @property
    def fO2buffer(cls):
        return cls._fO2buffer

    @fO2buffer.setter
    @_check_setter(fO2_buffers)
    def fO2buffer(cls, value: str):
        cls._fO2buffer = value

    @property
    def dfO2(cls):
        return cls._dfO2

    @dfO2.setter
    def dfO2(cls, value: float | int):
        cls._dfO2 = value

    @property
    def Kd_model(cls):
        return cls._Kd_model

    @Kd_model.setter
    @_check_setter(Kd_ol_FeMg_models)
    def Kd_model(cls, model: str):
        cls._Kd_model = model

    @property
    def Fe3Fe2_model(cls):
        return cls._Fe3Fe2_model

    @Fe3Fe2_model.setter
    @_check_setter(Fe3Fe2_models)
    def Fe3Fe2_model(cls, model: str):
        cls._Fe3Fe2_model = model

    @property
    def melt_thermometer(cls):
        return cls._melt_thermometer

    @melt_thermometer.setter
    @_check_setter(melt_thermometers)
    def melt_thermometer(cls, model: str):
        cls._melt_thermometer = model

    @property
    def volatile_solubility(cls):
        return cls._volatile_solubility

    @volatile_solubility.setter
    @_check_setter(volatile_solubility_models)
    def volatile_solubility(cls, model: str):
        cls._volatile_solubility = model

    @property
    def volatile_species(cls):
        return cls._volatile_species

    @volatile_species.setter
    @_check_setter(volatile_species_options)
    def volatile_species(cls, model: str):
        cls._volatile_species = model

    def __str__(cls):

        variables = {
            "fO2 buffer": "_fO2buffer",
            "\u0394 fO2": "_dfO2",
            "Melt Fe3+/Fe2+": "_Fe3Fe2_model",
            "Kd Fe-Mg ol-melt": "_Kd_model",
            "Melt thermometer": "_melt_thermometer",
            "Volatile solubility model": "_volatile_solubility",
            "Volatile species": "_volatile_species",
        }

        names_length = max([len(i) for i in variables.keys()]) + 5
        pad_right = 15
        pad_total = names_length + pad_right
        new_line = "\n"

        message = (
            f"{new_line}{' MagmaPandas ':#^{pad_total}}"
            f"{new_line}{'':#^{pad_total}}"
            f"{new_line}{'General settings':_<{pad_total}}"
            # f"{new_line}{'fO2 buffer':.<{names_length}}{'QFM':.>{pad_right}}"
        )

        parameter_settings = ""
        for param, value in variables.items():
            parameter_settings += (
                f"{new_line}{param:.<{names_length}}{getattr(cls, value):.>{pad_right}}"
            )

        return message + parameter_settings


class configuration(metaclass=_meta_configuration):
    """
    Class for configuring global settings in MagmaPandas.

    Attributes
    ----------
    dQFM    : int, float
        Log units shift of the QFM |fO2| buffer. Default value: 1
    Kd_model : str
        Olivine-melt Fe-Mg partitioning model. Available models: 'toplis'\ [10]_ and 'blundy'\ [11]_. Default value: 'toplis'
    Fe3Fe2_model : str
        Melt |Fe3Fe2| model. Available models: 'borisov'\ [1]_ and 'kressCarmichael'\ [2]_. Default value: 'borisov'
    melt_thermometer : str
        Melt-only thermometer. Available models: 'putirka2008_14'\ [15]_, 'putirka2008_15'\ [15]_ and 'putirka2008_16'\ [15]_. Default value: 'putirka2008_15'
    volatily_solubility : str
        |CO2|-|H2O| solubility model. Available models: 'IaconoMarziano'\ [18]_, 'Allison2022'\ [17]_, 'Shiskina'\ [19]_. Default value: 'IaconoMarziano'
    volatile_species : str
        Fluid phase species. Options: 'h2o', 'co2' or 'mixed'. Default value: 'mixed'
    """

    @classmethod
    def reset(cls):
        """
        Reset to default values
        """
        cls._fO2buffer = "QFM"
        cls._dfO2 = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "Sun2024"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"
        cls._volatile_species = "mixed"
