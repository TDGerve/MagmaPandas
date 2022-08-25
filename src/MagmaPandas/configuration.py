from .parse.validate import _check_setter


Fe3Fe2_models = ["borisov", "kressCarmichael"]
Kd_ol_FeMg_models = ["toplis", "blundy"]
melt_thermometers = ["putirka2008_14", "putirka2008_15", "putirka2008_16"]
volatile_solubility_models = ["IaconoMarziano", "Allison2022", "Shiskina"]
volatile_species_options = ["co2", "h2o", "mixed"]


class _meta_configuration(type):
    def __init__(cls, *args, **kwargs):
        cls._dQFM = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "borisov"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"
        cls._volatile_species = "co2"

    @property
    def dQFM(cls):
        return cls._dQFM

    @dQFM.setter
    def dQFM(cls, value):
        cls._dQFM = value

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


class configuration(metaclass=_meta_configuration):
    @classmethod
    def reset(cls):
        cls._dQFM = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "borisov"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"
        cls._volatile_species = "co2"

    @classmethod
    def print(cls):
        """
        Docstring
        """

        variables = {
            "\u0394QFM": "_dQFM",
            "Melt Fe3+/Fe2+": "_Fe3Fe2_model",
            "Kd Fe-Mg ol-melt": "_Kd_model",
            "Melt thermometer": "_melt_thermometer",
            "Volatile solubility model": "_volatile_solubility",
            "Volatile species": "_volatile_species",
        }

        names_length = max([len(i) for i in variables.keys()]) + 5
        pad_right = 15
        pad_total = names_length + pad_right

        print(" MagmaPandas ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("\nGeneral settings".ljust(pad_total, "_"))
        print(f"{'fO2 buffer':.<{names_length}}{'QFM':.>{pad_right}}")
        for param, value in variables.items():
            # model_attr = f"_configuration{model}"
            print(f"{param:.<{names_length}}{getattr(cls, value):.>{pad_right}}")
        print("\n")
