from .parse.validate import _check_setter
from .parse.validate import _check_argument


"""
ALL CONFIGURATION OBJECTS NEED TO BE REWRITTEN BECAUSE CLASS PROPERTIES DON'T WORK LIKE THIS.
IT EITHER NEEDS:
    - A METACLASS WHERE ALL CLASS VARIABLES AND PROPERTIES ARE DEFINED
    - A CUSTOM CLASSPROPERTY DECORATOR
    - ALL CONFIGURATION OBJECTS INITIALISED INSIDE THE MODULE __INIT__ WITH ALL ATTRIBUTES AS INSTANCE ATTRIBUTES
"""

Fe3Fe2_models = ["borisov", "kressCarmichael"]
Kd_ol_FeMg_models = ["toplis", "blundy"]
melt_thermometers = ["putirka2008_14", "putirka2008_15", "putirka2008_16"]
volatile_solubility_models = ["IaconoMarziano", "Allison2022"]


class _meta_configuration(type):

    def __init__(cls, *args, **kwargs):
        cls._dQFM = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "borisov"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"

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



class configuration(metaclass=_meta_configuration):

    @classmethod
    def reset(cls):
        cls._dQFM = 1
        cls._Kd_model = "toplis"
        cls._Fe3Fe2_model = "borisov"
        cls._melt_thermometer = "putirka2008_15"
        cls._volatile_solubility = "IaconoMarziano"

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
        }

        names_length = max([len(i) for i in variables.keys()]) + 5
        pad_right = 15
        pad_total = names_length + pad_right

        print(" MagmaPandas ".center(pad_total, "#"))
        print("".ljust(pad_total, "#"))
        print("\nGeneral settings".ljust(pad_total, "_"))
        print(f"{'fO2 buffer':.<{names_length}}{'QFM':.>{pad_right}}")
        for param, model in variables.items():
            # model_attr = f"_configuration{model}"
            print(f"{param:.<{names_length}}{getattr(cls, model):.>{pad_right}}")



    


# class configuration:

#     __QFMlogshift = 1
#     __Fe3Fe2_model = "borisov"
#     __Kd_model = "toplis"
#     __melt_thermometer = "putirka2008_14"
#     __volatile_solubility = "IaconoMarziano"

#     @property
#     def QFMlogshift(cls):
#         return cls.__QFMlogshift

#     @QFMlogshift.setter
#     def QFMlogshift(cls, value):
#         cls.__QFMlogshift = value

#     @property
#     def Kd_model(cls):
#         return cls.__Kd_model

#     @Kd_model.setter
#     @_check_setter(Kd_ol_FeMg_models)
#     def Kd_model(cls, model: str):
#         cls.__Kd_model = model

#     @property
#     def Fe3Fe2_model(cls):
#         return cls.__Fe3Fe2_model

#     @Fe3Fe2_model.setter
#     @_check_setter(Fe3Fe2_models)
#     def Fe3Fe2_model(cls, model: str):
#         cls.__Fe3Fe2_model = model

#     @property
#     def melt_thermometer(cls):
#         return cls.__melt_thermometer

#     @melt_thermometer.setter
#     @_check_setter(melt_thermometers)
#     def melt_thermometer(cls, model: str):
#         cls.__melt_thermometer = model

#     @property
#     def volatile_solubility(cls):
#         return cls.__volatile_solubility

#     @volatile_solubility.setter
#     @_check_setter(volatile_solubility_models)
#     def volatile_solubility(cls, model: str):
#         cls.__volatile_solubility = model

#     @classmethod
#     def reset(cls):
#         cls.__Kd_model = "toplis"
#         cls.__Fe3Fe2_model = "borisov"
#         cls.__melt_thermometer = "putirka2008_14"
#         cls.__volatile_solubility = "IaconoMarziano"

#     @classmethod
#     def print(cls):
#         """ 
#         Docstring
#         """

#         variables = {
#             "\u0394QFM": "__QFMlogshift",
#             "Melt Fe3+/Fe2+": "__Fe3Fe2_model",
#             "Kd Fe-Mg ol-melt": "__Kd_model",
#             "Melt thermometer": "__melt_thermometer",
#             "Volatile solubility model": "__volatile_solubility",
#         }

#         names_length = max([len(i) for i in variables.keys()]) + 5
#         pad_right = 15
#         pad_total = names_length + pad_right

#         print(" MagmaPandas ".center(pad_total, "#"))
#         print("".ljust(pad_total, "#"))
#         print("\nGeneral settings".ljust(pad_total, "_"))
#         print(f"{'fO2 buffer':.<{names_length}}{'QFM':.>{pad_right}}")
#         for param, model in variables.items():
#             model_attr = f"_configuration{model}"
#             print(f"{param:.<{names_length}}{getattr(cls, model_attr):.>{pad_right}}")
