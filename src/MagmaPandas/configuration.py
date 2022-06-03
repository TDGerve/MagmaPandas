from .parse.validate import _check_setter



Fe3Fe2_models = ["borisov", "kressCarmichael"]
Kd_ol_FeMg_models = ["toplis", "blundy"]
melt_thermometers = ["putirka2008_14", "putirka2008_16"]


class configuration:

    __Kd_model = "toplis"
    __Fe3Fe2_model = "borisov"
    __melt_thermometer = "putirka2008_14"

    @property
    def Kd_model(cls):
        return configuration.__Kd_model

    @Kd_model.setter
    @_check_setter(Kd_ol_FeMg_models)
    def Kd_model(cls, model: str):
        configuration.__Kd_model = model

    @property
    def Fe3Fe2_model(cls):
        return configuration.__Fe3Fe2_model

    @Fe3Fe2_model.setter
    @_check_setter(Fe3Fe2_models)
    def Fe3Fe2_model(cls, model: str):
        configuration.__Fe3Fe2_model = model
    
    @property
    def melt_thermometer(cls):
        return configuration.__melt_thermometer
        
    @melt_thermometer.setter
    @_check_setter(melt_thermometers)
    def melt_thermometer(cls, model: str):
        configuration.__melt_thermometer = model


    @staticmethod
    def reset():
        configuration.__Kd_model = "toplis"
        configuration.__Fe3Fe2_model = "borisov"
        configuration.__melt_thermometer = "putirka2008_14"

    @staticmethod
    def print():
        """ """
        parameters = ["Kd Fe-Mg ol-melt", "melt Fe3+/Fe2+", "melt thermometer"]
        models = [f"_configuration{i}" for i in ["__Kd_model", "__Fe3Fe2_model", "__melt_thermometer"]]
        for param, model in zip(parameters, models):
            print(f"{param:.<20}{getattr(configuration, model):.>15}")
