from .parse.validate import _check_setter



Fe3Fe2_models = ["borisov", "kressCarmichael"]
Kd_ol_FeMg_models = ["toplis", "blundy"]


class configuration:

    __Kd_model = "toplis"
    __Fe3Fe2_model = "borisov"

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

    @staticmethod
    def reset():
        configuration.__Kd_model = "toplis"
        configuration.__Fe3Fe2_model = "borisov"

    @staticmethod
    def print():
        """ """
        parameters = ["Kd Fe-Mg ol-melt", "melt Fe3+/Fe2+"]
        models = [f"_configuration{i}" for i in ["__Kd_model", "__Fe3Fe2_model"]]
        for param, model in zip(parameters, models):
            print(f"{param:.<20}{getattr(configuration, model):.>15}")
