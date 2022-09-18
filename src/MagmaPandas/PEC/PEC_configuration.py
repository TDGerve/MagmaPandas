from ..parse_io.validate import _check_setter, _check_value


Fe2_options = ["buffered", "incompatible"]


class _meta_PEC_configuration(type):
    """
    Metaclass for setting class properties
    """

    def __init__(cls, *args, **kwargs):
        cls._Fe2_behaviour = "buffered"
        # Initial stepsize
        cls._stepsize_equilibration = 0.002
        cls._stepsize_crystallisation = 0.05
        # Reduction factor for stepsize after overstepping
        cls._decrease_factor = 5
        # Convergence values
        cls.FeO_converge = 0.05
        cls.Kd_converge = 1e-3

    @property
    def Fe2_behaviour(cls):
        return cls._Fe2_behaviour

    @Fe2_behaviour.setter
    @_check_setter(Fe2_options)
    def Fe2_behaviour(cls, value):
        cls._Fe2_behaviour = value

    @property
    def stepsize_equilibration(cls):
        return cls._stepsize_equilibration

    @stepsize_equilibration.setter
    @_check_value("value", [0, 1])
    def stepsize_equilibration(cls, value):
        cls._stepsize_equilibration = value

    @property
    def stepsize_crystallisation(cls):
        return cls._stepsize_crystallisation

    @stepsize_crystallisation.setter
    @_check_value("value", [0, 1])
    def stepsize_crystallisation(cls, value):
        cls._stepsize_crystallisation = value

    @property
    def decrease_factor(cls):
        return cls._decrease_factor

    @decrease_factor.setter
    @_check_value("value", [1, 50])
    def decrease_factor(cls, value):
        cls._decrease_factor = value

    def __str__(cls):

        variables = {
            "Fe2+ behaviour": "_Fe2_behaviour",
            "Stepsize equilibration (moles)": "_stepsize_equilibration",
            "Stepsize crystallisation (moles)": "_stepsize_crystallisation",
            "Decrease factor": "_decrease_factor",
            "FeO convergence (wt. %)": "FeO_converge",
            "Kd convergence": "Kd_converge",
        }

        names_length = max([len(i) for i in variables.keys()]) + 5
        values_length = max([len(str(getattr(cls, i))) for i in variables.values()])

        pad_right = 20
        pad_total = names_length + pad_right
        new_line = "\n"

        message = (
            f"{new_line}{' Post-entrapment modification ':#^{pad_total}}"
            f"{new_line}{' correction model ':#^{pad_total}}"
            f"{new_line}{'Settings':_<{pad_total}}"
        )

        parameter_settings = ""
        for param, value in variables.items():
            value_str = f"{getattr(cls, value):<{values_length}}"
            parameter_settings += (
                f"{new_line}{param:.<{names_length}}{value_str:.>{pad_right}}"
            )

        return message + parameter_settings


class PEC_configuration(metaclass=_meta_PEC_configuration):
    @classmethod
    def reset(cls):
        cls._Fe2_behaviour = "buffered"
        cls._stepsize_equilibration = 0.002
        cls._stepsize_crystallisation = 0.05
        cls._decrease_factor = 5
        cls.FeO_converge = 0.05
        cls.Kd_converge = 1e-3
