import inspect
import sys

from MagmaPandas.volatile_solubility.volatile_solubility_models import (
    Allison2022,
    IaconoMarziano,
    Shiskina,
)

_members = inspect.getmembers(sys.modules[__name__], inspect.ismodule)


volatile_solubility_models_dict = {
    cls[0]: cls[1] for cls in _members if cls[1].__package__ == __name__
}
