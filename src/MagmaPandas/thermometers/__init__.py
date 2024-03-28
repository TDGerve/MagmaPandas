"""
Module with melt-only and melt-mineral thermometers
"""

import inspect

from MagmaPandas.thermometers import melt, ol_melt

melt_thermometers = {
    func[0]: func[1]
    for func in inspect.getmembers(melt, inspect.isfunction)
    if (func[0][0] != "_") and (func[1].__module__ == melt.__name__)
}

ol_melt_thermometers = {
    func[0]: func[1]
    for func in inspect.getmembers(ol_melt, inspect.isfunction)
    if (func[0][0] != "_") and (func[1].__module__ == ol_melt.__name__)
}
