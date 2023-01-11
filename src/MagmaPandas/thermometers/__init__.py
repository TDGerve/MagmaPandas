import inspect

from . import melt

melt_thermometers = {
    func[0]: func[1] for func in inspect.getmembers(melt, inspect.isfunction)
}
