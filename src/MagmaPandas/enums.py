from enum import Enum


class Unit(Enum):
    MOL_FRACTIONS = "mol fraction"
    WT_PERCENT = "wt. %"
    PPM = "ppm"
    UNKNOWN = None


class Datatype(Enum):
    CATION = "cation"
    OXIDE = "oxide"
    UNKNOWN = None
