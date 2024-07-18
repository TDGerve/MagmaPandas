import warnings as w
from fractions import Fraction
from typing import Dict, List

import elementMass as e
import numpy as np
import pandas as pd

from MagmaPandas.MagmaFrames.protocols import Magma


def _check_calibration_range_Series(melt: pd.Series, calibration_range: Dict) -> None:

    checks = np.array([], dtype=bool)
    for element, *range in calibration_range:
        row = [element] if isinstance(element, str) else element
        np.append(checks, not (range[0] < melt[row].sum() < range[1]))

    if sum(checks) < 1:
        return

    w.warn(f"sample has composition outside the thermometer calibration range")

    # SiO2_range = not (31 < melt["SiO2"] < 73.64)
    # alkali_range = melt[["Na2O", "K2O"]].sum() > 14.3
    # H2O_range = melt["H2O"] > 18.6

    # outside_range = SiO2_range | alkali_range | H2O_range

    # if not outside_range:
    #     return
    # w.warn(f"sample has composition outside the thermometer calibration range")


def _check_calibration_range_dataframe(melt: Magma, calibration_range: Dict) -> None:

    checks = pd.DataFrame(dtype=bool)

    for element, *range in calibration_range:
        col = [element] if isinstance(element, str) else element
        checks[str(element)] = ~melt[col].sum(axis=1).between(range[0], range[1])

    outside_range = checks.any(axis=1)

    if outside_range.sum() < 1:
        return

    w.warn(
        f"samples {*melt.index[outside_range],} have compositions outside the thermometer calibration range"
    )

    # SiO2_range = ~melt["SiO2"].between(31, 73.64)
    # alkali_range = melt[["Na2O", "K2O"]].sum(axis=1) > 14.3
    # H2O_range = melt["H2O"] > 18.6

    # outside_range = SiO2_range | alkali_range | H2O_range

    # if outside_range.sum() < 1:
    #     return

    # w.warn(
    #     f"samples {*melt.index[outside_range],} have compositions outside the thermometer calibration range"
    # )


def _check_calibration_range(melt: pd.Series | Magma, calibration_range) -> None:

    if isinstance(melt, pd.Series):
        _check_calibration_range_Series(melt, calibration_range)
    elif isinstance(melt, pd.DataFrame):
        _check_calibration_range_dataframe(melt, calibration_range)


def _get_elements(composition):

    import MagmaPandas as mp

    if isinstance(composition, mp.MagmaFrame):
        return composition.columns

    if isinstance(composition, mp.MagmaSeries):
        return composition.index


def _anhydrous_composition(composition):

    composition_H2O = composition.copy()

    try:
        composition_H2O = composition_H2O.drop("H2O")
    except KeyError:
        composition_H2O = composition_H2O.drop(columns=["H2O"])

    return composition_H2O.recalculate()


def _remove_elements(composition, drop: List[str]):

    composition_new = composition.copy()

    for e in drop:

        if e in composition_new.index:
            composition_new = composition_new.drop(e)
        try:
            if e in composition_new.columns:
                composition_new = composition_new.drop(columns=[e])
        except AttributeError:
            pass

        composition_new.recalculate(inplace=True)

    return composition_new


def moles_per_oxygen(moles):

    axis = [0, 1][isinstance(moles, pd.DataFrame)]

    O_number = e.oxygen_numbers(moles.elements)

    no_oxygen = O_number.index[O_number == 0]
    composition_new = _remove_elements(composition=moles, drop=no_oxygen)

    O_number = O_number.drop(no_oxygen)
    cations = e.cation_names(composition_new.elements)
    cation_number = e.cation_numbers(composition_new.elements)

    cations_normalised = cation_number.div(O_number)

    names_new = {
        ox: f"{cat}{Fraction(num).limit_denominator()}O"
        for cat, num, ox in zip(cations, cations_normalised, composition_new.elements)
    }
    names_new["total"] = "O_total"

    composition_new = (
        composition_new[composition_new.elements]
        .mul(O_number, axis=axis)
        .normalise(to=1)
    )

    return composition_new.rename(names_new, axis=axis)
