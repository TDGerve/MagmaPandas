from fractions import Fraction
from typing import List

import elementMass as e
import pandas as pd


def _get_elements(composition):

    import MagmaPandas as mp

    if isinstance(composition, mp.MagmaFrame):
        return composition.columns

    if isinstance(composition, mp.MagmaSeries):
        return composition.index


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


def _anhydrous_composition(composition):

    composition_H2O = composition.copy()

    try:
        composition_H2O = composition_H2O.drop("H2O")
    except KeyError:
        composition_H2O = composition_H2O.drop(columns=["H2O"])

    return composition_H2O.recalculate().normalise()
