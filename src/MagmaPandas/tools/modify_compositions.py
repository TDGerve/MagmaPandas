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

    for element in drop:

        if element in composition_new.index:
            composition_new = composition_new.drop(element)
        try:
            if element in composition_new.columns:
                composition_new = composition_new.drop(columns=[element])
        except AttributeError:
            pass

        composition_new.recalculate(inplace=True)

    return composition_new


def cation_moles_per_oxygen(moles):
    """
    Calculates cations moles per 1 mole total oxygen. This conversion is used by the thermobarometers and Kd model from :cite:t:`Sun2020a`.

    In the original paper this is described as 'oxide molar contents per unit oxygen'. However, Sun explained in emails that they are actually calculated as below, which is closer to cation moles per 1 mole total oxygen.
    """

    axis = [0, 1][isinstance(moles, pd.DataFrame)]
    # get number of oxygen per oxide
    O_number = e.oxygen_numbers(moles.elements)
    # remove non-oxides
    no_oxygen = O_number.index[O_number == 0]
    composition_new = _remove_elements(composition=moles, drop=no_oxygen)
    O_number = O_number.drop(no_oxygen)
    # get number of cations per oxide
    cation_number = e.cation_numbers(composition_new.elements)
    cations_per_oxygen = cation_number.div(O_number)

    # calculate oxide moles per 1 mole total oxygen
    composition_new = (
        composition_new[composition_new.elements]
        .mul(O_number, axis=axis)
        .normalise(to=1)
    )
    # convert to cations moles per 1 mole total oxygen
    composition_new = composition_new.mul(cations_per_oxygen, axis=axis)

    cations_names = e.cation_names(composition_new.elements)
    names_new = {
        ox: f"{cat}{Fraction(num).limit_denominator()}O"
        for cat, num, ox in zip(
            cations_names, cations_per_oxygen, composition_new.elements
        )
    }
    names_new["total"] = "O_total"

    return composition_new.rename(names_new, axis=axis)


def _anhydrous_composition(composition, normalise=True):

    axis = [0, 1][isinstance(composition, pd.DataFrame)]

    composition_H2O = composition.copy()

    try:
        composition_H2O = composition_H2O.drop("H2O", axis=axis)
    except KeyError:
        return composition.copy()

    if not normalise:
        return composition_H2O.recalculate()

    return composition_H2O.recalculate().normalise()
