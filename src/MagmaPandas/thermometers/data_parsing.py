def _check_oxides(composition, oxides_needed):

    oxides = _get_oxides(composition)

    absentOxides = oxides_needed.difference(oxides)

    if len(absentOxides) > 0:
        raise KeyError(f"{absentOxides} not found in melt")


def _fill_data(composition):

    composition_filled = composition.copy(deep=True)
    composition_filled = composition_filled.fillna(0.0)

    return composition_filled


def _get_oxides(composition):

    import MagmaPandas as mp

    if isinstance(composition, mp.MagmaFrame):
        return composition.columns

    if isinstance(composition, mp.MagmaSeries):
        return composition.index


def parse_data(composition, oxides_needed):

    _check_oxides(composition, oxides_needed)

    return _fill_data(composition)


def _anhydrous_composition(composition):

    composition_H2O = composition.copy()

    try:
        composition_H2O = composition_H2O.drop("H2O")
    except KeyError:
        composition_H2O = composition_H2O.drop(columns=["H2O"])

    return composition_H2O.recalculate()
