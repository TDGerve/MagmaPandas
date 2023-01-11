import pandas as pd
from multiprocessing import Pool
from scipy.optimize import root_scalar
from alive_progress import alive_bar

from MagmaPandas.MagmaFrames import Olivine, Melt
from MagmaPandas.Kd.ol_melt import calculate_olivine_Kd


def correct_olivine_accumulation(liquids: Melt, olivines: Olivine, pressures, **kwargs):
    """
    Subtract olivine from melt compositions until observed and model Mg-Fe Kds match.
    Only corrects compositions where Kd(observed) > Kd(equilibrium)
    """

    converge = kwargs.get("converge", 0.05)

    melts = liquids[liquids.elements].copy()

    olivines = olivines.reindex(columns=melts.columns).recalculate()
    olivines = olivines.fillna(0.0)

    forsterite = olivines.forsterite
    Kd_eq, Kd_observed = calculate_olivine_Kd(melts, forsterite, pressures)

    Kd_delta = abs(Kd_observed - Kd_eq)

    names = Kd_delta.index[Kd_delta > converge]

    total = len(names)
    samples = [
        (name, melts.loc[[name]], olivines.loc[[name]], pressures[name])
        for name in names
    ]

    olivine_accumulated = pd.Series(0.0, index=melts.index, dtype=float)

    with Pool() as pool, alive_bar(
        total,
        spinner=None,
        length=30,
        manual=True,
        theme="smooth",
        force_tty=True,
    ) as bar:

        results = pool.imap_unordered(_correct_accumulation_multicore, samples)
        pool.close()

        finished = 0.0
        bar(finished / total)

        for name, olivine_amount in results:
            olivine_accumulated[name] = olivine_amount
            finished += 1
            bar(finished / total)

    melts_corrected = melts.moles.loc[names] + olivines.moles.loc[names].mul(
        olivine_accumulated[names], axis=0
    )
    melts_corrected = melts_corrected.convert_moles_wtPercent

    difference = melts.index.difference(melts_corrected.index)
    melts_corrected = pd.concat(
        [melts_corrected, melts.normalise().loc[difference]], axis=0
    )

    return olivine_accumulated, melts_corrected


def _root_Kd(exchange_amount, melt, olivine, pressure):

    melt_moles = melt.moles
    olivine_moles = olivine.moles

    forsterite = olivine_moles.forsterite

    melt_new = melt_moles + olivine_moles.mul(exchange_amount)
    melt_new = melt_new.normalise()
    melt_new = melt_new.convert_moles_wtPercent
    Kd_equilibrium, Kd_real = calculate_olivine_Kd(melt_new, forsterite, pressure)

    return (Kd_equilibrium - Kd_real).iloc[0]


def _calculate_accumulation(melt, olivine, P_bar):

    olivine_accumulated = root_scalar(
        _root_Kd, args=(melt, olivine, P_bar), x0=-0.01, x1=-0.05
    ).root

    return olivine_accumulated


def _correct_accumulation_multicore(sample):
    """
    Refactor of calculate_accumulation for multiprocess calling
    """

    name, melt, olivine, P_bar = sample

    olivine_accumulated = _calculate_accumulation(melt, olivine, P_bar)

    return name, olivine_accumulated
