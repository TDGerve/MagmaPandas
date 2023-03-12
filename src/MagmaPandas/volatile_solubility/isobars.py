import itertools
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import MagmaPandas.volatile_solubility as vs
from MagmaPandas.MagmaSeries import MagmaSeries


def CO2H2O_isobar_data(
    oxide_wtPercents: MagmaSeries,
    temperature: float,
    isobars=None,
    isopleths=None,
    interpolate=True,
    n_points=10,
):

    species = "mixed"

    if isobars is None:
        isobars = np.arange(1e3, 7e3, 1e3)

    if isopleths is None:
        isopleths = np.arange(0.0, 1.1, 0.1)

    isobar_data = pd.DataFrame(
        index=np.arange(0, 1 + 1 / n_points, 1 / n_points),
        columns=pd.MultiIndex.from_tuples(itertools.product(isobars, ["H2O", "CO2"])),
        dtype=float,
    )
    isobar_data.index.name = "xfl"
    isobar_data.columns.names = ("P_bar", "species")

    isopleth_data = pd.DataFrame(
        index=np.arange(0, max(isobars) + 1 / n_points, max(isobars) / n_points),
        columns=pd.MultiIndex.from_tuples(itertools.product(isopleths, ["H2O", "CO2"])),
        dtype=float,
    )
    isopleth_data.index.name = "P_bar"
    isopleth_data.columns.names = ("xfl", "species")

    samples_isobars = [
        pair
        for pair in itertools.product(
            [oxide_wtPercents],
            isobars,
            [temperature],
            isobar_data.index.values,
            [species],
        )
    ]
    samples_isopleths = [
        pair
        for pair in itertools.product(
            [oxide_wtPercents],
            isopleth_data.index.values,
            [temperature],
            isopleths,
            [species],
        )
    ]

    with Pool() as pool:

        results_isobars = pool.imap_unordered(_solubility_multicore, samples_isobars)

        for (xfl, P_bar), volatiles in results_isobars:
            isobar_data.loc[xfl, P_bar] = volatiles

        results_isopleths = pool.imap_unordered(
            _solubility_multicore, samples_isopleths
        )

        for (xfl, P_bar), volatiles in results_isopleths:
            isopleth_data.loc[P_bar, xfl] = volatiles

    if interpolate:
        return _isobar_interpolate(isobar_data, isopleth_data)

    return isobar_data, isopleth_data


def _solubility_multicore(sample):

    oxide_wtPercents, P_bar, temperature, xfl, species = sample

    volatiles = vs.calculate_solubility(
        oxide_wtPercents, P_bar=P_bar, T_K=temperature, x_fluid=xfl, species=species
    )

    return (xfl, P_bar), volatiles


def _isobar_interpolate(isobar_data, isopleth_data):

    isobars_smooth = pd.DataFrame(dtype=float, columns=["H2O", "CO2", "P_bar"])

    for p in isobar_data.columns.levels[0]:
        x_data = isobar_data[(p, "H2O")]
        y_data = isobar_data[(p, "CO2")]
        f = interp1d(x=x_data, y=y_data, kind="cubic")
        x_interp = np.linspace(min(x_data), max(x_data), len(x_data) * 10)
        data_interp = pd.DataFrame({"H2O": x_interp, "CO2": f(x_interp), "P_bar": p})
        isobars_smooth = pd.concat([isobars_smooth, data_interp])

    isopleths_smooth = pd.DataFrame(dtype=float, columns=["H2O", "CO2", "xfl"])

    trim = (0 < isopleth_data.columns.levels[0]) & (isopleth_data.columns.levels[0] < 1)

    for x in isopleth_data.columns.levels[0][trim]:
        x_data = isopleth_data[(x, "H2O")]
        y_data = isopleth_data[(x, "CO2")]
        f = interp1d(x=x_data, y=y_data, kind="quadratic")
        x_interp = np.linspace(min(x_data), max(x_data), len(x_data) * 10)
        data_interp = pd.DataFrame({"H2O": x_interp, "CO2": f(x_interp), "xfl": x})
        isopleths_smooth = pd.concat([isopleths_smooth, data_interp])

    if 0.0 in isopleth_data.columns.levels[0]:
        CO2max = isobars_smooth.loc[isobars_smooth["H2O"] == 0.0, "CO2"].max()
        xfl0 = pd.DataFrame({"H2O": (0, 0), "CO2": (0, CO2max), "xfl": (0.0, 0.0)})
        isopleths_smooth = pd.concat([isopleths_smooth, xfl0], axis=0)

    if 1.0 in isopleth_data.columns.levels[0]:
        H2Omax = isobars_smooth.loc[isobars_smooth["CO2"] == 0.0, "H2O"].max()
        xfl1 = pd.DataFrame({"H2O": (0, H2Omax), "CO2": (0.0, 0.0), "xfl": (1.0, 1.0)})
        isopleths_smooth = pd.concat([isopleths_smooth, xfl1], axis=0)

    return isobars_smooth, isopleths_smooth
