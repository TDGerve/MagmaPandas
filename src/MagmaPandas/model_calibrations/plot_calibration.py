import json
import os
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes

import geoplot as gp
from geoplot.plot_layout import _normalise_kde
from MagmaPandas.configuration import (
    Fe3Fe2_models,
    Kd_ol_FeMg_models,
    configuration,
    melt_thermometers,
    volatile_solubility_models,
)
from MagmaPandas.parse_io.validate import _check_argument, _ensure_str_list

_models_dict = {
    "Fe3Fe2": Fe3Fe2_models,
    "Kd": Kd_ol_FeMg_models,
    "melt_thermometer": melt_thermometers,
    "volatile_solubility": volatile_solubility_models,
}

_here = os.path.dirname(__file__)
_calibration_datasets_path = os.path.join(_here, "data", "calibration_datasets.json")


@_check_argument(
    var_name="parameter",
    allowed_values=list(_models_dict.keys()),
)
def _read_calibration_data(parameter):
    """ """

    calibration_data_path = os.path.join(
        _here, "data", f"{parameter}_calibration_data.csv"
    )
    if not os.path.exists(calibration_data_path):
        raise RuntimeError(f"This function is not yet implemented for {parameter}.")

    with open(_calibration_datasets_path, "r") as f:
        calibration_datasets = json.load(f)

    calibration_data = pd.read_csv(calibration_data_path)

    return calibration_data, calibration_datasets[parameter]


def get_calibration_data(parameter, model=None):

    if model is None:
        parameter_name = (
            f"{parameter}_model" if parameter in ("Fe3Fe2", "Kd") else parameter
        )

        model = getattr(configuration, parameter_name)

    calibration_data, calibration_datasets = _read_calibration_data(parameter=parameter)
    datasets = calibration_datasets[model]

    if datasets is None:
        raise RuntimeError(
            f"The calibration dataset is not available for {model}. Please refer to the original publication."
        )

    data = calibration_data.query("ref in @datasets")

    return data


def plot_calibration(
    parameter: str,
    x_elements: str | List[str],
    y_elements: str | List[str],
    models: None | str = None,
    ax: Axes = None,
    melts: None | pd.DataFrame = None,
    sidekde: bool = True,
    calibration_props: Dict = {},
    melt_props: Dict = {},
):
    """"""

    if ax is not None:
        return _prepare_plot_calibration(
            parameter=parameter,
            x_elements=x_elements,
            y_elements=y_elements,
            ax=ax,
            melts=melts,
            models=models,
            sidekde=sidekde,
            calibration_props=calibration_props,
            melt_props=melt_props,
        )

    if isinstance(x_elements, str):
        x_elements = [x_elements]
    if isinstance(y_elements, str):
        y_elements = [y_elements]

    mm = 1 / 25.4
    gp.layout(colors=gp.colors.campfire)

    fig, ax = plt.subplots(figsize=(85 * mm, 80 * mm))

    _prepare_plot_calibration(
        parameter=parameter,
        x_elements=x_elements,
        y_elements=y_elements,
        ax=ax,
        melts=melts,
        models=models,
        sidekde=sidekde,
        calibration_props=calibration_props,
        melt_props=melt_props,
    )

    ax.set_xlabel(f"wt.% {'+'.join(gp.subscript_numbers(i) for i in x_elements)}")
    ax.set_ylabel(f"wt.% {'+'.join(gp.subscript_numbers(i) for i in y_elements)}")

    ax.legend(frameon=True, fancybox=False)

    plt.show()


def plot_calibration_PT(
    parameter: str,
    models: None | List[str] = None,
    ax: None | Axes = None,
    sidekde=True,
    **kwargs,
):
    """ """

    calibration_data, calibration_datasets = _read_calibration_data(parameter=parameter)

    if models is None:
        models = _models_dict[parameter]

    if isinstance(models, str):
        models = [models]
    _ensure_str_list(names=["models"], vals=[models])

    mm = 1 / 25.4
    plot = False

    if ax is None:
        plot = True
        gp.layout(colors=gp.colors.bright)
        fig, ax = plt.subplots(figsize=(90 * mm, 85 * mm))

    if sidekde:
        sax_x, sax_y = gp.side_plots(ax=ax, y_axis=True)

    for model in models:

        if model == "fixed":
            continue

        datasets = calibration_datasets[model]
        data = calibration_data.query("ref in @datasets")
        points = ax.plot(
            data["T_K"] - 273.15,
            data["P_bar"] / 1e3,
            lw=0.0,
            label=model,
            **kwargs,
        )
        color = points[0].get_markerfacecolor()

        if sidekde:
            kde_x = sns.kdeplot(
                ax=sax_x,
                x=data["T_K"] - 273.15,
                fill=True,
                color=color,
                multiple="stack",
                warn_singular=False,
            )

            # _normalise_kde(kde=kde_x, axis="x"

            kde_y = sns.kdeplot(
                ax=sax_y,
                y=data["P_bar"] / 1e3,
                fill=True,
                color=color,
                multiple="stack",
                warn_singular=False,
            )
            # _normalise_kde(kde=kde_y, axis="y")

    if not plot:
        return ax, (sax_x, sax_y)

    ax.invert_yaxis()

    ax.legend(frameon=True, fancybox=False)
    ax.set_xlabel("Temperature ($\degree$C)")
    ax.set_ylabel("Pressure (kbar)")

    plt.show()


def _prepare_plot_calibration(
    parameter: str,
    x_elements: str | List[str],
    y_elements: str | List[str],
    ax: Axes,
    models: None | List[str] = None,
    melts: None | pd.DataFrame = None,
    sidekde: bool = True,
    calibration_props: Dict = {},
    melt_props: Dict = {},
):
    """"""

    calibration_data, calibration_datasets = _read_calibration_data(parameter=parameter)

    if models is None:
        parameter_name = (
            f"{parameter}_model" if parameter in ("Fe3Fe2", "Kd") else parameter
        )

        models = getattr(configuration, parameter_name)

    _ensure_str_list(names=["x", "y", "models"], vals=[x_elements, y_elements, models])

    if isinstance(x_elements, str):
        x_elements = [x_elements]
    if isinstance(y_elements, str):
        y_elements = [y_elements]
    if isinstance(models, str):
        models = [models]

    results = [ax]

    if sidekde:
        (sax_x, sax_y) = gp.side_plots(ax=ax, x_axis=True, y_axis=True)
        results.append((sax_x, sax_y))

    for model in models:

        if model == "fixed":
            continue

        datasets = calibration_datasets[model]
        data = calibration_data.query("ref in @datasets")

        x_calibration_data = data[x_elements].sum(axis=1)
        y_calibration_data = data[y_elements].sum(axis=1)

        calibration_plot = ax.plot(
            x_calibration_data,
            y_calibration_data,
            # "D",
            lw=0.0,
            label=f"{model} dataset\n($n={data.shape[0]:d}$)",
            **calibration_props,
        )

        if sidekde:
            color_calibration = calibration_plot[0].get_markerfacecolor()

            sns.kdeplot(
                ax=sax_x,
                x=x_calibration_data,
                fill=True,
                color=color_calibration,
                multiple="stack",
            )
            sns.kdeplot(
                ax=sax_y,
                y=y_calibration_data,
                fill=True,
                color=color_calibration,
                multiple="stack",
            )

    if melts is None:
        return results

    x_melt_data = melts[x_elements].sum(axis=1)
    y_melt_data = melts[y_elements].sum(axis=1)

    melts_plot = ax.plot(x_melt_data, y_melt_data, lw=0.0, label="melts", **melt_props)

    results.append(melts_plot)

    if not sidekde:
        return results

    color_melt = melts_plot[0].get_markerfacecolor()

    sns.kdeplot(ax=sax_x, x=x_melt_data, fill=True, color=color_melt, multiple="stack")
    sns.kdeplot(ax=sax_y, y=y_melt_data, fill=True, color=color_melt, multiple="stack")

    return results
