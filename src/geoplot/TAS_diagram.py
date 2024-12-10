from importlib import resources
from typing import Dict

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd


def TAS(
    ax: plt.axis,
    labels=False,
    fontsize="medium",
    abbreviate=False,
    textkwargs: Dict = {},
    **kwargs,
):
    """Returns a line plot element of classification of volcanic rocks
    in total-alkali vs silica plots
    """

    # resources.open_text("geoplot.data", "TAS.csv")
    with resources.files("geoplot.data").joinpath("TAS.csv").open("r") as df:
        TAS = pd.read_csv(df)

    rock_labels = {
        "Picro-basalt": ["Picro\nbasalt", "Pi-Ba", [43, 1.5]],
        "Basalt": ["Basalt", "Ba", [49, 2.5]],
        "Basaltic andesite": ["Basaltic\nandesite", "Ba-An", [54.5, 2.5]],
        "Andesite": ["Andesite", "An", [60, 2.5]],
        "Dacite": ["Dacite", "Da", [68, 4]],
        "Trachy-basalt": ["Trachy-\nbasalt", "Tr-Ba", [49, 5.5]],
        "Basaltic trachy-andesite": [
            "Basaltic\ntrachy-\nandesite",
            "B\nTr-An",
            [52.5, 6.5],
        ],
        "Trachy-andesite": ["Trachy-\nandesite", "TrAn", [58.5, 8]],
        "Trachyte": ["Trachyte", "Tr", [64, 11]],
        "Tephrite": ["Tephrite", "Te", [45, 7]],
        "Phono-tephrite": ["Phono-\ntephrite", "Ph-Te", [49, 9.0]],
        "Tephri-phonolite": ["Tephri-\nphonolite", "Te-Ph", [53, 11]],
        "Phonolite": ["Phonolite", "Ph", [57, 15]],
        "Foidite": ["Foidite", "Fo", [45, 14]],
        "Rhyolite": ["Rhyolite", "Rh", [75, 8.5]],
    }

    if labels:
        for _, (fullname, abbreviation, coords) in rock_labels.items():
            label = [fullname, abbreviation][abbreviate]
            ax.text(
                *coords,
                label,
                fontsize=fontsize,
                fontfamily="monospace",
                clip_on=True,
                horizontalalignment="center",
                **textkwargs,
            )

    for id in TAS.id.unique():
        ax.plot(
            TAS.loc[TAS.id == id, "x"],
            TAS.loc[TAS.id == id, "y"],
            "-",
            color="k",
            **kwargs,
        )
