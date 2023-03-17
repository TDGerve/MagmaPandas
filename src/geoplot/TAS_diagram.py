from importlib import resources

import matplotlib.pyplot as plt
import pandas as pd


def TAS(ax, labels=False, fontsize="medium", **kwargs):
    """Returns a line plot element of classification of volcanic rocks
    in total-alkali vs silica plots
    """

    with resources.open_text("geoplot.data", "TAS.csv") as df:
        TAS = pd.read_csv(df)

    rock_labels = {
        "Picro-basalt": ["Picro\nbasalt", [41.7, 1.5]],
        "Basalt": ["Basalt", [47, 2.5]],
        "Basaltic andesite": ["Basaltic\nandesite", [53, 2.5]],
        "Andesite": ["Andesite", [58, 2.5]],
        "Dacite": ["Dacite", [65.5, 4]],
        "Trachy-basalt": ["Trachy-\nbasalt", [47.5, 5.5]],
        "Basaltic trachy-andesite": ["Basaltic\ntrachy-\nandesite", [51.6, 6.5]],
        "Trachy-andesite": ["Trachy-\nandesite", [56, 8]],
        "Trachyte": ["Trachyte", [64, 11]],
        "Tephrite": ["Tephrite", [43.5, 7]],
        "Phono-tephrite": ["Phono-\ntephrite", [47, 9.0]],
        "Tephri-phonolite": ["Tephri-\nphonolite", [51, 11]],
        "Phonolite": ["Phonolite", [55, 15]],
        "Foidite": ["Foidite", [45, 14]],
        "Rhyolite": ["Rhyolite", [72, 8.5]],
    }

    if labels:
        for _, rock in rock_labels.items():
            ax.text(
                *rock[1],
                rock[0],
                fontsize=fontsize,
                fontfamily="monospace",
                clip_on=True
            )

    for id in TAS.id.unique():
        ax.plot(
            TAS.loc[TAS.id == id, "x"],
            TAS.loc[TAS.id == id, "y"],
            "-",
            color="k",
            **kwargs
        )
