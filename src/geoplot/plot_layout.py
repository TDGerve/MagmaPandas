import matplotlib.pyplot as plt
import elements as e

# Color palettes
class colors:
    """
    Color palettes for plots
    """

    flatDesign = plt.cycler(
        color=["#e27a3d", "#344d5c", "#df5a49", "#43b29d", "#efc94d"]
    )

    firenze = plt.cycler(color=["#8E2800", "#468966", "#B64926", "#FFF0A5", "#FFB03B"])

    vitaminC = plt.cycler(color=["#FD7400", "#004358", "#FFE11A", "#1F8A70", "#BEDB39"])

    bella = plt.cycler(color=["#801637", "#047878", "#FFB733", "#F57336", "#C22121"])

    buddha = plt.cycler(color=["#192B33", "#FF8000", "#8FB359", "#FFD933", "#CCCC52"])

    elemental = plt.cycler(
        color=["#E64661", "#FFA644", "#998A2F", "#2C594F", "#002D40"]
    )

    carolina = plt.cycler(color=["#73839C", "#2E4569", "#AECCCF", "#D5957D", "#9C7873"])

    fourtyTwo = plt.cycler(
        color=["#2469A6", "#C4E1F2", "#F2E205", "#F2D22E", "#D9653B"]
    )

    terrazaverde = plt.cycler(
        color=["#DFE2F2", "#88ABF2", "#4384D9", "#56BFAC", "#D9B341"]
    )


markers = plt.cycler(marker=["^", "D", "s", "o", "p"])


def layout(colors=colors.firenze, fontsize=14, **kwargs):
    axTitleSize = int(fontsize / 1.2)
    axLabelSize = int(fontsize)
    tickLabelSize = int(fontsize / 1.2)
    markersize = kwargs.get("markersize", 8)
    linewidth = kwargs.get("linewidth", 2)

    plt.rcParams["figure.constrained_layout.use"] = True
    plt.rcParams["savefig.dpi"] = 300

    plt.rc("figure", figsize=(8, 7), facecolor="white")

    # Text
    plt.rc("font", family="sans-serif", size=fontsize)

    # Legend
    plt.rc("legend", fontsize=fontsize / 1.5, fancybox=False, facecolor="white")

    # Axes
    plt.rc("xtick", direction="in", labelsize=tickLabelSize)
    plt.rc("ytick", direction="in", labelsize=tickLabelSize)
    plt.rc(
        "axes",
        grid=True,
        titlesize=axTitleSize,
        labelsize=axLabelSize,
        axisbelow=True,
        prop_cycle=colors + markers,
        facecolor="whitesmoke",
        linewidth=1.2,
    )
    plt.rc("grid", color="snow")

    # Lines
    plt.rc("lines", markersize=markersize, linewidth=linewidth)


def side_plots(
    ax, x_axis: bool = False, y_axis: bool = False, side=0.1, spacing=0.03, **kwargs
):

    axes = []

    if x_axis:
        ax_kde_x = ax.inset_axes([0, 1.0 + spacing, 1, side])
        axes.append(ax_kde_x)

        ax_kde_x.spines["left"].set_visible(False)
        ax_kde_x.yaxis.set_visible(False)
        ax_kde_x.tick_params(axis="x", labelbottom=False, direction="in", **kwargs)
        ax_kde_x.set_xlim(ax.get_xlim())

    if y_axis:
        ax_kde_y = ax.inset_axes([1.0 + spacing, 0, side, 1])
        axes.append(ax_kde_y)

        ax_kde_y.xaxis.set_visible(False)
        ax_kde_y.spines["bottom"].set_visible(False)
        ax_kde_y.tick_params(axis="y", labelleft=False, direction="in", **kwargs)
        ax_kde_y.set_ylim(ax.get_ylim())

    for axis in axes:
        # axis.set_frame_on(False)
        axis.grid(False)
        axis.set_facecolor(color="white")

        for spine in ["right", "top"]:
            axis.spines[spine].set_visible(False)

    return axes


def subscript_numbers(compound: str):
    string = ""
    elements = e.find_elements(compound)
    for element in elements:
        name, quantity = e.find_quantity(element)
        string = string + name
        if float(quantity) > 1:
            string = string + f"$_{quantity}$"

    return string
