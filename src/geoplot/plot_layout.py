import elementMass as e
import matplotlib.pyplot as plt

markers = ["^", "D", "s", "o", "p", "8", "h", "d", "4", "5"]


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

    hollywood = plt.cycler(
        color=[
            "#8ECAE6",
            "#FD9E02",
            "#219EBC",
            "#FB8500",
            "#126782",
            "#FFB703",
            "#023047",
        ],
    )

    campfire = plt.cycler(color=["#588C7E", "#F2E394", "#F2AE72", "#D96459", "#8C4646"])

    pastel = plt.cycler(
        color=[
            "#FAD2E1",
            "#BCD4E6",
            "#C5DEDD",
            "#99C1DE",
            "#EDDCD2",
            "#DBE7E4",
            "#FFF1E6",
            "#F0EFEB",
            "#FDE2E4",
            "#D6E2E9",
        ],
    )

    autumn = plt.cycler(
        color=[
            "#797D62",
            "#FFCB69",
            "#9B9B7A",
            "#E8AC65",
            "#BAA587",
            "#D08C60",
            "#D9AE94",
            "#B58463",
            "#997B66",
            "#F1DCA7",
        ],
    )

    rainbow = plt.cycler(
        color=[
            "#54478C",
            "#83E377",
            "#2C699A",
            "#B9E769",
            "#048BA8",
            "#EFEA5A",
            "#0DB39E",
            "#F1C453",
            "#16DB93",
            "#F29E4C",
        ]
    )

    matteblue = plt.cycler(
        color=["#666A86", "#788AA3", "#92B6B1", "#B2C9AB", "#E8DDB5"],
    )


def layout(colors=colors.hollywood, fontsize=8, **kwargs):
    axTitleSize = int(fontsize)
    axLabelSize = int(fontsize)
    tickLabelSize = int(fontsize)
    markersize = kwargs.get("markersize", 6)
    linewidth = kwargs.get("linewidth", 1)

    plt.rcParams["figure.constrained_layout.use"] = True

    plt.rc("figure", figsize=(8, 7), facecolor="white")

    # Text
    plt.rc("font", family="sans-serif", size=fontsize)

    # Legend
    plt.rc(
        "legend",
        fontsize=int(fontsize * 0.8),
        frameon=False,
        # facecolor="white",
    )

    # Axes
    plt.rc("xtick", direction="in", labelsize=tickLabelSize)
    plt.rc("ytick", direction="in", labelsize=tickLabelSize)
    plt.rc(
        "axes",
        grid=True,
        titlesize=axTitleSize,
        labelsize=axLabelSize,
        axisbelow=True,
        prop_cycle=colors + plt.cycler(marker=markers[: len(colors)]),
        facecolor="whitesmoke",
        linewidth=0.8,
    )
    plt.rc("grid", color="snow")

    # Lines
    plt.rc(
        "lines",
        markersize=markersize,
        linewidth=linewidth,
        markeredgewidth=0.5,
        mec="k",
    )


def side_plots(
    ax,
    x_axis: bool = False,
    y_axis: bool = False,
    side=0.1,
    spacing=1.5e-2,
    y_position="right",
    **kwargs,
):
    axes = []

    if x_axis:
        ax_kde_x = ax.inset_axes([0, 1.0 + spacing, 1, side])
        axes.append(ax_kde_x)

        ax_kde_x.spines["left"].set_visible(False)
        ax_kde_x.yaxis.set_visible(False)
        ax_kde_x.tick_params(axis="x", labelbottom=False, direction="in", **kwargs)
        ax_kde_x.set_xlim(ax.get_xlim())
        ax_kde_x.spines["right"].set_visible(False)

    if y_axis:
        if y_position == "right":
            ax_kde_y = ax.inset_axes([1.0 + spacing, 0, side, 1])
            ax_kde_y.spines["right"].set_visible(False)
            ax_kde_y.tick_params(axis="y", labelright=False, direction="in", **kwargs)

        if y_position == "left":
            ax_kde_y = ax.inset_axes([0 - spacing - side, 0, side, 1])
            ax_kde_y.spines["left"].set_visible(False)
            ax_kde_y.tick_params(axis="y", labelleft=False, direction="in", **kwargs)
            ax_kde_y.invert_xaxis()
            ax_kde_y.get_yaxis().tick_right()

        axes.append(ax_kde_y)

        ax_kde_y.tick_params(axis="y", labelleft=False, direction="in", **kwargs)
        ax_kde_y.xaxis.set_visible(False)
        ax_kde_y.spines["bottom"].set_visible(False)
        # ax_kde_y.tick_params(axis="y", labelleft=False, direction="in", **kwargs)
        ax_kde_y.set_ylim(ax.get_ylim())

    for axis in axes:
        # axis.set_frame_on(False)
        axis.grid(False)
        axis.set_facecolor(color="white")

        for spine in ["top"]:
            axis.spines[spine].set_visible(False)

    return axes


def subscript_numbers(compound: str):
    string = ""
    elements = e.elements._find_elements(compound)
    for element in elements:
        name, quantity = e.elements._find_quantity(element)
        string = string + name
        if float(quantity) > 1:
            string = string + f"$_{quantity}$"

    return string
