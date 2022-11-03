import pandas as pd
from importlib import resources

REE_radii = pd.Series(
    {
        "La": 1.16,
        "Ce": 1.143,
        "Pr": 1.126,
        "Nd": 1.109,
        "Sm": 1.079,
        "Eu": 1.066,
        "Gd": 1.053,
        "Tb": 1.040,
        "Dy": 1.027,
        "Y": 1.019,
        "Ho": 1.015,
        "Er": 1.004,
        "Tm": 0.994,
        "Yb": 0.985,
        "Lu": 0.977,
    }
)

REE = REE_radii.index.values

divalent_radii = pd.Series({"Mg": 0.89, "Ba": 1.42, "Ca": 1.12, "Eu": 1.25, "Sr": 1.26})


def C1chondrite() -> pd.Series:
    """Returns C1 chondrite composition from McDonough & Sun (1995)"""

    with resources.open_text("geoplot.data", "Mcdonough_sun_1995.csv") as df:
        C1 = pd.read_csv(df, index_col=["reservoir"]).loc["C1"]

    return C1


def primitiveMantle() -> pd.Series:
    """Returns primitive mantle composition from McDonough & Sun (1995)"""

    with resources.open_text("geoplot.data", "Mcdonough_sun_1995.csv") as df:
        PM = pd.read_csv(df, index_col=["reservoir"]).loc["Pyrolite"]

    return PM


def MORB(type: str, add_majors=False, conf_interval=False) -> pd.Series:
    """Returns primitive mantle composition from McDonough & Sun (1995)"""

    with resources.open_text("geoplot.data", "MORB_Gale2013.csv") as df:
        MORB = pd.read_csv(df, index_col=["reservoir"])

    major_elements = [
        "MgO",
        "SiO2",
        "FeO",
        "CaO",
        "Na2O",
        "Al2O3",
        "TiO2",
        "K2O",
        "P2O5",
        "MnO",
    ]

    if not add_majors:
        MORB = MORB.drop(columns=major_elements)

    if not conf_interval:
        return MORB.loc[type]

    confidence95 = MORB.loc[f"{type}_conf95"]

    return MORB.loc[type], confidence95


def radii(valency):
    """Returns a dictionary with radii for common trace elements

    Parameters
    ----------
    valency : {2, 3, 'REE'}
        valency of requisted trace elements

    Returns
    -------
    dictionary
        dictionary with ionic radii in Angstrom
    """

    if valency not in [2, 3, "REE"]:
        raise Exception("valency should be 2, 3 or 'REE'")

    valency = {2: 0, 3: 1, "REE": 1}[valency]

    divalent = pd.Series({"Mg": 0.89, "Ba": 1.42, "Ca": 1.12, "Eu": 1.25, "Sr": 1.26})

    REE = pd.Series(
        {
            "La": 1.16,
            "Ce": 1.143,
            "Pr": 1.126,
            "Nd": 1.109,
            "Sm": 1.079,
            "Eu": 1.066,
            "Gd": 1.053,
            "Tb": 1.040,
            "Dy": 1.027,
            "Y": 1.019,
            "Ho": 1.015,
            "Er": 1.004,
            "Tm": 0.994,
            "Yb": 0.985,
            "Lu": 0.977,
        }
    )

    total = [divalent, REE]

    return total[valency]
