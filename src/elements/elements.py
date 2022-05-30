import re
from typing import List
from importlib import resources
import pandas as pd


def periodic_table(masses_only: bool = True):
    """
    Docstring
    """

    if not masses_only:
        with resources.open_text("elements.data", "PeriodicTable.csv") as df:
            return pd.read_csv(df, index_col=["AtomicNumber"])

    with resources.open_text("elements.data", "PeriodicTable.csv") as df:
        return pd.read_csv(
            df, usecols=["Symbol", "AtomicMass"], index_col=["Symbol"]
        ).squeeze("columns")


def find_elements(compound: str):
    """
    Docstring
    """
    elements = re.findall("([A-Z][^A-Z]*)", compound)

    # Raise an error if no elements are found
    if len(elements) < 1:
        raise ValueError(f"'{compound}' does not contain valid elements")

    for element in elements:
    # Raise an error for invalid elements with more than 1 lower case character
        if sum(c.islower() for c in element) > 1:
            raise ValueError(f"'{element}' is not a valid element")
    # Check for non-word characters
        elif len(re.findall(r"\W", element)) > 0:
            raise ValueError(f"'{element}' contains an invalid character")

    # Raise an error if there are any leftover characters
    length_elements = sum(len(s) for s in elements)
    if len(compound) != length_elements:
        raise ValueError(f"There are leftover characters in '{compound}'; elements found: {elements}")

    return elements


def find_quantity(element: str):
    """
    Docstring
    """

    element_quantity = re.findall(r"(\D+|\d[^A-Z]*)", element)

    if len(element_quantity) < 2:
        element_quantity.append(1)

    return element_quantity


def decompose(compound: str):
    """
    Docstring
    """

    elements = [find_quantity(i) for i in find_elements(compound)]

    elements_pd = pd.DataFrame(elements, columns=["element", "quantity"]).set_index(
        "element"
    )

    return elements_pd.astype(float).squeeze("columns")


def calculate_weight(compound: str):
    """
    Docstring
    """

    atomic_weights = periodic_table()

    elements = decompose(compound)

    return (atomic_weights[elements.index] * elements).sum()


def compound_weights(compounds: List[str]):
    """
    Docstring
    """

    weights = pd.Series(index=compounds, name="weights", dtype="float64")

    for i in weights.index:
        weights[i] = calculate_weight(i)

    return weights


def cation_numbers(compounds: List[str]):
    """
    Docstring
    """

    cations = pd.Series(index=compounds, name='cations', dtype=int)

    for i in cations.index:
        cations[i] = decompose(i)[0]

    return cations


def oxygen_numbers(compounds: List[str]):
    """
    Docstrings
    """

    oxygen = pd.Series(index=compounds, name='oxygen', dtype=int)

    for i in oxygen.index:
        oxygen[i] = decompose(i)['O']

    return oxygen


def cation_names(compounds: List[str]):
    """
    Docstrings
    """

    return [decompose(oxide).index[0] for oxide in compounds]
