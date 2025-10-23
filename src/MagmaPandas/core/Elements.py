import re
from typing import List, Optional

import elementMass as e
import numpy as np
import pandas as pd

from MagmaPandas.parse_io.validate import _check_argument


class Oxide_compositions:
    def __init__(self):
        self.names = np.array([], dtype=str)
        self._oxide_names = np.array([], dtype=str)
        self._cation_names = np.array([], dtype=str)
        self._cation_amount = np.array([], dtype=int)
        self._oxygen_amount = np.array([], dtype=int)

    def oxygen_amount(self, names, type="oxide"):
        if type == "oxide":
            self._process_names(names)
        idx = self._get_indeces(names, type=type)
        return self._oxygen_amount[idx]

    def cation_amount(self, names, type="oxide"):
        if type == "oxide":
            self._process_names(names)
        idx = self._get_indeces(names, type=type)
        return self._cation_amount[idx]

    def cation_names(self, oxides):
        self._process_names(oxides)
        idx = self._get_indeces(oxides)
        return self._cation_names[idx]

    def calculate_element_numbers(self, oxides):
        oxides_list = list(oxides)
        # new_oxides = oxides[~np.in1d(oxides, self.names)]
        # if len(new_oxides) == 0:
        #     return

        self.names = np.append(self.names, oxides_list)
        self._cation_names = np.append(self._cation_names, e.cation_names(oxides_list))
        self._cation_amount = np.append(
            self._cation_amount, e.cation_numbers(oxides_list).astype(int).values
        )
        self._oxygen_amount = np.append(
            self._oxygen_amount, e.oxygen_numbers(oxides_list).astype(int).values
        )

    def _process_names(self, oxides: List[str]) -> None:
        # cations = [c for c in oxides if "O" not in c]
        # oxides = list(set(oxides).difference(set(cations)))

        self._process_oxides(oxides=oxides)

        # difference = list(oxides_new.difference(set(self.names)))
        # if len(difference) > 0:
        #     self.calculate_element_numbers(difference)

        # idx = self._get_indeces(oxides)
        # [
        #     int(np.where(self.names == o)[0])
        #     for o in oxides[np.isin(oxides, self.names)]
        # ]

        # return idx

    def _process_oxides(self, oxides: List[str]) -> None:
        oxides_new = set(oxides)

        difference = list(oxides_new.difference(set(self.names)))
        if len(difference) > 0:
            self.calculate_element_numbers(difference)

    @_check_argument("type", [None, "oxide", "cation"])
    def _get_indeces(self, names: List[str], type: Optional[str] = None) -> List[int]:
        if type is None:
            type = "oxide"

        names_total = {"oxide": self.names, "cation": self._cation_names}[type]

        names_arr = np.array(names)

        return [
            int(np.where(names_total == name)[0][0])
            for name in names_arr[np.isin(names_arr, names_total)]
        ]


class Element_weights:
    def __init__(self):
        self.element_names = np.array([], dtype=str)
        self.all_weights = np.array([], dtype=float)
        self._not_an_element = np.array([], dtype=str)

    def weight(self, elements):
        idx = self._process_names(elements)
        return self.all_weights[idx].copy()

    def names(self, elements):
        idx = self._process_names(elements)
        return self.element_names[idx].copy()

    def calculate_weights(self, elements):
        new_elements = np.array(elements)[~np.in1d(elements, self.element_names)]
        if len(new_elements) == 0:
            return
        for element in new_elements:
            try:

                element_name = (
                    re.sub(r"\d+", "", element) if "O" not in element else element
                )  # strip numbers

                # Calculate element/oxide weight
                self.all_weights = np.append(
                    self.all_weights, e.calculate_weight(element_name)
                )
                self.element_names = np.append(self.element_names, element)
            except (ValueError, KeyError):
                self._not_an_element = np.append(self._not_an_element, element)

    def weights_as_series(self, elements):
        idx = self._process_names(elements)
        return pd.Series(self.all_weights[idx], index=self.element_names[idx])

    def _process_names(self, elements):
        elements = np.array(elements)

        difference = np.setdiff1d(
            elements, np.concatenate([self.element_names, self._not_an_element])
        )
        if len(difference) > 0:
            self.calculate_weights(difference)

        idx = [
            int(np.where(self.element_names == el)[0][0])
            for el in elements[np.isin(elements, self.element_names)]
        ]

        return idx


oxide_compositions = Oxide_compositions()
element_weights = Element_weights()
