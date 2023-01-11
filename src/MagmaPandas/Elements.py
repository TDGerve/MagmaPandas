import numpy as np
import elements as e
import pandas as pd


class Oxide_compositions:
    def __init__(self):

        self.oxide_names = np.array([], dtype=str)
        self.cation_names = np.array([], dtype=str)
        self.cation_amount = np.array([], dtype=int)

    def amount(self, oxides):
        idx = self._process_names(oxides)
        return self.cation_amount[idx]

    def names(self, oxides):
        idx = self._process_names(oxides)
        return self.cation_names[idx]

    def calculate_cations_numbers(self, oxides):
        oxides = np.array(oxides)
        new_oxides = oxides[~np.in1d(oxides, self.oxide_names)]
        if len(new_oxides) == 0:
            return

        self.oxide_names = np.append(self.oxide_names, oxides)
        self.cation_names = np.append(self.cation_names, e.cation_names(oxides))
        self.cation_amount = np.append(
            self.cation_amount, e.cation_numbers(oxides).astype(int).values
        )

    def _process_names(self, oxides):
        oxides = np.array(oxides)

        difference = np.setdiff1d(oxides, self.oxide_names)
        if len(difference) > 0:
            self.calculate_cations_numbers(difference)

        idx = [
            int(np.where(self.oxide_names == o)[0])
            for o in oxides[np.isin(oxides, self.oxide_names)]
        ]

        return idx


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
                # Calculate element/oxide weight
                self.all_weights = np.append(
                    self.all_weights, e.calculate_weight(element)
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
            int(np.where(self.element_names == e)[0])
            for e in elements[np.isin(elements, self.element_names)]
        ]

        return idx


oxide_compositions = Oxide_compositions()
element_weights = Element_weights()
