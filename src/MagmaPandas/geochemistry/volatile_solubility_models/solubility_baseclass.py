from abc import ABC, abstractmethod

class Solubility_model(ABC):

    @abstractmethod
    def calculate_saturation(self):
        pass

    @abstractmethod
    def calculate_solubility(self):
        pass