from abc import ABC, abstractmethod

class Solubility_model(ABC):

    @classmethod
    @abstractmethod
    def calculate_saturation(self):
        pass
    
    @classmethod        
    @abstractmethod
    def calculate_solubility(self):
        pass