from abc import ABC, abstractmethod

class Kd_model(ABC):

    @abstractmethod
    def calculate_Kd(self):
        pass