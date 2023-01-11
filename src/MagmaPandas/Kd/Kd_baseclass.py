from abc import ABC, abstractmethod


class Kd_model(ABC):
    @abstractmethod
    def calculate_Kd(
        cls, Melt_mol_fractions, forsterite_initial, T_K, P_bar, *args, **kwargs
    ):
        pass
