from abc import ABC, abstractmethod


class Fe3Fe2_model(ABC):
    @abstractmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions, T_K, fO2, *args, **kwargs
    ):
        pass

    @abstractmethod
    def get_error(cls, *args, **kwargs):
        pass

    @abstractmethod
    def get_offset_parameters(cls, n: int, *args, **kwargs):
        pass

    @abstractmethod
    def get_offset(cls, melt_composition, offset_parameters, *args, **kwargs):
        pass
