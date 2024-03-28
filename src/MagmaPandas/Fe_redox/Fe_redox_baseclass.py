from abc import ABC, abstractmethod

import numpy as np


class Fe3Fe2_model(ABC):
    @abstractmethod
    def calculate_Fe3Fe2(
        cls, Melt_mol_fractions, T_K, fO2, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate melt |Fe3Fe2| ratios

        Parameters
        ----------
        Melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        fO2 :   float, array-like
            Oxygen fugacity

        Returns
        -------
        float, array-like
            melt |Fe3Fe2| ratio
        """
        pass

    @abstractmethod
    def get_error(cls, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate one standard deviation errors on |Fe3Fe2| ratios.

        Returns
        -------
        float, array-like
            |Fe3Fe2| error
        """
        pass

    @abstractmethod
    def get_offset_parameters(cls, n: int, *args, **kwargs) -> float | np.ndarray:
        """
        Randomly sample within a standard normal distribution.
        """
        pass

    @abstractmethod
    def get_offset(cls, offset_parameters, *args, **kwargs) -> float | np.ndarray:
        """
        Calculate random samples of |Fe3Fe2| errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """

        pass
