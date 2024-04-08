from abc import ABC, abstractmethod

import numpy as np


class Kd_model(ABC):
    @abstractmethod
    def calculate_Kd(
        cls, Melt_mol_fractions, T_K, P_bar, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate mineral-melt partition coefficients

        Parameters
        ----------
        Melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
            Melt composition in oxide mol fractions
        T_K :   float, array-like
            temperature in Kelvin
        P_bar :   float, array-like
            pressure in bar

        Returns
        -------
        float, array-like
            mineral-melt partition coefficients
        """
        pass

    @abstractmethod
    def get_error(cls, *args, **kwargs):
        """
        Return one standard deviation errors on partition coefficients.

        Returns
        -------
        float, array-like
            partition coefficient errors
        """
        pass

    @abstractmethod
    def get_offset_parameters(cls, n: int, *args, **kwargs):
        """
        Randomly sample a standard normal distribution *n* times.

        n   : int
            sample amount.
        """
        pass

    @abstractmethod
    def get_offset(cls, offset_parameters, *args, **kwargs):
        """
        Calculate random samples of partition coefficient errors

        Parameters
        ----------
        offset_parameters : float, array-like
            random samples of a standard normal distribution.
        """
        pass
