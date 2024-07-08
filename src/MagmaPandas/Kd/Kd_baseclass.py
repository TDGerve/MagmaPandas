from abc import ABC, abstractmethod

import numpy as np


class Kd_model(ABC):
    @abstractmethod
    def calculate_Kd(
        cls, melt_mol_fractions, T_K, P_bar, *args, **kwargs
    ) -> float | np.ndarray:
        """
        Calculate mineral-melt partition coefficients

        Parameters
        ----------
        melt_mol_fractions   :   :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
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

    @classmethod
    def _calculate_Kd_(
        cls, melt_mol_fractions, T_K, P_bar, offset_parameters=0.0, *args, **kwargs
    ):

        Kd = cls.calculate_Kd(
            melt_mol_fractions=melt_mol_fractions, T_K=T_K, P_bar=P_bar, *args, **kwargs
        )

        if offset_parameters == 0.0:
            return Kd

        offset = cls.get_offset(
            melt_composition=melt_mol_fractions,
            offset_parameters=offset_parameters,
        )
        Kd = Kd + offset

        # Make sure no negative values are returned
        try:
            Kd[Kd <= 0.0] = 1e-6
        except TypeError:
            Kd = 1e-6 if Kd <= 0.0 else Kd

        # Make sure arrays of length 1 are returned as floats instead
        try:
            if len(Kd) == 1:
                return np.array(Kd)[0]
            return Kd

        except TypeError:
            return Kd

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
