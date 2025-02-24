import pandas as pd
from typing_extensions import Self

from MagmaPandas.configuration import configuration
from MagmaPandas.enums import Unit
from MagmaPandas.Fe_redox.Fe3Fe2_models import Fe3Fe2_models_dict
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.MagmaFrames.melt import Melt
from MagmaPandas.MagmaSeries import MagmaSeries


class Olivine(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with olivine specific methods.
    """

    # @property
    # def _constructor(self):
    #     """This is the key to letting Pandas know how to keep
    #     derivatives of `MagmaBase` the same type as yours.  It should
    #     be enough to return the name of the Class.  However, in
    #     some cases, `__finalize__` is not called and `new attributes` are
    #     not carried over.  We can fix that by constructing a callable
    #     that makes sure to call `__finalize__` every time."""

    #     def _c(*args, weights=None, **kwargs):
    #         if weights is None:
    #             weights = self._weights.copy(deep=True)
    #         return Olivine(*args, weights=weights, **kwargs).__finalize__(self)

    #     return _c

    @property
    def forsterite(self) -> pd.Series:
        """
        Forsterite contents
        """
        if self._units == Unit.WT_PERCENT:
            moles = self.moles()
        else:
            moles = self

        type = self._datatype.value

        Mg = {"oxide": "MgO", "cation": "Mg"}
        Fe = {"oxide": "FeO", "cation": "Fe"}
        self.recalculate(inplace=True)
        return pd.Series(
            moles[Mg[type]] / (moles[Fe[type]] + moles[Mg[type]]),
            name="Fo#",
        )

    @property
    def formula(self) -> Self:
        """
        Mineral formulas normalised to 4 O p.f.u.
        """
        return self.mineral_formula(O=4)

    def calculate_FeMg_Kd(
        self, melt_wtpc: Melt | MagmaSeries, T_K, P_bar=1, **kwargs
    ) -> pd.Series:
        """
        Calculate Fe-Mg exchange coefficients (Kd) based on measured olivine and melt compositions as (Fe2+/Mg)olivine / (Fe2+/Mg)liquid
        """

        if (not self.index.equals(melt_wtpc.index)) & (
            melt_wtpc.isinstance(pd.DataFrame)
        ):
            raise AttributeError("olivine and melt indeces do not match")

        dfO2 = kwargs.get("dfO2", configuration.dfO2)
        Fe3Fe2_model_name = kwargs.get("Fe3Fe2_model", configuration.Fe3Fe2_model)
        Fe3Fe2_model = Fe3Fe2_models_dict[Fe3Fe2_model_name]

        melt_mol_fractions = melt_wtpc.moles()
        olivine_mol_fractions = self.moles()

        Fe3Fe2 = Fe3Fe2_model.calculate_Fe3Fe2(
            mol_fractions=melt_mol_fractions,
            T_K=T_K,
            P_bar=P_bar,
            dfO2=dfO2,
            **kwargs,
        )

        Fe_melt_total = 1 + Fe3Fe2
        Fe2_FeTotal_melt = 1 / Fe_melt_total
        Fe2_melt = melt_mol_fractions["Fe"] * Fe2_FeTotal_melt

        Kd = (olivine_mol_fractions["Fe"] / olivine_mol_fractions["Mg"]) / (
            Fe2_melt / melt_mol_fractions["Mg"]
        )

        return Kd
