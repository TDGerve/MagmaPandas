import pandas as pd
from typing_extensions import Self

from MagmaPandas.enums import Unit
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame


class Clinopyroxene(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with clinopyroxene specific methods.
    """

    @property
    def formula(self) -> Self:
        """
        Mineral formulas normalised to 6 O p.f.u.
        """
        return self.mineral_formula(O=6)

    @property
    def endmembers(self) -> pd.DataFrame:
        """
        endmember components
        """

        cpx_formula = self.formula

        # TODO add code, check thermobar core.py calculate_clinopyroxene_components

        pass

    @property
    def mg_no(self) -> pd.Series:
        """
        Mg numbers
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
            name="Mg_no",
        )
