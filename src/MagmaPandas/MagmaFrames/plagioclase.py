import pandas as pd
from typing_extensions import Self

from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame


class Plagioclase(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with plagioclase specific methods.
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
    #         return Plagioclase(*args, weights=weights, **kwargs).__finalize__(self)

    #     return _c

    @property
    def anorthite(self) -> pd.Series:
        """
        Anorthite contents.
        """
        cations = self.cations()
        return pd.Series(
            cations["Ca"] * 100 / (cations["Ca"] + cations["Na"]), name="An"
        )

    @property
    def endmembers(self) -> pd.DataFrame:
        """
        endmember componenents
        """

        cations = self.cations()

        anorthite = cations["Ca"] * 100 / (cations["Ca"] + cations["Na"] + cations["K"])
        albite = cations["Na"] * 100 / (cations["Ca"] + cations["Na"] + cations["K"])
        orthoclase = cations["K"] * 100 / (cations["Ca"] + cations["Na"] + cations["K"])

        return pd.DataFrame(
            {"anorthite": anorthite, "albite": albite, "orthoclase": orthoclase}
        )

    @property
    def formula(self) -> Self:
        """
        Mineral formulas normalised to 8 O p.f.u.
        """

        return self.mineral_formula(O=8)
