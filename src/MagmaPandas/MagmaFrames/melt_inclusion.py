from typing import List, Union
from ..parse.readers import _read_file
from .melt import melt
import pandas as pd
from ..geochemistry.fO2 import fO2_QFM
from ..configuration import configuration




def read_melt_inclusion(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> "melt_inclusion":
    """
    Read olivine compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="melt_inclusion",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs,
    )


class melt_inclusion(melt):
    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, **kwargs):
            return melt_inclusion(*args, **kwargs).__finalize__(self)

        return _c

    def Fe_loss_sample(
        self,
        sample,
        FeO_initial: Union[int, float],
        Kd: float,
        P_bar: Union[float, int],
        **kwargs,
    ):
        """
        Docstrings

        """
        inclusion = self.loc[[sample], :]
        if inclusion.loc[sample, "FeO"] > FeO_initial:
            raise ValueError("Inclusion FeO is higher than initial FeO")

        stepsize = kwargs.get("stepsize", 0.001)
        forsterite = kwargs.get("forsterite", 0.8)
        QFM_logshift = kwargs.get("QFM_logshift", 1)
        temperature = inclusion.temperature(P_bar=1e3)


        fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)
        Fe3Fe2_model = configuration().Fe3Fe2
        Kd_model = configuration().Kd

        
        moles = inclusion.moles[inclusion.elements]
        moles.index = [0]
        equilibration_step = pd.Series(0, index=moles.columns)
        equilibration_step.loc[["FeO", "MgO"]] = stepsize, -stepsize

        wtPercent = moles.convert_moles_wtPercent
        FeO = wtPercent.loc[0, 'FeO']
        step = 1

        for step in range(4):

            moles.loc[stepsize*step] = moles.iloc[-1] + equilibration_step

            wtPercent = moles.convert_moles_wtPercent
            FeO = wtPercent.loc[stepsize*step, 'FeO']

            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(moles.iloc[-1], temperature, fO2)
            # FeMg ol-melt Kd
            Kd = Kd_model(moles.iloc[-1], forsterite, temperature, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            Fe2Mg = moles.loc[stepsize * step, "FeO"] * Fe2_FeTotal / moles.loc[stepsize * step, "MgO"]
            # Equilibrium forsterite content
            Fo_EQ = 1 / (1 + Kd * Fe2Mg)
            # Equilibrium olivine composition
            olivine = pd.Series({"MgO": Fo_EQ * 2, "FeO": (1 - Fo_EQ) * 2, "SiO2": 1}, index=[0])
            moles.iloc[-1] = moles.iloc[-1] - olivine * (stepsize / 4)

            wtPercent = moles.convert_moles_wtPercent
            FeO = wtPercent.loc[stepsize*step, 'FeO']

            temperature_new = wtPercent.temperature       


            # step += 1

        return wtPercent, temperature_new, temperature
