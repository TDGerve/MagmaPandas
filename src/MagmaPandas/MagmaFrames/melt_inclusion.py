from typing import List, Type, Union
from ..parse.readers import _read_file
from .melt import melt
import pandas as pd
from ..geochemistry.fO2 import fO2_QFM
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg
from ..geochemistry.Fe_redox import Fe_redox


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

        def _c(*args, weights=self._weights, **kwargs):
            return melt_inclusion(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

        return _c

    def Fe_loss_sample(
        self,
        sample,
        FeO_initial: Union[int, float],
        Kd: float,
        P_bar: float,
        **kwargs,
    ):
        """
        Docstrings

        """
        inclusion = self.loc[sample, :]
        if inclusion["FeO"] > FeO_initial:
            raise ValueError("Inclusion FeO is higher than initial FeO")

        stepsize = kwargs.get("stepsize", 0.001)
        olivine_stepsize = stepsize / 4
        forsterite = kwargs.get("forsterite", 0.8)
        QFM_logshift = kwargs.get("QFM_logshift", 1)
        temperature = inclusion.melt_temperature(P_bar=P_bar)

        fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

        moles = melt_inclusion(
            columns=self.elements, units="mol fraction", datatype="oxide"
        )
        moles.loc[0] = inclusion.moles[inclusion.elements].values
        equilibration_step = pd.Series(0, index=moles.columns)
        equilibration_step.loc[["FeO", "MgO"]] = stepsize, -stepsize

        FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

        step = 1
        olivine_crystallised = 0
        decrease_stepsize = True

        while FeO < FeO_initial:

            idx = moles.index[-1] + stepsize

            moles.loc[idx] = (moles.iloc[-1] + equilibration_step).values

            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(moles.iloc[-1], temperature, fO2)
            # FeMg ol-melt Kd
            Kd = Kd_model(moles.iloc[-1], forsterite, temperature, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            Fe2Mg = (
                moles.loc[idx, "FeO"]
                * Fe2_FeTotal
                / moles.loc[idx, "MgO"]
            )
            # Equilibrium forsterite content
            Fo_EQ = 1 / (1 + Kd * Fe2Mg)
            # Equilibrium olivine composition
            olivine = pd.Series(
                {"MgO": Fo_EQ * 2, "FeO": (1 - Fo_EQ) * 2, "SiO2": 1},
                index=moles.columns,
            ).fillna(0.0)

            temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(P_bar=P_bar)
            add_olivine = moles.iloc[-1]
            while temperature_new < temperature:
                add_olivine = add_olivine + olivine * (olivine_stepsize)
                temperature_new = add_olivine.convert_moles_wtPercent.melt_temperature(P_bar=P_bar)
                olivine_crystallised += olivine_stepsize
            moles.iloc[-1] = add_olivine.values

            FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

            if FeO > FeO_initial:
                if decrease_stepsize:
                    moles.drop([idx], inplace=True)
                    FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]
                    stepsize = stepsize / 10
                    equilibration_step.loc[["FeO", "MgO"]] = stepsize, -stepsize
                    decrease_stepsize = False

            step +=1        

        wtPercent = moles.convert_moles_wtPercent

        return wtPercent, temperature, temperature_new, olivine_crystallised
