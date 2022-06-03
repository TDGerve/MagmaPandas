from typing import List, Type, Union
from ..parse.readers import _read_file
from .melt import melt
import pandas as pd
from ..geochemistry.fO2 import fO2_QFM
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg
from ..geochemistry.Fe_redox import Fe_redox
from functools import partial


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

    @staticmethod
    def Fe_restore(
        inclusion: pd.Series, FeO_initial: Union[int, float], P_bar: float, **kwargs
    ):
        """
        Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
        Isothermal and isobaric.
        """
        # Grab model parameters
        converge = kwargs.get("converge", 0.02)
        temperature_converge = kwargs.get("temperature_converge", 1)
        stepsize = kwargs.get("stepsize", 0.001)
        # Select Fe loss or gain
        if inclusion["FeO"] > FeO_initial:
            stepsize = -stepsize
        olivine_stepsize = stepsize / 4
        forsterite = kwargs.get("forsterite", 0.8)
        QFM_logshift = kwargs.get("QFM_logshift", 1)

        # Collect configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

        # Calculate temperature and fO2
        temperature = inclusion.melt_temperature(P_bar=P_bar)
        fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)

        # Set up initial data
        moles = melt_inclusion(
            columns=inclusion.index, units="mol fraction", datatype="oxide"
        )
        moles.loc[0] = inclusion.moles[inclusion.elements].values
        # Fe-Mg exchange vector
        FeMg_exchange = pd.Series(0, index=moles.columns)
        FeMg_exchange.loc[["FeO", "MgO"]] = stepsize, -stepsize
        # Inclusion starting FeO
        FeO = inclusion["FeO"]

        # Parameters for the while loop
        olivine_melted = 0
        decrease_stepsize = True
        decrease_factor = 10

        while abs(FeO - FeO_initial) > converge:
            # Exchange Fe and Mg
            idx = moles.index[-1] + stepsize
            moles.loc[idx] = (moles.iloc[-1] + FeMg_exchange).values
            # New model temperature after Fe-Mg exchange
            temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
                P_bar=P_bar
            )

            # Calculate equilibrium olivine Fo#
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(moles.iloc[-1], temperature, fO2)
            # FeMg Kd
            Kd = Kd_model(moles.iloc[-1], forsterite, temperature, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            Fe2Mg = moles.loc[idx, "FeO"] * Fe2_FeTotal / moles.loc[idx, "MgO"]
            # Equilibrium forsterite content
            forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
            # Equilibrium olivine composition in oxide mol fraction
            olivine = pd.Series(
                {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1},
                index=moles.columns,
            ).fillna(0.0)

            # Crystallise or melt olivine to remain isothermal
            add_olivine = moles.iloc[-1]
            while abs(temperature_new - temperature) > temperature_converge:
                # crystallise olivine until calculated temperature is back at initial
                add_olivine = add_olivine + olivine * (olivine_stepsize)
                temperature_new = add_olivine.convert_moles_wtPercent.melt_temperature(
                    P_bar=P_bar
                )
                # Record crystallised olivine amount
                olivine_melted += olivine_stepsize
            # Copy olivine corrected composition
            moles.iloc[-1] = add_olivine.values

            # New inclusion FeO
            FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

            # Reduce stepsize (once) when approaching initial FeO
            if (abs(FeO - FeO_initial) < (converge * 5)) and decrease_stepsize:
                moles.drop([idx], inplace=True)
                FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]
                stepsize = stepsize / decrease_factor
                olivine_stepsize = stepsize / 4
                FeMg_exchange = FeMg_exchange.div(decrease_factor)
                decrease_stepsize = False

        # Recalculate compositions to oxide wt. %
        wtPercent = moles.convert_moles_wtPercent

        return wtPercent, temperature, temperature_new, olivine_melted

    @staticmethod
    def Fe_equilibrate(
        inclusion: pd.Series, forsterite_host: float, P_bar: float, **kwargs
    ):
        """
        Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
        Isothermal and isobaric.
        """
        # Grab model parameters
        converge = kwargs.get("converge", 0.002)
        temperature_converge = kwargs.get("temperature_converge", 1)
        QFM_logshift = kwargs.get("QFM_logshift", 1)

        # Calculate temperature and fO2
        temperature = inclusion.melt_temperature(P_bar=P_bar)
        fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)

        # Collect configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

        def calc_forsterite(mol_fractions, T_K, P_bar, oxygen_fugacity, Fo_initial):
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
            # FeMg Kd
            Kd = Kd_model(mol_fractions, Fo_initial, T_K, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            Fe2Mg = mol_fractions["FeO"] * Fe2_FeTotal / mol_fractions["MgO"]
            # Equilibrium forsterite content
            return 1 / (1 + Kd * Fe2Mg)

        # Fix temperature, pressure, fO2 and initial forsterite
        calc_forsterite_EQ = partial(
            calc_forsterite,
            T_K=temperature,
            P_bar=P_bar,
            oxygen_fugacity=fO2,
            Fo_initial=forsterite_host,
        )

        # Equilibrium forsterite content
        forsterite_EQ = calc_forsterite_EQ(inclusion.moles)

        stepsize = kwargs.get("stepsize", 0.001)
        # Model Fe loss or gain
        if forsterite_EQ > forsterite_host:
            stepsize = -stepsize
        olivine_stepsize = stepsize / 10

        # Set up initial data
        moles = melt_inclusion(
            columns=inclusion.index, units="mol fraction", datatype="oxide"
        )
        moles.loc[0] = inclusion.moles[inclusion.elements].values
        # Fe-Mg exchange vector
        FeMg_exchange = pd.Series(0, index=moles.columns)
        FeMg_exchange.loc[["FeO", "MgO"]] = -stepsize, stepsize

        # Parameters for the while loop
        olivine_crystallised = 0
        decrease_stepsize = True
        decrease_factor = 10

        while abs(forsterite_host - forsterite_EQ) > converge:
            # Exchange Fe-Mg
            idx = moles.index[-1] + stepsize
            moles.loc[idx] = (moles.iloc[-1] + FeMg_exchange).values
            # New model temperature after Fe-Mg exchange
            temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
                P_bar=P_bar
            )

            # Equilibrium forsterite content for new composition
            forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])
            # Equilibrium olivine composition in oxide mol fractions
            olivine = pd.Series(
                {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1},
                index=moles.columns,
            ).fillna(0.0)

            # Crystallise or melt olivine to remain isothermal
            remove_olivine = moles.iloc[-1]
            while abs(temperature_new - temperature) > temperature_converge:
                # Melt or crystallise olivine until the temperature is back to original.
                remove_olivine = remove_olivine - olivine * olivine_stepsize
                temperature_new = (
                    remove_olivine.convert_moles_wtPercent.melt_temperature(P_bar=P_bar)
                )
                # Record crystallised olivine amount
                olivine_crystallised += olivine_stepsize
            # Copy olivine corrected composition
            moles.iloc[-1] = remove_olivine.values

            # New equilibrium forsterite content
            forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])

            # Decrease stepsize (once) is Fo# is aproaching host olivine Fo#
            if (
                abs(forsterite_host - forsterite_EQ) < (converge * 5)
            ) and decrease_stepsize:
                moles.drop([idx], inplace=True)
                forsterite_EQ = calc_forsterite_EQ(moles.iloc[-1])
                stepsize = stepsize / decrease_factor
                olivine_stepsize = stepsize / 4
                FeMg_exchange = FeMg_exchange.div(decrease_factor)
                decrease_stepsize = False

        # Recalculate compositions to oxide wt. %
        wtPercent = moles.convert_moles_wtPercent

        return wtPercent, forsterite_host, forsterite_EQ, olivine_crystallised
