from typing import List, Type, Union

from sympy import rem
from ..parse.readers import _read_file
from .melt import melt
import pandas as pd
import numpy as np
from ..geochemistry.fO2 import fO2_QFM
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg_vectorised
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

    def Fe_restore(self, FeO_initial: Union[int, float], P_bar: float, **kwargs):
        """
        Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
        Isothermal and isobaric.
        """

        return

    def Fe_equilibrate(self, forsterite_host: float, P_bar: float, **kwargs):
        """
        Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
        Isothermal and isobaric.
        """
        # Data checks
        try:
            if len(forsterite_host) != self.shape[0]:
                raise ValueError(
                    "Host forsterite array length does not match the amount of melt inclusions"
                )
        except TypeError:
            forsterite_host = pd.Series(forsterite_host, index=self.index)
        try:
            if len(P_bar) != self.shape[0]:
                raise ValueError(
                    "Pressure array length does not match the amount of melt inclusions"
                )
        except TypeError:
            P_bar = pd.Series(P_bar, index=self.index)

        # Grab model parameters
        stepsize = kwargs.get("stepsize", -0.001)  # In molar fraction
        converge = kwargs.get("converge", 0.002)  # In molar fraction
        temperature_converge = kwargs.get("temperature_converge", 0.1)  # In degrees
        QFM_logshift = kwargs.get("QFM_logshift", 1)
        # Parameters for the while loop
        olivine_crystallised = pd.Series(0, index=self.index)
        olivine_stepsize_reduction = 2
        decrease_factor = 5
        # Calculate temperatures and fO2
        inclusions = self[self.elements].copy()
        inclusions = inclusions.normalise()
        temperatures = inclusions.temperature(P_bar=P_bar)
        fO2 = fO2_QFM(QFM_logshift, temperatures, P_bar)
        # Collect configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)

        # Function for calculating equilibrium forsterite content
        def calc_forsterite_EQ(
            mol_fractions,
            T_K=temperatures.copy(),
            P_bar=P_bar.copy(),
            oxygen_fugacity=fO2.copy(),
            Fo_initial=forsterite_host.copy(),
        ):

            Fo = Fo_initial.copy()
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
            # FeMg Kd
            Kd = Kd_model(mol_fractions, Fo, T_K, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)
            Fe2Mg = mol_fractions["FeO"] * Fe2_FeTotal / mol_fractions["MgO"]
            # Equilibrium forsterite content
            return 1 / (1 + Kd * Fe2Mg)

        # Set up initial data
        moles = inclusions.moles
        # Calculate equilibrium forsterite content
        forsterite_EQ = calc_forsterite_EQ(moles)
        # Select Fe removal or addition
        disequilibrium = ~np.isclose(forsterite_EQ, forsterite_host, atol=converge)
        stepsize = pd.Series(stepsize, index=inclusions.index)
        stepsize.loc[forsterite_EQ > forsterite_host] = -stepsize
        # Fe-Mg exchange vector
        FeMg_exchange = pd.DataFrame(0, index=moles.index, columns=moles.columns)
        FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]

        moles_new = moles.copy()

        while sum(disequilibrium) > 0:

            stepsize = stepsize.loc[disequilibrium]
            temperatures = temperatures.loc[disequilibrium]
            P_bar = P_bar.loc[disequilibrium]
            FeMg_exchange = FeMg_exchange.loc[disequilibrium]
            moles_new = moles_new.loc[disequilibrium]
            fO2 = fO2.loc[disequilibrium]
            forsterite_host = forsterite_host.loc[disequilibrium]

            moles_new = moles_new + FeMg_exchange.mul(stepsize, axis=0)
            temperatures_new = moles_new.convert_moles_wtPercent.temperature(
                P_bar=P_bar
            )

            forsterite_EQ_new = calc_forsterite_EQ(
                moles_new,
                T_K=temperatures,
                P_bar=P_bar,
                oxygen_fugacity=fO2,
                Fo_initial=forsterite_host,
            )

            olivine = pd.DataFrame(
                {
                    "MgO": forsterite_EQ_new * 2,
                    "FeO": (1 - forsterite_EQ_new) * 2,
                    "SiO2": 1,
                },
                columns=moles_new.columns,
                index=moles_new.index,
            ).fillna(0.0)

            olivine_stepsize = stepsize / olivine_stepsize_reduction
            temperature_mismatch = ~np.isclose(
                temperatures_new, temperatures, atol=temperature_converge
            )
            T = 1
            remove_olivine = moles_new.loc[temperature_mismatch].copy()
            while sum(temperature_mismatch) > 1:
                print(f"\rT correction {T:02}", end="")
                remove_olivine = remove_olivine.loc[temperature_mismatch]
                olivine = olivine.loc[temperature_mismatch]
                olivine_stepsize = olivine_stepsize.loc[temperature_mismatch]

                remove_olivine = remove_olivine + olivine.mul(olivine_stepsize, axis=0)
                temperatures_new = remove_olivine.convert_moles_wtPercent.temperature(
                    P_bar=P_bar.loc[remove_olivine.index]
                )
                olivine_crystallised.loc[olivine_stepsize.index] += olivine_stepsize
                temperature_mismatch = ~ np.isclose(
                    temperatures_new,
                    temperatures.loc[temperatures_new.index],
                    atol=temperature_converge,
                )

                T_overstepped = ~ np.equal(
                    np.sign(
                        temperatures.loc[temperatures_new.index] - temperatures_new
                    ),
                    np.sign(stepsize.loc[temperatures_new.index]),
                )
                decrease_stepsize_T = np.logical_and(
                    T_overstepped, temperature_mismatch
                )

                if sum(decrease_stepsize_T) > 0:
                    remove_olivine.loc[decrease_stepsize_T] = remove_olivine.loc[
                        decrease_stepsize_T
                    ] - olivine.loc[decrease_stepsize_T].mul(
                        olivine_stepsize.loc[decrease_stepsize_T], axis=0
                    )
                    olivine_crystallised.loc[remove_olivine.index] = (
                        olivine_crystallised.loc[remove_olivine.index]
                        - olivine_stepsize.loc[remove_olivine.index]
                    )
                    olivine_stepsize.loc[decrease_stepsize_T] = (
                        olivine_stepsize.loc[decrease_stepsize_T] / decrease_factor
                    )
                T += 1

            moles_new.loc[remove_olivine.index] = remove_olivine.copy()

            forsterite_EQ_new = calc_forsterite_EQ(
                moles_new,
                T_K=temperatures,
                P_bar=P_bar,
                oxygen_fugacity=fO2,
                Fo_initial=forsterite_host,
            )

            disequilibrium_new = ~np.isclose(
                forsterite_EQ_new, forsterite_host, atol=converge
            )
            overstepped = ~np.equal(
                np.sign(forsterite_EQ_new - forsterite_host), np.sign(stepsize)
            )
            decrease_stepsize = np.logical_and(overstepped, disequilibrium_new)

            if sum(decrease_stepsize) > 0:
                moles_new.loc[decrease_stepsize] = moles.loc[decrease_stepsize].copy()
                forsterite_EQ_new.loc[decrease_stepsize] = forsterite_EQ.loc[
                    decrease_stepsize
                ].copy()
                stepsize.loc[decrease_stepsize] = (
                    stepsize.loc[decrease_stepsize] / decrease_factor
                )
                FeMg_exchange.loc[decrease_stepsize] = FeMg_exchange.loc[
                    decrease_stepsize
                ].div(decrease_factor)

            moles.loc[moles_new.index] = moles_new.copy()
            forsterite_EQ.loc[forsterite_EQ_new.index] = forsterite_EQ_new.copy()
            disequilibrium = ~np.isclose(
                forsterite_EQ_new, forsterite_host, atol=converge
            )

        wtPercent = moles.convert_moles_wtPercent
        # Reverse overstepped inclusions and recalculate mismatch

        return forsterite_EQ, forsterite_host, olivine_crystallised
