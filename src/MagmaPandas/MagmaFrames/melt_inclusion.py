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

        # Grab default values
        stepsize = kwargs.get("stepsize", 0.005)  # In molar fraction
        converge = kwargs.get("converge", 0.002)  # In molar fraction
        olivine_stepsize_reduction = kwargs.get("olivine_stepsize_reduction", 4)
        decrease_factor = kwargs.get("decrease_factor", 5)
        temperature_converge = kwargs.get("temperature_converge", 0.1)  # In degrees
        QFM_logshift = kwargs.get("QFM_logshift", 1)
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
        olivine_crystallised = pd.Series(0, index=moles.index)
        # Calculate equilibrium forsterite content
        forsterite_EQ = calc_forsterite_EQ(moles)
        # Select Fe removal or addition
        disequilibrium = ~np.isclose(forsterite_EQ, forsterite_host, atol=converge)
        stepsize = pd.Series(stepsize, index=inclusions.index)
        stepsize.loc[forsterite_EQ < forsterite_host] = -stepsize
        # Fe-Mg exchange vector
        FeMg_exchange = pd.DataFrame(0, index=moles.index, columns=moles.columns)
        FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]

        ##### Main Fe-Mg exhange loop #####
        ###################################
        while sum(disequilibrium) > 0:

            stepsize_loop = stepsize.loc[disequilibrium]
            temperatures_loop = temperatures.loc[disequilibrium]
            P_loop = P_bar.loc[disequilibrium]
            FeMg_loop = FeMg_exchange.loc[disequilibrium, :]
            moles_loop = moles.loc[disequilibrium, :]
            fO2_loop = fO2.loc[disequilibrium]
            Fo_host_loop = forsterite_host.loc[disequilibrium]

            # Exchange Fe and Mg
            moles_loop = moles_loop + FeMg_loop.mul(stepsize_loop, axis=0)            
            # Calculate new liquidus temperature
            temperatures_FeMg = moles_loop.convert_moles_wtPercent.temperature(
                P_bar=P_bar
            )
            # Find inclusions outside the equilibrium temperature range
            temperature_mismatch = ~np.isclose(
                temperatures_FeMg, temperatures_loop, atol=temperature_converge
            )
            # Calculate new equilibrium forsterite content
            forsterite_EQ_new = calc_forsterite_EQ(
                moles_loop,
                T_K=temperatures_loop,
                P_bar=P_loop,
                oxygen_fugacity=fO2_loop,
                Fo_initial=Fo_host_loop,
            )
            # Create olivine equilibrium compositions in oxide mol fractions
            olivine = pd.DataFrame(
                {
                    "MgO": forsterite_EQ_new * 2,
                    "FeO": (1 - forsterite_EQ_new) * 2,
                    "SiO2": 1,
                },
                columns=moles_loop.columns,
                index=moles_loop.index,
            ).fillna(0.0)

            # Initialise data for the olivine melting/crystallisation loop
            olivine_stepsize = stepsize_loop.div(olivine_stepsize_reduction)
            olivine_correction = moles_loop.loc[temperature_mismatch, :].copy()
            olivine_correction_loop = olivine_correction.copy()
            T = 1
            ############################# Olivine melting/crystallisation loop ############################
            # Melt or crystallise olivine until the liquidus temperature is back in the equilibrium range #
            while sum(temperature_mismatch) > 1:
                print(f"\rT correction {T:02}", end="")
                # Gather all loop data
                olivine_correction_loop = olivine_correction_loop.loc[
                    temperature_mismatch, :
                ].copy()
                olivine = olivine.loc[temperature_mismatch, :]
                olivine_stepsize = olivine_stepsize.loc[temperature_mismatch]
                # Crystallise or melt olivine
                olivine_correction_loop = olivine_correction_loop + olivine.mul(
                    olivine_stepsize, axis=0
                )
                # Calculate the new liquidus temperature
                temperatures_olivine = (
                    olivine_correction_loop.convert_moles_wtPercent.temperature(
                        P_bar=P_bar.loc[olivine_correction_loop.index]
                    )
                )
                # Keep track of crystllised or melted olivine; negative values for crystallisation
                olivine_crystallised.loc[olivine_stepsize.index] += olivine_stepsize
                # Find inclusions outside the equilibrium temperature range
                temperature_mismatch = ~np.isclose(
                    temperatures_olivine,
                    temperatures_loop.loc[temperatures_olivine.index],
                    atol=temperature_converge,
                )
                # Find overcorrected inclusions
                T_overstepped = ~np.equal(
                    np.sign(
                        temperatures_loop.loc[temperatures_olivine.index] - temperatures_olivine
                    ),
                    np.sign(olivine_stepsize),
                )
                # Reverse one iteration and reduce stepsize for overcorrected inclusions
                decrease_stepsize_T = np.logical_and(
                    T_overstepped, temperature_mismatch
                )
                if sum(decrease_stepsize_T) > 0:
                    idx_1 = olivine_correction_loop.loc[decrease_stepsize_T].index
                    # Reverse meltind/crystallisation
                    olivine_correction_loop.loc[idx_1, :] = olivine_correction_loop.loc[
                        idx_1, :
                    ] - olivine.loc[idx_1, :].mul(olivine_stepsize.loc[idx_1], axis=0)
                    olivine_crystallised.loc[idx_1] -= olivine_stepsize.loc[idx_1]
                    # Decrease stepsize
                    olivine_stepsize.loc[idx_1] = olivine_stepsize.loc[idx_1].div(
                        decrease_factor
                    )
                # Save the corrected compositions
                olivine_correction.loc[
                    olivine_correction_loop.index, :
                ] = olivine_correction_loop.copy()
                T += 1
            # Copy the corrected compositions to the main loop
            moles_loop.loc[olivine_correction.index, :] = olivine_correction.copy()
            # Recalculate equilibrium forsterite content
            forsterite_EQ_new = calc_forsterite_EQ(
                moles_loop,
                T_K=temperatures_loop,
                P_bar=P_loop,
                oxygen_fugacity=fO2_loop,
                Fo_initial=Fo_host_loop,
            )
            # Find inclusions outside the equilibrium forsterite range
            disequilibrium_loop = ~np.isclose(
                forsterite_EQ_new, Fo_host_loop, atol=converge
            )
            # Find overcorrected inclusions
            overstepped = ~np.equal(
                np.sign(forsterite_EQ_new - Fo_host_loop), np.sign(stepsize_loop)
            )
            # Reverse one iteration and decrease stepsize for overcorrected inclusions
            decrease_stepsize = np.logical_and(overstepped, disequilibrium_loop)
            if sum(decrease_stepsize) > 0:
                idx_2 = moles_loop.loc[decrease_stepsize].index
                moles_loop.drop(labels=idx_2, axis=0, inplace=True)
                forsterite_EQ_new.drop(labels=idx_2, axis=0, inplace=True)
                stepsize.loc[idx_2] = stepsize_loop.loc[idx_2].div(decrease_factor)

            moles.loc[moles_loop.index, :] = moles_loop.copy()
            forsterite_EQ[forsterite_EQ_new.index] = forsterite_EQ_new.copy()
            disequilibrium = ~np.isclose(forsterite_EQ, forsterite_host, atol=converge)

        # Convert corrected compositions to oxide wt. %
        corrected_compositions = moles.convert_moles_wtPercent

        return corrected_compositions, olivine_crystallised, forsterite_EQ
