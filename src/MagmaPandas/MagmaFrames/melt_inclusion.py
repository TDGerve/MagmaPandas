from typing import List, Union
from ..parse.readers import _read_file
import pandas as pd
import numpy as np
from ..geochemistry.fO2 import fO2_QFM
from ..configuration import configuration
from ..geochemistry.Kd_ol_melt import Kd_FeMg_vectorised
from ..geochemistry.Fe_redox import Fe_redox
from functools import partial
from .olivine import Olivine
from .melt import Melt
from ..configuration import configuration


def read_melt_inclusion(
    filepath: str,
    *args,
    index_col: List[str],
    total_col: str = None,
    keep_columns: List[str] = None,
    **kwargs,
) -> "Melt_inclusion":
    """
    Read melt inclusion compositions in wt. % oxide from a .csv file

    """

    return _read_file(
        filepath=filepath,
        *args,
        phase="Melt_inclusion",
        index_col=index_col,
        total_col=total_col,
        keep_columns=keep_columns,
        units="wt. %",
        Type="oxide",
        **kwargs,
    )


class Melt_inclusion(Melt):
    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, weights=self._weights, **kwargs):
            return Melt_inclusion(*args, weights=weights, **kwargs).__finalize__(self)

        return _c

    def Fe_restore(self, FeO_target: Union[int, float], P_bar: float, **kwargs):
        """
        Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
        Isothermal and isobaric.
        """

        return

    def _Fe_equilibrate(self, forsterite_host: float, P_bar: float, **kwargs):
        """
        Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
        Exchange continues until the observed and modelled Fe-Mg Kd converge within
        the model convergence value. After each step of Fe-Mg exchange olivine is either melted
        or crystallised until the liquidus temperature of the inclusion is back to the
        original temperature.

        Isobaric.
        """
        # Data checks
        #
        try:
            if len(forsterite_host) != self.shape[0]:
                raise ValueError("Amount of olivines and inclusions does not match")
        except TypeError:
            pass
        if hasattr(forsterite_host, "index"):
            if not forsterite_host.index.equals(self.index):
                raise ValueError("Inclusions and olivines indeces don't match")
        else:
            forsterite_host = pd.Series(forsterite_host, index=self.index)

        try:
            if len(P_bar) != self.shape[0]:
                raise ValueError(
                    "Pressure array length does not match the amount of melt inclusions"
                )
        except TypeError:
            pass
        if hasattr(P_bar, "index"):
            if not P_bar.index.equals(self.index):
                raise ValueError("Pressure inputs and olivines indeces don't match")
        else:
            P_bar = pd.Series(P_bar, index=self.index)

        # Grab default values
        stepsize = kwargs.get("stepsize", 0.002)  # In molar fraction
        converge = kwargs.get("converge", 0.001)  # Kd converge
        temperature_converge = kwargs.get("temperature_converge", 0.1)  # In degrees
        QFM_logshift = kwargs.get("QFM_logshift", 1)
        # Parameters for while loop
        olivine_stepsize_reduction = 4
        decrease_factor = 5
        olivine_crystallised = pd.Series(0, index=self.index)
        stepsize = pd.Series(stepsize, index=self.index)
        # Normalise inclusion compositions
        inclusions = self[self.elements].copy()
        inclusions = inclusions.normalise()
        # Calculate temperatures and fO2
        temperatures = inclusions.temperature(P_bar=P_bar)
        fO2 = fO2_QFM(QFM_logshift, temperatures, P_bar)
        # Collect configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)

        # Function for calculating equilibrium forsterite content and Fe2+/Fe(total)
        def calculate_Kd(
            mol_fractions,
            T_K,
            P_bar,
            oxygen_fugacity,
            forsterite_host,
        ):
            forsterite = forsterite_host.copy()
            # melt Fe3+/Fe2+
            Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
            # FeMg Kd
            Kd = Kd_model(mol_fractions, forsterite, T_K, P_bar, Fe3Fe2)
            Fe2_FeTotal = 1 / (1 + Fe3Fe2)

            return Kd, Fe2_FeTotal

        # Set up initial data
        mi_moles = inclusions.moles
        # Calculate equilibrium Kd
        Kd_equilibrium, Fe2_FeTotal = calculate_Kd(
            mi_moles, temperatures, P_bar, fO2, forsterite_host
        )
        # Real Kd
        olivine_MgFe = forsterite_host / (1 - forsterite_host)
        melt_MgFe = mi_moles["MgO"] / (mi_moles["FeO"] * Fe2_FeTotal)
        Kd_real = melt_MgFe / olivine_MgFe
        # Fe-Mg exchange vectors
        FeMg_exchange = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
        FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]
        # Find disequilibrium inclusions
        disequilibrium = ~np.isclose(Kd_real, Kd_equilibrium, atol=converge)
        # Set stepsizes acoording to Kd disequilibrium
        stepsize.loc[Kd_real < Kd_equilibrium] = -stepsize

        ##### Main Fe-Mg exhange loop #####
        ###################################
        # Only reslice al the dataframes and series if the amount of disequilibrium inclusions has changed
        reslice = True
        while sum(disequilibrium) > 0:
            if reslice:
                stepsize_loop = stepsize.loc[disequilibrium]
                temperatures_loop = temperatures.loc[disequilibrium]
                P_loop = P_bar.loc[disequilibrium]
                FeMg_loop = FeMg_exchange.loc[disequilibrium, :]
                mi_moles_loop = mi_moles.loc[disequilibrium, :]
                fO2_loop = fO2.loc[disequilibrium]
                Fo_host_loop = forsterite_host.loc[disequilibrium]
                Kd_eq_loop = Kd_equilibrium.loc[disequilibrium]
                Kd_real_loop = Kd_real.loc[disequilibrium]

            # Exchange Fe and Mg
            mi_moles_loop = mi_moles_loop + FeMg_loop.mul(stepsize_loop, axis=0)
            mi_moles_loop = mi_moles_loop.normalise()
            # Calculate new liquidus temperature
            temperatures_new = mi_moles_loop.convert_moles_wtPercent.temperature(
                P_bar=P_bar
            )

            # Calculate new equilibrium Kd and Fe speciation
            Kd_eq_loop, Fe2_FeTotal = calculate_Kd(
                mi_moles_loop, temperatures_loop, P_loop, fO2_loop, Fo_host_loop
            )
            melt_FeMg = (mi_moles_loop["FeO"] * Fe2_FeTotal) / mi_moles_loop["MgO"]
            # equilibrium olivine composition in oxide mol fractions
            Fo_equilibrium = 1 / (1 + Kd_eq_loop * melt_FeMg)
            olivine = Olivine(
                {
                    "MgO": Fo_equilibrium * 2,
                    "FeO": (1 - Fo_equilibrium) * 2,
                    "SiO2": 1,
                },
                index=mi_moles_loop.index,
                units="mol fraction",
                datatype="oxide",
            )
            olivine = olivine.normalise().reindex(
                columns=mi_moles_loop.columns, fill_value=0.0
            )

            # Find inclusions outside the equilibrium temperature range
            temperature_mismatch = ~np.isclose(
                temperatures_new, temperatures_loop, atol=temperature_converge
            )
            ############################# Olivine melting/crystallisation loop ############################
            # Melt or crystallise olivine until the liquidus temperature is back in the equilibrium range #
            # Initialise data for the olivine melting/crystallisation loop
            olivine_stepsize = stepsize_loop.div(olivine_stepsize_reduction)
            # olivine_correction = mi_moles_loop.loc[temperature_mismatch, :].copy()
            reslice_olivine = True
            T = 1
            while sum(temperature_mismatch) > 0:
                print(f"\rT correction {T:02}", end="")
                # Gather all loop data
                if reslice_olivine:
                    olivine_correction = mi_moles_loop.loc[temperature_mismatch, :]
                    olivine_loop = olivine.loc[temperature_mismatch, :]
                    ol_stepsize_loop = olivine_stepsize.loc[temperature_mismatch]
                    idx_olivine = mi_moles_loop.index[temperature_mismatch]

                # Crystallise or melt olivine
                mi_moles_loop.loc[
                    temperature_mismatch, :
                ] = olivine_correction + olivine_loop.mul(ol_stepsize_loop, axis=0)
                # Keep track of crystllised or melted olivine; negative values for crystallisation
                olivine_crystallised.loc[idx_olivine] += olivine_stepsize
                # Calculate the new liquidus temperature
                temperatures_new = mi_moles_loop.convert_moles_wtPercent.temperature(
                    P_bar=P_loop
                )
                # Find inclusions outside the equilibrium temperature range
                temperature_mismatch_new = ~np.isclose(
                    temperatures_new,
                    temperatures_loop,
                    atol=temperature_converge,
                )
                # Find overcorrected inclusions
                T_overstepped = ~np.equal(
                    np.sign(temperatures_loop - temperatures_new),
                    np.sign(olivine_stepsize),
                )
                # Reverse one iteration and reduce stepsize for overcorrected inclusions
                decrease_stepsize_T = np.logical_and(
                    T_overstepped, temperature_mismatch_new
                )
                if sum(decrease_stepsize_T) > 0:
                    idx_reverse = mi_moles_loop.index[decrease_stepsize_T]
                    # Reverse meltind/crystallisation
                    reverse_olivine = olivine_loop.loc[idx_reverse, :].mul(
                        ol_stepsize_loop.loc[idx_reverse], axis=0
                    )
                    mi_moles_loop.loc[idx_reverse, :] = (
                        mi_moles_loop.loc[idx_reverse, :] - reverse_olivine
                    )
                    olivine_crystallised.loc[idx_reverse] -= olivine_stepsize.loc[
                        idx_reverse
                    ]
                    # Decrease stepsize
                    olivine_stepsize.loc[idx_reverse] = olivine_stepsize.loc[
                        idx_reverse
                    ].div(decrease_factor)

                reslice_olivine = ~np.array_equal(
                    temperature_mismatch_new, temperature_mismatch
                )
                temperature_mismatch = temperature_mismatch_new
                T += 1

            # Renormalise after olivine correction
            mi_moles_loop = mi_moles_loop.normalise()

            # Recalculate equilibrium Kd and Fe speciation
            Kd_eq_loop, Fe2_FeTotal = calculate_Kd(
                mi_moles_loop, temperatures_loop, P_loop, fO2_loop, Fo_host_loop
            )
            # Real Kd
            melt_MgFe = mi_moles_loop["MgO"] / (mi_moles_loop["FeO"] * (Fe2_FeTotal))
            Kd_real_loop = melt_MgFe / olivine_MgFe.loc[disequilibrium]

            # Find inclusions outside the equilibrium forsterite range
            disequilibrium_loop = ~np.isclose(Kd_eq_loop, Kd_real_loop, atol=converge)
            # Find overcorrected inclusions
            overstepped = ~np.equal(
                np.sign(Kd_real_loop - Kd_eq_loop), np.sign(stepsize_loop)
            )
            # Reverse one iteration and decrease stepsize for overcorrected inclusions
            decrease_stepsize = np.logical_and(overstepped, disequilibrium_loop)
            if sum(decrease_stepsize) > 0:
                idx_FeMg = mi_moles_loop.index[decrease_stepsize]
                mi_moles_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                Kd_eq_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                Kd_real_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
                stepsize.loc[idx_FeMg] = stepsize_loop.loc[idx_FeMg].div(
                    decrease_factor
                )

            mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
            Kd_equilibrium[Kd_eq_loop.index] = Kd_eq_loop.copy()
            Kd_real[Kd_real_loop.index] = Kd_real_loop.copy()
            disequilibrium_new = ~np.isclose(Kd_equilibrium, Kd_real, atol=converge)
            reslice = ~np.array_equal(disequilibrium_new, disequilibrium)
            disequilibrium = disequilibrium_new

        # Convert corrected compositions to oxide wt. %
        corrected_compositions = mi_moles.convert_moles_wtPercent

        return (
            corrected_compositions,
            olivine_crystallised,
            {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
        )

    def _crystallisation_correction(self, olivine_host, FeO_target, P_bar, **kwargs):
        """
        Correct olivine hosted melt inclusions for post entrapment crystallisation
        or melting by respectively melting or crystallising host olivine.
        Expects the melt inclusion is completely equilibrated with the host crystal.
        The models exits when the user input original melt inclusion FeO content is reached.
        Loosely based on the postentrapment reequilibration procedure in Petrolog:

        L. V. Danyushesky and P. Plechov (2011)
        Petrolog3: Integrated software for modeling crystallization processes
        Geochemistry, Geophysics, Geosystems, vol 12
        """
        # Model parameters
        FeO_as_function = False
        stepsize = kwargs.get(
            "stepsize", 0.001
        )  # In molar fraction, 0.001 is the maximum recommended value
        stepsize = pd.Series(stepsize, index=self.index)
        converge = kwargs.get("converge", 0.05)  # FeO convergence
        Kd_converge = kwargs.get("Kd_converge", 0.001)  # Kd converge
        QFM_logshift = kwargs.get("QFM_logshift", 1)
        # Parameters for the while loop
        olivine_corrected = pd.Series(0, self.index)
        FeMg_exchange_reduction = 4
        decrease_factor = 5
        # Collect configured models
        Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
        Kd_model = getattr(Kd_FeMg_vectorised, configuration().Kd_model)

        # Data checks
        # For olivine
        if hasattr(olivine_host, "index"):
            if not olivine_host.index.equals(self.index):
                raise ValueError("Inclusions and olivines indeces don't match")
        if not isinstance(olivine_host, Olivine):
            try:
                if len(olivine_host) != self.shape[0]:
                    raise ValueError("Number of olivines and inclusions does not match")
            except TypeError:
                pass
            forsterite = pd.Series(olivine_host, index=self.index)
            if (~forsterite.between(0, 1)).any:
                raise ValueError(
                    "olivine host forsterite contents are not all between 0 and 1"
                )
            olivine = Olivine(
                {"MgO": forsterite * 2, "FeO": (1 - forsterite) * 2, "SiO2": 1},
                index=self.index,
                units="mol fraction",
                datatype="oxide",
            )
            olivine = olivine.normalise()
        else:
            olivine = olivine_host.moles
        olivine = olivine.reindex(columns=self.columns, fill_value=0.0)
        
        # For pressure
        try:
            if len(P_bar) != self.shape[0]:
                raise ValueError(
                    "Number of pressure inputs and melt inclusions does not match"
                )
        except TypeError:
            pass
        if hasattr(P_bar, "index"):
            if not P_bar.index.equals(self.index):
                raise ValueError("Pressure inputs and inclusion indeces don't match")
        else:
            P_bar = pd.Series(P_bar, index=self.index)
        # For FeO
        try:
            if len(FeO_target) != len(self):
                raise ValueError("Number of initial FeO inputs and inclusions does not match")
            elif hasattr(FeO_target, "index"):
                if not FeO_target.equals(self.index):
                    raise ValueError("FeO target inputs and inclusion indeces don't match")
            else:
                FeO_target = pd.Series(FeO_target, index=self.index)
        except TypeError:
            if isinstance(FeO_target, (int, float)):
                FeO_target = pd.Series(FeO_target, index=self.index)
            elif hasattr(FeO_target, "__call__"):
                FeO_as_function = True
                FeO_function = FeO_target
                FeO_target = FeO_function(self)
