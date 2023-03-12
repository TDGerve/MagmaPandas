# from functools import partial
# import pandas as pd
# import numpy as np
# from typing import Union
# from scipy.optimize import root_scalar


# # Progress bar for lengthy calculations
# from alive_progress import alive_bar
# from alive_progress import config_handler

# config_handler.set_global(
#     title_length=17,
#     manual=True,
#     theme="smooth",
#     spinner=None,
#     stats=False,
#     length=30,
#     force_tty=True,
# )

# import elements as e

# from MagmaPandas.Kd.Ol_melt import calculate_FeMg_Kd
# from MagmaPandas.Fe_redox import Fe_redox
# from MagmaPandas.fO2 import fO2_QFM

# from MagmaPandas.MagmaSeries import MagmaSeries
# from MagmaPandas.MagmaFrames import Olivine, Melt

# from MagmaPandas.configuration import configuration
# from MagmaPandas.PEC.PEC_configuration import PEC_configuration


# def Fe_equilibrate(
#     inclusion: Melt,
#     olivine: Union[float, MagmaSeries],
#     P_bar: float,
#     **kwargs,
# ):
#     """
#     Equilibrate a melt inclusion with it's host olivine through Fe-Mg exchange.
#     Isothermal and isobaric.
#     """
#     # Grab model parameters
#     stepsize = kwargs.get(
#         "stepsize", getattr(PEC_configuration, "stepsize_equilibration")
#     )  # In molar fraction, this is the maximum recommende stepsize.
#     Kd_converge = kwargs.get(
#         "Kd_converge", getattr(PEC_configuration, "Kd_converge")
#     )  # Kd converge
#     Fe2_behaviour = getattr(PEC_configuration, "_Fe2_behaviour")
#     dQFM = kwargs.get("dQFM", configuration.dQFM)
#     # Parameters for the while loop
#     olivine_crystallised = np.array([0.0])
#     decrease_factor = getattr(PEC_configuration, "decrease_factor")
#     # Normalise inclusion composition
#     inclusion = inclusion.fillna(0.0)
#     inclusion = inclusion[inclusion.elements].copy()
#     inclusion = inclusion.normalise()
#     # Calculate temperature and fO2
#     temperature = inclusion.melt_temperature(P_bar=P_bar)
#     fO2 = fO2_QFM(dQFM, temperature, P_bar)
#     # Collect configured models
#     Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
#     Kd_model = calculate_FeMg_Kd
#     # Get olivine molar oxide fractions
#     if isinstance(olivine, float):
#         if olivine < 0 or olivine > 1:
#             raise ValueError(
#                 f"olivine host forsterite: {olivine:.3f} number is not between 0 and 1"
#             )
#         forsterite = olivine
#     elif not isinstance(olivine, MagmaSeries):
#         raise TypeError(
#             f"Olivine host should be forsterite number as float, or full composition as MagmaSeries, not {type(olivine)}"
#         )
#     else:
#         olivine = olivine.moles
#         forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])

#     def calculate_Fe2(
#         mol_fractions, Fe2_behaviour=Fe2_behaviour, temperature=temperature, fO2=fO2
#     ):
#         m_fractions = mol_fractions.copy()
#         m_fractions = m_fractions.normalise()
#         if Fe2_behaviour == "closed system":
#             Fe2_FeTotal = 1
#             Fe3Fe2 = mol_fractions.loc["Fe2O3"] * 2 / mol_fractions.loc["FeO"]
#         elif Fe2_behaviour == "buffered":
#             Fe3Fe2 = Fe3Fe2_model(m_fractions, temperature, fO2)
#             Fe2_FeTotal = 1 / (1 + Fe3Fe2)

#         return Fe3Fe2, Fe2_FeTotal

#     # Fix some parameters for Kd calculation
#     calculate_Kd = partial(
#         Kd_model, forsterite_initial=forsterite, T_K=temperature, P_bar=P_bar
#     )
#     # Calculate initial Fe speciation
#     Fe3Fe2, Fe2_FeTotal = calculate_Fe2(inclusion.moles, Fe2_behaviour="buffered")
#     # For non-buffered Fe2+
#     if Fe2_behaviour == "closed system":
#         # Calculate Fe speciation
#         Fe3_FeTotal = 1 - Fe2_FeTotal
#         # Calculate Fe2O3 wt. %
#         Fe2O3 = inclusion["FeO"] * Fe3_FeTotal
#         FeO_mass, Fe2O3_mass = e.compound_weights(["FeO", "Fe2O3"])
#         inclusion["Fe2O3"] = Fe2O3 * (FeO_mass * 2 / Fe2O3_mass)
#         # Recalculate FeO
#         inclusion["FeO"] = inclusion["FeO"] * (1 - Fe3_FeTotal)
#         inclusion.recalculate(inplace=True)
#     # Calculate moles
#     mi_moles = Melt(columns=inclusion.elements, units="mol fraction", datatype="oxide")
#     mi_moles.loc[0] = inclusion.moles[inclusion.elements].values
#     mi_moles = mi_moles.normalise()
#     # Equilibrium Kd
#     Kd_equilibrium = calculate_Kd(Melt_mol_fractions=mi_moles.iloc[-1], Fe3Fe2=Fe3Fe2)
#     # Real Kd
#     olivine_MgFe = forsterite / (1 - forsterite)
#     melt_MgFe = mi_moles.loc[0, "MgO"] / (mi_moles.loc[0, "FeO"] * Fe2_FeTotal)
#     Kd_real = melt_MgFe / olivine_MgFe
#     # Fe-Mg exchange vector
#     FeMg_vector = pd.Series(0, index=mi_moles.columns, name="FeMg_exchange")
#     FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
#     # Select Fe removal or addition
#     if Kd_real < Kd_equilibrium:
#         stepsize = -stepsize

#     ##### MAIN LOOP #####
#     #####################
#     while not np.isclose(Kd_real, Kd_equilibrium, atol=Kd_converge, rtol=0):

#         # Exchange Fe-Mg
#         idx = mi_moles.index[-1] + stepsize
#         mi_moles.loc[idx] = (mi_moles.iloc[-1] + FeMg_vector.mul(stepsize)).values
#         # Equilibrium Kd and Fe speciation for new composition
#         Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
#         Kd_equilibrium = calculate_Kd(
#             Melt_mol_fractions=mi_moles.normalise().loc[idx], Fe3Fe2=Fe3Fe2
#         )
#         melt_FeMg = (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal) / mi_moles.loc[idx, "MgO"]
#         # Equilibrium olivine composition in oxide mol fractions
#         Fo_equilibrium = 1 / (1 + Kd_equilibrium * melt_FeMg)
#         olivine = MagmaSeries(
#             {"MgO": Fo_equilibrium * 2, "FeO": (1 - Fo_equilibrium) * 2, "SiO2": 1},
#             index=mi_moles.columns,
#         )
#         olivine = olivine.fillna(0.0).normalise()

#         ######################################################
#         # Add or remove olivine to keep temperature constant #
#         olivine_amount = root_scalar(
#             _root_temperature,
#             args=(mi_moles.loc[idx], olivine, temperature, P_bar),
#             x0=0,
#             x1=0.05,
#         ).root
#         mi_moles.loc[idx] = mi_moles.loc[idx] + olivine.mul(olivine_amount)
#         olivine_crystallised = np.append(
#             olivine_crystallised, olivine_crystallised[-1] + olivine_amount
#         )
#         # mi_moles = mi_moles.normalise()
#         temperature_new = mi_moles.loc[idx].convert_moles_wtPercent.melt_temperature(
#             P_bar=P_bar
#         )
#         ######################################################
#         # New equilibrium Kd and Fe speciation
#         Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
#         Kd_equilibrium = calculate_Kd(
#             Melt_mol_fractions=mi_moles.normalise().loc[idx], Fe3Fe2=Fe3Fe2
#         )
#         # Real Kd
#         melt_MgFe = mi_moles.loc[idx, "MgO"] / (mi_moles.loc[idx, "FeO"] * Fe2_FeTotal)
#         Kd_real = melt_MgFe / olivine_MgFe
#         # Assess equilibrium
#         disequilibrium = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0)
#         overstepped = np.sign(Kd_real - Kd_equilibrium) != np.sign(stepsize)
#         decrease_stepsize = np.logical_and(disequilibrium, overstepped)
#         # Reverse one iteration and reduce stepsize if Kd
#         # gets oversteppend by more than the convergence value
#         if decrease_stepsize:
#             mi_moles.drop(index=idx, inplace=True)
#             olivine_crystallised = olivine_crystallised[:-1]
#             # Reset equilibrium and real Kd
#             Fe3Fe2, Fe2_FeTotal = calculate_Fe2(mi_moles.iloc[-1])
#             Kd_equilibrium = calculate_Kd(
#                 Melt_mol_fractions=mi_moles.normalise().iloc[-1], Fe3Fe2=Fe3Fe2
#             )
#             idx = mi_moles.index[-1]
#             melt_MgFe = mi_moles.loc[idx, "MgO"] / (
#                 mi_moles.loc[idx, "FeO"] * Fe2_FeTotal
#             )
#             Kd_real = melt_MgFe / olivine_MgFe
#             stepsize = stepsize / decrease_factor
#     # Recalculate compositions to oxide wt. %
#     equilibrated_composition = mi_moles.convert_moles_wtPercent

#     if len(olivine_crystallised) == 1:
#         olivine_crystallised = np.array([0])
#         temperature_new = temperature

#     equilibrated_composition.index = olivine_crystallised

#     return (
#         equilibrated_composition,
#         olivine_crystallised,
#         {"old": temperature, "new": temperature_new},
#         {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
#     )


# def crystallisation_correction(
#     inclusion: Melt,
#     olivine_host: Union[float, MagmaSeries],
#     FeO_target: Union[int, float, callable],
#     P_bar: float,
#     **kwargs,
# ):

#     """
#     Correct an olivine hosted melt inclusion for post entrapment crystallisation or melting by
#     respectively melting or crystallising host olivine.
#     Expects the melt inclusion is completely equilibrated with the host crystal.
#     The models exits when the user input original melt inclusion FeO content is reached.
#     Loosely based on the postentrapment reequilibration procedure in Petrolog:

#     L. V. Danyushesky and P. Plechov (2011)
#     Petrolog3: Integrated software for modeling crystallization processes
#     Geochemistry, Geophysics, Geosystems, vol 12
#     """
#     equilibration_crystallisation = kwargs.get("eq_crystallised", 0.0)
#     # Grab model parameters
#     stepsize = kwargs.get(
#         "stepsize", getattr(PEC_configuration, "stepsize_crystallisation")
#     )  # In molar fraction, 0.005 is the maximum recommended value
#     converge = kwargs.get(
#         "converge", getattr(PEC_configuration, "FeO_converge")
#     )  # FeO convergence
#     dQFM = kwargs.get("dQFM", configuration.dQFM)
#     calculate_FeO_target = False
#     # Parameters for the while loop
#     decrease_factor = getattr(PEC_configuration, "decrease_factor")
#     # Normalise inclusion composition
#     inclusion = inclusion[inclusion.elements].copy()
#     inclusion = inclusion.fillna(0.0)
#     inclusion = inclusion.normalise()
#     # SET UP INITIAL DATA
#     # Dataframe with new compositions
#     mi_moles = Melt(columns=inclusion.elements, units="mol fraction", datatype="oxide")
#     mi_moles.loc[equilibration_crystallisation] = inclusion.moles[
#         inclusion.elements
#     ].values
#     # Covert to moles left after equilibration
#     mi_moles = mi_moles.mul((1 + equilibration_crystallisation), axis=0)
#     # Get olivine molar oxide fractions
#     if isinstance(olivine_host, float):
#         if olivine_host < 0 or olivine_host > 1:
#             raise ValueError(
#                 f"olivine host forsterite: {olivine_host:.3f} number is not between 0 and 1"
#             )
#         olivine = MagmaSeries(
#             {"MgO": olivine_host * 2, "FeO": (1 - olivine_host) * 2, "SiO2": 1},
#             units="mol fraction",
#             datatype="oxide",
#         )
#         olivine = olivine.normalise().reindex(mi_moles.columns, fill_value=0.0)
#         forsterite = olivine_host
#     elif not isinstance(olivine_host, MagmaSeries):
#         raise TypeError(
#             f"Olivine host should be forsterite number as float, or full composition as MagmaSeries, not {type(olivine_host)}"
#         )
#     else:
#         olivine = olivine_host.moles
#         olivine = olivine.reindex(mi_moles.columns)
#         olivine = olivine.fillna(0.0)
#         olivine.recalculate(inplace=True)
#         forsterite = olivine["MgO"] / (olivine["MgO"] + olivine["FeO"])
#     # Fe-Mg exchange vector
#     FeMg_vector = pd.Series(0, index=mi_moles.columns, name="FeMg_exchange")
#     FeMg_vector.loc[["FeO", "MgO"]] = 1, -1
#     # Inclusion starting FeO
#     FeO = inclusion["FeO"]
#     temperature_old = mi_moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
#         P_bar=P_bar
#     )
#     # Function to calculate Kds
#     calculate_Kd = partial(calculate_Kds, P_bar=P_bar, forsterite=forsterite)

#     if hasattr(FeO_target, "__call__"):
#         calculate_FeO_target = FeO_target
#         FeO_target = calculate_FeO_target(inclusion)

#     if FeO > FeO_target:
#         stepsize = -stepsize

#     ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
#     while not np.isclose(FeO, FeO_target, atol=converge, rtol=0):

#         idx = mi_moles.index[-1] + stepsize
#         mi_moles.loc[idx] = (mi_moles.iloc[-1] + olivine.mul(stepsize)).values
#         # mi_moles = mi_moles.normalise()
#         ###### FE-MG EXCHANGE UNTIL KD EQUILIBRIIUM #####
#         exchange_amount = root_scalar(
#             _root_Kd,
#             args=(mi_moles.loc[idx], FeMg_vector, forsterite, P_bar, {"dQFM": dQFM}),
#             x0=0,
#             x1=0.05,
#             # bracket=[-0.5,0.5],
#         ).root
#         mi_moles.loc[idx] = mi_moles.loc[idx] + FeMg_vector.mul(exchange_amount)
#         # mi_moles = mi_moles.normalise()
#         #################################################
#         mi_wtPercent = mi_moles.convert_moles_wtPercent
#         FeO = mi_wtPercent.loc[idx, "FeO"]

#         if calculate_FeO_target:
#             FeO_target = calculate_FeO_target(mi_wtPercent.loc[idx])

#         disequilibrium = ~np.isclose(FeO, FeO_target, atol=converge, rtol=0)
#         overstepped = np.sign(FeO_target - FeO) != np.sign(stepsize)
#         decrease_stepsize = np.logical_and(disequilibrium, overstepped)
#         # Reverse one iteration and reduce stepsize if FeO
#         # gets oversteppend by more than the convergence value
#         if decrease_stepsize:
#             mi_moles.drop(index=idx, inplace=True)
#             # olivine_corrected.drop(index=idx, inplace=True)
#             mi_wtPercent = mi_moles.convert_moles_wtPercent
#             idx = mi_wtPercent.index[-1]
#             FeO = mi_wtPercent.loc[idx, "FeO"]
#             stepsize = stepsize / decrease_factor

#     Kd_equilibrium, Kd_real = calculate_Kd(mi_moles.iloc[-1])
#     temperature_new = mi_wtPercent.iloc[-1].melt_temperature(P_bar)
#     olivine_corrected = mi_moles.index.values
#     if olivine_corrected.max() == 0:
#         mi_wtPercent = mi_moles.convert_moles_wtPercent

#     return (
#         mi_wtPercent,
#         olivine_corrected,
#         {"equilibrium": Kd_equilibrium, "real": Kd_real},
#         {"old": temperature_old, "new": temperature_new},
#     )


# def _root_temperature(olivine_amount, melt_x_moles, olivine_x_moles, T_K, P_bar):

#     melt_x_new = melt_x_moles + olivine_x_moles.mul(olivine_amount)
#     melt_x_new = melt_x_new.normalise()
#     temperature_new = melt_x_new.convert_moles_wtPercent.melt_temperature(P_bar=P_bar)

#     return T_K - temperature_new


# def _root_Kd(exchange_amount, melt_x_moles, exchange_vector, forsterite, P_bar, kwargs):

#     melt_x_new = melt_x_moles + exchange_vector.mul(exchange_amount)
#     melt_x_new = melt_x_new.normalise()
#     Kd_equilibrium, Kd_real = calculate_Kds(melt_x_new, P_bar, forsterite, **kwargs)

#     return Kd_equilibrium - Kd_real


# def calculate_Kds(melt_x_moles, P_bar, forsterite, **kwargs):

#     Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
#     Kd_model = calculate_FeMg_Kd
#     dQFM = kwargs.get("dQFM", configuration.dQFM)

#     T_K = melt_x_moles.convert_moles_wtPercent.melt_temperature(P_bar)
#     fO2 = fO2_QFM(dQFM, T_K, P_bar)
#     Fe3Fe2 = Fe3Fe2_model(melt_x_moles, T_K, fO2)

#     Kd_observed = calculate_observed_Kd(melt_x_moles, Fe3Fe2, forsterite)

#     Kd_eq = Kd_model(melt_x_moles, forsterite, T_K, Fe3Fe2, P_bar=P_bar)

#     return Kd_eq, Kd_observed


# def calculate_observed_Kd(melt_x_moles, Fe3Fe2, forsterite):

#     melt_x_moles = melt_x_moles.normalise()
#     Fe2_FeTotal = 1 / (1 + Fe3Fe2)
#     melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
#     olivine_MgFe = forsterite / (1 - forsterite)

#     return melt_MgFe / olivine_MgFe


# class PEC_olivine:
#     def __init__(
#         self,
#         inclusions: Melt,
#         olivines: Union[Olivine, pd.Series, float],
#         P_bar: Union[float, int, pd.Series],
#         FeO_target: Union[float, pd.Series, callable],
#         **kwargs,
#     ):

#         self.olivine_corrected = pd.DataFrame(
#             0.0,
#             columns=[
#                 "equilibration_crystallisation",
#                 "PE_crystallisation",
#                 "total_crystallisation",
#             ],
#             index=inclusions.index,
#         )
#         self._FeO_as_function = False

#         # Process attributes
#         ######################

#         # For inclusions
#         if not isinstance(inclusions, Melt):
#             raise TypeError("Inclusions is not a Melt MagmaFrame")
#         else:
#             inclusions = inclusions.fillna(0.0)
#             self.inclusions = inclusions.normalise()
#             self.inclusions_uncorrected = self.inclusions.copy()
#         # For olivine
#         if hasattr(olivines, "index"):
#             try:
#                 olivines = olivines.loc[inclusions.index]
#             except KeyError as err:
#                 print("Inclusion and olivine indeces don't match")
#                 raise err

#         # For olivines
#         if not isinstance(olivines, Olivine):
#             try:
#                 if len(olivines) != self.inclusions.shape[0]:
#                     raise ValueError("Number of olivines and inclusions does not match")
#             except TypeError:
#                 pass
#             forsterite = pd.Series(
#                 olivines, index=self.inclusions.index, name="forsterite"
#             )
#             if (~forsterite.between(0, 1)).any():
#                 raise ValueError(
#                     "olivine host forsterite contents are not all between 0 and 1"
#                 )
#             olivine = Olivine(
#                 {"MgO": forsterite * 2, "FeO": (1 - forsterite) * 2, "SiO2": 1},
#                 index=self.inclusions.index,
#                 units="mol fraction",
#                 datatype="oxide",
#             )
#             self._olivine = olivine.normalise()
#         else:
#             olivines = olivines.fillna(0.0)
#             self._olivine = olivines.moles
#         self._olivine = self._olivine.reindex(
#             columns=self.inclusions.columns, fill_value=0.0
#         )

#         # For pressure
#         try:
#             if len(P_bar) != self.inclusions.shape[0]:
#                 raise ValueError(
#                     "Number of pressure inputs and melt inclusions does not match"
#                 )
#         except TypeError:
#             pass
#         if hasattr(P_bar, "index"):
#             try:
#                 P_bar = P_bar.loc[inclusions.index]
#             except KeyError as err:
#                 print("Inclusion and P_bar indeces don't match")
#                 raise err
#             self.P_bar = P_bar
#         else:
#             self.P_bar = pd.Series(P_bar, index=self.inclusions.index, name=P_bar)

#         # For FeO
#         try:
#             if len(FeO_target) != len(self.inclusions):
#                 raise ValueError(
#                     "Number of initial FeO inputs and inclusions does not match"
#                 )
#             if hasattr(FeO_target, "index"):
#                 if not FeO_target.equals(self.inclusions.index):
#                     raise ValueError(
#                         "FeO target inputs and inclusion indeces don't match"
#                     )
#                 self.FeO_target = FeO_target
#             else:
#                 self.FeO_target = pd.Series(
#                     FeO_target, index=self.inclusions.index, name="FeO_target"
#                 )
#         except TypeError:
#             if isinstance(FeO_target, (int, float)):
#                 self.FeO_target = pd.Series(
#                     FeO_target, index=self.inclusions.index, name="FeO_target"
#                 )
#             elif hasattr(FeO_target, "__call__"):
#                 self._FeO_as_function = True
#                 self.FeO_function = FeO_target
#                 self.FeO_target = self.FeO_function(self.inclusions)

#     def reset(self):
#         self.inclusions = self.inclusions_uncorrected.copy()
#         self.olivine_corrected = pd.DataFrame(
#             0.0,
#             columns=[
#                 "equilibration_crystallisation",
#                 "PE_crystallisation",
#                 "total_crystallisation",
#             ],
#             index=self.inclusions.index,
#         )
#         if self._FeO_as_function:
#             self.FeO_target = self.FeO_function(self.inclusions)

#     @property
#     def olivine(self):
#         return self._olivine.convert_moles_wtPercent

#     @olivine.setter
#     def olivine(self, value):
#         print("Olivine is read only")

#     @property
#     def equilibrated(self):
#         """
#         Check which inclusions are in equilibrium with their host olivine
#         based on modelled and observed Fe-Mg exchange Kd's.
#         """
#         Kd_converge = getattr(PEC_configuration, "Kd_converge") * 1.5
#         Kd_equilibrium, Kd_real = self.calculate_Kds()

#         return pd.Series(
#             np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0),
#             index=Kd_equilibrium.index,
#             name="equilibrated",
#         )

#     @property
#     def Fe_loss(self):
#         """
#         Check which inclusions have experienced Fe loss
#         """
#         FeO_converge = getattr(PEC_configuration, "FeO_converge")
#         if self._FeO_as_function:
#             FeO_target = self.FeO_function(self.inclusions)
#         else:
#             FeO_target = self.FeO_target
#         return pd.Series(
#             ~np.isclose(
#                 FeO_target,
#                 self.inclusions["FeO"],
#                 atol=FeO_converge,
#                 rtol=0,
#             ),
#             index=self.inclusions.index,
#             name="Fe_loss",
#         )

#     def calculate_Kds(self, **kwargs):
#         """
#         Calculate observed and modelled Kds
#         """
#         Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
#         Kd_model = calculate_FeMg_Kd

#         dQFM = kwargs.get("dQFM", configuration.dQFM)

#         melt = kwargs.get("melt", self.inclusions.moles)
#         pressure = kwargs.get("P_bar", self.P_bar)
#         forsterite = kwargs.get("forsterite", self._olivine.forsterite)

#         T_K = kwargs.get("T_K", melt.convert_moles_wtPercent.temperature(pressure))
#         fO2 = kwargs.get("fO2", fO2_QFM(dQFM, T_K, pressure))
#         Fe3Fe2 = Fe3Fe2_model(melt, T_K, fO2)

#         Fe2_FeTotal = 1 / (1 + Fe3Fe2)
#         melt_MgFe = melt["MgO"] / (melt["FeO"] * Fe2_FeTotal)
#         olivine_MgFe = forsterite / (1 - forsterite)
#         Kd_observed = melt_MgFe / olivine_MgFe
#         Kd_observed.rename("real", inplace=True)
#         Kd_observed.index.name = "name"

#         Kd_equilibrium = Kd_model(melt, forsterite, T_K, Fe3Fe2, P_bar=pressure)
#         if isinstance(Kd_equilibrium, (float, int)):
#             Kd_equilibrium = pd.Series(Kd_equilibrium, index=melt.index)
#         Kd_equilibrium.rename("equilibrium", inplace=True)
#         Kd_equilibrium.index.name = "name"

#         return Kd_equilibrium, Kd_observed

#     def _crystallise_multicore(self, sample):
#         """
#         Doesn't seem to work
#         """

#         name, inclusion, olivine, temperature, pressure = sample

#         inclusion_wtPercent = inclusion.convert_moles_wtPercent
#         temperature_new = inclusion_wtPercent.melt_temperature(pressure)

#         sign = np.sign(temperature - temperature_new)

#         olivine_amount = root_scalar(
#             _root_temperature,
#             args=(
#                 inclusion,
#                 olivine,
#                 temperature,
#                 pressure,
#             ),
#             x0=0,
#             x1=sign * 0.01,
#         ).root

#         return name, olivine_amount

#     def Fe_equilibrate(self, inplace=False, **kwargs):
#         """ """

#         # Get settings
#         stepsize = kwargs.get(
#             "stepsize", getattr(PEC_configuration, "stepsize_equilibration")
#         )
#         stepsize = pd.Series(stepsize, index=self.inclusions.index, name="stepsize")
#         decrease_factor = getattr(PEC_configuration, "decrease_factor")
#         Kd_converge = getattr(PEC_configuration, "Kd_converge")
#         dQFM = configuration.dQFM
#         P_bar = self.P_bar
#         # Get configured models
#         Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
#         # Calculate temperature and fO2
#         temperatures = self.inclusions.temperature(P_bar=P_bar)
#         fO2 = fO2_QFM(dQFM, temperatures, P_bar)
#         # Get olivine forsterite contents
#         forsterite_host = self._olivine.forsterite

#         def calculate_Fe2FeTotal(mol_fractions, T_K, oxygen_fugacity):

#             Fe3Fe2 = Fe3Fe2_model(mol_fractions, T_K, oxygen_fugacity)
#             return 1 / (1 + Fe3Fe2)

#         # Set up initial data
#         mi_moles = self.inclusions.moles
#         Kd_equilibrium, Kd_real = self.calculate_Kds(
#             melt=mi_moles,
#             T_K=temperatures,
#             P_bar=P_bar,
#             fO2=fO2,
#             forsterite=forsterite_host,
#         )
#         # Fe-Mg exchange vectors
#         FeMg_exchange = pd.DataFrame(0, index=mi_moles.index, columns=mi_moles.columns)
#         FeMg_exchange.loc[:, ["FeO", "MgO"]] = [1, -1]
#         # Find disequilibrium inclusions
#         disequilibrium = ~np.isclose(Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0)
#         # Set stepsizes acoording to Kd disequilibrium
#         stepsize.loc[Kd_real < Kd_equilibrium] = -stepsize.loc[
#             Kd_real < Kd_equilibrium
#         ].copy()
#         # Store olivine correction for the current itaration
#         olivine_corrected_loop = pd.Series(
#             0, index=mi_moles.index, name="olivine_corrected"
#         )

#         ##### Main Fe-Mg exhange loop #####
#         ###################################
#         total_inclusions = mi_moles.shape[0]
#         reslice = True

#         with alive_bar(
#             total=total_inclusions,
#             title=f"{'Equilibrating': <13} ...",
#         ) as bar:
#             # , Pool() as pool
#             bar(sum(~disequilibrium) / total_inclusions)
#             while sum(disequilibrium) > 0:

#                 if reslice:
#                     # Only reslice al the dataframes and series if the amount of disequilibrium inclusions has changed
#                     stepsize_loop = stepsize.loc[disequilibrium]
#                     temperatures_loop = temperatures.loc[disequilibrium]
#                     P_loop = P_bar.loc[disequilibrium]
#                     FeMg_loop = FeMg_exchange.loc[disequilibrium, :]
#                     mi_moles_loop = mi_moles.loc[disequilibrium, :].copy()
#                     fO2_loop = fO2.loc[disequilibrium]
#                     Fo_host_loop = forsterite_host.loc[disequilibrium]
#                     Kd_eq_loop = Kd_equilibrium.loc[disequilibrium].copy()
#                     Kd_real_loop = Kd_real.loc[disequilibrium].copy()

#                 # Exchange Fe and Mg
#                 mi_moles_loop = mi_moles_loop + FeMg_loop.mul(stepsize_loop, axis=0)
#                 # Calculate new equilibrium Kd and Fe speciation
#                 Kd_eq_loop, _ = self.calculate_Kds(
#                     melt=mi_moles_loop.normalise(),
#                     T_K=temperatures_loop,
#                     fO2=fO2_loop,
#                     P_bar=P_loop,
#                     forsterite=Fo_host_loop,
#                 )
#                 Fe2_FeTotal = calculate_Fe2FeTotal(
#                     mi_moles_loop.normalise(), temperatures_loop, fO2_loop
#                 )
#                 melt_FeMg = (mi_moles_loop["FeO"] * Fe2_FeTotal) / mi_moles_loop["MgO"]
#                 # equilibrium olivine composition in oxide mol fractions
#                 Fo_equilibrium = 1 / (1 + Kd_eq_loop * melt_FeMg)
#                 olivine = Olivine(
#                     {
#                         "MgO": Fo_equilibrium * 2,
#                         "FeO": (1 - Fo_equilibrium) * 2,
#                         "SiO2": 1,
#                     },
#                     index=mi_moles_loop.index,
#                     columns=mi_moles_loop.columns,
#                     units="mol fraction",
#                     datatype="oxide",
#                 )
#                 olivine = olivine.fillna(0.0).normalise()

#                 # Single core
#                 for sample in mi_moles_loop.index:
#                     sample_wtPercent = mi_moles_loop.loc[sample].convert_moles_wtPercent
#                     temperature_new = sample_wtPercent.melt_temperature(
#                         P_bar=P_loop.loc[sample]
#                     )
#                     sign = np.sign(temperatures_loop.loc[sample] - temperature_new)

#                     olivine_amount = root_scalar(
#                         _root_temperature,
#                         args=(
#                             mi_moles_loop.loc[sample],
#                             olivine.loc[sample],
#                             temperatures_loop.loc[sample],
#                             P_loop.loc[sample],
#                         ),
#                         x0=0,
#                         x1=sign * 0.01,
#                         # bracket=[-0.5,0.5],
#                     ).root
#                     mi_moles_loop.loc[sample] = mi_moles_loop.loc[sample] + olivine.loc[
#                         sample
#                     ].mul(olivine_amount)
#                     # current iteration
#                     olivine_corrected_loop.loc[sample] = olivine_amount
#                     # Running total
#                     self.olivine_corrected.loc[
#                         sample, "equilibration_crystallisation"
#                     ] += olivine_amount
#                 # mi_moles_loop = mi_moles_loop.normalise()
#                 ######################################################################
#                 # Recalculate Kds
#                 Kd_eq_loop, Kd_real_loop = self.calculate_Kds(
#                     melt=mi_moles_loop.normalise(),
#                     T_K=temperatures_loop,
#                     P_bar=P_loop,
#                     fO2=fO2_loop,
#                     forsterite=Fo_host_loop,
#                 )
#                 # Find inclusions outside the equilibrium forsterite range
#                 disequilibrium_loop = ~np.isclose(
#                     Kd_eq_loop, Kd_real_loop, atol=Kd_converge, rtol=0
#                 )
#                 # Find overcorrected inclusions
#                 overstepped = ~np.equal(
#                     np.sign(Kd_real_loop - Kd_eq_loop), np.sign(stepsize_loop)
#                 )
#                 # Reverse one iteration and decrease stepsize for overcorrected inclusions
#                 decrease_stepsize = np.logical_and(overstepped, disequilibrium_loop)
#                 if sum(decrease_stepsize) > 0:
#                     idx_FeMg = mi_moles_loop.index[decrease_stepsize]
#                     mi_moles_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
#                     Kd_eq_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
#                     Kd_real_loop.drop(labels=idx_FeMg, axis=0, inplace=True)
#                     self.olivine_corrected.loc[
#                         idx_FeMg, "equilibration_crystallisation"
#                     ] -= olivine_corrected_loop.loc[idx_FeMg]
#                     stepsize.loc[idx_FeMg] = stepsize_loop.loc[idx_FeMg].div(
#                         decrease_factor
#                     )

#                 mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
#                 Kd_equilibrium[Kd_eq_loop.index] = Kd_eq_loop.copy()
#                 Kd_real[Kd_real_loop.index] = Kd_real_loop.copy()
#                 disequilibrium_new = ~np.isclose(
#                     Kd_equilibrium, Kd_real, atol=Kd_converge, rtol=0
#                 )
#                 reslice = ~np.array_equal(disequilibrium_new, disequilibrium)
#                 disequilibrium = disequilibrium_new
#                 # Make sure that disequilibrium stays list-like
#                 # so that any sliced dataframe remains a dataframe
#                 try:
#                     len(disequilibrium)
#                 except TypeError:
#                     disequilibrium = [disequilibrium]

#                 bar(sum(~disequilibrium) / total_inclusions)

#         corrected_compositions = mi_moles.convert_moles_wtPercent
#         temperatures_new = corrected_compositions.temperature(P_bar=P_bar)

#         self.inclusions = corrected_compositions

#         if not inplace:
#             return (
#                 corrected_compositions,
#                 self.olivine_corrected["equilibration_crystallisation"].copy(),
#                 {"Equilibrium": Kd_equilibrium, "Real": Kd_real},
#             )

#     def correct_olivine(self, inplace=False, **kwargs):
#         """
#         Correct an olivine hosted melt inclusion for post entrapment crystallisation or melting by
#         respectively melting or crystallising host olivine.
#         Expects the melt inclusion is completely equilibrated with the host crystal.
#         The models exits when the user input original melt inclusion FeO content is reached.
#         Loosely based on the postentrapment reequilibration procedure in Petrolog:

#         L. V. Danyushesky and P. Plechov (2011)
#         Petrolog3: Integrated software for modeling crystallization processes
#         Geochemistry, Geophysics, Geosystems, vol 12
#         """

#         if not self.equilibrated.all():
#             raise RuntimeError("Inclusion compositions have not all been equilibrated")

#         # Get settings
#         stepsize = kwargs.get(
#             "stepsize", getattr(PEC_configuration, "stepsize_crystallisation")
#         )
#         stepsize = pd.Series(stepsize, index=self.inclusions.index, name="stepsize")
#         decrease_factor = getattr(PEC_configuration, "decrease_factor")
#         FeO_converge = kwargs.get(
#             "FeO_converge", getattr(PEC_configuration, "FeO_converge")
#         )
#         dQFM = getattr(configuration, "dQFM")
#         P_bar = self.P_bar
#         # Inclusion compositions in oxide mol fractions
#         mi_moles = self.inclusions.moles
#         # Convert to the total amount of moles after equilibration
#         mi_moles = mi_moles.mul(
#             (1 + self.olivine_corrected["equilibration_crystallisation"]), axis=0
#         )
#         # Olivine forsterite and Mg/Fe
#         forsterite = self._olivine.forsterite
#         # Fe-Mg exchange vectors
#         FeMg_vector = pd.Series(0, index=mi_moles.columns, name="FeMg_exchange")
#         FeMg_vector.loc[["FeO", "MgO"]] = [1, -1]
#         # Starting FeO and temperature
#         FeO = self.inclusions["FeO"].copy()
#         FeO_target = self.FeO_target.copy()
#         if self._FeO_as_function:
#             FeO_target = self.FeO_function(self.inclusions)
#         else:
#             FeO_target = self.FeO_target.copy()

#         temperature_old = self.inclusions.temperature(P_bar=P_bar)

#         stepsize.loc[FeO > FeO_target] = -stepsize.loc[FeO > FeO_target]
#         FeO_mismatch = ~np.isclose(FeO, FeO_target, atol=FeO_converge, rtol=0)

#         ##### OLIVINE MELTING/CRYSTALLISATION LOOP #####
#         reslice = True
#         total_inclusions = mi_moles.shape[0]
#         with alive_bar(
#             title=f"{'Correcting':<13} ...",
#             total=total_inclusions,
#         ) as bar:

#             while sum(FeO_mismatch) > 0:
#                 bar(sum(~FeO_mismatch) / total_inclusions)

#                 if reslice:
#                     stepsize_loop = stepsize.loc[FeO_mismatch]
#                     mi_moles_loop = mi_moles.loc[FeO_mismatch, :].copy()
#                     idx_olivine = mi_moles_loop.index
#                     olivine_loop = self._olivine.loc[FeO_mismatch, :]

#                 mi_moles_loop = mi_moles_loop + olivine_loop.mul(stepsize_loop, axis=0)
#                 # mi_moles_loop = mi_moles_loop.normalise()
#                 self.olivine_corrected.loc[
#                     idx_olivine, "PE_crystallisation"
#                 ] += stepsize_loop
#                 #################################################
#                 ##### Exchange Fe-Mg to keep Kd equilibrium #####
#                 for sample in mi_moles_loop.index:
#                     exchange_amount = root_scalar(
#                         _root_Kd,
#                         args=(
#                             mi_moles_loop.loc[sample],
#                             FeMg_vector,
#                             forsterite.loc[sample],
#                             P_bar.loc[sample],
#                             {"dQFM": dQFM},
#                         ),
#                         x0=0,
#                         x1=0.05,
#                     ).root
#                     mi_moles_loop.loc[sample] = mi_moles_loop.loc[
#                         sample
#                     ] + FeMg_vector.mul(exchange_amount)
#                 # mi_moles_loop = mi_moles_loop.normalise()
#                 #################################################
#                 # Recalculate FeO
#                 mi_wtPercent = mi_moles_loop.convert_moles_wtPercent
#                 FeO = mi_wtPercent["FeO"]
#                 if self._FeO_as_function:
#                     FeO_target_loop = self.FeO_function(mi_wtPercent)
#                     self.FeO_target.loc[FeO_mismatch] = FeO_target_loop
#                 else:
#                     FeO_target_loop = self.FeO_target.loc[FeO_mismatch]
#                 # Find FeO mismatched and overcorrected inclusions
#                 FeO_mismatch_loop = ~np.isclose(
#                     FeO, FeO_target_loop, atol=FeO_converge, rtol=0
#                 )
#                 FeO_overstepped = ~np.equal(
#                     np.sign(FeO_target_loop - FeO), np.sign(stepsize_loop)
#                 )
#                 decrease_stepsize_FeO = np.logical_and(
#                     FeO_overstepped, FeO_mismatch_loop
#                 )

#                 if sum(decrease_stepsize_FeO) > 0:
#                     # Reverse one step and decrease stepsize
#                     reverse_FeO = mi_moles_loop.index[decrease_stepsize_FeO]
#                     mi_moles_loop.drop(labels=reverse_FeO, axis=0, inplace=True)
#                     self.olivine_corrected.loc[
#                         reverse_FeO, "PE_crystallisation"
#                     ] -= stepsize_loop.loc[reverse_FeO]
#                     stepsize.loc[reverse_FeO] = stepsize.loc[reverse_FeO].div(
#                         decrease_factor
#                     )
#                 # Copy loop data to the main variables
#                 mi_moles.loc[mi_moles_loop.index, :] = mi_moles_loop.copy()
#                 mi_wtPercent = mi_moles.convert_moles_wtPercent
#                 FeO = mi_wtPercent["FeO"]
#                 if self._FeO_as_function:
#                     FeO_target = self.FeO_function(mi_wtPercent)
#                     self.FeO_target = FeO_target
#                 else:
#                     FeO_target = self.FeO_target
#                 FeO_mismatch_new = ~np.isclose(
#                     FeO, FeO_target, atol=FeO_converge, rtol=0
#                 )
#                 reslice = ~np.array_equal(FeO_mismatch_new, FeO_mismatch)
#                 FeO_mismatch = FeO_mismatch_new
#                 bar(sum(~FeO_mismatch) / total_inclusions)

#         corrected_compositions = mi_moles.convert_moles_wtPercent

#         self.olivine_corrected["total_crystallisation"] = self.olivine_corrected[
#             ["equilibration_crystallisation", "PE_crystallisation"]
#         ].sum(axis=1)
#         self.inclusions = corrected_compositions
#         if not inplace:
#             return (
#                 corrected_compositions.copy(),
#                 self.olivine_corrected.copy(),
#                 corrected_compositions.temperature(P_bar),
#             )

#     def correct(self, **kwargs):

#         self.Fe_equilibrate(inplace=True, **kwargs)
#         corrected_compositions, olivine_corrected, temperatures = self.correct_olivine(
#             inplace=False, **kwargs
#         )
#         return corrected_compositions, olivine_corrected, temperatures

#     def correct_inclusion(self, index, plot=True, **kwargs):

#         if type(index) == int:
#             index = self.inclusions_uncorrected.index[index]

#         inclusion = self.inclusions_uncorrected.loc[index].copy()
#         olivine = self._olivine.loc[index].copy()
#         FeO_target = self.FeO_target.loc[index]
#         P_bar = self.P_bar.loc[index]

#         if self._FeO_as_function:
#             FeO_target = self.FeO_function

#         equilibrated, olivine_equilibrated, *_ = Fe_equilibrate(
#             inclusion, olivine, P_bar, **kwargs
#         )
#         corrected, olivine_corrected, *_ = crystallisation_correction(
#             equilibrated.iloc[-1].copy(),
#             olivine,
#             FeO_target,
#             P_bar,
#             eq_crystallised=olivine_equilibrated[-1],
#             **kwargs,
#         )
#         total_corrected = olivine_corrected[-1]

#         equilibrated["correction"] = "equilibration"
#         corrected["correction"] = "correction"

#         total_inclusion = pd.concat([equilibrated, corrected.iloc[1:]], axis=0)

#         if plot:
#             import matplotlib.pyplot as plt
#             import matplotlib.lines as l
#             from labellines import labelLine

#             set_markers = kwargs.get("markers", True)

#             import geoplot as p

#             fontsize = 14

#             fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=False)

#             colors = p.colors.bella.by_key()["color"]

#             linewidth = 5
#             markersize = 90

#             FeO_color = tuple(np.repeat(0.25, 3))

#             plt.plot(
#                 equilibrated["MgO"],
#                 equilibrated["FeO"],
#                 ["-", ".-"][set_markers],
#                 color=colors[1],
#                 # label="equilibration",
#                 linewidth=linewidth,
#                 alpha=0.7,
#             )
#             plt.plot(
#                 corrected["MgO"],
#                 corrected["FeO"],
#                 ["-", ".-"][set_markers],
#                 color=colors[2],
#                 linewidth=linewidth,
#                 alpha=0.7,
#             )
#             ax.scatter(
#                 equilibrated.loc[equilibrated.index[0], "MgO"],
#                 equilibrated.loc[equilibrated.index[0], "FeO"],
#                 marker="^",
#                 color=colors[1],
#                 edgecolors="k",
#                 s=markersize,
#                 zorder=10,
#                 label="Glass",
#             )
#             ax.scatter(
#                 equilibrated.loc[equilibrated.index[-1], "MgO"],
#                 equilibrated.loc[equilibrated.index[-1], "FeO"],
#                 marker="o",
#                 edgecolors="k",
#                 color=colors[3],
#                 s=markersize,
#                 zorder=10,
#                 label="Equilibrated",
#             )
#             ax.scatter(
#                 corrected.loc[corrected.index[-1], "MgO"],
#                 corrected.loc[corrected.index[-1], "FeO"],
#                 marker="s",
#                 color=colors[2],
#                 edgecolors="k",
#                 s=markersize,
#                 zorder=10,
#                 label="Corrected",
#             )

#             middle = sum(ax.get_xlim()) / 2

#             if self._FeO_as_function:
#                 FeO_inital = self.FeO_function(corrected)
#                 ax.plot(
#                     corrected["MgO"], FeO_inital, "-", color=FeO_color, linestyle="-"
#                 )
#                 FeO_target = sum((min(FeO_inital), max(FeO_inital))) / 2
#             else:
#                 ax.axhline(FeO_target, linestyle="-", color=FeO_color, linewidth=1.5)

#             FeO_line = ax.get_lines()[-1]
#             labelLine(
#                 FeO_line,
#                 x=middle,
#                 label="initial FeO",
#                 size=fontsize * 0.8,
#                 color=FeO_color,
#             )

#             ax.set_ylim(ax.get_ylim()[0], max((FeO_target * 1.03, ax.get_ylim()[1])))

#             ax.set_xlabel("MgO (wt. %)")
#             ax.set_ylabel("FeO$^T$\n(wt. %)", rotation=0, labelpad=30)

#             handles, labels = ax.get_legend_handles_labels()

#             handles = handles + [l.Line2D([0], [0], linewidth=0)]

#             labels = labels + [f"{total_corrected * 100:.2f} mol. %\nPEC correction"]

#             ax.legend(
#                 handles,
#                 labels,
#                 title=index,
#                 prop={"family": "monospace", "size": fontsize / 1.5},
#                 fancybox=False,
#                 facecolor="white",
#             )

#         return total_inclusion
