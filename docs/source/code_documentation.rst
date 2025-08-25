.. include:: ./substitutions.rst

===================
Code documentation
===================


.. contents::
    :depth: 4

Configuration
#############
.. autoclass:: MagmaPandas.configuration.configuration
    :members:

MagmaFrames
###########
.. automodule:: MagmaPandas.MagmaFrames
    :members:

        .. autoclass:: MagmaFrame
            :members:
        .. autoclass:: Melt
            :members:
        .. autoclass:: Olivine
            :members:
        .. autoclass:: Clinopyroxene
            :members:
        .. autoclass:: Plagioclase
            :members:
        .. autoclass:: Magnetite
            :members:
        .. autoclass:: Ilmenite
            :members:

MagmaSeries
###########
.. automodule:: MagmaPandas.MagmaSeries
    :members:

        .. autoclass:: MagmaSeries
            :members:

melt |Fe3Fe2|
#############
.. automodule:: MagmaPandas.Fe_redox.Fe3Fe2_calculate
        
        .. autofunction:: calculate_Fe3Fe2 
.. automodule:: MagmaPandas.Fe_redox.Fe3Fe2_models
    :members:

|fO2|
#####
.. automodule:: MagmaPandas.fO2
    :members: 

        .. autofunction:: MagmaPandas.fO2.QFM.calculate_fO2
        .. autofunction:: MagmaPandas.fO2.IW.calculate_fO2

Olivine-melt Fe-Mg Kd
#####################
.. automodule:: MagmaPandas.Kd.Ol_melt
    :members: 
        
        .. autofunction:: MagmaPandas.Kd.Ol_melt.FeMg.Kd_calculate.calculate_FeMg_Kd
        .. autofunction:: MagmaPandas.Kd.Ol_melt.FeMg.Kd_calculate.observed_FeMg_Kd

    .. automodule:: MagmaPandas.Kd.Ol_melt.FeMg.Kd_models
        :members:

Silicate liquid density
#######################
.. automodule:: MagmaPandas.rheology.density
    
        .. autofunction:: calculate_density
.. automodule:: MagmaPandas.rheology.viscosity
    
        .. autofunction:: calculate_viscosity

Thermometers
#############
.. automodule:: MagmaPandas.thermometers
    :members:

    melt-only
    *********    

    .. automodule:: MagmaPandas.thermometers.melt
        :members:

    Olivine-melt
    ************

    .. automodule:: MagmaPandas.thermometers.ol_melt
        :members:   

|CO2|-|H2O| solubility models
#############################
.. automodule:: MagmaPandas.volatile_solubility
    :members:

    Allison 2022
    ************
    .. automodule:: MagmaPandas.volatile_solubility.volatile_solubility_models.allison2022
        :members:

    Iacono-Marziano 2012
    ********************
    .. automodule:: MagmaPandas.volatile_solubility.volatile_solubility_models.iaconomarziano2012
        :members:
        :exclude-members: NBO_O_calculate

    Shiskina 2014
    *************
    .. automodule:: MagmaPandas.volatile_solubility.volatile_solubility_models.shiskina2014
        :members:

.. _references:

References
##########

melt |Fe3Fe2| models
********************

.. [1] Borisov A., Behrens H., Holtz F. (2018) Ferric/ferrous ratio in silicate melts: a new model for 1 atm data with special emphasis on the effects of melt composition. Contributions to Mineralogy and Petrology. 173(98).
.. [2] Kress V. C., Carmichael I. S. E. (1991) The compressibility of silicate liquids containing |Fe2O3| and the effect of composition, temperature, oxygen fugacity and pressure on their redox states. Contributions to Mineralogy and Petrology. 108.
.. [3] Jayasuriya K. D., O'Neill H. St.C., Berry A. J., Campbell S. J. (2004) A Mössbauer study of the oxidation state of Fe in silicate melts. American mineralogist. Vol. 89, pp. 1597-1609
.. [4] K. Putirka (2016) Rates and styles of planetary cooling on earth, moon, mars and vesta, using new models for oxygen fugacity, ferric-ferrous ratios, olivine-liquid Fe-Mg exchange, and mantle potential temperature. American mineralogist. Vol. 101, pp 819-840
.. [5] Deng J., Du Z., Karki B. B., Ghosh D. B., Lee K. K. M. (2020) A magma ocean origin to divergent redox evolutions of rocky planetary bodies and early atmospheres. Nature communications. Vol. 11
.. [6] O'neill H. St.C., Berry A. J., McCammon C. C., Jayasuriya K. D., Campbell S. J., Foran G. (2006) An experimental determination of the effect of pressure on the Fe3+/ΣFe ratio of an anhydrous silicate melt to 3.0 GPa. American Mineralogist. Vol. 91, pp. 404-412.
.. [7] Oneill H. St.C., Berry A. J., Mallmann G. (2018) The oxidation state of iron in Mid-Ocean Ridge Basaltic (MORB) glasses: Implications for their petrogenesis and oxygen fugacities. Earth and planetary science letters. Vol. 504, pp. 152-162.
.. [8] Armstrong K. Frost D. J., McCammon C. A., Rubie D. C., Ballaran T. B. (2019) Deep magma formation set the oxidation state of earth's mantle. Science. Vol. 365, pp. 903-906.
.. [9] Zhang H. L., Hirschmann M. M., Cottrell E., Wither A. C. (2017) Effect of pressure on Fe3+/ΣFe ratio in a mafic magma and consequences for magma ocean redox gradients. Geochimica et Cosmochimica Acta. pp. 83-103
.. [10] Hirschmann M. M. (2022) Magma oceans, iron and chromium redox, and the origin of comparatively oxidized planetary mantles. Geochimica et Cosmochimica Acta. Vol. 328, pp. 221-241
.. [11] Sun C., Yao L. (2024) Redox equilibria of iron in low- to high-silica melts: A simple model and its applications to C-H-O-S degassing. Earth and planetary science letters. Vol. 638.


NBO parameterisation
********************

.. [12] Mysen, B. O. (1983) The structure of silicate melts. Ann. Rev. Earth Planet. Sci. 11. 75-97

|fO2| models
************

.. [13] O'Neill H. St. C. (1987) Quartz-fayalite-iron and quartz-fayalite-magnetite equilibria and the free energy of formation of fayalite (Fe2SiO4) and magnetite (Fe3O4). American Mineralogist. 72.
.. [14] Holland T. J. B., Powell R. (2011) An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving new equations of state for solids. Journal of metamorphic geology. 16.
.. [15] Holland T. J. B., Powell R. (1990). An enlarged and updated internally consistent thermodynamic dataset with uncertainties and correlations: the system K2O–Na2O–CaO–MgO–MnO–FeO–Fe2O3–Al2O3–TiO2–SiO2–C–H2–O2. Journal of Metamorphic Geology. 8.
.. [16] Holland T. J. B., Powell R. (1998). An internally consistent thermodynamic data set for phases of petrological interest. Journal of Metamorphic Geology. 16.
.. [17] Jennings E. S., Holland T. J. B. (2015). A simple thermodynamic model for melting of peridotite in the system NCFMASOCr. Journal of Petrology. 56.
.. [18] Hirschmann M. M. (2021) Iron-wüstite revisited: A revised calibration accounting for variable stoichiometry and the effects of pressure. Geochimica et Cosmochimica Acta. Vol. 313, pp. 74-84

Kd models
*********

.. [19] Toplis M. J. (2005) The thermodynamics of iron and magnesium partitioning between olivine and liquid: criteria for assessing and predicting equilibrium in natural and experimental systems. Contributions to mineralogy and petrology. 149. 
.. [20] Blundy (2020) Effect of redox on Fe–Mg–Mn exchange between olivine and melt and an oxybarometer for basalts. Contributions to mineralogy and petrology. 175
.. [21] K. Putirka (2016) Rates and styles of planetary cooling on earth, moon, mars and vesta, using new models for oxygen fugacity, ferric-ferrous ratios, olivine-liquid Fe-Mg exchange, and mantle potential temperature. American mineralogist. Vol. 101, pp 819-840
.. [22] Sun, C., Dasgupta, R. (2020) Thermobarometry of CO2-rich, silica-undersaturated melts constrains cratonic lithosphere thinning through time in areas of kimberlitic magmatism. Earth and Planetary Sience Letters. 550

Silicate melt rheology
**********************

.. [23] Iacovino, K, Till, C. B. (2019) DensityX: a program for calculating the densities of magmatic liquids up to 1627C and 30 kbar. Volcanica 2.1
.. [24] Lange R., Carmichael I. S. (1987) Densities of Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-TiO2-SiO2 liquids: New measurements and derived partial molar properties. Geochimica et cosmochimica acta. 51(11)
.. [25] Lange R. (1997) A revised model for the density and thermal expansivity of K2O-Na2O-CaO-MgO-Al2O3-SiO2 liquids from 700 to 1900 K: extension to crustal magmatic temperatures. Contributions to mineralogy and petrology. 130(1)
.. [26] Ochs III F. A., Lange R. (1999) The Density of Hydrous Magmatic Liquids. Science. 283(5406)
.. [27] Giordano, D., Russel, J. & Dingwell, D. (2008) Viscosity of magmatic liquids: a model. Earth and planetary science letters. Vol 271, issue 1-4, pp. 123-135.

Thermometers
************

.. [28] Putirka K. D. (2008) Thermometers and barometers for volcanic systems. Reviews in mineralogy and geochemistry. 69.
.. [29] Putirka K. D. (2007) Ambient and excess mantle temperatures, olivine thermometry, and active vs. passive upwelling. Chemical geology. 
.. [30] Shea, T., Matz, A. K., Mourey, A. J. (2022) Experimental study of Fe–Mg partitioning and zoning during rapid growth of olivine in Hawaiian tholeiites. Contributions to mineralogy and petrology. 177(12)
.. [31] Sugawara T. (2000) Empirical relationships between temperature, pressure, and MgO content in olivine and pyroxene saturated liquid. Journal of Geophysical Research: Solid earth. Vol. 105, issue B4, pp. 8457-8472


|CO2|-|H2O| solubilty models
****************************

.. [32] Allison C. M., Roggensack K., Clarke A. B. (2022) MafiCH: a general model for |H2O|–|CO2| solubility in mafic magmas. Contributions of mineralogy and petrology. 177:40.
.. [33] Iacono-Marziano G., Morizet Y., Le Trong E., Gaillard F. (2012) New experimental data and semi-empirical parameterization of |H2O|–|CO2| solubility in mafic melts. Geochimica et cosmochimica Acta. 97.
.. [34] Shiskina T. A., Botcharnikov R. E., Holtz F., Almeev R. R., Portnyagin M. V. (2010) Solubility of |H2O|- and |CO2|-bearing fluids in tholeiitic basalts at pressures up to 500 MPa. Chemical Geology. 277.









