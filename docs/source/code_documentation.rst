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

MagmaSeries
###########
.. automodule:: MagmaPandas.MagmaSeries
    :members:

        .. autoclass:: MagmaSeries
            :members:

melt |Fe3Fe2|
#############
.. automodule:: MagmaPandas.Fe_redox
    :members: 
        
        .. autofunction:: calculate_Fe3Fe2 
.. automodule:: MagmaPandas.Fe_redox.models
    :members:

|fO2|
#####
.. automodule:: MagmaPandas.fO2
    :members: 

        .. autofunction:: fO2_QFM

Olivine-melt Fe-Mg Kd
#####################
.. automodule:: MagmaPandas.Kd.Ol_melt
    :members: 
        
        .. autofunction:: calculate_FeMg_Kd
        .. autofunction:: observed_FeMg_Kd
        .. autoclass:: FeMg_blundy
            :members:
        .. autoclass:: FeMg_Toplis
            :members:

Silicate liquid density
#######################
.. automodule:: MagmaPandas.liquid_density
    :members: 
    
        .. autofunction:: calculate_density

Thermometers
#############
.. automodule:: MagmaPandas.thermometers
    :members:

    melt
    ****    

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
    .. automodule:: MagmaPandas.volatile_solubility.models.Allison2022
        :members:

    Iacono-Marziano 2012
    ********************
    .. automodule:: MagmaPandas.volatile_solubility.models.IaconoMarziano
        :members:
        :exclude-members: NBO_O_calculate

    Shiskina 2014
    *************
    .. automodule:: MagmaPandas.volatile_solubility.models.Shiskina
        :members:

.. _references:

References
##########

melt |Fe3Fe2| models
********************

.. [1] Borisov A., Behrens H., Holtz F. (2018) Ferric/ferrous ratio in silicate melts: a new model for 1 atm data with special emphasis on the effects of melt composition. Contributions to Mineralogy and Petrology. 173(98).
.. [2] Kress V. C., Carmichael I. S. E. (1991) The compressibility of silicate liquids containing |Fe2O3| and the effect of composition, temperature, oxygen fugacity and pressure on their redox states. Contributions to Mineralogy and Petrology. 108.


NBO parameterisation
********************

.. [4] Mysen, B. O. (1983) The structure of silicate melts. Ann. Rev. Earth Planet. Sci. 11. 75-97

|fO2| models
************

.. [5] O'Neill H. St. C. (1987) Quartz-fayalite-iron and quartz-fayalite-magnetite equilibria and the free energy of formation of fayalite (Fe2SiO4) and magnetite (Fe3O4). American Mineralogist. 72.
.. [6] Holland T. J. B., Powell R. (2011) An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving new equations of state for solids. Journal of metamorphic geology. 16.
.. [7] Holland T. J. B., Powell R. (1990). An enlarged and updated internally consistent thermodynamic dataset with uncertainties and correlations: the system K2O–Na2O–CaO–MgO–MnO–FeO–Fe2O3–Al2O3–TiO2–SiO2–C–H2–O2. Journal of Metamorphic Geology. 8.
.. [8] Holland T. J. B., Powell R. (1998). An internally consistent thermodynamic data set for phases of petrological interest. Journal of Metamorphic Geology. 16.
.. [9] Jennings E. S., Holland T. J. B. (2015). A simple thermodynamic model for melting of peridotite in the system NCFMASOCr. Journal of Petrology. 56.

Kd models
*********

.. [10] Toplis M. J. (2005) The thermodynamics of iron and magnesium partitioning between olivine and liquid: criteria for assessing and predicting equilibrium in natural and experimental systems. Contributions to mineralogy and petrology. 149. 
.. [11] Blundy (2020) Effect of redox on Fe–Mg–Mn exchange between olivine and melt and an oxybarometer for basalts. Contributions to mineralogy and petrology. 175


Silicate liquid density
***********************

.. [3] Iacovino, K, Till, C. B. (2019) DensityX: a program for calculating the densities of magmatic liquids up to 1627C and 30 kbar. Volcanica 2.1
.. [12] Lange R., Carmichael I. S. (1987) Densities of Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-TiO2-SiO2 liquids: New measurements and derived partial molar properties. Geochimica et cosmochimica acta. 51(11)
.. [13] Lange R. (1997) A revised model for the density and thermal expansivity of K2O-Na2O-CaO-MgO-Al2O3-SiO2 liquids from 700 to 1900 K: extension to crustal magmatic temperatures. Contributions to mineralogy and petrology. 130(1)
.. [14] Ochs III F. A., Lange R. (1999) The Density of Hydrous Magmatic Liquids. Science. 283(5406)

Thermometers
************

.. [15] Putirka K. D. (2008) Thermometers and barometers for volcanic systems. Reviews in mineralogy and geochemistry. 69.
.. [16] Putirka K. D. (2007) Ambient and excess mantle temperatures, olivine thermometry, and active vs. passive upwelling. Chemical geology. 

|CO2|-|H2O| solubilty models
****************************

.. [17] Allison C. M., Roggensack K., Clarke A. B. (2022) MafiCH: a general model for |H2O|–|CO2| solubility in mafic magmas. Contributions of mineralogy and petrology. 177:40.
.. [18] Iacono-Marziano G., Morizet Y., Le Trong E., Gaillard F. (2012) New experimental data and semi-empirical parameterization of |H2O|–|CO2| solubility in mafic melts. Geochimica et cosmochimica Acta. 97.
.. [19] Shiskina T. A., Botcharnikov R. E., Holtz F., Almeev R. R., Portnyagin M. V. (2010) Solubility of |H2O|- and |CO2|-bearing fluids in tholeiitic basalts at pressures up to 500 MPa. Chemical Geology. 277.

Equations of state
******************
.. [20] Katsura T., Tange Y. (2019) A simple derivation of the Birch-Murnaghan equations of state (EOSs) and comparions with EOSs derived from other definitions of finite strain. Minerals. 9:745

        




