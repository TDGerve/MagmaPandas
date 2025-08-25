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



