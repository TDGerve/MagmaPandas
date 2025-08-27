.. include:: ./substitutions.rst

================
Available models
================

.. list-table:: Models available in MagmaPandas.
   :header-rows: 1
   :widths: 10 10 10 10 10

   * - Parameter
     - Model
     - Equation
     - Error
     - MagmaPandas name
   * - |fO2| buffer
     - 
     - 
     -
     - 
   * - QFM
     - :cite:t:`oneill_quartz-fayalite-iron_1987`
     - method from 
     -
     - QFM
   * - 
     - :cite:t:`holland_enlarged_1990` 
     - :cite:t:`VanGerve2024`
     -
     -
   * - 
     - :cite:t:`holland1998`
     -
     -
     -
   * - 
     - :cite:t:`jennings_simple_2015`
     -
     -
     -
   * - IW
     - :cite:t:`Hirschmann2021` 
     - eqs. 4 \& 17
     -
     - IW
   * -
     -
     -
     -
     -
   * - |Kd|
     - :cite:t:`toplis_thermodynamics_2005`
     - eq. 10
     - 0.02
     - toplis2005
   * - 
     - :cite:t:`blundy_effect_2020` 
     - eq. 8
     - Na\ :sub:`2`\O + K\ :sub:`2`\O < 6 = 0.019
     - blundy2020
   * - 
     -
     -
     - Na\ :sub:`2`\O + K\ :sub:`2`\O ≥ 6 = 0.04
     -
   * - 
     -
     -
     - Na\ :sub:`2`\O + K\ :sub:`2`\O ≥ 9 = 0.063
     -
   * - 
     - :cite:t:`Putirka2016` 
     - eq. 8a
     - 0.04
     - putirka2016_8a
   * - 
     - 
     - eq. 8b
     - 0.04
     - putirka2016_8b
   * - 
     - 
     - eq. 8c
     - 0.04
     - putirka2016_8c
   * - 
     - 
     - eq. 8d
     - 0.04
     - putirka2016_8d
   * - 
     - :cite:t:`Sun2020a`
     - eq. 7
     - 0.03
     - sun2020
   * -
     - :cite:t:`Saper2022`
     - eq. 10
     - 0.0141
     - saper2022
   * - 
     - fixed value
     - user set
     - user set
     - fixed
   * - 
     -
     -
     -
     -
   * - melt |Fe3Fe2|
     - :cite:t:`borisov_ferricferrous_2018`
     - eq. 4
     - fitted
     - borisov2018
   * - 
     - :cite:t:`kress_compressibility_1991`
     - eq. 7
     - fitted
     - kress_carmichael1991
   * - 
     - :cite:t:`Jayasuriya2004` 
     - eq. 6b
     - fitted
     - jayasuriya2004
   * - 
     - :cite:t:`Putirka2016` 
     - eq. 6b
     - fitted
     - putirka2016_6b
   * - 
     - 
     - eq. 6c
     - fitted
     - putirka2016_6c
   * - 
     - :cite:t:`Deng2020` 
     - eq. 3
     - fitted
     - deng2020
   * - 
     - :cite:t:`ONeill2006` 
     - eq. 10
     - fitted
     - oneill2006
   * - 
     - :cite:t:`oneill_oxidation_2018`
     - eq. 9a
     - fitted
     - oneill2016
   * - 
     - :cite:t:`Armstrong2019` 
     - eq. S12
     - fitted
     - armstrong2019
   * - 
     - :cite:t:`Zhang2017`
     - eq. 11
     - fitted
     - zhang2017
   * - 
     - :cite:t:`Hirschmann2022`
     - eq. 21
     - fitted
     - hirschmann2022
   * - 
     - :cite:t:`Sun2024`  
     - eq. 9 
     - fitted
     - sun2024
   * - 
     - fixed value
     - user set
     - user set
     - fixed  
   * - 
     -
     -
     -
     -
   * - Olivine liquidus
     - :cite:t:`Putirka2008a`
     - eq. 13
     - 71°
     - putirka2008_13
   * - temperature
     - 
     - eq. 14
     - 58°
     - putirka2008_14
   * - 
     - 
     - eq. 15
     - 46°
     - putirka2008_15
   * - 
     - 
     - eq. 16
     - 26°
     - putirka2008_16
   * - 
     - :cite:t:`Sun2020a` 
     - eq. 9
     - 20°
     - sun2020
   * - 
     - :cite:t:`Shea2022`
     - eq. 1
     - 13°
     - shea2022
   * - 
     - :cite:t:`Sugawara2000`
     - eq. 3
     - 33°
     - sugawara2000_3
   * - 
     - 
     - eq. 6a
     - 30°
     - sugawara2000_6a
   * - 
     -
     -
     -
     -
   * - |H2O|-|CO2| saturation
     - :cite:t:`allison_mafich_2022` 
     - eqs. 5 \& 6
     - 
     - allison2022
   * - pressure
     - :cite:t:`iacono-marziano_new_2012`
     - eqs. 12 \& 13
     - 
     - iaconomarziano2012
   * - 
     - :cite:t:`shishkina_solubility_2010`
     - eqs. 9 \& 13
     -
     - shishkina2010


Errors are uncertainties on calibrations, as reported in their original publications. For fitted Fe3Fe2 errors, please see the `Fe3Fe2 error documentation <https://magmapandas.readthedocs.io/en/latest/notebooks/Fe3Fe2_errors.html>`_. For |fo2| and volatile solubility models, errors are not available

     

.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | **Parameter**              | **Model**                                | **Equation**                                            | **name**            |
.. +============================+==========================================+=========================================================+=====================+
.. | \fo2 buffer                |                                          |                                                         |                     |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | QFM                        | :cite:t:`oneill_quartz-fayalite-iron_1987` | Calculated according to method from :cite:t:`VanGerve2024`| QFM                 |
.. |                            | :cite:t:`holland_enlarged_1990`            |                                                         |                     |
.. |                            | :cite:t:`holland1998`                      |                                                         |                     |
.. |                            | :cite:t:`jennings_simple_2015`             |                                                         |                     |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | IW                         | :cite:t:`Hirschmann2021`                   | eqs. 4 \& 17                                            | IW                  |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | \kd                        | :cite:t:`toplis_thermodynamics_2005`       | eq. 10                                                  | toplis2005          |
.. |                            | :cite:t:`blundy_effect_2020`               | eq. 8                                                   | blundy2020          |
.. |                            | :cite:t:`Putirka2016`                      | eq. 8a                                                  | putirka2016_8a      |
.. |                            |                                          | eq. 8b                                                  | putirka2016_8b      |
.. |                            |                                          | eq. 8c                                                  | putirka2016_8c      |
.. |                            |                                          | eq. 8d                                                  | putirka2016_8d      |
.. |                            | :cite:t:`Sun2020a`                         | eq. 7                                                   | sun2020             |
.. |                            | :cite:t:`Saper2022`                        | eq. 10                                                  | saper2022           |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | \fe3fe2                    | :cite:t:`borisov_ferricferrous_2018`       | eq. 4                                                   | borisov2018         |
.. |                            | :cite:t:`kress_compressibility_1991`       | eq. 7                                                   | kress_carmichael1991|
.. |                            | :cite:t:`Jayasuriya2004`                   | eq. 12                                                  | jayasuriya2004      |
.. |                            | :cite:t:`Putirka2016`                      | eq. 6b                                                  | putirka2016_6b      |
.. |                            |                                          | eq. 6c                                                  | putirka2016_6c      |
.. |                            | :cite:t:`Deng2020`                         | eq. 3                                                   | deng2020            |
.. |                            | :cite:t:`ONeill2016`                       | eq. 10                                                  | oneill2016          |
.. |                            | :cite:t:`oneill_oxidation_2018`            | eq. 9a                                                  | oneill2018          |
.. |                            | :cite:t:`Armstrong2019`                    | eq. S12                                                 | armstrong2019       |
.. |                            | :cite:t:`Zhang2017`                        | eq. 11                                                  | zhang2017           |
.. |                            | :cite:t:`Hirschmann2022`                   | eq. 21                                                  | hirschmann2022      |
.. |                            | :cite:t:`Sun2024`                          | eq. 9                                                   | sun2024             |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | Olivine liquidus           | :cite:t:`Putirka2008a`                     | eq. 13                                                  | putirka2008_13      |
.. | temperature                |                                          | eq. 14                                                  | putirka2008_14      |
.. |                            |                                          | eq. 15                                                  | putirka2008_15      |
.. |                            |                                          | eq. 16                                                  | putirka2008_16      |
.. |                            | :cite:t:`Sun2020a`                         | eq. 5                                                   | sun2020             |
.. |                            | :cite:t:`Shea2022`                         | eq. 1                                                   | shea2022            |
.. |                            | :cite:t:`Sugawara2000`                     | eq. 3                                                   | sugawara2000_3      |
.. |                            |                                          | eq. 6a                                                  | sugawara2000_6a     |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
.. | H2O-CO2 saturation pressure| :cite:t:`allison_mafich_2022`              | eqs. 5 \& 8                                             | allison2022         |
.. |                            | :cite:t:`iacono-marziano_new_2012`         | eqs. 12 \& 13                                           | iaconomarziano2012  |
.. |                            | :cite:t:`shishkina_solubility_2010`        | eqs. 9 \& 13                                            | shishkina2010       |
.. +----------------------------+------------------------------------------+---------------------------------------------------------+---------------------+
