.. # MagmaPandas

.. Pandas based tool set for geochemical processing and modelling.

.. Includes algorithms for:

.. - post-entrapment crystallisation correction of olivine hosted melt inclusions
.. - calculating melt volatile (H<sub>2</sub>O, CO<sub>2</sub>) saturation pressures according to various solubility models
.. - calculating *f*O2 at given P, T and QFM buffer offset

.. and many more!

MagmaPandas is a `Pandas <https://pandas.pydata.org/>`_ based tool set for geochemical calculations and modelling.
It makes working with geochemical data easier by extending the Pandas DataFrame
with methods for common calculations, including:

   * mole and cation conversion,
   * stoichiometric mineral composition,
   * melt thermometry,
   * melt Fe speciation,
   * melt volatile saturation pressure
   * mineral-melt element partition coefficients,


MagmaPandas includes a model for post-entrapment crystallisation correction of olivine hosted melt inclusions.

Installation
------------
MagmaPandas can be installed with pip by running:

.. code-block:: bash

    pip install git+https://github.com/TDGerve/ramCOH.git

in a terminal.
If you want to install from a specific git branch or release, use:

.. code-block:: bash

    pip install git+https://github.com/TDGerve/ramCOH.git@tag

where *tag* should be repleaced by the release tag or branch name (e.g. *v1.0* or *development*)
