===========
MagmaPandas
===========

MagmaPandas is a `Pandas <https://pandas.pydata.org/>`_ based tool set for geochemical calculations and modelling.
It makes working with geochemical data easier by extending the Pandas DataFrame
with methods for common calculations, including those for:

   * mole and cation conversion,
   * mineral composition by stoichiometry,
   * melt thermometry,
   * melt Fe speciation,
   * \ *f*O\ :sub:`2`,
   * melt volatile (CO\ :sub:`2`\-H\ :sub:`2`\O) saturation pressure
   * mineral-melt element partition coefficients,


MagmaPandas can be combined with `MagmaPEC <https://github.com/TDGerve/MagmaPEC>`_ for post-entrapment crystallisation correction of compositions of olivine hosted melt inclusions.

Documentation
-------------
Code documentation is currently being worked on at `magmapandas.readthedocs.io <https://magmapandas.readthedocs.io>`_


How to cite MagmaPandas
------------------------------
If have used MagmaPandas in your research, please reference  `this <https://doi.org/10.1093/petrology/egae006>`_ paper published in Journal of Petrology. To ensure reproducibility, please also mention the release version of MagmaPandas and reference the specific models you used (see `documentation <https://magmapandas.readthedocs.io/en/latest/code_documentation.html#references>`_).



Installation
------------
MagmaPandas can be installed with pip by running:

.. code-block:: bash

    pip install magmapandas

in a terminal.
If you want to install from a specific git branch or release, use:

.. code-block:: bash

    pip install git+https://github.com/TDGerve/MagmaPandas.git@tag

where *tag* should be repleaced by the release tag or branch name (e.g. *v1.0.0* or *development*)


