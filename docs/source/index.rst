==============================
MagmaPandas quickstart
==============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   getting_started
   code_documentation
   support
   license

About
-----
MagmaPandas is a :py:mod:`Pandas <pandas:pandas>` based tool set for geochemical calculations and modelling.
It makes working with geochemical data easier by extending the :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
with methods for common calculations, including those for:

   * mole and cation conversion,
   * stoichiometric mineral composition,
   * melt thermometry,
   * melt Fe speciation,
   * melt volatile (CO\ :sub:`2`\-H\ :sub:`2`\O) saturation pressure
   * mineral-melt element partition coefficients,


MagmaPandas includes a model for post-entrapment crystallisation correction of compositions of olivine hosted melt inclusions.


Installation
------------
MagmaPandas can be installed with pip by running:

.. code-block:: bash

    pip install git+https://github.com/TDGerve/MagmaPandas.git

in a terminal.
If you want to install from a specific git branch or release, use:

.. code-block:: bash

    pip install git+https://github.com/TDGerve/MagmaPandas.git@tag

where *tag* should be repleaced by the release tag or branch name (e.g. *v1.0.0* or *development*)

Indices and tables
------------------


* :ref:`genindex`
* :ref:`search`
