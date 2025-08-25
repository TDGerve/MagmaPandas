.. include:: ./substitutions.rst

==============================
MagmaPandas quickstart
==============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   getting_started
   notebooks/Fe3Fe2_errors
   examples
   models
   bibliography
   extending
   code_documentation
   changelog
   support
   license

About
-----
MagmaPandas is a :py:mod:`Pandas <pandas:pandas>` based tool set for geochemical calculations and modelling.
It makes working with geochemical data easier by extending the :py:class:`Pandas DataFrame <pandas:pandas.DataFrame>`
with methods for common calculations, including those for:

   * mole and cation conversion,
   * mineral composition by stoichiometry,
   * melt thermometry,
   * melt Fe speciation,
   * melt |CO2|-|H2O| saturation pressure
   * mineral-melt element partition coefficients,

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

Indices and tables
------------------


* :ref:`genindex`
* :ref:`search`
