===============
Getting started
===============

Installation
------------
MagmaPandas versions 2.0.0 and above can be installed with pip by running

.. code-block:: bash

    pip install magmapandas

in a terminal.

If you would like to install from a specific git branch or release instead run

.. code-block:: bash

    pip install git+https://github.com/TDGerve/magmapandas.git@tag

where *tag* should be repleaced by the release tag or branch name (e.g. *v1.0.0* or *development*)

FeO and Fe2O3 in MagmaPandas
----------------------------

.. important::

    In MagmaPandas, if only FeO or Fe2O3 is included in a composition, it represents total Fe. They only represent :math:`\mathrm{Fe}^{2+}O` and :math:`\mathrm{Fe}^{3+}_{2}\mathrm{O}_{3}`, if both FeO and Fe2O3 are included in the composition.