.. include:: ./substitutions.rst

=========
Changelog
=========

v2.1.0
------
The following has been added:

- Magnetite and ilmenite subclasses for MagmaFrames.
- more melt Fe3Fe2 models
- more olivine-melt Fe-Mg Kd models
- more melt-only thermometers
- Iron-Wüstite fO2 buffer
- a melt viscosity model

The volatile_solubility module has been restructured and renamed for better easy of use. MagmaFrame properties moles and cations are now class methods.

v2.0.0
------
The PEC module has been removed and moved to it's own MagmaPEC package. Kd and |Fe3Fe2| model classes have been updated to accomodate error propagation in MagmaPEC.

v1.0.0
------
Initial release