starterDD Documentation
========================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api/index

Overview
--------

**starterDD** is a package for building Dragon & Donjon reactor physics models with optional GLOW geometry integration.

It provides tools for:

- **MaterialProperties** — Material composition and mixture handling
- **GeometryBuilder** — Geometry building with optional GLOW/SALOME integration
- **DDModel** — Dragon model definition (assembly, fuel pin, water rod models)
- **GeometryAnalysis** — TDT file parsing and geometry analysis
- **InterfaceToDD** — Dragon module call generation (LIB, MAC)

The package works in two modes:

1. **Standalone**: Works without GLOW, uses abstract classes for Dragon model definition
2. **With GLOW**: Full integration with GLOW for SALOME geometry generation and TDT export


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
