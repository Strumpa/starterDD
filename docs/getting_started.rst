Getting Started
===============

Installation
------------

Install the package in development mode:

.. code-block:: bash

   pip install -e .

To enable GLOW geometry generation features:

.. code-block:: bash

   pip install -e ".[glow]"

Quick Example
-------------

Here is a minimal example of defining a material composition and a material mixture:

.. code-block:: python

   from starterDD.MaterialProperties import Composition, MaterialMixture

   # Define a UOX composition
   uox = Composition("UOX_4.5wt%", {
       "U235": 1.05e-3,
       "U238": 2.17e-2,
       "O16":  4.55e-2,
   })

   # Create a material mixture with an index and temperature
   mix = MaterialMixture(
       material_name="UOX_4.5wt%",
       material_mixture_index=1,
       composition=uox,
       temperature=900.0,
       isdepletable=True,
   )

Building Documentation
----------------------

To rebuild this documentation locally:

.. code-block:: bash

   cd docs/
   make html

Then open ``_build/html/index.html`` in your browser.
