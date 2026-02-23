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

   from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
   from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
   from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
   from starterDD.DDModel.helpers import associate_material_to_rod_ID
   from starterDD.InterfaceToDD.dragon_module_calls import LIB

   # Goal : Create an assembly model and instantiate the DDModel builder
   # Model instantiation steps :
   flux_tracking_option = "TISO"  # Tracking type for the assembly model
   # import macros : true
   include_macros = True
   path_to_yaml_compositions = "../data/BWRProgressionProblems/GE14/inputs/material_compositions.yaml"
   path_to_yaml_geometry = "../data/BWRProgressionProblems/GE14/inputs/GEOM_GE14_DOM.yaml"
   path_to_tdt = "../../glow_data/tdt_data"
   tdt_file_name = "GE14_DOM_SSH_IC" #produced by glow.


   # Create the assembly model with the given tdt file and lattice description.
   # use CartesianAssemblyModel class method to create the zone_assignment properties.
   ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                path_to_yaml_geometry)

   BWR_assembly = CartesianAssemblyModel(name="BWR_assembly",
                                       tdt_file=path_to_tdt + "/" + tdt_file_name + ".tdt",
                                       geometry_description_yaml=path_to_yaml_geometry)
   BWR_assembly.set_rod_ID_to_material_mapping(ROD_to_material)
   BWR_assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0, coolant_temperature=600.0,moderator_temperature=600.0, structural_temperature=600.0)
   # Parses information from the geometry description YAML file.
   BWR_assembly.analyze_lattice_description(build_pins=True)
   # Parse compositions from the YAML and attach them to the assembly model
   compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
   BWR_assembly.set_material_compositions(compositions)

   # Run the by_pin numbering
   BWR_assembly.number_fuel_material_mixtures_by_pin()

   # Build a call to the DRAGON LIB: module
   BWR_assembly.identify_generating_and_daughter_mixes()
   lib = LIB(BWR_assembly)
   output_path = lib.write_to_c2m("outputs", "BWR_LIB_definition")

Building Documentation
----------------------

To rebuild this documentation locally:

.. code-block:: bash

   cd docs/
   make html

Then open ``_build/html/index.html`` in your browser or use the following command:

.. code-block:: bash

   sphinx-autobuild . ./_build/html
