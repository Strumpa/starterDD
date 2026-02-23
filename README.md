### starterDD
## A starter pack for Donjon&Dragon : deterministic neutron transport open-source codes
## Dragon and Donjon are part of the Version5 environment developped at Polytechnique Montréal and hosted on the OECD/NEA gitlab : https://git.oecd-nea.org/dragon/5.1 

## StarterDD is :
a python package aiming to help get a more user friendly experience with using Dragon and Donjon. 
## The main modules are : 
- DDModel : allows users to define an Donjon&Dragon model based on an input yaml file
    - The Dragon model handles fuel regions sub-divisions for self-shielding problems, physical property assignment and isotopic library definition in CLE-2000 format.
    - The Donjon model allows to build 3D cartesian core goemtries from a 2D core layout and associated axial layout descriptors provided in YAML format.
    - The DragonCalculationScheme model allows to create CalculationSteps and associate a SECTORIZED geometry to the DragonModel's TECHNOLOGICAL geometry based on selected dicretization options.
- GeometryBuilder : uses the GLOW application developped by newcleo, available at https://github.com/newcleo-dev-team/glow for surface element geometry definition of complex lattice geometries.


## To be implemented [23/02/2026] : 
- full "abstract" reactor model to support multi-physics.
- Calls to EDI: and COMPO: modules to assist with post treatment. Equivalent Serpent2 detector creation.
- adapting the meshing of neutronics calculation geometries to assembly geometries with control crosses.
- python to CLE-2000 generator : create a full Dragon&Donjon case from the python model and generate equivalent CLE-2000 instructions to run. 
