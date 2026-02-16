### starterDD
## A starter pack for Donjon&Dragon : deterministic neutron transport open-source codes
## Dragon and Donjon are part of the Version5 environment developped at Polytechnique Montréal and hosted on the OECD/NEA gitlab : https://git.oecd-nea.org/dragon/5.1 

## StarterDD is :
a python package aiming to help get a more user friendly experience with using Dragon and Donjon. 
## The main modules are : 
- DDModel : allows users to define an Donjon&Dragon model based on an input yaml file. The implementation is limited to a model for the Dragon lattice code for now [16/02/2026]. The Dragon model handles fuel regions sub-divisions for self-shielding problems, physical property assignment and iotopic library definition in CLE-2000 format.
- GeometryBuilder : uses the GLOW application developped by newcleo, avaialble at https://github.com/newcleo-dev-team/glow for surface element geometry definition of complex lattice geometries.


## To be implemented [16/02/2026] : 
- full "abstract" reactor model to support multi-physics
- meshing parameters in yaml input to create "sectorized" geometries for different SALT: tracking methods.
- DonjonModel to create reactor model and 3D geometry.
- python to CLE-2000 generator : create a full Dragon&Donjon case from the python model and generate equivalent CLE-2000 instructions to run. 
