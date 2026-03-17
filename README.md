### starterDD
## A starter pack for Donjon&Dragon : deterministic neutron transport open-source codes
## Dragon and Donjon are part of the Version5 environment developed at Polytechnique Montréal and hosted on the OECD/NEA gitlab : https://git.oecd-nea.org/dragon/5.1 [1]


## StarterDD is :
A python package aiming to provide a user-friendly experience for Dragon and Donjon users. 
## The main modules are : 
- DDModel : allows users to define an Donjon&Dragon model based on an input yaml file
    - The Dragon model handles fuel regions sub-divisions for self-shielding problems, physical property assignment and isotopic library definition in CLE-2000 format.
    - The Donjon model allows to build 3D cartesian core geometries from a 2D core layout and associated axial layout descriptors provided in YAML format.
    - The DragonCalculationScheme model allows to create CalculationSteps and associate a SECTORIZED geometry to the DragonModel's TECHNOLOGICAL geometry based on selected discretization options.
- GeometryBuilder : uses the GLOW application developed by newcleo, available at https://github.com/newcleo-dev-team/glow for surface element geometry definition of complex lattice geometries [2].
    - Support for rounded corner BWR assemblies.
    - Support for control crosses.

- InterfaceToDD : allows users to interface the starterDD functionalities to the Version5 (Donjon&Dragon) environment. As of 16/03/2026 Dragon can be executed directly from starterDD. 
    - Specify a path to a Dragon executable in file or as a "dragon_exec" environment variable,
    - Specify a path to draglibs in zipped or unzipped in file or as a "DRAGLIBS_DIR" envionment variable,
    - Execute the Dragon case and save time stamped results.
    - Only a direct self-shielding + MOC calculation scheme is supported for now, more schemes will be implemented.

# To install starterDD 

- git clone https://github.com/Strumpa/starterDD.git
- cd starterDD
- pip install .

# To build the documentation 
- cd docs
- make html

# To open the documentation : 
open ``_build/html/index.html`` in your browser
or use one of the following commands: 
- python -m http.server --directory _build/html/
- sphinx-autobuild . ./_build/html

# To run tests
pytest

# Examples

Examples are available in the tutorials/ folder
- simple_4x4_minicore : simple CoreModel with 4x4 2D layout consisting of 3 different assembly types, defined by 3 axial layouts.
- simple_AT10_1/2L_pincell_compo : single AT10 fuel pin example for Dragon COMPO generation with 1 or 2 levels schemes. 
- simple_glow_GE14 : simplified GE14 model from BWR Progression Problems [3], call to glow through the assembly geometry creation procedure and automated mix and output definitions in cle2000 format.
- simplified_GE14_model : simplified GE14 assembly model from 

# References
- [1] DRAGON5 and DONJON5, the contribution of École Polytechnique de Montréal to the SALOME platform, A. Hébert (2014). Available at  https://git.oecd-nea.org/dragon/5.1
- [2] GLOW by Davide Manzione and Daniele Tomatis : Open source Geometry Layout Oriented Workflow developed by newcleo. Available at https://github.com/newcleo-dev-team/glow 
- [3] Chase Lawing, Scott Palmtag, and Mehdi Asgari, "BWR Progression Problems," Oak Ridge National Laboratory, ORNL/TM-2020/1792 (Sept 2021). Data available at https://github.com/cdlawing1/BWRProgressionProblemInputs 


## To be implemented [16/03/2026] : 
- full "abstract" reactor model to support multi-physics.
- full reactor database creation with DRAGON to support DONJON full core calculations : need to implement support for burnup calculations.
- Support for colorsets / minicore 2D geometries for Dragon.
- python to CLE-2000 generator : create a full Dragon&Donjon case from the python model and generate equivalent CLE-2000 instructions to run. 
