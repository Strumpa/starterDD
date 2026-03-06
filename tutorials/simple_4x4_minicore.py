from starterDD.DDModel.DonjonModel import CoreModel



# Minimal example to instantiate a core model with 4x4 layout of GE14 assemblies with axial variation in the geometry description. 
# The core description YAML file describes the core layout and the axial variation in the geometry description for each assembly type. The geometry description YAML files pointed at in the core description YAML file describe the 2D lattice geometry for each axial region of the assemblies.
path_to_configs = "../data/BWRProgressionProblems/GE14/inputs/"
GE14_core_description_yaml = "GEOM_4x4_mini_CORE.yaml"
minicore_model = CoreModel(name="GE14_4x4_minicore", path_to_yaml_configs=path_to_configs, core_description_yaml=GE14_core_description_yaml)