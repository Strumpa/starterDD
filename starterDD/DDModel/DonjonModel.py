# DonjonModel classes : used to create a DONJON core model from a collection of DRAGON assembly models.

import yaml
from starterDD.DDModel.DragonModel import CartesianAssemblyModel
from ..MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from .helpers import associate_material_to_rod_ID

class CoreModel:
    """
    Docstring for CoreModel
    """
    def __init__(self, name, path_to_yaml_configs, core_description_yaml):
        self.name = name
        # extract the path 
        self.path_to_configs = path_to_yaml_configs
        self.parse_core_description_yaml(f"{self.path_to_configs}/{core_description_yaml}")

    def parse_core_description_yaml(self, core_description_yaml):
        """
        Parse the core description YAML file to extract the necessary information for building the core model.
        This includes the list of assemblies, their positions in the core, and any other relevant information.
        """
        with open(core_description_yaml, 'r') as file:
            yaml_data = yaml.safe_load(file)
        geometric_data = yaml_data.get("CORE_GEOMETRY", {})
        self.geometry_type = geometric_data.get("GEOMETRY_TYPE", "cartesian") # by default, assume cartesian geometry for the core
        self.reactor_type = geometric_data.get("REACTOR_TYPE", "BWR") # by default, assume BWR type for the core, which can be used to define default
        self.core_2D_layout = geometric_data.get("core_2D_layout", []) # list of lists representing the 2D layout of the core, with each entry corresponding to an assembly ID
        self.assembly_axial_layouts = geometric_data.get("assembly_axial_layouts", []) # list of lists representing the axial layout of the core, with each entry corresponding to an assembly ID and
        # the axial extent of the assembly in the core (if axial variation is present in the core description)
        print(self.core_2D_layout)
        number_of_assemblies = sum(len(row) for row in self.core_2D_layout)
        print(f"CoreModel: Parsed core description YAML for core '{self.name}' with {number_of_assemblies} assemblies.")
        self.assembly_axial_layouts = geometric_data.get("assembly_axial_layouts", {}) # dictionary mapping assembly IDs to their axial layouts, if axial variation is present in the core description
        print(f"CoreModel: Parsed assembly axial layouts for core '{self.name}' with {len(self.assembly_axial_layouts)} assemblies having axial variation.")
        # Create AxiallyExtrudedAssemblyModel instances for each assembly with axial variation in the core description
        self.assemblies = {}
        #for assembly_id, axial_layout in self.assembly_axial_layouts.items():
        count = 0
        for row_idx, row in enumerate(self.core_2D_layout):
            for col_idx, assembly_id in enumerate(row):
                print(f"CoreModel: Processing assembly '{assembly_id}' in row {row_idx} index {col_idx} for core '{self.name}'.")
                axial_layout = self.assembly_axial_layouts.get(assembly_id, [])
                slices_2D = []
                z_bounds = []
                axial_regions = []
                slice_to_geometry_dict = {}
                for slice_num, region in enumerate(axial_layout):
                    axial_region = region.get("axial_region", "unknown")
                    lower_bound, upper_bound = region.get("axial_bounds", [0.0, 0.0])
                    assembly_geometry_file = region.get("assembly_geometry_file", "")
                    # Here we would load the 2D slice geometry from the specified file (e.g., using a function like load_2D_slice_geometry(assembly_geometry_file))
                    # For this example, we will just create a placeholder for the 2D slice geometry
                    slice_2D = f"2D_slice_geometry_for_{axial_region}_region_of_{assembly_id}"
                    slice_to_geometry_dict[slice_2D] = assembly_geometry_file
                    slices_2D.append(slice_2D)
                    axial_regions.append(axial_region)
                    if slice_num == 0:
                        z_bounds.append(lower_bound) # add the lower bound of the first slice to the z bounds list
                    z_bounds.append(upper_bound) # add the upper bound of each slice to the z bounds list, which will define the axial extent of each slice in the assembly
                    # Make sure bounds are not repeated and in ascending order
                    z_bounds_sorted = sorted(set(z_bounds))
                    # test that z_bounds are consistent with z_bounds_sorted
                    if z_bounds != z_bounds_sorted:
                        raise ValueError(f"Z bounds for assembly '{assembly_id}' are not in ascending order. Z bounds: {z_bounds}, sorted unique Z bounds: {z_bounds_sorted}. Please check the axial bounds in the core description YAML file. Regions should be defined in z- to z+ order and in a way that the upper bound of one region corresponds to the lower bound of the next region in the axial layout.")

                count += 1    
                self.assemblies[(col_idx, row_idx, assembly_id)] = AxiallyExtrudedAssemblyModel(assembly_id, slices_2D, z_bounds_sorted)
                self.assemblies[(col_idx, row_idx, assembly_id)].set_slice_to_geometry_mapping(slice_to_geometry_dict)
            print(f"CoreModel: Created AxiallyExtrudedAssemblyModel for assembly '{assembly_id}' with axial regions {axial_regions} and axial bounds {z_bounds_sorted}.")
        print(f"CoreModel: Created AxiallyExtrudedAssemblyModel instances for core '{self.name}' with {len(self.assemblies)} assemblies having axial variation.")
        print(f"Assemblies in core '{self.name}': {list(self.assemblies.keys())}, count: {count}.")

    def createAssemblyModels(self):
        """
        Create DRAGON assembly models for each assembly in the core based on the geometry descriptions specified in the core description YAML file.
        This function can be used to build the assembly models for each axial region of the assemblies with axial variation in the core description.
        """
        self.assembly_models = {}
        for assembly_id, assembly_model in self.assemblies.items():
            for slice_2D in assembly_model.slices_2D:
                geometry_file = assembly_model.slice_to_geometry_dict.get(slice_2D, "")
                if f"{self.path_to_configs}/material_compositions.yaml":
                    rod_id_to_material = associate_material_to_rod_ID(f"{self.path_to_configs}/material_compositions.yaml", f"{self.path_to_configs}/{geometry_file}")
                    D5_assembly_model = CartesianAssemblyModel(name=slice_2D, tdt_file=None, geometry_description_yaml=f"{self.path_to_configs}/{geometry_file}")
                    D5_assembly_model.set_rod_ID_to_material_mapping(rod_id_to_material)
                else:
                    D5_assembly_model = CartesianAssemblyModel(name=slice_2D, tdt_file=None, geometry_description_yaml=f"{self.path_to_configs}/{geometry_file}")
                D5_assembly_model.analyze_lattice_description(build_pins=True)
                print(f"CoreModel: Created CartesianAssemblyModel for slice '{slice_2D}' of assembly '{assembly_id}' with geometry description from file '{geometry_file}' without material composition information.")
                self.assembly_models[(assembly_id, slice_2D)] = D5_assembly_model
                print(f"CoreModel: after AssemblyModel geometry analysis for slice '{slice_2D}' of assembly '{assembly_id}', the model has {D5_assembly_model.count_number_of_pins()} pins.")

class AxiallyExtrudedAssemblyModel:
    """
    Docstring for AxiallyExtrudedAssemblyModel
    """
    def __init__(self, name, slices_2D, z_bounds):
        self.name = name
        self.slices_2D = slices_2D # list of 2D slices representing the assembly geometry at different axial levels
        self.z_bounds = z_bounds # list of axial bounds corresponding to the slices, defining the axial extent of each slice in the assembly
        print(f"AxiallyExtrudedAssemblyModel: Created axially extruded assembly model '{self.name}' with {len(self.slices_2D)} slices and axial bounds {self.z_bounds}.")

    def set_slice_to_geometry_mapping(self, slice_to_geometry_dict):
        """
        Set the mapping between each 2D slice and its corresponding geometry description, which can be used to build the assembly model for each axial region.
        """
        self.slice_to_geometry_dict = slice_to_geometry_dict
        print(f"AxiallyExtrudedAssemblyModel: Set slice to geometry mapping for assembly '{self.name}' with {len(self.slice_to_geometry_dict)} slices.")