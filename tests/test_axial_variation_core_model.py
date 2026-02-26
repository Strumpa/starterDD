
from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.DonjonModel import CoreModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID

from conftest import GE14_INPUTS_DIR_STR

## Test the DonjonModel with simplified axially varying geometry.
# The GE14 assembly refered to in BWR Progression Problems.
# The "DOM" power dominated zone is modeled form 0.0 to 222.0595 cm, and the "VAN" void dominated zone is modeled from 222.0595 to 347.1291 cm.

def test_axial_variation_single_assembly_core_model():
    path_to_configs = GE14_INPUTS_DIR_STR
    GE14_core_description_yaml = "GEOM_single_assembly_CORE.yaml"
    single_assembly_core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=path_to_configs, core_description_yaml=GE14_core_description_yaml)
    assert single_assembly_core_model.geometry_type == "cartesian"
    assert single_assembly_core_model.reactor_type == "BWR"
    assert single_assembly_core_model.core_2D_layout == [["single_GE14_assembly"]]
    #assert single_assembly_core_model.assembly_axial_layouts == [["single_GE14_assembly"]]
    print(single_assembly_core_model.assembly_axial_layouts)
    assert single_assembly_core_model.assembly_axial_layouts == {'single_GE14_assembly': [{'axial_region': 'DOM', 'axial_bounds': [0.0, 222.0595], 'assembly_geometry_file': 'GEOM_GE14_DOM.yaml'}, {'axial_region': 'VAN', 'axial_bounds': [222.0595, 347.1291], 'assembly_geometry_file': 'GEOM_GE14_VAN.yaml'}]}
    assert len(single_assembly_core_model.assemblies) == 1
    assert (0, 0, "single_GE14_assembly") in single_assembly_core_model.assemblies
    assert single_assembly_core_model.assemblies[(0, 0, "single_GE14_assembly")].slices_2D == ["2D_slice_geometry_for_DOM_region_of_single_GE14_assembly", "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly"]
    assert single_assembly_core_model.assemblies[(0, 0, "single_GE14_assembly")].z_bounds == [0.0, 222.0595, 347.1291]
    assert single_assembly_core_model.assemblies[(0, 0, "single_GE14_assembly")].slice_to_geometry_dict == {"2D_slice_geometry_for_DOM_region_of_single_GE14_assembly": "GEOM_GE14_DOM.yaml", "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly": "GEOM_GE14_VAN.yaml"}

    single_assembly_core_model.createAssemblyModels()
    assert len(single_assembly_core_model.assembly_models) == 2
    # assert number of pins in the DOM region model is consistent with the geometry description pointed at in the core description YAML file
    assert single_assembly_core_model.assembly_models[((0, 0, "single_GE14_assembly"), "2D_slice_geometry_for_DOM_region_of_single_GE14_assembly")].count_number_of_pins() == 92
    # assert number of pins in the VAN region model is consistent
    assert single_assembly_core_model.assembly_models[((0, 0, "single_GE14_assembly"), "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly")].count_number_of_pins() == 78

def test_axial_variation_4x4_minicore():
    path_to_configs = GE14_INPUTS_DIR_STR
    GE14_core_description_yaml = "GEOM_4x4_mini_CORE.yaml"
    minicore_model = CoreModel(name="GE14_4x4_minicore", path_to_yaml_configs=path_to_configs, core_description_yaml=GE14_core_description_yaml)
    assert minicore_model.geometry_type == "cartesian"
    assert minicore_model.reactor_type == "BWR"
    assert minicore_model.core_2D_layout == [["GE14_10%", "GE14_20%", "GE14_20%", "GE14_10%"], ["GE14_20%", "GE14_30%", "GE14_30%", "GE14_20%"], ["GE14_20%", "GE14_30%", "GE14_30%", "GE14_20%"], ["GE14_10%", "GE14_20%", "GE14_20%", "GE14_10%"]]
    assert minicore_model.assembly_axial_layouts['GE14_10%'] == [{'axial_region': 'DOM', 'axial_bounds': [0.0, 222.0595], 'assembly_geometry_file': 'GEOM_GE14_DOM.yaml'}, {'axial_region': 'VAN', 'axial_bounds': [222.0595, 347.1291], 'assembly_geometry_file': 'GEOM_GE14_VAN.yaml'}]
    assert minicore_model.assembly_axial_layouts['GE14_20%'] == [{'axial_region': 'DOM', 'axial_bounds': [0.0, 222.0595], 'assembly_geometry_file': 'GEOM_GE14_DOM.yaml'}, {'axial_region': 'VAN', 'axial_bounds': [222.0595, 347.1291], 'assembly_geometry_file': 'GEOM_GE14_VAN.yaml'}]
    assert minicore_model.assembly_axial_layouts['GE14_30%'] == [{'axial_region': 'DOM', 'axial_bounds': [0.0, 222.0595], 'assembly_geometry_file': 'GEOM_GE14_DOM.yaml'}, {'axial_region': 'VAN', 'axial_bounds': [222.0595, 347.1291], 'assembly_geometry_file': 'GEOM_GE14_VAN.yaml'}]
    assert len(minicore_model.assemblies) == 16
    assert (0, 0, "GE14_10%") in minicore_model.assemblies
    assert minicore_model.assemblies[(0, 0, "GE14_10%")].slices_2D == ["2D_slice_geometry_for_DOM_region_of_GE14_10%", "2D_slice_geometry_for_VAN_region_of_GE14_10%"]
    assert minicore_model.assemblies[(0, 0, "GE14_10%")].z_bounds == [0.0, 222.0595, 347.1291]
    assert minicore_model.assemblies[(0, 0, "GE14_10%")].slice_to_geometry_dict == {"2D_slice_geometry_for_DOM_region_of_GE14_10%": "GEOM_GE14_DOM.yaml", "2D_slice_geometry_for_VAN_region_of_GE14_10%": "GEOM_GE14_VAN.yaml"}

if __name__ == "__main__":
    test_axial_variation_single_assembly_core_model()
    print("Test passed successfully.")
    test_axial_variation_4x4_minicore()
    print("Test passed successfully.")