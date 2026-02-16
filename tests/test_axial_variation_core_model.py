
from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.DonjonModel import CoreModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID

## Test the DonjonModel with simplified axially varying geometry.
# The GE14 assembly refered to in BWR Progression Problems.
# The "DOM" power dominated zone is modeled form 0.0 to 222.0595 cm, and the "VAN" void dominated zone is modeled from 222.0595 to 347.1291 cm.

def test_axial_variation_simple_core_model():
    path_to_configs = "../data/BWRProgressionProblems/GE14/inputs"
    GE14_core_description_yaml = "GEOM_equiv_CORE.yaml"
    single_assembly_core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=path_to_configs, core_description_yaml=GE14_core_description_yaml)
    assert single_assembly_core_model.geometry_type == "cartesian"
    assert single_assembly_core_model.reactor_type == "BWR"
    assert single_assembly_core_model.core_2D_layout == [["single_GE14_assembly"]]
    #assert single_assembly_core_model.assembly_axial_layouts == [["single_GE14_assembly"]]
    print(single_assembly_core_model.assembly_axial_layouts)
    assert single_assembly_core_model.assembly_axial_layouts == {'single_GE14_assembly': [{'axial_region': 'DOM', 'axial_bounds': [0.0, 222.0595], 'assembly_geometry_file': 'GEOM_GE14_DOM.yaml'}, {'axial_region': 'VAN', 'axial_bounds': [222.0595, 347.1291], 'assembly_geometry_file': 'GEOM_GE14_VAN.yaml'}]}
    assert len(single_assembly_core_model.assemblies) == 1
    assert "single_GE14_assembly" in single_assembly_core_model.assemblies
    assert single_assembly_core_model.assemblies["single_GE14_assembly"].slices_2D == ["2D_slice_geometry_for_DOM_region_of_single_GE14_assembly", "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly"]
    assert single_assembly_core_model.assemblies["single_GE14_assembly"].z_bounds == [0.0, 222.0595, 347.1291]
    assert single_assembly_core_model.assemblies["single_GE14_assembly"].slice_to_geometry_dict == {"2D_slice_geometry_for_DOM_region_of_single_GE14_assembly": "GEOM_GE14_DOM.yaml", "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly": "GEOM_GE14_VAN.yaml"}

    single_assembly_core_model.createAssemblyModels()
    assert len(single_assembly_core_model.assembly_models) == 2
    # assert number of pins in the DOM region model is consistent with the geometry description pointed at in the core description YAML file
    assert single_assembly_core_model.assembly_models[("single_GE14_assembly", "2D_slice_geometry_for_DOM_region_of_single_GE14_assembly")].count_number_of_pins() == 92
    # assert number of pins in the VAN region model is consistent
    assert single_assembly_core_model.assembly_models[("single_GE14_assembly", "2D_slice_geometry_for_VAN_region_of_single_GE14_assembly")].count_number_of_pins() == 78

if __name__ == "__main__":
    test_axial_variation_simple_core_model()
    print("Test passed successfully.")