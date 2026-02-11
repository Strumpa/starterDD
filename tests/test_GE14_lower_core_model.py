## test to ensure that the assembly model is working as expected


from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID


# write a test that checks that the assembly model is correctly created and that the lattice description is correctly analyzed.
def test_full_assembly_model_creation():
    # Goal : Create an assembly model and instantiate the DDModel builder
    # Model instantiation steps :
    flux_tracking_option = "TISO"  # Tracking type for the assembly model
    # import macros : true
    include_macros = True
    path_to_yaml_compositions = "../data/BWRProgressionProblems/GE14/inputs/material_compositions.yaml"
    path_to_yaml_geometry = "../data/BWRProgressionProblems/GE14/inputs/GE14_lower_core_geometry.yaml"
    path_to_tdt = "../../glow_data/tdt_data"
    tdt_file_name = "GE14_lower_core" #to be produced by glow


    # Create the assembly model with the given tdt file and lattice description.
    # use CartesianAssemblyModel class method to create the zone_assignment properties.
    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)

    GE14_assembly = CartesianAssemblyModel(name="GE14_assembly",
                                        tdt_file=path_to_tdt + "/" + tdt_file_name + ".tdt",
                                        geometry_description_yaml=path_to_yaml_geometry)
    GE14_assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    GE14_assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0, coolant_temperature=600.0,moderator_temperature=600.0, structural_temperature=600.0)

    GE14_assembly.analyze_lattice_description(build_pins=True)

    assert GE14_assembly.lattice_description is not None
    assert GE14_assembly.lattice is not None
    assert len(GE14_assembly.lattice) == 10
    assert len(GE14_assembly.lattice[0]) == 10
    assert GE14_assembly.lattice[0][0].rod_ID == "ROD1"
    assert GE14_assembly.lattice[0][0].fuel_material_name == "UOX16"
    assert GE14_assembly.lattice[1][2].fuel_material_name == "UOX40Gd8"
    assert GE14_assembly.lattice[3][3].rod_ID == "WROD"

    # ------------------------------------------------------------------
    # Test: number_fuel_material_mixtures with "by_material" strategy
    # ------------------------------------------------------------------
    # Parse compositions from the YAML and attach them to the assembly model
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
    GE14_assembly.set_material_compositions(compositions)

    # Run the by_material numbering
    GE14_assembly.number_fuel_material_mixtures_by_material()

    # --- Check that assembly-level attributes are populated ---
    mixture_names = GE14_assembly.get_fuel_material_mixture_names()
    mixture_indices = GE14_assembly.get_fuel_material_mixture_indices()
    mixture_objects = GE14_assembly.get_fuel_material_mixtures()

    assert len(mixture_names) > 0, "Expected at least one fuel material mixture name."
    assert len(mixture_names) == len(mixture_indices), "Names and indices lists must have the same length."
    assert len(mixture_names) == len(mixture_objects), "Names and objects lists must have the same length."

    # Indices should be sequential starting from 1
    assert mixture_indices == list(range(1, len(mixture_indices) + 1)), "Indices must be sequential starting from 1."

    # Each MaterialMixture must carry a non-None Composition
    for mix in mixture_objects:
        assert mix.composition is not None, f"MaterialMixture '{mix.unique_material_mixture_name}' has no Composition."
        assert mix.material_mixture_index in mixture_indices

    # --- Check that the same-material pins share the same mixture indices ---
    first_uox16_pin = GE14_assembly.lattice[0][0]
    assert hasattr(first_uox16_pin, 'fuel_material_mixture_indices'), "FuelPinModel should have fuel_material_mixture_indices after numbering."
    assert hasattr(first_uox16_pin, 'fuel_material_mixture_names'), "FuelPinModel should have fuel_material_mixture_names after numbering."

    # All UOX16 pins must share identical zone indices / names
    for row in GE14_assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel) and pin.fuel_material_name == "UOX16":
                assert pin.fuel_material_mixture_indices == first_uox16_pin.fuel_material_mixture_indices, \
                    "All pins with the same material must share the same mixture indices (by_material rule)."
                assert pin.fuel_material_mixture_names == first_uox16_pin.fuel_material_mixture_names, \
                    "All pins with the same material must share the same mixture names (by_material rule)."

    # UOX40Gd8 pins must NOT share indices with UOX16 pins
    gd_pin = GE14_assembly.lattice[1][2]
    assert gd_pin.fuel_material_mixture_indices != first_uox16_pin.fuel_material_mixture_indices, \
        "Pins with different materials must have different mixture indices."

    # Temperature should reflect the uniform fuel temperature set earlier
    for mix in mixture_objects:
        assert mix.temperature == 900.0, f"Expected fuel temperature 900.0 but got {mix.temperature}."

    print("  -> number_fuel_material_mixtures (by_material) test PASSED")
    print(f"     Created {len(mixture_names)} fuel material mixtures:")
    for name, idx in zip(mixture_names, mixture_indices):
        print(f"       index {idx}: {name}")

    # ------------------------------------------------------------------
    # Test: read_material_mixture_indices_from_tdt_file parser
    # ------------------------------------------------------------------
    tdt_indices = read_material_mixture_indices_from_tdt_file(
        path_to_tdt,
        tdt_file_name=tdt_file_name,
        tracking_option=flux_tracking_option,
        include_macros=include_macros,
        material_names=None  # get ALL entries
    )

    # The TDT file should contain at least the fuel mixture names
    for fuel_name in mixture_names:
        assert fuel_name in tdt_indices, f"Fuel mixture name '{fuel_name}' not found in TDT file."
    # Indices should all be positive integers
    for name, idx in tdt_indices.items():
        assert isinstance(idx, int) and idx > 0, f"Expected positive int index for '{name}', got {idx}"

    print(f"  -> TDT parser test PASSED, found {len(tdt_indices)} entries: {tdt_indices}")

    # ------------------------------------------------------------------
    # Test: enforce_material_mixture_indices_from_tdt
    # ------------------------------------------------------------------
    # Before enforcement, indices were sequential 1, 2, ...
    old_indices = list(mixture_indices)

    GE14_assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)

    new_indices = GE14_assembly.get_fuel_material_mixture_indices()
    new_mixtures = GE14_assembly.get_fuel_material_mixtures()

    # After enforcement, each fuel mixture index must match the TDT value
    for name, mix in zip(mixture_names, new_mixtures):
        expected_idx = tdt_indices[name]
        assert mix.material_mixture_index == expected_idx, \
            f"MaterialMixture '{name}' index should be {expected_idx}, got {mix.material_mixture_index}"

    # Pin-level indices must also be updated
    first_uox16_pin = GE14_assembly.lattice[0][0]
    expected_uox16_idx = tdt_indices["UOX16_zone_1"]
    assert first_uox16_pin.fuel_material_mixture_indices[0] == expected_uox16_idx, \
        f"Pin UOX16 zone 1 index should be {expected_uox16_idx}, got {first_uox16_pin.fuel_material_mixture_indices[0]}"

    gd_pin = GE14_assembly.lattice[1][2]
    expected_gd_idx = tdt_indices["UOX40Gd8_zone_1"]
    assert gd_pin.fuel_material_mixture_indices[0] == expected_gd_idx, \
        f"Pin UOX40Gd8 zone 1 index should be {expected_gd_idx}, got {gd_pin.fuel_material_mixture_indices[0]}"

    # Non-fuel entries should be stored
    assert hasattr(GE14_assembly, 'non_fuel_material_mixture_indices')
    assert "COOLANT" in GE14_assembly.non_fuel_material_mixture_indices
    assert "CLAD" in GE14_assembly.non_fuel_material_mixture_indices

    print(f"  -> enforce_material_mixture_indices_from_tdt test PASSED")
    print(f"     Old indices: {old_indices}")
    print(f"     New (TDT) indices: {new_indices}")
    print(f"     Non-fuel entries: {GE14_assembly.non_fuel_material_mixture_indices}")


if __name__ == "__main__":
    test_full_assembly_model_creation()
    print("Test passed successfully!")