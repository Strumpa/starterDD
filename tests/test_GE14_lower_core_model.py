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
    tdt_file_name = "GE14_lower_core_IC" #to be produced by glow


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
    assert GE14_assembly.lattice[0][0].rod_ID == "ROD2"
    assert GE14_assembly.lattice[0][0].fuel_material_name == "UOX28"
    assert GE14_assembly.lattice[1][2].fuel_material_name == "UOX44Gd6"
    assert GE14_assembly.lattice[3][3].rod_ID == "WROD"

    # ------------------------------------------------------------------
    # Test: number_fuel_material_mixtures with "by_pin" strategy
    # ------------------------------------------------------------------
    # Parse compositions from the YAML and attach them to the assembly model
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
    GE14_assembly.set_material_compositions(compositions)

    # Run the by_pin numbering
    GE14_assembly.number_fuel_material_mixtures_by_pin()

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

    # --- Check that mixture names follow the by_pin convention ---
    # Names should follow the pattern: <material>_zone<N>_pin<M>
    for name in mixture_names:
        assert "_pin" in name, f"by_pin mixture name '{name}' should contain '_pin'"
        assert "_zone" in name, f"by_pin mixture name '{name}' should contain '_zone'"

    # --- Check that pins have pin_idx attribute assigned ---
    first_fuel_pin = None
    for row in GE14_assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel):
                assert hasattr(pin, 'pin_idx'), "FuelPinModel should have pin_idx after by_pin numbering."
                assert hasattr(pin, 'fuel_material_mixture_indices'), "FuelPinModel should have fuel_material_mixture_indices."
                assert hasattr(pin, 'fuel_material_mixture_names'), "FuelPinModel should have fuel_material_mixture_names."
                if first_fuel_pin is None:
                    first_fuel_pin = pin
                break
        if first_fuel_pin is not None:
            break

    # --- Check diagonal symmetry detection ---
    if hasattr(GE14_assembly, 'lattice_has_diagonal_symmetry'):
        symmetry_type = GE14_assembly.lattice_has_diagonal_symmetry
        print(f"  -> Detected diagonal symmetry: {symmetry_type}")
        
        if symmetry_type is not None:
            # If symmetric, find mirror partners and verify they share pin_idx
            n_rows = len(GE14_assembly.lattice)
            for i in range(n_rows):
                for j in range(len(GE14_assembly.lattice[i])):
                    pin = GE14_assembly.lattice[i][j]
                    if isinstance(pin, FuelPinModel):
                        # Calculate mirror position
                        if symmetry_type == "anti-diagonal":
                            mirror = (n_rows - 1 - j, n_rows - 1 - i)
                        else:  # "main-diagonal"
                            mirror = (j, i)
                        
                        if mirror != (i, j):
                            mirror_pin = GE14_assembly.lattice[mirror[0]][mirror[1]]
                            if isinstance(mirror_pin, FuelPinModel):
                                assert pin.pin_idx == mirror_pin.pin_idx, \
                                    f"Mirror pins at ({i},{j}) and {mirror} should share the same pin_idx"
                                assert pin.fuel_material_mixture_indices == mirror_pin.fuel_material_mixture_indices, \
                                    f"Mirror pins should share the same mixture indices"

    # --- Check that different pins with the same material have DIFFERENT mixture indices ---
    # (This is the key difference from by_material strategy)
    uox28_pins = []
    for row in GE14_assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel) and pin.fuel_material_name == "UOX28":
                uox28_pins.append(pin)
    
    if len(uox28_pins) > 1:
        # Check if any two pins have different pin_idx (i.e., are not symmetric partners)
        non_symmetric_pairs_found = False
        for i in range(len(uox28_pins)):
            for j in range(i + 1, len(uox28_pins)):
                if uox28_pins[i].pin_idx != uox28_pins[j].pin_idx:
                    # These are different pins - they should have different mixture indices
                    assert uox28_pins[i].fuel_material_mixture_indices != uox28_pins[j].fuel_material_mixture_indices, \
                        f"Different pins (pin_idx {uox28_pins[i].pin_idx} vs {uox28_pins[j].pin_idx}) with same material " \
                        f"should have different mixture indices (by_pin rule)"
                    non_symmetric_pairs_found = True
                    break
            if non_symmetric_pairs_found:
                break

    # --- Verify that different materials have different mixture indices ---
    uox28_pin = None
    gd_pin = None
    for row in GE14_assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel):
                if pin.fuel_material_name == "UOX28" and uox28_pin is None:
                    uox28_pin = pin
                elif "Gd" in pin.fuel_material_name and gd_pin is None:
                    gd_pin = pin
        if uox28_pin and gd_pin:
            break

    if uox28_pin and gd_pin:
        assert uox28_pin.fuel_material_mixture_indices != gd_pin.fuel_material_mixture_indices, \
            "Pins with different materials must have different mixture indices."

    # Temperature should reflect the uniform fuel temperature set earlier
    for mix in mixture_objects:
        assert mix.temperature == 900.0, f"Expected fuel temperature 900.0 but got {mix.temperature}."

    print("  -> number_fuel_material_mixtures_by_pin test PASSED")
    print(f"     Created {len(mixture_names)} fuel material mixtures:")
    for name, idx in zip(mixture_names[:10], mixture_indices[:10]):  # Show first 10
        print(f"       index {idx}: {name}")
    if len(mixture_names) > 10:
        print(f"       ... and {len(mixture_names) - 10} more")

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

    print(f"  -> TDT parser test PASSED, found {len(tdt_indices)} entries")

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
    if first_fuel_pin:
        first_zone_name = first_fuel_pin.fuel_material_mixture_names[0]
        expected_idx = tdt_indices[first_zone_name]
        assert first_fuel_pin.fuel_material_mixture_indices[0] == expected_idx, \
            f"Pin zone 1 index should be {expected_idx}, got {first_fuel_pin.fuel_material_mixture_indices[0]}"

    if gd_pin:
        gd_zone_name = gd_pin.fuel_material_mixture_names[0]
        expected_gd_idx = tdt_indices[gd_zone_name]
        assert gd_pin.fuel_material_mixture_indices[0] == expected_gd_idx, \
            f"Gd pin zone 1 index should be {expected_gd_idx}, got {gd_pin.fuel_material_mixture_indices[0]}"

    # Non-fuel entries should be stored
    assert hasattr(GE14_assembly, 'non_fuel_material_mixture_indices')
    assert "COOLANT" in GE14_assembly.non_fuel_material_mixture_indices
    assert "CLAD" in GE14_assembly.non_fuel_material_mixture_indices

    print(f"  -> enforce_material_mixture_indices_from_tdt test PASSED")
    print(f"     Old indices: {old_indices[:10]}..." if len(old_indices) > 10 else f"     Old indices: {old_indices}")
    print(f"     New (TDT) indices: {new_indices[:10]}..." if len(new_indices) > 10 else f"     New (TDT) indices: {new_indices}")
    print(f"     Non-fuel entries: {GE14_assembly.non_fuel_material_mixture_indices}")


if __name__ == "__main__":
    test_full_assembly_model_creation()
    print("Test passed successfully!")