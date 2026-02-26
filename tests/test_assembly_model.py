## test to ensure that the assembly model is working as expected


from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.InterfaceToDD.dragon_module_calls import LIB
from starterDD.InterfaceToDD.dragon_module_calls import EDI, COMPO, EDI_COMPO

from conftest import (
    GE14_COMPOSITIONS_YAML,
    GE14_SIMPLE_GEOMETRY_YAML,
    GE14_TDT_DIR,
    OUTPUTS_DIR,
)


# write a test that checks that the assembly model is correctly created and that the lattice description is correctly analyzed.
def test_assembly_model_creation():
    # Goal : Create an assembly model and instantiate the DDModel builder
    # Model instantiation steps :
    flux_tracking_option = "TISO"  # Tracking type for the assembly model
    # import macros : true
    include_macros = True
    path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
    path_to_yaml_geometry = GE14_SIMPLE_GEOMETRY_YAML
    path_to_tdt = GE14_TDT_DIR
    tdt_file_name = "GE14_simplified"


    # Create the assembly model with the given tdt file and lattice description.
    # use CartesianAssemblyModel class method to create the zone_assignment properties.
    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)

    GE14_simple_assembly = CartesianAssemblyModel(name="GE14_simple_assembly",
                                        tdt_file=path_to_tdt + "/" + tdt_file_name + ".tdt",
                                        geometry_description_yaml=path_to_yaml_geometry)
    GE14_simple_assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    GE14_simple_assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0, coolant_temperature=600.0,moderator_temperature=600.0, structural_temperature=600.0)

    GE14_simple_assembly.analyze_lattice_description(build_pins=True)

    assert GE14_simple_assembly.lattice_description is not None
    assert GE14_simple_assembly.lattice is not None
    assert len(GE14_simple_assembly.lattice) == 10
    assert len(GE14_simple_assembly.lattice[0]) == 10
    assert GE14_simple_assembly.lattice[0][0].rod_ID == "ROD1"
    assert GE14_simple_assembly.lattice[0][0].fuel_material_name == "UOX16"
    assert GE14_simple_assembly.lattice[1][2].fuel_material_name == "UOX40Gd8"
    assert GE14_simple_assembly.lattice[3][3].rod_ID == "WROD"

    # ------------------------------------------------------------------
    # Test: number_fuel_material_mixtures with "by_material" strategy
    # ------------------------------------------------------------------
    # Parse compositions from the YAML and attach them to the assembly model
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
    GE14_simple_assembly.set_material_compositions(compositions)

    # Run the by_material numbering
    GE14_simple_assembly.number_fuel_material_mixtures_by_material()

    # --- Check that assembly-level attributes are populated ---
    mixture_names = GE14_simple_assembly.get_fuel_material_mixture_names()
    mixture_indices = GE14_simple_assembly.get_fuel_material_mixture_indices()
    mixture_objects = GE14_simple_assembly.get_fuel_material_mixtures()

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
    first_uox16_pin = GE14_simple_assembly.lattice[0][0]
    assert hasattr(first_uox16_pin, 'fuel_material_mixture_indices'), "FuelPinModel should have fuel_material_mixture_indices after numbering."
    assert hasattr(first_uox16_pin, 'fuel_material_mixture_names'), "FuelPinModel should have fuel_material_mixture_names after numbering."

    # All UOX16 pins must share identical zone indices / names
    for row in GE14_simple_assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel) and pin.fuel_material_name == "UOX16":
                assert pin.fuel_material_mixture_indices == first_uox16_pin.fuel_material_mixture_indices, \
                    "All pins with the same material must share the same mixture indices (by_material rule)."
                assert pin.fuel_material_mixture_names == first_uox16_pin.fuel_material_mixture_names, \
                    "All pins with the same material must share the same mixture names (by_material rule)."

    # UOX40Gd8 pins must NOT share indices with UOX16 pins
    gd_pin = GE14_simple_assembly.lattice[1][2]
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

    GE14_simple_assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)

    new_indices = GE14_simple_assembly.get_fuel_material_mixture_indices()
    new_mixtures = GE14_simple_assembly.get_fuel_material_mixtures()

    # After enforcement, each fuel mixture index must match the TDT value
    for name, mix in zip(mixture_names, new_mixtures):
        expected_idx = tdt_indices[name]
        assert mix.material_mixture_index == expected_idx, \
            f"MaterialMixture '{name}' index should be {expected_idx}, got {mix.material_mixture_index}"

    # Pin-level indices must also be updated
    first_uox16_pin = GE14_simple_assembly.lattice[0][0]
    expected_uox16_idx = tdt_indices["UOX16_zone_1"]
    assert first_uox16_pin.fuel_material_mixture_indices[0] == expected_uox16_idx, \
        f"Pin UOX16 zone 1 index should be {expected_uox16_idx}, got {first_uox16_pin.fuel_material_mixture_indices[0]}"

    gd_pin = GE14_simple_assembly.lattice[1][2]
    expected_gd_idx = tdt_indices["UOX40Gd8_zone_1"]
    assert gd_pin.fuel_material_mixture_indices[0] == expected_gd_idx, \
        f"Pin UOX40Gd8 zone 1 index should be {expected_gd_idx}, got {gd_pin.fuel_material_mixture_indices[0]}"

    # Non-fuel entries should be stored
    assert hasattr(GE14_simple_assembly, 'non_fuel_material_mixture_indices')
    assert "COOLANT" in GE14_simple_assembly.non_fuel_material_mixture_indices
    assert "CLAD" in GE14_simple_assembly.non_fuel_material_mixture_indices

    print(f"  -> enforce_material_mixture_indices_from_tdt test PASSED")
    print(f"     Old indices: {old_indices}")
    print(f"     New (TDT) indices: {new_indices}")
    print(f"     Non-fuel entries: {GE14_simple_assembly.non_fuel_material_mixture_indices}")


def test_number_fuel_material_mixtures_by_pin():
    """
    Test the by_pin numbering strategy on the GE14 simplified assembly.
    The GE14 lattice IS anti-diagonally symmetric (axis from top-left to
    bottom-right in physical space), so symmetric partner fuel pins should
    share the same pin_idx, reducing 92 fuel pins to 51 unique positions.
    """
    path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
    path_to_yaml_geometry = GE14_SIMPLE_GEOMETRY_YAML

    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)

    assembly = CartesianAssemblyModel(name="GE14_by_pin",
                                      tdt_file="dummy.tdt",
                                      geometry_description_yaml=path_to_yaml_geometry)
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                      coolant_temperature=600.0, moderator_temperature=600.0,
                                      structural_temperature=600.0)
    assembly.analyze_lattice_description(build_pins=True)

    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
    assembly.set_material_compositions(compositions)

    # --- Verify anti-diagonal symmetry detection ---
    assert assembly.check_diagonal_symmetry() == "anti-diagonal", \
        "The GE14 simplified lattice should be anti-diagonally symmetric."
    assert assembly.check_anti_diagonal_symmetry() is True
    assert assembly.check_main_diagonal_symmetry() is False

    # --- Run by_pin numbering ---
    assembly.number_fuel_material_mixtures_by_pin()

    mixture_names = assembly.get_fuel_material_mixture_names()
    mixture_indices = assembly.get_fuel_material_mixture_indices()
    mixture_objects = assembly.get_fuel_material_mixtures()

    # Count fuel pins in the lattice
    n_fuel_pins = 0
    for row in assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel):
                n_fuel_pins += 1
    # GE14: 100 positions - 8 water rods = 92 fuel pins
    assert n_fuel_pins == 92, f"Expected 92 fuel pins, got {n_fuel_pins}"

    # With anti-diagonal symmetry: 10 on-anti-diagonal + 41 off-diagonal pairs = 51
    n_unique = 51
    assert len(mixture_names) == n_unique, f"Expected {n_unique} mixture names, got {len(mixture_names)}"
    assert len(mixture_indices) == n_unique
    assert len(mixture_objects) == n_unique

    # Indices must be sequential starting from 1
    assert mixture_indices == list(range(1, n_unique + 1)), f"Indices must be sequential 1..{n_unique}"

    # Every name must follow the convention <material>_zone<N>_pin<M>
    for name in mixture_names:
        assert "_zone" in name and "_pin" in name, \
            f"Name '{name}' does not follow the <material>_zone<N>_pin<M> convention."

    # All names must be unique
    assert len(set(mixture_names)) == len(mixture_names), "Mixture names must be unique."

    # Each pin should have a pin_idx; symmetric partners should share the same one
    pin_indices_seen = set()
    for row in assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel):
                assert hasattr(pin, 'pin_idx'), "Each fuel pin must have a pin_idx attribute."
                pin_indices_seen.add(pin.pin_idx)
                # Check that the pin's mixture names contain its pin_idx
                for mname in pin.fuel_material_mixture_names:
                    assert f"_pin{pin.pin_idx}" in mname, \
                        f"Pin mixture name '{mname}' should contain '_pin{pin.pin_idx}'."

    assert len(pin_indices_seen) == n_unique, \
        f"Expected {n_unique} unique pin indices (anti-diagonal symmetry), got {len(pin_indices_seen)}"

    # Verify specific symmetric pair: (i,j) and (n-1-j, n-1-i) share pin_idx
    # e.g. (1,2) and (7,8) should share the same pin_idx
    pin_12 = assembly.lattice[1][2]  # ROD5G
    pin_78 = assembly.lattice[7][8]  # ROD5G
    assert pin_12.pin_idx == pin_78.pin_idx, \
        f"Anti-diagonal partners (1,2) and (7,8) should share pin_idx, got {pin_12.pin_idx} vs {pin_78.pin_idx}"
    assert pin_12.fuel_material_mixture_indices == pin_78.fuel_material_mixture_indices

    # Verify an on-anti-diagonal fuel pin has a unique idx
    # Anti-diagonal: i+j=9, e.g. (0,9)
    pin_09 = assembly.lattice[0][9]
    pin_90 = assembly.lattice[9][0]
    # (0,9) mirrors to (n-1-9, n-1-0) = (0,9) itself → self-paired
    assert pin_09.pin_idx != pin_90.pin_idx or (0, 9) == (9-9, 9-0), \
        "On-anti-diagonal pins should be self-paired."

    # Compositions must be non-None
    for mix in mixture_objects:
        assert mix.composition is not None
        assert mix.temperature == 900.0

    # The lattice_has_diagonal_symmetry flag must be stored
    assert assembly.lattice_has_diagonal_symmetry == "anti-diagonal"

    print("  -> number_fuel_material_mixtures_by_pin test PASSED")
    print(f"     {len(mixture_names)} unique mixtures for {n_fuel_pins} fuel pins (anti-diagonal symmetry)")
    print(f"     First 5 names: {mixture_names[:5]}")
    print(f"     Last 5 names:  {mixture_names[-5:]}")


def test_by_pin_diagonal_symmetry():
    """
    Test the by_pin numbering with a small anti-diagonally symmetric lattice.
    A 3x3 lattice (row 0 = bottom):
        Row 0: ROD1   ROD5G  ROD1
        Row 1: ROD1   WROD   ROD5G
        Row 2: ROD1   ROD1   ROD1
    Anti-diagonal symmetry: lattice[i][j] == lattice[n-1-j][n-1-i].
    Anti-diagonal partners for n=3:
        (0,0) <-> (2,2),  (0,1) <-> (1,2),  (1,0) <-> (2,1)
        On anti-diagonal (i+j=2): (0,2), (1,1), (2,0) are self-paired.
    """
    import tempfile, os, yaml

    # Create a small anti-diagonally symmetric geometry YAML
    geometry_data = {
        "PIN_GEOMETRY": {
            "fuel_radius": 0.438,
            "gap_radius": 0.447,
            "clad_radius": 0.515,
            "pin_pitch": 1.3,
            "self_shielding_option": "automatic",
            "options_dict": {"num_radial_zones": 1},
        },
        "ASSEMBLY_GEOMETRY": {
            "lattice_type": "Cartesian",
            "reactor_type": "BWR",
            "lattice_description": [
                ["ROD1", "ROD5G", "ROD1"],
                ["ROD1", "WROD", "ROD5G"],
                ["ROD1", "ROD1", "ROD1"],
            ],
            "assembly_pitch": 3.9,
            "gap_wide": 0.0,
            "channel_box_thickness": 0.0,
            "Gd_rod_ids": ["ROD5G"],
            "non_fuel_rod_ids": ["WROD"],
        },
    }

    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(geometry_data, f)
        geometry_yaml_path = f.name

    try:
        path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
        ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                       geometry_yaml_path)

        assembly = CartesianAssemblyModel(name="sym_test",
                                          tdt_file="dummy.tdt",
                                          geometry_description_yaml=geometry_yaml_path)
        assembly.set_rod_ID_to_material_mapping(ROD_to_material)
        assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                          coolant_temperature=600.0, moderator_temperature=600.0,
                                          structural_temperature=600.0)
        assembly.analyze_lattice_description(build_pins=True)

        compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
        assembly.set_material_compositions(compositions)

        # --- Verify anti-diagonal symmetry ---
        assert assembly.check_diagonal_symmetry() == "anti-diagonal", \
            "This 3x3 lattice should be anti-diagonally symmetric."

        # --- Run by_pin ---
        assembly.number_fuel_material_mixtures_by_pin()

        assert assembly.lattice_has_diagonal_symmetry == "anti-diagonal"

        # 3x3 = 9 positions, 1 WROD => 8 fuel pins
        n_fuel_pins = 0
        for row in assembly.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    n_fuel_pins += 1
        assert n_fuel_pins == 8, f"Expected 8 fuel pins, got {n_fuel_pins}"

        # Anti-diagonal pairs (n=3, mirror = (2-j, 2-i)):
        #   (0,0) <-> (2,2): pair
        #   (0,1) <-> (1,2): pair (both ROD5G)
        #   (1,0) <-> (2,1): pair
        # On anti-diagonal (self-paired fuel pins):
        #   (0,2): ROD1
        #   (2,0): ROD1
        #   (1,1): WROD (not fuel)
        # Total unique fuel pin positions = 3 pairs + 2 self = 5
        mixture_names = assembly.get_fuel_material_mixture_names()
        assert len(mixture_names) == 5, \
            f"Expected 5 unique mixtures (anti-diagonal symmetry), got {len(mixture_names)}: {mixture_names}"

        # Pins at (0,1) and (1,2) must share the same pin_idx (anti-diag partners)
        pin_01 = assembly.lattice[0][1]
        pin_12 = assembly.lattice[1][2]
        assert pin_01.pin_idx == pin_12.pin_idx, \
            f"Anti-diag partners (0,1) and (1,2) should share pin_idx, got {pin_01.pin_idx} vs {pin_12.pin_idx}"
        assert pin_01.fuel_material_mixture_indices == pin_12.fuel_material_mixture_indices

        # Pins at (0,0) and (2,2) must share the same pin_idx
        pin_00 = assembly.lattice[0][0]
        pin_22 = assembly.lattice[2][2]
        assert pin_00.pin_idx == pin_22.pin_idx, \
            f"Anti-diag partners (0,0) and (2,2) should share pin_idx, got {pin_00.pin_idx} vs {pin_22.pin_idx}"
        assert pin_00.fuel_material_mixture_indices == pin_22.fuel_material_mixture_indices

        # Pins at (1,0) and (2,1) must share the same pin_idx
        pin_10 = assembly.lattice[1][0]
        pin_21 = assembly.lattice[2][1]
        assert pin_10.pin_idx == pin_21.pin_idx, \
            f"Anti-diag partners (1,0) and (2,1) should share pin_idx, got {pin_10.pin_idx} vs {pin_21.pin_idx}"
        assert pin_10.fuel_material_mixture_indices == pin_21.fuel_material_mixture_indices

        # On-anti-diagonal fuel pins (0,2) and (2,0) are each self-paired
        # and must have different pin_idx from each other
        pin_02 = assembly.lattice[0][2]
        pin_20 = assembly.lattice[2][0]
        assert pin_02.pin_idx != pin_20.pin_idx, \
            "Self-paired anti-diagonal pins at different positions should have different pin_idx."

        # All names must follow convention
        for name in mixture_names:
            assert "_zone" in name and "_pin" in name

        # Indices sequential from 1
        mixture_indices = assembly.get_fuel_material_mixture_indices()
        assert mixture_indices == list(range(1, 6))

        print("  -> by_pin anti-diagonal symmetry test PASSED")
        print(f"     {len(mixture_names)} unique mixtures for {n_fuel_pins} fuel pins (3 anti-diag pairs)")
        print(f"     Names: {mixture_names}")
        print(f"     Pin (0,1) idx={pin_01.pin_idx}, Pin (1,2) idx={pin_12.pin_idx}")
    finally:
        os.unlink(geometry_yaml_path)


def test_by_pin_main_diagonal_symmetry():
    """
    Test the by_pin numbering with a small main-diagonally symmetric lattice.
    A 3x3 lattice (transpose symmetric: lattice[i][j] == lattice[j][i]):
        Row 0: ROD1   ROD5G  ROD1
        Row 1: ROD5G  WROD   ROD1
        Row 2: ROD1   ROD1   ROD1
    Main-diagonal partners: (0,1)<->(1,0), (0,2)<->(2,0), (1,2)<->(2,1).
    On main diagonal: (0,0), (2,2) are self-paired, (1,1) is WROD.
    """
    import tempfile, os, yaml

    geometry_data = {
        "PIN_GEOMETRY": {
            "fuel_radius": 0.438,
            "gap_radius": 0.447,
            "clad_radius": 0.515,
            "pin_pitch": 1.3,
            "self_shielding_option": "automatic",
            "options_dict": {"num_radial_zones": 1},
        },
        "ASSEMBLY_GEOMETRY": {
            "lattice_type": "Cartesian",
            "reactor_type": "BWR",
            "lattice_description": [
                ["ROD1", "ROD5G", "ROD1"],
                ["ROD5G", "WROD", "ROD1"],
                ["ROD1", "ROD1", "ROD1"],
            ],
            "assembly_pitch": 3.9,
            "gap_wide": 0.0,
            "channel_box_thickness": 0.0,
            "Gd_rod_ids": ["ROD5G"],
            "non_fuel_rod_ids": ["WROD"],
        },
    }

    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(geometry_data, f)
        geometry_yaml_path = f.name

    try:
        path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
        ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                       geometry_yaml_path)

        assembly = CartesianAssemblyModel(name="main_diag_test",
                                          tdt_file="dummy.tdt",
                                          geometry_description_yaml=geometry_yaml_path)
        assembly.set_rod_ID_to_material_mapping(ROD_to_material)
        assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                          coolant_temperature=600.0, moderator_temperature=600.0,
                                          structural_temperature=600.0)
        assembly.analyze_lattice_description(build_pins=True)

        compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
        assembly.set_material_compositions(compositions)

        # This lattice is main-diag symmetric but NOT anti-diag symmetric
        assert assembly.check_main_diagonal_symmetry() is True
        assert assembly.check_anti_diagonal_symmetry() is False
        assert assembly.check_diagonal_symmetry() == "main-diagonal"

        assembly.number_fuel_material_mixtures_by_pin()

        assert assembly.lattice_has_diagonal_symmetry == "main-diagonal"

        # 8 fuel pins, main-diag pairs: (0,1)<->(1,0), (0,2)<->(2,0), (1,2)<->(2,1)
        # Self-paired: (0,0), (2,2)  => 5 unique
        mixture_names = assembly.get_fuel_material_mixture_names()
        assert len(mixture_names) == 5, \
            f"Expected 5 unique mixtures (main-diagonal symmetry), got {len(mixture_names)}: {mixture_names}"

        # Main-diag partners (0,1) and (1,0) should share pin_idx
        pin_01 = assembly.lattice[0][1]
        pin_10 = assembly.lattice[1][0]
        assert pin_01.pin_idx == pin_10.pin_idx, \
            f"Main-diag partners (0,1) and (1,0) should share pin_idx, got {pin_01.pin_idx} vs {pin_10.pin_idx}"
        assert pin_01.fuel_material_mixture_indices == pin_10.fuel_material_mixture_indices

        # (0,2) and (2,0)
        pin_02 = assembly.lattice[0][2]
        pin_20 = assembly.lattice[2][0]
        assert pin_02.pin_idx == pin_20.pin_idx

        # (1,2) and (2,1)
        pin_12 = assembly.lattice[1][2]
        pin_21 = assembly.lattice[2][1]
        assert pin_12.pin_idx == pin_21.pin_idx

        # On-diagonal (0,0) and (2,2) are self-paired, different from each other
        pin_00 = assembly.lattice[0][0]
        pin_22 = assembly.lattice[2][2]
        assert pin_00.pin_idx != pin_22.pin_idx

        mixture_indices = assembly.get_fuel_material_mixture_indices()
        assert mixture_indices == list(range(1, 6))

        print("  -> by_pin main-diagonal symmetry test PASSED")
        print(f"     {len(mixture_names)} unique mixtures for 8 fuel pins")
        print(f"     Names: {mixture_names}")
    finally:
        os.unlink(geometry_yaml_path)


def test_generating_and_daughter_mixes():
    """
    Test identify_generating_and_daughter_mixes on the GE14 simplified
    assembly using by_material numbering.

    With by_material numbering and 1 zone per pin, there are 2 unique fuel
    materials (UOX16, UOX40Gd8) → 2 material mixtures → 2 generating mixes
    and 0 daughter mixes.

    Then test with by_pin numbering: 51 unique pin positions with 1 zone each.
    2 generating mixes (first UOX16 pin, first UOX40Gd8 pin) and 49 daughters.
    """
    path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
    path_to_yaml_geometry = GE14_SIMPLE_GEOMETRY_YAML

    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)

    # === Test with by_material numbering ===
    assembly_mat = CartesianAssemblyModel(name="GE14_gen_mat",
                                          tdt_file="dummy.tdt",
                                          geometry_description_yaml=path_to_yaml_geometry)
    assembly_mat.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly_mat.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                          coolant_temperature=600.0, moderator_temperature=600.0,
                                          structural_temperature=600.0)
    assembly_mat.analyze_lattice_description(build_pins=True)
    assembly_mat.set_material_compositions(compositions)
    assembly_mat.number_fuel_material_mixtures_by_material()
    assembly_mat.identify_generating_and_daughter_mixes()

    # 2 materials × 1 zone = 2 mixtures, all generating, no daughters
    assert len(assembly_mat.generating_mixes) == 2
    assert len(assembly_mat.daughter_mixes) == 0

    for mix in assembly_mat.generating_mixes:
        assert mix.is_generating is True
        assert mix.generating_mix is None

    print("  -> by_material generating/daughter test PASSED "
          f"({len(assembly_mat.generating_mixes)} generating, "
          f"{len(assembly_mat.daughter_mixes)} daughter)")

    # === Test with by_pin numbering ===
    assembly_pin = CartesianAssemblyModel(name="GE14_gen_pin",
                                          tdt_file="dummy.tdt",
                                          geometry_description_yaml=path_to_yaml_geometry)
    assembly_pin.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly_pin.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                          coolant_temperature=600.0, moderator_temperature=600.0,
                                          structural_temperature=600.0)
    assembly_pin.analyze_lattice_description(build_pins=True)
    assembly_pin.set_material_compositions(compositions)
    assembly_pin.number_fuel_material_mixtures_by_pin()
    assembly_pin.identify_generating_and_daughter_mixes()

    # 51 mixtures total (anti-diagonal symmetry), 2 unique materials
    assert len(assembly_pin.generating_mixes) == 2
    assert len(assembly_pin.daughter_mixes) == 51 - 2

    # Every generating mix must be first appearance of its material
    gen_materials = [m.material_name for m in assembly_pin.generating_mixes]
    assert "UOX16" in gen_materials
    assert "UOX40Gd8" in gen_materials

    # Every daughter must point to a valid generating mix of the same material
    for mix in assembly_pin.daughter_mixes:
        assert mix.is_generating is False
        assert mix.generating_mix is not None
        assert mix.generating_mix.is_generating is True
        assert mix.generating_mix.material_name == mix.material_name

    print("  -> by_pin generating/daughter test PASSED "
          f"({len(assembly_pin.generating_mixes)} generating, "
          f"{len(assembly_pin.daughter_mixes)} daughter)")


def test_lib_write_to_c2m():
    """
    Test the full pipeline: by_pin numbering → TDT enforcement → generating/
    daughter identification → LIB .c2m generation.

    Uses the real TDT file ``GE14_simplified_pin_numbering_TISO_MACRO.dat``
    produced by glow/SALOME so that material mixture indices in the generated
    .c2m are coherent with the SALOME-internal numbering.
    """
    import tempfile, shutil, os

    path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
    path_to_yaml_geometry = GE14_SIMPLE_GEOMETRY_YAML
    path_to_tdt = GE14_TDT_DIR
    tdt_file_name = "GE14_simplified_pin_numbering"
    flux_tracking_option = "TISO"

    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)

    assembly = CartesianAssemblyModel(name="GE14_lib_test",
                                      tdt_file=path_to_tdt + "/" + tdt_file_name + ".tdt",
                                      geometry_description_yaml=path_to_yaml_geometry)
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                      coolant_temperature=600.0, moderator_temperature=600.0,
                                      structural_temperature=600.0)
    assembly.analyze_lattice_description(build_pins=True)
    assembly.set_material_compositions(compositions)

    # ------------------------------------------------------------------
    # Step 1: Number fuel mixtures by pin (sequential 1 .. 51)
    # ------------------------------------------------------------------
    assembly.number_fuel_material_mixtures_by_pin()

    mixture_names_before = list(assembly.get_fuel_material_mixture_names())
    indices_before = list(assembly.get_fuel_material_mixture_indices())
    assert indices_before == list(range(1, 52)), "Pre-TDT indices must be sequential 1..51"

    # ------------------------------------------------------------------
    # Step 2: Read TDT indices from the real glow-produced .dat file
    # ------------------------------------------------------------------
    tdt_indices = read_material_mixture_indices_from_tdt_file(
        path_to_tdt,
        tdt_file_name=tdt_file_name,
        tracking_option=flux_tracking_option,
        include_macros=True,
        material_names=None  # get ALL entries (fuel + non-fuel)
    )

    # The TDT file must contain every fuel mixture name
    for fuel_name in mixture_names_before:
        assert fuel_name in tdt_indices, \
            f"Fuel mixture name '{fuel_name}' not found in TDT file."

    # Non-fuel entries must be present
    for nf_name in ["COOLANT", "CLAD", "MODERATOR", "GAP", "CHANNEL_BOX"]:
        assert nf_name in tdt_indices, f"Non-fuel '{nf_name}' not found in TDT file."

    print(f"  -> TDT parser (by_pin) PASSED: {len(tdt_indices)} entries")

    # ------------------------------------------------------------------
    # Step 3: Enforce TDT indices
    # ------------------------------------------------------------------
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)

    indices_after = assembly.get_fuel_material_mixture_indices()
    # Indices should now match the TDT file (not sequential)
    for name, mix in zip(mixture_names_before, assembly.get_fuel_material_mixtures()):
        expected = tdt_indices[name]
        assert mix.material_mixture_index == expected, \
            f"After enforcement, '{name}' should have index {expected}, got {mix.material_mixture_index}"

    # Pin-level indices must also be updated
    for row in assembly.lattice:
        for pin in row:
            if isinstance(pin, FuelPinModel):
                for zone_name, pin_idx in zip(pin.fuel_material_mixture_names,
                                               pin.fuel_material_mixture_indices):
                    assert pin_idx == tdt_indices[zone_name], \
                        f"Pin zone '{zone_name}': expected {tdt_indices[zone_name]}, got {pin_idx}"

    # Non-fuel indices stored
    assert "COOLANT" in assembly.non_fuel_material_mixture_indices
    assert "CLAD" in assembly.non_fuel_material_mixture_indices
    assert "MODERATOR" in assembly.non_fuel_material_mixture_indices
    assert "GAP" in assembly.non_fuel_material_mixture_indices
    assert "CHANNEL_BOX" in assembly.non_fuel_material_mixture_indices

    print(f"  -> enforce_material_mixture_indices_from_tdt (by_pin) PASSED")
    print(f"     Non-fuel: {assembly.non_fuel_material_mixture_indices}")

    # ------------------------------------------------------------------
    # Step 4: Identify generating / daughter mixes
    # ------------------------------------------------------------------
    assembly.identify_generating_and_daughter_mixes()

    assert len(assembly.generating_mixes) == 2  # UOX16 + UOX40Gd8
    assert len(assembly.daughter_mixes) == 49

    # The generating mix for UOX16 should be the one with TDT index 5
    # (UOX16_zone1_pin1 is the first UOX16 in fuel_material_mixtures)
    gen_uox16 = [m for m in assembly.generating_mixes if m.material_name == "UOX16"][0]
    assert gen_uox16.material_mixture_index == tdt_indices["UOX16_zone1_pin1"]

    gen_gd = [m for m in assembly.generating_mixes if m.material_name == "UOX40Gd8"][0]
    assert gen_gd.material_mixture_index == tdt_indices["UOX40Gd8_zone1_pin13"]

    # ------------------------------------------------------------------
    # Step 5: Build LIB and write .c2m
    # ------------------------------------------------------------------
    lib = LIB(assembly)
    lib.set_isotope_alias("MODERATOR", "H1", "H1_H2O")
    lib.set_isotope_alias("COOLANT", "H1", "H1_H2O")

    # --- Generating mix lines use SALOME indices ---
    gen_lines = lib.build_generating_mix_lines()
    gen_uox16_idx = tdt_indices["UOX16_zone1_pin1"]
    gen_gd_idx = tdt_indices["UOX40Gd8_zone1_pin13"]
    assert f"MIX {gen_uox16_idx} <<TFUEL>>" in gen_lines
    assert f"MIX {gen_gd_idx} <<TFUEL>>" in gen_lines
    assert "U235" in gen_lines
    assert "Gd155" in gen_lines

    # --- Daughter mix lines reference the generating SALOME index ---
    daughter_lines = lib.build_daughter_mix_lines()
    assert daughter_lines.count("COMB") == 49
    # Spot-check: UOX16_zone1_pin2 (TDT idx 15) should be COMB of gen UOX16 (TDT idx 5)
    assert f"MIX {tdt_indices['UOX16_zone1_pin2']} COMB {gen_uox16_idx} 1.0" in daughter_lines
    # UOX40Gd8_zone1_pin21 (TDT idx 17) should be COMB of gen UOX40Gd8 (TDT idx 25)
    assert f"MIX {tdt_indices['UOX40Gd8_zone1_pin21']} COMB {gen_gd_idx} 1.0" in daughter_lines

    # --- Non-fuel lines from TDT enforcement (automatic, not manual) ---
    non_fuel_lines = lib.build_non_fuel_mix_lines()
    assert "NOEV" in non_fuel_lines
    assert f"MIX {tdt_indices['COOLANT']} <<TCOOL>> NOEV" in non_fuel_lines
    assert f"MIX {tdt_indices['MODERATOR']} <<TMODE>> NOEV" in non_fuel_lines
    assert f"MIX {tdt_indices['CLAD']} <<TCLAD>> NOEV" in non_fuel_lines
    assert f"MIX {tdt_indices['CHANNEL_BOX']} <<TBOX>> NOEV" in non_fuel_lines
    assert "H1_H2O" in non_fuel_lines

    # --- NMIX = max index across all (fuel + non-fuel) ---
    max_idx = lib._get_max_mix_index()
    assert max_idx == max(tdt_indices.values()), \
        f"NMIX should be {max(tdt_indices.values())}, got {max_idx}"

    # --- Full lib call ---
    lib_call = lib.build_lib_module_call()
    assert "LIBRARY := LIB: ::" in lib_call
    assert f"NMIX {max_idx}" in lib_call

    # --- Comment block ---
    comment_block = lib.build_mix_index_comment_block()
    assert "(generating)" in comment_block
    assert "(daughter of mix" in comment_block
    # SALOME indices should appear in the comment
    assert f"{gen_uox16_idx:4d} : UOX16_zone1_pin1 (generating)" in comment_block

    # --- Write to file ---
    tmpdir = tempfile.mkdtemp()
    try:
        filepath = lib.write_to_c2m(tmpdir, "MIX_LIB_GE14")
        assert os.path.isfile(filepath)

        with open(filepath, 'r') as f:
            content = f.read()

        # Structure checks
        assert "*PROCEDURE MIX_LIB_GE14.c2m" in content
        assert "PARAMETER LIBRARY ::" in content
        assert "REAL TFUEL := DTFUEL D_TO_R" in content
        assert "MODULE  LIB:" in content
        assert "LIBRARY := LIB: ::" in content
        assert "END: ;" in content
        assert "QUIT ." in content

        # No CLE-2000 FMIX variable declarations
        assert "INTEGER FMIX" not in content

        # Uses SALOME indices (not sequential 1, 2, 3 ...)
        assert f"MIX {gen_uox16_idx} <<TFUEL>>" in content
        assert f"MIX {gen_gd_idx} <<TFUEL>>" in content
        assert f"COMB {gen_uox16_idx} 1.0" in content
        assert f"COMB {gen_gd_idx} 1.0" in content

        # Non-fuel with NOEV and correct SALOME indices
        assert f"MIX {tdt_indices['COOLANT']} <<TCOOL>> NOEV" in content
        assert f"MIX {tdt_indices['CHANNEL_BOX']} <<TBOX>> NOEV" in content

        print(f"  -> LIB write_to_c2m (with TDT enforcement) test PASSED")
        print(f"     Generated: {filepath}")
        print(f"     File size: {os.path.getsize(filepath)} bytes")
        print(f"     NMIX: {max_idx}")
        print(f"     Generating UOX16 idx: {gen_uox16_idx}, UOX40Gd8 idx: {gen_gd_idx}")

        # Also write to the test outputs directory for inspection
        outputs_path = lib.write_to_c2m(str(OUTPUTS_DIR), "MIX_LIB_simple_GE14_by_pin")
        print(f"     Also saved to: {outputs_path}")

    finally:
        shutil.rmtree(tmpdir)

# ---------------------------------------------------------------------------
# Helper: build a by_pin assembly with TDT-enforced indices for EDI tests
# ---------------------------------------------------------------------------
def _build_by_pin_assembly_with_tdt():
    """
    Build and return a GE14 simplified assembly with by_pin numbering and
    TDT-enforced indices.  Shared setup for several EDI/COMPO tests.
    """
    path_to_yaml_compositions = GE14_COMPOSITIONS_YAML
    path_to_yaml_geometry = GE14_SIMPLE_GEOMETRY_YAML
    path_to_tdt = GE14_TDT_DIR
    tdt_file_name = "GE14_simplified_pin_numbering"
    flux_tracking_option = "TISO"

    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                                   path_to_yaml_geometry)
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)

    assembly = CartesianAssemblyModel(name="GE14_edi_test",
                                      tdt_file=path_to_tdt + "/" + tdt_file_name + ".tdt",
                                      geometry_description_yaml=path_to_yaml_geometry)
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                      coolant_temperature=600.0, moderator_temperature=600.0,
                                      structural_temperature=600.0)
    assembly.analyze_lattice_description(build_pins=True)
    assembly.set_material_compositions(compositions)
    assembly.number_fuel_material_mixtures_by_pin()

    tdt_indices = read_material_mixture_indices_from_tdt_file(
        path_to_tdt,
        tdt_file_name=tdt_file_name,
        tracking_option=flux_tracking_option,
        include_macros=True,
        material_names=None
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)
    assembly.identify_generating_and_daughter_mixes()
    return assembly, tdt_indices


# ---------------------------------------------------------------------------
# Tests: EDI MERG MIX vector construction
# ---------------------------------------------------------------------------
def test_edi_merg_mix_by_material():
    """
    Test by_material spatial mode: all mixes sharing the same base composition
    get the same output region index.
    """
    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    edi = EDI("test_by_material", assembly)
    edi.set_isotopes(["U235", "U238"])
    edi.set_spatial_homogenization("by_material")

    vector = edi.build_merg_mix_vector()
    max_mix = edi._get_max_mix_index()

    assert len(vector) == max_mix, f"Vector length {len(vector)} != max_mix {max_mix}"

    # All fuel mixes should have a value > 0
    # All non-fuel mixes should have value 0
    non_fuel_idxs = set(assembly.non_fuel_material_mixture_indices.values())

    # Collect unique family IDs assigned to fuel mixes
    family_ids = set()
    for mix in assembly.fuel_material_mixtures:
        idx = mix.material_mixture_index
        val = vector[idx - 1]
        assert val > 0, f"Fuel mix {mix.unique_material_mixture_name} (idx {idx}) should be > 0, got {val}"
        family_ids.add(val)

    # Two unique fuel materials (UOX16, UOX40Gd8) → 2 family IDs
    assert len(family_ids) == 2, f"Expected 2 family IDs, got {len(family_ids)}: {family_ids}"

    # Non-fuel should be 0
    for nf_idx in non_fuel_idxs:
        assert vector[nf_idx - 1] == 0, f"Non-fuel idx {nf_idx} should be 0, got {vector[nf_idx - 1]}"

    # Mixes with same material_name must have same value
    material_to_value = {}
    for mix in assembly.fuel_material_mixtures:
        val = vector[mix.material_mixture_index - 1]
        if mix.material_name in material_to_value:
            assert material_to_value[mix.material_name] == val, \
                f"All mixes of {mix.material_name} should have same family id"
        else:
            material_to_value[mix.material_name] = val

    print("  -> test_edi_merg_mix_by_material PASSED")


def test_edi_merg_mix_by_pin():
    """
    Test by_pin spatial mode: each unique pin position gets its own output
    region, zones within the same pin are merged.
    """
    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    edi = EDI("test_by_pin", assembly)
    edi.set_isotopes(["U235", "U238"])
    edi.set_spatial_homogenization("by_pin")

    vector = edi.build_merg_mix_vector()
    max_mix = edi._get_max_mix_index()

    assert len(vector) == max_mix

    # Non-fuel → 0
    non_fuel_idxs = set(assembly.non_fuel_material_mixture_indices.values())
    for nf_idx in non_fuel_idxs:
        assert vector[nf_idx - 1] == 0

    # Check that fuel mixes carry the pin_idx from their name
    for mix in assembly.fuel_material_mixtures:
        expected_pin_idx = int(mix.unique_material_mixture_name.rsplit("_pin", 1)[1])
        actual = vector[mix.material_mixture_index - 1]
        assert actual == expected_pin_idx, \
            f"Mix '{mix.unique_material_mixture_name}' should have pin_idx {expected_pin_idx}, got {actual}"

    # Unique output values = number of unique pin positions = 51
    unique_output_ids = set(v for v in vector if v > 0)
    assert len(unique_output_ids) == 51, f"Expected 51 unique pin ids, got {len(unique_output_ids)}"

    print("  -> test_edi_merg_mix_by_pin PASSED")


def test_edi_merg_mix_fuel_and_all():
    """
    Test FUEL and ALL spatial modes.
    """
    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    # --- FUEL mode ---
    edi_fuel = EDI("test_fuel", assembly)
    edi_fuel.set_isotopes(["U235"])
    edi_fuel.set_spatial_homogenization("FUEL")

    vector_fuel = edi_fuel.build_merg_mix_vector()

    # All fuel positions → 1
    for mix in assembly.fuel_material_mixtures:
        assert vector_fuel[mix.material_mixture_index - 1] == 1

    # Non-fuel → 0
    for nf_idx in assembly.non_fuel_material_mixture_indices.values():
        assert vector_fuel[nf_idx - 1] == 0

    # Only one unique output region
    assert set(v for v in vector_fuel if v > 0) == {1}

    # --- ALL mode ---
    edi_all = EDI("test_all", assembly)
    edi_all.set_isotopes(["U235"])
    edi_all.set_spatial_homogenization("ALL")

    vector_all = edi_all.build_merg_mix_vector()

    # Every position → 1
    assert all(v == 1 for v in vector_all), "ALL mode: every position should be 1"

    print("  -> test_edi_merg_mix_fuel_and_all PASSED")


def test_edi_merg_mix_by_mix():
    """
    Test by_mix spatial mode: each MaterialMixture gets its own sequential
    output region.
    """
    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    edi = EDI("test_by_mix", assembly)
    edi.set_isotopes(["U235"])
    edi.set_spatial_homogenization("by_mix")

    vector = edi.build_merg_mix_vector()

    # Each fuel mix gets a unique sequential rank
    ranks_seen = set()
    for rank, mix in enumerate(assembly.fuel_material_mixtures, start=1):
        val = vector[mix.material_mixture_index - 1]
        assert val == rank, f"Mix rank {rank} should have value {rank}, got {val}"
        ranks_seen.add(val)

    assert len(ranks_seen) == len(assembly.fuel_material_mixtures)

    # Non-fuel → 0
    for nf_idx in assembly.non_fuel_material_mixture_indices.values():
        assert vector[nf_idx - 1] == 0

    print("  -> test_edi_merg_mix_by_mix PASSED")


def test_edi_merg_mix_custom():
    """
    Test custom spatial mode: user supplies an explicit mapping.
    """
    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    # Build a custom map: first 3 fuel mixes → output regions 10, 20, 30
    first_three = assembly.fuel_material_mixtures[:3]
    custom_map = {
        first_three[0].material_mixture_index: 10,
        first_three[1].material_mixture_index: 20,
        first_three[2].material_mixture_index: 30,
    }

    edi = EDI("test_custom", assembly)
    edi.set_isotopes(["U235"])
    edi.set_spatial_homogenization("custom", custom_map=custom_map)

    vector = edi.build_merg_mix_vector()

    for tdt_idx, output_idx in custom_map.items():
        assert vector[tdt_idx - 1] == output_idx, \
            f"Custom map: idx {tdt_idx} should be {output_idx}, got {vector[tdt_idx - 1]}"

    # Positions not in custom_map should be 0
    mapped_positions = set(k - 1 for k in custom_map.keys())
    for i, v in enumerate(vector):
        if i not in mapped_positions:
            assert v == 0, f"Unmapped position {i+1} should be 0, got {v}"

    print("  -> test_edi_merg_mix_custom PASSED")


def test_edi_energy_condensation_options():
    """
    Test EDI energy condensation keyword generation.
    """
    assembly, _ = _build_by_pin_assembly_with_tdt()

    # None → no COND
    edi = EDI("test", assembly)
    edi.set_energy_condensation(None)
    assert edi.build_cond_option() == ""

    # [] → COND
    edi.set_energy_condensation([])
    assert edi.build_cond_option() == "COND"

    # [0.625] → COND 0.625
    edi.set_energy_condensation([0.625])
    assert edi.build_cond_option() == "COND 0.625"

    # Multiple bounds
    edi.set_energy_condensation([0.625, 100.0])
    assert edi.build_cond_option() == "COND 0.625 100.0"

    # Verify COND appears/does not appear in full EDI call
    edi.set_isotopes(["U235"])
    edi.set_spatial_homogenization("FUEL")

    edi.set_energy_condensation(None)
    call_no_cond = edi.build_edi_call()
    assert "COND" not in call_no_cond

    edi.set_energy_condensation([0.625])
    call_with_cond = edi.build_edi_call()
    assert "COND 0.625" in call_with_cond

    print("  -> test_edi_energy_condensation_options PASSED")


def test_edi_spatial_mode_validation():
    """
    Test that invalid spatial modes raise ValueError.
    """
    assembly, _ = _build_by_pin_assembly_with_tdt()

    edi = EDI("test", assembly)
    import pytest
    with pytest.raises(ValueError, match="Invalid spatial mode"):
        edi.set_spatial_homogenization("invalid_mode")

    with pytest.raises(ValueError, match="custom_map is required"):
        edi.set_spatial_homogenization("custom")

    print("  -> test_edi_spatial_mode_validation PASSED")


def test_compo_init_and_store():
    """
    Test COMPO initialization and store block generation.
    """
    compo = COMPO()
    compo.add_directory("EDIHOM_COND",
                        "Condensed, Homogenized over all fuel cells",
                        ["U235", "U238", "Gd155"])
    compo.add_directory("EDIHOM_295",
                        "Homogenized, 295g",
                        ["U235", "U238"])

    init_block = compo.build_compo_init()

    # Structural checks
    assert "COMPO := COMPO: ::" in init_block
    assert "STEP UP 'EDIHOM_COND'" in init_block
    assert "STEP UP 'EDIHOM_295'" in init_block
    assert "ISOT 3 U235 U238 Gd155" in init_block
    assert "ISOT 2 U235 U238" in init_block
    assert "COMM 'Condensed, Homogenized over all fuel cells' ENDC" in init_block
    assert "COMM 'Homogenized, 295g' ENDC" in init_block
    assert init_block.count("INIT") == 2

    # Store block
    store = compo.build_compo_store("EDIHOM_COND")
    assert "COMPO := COMPO: COMPO EDIRATES LIBRARY2 ::" in store
    assert "STEP UP EDIHOM_COND" in store

    print("  -> test_compo_init_and_store PASSED")


def test_edi_compo_write_to_c2m():
    """
    Test full EDI_COMPO pipeline: add multiple editions, write .c2m,
    verify structural correctness.
    """
    import tempfile, shutil, os

    assembly, tdt_indices = _build_by_pin_assembly_with_tdt()

    edi_compo = EDI_COMPO(assembly)

    edi_compo.add_edition(
        name="EDIHOM_COND",
        comment="Condensed, Homogenized over all fuel cells",
        isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
        spatial_mode="FUEL",
        energy_bounds=[],
    )
    edi_compo.add_edition(
        name="EDIHOM_295",
        comment="Homogenized, 295g",
        isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
        spatial_mode="ALL",
        energy_bounds=None,
    )
    edi_compo.add_edition(
        name="H_EDI_REGI_1g",
        comment="Condensed, per unique pin",
        isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
        spatial_mode="by_pin",
        energy_bounds=[],
    )
    edi_compo.add_edition(
        name="H_EDI_REGI_2g",
        comment="Condensed to 2 groups, per unique pin",
        isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
        spatial_mode="by_pin",
        energy_bounds=[0.625],
    )

    tmpdir = tempfile.mkdtemp()
    try:
        filepath = edi_compo.write_to_c2m(tmpdir, "EDICPO_GE14_test")
        assert os.path.isfile(filepath)

        with open(filepath, 'r') as f:
            content = f.read()

        # Header checks
        assert "PARAMETER COMPO FLUX LIBRARY2 TRACK" in content
        assert "MODULE EDI: COMPO: DELETE: END:" in content
        assert "LINKED_LIST EDIRATES" in content
        assert "SEQ_ASCII _COMPO :: FILE <<name_cpo>>" in content

        # COMPO initialization
        assert "COMPO := COMPO: ::" in content
        assert "STEP UP 'EDIHOM_COND'" in content
        assert "STEP UP 'EDIHOM_295'" in content
        assert "STEP UP 'H_EDI_REGI_1g'" in content
        assert "STEP UP 'H_EDI_REGI_2g'" in content
        assert content.count("\n   INIT\n") == 4  # 4 directories

        # EDI calls
        assert content.count("EDIRATES := EDI: FLUX LIBRARY2 TRACK") == 4
        assert "SAVE ON EDIHOM_COND" in content
        assert "SAVE ON EDIHOM_295" in content
        assert "SAVE ON H_EDI_REGI_1g" in content
        assert "SAVE ON H_EDI_REGI_2g" in content

        # COMPO store calls
        assert content.count("COMPO := COMPO: COMPO EDIRATES LIBRARY2") == 4

        # Cleanup
        assert content.count("EDIRATES := DELETE: EDIRATES") == 4

        # Footer
        assert "IF save_opt 'SAVE' = THEN" in content
        assert "_COMPO := COMPO" in content
        assert "END: ;" in content

        # COND keyword presence
        assert "COND 0.625" in content  # H_EDI_REGI_2g
        # EDIHOM_295 has no COND
        # EDIHOM_COND and H_EDI_REGI_1g have "COND" (1 group)

        # MERG MIX vectors should be present
        assert "MERG MIX" in content

        # Also write to outputs for inspection
        filepath_out = edi_compo.write_to_c2m(str(OUTPUTS_DIR), "EDICPO_GE14_simple_by_pin")

        print(f"  -> test_edi_compo_write_to_c2m PASSED")
        print(f"     Generated: {filepath}")
        print(f"     Also saved to: {filepath_out}")

    finally:
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_assembly_model_creation()
    test_number_fuel_material_mixtures_by_pin()
    test_by_pin_diagonal_symmetry()
    test_by_pin_main_diagonal_symmetry()
    test_generating_and_daughter_mixes()
    test_lib_write_to_c2m()
    test_edi_merg_mix_by_material()
    test_edi_merg_mix_by_pin()
    test_edi_merg_mix_fuel_and_all()
    test_edi_merg_mix_by_mix()
    test_edi_merg_mix_custom()
    test_edi_energy_condensation_options()
    test_edi_spatial_mode_validation()
    test_compo_init_and_store()
    test_edi_compo_write_to_c2m()
    print("All tests passed successfully!")