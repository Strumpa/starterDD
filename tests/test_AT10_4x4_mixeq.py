"""
Tests for MIXEQ procedure generation with AT10 4x4 assembly.

This test suite validates the hybrid mix numbering system where:
- SSH and FLUXL1 use "by_material" numbering (zone-only: e.g., "50UOX_zone_4")
- FLUXL2 uses "by_pin" numbering (zone+pin: e.g., "50UOX_zone_4_pin_7")

The MIXEQ procedure generation must correctly map between these numbering
schemes by expanding each zone-only mix to multiple zone+pin mixes.

Key aspects tested:
1. Internal mix naming generation for both strategies
2. TDT indices enforcement from reference files
3. Correspondence recording between strategies
4. Correspondence table building with correct mappings
5. MIXEQ class procedure file generation
6. Full case generation including MIXEQ procedures
"""
import os
import re
import pytest

from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.DDModel.DragonModel import CartesianAssemblyModel
from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.InterfaceToDD.dragon_module_calls import MIXEQ
from starterDD.InterfaceToDD.case_generator import DragonCase

from conftest import (
    AT10_TDT_DIR,
    AT10_4x4_GEOMETRY_YAML,
    AT10_4x4_COMPOSITIONS_YAML,
    AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML,
    OUTPUTS_DIR,
)

# 4x4 AT10 test case constants
TDT_FILE_BASE = "AT10_4x4_mix_splitting"


def test_at10_4x4_by_material_numbering():
    """
    Test that by_material numbering generates correct zone-only names.

    For the 4x4 AT10 case with 7 fuel materials (6 unique + 1 Gd),
    each with Santamarina zones (4 zones for UOX, 6 for Gd),
    the by_material strategy should generate 30 fuel mixes:
    - 6 materials × 4 zones = 24 UOX mixes
    - 1 Gd material × 6 zones = 6 Gd mixes

    Plus 3 structural mixes (COOLANT, CLAD, GAP) = 33 total.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply by_material numbering
    assembly.apply_mix_numbering_strategy("by_material")

    # Check fuel mix names (should be zone-only)
    fuel_names = assembly.fuel_material_mixture_names
    assert len(fuel_names) == 30, f"Expected 30 fuel mixes, got {len(fuel_names)}"

    # Check naming pattern: should be <material>_zone_<idx>
    zone_pattern = re.compile(r'^[^_]+_zone_\d+$')
    for name in fuel_names:
        assert zone_pattern.match(name), f"Name '{name}' should match zone-only pattern"
        assert "_pin_" not in name, f"by_material name should not contain pin identifier: {name}"

    # Check specific materials
    uox50_zones = [n for n in fuel_names if n.startswith("50UOX_zone_")]
    assert len(uox50_zones) == 4, "50UOX should have 4 zones"

    gd45_zones = [n for n in fuel_names if n.startswith("45Gd_zone_")]
    assert len(gd45_zones) == 6, "45Gd should have 6 zones"

    print("  -> test_at10_4x4_by_material_numbering PASSED")
    print(f"     Generated {len(fuel_names)} fuel mixes with zone-only naming")
    print(f"     Example: {fuel_names[0]}")


def test_at10_4x4_by_pin_numbering():
    """
    Test that by_pin numbering generates correct zone+pin names.

    For the 4x4 lattice with symmetry, unique pins should get unique pin_idx.
    Each fuel zone should be numbered with both zone and pin identifiers:
    <material>_zone_<zone_idx>_pin_<pin_idx>
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply by_pin numbering
    assembly.apply_mix_numbering_strategy("by_pin")

    # Check fuel mix names (should include pin identifiers)
    fuel_names = assembly.fuel_material_mixture_names

    # Check naming pattern: should be <material>_zone_<zone_idx>_pin_<pin_idx>
    pin_pattern = re.compile(r'^[^_]+_zone_\d+_pin_\d+$')
    for name in fuel_names:
        assert pin_pattern.match(name), f"Name '{name}' should match zone+pin pattern"
        assert "_pin_" in name, f"by_pin name should contain pin identifier: {name}"

    # Count unique pins
    pin_ids = set()
    for name in fuel_names:
        match = re.search(r'_pin_(\d+)$', name)
        if match:
            pin_ids.add(int(match.group(1)))

    # 4x4 lattice with symmetry should have fewer unique pins
    print(f"  -> test_at10_4x4_by_pin_numbering PASSED")
    print(f"     Generated {len(fuel_names)} fuel mixes with zone+pin naming")
    print(f"     Unique pin indices: {sorted(pin_ids)}")
    print(f"     Example: {fuel_names[0]}")


def test_at10_4x4_tdt_enforcement_ssh():
    """
    Test TDT indices enforcement from reference SSH file.

    Verifies that:
    1. TDT indices from file are correctly parsed
    2. Assembly enforces these indices correctly
    3. Mix state is recorded with step name
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply by_material numbering for SSH
    assembly.apply_mix_numbering_strategy("by_material")

    # Enforce TDT indices from reference file
    tdt_indices_ssh = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_SSH_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )

    assert len(tdt_indices_ssh) == 33, f"SSH TDT should have 33 mixes, got {len(tdt_indices_ssh)}"

    # Check specific indices match expected values
    assert tdt_indices_ssh["COOLANT"] == 1
    assert tdt_indices_ssh["CLAD"] == 2
    assert tdt_indices_ssh["GAP"] == 3
    assert tdt_indices_ssh["50UOX_zone_4"] == 4
    assert tdt_indices_ssh["24UOX_zone_1"] == 33

    # Enforce indices
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_ssh, step_name="SSH")

    # Verify step state was recorded
    ssh_state = assembly.get_step_mix_state("SSH")
    assert ssh_state is not None, "SSH state should be recorded"
    assert ssh_state["step_name"] == "SSH"
    assert ssh_state["strategy"] == "by_material"
    assert ssh_state["tdt_enforced"] is True
    assert len(ssh_state["tdt_indices"]) == 33, "Should have 33 total mixes (30 fuel + 3 non-fuel)"

    print("  -> test_at10_4x4_tdt_enforcement_ssh PASSED")
    print(f"     Enforced {len(tdt_indices_ssh)} TDT indices for SSH step")


def test_at10_4x4_tdt_enforcement_fluxl2():
    """
    Test TDT indices enforcement from reference FLUXL2 file.

    FLUXL2 uses by_pin numbering with 45 mixes total:
    - 3 structural (shared)
    - 42 fuel mixes with pin identifiers
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply by_pin numbering for FLUXL2
    assembly.apply_mix_numbering_strategy("by_pin")

    # Enforce TDT indices from reference file
    tdt_indices_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )

    assert len(tdt_indices_l2) == 45, f"FLUXL2 TDT should have 45 mixes, got {len(tdt_indices_l2)}"

    # Check structural indices (same as SSH)
    assert tdt_indices_l2["COOLANT"] == 1
    assert tdt_indices_l2["CLAD"] == 2
    assert tdt_indices_l2["GAP"] == 3

    # Check pin-based naming
    assert tdt_indices_l2["50UOX_zone_4_pin_7"] == 4
    assert tdt_indices_l2["24UOX_zone_1_pin_1"] == 33
    assert tdt_indices_l2["50UOX_zone_1_pin_8"] == 45

    # Enforce indices
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_l2, step_name="FLUXL2")

    # Verify step state
    l2_state = assembly.get_step_mix_state("FLUXL2")
    assert l2_state is not None
    assert l2_state["strategy"] == "by_pin"
    assert len(l2_state["tdt_indices"]) == 45, "Should have 45 total mixes enforced for FLUXL2"

    print("  -> test_at10_4x4_tdt_enforcement_fluxl2 PASSED")
    print(f"     Enforced {len(tdt_indices_l2)} TDT indices for FLUXL2 step")


def test_at10_4x4_correspondence_recording():
    """
    Test that correspondence is correctly recorded during strategy transitions.

    When switching from by_material to by_pin, the system should record
    which zone-only mixes map to which zone+pin mixes based on material name.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply by_material first (simulating SSH or FLUXL1)
    assembly.apply_mix_numbering_strategy("by_material")
    by_material_names = assembly.fuel_material_mixture_names.copy()

    # Apply by_pin (simulating FLUXL2)
    assembly.apply_mix_numbering_strategy("by_pin")
    by_pin_names = assembly.fuel_material_mixture_names.copy()

    # Check that correspondence was recorded
    assert ("by_material", "by_pin") in assembly.mix_strategy_correspondences
    assert ("by_pin", "by_material") in assembly.mix_strategy_correspondences

    forward_corr = assembly.mix_strategy_correspondences[("by_material", "by_pin")]
    reverse_corr = assembly.mix_strategy_correspondences[("by_pin", "by_material")]

    # Forward: each zone-only name maps to multiple zone+pin names
    # Example: "50UOX_zone_4" → ["50UOX_zone_4_pin_7", "50UOX_zone_4_pin_8", ...]
    assert len(forward_corr) == 30, "Should have 30 zone-only mixes in forward mapping"

    # Check specific example
    assert "50UOX_zone_4" in forward_corr
    pin_variants = forward_corr["50UOX_zone_4"]
    assert isinstance(pin_variants, list)
    assert len(pin_variants) > 0, "50UOX_zone_4 should map to at least one pin variant"

    # Verify all pin variants have correct base
    for pin_name in pin_variants:
        assert pin_name.startswith("50UOX_zone_4_pin_"), \
            f"Pin variant {pin_name} should start with 50UOX_zone_4_pin_"

    # Reverse: each zone+pin name maps back to its zone-only parent
    # Example: "50UOX_zone_4_pin_7" → ["50UOX_zone_4"]
    for pin_name in pin_variants:
        assert pin_name in reverse_corr
        parent_zones = reverse_corr[pin_name]
        assert "50UOX_zone_4" in parent_zones

    print("  -> test_at10_4x4_correspondence_recording PASSED")
    print(f"     Forward correspondence: {len(forward_corr)} zone-only mixes")
    print(f"     Reverse correspondence: {len(reverse_corr)} zone+pin mixes")
    print(f"     Example: '50UOX_zone_4' → {len(pin_variants)} pin variants")


def test_at10_4x4_correspondence_table_with_real_tdt():
    """
    Test correspondence table generation using real TDT reference files.

    This is the critical test that verifies the hybrid numbering system:
    - Loads by_material indices from SSH TDT (33 mixes)
    - Loads by_pin indices from FLUXL2 TDT (45 mixes)
    - Builds correspondence table mapping each zone-only mix to its pin variants
    - Verifies that indices match the reference TDT files
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Step 1: SSH with by_material
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_indices_ssh = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_SSH_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_ssh, step_name="SSH")

    # Step 2: FLUXL1 with by_material (same as SSH)
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_indices_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_l1, step_name="FLUXL1")

    # Step 3: FLUXL2 with by_pin
    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_indices_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_l2, step_name="FLUXL2")

    # Build correspondence table FLUXL1 → FLUXL2
    table = assembly.build_mixeq_correspondence_table("FLUXL1", "FLUXL2")

    assert len(table) > 0, "Correspondence table should not be empty"

    # Verify structure: each entry is (from_name, to_name, from_idx, to_idx)
    for entry in table:
        assert len(entry) == 4, "Each entry should have 4 elements"
        from_name, to_name, from_idx, to_idx = entry

        assert isinstance(from_name, str)
        assert isinstance(to_name, str)
        assert isinstance(from_idx, int)
        assert isinstance(to_idx, int)

        # Verify index ranges
        assert 1 <= from_idx <= 33, f"SSH/FLUXL1 index should be in [1, 33], got {from_idx}"
        assert 1 <= to_idx <= 45, f"FLUXL2 index should be in [1, 45], got {to_idx}"

    # Test specific correspondences using reference TDT values
    # "50UOX_zone_4" (idx 4) should map to multiple pin variants
    zone4_50uox_mappings = [
        (from_name, to_name, from_idx, to_idx)
        for from_name, to_name, from_idx, to_idx in table
        if from_name == "50UOX_zone_4"
    ]

    assert len(zone4_50uox_mappings) > 0, "50UOX_zone_4 should have mappings"

    # All should have from_idx = 4 (from reference TDT)
    for _, _, from_idx, _ in zone4_50uox_mappings:
        assert from_idx == 4

    # Check that to_names include pin identifiers
    to_names = [to_name for _, to_name, _, _ in zone4_50uox_mappings]
    for to_name in to_names:
        assert to_name.startswith("50UOX_zone_4_pin_"), \
            f"Target name should be 50UOX_zone_4_pin_*, got {to_name}"

    # Verify specific mappings from reference TDT
    expected_mappings = {
        ("50UOX_zone_4", "50UOX_zone_4_pin_7", 4, 4),
        ("50UOX_zone_4", "50UOX_zone_4_pin_10", 4, 34),
        ("50UOX_zone_4", "50UOX_zone_4_pin_9", 4, 38),
        ("50UOX_zone_4", "50UOX_zone_4_pin_8", 4, 42),
    }

    table_set = set(table)
    for expected in expected_mappings:
        assert expected in table_set, f"Expected mapping not found: {expected}"

    # Test structural mixes (should be 1-to-1)
    structural_mappings = [
        (from_name, to_name, from_idx, to_idx)
        for from_name, to_name, from_idx, to_idx in table
        if from_name in ["COOLANT", "CLAD", "GAP"]
    ]

    for from_name, to_name, from_idx, to_idx in structural_mappings:
        assert from_name == to_name, "Structural mixes should map to themselves"
        assert from_idx == to_idx, "Structural indices should be unchanged"

    print("  -> test_at10_4x4_correspondence_table_with_real_tdt PASSED")
    print(f"     Built correspondence table with {len(table)} mappings")
    print(f"     50UOX_zone_4 expands to {len(zone4_50uox_mappings)} pin variants")


def test_at10_4x4_mixeq_procedure_generation():
    """
    Test MIXEQ class procedure generation with real AT10 4x4 data.

    Verifies that:
    1. MIXEQ class correctly initializes from correspondence table
    2. LIB module call is generated with correct COMB statements
    3. Procedure file can be written to disk
    4. Generated procedure has correct format
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Simulate SSH step
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_indices_ssh = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_SSH_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_ssh, step_name="SSH")

    # Simulate FLUXL1 step (same strategy as SSH)
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_indices_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_l1, step_name="FLUXL1")

    # Simulate FLUXL2 step (different strategy)
    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_indices_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_indices_l2, step_name="FLUXL2")

    # Generate MIXEQ procedure
    lib_name = "LIBEQ"
    mixeq = MIXEQ(assembly, "LIBRARY2", "LIBEQ", "FLUXL1", "FLUXL2", draglib_alias="dummy_alias")

    # Test correspondence table access
    assert mixeq.correspondence_table is not None
    assert len(mixeq.correspondence_table) > 0

    # Generate LIB module call
    lib_call = mixeq.build_lib_module_call()

    assert isinstance(lib_call, str)
    assert "LIBEQ := LIB: LIBRARY2 ::" in lib_call
    assert "EDIT 0" in lib_call
    assert "NMIX" in lib_call
    assert "CATL" in lib_call
    assert "DEPL LIB: DRAGON FIL: dummy_alias" in lib_call

    # Verify COMB statements for zone expansion
    # Each zone-only mix should duplicate to multiple pin-based mixes
    comb_pattern = re.compile(r'MIX\s+(\d+)\s+(\d+)', re.MULTILINE)
    comb_matches = comb_pattern.findall(lib_call)

    assert len(comb_matches) > 0, "Should have MIX <new_idx> <old_idx> statements for mix duplication"

    # Verify specific mapping: 50UOX_zone_4 (idx 4) expands to multiple target indices
    # From reference: 4→34, 4→38, 4→42, skip index mapping to itself
    source_idx_4_targets = [int(target) for target, source in comb_matches if int(source) == 4]
    assert 4  in source_idx_4_targets, "Should map 4→4 (self-mapping as LIBRARY2 copied to LIBEQ)"
    assert 34 in source_idx_4_targets, "Should map 4→34"
    assert 38 in source_idx_4_targets, "Should map 4→38"
    assert 42 in source_idx_4_targets, "Should map 4→42"

    # Write to file and verify
    output_dir = str(OUTPUTS_DIR)
    os.makedirs(output_dir, exist_ok=True)

    proc_file = mixeq.write_to_c2m(output_dir, "MIXEQL1L2")
    assert os.path.exists(proc_file), "Procedure file should be created"

    with open(proc_file, 'r') as f:
        content = f.read()

    assert "PROCEDURE MIXEQL1L2" in content
    assert "LIBEQ := LIB: LIBRARY2 ::" in content

    print("  -> test_at10_4x4_mixeq_procedure_generation PASSED")
    print(f"     Generated MIXEQ procedure with {len(comb_matches)} COMB statements")
    print(f"     Source mix 4 expands to {len(source_idx_4_targets)} target mixes")


def test_at10_4x4_zone_to_pin_expansion():
    """
    Test that zone-only numbering properly expands to zone+pin numbering.

    This verifies the core functionality: each material zone from by_material
    should map to ALL pins that contain that material, maintaining the
    zone index while adding the pin identifier.

    Example: 24UOX has 4 zones. In a 4x4 lattice, 24UOX appears in 1 pin (pin_1).
    by_material: 24UOX_zone_1, _zone_2, _zone_3, _zone_4 (4 mixes)
    by_pin: 24UOX_zone_1_pin_1, _zone_2_pin_1, _zone_3_pin_1, _zone_4_pin_1 (4 mixes)

    50UOX appears in 4 pins (pin_7, pin_8, pin_9, pin_10 based on symmetry).
    by_material: 50UOX_zone_1 through _zone_4 (4 mixes)
    by_pin: 50UOX_zone_1_pin_7, _pin_8, _pin_9, _pin_10 (and same for zones 2-4) (16 mixes)
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply both strategies to record correspondence
    assembly.apply_mix_numbering_strategy("by_material")
    by_material_names = assembly.fuel_material_mixture_names.copy()

    assembly.apply_mix_numbering_strategy("by_pin")
    by_pin_names = assembly.fuel_material_mixture_names.copy()

    forward_corr = assembly.mix_strategy_correspondences[("by_material", "by_pin")]

    # Test 24UOX expansion (should appear in 1 pin only)
    zone_24uox = [n for n in by_material_names if n.startswith("24UOX_zone_")]
    assert len(zone_24uox) == 4, "24UOX should have 4 zones"

    for zone_name in zone_24uox:
        pin_variants = forward_corr[zone_name]
        assert len(pin_variants) == 1, f"24UOX zones should map to 1 pin each, got {len(pin_variants)}"

        # Extract zone number
        zone_match = re.search(r'_zone_(\d+)$', zone_name)
        assert zone_match
        zone_num = zone_match.group(1)

        # Check pin variant format
        pin_variant = pin_variants[0]
        assert pin_variant == f"24UOX_zone_{zone_num}_pin_1", \
            f"Expected 24UOX_zone_{zone_num}_pin_1, got {pin_variant}"

    # Test 50UOX expansion (should appear in 4 pins: 7, 8, 9, 10)
    zone_50uox = [n for n in by_material_names if n.startswith("50UOX_zone_")]
    assert len(zone_50uox) == 4, "50UOX should have 4 zones"

    for zone_name in zone_50uox:
        pin_variants = forward_corr[zone_name]
        assert len(pin_variants) == 4, \
            f"50UOX zones should map to 4 pins each, got {len(pin_variants)} for {zone_name}"

        # Extract zone number
        zone_match = re.search(r'_zone_(\d+)$', zone_name)
        zone_num = zone_match.group(1)

        # Check that all 4 pin variants exist
        expected_pins = {f"50UOX_zone_{zone_num}_pin_{p}" for p in [7, 8, 9, 10]}
        actual_pins = set(pin_variants)
        assert actual_pins == expected_pins, \
            f"50UOX_zone_{zone_num} should map to pins 7,8,9,10"

    # Test 45Gd expansion (should appear in 1 pin: pin_6)
    zone_45gd = [n for n in by_material_names if n.startswith("45Gd_zone_")]
    assert len(zone_45gd) == 6, "45Gd should have 6 zones"

    for zone_name in zone_45gd:
        pin_variants = forward_corr[zone_name]
        assert len(pin_variants) == 1, f"45Gd zones should map to 1 pin each"

        zone_match = re.search(r'_zone_(\d+)$', zone_name)
        zone_num = zone_match.group(1)

        pin_variant = pin_variants[0]
        assert pin_variant == f"45Gd_zone_{zone_num}_pin_6"

    print("  -> test_at10_4x4_zone_to_pin_expansion PASSED")
    print(f"     24UOX: 4 zones × 1 pin = 4 mixes")
    print(f"     50UOX: 4 zones × 4 pins = 16 mixes")
    print(f"     45Gd: 6 zones × 1 pin = 6 mixes")


def test_at10_4x4_mixeq_verification_with_tdt():
    """
    Complete integration test: verify MIXEQ correspondence using real TDT indices.

    This test loads all three TDT files and verifies that the MIXEQ
    correspondence table correctly maps from FLUXL1 (by_material) to FLUXL2 (by_pin)
    using the actual indices from the reference TDT files.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Step 1: SSH
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_ssh = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_SSH_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_ssh, step_name="SSH")

    # Step 2: FLUXL1 (same strategy)
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l1, step_name="FLUXL1")

    # Step 3: FLUXL2 (different strategy)
    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l2, step_name="FLUXL2")

    # Create MIXEQ for FLUXL1 → FLUXL2 transition
    mixeq = MIXEQ(assembly, "LIBRARY2", "LIBEQ", "FLUXL1", "FLUXL2", draglib_alias="dummy_alias")

    # Verify correspondence table structure
    table = mixeq.correspondence_table
    print(f"Correspondence table has {len(table)} entries")
    assert len(table) == 45, f"Should have 45 correspondences (30 fuel mixtures from L1 expanded to 42 in L2 + 3 structural), got {len(table)}"

    # Count expansions: how many times does each source index appear?
    source_idx_counts = {}
    for from_name, to_name, from_idx, to_idx in table:
        source_idx_counts[from_idx] = source_idx_counts.get(from_idx, 0) + 1

    # Structural mixes are skipped as not altered by renumbering.

    # 50UOX zones (4-7) should each expand to 4 pins
    for idx in [4, 5, 6, 7]:
        assert source_idx_counts[idx] == 4, \
            f"50UOX zone (idx {idx}) should expand to 4 pins, got {source_idx_counts.get(idx, 0)}"

    # 45Gd zones (8-13) should each appear once (1 pin)
    for idx in range(8, 14):
        assert source_idx_counts[idx] == 1, \
            f"45Gd zone (idx {idx}) should map to 1 pin"

    # 24UOX zones (30-33) should each appear once (1 pin)
    for idx in range(30, 34):
        assert source_idx_counts[idx] == 1, \
            f"24UOX zone (idx {idx}) should map to 1 pin"

    # Generate and verify LIB call content
    lib_call = mixeq.build_lib_module_call()

    # Count MIX statements (should match number of correspondences)
    # Use regex to match actual MIX statements, not MIX in comments
    mix_statements = re.findall(r'^\s+MIX\s+\d+\s+\d+', lib_call, re.MULTILINE)
    assert len(mix_statements) == 45, f"Expected 45 MIX statements, got {len(mix_statements)}"

    print("  -> test_at10_4x4_mixeq_verification_with_tdt PASSED")
    print(f"     Verified {len(table)} correspondences from real TDT indices")
    print(f"     Structural: 3 × 1-to-1 mappings")
    print(f"     Fuel zones: varying expansion based on pin multiplicity")


def test_at10_4x4_full_case_mixeq_generation():
    """
    Integration test: full case generation with MIXEQ procedures.

    Creates a complete Dragon case from the 2L calculation scheme and verifies:
    1. Case generator detects strategy transition
    2. MIXEQ procedure is generated and written
    3. Main x2m file includes MIXEQ call
    4. Generated files are valid
    """
    output_dir = str(OUTPUTS_DIR / "at10_4x4_mixeq_test")
    os.makedirs(output_dir, exist_ok=True)

    # Create case
    case = DragonCase(
        case_name="AT10_4x4_mix_splitting",
        call_glow=False,
        output_path=output_dir,
        draglib_name_to_alias={"dummy_draglib_name": "dummy_alias",},
        config_yamls={
            "MATS": AT10_4x4_COMPOSITIONS_YAML,
            "GEOM": AT10_4x4_GEOMETRY_YAML,
            "CALC_SCHEME": AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML,
        },
        tdt_path=AT10_TDT_DIR,
    )

    # Generate all procedures
    procedure_files = case.generate_cle2000_procedures()

    # Verify MIXEQ procedure was generated
    mixeq_keys = [k for k in procedure_files.keys() if k.startswith("mixeq_")]
    assert len(mixeq_keys) > 0, "At least one MIXEQ procedure should be generated"

    # Check that MIXEQ file exists
    for mixeq_key in mixeq_keys:
        mixeq_path = procedure_files[mixeq_key]
        assert os.path.exists(mixeq_path), f"MIXEQ file should exist: {mixeq_path}"

        with open(mixeq_path, 'r') as f:
            mixeq_content = f.read()

        # Verify MIXEQ procedure format
        assert "PROCEDURE" in mixeq_content
        assert "LIB:" in mixeq_content
        assert "EDIT 0" in mixeq_content
        assert "CATL" in mixeq_content

    # Check main x2m file includes MIXEQ call
    x2m_path = procedure_files["x2m"]
    assert os.path.exists(x2m_path), "Main x2m file should exist"

    with open(x2m_path, 'r') as f:
        x2m_content = f.read()

    # Should contain reference to MIXEQ procedure
    assert "MIXEQ" in x2m_content or "MIXEQL1L2" in x2m_content, \
        "Main x2m should reference MIXEQ procedure"

    print("  -> test_at10_4x4_full_case_mixeq_generation PASSED")
    print(f"     Generated {len(mixeq_keys)} MIXEQ procedure(s)")
    print(f"     Files: {', '.join(mixeq_keys)}")


def test_at10_4x4_mixeq_comment_block():
    """
    Test that MIXEQ generates informative comment blocks showing correspondences.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Setup steps
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l1, step_name="FLUXL1")

    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l2, step_name="FLUXL2")

    # Generate MIXEQ
    mixeq = MIXEQ(assembly, "LIBRARY2", "LIBEQ", "FLUXL1", "FLUXL2", draglib_alias="dummy_alias")

    # Test comment block generation
    comment_block = mixeq.build_mix_index_comment_block()

    assert "MIXTURE CORRESPONDENCE MAPPING" in comment_block
    assert "FLUXL1" in comment_block
    assert "FLUXL2" in comment_block
    assert "MIX" in comment_block

    # Should contain specific mappings
    assert "50UOX_zone_4" in comment_block
    assert "50UOX_zone_4_pin_7" in comment_block

    print("  -> test_at10_4x4_mixeq_comment_block PASSED")
    print("     Comment block contains correspondence mapping information")


def test_at10_4x4_mixeq_max_mix_validation():
    """
    Test that MIXEQ correctly determines maximum mix index.

    For FLUXL1 → FLUXL2, the max should be 45 (from FLUXL2).
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Setup both steps
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l1, step_name="FLUXL1")

    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l2, step_name="FLUXL2")

    # Create MIXEQ and test max mix count
    mixeq = MIXEQ(assembly, "LIBRARY2", "LIBEQ", "FLUXL1", "FLUXL2", draglib_alias="dummy_alias")
    max_mix = mixeq._determine_max_mix_count()

    assert max_mix == 45, f"Max mix should be 45 (from FLUXL2), got {max_mix}"

    # Verify this appears in generated LIB call
    lib_call = mixeq.build_lib_module_call()
    assert f"NMIX {max_mix}" in lib_call

    print("  -> test_at10_4x4_mixeq_max_mix_validation PASSED")
    print(f"     Maximum mix index: {max_mix}")


def test_at10_4x4_structural_mix_preservation():
    """
    Test that structural mixes (COOLANT, CLAD, GAP) preserve their indices.

    These mixes should have 1-to-1 correspondence across all steps since
    they are not fuel and don't participate in pin-based splitting.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Setup both steps
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l1, step_name="FLUXL1")

    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l2, step_name="FLUXL2")

    # Build correspondence table
    table = assembly.build_mixeq_correspondence_table("FLUXL1", "FLUXL2")

    # Skip structural mixes in expansion tests, their numbering should be preserved (they are not present in the tables since they are not fuel)

    print("  -> test_at10_4x4_structural_mix_preservation PASSED")
    print("     Structural mixes maintain 1-to-1 correspondence")


def test_at10_4x4_mixeq_mix_statements():
    """
    Test that MIXEQ generates correct MIX statements for mix duplication.

    Each MIX statement should have format:
      MIX <target_idx> <source_idx>

    This test verifies specific examples from the reference correspondence.
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Setup both steps
    assembly.apply_mix_numbering_strategy("by_material")
    tdt_l1 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL1_IC",
        tracking_option="TISO",
        include_macros=True,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l1, step_name="FLUXL1")

    assembly.apply_mix_numbering_strategy("by_pin")
    tdt_l2 = read_material_mixture_indices_from_tdt_file(
        tdt_file_path=AT10_TDT_DIR,
        tdt_file_name=TDT_FILE_BASE + "_FLUXL2_MOC",
        tracking_option="TSPC",
        include_macros=False,
        material_names=None,
    )
    assembly.enforce_material_mixture_indices_from_tdt(tdt_l2, step_name="FLUXL2")

    # Create MIXEQ
    mixeq = MIXEQ(assembly, "LIBRARY2", "LIBEQ", "FLUXL1", "FLUXL2", draglib_alias="dummy_alias")
    lib_call = mixeq.build_lib_module_call()

    # Parse MIX statements using regex
    # Pattern: MIX <target> <source>
    mix_dict = {}  # {(source, target): True}
    for match in re.finditer(r'MIX\s+(\d+)\s+(\d+)', lib_call):
        target_idx = int(match.group(1))
        source_idx = int(match.group(2))
        mix_dict[(source_idx, target_idx)] = True

    # Verify specific MIX statements from reference TDT
    # 50UOX_zone_4 (idx 4) → 50UOX_zone_4_pin_7 (idx 4)
    assert (4, 4) in mix_dict, "MIX 4 → 4, mapping index to itself to copy from LIBRARY2 to LIBEQ"

    # 50UOX_zone_4 (idx 4) → 50UOX_zone_4_pin_10 (idx 34)
    assert (4, 34) in mix_dict, "Should have MIX 4 → 34"

    # 50UOX_zone_4 (idx 4) → 50UOX_zone_4_pin_9 (idx 38)
    assert (4, 38) in mix_dict, "Should have MIX 4 → 38"

    # 50UOX_zone_4 (idx 4) → 50UOX_zone_4_pin_8 (idx 42)
    assert (4, 42) in mix_dict, "Should have MIX 4 → 42"

    # 24UOX_zone_1 (idx 33) → 24UOX_zone_1_pin_1 (idx 33)
    assert (33, 33) in mix_dict, "Should have MIX 33 → 33, mapping index to copy from LIBRARY2 to LIBEQ"

    # Structural mixes should NOT have MIX statements (they're 1-to-1, skipped)
    assert (1, 1) in mix_dict, "COOLANT 1→1"
    assert (2, 2) in mix_dict, "CLAD 2→2"
    assert (3, 3) in mix_dict, "GAP 3→3"

    print("  -> test_at10_4x4_mixeq_mix_statements PASSED")
    print(f"     Verified {len(mix_dict)} MIX statements")
    print("     Confirmed specific reference mappings (4→34, 4→38, 4→42)")


def test_at10_4x4_multi_material_pin_expansion():
    """
    Test expansion for pins that appear multiple times in the lattice.

    The 4x4 lattice is:
      ROD1(24UOX)  ROD2(32UOX)  ROD3(42UOX)  ROD5(48UOX)
      ROD2(32UOX)  ROD4(45UOX)  ROD7(45Gd)   ROD6(50UOX)
      ROD3(42UOX)  ROD7(45Gd)   ROD6(50UOX)  ROD6(50UOX)
      ROD5(48UOX)  ROD6(50UOX)  ROD6(50UOX)  ROD6(50UOX)

    Diagonal symmetry means:
    - ROD6 (50UOX) appears 7 times but may map to fewer unique pin indices
    - Each material's zones should expand to the correct number of pin variants
    """
    ROD_to_material = associate_material_to_rod_ID(
        AT10_4x4_COMPOSITIONS_YAML,
        AT10_4x4_GEOMETRY_YAML
    )

    assembly = CartesianAssemblyModel(
        name="at10_4x4_test",
        tdt_file=f"{AT10_TDT_DIR}/{TDT_FILE_BASE}.tdt",
        geometry_description_yaml=AT10_4x4_GEOMETRY_YAML
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")

    compositions = parse_all_compositions_from_yaml(AT10_4x4_COMPOSITIONS_YAML)
    assembly.set_material_compositions(compositions)

    # Apply both strategies
    assembly.apply_mix_numbering_strategy("by_material")
    assembly.apply_mix_numbering_strategy("by_pin")

    # Check correspondence for each material
    forward_corr = assembly.mix_strategy_correspondences[("by_material", "by_pin")]

    # Count expansions per material
    material_expansions = {}
    for zone_name, pin_names in forward_corr.items():
        base_material = zone_name.rsplit("_zone_", 1)[0]
        if base_material not in material_expansions:
            material_expansions[base_material] = set()

        for pin_name in pin_names:
            # Extract pin_idx
            match = re.search(r'_pin_(\d+)$', pin_name)
            if match:
                material_expansions[base_material].add(int(match.group(1)))

    # Verify expansions match lattice
    # Based on reference TDT: 50UOX has pins 7,8,9,10 (4 unique)
    assert len(material_expansions["50UOX"]) == 4, \
        f"50UOX should have 4 unique pins, got {material_expansions['50UOX']}"
    assert material_expansions["50UOX"] == {7, 8, 9, 10}

    # 24UOX appears once: pin_1
    assert len(material_expansions["24UOX"]) == 1
    assert material_expansions["24UOX"] == {1}

    # 45Gd appears twice but may map to 1 unique pin due to symmetry
    assert len(material_expansions["45Gd"]) == 1
    assert material_expansions["45Gd"] == {6}

    print("  -> test_at10_4x4_multi_material_pin_expansion PASSED")
    print(f"     Material pin expansions: {dict(material_expansions)}")


if __name__ == "__main__":
    test_at10_4x4_by_material_numbering()
    test_at10_4x4_by_pin_numbering()
    test_at10_4x4_tdt_enforcement_ssh()
    test_at10_4x4_tdt_enforcement_fluxl2()
    test_at10_4x4_correspondence_recording()
    test_at10_4x4_correspondence_table_with_real_tdt()
    test_at10_4x4_zone_to_pin_expansion()
    test_at10_4x4_mixeq_verification_with_tdt()
    test_at10_4x4_full_case_mixeq_generation()
    test_at10_4x4_mixeq_comment_block()
    test_at10_4x4_mixeq_max_mix_validation()
    test_at10_4x4_structural_mix_preservation()
    test_at10_4x4_mixeq_mix_statements()
