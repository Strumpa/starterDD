"""
Tests for the GE14 lower core assembly model.

These tests verify:
- Assembly model creation and lattice analysis
- Material mixture numbering strategies (by_pin)
- TDT file parsing and index enforcement
- LIB module output generation
"""
import os
import pytest

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel, DummyPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.InterfaceToDD.dragon_module_calls import LIB


# ---------------------------------------------------------------------------
# Test configuration constants
# ---------------------------------------------------------------------------
FLUX_TRACKING_OPTION = "TISO"
INCLUDE_MACROS = True
PATH_TO_YAML_COMPOSITIONS = "../data/BWRProgressionProblems/GE14/inputs/material_compositions.yaml"
PATH_TO_YAML_GEOMETRY = "../data/BWRProgressionProblems/GE14/inputs/GEOM_GE14_DOM.yaml"
PATH_TO_TDT = "../../glow_data/tdt_data"
TDT_FILE_NAME = "GE14_DOM_SSH_IC"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ge14_assembly_base():
    """
    Create a base GE14 assembly model with lattice analyzed but without
    material mixture numbering. This fixture is shared across tests.
    """
    ROD_to_material = associate_material_to_rod_ID(
        PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
    )

    assembly = CartesianAssemblyModel(
        name="GE14_assembly",
        tdt_file=f"{PATH_TO_TDT}/{TDT_FILE_NAME}.tdt",
        geometry_description_yaml=PATH_TO_YAML_GEOMETRY
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(
        fuel_temperature=900.0,
        gap_temperature=600.0,
        coolant_temperature=600.0,
        moderator_temperature=600.0,
        structural_temperature=600.0
    )
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
    return assembly


@pytest.fixture(scope="module")
def ge14_assembly_numbered(ge14_assembly_base):
    """
    GE14 assembly model with material compositions set and by_pin numbering applied.
    """
    compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)
    ge14_assembly_base.set_material_compositions(compositions)
    ge14_assembly_base.number_fuel_material_mixtures_by_pin()
    return ge14_assembly_base


@pytest.fixture(scope="module")
def ge14_assembly_with_tdt_indices(ge14_assembly_numbered):
    """
    GE14 assembly model with TDT indices enforced.
    """
    tdt_indices = read_material_mixture_indices_from_tdt_file(
        PATH_TO_TDT,
        tdt_file_name=TDT_FILE_NAME,
        tracking_option=FLUX_TRACKING_OPTION,
        include_macros=INCLUDE_MACROS,
        material_names=None
    )
    ge14_assembly_numbered.enforce_material_mixture_indices_from_tdt(tdt_indices)
    return ge14_assembly_numbered


@pytest.fixture(scope="module")
def ge14_assembly_for_lib(ge14_assembly_with_tdt_indices):
    """
    GE14 assembly model ready for LIB generation (with generating/daughter mixes identified).
    """
    ge14_assembly_with_tdt_indices.identify_generating_and_daughter_mixes()
    return ge14_assembly_with_tdt_indices


# ---------------------------------------------------------------------------
# Tests: Assembly model creation and lattice analysis
# ---------------------------------------------------------------------------
class TestAssemblyModelCreation:
    """Tests for basic assembly model creation and lattice analysis."""

    def test_lattice_description_exists(self, ge14_assembly_base):
        """Verify lattice description is parsed from YAML."""
        assert ge14_assembly_base.lattice_description is not None

    def test_lattice_structure(self, ge14_assembly_base):
        """Verify lattice is built with correct dimensions."""
        assert ge14_assembly_base.lattice is not None
        assert len(ge14_assembly_base.lattice) == 10
        assert len(ge14_assembly_base.lattice[0]) == 10

    def test_lattice_pin_assignments(self, ge14_assembly_base):
        """Verify specific pins have correct rod IDs and materials."""
        assert ge14_assembly_base.lattice[0][0].rod_ID == "ROD2"
        assert ge14_assembly_base.lattice[0][0].fuel_material_name == "UOX28"
        assert ge14_assembly_base.lattice[1][2].fuel_material_name == "UOX44Gd6"
        assert ge14_assembly_base.lattice[3][3].rod_ID == "WROD"


# ---------------------------------------------------------------------------
# Tests: Pin geometry, translation offset, and center coordinates
# ---------------------------------------------------------------------------
class TestPinGeometryAndCenters:
    """Tests for pin center coordinates and assembly translation offset computation."""

    def test_channel_box_inner_side_computed(self, ge14_assembly_base):
        """Verify channel_box_inner_side is computed correctly."""
        assembly = ge14_assembly_base
        assert assembly.channel_box_inner_side is not None
        # channel_box_inner_side = assembly_pitch - 2 * channel_box_thickness - 2 * gap_wide
        expected = assembly.assembly_pitch - 2 * assembly.channel_box_thickness - 2 * assembly.gap_wide
        assert assembly.channel_box_inner_side == pytest.approx(expected)

    def test_intra_assembly_coolant_width_computed(self, ge14_assembly_base):
        """Verify intra_assembly_coolant_width is computed correctly with /2.0 division."""
        assembly = ge14_assembly_base
        assert assembly.intra_assembly_coolant_width is not None
        n_cols = len(assembly.lattice_description[0])
        pin_pitch = assembly.pin_geometry_dict.get("pin_pitch", 0)
        # intra_assembly_coolant_width = (channel_box_inner_side - n_cols * pin_pitch) / 2.0
        expected = (assembly.channel_box_inner_side - n_cols * pin_pitch) / 2.0
        assert assembly.intra_assembly_coolant_width == pytest.approx(expected)

    def test_translation_offset_computed(self, ge14_assembly_base):
        """Verify translation_offset is computed correctly."""
        assembly = ge14_assembly_base
        assert assembly.translation_offset is not None
        # translation_offset = gap_wide + channel_box_thickness + intra_assembly_coolant_width
        expected = assembly.gap_wide + assembly.channel_box_thickness + assembly.intra_assembly_coolant_width
        assert assembly.translation_offset == pytest.approx(expected)

    def test_fuel_pins_have_center_attribute(self, ge14_assembly_base):
        """Verify all FuelPinModel objects have center coordinates set."""
        for row in ge14_assembly_base.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    assert hasattr(pin, 'center'), "FuelPinModel should have a center attribute."
                    assert hasattr(pin, 'center_x'), "FuelPinModel should have a center_x attribute."
                    assert hasattr(pin, 'center_y'), "FuelPinModel should have a center_y attribute."
                    assert pin.center is not None
                    assert isinstance(pin.center, tuple)
                    assert len(pin.center) == 2

    def test_dummy_pins_have_center_attribute(self, ge14_assembly_base):
        """Verify DummyPinModel objects (water rod placeholders) have center coordinates set."""
        found_dummy = False
        for row in ge14_assembly_base.lattice:
            for pin in row:
                if isinstance(pin, DummyPinModel):
                    found_dummy = True
                    assert hasattr(pin, 'center'), "DummyPinModel should have a center attribute."
                    assert hasattr(pin, 'center_x'), "DummyPinModel should have a center_x attribute."
                    assert hasattr(pin, 'center_y'), "DummyPinModel should have a center_y attribute."
                    assert pin.center is not None
        assert found_dummy, "Expected at least one DummyPinModel in the lattice."

    def test_pin_center_coordinates_formula(self, ge14_assembly_base):
        """Verify pin center coordinates follow the expected formula."""
        assembly = ge14_assembly_base
        pin_pitch = assembly.pin_geometry_dict.get("pin_pitch", 0)
        translation = assembly.translation_offset

        for row in assembly.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    x_index = pin.x_index
                    y_index = pin.y_index
                    # Expected: center_x = translation_offset + x_index * pin_pitch + pin_pitch / 2.0
                    expected_x = translation + x_index * pin_pitch + pin_pitch / 2.0
                    expected_y = translation + y_index * pin_pitch + pin_pitch / 2.0
                    assert pin.center_x == pytest.approx(expected_x), \
                        f"Pin at ({x_index}, {y_index}) has wrong center_x: {pin.center_x} != {expected_x}"
                    assert pin.center_y == pytest.approx(expected_y), \
                        f"Pin at ({x_index}, {y_index}) has wrong center_y: {pin.center_y} != {expected_y}"

    def test_corner_pin_centers(self, ge14_assembly_base):
        """Verify corner pin centers are within assembly bounds."""
        assembly = ge14_assembly_base
        pin_pitch = assembly.pin_geometry_dict.get("pin_pitch", 0)
        
        # Check pin at (0, 0) - should be near the translation offset
        pin_00 = assembly.lattice[0][0]
        if isinstance(pin_00, FuelPinModel):
            expected_center = assembly.translation_offset + pin_pitch / 2.0
            assert pin_00.center_x == pytest.approx(expected_center)
            assert pin_00.center_y == pytest.approx(expected_center)
            # Center should be positive and within assembly pitch
            assert pin_00.center_x > 0
            assert pin_00.center_y > 0
            assert pin_00.center_x < assembly.assembly_pitch
            assert pin_00.center_y < assembly.assembly_pitch

    def test_pin_indices_match_lattice_position(self, ge14_assembly_base):
        """Verify pin x_index and y_index match their position in the lattice."""
        for y_index, row in enumerate(ge14_assembly_base.lattice):
            for x_index, pin in enumerate(row):
                if isinstance(pin, FuelPinModel):
                    assert pin.x_index == x_index, \
                        f"Pin x_index {pin.x_index} doesn't match lattice column {x_index}"
                    assert pin.y_index == y_index, \
                        f"Pin y_index {pin.y_index} doesn't match lattice row {y_index}"


# ---------------------------------------------------------------------------
# Tests: Material mixture numbering (by_pin strategy)
# ---------------------------------------------------------------------------
class TestMixNumberingByPin:
    """Tests for number_fuel_material_mixtures_by_pin strategy."""

    def test_assembly_level_attributes_populated(self, ge14_assembly_numbered):
        """Verify assembly-level mixture attributes are populated."""
        mixture_names = ge14_assembly_numbered.get_fuel_material_mixture_names()
        mixture_indices = ge14_assembly_numbered.get_fuel_material_mixture_indices()
        mixture_objects = ge14_assembly_numbered.get_fuel_material_mixtures()

        assert len(mixture_names) > 0, "Expected at least one fuel material mixture name."
        assert len(mixture_names) == len(mixture_indices), "Names and indices lists must have the same length."
        assert len(mixture_names) == len(mixture_objects), "Names and objects lists must have the same length."

    def test_indices_sequential(self, ge14_assembly_numbered):
        """Verify mixture indices are sequential starting from 1."""
        mixture_indices = ge14_assembly_numbered.get_fuel_material_mixture_indices()
        assert mixture_indices == list(range(1, len(mixture_indices) + 1)), "Indices must be sequential starting from 1."

    def test_mixtures_have_compositions(self, ge14_assembly_numbered):
        """Verify each MaterialMixture has a non-None Composition."""
        mixture_indices = ge14_assembly_numbered.get_fuel_material_mixture_indices()
        for mix in ge14_assembly_numbered.get_fuel_material_mixtures():
            assert mix.composition is not None, f"MaterialMixture '{mix.unique_material_mixture_name}' has no Composition."
            assert mix.material_mixture_index in mixture_indices

    def test_mixture_names_follow_by_pin_convention(self, ge14_assembly_numbered):
        """Verify mixture names follow the by_pin naming convention."""
        for name in ge14_assembly_numbered.get_fuel_material_mixture_names():
            assert "_pin" in name, f"by_pin mixture name '{name}' should contain '_pin'"
            assert "_zone" in name, f"by_pin mixture name '{name}' should contain '_zone'"

    def test_pins_have_required_attributes(self, ge14_assembly_numbered):
        """Verify FuelPinModel objects have pin_idx and mixture attributes after by_pin numbering."""
        for row in ge14_assembly_numbered.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    assert hasattr(pin, 'pin_idx'), "FuelPinModel should have pin_idx after by_pin numbering."
                    assert hasattr(pin, 'fuel_material_mixture_indices'), "FuelPinModel should have fuel_material_mixture_indices."
                    assert hasattr(pin, 'fuel_material_mixture_names'), "FuelPinModel should have fuel_material_mixture_names."

    def test_diagonal_symmetry_detection(self, ge14_assembly_numbered):
        """Verify symmetric pins share pin_idx when diagonal symmetry is detected."""
        if not hasattr(ge14_assembly_numbered, 'lattice_has_diagonal_symmetry'):
            pytest.skip("Symmetry detection not available")

        symmetry_type = ge14_assembly_numbered.lattice_has_diagonal_symmetry
        if symmetry_type is None:
            pytest.skip("No diagonal symmetry detected")

        n_rows = len(ge14_assembly_numbered.lattice)
        for i in range(n_rows):
            for j in range(len(ge14_assembly_numbered.lattice[i])):
                pin = ge14_assembly_numbered.lattice[i][j]
                if isinstance(pin, FuelPinModel):
                    if symmetry_type == "anti-diagonal":
                        mirror = (n_rows - 1 - j, n_rows - 1 - i)
                    else:  # "main-diagonal"
                        mirror = (j, i)

                    if mirror != (i, j):
                        mirror_pin = ge14_assembly_numbered.lattice[mirror[0]][mirror[1]]
                        if isinstance(mirror_pin, FuelPinModel):
                            assert pin.pin_idx == mirror_pin.pin_idx, \
                                f"Mirror pins at ({i},{j}) and {mirror} should share the same pin_idx"
                            assert pin.fuel_material_mixture_indices == mirror_pin.fuel_material_mixture_indices, \
                                "Mirror pins should share the same mixture indices"

    def test_different_pins_same_material_different_indices(self, ge14_assembly_numbered):
        """Verify different pins with the same material have DIFFERENT mixture indices (key by_pin rule)."""
        uox28_pins = []
        for row in ge14_assembly_numbered.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel) and pin.fuel_material_name == "UOX28":
                    uox28_pins.append(pin)

        if len(uox28_pins) <= 1:
            pytest.skip("Not enough UOX28 pins to test")

        # Find two pins with different pin_idx (non-symmetric partners)
        for i in range(len(uox28_pins)):
            for j in range(i + 1, len(uox28_pins)):
                if uox28_pins[i].pin_idx != uox28_pins[j].pin_idx:
                    assert uox28_pins[i].fuel_material_mixture_indices != uox28_pins[j].fuel_material_mixture_indices, \
                        f"Different pins (pin_idx {uox28_pins[i].pin_idx} vs {uox28_pins[j].pin_idx}) with same material " \
                        f"should have different mixture indices (by_pin rule)"
                    return  # Found a valid pair, test passes

    def test_different_materials_different_indices(self, ge14_assembly_numbered):
        """Verify pins with different materials have different mixture indices."""
        uox28_pin = None
        gd_pin = None
        for row in ge14_assembly_numbered.lattice:
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

    def test_mixture_temperatures(self, ge14_assembly_numbered):
        """Verify mixture temperatures match the uniform fuel temperature."""
        for mix in ge14_assembly_numbered.get_fuel_material_mixtures():
            assert mix.temperature == 900.0, f"Expected fuel temperature 900.0 but got {mix.temperature}."


# ---------------------------------------------------------------------------
# Tests: TDT file parsing
# ---------------------------------------------------------------------------
class TestTDTParser:
    """Tests for TDT file parsing functionality."""

    def test_tdt_contains_fuel_mixture_names(self, ge14_assembly_numbered):
        """Verify TDT file contains all fuel mixture names."""
        tdt_indices = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=None
        )

        mixture_names = ge14_assembly_numbered.get_fuel_material_mixture_names()
        for fuel_name in mixture_names:
            assert fuel_name in tdt_indices, f"Fuel mixture name '{fuel_name}' not found in TDT file."

    def test_tdt_indices_positive_integers(self):
        """Verify all TDT indices are positive integers."""
        tdt_indices = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=None
        )

        for name, idx in tdt_indices.items():
            assert isinstance(idx, int) and idx > 0, f"Expected positive int index for '{name}', got {idx}"


# ---------------------------------------------------------------------------
# Tests: TDT index enforcement
# ---------------------------------------------------------------------------
class TestTDTIndexEnforcement:
    """Tests for enforce_material_mixture_indices_from_tdt functionality."""

    def test_mixture_indices_match_tdt(self, ge14_assembly_with_tdt_indices):
        """Verify mixture indices match TDT values after enforcement."""
        tdt_indices = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=None
        )

        mixture_names = ge14_assembly_with_tdt_indices.get_fuel_material_mixture_names()
        for name, mix in zip(mixture_names, ge14_assembly_with_tdt_indices.get_fuel_material_mixtures()):
            expected_idx = tdt_indices[name]
            assert mix.material_mixture_index == expected_idx, \
                f"MaterialMixture '{name}' index should be {expected_idx}, got {mix.material_mixture_index}"

    def test_pin_level_indices_updated(self, ge14_assembly_with_tdt_indices):
        """Verify pin-level indices are also updated after TDT enforcement."""
        tdt_indices = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=None
        )

        # Find first fuel pin
        for row in ge14_assembly_with_tdt_indices.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    first_zone_name = pin.fuel_material_mixture_names[0]
                    expected_idx = tdt_indices[first_zone_name]
                    assert pin.fuel_material_mixture_indices[0] == expected_idx, \
                        f"Pin zone 1 index should be {expected_idx}, got {pin.fuel_material_mixture_indices[0]}"
                    return

    def test_non_fuel_entries_stored(self, ge14_assembly_with_tdt_indices):
        """Verify non-fuel material entries are stored separately."""
        assert hasattr(ge14_assembly_with_tdt_indices, 'non_fuel_material_mixture_indices')
        assert "COOLANT" in ge14_assembly_with_tdt_indices.non_fuel_material_mixture_indices
        assert "CLAD" in ge14_assembly_with_tdt_indices.non_fuel_material_mixture_indices


# ---------------------------------------------------------------------------
# Tests: LIB module output generation
# ---------------------------------------------------------------------------
class TestLIBCreation:
    """Tests for LIB module .c2m file generation."""

    @pytest.fixture
    def lib_output_path(self, ge14_assembly_for_lib):
        """Generate LIB output file and return its path."""
        lib = LIB(ge14_assembly_for_lib)
        return lib.write_to_c2m("outputs", "GE14_DOM_LIB")

    @pytest.fixture
    def lib_content(self, lib_output_path):
        """Read and return LIB output file content."""
        with open(lib_output_path, 'r') as f:
            return f.read()

    def test_output_file_exists(self, lib_output_path):
        """Verify output file is created."""
        assert os.path.isfile(lib_output_path), f"Expected output file at {lib_output_path}"

    def test_output_not_empty(self, lib_content):
        """Verify output file is not empty."""
        assert len(lib_content) > 0, "Generated .c2m file should not be empty."

    def test_cle2000_structural_keywords(self, lib_content):
        """Verify CLE-2000 structural keywords are present."""
        assert "LIBRARY := LIB:" in lib_content, "Missing 'LIBRARY := LIB:' module call."
        assert "NMIX" in lib_content, "Missing NMIX keyword."
        assert "DEPL LIB:" in lib_content, "Missing DEPL directive."
        assert "MIXS LIB:" in lib_content, "Missing MIXS directive."
        assert "END: ;" in lib_content, "Missing END: statement."
        assert "QUIT ." in lib_content, "Missing QUIT statement."

    def test_procedure_header_parameters(self, lib_content):
        """Verify procedure header parameters are present."""
        assert "PARAMETER LIBRARY" in lib_content, "Missing PARAMETER declaration."
        assert "<<Library>>" in lib_content, "Missing Library placeholder."
        assert "<<ssh_method>>" in lib_content, "Missing ssh_method placeholder."
        assert "<<anis_level>>" in lib_content, "Missing anis_level placeholder."
        assert "<<tran_correc>>" in lib_content, "Missing tran_correc placeholder."
        assert "DTFUEL" in lib_content, "Missing DTFUEL temperature parameter."

    def test_generating_mix_lines(self, lib_content, ge14_assembly_for_lib):
        """Verify generating mixes have full isotopic composition."""
        for gen_mix in ge14_assembly_for_lib.generating_mixes:
            idx = gen_mix.material_mixture_index
            assert f"MIX {idx}" in lib_content, \
                f"Generating mix index {idx} not found in LIB output."

            iso_comp = gen_mix.composition.get_isotope_name_composition()
            found_isotope = any(isotope in lib_content for isotope in iso_comp)
            assert found_isotope, \
                f"No isotope from generating mix '{gen_mix.unique_material_mixture_name}' found in output."

    def test_daughter_mix_comb_lines(self, lib_content, ge14_assembly_for_lib):
        """Verify daughter mixes use COMB duplication."""
        for daughter in ge14_assembly_for_lib.daughter_mixes:
            idx = daughter.material_mixture_index
            gen_idx = daughter.generating_mix.material_mixture_index
            expected_comb = f"MIX {idx} COMB {gen_idx} 1.0"
            assert expected_comb in lib_content, \
                f"Expected daughter COMB line '{expected_comb}' not found in output."

    def test_noev_keyword_present(self, lib_content):
        """Verify NOEV keyword is present for non-fuel materials."""
        assert "NOEV" in lib_content, "Missing NOEV keyword for non-fuel materials."

    def test_non_fuel_materials_present(self, lib_content, ge14_assembly_for_lib):
        """Verify known non-fuel materials are present in output."""
        non_fuel_indices = ge14_assembly_for_lib.non_fuel_material_mixture_indices
        for mat_name, idx in non_fuel_indices.items():
            if mat_name in ge14_assembly_for_lib.composition_lookup:
                assert f"MIX {idx}" in lib_content, \
                    f"Non-fuel material '{mat_name}' (index {idx}) not found in output."

    def test_mix_index_comment_block(self, lib_content):
        """Verify mix index comment block is present."""
        assert "MIX INDEX ASSIGNMENTS" in lib_content, "Missing mix index comment block."
        assert "Fuel mixes" in lib_content, "Missing fuel mixes section in comments."

    def test_nmix_value(self, lib_content, ge14_assembly_for_lib):
        """Verify NMIX value is the maximum index."""
        all_indices = list(ge14_assembly_for_lib.fuel_material_mixture_indices)
        all_indices.extend(ge14_assembly_for_lib.non_fuel_material_mixture_indices.values())
        expected_max = max(all_indices)
        assert f"NMIX {expected_max}" in lib_content, \
            f"NMIX should be {expected_max}, not found in output."