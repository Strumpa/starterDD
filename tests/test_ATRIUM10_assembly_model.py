"""
Tests for the ATRIUM-10 assembly model with a single square water rod.

These tests verify:
- Assembly model creation with square water rod geometry
- Material mixture numbering (by_pin strategy)
- Diagonal symmetry detection for the ATRIUM-10 lattice
- Water rod model attributes (SquareWaterRodModel)
- Square water rod corner_radius support
- Calculation scheme loading with 'splits' discretization for square water rods
- LIB module output generation
- EDI / COMPO module output generation
- Serpent2 model generation from the same assembly model
"""
import os
import tempfile
import pytest
import yaml

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.DDModel.DragonModel import (
    CartesianAssemblyModel, FuelPinModel, DummyPinModel, SquareWaterRodModel,
)
from starterDD.DDModel.DragonCalculationScheme import (
    DragonCalculationScheme, CalculationStep, SectorConfig,
)
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.InterfaceToDD.dragon_module_calls import LIB, EDI, COMPO, EDI_COMPO
from starterDD.InterfaceToDD.serpent2_cards import (
    Serpent2Model, S2_Settings, S2_Material, S2_PinUniverse, S2_Lattice,
    S2_ChannelGeometry, S2_EnergyGrid, S2_Detector,
    REACTION_TO_MT_NUMBER, _reaction_name_to_mt,
)


# ---------------------------------------------------------------------------
# Test configuration constants
# ---------------------------------------------------------------------------
PATH_TO_YAML_COMPOSITIONS = "../data/ATRIUM10/inputs/material_compositions.yaml"
PATH_TO_YAML_GEOMETRY = "../data/ATRIUM10/inputs/GEOM_ATRIUM10.yaml"
PATH_TO_YAML_CALC_SCHEME = "../data/ATRIUM10/inputs/CALC_SCHEME_AT10.yaml"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def at10_assembly_base():
    """
    Create a base ATRIUM-10 assembly model with lattice analyzed but without
    material mixture numbering. Uses the square water rod geometry.
    """
    ROD_to_material = associate_material_to_rod_ID(
        PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
    )

    assembly = CartesianAssemblyModel(
        name="ATRIUM10_assembly",
        tdt_file="dummy.tdt",
        geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
    )
    assembly.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly.set_uniform_temperatures(
        fuel_temperature=900.0,
        gap_temperature=600.0,
        coolant_temperature=600.0,
        moderator_temperature=600.0,
        structural_temperature=600.0,
    )
    assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
    return assembly


@pytest.fixture(scope="module")
def at10_assembly_numbered(at10_assembly_base):
    """
    ATRIUM-10 assembly model with material compositions set and by_pin
    numbering applied.
    """
    compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)
    at10_assembly_base.set_material_compositions(compositions)
    at10_assembly_base.number_fuel_material_mixtures_by_pin()
    return at10_assembly_base


@pytest.fixture(scope="module")
def at10_assembly_for_lib(at10_assembly_numbered):
    """
    ATRIUM-10 assembly with generating / daughter mixes identified and
    non-fuel mixture indices set manually (no TDT file for ATRIUM-10).
    """
    at10_assembly_numbered.identify_generating_and_daughter_mixes()

    # Since we don't have a TDT file for ATRIUM-10, create non-fuel indices
    # manually.  Use indices starting after the max fuel index.
    max_fuel_idx = max(at10_assembly_numbered.fuel_material_mixture_indices)
    at10_assembly_numbered.non_fuel_material_mixture_indices = {
        "COOLANT": max_fuel_idx + 1,
        "CLAD": max_fuel_idx + 2,
        "GAP": max_fuel_idx + 3,
        "MODERATOR": max_fuel_idx + 4,
        "CHANNEL_BOX": max_fuel_idx + 5,
    }
    return at10_assembly_numbered


# ---------------------------------------------------------------------------
# Tests: Assembly model creation and lattice analysis
# ---------------------------------------------------------------------------
class TestAT10AssemblyCreation:
    """Tests for basic ATRIUM-10 assembly creation and lattice analysis."""

    def test_lattice_description_exists(self, at10_assembly_base):
        """Verify lattice description is parsed from YAML."""
        assert at10_assembly_base.lattice_description is not None

    def test_lattice_structure(self, at10_assembly_base):
        """Verify lattice is built with correct 10x10 dimensions."""
        assert at10_assembly_base.lattice is not None
        assert len(at10_assembly_base.lattice) == 10
        assert len(at10_assembly_base.lattice[0]) == 10

    def test_lattice_corner_pin(self, at10_assembly_base):
        """Verify corner pin has correct rod ID and material."""
        pin_00 = at10_assembly_base.lattice[0][0]
        assert isinstance(pin_00, FuelPinModel)
        assert pin_00.rod_ID == "ROD1"
        assert pin_00.fuel_material_name == "24UOX"

    def test_lattice_gd_pin(self, at10_assembly_base):
        """Verify a Gd-bearing pin is correctly identified."""
        # ROD7 is in Gd_rod_ids; row 0, col 2 is ROD3 (not Gd)
        # row 1, col 2 is ROD7 (Gd)
        pin = at10_assembly_base.lattice[1][2]
        assert isinstance(pin, FuelPinModel)
        assert pin.rod_ID == "ROD7"
        assert pin.isGd is True

    def test_water_rod_placeholder(self, at10_assembly_base):
        """Verify WROD entries become DummyPinModel in the lattice."""
        # ATRIUM-10 lattice: row 4 col 4 should be WROD
        pin = at10_assembly_base.lattice[4][4]
        assert isinstance(pin, DummyPinModel)
        assert pin.rod_ID == "WROD"

    def test_number_of_water_rods(self, at10_assembly_base):
        """Verify a single square water rod is created."""
        assert len(at10_assembly_base.water_rods) == 1

    def test_water_rod_is_square(self, at10_assembly_base):
        """Verify the water rod model is SquareWaterRodModel."""
        wr = at10_assembly_base.water_rods[0]
        assert isinstance(wr, SquareWaterRodModel)

    def test_water_rod_dimensions(self, at10_assembly_base):
        """Verify the square water rod inner/outer side dimensions."""
        wr = at10_assembly_base.water_rods[0]
        assert wr.moderator_box_inner_side == pytest.approx(3.34)
        assert wr.moderator_box_outer_side == pytest.approx(3.5)

    def test_water_rod_center(self, at10_assembly_base):
        """Verify the water rod center coordinates are stored."""
        wr = at10_assembly_base.water_rods[0]
        assert wr.center is not None
        # From YAML: centers: [[7.1225, 7.1225]]
        assert wr.center[0] == pytest.approx(7.1225)
        assert wr.center[1] == pytest.approx(7.1225)

    def test_water_rod_materials(self, at10_assembly_base):
        """Verify water rod material assignments."""
        wr = at10_assembly_base.water_rods[0]
        assert wr.moderator_material_name == "MODERATOR"
        assert wr.cladding_material_name == "CLAD"
        assert wr.coolant_material_name == "COOLANT"

    def test_assembly_pitch(self, at10_assembly_base):
        """Verify assembly pitch is correctly parsed."""
        assert at10_assembly_base.assembly_pitch == pytest.approx(15.24)

    def test_channel_box_thickness(self, at10_assembly_base):
        """Verify channel box thickness."""
        assert at10_assembly_base.channel_box_thickness == pytest.approx(0.17)

    def test_pin_pitch(self, at10_assembly_base):
        """Verify pin pitch from geometry dict."""
        assert at10_assembly_base.pin_geometry_dict["pin_pitch"] == pytest.approx(1.295)


# ---------------------------------------------------------------------------
# Tests: Pin geometry and center coordinates
# ---------------------------------------------------------------------------
class TestAT10PinGeometry:
    """Tests for pin center coordinates in the ATRIUM-10 assembly."""

    def test_fuel_pins_have_centers(self, at10_assembly_base):
        """Verify all FuelPinModel objects have center coordinates set."""
        for row in at10_assembly_base.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    assert pin.center is not None
                    assert isinstance(pin.center, tuple)
                    assert len(pin.center) == 2

    def test_dummy_pins_have_centers(self, at10_assembly_base):
        """Verify DummyPinModel objects have center coordinates set."""
        found_dummy = False
        for row in at10_assembly_base.lattice:
            for pin in row:
                if isinstance(pin, DummyPinModel):
                    found_dummy = True
                    assert pin.center is not None
        assert found_dummy

    def test_pin_count(self, at10_assembly_base):
        """Count fuel pins: 100 lattice positions - 9 WROD placeholders = 91 fuel pins."""
        n_fuel = 0
        n_dummy = 0
        for row in at10_assembly_base.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    n_fuel += 1
                elif isinstance(pin, DummyPinModel):
                    n_dummy += 1
        assert n_fuel + n_dummy == 100
        assert n_dummy == 9, f"Expected 9 WROD placeholders, got {n_dummy}"
        assert n_fuel == 91, f"Expected 91 fuel pins, got {n_fuel}"


# ---------------------------------------------------------------------------
# Tests: Material mixture numbering (by_pin strategy)
# ---------------------------------------------------------------------------
class TestAT10MixNumberingByPin:
    """Tests for number_fuel_material_mixtures_by_pin on the ATRIUM-10 assembly."""

    def test_assembly_level_attributes(self, at10_assembly_numbered):
        """Verify mixture attributes are populated after numbering."""
        names = at10_assembly_numbered.get_fuel_material_mixture_names()
        indices = at10_assembly_numbered.get_fuel_material_mixture_indices()
        objects = at10_assembly_numbered.get_fuel_material_mixtures()

        assert len(names) > 0
        assert len(names) == len(indices)
        assert len(names) == len(objects)

    def test_indices_are_sequential(self, at10_assembly_numbered):
        """Verify mixture indices are sequential starting from 1."""
        indices = at10_assembly_numbered.get_fuel_material_mixture_indices()
        assert indices == list(range(1, len(indices) + 1))

    def test_mixture_names_follow_by_pin_convention(self, at10_assembly_numbered):
        """Verify names follow <material>_zone<N>_pin<M> convention."""
        for name in at10_assembly_numbered.get_fuel_material_mixture_names():
            assert "_zone" in name, f"Name '{name}' should contain '_zone'"
            assert "_pin" in name, f"Name '{name}' should contain '_pin'"

    def test_mixtures_have_compositions(self, at10_assembly_numbered):
        """Verify each MaterialMixture has a Composition."""
        for mix in at10_assembly_numbered.get_fuel_material_mixtures():
            assert mix.composition is not None

    def test_pins_have_pin_idx(self, at10_assembly_numbered):
        """Verify each FuelPinModel has a pin_idx after numbering."""
        for row in at10_assembly_numbered.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    assert hasattr(pin, 'pin_idx')
                    assert pin.pin_idx >= 1

    def test_symmetric_pins_share_indices(self, at10_assembly_numbered):
        """If diagonal symmetry is detected, verify symmetric pins share pin_idx."""
        symmetry = at10_assembly_numbered.lattice_has_diagonal_symmetry
        if symmetry is None:
            pytest.skip("No diagonal symmetry detected for ATRIUM-10")

        n = len(at10_assembly_numbered.lattice)
        for i in range(n):
            for j in range(len(at10_assembly_numbered.lattice[i])):
                pin = at10_assembly_numbered.lattice[i][j]
                if isinstance(pin, FuelPinModel):
                    if symmetry == "anti-diagonal":
                        mirror = (n - 1 - j, n - 1 - i)
                    else:
                        mirror = (j, i)
                    if mirror != (i, j):
                        mirror_pin = at10_assembly_numbered.lattice[mirror[0]][mirror[1]]
                        if isinstance(mirror_pin, FuelPinModel):
                            assert pin.pin_idx == mirror_pin.pin_idx

    def test_mixture_temperatures(self, at10_assembly_numbered):
        """Verify fuel temperatures."""
        for mix in at10_assembly_numbered.get_fuel_material_mixtures():
            assert mix.temperature == 900.0

    def test_all_names_unique(self, at10_assembly_numbered):
        """Verify all mixture names are unique."""
        names = at10_assembly_numbered.get_fuel_material_mixture_names()
        assert len(set(names)) == len(names)


# ---------------------------------------------------------------------------
# Tests: Generating and daughter mixes
# ---------------------------------------------------------------------------
class TestAT10GeneratingDaughterMixes:
    """Tests for generating/daughter mix identification on ATRIUM-10."""

    def test_generating_mixes_exist(self, at10_assembly_for_lib):
        """Verify generating mixes are identified."""
        assert len(at10_assembly_for_lib.generating_mixes) > 0

    def test_daughter_mixes_exist(self, at10_assembly_for_lib):
        """With by_pin numbering and 8 unique materials, daughters should exist."""
        assert len(at10_assembly_for_lib.daughter_mixes) > 0

    def test_generating_count_matches_materials(self, at10_assembly_for_lib):
        """Number of generating mixes should equal number of unique fuel materials."""
        unique_materials = set()
        for row in at10_assembly_for_lib.lattice:
            for pin in row:
                if isinstance(pin, FuelPinModel):
                    unique_materials.add(pin.fuel_material_name)
        assert len(at10_assembly_for_lib.generating_mixes) == len(unique_materials)

    def test_total_mixes_consistent(self, at10_assembly_for_lib):
        """Total = generating + daughter."""
        total = len(at10_assembly_for_lib.generating_mixes) + len(at10_assembly_for_lib.daughter_mixes)
        assert total == len(at10_assembly_for_lib.fuel_material_mixtures)

    def test_daughter_points_to_valid_generating(self, at10_assembly_for_lib):
        """Every daughter mix must reference a valid generating mix of the same material."""
        for mix in at10_assembly_for_lib.daughter_mixes:
            assert mix.is_generating is False
            assert mix.generating_mix is not None
            assert mix.generating_mix.is_generating is True
            assert mix.generating_mix.material_name == mix.material_name


# ---------------------------------------------------------------------------
# Tests: LIB module output generation
# ---------------------------------------------------------------------------
class TestAT10LIBCreation:
    """Tests for LIB module .c2m file generation for ATRIUM-10."""

    @pytest.fixture
    def lib_output_path(self, at10_assembly_for_lib):
        """Generate LIB output file and return its path."""
        lib = LIB(at10_assembly_for_lib)
        return lib.write_to_c2m("outputs", "ATRIUM10_LIB")

    @pytest.fixture
    def lib_content(self, lib_output_path):
        """Read and return LIB output file content."""
        with open(lib_output_path, 'r') as f:
            return f.read()

    def test_output_file_exists(self, lib_output_path):
        """Verify output file is created."""
        assert os.path.isfile(lib_output_path)

    def test_output_not_empty(self, lib_content):
        """Verify output file is not empty."""
        assert len(lib_content) > 0

    def test_cle2000_keywords(self, lib_content):
        """Verify CLE-2000 structural keywords are present."""
        assert "LIBRARY := LIB:" in lib_content
        assert "NMIX" in lib_content
        assert "DEPL LIB:" in lib_content
        assert "MIXS LIB:" in lib_content
        assert "END: ;" in lib_content
        assert "QUIT ." in lib_content

    def test_procedure_header(self, lib_content):
        """Verify procedure header parameters."""
        assert "PARAMETER LIBRARY" in lib_content

    def test_generating_mix_lines(self, lib_content, at10_assembly_for_lib):
        """Verify generating mixes have full isotopic composition."""
        for gen_mix in at10_assembly_for_lib.generating_mixes:
            idx = gen_mix.material_mixture_index
            assert f"MIX {idx}" in lib_content

    def test_daughter_mix_comb_lines(self, lib_content, at10_assembly_for_lib):
        """Verify daughter mixes use COMB duplication."""
        for daughter in at10_assembly_for_lib.daughter_mixes:
            idx = daughter.material_mixture_index
            gen_idx = daughter.generating_mix.material_mixture_index
            expected_comb = f"MIX {idx} COMB {gen_idx} 1.0"
            assert expected_comb in lib_content

    def test_mix_index_comment_block(self, lib_content):
        """Verify mix index comment block is present."""
        assert "MIX INDEX ASSIGNMENTS" in lib_content
        assert "Fuel mixes" in lib_content


# ---------------------------------------------------------------------------
# Tests: EDI / COMPO module output generation
# ---------------------------------------------------------------------------
class TestAT10EDICOMPOCreation:
    """Tests for EDI_COMPO orchestrator .c2m file generation for ATRIUM-10."""

    @pytest.fixture
    def edi_compo_output_path(self, at10_assembly_for_lib):
        """Generate EDI_COMPO output file and return its path."""
        edi_compo = EDI_COMPO(at10_assembly_for_lib)

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

        return edi_compo.write_to_c2m("outputs", "EDICPO_ATRIUM10")

    @pytest.fixture
    def edi_compo_content(self, edi_compo_output_path):
        """Read and return EDI_COMPO output file content."""
        with open(edi_compo_output_path, 'r') as f:
            return f.read()

    def test_output_file_exists(self, edi_compo_output_path):
        """Verify output file is created."""
        assert os.path.isfile(edi_compo_output_path)

    def test_output_not_empty(self, edi_compo_content):
        """Verify output file is not empty."""
        assert len(edi_compo_content) > 0

    def test_procedure_header(self, edi_compo_content):
        """Verify CLE-2000 procedure header."""
        assert "PARAMETER COMPO FLUX LIBRARY2 TRACK" in edi_compo_content
        assert "MODULE EDI: COMPO: DELETE: END:" in edi_compo_content

    def test_compo_initialization(self, edi_compo_content):
        """Verify COMPO: initialization with all directories."""
        assert "COMPO := COMPO: ::" in edi_compo_content
        assert "STEP UP 'EDIHOM_COND'" in edi_compo_content
        assert "STEP UP 'EDIHOM_295'" in edi_compo_content
        assert "STEP UP 'H_EDI_REGI_1g'" in edi_compo_content
        assert "STEP UP 'H_EDI_REGI_2g'" in edi_compo_content

    def test_edi_calls_present(self, edi_compo_content):
        """Verify all 4 EDI calls are present."""
        assert edi_compo_content.count("EDIRATES := EDI: FLUX LIBRARY2 TRACK") == 4
        assert "SAVE ON EDIHOM_COND" in edi_compo_content
        assert "SAVE ON EDIHOM_295" in edi_compo_content
        assert "SAVE ON H_EDI_REGI_1g" in edi_compo_content
        assert "SAVE ON H_EDI_REGI_2g" in edi_compo_content

    def test_compo_store_calls(self, edi_compo_content):
        """Verify all 4 COMPO: store calls are present."""
        assert edi_compo_content.count("COMPO := COMPO: COMPO EDIRATES LIBRARY2") == 4

    def test_edirates_cleanup(self, edi_compo_content):
        """Verify EDIRATES cleanup after each edition."""
        assert edi_compo_content.count("EDIRATES := DELETE: EDIRATES") == 4

    def test_energy_condensation_keywords(self, edi_compo_content):
        """Verify COND keywords for editions with condensation."""
        assert "COND 0.625" in edi_compo_content

    def test_save_and_end_footer(self, edi_compo_content):
        """Verify footer with save conditional and END."""
        assert "IF save_opt 'SAVE' = THEN" in edi_compo_content
        assert "END: ;" in edi_compo_content

    def test_merg_mix_vectors_present(self, edi_compo_content):
        """Verify MERG MIX keyword is present."""
        assert edi_compo_content.count("MERG MIX") == 4

    def test_micr_keyword(self, edi_compo_content):
        """Verify MICR keyword with isotope lists."""
        assert "MICR 5 U235 U238 U234 Gd155 Gd157" in edi_compo_content


# ---------------------------------------------------------------------------
# Tests: Serpent2 model generation
# ---------------------------------------------------------------------------
class TestAT10Serpent2Model:
    """Tests for Serpent2 model generation from the ATRIUM-10 assembly."""

    @pytest.fixture(scope="class")
    def at10_serpent2_model(self):
        """Build and return a Serpent2Model from the ATRIUM-10 assembly."""
        ROD_to_material = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)

        assembly = CartesianAssemblyModel(
            name="ATRIUM10_S2",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
        )
        assembly.set_rod_ID_to_material_mapping(ROD_to_material)
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0,
            gap_temperature=600.0,
            coolant_temperature=600.0,
            moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
        assembly.set_material_compositions(compositions)
        assembly.number_fuel_material_mixtures_by_pin()
        assembly.identify_generating_and_daughter_mixes()

        # Set non-fuel indices manually
        max_fuel_idx = max(assembly.fuel_material_mixture_indices)
        assembly.non_fuel_material_mixture_indices = {
            "COOLANT": max_fuel_idx + 1,
            "CLAD": max_fuel_idx + 2,
            "GAP": max_fuel_idx + 3,
            "MODERATOR": max_fuel_idx + 4,
            "CHANNEL_BOX": max_fuel_idx + 5,
        }

        settings = S2_Settings()
        settings.title = "ATRIUM-10 BWR assembly - starterDD generated"
        settings.bc = 2
        settings.neutrons_per_cycle = 10000
        settings.active_cycles = 200
        settings.inactive_cycles = 50

        model = Serpent2Model(assembly_model=assembly, settings=settings)
        model.build(
            gap_material_name="gap",
            clad_material_name="clad",
            coolant_material_name="cool",
            outer_water_material_name="cool_outer",
            channel_box_material_name="zr4",
        )
        # Add structural materials from the assembly composition lookup
        model.build_structural_materials_from_assembly(
            name_map={
                "COOLANT": "cool",
                "CLAD": "clad",
                "GAP": "gap",
                "MODERATOR": "moderator",
                "CHANNEL_BOX": "zr4",
            },
            temperature_map={
                "COOLANT": 600.0,
                "CLAD": 600.0,
                "GAP": 600.0,
                "MODERATOR": 600.0,
                "CHANNEL_BOX": 600.0,
            },
        )
        return model, assembly

    @pytest.fixture(scope="class")
    def at10_serpent2_output_path(self, at10_serpent2_model):
        """Write Serpent2 model and return the file path."""
        model, _ = at10_serpent2_model
        filepath = "outputs/ATRIUM10_assembly.serp"
        model.write(filepath)
        return filepath

    @pytest.fixture(scope="class")
    def at10_serpent2_content(self, at10_serpent2_output_path):
        """Read and return Serpent2 output file content."""
        with open(at10_serpent2_output_path, 'r') as f:
            return f.read()

    def test_output_file_exists(self, at10_serpent2_output_path):
        """Verify Serpent2 output file is created."""
        assert os.path.isfile(at10_serpent2_output_path)

    def test_output_not_empty(self, at10_serpent2_content):
        """Verify output file is not empty."""
        assert len(at10_serpent2_content) > 0

    def test_title_present(self, at10_serpent2_content):
        """Verify title is in the output."""
        assert "ATRIUM-10" in at10_serpent2_content

    def test_fuel_materials_present(self, at10_serpent2_content, at10_serpent2_model):
        """Verify fuel material cards are present."""
        _, assembly = at10_serpent2_model
        for mix in assembly.fuel_material_mixtures:
            mat_name = mix.unique_material_mixture_name
            assert f"mat {mat_name}" in at10_serpent2_content, \
                f"Fuel material '{mat_name}' not found in Serpent2 output"

    def test_pin_universes_present(self, at10_serpent2_content, at10_serpent2_model):
        """Verify pin universe cards are present."""
        model, _ = at10_serpent2_model
        for pin_univ in model.pin_universes:
            assert f"pin {pin_univ.universe_name}" in at10_serpent2_content

    def test_lattice_present(self, at10_serpent2_content):
        """Verify lattice card is present."""
        assert "lat " in at10_serpent2_content

    def test_channel_geometry_present(self, at10_serpent2_content):
        """Verify channel box and assembly boundary surfaces are present."""
        assert "surf" in at10_serpent2_content
        assert "cell" in at10_serpent2_content
        assert "sqc" in at10_serpent2_content

    def test_square_water_rod_surfaces(self, at10_serpent2_content):
        """Verify square moderator box surfaces are generated."""
        # Square water rods use sqc surface type
        assert "Moderator box" in at10_serpent2_content

    def test_boundary_conditions(self, at10_serpent2_content):
        """Verify reflective boundary conditions."""
        assert "set bc 2" in at10_serpent2_content

    def test_materials_have_temperature(self, at10_serpent2_content):
        """Verify materials have temperature specification."""
        assert "tmp 900.0" in at10_serpent2_content  # fuel
        assert "tmp 600.0" in at10_serpent2_content  # structural

    def test_burnable_flag_on_fuel(self, at10_serpent2_content, at10_serpent2_model):
        """Verify depletable fuel materials have burn flag."""
        _, assembly = at10_serpent2_model
        for mix in assembly.fuel_material_mixtures:
            if mix.isdepletable:
                mat_name = mix.unique_material_mixture_name
                # Find the mat card for this material and check for 'burn 1'
                idx = at10_serpent2_content.find(f"mat {mat_name}")
                if idx >= 0:
                    card_end = at10_serpent2_content.find("\n\n", idx)
                    card_text = at10_serpent2_content[idx:card_end] if card_end > idx else at10_serpent2_content[idx:]
                    assert "burn 1" in card_text, \
                        f"Depletable material '{mat_name}' should have 'burn 1'"
                    break  # test at least one

    def test_model_summary(self, at10_serpent2_model):
        """Verify model summary is generated."""
        model, _ = at10_serpent2_model
        summary = model.summary()
        assert "Serpent2 Model Summary" in summary
        assert "Materials" in summary
        assert "Pin universes" in summary

    def test_empty_universe_in_lattice(self, at10_serpent2_content):
        """Verify empty universe for water rod placeholders is defined."""
        assert "pin empty" in at10_serpent2_content

    def test_number_of_unique_pin_universes(self, at10_serpent2_model):
        """Verify number of unique pin universes matches expectations."""
        model, assembly = at10_serpent2_model
        # Pin universes = unique fuel pin positions + 1 (empty)
        n_unique_fuel = len(set(
            pin.pin_idx for row in assembly.lattice for pin in row
            if isinstance(pin, FuelPinModel)
        ))
        # model.pin_universes includes fuel + empty
        assert len(model.pin_universes) == n_unique_fuel + 1


# ---------------------------------------------------------------------------
# Tests: SquareWaterRodModel corner_radius attribute
# ---------------------------------------------------------------------------
class TestAT10SquareWaterRodCornerRadius:
    """Tests for the corner_radius attribute on SquareWaterRodModel."""

    def test_corner_radius_default_is_none(self, at10_assembly_base):
        """ATRIUM-10 YAML has no corner_radius — should default to None."""
        wr = at10_assembly_base.water_rods[0]
        assert isinstance(wr, SquareWaterRodModel)
        assert wr.corner_radius is None

    def test_corner_radius_stored_on_model(self):
        """Directly construct a SquareWaterRodModel with corner_radius."""
        wr = SquareWaterRodModel(
            bounding_box_side_length=4.0,
            moderator_box_outer_side=3.5,
            moderator_box_inner_side=3.34,
            center=(0.0, 0.0),
            rod_ID="WROD_TEST",
            corner_radius=0.5,
        )
        assert wr.corner_radius == pytest.approx(0.5)

    def test_corner_radius_zero_treated_as_sharp(self):
        """corner_radius=0.0 should be stored as 0.0 (sharp corners)."""
        wr = SquareWaterRodModel(
            bounding_box_side_length=4.0,
            moderator_box_outer_side=3.5,
            moderator_box_inner_side=3.34,
            center=(0.0, 0.0),
            rod_ID="WROD_TEST",
            corner_radius=0.0,
        )
        assert wr.corner_radius == 0.0

    def test_corner_radius_from_yaml(self):
        """Verify corner_radius is parsed from YAML when present."""
        geometry_data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.295,
            },
            "ASSEMBLY_GEOMETRY": {
                "lattice_type": "Cartesian",
                "lattice_description": [
                    ["ROD1", "ROD1", "WROD"],
                    ["ROD1", "WROD", "WROD"],
                    ["WROD", "WROD", "WROD"],
                ],
                "assembly_pitch": 3.885,
                "channel_box_thickness": 0.17,
                "non_fuel_rod_ids": ["WROD"],
            },
            "WATER_ROD_GEOMETRY": {
                "type": "square",
                "inner_side": 1.5,
                "outer_side": 1.7,
                "corner_radius": 0.25,
                "centers": [[2.59, 2.59]],
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            assembly = CartesianAssemblyModel(
                name="corner_radius_test",
                tdt_file="dummy.tdt",
                geometry_description_yaml=tmp_path,
            )
            assert assembly.water_box_corner_radius == pytest.approx(0.25)
        finally:
            os.unlink(tmp_path)

    def test_corner_radius_key_accepted_in_yaml_validation(self):
        """The YAML key 'corner_radius' should be accepted without error."""
        geometry_data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.295,
            },
            "ASSEMBLY_GEOMETRY": {
                "lattice_type": "Cartesian",
                "lattice_description": [["ROD1"]],
                "assembly_pitch": 1.295,
            },
            "WATER_ROD_GEOMETRY": {
                "type": "square",
                "inner_side": 1.0,
                "outer_side": 1.1,
                "corner_radius": 0.1,
                "centers": [[0.0, 0.0]],
            },
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(geometry_data, f)
            tmp_path = f.name

        try:
            # Should not raise ValueError for unrecognised key
            assembly = CartesianAssemblyModel(
                name="key_validation_test",
                tdt_file="dummy.tdt",
                geometry_description_yaml=tmp_path,
            )
        finally:
            os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Tests: ATRIUM-10 Calculation Scheme with splits for square water rods
# ---------------------------------------------------------------------------
class TestAT10CalcSchemeWithSplits:
    """Tests for loading the ATRIUM-10 calculation scheme YAML
    and verifying the 'splits' water rod discretization."""

    @pytest.fixture(scope="class")
    def at10_calc_scheme(self):
        """Load the ATRIUM-10 calculation scheme from YAML."""
        return DragonCalculationScheme.from_yaml(PATH_TO_YAML_CALC_SCHEME)

    def test_scheme_loads(self, at10_calc_scheme):
        """Verify the calc scheme loads without error."""
        assert at10_calc_scheme is not None
        assert at10_calc_scheme.name == "AT10_standard"

    def test_scheme_has_ssh_step(self, at10_calc_scheme):
        """Verify the SSH step exists."""
        ssh = at10_calc_scheme.get_step("SSH")
        assert ssh is not None
        assert ssh.step_type == "self_shielding"

    def test_ssh_water_rod_splits(self, at10_calc_scheme):
        """Verify the SSH step water_rod_sectors has splits=(3,3)."""
        ssh = at10_calc_scheme.get_step("SSH")
        wr_cfg = ssh.get_water_rod_sectorization()
        # Sectorization is disabled for SSH, so wr_cfg may be None
        # depending on whether the YAML has it under enabled: False
        if wr_cfg is not None:
            # When present under disabled sectorization, splits should
            # still be parsed
            assert wr_cfg.splits == (3, 3)

    def test_water_rod_sectors_splits_attribute(self, at10_calc_scheme):
        """Verify water_rod_sectors.splits is a tuple from YAML parsing."""
        for step in at10_calc_scheme.steps:
            wr_cfg = step.get_water_rod_sectorization()
            if wr_cfg is not None and wr_cfg.splits is not None:
                assert isinstance(wr_cfg.splits, tuple)
                assert len(wr_cfg.splits) == 2
                assert all(isinstance(x, int) for x in wr_cfg.splits)

    def test_water_rod_sectors_splits_values(self, at10_calc_scheme):
        """Verify the splits values match the YAML (3, 3)."""
        # Find the step that has water rod splits defined
        found = False
        for step in at10_calc_scheme.steps:
            wr_cfg = step.get_water_rod_sectorization()
            if wr_cfg is not None and wr_cfg.splits is not None:
                assert wr_cfg.splits == (3, 3)
                found = True
        assert found, "No step with water_rod splits found in CALC_SCHEME_AT10.yaml"


# ---------------------------------------------------------------------------
# Tests: Serpent2 Detector API (reaction_isotope_map)
# ---------------------------------------------------------------------------
class TestSerpent2DetectorAPI:
    """Tests for the updated Serpent2 detector API with reaction_isotope_map."""

    def test_reaction_name_to_mt_conversion(self):
        """Verify reaction name to MT number conversion works."""
        assert _reaction_name_to_mt('Fission') == 18
        assert _reaction_name_to_mt('absorption') == 27
        assert _reaction_name_to_mt('n,gamma') == 102
        # Case-insensitive
        assert _reaction_name_to_mt('fission') == 18
        assert _reaction_name_to_mt('ABSORPTION') == 27
        # Already an int should pass through
        assert _reaction_name_to_mt(18) == 18
        assert _reaction_name_to_mt(102) == 102

    def test_reaction_name_invalid_raises(self):
        """Verify invalid reaction name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown reaction name"):
            _reaction_name_to_mt('invalid_reaction')

    def test_energy_grid_full_range(self):
        """Verify S2_EnergyGrid.full_range() creates single-bin grid."""
        eg = S2_EnergyGrid.full_range(name="full_test")
        assert eg.name == "full_test"
        assert eg.n_groups == 1
        assert len(eg.boundaries) == 2
        assert eg.boundaries[0] < eg.boundaries[1]

    def test_energy_grid_two_group(self):
        """Verify S2_EnergyGrid.two_group() creates 2-bin grid."""
        eg = S2_EnergyGrid.two_group(name="2g_test")
        assert eg.name == "2g_test"
        assert eg.n_groups == 2
        assert len(eg.boundaries) == 3

    def test_detector_for_pin_with_reaction_isotope_map(self):
        """Verify S2_Detector.for_pin() works with reaction_isotope_map."""
        # Create a mock pin universe
        pin_univ = S2_PinUniverse(
            universe_name="UOX_test_1",
            material_names=["UOX_zone1_pin1", "UOX_zone2_pin1", "gap", "clad", "cool"],
            radii=[0.3, 0.4, 0.45, 0.51],
        )
        
        reaction_map = {
            'Fission': ['U235', 'U238'],
            'absorption': ['U235', 'U238', 'Gd155'],
        }
        
        det = S2_Detector.for_pin(
            pin_universe=pin_univ,
            reaction_isotope_map=reaction_map,
            energy_grid_name="full",
            detector_type=-4,
        )
        
        # Check detector name
        assert det.name == "det_UOX_test_1"
        
        # Check domain materials (fuel zones only)
        assert "UOX_zone1_pin1" in det.domain_materials
        assert "UOX_zone2_pin1" in det.domain_materials
        assert "gap" not in det.domain_materials
        assert "clad" not in det.domain_materials
        
        # Check responses: should have (2 reactions × isotopes per reaction)
        # Fission: 2 isotopes, absorption: 3 isotopes = 5 total responses
        assert len(det.responses) == 5
        
        # Verify Fission responses use MT=18
        fission_responses = [r for r in det.responses if r[0] == 18]
        assert len(fission_responses) == 2
        fission_isotopes = {r[1] for r in fission_responses}
        assert fission_isotopes == {'U235', 'U238'}
        
        # Verify absorption responses use MT=27
        absorption_responses = [r for r in det.responses if r[0] == 27]
        assert len(absorption_responses) == 3
        absorption_isotopes = {r[1] for r in absorption_responses}
        assert absorption_isotopes == {'U235', 'U238', 'Gd155'}

    def test_detector_format_card_with_reaction_isotope_map(self):
        """Verify detector card formatting with reaction_isotope_map."""
        pin_univ = S2_PinUniverse(
            universe_name="test_pin",
            material_names=["fuel_z1", "fuel_z2", "gap", "clad", "cool"],
            radii=[0.3, 0.4, 0.45, 0.51],
        )
        
        det = S2_Detector.for_pin(
            pin_universe=pin_univ,
            reaction_isotope_map={'Fission': ['U235']},
            energy_grid_name="full",
            detector_type=-4,
        )
        
        card = det.format_card()
        
        assert "det det_test_pin" in card
        assert "de full" in card
        assert "dt -4" in card
        assert "dm fuel_z1" in card
        assert "dm fuel_z2" in card
        assert "dr 18  U235" in card

    def test_deprecated_api_issues_warning(self):
        """Verify deprecated isotopes/reactions API issues warning."""
        pin_univ = S2_PinUniverse(
            universe_name="deprecated_test",
            material_names=["fuel_z1", "gap", "clad", "cool"],
            radii=[0.4, 0.45, 0.51],
        )
        
        with pytest.warns(DeprecationWarning, match="reaction_isotope_map"):
            det = S2_Detector.for_pin(
                pin_universe=pin_univ,
                isotopes=['U235', 'U238'],
                reactions=[18, 27],
            )
        
        # Should still work
        assert len(det.responses) == 4  # 2 isotopes × 2 reactions

    def test_serpent2_model_add_detector_config_with_map(self):
        """Verify Serpent2Model.add_detector_config() with reaction_isotope_map."""
        # Load assembly model
        ROD_to_material = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)

        assembly = CartesianAssemblyModel(
            name="detector_test",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
        )
        assembly.set_rod_ID_to_material_mapping(ROD_to_material)
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0,
            gap_temperature=600.0,
            coolant_temperature=600.0,
            moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
        assembly.set_material_compositions(compositions)
        assembly.number_fuel_material_mixtures_by_pin()

        model = Serpent2Model(assembly_model=assembly)
        model.build(
            gap_material_name="gap",
            clad_material_name="clad",
            coolant_material_name="cool",
        )
        
        # Use new API
        model.add_detector_config(
            reaction_isotope_map={
                'Fission': ['U235', 'U238'],
                'absorption': ['U235', 'U238', 'Gd155'],
            },
        )
        
        # Check energy grid is full range by default
        assert len(model.energy_grids) == 1
        assert model.energy_grids[0].n_groups == 1
        
        # Check isotope response materials are created for unique isotopes
        response_isotopes = {rm.isotope_name for rm in model.isotope_response_materials}
        assert response_isotopes == {'U235', 'U238', 'Gd155'}
        
        # Check detectors are created
        assert len(model.detectors) > 0

    def test_energy_grid_full_range_is_default(self):
        """Verify full_range energy grid is default when energy_bounds=None."""
        # Create a simple assembly model
        ROD_to_material = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)

        assembly = CartesianAssemblyModel(
            name="energy_grid_test",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
        )
        assembly.set_rod_ID_to_material_mapping(ROD_to_material)
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0,
            gap_temperature=600.0,
            coolant_temperature=600.0,
            moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
        assembly.set_material_compositions(compositions)
        assembly.number_fuel_material_mixtures_by_pin()

        model = Serpent2Model(assembly_model=assembly)
        model.build(
            gap_material_name="gap",
            clad_material_name="clad",
            coolant_material_name="cool",
        )
        
        # energy_bounds=None should create full_range grid
        model.add_detector_config(
            reaction_isotope_map={'Fission': ['U235']},
            energy_bounds=None,
        )
        
        eg = model.energy_grids[0]
        assert eg.n_groups == 1  # Single bin = full range
        assert eg.name == "full"

