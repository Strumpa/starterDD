"""Tests for internal functionalities, called and implicitly tested by the remaining test suite.

Covers gaps identified in the comprehensive review, grouped by module:
  A. MaterialProperties — get_element_symbol, zaid_to_isotope, XSData edges, HM_isotopes
  B. DDModel — count helpers, Santamarina radii, helpers
  C. GeometryAnalysis — tdt_parser material_names filtering
  D. InterfaceToDD — MAC edge cases, EDI edit_level
  E. Serpent2_exports — S2_Material unit tests, S2_PinUniverse.water_rod,
     S2_EnergyGrid presets/validation, S2_Detector.for_assembly,
     S2_IsotopeResponseMaterial, S2_Settings advanced, Serpent2Model extras,
     parse_isotope_name, isotope_name_to_zaid_str
"""
import os
import pytest
import tempfile
import warnings
import yaml
from pathlib import Path

# ---------------------------------------------------------------------------
# Imports from starterDD
# ---------------------------------------------------------------------------
from starterDD.MaterialProperties.material_mixture import (
    Composition,
    MaterialMixture,
    XSData,
    get_element_symbol,
    get_isotope_atomic_mass,
    fractions_to_iso_densities,
    parse_all_compositions_from_yaml,
    HM_isotopes,
)
from starterDD.DDModel.DragonModel import (
    CartesianAssemblyModel,
    FuelPinModel,
    CircularWaterRodModel,
    SquareWaterRodModel,
    DummyPinModel,
)
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.GeometryBuilder.helpers import computeSantamarinaradii
from starterDD.InterfaceToDD.dragon_module_calls import (
    LIB,
    MAC,
    EDI,
    COMPO,
    EDI_COMPO,
)
from starterDD.InterfaceToDD.Serpent2_exports import (
    S2_Material,
    S2_PinUniverse,
    S2_Lattice,
    S2_EnergyGrid,
    S2_Detector,
    S2_IsotopeResponseMaterial,
    S2_ThermalScattering,
    S2_Settings,
    S2_ChannelGeometry,
    Serpent2Model,
    parse_isotope_name,
    isotope_name_to_zaid_str,
    get_xs_suffix,
    SUPPORTED_ENERGY_MESHES,
)

# ---------------------------------------------------------------------------
# Paths — use conftest constants where available
# ---------------------------------------------------------------------------
from conftest import (
    TESTS_DIR, OUTPUTS_DIR, DATA_DIR,
    GE14_COMPOSITIONS_YAML, GE14_DOM_GEOMETRY_YAML,
    GE14_TDT_DIR, GE14_CALC_SCHEME_YAML,
)

PATH_TO_YAML_COMPOSITIONS = GE14_COMPOSITIONS_YAML
PATH_TO_YAML_GEOMETRY = GE14_DOM_GEOMETRY_YAML
PATH_TO_TDT = GE14_TDT_DIR
TDT_FILE_NAME = "GE14_DOM_SSH_IC"
FLUX_TRACKING_OPTION = "TISO"
INCLUDE_MACROS = True


# ═══════════════════════════════════════════════════════════════
#  A. MaterialProperties / material_mixture.py
# ═══════════════════════════════════════════════════════════════

class TestGetElementSymbol:
    """Unit tests for get_element_symbol(Z) which takes an atomic number (int)."""

    def test_uranium(self):
        """Z=92 → 'U'."""
        assert get_element_symbol(92) == "U"

    def test_gadolinium(self):
        """Z=64 → 'Gd'."""
        assert get_element_symbol(64) == "Gd"

    def test_actinides(self):
        """Test actinide elements by atomic number."""
        assert get_element_symbol(94) == "Pu"
        assert get_element_symbol(95) == "Am"
        assert get_element_symbol(96) == "Cm"

    def test_hydrogen(self):
        """Z=1 → 'H'."""
        assert get_element_symbol(1) == "H"

    def test_oxygen(self):
        """Z=8 → 'O'."""
        assert get_element_symbol(8) == "O"


class TestCompositionZaidConversion:
    """Tests for Composition.zaid_to_isotope() and ZAID-based init."""

    def test_zaid_keys_auto_converted(self):
        """zaid_to_isotope() creates isotope_name_composition attribute."""
        comp = Composition("test_zaid", {"92235": 1.0e-3, "92238": 2.0e-2})
        comp.zaid_to_isotope()
        # The result is stored in isotope_name_composition (not modifying isotopic_composition)
        result = comp.isotope_name_composition
        assert len(result) > 0
        for key in result:
            assert not key.isdigit(), f"Key '{key}' should have been converted from ZAID"

    def test_mixed_keys_not_converted(self):
        """When keys are isotope names (not all digits), no conversion."""
        comp = Composition("test_names", {"U235": 1.0e-3, "U238": 2.0e-2})
        assert "U235" in comp.isotopic_composition
        assert "U238" in comp.isotopic_composition


class TestCompositionGetIsotopeNameComposition:
    """Tests for Composition.get_isotope_name_composition()."""

    def test_returns_isotope_name_dict(self):
        """Verify get_isotope_name_composition returns correct dict."""
        comp = Composition("test", {"U235": 1.0e-3, "O16": 4.0e-2})
        result = comp.get_isotope_name_composition()
        assert isinstance(result, dict)
        assert "U235" in result or any("U" in k and "235" in k for k in result)


class TestHMIsotopes:
    """Tests for HM_isotopes constant."""

    def test_contains_major_actinides(self):
        """Verify HM_isotopes includes key actinides."""
        assert "U235" in HM_isotopes
        assert "U238" in HM_isotopes
        assert "Pu239" in HM_isotopes

    def test_no_duplicates(self):
        """HM_isotopes list should have no duplicates."""
        assert len(HM_isotopes) == len(set(HM_isotopes))

    def test_all_strings(self):
        """All entries should be strings."""
        assert all(isinstance(iso, str) for iso in HM_isotopes)


class TestXSDataEdgeCases:
    """Tests for XSData edge cases."""

    def test_check_completeness_missing_field(self):
        """XSData init should raise ValueError when required keys are missing."""
        with pytest.raises(ValueError):
            XSData(1, "macroscopic", {"total": [1.0]})  # missing "scattering"

    def test_set_and_get_roundtrip(self):
        """Set then get should return same values."""
        xs = XSData(1, "macroscopic", {"total": [1.0, 2.0], "scattering": [[0.5, 0.1]]})
        xs.set_cross_sections("total", [3.0, 4.0])
        result = xs.get_cross_sections("total")
        assert result == [3.0, 4.0]

    def test_get_nonexistent_returns_none(self):
        """Getting a field that was never set returns None."""
        xs = XSData(1, "macroscopic", {"total": [1.0], "scattering": [[0.5]]})
        result = xs.get_cross_sections("absorption")
        assert result is None


# ═══════════════════════════════════════════════════════════════
#  B. DDModel
# ═══════════════════════════════════════════════════════════════

class TestSantamarinaRadii:
    """Tests for computeSantamarinaradii() in GeometryBuilder/helpers.py."""

    def test_uox_returns_list(self):
        """Non-gadolinium mode should return a list of radii."""
        radii = computeSantamarinaradii(fuel_radius=0.438, gap_radius=0.45, clad_radius=0.51)
        assert isinstance(radii, list)
        assert len(radii) > 0

    def test_gd_returns_more_zones(self):
        """Gadolinium mode should return more radial zones than non-Gd."""
        radii_uox = computeSantamarinaradii(fuel_radius=0.438, gap_radius=0.45, clad_radius=0.51, gadolinium=False)
        radii_gd = computeSantamarinaradii(fuel_radius=0.438, gap_radius=0.45, clad_radius=0.51, gadolinium=True)
        assert len(radii_gd) >= len(radii_uox)

    def test_radii_are_increasing(self):
        """All radii should be strictly increasing."""
        radii = computeSantamarinaradii(fuel_radius=0.438, gap_radius=0.45, clad_radius=0.51)
        for i in range(1, len(radii)):
            assert radii[i] > radii[i - 1]

    def test_fuel_radius_in_output(self):
        """The fuel pellet radius should appear in the radii list."""
        fuel_r = 0.438
        radii = computeSantamarinaradii(fuel_radius=fuel_r, gap_radius=0.45, clad_radius=0.51)
        assert fuel_r in radii or any(
            abs(r - fuel_r) < 1e-6 for r in radii
        ), f"Fuel radius {fuel_r} not found in radii {radii}"


class TestAssociateMaterialToRodID:
    """Tests for helpers.associate_material_to_rod_ID()."""

    def test_basic_mapping(self):
        """Verify basic rod_id → material mapping from real YAML files."""
        mapping = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        assert isinstance(mapping, dict)
        assert len(mapping) > 0

    def test_missing_rod_id_raises(self):
        """A lattice rod_id with no material should raise ValueError."""
        mat_data = {
            "MIX_COMPOSITIONS": [
                {"name": "MAT_A", "rod_id": "ROD_A"},
            ]
        }
        geom_data = {
            "lattice_description": [["ROD_A", "ROD_B"]]  # ROD_B has no material
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f_mat:
            yaml.dump(mat_data, f_mat)
            mat_path = f_mat.name
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f_geo:
            yaml.dump(geom_data, f_geo)
            geo_path = f_geo.name
        try:
            with pytest.raises(ValueError, match="No material associated"):
                associate_material_to_rod_ID(mat_path, geo_path)
        finally:
            os.unlink(mat_path)
            os.unlink(geo_path)

    def test_invalid_materials_file_raises(self):
        """A non-existent materials file should raise ValueError."""
        with pytest.raises(ValueError, match="Error reading"):
            associate_material_to_rod_ID("/nonexistent/file.yaml", PATH_TO_YAML_GEOMETRY)


class TestCountHelpers:
    """Tests for CartesianAssemblyModel count helper methods."""

    @pytest.fixture(scope="class")
    def ge14_assembly(self):
        """Build a GE14 assembly with pins."""
        rod_map = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        assembly = CartesianAssemblyModel(
            name="GE14_count_test",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
        )
        assembly.set_rod_ID_to_material_mapping(rod_map)
        assembly.analyze_lattice_description(build_pins=True)
        return assembly

    def test_count_number_of_pins(self, ge14_assembly):
        """count_number_of_pins() should return a positive integer."""
        n = ge14_assembly.count_number_of_pins()
        assert isinstance(n, int)
        assert n > 0

    def test_count_number_of_unique_fuel_materials(self, ge14_assembly):
        """count_number_of_unique_fuel_materials() should return > 0."""
        n = ge14_assembly.count_number_of_unique_fuel_materials()
        assert isinstance(n, int)
        assert n > 0

    def test_count_pins_equals_fuel_lattice_positions(self, ge14_assembly):
        """Pin count should match fuel pin positions in the lattice."""
        expected = sum(
            1 for row in ge14_assembly.lattice for pin in row
            if isinstance(pin, FuelPinModel)
        )
        assert ge14_assembly.count_number_of_pins() == expected


# ═══════════════════════════════════════════════════════════════
#  C. GeometryAnalysis / tdt_parser
# ═══════════════════════════════════════════════════════════════

class TestTDTParserFiltering:
    """Tests for read_material_mixture_indices_from_tdt_file with material_names filtering."""

    def test_filter_by_material_names(self):
        """Only requested material names should be returned."""
        # First read all to know what names exist
        all_indices = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=None,
        )
        # Pick a subset
        sample_names = list(all_indices.keys())[:3]
        filtered = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=sample_names,
        )
        assert set(filtered.keys()).issubset(set(sample_names))
        assert len(filtered) <= len(sample_names)

    def test_empty_filter_returns_empty(self):
        """Filtering by an empty list should return empty dict."""
        filtered = read_material_mixture_indices_from_tdt_file(
            PATH_TO_TDT,
            tdt_file_name=TDT_FILE_NAME,
            tracking_option=FLUX_TRACKING_OPTION,
            include_macros=INCLUDE_MACROS,
            material_names=[],
        )
        assert len(filtered) == 0


# ═══════════════════════════════════════════════════════════════
#  D. InterfaceToDD / dragon_module_calls
# ═══════════════════════════════════════════════════════════════

class TestMACEdgeCases:
    """Tests for MAC class edge cases."""

    def test_empty_mixtures_raises(self):
        """write_to_c2m with no mixtures should raise ValueError."""
        mac = MAC("TESTLIB")
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="No material mixtures"):
                mac.write_to_c2m(tmpdir, "empty_test")

    def test_scattering_matrix_output(self):
        """MAC should write SCAT lines for scattering matrix data."""
        mac = MAC("TESTLIB")
        mac.ngroup = 2

        xs = XSData(1, "macroscopic", {
            "total": [1.0, 2.0],
            "absorption": [0.1, 0.2],
            "nu*fission": [0.0, 0.0],
            "chi": [1.0, 0.0],
            "scattering": [[0.5, 0.1], [0.0, 0.8]],
        })

        mix = MaterialMixture(
            material_name="test_scat",
            material_mixture_index=1,
            composition=Composition("test_scat", {"U235": 1.0e-3}),
            temperature=900.0,
        )
        mix.xs_data = xs
        mac.add_material_mixture(mix)

        with tempfile.TemporaryDirectory() as tmpdir:
            mac.write_to_c2m(tmpdir, "scat_test")
            filepath = os.path.join(tmpdir, "scat_test.c2m")
            with open(filepath, 'r') as f:
                content = f.read()
            assert "SCAT" in content
            assert "MIX 1" in content


class TestEDISetEditLevel:
    """Tests for EDI.set_edit_level() method."""

    def test_custom_edit_level_in_output(self):
        """Verify custom EDIT level appears in the EDI call."""
        # Need a minimal assembly-like object
        rod_map = associate_material_to_rod_ID(
            PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY
        )
        compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)
        assembly = CartesianAssemblyModel(
            name="edit_level_test",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY,
        )
        assembly.set_rod_ID_to_material_mapping(rod_map)
        assembly.analyze_lattice_description(build_pins=True)
        assembly.set_material_compositions(compositions)
        assembly.number_fuel_material_mixtures_by_pin()

        edi = EDI("TEST_EDI", assembly)
        edi.set_isotopes(["U235"])
        edi.set_spatial_homogenization("FUEL")
        edi.set_edit_level(3)

        call = edi.build_edi_call()
        assert "EDIT 3" in call


# ═══════════════════════════════════════════════════════════════
#  E. Serpent2_exports — Unit tests for individual S2_ classes
# ═══════════════════════════════════════════════════════════════

class TestParseIsotopeName:
    """Tests for parse_isotope_name()."""

    def test_standard_isotope(self):
        """Parse 'U235' into element and mass number."""
        result = parse_isotope_name("U235")
        assert result is not None
        # Should return something with element='U' and A=235

    def test_gadolinium_isotope(self):
        """Parse 'Gd155'."""
        result = parse_isotope_name("Gd155")
        assert result is not None

    def test_hydrogen(self):
        """Parse 'H1'."""
        result = parse_isotope_name("H1")
        assert result is not None


class TestIsotopeNameToZaidStr:
    """Tests for isotope_name_to_zaid_str()."""

    def test_uranium235(self):
        """U235 → '92235'."""
        zaid = isotope_name_to_zaid_str("U235")
        assert zaid == "92235"

    def test_hydrogen1(self):
        """H1 → '1001'."""
        zaid = isotope_name_to_zaid_str("H1")
        assert zaid == "1001"

    def test_oxygen16(self):
        """O16 → '8016'."""
        zaid = isotope_name_to_zaid_str("O16")
        assert zaid == "8016"

    def test_gadolinium155(self):
        """Gd155 → '64155'."""
        zaid = isotope_name_to_zaid_str("Gd155")
        assert zaid == "64155"

    def test_plutonium239(self):
        """Pu239 → '94239'."""
        zaid = isotope_name_to_zaid_str("Pu239")
        assert zaid == "94239"


class TestS2MaterialUnit:
    """Unit tests for S2_Material construction and formatting."""

    def test_from_raw_sum_density(self):
        """from_raw with density_type='sum' should produce 'sum' card."""
        mat = S2_Material.from_raw(
            name="test_sum",
            isotope_densities={"U235": 1.0e-3, "O16": 4.0e-2},
            temperature=900.0,
        )
        card = mat.format_card()
        assert "mat test_sum" in card
        assert "sum" in card
        assert "tmp 900.0" in card

    def test_from_mass_fractions(self):
        """from_mass_fractions should produce negative density card."""
        mat = S2_Material.from_mass_fractions(
            name="test_mass",
            mass_fractions={"92235": 0.04, "U238": 0.96},
            mass_density=10.0,
            temperature=600.0,
        )
        card = mat.format_card()
        assert "mat test_mass" in card
        assert "-10.0" in card
        assert "tmp 600.0" in card

    def test_repr(self):
        """Verify __repr__ returns informative string."""
        mat = S2_Material.from_raw(
            name="repr_test",
            isotope_densities={"U235": 1.0e-3},
            temperature=900.0,
        )
        r = repr(mat)
        assert "repr_test" in r
        assert "900" in r


class TestS2PinUniverseUnit:
    """Unit tests for S2_PinUniverse."""

    def test_water_rod_constructor(self):
        """S2_PinUniverse.water_rod() should create 3-region pin."""
        wr = S2_PinUniverse.water_rod(
            universe_name="WR_1",
            inner_radius=0.5,
            outer_radius=0.6,
        )
        assert wr.universe_name == "WR_1"
        assert len(wr.material_names) == 3
        assert len(wr.radii) == 2

    def test_empty_constructor(self):
        """S2_PinUniverse.empty() should create single-region pin."""
        e = S2_PinUniverse.empty(universe_name="void")
        assert len(e.material_names) == 1
        assert len(e.radii) == 0

    def test_get_fuel_material_names_fuel_pin(self):
        """Fuel pin with 5 regions: 2 fuel + gap + clad + cool → 2 fuel names."""
        pin = S2_PinUniverse(
            universe_name="fuel_1",
            material_names=["fuel_z1", "fuel_z2", "gap", "clad", "cool"],
            radii=[0.3, 0.4, 0.45, 0.51],
        )
        fuel_names = pin.get_fuel_material_names()
        assert fuel_names == ["fuel_z1", "fuel_z2"]

    def test_get_fuel_material_names_water_rod(self):
        """Water rod with 3 regions has no fuel materials."""
        wr = S2_PinUniverse.water_rod("WR", 0.5, 0.6)
        assert wr.get_fuel_material_names() == []

    def test_format_card(self):
        """Verify pin card formatting."""
        pin = S2_PinUniverse(
            universe_name="test_pin",
            material_names=["fuel", "gap", "clad", "cool"],
            radii=[0.4, 0.45, 0.51],
        )
        card = pin.format_card()
        assert "pin test_pin" in card
        assert "fuel" in card
        assert "cool" in card

    def test_invalid_material_count_raises(self):
        """len(materials) != len(radii)+1 should raise ValueError."""
        with pytest.raises(ValueError, match="must equal"):
            S2_PinUniverse("bad", ["a", "b"], [0.5, 0.6, 0.7])


class TestS2EnergyGridUnit:
    """Unit tests for S2_EnergyGrid."""

    def test_from_preset_26g(self):
        """26-group preset should load correctly."""
        eg = S2_EnergyGrid.from_preset("26g")
        assert eg.n_groups == 26
        assert len(eg.boundaries) == 27

    def test_from_preset_295g(self):
        """295-group preset should load correctly."""
        eg = S2_EnergyGrid.from_preset("295g")
        assert eg.n_groups == 295

    def test_from_preset_invalid_raises(self):
        """Invalid preset name should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown energy mesh preset"):
            S2_EnergyGrid.from_preset("999g")

    def test_insufficient_boundaries_raises(self):
        """Fewer than 2 boundaries should raise ValueError."""
        with pytest.raises(ValueError, match="at least 2 boundaries"):
            S2_EnergyGrid("bad", [1.0])

    def test_from_dragon_energy_mesh(self):
        """Convert from Dragon (eV, descending) to Serpent2 (MeV, ascending)."""
        # Dragon convention: high → low in eV
        dragon_bounds = [20.0e6, 0.625, 1.0e-5]
        eg = S2_EnergyGrid.from_dragon_energy_mesh("dragon_test", dragon_bounds)
        assert eg.n_groups == 2
        # Boundaries should be in MeV, ascending
        assert eg.boundaries[0] < eg.boundaries[1]
        assert eg.boundaries[0] == pytest.approx(1.0e-5 * 1.0e-6)

    def test_format_card(self):
        """Verify ene card formatting."""
        eg = S2_EnergyGrid("test_grid", [1e-11, 0.625e-6, 20.0])
        card = eg.format_card()
        assert "ene test_grid  1" in card

    def test_two_group_with_custom_cutoff(self):
        """two_group() with custom cutoff should create correct grid."""
        eg = S2_EnergyGrid.two_group(
            name="custom_2g",
            thermal_cutoff=1.0e-6,
            e_min=1.0e-11,
            e_max=20.0,
        )
        assert eg.n_groups == 2
        # Middle boundary should be the cutoff
        assert eg.boundaries[1] == pytest.approx(1.0e-6)


class TestS2DetectorForAssembly:
    """Tests for S2_Detector.for_assembly()."""

    def test_creates_assembly_detector(self):
        """Verify for_assembly creates a single detector over all fuel materials."""
        fuel_mats = ["fuel_z1_pin1", "fuel_z2_pin1", "fuel_z1_pin2"]
        det = S2_Detector.for_assembly(
            all_fuel_material_names=fuel_mats,
            reaction_isotope_map={"fission": ["U235"], "absorption": ["U235", "Gd155"]},
            energy_grid_name="26g",
        )
        assert det.name == "det_assembly_26g"
        assert len(det.domain_materials) == 3
        assert len(det.responses) == 3  # 1 fission + 2 absorption
        assert det.detector_type == -4

    def test_custom_name(self):
        """Verify custom detector name."""
        det = S2_Detector.for_assembly(
            all_fuel_material_names=["fuel"],
            reaction_isotope_map={"fission": ["U235"]},
            name="my_detector",
        )
        assert det.name == "my_detector"


class TestS2IsotopeResponseMaterial:
    """Tests for S2_IsotopeResponseMaterial."""

    def test_uranium_response_material(self):
        """U235 response material should have ZAID 92235."""
        irm = S2_IsotopeResponseMaterial("U235", temperature=900.0)
        assert irm.zaid == "92235"
        assert irm.material_name == "U235"

    def test_format_card(self):
        """Verify response material card format."""
        irm = S2_IsotopeResponseMaterial("Gd155", temperature=600.0)
        card = irm.format_card()
        assert "mat Gd155" in card
        assert "1.0" in card
        assert "tmp 600.0" in card

    def test_xs_suffix_property(self):
        """xs_suffix should return a valid suffix string."""
        irm = S2_IsotopeResponseMaterial("U238", temperature=900.0)
        suffix = irm.xs_suffix
        assert isinstance(suffix, str)
        assert suffix.startswith(".")


class TestS2SettingsAdvanced:
    """Tests for S2_Settings advanced features."""

    def test_set_endfb8r1_libraries(self):
        """Verify ENDF/B-VIII paths are set correctly."""
        settings = S2_Settings()
        settings.set_endfb8r1_libraries("/data/nuclear")
        assert settings.acelib == "/data/nuclear/endfb8r1.xsdata"
        assert settings.declib == "/data/nuclear/endfb8r1.dec"
        assert settings.nfylib == "/data/nuclear/endfb8r1.nfy"

    def test_set_jeff311_libraries(self):
        """Verify JEFF-3.1.1 paths are set correctly."""
        settings = S2_Settings()
        settings.set_jeff311_libraries("/data/nuclear")
        assert settings.acelib == "/data/nuclear/JEFF-311.xsdata"

    def test_burnup_settings_in_output(self):
        """Verify burnup settings appear in formatted output."""
        settings = S2_Settings()
        settings.power_density = 40.0
        settings.burnup_steps = [0.0, 0.1, 0.5, 1.0, 5.0]
        settings.bumode = 2
        settings.pcc = 1
        cards = settings.format_cards()
        assert "powdens" in cards
        assert "bumode 2" in cards
        assert "pcc 1" in cards
        assert "dep butot" in cards

    def test_gcu_setting(self):
        """Verify GCU setting appears in output."""
        settings = S2_Settings()
        settings.gcu = [0]
        cards = settings.format_cards()
        assert "set gcu 0" in cards

    def test_nfg_setting(self):
        """Verify NFG setting appears in output."""
        settings = S2_Settings()
        settings.nfg = [0.625e-6]
        cards = settings.format_cards()
        assert "set nfg" in cards

    def test_plot_commands(self):
        """Verify plot commands are appended."""
        settings = S2_Settings()
        settings.add_plot(plot_type=3, x_pixels=500, y_pixels=500)
        cards = settings.format_cards()
        assert "plot 3 500 500" in cards

    def test_nuclear_data_evaluation(self):
        """Verify nuclear data evaluation selection."""
        settings = S2_Settings()
        settings.set_nuclear_data_evaluation("endfb8r1")
        assert settings.nuclear_data_evaluation == "endfb8r1"

    def test_unknown_evaluation_warns(self):
        """Unknown evaluation should warn."""
        settings = S2_Settings()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            settings.set_nuclear_data_evaluation("nonexistent_eval")
            assert len(w) >= 1


class TestSerpent2ModelExtras:
    """Tests for Serpent2Model methods not covered by integration tests."""

    def test_add_structural_material(self):
        """Verify add_structural_material adds to materials list."""
        model = Serpent2Model()
        model.add_structural_material(
            name="zr4",
            isotope_densities={"Zr90": 2.0e-2, "Zr91": 5.0e-3},
            temperature=600.0,
        )
        assert len(model.materials) == 1
        assert model.materials[0].name == "zr4"

    def test_build_without_assembly_raises(self):
        """build() without assembly_model should raise RuntimeError."""
        model = Serpent2Model()
        with pytest.raises(RuntimeError, match="No assembly model"):
            model.build()

    def test_add_flux_detector(self):
        """Verify add_flux_detector creates energy grid and detector."""
        model = Serpent2Model()
        model.add_flux_detector(energy_grid_name="2g", name="flux_det")
        assert len(model.energy_grids) == 1
        assert model.energy_grids[0].n_groups == 2
        assert len(model.detectors) == 1
        assert model.detectors[0].name == "flux_det"

    def test_add_thermal_scattering_explicit(self):
        """Verify explicit library thermal scattering."""
        model = Serpent2Model()
        model.add_thermal_scattering(name="lwtr", library_name="lwtr.16t")
        assert len(model.thermal_scattering_laws) == 1
        assert model.thermal_scattering_laws[0].library_name == "lwtr.16t"

    def test_add_thermal_scattering_temperature(self):
        """Verify temperature-resolved thermal scattering."""
        model = Serpent2Model()
        model.add_thermal_scattering(name="lwtr", temperature=600.0)
        assert len(model.thermal_scattering_laws) == 1

    def test_add_thermal_scattering_dedup(self):
        """Adding same therm name twice should not duplicate."""
        model = Serpent2Model()
        model.add_thermal_scattering(name="lwtr", temperature=600.0)
        model.add_thermal_scattering(name="lwtr", temperature=600.0)
        assert len(model.thermal_scattering_laws) == 1

    def test_summary_method(self):
        """Verify summary string contains expected sections."""
        model = Serpent2Model()
        s = model.summary()
        assert "Serpent2 Model Summary" in s
        assert "Materials" in s
        assert "Pin universes" in s
        assert "Energy grids" in s


class TestS2LatticeUnit:
    """Unit tests for S2_Lattice."""

    def test_constructor_validates_rows(self):
        """Mismatched row count should raise ValueError."""
        with pytest.raises(ValueError, match="rows"):
            S2_Lattice("bad", 0.0, 0.0, 2, 3, 1.0,
                        universe_map=[["a", "b"], ["c", "d"]])

    def test_constructor_validates_columns(self):
        """Mismatched column count should raise ValueError."""
        with pytest.raises(ValueError, match="elements"):
            S2_Lattice("bad", 0.0, 0.0, 3, 2, 1.0,
                        universe_map=[["a", "b"], ["c", "d"]])

    def test_format_card(self):
        """Verify lattice card formatting."""
        lat = S2_Lattice("test_lat", 5.0, 5.0, 2, 2, 1.3,
                          universe_map=[["A", "B"], ["C", "D"]])
        card = lat.format_card()
        assert "lat test_lat  1" in card
        assert "2 2" in card


# ═══════════════════════════════════════════════════════════════
#  F. DragonCalculationScheme — summary() method
# ═══════════════════════════════════════════════════════════════

class TestDragonCalcSchemeSummary:
    """Test for DragonCalculationScheme.summary() method."""

    def test_summary_returns_string(self):
        """Verify summary() returns a non-empty string."""
        from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
        PATH_TO_CALC_SCHEME = str(
            DATA_DIR / "BWRProgressionProblems" / "GE14_inputs" / "ASSEMBLY" / "CALC_SCHEME_GE14.yaml"
        )
        scheme = DragonCalculationScheme.from_yaml(PATH_TO_CALC_SCHEME)
        s = scheme.summary()
        assert isinstance(s, str)
        assert len(s) > 0
