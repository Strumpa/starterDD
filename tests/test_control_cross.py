"""
Tests for BWR control cross assembly model support.

These tests verify:
- Parsing of CONTROL_CROSS_GEOMETRY from the YAML geometry file
- ControlCrossModel creation with correct attributes
- Control cross absorber tube numbering (lumped and by_tube strategies)
- Generating/daughter mix identification for per-tube absorber mixtures
- LIB module output generation with control cross mix lines
- YAML validation rejects unknown keys in CONTROL_CROSS_GEOMETRY
"""
import os
import pytest
import yaml
import tempfile

from starterDD.MaterialProperties.material_mixture import (
    MaterialMixture,
    Composition,
    parse_all_compositions_from_yaml,
)
from starterDD.DDModel.DragonModel import (
    CartesianAssemblyModel,
    ControlCrossModel,
    FuelPinModel,
    DummyPinModel,
)
from starterDD.DDModel.helpers import associate_material_to_rod_ID
from starterDD.InterfaceToDD.dragon_module_calls import LIB

from conftest import (
    GE14_COMPOSITIONS_YAML,
    GE14_DOM_C_GEOMETRY_YAML,
    OUTPUTS_DIR,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PATH_TO_YAML_COMPOSITIONS = GE14_COMPOSITIONS_YAML
PATH_TO_YAML_GEOMETRY_WITH_CTRL = GE14_DOM_C_GEOMETRY_YAML


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ge14_ctrl_assembly():
    """
    Create a GE14 assembly model from the DOM-C YAML that contains a
    CONTROL_CROSS_GEOMETRY section.  Lattice is analysed and material
    compositions are set but fuel mix numbering is NOT yet done.
    """
    ROD_to_material = associate_material_to_rod_ID(
        PATH_TO_YAML_COMPOSITIONS, PATH_TO_YAML_GEOMETRY_WITH_CTRL
    )
    assembly = CartesianAssemblyModel(
        name="GE14_ctrl",
        tdt_file="dummy.tdt",
        geometry_description_yaml=PATH_TO_YAML_GEOMETRY_WITH_CTRL,
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

    compositions = parse_all_compositions_from_yaml(PATH_TO_YAML_COMPOSITIONS)
    assembly.set_material_compositions(compositions)
    return assembly


# ---------------------------------------------------------------------------
# Tests: YAML parsing and ControlCrossModel creation
# ---------------------------------------------------------------------------
class TestControlCrossParsing:
    """Tests for parsing the CONTROL_CROSS_GEOMETRY section from YAML."""

    def test_has_control_cross_flag(self, ge14_ctrl_assembly):
        """Assembly parsed from DOM-C YAML should have has_control_cross = True."""
        assert ge14_ctrl_assembly.has_control_cross is True

    def test_control_cross_object_created(self, ge14_ctrl_assembly):
        """A ControlCrossModel instance should be stored on the assembly."""
        assert ge14_ctrl_assembly.control_cross is not None
        assert isinstance(ge14_ctrl_assembly.control_cross, ControlCrossModel)

    def test_control_cross_center(self, ge14_ctrl_assembly):
        """The cross should sit at the north-west corner."""
        assert ge14_ctrl_assembly.control_cross.center == "north-west"

    def test_control_cross_dimensions(self, ge14_ctrl_assembly):
        """Verify geometric parameters are parsed correctly from YAML."""
        ctrl = ge14_ctrl_assembly.control_cross
        assert ctrl.number_tubes_per_wing == 21
        assert ctrl.blade_half_span == pytest.approx(12.3825)
        assert ctrl.blade_thickness == pytest.approx(0.79248)
        assert ctrl.tip_radius == pytest.approx(0.39624)
        assert ctrl.central_structure_half_span == pytest.approx(1.98501)
        assert ctrl.sheath_thickness == pytest.approx(0.14224)
        assert ctrl.absorber_tube_outer_radius == pytest.approx(0.23876)
        assert ctrl.absorber_tube_inner_radius == pytest.approx(0.17526)

    def test_control_cross_material_names(self, ge14_ctrl_assembly):
        """Verify custom material names from YAML are stored."""
        ctrl = ge14_ctrl_assembly.control_cross
        assert ctrl.absorber_material == "ABS_B4C"
        assert ctrl.sheath_material == "SHEATH_SS304"

    def test_control_cross_derived_dimensions(self, ge14_ctrl_assembly):
        """Verify derived dimensions (inner_sheath_width, wing_length)."""
        ctrl = ge14_ctrl_assembly.control_cross
        expected_inner_sheath = ctrl.blade_thickness - 2 * ctrl.sheath_thickness
        expected_wing_length = ctrl.blade_half_span - ctrl.central_structure_half_span
        assert ctrl.inner_sheath_width == pytest.approx(expected_inner_sheath)
        assert ctrl.wing_length == pytest.approx(expected_wing_length)

    def test_assembly_without_control_cross(self):
        """An assembly parsed from a YAML without CONTROL_CROSS_GEOMETRY
        should have has_control_cross = False and control_cross = None."""
        from conftest import GE14_DOM_GEOMETRY_YAML
        assembly = CartesianAssemblyModel(
            name="no_ctrl",
            tdt_file="dummy.tdt",
            geometry_description_yaml=GE14_DOM_GEOMETRY_YAML,
        )
        assert assembly.has_control_cross is False
        assert assembly.control_cross is None


# ---------------------------------------------------------------------------
# Tests: ControlCrossModel standalone
# ---------------------------------------------------------------------------
class TestControlCrossModelStandalone:
    """Unit tests for the ControlCrossModel class itself."""

    def test_valid_centers(self):
        """All four valid centers should be accepted."""
        for center in ("north-west", "north-east", "south-west", "south-east"):
            model = ControlCrossModel(
                center=center,
                number_tubes_per_wing=5,
                blade_half_span=10.0,
                blade_thickness=0.8,
                tip_radius=0.4,
                central_structure_half_span=2.0,
                sheath_thickness=0.15,
                absorber_tube_outer_radius=0.25,
                absorber_tube_inner_radius=0.18,
            )
            assert model.center == center

    def test_invalid_center_raises(self):
        """An invalid center string should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid control cross center"):
            ControlCrossModel(
                center="east",
                number_tubes_per_wing=5,
                blade_half_span=10.0,
                blade_thickness=0.8,
                tip_radius=0.4,
                central_structure_half_span=2.0,
                sheath_thickness=0.15,
                absorber_tube_outer_radius=0.25,
                absorber_tube_inner_radius=0.18,
            )

    def test_default_materials(self):
        """Default absorber and sheath material names should be 'ABS' and 'SHEATH'."""
        model = ControlCrossModel(
            center="north-west",
            number_tubes_per_wing=5,
            blade_half_span=10.0,
            blade_thickness=0.8,
            tip_radius=0.4,
            central_structure_half_span=2.0,
            sheath_thickness=0.15,
            absorber_tube_outer_radius=0.25,
            absorber_tube_inner_radius=0.18,
        )
        assert model.absorber_material == "ABS"
        assert model.sheath_material == "SHEATH"

    def test_set_temperatures(self):
        """set_temperatures should store absorber and sheath temperatures."""
        model = ControlCrossModel(
            center="south-east",
            number_tubes_per_wing=5,
            blade_half_span=10.0,
            blade_thickness=0.8,
            tip_radius=0.4,
            central_structure_half_span=2.0,
            sheath_thickness=0.15,
            absorber_tube_outer_radius=0.25,
            absorber_tube_inner_radius=0.18,
        )
        model.set_temperatures(absorber_temperature=750.0, sheath_temperature=650.0)
        assert model.absorber_temperature == 750.0
        assert model.sheath_temperature == 650.0

    def test_repr(self):
        """__repr__ should include center, tubes_per_wing, and materials."""
        model = ControlCrossModel(
            center="north-west",
            number_tubes_per_wing=21,
            blade_half_span=12.0,
            blade_thickness=0.8,
            tip_radius=0.4,
            central_structure_half_span=2.0,
            sheath_thickness=0.15,
            absorber_tube_outer_radius=0.25,
            absorber_tube_inner_radius=0.18,
            absorber_material="B4C",
            sheath_material="SS304",
        )
        r = repr(model)
        assert "north-west" in r
        assert "21" in r
        assert "B4C" in r
        assert "SS304" in r

    def test_tube_material_names_initially_none(self):
        """Before numbering, tube_material_names should be None."""
        model = ControlCrossModel(
            center="north-west",
            number_tubes_per_wing=5,
            blade_half_span=10.0,
            blade_thickness=0.8,
            tip_radius=0.4,
            central_structure_half_span=2.0,
            sheath_thickness=0.15,
            absorber_tube_outer_radius=0.25,
            absorber_tube_inner_radius=0.18,
        )
        assert model.tube_material_names is None


# ---------------------------------------------------------------------------
# Tests: Lumped numbering strategy
# ---------------------------------------------------------------------------
class TestControlCrossLumpedNumbering:
    """Tests for number_control_cross_absorber_tubes(strategy='lumped')."""

    def test_lumped_sets_no_tube_names(self, ge14_ctrl_assembly):
        """Lumped strategy should leave tube_material_names as None."""
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="lumped")
        assert ge14_ctrl_assembly.control_cross.tube_material_names is None

    def test_lumped_creates_no_mixtures(self, ge14_ctrl_assembly):
        """Lumped strategy should not create any MaterialMixture objects."""
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="lumped")
        assert ge14_ctrl_assembly.control_cross_absorber_mixtures == []

    def test_lumped_no_generating_mixes_needed(self, ge14_ctrl_assembly):
        """After lumped numbering, calling identify_generating_and_daughter_control_cross_mixes
        should raise because there are no per-tube mixtures."""
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="lumped")
        with pytest.raises(RuntimeError, match="Control cross absorber mixtures have not been created"):
            ge14_ctrl_assembly.identify_generating_and_daughter_control_cross_mixes()


# ---------------------------------------------------------------------------
# Tests: By-tube numbering strategy
# ---------------------------------------------------------------------------
class TestControlCrossByTubeNumbering:
    """Tests for number_control_cross_absorber_tubes(strategy='by_tube')."""

    @pytest.fixture(autouse=True)
    def _setup_by_tube(self, ge14_ctrl_assembly):
        """Number tubes by_tube before each test in this class."""
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="by_tube")
        self.assembly = ge14_ctrl_assembly

    def test_unique_tube_names_created(self):
        """by_tube should create one unique name per tube in a wing."""
        ctrl = self.assembly.control_cross
        n = ctrl.number_tubes_per_wing
        assert ctrl.tube_material_names is not None
        assert len(ctrl.tube_material_names["wing_1"]) == n
        assert len(ctrl.tube_material_names["wing_2"]) == n

    def test_tube_naming_convention(self):
        """Unique names should follow <absorber_material>_tube_<k> pattern."""
        ctrl = self.assembly.control_cross
        base = ctrl.absorber_material
        for k in range(1, ctrl.number_tubes_per_wing + 1):
            expected_name = f"{base}_tube_{k}"
            assert expected_name in ctrl.tube_material_names["wing_1"]

    def test_wing2_is_reversed_wing1(self):
        """Wing 2 tube names should be wing 1 in reversed order (symmetry)."""
        ctrl = self.assembly.control_cross
        assert ctrl.tube_material_names["wing_2"] == list(
            reversed(ctrl.tube_material_names["wing_1"])
        )

    def test_material_mixture_objects_created(self):
        """by_tube should create MaterialMixture objects for each unique tube."""
        n = self.assembly.control_cross.number_tubes_per_wing
        mixtures = self.assembly.control_cross_absorber_mixtures
        assert len(mixtures) == n

    def test_mixture_names_match_tube_names(self):
        """Each MaterialMixture.unique_material_mixture_name should match a tube name."""
        ctrl = self.assembly.control_cross
        unique_names = set(ctrl.tube_material_names["wing_1"])
        mixture_names = {m.unique_material_mixture_name for m in self.assembly.control_cross_absorber_mixtures}
        assert mixture_names == unique_names

    def test_mixture_base_material_name(self):
        """All tube MaterialMixture.material_name should be the base absorber material."""
        base = self.assembly.control_cross.absorber_material
        for mix in self.assembly.control_cross_absorber_mixtures:
            assert mix.material_name == base

    def test_mixture_composition_is_absorber(self):
        """Each tube MaterialMixture should reference the absorber Composition."""
        base = self.assembly.control_cross.absorber_material
        expected_composition = self.assembly.composition_lookup[base]
        for mix in self.assembly.control_cross_absorber_mixtures:
            assert mix.composition is expected_composition

    def test_mixture_depletable_flag(self):
        """Depletable flag on tube mixtures should match absorber composition."""
        base = self.assembly.control_cross.absorber_material
        expected = self.assembly.composition_lookup[base].depletable
        for mix in self.assembly.control_cross_absorber_mixtures:
            assert mix.isdepletable == expected

    def test_temporary_indices_high_offset(self):
        """Temporary indices should start from 9001 to avoid collision with fuel."""
        for i, mix in enumerate(self.assembly.control_cross_absorber_mixtures):
            assert mix.material_mixture_index == 9001 + i

    def test_invalid_strategy_raises(self):
        """Passing an unknown strategy string should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid control cross numbering strategy"):
            self.assembly.number_control_cross_absorber_tubes(strategy="by_wing")


# ---------------------------------------------------------------------------
# Tests: Generating / daughter identification for per-tube mixes
# ---------------------------------------------------------------------------
class TestControlCrossGeneratingDaughterMixes:
    """Tests for identify_generating_and_daughter_control_cross_mixes."""

    @pytest.fixture(autouse=True)
    def _setup(self, ge14_ctrl_assembly):
        """Number tubes by_tube and identify generating/daughter mixes."""
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="by_tube")
        ge14_ctrl_assembly.identify_generating_and_daughter_control_cross_mixes()
        self.assembly = ge14_ctrl_assembly

    def test_one_generating_mix(self):
        """There should be exactly one generating control cross absorber mix."""
        assert len(self.assembly.control_cross_generating_mixes) == 1

    def test_daughters_count(self):
        """Number of daughters = total tubes - 1 (the generating one)."""
        n = self.assembly.control_cross.number_tubes_per_wing
        assert len(self.assembly.control_cross_daughter_mixes) == n - 1

    def test_generating_is_first_tube(self):
        """The generating mix should be the first tube (tube_1)."""
        gen = self.assembly.control_cross_generating_mixes[0]
        assert gen.is_generating is True
        assert gen.generating_mix is None
        base = self.assembly.control_cross.absorber_material
        assert gen.unique_material_mixture_name == f"{base}_tube_1"

    def test_daughters_reference_generating(self):
        """All daughter mixes should reference the single generating mix."""
        gen = self.assembly.control_cross_generating_mixes[0]
        for d in self.assembly.control_cross_daughter_mixes:
            assert d.is_generating is False
            assert d.generating_mix is gen

    def test_all_tubes_accounted_for(self):
        """Total generating + daughter should equal number_tubes_per_wing."""
        n = self.assembly.control_cross.number_tubes_per_wing
        total = (len(self.assembly.control_cross_generating_mixes) +
                 len(self.assembly.control_cross_daughter_mixes))
        assert total == n


# ---------------------------------------------------------------------------
# Tests: LIB module output with control cross
# ---------------------------------------------------------------------------
class TestLIBWithControlCross:
    """Tests for LIB class integration with control cross mixes."""

    @pytest.fixture()
    def lib_lumped(self, ge14_ctrl_assembly):
        """LIB object created from assembly with lumped control cross."""
        ge14_ctrl_assembly.number_fuel_material_mixtures_by_pin()
        ge14_ctrl_assembly.identify_generating_and_daughter_mixes()
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="lumped")
        return LIB(ge14_ctrl_assembly)

    @pytest.fixture()
    def lib_by_tube(self, ge14_ctrl_assembly):
        """LIB object created from assembly with by_tube control cross."""
        ge14_ctrl_assembly.number_fuel_material_mixtures_by_pin()
        ge14_ctrl_assembly.identify_generating_and_daughter_mixes()
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="by_tube")
        ge14_ctrl_assembly.identify_generating_and_daughter_control_cross_mixes()
        return LIB(ge14_ctrl_assembly)

    def test_lumped_no_ctrl_generating_lines(self, lib_lumped):
        """Lumped strategy: build_control_cross_generating_mix_lines should return ''."""
        assert lib_lumped.build_control_cross_generating_mix_lines() == ""

    def test_lumped_no_ctrl_daughter_lines(self, lib_lumped):
        """Lumped strategy: build_control_cross_daughter_mix_lines should return ''."""
        assert lib_lumped.build_control_cross_daughter_mix_lines() == ""

    def test_by_tube_ctrl_generating_lines_not_empty(self, lib_by_tube):
        """by_tube strategy: control cross generating lines should not be empty."""
        lines = lib_by_tube.build_control_cross_generating_mix_lines()
        assert len(lines) > 0
        assert "MIX" in lines

    def test_by_tube_ctrl_generating_lines_contain_isotopes(self, lib_by_tube):
        """by_tube generating lines should contain absorber isotopes (B10, B11, C12)."""
        lines = lib_by_tube.build_control_cross_generating_mix_lines()
        assert "B10" in lines
        assert "B11" in lines
        assert "C12" in lines or "C0" in lines  # depends on naming convention

    def test_by_tube_ctrl_daughter_lines_have_comb(self, lib_by_tube):
        """by_tube daughter lines should contain COMB keyword."""
        lines = lib_by_tube.build_control_cross_daughter_mix_lines()
        assert "COMB" in lines

    def test_by_tube_daughter_lines_count(self, lib_by_tube):
        """Number of COMB lines should equal number_tubes_per_wing - 1."""
        lines = lib_by_tube.build_control_cross_daughter_mix_lines()
        comb_count = lines.count("COMB")
        n = lib_by_tube.assembly.control_cross.number_tubes_per_wing
        assert comb_count == n - 1

    def test_lib_module_call_includes_ctrl_lines(self, lib_by_tube):
        """The full LIB module call should include control cross mix lines."""
        call = lib_by_tube.build_lib_module_call()
        assert "LIBRARY := LIB:" in call
        # Check that the control cross generating mix index appears
        gen_idx = lib_by_tube.assembly.control_cross_generating_mixes[0].material_mixture_index
        assert f"MIX {gen_idx}" in call

    def test_mix_index_comment_block_with_ctrl(self, lib_by_tube):
        """Comment block should list control cross absorber tube mixes."""
        block = lib_by_tube.build_mix_index_comment_block()
        assert "Control cross absorber tubes" in block
        base = lib_by_tube.assembly.control_cross.absorber_material
        assert f"{base}_tube_1" in block

    def test_max_mix_index_includes_ctrl(self, lib_by_tube):
        """_get_max_mix_index should account for control cross tube indices."""
        max_idx = lib_by_tube._get_max_mix_index()
        n = lib_by_tube.assembly.control_cross.number_tubes_per_wing
        # Temporary indices start at 9001
        assert max_idx >= 9000 + n

    def test_ctrl_temp_var_registered(self, lib_by_tube):
        """Custom material names (ABS_B4C, SHEATH_SS304) should be in
        non_fuel_temperature_map with 'TCTRL'."""
        assert "ABS_B4C" in lib_by_tube.non_fuel_temperature_map
        assert "SHEATH_SS304" in lib_by_tube.non_fuel_temperature_map
        assert lib_by_tube.non_fuel_temperature_map["ABS_B4C"] == "TCTRL"
        assert lib_by_tube.non_fuel_temperature_map["SHEATH_SS304"] == "TCTRL"


# ---------------------------------------------------------------------------
# Tests: YAML validation for CONTROL_CROSS_GEOMETRY
# ---------------------------------------------------------------------------
class TestControlCrossYAMLValidation:
    """Tests for YAML key validation in CONTROL_CROSS_GEOMETRY section."""

    def _write_yaml(self, tmp_path, ctrl_cross_data):
        """Helper: write a minimal valid YAML with a custom CONTROL_CROSS_GEOMETRY."""
        data = {
            "PIN_GEOMETRY": {
                "fuel_radius": 0.438,
                "pin_pitch": 1.3,
            },
            "ASSEMBLY_GEOMETRY": {
                "lattice_description": [["ROD1"]],
                "assembly_pitch": 15.0,
            },
        }
        if ctrl_cross_data is not None:
            data["CONTROL_CROSS_GEOMETRY"] = ctrl_cross_data

        path = os.path.join(tmp_path, "test_geom.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return path

    def test_unrecognised_key_raises(self, tmp_path):
        """An unrecognised key in CONTROL_CROSS_GEOMETRY should raise ValueError."""
        ctrl_data = {
            "center": "north-west",
            "number_tubes_per_wing": 5,
            "blade_half_span": 10.0,
            "blade_thickness": 0.8,
            "tip_radius": 0.4,
            "central_structure_half_span": 2.0,
            "sheath_thickness": 0.15,
            "absorber_tube_outer_radius": 0.25,
            "absorber_tube_inner_radius": 0.18,
            "bogus_key": 42,  # invalid
        }
        path = self._write_yaml(str(tmp_path), ctrl_data)
        with pytest.raises(ValueError, match="Unrecognised key"):
            CartesianAssemblyModel("test", "dummy.tdt", path)

    def test_missing_required_key_raises(self, tmp_path):
        """Missing a required key should raise ValueError."""
        ctrl_data = {
            "center": "north-west",
            # "number_tubes_per_wing" is missing
            "blade_half_span": 10.0,
            "blade_thickness": 0.8,
            "tip_radius": 0.4,
            "central_structure_half_span": 2.0,
            "sheath_thickness": 0.15,
            "absorber_tube_outer_radius": 0.25,
            "absorber_tube_inner_radius": 0.18,
        }
        path = self._write_yaml(str(tmp_path), ctrl_data)
        with pytest.raises(ValueError, match="Missing required key"):
            CartesianAssemblyModel("test", "dummy.tdt", path)


# ---------------------------------------------------------------------------
# Tests: no control cross guard
# ---------------------------------------------------------------------------
class TestNoControlCrossGuards:
    """Tests that methods raise when called on an assembly without a control cross."""

    def test_number_absorber_tubes_raises_without_cross(self):
        """number_control_cross_absorber_tubes should raise if no cross defined."""
        from conftest import GE14_DOM_GEOMETRY_YAML
        assembly = CartesianAssemblyModel(
            name="no_ctrl",
            tdt_file="dummy.tdt",
            geometry_description_yaml=GE14_DOM_GEOMETRY_YAML,
        )
        with pytest.raises(RuntimeError, match="No control cross defined"):
            assembly.number_control_cross_absorber_tubes(strategy="lumped")

    def test_number_absorber_tubes_by_tube_requires_compositions(self, ge14_ctrl_assembly):
        """by_tube strategy should raise if compositions are not set.
        (We test with a fresh assembly that has no compositions.)"""
        assembly = CartesianAssemblyModel(
            name="fresh_ctrl",
            tdt_file="dummy.tdt",
            geometry_description_yaml=PATH_TO_YAML_GEOMETRY_WITH_CTRL,
        )
        # has_control_cross is True but compositions not set
        with pytest.raises(RuntimeError, match="Material compositions have not been set"):
            assembly.number_control_cross_absorber_tubes(strategy="by_tube")


# ---------------------------------------------------------------------------
# Tests: enforce_material_mixture_indices_from_tdt with control cross
# ---------------------------------------------------------------------------
class TestEnforceTDTWithControlCross:
    """Tests for TDT index enforcement when per-tube control cross mixes exist."""

    @pytest.fixture()
    def assembly_with_by_tube(self, ge14_ctrl_assembly):
        """Assembly with by_pin fuel numbering + by_tube control cross."""
        ge14_ctrl_assembly.number_fuel_material_mixtures_by_pin()
        ge14_ctrl_assembly.number_control_cross_absorber_tubes(strategy="by_tube")
        ge14_ctrl_assembly.identify_generating_and_daughter_control_cross_mixes()
        return ge14_ctrl_assembly

    def test_ctrl_tube_names_excluded_from_non_fuel(self, assembly_with_by_tube):
        """When enforcing TDT indices, per-tube absorber names should NOT
        be placed in non_fuel_material_mixture_indices."""
        # Build a fake TDT mapping: fuel names + non-fuel + ctrl tube names
        tdt = {}
        for i, name in enumerate(assembly_with_by_tube.fuel_material_mixture_names, start=1):
            tdt[name] = i
        # Add usual non-fuel materials
        offset = len(tdt)
        for nf_name in ("CLAD", "COOLANT", "MODERATOR", "GAP", "CHANNEL_BOX"):
            offset += 1
            tdt[nf_name] = offset
        # Add ctrl tube names
        for mix in assembly_with_by_tube.control_cross_absorber_mixtures:
            offset += 1
            tdt[mix.unique_material_mixture_name] = offset

        assembly_with_by_tube.enforce_material_mixture_indices_from_tdt(tdt)

        # Verify: ctrl tube names should NOT be in non_fuel dict
        ctrl_names = {m.unique_material_mixture_name for m in assembly_with_by_tube.control_cross_absorber_mixtures}
        for name in ctrl_names:
            assert name not in assembly_with_by_tube.non_fuel_material_mixture_indices

    def test_ctrl_tube_indices_updated(self, assembly_with_by_tube):
        """After TDT enforcement, control cross mixture indices should be updated."""
        tdt = {}
        for i, name in enumerate(assembly_with_by_tube.fuel_material_mixture_names, start=1):
            tdt[name] = i
        offset = len(tdt)
        for nf_name in ("CLAD", "COOLANT", "MODERATOR", "GAP", "CHANNEL_BOX"):
            offset += 1
            tdt[nf_name] = offset
        # Assign specific indices to ctrl tubes
        expected_ctrl_indices = {}
        for mix in assembly_with_by_tube.control_cross_absorber_mixtures:
            offset += 1
            tdt[mix.unique_material_mixture_name] = offset
            expected_ctrl_indices[mix.unique_material_mixture_name] = offset

        assembly_with_by_tube.enforce_material_mixture_indices_from_tdt(tdt)

        for mix in assembly_with_by_tube.control_cross_absorber_mixtures:
            assert mix.material_mixture_index == expected_ctrl_indices[mix.unique_material_mixture_name]
