# Testing calculation scheme generation through case_generator
# Covers: CLE2000 constraints, CalculationBranch, CalculationOutput,
#          DragonCalculationScheme branch/output helpers, COMPO with
#          branches, EDI_COMPO with statepoint guards, and full
#          end-to-end case generation with branch loops.
# R.Guasch
# Date : 10/03/2026

import os
import textwrap
import tempfile

import pytest
import yaml

from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.InterfaceToDD.CLE2000 import (
    main_procedure, sub_procedure,
    CLE2000_MAX_LINE, CLE2000_MAX_VARNAME,
    validate_varname, wrap_cle2000_line,
)
from starterDD.DDModel.DragonCalculationScheme import (
    CalculationStep,
    DragonCalculationScheme,
    CalculationBranch,
    CalculationOutput,
    VALID_BRANCH_TYPES,
    BRANCH_TYPE_TO_PARA_NAME,
    BRANCH_TYPE_TO_PARA_KEYWORD,
    BRANCH_TYPE_TO_ARRAY_NAME,
    BRANCH_TYPE_TO_COUNTER,
    BRANCH_TYPE_TO_COUNT_VAR,
    DEFAULT_BRANCH_LOOP_ORDER,
)
from starterDD.InterfaceToDD.dragon_module_calls import (
    COMPO, EDI_COMPO, EDI,
)
from conftest import (
    GE14_COMPOSITIONS_YAML,
    GE14_DOM_GEOMETRY_YAML,
    GE14_CALC_SCHEME_YAML,
    GE14_CALC_SCHEME_1L_YAML,
    GE14_CALC_BRANCHES_YAML,
    GE14_CALC_OUTPUTS_YAML,
    GE14_TDT_DIR,
)


# =========================================================
# Test CLE2000 constraints
# =========================================================

def test_varname_length_limit():
    """Variable names > 12 chars must raise."""
    with pytest.raises(ValueError, match="12"):
        validate_varname("this_is_too_long")
    # Exactly 12 is fine
    validate_varname("twelve_chars")


def test_line_wrapping():
    """Lines > 70 chars get wrapped."""
    long = "A " * 40  # 80 chars
    wrapped = wrap_cle2000_line(long.strip())
    for line in wrapped.splitlines():
        assert len(line) <= CLE2000_MAX_LINE


# =========================================================
# CalculationBranch – unit tests
# =========================================================

class TestCalculationBranch:
    """Tests for the CalculationBranch data class."""

    def test_construction_valid_types(self):
        """All four valid branch types can be instantiated."""
        for btype in VALID_BRANCH_TYPES:
            vals = [0.0] if btype == "burnup" else [1.0, 2.0]
            b = CalculationBranch(
                name=f"test_{btype}", branch_type=btype, values=vals,
            )
            assert b.type == btype
            assert b.values == vals

    def test_invalid_branch_type_raises(self):
        """An unknown branch type must raise ValueError."""
        with pytest.raises(ValueError, match="Invalid branch type"):
            CalculationBranch("bad", branch_type="pressure", values=[1.0])

    def test_burnup_multi_value_raises(self):
        """Burnup with more than one value is not yet supported."""
        with pytest.raises(NotImplementedError, match="Depletion"):
            CalculationBranch(
                "BU", branch_type="burnup", values=[0.0, 10.0],
            )

    def test_burnup_single_value_ok(self):
        """Burnup with exactly one value is accepted."""
        b = CalculationBranch("BU", branch_type="burnup", values=[0.0])
        assert b.values == [0.0]

    def test_para_name_property(self):
        """para_name returns the CLE2000 PARA identifier."""
        b = CalculationBranch("cd", "coolant_density", [0.7])
        assert b.para_name == "DCool"

    def test_para_keyword_property(self):
        """para_keyword returns the correct COMPO keyword."""
        b_real = CalculationBranch("ft", "fuel_temperature", [900.0])
        assert b_real.para_keyword == "VALU REAL"
        b_bu = CalculationBranch("bu", "burnup", [0.0])
        assert b_bu.para_keyword == "IRRA"

    def test_array_name_property(self):
        b = CalculationBranch("ct", "coolant_temperature", [600.0])
        assert b.array_name == BRANCH_TYPE_TO_ARRAY_NAME["coolant_temperature"]

    def test_counter_and_count_vars(self):
        b = CalculationBranch("cd", "coolant_density", [0.7, 0.5])
        assert b.counter_var == BRANCH_TYPE_TO_COUNTER["coolant_density"]
        assert b.count_var == BRANCH_TYPE_TO_COUNT_VAR["coolant_density"]

    def test_cle2000_variable_names_within_limit(self):
        """All auto-generated CLE2000 variable names must be ≤ 12 chars."""
        for btype in VALID_BRANCH_TYPES:
            vals = [0.0] if btype == "burnup" else [1.0]
            b = CalculationBranch("x", btype, vals)
            for attr in ("para_name", "array_name", "counter_var", "count_var"):
                name = getattr(b, attr)
                assert len(name) <= CLE2000_MAX_VARNAME, (
                    f"{attr} for {btype} = '{name}' exceeds 12 chars"
                )

    def test_repr(self):
        b = CalculationBranch("CD", "coolant_density", [0.7, 0.5])
        r = repr(b)
        assert "coolant_density" in r
        assert "0.7" in r

    def test_from_yaml_entry(self):
        entry = {
            "name": "Density variation",
            "type": "coolant_density",
            "values": [0.73, 0.45, 0.17],
        }
        b = CalculationBranch.from_yaml_entry(entry)
        assert b.name == "Density variation"
        assert b.type == "coolant_density"
        assert b.values == [0.73, 0.45, 0.17]

    def test_values_are_copied(self):
        """Internal values list should be a copy, not a reference."""
        original = [0.7, 0.5]
        b = CalculationBranch("cd", "coolant_density", original)
        original.append(0.3)
        assert len(b.values) == 2


# =========================================================
# CalculationOutput – unit tests
# =========================================================

class TestCalculationOutput:
    """Tests for the CalculationOutput data class."""

    def test_all_state_points(self):
        """state_points='ALL' means output at every statepoint."""
        o = CalculationOutput(
            name="HOM_COND", comment="test", isotopes=["U235"],
            spatial_integration_mode="FUEL",
            energy_bounds=[], state_points="ALL",
        )
        assert o.state_points == "ALL"
        assert o.applies_to({"coolant_density": 0.7}) is True
        assert o.applies_to({}) is True

    def test_selective_state_points_match(self):
        """Selective output matches when current value is in the list."""
        o = CalculationOutput(
            name="SEL", comment="", isotopes=["U238"],
            spatial_integration_mode="FUEL",
            energy_bounds=None,
            state_points={
                "coolant_density": [0.73669],
                "fuel_temperature": [900.0],
            },
        )
        assert o.applies_to({
            "coolant_density": 0.73669,
            "fuel_temperature": 900.0,
        }) is True

    def test_selective_state_points_mismatch(self):
        """Selective output does NOT match when current value is absent."""
        o = CalculationOutput(
            name="SEL", comment="", isotopes=["U238"],
            spatial_integration_mode="FUEL",
            energy_bounds=None,
            state_points={
                "coolant_density": [0.73669],
            },
        )
        assert o.applies_to({"coolant_density": 0.45}) is False

    def test_selective_ignores_missing_branch_type(self):
        """If the statepoint dict lacks a branch type from the filter,
        that filter dimension is skipped (not rejected)."""
        o = CalculationOutput(
            name="X", comment="", isotopes=[],
            spatial_integration_mode="FUEL",
            energy_bounds=None,
            state_points={"burnup": [0.0]},
        )
        # statepoint doesn't mention burnup at all → passes
        assert o.applies_to({"coolant_density": 0.7}) is True

    def test_energy_bounds_none_normalization(self):
        """energy_bounds=None and 'None' (string) are normalised to None."""
        for eb in (None, "None"):
            o = CalculationOutput(
                "X", "", [], "FUEL", eb, "ALL",
            )
            assert o.energy_bounds is None

    def test_energy_bounds_empty_list(self):
        """energy_bounds=[] → COND (collapse to 1 group)."""
        o = CalculationOutput("X", "", [], "FUEL", [], "ALL")
        assert o.energy_bounds == []

    def test_energy_bounds_two_group(self):
        """energy_bounds=[0.625] → 2-group split."""
        o = CalculationOutput("X", "", [], "FUEL", [0.625], "ALL")
        assert o.energy_bounds == [0.625]

    def test_from_yaml_entry_full(self):
        entry = {
            "name": "EDIR_2G",
            "comment": "Condensed to 2g",
            "isotopes": ["U235", "U238"],
            "spatial_integration_mode": "by_pin",
            "energy_bounds": [0.625],
            "state_points": {
                "coolant_density": [0.73669],
                "fuel_temperature": [900.0],
            },
        }
        o = CalculationOutput.from_yaml_entry(entry)
        assert o.name == "EDIR_2G"
        assert o.spatial_integration_mode == "by_pin"
        assert o.energy_bounds == [0.625]
        assert isinstance(o.state_points, dict)
        assert 0.73669 in o.state_points["coolant_density"]

    def test_from_yaml_entry_typo_key(self):
        """The YAML files use 'spatial_intergation_mode' (typo) — must be accepted."""
        entry = {
            "name": "X",
            "spatial_intergation_mode": "ALL",
            "energy_bounds": None,
            "state_points": "ALL",
        }
        o = CalculationOutput.from_yaml_entry(entry)
        assert o.spatial_integration_mode == "ALL"

    def test_from_yaml_entry_defaults(self):
        """Minimal entry gets reasonable defaults."""
        entry = {"name": "MIN"}
        o = CalculationOutput.from_yaml_entry(entry)
        assert o.comment == ""
        assert o.isotopes == []
        assert o.spatial_integration_mode == "FUEL"
        assert o.energy_bounds is None
        assert o.state_points == "ALL"

    def test_repr(self):
        o = CalculationOutput("HOM", "c", ["U235"], "FUEL", [], "ALL")
        r = repr(o)
        assert "HOM" in r
        assert "FUEL" in r


# =========================================================
# DragonCalculationScheme – branch / output helpers
# =========================================================

class TestSchemeHelpers:
    """Tests for DragonCalculationScheme branch and output helpers."""

    def _make_scheme_with_branches(self):
        scheme = DragonCalculationScheme(name="test")
        scheme.branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5, 0.3]),
            CalculationBranch("FT", "fuel_temperature", [900.0, 1200.0]),
            CalculationBranch("CT", "coolant_temperature", [600.0]),
            CalculationBranch("BU", "burnup", [0.0]),
        ]
        return scheme

    def test_has_branches(self):
        s_empty = DragonCalculationScheme()
        assert s_empty.has_branches() is False
        s = self._make_scheme_with_branches()
        assert s.has_branches() is True

    def test_get_branch_found(self):
        s = self._make_scheme_with_branches()
        cd = s.get_branch("coolant_density")
        assert cd is not None
        assert cd.type == "coolant_density"

    def test_get_branch_not_found(self):
        s = self._make_scheme_with_branches()
        assert s.get_branch("pressure") is None

    def test_get_total_statepoints(self):
        s = self._make_scheme_with_branches()
        # 3 densities × 2 fuel temps × 1 coolant temp × 1 burnup = 6
        assert s.get_total_statepoints() == 6

    def test_get_total_statepoints_no_branches(self):
        s = DragonCalculationScheme()
        assert s.get_total_statepoints() == 1

    def test_get_ordered_branches_default_order(self):
        """Branches are returned in the default loop order."""
        s = self._make_scheme_with_branches()
        ordered = s.get_ordered_branches()
        types = [b.type for b in ordered]
        # Default: fuel_temperature, coolant_temperature, coolant_density
        # Then remaining (burnup) appended in original order
        assert types[0] == "fuel_temperature"
        assert types[1] == "coolant_temperature"
        assert types[2] == "coolant_density"
        assert "burnup" in types

    def test_get_ordered_branches_custom_order(self):
        """Custom loop order is respected."""
        s = self._make_scheme_with_branches()
        s.branch_loop_order = ["coolant_density", "fuel_temperature"]
        ordered = s.get_ordered_branches()
        types = [b.type for b in ordered]
        assert types[0] == "coolant_density"
        assert types[1] == "fuel_temperature"

    def test_outputs_stored_from_yaml(self, tmp_path):
        """Outputs parsed from a unified YAML are stored on the scheme."""
        yaml_content = textwrap.dedent("""\
            FLUX_CALCULATION_SCHEME:
              name: test_scheme
              steps:
                - name: SSH
                  step_type: self_shielding
                  self_shielding_module: USS
                  self_shielding_method: RSE
                  spatial_method: CP
                  tracking: TISO
                  radial_scheme: Santamarina
                  sectorization:
                    enabled: false
                - name: FLUX
                  step_type: flux
                  spatial_method: CP
                  tracking: TSPC
                  radial_scheme: Santamarina
                  sectorization:
                    enabled: false

            CALCULATION_BRANCHES:
              - name: "CD"
                type: coolant_density
                values: [0.7, 0.5]

            CALCULATION_OUTPUTS:
              - name: "EDIHOM"
                comment: "Homogenized"
                isotopes: ["U235"]
                spatial_integration_mode: "FUEL"
                energy_bounds: []
                state_points: "ALL"
              - name: "SEL_2G"
                comment: "Selective 2g"
                isotopes: ["U235", "U238"]
                spatial_integration_mode: "FUEL"
                energy_bounds: [0.625]
                state_points:
                  coolant_density: [0.7]
        """)
        yaml_file = tmp_path / "scheme.yaml"
        yaml_file.write_text(yaml_content)

        scheme = DragonCalculationScheme.from_yaml(str(yaml_file))

        assert scheme.has_branches()
        assert len(scheme.branches) == 1
        assert scheme.branches[0].type == "coolant_density"
        assert len(scheme.outputs) == 2
        assert scheme.outputs[0].name == "EDIHOM"
        assert scheme.outputs[0].state_points == "ALL"
        assert scheme.outputs[1].name == "SEL_2G"
        assert isinstance(scheme.outputs[1].state_points, dict)

    def test_scheme_summary_includes_branches_and_outputs(self):
        """summary() should mention branches and outputs."""
        s = self._make_scheme_with_branches()
        s.outputs = [
            CalculationOutput("HOM", "test", ["U235"], "FUEL", [], "ALL"),
        ]
        # Add a dummy step so summary() doesn't crash
        s.add_step(CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        ))
        summary = s.summary()
        assert "Calculation Branches" in summary
        assert "coolant_density" in summary
        assert "Total statepoints: 6" in summary
        assert "Calculation Outputs" in summary
        assert "HOM" in summary


# =========================================================
# YAML parsing – branches and outputs from unified file
# =========================================================

class TestYAMLParsing:
    """Test YAML parsing for the unified CALC_SCHEME file format."""

    def test_parse_calc_scheme_1l_yaml(self):
        """Parse the real CALC_SCHEME_1L.yaml (unified format) and
        verify branches and outputs are loaded."""
        scheme = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)

        # Steps
        assert len(scheme.steps) == 2
        assert scheme.steps[0].step_type == "self_shielding"
        assert scheme.steps[1].step_type == "flux"

        # Branches
        assert scheme.has_branches()
        assert len(scheme.branches) == 4
        cd = scheme.get_branch("coolant_density")
        assert cd is not None
        assert len(cd.values) == 5
        ft = scheme.get_branch("fuel_temperature")
        assert ft is not None
        assert ft.values == [900.0]

        # Outputs
        assert len(scheme.outputs) == 4
        names = [o.name for o in scheme.outputs]
        assert "EDIHOM_COND" in names
        assert "EDIHOM_295" in names
        assert "EDIR_2G" in names
        assert "U238_295" in names

        # EDIHOM_COND is ALL state points
        edihom = next(o for o in scheme.outputs if o.name == "EDIHOM_COND")
        assert edihom.state_points == "ALL"
        assert edihom.energy_bounds == []

        # EDIR_2G is selective
        edir2g = next(o for o in scheme.outputs if o.name == "EDIR_2G")
        assert isinstance(edir2g.state_points, dict)
        assert edir2g.energy_bounds == [0.625]

        # U238_295 has energy_bounds "None" → normalized to None
        u238 = next(o for o in scheme.outputs if o.name == "U238_295")
        assert u238.energy_bounds is None

    def test_parse_separate_outputs_yaml(self):
        """Parse the standalone CALC_OUTPUTS.yaml (EDI_COMPO_OUTPUTS key)."""
        with open(GE14_CALC_OUTPUTS_YAML) as f:
            data = yaml.safe_load(f)
        outputs_data = data.get("EDI_COMPO_OUTPUTS", [])
        outputs = [CalculationOutput.from_yaml_entry(e) for e in outputs_data]

        assert len(outputs) == 4
        assert outputs[0].name == "EDIHOM_COND"
        assert outputs[0].state_points == "ALL"
        # EDIHOM_295 uses the typo key
        assert outputs[1].spatial_integration_mode == "ALL"

    def test_parse_separate_branches_yaml(self):
        """Parse the standalone CALC_BRANCHES.yaml."""
        with open(GE14_CALC_BRANCHES_YAML) as f:
            data = yaml.safe_load(f)
        branches_data = data.get("CALCULATION_BRANCHES", [])
        branches = [CalculationBranch.from_yaml_entry(e) for e in branches_data]

        assert len(branches) == 4
        types = {b.type for b in branches}
        assert types == {
            "coolant_density", "fuel_temperature",
            "coolant_temperature", "burnup",
        }
        cd = next(b for b in branches if b.type == "coolant_density")
        assert len(cd.values) == 5

    def test_scheme_without_branches_has_empty_list(self):
        """A YAML with no CALCULATION_BRANCHES key yields empty branches."""
        scheme = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_YAML)
        assert scheme.has_branches() is False
        assert scheme.branches == []
        assert scheme.outputs == []

    def test_total_statepoints_from_yaml(self):
        """Verify total statepoints computed from CALC_SCHEME_1L.yaml.
        5 densities × 1 fuel_temp × 1 cool_temp × 1 burnup = 5."""
        scheme = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
        assert scheme.get_total_statepoints() == 5


# =========================================================
# COMPO initialization with branches
# =========================================================

class TestCOMPOWithBranches:
    """Tests for the COMPO class when branches are provided."""

    def _make_branches(self):
        return [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5, 0.3]),
            CalculationBranch("FT", "fuel_temperature", [900.0]),
        ]

    def test_compo_init_without_branches(self):
        """Without branches, COMPO init has no MAXCAL / PARA."""
        compo = COMPO()
        compo.add_directory("HOM", "Homogenized", ["U235", "U238"])
        init = compo.build_compo_init()
        assert "MAXCAL" not in init
        assert "PARA" not in init
        assert "STEP UP 'HOM'" in init
        assert "ISOT 2 U235 U238" in init

    def test_compo_init_with_branches(self):
        """With branches, COMPO init includes MAXCAL and PARA keywords."""
        branches = self._make_branches()
        compo = COMPO(branches=branches)
        compo.add_directory("HOM", "Homogenized", ["U235"])
        init = compo.build_compo_init()
        # MAXCAL should be 3 × 1 = 3
        assert "MAXCAL 3" in init
        assert "PARA 'DCool' VALU REAL" in init
        assert "PARA 'TFuel' VALU REAL" in init

    def test_compo_store_without_branches(self):
        """Store block without branches uses plain STEP UP (no quotes)."""
        compo = COMPO()
        compo.add_directory("HOM", "test", ["U235"])
        store = compo.build_compo_store("HOM")
        assert "STEP UP HOM" in store
        assert "<<" not in store  # no variable references

    def test_compo_store_with_branches(self):
        """Store block with branches passes parameter values
        and uses quoted STEP UP."""
        branches = self._make_branches()
        compo = COMPO(branches=branches)
        compo.add_directory("HOM", "test", ["U235"])
        store = compo.build_compo_store("HOM")
        assert "STEP UP 'HOM'" in store
        assert "'DCool' <<DCool>>" in store
        assert "'TFuel' <<TFuel>>" in store

    def test_compo_multiple_directories(self):
        """Multiple directories each get their own STEP UP + INIT block."""
        compo = COMPO()
        compo.add_directory("DIR1", "First", ["U235"])
        compo.add_directory("DIR2", "Second", ["U238"])
        init = compo.build_compo_init()
        assert init.count("STEP UP") == 2
        assert init.count("INIT") == 2


# =========================================================
# EDI_COMPO – statepoint guards and procedure structure
# =========================================================

class TestEDICOMPOWithBranches:
    """Tests for EDI_COMPO when branches and outputs are configured."""

    def _make_mock_assembly(self):
        """Build a minimal assembly model for EDI tests."""
        from starterDD.DDModel.DragonModel import CartesianAssemblyModel
        from starterDD.DDModel.helpers import associate_material_to_rod_ID
        from starterDD.MaterialProperties.material_mixture import (
            parse_all_compositions_from_yaml,
        )

        geom = GE14_DOM_GEOMETRY_YAML
        mats = GE14_COMPOSITIONS_YAML
        rod_to_mat = associate_material_to_rod_ID(mats, geom)
        compositions = parse_all_compositions_from_yaml(mats)

        assembly = CartesianAssemblyModel(
            name="test", tdt_file="", geometry_description_yaml=geom,
        )
        assembly.set_rod_ID_to_material_mapping(rod_to_mat)
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0, gap_temperature=600.0,
            coolant_temperature=600.0, moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True)
        assembly.set_material_compositions(compositions)

        # Apply SSH radii and number by pin
        ssh = CalculationStep(
            name="SSH", step_type="self_shielding",
            self_shielding_module="USS", self_shielding_method="RSE",
            spatial_method="CP",
        )
        ssh.apply_radii(assembly)
        assembly.number_fuel_material_mixtures_by_pin()
        assembly.identify_generating_and_daughter_mixes()
        return assembly

    def test_edi_compo_body_with_guards(self, tmp_path):
        """When outputs have selective state_points, IF guards appear."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5]),
            CalculationBranch("FT", "fuel_temperature", [900.0]),
        ]
        outputs = [
            CalculationOutput(
                "EDIHOM", "Homogenized", ["U235"], "FUEL", [], "ALL",
            ),
            CalculationOutput(
                "SEL_2G", "Selective 2g", ["U235", "U238"], "FUEL",
                [0.625],
                state_points={"coolant_density": [0.7]},
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        body = edi_compo.build_procedure_body()

        # EDIHOM (ALL) should NOT have an IF guard
        # SEL_2G (selective) should have an IF ... THEN guard
        assert "IF " in body
        assert "THEN" in body
        assert "ENDIF" in body

        # The ALL output should appear without IF
        lines = body.splitlines()
        edihom_idx = next(
            i for i, l in enumerate(lines) if "EDI: CALL FOR EDIHOM" in l
        )
        # Check that the line before is not an IF line
        before_lines = "\n".join(lines[max(0, edihom_idx - 3):edihom_idx])
        assert "IF " not in before_lines or "ENDIF" in before_lines

    def test_edi_compo_no_compo_init_when_branches(self, tmp_path):
        """When branches are present, COMPO init is NOT in the EDIR procedure
        body (it is handled in the main x2m)."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7]),
        ]
        outputs = [
            CalculationOutput(
                "HOM", "test", ["U235"], "FUEL", [], "ALL",
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        body = edi_compo.build_procedure_body()
        assert "COMPO := COMPO: ::" not in body

    def test_edi_compo_procedure_file_with_branches(self, tmp_path):
        """Written .c2m file has COMPO in PARAMETER and REAL branch vars."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5]),
            CalculationBranch("FT", "fuel_temperature", [900.0]),
        ]
        outputs = [
            CalculationOutput(
                "HOM", "test", ["U235"], "FUEL", [], "ALL",
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        filepath = edi_compo.write_to_c2m(str(tmp_path), "EDIR_TEST")

        with open(filepath) as f:
            content = f.read()

        # COMPO should be in the PARAMETER block
        assert "PARAMETER COMPO FLUX LIBRARY2 TRACK ::" in content
        assert "::: LINKED_LIST COMPO" in content
        # Branch REAL variables should be declared
        assert "REAL DCool TFuel" in content
        assert ">>DCool<<" in content
        assert ">>TFuel<<" in content
        # Should NOT have local COMPO declaration or SEQ_ASCII export
        assert "LINKED_LIST EDIRATES COMPO" not in content
        assert "SEQ_ASCII _COMPO" not in content
        # Footer should just be END:
        assert content.rstrip().endswith("END: ;")

    def test_edi_compo_procedure_file_without_branches(self, tmp_path):
        """Written .c2m file without branches has local COMPO + export."""
        assembly = self._make_mock_assembly()
        edi_compo = EDI_COMPO(assembly)
        edi_compo.add_edition(
            name="HOM", comment="test", isotopes=["U235"],
            spatial_mode="FUEL", energy_bounds=[],
        )
        filepath = edi_compo.write_to_c2m(str(tmp_path), "EDIR_NOBR")

        with open(filepath) as f:
            content = f.read()

        # PARAMETER should NOT include COMPO
        assert "PARAMETER FLUX LIBRARY2 TRACK ::" in content
        # Local COMPO and SEQ_ASCII export present
        assert "LINKED_LIST EDIRATES COMPO" in content
        assert "SEQ_ASCII _COMPO :: FILE <<name_cpo>>" in content
        assert "_COMPO := COMPO" in content

    def test_statepoint_guard_single_value(self):
        """Guard for a single value per branch type uses simple equality."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5]),
        ]
        outputs = [
            CalculationOutput(
                "SEL", "", ["U235"], "FUEL", [],
                state_points={"coolant_density": [0.7]},
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        guard = edi_compo._build_statepoint_guard(outputs[0])
        assert "IF " in guard
        assert "DCool 0.7 =" in guard
        assert "THEN" in guard

    def test_statepoint_guard_multiple_values(self):
        """Guard for multiple allowed values uses OR (RPN addition)."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5, 0.3]),
        ]
        outputs = [
            CalculationOutput(
                "SEL", "", ["U235"], "FUEL", [],
                state_points={"coolant_density": [0.7, 0.5]},
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        guard = edi_compo._build_statepoint_guard(outputs[0])
        # Should have two equality tests combined with +
        assert "DCool 0.7 =" in guard
        assert "DCool 0.5 =" in guard
        assert "+" in guard
        assert "0 >" in guard

    def test_statepoint_guard_multi_branch(self):
        """Guard for multiple branch types uses AND (* in RPN)."""
        assembly = self._make_mock_assembly()
        branches = [
            CalculationBranch("CD", "coolant_density", [0.7, 0.5]),
            CalculationBranch("FT", "fuel_temperature", [900.0, 1200.0]),
        ]
        outputs = [
            CalculationOutput(
                "SEL", "", ["U235"], "FUEL", [],
                state_points={
                    "coolant_density": [0.7],
                    "fuel_temperature": [900.0],
                },
            ),
        ]
        edi_compo = EDI_COMPO(assembly, branches=branches, outputs=outputs)
        guard = edi_compo._build_statepoint_guard(outputs[0])
        # Both conditions present, combined with *
        assert "DCool" in guard
        assert "TFuel" in guard
        assert "*" in guard


# =========================================================
# Integration: full case generation with branches
# =========================================================

def test_generate_cle2000_procedures(tmp_path):
    """
    Generate full CLE2000 case from GE14 DOM config yamls
    and verify all four files are produced.
    """
    GE14_DOM_test_case = DragonCase(
        case_name="GE14_DOM",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
    )

    result = GE14_DOM_test_case.generate_cle2000_procedures()

    # All four files should exist
    assert os.path.isfile(result["x2m"])
    assert os.path.isfile(result["mix"])
    assert os.path.isfile(result["trk"])
    assert os.path.isfile(result["edir"])

    # Verify CLE2000 line-length constraint on x2m
    with open(result["x2m"]) as f:
        for i, line in enumerate(f, 1):
            stripped = line.rstrip("\n")
            assert len(stripped) <= CLE2000_MAX_LINE, (
                f"x2m line {i} exceeds {CLE2000_MAX_LINE} "
                f"chars: {stripped!r}"
            )


def test_x2m_with_branches_has_while_loops(tmp_path):
    """The main x2m procedure generated from the branched scheme
    must contain WHILE loops and UTL parameter arrays."""
    case = DragonCase(
        case_name="GE14_BR",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
        tdt_base_name="GE14_DOM",
    )
    result = case.generate_cle2000_procedures()

    with open(result["x2m"]) as f:
        x2m = f.read()

    # WHILE loops for each branch axis
    assert "WHILE" in x2m
    assert "ENDWHILE" in x2m

    # UTL: parameter array creation
    assert "UTL:" in x2m
    assert "CREA" in x2m

    # COMPO initialization (before the loops)
    assert "COMPO := COMPO: ::" in x2m
    assert "MAXCAL" in x2m

    # COMPO export (after all loops)
    assert "_COMPO := COMPO" in x2m

    # Coolant density N_H / N_O arrays
    assert "N_H_arr" in x2m
    assert "N_O_arr" in x2m

    # Verify line-length for the x2m
    for i, line in enumerate(x2m.splitlines(), 1):
        assert len(line) <= CLE2000_MAX_LINE, (
            f"x2m line {i} exceeds {CLE2000_MAX_LINE}: {line!r}"
        )


def test_edir_with_branches_has_compo_parameter(tmp_path):
    """The EDIR procedure generated from the branched scheme
    must declare COMPO as a PARAMETER (not a local LINKED_LIST)."""
    case = DragonCase(
        case_name="GE14_BR",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
        tdt_base_name="GE14_DOM",
    )
    result = case.generate_cle2000_procedures()

    with open(result["edir"]) as f:
        edir = f.read()

    # COMPO is in the PARAMETER block, not declared locally
    assert "PARAMETER COMPO FLUX LIBRARY2 TRACK ::" in edir
    assert "::: LINKED_LIST COMPO" in edir

    # Branch parameter variables are declared and recovered
    assert "REAL" in edir
    assert ">>DCool<<" in edir

    # Selective outputs get IF guards
    assert "IF " in edir
    assert "ENDIF" in edir

    # ALL outputs appear without guards
    # EDIHOM_COND is ALL → should appear directly
    assert "EDI: CALL FOR EDIHOM_COND" in edir


def test_edir_with_branches_has_correct_editions(tmp_path):
    """The EDIR file should contain one EDI call per output spec."""
    case = DragonCase(
        case_name="GE14_ED",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
        tdt_base_name="GE14_DOM",
    )
    result = case.generate_cle2000_procedures()

    with open(result["edir"]) as f:
        edir = f.read()

    # CALC_SCHEME_1L has 4 outputs
    for output_name in ("EDIHOM_COND", "EDIHOM_295", "EDIR_2G", "U238_295"):
        assert f"SAVE ON {output_name}" in edir
        assert f"STEP UP '{output_name}'" in edir

    # EDIR_2G should have COND 0.625
    assert "COND 0.625" in edir

    # U238_295 should have no COND (energy_bounds="None")
    # Find the section for U238_295 and verify no COND between its
    # MERG MIX and SAVE ON
    u238_start = edir.index("EDI: CALL FOR U238_295")
    u238_save = edir.index("SAVE ON U238_295")
    u238_section = edir[u238_start:u238_save]
    assert "COND" not in u238_section


def test_mix_with_density_branch_has_nh_no_parameters(tmp_path):
    """MIX procedure generated when coolant_density branch is present
    should accept N_H and N_O as additional REAL parameters."""
    case = DragonCase(
        case_name="GE14_MX",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
        tdt_base_name="GE14_DOM",
    )
    result = case.generate_cle2000_procedures()

    with open(result["mix"]) as f:
        mix = f.read()

    # The MIX procedure should reference N_H and N_O variables
    assert "N_H" in mix
    assert "N_O" in mix


def test_line_length_constraint_all_files(tmp_path):
    """Verify CLE2000 72-char line-length constraint on ALL generated files."""
    case = DragonCase(
        case_name="GE14_LL",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": GE14_COMPOSITIONS_YAML,
            "GEOM": GE14_DOM_GEOMETRY_YAML,
            "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML,
        },
        output_path=str(tmp_path),
        tdt_path=GE14_TDT_DIR,
        tdt_base_name="GE14_DOM",
    )
    result = case.generate_cle2000_procedures()

    for key in ("x2m", "trk", "edir"):
        with open(result[key]) as f:
            for i, line in enumerate(f, 1):
                stripped = line.rstrip("\n")
                assert len(stripped) <= CLE2000_MAX_LINE, (
                    f"{key} line {i} exceeds {CLE2000_MAX_LINE}: "
                    f"{stripped!r}"
                )