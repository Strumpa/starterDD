# Tests for the programmatic DragonCalculationScheme API (plan P2):
# to_dict / to_yaml / from_dict round trips, copy, set_step, with_branch_values,
# validate, variant_tag, and identical CLE2000 procedures after a round trip.

import glob
import os
from pathlib import Path

import pytest
import yaml

from starterDD.DDModel.DragonCalculationScheme import (
    DragonCalculationScheme,
    values_match,
)
from starterDD.InterfaceToDD.case_generator import DragonCase

from conftest import (
    AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML,
    AT10_4x4_COMPOSITIONS_YAML,
    AT10_4x4_GEOMETRY_YAML,
    AT10_TDT_DIR,
    GE14_CALC_SCHEME_1L_YAML,
    GE14_COMPOSITIONS_YAML,
    GE14_DOM_GEOMETRY_YAML,
    GE14_TDT_DIR,
)

REPO = Path(__file__).resolve().parents[1]
ALL_SCHEME_YAMLS = sorted(glob.glob(str(REPO / "data" / "**" / "CALC_SCHEME*.yaml"),
                                    recursive=True))

GE14_CASE = dict(case_name="GE14_DOM", mats=GE14_COMPOSITIONS_YAML,
                 geom=GE14_DOM_GEOMETRY_YAML, tdt=GE14_TDT_DIR,
                 alias={"draglibendfb8r1SHEM295": "endfb8r1_295"})
AT10_2L_CASE = dict(case_name="AT10_4x4_mix_splitting", mats=AT10_4x4_COMPOSITIONS_YAML,
                    geom=AT10_4x4_GEOMETRY_YAML, tdt=AT10_TDT_DIR,
                    alias={"dummy": "dummy_alias"})


def _generate(out_dir, scheme, case):
    """Generate the CLE2000 procedures of *case* with *scheme* (a YAML path
    or an object) and return {file name: content}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dc = DragonCase(
        case_name=case["case_name"], call_glow=False,
        draglib_name_to_alias=case["alias"],
        config_yamls={"MATS": case["mats"], "GEOM": case["geom"], "CALC_SCHEME": scheme},
        output_path=str(out_dir), tdt_path=case["tdt"],
    )
    dc.generate_cle2000_procedures()
    return {p.name: p.read_text() for p in sorted(out_dir.iterdir()) if p.is_file()}


# ─────────────────────────────────────────────────────────────────────────────
# Round trips
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("path", ALL_SCHEME_YAMLS,
                         ids=[Path(p).parent.name + "/" + Path(p).name for p in ALL_SCHEME_YAMLS])
def test_dict_and_yaml_round_trip(path, tmp_path):
    s = DragonCalculationScheme.from_yaml(path)
    d = s.to_dict()
    assert DragonCalculationScheme.from_dict(d).to_dict() == d
    out = tmp_path / "resolved.yaml"
    s.to_yaml(out)
    again = DragonCalculationScheme.from_yaml(out)
    assert again.to_dict() == d
    assert again.variant_tag == s.variant_tag
    # The YAML is plain (safe_load) and self-describing.
    assert yaml.safe_load(out.read_text())["DRAGON_CALCULATION_SCHEME"]["steps"]


@pytest.mark.parametrize("case,yaml_path", [
    (GE14_CASE, GE14_CALC_SCHEME_1L_YAML),
    (AT10_2L_CASE, AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML),
], ids=["GE14_1L_branches", "AT10_2L"])
def test_resolved_yaml_generates_identical_procedures(case, yaml_path, tmp_path):
    """Plan P2 'done when': a YAML round trip gives identical .x2m/.c2m."""
    original = _generate(tmp_path / "orig", yaml_path, case)
    resolved = tmp_path / "resolved.yaml"
    DragonCalculationScheme.from_yaml(yaml_path).to_yaml(resolved)
    from_file = _generate(tmp_path / "file", str(resolved), case)
    from_obj = _generate(tmp_path / "obj",
                         DragonCalculationScheme.from_yaml(yaml_path).copy(), case)
    assert any(n.endswith(".x2m") for n in original)
    assert original == from_file == from_obj


def test_copy_is_independent():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    c = s.copy()
    c.get_step(c.steps[0].name).num_angles_2d = 99
    c.branches[0].values.append(0.1)
    c.outputs[1].state_points["coolant_density"].append(0.1)
    assert s.to_dict() == DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).to_dict()


# ─────────────────────────────────────────────────────────────────────────────
# set_step / with_step
# ─────────────────────────────────────────────────────────────────────────────

def test_set_step_changes_parameters_and_revalidates():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    flux = s.get_flux_steps()[-1].name
    ref = s.to_dict()
    v = s.with_step(flux, num_angles_2d=24, line_density=60.0)
    st = v.get_step(flux)
    assert (st.num_angles_2d, st.line_density) == (24, 60.0)
    assert s.to_dict() == ref                      # with_step copies
    assert v.content_hash() != s.content_hash()
    with pytest.raises(ValueError, match="MOC"):    # parser validation applies
        s.with_step(s.get_self_shielding_steps()[0].name, spatial_method="MOC")
    with pytest.raises(KeyError, match="num_angle_2d"):
        s.with_step(flux, num_angle_2d=12)          # typo
    with pytest.raises(KeyError, match="No step"):
        s.with_step("NOPE", num_angles_2d=12)


def test_set_step_merges_nested_dicts():
    s = DragonCalculationScheme.from_yaml(AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML)
    step = next(st for st in s.get_flux_steps() if st.sectorization_enabled)
    sect = s.get_step(step.name).to_dict()["sectorization"]
    s.set_step(step.name, sectorization={"windmill": not sect["windmill"]})
    new = s.get_step(step.name).to_dict()["sectorization"]
    assert new["windmill"] is (not sect["windmill"])
    assert {k: v for k, v in new.items() if k != "windmill"} == \
           {k: v for k, v in sect.items() if k != "windmill"}


# ─────────────────────────────────────────────────────────────────────────────
# with_branch_values / validate
# ─────────────────────────────────────────────────────────────────────────────

def test_with_branch_values_keeps_selective_outputs_in_sync():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    v = s.with_branch_values("coolant_density", [0.7366900001, 0.45702320],
                             labels=["00", "40"])
    br = v.get_branch("coolant_density")
    assert br.values == [0.7366900001, 0.4570232] and br.labels == ["00", "40"]
    sp = {o.name: o.state_points for o in v.outputs}
    # Kept values are rewritten to the exact new branch values.
    assert sp["EDIR_2G"]["coolant_density"] == [0.7366900001]
    assert sp["EDIHOM_295"]["coolant_density"] == [0.7366900001, 0.4570232]
    assert sp["EDIHOM_COND"] == "ALL"
    assert s.get_branch("coolant_density").values[0] == 0.73669   # original untouched


def test_with_branch_values_lost_output_raises_or_is_redirected():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    with pytest.raises(ValueError, match="'EDIR_2G'.*'U238_295'") as exc:
        s.with_branch_values("coolant_density", [0.5, 0.3])
    assert "EDIHOM_295" not in str(exc.value)     # it asked for every value
    v = s.with_branch_values("coolant_density", [0.5, 0.3],
                             selective_outputs={"EDIR_2G": [0.5], "U238_295": [0.3]})
    sp = {o.name: o.state_points for o in v.outputs}
    assert sp["EDIR_2G"]["coolant_density"] == [0.5]
    assert sp["U238_295"]["coolant_density"] == [0.3]
    assert sp["EDIHOM_295"]["coolant_density"] == [0.5, 0.3]
    a = s.with_branch_values("coolant_density", [0.5, 0.3], selective_outputs="all")
    assert all("coolant_density" not in o.state_points
               for o in a.outputs if isinstance(o.state_points, dict))
    with pytest.raises(ValueError, match="not in the branch values"):
        s.with_branch_values("coolant_density", [0.5],
                             selective_outputs={"EDIR_2G": [0.6], "U238_295": [0.5]})


def test_branch_override_keeps_every_output_in_the_procedure(tmp_path):
    """Plan P2 'done when': after a branch override every output is still
    present in the generated EDI/COMPO procedure, with a guard that fires."""
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    v = s.with_branch_values("coolant_density", [0.5, 0.3],
                             selective_outputs={"EDIR_2G": [0.5], "U238_295": [0.5]})
    files = _generate(tmp_path, v, GE14_CASE)
    edir = next(c for n, c in files.items() if n.startswith("EDIR") and n.endswith(".c2m"))
    for out in v.outputs:
        assert f"EDI: CALL FOR {out.name}" in edir
    assert "DCool 0.7366" not in edir                 # no stale guard value
    assert "DCool 0.5 - ABS" in edir


def test_validate_reports_inconsistent_outputs():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    s.outputs[2].state_points["coolant_density"] = [0.6]
    with pytest.raises(ValueError, match=r"EDIR_2G.*\[0\.6\]"):
        s.validate()
    s.branches = []                                   # no branches: no guards
    assert s.validate() is s


def test_values_match_tolerance():
    assert values_match(0.73669, 0.736690001)
    assert not values_match(0.73669, 0.7367)
    assert values_match(0.0, 0.0) and not values_match(0.0, 1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# variant_tag
# ─────────────────────────────────────────────────────────────────────────────

def test_variant_tag():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    assert s.variant_tag.startswith(s.descriptor() + "_")
    renamed = s.copy()
    renamed.name = "another name"
    assert renamed.variant_tag == s.variant_tag       # name is not content
    other = s.with_step(s.get_flux_steps()[-1].name, num_angles_2d=30)
    assert other.variant_tag != s.variant_tag
    assert s.copy().set_tag("GE14_1L_ref").variant_tag == "GE14_1L_ref"
    with pytest.raises(ValueError, match="Invalid tag"):
        s.copy().set_tag("bad tag/with spaces")
    tagged = s.copy().set_tag("v1")
    assert DragonCalculationScheme.from_dict(tagged.to_dict()).variant_tag == "v1"


# ─────────────────────────────────────────────────────────────────────────────
# P3: exact coolant statepoints
# ─────────────────────────────────────────────────────────────────────────────

GE14_ASSEMBLY = REPO / "data" / "BWRProgressionProblems" / "GE14_inputs" / "ASSEMBLY"
BENCH = {k: str(GE14_ASSEMBLY / f"material_compos_{k}.yaml") for k in ("00", "40", "80")}


def _bench_coolant(label):
    data = yaml.safe_load(open(BENCH[label]))
    entry = next(m for m in data["MIX_COMPOSITIONS"] if m["name"] == "COOLANT")
    return {"H1": float(entry["isotopic_composition"]["1001"]),
            "O16": float(entry["isotopic_composition"]["8016"])}


def _utl_array(x2m, name):
    """Values of the UTL: CREA array *name* in a generated procedure."""
    lines = x2m.splitlines()
    i = next(k for k, l in enumerate(lines) if l.split()[:1] == [name])
    vals = []
    for l in lines[i + 1:]:
        if l.strip() in ("CREA", ";") or not l.strip():
            break
        vals += [float(x) for x in l.split()]
    return vals


def test_composition_helpers():
    from starterDD.DDModel.DragonCalculationScheme import (
        composition_mass_density, normalise_composition, read_material_composition)
    assert normalise_composition({"1001": "1e-2", 8016: 5e-3, "O17": 1e-6}) == \
        {"H1": 1e-2, "O16": 5e-3, "O17": 1e-6}
    with pytest.raises(ValueError):
        normalise_composition({"8000": 1.0})          # natural element
    with pytest.raises(ValueError):
        normalise_composition({"H1": -1.0})
    comp = read_material_composition(BENCH["40"])
    assert comp == _bench_coolant("40")
    # Benchmark isotopics <-> the 0/40/80 % void densities of the GE14 YAMLs.
    for label, rho in (("00", 0.73669), ("40", 0.4570232), ("80", 0.1773564)):
        assert composition_mass_density(read_material_composition(BENCH[label])) == \
            pytest.approx(rho, rel=1e-6)


def test_with_coolant_statepoints_branch_and_outputs():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    v = s.with_coolant_statepoints(BENCH)
    br = v.get_branch("coolant_density")
    assert br.labels == ["00", "40", "80"]
    assert br.values == pytest.approx([0.73669, 0.4570232, 0.1773564], rel=1e-6)
    assert br.compositions == [_bench_coolant(k) for k in ("00", "40", "80")]
    assert br.label_of(br.values[1]) == "40"
    # EDIR_2G asked for 0.73669: kept, rewritten to the isotopics-derived value.
    sp = {o.name: o.state_points for o in v.outputs}
    assert sp["EDIR_2G"]["coolant_density"] == [br.values[0]]
    assert sp["EDIHOM_295"]["coolant_density"] == br.values   # "every statepoint"
    # Explicit PARA values are allowed.
    w = s.with_coolant_statepoints(BENCH, densities={"00": 0.73669, "40": 0.4570232,
                                                     "80": 0.1773564})
    assert w.get_branch("coolant_density").values == [0.73669, 0.4570232, 0.1773564]
    with pytest.raises(ValueError, match="same coolant density"):
        s.with_coolant_statepoints({"a": BENCH["00"], "b": BENCH["00"]})


def test_generated_mix_has_exact_benchmark_isotopics(tmp_path):
    """Plan P3 'done when': the MIX inputs for statepoint '40' are the
    isotopics of material_compos_40.yaml."""
    v = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).with_coolant_statepoints(BENCH)
    files = _generate(tmp_path, v, GE14_CASE)
    x2m = files["GE14_DOM.x2m"]
    br = v.get_branch("coolant_density")
    i40 = br.labels.index("40")
    for iso, ref in _bench_coolant("40").items():
        arr = _utl_array(x2m, f"N_{iso}_arr")
        assert len(arr) == 3
        assert arr[i40] == pytest.approx(ref, rel=1e-11)     # written with 11 digits
    assert _utl_array(x2m, "cool_dens") == pytest.approx(br.values, rel=1e-12)
    flat = " ".join(x2m.split())
    assert "GETVAL 'N_H1_arr' <<i_DCool>> >>N_H1<<" in flat
    mix = " ".join(files["MIX.c2m"].split())
    assert "H1 = H1_H2O <<N_H1>>" in mix and "O16 = O16 <<N_O16>>" in mix
    for line in x2m.splitlines():
        assert len(line) <= 72


def test_extra_coolant_isotope_is_carried_to_mix(tmp_path):
    """An isotope absent from the material YAML (here O17) becomes one more
    statepoint input; a statepoint without it gets 0."""
    comps = {"A": {"1001": 4.9e-2, "8016": 2.45e-2, "8017": 9.3e-6},
             "B": {"1001": 1.2e-2, "8016": 6.0e-3}}
    v = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML) \
        .with_coolant_statepoints(comps, selective_outputs="all")
    br = v.get_branch("coolant_density")
    assert br.composition_isotopes() == ["H1", "O16", "O17"]
    assert br.coolant_compositions()[1]["O17"] == 0.0
    files = _generate(tmp_path, v, GE14_CASE)
    assert _utl_array(files["GE14_DOM.x2m"], "N_O17_arr") == [9.3e-6, 0.0]
    mix = " ".join(files["MIX.c2m"].split())
    assert "REAL N_H1 N_O16 N_O17 ;" in mix
    assert "O17 = O17 <<N_O17>>" in mix
    assert "<<N_H1>> <<N_O16>> <<N_O17>> ;" in " ".join(files["GE14_DOM.x2m"].split())


def test_statepoint_compositions_round_trip_and_identical_procedures(tmp_path):
    v = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).with_coolant_statepoints(BENCH)
    out = tmp_path / "resolved.yaml"
    v.to_yaml(out)
    again = DragonCalculationScheme.from_yaml(out)
    assert again.to_dict() == v.to_dict()
    assert again.get_branch("coolant_density").compositions == \
        v.get_branch("coolant_density").compositions
    assert _generate(tmp_path / "a", v, GE14_CASE) == _generate(tmp_path / "b", str(out), GE14_CASE)


def test_composition_validation():
    from starterDD.DDModel.DragonCalculationScheme import CalculationBranch
    with pytest.raises(ValueError, match="only supported for coolant_density"):
        CalculationBranch("T", "fuel_temperature", [900.0], compositions=[{"H1": 1.0}])
    with pytest.raises(ValueError, match="2 compositions for 1 values"):
        CalculationBranch("D", "coolant_density", [0.7], compositions=[{"H1": 1.0}] * 2)
    with pytest.raises(ValueError, match="Cannot read isotope"):
        CalculationBranch("D", "coolant_density", [0.7], compositions=[{"water": 1.0}])
    # New values drop the old compositions (they belonged to the old values).
    v = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).with_coolant_statepoints(BENCH)
    w = v.with_branch_values("coolant_density", [0.7, 0.5], selective_outputs="all")
    assert w.get_branch("coolant_density").compositions is None
