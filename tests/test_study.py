# Tests for the code-to-code study layer (plan §3.9, P5):
# configuration and YAML input, statepoint labels, the study tree written by a
# dry run, case records and the re-run policy, Serpent2 reference generation,
# the manifest, and the committed GE14 HZP/HFP study input.
import os
from pathlib import Path

import pytest
import yaml

from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.InterfaceToDD.dragon_runner import DragonRunner
from starterDD.PostTreatment.references import ReferenceRegistry
from starterDD.Workflow import CodeToCodeStudy, StudyCondition
from starterDD.Workflow.study import CASE_FILE, MANIFEST_FILE, case_statepoints

from conftest import (GE14_CALC_SCHEME_1L_YAML, GE14_COMPOSITIONS_YAML,
                      GE14_DOM_GEOMETRY_YAML, GE14_INPUTS_DIR, GE14_TDT_DIR)
from test_scheme_api import BENCH, _bench_coolant

GE14_STUDY_YAML = GE14_INPUTS_DIR / "STUDIES" / "GE14_2L_PT3x3_IC3x3_6x6_HZP_HFP.yaml"
ISOTOPES = ["U235", "U238", "Gd155", "Gd157"]


def _study(root, **kw):
    """GE14 DOM 1L (committed reference TDTs), HZP (00) and HFP (00/40/80)."""
    args = dict(
        name="GE14_test", root=root,
        assemblies={"GE14_DOM": dict(geometry=GE14_DOM_GEOMETRY_YAML,
                                     materials=GE14_COMPOSITIONS_YAML,
                                     tdt_path=GE14_TDT_DIR, scheme=GE14_CALC_SCHEME_1L_YAML)},
        conditions={"HZP": dict(coolant={"00": BENCH["00"]}, fuel_temperature=[600.0],
                                coolant_temperature=[600.0]),
                    "HFP": dict(coolant=dict(BENCH), fuel_temperature=[900.0],
                                coolant_temperature=[600.0])},
        variants={"1L": {}},
        spec=dict(preset="pin_rates", isotopes=ISOTOPES),
        draglib={"draglibendfb8r1SHEM295": "endfb8r1_295"},
    )
    args.update(kw)
    return CodeToCodeStudy(**args)


@pytest.fixture
def study(tmp_path):
    return _study(tmp_path)


# ── cases and statepoints ─────────────────────────────────────────────────


def test_cases_statepoints_and_resolved_schemes(study):
    cases = study.cases()
    assert [c.key for c in cases] == ["GE14_DOM/1L/HZP", "GE14_DOM/1L/HFP"]
    hzp, hfp = cases
    assert [sp.label for sp in hzp.statepoints] == ["HZP_00"]
    assert [sp.label for sp in hfp.statepoints] == ["HFP_00", "HFP_40", "HFP_80"]
    assert hzp.scheme.get_branch("fuel_temperature").values == [600.0]
    assert hfp.scheme.get_branch("fuel_temperature").values == [900.0]
    # exact benchmark isotopics carried by the coolant branch
    br = hfp.scheme.get_branch("coolant_density")
    for comp, label in zip(br.coolant_compositions(), ("00", "40", "80")):
        bench = _bench_coolant(label)
        assert all(abs(comp[i] / bench[i] - 1) < 1e-12 for i in bench)
    sp40 = hfp.statepoints[1]
    assert sp40.params["TFuel"] == 900.0 and sp40.params["DCool"] == br.values[1]
    # the spec sets the DRAGON outputs
    assert [o.name for o in hfp.scheme.outputs] == ["FUEL_1G", "PIN_2G"]
    assert hzp.directory == study.directory / "GE14_DOM" / "dragon" / "1L" / "HZP"


def test_unlabelled_multivalued_branch_gets_value_labels(study):
    s = StudyCondition("HFP", coolant={"00": BENCH["00"]}, fuel_temperature=[600, 900]).apply(
        study._base_scheme(study.assemblies["GE14_DOM"]))
    assert [sp.label for sp in case_statepoints(s, "HFP")] == ["HFP_TF600_00", "HFP_TF900_00"]


def test_variant_steps_change_the_scheme_and_the_fingerprint(tmp_path):
    st = _study(tmp_path, variants={"a": {}, "b": {"steps": {"FLUX": {"num_angles_2d": 6}}}})
    a, b = st.case("GE14_DOM", "a", "HFP"), st.case("GE14_DOM", "b", "HFP")
    assert b.scheme.get_step("FLUX").num_angles_2d == 6
    assert a.fingerprint != b.fingerprint
    assert a.scheme.variant_tag != b.scheme.variant_tag


# ── DRAGON: dry run, records, re-run policy ───────────────────────────────


def test_dry_run_creates_the_full_tree(study):
    res = study.run_dragon(dry_run=True)
    assert res == {"GE14_DOM/1L/HZP": "dry_run", "GE14_DOM/1L/HFP": "dry_run"}
    for case in study.cases():
        d = case.directory
        rec = yaml.safe_load((d / CASE_FILE).read_text())
        assert rec["status"] == "dry_run" and rec["fingerprint"] == case.fingerprint
        assert [s["label"] for s in rec["statepoints"]] == [sp.label for sp in case.statepoints]
        assert (d / "scheme.yaml").exists()
        assert (d / "procedures" / "GE14_DOM.x2m").exists()
        run = d / rec["runs"][-1]["run_directory"]
        assert run.parent == d                      # runs directly in the case folder
        assert (run / "run_manifest.yaml").exists() and (run / "_staging_contents.txt").exists()
        assert (run / "inputs" / "CALC_SCHEME_resolved.yaml").exists()
    man = yaml.safe_load((study.directory / MANIFEST_FILE).read_text())
    assert man["spec"] == ["FUEL_1G", "PIN_2G"]
    sps = [sp for c in man["cases"] for sp in c["statepoints"]]
    assert [sp["label"] for sp in sps] == ["HZP_00", "HFP_00", "HFP_40", "HFP_80"]
    assert sps[2]["reference"] == "references/GE14_DOM/endfb8r1/HFP_40"
    assert sps[2]["compo_parameters"]["TFuel"] == 900.0
    assert (study.directory / "study_config.yaml").exists()


def test_completed_cases_are_skipped_and_changed_inputs_are_stale(tmp_path):
    st = _study(tmp_path)
    case = st.case("GE14_DOM", "1L", "HZP")
    case.write_record({"key": case.key, "fingerprint": case.fingerprint, "status": "completed",
                       "current_run": "somewhere"})
    assert case.status == "completed"
    assert st.run_dragon(dry_run=True, conditions="HZP") == {case.key: "skipped"}
    # same study folder, another fuel temperature: the record is stale
    st2 = _study(tmp_path, conditions={"HZP": dict(coolant={"00": BENCH["00"]},
                                                   fuel_temperature=[550.0])})
    assert st2.case("GE14_DOM", "1L", "HZP").status == "stale"


def test_real_run_with_a_stand_in_executable(tmp_path):
    """Non-dry path: a fake DRAGON that ends normally and writes a COMPO."""
    exe = tmp_path / "fake_dragon"
    exe.write_text("#!/bin/sh\ncat > /dev/null\necho 'kinf= 1.012345'\n"
                   "touch _CPO_GE14_DOM\necho 'normal end'\n")
    exe.chmod(0o755)
    lib = tmp_path / "draglibendfb8r1SHEM295"
    lib.write_text("dummy")
    st = _study(tmp_path, dragon_executable=exe, draglib_paths={lib.name: lib})
    assert st.run_dragon(conditions="HZP") == {"GE14_DOM/1L/HZP": "completed"}
    case = st.case("GE14_DOM", "1L", "HZP")
    rec = case.record
    assert rec["status"] == "completed" and rec["runs"][-1]["kinf"] == pytest.approx(1.012345)
    assert case.run_directory == case.directory / rec["current_run"]
    assert case.compo_path == case.run_directory / "_CPO_GE14_DOM"
    assert (case.directory / "latest").resolve() == case.run_directory.resolve()
    assert st.run_dragon(conditions="HZP") == {case.key: "skipped"}
    assert st.run_dragon(conditions="HZP", force=True) == {case.key: "completed"}
    man = yaml.safe_load(st.manifest_path.read_text())
    assert man["cases"][0]["compo"].endswith("_CPO_GE14_DOM")
    st.write_serpent_references(conditions="HZP")
    rows = [r for r in st.post_process() if r["condition"] == "HZP"]
    assert [(r["statepoint"], r["status"]) for r in rows] == [("HZP_00", "reference not run yet")]


def test_failed_case_is_recorded_and_the_next_one_runs(tmp_path):
    st = _study(tmp_path, assemblies={"GE14_DOM": dict(
        geometry=GE14_DOM_GEOMETRY_YAML, materials=GE14_COMPOSITIONS_YAML,
        tdt_path=tmp_path / "no_tdt", scheme=GE14_CALC_SCHEME_1L_YAML)})
    res = st.run_dragon(dry_run=True)
    assert set(res.values()) == {"failed"}
    rec = st.case("GE14_DOM", "1L", "HFP").record
    assert "TDT file not found" in rec["runs"][-1]["error"]
    with pytest.raises(FileNotFoundError):
        st.run_dragon(dry_run=True, stop_on_error=True)


def test_runner_run_subdirs(tmp_path):
    dc = DragonCase(case_name="GE14_DOM", call_glow=False,
                    draglib_name_to_alias={"draglibendfb8r1SHEM295": "endfb8r1_295"},
                    config_yamls={"MATS": GE14_COMPOSITIONS_YAML, "GEOM": GE14_DOM_GEOMETRY_YAML,
                                  "CALC_SCHEME": GE14_CALC_SCHEME_1L_YAML},
                    output_path=str(tmp_path / "procs"), tdt_path=GE14_TDT_DIR)
    res = dc.run(results_root=str(tmp_path / "flat"), dry_run=True, run_subdirs=[])
    assert Path(res.run_directory).parent == tmp_path / "flat"
    res = DragonRunner(dc, results_root=str(tmp_path / "default")).run(dry_run=True)
    assert Path(res.run_directory).parent.parent.name == "GE14_DOM"     # unchanged default


# ── Serpent references ────────────────────────────────────────────────────


def test_serpent_references_one_per_statepoint(study):
    refs = study.write_serpent_references()
    assert [r.meta["sp_label"] for r in refs] == ["HZP_00", "HFP_00", "HFP_40", "HFP_80"]
    ref40 = ReferenceRegistry(study.references_root).find("GE14_DOM", "endfb8r1", "HFP_40")
    assert ref40.meta["statepoint"]["fuel_temperature"] == 900.0
    bench = _bench_coolant("40")
    assert all(abs(ref40.meta["coolant_composition"][i] / bench[i] - 1) < 1e-10 for i in bench)
    hzp = ReferenceRegistry(study.references_root).find("GE14_DOM", "endfb8r1", "HZP_00")
    assert hzp.meta["statepoint"]["fuel_temperature"] == 600.0
    assert study.check_references() == {}           # every case covered
    assert study.write_serpent_references() == []   # existing references are kept
    rows = {r["case"]: r for r in study.status()}
    assert rows["GE14_DOM/1L/HFP"]["references"] == "3/3"
    assert rows["GE14_DOM/1L/HFP"]["serpent_results"] == "0/3"


def test_reference_spec_finer_than_the_dragon_spec(tmp_path):
    st = _study(tmp_path, conditions={"HZP": dict(coolant={"00": BENCH["00"]})},
                reference_spec=dict(preset="pin_rates_and_spectra", isotopes=ISOTOPES))
    assert "PIN_295G" in st.reference_spec().names()
    st.write_serpent_references()
    assert st.check_references() == {}


def test_post_process_without_results(study):
    study.write_serpent_references(conditions="HZP")
    rows = study.post_process()
    assert {r["status"] for r in rows} == {"no COMPO (planned)"}
    assert (study.directory / "post" / "summary.csv").exists()


# ── imported DRAGON run, normalisation, figures ───────────────────────────


def _plain_compo(tmp_path):
    """The committed mini COMPO (GE14 DOM 1L, HFP 00) as a plain run folder."""
    import gzip
    import shutil
    from conftest import GE14_MINI_DIR
    run = tmp_path / "old_run"
    run.mkdir()
    with gzip.open(GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz", "rb") as f, \
            open(run / "_CPO_GE14_DOM", "wb") as g:
        shutil.copyfileobj(f, g)
    return run


def test_import_dragon_run_post_process_and_figures(tmp_path):
    """An existing COMPO is imported without running DRAGON; the reference is
    a synthetic Serpent result (DRAGON x 3.7, U238 x 1.02). fission_source
    recovers 3.7, U235 rates agree, U238 shows its bias, and the figures are
    made from study.h5 alone."""
    import h5py
    import numpy as np
    from starterDD.PostTreatment import SerpentRun, physics_data as phys
    from starterDD.PostTreatment.compo_reader import CompoReader
    from test_macro_rates import _write_det0, synthetic_pin_blocks
    hfp00 = {"HFP": dict(coolant={"00": BENCH["00"]}, fuel_temperature=[900.0],
                         coolant_temperature=[600.0])}
    st = _study(tmp_path, conditions=hfp00)
    run = _plain_compo(tmp_path)
    case = st.import_dragon_run("GE14_DOM", "1L", "HFP", run)
    assert case.status == "completed" and case.record["origin"] == "imported"
    assert case.compo_path.is_symlink() and case.compo_path.resolve() == run / "_CPO_GE14_DOM"
    assert st.run_dragon(dry_run=True) == {case.key: "skipped"}      # never re-run
    assert st.import_dragon_run("GE14_DOM", "1L", "HFP", run) is case  # kept

    # synthetic Serpent result for the generated reference
    (ref,) = st.write_serpent_references()
    lat = st.lattice("GE14_DOM")
    with CompoReader(case.compo_path, cache=False) as c:
        blocks = synthetic_pin_blocks(c, "EDIR_2G", SerpentRun(ref.input_path), "PIN_2G", lat,
                                      3.7, {"U238": 1.02})
    _write_det0(str(ref.input_path) + "_det0.m", blocks, phys.G2_EV)
    Path(str(ref.input_path) + "_res.m").write_text(
        "IMP_KEFF                  (idx, [1:   2]) = [  9.90261E-01 7.2E-06 ];\n")

    rows = st.post_process(isotope_sets=[["U235"], ["U238"]])
    by = {(r["item"], r.get("isotopes"), r.get("quantity"), r.get("group")): r for r in rows}
    assert all(r["status"] == "ok" for r in rows), [r for r in rows if r["status"] != "ok"]
    assert by[("keff", None, None, None)]["delta"] == pytest.approx(-117.1511, abs=1e-3)
    norm = by[("normalisation", None, "fission_source", None)]
    assert norm["factor"] == pytest.approx(3.7, rel=1e-9) and norm["normalisation"] == "fission_source"
    for g in (1, 2):
        u5 = by[("pin_rates", "U235", "absorption", g)]
        u8 = by[("pin_rates", "U238", "absorption", g)]
        assert u5["rms"] == pytest.approx(0.0, abs=1e-6)
        assert u8["mean"] == pytest.approx(100 * (1 / 1.02 - 1), abs=1e-6)
        assert by[("macro_rates", "all", "nufission", g)]["rms"] == pytest.approx(0.0, abs=1e-6)
    h5p = st.directory / "post" / "study.h5"
    with h5py.File(h5p, "r") as h5:
        g = h5[f"{case.key}/HFP_00"]
        assert g["normalisation"].attrs["method"] == "fission_source"
        assert g["pin_rates/PIN_2G/U238/absorption"].attrs["factor"] == pytest.approx(3.7)
        assert h5[case.key].attrs["origin"] == "imported"
        assert list(h5["GE14_DOM/lattice/ordered_pins"]) == lat.ordered_pins
    # subset normalisation: U235 alone is still exact, but the factor changes per set
    rows_ps = st.post_process(isotope_sets=[["U235", "U238"]], normalisation="power_subset")
    r = next(r for r in rows_ps if r["item"] == "pin_rates" and r["quantity"] == "absorption")
    assert r["factor"] != pytest.approx(3.7, rel=1e-4)

    st.post_process(isotope_sets=[["U235"], ["U238"]])
    figs = st.figures()
    names = sorted(f.name for f in figs)
    assert "keff_GE14_DOM.png" in names
    assert "maps_pin_rates_PIN_2G_U238_absorption.png" in names
    assert "maps_macro_rates_PIN_2G_nufission.png" in names
    assert all(f.exists() and f.stat().st_size > 1000 for f in figs)


def test_import_dragon_run_checks_the_statepoints(tmp_path):
    """The single-statepoint COMPO cannot be the HFP 00/40/80 case."""
    st = _study(tmp_path)
    run = _plain_compo(tmp_path)
    with pytest.raises(ValueError, match="HFP_40"):
        st.import_dragon_run("GE14_DOM", "1L", "HFP", run)
    case = st.case("GE14_DOM", "1L", "HFP")
    assert case.status == "planned" and not list(case.directory.glob("imported_*"))
    st2 = _study(tmp_path / "b", conditions={"HFP": dict(coolant={"00": BENCH["00"]},
                                                         fuel_temperature=[600.0])})
    with pytest.raises(ValueError, match="match"):
        st2.import_dragon_run("GE14_DOM", "1L", "HFP", run)   # TFuel 600 vs 900


def test_spec_config_renamed_keys(tmp_path):
    from starterDD.Workflow.study import spec_from_config
    with pytest.raises(KeyError, match="fuel_isotopes_inventory"):
        spec_from_config(dict(isotopes=ISOTOPES, keff_isotopes=ISOTOPES))
    with pytest.raises(KeyError, match="serpent_fine_pins"):
        spec_from_config(dict(isotopes=ISOTOPES, pins_fine=[1]))
    spec = spec_from_config(dict(preset="pin_rates_and_spectra", isotopes=ISOTOPES,
                                 serpent_fine_pins="diagonal"))
    st = _study(tmp_path, spec=dict(preset="pin_rates_and_spectra", isotopes=ISOTOPES,
                                    serpent_fine_pins="diagonal"),
                conditions={"HZP": dict(coolant={"00": BENCH["00"]})})
    case = st.cases()[0]
    resolved = st.case_spec(case)
    assert spec.pins_unresolved and not resolved.pins_unresolved
    assert resolved["PIN_295G"].pins == sorted(st.lattice("GE14_DOM").diagonal_pins("axis"))
    st.write_serpent_references()
    assert st.check_references() == {}


# ── YAML input ────────────────────────────────────────────────────────────


def _write_yaml(tmp_path, text):
    p = tmp_path / "study.yaml"
    p.write_text(text)
    return p


STUDY_TEXT = f"""
variables:
  ROOT: default_root
  TDT: {GE14_TDT_DIR}
study:
  name: mini
  root: ${{ROOT}}
dragon: {{draglib: {{draglibendfb8r1SHEM295: endfb8r1_295}}, num_threads: 4}}
spec: {{preset: pin_rates, isotopes: [U235, U238]}}
assemblies:
  GE14_DOM:
    geometry: {GE14_DOM_GEOMETRY_YAML}
    materials: {GE14_COMPOSITIONS_YAML}
    scheme: {GE14_CALC_SCHEME_1L_YAML}
    tdt_path: ${{TDT}}
conditions:
  HFP:
    coolant: {{"00": {BENCH["00"]}, "80": {BENCH["80"]}}}
    fuel_temperature: 900
"""


def test_from_yaml_variables_and_relative_paths(tmp_path, monkeypatch):
    p = _write_yaml(tmp_path, STUDY_TEXT)
    st = CodeToCodeStudy.from_yaml(p)
    assert st.root == (tmp_path / "default_root").resolve()   # relative to the file
    assert st.num_threads == 4 and st.variants.keys() == {"base"}
    monkeypatch.setenv("ROOT", str(tmp_path / "from_env"))
    assert CodeToCodeStudy.from_yaml(p).root == tmp_path / "from_env"
    st = CodeToCodeStudy.from_yaml(p, variables={"ROOT": str(tmp_path / "explicit")})
    assert st.root == tmp_path / "explicit"
    assert [sp.label for sp in st.case("GE14_DOM", "base", "HFP").statepoints] == \
        ["HFP_00", "HFP_80"]


def test_from_yaml_errors(tmp_path):
    with pytest.raises(KeyError, match="Undefined variable"):
        CodeToCodeStudy.from_yaml(_write_yaml(tmp_path, STUDY_TEXT.replace("${TDT}", "${NOPE}")))
    with pytest.raises(KeyError, match="unknown section"):
        CodeToCodeStudy.from_yaml(_write_yaml(tmp_path, STUDY_TEXT + "typo_section: 1\n"))
    with pytest.raises(ValueError, match="folder names"):
        _study(tmp_path, name="bad name")


def test_ge14_study_input(tmp_path):
    """The committed GE14 2L PT 3x3 / IC 3x3_6x6 study: 3 lattices x HZP/HFP,
    12 statepoints, spec outputs on 26 groups, references on SHEM295."""
    st = CodeToCodeStudy.from_yaml(GE14_STUDY_YAML, variables={"STUDY_ROOT": str(tmp_path)})
    assert list(st.assemblies) == ["GE14_DOM", "GE14_DOM-C", "GE14_VAN"]
    cases = st.cases()
    assert len(cases) == 6
    labels = {c.key: [sp.label for sp in c.statepoints] for c in cases}
    assert labels["GE14_DOM-C/2L_PT3x3_IC3x3_6x6/HZP"] == ["HZP_00"]
    assert labels["GE14_VAN/2L_PT3x3_IC3x3_6x6/HFP"] == ["HFP_00", "HFP_40", "HFP_80"]
    for c in cases:
        tf = c.scheme.get_branch("fuel_temperature").values
        assert tf == ([600.0] if c.condition.name == "HZP" else [900.0])
        assert c.scheme.get_branch("coolant_temperature").values == [600.0]
        assert c.scheme.is_two_level_scheme()
        assert [o.name for o in c.scheme.outputs] == [
            "FUEL_1G", "PIN_2G", "PIN_26G", "FUEL_26G", "ALL_26G"]
        assert c.assembly.geometry.exists() and c.assembly.materials.exists()
    assert st.reference_spec().names() == ["FUEL_1G", "PIN_2G", "PIN_295G", "FUEL_295G",
                                           "ALL_295G"]
    assert st.serpent_settings["neutrons_per_cycle"] == 2000000
    assert set(st.serpent_imports) == {"GE14_DOM", "GE14_DOM-C", "GE14_VAN"}
    assert st.case("GE14_DOM-C", "2L_PT3x3_IC3x3_6x6", "HFP").scheme.get_step(
        "FLUXL2").tdt_file_id == "FLUXL2_wr_8_cross_submesh_extend_moder308"
