"""Serpent2 runner (plan P6): executable resolution, input checks, runs with
a stand-in ``sss2``, reference records and the study ``run_serpent``."""

import os
import stat

import pytest
import yaml

from starterDD.InterfaceToDD.serpent_runner import (Serpent2Runner, SerpentExecutableError,
                                                    check_input, data_library_paths,
                                                    resolve_serpent_executable)
from starterDD.PostTreatment.references import SerpentReference

FAKE_SSS2 = """#!/bin/sh
# stand-in for sss2: sss2 -omp N <input>
echo "args: $@"
echo "omp_env: $OMP_NUM_THREADS"
[ -n "$FAKE_SSS2_FAIL" ] && { echo "fatal error"; exit 3; }
inp="$3"
if [ -z "$FAKE_SSS2_NO_DET" ] && grep -q '^det' "$inp"; then
  echo "DET1 = [ 1 1 1 1 1 1 1 1 1.0 0.01 ];" > "${inp}_det0.m"
fi
echo "IMP_KEFF (idx, [1:   2]) = [  1.01234E+00  0.00050 ];" >> "${inp}_res.m"
echo "SIMULATED_HISTORIES (idx, 1) = 12345 ;" >> "${inp}_res.m"
"""


@pytest.fixture
def sss2(tmp_path, monkeypatch):
    exe = tmp_path / "bin" / "sss2"
    exe.parent.mkdir()
    exe.write_text(FAKE_SSS2)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    for v in ("SERPENT_EXEC", "sss2", "FAKE_SSS2_FAIL", "FAKE_SSS2_NO_DET", "SERPENT_ACELIB"):
        monkeypatch.delenv(v, raising=False)
    return exe


def _deck(folder, name="case.serp", detectors=True, acelib=None):
    folder.mkdir(parents=True, exist_ok=True)
    ace = acelib or (folder / "lib.xsdata")
    if acelib is None:
        ace.write_text("")
    lines = [f'set acelib "{ace}"  % nuclear data', "set pop 100 10 5"]
    if detectors:
        lines.append("det D1 de 2g")
    inp = folder / name
    inp.write_text("\n".join(lines) + "\n")
    return inp


def test_executable_resolution(sss2, monkeypatch, tmp_path):
    assert resolve_serpent_executable(sss2) == str(sss2.resolve())
    monkeypatch.setenv("SERPENT_EXEC", str(sss2))
    assert resolve_serpent_executable() == str(sss2.resolve())
    monkeypatch.delenv("SERPENT_EXEC")
    monkeypatch.setenv("PATH", str(sss2.parent) + os.pathsep + os.environ.get("PATH", ""))
    assert resolve_serpent_executable() == str(sss2)
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    with pytest.raises(SerpentExecutableError, match="SERPENT_EXEC"):
        resolve_serpent_executable(tmp_path / "missing_sss2")


def test_input_checks(tmp_path):
    inp = _deck(tmp_path / "a")
    assert data_library_paths(inp) == [("acelib", str(tmp_path / "a" / "lib.xsdata"))]
    assert check_input(inp) == []
    bad = _deck(tmp_path / "b", acelib=tmp_path / "nowhere.xsdata")
    assert check_input(bad) == [f"set acelib: {tmp_path / 'nowhere.xsdata'} not found"]
    assert check_input(tmp_path / "none.serp") == [f"input {tmp_path / 'none.serp'} not found"]


def test_run_skip_force_and_failures(sss2, tmp_path, monkeypatch):
    inp = _deck(tmp_path / "ref")
    runner = Serpent2Runner(sss2, num_threads=4)
    dry = runner.run(inp, dry_run=True)
    assert dry.status == "dry_run" and "-omp 4 case.serp" in dry.message
    assert not (tmp_path / "ref" / "case.serp_res.m").exists()

    res = runner.run(inp)
    assert res.success and res.return_code == 0
    assert res.keff == pytest.approx((1.01234, 0.0005))
    log = (tmp_path / "ref" / "case.serp.log").read_text()
    assert "args: -omp 4 case.serp" in log and "omp_env: 4" in log      # run in the folder
    assert runner.run(inp).status == "skipped"
    # force: old outputs removed first (the stand-in appends to _res.m)
    assert runner.run(inp, force=True).success
    assert (tmp_path / "ref" / "case.serp_res.m").read_text().count("IMP_KEFF") == 1

    monkeypatch.setenv("FAKE_SSS2_FAIL", "1")
    failed = runner.run(inp, force=True)
    assert failed.status == "failed" and "return code 3" in failed.message
    monkeypatch.delenv("FAKE_SSS2_FAIL")
    monkeypatch.setenv("FAKE_SSS2_NO_DET", "1")
    nodet = runner.run(inp, force=True)
    assert nodet.status == "failed" and "case.serp_det0.m" in nodet.message


def test_bad_data_path_and_symlinked_outputs_are_refused(sss2, tmp_path):
    runner = Serpent2Runner(sss2)
    bad = _deck(tmp_path / "bad", acelib=tmp_path / "nowhere.xsdata")
    assert runner.run(bad, dry_run=True).status == "not_ready"
    r = runner.run(bad)
    assert r.status == "failed" and "not found" in r.message
    assert not (tmp_path / "bad" / "case.serp.log").exists()          # never started
    inp = _deck(tmp_path / "linked")
    original = tmp_path / "original_res.m"
    original.write_text("ORIGINAL")
    (tmp_path / "linked" / "case.serp_res.m").symlink_to(original)
    r = runner.run(inp, force=True)
    assert r.status == "failed" and "symlink" in r.message
    assert original.read_text() == "ORIGINAL"


def test_reference_run_serpent_records_and_refuses_imported(sss2, tmp_path):
    folder = tmp_path / "refs" / "A" / "endfb8r1" / "HFP_00"
    inp = _deck(folder, "A_HFP_00.serp")
    (folder / "reference.yaml").write_text(yaml.safe_dump(
        {"assembly": "A", "library": "endfb8r1", "sp_label": "HFP_00",
         "origin": "generated", "input": inp.name, "settings": {}}))
    ref = SerpentReference.open(folder)
    assert not ref.has_results
    res = ref.run_serpent(sss2, num_threads=2)
    assert res.success and ref.has_results
    meta = yaml.safe_load((folder / "reference.yaml").read_text())
    assert meta["runs"][-1]["status"] == "completed" and meta["runs"][-1]["num_threads"] == 2
    assert meta["runs"][-1]["keff"] == pytest.approx([1.01234, 0.0005])
    assert meta["settings"]["simulated_histories"] == 12345
    assert ref.last_run["log"] == "A_HFP_00.serp.log"
    meta["origin"] = "imported"
    (folder / "reference.yaml").write_text(yaml.safe_dump(meta))
    r = SerpentReference.open(folder).run_serpent(sss2, force=True)
    assert r.status == "skipped" and "imported" in r.message


def test_study_run_serpent(sss2, tmp_path, capsys):
    from test_study import _study
    ace = tmp_path / "endfb8r1.xsdata"
    ace.write_text("")
    st = _study(tmp_path / "s", serpent_executable=sss2, serpent_num_threads=3,
                serpent_settings={"acelib": str(ace)})
    assert st.run_serpent(conditions="HZP") == {"GE14_DOM/HZP_00": "no reference"}
    st.write_serpent_references()
    assert st.run_serpent(dry_run=True, statepoints=["HFP_40"]) == {"GE14_DOM/HFP_40": "dry_run"}
    assert st.run_serpent(conditions="HFP", statepoints=["HFP_40"]) == {"GE14_DOM/HFP_40": "completed"}
    assert "keff = 1.012340" in capsys.readouterr().out
    rows = {r["case"]: r for r in st.status()}
    assert rows["GE14_DOM/1L/HFP"]["serpent_results"] == "1/3"
    out = st.run_serpent()
    assert out["GE14_DOM/HFP_40"] == "skipped" and out["GE14_DOM/HZP_00"] == "completed"
    assert st.to_dict()["serpent"]["num_threads"] == 3
    # post-processing now finds results (no DRAGON COMPO yet, so no comparison)
    assert {r["status"] for r in st.post_process()} == {"no COMPO (planned)"}


def test_cli_run_serpent_and_executable_variable(sss2, tmp_path, monkeypatch, capsys):
    from starterDD.Workflow.study import main
    from test_study import STUDY_TEXT, _write_yaml
    monkeypatch.setenv("SERPENT_EXEC", str(sss2))
    study_yaml = _write_yaml(tmp_path, STUDY_TEXT)
    from starterDD.Workflow import CodeToCodeStudy
    st = CodeToCodeStudy.from_yaml(study_yaml, site=False)
    assert st.serpent_executable == str(sss2)
    assert st.variables["SERPENT_EXEC"]["source"] == "environment"
    assert main([str(study_yaml), "run-serpent", "--site", "none", "--dry-run"]) == 1
    assert "no reference" in capsys.readouterr().out
