# Tests for PostTreatment.compo_repair: diagnosis of EDIR.c2m / x2m argument
# order and exact relabelling of COMPO global parameters (bug fixed 2026-11-15,
# see the module docstring).

import gzip
from pathlib import Path

import numpy as np
import pytest

from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
from starterDD.PostTreatment import CompoReader
from starterDD.PostTreatment.compo_repair import (
    diagnose_procedures, relabel_compo_parameters, scheme_parameter_mapping,
)

from conftest import GE14_CALC_SCHEME_1L_YAML, GE14_MINI_DIR

MINI_CPO = GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz"
# What the bug did with the GE14 1L YAML (DCool declared first, loop order
# TFuel, TCool, DCool): 'DCool' received TFuel, 'TFuel' TCool, 'TCool' DCool.
BUG = {"DCool": "TFuel", "TFuel": "TCool", "TCool": "DCool"}
INVERSE = {v: k for k, v in BUG.items()}


def _params(path):
    with CompoReader(path, cache=False) as c:
        return {o: c.parameters(o) for o in c.outputs()}


def test_relabel_permutes_names_only(tmp_path):
    out = tmp_path / "permuted"
    assert relabel_compo_parameters(MINI_CPO, out, INVERSE) == 2      # one per output
    src = gzip.decompress(MINI_CPO.read_bytes())
    new = out.read_bytes()
    assert len(new) == len(src)
    changed = [a for a, b in zip(src.splitlines(), new.splitlines()) if a != b]
    assert changed and all(b"TFuel" in ln for ln in changed)          # PARKEY lines only
    before, after = _params(MINI_CPO), _params(out)
    for o in before:
        assert after[o] == {INVERSE[k]: v for k, v in before[o].items()}
    with CompoReader(MINI_CPO, cache=False) as a, CompoReader(out, cache=False) as b:
        assert a.keff("EDIHOM_COND", 1) == b.keff("EDIHOM_COND", 1)
        np.testing.assert_array_equal(a.flux("EDIR_2G", 1), b.flux("EDIR_2G", 1))
    # Repairing with the bug mapping restores the original bytes.
    back = tmp_path / "back.gz"
    relabel_compo_parameters(out, back, BUG)
    assert gzip.decompress(back.read_bytes()) == src


def test_relabel_identity_and_validation(tmp_path):
    relabel_compo_parameters(MINI_CPO, tmp_path / "same", {})
    assert (tmp_path / "same").read_bytes() == gzip.decompress(MINI_CPO.read_bytes())
    with pytest.raises(ValueError, match="permutation"):
        relabel_compo_parameters(MINI_CPO, tmp_path / "x", {"DCool": "TFuel"})
    with pytest.raises(ValueError, match="does not contain"):
        relabel_compo_parameters(MINI_CPO, tmp_path / "x", {"Burnup": "Boron", "Boron": "Burnup"})


def test_scheme_mapping_of_pre_fix_procedures():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    m = scheme_parameter_mapping(s)
    assert {k: v for k, v in m.items() if k != v} == BUG


def _write_procs(d, reads, passes, guard=False):
    d.mkdir()
    guard_txt = "IF DCool 0.5 - ABS 1.0E-05 < THEN\n" if guard else ""
    (d / "EDIR.c2m").write_text(
        "PARAMETER COMPO FLUX LIBRARY2 TRACK ::\n"
        f"REAL {' '.join(reads)} ;\n:: {' '.join(f'>>{r}<<' for r in reads)} ;\n"
        "* --------------------------------\n*    EDI: CALL FOR EDIR_2G\n"
        f"{guard_txt}EDIRATES := EDI: FLUX ;\n")
    (d / "CASE.x2m").write_text(
        "COMPO := EDIR FLUXL2 LIBEQL1 TRKFLUXL2 COMPO :: <<name_compo>> "
        + " ".join(f"<<{p}>>" for p in passes[:2]) + "\n  "
        + " ".join(f"<<{p}>>" for p in passes[2:]) + " ;\n")


def test_diagnose_procedures(tmp_path):
    ok = ["TFuel", "TCool", "DCool"]
    _write_procs(tmp_path / "ok", ok, ok)
    d = diagnose_procedures(tmp_path / "ok")
    assert not d.permuted and d.repairable and d.summary().startswith("OK")

    _write_procs(tmp_path / "bad", ["DCool", "TFuel", "TCool"], ok)
    d = diagnose_procedures(tmp_path / "bad")
    assert d.permuted and d.mapping == BUG and d.repairable

    _write_procs(tmp_path / "guard", ["DCool", "TFuel", "TCool"], ok, guard=True)
    d = diagnose_procedures(tmp_path / "guard")
    assert d.guarded_outputs == ["EDIR_2G"] and not d.repairable
    assert "cannot recover" in d.summary()


def test_diagnose_generated_procedures_are_consistent(tmp_path):
    """After the fix, freshly generated procedures diagnose as OK even with
    the GE14 YAML that declares DCool first."""
    from test_scheme_api import GE14_CASE, _generate
    _generate(tmp_path / "p", GE14_CALC_SCHEME_1L_YAML, GE14_CASE)
    assert not diagnose_procedures(tmp_path / "p").permuted


@pytest.mark.fulldata
def test_repairs_the_multi_statepoint_run(tmp_path):
    """The 2026-11-14 multi-statepoint GE14 run: relabelled COMPO has DCool
    = the three densities and identical keff."""
    import os
    root = Path(os.environ.get("STARTERDD_TEST_DATA",
                               Path(__file__).resolve().parents[1] / "temp_examples"))
    run = root / "GE14_DOM" / "results" / "PT_IC_2L_IC_MOC_3sp" / "2026-11-14T22-02-26"
    if not run.exists():
        pytest.skip("multi-statepoint run not available")
    d = diagnose_procedures(run / "procedures")
    assert d.mapping == BUG and d.repairable
    out = tmp_path / "_CPO"
    assert relabel_compo_parameters(run / "_CPO_GE14_DOM_DIAG", out, d.mapping) == 4
    with CompoReader(out, cache=False) as c:
        assert c.parameters("EDIR_2G")["TFuel"] == [900.0]
        for rho, k in [(0.73669, 0.988434136), (0.4570232, 0.971813202), (0.1773564, 0.94893986)]:
            assert c.keff("EDIHOM_COND", {"DCool": rho}) == pytest.approx(k, abs=1e-9)
