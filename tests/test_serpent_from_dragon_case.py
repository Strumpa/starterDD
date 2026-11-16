# Plan P4 'done when' (second half): the Serpent2 reference generated for a
# statepoint of a DRAGON case has exactly the coolant isotopics that the
# DRAGON procedure uses at that statepoint, and the temperatures of its MIX.

import pytest

from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.InterfaceToDD.Serpent2_exports import Serpent2Model
from starterDD.PostTreatment import DetectorManifest, parse_serpent_input
from starterDD.PostTreatment.spec import ComparisonSpec

from conftest import GE14_CALC_SCHEME_1L_YAML
from test_scheme_api import BENCH, GE14_CASE, _bench_coolant, _utl_array

ISOS = ["U235", "U238", "Gd155", "Gd157"]


@pytest.fixture(scope="module")
def ge14_case(tmp_path_factory):
    """GE14 1L case with the benchmark 00/40/80 statepoints, procedures generated."""
    out = tmp_path_factory.mktemp("case")
    scheme = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).with_coolant_statepoints(BENCH)
    dc = DragonCase(case_name=GE14_CASE["case_name"], call_glow=False,
                    draglib_name_to_alias=GE14_CASE["alias"],
                    config_yamls={"MATS": GE14_CASE["mats"], "GEOM": GE14_CASE["geom"],
                                  "CALC_SCHEME": scheme},
                    output_path=str(out), tdt_path=GE14_CASE["tdt"])
    dc.generate_cle2000_procedures()
    x2m = (out / "GE14_DOM.x2m").read_text()
    return dc, scheme, x2m


@pytest.mark.parametrize("label", ["00", "40", "80"])
def test_coolant_isotopics_identical_in_both_codes(ge14_case, label, tmp_path):
    dc, scheme, x2m = ge14_case
    br = scheme.get_branch("coolant_density")
    i = br.labels.index(label)
    model = Serpent2Model.from_dragon_case(dc, label)
    assert model.statepoint["coolant_density"] == br.values[i]
    serp = tmp_path / f"GE14_{label}.serp"
    model.write(str(serp))
    inp = parse_serpent_input(serp)
    cool = inp.materials["coolant"]
    written = cool.composition()
    assert set(written) == set(br.composition_isotopes())
    for iso, n_serp in written.items():
        n_dragon = _utl_array(x2m, f"N_{iso}_arr")[i]          # what DRAGON reads
        assert n_serp == pytest.approx(n_dragon, rel=1e-10)     # same 11 digits
        assert n_serp == pytest.approx(_bench_coolant(label)[iso], rel=1e-10)
    assert cool.temperature == 600.0
    assert inp.materials["moderator"].temperature == 600.0
    assert inp.materials["gap"].temperature == 900.0           # MIX 3 <<TFUEL>>
    fuel = [m for n, m in inp.materials.items() if "_pin_" in n]
    assert fuel and all(m.temperature == 900.0 for m in fuel)


def test_reference_with_spec_detectors_and_manifest(ge14_case, tmp_path):
    dc, _, _ = ge14_case
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    model = Serpent2Model.from_dragon_case(dc, "40", spec=spec)
    serp = tmp_path / "GE14_40.serp"
    model.write(str(serp))
    man = DetectorManifest.read(tmp_path / "GE14_40.serp.detectors.json")
    assert man.source["inferred"] is False and spec.coverage(man) == {}
    assert man.summary() == DetectorManifest.infer(serp).summary()
    assert "coolant_density=0.457" in model.settings.title


def _legacy_manifest(which):
    import os
    from pathlib import Path
    from conftest import GE14_MINI_DIR
    if which == "mini":
        return DetectorManifest.read(GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.detectors.json")
    root = Path(os.environ.get("STARTERDD_TEST_DATA",
                               Path(__file__).resolve().parents[1] / "temp_examples"))
    serp = root / "GE14_DOM" / "Serpent" / "GE14_DOM_00_endfb8r1.serp"
    if not serp.exists():
        pytest.skip("full GE14 reference not available")
    return DetectorManifest.infer(serp)


@pytest.mark.parametrize("which", ["mini", pytest.param("full", marks=pytest.mark.fulldata)])
def test_generated_reference_covers_the_legacy_detector_set(ge14_case, tmp_path, which):
    """Plan P4 'done when': a generated reference covers the detector set of
    the existing GE14 references, detector by detector (same observable,
    kind and pin; every (MT, nuclide) response except the MT 101 summary;
    the DRAGON SHEM295 grid instead of the legacy one)."""
    dc, _, _ = ge14_case
    legacy = _legacy_manifest(which)
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    model = Serpent2Model.from_dragon_case(dc, "00", spec=spec)
    model.write(str(tmp_path / "GE14_00.serp"))
    new = DetectorManifest.read(tmp_path / "GE14_00.serp.detectors.json")
    for e in legacy:
        g = new.find(e.observable, e.kind, e.pin_idx)
        assert g.materials == e.materials
        want = {(r.mt, r.nuclide) for r in e.responses if r.mt != 101}
        assert want <= {(r.mt, r.nuclide) for r in g.responses}, e.detector
        assert len(g.bounds_eV) == len(e.bounds_eV) == 296


def test_generated_reference_folder(ge14_case, tmp_path):
    from starterDD.PostTreatment.references import (
        ReferenceCoverageError, ReferenceRegistry, SerpentReference)
    dc, scheme, _ = ge14_case
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    model = Serpent2Model.from_dragon_case(dc, "80", spec=spec)
    ref = SerpentReference.create(tmp_path, "GE14_DOM", "endfb8r1", "80", model, spec=spec)
    assert ref.path == tmp_path / "GE14_DOM" / "endfb8r1" / "80"
    assert ref.meta["origin"] == "generated" and ref.meta["spec"] == spec.names()
    assert ref.meta["coolant_composition"] == _bench_coolant("80")
    assert ref.meta["statepoint"]["coolant_density"] == scheme.get_branch("coolant_density").values[2]
    assert not ref.has_results
    with pytest.raises(FileNotFoundError):
        ref.run()
    ref.require(spec, n_pins=51)
    with pytest.raises(ReferenceCoverageError, match="PIN_26G"):
        ref.require(ComparisonSpec.preset("pin_rates_and_spectra", ["U235", "Pu239"],
                                          fine_mesh="26g"), n_pins=51)
    with pytest.raises(FileExistsError):
        SerpentReference.create(tmp_path, "GE14_DOM", "endfb8r1", "80", model)
    assert ReferenceRegistry(tmp_path).find("GE14_DOM", "endfb8r1", "80").meta == ref.meta


def test_import_existing_legacy_reference(tmp_path):
    """Decision 8: an existing run (committed GE14 subset, gzipped, legacy
    names and grid) is registered without re-running and post-processed."""
    from conftest import GE14_MINI_DIR
    from starterDD.PostTreatment.references import (
        ReferenceCoverageError, ReferenceRegistry, SerpentReference)
    serp = GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.gz"
    ref = SerpentReference.import_existing(serp, tmp_path, "GE14_DOM", "endfb8r1", "00")
    assert ref.meta["origin"] == "imported" and ref.has_results
    assert (ref.path / serp.name).is_symlink()
    assert ref.meta["settings"]["neutrons_per_cycle"] == 2_000_000
    k, _ = ref.run().keff()
    assert k == pytest.approx(0.990261)
    pins = ref.manifest.pins("PIN_295G", "rr")                   # the 4 committed pins
    # (legacy runs have no 'void' responses: no whole-assembly totals)
    spec26 = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, fine_mesh="26g",
                                   global_totals=False)
    ref.require(spec26, pins=pins)
    with pytest.raises(ReferenceCoverageError, match="void"):
        ref.require(ComparisonSpec.preset("pin_rates_and_spectra", ISOS, fine_mesh="26g"),
                    pins=pins)
    # A generated reference would score PIN_26G; this one provides it (and
    # PIN_2G) through its legacy per-pin 295g detectors.
    assert spec26.serpent_source("PIN_2G") == "PIN_26G"
    assert ref.source(spec26, "PIN_2G", pins=pins) == "PIN_295G"
    assert ref.source(spec26, "ALL_26G", "flux") == "ALL_295G"
    with pytest.raises(ReferenceCoverageError, match="legacy 295g grid"):
        ref.require(ComparisonSpec.preset("pin_rates_and_spectra", ISOS), pins=pins)
    assert [r.meta["sp_label"] for r in ReferenceRegistry(tmp_path).list()] == ["00"]


def test_statepoint_values():
    s = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML).with_coolant_statepoints(BENCH)
    v = s.statepoint_values("80")
    assert v["coolant_density"] == s.get_branch("coolant_density").values[2]
    assert v["fuel_temperature"] == 900.0
    assert s.statepoint_values({"coolant_density": 0.4570232})["coolant_density"] == \
        pytest.approx(0.4570232, rel=1e-6)
    with pytest.raises(ValueError, match="must give a value"):
        s.statepoint_values()
    with pytest.raises(KeyError):
        s.statepoint_values("99")
