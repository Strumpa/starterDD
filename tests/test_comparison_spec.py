# Tests for PostTreatment.spec (plan P4): observables and naming, DRAGON
# outputs (COND as group indices of the scheme output mesh), Serpent detector
# generation with a written manifest, request deduplication and coverage.

from pathlib import Path

import numpy as np
import pytest

from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
from starterDD.PostTreatment import DetectorManifest, PinLattice, physics_data as phys
from starterDD.PostTreatment.spec import (
    ComparisonSpec, Observable, dragon_cond, find_output, output_observable,
    scheme_output_mesh,
)

from conftest import (GE14_CALC_SCHEME_1L_YAML, GE14_COMPOSITIONS_YAML,
                      GE14_DOM_GEOMETRY_YAML, GE14_MINI_DIR)

GE14_2L = str(Path(GE14_CALC_SCHEME_1L_YAML).with_name("CALC_SCHEME_2L.yaml"))
ISOS = ["U235", "U238", "Gd155", "Gd157"]


# ── naming and meshes ───────────────────────────────────────────────────────

def test_observable_names_and_validation():
    assert Observable("by_pin", "2g", ISOS).name == "PIN_2G"
    assert Observable("FUEL", "1g", ISOS, serpent=()).name == "FUEL_1G"
    assert Observable("ALL", "295g", ISOS, serpent=("flux",)).name == "ALL_295G"
    custom = Observable("PIN", [1e-5, 0.625, 1e3, 2e7], ISOS, tag="X")
    assert custom.name == "PIN_C3GX" and custom.n_groups == 3
    with pytest.raises(ValueError, match="12 characters"):
        Observable("PIN", "2g", ISOS, tag="_TOO_LONG_TAG")
    with pytest.raises(ValueError, match="macroscopic rates"):
        Observable("ALL", "2g", ISOS)                     # ALL per-isotope rr not scorable
    with pytest.raises(ValueError, match="pins="):
        Observable("FUEL", "2g", ISOS, pins=[1])
    with pytest.raises(ValueError, match="Duplicate"):
        ComparisonSpec([Observable("PIN", "2g", ISOS), Observable("PIN", "2g", ["U235"])])


def test_scheme_output_mesh_1l_and_2l():
    s1 = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    s2 = DragonCalculationScheme.from_yaml(GE14_2L)
    assert phys.bounds_match(scheme_output_mesh(s1), phys.G295_EV)
    assert phys.bounds_match(scheme_output_mesh(s2), phys.G26_EV)   # COND indices -> 26g preset


def _edione_group(e, mesh_asc):
    """Group index DRAGON EDIENE.f gives to a real COND energy e (float32
    library energies, first group whose lower bound is <= e)."""
    energy = np.asarray(mesh_asc, dtype=np.float32)[::-1]
    for j in range(len(energy) - 1):
        if np.float32(e) >= energy[j + 1]:
            return j + 1
    return len(energy) - 1


def test_dragon_cond_indices_match_edi_energy_mapping():
    out = phys.G295_EV
    assert dragon_cond(out, out) is None
    assert dragon_cond(phys.G2_EV[[0, -1]], out) == []
    idx = dragon_cond(phys.G2_EV, out)
    assert idx == [_edione_group(0.625, out)] == [261]       # legacy COND 0.625
    assert len(phys.G295_EV) - 1 - idx[0] == 34                 # 34 thermal groups
    idx26 = dragon_cond(phys.G26_EV, out)
    s2 = DragonCalculationScheme.from_yaml(GE14_2L)
    assert idx26 == s2.steps[2].energy_groups_bounds           # the scheme's own COND
    with pytest.raises(phys.MeshMismatchError):
        dragon_cond(phys.G295_EV, phys.G26_EV)


def test_to_dragon_outputs_and_generated_procedure(tmp_path):
    from test_scheme_api import GE14_CASE, _generate
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    s1 = spec.apply_to(DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML))
    outs = {o.name: o for o in s1.outputs}
    assert list(outs) == spec.names()
    assert outs["PIN_2G"].energy_bounds == [261]
    assert outs["PIN_295G"].energy_bounds is None and outs["FUEL_1G"].energy_bounds == []
    assert outs["ALL_295G"].spatial_integration_mode == "ALL"
    files = _generate(tmp_path / "p", s1, GE14_CASE)
    edir = next(c for n, c in files.items() if n.startswith("EDIR"))
    assert "COND 261" in edir and "SAVE ON PIN_2G" in edir
    # YAML round trip keeps the integer indices
    s1.to_yaml(tmp_path / "s.yaml")
    again = DragonCalculationScheme.from_yaml(tmp_path / "s.yaml")
    assert [o.energy_bounds for o in again.outputs] == [o.energy_bounds for o in s1.outputs]


def test_2l_scheme_cannot_edit_295_groups():
    s2 = DragonCalculationScheme.from_yaml(GE14_2L)
    fine = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    with pytest.raises(phys.MeshMismatchError, match="26-group"):
        fine.to_dragon_outputs(s2)
    with pytest.warns(UserWarning):
        assert [o.name for o in fine.to_dragon_outputs(s2, skip_unavailable=True)] == \
            ["FUEL_1G", "PIN_2G"]
    coarse = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, fine_mesh="26g")
    outs = {o.name: o.energy_bounds for o in coarse.to_dragon_outputs(s2)}
    assert outs == {"FUEL_1G": [], "PIN_2G": [19], "PIN_26G": None,
                    "FUEL_26G": None, "ALL_26G": None}


# ── Serpent side ────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def ge14_model_factory():
    from starterDD.DDModel.DragonModel import CartesianAssemblyModel
    from starterDD.DDModel.helpers import associate_material_to_rod_ID
    from starterDD.InterfaceToDD.Serpent2_exports import S2_Settings, Serpent2Model
    from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml

    def make():
        a = CartesianAssemblyModel(name="GE14_S2", tdt_file="dummy.tdt",
                                   geometry_description_yaml=GE14_DOM_GEOMETRY_YAML)
        a.set_rod_ID_to_material_mapping(
            associate_material_to_rod_ID(GE14_COMPOSITIONS_YAML, GE14_DOM_GEOMETRY_YAML))
        a.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                   coolant_temperature=600.0, moderator_temperature=600.0,
                                   structural_temperature=600.0)
        a.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
        a.set_material_compositions(parse_all_compositions_from_yaml(GE14_COMPOSITIONS_YAML))
        a.number_fuel_material_mixtures_by_pin()
        m = Serpent2Model(a, S2_Settings())
        m.build()
        return m
    return make


def test_serpent_requests_are_deduplicated():
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, global_totals=False)
    assert [(r.observable, r.kind) for r in spec.serpent_requests()] == [
        ("PIN_295G", "rr"), ("PIN_295G", "flux"), ("FUEL_295G", "rr"), ("ALL_295G", "flux")]
    assert spec.serpent_source("PIN_2G") == "PIN_295G"
    # Fine pin detectors restricted to some pins no longer provide PIN_2G.
    lim = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, serpent_fine_pins=[1, 11, 20],
                                global_totals=False)
    assert lim.serpent_source("PIN_2G") == "PIN_2G"
    # With the whole-assembly totals (default), PIN_2G has its own 2-group
    # detectors (isotopes + fission/nu-fission 'void'), and ALL scores the absorption.
    full = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    reqs = {(r.observable, r.kind): r for r in full.serpent_requests()}
    assert list(reqs) == [("PIN_2G", "rr"), ("PIN_295G", "rr"), ("PIN_295G", "flux"),
                          ("FUEL_295G", "rr"), ("ALL_295G", "flux"), ("ALL_295G", "rr")]
    assert reqs[("PIN_2G", "rr")].isotopes == ISOS and reqs[("PIN_2G", "rr")].macro_mts == [-6, -7]
    assert reqs[("ALL_295G", "rr")].macro_mts == [-2, -6]


def test_fuel_isotopes_inventory_and_serpent_fine_pins():
    inv = ["U235", "U238", "Pu239", "Xe135"]
    spec = ComparisonSpec.preset("pin_rates", ISOS, fuel_isotopes_inventory=inv)
    assert spec["FUEL_1G"].isotopes == inv and spec["FUEL_1G"].serpent == ()
    assert spec["PIN_2G"].isotopes == ISOS
    assert ComparisonSpec.preset("pin_rates", ISOS)["FUEL_1G"].isotopes == ISOS   # default
    with pytest.raises(TypeError):
        ComparisonSpec.preset("pin_rates", ISOS, keff_isotopes=inv)
    with pytest.raises(ValueError, match="only applies"):
        ComparisonSpec.preset("integrated", ISOS, serpent_fine_pins="diagonal")
    with pytest.raises(ValueError, match="Unknown pin selector"):
        Observable("PIN", "26g", ISOS, pins="edge")
    diag = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, serpent_fine_pins="diagonal")
    assert diag["PIN_295G"].pins == "axis" and diag.pins_unresolved
    with pytest.raises(ValueError, match="resolve_pins"):
        diag.serpent_requests()
    lat = PinLattice.from_inputs(GE14_DOM_GEOMETRY_YAML, GE14_COMPOSITIONS_YAML)
    res = diag.resolve_pins(lat)
    assert res["PIN_295G"].pins == sorted(lat.diagonal_pins("axis"))
    assert all(lat.on_axis(p) for p in res["PIN_295G"].pins)
    assert res["PIN_2G"].pins is None and not res.pins_unresolved
    anti = Observable("PIN", "26g", ISOS, pins="anti").resolve_pins(lat)
    assert anti.pins == sorted(lat.diagonal_pins("anti"))


def test_configure_serpent_writes_a_manifest_equal_to_the_inferred_one(ge14_model_factory, tmp_path):
    model = ge14_model_factory()
    spec = ComparisonSpec.preset("pin_rates_and_spectra", ISOS)
    manifest = spec.configure_serpent(model, "endfb8r1")
    serp = tmp_path / "GE14.serp"
    model.write(str(serp))
    text = serp.read_text()
    assert text.count("\nene 295g ") == 1                       # one grid card
    assert "det PIN_295G_RR_p15" in text and "det ALL_295G_FLX" in text
    assert "det PIN_2G_RR_p15" in text and "det ALL_295G_RR" in text   # global totals
    written = DetectorManifest.read(tmp_path / "GE14.serp.detectors.json")
    inferred = DetectorManifest.infer(serp)
    assert len(written) == len(inferred) == len(manifest) == 51 + 51 + 51 + 1 + 1 + 1
    for e in written:
        i = inferred.by_detector(e.detector)
        assert (i.observable, i.kind, i.spatial, i.pin_idx, i.materials, i.responses, i.grid) == \
            (e.observable, e.kind, e.spatial, e.pin_idx, e.materials, e.responses, e.grid)
        np.testing.assert_allclose(i.bounds_eV, e.bounds_eV, rtol=1e-8)
        ref_mesh = phys.G2_EV if e.observable == "PIN_2G" else phys.G295_EV
        np.testing.assert_allclose(e.bounds_eV, ref_mesh, rtol=1e-12)        # DRAGON meshes
    u238 = [r.mt for r in written.find("PIN_295G", "rr", 1).responses if r.nuclide == "U238"]
    assert tuple(u238) == phys.serpent_detector_mts("endfb8r1", "U238")
    assert spec.coverage(written) == {}


def test_pins_option_and_flux_helper(ge14_model_factory):
    model = ge14_model_factory()
    spec = ComparisonSpec([Observable("PIN", "26g", ISOS, serpent=("rr", "flux"), pins=[1, 28])])
    man = spec.configure_serpent(model, "endfb8r1")
    assert man.pins("PIN_26G", "rr") == [1, 28]
    model.add_pin_flux_detectors("26g", pins=[2])
    assert [g.name for g in model.energy_grids] == ["26g"]
    with pytest.raises(ValueError, match="different groups"):
        model.ensure_energy_grid("26g", [1e-10, 1.0, 20.0])
    # a selector is resolved on the model's lattice when the detectors are written
    model2 = ge14_model_factory()
    diag = ComparisonSpec([Observable("PIN", "26g", ISOS, serpent=("flux",), pins="diagonal")])
    man2 = diag.configure_serpent(model2, "endfb8r1")
    lat = PinLattice.from_assembly(model2.assembly)
    assert man2.pins("PIN_26G", "flux") == sorted(lat.diagonal_pins("axis"))


def test_coverage_of_the_legacy_reference():
    """The committed legacy run (4 pins, legacy 295g grid): 2G and 26G are
    covered on its pins, 295G is not (wrong grid), missing pins reported."""
    serp = GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.gz"
    man = DetectorManifest.read(GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.detectors.json")
    pins = man.pins("PIN_295G", "rr")
    spec26 = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, fine_mesh="26g",
                                   global_totals=False)
    assert spec26.coverage(man, pins=pins) == {}
    gaps = spec26.coverage(man, n_pins=51)
    assert "no detector for pins" in gaps["PIN_26G/rr"][0]
    # the legacy run has no 'void' responses: the whole-assembly totals are missing
    tot = ComparisonSpec.preset("pin_rates_and_spectra", ISOS, fine_mesh="26g").coverage(
        man, pins=pins)
    assert set(tot) == {"PIN_2G/rr", "ALL_26G/rr"}
    assert "no 'void' response" in tot["PIN_2G/rr"][0]
    fine = ComparisonSpec.preset("pin_rates_and_spectra", ISOS,
                                 global_totals=False).coverage(man, pins=pins)
    assert set(fine) == {"PIN_295G/rr", "PIN_295G/flux", "FUEL_295G/rr", "ALL_295G/flux"}
    assert "legacy 295g grid" in fine["PIN_295G/rr"][0]
    assert serp.exists()


def test_output_observable_and_legacy_aliases():
    from starterDD.PostTreatment import CompoReader
    with CompoReader(GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz", cache=False) as c:
        assert output_observable(c, "EDIR_2G") == "PIN_2G"
        assert output_observable(c, "EDIHOM_COND") == "FUEL_1G"
        assert find_output(c, "PIN_2G") == "EDIR_2G"
        with pytest.raises(KeyError):
            find_output(c, "PIN_26G")


def test_union_merges_isotopes():
    a = ComparisonSpec([Observable("PIN", "2g", ["U235"])])
    b = ComparisonSpec([Observable("PIN", "2g", ["U238"], serpent=("rr", "flux")),
                        Observable("FUEL", "1g", ISOS, serpent=())])
    u = ComparisonSpec.union(a, b)
    assert u.names() == ["PIN_2G", "FUEL_1G"]
    assert u["PIN_2G"].isotopes == ["U235", "U238"] and u["PIN_2G"].serpent == ("rr", "flux")
