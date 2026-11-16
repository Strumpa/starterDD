# Macroscopic (all-nuclide) observables: Serpent 'dr <mt> void' responses,
# the "integrated" default preset, DRAGON-side sums over the edited isotopes
# plus *MAC*RES, from_scheme(), and the DRAGON gap temperature in Serpent.

import warnings
from pathlib import Path

import numpy as np
import pytest

from starterDD.DDModel.DragonCalculationScheme import DragonCalculationScheme
from starterDD.PostTreatment import (CompoReader, DetectorManifest, PinLattice, SerpentRun,
                                     compute_macro_rates, physics_data as phys)
from starterDD.PostTreatment.spec import ComparisonSpec, Observable, cond_to_bounds

from conftest import (GE14_CALC_SCHEME_1L_YAML, GE14_COMPOSITIONS_YAML, GE14_MINI_DIR)

GE14_ASSEMBLY = Path(GE14_CALC_SCHEME_1L_YAML).parent
GE14_2L = str(GE14_ASSEMBLY / "CALC_SCHEME_2L.yaml")
ISOS = ["U234", "U235", "U236", "U238", "Gd155", "Gd157"]
MINI_CPO = GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz"


# ── spec ────────────────────────────────────────────────────────────────────

def test_integrated_preset_is_the_default_and_small():
    spec = ComparisonSpec.default(ISOS)
    assert spec.names() == ["FUEL_1G", "PIN_2G", "FUEL_295G", "ALL_295G"]
    reqs = {(r.observable, r.kind): r for r in spec.serpent_requests()}
    assert set(reqs) == {("PIN_2G", "rr"), ("FUEL_295G", "rr"), ("FUEL_295G", "flux"),
                         ("ALL_295G", "flux"), ("ALL_295G", "rr")}
    # fission + nu-fission per pin, assembly absorption: whole-assembly normalisations
    assert reqs[("PIN_2G", "rr")].isotopes == [] and reqs[("PIN_2G", "rr")].macro_mts == [-6, -7]
    assert reqs[("FUEL_295G", "rr")].macro_mts == [-2, -6]
    assert reqs[("ALL_295G", "rr")].isotopes == [] and reqs[("ALL_295G", "rr")].macro_mts == [-2, -6]
    small = ComparisonSpec.default(ISOS, global_totals=False)
    assert [r.macro_mts for r in small.serpent_requests() if r.observable == "PIN_2G"] == [[-6]]
    assert ("ALL_295G", "rr") not in {(r.observable, r.kind) for r in small.serpent_requests()}
    assert spec["PIN_2G"].isotopes == ISOS                     # DRAGON edit (MICR)
    with pytest.raises(ValueError, match="Unknown macroscopic"):
        Observable("PIN", "2g", ISOS, macro=("heat",))
    with pytest.raises(ValueError, match="needs isotopes"):
        Observable("PIN", "2g", ISOS, serpent_isotopes=())
    # ALL may score macroscopic rates (no domain), not per-isotope ones
    Observable("ALL", "2g", ISOS, macro=("fission",), serpent_isotopes=())


def test_integrated_detectors_and_manifest(ge14_serpent_model, tmp_path):
    from copy import deepcopy
    model = deepcopy(ge14_serpent_model)
    spec = ComparisonSpec.default(["U235", "U238"])
    spec.configure_serpent(model, "endfb8r1")
    serp = tmp_path / "GE14.serp"
    model.write(str(serp))
    text = serp.read_text()
    assert "    dr -6  void" in text and "    dr -2  void" in text
    assert "mat U235" not in text                               # no response materials needed
    written = DetectorManifest.read(tmp_path / "GE14.serp.detectors.json")
    inferred = DetectorManifest.infer(serp)
    assert len(written) == 51 + 1 + 1 + 1 + 1          # + ALL_295G_RR (assembly absorption)
    for e in written:
        assert inferred.by_detector(e.detector).responses == e.responses
    r = written.find("PIN_2G", "rr", 15).responses
    assert [(x.mt, x.material, x.nuclide) for x in r] == [(-6, "void", None), (-7, "void", None)]
    r = written.find("ALL_295G", "rr", None)
    assert r.materials == [] and [x.mt for x in r.responses] == [-2, -6]
    assert spec.coverage(written) == {}
    legacy = DetectorManifest.read(GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.detectors.json")
    gaps = spec.coverage(legacy)
    assert "no 'void' response" in gaps["PIN_2G/rr"][0]


def test_cond_to_bounds_and_from_scheme():
    s1 = DragonCalculationScheme.from_yaml(GE14_CALC_SCHEME_1L_YAML)
    out = phys.G295_EV
    assert np.allclose(cond_to_bounds([0.625], out), phys.G2_EV, rtol=1e-6)   # EDIENE rule
    assert np.allclose(cond_to_bounds([261], out), phys.G2_EV, rtol=1e-6)
    assert len(cond_to_bounds([], out)) == 2 and cond_to_bounds(None, out) is not out
    spec = ComparisonSpec.from_scheme(s1, macro=("fission",), global_totals=False)
    by_out = {o.dragon_name: o for o in spec}
    assert {k: v.name for k, v in by_out.items()} == {
        "EDIHOM_COND": "FUEL_1G", "EDIHOM_295": "ALL_295G",
        "EDIR_2G": "PIN_2G", "U238_295": "FUEL_295G"}
    assert by_out["EDIHOM_COND"].serpent == () and by_out["EDIHOM_295"].serpent == ("flux",)
    assert by_out["EDIR_2G"].macro == ("fission",) and by_out["U238_295"].serpent == ("rr",)
    assert by_out["EDIR_2G"].micro_isotopes == by_out["EDIR_2G"].isotopes
    # default: whole-assembly totals on the coarsest PIN output and on ALL (Serpent only)
    tot = {o.dragon_name: o for o in ComparisonSpec.from_scheme(s1)}
    assert tot["EDIR_2G"].macro == ("fission", "nufission")
    assert tot["EDIHOM_295"].serpent == ("flux", "rr") and tot["EDIHOM_295"].macro == ("absorption",)
    assert tot["EDIHOM_295"].micro_isotopes == [] and tot["U238_295"].macro == ()
    # The spec reproduces the scheme's outputs (names, modes, isotopes).
    regen = {o.name: o for o in spec.to_dragon_outputs(s1)}
    for o in s1.outputs:
        assert regen[o.name].spatial_integration_mode == o.spatial_integration_mode
        assert regen[o.name].isotopes == o.isotopes
    s2 = DragonCalculationScheme.from_yaml(GE14_2L)
    names2 = {o.dragon_name: o.name for o in ComparisonSpec.from_scheme(s2)}
    assert names2 == {"EDIHOM_COND": "FUEL_1G", "EDIHOM_295": "ALL_26G",   # 26 groups in 2L
                      "EDIR_2G": "PIN_2G", "U238_295": "FUEL_26G"}


# ── DRAGON + Serpent macroscopic rates ──────────────────────────────────────

@pytest.fixture(scope="module")
def ge14_serpent_model():
    from starterDD.DDModel.DragonModel import CartesianAssemblyModel
    from starterDD.DDModel.helpers import associate_material_to_rod_ID
    from starterDD.InterfaceToDD.Serpent2_exports import S2_Settings, Serpent2Model
    from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
    geom = str(GE14_ASSEMBLY / "GEOM_GE14_DOM_DIAG.yaml")
    a = CartesianAssemblyModel(name="GE14", tdt_file="", geometry_description_yaml=geom)
    a.set_rod_ID_to_material_mapping(associate_material_to_rod_ID(GE14_COMPOSITIONS_YAML, geom))
    a.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=900.0,
                               coolant_temperature=600.0, moderator_temperature=600.0,
                               structural_temperature=600.0)
    a.analyze_lattice_description(build_pins=True, apply_self_shielding="from_yaml")
    a.set_material_compositions(parse_all_compositions_from_yaml(GE14_COMPOSITIONS_YAML))
    a.number_fuel_material_mixtures_by_pin()
    m = Serpent2Model(a, S2_Settings())
    m.build()
    return m


def _dragon_macro(c, output, mix, keys):
    """Independent sum: ISOTOPESDENS x records, read directly."""
    p = f"{output}/MIXTURES/{mix}/CALCULATIONS/1"
    dens = np.atleast_1d(c.get(p + "/ISOTOPESDENS"))
    tot = 0.0
    for i, n in enumerate(dens, start=1):
        for k in keys:
            v = c.get(f"{p}/ISOTOPESLIST/{i}/{k}", None)
            if v is not None:
                tot = tot + n * np.atleast_1d(np.asarray(v, dtype=float))[::-1]
    phi = np.atleast_1d(np.asarray(c.get(f"{p}/ISOTOPESLIST/1/NWT0"), dtype=float))[::-1]
    return tot * phi * float(c.get(p + "/MIXTURESVOL"))


def _write_det0(path, blocks, bounds_eV, grids=None):
    """Minimal det0.m: {detector: values (n_reactions, n_groups)}, rel err 1e-4.
    *grids*: ``{detector: bounds_eV}`` for detectors on another grid."""
    lines = []
    for name, vals in blocks.items():
        e = np.asarray((grids or {}).get(name, bounds_eV)) * 1e-6
        lines.append(f"DET{name} = [")
        i = 0
        for r, row in enumerate(vals):
            for g, v in enumerate(row):
                i += 1
                lines.append(f" {i} {g + 1} 1 1 1 1 {r + 1} 1 1 1 {v:.10E} 0.00010")
        lines.append("];")
        lines.append(f"DET{name}E = [")
        lines += [f" {e[g]:.8E} {e[g + 1]:.8E} {0.5 * (e[g] + e[g + 1]):.8E}"
                  for g in range(len(e) - 1)]
        lines.append("];")
    Path(path).write_text("\n".join(lines) + "\n")


def test_macro_rates_pin_fission_end_to_end(ge14_serpent_model, tmp_path):
    """DRAGON side (isotopes + *MAC*RES) against an independent sum, and the
    Serpent side through a written input and a det0.m with known values."""
    lat = PinLattice.from_inputs(GE14_ASSEMBLY / "GEOM_GE14_DOM_DIAG.yaml",
                                 GE14_COMPOSITIONS_YAML)
    from copy import deepcopy
    model = deepcopy(ge14_serpent_model)
    spec = ComparisonSpec([Observable("PIN", "2g", ISOS, macro=("fission",),
                                      serpent_isotopes=())])
    spec.configure_serpent(model, "endfb8r1")
    serp = tmp_path / "GE14.serp"
    model.write(str(serp))
    with CompoReader(MINI_CPO, cache=False) as c:
        d = {p: _dragon_macro(c, "EDIR_2G", lat.mixture_of(p), ["NFTOT"]) / lat.dragon_fraction(p)
             for p in lat.ordered_pins}
        assert "*MAC*RES" in c.isotopes("EDIR_2G", 1)
        delta = {p: (0.01 if p == 15 else 0.0) for p in lat.ordered_pins}
        scale = 3.7
        _write_det0(str(serp) + "_det0.m",
                    {f"PIN_2G_RR_p{p}": (d[p] * lat.multiplicity(p) * scale * (1 + delta[p]))[None, :]
                     for p in lat.ordered_pins}, phys.G2_EV)
        run = SerpentRun(serp)                                   # exporter manifest found
        cmp = compute_macro_rates(c, "EDIR_2G", run, "PIN_2G", lattice=lat)
    assert cmp.quantities == ["fission"]
    for k, p in enumerate(cmp.regions):
        np.testing.assert_allclose(cmp.dragon["fission"][k], d[p], rtol=1e-12)
    m = np.array([lat.multiplicity(p) for p in cmp.regions], dtype=float)
    dd = np.array([d[p] for p in cmp.regions])
    dl = np.array([delta[p] for p in cmp.regions])
    f = (m[:, None] * dd * (1 + dl[:, None])).sum() * scale / (m[:, None] * dd).sum()
    expected = np.broadcast_to(100 * (f / (scale * (1 + dl[:, None])) - 1), (len(dl), 2))
    np.testing.assert_allclose(cmp.errors("fission"), expected, atol=1e-6)
    k15 = cmp.regions.index(15)
    assert cmp.errors("fission")[k15, 0] == pytest.approx(100 * (f / scale / 1.01 - 1), abs=1e-6)
    np.testing.assert_allclose(cmp.rel_sigma("fission"), 0.01, rtol=1e-3)


def synthetic_pin_blocks(c, output, run0, observable, lat, scale, bias=None):
    """det0 blocks ``{detector: (n_responses, n_groups)}`` of the per-pin
    detectors of *observable*: Serpent = DRAGON *output* x *scale* (per
    isotope x *bias*), for every per-isotope and ``void`` (-6, -7) response."""
    from starterDD.PostTreatment.comparison import _domain_composition
    bias = bias or {}
    key_of = {phys.FISSION.mt: phys.FISSION.dragon_key,
              **{r.mt: r.dragon_key for r in phys.ABSORPTION_REACTIONS}}
    ng = c.n_groups(output, 1)
    blocks = {}
    for p in lat.ordered_pins:
        e = run0.manifest.find(observable, "rr", p)
        mix, m, fr = lat.mixture_of(p), lat.multiplicity(p), lat.dragon_fraction(p)
        n_s = _domain_composition(run0, e)
        phi = c.micro_xs(output, "U235", "NWT0", 1, mixtures=mix)
        dens = c.densities(output, 1, mix)
        vol = c.isotope_volumes(output, 1, mix)
        rows = []
        for r in e.responses:
            if r.material == phys.VOID:
                key = {-6: ["NFTOT"], -7: ["NUSIGF"]}[r.mt]
                rows.append(_dragon_macro(c, output, mix, key) / fr * m * scale)
                continue
            xs = c.micro_xs(output, r.nuclide, key_of.get(r.mt, "NONE"), 1, mixtures=mix)
            if xs is None or n_s.get(r.nuclide, 0.0) == 0.0:
                rows.append(np.zeros(ng))
                continue
            rate = xs * phi * dens[r.nuclide] * vol[r.nuclide] / fr      # per physical pin
            rows.append(rate * m * scale * bias.get(r.nuclide, 1.0) / n_s[r.nuclide])
        blocks[e.detector] = np.array(rows)
    return blocks


def test_global_normalisation_does_not_depend_on_the_compared_isotopes(ge14_serpent_model,
                                                                       tmp_path):
    """Serpent = DRAGON x a known scale, U238 responses x 1.02. The
    whole-assembly factors (fission_source from PIN_2G nu-fission, absorption
    from a 1-mixture output with the DIAG geometry fraction) recover the scale
    exactly, and the per-isotope errors with that factor do not depend on the
    isotope set; power_subset does."""
    from copy import deepcopy
    from starterDD.PostTreatment import compute_pin_rates, global_normalisation
    lat =PinLattice.from_inputs(GE14_ASSEMBLY / "GEOM_GE14_DOM_DIAG.yaml",
                                 GE14_COMPOSITIONS_YAML)
    assert lat.dragon_geometry_fraction == 0.5
    model = deepcopy(ge14_serpent_model)
    # ALL_1G stands for an assembly output: EDIHOM_COND is the only 1-mixture
    # output of the mini COMPO (the geometry-fraction logic is the same).
    spec = ComparisonSpec([
        Observable("PIN", "2g", ISOS, macro=("fission", "nufission")),
        Observable("ALL", "1g", ISOS, serpent=("rr",), macro=("absorption",),
                   serpent_isotopes=(), dragon=False)])
    with pytest.raises(phys.UnknownReactionsError, match="U236"):
        spec.configure_serpent(deepcopy(ge14_serpent_model), "endfb8r1")
    # U236 is not in physics_data.SERPENT_DETECTOR_MTS: its reactions are given
    spec.configure_serpent(model, "endfb8r1", reaction_mts={"U236": [18, 102, 16, 17]})
    u236 = [r.mt for r in model.detector_manifest_entries[0].responses if r.nuclide == "U236"]
    assert u236 == [18, 102, 16, 17]
    serp = tmp_path / "GE14.serp"
    model.write(str(serp))
    scale, bias = 3.7, {"U238": 1.02}
    abs_keys = phys.MACRO_REACTIONS["absorption"].dragon_keys
    with CompoReader(MINI_CPO, cache=False) as c:
        blocks = synthetic_pin_blocks(c, "EDIR_2G", SerpentRun(serp), "PIN_2G", lat,
                                      scale, bias)
        all_abs = _dragon_macro(c, "EDIHOM_COND", 1, abs_keys)
        cap = _dragon_macro(c, "EDIHOM_COND", 1, [k for k in abs_keys if k != "NFTOT"])
        full = scale / lat.dragon_geometry_fraction        # DIAG half -> full assembly
        blocks["ALL_1G_RR"] = np.array([cap * full, (all_abs - cap) * full])   # -2, -6
        _write_det0(str(serp) + "_det0.m", blocks, phys.G2_EV,
                    grids={"ALL_1G_RR": phys.G2_EV[[0, -1]]})
        run = SerpentRun(serp)
        pin = compute_macro_rates(c, "EDIR_2G", run, "PIN_2G", ["fission", "nufission"],
                                  lattice=lat)
        n_fs = global_normalisation("fission_source", pin)
        assert n_fs.factor == pytest.approx(scale, rel=1e-9)
        assert global_normalisation("fission", pin).factor == pytest.approx(scale, rel=1e-9)
        tot = compute_macro_rates(c, "EDIHOM_COND", run, "ALL_1G", ["absorption"])
        n_abs = global_normalisation("absorption", tot, lat.dragon_geometry_fraction)
        assert n_abs.factor == pytest.approx(scale, rel=1e-9)
        assert n_abs.dragon_total == pytest.approx(all_abs.sum() / 0.5, rel=1e-12)
        with pytest.raises(Exception, match="needs a ALL"):
            global_normalisation("absorption", pin)
        u5 = compute_pin_rates(c, "EDIR_2G", run, "PIN_2G", ["U235"], lat, "endfb8r1")
        u58 = compute_pin_rates(c, "EDIR_2G", run, "PIN_2G", ["U235", "U238"], lat, "endfb8r1")
        u8 = compute_pin_rates(c, "EDIR_2G", run, "PIN_2G", ["U238"], lat, "endfb8r1")
    # whole-assembly factor: U235 exact, U238 off by its bias, whatever the set
    np.testing.assert_allclose(u5.errors("absorption", factor=n_fs.factor), 0.0, atol=1e-6)
    np.testing.assert_allclose(u8.errors("absorption", factor=n_fs.factor),
                               100 * (1 / 1.02 - 1), atol=1e-6)
    # power_subset depends on the set: U235 alone is exact, U235+U238 is not
    np.testing.assert_allclose(u5.errors("absorption", "power_subset"), 0.0, atol=1e-6)
    assert u58.normalisation_factor("power_subset") != pytest.approx(scale, rel=1e-4)


# ── gap temperature (DRAGON: MIX <<TFUEL>>) ─────────────────────────────────

def test_structural_temperatures_follow_dragon(ge14_serpent_model):
    from copy import deepcopy
    m = deepcopy(ge14_serpent_model)
    m.materials = [x for x in m.materials if "_pin_" in x.name]
    names = {"COOLANT": "coolant", "CLAD": "zr2", "GAP": "gap", "MODERATOR": "moderator",
             "CHANNEL_BOX": "zr4"}
    m.build_structural_materials_from_assembly(name_map=names)
    t = {x.name: x.temperature for x in m.materials}
    assert t["gap"] == 900.0 and t["coolant"] == 600.0 and t["zr2"] == 600.0
    m2 = deepcopy(ge14_serpent_model)
    m2.materials = [x for x in m2.materials if "_pin_" in x.name]
    with pytest.warns(UserWarning, match="GAP at 600.0 K"):
        m2.build_structural_materials_from_assembly(name_map=names, temperature_map={"GAP": 600.0})
