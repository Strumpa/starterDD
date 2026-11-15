"""P1 comparison core: lattice_maps, comparison, plotting.

Tiers as in ``tests/data/README.md``: synthetic / committed inputs (always),
``tests/data/ge14_mini`` (always), full local GE14 data (``fulldata``).

The full-data regression rebuilds the outputs of the original GE14 analysis
(``2G_data.npz``, ``density_00_<k>.txt``, summary CSV rows) for DOM/VAN x
1L/2L at 0 % void with the per-physical-pin convention of
:mod:`starterDD.PostTreatment.comparison`.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from starterDD.PostTreatment import physics_data as phys
from starterDD.PostTreatment.compo_reader import CompoReader
from starterDD.PostTreatment.comparison import (
    CodeRates, NormalisationError, PinRateComparison, ReactionCoverage,
    compare_keff, compare_spectra, compute_pin_rates,
)
from starterDD.PostTreatment.lattice_maps import PinLattice
from starterDD.PostTreatment.pylcm import Pylcm
from starterDD.PostTreatment.serpent_reader import DetectorManifest, SerpentRun

from conftest import (
    FULL_DATA_DIR, GE14_INPUTS_DIR, GE14_MATERIAL_COMPOS_00_YAML, GE14_MINI_DIR,
)

DOM_DIAG_GEOM = GE14_INPUTS_DIR / "ASSEMBLY" / "GEOM_GE14_DOM_DIAG.yaml"
DOM_AXIS = [1, 11, 20, 28, 33, 37, 42, 46, 49, 51]
ISOTOPE_SETS = [["U235"], ["U238"], ["Gd155", "Gd157"], ["U235", "U238", "Gd155", "Gd157"]]


@pytest.fixture(scope="module")
def dom_lattice():
    return PinLattice.from_inputs(DOM_DIAG_GEOM, GE14_MATERIAL_COMPOS_00_YAML,
                                  name="GE14_DOM_DIAG")


# ═════════════════════════════════════════════════════════════════════════════
# physics_data additions
# ═════════════════════════════════════════════════════════════════════════════

def test_legacy_exporter_shem295_differs_only_in_bounds_235_238():
    legacy = phys.G295_STARTERDD_LEGACY_EV
    diff = np.where(np.abs(legacy / phys.G295_EV - 1) > phys.BOUNDS_RTOL)[0]
    assert list(diff) == [235, 236, 237, 238]
    common = phys.common_mesh(phys.G295_EV, legacy)
    assert len(common) - 1 == 291
    assert phys.merged_groups(phys.G295_EV, common) == [(234, 238)]
    for coarse in (phys.G2_EV, phys.G26_EV):           # unaffected meshes
        phys.condensation_indices(phys.G295_EV, coarse)
        phys.condensation_indices(legacy, coarse)
    with pytest.raises(phys.MeshMismatchError):
        phys.condensation_indices(legacy, phys.G295_EV)
    # Legacy references keep their "295G" observable names.
    assert phys.mesh_code(legacy) == "295G" and phys.is_legacy_mesh(legacy)
    assert not phys.is_legacy_mesh(phys.G295_EV)
    assert phys.G295_DRAGON_EV is phys.G295_EV        # deprecated alias


def test_target_mesh_dragon_default_and_common_merged_ranges():
    from starterDD.PostTreatment.comparison import _target_mesh
    b, merged = _target_mesh(phys.G295_EV, phys.G295_EV, "dragon", "x")
    assert b is phys.G295_EV and merged == []
    with pytest.raises(phys.MeshMismatchError, match="legacy starterDD 295g grid"):
        _target_mesh(phys.G295_EV, phys.G295_STARTERDD_LEGACY_EV, "dragon", "x")
    with pytest.warns(UserWarning, match="57-61"):
        b, merged = _target_mesh(phys.G295_EV, phys.G295_STARTERDD_LEGACY_EV, "common", "x")
    assert len(b) - 1 == 291
    assert np.allclose(merged, [(3481.07, 9118.81)], rtol=1e-5)


# ═════════════════════════════════════════════════════════════════════════════
# PinLattice
# ═════════════════════════════════════════════════════════════════════════════

class TestPinLattice:
    def test_counts_and_axis(self, dom_lattice):
        lat = dom_lattice
        assert lat.symmetry == "main-diagonal"
        assert lat.n_pins == 51 and lat.n_physical_pins == 92
        assert lat.axis_pins == DOM_AXIS
        assert lat.diagonal_pins("axis") == DOM_AXIS == lat.diagonal_pins("main")
        assert lat.ordered_pins == list(range(1, 52))
        assert lat.mixture_of(28) == 28

    def test_multiplicity_and_dragon_fraction(self, dom_lattice):
        lat = dom_lattice
        for p in lat.ordered_pins:
            on = p in DOM_AXIS
            assert lat.multiplicity(p) == (1 if on else 2)
            assert lat.dragon_fraction(p) == (0.5 if on else 1.0)
            # per-physical-pin convention reproduces the legacy factor:
            # multiplicity * legacy factor == 2 for every pin
            assert lat.multiplicity(p) * lat.legacy_symmetry_factors()[p - 1] == 2

    def test_to_map(self, dom_lattice):
        lat = dom_lattice
        grid = lat.to_map(np.array(lat.ordered_pins, float))
        assert grid.shape == lat.shape
        assert np.count_nonzero(np.isfinite(grid)) == 92
        assert np.array_equal(grid, grid.T, equal_nan=True)      # main diagonal
        g2 = lat.to_map(np.ones((51, 2)))
        assert g2.shape == (2,) + lat.shape
        assert np.nansum(lat.to_map({28: 5.0}, fill=0.0)) == 5.0

    def test_broken_symmetry_is_detected(self, dom_lattice):
        lat = PinLattice(dom_lattice.shape, {k: list(v) for k, v in dom_lattice.positions.items()},
                         dom_lattice.ordered_pins, dom_lattice.symmetry)
        lat.positions[2] = lat.positions[2][:1]               # drop the mirror
        with pytest.raises(ValueError, match="not closed"):
            lat._check_symmetry()


# ═════════════════════════════════════════════════════════════════════════════
# Normalisation and statistics on synthetic rates
# ═════════════════════════════════════════════════════════════════════════════

def _synthetic(lat: PinLattice, seed=0, ng=2) -> PinRateComparison:
    rng = np.random.default_rng(seed)
    pins = lat.ordered_pins
    mult = lat.multiplicities().astype(float)

    # Serpent rates, and DRAGON rates = 1 % noise and an arbitrary scale (37)
    sf = rng.uniform(0.5, 1.5, (len(pins), ng))
    sa = sf * rng.uniform(1.5, 2.5, (len(pins), ng))
    spw = sf * 200.0
    df = sf * rng.normal(1, 0.01, sf.shape) * 37.0
    da = sa * rng.normal(1, 0.01, sa.shape) * 37.0
    dpw = df * 200.0
    bounds = phys.G2_EV if ng == 2 else np.geomspace(1e-5, 2e7, ng + 1)
    s = CodeRates("serpent", pins, bounds, sf, sa, spw, mult, sf * 1e-4, sa * 1e-4)
    d = CodeRates("dragon", pins, bounds, df, da, dpw, mult)
    cov = {"X": ReactionCoverage("X", 1, fission=True)}
    return PinRateComparison(d, s, lat, ["X"], cov)


def _legacy_map_stats(lat, err, dist):
    """Statistics as computed by the original analysis on unfolded 2D maps."""
    out = []
    for g in range(err.shape[1]):
        e = lat.to_map(err[:, g])
        d = lat.to_map(dist[:, g])
        out.append({"min_abs": np.nanmin(np.abs(e)), "min": np.nanmin(e),
                    "max": np.nanmax(e), "mean": np.nanmean(e),
                    "weighted_mean": np.nanmean(e * d),
                    "rms": np.sqrt(np.nanmean(e ** 2)),
                    "weighted_rms": np.sqrt(np.nanmean(e ** 2 * d))})
    return out


class TestNormalisation:
    def test_errors_invariant_to_code_scale(self, dom_lattice):
        c = _synthetic(dom_lattice)
        e0 = c.errors("absorption")
        c.serpent.absorption *= 7.0; c.serpent.fission *= 7.0; c.serpent.power *= 7.0
        assert np.allclose(c.errors("absorption"), e0, rtol=1e-12, atol=1e-12)
        assert np.allclose(c.errors("fission_over_absorption"),
                           _synthetic(dom_lattice).errors("fission_over_absorption"))

    def test_power_normalisation_matches_totals(self, dom_lattice):
        c = _synthetic(dom_lattice)
        k = c.normalisation_factor("power")
        assert (c.dragon.total("power").sum() * k) == pytest.approx(c.serpent.total("power").sum())
        assert c.normalisation_factor("absorption") != pytest.approx(k)

    def test_distribution_mean_one_over_physical_pins(self, dom_lattice):
        c = _synthetic(dom_lattice)
        dist = c.distribution("absorption")
        grid = dom_lattice.to_map(dist[:, 0])
        assert np.nanmean(grid) == pytest.approx(1.0)

    def test_stats_equal_legacy_map_formulas(self, dom_lattice):
        c = _synthetic(dom_lattice, seed=3)
        for q in ("absorption", "fission"):
            new = c.stats(q)
            old = _legacy_map_stats(dom_lattice, c.errors(q), c.distribution(q))
            for a, b in zip(new, old):
                for key in b:
                    assert a[key] == pytest.approx(b[key], rel=1e-12, abs=1e-14), (q, key)

    def test_partial_coverage_cannot_be_normalised(self, dom_lattice):
        c = _synthetic(dom_lattice)
        sub = PinRateComparison(
            *[CodeRates(r.code, r.pins[:3], r.bounds_eV, r.fission[:3], r.absorption[:3],
                        r.power[:3], r.multiplicity[:3]) for r in (c.dragon, c.serpent)],
            dom_lattice, c.isotopes, c.coverage)
        with pytest.raises(NormalisationError):
            sub.errors()
        assert sub.errors("fission_over_absorption").shape == (3, 2)

    def test_plots_smoke(self, dom_lattice, tmp_path):
        pytest.importorskip("matplotlib")
        from starterDD.PostTreatment import plotting
        c = _synthetic(dom_lattice)
        plotting.plot_rate_maps(c, tmp_path / "maps.png")
        plotting.plot_profile(c, dom_lattice.diagonal_pins("axis"), path=tmp_path / "p.pdf")
        assert (tmp_path / "maps.png").stat().st_size > 0
        assert (tmp_path / "p.pdf").stat().st_size > 0


# ═════════════════════════════════════════════════════════════════════════════
# Tier 2: ge14_mini
# ═════════════════════════════════════════════════════════════════════════════

MINI_PINS = [1, 2, 15, 28]


@pytest.fixture(scope="module")
def mini():
    compo = CompoReader(GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz", cache=False)
    serp = SerpentRun(GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.gz")
    golden = np.load(GE14_MINI_DIR / "golden" / "2G_data.npz")
    return compo, serp, golden


class TestGE14MiniComparison:
    def test_reproduces_golden_npz_columns(self, mini, dom_lattice):
        compo, serp, golden = mini
        c = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ISOTOPE_SETS[3],
                              dom_lattice, "endfb8r1", pins=MINI_PINS,
                              serpent_densities=dom_lattice.composition,
                              error_mode="legacy")
        assert not c.is_complete
        new = c.to_legacy_npz()
        assert "serp_pw" not in new                        # needs the whole assembly
        cols = np.array(MINI_PINS) - 1
        for key in ("serp_groups_fiss_rate_tot", "serp_groups_fiss_rate_err_tot",
                    "serp_groups_abs_rate_tot", "serp_groups_abs_rate_err_tot",
                    "drag_groups_fiss_rate_tot", "drag_groups_abs_rate_tot"):
            assert np.allclose(new[key][:, cols], golden[key][:, cols], rtol=1e-12, atol=0), key
        assert np.array_equal(new["symmetry_factors"], golden["symmetry_factors"])
        with pytest.raises(NormalisationError):
            c.errors("absorption")

    def test_serpent_mat_densities_are_equivalent(self, mini, dom_lattice):
        compo, serp, _ = mini
        kw = dict(pins=MINI_PINS)
        a = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ISOTOPE_SETS[3],
                              dom_lattice, "endfb8r1", **kw)
        b = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ISOTOPE_SETS[3],
                              dom_lattice, "endfb8r1", serpent_densities=dom_lattice.composition, **kw)
        # Serpent mat cards carry 7 significant digits of the benchmark values
        assert np.allclose(a.serpent.absorption, b.serpent.absorption, rtol=1e-6)
        assert np.array_equal(a.dragon.absorption, b.dragon.absorption)

    def test_coverage_report(self, mini, dom_lattice):
        compo, serp, _ = mini
        c = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ISOTOPE_SETS[3],
                              dom_lattice, "endfb8r1", pins=MINI_PINS)
        cov = c.coverage
        assert cov["U235"].fission and cov["U235"].compared == [102, 107, 16, 17]
        assert 103 in cov["U235"].dragon_only                # (n,p) not scored in Serpent
        assert cov["U238"].compared == [102, 103, 107, 16, 17]
        assert not cov["Gd155"].fission and cov["Gd155"].reference_pin == 15
        assert c.serpent.multiplicity.tolist() == [1, 2, 2, 1]

    def test_quadrature_errors_smaller_than_linear(self, mini, dom_lattice):
        compo, serp, _ = mini
        args = (compo, "EDIR_2G", serp, "PIN_295G", ISOTOPE_SETS[1], dom_lattice, "endfb8r1")
        q = compute_pin_rates(*args, pins=MINI_PINS)
        lin = compute_pin_rates(*args, pins=MINI_PINS, error_mode="linear")
        assert np.all(q.serpent.absorption_err <= lin.serpent.absorption_err)
        assert np.all(q.serpent.absorption_err > 0)

    def test_pin_spectra(self, mini, dom_lattice):
        compo, serp, _ = mini
        s = compare_spectra(compo, "EDIR_2G", serp, "PIN_295G", pins=MINI_PINS, lattice=dom_lattice)
        assert s.dragon.shape == s.serpent.shape == (4, 2)
        assert np.allclose(s.dragon.sum(1), 1) and np.allclose(s.serpent.sum(1), 1)
        assert np.all(np.abs(s.errors()) < 5.0)

    def test_keff(self, mini):
        compo, serp, _ = mini
        k = compare_keff(compo, "EDIHOM_COND", serp)
        assert k.delta_pcm == pytest.approx(-117.1511, abs=1e-3)
        assert k.sigma_pcm == pytest.approx(0.990261 * 7.2e-6 * 1e5)

    def test_multiplicity_mismatch_is_detected(self, mini, dom_lattice):
        compo, serp, _ = mini
        lat = PinLattice(dom_lattice.shape, {k: list(v) for k, v in dom_lattice.positions.items()},
                         dom_lattice.ordered_pins, dom_lattice.symmetry, dom_lattice.material)
        lat.positions[1] = lat.positions[1] * 2   # model claims 2 positions for pin 1
        with pytest.raises(ValueError, match="Serpent detector covers 1"):
            compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ["U235"], lat,
                              "endfb8r1", pins=[1])


# ═════════════════════════════════════════════════════════════════════════════
# Tier 3: full regression against the original analysis
# ═════════════════════════════════════════════════════════════════════════════

CASES = {
    "DOM_1L": ("DOM", "1L_PT_3x3_HFP", "PT_IC_1L_MOC", "2026-07-19T10-46-27"),
    "DOM_2L": ("DOM", "2L_PT_3x3_IC_3x3_6x6_HFP", "PT_IC_2L_IC_MOC", "2026-07-21T16-23-54"),
    "VAN_1L": ("VAN", "1L_PT_3x3_HFP", "PT_IC_1L_MOC", "2026-07-19T10-47-53"),
    "VAN_2L": ("VAN", "2L_PT_3x3_IC_3x3_6x6_HFP", "PT_IC_2L_IC_MOC", "2026-07-21T16-25-56"),
}


def _case_paths(asm, case, desc, ts):
    root = FULL_DATA_DIR / f"GE14_{asm}"
    base = root / f"results/coolant_density_{case}/GE14_{asm}_DIAG/{desc}"
    return root, base, base / ts


def _numbers(txt):
    return np.array([float(x) for x in re.findall(r"(-?\d+\.\d+(?:e[-+]?\d+)?) %", txt)])


@pytest.mark.fulldata
@pytest.mark.parametrize("case", list(CASES))
def test_full_regression_against_legacy_analysis(case):
    asm, cname, desc, ts = CASES[case]
    root, base, run = _case_paths(asm, cname, desc, ts)
    serp_path = root / "Serpent" / f"GE14_{asm}_00_endfb8r1.serp"
    if not (run / f"_CPO_GE14_{asm}_DIAG").exists() or not serp_path.exists():
        pytest.skip(f"GE14 {case} data not present")
    lat = PinLattice.from_inputs(run / "inputs" / "GEOM_DIAG.yaml",
                                 next((run / "inputs").glob("material_compos_*.yaml")))
    compo = CompoReader(Pylcm.from_ascii(run / f"_CPO_GE14_{asm}_DIAG",
                                         only=["EDIHOM_COND", "EDIR_2G"]))
    serp = SerpentRun(serp_path, manifest=DetectorManifest.infer(serp_path))
    k = compare_keff(compo, "EDIHOM_COND", serp)
    csv_rows = [ln.split(";") for ln in next(base.glob("sortie_*.csv")).read_text().splitlines()
                if ln and abs(float(ln.split(";")[1]) - k.dragon) < 1e-12]
    assert len(csv_rows) == 4
    for kset, isos in enumerate(ISOTOPE_SETS):
        c = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", isos, lat, "endfb8r1",
                              serpent_densities=lat.composition, error_mode="legacy")
        assert c.is_complete
        old = (run / f"density_00_{kset}.txt").read_text()
        new = c.legacy_stats_text()
        assert [l.split(":")[0] for l in new.splitlines()] == [l.split(":")[0] for l in old.splitlines()]
        assert np.allclose(_numbers(new), _numbers(old), rtol=1e-10, atol=1e-11)
        row = csv_rows[kset]
        assert float(row[3]) == pytest.approx(k.delta_pcm, abs=1e-8)
        assert np.allclose(c.legacy_csv_fields(), [float(x) for x in row[4:]],
                           rtol=1e-10, atol=1e-11)
        # Same % with the physical (mat-card) densities, within 1e-4 points
        if kset == 3:
            golden = np.load(run / "2G_data.npz")
            new_npz = c.to_legacy_npz()
            for key in golden.files:
                g, n = golden[key], new_npz[key]
                if np.asarray(g).dtype.kind == "f":
                    assert np.allclose(n, g, rtol=1e-12, atol=0), key
                else:
                    assert np.array_equal(n, g), key
            phys_dens = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", isos, lat,
                                          "endfb8r1")
            assert np.nanmax(np.abs(phys_dens.errors() - c.errors())) < 1e-4


@pytest.mark.fulldata
def test_full_fine_group_mesh_mismatch_is_handled():
    root, _, run = _case_paths(*CASES["DOM_1L"])
    serp_path = root / "Serpent" / "GE14_DOM_00_endfb8r1.serp"
    if not (run / "_CPO_GE14_DOM_DIAG").exists() or not serp_path.exists():
        pytest.skip("GE14 DOM data not present")
    compo = CompoReader(Pylcm.from_ascii(run / "_CPO_GE14_DOM_DIAG", only=["EDIHOM_295"]))
    serp = SerpentRun(serp_path, manifest=DetectorManifest.infer(serp_path))
    # DRAGON library mesh == canonical 295g, to the 9 digits of the COMPO.
    assert np.allclose(compo.energy_bounds("EDIHOM_295"), phys.G295_EV, rtol=1e-8, atol=0)
    with pytest.raises(phys.MeshMismatchError, match="legacy"):
        compare_spectra(compo, "EDIHOM_295", serp, "ALL_295G")
    with pytest.warns(UserWarning, match="57-61"):
        s = compare_spectra(compo, "EDIHOM_295", serp, "ALL_295G", mesh="common")
    assert s.dragon.shape == (1, 291)
    assert np.allclose(s.merged_eV, [(3481.07, 9118.81)], rtol=1e-5)


# ═════════════════════════════════════════════════════════════════════════════
# VAN geometry (committed) and the legacy report wrapper
# ═════════════════════════════════════════════════════════════════════════════

VAN_DIAG_GEOM = GE14_INPUTS_DIR / "ASSEMBLY" / "GEOM_GE14_VAN_DIAG.yaml"


def test_van_lattice_from_committed_geometry():
    lat = PinLattice.from_inputs(VAN_DIAG_GEOM, GE14_MATERIAL_COMPOS_00_YAML, name="GE14_VAN_DIAG")
    assert lat.symmetry == "main-diagonal"
    assert lat.n_pins == 42 and lat.n_physical_pins == 78
    # Same list as the original driver's hardcoded VAN "anti_diag" pins
    assert lat.diagonal_pins("axis") == lat.axis_pins == [1, 16, 24, 35, 38, 42]
    assert sum(lat.multiplicity(p) for p in lat.ordered_pins) == 78


def test_report_output_selection_by_content(tmp_path):
    from test_post_treatment_readers import N_MIX, write_synthetic_compo
    from starterDD.PostTreatment.reports import _choose_outputs
    compo = CompoReader(write_synthetic_compo(tmp_path / "_CPO"), cache=False)
    # HOM_1G: 1 mixture, OUT_2G: one mixture per pin
    assert _choose_outputs(compo, N_MIX, sp=1) == ("HOM_1G", None, "OUT_2G", None)


@pytest.mark.fulldata
@pytest.mark.parametrize("case", ["DOM_1L", "VAN_2L"])
def test_full_legacy_report_wrapper(case, tmp_path):
    from starterDD.PostTreatment.reports import legacy_report_for_run
    asm, cname, desc, ts = CASES[case]
    root, base, run = _case_paths(asm, cname, desc, ts)
    serp_path = root / "Serpent" / f"GE14_{asm}_00_endfb8r1.serp"
    if not (run / f"_CPO_GE14_{asm}_DIAG").exists() or not serp_path.exists():
        pytest.skip(f"GE14 {case} data not present")
    csv = tmp_path / "summary.csv"
    rep = legacy_report_for_run(run, serp_path, case_name="density_00", out_dir=tmp_path,
                                csv_path=csv, cache=tmp_path / "cpo.h5",
                                figures=(case == "VAN_2L"))
    assert rep.fine_rates is not None
    assert rep.fine_rates.output == ("EDIR_295" if "1L" in case else "EDIR_26")
    for k in range(4):
        assert np.allclose(_numbers((tmp_path / f"density_00_{k}.txt").read_text()),
                           _numbers((run / f"density_00_{k}.txt").read_text()),
                           rtol=1e-10, atol=1e-11)
    old = [l.split(";") for l in next(base.glob("sortie_*.csv")).read_text().splitlines()
           if l.split(";")[1] == str(rep.keff.dragon)]
    new = [l.split(";") for l in csv.read_text().splitlines()]
    assert len(old) == len(new) == 4
    for a, b in zip(old, new):
        assert a[1:4] == b[1:4]                                  # keff, keff, pcm: same text
        assert float(a[0]) == pytest.approx(float(b[0]), abs=0.01)   # wall time (rounded in manifest)
        assert np.allclose([float(x) for x in a[4:]], [float(x) for x in b[4:]],
                           rtol=1e-10, atol=1e-11)
    g, n = np.load(run / "2G_data.npz"), np.load(tmp_path / "2G_data.npz")
    for key in g.files:
        if g[key].dtype.kind == "f":
            assert np.allclose(n[key], g[key], rtol=1e-12, atol=0), key
        else:
            assert np.array_equal(n[key], g[key]), key
    if case == "VAN_2L":
        for name in ("flux", "anti_diag", "flux_anti_diag", "anti_diag_abs_rates",
                     "anti_diag_fiss_rates", "density_00_3"):
            assert (tmp_path / f"{name}.png").stat().st_size > 0
