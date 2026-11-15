# Tutorial: DRAGON5 / Serpent2 code-to-code post-treatment with starterDD
#
# Shows the PostTreatment API (readers, detector manifests, pin bookkeeping,
# comparisons and figures) on the GE14 DOM assembly at 0 % void.
#
#   Part A  runs anywhere: it only uses data committed in the repository
#           (tests/data/ge14_mini: a real subset of a DRAGON COMPO and of an
#           existing Serpent2 reference, 4 pins out of 51).
#   Part B  runs when the full local GE14 results are available
#           ($STARTERDD_TEST_DATA, default <repo>/temp_examples): normalised
#           full-assembly maps, axis profiles, fine-group spectra and the
#           one-call legacy report.
#
# Every file is written to tutorials/results/post_treatment/ (git-ignored).
# Requires the post-treatment extra:  pip install -e ".[post]"
#
# Run:  python tutorials/post_treatment_code_to_code.py
#
# Date : 14/11/2026

import os
import warnings
from pathlib import Path

import numpy as np

from starterDD.PostTreatment import (
    CompoReader, DetectorManifest, PinLattice, SerpentRun,
    compare_keff, compare_spectra, compute_pin_rates, physics_data,
)
from starterDD.PostTreatment import plotting
from starterDD.PostTreatment.reports import legacy_report_for_run

# =====================================================================
# Paths
# =====================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
GE14_INPUTS = PROJECT_ROOT / "data" / "BWRProgressionProblems" / "GE14_inputs" / "ASSEMBLY"
MINI = PROJECT_ROOT / "tests" / "data" / "ge14_mini"
FULL_DATA = Path(os.environ.get("STARTERDD_TEST_DATA", PROJECT_ROOT / "temp_examples"))

OUT = PROJECT_ROOT / "tutorials" / "results" / "post_treatment"
OUT_A = OUT / "A_committed_subset"
OUT_B = OUT / "B_full_assembly"
OUT_A.mkdir(parents=True, exist_ok=True)

LIBRARY = "endfb8r1"                       # Q-values used by the power normalisation
ISOTOPES = ["U235", "U238", "Gd155", "Gd157"]


def title(text):
    print("\n" + f" {text} ".center(72, "="))


# =====================================================================
# PART A - committed data only
# =====================================================================

# ---------------------------------------------------------------------
# A.1  DRAGON COMPO: outputs, parameters, statepoints, keff
# ---------------------------------------------------------------------
title("A.1 DRAGON COMPO")
# CompoReader accepts the ASCII export (here gzipped). The first call parses
# it and writes an HDF5 cache; later calls read only what they need from the
# cache. The cache is put in the output folder (default: next to the COMPO).
compo = CompoReader(MINI / "_CPO_GE14_DOM_DIAG_mini.gz",
                    cache=OUT_A / "_CPO_GE14_DOM_DIAG_mini.h5",
                    labels={"00": {"DCool": 0.73669}})   # optional statepoint labels
print("outputs      :", compo.outputs())
print("parameters   :", compo.parameters("EDIR_2G"))
for sp in compo.statepoints("EDIR_2G"):
    print("statepoint   :", sp.calc, sp.values, "label:", sp.label)
print("keff (00)    :", compo.keff("EDIHOM_COND", sp="00"))
print("energy bounds:", compo.energy_bounds("EDIR_2G"), "eV (ascending)")
# Multi-statepoint COMPOs: select with sp=<calc index>, sp={"DCool": 0.457},
# or a label. Group-wise accessors return (n_mixtures, n_groups) arrays:
flux_2g = compo.flux("EDIR_2G", sp="00")
print("NWT0 flux    :", flux_2g.shape, "(pins x groups)")

# ---------------------------------------------------------------------
# A.2  An existing Serpent2 reference with non-standard detector names
# ---------------------------------------------------------------------
title("A.2 Serpent2 reference and its detector manifest")
serp_input = MINI / "GE14_DOM_00_endfb8r1.serp.gz"
# The manifest says what each detector measures (observable, pin, ordered
# (MT, nuclide) responses, energy grid). For an existing run it is inferred
# from the content of the det cards, never from the detector names:
inferred = DetectorManifest.infer(serp_input)
print(inferred.summary())
print("pin 15 detector:", inferred.find("PIN_295G", "rr", pin=15).detector)
# Freeze it as JSON (normally next to the run as <input>.detectors.json, where
# SerpentRun picks it up automatically; here in the output folder).
inferred.write(OUT_A / "GE14_DOM_00_endfb8r1.serp.detectors.json")

# This subset only contains the det0.m blocks of 4 pins, so we use the frozen
# manifest committed with it (SerpentRun finds it next to the input).
serp = SerpentRun(serp_input)
print("detectors in det0.m:", len(serp.detectors), "| consistency problems:", serp.check())
pins_available = serp.manifest.pins("PIN_295G", "rr")
print("pins with detectors:", pins_available)

k = compare_keff(compo, "EDIHOM_COND", serp)
print(f"keff DRAGON {k.dragon:.6f}  Serpent {k.serpent:.6f} +- {k.sigma_pcm:.1f} pcm"
      f"  ->  {k.delta_pcm:+.1f} pcm")

# ---------------------------------------------------------------------
# A.3  Pin bookkeeping and a full-lattice DRAGON map
# ---------------------------------------------------------------------
title("A.3 Pin lattice")
# Built from the same geometry / composition YAMLs as the DRAGON case.
lattice = PinLattice.from_inputs(GE14_INPUTS / "GEOM_GE14_DOM_DIAG.yaml",
                                 GE14_INPUTS / "material_compos_00.yaml",
                                 name="GE14_DOM_DIAG")
print(f"{lattice.n_pins} unique pins, {lattice.n_physical_pins} physical pins, "
      f"symmetry: {lattice.symmetry}")
print("pins on the symmetry axis:", lattice.axis_pins)
print("pin 2: multiplicity", lattice.multiplicity(2), "| pin 1: DRAGON fraction",
      lattice.dragon_fraction(1), "(cut in half by the diagonal)")

# DRAGON alone covers every pin: map the thermal / fast flux of each pin.
# NWT0 is a volume-averaged flux, so no symmetry correction is needed here.
fig = plotting.new_figure(figsize=(11, 5))
axs = fig.subplots(1, 2)
for g, name in enumerate(plotting.group_names(2)):
    rel = flux_2g[:, g] / np.average(flux_2g[:, g], weights=lattice.multiplicities())
    plotting.draw_pin_map(axs[g], lattice, rel, cmap="viridis",
                          title=f"DRAGON5 {name} flux / assembly mean")
fig.suptitle("GE14 DOM, 0 % void - pin-wise flux (EDIR_2G)", fontweight="bold")
fig.savefig(OUT_A / "dragon_flux_map.png", dpi=150)

# ---------------------------------------------------------------------
# A.4  Code-to-code on the 4 committed pins
# ---------------------------------------------------------------------
title("A.4 DRAGON5 vs Serpent2 pin rates (4 pins)")
# Rates are per physical pin in both codes (Serpent detector / multiplicity,
# DRAGON / dragon_fraction). The Serpent 295g tallies are condensed onto the
# DRAGON 2g mesh, errors in quadrature.
cmp4 = compute_pin_rates(compo, "EDIR_2G", serp, "PIN_295G", ISOTOPES, lattice,
                         LIBRARY, sp="00", pins=pins_available)
for iso, cov in cmp4.coverage.items():
    print(f"{iso:6s} compared MTs {cov.compared}  fission={cov.fission}  "
          f"DRAGON-only {list(cov.dragon_only)}")
# A global normalisation needs every pin (cmp4.errors() would raise
# NormalisationError). The fission / absorption ratio needs none:
f_over_a = cmp4.errors("fission_over_absorption")
rel_sigma = cmp4.serpent.absorption_err / cmp4.serpent.absorption
for p, e, s in zip(cmp4.pins, f_over_a, rel_sigma):
    print(f"pin {p:2d} ({lattice.material[p]:9s}) fission/absorption D5-S2: "
          f"{e[0]:+.2f} % (thermal) {e[1]:+.2f} % (fast) | Serpent 1-sigma "
          f"{100 * s[0]:.4f} %")

fig = plotting.new_figure(figsize=(11, 5))
axs = fig.subplots(1, 2)
for g, name in enumerate(plotting.group_names(2)):
    plotting.draw_pin_map(axs[g], lattice, f_over_a[:, g], pins=cmp4.pins,
                          symmetric=True, title=f"fission/absorption, {name}: D5-S2 (%)")
fig.suptitle("4 committed pins (other pins not in the subset)", fontweight="bold")
fig.savefig(OUT_A / "fission_over_absorption_4pins.png", dpi=150)

# ---------------------------------------------------------------------
# A.5  Spectra of the 4 pins
# ---------------------------------------------------------------------
spec4 = compare_spectra(compo, "EDIR_2G", serp, "PIN_295G", sp="00",
                        pins=pins_available, lattice=lattice)
plotting.plot_spectrum(spec4, OUT_A / "pin_spectra_2g.png")
print("\nPart A figures written to", OUT_A)


# =====================================================================
# PART B - full local results (skipped when absent)
# =====================================================================
RUN = (FULL_DATA / "GE14_DOM" / "results" / "coolant_density_1L_PT_3x3_HFP"
       / "GE14_DOM_DIAG" / "PT_IC_1L_MOC" / "2026-07-19T10-46-27")
SERP_FULL = FULL_DATA / "GE14_DOM" / "Serpent" / "GE14_DOM_00_endfb8r1.serp"

if not ((RUN / "_CPO_GE14_DOM_DIAG").exists() and SERP_FULL.exists()):
    title("PART B skipped")
    print(f"Full GE14 results not found under {FULL_DATA}.")
    print("Set STARTERDD_TEST_DATA to the folder holding GE14_DOM/ to run it.")
else:
    OUT_B.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    # B.1  Full-assembly pin-rate comparison
    # -----------------------------------------------------------------
    title("B.1 Full assembly (51 pins)")
    compo_f = CompoReader(RUN / "_CPO_GE14_DOM_DIAG", cache=OUT_B / "_CPO_GE14_DOM_DIAG.h5")
    serp_f = SerpentRun(SERP_FULL, manifest=DetectorManifest.infer(SERP_FULL))
    cmp = compute_pin_rates(compo_f, "EDIR_2G", serp_f, "PIN_295G", ISOTOPES, lattice, LIBRARY)
    print("power normalisation factor (DRAGON -> Serpent):", cmp.normalisation_factor())
    for q in ("absorption", "fission"):
        for g, s in zip(plotting.group_names(2), cmp.stats(q)):
            print(f"{q:10s} {g:7s} RMS {s['rms']:.3f} %  max |err| "
                  f"{max(abs(s['min']), abs(s['max'])):.3f} %  weighted mean {s['weighted_mean']:+.3f} %")
    plotting.plot_rate_maps(cmp, OUT_B / "rate_maps_2g.png")
    plotting.plot_profile(cmp, lattice.diagonal_pins("axis"), "absorption",
                          path=OUT_B / "axis_profile_absorption.png")
    plotting.plot_profile(cmp, lattice.diagonal_pins("axis"), "fission",
                          path=OUT_B / "axis_profile_fission.pdf")   # any format

    # -----------------------------------------------------------------
    # B.2  Fine-group comparisons: the 295-group meshes differ
    # -----------------------------------------------------------------
    title("B.2 Fine groups")
    # Energy meshes: DRAGON is the code of interest, so comparisons are made
    # on the DRAGON output mesh by default (mesh="dragon"), and the Serpent
    # grid must refine it. The exporter writes the DRAGON SHEM295 bounds
    # exactly (physics_data.G295_EV), so new references compare group by group.
    #
    # This reference predates the 2026-11-14 grid fix: its 295g grid has four
    # wrong bounds (indices 235-238, 3.9-7.5 keV), so DRAGON groups 57-61
    # (3481-9119 eV) do not line up with the Serpent bins. The strict default
    # refuses; mesh="common" merges these 5 groups into one bin (exact, the
    # outer bounds are shared) and compares on 291 groups. The merged range is
    # recorded in .merged_eV and shaded grey in the figures. Rates are plotted
    # per unit lethargy, so the wide merged bin is not drawn as a bump.
    print("Serpent grid is the legacy one:",
          physics_data.is_legacy_mesh(serp_f.tally("ALL_295G", "flux").bounds_eV))
    try:
        compare_spectra(compo_f, "EDIHOM_295", serp_f, "ALL_295G")
    except physics_data.MeshMismatchError as exc:
        print("strict mesh check:", str(exc)[-150:])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        spec = compare_spectra(compo_f, "EDIHOM_295", serp_f, "ALL_295G", mesh="common")
        fine = compute_pin_rates(compo_f, "EDIR_295", serp_f, "PIN_295G", ["U238"],
                                 lattice, LIBRARY, mesh="common")
    print(f"homogenised spectrum on {len(spec.bounds_eV) - 1} common groups, "
          f"merged {[(round(a), round(b)) for a, b in spec.merged_eV]} eV, "
          f"max |D5-S2| {np.nanmax(np.abs(spec.errors())):.1f} %")
    plotting.plot_spectrum(spec, OUT_B / "assembly_spectrum.png", per_lethargy=True)
    plotting.plot_pin_spectra(fine, [1, 28], "absorption",
                              path=OUT_B / "U238_absorption_pins_1_28.png")

    # -----------------------------------------------------------------
    # B.3  One call: the outputs of the original analysis scripts
    # -----------------------------------------------------------------
    title("B.3 Legacy report")
    report = legacy_report_for_run(RUN, SERP_FULL, case_name="density_00",
                                   out_dir=OUT_B / "legacy_report",
                                   csv_path=OUT_B / "legacy_report" / "summary.csv",
                                   cache=OUT_B / "_CPO_GE14_DOM_DIAG.h5")
    print(f"{len(report.files)} files, keff {report.keff.delta_pcm:+.2f} pcm")
    compo_f.close()
    print("\nPart B figures written to", OUT_B)

compo.close()
