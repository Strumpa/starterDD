# -*- coding: utf-8 -*-
"""Drop-in replacement for the original GE14 analysis driver.

:func:`write_legacy_report` produces, from the new comparison core, the files
that the original ``results_analysis_GE14_functions`` calls wrote into a run
directory:

======================================  ==========================================
``<case>_<k>.txt``                      statistics per isotope set (same format)
``<case>_<k>.png``                      error / distribution maps per isotope set
``<n>G_data.npz``                       pin-wise arrays of the last isotope set
summary CSV row(s)                      ``wall;k_D5;k_S2;pcm;<statistics>``
``flux.png``                            homogenised spectrum
``anti_diag.png``                       coarse rates along the symmetry axis
``flux_anti_diag.png``                  pin spectra along the axis
``anti_diag_{abs,fiss}_rates.png``      fine-group pin rates along the axis
======================================  ==========================================

With ``legacy_numbers=True`` (default) the txt / csv / npz content is the same
as the original scripts (benchmark densities, legacy error combination);
``False`` uses the physical defaults of :func:`compute_pin_rates` (Serpent
``mat`` densities, quadrature errors), which change the % values by <1e-4
points. Figures are new (same content, different layout).

:func:`legacy_report_for_run` does the whole job for one DRAGON run directory
(the ``dragon_runner`` layout with ``run_manifest.yaml`` and ``inputs/``),
choosing the COMPO outputs by content rather than by name.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from . import physics_data as phys
from .comparison import (
    KeffComparison, PinRateComparison, compare_keff, compare_spectra,
    compute_pin_rates,
)
from .compo_reader import CompoReader, StatepointSel
from .lattice_maps import PinLattice
from .serpent_reader import DetectorManifest, PinResolver, SerpentRun

PathLike = Union[str, Path]

#: Isotope sets of the original GE14 analysis (one txt/png/csv row each).
DEFAULT_ISOTOPE_SETS: Tuple[Tuple[str, ...], ...] = (
    ("U235",), ("U238",), ("Gd155", "Gd157"), ("U235", "U238", "Gd155", "Gd157"))


@dataclass
class LegacyReport:
    """What :func:`write_legacy_report` computed and wrote."""

    out_dir: Path
    keff: KeffComparison
    rates: List[PinRateComparison]
    fine_rates: Optional[PinRateComparison] = None
    files: List[Path] = field(default_factory=list)
    csv_rows: List[str] = field(default_factory=list)


def legacy_csv_row(keff: KeffComparison, cmp: PinRateComparison,
                   wall_time: Optional[float] = None, method: str = "power") -> str:
    """One summary-CSV line in the original format (no trailing newline)."""
    head = [("" if wall_time is None else wall_time), keff.dragon, keff.serpent,
            keff.delta_pcm]
    return ";".join(str(x) for x in head + cmp.legacy_csv_fields(method))


def write_legacy_report(compo: CompoReader, serpent: SerpentRun, lattice: PinLattice,
                        out_dir: PathLike, *, library: str,
                        case_name: str = "comparison",
                        isotope_sets: Sequence[Sequence[str]] = DEFAULT_ISOTOPE_SETS,
                        sp: StatepointSel = None,
                        keff_output: str = "EDIHOM_COND",
                        rate_output: str = "EDIR_2G",
                        spectrum_output: Optional[str] = None,
                        fine_rate_output: Optional[str] = None,
                        rate_observable: str = "PIN_295G",
                        spectrum_observable: str = "ALL_295G",
                        profile_pins: Optional[Sequence[int]] = None,
                        csv_path: Optional[PathLike] = None,
                        wall_time: Optional[float] = None,
                        legacy_numbers: bool = True,
                        fine_mesh: str = "common",
                        figures: bool = True,
                        figure_format: str = "png") -> LegacyReport:
    """Write the original analysis outputs for one statepoint.

    Args:
        compo, serpent, lattice: Inputs of :func:`compute_pin_rates`.
        out_dir: Directory for the files (created if needed).
        library: Q-value library (e.g. ``"endfb8r1"``).
        case_name: Prefix of the per-isotope-set files (the original driver
            used ``density_<void label>``).
        isotope_sets: One statistics file / figure / CSV row per set; the
            npz is written for the last set.
        keff_output: Output holding the keff to compare.
        rate_output: ``by_pin`` output of the pin-rate maps (coarse mesh).
        spectrum_output: Homogenised output of ``flux.png`` (skipped if None).
        fine_rate_output: Fine-mesh ``by_pin`` output of the axis spectra and
            fine-group rates (skipped if None).
        rate_observable, spectrum_observable: Serpent observables (rates and
            per-pin flux use *rate_observable*).
        profile_pins: Pins of the profile figures (default: the symmetry
            axis, as the original ``anti_diag`` lists).
        csv_path: Summary CSV to append rows to (skipped if None).
        wall_time: First CSV column (DRAGON wall time).
        legacy_numbers: Reproduce the original numbers exactly (see module
            docstring).
        fine_mesh: Mesh option of the fine-group comparisons. ``"common"``
            (default) also handles references made on the legacy 295g grid
            (:data:`physics_data.G295_STARTERDD_LEGACY_EV`); on a reference
            generated after the 2026-11-14 grid fix it equals ``"dragon"``.
        figures: Write figures (needs matplotlib).
        figure_format: Extension of the figures.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rate_kw = dict(library=library, sp=sp)
    if legacy_numbers:
        rate_kw.update(serpent_densities=lattice.composition, error_mode="legacy")
    if profile_pins is None:
        profile_pins = lattice.diagonal_pins("axis") if lattice.symmetry else []
    plotting = None
    if figures:
        from . import plotting as _plotting   # matplotlib needed only here
        plotting = _plotting
    ext = figure_format.lstrip(".")

    keff = compare_keff(compo, keff_output, serpent, sp)
    report = LegacyReport(out, keff, [])

    for k, isos in enumerate(isotope_sets):
        cmp = compute_pin_rates(compo, rate_output, serpent, rate_observable,
                                list(isos), lattice, **rate_kw)
        report.rates.append(cmp)
        txt = out / f"{case_name}_{k}.txt"
        txt.write_text(cmp.legacy_stats_text())
        report.files.append(txt)
        report.csv_rows.append(legacy_csv_row(keff, cmp, wall_time))
        if plotting:
            png = out / f"{case_name}_{k}.{ext}"
            plotting.plot_rate_maps(cmp, png)
            report.files.append(png)
    last = report.rates[-1]
    npz = out / f"{last.n_groups}G_data.npz"
    last.to_legacy_npz(npz)
    report.files.append(npz)

    if csv_path is not None:
        with open(csv_path, "a") as f:
            for row in report.csv_rows:
                f.write(row + "\n")

    if plotting:
        def fig(fn, name, *a, **kw):
            path = out / f"{name}.{ext}"
            fn(*a, path=path, **kw)
            report.files.append(path)

        if spectrum_output:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                spec = compare_spectra(compo, spectrum_output, serpent, spectrum_observable,
                                       sp, mesh=fine_mesh)
            fig(plotting.plot_spectrum, "flux", spec, per_lethargy=True)
        if profile_pins:
            fig(plotting.plot_profile, "anti_diag", last, profile_pins)
    if fine_rate_output:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)   # merged groups, see mesh docs
            fine = compute_pin_rates(compo, fine_rate_output, serpent, rate_observable,
                                     list(isotope_sets[-1]), lattice, mesh=fine_mesh,
                                     **rate_kw)
            report.fine_rates = fine
            if plotting and profile_pins:
                spec = compare_spectra(compo, fine_rate_output, serpent, rate_observable,
                                       sp, pins=profile_pins, lattice=lattice, mesh=fine_mesh)
                fig(plotting.plot_spectrum, "flux_anti_diag", spec)
                fig(plotting.plot_pin_spectra, "anti_diag_abs_rates", fine, profile_pins,
                    "absorption")
                if fine.has_fission:
                    fig(plotting.plot_pin_spectra, "anti_diag_fiss_rates", fine,
                        profile_pins, "fission")
    return report


def _choose_outputs(compo: CompoReader, n_pins: int, sp: StatepointSel):
    """(keff, spectrum, rate, fine_rate) outputs, by number of mixtures and
    groups: homogenised outputs have 1 mixture, by_pin outputs one per pin."""
    homog, by_pin = [], []
    for o in compo.outputs():
        n_mix = compo.n_mixtures(o)
        ng = compo.n_groups(o, sp)
        if n_mix == 1:
            homog.append((ng, o))
        elif n_mix == n_pins:
            by_pin.append((ng, o))
    if not homog or not by_pin:
        raise ValueError(f"Cannot identify homogenised and by_pin outputs among "
                         f"{compo.outputs()}.")
    homog.sort()
    by_pin.sort()
    spectrum = homog[-1][1] if len(homog) > 1 else None
    fine = by_pin[-1][1] if len(by_pin) > 1 else None
    return homog[0][1], spectrum, by_pin[0][1], fine


def legacy_report_for_run(run_dir: PathLike, serpent_input: PathLike, *,
                          case_name: str = "comparison",
                          library: Optional[str] = None,
                          sp: StatepointSel = None,
                          csv_path: Optional[PathLike] = None,
                          out_dir: Optional[PathLike] = None,
                          cache: Union[bool, PathLike] = True,
                          serpent_manifest: Union[None, PathLike, DetectorManifest] = None,
                          pin_resolver: Optional[PinResolver] = None,
                          outputs: Optional[Mapping[str, Optional[str]]] = None,
                          **kwargs) -> LegacyReport:
    """Write the original analysis outputs for a ``dragon_runner`` run dir.

    Reads ``run_manifest.yaml`` (case name, inputs, library, wall time),
    builds the :class:`PinLattice` from the archived ``inputs/`` YAMLs, opens
    the COMPO (with its HDF5 cache unless ``cache=False``) and the Serpent
    run, then calls :func:`write_legacy_report` into *out_dir* (default: the
    run directory, like the original driver).

    Args:
        outputs: Override the automatic choice, keys ``keff_output``,
            ``spectrum_output``, ``rate_output``, ``fine_rate_output``.
        kwargs: Passed to :func:`write_legacy_report` (``isotope_sets``,
            ``legacy_numbers``, ``figures``...).
    """
    import yaml

    run = Path(run_dir)
    man = yaml.safe_load((run / "run_manifest.yaml").read_text())
    inputs = run / "inputs"
    cfg = man.get("config_yamls", {})
    geom = inputs / Path(cfg.get("GEOM", "GEOM.yaml")).name
    mats = inputs / Path(cfg.get("MATS", "MATS.yaml")).name
    lattice = PinLattice.from_inputs(geom, mats, name=man.get("case_name", "assembly"))
    if library is None:
        aliases = list((man.get("draglibs") or {}).values())
        if not aliases:
            raise ValueError("library not given and not found in run_manifest.yaml")
        library = phys.canonical_library(aliases[0])
    compo = CompoReader(run / f"_CPO_{man['case_name']}", cache=cache)
    serpent = SerpentRun(serpent_input, manifest=serpent_manifest, pin_resolver=pin_resolver)
    keff_o, spec_o, rate_o, fine_o = _choose_outputs(compo, lattice.n_pins, sp)
    chosen = dict(keff_output=keff_o, spectrum_output=spec_o, rate_output=rate_o,
                  fine_rate_output=fine_o)
    chosen.update(outputs or {})
    wall = (man.get("execution") or {}).get("wall_time_seconds")
    try:
        return write_legacy_report(compo, serpent, lattice, out_dir or run,
                                   library=library, case_name=case_name, sp=sp,
                                   csv_path=csv_path, wall_time=wall, **chosen, **kwargs)
    finally:
        compo.close()


__all__ = ["DEFAULT_ISOTOPE_SETS", "LegacyReport", "legacy_csv_row",
           "write_legacy_report", "legacy_report_for_run"]
