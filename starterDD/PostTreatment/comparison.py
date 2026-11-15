# -*- coding: utf-8 -*-
"""DRAGON5 / Serpent2 code-to-code comparison core.

Pairing goes through the COMPO statepoint API (:class:`CompoReader`) and the
Serpent detector manifest (:class:`SerpentRun`), never through names.

Conventions
-----------
**Per physical pin.** Every pin-wise quantity is a volume-integrated rate for
*one physical pin*:

- Serpent2: a detector scoring ``pin_idx`` covers every lattice position of
  that pin (the universe is reused), so ``rate = N * det / multiplicity``. The
  multiplicity is read from the Serpent input (``pin``/``lat`` cards) and
  checked against the model. Detectors are volume-integrated (the per-pin
  detectors add up to the assembly detector), so no volume factor is applied.
- DRAGON5: ``rate = sigma * phi * N * V / dragon_fraction`` with ``phi`` the
  ``NWT0`` flux, ``N``/``V`` the COMPO ``ISOTOPESDENS``/``ISOTOPESVOL``, and
  ``dragon_fraction`` = 1/2 for pins cut by the diagonal symmetry axis.

**Normalisation** (method ``"power"``, the default): both codes are scaled to
the same full-assembly total ``sum_p multiplicity(p) * power(p)``, where
``power`` is the Q-weighted sum of the compared reactions of the studied
isotopes (the ``power_subset`` of the plan; the legacy analysis does the same).

These conventions give exactly the same relative differences as the original
GE14 analysis, which carried a uniform hidden factor ``2*pi*r^2`` on every
Serpent value (see :meth:`PinRateComparison.to_legacy_npz`).

**Errors.** Serpent statistical errors are absolute 1-sigma values. They are
condensed in energy and combined over reactions/isotopes according to
``error_mode``: ``"quadrature"`` (default, independent tallies),
``"linear"`` (sum of |contributions|, conservative) or ``"legacy"`` (linear
with the signed absorption factors, as in the original scripts; only for
regression).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from . import physics_data as phys
from .compo_reader import CompoReader, StatepointSel
from .lattice_maps import PinLattice
from .serpent_reader import ManifestEntry, SerpentRun

QUANTITIES = ("absorption", "fission", "fission_over_absorption")


# ═════════════════════════════════════════════════════════════════════════════
# keff
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class KeffComparison:
    dragon: float
    serpent: float
    serpent_rel_err: float
    estimator: str

    @property
    def delta_pcm(self) -> float:
        """DRAGON - Serpent, in pcm (1e-5)."""
        return (self.dragon - self.serpent) * 1e5

    @property
    def sigma_pcm(self) -> float:
        """Serpent 1-sigma, in pcm."""
        return self.serpent * self.serpent_rel_err * 1e5


def compare_keff(compo: CompoReader, output: str, serpent: SerpentRun,
                 sp: StatepointSel = None, estimator: str = "IMP_KEFF") -> KeffComparison:
    k_s, rel = serpent.keff(estimator)
    return KeffComparison(compo.keff(output, sp), k_s, rel, estimator)


# ═════════════════════════════════════════════════════════════════════════════
# Pin-wise reaction rates
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class CodeRates:
    """Pin-wise rates of one code (per physical pin, ascending energy).

    Arrays are shaped ``(n_pins, n_groups)``, rows in the order of *pins*.
    """

    code: str
    pins: List[int]
    bounds_eV: np.ndarray
    fission: np.ndarray
    absorption: np.ndarray
    power: np.ndarray
    multiplicity: np.ndarray
    fission_err: Optional[np.ndarray] = None
    absorption_err: Optional[np.ndarray] = None

    def get(self, quantity: str) -> np.ndarray:
        if quantity == "fission_over_absorption":
            with np.errstate(divide="ignore", invalid="ignore"):
                return self.fission / self.absorption
        return getattr(self, quantity)

    def total(self, quantity: str = "power") -> np.ndarray:
        """Full-assembly total per group: ``sum_p multiplicity(p) * q(p)``."""
        return (self.multiplicity[:, None] * self.get(quantity)).sum(axis=0)


@dataclass
class ReactionCoverage:
    """Which reactions of one isotope were compared, and the share of the
    total absorption of those scored by one code only (max over pins and
    groups)."""

    isotope: str
    reference_pin: Optional[int]
    fission: bool = False
    compared: List[int] = field(default_factory=list)
    dragon_only: Dict[int, float] = field(default_factory=dict)
    serpent_only: Dict[int, float] = field(default_factory=dict)


class NormalisationError(ValueError):
    """The rates do not cover the whole assembly, so no global
    normalisation can be computed."""


@dataclass
class PinRateComparison:
    """Pin-wise fission/absorption rates of both codes for a set of isotopes."""

    dragon: CodeRates
    serpent: CodeRates
    lattice: PinLattice
    isotopes: List[str]
    coverage: Dict[str, ReactionCoverage]
    output: str = ""
    observable: str = ""
    error_mode: str = "quadrature"
    #: ``(E_lo, E_hi)`` eV ranges where DRAGON groups were merged
    #: (``mesh="common"``); empty on the DRAGON mesh.
    merged_eV: List[Tuple[float, float]] = field(default_factory=list)

    @property
    def pins(self) -> List[int]:
        return self.dragon.pins

    @property
    def bounds_eV(self) -> np.ndarray:
        return self.dragon.bounds_eV

    @property
    def n_groups(self) -> int:
        return len(self.bounds_eV) - 1

    @property
    def has_fission(self) -> bool:
        return any(c.fission for c in self.coverage.values())

    @property
    def is_complete(self) -> bool:
        return sorted(self.pins) == sorted(self.lattice.ordered_pins)

    # ── normalisation ──────────────────────────────────────────────────────
    def normalisation_factor(self, method: str = "power") -> float:
        """Factor applied to DRAGON rates so that both codes share the same
        full-assembly total."""
        if not self.is_complete:
            raise NormalisationError(
                f"Rates cover {len(self.pins)} of {self.lattice.n_pins} pins; "
                "a global normalisation needs the whole assembly.")
        if method == "power":
            return float(self.serpent.total("power").sum() / self.dragon.total("power").sum())
        if method == "absorption":
            return float(self.serpent.total("absorption").sum()
                         / self.dragon.total("absorption").sum())
        if method == "fission":
            return float(self.serpent.total("fission").sum()
                         / self.dragon.total("fission").sum())
        raise ValueError(f"Unknown normalisation method '{method}'.")

    def errors(self, quantity: str = "absorption", method: str = "power") -> np.ndarray:
        """Relative difference DRAGON vs Serpent, in %, ``(n_pins, n_groups)``.

        ``fission_over_absorption`` is normalisation-free. NaN where Serpent
        scores zero.
        """
        s = self.serpent.get(quantity)
        if quantity == "fission_over_absorption":
            d = self.dragon.get(quantity)
        else:
            d = self.dragon.get(quantity) * self.normalisation_factor(method)
        with np.errstate(divide="ignore", invalid="ignore"):
            e = (d - s) / s * 100.0
        e[~np.isfinite(e)] = np.nan
        return e

    def distribution(self, quantity: str = "absorption", code: str = "serpent") -> np.ndarray:
        """Pin-wise relative rate, mean 1 over the physical pins where the
        rate is non-zero (NaN elsewhere)."""
        rates = (self.serpent if code == "serpent" else self.dragon)
        q = rates.get(quantity).astype(float)
        m = rates.multiplicity[:, None]
        nonzero = q != 0
        total = (m * q).sum(axis=0)
        n_phys = (m * nonzero).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = q / total * n_phys
        out[~nonzero] = np.nan
        return out

    # ── statistics ─────────────────────────────────────────────────────────
    def stats(self, quantity: str = "absorption", method: str = "power") -> List[Dict[str, float]]:
        """Per-group statistics over **physical** pins (each unique pin
        weighted by its multiplicity). Same definitions as the original
        analysis: min |err|, min, max, mean, mean weighted by the Serpent
        distribution, RMS and weighted RMS (no weighted values for
        ``fission_over_absorption``)."""
        err = self.errors(quantity, method)
        m = self.serpent.multiplicity.astype(float)
        dist = None if quantity == "fission_over_absorption" else \
            self.distribution(quantity, "serpent")
        out = []
        for g in range(self.n_groups):
            e = err[:, g]
            ok = np.isfinite(e) if dist is None else np.isfinite(e) & np.isfinite(dist[:, g])
            w = m[ok]
            eg = e[ok]
            row = {
                "min_abs": float(np.min(np.abs(eg))),
                "min": float(np.min(eg)),
                "max": float(np.max(eg)),
                "mean": float(np.sum(w * eg) / w.sum()),
                "rms": float(np.sqrt(np.sum(w * eg ** 2) / w.sum())),
            }
            if dist is not None:
                d = dist[ok, g]
                row["weighted_mean"] = float(np.sum(w * eg * d) / w.sum())
                row["weighted_rms"] = float(np.sqrt(np.sum(w * eg ** 2 * d) / w.sum()))
            out.append(row)
        return out

    # ── legacy outputs ─────────────────────────────────────────────────────
    def _legacy_arrays(self, fuel_radius: float) -> Dict[str, np.ndarray]:
        n = max(self.lattice.ordered_pins)
        idx = np.array(self.pins) - 1
        k = 2.0 * np.pi * fuel_radius ** 2   # hidden legacy factor on Serpent

        def unfold(a, scale=1.0):
            out = np.zeros((self.n_groups, n))
            out[:, idx] = (a * scale).T
            return out

        s, d = self.serpent, self.dragon
        return {
            "serp_groups_fiss_rate_tot": unfold(s.fission, k),
            "serp_groups_fiss_rate_err_tot": unfold(s.fission_err, k),
            "serp_groups_abs_rate_tot": unfold(s.absorption, k),
            "serp_groups_abs_rate_err_tot": unfold(s.absorption_err, k),
            "drag_groups_fiss_rate_tot": unfold(d.fission),
            "drag_groups_abs_rate_tot": unfold(d.absorption),
        }

    def to_legacy_npz(self, path: Union[str, Path, None] = None,
                      fuel_radius: Optional[float] = None) -> Dict[str, np.ndarray]:
        """Arrays of the original ``2G_data.npz`` (written to *path* if
        given). Requires ``error_mode="legacy"`` for the error arrays to be
        identical, and the whole assembly for ``serp_pw``/``drag_pw``."""
        r = fuel_radius if fuel_radius is not None else self.lattice.fuel_radius
        if r is None:
            raise ValueError("fuel_radius is required for legacy outputs.")
        data = self._legacy_arrays(r)
        if self.is_complete:
            data["serp_pw"] = np.float64(2.0 * np.pi * r ** 2 * self.serpent.total("power").sum())
            data["drag_pw"] = np.float64(self.dragon.total("power").sum())
        data.update(
            n_energy_groups=np.int64(self.n_groups),
            symmetry=np.str_(self.lattice.symmetry),
            symmetry_factors=self.lattice.legacy_symmetry_factors(),
            is_fission_tot=np.bool_(self.has_fission),
        )
        if path is not None:
            np.savez(path, **data)
        return data

    def legacy_stats_text(self, method: str = "power") -> str:
        """Text of the original ``<name>_<k>.txt`` statistics file."""
        lines = []
        blocks = [("absorption", "absorption rate", True)]
        if self.has_fission:
            blocks += [("fission", "fission rate", True),
                       ("fission_over_absorption", "fission over absorption rate", False)]
        for q, label, weighted in blocks:
            for g, s in enumerate(self.stats(q, method), start=1):
                lines.append(f"For the {label} in the energy group {g}:")
                lines.append(f"    - Minimal absolute error: {s['min_abs']} %")
                lines.append(f"    - Maximum errors (bounds): ({s['min']} %, {s['max']} %)")
                lines.append(f"    - Average error: {s['mean']} %")
                if weighted:
                    lines.append(f"    - Weighted average error: {s['weighted_mean']} %")
                lines.append(f"    - RMS error: {s['rms']} %")
                if weighted:
                    lines.append(f"    - Weighted RMS error: {s['weighted_rms']} %")
                lines.append("")
        return "\n".join(lines) + "\n"

    def legacy_csv_fields(self, method: str = "power") -> List[float]:
        """Statistics columns appended by the original analysis to the
        summary CSV (after wall time, keff values and pcm)."""
        keys_w = ["min_abs", "min", "max", "mean", "weighted_mean", "rms", "weighted_rms"]
        keys = ["min_abs", "min", "max", "mean", "rms"]
        out: List[float] = []
        for s in self.stats("absorption", method):
            out += [s[k] for k in keys_w]
        if self.has_fission:
            for s in self.stats("fission", method):
                out += [s[k] for k in keys_w]
            for s in self.stats("fission_over_absorption", method):
                out += [s[k] for k in keys]
        return out


def _target_mesh(dragon_bounds: np.ndarray, serpent_bounds: Optional[Sequence[float]],
                 mesh: str, what: str) -> Tuple[np.ndarray, List[Tuple[float, float]]]:
    """Comparison mesh and merged energy ranges.

    ``"dragon"`` (default): the DRAGON output mesh, DRAGON being the code of
    interest; the Serpent grid must refine it. ``"common"``: the bounds common
    to both codes, for references whose grid does not refine DRAGON's (e.g.
    the pre-2026-11-14 ``295g`` grid, :data:`physics_data.G295_STARTERDD_LEGACY_EV`).

    Returns ``(bounds, merged)`` where *merged* lists the ``(E_lo, E_hi)`` eV
    ranges in which several DRAGON groups were merged (empty for ``"dragon"``).
    """
    if mesh not in ("dragon", "common"):
        raise ValueError(f"Unknown mesh option '{mesh}' (use 'dragon' or 'common').")
    if mesh == "dragon" or serpent_bounds is None:
        if serpent_bounds is not None:
            try:
                phys.condensation_indices(serpent_bounds, dragon_bounds)
            except phys.MeshMismatchError as exc:
                hint = (" The Serpent reference uses the legacy starterDD 295g grid"
                        " (before the 2026-11-14 fix); compare with mesh='common',"
                        " or re-run Serpent with the current exporter."
                        if phys.is_legacy_mesh(serpent_bounds) else
                        " Use mesh='common' to compare on the bounds shared by both codes.")
                raise phys.MeshMismatchError(f"{what}: {exc}{hint}") from None
        return dragon_bounds, []
    common = phys.common_mesh(dragon_bounds, serpent_bounds)
    merged: List[Tuple[float, float]] = []
    if len(common) != len(dragon_bounds):
        groups = phys.merged_groups(dragon_bounds, common)
        d = phys.to_ascending(dragon_bounds)
        merged = [(float(d[a]), float(d[b + 1])) for a, b in groups]
        n = len(d) - 1
        warnings.warn(
            f"{what}: DRAGON and Serpent meshes differ; comparing on the "
            f"{len(common) - 1}-group common mesh. Merged DRAGON groups "
            + ", ".join(f"{n - b}-{n - a} ({d[a]:.6g}-{d[b + 1]:.6g} eV)" for a, b in groups)
            + ".", stacklevel=3)
    return common, merged


def _domain_composition(run: SerpentRun, entry: ManifestEntry,
                        rtol: float = 1e-6) -> Dict[str, float]:
    """Nuclide densities of a detector domain, from the Serpent ``mat``
    cards. All domain materials must share the same composition (e.g. the
    zones of a fresh fuel pin)."""
    comps = []
    for name in entry.materials:
        mat = run.input.materials.get(name)
        if mat is None:
            raise KeyError(f"Detector '{entry.detector}': material '{name}' not "
                           "found in the Serpent input.")
        if any(v < 0 for v in mat.nuclides.values()):
            raise ValueError(f"Material '{name}' is given in mass fractions; "
                             "pass serpent_densities explicitly.")
        comps.append(mat.composition())
    ref = comps[0]
    for c in comps[1:]:
        if set(c) != set(ref) or any(abs(c[k] - ref[k]) > rtol * abs(ref[k]) for k in ref):
            raise NotImplementedError(
                f"Detector '{entry.detector}': domain materials have different "
                "compositions; zone-wise detectors are needed for exact rates.")
    return ref


def compute_pin_rates(compo: CompoReader, output: str, serpent: SerpentRun,
                      observable: str, isotopes: Sequence[str], lattice: PinLattice,
                      library: str, sp: StatepointSel = None,
                      pins: Optional[Sequence[int]] = None,
                      serpent_densities: Optional[Mapping[int, Mapping[str, float]]] = None,
                      error_mode: str = "quadrature",
                      check_multiplicity: bool = True,
                      mesh: str = "dragon") -> PinRateComparison:
    """Pin-wise fission and absorption rates of *isotopes* in both codes.

    Args:
        compo: DRAGON COMPO reader.
        output: ``by_pin`` COMPO output (e.g. ``"EDIR_2G"``); its energy mesh
            defines the comparison groups.
        serpent: Serpent run; its manifest must have ``(observable, "rr",
            pin)`` entries.
        observable: Serpent observable (e.g. ``"PIN_295G"``); its mesh must be
            a refinement of the DRAGON mesh.
        isotopes: Isotopes summed in the rates.
        lattice: Pin layout of the DRAGON model.
        library: Q-value table (``physics_data.Q_VALUES``).
        sp: Statepoint selection in the COMPO.
        pins: Restrict to these pins (default: every pin of the lattice).
            Global normalisation then needs the whole assembly.
        serpent_densities: ``{pin: {isotope: N}}`` used with the Serpent
            detectors. Default: the Serpent ``mat`` cards of each detector
            domain (what Serpent actually used). The original analysis used
            the benchmark compositions (``lattice.composition``).
        error_mode: ``"quadrature"``, ``"linear"`` or ``"legacy"``.
        check_multiplicity: Check the Serpent detector multiplicity against
            the model.
        mesh: ``"dragon"`` (default): compare on the DRAGON output mesh, which
            the Serpent grid must refine. ``"common"``: condense both codes
            onto their common bounds (needed when the meshes differ, e.g.
            fine-group comparisons with references made on the legacy
            ``295g`` grid, :data:`physics_data.G295_STARTERDD_LEGACY_EV`).
            The merged ranges are kept in ``merged_eV`` and shaded by the
            plotting functions.

    The reactions compared for an isotope are fission (MT 18 / ``NFTOT``)
    and the absorption reactions of :data:`physics_data.ABSORPTION_REACTIONS`
    present in **both** codes, decided on the first pin that contains the
    isotope. The others are reported in :attr:`PinRateComparison.coverage`.
    """
    if error_mode not in ("quadrature", "linear", "legacy"):
        raise ValueError(f"Unknown error_mode '{error_mode}'.")
    pins = list(pins or lattice.ordered_pins)
    n_mix = compo.n_mixtures(output)
    if n_mix != lattice.n_pins:
        raise ValueError(f"{output} has {n_mix} mixtures but the lattice has "
                         f"{lattice.n_pins} unique pins (is it a by_pin output?).")
    d_bounds = compo.energy_bounds(output, sp)
    ngd = len(d_bounds) - 1
    cond_mode = "quadrature" if error_mode == "quadrature" else "linear"

    entries = {p: serpent.manifest.find(observable, "rr", p) for p in pins}
    raw = {p: serpent.tally_of(e) for p, e in entries.items()}
    bounds, merged = _target_mesh(d_bounds, raw[pins[0]].bounds_eV, mesh,
                                  f"{output} vs {observable}")
    ng = len(bounds) - 1
    d_idx = phys.condensation_indices(d_bounds, bounds)
    tallies = {p: t.condense(bounds, error_mode=cond_mode) for p, t in raw.items()}

    mult = np.array([lattice.multiplicity(p) for p in pins], dtype=float)
    if check_multiplicity:
        for k, p in enumerate(pins):
            m_s = serpent.input.domain_multiplicity(entries[p].materials)
            if m_s is None:
                warnings.warn(f"Cannot read the multiplicity of pin {p} from the "
                              "Serpent lattice; using the model's.", stacklevel=2)
            elif m_s != mult[k]:
                raise ValueError(
                    f"Pin {p}: the Serpent detector covers {m_s} lattice positions "
                    f"but the model has {int(mult[k])} positions for this pin.")
    frac = np.array([lattice.dragon_fraction(p) for p in pins])
    _check_dragon_fractions(compo, output, sp, lattice)

    if serpent_densities is None:
        n_serp = {p: _domain_composition(serpent, entries[p]) for p in pins}
    else:
        n_serp = {p: dict(serpent_densities[p]) for p in pins}

    shape = (len(pins), ng)
    S = {k: np.zeros(shape) for k in ("fission", "absorption", "power",
                                     "fission_err", "absorption_err")}
    # DRAGON rates are built on the native DRAGON mesh, condensed at the end.
    D = {k: np.zeros((len(pins), ngd)) for k in ("fission", "absorption", "power")}
    coverage: Dict[str, ReactionCoverage] = {}

    for iso in isotopes:
        with_iso = [p for p in pins if n_serp[p].get(iso, 0.0) > 0.0]
        cov = ReactionCoverage(iso, with_iso[0] if with_iso else None)
        coverage[iso] = cov
        if not with_iso:
            continue
        ref = with_iso[0]
        ref_entry = entries[ref]
        ref_mix = lattice.mixture_of(ref)

        def in_dragon(key: str) -> bool:
            return compo.isotope_record(output, iso, key, sp, ref_mix, default=None) is not None

        cov.fission = (ref_entry.response_index(phys.FISSION.mt, iso) is not None
                       and in_dragon(phys.FISSION.dragon_key))
        for r in phys.ABSORPTION_REACTIONS:
            in_s = ref_entry.response_index(r.mt, iso) is not None
            in_d = in_dragon(r.dragon_key)
            if in_s and in_d:
                cov.compared.append(r.mt)
        for r in phys.ABSORPTION_REACTIONS:
            in_s = ref_entry.response_index(r.mt, iso) is not None
            in_d = in_dragon(r.dragon_key)
            if in_d and not in_s:
                cov.dragon_only[r.mt] = 0.0
            elif in_s and not in_d:
                cov.serpent_only[r.mt] = 0.0
        q18 = phys.q_value(library, iso, phys.FISSION.mt)

        for k, p in enumerate(pins):
            if p not in with_iso:
                continue
            mix = lattice.mixture_of(p)
            t = tallies[p]
            # ── Serpent (micro responses, condensed) ──
            s_f = np.zeros(ng); s_ferr2 = np.zeros(ng); s_ferr_lin = np.zeros(ng)
            s_a = np.zeros(ng); s_aerr2 = np.zeros(ng); s_aerr_lin = np.zeros(ng)
            s_pw = np.zeros(ng); s_waste = []
            if cov.fission:
                v, rel = t.response(phys.FISSION.mt, iso)
                a = v * rel
                s_f += v; s_a += v; s_pw += v * q18
                s_ferr2 += a ** 2; s_aerr2 += a ** 2
                s_ferr_lin += a; s_aerr_lin += a
            for r in phys.ABSORPTION_REACTIONS:
                resp = t.response(r.mt, iso)
                if resp is None:
                    continue
                v, rel = resp
                a = v * rel
                if r.mt in cov.compared:
                    s_a += r.absorption_factor * v
                    s_pw += v * phys.q_value(library, iso, r.mt)
                    s_aerr2 += (r.absorption_factor * a) ** 2
                    s_aerr_lin += (r.absorption_factor if error_mode == "legacy"
                                   else abs(r.absorption_factor)) * a
                elif r.mt in cov.serpent_only:
                    s_waste.append((r.mt, r.absorption_factor * v))
            n_s = n_serp[p][iso] / mult[k]
            S["fission"][k] += s_f * n_s
            S["absorption"][k] += s_a * n_s
            S["power"][k] += s_pw * n_s
            if error_mode == "quadrature":
                S["fission_err"][k] = np.sqrt(S["fission_err"][k] ** 2 + s_ferr2 * n_s ** 2)
                S["absorption_err"][k] = np.sqrt(S["absorption_err"][k] ** 2 + s_aerr2 * n_s ** 2)
            else:
                S["fission_err"][k] += s_ferr_lin * n_s
                S["absorption_err"][k] += s_aerr_lin * n_s

            # ── DRAGON ──
            def xs(key):
                return compo.micro_xs(output, iso, key, sp, mixtures=mix)
            phi = compo.micro_xs(output, iso, "NWT0", sp, mixtures=mix, missing="raise")
            dens = compo.densities(output, sp, mix)[iso]
            vol = compo.isotope_volumes(output, sp, mix)[iso]
            d_f = np.zeros(ngd); d_a = np.zeros(ngd); d_pw = np.zeros(ngd); d_waste = []
            if cov.fission:
                sf = xs(phys.FISSION.dragon_key)
                d_f += sf; d_a += sf; d_pw += sf * q18
            for r in phys.ABSORPTION_REACTIONS:
                sx = xs(r.dragon_key)
                if sx is None:
                    continue
                if r.mt in cov.compared:
                    d_a += r.absorption_factor * sx
                    d_pw += sx * phys.q_value(library, iso, r.mt)
                elif r.mt in cov.dragon_only:
                    d_waste.append((r.mt, r.absorption_factor * sx))
            scale = dens * vol / frac[k]
            D["fission"][k] += d_f * phi * scale
            D["absorption"][k] += d_a * phi * scale
            D["power"][k] += d_pw * phi * scale

            # ── coverage: share of the total absorption missing in one code ──
            for waste, tot_base, store in ((d_waste, d_a, cov.dragon_only),
                                           (s_waste, s_a, cov.serpent_only)):
                if not waste:
                    continue
                tot = tot_base + sum(w for _, w in waste)
                with np.errstate(divide="ignore", invalid="ignore"):
                    for mt, w in waste:
                        share = float(np.nanmax(np.abs(w) / np.abs(tot)))
                        store[mt] = max(store[mt], share)

    if ngd != ng:
        D = {k: phys.condense(v, d_idx, axis=1) for k, v in D.items()}
    serp = CodeRates("serpent", pins, bounds, S["fission"], S["absorption"], S["power"],
                     mult, S["fission_err"], S["absorption_err"])
    drag = CodeRates("dragon", pins, bounds, D["fission"], D["absorption"], D["power"], mult)
    return PinRateComparison(drag, serp, lattice, list(isotopes), coverage,
                             output, observable, error_mode, merged)


def _check_dragon_fractions(compo: CompoReader, output: str, sp: StatepointSel,
                            lattice: PinLattice, rtol: float = 1e-3) -> None:
    """Warn if COMPO mixture volumes contradict the assumed DRAGON fractions
    (compares on-axis and off-axis pins of the same fuel material)."""
    if not lattice.axis_pins:
        return
    vol = {p: compo.mixture_volume(output, sp, lattice.mixture_of(p))
           for p in lattice.ordered_pins}
    for p in lattice.axis_pins:
        same = [q for q in lattice.ordered_pins if not lattice.on_axis(q)
                and lattice.material.get(q) == lattice.material.get(p)]
        if not same:
            continue
        ratio = vol[p] / vol[same[0]]
        if abs(ratio - lattice.dragon_fraction(p)) > rtol:
            warnings.warn(
                f"{output}: volume of on-axis pin {p} is {ratio:.4f} x that of pin "
                f"{same[0]} (expected {lattice.dragon_fraction(p)}); check the "
                "DRAGON symmetry convention.", stacklevel=3)
            return


# ═════════════════════════════════════════════════════════════════════════════
# Spectra
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class SpectrumComparison:
    """Flux spectra normalised to unit sum (rows = regions)."""

    bounds_eV: np.ndarray
    labels: List[str]
    dragon: np.ndarray
    serpent: np.ndarray
    serpent_rel_err: np.ndarray
    #: As :attr:`PinRateComparison.merged_eV`.
    merged_eV: List[Tuple[float, float]] = field(default_factory=list)

    def per_lethargy(self) -> Tuple[np.ndarray, np.ndarray]:
        du = phys.lethargy_widths(self.bounds_eV)
        return self.dragon / du, self.serpent / du

    def errors(self) -> np.ndarray:
        """(DRAGON - Serpent) / Serpent in %, per group."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return (self.dragon - self.serpent) / self.serpent * 100.0


def compare_spectra(compo: CompoReader, output: str, serpent: SerpentRun,
                    observable: str, sp: StatepointSel = None,
                    pins: Optional[Sequence[int]] = None,
                    lattice: Optional[PinLattice] = None,
                    error_mode: str = "quadrature",
                    mesh: str = "dragon") -> SpectrumComparison:
    """Compare ``NWT0`` spectra with Serpent flux detectors.

    Without *pins*, the COMPO output must have one mixture (e.g. ``ALL``
    homogenisation) and the Serpent observable one flux detector. With
    *pins* (and *lattice*), each pin's mixture is compared with the
    ``(observable, "flux", pin)`` detector. *mesh* as in
    :func:`compute_pin_rates`.
    """
    d_bounds = compo.energy_bounds(output, sp)
    if pins is None:
        d = compo.flux(output, sp, mixtures=1)[None, :]
        raw = [serpent.tally(observable, "flux")]
        labels = [output]
    else:
        if lattice is None:
            raise ValueError("lattice is required with pins.")
        d = np.vstack([compo.flux(output, sp, mixtures=lattice.mixture_of(p)) for p in pins])
        raw = [serpent.tally(observable, "flux", p) for p in pins]
        labels = [f"{p}: {lattice.material.get(p, '')}" for p in pins]
    bounds, merged = _target_mesh(d_bounds, raw[0].bounds_eV, mesh, f"{output} vs {observable}")
    if len(bounds) != len(d_bounds):
        d = phys.condense(d, phys.condensation_indices(d_bounds, bounds), axis=1)
    ts = [t.condense(bounds, error_mode=error_mode) for t in raw]
    s = np.vstack([t.values[0] for t in ts])
    se = np.vstack([t.rel_err[0] for t in ts])
    return SpectrumComparison(bounds, labels, d / d.sum(axis=1, keepdims=True),
                              s / s.sum(axis=1, keepdims=True), se, merged)


__all__ = [
    "KeffComparison", "compare_keff", "CodeRates", "ReactionCoverage",
    "PinRateComparison", "NormalisationError", "compute_pin_rates",
    "SpectrumComparison", "compare_spectra", "QUANTITIES",
]
