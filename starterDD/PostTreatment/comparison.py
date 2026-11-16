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

**Normalisation.** DRAGON is scaled by one factor ``f = T_serpent / T_dragon``
where ``T`` is a full-assembly total. Two families (plan §3.6):

- *whole-assembly* methods (:data:`GLOBAL_NORMALISATIONS`), independent of
  the compared isotopes, computed once per statepoint from macroscopic
  (all-nuclide) rates by :func:`global_normalisation` and passed as
  ``factor=``: ``fission_source`` (Σ ν·fission over the fuel, the usual
  k-eigenvalue normalisation), ``fission`` and ``absorption`` (whole
  assembly, ALL output);
- *subset* methods, computed by each comparison from its own rates:
  ``power_subset`` (Q-weighted rates of the studied isotopes; the legacy
  analysis, also accepted as ``"power"``), ``absorption_subset``,
  ``fission_subset`` and ``pin_average`` (same total of the compared
  quantity: shape only). They depend on which isotopes are compared.

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


#: Whole-assembly normalisations: ``{method: (spatial, macroscopic quantity)}``
#: of the macroscopic rates that give the total (:func:`global_normalisation`).
GLOBAL_NORMALISATIONS: Dict[str, Tuple[str, str]] = {
    "fission_source": ("PIN", "nufission"),
    "fission": ("PIN", "fission"),
    "absorption": ("ALL", "absorption"),
}
#: Normalisations computed by a comparison from its own rates.
SUBSET_NORMALISATIONS = ("power_subset", "absorption_subset", "fission_subset", "pin_average")
#: Legacy names accepted by :class:`PinRateComparison`.
_SUBSET_ALIASES = {"power": "power_subset"}
NORMALISATIONS = tuple(GLOBAL_NORMALISATIONS) + SUBSET_NORMALISATIONS


@dataclass
class Normalisation:
    """A whole-assembly normalisation of one statepoint.

    ``factor`` multiplies DRAGON values; ``dragon_total`` and
    ``serpent_total`` are full-assembly totals (DRAGON already scaled from
    the modelled fraction of the assembly to the full one)."""

    method: str
    factor: float
    dragon_total: float
    serpent_total: float
    serpent_rel_err: float
    observable: str = ""
    output: str = ""
    note: str = ""


def global_normalisation(method: str, cmp: "MacroRateComparison",
                         geometry_fraction: float = 1.0) -> Normalisation:
    """Whole-assembly normalisation from macroscopic rates.

    Args:
        method: A :data:`GLOBAL_NORMALISATIONS` key.
        cmp: Macroscopic rates of the matching observable: per pin (``PIN``,
            ``nufission`` or ``fission``) for ``fission_source`` / ``fission``,
            the assembly (``ALL``, ``absorption``) for ``absorption``.
        geometry_fraction: Fraction of the full assembly modelled by DRAGON
            (:attr:`PinLattice.dragon_geometry_fraction`, 1/2 for a diagonal
            symmetry). Used for ``FUEL``/``ALL`` outputs only: per-pin rates
            are already per physical pin in both codes.

    Physics: in an infinite lattice, production / k∞ = absorption. With
    ``fission_source`` both codes describe one source neutron, so their total
    absorptions differ by about -Δk/k (the keff bias appears in the
    absorption rates). With ``absorption`` the balance is matched and the
    bias goes to production. A DRAGON flux with a leakage model and a
    non-zero buckling adds a leakage term that Serpent's infinite lattice
    does not have, which biases ``absorption``; prefer ``fission_source``.
    """
    if method not in GLOBAL_NORMALISATIONS:
        raise ValueError(f"'{method}' is not a whole-assembly normalisation "
                         f"({sorted(GLOBAL_NORMALISATIONS)}).")
    spatial, q = GLOBAL_NORMALISATIONS[method]
    if cmp.spatial != spatial:
        raise NormalisationError(f"{method} needs a {spatial} observable, got "
                                 f"{cmp.observable} ({cmp.spatial}).")
    if q not in cmp.dragon:
        raise NormalisationError(f"{method}: '{q}' is not scored by {cmp.observable} "
                                 f"({cmp.quantities}).")
    scale = 1.0 if spatial == "PIN" else 1.0 / float(geometry_fraction)
    d = float(cmp.total("dragon", q).sum()) * scale
    s = float(cmp.total("serpent", q).sum())
    err = float(np.sqrt(np.sum((cmp.multiplicity[:, None] * cmp.serpent_err[q]) ** 2)))
    if d <= 0 or s <= 0:
        raise NormalisationError(f"{method}: non-positive total (DRAGON {d}, Serpent {s}).")
    return Normalisation(method, s / d, d, s, err / s, cmp.observable, cmp.output)


def _subset_method(method: str) -> str:
    return _SUBSET_ALIASES.get(method, method)


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
    def normalisation_factor(self, method: str = "power_subset",
                             quantity: str = "absorption") -> float:
        """Factor applied to DRAGON rates so that both codes share the same
        full-assembly total of the studied isotopes: ``power_subset`` (alias
        ``power``), ``absorption_subset``, ``fission_subset``, or
        ``pin_average`` (the total of *quantity*: shape only). The
        whole-assembly methods (:data:`GLOBAL_NORMALISATIONS`) do not depend
        on the isotopes: compute them with :func:`global_normalisation` and
        pass ``factor=`` to :meth:`errors` / :meth:`stats`."""
        method = _subset_method(method)
        if method in GLOBAL_NORMALISATIONS:
            raise NormalisationError(
                f"'{method}' is a whole-assembly normalisation: compute it with "
                "global_normalisation() from macroscopic rates and pass factor=.")
        if not self.is_complete:
            raise NormalisationError(
                f"Rates cover {len(self.pins)} of {self.lattice.n_pins} pins; "
                "a global normalisation needs the whole assembly.")
        q = {"power_subset": "power", "absorption_subset": "absorption",
             "fission_subset": "fission", "pin_average": quantity}.get(method)
        if q is None:
            raise ValueError(f"Unknown normalisation method '{method}' "
                             f"({', '.join(NORMALISATIONS)}).")
        return float(self.serpent.total(q).sum() / self.dragon.total(q).sum())

    def errors(self, quantity: str = "absorption", method: str = "power_subset",
               factor: Optional[float] = None) -> np.ndarray:
        """Relative difference DRAGON vs Serpent, in %, ``(n_pins, n_groups)``.

        *factor* (a whole-assembly :class:`Normalisation`'s ``factor``)
        overrides *method*. ``fission_over_absorption`` is
        normalisation-free. NaN where Serpent scores zero.
        """
        s = self.serpent.get(quantity)
        if quantity == "fission_over_absorption":
            d = self.dragon.get(quantity)
        else:
            k = factor if factor is not None else self.normalisation_factor(method, quantity)
            d = self.dragon.get(quantity) * k
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
    def stats(self, quantity: str = "absorption", method: str = "power_subset",
              factor: Optional[float] = None) -> List[Dict[str, float]]:
        """Per-group statistics over **physical** pins (each unique pin
        weighted by its multiplicity). Same definitions as the original
        analysis: min |err|, min, max, mean, mean weighted by the Serpent
        distribution, RMS and weighted RMS (no weighted values for
        ``fission_over_absorption``)."""
        err = self.errors(quantity, method, factor)
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
# Macroscopic (all-nuclide) rates: Serpent ``void`` responses
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class MacroRateComparison:
    """Macroscopic rates of both codes, per region (pins, or one region for
    FUEL/ALL observables), ascending energy, arrays ``(n_regions, n_groups)``.

    Per physical pin for PIN observables (same convention as
    :class:`PinRateComparison`). For FUEL/ALL the codes cover different
    geometries (DRAGON's DIAG half vs Serpent's full assembly), so only
    normalised values are compared.
    """

    observable: str
    output: str
    spatial: str
    regions: List[Optional[int]]
    bounds_eV: np.ndarray
    multiplicity: np.ndarray
    dragon: Dict[str, np.ndarray]
    serpent: Dict[str, np.ndarray]
    serpent_err: Dict[str, np.ndarray]
    merged_eV: List[Tuple[float, float]] = field(default_factory=list)

    @property
    def quantities(self) -> List[str]:
        return list(self.dragon)

    @property
    def n_groups(self) -> int:
        return len(self.bounds_eV) - 1

    def total(self, code: str, quantity: str) -> np.ndarray:
        """Per-group total over regions (``sum_r multiplicity * rate``)."""
        a = (self.dragon if code == "dragon" else self.serpent)[quantity]
        return (self.multiplicity[:, None] * a).sum(axis=0)

    def normalisation_factor(self, method: Optional[str] = None,
                             quantity: Optional[str] = None) -> float:
        """Factor applied to DRAGON so that both codes have the same total of
        *method* over these regions: a compared quantity name, or
        ``pin_average`` (*quantity* itself: shape only). Default: ``fission``
        if compared, else *quantity*. For a whole-assembly normalisation pass
        ``factor=`` to :meth:`errors` (:func:`global_normalisation`)."""
        if method == "pin_average":
            method = quantity
        method = method or ("fission" if "fission" in self.dragon else quantity)
        if method not in self.dragon:
            raise NormalisationError(f"'{method}' is not compared ({self.quantities}).")
        return float(self.total("serpent", method).sum() / self.total("dragon", method).sum())

    def errors(self, quantity: str, method: Optional[str] = None,
               factor: Optional[float] = None) -> np.ndarray:
        """DRAGON vs Serpent in %, ``(n_regions, n_groups)``. *factor*
        overrides *method*."""
        s = self.serpent[quantity]
        k = factor if factor is not None else self.normalisation_factor(method, quantity)
        d = self.dragon[quantity] * k
        with np.errstate(divide="ignore", invalid="ignore"):
            e = (d - s) / s * 100.0
        e[~np.isfinite(e)] = np.nan
        return e

    def rel_sigma(self, quantity: str) -> np.ndarray:
        """Serpent 1σ in %."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return 100.0 * self.serpent_err[quantity] / self.serpent[quantity]

    def stats(self, quantity: str, method: Optional[str] = None,
              factor: Optional[float] = None) -> List[Dict[str, float]]:
        """Per-group min/max/mean/RMS over physical regions."""
        err = self.errors(quantity, method, factor)
        out = []
        for g in range(self.n_groups):
            e = err[:, g]
            ok = np.isfinite(e)
            w, eg = self.multiplicity[ok], e[ok]
            out.append({"min": float(eg.min()), "max": float(eg.max()),
                        "mean": float(np.sum(w * eg) / w.sum()),
                        "rms": float(np.sqrt(np.sum(w * eg ** 2) / w.sum()))})
        return out


def compute_macro_rates(compo: CompoReader, output: str, serpent: SerpentRun,
                        observable: str, quantities: Optional[Sequence[str]] = None,
                        lattice: Optional[PinLattice] = None, sp: StatepointSel = None,
                        pins: Optional[Sequence[int]] = None,
                        error_mode: str = "quadrature", mesh: str = "dragon",
                        check_multiplicity: bool = True) -> MacroRateComparison:
    """Macroscopic rates (:data:`physics_data.MACRO_REACTIONS`) of a DRAGON
    output against Serpent ``dr <mt> void`` responses.

    DRAGON: ``V · φ_g · Σ_i N_i Σ_k σ_i,k,g`` over the edited isotopes **and**
    ``*MAC*RES`` (density 1, the macroscopic residual of every other
    nuclide), so every nuclide of the region is included, as in Serpent's
    collision-material response. PIN rates are per physical pin
    (DRAGON ÷ ``dragon_fraction``, Serpent ÷ multiplicity).

    Args:
        output: DRAGON output (``by_pin`` for a PIN observable, one mixture
            for FUEL/ALL).
        observable: Serpent observable with ``void`` responses
            (``ComparisonSpec.serpent_source`` / ``SerpentReference.source``).
        quantities: Default: every macroscopic quantity whose MTs the
            Serpent detectors score.
        lattice: Required for PIN observables.
        pins: Restrict PIN observables to these pins.
    """
    if error_mode not in ("quadrature", "linear"):
        raise ValueError(f"Unknown error_mode '{error_mode}'.")
    man = serpent.manifest
    spatial = man.select(observable, "rr")[0].spatial if man.select(observable, "rr") else None
    if spatial is None:
        raise KeyError(f"No '{observable}' reaction-rate detector in the Serpent manifest.")
    if spatial == "PIN":
        if lattice is None:
            raise ValueError("A PIN observable needs the lattice.")
        regions: List[Optional[int]] = list(pins or lattice.ordered_pins)
        entries = {p: man.find(observable, "rr", p) for p in regions}
        mixtures = {p: lattice.mixture_of(p) for p in regions}
        frac = np.array([lattice.dragon_fraction(p) for p in regions])
        mult = np.array([lattice.multiplicity(p) for p in regions], dtype=float)
        if compo.n_mixtures(output) != lattice.n_pins:
            raise ValueError(f"{output} is not a by_pin output of this lattice.")
        _check_dragon_fractions(compo, output, sp, lattice)
    else:
        if compo.n_mixtures(output) != 1:
            raise ValueError(f"{output} has {compo.n_mixtures(output)} mixtures; a "
                             f"{spatial} observable needs a homogenised output.")
        regions = [None]
        entries = {None: man.find(observable, "rr", None)}
        mixtures = {None: 1}
        frac = np.ones(1)
        mult = np.ones(1)

    void_mts = {r.mt for r in entries[regions[0]].responses if r.material == phys.VOID}
    if quantities is None:
        quantities = [q for q, m in phys.MACRO_REACTIONS.items()
                      if set(m.serpent_mts) <= void_mts]
    missing = [q for q in quantities if not set(phys.MACRO_REACTIONS[q].serpent_mts) <= void_mts]
    if not quantities or missing:
        raise KeyError(f"{observable}: no Serpent 'void' responses for {missing or 'any rate'} "
                       f"(scored MTs: {sorted(void_mts)}).")

    raw = {r: serpent.tally_of(e) for r, e in entries.items()}
    d_bounds = compo.energy_bounds(output, sp)
    bounds, merged = _target_mesh(d_bounds, raw[regions[0]].bounds_eV, mesh,
                                  f"{output} vs {observable}")
    d_idx = phys.condensation_indices(d_bounds, bounds)
    tallies = {r: t.condense(bounds, error_mode=error_mode) for r, t in raw.items()}
    if check_multiplicity and spatial == "PIN":
        for k, p in enumerate(regions):
            m_s = serpent.input.domain_multiplicity(entries[p].materials)
            if m_s is not None and m_s != mult[k]:
                raise ValueError(f"Pin {p}: Serpent detector covers {m_s} positions, "
                                 f"the model {int(mult[k])}.")

    ng, ngd = len(bounds) - 1, len(d_bounds) - 1
    D = {q: np.zeros((len(regions), ngd)) for q in quantities}
    S = {q: np.zeros((len(regions), ng)) for q in quantities}
    E = {q: np.zeros((len(regions), ng)) for q in quantities}
    for k, r in enumerate(regions):
        mix = mixtures[r]
        dens = compo.densities(output, sp, mix)
        vol = compo.mixture_volume(output, sp, mix)
        phi = compo.flux(output, sp, mixtures=mix)
        xs_cache: Dict[Tuple[str, str], np.ndarray] = {}
        for q in quantities:
            macro = np.zeros(ngd)
            for iso, n in dens.items():
                for key in phys.MACRO_REACTIONS[q].dragon_keys:
                    if (iso, key) not in xs_cache:
                        v = compo.micro_xs(output, iso, key, sp, mixtures=mix)
                        xs_cache[(iso, key)] = np.zeros(ngd) if v is None else v
                    macro += n * xs_cache[(iso, key)]
            D[q][k] = macro * phi * vol / frac[k]
            t = tallies[r]
            v_sum = np.zeros(ng)
            a2 = np.zeros(ng)
            a_lin = np.zeros(ng)
            for mt in phys.MACRO_REACTIONS[q].serpent_mts:
                v, rel = t.response(mt, None)
                v_sum += v
                a2 += (v * rel) ** 2
                a_lin += np.abs(v * rel)
            S[q][k] = v_sum / mult[k]
            E[q][k] = (np.sqrt(a2) if error_mode == "quadrature" else a_lin) / mult[k]
    if ngd != ng:
        D = {q: phys.condense(v, d_idx, axis=1) for q, v in D.items()}
    return MacroRateComparison(observable, output, spatial, regions, bounds, mult,
                               D, S, E, merged)


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
    "Normalisation", "global_normalisation", "GLOBAL_NORMALISATIONS",
    "SUBSET_NORMALISATIONS", "NORMALISATIONS",
    "MacroRateComparison", "compute_macro_rates",
    "SpectrumComparison", "compare_spectra", "QUANTITIES",
]
