# -*- coding: utf-8 -*-
"""Comparison specification: one declaration of the observables, from which
the DRAGON outputs, the Serpent2 detectors and the post-processing pairings
are all derived (plan §3.1-§3.3).

An :class:`Observable` is named ``{SPATIAL}_{MESH}[{TAG}]`` (≤ 12 characters,
also the EDI ``SAVE ON`` / COMPO ``STEP UP`` name):

==========  ==========================  ==========================================
spatial     DRAGON EDI mode             Serpent2 domain
==========  ==========================  ==========================================
``PIN``     ``by_pin``                  one detector per pin, ``dm`` its fuel zones
``FUEL``    ``FUEL``                    one detector, ``dm`` every fuel zone
``ALL``     ``ALL``                     no domain (flux only)
==========  ==========================  ==========================================

Meshes are preset names (``1g``, ``2g``, ``26g``, ``295g``: the DRAGON library
meshes, :mod:`physics_data`) or explicit bounds in eV.

DRAGON side (:meth:`ComparisonSpec.to_dragon_outputs`). The ``COND`` of each
output is written as **group indices** of the scheme's output mesh (SHEM295
for a 1L scheme, the L1-condensed library for a 2L scheme). ``EDI:`` maps a
real energy to the library bound *at or below* it (``EDIENE.f``) and silently
ignores a real that is not decreasing (``EDIGET.f``), so indices are the only
exact form. An observable whose mesh is not a coarsening of the output mesh
cannot be edited (e.g. ``PIN_295G`` in a 2L scheme) and raises.

Serpent side (:meth:`ComparisonSpec.configure_serpent`). Each observable is
scored on its own mesh, unless another observable of the spec already
provides it: same spatial and kind, a mesh that refines it, at least the same
isotopes and pins (:meth:`ComparisonSpec.serpent_source`). Detectors follow the
naming convention (``PIN_295G_RR_p12``, ``FUEL_295G_RR``, ``ALL_295G_FLX``)
and every detector gets a manifest entry (``source="exporter"``), written next
to the input by :meth:`Serpent2Model.write`.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from . import physics_data as phys
from .serpent_reader import DetectorManifest, ManifestEntry, Response

MeshLike = Union[str, Sequence[float]]

SPATIAL_TO_DRAGON = {"PIN": "by_pin", "FUEL": "FUEL", "ALL": "ALL"}
DRAGON_TO_SPATIAL = {v: k for k, v in SPATIAL_TO_DRAGON.items()}
SERPENT_KINDS = ("rr", "flux")
NAME_MAX = 12
#: Named pin subsets of :attr:`Observable.pins`, resolved on the lattice
#: (:meth:`PinLattice.diagonal_pins`): ``"axis"`` (alias ``"diagonal"``) is the
#: symmetry axis of the DRAGON geometry, ``"main"`` / ``"anti"`` the main
#: (row == col) / anti (row + col == n-1) diagonal.
PIN_SELECTORS = {"axis": "axis", "diagonal": "axis", "main": "main", "anti": "anti"}
#: Macroscopic rates added by the presets (``global_totals=True``) for the
#: whole-assembly normalisations of the post-processing:
#: ``fission_source`` / ``fission`` need the per-pin totals, ``absorption``
#: the assembly (``ALL``) total.
GLOBAL_PIN_MACRO = ("fission", "nufission")
GLOBAL_ALL_MACRO = ("absorption",)
_NAME_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")

#: DRAGON output names of the GE14 legacy schemes -> spatial code. The mesh
#: code is taken from the COMPO ENERGY record (``EDIR_295`` of a 2L scheme
#: holds 26 groups, so it is ``PIN_26G``).
LEGACY_OUTPUT_SPATIAL = {"EDIR_2G": "PIN", "EDIR_295": "PIN",
                         "EDIHOM_COND": "FUEL", "EDIHOM_295": "ALL"}


# ─────────────────────────────────────────────────────────────────────────────
# Meshes
# ─────────────────────────────────────────────────────────────────────────────


def resolve_mesh(mesh: MeshLike) -> np.ndarray:
    """Ascending bounds in eV of a preset name or explicit bounds."""
    if isinstance(mesh, str):
        if mesh.lower() == "1g":
            return phys.G2_EV[[0, -1]].copy()
        return phys.mesh(mesh)
    b = phys.to_ascending(mesh)
    if len(b) < 2 or np.any(np.diff(b) <= 0):
        raise ValueError(f"Invalid energy bounds {list(mesh)}.")
    return b


def mesh_code(bounds: Sequence[float]) -> str:
    """``1G``, a preset code (``2G``, ``26G``, ``295G``) or ``C{n}G``."""
    return "1G" if len(bounds) == 2 else phys.mesh_code(bounds)


def refines(fine: Sequence[float], coarse: Sequence[float],
            rtol: float = phys.BOUNDS_RTOL) -> bool:
    """True if every bound of *coarse* is a bound of *fine*."""
    f = np.asarray(fine, dtype=float)
    return all(np.any(np.isclose(f, c, rtol=rtol, atol=0.0)) for c in coarse)


def scheme_output_mesh(scheme, library_mesh: Sequence[float] = phys.G295_EV) -> np.ndarray:
    """Ascending eV bounds of the library on which the scheme's outputs are
    edited: *library_mesh* for a 1L scheme, condensed by every
    ``edition_between_levels`` step (``COND`` group indices) otherwise."""
    from ..DDModel.DragonCalculationScheme import EditionBetweenLevelsStep
    desc = np.asarray(library_mesh, dtype=float)[::-1]
    for step in scheme.steps:
        if isinstance(step, EditionBetweenLevelsStep):
            idx = [int(i) for i in step.energy_groups_bounds if 0 < int(i) < len(desc) - 1]
            desc = np.concatenate([desc[:1], desc[idx], desc[-1:]])
    return desc[::-1].copy()


def dragon_cond(bounds: Sequence[float], output_mesh: Sequence[float]):
    """``energy_bounds`` of a :class:`CalculationOutput` editing *bounds*
    from *output_mesh* (both ascending eV): ``None`` (no ``COND``), ``[]``
    (one group) or DRAGON group indices (``IGCR``: last fine group, counted
    from the highest energy, of every coarse group but the last)."""
    b = np.asarray(bounds, dtype=float)
    out = np.asarray(output_mesh, dtype=float)
    if len(b) == len(out) and phys.bounds_match(b, out):
        return None
    if not (np.isclose(b[0], out[0], rtol=phys.BOUNDS_RTOL)
            and np.isclose(b[-1], out[-1], rtol=phys.BOUNDS_RTOL)):
        raise phys.MeshMismatchError(
            f"Mesh [{b[0]:.6g}, {b[-1]:.6g}] eV does not span the output mesh "
            f"[{out[0]:.6g}, {out[-1]:.6g}] eV.")
    if len(b) == 2:
        return []
    desc = out[::-1]
    idx, missing = [], []
    for e in b[-2:0:-1]:                       # inner bounds, descending
        hit = np.nonzero(np.isclose(desc, e, rtol=phys.BOUNDS_RTOL, atol=0.0))[0]
        if len(hit) == 0:
            missing.append(e)
        else:
            idx.append(int(hit[0]))
    if missing:
        raise phys.MeshMismatchError(
            f"{len(missing)} bound(s) are not bounds of the {len(out) - 1}-group output "
            f"mesh, e.g. {[float(f'{m:.6g}') for m in missing[:5]]} eV.")
    return idx


def cond_to_bounds(energy_bounds, output_mesh: Sequence[float]) -> np.ndarray:
    """Ascending eV bounds produced by a ``CalculationOutput.energy_bounds``
    on *output_mesh*: ``None`` (no ``COND``), ``[]`` (1 group), integer group
    indices, or real energies in eV mapped as ``EDIENE.f`` does (each energy
    to the library bound at or below it)."""
    out = np.asarray(output_mesh, dtype=float)
    if energy_bounds is None:
        return out.copy()
    if len(energy_bounds) == 0:
        return out[[0, -1]].copy()
    desc = out[::-1]
    idx = []
    for b in energy_bounds:
        if isinstance(b, (int, np.integer)):
            j = int(b)
        else:
            below = np.nonzero(np.asarray(desc, dtype=np.float32) <= np.float32(b))[0]
            j = int(below[0]) if len(below) else len(desc) - 1
        if 0 < j < len(desc) - 1 and j not in idx:
            idx.append(j)
    return np.concatenate([desc[:1], desc[sorted(idx)], desc[-1:]])[::-1].copy()


# ─────────────────────────────────────────────────────────────────────────────
# Observables
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class Observable:
    """One quantity compared between the codes.

    Args:
        spatial: ``PIN``, ``FUEL`` or ``ALL`` (DRAGON modes ``by_pin`` /
            ``FUEL`` / ``ALL`` accepted).
        mesh: Preset name or ascending/descending bounds in eV.
        isotopes: Isotopes edited by DRAGON (``MICR``) and scored by Serpent.
        serpent: Serpent quantities, any of ``"rr"`` and ``"flux"``; empty for
            a DRAGON-only output (e.g. the 1-group fuel inventory edition).
        dragon: Produce a DRAGON output.
        pins: Serpent detector subset (``PIN`` only): a list of ``pin_idx``
            or a :data:`PIN_SELECTORS` name (``"diagonal"``: the symmetry
            axis), resolved on the lattice (:meth:`resolve_pins`). DRAGON
            always edits every pin; only the Serpent detectors are limited,
            to keep large fine-mesh outputs small.
        statepoints: ``"ALL"`` or ``{branch_type: [values]}`` (DRAGON).
        tag: Suffix of the name, ``[A-Z0-9_]``.
        reactions: ``{isotope: [MT, ...]}`` overriding
            :func:`physics_data.serpent_detector_mts` for the Serpent responses.
        comment: COMPO comment.
        macro: Macroscopic (all-nuclide) rates scored by Serpent ``dr <mt>
            void`` responses, any of :data:`physics_data.MACRO_REACTIONS`
            (``fission``, ``capture``, ``absorption`` = capture + fission,
            ``nufission``, ``total``). DRAGON rebuilds them from the edited
            isotopes plus ``*MAC*RES`` (``compute_macro_rates``).
        serpent_isotopes: Isotopes given per-isotope Serpent responses
            (default: *isotopes*). ``()`` keeps only the macroscopic ones,
            which makes small detectors.
        output: DRAGON output name (EDI ``SAVE ON`` / COMPO ``STEP UP``),
            default: the observable name. Lets a spec describe the outputs of
            an existing scheme (``from_scheme``).
    """

    spatial: str
    mesh: MeshLike
    isotopes: Sequence[str] = ()
    serpent: Sequence[str] = ("rr",)
    dragon: bool = True
    pins: Optional[Union[str, Sequence[int]]] = None
    statepoints: Union[str, Mapping] = "ALL"
    tag: str = ""
    reactions: Optional[Mapping[str, Sequence[int]]] = None
    comment: str = ""
    macro: Sequence[str] = ()
    serpent_isotopes: Optional[Sequence[str]] = None
    output: Optional[str] = None

    def __post_init__(self):
        sp = DRAGON_TO_SPATIAL.get(self.spatial, str(self.spatial).upper())
        if sp not in SPATIAL_TO_DRAGON:
            raise ValueError(f"Unknown spatial '{self.spatial}' (PIN, FUEL, ALL).")
        self.spatial = sp
        self.isotopes = list(self.isotopes)
        self.serpent = tuple(self.serpent)
        bad = set(self.serpent) - set(SERPENT_KINDS)
        if bad:
            raise ValueError(f"Unknown Serpent kind(s) {sorted(bad)}; use {SERPENT_KINDS}.")
        if self.pins is not None:
            if sp != "PIN":
                raise ValueError("pins= only applies to PIN observables.")
            if isinstance(self.pins, str):
                key = self.pins.lower()
                if key not in PIN_SELECTORS:
                    raise ValueError(f"Unknown pin selector '{self.pins}'; use a list of "
                                     f"pin_idx or one of {sorted(PIN_SELECTORS)}.")
                self.pins = PIN_SELECTORS[key]
            else:
                self.pins = sorted(int(p) for p in self.pins)
        self.macro = tuple(dict.fromkeys(self.macro))
        bad = set(self.macro) - set(phys.MACRO_REACTIONS)
        if bad:
            raise ValueError(f"Unknown macroscopic rate(s) {sorted(bad)}; use "
                             f"{sorted(phys.MACRO_REACTIONS)}.")
        if self.serpent_isotopes is not None:
            self.serpent_isotopes = list(self.serpent_isotopes)
        if sp == "ALL" and "rr" in self.serpent and self.micro_isotopes:
            raise ValueError("ALL observables can only score macroscopic rates (macro=, "
                             "serpent_isotopes=()) or the flux in Serpent; use FUEL for "
                             "per-isotope assembly rates.")
        if "rr" in self.serpent and not (self.micro_isotopes or self.macro):
            raise ValueError("Serpent 'rr' needs isotopes (serpent_isotopes) or macro rates.")
        self.bounds_eV = resolve_mesh(self.mesh)
        for label, name in (("Observable", self.name), ("Output", self.dragon_name)):
            if len(name) > NAME_MAX or not _NAME_RE.match(name):
                raise ValueError(f"{label} name '{name}' must match [A-Z][A-Z0-9_]* and have "
                                 f"≤ {NAME_MAX} characters (EDI/COMPO directory name).")
        if self.dragon and not self.isotopes:
            raise ValueError(f"{self.name}: a DRAGON output needs isotopes (EDI MICR); the "
                             "others are kept in *MAC*RES.")

    @property
    def pins_unresolved(self) -> bool:
        """True while :attr:`pins` is a selector name (see :meth:`resolve_pins`)."""
        return isinstance(self.pins, str)

    def resolve_pins(self, lattice) -> "Observable":
        """Copy whose selector :attr:`pins` (``"axis"``, ``"main"``,
        ``"anti"``) is replaced by the ``pin_idx`` list of *lattice*
        (a :class:`PinLattice`). Unchanged copy otherwise."""
        kw = _obs_kwargs(self)
        if self.pins_unresolved:
            pins = lattice.diagonal_pins(self.pins)
            if not pins:
                raise ValueError(f"{self.name}: no fuel pin on the '{self.pins}' diagonal.")
            kw["pins"] = pins
        return Observable(**kw)

    @property
    def micro_isotopes(self) -> List[str]:
        """Isotopes with per-isotope Serpent responses."""
        return list(self.isotopes if self.serpent_isotopes is None else self.serpent_isotopes)

    @property
    def dragon_name(self) -> str:
        """DRAGON output name."""
        return self.output or self.name

    @property
    def mesh_code(self) -> str:
        return mesh_code(self.bounds_eV)

    @property
    def name(self) -> str:
        return f"{self.spatial}_{self.mesh_code}{self.tag}"

    @property
    def n_groups(self) -> int:
        return len(self.bounds_eV) - 1


@dataclass
class SerpentRequest:
    """Detectors of one ``(observable, kind)`` to generate."""

    observable: str
    kind: str
    spatial: str
    bounds_eV: np.ndarray
    isotopes: List[str]
    pins: Optional[List[int]]
    reactions: Optional[Mapping[str, Sequence[int]]] = None
    macro: Tuple[str, ...] = ()

    @property
    def grid_name(self) -> Optional[str]:
        return None if len(self.bounds_eV) == 2 else mesh_code(self.bounds_eV).lower()

    @property
    def macro_mts(self) -> List[int]:
        """Serpent MTs of the ``void`` responses, ordered, unique."""
        return list(dict.fromkeys(mt for q in self.macro
                                  for mt in phys.MACRO_REACTIONS[q].serpent_mts))

    def provides(self, other: "SerpentRequest") -> bool:
        """This request's detectors also give *other* (by condensation)."""
        return (self.spatial == other.spatial and self.kind == other.kind
                and refines(self.bounds_eV, other.bounds_eV)
                and (self.kind == "flux" or set(other.isotopes) <= set(self.isotopes))
                and (self.kind == "flux" or set(other.macro_mts) <= set(self.macro_mts))
                and (self.pins is None or (other.pins is not None
                                           and set(other.pins) <= set(self.pins)))
                and (other.reactions is None or other.reactions == self.reactions))


# ─────────────────────────────────────────────────────────────────────────────
# Spec
# ─────────────────────────────────────────────────────────────────────────────


class ComparisonSpec:
    """Ordered, uniquely named set of :class:`Observable`."""

    def __init__(self, observables: Sequence[Observable]):
        self.observables: List[Observable] = []
        for o in observables:
            self.add(o)

    def add(self, obs: Observable) -> "ComparisonSpec":
        if obs.name in self.names():
            raise ValueError(f"Duplicate observable '{obs.name}' (use tag=).")
        self.observables.append(obs)
        return self

    def names(self) -> List[str]:
        return [o.name for o in self.observables]

    def __getitem__(self, name: str) -> Observable:
        for o in self.observables:
            if o.name == name:
                return o
        raise KeyError(f"No observable '{name}'; known: {self.names()}")

    def __iter__(self):
        return iter(self.observables)

    def __len__(self):
        return len(self.observables)

    @property
    def pins_unresolved(self) -> bool:
        return any(o.pins_unresolved for o in self.observables)

    def resolve_pins(self, lattice) -> "ComparisonSpec":
        """Copy with every pin selector (``"diagonal"``...) resolved on
        *lattice*; needed before generating or matching Serpent detectors."""
        return ComparisonSpec([o.resolve_pins(lattice) for o in self.observables])

    @classmethod
    def union(cls, *specs: "ComparisonSpec") -> "ComparisonSpec":
        """Union of several specs (e.g. over DRAGON variants). Observables
        with the same name must be identical; isotopes are merged."""
        merged: Dict[str, Observable] = {}
        for s in specs:
            for o in s:
                if o.name not in merged:
                    merged[o.name] = Observable(**_obs_kwargs(o))
                    continue
                m = merged[o.name]
                if (m.spatial, m.mesh_code, m.pins, m.statepoints) != \
                        (o.spatial, o.mesh_code, o.pins, o.statepoints):
                    raise ValueError(f"Observable '{o.name}' differs between specs.")
                m.isotopes = m.isotopes + [i for i in o.isotopes if i not in m.isotopes]
                m.serpent = tuple(dict.fromkeys(m.serpent + o.serpent))
                m.dragon = m.dragon or o.dragon
        return cls(list(merged.values()))

    # ── presets ────────────────────────────────────────────────────────────
    @classmethod
    def default(cls, isotopes: Sequence[str], fine_mesh: str = "295g",
                **kwargs) -> "ComparisonSpec":
        """The default observable set: :meth:`preset` ``"integrated"``."""
        return cls.preset("integrated", isotopes, fine_mesh=fine_mesh, **kwargs)

    @classmethod
    def preset(cls, name: str, isotopes: Sequence[str],
               fine_mesh: str = "295g",
               serpent_fine_pins: Optional[Union[str, Sequence[int]]] = None,
               fuel_isotopes_inventory: Optional[Sequence[str]] = None,
               global_totals: bool = True) -> "ComparisonSpec":
        """Named observable sets.

        Every preset starts with ``FUEL_1G``, the **fuel inventory** edition:
        DRAGON only, one group, homogenised over the fuel. Its ``MICR`` list
        is *fuel_isotopes_inventory* (default *isotopes*), so the COMPO holds
        their fuel-averaged number densities and 1-group cross sections
        (``CompoReader.densities``; × ``mixture_volume`` for the total). The
        keff of the comparison is read from this output, but keff does not
        depend on the list (it is a property of the flux calculation).

        ``"integrated"`` (the default, small Serpent outputs): macroscopic
        rates from ``void`` responses, fine meshes only where spatially
        integrated.

        - ``PIN_2G``: per-pin 2-group macroscopic fission rate (pin-power
          shape; ``dr -6 void`` per pin, plus ``-7`` with *global_totals*);
        - ``FUEL_{fine}``: fuel-integrated macroscopic capture and fission
          (``-2``, ``-6``, so absorption = capture + fission) and flux;
        - ``ALL_{fine}``: assembly flux spectrum (plus the assembly
          absorption with *global_totals*).

        *isotopes* are the DRAGON edit (``MICR``); every other nuclide is in
        ``*MAC*RES``, so the DRAGON macroscopic rates are complete.

        ``"pin_rates"``: ``FUEL_1G`` and ``PIN_2G`` with per-isotope rates.

        ``"pin_rates_and_spectra"`` (the GE14 DOM/DOM-C/VAN analysis):
        per-isotope rates **on every pin and on the fine mesh**, i.e. the
        DRAGON fine-mesh edition pin by pin (large ``det0.m``):

        - ``PIN_2G``: per-isotope pin rates, 2 groups;
        - ``PIN_{fine}``: per-isotope pin rates and pin flux on *fine_mesh*.
          DRAGON edits every pin; *serpent_fine_pins* limits the Serpent
          detectors to a list of ``pin_idx`` or to a diagonal (``"diagonal"``
          = the symmetry axis, ``"main"``, ``"anti"``), which is what the
          spectra along the diagonal need and keeps ``det0.m`` small;
        - ``FUEL_{fine}``: fuel-integrated rates, every partial reaction
          (:data:`physics_data.DEFAULT_DETECTOR_MTS`);
        - ``ALL_{fine}``: assembly spectrum.

        *global_totals* (default) adds the Serpent ``void`` responses of the
        whole-assembly normalisations (``fission_source``, ``fission``,
        ``absorption``): macroscopic fission and ν-fission on ``PIN_2G``
        (:data:`GLOBAL_PIN_MACRO`), absorption on ``ALL_{fine}``
        (:data:`GLOBAL_ALL_MACRO`). They cost one or two responses per
        detector.

        Use ``fine_mesh="26g"`` for a 2L scheme (its outputs are on 26 groups).
        """
        isotopes = list(isotopes)
        inventory = Observable(
            "FUEL", "1g", list(fuel_isotopes_inventory or isotopes), serpent=(),
            comment="Fuel inventory: 1 group, homogenised over the fuel (keff)")
        f = fine_mesh
        pin_macro = GLOBAL_PIN_MACRO if global_totals else ()
        all_obs = Observable(
            "ALL", f, isotopes,
            serpent=("flux", "rr") if global_totals else ("flux",),
            macro=GLOBAL_ALL_MACRO if global_totals else (),
            serpent_isotopes=() if global_totals else None,
            comment=f"{f}, homogenised over the assembly")
        if serpent_fine_pins is not None and name != "pin_rates_and_spectra":
            raise ValueError("serpent_fine_pins only applies to the 'pin_rates_and_spectra' "
                             "preset (its PIN_{fine} observable).")
        if name == "integrated":
            return cls([
                inventory,
                Observable("PIN", "2g", isotopes, serpent_isotopes=(),
                           macro=tuple(dict.fromkeys(("fission",) + pin_macro)),
                           comment="2 groups, per pin (macroscopic fission)"),
                Observable("FUEL", f, isotopes, macro=("capture", "fission"),
                           serpent=("rr", "flux"), serpent_isotopes=(),
                           comment=f"{f}, homogenised over the fuel"),
                all_obs,
            ])
        base = [inventory, Observable("PIN", "2g", isotopes, macro=pin_macro,
                                      comment="2 groups, per pin")]
        if name == "pin_rates":
            return cls(base)
        if name != "pin_rates_and_spectra":
            raise KeyError(f"Unknown preset '{name}' (integrated, pin_rates, "
                           "pin_rates_and_spectra).")
        return cls(base + [
            Observable("PIN", f, isotopes, serpent=("rr", "flux"), pins=serpent_fine_pins,
                       comment=f"{f}, per pin"),
            # One detector: every partial reaction is cheap here (as in the
            # GE14 references' det_assembly_295g).
            Observable("FUEL", f, isotopes, comment=f"{f}, homogenised over the fuel",
                       reactions={i: phys.DEFAULT_DETECTOR_MTS for i in isotopes}),
            all_obs,
        ])

    @classmethod
    def from_scheme(cls, scheme, macro: Sequence[str] = (),
                    serpent_isotopes: Optional[Sequence[str]] = None,
                    serpent: Optional[Mapping[str, Sequence[str]]] = None,
                    library_mesh: Sequence[float] = phys.G295_EV,
                    global_totals: bool = True) -> "ComparisonSpec":
        """Observables matching the DRAGON outputs of *scheme*
        (``CALCULATION_OUTPUTS``), keeping their names (``output=``), so that
        the Serpent detectors match what the DRAGON run edits.

        Args:
            macro: Macroscopic rates added to every PIN/FUEL observable.
            serpent_isotopes: Per-isotope Serpent responses (default: each
                output's isotopes; ``()`` for macroscopic only).
            serpent: ``{output_name: kinds}`` overrides of the Serpent kinds.
                Defaults: ``()`` for 1-group outputs (fuel inventory, keff), ``("flux",)`` for
                ``ALL``, ``("rr",)`` otherwise.
            global_totals: Add the Serpent ``void`` responses of the
                whole-assembly normalisations (as :meth:`preset`): fission
                and ν-fission on the coarsest PIN output, absorption on every
                ALL output. Serpent side only; DRAGON rebuilds them from its
                output (edited isotopes + ``*MAC*RES``).
        """
        out_mesh = scheme_output_mesh(scheme, library_mesh)
        serpent = dict(serpent or {})
        coarsest_pin = None
        if global_totals:
            pins_out = [(len(cond_to_bounds(o.energy_bounds, out_mesh)), o.name)
                        for o in scheme.outputs
                        if DRAGON_TO_SPATIAL.get(o.spatial_integration_mode) == "PIN"]
            coarsest_pin = min(pins_out)[1] if pins_out else None
        obs: List[Observable] = []
        for out in scheme.outputs:
            sp = DRAGON_TO_SPATIAL.get(out.spatial_integration_mode)
            if sp is None:
                raise ValueError(f"Output {out.name}: spatial mode "
                                 f"'{out.spatial_integration_mode}' has no Serpent equivalent.")
            bounds = cond_to_bounds(out.energy_bounds, out_mesh)
            kinds = serpent.get(out.name)
            if kinds is None:
                kinds = () if len(bounds) == 2 else (("flux",) if sp == "ALL" else ("rr",))
            m = () if sp == "ALL" else tuple(macro)
            s_iso = None if serpent_isotopes is None else list(serpent_isotopes)
            if global_totals and out.name == coarsest_pin:
                m = tuple(dict.fromkeys(m + GLOBAL_PIN_MACRO))
                kinds = tuple(dict.fromkeys(tuple(kinds) + ("rr",)))
            if global_totals and sp == "ALL":
                m = GLOBAL_ALL_MACRO
                kinds = tuple(dict.fromkeys(tuple(kinds) + ("rr",)))
                s_iso = []                       # ALL: macroscopic responses only
            kw = dict(serpent=tuple(kinds), statepoints=out.state_points, output=out.name,
                      comment=out.comment or "", macro=m, serpent_isotopes=s_iso)
            o = Observable(sp, bounds, list(out.isotopes), **kw)
            k = 2
            while o.name in [x.name for x in obs]:
                o = Observable(sp, bounds, list(out.isotopes), tag=f"_{k}", **kw)
                k += 1
            obs.append(o)
        return cls(obs)

    def dragon_output(self, observable: str) -> str:
        """DRAGON output name of *observable*."""
        return self[observable].dragon_name

    # ── DRAGON side ────────────────────────────────────────────────────────
    def to_dragon_outputs(self, scheme, skip_unavailable: bool = False,
                          library_mesh: Sequence[float] = phys.G295_EV) -> list:
        """:class:`CalculationOutput` list for *scheme*. Observables whose mesh
        the scheme cannot edit raise (all listed), or are skipped with a
        warning if *skip_unavailable*."""
        from ..DDModel.DragonCalculationScheme import CalculationOutput
        out_mesh = scheme_output_mesh(scheme, library_mesh)
        outputs, errors = [], []
        for o in self.observables:
            if not o.dragon:
                continue
            try:
                cond = dragon_cond(o.bounds_eV, out_mesh)
            except phys.MeshMismatchError as exc:
                errors.append(f"{o.name}: {exc}")
                continue
            outputs.append(CalculationOutput(
                name=o.dragon_name, comment=o.comment or o.name, isotopes=list(o.isotopes),
                spatial_integration_mode=SPATIAL_TO_DRAGON[o.spatial],
                energy_bounds=cond,
                state_points=o.statepoints if isinstance(o.statepoints, str)
                else {k: list(v) for k, v in o.statepoints.items()}))
        if errors:
            msg = (f"The scheme edits on a {len(out_mesh) - 1}-group mesh; these "
                   "observables cannot be produced:\n  " + "\n  ".join(errors))
            if not skip_unavailable:
                raise phys.MeshMismatchError(msg)
            warnings.warn(msg, stacklevel=2)
        return outputs

    def apply_to(self, scheme, **kwargs):
        """Copy of *scheme* whose outputs are this spec's (validated)."""
        s = scheme.copy()
        s.outputs = self.to_dragon_outputs(s, **kwargs)
        s.validate()
        return s

    # ── Serpent side ───────────────────────────────────────────────────────
    def _raw_requests(self) -> List[SerpentRequest]:
        if self.pins_unresolved:
            bad = [f"{o.name} (pins='{o.pins}')" for o in self if o.pins_unresolved]
            raise ValueError(f"Pin selectors of {bad} are not resolved: call "
                             "spec.resolve_pins(lattice) first.")
        return [SerpentRequest(o.name, k, o.spatial, o.bounds_eV,
                               [] if k == "flux" else o.micro_isotopes,
                               o.pins, None if o.reactions is None else dict(o.reactions),
                               () if k == "flux" else tuple(o.macro))
                for o in self.observables for k in o.serpent]

    def serpent_requests(self) -> List[SerpentRequest]:
        """Requests to generate, after removing those another request
        provides (finest mesh kept)."""
        raw = sorted(self._raw_requests(), key=lambda r: -len(r.bounds_eV))
        kept: List[SerpentRequest] = []
        for r in raw:
            if not any(q.provides(r) for q in kept):
                kept.append(r)
        order = {(r.observable, r.kind): i for i, r in enumerate(self._raw_requests())}
        return sorted(kept, key=lambda r: order[(r.observable, r.kind)])

    def serpent_source(self, observable: str, kind: str = "rr",
                       manifest: Optional[DetectorManifest] = None,
                       n_pins: Optional[int] = None,
                       pins: Optional[Sequence[int]] = None) -> str:
        """Serpent observable whose detectors give *observable*/*kind*: in
        the detectors this spec generates, or, with *manifest*, in an
        existing reference (e.g. ``PIN_2G`` from the legacy ``PIN_295G``).
        Pass the result as ``observable`` to ``compute_pin_rates``."""
        target = [r for r in self._raw_requests() if (r.observable, r.kind) == (observable, kind)]
        if not target:
            raise KeyError(f"'{observable}' has no Serpent {kind} in this spec.")
        if manifest is not None:
            src, msgs = _source_in(manifest, target[0], n_pins, pins)
            if src is None:
                raise KeyError(f"{observable}/{kind} not provided by the reference: "
                               + "; ".join(msgs))
            return src
        for q in self.serpent_requests():
            if q.provides(target[0]):
                return q.observable
        raise AssertionError("unreachable")   # every request is kept or provided

    def configure_serpent(self, model, library: str, response_temperature: Optional[float] = None,
                          xs_suffix: Optional[str] = None) -> DetectorManifest:
        """Add energy grids, response materials and detectors to a built
        :class:`Serpent2Model`, and return (and attach) their manifest."""
        from ..InterfaceToDD.Serpent2_exports import S2_Detector
        if not getattr(model, "_built", False):
            raise RuntimeError("Build the Serpent2Model (model.build()) before configuring detectors.")
        if self.pins_unresolved:              # e.g. serpent_fine_pins="diagonal"
            from .lattice_maps import PinLattice
            return self.resolve_pins(PinLattice.from_assembly(model.assembly)).configure_serpent(
                model, library, response_temperature=response_temperature, xs_suffix=xs_suffix)
        if response_temperature is None:
            response_temperature = getattr(model.assembly, "default_fuel_temperature", 900.0)
        pins_map = dict(sorted(model._pin_universe_map.items()))
        entries: List[ManifestEntry] = []
        for req in self.serpent_requests():
            grid = req.grid_name
            if grid is not None:
                model.ensure_energy_grid(grid, list(req.bounds_eV * 1e-6))
            responses: List[Tuple[int, str]] = []
            for iso in req.isotopes:
                mts = (req.reactions or {}).get(iso) or phys.serpent_detector_mts(library, iso)
                responses += [(int(mt), iso) for mt in mts]
                model.ensure_response_material(iso, response_temperature, xs_suffix)
            responses += [(mt, phys.VOID) for mt in req.macro_mts]   # collision material
            if req.spatial == "PIN":
                pins = req.pins or list(pins_map)
                missing = sorted(set(pins) - set(pins_map))
                if missing:
                    raise KeyError(f"{req.observable}: no fuel pin universe for pins {missing}.")
                domains = [(p, pins_map[p].get_fuel_material_names()) for p in pins]
            elif req.spatial == "FUEL":
                domains = [(None, [m for u in pins_map.values() for m in u.get_fuel_material_names()])]
            else:
                domains = [(None, [])]
            suffix = "RR" if req.kind == "rr" else "FLX"
            for pin, mats in domains:
                name = f"{req.observable}_{suffix}" + ("" if pin is None else f"_p{pin}")
                if any(d.name == name for d in model.detectors):
                    raise ValueError(f"Detector '{name}' already exists in the model.")
                det = S2_Detector(name, energy_grid_name=grid, domain_materials=list(mats),
                                  responses=list(responses),
                                  detector_type=-4 if mats else None)
                model.detectors.append(det)
                entries.append(ManifestEntry(
                    detector=name, observable=req.observable, kind=req.kind,
                    spatial=req.spatial, pin_idx=pin,
                    pins=[pin] if pin is not None else (list(pins_map) if req.spatial == "FUEL" else []),
                    materials=list(mats),
                    responses=[Response(mt, mat, None if mat == phys.VOID else mat)
                               for mt, mat in responses],
                    grid=grid, n_groups=len(req.bounds_eV) - 1,
                    bounds_eV=[float(x) for x in req.bounds_eV],
                    detector_type=["-4"] if mats else [], source="exporter"))
        manifest = DetectorManifest(entries, source={"exporter": "starterDD",
                                                     "library": library,
                                                     "spec": self.names()})
        manifest.validate()
        model.detector_manifest_entries = list(getattr(model, "detector_manifest_entries", [])) + entries
        return manifest

    # ── coverage of a Serpent reference ────────────────────────────────────
    def coverage(self, manifest: DetectorManifest, n_pins: Optional[int] = None,
                 pins: Optional[Sequence[int]] = None) -> Dict[str, List[str]]:
        """Problems per ``observable/kind`` for which *manifest* (an existing
        or generated reference) cannot provide the spec. Empty dict = covered.

        A Serpent observable covers a request when it has the same spatial
        and kind, a mesh that refines the request's mesh, a response for
        every requested isotope and, for ``PIN``, a detector for every pin
        (*pins*, or the request's pins, or all pins seen in the manifest).
        """
        problems: Dict[str, List[str]] = {}
        for req in self._raw_requests():
            src, msgs = _source_in(manifest, req, n_pins, pins)
            if src is None:
                problems[f"{req.observable}/{req.kind}"] = msgs
        return problems


def _source_in(manifest: DetectorManifest, req: SerpentRequest, n_pins: Optional[int],
               pins: Optional[Sequence[int]]) -> Tuple[Optional[str], List[str]]:
    """First manifest observable that fully provides *req* (finest mesh
    first), or ``(None, reasons)``."""
    cands = [e for e in manifest
             if e.spatial == req.spatial and e.kind == req.kind and e.bounds_eV]
    if not cands:
        return None, [f"no {req.spatial} {req.kind} detector"]
    fine = [e for e in cands if refines(e.bounds_eV, req.bounds_eV)]
    if not fine:
        meshes = sorted({mesh_code(e.bounds_eV) for e in cands})
        legacy = any(phys.is_legacy_mesh(e.bounds_eV) for e in cands)
        return None, [f"no detector mesh refines {mesh_code(req.bounds_eV)} "
                      f"(available: {meshes})"
                      + (" - legacy 295g grid, compare with mesh='common'" if legacy else "")]
    by_obs: Dict[str, List[ManifestEntry]] = {}
    for e in sorted(fine, key=lambda e: -(e.n_groups or 0)):
        by_obs.setdefault(e.observable, []).append(e)
    msgs: List[str] = []
    for obs, es in by_obs.items():
        m = []
        if req.isotopes:
            have = {r.nuclide for e in es for r in e.responses}
            lack = [i for i in req.isotopes if i not in have]
            if lack:
                m.append(f"{obs}: no response for {lack}")
        if req.macro_mts:
            have_v = {r.mt for e in es for r in e.responses if r.material == phys.VOID}
            lack_v = [mt for mt in req.macro_mts if mt not in have_v]
            if lack_v:
                m.append(f"{obs}: no 'void' response for MT {lack_v} ({', '.join(req.macro)})")
        if req.spatial == "PIN":
            want = req.pins or pins or (list(range(1, n_pins + 1)) if n_pins else None)
            have_p = {e.pin_idx for e in es}
            want = want or sorted(p for p in have_p if p is not None)
            lack_p = sorted(set(want) - have_p)
            if lack_p:
                m.append(f"{obs}: no detector for pins {lack_p}")
        if not m:
            return obs, []
        msgs += m
    return None, msgs


def _obs_kwargs(o: Observable) -> dict:
    return dict(spatial=o.spatial, mesh=list(o.bounds_eV), isotopes=list(o.isotopes),
                serpent=tuple(o.serpent), dragon=o.dragon, pins=o.pins,
                statepoints=o.statepoints, tag=o.tag, reactions=o.reactions,
                comment=o.comment, macro=tuple(o.macro),
                serpent_isotopes=o.serpent_isotopes, output=o.output)


# ─────────────────────────────────────────────────────────────────────────────
# DRAGON outputs <-> observables (including legacy names)
# ─────────────────────────────────────────────────────────────────────────────


def output_observable(compo, output: str, sp=None) -> str:
    """Observable name of a COMPO output, from its name (naming convention or
    :data:`LEGACY_OUTPUT_SPATIAL`) and its ENERGY record."""
    spatial = output.split("_")[0] if output.split("_")[0] in SPATIAL_TO_DRAGON \
        else LEGACY_OUTPUT_SPATIAL.get(output)
    if spatial is None:
        raise KeyError(f"Cannot tell the spatial mode of output '{output}'.")
    if sp is None:
        sp = compo.statepoints(output)[0].calc
    return f"{spatial}_{mesh_code(compo.energy_bounds(output, sp))}"


def find_output(compo, observable: str) -> str:
    """COMPO output holding *observable*: same name, or a legacy output whose
    spatial mode and energy mesh give that observable."""
    outs = compo.outputs()
    if observable in outs:
        return observable
    hits = []
    for o in outs:
        try:
            if output_observable(compo, o) == observable:
                hits.append(o)
        except KeyError:
            continue
    if len(hits) != 1:
        raise KeyError(f"{'No' if not hits else 'Several'} output(s) for '{observable}' "
                       f"among {outs} {hits or ''}".rstrip())
    return hits[0]


__all__ = ["Observable", "ComparisonSpec", "SerpentRequest", "resolve_mesh", "mesh_code",
           "refines", "scheme_output_mesh", "dragon_cond", "cond_to_bounds", "output_observable",
           "find_output", "LEGACY_OUTPUT_SPATIAL", "SPATIAL_TO_DRAGON", "PIN_SELECTORS",
           "GLOBAL_PIN_MACRO", "GLOBAL_ALL_MACRO"]
