# -*- coding: utf-8 -*-
"""Serpent2 output readers and detector manifests.

Three layers:

1. **Raw readers** (no knowledge of starterDD conventions):
   :func:`read_det` (``*_det0.m``), :func:`read_res` (``*_res.m``) and
   :func:`parse_serpent_input` (``.serp`` input: ``mat``, ``ene``, ``det``
   cards, ``include`` files).

2. **Detector manifest** (:class:`DetectorManifest`): for every detector, the
   observable it contributes to, its kind (reaction rates or flux), the pin it
   covers, its domain, its ordered ``(MT, nuclide)`` responses and its energy
   grid. Post-processing only uses the manifest, never detector names.

   A manifest is either written by the starterDD exporter (P4) or **inferred
   from an existing input** with :meth:`DetectorManifest.infer`. Inference
   uses the *content* of each ``det`` card (``de``/``dm``/``du``/``dc``/``dr``
   and the response materials' ``mat`` cards), so legacy references with
   arbitrary detector names (``det_295g_UOX16_1``, ``det_assembly_295g``,
   ``flux_295g``...) can be reused without re-running Serpent. The inferred
   manifest can be frozen to JSON next to the run (:meth:`DetectorManifest.write`)
   and edited by hand if needed.

3. **Run access** (:class:`SerpentRun`): input + det0 + res + manifest, with
   consistency checks (response counts, energy bounds) and accessors that
   return group-wise arrays in ascending energy.

The legacy functions of ``temp_examples/serpent_reader.py``
(:func:`serpent_simple_det_reader`, :func:`serpent_simple_serp_reader`,
:func:`serpent_simple_res_reader`) are kept with identical outputs so that the
existing analysis scripts can switch imports before being refactored.
"""

from __future__ import annotations

import hashlib
import json
import re
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from . import physics_data as phys
from ._io import existing_variant, read_text, strip_gz

PathLike = Union[str, Path]

MANIFEST_FORMAT = "starterDD-serpent-detectors"
MANIFEST_VERSION = 1

# ═════════════════════════════════════════════════════════════════════════════
# 1. Raw readers
# ═════════════════════════════════════════════════════════════════════════════

# ── det0.m ──────────────────────────────────────────────────────────────────

#: Column layout of a 12-column detector array (Serpent 2).
DET_COLUMNS_12 = ("index", "energy", "universe", "cell", "material", "lattice",
                  "reaction", "z", "y", "x", "value", "rel_err")
#: 13-column layout written by Serpent versions with time binning.
DET_COLUMNS_13 = ("index", "time", "energy", "universe", "cell", "material",
                  "lattice", "reaction", "z", "y", "x", "value", "rel_err")

#: Suffixes of auxiliary arrays written after a detector (energy grid,
#: mesh coordinates...).
_DET_AUX_SUFFIXES = ("COORD", "PHI", "THETA", "E", "X", "Y", "Z", "R", "T")

_DET_BLOCK_RE = re.compile(r"^DET(\S+?)\s*=\s*\[(.*?)\];", re.M | re.S)


@dataclass
class DetectorResult:
    """One detector from a ``det0.m`` file.

    Attributes:
        name: Detector name (as in the ``det`` card).
        raw: The raw array, one row per bin (12 or 13 columns).
        aux: Auxiliary arrays (``"E"`` for the energy grid, ``"X"``...).
    """

    name: str
    raw: np.ndarray
    aux: Dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def columns(self) -> Tuple[str, ...]:
        n = self.raw.shape[1]
        if n == 12:
            return DET_COLUMNS_12
        if n == 13:
            return DET_COLUMNS_13
        raise ValueError(f"Detector '{self.name}': unsupported number of "
                         f"columns ({n}) in det0.m.")

    def _col(self, name: str) -> np.ndarray:
        return self.raw[:, self.columns.index(name)]

    @property
    def bin_axes(self) -> Tuple[str, ...]:
        """Names of the bin axes, in the order of :meth:`values`."""
        return tuple(c for c in self.columns
                     if c not in ("index", "value", "rel_err"))

    @property
    def shape(self) -> Tuple[int, ...]:
        """Number of bins along each axis of :attr:`bin_axes`."""
        return tuple(int(self._col(a).max()) for a in self.bin_axes)

    def _scatter(self, column: str) -> np.ndarray:
        out = np.full(self.shape, np.nan)
        idx = tuple(self._col(a).astype(int) - 1 for a in self.bin_axes)
        out[idx] = self._col(column)
        return out

    def values(self) -> np.ndarray:
        """Tally means as an array indexed by :attr:`bin_axes` (1 bin per
        unused axis). Energy is ascending (Serpent order)."""
        return self._scatter("value")

    def rel_err(self) -> np.ndarray:
        """Relative 1-sigma errors, same layout as :meth:`values`."""
        return self._scatter("rel_err")

    @property
    def n_reactions(self) -> int:
        return int(self._col("reaction").max())

    @property
    def n_energy(self) -> int:
        return int(self._col("energy").max())

    def response_by_energy(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return ``(values, rel_err)`` shaped ``(n_reactions, n_energy)``.

        Raises ``ValueError`` if any bin axis other than reaction / energy
        has more than one bin (mesh or multi-material detectors need the
        generic :meth:`values`).
        """
        axes = self.bin_axes
        shape = self.shape
        extra = {a: n for a, n in zip(axes, shape)
                 if a not in ("reaction", "energy") and n > 1}
        if extra:
            raise ValueError(
                f"Detector '{self.name}' has extra bins {extra}; use values()."
            )
        ie, ir = axes.index("energy"), axes.index("reaction")
        v = np.moveaxis(self.values(), (ir, ie), (0, 1)).reshape(
            shape[ir], shape[ie])
        e = np.moveaxis(self.rel_err(), (ir, ie), (0, 1)).reshape(
            shape[ir], shape[ie])
        return v, e

    @property
    def energy_bounds_eV(self) -> Optional[np.ndarray]:
        """Ascending bounds (eV) from the ``E`` auxiliary array, or ``None``
        when the detector has no energy binning."""
        e = self.aux.get("E")
        if e is None:
            return None
        e = np.atleast_2d(e)
        return np.append(e[:, 0], e[-1, 1]) * 1.0e6


def _parse_numeric_block(body: str) -> np.ndarray:
    lines = [ln for ln in body.strip().splitlines() if ln.strip()]
    if not lines:
        return np.zeros((0, 0))
    ncol = len(lines[0].split())
    flat = np.array(" ".join(lines).split(), dtype=float)
    if flat.size % ncol:
        raise ValueError("Ragged detector block in det0.m")
    return flat.reshape(-1, ncol)


def read_det(path: PathLike) -> Dict[str, DetectorResult]:
    """Read a Serpent2 ``*_det<N>.m`` file.

    Returns:
        ``{detector_name: DetectorResult}``. Auxiliary arrays
        (``DET<name>E`` energy grids, mesh coordinates) are attached to their
        detector instead of being returned as separate entries.
    """
    text = read_text(path)
    blocks: Dict[str, np.ndarray] = {}
    for m in _DET_BLOCK_RE.finditer(text):
        blocks[m.group(1)] = _parse_numeric_block(m.group(2))

    names = set(blocks)
    aux_of: Dict[str, Tuple[str, str]] = {}
    for name in names:
        for suf in _DET_AUX_SUFFIXES:
            base = name[: -len(suf)]
            if name.endswith(suf) and base in names and base != name:
                # Prefer the longest base (e.g. "detE" vs "det" + "E").
                prev = aux_of.get(name)
                if prev is None or len(base) > len(prev[0]):
                    aux_of[name] = (base, suf)
    # An array that is somebody's aux AND has its own aux is a detector
    # (e.g. detectors "a" and "aE" both with energy grids): keep it.
    has_aux = {b for b, _ in aux_of.values()}
    for name in list(aux_of):
        if name in has_aux and blocks[name].shape[1] in (12, 13):
            del aux_of[name]

    out: Dict[str, DetectorResult] = {}
    for name, arr in blocks.items():
        if name not in aux_of:
            out[name] = DetectorResult(name, arr)
    for name, (base, suf) in aux_of.items():
        if base in out:
            out[base].aux[suf] = blocks[name]
    return out


# ── res.m ───────────────────────────────────────────────────────────────────

_RES_LINE_RE = re.compile(
    r"^(?P<key>[A-Z][A-Z0-9_]*)\s+\(idx,\s*(?P<shape>\[[^\]]*\]|\d+)\)\s*=\s*"
    r"(?P<val>.*?)\s*;\s*$"
)


def _parse_res_value(raw: str):
    raw = raw.strip()
    if raw.startswith("["):
        body = raw[1:-1].strip() if raw.endswith("]") else raw[1:]
        return np.array(body.split(), dtype=float) if body else np.zeros(0)
    if raw[:1] in ("'", '"'):
        return raw[1:-1] if len(raw) >= 2 else ""
    try:
        v = float(raw)
    except ValueError:
        return raw
    # Same rule as the legacy reader: integral values become int.
    return int(v) if v.is_integer() else v


@dataclass
class ResFile:
    """Content of a ``*_res.m`` file.

    ``values[key]`` holds the value of the **last** step (the legacy reader
    behaviour); ``steps[key]`` holds one value per step for keys written at
    several burnup / branch steps.
    """

    values: Dict[str, object]
    steps: Dict[str, List[object]]

    def __getitem__(self, key: str):
        return self.values[key]

    def __contains__(self, key: str) -> bool:
        return key in self.values

    def get(self, key: str, default=None):
        return self.values.get(key, default)

    def keff(self, estimator: str = "IMP_KEFF") -> Tuple[float, float]:
        """Return ``(k, relative 1-sigma)`` for *estimator*
        (``IMP_KEFF``, ``ANA_KEFF``, ``ABS_KEFF``, ``COL_KEFF``,
        ``INF_KINF``...)."""
        v = np.asarray(self.values[estimator], dtype=float)
        return float(v[0]), float(v[1])


def read_res(path: PathLike) -> ResFile:
    """Read a Serpent2 ``*_res.m`` file (all steps)."""
    values: Dict[str, object] = {}
    steps: Dict[str, List[object]] = {}
    for line in read_text(path).splitlines():
        if not line or line[0] == "%":
            continue
        m = _RES_LINE_RE.match(line)
        if not m:
            continue
        val = _parse_res_value(m.group("val"))
        key = m.group("key")
        values[key] = val
        steps.setdefault(key, []).append(val)
    return ResFile(values, steps)


# ── .serp input ─────────────────────────────────────────────────────────────

#: Top-level Serpent2 card keywords (anything else continues the current card).
SERPENT_CARDS = frozenset({
    "set", "mat", "mix", "therm", "thermstoch", "det", "ene", "surf", "cell",
    "pin", "lat", "include", "plot", "mesh", "src", "dep", "div", "trans",
    "utrans", "strans", "ftrans", "ltrans", "particle", "nest", "solid",
    "umsh", "ifc", "datamesh", "sens", "wwgen", "wwin", "branch", "coef",
    "casematrix", "hisv", "rep", "sample", "tme", "dtrans", "disperse",
    "pbed", "usym", "fun", "phb",
})

#: Sub-options of the ``det`` card.
DET_OPTIONS = frozenset({
    "de", "dt", "dm", "du", "dc", "dl", "dr", "dx", "dy", "dz", "dn", "dh",
    "dumsh", "dfl", "dtl", "dir", "ds", "dhis", "di", "dv", "dphb", "dfet",
    "da", "dtt",
})

#: Arity of ``mat`` card options.
_MAT_OPTIONS = {"tmp": 1, "tms": 1, "tft": 2, "rgb": 3, "vol": 1, "mass": 1,
                "burn": 1, "fix": 2, "moder": 2, "cdom": 1}

_TOKEN_RE = re.compile(r'"[^"]*"|\'[^\']*\'|%.*$|[^\s%]+')


@dataclass
class Token:
    text: str
    file: str
    line: int


def _tokenize(path: Path, follow_includes: bool = True,
              _seen: Optional[set] = None) -> List[Token]:
    _seen = set() if _seen is None else _seen
    path = path.resolve()
    if path in _seen:
        raise ValueError(f"Recursive include of {path}")
    _seen.add(path)
    text = read_text(path)
    # Block comments, keeping line numbers.
    text = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"),
                  text, flags=re.S)
    tokens: List[Token] = []
    pending_include = False
    for lineno, line in enumerate(text.splitlines(), start=1):
        for m in _TOKEN_RE.finditer(line):
            tok = m.group(0)
            if tok.startswith("%"):
                break
            if pending_include:
                pending_include = False
                inc = Path(tok.strip("\"'"))
                if not inc.is_absolute():
                    inc = path.parent / inc
                tokens.extend(_tokenize(inc, follow_includes, _seen))
                continue
            if tok == "include" and follow_includes:
                pending_include = True
                continue
            tokens.append(Token(tok, str(path), lineno))
    return tokens


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


@dataclass
class SerpentMaterial:
    """A ``mat`` card.

    Attributes:
        name: Material name.
        density: Density as written (>0 atomic b^-1cm^-1, <0 mass g/cm3),
            ``None`` for ``sum`` (sum of the nuclide fractions).
        options: ``mat`` options (``tmp``, ``moder``, ``vol``...), values as
            token lists.
        nuclides: Ordered ``{nuclide token: fraction}`` as written.
    """

    name: str
    density: Optional[float]
    options: Dict[str, List[str]] = field(default_factory=dict)
    nuclides: Dict[str, float] = field(default_factory=dict)

    @property
    def temperature(self) -> Optional[float]:
        for key in ("tmp", "tms"):
            if key in self.options:
                return float(self.options[key][0])
        return None

    def composition(self) -> Dict[str, float]:
        """``{DRAGON-style nuclide name: fraction}``, summing duplicate
        nuclides (different library suffixes); unparseable tokens kept."""
        out: Dict[str, float] = {}
        for tok, frac in self.nuclides.items():
            key = phys.nuclide_name(tok) or tok
            out[key] = out.get(key, 0.0) + frac
        return out

    @property
    def single_nuclide(self) -> Optional[str]:
        """Nuclide name if the material holds exactly one nuclide (response
        materials), else ``None``."""
        comp = self.composition()
        if len(comp) == 1:
            name = next(iter(comp))
            return name if phys.nuclide_name(name) else None
        return None


@dataclass
class EnergyGridDef:
    """An ``ene`` card. ``bounds_eV`` is ``None`` for predefined structures
    (type 3), whose bounds are only known from the det0.m output."""

    name: str
    grid_type: int
    args: List[str]
    bounds_eV: Optional[np.ndarray]

    @classmethod
    def from_tokens(cls, name: str, args: List[str]) -> "EnergyGridDef":
        gtype = int(args[0])
        rest = args[1:]
        bounds = None
        if gtype == 1:
            bounds = np.array(rest, dtype=float) * 1.0e6
        elif gtype == 2:
            n, emin, emax = int(rest[0]), float(rest[1]), float(rest[2])
            bounds = phys.equal_lethargy_mesh(n, emin * 1e6, emax * 1e6)
        elif gtype == 4:
            n, emin, emax = int(rest[0]), float(rest[1]), float(rest[2])
            bounds = phys.uniform_mesh(n, emin * 1e6, emax * 1e6)
        return cls(name, gtype, args, None if bounds is None
                   else phys.to_ascending(bounds))


@dataclass
class DetectorDef:
    """A ``det`` card as written in the input."""

    name: str
    grid: Optional[str] = None
    dtype: List[str] = field(default_factory=list)
    materials: List[str] = field(default_factory=list)
    universes: List[str] = field(default_factory=list)
    cells: List[str] = field(default_factory=list)
    lattices: List[str] = field(default_factory=list)
    responses: List[Tuple[int, str]] = field(default_factory=list)
    other: Dict[str, List[List[str]]] = field(default_factory=dict)
    file: str = ""
    line: int = 0

    @property
    def is_mesh(self) -> bool:
        return any(k in self.other for k in ("dx", "dy", "dz", "dn", "dh",
                                             "dumsh"))


@dataclass
class SerpentInput:
    """Parsed Serpent2 input (the parts needed for post-processing)."""

    path: Path
    materials: Dict[str, SerpentMaterial]
    grids: Dict[str, EnergyGridDef]
    detectors: Dict[str, DetectorDef]
    settings: Dict[str, List[str]]
    #: ``pin`` cards: universe -> material names (radii and ``fill`` dropped).
    pins: Dict[str, List[str]] = field(default_factory=dict)
    #: ``lat`` cards: lattice name -> all tokens after the name.
    lattices: Dict[str, List[str]] = field(default_factory=dict)

    @property
    def title(self) -> Optional[str]:
        t = self.settings.get("title")
        return " ".join(x.strip("\"'") for x in t) if t else None

    def universe_counts(self) -> Dict[str, int]:
        """Number of lattice positions filled by each ``pin`` universe
        (summed over all ``lat`` cards; nested lattices are not expanded)."""
        counts: Dict[str, int] = {}
        for toks in self.lattices.values():
            for t in toks:
                if t in self.pins:
                    counts[t] = counts.get(t, 0) + 1
        return counts

    def domain_multiplicity(self, materials: Sequence[str]) -> Optional[int]:
        """Number of lattice positions covered by a detector domain made of
        *materials*: the positions of every ``pin`` universe containing one
        of them. ``None`` if the materials are in no lattice universe."""
        mats = set(materials)
        counts = self.universe_counts()
        unis = [u for u, ms in self.pins.items() if mats & set(ms)]
        if not unis or not any(u in counts for u in unis):
            return None
        return sum(counts.get(u, 0) for u in unis)


def _split_cards(tokens: List[Token]) -> List[List[Token]]:
    cards: List[List[Token]] = []
    for tok in tokens:
        if tok.text in SERPENT_CARDS:
            cards.append([tok])
        elif cards:
            cards[-1].append(tok)
    return cards


def _parse_mat(card: List[Token]) -> SerpentMaterial:
    words = [t.text for t in card[1:]]
    # Density "sum" = total of the nuclide fractions (stored as None).
    density = None if words[1].lower() == "sum" else float(words[1])
    mat = SerpentMaterial(name=words[0], density=density)
    i = 2
    while i < len(words):
        w = words[i]
        if w in _MAT_OPTIONS:
            n = _MAT_OPTIONS[w]
            mat.options[w] = words[i + 1:i + 1 + n]
            i += 1 + n
        elif i + 1 < len(words) and _is_number(words[i + 1]):
            mat.nuclides[w] = mat.nuclides.get(w, 0.0) + float(words[i + 1])
            i += 2
        else:
            # Unknown option without arity: store and move on.
            mat.options.setdefault(w, [])
            i += 1
    return mat


def _parse_det(card: List[Token]) -> DetectorDef:
    head = card[1]
    det = DetectorDef(name=head.text, file=head.file, line=head.line)
    groups: List[Tuple[str, List[str]]] = []
    for tok in card[2:]:
        if tok.text in DET_OPTIONS:
            groups.append((tok.text, []))
        elif groups:
            groups[-1][1].append(tok.text)
        else:
            raise ValueError(f"{tok.file}:{tok.line}: unexpected token "
                             f"'{tok.text}' in det card '{det.name}'")
    for opt, args in groups:
        if opt == "de":
            det.grid = args[0]
        elif opt == "dt":
            det.dtype = args
        elif opt == "dm":
            det.materials.extend(args)
        elif opt == "du":
            det.universes.extend(args)
        elif opt == "dc":
            det.cells.extend(args)
        elif opt == "dl":
            det.lattices.extend(args)
        elif opt == "dr":
            if len(args) % 2:
                raise ValueError(f"{head.file}:{head.line}: odd number of "
                                 f"dr arguments in det '{det.name}'")
            for k in range(0, len(args), 2):
                det.responses.append((int(args[k]), args[k + 1]))
        else:
            det.other.setdefault(opt, []).append(args)
    return det


def parse_serpent_input(path: PathLike, follow_includes: bool = True) -> SerpentInput:
    """Parse the ``mat``, ``ene``, ``det`` and ``set`` cards of a Serpent2
    input (following ``include`` files relative to the including file)."""
    path = Path(path)
    cards = _split_cards(_tokenize(path, follow_includes))
    materials: Dict[str, SerpentMaterial] = {}
    grids: Dict[str, EnergyGridDef] = {}
    detectors: Dict[str, DetectorDef] = {}
    settings: Dict[str, List[str]] = {}
    pins: Dict[str, List[str]] = {}
    lattices: Dict[str, List[str]] = {}
    for card in cards:
        kw = card[0].text
        if len(card) < 2:
            continue
        if kw == "mat":
            mat = _parse_mat(card)
            materials[mat.name] = mat
        elif kw == "ene":
            g = EnergyGridDef.from_tokens(card[1].text,
                                          [t.text for t in card[2:]])
            grids[g.name] = g
        elif kw == "det":
            det = _parse_det(card)
            if det.name in detectors:
                raise ValueError(f"Detector '{det.name}' defined twice "
                                 f"({card[1].file}:{card[1].line}).")
            detectors[det.name] = det
        elif kw == "set":
            settings[card[1].text] = [t.text for t in card[2:]]
        elif kw == "pin":
            # pin NAME  mat1 r1  mat2 r2 ... matN   (or "fill UNI r")
            words = [t.text for t in card[2:]]
            pin_mats = [w for w in words if not _is_number(w) and w != "fill"]
            pins[card[1].text] = pin_mats
        elif kw == "lat":
            lattices[card[1].text] = [t.text for t in card[2:]]
    return SerpentInput(path, materials, grids, detectors, settings, pins, lattices)


# ═════════════════════════════════════════════════════════════════════════════
# 2. Detector manifest
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class PinResolver:
    """Maps detector domains (materials / universes / cells) to pin indices.

    Resolution order for one domain name: explicit maps first
    (``material_to_pin``, ``universe_to_pin``, ``cell_to_pin``), then the
    regular expressions (first capture group = pin index).

    The defaults match the starterDD exporter: fuel materials are named
    ``<MAT>_zone_<z>_pin_<p>`` and pin universes ``<MAT>_<p>``.
    """

    material_patterns: Sequence[str] = (r"_pin_(\d+)$",)
    universe_patterns: Sequence[str] = (r"_(\d+)$",)
    cell_patterns: Sequence[str] = ()
    material_to_pin: Mapping[str, int] = field(default_factory=dict)
    universe_to_pin: Mapping[str, int] = field(default_factory=dict)
    cell_to_pin: Mapping[str, int] = field(default_factory=dict)

    def _resolve(self, name: str, mapping: Mapping[str, int],
                 patterns: Sequence[str]) -> Optional[int]:
        if name in mapping:
            return int(mapping[name])
        for pat in patterns:
            m = re.search(pat, name)
            if m:
                return int(m.group(1))
        return None

    def pins_of(self, det: DetectorDef) -> List[Optional[int]]:
        """Pin index of every domain entry (``None`` if unresolved)."""
        out = [self._resolve(n, self.material_to_pin, self.material_patterns)
               for n in det.materials]
        out += [self._resolve(n, self.universe_to_pin, self.universe_patterns)
                for n in det.universes]
        out += [self._resolve(n, self.cell_to_pin, self.cell_patterns)
                for n in det.cells]
        return out

    def to_dict(self) -> dict:
        return {
            "material_patterns": list(self.material_patterns),
            "universe_patterns": list(self.universe_patterns),
            "cell_patterns": list(self.cell_patterns),
            "material_to_pin": dict(self.material_to_pin),
            "universe_to_pin": dict(self.universe_to_pin),
            "cell_to_pin": dict(self.cell_to_pin),
        }

    @classmethod
    def from_dict(cls, d: Mapping) -> "PinResolver":
        return cls(**{k: d[k] for k in (
            "material_patterns", "universe_patterns", "cell_patterns",
            "material_to_pin", "universe_to_pin", "cell_to_pin") if k in d})


@dataclass
class Response:
    """One ``dr`` response: MT number, response material and the nuclide it
    represents (``None`` for macroscopic / unresolved responses)."""

    mt: int
    material: str
    nuclide: Optional[str]


@dataclass
class ManifestEntry:
    """Semantic description of one Serpent detector.

    Attributes:
        detector: Detector name in the Serpent input / det0.m.
        observable: Observable name, ``{SPATIAL}_{MESH}`` (e.g. ``PIN_295G``).
        kind: ``"rr"`` (reaction rates, has responses) or ``"flux"``.
        spatial: ``PIN`` | ``FUEL`` | ``GROUP`` | ``ALL`` | ``MAT`` | ``MESH``.
        pin_idx: Pin index for ``PIN`` detectors, else ``None``.
        pins: All pin indices covered (sorted).
        materials, universes, cells: The detector domain as written.
        responses: Ordered responses; position ``k`` is reaction bin ``k+1``.
        grid: Energy grid name (``None`` = no energy binning).
        n_groups: Number of energy groups (``None`` if unknown before reading
            det0.m, i.e. predefined grids).
        bounds_eV: Ascending energy bounds, if known.
        detector_type: ``dt`` arguments as written.
        source: ``"inferred"``, ``"exporter"`` or ``"override"``.
    """

    detector: str
    observable: str
    kind: str
    spatial: str
    pin_idx: Optional[int] = None
    pins: List[int] = field(default_factory=list)
    materials: List[str] = field(default_factory=list)
    universes: List[str] = field(default_factory=list)
    cells: List[str] = field(default_factory=list)
    responses: List[Response] = field(default_factory=list)
    grid: Optional[str] = None
    n_groups: Optional[int] = None
    bounds_eV: Optional[List[float]] = None
    detector_type: List[str] = field(default_factory=list)
    source: str = "inferred"

    @property
    def key(self) -> Tuple[str, str, Optional[int]]:
        return (self.observable, self.kind, self.pin_idx)

    def response_index(self, mt: int, nuclide: str) -> Optional[int]:
        """0-based reaction bin of ``(mt, nuclide)`` or ``None``."""
        for k, r in enumerate(self.responses):
            if r.mt == mt and r.nuclide == nuclide:
                return k
        return None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["responses"] = [asdict(r) for r in self.responses]
        return d

    @classmethod
    def from_dict(cls, d: Mapping) -> "ManifestEntry":
        d = dict(d)
        d["responses"] = [Response(**r) for r in d.get("responses", [])]
        return cls(**d)


class ManifestConflictError(ValueError):
    """Two detectors map to the same ``(observable, kind, pin)``."""


@dataclass
class DetectorManifest:
    """Ordered collection of :class:`ManifestEntry`, with JSON round trip."""

    entries: List[ManifestEntry]
    source: Dict[str, object] = field(default_factory=dict)
    pin_resolver: Optional[PinResolver] = None

    # ── queries ────────────────────────────────────────────────────────────
    def __iter__(self):
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def by_detector(self, name: str) -> ManifestEntry:
        for e in self.entries:
            if e.detector == name:
                return e
        raise KeyError(f"No detector '{name}' in manifest.")

    def observables(self) -> List[str]:
        seen: List[str] = []
        for e in self.entries:
            if e.observable not in seen:
                seen.append(e.observable)
        return seen

    def select(self, observable: Optional[str] = None,
               kind: Optional[str] = None,
               pin: Optional[int] = None) -> List[ManifestEntry]:
        return [e for e in self.entries
                if (observable is None or e.observable == observable)
                and (kind is None or e.kind == kind)
                and (pin is None or e.pin_idx == pin)]

    def find(self, observable: str, kind: str,
             pin: Optional[int] = None) -> ManifestEntry:
        hits = [e for e in self.entries if e.key == (observable, kind, pin)]
        if not hits:
            raise KeyError(
                f"No detector for observable={observable!r}, kind={kind!r}, "
                f"pin={pin!r}. Available observables: {self.observables()}"
            )
        return hits[0]

    def pins(self, observable: str, kind: str) -> List[int]:
        return sorted(e.pin_idx for e in self.entries
                      if e.observable == observable and e.kind == kind
                      and e.pin_idx is not None)

    def validate(self) -> None:
        """Raise :class:`ManifestConflictError` on duplicate keys."""
        seen: Dict[Tuple, str] = {}
        clashes = []
        for e in self.entries:
            if e.key in seen:
                clashes.append(f"{seen[e.key]} / {e.detector} -> {e.key}")
            else:
                seen[e.key] = e.detector
        if clashes:
            raise ManifestConflictError(
                "Several detectors map to the same (observable, kind, pin):\n  "
                + "\n  ".join(clashes)
                + "\nUse `exclude=[...]` or `overrides={det: {...}}` "
                  "to disambiguate."
            )

    def summary(self) -> str:
        """One line per (observable, kind) with the number of detectors."""
        rows: Dict[Tuple[str, str], List[ManifestEntry]] = {}
        for e in self.entries:
            rows.setdefault((e.observable, e.kind), []).append(e)
        lines = []
        for (obs, kind), es in rows.items():
            pins = {p for e in es for p in e.pins}   # distinct pins covered
            nres = sorted({len(e.responses) for e in es})
            ng = sorted({e.n_groups for e in es if e.n_groups is not None})
            lines.append(
                f"{obs:<12} {kind:<4} detectors={len(es):<4} "
                f"pins={len(pins):<4} responses={nres} groups={ng}"
            )
        return "\n".join(lines)

    # ── inference from an existing input ───────────────────────────────────
    @classmethod
    def infer(cls, serp: Union[PathLike, SerpentInput],
              pin_resolver: Optional[PinResolver] = None,
              overrides: Optional[Mapping[str, Mapping]] = None,
              exclude: Iterable[str] = (),
              validate: bool = True) -> "DetectorManifest":
        """Build a manifest from the ``det`` cards of an existing input.

        Classification uses only the card content:

        - ``kind`` = ``"rr"`` if the card has ``dr`` responses, else ``"flux"``;
        - ``spatial`` = ``MESH`` (``dx``/``dy``/... present), ``ALL`` (no
          domain), ``PIN`` (every domain entry resolves to the same pin),
          ``FUEL`` (domain spans every pin seen in the input), ``GROUP``
          (several pins, not all) or ``MAT`` (domain not resolvable to pins);
        - mesh code from the energy grid (``2G``, ``295G``, ``C{n}G``; ``1G``
          without ``de``);
        - each response's nuclide from its response material's ``mat`` card,
          falling back to the material name itself (``U235``).

        Args:
            serp: Input path or a parsed :class:`SerpentInput`.
            pin_resolver: How domain names map to pins (defaults to the
                starterDD naming).
            overrides: ``{detector_name: {field: value}}`` applied after
                inference (e.g. ``{"det_x": {"observable": "FUEL_26G"}}``).
            exclude: Detector names to leave out.
            validate: Check that no two detectors share a key.
        """
        inp = serp if isinstance(serp, SerpentInput) else parse_serpent_input(serp)
        resolver = pin_resolver or PinResolver()
        overrides = dict(overrides or {})
        exclude = set(exclude)
        unknown = (set(overrides) | exclude) - set(inp.detectors)
        if unknown:
            raise KeyError(f"overrides/exclude name unknown detectors: "
                           f"{sorted(unknown)}")

        pins_by_det = {n: resolver.pins_of(d) for n, d in inp.detectors.items()
                       if n not in exclude}
        all_pins = sorted({p for ps in pins_by_det.values() for p in ps
                           if p is not None})

        entries: List[ManifestEntry] = []
        for name, det in inp.detectors.items():
            if name in exclude:
                continue
            pins = pins_by_det[name]
            resolved = sorted({p for p in pins if p is not None})
            has_domain = bool(det.materials or det.universes or det.cells
                              or det.lattices)
            if det.is_mesh:
                spatial = "MESH"
            elif not has_domain:
                spatial = "ALL"
            elif not pins or any(p is None for p in pins) or det.lattices:
                spatial = "MAT"
            elif len(resolved) == 1:
                spatial = "PIN"
            elif resolved == all_pins:
                spatial = "FUEL"
            else:
                spatial = "GROUP"

            bounds = None
            n_groups: Optional[int] = 1
            if det.grid is not None:
                g = inp.grids.get(det.grid)
                if g is None:
                    raise KeyError(f"{det.file}:{det.line}: det '{name}' uses "
                                   f"undefined energy grid '{det.grid}'.")
                bounds = g.bounds_eV
                n_groups = None if bounds is None else len(bounds) - 1
            if bounds is not None:
                mesh = phys.mesh_code(bounds)
            elif det.grid is None:
                mesh = "1G"
            else:
                mesh = det.grid.upper()

            responses = []
            for mt, rmat in det.responses:
                m = inp.materials.get(rmat)
                nuc = m.single_nuclide if m is not None else None
                if nuc is None:
                    nuc = phys.nuclide_name(rmat)
                responses.append(Response(mt, rmat, nuc))

            entry = ManifestEntry(
                detector=name,
                observable=f"{spatial}_{mesh}",
                kind="rr" if det.responses else "flux",
                spatial=spatial,
                pin_idx=resolved[0] if spatial == "PIN" else None,
                pins=resolved,
                materials=list(det.materials),
                universes=list(det.universes),
                cells=list(det.cells),
                responses=responses,
                grid=det.grid,
                n_groups=n_groups,
                bounds_eV=None if bounds is None else bounds.tolist(),
                detector_type=list(det.dtype),
            )
            if name in overrides:
                for k, v in overrides[name].items():
                    if not hasattr(entry, k):
                        raise AttributeError(f"Unknown manifest field '{k}'")
                    setattr(entry, k, v)
                entry.source = "override"
            entries.append(entry)

        manifest = cls(entries, source={
            "input": inp.path.name,           # the manifest sits next to it; no machine path
            "input_sha256": _sha256(inp.path),
            "title": inp.title,
            "inferred": True,
        }, pin_resolver=resolver)
        if validate:
            manifest.validate()
        return manifest

    # ── JSON ───────────────────────────────────────────────────────────────
    def to_dict(self) -> dict:
        return {
            "format": MANIFEST_FORMAT,
            "version": MANIFEST_VERSION,
            "source": self.source,
            "pin_resolver": None if self.pin_resolver is None
            else self.pin_resolver.to_dict(),
            "detectors": [e.to_dict() for e in self.entries],
        }

    def write(self, path: PathLike) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=1))
        return path

    @classmethod
    def from_dict(cls, d: Mapping) -> "DetectorManifest":
        if d.get("format") != MANIFEST_FORMAT:
            raise ValueError("Not a starterDD detector manifest.")
        if int(d.get("version", 0)) > MANIFEST_VERSION:
            raise ValueError(f"Manifest version {d['version']} is newer than "
                             f"supported ({MANIFEST_VERSION}).")
        pr = d.get("pin_resolver")
        return cls([ManifestEntry.from_dict(e) for e in d["detectors"]],
                   source=dict(d.get("source", {})),
                   pin_resolver=None if pr is None else PinResolver.from_dict(pr))

    @classmethod
    def read(cls, path: PathLike) -> "DetectorManifest":
        return cls.from_dict(json.loads(Path(path).read_text()))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def default_manifest_path(input_path: PathLike) -> Path:
    """``<input>.detectors.json`` next to the input (``.gz`` suffix of a
    compressed input is dropped)."""
    p = strip_gz(input_path)
    return p.with_name(p.name + ".detectors.json")


# ═════════════════════════════════════════════════════════════════════════════
# 3. Run access
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class Tally:
    """Group-wise data of one detector, ascending energy.

    Attributes:
        entry: The manifest entry.
        values: ``(n_responses, n_groups)`` for ``rr``, ``(1, n_groups)``
            for ``flux``.
        rel_err: Relative 1-sigma errors, same shape.
        bounds_eV: Ascending energy bounds (``None`` without energy binning).
    """

    entry: ManifestEntry
    values: np.ndarray
    rel_err: np.ndarray
    bounds_eV: Optional[np.ndarray]

    @property
    def abs_err(self) -> np.ndarray:
        return self.values * self.rel_err

    def response(self, mt: int, nuclide: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """``(values, rel_err)`` of one response, ``None`` if not scored."""
        k = self.entry.response_index(mt, nuclide)
        if k is None:
            return None
        return self.values[k], self.rel_err[k]

    def condense(self, coarse_bounds: Sequence[float],
                 error_mode: str = "quadrature") -> "Tally":
        """Condense onto a coarser mesh (must be a sub-mesh, see
        :func:`physics_data.condensation_indices`)."""
        if self.bounds_eV is None:
            raise ValueError(f"Detector '{self.entry.detector}' has no "
                             "energy binning.")
        idx = phys.condensation_indices(self.bounds_eV, coarse_bounds)
        v = phys.condense(self.values, idx)
        a = phys.condense_abs_error(self.abs_err, idx, mode=error_mode)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(v != 0, a / np.abs(v), 0.0)
        return Tally(self.entry, v, r, phys.to_ascending(coarse_bounds))


class SerpentRun:
    """A finished Serpent2 run: input, detectors, results and manifest.

    Typical legacy use::

        run = SerpentRun("GE14_DOM_00_endfb8r1.serp")    # infers the manifest
        run.manifest.write(default_manifest_path(run.input_path))  # freeze it
        k, dk = run.keff()
        t = run.tally("PIN_295G", "rr", pin=12).condense(phys.G2_EV)

    The manifest is taken, in order, from the *manifest* argument, the
    ``<input>.detectors.json`` file next to the input, or inferred from the
    input. The det0.m and res.m files are read lazily.
    """

    def __init__(self, input_path: PathLike,
                 manifest: Union[None, PathLike, DetectorManifest] = None,
                 det_path: Optional[PathLike] = None,
                 res_path: Optional[PathLike] = None,
                 pin_resolver: Optional[PinResolver] = None,
                 overrides: Optional[Mapping[str, Mapping]] = None,
                 exclude: Iterable[str] = ()):
        self.input_path = Path(input_path)
        base = strip_gz(self.input_path)   # outputs are named after the plain input
        self.det_path = Path(det_path) if det_path else \
            existing_variant(base.with_name(base.name + "_det0.m"))
        self.res_path = Path(res_path) if res_path else \
            existing_variant(base.with_name(base.name + "_res.m"))
        self._input: Optional[SerpentInput] = None
        self._det: Optional[Dict[str, DetectorResult]] = None
        self._res: Optional[ResFile] = None

        if isinstance(manifest, DetectorManifest):
            self.manifest = manifest
        elif manifest is not None:
            self.manifest = DetectorManifest.read(manifest)
        elif default_manifest_path(self.input_path).exists() and not (
                pin_resolver or overrides or exclude):
            self.manifest = DetectorManifest.read(
                default_manifest_path(self.input_path))
        else:
            self.manifest = DetectorManifest.infer(
                self.input, pin_resolver=pin_resolver, overrides=overrides,
                exclude=exclude)
        self._check_input_hash()

    # ── lazy members ───────────────────────────────────────────────────────
    @property
    def input(self) -> SerpentInput:
        if self._input is None:
            self._input = parse_serpent_input(self.input_path)
        return self._input

    @property
    def detectors(self) -> Dict[str, DetectorResult]:
        if self._det is None:
            self._det = read_det(self.det_path)
        return self._det

    @property
    def results(self) -> ResFile:
        if self._res is None:
            self._res = read_res(self.res_path)
        return self._res

    def _check_input_hash(self) -> None:
        ref = self.manifest.source.get("input_sha256")
        if ref and self.input_path.exists() and _sha256(self.input_path) != ref:
            warnings.warn(
                f"{self.input_path.name} changed since its manifest was "
                "written; check that the manifest still matches the outputs.",
                stacklevel=3,
            )

    # ── accessors ──────────────────────────────────────────────────────────
    def keff(self, estimator: str = "IMP_KEFF") -> Tuple[float, float]:
        """``(k, relative 1-sigma)``."""
        return self.results.keff(estimator)

    def detector(self, name: str) -> DetectorResult:
        try:
            return self.detectors[name]
        except KeyError:
            raise KeyError(f"Detector '{name}' not found in "
                           f"{self.det_path.name}") from None

    def tally_of(self, entry: ManifestEntry) -> Tally:
        """Read the detector of *entry*, checking it against the manifest."""
        det = self.detector(entry.detector)
        values, rel = det.response_by_energy()
        n_expected = max(len(entry.responses), 1)
        if values.shape[0] != n_expected:
            raise ValueError(
                f"Detector '{entry.detector}': det0.m has {values.shape[0]} "
                f"reaction bins but the manifest lists {n_expected} responses."
            )
        bounds = det.energy_bounds_eV
        if bounds is not None and entry.bounds_eV is not None and \
                not phys.bounds_match(bounds, entry.bounds_eV):
            raise ValueError(
                f"Detector '{entry.detector}': energy bounds in det0.m differ "
                f"from grid '{entry.grid}' in the input/manifest."
            )
        if bounds is None and entry.bounds_eV is not None:
            bounds = np.asarray(entry.bounds_eV)
        if bounds is not None and len(bounds) - 1 != values.shape[1]:
            raise ValueError(f"Detector '{entry.detector}': {values.shape[1]} "
                             f"energy bins but {len(bounds) - 1} groups.")
        return Tally(entry, values, rel, bounds)

    def tally(self, observable: str, kind: str,
              pin: Optional[int] = None) -> Tally:
        return self.tally_of(self.manifest.find(observable, kind, pin))

    def check(self) -> List[str]:
        """Read every manifest detector and return a list of problems
        (empty if the manifest and outputs are consistent)."""
        problems = []
        for e in self.manifest:
            try:
                self.tally_of(e)
            except (KeyError, ValueError) as exc:
                problems.append(str(exc))
        return problems


# ═════════════════════════════════════════════════════════════════════════════
# Legacy API (temp_examples/serpent_reader.py, B. Daré)
# ═════════════════════════════════════════════════════════════════════════════


def serpent_simple_det_reader(path: str) -> Dict[str, np.ndarray]:
    """Legacy: ``{name: raw array}`` including ``<name>E`` energy arrays."""
    text = read_text(path)
    return {m.group(1): _parse_numeric_block(m.group(2))
            for m in _DET_BLOCK_RE.finditer(text)}


def serpent_simple_serp_reader(path: str) -> Dict[str, dict]:
    """Legacy: ``{name: {"energy": n_groups, (mt, mat): reaction_bin}}``.

    Reaction bins are 1-based, as in column 7 of the det0.m arrays. Unlike
    the original, multiple ``dr`` per line and ``include`` files are handled.
    """
    inp = parse_serpent_input(path)
    out: Dict[str, dict] = {}
    for name, det in inp.detectors.items():
        d: dict = {}
        if det.grid is None:
            d["energy"] = 1
        else:
            g = inp.grids.get(det.grid)
            # The original reader stored int(grid_name[:-1]) ("295g" -> 295).
            if g is not None and g.bounds_eV is not None:
                d["energy"] = len(g.bounds_eV) - 1
            else:
                d["energy"] = int(re.sub(r"\D", "", det.grid) or 1)
        for k, (mt, mat) in enumerate(det.responses, start=1):
            d[(mt, mat)] = k
        out[name] = d
    return out


def serpent_simple_res_reader(path: str) -> Dict[str, object]:
    """Legacy: ``{key: value}`` from res.m (last step)."""
    return dict(read_res(path).values)


__all__ = [
    "DetectorResult", "read_det", "ResFile", "read_res",
    "SerpentMaterial", "EnergyGridDef", "DetectorDef", "SerpentInput",
    "parse_serpent_input", "PinResolver", "Response", "ManifestEntry",
    "ManifestConflictError", "DetectorManifest", "default_manifest_path",
    "Tally", "SerpentRun",
    "serpent_simple_det_reader", "serpent_simple_serp_reader",
    "serpent_simple_res_reader",
]
