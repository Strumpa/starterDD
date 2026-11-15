# -*- coding: utf-8 -*-
"""Statepoint-aware reader for DRAGON MULTICOMPO objects.

A MULTICOMPO export contains one top-level directory per EDI/COMPO output
(``EDIR_2G``, ``EDIHOM_COND``...). Each directory holds:

- ``GLOBAL``: the global parameters (``PARKEY``, ``PARFMT``, ``NVALUE`` and
  one ``pval%08d`` record of sorted values per parameter);
- ``MIXTURES``: a list with, for each mixture, a ``CALCULATIONS`` list (one
  item per elementary calculation = statepoint) and a ``TREE`` table that maps
  parameter-value tuples to calculation indices.

Tree decoding follows DRAGON5 ``NCRCAL`` / ``COMARB``: nodes are 1-based, the
root is node 1, the children of node ``i`` are ``DEBARB(i) .. DEBARB(i+1)-1``
and ``ARBVAL(j)`` is the (1-based) index of node ``j``'s value in the
parameter of its depth. At depth ``NPAR`` (a leaf), the calculation index is
``DEBARB(i+1)``. ``DEBARB``/``ARBVAL`` are allocated with spare room: only the
first ``NVP(1)+1`` / ``NVP(1)`` entries are meaningful.

Conventions of all accessors: energies in eV, **ascending**; group-wise arrays
shaped ``(n_mixtures, n_groups)`` (or ``(n_groups,)`` for a single mixture).
"""

from __future__ import annotations

import hashlib
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .pylcm import LcmChars, LcmH5, Pylcm, h5py_available

PathLike = Union[str, Path]
Node = Union[Pylcm, LcmH5]
StatepointSel = Union[None, int, Mapping[str, object]]

CACHE_SUFFIX = ".h5"
CACHE_FORMAT_VERSION = 1


# ═════════════════════════════════════════════════════════════════════════════
# HDF5 cache
# ═════════════════════════════════════════════════════════════════════════════


def _file_signature(path: Path) -> Dict[str, object]:
    st = path.stat()
    return {"source_size": int(st.st_size), "source_mtime_ns": int(st.st_mtime_ns)}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def default_cache_path(compo_path: PathLike) -> Path:
    """``<compo>.h5`` next to the ASCII COMPO."""
    p = Path(compo_path)
    return p.with_name(p.name + CACHE_SUFFIX)


def cache_is_valid(compo_path: PathLike, cache_path: PathLike,
                   check_hash: bool = True) -> bool:
    """True if *cache_path* was built from the current *compo_path*.

    Size must match. If the modification time also matches, the cache is
    accepted; otherwise (e.g. the files were copied) the SHA-256 of the
    source is compared when *check_hash* is true.
    """
    compo_path, cache_path = Path(compo_path), Path(cache_path)
    if not cache_path.exists():
        return False
    import h5py
    try:
        with h5py.File(cache_path, "r") as h5:
            a = dict(h5.attrs)
    except OSError:
        return False
    if int(a.get("cache_format_version", -1)) != CACHE_FORMAT_VERSION:
        return False
    sig = _file_signature(compo_path)
    if int(a.get("source_size", -1)) != sig["source_size"]:
        return False
    if int(a.get("source_mtime_ns", -1)) == sig["source_mtime_ns"]:
        return True
    return bool(check_hash) and a.get("source_sha256") == _sha256(compo_path)


def build_cache(compo_path: PathLike, cache_path: Optional[PathLike] = None,
                only: Optional[Iterable[str]] = None) -> Path:
    """Parse an ASCII COMPO and write its HDF5 cache."""
    compo_path = Path(compo_path)
    cache_path = Path(cache_path) if cache_path else default_cache_path(compo_path)
    lcm = Pylcm.from_ascii(compo_path, name="MULTICOMPO", only=only)
    attrs = {
        "cache_format_version": CACHE_FORMAT_VERSION,
        "source_path": str(compo_path.resolve()),
        "source_sha256": _sha256(compo_path),
        **_file_signature(compo_path),
    }
    if only is not None:
        attrs["only"] = ",".join(only)
    return lcm.to_hdf5(cache_path, attrs=attrs)


def open_compo(path: PathLike, cache: Union[bool, PathLike] = True,
               only: Optional[Iterable[str]] = None) -> Node:
    """Open a MULTICOMPO, from its HDF5 cache when possible.

    Args:
        path: ASCII COMPO export, or an ``.h5`` file written by
            :meth:`Pylcm.to_hdf5` / :func:`build_cache`.
        cache: ``True`` = use/create ``<path>.h5``; a path = use/create that
            cache file; ``False`` = parse the ASCII file into memory.
        only: Top-level outputs to load when parsing into memory
            (ignored when a full cache is used).

    Returns:
        A lazy :class:`LcmH5` view (cached) or an in-memory :class:`Pylcm`.
    """
    path = Path(path)
    if path.suffix == CACHE_SUFFIX:
        return LcmH5(path)
    if cache and not h5py_available():
        warnings.warn("h5py not installed: parsing the ASCII COMPO in memory "
                      "(install with pip install -e '.[post]').", stacklevel=2)
        cache = False
    if not cache:
        return Pylcm.from_ascii(path, name="MULTICOMPO", only=only)
    cache_path = default_cache_path(path) if cache is True else Path(cache)
    if not cache_is_valid(path, cache_path):
        try:
            build_cache(path, cache_path)
        except OSError as exc:  # read-only results directory, full disk...
            warnings.warn(f"Could not write COMPO cache {cache_path} ({exc}); "
                          "parsing in memory.", stacklevel=2)
            return Pylcm.from_ascii(path, name="MULTICOMPO", only=only)
    return LcmH5(cache_path)


# ═════════════════════════════════════════════════════════════════════════════
# Parameter tree
# ═════════════════════════════════════════════════════════════════════════════


def decode_tree(nvp: int, debarb: Sequence[int], arbval: Sequence[int],
                n_par: int) -> List[Tuple[Tuple[int, ...], int]]:
    """Enumerate the leaves of a COMPO parameter tree.

    Args:
        nvp: Number of nodes (``NVP(1)``).
        debarb: ``DEBARB`` record (at least ``nvp+1`` entries).
        arbval: ``ARBVAL`` record (at least ``nvp`` entries).
        n_par: Number of parameters (tree depth).

    Returns:
        ``[(value_indices, calc_index), ...]`` in tree order, where
        ``value_indices`` are 1-based indices into each parameter's sorted
        ``pval`` values and ``calc_index`` is the 1-based CALCULATIONS index.
    """
    deb = np.atleast_1d(np.asarray(debarb, dtype=int))[: nvp + 1]
    arb = np.atleast_1d(np.asarray(arbval, dtype=int))[:nvp]
    if deb.size < nvp + 1 or arb.size < nvp:
        raise ValueError(f"Parameter tree records too short for NVP={nvp}.")
    out: List[Tuple[Tuple[int, ...], int]] = []
    # Iterative DFS: (node, depth, tuple)
    stack = [(1, 0, ())]
    while stack:
        node, depth, tup = stack.pop()
        if depth == n_par:
            out.append((tup, int(deb[node])))   # DEBARB(node+1), 1-based
            continue
        first, last = int(deb[node - 1]), int(deb[node]) - 1
        if not (1 <= first <= nvp + 1) or last > nvp:
            raise ValueError(f"Corrupt parameter tree at node {node}.")
        # Push in reverse so that children are visited in order.
        for child in range(last, first - 1, -1):
            stack.append((child, depth + 1, tup + (int(arb[child - 1]),)))
    return out


@dataclass(frozen=True)
class Statepoint:
    """One elementary calculation of a COMPO output.

    Attributes:
        calc: 1-based index in ``MIXTURES/i/CALCULATIONS``.
        values: ``{parameter: value}``.
        indices: ``{parameter: 1-based index in the sorted pval list}``.
        label: Optional label from a study manifest (``"00"``, ``"40"``...).
    """

    calc: int
    values: Dict[str, object]
    indices: Dict[str, int] = field(default_factory=dict)
    label: Optional[str] = None


class StatepointError(LookupError):
    """No (or more than one) calculation matches a statepoint selection."""


# ═════════════════════════════════════════════════════════════════════════════
# Reader
# ═════════════════════════════════════════════════════════════════════════════


def _to_list(x) -> list:
    if x is None:
        return []
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


def _names(chars, width: int) -> List[str]:
    if isinstance(chars, LcmChars):
        return [n for n in chars.chunks(width) if n]
    return str(chars).split()


class CompoReader:
    """Read statepoints and isotope-wise data from a MULTICOMPO.

    Args:
        source: ASCII COMPO path, ``.h5`` cache path, or an already opened
            :class:`Pylcm` / :class:`LcmH5`.
        cache: See :func:`open_compo`.
        labels: Optional ``{label: {parameter: value}}`` so that statepoints
            can be selected by label (``reader.keff("EDIR_2G", sp="40")``).
        rtol: Relative tolerance when matching REAL parameter values.

    Statepoint selection (argument ``sp`` of the accessors):

    - ``None``: the output must contain a single calculation;
    - ``int``: a CALCULATIONS index (1-based);
    - ``str``: a label from *labels*;
    - mapping ``{parameter: value}``: matched with tolerance; parameters
      that are omitted must be unambiguous.
    """

    def __init__(self, source: Union[PathLike, Node],
                 cache: Union[bool, PathLike] = True,
                 labels: Optional[Mapping[str, Mapping[str, object]]] = None,
                 rtol: float = 1.0e-5):
        if isinstance(source, (Pylcm, LcmH5)):
            self.root: Node = source
            self.path: Optional[Path] = None
        else:
            self.path = Path(source)
            self.root = open_compo(self.path, cache=cache)
        self.labels = dict(labels or {})
        self.rtol = rtol
        self._sp_cache: Dict[Tuple[str, int], List[Statepoint]] = {}
        self._par_cache: Dict[str, Dict[str, list]] = {}

    # ── context manager ────────────────────────────────────────────────────
    def close(self) -> None:
        if isinstance(self.root, LcmH5):
            self.root.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ── raw access ─────────────────────────────────────────────────────────
    def get(self, path: str, default=KeyError):
        return self.root.get(path, default)

    def _keys(self, path: str) -> list:
        node = self.root.get(path)
        return node.keys() if isinstance(node, LcmH5) else list(node.keys())

    # ── outputs and parameters ─────────────────────────────────────────────
    def outputs(self) -> List[str]:
        """Top-level directories that are COMPO outputs (have ``MIXTURES``)."""
        out = []
        for k in self._keys(""):
            node = self.root.get(str(k))
            if isinstance(node, (Pylcm, LcmH5)) and \
                    self.root.get(f"{k}/MIXTURES", None) is not None:
                out.append(str(k))
        return out

    def _check_output(self, output: str) -> None:
        if self.root.get(f"{output}/MIXTURES", None) is None:
            raise KeyError(f"'{output}' is not a COMPO output. "
                           f"Available: {self.outputs()}")

    def comment(self, output: str) -> Optional[str]:
        c = self.root.get(f"{output}/COMMENT", None)
        return None if c is None else str(c)

    def parameters(self, output: str) -> Dict[str, list]:
        """``{parameter: sorted values}`` from ``GLOBAL`` (PARKEY order).

        REAL values are floats, INTEGER ints, STRING stripped strings.
        """
        if output in self._par_cache:
            return self._par_cache[output]
        self._check_output(output)
        g = f"{output}/GLOBAL"
        if self.root.get(f"{output}/LOCAL", None) is not None:
            raise NotImplementedError(
                f"'{output}' has local parameters (LOCAL directory); only "
                "global parameters are supported.")
        parkey = self.root.get(f"{g}/PARKEY", None)
        if parkey is None:
            self._par_cache[output] = {}
            return {}
        keys = _names(parkey, 12)
        fmts = _names(self.root.get(f"{g}/PARFMT"), 8)
        nval = _to_list(self.root.get(f"{g}/NVALUE"))
        params: Dict[str, list] = {}
        for i, (key, fmt) in enumerate(zip(keys, fmts), start=1):
            raw = self.root.get(f"{g}/pval{i:08d}", None)
            if raw is None:
                vals: list = []
            elif fmt == "STRING":
                vals = raw.chunks(12) if isinstance(raw, LcmChars) else str(raw).split()
            elif fmt == "INTEGER":
                vals = [int(v) for v in _to_list(raw)]
            else:
                vals = [float(v) for v in _to_list(raw)]
            if i - 1 < len(nval) and len(vals) != int(nval[i - 1]):
                vals = vals[: int(nval[i - 1])]
            params[key] = vals
        self._par_cache[output] = params
        return params

    def n_mixtures(self, output: str) -> int:
        self._check_output(output)
        return len(self._keys(f"{output}/MIXTURES"))

    def n_calculations(self, output: str, mixture: int = 1) -> int:
        ncals = self.root.get(f"{output}/MIXTURES/{mixture}/TREE/NCALS", None)
        if ncals is not None:
            return int(np.atleast_1d(ncals)[0])
        return len(self._keys(f"{output}/MIXTURES/{mixture}/CALCULATIONS"))

    # ── statepoints ────────────────────────────────────────────────────────
    def statepoints(self, output: str, mixture: int = 1) -> List[Statepoint]:
        """All statepoints of *output* (decoded from ``MIXTURES/m/TREE``),
        sorted by calculation index."""
        ck = (output, mixture)
        if ck in self._sp_cache:
            return self._sp_cache[ck]
        params = self.parameters(output)
        names = list(params)
        tree = f"{output}/MIXTURES/{mixture}/TREE"
        ncals = self.n_calculations(output, mixture)
        if not names:
            sps = [Statepoint(c, {}, {}) for c in range(1, ncals + 1)]
        else:
            nvp = int(np.atleast_1d(self.root.get(f"{tree}/NVP"))[0])
            leaves = decode_tree(nvp, self.root.get(f"{tree}/DEBARB"),
                                 self.root.get(f"{tree}/ARBVAL"), len(names))
            sps = []
            for tup, calc in leaves:
                vals = {n: params[n][i - 1] for n, i in zip(names, tup)}
                sps.append(Statepoint(calc, vals, dict(zip(names, tup)),
                                      self._label_of(vals)))
            sps.sort(key=lambda s: s.calc)
            calcs = [s.calc for s in sps]
            if sorted(set(calcs)) != list(range(1, ncals + 1)) or len(calcs) != ncals:
                raise ValueError(
                    f"{output}: parameter tree gives calculations {calcs} but "
                    f"NCALS={ncals}.")
        self._sp_cache[ck] = sps
        return sps

    def _match(self, value, target) -> bool:
        if isinstance(value, str) or isinstance(target, str):
            return str(value).strip() == str(target).strip()
        a, b = float(value), float(target)
        return abs(a - b) <= self.rtol * max(abs(a), abs(b), 1e-300)

    def _label_of(self, values: Mapping[str, object]) -> Optional[str]:
        for label, sel in self.labels.items():
            if all(k in values and self._match(values[k], v) for k, v in sel.items()):
                return label
        return None

    def find_statepoints(self, output: str, mixture: int = 1,
                         **params) -> List[Statepoint]:
        """Statepoints whose parameters match *params* (with tolerance)."""
        known = self.parameters(output)
        unknown = set(params) - set(known)
        if unknown:
            raise KeyError(f"{output}: unknown parameter(s) {sorted(unknown)}; "
                           f"available {list(known)}")
        return [s for s in self.statepoints(output, mixture)
                if all(self._match(s.values[k], v) for k, v in params.items())]

    def calc_index(self, output: str, sp: StatepointSel = None,
                   mixture: int = 1) -> int:
        """Resolve a statepoint selection to a CALCULATIONS index."""
        if isinstance(sp, (int, np.integer)) and not isinstance(sp, bool):
            n = self.n_calculations(output, mixture)
            if not 1 <= int(sp) <= n:
                raise StatepointError(f"{output}: calculation {sp} out of "
                                      f"range 1..{n}")
            return int(sp)
        if isinstance(sp, str):
            if sp not in self.labels:
                raise StatepointError(f"Unknown statepoint label '{sp}'. "
                                      f"Known: {list(self.labels)}")
            sp = self.labels[sp]
        hits = self.find_statepoints(output, mixture, **dict(sp or {}))
        if len(hits) != 1:
            avail = [s.values for s in self.statepoints(output, mixture)]
            what = "No" if not hits else f"{len(hits)}"
            raise StatepointError(
                f"{output}: {what} calculation(s) match {dict(sp or {})}. "
                f"Statepoints: {avail[:10]}{' ...' if len(avail) > 10 else ''}"
            )
        return hits[0].calc

    def _calc_path(self, output: str, mixture: int, sp: StatepointSel) -> str:
        c = self.calc_index(output, sp, mixture)
        return f"{output}/MIXTURES/{mixture}/CALCULATIONS/{c}"

    def _mixtures(self, output: str, mixtures) -> List[int]:
        if mixtures is None:
            return list(range(1, self.n_mixtures(output) + 1))
        if isinstance(mixtures, (int, np.integer)):
            return [int(mixtures)]
        return [int(m) for m in mixtures]

    # ── scalar results ─────────────────────────────────────────────────────
    def keff(self, output: str, sp: StatepointSel = None, mixture: int = 1) -> float:
        return float(self.get(f"{self._calc_path(output, mixture, sp)}/K-EFFECTIVE"))

    def kinf(self, output: str, sp: StatepointSel = None, mixture: int = 1) -> float:
        return float(self.get(f"{self._calc_path(output, mixture, sp)}/K-INFINITY"))

    def energy_bounds(self, output: str, sp: StatepointSel = None,
                      mixture: int = 1) -> np.ndarray:
        """Energy bounds (eV, ascending) from the ``ENERGY`` record."""
        e = np.atleast_1d(np.asarray(
            self.get(f"{self._calc_path(output, mixture, sp)}/ENERGY"), float))
        return e[::-1].copy() if e.size > 1 and e[0] > e[-1] else e

    def n_groups(self, output: str, sp: StatepointSel = None) -> int:
        return len(self.energy_bounds(output, sp)) - 1

    # ── isotopes ───────────────────────────────────────────────────────────
    def isotopes(self, output: str, sp: StatepointSel = None,
                 mixture: int = 1) -> List[str]:
        """Isotope names in ``ISOTOPESLIST`` order (``ISOTOPESUSED``)."""
        used = self.get(f"{self._calc_path(output, mixture, sp)}/ISOTOPESUSED")
        return _names(used, 12)

    def _iso_index(self, output: str, isotope: str, sp, mixture: int) -> int:
        isos = self.isotopes(output, sp, mixture)
        try:
            return isos.index(isotope) + 1
        except ValueError:
            raise KeyError(f"{output} mixture {mixture}: no isotope "
                           f"'{isotope}' (have {isos})") from None

    def densities(self, output: str, sp: StatepointSel = None,
                  mixture: int = 1) -> Dict[str, float]:
        """``{isotope: number density (1e24/cm3)}``."""
        p = self._calc_path(output, mixture, sp)
        d = np.atleast_1d(self.get(f"{p}/ISOTOPESDENS"))
        return dict(zip(self.isotopes(output, sp, mixture), map(float, d)))

    def isotope_volumes(self, output: str, sp: StatepointSel = None,
                        mixture: int = 1) -> Dict[str, float]:
        """``{isotope: volume (cm3)}`` (``ISOTOPESVOL``)."""
        p = self._calc_path(output, mixture, sp)
        v = np.atleast_1d(self.get(f"{p}/ISOTOPESVOL"))
        return dict(zip(self.isotopes(output, sp, mixture), map(float, v)))

    def mixture_volume(self, output: str, sp: StatepointSel = None,
                       mixture: int = 1) -> float:
        return float(self.get(f"{self._calc_path(output, mixture, sp)}/MIXTURESVOL"))

    def isotope_record(self, output: str, isotope: str, key: str,
                       sp: StatepointSel = None, mixture: int = 1,
                       default=KeyError):
        """Raw record of an isotope directory (no reordering)."""
        p = self._calc_path(output, mixture, sp)
        i = self._iso_index(output, isotope, sp, mixture)
        return self.get(f"{p}/ISOTOPESLIST/{i}/{key}", default)

    def _groupwise(self, output: str, isotope: str, key: str,
                   sp: StatepointSel, mixtures, missing: str) -> Optional[np.ndarray]:
        mix = self._mixtures(output, mixtures)
        ng = self.n_groups(output, sp)
        rows = []
        found = False
        for m in mix:
            try:
                v = self.isotope_record(output, isotope, key, sp, m, default=None)
            except KeyError:
                v = None  # isotope absent from this mixture
            if v is None:
                rows.append(np.zeros(ng))
                continue
            found = True
            v = np.atleast_1d(np.asarray(v, dtype=float))
            if v.size != ng:
                raise ValueError(f"{output}/{isotope}/{key}: {v.size} values "
                                 f"for {ng} groups.")
            rows.append(v[::-1])  # DRAGON group 1 = highest energy
        if not found:
            if missing == "raise":
                raise KeyError(f"{output}: record '{key}' of '{isotope}' not "
                               "found in any requested mixture.")
            return None
        out = np.vstack(rows)
        return out[0] if isinstance(mixtures, (int, np.integer)) else out

    def micro_xs(self, output: str, isotope: str, reaction: str,
                 sp: StatepointSel = None, mixtures=None,
                 missing: str = "none") -> Optional[np.ndarray]:
        """Microscopic cross section *reaction* (DRAGON key, e.g. ``NG``,
        ``NFTOT``) of *isotope*, ascending energy.

        Returns ``(n_mixtures, n_groups)`` (``(n_groups,)`` for an int
        *mixtures*); mixtures without the isotope/reaction are zero rows.
        Returns ``None`` if the reaction exists in no mixture, unless
        ``missing="raise"``.
        """
        return self._groupwise(output, isotope, reaction, sp, mixtures, missing)

    def flux(self, output: str, sp: StatepointSel = None, mixtures=None,
             isotope: Optional[str] = None) -> np.ndarray:
        """Integrated flux (``NWT0``) per mixture, ascending energy.

        ``NWT0`` is stored per isotope and is identical for all isotopes of
        a mixture; by default the first isotope of each mixture is used.
        """
        mix = self._mixtures(output, mixtures)
        rows = []
        for m in mix:
            iso = isotope or self.isotopes(output, sp, m)[0]
            rows.append(self._groupwise(output, iso, "NWT0", sp, m, "raise"))
        out = np.vstack(rows)
        return out[0] if isinstance(mixtures, (int, np.integer)) else out


__all__ = [
    "CompoReader", "Statepoint", "StatepointError", "decode_tree",
    "open_compo", "build_cache", "cache_is_valid", "default_cache_path",
]
