# -*- coding: utf-8 -*-
"""Registry of Serpent2 references (plan §3.7).

One reference per (assembly, library, statepoint label), many DRAGON variants
compared against it::

    <root>/<assembly>/<library>/<sp_label>/
        <input>                       Serpent input (written or linked)
        <input>.detectors.json        detector manifest (written or frozen)
        <input>_det0.m, _res.m, ...   outputs, once run (or linked)
        reference.yaml                what this reference is

``reference.yaml`` records the origin (``generated`` or ``imported``), the
statepoint values and coolant composition (generated), the Serpent settings
(population, cycles, histories and seed when known), the spec observables and
the SHA-256 of the input. :meth:`SerpentReference.require` checks that a
reference covers a :class:`ComparisonSpec` before post-processing.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import yaml

from ._io import read_text, strip_gz
from .serpent_reader import (DetectorManifest, PinResolver, SerpentRun, _sha256,
                             default_manifest_path, read_res)

PathLike = Union[str, Path]
REFERENCE_FILE = "reference.yaml"
_LABEL_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")
#: Output files of a Serpent run, by suffix appended to the input name.
OUTPUT_SUFFIXES = ("_det0.m", "_res.m", ".out", ".seed")


class ReferenceCoverageError(LookupError):
    """The reference cannot provide some observables of a spec."""


def _check_label(*parts: str) -> None:
    for p in parts:
        if not _LABEL_RE.match(str(p)):
            raise ValueError(f"'{p}' is not a valid folder label ([A-Za-z0-9_.+-]).")


def _serpent_settings(text: str) -> Dict[str, object]:
    """``set pop`` / ``set seed`` of an input (comments stripped)."""
    out: Dict[str, object] = {}
    for line in text.splitlines():
        line = line.split("%")[0].strip()
        w = line.split()
        if len(w) >= 3 and w[0] == "set" and w[1] == "pop":
            out.update(neutrons_per_cycle=int(float(w[2])))
            if len(w) >= 5:
                out.update(active_cycles=int(float(w[3])), inactive_cycles=int(float(w[4])))
        elif len(w) >= 3 and w[0] == "set" and w[1] == "seed":
            out["seed"] = int(float(w[2]))
    return out


@dataclass
class SerpentReference:
    """One reference folder."""

    path: Path
    meta: Dict[str, object]

    # ── creation ───────────────────────────────────────────────────────────
    @staticmethod
    def folder(root: PathLike, assembly: str, library: str, sp_label: str) -> Path:
        _check_label(assembly, library, sp_label)
        return Path(root) / assembly / library / sp_label

    @classmethod
    def create(cls, root: PathLike, assembly: str, library: str, sp_label: str,
               model, spec=None, input_name: Optional[str] = None,
               overwrite: bool = False) -> "SerpentReference":
        """Write the input of *model* (``Serpent2Model``, e.g. from
        ``from_dragon_case``) and its manifest into a new reference folder."""
        d = cls.folder(root, assembly, library, sp_label)
        if (d / REFERENCE_FILE).exists() and not overwrite:
            raise FileExistsError(f"Reference {d} exists (overwrite=True to replace it).")
        d.mkdir(parents=True, exist_ok=True)
        name = input_name or f"{assembly}_{sp_label}_{library}.serp"
        inp = d / name
        model.write(str(inp))
        s = model.settings
        meta = {
            "assembly": assembly, "library": library, "sp_label": sp_label,
            "origin": "generated",
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "input": name, "input_sha256": _sha256(inp),
            "statepoint": {k: float(v) for k, v in (getattr(model, "statepoint", None) or {}).items()},
            "coolant_composition": getattr(model, "coolant_composition", None),
            "settings": {k: getattr(s, k) for k in ("neutrons_per_cycle", "active_cycles",
                                                     "inactive_cycles") if hasattr(s, k)},
            "spec": None if spec is None else spec.names(),
            "manifest": default_manifest_path(inp).name
            if default_manifest_path(inp).exists() else None,
        }
        (d / REFERENCE_FILE).write_text(yaml.safe_dump(meta, sort_keys=False))
        return cls(d, meta)

    @classmethod
    def import_existing(cls, serp_path: PathLike, root: PathLike, assembly: str,
                        library: str, sp_label: str, link: bool = True,
                        pin_resolver: Optional[PinResolver] = None,
                        manifest: Optional[DetectorManifest] = None,
                        overwrite: bool = False) -> "SerpentReference":
        """Register an existing run (decision 8: no re-run). The input and
        its outputs are symlinked (*link*) or copied; the manifest is the
        frozen one next to the run if any, else *manifest*, else inferred
        from the card content and frozen in the reference folder."""
        src = Path(serp_path).resolve()
        if not src.exists():
            raise FileNotFoundError(src)
        d = cls.folder(root, assembly, library, sp_label)
        if (d / REFERENCE_FILE).exists() and not overwrite:
            raise FileExistsError(f"Reference {d} exists (overwrite=True to replace it).")
        d.mkdir(parents=True, exist_ok=True)
        base = strip_gz(src)
        files = [src] + [f for suf in OUTPUT_SUFFIXES
                         for f in (base.with_name(base.name + suf),
                                   base.with_name(base.name + suf + ".gz")) if f.exists()]
        for f in files:
            target = d / f.name
            if target.exists() or target.is_symlink():
                target.unlink()
            if link:
                os.symlink(f, target)
            else:
                shutil.copy2(f, target)
        inp = d / src.name
        frozen = default_manifest_path(src)
        if manifest is None and frozen.exists():
            manifest = DetectorManifest.read(frozen)
        if manifest is None:
            manifest = DetectorManifest.infer(src, pin_resolver=pin_resolver)
        manifest.write(default_manifest_path(inp))
        settings = _serpent_settings(read_text(src))
        res_path = next((f for f in files if f.name.startswith(base.name + "_res.m")), None)
        if res_path is not None:
            res = read_res(res_path)
            for key, name in (("SIMULATED_HISTORIES", "simulated_histories"),
                              ("CYCLES", "cycles")):
                v = res.get(key)
                if v is not None:
                    settings[name] = int(float(v if not hasattr(v, "__len__") else v[0]))
        meta = {
            "assembly": assembly, "library": library, "sp_label": sp_label,
            "origin": "imported",
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "source": str(src), "linked": link,
            "input": src.name, "input_sha256": _sha256(src),
            "settings": settings,
            "files": [f.name for f in files],
            "manifest": default_manifest_path(inp).name,
        }
        (d / REFERENCE_FILE).write_text(yaml.safe_dump(meta, sort_keys=False))
        return cls(d, meta)

    # ── access ─────────────────────────────────────────────────────────────
    @classmethod
    def open(cls, path: PathLike) -> "SerpentReference":
        d = Path(path)
        return cls(d, yaml.safe_load((d / REFERENCE_FILE).read_text()))

    @property
    def input_path(self) -> Path:
        return self.path / self.meta["input"]

    @property
    def manifest(self) -> DetectorManifest:
        return DetectorManifest.read(default_manifest_path(self.input_path))

    @property
    def has_results(self) -> bool:
        base = strip_gz(self.input_path)
        return all(any(base.with_name(base.name + suf + gz).exists() for gz in ("", ".gz"))
                   for suf in ("_det0.m", "_res.m"))

    def run(self) -> SerpentRun:
        """:class:`SerpentRun` on this reference (needs the outputs)."""
        if not self.has_results:
            raise FileNotFoundError(f"{self.path}: no _det0.m/_res.m yet (run Serpent2).")
        return SerpentRun(self.input_path, manifest=self.manifest)

    def source(self, spec, observable: str, kind: str = "rr",
               n_pins: Optional[int] = None, pins: Optional[Sequence[int]] = None) -> str:
        """Observable of this reference's detectors that provides
        *observable*/*kind* of *spec* (the ``observable`` argument of
        ``compute_pin_rates`` / ``compare_spectra``)."""
        return spec.serpent_source(observable, kind, manifest=self.manifest,
                                   n_pins=n_pins, pins=pins)

    def coverage(self, spec, n_pins: Optional[int] = None,
                 pins: Optional[Sequence[int]] = None) -> Dict[str, List[str]]:
        return spec.coverage(self.manifest, n_pins=n_pins, pins=pins)

    def require(self, spec, n_pins: Optional[int] = None,
                pins: Optional[Sequence[int]] = None) -> None:
        """Raise :class:`ReferenceCoverageError` unless the reference covers
        *spec* (every observable, every pin)."""
        gaps = self.coverage(spec, n_pins=n_pins, pins=pins)
        if gaps:
            lines = [f"{k}: {'; '.join(v)}" for k, v in gaps.items()]
            raise ReferenceCoverageError(
                f"Reference {self.path} does not cover the spec:\n  " + "\n  ".join(lines))


class ReferenceRegistry:
    """All references under *root*."""

    def __init__(self, root: PathLike):
        self.root = Path(root)

    def list(self) -> List[SerpentReference]:
        return [SerpentReference.open(p.parent)
                for p in sorted(self.root.glob(f"*/*/*/{REFERENCE_FILE}"))]

    def find(self, assembly: str, library: str, sp_label: str) -> SerpentReference:
        d = SerpentReference.folder(self.root, assembly, library, sp_label)
        if not (d / REFERENCE_FILE).exists():
            have = [f"{r.meta['assembly']}/{r.meta['library']}/{r.meta['sp_label']}"
                    for r in self.list()]
            raise KeyError(f"No reference {assembly}/{library}/{sp_label}; registered: {have}")
        return SerpentReference.open(d)


__all__ = ["SerpentReference", "ReferenceRegistry", "ReferenceCoverageError",
           "REFERENCE_FILE"]
