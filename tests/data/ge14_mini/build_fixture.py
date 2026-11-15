#!/usr/bin/env python
"""Build the committed GE14 "mini" test fixture from the full local samples.

The full GE14 DOM results (~1.5 GB, ``temp_examples/GE14_DOM``) stay local.
This script extracts a small, real subset (~0.5 MB) used by the always-run
tests, and verifies that the subset reads back identically to the source:

- COMPO: ``SIGNATURE`` + ``EDIHOM_COND`` + ``EDIR_2G`` top-level entries,
  copied verbatim (byte ranges of the ASCII export), gzipped.
- Serpent2: the unmodified legacy ``.serp`` input and ``res.m`` (gzipped);
  ``det0.m`` reduced to the per-pin detectors of ``--pins`` plus the
  assembly-integrated and global flux detectors (blocks copied verbatim);
  a frozen detector manifest restricted to those detectors.
- Golden values: ``2G_data.npz``, ``density_00_*.txt`` and the 0 % void rows
  of ``sortie_DOM_1L_PT_3x3_HFP.csv`` produced by the original analysis
  scripts, plus ``reference_values.json`` (keff, pcm).
- ``SOURCES.json``: provenance (relative source paths, sizes, SHA-256).

Usage (from the repository root)::

    python tests/data/ge14_mini/build_fixture.py [--src temp_examples/GE14_DOM]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
from datetime import date
from pathlib import Path

import numpy as np

from starterDD.PostTreatment.compo_reader import CompoReader
from starterDD.PostTreatment.pylcm import Pylcm
from starterDD.PostTreatment.serpent_reader import (
    DetectorManifest, SerpentRun, default_manifest_path,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

RUN_REL = Path("results/coolant_density_1L_PT_3x3_HFP/GE14_DOM_DIAG/PT_IC_1L_MOC")
RUN_00 = "2026-07-19T10-46-27"
CPO_NAME = "_CPO_GE14_DOM_DIAG"
SERP_NAME = "GE14_DOM_00_endfb8r1.serp"
CSV_NAME = "sortie_DOM_1L_PT_3x3_HFP.csv"
COMPO_KEEP = ("SIGNATURE", "EDIHOM_COND", "EDIR_2G")
# UOX16 (on axis), UOX28 (off axis), UOX40Gd8 (Gd, off axis), UOX49Gd8 (Gd, on axis)
DEFAULT_PINS = (1, 2, 15, 28)

_TOP_RE = re.compile(rb"^->[ \t]*1[ \t]+(\d+)[ \t]+-?\d+[ \t]+-?\d+[ \t]*<-[^\n]*\n([^\n]*)", re.M)
_ROOT_CLOSE_RE = re.compile(rb"^->[ \t]*-1[ \t]+0[ \t]+0[ \t]+0[ \t]*<-[^\n]*\n?", re.M)
_DET_RE = re.compile(r"^DET(\S+?)\s*=\s*\[.*?\];\n", re.M | re.S)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def write_gz(path: Path, data: bytes) -> Path:
    """Deterministic gzip (no timestamp / file name in the header)."""
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb",
                                                 compresslevel=9, mtime=0,
                                                 filename="") as gz:
        gz.write(data)
    return path


def trim_compo(src: Path, keep) -> bytes:
    buf = src.read_bytes()
    tops = [(m.start(), m.group(2).strip().decode()) for m in _TOP_RE.finditer(buf)]
    close = list(_ROOT_CLOSE_RE.finditer(buf))
    end_all = close[-1].start() if close else len(buf)
    out = []
    for i, (start, key) in enumerate(tops):
        stop = tops[i + 1][0] if i + 1 < len(tops) else end_all
        if key in keep:
            out.append(buf[start:stop])
    missing = set(keep) - {k for _, k in tops}
    if missing:
        raise KeyError(f"COMPO has no top-level entries {sorted(missing)}")
    out.append(close[-1].group(0) if close else b"")
    return b"".join(out)


def same_records(a: Pylcm, b: Pylcm) -> bool:
    wa, wb = list(a.walk()), list(b.walk())
    if [p for p, _ in wa] != [p for p, _ in wb]:
        return False
    for (_, x), (_, y) in zip(wa, wb):
        if isinstance(x, np.ndarray):
            if not np.array_equal(x, y):
                return False
        elif x != y:
            return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", type=Path, default=REPO / "temp_examples" / "GE14_DOM")
    ap.add_argument("--out", type=Path, default=HERE)
    ap.add_argument("--pins", type=int, nargs="+", default=list(DEFAULT_PINS))
    args = ap.parse_args()
    src, out = args.src.resolve(), args.out.resolve()
    run_dir = src / RUN_REL / RUN_00
    serp_dir = src / "Serpent"
    (out / "golden").mkdir(parents=True, exist_ok=True)
    sources = {}

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(REPO))
        except ValueError:
            return str(p)

    def record(p: Path):
        sources[rel(p)] = {"size": p.stat().st_size, "sha256": sha256(p)}

    # ── COMPO ──────────────────────────────────────────────────────────────
    cpo_src = run_dir / CPO_NAME
    record(cpo_src)
    cpo_out = write_gz(out / f"{CPO_NAME}_mini.gz", trim_compo(cpo_src, COMPO_KEEP))
    full = Pylcm.from_ascii(cpo_src, only=COMPO_KEEP)
    mini = Pylcm.from_ascii(cpo_out)
    assert mini.keys() == list(COMPO_KEEP), mini.keys()
    assert same_records(full, mini), "trimmed COMPO differs from the source"
    print(f"COMPO   {cpo_out.name}: {cpo_out.stat().st_size / 1024:.0f} KB, records identical")

    # ── Serpent ────────────────────────────────────────────────────────────
    serp_src = serp_dir / SERP_NAME
    det_src = serp_dir / f"{SERP_NAME}_det0.m"
    res_src = serp_dir / f"{SERP_NAME}_res.m"
    for p in (serp_src, det_src, res_src):
        record(p)
    serp_out = write_gz(out / f"{SERP_NAME}.gz", serp_src.read_bytes())
    write_gz(out / f"{SERP_NAME}_res.m.gz", res_src.read_bytes())

    full_man = DetectorManifest.infer(serp_src)
    keep = [e for e in full_man
            if (e.spatial == "PIN" and e.pin_idx in args.pins)
            or e.spatial in ("FUEL", "ALL")]
    keep_names = {e.detector for e in keep}
    text = det_src.read_text()
    blocks = []
    for m in _DET_RE.finditer(text):
        name = m.group(1)
        base = name[:-1] if name.endswith("E") and name[:-1] in keep_names else name
        if base in keep_names:
            blocks.append(m.group(0))
    det_out = write_gz(out / f"{SERP_NAME}_det0.m.gz", ("\n" + "".join(blocks)).encode())

    man = DetectorManifest.infer(serp_out)          # hash of the committed input
    man.entries = [e for e in man if e.detector in keep_names]
    man.source.update({"subset_of_det0_sha256": sources[rel(det_src)]["sha256"],
                       "subset_pins": sorted(args.pins)})
    man_path = man.write(default_manifest_path(serp_out))

    run_full = SerpentRun(serp_src, manifest=full_man)
    run_mini = SerpentRun(serp_out)
    assert run_mini.manifest.source.get("subset_pins") == sorted(args.pins)
    assert run_mini.check() == [], run_mini.check()
    for e in run_mini.manifest:
        a, b = run_full.tally_of(full_man.by_detector(e.detector)), run_mini.tally_of(e)
        assert np.array_equal(a.values, b.values) and np.array_equal(a.rel_err, b.rel_err)
    assert run_mini.keff() == run_full.keff()
    print(f"Serpent {det_out.name}: {det_out.stat().st_size / 1024:.0f} KB, "
          f"{len(keep)} detectors identical; manifest {man_path.name}")

    # ── Golden values from the original analysis ───────────────────────────
    for name in ["2G_data.npz"] + sorted(p.name for p in run_dir.glob("density_00_*.txt")):
        p = run_dir / name
        record(p)
        shutil.copyfile(p, out / "golden" / name)
    shutil.copyfile(run_dir / "run_manifest.yaml", out / "golden" / "dragon_run_manifest.yaml")
    record(run_dir / "run_manifest.yaml")
    csv = src / RUN_REL / CSV_NAME
    record(csv)
    k_dragon = CompoReader(mini).keff("EDIHOM_COND")
    k_serp, k_serp_rel = run_mini.keff()
    rows = [ln for ln in csv.read_text().splitlines()
            if ln and abs(float(ln.split(";")[1]) - k_dragon) < 1e-12]
    (out / "golden" / "sortie_DOM_1L_00.csv").write_text("\n".join(rows) + "\n")
    ref = {
        "keff_dragon_EDIHOM_COND": k_dragon,
        "keff_serpent_IMP_KEFF": k_serp,
        "keff_serpent_IMP_KEFF_rel_err": k_serp_rel,
        "delta_keff_pcm": (k_dragon - k_serp) * 1e5,
        "statepoint": {"TFuel": 900.0, "TCool": 600.0, "DCool": 0.73669, "label": "00"},
    }
    (out / "golden" / "reference_values.json").write_text(json.dumps(ref, indent=1) + "\n")
    print(f"Golden  keff {k_dragon} vs {k_serp} -> {ref['delta_keff_pcm']:.2f} pcm, "
          f"{len(rows)} csv rows")

    (out / "SOURCES.json").write_text(json.dumps({
        "built": str(date.today()),
        "builder": rel(Path(__file__).resolve()),
        "pins": sorted(args.pins),
        "compo_entries": list(COMPO_KEEP),
        "sources": sources,
    }, indent=1) + "\n")
    total = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"Fixture total: {total / 1024:.0f} KB in {out}")


if __name__ == "__main__":
    main()
