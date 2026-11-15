# -*- coding: utf-8 -*-
"""Small I/O helpers shared by the readers: transparent ``.gz`` support.

Test fixtures (``tests/data/ge14_mini``) store DRAGON and Serpent2 outputs
gzipped; every reader accepts either form.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]


def is_gz(path: PathLike) -> bool:
    return str(path).endswith(".gz")


def read_bytes(path: PathLike) -> bytes:
    path = Path(path)
    if is_gz(path):
        with gzip.open(path, "rb") as f:
            return f.read()
    return path.read_bytes()


def read_text(path: PathLike) -> str:
    return read_bytes(path).decode("utf-8", "replace")


def strip_gz(path: PathLike) -> Path:
    """``a/b.serp.gz`` -> ``a/b.serp``; other paths unchanged."""
    p = Path(path)
    return p.with_name(p.name[:-3]) if is_gz(p) else p


def existing_variant(path: PathLike) -> Path:
    """Return *path* if it exists, else ``<path>.gz`` if that exists, else
    *path* (so that the caller's error names the expected file)."""
    p = Path(path)
    if p.exists():
        return p
    gz = p.with_name(p.name + ".gz")
    return gz if gz.exists() else p
