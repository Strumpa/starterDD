# -*- coding: utf-8 -*-
# __author__ = B. Daré (original Pylcm), extended for starterDD post-treatment
# __date__ = 2026-04-17

"""Python representation of DRAGON/DONJON LCM data structures.

:class:`Pylcm` mirrors an LCM object (associative table or list) read from a
DRAGON ASCII export (``_CPO_*`` files and similar). On top of the original
API (``content``, ``add``, ``get``, ``print_tree``, ``generate_from_ascii``),
this module adds:

- ``get("A/B/3/C")`` path access (list items addressed by 1-based index),
  ``keys()``, ``in`` and ``[]``;
- a single-pass iterative ASCII parser working on a memory map, with exact
  decoding of character records and support for logical / complex records;
- **selective loading** of top-level directories (``only=[...]``) that skips
  the other directories without parsing them;
- :meth:`Pylcm.to_hdf5` / :meth:`Pylcm.from_hdf5` and the lazy
  :class:`LcmH5` view, so that a large COMPO is parsed once and then read
  partially (one output, one record) from an HDF5 cache.

ASCII format (one node)::

    ->   level  name_type  content_type  content_size   <-   [list index]
    KEY (only when name_type > 0, on its own line)
    data...

``content_type``: 0 table, 1 integer, 2 real, 3 character (4-char words),
4 double, 5 logical, 6 complex, 10 list, 99 undefined. ``content_size`` is -1
for a table, the number of items for a list, the number of elements (words
for characters) otherwise. A table whose children are at level ``L`` is
closed by ``-> -L 0 0 0 <-``; lists are closed implicitly after their items.
"""

from __future__ import annotations

import mmap
import os
import re
import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Union

import numpy as np

from ._io import is_gz, read_bytes

PathLike = Union[str, Path]

_MISSING = object()

# LCM content types
CT_TABLE, CT_INT, CT_REAL, CT_CHAR, CT_DOUBLE, CT_LOGICAL, CT_COMPLEX = 0, 1, 2, 3, 4, 5, 6
CT_LIST, CT_UNDEF = 10, 99


class LcmChars(str):
    """A character record.

    Behaves as the stripped string (backward compatible with the original
    parser), and keeps the exact fixed-width content in :attr:`raw` so that
    arrays of names can be split by width with :meth:`chunks`.
    """

    raw: str

    def __new__(cls, raw: str):
        obj = super().__new__(cls, raw.strip())
        obj.raw = raw
        return obj

    def chunks(self, width: int) -> List[str]:
        """Split :attr:`raw` into fixed-width fields, stripped."""
        r = self.raw
        return [r[i:i + width].strip() for i in range(0, len(r), width)]

    def __reduce__(self):
        return (LcmChars, (self.raw,))


def _split_path(path: Union[str, int, Sequence]) -> List[str]:
    if isinstance(path, int):
        return [str(path)]
    if isinstance(path, (list, tuple)):
        return [str(p) for p in path]
    return [p for p in str(path).split("/") if p != ""]


class Pylcm:
    """
    This class allows the user to represent LCM data structures, used in the
    codes DRAGON and DONJON. A LCM data structure can be interpreted as an
    associative table, containing heterogeneous elements by associating a
    12-character key and a content, or as a list, containing heterogeneous
    contents identified by their index.

    A LCM data structure can contain other LCM data structures (associative
    tables or lists), or elementary objects:

    - integer
    - single or double precision real number
    - character string
    - logical number (0 or 1)
    - complex number

    Note that an elementary content in a LCM data structure can be an
    elementary object alone, or a set of several elementary objects of the
    same type.

    Those properties allow consequently tree structures.

    List items are stored with integer keys (1-based, as in DRAGON).
    Elementary records of size 1 are stored as Python scalars, larger ones as
    numpy arrays, characters as :class:`LcmChars`.
    """

    def __init__(self, name='Pylcm_object', level=0, type='table', size=np.inf, parent=None):
        if len(name) > 12:
            raise ValueError("Name must be a 12-character string")
        self.name = name # 12-character name
        self.content = {}
        self.level = level # level in the global tree
        self.type = type # 'table' or 'list'
        self.size = size # np.inf if table, len(list) if list
        self.parent = parent # Pylcm parent, None if root
        self.ctypes = {} # key -> LCM content type of elementary records

    # ── basic access ────────────────────────────────────────────────────────
    def add(self, key: str, object, ctype: Optional[int] = None):
        if self.type == 'list':
            key = int(key)
        self.content[key] = object
        if ctype is not None:
            self.ctypes[key] = ctype

    def _get1(self, key):
        if self.type == 'list':
            return self.content[int(key)]
        return self.content[key]

    def get(self, key, default=KeyError):
        """Return a child by key, or by path ``"A/B/3/C"``.

        Keys containing ``/`` are first tried verbatim. Raises ``KeyError``
        unless *default* is given.
        """
        if isinstance(key, str) and key.strip("/") == "" and key not in self.content:
            return self
        try:
            if isinstance(key, str) and "/" in key and key not in self.content:
                node = self
                for part in _split_path(key):
                    if not isinstance(node, Pylcm):
                        raise KeyError(key)
                    node = node._get1(part)
                return node
            return self._get1(key)
        except (KeyError, ValueError):
            if default is KeyError:
                raise KeyError(f"{key!r} not found in {self.path or '<root>'}") from None
            return default

    def __getitem__(self, key):
        return self.get(key)

    def __contains__(self, key) -> bool:
        return self.get(key, _MISSING) is not _MISSING

    def keys(self) -> List:
        return list(self.content.keys())

    def items(self):
        return self.content.items()

    def __len__(self) -> int:
        return len(self.content)

    def __iter__(self) -> Iterator:
        return iter(self.content)

    @property
    def path(self) -> str:
        """Path of this node from the root (``""`` for the root)."""
        parts = []
        node = self
        while node.parent is not None:
            parent = node.parent
            for k, v in parent.content.items():
                if v is node:
                    parts.append(str(k))
                    break
            node = parent
        return "/".join(reversed(parts))

    def walk(self, prefix: str = "") -> Iterator[tuple]:
        """Yield ``(path, value)`` for every elementary record below."""
        for k, v in self.content.items():
            p = f"{prefix}/{k}" if prefix else str(k)
            if isinstance(v, Pylcm):
                yield from v.walk(p)
            else:
                yield p, v

    # ── display ─────────────────────────────────────────────────────────────
    def print_tree(self, path: str, name: str=None):
        """

        This method creates a txt file at the specified *path* representing the
        tree structure of a ``Pylcm`` object (and that of any other child
        ``Pylcm`` object contained).

        Args:
            *path* (str): The path to the folder where the output file will
                          be created.
            *name* (str): Optional, name of the output file. If None, the
                          function will give the name of the ``Pylcm`` object.

        Returns:
            ``None``

        """

        # Define the name of the file
        if name is None:
            name = self.name

        # Create the output file
        with open(f"{path}/{name}.txt", 'w') as f:
            f.write(self.name.ljust(14,'_'))
            self._print_level(f, self.level, [])


    def _print_level(self, f, root_level, state_list):
        """

        Internal method used by the method ``print_tree`` to recursively
        generate the tree structure in the output file by browsing the ``Pylcm``
        object.

        Args:
            *f*: The file object where to write.
            *root_level*: The level of the root, in other words the level of the
                          first ``Pylcm`` object on which the function
                          ``print_tree`` is called.
            *state_list*: List of boolean values used to determine if at a
                          specified level (linked to the index in the list) all
                          the branches of the current node have been explored
                          (True) or not (False). It is used for display purpuse
                          only.

        Returns:
            ``None``

        """

        if (self.level - root_level) < len(state_list):
            state_list[self.level - root_level] = False
        else:
            state_list.append(False)
        children = list(self.content.keys())
        for i in range(len(children)):
            if i == 0:
                ch = '_'
            else:
                ch = '|'
                f.write(' '*14)
                for j in range(self.level - root_level):
                    if not state_list[j]:
                        f.write('|'.ljust(17,' '))
                    else:
                        f.write(' '*17)
            if i == len(children)-1:
                state_list[self.level - root_level] = True
            if isinstance(self.get(children[i]), Pylcm):
                f.write(ch + f"__{children[i]}".ljust(16,'_'))
                self.get(children[i])._print_level(f, root_level, state_list)
            else:
                f.write(ch + f"__{children[i]}\n")

    # ── ASCII import ────────────────────────────────────────────────────────
    @classmethod
    def from_ascii(cls, path: PathLike, name: str = "Pylcm_object",
                   only: Optional[Iterable[str]] = None) -> "Pylcm":
        """Parse a DRAGON ASCII LCM export.

        Args:
            path: File to read (``.gz`` files are decompressed in memory).
            name: Name of the returned root object (≤ 12 chars).
            only: Names of the top-level entries to load. Other top-level
                directories are skipped without being parsed. ``None`` loads
                everything.
        """
        root = cls(name=name, level=0)
        if is_gz(path):
            _AsciiParser(read_bytes(path), root, only).run()
            return root
        with open(path, "rb") as f:
            size = os.fstat(f.fileno()).st_size
            if size == 0:
                return root
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buf:
                _AsciiParser(buf, root, only).run()
        return root

    @classmethod
    def generate_from_ascii(cls, path: str, **kwargs):
        """Original entry point, kept for compatibility. Accepts ``name``
        and ``only`` keyword arguments (see :meth:`from_ascii`)."""
        name = kwargs.get('name', 'Pylcm_object')
        if len(name) > 12:
            raise ValueError("Name must be a 12-character string")
        return cls.from_ascii(path, name=name, only=kwargs.get('only'))

    # ── HDF5 ────────────────────────────────────────────────────────────────
    def to_hdf5(self, path: PathLike, compression: Optional[str] = "gzip",
                attrs: Optional[Dict[str, Any]] = None) -> Path:
        """Write this object to an HDF5 file.

        Layout: tables and lists are groups (attribute ``lcm_type``), list
        items are named ``00000001``, ``00000002``..., elementary records are
        datasets with an ``lcm_ctype`` attribute. Root attributes from
        *attrs* are added (e.g. source file metadata).
        """
        h5py = _require_h5py()
        path = Path(path)
        tmp = path.with_name(path.name + ".part")
        with h5py.File(tmp, "w", track_order=True) as h5:
            h5.attrs["lcm_format"] = "starterDD-pylcm"
            h5.attrs["lcm_format_version"] = 1
            h5.attrs["name"] = self.name
            for k, v in (attrs or {}).items():
                h5.attrs[k] = v
            _write_group(h5, self, compression)
        os.replace(tmp, path)
        return path

    @classmethod
    def from_hdf5(cls, path: PathLike, only: Optional[Iterable[str]] = None,
                  name: Optional[str] = None) -> "Pylcm":
        """Load an HDF5 file written by :meth:`to_hdf5` into memory
        (optionally only some top-level entries)."""
        h5py = _require_h5py()
        with h5py.File(path, "r") as h5:
            root = cls(name=name or str(h5.attrs.get("name", "Pylcm_object")))
            only = None if only is None else set(only)
            for key in h5:
                if only is not None and _h5_unescape(key) not in only:
                    continue
                _read_into(root, key, h5[key], level=1)
        return root

    def __repr__(self):
        if self.type == 'list':
            return f"Pylcm(name={self.name!r}, type={self.type!r}, size={self.size!r}, parent={self.parent!r})"
        else:
            return f"Pylcm(name={self.name!r}, type={self.type!r}, parent={self.parent!r})"


# ═════════════════════════════════════════════════════════════════════════════
# ASCII parser
# ═════════════════════════════════════════════════════════════════════════════

_HEADER_RE = re.compile(
    rb"^->[ \t]*(-?\d+)[ \t]+(-?\d+)[ \t]+(-?\d+)[ \t]+(-?\d+)[ \t]*<-([^\r\n]*)\r?$",
    re.M,
)
_TOP_HEADER_RE = re.compile(rb"^->[ \t]*1[ \t]+", re.M)
_FORTRAN_EXP_RE = re.compile(r"(?<=[0-9.])([+-]\d{2,3})(?=\s|$)")


def _to_floats(text: str) -> np.ndarray:
    tokens = text.split()
    try:
        return np.array(tokens, dtype=float)
    except ValueError:
        # Fortran drops the "E" for 3-digit exponents ("1.0-100") and may
        # use "D" for double precision.
        fixed = _FORTRAN_EXP_RE.sub(r"E\1", text.replace("D", "E").replace("d", "E"))
        return np.array(fixed.split(), dtype=float)


class _AsciiParser:
    """Single-pass, stack-based parser (no recursion)."""

    def __init__(self, buf, root: Pylcm, only: Optional[Iterable[str]]):
        self.buf = buf
        self.root = root
        self.only = None if only is None else set(only)
        # Stack of containers; a container's children are at level+1.
        self.stack: List[Pylcm] = [root]

    def run(self) -> None:
        buf = self.buf
        m = _HEADER_RE.search(buf, 0)
        while m is not None:
            level, ntype, ctype, csize = (int(x) for x in m.group(1, 2, 3, 4))
            rest = m.group(5)
            nxt = _HEADER_RE.search(buf, m.end())
            payload_end = nxt.start() if nxt is not None else len(buf)

            if level < 0:
                self._close(-level)
                m = nxt
                continue

            self._unwind_to(level)
            payload = buf[m.end():payload_end].decode("ascii", "replace")
            lines = payload.split("\n")
            if lines and lines[0].strip() == "":
                lines = lines[1:]  # remainder of the header line
            if ntype > 0:
                key = lines[0][:ntype].strip() if lines else ""
                lines = lines[1:]
            else:
                key = rest.decode("ascii", "replace").strip()

            parent = self.stack[-1]
            if level == 1 and self.only is not None and key not in self.only:
                # Skip this whole top-level entry without parsing it.
                skip = _TOP_HEADER_RE.search(buf, m.end())
                m = _HEADER_RE.search(buf, skip.start()) if skip else None
                continue

            if ctype == CT_TABLE:
                node = Pylcm(name=key[:12] or "item", level=level, type="table",
                             parent=parent)
                parent.add(key, node)
                self.stack.append(node)
            elif ctype == CT_LIST:
                node = Pylcm(name=key[:12] or "item", level=level, type="list",
                             size=csize, parent=parent)
                parent.add(key, node)
                self.stack.append(node)
            else:
                parent.add(key, self._decode(ctype, csize, lines), ctype)
            m = nxt

    def _unwind_to(self, level: int) -> None:
        # Pop finished containers (lists have no explicit close).
        while len(self.stack) > 1 and self.stack[-1].level + 1 > level:
            self.stack.pop()

    def _close(self, child_level: int) -> None:
        self._unwind_to(child_level)
        if len(self.stack) > 1 and self.stack[-1].level + 1 == child_level:
            self.stack.pop()

    @staticmethod
    def _decode(ctype: int, size: int, lines: List[str]):
        if ctype == CT_UNDEF:
            return None
        if size <= 0:
            return LcmChars("") if ctype == CT_CHAR else np.zeros(0)
        if ctype == CT_CHAR:
            # Word lengths first (8 per line), then the text in 80-col lines.
            n_len_lines = -(-size // 8)
            text_lines = [ln.rstrip("\r") for ln in lines[n_len_lines:]]
            while text_lines and text_lines[-1] == "":
                text_lines.pop()
            nchar = 4 * size
            raw = "".join(ln.ljust(80) for ln in text_lines[:-1])
            raw += text_lines[-1] if text_lines else ""
            return LcmChars(raw[:nchar].ljust(nchar))
        text = " ".join(lines)
        if ctype in (CT_INT, CT_LOGICAL):
            arr = np.array(text.split(), dtype=np.int64)
            if ctype == CT_LOGICAL:
                arr = arr.astype(bool)
            return arr.item() if size == 1 else arr
        if ctype in (CT_REAL, CT_DOUBLE):
            arr = _to_floats(text)
            return float(arr[0]) if size == 1 else arr
        if ctype == CT_COMPLEX:
            arr = _to_floats(text)
            arr = arr[0::2] + 1j * arr[1::2]
            return complex(arr[0]) if size == 1 else arr
        warnings.warn(f"Unknown LCM content type {ctype}; record kept as text.")
        return text


# ═════════════════════════════════════════════════════════════════════════════
# HDF5 backend
# ═════════════════════════════════════════════════════════════════════════════


def _require_h5py():
    try:
        import h5py  # noqa: F401
    except ImportError as exc:  # pragma: no cover - depends on env
        raise ImportError(
            "h5py is required for HDF5 support: pip install -e '.[post]'"
        ) from exc
    return h5py


def h5py_available() -> bool:
    try:
        import h5py  # noqa: F401
        return True
    except ImportError:
        return False


def _h5_escape(key) -> str:
    if isinstance(key, int):
        return f"{key:08d}"
    s = str(key).replace("%", "%25").replace("/", "%2F")
    return "%2E" if s == "." else s


def _h5_unescape(key: str) -> str:
    return key.replace("%2F", "/").replace("%2E", ".").replace("%25", "%")


def _write_group(group, node: Pylcm, compression) -> None:
    group.attrs["lcm_type"] = node.type
    if node.type == "list":
        group.attrs["size"] = int(node.size) if np.isfinite(node.size) else -1
    for key, val in node.content.items():
        name = _h5_escape(key)
        if isinstance(val, Pylcm):
            # track_order keeps the LCM key order (h5py sorts by name otherwise)
            _write_group(group.create_group(name, track_order=True), val,
                         compression)
            continue
        ctype = node.ctypes.get(key)
        if val is None:
            ds = group.create_dataset(name, data=np.zeros(0, dtype=np.int8))
            ctype = CT_UNDEF
        elif isinstance(val, str):
            raw = val.raw if isinstance(val, LcmChars) else val
            ds = group.create_dataset(name, data=np.bytes_(raw.encode("ascii", "replace")))
            ctype = CT_CHAR
        else:
            arr = np.asarray(val)
            kw = {}
            if compression and arr.size > 64:
                kw = {"compression": compression, "shuffle": True}
            ds = group.create_dataset(name, data=arr, **kw)
            if ctype is None:
                ctype = (CT_LOGICAL if arr.dtype == bool else CT_INT
                         if np.issubdtype(arr.dtype, np.integer) else CT_COMPLEX
                         if np.iscomplexobj(arr) else CT_REAL)
        ds.attrs["lcm_ctype"] = int(ctype)


def _decode_dataset(ds):
    ctype = int(ds.attrs.get("lcm_ctype", CT_REAL))
    if ctype == CT_UNDEF:
        return None
    if ctype == CT_CHAR:
        v = ds[()]
        return LcmChars(v.decode("ascii", "replace") if isinstance(v, bytes) else str(v))
    v = ds[()]
    if np.ndim(v) == 0:
        v = v.item()
    return v


def _read_into(parent: Pylcm, h5key: str, obj, level: int) -> None:
    import h5py
    key = _h5_unescape(h5key)
    if isinstance(obj, h5py.Group):
        ltype = obj.attrs.get("lcm_type", "table")
        ltype = ltype.decode() if isinstance(ltype, bytes) else str(ltype)
        size = int(obj.attrs.get("size", -1)) if ltype == "list" else np.inf
        node = Pylcm(name=key[:12] if parent.type != "list" else "item",
                     level=level, type=ltype, size=size, parent=parent)
        parent.add(key, node)
        for k in obj:
            _read_into(node, k, obj[k], level + 1)
    else:
        parent.add(key, _decode_dataset(obj), int(obj.attrs.get("lcm_ctype", CT_REAL)))


class LcmH5:
    """Lazy, read-only view of an HDF5 file written by :meth:`Pylcm.to_hdf5`.

    Offers the same ``get(path)`` / ``keys(path)`` access as :class:`Pylcm`
    but only reads the requested records. Groups are returned as
    :class:`LcmH5` sub-views. Use as a context manager or call :meth:`close`.
    """

    def __init__(self, source, _root=None):
        h5py = _require_h5py()
        if isinstance(source, (str, Path)):
            self._file = h5py.File(source, "r")
            self._group = self._file
        else:
            self._file = _root
            self._group = source

    @property
    def attrs(self):
        return self._group.attrs

    @property
    def type(self) -> str:
        t = self._group.attrs.get("lcm_type", "table")
        return t.decode() if isinstance(t, bytes) else str(t)

    def _h5name(self, part: str, group) -> str:
        t = group.attrs.get("lcm_type", "table")
        t = t.decode() if isinstance(t, bytes) else str(t)
        if t == "list":
            return f"{int(part):08d}"
        return _h5_escape(part)

    def get(self, path="", default=KeyError):
        import h5py
        obj = self._group
        try:
            for part in _split_path(path):
                if not isinstance(obj, h5py.Group):
                    raise KeyError(path)
                obj = obj[self._h5name(part, obj)]
        except (KeyError, ValueError):
            if default is KeyError:
                raise KeyError(f"{path!r} not found in HDF5 LCM view") from None
            return default
        if isinstance(obj, h5py.Group):
            return LcmH5(obj, _root=self._file)
        return _decode_dataset(obj)

    def __getitem__(self, path):
        return self.get(path)

    def __contains__(self, path) -> bool:
        return self.get(path, _MISSING) is not _MISSING

    def keys(self, path="") -> List:
        node = self.get(path) if path else self
        if not isinstance(node, LcmH5):
            raise KeyError(f"{path!r} is not a table or list")
        if node.type == "list":
            return [int(k) for k in node._group]
        return [_h5_unescape(k) for k in node._group]

    def __len__(self) -> int:
        return len(self._group)

    def to_pylcm(self, name: Optional[str] = None) -> Pylcm:
        root = Pylcm(name=name or "Pylcm_object", type=self.type)
        for k in self._group:
            _read_into(root, k, self._group[k], level=1)
        return root

    def close(self) -> None:
        if self._file is not None and self._group is self._file:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


__all__ = ["Pylcm", "LcmChars", "LcmH5", "h5py_available"]
