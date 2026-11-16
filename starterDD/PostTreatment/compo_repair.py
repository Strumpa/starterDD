# -*- coding: utf-8 -*-
"""Diagnose and repair COMPOs whose global parameter names are permuted.

The bug (fixed 2026-11-15)
--------------------------
From the introduction of branch calculations (commit 358e43d, 2026-03-12) to
the fix, the generated ``EDIR.c2m`` read the statepoint values in the YAML
declaration order of ``CALCULATION_BRANCHES``::

    :: >>DCool<< >>TFuel<< >>TCool<< ;

while the main ``.x2m`` passed them in the loop order
(``BRANCH_LOOP_ORDER``, default ``fuel_temperature, coolant_temperature,
coolant_density``)::

    COMPO := EDIR ... :: <<name_compo>> <<TFuel>> <<TCool>> <<DCool>> ;

CLE-2000 binds procedure arguments by position, so each ``COMPO:`` call
stored the value of one parameter under another one's name. The transport
calculation itself was correct (MIX/LIB get their values from the main
loop); only the labels of the COMPO parameter tree are wrong.

Why relabelling is an exact repair
----------------------------------
In a COMPO, everything about a global parameter is indexed by its
*position* ``i``: ``GLOBAL/pval{i:08d}`` (sorted values), ``NVALUE(i)``,
``PARTYP(i)``/``PARFMT(i)``, and the tree records ``DEBARB``/``ARBVAL``
of every mixture. The name appears only in ``GLOBAL/PARKEY(i)``. DONJON
``NCR:`` (``NCRDRV.f``) and DRAGON ``COMPO:`` (``COMPO.f``) map a name to
its position through ``PARKEY`` and use the position afterwards. Renaming the
``PARKEY`` entries therefore puts every value back under its true name, with
no other change and no new calculation. All starterDD branch parameters are
``VALU REAL``, so their ``PARTYP``/``PARFMT`` are identical and need no change.

What relabelling cannot repair: outputs restricted to some statepoints
(selective ``CalculationOutput.state_points``). Their ``IF`` guards in
``EDIR.c2m`` tested the permuted variables, so those outputs were written at
the wrong statepoints or not at all. :func:`diagnose_procedures` reports them.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

from ._io import is_gz, read_bytes, read_text
from .pylcm import CT_CHAR, _HEADER_RE

PathLike = Union[str, Path]
NAME_WIDTH = 12


@dataclass
class ParameterDiagnosis:
    """Result of :func:`diagnose_procedures`.

    ``mapping[stored] = actual``: the values stored in the COMPO under the
    name ``stored`` are those of the parameter ``actual``.
    """

    reads: List[str]
    passes: List[str]
    guarded_outputs: List[str] = field(default_factory=list)

    @property
    def mapping(self) -> Dict[str, str]:
        return dict(zip(self.reads, self.passes))

    @property
    def permuted(self) -> bool:
        return self.reads != self.passes

    @property
    def repairable(self) -> bool:
        """True if :func:`relabel_compo_parameters` gives a fully correct
        COMPO: no permutation, or a permutation with no statepoint guard."""
        return not (self.permuted and self.guarded_outputs)

    def summary(self) -> str:
        if not self.permuted:
            return f"OK: EDIR reads {self.reads} in the order they are passed."
        moved = ", ".join(f"'{s}' holds {a}" for s, a in self.mapping.items() if s != a)
        txt = f"PERMUTED: {moved}."
        if self.guarded_outputs:
            txt += (f" Outputs with statepoint guards were evaluated on the wrong "
                    f"variables: {self.guarded_outputs}; relabelling cannot recover them.")
        return txt


def _main_x2m(procedures_dir: Path) -> Path:
    x2m = sorted(procedures_dir.glob("*.x2m"))
    if len(x2m) != 1:
        raise FileNotFoundError(f"Expected one .x2m in {procedures_dir}, found {len(x2m)}.")
    return x2m[0]


def diagnose_procedures(procedures_dir: PathLike, edir_name: str = "EDIR.c2m") -> ParameterDiagnosis:
    """Compare the statepoint arguments read by ``EDIR.c2m`` with those
    passed by the main ``.x2m`` (e.g. a run's archived ``procedures/``)."""
    d = Path(procedures_dir)
    edir = read_text(d / edir_name)
    x2m = read_text(_main_x2m(d))
    m = re.search(r"^REAL ([^;\n]*);\s*\n::((?:\s*>>\w+<<)+)\s*;", edir, re.M)
    reads = re.findall(r">>(\w+)<<", m.group(2)) if m else []
    stem = Path(edir_name).stem
    c = re.search(rf"{stem}\b[^;]*?COMPO\s*::\s*<<\w+>>([^;]*);", x2m, re.S)
    passes = re.findall(r"<<(\w+)>>", c.group(1)) if c else []
    if len(reads) != len(passes):
        raise ValueError(f"EDIR reads {reads} but the x2m passes {passes}.")
    guarded: List[str] = []
    current = None
    for line in edir.splitlines():
        h = re.match(r"\*\s+EDI: CALL FOR (\S+)", line)
        if h:
            current = h.group(1)
        elif current and re.match(r"\s*IF\s", line) and current not in guarded:
            guarded.append(current)
    return ParameterDiagnosis(reads, passes, guarded)


def scheme_parameter_mapping(scheme) -> Dict[str, str]:
    """``{stored: actual}`` implied by a scheme for procedures generated
    **before** the 2026-11-15 fix (YAML order read, loop order passed). Prefer
    :func:`diagnose_procedures` on the archived procedures when available."""
    reads = [b.para_name for b in scheme.branches]
    passes = [b.para_name for b in scheme.get_ordered_branches()]
    return dict(zip(reads, passes))


def _iter_global_parkey(buf: bytes):
    """Yield ``(text_start, text_end, n_words)`` byte spans of the text of
    every ``GLOBAL/PARKEY`` character record."""
    stack: List[tuple] = []          # (level, key) of open tables/lists
    for m in _HEADER_RE.finditer(buf):
        level, ntype, ctype, csize = (int(x) for x in m.group(1, 2, 3, 4))
        if level < 0:
            while stack and stack[-1][0] >= -level:
                stack.pop()
            continue
        while stack and stack[-1][0] >= level:
            stack.pop()
        line_end = buf.index(b"\n", m.end()) + 1
        key_end = buf.find(b"\n", line_end) + 1 if ntype > 0 else line_end
        key = buf[line_end:key_end].decode("ascii", "replace")[:max(ntype, 0)].strip()
        if ctype in (0, 10):                     # table / list
            stack.append((level, key))
            continue
        if ctype == CT_CHAR and key == "PARKEY" and stack and stack[-1][1] == "GLOBAL" and csize > 0:
            pos = key_end
            for _ in range(-(-csize // 8)):      # word-length lines
                pos = buf.index(b"\n", pos) + 1
            nxt = _HEADER_RE.search(buf, pos)
            yield pos, (nxt.start() if nxt else len(buf)), csize


def _rename_text(text: str, n_words: int, mapping: Mapping[str, str]) -> str:
    nchar = 4 * n_words
    lines = text.split("\n")
    raw = "".join(ln.rstrip("\r").ljust(80) if i < len(lines) - 1 else ln.rstrip("\r")
                  for i, ln in enumerate(lines))
    names = [raw[i:i + NAME_WIDTH] for i in range(0, nchar, NAME_WIDTH)]
    new = []
    for n in names:
        s = n.strip()
        new.append(mapping.get(s, s).ljust(NAME_WIDTH))
    new_raw = "".join(new)[:nchar]
    # Same line layout as the original (80 columns, original line lengths).
    out, pos = [], 0
    for ln in lines:
        width = len(ln.rstrip("\r"))
        out.append(new_raw[pos:pos + width] + ln[len(ln.rstrip("\r")):])
        pos += 80 if width == 80 else width
    return "\n".join(out)


def relabel_compo_parameters(src: PathLike, dst: PathLike, mapping: Mapping[str, str]) -> int:
    """Write *dst*, a copy of the ASCII COMPO *src* in which every
    ``GLOBAL/PARKEY`` entry ``stored`` is renamed ``mapping[stored]``.

    *mapping* must be a permutation of the names it touches. Only the
    ``PARKEY`` text bytes change (names keep their 12-character width), so
    *dst* has the same size as *src*. ``.gz`` is supported on both sides.
    Returns the number of ``PARKEY`` records rewritten (one per output).
    """
    mapping = {k: v for k, v in mapping.items() if k != v}
    if sorted(mapping) != sorted(mapping.values()):
        raise ValueError(f"mapping {mapping} is not a permutation of parameter names.")
    if any(len(n) > NAME_WIDTH for n in mapping):
        raise ValueError("parameter names are limited to 12 characters.")
    buf = read_bytes(src)
    pieces, last, n = [], 0, 0
    for start, end, n_words in _iter_global_parkey(buf):
        text = buf[start:end].decode("ascii")
        names = {text.replace("\n", "")[i:i + NAME_WIDTH].strip()
                 for i in range(0, 4 * n_words, NAME_WIDTH)}
        if mapping and not set(mapping) <= names:
            raise ValueError(f"PARKEY {sorted(names)} does not contain {sorted(mapping)}.")
        pieces += [buf[last:start], _rename_text(text, n_words, mapping).encode("ascii")]
        last, n = end, n + 1
    pieces.append(buf[last:])
    out = b"".join(pieces)
    assert len(out) == len(buf)
    if is_gz(dst):
        with gzip.open(dst, "wb") as f:
            f.write(out)
    else:
        Path(dst).write_bytes(out)
    return n


__all__ = ["ParameterDiagnosis", "diagnose_procedures", "scheme_parameter_mapping",
           "relabel_compo_parameters"]
