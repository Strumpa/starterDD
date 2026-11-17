# -*- coding: utf-8 -*-
"""Run Serpent2 on an input deck with OpenMP threads (plan P6).

The Python counterpart of the ``rserpent_para`` shell script::

    sss2 -omp <num_threads> <input>  > <input>.log

Differences with the script:

- Serpent runs **in the folder of the input** (no copy into a ``<machine>/``
  folder), so ``<input>_res.m`` / ``_det0.m`` land next to the input, where
  ``SerpentReference`` and ``SerpentRun`` look for them.
- No ``xs`` link: the nuclear data paths come from the ``set acelib`` /
  ``declib`` / ``nfylib`` cards of the deck. They are checked **before** the
  run, so a wrong path fails at once instead of in a cluster job.
- stdout and stderr both go to ``<input>.log``; the return code, the wall
  time and the produced outputs decide success, and ``IMP_KEFF`` is read
  from ``_res.m``.

Executable: the ``executable`` argument, else ``$SERPENT_EXEC``, else
``$sss2``, else ``sss2`` on the ``PATH``.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

PathLike = Union[str, Path]

#: Environment variables naming the Serpent2 executable, in order.
SERPENT_EXEC_ENV = ("SERPENT_EXEC", "sss2")
#: Executable name looked up on the PATH last.
SERPENT_EXEC_NAME = "sss2"
#: Suffix of the log (stdout + stderr) written next to the input.
LOG_SUFFIX = ".log"
#: Outputs removed before a forced re-run (Serpent appends to some of them).
RUN_OUTPUT_SUFFIXES = ("_res.m", "_det0.m", ".out", ".seed", LOG_SUFFIX)

_DATA_CARD_RE = re.compile(r'^\s*set\s+(acelib|declib|nfylib)\s+(.+)$', re.IGNORECASE)


class SerpentExecutableError(FileNotFoundError):
    """No usable Serpent2 executable."""


@dataclass
class SerpentRunResult:
    """Outcome of one Serpent2 run (or dry run)."""

    input_path: Path
    status: str                      # completed | failed | skipped | dry_run | not_ready
    return_code: Optional[int] = None
    wall_time_seconds: float = 0.0
    log_path: Optional[Path] = None
    keff: Optional[Tuple[float, float]] = None     # IMP_KEFF, relative 1-sigma
    executable: Optional[str] = None
    num_threads: Optional[int] = None
    message: str = ""
    problems: List[str] = field(default_factory=list)

    @property
    def success(self) -> bool:
        return self.status == "completed"

    def to_dict(self) -> dict:
        """Entry recorded in ``reference.yaml`` (``runs:``)."""
        return {
            "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
            "status": self.status, "return_code": self.return_code,
            "wall_time_seconds": float(self.wall_time_seconds),
            "executable": self.executable, "num_threads": self.num_threads,
            "keff": None if self.keff is None else [float(self.keff[0]), float(self.keff[1])],
            "log": None if self.log_path is None else self.log_path.name,
            "message": self.message,
        }


def resolve_serpent_executable(explicit: Optional[PathLike] = None) -> str:
    """The Serpent2 executable: *explicit*, ``$SERPENT_EXEC``, ``$sss2``,
    then ``sss2`` on the PATH. Raises :class:`SerpentExecutableError`."""
    tried = []
    candidates = [("argument / study file", explicit)] + \
                 [(f"${n}", os.environ.get(n)) for n in SERPENT_EXEC_ENV]
    for origin, c in candidates:
        if not c:
            continue
        p = Path(str(c)).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p.resolve())
        tried.append(f"{origin} = {c} (not an executable file)")
    found = shutil.which(SERPENT_EXEC_NAME)
    if found:
        return found
    tried.append(f"'{SERPENT_EXEC_NAME}' not on the PATH")
    raise SerpentExecutableError(
        "No Serpent2 executable: " + "; ".join(tried) + ". Give executable= "
        "(study file 'serpent: executable:' or ${SERPENT_EXEC} in the site file), "
        "or set $SERPENT_EXEC.")


def data_library_paths(input_path: PathLike) -> List[Tuple[str, str]]:
    """``[(card, path)]`` of the ``set acelib/declib/nfylib`` cards of a
    deck (``%`` comments stripped; an acelib card may list several files)."""
    out = []
    for line in Path(input_path).read_text(errors="replace").splitlines():
        line = line.split("%", 1)[0]
        m = _DATA_CARD_RE.match(line)
        if m:
            for word in shlex.split(m.group(2)):
                out.append((m.group(1).lower(), word))
    return out


def check_input(input_path: PathLike) -> List[str]:
    """Problems that would stop a run before it starts: missing input, or
    a nuclear data file of the deck that does not exist. Relative data paths
    are resolved from the input folder (Serpent also tries
    ``$SERPENT_DATA``)."""
    inp = Path(input_path)
    if not inp.is_file():
        return [f"input {inp} not found"]
    problems = []
    cards = data_library_paths(inp)
    if not any(c == "acelib" for c, _ in cards) and not os.environ.get("SERPENT_ACELIB"):
        problems.append("no 'set acelib' card and no $SERPENT_ACELIB")
    data_dir = os.environ.get("SERPENT_DATA")
    for card, path in cards:
        p = Path(path).expanduser()
        options = [p] if p.is_absolute() else [inp.parent / p] + (
            [Path(data_dir) / p] if data_dir else [])
        if not any(o.exists() for o in options):
            problems.append(f"set {card}: {path} not found")
    return problems


def _has_detectors(input_path: Path) -> bool:
    for line in input_path.read_text(errors="replace").splitlines():
        if re.match(r"^\s*det\s+\S+", line.split("%", 1)[0]):
            return True
    return False


class Serpent2Runner:
    """Runs ``sss2 -omp N <input>`` in the input's folder.

    Args:
        executable: Serpent2 executable (see :func:`resolve_serpent_executable`).
        num_threads: OpenMP threads (``-omp``); also sets ``OMP_NUM_THREADS``.
        extra_args: More command-line options, before the input name.
        env: Extra environment variables for the run.
        timeout: Seconds before the run is killed (``None``: no limit).
    """

    def __init__(self, executable: Optional[PathLike] = None, num_threads: int = 1,
                 extra_args: Sequence[str] = (), env: Optional[dict] = None,
                 timeout: Optional[float] = None):
        self._raw_executable = executable
        self.num_threads = int(num_threads)
        if self.num_threads < 1:
            raise ValueError("num_threads must be >= 1.")
        self.extra_args = [str(a) for a in extra_args]
        self.env = dict(env or {})
        self.timeout = timeout

    @property
    def executable(self) -> str:
        return resolve_serpent_executable(self._raw_executable)

    def command(self, input_path: PathLike) -> List[str]:
        return [self.executable, "-omp", str(self.num_threads), *self.extra_args,
                Path(input_path).name]

    def run(self, input_path: PathLike, dry_run: bool = False,
            force: bool = False) -> SerpentRunResult:
        """Run Serpent2 on *input_path*.

        An input that already has ``_res.m`` (and ``_det0.m`` when it has
        detectors) is skipped unless *force*; with *force* the old outputs
        are removed first (never symlinks: those belong to imported runs).
        *dry_run* checks the executable and the deck without running:
        ``dry_run`` when the run could start, ``not_ready`` otherwise (the
        reasons in ``message``).
        """
        inp = Path(input_path).resolve()
        base = inp
        res, det = Path(f"{base}_res.m"), Path(f"{base}_det0.m")
        log = Path(f"{base}{LOG_SUFFIX}")
        need_det = inp.is_file() and _has_detectors(inp)
        result = SerpentRunResult(inp, "failed", log_path=log, num_threads=self.num_threads)

        problems = check_input(inp)
        try:
            result.executable = self.executable
        except SerpentExecutableError as exc:
            problems.append(str(exc))
        result.problems = problems
        if problems:
            result.message = "; ".join(problems)
            if dry_run:
                result.status = "not_ready"
            return result

        done = res.exists() and (det.exists() or not need_det)
        if done and not force:
            result.status, result.message = "skipped", "has results (force=True to re-run)"
            result.log_path = log if log.exists() else None
            return result
        if dry_run:
            result.status = "dry_run"
            result.message = "would run: " + " ".join(shlex.quote(c) for c in self.command(inp))
            return result

        links = [p for suf in RUN_OUTPUT_SUFFIXES for p in [Path(f"{base}{suf}")] if p.is_symlink()]
        if links:
            result.message = (f"refusing to run: {[p.name for p in links]} are symlinks "
                              "(outputs of an imported run would be overwritten)")
            return result
        for suf in RUN_OUTPUT_SUFFIXES:
            Path(f"{base}{suf}").unlink(missing_ok=True)

        env = os.environ.copy()
        env.update({k: str(v) for k, v in self.env.items()})
        env["OMP_NUM_THREADS"] = str(self.num_threads)
        t0 = time.time()
        try:
            with open(log, "wb") as out:
                proc = subprocess.run(self.command(inp), cwd=inp.parent, env=env,
                                      stdout=out, stderr=subprocess.STDOUT,
                                      stdin=subprocess.DEVNULL, timeout=self.timeout)
            result.return_code = proc.returncode
        except subprocess.TimeoutExpired:
            result.message = f"killed after {self.timeout} s"
        result.wall_time_seconds = time.time() - t0

        missing = [p.name for p, needed in ((res, True), (det, need_det))
                   if needed and not p.exists()]
        if result.return_code == 0 and not missing:
            result.status = "completed"
            try:
                from ..PostTreatment.serpent_reader import read_res
                result.keff = read_res(res).keff("IMP_KEFF")
            except Exception as exc:                 # outputs exist; keff is a bonus
                result.message = f"IMP_KEFF not read: {exc}"
        elif not result.message:
            result.message = (f"return code {result.return_code}"
                              + (f", missing {missing}" if missing else "")
                              + f"; see {log.name}")
        return result


__all__ = ["Serpent2Runner", "SerpentRunResult", "SerpentExecutableError",
           "resolve_serpent_executable", "check_input", "data_library_paths",
           "SERPENT_EXEC_ENV", "LOG_SUFFIX"]
