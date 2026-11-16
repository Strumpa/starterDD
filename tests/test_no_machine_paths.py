# Guard: no machine-specific path (cluster account, local home folder) may be
# committed. Machine paths go in the site file (~/.config/starterDD/site.yaml,
# template docs/site.example.yaml) or in a git-ignored <STUDY>.local.yaml, and
# committed inputs reference them as ${NAME}.
import re
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"

#: Home / cluster scratch prefixes. /home/user/ is the generic SALOME
#: container path of the glow tutorials, not a machine path.
MACHINE_PATH_RE = re.compile(
    r"(?<![\w.])(/home/(?!user/)[\w.-]+|/Users/[\w.-]+|/scratch/|/lustre/|/gpfs/|"
    r"[A-Za-z]:\\\\?Users\\\\?)")
TEXT_SUFFIXES = {".py", ".yaml", ".yml", ".json", ".md", ".rst", ".txt", ".c2m",
                 ".x2m", ".cfg", ".toml", ".ini", ".sh", ".serp"}


def _tracked_files():
    try:
        out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO, capture_output=True,
                             check=True).stdout.decode()
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    return [REPO / f for f in out.split("\0") if f]


def test_no_machine_path_in_tracked_text_files():
    offenders = []
    for f in _tracked_files():
        if f.suffix not in TEXT_SUFFIXES or not f.is_file() or f.name == Path(__file__).name:
            continue
        try:
            text = f.read_text()
        except UnicodeDecodeError:
            continue
        for n, line in enumerate(text.splitlines(), 1):
            if MACHINE_PATH_RE.search(line):
                offenders.append(f"{f.relative_to(REPO)}:{n}: {line.strip()[:120]}")
    assert not offenders, ("Machine-specific paths are committed; move them to the site "
                           "file / <STUDY>.local.yaml and use ${NAME}:\n" + "\n".join(offenders))


def _absolute_values(obj, where=""):
    if isinstance(obj, str):
        s = obj.strip()
        if s.startswith(("/", "~")) or re.match(r"^[A-Za-z]:[\\/]", s):
            yield where, s
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from _absolute_values(v, f"{where}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _absolute_values(v, f"{where}[{i}]")


def test_data_yaml_files_hold_no_absolute_path():
    """Inputs under data/ are portable: relative paths or ${NAME} only."""
    offenders = []
    for f in sorted(DATA.rglob("*.y*ml")):
        if f.name.endswith(".local.yaml"):
            continue
        for doc in yaml.safe_load_all(f.read_text()):
            for where, value in _absolute_values(doc):
                offenders.append(f"{f.relative_to(REPO)}: {where} = {value}")
    assert not offenders, "\n".join(offenders)


def test_machine_path_pattern():
    assert MACHINE_PATH_RE.search("GE14_WORK: /home/someone/work")
    assert MACHINE_PATH_RE.search('"/Users/me/x"')
    assert not MACHINE_PATH_RE.search("PROJECT_ROOT = Path('/home/user/data/starterDD')")
    assert not MACHINE_PATH_RE.search("acelib: ${SERPENT_XS}/endfb8r1.xsdata")
