# -*- coding: utf-8 -*-
"""DRAGON5 / Serpent2 code-to-code studies (plan §3.9, phase P5).

A :class:`CodeToCodeStudy` is the product

    assemblies × scheme variants × operating conditions

- an **assembly** is a geometry + materials + TDT folder + base scheme (the
  base schemes of GE14 DOM, DOM-C and VAN differ by their meshing);
- a **variant** is a set of step changes applied to every assembly's base
  scheme (``DragonCalculationScheme.set_step``), e.g. another SSH macro
  grouping or another L1 method;
- a **condition** (``HZP``, ``HFP``) is a set of branch values: the coolant
  statepoints with their exact benchmark isotopics
  (``with_coolant_statepoints``), the fuel and coolant temperatures.

Each (assembly, variant, condition) is **one DRAGON case**, i.e. one
CLE-2000 procedure set and one multi-statepoint COMPO. Each statepoint gets a
label ``<condition>_<coolant label>`` (``HFP_40``) that names its Serpent2
reference, which is shared by every variant::

    <root>/<study>/
        study_manifest.yaml        rebuilt from the case records by every command
        study_config.yaml          the study definition (absolute paths)
        references/<assembly>/<library>/<sp_label>/     Serpent2 references
        <assembly>/dragon/<variant>/<condition>/
            case.yaml              case record (status, runs, statepoints)
            scheme.yaml            resolved scheme
            procedures/            generated CLE-2000 files
            <timestamp>/  latest   DragonRunner run folders
        post/                      study.h5, summary.csv, summary.json

Each case writes **only its own** ``case.yaml``; ``study_manifest.yaml`` is an
aggregate rewritten atomically, so several cluster jobs can run different
cases of one study at the same time. A completed case is skipped unless
``force=True`` or its inputs changed (``fingerprint``: resolved scheme hash,
geometry and material file hashes, library, case name).

The study writes Serpent2 inputs; it does not run Serpent2 (P6). Existing
Serpent2 runs are registered with :meth:`CodeToCodeStudy.import_serpent_references`,
existing DRAGON COMPOs with :meth:`CodeToCodeStudy.import_dragon_run` (no
re-run; ``<case>/imported_<time>/``, ``origin: imported``).

Post-processing (:meth:`CodeToCodeStudy.post_process`) writes ``post/study.h5``;
figures are made from that file alone (:mod:`.study_figures`).

Command line (one case per cluster job)::

    python -m starterDD.Workflow STUDY.yaml status
    python -m starterDD.Workflow STUDY.yaml run --assembly GE14_DOM --condition HFP --threads 10
    python -m starterDD.Workflow STUDY.yaml import-dragon     # dragon: import: block
    python -m starterDD.Workflow STUDY.yaml post
    python -m starterDD.Workflow STUDY.yaml figures
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import itertools
import json
import os
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import yaml

from ..DDModel.DragonCalculationScheme import DragonCalculationScheme

PathLike = Union[str, Path]

MANIFEST_FILE = "study_manifest.yaml"
CONFIG_FILE = "study_config.yaml"
CASE_FILE = "case.yaml"
SCHEME_FILE = "scheme.yaml"
PROCEDURES_DIR = "procedures"

#: Draglib file → CLE-2000 alias, per library (the GE14 runs' names).
DEFAULT_DRAGLIBS = {
    "endfb8r1": {"draglibendfb8r1SHEM295_v5p1": "endfb8r1295"},
    "jeff311": {"draglibJeff3p1p1SHEM295_v5p1": "J311_295"},
}

#: Label prefix of a multi-valued branch without labels (``TF900``).
_BRANCH_SHORT = {"fuel_temperature": "TF", "coolant_temperature": "TC",
                 "coolant_density": "D", "burnup": "BU"}
_BRANCH_NAMES = {"fuel_temperature": "Fuel temperature variation",
                 "coolant_temperature": "Coolant temperature variation"}
_LABEL_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")
_VAR_RE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")


def _check_label(kind: str, label: str) -> str:
    label = str(label)
    if not _LABEL_RE.match(label):
        raise ValueError(f"{kind} '{label}' is used in folder names: use letters, digits "
                         "and _ - . + only.")
    return label


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _sha256(path: PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _atomic_write_yaml(path: Path, data: Any) -> None:
    """Write YAML through a temporary file and ``os.replace``, so readers
    never see a half-written file (concurrent cluster jobs)."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(yaml.safe_dump(data, sort_keys=False, default_flow_style=None, width=100))
    os.replace(tmp, path)


def _rel(path: Optional[PathLike], base: Path) -> Optional[str]:
    if path is None:
        return None
    try:
        return os.path.relpath(path, base)
    except ValueError:          # other drive (Windows)
        return str(path)


def _plain(obj):
    """YAML-safe copy (numpy scalars, tuples, paths)."""
    if isinstance(obj, Mapping):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        return obj.item()
    return obj


# ─────────────────────────────────────────────────────────────────────────────
# Study building blocks
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class StudyAssembly:
    """One lattice of the study.

    Args:
        name: Label (folder names, Serpent reference ``assembly``), e.g.
            ``"GE14_DOM"``.
        geometry: GEOM YAML (``DragonCase`` ``GEOM``).
        materials: MATS YAML. The coolant isotopics of each statepoint come
            from the condition, the rest from this file.
        tdt_path: Folder holding the TDT files of this lattice.
        scheme: Base scheme (YAML path or object). ``None``: the study's.
        case_name: DRAGON case name (x2m, ``_CPO_<case_name>``, TDT prefix).
            Default: *name*.
        tdt_base_name: TDT file prefix if it differs from *case_name*.
        steps: ``{step: {yaml_key: value}}`` changes applied to the base
            scheme of this lattice only (``set_step``).
    """

    name: str
    geometry: PathLike
    materials: PathLike
    tdt_path: PathLike
    scheme: Optional[Union[PathLike, DragonCalculationScheme]] = None
    case_name: Optional[str] = None
    tdt_base_name: Optional[str] = None
    steps: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self):
        self.name = _check_label("Assembly", self.name)
        self.geometry = Path(self.geometry)
        self.materials = Path(self.materials)
        self.tdt_path = Path(self.tdt_path)
        if self.scheme is not None and not isinstance(self.scheme, DragonCalculationScheme):
            self.scheme = Path(self.scheme)
        self.case_name = self.case_name or self.name
        self.steps = {k: dict(v) for k, v in (self.steps or {}).items()}

    def to_dict(self) -> dict:
        d = {"geometry": str(self.geometry), "materials": str(self.materials),
             "tdt_path": str(self.tdt_path), "case_name": self.case_name}
        if isinstance(self.scheme, Path):
            d["scheme"] = str(self.scheme)
        elif self.scheme is not None:
            d["scheme"] = f"<object: {self.scheme.name}>"
        if self.tdt_base_name:
            d["tdt_base_name"] = self.tdt_base_name
        if self.steps:
            d["steps"] = _plain(self.steps)
        return d


@dataclass
class StudyVariant:
    """Changes of the scheme steps shared by every assembly.

    Args:
        name: Folder label, e.g. ``"2L_PT3x3_IC3x3_6x6"``.
        steps: ``{step: {yaml_key: value}}`` (``set_step``); empty = the
            base schemes as they are.
        description: Free text.
    """

    name: str
    steps: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    description: str = ""

    def __post_init__(self):
        self.name = _check_label("Variant", self.name)
        self.steps = {k: dict(v) for k, v in (self.steps or {}).items()}

    def to_dict(self) -> dict:
        d: Dict[str, Any] = {"description": self.description}
        if self.steps:
            d["steps"] = _plain(self.steps)
        return d


@dataclass
class StudyCondition:
    """Operating condition: the branch values of one DRAGON case.

    Args:
        name: Label, e.g. ``"HFP"``; it prefixes the statepoint labels.
        coolant: ``{label: source}``, a material YAML (its ``COOLANT``
            entry) or a composition dict, one per coolant statepoint
            (``with_coolant_statepoints``). ``None``: keep the scheme's
            coolant density branch.
        fuel_temperature: Fuel temperatures (K). ``None``: keep the scheme's.
        coolant_temperature: Coolant (and moderator) temperatures (K).
        densities: ``{label: rho}`` PARA values, default: from the isotopics.
        description: Free text.
    """

    name: str
    coolant: Optional[Mapping[str, Union[PathLike, Mapping[str, float]]]] = None
    fuel_temperature: Optional[Sequence[float]] = None
    coolant_temperature: Optional[Sequence[float]] = None
    densities: Optional[Mapping[str, float]] = None
    description: str = ""

    def __post_init__(self):
        self.name = _check_label("Condition", self.name)
        if self.coolant is not None:
            self.coolant = {_check_label("Coolant statepoint", k):
                            (v if isinstance(v, Mapping) else str(v))
                            for k, v in self.coolant.items()}
        for attr in ("fuel_temperature", "coolant_temperature"):
            v = getattr(self, attr)
            if v is not None:
                setattr(self, attr, [float(x) for x in (v if isinstance(v, (list, tuple)) else [v])])

    def apply(self, scheme: DragonCalculationScheme) -> DragonCalculationScheme:
        """Copy of *scheme* with this condition's branch values."""
        s = scheme
        if self.coolant is not None:
            s = s.with_coolant_statepoints(dict(self.coolant), densities=self.densities)
        for btype in ("fuel_temperature", "coolant_temperature"):
            vals = getattr(self, btype)
            if vals is not None:
                s = s.with_branch_values(btype, vals, name=_BRANCH_NAMES[btype])
        return s if s is not scheme else scheme.copy()

    def to_dict(self) -> dict:
        d: Dict[str, Any] = {"description": self.description}
        if self.coolant is not None:
            d["coolant"] = _plain(self.coolant)
        for attr in ("fuel_temperature", "coolant_temperature", "densities"):
            if getattr(self, attr) is not None:
                d[attr] = _plain(getattr(self, attr))
        return d


@dataclass
class StudyStatepoint:
    """One statepoint of a case."""

    label: str
    values: Dict[str, float]           # {branch_type: value}
    params: Dict[str, float]           # {COMPO PARA name: value}

    def to_dict(self) -> dict:
        return {"label": self.label, "values": _plain(self.values),
                "compo_parameters": _plain(self.params)}


def case_statepoints(scheme: DragonCalculationScheme, prefix: str) -> List[StudyStatepoint]:
    """Statepoints of the branch grid of *scheme* (loop order), labelled
    ``<prefix>_<branch labels>``. A labelled branch contributes its label; an
    unlabelled one contributes ``TF900``-like parts only if it has several
    values."""
    branches = scheme.get_ordered_branches()
    out = []
    for idx in itertools.product(*[range(len(b.values)) for b in branches]):
        values, params, parts = {}, {}, [prefix]
        for b, i in zip(branches, idx):
            v = b.values[i]
            values[b.type] = v
            params[b.para_name] = v
            if b.labels:
                parts.append(b.labels[i])
            elif len(b.values) > 1:
                parts.append(f"{_BRANCH_SHORT.get(b.type, b.type)}{v:g}")
        out.append(StudyStatepoint(_check_label("Statepoint", "_".join(parts)), values, params))
    if len({sp.label for sp in out}) != len(out):
        raise ValueError(f"Statepoint labels are not unique: {[sp.label for sp in out]}")
    return out


@dataclass
class StudyCase:
    """One DRAGON case: (assembly, variant, condition) and its resolved scheme."""

    study: "CodeToCodeStudy" = field(repr=False)
    assembly: StudyAssembly
    variant: StudyVariant
    condition: StudyCondition
    scheme: DragonCalculationScheme = field(repr=False)

    @property
    def key(self) -> str:
        return f"{self.assembly.name}/{self.variant.name}/{self.condition.name}"

    @property
    def case_name(self) -> str:
        return self.assembly.case_name

    @property
    def directory(self) -> Path:
        return (self.study.directory / self.assembly.name / "dragon"
                / self.variant.name / self.condition.name)

    @property
    def statepoints(self) -> List[StudyStatepoint]:
        return case_statepoints(self.scheme, self.condition.name)

    @property
    def fingerprint(self) -> str:
        """Hash of everything that changes the DRAGON results of this case."""
        d = {"scheme": self.scheme.content_hash(16),
             "geometry": self.study._file_hash(self.assembly.geometry),
             "materials": self.study._file_hash(self.assembly.materials),
             "case_name": self.case_name, "tdt_base_name": self.assembly.tdt_base_name,
             "draglib": self.study.draglib,
             "structural_temperature": self.study.structural_temperature}
        blob = json.dumps(d, sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:16]

    # ── record ─────────────────────────────────────────────────────────────
    @property
    def record_path(self) -> Path:
        return self.directory / CASE_FILE

    @property
    def record(self) -> dict:
        p = self.record_path
        return (yaml.safe_load(p.read_text()) or {}) if p.exists() else {}

    def write_record(self, rec: dict) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        rec["updated"] = _now()
        _atomic_write_yaml(self.record_path, _plain(rec))

    @property
    def status(self) -> str:
        """``planned``, ``generated``, ``dry_run``, ``completed`` or
        ``failed``; ``stale`` when the inputs changed since then."""
        rec = self.record
        if not rec:
            return "planned"
        if rec.get("fingerprint") != self.fingerprint:
            return "stale"
        return rec.get("status", "planned")

    @property
    def run_directory(self) -> Optional[Path]:
        """Folder of the last completed DRAGON run (``None`` if none)."""
        cur = self.record.get("current_run")
        return None if cur is None else self.directory / cur

    @property
    def compo_path(self) -> Optional[Path]:
        run = self.run_directory
        if run is None:
            return None
        p = run / f"_CPO_{self.case_name}"
        if p.exists():
            return p
        hits = sorted(q for q in run.glob("_CPO_*") if q.suffix != ".h5")
        return hits[0] if len(hits) == 1 else None

    def reference_folder(self, sp: StudyStatepoint) -> Path:
        from ..PostTreatment.references import SerpentReference
        return SerpentReference.folder(self.study.references_root, self.assembly.name,
                                       self.study.library, sp.label)


# ─────────────────────────────────────────────────────────────────────────────
# Config helpers
# ─────────────────────────────────────────────────────────────────────────────


def _expand(obj, variables: Mapping[str, str], where: str = ""):
    """Substitute ``${NAME}`` in every string of *obj* (recursively)."""
    if isinstance(obj, str):
        def sub(m):
            name = m.group(1)
            if name not in variables:
                raise KeyError(f"Undefined variable ${{{name}}}{' in ' + where if where else ''}; "
                               "define it under 'variables:', in the environment, or pass "
                               "variables={...}.")
            return str(variables[name])
        for _ in range(10):           # variables may reference variables
            new = _VAR_RE.sub(sub, obj)
            if new == obj:
                return new
            obj = new
        raise ValueError(f"Variable expansion does not terminate: {obj}")
    if isinstance(obj, Mapping):
        return {k: _expand(v, variables, where or str(k)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand(v, variables, where) for v in obj]
    return obj


#: Spec keys renamed on 2026-11-15 -> their new name.
_RENAMED_SPEC_KEYS = {"keff_isotopes": "fuel_isotopes_inventory",
                      "pins_fine": "serpent_fine_pins"}
_SPEC_KEYS = ("preset", "isotopes", "fine_mesh", "fuel_isotopes_inventory",
              "serpent_fine_pins", "global_totals")


def spec_from_config(cfg):
    """``ComparisonSpec`` from a config value: ``None`` (keep the scheme
    outputs), a spec, or ``{preset, isotopes, fine_mesh,
    fuel_isotopes_inventory, serpent_fine_pins, global_totals}``
    (:meth:`ComparisonSpec.preset`)."""
    from ..PostTreatment.spec import ComparisonSpec
    if cfg is None or cfg == "from_scheme":
        return None
    if isinstance(cfg, ComparisonSpec):
        return cfg
    cfg = dict(cfg)
    old = [k for k in cfg if k in _RENAMED_SPEC_KEYS]
    if old:
        raise KeyError("Renamed spec key(s): " + ", ".join(
            f"'{k}' -> '{_RENAMED_SPEC_KEYS[k]}'" for k in old) + ".")
    unknown = set(cfg) - set(_SPEC_KEYS)
    if unknown:
        raise KeyError(f"Unknown spec key(s) {sorted(unknown)}; known {list(_SPEC_KEYS)}.")
    preset = cfg.pop("preset", "integrated")
    if "isotopes" not in cfg:
        raise KeyError("A spec needs 'isotopes'.")
    isotopes = cfg.pop("isotopes")
    return ComparisonSpec.preset(preset, isotopes, **cfg)


# ─────────────────────────────────────────────────────────────────────────────
# Study
# ─────────────────────────────────────────────────────────────────────────────


class CodeToCodeStudy:
    """DRAGON5 cases over assemblies × variants × conditions, their Serpent2
    references, the manifest and the post-processing.

    Args:
        name: Study label (folder name).
        root: Parent folder of the study folder.
        assemblies: ``{name: StudyAssembly or dict of its arguments}``.
        conditions: ``{name: StudyCondition or dict}``.
        variants: ``{name: StudyVariant or dict}``; default one variant
            ``"base"`` without changes.
        scheme: Default base scheme for assemblies without their own.
        spec: ``ComparisonSpec`` (or its config dict) that sets the DRAGON
            outputs of every case (``spec.apply_to``). ``None``: keep the
            schemes' ``CALCULATION_OUTPUTS``.
        reference_spec: Spec of the Serpent2 detectors. Default: the union of
            the case specs. Give a finer one (e.g. 295g) to make references
            reusable by 1L variants.
        library: Nuclear-data evaluation (Serpent suffixes, Q-values,
            reference folders).
        draglib: ``{draglib file: CLE-2000 alias}``; default
            :data:`DEFAULT_DRAGLIBS` of *library*.
        draglib_paths: ``{draglib file: path}``; default ``$DRAGLIB_DIR``...
            (``DragonRunner``).
        dragon_executable: Default ``$dragon_exec`` / ``$DRAGON_EXEC``.
        num_threads: OpenMP threads of each DRAGON run.
        references_root: Serpent reference registry, default
            ``<study>/references``. Several studies may share one.
        serpent_settings: ``S2_Settings`` attributes of generated references
            (``acelib``, ``neutrons_per_cycle``, ``active_cycles``...).
        structural_temperature: Clad, box and other structures (K).
        post_processing: Defaults of :meth:`post_process`.
        source: YAML file the study was read from (recorded only).
    """

    def __init__(self, name: str, root: PathLike,
                 assemblies: Mapping[str, Union[StudyAssembly, Mapping]],
                 conditions: Mapping[str, Union[StudyCondition, Mapping]],
                 variants: Optional[Mapping[str, Union[StudyVariant, Mapping]]] = None,
                 scheme: Optional[Union[PathLike, DragonCalculationScheme]] = None,
                 spec=None, reference_spec=None, library: str = "endfb8r1",
                 draglib: Optional[Mapping[str, str]] = None,
                 draglib_paths: Optional[Mapping[str, PathLike]] = None,
                 dragon_executable: Optional[PathLike] = None, num_threads: int = 1,
                 references_root: Optional[PathLike] = None,
                 serpent_settings: Optional[Mapping[str, Any]] = None,
                 structural_temperature: float = 600.0,
                 post_processing: Optional[Mapping[str, Any]] = None,
                 source: Optional[PathLike] = None):
        self.name = _check_label("Study", name)
        self.root = Path(root).expanduser().resolve()
        self.assemblies = self._build(assemblies, StudyAssembly)
        self.conditions = self._build(conditions, StudyCondition)
        self.variants = self._build(variants or {"base": {}}, StudyVariant)
        if not self.assemblies or not self.conditions:
            raise ValueError("A study needs at least one assembly and one condition.")
        self.scheme = (Path(scheme) if scheme is not None
                       and not isinstance(scheme, DragonCalculationScheme) else scheme)
        self.spec_config = spec if isinstance(spec, Mapping) else None
        self.reference_spec_config = reference_spec if isinstance(reference_spec, Mapping) else None
        self.spec = spec_from_config(spec)
        self._reference_spec = spec_from_config(reference_spec)
        self.library = library
        if draglib is None:
            if library not in DEFAULT_DRAGLIBS:
                raise KeyError(f"No default draglib for library '{library}'; pass draglib=.")
            draglib = DEFAULT_DRAGLIBS[library]
        self.draglib = dict(draglib)
        self.draglib_paths = {k: str(v) for k, v in (draglib_paths or {}).items()}
        self.dragon_executable = None if dragon_executable is None else str(dragon_executable)
        self.num_threads = int(num_threads)
        self.references_root = (Path(references_root).expanduser().resolve()
                                if references_root else self.directory / "references")
        self.serpent_settings = dict(serpent_settings or {})
        self.structural_temperature = float(structural_temperature)
        self.post_processing = dict(post_processing or {})
        self.source = None if source is None else Path(source).resolve()
        #: ``{assembly: {label: .serp}}`` of the study file's ``serpent: import:``.
        self.serpent_imports: Dict[str, Dict[str, Path]] = {}
        #: ``{"<assembly>/<variant>/<condition>": COMPO or run folder}`` of the
        #: study file's ``dragon: import:`` (:meth:`import_dragon_runs`).
        self.dragon_imports: Dict[str, Path] = {}
        self._schemes: Dict[str, DragonCalculationScheme] = {}
        self._cases: Dict[str, StudyCase] = {}
        self._hashes: Dict[str, str] = {}
        self._lattices: Dict[str, Any] = {}

    @staticmethod
    def _build(items, cls):
        out = {}
        for k, v in (items or {}).items():
            obj = v if isinstance(v, cls) else cls(name=k, **dict(v or {}))
            if obj.name != k:
                raise ValueError(f"{cls.__name__} key '{k}' differs from its name '{obj.name}'.")
            out[k] = obj
        return out

    # ── construction from YAML ─────────────────────────────────────────────
    @classmethod
    def from_yaml(cls, path: PathLike, variables: Optional[Mapping[str, str]] = None,
                  **overrides) -> "CodeToCodeStudy":
        """Read a study file (see ``data/BWRProgressionProblems/GE14_inputs/
        STUDIES`` for an example).

        ``${NAME}`` is replaced by, in order of precedence, *variables*, the
        environment, then the file's ``variables:`` block. Relative paths are
        relative to the YAML file. *overrides* replace constructor arguments
        (e.g. ``root=``, ``num_threads=``).
        """
        path = Path(path).resolve()
        raw = yaml.safe_load(path.read_text()) or {}
        file_vars = {k: str(v) for k, v in (raw.pop("variables", None) or {}).items()}
        env = {k: v for k, v in os.environ.items() if _VAR_RE.fullmatch("${" + k + "}")}
        allvars = {**file_vars, **env, **{k: str(v) for k, v in (variables or {}).items()}}
        cfg = _expand(raw, allvars)
        base = path.parent

        def p(value):
            if value is None:
                return None
            q = Path(str(value)).expanduser()
            return q if q.is_absolute() else (base / q).resolve()

        known = {"study", "dragon", "serpent", "spec", "reference_spec", "scheme",
                 "assemblies", "variants", "conditions", "post_processing"}
        unknown = set(cfg) - known
        if unknown:
            raise KeyError(f"{path.name}: unknown section(s) {sorted(unknown)}; known {sorted(known)}.")
        st = dict(cfg.get("study") or {})
        dr = dict(cfg.get("dragon") or {})
        se = dict(cfg.get("serpent") or {})

        assemblies = {}
        for name, a in (cfg.get("assemblies") or {}).items():
            a = dict(a)
            for k in ("geometry", "materials", "tdt_path", "scheme"):
                if k in a:
                    a[k] = p(a[k])
            assemblies[name] = a
        conditions = {}
        for name, c in (cfg.get("conditions") or {}).items():
            c = dict(c)
            if c.get("coolant") is not None:
                c["coolant"] = {str(k): (v if isinstance(v, Mapping) else p(v))
                                for k, v in c["coolant"].items()}
            conditions[name] = c

        draglib_paths = dict(dr.get("draglib_paths") or {})
        if dr.get("draglib_dir"):
            for lib in (dr.get("draglib") or DEFAULT_DRAGLIBS.get(st.get("library", "endfb8r1"), {})):
                draglib_paths.setdefault(lib, p(dr["draglib_dir"]) / lib)

        kwargs = dict(
            name=st.get("name", path.stem), root=p(st.get("root", ".")),
            assemblies=assemblies, conditions=conditions, variants=cfg.get("variants"),
            scheme=p(cfg.get("scheme")), spec=cfg.get("spec"),
            reference_spec=cfg.get("reference_spec"),
            library=st.get("library", "endfb8r1"), draglib=dr.get("draglib"),
            draglib_paths={k: p(v) for k, v in draglib_paths.items()},
            dragon_executable=p(dr.get("executable")),
            num_threads=dr.get("num_threads", 1),
            references_root=p(se.get("references_root")),
            serpent_settings=se.get("settings"),
            structural_temperature=dr.get("structural_temperature", 600.0),
            post_processing=cfg.get("post_processing"), source=path)
        kwargs.update(overrides)
        study = cls(**kwargs)
        study.serpent_imports = {a: {str(lab): p(f) for lab, f in sps.items()}
                                 for a, sps in (se.get("import") or {}).items()}
        study.dragon_imports = {str(k): p(v) for k, v in (dr.get("import") or {}).items()}
        return study

    # ── paths, cached inputs ───────────────────────────────────────────────
    @property
    def directory(self) -> Path:
        return self.root / self.name

    @property
    def manifest_path(self) -> Path:
        return self.directory / MANIFEST_FILE

    @property
    def registry(self):
        from ..PostTreatment.references import ReferenceRegistry
        return ReferenceRegistry(self.references_root)

    def _file_hash(self, path: PathLike) -> str:
        key = str(Path(path).resolve())
        if key not in self._hashes:
            self._hashes[key] = _sha256(path)[:16]
        return self._hashes[key]

    def _base_scheme(self, asm: StudyAssembly) -> DragonCalculationScheme:
        src = asm.scheme if asm.scheme is not None else self.scheme
        if src is None:
            raise ValueError(f"Assembly '{asm.name}' has no scheme and the study has no default.")
        if isinstance(src, DragonCalculationScheme):
            return src
        key = str(src)
        if key not in self._schemes:
            self._schemes[key] = DragonCalculationScheme.from_yaml(key)
        return self._schemes[key]

    def lattice(self, assembly: Union[str, StudyAssembly]):
        """``PinLattice`` of an assembly (pin numbering of the DRAGON case)."""
        from ..PostTreatment.lattice_maps import PinLattice
        asm = self.assemblies[assembly] if isinstance(assembly, str) else assembly
        if asm.name not in self._lattices:
            self._lattices[asm.name] = PinLattice.from_inputs(asm.geometry, asm.materials,
                                                              name=asm.case_name)
        return self._lattices[asm.name]

    # ── cases ──────────────────────────────────────────────────────────────
    def _resolve_scheme(self, asm: StudyAssembly, var: StudyVariant,
                        cond: StudyCondition) -> DragonCalculationScheme:
        s = self._base_scheme(asm).copy()
        for step, params in asm.steps.items():
            s.set_step(step, **params)
        for step, params in var.steps.items():
            s.set_step(step, **params)
        s = cond.apply(s)
        if self.spec is not None:
            s = self.spec.apply_to(s)
        return s.validate()

    def case(self, assembly: str, variant: str, condition: str) -> StudyCase:
        key = f"{assembly}/{variant}/{condition}"
        if key not in self._cases:
            asm, var, cond = (self.assemblies[assembly], self.variants[variant],
                              self.conditions[condition])
            self._cases[key] = StudyCase(self, asm, var, cond,
                                         self._resolve_scheme(asm, var, cond))
        return self._cases[key]

    @staticmethod
    def _pick(names: Iterable[str], selected) -> List[str]:
        names = list(names)
        if selected is None:
            return names
        selected = [selected] if isinstance(selected, str) else list(selected)
        unknown = set(selected) - set(names)
        if unknown:
            raise KeyError(f"Unknown {sorted(unknown)}; known {names}.")
        return [n for n in names if n in selected]

    def cases(self, assemblies=None, variants=None, conditions=None) -> List[StudyCase]:
        """Cases in assembly, variant, condition order (optionally filtered by
        name or list of names)."""
        return [self.case(a, v, c)
                for a in self._pick(self.assemblies, assemblies)
                for v in self._pick(self.variants, variants)
                for c in self._pick(self.conditions, conditions)]

    def case_spec(self, case: StudyCase, resolve_pins: bool = True):
        """Observables of one case: the study spec, or the scheme outputs,
        with pin selectors (``serpent_fine_pins: diagonal``) resolved on the
        assembly lattice unless *resolve_pins* is False."""
        from ..PostTreatment.spec import ComparisonSpec
        spec = self.spec if self.spec is not None else ComparisonSpec.from_scheme(case.scheme)
        if resolve_pins and spec.pins_unresolved:
            spec = spec.resolve_pins(self.lattice(case.assembly))
        return spec

    def reference_spec(self):
        """Spec of the Serpent detectors: explicit, or the union of the case
        specs (a reference must serve every variant)."""
        from ..PostTreatment.spec import ComparisonSpec
        if self._reference_spec is not None:
            return self._reference_spec
        if self.spec is not None:
            return self.spec
        return ComparisonSpec.union(*[self.case_spec(c, resolve_pins=False)
                                      for c in self.cases()])

    # ── DRAGON ─────────────────────────────────────────────────────────────
    def dragon_case(self, case: StudyCase):
        """``DragonCase`` of a study case (procedures not generated yet)."""
        from ..InterfaceToDD.case_generator import DragonCase
        asm = case.assembly
        dc = DragonCase(case_name=asm.case_name, call_glow=False,
                        draglib_name_to_alias=dict(self.draglib),
                        config_yamls={"MATS": str(asm.materials), "GEOM": str(asm.geometry),
                                      "CALC_SCHEME": case.scheme.copy()},
                        output_path=str(case.directory / PROCEDURES_DIR),
                        tdt_path=str(asm.tdt_path), tdt_base_name=asm.tdt_base_name)
        # Fuel temperature, coolant/moderator temperature and density come
        # from the branch loop of the x2m; this sets the structures.
        dc.set_non_fuel_temperatures(structural_temperature=self.structural_temperature)
        return dc

    def _generate_case(self, case: StudyCase):
        case.directory.mkdir(parents=True, exist_ok=True)
        dc = self.dragon_case(case)
        dc.generate_cle2000_procedures()
        case.scheme.to_yaml(case.directory / SCHEME_FILE,
                            header=f"study: {self.name}\ncase: {case.key}")
        rec = case.record
        keep_status = rec.get("status") == "completed" and rec.get("fingerprint") == case.fingerprint
        rec.update({
            "key": case.key, "assembly": case.assembly.name, "variant": case.variant.name,
            "condition": case.condition.name, "case_name": case.case_name,
            "fingerprint": case.fingerprint, "variant_tag": case.scheme.variant_tag,
            "descriptor": case.scheme.descriptor(), "scheme": SCHEME_FILE,
            "statepoints": [sp.to_dict() for sp in case.statepoints],
            "generated": {"time": _now(), "procedures": PROCEDURES_DIR,
                          "files": sorted(os.listdir(case.directory / PROCEDURES_DIR)),
                          "tdt_files": dict(dc.tdt_files_used)},
            "status": rec["status"] if keep_status else "generated",
        })
        case.write_record(rec)
        return dc

    def generate(self, assemblies=None, variants=None, conditions=None,
                 force: bool = False) -> List[StudyCase]:
        """Write the CLE-2000 procedures of the selected cases (needs the TDT
        files, not DRAGON). Up-to-date cases are skipped unless *force*."""
        done = []
        for case in self.cases(assemblies, variants, conditions):
            if not force and case.status in ("generated", "dry_run", "completed"):
                print(f"[study] {case.key}: up to date ({case.status})")
                continue
            print(f"[study] {case.key}: generating procedures")
            self._generate_case(case)
            done.append(case)
        self.write_manifest()
        return done

    def run_dragon(self, num_threads: Optional[int] = None, dry_run: bool = False,
                   force: bool = False, assemblies=None, variants=None, conditions=None,
                   stop_on_error: bool = False) -> Dict[str, str]:
        """Generate and run the selected cases, one after the other.

        Completed, up-to-date cases are skipped unless *force*. *dry_run*
        stages everything (procedures, TDT links, manifest) without running
        DRAGON, so the whole tree can be checked first. A failing case is
        recorded and the next one starts, unless *stop_on_error*.

        Returns ``{case key: status}``.
        """
        out: Dict[str, str] = {}
        for case in self.cases(assemblies, variants, conditions):
            if not force and case.status == "completed":
                print(f"[study] {case.key}: completed, skipped (force=True to re-run)")
                out[case.key] = "skipped"
                continue
            print(f"[study] {case.key}: {'dry run' if dry_run else 'running DRAGON'} "
                  f"({len(case.statepoints)} statepoints)")
            error = None
            try:
                dc = self._generate_case(case)
                res = dc.run(dragon_executable=self.dragon_executable,
                             draglib_paths=self.draglib_paths or None,
                             results_root=str(case.directory),
                             num_threads=num_threads or self.num_threads,
                             dry_run=dry_run, run_subdirs=[])
                entry = {"timestamp": Path(res.run_directory).name,
                         "run_directory": _rel(res.run_directory, case.directory),
                         "dry_run": dry_run, "success": bool(res.success),
                         "return_code": res.return_code,
                         "wall_time_seconds": float(res.wall_time_seconds),
                         "kinf": res.kinf}
                status = "dry_run" if dry_run else ("completed" if res.success else "failed")
            except Exception as exc:            # recorded, then next case
                error = exc
                entry = {"timestamp": _dt.datetime.now().strftime("%Y-%m-%dT%H-%M-%S"),
                         "dry_run": dry_run, "success": False,
                         "error": f"{type(exc).__name__}: {exc}"}
                status = "failed"
            entry["fingerprint"] = case.fingerprint
            rec = case.record
            rec.setdefault("runs", []).append(entry)
            rec["status"] = status
            if status == "completed":
                rec["current_run"] = entry["run_directory"]
            case.write_record(rec)
            self.write_manifest()
            out[case.key] = status
            print(f"[study] {case.key}: {status}"
                  + (f" ({entry['error']})" if "error" in entry else "")
                  + (f", kinf = {entry['kinf']}" if entry.get("kinf") else ""))
            if error is not None and stop_on_error:
                raise error
        return out

    # ── Serpent2 references ────────────────────────────────────────────────
    def _serpent_settings(self, case: StudyCase, sp: StudyStatepoint):
        from ..InterfaceToDD.Serpent2_exports import S2_Settings
        s = S2_Settings()
        s.set_nuclear_data_evaluation(self.library)
        for k, v in self.serpent_settings.items():
            if not hasattr(s, k):
                raise KeyError(f"Unknown Serpent setting '{k}' (an S2_Settings attribute).")
            setattr(s, k, v)
        vals = ", ".join(f"{k}={v:g}" for k, v in sp.values.items())
        s.title = f"{case.case_name} {sp.label} - starterDD study {self.name} ({vals})"
        return s

    def _reference_cases(self, assemblies=None, conditions=None) -> List[StudyCase]:
        """One case per (assembly, condition): references do not depend on
        the DRAGON variant."""
        first = next(iter(self.variants))
        return [self.case(a, first, c) for a in self._pick(self.assemblies, assemblies)
                for c in self._pick(self.conditions, conditions)]

    def write_serpent_references(self, assemblies=None, conditions=None,
                                 overwrite: bool = False) -> list:
        """Write one Serpent2 input per (assembly, statepoint) in the
        registry, with the detectors of :meth:`reference_spec` and the
        statepoint applied as in DRAGON (``Serpent2Model.from_dragon_case``).
        Existing references, generated or imported, are kept unless
        *overwrite*. Serpent2 itself is not run."""
        from ..InterfaceToDD.Serpent2_exports import Serpent2Model
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        spec = self.reference_spec()
        created = []
        for case in self._reference_cases(assemblies, conditions):
            dc = None
            for sp in case.statepoints:
                folder = case.reference_folder(sp)
                if (folder / REFERENCE_FILE).exists() and not overwrite:
                    print(f"[study] reference {_rel(folder, self.references_root)}: exists, kept")
                    continue
                dc = dc or self.dragon_case(case)
                model = Serpent2Model.from_dragon_case(
                    dc, statepoint=dict(sp.values), spec=spec, library=self.library,
                    settings=self._serpent_settings(case, sp),
                    structural_temperature=self.structural_temperature)
                ref = SerpentReference.create(self.references_root, case.assembly.name,
                                              self.library, sp.label, model, spec=spec,
                                              overwrite=overwrite)
                print(f"[study] reference {_rel(ref.input_path, self.references_root)}: written")
                created.append(ref)
        gaps = self.check_references(assemblies=assemblies, conditions=conditions)
        for k, v in gaps.items():
            warnings.warn(f"Serpent reference {k}: {v}", stacklevel=2)
        self.write_manifest()
        return created

    def import_serpent_references(self, mapping: Optional[Mapping[str, Mapping[str, PathLike]]] = None,
                                  link: bool = True, overwrite: bool = False,
                                  pin_resolver=None) -> list:
        """Register existing Serpent2 runs (decision 8, no re-run):
        ``{assembly: {statepoint label: .serp path}}``, default the
        ``serpent: import:`` block of the study file."""
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        mapping = mapping if mapping is not None else self.serpent_imports
        labels = {a: {sp.label for c in self.cases(assemblies=a) for sp in c.statepoints}
                  for a in self.assemblies}
        out = []
        for asm, sps in mapping.items():
            if asm not in self.assemblies:
                raise KeyError(f"Unknown assembly '{asm}' in the Serpent imports.")
            for label, serp in sps.items():
                if label not in labels[asm]:
                    warnings.warn(f"{asm}: '{label}' is not a statepoint of the study "
                                  f"({sorted(labels[asm])}); imported anyway.", stacklevel=2)
                folder = SerpentReference.folder(self.references_root, asm, self.library, label)
                if (folder / REFERENCE_FILE).exists() and not overwrite:
                    print(f"[study] reference {asm}/{self.library}/{label}: exists, kept")
                    continue
                out.append(SerpentReference.import_existing(
                    serp, self.references_root, asm, self.library, label, link=link,
                    pin_resolver=pin_resolver, overwrite=overwrite))
                print(f"[study] reference {asm}/{self.library}/{label}: imported from {serp}")
        self.write_manifest()
        return out

    # ── existing DRAGON runs ───────────────────────────────────────────────
    def check_dragon_compo(self, case: StudyCase, compo_path: PathLike,
                           cache: bool = True) -> Tuple[List[str], List[str]]:
        """``(problems, notes)`` of using *compo_path* as the COMPO of *case*.

        Problems: a DRAGON output that :meth:`post_process` reads is missing
        (the fuel 1-group output, the coarsest PIN output, the ALL flux
        outputs), a statepoint of the case matches no (or several)
        calculation, or the COMPO lacks a parameter needed to tell the
        statepoints apart. Permuted parameter names (issue 15, COMPOs written
        before 2026-11-15) show up as unmatched statepoints: repair them with
        ``PostTreatment.compo_repair`` first. Notes: other missing outputs,
        parameters that cannot be checked (constant in the case).
        """
        from ..PostTreatment.compo_reader import CompoReader
        problems: List[str] = []
        notes: List[str] = []
        spec = self.case_spec(case)
        used = set(_post_processed_observables(spec))
        with CompoReader(compo_path, cache=cache) as compo:
            outs = []
            for o in spec:
                if not o.dragon:
                    continue
                try:
                    outs.append(_dragon_output(compo, o))
                except KeyError as exc:
                    (problems if o.name in used else notes).append(
                        f"{o.name}: {exc}" + ("" if o.name in used else " (not post-processed)"))
            sps = case.statepoints
            for out in outs:
                have = compo.parameters(out)
                missing = sorted({k for sp in sps for k in sp.params} - set(have))
                varying = [k for k in missing if len({sp.params[k] for sp in sps}) > 1]
                if varying:
                    problems.append(f"{out}: no COMPO parameter(s) {varying}: cannot tell "
                                    f"the {len(sps)} statepoints apart")
                    continue
                if missing:
                    notes.append(f"{out}: no COMPO parameter(s) {missing} (constant in the "
                                 "case, values not checked)")
                calcs = {}
                for sp in sps:
                    try:
                        calcs[sp.label] = compo.calc_index(out, _selection(compo, out, sp))
                    except Exception as exc:
                        problems.append(f"{out} {sp.label}: {exc}"
                                        + (" (permuted parameter names? see compo_repair)"
                                           if "match" in str(exc) else ""))
                if len(set(calcs.values())) < len(calcs):
                    problems.append(f"{out}: statepoints {calcs} share calculations")
        return problems, notes

    def import_dragon_run(self, assembly: str, variant: str, condition: str,
                          source: PathLike, link: bool = True, force: bool = False,
                          check: bool = True) -> StudyCase:
        """Register an existing DRAGON COMPO as the result of a case, without
        running DRAGON (post-processing only).

        *source* is a COMPO file or a run folder (``_CPO_<case_name>``, or its
        only ``_CPO_*``). The COMPO (and its ``.h5`` cache when present) is
        linked (copied if not *link*) into ``<case>/imported_<time>/``, so
        the case looks like a completed run (``origin: imported`` in
        ``case.yaml``). With *check*, :meth:`check_dragon_compo` must find no
        problem. If the run folder has an archived scheme
        (``inputs/CALC_SCHEME_resolved.yaml``), its hash is compared with the
        case scheme and recorded as ``scheme_match``. A completed case is
        kept unless *force*.
        """
        import shutil
        from ..PostTreatment.compo_reader import default_cache_path
        case = self.case(assembly, variant, condition)
        src = Path(source).expanduser().resolve()
        if src.is_dir():
            run_dir = src
            compo = src / f"_CPO_{case.case_name}"
            if not compo.exists():
                hits = sorted(q for q in src.glob("_CPO_*") if q.suffix != ".h5")
                if len(hits) != 1:
                    raise FileNotFoundError(f"{src}: {'no' if not hits else 'several'} "
                                            f"_CPO_* file(s) {[h.name for h in hits]}; "
                                            "give the COMPO file.")
                compo = hits[0]
        else:
            run_dir, compo = src.parent, src
        if not compo.exists():
            raise FileNotFoundError(f"COMPO not found: {compo}")
        if case.status == "completed" and not force:
            print(f"[study] {case.key}: completed, kept (force=True to import {compo})")
            return case

        stamp =_dt.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        folder = case.directory / f"imported_{stamp}"
        folder.mkdir(parents=True, exist_ok=False)
        dest = folder / f"_CPO_{case.case_name}"
        pairs = [(compo, dest)]
        if default_cache_path(compo).exists():
            pairs.append((default_cache_path(compo), default_cache_path(dest)))
        for a, b in pairs:
            if link:
                b.symlink_to(a)
            else:
                shutil.copy2(a, b)
        notes: List[str] = []
        if check:
            problems, notes = self.check_dragon_compo(case, dest)
            if problems:
                shutil.rmtree(folder)
                raise ValueError(f"{compo} cannot be the COMPO of {case.key}:\n  "
                                 + "\n  ".join(problems))
        scheme_match = None
        archived = run_dir / "inputs" / "CALC_SCHEME_resolved.yaml"
        if archived.exists():
            try:
                scheme_match = (DragonCalculationScheme.from_yaml(str(archived)).content_hash(16)
                                == case.scheme.content_hash(16))
            except Exception as exc:          # recorded, not fatal
                warnings.warn(f"{archived}: cannot compare schemes ({exc}).", stacklevel=2)
            if scheme_match is False:
                warnings.warn(f"{case.key}: the imported run's scheme differs from the "
                              "case scheme (scheme_match: false).", stacklevel=2)
        info = {"source": str(compo), "linked": link, "imported": _now(),
                "source_size": compo.stat().st_size, "run_directory": str(run_dir),
                "scheme_match": scheme_match, "notes": notes}
        _atomic_write_yaml(folder / "import.yaml", info)
        case.scheme.to_yaml(case.directory / SCHEME_FILE,
                            header=f"study: {self.name}\ncase: {case.key}")
        rec = case.record
        rec.update({
            "key": case.key, "assembly": case.assembly.name, "variant": case.variant.name,
            "condition": case.condition.name, "case_name": case.case_name,
            "fingerprint": case.fingerprint, "variant_tag": case.scheme.variant_tag,
            "descriptor": case.scheme.descriptor(), "scheme": SCHEME_FILE,
            "statepoints": [sp.to_dict() for sp in case.statepoints],
            "status": "completed", "origin": "imported", "current_run": folder.name,
        })
        rec.setdefault("runs", []).append({
            "timestamp": stamp, "run_directory": folder.name, "dry_run": False,
            "success": True, "imported_from": str(compo), "scheme_match": scheme_match,
            "fingerprint": case.fingerprint})
        case.write_record(rec)
        self.write_manifest()
        print(f"[study] {case.key}: imported {compo}")
        return case

    def import_dragon_runs(self, mapping: Optional[Mapping[str, PathLike]] = None,
                           **kwargs) -> List[StudyCase]:
        """:meth:`import_dragon_run` for ``{"<assembly>/<variant>/<condition>":
        COMPO or run folder}``, default the study file's ``dragon: import:``
        block."""
        mapping = mapping if mapping is not None else self.dragon_imports
        out = []
        for key, src in mapping.items():
            parts = str(key).split("/")
            if len(parts) != 3:
                raise KeyError(f"DRAGON import key '{key}' must be <assembly>/<variant>/<condition>.")
            out.append(self.import_dragon_run(*parts, source=src, **kwargs))
        return out

    def check_references(self, assemblies=None, variants=None, conditions=None) -> Dict[str, Any]:
        """``{"<assembly>/<sp>/<variant>": problem}`` for every statepoint
        whose reference is missing or does not cover the case spec. Empty =
        every case can be post-processed once Serpent has run."""
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        problems: Dict[str, Any] = {}
        for case in self.cases(assemblies, variants, conditions):
            spec = self.case_spec(case)
            n_pins = self.lattice(case.assembly).n_pins
            for sp in case.statepoints:
                key = f"{case.assembly.name}/{sp.label}/{case.variant.name}"
                folder = case.reference_folder(sp)
                if not (folder / REFERENCE_FILE).exists():
                    problems[key] = "missing"
                    continue
                gaps = SerpentReference.open(folder).coverage(spec, n_pins=n_pins)
                if gaps:
                    problems[key] = {k: "; ".join(v) for k, v in gaps.items()}
        return problems

    # ── manifest and status ────────────────────────────────────────────────
    def to_dict(self) -> dict:
        """The study definition (absolute paths)."""
        return _plain({
            "study": {"name": self.name, "root": str(self.root), "library": self.library,
                      "source": None if self.source is None else str(self.source)},
            "dragon": {"draglib": self.draglib, "draglib_paths": self.draglib_paths,
                       "executable": self.dragon_executable, "num_threads": self.num_threads,
                       "structural_temperature": self.structural_temperature},
            "serpent": {"references_root": str(self.references_root),
                        "settings": self.serpent_settings},
            "spec": self.spec_config if self.spec_config is not None else
            (None if self.spec is None else self.spec.names()),
            "reference_spec": self.reference_spec_config,
            "scheme": None if self.scheme is None else str(self.scheme),
            "assemblies": {k: a.to_dict() for k, a in self.assemblies.items()},
            "variants": {k: v.to_dict() for k, v in self.variants.items()},
            "conditions": {k: c.to_dict() for k, c in self.conditions.items()},
            "post_processing": self.post_processing,
        })

    def manifest_dict(self) -> dict:
        """Everything post-processing needs: per case, its status, run folder,
        COMPO and, per statepoint, the COMPO parameters and the reference."""
        from .. import __version__
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        cases = []
        for case in self.cases():
            sps = []
            for sp in case.statepoints:
                folder = case.reference_folder(sp)
                exists = (folder / REFERENCE_FILE).exists()
                d = sp.to_dict()
                d.update(reference=_rel(folder, self.directory), reference_exists=exists,
                         reference_has_results=exists and SerpentReference.open(folder).has_results)
                sps.append(d)
            cases.append({
                "key": case.key, "assembly": case.assembly.name, "variant": case.variant.name,
                "condition": case.condition.name, "case_name": case.case_name,
                "directory": _rel(case.directory, self.directory), "status": case.status,
                "variant_tag": case.scheme.variant_tag, "fingerprint": case.fingerprint,
                "run_directory": _rel(case.run_directory, self.directory),
                "compo": _rel(case.compo_path, self.directory),
                "statepoints": sps,
            })
        try:
            ref_spec = self.reference_spec().names()
        except Exception as exc:              # reported, not fatal for the manifest
            ref_spec = f"<error: {exc}>"
        return _plain({
            "study": self.name, "updated": _now(), "starterDD_version": __version__,
            "library": self.library, "draglib": self.draglib,
            "references_root": _rel(self.references_root, self.directory),
            "spec": None if self.spec is None else self.spec.names(),
            "reference_spec": ref_spec, "cases": cases,
        })

    def write_manifest(self) -> Path:
        self.directory.mkdir(parents=True, exist_ok=True)
        _atomic_write_yaml(self.directory / CONFIG_FILE, self.to_dict())
        _atomic_write_yaml(self.manifest_path, self.manifest_dict())
        return self.manifest_path

    def status(self) -> List[Dict[str, Any]]:
        """One row per case: status, statepoints, references written / with
        results, last run."""
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        rows = []
        for case in self.cases():
            refs = [case.reference_folder(sp) for sp in case.statepoints]
            have = [r for r in refs if (r / REFERENCE_FILE).exists()]
            runs = case.record.get("runs") or [{}]
            rows.append({"case": case.key, "status": case.status,
                         "statepoints": ",".join(sp.label for sp in case.statepoints),
                         "references": f"{len(have)}/{len(refs)}",
                         "serpent_results": f"{sum(SerpentReference.open(r).has_results for r in have)}"
                                            f"/{len(refs)}",
                         "last_run": runs[-1].get("timestamp", "")})
        return rows

    def print_status(self) -> None:
        rows = self.status()
        cols = list(rows[0])
        width = {c: max(len(c), *(len(str(r[c])) for r in rows)) for c in cols}
        print("  ".join(c.ljust(width[c]) for c in cols))
        for r in rows:
            print("  ".join(str(r[c]).ljust(width[c]) for c in cols))

    # ── post-processing ────────────────────────────────────────────────────
    def post_process(self, isotope_sets: Optional[Sequence[Sequence[str]]] = None,
                     normalisation: Optional[str] = None, assemblies=None, variants=None,
                     conditions=None, out_dir: Optional[PathLike] = None) -> List[dict]:
        """Code-to-code comparison of every completed case against its
        references that have results. Writes ``post/study.h5`` (needs h5py),
        ``summary.csv`` and ``summary.json``; returns the summary rows.

        Per statepoint: keff (the FUEL 1-group output vs ``IMP_KEFF``); the
        normalisation; pin rates on the coarsest PIN observable (macroscopic
        rates if the spec scores them, and per isotope set if it scores
        isotopes); spectra of the ALL flux observables. Statepoints that
        cannot be compared get a row with the reason in ``status``.

        *normalisation* (default ``fission_source``): a whole-assembly method
        (``fission_source``, ``fission``, ``absorption``), computed once per
        statepoint from the macroscopic rates (``global_totals`` of the spec)
        and applied to every rate comparison of that statepoint, or a subset
        method (``power_subset``, ``absorption_subset``, ``fission_subset``,
        ``pin_average``). A reference without the needed ``void`` responses
        (legacy runs) falls back to ``power_subset``, reported in the
        ``normalisation`` row and in ``study.h5``.
        """
        from ..PostTreatment.comparison import NORMALISATIONS
        from ..PostTreatment.compo_reader import CompoReader
        from ..PostTreatment.references import REFERENCE_FILE, SerpentReference
        from ..PostTreatment.reports import DEFAULT_ISOTOPE_SETS

        pp = self.post_processing
        isotope_sets = [list(s) for s in (isotope_sets or pp.get("isotope_sets")
                                          or DEFAULT_ISOTOPE_SETS)]
        normalisation = normalisation or pp.get("normalisation", "fission_source")
        if normalisation == "power":
            normalisation = "power_subset"
        if normalisation not in NORMALISATIONS:
            raise ValueError(f"Unknown normalisation '{normalisation}' ({', '.join(NORMALISATIONS)}).")
        out = Path(out_dir) if out_dir is not None else self.directory / "post"
        out.mkdir(parents=True, exist_ok=True)
        try:
            import h5py
        except ImportError:
            h5py = None
            warnings.warn("h5py not installed: study.h5 is not written.", stacklevel=2)
        h5 = h5py.File(out / "study.h5", "w") if h5py is not None else None
        rows: List[dict] = []
        try:
            if h5 is not None:
                h5.attrs["study"] = self.name
                h5.attrs["library"] = self.library
                h5.attrs["normalisation"] = normalisation
                h5.attrs["created"] = _now()
            for case in self.cases(assemblies, variants, conditions):
                base = {"assembly": case.assembly.name, "variant": case.variant.name,
                        "condition": case.condition.name}
                compo_path = case.compo_path
                if compo_path is None:
                    rows.append({**base, "item": "case", "status": f"no COMPO ({case.status})"})
                    continue
                compo = CompoReader(compo_path)
                spec = self.case_spec(case)
                lattice = self.lattice(case.assembly)
                leakage = sorted({s.leakage_model for s in case.scheme.get_trackable_steps()
                                  if getattr(s, "leakage_model", None)})
                if h5 is not None:
                    _write_lattice(h5, case.assembly.name, lattice)
                    cg = h5.require_group(case.key)
                    cg.attrs["spec"] = json.dumps(_spec_description(spec))
                    cg.attrs["leakage_model"] = ",".join(leakage)
                    cg.attrs["status"] = case.status
                    cg.attrs["origin"] = case.record.get("origin", "run")
                try:
                    for sp in case.statepoints:
                        row0 = {**base, "statepoint": sp.label}
                        folder = case.reference_folder(sp)
                        if not (folder / REFERENCE_FILE).exists():
                            rows.append({**row0, "item": "statepoint", "status": "no reference"})
                            continue
                        ref = SerpentReference.open(folder)
                        if not ref.has_results:
                            rows.append({**row0, "item": "statepoint",
                                         "status": "reference not run yet"})
                            continue
                        grp = None
                        if h5 is not None:
                            grp = h5.require_group(f"{case.key}/{sp.label}")
                            grp.attrs["compo"] = str(compo_path)
                            grp.attrs["reference"] = str(folder)
                            for k, v in sp.params.items():
                                grp.attrs[f"param_{k}"] = v
                        rows += _compare_statepoint(compo, ref, spec, lattice, sp, self.library,
                                                    isotope_sets, normalisation, row0, grp,
                                                    leakage=leakage)
                finally:
                    compo.close()
        finally:
            if h5 is not None:
                h5.close()
        _write_summary(rows, out)
        print(f"[study] post-processing: {len(rows)} rows -> {out}")
        return rows

    def figures(self, post_dir: Optional[PathLike] = None, **kwargs) -> List[Path]:
        """Figures from ``<post_dir>/study.h5`` (default ``<study>/post``),
        see :func:`study_figures.study_figures` (``items=``, ``quantities=``)."""
        from .study_figures import study_figures
        post = Path(post_dir) if post_dir is not None else self.directory / "post"
        return study_figures(post / "study.h5", **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing helpers
# ─────────────────────────────────────────────────────────────────────────────

SUMMARY_COLUMNS = ["assembly", "variant", "condition", "statepoint", "item", "observable",
                   "isotopes", "quantity", "group", "dragon", "serpent", "sigma", "delta",
                   "mean", "rms", "min", "max", "max_abs", "pin_max_abs", "mesh",
                   "normalisation", "factor", "status"]


def _spec_description(spec) -> List[dict]:
    """JSON-able description of the observables (``study.h5`` attribute):
    which pins the Serpent detectors cover (``serpent_pins``, ``None`` =
    all), isotopes, macroscopic rates."""
    return [{"name": o.name, "output": o.dragon_name, "spatial": o.spatial,
             "n_groups": o.n_groups, "isotopes": list(o.isotopes),
             "serpent": list(o.serpent), "serpent_isotopes": list(o.micro_isotopes),
             "macro": list(o.macro), "dragon": o.dragon,
             "serpent_pins": None if o.pins is None else [int(p) for p in o.pins]}
            for o in spec]


def _write_lattice(h5, assembly: str, lattice) -> None:
    """``/<assembly>/lattice``: pin positions, so that maps can be drawn
    from ``study.h5`` alone."""
    import numpy as np
    g = h5.require_group(assembly)
    if "lattice" in g:
        return
    lg = g.create_group("lattice")
    lg["ordered_pins"] = np.array(lattice.ordered_pins, dtype=int)
    lg["positions"] = np.array([(p, i, j) for p, pos in sorted(lattice.positions.items())
                                for (i, j) in pos], dtype=int)
    lg.attrs["shape"] = list(lattice.shape)
    lg.attrs["symmetry"] = lattice.symmetry or ""
    lg.attrs["material"] = json.dumps({str(k): v for k, v in lattice.material.items()})


def _dragon_output(compo, obs) -> str:
    from ..PostTreatment.spec import find_output
    return obs.dragon_name if obs.dragon_name in compo.outputs() else find_output(compo, obs.name)


def _selection(compo, output: str, sp: StudyStatepoint) -> Dict[str, float]:
    """COMPO selection of *sp*: the parameters this output has."""
    have = compo.parameters(output)
    return {k: v for k, v in sp.params.items() if k in have}


def _with_mesh_fallback(fn, *args, **kwargs):
    """Compare on the DRAGON mesh; on the common mesh if the reference grid
    does not refine it (legacy SHEM295 references)."""
    from ..PostTreatment.physics_data import MeshMismatchError
    try:
        return fn(*args, mesh="dragon", **kwargs), "dragon"
    except MeshMismatchError:
        return fn(*args, mesh="common", **kwargs), "common"


def _isotope_source(ref, obs, isotopes, n_pins: int) -> str:
    """Serpent observable of *ref* giving *obs* for these *isotopes* only, so
    that a reference scoring fewer isotopes than the DRAGON edit (e.g. the
    legacy GE14 ones, without U234/U236) still serves the sets it covers."""
    from ..PostTreatment.spec import ComparisonSpec, Observable
    sub = ComparisonSpec([Observable(obs.spatial, obs.bounds_eV, list(isotopes), serpent=("rr",),
                                     dragon=False, pins=obs.pins, tag=obs.tag,
                                     reactions=obs.reactions)])
    return ref.source(sub, obs.name, "rr", n_pins=n_pins)


def _group_rows(row0: dict, stats: List[dict], err, pins, mesh: str) -> List[dict]:
    import numpy as np
    out = []
    for g, st in enumerate(stats):
        e = err[:, g]
        ok = np.isfinite(e)
        i = int(np.nanargmax(np.abs(np.where(ok, e, np.nan)))) if ok.any() else None
        out.append({**row0, "group": g + 1, "mean": st.get("mean"), "rms": st.get("rms"),
                    "min": st.get("min"), "max": st.get("max"),
                    "max_abs": None if i is None else float(abs(e[i])),
                    "pin_max_abs": None if i is None or pins is None else pins[i],
                    "mesh": mesh, "status": "ok"})
    return out


def _keff_observables(spec) -> list:
    return [o for o in spec if o.spatial == "FUEL" and o.n_groups == 1 and o.dragon][:1]


def _pin_observables(spec) -> list:
    """The coarsest PIN observable with Serpent reaction rates (compared)."""
    return sorted((o for o in spec if o.spatial == "PIN" and "rr" in o.serpent and o.dragon),
                  key=lambda o: o.n_groups)[:1]


def _spectrum_observables(spec) -> list:
    return [o for o in spec if o.spatial == "ALL" and "flux" in o.serpent and o.dragon]


def _post_processed_observables(spec) -> List[str]:
    """Names of the observables whose DRAGON outputs post_process reads."""
    return [o.name for o in _keff_observables(spec) + _pin_observables(spec)
            + _spectrum_observables(spec)]


def _macro_source(ref, obs, quantities, n_pins: int) -> str:
    """Serpent observable of *ref* giving the macroscopic *quantities* of
    *obs* (whatever per-isotope responses it has)."""
    from ..PostTreatment.spec import ComparisonSpec, Observable
    sub = ComparisonSpec([Observable(obs.spatial, obs.bounds_eV, list(obs.isotopes),
                                     serpent=("rr",), dragon=False, pins=obs.pins, tag=obs.tag,
                                     macro=tuple(quantities), serpent_isotopes=())])
    return ref.source(sub, obs.name, "rr", n_pins=n_pins)


def _statepoint_normalisation(method: str, compo, ref, serp, spec, lattice,
                              sp: StudyStatepoint, leakage: Sequence[str]):
    """``(Normalisation or None, subset method, note)``. A whole-assembly
    *method* that the spec or the reference cannot provide falls back to
    ``power_subset`` (with the reason in *note*)."""
    from ..PostTreatment.comparison import (GLOBAL_NORMALISATIONS, compute_macro_rates,
                                            global_normalisation)
    if method not in GLOBAL_NORMALISATIONS:
        return None, method, ""
    spatial, q = GLOBAL_NORMALISATIONS[method]
    cands = sorted((o for o in spec if o.spatial == spatial and q in o.macro
                    and "rr" in o.serpent and o.dragon), key=lambda o: o.n_groups)
    reasons = []
    if not cands:
        reasons.append(f"the spec scores no {spatial} '{q}' (global_totals)")
    for o in cands:
        try:
            out = _dragon_output(compo, o)
            src = _macro_source(ref, o, [q], lattice.n_pins)
            cmp, _ = _with_mesh_fallback(compute_macro_rates, compo, out, serp, src,
                                         quantities=[q], lattice=lattice,
                                         sp=_selection(compo, out, sp))
            norm = global_normalisation(method, cmp, lattice.dragon_geometry_fraction)
            if method == "absorption" and leakage:
                norm.note = (f"DRAGON flux with leakage model {','.join(leakage)}: the "
                             "absorption balance includes leakage if the buckling is non-zero")
                warnings.warn(f"{sp.label}: {norm.note}.", stacklevel=3)
            return norm, method, norm.note
        except Exception as exc:          # next candidate, then fallback
            reasons.append(f"{o.name}: {type(exc).__name__}: {exc}")
    return None, "power_subset", f"fallback to power_subset ({'; '.join(reasons)})"


def _compare_statepoint(compo, ref, spec, lattice, sp: StudyStatepoint, library: str,
                        isotope_sets, normalisation: str, row0: dict, grp,
                        leakage: Sequence[str] = ()) -> List[dict]:
    import numpy as np
    from ..PostTreatment.comparison import (compare_keff, compare_spectra, compute_macro_rates,
                                            compute_pin_rates)
    rows: List[dict] = []
    serp = ref.run()
    n_pins = lattice.n_pins

    def failed(item, obs, exc):
        rows.append({**row0, "item": item, "observable": obs,
                     "status": f"error: {type(exc).__name__}: {exc}"})

    # keff
    for o in _keff_observables(spec):
        try:
            out = _dragon_output(compo, o)
            k = compare_keff(compo, out, serp, _selection(compo, out, sp))
            rows.append({**row0, "item": "keff", "observable": out, "dragon": k.dragon,
                         "serpent": k.serpent, "sigma": k.sigma_pcm, "delta": k.delta_pcm,
                         "status": "ok"})
            if grp is not None:
                g = grp.require_group("keff")
                for a in ("dragon", "serpent", "sigma_pcm", "delta_pcm"):
                    g.attrs[a] = getattr(k, a)
                g.attrs["estimator"] = k.estimator
                g.attrs["output"] = out
        except Exception as exc:
            failed("keff", o.name, exc)

    # normalisation: one factor per statepoint for every rate comparison
    norm, subset, note = _statepoint_normalisation(normalisation, compo, ref, serp, spec,
                                                   lattice, sp, leakage)
    factor = None if norm is None else norm.factor
    used = normalisation if norm is not None else subset
    nrow = {**row0, "item": "normalisation", "quantity": normalisation, "normalisation": used,
            "status": "ok" if not note else note}
    if norm is not None:
        nrow.update(observable=norm.output, dragon=norm.dragon_total, serpent=norm.serpent_total,
                    sigma=100.0 * norm.serpent_rel_err, factor=norm.factor)
    rows.append(nrow)
    if grp is not None:
        g = grp.require_group("normalisation")
        g.attrs["requested"] = normalisation
        g.attrs["method"] = used
        g.attrs["note"] = note
        if norm is not None:
            for a in ("factor", "dragon_total", "serpent_total", "serpent_rel_err"):
                g.attrs[a] = getattr(norm, a)
            g.attrs["output"] = norm.output
            g.attrs["observable"] = norm.observable

    # macroscopic rates without a whole-assembly factor: shape only
    # (pin_average) or the total fission of the compared pins (default)
    m_method = "pin_average" if used == "pin_average" else None

    # pin rates on the coarsest PIN observable with Serpent reaction rates
    for o in _pin_observables(spec):
        try:
            out = _dragon_output(compo, o)
            sel = _selection(compo, out, sp)
        except Exception as exc:
            failed("pin_rates", o.name, exc)
            continue
        src = None
        if o.macro:
            try:
                src = _macro_source(ref, o, o.macro, n_pins)
            except KeyError as exc:          # e.g. legacy references: no 'void' responses
                rows.append({**row0, "item": "macro_rates", "observable": out,
                             "status": f"not scored by the reference: {exc}"})
        if src is not None:
            try:
                cmp, mesh = _with_mesh_fallback(compute_macro_rates, compo, out, serp, src,
                                                quantities=list(o.macro), lattice=lattice, sp=sel)
                for q in cmp.quantities:
                    err = cmp.errors(q, m_method, factor)
                    k = factor if factor is not None else cmp.normalisation_factor(m_method, q)
                    r0 = {**row0, "item": "macro_rates", "observable": out, "isotopes": "all",
                          "quantity": q, "normalisation": used, "factor": k}
                    rows += _group_rows(r0, cmp.stats(q, m_method, factor), err,
                                        cmp.regions, mesh)
                    if grp is not None:
                        g = grp.require_group(f"macro_rates/{o.name}/{q}")
                        g["dragon"], g["serpent"] = cmp.dragon[q], cmp.serpent[q]
                        g["serpent_err"], g["errors"] = cmp.serpent_err[q], err
                        g.attrs["pins"] = [int(p) for p in cmp.regions]
                        g.attrs["bounds_eV"] = cmp.bounds_eV
                        g.attrs["normalisation"] = used
                        g.attrs["factor"] = k
            except Exception as exc:
                failed("macro_rates", o.name, exc)
        if not o.micro_isotopes:
            continue
        for isos in isotope_sets:
            name = "+".join(isos)
            try:
                src = _isotope_source(ref, o, isos, n_pins)
                cmp, mesh = _with_mesh_fallback(compute_pin_rates, compo, out, serp, src,
                                                list(isos), lattice, library, sp=sel)
                quantities = ["absorption"] + (["fission", "fission_over_absorption"]
                                               if cmp.has_fission else [])
                for q in quantities:
                    err = cmp.errors(q, subset, factor)
                    k = (None if q == "fission_over_absorption" else
                         factor if factor is not None else cmp.normalisation_factor(subset, q))
                    r0 = {**row0, "item": "pin_rates", "observable": out, "isotopes": name,
                          "quantity": q, "normalisation": used, "factor": k}
                    rows += _group_rows(r0, cmp.stats(q, subset, factor), err, cmp.pins, mesh)
                    if grp is not None:
                        g = grp.require_group(f"pin_rates/{o.name}/{name}/{q}")
                        g["dragon"], g["serpent"] = cmp.dragon.get(q), cmp.serpent.get(q)
                        g["errors"] = err
                        g.attrs["pins"] = list(cmp.pins)
                        g.attrs["bounds_eV"] = cmp.bounds_eV
                        g.attrs["normalisation"] = used
                        if k is not None:
                            g.attrs["factor"] = k
            except Exception as exc:
                failed("pin_rates", f"{o.name} {name}", exc)

    # spectra of the ALL flux observables
    for o in _spectrum_observables(spec):
        try:
            out = _dragon_output(compo, o)
            src = ref.source(spec, o.name, "flux")
            cmp, mesh = _with_mesh_fallback(compare_spectra, compo, out, serp, src,
                                            sp=_selection(compo, out, sp))
            e = cmp.errors()[0]
            ok = np.isfinite(e)
            rows.append({**row0, "item": "spectrum", "observable": out, "group": "all",
                         "mean": float(np.mean(e[ok])), "rms": float(np.sqrt(np.mean(e[ok] ** 2))),
                         "min": float(np.min(e[ok])), "max": float(np.max(e[ok])),
                         "max_abs": float(np.max(np.abs(e[ok]))), "mesh": mesh, "status": "ok"})
            if grp is not None:
                g = grp.require_group(f"spectra/{o.name}")
                g["dragon"], g["serpent"] = cmp.dragon[0], cmp.serpent[0]
                g["serpent_rel_err"], g["errors"] = cmp.serpent_rel_err[0], e
                g.attrs["bounds_eV"] = cmp.bounds_eV
        except Exception as exc:
            failed("spectrum", o.name, exc)
    return rows


def _write_summary(rows: List[dict], out: Path) -> None:
    import csv
    with open(out / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLUMNS, delimiter=";", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in SUMMARY_COLUMNS})
    (out / "summary.json").write_text(json.dumps(_plain(rows), indent=1, default=str))


# ─────────────────────────────────────────────────────────────────────────────
# Command line
# ─────────────────────────────────────────────────────────────────────────────


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(
        prog="python -m starterDD.Workflow",
        description="Generate / run the DRAGON cases of a code-to-code study, write its "
                    "Serpent2 references, post-process.")
    ap.add_argument("study_yaml")
    ap.add_argument("command", choices=["status", "generate", "run", "serpent",
                                        "import-serpent", "import-dragon", "check", "post",
                                        "figures", "manifest"])
    ap.add_argument("--assembly", action="append", help="repeat to select several")
    ap.add_argument("--variant", action="append")
    ap.add_argument("--condition", action="append")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="re-generate / re-run / overwrite")
    ap.add_argument("--set", action="append", default=[], metavar="VAR=VALUE",
                    help="value of a ${VAR} of the study file")
    a = ap.parse_args(argv)
    variables = dict(s.split("=", 1) for s in a.set)
    study = CodeToCodeStudy.from_yaml(a.study_yaml, variables=variables)
    sel = dict(assemblies=a.assembly, variants=a.variant, conditions=a.condition)
    if a.command == "status":
        study.print_status()
    elif a.command == "generate":
        study.generate(force=a.force, **sel)
    elif a.command == "run":
        res = study.run_dragon(num_threads=a.threads, dry_run=a.dry_run, force=a.force, **sel)
        return 0 if all(v in ("completed", "dry_run", "skipped") for v in res.values()) else 1
    elif a.command == "serpent":
        study.write_serpent_references(assemblies=a.assembly, conditions=a.condition,
                                       overwrite=a.force)
    elif a.command == "import-serpent":
        study.import_serpent_references(overwrite=a.force)
    elif a.command == "import-dragon":
        study.import_dragon_runs(force=a.force)
    elif a.command == "figures":
        for f in study.figures():
            print(f)
    elif a.command == "check":
        problems = study.check_references(**sel)
        for k, v in problems.items():
            print(f"{k}: {v}")
        print("all references cover their cases" if not problems else f"{len(problems)} problem(s)")
    elif a.command == "post":
        study.post_process(**sel)
    elif a.command == "manifest":
        print(study.write_manifest())
    return 0


__all__ = ["CodeToCodeStudy", "StudyAssembly", "StudyVariant", "StudyCondition",
           "StudyStatepoint", "StudyCase", "case_statepoints", "spec_from_config",
           "DEFAULT_DRAGLIBS", "MANIFEST_FILE", "CASE_FILE", "SUMMARY_COLUMNS", "main"]
