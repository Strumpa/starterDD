# -*- coding: utf-8 -*-
"""Figures of a code-to-code study, from ``post/study.h5`` alone (no COMPO,
no Serpent output needed).

:func:`study_figures` writes, under ``<post>/figures/``:

- ``keff_<assembly>.png``: DRAGON - Serpent keff (pcm) per statepoint, one
  marker per variant, Serpent 1σ as error bars;
- ``<assembly>/<variant>/<condition>/<sp>/maps_<item>_<obs>_<set>_<q>.png``:
  per-pin error maps (%), one panel per group, with the normalisation and
  its factor in the title;
- ``.../spectrum_<obs>.png``: assembly spectra per unit lethargy and their
  relative difference.

Maps are drawn on the lattice stored in ``/<assembly>/lattice`` by
:meth:`CodeToCodeStudy.post_process`. Every function builds detached
:class:`matplotlib.figure.Figure` objects (safe in batch / cluster runs).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

PathLike = Union[str, Path]

#: Categorical slots 1-3 (reference palette, light mode): variants of the
#: keff figure, and DRAGON / Serpent in the spectra. Markers double the
#: encoding so identity never relies on color alone.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
MARKERS = ("o", "s", "D", "^", "v", "P", "X", "*")
#: Diverging map (two hues, neutral midpoint) for signed errors.
ERROR_CMAP = "RdBu_r"


def lattice_from_h5(h5, assembly: str):
    """:class:`PinLattice` stored by the post-processing."""
    from ..PostTreatment.lattice_maps import PinLattice
    g = h5[f"{assembly}/lattice"]
    positions: Dict[int, List[Tuple[int, int]]] = {}
    for p, i, j in np.asarray(g["positions"]):
        positions.setdefault(int(p), []).append((int(i), int(j)))
    mat = json.loads(g.attrs.get("material", "{}"))
    return PinLattice(shape=tuple(int(x) for x in g.attrs["shape"]), positions=positions,
                      ordered_pins=[int(p) for p in np.asarray(g["ordered_pins"])],
                      symmetry=g.attrs.get("symmetry") or None,
                      material={int(k): v for k, v in mat.items()})


def statepoint_groups(h5) -> Iterator[Tuple[str, str, str, str, object]]:
    """``(assembly, variant, condition, statepoint, group)`` of every
    compared statepoint of *h5*."""
    for asm in h5:
        for var in h5[asm]:
            if var == "lattice":
                continue
            for cond in h5[f"{asm}/{var}"]:
                for sp in h5[f"{asm}/{var}/{cond}"]:
                    yield asm, var, cond, sp, h5[f"{asm}/{var}/{cond}/{sp}"]


def _norm_label(g) -> str:
    method = g.attrs.get("normalisation", "")
    k = g.attrs.get("factor", None)
    return f"{method} normalisation" + (f", factor {k:.5g}" if k is not None else "")


def plot_keff(h5, assembly: str, path: Optional[PathLike] = None):
    """DRAGON - Serpent keff (pcm) per statepoint of *assembly*, one series
    per variant (statepoints on the x axis, in study order)."""
    from ..PostTreatment.plotting import _figure, _save
    data: Dict[str, Dict[str, Tuple[float, float]]] = {}
    labels: List[str] = []
    for asm, var, cond, sp, g in statepoint_groups(h5):
        if asm != assembly or "keff" not in g:
            continue
        k = g["keff"].attrs
        data.setdefault(var, {})[sp] = (float(k["delta_pcm"]), float(k["sigma_pcm"]))
        if sp not in labels:
            labels.append(sp)
    fig = _figure(figsize=(max(4.0, 1.1 * len(labels) + 2.5), 3.6))
    ax = fig.subplots()
    x = np.arange(len(labels))
    n = len(data)
    width = 0.5 / max(n, 1)
    for i, (var, vals) in enumerate(data.items()):
        xs = [x[labels.index(sp)] + (i - (n - 1) / 2) * width for sp in vals]
        ys = [v[0] for v in vals.values()]
        es = [v[1] for v in vals.values()]
        ax.errorbar(xs, ys, yerr=es, fmt=MARKERS[i % len(MARKERS)], color=SERIES[i % len(SERIES)],
                    markersize=8, elinewidth=1.5, capsize=3, linestyle="none", label=var)
    ax.axhline(0.0, color="0.45", linewidth=1)
    ax.set_xticks(x, labels)
    ax.set_ylabel("keff DRAGON5 - Serpent2 (pcm)")
    ax.set_title(f"{assembly}: keff difference per statepoint")
    ax.grid(True, axis="y", color="0.9")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if n > 1:
        ax.legend(frameon=False, title="variant")
    return _save(fig, path)


def plot_error_maps(h5, sp_path: str, rel: str, path: Optional[PathLike] = None):
    """Per-pin error maps (%) of ``<sp_path>/<rel>`` (a ``pin_rates/...`` or
    ``macro_rates/...`` group), one panel per group."""
    from ..PostTreatment.plotting import _figure, _save, draw_pin_map, group_names
    g = h5[f"{sp_path}/{rel}"]
    lattice = lattice_from_h5(h5, sp_path.split("/")[0])
    err = np.asarray(g["errors"])
    pins = [int(p) for p in g.attrs["pins"]]
    ng = err.shape[1]
    names = group_names(ng)
    fig = _figure(figsize=(5.0 * ng, 4.8))
    axs = np.atleast_1d(fig.subplots(1, ng, squeeze=False)[0])
    for k in range(ng):
        draw_pin_map(axs[k], lattice, err[:, k], pins, cmap=ERROR_CMAP, symmetric=True,
                     title=f"group {names[k]}: DRAGON5 vs Serpent2 (%)")
    parts = rel.split("/")
    fig.suptitle(f"{sp_path.replace('/', ' / ')}\n{' '.join(parts[1:])} ({parts[0]}), "
                 f"{_norm_label(g)}", fontsize=12)
    return _save(fig, path)


def plot_study_spectrum(h5, sp_path: str, obs: str, path: Optional[PathLike] = None):
    """Assembly spectra per unit lethargy and relative difference."""
    from ..PostTreatment.comparison import SpectrumComparison
    from ..PostTreatment.plotting import plot_spectrum
    g = h5[f"{sp_path}/spectra/{obs}"]
    spec = SpectrumComparison(np.asarray(g.attrs["bounds_eV"]), [obs],
                              np.asarray(g["dragon"])[None, :], np.asarray(g["serpent"])[None, :],
                              np.asarray(g["serpent_rel_err"])[None, :])
    return plot_spectrum(spec, path, per_lethargy=True,
                         title=f"{sp_path.replace('/', ' / ')}: {obs}")


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def study_figures(h5_path: PathLike, out_dir: Optional[PathLike] = None,
                  items: Sequence[str] = ("keff", "pin_rates", "macro_rates", "spectra"),
                  quantities: Optional[Sequence[str]] = None) -> List[Path]:
    """Write the figures of a ``study.h5`` (see the module docstring).

    Args:
        out_dir: Default ``<folder of study.h5>/figures``.
        items: Which figures: ``keff``, ``pin_rates``, ``macro_rates``,
            ``spectra``.
        quantities: Restrict the maps to these quantities (e.g.
            ``["absorption"]``).

    Returns the written paths.
    """
    import h5py
    h5_path = Path(h5_path)
    out = Path(out_dir) if out_dir is not None else h5_path.parent / "figures"
    written: List[Path] = []

    def done(p):
        written.append(p)

    with h5py.File(h5_path, "r") as h5:
        if "keff" in items:
            for asm in h5:
                if any("keff" in g for a, *_, g in statepoint_groups(h5) if a == asm):
                    p = out / f"keff_{_safe(asm)}.png"
                    plot_keff(h5, asm, p)
                    done(p)
        for asm, var, cond, sp, g in statepoint_groups(h5):
            sp_path = f"{asm}/{var}/{cond}/{sp}"
            folder = out / _safe(asm) / _safe(var) / _safe(cond) / _safe(sp)
            for item in ("pin_rates", "macro_rates"):
                if item not in items or item not in g:
                    continue
                leaves: List[str] = []
                g[item].visititems(lambda name, obj: leaves.append(name)
                                   if isinstance(obj, h5py.Group) and "errors" in obj else None)
                for leaf in leaves:
                    if quantities and leaf.split("/")[-1] not in quantities:
                        continue
                    p = folder / f"maps_{item}_{_safe(leaf.replace('/', '_'))}.png"
                    plot_error_maps(h5, sp_path, f"{item}/{leaf}", p)
                    done(p)
            if "spectra" in items and "spectra" in g:
                for obs in g["spectra"]:
                    p = folder / f"spectrum_{_safe(obs)}.png"
                    plot_study_spectrum(h5, sp_path, obs, p)
                    done(p)
    return written


__all__ = ["study_figures", "plot_keff", "plot_error_maps", "plot_study_spectrum",
           "lattice_from_h5", "statepoint_groups"]
