# -*- coding: utf-8 -*-
"""Figures for DRAGON5 / Serpent2 comparisons.

All functions build a :class:`matplotlib.figure.Figure` without touching the
``pyplot`` state (no ``plt.show()``, safe in batch / cluster runs) and save it
when *path* is given; the format follows the file extension (png, pdf, eps,
svg). Requires the ``[post]`` extra (matplotlib).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from . import physics_data as phys
from .comparison import PinRateComparison, SpectrumComparison
from .lattice_maps import PinLattice

PathLike = Union[str, Path]

_QUANTITY_LABEL = {
    "absorption": "Absorption rate",
    "fission": "Fission rate",
    "fission_over_absorption": "Fission / absorption",
}


def _figure(**kw):
    try:
        from matplotlib.figure import Figure
    except ImportError as exc:  # pragma: no cover - depends on env
        raise ImportError("matplotlib is required: pip install -e '.[post]'") from exc
    kw.setdefault("layout", "constrained")
    return Figure(**kw)


def new_figure(**kw):
    """A :class:`matplotlib.figure.Figure` (constrained layout) detached from
    ``pyplot``: build custom panels with ``fig.subplots(...)`` and
    :func:`draw_pin_map`, then ``fig.savefig(path)``."""
    return _figure(**kw)


def _save(fig, path: Optional[PathLike], dpi: int = 150):
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return fig


def group_names(n_groups: int) -> list:
    """``["thermal", "fast"]`` for 2 groups, else ``["1", "2", ...]``
    (ascending energy)."""
    return ["thermal", "fast"] if n_groups == 2 else [str(g + 1) for g in range(n_groups)]


def draw_pin_map(ax, lattice: PinLattice, values, pins: Optional[Sequence[int]] = None,
                 cmap: str = "seismic", symmetric: bool = False, fmt: str = "{:.2f}",
                 annotate: bool = True, colorbar: bool = True, title: str = ""):
    """Draw per-pin *values* (aligned with *pins*) on the full lattice."""
    import matplotlib as mpl

    grid = lattice.to_map(values, pins)
    finite = grid[np.isfinite(grid)]
    vmin, vmax = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)
    if symmetric:
        lim = max(abs(vmin), abs(vmax)) or 1.0
        vmin, vmax = -lim, lim
    cm = mpl.colormaps[cmap].with_extremes(bad="white")
    im = ax.imshow(np.ma.masked_invalid(grid), cmap=cm, vmin=vmin, vmax=vmax,
                   interpolation="nearest", origin="upper")
    if annotate:
        mid, span = 0.5 * (vmax + vmin), (vmax - vmin) or 1.0
        for (i, j), v in np.ndenumerate(grid):
            if np.isfinite(v):
                color = "black" if abs(v - mid) < 0.25 * span else "white"
                ax.text(j, i, fmt.format(v), ha="center", va="center",
                        color=color, fontsize=7)
    n_r, n_c = lattice.shape
    ax.set_xticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_r, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", bottom=False, left=False, labelbottom=False,
                   labelleft=False)
    ax.set_title(title)
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def plot_rate_maps(cmp: PinRateComparison, path: Optional[PathLike] = None,
                   method: str = "power", title: Optional[str] = None):
    """Error and distribution maps for every group (the figure of the
    original lattice analysis): absorption error / distribution, and, when
    fission is compared, fission error / distribution and fission-over-
    absorption error."""
    rows = [("absorption", "error"), ("absorption", "distribution")]
    if cmp.has_fission:
        rows += [("fission", "error"), ("fission", "distribution"),
                 ("fission_over_absorption", "error")]
    ng = cmp.n_groups
    names = group_names(ng)
    fig = _figure(figsize=(5.5 * ng, 5 * len(rows)))
    axs = np.atleast_2d(fig.subplots(len(rows), ng, squeeze=False))
    for r, (q, kind) in enumerate(rows):
        data = cmp.errors(q, method) if kind == "error" else cmp.distribution(q)
        for g in range(ng):
            label = "error (%)" if kind == "error" else "distribution"
            draw_pin_map(axs[r, g], cmp.lattice, data[:, g], cmp.pins,
                         symmetric=(kind == "error"),
                         title=f"{_QUANTITY_LABEL[q]} - group {names[g]} - {label}")
    fig.suptitle(title or f"{', '.join(cmp.isotopes)} - DRAGON5 vs Serpent2 "
                 f"({ng} groups, {method} normalisation)", fontsize=16, fontweight="bold")
    return _save(fig, path)


def _shade_merged(axes, merged_eV, label: bool = True) -> None:
    """Grey band over each energy range where DRAGON groups were merged
    (``mesh="common"`` comparisons), so a wide merged bin is not mistaken
    for a physical feature."""
    for ax in axes:
        for n, (lo, hi) in enumerate(merged_eV):
            ax.axvspan(lo, hi, color="0.85", zorder=0,
                       label="merged groups" if (label and n == 0 and ax is axes[0]) else None)


def plot_spectrum(spec: SpectrumComparison, path: Optional[PathLike] = None,
                  per_lethargy: bool = False, with_errors: bool = True,
                  title: Optional[str] = None):
    """Step plot of normalised spectra (one panel per region), with the
    relative difference underneath."""
    n = len(spec.labels)
    rows = 2 if with_errors else 1
    fig = _figure(figsize=(8, 4.5 * n * (1.4 if with_errors else 1)))
    axs = fig.subplots(rows * n, 1, squeeze=False,
                       gridspec_kw={"height_ratios": [3, 1] * n if with_errors else [1] * n})
    d, s = spec.per_lethargy() if per_lethargy else (spec.dragon, spec.serpent)
    e = spec.bounds_eV
    err = spec.errors()
    for k, label in enumerate(spec.labels):
        ax = axs[rows * k, 0]
        ax.step(e, np.append(d[k], d[k][-1]), where="post", label="DRAGON5")
        ax.step(e, np.append(s[k], s[k][-1]), where="post", label="Serpent2")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(e[0], e[-1])
        ax.set_ylabel("flux per unit lethargy" if per_lethargy else "normalised flux")
        ax.set_title(label)
        ax.grid(True, which="both", alpha=0.3)
        _shade_merged([ax], spec.merged_eV)
        ax.legend()
        if with_errors:
            ax2 = axs[rows * k + 1, 0]
            _shade_merged([ax2], spec.merged_eV, label=False)
            ax2.step(e, np.append(err[k], err[k][-1]), where="post", color="tab:green")
            ax2.axhline(0.0, color="black", linewidth=0.8)
            ax2.set_xscale("log")
            ax2.set_xlim(e[0], e[-1])
            ax2.set_ylabel("D5 - S2 (%)")
            ax2.set_xlabel("Energy (eV)")
    fig.suptitle(title or f"Flux spectrum ({len(e) - 1} groups)", fontweight="bold")
    return _save(fig, path)


def plot_profile(cmp: PinRateComparison, pins: Sequence[int], quantity: str = "absorption",
                 path: Optional[PathLike] = None, method: str = "power",
                 title: Optional[str] = None):
    """Rates and relative differences along a list of pins (e.g. the
    symmetry axis, ``lattice.diagonal_pins("axis")``), one column per group.

    Rates are shown normalised to the full-assembly total of each code.
    """
    idx = [cmp.pins.index(p) for p in pins]
    ng = cmp.n_groups
    names = group_names(ng)
    k = cmp.normalisation_factor(method)
    d = cmp.dragon.get(quantity)[idx] * k
    s = cmp.serpent.get(quantity)[idx]
    tot = cmp.serpent.total(quantity) if quantity != "fission_over_absorption" else np.ones(ng)
    err = cmp.errors(quantity, method)[idx]
    labels = [f"{p}\n{cmp.lattice.material.get(p, '')}" for p in pins]
    x = np.arange(len(pins))
    fig = _figure(figsize=(max(6, 0.9 * len(pins)) * ng, 7))
    axs = fig.subplots(2, ng, squeeze=False, sharex="col")
    for g in range(ng):
        ax = axs[0, g]
        ax.bar(x - 0.2, d[:, g] / tot[g], width=0.4, label="DRAGON5")
        ax.bar(x + 0.2, s[:, g] / tot[g], width=0.4, label="Serpent2")
        ax.set_yscale("log")
        ax.set_ylabel(f"{_QUANTITY_LABEL[quantity]} (normalised)")
        ax.set_title(f"group {names[g]}")
        ax.legend()
        ax2 = axs[1, g]
        ax2.bar(x, err[:, g], width=0.6, color="tab:green")
        ax2.axhline(0.0, color="black", linewidth=0.8)
        ax2.set_ylabel("D5 - S2 (%)")
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels, fontsize=8)
    fig.suptitle(title or f"{_QUANTITY_LABEL[quantity]} along pins {list(pins)}",
                 fontweight="bold")
    return _save(fig, path)


def plot_pin_spectra(cmp: PinRateComparison, pins: Sequence[int],
                     quantity: str = "absorption", path: Optional[PathLike] = None,
                     method: str = "power", title: Optional[str] = None,
                     per_lethargy: bool = True):
    """Group-wise rates of selected pins (one panel per pin), each code
    normalised to its full-assembly power, with the relative difference
    underneath. Fine-group counterpart of :func:`plot_profile`.

    Args:
        per_lethargy: Divide the rates by the group lethargy width (default).
            On a non-uniform mesh, and in particular on a ``mesh="common"``
            comparison where several groups were merged into one wide bin,
            per-group rates show bumps that are only bin widths. Merged
            ranges are shaded in grey either way.
    """
    idx = [cmp.pins.index(p) for p in pins]
    k = cmp.normalisation_factor(method)
    p_s = cmp.serpent.total("power").sum()
    d = cmp.dragon.get(quantity)[idx] * k / p_s
    s = cmp.serpent.get(quantity)[idx] / p_s
    e = cmp.bounds_eV
    if per_lethargy:
        du = phys.lethargy_widths(e)
        d, s = d / du, s / du
    err = cmp.errors(quantity, method)[idx]
    merged = getattr(cmp, "merged_eV", [])
    fig = _figure(figsize=(8, 5.5 * len(pins)))
    axs = fig.subplots(2 * len(pins), 1, squeeze=False,
                       gridspec_kw={"height_ratios": [3, 1] * len(pins)})
    ylabel = f"{_QUANTITY_LABEL[quantity]} / power"
    if per_lethargy:
        ylabel += " / unit lethargy"
    for n, p in enumerate(pins):
        ax = axs[2 * n, 0]
        ax.step(e, np.append(d[n], d[n][-1]), where="post", label="DRAGON5")
        ax.step(e, np.append(s[n], s[n][-1]), where="post", label="Serpent2")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(e[0], e[-1])
        ax.set_ylabel(ylabel)
        ax.set_title(f"Pin {p}: {cmp.lattice.material.get(p, '')}")
        ax.grid(True, which="both", alpha=0.3)
        _shade_merged([ax], merged)
        ax.legend()
        ax2 = axs[2 * n + 1, 0]
        _shade_merged([ax2], merged, label=False)
        ax2.step(e, np.append(err[n], err[n][-1]), where="post", color="tab:green")
        ax2.axhline(0.0, color="black", linewidth=0.8)
        ax2.set_xscale("log")
        ax2.set_xlim(e[0], e[-1])
        ax2.set_ylabel("D5 - S2 (%)")
        ax2.set_xlabel("Energy (eV)")
    fig.suptitle(title or f"{_QUANTITY_LABEL[quantity]} per pin ({cmp.n_groups} groups)",
                 fontweight="bold")
    return _save(fig, path)


__all__ = ["draw_pin_map", "plot_rate_maps", "plot_spectrum", "plot_profile",
           "plot_pin_spectra", "group_names", "new_figure"]
