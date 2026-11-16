# -*- coding: utf-8 -*-
"""Pin bookkeeping for code-to-code comparisons.

:class:`PinLattice` links the ``pin_idx`` numbering of a starterDD assembly
model to:

- the **physical positions** of each unique pin in the full lattice (several
  positions share a ``pin_idx`` when the model uses a diagonal symmetry);
- the **multiplicity** of each pin (number of physical positions);
- the **fraction of a physical pin** that a DRAGON mixture represents (1/2 for
  pins cut by a diagonal symmetry axis, 1 otherwise);
- the order of the DRAGON ``MIXTURES`` of a ``by_pin`` EDI output.

Normalisation convention used by the comparison (see
``docs/dev_plans/post_treatment_and_code_to_code_plan.md`` §3.6):

- rates are reported **per physical pin** in both codes;
- a Serpent2 detector that scores a ``pin_idx`` covers every physical
  position of that pin (the universe is reused), so its score is divided by
  :meth:`PinLattice.multiplicity`;
- a DRAGON mixture of a pin on the symmetry axis covers half a pin, so its
  score is divided by :meth:`PinLattice.dragon_fraction`;
- assembly totals are ``sum_p multiplicity(p) * rate(p)``.

This replaces the hardcoded ``DOM_LATTICE`` / ``VAN_LATTICE`` maps,
``rot90`` calls and ``range(1, 11)`` loops of the original GE14 analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


@dataclass
class PinLattice:
    """Physical layout of the unique pins of an assembly.

    Attributes:
        shape: ``(n_rows, n_cols)`` of the full lattice.
        positions: ``{pin_idx: [(row, col), ...]}`` physical positions.
        ordered_pins: Unique ``pin_idx`` in row-major traversal order; the
            DRAGON ``MIXTURES`` of a ``by_pin`` output follow this order.
        symmetry: ``"main-diagonal"``, ``"anti-diagonal"`` or ``None``:
            the symmetry used by the DRAGON geometry.
        material: ``{pin_idx: fuel material name}`` (labels only).
        composition: ``{pin_idx: {isotope: number density}}`` from the model
            (benchmark compositions, used by the legacy analysis).
        fuel_radius: Fuel pellet radius (cm), if known.
    """

    shape: Tuple[int, int]
    positions: Dict[int, List[Tuple[int, int]]]
    ordered_pins: List[int]
    symmetry: Optional[str] = None
    material: Dict[int, str] = field(default_factory=dict)
    composition: Dict[int, Dict[str, float]] = field(default_factory=dict)
    fuel_radius: Optional[float] = None

    # ── construction ───────────────────────────────────────────────────────
    @classmethod
    def from_assembly(cls, assembly) -> "PinLattice":
        """Build from a starterDD ``CartesianAssemblyModel``.

        Prerequisites: ``analyze_lattice_description(build_pins=True)`` and
        ``number_fuel_material_mixtures_by_pin()`` (every fuel pin carries a
        ``pin_idx``).
        """
        from starterDD.DDModel.DragonModel import FuelPinModel

        if assembly.lattice is None:
            raise RuntimeError("Lattice not built: call "
                               "analyze_lattice_description(build_pins=True).")
        info = assembly.get_postprocessing_lattice_info()
        positions: Dict[int, List[Tuple[int, int]]] = {}
        n_rows = len(assembly.lattice)
        n_cols = max(len(r) for r in assembly.lattice)
        for i, row in enumerate(assembly.lattice):
            for j, pin in enumerate(row):
                if isinstance(pin, FuelPinModel):
                    if getattr(pin, "pin_idx", None) is None:
                        raise RuntimeError("Fuel pins have no pin_idx: call "
                                           "number_fuel_material_mixtures_by_pin().")
                    positions.setdefault(int(pin.pin_idx), []).append((i, j))
        lat = cls(shape=(n_rows, n_cols), positions=positions,
                  ordered_pins=list(info["ordered_pin_indices"]),
                  symmetry=info["symmetry_type"],
                  material=dict(info["pin_idx_to_material_name"]),
                  composition={k: dict(v) for k, v in info["pin_idx_to_composition"].items()},
                  fuel_radius=info.get("fuel_radius"))
        lat._check_symmetry()
        return lat

    @classmethod
    def from_inputs(cls, geometry_yaml, compositions_yaml, name: str = "assembly",
                    apply_self_shielding: Optional[str] = "from_yaml",
                    quiet: bool = True) -> "PinLattice":
        """Build from the geometry and composition YAMLs of a run (as
        archived in its ``inputs/`` folder), without a ``DragonCase``.

        The pin numbering is the ``by_pin`` numbering used by the DRAGON
        case. Temperatures do not affect it and are set to placeholders.
        """
        import contextlib
        import io
        from starterDD.DDModel.DragonModel import CartesianAssemblyModel
        from starterDD.DDModel.helpers import associate_material_to_rod_ID
        from starterDD.MaterialProperties.material_mixture import (
            parse_all_compositions_from_yaml,
        )

        sink = io.StringIO() if quiet else None
        with contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext():
            asm = CartesianAssemblyModel(name=name, tdt_file=None,
                                         geometry_description_yaml=str(geometry_yaml))
            asm.set_rod_ID_to_material_mapping(
                associate_material_to_rod_ID(str(compositions_yaml), str(geometry_yaml)))
            asm.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0,
                                         coolant_temperature=600.0,
                                         moderator_temperature=600.0,
                                         structural_temperature=600.0)
            asm.analyze_lattice_description(build_pins=True,
                                            apply_self_shielding=apply_self_shielding)
            asm.set_material_compositions(
                parse_all_compositions_from_yaml(str(compositions_yaml)))
            asm.number_fuel_material_mixtures_by_pin()
        return cls.from_assembly(asm)

    def _check_symmetry(self) -> None:
        """Every off-axis pin must have its mirror; on-axis pins one position."""
        if self.symmetry is None:
            return
        for p, pos in self.positions.items():
            mirrors = {self._mirror(rc) for rc in pos}
            if mirrors != set(pos):
                raise ValueError(f"pin {p}: positions {pos} are not closed under "
                                 f"the {self.symmetry} symmetry.")

    def _mirror(self, rc: Tuple[int, int]) -> Tuple[int, int]:
        i, j = rc
        n = self.shape[0]
        if self.symmetry == "main-diagonal":
            return (j, i)
        if self.symmetry == "anti-diagonal":
            return (n - 1 - j, n - 1 - i)
        return (i, j)

    # ── queries ────────────────────────────────────────────────────────────
    @property
    def pins(self) -> List[int]:
        """Unique pins in mixture order."""
        return list(self.ordered_pins)

    @property
    def n_pins(self) -> int:
        return len(self.ordered_pins)

    @property
    def n_physical_pins(self) -> int:
        return sum(len(v) for v in self.positions.values())

    def multiplicity(self, pin: int) -> int:
        """Number of physical positions of *pin*."""
        return len(self.positions[pin])

    def multiplicities(self, pins: Optional[Sequence[int]] = None) -> np.ndarray:
        return np.array([self.multiplicity(p) for p in (pins or self.ordered_pins)])

    def on_axis(self, pin: int) -> bool:
        """True if *pin* lies on the symmetry axis (cut in half in DRAGON)."""
        return self.symmetry is not None and any(
            self._mirror(rc) == rc for rc in self.positions[pin])

    @property
    def axis_pins(self) -> List[int]:
        return [p for p in self.ordered_pins if self.on_axis(p)]

    def dragon_fraction(self, pin: int) -> float:
        """Fraction of one physical pin represented by the DRAGON mixture of
        *pin* (1/2 on a diagonal symmetry axis, 1 otherwise)."""
        return 0.5 if self.on_axis(pin) else 1.0

    @property
    def dragon_geometry_fraction(self) -> float:
        """Fraction of the full assembly modelled by the DRAGON geometry
        (1/2 with a diagonal symmetry, 1 otherwise). Converts the totals of
        FUEL/ALL outputs to the full assembly modelled by Serpent."""
        return 0.5 if self.symmetry is not None else 1.0

    def dragon_fractions(self, pins: Optional[Sequence[int]] = None) -> np.ndarray:
        return np.array([self.dragon_fraction(p) for p in (pins or self.ordered_pins)])

    def mixture_of(self, pin: int) -> int:
        """1-based DRAGON mixture index of *pin* in a ``by_pin`` output."""
        return self.ordered_pins.index(pin) + 1

    def legacy_symmetry_factors(self) -> np.ndarray:
        """Factors of the original analysis (2 on the axis, 1 elsewhere),
        indexed by ``pin_idx - 1``. Only for reproducing legacy outputs."""
        n = max(self.ordered_pins)
        out = np.ones(n)
        for p in self.axis_pins:
            out[p - 1] = 2.0
        return out

    def diagonal_pins(self, which: str = "axis") -> List[int]:
        """Pins along a diagonal, in traversal order.

        Args:
            which: ``"main"`` (row == col), ``"anti"`` (row + col == n-1) or
                ``"axis"`` (the symmetry axis of the model).
        """
        if which == "axis":
            if self.symmetry is None:
                raise ValueError("The lattice has no symmetry axis.")
            which = "main" if self.symmetry == "main-diagonal" else "anti"
        n = self.shape[0]
        cells = [(i, i) for i in range(n)] if which == "main" else \
            [(i, n - 1 - i) for i in range(n)]
        at = {rc: p for p, pos in self.positions.items() for rc in pos}
        return [at[rc] for rc in cells if rc in at]

    # ── maps ───────────────────────────────────────────────────────────────
    def to_map(self, values: Mapping[int, float] | Sequence[float] | np.ndarray,
               pins: Optional[Sequence[int]] = None, fill: float = np.nan) -> np.ndarray:
        """Unfold per-pin values onto the full lattice.

        Args:
            values: ``{pin: value}`` or an array aligned with *pins*
                (default: :attr:`ordered_pins`). A trailing axis (e.g. energy
                groups) is allowed for arrays shaped ``(n_pins, n_groups)``.
            fill: Value of non-fuel cells.

        Returns:
            ``(n_rows, n_cols)`` or ``(n_groups, n_rows, n_cols)`` array, with
            row 0 the first row of the geometry description.
        """
        pins = list(pins or self.ordered_pins)
        if isinstance(values, Mapping):
            vals = {int(k): np.asarray(v, float) for k, v in values.items()}
        else:
            arr = np.asarray(values, float)
            if arr.shape[0] != len(pins):
                raise ValueError(f"{arr.shape[0]} values for {len(pins)} pins.")
            vals = dict(zip(pins, arr))
        extra = np.shape(next(iter(vals.values())))
        out = np.full(extra + self.shape, fill, dtype=float)
        for p, v in vals.items():
            for (i, j) in self.positions[p]:
                out[(..., i, j)] = v
        return out

    def physical_weights(self, pins: Optional[Sequence[int]] = None) -> np.ndarray:
        """Weights to average per-pin values over physical pins
        (multiplicity / number of physical pins)."""
        m = self.multiplicities(pins)
        return m / m.sum()


__all__ = ["PinLattice"]
