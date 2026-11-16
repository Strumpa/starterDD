# Classes defining a DRAGON calculation scheme for starterDD.
# A calculation scheme orchestrates the sequence of DRAGON calculation steps
# (self-shielding, flux) and their respective geometry
# discretization options (radial subdivision, azimuthal sectorization).
#
# R.Guasch
# Date : 19/02/2026 [created]

import logging
import re
import yaml
from copy import deepcopy

log = logging.getLogger(__name__)


def _plain(obj):
    """Convert *obj* to plain YAML-safe Python types (tuples -> lists,
    numpy scalars/arrays -> floats/lists), recursively."""
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if hasattr(obj, "tolist"):          # numpy scalar or array
        return _plain(obj.tolist())
    return obj


def _deep_merge(base, updates):
    """Return a copy of dict *base* with *updates* merged in. Nested dicts
    are merged key by key; every other value replaces the base value."""
    out = deepcopy(base)
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


#: Relative tolerance used to match branch values (output state points,
#: CLE2000 statepoint guards). CLE2000 REALs are single precision (~6e-8).
STATEPOINT_RTOL = 1.0e-5


def values_match(a, b, rtol=STATEPOINT_RTOL):
    """True if two branch values are equal within *rtol* (absolute *rtol*
    when both are zero)."""
    a, b = float(a), float(b)
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1.0e-30) or (a == b)


# ---------------------------------------------------------------------------
# Statepoint compositions (plan §3.4, P3)
# ---------------------------------------------------------------------------

_ISO_NAME_RE = re.compile(r"^([A-Z][a-z]?)(\d{1,3})$")


def normalise_composition(composition):
    """``{ZAID or isotope name: N}`` -> ``{isotope name: float}``.

    ZAIDs (``"1001"``, ``8016``) become DRAGON/starterDD isotope names
    (``"H1"``, ``"O16"``), as in ``Composition.zaid_to_isotope``. Number
    densities are in atoms/(barn·cm) and must be non-negative.
    """
    from ..MaterialProperties.material_mixture import get_element_symbol
    out = {}
    for key, n in composition.items():
        k = str(key).strip()
        if k.isdigit():
            z, a = divmod(int(k), 1000)
            if a == 0:
                raise ValueError(f"Natural-element ZAID '{k}' is not supported.")
            name = f"{get_element_symbol(z)}{a}"
        elif _ISO_NAME_RE.match(k):
            name = k
        else:
            raise ValueError(f"Cannot read isotope '{k}' (use a ZAID or e.g. 'H1').")
        if name in out:
            raise ValueError(f"Isotope {name} given twice in a composition.")
        n = float(n)
        if n < 0.0:
            raise ValueError(f"Negative number density for {name}: {n}")
        out[name] = n
    return out


def composition_mass_density(composition):
    """Mass density (g/cm³) of a composition in atoms/(barn·cm), from the
    starterDD isotopic masses (``material_mixture.ISOTOPIC_MASSES``)."""
    from ..MaterialProperties.material_mixture import (
        AVOGADRO, CM2_TO_BARN, get_element_symbol, get_isotope_atomic_mass,
    )
    symbol_to_z = {get_element_symbol(z): z for z in range(1, 104)}
    total = 0.0
    for name, n in normalise_composition(composition).items():
        sym, a = _ISO_NAME_RE.match(name).groups()
        total += n * get_isotope_atomic_mass(str(symbol_to_z[sym] * 1000 + int(a)))
    return total * CM2_TO_BARN / AVOGADRO


def read_material_composition(materials_yaml, material="COOLANT"):
    """Isotopic composition of *material* in a starterDD material YAML
    (``MIX_COMPOSITIONS`` entry with an ``isotopic_composition`` block),
    normalised to isotope names."""
    with open(materials_yaml) as f:
        data = yaml.safe_load(f)
    for entry in data.get("MIX_COMPOSITIONS", []):
        if entry.get("name") == material:
            comp = entry.get("isotopic_composition")
            if not comp:
                raise ValueError(
                    f"{materials_yaml}: '{material}' has no isotopic_composition "
                    f"(number densities are required for statepoint compositions).")
            return normalise_composition(comp)
    raise KeyError(f"No material '{material}' in {materials_yaml}.")


def coolant_density_var(isotope):
    """CLE2000 REAL variable holding the coolant number density of
    *isotope* at the current statepoint (``N_H1``, ``N_O16``)."""
    name = f"N_{isotope}"
    if len(f"{name}_arr") > 12:
        raise ValueError(f"CLE2000 name '{name}_arr' for {isotope} exceeds 12 characters.")
    return name


def coolant_density_array(isotope):
    """CLE2000 ``UTL:`` array of the coolant number densities of *isotope*
    over the coolant_density branch (``N_H1_arr``)."""
    return f"{coolant_density_var(isotope)}_arr"


# ---------------------------------------------------------------------------
# Valid option constants
# ---------------------------------------------------------------------------

VALID_STEP_TYPES = ("self_shielding", "flux", "edition_between_levels")
VALID_SELF_SHIELDING_MODULES = ("USS", "SHI")  # USS: unresolved resonance, SHI: Stammler
VALID_SELF_SHIELDING_METHODS = ("RSE", "PT")  # RSE: subgroup+equivalence, PT: probability tables
VALID_SPATIAL_METHODS = ("CP", "IC", "MOC")
VALID_TRACKING_OPTIONS = ("TISO", "TSPC")
VALID_RADIAL_SCHEMES = ("Santamarina", "Santamarina_fine_Gd", "automatic", "user_defined")
VALID_MIX_NUMBERING_STRATEGIES = ("by_material", "by_pin")
# FLU: leakage options (IGE335 §3.7, descleak structure). Only the TYPE K + B1 + SIGS
# combination is enabled for now; extend these tuples to enable more models.
VALID_LEAKAGE_MODELS = ("B1",)        # descleak group 1: leakage coefficients model
VALID_LEAKAGE_CORRECTIONS = ("SIGS",) # descleak group 2: how the leakage correction is applied


# ---------------------------------------------------------------------------
# SectorConfig – lightweight container for sectorization parameters
# ---------------------------------------------------------------------------

class SectorConfig:
    """
    Container for discretization parameters for cell regions.

    For **circular** water rods and fuel pins this stores azimuthal
    sectorization parameters used by glow's ``RectCell.sectorize()``
    method (``sectors``, ``angles``, ``windmill``).

    For **square** water rods this stores a Cartesian grid subdivision
    specification (``splits``) applied over the bounding box, analogous
    to the assembly-box discretization.

    Attributes
    ----------
    sectors : list[int]
        Number of sectors per radial ring (from innermost to outermost
        region, including coolant).  Used for circular geometries.
    angles : list[float]
        Starting angle offset for each ring (degrees).  Used for circular
        geometries.
    windmill : bool
        Whether to apply windmill sectorization.  Used for circular
        geometries.
    splits : tuple[int, int] or None
        ``(nx, ny)`` Cartesian grid subdivisions applied to the whole
        bounding box of a square water rod.  Material assignment to each
        resulting sub-face is performed by geometric containment.
        ``None`` means no grid sub-meshing (the cell keeps its base
        3-region technological geometry).
    additional_radial_splits_in_moderator : int or list[float]
        Controls radial subdivision of the innermost moderator region
        of a **circular** water rod (the region with ``r < inner_radius``).

        * ``int`` **N** (default ``1``): N = 1 means no extra splits.
          N ≥ 2 produces N − 1 equally-spaced intermediate radii at
          ``inner_radius × k / N`` for ``k = 1 … N − 1``.
        * ``list[float]``: explicit radii (must all be > 0 and
          < ``inner_radius``; validated at geometry-build time).

        When extra rings are added, ``sectors[0]`` and ``angles[0]``
        are automatically replicated for every new sub-ring so that
        the user only needs to specify sectors/angles for the **base**
        regions (moderator, clad, coolant).
    """

    VALID_SECTOR_NO_TO_ANGLE = {1: [0], 4: [0, 45], 8: [0, 22.5], 16: [0]}

    def __init__(self, sectors=None, angles=None, windmill=False,
                 splits=None, additional_radial_splits_in_moderator=None,
                 subdivisions_coolant_corners=None):
        self.sectors = sectors or []
        self.angles = angles or []
        self.windmill = windmill
        self.splits = tuple(splits) if splits else None
        self.additional_radial_splits_in_moderator = (
            additional_radial_splits_in_moderator
        )
        self.subdivisions_coolant_corners = subdivisions_coolant_corners
        if sectors is not None and angles is not None:
            self._validate_sector_to_angles()

    def _validate_sector_to_angles(self):
        """
        Validate that the number of sectors matches the number of angles
        and that the combinations of number of sectors and angles are supported by GLOW.

        Raises
        ------
        ValueError
            If the number of sectors does not match the number of angles
            and is not in the predefined mapping.
        """
        if len(self.sectors) != len(self.angles):
            raise ValueError(
                f"Number of sector rings {len(self.sectors)} does not match "
                f"number of angle entries {len(self.angles)}. Either they must "
                f"match or the number of sectors must be in the predefined mapping."
            )
        for i, n in enumerate(self.sectors):
            if n not in self.VALID_SECTOR_NO_TO_ANGLE:
                raise ValueError(
                    f"Number of sectors {n} does not match the number of sectors supported by GLOW"
                    f"The valid sector to angles mapping is: {self.VALID_SECTOR_NO_TO_ANGLE}"
                )
            elif self.angles[i] not in self.VALID_SECTOR_NO_TO_ANGLE[n]:
                raise ValueError(
                    f"Angle {self.angles[i]} is not valid for {n} sectors. "
                    f"The valid angles for {n} sectors are: {self.VALID_SECTOR_NO_TO_ANGLE[n]}"
                )


    # ------------------------------------------------------------------
    # Radial-split helpers (circular water rods)
    # ------------------------------------------------------------------

    def resolve_water_rod_radii(self, inner_radius):
        """
        Compute intermediate radii for moderator sub-division.

        Parameters
        ----------
        inner_radius : float
            The base inner radius of the circular water rod (boundary
            between moderator and cladding).

        Returns
        -------
        list[float]
            Sorted list of intermediate radii (all < ``inner_radius``).
            Empty list when no additional splits are requested.

        Raises
        ------
        ValueError
            If user-supplied radii are outside ``(0, inner_radius)``.
        """
        spec = self.additional_radial_splits_in_moderator

        if spec is None:
            return []

        # Integer mode
        if isinstance(spec, int):
            if spec <= 1:
                return []
            return [inner_radius * k / spec for k in range(1, spec)]

        # List-of-radii mode
        if isinstance(spec, (list, tuple)):
            if len(spec) == 0:
                return []
            radii = sorted(float(r) for r in spec)
            for r in radii:
                if r <= 0.0 or r >= inner_radius:
                    raise ValueError(
                        f"additional_radial_splits_in_moderator: radius "
                        f"{r} is out of range (0, {inner_radius})."
                    )
            return radii

        raise TypeError(
            f"additional_radial_splits_in_moderator must be int or "
            f"list[float], got {type(spec).__name__}"
        )

    def expanded_sectors_and_angles(self, inner_radius):
        """
        Return ``(sectors, angles)`` lists expanded to account for
        additional moderator sub-rings.

        The first entry in ``self.sectors`` / ``self.angles`` corresponds
        to the original moderator region.  For each extra sub-ring the
        same value is prepended so that the resulting lists match the
        total number of radial regions in the geometry.

        Parameters
        ----------
        inner_radius : float
            Base inner radius (used to compute the number of extra rings).

        Returns
        -------
        (list[int], list[float])
            Expanded sectors and angles lists.
        """
        extra_radii = self.resolve_water_rod_radii(inner_radius)
        n_extra = len(extra_radii)

        if n_extra == 0 or len(self.sectors) == 0:
            return list(self.sectors), list(self.angles)

        base_sector = self.sectors[0]
        base_angle = self.angles[0] if self.angles else 0.0

        expanded_s = [base_sector] * n_extra + list(self.sectors)
        expanded_a = [base_angle] * n_extra + list(self.angles)

        return expanded_s, expanded_a

    def __repr__(self):
        if self.splits is not None:
            return f"SectorConfig(splits={self.splits})"
        extra = ""
        spec = self.additional_radial_splits_in_moderator
        if spec is not None and spec != 1:
            extra = f", additional_radial_splits={spec}"
        return (f"SectorConfig(sectors={self.sectors}, "
                f"angles={self.angles}, windmill={self.windmill}{extra})")

    def to_dict(self):
        """YAML entry of a ``sectorization`` sub-block (``fuel_pins``,
        ``Gd_pins``, ``water_rods``). ``windmill`` lives one level up."""
        d = {"sectors": list(self.sectors), "angles": list(self.angles)}
        if self.splits is not None:
            d["splits"] = list(self.splits)
        if self.additional_radial_splits_in_moderator is not None:
            d["additional_radial_splits_in_moderator"] = (
                self.additional_radial_splits_in_moderator)
        if self.subdivisions_coolant_corners is not None:
            d["subdivisions_coolant_corners"] = self.subdivisions_coolant_corners
        return _plain(d)


# ---------------------------------------------------------------------------
# ControlCrossSubmeshConfig – sub-meshing options for the control cross wings
# ---------------------------------------------------------------------------

class ControlCrossSubmeshConfig:
    """
    Configuration for sub-meshing the control cross wing regions.

    Each wing arm is decomposed into three axial zones (perpendicular to
    the arm axis):

    1. **Control-cross corner zone** — the ``bt/2 × bt/2`` square where
       both arms overlap at the cross centre.
    2. **Central-structure-to-absorber zone** — the span from the corner
       zone edge to the first absorber tube boundary.
    3. **Absorber-pin zone** — the span containing the absorber tubes
       up to the wing tip.

    Additional options control whether tube bounding surfaces are
    extended to the sheath border and whether tubes are bisected along
    the arm axis.

    Attributes
    ----------
    enabled : bool
        Whether control-cross sub-meshing is active.
    control_cross_corner_splits : tuple[int, int] or None
        ``(n_along_arm, n_across_arm)`` grid for the control-cross
        corner zone (the ``bt/2 × bt/2`` square at the cross centre).
        ``None`` → ``(1, 1)`` (no sub-division).
    central_structure_splits : tuple[int, int] or None
        ``(n_along_arm, n_across_arm)`` grid for the CS-to-absorber zone.
        ``None`` → ``(1, 1)``.
    extend_splits_at_tube_boundaries : bool
        When ``True``, generate ``bt``-wide splitting faces at each
        absorber tube centre so that tube bounding surfaces extend to
        the sheath border.  Default ``True``.
    split_tubes_in_half : bool
        When ``True``, add a bisecting split at each tube centre along
        the arm axis.  Default ``False``.
    """

    # Default values used when ``control_cross_submesh: true`` (bool shorthand)
    _DEFAULTS = dict(
        enabled=True,
        control_cross_corner_splits=None,
        central_structure_splits=None,
        extend_splits_at_tube_boundaries=True,
        split_tubes_in_half=False,
    )

    def __init__(self, enabled=True, control_cross_corner_splits=None,
                 central_structure_splits=None,
                 extend_splits_at_tube_boundaries=True,
                 split_tubes_in_half=False):
        self.enabled = enabled
        self.control_cross_corner_splits = (
            tuple(control_cross_corner_splits)
            if control_cross_corner_splits else None
        )
        self.central_structure_splits = (
            tuple(central_structure_splits)
            if central_structure_splits else None
        )
        self.extend_splits_at_tube_boundaries = extend_splits_at_tube_boundaries
        self.split_tubes_in_half = split_tubes_in_half

    @classmethod
    def from_yaml(cls, raw):
        """
        Parse a ``control_cross_submesh`` YAML value.

        Accepts:

        - ``bool``: ``True`` → default config, ``False`` → disabled.
        - ``dict``: explicit configuration keys.

        Parameters
        ----------
        raw : bool or dict
            The raw YAML value.

        Returns
        -------
        ControlCrossSubmeshConfig or None
            ``None`` when disabled.
        """
        if isinstance(raw, bool):
            return cls(**cls._DEFAULTS) if raw else None
        if isinstance(raw, dict):
            return cls(
                enabled=raw.get("enabled", True),
                control_cross_corner_splits=raw.get(
                    "control_cross_corner_splits", None
                ),
                central_structure_splits=raw.get(
                    "central_structure_splits", None
                ),
                extend_splits_at_tube_boundaries=raw.get(
                    "extend_splits_at_tube_boundaries", True
                ),
                split_tubes_in_half=raw.get("split_tubes_in_half", False),
            )
        raise TypeError(
            f"control_cross_submesh must be a bool or dict, got "
            f"{type(raw).__name__}"
        )

    def __repr__(self):
        return (
            f"ControlCrossSubmeshConfig(enabled={self.enabled}, "
            f"control_cross_corner_splits="
            f"{self.control_cross_corner_splits}, "
            f"central_structure_splits={self.central_structure_splits}, "
            f"extend_tube_boundaries={self.extend_splits_at_tube_boundaries}, "
            f"split_tubes_in_half={self.split_tubes_in_half})"
        )

    def to_dict(self):
        """YAML ``control_cross_submesh`` dict (inverse of :meth:`from_yaml`)."""
        return _plain({
            "enabled": self.enabled,
            "control_cross_corner_splits": self.control_cross_corner_splits,
            "central_structure_splits": self.central_structure_splits,
            "extend_splits_at_tube_boundaries": self.extend_splits_at_tube_boundaries,
            "split_tubes_in_half": self.split_tubes_in_half,
        })


# ---------------------------------------------------------------------------
# CrossModeratorDiscretizationConfig – sub-meshing options for moderator
# regions surrounding the control cross
# ---------------------------------------------------------------------------

class CrossModeratorDiscretizationConfig:
    """
    Sub-meshing options for the moderator regions surrounding a control
    cross device, used within ``BoxDiscretizationConfig``.

    All split fields are optional.  When ``None``, splits are computed
    automatically from the ``gap_splits`` density applied to the actual
    dimensions of each region.

    Attributes
    ----------
    cross_side_gap_splits : tuple[int, int] or None
        ``(n_parallel, n_perpendicular)`` grid for the moderator gap
        on the sides affected by the control cross (between the blade edge
        and the pin-lattice footprint edge). This overrides general
        ``gap_narrow_splits`` on affected sides only.
        Orientation follows the same convention as ``gap_splits``:
        ``n_parallel`` is along the lattice edge, ``n_perpendicular``
        is across the gap width.  Auto-permuted for horizontal vs
        vertical affected sides.
    moderator_at_cross_corner_splits : tuple[int, int] or None
        ``(nx, ny)`` grid for the small moderator gap rectangle between
        the two cross arms and the pin-lattice corner.
    stub_splits : tuple[int, int] or None
        ``(n_parallel, n_perpendicular)`` grid for the moderator stubs
        that sit in the blade-thickness zone beyond the arm tip extent.
        Auto-permuted for horizontal vs vertical stubs.
    """

    def __init__(self, cross_side_gap_splits=None,
                 moderator_at_cross_corner_splits=None,
                 stub_splits=None,
                 # Backward compatibility
                 narrow_gap_splits=None):
        import warnings

        # Handle backward compatibility for narrow_gap_splits
        if narrow_gap_splits is not None and cross_side_gap_splits is None:
            warnings.warn(
                "Parameter 'narrow_gap_splits' is deprecated. "
                "Use 'cross_side_gap_splits' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            cross_side_gap_splits = narrow_gap_splits

        self.cross_side_gap_splits = (
            tuple(cross_side_gap_splits) if cross_side_gap_splits else None
        )
        self.moderator_at_cross_corner_splits = (
            tuple(moderator_at_cross_corner_splits)
            if moderator_at_cross_corner_splits else None
        )
        self.stub_splits = (
            tuple(stub_splits) if stub_splits else None
        )

    def resolve(self, gap_splits, lattice_pitch, wide_gap_width,
                narrow_gap_width, cross_corner_dims, stub_dims):
        """
        Return resolved ``(cross_side_gap, moderator_corner, stub)`` split
        tuples, auto-computing from mesh density when a field is ``None``.

        The reference density is derived from ``gap_splits`` applied to a
        strip of size ``lattice_pitch × wide_gap_width``:

        - ``density_parallel  = gap_splits[0] / lattice_pitch``
        - ``density_perpendicular = gap_splits[1] / wide_gap_width``

        Parameters
        ----------
        gap_splits : tuple[int, int]
            ``(n_parallel, n_perpendicular)`` for unaffected gap strips.
        lattice_pitch : float
            Parallel extent of an unaffected gap strip (lattice edge).
        wide_gap_width : float
            Perpendicular extent of an unaffected gap strip.
        narrow_gap_width : float
            Perpendicular extent of the narrow gap (blade edge to
            lattice edge).
        cross_corner_dims : tuple[float, float]
            ``(width, height)`` of the moderator corner gap rectangle.
        stub_dims : tuple[float, float]
            ``(parallel_length, perpendicular_length)`` of a typical
            stub region.

        Returns
        -------
        cross_side_gap : tuple[int, int]
            ``(n_parallel, n_perpendicular)`` for the affected side gap.
        moderator_corner : tuple[int, int]
            ``(nx, ny)`` for the moderator corner rectangle.
        stub : tuple[int, int]
            ``(n_parallel, n_perpendicular)`` for stub regions.
        """
        d_par = gap_splits[0] / lattice_pitch if lattice_pitch > 0 else 1.0
        d_perp = (
            gap_splits[1] / wide_gap_width if wide_gap_width > 0 else 1.0
        )

        # Cross-affected side gap
        if self.cross_side_gap_splits is not None:
            cross_side_gap = self.cross_side_gap_splits
        else:
            cross_side_gap = (
                max(1, round(d_par * lattice_pitch)),
                max(1, round(d_perp * narrow_gap_width)),
            )

        # Moderator corner
        if self.moderator_at_cross_corner_splits is not None:
            moderator_corner = self.moderator_at_cross_corner_splits
        else:
            moderator_corner = (
                max(1, round(d_perp * cross_corner_dims[0])),
                max(1, round(d_perp * cross_corner_dims[1])),
            )

        # Stubs
        if self.stub_splits is not None:
            stub = self.stub_splits
        else:
            stub = (
                max(1, round(d_par * stub_dims[0])),
                max(1, round(d_perp * stub_dims[1])),
            )

        return cross_side_gap, moderator_corner, stub

    def __repr__(self):
        return (
            f"CrossModeratorDiscretizationConfig("
            f"cross_side_gap_splits={self.cross_side_gap_splits}, "
            f"moderator_at_cross_corner_splits="
            f"{self.moderator_at_cross_corner_splits}, "
            f"stub_splits={self.stub_splits})"
        )

    def to_dict(self):
        """YAML ``cross_moderator_discretization`` dict."""
        return _plain({
            "cross_side_gap_splits": self.cross_side_gap_splits,
            "moderator_at_cross_corner_splits": self.moderator_at_cross_corner_splits,
            "stub_splits": self.stub_splits,
        })


# ---------------------------------------------------------------------------
# BoxDiscretizationConfig – sub-meshing of the assembly box for MOC
# ---------------------------------------------------------------------------

class BoxDiscretizationConfig:
    """
    Configuration for subdividing the assembly-box peripheral regions
    into a grid of sub-faces, typically needed for MOC tracking.

    The assembly box is divided into 8 strips surrounding the pin-lattice
    footprint (4 corners + 4 sides).  Each strip category can be
    independently sub-meshed:

    * **corner_splits** ``(nx, ny)`` — grid size for each corner region.
    * **gap_splits** ``(n_parallel, n_perpendicular)`` — grid size for
      side strips.  ``n_parallel`` is along the lattice edge,
      ``n_perpendicular`` is across the gap width.  The tuple is
      automatically permuted for horizontal vs vertical strips:
      horizontal strips use ``(n_par, n_perp)`` as ``(nx, ny)``,
      vertical strips use ``(n_perp, n_par)``.

    For asymmetric gap configurations (gap_wide ≠ gap_narrow), the following
    optional parameters allow region-specific discretization:

    * **gap_wide_splits** ``(n_parallel, n_perpendicular)`` — grid size for
      regions bordering the wide gap (applies to sides and wide-gap corners).
    * **gap_narrow_splits** ``(n_parallel, n_perpendicular)`` — grid size for
      regions bordering the narrow gap (applies to sides and narrow-gap corners).
    * **wide_wide_corner_splits**, **narrow_narrow_corner_splits**,
      **mixed_corner_splits** — individual corner-specific splits for precise control.

    When a control cross is present, two optional sibling configs control
    the sub-meshing of the surrounding moderator regions and the cross
    structure itself:

    * ``cross_moderator_discretization`` — split parameters for the
      moderator regions affected by the cross (narrow gap, moderator
      corner, stubs).
    * ``control_cross_submesh`` — sub-meshing of the control cross
      wings and corner zone.

    Attributes
    ----------
    enabled : bool
        Whether box discretization is active.
    corner_splits : tuple[int, int]
        ``(nx, ny)`` grid subdivisions for each corner region (default when
        specific corner splits are not set).
    gap_splits : tuple[int, int] or None
        ``(n_parallel, n_perpendicular)`` grid subdivisions for
        gap side strips.  Defaults to ``(n_cols, 1)`` when ``None``.
    gap_wide_splits : tuple[int, int] or None
        Grid subdivisions for wide-gap regions. Overrides gap_splits for
        regions bordering gap_wide. Only applicable with asymmetric gaps.
    gap_narrow_splits : tuple[int, int] or None
        Grid subdivisions for narrow-gap regions. Overrides gap_splits for
        regions bordering gap_narrow. Only applicable with asymmetric gaps.
    wide_wide_corner_splits : tuple[int, int] or None
        Grid subdivisions for the corner with (gap_wide, gap_wide) on both axes.
        Overrides gap_wide_splits if set.
    narrow_narrow_corner_splits : tuple[int, int] or None
        Grid subdivisions for the corner with (gap_narrow, gap_narrow) on both axes.
        Overrides gap_narrow_splits if set.
    mixed_corner_splits : tuple[int, int] or None
        Grid subdivisions for mixed corners (one axis wide, one axis narrow).
    cross_moderator_discretization : CrossModeratorDiscretizationConfig or None
        Optional sub-meshing config for moderator regions surrounding
        the control cross.
    control_cross_submesh : ControlCrossSubmeshConfig or None
        Optional sub-meshing config for the control cross wings and
        corner zone.
    reassign_materials : bool
        When ``True`` (default), after the discretization partition an
        explicit material re-assignment pass is performed using
        geometric containment against the original reference shapes
        (coolant boundary, channel box boundary, control cross shapes).
        This guarantees correct materials on every sub-face regardless
        of how glow's internal ``update_geometry_from_face`` propagates
        properties across successive partitions.  Set to ``False`` to
        skip this step and rely solely on property propagation (faster
        but may produce incorrect materials on control-cross sub-faces).
    """

    def __init__(self, enabled=False, corner_splits=None,
                 gap_splits=None, gap_wide_splits=None, gap_narrow_splits=None,
                 wide_wide_corner_splits=None, narrow_narrow_corner_splits=None,
                 mixed_corner_splits=None,
                 cross_moderator_discretization=None,
                 control_cross_submesh=None,
                 reassign_materials=True,
                 # deprecated aliases kept for backward compatibility
                 side_x_splits=None, side_y_splits=None):
        import warnings

        self.enabled = enabled
        self.corner_splits = tuple(corner_splits) if corner_splits else (4, 4)

        # Asymmetric gap split parameters (new)
        self.gap_wide_splits = tuple(gap_wide_splits) if gap_wide_splits else None
        self.gap_narrow_splits = tuple(gap_narrow_splits) if gap_narrow_splits else None
        self.wide_wide_corner_splits = tuple(wide_wide_corner_splits) if wide_wide_corner_splits else None
        self.narrow_narrow_corner_splits = tuple(narrow_narrow_corner_splits) if narrow_narrow_corner_splits else None
        self.mixed_corner_splits = tuple(mixed_corner_splits) if mixed_corner_splits else None

        self.cross_moderator_discretization = cross_moderator_discretization
        self.control_cross_submesh = control_cross_submesh
        self.reassign_materials = reassign_materials

        # Handle deprecated aliases
        if gap_splits is not None:
            self.gap_splits = tuple(gap_splits)
        elif side_x_splits is not None:
            warnings.warn(
                "'side_x_splits' / 'side_y_splits' are deprecated; "
                "use 'gap_splits' instead.  'side_x_splits' is being "
                "used as gap_splits = (n_parallel, n_perpendicular).",
                DeprecationWarning,
                stacklevel=2,
            )
            self.gap_splits = tuple(side_x_splits)
        else:
            self.gap_splits = None  # will be defaulted in resolve_splits

    def resolve_splits(self, n_cols, n_rows):
        """
        Return the resolved ``(corner, side_h, side_v)`` split tuples,
        filling in lattice-dimension defaults for ``None`` values.

        ``gap_splits`` is auto-permuted:

        - Horizontal strips (top / bottom): ``(n_par, n_perp)``
        - Vertical strips (left / right): ``(n_perp, n_par)``

        Parameters
        ----------
        n_cols : int
            Number of pin columns in the lattice.
        n_rows : int
            Number of pin rows in the lattice.

        Returns
        -------
        corner : tuple[int, int]
        side_h : tuple[int, int]
            Splits for horizontal (top/bottom) strips.
        side_v : tuple[int, int]
            Splits for vertical (left/right) strips (permuted).
        """
        corner = self.corner_splits
        gap = self.gap_splits if self.gap_splits else (n_cols, 1)
        side_h = gap                   # (n_parallel, n_perpendicular)
        side_v = (gap[1], gap[0])      # permuted for vertical strips
        return corner, side_h, side_v

    def resolve_splits_with_symmetry(self, n_cols, n_rows, assembly_model):
        """
        Resolve region-specific splits for all 8 box regions (4 corners + 4 sides)
        based on detected lattice symmetry and gap_wide vs gap_narrow configuration.

        For assemblies with asymmetric gaps (gap_wide ≠ gap_narrow), corners and
        sides are classified based on symmetry type:

        - **Main-diagonal symmetry** (BL=wide-wide, TR=narrow-narrow):
          Grid structure: left=wide, right=narrow, bottom=wide, top=narrow

        - **Anti-diagonal symmetry** (TL=wide-wide, BR=narrow-narrow):
          Grid structure: left=wide, right=narrow, bottom=narrow, top=wide

        - **No symmetry**: All regions treated uniformly (backward compatible)

        **Mesh Alignment Convention:**

        - **wide_wide and narrow_narrow corners** keep fixed splits (nx, ny) — their mesh
          orientation is determined by which gaps they border, not affected by symmetry.

        - **Mixed corners** have their splits transposed based on symmetry to maintain
          consistent mesh alignment relative to the gap configuration:

          * Main-diagonal: TL mixed corner gets transposed splits for consistent orientation
          * Anti-diagonal: BL mixed corner gets transposed splits for consistent orientation

        This ensures the mesh refinement pattern aligns consistently with the detected
        gap structure and maintains geometric symmetry.

        Parameters
        ----------
        n_cols : int
            Number of pin columns in the lattice.
        n_rows : int
            Number of pin rows in the lattice.
        assembly_model : CartesianAssemblyModel
            Assembly model providing symmetry information via check_diagonal_symmetry().

        Returns
        -------
        dict
            A dictionary with keys for the 8 regions:
            - 'corner_bl', 'corner_br', 'corner_tl', 'corner_tr': corners
            - 'side_top', 'side_bottom', 'side_left', 'side_right': sides
            Each maps to a tuple (nx, ny) of grid subdivisions.
        """
        # Get symmetry from assembly model
        symmetry = assembly_model.check_diagonal_symmetry()

        # Helper to resolve a single region's splits
        def resolve_region_splits(gap_type):
            """
            Resolve splits for a region based on its gap type:
            - 'wide_wide': use wide_wide_corner_splits else gap_wide_splits else gap_splits
            - 'narrow_narrow': use narrow_narrow_corner_splits else gap_narrow_splits else gap_splits
            - 'mixed': use mixed_corner_splits else gap_splits
            - 'side_wide': use gap_wide_splits else gap_splits
            - 'side_narrow': use gap_narrow_splits else gap_splits
            """
            if gap_type == 'wide_wide' and self.wide_wide_corner_splits:
                return self.wide_wide_corner_splits
            elif gap_type == 'narrow_narrow' and self.narrow_narrow_corner_splits:
                return self.narrow_narrow_corner_splits
            elif gap_type == 'mixed' and self.mixed_corner_splits:
                return self.mixed_corner_splits
            elif gap_type in ('wide_wide', 'side_wide') and self.gap_wide_splits:
                return self.gap_wide_splits
            elif gap_type in ('narrow_narrow', 'side_narrow') and self.gap_narrow_splits:
                return self.gap_narrow_splits
            else:
                # Fall back to gap_splits
                gap = self.gap_splits if self.gap_splits else (n_cols, 1)
                if gap_type in ('side_wide', 'side_narrow'):
                    # For sides, return horizontally oriented (will be permuted by caller if needed)
                    return gap
                elif gap_type in ('mixed', 'wide_wide', 'narrow_narrow'):
                    # For corners, return as-is
                    return self.corner_splits
                return gap

        # Initialize result dictionary
        regions = {}

        # Helper to transpose (nx, ny) → (ny, nx)
        def transpose_splits(splits):
            """Transpose a split tuple (nx, ny) → (ny, nx) for symmetry."""
            if isinstance(splits, (tuple, list)) and len(splits) == 2:
                return (splits[1], splits[0])
            return splits

        # Classify corners and sides based on symmetry
        if symmetry == "main-diagonal":
            # Main-diagonal: gap_wide on left & bottom, gap_narrow on right & top
            # BL (bottom-left): both wide — no transpose (fixed by gap configuration)
            regions['corner_bl'] = resolve_region_splits('wide_wide')
            # BR (bottom-right): mixed (bottom=wide, right=narrow) — no transpose
            regions['corner_br'] = resolve_region_splits('mixed')
            # TL (top-left): mixed (left=wide, top=narrow) — transpose for consistent mesh orientation
            regions['corner_tl'] = transpose_splits(resolve_region_splits('mixed'))
            # TR (top-right): both narrow — no transpose (fixed by gap configuration)
            regions['corner_tr'] = resolve_region_splits('narrow_narrow')
            # Sides
            regions['side_bottom'] = resolve_region_splits('side_wide')
            regions['side_top'] = resolve_region_splits('side_narrow')
            regions['side_left'] = resolve_region_splits('side_wide')
            regions['side_right'] = resolve_region_splits('side_narrow')

        elif symmetry == "anti-diagonal":
            # Anti-diagonal: gap_wide on left & top, gap_narrow on right & bottom
            # TL (top-left): both wide — no transpose (fixed by gap configuration)
            regions['corner_tl'] = resolve_region_splits('wide_wide')
            # TR (top-right): mixed (top=wide, right=narrow) — no transpose
            regions['corner_tr'] = resolve_region_splits('mixed')
            # BL (bottom-left): mixed (left=wide, bottom=narrow) — transpose for consistent mesh orientation
            regions['corner_bl'] = transpose_splits(resolve_region_splits('mixed'))
            # BR (bottom-right): both narrow — no transpose (fixed by gap configuration)
            regions['corner_br'] = resolve_region_splits('narrow_narrow')
            # Sides
            regions['side_top'] = resolve_region_splits('side_wide')
            regions['side_bottom'] = resolve_region_splits('side_narrow')
            regions['side_left'] = resolve_region_splits('side_wide')
            regions['side_right'] = resolve_region_splits('side_narrow')

        else:
            # No symmetry detected: use uniform splits (backward compatible)
            uniform_corner = self.corner_splits
            uniform_gap = self.gap_splits if self.gap_splits else (n_cols, 1)
            for corner_key in ['corner_bl', 'corner_br', 'corner_tl', 'corner_tr']:
                regions[corner_key] = uniform_corner
            for side_key in ['side_top', 'side_bottom', 'side_left', 'side_right']:
                regions[side_key] = uniform_gap

        # Apply horizontal/vertical permutation for sides
        # Horizontal sides (top/bottom): use (n_parallel, n_perpendicular) as (nx, ny)
        # Vertical sides (left/right): permute to (n_perpendicular, n_parallel)
        side_top_h = regions['side_top']  # horizontal strip
        side_bottom_h = regions['side_bottom']  # horizontal strip
        side_left_v = regions['side_left']  # vertical strip - permute for vertical
        side_right_v = regions['side_right']  # vertical strip - permute for vertical

        # Apply permutation to vertical sides
        regions['side_top'] = side_top_h
        regions['side_bottom'] = side_bottom_h
        regions['side_left'] = (side_left_v[1], side_left_v[0]) if len(side_left_v) == 2 else side_left_v
        regions['side_right'] = (side_right_v[1], side_right_v[0]) if len(side_right_v) == 2 else side_right_v

        return regions

    def __repr__(self):
        return (f"BoxDiscretizationConfig(enabled={self.enabled}, "
                f"corner_splits={self.corner_splits}, "
                f"gap_splits={self.gap_splits}, "
                f"gap_wide_splits={self.gap_wide_splits}, "
                f"gap_narrow_splits={self.gap_narrow_splits}, "
                f"wide_wide_corner_splits={self.wide_wide_corner_splits}, "
                f"narrow_narrow_corner_splits={self.narrow_narrow_corner_splits}, "
                f"mixed_corner_splits={self.mixed_corner_splits}, "
                f"cross_moderator_discretization="
                f"{self.cross_moderator_discretization}, "
                f"control_cross_submesh={self.control_cross_submesh}, "
                f"reassign_materials={self.reassign_materials})")

    def to_dict(self):
        """YAML ``box_discretization`` dict."""
        d = {
            "enabled": self.enabled,
            "corner_splits": self.corner_splits,
            "gap_splits": self.gap_splits,
            "gap_wide_splits": self.gap_wide_splits,
            "gap_narrow_splits": self.gap_narrow_splits,
            "wide_wide_corner_splits": self.wide_wide_corner_splits,
            "narrow_narrow_corner_splits": self.narrow_narrow_corner_splits,
            "mixed_corner_splits": self.mixed_corner_splits,
            "reassign_materials": self.reassign_materials,
            "control_cross_submesh": (self.control_cross_submesh.to_dict()
                                      if self.control_cross_submesh is not None
                                      else False),
        }
        if self.cross_moderator_discretization is not None:
            d["cross_moderator_discretization"] = (
                self.cross_moderator_discretization.to_dict())
        return _plain(d)

# -------------------------------------------------------------------------------
# VanishedRodDiscretizationConfig – sub-meshing options for vanished rod regions
# -------------------------------------------------------------------------------


class VanishedRodDiscretizationConfig:
    """
    Sub-meshing options for the vanished rod regions in the lattice
    
    When a vanished rod is present, the region originally occupied by the rod is filled with coolant,
    proper discretization of this region is required to capture thermalization effects and ensure an accurate flux distribution is obtained.

    Attributes
    ----------
    base_radius : float (optional)
        Base radius, used to determine the extent of the sub-meshing region, if None, defaults to outer clad radius from PIN_GEOMETRY.
    radial_splits : int
        Number of radial splits to apply within the vanished rod region, default is 1.
    sector : int
        Number of azimuthal sector splits to apply within the vanished rod region, default is 1.
    angles : list[float] or None
        Optional list of angles for sector splits, in degrees.  If None, sectors are evenly spaced.
    """

    def __init__(self, base_radius=None, radial_splits=1, sectors=None, angles=None, windmill=False):
        self.base_radius = base_radius
        if radial_splits < 1:
            raise ValueError("radial_splits must be at least 1")
        self.radial_splits = radial_splits
        self.sector_config = SectorConfig(sectors=sectors, angles=angles, windmill=windmill)
        self.windmill = windmill

    def __repr__(self):
        return (f"VanishedRodDiscretizationConfig(base_radius={self.base_radius}, "
                f"radial_splits={self.radial_splits}, sector_config={self.sector_config.__repr__()})")

    def to_dict(self):
        """YAML ``vanished_rods`` entry of a ``sectorization`` block."""
        d = {"base_radius": self.base_radius, "radial_splits": self.radial_splits,
             "sectors": list(self.sector_config.sectors)}
        if self.sector_config.angles:
            d["angles"] = list(self.sector_config.angles)
        return _plain(d)

    def resolve_radii_and_sectors(self, radius=None):
        """
        Return the resolved list of radial splits for the vanished rod region, auto-computing from base_radius if needed.
        If base_radius is provided (set from the calculation scheme yaml), radial splits are computed as equally spaced radii between 0. and base_radius.
        If base_radius is None, and radius is provided, radial splits are computed as equally spaced radii between 0. and radius.
        If both base_radius and radius are None, radial splits are set to 0. This is 
        """

        if radius is not None:
            # override the user input radius to ensure consistency with the provided radius
            self.base_radius = radius
        resolved_radius = self.base_radius
        if resolved_radius is None:
            self.radial_split_points = []
        else:
            self.radial_split_points = [round(resolved_radius * (i + 1) / self.radial_splits, 6) for i in range(self.radial_splits)]
        # extend the sector config to match the number of radial splits if needed
        resolved_radial_splits = len(self.radial_split_points)
        total_number_of_regions = resolved_radial_splits + 1  # radial splits create n + 1 regions
        if self.sector_config and self.sector_config.sectors:
            sectors = [self.sector_config.sectors[0]] * resolved_radial_splits
            if self.sector_config.angles:
                angles = [self.sector_config.angles[0]] * resolved_radial_splits
                angles.extend([self.sector_config.angles[-1]] * (total_number_of_regions - len(angles))) # radial splits are within r<radius, last entry is the outermost region with r>radius
            else:
                angles = None
            # same logic for sectors
            sectors.extend([self.sector_config.sectors[-1]] * (total_number_of_regions - len(sectors)))
            angle_count = len(self.sector_config.angles) if self.sector_config.angles else 0
            if len(self.sector_config.sectors) != total_number_of_regions or angle_count != total_number_of_regions:
                    # over ride the user input sectorization to ensure consistency with the number of radial splits
                    self.sector_config = SectorConfig(sectors=sectors, angles=angles, windmill=self.windmill)




# ---------------------------------------------------------------------------
# CalculationStep – one step inside a calculation scheme
# ---------------------------------------------------------------------------

class CalculationStep:
    """
    Describes a single step in a DRAGON calculation scheme.

    A step combines:

    * **Step type** (``"self_shielding"`` or ``"flux"``).
    * **Self-shielding module** (only for ``step_type="self_shielding"``):
      ``"USS"`` or ``"SHI"`` DRAGON5 module to be used at the self-shielding step.
    * **Self-shielding method** (only for ``step_type="self_shielding"``):
      ``"RSE"`` (Resoncance Spectrum Expansion) or ``"PT"`` (CALENDF Probability Tables).
    * **Spatial solution method** (``"CP"``, ``"IC"``, ``"MOC"``).
      ``"MOC"`` is only valid for flux steps.
    * **Tracking option** (``"TISO"`` or ``"TSPC"``).  For the ``IC`` spatial
      method only ``"TISO"`` is valid; ``CP`` and ``MOC`` support both.
    * **Radial discretization** (self-shielding radii strategy) applied to
      fuel pins before this step.
    * **Azimuthal sectorization** config per pin category.

    For flux calculations the step can optionally carry a ``flux_level``
    integer (1, 2, …) so that multi-level flux schemes (e.g. coarse then
    fine) can be represented as separate steps.

    Parameters
    ----------
    name : str
        Human-readable step identifier (e.g. ``"SSH"``, ``"FLUX_L1"``).
    step_type : str
        ``"self_shielding"`` or ``"flux"``.
    self_shielding_module : str or None
        DRAGON self-shielding module: ``"USS"`` or ``"SHI"``.
        Only applicable for ``step_type="self_shielding"``, ignored otherwise.
    self_shielding_method : str or None
        Energy self-shielding method: ``"RSE"`` or ``"PT"``.
        Only applicable for ``step_type="self_shielding"``, ignored otherwise.
    fuel_self_shielded_isotopes : list[str] or None
        Optional list of fuel isotopes to be self-shielded (e.g. ``["U235", "Gd155"]``).
        When ``None`` (default), default fuel isotopes are self-shielded (see LIB).
    clad_self_shielded_isotopes : list[str] or None
        Optional list of cladding isotopes to be self-shielded (e.g. ``["Zr90"]``).
        When ``None`` (default), default cladding isotopes are self-shielded (see LIB).
    spatial_method : str
        Spatial solution method: ``"CP"``, ``"IC"``, or ``"MOC"``.
        ``"MOC"`` is only valid for flux steps.
    tracking : str
        ``"TISO"`` or ``"TSPC"``.  Automatically validated against the
        spatial method (``IC`` → ``TISO`` only).
    flux_level : int or None
        For flux steps in a multi-level scheme, the level number (1, 2, …).
        ``None`` for single-level flux or self-shielding steps.
    radial_scheme : str
        Radial discretization strategy: ``"Santamarina"``, ``"automatic"``,
        or ``"user_defined"``.
    radial_params : dict or None
        Extra parameters for the radial scheme.  E.g.
        ``{"num_radial_zones": 4}`` for ``"automatic"``.
    radial_overrides : dict or None
        Per-rod-type overrides.  Keys are rod-type identifiers
        (``"Gd"``, ``"fuel"``, or specific rod IDs); values are dicts with
        ``"scheme"`` and optional ``"params"`` keys.
    sectorization_enabled : bool
        Whether azimuthal sectorization is applied for this step.
    fuel_sectors : SectorConfig or None
        Sectorization config for standard fuel pins.
    gd_sectors : SectorConfig or None
        Sectorization config for gadolinium-bearing fuel pins.
    water_rod_sectors : SectorConfig or None
        Sectorization config for water rods.
    vanished_rods_sectors : VanishedRodDiscretizationConfig or None
        Optional sub-meshing config for the vanished rod regions.
    export_macros : bool
        Whether to export MACRO properties in the geometry.
    box_discretization : BoxDiscretizationConfig or None
        Configuration for sub-meshing the assembly-box peripheral
        regions into a grid of sub-faces (for MOC tracking).  When
        ``None``, no box discretization is applied.
    macro_assignment : 
        List of list to assign MACRO properties to the glow model.
        n_rows by n_cols table ordered in x-increasing, y-increasing logic with GEOM.yaml::lattice_description dimensions 
        Entries should be strings that will be used as glow's PropertyType.MACRO property.
        If None : by default each fuel cell, water rod and strip rectangle gets an independent macro assigned.
    leakage_model : str or None
        Flux steps only.  Model used by ``FLU:`` to compute the leakage
        coefficients in a fundamental mode (``TYPE K``) calculation with
        conservative boundary conditions.  Supported: ``"B1"``.
        When ``None`` (default) the flux calculation is a plain ``TYPE K``
        (infinite lattice, no leakage).
    leakage_correction : str or None
        Flux steps only.  How the leakage correction is applied
        (supported: ``"SIGS"``).  Requires ``leakage_model``.  When ``None``
        with a ``leakage_model`` set, DRAGON applies ``SIGS`` by default.
    """

    VALID_POLAR_QUADRATURES = ("GAUS", "CACA", "CACB", "LCMD", "OPP1", "OGAU")
    VALID_TRANSPORT_CORRECTIONS = ("NONE", "APOL") # add more as needed

    def __init__(
        self,
        name,
        step_type,
        spatial_method,
        self_shielding_module=None,
        self_shielding_method=None,
        fuel_self_shielded_isotopes=None,
        clad_self_shielded_isotopes=None,
        tracking="TISO",
        flux_level=None,
        radial_scheme="Santamarina",
        radial_params=None,
        radial_overrides=None,
        sectorization_enabled=False,
        fuel_sectors=None,
        gd_sectors=None,
        water_rod_sectors=None,
        vanished_rods_sectors=None,
        export_macros=False,
        box_discretization=None,
        num_angles_2d=8,
        line_density=25.0,
        anisotropy_level=1,
        batch_size=None,
        transport_correction=None,
        polar_angles_quadrature=None,
        number_of_polar_angles=None,
        mix_numbering_strategy="by_material",
        tdt_file_id=None,
        macro_assignment=None,
        leakage_model=None,
        leakage_correction=None,
    ):
        # --- Validate step type ---
        if step_type not in VALID_STEP_TYPES:
            raise ValueError(
                f"Invalid step_type '{step_type}'. "
                f"Valid options: {VALID_STEP_TYPES}"
            )

        # --- Validate spatial method ---
        if spatial_method not in VALID_SPATIAL_METHODS:
            raise ValueError(
                f"Invalid spatial_method '{spatial_method}'. "
                f"Valid options: {VALID_SPATIAL_METHODS}"
            )
        if spatial_method == "MOC" and step_type != "flux":
            raise ValueError(
                "MOC spatial method is only available for flux steps, "
                "not for self-shielding steps."
            )

        # --- Validate self-shielding module/method (only for self_shielding steps) ---
        if step_type == "self_shielding":
            if self_shielding_module is None:
                raise ValueError(
                    "self_shielding_module is required for self_shielding steps. "
                    f"Valid options: {VALID_SELF_SHIELDING_MODULES}"
                )
            if self_shielding_module not in VALID_SELF_SHIELDING_MODULES:
                raise ValueError(
                    f"Invalid self_shielding_module '{self_shielding_module}'. "
                    f"Valid options: {VALID_SELF_SHIELDING_MODULES}"
                )
            if self_shielding_method is None:
                raise ValueError(
                    "self_shielding_method is required for self_shielding steps. "
                    f"Valid options: {VALID_SELF_SHIELDING_METHODS}"
                )
            if self_shielding_method not in VALID_SELF_SHIELDING_METHODS:
                raise ValueError(
                    f"Invalid self_shielding_method '{self_shielding_method}'. "
                    f"Valid options: {VALID_SELF_SHIELDING_METHODS}"
                )
            if transport_correction is not None and transport_correction not in self.VALID_TRANSPORT_CORRECTIONS:
                raise ValueError(
                    f"Invalid transport_correction '{transport_correction}'. "
                    f"Valid options: {self.VALID_TRANSPORT_CORRECTIONS}"
                )
        else:
            # For non-self-shielding steps, these should be None
            if self_shielding_module is not None:
                raise ValueError(
                    "self_shielding_module should not be set for flux steps."
                )
            if self_shielding_method is not None:
                raise ValueError(
                    "self_shielding_method should not be set for flux steps."
                )

        # --- Validate tracking vs spatial method ---
        if tracking not in VALID_TRACKING_OPTIONS:
            raise ValueError(
                f"Invalid tracking '{tracking}'. "
                f"Valid options: {VALID_TRACKING_OPTIONS}"
            )
        if spatial_method == "IC" and tracking != "TISO":
            raise ValueError(
                f"Spatial method 'IC' only supports 'TISO' tracking, "
                f"but '{tracking}' was requested."
            )

        # --- Validate radial scheme ---
        if radial_scheme not in VALID_RADIAL_SCHEMES:
            raise ValueError(
                f"Invalid radial_scheme '{radial_scheme}'. "
                f"Valid options: {VALID_RADIAL_SCHEMES}"
            )

        # --- Validate mix numbering strategy ---
        if mix_numbering_strategy not in VALID_MIX_NUMBERING_STRATEGIES:
            raise ValueError(
                f"Invalid mix_numbering_strategy '{mix_numbering_strategy}'. "
                f"Valid options: {VALID_MIX_NUMBERING_STRATEGIES}"
            )

        # --- Validate tracking parameters ---
        if spatial_method == "MOC":
            if polar_angles_quadrature is None:
                raise ValueError(
                    "polar_angles_quadrature is required for MOC steps. "
                    f"Valid options: {self.VALID_POLAR_QUADRATURES}"
                )
            if polar_angles_quadrature not in self.VALID_POLAR_QUADRATURES:
                raise ValueError(
                    f"Invalid polar_angles_quadrature '{polar_angles_quadrature}'. "
                    f"Valid options: {self.VALID_POLAR_QUADRATURES}"
                )
            if number_of_polar_angles is None:
                raise ValueError(
                    "number_of_polar_angles is required for MOC steps."
                )

        # --- Validate FLU: leakage options (flux steps only) ---
        if step_type != "flux" and (
            leakage_model is not None or leakage_correction is not None
        ):
            raise ValueError(
                f"leakage_model / leakage_correction are only valid for "
                f"flux steps, but step '{name}' has step_type '{step_type}'."
            )
        if leakage_model is not None and leakage_model not in VALID_LEAKAGE_MODELS:
            raise ValueError(
                f"Invalid leakage_model '{leakage_model}' in step '{name}'. "
                f"Valid options: {VALID_LEAKAGE_MODELS}"
            )
        if leakage_correction is not None:
            if leakage_model is None:
                raise ValueError(
                    f"leakage_correction '{leakage_correction}' in step "
                    f"'{name}' requires a leakage_model (the FLU: descleak "
                    f"structure needs the leakage coefficients model first). "
                    f"Valid leakage_model options: {VALID_LEAKAGE_MODELS}"
                )
            if leakage_correction not in VALID_LEAKAGE_CORRECTIONS:
                raise ValueError(
                    f"Invalid leakage_correction '{leakage_correction}' in "
                    f"step '{name}'. "
                    f"Valid options: {VALID_LEAKAGE_CORRECTIONS}"
                )
        elif leakage_model is not None:
            log.info(
                f"FLU step '{name}': leakage_model {leakage_model} set "
                f"without leakage_correction; DRAGON will apply SIGS by "
                f"default (IGE335 §3.7)."
            )

        self.name = name
        self.step_type = step_type
        self.self_shielding_module = self_shielding_module
        self.self_shielding_method = self_shielding_method
        self.fuel_self_shielded_isotopes = fuel_self_shielded_isotopes
        self.clad_self_shielded_isotopes = clad_self_shielded_isotopes
        self.spatial_method = spatial_method
        self.tracking = tracking
        self.flux_level = flux_level
        self.radial_scheme = radial_scheme
        self.radial_params = radial_params or {}
        self.radial_overrides = radial_overrides or {}
        self.sectorization_enabled = sectorization_enabled
        self.fuel_sectors = fuel_sectors
        self.gd_sectors = gd_sectors
        self.water_rod_sectors = water_rod_sectors
        self.vanished_rods_sectors = vanished_rods_sectors
        self.export_macros = export_macros
        self.box_discretization = box_discretization
        self.num_angles_2d = num_angles_2d
        self.line_density = line_density
        self.anisotropy_level = anisotropy_level
        self.batch_size = batch_size
        self.transport_correction = transport_correction
        self.polar_angles_quadrature = polar_angles_quadrature
        self.number_of_polar_angles = number_of_polar_angles
        self.mix_numbering_strategy = mix_numbering_strategy
        self.tdt_file_id = tdt_file_id
        self.macro_assignment = macro_assignment
        self.leakage_model = leakage_model
        self.leakage_correction = leakage_correction

    # ------------------------------------------------------------------
    # Radii application
    # ------------------------------------------------------------------

    def apply_radii(self, assembly_model):
        """
        Update ``FuelPinModel.radii`` for every fuel pin in the assembly
        lattice according to this step's radial discretization config.

        Must be called before ``number_fuel_material_mixtures_by_*`` so
        that the zone count is correct.

        Parameters
        ----------
        assembly_model : CartesianAssemblyModel
            Assembly model whose lattice has been built
            (``analyze_lattice_description(build_pins=True)``).
        """
        from .DragonModel import FuelPinModel

        if assembly_model.lattice is None:
            raise RuntimeError(
                "Assembly lattice has not been built yet. "
                "Call analyze_lattice_description(build_pins=True) first."
            )

        for row in assembly_model.lattice:
            for pin in row:
                if not isinstance(pin, FuelPinModel):
                    continue

                # Determine which radial scheme to use for this pin
                scheme, params = self._resolve_radial_config(pin)

                if scheme == "Santamarina":
                    pin.subdivide_into_Santamarina_radii()
                elif scheme == "Santamarina_fine_Gd":
                    pin.subdivide_into_fine_Gd_radii()
                elif scheme == "automatic":
                    n_zones = params.get("num_radial_zones", 1)
                    pin.subdivide_into_radial_zones(n_zones)
                elif scheme == "user_defined":
                    user_radii = params.get("user_defined_radii", None)
                    if user_radii is None:
                        raise ValueError(
                            f"'user_defined' radial scheme for pin "
                            f"'{pin.fuel_material_name}' requires "
                            f"'user_defined_radii' in params."
                        )
                    pin.subdivide_into_user_defined_radii(user_radii)

                pin.self_shielding_option = scheme

    def _resolve_radial_config(self, pin):
        """
        Return ``(scheme, params)`` for a given pin, checking overrides first.

        Override keys checked (in priority order):
        1. pin.rod_ID (e.g. ``"ROD5G"``)
        2. ``"Gd"`` if pin.isGd is True
        3. ``"fuel"`` (generic fuel override)
        4. Fall back to step-level defaults.
        """
        # Check specific rod_ID override
        rod_id = getattr(pin, "rod_ID", None)
        if rod_id and rod_id in self.radial_overrides:
            ovr = self.radial_overrides[rod_id]
            return ovr.get("scheme", self.radial_scheme), ovr.get("params", self.radial_params)

        # Check Gd category override
        if pin.isGd and "Gd" in self.radial_overrides:
            ovr = self.radial_overrides["Gd"]
            return ovr.get("scheme", self.radial_scheme), ovr.get("params", self.radial_params)

        # Check generic fuel override
        if "fuel" in self.radial_overrides:
            ovr = self.radial_overrides["fuel"]
            return ovr.get("scheme", self.radial_scheme), ovr.get("params", self.radial_params)

        # Step-level defaults
        return self.radial_scheme, self.radial_params

    # ------------------------------------------------------------------
    # Sectorization query
    # ------------------------------------------------------------------

    def get_sectorization_for_pin(self, isGd=False):
        """
        Return the ``SectorConfig`` applicable to a given pin.

        Parameters
        ----------
        isGd : bool
            Whether the pin is a gadolinium-bearing pin or not.

        Returns
        -------
        SectorConfig or None
            ``None`` if sectorization is disabled for this step.
        """
        if not self.sectorization_enabled:
            return None

        if isGd and self.gd_sectors is not None:
            return self.gd_sectors
        if self.fuel_sectors is not None:
            return self.fuel_sectors
        return None

    def get_water_rod_sectorization(self):
        """
        Return the ``SectorConfig`` for water rods, or ``None`` if
        sectorization is disabled.
        """
        if not self.sectorization_enabled:
            return None
        return self.water_rod_sectors

    def get_vanished_rod_sectorization(self):
        """
        Return the ``VanishedRodDiscretizationConfig`` for vanished rods, or ``None`` if
        sectorization is disabled.
        """
        if not self.sectorization_enabled:
            return None
        return self.vanished_rods_sectors

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self):
        level_str = f", flux_level={self.flux_level}" if self.flux_level else ""
        if self.step_type == "self_shielding":
            method_str = f"module={self.self_shielding_module}, method={self.self_shielding_method}, spatial={self.spatial_method}"
        else:
            method_str = f"spatial={self.spatial_method}"
            if self.leakage_model:
                method_str += f", leakage={self.leakage_model}"
                if self.leakage_correction:
                    method_str += f" {self.leakage_correction}"
        return (
            f"CalculationStep(name='{self.name}', "
            f"type='{self.step_type}', "
            f"{method_str}, "
            f"tracking='{self.tracking}', "
            f"mix_numbering='{self.mix_numbering_strategy}'{level_str})"
        )

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self):
        """YAML step entry; ``DragonCalculationScheme._parse_step`` of it
        rebuilds an identical step (every value is written explicitly, so
        the parser defaults never apply)."""
        d = {
            "name": self.name,
            "step_type": self.step_type,
            "spatial_method": self.spatial_method,
            "tracking": self.tracking,
            "mix_numbering_strategy": self.mix_numbering_strategy,
            "tdt_file_id": self.tdt_file_id,
            "flux_level": self.flux_level,
            "num_angles_2d": self.num_angles_2d,
            "line_density": self.line_density,
            "anisotropy_level": self.anisotropy_level,
            "polar_angles_quadrature": self.polar_angles_quadrature,
            "number_of_polar_angles": self.number_of_polar_angles,
            "batch_size": self.batch_size,
            "radial_scheme": self.radial_scheme,
            "radial_params": self.radial_params,
            "radial_overrides": self.radial_overrides,
            "export_macros": self.export_macros,
        }
        if self.step_type == "self_shielding":
            d["self_shielding_module"] = self.self_shielding_module
            d["self_shielding_method"] = self.self_shielding_method
            d["transport_correction"] = self.transport_correction
            d["fuel_self_shielded_isotopes"] = self.fuel_self_shielded_isotopes
            d["clad_self_shielded_isotopes"] = self.clad_self_shielded_isotopes
        else:
            d["leakage_model"] = self.leakage_model
            d["leakage_correction"] = self.leakage_correction
        sect = {"enabled": bool(self.sectorization_enabled)}
        if self.sectorization_enabled:
            cfgs = {"fuel_pins": self.fuel_sectors, "Gd_pins": self.gd_sectors,
                    "water_rods": self.water_rod_sectors,
                    "vanished_rods": self.vanished_rods_sectors}
            windmills = [c.windmill for c in cfgs.values() if c is not None]
            sect["windmill"] = bool(windmills[0]) if windmills else False
            for key, cfg in cfgs.items():
                if cfg is not None:
                    sect[key] = cfg.to_dict()
        d["sectorization"] = sect
        if self.box_discretization is not None:
            d["box_discretization"] = self.box_discretization.to_dict()
        if self.macro_assignment is not None:
            d["macro_assignment"] = self.macro_assignment
        return _plain(d)


# ---------------------------------------------------------------------------
# CalculationBranch – parameter variation for branch calculations
# ---------------------------------------------------------------------------

VALID_BRANCH_TYPES = (
    "coolant_density", "fuel_temperature",
    "coolant_temperature", "burnup",
)

# Default ordering of branch loops (outermost → innermost)
DEFAULT_BRANCH_LOOP_ORDER = [
    "fuel_temperature", "coolant_temperature", "coolant_density",
]

# Map branch type → COMPO PARA name (CLE2000 ≤ 12 chars)
BRANCH_TYPE_TO_PARA_NAME = {
    "fuel_temperature":    "TFuel",
    "coolant_temperature": "TCool",
    "coolant_density":     "DCool",
    "burnup":              "Burnup",
}

# Map branch type → COMPO PARA keyword
BRANCH_TYPE_TO_PARA_KEYWORD = {
    "fuel_temperature":    "VALU REAL",
    "coolant_temperature": "VALU REAL",
    "coolant_density":     "VALU REAL",
    "burnup":              "IRRA",
}

# Map branch type → CLE2000 UTL: array name (≤ 12 chars)
BRANCH_TYPE_TO_ARRAY_NAME = {
    "fuel_temperature":    "fuel_temp",
    "coolant_temperature": "cool_temp",
    "coolant_density":     "cool_dens",
    "burnup":              "burnup_pts",
}

# Map branch type → CLE2000 counter variable (≤ 12 chars)
BRANCH_TYPE_TO_COUNTER = {
    "fuel_temperature":    "i_Tfuel",
    "coolant_temperature": "i_Tcool",
    "coolant_density":     "i_DCool",
    "burnup":              "i_BU",
}

# Map branch type → CLE2000 count variable (≤ 12 chars)
BRANCH_TYPE_TO_COUNT_VAR = {
    "fuel_temperature":    "n_Tfuel",
    "coolant_temperature": "n_Tcool",
    "coolant_density":     "n_DCool",
    "burnup":              "n_BU",
}


class CalculationBranch:
    """
    A parameter variation axis for branch calculations.

    Represents one dimension in the multi-parameter space over which
    DRAGON calculations are tabulated (e.g. coolant density at 5 void
    fractions, fuel temperature at 3 values, etc.).

    Attributes
    ----------
    name : str
        Human-readable label, e.g. ``"Coolant density variation"``.
    type : str
        One of ``"coolant_density"``, ``"fuel_temperature"``,
        ``"coolant_temperature"``, ``"burnup"``.
    values : list[float]
        Parameter values to iterate over (stored as floats: they are
        written to CLE2000 REAL arrays).
    labels : list[str] or None
        Optional statepoint labels, one per value (e.g. ``["00", "40",
        "80"]``), used in folder names, HDF5 groups and plots.
    compositions : list[dict] or None
        ``coolant_density`` only. One coolant composition per value,
        ``{ZAID or isotope name: N [atoms/(barn·cm)]}``, written verbatim to
        the MIX procedure (e.g. the exact benchmark isotopics). Without it the
        H1/O16 densities are computed from the mass density (pure H2O).
    """

    def __init__(self, name, branch_type, values, labels=None, compositions=None):
        if branch_type not in VALID_BRANCH_TYPES:
            raise ValueError(
                f"Invalid branch type '{branch_type}'. "
                f"Valid types: {VALID_BRANCH_TYPES}"
            )
        if branch_type == "burnup" and len(values) > 1:
            raise NotImplementedError(
                "Depletion branches (burnup with multiple values) "
                "are not yet supported. Use a single burnup value "
                "(e.g. [0.0]) for fresh-fuel calculations."
            )
        self.name = name
        self.type = branch_type
        self.values = [float(v) for v in values]
        if labels is not None:
            labels = [str(l) for l in labels]
            if len(labels) != len(self.values):
                raise ValueError(
                    f"Branch '{name}': {len(labels)} labels for "
                    f"{len(self.values)} values.")
            if len(set(labels)) != len(labels):
                raise ValueError(f"Branch '{name}': duplicate labels {labels}.")
        self.labels = labels
        if compositions is not None:
            if branch_type != "coolant_density":
                raise ValueError(
                    f"Branch '{name}': compositions are only supported for "
                    f"coolant_density branches, not '{branch_type}'.")
            if len(compositions) != len(self.values):
                raise ValueError(
                    f"Branch '{name}': {len(compositions)} compositions for "
                    f"{len(self.values)} values.")
            compositions = [normalise_composition(c) for c in compositions]
            for name_iso in {i for c in compositions for i in c}:
                coolant_density_var(name_iso)          # CLE2000 name length
        self.compositions = compositions

    def coolant_compositions(self):
        """Coolant composition at each value, ``[{isotope: N}, ...]``: the
        explicit :attr:`compositions`, or H1/O16 computed from the mass
        density (``compute_water_iso_densities_at_densities``, pure H2O).
        Isotopes absent from one statepoint get 0."""
        if self.type != "coolant_density":
            raise ValueError(f"Branch '{self.name}' is not a coolant_density branch.")
        if self.compositions is None:
            from ..MaterialProperties.water_properties import (
                compute_water_iso_densities_at_densities,
            )
            comps = compute_water_iso_densities_at_densities(self.values)
        else:
            comps = self.compositions
        isotopes = self.composition_isotopes()
        return [{iso: float(c.get(iso, 0.0)) for iso in isotopes} for c in comps]

    def composition_isotopes(self):
        """Ordered isotopes whose coolant density varies with this branch."""
        if self.compositions is None:
            return ["H1", "O16"]
        seen = []
        for c in self.compositions:
            seen += [i for i in c if i not in seen]
        return seen

    def label_of(self, value):
        """Label of the statepoint *value* (``None`` without labels)."""
        i = self.index_of(value)
        return None if (i is None or self.labels is None) else self.labels[i]

    def index_of(self, value, rtol=STATEPOINT_RTOL):
        """Index of *value* in :attr:`values` (within *rtol*), or ``None``."""
        for i, v in enumerate(self.values):
            if values_match(v, value, rtol):
                return i
        return None

    def to_dict(self):
        """YAML ``CALCULATION_BRANCHES`` entry."""
        d = {"name": self.name, "type": self.type, "values": list(self.values)}
        if self.labels is not None:
            d["labels"] = list(self.labels)
        if self.compositions is not None:
            d["compositions"] = [dict(c) for c in self.compositions]
        return _plain(d)

    @property
    def para_name(self):
        """CLE2000 COMPO PARA name for this branch type."""
        return BRANCH_TYPE_TO_PARA_NAME[self.type]

    @property
    def para_keyword(self):
        """CLE2000 COMPO PARA keyword for this branch type."""
        return BRANCH_TYPE_TO_PARA_KEYWORD[self.type]

    @property
    def array_name(self):
        """CLE2000 UTL: CREA array name."""
        return BRANCH_TYPE_TO_ARRAY_NAME[self.type]

    @property
    def counter_var(self):
        """CLE2000 WHILE loop counter variable name."""
        return BRANCH_TYPE_TO_COUNTER[self.type]

    @property
    def count_var(self):
        """CLE2000 INTEGER variable holding the number of points."""
        return BRANCH_TYPE_TO_COUNT_VAR[self.type]

    def __repr__(self):
        return (
            f"CalculationBranch(name='{self.name}', "
            f"type='{self.type}', "
            f"values={self.values})"
        )

    @classmethod
    def from_yaml_entry(cls, entry):
        """Parse a single branch entry from YAML.

        Parameters
        ----------
        entry : dict
            ``{"name": ..., "type": ..., "values": [...]}``.

        Returns
        -------
        CalculationBranch
        """
        return cls(
            name=entry["name"],
            branch_type=entry["type"],
            values=entry["values"],
            labels=entry.get("labels"),
            compositions=entry.get("compositions"),
        )


# ---------------------------------------------------------------------------
# CalculationOutput – specification for an EDI/COMPO output
# ---------------------------------------------------------------------------

class CalculationOutput:
    """
    Specification for a single EDI → COMPO output volume.

    Each output defines the spatial integration, energy condensation,
    isotope list, and the subset of state points at which to generate
    the edition.

    Attributes
    ----------
    name : str
        COMPO directory / EDI ``SAVE ON`` name (≤ 12 chars for CLE2000).
    comment : str
        Human-readable description.
    isotopes : list[str]
        Isotopes tracked (``MICR`` / ``ISOT`` keywords).
    spatial_integration_mode : str
        ``"FUEL"``, ``"ALL"``, ``"by_pin"``, etc.
    energy_bounds : list[float] | None | str
        ``None`` or ``"None"`` → no COND (keep all groups).
        ``[]`` → COND (collapse to 1 group).
        ``[0.625]`` → COND 0.625 (2-group).
    state_points : str | dict
        ``"ALL"`` — tabulate at every branch combination.
        ``dict`` — selective: ``{branch_type: [value, ...]}`` listing
        the specific values at which this output is generated.
    """

    def __init__(self, name, comment, isotopes,
                 spatial_integration_mode, energy_bounds,
                 state_points):
        self.name = name
        self.comment = comment
        self.isotopes = list(isotopes)
        self.spatial_integration_mode = spatial_integration_mode

        # Normalise energy_bounds
        if energy_bounds is None or energy_bounds == "None":
            self.energy_bounds = None
        else:
            self.energy_bounds = list(energy_bounds)

        # Normalise state_points
        if isinstance(state_points, str) and state_points.upper() == "ALL":
            self.state_points = "ALL"
        elif isinstance(state_points, dict):
            self.state_points = {
                k: [float(x) for x in v] for k, v in state_points.items()
            }
        else:
            self.state_points = "ALL"

    def applies_to(self, statepoint_dict):
        """Check whether this output should be generated at a given statepoint.

        Parameters
        ----------
        statepoint_dict : dict
            ``{branch_type: current_value}`` for each active branch,
            e.g. ``{"coolant_density": 0.736, "fuel_temperature": 900.0, ...}``.

        Returns
        -------
        bool
        """
        if self.state_points == "ALL":
            return True
        for branch_type, required_values in self.state_points.items():
            current = statepoint_dict.get(branch_type)
            if current is None:
                continue
            if not any(values_match(current, v) for v in required_values):
                return False
        return True

    def to_dict(self):
        """YAML ``CALCULATION_OUTPUTS`` entry."""
        return _plain({
            "name": self.name,
            "comment": self.comment,
            "isotopes": list(self.isotopes),
            "spatial_integration_mode": self.spatial_integration_mode,
            "energy_bounds": self.energy_bounds,
            "state_points": self.state_points,
        })

    def __repr__(self):
        return (
            f"CalculationOutput(name='{self.name}', "
            f"mode='{self.spatial_integration_mode}', "
            f"bounds={self.energy_bounds}, "
            f"state_points={self.state_points})"
        )

    @classmethod
    def from_yaml_entry(cls, entry):
        """Parse a single output entry from YAML.

        Parameters
        ----------
        entry : dict

        Returns
        -------
        CalculationOutput
        """
        return cls(
            name=entry["name"],
            comment=entry.get("comment", ""),
            isotopes=entry.get("isotopes", []),
            spatial_integration_mode=entry.get(
                "spatial_integration_mode",
                entry.get("spatial_intergation_mode", "FUEL"),
            ),
            energy_bounds=entry.get("energy_bounds", None),
            state_points=entry.get("state_points", "ALL"),
        )


# ---------------------------------------------------------------------------
# EditionBetweenLevelsStep – energy condensation step between flux levels
# ---------------------------------------------------------------------------

class EditionBetweenLevelsStep:
    """
    Describes the energy condensation step between two flux levels in a
    multi-level DRAGON calculation scheme.

    This step does NOT involve geometry tracking of its own; it operates
    on the flux and library structures produced by the preceding flux step.
    It therefore does not inherit from ``CalculationStep`` and has no
    spatial_method, tracking, or sectorization attributes.

    The step drives:

    * An ``EDI:`` module call to condense cross sections from the fine
      group mesh to ``number_of_macro_groups`` coarse groups using the
      boundaries listed in ``energy_groups_bounds``.
    * An optional ``SPH:`` equivalence correction applied to the condensed
      library up to group ``max_sph_group``.

    Parameters
    ----------
    name : str
        Human-readable step identifier.
    number_of_macro_groups : int
        Number of coarse energy groups to condense to (e.g. 26).
    energy_groups_bounds : list[int]
        Upper-group indices of each coarse energy group boundary used in
        the ``COND`` keyword of ``EDI:``.  Typically ``number_of_macro_groups - 1``
        values (the last group boundary is implicit).
    sph_correction : bool
        Whether to apply an SPH equivalence correction after condensation.
        Default ``False``.
    max_sph_group : int or None
        Maximum coarse group index up to which SPH is applied
        (``GRMAX`` keyword).  Required when ``sph_correction=True``.
    max_iterations : int or None
        Maximum number of fixed point SPH iterations to be performed.
        ("ITER" keyword). By default, SPH: will iterate until sph factors are converged to tolerance 
        or reaches the maximal value set in DRAGON (200).
    tolerance : float 
        Tolerance to test SPH factors' convergence, by default this takes the value of 1e-4.
    """

    step_type = "edition_between_levels"

    def __init__(self, name, number_of_macro_groups,
                 energy_groups_bounds, sph_correction=False,
                 max_sph_group=None, max_iterations = 200, tolerance = 1e-4):
        self.name = name
        self.number_of_macro_groups = number_of_macro_groups
        self.energy_groups_bounds = list(energy_groups_bounds)
        self.sph_correction = bool(sph_correction)
        self.max_sph_group = max_sph_group
        self.max_iterations = max_iterations
        self.tolerance = tolerance

        if self.sph_correction and self.max_sph_group is None:
            raise ValueError(
                f"EditionBetweenLevelsStep '{name}': "
                "max_sph_group is required when sph_correction=True."
            )

    def __repr__(self):
        return (
            f"EditionBetweenLevelsStep("
            f"name='{self.name}', "
            f"n_groups={self.number_of_macro_groups}, "
            f"sph={self.sph_correction}, "
            f"max_iterations={self.max_iterations}, "
            f"tolerance={self.tolerance})"
        )

    def to_dict(self):
        """YAML step entry (inverse of ``DragonCalculationScheme._parse_step``)."""
        return _plain({
            "name": self.name,
            "step_type": self.step_type,
            "number_of_macro_groups": self.number_of_macro_groups,
            "energy_groups_bounds": self.energy_groups_bounds,
            "SPH_correction": self.sph_correction,
            "max_SPH_group": self.max_sph_group,
            "max_iterations": self.max_iterations,
            "tolerance": self.tolerance,
        })


# ---------------------------------------------------------------------------
# DragonCalculationScheme – ordered collection of CalculationSteps
# ---------------------------------------------------------------------------

class DragonCalculationScheme:
    """
    Orchestrates the sequence of DRAGON calculation steps for an assembly.

    A typical scheme contains at minimum:

    1. A **self-shielding** step (SSH) – computes resonance integrals on a
       coarse geometry.
    2. One or more **flux** steps – solves the transport equation on a
       (possibly finer) geometry.

    Each step carries its own radial discretization and sectorization
    config, so different geometries can be used at each stage.

    Construction
    ------------
    * Programmatic: instantiate and call ``add_step()``.
    * From YAML:    ``DragonCalculationScheme.from_yaml(path)``.
    * From preset:  ``DragonCalculationScheme.preset("BWR_standard")``.

    Parameters
    ----------
    name : str
        Scheme identifier.
    """

    def __init__(self, name="default", tag=None):
        self.name = name
        self.steps = []
        self.branches = []   # list[CalculationBranch]
        self.outputs = []    # list[CalculationOutput]
        self.branch_loop_order = list(DEFAULT_BRANCH_LOOP_ORDER)
        self.tag = None
        if tag is not None:
            self.set_tag(tag)

    # ------------------------------------------------------------------
    # Step management
    # ------------------------------------------------------------------

    def add_step(self, step):
        """
        Append a ``CalculationStep`` or ``EditionBetweenLevelsStep`` to the scheme.

        Parameters
        ----------
        step : CalculationStep or EditionBetweenLevelsStep
        """
        if not isinstance(step, (CalculationStep, EditionBetweenLevelsStep)):
            raise TypeError(
                f"Expected CalculationStep or EditionBetweenLevelsStep, "
                f"got {type(step).__name__}"
            )
        self.steps.append(step)

    def get_step(self, name):
        """
        Retrieve a step by name.

        Parameters
        ----------
        name : str

        Returns
        -------
        CalculationStep

        Raises
        ------
        KeyError
            If no step with the given name exists.
        """
        for step in self.steps:
            if step.name == name:
                return step
        raise KeyError(f"No step named '{name}' in scheme '{self.name}'.")

    def get_self_shielding_steps(self):
        """Return all self-shielding steps."""
        return [s for s in self.steps if s.step_type == "self_shielding"]

    def get_flux_steps(self):
        """Return all flux steps, ordered by flux_level if set."""
        flux = [s for s in self.steps if s.step_type == "flux"]
        flux.sort(key=lambda s: (s.flux_level or 0))
        return flux

    def get_edition_between_levels_steps(self):
        """Return all edition-between-levels steps (energy condensation)."""
        return [s for s in self.steps
                if s.step_type == "edition_between_levels"]

    def is_two_level_scheme(self):
        """Return ``True`` if the scheme contains an edition-between-levels step."""
        return len(self.get_edition_between_levels_steps()) > 0

    def get_trackable_steps(self):
        """Return steps that require TDT geometry and tracking.

        Excludes ``edition_between_levels`` steps, which operate on
        existing tracking structures rather than their own geometry.
        """
        return [s for s in self.steps
                if s.step_type != "edition_between_levels"]

    # ------------------------------------------------------------------
    # Branch / output helpers
    # ------------------------------------------------------------------

    def has_branches(self):
        """Return ``True`` if calculation branches are defined."""
        return len(self.branches) > 0

    def get_branch(self, branch_type):
        """Return the ``CalculationBranch`` for *branch_type*, or ``None``."""
        for b in self.branches:
            if b.type == branch_type:
                return b
        return None

    def statepoint_values(self, statepoint=None):
        """``{branch_type: value}`` of one statepoint of the branch grid.

        *statepoint* is ``None`` (only if every branch has a single value),
        a label of one labelled branch (e.g. ``"40"``; the other branches
        must be single-valued), or ``{branch_type: value or label}`` (missing
        single-valued branches are filled in). Values are matched within
        ``STATEPOINT_RTOL``; the stored value is returned.
        """
        chosen = {}
        if isinstance(statepoint, str):
            owners = [b for b in self.branches if b.labels and statepoint in b.labels]
            if len(owners) != 1:
                raise KeyError(f"Label '{statepoint}' found in {len(owners)} branches.")
            chosen[owners[0].type] = statepoint
        elif statepoint is not None:
            chosen = dict(statepoint)
        unknown = set(chosen) - {b.type for b in self.branches}
        if unknown:
            raise KeyError(f"No branch of type {sorted(unknown)}.")
        values = {}
        for b in self.branches:
            if b.type not in chosen:
                if len(b.values) != 1:
                    raise ValueError(f"Statepoint must give a value for '{b.type}' "
                                     f"({len(b.values)} values).")
                values[b.type] = b.values[0]
                continue
            sel = chosen[b.type]
            if isinstance(sel, str):
                if not b.labels or sel not in b.labels:
                    raise KeyError(f"'{sel}' is not a label of branch '{b.type}'.")
                values[b.type] = b.values[b.labels.index(sel)]
            else:
                i = b.index_of(float(sel))
                if i is None:
                    raise ValueError(f"{sel} is not a value of branch '{b.type}': {b.values}")
                values[b.type] = b.values[i]
        return values

    def get_ordered_branches(self):
        """Return branches in the configured loop order (outermost first).

        Branches whose type is not in ``branch_loop_order`` are appended
        at the end in their original order.
        """
        by_type = {b.type: b for b in self.branches}
        ordered = []
        for btype in self.branch_loop_order:
            if btype in by_type:
                ordered.append(by_type.pop(btype))
        # Append any remaining branches not in the configured order
        for b in self.branches:
            if b.type in by_type:
                ordered.append(by_type.pop(b.type))
        return ordered

    def get_total_statepoints(self):
        """Return the total number of statepoint combinations.

        This is the product of the number of values in each branch.
        Returns 1 when no branches are defined (single statepoint).
        """
        total = 1
        for b in self.branches:
            total *= len(b.values)
        return total

    # ------------------------------------------------------------------
    # Serialisation, copies and programmatic variants (plan §3.8, P2)
    # ------------------------------------------------------------------

    def to_dict(self):
        """Plain dict in the YAML layout read by :meth:`from_dict` /
        :meth:`from_yaml`. ``from_dict(to_dict())`` rebuilds a scheme that
        generates identical CLE2000 procedures."""
        head = {"name": self.name}
        if self.tag is not None:
            head["tag"] = self.tag
        head["steps"] = [st.to_dict() for st in self.steps]
        return _plain({
            "DRAGON_CALCULATION_SCHEME": head,
            "CALCULATION_BRANCHES": [b.to_dict() for b in self.branches],
            "CALCULATION_OUTPUTS": [o.to_dict() for o in self.outputs],
            "BRANCH_LOOP_ORDER": list(self.branch_loop_order),
        })

    def to_yaml(self, path=None, header=None):
        """YAML text of :meth:`to_dict`; also written to *path* if given.

        Every parameter is written explicitly (parser defaults included), so
        the file is a complete, resolved description of the scheme.
        """
        lines = [
            "# Resolved DRAGON calculation scheme written by starterDD",
            f"# scheme: {self.name}   variant_tag: {self.variant_tag}",
        ]
        if header:
            lines += [f"# {l}" for l in str(header).splitlines()]
        text = "\n".join(lines) + "\n" + yaml.safe_dump(
            self.to_dict(), sort_keys=False, default_flow_style=None, width=100)
        if path is not None:
            with open(path, "w") as f:
                f.write(text)
        return text

    def copy(self):
        """Independent deep copy."""
        return deepcopy(self)

    #: YAML step keys accepted by :meth:`set_step` (everything
    #: :meth:`_parse_step` reads).
    STEP_KEYS = frozenset({
        "name", "step_type", "spatial_method", "tracking", "mix_numbering_strategy",
        "tdt_file_id", "flux_level", "num_angles_2d", "line_density",
        "anisotropy_level", "polar_angles_quadrature", "number_of_polar_angles",
        "batch_size", "radial_scheme", "radial_params", "radial_overrides",
        "export_macros", "self_shielding_module", "self_shielding_method",
        "transport_correction", "fuel_self_shielded_isotopes",
        "clad_self_shielded_isotopes", "leakage_model", "leakage_correction",
        "sectorization", "box_discretization", "macro_assignment",
        "number_of_macro_groups", "energy_groups_bounds", "SPH_correction",
        "max_SPH_group", "max_iterations", "tolerance",
    })

    def set_step(self, name, **params):
        """Change parameters of step *name* **in place** and return ``self``.

        *params* use the YAML step keys (e.g. ``num_angles_2d=12``,
        ``spatial_method="CP"``, ``macro_assignment=[[...]]``). Nested dicts
        such as ``sectorization`` are merged key by key, so
        ``sectorization={"enabled": False}`` keeps the sub-blocks. The step is
        rebuilt by the YAML parser, so every validation applies.

        Raises:
            KeyError: unknown step or parameter name.
        """
        unknown = set(params) - self.STEP_KEYS
        if unknown:
            raise KeyError(
                f"Unknown step parameter(s) {sorted(unknown)}. "
                f"Known: {sorted(self.STEP_KEYS)}")
        for i, st in enumerate(self.steps):
            if st.name == name:
                self.steps[i] = self._parse_step(_deep_merge(st.to_dict(), params))
                return self
        raise KeyError(f"No step named '{name}' in scheme '{self.name}'.")

    def with_step(self, name, **params):
        """Copy of the scheme with :meth:`set_step` applied."""
        return self.copy().set_step(name, **params)

    def with_branch_values(self, branch_type, values, labels=None,
                           selective_outputs="keep", name=None,
                           compositions=None):
        """Copy of the scheme with new values for one branch axis.

        Selective outputs (``state_points`` dicts) are kept **in sync**, so an
        output can no longer disappear silently because its CLE2000 guard
        never fires.

        Parameters
        ----------
        branch_type : str
            One of ``VALID_BRANCH_TYPES``. The branch is created if absent.
        values : list[float]
        labels : list[str] or None
            Statepoint labels, one per value.
        selective_outputs : ``"keep"`` | ``"all"`` | dict
            ``"keep"`` (default): an output that requested **every** old
            value of the axis gets every new value; otherwise it keeps the
            requested values that are still present (matched within
            ``STATEPOINT_RTOL``). Outputs left with none raise
            ``ValueError`` (all listed). ``"all"``: drop the restriction on
            this axis.
            ``{output_name: [values]}``: explicit new values for those
            outputs (the others follow ``"keep"``).
        name : str or None
            Name of a new branch (default: *branch_type*).
        compositions : list[dict] or None
            ``coolant_density`` only: one coolant composition per value (see
            :class:`CalculationBranch`). Compositions of the replaced branch
            are **not** carried over, since they belong to the old values.
        """
        if not isinstance(selective_outputs, dict) and selective_outputs not in ("keep", "all"):
            raise ValueError(
                f"selective_outputs must be 'keep', 'all' or a dict, "
                f"got {selective_outputs!r}.")
        new = self.copy()
        br = new.get_branch(branch_type)
        old_values = list(br.values) if br is not None else []
        if br is None:
            br = CalculationBranch(name or branch_type, branch_type, values, labels,
                                   compositions)
            new.branches.append(br)
        else:
            idx = new.branches.index(br)
            br = new.branches[idx] = CalculationBranch(br.name, branch_type, values,
                                                       labels, compositions)
        explicit = selective_outputs if isinstance(selective_outputs, dict) else {}
        unknown = set(explicit) - {o.name for o in new.outputs}
        if unknown:
            raise KeyError(f"Unknown output(s) in selective_outputs: {sorted(unknown)}")
        lost = []
        for out in new.outputs:
            if out.name in explicit:
                if not isinstance(out.state_points, dict):
                    out.state_points = {}
                out.state_points[branch_type] = [float(v) for v in explicit[out.name]]
                continue
            sp = out.state_points
            if not isinstance(sp, dict) or branch_type not in sp:
                continue
            if selective_outputs == "all":
                del sp[branch_type]
                if not sp:
                    out.state_points = "ALL"
                continue
            requested = sp[branch_type]
            if old_values and all(any(values_match(o, r) for r in requested)
                                  for o in old_values):
                sp[branch_type] = list(br.values)    # meant "every statepoint"
                continue
            kept = []
            for v in requested:
                j = br.index_of(v)
                if j is not None and br.values[j] not in kept:
                    kept.append(br.values[j])
            if not kept:
                lost.append(f"'{out.name}' ({branch_type} = {requested})")
                continue
            sp[branch_type] = kept
        if lost:
            raise ValueError(
                f"These outputs request no value of the new {branch_type} "
                f"values {br.values} and would never be generated: "
                f"{', '.join(lost)}. Pass selective_outputs={{name: [values]}} "
                f"for them, or selective_outputs='all'.")
        return new.validate()

    def with_coolant_statepoints(self, statepoints, densities=None,
                                 material="COOLANT", selective_outputs="keep"):
        """Copy of the scheme whose ``coolant_density`` branch carries the
        **exact** coolant isotopics of each statepoint (plan §3.4).

        Parameters
        ----------
        statepoints : dict
            ``{label: source}`` in statepoint order, e.g. ``{"00":
            "material_compos_00.yaml", "40": ..., "80": ...}``. A source is a
            starterDD material YAML (its *material* entry is read) or a
            composition dict ``{ZAID or isotope: N}``.
        densities : dict or None
            ``{label: rho [g/cm³]}`` used as the ``DCool`` PARA value. By
            default rho is computed from the isotopics
            (:func:`composition_mass_density`), so the PARA value and the
            isotopics are consistent by construction.
        material : str
            Material name read from the YAML sources.
        selective_outputs
            As in :meth:`with_branch_values`.

        The compositions are written verbatim to the MIX procedure, one
        ``UTL:`` array per isotope, so a single multi-statepoint COMPO uses
        the same isotopics as separate per-statepoint runs (and as Serpent).
        """
        if not statepoints:
            raise ValueError("No statepoints given.")
        labels, comps, values = [], [], []
        for label, src in statepoints.items():
            if isinstance(src, dict):
                comp = normalise_composition(src)
            else:
                comp = read_material_composition(src, material)
            labels.append(str(label))
            comps.append(comp)
            if densities is not None and label in densities:
                values.append(float(densities[label]))
            else:
                values.append(composition_mass_density(comp))
        for i in range(len(values)):
            for j in range(i):
                if values_match(values[i], values[j]):
                    raise ValueError(
                        f"Statepoints '{labels[j]}' and '{labels[i]}' have the same "
                        f"coolant density {values[i]:.7g} g/cm3 (within "
                        f"{STATEPOINT_RTOL:g}); they cannot be told apart in the COMPO.")
        return self.with_branch_values(
            "coolant_density", values, labels=labels, compositions=comps,
            selective_outputs=selective_outputs, name="Coolant density variation")

    def validate(self):
        """Check that branches and selective outputs are consistent.

        Every ``state_points`` entry of an output must name an existing
        branch and list values present in that branch (within
        ``STATEPOINT_RTOL``). Otherwise the generated CLE2000 guard never
        fires (or uses an undeclared parameter) and the output is missing
        from the COMPO. Without branches no guard is generated (single
        statepoint, every output produced), so nothing is checked. Returns
        ``self``.
        """
        if not self.branches:
            return self
        problems = []
        for out in self.outputs:
            if not isinstance(out.state_points, dict):
                continue
            for btype, vals in out.state_points.items():
                br = self.get_branch(btype)
                if br is None:
                    problems.append(
                        f"output '{out.name}': state_points on '{btype}', "
                        f"which is not a branch of the scheme")
                    continue
                if not vals:
                    problems.append(f"output '{out.name}': empty {btype} list")
                missing = [v for v in vals if br.index_of(v) is None]
                if missing:
                    problems.append(
                        f"output '{out.name}': {btype} value(s) {missing} not "
                        f"in the branch values {br.values}")
        if problems:
            raise ValueError(
                "Inconsistent calculation outputs (they would never be "
                "generated):\n  - " + "\n  - ".join(problems))
        return self

    def descriptor(self):
        """Short scheme descriptor ``<SSH_method>_<SSH_spatial>_<N>L_<flux
        methods>[_<n>sp]``, e.g. ``PT_IC_2L_IC_MOC`` (run directory name)."""
        ssh_steps = self.get_self_shielding_steps()
        flux_steps = self.get_flux_steps()
        if ssh_steps:
            ssh_method = (ssh_steps[0].self_shielding_method or "RSE").upper()
            ssh_spatial = (ssh_steps[0].spatial_method or "CP").upper()
        else:
            ssh_method = ssh_spatial = "UNK"
        flux_methods = ("_".join(st.spatial_method.upper() for st in flux_steps)
                        if flux_steps else "NONE")
        parts = [ssh_method, ssh_spatial, f"{len(flux_steps)}L", flux_methods]
        n_sp = self.get_total_statepoints()
        if n_sp > 1:
            parts.append(f"{n_sp}sp")
        return "_".join(parts)

    def content_hash(self, length=8):
        """Hash of everything that changes the generated procedures (steps,
        branches, outputs, loop order), not of the scheme name or tag."""
        import hashlib
        import json
        d = self.to_dict()
        d["DRAGON_CALCULATION_SCHEME"].pop("name", None)
        d["DRAGON_CALCULATION_SCHEME"].pop("tag", None)
        blob = json.dumps(d, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode()).hexdigest()[:length]

    @property
    def variant_tag(self):
        """Identifier of this scheme variant: the user :attr:`tag` if set,
        otherwise ``<descriptor>_<content hash>`` (e.g.
        ``PT_IC_2L_IC_MOC_3sp_1a2b3c4d``). Two schemes share an automatic tag
        only if they generate the same procedures."""
        if self.tag is not None:
            return self.tag
        return f"{self.descriptor()}_{self.content_hash()}"

    def set_tag(self, tag):
        """Set a user variant tag (letters, digits, ``_ - . +``; it is used in
        directory names). ``None`` restores the automatic tag. Returns
        ``self``."""
        import re
        if tag is not None and not re.fullmatch(r"[A-Za-z0-9_.+-]+", str(tag)):
            raise ValueError(
                f"Invalid tag {tag!r}: use letters, digits and _ - . + only.")
        self.tag = None if tag is None else str(tag)
        return self

    # ------------------------------------------------------------------
    # YAML construction
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, yaml_path):
        """
        Parse a DRAGON calculation scheme from a YAML file.

        Expected YAML structure::

            DRAGON_CALCULATION_SCHEME:
              name: "MyScheme" (human-readable identifier)
              steps:
                - name: "name_of_step" (human-readable identifier)
                  step_type: "flux" or "self_shielding"
                  mix_numbering_strategy: "by_material" or "by_pin" (optional, default "by_material")
                  self_shielding_module: "USS" or "SHI" (only for self_shielding steps)
                  self_shielding_method: "RSE" or "PT"  (only for self_shielding steps)
                  spatial_method: "IC" or "MOC" or "CP"
                  tracking: "TISO" or "TSPC" (TISO only for IC)
                  flux_level: <n>  # optional, for multi-level flux schemes
                  radial_scheme: "Santamarina"
                  export_macros: boolean flag to export MACRO properties defined in glow_builder, used in IC tracking
                  
                  sectorization:
                    enabled: boolean flag to enable azimuthal sectorization for this step
                    windmill: boolean flag to apply windmill sectorization (in outermost region only)
                    fuel_pins:
                      sectors: [(list of sector counts per region, e.g. [1, 1, 1, 1, 1, 1, 8] for UOX (4 regions)+gap+clad with subdivided coolant)]
                      angles:  [(list of angles to apply to each region sector, e.g. [0, 0, 0, 0, 0, 0, 0, 0, 22.5] for 22.5° offset on coolant only)]
                    Gd_pins:
                      sectors: [(same as fuel_pins but for Gd-bearing pins, e.g. [1, 1, 1, 1, 1, 1, 1, 1, 8] for Santamarina sectorization on all rings and subdivided coolant)]
                      angles:  [(same as fuel_pins but for Gd-bearing pins, e.g. [0, 0, 0, 0, 0, 0, 0, 0, 22.5] for 22.5° offset on coolant only)]
                    water_rods:
                      # For circular water rods (azimuthal sectorization):
                      sectors: [(list of sector counts per ring for water rods, e.g. [1, 1, 8] for 1 moderator, 1 clad, and subdivided coolant)]
                      angles:  [(list of angles for water rods, e.g. [0, 0, 22.5] for 22.5° offset on coolant only)]
                      additional_radial_splits_in_moderator: <int or list[float]>
                      # int N (default 1): N=1 no extra splits; N>=2 produces N-1 evenly spaced radii in the moderator region (r < inner_radius).
                      # list[float]: explicit user-defined radii (must be > 0 and < inner_radius).
                      # When extra rings are added, sectors[0]/angles[0] are automatically replicated for each new sub-ring.
                      subdivisions_coolant_corners: <int> : number of subdivisions to apply to the corner regions of the coolant corners in circular water rods
                      allows for the submeshing of regions with r>outer_radius, but within the water rod cell.
                      # For square water rods (Cartesian grid sub-meshing):
                      splits: [nx, ny]  # grid subdivisions applied to the whole bounding box; material is reassigned by geometric containment
                    vanished_rods:
                      base_radius: <float> (optional)  # base radius for sectorization : if not provided use the default_sectorization_radius from the vanished rod model.
                      radial_splits : int (optional) # number of radial splits used to sundivide the 0<r<base_radius or 0<r<default_sectorization_radius region of the vanished rod cell; default 1 (no subdivision))
                      sectors: [(list of sector counts per ring for vanished rods, e.g. [1, 1, 8] for 2 inner rings with 1 sector each and subdivided outer region)]
                      angles: [(list of angles for vanished rods, e.g. [0, 0, 22.5] for 22.5° offset on outer region only)]
                
                  box_discretization: 
                    enabled: boolean flag to enable sub-meshing of the assembly box peripheral regions into a grid of sub-faces (for MOC tracking)
                    corner_splits: [nx, ny]  # grid size for corner regions, default [4, 4]
                    gap_splits: [n_parallel, n_perpendicular]  # grid size for side strips, auto-permuted for H/V
                    # Deprecated aliases (still accepted):
                    # side_x_splits / side_y_splits
                    cross_moderator_discretization:  # optional, moderator regions surrounding control cross
                      narrow_gap_splits: [n_par, n_perp]  # null = auto from gap_splits density
                      moderator_at_cross_corner_splits: [nx, ny]  # null = auto
                      stub_splits: [n_par, n_perp]        # null = auto
                    control_cross_submesh: false  # or true (defaults) or dict:
                    #   control_cross_submesh:
                    #     enabled: true
                    #     control_cross_corner_splits: [n_along, n_across]  # null = (1,1)
                    #     central_structure_splits: [n_along, n_across]  # null = (1,1)
                    #     extend_splits_at_tube_boundaries: true
                    #     split_tubes_in_half: false
        
        Note : the sectorization config matches sectorization options in glow's ``RectCell.sectorize()`` method.
                      
        Parameters
        ----------
        yaml_path : str
            Path to the YAML file.

        Returns
        -------
        DragonCalculationScheme
        """
        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data):
        """Build a scheme from a dict with the :meth:`from_yaml` layout
        (e.g. the output of :meth:`to_dict`)."""
        data = deepcopy(data)
        # Accept multiple top-level key names for the flux scheme
        scheme_data = (
            data.get("DRAGON_CALCULATION_SCHEME")
            or data.get("FLUX_CALCULATION_SCHEME")
            or data
        )
        scheme_name = scheme_data.get("name", "from_yaml")
        scheme = cls(name=scheme_name, tag=scheme_data.get("tag"))

        for step_data in scheme_data.get("steps", []):
            step = cls._parse_step(step_data)
            scheme.add_step(step)

        # --- Parse CALCULATION_BRANCHES ---
        branches_data = (
            data.get("CALCULATION_BRANCHES")
            or scheme_data.get("CALCULATION_BRANCHES")
            or []
        )
        # branches_data may be at the same level as FLUX_CALCULATION_SCHEME
        # (top-level) or nested inside it (unified YAML).
        if isinstance(branches_data, list):
            for entry in branches_data:
                scheme.branches.append(
                    CalculationBranch.from_yaml_entry(entry)
                )

        # --- Parse CALCULATION_OUTPUTS ---
        outputs_data = (
            data.get("CALCULATION_OUTPUTS")
            or scheme_data.get("CALCULATION_OUTPUTS")
            or data.get("EDI_COMPO_OUTPUTS")
            or []
        )
        if isinstance(outputs_data, list):
            for entry in outputs_data:
                scheme.outputs.append(
                    CalculationOutput.from_yaml_entry(entry)
                )

        # --- Optional loop ordering ---
        loop_order = (
            data.get("BRANCH_LOOP_ORDER")
            or scheme_data.get("BRANCH_LOOP_ORDER")
        )
        if loop_order:
            scheme.branch_loop_order = list(loop_order)

        return scheme

    @staticmethod
    def _parse_step(d):
        """Parse a single step dict into a CalculationStep or EditionBetweenLevelsStep."""
        # --- Early exit for edition_between_levels (no geometry/tracking attributes) ---
        step_type = d.get("step_type")
        if step_type == "edition_between_levels":
            return EditionBetweenLevelsStep(
                name=d["name"],
                number_of_macro_groups=d["number_of_macro_groups"],
                energy_groups_bounds=d["energy_groups_bounds"],
                sph_correction=d.get("SPH_correction", False),
                max_sph_group=d.get("max_SPH_group", None),
                max_iterations=d.get("max_iterations", 200),
                tolerance=d.get("tolerance", 1e-4),
            )

        # --- Sectorization ---
        sect = d.get("sectorization", {})
        sect_enabled = sect.get("enabled", False)

        fuel_sectors = None
        gd_sectors = None
        water_rod_sectors = None
        vanished_rods_sectors = None

        if sect_enabled:
            fp = sect.get("fuel_pins", {})
            if fp:
                fuel_sectors = SectorConfig(
                    sectors=fp.get("sectors", []),
                    angles=fp.get("angles", []),
                    windmill=sect.get("windmill", False),
                )
            gp = sect.get("Gd_pins", {})
            if gp:
                gd_sectors = SectorConfig(
                    sectors=gp.get("sectors", []),
                    angles=gp.get("angles", []),
                    windmill=sect.get("windmill", False),
                )
            wr = sect.get("water_rods", {})
            if wr:
                wr_splits = wr.get("splits", None)
                wr_sectors = wr.get("sectors", [])
                wr_angles = wr.get("angles", [])
                wr_additional = wr.get(
                    "additional_radial_splits_in_moderator", None
                )
                wr_corners = wr.get("subdivisions_coolant_corners", None)

                # Warn if both circular (sectors/angles) and square
                # (splits) keys are present — user should use one or
                # the other depending on the geometry type.
                import warnings
                if wr_splits and (wr_sectors or wr_angles):
                    warnings.warn(
                        "water_rods block contains both 'splits' and "
                        "'sectors'/'angles'.  For circular water rods "
                        "only sectors/angles are used; for square water "
                        "rods only splits is used.  The unused keys will "
                        "be ignored at geometry-build time.",
                        stacklevel=2,
                    )

                water_rod_sectors = SectorConfig(
                    sectors=wr_sectors,
                    angles=wr_angles,
                    windmill=sect.get("windmill", False),
                    splits=wr_splits,
                    additional_radial_splits_in_moderator=wr_additional,
                    subdivisions_coolant_corners=wr_corners,
                )
            van_rod = sect.get("vanished_rods", {})
            if van_rod:
                vanished_rods_sectors = VanishedRodDiscretizationConfig(
                    base_radius=van_rod.get("base_radius", None),
                    radial_splits=van_rod.get("radial_splits", 1),
                    sectors=van_rod.get("sectors", []),
                    angles=van_rod.get("angles"),
                    windmill=sect.get("windmill", False),
                )

        # --- Radial overrides ---
        radial_overrides_raw = d.get("radial_overrides", {})
        radial_overrides = {}
        for key, val in radial_overrides_raw.items():
            radial_overrides[key] = {
                "scheme": val.get("scheme", d.get("radial_scheme", "Santamarina")),
                "params": val.get("params", {}),
            }

        # --- Box discretization (MOC sub-meshing) ---
        box_disc_raw = d.get("box_discretization", {})
        box_disc = None
        if box_disc_raw.get("enabled", False):
            # Parse optional cross_moderator_discretization sub-block
            cross_mod_disc = None
            cross_raw = box_disc_raw.get(
                "cross_moderator_discretization", {}
            )
            if cross_raw:
                # Accept both old "narrow_gap_splits" and new "cross_side_gap_splits" with deprecation
                cross_side_gap_splits_val = cross_raw.get("cross_side_gap_splits")
                narrow_gap_splits_val = cross_raw.get("narrow_gap_splits")

                if narrow_gap_splits_val is not None and cross_side_gap_splits_val is None:
                    import warnings
                    warnings.warn(
                        "YAML key 'cross_moderator_discretization.narrow_gap_splits' is deprecated. "
                        "Use 'cross_side_gap_splits' instead.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                    cross_side_gap_splits_val = narrow_gap_splits_val

                cross_mod_disc = CrossModeratorDiscretizationConfig(
                    cross_side_gap_splits=cross_side_gap_splits_val,
                    moderator_at_cross_corner_splits=cross_raw.get(
                        "moderator_at_cross_corner_splits"
                    ),
                    stub_splits=cross_raw.get("stub_splits"),
                )

            # Parse optional control_cross_submesh sub-block (sibling)
            ctrl_cross_raw = box_disc_raw.get(
                "control_cross_submesh", False
            )
            if isinstance(ctrl_cross_raw, (bool, dict)):
                ctrl_cross_cfg = ControlCrossSubmeshConfig.from_yaml(
                    ctrl_cross_raw
                )
            else:
                ctrl_cross_cfg = None

            box_disc = BoxDiscretizationConfig(
                enabled=True,
                corner_splits=box_disc_raw.get("corner_splits", None),
                gap_splits=box_disc_raw.get("gap_splits", None),
                gap_wide_splits=box_disc_raw.get("gap_wide_splits", None),
                gap_narrow_splits=box_disc_raw.get("gap_narrow_splits", None),
                wide_wide_corner_splits=box_disc_raw.get("wide_wide_corner_splits", None),
                narrow_narrow_corner_splits=box_disc_raw.get("narrow_narrow_corner_splits", None),
                mixed_corner_splits=box_disc_raw.get("mixed_corner_splits", None),
                cross_moderator_discretization=cross_mod_disc,
                control_cross_submesh=ctrl_cross_cfg,
                reassign_materials=box_disc_raw.get(
                    "reassign_materials", True),
                # deprecated aliases (backward compatibility)
                side_x_splits=box_disc_raw.get("side_x_splits", None),
                side_y_splits=box_disc_raw.get("side_y_splits", None),
            )

        # --- Tracking parameters (common to all flux step types) ---
        tracking_kwargs = dict(
            num_angles_2d=d.get("num_angles_2d", 8),
            line_density=d.get("line_density", 10.0),
            anisotropy_level=d.get("anisotropy_level", 1),
            polar_angles_quadrature=d.get(
                "polar_angles_quadrature", "GAUS"),
            number_of_polar_angles=d.get(
                "number_of_polar_angles", 1),
        )

        # --- Build CalculationStep with appropriate parameters ---
        step_type = d["step_type"]

        # Common kwargs shared by both branches
        common_kwargs = dict(
            tracking=d.get("tracking", "TISO"),
            batch_size=d.get("batch_size", None),
            flux_level=d.get("flux_level", None),
            radial_scheme=d.get("radial_scheme", "Santamarina"),
            radial_params=d.get("radial_params", {}),
            radial_overrides=radial_overrides,
            sectorization_enabled=sect_enabled,
            fuel_sectors=fuel_sectors,
            gd_sectors=gd_sectors,
            water_rod_sectors=water_rod_sectors,
            vanished_rods_sectors=vanished_rods_sectors,
            export_macros=d.get("export_macros", False),
            box_discretization=box_disc,
            mix_numbering_strategy=d.get("mix_numbering_strategy", "by_material"),
            tdt_file_id=d.get("tdt_file_id", None),
            macro_assignment = d.get("macro_assignment", None),
            **tracking_kwargs,
        )

        if step_type == "self_shielding":
            for key in ("leakage_model", "leakage_correction"):
                if key in d:
                    raise ValueError(
                        f"'{key}' is only valid for flux steps, but was "
                        f"given for self_shielding step '{d['name']}'."
                    )
            return CalculationStep(
                name=d["name"],
                step_type=step_type,
                self_shielding_module=d.get(
                    "self_shielding_module", "USS"),
                self_shielding_method=d.get(
                    "self_shielding_method", "RSE"),
                spatial_method=d.get("spatial_method", "CP"),
                transport_correction=d.get("transport_correction", None),
                fuel_self_shielded_isotopes=d.get("fuel_self_shielded_isotopes", None),
                clad_self_shielded_isotopes=d.get("clad_self_shielded_isotopes", None),
                **common_kwargs,
            )
        else:
            return CalculationStep(
                name=d["name"],
                step_type=step_type,
                spatial_method=d.get("spatial_method", "CP"),
                leakage_model=d.get("leakage_model", None),
                leakage_correction=d.get("leakage_correction", None),
                **common_kwargs,
            )

    # ------------------------------------------------------------------
    # Presets
    # ------------------------------------------------------------------

    @classmethod
    def preset(cls, preset_name):
        """
        Return a pre-built calculation scheme.

        Available presets:

        ``"BWR_fine_1L"``
            Two-step scheme for BWR assemblies with fine MOC flux:

            1. **SSH**: USS module, RSE self-shielding method, IC spatial
               method, Santamarina radii, no sectorization, TISO tracking,
               macro export enabled.
            2. **FLUX**: MOC spatial method, TSPC tracking, Santamarina radii,
               full windmill sectorization (8 sectors on all regions),
               box discretization enabled.

        ``"BWR_2L"``
            Three-step scheme with two flux levels:

            1. **SSH**: USS module, RSE self-shielding method, IC spatial
               method, Santamarina radii, no sectorization, TISO tracking,
               macro export enabled.
            2. **FLUX_L1**: IC spatial method, TISO tracking, flux_level=1,
               Santamarina radii, windmill sectorization (8 sectors on
               coolant only), macro export enabled.
            3. **FLUX_L2**: MOC spatial method, TSPC tracking, flux_level=2,
               Santamarina radii, full windmill sectorization (8 sectors
               on all regions), box discretization enabled, macro export.

        ``"BWR_CP"``
            Two-step CP-based scheme:

            1. **SSH**: USS module, PT self-shielding method, CP spatial
               method, Santamarina radii, no sectorization, TISO tracking.
            2. **FLUX**: CP spatial method, TSPC tracking, Santamarina radii,
               no sectorization.

        Parameters
        ----------
        preset_name : str

        Returns
        -------
        DragonCalculationScheme
        """
        builders = {
            "BWR_fine_1L": cls._preset_BWR_fine_1L,
            "BWR_2L": cls._preset_BWR_2L,
            "BWR_CP": cls._preset_BWR_CP,
        }
        if preset_name not in builders:
            raise ValueError(
                f"Unknown preset '{preset_name}'. "
                f"Available presets: {list(builders.keys())}"
            )
        return builders[preset_name]()

    @classmethod
    def _preset_BWR_fine_1L(cls):
        scheme = cls(name="BWR_fine_1L")

        # Step 1: Self-shielding (USS + RSE + IC, no sectorization)
        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
            export_macros=True,
        ))

        # Step 2: Flux (MOC, windmill sectorization)
        scheme.add_step(CalculationStep(
            name="FLUX",
            step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            radial_scheme="Santamarina",
            sectorization_enabled=True,
            export_macros=False,
            num_angles_2d=24,
            line_density=40.0,
            anisotropy_level=4,
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
            fuel_sectors=SectorConfig(
                sectors=[8, 8, 8, 8, 8, 8, 8],
                angles=[22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5],
                windmill=True,
            ),
            gd_sectors=SectorConfig(
                sectors=[8, 8, 8, 8, 8, 8, 8, 8, 8],
                angles=[22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5],
                windmill=True,
            ),
            water_rod_sectors=SectorConfig(
                sectors=[8, 8, 8],
                angles=[22.5, 22.5, 22.5],
                windmill=True,
            ),
            box_discretization=BoxDiscretizationConfig(
                enabled=True,
                corner_splits=(8, 8),
                gap_splits=(30, 8),
            ),
        ))

        return scheme

    @classmethod
    def _preset_BWR_2L(cls):
        scheme = cls(name="BWR_2L")

        # Step 1: Self-shielding
        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="RSE",
            spatial_method="IC",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
            export_macros=True,
        ))

        # Step 2: Flux level 1 (IC, windmill)
        scheme.add_step(CalculationStep(
            name="FLUX_L1",
            step_type="flux",
            spatial_method="IC",
            tracking="TISO",
            flux_level=1,
            radial_scheme="Santamarina",
            sectorization_enabled=True,
            export_macros=True,
            fuel_sectors=SectorConfig(
                sectors=[1, 1, 1, 1, 1, 1, 8],
                angles=[0, 0, 0, 0, 0, 0, 22.5],
                windmill=True,
            ),
            gd_sectors=SectorConfig(
                sectors=[1, 1, 1, 1, 1, 1, 1, 1, 8],
                angles=[0, 0, 0, 0, 0, 0, 0, 0, 22.5],
                windmill=True,
            ),
            water_rod_sectors=SectorConfig(
                sectors=[1, 1, 8],
                angles=[0, 0, 0],
                windmill=True,
            ),
        ))

        # Step 3: Flux level 2 (MOC, finer sectorization, TSPC)
        scheme.add_step(CalculationStep(
            name="FLUX_L2",
            step_type="flux",
            spatial_method="MOC",
            tracking="TSPC",
            flux_level=2,
            radial_scheme="Santamarina",
            sectorization_enabled=True,
            export_macros=True,
            num_angles_2d=24,
            line_density=40.0,
            anisotropy_level=4,
            polar_angles_quadrature="GAUS",
            number_of_polar_angles=4,
            fuel_sectors=SectorConfig(
                sectors=[8, 8, 8, 8, 8, 8, 8],
                angles=[22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5],
                windmill=True,
            ),
            gd_sectors=SectorConfig(
                sectors=[8, 8, 8, 8, 8, 8, 8, 8, 8],
                angles=[22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5],
                windmill=True,
            ),
            water_rod_sectors=SectorConfig(
                sectors=[8, 8, 8],
                angles=[22.5, 22.5, 22.5],
                windmill=True,
            ),
            box_discretization=BoxDiscretizationConfig(
                enabled=True,
                corner_splits=(8, 8),
                gap_splits=(30, 8),
            ),
        ))

        return scheme

    @classmethod
    def _preset_BWR_CP(cls):
        scheme = cls(name="BWR_CP")

        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_module="USS",
            self_shielding_method="PT",
            spatial_method="CP",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
        ))

        scheme.add_step(CalculationStep(
            name="FLUX",
            step_type="flux",
            spatial_method="CP",
            tracking="TSPC",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
        ))

        return scheme

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self):
        step_names = [s.name for s in self.steps]
        return f"DragonCalculationScheme(name='{self.name}', steps={step_names})"

    def summary(self):
        """Return a human-readable multi-line summary of the scheme."""
        lines = [f"Calculation Scheme: {self.name}", "=" * 50]
        for i, step in enumerate(self.steps, 1):
            lines.append(f"\nStep {i}: {step.name}")
            lines.append(f"  Type:       {step.step_type}")
            if step.step_type == "edition_between_levels":
                lines.append(f"  Groups:     {step.number_of_macro_groups}")
                lines.append(f"  Group bounds: {step.energy_groups_bounds}")
                lines.append(f"  SPH:        {step.sph_correction}")
                if step.sph_correction:
                    lines.append(f"  GRMAX:      {step.max_sph_group}")
                continue
            if step.step_type == "self_shielding":
                lines.append(f"  SSH module:  {step.self_shielding_module}")
                lines.append(f"  SSH method:  {step.self_shielding_method}")
            lines.append(f"  Method:     {step.spatial_method}")
            if step.step_type == "flux":
                leak = step.leakage_model or "none"
                if step.leakage_model and step.leakage_correction:
                    leak += f" {step.leakage_correction}"
                lines.append(f"  Leakage:    {leak}")
            lines.append(f"  Tracking:   {step.tracking}")
            if step.flux_level is not None:
                lines.append(f"  Flux level: {step.flux_level}")
            lines.append(f"  Radial:     {step.radial_scheme}")
            if step.radial_params:
                lines.append(f"  Radial params: {step.radial_params}")
            if step.radial_overrides:
                lines.append(f"  Radial overrides: {list(step.radial_overrides.keys())}")
            lines.append(f"  Sectors:    {'enabled' if step.sectorization_enabled else 'disabled'}")
            if step.sectorization_enabled:
                if step.fuel_sectors:
                    lines.append(f"    fuel:  {step.fuel_sectors}")
                if step.gd_sectors:
                    lines.append(f"    Gd:    {step.gd_sectors}")
                if step.water_rod_sectors:
                    lines.append(f"    WRod:  {step.water_rod_sectors}")
            lines.append(f"  Macros:     {step.export_macros}")
            lines.append(f"  Box disc:   {'enabled' if step.box_discretization and step.box_discretization.enabled else 'disabled'}")
            if step.box_discretization and step.box_discretization.enabled:
                lines.append(f"    {step.box_discretization}")

        if self.branches:
            lines.append("\nCalculation Branches")
            lines.append("-" * 50)
            for b in self.branches:
                lines.append(f"  {b.type}: {b.values}")
            lines.append(f"  Total statepoints: {self.get_total_statepoints()}")
            lines.append(f"  Loop order: {self.branch_loop_order}")

        if self.outputs:
            lines.append("\nCalculation Outputs")
            lines.append("-" * 50)
            for o in self.outputs:
                lines.append(f"  {o.name}: mode={o.spatial_integration_mode}, "
                             f"bounds={o.energy_bounds}, "
                             f"state_points={o.state_points}")

        return "\n".join(lines)
