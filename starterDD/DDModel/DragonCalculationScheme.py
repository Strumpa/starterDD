# Classes defining a DRAGON calculation scheme for starterDD.
# A calculation scheme orchestrates the sequence of DRAGON calculation steps
# (self-shielding, flux) and their respective geometry
# discretization options (radial subdivision, azimuthal sectorization).
#
# R.Guasch
# Date : 19/02/2026 [created]

import yaml
from copy import deepcopy


# ---------------------------------------------------------------------------
# Valid option constants
# ---------------------------------------------------------------------------

VALID_STEP_TYPES = ("self_shielding", "flux", "edition_between_levels")
VALID_SELF_SHIELDING_MODULES = ("USS", "SHI")  # USS: unresolved resonance, SHI: Stammler
VALID_SELF_SHIELDING_METHODS = ("RSE", "PT")  # RSE: subgroup+equivalence, PT: probability tables
VALID_SPATIAL_METHODS = ("CP", "IC", "MOC")
VALID_TRACKING_OPTIONS = ("TISO", "TSPC")
VALID_RADIAL_SCHEMES = ("Santamarina", "automatic", "user_defined")
VALID_MIX_NUMBERING_STRATEGIES = ("by_material", "by_pin")


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
        Default sectorization config for standard fuel pins.
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

    def get_sectorization_for_pin(self, pin_or_rod_type, isGd=False):
        """
        Return the ``SectorConfig`` applicable to a given pin.

        Parameters
        ----------
        pin_or_rod_type : str or FuelPinModel
            Either a rod type string or a pin model object.
        isGd : bool
            Whether the pin is a gadolinium-bearing pin (used when
            ``pin_or_rod_type`` is a string).

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
        return (
            f"CalculationStep(name='{self.name}', "
            f"type='{self.step_type}', "
            f"{method_str}, "
            f"tracking='{self.tracking}', "
            f"mix_numbering='{self.mix_numbering_strategy}'{level_str})"
        )


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
        Parameter values to iterate over.
    """

    def __init__(self, name, branch_type, values):
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
        self.values = list(values)

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
                k: list(v) for k, v in state_points.items()
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
            if current not in required_values:
                return False
        return True

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
    """

    step_type = "edition_between_levels"

    def __init__(self, name, number_of_macro_groups,
                 energy_groups_bounds, sph_correction=False,
                 max_sph_group=None):
        self.name = name
        self.number_of_macro_groups = number_of_macro_groups
        self.energy_groups_bounds = list(energy_groups_bounds)
        self.sph_correction = bool(sph_correction)
        self.max_sph_group = max_sph_group

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
            f"sph={self.sph_correction})"
        )


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

    def __init__(self, name="default"):
        self.name = name
        self.steps = []
        self.branches = []   # list[CalculationBranch]
        self.outputs = []    # list[CalculationOutput]
        self.branch_loop_order = list(DEFAULT_BRANCH_LOOP_ORDER)

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

        # Accept multiple top-level key names for the flux scheme
        scheme_data = (
            data.get("DRAGON_CALCULATION_SCHEME")
            or data.get("FLUX_CALCULATION_SCHEME")
            or data
        )
        scheme_name = scheme_data.get("name", "from_yaml")
        scheme = cls(name=scheme_name)

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
            **tracking_kwargs,
        )

        if step_type == "self_shielding":
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
                lines.append(f"  SH module:  {step.self_shielding_module}")
                lines.append(f"  SH method:  {step.self_shielding_method}")
            lines.append(f"  Method:     {step.spatial_method}")
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
