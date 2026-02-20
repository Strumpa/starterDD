# Classes defining a DRAGON calculation scheme for starterDD.
# A calculation scheme orchestrates the sequence of DRAGON calculation steps
# (self-shielding, flux, depletion) and their respective geometry
# discretization options (radial subdivision, azimuthal sectorization).
#
# R.Guasch
# Date : 19/02/2026 [created]

import yaml
from copy import deepcopy


# ---------------------------------------------------------------------------
# Valid option constants
# ---------------------------------------------------------------------------

VALID_STEP_TYPES = ("self_shielding", "flux")
VALID_SELF_SHIELDING_METHODS = ("USS", "SHI", "RSE", "None")  # USS: subgroup, SHI: Stammler, RSE: subgroup+equivalence, None: no self-shielding
VALID_SPATIAL_METHODS = ("CP", "IC", "MOC")
VALID_TRACKING_OPTIONS = ("TISO", "TSPC")
VALID_RADIAL_SCHEMES = ("Santamarina", "automatic", "user_defined")


# ---------------------------------------------------------------------------
# SectorConfig – lightweight container for sectorization parameters
# ---------------------------------------------------------------------------

class SectorConfig:
    """
    Container for azimuthal sectorization parameters used by glow's
    ``RectCell.sectorize()`` method.

    Attributes
    ----------
    sectors : list[int]
        Number of sectors per radial ring (from innermost to outermost
        region, including coolant).
    angles : list[float]
        Starting angle offset for each ring (degrees).
    windmill : bool
        Whether to apply windmill sectorization.
    """

    def __init__(self, sectors=None, angles=None, windmill=False):
        self.sectors = sectors or []
        self.angles = angles or []
        self.windmill = windmill

    def __repr__(self):
        return (f"SectorConfig(sectors={self.sectors}, "
                f"angles={self.angles}, windmill={self.windmill})")


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
    * **side_x_splits** ``(nx, ny)`` — grid size for horizontal (top /
      bottom) side strips.  ``nx`` defaults to the number of pin columns.
    * **side_y_splits** ``(nx, ny)`` — grid size for vertical (left /
      right) side strips.  ``ny`` defaults to the number of pin rows.

    Attributes
    ----------
    enabled : bool
        Whether box discretization is active.
    corner_splits : tuple[int, int]
        ``(nx, ny)`` grid subdivisions for each corner region.
    side_x_splits : tuple[int, int]
        ``(nx, ny)`` grid subdivisions for horizontal side strips.
    side_y_splits : tuple[int, int]
        ``(nx, ny)`` grid subdivisions for vertical side strips.
    """

    def __init__(self, enabled=False, corner_splits=None,
                 side_x_splits=None, side_y_splits=None):
        self.enabled = enabled
        self.corner_splits = tuple(corner_splits) if corner_splits else (4, 4)
        self.side_x_splits = tuple(side_x_splits) if side_x_splits else None
        self.side_y_splits = tuple(side_y_splits) if side_y_splits else None

    def resolve_splits(self, n_cols, n_rows):
        """
        Return the resolved ``(corner, side_x, side_y)`` split tuples,
        filling in lattice-dimension defaults for ``None`` values.

        Parameters
        ----------
        n_cols : int
            Number of pin columns in the lattice.
        n_rows : int
            Number of pin rows in the lattice.

        Returns
        -------
        corner : tuple[int, int]
        side_x : tuple[int, int]
            Splits for top/bottom strips (default ``(n_cols, 1)``).
        side_y : tuple[int, int]
            Splits for left/right strips (default ``(1, n_rows)``).
        """
        corner = self.corner_splits
        side_x = self.side_x_splits if self.side_x_splits else (n_cols, 1)
        side_y = self.side_y_splits if self.side_y_splits else (1, n_rows)
        return corner, side_x, side_y

    def __repr__(self):
        return (f"BoxDiscretizationConfig(enabled={self.enabled}, "
                f"corner_splits={self.corner_splits}, "
                f"side_x_splits={self.side_x_splits}, "
                f"side_y_splits={self.side_y_splits})")


# ---------------------------------------------------------------------------
# CalculationStep – one step inside a calculation scheme
# ---------------------------------------------------------------------------

class CalculationStep:
    """
    Describes a single step in a DRAGON calculation scheme.

    A step combines:

    * **Step type** (``"self_shielding"`` or ``"flux"``).
    * **Method pair**: a self-shielding/energy method (``"USS"`` or ``"SHI"``)
      coupled with a spatial solution method (``"CP"``, ``"IC"``, ``"MOC"``).
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
    self_shielding_method : str
        Energy self-shielding method: ``"USS"`` or ``"SHI"``.
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
    export_macros : bool
        Whether to export MACRO properties in the geometry.
    box_discretization : BoxDiscretizationConfig or None
        Configuration for sub-meshing the assembly-box peripheral
        regions into a grid of sub-faces (for MOC tracking).  When
        ``None``, no box discretization is applied.
    """

    def __init__(
        self,
        name,
        step_type,
        self_shielding_method,
        spatial_method,
        tracking="TISO",
        flux_level=None,
        radial_scheme="Santamarina",
        radial_params=None,
        radial_overrides=None,
        sectorization_enabled=False,
        fuel_sectors=None,
        gd_sectors=None,
        water_rod_sectors=None,
        export_macros=False,
        box_discretization=None,
    ):
        # --- Validate step type ---
        if step_type not in VALID_STEP_TYPES:
            raise ValueError(
                f"Invalid step_type '{step_type}'. "
                f"Valid options: {VALID_STEP_TYPES}"
            )

        # --- Validate method pair ---
        if self_shielding_method not in VALID_SELF_SHIELDING_METHODS:
            raise ValueError(
                f"Invalid self_shielding_method '{self_shielding_method}'. "
                f"Valid options: {VALID_SELF_SHIELDING_METHODS}"
            )
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

        self.name = name
        self.step_type = step_type
        self.self_shielding_method = self_shielding_method
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
        self.export_macros = export_macros
        self.box_discretization = box_discretization

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

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self):
        level_str = f", flux_level={self.flux_level}" if self.flux_level else ""
        return (
            f"CalculationStep(name='{self.name}', "
            f"type='{self.step_type}', "
            f"method=({self.self_shielding_method}, {self.spatial_method}), "
            f"tracking='{self.tracking}'{level_str})"
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

    # ------------------------------------------------------------------
    # Step management
    # ------------------------------------------------------------------

    def add_step(self, step):
        """
        Append a ``CalculationStep`` to the scheme.

        Parameters
        ----------
        step : CalculationStep
        """
        if not isinstance(step, CalculationStep):
            raise TypeError(f"Expected CalculationStep, got {type(step).__name__}")
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

    # ------------------------------------------------------------------
    # YAML construction
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, yaml_path):
        """
        Parse a DRAGON calculation scheme from a YAML file.

        Expected YAML structure::

            DRAGON_CALCULATION_SCHEME:
              name: "GE14_standard"
              steps:
                - name: "SSH"
                  step_type: "self_shielding"
                  self_shielding_method: "USS"
                  spatial_method: "CP"
                  tracking: "TISO"
                  radial_scheme: "Santamarina"
                  sectorization:
                    enabled: false

                - name: "FLUX"
                  step_type: "flux"
                  self_shielding_method: "USS"
                  spatial_method: "IC"
                  tracking: "TISO"
                  flux_level: null
                  radial_scheme: "Santamarina"
                  export_macros: true
                  sectorization:
                    enabled: true
                    windmill: true
                    fuel_pins:
                      sectors: [1, 1, 1, 1, 1, 1, 8]
                      angles:  [0, 0, 0, 0, 0, 0, 22.5]
                    Gd_pins:
                      sectors: [1, 1, 1, 1, 1, 1, 1, 1, 8]
                      angles:  [0, 0, 0, 0, 0, 0, 0, 0, 22.5]
                    water_rods:
                      sectors: [1, 1, 8]
                      angles:  [0, 0, 0]

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

        scheme_data = data.get("DRAGON_CALCULATION_SCHEME", data)
        scheme_name = scheme_data.get("name", "from_yaml")
        scheme = cls(name=scheme_name)

        for step_data in scheme_data.get("steps", []):
            step = cls._parse_step(step_data)
            scheme.add_step(step)

        return scheme

    @staticmethod
    def _parse_step(d):
        """Parse a single step dict into a CalculationStep."""
        # --- Sectorization ---
        sect = d.get("sectorization", {})
        sect_enabled = sect.get("enabled", False)

        fuel_sectors = None
        gd_sectors = None
        water_rod_sectors = None

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
                water_rod_sectors = SectorConfig(
                    sectors=wr.get("sectors", []),
                    angles=wr.get("angles", []),
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
            box_disc = BoxDiscretizationConfig(
                enabled=True,
                corner_splits=box_disc_raw.get("corner_splits", None),
                side_x_splits=box_disc_raw.get("side_x_splits", None),
                side_y_splits=box_disc_raw.get("side_y_splits", None),
            )

        return CalculationStep(
            name=d["name"],
            step_type=d["step_type"],
            self_shielding_method=d.get("self_shielding_method", "USS"),
            spatial_method=d.get("spatial_method", "CP"),
            tracking=d.get("tracking", "TISO"),
            flux_level=d.get("flux_level", None),
            radial_scheme=d.get("radial_scheme", "Santamarina"),
            radial_params=d.get("radial_params", {}),
            radial_overrides=radial_overrides,
            sectorization_enabled=sect_enabled,
            fuel_sectors=fuel_sectors,
            gd_sectors=gd_sectors,
            water_rod_sectors=water_rod_sectors,
            export_macros=d.get("export_macros", False),
            box_discretization=box_disc,
        )

    # ------------------------------------------------------------------
    # Presets
    # ------------------------------------------------------------------

    @classmethod
    def preset(cls, preset_name):
        """
        Return a pre-built calculation scheme.

        Available presets:

        ``"BWR_standard"``
            Two-step scheme for BWR assemblies:

            1. **SSH**: USS self-shielding with CP spatial method,
               Santamarina radii, no sectorization, TISO tracking.
            2. **FLUX**: USS self-shielding with IC spatial method,
               Santamarina radii, windmill sectorization
               (8 coolant sectors), TISO tracking, macro export.

        ``"BWR_fine"``
            Three-step scheme with two flux levels:

            1. **SSH**: same as BWR_standard.
            2. **FLUX_L1**: IC flux on Santamarina radii with windmill.
            3. **FLUX_L2**: MOC flux on Santamarina radii with finer
               sectorization, TSPC tracking.

        ``"BWR_CP"``
            Two-step CP-based scheme:

            1. **SSH**: USS + CP, Santamarina, no sectors, TISO.
            2. **FLUX**: USS + CP, Santamarina, no sectors, TISO.

        Parameters
        ----------
        preset_name : str

        Returns
        -------
        DragonCalculationScheme
        """
        builders = {
            "BWR_standard": cls._preset_BWR_standard,
            "BWR_fine": cls._preset_BWR_fine,
            "BWR_CP": cls._preset_BWR_CP,
        }
        if preset_name not in builders:
            raise ValueError(
                f"Unknown preset '{preset_name}'. "
                f"Available presets: {list(builders.keys())}"
            )
        return builders[preset_name]()

    @classmethod
    def _preset_BWR_standard(cls):
        scheme = cls(name="BWR_standard")

        # Step 1: Self-shielding (USS + CP, no sectorization)
        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_method="USS",
            spatial_method="CP",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
        ))

        # Step 2: Flux (USS + IC, windmill sectorization)
        scheme.add_step(CalculationStep(
            name="FLUX",
            step_type="flux",
            self_shielding_method="USS",
            spatial_method="IC",
            tracking="TISO",
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

        return scheme

    @classmethod
    def _preset_BWR_fine(cls):
        scheme = cls(name="BWR_fine")

        # Step 1: Self-shielding
        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_method="USS",
            spatial_method="CP",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
        ))

        # Step 2: Flux level 1 (IC, windmill)
        scheme.add_step(CalculationStep(
            name="FLUX_L1",
            step_type="flux",
            self_shielding_method="USS",
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
            self_shielding_method="USS",
            spatial_method="MOC",
            tracking="TSPC",
            flux_level=2,
            radial_scheme="Santamarina",
            sectorization_enabled=True,
            export_macros=True,
            fuel_sectors=SectorConfig(
                sectors=[1, 1, 1, 1, 1, 1, 12],
                angles=[0, 0, 0, 0, 0, 0, 15.0],
                windmill=True,
            ),
            gd_sectors=SectorConfig(
                sectors=[1, 1, 1, 1, 1, 1, 1, 1, 12],
                angles=[0, 0, 0, 0, 0, 0, 0, 0, 15.0],
                windmill=True,
            ),
            water_rod_sectors=SectorConfig(
                sectors=[1, 1, 12],
                angles=[0, 0, 0],
                windmill=True,
            ),
            box_discretization=BoxDiscretizationConfig(
                enabled=True,
                corner_splits=(4, 4),
                side_x_splits=None,  # defaults to (n_cols, 1)
                side_y_splits=None,  # defaults to (1, n_rows)
            ),
        ))

        return scheme

    @classmethod
    def _preset_BWR_CP(cls):
        scheme = cls(name="BWR_CP")

        scheme.add_step(CalculationStep(
            name="SSH",
            step_type="self_shielding",
            self_shielding_method="USS",
            spatial_method="CP",
            tracking="TISO",
            radial_scheme="Santamarina",
            sectorization_enabled=False,
        ))

        scheme.add_step(CalculationStep(
            name="FLUX",
            step_type="flux",
            self_shielding_method="USS",
            spatial_method="CP",
            tracking="TISO",
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
            if step.step_type == "self_shielding":
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
            print(f"  Box disc:   {'enabled' if step.box_discretization and step.box_discretization.enabled else 'disabled'}")
            if step.box_discretization and step.box_discretization.enabled:
                lines.append(f"  Box disc:   {step.box_discretization}")
        return "\n".join(lines)
