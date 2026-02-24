## Some helpful functions to deal with GLOW geometry building
# Author : R. Guasch
# Date : 04/02/2026
# ----------------------------------------------------------------------------


from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.geometries import Rectangle
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.main import TdtSetup, analyse_and_generate_tdt
from glow.interface.geom_interface import *
from glow.support.types import *
from .helpers import computeSantamarinaradii
# Note: CartesianAssemblyModel and FuelPinModel are imported inside functions to avoid circular imports


def make_grid_faces(parent: Rectangle, nx: int, ny: int):
    """
    Create a regular nx by ny grid of faces within the given parent rectangle, returning a list of the created faces.
    """
    lx = float(parent.lx) 
    ly = float(parent.ly)  
    dx = lx / nx
    dy = ly / ny
    cx_parent, cy_parent = float(parent.o.GetParameters().split(":")[0]), float(parent.o.GetParameters().split(":")[1])
    x0 = cx_parent - lx / 2.0
    y0 = cy_parent - ly / 2.0

    # create a (nx+1) x (ny+1) grid of vertices
    verts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = x0 + i * dx
            y = y0 + j * dy
            v = make_vertex((x, y, 0.0))
            verts.append(v)

    # helper to index vertex at grid (i,j)
    def v_at(i, j):
        return verts[j * (nx + 1) + i]

    faces = []
    for j in range(ny):
        for i in range(nx):
            # rectangle corners (lower-left, lower-right, upper-right, upper-left)
            v00 = v_at(i, j)
            v10 = v_at(i + 1, j)
            v11 = v_at(i + 1, j + 1)
            v01 = v_at(i, j + 1)

            # create four edges for the rectangle (order matters for consistent orientation)
            e_bottom = make_edge(v00, v10)
            e_right = make_edge(v10, v11)
            e_top = make_edge(v11, v01)
            e_left = make_edge(v01, v00)

            # assemble a face from the four edges
            face = make_face([e_bottom, e_right, e_top, e_left])

            faces.append(face)

    return faces


def generate_fuel_cells(assemblyModel, calculation_step=None):
    """
    Generate RectCell objects for each individual subgeometry in the lattice.

    Parameters
    ----------
    assemblyModel : CartesianAssemblyModel
        The assembly model containing the lattice description and rod ID to material mapping.
        Additionally contains the pin geometry parameters needed to define the fuel cells
        and the material mixture unique names to assign to the cell materials.
    calculation_step : CalculationStep or None
        Optional calculation step providing sectorization configuration.
        When provided, each fuel cell is sectorized according to the step's
        ``SectorConfig`` for fuel / Gd / water-rod pins.  When ``None``,
        no sectorization is applied (backward-compatible behaviour).

    Returns
    -------
    lattice_components : list of list of RectCell
        2D list of RectCell objects representing the lattice layout, with properties set
        according to assemblyModel information. To be used to build lattice geometry with
        glow and export to TDT file.
    """
    # Import here to avoid circular import issues
    from ..DDModel.DragonModel import FuelPinModel
    
    lattice_components = []
    pitch = assemblyModel.pin_geometry_dict["pin_pitch"]
    
    fuel_material_mixtures = assemblyModel.fuel_material_mixtures
    row_idx = -1
    for row in assemblyModel.lattice:
        row_idx += 1 
        row_of_cells = []
        cell_idx = -1
        for pin in row:
            cell_idx += 1
            if isinstance(pin, FuelPinModel):
                fuel_material_mixtures = pin.fuel_material_mixtures
                radii = pin.radii
                tmp_cell = RectCell(
                    name=pin.fuel_material_name,
                    height_x_width=(pitch, pitch),
                    center=(0.0, 0.0, 0.0),
                )
                for radius in radii:
                    tmp_cell.add_circle(radius)
                list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["GAP", "CLAD", "COOLANT"]

                # Apply sectorization from calculation step if provided
                if calculation_step is not None:
                    sector_cfg = calculation_step.get_sectorization_for_pin(pin, isGd=pin.isGd)
                    if sector_cfg is not None:
                        tmp_cell.sectorize(sector_cfg.sectors, sector_cfg.angles, windmill=sector_cfg.windmill)

                tmp_cell.set_properties({
                    PropertyType.MATERIAL: list_of_cell_mats,
                    PropertyType.MACRO: [f"MACRO{row_idx}{cell_idx}"] * len(list_of_cell_mats)
                })
            else:  # Water rod placeholder or other non-fuel cell
                tmp_cell = RectCell(
                    name=pin.rod_ID,
                    height_x_width=(pitch, pitch),
                    center=(0.0, 0.0, 0.0),
                )
                tmp_cell.set_properties({
                    PropertyType.MATERIAL: ["MODERATOR"],
                    PropertyType.MACRO: [f"MACRO_{row_idx}{cell_idx}"]
                })
            print(f"Generated cell {tmp_cell.name} at position ({cell_idx}, {row_idx})")
            row_of_cells.append(tmp_cell)
        lattice_components.append(row_of_cells)
    return lattice_components


def add_cells_to_regular_lattice(lattice, ordered_cells, cell_pitch, translation=0.0):
    """
    Add fuel cells to the lattice, skipping water rod placeholders.

    Parameters
    ----------
    lattice : Lattice
        The lattice to which cells will be added
    ordered_cells : list of list of RectCell
        2D list of RectCell objects representing the lattice layout
    cell_pitch : float
        Pitch of each cell in the lattice
    translation : float
        Translation to apply to cell positions
    """
    for row_idx in range(len(ordered_cells)):
        row_of_cells = ordered_cells[row_idx]
        for cell_idx in range(len(row_of_cells)):
            cell = row_of_cells[cell_idx]
            if "W" in cell.name:  # Skip water rod placeholders
                continue
            else:
                print(f"Adding cell {cell.name} at position ({cell_idx}, {row_idx})")
                print(f"cell_pitch: {cell_pitch}, translation: {translation}")
                lattice.add_cell(
                    cell, ((cell_idx + 0.5) * cell_pitch + translation,
                           (row_idx + 0.5) * cell_pitch + translation,
                           0.0)
                )
    return lattice


def _build_square_water_rod_cell(water_rod_model, calculation_step=None):
    """
    Build a ``RectCell`` for a square water rod with 3 concentric
    rectangular regions (moderator / cladding / coolant).

    The construction follows the same pattern as ``build_assembly_box``:
    two inner ``Rectangle`` boundaries are partitioned into the bounding
    box cell, then materials are assigned from innermost to outermost.

    If the ``calculation_step`` provides a ``SectorConfig`` with a
    ``splits`` attribute, the cell is further sub-meshed into an
    ``(nx, ny)`` Cartesian grid.  Material is then reassigned to each
    sub-face by geometric containment against the two boundary
    rectangles.

    Parameters
    ----------
    water_rod_model : SquareWaterRodModel
        The square water rod model with ``bounding_box_side_length``,
        ``moderator_box_inner_side``, ``moderator_box_outer_side``,
        ``center``, ``rod_ID``, and material names.
    calculation_step : CalculationStep or None
        Optional calculation step providing discretization config via
        ``get_water_rod_sectorization().splits``.

    Returns
    -------
    tmp_cell : RectCell
        The constructed cell ready to be added to the lattice.
    """
    import warnings

    bb = water_rod_model.bounding_box_side_length
    inner_side = water_rod_model.moderator_box_inner_side
    outer_side = water_rod_model.moderator_box_outer_side
    corner_radius = getattr(water_rod_model, 'corner_radius', None)
    center = (0.0, 0.0, 0.0)

    # Build rounded-corner specs if a corner radius is provided
    if corner_radius is not None and corner_radius > 0.0:
        wall_thickness = (outer_side - inner_side) / 2.0
        inner_rc = [(i, corner_radius) for i in range(4)]
        outer_rc = [(i, corner_radius + wall_thickness) for i in range(4)]
    else:
        inner_rc = None
        outer_rc = None

    # --- 1. Create bounding-box cell ---
    tmp_cell = RectCell(
        name=water_rod_model.rod_ID,
        height_x_width=(bb, bb),
        center=center,
    )

    # --- 2. Create inner boundary rectangles ---
    inner_rect = Rectangle(
        name=f"{water_rod_model.rod_ID}_inner",
        height=inner_side,
        width=inner_side,
        center=center,
        rounded_corners=inner_rc,
    )
    outer_rect = Rectangle(
        name=f"{water_rod_model.rod_ID}_outer",
        height=outer_side,
        width=outer_side,
        center=center,
        rounded_corners=outer_rc,
    )

    # --- 3. Partition the cell face with the two boundaries ---
    partitioned_face = make_partition(
        [tmp_cell.face],
        [inner_rect.face, outer_rect.face],
        shape_type=ShapeType.COMPOUND,
    )
    tmp_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face,
    )

    # --- 4. Base material and MACRO assignment (3 regions) ---
    tmp_cell.set_properties({
        PropertyType.MATERIAL: [
            water_rod_model.moderator_material_name,
            water_rod_model.cladding_material_name,
            water_rod_model.coolant_material_name,
        ],
        PropertyType.MACRO: [f"MACRO_{water_rod_model.rod_ID}"] * 3,
    })

    # --- 5. Optional Cartesian grid sub-meshing ---
    splits = None
    if calculation_step is not None:
        wr_cfg = calculation_step.get_water_rod_sectorization()
        if wr_cfg is not None:
            if wr_cfg.splits is not None:
                splits = wr_cfg.splits
                # Warn if circular-only keys are also populated
                if wr_cfg.sectors:
                    warnings.warn(
                        "Square water rod: 'sectors'/'angles' in the "
                        "water_rods config are ignored; only 'splits' "
                        "is used for square water rods.",
                        stacklevel=2,
                    )
            elif wr_cfg.sectors:
                # Warn if user supplied sectors for a square water rod
                warnings.warn(
                    "Square water rod: 'sectors'/'angles' sectorization "
                    "is not applicable to square water rods.  Use "
                    "'splits: [nx, ny]' instead.  No discretization "
                    "will be applied.",
                    stacklevel=2,
                )

    if splits is not None:
        nx, ny = splits
        # Build the grid of splitting faces over the bounding box
        bb_rect = Rectangle(
            name=f"{water_rod_model.rod_ID}_grid",
            height=bb,
            width=bb,
            center=center,
        )
        splitting_faces = make_grid_faces(bb_rect, nx, ny)

        # Re-partition the (already 3-region) cell face
        re_partitioned = make_partition(
            [tmp_cell.face],
            splitting_faces,
            shape_type=ShapeType.COMPOUND,
        )
        tmp_cell.update_geometry_from_face(
            GeometryType.TECHNOLOGICAL, re_partitioned,
        )

        # Reassign materials by geometric containment
        subfaces = tmp_cell.extract_subfaces()
        n_regions = len(subfaces)
        materials_list = [""] * n_regions
        macros_list = [f"MACRO_{water_rod_model.rod_ID}"] * n_regions

        for i, subface in enumerate(subfaces):
            pt = make_vertex_inside_face(subface)
            if is_point_inside_shape(pt, inner_rect.face):
                materials_list[i] = water_rod_model.moderator_material_name
            elif is_point_inside_shape(pt, outer_rect.face):
                materials_list[i] = water_rod_model.cladding_material_name
            else:
                materials_list[i] = water_rod_model.coolant_material_name

        tmp_cell.set_properties({
            PropertyType.MATERIAL: materials_list,
            PropertyType.MACRO: macros_list,
        })

        print(f"_build_square_water_rod_cell: sub-meshed "
              f"'{water_rod_model.rod_ID}' into {n_regions} sub-regions "
              f"(splits={splits}).")
    else:
        print(f"_build_square_water_rod_cell: built "
              f"'{water_rod_model.rod_ID}' with 3 base regions "
              f"(no grid sub-meshing).")

    return tmp_cell


def create_and_add_water_rods_to_lattice(lattice, assembly_model, translation=0.0, windmill=False, calculation_step=None):
    """
    Create water rod cells from the assembly model and add them to the lattice at their centers.

    Parameters:
    -----------
    lattice : Lattice
        The lattice to which water rod cells will be added
    assembly_model : CartesianAssemblyModel
        The assembly model containing the water rod geometry parameters
        (water_rod_type, water_rods list with center, radii, materials, etc.)
    translation : float
        Translation offset to apply to water rod center positions
    windmill : bool
        Whether to apply windmill sectorization to the water rod coolant region.
        Ignored if ``calculation_step`` is provided.
    calculation_step : CalculationStep or None
        Optional calculation step providing water-rod sectorization config.
        When provided, overrides the ``windmill`` parameter.
    """
    from ..DDModel.DragonModel import CircularWaterRodModel, SquareWaterRodModel

    if assembly_model.water_rod_type not in ("circular", "square"):
        raise ValueError(
            f"Unsupported water rod type: {assembly_model.water_rod_type}. "
            "Supported types are 'circular' and 'square'."
        )

    for water_rod_model in assembly_model.water_rods:
        if assembly_model.water_rod_type == "circular":
            tmp_cell = RectCell(
                name=water_rod_model.rod_ID,
                height_x_width=(
                    water_rod_model.bounding_box_side_length,
                    water_rod_model.bounding_box_side_length,
                ),
                center=(0.0, 0.0, 0.0),
            )
            tmp_cell.add_circle(water_rod_model.inner_radius)
            tmp_cell.add_circle(water_rod_model.outer_radius)
            tmp_cell.set_properties({
                PropertyType.MATERIAL: [
                    water_rod_model.moderator_material_name,
                    water_rod_model.cladding_material_name,
                    water_rod_model.coolant_material_name,
                ],
                PropertyType.MACRO: [f"MACRO_{water_rod_model.rod_ID}"] * 3,
            })
            # Apply sectorization: prefer calculation_step config, fall back to windmill flag
            if calculation_step is not None:
                wr_sectors = calculation_step.get_water_rod_sectorization()
                if wr_sectors is not None:
                    if wr_sectors.splits is not None:
                        import warnings
                        warnings.warn(
                            "Circular water rod: 'splits' in the "
                            "water_rods config is ignored; only "
                            "'sectors'/'angles' are used for circular "
                            "water rods.",
                            stacklevel=2,
                        )
                    tmp_cell.sectorize(wr_sectors.sectors, wr_sectors.angles, windmill=wr_sectors.windmill)
            elif windmill:
                tmp_cell.sectorize([1, 1, 8], [0, 0, 0], windmill=True)
        elif assembly_model.water_rod_type == "square":
            tmp_cell = _build_square_water_rod_cell(
                water_rod_model, calculation_step=calculation_step,
            )

        # water_rod_model.center is set by analyze_lattice_description and
        # expressed in lattice coordinates (pin-pitch units already resolved)
        cx, cy = water_rod_model.center
        lattice.add_cell(
            tmp_cell,
            (cx + translation, cy + translation, 0.0),
        )

    return lattice


def export_glow_geom(output_path, output_file_name, lattice, tracking_option, export_macro=False):
    """
    Export the geometry of the lattice to a TDT file for GLOW simulation.

    Parameters
    ----------
    output_path : str
        Path to save the exported TDT file
    output_file_name : str
        Name of the exported TDT file
    lattice : Lattice
        The lattice whose geometry is to be exported
    tracking_option : str
        Tracking option, either ``"TISO"`` or ``"TSPC"``
    """
    if export_macro:
        properties_to_export = [PropertyType.MATERIAL, PropertyType.MACRO]
        output_file_name = f"{output_file_name}_{tracking_option}_MACRO"
    else:
        properties_to_export = [PropertyType.MATERIAL]
        output_file_name = f"{output_file_name}_{tracking_option}"

    if tracking_option == "TISO":
        lattice.type_geo = LatticeGeometryType.ISOTROPIC
        analyse_and_generate_tdt(
            [lattice], f"{output_path}/{output_file_name}", TdtSetup(GeometryType.SECTORIZED,
                                             property_types=properties_to_export,
                                             type_geo=LatticeGeometryType.ISOTROPIC,
                                             symmetry_type=SymmetryType.FULL))
    elif tracking_option == "TSPC":
        lattice.type_geo = LatticeGeometryType.RECTANGLE_SYM    
        analyse_and_generate_tdt(
            [lattice], f"{output_path}/{output_file_name}", TdtSetup(GeometryType.SECTORIZED,
                                             property_types=properties_to_export,
                                             type_geo=LatticeGeometryType.RECTANGLE_SYM,
                                             symmetry_type=SymmetryType.FULL))


def build_assembly_box(assembly_model, center=None):
    """
    Build the 3-layer assembly box cell (intra-assembly coolant, channel
    box, inter-assembly moderator) from the assembly model dimensions.

    The resulting ``RectCell`` has three concentric regions whose materials
    are set to ``["COOLANT", "CHANNEL_BOX", "MODERATOR"]`` (innermost to
    outermost).

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model providing ``assembly_pitch``, ``gap_wide``,
        ``channel_box_thickness``, and ``corner_inner_radius_of_curvature``.
    center : tuple or None
        ``(x, y, z)`` centre of the box.  Defaults to
        ``(assembly_pitch / 2, assembly_pitch / 2, 0)``.

    Returns
    -------
    assembly_box_cell : RectCell
        The partitioned assembly box cell with 3 material regions.
    """
    ap = assembly_model.assembly_pitch
    gap = assembly_model.gap_wide
    cbt = assembly_model.channel_box_thickness
    corner_r_inner = assembly_model.corner_inner_radius_of_curvature

    if center is None:
        center = (ap / 2.0, ap / 2.0, 0.0)

    channel_box_inner_side = ap - 2.0 * cbt - 2.0 * gap
    channel_box_outer_side = ap - 2.0 * gap
    if corner_r_inner > 0.0:
        corner_r_outer = corner_r_inner + cbt
    else:
        corner_r_outer = 0.0
    print(f"Building assembly box with parameters:\n")
    print(f"corner_r_inner: {corner_r_inner}")
    print(f"corner_r_outer: {corner_r_outer}")
    print(f"channel_box_inner_side: {channel_box_inner_side}")
    print(f"channel_box_outer_side: {channel_box_outer_side}")

    if corner_r_inner > 0.0 and corner_r_outer > 0.0:
        print("Using rounded corners for assembly box boundaries.")
        rounded_corners_coolant = [
            (0, corner_r_inner),
            (1, corner_r_inner),
            (2, corner_r_inner),
            (3, corner_r_inner),
        ]
        rounded_corners_chanbox = [
            (0, corner_r_outer),
            (1, corner_r_outer),
            (2, corner_r_outer),
            (3, corner_r_outer),
        ]
    else:
        rounded_corners_coolant = None
        rounded_corners_chanbox = None
        
    # Inner coolant boundary (rounded-corner rectangle)
    coolant_rect = Rectangle(
        name="intra_assembly_coolant",
        height=channel_box_inner_side,
        width=channel_box_inner_side,
        center=center,
        rounded_corners=rounded_corners_coolant
    )
    # Channel box boundary (rounded-corner rectangle)
    channel_box_rect = Rectangle(
        name="channel_box",
        height=channel_box_outer_side,
        width=channel_box_outer_side,
        center=center,
        rounded_corners=rounded_corners_chanbox,
    )
    # Outer moderator cell
    assembly_box_cell = RectCell(
        name="assembly_box",
        height_x_width=(ap, ap),
        center=center,
    )
    # Partition the box with the two inner boundaries
    partitioned_face = make_partition(
        [assembly_box_cell.face],
        [channel_box_rect.face, coolant_rect.face],
        shape_type=ShapeType.COMPOUND,
    )
    assembly_box_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face
    )
    # Assign materials: innermost (coolant) → channel box → moderator
    assembly_box_cell.set_properties({
        PropertyType.MATERIAL: ["COOLANT", "CHANNEL_BOX", "MODERATOR"],
    })

    return assembly_box_cell


def _classify_point_to_macro(x, y, x0, y0, x1, y1, pin_pitch, n_cols, n_rows):
    """
    Classify an (x, y) coordinate into a MACRO name based on its position
    relative to the pin-lattice footprint.

    The pin lattice occupies the rectangle ``[x0, x1] × [y0, y1]``.
    Side strips are subdivided into per-pin-row/column regions.

    Parameters
    ----------
    x, y : float
        Coordinates of the point.
    x0, y0 : float
        Lower-left corner of the pin-lattice footprint.
    x1, y1 : float
        Upper-right corner of the pin-lattice footprint.
    pin_pitch : float
        Pin pitch (used to identify column/row index).
    n_cols, n_rows : int
        Number of pin columns and rows.

    Returns
    -------
    str
        MACRO name (e.g. ``"LEFT_3"``, ``"CORNER_BL"``, ``"BASE_CELL"``).
    """
    eps = 1e-6  # tolerance for boundary checks

    in_x_band = (x0 - eps) <= x <= (x1 + eps)
    in_y_band = (y0 - eps) <= y <= (y1 + eps)

    if in_x_band and in_y_band:
        # Inside the pin-lattice footprint (residual coolant region)
        return "BASE_CELL"

    if in_x_band:
        # Top or bottom strip — index by column
        col = int((x - x0) / pin_pitch)
        col = max(0, min(col, n_cols - 1))
        col_label = col + 1  # 1-based
        if y < y0:
            return f"BOT_{col_label}"
        else:
            return f"TOP_{col_label}"

    if in_y_band:
        # Left or right strip — index by row
        row = int((y - y0) / pin_pitch)
        row = max(0, min(row, n_rows - 1))
        row_label = row + 1  # 1-based
        if x < x0:
            return f"LEFT_{row_label}"
        else:
            return f"RIGHT_{row_label}"

    # Corner regions
    if x < x0 and y < y0:
        return "CORNER_BL"
    elif x >= x1 and y < y0:
        return "CORNER_BR"
    elif x < x0 and y >= y1:
        return "CORNER_TL"
    else:
        return "CORNER_TR"


def subdivide_box_into_macros(assembly_box_cell, assembly_model):
    """
    Partition the assembly box cell into per-pin-row/column MACRO regions
    required for the IC spatial method in DRAGON.

    The algorithm:

    1. Compute the pin-lattice footprint from the assembly dimensions.
    2. Create splitting rectangles for the 8 strips surrounding the
       lattice (4 sides + 4 corners), with side strips further divided
       into ``n_cols`` or ``n_rows`` sub-strips.
    3. Partition the box cell face with these splitting faces.
    4. Update the technological geometry (``update_geometry_from_face``),
       which automatically inherits MATERIAL properties from the
       pre-partition 3-region cell via spatial containment.
    5. Query each resulting sub-face's centroid position to assign a
       MACRO name programmatically (no trial-and-error ordering).

    MACRO naming convention (all 3 material layers in a strip share the
    same MACRO):

    - ``BOT_k``  / ``TOP_k``  — bottom / top strip at pin column *k*
    - ``LEFT_k`` / ``RIGHT_k`` — left / right strip at pin row *k*
    - ``CORNER_BL``, ``CORNER_BR``, ``CORNER_TL``, ``CORNER_TR``
    - ``BASE_CELL`` — residual intra-assembly coolant inside the
      pin-lattice footprint

    Parameters
    ----------
    assembly_box_cell : RectCell
        The 3-region assembly box cell (as returned by
        ``build_assembly_box``).
    assembly_model : CartesianAssemblyModel
        Assembly model providing dimensional information.

    Returns
    -------
    assembly_box_cell : RectCell
        The updated cell with MACRO and MATERIAL properties set on all
        sub-regions.
    """
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_rows = len(assembly_model.lattice_description)
    n_cols = len(assembly_model.lattice_description[0])

    lattice_pitch_x = n_cols * pin_pitch
    lattice_pitch_y = n_rows * pin_pitch

    # Pin-lattice footprint corners
    x0 = (ap - lattice_pitch_x) / 2.0
    y0 = (ap - lattice_pitch_y) / 2.0
    x1 = x0 + lattice_pitch_x
    y1 = y0 + lattice_pitch_y

    # ------------------------------------------------------------------
    # Build splitting rectangles (8 strips around the lattice, no centre)
    # ------------------------------------------------------------------
    rectangles_and_splits = [
        # (Rectangle, (nx, ny))
        # Bottom-left corner
        (Rectangle(height=y0, width=x0,
                   center=(x0 / 2.0, y0 / 2.0, 0.0)),
         (1, 1)),
        # Bottom-middle strip
        (Rectangle(height=y0, width=lattice_pitch_x,
                   center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
         (n_cols, 1)),
        # Bottom-right corner
        (Rectangle(height=y0, width=(ap - x1),
                   center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
         (1, 1)),
        # Middle-left strip
        (Rectangle(height=lattice_pitch_y, width=x0,
                   center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
         (1, n_rows)),
        # Middle-right strip
        (Rectangle(height=lattice_pitch_y, width=(ap - x1),
                   center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
         (1, n_rows)),
        # Top-left corner
        (Rectangle(height=(ap - y1), width=x0,
                   center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
         (1, 1)),
        # Top-middle strip
        (Rectangle(height=(ap - y1), width=lattice_pitch_x,
                   center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
         (n_cols, 1)),
        # Top-right corner
        (Rectangle(height=(ap - y1), width=(ap - x1),
                   center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
         (1, 1)),
    ]

    splitting_faces = []
    for rect, (nx, ny) in rectangles_and_splits:
        splitting_faces.extend(make_grid_faces(rect, nx, ny))

    # ------------------------------------------------------------------
    # Partition the box cell face
    # ------------------------------------------------------------------
    partitioned_face = make_partition(
        [assembly_box_cell.face],
        splitting_faces,
        shape_type=ShapeType.COMPOUND,
    )
    # update_geometry_from_face inherits MATERIAL from the 3-region cell
    assembly_box_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face
    )

    # ------------------------------------------------------------------
    # Build reference boundary faces for material classification.
    # We cannot rely on tech_geom_props inheritance from
    # update_geometry_from_face because deepcopy corrupts SALOME CORBA
    # proxy keys.  Instead, recreate the coolant / channel-box boundary
    # shapes (same parameters as build_assembly_box) and classify each
    # sub-face by geometric containment.
    # ------------------------------------------------------------------
    gap = assembly_model.gap_wide
    cbt = assembly_model.channel_box_thickness
    corner_r_inner = assembly_model.corner_inner_radius_of_curvature
    center = (ap / 2.0, ap / 2.0, 0.0)

    channel_box_inner_side = ap - 2.0 * cbt - 2.0 * gap
    channel_box_outer_side = ap - 2.0 * gap
    if corner_r_inner > 0.0:
        corner_r_outer = corner_r_inner + cbt
    else:
        corner_r_outer = 0.0
        
    if corner_r_inner > 0.0 and corner_r_outer > 0.0:
        print("Using rounded corners for assembly box boundaries.")
        coolant_corners = [
            (0, corner_r_inner),
            (1, corner_r_inner),
            (2, corner_r_inner),
            (3, corner_r_inner),
        ]
        chanbox_corners = [
            (0, corner_r_outer),
            (1, corner_r_outer),
            (2, corner_r_outer),
            (3, corner_r_outer),
        ]
    else:
        coolant_corners = None
        chanbox_corners = None
        
    print(f"Rebuilding reference boundaries for material classification with parameters:\n")
    print(f"corner_r_inner: {corner_r_inner}")
    print(f"corner_r_outer: {corner_r_outer}")
    print(f"channel_box_inner_side: {channel_box_inner_side}")
    print(f"channel_box_outer_side: {channel_box_outer_side}")
    print(f"center: {center}")
    print(f"coolant_corners: {coolant_corners}")
    print(f"chanbox_corners: {chanbox_corners}")

    coolant_boundary = Rectangle(
        name="_coolant_ref",
        height=channel_box_inner_side,
        width=channel_box_inner_side,
        center=center,
        rounded_corners=coolant_corners,
    )
    channel_box_boundary = Rectangle(
        name="_chanbox_ref",
        height=channel_box_outer_side,
        width=channel_box_outer_side,
        center=center,
        rounded_corners=chanbox_corners,
    )

    # ------------------------------------------------------------------
    # Query each sub-face to assign MATERIAL and MACRO
    # ------------------------------------------------------------------
    subfaces = assembly_box_cell.extract_subfaces()
    n_regions = len(subfaces)

    materials_list = [""] * n_regions
    macros_list = [""] * n_regions

    for i, subface in enumerate(subfaces):
        pt = make_vertex_inside_face(subface)
        px, py, _ = get_point_coordinates(pt)

        # Classify material by geometric containment against the
        # concentric boundary shapes (innermost test first)
        if is_point_inside_shape(pt, coolant_boundary.face):
            mat = "COOLANT"
        elif is_point_inside_shape(pt, channel_box_boundary.face):
            mat = "CHANNEL_BOX"
        else:
            mat = "MODERATOR"
        materials_list[i] = mat

        # Classify position → MACRO name
        macro = _classify_point_to_macro(
            px, py, x0, y0, x1, y1, pin_pitch, n_cols, n_rows
        )
        macros_list[i] = macro
    
    # handle corner cases : if the point is classified as "BASE_CELL" but material is "CHANNEL_BOX", it means it's a corner region that should be classified as a corner macro, not base cell
    for i in range(n_regions):
        if materials_list[i] == "CHANNEL_BOX" and macros_list[i] == "BASE_CELL":
            pt = make_vertex_inside_face(subfaces[i])
            px, py, _ = get_point_coordinates(pt)
            mid_x = ap / 2.0
            mid_y = ap / 2.0
            # test which corner it is closest to and map to corner pin MACRO{row}{col}
            row = 0 if py < mid_y else (n_rows - 1) 
            col = 0 if px < mid_x else (n_cols - 1)
            macros_list[i] = f"MACRO{row}{col}"

    # Sanity check
    unknowns = [i for i, m in enumerate(materials_list) if m == "UNKNOWN"]
    if unknowns:
        print(f"WARNING: {len(unknowns)} sub-faces could not be assigned a "
              f"material after box subdivision.")

    # Set both MATERIAL and MACRO properties in one call
    assembly_box_cell.set_properties({
        PropertyType.MATERIAL: materials_list,
        PropertyType.MACRO: macros_list,
    })

    print(f"subdivide_box_into_macros: split assembly box into "
          f"{n_regions} sub-regions with "
          f"{len(set(macros_list))} unique MACROs.")

    return assembly_box_cell


def discretize_box(assembly_box_cell, assembly_model, box_discretization_config):
    """
    Subdivide the assembly-box peripheral regions into a grid of
    sub-faces for MOC tracking.

    The pin-lattice footprint is computed from the assembly dimensions.
    The surrounding area is split into 8 rectangular strips (4 corners
    + 4 sides) and each strip is further gridded into sub-faces using
    ``make_grid_faces``.

    Unlike ``subdivide_box_into_macros`` (used for the IC method), this
    function does **not** assign MACRO properties by default.  Material
    properties are inherited from the 3-region box cell via
    ``update_geometry_from_face``.

    Parameters
    ----------
    assembly_box_cell : RectCell
        The 3-region assembly box cell (as returned by
        ``build_assembly_box``).
    assembly_model : CartesianAssemblyModel
        Assembly model providing dimensional information.
    box_discretization_config : BoxDiscretizationConfig
        Configuration specifying the grid splits for corner and side
        strips.

    Returns
    -------
    assembly_box_cell : RectCell
        The updated cell with a finer technological geometry.
    """
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_rows = len(assembly_model.lattice_description)
    n_cols = len(assembly_model.lattice_description[0])

    lattice_pitch_x = n_cols * pin_pitch
    lattice_pitch_y = n_rows * pin_pitch

    # Pin-lattice footprint corners
    x0 = (ap - lattice_pitch_x) / 2.0
    y0 = (ap - lattice_pitch_y) / 2.0
    x1 = x0 + lattice_pitch_x
    y1 = y0 + lattice_pitch_y

    # Resolve split counts (uses lattice dims as defaults for sides)
    corner, side_x, side_y = box_discretization_config.resolve_splits(
        n_cols, n_rows
    )

    # ------------------------------------------------------------------
    # 8 peripheral rectangles (same layout as subdivide_box_into_macros)
    # ------------------------------------------------------------------
    rectangles_and_splits = [
        # Bottom-left corner
        (Rectangle(height=y0, width=x0,
                   center=(x0 / 2.0, y0 / 2.0, 0.0)),
         corner),
        # Bottom-middle strip
        (Rectangle(height=y0, width=lattice_pitch_x,
                   center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
         side_x),
        # Bottom-right corner
        (Rectangle(height=y0, width=(ap - x1),
                   center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
         corner),
        # Middle-left strip
        (Rectangle(height=lattice_pitch_y, width=x0,
                   center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
         side_y),
        # Middle-right strip
        (Rectangle(height=lattice_pitch_y, width=(ap - x1),
                   center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
         side_y),
        # Top-left corner
        (Rectangle(height=(ap - y1), width=x0,
                   center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
         corner),
        # Top-middle strip
        (Rectangle(height=(ap - y1), width=lattice_pitch_x,
                   center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
         side_x),
        # Top-right corner
        (Rectangle(height=(ap - y1), width=(ap - x1),
                   center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
         corner),
    ]

    splitting_faces = []
    for rect, (nx, ny) in rectangles_and_splits:
        splitting_faces.extend(make_grid_faces(rect, nx, ny))

    # ------------------------------------------------------------------
    # Partition the box cell face
    # ------------------------------------------------------------------
    partitioned_face = make_partition(
        [assembly_box_cell.face],
        splitting_faces,
        shape_type=ShapeType.COMPOUND,
    )
    assembly_box_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face
    )

    n_subfaces = len(assembly_box_cell.extract_subfaces())
    print(f"discretize_box: split assembly box into "
          f"{n_subfaces} sub-regions (corner={corner}, "
          f"side_x={side_x}, side_y={side_y}).")

    return assembly_box_cell


def build_full_assembly_geometry(assembly_model, calculation_step,
                                 output_path, output_file_name,
                                 translation=None):
    """
    High-level function that builds a complete assembly geometry —
    including the assembly box with automatic MACRO subdivision for
    IC-method compatibility — and exports it to a TDT file.

    This function orchestrates the full pipeline:

    1. Apply radial discretization from the ``calculation_step``.
    2. Generate fuel ``RectCell`` objects (with sectorization).
    3. Build the 3-layer assembly box (coolant / channel box / moderator).
    4. If the step uses the IC spatial method with macro export, subdivide
       the assembly box into per-pin-row/column MACRO regions
       automatically.
    5. Assemble the ``Lattice``, add fuel cells, water rods, and box.
    6. Export the TDT file.

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model with lattice built and material mixtures numbered.
    calculation_step : CalculationStep
        The calculation step whose discretization config drives the geometry.
    output_path : str
        Directory to write the TDT file.
    output_file_name : str
        Base name for the output TDT file (tracking suffix is appended
        automatically by ``export_glow_geom``).
    translation : float or None
        Translation offset applied to cell positions in the lattice.
        If ``None``, computed automatically from assembly dimensions as
        ``gap_wide + channel_box_thickness + coolant_intra_assembly_width``.

    Returns
    -------
    lattice : Lattice
        The assembled glow ``Lattice`` object, already exported.
    assembly_box_cell : RectCell
        The assembly box cell (with MACRO properties if applicable).
    """
    from ..DDModel.DragonModel import CartesianAssemblyModel  # type check only

    # ----- Derived dimensions -----
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    gap = assembly_model.gap_wide
    cbt = assembly_model.channel_box_thickness
    n_cols = len(assembly_model.lattice_description[0])
    channel_box_inner_side = ap - 2.0 * cbt - 2.0 * gap
    coolant_intra_assembly_width = (channel_box_inner_side - n_cols * pin_pitch) / 2.0

    if translation is None:
        translation = gap + cbt + coolant_intra_assembly_width

    center = (ap / 2.0, ap / 2.0, 0.0)

    # Step 1: Apply radii from the calculation step
    calculation_step.apply_radii(assembly_model)

    # Step 2: Generate fuel cells with sectorization
    ordered_cells = generate_fuel_cells(
        assembly_model, calculation_step=calculation_step
    )

    # Step 3: Build assembly box
    assembly_box_cell = build_assembly_box(assembly_model, center=center)

    # Step 4: Subdivide assembly box if needed
    #   - IC + macros  → per-pin-row/column MACRO regions
    #   - MOC + box_discretization enabled → fine grid for MOC tracking
    if (calculation_step.spatial_method == "IC"
            and calculation_step.export_macros):
        assembly_box_cell = subdivide_box_into_macros(
            assembly_box_cell, assembly_model
        )
    elif (calculation_step.box_discretization is not None
              and calculation_step.box_discretization.enabled):
        assembly_box_cell = discretize_box(
            assembly_box_cell, assembly_model,
            calculation_step.box_discretization,
        )

    # Step 5: Build lattice and add cells
    lattice = Lattice(
        name=f"{assembly_model.name}_{calculation_step.name}",
        center=center,
    )

    lattice = add_cells_to_regular_lattice(
        lattice, ordered_cells, pin_pitch, translation=translation
    )

    if hasattr(assembly_model, "water_rods") and assembly_model.water_rods:
        lattice = create_and_add_water_rods_to_lattice(
            lattice, assembly_model,
            translation=translation,
            calculation_step=calculation_step,
        )

    lattice.lattice_box = assembly_box_cell

    # Step 6: Export TDT
    export_glow_geom(
        output_path,
        output_file_name,
        lattice,
        tracking_option=calculation_step.tracking,
        export_macro=calculation_step.export_macros,
    )

    return lattice, assembly_box_cell

