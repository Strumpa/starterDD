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


def _corner_transform(corner, x, y, ap):
    """
    Map canonical north-west coordinates ``(x, y)`` to the actual
    assembly corner.

    In the canonical (north-west) system the cross centre sits at the
    top-left corner ``(0, ap)``.  This helper mirrors the coordinates
    for the other three corners.

    Parameters
    ----------
    corner : str
        ``"north-west"``, ``"north-east"``, ``"south-west"``, or
        ``"south-east"``.
    x, y : float
        Coordinates in the canonical north-west system.
    ap : float
        Assembly pitch.

    Returns
    -------
    (float, float)
        Transformed ``(x, y)``.
    """
    if corner == "north-west":
        return (x, y)
    elif corner == "north-east":
        return (ap - x, y)
    elif corner == "south-west":
        return (x, ap - y)
    elif corner == "south-east":
        return (ap - x, ap - y)
    else:
        raise ValueError(f"Unknown corner '{corner}'.")


def _remap_rounded_corner_indices(corner_indices, cross_corner):
    """
    Remap glow rounded-corner indices for a wing-tip rectangle after
    applying the corner transform.

    In the canonical north-west system the wing tips have a rounded
    corner at glow index 1 (= bottom-right of the rectangle in glow
    convention).  When the cross is placed at a different assembly
    corner the rectangle is mirrored and the rounded-corner index
    must change accordingly.

    Parameters
    ----------
    corner_indices : list of (int, float)
        List of ``(glow_corner_index, radius)`` pairs in the canonical
        system (north-west).
    cross_corner : str
        The actual assembly corner.

    Returns
    -------
    list of (int, float)
        Remapped corner/radius pairs.
    """
    # Mapping: NW is identity.  Mirror in x flips left↔right (0↔1, 3↔2).
    # Mirror in y flips top↔bottom (0↔3, 1↔2).
    _mirror_x = {0: 1, 1: 0, 2: 3, 3: 2}
    _mirror_y = {0: 3, 1: 2, 2: 1, 3: 0}

    result = list(corner_indices)
    if cross_corner in ("north-east", "south-east"):
        result = [(_mirror_x[idx], r) for idx, r in result]
    if cross_corner in ("south-west", "south-east"):
        result = [(_mirror_y[idx], r) for idx, r in result]
    return result


def _build_control_cross_shapes(ctrl, ap):
    """
    Build all glow geometry shapes for a control cross and return them
    in a structured dict.

    Shapes are built at the correct position for ``ctrl.center`` using
    ``_corner_transform``.

    Parameters
    ----------
    ctrl : ControlCrossModel
        The control cross model with all geometric dimensions.
    ap : float
        Assembly pitch.

    Returns
    -------
    dict with keys:

    - ``"sheath_rectangles"`` : list of ``Rectangle``
        Outer sheath/structural boundary shapes (5 rectangles:
        quarter centre, right central structure half, bottom central
        structure half, right wing half, bottom wing half).
    - ``"inner_sheath_rectangles"`` : list of ``Rectangle``
        Inner sheath boundary shapes for each wing (2 rectangles).
    - ``"absorber_tubes"`` : list of ``RectCell``
        Absorber tube cells (2 × ``number_tubes_per_wing``), each with
        inner and outer circle radii.
    - ``"splitting_rectangles"`` : list of ``Rectangle``
        Rectangles that split each wing bounding box at the boundary
        of the last absorber tube, so the tip region is separated.
    - ``"wing_footprints"`` : dict
        ``{"wing_1": (x_min, y_min, x_max, y_max),
          "wing_2": (x_min, y_min, x_max, y_max),
          "center": (x_min, y_min, x_max, y_max)}``
        Axis-aligned bounding boxes for the two wings and central
        structure (in assembly coordinates), used for MACRO
        classification.
    """
    corner = ctrl.center
    bt = ctrl.blade_thickness
    bhs = ctrl.blade_half_span
    cshs = ctrl.central_structure_half_span
    st = ctrl.sheath_thickness
    tr = ctrl.tip_radius
    n_tubes = ctrl.number_tubes_per_wing
    r_inner = ctrl.absorber_tube_inner_radius
    r_outer = ctrl.absorber_tube_outer_radius
    inner_w = ctrl.inner_sheath_width  # = bt - 2*st
    delta = ctrl.tube_spacing
    first_offset = ctrl.first_tube_offset

    def ct(x, y):
        """Shorthand for corner transform."""
        return _corner_transform(corner, x, y, ap)

    # ------------------------------------------------------------------
    # 1. Sheath / structural rectangles (canonical NW coords)
    # ------------------------------------------------------------------
    sheath_rects = []

    # Quarter centre
    cx, cy = ct(bt / 4.0, ap - bt / 4.0)
    sheath_rects.append(Rectangle(
        name="CTRL_CROSS_CENTER_QTR",
        height=bt / 2.0, width=bt / 2.0,
        center=(cx, cy, 0.0),
    ))

    # Right central structure half
    w_rcs = cshs - bt / 2.0
    cx, cy = ct(w_rcs / 2.0 + bt / 2.0, ap - bt / 4.0)
    sheath_rects.append(Rectangle(
        name="CTRL_CROSS_RIGHT_CS_HALF",
        height=bt / 2.0, width=w_rcs,
        center=(cx, cy, 0.0),
    ))

    # Bottom central structure half
    h_bcs = cshs - bt / 2.0
    cx, cy = ct(bt / 4.0, ap - h_bcs / 2.0 - bt / 2.0)
    sheath_rects.append(Rectangle(
        name="CTRL_CROSS_BOT_CS_HALF",
        height=h_bcs, width=bt / 2.0,
        center=(cx, cy, 0.0),
    ))

    # Right wing half (horizontal arm) — with tip radius
    wing_len = bhs - cshs
    cx, cy = ct(wing_len / 2.0 + cshs, ap)
    # In canonical NW, tip rounded corner is at index 1
    if tr > 0.0:
        rc_wing_h = _remap_rounded_corner_indices([(1, tr)], corner)
    else:
        rc_wing_h = None
    sheath_rects.append(Rectangle(
        name="CTRL_CROSS_WING_H",
        height=bt, width=wing_len,
        center=(cx, cy, 0.0),
        rounded_corners=rc_wing_h,
    ))

    # Bottom wing half (vertical arm) — with tip radius
    cx, cy = ct(0.0, ap - wing_len / 2.0 - cshs)
    if tr > 0.0:
        rc_wing_v = _remap_rounded_corner_indices([(1, tr)], corner)
    else:
        rc_wing_v = None
    sheath_rects.append(Rectangle(
        name="CTRL_CROSS_WING_V",
        height=wing_len, width=bt,
        center=(cx, cy, 0.0),
        rounded_corners=rc_wing_v,
    ))

    # ------------------------------------------------------------------
    # 2. Inner sheath rectangles (absorber cavity boundary)
    #    Skipped for solid crosses (st == 0): there is no hollow
    #    cavity inside the blade — the blade material fills
    #    everything between absorber rods.
    # ------------------------------------------------------------------
    inner_rects = []

    if st > 0:
        inner_wing_w = wing_len - st  # inner wing length
        cx, cy = ct(inner_wing_w / 2.0 + cshs, ap)
        if tr > 0.0:
            rc_inner_h = _remap_rounded_corner_indices([(1, tr - st)], corner)
        else:
            rc_inner_h = None
        inner_rects.append(Rectangle(
            name="CTRL_CROSS_WING_H_INNER",
            height=inner_w, width=inner_wing_w,
            center=(cx, cy, 0.0),
            rounded_corners=rc_inner_h,
        ))

        cx, cy = ct(0.0, ap - inner_wing_w / 2.0 - cshs)
        if tr > 0.0:
            rc_inner_v = _remap_rounded_corner_indices([(1, tr - st)], corner)
        else:
            rc_inner_v = None
        inner_rects.append(Rectangle(
            name="CTRL_CROSS_WING_V_INNER",
            height=inner_wing_w, width=inner_w,
            center=(cx, cy, 0.0),
            rounded_corners=rc_inner_v,
        ))

    # ------------------------------------------------------------------
    # 3. Absorber tubes
    # ------------------------------------------------------------------
    tubes = []
    for i in range(n_tubes):
        offset = first_offset + i * delta

        # Horizontal wing tube (along x)
        tx, ty = ct(offset, ap)
        tube_h = RectCell(
            name=f"CTRL_TUBE_H_{i}",
            height_x_width=(inner_w, delta),
            center=(tx, ty, 0.0),
        )
        if r_inner < r_outer:
            # Hollow tube (GE-14 style): inner absorber + outer cladding
            tube_h.add_circle(r_inner)
            tube_h.add_circle(r_outer)
        else:
            # Solid rod (AT10 style): single absorber circle
            tube_h.add_circle(r_outer)
        tubes.append(tube_h)

        # Vertical wing tube (along y)
        tx, ty = ct(0.0, ap - offset)
        tube_v = RectCell(
            name=f"CTRL_TUBE_V_{i}",
            height_x_width=(delta, inner_w),
            center=(tx, ty, 0.0),
        )
        if r_inner < r_outer:
            tube_v.add_circle(r_inner)
            tube_v.add_circle(r_outer)
        else:
            tube_v.add_circle(r_outer)
        tubes.append(tube_v)

    # ------------------------------------------------------------------
    # 4. Splitting rectangles at last-tube boundary
    # ------------------------------------------------------------------
    last_tube_offset = first_offset + (n_tubes - 1) * delta
    last_boundary = last_tube_offset + delta / 2.0

    # Horizontal wing split: full-width rectangle from x=0 to
    # blade_half_span, height = bt/2, centered on the top edge.
    # Split at x = last_boundary.
    split_rects = []

    # Left part of north arm (covers tubes region)
    lw = last_boundary
    if lw > 1e-6:
        cx, cy = ct(lw / 2.0, ap - bt / 4.0)
        split_rects.append(Rectangle(
            name="CTRL_SPLIT_H_LEFT",
            height=bt / 2.0, width=lw,
            center=(cx, cy, 0.0),
        ))

    # Right part of north arm (tip region beyond last tube)
    rw = bhs - last_boundary
    cx, cy = ct(last_boundary + rw / 2.0, ap - bt / 4.0)
    print(f"building rect at cx={cx}, cy={cy} with rw={rw}, bt={bt}")
    if rw>1e-6:
        split_rects.append(Rectangle(
            name="CTRL_SPLIT_H_RIGHT",
            height=bt / 2.0, width=rw,
            center=(cx, cy, 0.0),
        ))

    # Bottom part of west arm (covers tubes region)
    bh = last_boundary
    cx, cy = ct(bt / 4.0, ap - lw / 2.0)
    print(f"building rect at cx={cx}, cy={cy} with bh={bh}, bt={bt}")
    if bh>1e-6:
        split_rects.append(Rectangle(
            name="CTRL_SPLIT_V_BOT",
            height=bh, width=bt / 2.0,
            center=(cx, cy, 0.0),
        ))

    # Top part of west arm (tip region beyond last tube)
    th = bhs - last_boundary
    cx, cy = ct(bt / 4.0, ap - last_boundary - th / 2.0)
    if th>1e-6:
        split_rects.append(Rectangle(
            name="CTRL_SPLIT_V_TOP",
            height=th, width=bt / 2.0,
            center=(cx, cy, 0.0),
        ))

    # ------------------------------------------------------------------
    # 5. Wing footprints (tight AABBs in assembly coordinates)
    # ------------------------------------------------------------------
    # Build one tight AABB per sheath structural region, clamped to the
    # assembly domain [0, ap]^2, so that the L-shaped central structure
    # does not produce an over-sized bounding box that would swallow
    # adjacent gap regions.
    def _tight_aabb(corners_canon):
        pts = [ct(cx_, cy_) for cx_, cy_ in corners_canon]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        return (max(0.0, min(xs)), max(0.0, min(ys)),
                min(ap, max(xs)), min(ap, max(ys)))

    # Each MACRO region must be convex.  The central structure is an
    # L-shape, so we decompose it: the quarter-centre square stays as
    # CROSS_CENTER; the right CS half is merged into CROSS_WING_1
    # (together they form the convex rectangle [bt/2, bhs] × blade-y);
    # the bottom CS half is merged into CROSS_WING_2.
    wing_footprints = [
        # Quarter centre — convex square
        (_tight_aabb([(0, ap - bt / 2.0), (bt / 2.0, ap - bt / 2.0),
                       (bt / 2.0, ap), (0, ap)]),
         "CROSS_CENTER"),
        # Horizontal wing + right CS half — convex rectangle
        (_tight_aabb([(bt / 2.0, ap - bt / 2.0), (bhs, ap - bt / 2.0),
                       (bhs, ap + bt / 2.0), (bt / 2.0, ap + bt / 2.0)]),
         "CROSS_WING_1"),
        # Vertical wing + bottom CS half — convex rectangle
        (_tight_aabb([(-bt / 2.0, ap - bhs), (bt / 2.0, ap - bhs),
                       (bt / 2.0, ap - bt / 2.0), (-bt / 2.0, ap - bt / 2.0)]),
         "CROSS_WING_2"),
    ]

    print(f"[control cross] Built {len(sheath_rects)} sheath rects, "
          f"{len(inner_rects)} inner rects, {len(tubes)} tubes, "
          f"{len(split_rects)} split rects for corner '{corner}'.")

    return {
        "sheath_rectangles": sheath_rects,
        "inner_sheath_rectangles": inner_rects,
        "absorber_tubes": tubes,
        "splitting_rectangles": split_rects,
        "wing_footprints": wing_footprints,
    }


def build_assembly_box(assembly_model, center=None):
    """
    Build the assembly box cell from the assembly model dimensions.

    Without a control cross the result is a 3-region cell (intra-assembly
    coolant / channel box / inter-assembly moderator).

    When ``assembly_model.has_control_cross`` is ``True`` the control
    cross shapes (sheath, inner cavity, absorber tubes) are also
    included in the partition and materials are assigned by geometric
    containment against all reference shapes.

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model providing ``assembly_pitch``, ``gap_wide``,
        ``channel_box_thickness``, ``corner_inner_radius_of_curvature``,
        and optionally ``control_cross``.
    center : tuple or None
        ``(x, y, z)`` centre of the box.  Defaults to
        ``(assembly_pitch / 2, assembly_pitch / 2, 0)``.

    Returns
    -------
    assembly_box_cell : RectCell
        The partitioned assembly box cell with material properties set
        on every sub-face.
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
        rounded_corners=rounded_corners_coolant,
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

    # ------------------------------------------------------------------
    # Build partition shapes list
    # ------------------------------------------------------------------
    partition_shapes = [channel_box_rect.face, coolant_rect.face]

    ctrl_shapes = None
    has_cross = getattr(assembly_model, "has_control_cross", False)
    if has_cross:
        ctrl = assembly_model.control_cross
        ctrl_shapes = _build_control_cross_shapes(ctrl, ap)
        for r in ctrl_shapes["sheath_rectangles"]:
            partition_shapes.append(r.face)
        for r in ctrl_shapes["inner_sheath_rectangles"]:
            partition_shapes.append(r.face)
        for t in ctrl_shapes["absorber_tubes"]:
            partition_shapes.append(t.face)
        for r in ctrl_shapes["splitting_rectangles"]:
            partition_shapes.append(r.face)

    # ------------------------------------------------------------------
    # Single partition call
    # ------------------------------------------------------------------
    partitioned_face = make_partition(
        [assembly_box_cell.face],
        partition_shapes,
        shape_type=ShapeType.COMPOUND,
    )
    assembly_box_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face,
    )

    # ------------------------------------------------------------------
    # Assign materials
    # ------------------------------------------------------------------
    if not has_cross:
        # Simple 3-region case (backward compatible)
        assembly_box_cell.set_properties({
            PropertyType.MATERIAL: ["COOLANT", "CHANNEL_BOX", "MODERATOR"],
        })
    else:
        # Classify every sub-face by geometric containment
        _assign_materials_with_control_cross(
            assembly_box_cell,
            coolant_rect,
            channel_box_rect,
            ctrl_shapes,
            ctrl,
        )
        # Attach cross metadata for downstream use by
        # subdivide_box_into_macros
        assembly_box_cell._ctrl_cross_shapes = ctrl_shapes

    return assembly_box_cell


def _assign_materials_with_control_cross(
    assembly_box_cell, coolant_rect, channel_box_rect, ctrl_shapes, ctrl
):
    """
    Assign MATERIAL properties to every sub-face of the assembly box
    cell using geometric containment against reference shapes.

    Priority order (innermost first):

    1. Absorber tube inner circle → absorber material
    2. Absorber tube outer circle → sheath material (tube cladding)
    3. Inside inner-sheath rectangle but outside tubes → MODERATOR
       (inter-tube gap)
    4. Inside outer-sheath / structure rectangle → sheath material
    5. Inside coolant boundary → COOLANT
    6. Inside channel-box boundary → CHANNEL_BOX
    7. Otherwise → MODERATOR
    """
    subfaces = assembly_box_cell.extract_subfaces()
    n = len(subfaces)
    materials = [""] * n

    # Collect reference faces for containment tests
    tube_cells = ctrl_shapes["absorber_tubes"]
    inner_rects = ctrl_shapes["inner_sheath_rectangles"]
    sheath_rects = ctrl_shapes["sheath_rectangles"]

    solid = ctrl.is_solid  # AT10-style (no sheath cavity)

    # Pre-extract absorber tube circle faces for containment.
    # Hollow tubes (GE-14): 2 circles per tube (inner absorber, outer cladding).
    # Solid rods  (AT10):   1 circle per tube (absorber only).
    if solid:
        tube_circle_faces = [t.inner_circles[0].face for t in tube_cells]
    else:
        tube_inner_faces = []
        tube_outer_faces = []
        for t in tube_cells:
            # inner_circles[0] is the smaller (absorber inner) circle
            # inner_circles[1] is the larger (absorber outer) circle
            tube_inner_faces.append(t.inner_circles[0].face)
            tube_outer_faces.append(t.inner_circles[1].face)

    inner_rect_faces = [r.face for r in inner_rects]
    sheath_rect_faces = [r.face for r in sheath_rects]

    absorber_mat = ctrl.absorber_material
    sheath_mat = ctrl.sheath_material

    for i, subface in enumerate(subfaces):
        pt = make_vertex_inside_face(subface)

        # 1–2. Absorber tubes
        found_tube = False
        if solid:
            # Solid rod: single circle → absorber material
            for j, cf in enumerate(tube_circle_faces):
                if is_point_inside_shape(pt, cf):
                    materials[i] = absorber_mat
                    found_tube = True
                    break
        else:
            # Hollow tube: inner circle → absorber, outer annulus → sheath
            for j in range(len(tube_cells)):
                if is_point_inside_shape(pt, tube_inner_faces[j]):
                    materials[i] = absorber_mat
                    found_tube = True
                    break
                if is_point_inside_shape(pt, tube_outer_faces[j]):
                    materials[i] = sheath_mat
                    found_tube = True
                    break
        if found_tube:
            continue

        # 3. Inner sheath cavity (inter-tube moderator)
        #    Only applies to sheathed crosses — for solid crosses the
        #    inner_rects list is empty and this block is skipped.
        in_inner = False
        for irf in inner_rect_faces:
            if is_point_inside_shape(pt, irf):
                materials[i] = "MODERATOR"
                in_inner = True
                break
        if in_inner:
            continue

        # 4. Outer sheath / structural rectangles
        in_sheath = False
        for srf in sheath_rect_faces:
            if is_point_inside_shape(pt, srf):
                materials[i] = sheath_mat
                in_sheath = True
                break
        if in_sheath:
            continue

        # 5–7. Standard box regions
        if is_point_inside_shape(pt, coolant_rect.face):
            materials[i] = "COOLANT"
        elif is_point_inside_shape(pt, channel_box_rect.face):
            materials[i] = "CHANNEL_BOX"
        else:
            materials[i] = "MODERATOR"

    assembly_box_cell.set_properties({
        PropertyType.MATERIAL: materials,
    })

    # Summary
    from collections import Counter
    counts = Counter(materials)
    print(f"[control cross] Assigned materials to {n} sub-faces: "
          f"{dict(counts)}")


def _classify_point_to_macro(x, y, x0, y0, x1, y1, pin_pitch, n_cols, n_rows,
                             wing_footprints=None):
    """
    Classify an (x, y) coordinate into a MACRO name based on its position
    relative to the pin-lattice footprint and, optionally, the control
    cross wing footprints.

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
    wing_footprints : dict or None
        If provided, a dict with keys ``"wing_1"``, ``"wing_2"``,
        ``"center"`` each mapping to an ``(x_min, y_min, x_max, y_max)``
        axis-aligned bounding box.  Points inside these boxes are
        classified as ``"CROSS_WING_1"``, ``"CROSS_WING_2"``, or
        ``"CROSS_CENTER"`` instead of the standard MACRO names.

    Returns
    -------
    str
        MACRO name (e.g. ``"LEFT_3"``, ``"CORNER_BL"``, ``"BASE_CELL"``,
        ``"CROSS_WING_1"``).
    """
    eps = 1e-6  # tolerance for boundary checks

    # --- Control cross test (highest priority) ---
    if wing_footprints is not None:
        for (xmin, ymin, xmax, ymax), label in wing_footprints:
            if (xmin - eps) <= x <= (xmax + eps) and (ymin - eps) <= y <= (ymax + eps):
                return label

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


def _build_cross_aware_splitting_rects(
    ap, x0, y0, x1, y1, lattice_pitch_x, lattice_pitch_y,
    n_cols, n_rows, cross_corner, ctrl,
):
    """
    Build the MACRO-splitting rectangles for the 8 peripheral strips
    around the pin lattice, trimming strips adjacent to control cross
    wings so that they stop at the wing boundary.

    Parameters
    ----------
    ap : float
        Assembly pitch.
    x0, y0, x1, y1 : float
        Pin-lattice footprint corners.
    lattice_pitch_x, lattice_pitch_y : float
        Pin-lattice extents.
    n_cols, n_rows : int
        Number of pin columns / rows.
    cross_corner : str
        Cross centre corner (``"north-west"``, etc.).
    ctrl : ControlCrossModel
        Control cross geometry.

    Returns
    -------
    list of (Rectangle, (nx, ny))
        Splitting rectangles with their grid split counts.
    """
    bt = ctrl.blade_thickness
    bhs = ctrl.blade_half_span
    bt2 = bt / 2.0  # half thickness

    eps = 1e-8

    # ------------------------------------------------------------------
    # Compute blade and arm extents in assembly coordinates
    # ------------------------------------------------------------------
    # Vertical blade half-thickness: x range
    xv0, _ = _corner_transform(cross_corner, 0.0, 0.0, ap)
    xv1, _ = _corner_transform(cross_corner, bt2, 0.0, ap)
    x_blade_lo = min(xv0, xv1)
    x_blade_hi = max(xv0, xv1)

    # Horizontal blade half-thickness: y range
    _, yh0 = _corner_transform(cross_corner, 0.0, ap - bt2, ap)
    _, yh1 = _corner_transform(cross_corner, 0.0, ap, ap)
    y_blade_lo = min(yh0, yh1)
    y_blade_hi = max(yh0, yh1)

    # Vertical arm full y extent (from wing tip to centre)
    _, yv_tip = _corner_transform(cross_corner, 0.0, ap - bhs, ap)
    _, yv_ctr = _corner_transform(cross_corner, 0.0, ap, ap)
    y_arm_lo = min(yv_tip, yv_ctr)
    y_arm_hi = max(yv_tip, yv_ctr)

    # Horizontal arm full x extent (from centre to wing tip)
    xh_ctr, _ = _corner_transform(cross_corner, 0.0, 0.0, ap)
    xh_tip, _ = _corner_transform(cross_corner, bhs, 0.0, ap)
    x_arm_lo = min(xh_ctr, xh_tip)
    x_arm_hi = max(xh_ctr, xh_tip)

    # Determine which sides are affected by the cross
    cross_on_left = cross_corner in ("north-west", "south-west")
    cross_on_right = cross_corner in ("north-east", "south-east")
    cross_on_top = cross_corner in ("north-west", "north-east")
    cross_on_bottom = cross_corner in ("south-west", "south-east")

    rects = []

    # ==================================================================
    # CORNERS -- the cross-affected corner gets a gap-only rectangle;
    # the other three corners are unchanged.
    # ==================================================================

    # ---- Bottom-left corner ----
    if cross_corner == "south-west":
        gw = x0 - x_blade_hi
        gh = y0 - y_blade_hi
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x_blade_hi + gw / 2.0,
                                  y_blade_hi + gh / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=x0,
                      center=(x0 / 2.0, y0 / 2.0, 0.0)),
            (1, 1),
        ))

    # ---- Bottom-right corner ----
    if cross_corner == "south-east":
        gw = x_blade_lo - x1
        gh = y0 - y_blade_hi
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x1 + gw / 2.0,
                                  y_blade_hi + gh / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=(ap - x1),
                      center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
            (1, 1),
        ))

    # ---- Top-left corner ----
    if cross_corner == "north-west":
        gw = x0 - x_blade_hi
        gh = y_blade_lo - y1
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x_blade_hi + gw / 2.0,
                                  y1 + gh / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=x0,
                      center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
            (1, 1),
        ))

    # ---- Top-right corner ----
    if cross_corner == "north-east":
        gw = x_blade_lo - x1
        gh = y_blade_lo - y1
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x1 + gw / 2.0,
                                  y1 + gh / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=(ap - x1),
                      center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
            (1, 1),
        ))

    # ==================================================================
    # SIDE STRIPS -- cross-affected sides get a narrow gap column/row
    # (with n_rows/n_cols subdivisions aligned to the lattice pitch)
    # plus optional stubs for blade-width regions free of cross
    # structure (e.g. below the wing tip).
    # ==================================================================

    # ---- Bottom-middle strip ----
    if cross_on_bottom:
        # Gap row below the blade, at full lattice x-width
        gap_h = y0 - y_blade_hi
        if gap_h > eps:
            rects.append((
                Rectangle(height=gap_h, width=lattice_pitch_x,
                          center=((x0 + x1) / 2.0,
                                  y_blade_hi + gap_h / 2.0, 0.0)),
                (n_cols, 1),
            ))
        # Stubs in the blade region [0, y_blade_hi]
        blade_h = y_blade_hi
        arm_lo_c = max(x_arm_lo, x0)
        arm_hi_c = min(x_arm_hi, x1)
        stub_w = arm_lo_c - x0
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(x0 + stub_w / 2.0,
                                  blade_h / 2.0, 0.0)),
                (1, 1),
            ))
        stub_w = x1 - arm_hi_c
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(arm_hi_c + stub_w / 2.0,
                                  blade_h / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=lattice_pitch_x,
                      center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
            (n_cols, 1),
        ))

    # ---- Top-middle strip ----
    if cross_on_top:
        # Gap row above the lattice, below the blade
        gap_h = y_blade_lo - y1
        if gap_h > eps:
            rects.append((
                Rectangle(height=gap_h, width=lattice_pitch_x,
                          center=((x0 + x1) / 2.0,
                                  y1 + gap_h / 2.0, 0.0)),
                (n_cols, 1),
            ))
        # Stubs in the blade region [y_blade_lo, ap]
        blade_h = ap - y_blade_lo
        arm_lo_c = max(x_arm_lo, x0)
        arm_hi_c = min(x_arm_hi, x1)
        stub_w = arm_lo_c - x0
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(x0 + stub_w / 2.0,
                                  y_blade_lo + blade_h / 2.0, 0.0)),
                (1, 1),
            ))
        stub_w = x1 - arm_hi_c
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(arm_hi_c + stub_w / 2.0,
                                  y_blade_lo + blade_h / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=lattice_pitch_x,
                      center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
            (n_cols, 1),
        ))

    # ---- Middle-left strip ----
    if cross_on_left:
        # Gap column: from blade edge to lattice edge
        gap_w = x0 - x_blade_hi
        if gap_w > eps:
            rects.append((
                Rectangle(height=lattice_pitch_y, width=gap_w,
                          center=(x_blade_hi + gap_w / 2.0,
                                  (y0 + y1) / 2.0, 0.0)),
                (1, n_rows),
            ))
        # Stubs in the blade x-extent [x_blade_lo, x_blade_hi]
        blade_w = x_blade_hi - x_blade_lo
        arm_lo_c = max(y_arm_lo, y0)
        arm_hi_c = min(y_arm_hi, y1)
        stub_h = arm_lo_c - y0
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  y0 + stub_h / 2.0, 0.0)),
                (1, 1),
            ))
        stub_h = y1 - arm_hi_c
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  arm_hi_c + stub_h / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=lattice_pitch_y, width=x0,
                      center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
            (1, n_rows),
        ))

    # ---- Middle-right strip ----
    if cross_on_right:
        # Gap column: from lattice edge to blade edge
        gap_w = x_blade_lo - x1
        if gap_w > eps:
            rects.append((
                Rectangle(height=lattice_pitch_y, width=gap_w,
                          center=(x1 + gap_w / 2.0,
                                  (y0 + y1) / 2.0, 0.0)),
                (1, n_rows),
            ))
        # Stubs in the blade x-extent [x_blade_lo, x_blade_hi]
        blade_w = x_blade_hi - x_blade_lo
        arm_lo_c = max(y_arm_lo, y0)
        arm_hi_c = min(y_arm_hi, y1)
        stub_h = arm_lo_c - y0
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  y0 + stub_h / 2.0, 0.0)),
                (1, 1),
            ))
        stub_h = y1 - arm_hi_c
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  arm_hi_c + stub_h / 2.0, 0.0)),
                (1, 1),
            ))
    else:
        rects.append((
            Rectangle(height=lattice_pitch_y, width=(ap - x1),
                      center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
            (1, n_rows),
        ))

    return rects


def subdivide_box_into_macros(assembly_box_cell, assembly_model):
    """
    Partition the assembly box cell into per-pin-row/column MACRO regions
    required for the IC spatial method in DRAGON.

    The algorithm:

    1. Compute the pin-lattice footprint from the assembly dimensions.
    2. Create splitting rectangles for the 8 strips surrounding the
       lattice (4 sides + 4 corners), with side strips further divided
       into ``n_cols`` or ``n_rows`` sub-strips.  When a control cross
       is present, strips adjacent to the cross wings are trimmed so
       they do not extend into the wing footprint.
    3. Partition the box cell face with these splitting faces.
    4. Rebuild material + MACRO properties for every resulting sub-face
       by geometric containment (materials) and coordinate classification
       (MACROs).

    MACRO naming convention:

    - ``BOT_k``  / ``TOP_k``  — bottom / top strip at pin column *k*
    - ``LEFT_k`` / ``RIGHT_k`` — left / right strip at pin row *k*
    - ``CORNER_BL``, ``CORNER_BR``, ``CORNER_TL``, ``CORNER_TR``
    - ``BASE_CELL`` — residual intra-assembly coolant inside the
      pin-lattice footprint
    - ``CROSS_WING_1`` / ``CROSS_WING_2`` / ``CROSS_CENTER`` — control
      cross regions (when present)

    Parameters
    ----------
    assembly_box_cell : RectCell
        The assembly box cell (as returned by ``build_assembly_box``),
        already partitioned and with MATERIAL properties assigned.
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
    # Build splitting rectangles
    # ------------------------------------------------------------------
    has_cross = getattr(assembly_model, "has_control_cross", False)
    wing_footprints = None

    if has_cross:
        ctrl = assembly_model.control_cross
        rectangles_and_splits = _build_cross_aware_splitting_rects(
            ap, x0, y0, x1, y1,
            lattice_pitch_x, lattice_pitch_y,
            n_cols, n_rows,
            ctrl.center, ctrl,
        )
        # Recover wing footprints from build_assembly_box
        ctrl_shapes = getattr(assembly_box_cell, "_ctrl_cross_shapes", None)
        if ctrl_shapes is not None:
            wing_footprints = list(ctrl_shapes["wing_footprints"])

        # Compute stub footprints — moderator rectangles in the blade
        # region but outside the arm extent.  These must get their own
        # MACRO so they are not merged with the gap column/row MACROs.
        bt = ctrl.blade_thickness
        bhs = ctrl.blade_half_span
        bt2 = bt / 2.0
        corner = ctrl.center

        # Blade edge positions (secondary axis)
        xv0, _ = _corner_transform(corner, 0.0, 0.0, ap)
        xv1, _ = _corner_transform(corner, bt2, 0.0, ap)
        x_blade_lo = min(xv0, xv1)
        x_blade_hi = max(xv0, xv1)
        _, yh0 = _corner_transform(corner, 0.0, ap - bt2, ap)
        _, yh1 = _corner_transform(corner, 0.0, ap, ap)
        y_blade_lo = min(yh0, yh1)
        y_blade_hi = max(yh0, yh1)

        # Arm tip positions (primary axis)
        _, yv_tip = _corner_transform(corner, 0.0, ap - bhs, ap)
        _, yv_ctr = _corner_transform(corner, 0.0, ap, ap)
        y_arm_lo = min(yv_tip, yv_ctr)
        y_arm_hi = max(yv_tip, yv_ctr)
        xh_ctr, _ = _corner_transform(corner, 0.0, 0.0, ap)
        xh_tip, _ = _corner_transform(corner, bhs, 0.0, ap)
        x_arm_lo = min(xh_ctr, xh_tip)
        x_arm_hi = max(xh_ctr, xh_tip)

        eps_s = 1e-8
        cross_on_left = corner in ("north-west", "south-west")
        cross_on_right = corner in ("north-east", "south-east")
        cross_on_top = corner in ("north-west", "north-east")
        cross_on_bottom = corner in ("south-west", "south-east")

        # Left-side stubs (below/above vertical arm tip)
        if cross_on_left:
            blade_w = x_blade_hi - x_blade_lo
            arm_lo_c = max(y_arm_lo, y0)
            arm_hi_c = min(y_arm_hi, y1)
            stub_h = arm_lo_c - y0
            if stub_h > eps_s and blade_w > eps_s:
                wing_footprints.append((
                    (x_blade_lo, y0, x_blade_hi, arm_lo_c),
                    "CROSS_WING_2_STUB",
                ))
            stub_h = y1 - arm_hi_c
            if stub_h > eps_s and blade_w > eps_s:
                wing_footprints.append((
                    (x_blade_lo, arm_hi_c, x_blade_hi, y1),
                    "CROSS_WING_2_STUB",
                ))

        # Right-side stubs
        if cross_on_right:
            blade_w = x_blade_hi - x_blade_lo
            arm_lo_c = max(y_arm_lo, y0)
            arm_hi_c = min(y_arm_hi, y1)
            stub_h = arm_lo_c - y0
            if stub_h > eps_s and blade_w > eps_s:
                wing_footprints.append((
                    (x_blade_lo, y0, x_blade_hi, arm_lo_c),
                    "CROSS_WING_2_STUB",
                ))
            stub_h = y1 - arm_hi_c
            if stub_h > eps_s and blade_w > eps_s:
                wing_footprints.append((
                    (x_blade_lo, arm_hi_c, x_blade_hi, y1),
                    "CROSS_WING_2_STUB",
                ))

        # Top-side stubs (left/right of horizontal arm tip)
        if cross_on_top:
            blade_h = y_blade_hi - y_blade_lo
            arm_lo_c = max(x_arm_lo, x0)
            arm_hi_c = min(x_arm_hi, x1)
            stub_w = arm_lo_c - x0
            if stub_w > eps_s and blade_h > eps_s:
                wing_footprints.append((
                    (x0, y_blade_lo, arm_lo_c, y_blade_hi),
                    "CROSS_WING_1_STUB",
                ))
            stub_w = x1 - arm_hi_c
            if stub_w > eps_s and blade_h > eps_s:
                wing_footprints.append((
                    (arm_hi_c, y_blade_lo, x1, y_blade_hi),
                    "CROSS_WING_1_STUB",
                ))

        # Bottom-side stubs
        if cross_on_bottom:
            blade_h = y_blade_hi - y_blade_lo
            arm_lo_c = max(x_arm_lo, x0)
            arm_hi_c = min(x_arm_hi, x1)
            stub_w = arm_lo_c - x0
            if stub_w > eps_s and blade_h > eps_s:
                wing_footprints.append((
                    (x0, y_blade_lo, arm_lo_c, y_blade_hi),
                    "CROSS_WING_1_STUB",
                ))
            stub_w = x1 - arm_hi_c
            if stub_w > eps_s and blade_h > eps_s:
                wing_footprints.append((
                    (arm_hi_c, y_blade_lo, x1, y_blade_hi),
                    "CROSS_WING_1_STUB",
                ))
    else:
        rectangles_and_splits = [
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
    assembly_box_cell.update_geometry_from_face(
        GeometryType.TECHNOLOGICAL, partitioned_face
    )

    # ------------------------------------------------------------------
    # Build reference boundary faces for material classification.
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

    # If control cross present, rebuild tube/sheath reference faces for
    # material classification after the MACRO partition.
    ctrl_tube_inner_faces = []
    ctrl_tube_outer_faces = []
    ctrl_inner_rect_faces = []
    ctrl_sheath_rect_faces = []
    ctrl_absorber_mat = None
    ctrl_sheath_mat = None

    if has_cross:
        ctrl_shapes_ref = getattr(assembly_box_cell, "_ctrl_cross_shapes", None)
        if ctrl_shapes_ref is None:
            # Rebuild if not attached
            ctrl_shapes_ref = _build_control_cross_shapes(ctrl, ap)
        ctrl_absorber_mat = ctrl.absorber_material
        ctrl_sheath_mat = ctrl.sheath_material
        solid = ctrl.is_solid  # AT10-style (no sheath cavity)

        for t in ctrl_shapes_ref["absorber_tubes"]:
            ctrl_tube_inner_faces.append(t.inner_circles[0].face)
            if not solid:
                ctrl_tube_outer_faces.append(t.inner_circles[1].face)
        ctrl_inner_rect_faces = [r.face for r in ctrl_shapes_ref["inner_sheath_rectangles"]]
        ctrl_sheath_rect_faces = [r.face for r in ctrl_shapes_ref["sheath_rectangles"]]

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

        # Classify material by geometric containment
        mat_assigned = False
        if has_cross:
            # 1–2. Absorber tubes
            for j in range(len(ctrl_tube_inner_faces)):
                if is_point_inside_shape(pt, ctrl_tube_inner_faces[j]):
                    materials_list[i] = ctrl_absorber_mat
                    mat_assigned = True
                    break
                if not solid:
                    if is_point_inside_shape(pt, ctrl_tube_outer_faces[j]):
                        materials_list[i] = ctrl_sheath_mat
                        mat_assigned = True
                        break
            if not mat_assigned:
                # 3. Inner sheath cavity
                for irf in ctrl_inner_rect_faces:
                    if is_point_inside_shape(pt, irf):
                        if solid:
                            materials_list[i] = ctrl_sheath_mat
                            mat_assigned = True
                            break
                        else:
                            materials_list[i] = "MODERATOR"
                            mat_assigned = True
                            break
            if not mat_assigned:
                # 4. Outer sheath / structural rectangles
                for srf in ctrl_sheath_rect_faces:
                    if is_point_inside_shape(pt, srf):
                        materials_list[i] = ctrl_sheath_mat
                        mat_assigned = True
                        break

        if not mat_assigned:
            # 5–7. Standard box regions
            if is_point_inside_shape(pt, coolant_boundary.face):
                materials_list[i] = "COOLANT"
            elif is_point_inside_shape(pt, channel_box_boundary.face):
                materials_list[i] = "CHANNEL_BOX"
            else:
                materials_list[i] = "MODERATOR"

        # Classify position → MACRO name
        macro = _classify_point_to_macro(
            px, py, x0, y0, x1, y1, pin_pitch, n_cols, n_rows,
            wing_footprints=wing_footprints,
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


def _build_cross_aware_discretization_rects(
    ap, x0, y0, x1, y1, lattice_pitch_x, lattice_pitch_y,
    n_cols, n_rows, cross_corner, ctrl,
    unaffected_side_h,
    unaffected_side_v,
    unaffected_corner,
    narrow_gap_splits_h,
    narrow_gap_splits_v,
    cross_corner_splits,
    stub_splits_h,
    stub_splits_v,
):
    """
    Build the discretization rectangles for the peripheral strips around
    the pin lattice when a control cross is present, with user-configurable
    split counts for each region type.

    The geometric decomposition follows the same logic as
    ``_build_cross_aware_splitting_rects`` (used for IC MACRO assignment):
    cross-affected corners are shrunk, affected side strips are split into
    a narrow gap and blade-width stubs, and unaffected sides/corners remain
    unchanged.  The difference is that split counts are configurable rather
    than hard-coded.

    Split tuples must be pre-permuted by the caller:

    - ``_h`` variants are ``(nx, ny)`` for horizontal strips/stubs.
    - ``_v`` variants are ``(nx, ny)`` for vertical strips/stubs.

    Parameters
    ----------
    ap : float
        Assembly pitch.
    x0, y0, x1, y1 : float
        Pin-lattice footprint corners.
    lattice_pitch_x, lattice_pitch_y : float
        Pin-lattice extents.
    n_cols, n_rows : int
        Number of pin columns / rows.
    cross_corner : str
        Cross centre corner (``"north-west"``, etc.).
    ctrl : ControlCrossModel
        Control cross geometry.
    unaffected_side_h : tuple[int, int]
        ``(nx, ny)`` for unaffected horizontal side strips.
    unaffected_side_v : tuple[int, int]
        ``(nx, ny)`` for unaffected vertical side strips.
    unaffected_corner : tuple[int, int]
        ``(nx, ny)`` for unaffected corner rectangles.
    narrow_gap_splits_h : tuple[int, int]
        ``(nx, ny)`` for horizontal narrow gap strips.
    narrow_gap_splits_v : tuple[int, int]
        ``(nx, ny)`` for vertical narrow gap strips.
    cross_corner_splits : tuple[int, int]
        ``(nx, ny)`` for the cross corner gap rectangle.
    stub_splits_h : tuple[int, int]
        ``(nx, ny)`` for horizontal stub regions.
    stub_splits_v : tuple[int, int]
        ``(nx, ny)`` for vertical stub regions.

    Returns
    -------
    list of (Rectangle, (nx, ny))
        Splitting rectangles with their grid split counts.
    """
    bt = ctrl.blade_thickness
    bhs = ctrl.blade_half_span
    bt2 = bt / 2.0

    eps = 1e-8

    # ------------------------------------------------------------------
    # Compute blade and arm extents in assembly coordinates
    # ------------------------------------------------------------------
    xv0, _ = _corner_transform(cross_corner, 0.0, 0.0, ap)
    xv1, _ = _corner_transform(cross_corner, bt2, 0.0, ap)
    x_blade_lo = min(xv0, xv1)
    x_blade_hi = max(xv0, xv1)

    _, yh0 = _corner_transform(cross_corner, 0.0, ap - bt2, ap)
    _, yh1 = _corner_transform(cross_corner, 0.0, ap, ap)
    y_blade_lo = min(yh0, yh1)
    y_blade_hi = max(yh0, yh1)

    _, yv_tip = _corner_transform(cross_corner, 0.0, ap - bhs, ap)
    _, yv_ctr = _corner_transform(cross_corner, 0.0, ap, ap)
    y_arm_lo = min(yv_tip, yv_ctr)
    y_arm_hi = max(yv_tip, yv_ctr)

    xh_ctr, _ = _corner_transform(cross_corner, 0.0, 0.0, ap)
    xh_tip, _ = _corner_transform(cross_corner, bhs, 0.0, ap)
    x_arm_lo = min(xh_ctr, xh_tip)
    x_arm_hi = max(xh_ctr, xh_tip)

    cross_on_left = cross_corner in ("north-west", "south-west")
    cross_on_right = cross_corner in ("north-east", "south-east")
    cross_on_top = cross_corner in ("north-west", "north-east")
    cross_on_bottom = cross_corner in ("south-west", "south-east")

    rects = []

    # ==================================================================
    # CORNERS
    # ==================================================================

    # ---- Bottom-left corner ----
    if cross_corner == "south-west":
        gw = x0 - x_blade_hi
        gh = y0 - y_blade_hi
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x_blade_hi + gw / 2.0,
                                  y_blade_hi + gh / 2.0, 0.0)),
                cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=x0,
                      center=(x0 / 2.0, y0 / 2.0, 0.0)),
            unaffected_corner,
        ))

    # ---- Bottom-right corner ----
    if cross_corner == "south-east":
        gw = x_blade_lo - x1
        gh = y0 - y_blade_hi
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x1 + gw / 2.0,
                                  y_blade_hi + gh / 2.0, 0.0)),
                cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=(ap - x1),
                      center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
            unaffected_corner,
        ))

    # ---- Top-left corner ----
    if cross_corner == "north-west":
        gw = x0 - x_blade_hi
        gh = y_blade_lo - y1
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x_blade_hi + gw / 2.0,
                                  y1 + gh / 2.0, 0.0)),
                cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=x0,
                      center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
            unaffected_corner,
        ))

    # ---- Top-right corner ----
    if cross_corner == "north-east":
        gw = x_blade_lo - x1
        gh = y_blade_lo - y1
        if gw > eps and gh > eps:
            rects.append((
                Rectangle(height=gh, width=gw,
                          center=(x1 + gw / 2.0,
                                  y1 + gh / 2.0, 0.0)),
                cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=(ap - x1),
                      center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
            unaffected_corner,
        ))

    # ==================================================================
    # SIDE STRIPS
    # ==================================================================

    # ---- Bottom-middle strip ----
    if cross_on_bottom:
        gap_h = y0 - y_blade_hi
        if gap_h > eps:
            rects.append((
                Rectangle(height=gap_h, width=lattice_pitch_x,
                          center=((x0 + x1) / 2.0,
                                  y_blade_hi + gap_h / 2.0, 0.0)),
                narrow_gap_splits_h,
            ))
        # Stubs in the blade region
        blade_h = y_blade_hi
        arm_lo_c = max(x_arm_lo, x0)
        arm_hi_c = min(x_arm_hi, x1)
        stub_w = arm_lo_c - x0
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(x0 + stub_w / 2.0,
                                  blade_h / 2.0, 0.0)),
                stub_splits_h,
            ))
        stub_w = x1 - arm_hi_c
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(arm_hi_c + stub_w / 2.0,
                                  blade_h / 2.0, 0.0)),
                stub_splits_h,
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=lattice_pitch_x,
                      center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
            unaffected_side_h,
        ))

    # ---- Top-middle strip ----
    if cross_on_top:
        gap_h = y_blade_lo - y1
        if gap_h > eps:
            rects.append((
                Rectangle(height=gap_h, width=lattice_pitch_x,
                          center=((x0 + x1) / 2.0,
                                  y1 + gap_h / 2.0, 0.0)),
                narrow_gap_splits_h,
            ))
        blade_h = ap - y_blade_lo
        arm_lo_c = max(x_arm_lo, x0)
        arm_hi_c = min(x_arm_hi, x1)
        stub_w = arm_lo_c - x0
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(x0 + stub_w / 2.0,
                                  y_blade_lo + blade_h / 2.0, 0.0)),
                stub_splits_h,
            ))
        stub_w = x1 - arm_hi_c
        if stub_w > eps and blade_h > eps:
            rects.append((
                Rectangle(height=blade_h, width=stub_w,
                          center=(arm_hi_c + stub_w / 2.0,
                                  y_blade_lo + blade_h / 2.0, 0.0)),
                stub_splits_h,
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=lattice_pitch_x,
                      center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
            unaffected_side_h,
        ))

    # ---- Middle-left strip ----
    if cross_on_left:
        gap_w = x0 - x_blade_hi
        if gap_w > eps:
            rects.append((
                Rectangle(height=lattice_pitch_y, width=gap_w,
                          center=(x_blade_hi + gap_w / 2.0,
                                  (y0 + y1) / 2.0, 0.0)),
                narrow_gap_splits_v,
            ))
        blade_w = x_blade_hi - x_blade_lo
        arm_lo_c = max(y_arm_lo, y0)
        arm_hi_c = min(y_arm_hi, y1)
        stub_h = arm_lo_c - y0
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  y0 + stub_h / 2.0, 0.0)),
                stub_splits_v,
            ))
        stub_h = y1 - arm_hi_c
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  arm_hi_c + stub_h / 2.0, 0.0)),
                stub_splits_v,
            ))
    else:
        rects.append((
            Rectangle(height=lattice_pitch_y, width=x0,
                      center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
            unaffected_side_v,
        ))

    # ---- Middle-right strip ----
    if cross_on_right:
        gap_w = x_blade_lo - x1
        if gap_w > eps:
            rects.append((
                Rectangle(height=lattice_pitch_y, width=gap_w,
                          center=(x1 + gap_w / 2.0,
                                  (y0 + y1) / 2.0, 0.0)),
                narrow_gap_splits_v,
            ))
        blade_w = x_blade_hi - x_blade_lo
        arm_lo_c = max(y_arm_lo, y0)
        arm_hi_c = min(y_arm_hi, y1)
        stub_h = arm_lo_c - y0
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  y0 + stub_h / 2.0, 0.0)),
                stub_splits_v,
            ))
        stub_h = y1 - arm_hi_c
        if stub_h > eps and blade_w > eps:
            rects.append((
                Rectangle(height=stub_h, width=blade_w,
                          center=(x_blade_lo + blade_w / 2.0,
                                  arm_hi_c + stub_h / 2.0, 0.0)),
                stub_splits_v,
            ))
    else:
        rects.append((
            Rectangle(height=lattice_pitch_y, width=(ap - x1),
                      center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
            unaffected_side_v,
        ))

    return rects


def discretize_box(assembly_box_cell, assembly_model, box_discretization_config):
    """
    Subdivide the assembly-box peripheral regions into a grid of
    sub-faces for MOC tracking.

    The pin-lattice footprint is computed from the assembly dimensions.
    The surrounding area is split into rectangular strips and each strip
    is further gridded into sub-faces using ``make_grid_faces``.

    When no control cross is present, the standard 8-strip layout is
    used (4 corners + 4 sides).

    When a control cross is present, the affected sides are decomposed
    into narrow gap strips, blade-width stubs, and a cross corner
    rectangle, each with independently configurable split counts
    (resolved from ``box_discretization_config.cross_discretization``).

    Unlike ``subdivide_box_into_macros`` (used for the IC method), this
    function does **not** assign MACRO properties.  Material properties
    are inherited from the already-partitioned box cell via
    ``update_geometry_from_face``.

    Parameters
    ----------
    assembly_box_cell : RectCell
        The assembly box cell as returned by ``build_assembly_box``
        (3-region for uncontrolled, multi-region for controlled).
    assembly_model : CartesianAssemblyModel
        Assembly model providing dimensional information.
    box_discretization_config : BoxDiscretizationConfig
        Configuration specifying the grid splits for corner, side,
        and (optionally) cross-affected strips.

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
    corner, side_h, side_v = box_discretization_config.resolve_splits(
        n_cols, n_rows
    )

    has_cross = getattr(assembly_model, "has_control_cross", False)

    if has_cross:
        # --------------------------------------------------------------
        # Cross-aware discretization
        # --------------------------------------------------------------
        ctrl = assembly_model.control_cross
        bt2 = ctrl.blade_thickness / 2.0

        # Wide gap width (unaffected side: lattice edge to assembly edge)
        wide_gap_width = x0  # symmetric assembly → x0 == y0

        # Narrow gap width (blade edge to lattice edge)
        narrow_gap_width = x0 - bt2

        # Cross corner gap rectangle dimensions
        cross_corner_dims = (narrow_gap_width, narrow_gap_width)

        # Typical stub dimensions:
        #   parallel = extent of lattice edge minus arm span within it
        #   perpendicular = blade half-thickness
        bhs = ctrl.blade_half_span
        # Horizontal stubs sit beside the horizontal arm
        stub_par_h = max(0.0, lattice_pitch_x - (bhs - x0))
        # Vertical stubs sit below/above the vertical arm
        stub_par_v = max(0.0, lattice_pitch_y - (bhs - y0))
        stub_perp = bt2

        # Resolve cross splits from density or explicit config
        cross_disc = box_discretization_config.cross_discretization
        gap_ref = box_discretization_config.gap_splits or (n_cols, 1)

        if cross_disc is not None:
            narrow_gap, cc_splits, stub = cross_disc.resolve(
                gap_splits=gap_ref,
                lattice_pitch=lattice_pitch_x,
                wide_gap_width=wide_gap_width,
                narrow_gap_width=narrow_gap_width,
                cross_corner_dims=cross_corner_dims,
                stub_dims=(stub_par_h, stub_perp),
            )
        else:
            # Auto-compute from gap density
            from ..DDModel.DragonCalculationScheme import CrossDiscretizationConfig
            _auto = CrossDiscretizationConfig()
            narrow_gap, cc_splits, stub = _auto.resolve(
                gap_splits=gap_ref,
                lattice_pitch=lattice_pitch_x,
                wide_gap_width=wide_gap_width,
                narrow_gap_width=narrow_gap_width,
                cross_corner_dims=cross_corner_dims,
                stub_dims=(stub_par_h, stub_perp),
            )

        # Permute for horizontal vs vertical orientation
        narrow_gap_splits_h = narrow_gap                     # (n_par, n_perp)
        narrow_gap_splits_v = (narrow_gap[1], narrow_gap[0]) # permuted
        stub_splits_h = stub
        stub_splits_v = (stub[1], stub[0])

        rectangles_and_splits = _build_cross_aware_discretization_rects(
            ap, x0, y0, x1, y1,
            lattice_pitch_x, lattice_pitch_y,
            n_cols, n_rows,
            ctrl.center, ctrl,
            unaffected_side_h=side_h,
            unaffected_side_v=side_v,
            unaffected_corner=corner,
            narrow_gap_splits_h=narrow_gap_splits_h,
            narrow_gap_splits_v=narrow_gap_splits_v,
            cross_corner_splits=cc_splits,
            stub_splits_h=stub_splits_h,
            stub_splits_v=stub_splits_v,
        )

        print(f"discretize_box [cross-aware]: narrow_gap={narrow_gap}, "
              f"cross_corner={cc_splits}, stub={stub}, "
              f"unaffected_side_h={side_h}, corner={corner}")
    else:
        # --------------------------------------------------------------
        # Standard 8 peripheral rectangles (no control cross)
        # --------------------------------------------------------------
        rectangles_and_splits = [
            # Bottom-left corner
            (Rectangle(height=y0, width=x0,
                       center=(x0 / 2.0, y0 / 2.0, 0.0)),
             corner),
            # Bottom-middle strip
            (Rectangle(height=y0, width=lattice_pitch_x,
                       center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
             side_h),
            # Bottom-right corner
            (Rectangle(height=y0, width=(ap - x1),
                       center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
             corner),
            # Middle-left strip
            (Rectangle(height=lattice_pitch_y, width=x0,
                       center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
             side_v),
            # Middle-right strip
            (Rectangle(height=lattice_pitch_y, width=(ap - x1),
                       center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
             side_v),
            # Top-left corner
            (Rectangle(height=(ap - y1), width=x0,
                       center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
             corner),
            # Top-middle strip
            (Rectangle(height=(ap - y1), width=lattice_pitch_x,
                       center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
             side_h),
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
    if not has_cross:
        print(f"discretize_box: split assembly box into "
              f"{n_subfaces} sub-regions (corner={corner}, "
              f"side_h={side_h}, side_v={side_v}).")
    else:
        print(f"discretize_box: split assembly box into "
              f"{n_subfaces} sub-regions (cross-aware).")

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

