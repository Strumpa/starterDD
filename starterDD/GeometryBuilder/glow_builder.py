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
import numpy as np
import os
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
    lattice_components : dictionary mapping (col_idx, row_idx) tuple to [RectCell, FuelPinModel]
        Allows to keep track of of the pin model describing each cell to ensure proper order of addition to the lattice in 
        ``add_cells_to_regular_lattice`` (non-generating cells are added first, then generating cells, to enforce correct mix numbering in SALOME).
    """
    # Import here to avoid circular import issues
    from ..DDModel.DragonModel import FuelPinModel
    
    lattice_components = {} # change to a dictionaty to store cells and their pin model, at a given position. This will allow to keep track of the pin models associated with each cell.
    pitch = assemblyModel.pin_geometry_dict["pin_pitch"]
    
    fuel_material_mixtures = assemblyModel.fuel_material_mixtures
    row_idx = -1
    for row in assemblyModel.lattice:
        row_idx += 1 
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
                # Recover FuelPinModel technologocal radii to assign materials in the correct order (from innermost to outermost regions)
                techo_radii = pin.technological_radii
                fuel_radius = techo_radii[0]
                gap_radius = techo_radii[1] if len(techo_radii) > 1 else None
                clad_radius = techo_radii[2] if len(techo_radii) > 2 else None
                if gap_radius is not None and gap_radius < fuel_radius:
                    # gap is inner most region, so the order of materials from innermost to outermost is : gap, fuel zones, clad (if clad radius provided and larger than fuel radius), coolant
                    list_of_cell_mats = ["GAP"] + [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures]
                    if clad_radius is not None and clad_radius > fuel_radius:
                        list_of_cell_mats.append("CLAD")
                    list_of_cell_mats.append("COOLANT")
                elif gap_radius is not None and gap_radius > fuel_radius and clad_radius is not None and clad_radius > fuel_radius:
                    # Fuel regions are inner most, then gap, then clad as outer most solid region, coolant is outside: fuel zones, gap, clad, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["GAP", "CLAD", "COOLANT"]
                elif gap_radius is None and clad_radius is not None and clad_radius > fuel_radius:
                    # No gap, clad is outer most solid region, so the order of materials from innermost to outermost is : fuel zones, clad, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["CLAD", "COOLANT"]
                elif gap_radius is not None and gap_radius > fuel_radius and clad_radius is None:
                    # No clad, gap is outer most solid region, so the order of materials from innermost to outermost is : fuel zones, gap, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["GAP", "COOLANT"]
                elif gap_radius is None and clad_radius is None:
                    # No gap, no clad, so only fuel zones and coolant, order of materials from innermost to outermost is : fuel zones, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["COOLANT"]
                else:
                    raise ValueError(
                        f"Invalid combination of radii: fuel_radius={fuel_radius}, gap_radius={gap_radius}, clad_radius={clad_radius}"
                    )

                # Apply sectorization from calculation step if provided
                if calculation_step is not None:
                    sector_cfg = calculation_step.get_sectorization_for_pin(pin, isGd=pin.isGd)
                    if sector_cfg is not None:
                        tmp_cell.sectorize(sector_cfg.sectors, sector_cfg.angles, windmill=sector_cfg.windmill)

                tmp_cell.set_properties({
                    PropertyType.MATERIAL: list_of_cell_mats,
                    PropertyType.MACRO: [f"MACRO{row_idx}{cell_idx}"] * len(list_of_cell_mats)
                })
                # cell_idx in row : column number, ie position along the x direction,
                # row_idx in lattice: row number, ie position along the y direction
                lattice_components[(cell_idx, row_idx)] = [tmp_cell, pin] # store both the cell and its associated pin model for later reference
    return lattice_components


def add_cells_to_regular_lattice(lattice, lattice_components, cell_pitch, translation_x=0.0, translation_y=0.0):
    """
    Add fuel cells to the lattice, skipping water rod placeholders.
    Generating cells are added last in order to enforce order of mix attribution in SALOME (cells added last are assigned first mix numbers)

    Parameters
    ----------
    lattice : Lattice
        The lattice to which cells will be added
    lattice_components : dict
        Dictionary mapping positions to lists of RectCell objects and their associated pin models
    cell_pitch : float
        Pitch of each cell in the lattice
    translation_x : float
        X-axis translation offset to apply to cell positions (supports asymmetric gaps)
    translation_y : float
        Y-axis translation offset to apply to cell positions (supports asymmetric gaps)
    """
    # Import here to avoid circular import issues
    from ..DDModel.DragonModel import FuelPinModel

    for pos, cell_and_pin in lattice_components.items():
        cell, pin = cell_and_pin
        if pin is not None and isinstance(pin, FuelPinModel):
            if pin.isGeneratingCell is False: # add all non generating cells first
                lattice.add_cell(
                    cell, ((pos[0] + 0.5) * cell_pitch + translation_x,
                            (pos[1] + 0.5) * cell_pitch + translation_y,
                            0.0)
                )
    for pos, cell_and_pin in lattice_components.items():
        cell, pin = cell_and_pin
        if pin is not None and isinstance(pin, FuelPinModel):
            if pin.isGeneratingCell: # add generating cells last
                lattice.add_cell(
                    cell, ((pos[0] + 0.5) * cell_pitch + translation_x,
                            (pos[1] + 0.5) * cell_pitch + translation_y,
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


def create_and_add_water_rods_to_lattice(lattice, assembly_model, translation_x=0.0, translation_y=0.0, windmill=False, calculation_step=None):
    """
    Create water rod cells from the assembly model and add them to the lattice at their centers.

    Parameters:
    -----------
    lattice : Lattice
        The lattice to which water rod cells will be added
    assembly_model : CartesianAssemblyModel
        The assembly model containing the water rod geometry parameters
        (water_rod_type, water_rods list with center, radii, materials, etc.)
    translation_x : float
        Unused. Water rod centers are already in assembly coordinates (from YAML).
        Kept for function signature consistency with pin positioning.
    translation_y : float
        Unused. Water rod centers are already in assembly coordinates (from YAML).
        Kept for function signature consistency with pin positioning.
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

            # --- Determine extra moderator radii from calculation step ---
            extra_radii = []
            wr_sectors = None
            if calculation_step is not None:
                wr_sectors = calculation_step.get_water_rod_sectorization()
                if wr_sectors is not None:
                    if wr_sectors.additional_radial_splits_in_moderator:
                        extra_radii = wr_sectors.resolve_water_rod_radii(
                            water_rod_model.inner_radius
                        )

            # Add circles: extra moderator sub-rings, then inner, then outer
            for r in extra_radii:
                tmp_cell.add_circle(r)
            tmp_cell.add_circle(water_rod_model.inner_radius)
            tmp_cell.add_circle(water_rod_model.outer_radius)

            # Build material list: one moderator entry per sub-ring + base 3
            n_extra = len(extra_radii)
            materials = (
                [water_rod_model.moderator_material_name] * (1 + n_extra)
                + [water_rod_model.cladding_material_name,
                   water_rod_model.coolant_material_name]
            )
            n_regions = len(materials)
            tmp_cell.set_properties({
                PropertyType.MATERIAL: materials,
                PropertyType.MACRO: [f"MACRO_{water_rod_model.rod_ID}"] * n_regions,
            })

            # Apply sectorization: prefer calculation_step config, fall back to windmill flag
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
                expanded_s, expanded_a = wr_sectors.expanded_sectors_and_angles(
                    water_rod_model.inner_radius
                )
                tmp_cell.sectorize(expanded_s, expanded_a, windmill=wr_sectors.windmill)
            elif windmill:
                tmp_cell.sectorize([1, 1, 8], [0, 0, 0], windmill=True)
            split_coolant_corners = wr_sectors.subdivisions_coolant_corners if wr_sectors is not None else False
            if split_coolant_corners:
                # circular water rods with sectorization : glow does not allow to sub mesh the coolant
                # at the square corners further than the 16 angles of the .sectorize method.
                # This leads to potentially large coolant regions in the coreners where no inner circle could be added.
                # Need to add extra splitting faces : 
                # compute base point where 16-sector splits intersects with the outer square boundary 
                # add n splitting faces that split base point to corner evenly in parallel splits ?
                # For top right corner :
                alpha = 360.0 / 16.0 # angle of each sector
                adj = water_rod_model.bounding_box_side_length / 2.0
                top_right_corner = (adj, adj, 0.0)
                opp = adj * np.tan(np.radians(alpha))
                base_pt_1 = (opp, adj, 0.0)
                distance_to_split = adj - opp
                # symmetric along y=x
                base_pt_1_sym = (adj, opp, 0.0)
                n_corner_splits = wr_sectors.subdivisions_coolant_corners # this would be retrieved from calculation step config in a more complete implementation
                delta_split = distance_to_split / n_corner_splits
                splitting_faces = []
                for i in range(n_corner_splits):
                    split_pt_1 = (base_pt_1[0] + i * delta_split, adj, 0.0)
                    split_pt_2 = (adj, base_pt_1_sym[1] + i * delta_split, 0.0)
                    splitting_face = make_edge(
                        make_vertex(split_pt_1),
                        make_vertex(split_pt_2),
                    )
                    splitting_faces.append(splitting_face)
                
                    # split the top left corner now : reflect the split points across y axis
                    split_pt_1 = (-base_pt_1[0] - i * delta_split, adj, 0.0)
                    split_pt_2 = (-adj, base_pt_1_sym[1] + i * delta_split, 0.0)
                    splitting_face = make_edge(
                        make_vertex(split_pt_1),
                        make_vertex(split_pt_2),
                    )
                    splitting_faces.append(splitting_face)
                    
                    # split the bottom right corner now : reflect the split points across x axis
                    split_pt_1 = (base_pt_1[0] + i * delta_split, -adj, 0.0)
                    split_pt_2 = (adj, -base_pt_1_sym[1] - i * delta_split, 0.0)
                    splitting_face = make_edge(
                        make_vertex(split_pt_1),
                        make_vertex(split_pt_2),
                    )
                    splitting_faces.append(splitting_face)
                    
                    # split the bottom left corner now : reflect the split points across both axis
                    split_pt_1 = (-base_pt_1[0] - i * delta_split, -adj, 0.0)
                    split_pt_2 = (-adj, -base_pt_1_sym[1] - i * delta_split, 0.0)
                    splitting_face = make_edge(
                        make_vertex(split_pt_1),
                        make_vertex(split_pt_2),
                    )
                    splitting_faces.append(splitting_face)
                    
                    
                    
                re_partitioned = make_partition(
                        [tmp_cell.face],
                        splitting_faces,
                        shape_type=ShapeType.COMPOUND,
                    )
                tmp_cell.update_geometry_from_face(
                        GeometryType.TECHNOLOGICAL, re_partitioned,
                    )
                    

        elif assembly_model.water_rod_type == "square":
            tmp_cell = _build_square_water_rod_cell(
                water_rod_model, calculation_step=calculation_step,
            )

        # water_rod_model.center is in assembly coordinates (from YAML).
        # No translation applied - centers are already positioned within the assembly frame [0, assembly_pitch].
        cx, cy = water_rod_model.center
        lattice.add_cell(
            tmp_cell,
            (cx, cy, 0.0),
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
    # check of output path exists and create if not
    cwd = os.getcwd()
    if not os.path.isabs(output_path):
        output_path = os.path.join(cwd, output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    if export_macro:
        properties_to_export = [PropertyType.MATERIAL, PropertyType.MACRO]
        output_file_name = f"{output_file_name}_{tracking_option}_MACRO"
    else:
        properties_to_export = [PropertyType.MATERIAL]
        output_file_name = f"{output_file_name}_{tracking_option}"

    full_tdt_path = os.path.join(output_path, output_file_name)

    if tracking_option == "TISO":
        lattice.type_geo = LatticeGeometryType.ISOTROPIC
        analyse_and_generate_tdt(
            [lattice], full_tdt_path, TdtSetup(GeometryType.SECTORIZED,
                                             property_types=properties_to_export,
                                             type_geo=LatticeGeometryType.ISOTROPIC,
                                             symmetry_type=SymmetryType.FULL))
    elif tracking_option == "TSPC":
        lattice.type_geo = LatticeGeometryType.RECTANGLE_SYM    
        analyse_and_generate_tdt(
            [lattice], full_tdt_path, TdtSetup(GeometryType.SECTORIZED,
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
    if rw>1e-6:
        split_rects.append(Rectangle(
            name="CTRL_SPLIT_H_RIGHT",
            height=bt / 2.0, width=rw,
            center=(cx, cy, 0.0),
        ))

    # Bottom part of west arm (covers tubes region)
    bh = last_boundary
    cx, cy = ct(bt / 4.0, ap - lw / 2.0)
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


def _compute_asymmetric_coolant_channel_box_rects(assembly_model, center):
    """
    Compute coolant and channel box rectangles with asymmetric gap support.

    With asymmetric gaps (gap_wide ≠ gap_narrow), the rectangles must be
    off-centered to account for different moderator widths on different sides.
    This function derives the correct dimensions and center positions based on
    the detected lattice symmetry.

    For **anti-diagonal symmetry** (top-left to bottom-right):
        - X-axis: gap_wide on left, gap_narrow on right → asymmetric
        - Y-axis: gap_narrow on bottom, gap_wide on top → asymmetric (reversed)
        - Rectangle dimensions differ on X and Y

    For **main-diagonal symmetry** (transpose):
        - Both axes use gap_wide (symmetric despite different gaps)
        - Rectangles remain centered

    For **no symmetry** or **quarter/eighth symmetry**:
        - Treat as isotropic (both axes use gap_wide)
        - Rectangles remain centered (backward compatible)

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model providing gap_wide, gap_narrow, channel_box_thickness,
        and corner_inner_radius_of_curvature
    center : tuple
        (x, y, z) center of the assembly box

    Returns
    -------
    coolant_rect : Rectangle
        Inner coolant boundary with proper centering for asymmetric gaps
    channel_box_rect : Rectangle
        Outer channel box boundary with proper centering for asymmetric gaps
    """
    from glow.geometry_layouts.geometries import Rectangle

    ap = assembly_model.assembly_pitch
    cbt = assembly_model.channel_box_thickness
    gap_wide = assembly_model.gap_wide
    gap_narrow = assembly_model.gap_narrow
    corner_r_inner = assembly_model.corner_inner_radius_of_curvature

    # Detect lattice symmetry to determine offset strategy
    sym_type = assembly_model.check_diagonal_symmetry()

    if sym_type == "anti-diagonal":
        # Asymmetric configuration: wide-wide corner on top-left, narrow-narrow on bottom-right
        # X-axis: gap_wide on left (low x), gap_narrow on right (high x)
        # Y-axis: gap_narrow on bottom (low y), gap_wide on top (high y)

        channel_box_outer_x = ap - gap_wide - gap_narrow
        channel_box_outer_y = ap - gap_narrow - gap_wide 

        channel_box_inner_x = channel_box_outer_x - 2.0 * cbt
        channel_box_inner_y = channel_box_outer_y - 2.0 * cbt

        # Center offset: shift to balance the asymmetric gaps
        # For anti-diagonal: gap_wide on left/top, gap_narrow on right/bottom
        # Channel box occupies:
        #   X range: [gap_wide + cbt, ap - gap_narrow - cbt]
        #   Y range: [gap_narrow + cbt, ap - gap_wide - cbt]
        # Center X = (gap_wide + ap - gap_narrow) / 2
        # Center Y = (gap_narrow + ap - gap_wide) / 2
        # Offsets from (ap/2, ap/2):
        offset_x = (gap_wide - gap_narrow) / 2.0
        offset_y = (gap_narrow - gap_wide) / 2.0

        rect_center = (center[0] + offset_x, center[1] + offset_y, center[2])

    elif sym_type == "main-diagonal":
        # Symmetric on both axes (both use gap_wide due to symmetry)
        channel_box_outer_x = ap - 2.0 * gap_wide
        channel_box_outer_y = ap - 2.0 * gap_wide

        channel_box_inner_x = channel_box_outer_x - 2.0 * cbt
        channel_box_inner_y = channel_box_outer_y - 2.0 * cbt

        rect_center = center

    else:
        # No symmetry or quarter/eighth symmetry
        # Treat as isotropic (both axes use gap_wide)
        channel_box_outer_x = ap - 2.0 * gap_wide
        channel_box_outer_y = ap - 2.0 * gap_wide

        channel_box_inner_x = channel_box_outer_x - 2.0 * cbt
        channel_box_inner_y = channel_box_outer_y - 2.0 * cbt

        rect_center = center

    # Configure rounded corners if applicable
    if corner_r_inner > 0.0:
        corner_r_outer = corner_r_inner + cbt
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

    # Build rectangles with potentially asymmetric dimensions
    coolant_rect = Rectangle(
        name="intra_assembly_coolant",
        height=channel_box_inner_y,
        width=channel_box_inner_x,
        center=rect_center,
        rounded_corners=rounded_corners_coolant,
    )

    channel_box_rect = Rectangle(
        name="channel_box",
        height=channel_box_outer_y,
        width=channel_box_outer_x,
        center=rect_center,
        rounded_corners=rounded_corners_chanbox,
    )

    return coolant_rect, channel_box_rect


def build_assembly_box(assembly_model, center=None):
    """
    Build the assembly box cell from the assembly model dimensions.

    Supports both symmetric and asymmetric gap configurations. With asymmetric
    gaps (gap_wide ≠ gap_narrow), the coolant and channel box rectangles are
    automatically off-centered to account for different moderator widths on
    different sides, based on the detected lattice symmetry.

    Without a control cross the result is a 3-region cell (intra-assembly
    coolant / channel box / inter-assembly moderator).

    When ``assembly_model.has_control_cross`` is ``True`` the control
    cross shapes (sheath, inner cavity, absorber tubes) are also
    included in the partition and materials are assigned by geometric
    containment against all reference shapes.

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model providing ``assembly_pitch``, ``gap_wide``, ``gap_narrow``,
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
    cbt = assembly_model.channel_box_thickness

    if center is None:
        center = (ap / 2.0, ap / 2.0, 0.0)

    # Compute asymmetry-aware coolant and channel box rectangles
    # This handles both symmetric (backward compatible) and asymmetric gap configurations
    coolant_rect, channel_box_rect = _compute_asymmetric_coolant_channel_box_rects(
        assembly_model, center
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


def _reassign_materials_by_containment(assembly_box_cell, assembly_model):
    """
    Reassign MATERIAL properties to every sub-face of the assembly box
    cell using geometric containment against freshly built reference
    shapes.

    This function is intended to be called **after** a discretization
    partition (e.g. in ``discretize_box``) to guarantee correct
    material assignments regardless of how glow's internal
    ``update_geometry_from_face`` propagates properties across
    successive partitions.

    The classification logic mirrors ``_assign_materials_with_control_cross``
    and the material loop in ``subdivide_box_into_macros``:

    Priority order (innermost first):

    1. Absorber tube inner circle → absorber material
    2. Absorber tube outer circle → sheath material (tube cladding)
    3. Inside inner-sheath rectangle but outside tubes → MODERATOR
       (inter-tube gap; sheath material for solid crosses)
    4. Inside outer-sheath / structure rectangle → sheath material
    5. Inside coolant boundary → COOLANT
    6. Inside channel-box boundary → CHANNEL_BOX
    7. Otherwise → MODERATOR

    Parameters
    ----------
    assembly_box_cell : RectCell
        The assembly box cell whose sub-faces need material
        re-classification.
    assembly_model : CartesianAssemblyModel
        Assembly model providing dimensional information and optional
        control cross data.
    """
    from collections import Counter

    ap = assembly_model.assembly_pitch
    center = (ap / 2.0, ap / 2.0, 0.0)

    # Use the asymmetry-aware helper to build reference rectangles
    # This ensures material classification matches the actual geometry
    # (which may be off-centered with asymmetric gaps)
    coolant_boundary, channel_box_boundary = _compute_asymmetric_coolant_channel_box_rects(
        assembly_model, center
    )

    # ------------------------------------------------------------------
    # Control cross reference shapes (if present)
    # ------------------------------------------------------------------
    has_cross = getattr(assembly_model, "has_control_cross", False)
    ctrl_tube_inner_faces = []
    ctrl_tube_outer_faces = []
    ctrl_inner_rect_faces = []
    ctrl_sheath_rect_faces = []
    ctrl_absorber_mat = None
    ctrl_sheath_mat = None
    solid = False

    if has_cross:
        ctrl = assembly_model.control_cross
        ctrl_shapes_ref = getattr(
            assembly_box_cell, "_ctrl_cross_shapes", None
        )
        if ctrl_shapes_ref is None:
            ctrl_shapes_ref = _build_control_cross_shapes(ctrl, ap)
        ctrl_absorber_mat = ctrl.absorber_material
        ctrl_sheath_mat = ctrl.sheath_material
        solid = ctrl.is_solid

        # Pre-extract absorber tube circle faces
        for t in ctrl_shapes_ref["absorber_tubes"]:
            ctrl_tube_inner_faces.append(t.inner_circles[0].face)
            if not solid:
                ctrl_tube_outer_faces.append(t.inner_circles[1].face)
        ctrl_inner_rect_faces = [
            r.face for r in ctrl_shapes_ref["inner_sheath_rectangles"]
        ]
        ctrl_sheath_rect_faces = [
            r.face for r in ctrl_shapes_ref["sheath_rectangles"]
        ]

    # ------------------------------------------------------------------
    # Classify every sub-face
    # ------------------------------------------------------------------
    subfaces = assembly_box_cell.extract_subfaces()
    n = len(subfaces)
    materials = [""] * n

    for i, subface in enumerate(subfaces):
        pt = make_vertex_inside_face(subface)

        mat_assigned = False
        if has_cross:
            # 1–2. Absorber tubes
            if solid:
                for cf in ctrl_tube_inner_faces:
                    if is_point_inside_shape(pt, cf):
                        materials[i] = ctrl_absorber_mat
                        mat_assigned = True
                        break
            else:
                for j in range(len(ctrl_tube_inner_faces)):
                    if is_point_inside_shape(pt, ctrl_tube_inner_faces[j]):
                        materials[i] = ctrl_absorber_mat
                        mat_assigned = True
                        break
                    if is_point_inside_shape(pt, ctrl_tube_outer_faces[j]):
                        materials[i] = ctrl_sheath_mat
                        mat_assigned = True
                        break

            if not mat_assigned:
                # 3. Inner sheath cavity (inter-tube moderator) or
                #    solid cross sheath material
                for irf in ctrl_inner_rect_faces:
                    if is_point_inside_shape(pt, irf):
                        if solid:
                            materials[i] = ctrl_sheath_mat
                        else:
                            materials[i] = "MODERATOR"
                        mat_assigned = True
                        break

            if not mat_assigned:
                # 4. Outer sheath / structural rectangles
                for srf in ctrl_sheath_rect_faces:
                    if is_point_inside_shape(pt, srf):
                        materials[i] = ctrl_sheath_mat
                        mat_assigned = True
                        break

        if not mat_assigned:
            # 5–7. Standard box regions
            if is_point_inside_shape(pt, coolant_boundary.face):
                materials[i] = "COOLANT"
            elif is_point_inside_shape(pt, channel_box_boundary.face):
                materials[i] = "CHANNEL_BOX"
            else:
                materials[i] = "MODERATOR"

    assembly_box_cell.set_properties({
        PropertyType.MATERIAL: materials,
    })

    counts = Counter(materials)
    print(f"[reassign_materials] Re-assigned materials to {n} sub-faces: "
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
    # Use translation offsets instead of centered assumption to support asymmetric gaps
    x0 = assembly_model.translation_offset_x if assembly_model.translation_offset_x is not None else (ap - lattice_pitch_x) / 2.0
    y0 = assembly_model.translation_offset_y if assembly_model.translation_offset_y is not None else (ap - lattice_pitch_y) / 2.0
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

        # Compute stub footprints — moderator rectangles in the blade region but outside the arm extent. 
        # These must get their own MACRO so they are not merged with the gap column/row MACROs.
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
    center = (ap / 2.0, ap / 2.0, 0.0)

    # Use the asymmetry-aware helper to build reference rectangles
    # This ensures material classification matches the actual geometry
    # (which may be off-centered with asymmetric gaps)
    coolant_boundary, channel_box_boundary = _compute_asymmetric_coolant_channel_box_rects(
        assembly_model, center
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
    corner_bl,
    corner_br,
    corner_tl,
    corner_tr,
    narrow_gap_splits_h,
    narrow_gap_splits_v,
    moderator_at_cross_corner_splits,
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
    corner_bl : tuple[int, int]
        ``(nx, ny)`` for bottom-left corner rectangle.
    corner_br : tuple[int, int]
        ``(nx, ny)`` for bottom-right corner rectangle.
    corner_tl : tuple[int, int]
        ``(nx, ny)`` for top-left corner rectangle.
    corner_tr : tuple[int, int]
        ``(nx, ny)`` for top-right corner rectangle.
    narrow_gap_splits_h : tuple[int, int]
        ``(nx, ny)`` for horizontal narrow gap strips.
    narrow_gap_splits_v : tuple[int, int]
        ``(nx, ny)`` for vertical narrow gap strips.
    moderator_at_cross_corner_splits : tuple[int, int]
        ``(nx, ny)`` for the moderator corner gap rectangle near the
        control cross corner.
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
                moderator_at_cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=x0,
                      center=(x0 / 2.0, y0 / 2.0, 0.0)),
            corner_bl,
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
                moderator_at_cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=y0, width=(ap - x1),
                      center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
            corner_br,
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
                moderator_at_cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=x0,
                      center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
            corner_tl,
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
                moderator_at_cross_corner_splits,
            ))
    else:
        rects.append((
            Rectangle(height=(ap - y1), width=(ap - x1),
                      center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
            corner_tr,
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


def _build_wing_submesh_rects(ctrl, ap, control_cross_submesh_config, ctrl_shapes):
    """
    Build splitting faces to sub-mesh the control cross wings into
    axial zones.

    Each wing arm is decomposed into three zones by cuts perpendicular
    to the arm axis:

    1. **Corner zone** — the ``bt/2 × bt/2`` square at the cross centre
       where both arms overlap.
    2. **Central-structure-to-absorber zone** — from the corner zone
       edge to the first absorber tube boundary.
    3. **Absorber-pin zone** — optionally split at each tube boundary
       (``extend_splits_at_tube_boundaries``) and optionally bisected
       at each tube centre (``split_tubes_in_half``).

    Zones 1 and 2 can be further gridded by user-specified split counts.
    All splitting faces are ``Rectangle`` objects suitable for passing to
    ``make_grid_faces`` or directly into ``make_partition``.

    Parameters
    ----------
    ctrl : ControlCrossModel
        The control cross model with all geometric dimensions.
    ap : float
        Assembly pitch.
    control_cross_submesh_config : ControlCrossSubmeshConfig
        Configuration specifying the grid splits and options.
    ctrl_shapes : dict
        The shapes dict as returned by ``_build_control_cross_shapes``.

    Returns
    -------
    list of (Rectangle, (nx, ny))
        Splitting rectangles with their grid split counts.
    """
    corner = ctrl.center
    bt = ctrl.blade_thickness
    bt2 = bt / 2.0
    cshs = ctrl.central_structure_half_span
    bhs = ctrl.blade_half_span
    n_tubes = ctrl.number_tubes_per_wing
    delta = ctrl.tube_spacing
    first_offset = ctrl.first_tube_offset
    st = ctrl.sheath_thickness
    inner_w = ctrl.inner_sheath_width  # = bt - 2*st

    # Resolve split counts with defaults
    ctrl_cross_corner_splits = control_cross_submesh_config.control_cross_corner_splits or (1, 1)
    cs_splits = control_cross_submesh_config.central_structure_splits or (1, 1)
    extend_tube = control_cross_submesh_config.extend_splits_at_tube_boundaries
    bisect_tube = control_cross_submesh_config.split_tubes_in_half

    def ct(x, y):
        """Shorthand for corner transform."""
        return _corner_transform(corner, x, y, ap)

    rects = []

    # ------------------------------------------------------------------
    # Helper: determine which arm axis is which.
    # The horizontal arm extends along x (canonical NW: positive x).
    # The vertical arm extends along y (canonical NW: negative y from ap).
    #
    # For the horizontal arm, "along arm" = x, "across arm" = y.
    #   corner_splits = (n_along, n_across) → (nx, ny) for H arm
    # For the vertical arm, "along arm" = y, "across arm" = x.
    #   corner_splits = (n_along, n_across) → (ny, nx) for V arm → permute
    # ------------------------------------------------------------------

    # ==================================================================
    # ZONE A — Corner zone: bt/2 × bt/2 at cross centre
    # Both arms share this zone, so we only create it once.
    # ==================================================================
    cx, cy = ct(bt2 / 2.0, ap - bt2 / 2.0)
    rects.append((
        Rectangle(
            name="WING_SUBMESH_CORNER",
            height=bt2, width=bt2,
            center=(cx, cy, 0.0),
        ),
        ctrl_cross_corner_splits,
    ))

    # ==================================================================
    # ZONE B — Central-structure-to-absorber zone
    # Span from bt/2 to first tube boundary along each arm.
    # ==================================================================
    # First tube boundary (lower edge of first tube bounding box)
    first_tube_boundary = first_offset - delta / 2.0

    # Horizontal arm zone B: x from bt/2 to first_tube_boundary
    zone_b_len_h = first_tube_boundary - bt2
    if zone_b_len_h > 1e-8:
        cx, cy = ct(bt2 + zone_b_len_h / 2.0, ap)
        rects.append((
            Rectangle(
                name="WING_SUBMESH_CS_H",
                height=bt, width=zone_b_len_h,
                center=(cx, cy, 0.0),
            ),
            cs_splits,
        ))

    # Vertical arm zone B: y from (ap - bt/2) to (ap - first_tube_boundary)
    zone_b_len_v = first_tube_boundary - bt2
    if zone_b_len_v > 1e-8:
        cx, cy = ct(0.0, ap - bt2 - zone_b_len_v / 2.0)
        rects.append((
            Rectangle(
                name="WING_SUBMESH_CS_V",
                height=zone_b_len_v, width=bt,
                center=(cx, cy, 0.0),
            ),
            (cs_splits[1], cs_splits[0]),  # permuted for vertical arm
        ))

    # ==================================================================
    # ZONE C — Absorber pin zone
    # Optionally extend tube bounding surfaces to sheath border and/or
    # bisect tubes.
    # ==================================================================
    tubes = ctrl_shapes["absorber_tubes"]

    if extend_tube or bisect_tube:
        for i in range(n_tubes):
            # Get tube centres from the Salome geometry objects.
            # Horizontal tube is at index 2*i, vertical at 2*i+1.
            tube_h = tubes[2 * i]
            tube_v = tubes[2 * i + 1]

            # Retrieve tube centres via Salome GetParameters
            tx_h = float(tube_h.inner_circles[0].o.GetParameters().split(":")[0])
            ty_h = float(tube_h.inner_circles[0].o.GetParameters().split(":")[1])
            tx_v = float(tube_v.inner_circles[0].o.GetParameters().split(":")[0])
            ty_v = float(tube_v.inner_circles[0].o.GetParameters().split(":")[1])

            if extend_tube:
                # Horizontal arm: full bt-wide rectangle at tube centre,
                # spanning the blade thickness.
                # The tube bounding box is inner_w × delta; we extend to
                # bt × delta by creating a full-width splitting face.
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_EXT_H_{i}",
                        height=bt, width=delta,
                        center=(tx_h, ty_h, 0.0),
                    ),
                    (1, 1),
                ))

                # Vertical arm: full bt-wide rectangle at tube centre
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_EXT_V_{i}",
                        height=delta, width=bt,
                        center=(tx_v, ty_v, 0.0),
                    ),
                    (1, 1),
                ))

            if bisect_tube:
                # Horizontal arm: bisect the tube at its centre along x
                # (perpendicular to arm axis = a thin horizontal cut)
                # A thin rectangle spanning bt across, half-delta wide,
                # placed at the tube centre.  We split by creating two
                # half-delta faces.
                half_delta = delta / 2.0
                # Left half
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_BISECT_H_{i}_L",
                        height=bt, width=half_delta,
                        center=(tx_h - half_delta / 2.0, ty_h, 0.0),
                    ),
                    (1, 1),
                ))
                # Right half
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_BISECT_H_{i}_R",
                        height=bt, width=half_delta,
                        center=(tx_h + half_delta / 2.0, ty_h, 0.0),
                    ),
                    (1, 1),
                ))

                # Vertical arm: bisect at tube centre along y
                # Top half
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_BISECT_V_{i}_T",
                        height=half_delta, width=bt,
                        center=(tx_v, ty_v + half_delta / 2.0, 0.0),
                    ),
                    (1, 1),
                ))
                # Bottom half
                rects.append((
                    Rectangle(
                        name=f"WING_TUBE_BISECT_V_{i}_B",
                        height=half_delta, width=bt,
                        center=(tx_v, ty_v - half_delta / 2.0, 0.0),
                    ),
                    (1, 1),
                ))

    n_faces = sum(nx * ny for _, (nx, ny) in rects)
    print(f"[control cross submesh] Built {len(rects)} splitting rects "
          f"({n_faces} total faces) for corner '{corner}', "
          f"ctrl_cross_corner_splits={ctrl_cross_corner_splits}, cs_splits={cs_splits}, "
          f"extend_tube={extend_tube}, bisect_tube={bisect_tube}.")

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
    (resolved from ``box_discretization_config.cross_moderator_discretization``).

    If ``control_cross_submesh`` is enabled in the box discretization
    config, the control cross wings are further subdivided into three
    axial zones per arm (corner, central-structure-to-absorber,
    absorber pin zone) with optional tube boundary extension and tube
    bisection.  See ``ControlCrossSubmeshConfig`` for details.

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
    # Use translation offsets instead of centered assumption to support asymmetric gaps
    x0 = assembly_model.translation_offset_x if assembly_model.translation_offset_x is not None else (ap - lattice_pitch_x) / 2.0
    y0 = assembly_model.translation_offset_y if assembly_model.translation_offset_y is not None else (ap - lattice_pitch_y) / 2.0
    x1 = x0 + lattice_pitch_x
    y1 = y0 + lattice_pitch_y

    # Resolve split counts
    # Check if asymmetric gap splits are configured
    has_asym_gap_splits = (
        box_discretization_config.gap_wide_splits or
        box_discretization_config.gap_narrow_splits or
        box_discretization_config.wide_wide_corner_splits or
        box_discretization_config.narrow_narrow_corner_splits or
        box_discretization_config.mixed_corner_splits
    )

    if has_asym_gap_splits:
        # Use symmetry-aware resolution for asymmetric gaps
        region_splits = box_discretization_config.resolve_splits_with_symmetry(
            n_cols, n_rows, assembly_model
        )
        # Extract splits for later use
        corner_bl = region_splits['corner_bl']
        corner_br = region_splits['corner_br']
        corner_tl = region_splits['corner_tl']
        corner_tr = region_splits['corner_tr']
        side_bottom = region_splits['side_bottom']
        side_top = region_splits['side_top']
        side_left = region_splits['side_left']
        side_right = region_splits['side_right']
    else:
        # Use traditional uniform resolution (backward compatible)
        corner, side_h, side_v = box_discretization_config.resolve_splits(
            n_cols, n_rows
        )
        # Apply uniform splits to all regions
        corner_bl = corner_br = corner_tl = corner_tr = corner
        side_bottom = side_top = side_h
        side_left = side_right = side_v

    has_cross = getattr(assembly_model, "has_control_cross", False)

    if has_cross:
        # Determine which sides are affected by cross placement
        cross_corner = assembly_model.control_cross.center
        cross_on_left = cross_corner in ("north-west", "south-west")
        cross_on_right = cross_corner in ("north-east", "south-east")
        cross_on_top = cross_corner in ("north-west", "north-east")
        cross_on_bottom = cross_corner in ("south-west", "south-east")

        # Get cross-specific discretization if available
        cross_mod_disc = box_discretization_config.cross_moderator_discretization
        gap_ref = box_discretization_config.gap_splits or (n_cols, 1)

        # Resolve cross-specific splits for affected sides
        cross_side_gap = None
        if cross_mod_disc is not None:
            cross_side_gap, _, _ = cross_mod_disc.resolve(
                gap_splits=gap_ref,
                lattice_pitch=lattice_pitch_x,
                wide_gap_width=x0,
                narrow_gap_width=x0 - (assembly_model.control_cross.blade_thickness / 2.0),
                cross_corner_dims=(x0 - assembly_model.control_cross.blade_thickness / 2.0,
                                   y0 - assembly_model.control_cross.blade_thickness / 2.0),
                stub_dims=(0, assembly_model.control_cross.blade_thickness / 2.0),
            )

        # Determine splits for each side (affected vs unaffected)
        # For affected sides, use cross-specific splits; for unaffected, use original splits
        if has_asym_gap_splits:
            # Use region-specific splits from resolve_splits_with_symmetry
            if cross_on_top:
                side_top_used = cross_side_gap if cross_side_gap else region_splits['side_top']
            else:
                side_top_used = region_splits['side_top']

            if cross_on_bottom:
                side_bottom_used = cross_side_gap if cross_side_gap else region_splits['side_bottom']
            else:
                side_bottom_used = region_splits['side_bottom']

            if cross_on_left:
                side_left_used = (cross_side_gap[1], cross_side_gap[0]) if cross_side_gap else region_splits['side_left']
            else:
                side_left_used = region_splits['side_left']

            if cross_on_right:
                side_right_used = (cross_side_gap[1], cross_side_gap[0]) if cross_side_gap else region_splits['side_right']
            else:
                side_right_used = region_splits['side_right']
        else:
            # Use uniform splits but distinguish affected/unaffected
            corner, side_h, side_v = box_discretization_config.resolve_splits(n_cols, n_rows)
            if cross_on_top:
                side_top_used = cross_side_gap if cross_side_gap else side_h
            else:
                side_top_used = side_h

            if cross_on_bottom:
                side_bottom_used = cross_side_gap if cross_side_gap else side_h
            else:
                side_bottom_used = side_h

            if cross_on_left:
                side_left_used = (cross_side_gap[1], cross_side_gap[0]) if cross_side_gap else side_v
            else:
                side_left_used = side_v

            if cross_on_right:
                side_right_used = (cross_side_gap[1], cross_side_gap[0]) if cross_side_gap else side_v
            else:
                side_right_used = side_v

        # For cross-aware building, pass unaffected splits (which could be affected by symmetry)
        # Use unaffected horizontal (TOP/BOTTOM unaffected) and vertical (LEFT/RIGHT unaffected)
        # The _build_cross_aware_discretization_rects function will use narrow_gap_splits for affected
        # and unaffected_side_* for unaffected

        # Choose reference sides for unaffected
        if cross_on_bottom:
            unaffected_side_h = side_top_used  # TOP is unaffected
        else:
            unaffected_side_h = side_bottom_used  # BOTTOM is unaffected

        if cross_on_right:
            unaffected_side_v = side_left_used  # LEFT is unaffected
        else:
            unaffected_side_v = side_right_used  # RIGHT is unaffected

        # Show debug info about which mode is used
        mode_str = "asymmetric gaps" if has_asym_gap_splits else "symmetric"
        print(f"discretize_box [cross-aware, {mode_str}, {cross_corner}]: "
              f"affected: top={cross_on_top}, bottom={cross_on_bottom}, left={cross_on_left}, right={cross_on_right}, "
              f"cross_side_gap={cross_side_gap}, unaffected_h={unaffected_side_h}, unaffected_v={unaffected_side_v}")
    else:
        # Standard 8 peripheral rectangles (no control cross)
        rectangles_and_splits = [
            # Bottom-left corner
            (Rectangle(height=y0, width=x0,
                       center=(x0 / 2.0, y0 / 2.0, 0.0)),
             corner_bl),
            # Bottom-middle strip
            (Rectangle(height=y0, width=lattice_pitch_x,
                       center=((x0 + x1) / 2.0, y0 / 2.0, 0.0)),
             side_bottom),
            # Bottom-right corner
            (Rectangle(height=y0, width=(ap - x1),
                       center=((x1 + ap) / 2.0, y0 / 2.0, 0.0)),
             corner_br),
            # Middle-left strip
            (Rectangle(height=lattice_pitch_y, width=x0,
                       center=(x0 / 2.0, (y0 + y1) / 2.0, 0.0)),
             side_left),
            # Middle-right strip
            (Rectangle(height=lattice_pitch_y, width=(ap - x1),
                       center=((x1 + ap) / 2.0, (y0 + y1) / 2.0, 0.0)),
             side_right),
            # Top-left corner
            (Rectangle(height=(ap - y1), width=x0,
                       center=(x0 / 2.0, (y1 + ap) / 2.0, 0.0)),
             corner_tl),
            # Top-middle strip
            (Rectangle(height=(ap - y1), width=lattice_pitch_x,
                       center=((x0 + x1) / 2.0, (y1 + ap) / 2.0, 0.0)),
             side_top),
            # Top-right corner
            (Rectangle(height=(ap - y1), width=(ap - x1),
                       center=((x1 + ap) / 2.0, (y1 + ap) / 2.0, 0.0)),
             corner_tr),
        ]

    # For control cross, build the cross-aware discretization
    if has_cross:
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

        # Get cross-specific discretization if available (already resolved above for affected sides)
        cross_mod_disc = box_discretization_config.cross_moderator_discretization
        gap_ref = box_discretization_config.gap_splits or (n_cols, 1)

        if cross_mod_disc is not None:
            narrow_gap, cc_splits, stub = cross_mod_disc.resolve(
                gap_splits=gap_ref,
                lattice_pitch=lattice_pitch_x,
                wide_gap_width=wide_gap_width,
                narrow_gap_width=narrow_gap_width,
                cross_corner_dims=cross_corner_dims,
                stub_dims=(stub_par_h, stub_perp),
            )
        else:
            # Auto-compute from gap density
            from ..DDModel.DragonCalculationScheme import CrossModeratorDiscretizationConfig
            _auto = CrossModeratorDiscretizationConfig()
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
            unaffected_side_h=unaffected_side_h,
            unaffected_side_v=unaffected_side_v,
            corner_bl=corner_bl,
            corner_br=corner_br,
            corner_tl=corner_tl,
            corner_tr=corner_tr,
            narrow_gap_splits_h=narrow_gap_splits_h,
            narrow_gap_splits_v=narrow_gap_splits_v,
            moderator_at_cross_corner_splits=cc_splits,
            stub_splits_h=stub_splits_h,
            stub_splits_v=stub_splits_v,
        )

        mode_str = "asymmetric gaps" if has_asym_gap_splits else "symmetric"
        print(f"discretize_box [cross-aware, {mode_str}]: narrow_gap={narrow_gap}, "
              f"moderator_at_cross_corner={cc_splits}, stub={stub}, "
              f"unaffected_side_h={unaffected_side_h}, "
              f"corner_bl={corner_bl}, corner_br={corner_br}, corner_tl={corner_tl}, corner_tr={corner_tr}")

    splitting_faces = []
    for rect, (nx, ny) in rectangles_and_splits:
        splitting_faces.extend(make_grid_faces(rect, nx, ny))

    # ------------------------------------------------------------------
    # Control cross sub-mesh (if enabled)
    # ------------------------------------------------------------------
    if has_cross:
        ctrl_cross_cfg = box_discretization_config.control_cross_submesh
        if ctrl_cross_cfg is not None and ctrl_cross_cfg.enabled:
            ctrl_shapes = getattr(assembly_box_cell, "_ctrl_cross_shapes", None)
            if ctrl_shapes is None:
                import warnings
                warnings.warn(
                    "Control cross sub-mesh requested but "
                    "_ctrl_cross_shapes not attached to "
                    "assembly_box_cell.  Skipping control cross "
                    "sub-meshing.",
                    stacklevel=2,
                )
            else:
                wing_rects = _build_wing_submesh_rects(
                    ctrl, ap, ctrl_cross_cfg, ctrl_shapes,
                )
                for rect, (nx, ny) in wing_rects:
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
        if has_asym_gap_splits:
            print(f"discretize_box: split assembly box into "
                  f"{n_subfaces} sub-regions (asymmetric gaps - "
                  f"corner_bl={corner_bl}, corner_br={corner_br}, "
                  f"corner_tl={corner_tl}, corner_tr={corner_tr}, "
                  f"side_bottom={side_bottom}, side_top={side_top}, "
                  f"side_left={side_left}, side_right={side_right}).")
        else:
            print(f"discretize_box: split assembly box into "
                  f"{n_subfaces} sub-regions (corner={corner_bl}, "
                  f"side_h={side_bottom}, side_v={side_left}).")
    else:
        print(f"discretize_box: split assembly box into "
              f"{n_subfaces} sub-regions (cross-aware).")

    # ------------------------------------------------------------------
    # Explicit material re-assignment by geometric containment
    # ------------------------------------------------------------------
    if box_discretization_config.reassign_materials:
        _reassign_materials_by_containment(assembly_box_cell, assembly_model)

    return assembly_box_cell


def build_full_assembly_geometry(assembly_model, calculation_step,
                                 output_path, output_file_name):
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
        Must have ``translation_offset_x`` and ``translation_offset_y`` computed
        from asymmetric gap configuration (done automatically in model init).
    calculation_step : CalculationStep
        The calculation step whose discretization config drives the geometry.
    output_path : str
        Directory to write the TDT file.
    output_file_name : str
        Base name for the output TDT file (tracking suffix is appended
        automatically by ``export_glow_geom``).

    Returns
    -------
    lattice : Lattice
        The assembled glow ``Lattice`` object, already exported.
    assembly_box_cell : RectCell
        The assembly box cell (with MACRO properties if applicable).
    """
    from ..DDModel.DragonModel import CartesianAssemblyModel  # type check only

    # ----- Extract asymmetric translation offsets from model -----
    translation_x = assembly_model.translation_offset_x if assembly_model.translation_offset_x is not None else 0.0
    translation_y = assembly_model.translation_offset_y if assembly_model.translation_offset_y is not None else 0.0

    # ----- Derived dimensions -----
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_cols = len(assembly_model.lattice_description[0])

    center = (ap / 2.0, ap / 2.0, 0.0)

    # Step 1: Apply radii from the calculation step
    calculation_step.apply_radii(assembly_model)

    # Step 2: Generate fuel cells with sectorization
    ordered_cells = generate_fuel_cells(
        assembly_model, calculation_step=calculation_step
    )
    if ap > pin_pitch * n_cols:
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
    else:        
        print(f"Assembly pitch {ap} is not larger than pin lattice "
            f"footprint {pin_pitch * n_cols}; skipping assembly box.")
        assembly_box_cell = None

    # Step 5: Build lattice and add cells
    lattice = Lattice(
        name=f"{assembly_model.name}_{calculation_step.name}",
        center=center,
    )

    lattice = add_cells_to_regular_lattice(
        lattice, ordered_cells, pin_pitch, translation_x=translation_x, translation_y=translation_y
    )
    # Build optional water rods if present in the model and add to lattice
    if hasattr(assembly_model, "water_rods") and assembly_model.water_rods:
        lattice = create_and_add_water_rods_to_lattice(
            lattice, assembly_model,
            translation_x=translation_x, translation_y=translation_y,
            calculation_step=calculation_step,
        )
    
    if assembly_box_cell is not None:
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

