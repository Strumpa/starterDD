## Some helpful functions to deal with GLOW geometry building
# Author : R. Guasch
# Date : 04/02/2026
# Updated on 14/04/2026 to support part-length (vanished) rods.
# ----------------------------------------------------------------------------


from curses.ascii import ctrl
from glow import *
from glow.geometry_layouts.layouts import associate_colors_to_regions, build_compound_regions
from glow.support.types import *
from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.layouts import Region
from glow.geometry_layouts.lattices import CartesianLattice
from glow.geometry_layouts.geometries import Rectangle, Circle
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.main import TdtSetup, export_layout_to_tdt
from glow.interface.geom_interface import *
from glow.support.types import GeometryType, LayoutGeometryType, PropertyType, \
    SymmetryType
from glow.interface.geom_entities import wrap_shape
import numpy as np
import os
# Note: CartesianAssemblyModel and FuelPinModel are imported inside functions to avoid circular imports


SYM_TO_LAYOUT_AND_BOUNDARY = {
    SymmetryType.FULL: {"TISO": {"layout": LayoutGeometryType.ISOTROPIC, "boundary": None},
                        "TSPC": {"layout": LayoutGeometryType.RECTANGLE_SYM, "boundary": "AXIAL_SYMMETRY"}},
    
    SymmetryType.HALF: {"TISO": {"layout": LayoutGeometryType.SYMMETRIES_TWO, "boundary": "AXIAL_SYMMETRY"},
                        "TSPC": {"layout": LayoutGeometryType.RECTANGLE_SYM, "boundary": "AXIAL_SYMMETRY"}},
    
    SymmetryType.DIAG: {"TISO": {"layout": LayoutGeometryType.SYMMETRIES_TWO, "boundary": "AXIAL_SYMMETRY"},
                        "TSPC": {"layout": LayoutGeometryType.RECTANGLE_EIGHT, "boundary": "AXIAL_SYMMETRY"}},
    
    SymmetryType.QUARTER: {"TISO": {"layout": LayoutGeometryType.SYMMETRIES_TWO, "boundary": "AXIAL_SYMMETRY"},
                            "TSPC": {"layout": LayoutGeometryType.RECTANGLE_SYM, "boundary": "AXIAL_SYMMETRY"}},

    SymmetryType.EIGHTH: {"TISO": {"layout": LayoutGeometryType.SYMMETRIES_TWO, "boundary": "AXIAL_SYMMETRY"},
                            "TSPC": {"layout": LayoutGeometryType.RECTANGLE_EIGHT, "boundary": "AXIAL_SYMMETRY"}},
}

def make_grid_faces(parent: Rectangle, nx: int, ny: int):
    """
    Create a regular nx by ny grid of faces within the given parent rectangle, returning a list of the created faces.
    """
    lx, ly = parent.dimensions
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
    Generate CartesianCell objects for each individual subgeometry in the lattice.

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
    lattice_components : dictionary mapping (col_idx, row_idx) tuple to [CartesianCell, FuelPinModel]
        Allows to keep track of of the pin model describing each cell to ensure proper order of addition to the lattice in 
        ``add_cells_to_regular_lattice`` (non-generating cells are added first, then generating cells, to enforce correct mix numbering in SALOME).
    """
    # Import here to avoid circular import issues
    from ..DDModel.DragonModel import FuelPinModel
    
    lattice_components = {} # change to a dictionary to store cells and their pin model, at a given position. This will allow to keep track of the pin models associated with each cell.
    pitch = assemblyModel.pin_geometry_dict["pin_pitch"]
    
    row_idx = -1
    for row in assemblyModel.lattice:
        row_idx += 1 
        cell_idx = -1
        for pin in row:
            cell_idx += 1
            if isinstance(pin, FuelPinModel):
                fuel_material_mixtures = pin.fuel_material_mixtures
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
                    last_mat = "COOLANT"
                elif gap_radius is not None and gap_radius > fuel_radius and clad_radius is not None and clad_radius > fuel_radius:
                    # Fuel regions are inner most, then gap, then clad as outer most solid region, coolant is outside: fuel zones, gap, clad, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["GAP", "CLAD"]
                    last_mat = "COOLANT"
                elif gap_radius is None and clad_radius is not None and clad_radius > fuel_radius:
                    # No gap, clad is outer most solid region, so the order of materials from innermost to outermost is : fuel zones, clad, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["CLAD"]
                    last_mat = "COOLANT"
                elif gap_radius is not None and gap_radius > fuel_radius and clad_radius is None:
                    # No clad, gap is outer most solid region, so the order of materials from innermost to outermost is : fuel zones, gap, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures] + ["GAP"]
                    last_mat = "COOLANT"
                elif gap_radius is None and clad_radius is None:
                    # No gap, no clad, so only fuel zones and coolant, order of materials from innermost to outermost is : fuel zones, coolant
                    list_of_cell_mats = [fuel_mat.unique_material_mixture_name for fuel_mat in fuel_material_mixtures]
                    last_mat = "COOLANT"
                else:
                    raise ValueError(
                        f"Invalid combination of radii: fuel_radius={fuel_radius}, gap_radius={gap_radius}, clad_radius={clad_radius}"
                    )
                radii = pin.radii
                tmp_cell = CartesianCell(
                    name=pin.fuel_material_name,
                    width_height=(pitch, pitch),
                    center=(0.0, 0.0, 0.0),
                    base_props={PropertyType.MATERIAL:last_mat,
                                PropertyType.MACRO: f"MACRO{row_idx}{cell_idx}"}
                )
                for radius, mat in zip(radii[::-1],list_of_cell_mats[::-1]):
                    tmp_cell.add(
                        Region(Circle(radius=radius), properties={PropertyType.MATERIAL:mat, 
                                                                 PropertyType.MACRO:f"MACRO{row_idx}{cell_idx}"})
                    )
                            
                # Apply sectorization from calculation step if provided
                if calculation_step is not None:
                    sector_cfg = calculation_step.get_sectorization_for_pin(pin, isGd=pin.isGd)
                    if sector_cfg is not None:
                        print(f"Applying sectorization to cell at position ({cell_idx}, {row_idx}): sectors={sector_cfg.sectors}, angles={sector_cfg.angles}, windmill={sector_cfg.windmill}")
                        tmp_cell.sectorize(sector_cfg.sectors, sector_cfg.angles, windmill=sector_cfg.windmill)
                    else:
                        print(f"Warning: No sectorization config found for pin at position ({cell_idx}, {row_idx}). No sectorization will be applied to this cell.")
                # cell_idx in row : column number, ie position along the x direction,
                # row_idx in lattice: row number, ie position along the y direction
                lattice_components[(cell_idx, row_idx)] = [tmp_cell, pin] # store both the cell and its associated pin model for later reference
    return lattice_components


def add_cells_to_cartesian_lattice(lattice, lattice_components, cell_pitch, translation_x=0.0, translation_y=0.0):
    """
    Add fuel cells to the lattice, skipping water rod placeholders.
    Generating cells are added last in order to enforce order of mix attribution in SALOME (cells added last are assigned first mix numbers)

    Parameters
    ----------
    lattice : CartesianLattice
        The lattice to which cells will be added
    lattice_components : dict
        Dictionary mapping positions to lists of CartesianCell objects and their associated pin models
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
                lattice.add(
                    cell, position=((pos[0] + 0.5) * cell_pitch + translation_x,
                            (pos[1] + 0.5) * cell_pitch + translation_y,
                            0.0)
                )
    for pos, cell_and_pin in lattice_components.items():
        cell, pin = cell_and_pin
        if pin is not None and isinstance(pin, FuelPinModel):
            if pin.isGeneratingCell: # add generating cells last
                lattice.add(
                    cell, position=((pos[0] + 0.5) * cell_pitch + translation_x,
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
    tmp_cell = CartesianCell(
        name=water_rod_model.rod_ID,
        width_height=(bb, bb),
        center=center,
        base_props={PropertyType.MATERIAL: water_rod_model.coolant_material_name,
                    PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"}
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

    water_box_cell = Region(
        name=f"{water_rod_model.rod_ID}_box",
        geom_obj = outer_rect - inner_rect,
        properties={PropertyType.MATERIAL: water_rod_model.cladding_material_name,
                    PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"},
    )
    inner_moderator_cell = Region(
        name=f"{water_rod_model.rod_ID}_moderator",
        geom_obj = inner_rect,
        properties={PropertyType.MATERIAL: water_rod_model.moderator_material_name,
                    PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"},
    )

    tmp_cell.add(inner_moderator_cell)
    tmp_cell.add(water_box_cell)

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
            [tmp_cell],
            splitting_faces,
            shape_type=ShapeType.COMPOUND,
        )
        tmp_cell.geometry_maps[GeometryType.SECTORIZED] = \
        tmp_cell.get_geometry_map(GeometryType.SECTORIZED) // wrap_shape(re_partitioned)
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
            tmp_cell = CartesianCell(
                name=water_rod_model.rod_ID,
                width_height=(
                    water_rod_model.bounding_box_side_length,
                    water_rod_model.bounding_box_side_length,
                ),
                center=(0.0, 0.0, 0.0),
                base_props={PropertyType.MATERIAL: water_rod_model.coolant_material_name,
                            PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"},
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

            tmp_cell.add(Region(Circle(radius=water_rod_model.outer_radius), properties={PropertyType.MATERIAL: water_rod_model.cladding_material_name,
                                                                 PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"}))
            tmp_cell.add(Region(Circle(radius=water_rod_model.inner_radius), properties={PropertyType.MATERIAL: water_rod_model.moderator_material_name,
                                                                 PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"}))
            # Add circles: extra moderator sub-rings, then inner, then outer
            for r in extra_radii[::-1]:  # add extra moderator radii from outermost to innermost
                tmp_cell.add(Region(Circle(radius=r), properties={PropertyType.MATERIAL: water_rod_model.moderator_material_name,
                                                 PropertyType.MACRO: f"MACRO_{water_rod_model.rod_ID}"}))



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
            else:
                tmp_cell.sectorize([1,1,1], [0,0,0], windmill=False)
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
        lattice.add(
            tmp_cell,
            position=(cx, cy, 0.0),
        )

    return lattice

def add_vanished_rods_to_lattice(lattice, assembly_model, translation_x=0.0, translation_y=0.0, calculation_step=None):
    """
    Create vanished rod cells from the assembly model and add them to the lattice at their centers.

    Parameters:
    -----------
    lattice : Lattice
        The lattice to which vanished rod cells will be added
    assembly_model : CartesianAssemblyModel
        The assembly model containing a ``vanished_rods`` list of ``VanishedRodModel``
        instances. This function uses each model's ``rod_ID``, ``center``, and
        ``default_sectorization_radius`` attributes; lattice indices may also be
        present on the model but are not used here. Vanished rod centers are
        expected to already be in assembly/YAML coordinates.
    translation_x : float
        Unused. Vanished rod centers are already in assembly coordinates if translation_offset_x has been defined 
        from YAML input geometry.
        Kept for function signature consistency with pin positioning.
    translation_y : float
        Unused. Vanished rod centers are already in assembly coordinates if translation_offset_x has been defined 
        from YAML input geometry.
        Kept for function signature consistency with pin positioning.
    calculation_step : CalculationStep or None
        CalculationStep object to retrieve vanished rod sectorization options from.

    If translation offsets are not provided, vanished rods will be positioned according to their position in the lattice.

    """
    from ..DDModel.DragonModel import VanishedRodModel

    lattice_pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]

    for rod_model in assembly_model.vanished_rods:
        if not isinstance(rod_model, VanishedRodModel):
            raise ValueError(
                f"Expected VanishedRodModel in assembly_model.vanished_rods, "
                f"but got {type(rod_model)}"
            )
        # RectCell replaced by CartesianCell with region and circular subregions :

        tmp_cell = CartesianCell(
            name=rod_model.rod_ID,
            width_height=(lattice_pin_pitch, lattice_pin_pitch),
            center=(0.0, 0.0, 0.0),
            base_props={PropertyType.MATERIAL:"COOLANT",
                        PropertyType.MACRO: f"MACRO_{rod_model.rod_ID}"}
        )

        n_regions = 1 # default number of regions if no sectorization provided
        
        if calculation_step is not None:
            vr_sector_cfg = calculation_step.get_vanished_rod_sectorization()
            if vr_sector_cfg is not None:
                if vr_sector_cfg.base_radius is not None:
                    vr_sector_cfg.resolve_radii_and_sectors()
                else:
                    # If no base_radius provided, use the default sectorization radius from the rod model (set to be equal to the cladding radius)
                    vr_sector_cfg.resolve_radii_and_sectors(rod_model.default_sectorization_radius)
                radii = vr_sector_cfg.radial_split_points
                for radius in radii:
                    tmp_cell.add(Region(Circle(radius=radius, properties={PropertyType.MATERIAL:"COOLANT",
                                                              PropertyType.MACRO:f"MACRO_{rod_model.rod_ID}"})))
                if vr_sector_cfg.sector_config:
                    tmp_cell.sectorize(vr_sector_cfg.sector_config.sectors, vr_sector_cfg.sector_config.angles, windmill=vr_sector_cfg.windmill)
                n_regions = len(radii) + 1 # number of regions is number of circles + 1 (the central region inside the innermost circle)

        tmp_cell.set_properties({
            PropertyType.MATERIAL: ["COOLANT"]*n_regions,
            PropertyType.MACRO: [f"MACRO_{rod_model.rod_ID}"]*n_regions,
        })

        if rod_model.center is None:
            # compute center from lattice indices if not provided
            posx, posy = rod_model.x_index, rod_model.y_index
            cx = (posx + 0.5) * lattice_pin_pitch
            cy = (posy + 0.5) * lattice_pin_pitch
        else:
            cx, cy = rod_model.center
        lattice.add(
            tmp_cell,
            position=(cx, cy, 0.0)
        )

    return lattice


def export_glow_geom(output_path, output_file_name, assembly_universe, symmetry_type, tracking_option, export_macro=False):
    """
    Export the geometry of the lattice to a TDT file for GLOW simulation.

    Parameters
    ----------
    output_path : str
        Path to save the exported TDT file
    output_file_name : str
        Name of the exported TDT file
    assembly_universe : CartesianCell filled with lattice and box components
        The CSG assembly universe whose geometry is to be exported
    symmetry_type : SymmetryType
        Symmetry type of the geometry, e.g. ``SymmetryType.FULL`` or ``SymmetryType.QUARTER``.  This is used to determine the geometry type to
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
        property_to_show = PropertyType.MACRO
        geometry_type_to_show = GeometryType.SECTORIZED
        print(f"Attempting to show macro properties on SECTORIZED geometry for export...")
        try:
            assembly_universe.show(property_to_show, geometry_type_to_show)
        except RuntimeError as e:
            print(f"Error occurred while showing macro properties on SECTORIZED geometry: {e}, trying to show properties on TECHNOLOGICAL geometry instead")
            geometry_type_to_show = GeometryType.TECHNOLOGICAL
            assembly_universe.show(property_to_show, geometry_type_to_show)
    else:
        properties_to_export = [PropertyType.MATERIAL]
        output_file_name = f"{output_file_name}_{tracking_option}"
        property_to_show = PropertyType.MATERIAL
        geometry_type_to_show = GeometryType.SECTORIZED
        try:
            assembly_universe.show(property_to_show, geometry_type_to_show)
        except RuntimeError as e:
            print(f"Error occurred while showing material properties on SECTORIZED geometry: {e}, trying to show properties on TECHNOLOGICAL geometry instead")
            geometry_type_to_show = GeometryType.TECHNOLOGICAL
            assembly_universe.show(property_to_show, geometry_type_to_show)

    full_tdt_path = os.path.join(output_path, output_file_name)

    geometry_type_to_export = geometry_type_to_show

    if tracking_option == "TISO":
        
        layout_geometry_type = SYM_TO_LAYOUT_AND_BOUNDARY[symmetry_type]["TISO"]["layout"]
        
        export_layout_to_tdt(
            assembly_universe, full_tdt_path, TdtSetup(geometry_type_to_export,
                                             property_types=properties_to_export,
                                             type_geo=layout_geometry_type,
                                             symmetry_type=symmetry_type))
    elif tracking_option == "TSPC":  
        
        layout_geometry_type = SYM_TO_LAYOUT_AND_BOUNDARY[symmetry_type]["TSPC"]["layout"]
        
        export_layout_to_tdt(
            assembly_universe, full_tdt_path, TdtSetup(geometry_type_to_export,
                                             property_types=properties_to_export,
                                             type_geo=layout_geometry_type,
                                             symmetry_type=symmetry_type))


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
    # Mapping: NW is identity.  Mirror in x flips left <-> right (0 <-> 1, 3 <-> 2).
    # Mirror in y flips top <-> bottom (0 <-> 3, 1 <-> 2).
    _mirror_x = {0: 1, 1: 0, 2: 3, 3: 2}
    _mirror_y = {0: 3, 1: 2, 2: 1, 3: 0}

    result = list(corner_indices)
    if cross_corner in ("north-east", "south-east"):
        result = [(_mirror_x[idx], r) for idx, r in result]
    if cross_corner in ("south-west", "south-east"):
        result = [(_mirror_y[idx], r) for idx, r in result]
    return result


def _build_control_cross_elements(ctrl, ap, assembly_center=(None, None, None)):
    """
    Build a control cross geometry in CSG appraoch using glow v1.1.0
    
    Parameters
    ----------
    ctrl : ControlCrossModel
        The control cross model with all geometric dimensions.
    ap : float
        Assembly pitch.
    assembly_center : tuple, optional
        (x, y, z) center of the assembly_universe. Defaults to (ap/2, ap/2, 0).
        Region centers are converted from absolute assembly coords to local
        (center-relative) coords.


    Returns 
    -------
    control_cross_universe : CartesianCell
        A CartesianCell containing all the Regions that define the control cross geometry, centered at 
        either the NW or SW corner depending on the identified lattice symmetry.

    """

    # retrieve dimensions and properties of the control cross model

    if assembly_center[0] is None:
        assembly_center = (ap / 2.0, ap / 2.0, 0.0)

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
    extra_moderator_gap = (delta - 2.0*r_outer) / 2.0
    half_square_side = r_outer + extra_moderator_gap
    distance_from_first_center_to_last_center = (n_tubes - 1) * 2.0 * half_square_side
    width_controlled_section = bhs - st - cshs
    print(f"Computed width of controlled section: {width_controlled_section} cm, distance from first center to last center: {distance_from_first_center_to_last_center} cm")

    extra_distance_from_last_center_to_edge = (width_controlled_section - distance_from_first_center_to_last_center) / 2.0

    first_tube_offset = cshs + r_outer + extra_moderator_gap # In assembly coordinates.

    print(f"Computed first tube offset: {first_tube_offset}, distance from first center to last center: {distance_from_first_center_to_last_center}, extra moderator gap: {extra_moderator_gap}, extra distance from last center to edge: {extra_distance_from_last_center_to_edge}")
    print(f"Check: computed first tube offset : {first_tube_offset} cm, provided first tube offset: {first_offset} cm, difference: {first_tube_offset - first_offset} cm")
    print(f"First offset: {first_offset} cm, expected first offset based on tube dimensions and extra moderator gap: {cshs + r_outer + extra_moderator_gap} cm, difference: {first_offset - (cshs + r_outer + extra_moderator_gap)} cm")

    print(f"Check provided delta: {delta} cm, expected delta based on tube dimensions and extra moderator gap: {2.0 * half_square_side} cm, difference: {delta - 2.0 * half_square_side} cm")

    is_solid = ctrl.is_solid
    absorber_mat = ctrl.absorber_material
    sheath_mat = ctrl.sheath_material

    assembly_bounding_rect = Rectangle(
        height=ap,
        width=ap,
        center=(ap/2,ap/2,0.0)
    )

    elements = {}
    elements["CTRL_H_SHEATH"] = {}
    elements["CTRL_V_SHEATH"] = {}
    # Create a universe cell for the control cross and add all regions
    control_cross_horizontal_wing = Rectangle(
        name="CTRL_CROSS_H",
        width=bhs,
        height=bt,
        center=(bhs/2.0, ap, 0.0)
        #rounded_corners=_remap_rounded_corner_indices([(1, tr), (2, tr)], corner) if tr > 0.0 else None
    )
    control_cross_horizontal_wing_rounded = Rectangle(
        name="CTRL_CROSS_H_rounded",
        width=bhs,
        height=bt,
        center=(bhs/2.0, ap, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(1, tr), (2, tr)], corner) if tr > 0.0 else None
    )

    # compute delta between centers of horizontal wing with rounded corners and horizontal wing without rounded corners
    delta_centers_x_horizontal_wing = get_point_coordinates(control_cross_horizontal_wing_rounded.o)[0] - get_point_coordinates(control_cross_horizontal_wing.o)[0]
    delta_centers_y_horizontal_wing = get_point_coordinates(control_cross_horizontal_wing_rounded.o)[1] - get_point_coordinates(control_cross_horizontal_wing.o)[1]
    #wing_elements[(bhs/2.0 + delta_centers_x_horizontal_wing, ap + delta_centers_y_horizontal_wing, 0.0)] = control_cross_horizontal_wing_rounded
    
    h_wing_in_assembly_footprint = control_cross_horizontal_wing_rounded * assembly_bounding_rect

    elements["CTRL_H_SHEATH"]["geometry"] = h_wing_in_assembly_footprint
    elements["CTRL_H_SHEATH"]["offset"] = (delta_centers_x_horizontal_wing, delta_centers_y_horizontal_wing, 0.0)
    elements["CTRL_H_SHEATH"]["macro"] = "CTRL_H"
    elements["CTRL_H_SHEATH"]["material"] = sheath_mat

    control_cross_vertical_wing = Rectangle(
        name="CTRL_CROSS_V",
        width=bt,
        height=bhs,
        center=(0.0, ap - bhs/2.0, 0.0)
    )
    control_cross_vertical_wing_rounded = Rectangle(
        name="CTRL_CROSS_V_rounded",
        width=bt,
        height=bhs,
        center=(0.0, ap - bhs/2.0, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(0, tr), (1, tr)], corner) if tr > 0.0 else None
    )

    delta_centers_x_vertical_wing = get_point_coordinates(control_cross_vertical_wing_rounded.o)[0] - get_point_coordinates(control_cross_vertical_wing.o)[0]
    delta_centers_y_vertical_wing = get_point_coordinates(control_cross_vertical_wing_rounded.o)[1] - get_point_coordinates(control_cross_vertical_wing.o)[1]
    
    v_wing_in_assembly_footprint = control_cross_vertical_wing_rounded * assembly_bounding_rect
    elements["CTRL_V_SHEATH"]["geometry"] = v_wing_in_assembly_footprint
    elements["CTRL_V_SHEATH"]["offset"] = (delta_centers_x_vertical_wing, delta_centers_y_vertical_wing, 0.0)
    elements["CTRL_V_SHEATH"]["macro"] = "CTRL_V"
    elements["CTRL_V_SHEATH"]["material"] = sheath_mat

    # build hollow sheath region and absorber pins 
    if not is_solid:
        base_material = "MODERATOR"
        elements["CTRL_H_HOLLOW"] = {}
        elements["CTRL_V_HOLLOW"] = {}
    # build inner region for hollow sheath (moderator-filled cavity inside the blade)
        inner_sheath_horizontal_wing = Rectangle(
            name="CTRL_CROSS_H_INNER",
            width=(bhs - st - cshs),
            height=(bt - 2.0*st),
            center=(cshs + (bhs - st - cshs)/2.0, ap, 0.0)
        )
        # get center for the horizontal wing inner sheath without rounded corners :
        inner_sheath_horizontal_wing_rounded = Rectangle(
            name="CTRL_CROSS_H_INNER_HOLLOW",
            width=(bhs - st - cshs),
            height=(bt - 2.0*st),
            center=(cshs + (bhs - st - cshs)/2.0, ap, 0.0),
            rounded_corners=_remap_rounded_corner_indices([(1, tr - st), (2, tr - st)], corner) if tr > 0.0 else None
        )
        delta_centers_x_inner_sheath = get_point_coordinates(inner_sheath_horizontal_wing_rounded.o)[0] - get_point_coordinates(inner_sheath_horizontal_wing.o)[0]
        delta_centers_y_inner_sheath = get_point_coordinates(inner_sheath_horizontal_wing_rounded.o)[1] - get_point_coordinates(inner_sheath_horizontal_wing.o)[1]
        
        intersection_with_assembly_footprint = inner_sheath_horizontal_wing_rounded * assembly_bounding_rect
        elements["CTRL_H_HOLLOW"]["geometry"] = intersection_with_assembly_footprint
        elements["CTRL_H_HOLLOW"]["offset"] = (delta_centers_x_inner_sheath, delta_centers_y_inner_sheath, 0.0)
        elements["CTRL_H_HOLLOW"]["macro"] = "CTRL_H"
        elements["CTRL_H_HOLLOW"]["material"] = base_material
    

        # inner region for hollow sheath in vertical wing :
        inner_sheath_vertical_wing = Rectangle(
            name="CTRL_CROSS_V_INNER",
            width=(bt - 2.0*st),
            height=(bhs - st - cshs),
            center=(0.0, ap - (cshs + (bhs - st - cshs)/2.0), 0.0)
        )
        inner_sheath_vertical_wing_rounded = Rectangle(
            name="CTRL_CROSS_V_INNER_HOLLOW",
            width=(bt - 2.0*st),
            height=(bhs - st - cshs),
            center=(0.0, ap - (cshs + (bhs - st - cshs)/2.0), 0.0),
            rounded_corners=_remap_rounded_corner_indices([(0, tr - st), (1, tr - st)], corner) if tr > 0.0 else None
        )

        delta_centers_x_inner_sheath_vertical = get_point_coordinates(inner_sheath_vertical_wing_rounded.o)[0] - get_point_coordinates(inner_sheath_vertical_wing.o)[0]
        delta_centers_y_inner_sheath_vertical = get_point_coordinates(inner_sheath_vertical_wing_rounded.o)[1] - get_point_coordinates(inner_sheath_vertical_wing.o)[1]
        
        intersection_with_assembly_footprint = inner_sheath_vertical_wing_rounded * assembly_bounding_rect
        elements["CTRL_V_HOLLOW"]["geometry"] = intersection_with_assembly_footprint
        elements["CTRL_V_HOLLOW"]["offset"] = (delta_centers_x_inner_sheath_vertical, delta_centers_y_inner_sheath_vertical, 0.0)
        elements["CTRL_V_HOLLOW"]["macro"] = "CTRL_V"
        elements["CTRL_V_HOLLOW"]["material"] = base_material
    
    for i in range(n_tubes):
        elements[f"ABS_TUBE_H_{i}"] = {}
        elements[f"ABS_TUBE_V_{i}"] = {}
        offset = first_tube_offset + float(i * delta)
        # Horizontal wing tube (along x)
        tx_h, ty_h = _corner_transform(corner, offset, ap, ap)
        tx_v, ty_v = _corner_transform(corner, 0.0, ap - offset, ap)
        if is_solid:
            # Solid rod (AT10 style): single absorber circle
            abs_tube_h = Circle(radius=r_outer, center=(tx_h, ty_h, 0.0))
            h_tube_intersection_with_af = abs_tube_h * assembly_bounding_rect
            elements[f"ABS_TUBE_H_{i}"]["geometry"] = h_tube_intersection_with_af
            elements[f"ABS_TUBE_H_{i}"]["macro"] = "CTRL_H"
            elements[f"ABS_TUBE_H_{i}"]["material"] = absorber_mat

            abs_tube_v = Circle(radius=r_outer, center=(tx_v, ty_v, 0.0))
            v_tube_intersection_with_af = abs_tube_v * assembly_bounding_rect
            elements[f"ABS_TUBE_V_{i}"]["geometry"] = v_tube_intersection_with_af
            elements[f"ABS_TUBE_V_{i}"]["macro"] = "CTRL_V"
            elements[f"ABS_TUBE_V_{i}"]["material"] = absorber_mat
            

        else:
            # Hollow tube (GE-14 style): inner absorber + outer cladding
            abs_tube_h = Circle(radius=r_inner, center=(tx_h, ty_h, 0.0))
            h_tube_intersection_with_af = abs_tube_h * assembly_bounding_rect
            elements[f"ABS_TUBE_H_{i}"]["geometry"] = h_tube_intersection_with_af
            elements[f"ABS_TUBE_H_{i}"]["macro"] = "CTRL_H"
            elements[f"ABS_TUBE_H_{i}"]["material"] = absorber_mat

            abs_tube_v = Circle(radius=r_inner, center=(tx_v, ty_v, 0.0))
            v_tube_intersection_with_af = abs_tube_v * assembly_bounding_rect
            elements[f"ABS_TUBE_V_{i}"]["geometry"] = v_tube_intersection_with_af
            elements[f"ABS_TUBE_V_{i}"]["macro"] = "CTRL_V"
            elements[f"ABS_TUBE_V_{i}"]["material"] = absorber_mat

            elements[f"SHEATH_TUBE_H_{i}"] = {}
            elements[f"SHEATH_TUBE_V_{i}"] = {}
            sheath_tube_h = Circle(radius=r_outer, center=(tx_h, ty_h, 0.0))
            h_tube_intersection_with_af = sheath_tube_h * assembly_bounding_rect - abs_tube_h
            elements[f"SHEATH_TUBE_H_{i}"]["geometry"] = h_tube_intersection_with_af
            elements[f"SHEATH_TUBE_H_{i}"]["macro"] = "CTRL_H"
            elements[f"SHEATH_TUBE_H_{i}"]["material"] = sheath_mat

            sheath_tube_v = Circle(radius=r_outer, center=(tx_v, ty_v, 0.0))
            v_tube_intersection_with_af = sheath_tube_v * assembly_bounding_rect - abs_tube_v
            elements[f"SHEATH_TUBE_V_{i}"]["geometry"] = v_tube_intersection_with_af
            elements[f"SHEATH_TUBE_V_{i}"]["macro"] = "CTRL_V"
            elements[f"SHEATH_TUBE_V_{i}"]["material"] = sheath_mat
    
    return elements

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
    from glow.geometry_layouts.layouts import Region

    ap = assembly_model.assembly_pitch
    cbt = assembly_model.channel_box_thickness
    gap_wide = assembly_model.gap_wide
    gap_narrow = assembly_model.gap_narrow
    corner_r_inner = assembly_model.corner_inner_radius_of_curvature
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_rows = len(assembly_model.lattice_description)
    n_cols = len(assembly_model.lattice_description[0])

    lattice_pitch_x = n_cols * pin_pitch
    lattice_pitch_y = n_rows * pin_pitch

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
    print(f"Created coolant rectangle with width {channel_box_inner_x:.3f} and height {channel_box_inner_y:.3f} at center {rect_center} with rounded corners {rounded_corners_coolant}")

    channel_box_rect = Rectangle(
        name="channel_box",
        height=channel_box_outer_y,
        width=channel_box_outer_x,
        center=rect_center,
        rounded_corners=rounded_corners_chanbox,
    )
    print(f"Created channel box rectangle with width {channel_box_outer_x:.3f} and height {channel_box_outer_y:.3f} at center {rect_center} with rounded corners {rounded_corners_chanbox}")
    
    if lattice_pitch_x < channel_box_inner_x and lattice_pitch_y < channel_box_inner_y:
        lattice_rect = Rectangle(
            name="lattice_bounding_rectangular_box",
            height=lattice_pitch_y,
            width=lattice_pitch_x,
            center=rect_center
        )
        print(f"Created lattice bounding rectangle with width {lattice_pitch_x:.3f} and height {lattice_pitch_y:.3f} at center {rect_center}")
    else:
        lattice_rect = None
        print(f"Warning: Lattice pitch ({lattice_pitch_x:.3f}, {lattice_pitch_y:.3f}) is greater or equal to channel box dimensions ({channel_box_inner_x:.3f}, {channel_box_inner_y:.3f}). Skipping lattice bounding rectangle.")

    return  lattice_rect, coolant_rect, channel_box_rect


# ==============================================================================
# PHASE 1: Generate MACRO Subdivision Rectangles
# ==============================================================================

def _generate_macro_subdivision_rectangles(
    assembly_model,
    x0, y0, x1, y1,
    center
):
    """
    Generate rectangles for all MACRO regions surrounding the pin lattice.

    Handles symmetric and asymmetric gap configurations. The lattice footprint
    (x0, y0, x1, y1) already incorporates asymmetric offsets from translation_offset_x/y.

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model with pin_pitch, lattice_description, etc.
    x0, y0, x1, y1 : float
        Lattice footprint bounds (already include asymmetric gap handling)
    center : tuple
        Assembly center (x_center, y_center, z_center)

    Returns
    -------
    Dict[str, tuple]
        Maps MACRO name → (name, height, width, (x, y, z), rounded_corners)

        Example:
        {
            "LEFT_1": ("LEFT_1", strip_width, lattice_pitch_y, (x, y, 0), False),
            "BOT_1": ("BOT_1", strip_height, pin_pitch, (x, y, 0), False),
            "CORNER_BL": ("CORNER_BL", y0, x0, (x0/2, y0/2, 0), False),
            ...
        }
    """
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_rows = len(assembly_model.lattice_description)
    n_cols = len(assembly_model.lattice_description[0])

    macros = {}
    z_center = center[2]

    # ======================================================================
    # BOTTOM STRIP: [x0, x1] × [0, y0]
    # Subdivided into n_cols regions (one per pin column)
    # ======================================================================
    strip_height = y0  # Distance from bottom edge (y=0) to lattice bottom (y=y0)

    for col_idx in range(n_cols):
        # Rectangle spans one pin column width
        col_x_min = x0 + col_idx * pin_pitch
        col_x_max = col_x_min + pin_pitch

        # Center of this bottom strip rectangle
        rect_center_x = (col_x_min + col_x_max) / 2.0
        rect_center_y = strip_height / 2.0  # Centered in the bottom strip

        macro_name = f"BOT_{col_idx + 1}"  # 1-based indexing
        print(f"Create rect_params for {macro_name} with center ({rect_center_x:.3f}, {rect_center_y:.3f}, {z_center:.3f}), width {pin_pitch:.3f}, height {strip_height:.3f}")
        rect_params = (
            macro_name,                                  # name
            strip_height,                                # height (y-direction)
            pin_pitch,                                   # width (x-direction)
            (rect_center_x, rect_center_y, z_center),   # center (x, y, z)
            None                                         # rounded_corners (None = no rounding)
        )
        macros[macro_name] = rect_params

    # ======================================================================
    # TOP STRIP: [x0, x1] × [y1, ap]
    # Subdivided into n_cols regions (one per pin column)
    # ======================================================================
    strip_height = ap - y1  # Distance from lattice top (y=y1) to top edge (y=ap)

    for col_idx in range(n_cols):
        col_x_min = x0 + col_idx * pin_pitch
        col_x_max = col_x_min + pin_pitch

        rect_center_x = (col_x_min + col_x_max) / 2.0
        rect_center_y = y1 + (ap - y1) / 2.0  # Centered in the top strip

        macro_name = f"TOP_{col_idx + 1}"
        print(f"Create rect_params for {macro_name} with center ({rect_center_x:.3f}, {rect_center_y:.3f}, {z_center:.3f}), width {pin_pitch:.3f}, height {strip_height:.3f}")
        rect_params = (
            macro_name,
            strip_height,
            pin_pitch,
            (rect_center_x, rect_center_y, z_center),
            None                                         # rounded_corners (None = no rounding)
        )
        macros[macro_name] = rect_params

    # ======================================================================
    # LEFT STRIP: [0, x0] × [y0, y1]
    # Subdivided into n_rows regions (one per pin row)
    # ======================================================================
    strip_width = x0  # Distance from left edge (x=0) to lattice left (x=x0)
    lattice_pitch_y = y1 - y0

    for row_idx in range(n_rows):
        row_y_min = y0 + row_idx * pin_pitch
        row_y_max = row_y_min + pin_pitch

        rect_center_x = strip_width / 2.0  # Centered in the left strip
        rect_center_y = (row_y_min + row_y_max) / 2.0

        macro_name = f"LEFT_{row_idx + 1}"
        print(f"Create rect_params for {macro_name} with center ({rect_center_x:.3f}, {rect_center_y:.3f}, {z_center:.3f}), width {pin_pitch:.3f}, height {strip_width:.3f}")
        rect_params = (
            macro_name,
            pin_pitch,                                   # height (y-direction)
            strip_width,                                 # width (x-direction)
            (rect_center_x, rect_center_y, z_center),
            None                                         # rounded_corners (None = no rounding)
        )
        macros[macro_name] = rect_params

    # ======================================================================
    # RIGHT STRIP: [x1, ap] × [y0, y1]
    # Subdivided into n_rows regions (one per pin row)
    # ======================================================================
    strip_width = ap - x1  # Distance from lattice right (x=x1) to right edge (x=ap)

    for row_idx in range(n_rows):
        row_y_min = y0 + row_idx * pin_pitch
        row_y_max = row_y_min + pin_pitch

        rect_center_x = x1 + (ap - x1) / 2.0  # Centered in the right strip
        rect_center_y = (row_y_min + row_y_max) / 2.0

        macro_name = f"RIGHT_{row_idx + 1}"
        print(f"Create rect_params for {macro_name} with center ({rect_center_x:.3f}, {rect_center_y:.3f}, {z_center:.3f}), width {pin_pitch:.3f}, height {strip_width:.3f}")
        rect_params = (
            macro_name,
            pin_pitch,
            strip_width,
            (rect_center_x, rect_center_y, z_center),
            None                                         # rounded_corners (None = no rounding)
        )
        macros[macro_name] = rect_params

    # ======================================================================
    # CORNERS
    # ======================================================================
    corner_width_bl = x0
    corner_height_bl = y0

    # Bottom-Left
    macros["CORNER_BL"] = (
        "CORNER_BL",
        corner_height_bl,
        corner_width_bl,
        (corner_width_bl / 2.0, corner_height_bl / 2.0, z_center),
        None                                             # rounded_corners (None = no rounding)
    )

    # Bottom-Right
    corner_width_br = ap - x1
    corner_height_br = y0
    macros["CORNER_BR"] = (
        "CORNER_BR",
        corner_height_br,
        corner_width_br,
        (x1 + corner_width_br / 2.0, corner_height_br / 2.0, z_center),
        None                                             # rounded_corners (None = no rounding)
    )

    # Top-Left
    corner_width_tl = x0
    corner_height_tl = ap - y1
    macros["CORNER_TL"] = (
        "CORNER_TL",
        corner_height_tl,
        corner_width_tl,
        (corner_width_tl / 2.0, y1 + corner_height_tl / 2.0, z_center),
        None                                             # rounded_corners (None = no rounding)
    )

    # Top-Right
    corner_width_tr = ap - x1
    corner_height_tr = ap - y1
    macros["CORNER_TR"] = (
        "CORNER_TR",
        corner_height_tr,
        corner_width_tr,
        (x1 + corner_width_tr / 2.0, y1 + corner_height_tr / 2.0, z_center),
        None                                             # rounded_corners (None = no rounding)
    )

    return macros


# ==============================================================================
# PHASE 1B: Generate Layer-Aware MACRO Region Subdivisions
# ==============================================================================

def is_degenerate(geometry):
    """
    Check if a geometry is degenerate (zero area, empty, or None).

    Parameters
    ----------
    geometry : any
        Geometry object (Rectangle, Circle, or compound geometry)

    Returns
    -------
    bool
        True if geometry is None, has zero area, or is empty
    """
    if geometry is None:
        return True

    # Check for zero-area geometries
    # Rectangles have dimensions property
    if hasattr(geometry, 'dimensions'):
        dims = geometry.dimensions
        if dims[0] < 1e-10 or dims[1] < 1e-10:
            return True

    # Check for area attribute
    if hasattr(geometry, 'area'):
        if geometry.area < 1e-10:
            return True

    return False


def _wrap_and_validate_geometry(geometry, layer_name=""):
    """
    Wrap geometry to ensure it's valid for Region creation.
    Handles COMPOUND geometries by wrapping with wrap_shape.

    Parameters
    ----------
    geometry : any
        Geometry object from set operations
    layer_name : str
        Name for logging purposes

    Returns
    -------
    geometry or None
        Wrapped geometry or None if invalid
    """
    if geometry is None:
        return None

    try:
        # Try wrapping the geometry to ensure it's a valid FACE
        wrapped = wrap_shape(geometry)
        return wrapped
    except Exception as e:
        print(f"    Warning: Could not wrap {layer_name} geometry: {e}")
        return None


def _generate_macro_region_layers(macro_rect, macro_name, coolant_rect, channel_box_rect):
    """
    Decompose a MACRO rectangle into layer-specific geometries using set operations.

    For each MACRO rectangle, this function computes the intersection with layer
    boundaries to produce layer-specific geometries. The layers are:
    - COOLANT: Inside coolant_rect
    - CHANNEL_BOX: Between coolant_rect and channel_box_rect
    - MODERATOR: Outside channel_box_rect

    Parameters
    ----------
    macro_rect : Rectangle
        Rectangle geometry for the MACRO region
    macro_name : str
        Identifier for the MACRO (e.g., "BOT_1", "CORNER_BL")
    coolant_rect : Rectangle
        Reference boundary for inner coolant region
    channel_box_rect : Rectangle
        Reference boundary for outer channel box region

    Returns
    -------
    dict
        Maps layer type to geometry:
        {
            'COOLANT': geometry or None,
            'CHANNEL_BOX': geometry or None,
            'MODERATOR': geometry or None
        }

        Each geometry may be:
        - Rectangle: if layer aligns with MACRO bounds
        - Wrapped compound geometry: if MACRO spans layer boundary
        - None: if layer doesn't exist for this MACRO
    """
    result = {}

    # Layer 1: COOLANT = macro_rect ∩ coolant_rect (intersection)
    try:
        coolant_geom = macro_rect * coolant_rect  # intersection operator
        # Wrap geometry to handle COMPOUND types
        coolant_geom = _wrap_and_validate_geometry(coolant_geom, f"{macro_name}/COOLANT")
        if coolant_geom is not None and not is_degenerate(coolant_geom):
            result['COOLANT'] = coolant_geom
            print(f"  {macro_name}: COOLANT geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute COOLANT layer for {macro_name}: {e}")

    # Layer 2: CHANNEL_BOX = (macro_rect ∩ channel_box_rect) - coolant_rect (intersection minus)
    try:
        temp_geom = macro_rect * channel_box_rect  # intersection
        channel_box_geom = temp_geom - coolant_rect  # difference
        # Wrap geometry to handle COMPOUND types
        channel_box_geom = _wrap_and_validate_geometry(channel_box_geom, f"{macro_name}/CHANNEL_BOX")
        if channel_box_geom is not None and not is_degenerate(channel_box_geom):
            result['CHANNEL_BOX'] = channel_box_geom
            print(f"  {macro_name}: CHANNEL_BOX geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute CHANNEL_BOX layer for {macro_name}: {e}")

    # Layer 3: MODERATOR = macro_rect - channel_box_rect (difference)
    try:
        moderator_geom = macro_rect - channel_box_rect  # difference operator
        # Wrap geometry to handle COMPOUND types
        moderator_geom = _wrap_and_validate_geometry(moderator_geom, f"{macro_name}/MODERATOR")
        if moderator_geom is not None and not is_degenerate(moderator_geom):
            result['MODERATOR'] = moderator_geom
            print(f"  {macro_name}: MODERATOR geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute MODERATOR layer for {macro_name}: {e}")

    return result


def _generate_corner_channel_box_regions(x0, y0, x1, y1, ap, n_rows, n_cols, coolant_rect, channel_box_rect, center):
    """
    Generate channel box material regions at assembly corners that intersect with lattice footprint.

    These are the 4 corner regions where:
    - The channel box extends into the lattice footprint
    - They are OUTSIDE coolant_rect, INSIDE channel_box_rect, INSIDE lattice footprint
    - They should be assigned to the corner fuel pin macros (MACRO00, MACRO0{n_cols-1}, etc.)

    Parameters
    ----------
    x0, y0, x1, y1 : float
        Lattice footprint bounds
    ap : float
        Assembly pitch
    n_rows, n_cols : int
        Number of rows and columns in lattice
    coolant_rect : Rectangle
        Reference boundary for inner coolant region
    channel_box_rect : Rectangle
        Reference boundary for outer channel box region
    center : tuple
        Assembly center (x_center, y_center, z_center)

    Returns
    -------
    dict
        Maps corner identifier to dict with geometry and macro name:
        {
            'BL': {'geometry': Rectangle or None, 'macro': 'MACRO00'},
            'BR': {'geometry': Rectangle or None, 'macro': 'MACRO0{n_cols-1}'},
            'TL': {'geometry': Rectangle or None, 'macro': 'MACRO{n_rows-1}0'},
            'TR': {'geometry': Rectangle or None, 'macro': 'MACRO{n_rows-1}{n_cols-1}'}
        }
    """
    z_center = center[2]
    result = {}
    pin_pitch_x = (x1 - x0) / n_cols  # Assuming uniform pin pitch in x
    pin_pitch_y = (y1 - y0) / n_rows  # Assuming uniform pin pitch in y
    # Define the 4 corner regions in assembly space
    corners = {
        'BL': {'x_min': x0, 'x_max': x0 + pin_pitch_x, 'y_min': y0, 'y_max': y0 + pin_pitch_y, 'macro_row': 0, 'macro_col': 0},
        'BR': {'x_min': x1 - pin_pitch_x, 'x_max': x1, 'y_min': y0, 'y_max': y0 + pin_pitch_y, 'macro_row': 0, 'macro_col': n_cols - 1},
        'TL': {'x_min': x0, 'x_max': x0 + pin_pitch_x, 'y_min': y1 - pin_pitch_y, 'y_max': y1, 'macro_row': n_rows - 1, 'macro_col': 0},
        'TR': {'x_min': x1 - pin_pitch_x, 'x_max': x1, 'y_min': y1 - pin_pitch_y, 'y_max': y1, 'macro_row': n_rows - 1, 'macro_col': n_cols - 1}
    }

    for corner_id, corner_def in corners.items():
        x_min, x_max = corner_def['x_min'], corner_def['x_max']
        y_min, y_max = corner_def['y_min'], corner_def['y_max']
        macro_row = corner_def['macro_row']
        macro_col = corner_def['macro_col']

        try:
            # Create corner rectangle in assembly coordinates
            width = x_max - x_min
            height = y_max - y_min
            corner_center_x = (x_min + x_max) / 2.0
            corner_center_y = (y_min + y_max) / 2.0

            corner_rect = Rectangle(
                name=f"CORNER_CB_{corner_id}",
                height=height,
                width=width,
                center=(corner_center_x, corner_center_y, z_center)
            )

            # Compute channel box material portion: corner_rect ∩ channel_box_rect - coolant_rect
            temp_geom = corner_rect * channel_box_rect  # intersection
            corner_cb_geom = temp_geom - coolant_rect   # difference (remove coolant)

            # Wrap and validate
            corner_cb_geom = _wrap_and_validate_geometry(corner_cb_geom, f"CORNER_CB_{corner_id}")

            if corner_cb_geom is not None and not is_degenerate(corner_cb_geom):
                macro_name = f"MACRO{macro_row}{macro_col}"
                result[corner_id] = {
                    'geometry': corner_cb_geom,
                    'macro': macro_name
                }
                print(f"  Corner {corner_id}: Channel box geometry computed for {macro_name}")
            else:
                result[corner_id] = {
                    'geometry': None,
                    'macro': f"MACRO{macro_row}{macro_col}"
                }
                print(f"  Corner {corner_id}: No channel box geometry (degenerate or absent)")

        except Exception as e:
            print(f"  Warning: Failed to compute channel box region for corner {corner_id}: {e}")
            macro_name = f"MACRO{macro_row}{macro_col}"
            result[corner_id] = {
                'geometry': None,
                'macro': macro_name
            }

    return result


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


def discretize_box(assembly_universe, assembly_model, box_discretization_config):
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
    assembly_universe : CartesianCell
        The assembly universe as returned by ``build_assembly_with_macros``.
        (3-region for uncontrolled, multi-region for controlled).
    assembly_model : CartesianAssemblyModel
        Assembly model providing dimensional information.
    box_discretization_config : BoxDiscretizationConfig
        Configuration specifying the grid splits for corner, side,
        and (optionally) cross-affected strips.

    Returns
    -------
    assembly_universe : CartesianCell
        The updated universe with a finer technological geometry.
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
            ctrl_shapes = getattr(assembly_universe, "_ctrl_cross_shapes", None)
            if ctrl_shapes is None:
                import warnings
                warnings.warn(
                    "Control cross sub-mesh requested but "
                    "_ctrl_cross_shapes not attached to "
                    "assembly_universe.  Skipping control cross "
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
        [assembly_universe],
        splitting_faces,
        shape_type=ShapeType.COMPOUND,
    )
    assembly_universe.geometry_maps[GeometryType.SECTORIZED] = \
    assembly_universe.get_geometry_map(GeometryType.SECTORIZED) // wrap_shape(partitioned_face)

    return assembly_universe


# ==============================================================================
# PHASE 3: Build Assembly with MACRO Regions
# ==============================================================================

def build_assembly_with_macros(assembly_model, calculation_step, center=None):
    """
    Build complete assembly with fuel lattice and layer-aware MACRO subdivision regions.

    This function orchestrates the assembly construction with a new layer-aware MACRO
    subdivision approach:

    1. Generate fuel cells with sectorization (if calculation_step provided)
    2. Create CartesianLattice and add fuel cells
    3. Add vanished cells and water rods to lattice
    4. Add lattice to assembly
    5. Compute lattice footprint (with asymmetric gap support)
    6. Generate MACRO rectangles around lattice
    7. Get reference rectangles for layer boundaries (coolant, channel_box)
    8. For each MACRO, generate layer-aware sub-regions (COOLANT, CHANNEL_BOX, MODERATOR)
    9. Create Region objects with both MATERIAL and MACRO properties
    10. Add all layer-aware MACRO regions to assembly

    CHANGES IN THIS VERSION:
    - Base regions (inner_box_coolant, channel_box, outer_moderator) are NOT created
    - Instead, MACRO rectangles are subdivided into layer-specific geometries
    - Each layer portion maintains the same MACRO property but different MATERIAL
    - Layer boundaries respect asymmetric gap configurations

    Parameters
    ----------
    assembly_model : CartesianAssemblyModel
        Assembly model with lattice description and dimensions
    calculation_step : CalculationStep or None
        Optional calculation step for fuel cell sectorization
    center : tuple or None
        Assembly center (x, y, z). Defaults to (ap/2, ap/2, 0)

    Returns
    -------
    assembly_universe : CartesianCell
        Complete assembly with lattice and layer-aware MACRO regions, all properties set
    """
    # Extract dimensions and translation offsets
    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_rows = len(assembly_model.lattice_description)
    n_cols = len(assembly_model.lattice_description[0])
    lattice_pitch_x = n_cols * pin_pitch
    lattice_pitch_y = n_rows * pin_pitch

    if center is None:
        center = (ap / 2.0, ap / 2.0, 0.0)

    translation_x = assembly_model.translation_offset_x if assembly_model.translation_offset_x is not None else 0.0
    translation_y = assembly_model.translation_offset_y if assembly_model.translation_offset_y is not None else 0.0

    # ======================================================================
    # STEP 1: Generate fuel cells with sectorization
    # ======================================================================
    if calculation_step:
        calculation_step.apply_radii(assembly_model)

    ordered_cells = generate_fuel_cells(
        assembly_model, calculation_step=calculation_step
    )

    # ======================================================================
    # STEP 2: Initialize assembly (WITHOUT base regions)
    # ======================================================================
    # Create assembly container directly without detailed base regions
    # Base regions will be replaced by layer-aware MACRO subdivisions
    # Use CartesianCell with width_height and center parameters to define bounds
    assembly_universe = CartesianCell(
        name=f"{assembly_model.name}_universe",
        width_height=(ap, ap),
        center=center,
        base_props={PropertyType.MATERIAL: "COOLANT",
                    PropertyType.MACRO:"BASE_CELL"},
    )

    # ======================================================================
    # STEP 3-5: Create lattice and add all cells
    # ======================================================================
    lattice = CartesianLattice(
        name=f"{assembly_model.name}_lattice",
        centre=center,
        cells=[],
    )

    # Add fuel cells to lattice
    lattice = add_cells_to_cartesian_lattice(
        lattice, ordered_cells, pin_pitch,
        translation_x=translation_x, translation_y=translation_y
    )

    # Add water rods if present
    if hasattr(assembly_model, "water_rods") and assembly_model.water_rods:
        lattice = create_and_add_water_rods_to_lattice(
            lattice, assembly_model,
            translation_x=translation_x, translation_y=translation_y,
            calculation_step=calculation_step,
        )

    # Add vanished rods if present
    if hasattr(assembly_model, "vanished_rods") and assembly_model.vanished_rods:
        lattice = add_vanished_rods_to_lattice(
            lattice, assembly_model,
            translation_x=translation_x, translation_y=translation_y,
            calculation_step=calculation_step,
        )

    # Add lattice to assembly
    assembly_universe.add(lattice)

    # ======================================================================
    # STEP 6: Compute lattice footprint (with asymmetric gap support)
    # ======================================================================
    x0 = translation_x if translation_x != 0.0 else (ap - lattice_pitch_x) / 2.0
    y0 = translation_y if translation_y != 0.0 else (ap - lattice_pitch_y) / 2.0
    x1 = x0 + lattice_pitch_x
    y1 = y0 + lattice_pitch_y

    # ======================================================================
    # STEP 7: Generate MACRO rectangles
    # ======================================================================
    macro_rects = _generate_macro_subdivision_rectangles(
        assembly_model,
        x0, y0, x1, y1,
        center
    )

    # ======================================================================
    # STEP 8: Get reference rectangles for layer boundary computation
    # ======================================================================
    lattice_rect, coolant_rect, channel_box_rect = _compute_asymmetric_coolant_channel_box_rects(
        assembly_model, center
    )

    # ======================================================================
    # STEP 9-10: Create layer-aware MACRO region subdivisions
    # ======================================================================
    print(f"\n=== Creating Layer-Aware MACRO Regions ===")
    for macro_name, rect_params in macro_rects.items():
        name, height, width, rect_center, rounded = rect_params

        # Create Rectangle geometry for this MACRO
        print(f"\nProcessing MACRO '{macro_name}' with rect center {rect_center}, height {height}, width {width}")
        macro_rect = Rectangle(
            name=name,
            height=height,
            width=width,
            center=rect_center,  # Keep center for boolean operations to work correctly
            rounded_corners=rounded
        )

        # Generate layer-specific sub-regions using set operations
        layer_geometries = _generate_macro_region_layers(
            macro_rect,
            macro_name,
            coolant_rect,
            channel_box_rect
        )

        # Create Region objects for each layer portion
        for material_type, geometry in layer_geometries.items():
            if geometry is None:  # Skip if layer doesn't exist for this MACRO
                continue

            try:
                # Create Region with layer-aware geometry
                # All portions of a MACRO (e.g., BOT_1_COOLANT and BOT_1_CHANNEL_BOX)
                # share the same MACRO property identifier (BOT_1)
                region = Region(
                    geometry,
                    properties={
                        PropertyType.MATERIAL: material_type,
                        PropertyType.MACRO: macro_name  # Keep original MACRO name
                    }
                )

                # Add region to assembly_universe at the CENTER of THIS LAYER'S GEOMETRY.
                # Boolean operations (intersection/difference) produce layer-specific geometries
                # with different centers than the original rect_center. We must use each layer's
                # actual center (computed as CDG by Region) to avoid stacking layers on top of each other.
                layer_center = get_point_coordinates(region.o)
                assembly_universe.add(region, position=layer_center)

            except Exception as e:
                print(f"  Warning: Failed to create Region for {macro_name}/{material_type}: {e}")

    print(f"\n=== Finished Creating Layer-Aware MACRO Regions ===\n")

    # ======================================================================
    # STEP 11: Generate corner channel box regions
    # ======================================================================
    # These regions bridge between lattice footprint and channel box at corners,
    # assigning them to the corner fuel pin macros (MACRO00, MACRO0{n_cols-1}, etc.)
    print(f"\n=== Creating Corner Channel Box Regions ===")
    corner_cb_regions = _generate_corner_channel_box_regions(
        x0, y0, x1, y1, ap,
        n_rows, n_cols,
        coolant_rect, channel_box_rect,
        center
    )

    for corner_id, corner_data in corner_cb_regions.items():
        geometry = corner_data['geometry']
        macro_name = corner_data['macro']

        if geometry is None:
            continue

        try:
            region = Region(
                geometry,
                properties={
                    PropertyType.MATERIAL: "CHANNEL_BOX",
                    PropertyType.MACRO: macro_name
                }
            )
            layer_center = get_point_coordinates(region.o)
            assembly_universe.add(region, position=layer_center)
            print(f"  Corner {corner_id}: Created CHANNEL_BOX region with {macro_name}")
        except Exception as e:
            print(f"  Warning: Failed to create corner CB region {corner_id}: {e}")

    print(f"=== Finished Creating Corner Channel Box Regions ===\n")

    # If control cross is present, build control cross shapes to assembly_universe for later use in box discretization
    if hasattr(assembly_model, "control_cross") and assembly_model.control_cross is not None:
        elements = _build_control_cross_elements(
            ctrl=assembly_model.control_cross,
            ap=ap,
            assembly_center=center
        )
        for element_name in elements.keys():
            geom_obj = elements[element_name]["geometry"]
            material = elements[element_name]["material"]
            macro = elements[element_name]["macro"]
            
            if "offset" in elements[element_name].keys():
                offset = elements[element_name]["offset"] 
            else:
                offset = (0.0, 0.0, 0.0)

            try: 
                region = Region(
                    geom_obj=geom_obj,
                    properties={
                        PropertyType.MATERIAL: material,
                        PropertyType.MACRO: macro
                    }
                )
                region_center = get_point_coordinates(region.o)
                offset_center = (region_center[0]+offset[0],region_center[1]+offset[1], region_center[2]+offset[2])
                assembly_universe.add(region, position=offset_center)
            except Exception as e:
                print(f"  Warning failed to create Region for control cross element : {element_name}")
           

    return assembly_universe



def build_full_assembly_geometry(assembly_model, calculation_step,
                                 output_path, output_file_name):
    """
    High-level function that builds a complete assembly geometry with fuel lattice,
    MACRO subdivision regions (if IC method), and optional MOC discretization.
    Exports result to a TDT file.

    This function orchestrates the full pipeline:

    1. Build assembly with fuel lattice and MACRO regions (all cells, water rods, vanished rods)
    2. Apply symmetry
    3. Optionally apply MOC box discretization if enabled
    4. Export the TDT file

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
    assembly_universe : CartesianCell
        The complete assembly geometry with all properties set and exported.
    """
    from ..DDModel.DragonModel import CartesianAssemblyModel  # type check only

    ap = assembly_model.assembly_pitch
    pin_pitch = assembly_model.pin_geometry_dict["pin_pitch"]
    n_cols = len(assembly_model.lattice_description[0])

    print(f"Building assembly geometry for {assembly_model.name} with {calculation_step.spatial_method} method")

    # ======================================================================
    # STEP 1: Build complete assembly with MACRO regions
    # ======================================================================
    # This function handles:
    # - Fuel cell generation and sectorization
    # - Assembly box creation (coolant/channel_box/moderator)
    # - Lattice creation and population with all cell types
    # - MACRO region generation with properties
    assembly_universe = build_assembly_with_macros(
        assembly_model,
        calculation_step=calculation_step
    )

    # ======================================================================
    # STEP 2: Apply symmetry
    # ======================================================================
    assembly_universe.update_hierarchical_structure(True)
    # check for symmetries
    is_anti_diag_symmetric = assembly_model.check_anti_diagonal_symmetry()
    is_main_diag_symmetric = assembly_model.check_main_diagonal_symmetry()
    is_quarter_symmetric = assembly_model.check_quarter_symmetry()
    is_half_symmetric = assembly_model.check_half_symmetry()
    if is_quarter_symmetric and not (is_anti_diag_symmetric or is_main_diag_symmetric):
        print("Assembly is quarter symmetric; applying quarter symmetry")
        assembly_universe.apply_symmetry(SymmetryType.QUARTER)
        symmetry = SymmetryType.QUARTER
    elif is_quarter_symmetric and is_anti_diag_symmetric and is_main_diag_symmetric:
        print("Assembly is quarter symmetric with diagonal symmetry; applying eigth symmetry")
        assembly_universe.apply_symmetry(SymmetryType.EIGHTH)
        symmetry = SymmetryType.EIGHTH
    elif is_half_symmetric:
        print("Assembly is half symmetric; applying half symmetry")
        assembly_universe.apply_symmetry(SymmetryType.HALF)
        symmetry = SymmetryType.HALF
    elif is_main_diag_symmetric:
        print("Assembly is symmetric across a diagonal; applying diagonal symmetry")
        assembly_universe.apply_symmetry(SymmetryType.DIAG)
        symmetry = SymmetryType.DIAG
    elif is_anti_diag_symmetric:
        print("Assembly is symmetric across anti-diagonal; applying anti-diagonal symmetry")
        #assembly_universe.rotate(90)  # Rotate 90 degrees to align anti-diagonal with main diagonal
        assembly_universe.apply_symmetry(SymmetryType.FULL)
        symmetry = SymmetryType.FULL
    else:
        print("Assembly has no symmetries; using full geometry")
        assembly_universe.apply_symmetry(SymmetryType.FULL)
        symmetry = SymmetryType.FULL

    


    # ======================================================================
    # STEP 3: Optionally apply MOC box discretization
    # ======================================================================
    # If MOC method with box discretization enabled, subdivide assembly box further
    # for fine MOC tracking grid
    if (calculation_step.box_discretization is not None
            and calculation_step.box_discretization.enabled):
        print(f"Box discretization enabled for MOC tracking; subdividing assembly box")
        assembly_universe = discretize_box(
            assembly_universe, assembly_model,
            calculation_step.box_discretization,
        )
    assembly_universe.show(PropertyType.MATERIAL, GeometryType.TECHNOLOGICAL)

    # ======================================================================
    # STEP 4: Export to TDT file
    # ======================================================================
    export_glow_geom(
        output_path,
        output_file_name,
        assembly_universe,
        symmetry_type=symmetry,
        tracking_option=calculation_step.tracking,
        export_macro=calculation_step.export_macros,
    )

    return assembly_universe

