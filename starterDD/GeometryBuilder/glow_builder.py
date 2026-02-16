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
    # parent width/height and lower-left corner
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


def generate_IC_cells(lattice_desc, Gd_cells, pitch, C_to_mat, fuel_rad, gap_rad, clad_rad, corner_radius=0, windmill=False):
    """
    Generate RectCell objects for each individual subgeometry in the lattice
    
    Parameters:
    -----------
    lattice_desc : list of list of str
        2D list describing the lattice layout with cell IDs
    Gd_cells : list of str
        List of cell IDs that correspond to gadolinium fuel cells
    pitch : float
        Pitch of each cell in the lattice
    C_to_mat : dict
        Mapping from cell IDs to material names
    fuel_rad : float
        Radius of the fuel region in fuel pins
    gap_rad : float
        Radius of the gap region in fuel pins
    clad_rad : float
        Radius of the cladding region in fuel pins
    corner_radius : float, optional
        Radius of curvature for the outer corner of corner fuel cells.
        Only the 4 corner cells (positions (0,0), (0,n-1), (n-1,0), (n-1,n-1)) 
        will have a single rounded corner to match the channel box inner boundary.
    windmill : bool, optional
        Whether to apply windmill sectorization to fuel pins.
    """
    lattice_components = []
    n_rows = len(lattice_desc)
    n_cols = len(lattice_desc[0]) if lattice_desc else 0
    
    for row_idx in range(n_rows):
        row = lattice_desc[row_idx]
        row_of_cells = []
        for cell_idx in range(len(row)):
            cell_id = row[cell_idx]
            
            # Check if this is a corner cell that needs a rounded corner
            # to match the channel box inner boundary
            rounded_corners = None
            if corner_radius > 0:
                # Corner 0 = bottom-left of cell (row=0, col=0 -> assembly corner)
                # Corner 1 = bottom-right of cell (row=0, col=n_cols-1 -> assembly corner)
                # Corner 2 = top-right of cell (row=n_rows-1, col=n_cols-1 -> assembly corner)
                # Corner 3 = top-left of cell (row=n_rows-1, col=0 -> assembly corner)
                if row_idx == 0 and cell_idx == 0:
                    rounded_corners = [(0, corner_radius)]
                elif row_idx == 0 and cell_idx == n_cols - 1:
                    rounded_corners = [(1, corner_radius)]
                elif row_idx == n_rows - 1 and cell_idx == n_cols - 1:
                    rounded_corners = [(2, corner_radius)]
                elif row_idx == n_rows - 1 and cell_idx == 0:
                    rounded_corners = [(3, corner_radius)]
            
            tmp_cell = RectCell(
                name=cell_id,
                height_x_width=(pitch, pitch),
                center=(0.0, 0.0, 0.0),
                rounded_corners=rounded_corners
            )
            if C_to_mat is None:
                mat_name = cell_id  # If no mapping provided, use cell ID as material name (assuming they match)
            else:
                mat_name = C_to_mat[cell_id]
            
            if cell_id in Gd_cells:
                radii = computeSantamarinaradii(fuel_rad, gap_rad, clad_rad, gadolinium=True)
            elif "ROD" in cell_id and "W" not in cell_id:
                radii = computeSantamarinaradii(fuel_rad, gap_rad, clad_rad, gadolinium=False)
            else:  # Water rod placeholder
                radii = []
                mat_name = "MODERATOR"
            
            for radius in radii:
                tmp_cell.add_circle(radius)

            if mat_name == "MODERATOR":
                tmp_cell.set_properties({
                    PropertyType.MATERIAL: ["MODERATOR"],
                    PropertyType.MACRO: [f"MACRO{row_idx}{cell_idx}"]
                })
            else:
                if cell_id in Gd_cells:
                    list_of_cell_mats = [mat_name] * 6
                    if windmill:
                        tmp_cell.sectorize([1]*8 + [8], [0]*8 + [22.5], windmill=True)
                else:
                    list_of_cell_mats = [mat_name] * 4
                    if windmill:
                        tmp_cell.sectorize([1]*6 + [8], [0]*6 + [22.5], windmill=True)
                list_of_cell_mats.extend(["GAP", "CLAD", "COOLANT"])
                tmp_cell.set_properties({
                    PropertyType.MATERIAL: list_of_cell_mats,
                    PropertyType.MACRO: [f"MACRO{row_idx}{cell_idx}"] * len(list_of_cell_mats)
                })
            row_of_cells.append(tmp_cell)
        lattice_components.append(row_of_cells)
    return lattice_components


def generate_fuel_cells(assemblyModel):
    """
    Generate RectCell objects for each individual subgeometry in the lattice
    
    Parameters:
    -----------
    assemblyModel : CartesianAssemblyModel
        The assembly model containing the lattice description and rod ID to material mapping
        additionally contains the pin geometry parameters needed to define the fuel cells and the material mixture unique names to assign to the cell materials.
    
    Returns:
    -----------
    lattice_components : list of list of RectCell
        2D list of RectCell objects representing the lattice layout, with properties set according to assemblyModel information.
        To be used to build lattice geometry with glow and export to TDT file.
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
                print(list_of_cell_mats)
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
    Add fuel cells to the lattice, skipping water rod placeholders
    Parameters:
    -----------
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


def create_and_add_water_rods_to_lattice(lattice, assembly_model, translation=0.0, windmill=False):
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
        Whether to apply windmill sectorization to the water rod coolant region
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
            if windmill:
                tmp_cell.sectorize([1, 1, 8], [0, 0, 0], windmill=True)
        elif assembly_model.water_rod_type == "square":
            raise NotImplementedError(
                "Square water rod geometry is not yet implemented in the glow builder."
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
    Export the geometry of the lattice to a TDT file for GLOW simulation
    Parameters:
    -----------
    output_path : str
        Path to save the exported TDT file
    output_file_name : str
        Name of the exported TDT file
    lattice : Lattice
        The lattice whose geometry is to be exported
    tracking_option : str
        Tracking option, either "TISO" or "TSPC"

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

