## Some helpful functiont to deal with GLOW geometry building
from glow.geometry_layouts.cells import RectCell
from glow.geometry_layouts.geometries import Rectangle
from glow.support.types import GeometryType, PropertyType, SymmetryType
from glow.geometry_layouts.lattices import Lattice
from glow.main import TdtSetup, analyse_and_generate_tdt
from glow.interface.geom_interface import *
from glow.support.types import *


def computeSantamarinaradii(fuel_radius, gap_radius, clad_radius, gadolinium=False):
    """
    Helper to define fuel region radii for fuel pins
    A. Santamarina recommendations:
    - UOX pins: 50%, 80%, 95% and 100% volume fractions
    - Gd2O3 pins: 20%, 40%, 60%, 80%, 95% and 100% volume fractions
    """
    if not gadolinium:
        pin_radii = [
            (0.5**0.5) * fuel_radius,
            (0.8**0.5) * fuel_radius,
            (0.95**0.5) * fuel_radius,
            fuel_radius,
            gap_radius,
            clad_radius
        ]
    else:
        pin_radii = [
            (0.2**0.5) * fuel_radius,
            (0.4**0.5) * fuel_radius,
            (0.6**0.5) * fuel_radius,
            (0.8**0.5) * fuel_radius,
            (0.95**0.5) * fuel_radius,
            fuel_radius,
            gap_radius,
            clad_radius
        ]
    return pin_radii


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
    #Gd_cells = ["ROD5G", "ROD6H", "ROD6K", "ROD7G", "ROD7H"]
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
                center=(pitch / 2, pitch / 2, 0.0),
                rounded_corners=rounded_corners
            )
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
                    cell,
                    (
                        (cell_idx + 0.5) * cell_pitch + translation,
                        (row_idx + 0.5) * cell_pitch + translation,
                        0.0
                    )
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
    lattice.apply_symmetry(SymmetryType.FULL)
    if export_macro:
        propertis_to_export = [PropertyType.MATERIAL, PropertyType.MACRO]
    else:
        propertis_to_export = [PropertyType.MATERIAL]

    if tracking_option == "TISO":
        lattice.type_geo = LatticeGeometryType.ISOTROPIC
        analyse_and_generate_tdt(
            [lattice], f"{output_path}/{output_file_name}", TdtSetup(GeometryType.SECTORIZED,
                                             property_types=propertis_to_export,
                                             type_geo=LatticeGeometryType.ISOTROPIC,
                                             symmetry_type=BoundaryType.AXIAL_SYMMETRY))
    elif tracking_option == "TSPC":
        lattice.type_geo = LatticeGeometryType.RECTANGLE_SYM
        analyse_and_generate_tdt(
            [lattice], f"{output_path}/{output_file_name}", TdtSetup(GeometryType.SECTORIZED,
                                             property_types=propertis_to_export,
                                             type_geo=LatticeGeometryType.RECTANGLE_SYM,
                                             symmetry_type=BoundaryType.AXIAL_SYMMETRY))
