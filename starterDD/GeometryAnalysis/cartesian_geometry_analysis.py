import os
import numpy as np
import yaml
from shapely.geometry import box, Point, MultiLineString
from shapely.ops import unary_union
import math

_GL_NODES   = np.array([-0.90617984593866399, -0.53846931010568309,  0.0,
                          0.53846931010568309,  0.90617984593866399])
_GL_WEIGHTS = np.array([ 0.23692688505618908,  0.47862867049936647,  0.56888888888888889,
                          0.47862867049936647,  0.23692688505618908])

def _extract_box_dimensions(dragon_assembly_model):
    """
    Extrait les dimensions du boîtier. Si elles sont commentées dans le YAML,
    on les recalcule dynamiquement à partir des surfaces (assembly_box_description).
    """
    L_ext = dragon_assembly_model.assembly_pitch
    gap_wide = getattr(dragon_assembly_model, 'gap_wide', None)
    cbt = getattr(dragon_assembly_model, 'channel_box_thickness', None)
    
    # 1. Si les paramètres sont définis explicitement
    if L_ext is not None and gap_wide is not None and cbt is not None:
        W_start = gap_wide + cbt
        L_int = L_ext - 2 * W_start
        return W_start, L_ext, L_int
        
    # 2. Sinon, on lit les surfaces brutes !
    ass_geo = getattr(dragon_assembly_model, 'raw_yaml_data', {}).get('ASSEMBLY_GEOMETRY', {})
    surfaces = ass_geo.get('assembly_box_description', {}).get('surfaces', [])
    
    pb_side = None
    icb_side = None
    for surf in surfaces:
        if surf['name'] == 'problem_boundaries':
            pb_side = surf['parameters']['side_length']
        elif surf['name'] == 'inner_channel_box':
            icb_side = surf['parameters']['side_length']
            
    if pb_side is not None and icb_side is not None:
        L_ext = float(pb_side)
        L_int = float(icb_side)
        W_start = (L_ext - L_int) / 2.0
        return W_start, L_ext, L_int
        
    raise ValueError(f"Impossible de déduire les dimensions du boîtier pour {dragon_assembly_model.name}. Vérifiez le YAML.")

def build_assembly_geometry(dragon_assembly_model, r_wr_override=None):
    data_ref = {}
    data_ref["ASSEMBLY_GEOMETRY"] = {}
    data_ref["PIN_GEOMETRY"] = {}
    data_ref["WATER_ROD_GEOMETRY"] = {}
    data_ref["WATER_ROD_GEOMETRY"]["centers"] = []

    pin_geo = dragon_assembly_model.pin_geometry_dict
    lattice = dragon_assembly_model.lattice_description
    exclusions = set(dragon_assembly_model.non_fuel_rod_ids)

    d = pin_geo['pin_pitch']
    data_ref["PIN_GEOMETRY"]["pin_pitch"] = d
    r_clad = pin_geo['clad_radius']
    data_ref["PIN_GEOMETRY"]["clad_radius"] = r_clad
    
    data_ref["ASSEMBLY_GEOMETRY"]["lattice_description"] = lattice
    data_ref["ASSEMBLY_GEOMETRY"]["corner_inner_radius_of_curvature"] = dragon_assembly_model.corner_inner_radius_of_curvature

    if r_wr_override is not None:
        r_wr = r_wr_override
    elif dragon_assembly_model.water_rod_type == "circular":
        r_wr = dragon_assembly_model.water_rod_outer_radius
    elif dragon_assembly_model.water_rod_type == "conical":
        r_wr = getattr(dragon_assembly_model, "water_rod_outer_radius_start", 0.0)
    else:
        raise ValueError(f"Analyzer does not support water rods with type {dragon_assembly_model.water_rod_type}")
    
    data_ref["WATER_ROD_GEOMETRY"]["outer_radius"] = r_wr
    
    grid_thickness = pin_geo.get('grid_thickness', 0.0)
    r_solid_clad = r_clad + grid_thickness
    r_solid_wr   = r_wr   + grid_thickness
    
    # --- APPEL DE NOTRE NOUVELLE FONCTION INTELLIGENTE ---
    W_start, L_ext, L_int = _extract_box_dimensions(dragon_assembly_model)
    R_c = dragon_assembly_model.corner_inner_radius_of_curvature
    W_end = L_ext - W_start

    # Sauvegarde dans le dictionnaire pour plus tard
    data_ref["ASSEMBLY_GEOMETRY"]["assembly_pitch"] = L_ext
    data_ref["ASSEMBLY_GEOMETRY"]["W_start"] = W_start
    data_ref["ASSEMBLY_GEOMETRY"]["L_int"] = L_int

    inner_box = box(W_start + R_c, W_start + R_c, W_end - R_c, W_end - R_c).buffer(R_c, quad_segs=64)

    formes_solides = []
    lignes_chauffantes = []
    lignes_water = []

    for water_rod in dragon_assembly_model.water_rods:
        center = water_rod.center
        data_ref["WATER_ROD_GEOMETRY"]["centers"].append(center)
        wr_circle = Point(center[0], center[1]).buffer(r_solid_wr, quad_segs=64)
        formes_solides.append(wr_circle)
        lignes_water.append(wr_circle.exterior)
    
    l_gap_int = (L_int - ((len(lattice[0]) - 1) * d) - 2 * r_clad) / 2.0
    start_x = W_start + l_gap_int + r_clad
    start_y = W_start + l_gap_int + r_clad

    for row_idx, row in enumerate(lattice):
        for col_idx, item in enumerate(row):
            cx = start_x + col_idx * d
            cy = start_y + row_idx * d

            if item == 'VROD':
                if grid_thickness > 0.0:
                    cercle_ext = Point(cx, cy).buffer(r_solid_clad, quad_segs=64)
                    cercle_int = Point(cx, cy).buffer(r_clad, quad_segs=64)
                    anneau = cercle_ext.difference(cercle_int)
                    formes_solides.append(anneau)
                    lignes_water.append(cercle_ext.exterior)
                    lignes_water.append(cercle_int.exterior)
                continue

            if item in exclusions:
                continue

            rod_circle = Point(cx, cy).buffer(r_solid_clad, quad_segs=64)
            formes_solides.append(rod_circle)
            lignes_chauffantes.append(Point(cx, cy).buffer(r_clad, quad_segs=64).exterior)

    solide_total = unary_union(formes_solides)
    multi_lignes_chauffantes = MultiLineString(lignes_chauffantes)
    multi_lignes_water = MultiLineString(lignes_water)

    return inner_box, solide_total, multi_lignes_chauffantes, multi_lignes_water, data_ref

def analyse_mesh(x1, y1, x2, y2, inner_box, solide_total, lignes_chauffantes, lignes_water):
    """
    Calculate geometric properties for a given mesh cell defined by (x1, y1, x2, y2).
    
    Parameters:
    - (x1, y1): Bottom-left corner of the cell
    - (x2, y2): Top-right corner of the cell
    - inner_box: Polygon of the internal region (excluding walls)
    - solide_total: Polygon of the solid region (fuel rods + water rods)
    - lignes_chauffantes: MultiLineString of heating perimeters (e.g. rod boundaries)
    - lignes_water: MultiLineString of water rod perimeters
    
    Returns:
    - s_tot_valide: Valid total area of the cell within the inner box
    - s_m: Area of the fluid (non-solid) part of the cell
    - p_m: Wetted perimeter of the fluid part (interface with solids)
    - p_h: Perimeter of the cell that is in contact with heating surfaces (rod boundaries)
    - p_box_in_cv: Perimeter of the cell that is in contact with the inner box boundary
    - p_wr: Perimeter of the cell that is in contact with water rod boundaries
    """
    cv = box(x1, y1, x2, y2)
    valid_cv = cv.intersection(inner_box)
    
    s_tot_valide = valid_cv.area 
    if s_tot_valide <= 1e-9:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    s_solid_rect = solide_total.intersection(valid_cv).area
    p_box_in_cv = inner_box.exterior.intersection(valid_cv).length
    p_m_rods = solide_total.boundary.intersection(valid_cv).length
    p_h = lignes_chauffantes.intersection(valid_cv).length
    p_wr = lignes_water.intersection(valid_cv).length

    s_m = s_tot_valide - s_solid_rect
    p_m = p_box_in_cv + p_m_rods
    
    return s_tot_valide, s_m, p_m, p_h, p_box_in_cv, p_wr

def analyse_3d_volume(tranches_axiales, x1, y1, x2, y2, z1, z2):
    """
    Calculate porosity, fluid area, hydraulic diameter, heating perimeter, inner box perimeter, and water rod perimeter for a 3D volume 
    defined by (x1, y1, z1) to (x2, y2, z2) across multiple axial slices.
    
    Parameters:
    - tranches_axiales: List of axial slices, 
    each with its own geometry (inner_box, solide, lignes_chauffantes, lignes_water) and z_start/z_end
    - (x1, y1): Bottom-left corner of the vertical prism
    - (x2, y2): Top-right corner of the vertical prism
    - z1: Starting axial coordinate of the prism
    - z2: Ending axial coordinate of the prism      
    
    Returns:
    - poro_3d: Porosity of the volume (fluid volume / total volume)
    - a_cool: Total area of fluid in the volume (sum of fluid area across slices * slice thickness)
    - dh_3d: Hydraulic diameter of the fluid region (4 * fluid volume / wetted perimeter)
    - ph_moyen: Average heating perimeter in contact with the fluid across the axial height
    - pbox_moyen: Average inner box perimeter in contact with the fluid across the axial height
    - pwr_moyen: Average water rod perimeter in contact with the fluid across the axial height
    """
    v_total = 0.0
    v_fluide = 0.0
    s_mouillee = 0.0
    s_chauffante = 0.0
    s_box = 0.0
    s_wr = 0.0
    
    hauteur_totale = z2 - z1
    if hauteur_totale <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    for tranche in tranches_axiales:
        z_min_overlap = max(z1, tranche['z_start'])
        z_max_overlap = min(z2, tranche['z_end'])
        dz = z_max_overlap - z_min_overlap

        if dz <= 0:
            continue

        if tranche.get('is_conical', False):
            # Gauss-Legendre quadrature: integrate exactly over the z-overlap interval.
            # r(z) is linear in z, so A(z) = π r(z)² is quadratic — GL-5 is exact.
            z_s    = tranche['z_start']
            z_e    = tranche['z_end']
            r_s    = tranche['r_cone_start']
            r_e    = tranche['r_cone_end']
            data   = tranche['dragon_assembly_model']
            z_mid_ov = (z_min_overlap + z_max_overlap) / 2.0
            half_dz  = dz / 2.0
            for xi, wi in zip(_GL_NODES, _GL_WEIGHTS):
                z_i  = z_mid_ov + half_dz * xi
                frac = (z_i - z_s) / (z_e - z_s)
                r_i  = r_s + (r_e - r_s) * frac
                ib_i, sol_i, lch_i, lwr_i, _ = build_assembly_geometry(data, r_wr_override=r_i)
                st, sm, pm, ph, pb, pw = analyse_mesh(x1, y1, x2, y2, ib_i, sol_i, lch_i, lwr_i)
                v_total      += wi * st * half_dz
                v_fluide     += wi * sm * half_dz
                s_mouillee   += wi * pm * half_dz
                s_chauffante += wi * ph * half_dz
                s_box        += wi * pb * half_dz
                s_wr         += wi * pw * half_dz
        else:
            s_tot, s_m, p_m, p_h, p_box, p_wr = analyse_mesh(
                x1, y1, x2, y2,
                tranche['inner_box'], tranche['solide'], tranche['lignes_chauffantes'], tranche['lignes_water']
            )
            v_total      += s_tot * dz
            v_fluide     += s_m   * dz
            s_mouillee   += p_m   * dz
            s_chauffante += p_h   * dz
            s_box        += p_box * dz
            s_wr         += p_wr  * dz

    if v_total <= 1e-9:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    poro_3d = v_fluide / v_total
    a_cool = v_fluide / (z2-z1)
    dh_3d = (4 * v_fluide) / s_mouillee if s_mouillee > 0 else 0.0
    ph_moyen = s_chauffante / hauteur_totale
    pbox_moyen = s_box / hauteur_totale
    pwr_moyen = s_wr / hauteur_totale
    
    return round(poro_3d, 5), round(a_cool, 5), round(dh_3d, 5), round(ph_moyen, 5), round(pbox_moyen, 5), round(pwr_moyen, 5)

class CartesianGeometricAnalyser:
    """
    Class responsible for analyzing the geometry of a nuclear reactor assembly 
    based on YAML input files. 
    It provides methods to compute porosity, hydraulic diameter, 
    and heating perimeter for specified volumes or profiles within the assembly.
    """
    def __init__(self, core_model, core_i=1, core_j=1):
        """
        Initialize the geometric analyzer by reading the core model,
        extracting the relevant axial slices for the specified assembly position.
        
        Parameters:
        - core_model: Object containing the overall geometry
        - core_i: Row index of the assembly in the core layout (1 is bottom row)
        - core_j: Column index of the assembly in the core layout (1 is leftmost column)
        """
        self.slices_data = []
        self.data_ref = None
        self.core_model = core_model 
        
        # 1. On lit la carte du cœur et on vérifie les indices
        core_map_2D = core_model.core_2D_layout
        max_i = len(core_map_2D)
        max_j = len(core_map_2D[0])
        
        if not (1 <= core_i <= max_i and 1 <= core_j <= max_j):
            raise ValueError(f"Erreur : Assemblage ({core_i},{core_j}) hors limites. Le cœur fait {max_i}x{max_j}.")
            
        assembly_id = core_map_2D[core_i - 1][core_j - 1]
        self.core_assembly_name = assembly_id 
        print(f"Info : Chargement de l'assemblage '{assembly_id}' à la position cœur ({core_i},{core_j})")
        
        # 3. On récupère l'empilement axial pour cet assemblage précis
        fuel_assembly_bundle = core_model.assemblies[(core_i-1, core_j-1, assembly_id)]
        axial_regions = fuel_assembly_bundle.slices_2D
        z_bounds = fuel_assembly_bundle.z_bounds
        
        # 4. On boucle sur les tranches (DOM, VAN, PLENUM_BOT, PLENUM_TOP, etc.)
        for region_idx in range(len(axial_regions)):
            z_s, z_e = z_bounds[region_idx], z_bounds[region_idx+1]
            slice_id = axial_regions[region_idx]           
                    
            dragon_assembly_model = core_model.assembly_models[((core_i-1, core_j-1, assembly_id), slice_id)]
            box_geom, solide, lignes_chauffantes, lignes_water, data_slice = build_assembly_geometry(dragon_assembly_model)
            cylinders_de_cette_tranche = self._get_cylinders(dragon_assembly_model)

            is_conical = getattr(dragon_assembly_model, "water_rod_type", "") == "conical"
            if is_conical:
                r_ws = float(dragon_assembly_model.water_rod_outer_radius_start)
                r_we = float(getattr(dragon_assembly_model, "water_rod_outer_radius_end",
                                    dragon_assembly_model.water_rod_outer_radius_start))
                r_mid = (r_ws + r_we) / 2.0
                box_geom, solide, lignes_chauffantes, lignes_water, data_slice = build_assembly_geometry(
                    dragon_assembly_model, r_wr_override=r_mid)
                cylinders_de_cette_tranche = self._get_cylinders(dragon_assembly_model, r_wr_override=r_mid)
            else:
                box_geom, solide, lignes_chauffantes, lignes_water, data_slice = build_assembly_geometry(
                    dragon_assembly_model)
                cylinders_de_cette_tranche = self._get_cylinders(dragon_assembly_model)
            n_wr = len(dragon_assembly_model.water_rods) if hasattr(dragon_assembly_model, "water_rods") else 1
            k_wall = float(getattr(dragon_assembly_model, "wall_conductivity", 18.0))

            tranche_dict = {
                'z_start': float(z_s),
                'z_end': float(z_e),
                'inner_box': box_geom,
                'solide': solide,
                'lignes_chauffantes': lignes_chauffantes,
                'lignes_water': lignes_water,
                'cylinders': cylinders_de_cette_tranche,
                'axial_region': slice_id,
                'is_conical': is_conical,
                'n_wr': n_wr,
                'k_wall': k_wall,
                'dragon_assembly_model': dragon_assembly_model,
            }

            if is_conical:
                # Conical region: store a SINGLE tranche entry.
                # Volume integrals are computed analytically via GL-5 quadrature in
                # analyse_3d_volume — no discretisation into sub-steps needed.
                tranche_dict.update({
                    'r_cone_start': float(getattr(dragon_assembly_model, "water_rod_outer_radius_start", 0.0)),
                    'r_cone_end': float(getattr(dragon_assembly_model, "water_rod_outer_radius_end", 0.0)),
                    'r_cone_inner_start': float(getattr(dragon_assembly_model, "water_rod_inner_radius_start", 0.0)),
                    'r_cone_inner_end': float(getattr(dragon_assembly_model, "water_rod_inner_radius_end", 0.0)),
                    'r_wr_inner': 0.0 
                })
            else:
                tranche_dict.update({
                    'r_wr_inner': float(getattr(dragon_assembly_model, "water_rod_inner_radius", 0.0)),
                    'r_wr_outer': float(getattr(dragon_assembly_model, "water_rod_outer_radius", 0.0))
                })

            # On utilise les z_s et z_e lus dans le CORE
            self.slices_data.append(tranche_dict)
            
            if self.data_ref is None:
                self.data_ref = data_slice

        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_rows = len(lattice)
        n_cols = len(lattice[0])
        self.U_nominal_x = np.zeros((n_rows + 1, n_cols + 1))
        self.U_nominal_y = np.zeros((n_rows + 1, n_cols + 1))

        # Initialisation du cache
        self._dernier_args = None
        self._dernier_resultats = None 

    def run_THM_analysis(self, nz, include_water_rods=False):
        """
        Wrapper method to perform geometry analysis and return the necessary data to run a pyTHM or THM: case
        
        nz : integer : number of axial mesh points
        include_water_rods : boolean : include analysis of water rods present in the geometry.
        """

        geometric_data = {}
        geometric_data["active_flow_data"] = {}
        geometric_data["fuel_data"] = {} 
        geometric_data["water_rod_data"] = {}
        # --- 1. Extract geometric information ---
        z_min, maxh = self.get_z_global_bounds() # in cm
        dz = (maxh - z_min) / nz
        _, pitch_cm = self.get_x_global_bounds()
        pitch_m = pitch_cm * 1E-2
        # Analyse the axial profile for in control volume mode ('cv'), recover active flow parameters
        geom_profiles = self.execute_profile_z(['cv', [0, 0, pitch_cm, pitch_cm]], dz, dz, z_min, maxh)
        
        # Recover porosities, coolant flow cross sectional areas, hydraulic diameters, 
        porosities_profile = geom_profiles[1]
        acool_profile = [a * 1e-4 for a in geom_profiles[2]] # Conversion cm² -> m²
        dhs_profile = [dh * 1e-2 for dh in geom_profiles[3]]    # Conversion cm -> m
        phs_profile = [pch * 1e-2 for pch in geom_profiles[4]]  # Conversion cm -> m
        kexp_profile = geom_profiles[5]
        kcon_profile = geom_profiles[6]
        rsin_profile = geom_profiles[7]

        geometric_data["active_flow_data"]["number_of_axial_meshes"] = nz
        geometric_data["active_flow_data"]["porosities"] = porosities_profile
        geometric_data["active_flow_data"]["coolant_cross_sectional_areas"] = acool_profile
        geometric_data["active_flow_data"]["hydraulic_diamters"] = dhs_profile
        geometric_data["active_flow_data"]["heated_perimeters"] = phs_profile
        geometric_data["active_flow_data"]["k_expansion"] = kexp_profile
        geometric_data["active_flow_data"]["k_contraction"] = kcon_profile
        geometric_data["active_flow_data"]["reference_coolant_cross_sectional_area"] = min(acool_profile)
        geometric_data["active_flow_data"]["singular_contraction_ratios"] = rsin_profile
        geometric_data["active_flow_data"]["pitch"] = pitch_m

        if include_water_rods:
            geom_profiles_wr = self.execute_profile_z(
                ('wr_tube',),
                dz, dz, z_min, maxh
            )

            acools_wr = [a * 1e-4 for a in geom_profiles_wr[2]] # Conversion cm² -> m²
            porosities_wr = geom_profiles_wr[1]
            dhs_wr = [dh * 1e-2 for dh in geom_profiles_wr[3]]    # Conversion cm -> m
            kexp_wr = geom_profiles_wr[5]

            
            p_wr = []
            rwall_wr   = []
            curr_z = z_min
            while curr_z + dz <= maxh + 1e-10:
                z1, z2 = curr_z, curr_z + dz
                p_wr.append(self.get_pch_wr_outer(z1, z2) * 1e-2)
                rwall_wr.append(self.get_rwall_wr_tube(z1, z2))
                curr_z += dz
            wr_holes, Idelchik_exit, Idelchik_enter = self.get_wr_hole_data()

            hole_z = [h['z'] * 1e-2 for h in wr_holes]
            hole_A = [math.pi * (h['D_hole'] * 0.5e-2) ** 2 for h in wr_holes]

            geometric_data["water_rod_data"]["moderator_cross_sectional_areas"] = acools_wr
            geometric_data["water_rod_data"]["porosities"] = porosities_wr
            geometric_data["water_rod_data"]["hydraulic_diamters"] = dhs_wr
            geometric_data["water_rod_data"]["k_expansion"] = kexp_wr
            geometric_data["water_rod_data"]["permieters"] = p_wr
            geometric_data["water_rod_data"]["thermal_resistances"] = rwall_wr
            geometric_data["water_rod_data"]["Idelchik_enter"] = Idelchik_enter
            geometric_data["water_rod_data"]["Idelchik_exit"] = Idelchik_exit
            geometric_data["water_rod_data"]["hole_A"] = hole_A
            geometric_data["water_rod_data"]["hole_Z"] = hole_z 


        # --- Recover fuel pin related parameters from the DRAGON assembly model
        # Convert from cm to m for downstream use
        pin_geom = self.slices_data[0]['dragon_assembly_model'].pin_geometry_dict
        fuel_radius = pin_geom['fuel_radius'] * 1E-2
        gap_radius = pin_geom['gap_radius'] * 1E-2
        clad_radius = pin_geom['clad_radius'] * 1E-2
        pin_pitch = pin_geom['pin_pitch'] * 1E-2
        fuel_rod_length = (maxh - z_min) * 1E-2
        geometric_data["fuel_data"]["fuel_radius"] = fuel_radius
        geometric_data["fuel_data"]["gap_radius"] = gap_radius
        geometric_data["fuel_data"]["clad_radius"] = clad_radius
        geometric_data["fuel_data"]["pin_pitch"] = pin_pitch
        geometric_data["fuel_data"]["max_rod_length"] = fuel_rod_length

        return geometric_data


    def _get_box_geom(self):
        """Helper mis à jour pour lire les valeurs extraites via les surfaces"""
        L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
        W_start = self.data_ref['ASSEMBLY_GEOMETRY']['W_start']
        L_int = self.data_ref['ASSEMBLY_GEOMETRY']['L_int']
        W_end = L_ext - W_start
        return W_start, W_end, L_ext, L_int

    def _obtenir_bornes_rod(self, i, j):
        """
        Calculate the bounding box of a fuel rod based on its lattice position (i, j).
        
        Parameters:
        - i: Row index of the fuel rod in the lattice
        - j: Column index of the fuel rod in the lattice
        
        Returns:
        - Tuple of (x_min, y_min, x_max, y_max) representing the bounding box of the fuel rod
        """
        W_start, W_end, L_ext, L_int = self._get_box_geom()
        d = self.data_ref['PIN_GEOMETRY']['pin_pitch']
        r_clad = self.data_ref['PIN_GEOMETRY']['clad_radius']
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        l_gap_int = (L_int - ((n_cols - 1) * d) - 2 * r_clad) / 2.0
        
        # Calcul du centre du crayon (j=col, i=row)
        cx = W_start + l_gap_int + r_clad + (j - 1) * d
        cy = W_start + l_gap_int + r_clad + (i - 1) * d
        return cx - d/2.0, cy - d/2.0, cx + d/2.0, cy + d/2.0

    def _obtenir_bornes_water(self, i, j):
        """
        Calculate the bounding box of a channel based on its lattice position (i, j).
        
        Parameters:
        - i: Row index of the channel in the lattice
        - j: Column index of the channel in the lattice
        
        Returns:
        - Tuple of (x_min, y_min, x_max, y_max) representing the bounding box of the channel
        """
        W_start, W_end, L_ext, L_int = self._get_box_geom()
        d = self.data_ref['PIN_GEOMETRY']['pin_pitch']
        r_clad = self.data_ref['PIN_GEOMETRY']['clad_radius']
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        l_gap_int = (L_int - ((n_cols - 1) * d) - 2 * r_clad) / 2.0

        cx_int, cy_int = L_ext / 2.0, L_ext / 2.0
        
        # Graduation des centres
        start_x = W_start + l_gap_int + r_clad
        centers_x = [start_x + k * d for k in range(n_cols)]
        start_y = W_start + l_gap_int + r_clad
        centers_y = [start_y + k * d for k in range(n_cols)]

        # Bornes incluant les parois internes du boîtier
        bounds_x = [cx_int - L_int/2.0] + centers_x + [cx_int + L_int/2.0]
        bounds_y = [cy_int - L_int/2.0] + centers_y + [cy_int + L_int/2.0]
        
        return bounds_x[j-1], bounds_y[i-1], bounds_x[j], bounds_y[i]

    def _calculer_et_mettre_en_cache(self, x1, y1, x2, y2, z1, z2):
        """
        Internal method to calculate porosity, fluid area,
        hydraulic diameter, heating perimeter, inner box perimeter, and water rod perimeter for a given volume.
        
        Parameters:
        - (x1, y1): Bottom-left corner of the vertical prism
        - (x2, y2): Top-right corner of the vertical prism
        - z1: Starting axial coordinate of the prism
        - z2: Ending axial coordinate of the prism
        
        Returns:
        - poro_3d: Porosity of the volume (fluid volume / total volume)
        - a_cool: Total area of fluid in the volume (sum of fluid area across slices * slice thickness)
        - dh_3d: Hydraulic diameter of the fluid region (4 * fluid volume / wetted perimeter)
        - ph_moyen: Average heating perimeter in contact with the fluid across the axial height
        - pbox_moyen: Average inner box perimeter in contact with the fluid across the axial height
        - pwr_moyen: Average water rod perimeter in contact with the fluid across the axial height.
        """
        args_actuels = (round(x1, 5), round(y1, 5), round(x2, 5), round(y2, 5), round(z1, 5), round(z2, 5))
        if self._dernier_args == args_actuels:
            return self._dernier_resultats
            
        p, a_cool, dh, ph, pbox, pwr = analyse_3d_volume(self.slices_data, x1, y1, x2, y2, z1, z2)
        self._dernier_args = args_actuels
        self._dernier_resultats = (p, a_cool, dh, ph, pbox, pwr)
        return p, a_cool, dh, ph, pbox, pwr
    
    def _get_cylinders(self, dragon_assembly_model, r_wr_override=None):
        """
        Extrait la liste des cylindres (x, y, R) pour UNE tranche axiale spécifique.

        Parameters:
        - data: assembly geometry dict
        - r_wr_override: If provided, use this as the water rod outer radius (required for
          conical tranches where WATER_ROD_GEOMETRY has no fixed outer_radius).
        """
        cylinders = []
        d = dragon_assembly_model.pin_geometry_dict['pin_pitch']
        r_clad = dragon_assembly_model.pin_geometry_dict['clad_radius']
        
        if r_wr_override is not None:
            r_wr = r_wr_override
        elif dragon_assembly_model.water_rod_type == "circular":
            r_wr = dragon_assembly_model.water_rod_outer_radius
        elif dragon_assembly_model.water_rod_type == "conical":
            r_wr = getattr(dragon_assembly_model, "water_rod_outer_radius_start", 0.0)
        else:
            r_wr = 0.0

        grid_thickness = dragon_assembly_model.pin_geometry_dict.get('grid_thickness', 0.0)
        r_solid_clad = r_clad + grid_thickness
        r_solid_wr   = r_wr   + grid_thickness

        W_start, L_ext, L_int = _extract_box_dimensions(dragon_assembly_model)

        lattice = dragon_assembly_model.lattice_description
        exclusions = set(dragon_assembly_model.non_fuel_rod_ids)
        l_gap_int = (L_int - ((len(lattice[0]) - 1) * d) - 2 * r_clad) / 2.0

        # 1. Crayons Combustibles (Fuel rods)
        for row_idx, row in enumerate(lattice):
            for col_idx, item in enumerate(row):
                cx = W_start + l_gap_int + r_clad + col_idx * d
                cy = W_start + l_gap_int + r_clad + row_idx * d
                
                if item == 'VROD':
                    if grid_thickness > 0.0:
                        cylinders.append((cx, cy, r_solid_clad))
                    continue

                if item not in exclusions:
                    cylinders.append((cx, cy, r_solid_clad))

        # 2. Tubes d'eau (Water rods)
        for water_rod in dragon_assembly_model.water_rods:
            center = water_rod.center
            cylinders.append((center[0], center[1], r_solid_wr))

        return cylinders
    
    def _get_active_y_bounds(self, x_slice):
        """
        Calculate the active Y bounds (y values within the inner box) for a given X slice in the YZ plane.

        Parameters:
        - x_slice: the x-coordinate of the slice in the YZ plane

        Returns:
        - (y_min, y_max): the minimum and maximum y-coordinates of the active region for this x_slice
        """
        W_start, W_end, L_ext, L_int = self._get_box_geom()
        R_c = self.data_ref['ASSEMBLY_GEOMETRY']['corner_inner_radius_of_curvature']
        
        # 1. Complètement en dehors du boîtier
        if x_slice <= W_start or x_slice >= W_end:
            return None, None
        # 2. Zone centrale (bords droits)
        if W_start + R_c <= x_slice <= W_end - R_c:
            return W_start, W_end
        # 3. Dans l'arrondi gauche
        if x_slice < W_start + R_c:
            # Distance au centre du coin
            dx = (W_start + R_c) - x_slice 
            dy = math.sqrt(max(0, R_c**2 - dx**2)) 
            return W_start + R_c - dy, W_end - R_c + dy
        # 4. Dans l'arrondi droit
        if x_slice > W_end - R_c:
            # Distance au centre du coin
            dx = x_slice - (W_end - R_c) 
            dy = math.sqrt(max(0, R_c**2 - dx**2)) 
            return W_start + R_c - dy, W_end - R_c + dy
        return None, None
    
    def _get_active_x_bounds(self, y_slice):
        """
        Calculate the active X bounds (x values within the inner box) for a given Y slice in the XZ plane.

        Parameters:
        - y_slice: the y-coordinate of the slice in the XZ plane

        Returns:
        - (x_min, x_max): the minimum and maximum x-coordinates of the active region for this y_slice
        """
        W_start, W_end, L_ext, L_int = self._get_box_geom()
        R_c = self.data_ref['ASSEMBLY_GEOMETRY']['corner_inner_radius_of_curvature']
        
        # 1. Complètement en dehors du boîtier
        if y_slice <= W_start or y_slice >= W_end:
            return None, None
        # 2. Zone centrale (bords droits)
        if W_start + R_c <= y_slice <= W_end - R_c:
            return W_start, W_end
        # 3. Dans l'arrondi bas
        if y_slice < W_start + R_c:
            # Distance au centre du coin en Y
            dy = (W_start + R_c) - y_slice 
            dx = math.sqrt(max(0, R_c**2 - dy**2)) 
            return W_start + R_c - dx, W_end - R_c + dx
        # 4. Dans l'arrondi haut
        if y_slice > W_end - R_c:
            # Distance au centre du coin en Y
            dy = y_slice - (W_end - R_c)
            dx = math.sqrt(max(0, R_c**2 - dy**2)) 
            return W_start + R_c - dx, W_end - R_c + dx
        return None, None
    
    def get_wr_hole_data(self):
        """
        Extract water rod hole positions/areas and Idelchik loss tables from the core YAML.

        Returns:
        - holes   : list of dicts {z (cm), A_hole (cm²)} sorted by z
        - f_table : list of [winf_over_w0, K] pairs — water exits WR (P_actif < P_WR)
        - g_table : list of [winf_over_w0, K] pairs — water enters WR (P_actif > P_WR)
        """
        model = getattr(self.core_model, 'wr_hole_model', {}) 
        f_table = model.get('IDELCHIK-EXIT', [])
        g_table = model.get('IDELCHIK-ENTER', [])

        # n_wr: number of water rods — read from the first tranche that has it
        n_wr = next((t.get('n_wr', 1) for t in self.slices_data if t.get('n_wr', 0) > 0), 1)
        holes = []
        layouts = getattr(self.core_model, 'assembly_axial_layouts', {})
        regions = layouts.get(self.core_assembly_name, [])
        for region in regions:
            for hole in region.get('wr_holes', []):
                # Each wr_holes entry describes one hole per WR at this altitude.
                # Total equivalent area = n_wr * pi*(D/2)^2, represented as
                # a single equivalent hole with D_eq = D * sqrt(n_wr).
                d_single = float(hole['D_hole'])
                d_eq = d_single * math.sqrt(n_wr)
                holes.append({'z': float(hole['z']), 'D_hole': d_eq})
        holes.sort(key=lambda h: h['z'])
        return holes, f_table, g_table

    def get_z_global_bounds(self):
        """ 
        Get the global minimum and maximum axial coordinates (z) across all slices.
        
        Parameters: None
        
        Returns:
        - z_min: Global minimum axial coordinate
        - z_max: Global maximum axial coordinate
        """
        z_min = min(t['z_start'] for t in self.slices_data)
        z_max = max(t['z_end'] for t in self.slices_data)
        return z_min, z_max
    
    def get_x_global_bounds(self):
        return 0.0, self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']

    def get_y_global_bounds(self):
        return 0.0, self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']

# --- Z-DIRECTION: PUBLIC METHODS FOR POROSITY, HYDRAULIC DIAMETER, HEATING PERIMETER, WATER ROD PERIMETER AND INNER BOX PERIMETER  ---

# Control Volume (CV)
    def get_porosity_z_cv(self, x1, y1, x2, y2, z1, z2): 
        """
        Public method to get porosity in the z-direction for a control volume 
        defined by (x1, y1, z1) to (x2, y2, z2).
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume 
        
        Returns:
        - porosity: Porosity in the z-direction of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[0]
    
    def get_a_cool_z_cv(self, x1, y1, x2, y2, z1, z2): 
        """
        Public method to get total fluid area in the z-direction for a control volume 
        defined by (x1, y1, z1) to (x2, y2, z2).
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume 
        
        Returns:
        - a_cool: Total area of fluid in the z-direction of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[1]
    
    def get_dh_z_cv(self, x1, y1, x2, y2, z1, z2): 
        """Public method to get hydraulic diameter in the z-direction for a control volume 
        defined by (x1, y1, z1) to (x2, y2, z2). 
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume 
        
        Returns:
        - hydraulic_diameter: Hydraulic diameter in the z-direction of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[2]
    
    def get_ph_cv(self, x1, y1, x2, y2, z1, z2): 
        """Public method to get average heating perimeter for a control volume 
        defined by (x1, y1, z1) to (x2, y2, z2). 
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume 
        
        Returns:
        - average_heating_perimeter: Average heating perimeter of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[3]
    
    def get_pbox_cv(self, x1, y1, x2, y2, z1, z2): 
        """
        Public method to get average inner box perimeter for a control volume
        defined by (x1, y1, z1) to (x2, y2, z2).
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume
        
        Returns:
        - average_inner_box_perimeter: Average inner box perimeter of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[4]
    
    def get_pwr_cv(self, x1, y1, x2, y2, z1, z2): 
        """
        Public method to get average water rod perimeter for a control volume
        defined by (x1, y1, z1) to (x2, y2, z2).
        
        Parameters:
        - (x1, y1): Bottom-left corner of the control volume
        - (x2, y2): Top-right corner of the control volume
        - z1: Starting axial coordinate of the control volume
        - z2: Ending axial coordinate of the control volume
        
        Returns:
        - average_water_rod_perimeter: Average water rod perimeter of the specified control volume
        """
        return self._calculer_et_mettre_en_cache(x1, y1, x2, y2, z1, z2)[5]

# Water Channel
    def get_porosity_z_water(self, i, j, z1, z2):
        """Public method to get porosity in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - porosity: Porosity in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_porosity_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_a_cool_z_water(self, i, j, z1, z2):
        """Public method to get total fluid area in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - a_cool: Total area of fluid in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_dh_z_water(self, i, j, z1, z2):
        """Public method to get hydraulic diameter in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - hydraulic_diameter: Hydraulic diameter in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_dh_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_ph_water(self, i, j, z1, z2):
        """Public method to get average heating perimeter for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - average_heating_perimeter: Average heating perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_ph_cv(x1, y1, x2, y2, z1, z2)
    
    def get_pbox_water(self, i, j, z1, z2):
        """
        Public method to get average inner box perimeter for a water channel
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - average_inner_box_perimeter: Average inner box perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_pbox_cv(x1, y1, x2, y2, z1, z2)
    
    def get_pwr_water(self, i, j, z1, z2):
        """
        Public method to get average water rod perimeter for a water channel
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the water channel
        - z2: Ending axial coordinate of the water channel
        
        Returns:
        - average_water_rod_perimeter: Average water rod perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_water(i, j)
        return self.get_pwr_cv(x1, y1, x2, y2, z1, z2)
    
# Rod (Fuel Rod)
    def get_porosity_z_rod(self, i, j, z1, z2):
        """
        Public method to get porosity in the z-direction for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod
        
        Returns:
        - porosity: Porosity in the z-direction of the specified fuel rod
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_porosity_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_a_cool_z_rod(self, i, j, z1, z2):
        """
        Public method to get total fluid area in the z-direction for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod

        Returns:
        - a_cool: Total area of fluid in the z-direction of the specified fuel rod
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_dh_z_rod(self, i, j, z1, z2):
        """
        Public method to get hydraulic diameter in the z-direction for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod

        Returns:
        - hydraulic_diameter: Hydraulic diameter in the z-direction of the specified fuel rod
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_dh_z_cv(x1, y1, x2, y2, z1, z2)

    def get_ph_rod(self, i, j, z1, z2):
        """
        Public method to get average heating perimeter for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod

        Returns:
        - average_heating_perimeter: Average heating perimeter of the specified fuel rod
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_ph_cv(x1, y1, x2, y2, z1, z2)
    
    def get_pbox_rod(self, i, j, z1, z2):
        """
        Public method to get average inner box perimeter for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod    

        Returns:
        - average_inner_box_perimeter: Average inner box perimeter of the specified fuel rod    
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_pbox_cv(x1, y1, x2, y2, z1, z2)
    
    def get_pwr_rod(self, i, j, z1, z2):
        """
        Public method to get average water rod perimeter for a fuel rod
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is left column)
        - z1: Starting axial coordinate of the rod
        - z2: Ending axial coordinate of the rod

        Returns:
        - average_water_rod_perimeter: Average water rod perimeter of the specified fuel rod
        """
        x1, y1, x2, y2 = self._obtenir_bornes_rod(i, j)
        return self.get_pwr_cv(x1, y1, x2, y2, z1, z2)

#--- AXIAL Z-PROFILE AND AXIAL MESH MATRICES ---#

    def _get_acool_at_z(self, section_type, tranche, at_z):
        """
        Returns the fluid cross-sectional area AT a specific axial position z = at_z.

        For cylindrical tranches the area is uniform in z, so this equals the average over
        the tranche.  For conical tranches the water-rod radius is interpolated to at_z
        and the cross-section is evaluated analytically at that single elevation.

        This is the correct quantity to use when computing KSING at a tranche boundary z_b:
        evaluate each adjacent tranche AT z_b (not averaged over the tranche) so that a
        smooth cone connecting to a cylinder gives KSING = 0 by construction.

        Parameters:
        - section_type: same format as execute_profile_z
        - tranche: a tranche dict from self.tranches
        - at_z: axial position at which to evaluate the area [same units as z_start/z_end]

        Returns:
        - a_cool: fluid cross-sectional area at z = at_z  (cm² if geometry is in cm)
        """
        if tranche.get('is_conical', False):
            z_s, z_e = tranche['z_start'], tranche['z_end']
            frac = max(0.0, min(1.0, (at_z - z_s) / (z_e - z_s))) if z_e > z_s else 0.5
            r = tranche['r_cone_start'] + (tranche['r_cone_end'] - tranche['r_cone_start']) * frac
            data = tranche['dragon_assembly_model']
            ib, sol, lch, lwr, _ = build_assembly_geometry(data, r_wr_override=r)
            if section_type[0] == 'cv':
                x1, y1, x2, y2 = section_type[1]
            elif section_type[0] == 'rod':
                x1, y1, x2, y2 = self._obtenir_bornes_rod(*section_type[1])
            elif section_type[0] == 'water':
                x1, y1, x2, y2 = self._obtenir_bornes_water(*section_type[1])
            else:
                return 0.0
            return analyse_mesh(x1, y1, x2, y2, ib, sol, lch, lwr)[1]
        else:
            # Cylindrical: area is constant → average over tranche = point value
            z1, z2 = tranche['z_start'], tranche['z_end']
            if section_type[0] == 'cv':
                x1, y1, x2, y2 = section_type[1]
                return self.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
            elif section_type[0] == 'rod':
                i, j = section_type[1]
                return self.get_a_cool_z_rod(i, j, z1, z2)
            elif section_type[0] == 'water':
                i, j = section_type[1]
                return self.get_a_cool_z_water(i, j, z1, z2)
            elif section_type[0] == 'wr_tube':
                r_inner = self._get_r_inner_at_z(tranche, at_z)
                n_wr = tranche.get('n_wr', 0)
                return math.pi * r_inner ** 2 * n_wr
            return 0.0

    def _get_r_inner_at_z(self, tranche, at_z):
        """
        Returns the water-rod flow radius at a specific axial position for 'wr_tube'.

        For cylindrical tranches: constant inner bore radius (from 'r_wr_inner').
        For conical tranches: the outer radius is the flow boundary (no cladding in
        plenum transitions), linearly interpolated between r_cone_start and r_cone_end.

        Parameters:
        - tranche: tranche dict from self.tranches
        - at_z: axial position [cm]

        Returns:
        - r: flow radius [cm]
        """
        if tranche.get('is_conical', False):
            z_s, z_e = tranche['z_start'], tranche['z_end']
            frac = max(0.0, min(1.0, (at_z - z_s) / (z_e - z_s))) if z_e > z_s else 0.5
            r_i_s = tranche.get('r_cone_inner_start', 0.0)
            r_i_e = tranche.get('r_cone_inner_end',   0.0)
            return r_i_s + (r_i_e - r_i_s) * frac
        return tranche.get('r_wr_inner', 0.0)

    def get_a_cool_wr_tube(self, z1, z2):
        """
        Average flow area inside all water rod tubes over the axial interval [z1, z2].

        A(z) = π * r(z)² * n_wr, where r(z) is the flow radius:
        - Cylindrical tranche: r = r_inner (constant)
        - Conical tranche: r(z) linearly interpolated from r_cone_start to r_cone_end;
          the integral ∫ π*r(z)² dz is evaluated analytically.

        Parameters:
        - z1: Lower axial bound [cm]
        - z2: Upper axial bound [cm]

        Returns:
        - a_cool: Average flow area [cm²]
        """
        span = z2 - z1
        if span <= 0:
            return 0.0
        total = 0.0
        for tranche in self.slices_data:
            z_s, z_e = tranche['z_start'], tranche['z_end']
            lo = max(z_s, z1)
            hi = min(z_e, z2)
            if hi <= lo:
                continue
            n_wr = tranche.get('n_wr', 0)
            if tranche.get('is_conical', False):
                # r(z) = r_s + m*(z - z_s), with m = (r_e - r_s)/(z_e - z_s)
                # Use inner radius for flow area
                r_s = tranche.get('r_cone_inner_start', 0.0)
                r_e = tranche.get('r_cone_inner_end',   0.0)
                dz_t = z_e - z_s
                m = (r_e - r_s) / dz_t if dz_t > 0 else 0.0
                u_lo, u_hi = lo - z_s, hi - z_s
                r_lo = r_s + m * u_lo
                r_hi = r_s + m * u_hi
                if abs(m) < 1e-12:
                    integral = math.pi * r_s ** 2 * (u_hi - u_lo)
                else:
                    integral = math.pi / (3.0 * m) * (r_hi ** 3 - r_lo ** 3)
                total += integral * n_wr
            else:
                r = tranche.get('r_wr_inner', 0.0)
                total += math.pi * r ** 2 * n_wr * (hi - lo)
        return total / span

    def get_dh_wr_tube(self, z1, z2):
        """
        Average hydraulic diameter of the water rod tube interior over [z1, z2].

        For a circular bore: DH(z) = 2 * r(z).  Averaged as ∫ 2*r(z) dz / (z2-z1).
        For conical tranches the average is simply the mean of the endpoint diameters
        over the overlap interval.

        Parameters:
        - z1: Lower axial bound [cm]
        - z2: Upper axial bound [cm]

        Returns:
        - dh: Average hydraulic diameter [cm]
        """
        span = z2 - z1
        if span <= 0:
            return 0.0
        total_r = 0.0
        for tranche in self.slices_data:
            z_s, z_e = tranche['z_start'], tranche['z_end']
            lo = max(z_s, z1)
            hi = min(z_e, z2)
            if hi <= lo:
                continue
            if tranche.get('is_conical', False):
                # ∫ r(z) dz = (r_lo + r_hi)/2 * (hi - lo), using inner radius
                r_s = tranche.get('r_cone_inner_start', 0.0)
                r_e = tranche.get('r_cone_inner_end',   0.0)
                dz_t = z_e - z_s
                m = (r_e - r_s) / dz_t if dz_t > 0 else 0.0
                r_lo = r_s + m * (lo - z_s)
                r_hi = r_s + m * (hi - z_s)
                total_r += (r_lo + r_hi) / 2.0 * (hi - lo)
            else:
                r = tranche.get('r_wr_inner', 0.0)
                total_r += r * (hi - lo)
        return 2.0 * total_r / span

    def get_pch_wr_outer(self, z1, z2):
        """
        Average outer perimeter of all water rod tubes over [z1, z2].

        P_outer(z) = N_WR * 2 * pi * r_outer(z)

        For conical tranches, r_outer(z) is interpolated from r_cone_start to r_cone_end.
        For cylindrical tranches, r_outer = r_wr_outer (constant).

        Parameters:
        - z1, z2: axial bounds [cm]

        Returns:
        - pch_outer: average outer perimeter [cm]
        """
        span = z2 - z1
        if span <= 0:
            return 0.0
        total_r = 0.0
        for tranche in self.slices_data:
            z_s, z_e = tranche['z_start'], tranche['z_end']
            lo = max(z_s, z1)
            hi = min(z_e, z2)
            if hi <= lo:
                continue
            n_wr = tranche.get('n_wr', 0)
            if tranche.get('is_conical', False):
                r_s = tranche.get('r_cone_start', 0.0)
                r_e = tranche.get('r_cone_end',   0.0)
                dz_t = z_e - z_s
                m = (r_e - r_s) / dz_t if dz_t > 0 else 0.0
                r_lo = r_s + m * (lo - z_s)
                r_hi = r_s + m * (hi - z_s)
                total_r += n_wr * (r_lo + r_hi) / 2.0 * (hi - lo)
            else:
                r = tranche.get('r_wr_outer', 0.0)
                total_r += n_wr * r * (hi - lo)
        return 2.0 * math.pi * total_r / span

    def get_rwall_wr_tube(self, z1, z2):
        """
        Average wall conduction resistance of all water rod tubes over [z1, z2].

        R_wall(z) = ln(r_outer(z) / r_inner(z)) / (2 * pi * k_wall * N_WR)  [K*m/W]

        wall_conductivity k_wall [W/m/K] is read from the tranche (stored at parse time
        from WATER_ROD_GEOMETRY.wall_conductivity in the assembly YAML).

        For conical tranches, both r_outer and r_inner are linearly interpolated;
        the midpoint is used as a representative value for the node.

        Parameters:
        - z1, z2: axial bounds [cm]

        Returns:
        - r_wall: average wall resistance [K*m/W]
        """
        span = z2 - z1
        if span <= 0:
            return 0.0
        total_r = 0.0
        for tranche in self.slices_data:
            z_s, z_e = tranche['z_start'], tranche['z_end']
            lo = max(z_s, z1)
            hi = min(z_e, z2)
            if hi <= lo:
                continue
            n_wr = tranche.get('n_wr', 0)
            k_wall = tranche.get('k_wall', 18.0) 
            if n_wr == 0 or k_wall <= 0:
                continue
            z_mid = (lo + hi) / 2.0
            if tranche.get('is_conical', False):
                dz_t = z_e - z_s
                frac = (z_mid - z_s) / dz_t if dz_t > 0 else 0.5
                r_o = tranche.get('r_cone_start', 0.0) + frac * (
                    tranche.get('r_cone_end', 0.0) - tranche.get('r_cone_start', 0.0))
                r_i = tranche.get('r_cone_inner_start', 0.0) + frac * (
                    tranche.get('r_cone_inner_end', 0.0) - tranche.get('r_cone_inner_start', 0.0))
            else:
                r_o = tranche.get('r_wr_outer', 0.0)
                r_i = tranche.get('r_wr_inner', 0.0)
            if r_i > 1e-12 and r_o > r_i:
                r_wall_node = math.log(r_o / r_i) / (2.0 * math.pi * k_wall * n_wr)
            else:
                r_wall_node = 0.0
            total_r += r_wall_node * (hi - lo)
        return total_r / span

    @staticmethod
    def _conical_expansion_ksing(r_inlet, r_outlet, cone_length):
        """
        Singular pressure-loss coefficient for a conical expansion of the water rod bore.

        Model (from hydraulic handbooks):
          b     = d1 / d2      (inlet/outlet diameter ratio, b < 1 for expansion)
          alpha = arctan((d2 - d1) / (2 * L))  [degrees, half-angle of the cone]

          For a contraction (r_outlet < r_inlet):  K = 0  (conical contraction is lossless)

          0° ≤ alpha ≤ 10°:
            K = 8.3 * tan(alpha)^1.75 * (1 - b²)²

          10° < alpha ≤ 30°:
            K_base = 1.366 * sin_deg(2 * sqrt(2*alpha - 15)) - 0.17
            if b < 0.5: K_base -= 3.28 * (0.0625 - b⁴) * sqrt((alpha-10)/20)
            K = max(0, K_base) * (1 - b²)²

          30° < alpha ≤ 90°:
            if b <= 0.5:
              K = max(0, 1.205 - 3.28*(0.0625-b⁴) - 12.8*b⁶*sqrt((alpha-30)/60)) * (1-b²)²
            else:
              K = max(0, 1.205 - 0.2*sqrt((alpha-30)/60)) * (1-b²)²

        Parameters:
        - r_inlet:      inlet radius of the cone [any consistent unit]
        - r_outlet:     outlet radius of the cone [same unit]
        - cone_length:  axial length of the cone [same unit]

        Returns:
        - K: singular loss coefficient (dimensionless, ≥ 0)
        """
        if r_outlet <= r_inlet:
            return 0.0  # conical contraction: no singular loss
        d1 = 2.0 * r_inlet
        d2 = 2.0 * r_outlet
        b = d1 / d2
        alpha = math.degrees(math.atan((d2 - d1) / (2.0 * cone_length)))
        factor = (1.0 - b ** 2) ** 2

        def sin_deg(x):
            """sin with argument in degrees."""
            return math.sin(math.radians(x))

        if alpha <= 20.0:
            K = 3.2 * math.tan(math.radians(alpha)) ** 1.25 * factor
        elif alpha <= 40.0:
            K=(0.905-0.295*math.cos(math.radians(9.0*alpha/2.0)))*factor
        else:  # 40° < alpha <= 90°
            K = (1+0.2*(1+math.cos(math.radians(18.0*(alpha-40.0)/5.0)))/2.0) * factor
        return K

    def execute_profile_z(self, section_type, h, p, z_min, z_max):
        """
        Execute an axial profile analysis for the specified control volume, rod, or water channel.
        It computes porosity, fluid area, hydraulic diameter,
        and heating perimeter in the z direction at multiple axial positions.
        It also computes the singular pressure loss coefficient KSING at each node,
        which is non-zero only when a sudden expansion or contraction (tranche boundary
        from the YAML geometry) falls within the node's sampling window [curr_z, curr_z+h).
        KSING[k] is computed from the exact ACOOL of the two adjacent tranches,
        independently of the sampling step h.

        Parameters:
        - section_type: Tuple indicating the type of profile and its parameters
            - For control volume: ('cv', (x1, y1, x2, y2))
            - For rod: ('rod', (i, j))
            - For water channel: ('water', (i, j))
            - For water rod tube interior: ('wr_tube',)
        - h: Axial step size for sampling the profile
        - p: Axial step size for moving to the next sample point (can be equal to h or different for overlapping samples)
        - z_min: Minimum axial coordinate to start the profile
        - z_max: Maximum axial coordinate to end the profile

        Returns:
        - z_coords: List of axial coordinates where the properties were computed
        - porosities: List of porosity values corresponding to z_coords
        - a_cools: List of fluid area values corresponding to z_coords
        - dhs: List of hydraulic diameter values corresponding to z_coords
        - phs: List of heating perimeter values corresponding to z_coords
        - kexp_list: Singular expansion loss coefficients (dimensionless) at each node.
            kexp_list[k] > 0 when a sudden expansion boundary falls within [curr_z, curr_z+h).
            Borda-Carnot: K_exp = (1 - A_before/A_after)^2, reference area = a_before.
        - kcon_list: Singular contraction loss coefficients (dimensionless) at each node.
            kcon_list[k] > 0 when a sudden contraction boundary falls within [curr_z, curr_z+h).
            Borda-Carnot: K_con = 0.42*(1-r) if sqrt(r)<=0.76, else (1-r)^2, r=A_after/A_before.
            Reference area = a_after.
        - rsin_list: Area ratio A_small/A_large for the contraction(s) in each node.
            Weighted average of r over all contractions in the node (weighted by k_phys).
            rsin_list[k] = 1.0 when there is no contraction in node k.
        - cible_str: String identifier of the target (e.g. "cv_x1_y1_x2_y2", "rod_i_j", "water_i_j") for labeling purposes
        """
        z_coords = []
        porosities, a_cools, dhs, phs, kexp_list, kcon_list, rsin_list = [], [], [], [], [], [], []

        # Build list of all inter-tranche boundaries within [z_min, z_max].
        # Each entry: (index of tranche before boundary, z_boundary).
        # Conical regions are now stored as single tranches, so there are no spurious
        # sub-tranche boundaries to filter.  The transition PLENUM→DOM automatically
        # yields KSING≈0 because _get_acool_at_z evaluates the area AT z_b (cone tip
        # radius = nominal radius) for both adjacent tranches.
        tranche_bounds = [
            (k, float(self.slices_data[k]['z_end']))
            for k in range(len(self.slices_data) - 1)
            if z_min < float(self.slices_data[k]['z_end']) < z_max
        ]

        curr_z = z_min
        while curr_z + h <= z_max:
            z1 = curr_z
            z2 = curr_z + h

            if section_type[0]=='cv':
                x1, y1, x2, y2 = section_type[1]
                phi = self.get_porosity_z_cv(x1, y1, x2, y2, z1, z2)
                a_cool = self.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
                dh = self.get_dh_z_cv(x1, y1, x2, y2, z1, z2)
                ph = self.get_ph_cv(x1, y1, x2, y2, z1, z2)
                cible_str = f"cv_{x1}_{y1}_{x2}_{y2}"
            elif section_type[0]=='rod':
                i, j = section_type[1]
                phi = self.get_porosity_z_rod(i, j, z1, z2)
                a_cool = self.get_a_cool_z_rod(i, j, z1, z2)
                dh = self.get_dh_z_rod(i, j, z1, z2)
                ph = self.get_ph_rod(i, j, z1, z2)
                cible_str = f"rod_{i}_{j}"
            elif section_type[0]=='water':
                i, j = section_type[1]
                phi = self.get_porosity_z_water(i, j, z1, z2)
                a_cool = self.get_a_cool_z_water(i, j, z1, z2)
                dh = self.get_dh_z_water(i, j, z1, z2)
                ph = self.get_ph_water(i, j, z1, z2)
                cible_str = f"water_{i}_{j}"
            elif section_type[0]=='wr_tube':
                # phi = 1.0 by definition (entire bore is fluid)
                phi = 1.0
                a_cool = self.get_a_cool_wr_tube(z1, z2)
                dh = self.get_dh_wr_tube(z1, z2)
                ph = 0.0  # water rod interior is unheated
                cible_str = "wr_tube"
            else:
                raise ValueError(f"Unknown section_type: {section_type[0]!r}")

            # Compute KEXP/KCON at each inter-tranche boundary within this node.
            ksing_exp = 0.0
            ksing_con = 0.0
            r_sing_num = 0.0  # weighted sum of r = A_small/A_large for contractions
            r_sing_den = 0.0  # sum of k_phys weights for r averaging

            if section_type[0] == 'wr_tube':
                # --- Conical expansion KSING: distributed uniformly over all nodes
                #     that overlap the cone (K_total * overlap / L_cone per node).
                #     This avoids assigning the full K to a single node at the cone start,
                #     which would give physically wrong pressure-drop profiles.
                for t in self.slices_data:
                    if t.get('is_conical', False):
                        z_cs = t['z_start']
                        z_ce = t['z_end']
                        overlap = min(z2, z_ce) - max(z1, z_cs)
                        if overlap > 1e-12:
                            L_cone = z_ce - z_cs
                            K_total = self._conical_expansion_ksing(
                                t.get('r_cone_inner_start', t.get('r_cone_start', 0.0)),
                                t.get('r_cone_inner_end',   t.get('r_cone_end', 0.0)),
                                L_cone)
                            ksing_exp += K_total * overlap / L_cone

                # --- Abrupt steps at boundaries (cone-end → cylinder, or cylinder → cylinder).
                #     Conical-start boundaries are skipped here (handled above).
                #     Use z1+p (not z2) so each boundary is counted in exactly one node,
                #     independent of the window size h.
                for (k_idx, z_b) in tranche_bounds:
                    if z1 <= z_b < z1 + p:
                        t_before = self.slices_data[k_idx]
                        t_after  = self.slices_data[k_idx + 1]
                        if t_after.get('is_conical', False):
                            continue  # distributed above, not a point loss
                        r_before = self._get_r_inner_at_z(t_before, z_b)
                        r_after  = t_after.get('r_wr_inner', 0.0)
                        if r_before > 1e-12 and r_after > 1e-12:
                            ratio_A = (r_after / r_before) ** 2
                            if abs(ratio_A - 1.0) > 1e-6:
                                if ratio_A > 1.0:
                                    ksing_exp += (1.0 - 1.0 / ratio_A) ** 2
                                else:
                                    k_c = (0.42 * (1.0 - ratio_A) if math.sqrt(ratio_A) <= 0.76
                                           else (1.0 - ratio_A) ** 2)
                                    ksing_con += k_c
                                    r_sing_num += ratio_A * k_c
                                    r_sing_den += k_c
            else:
                # Generic Borda-Carnot (cv, rod, water).
                # _get_acool_at_z evaluates area AT z_b so conical→cylindrical
                # transitions give the correct ratio automatically.
                # Use z1+p (not z2) so each boundary is counted in exactly one node,
                # independent of the window size h.
                #
                # Mesh-independence correction: THMPV (LKSING=.TRUE.) computes
                #   DP_SING = rho * (Q/ACOOL_K)^2 * KSING
                # where ACOOL_K is the volume-average area of the node.  To ensure
                # the physical pressure drop is the same regardless of how a spacer
                # aligns with the mesh, each K contribution is scaled by
                # (a_ref / a_cool)^2, where a_ref is the Borda-Carnot reference area:
                #   expansion   -> a_ref = a_before  (upstream)
                #   contraction -> a_ref = a_after   (downstream)
                # This guarantees DP = rho*(Q/a_ref)^2 * K_phys independent of a_cool.
                for (k_idx, z_b) in tranche_bounds:
                    if z1 <= z_b < z1 + p:
                        t_before = self.slices_data[k_idx]
                        t_after  = self.slices_data[k_idx + 1]
                        a_before = self._get_acool_at_z(section_type, t_before, z_b)
                        a_after  = self._get_acool_at_z(section_type, t_after,  z_b)
                        if a_before > 1e-12 and a_cool > 1e-12:
                            ratio_A = a_after / a_before
                            if abs(ratio_A - 1.0) > 1e-6:
                                if ratio_A > 1.0: # expansion: reference = upstream = a_before
                                    k_phys = (1.0 - 1.0 / ratio_A) ** 2
                                    a_ref = a_before
                                    ksing_exp += k_phys * (a_cool / a_ref) ** 2
                                else: # contraction: reference = downstream = a_after
                                    if math.sqrt(ratio_A) <= 0.76:
                                        k_phys = 0.42 * (1.0 - ratio_A)
                                    else:
                                        k_phys = (1.0 - ratio_A) ** 2
                                    a_ref = a_after
                                    ksing_con += k_phys * (a_cool / a_ref) ** 2
                                    r_sing_num += ratio_A * k_phys
                                    r_sing_den += k_phys

            z_coords.append(curr_z + h/2.0)
            porosities.append(phi)
            a_cools.append(a_cool)
            dhs.append(dh)
            phs.append(ph)
            r_sing = r_sing_num / r_sing_den if r_sing_den > 1e-12 else 1.0
            kexp_list.append(ksing_exp)
            kcon_list.append(ksing_con)
            rsin_list.append(r_sing)
            curr_z += p
        return z_coords, porosities, a_cools, dhs, phs, kexp_list, kcon_list, rsin_list, cible_str
    
    def execute_mesh_z(self, mesh_type, z1, z2):
        """
        Execute a mesh analysis for the entire assembly based on the 
        specified mesh type (rod, water, or regular).
        It computes porosity, fluid area, hydraulic diameter, and heating perimeter in the z direction
        for each cell in the mesh and organizes the results in matrices of x and y coordinates.
        
        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For rod: ('rod', None)
            - For water: ('water', None)
            - For regular: ('regular', (n))
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        
        Returns:
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        - mat_ph: 2D list (matrix) of heating perimeter values for each cell
        
        The structure of the matrices depends on the mesh type:
        - For 'rod': mat_p[i][j] corresponds to the rod at lattice position (i+1, j+1)
        - For 'water': mat_p[i][j] corresponds to the canal at lattice position (i+1, j+1)
        - For 'regular': mat_p[i][j] corresponds to the control volume defined by the regular grid cell
        """
        mat_p, mat_dh, mat_ph = [], [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])

        if mesh_type[0] == 'rod':
            for i in range(1, n_cols+1):
                row_p, row_dh, row_ph = [], [], []
                for j in range(1, n_cols + 1):
                    row_p.append(self.get_porosity_z_rod(i, j, z1, z2))
                    row_dh.append(self.get_dh_z_rod(i, j, z1, z2))
                    row_ph.append(self.get_ph_rod(i, j, z1, z2))
                mat_p.append(row_p); mat_dh.append(row_dh); mat_ph.append(row_ph)
                
        elif mesh_type[0] == 'water':
            for i in range(1, n_cols+2):
                row_p, row_dh, row_ph = [], [], []
                for j in range(1, n_cols + 2):
                    row_p.append(self.get_porosity_z_water(i, j, z1, z2))
                    row_dh.append(self.get_dh_z_water(i, j, z1, z2))
                    row_ph.append(self.get_ph_water(i, j, z1, z2))
                mat_p.append(row_p); mat_dh.append(row_dh); mat_ph.append(row_ph)
                
        elif mesh_type[0] == 'regular':
            W_start, _, _, L_int = self._get_box_geom()
            bounds = [W_start + k * L_int / mesh_type[1] for k in range(mesh_type[1] + 1)]
            
            for i in range(len(bounds)-1):
                row_p, row_dh, row_ph = [], [], []
                idx_y = i
                y1, y2 = bounds[idx_y], bounds[idx_y+1]
                for j in range(len(bounds)-1):
                    x1, x2 = bounds[j], bounds[j+1]
                    row_p.append(self.get_porosity_z_cv(x1, y1, x2, y2, z1, z2))
                    row_dh.append(self.get_dh_z_cv(x1, y1, x2, y2, z1, z2))
                    row_ph.append(self.get_ph_cv(x1, y1, x2, y2, z1, z2))
                mat_p.append(row_p); mat_dh.append(row_dh); mat_ph.append(row_ph)
        
        return mat_p, mat_dh, mat_ph
    
    def get_porosity_x_cv(self, x_slice, y1, y2, z1, z2):
        """
        Calculates the porosity in the x-direction for a slice defined by 
        x_slice and vertical bounds [y1, y2] across the axial range [z1, z2].

        Parameters:
        - x_slice: The x-coordinate of the slice in the YZ plane
        - y1: Minimum y-coordinate of the slice
        - y2: Maximum y-coordinate of the slice
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - porosity: The porosity in the x-direction for the specified slice
        """
        v_totale = 0.0
        v_fluide = 0.0
        
        hauteur_totale = z2 - z1
        if hauteur_totale <= 0 or y2 <= y1:
            return 0.0

        # 1. Calcul les limites de la boîte à cette coordonnée X
        box_y_min, box_y_max = self._get_active_y_bounds(x_slice)
        
        # Si la coupe est dans le by-pass
        if box_y_min is None:
            return 0.0

        # 2. On restreint notre fenêtre d'étude [y1, y2] à l'intérieur du boîtier
        min_y = max(y1, box_y_min)
        max_y = min(y2, box_y_max)
        
        w_tot_valide = max_y - min_y
        if w_tot_valide <= 1e-9:
            return 0.0 

        # 3. Boucle d'intégration axiale
        for tranche in self.slices_data:
            z_min_overlap = max(z1, tranche['z_start'])
            z_max_overlap = min(z2, tranche['z_end'])
            dz = z_max_overlap - z_min_overlap
            
            if dz > 0:
                cylinders = tranche['cylinders']
                w_bloque = 0.0
                
                # Calcul analytique des cordes
                for (xc, yc, R) in cylinders:
                    dx = abs(x_slice - xc)
                    if dx < R:
                        corde_y = 2.0 * math.sqrt(R**2 - dx**2)
                        
                        # Intersection de la corde avec l'espace actif [min_y, max_y]
                        inter_y_start = max(min_y, yc - (corde_y / 2.0))
                        inter_y_end = min(max_y, yc + (corde_y / 2.0))
                        
                        if inter_y_start < inter_y_end:
                            w_bloque += (inter_y_end - inter_y_start)
                
                v_totale += w_tot_valide * dz
                v_fluide += (w_tot_valide - w_bloque) * dz

        if v_totale <= 1e-9:
            return 0.0

        porosite = v_fluide / v_totale
        return round(porosite, 5)
    
    def get_dh_x_cv(self, x_slice, y1, y2, z1, z2):
        """
        Calculates the hydraulic diameter in the x-direction for a slice defined by 
        x_slice and vertical bounds [y1, y2] across the axial range [z1, z2].

        Parameters:
        - x_slice: The x-coordinate of the slice in the YZ plane
        - y1: Minimum y-coordinate of the slice
        - y2: Maximum y-coordinate of the slice
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - hydraulic_diameter: The hydraulic diameter in the x-direction for the specified slice
        """
        v_fluide = 0.0
        s_mouillee = 0.0
        
        hauteur_totale = z2 - z1
        if hauteur_totale <= 0 or y2 <= y1:
            return 0.0

        # 1. Calcul les limites de la boîte à cette coordonnée X
        box_y_min, box_y_max = self._get_active_y_bounds(x_slice)
        
        # Si la coupe est dans le by-pass
        if box_y_min is None:
            return 0.0

        # 2. On restreint notre fenêtre d'étude [y1, y2] à l'intérieur du boîtier
        min_y = max(y1, box_y_min)
        max_y = min(y2, box_y_max)
        
        w_tot_valide = max_y - min_y
        if w_tot_valide <= 1e-9:
            return 0.0

        # 3. Boucle d'intégration axiale
        for tranche in self.slices_data:
            z_min_overlap = max(z1, tranche['z_start'])
            z_max_overlap = min(z2, tranche['z_end'])
            dz = z_max_overlap - z_min_overlap
            
            if dz > 0:
                w_bloque = 0.0
                p_m_tranche = 0.0
                cylinders = tranche['cylinders']
                
                # Ajout des murs du boîtier interne dans le périmètre
                if box_y_min >= y1: p_m_tranche += 1.0
                if box_y_max <= y2: p_m_tranche += 1.0
                
                # Calcul analytique des cordes et des parois de crayons
                for (xc, yc, R) in cylinders:
                    dx = abs(x_slice - xc)
                    if dx < R:
                        corde_y = 2.0 * math.sqrt(R**2 - dx**2)
                        y_min_rod = yc - (corde_y / 2.0)
                        y_max_rod = yc + (corde_y / 2.0)
                        
                        inter_y_start = max(min_y, y_min_rod)
                        inter_y_end = min(max_y, y_max_rod)
                        
                        if inter_y_start < inter_y_end:
                            w_bloque += (inter_y_end - inter_y_start)
                            # Ajout des bords du crayon au périmètre mouillé
                            if y_min_rod >= min_y: p_m_tranche += 1.0
                            if y_max_rod <= max_y: p_m_tranche += 1.0
                
                v_fluide += (w_tot_valide - w_bloque) * dz
                s_mouillee += p_m_tranche * dz

        # Pour éviter la division par zéro
        if s_mouillee <= 1e-6:
            return 0.0

        dh = 4.0 * v_fluide / s_mouillee
        return round(dh, 5)

    def get_porosity_x_rod(self, i, x_slice, z1, z2):
        """
        Calculate the porosity in the x-direction for the horizontal band corresponding to the fuel rod line 'i'.

        Parameters:
        - i: Row index of the rod line (1 is bottom row)
        - x_slice: The x-coordinate of the slice in the YZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - porosity: The porosity in the x-direction for the specified rod line and slice
        """
        _, y1, _, y2 = self._obtenir_bornes_rod(i, 1)
        return self.get_porosity_x_cv(x_slice, y1, y2, z1, z2)

    def get_dh_x_rod(self, i, x_slice, z1, z2):
        """
        Calculate the hydraulic diameter in the x-direction for the horizontal band corresponding to the fuel rod line 'i'.

        Parameters:
        - i: Row index of the rod line (1 is bottom row)
        - x_slice: The x-coordinate of the slice in the YZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - hydraulic_diameter: The hydraulic diameter in the x-direction for the specified rod line and slice
        """
        _, y1, _, y2 = self._obtenir_bornes_rod(i, 1)
        return self.get_dh_x_cv(x_slice, y1, y2, z1, z2)

    def get_porosity_x_water(self, i, x_slice, z1, z2):
        """
        Calculate the porosity in the x-direction for the horizontal band corresponding to the water channel 'i'.

        Parameters:
        - i: Row index of the channel (1 is bottom row)
        - x_slice: The x-coordinate of the slice in the YZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - porosity: The porosity in the x-direction for the specified channel and slice
        """
        _, y1, _, y2 = self._obtenir_bornes_water(i, 1)
        return self.get_porosity_x_cv(x_slice, y1, y2, z1, z2)

    def get_dh_x_water(self, i, x_slice, z1, z2):
        """
        Calculate the hydraulic diameter in the x-direction for the horizontal band corresponding to the water channel 'i'.
        
        Parameters:
        - i: Row index of the channel (1 is bottom row)
        - x_slice: The x-coordinate of the slice in the YZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        
        Returns:
        - hydraulic_diameter: The hydraulic diameter in the x-direction for the specified channel and slice
        """
        _, y1, _, y2 = self._obtenir_bornes_water(i, 1)
        return self.get_dh_x_cv(x_slice, y1, y2, z1, z2)

    def execute_profile_x(self, section_type, z1, z2, p):
        """
        Execute a lateral x-profile (Cross-flow in X) for a given section.
        The YZ slice moves along the X-axis with a step size p.

        Parameters:
        - section_type: Tuple indicating the type of section and its parameters
            - For control volume: ('cv', (y1, z1, y2, z2))
            - For rod: ('rod', (i, j))
            - For water channel: ('water', (i, j))
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        - p: Step size for moving the slice along the X-axis

        Returns:
        - x_coords: List of x-coordinates where the properties were computed
        - porosities: List of porosity values corresponding to x_coords
        - dhs: List of hydraulic diameter values corresponding to x_coords
        - cible_str: String identifier of the target (e.g. "latX_cv_y1_z1_y2_z2", "latX_rod_i_j", "latX_water_i_j") for labeling purposes
        """
        x_coords = []
        porosities, dhs = [], []
        
        # 1. Détermination des bornes Y selon le type de section
        if section_type[0] == 'cv':
            y1, _, y2, _ = section_type[1]
            cible_str = f"cv_Y{round(y1,2)}_Y{round(y2,2)}"
        elif section_type[0] == 'rod':
            i, j = section_type[1]
            # On extrait uniquement les coordonnées Y de cette ligne
            _, y1, _, y2 = self._obtenir_bornes_rod(i, 1) 
            cible_str = f"rod_{i}_{j}"
        elif section_type[0] == 'water':
            i, j = section_type[1]
            # On extrait uniquement les coordonnées Y de ce canal
            _, y1, _, y2 = self._obtenir_bornes_water(i, 1)
            cible_str = f"water_{i}_{j}"
        else:
            raise ValueError(f"Type de section inconnu : {section_type[0]}")

        # 2. Balayage analytique sur l'axe X
        x_min, x_max = self.get_x_global_bounds()
        curr_x = x_min
        
        # On ajoute une petite tolérance (1e-9) pour les erreurs d'arrondi des floats
        while curr_x <= x_max + 1e-9:
            phi = self.get_porosity_x_cv(curr_x, y1, y2, z1, z2)
            dh = self.get_dh_x_cv(curr_x, y1, y2, z1, z2)
            
            x_coords.append(curr_x)
            porosities.append(phi)
            dhs.append(dh)
            
            curr_x += p
            
        return x_coords, porosities, dhs, cible_str
    
    def execute_mesh_x(self, mesh_type, x_slice, z1, z2, h, p):
        """
        Execute a mesh analysis in the x-direction for the entire assembly based on the
        specified mesh type (rod, water, or regular). It computes porosity and hydraulic diameter
        for each cell in the mesh and organizes the results in matrices of y and z coordinates.

        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For rod: ('rod', None)
            - For water: ('water', None)
            - For regular: ('regular', n_y) where n_y is the number of divisions in the y-direction
        - x_slice: The x-coordinate of the slice in the YZ plane where the properties are computed
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        - h: Height of each axial window (cm)
        - p: Step between successive windows (cm)

        Returns:         
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        """
        mat_p, mat_dh = [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_rows = len(lattice)
        
        # 1. Création des fenêtres axiales glissantes
        z_windows = []
        curr_z = z1
        while curr_z + h <= z2 + 1e-9:
            z_windows.append((curr_z, curr_z + h))
            curr_z += p
        
        # 2. Détermination des intervalles sur l'axe Y (les bandes horizontales)
        y_intervals = []
        
        if mesh_type[0] == 'regular':
            W_start, _, _, L_int = self._get_box_geom()
            # On découpe l'espace interne en mesh_type[1] parties égales
            n_y = mesh_type[1]
            y_bounds = [W_start + k * L_int / n_y for k in range(n_y + 1)]
            for j in range(n_y):
                y_intervals.append((y_bounds[j], y_bounds[j+1]))
        elif mesh_type[0] == 'rod':
            # On parcourt les lignes de crayons de haut en bas (i de n_rows à 1)
            for i in range(n_rows, 0, -1):
                _, y1_rod, _, y2_rod = self._obtenir_bornes_rod(i, 1) # j=1 suffit pour avoir Y
                y_intervals.append((y1_rod, y2_rod))
        elif mesh_type[0] == 'water':
            # On parcourt les canaux d'eau de haut en bas (il y a n_rows + 1 canaux Y)
            for i in range(n_rows + 1, 0, -1):
                _, y1_water, _, y2_water = self._obtenir_bornes_water(i, 1)
                y_intervals.append((y1_water, y2_water))
        else:
            raise ValueError(f"Type de maillage inconnu : {mesh_type[0]}")

        # 3. La double boucle de construction des matrices
        for (z1_cell, z2_cell) in z_windows:
            row_p, row_dh = [], []
            for (y1_cell, y2_cell) in y_intervals:
                # Appel de nos fonctions analytiques ultra-rapides
                phi = self.get_porosity_x_cv(x_slice, y1_cell, y2_cell, z1_cell, z2_cell)
                dh = self.get_dh_x_cv(x_slice, y1_cell, y2_cell, z1_cell, z2_cell)
                row_p.append(phi)
                row_dh.append(dh)
            mat_p.append(row_p)
            mat_dh.append(row_dh)
            
        return mat_p, mat_dh

    def get_porosity_y_cv(self, y_slice, x1, x2, z1, z2):
        """
        Calculates the porosity in the y-direction for a slice defined by
        y_slice and horizontal bounds [x1, x2] across the axial range [z1, z2].

        Parameters:
        - y_slice: The y-coordinate of the slice in the XZ plane
        - x1: Minimum x-coordinate of the slice
        - x2: Maximum x-coordinate of the slice
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - porosity: The porosity in the y-direction for the specified slice
        """
        v_totale = 0.0
        v_fluide = 0.0
        
        hauteur_totale = z2 - z1
        if hauteur_totale <= 0 or x2 <= x1:
            return 0.0

        box_x_min, box_x_max = self._get_active_x_bounds(y_slice)
        if box_x_min is None:
            return 0.0

        min_x = max(x1, box_x_min)
        max_x = min(x2, box_x_max)
        
        w_tot_valide = max_x - min_x
        if w_tot_valide <= 1e-9:
            return 0.0

        for tranche in self.slices_data:
            z_min_overlap = max(z1, tranche['z_start'])
            z_max_overlap = min(z2, tranche['z_end'])
            dz = z_max_overlap - z_min_overlap
            
            if dz > 0:
                w_bloque = 0.0
                cylinders = tranche['cylinders']
                for (xc, yc, R) in cylinders:
                    dy = abs(y_slice - yc) # Distance au centre du cercle selon Y
                    if dy < R:
                        corde_x = 2.0 * math.sqrt(R**2 - dy**2)
                        inter_x_start = max(min_x, xc - (corde_x / 2.0))
                        inter_x_end = min(max_x, xc + (corde_x / 2.0))
                        
                        if inter_x_start < inter_x_end:
                            w_bloque += (inter_x_end - inter_x_start)
                
                v_totale += w_tot_valide * dz
                v_fluide += (w_tot_valide - w_bloque) * dz

        if v_totale <= 1e-9: return 0.0
        return round(v_fluide / v_totale, 5)

    def get_dh_y_cv(self, y_slice, x1, x2, z1, z2):
        """
        Calculates the hydraulic diameter in the y-direction for a slice defined by
        y_slice and horizontal bounds [x1, x2] across the axial range [z1, z2]. 

        Parameters:
        - y_slice: The y-coordinate of the slice in the XZ plane
        - x1: Minimum x-coordinate of the slice
        - x2: Maximum x-coordinate of the slice
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice

        Returns:
        - hydraulic_diameter: The hydraulic diameter in the y-direction for the specified slice
        """
        v_fluide = 0.0
        s_mouillee = 0.0
        
        hauteur_totale = z2 - z1
        if hauteur_totale <= 0 or x2 <= x1:
            return 0.0

        box_x_min, box_x_max = self._get_active_x_bounds(y_slice)
        if box_x_min is None: return 0.0

        min_x = max(x1, box_x_min)
        max_x = min(x2, box_x_max)
        
        w_tot_valide = max_x - min_x
        if w_tot_valide <= 1e-9: return 0.0

        for tranche in self.slices_data:
            z_min_overlap = max(z1, tranche['z_start'])
            z_max_overlap = min(z2, tranche['z_end'])
            dz = z_max_overlap - z_min_overlap
            
            if dz > 0:
                w_bloque = 0.0
                p_m_tranche = 0.0
                cylinders = tranche['cylinders']
                
                # Ajout des murs Ouest/Est du boîtier interne
                if box_x_min >= x1: p_m_tranche += 1.0
                if box_x_max <= x2: p_m_tranche += 1.0
                
                for (xc, yc, R) in cylinders:
                    dy = abs(y_slice - yc)
                    if dy < R:
                        corde_x = 2.0 * math.sqrt(R**2 - dy**2)
                        x_min_rod = xc - (corde_x / 2.0)
                        x_max_rod = xc + (corde_x / 2.0)
                        
                        inter_x_start = max(min_x, x_min_rod)
                        inter_x_end = min(max_x, x_max_rod)
                        
                        if inter_x_start < inter_x_end:
                            w_bloque += (inter_x_end - inter_x_start)
                            # Parois Nord/Sud du crayon
                            if x_min_rod >= min_x: p_m_tranche += 1.0
                            if x_max_rod <= max_x: p_m_tranche += 1.0
                
                v_fluide += (w_tot_valide - w_bloque) * dz
                s_mouillee += p_m_tranche * dz

        if s_mouillee <= 1e-6: return 0.0
        return round(4.0 * v_fluide / s_mouillee, 5)
    
    def get_porosity_y_rod(self, j, y_slice, z1, z2):
        """
        Porosity in the y-direction for the vertical band corresponding to the fuel rod column 'j'.
        
        Parameters:
        - j: Column index of the rod line (1 is leftmost column)
        - y_slice: The y-coordinate of the slice in the XZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        
        Returns:
        - porosity: The porosity in the y-direction for the specified rod column and slice
        """
        x1, _, x2, _ = self._obtenir_bornes_rod(1, j) 
        return self.get_porosity_y_cv(y_slice, x1, x2, z1, z2)

    def get_dh_y_rod(self, j, y_slice, z1, z2):
        """
        Hydraulic diameter in the y-direction for the vertical band corresponding to the fuel rod column 'j'.
        
        Parameters:
        - j: Column index of the rod line (1 is leftmost column)
        - y_slice: The y-coordinate of the slice in the XZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        
        Returns:
        - hydraulic_diameter: The hydraulic diameter in the y-direction for the specified rod column and slice
        """
        x1, _, x2, _ = self._obtenir_bornes_rod(1, j)
        return self.get_dh_y_cv(y_slice, x1, x2, z1, z2)

    def get_porosity_y_water(self, j, y_slice, z1, z2):
        """
        Porosity in the y-direction for the vertical band corresponding to the water channel column 'j'.
        
        Parameters:
        - j: Column index of the water channel line (1 is leftmost column)
        - y_slice: The y-coordinate of the slice in the XZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        
        Returns:
        - porosity: The porosity in the y-direction for the specified water channel column and slice
        """
        x1, _, x2, _ = self._obtenir_bornes_water(1, j)
        return self.get_porosity_y_cv(y_slice, x1, x2, z1, z2)

    def get_dh_y_water(self, j, y_slice, z1, z2):
        """
        Hydraulic diameter in the y-direction for the vertical band corresponding to the water channel column 'j'.
        
        Parameters:
        - j: Column index of the water channel line (1 is leftmost column)
        - y_slice: The y-coordinate of the slice in the XZ plane
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        
        Returns:
        - hydraulic_diameter: The hydraulic diameter in the y-direction for the specified water channel column and slice
        """
        x1, _, x2, _ = self._obtenir_bornes_water(1, j)
        return self.get_dh_y_cv(y_slice, x1, x2, z1, z2)

    def execute_profile_y(self, section_type, z1, z2, p):
        """
        Execute a lateral y-profile for a given section.
        The XZ slice moves along the Y-axis with a step size p.

        Parameters:
        - section_type: Tuple indicating the type of section and its parameters
            - For control volume: ('cv', (x1, z1, x2, z2))
            - For rod: ('rod', (i, j))
            - For water channel: ('water', (i, j))
        - z1: Starting axial coordinate of the slice
        - z2: Ending axial coordinate of the slice
        - p: Step size for moving the slice along the Y-axis

        Returns:
        - y_coords: List of y-coordinates where the properties were computed
        - porosities: List of porosity values corresponding to y_coords
        - dhs: List of hydraulic diameter values corresponding to y_coords
        - cible_str: String identifier of the target (e.g. "latY_cv_x1_z1_x2_z2", "latY_rod_i_j", "latY_water_i_j") for labeling purposes
        """
        y_coords = []
        porosities, dhs = [], []
        
        if section_type[0] == 'cv':
            x1, _, x2, _ = section_type[1] # On extrait x1 et x2
            cible_str = f"cv_X{round(x1,2)}_X{round(x2,2)}"
        elif section_type[0] == 'rod':
            i, j = section_type[1]
            x1, _, x2, _ = self._obtenir_bornes_rod(1, j) 
            cible_str = f"rod_{i}_{j}"
        elif section_type[0] == 'water':
            i, j = section_type[1]
            x1, _, x2, _ = self._obtenir_bornes_water(1, j)
            cible_str = f"water_{i}_{j}"
        else:
            raise ValueError(f"Type de section inconnu : {section_type[0]}")

        _, y_max = self.get_y_global_bounds()
        curr_y = 0.0 # ou get_y_global_bounds()[0]
        
        while curr_y <= y_max + 1e-9:
            phi = self.get_porosity_y_cv(curr_y, x1, x2, z1, z2)
            dh = self.get_dh_y_cv(curr_y, x1, x2, z1, z2)
            
            y_coords.append(curr_y)
            porosities.append(phi)
            dhs.append(dh)
            curr_y += p
            
        return y_coords, porosities, dhs, cible_str

    def execute_mesh_y(self, mesh_type, y_slice, z1, z2, h, p):
        """
        Execute a 2D mesh on the XZ plane (Cross-Flow in Y) at a fixed y_slice coordinate.
        The mesh is constructed based on the specified mesh type (rod, water, or regular) and the assembly geometry.
        It computes porosity and hydraulic diameter for each cell in the mesh 
        and organizes the results in matrices of x and z coordinates.

        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For rod: ('rod', None)
            - For water: ('water', None)
            - For regular: ('regular', n_x) where n_x is the number of divisions in the x-direction
        - y_slice: The y-coordinate of the slice in the XZ plane where the properties are computed
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        - h: Height of each axial window (cm)
        - p: Step between successive windows (cm)

        Returns:
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        """
        mat_p, mat_dh = [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        
        # 1. Création des fenêtres axiales glissantes
        z_windows = []
        curr_z = z1
        while curr_z + h <= z2 + 1e-9:
            z_windows.append((curr_z, curr_z + h))
            curr_z += p

        x_intervals = []
        
        if mesh_type[0] == 'regular':
            W_start, _, _, L_int = self._get_box_geom()
            n_x = mesh_type[1]
            x_bounds = [W_start + k * L_int / n_x for k in range(n_x + 1)]
            for j in range(n_x):
                x_intervals.append((x_bounds[j], x_bounds[j+1]))
        elif mesh_type[0] == 'rod':
            for j in range(1, n_cols + 1):
                x1_rod, _, x2_rod, _ = self._obtenir_bornes_rod(1, j)
                x_intervals.append((x1_rod, x2_rod))
        elif mesh_type[0] == 'water':
            for j in range(1, n_cols + 2):
                x1_water, _, x2_water, _ = self._obtenir_bornes_water(1, j)
                x_intervals.append((x1_water, x2_water))

        for (z1_cell, z2_cell) in z_windows:
            row_p, row_dh = [], []
            for (x1_cell, x2_cell) in x_intervals:
                phi = self.get_porosity_y_cv(y_slice, x1_cell, x2_cell, z1_cell, z2_cell)
                dh = self.get_dh_y_cv(y_slice, x1_cell, x2_cell, z1_cell, z2_cell)
                row_p.append(phi)
                row_dh.append(dh)
            mat_p.append(row_p)
            mat_dh.append(row_dh)
            
        return mat_p, mat_dh

    def get_CD(self, Re, i, j):
        """"""
        p = self.data_ref['PIN_GEOMETRY']['pin_pitch']
        D = 2*self.data_ref['PIN_GEOMETRY']['clad_radius']
        A = 1.8 * (p/D -1)**(-0.5)
        U_x  = self.U_nominal_x[i][j - 1]
        U_y = self.U_nominal_y[i][j - 1]
        U = math.sqrt(U_x**2 + U_y**2)
        return p/D * U/U_x * A * Re**(-0.2) *(p/(p-D))**2
    
    def get_f(self, Re):
        """"""
        if Re < 2000:
            return 64/Re
        elif 2000 <= Re < 3000:
            return 1.07*10**(-5)*Re + 1.06*10**(-2)
        else:
            return 0.3164*Re**(-0.25)
        
    def get_dh_x_water_meca_flu(self, i, j):
        """
        Compute the hydraulic diameter for cross-flow in the x-direction for subchannel (i, j),
        using a momentum / friction-force balance:

            deltaP * py * pz = 0.5 * CD * rho * U^2 * D * pz   (drag on rod)
            deltaP = (px / 2Dh) * f * rho * U^2                 (Darcy-Weisbach)

        Solving: Dh = f * px * py / (CD * D)

        For interior subchannels px = py = p (pin pitch).
        For border subchannels px and/or py are smaller (half-gap to box wall);
        these are read directly from _obtenir_bornes_water so no special-casing is needed.

        Indexing convention: i is 0-based (0 = bottom border, n_rows = top border),
                             j is 0-based (0 = left border,  n_cols = right border).
        U_nominal_x[i][k] stores the cross-flow gap velocity between subchannel k and k+1,
        so the matrix shape is (n_rows+1) x n_cols.
        """
        D  = 2 * self.data_ref['PIN_GEOMETRY']['clad_radius']
        # U_nominal_x[i][k] is the cross-flow gap velocity between subchannel k and k+1.
        # The gap entering subchannel j from the left is therefore at index j-1.
        U_x  = self.U_nominal_x[i][j - 1]
        U_y = self.U_nominal_y[i][j - 1]
        U = math.sqrt(U_x**2 + U_y**2)
        
        mu = self.data_ref.get('FLUID_PROPERTIES', {}).get('dynamic_viscosity', 8.5e-5)
        rho = self.data_ref.get('FLUID_PROPERTIES', {}).get('density', 740.0)
        
        Re = rho * U * D / mu
        f  = self.get_f(Re)
        CD = self.get_CD(Re, i, j)
        
        # _obtenir_bornes_water uses 1-based indices; i and j are 0-based so +1
        x1_c, y1_c, x2_c, y2_c = self._obtenir_bornes_water(i + 1, j + 1)
        px = x2_c - x1_c   # width  of subchannel in x (= p for interior, < p on x-borders)
        py = y2_c - y1_c   # height of subchannel in y (= p for interior, < p on y-borders)
        return round(px * py * f / (CD * D), 5)