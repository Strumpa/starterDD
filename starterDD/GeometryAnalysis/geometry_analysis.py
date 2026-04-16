import os
import numpy as np
import yaml
from shapely.geometry import box, Point, MultiLineString
from shapely.ops import unary_union
import math

def build_global_geometry(data):
    """
    Build the full 2D reactor geometry for a single axial slice.

    Parameters:
    - data: Dictionary containing the geometry information for the assembly (from YAML)

    Returns: (inner_box, solide_total, multi_lignes_chauffantes)
    - inner_box: Polygon of the internal region (excluding the walls)
    - solide_total: Polygon of the solid region (fuel rods + water rods)
    - multi_lignes_chauffantes: MultiLineString of heating perimeters (e.g. rod boundaries)
    - multi_lignes_water: MultiLineString of water rod perimeters

    Used during AnalyseurGeometrique initialization for each axial slice.
    """

    pin_geo = data['PIN_GEOMETRY']
    ass_geo = data['ASSEMBLY_GEOMETRY']
    wr_geo = data['WATER_ROD_GEOMETRY']
    lattice = ass_geo['lattice_description']
    exclusions = set(ass_geo['non_fuel_rod_ids'])

    d = pin_geo['pin_pitch']
    r_clad = pin_geo['clad_radius']
    r_wr = wr_geo['outer_radius']
    
    L_ext = ass_geo['assembly_pitch']
    W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
    L_int = L_ext - 2 * W_start
    R_c = ass_geo['corner_inner_radius_of_curvature']

    l_gap_int = (L_int - ((len(lattice[0]) - 1) * d) - 2 * r_clad) / 2.0

    # Le Boitier Interne
    W_end = L_ext - W_start
    inner_box = box(W_start + R_c, W_start + R_c, W_end - R_c, W_end - R_c).buffer(R_c, resolution=64)

    formes_solides = []
    lignes_chauffantes = []
    lignes_water = []

    # 2. Water Rods
    for center in wr_geo['centers']:
        wr_circle = Point(center[0], center[1]).buffer(r_wr, resolution=64)
        formes_solides.append(wr_circle)
        lignes_water.append(wr_circle.exterior)
    
    # 3. Fuel Rods
    for row_idx, row in enumerate(lattice):
        for col_idx, item in enumerate(row):
            if item not in exclusions:
                cx = W_start + l_gap_int + r_clad + col_idx * d
                cy = W_start + l_gap_int + r_clad + row_idx * d
                
                rod_circle = Point(cx, cy).buffer(r_clad, resolution=64)
                formes_solides.append(rod_circle)
                lignes_chauffantes.append(rod_circle.exterior)

    # 4. Fusion
    solide_total = unary_union(formes_solides)
    multi_lignes_chauffantes = MultiLineString(lignes_chauffantes)
    multi_lignes_water = MultiLineString(lignes_water)

    return inner_box, solide_total, multi_lignes_chauffantes, multi_lignes_water

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
        
        if dz > 0:
            s_tot, s_m, p_m, p_h, p_box, p_wr = analyse_mesh(
                x1, y1, x2, y2, 
                tranche['inner_box'], tranche['solide'], tranche['lignes_chauffantes'], tranche['lignes_water']
            )
            v_total += s_tot * dz
            v_fluide += s_m * dz
            s_mouillee += p_m * dz
            s_chauffante += p_h * dz
            s_box += p_box * dz
            s_wr += p_wr * dz

    if v_total <= 1e-9:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    poro_3d = v_fluide / v_total
    a_cool = v_fluide / (z2-z1)
    dh_3d = (4 * v_fluide) / s_mouillee if s_mouillee > 0 else 0.0
    ph_moyen = s_chauffante / hauteur_totale
    pbox_moyen = s_box / hauteur_totale
    pwr_moyen = s_wr / hauteur_totale
    
    return round(poro_3d, 5), round(a_cool, 5), round(dh_3d, 5), round(ph_moyen, 5), round(pbox_moyen, 5), round(pwr_moyen, 5)


class GeometricAnalyser:
    """
    Class responsible for analyzing the geometry of a nuclear reactor assembly 
    based on YAML input files. 
    It provides methods to compute porosity, hydraulic diameter, 
    and heating perimeter for specified volumes or profiles within the assembly.
    """
    def __init__(self, core_yaml_path, core_i=1, core_j=1):
        """
        Initialize the geometric analyzer by reading the core YAML file,
        extracting the relevant axial slices for the specified assembly position.
        
        Parameters:
        - core_yaml_path: Path to the CORE YAML file containing the overall geometry
        - core_i: Row index of the assembly in the core layout (1 is bottom row)
        - core_j: Column index of the assembly in the core layout (1 is leftmost column)
        """
        self.tranches = []
        self.data_ref = None
        
        # 1. On lit le fichier CORE maître
        with open(core_yaml_path, 'r', encoding='utf-8') as f:
            core_data = yaml.safe_load(f)
            
        # 2. On lit la carte du cœur et on vérifie les indices
        carte_coeur = core_data['CORE_GEOMETRY']['core_2D_layout']
        max_i = len(carte_coeur)
        max_j = len(carte_coeur[0])
        
        if not (1 <= core_i <= max_i and 1 <= core_j <= max_j):
            raise ValueError(f"Erreur : Assemblage ({core_i},{core_j}) hors limites. Le cœur fait {max_i}x{max_j}.")
            
        nom_assemblage = carte_coeur[core_i - 1][core_j - 1]
        print(f"Info : Chargement de l'assemblage '{nom_assemblage}' à la position cœur ({core_i},{core_j})")
        
        # 3. On récupère l'empilement axial pour cet assemblage précis
        layouts = core_data['CORE_GEOMETRY']['assembly_axial_layouts']
        regions_axiales = layouts[nom_assemblage]
        
        dossier_core = os.path.dirname(os.path.abspath(core_yaml_path))
        
        # 4. On boucle sur les tranches (DOM, VAN, etc.)
        for region in regions_axiales:
            z_s, z_e = region['axial_bounds']
            chemin_relatif = region['assembly_geometry_file']
            
            # Reconstitution du chemin absolu
            chemin_absolu = os.path.abspath(os.path.join(dossier_core, chemin_relatif))
            
            with open(chemin_absolu, 'r', encoding='utf-8') as f:
                data_assemblage = yaml.safe_load(f)
                
                if self.data_ref is None: 
                    self.data_ref = data_assemblage
                    
                box_geom, solide, lignes_chauffantes, lignes_water = build_global_geometry(data_assemblage)
                cylinders_de_cette_tranche = self._get_cylinders(data_assemblage)
                
                # On utilise les z_s et z_e lus dans le CORE (ignore AXIAL_GEOMETRY de l'assemblage)
                self.tranches.append({
                    'z_start': float(z_s),
                    'z_end': float(z_e),
                    'inner_box': box_geom,
                    'solide': solide,
                    'lignes_chauffantes': lignes_chauffantes,
                    'lignes_water': lignes_water,
                    'cylinders': cylinders_de_cette_tranche
                })
        
        # Initialisation du cache
        self._dernier_args = None
        self._dernier_resultats = None

    def _obtenir_bornes_canal(self, i, j):
        """
        Calculate the bounding box of the water channel at position (i, j) in the lattice.
        This is used to define the control volume for canal-based queries.
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is leftmost column)
        Canals are defined on the border of the lattice too.
        
        Returns:
        - (x1, y1, x2, y2): Coordinates of the bounding box of the canal
        """
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        d = self.data_ref['PIN_GEOMETRY']['pin_pitch']
        r_clad = self.data_ref['PIN_GEOMETRY']['clad_radius']
        
        W_start = self.data_ref['ASSEMBLY_GEOMETRY']['gap_wide'] + self.data_ref['ASSEMBLY_GEOMETRY']['channel_box_thickness']
        L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
        W_end = L_ext - W_start
        L_int = L_ext - 2 * W_start
        l_gap_int = (L_int - ((n_cols - 1) * d) - 2 * r_clad) / 2.0

        first_center = W_start + l_gap_int + r_clad
        centers = [first_center + k * d for k in range(n_cols)]
        bounds = [W_start] + centers + [W_end]
        
        return bounds[j-1], bounds[i-1], bounds[j], bounds[i]

    def _obtenir_bornes_rod(self, i, j):
        """
        Calculate the bounding box of the fuel rod at position (i, j) in the lattice.
        
        Parameters:
        - i: Row index of the rod (1 is bottom row)
        - j: Column index of the rod (1 is leftmost column)
        
        Returns:
        - (x1, y1, x2, y2): Coordinates of the bounding box of the rod
        """

        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        d = self.data_ref['PIN_GEOMETRY']['pin_pitch']
        r_clad = self.data_ref['PIN_GEOMETRY']['clad_radius']
        
        W_start = self.data_ref['ASSEMBLY_GEOMETRY']['gap_wide'] + self.data_ref['ASSEMBLY_GEOMETRY']['channel_box_thickness']
        L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
        L_int = L_ext - 2 * W_start
        l_gap_int = (L_int - ((n_cols - 1) * d) - 2 * r_clad) / 2.0

        cx = W_start + l_gap_int + r_clad + (j - 1) * d
        cy = W_start + l_gap_int + r_clad + (i - 1) * d
        return cx - d/2.0, cy - d/2.0, cx + d/2.0, cy + d/2.0

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
            
        p, a_cool, dh, ph, pbox, pwr = analyse_3d_volume(self.tranches, x1, y1, x2, y2, z1, z2)
        self._dernier_args = args_actuels
        self._dernier_resultats = (p, a_cool, dh, ph, pbox, pwr)
        return p, a_cool, dh, ph, pbox, pwr
    
    def _get_cylinders(self, data):
        """
        Extrait la liste des cylindres (x, y, R) pour UNE tranche axiale spécifique.
        """
        cylinders = []
        pin_geo = data['PIN_GEOMETRY']
        ass_geo = data['ASSEMBLY_GEOMETRY']
        wr_geo = data['WATER_ROD_GEOMETRY']
        lattice = ass_geo['lattice_description']
        exclusions = set(ass_geo['non_fuel_rod_ids'])

        d = pin_geo['pin_pitch']
        r_clad = pin_geo['clad_radius']
        r_wr = wr_geo['outer_radius']
        
        W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
        L_ext = ass_geo['assembly_pitch']
        L_int = L_ext - 2 * W_start
        l_gap_int = (L_int - ((len(lattice[0]) - 1) * d) - 2 * r_clad) / 2.0

        # 1. Crayons Combustibles (Fuel rods)
        for row_idx, row in enumerate(lattice):
            for col_idx, item in enumerate(row):
                if item not in exclusions:
                    cx = W_start + l_gap_int + r_clad + col_idx * d
                    cy = W_start + l_gap_int + r_clad + row_idx * d
                    cylinders.append((cx, cy, r_clad))

        # 2. Tubes d'eau (Water rods)
        for center in wr_geo['centers']:
            cylinders.append((center[0], center[1], r_wr))

        return cylinders
    
    def _get_active_y_bounds(self, x_slice):
        """
        Calculate the active Y bounds (y values within the inner box) for a given X slice in the YZ plane.

        Parameters:
        - x_slice: the x-coordinate of the slice in the YZ plane

        Returns:
        - (y_min, y_max): the minimum and maximum y-coordinates of the active region for this x_slice
        """
        ass_geo = self.data_ref['ASSEMBLY_GEOMETRY']
        W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
        W_end = ass_geo['assembly_pitch'] - W_start
        R_c = ass_geo['corner_inner_radius_of_curvature']
        
        # 1. Complètement en dehors du boîtier
        if x_slice <= W_start or x_slice >= W_end:
            return None, None
            
        # 2. Zone centrale (bords droits)
        if W_start + R_c <= x_slice <= W_end - R_c:
            return W_start, W_end
            
        # 3. Dans l'arrondi gauche
        if x_slice < W_start + R_c:
            dx = (W_start + R_c) - x_slice # Distance au centre du coin
            dy = math.sqrt(max(0, R_c**2 - dx**2)) 
            return W_start + R_c - dy, W_end - R_c + dy
            
        # 4. Dans l'arrondi droit
        if x_slice > W_end - R_c:
            dx = x_slice - (W_end - R_c) # Distance au centre du coin
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
        ass_geo = self.data_ref['ASSEMBLY_GEOMETRY']
        W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
        W_end = ass_geo['assembly_pitch'] - W_start
        R_c = ass_geo['corner_inner_radius_of_curvature']
        
        # 1. Complètement en dehors du boîtier
        if y_slice <= W_start or y_slice >= W_end:
            return None, None
            
        # 2. Zone centrale (bords droits)
        if W_start + R_c <= y_slice <= W_end - R_c:
            return W_start, W_end
            
        # 3. Dans l'arrondi bas
        if y_slice < W_start + R_c:
            dy = (W_start + R_c) - y_slice # Distance au centre du coin en Y
            dx = math.sqrt(max(0, R_c**2 - dy**2)) 
            return W_start + R_c - dx, W_end - R_c + dx
            
        # 4. Dans l'arrondi haut
        if y_slice > W_end - R_c:
            dy = y_slice - (W_end - R_c) # Distance au centre du coin en Y
            dx = math.sqrt(max(0, R_c**2 - dy**2)) 
            return W_start + R_c - dx, W_end - R_c + dx
            
        return None, None
    
    def get_z_global_bounds(self):
        """ 
        Get the global minimum and maximum axial coordinates (z) across all slices.
        
        Parameters: None
        
        Returns:
        - z_min: Global minimum axial coordinate
        - z_max: Global maximum axial coordinate
        """
        z_min = min(t['z_start'] for t in self.tranches)
        z_max = max(t['z_end'] for t in self.tranches)
        return z_min, z_max
    
    def get_x_global_bounds(self):
        """ 
        Get the global minimum and maximum x coordinates of the geometry.
        
        Parameters: None
        
        Returns:
        - x_min: Global minimum x coordinate
        - x_max: Global maximum x coordinate
        """
        x_min = 0.0
        x_max = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
        return x_min, x_max 
     
    def get_y_global_bounds(self):
        """ 
        Get the global minimum and maximum y coordinates of the geometry.
        
        Parameters: None
        
        Returns:
        - y_min: Global minimum y coordinate
        - y_max: Global maximum y coordinate
        """
        y_min = 0.0
        y_max = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
        return y_min, y_max

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

# Canal (Water Channel)
    def get_porosity_z_canal(self, i, j, z1, z2):
        """Public method to get porosity in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - porosity: Porosity in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
        return self.get_porosity_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_a_cool_z_canal(self, i, j, z1, z2):
        """Public method to get total fluid area in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - a_cool: Total area of fluid in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
        return self.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_dh_z_canal(self, i, j, z1, z2):
        """Public method to get hydraulic diameter in the z-direction for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - hydraulic_diameter: Hydraulic diameter in the z-direction of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
        return self.get_dh_z_cv(x1, y1, x2, y2, z1, z2)
    
    def get_ph_canal(self, i, j, z1, z2):
        """Public method to get average heating perimeter for a water channel 
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - average_heating_perimeter: Average heating perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
        return self.get_ph_cv(x1, y1, x2, y2, z1, z2)

    def get_pbox_canal(self, i, j, z1, z2):
        """
        Public method to get average inner box perimeter for a water channel
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - average_inner_box_perimeter: Average inner box perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
        return self.get_pbox_cv(x1, y1, x2, y2, z1, z2)
    
    def get_pwr_canal(self, i, j, z1, z2):
        """
        Public method to get average water rod perimeter for a water channel
        defined by its lattice position (i, j) and axial bounds (z1, z2).
        
        Parameters:
        - i: Row index of the canal (1 is bottom row)
        - j: Column index of the canal (1 is left column)
        - z1: Starting axial coordinate of the canal
        - z2: Ending axial coordinate of the canal
        
        Returns:
        - average_water_rod_perimeter: Average water rod perimeter of the specified water channel
        """
        x1, y1, x2, y2 = self._obtenir_bornes_canal(i, j)
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

    def execute_profile_z(self, section_type, h, p, z_min, z_max):
        """
        Execute an axial profile analysis for the specified control volume, rod, or canal.
        It computes porosity, fluid area, hydraulic diameter,
        and heating perimeter in the z direction at multiple axial positions.

        Parameters:
        - section_type: Tuple indicating the type of profile and its parameters
            - For control volume: ('cv', (x1, y1, x2, y2))
            - For rod: ('rod', (i, j))
            - For water channel: ('water', (i, j))
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
        - cible_str: String identifier of the target (e.g. "cv_x1_y1_x2_y2", "rod_i_j", "water_i_j") for labeling purposes
        """
        z_coords = []
        porosities, a_cools, dhs, phs = [], [], [], []
        
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
                phi = self.get_porosity_z_canal(i, j, z1, z2)
                a_cool = self.get_a_cool_z_canal(i, j, z1, z2)
                dh = self.get_dh_z_canal(i, j, z1, z2)
                ph = self.get_ph_canal(i, j, z1, z2)
                cible_str = f"water_{i}_{j}"
                
            z_coords.append(curr_z + h/2.0) 
            porosities.append(phi)
            a_cools.append(a_cool)
            dhs.append(dh)
            phs.append(ph)
            curr_z += p
        return z_coords, porosities, a_cools, dhs, phs, cible_str
    
    def execute_mesh_z(self, mesh_type, z1, z2):
        """
        Execute a mesh analysis for the entire assembly based on the 
        specified mesh type (fuel, water, or regular).
        It computes porosity, fluid area, hydraulic diameter, and heating perimeter in the z direction
        for each cell in the mesh and organizes the results in matrices of x and y coordinates.
        
        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For fuel: ('fuel', None)
            - For water: ('water', None)
            - For regular: ('regular', (n))
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        
        Returns:
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        - mat_ph: 2D list (matrix) of heating perimeter values for each cell
        
        The structure of the matrices depends on the mesh type:
        - For 'fuel': mat_p[i][j] corresponds to the rod at lattice position (i+1, j+1)
        - For 'water': mat_p[i][j] corresponds to the canal at lattice position (i+1, j+1)
        - For 'regular': mat_p[i][j] corresponds to the control volume defined by the regular grid cell
        """
        mat_p, mat_dh, mat_ph = [], [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])

        if mesh_type[0] == 'fuel':
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
                    row_p.append(self.get_porosity_z_canal(i, j, z1, z2))
                    row_dh.append(self.get_dh_z_canal(i, j, z1, z2))
                    row_ph.append(self.get_ph_canal(i, j, z1, z2))
                mat_p.append(row_p); mat_dh.append(row_dh); mat_ph.append(row_ph)
                
        elif mesh_type[0] == 'regular':
            L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
            W_start = self.data_ref['ASSEMBLY_GEOMETRY']['gap_wide'] + self.data_ref['ASSEMBLY_GEOMETRY']['channel_box_thickness']
            L_int = L_ext - 2 * W_start
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
    

#--- X-DIRECTION: PUBLIC METHODS FOR CROSS-FLOW POROSITY AND HYDRAULIC DIAMETER ---

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
        for tranche in self.tranches:
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
        
        if box_y_min is None:
            return 0.0

        # 2. On restreint notre fenêtre d'étude [y1, y2] à l'intérieur du boîtier
        min_y = max(y1, box_y_min)
        max_y = min(y2, box_y_max)
        
        w_tot_valide = max_y - min_y
        if w_tot_valide <= 1e-9:
            return 0.0

        # 3. Boucle d'intégration axiale
        for tranche in self.tranches:
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

# Rods
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


# Canals
    def get_porosity_x_canal(self, i, x_slice, z1, z2):
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
        _, y1, _, y2 = self._obtenir_bornes_canal(i, 1)
        return self.get_porosity_x_cv(x_slice, y1, y2, z1, z2)

    def get_dh_x_canal(self, i, x_slice, z1, z2):
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
        _, y1, _, y2 = self._obtenir_bornes_canal(i, 1)
        return self.get_dh_x_cv(x_slice, y1, y2, z1, z2)

# LATERAL X-PROFILE AND YZ MESH ---

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
            _, y1, _, y2 = self._obtenir_bornes_canal(i, 1)
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
    
    def execute_mesh_x(self, mesh_type, x_slice, z1, z2, n_z):
        """
        Execute a mesh analysis in the x-direction for the entire assembly based on the
        specified mesh type (fuel, water, or regular). It computes porosity and hydraulic diameter
        for each cell in the mesh and organizes the results in matrices of y and z coordinates.

        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For fuel: ('fuel', None)
            - For water: ('water', None)
            - For regular: ('regular', n_y) where n_y is the number of divisions in the y-direction
        - x_slice: The x-coordinate of the slice in the YZ plane where the properties are computed
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        - n_z: Number of divisions in the z-direction for the mesh

        Returns:         
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        """
        mat_p, mat_dh = [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_rows = len(lattice)
        
        # 1. Création des graduations pour Z (découpage régulier de la hauteur)
        z_bounds = [z1 + k * (z2 - z1) / n_z for k in range(n_z + 1)]
        
        # 2. Détermination des intervalles sur l'axe Y (les bandes horizontales)
        y_intervals = []
        
        if mesh_type[0] == 'regular':
            L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
            W_start = self.data_ref['ASSEMBLY_GEOMETRY']['gap_wide'] + self.data_ref['ASSEMBLY_GEOMETRY']['channel_box_thickness']
            L_int = L_ext - 2 * W_start
            
            # On découpe l'espace interne en mesh_type[1] parties égales
            n_y = mesh_type[1]
            y_bounds = [W_start + k * L_int / n_y for k in range(n_y + 1)]
            for j in range(n_y):
                y_intervals.append((y_bounds[j], y_bounds[j+1]))
                
        elif mesh_type[0] == 'fuel':
            # On parcourt les lignes de crayons de haut en bas (i de n_rows à 1)
            for i in range(n_rows, 0, -1):
                _, y1_rod, _, y2_rod = self._obtenir_bornes_rod(i, 1) # j=1 suffit pour avoir Y
                y_intervals.append((y1_rod, y2_rod))
                
        elif mesh_type[0] == 'water':
            # On parcourt les canaux d'eau de haut en bas (il y a n_rows + 1 canaux Y)
            for i in range(n_rows + 1, 0, -1):
                _, y1_canal, _, y2_canal = self._obtenir_bornes_canal(i, 1)
                y_intervals.append((y1_canal, y2_canal))
                
        else:
            raise ValueError(f"Type de maillage inconnu : {mesh_type[0]}")

        # 3. La double boucle de construction des matrices
        # Convention : on lit Z de haut en bas pour un affichage Heatmap (Matplotlib) logique
        for idx_z in range(n_z):
            z1_cell, z2_cell = z_bounds[idx_z], z_bounds[idx_z+1]
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

#--- Y-DIRECTION: PUBLIC METHODS FOR CROSS-FLOW POROSITY AND HYDRAULIC DIAMETER ---
    
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

        for tranche in self.tranches:
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

        for tranche in self.tranches:
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
    
# Rods
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
# Canals
    def get_porosity_y_canal(self, j, y_slice, z1, z2):
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
        x1, _, x2, _ = self._obtenir_bornes_canal(1, j)
        return self.get_porosity_y_cv(y_slice, x1, x2, z1, z2)

    def get_dh_y_canal(self, j, y_slice, z1, z2):
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
        x1, _, x2, _ = self._obtenir_bornes_canal(1, j)
        return self.get_dh_y_cv(y_slice, x1, x2, z1, z2)

#--- LATERAL Y-PROFILE AND XZ MESH ---

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
            x1, _, x2, _ = self._obtenir_bornes_canal(1, j)
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

    def execute_mesh_y(self, mesh_type, y_slice, z1, z2, n_z):
        """
        Execute a 2D mesh on the XZ plane (Cross-Flow in Y) at a fixed y_slice coordinate.
        The mesh is constructed based on the specified mesh type (fuel, water, or regular) and the assembly geometry.
        It computes porosity and hydraulic diameter for each cell in the mesh 
        and organizes the results in matrices of x and z coordinates.

        Parameters:
        - mesh_type: Tuple indicating the type of mesh and its parameters
            - For fuel: ('fuel', None)
            - For water: ('water', None)
            - For regular: ('regular', n_x) where n_x is the number of divisions in the x-direction
        - y_slice: The y-coordinate of the slice in the XZ plane where the properties are computed
        - z1: Starting axial coordinate of the mesh
        - z2: Ending axial coordinate of the mesh
        - n_z: Number of divisions in the z-direction for the mesh

        Returns:
        - mat_p: 2D list (matrix) of porosity values for each cell
        - mat_dh: 2D list (matrix) of hydraulic diameter values for each cell
        """
        mat_p, mat_dh = [], []
        lattice = self.data_ref['ASSEMBLY_GEOMETRY']['lattice_description']
        n_cols = len(lattice[0])
        
        z_bounds = [z1 + k * (z2 - z1) / n_z for k in range(n_z + 1)]
        x_intervals = []
        
        if mesh_type[0] == 'regular':
            L_ext = self.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch']
            W_start = self.data_ref['ASSEMBLY_GEOMETRY']['gap_wide'] + self.data_ref['ASSEMBLY_GEOMETRY']['channel_box_thickness']
            L_int = L_ext - 2 * W_start
            
            n_x = mesh_type[1]
            x_bounds = [W_start + k * L_int / n_x for k in range(n_x + 1)]
            for j in range(n_x):
                x_intervals.append((x_bounds[j], x_bounds[j+1]))
                
        elif mesh_type[0] == 'fuel':
            for j in range(1, n_cols + 1):
                x1_rod, _, x2_rod, _ = self._obtenir_bornes_rod(1, j)
                x_intervals.append((x1_rod, x2_rod))
                
        elif mesh_type[0] == 'water':
            for j in range(1, n_cols + 2):
                x1_canal, _, x2_canal, _ = self._obtenir_bornes_canal(1, j)
                x_intervals.append((x1_canal, x2_canal))

        for idx_z in range(n_z):
            z1_cell, z2_cell = z_bounds[idx_z], z_bounds[idx_z+1]
            row_p, row_dh = [], []
            
            for (x1_cell, x2_cell) in x_intervals:
                phi = self.get_porosity_y_cv(y_slice, x1_cell, x2_cell, z1_cell, z2_cell)
                dh = self.get_dh_y_cv(y_slice, x1_cell, x2_cell, z1_cell, z2_cell)
                
                row_p.append(phi)
                row_dh.append(dh)
                
            mat_p.append(row_p)
            mat_dh.append(row_dh)
            
        return mat_p, mat_dh