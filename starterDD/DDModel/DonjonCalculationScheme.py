### Initial implementation of a Donjon Calculation scheme.

# Date : 18/08/2026
# Author : R. Guasch
# Purpose : provide a structure to hold steps for a Donjon calculation, 

# Minimal steps :
# - Geometry creation and material indexation
# - Fuel Map creation
# - Finite element analysis for Diffusion or SPn or Sn tracking with SNT: for transport.

import numpy as np

#_DONJON_CALCULATION_STEPS = [""]

class DonjonCalculationScheme:
    """
    Minimal implementation : 
    Choice of solver / operator : Diffusion, SPn or Sn.
    Choice of mesh refinement : number of axial meshes, radial spatial homogenisation scheme

    Steps : 
      - ModelInitialisation : Geometry, Matex, Fmap creation
      - Traking : finite element analysis of the geometry
      - THsolve : solver thermal-hydraulics problem (THM: call)
      - NeutronicsSolve : Interpolate cross sections and solve for neutron flux -> power density.
    """

    def __init__(self, name="donjon_scheme"):

        self.name = name
        self.steps = []

    def add_initialisation_step(self, step):
        if not isinstance(step, (ModelInitialisation)):
            raise TypeError(
                f"Expected ModelInitialisation, got {type(step).__name__}"
            )
        self.initialisation_step = step
        self.steps.append(step)

    def add_neutronics_step(self, step):
        if not isinstance(step, (NeutronicsSolve)):
            raise TypeError(
                f"Expected NeutronicsSolve, got {type(step).__name__}"
            )
        self.neutronics_step = step
        self.steps.append(step)

    def add_step(self, step):
        """
        Append a ``ModelInitialisation``, ``ThermalHydraulicsSolve`` or ``NeutronicsSolve`` to the scheme.

        Parameters
        ----------
        step : ModelInitialisation, ThermalHydraulicsSolve or NeutronicsSolve
        """
        if not isinstance(step, (ModelInitialisation, ThermalHydraulicsSolve, NeutronicsSolve)):
            raise TypeError(
                f"Expected ModelInitialisation, ThermalHydraulicsSolve or NeutronicsSolve "
                f"got {type(step).__name__}"
            )
        self.steps.append(step)

    def get_step(self, step_type):
        """
        Retrieve a step by step_type.

        Parameters
        ----------
        step_type : str

        Returns
        -------
        step with given step_type

        Raises
        ------
        KeyError
            If no step with the given step_type exists.
        """
        for step in self.steps:
            if isinstance(step, step_type):
                return step
        raise KeyError(f"No step with type '{step_type}' in scheme '{self.name}'.")

class ModelInitialisation:
    """
    Instantiate a Donjon model by defining geometry, meshing, material indexation, FuelMap parameters, 
    """

    def __init__(self, core_model, radial_homogenization_strategy , boundary_conditions_dict, number_of_axial_materials_per_slice, number_of_energy_groups, interpolation_variables):
        """
        core_model : DDModel::DonjonModel::CoreModel object holding the geometry definition information
        radial_homogenization_strategy  (str) : "by_pin" or "by_assembly" to define radial, sub-assembly meshing strategy.
        boundary_conditions_dict (dictionary): associate a x/y/z direction to a type of boundary condition.
            ex : {x: "reflective", y: "reflective", z: "void"}
        number_of_axial_materials_per_slice list of (int) : number of material subdvisions in each slice : should have the same length as number of slices.
        interpolation_variables : list of tuples representing the variable parameter name used in the compo, whether its is local or global, and the variable type.
            Variable type is recommended to be one of the following types :
            "boron_concentration", "effective_fuel_temperature", "fuel_surface_temperature", "coolant_temperature", "coolant_density", "coolant_pressure", "moderator_temperature", "moderator_density"

        idea: add support for SPLITX/Y/Z to include more calculation nodes than mixes?     
        """

        self.core_model = core_model
        self.radial_homogenization_strategy = radial_homogenization_strategy
        self.boundary_conditions = boundary_conditions_dict
        self.number_of_axial_materials_per_slice = number_of_axial_materials_per_slice
        self.number_of_energy_groups = number_of_energy_groups
        self.slice_to_compodir_correspondance = {}
        self.fuel_mixtures = []
        self.validate_interpolation_variables(interpolation_variables)


        self.core_model.createAssemblyModels()
        self.resolve_2D_fuel_map()

    def validate_interpolation_variables(self, interpolation_variables):
        AVAILABLE_PARAMETERS = {"fuel_temperature": "T-FUEL", 
                            "coolant_temperature": "T-COOL", 
                            "coolant_density": "D-COOL", 
                            "boron_concentration": "C-BORE", 
                            "fuel_surface_temperature": "T-SURF", 
                            "coolant_pressure": "P-COOL", 
                            "moderatore_temperature": "T-MODE", 
                            "moderator_density": "D-MODE"
                        }
        var_name_flag_pname = {}
        for var in interpolation_variables:
            print(var)
            print(var[0])
            print(var[1])
            if var[1] not in ["global", "local"]:
                raise ValueError(f"Only local or global parameters are supported, got {var[1]}.")
            if var[2] not in AVAILABLE_PARAMETERS.keys():
                raise ValueError(f"Parameter names in {AVAILABLE_PARAMETERS.keys()} are supported, got {var[2]}.")
            var_name_flag_pname[var[0]] = (var[1].upper(), AVAILABLE_PARAMETERS[var[2]])    



        self.interpolation_variables = var_name_flag_pname

    def resolve_mesh(self):
        meshx = [0.0]
        meshy = [0.0]
        meshz = [0.0]
        ny = len(self.core_model.core_2D_layout)
        nx = len(self.core_model.core_2D_layout[0])
        for idx in range(nx):
            # set y = 0, traverse x direction to retrieve xmesh
            idy = 0
            assembly_id = self.core_model.core_2D_layout[idy][idx]
            extruded_assembly = self.core_model.assemblies[(idx, idy, assembly_id)]
            first_slice = extruded_assembly.slices_2D[0]
            if assembly_id == "reflector":
                reflector_model = self.core_model.reflector_models[((idx, idy, assembly_id), first_slice)]
                dimensions = reflector_model.dimensions
                meshx.append(meshx[-1] + dimensions[1])
            else: 
                dragon_model_ij0 = self.core_model.assembly_models[((idx, idy, assembly_id), first_slice)]
                assembly_pitch = dragon_model_ij0.assembly_pitch

                if self.radial_homogenization_strategy == "by_assembly":
                    meshx.append(meshx[-1] + assembly_pitch)
                elif self.radial_homogenization_strategy == "by_pin":
                    nx_lat = len(dragon_model_ij0.lattice[0])
                    lattice_meshx = np.zeros(nx_lat+2)
                    ini_delx = dragon_model_ij0.translation_offset_x

                    # sweep x direction to reconstruct pin-by-pin lattice_meshx
                    lattice_meshx[0] = ini_delx
                    for i in range(1, nx_lat+2):
                        lattice_meshx[i] = ini_delx + i * dragon_model_ij0.pin_geometry_dict["pin_pitch"]
                    lattice_meshx[-1] = dragon_model_ij0.assembly_pitch
                    last_meshx = meshx[-1]
                    for mesh_pt in lattice_meshx:
                        meshx.append(last_meshx + mesh_pt)
                else:
                    raise ValueError(f"Invalid spatial homogenization strategy : {self.radial_homogenization_strategy}")
                
        for idy in range(ny):
            # set x = 0, traverse y direction to retrieve xmesh
            idx = 0
            assembly_id = self.core_model.core_2D_layout[idy][idx]
            extruded_assembly = self.core_model.assemblies[(idx, idy, assembly_id)]
            first_slice = extruded_assembly.slices_2D[0]
            if assembly_id == "reflector":
                reflector_model = self.core_model.reflector_models[((idx, idy, assembly_id), first_slice)]
                dimensions = reflector_model.dimensions
                meshy.append(meshy[-1] + dimensions[1])
            else: 
                dragon_model_ij0 = self.core_model.assembly_models[((idx, idy, assembly_id), first_slice)]
                assembly_pitch = dragon_model_ij0.assembly_pitch

                if self.radial_homogenization_strategy == "by_assembly":
                    meshy.append(meshy[-1] + assembly_pitch)
                elif self.radial_homogenization_strategy == "by_pin":
                    ny_lat = len(dragon_model_ij0.lattice)
                    lattice_meshy = np.zeros(ny_lat+2)
                    ini_dely = dragon_model_ij0.translation_offset_y

                    # sweep x direction to reconstruct pin-by-pin lattice_meshx
                    lattice_meshy[0] = ini_dely
                    for i in range(1, ny_lat+2):
                        lattice_meshy[i] = ini_delx + i * dragon_model_ij0.pin_geometry_dict["pin_pitch"]
                    lattice_meshy[-1] = dragon_model_ij0.assembly_pitch
                    last_meshy = meshy[-1]
                    for mesh_pt in lattice_meshy:
                        meshy.append(last_meshy + mesh_pt)
                else:
                    raise ValueError(f"Invalid spatial homogenization strategy : {self.radial_homogenization_strategy}")

        # Analyse zmesh, assume central assembly is not a reflector and all assemblies on x-y plane share the same axial bounds
        # TODO : generalise to finding assemblies and merging meshes

        idx, idy = int(nx/2), int(ny/2)
        assembly_id = self.core_model.core_2D_layout[idy][idx]
        extruded_assembly = self.core_model.assemblies[(idx, idy, assembly_id)]
        nz_bounds = len(extruded_assembly.z_bounds)
        n_slices = nz_bounds - 1
        total_axial_meshes = 0 
        for nslice in range(n_slices):
            last_meshz = meshz[-1]
            slice_thickness =  extruded_assembly.z_bounds[nslice+1] -  extruded_assembly.z_bounds[nslice] 
            number_of_axial_mesh_points_in_slice = self.number_of_axial_materials_per_slice[nslice]
            total_axial_meshes += number_of_axial_mesh_points_in_slice
            delta_z_in_slice = slice_thickness / number_of_axial_mesh_points_in_slice
            for idz in range(number_of_axial_mesh_points_in_slice):
                meshz.append(last_meshz + (idz + 1) * delta_z_in_slice) # issue here with meshz construction : should end up with nz+1 bounds
        

        self.meshx = meshx
        self.meshy = meshy
        self.meshz = meshz

        self._build_fuel_index_mapping()

    def resolve_2D_fuel_map(self):
        """
        The goal is to analyse the core layout and identify which elements are reflectors.
        """

        reflector_pairs = []
        fuel_channel_pairs = []
        ny = len(self.core_model.core_2D_layout)
        nx = len(self.core_model.core_2D_layout[0])
        for j in range(ny):
            for i in range(nx):
                assembly_identifier = self.core_model.core_2D_layout[j][i]
                if assembly_identifier == "reflector":
                    reflector_pairs.append((i,j))
                else:
                    fuel_channel_pairs.append((i,j))

        self.fuel_channel_pairs = fuel_channel_pairs
        self.reflector_pairs = reflector_pairs
        self.number_fuel_channels = len(self.fuel_channel_pairs)

    def _build_fuel_index_mapping(self):
        ny = len(self.core_model.core_2D_layout)
        nx = len(self.core_model.core_2D_layout[0])
        nz = len(self.meshz) -1
        mapping = {}
        m = 1

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if (i, j) in self.reflector_pairs:
                        mapping[(i, j, k)] = 0
                    else:
                        mapping[(i, j, k)] = m
                        m += 1
        self.fuel_index_mapping = mapping


    def _get_fuel_index_from_ijk(self, i, j, k):
        """
        return fuel index from position (i,j,k),
        if (i,j) is not a reflector, return the unique fuel index,
        if (i,j) is a refelctor, return 0 
        """
        return self.fuel_index_mapping[(i,j,k)]


    def set_slice_to_compodir_correspondance(self, xpos, ypos, z_bounds, compo_name, dir_name):
        """
        Helper to set COMPONAM and DIRNAM in corresponding to a slice in a 3D model.
        Provide flexibility for users to specify what compo LCM object is linked to what set of mesh points, belonging to a slice.
        Provide flexibility for users to specify which directory is to be used in the specified LCM object.

        Slice is identified by :
        xpos, ypos (int) : x/y position of the slice in the 2D core layout according to x-increasing, y-increasing numbering
        z_bounds (list of float) : list of 2 axial coordinates specifying the axial bounds.  
        and associated to :
        compo_name (str) : name of the LCM compo object to be associated with a slice
        dir_name (str) : name of the directory to read from the compo
        """

        self.slice_to_compodir_correspondance[(xpos, ypos, z_bounds[0], z_bounds[-1])] = (compo_name, dir_name)

    def _get_slice_info_from_k(self, i, j, k):
        """
        Helper to get slice physical bounds from an index k in the axial mesh
        i, j (int) : assembly id indices on the 2D core layout
        k (int) : axial plane index
        """
        k_plane_min_bound = self.meshz[k]
        k_plane_max_bound = self.meshz[k+1]
        assembly_id = self.core_model.core_2D_layout[j][i]
        extruded_assembly = self.core_model.assemblies[(i, j, assembly_id)]
        z_bounds = extruded_assembly.z_bounds
        n_slices = len(z_bounds) - 1

        for slice_number in range(n_slices):
            z_min_slice = z_bounds[slice_number]
            z_max_slice = z_bounds[slice_number+1]
            if k_plane_min_bound >= z_min_slice and k_plane_max_bound <= z_max_slice:
                slice_info = z_min_slice, z_max_slice
            else: 
                slice_info = None
                raise ValueError(f"No intersection between meshz plane k={k} and slices data were found for positions (i,j) = ({i},{j})")
        return slice_info

    def _get_ks_from_slice_info(self, slice_z_bounds):
        """
        Return the plane index k from x-y indexing and slice axial bounds info
        Limited to the assumption that all fuel channels have the same axial geometry / bounds
        """

        k_indices = []
        for k in range(len(self.meshz)-1):
            if self.meshz[k] >= min(slice_z_bounds) and self.meshz[k+1] <= max(slice_z_bounds):
                # the kth zmesh plane was found to intersect with the given slice data
                k_indices.append(k)
        if len(k_indices) == 0:
            print(f"Warning : no intersection found between axial mesh {self.meshz} and the specified slice bounds {slice_z_bounds}.")
        return k_indices

    
    def _bc_keyword_handling(self, direction):
        """
        direction : "x", "y" or "z"

        returns the GEO: keyword to treat the type of user facing boundary condition specified
        """

        if self.boundary_conditions[direction] == "reflective":
            return "REFL"
        elif self.boundary_conditions[direction] == "void":
            return "VOID"
        else:
            raise ValueError(f"The type of boundary condition {self.boundary_conditions[direction]} is not implemented yet!")


    def _set_number_of_fuel_mixtures(self, nfuel):
        """
        Set the total number of fuel mixtures
        """
        self.nfuel = nfuel

    def _add_fuel_mix_index(self, index):
        """
        Add a mixture index to the set of indices used for fuel mixtures.
        """
        self.fuel_mixtures.append(index)

    def set_reactor_power(self, power):
        """
        set the total thermal power which should be used in the full core problem
        
        power (float) : total thermal power of the reaxctor in W
        """
        self.reactor_power = power

    def set_initial_axial_power_form(self, power_shape):
        """
        Initialize the axial power profile based on power_shape.
        power_shape :: (str) | (list) | (np.ndarray) : 'cosine', 'sine', 'uniform' or a list/array of length nz representing the power distribution along the axial dimension of the fuel assembly.
        """
        profile = []
        nz_meshes = len(self.meshz) - 1
        for i in range(nz_meshes):
            z_norm = (i + 0.5) / nz_meshes
            if isinstance(power_shape, (list, np.ndarray)) and len(power_shape) == nz_meshes:
                val = power_shape[i]
            elif power_shape == 'cosine':
                val = np.cos((np.pi / 2.0) * z_norm)
            elif power_shape == 'sine':
                val = np.sin(np.pi * z_norm)
            elif power_shape == 'uniform':
                val = 1.0
            else: 
                raise ValueError("Invalid power_distribution input. Must be 'cosine', 'sine', 'uniform', or a list/array of length nz.")
            profile.append(val)
        
        self.axial_power_form = np.array(profile / np.mean(profile)) 

    def set_initial_parameters(self, parameter_to_value_dict):
        """
        When initializing a model, set initial values for desired parameters :
        parameter_to_value_dict (dictionnary) : key variable name used as parameter in the compo / FuelMap, value : initial value guess for that variable
        """

        self.initial_parameters = parameter_to_value_dict

    def get_initial_parameter(self, parameter_name):
        """
        Return the value of an initial parameter
        """
        if parameter_name in self.initial_parameters.keys():
            return self.initial_parameters[parameter_name]
        else:
            raise ValueError(f"Parameter with name {parameter_name} not available in the set of initial parameters {self.initial_parameters}") 

    def set_fuel_mass(self, total_fuel_mass):
        """
        set the total mass of fuel present in the problem in kg.
        """
        self.total_fuel_mass = total_fuel_mass




class ThermalHydraulicsSolve:
    """
    placeholder for now
    """

    def __init__(self, name="THSolve"):
        self.name = name


class NeutronicsSolve:
    """
    Class representing a Neutronics Solve step in a Donjon Calculation Scheme
    Holds the necessary information to compute the neutron flux in Donjon based on :
    initial_model (ModelInitialisation) :  object carying the information related to the calculation scheme's initialisation step,
    operator (str) : diffusion / transport operator, option to select the type of neutronics solution
    interpolation_type (str) : type of interpolation to be used when evaluating the neutron cross sections 
    """

    def __init__(self, initial_model, operator, interpolation_type, discretization_type):
        self.model = initial_model
        self.operator = operator
        self.interpolation_type = interpolation_type
        self.parameter_values_dict = {}
        self.read_from_lcm = {}
        self.discretization_type = discretization_type

    def set_interpolation_parameter(self, parameter_key, values, from_lcm=False):
        """
        Set the values a parameter should take before interpolation of cross sections
        parameter_key (str) : key to the parameter to set interpolation values from,
        values : list of values to set the corresponding parameter to,
        from_lcm (boolean) : default set to False -> set parameters to values stored in the list of values, if set to True : 
            define a list of variables to recover data from a seperate   
        """
        nbChannels = self.model.number_fuel_channels
        nz = len(self.model.meshz) - 1
        if parameter_key not in self.model.interpolation_variables.keys():
            raise ValueError(f"Invalid parameter key : got {parameter_key}, model was initialised with the following parameters : {self.model.interpolation_variables.keys()}.")
        if len(values) != nbChannels * nz and self.model.interpolation_variables[parameter_key][0] == "LOCAL":
            raise ValueError (f"For local parameter {parameter_key} : Expected {nbChannels*nz} local parameters, got {len(values)}.")

        self.parameter_values_dict[parameter_key] = values
        self.read_from_lcm[parameter_key] = from_lcm

    def resolve_cpo_dir_to_mix(self):
        """
        Resolve the cpo/compodir to i,j,k association to find interpolation points and define mixes
        Limited to non-reflector assembly ids
        """

        print(self.model.slice_to_compodir_correspondance)
        self.compo_dir_to_mix = {}
        nx = len(self.model.meshx) - 1
        ny = len(self.model.meshy) - 1
        for slice_data in self.model.slice_to_compodir_correspondance.keys():
            slice_z_min = slice_data[2]
            slice_z_max = slice_data[3]
            k_indices = self.model._get_ks_from_slice_info([slice_z_min, slice_z_max])
            if self.model.slice_to_compodir_correspondance[slice_data] not in self.compo_dir_to_mix.keys():
                self.compo_dir_to_mix[self.model.slice_to_compodir_correspondance[slice_data]] = []
            for idz in k_indices:
                for idy in range(ny):
                    for idx in range(nx):
                        if self.model.fuel_index_mapping[(idx, idy, idz)] != 0:
                           self.compo_dir_to_mix[self.model.slice_to_compodir_correspondance[slice_data]].append(self.model.fuel_index_mapping[(idx, idy, idz)])  


        

        
        





    
                

