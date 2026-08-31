## collection of classes to handle DONJON module calls
# Authors: B. Godard, R. Guasch
# Date: 16/04/2026 (creation)
# Updated : 18/08/2026
# Purpose : Define class structures to handle DONJON module calls for starterDD package
# ----------------------------------------------------------------------------------------------- 

import os
import numpy as np
from .CLE2000 import (
    main_procedure, sub_procedure,
    validate_varname, wrap_cle2000_line,
)


def format_cle2000_float_array(vector, max_per_line=10):
    """
    List formating of a python list into a CLE2000 block
    vector : list to format,
    max_per_line  (optional / defaults to 10) : maximum values per line.
    """
    lines = []
    current_line = []
    for val in vector:
        current_line.append(f"{val:.5E}")
        if len(current_line) == max_per_line:
            lines.append(" ".join(current_line))
            current_line = []
    if current_line:
        lines.append(" ".join(current_line))
    return "\n".join(lines)

def format_cle2000_integer_array(array, max_per_line = 10):
    """
    List formating of a python list into a CLE2000 block
    vector : list to format,
    max_per_line  (optional / defaults to 10) : maximum values per line.
    """
    lines = []
    current_line = []
    for val in array:
        current_line.append(f"{int(val)}")
        if len(current_line) == max_per_line:
            lines.append(" ".join(current_line))
            current_line = []
    if current_line:
        lines.append(" ".join(current_line))
    return "\n".join(lines)


# -----------------------------------------------------------------------------------------------
# Classes representing DONJON modules
# -----------------------------------------------------------------------------------------------

class GEO:
    """
        Generate calls to the GEO: module, defining the reactor geometry.
        model_initialisation : ModelInitialisation object storing meshing information.

    """
    def __init__(self, lcm_geometry_name, model_initialisation):

        self.geometry_lcm_name = lcm_geometry_name
        self.nx = len(model_initialisation.meshx) - 1
        self.ny = len(model_initialisation.meshy) - 1
        self.nz = len(model_initialisation.meshz) - 1
        self.radial_homogenisation_option = model_initialisation.radial_homogenization_strategy
        self.x_boundary_condition = model_initialisation._bc_keyword_handling("x")
        self.y_boundary_condition = model_initialisation._bc_keyword_handling("y")
        self.z_boundary_condition = model_initialisation._bc_keyword_handling("z")
        self.initModel = model_initialisation


    def write_c2m(self, embedded_call=False):
        """
        embedded_call (boolean, optional) : option to exclude "GEONAM :=" section of the call
        to support embedded GEO: calls 
        """
        if embedded_call:
            call_option = f"  :::  GEO:"
        else:
            call_option = f"{self.geometry_lcm_name} := GEO:"
        block = (
            f"{call_option} :: CAR3D {self.nx} {self.ny} {self.nz}\n"
            f" X- {self.x_boundary_condition} X+ {self.x_boundary_condition} \n"
            f" Y- {self.y_boundary_condition} Y+ {self.y_boundary_condition} \n"
            f" Z- {self.z_boundary_condition} Z+ {self.z_boundary_condition} \n"
        )
        block += self.format_mesh()

        block += " MIX\n"
        for k in range(1, self.nz + 1):
            block += f"  PLANE {k}  {self.format_plane(k)}\n"
        block += " ;\n"
        if embedded_call is False:
            self.initModel._set_number_of_fuel_mixtures(len(self.initModel.fuel_mixtures))
        return block

    def format_plane(self, k):
        """
        Format definition of PLANE keyword
        k (int) : axial plane index
        """

        lines = []
        if self.radial_homogenisation_option == "by_assembly":
            # each assembly in a core 2D slice gets a unique number 
            # escape keyword "reflector" sets to 0 ?
            for j in range(self.ny):
                lines.append("\n")
                for i in range(self.nx):
                    unique_index = self.initModel._get_fuel_index_from_ijk(i,j,k-1)
                    lines.append(f"{unique_index}")
                    if unique_index != 0:
                        self.initModel._add_fuel_mix_index(unique_index)
        else: 
            raise ValueError(f"This type of radial homogenisation option is not implemented yet. Expected 'by_assembly', got {self.radial_homogenisation_option}")

        return " ".join(lines)


    def format_mesh(self, values_per_line=8):
        """
        Format a list of mesh bounds into CLE-2000 lines to be used within the MESHZ GEO: keyword.
        """
        xbounds = self.initModel.meshx
        ybounds = self.initModel.meshy
        zbounds = self.initModel.meshz

        tokens = [f"{x:.4f}" for x in xbounds]
        chunks = [tokens[i:i + values_per_line] for i in range(0, len(tokens), values_per_line)]

        firstx = " MESHX " + " ".join(chunks[0])
        restx = [" " * 7 + " ".join(chunk) for chunk in chunks[1:]]

        tokens = [f"{y:.4f}" for y in ybounds]
        chunks = [tokens[i:i + values_per_line] for i in range(0, len(tokens), values_per_line)]

        firsty = "\n MESHY " + " ".join(chunks[0])
        resty = [" " * 7 + " ".join(chunk) for chunk in chunks[1:]]

        tokens = [f"{z:.4f}" for z in zbounds]
        chunks = [tokens[i:i + values_per_line] for i in range(0, len(tokens), values_per_line)]

        firstz = "\n MESHZ " + " ".join(chunks[0])
        restz = [" " * 7 + " ".join(chunk) for chunk in chunks[1:]]
    

        lines = ( 
            "".join([firstx] + restx) + "\n"
            "".join([firsty] + resty) + "\n"
            "".join([firstz] + restz) + "\n"
        )
        return lines

class USPLIT:
    """
        Class handling calls to the USPLIT Donjon module generating a material indexation.
        geometry_name (str) : name of the geometry CLE-2000 variable.
        matex_name (str) : name of the material indexation object to be created by USPLIT:
        model_initialisation : ModelInitialisation object to hold information about the calculation to be performed.
    """

    def __init__(self, geometry_name, matex_name, model_initialisation):
        self.geometry_name = geometry_name
        self.matex_name = matex_name
        self.initModel = model_initialisation
        self.max_number_of_regions = 100000

    def write_c2m(self):
        block = ""
        block += f"{self.geometry_name} {self.matex_name} := USPLIT: {self.geometry_name} ::\n"
        block += f"NGRP {self.initModel.number_of_energy_groups} MAXR {self.max_number_of_regions} \n"
        block += f"NFUEL {self.initModel.nfuel} FMIX \n"
        block += f"{format_cle2000_integer_array(self.initModel.fuel_mixtures)} \n"
        block += ";\n"

        return block



class RESINI:
    """
        Calls handling calls to module RESINI: to create a Fuel Map.
        Minimal implementation : assume FMap and calculation geometries are the same.

        fuel_map_name (str)
        matex_name (str)
        model_initialiser (ModelInitialisation object)
        model_updater (NeutronicsSolve object) : if provided, the Fuel Map is updated accoding to values associated with the desired neutronics solution. 
    """
    def __init__(self, fuel_map_name, matex_name, model_initialiser, model_updater=None):

        self.fuel_map_name = fuel_map_name
        self.matex_name = matex_name
        self.model_initialiser = model_initialiser
        self.model_updater = model_updater

    def write_c2m(self):
        if self.model_updater == None: 
            initialisation = True
        else:
            initialisation = False

        if initialisation:
            geo_inner = GEO(lcm_geometry_name="Geom", model_initialisation=self.model_initialiser)
            geo_def = geo_inner.write_c2m(embedded_call=True)

            power_mw = self.model_initialiser.reactor_power*1e-6
            pform_str = format_cle2000_float_array(self.model_initialiser.axial_power_form, max_per_line=5)
            nx_names, ny_names = self.format_nx_ny_names()
            b_zones = self.format_burnup_zones()
            parameters_definition = self.format_parameters_definition()
            set_parameters_bloc = self.format_set_parameters_bloc(initialisation)
            fuel_weight_bloc = self.format_fuel_weight_bloc()

            block = (
                "*--------------------------------------------------------\n"
                "* Fuel map definition\n"
                "*--------------------------------------------------------\n"
                f"{self.fuel_map_name} {self.matex_name} := RESINI: {self.matex_name} ::\n"
                f"{geo_def}"
                f"NXNAME {nx_names} \nNYNAME {ny_names}\n" 
                f"NCOMB {self.model_initialiser.nfuel}\n"
                "B-ZONE \n"
                f"{b_zones}\n"
                f"{parameters_definition}\n"
                "    BTYPE INST-BURN INST-BVAL CHAN 0.0\n"
                f" REACTOR-POW {max(power_mw, 1e-10):.5E}\n"
                "  AXIAL-PFORM\n"
                f"{pform_str}\n"
                f"{set_parameters_bloc}\n"
                f"{fuel_weight_bloc}"
                ";\n"
            )
            return block

        else:
            # Create a precodure updating the fuel map by setting local parameters

            set_parameters_bloc = self.format_set_parameters_bloc(initialisation)

            block = (
                "*--------------------------------------------------------\n"
                "* Update Fuel Map \n"
                "*--------------------------------------------------------\n"
                f"{self.fuel_map_name} := RESINI: {self.fuel_map_name} ::\n"
                " BTYPE INST-BURN INST-BVAL CHAN 0.0\n"
                f"{set_parameters_bloc}\n"
                ";\n"
            )
            return block


    def format_nx_ny_names(self):
        """
        Format the NXNAME and NYNAME keyword entries
        This assumes that each entry in the cartesian nx by ny grid is a fuel assembly
        TODO : implement escape parameter if a channel is declared as non-fuel, denoted by '-'
        """
        NX_NAMES = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                    '21', '22', '23', '24', '25', '26']
        NY_NAMES = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
                    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
                    'U', 'V', 'W', 'X', 'Y', 'Z']
        nx = len(self.model_initialiser.meshx) - 1
        ny = len(self.model_initialiser.meshy) - 1
        if nx > len(NX_NAMES) or ny > len(NY_NAMES):
            raise ValueError(f"RESINI: Maximal core dimensions supported are 28 by 28, got {nx} by {ny}.")
        nxname_str = ""
        nyname_str = ""
        x_pos_with_refl = []
        y_pos_with_refl = []
        for pair in self.model_initialiser.reflector_pairs:
            x_pos_with_refl.append(pair[0])
            y_pos_with_refl.append(pair[1])
        for idx in range(nx):
            if len(x_pos_with_refl)>0:
                if idx == max(x_pos_with_refl) or idx == min(x_pos_with_refl):
                    nxname_str += "'-' "
                else:    
                    nxname_str += f"'{NX_NAMES[idx-1]}' "
            else:
                nxname_str += f"'{NX_NAMES[idx]}' "
        for idy in range(ny):
            if len(y_pos_with_refl)>0:
                if idy == max(y_pos_with_refl) or idy == min(y_pos_with_refl) and len(y_pos_with_refl)>0:
                    nyname_str += "'-' "
                else:
                    nyname_str += f"'{NY_NAMES[idy-1]}' "
            else:
                nyname_str += f"'{NY_NAMES[idy]}' "

        return nxname_str, nyname_str

    def format_burnup_zones(self):
        """
        Format B-ZONE entry
        For now assume one burnup zone per fuel channel and per height. 
        TODO : implement grouping of zones for more flexible burnup zone assignment.
        """

        nx = len(self.model_initialiser.meshx) - 1
        ny = len(self.model_initialiser.meshy) - 1
        nz = len(self.model_initialiser.meshz) - 1

        b_zones_str = ""
        zone_index = 1
        for idz in range(nz):
            for idy in range(ny):
                b_zones_str += "\n"
                for idx in range(nx):
                    if (idx, idy) in self.model_initialiser.fuel_channel_pairs:
                        b_zones_str += f"{zone_index} "
                        zone_index += 1

        return b_zones_str

    def format_parameters_definition(self):
        """
        Format paramters definition that should be included in the fuel map.
        """
        params_str = ""
        parameters = self.model_initialiser.interpolation_variables
        for parameter in parameters.keys():
            params_str += f"ADD-PARAM PNAME '{parameters[parameter][1]}' PARKEY '{parameter}' {parameters[parameter][0]} \n"
        return params_str

    def format_set_parameters_bloc(self, isInitialisation):
        """
        Format the SET-PARAM block 
        """
        set_params_block = ""
        parameters = self.model_initialiser.interpolation_variables
        for param in parameters.keys():
            if isInitialisation:
                # recover guess initial value for parameter
                inital_pvalue = self.model_initialiser.get_initial_parameter(param) 
                if parameters[param][0] == "GLOBAL":
                    set_params_block += f"SET-PARAM '{parameters[param][1]}' {inital_pvalue:.5E}\n"
                else: 
                    # set bundle-wise values : 1 value per channel, 1 value per axial mesh
                    num_chan = self.model_initialiser.number_fuel_channels
                    num_bundles = num_chan * (len(self.model_initialiser.meshz) - 1) 
                    param_values = [inital_pvalue] * num_bundles
                    set_params_block += f"SET-PARAM '{parameters[param][1]}' BUND \n"
                    set_params_block += f"{format_cle2000_float_array(param_values, max_per_line=5)}\n"
            else:
                param_values = self.model_updater.parameter_values_dict[param]
                if parameters[param][0] == "GLOBAL" and len(param_values)==1: 
                    set_params_block += f"SET-PARAM '{parameters[param][1]}' {param_values[0]:.5E}\n"
                else:
                    set_params_block += f"SET-PARAM '{parameters[param][1]}' BUND \n"
                    set_params_block += f"{format_cle2000_float_array(param_values, max_per_line=5)}\n"
        return set_params_block

    def format_fuel_weight_bloc(self):
        """
        This function will not work in this current implementation with vanished rods as it assumes that the fuel mass is evenly 
        distributed along the core's height.
        """
        fuel_weight_bloc = ""
        num_chan = self.model_initialiser.number_fuel_channels
        nz = (len(self.model_initialiser.meshz) - 1) 
        num_bundles = num_chan * nz
        fuel_weight_bloc += f"FUEL WEIGHT\n" 
        fuel_weight_values = [self.model_initialiser.total_fuel_mass / self.model_initialiser.nfuel] * num_bundles
        fuel_weight_bloc += f"{format_cle2000_float_array(fuel_weight_values, max_per_line=5)}\n"

        return fuel_weight_bloc




class NCR:
    """
        Call to the NCR: module performing interpolation of cross sections based on a set of parameters.
    """
    def __init__(self, libname, fuelmap_name, initialModel, neutronicsSolve):
        """
        libname (str) : CLE2000 variable name to host the output cross section library.
        fuelmap_name (str) : CLE2000 variable name hosting the fuel map to recover parameters from.
        initial_model : ModelInitialisation object carrying the variables to be included in the interpolation.
        neutronicsSolve : NeutronicsSolve object carrying the parameters used to interpolate cross sections on
        compo1 : CLE2000 variable name hosting the 
        """

        self.libname = libname
        self.fuelmap_name = fuelmap_name
        self.initial_model = initialModel
        self.neutronics_model = neutronicsSolve


    def write_c2m(self):
        # include second compo ? 
        # call NCR: for as many compos as there are and pair mixes to compos + DIR to interpolate from
        compos_dir_keys = self.neutronics_model.compo_dir_to_mix.keys()
        compo_names = [cpo_dir_pair[0] for cpo_dir_pair in compos_dir_keys]
        updateLib = False

        ncr_call_block = ""
        for compo_dir_pair in compos_dir_keys:
            cpo_lcm_name = compo_dir_pair[0]
            cpodir = compo_dir_pair[1]
            if updateLib:
                ncr_call_block += f"{self.libname} := NCR: {self.fuelmap_name} {self.libname} {cpo_lcm_name} ::\n"
                ncr_call_block += "    EDIT 0 \n"
            else:
                ncr_call_block += f"{self.libname} := NCR: {self.fuelmap_name} {cpo_lcm_name} ::\n"
                ncr_call_block += "    EDIT 0 \n"
                ncr_call_block += f"    NMIX {self.initial_model.nfuel}"

            updateLib = True
            ncr_call_block += f"    MICRO {self.neutronics_model.interpolation_type.upper()} \n"

            mix_def = self.format_compo_dir_pair_mix_def(cpo_lcm_name, cpodir)
            ncr_call_block += f"{mix_def}\n"
            ncr_call_block += ";\n"

        return ncr_call_block


    def format_compo_dir_pair_mix_def(self, compo_lcm, compo_dir):
        """
        Format MIX definitions from a CPO + CPODIR and a set of parameters
        """
        mix_definitions = ""
        mix_indices = self.neutronics_model.compo_dir_to_mix[(compo_lcm, compo_dir)]
        for mix_idx in mix_indices:
            mix_definitions += f"    COMPO {compo_lcm} {compo_dir}\n"
            mix_definitions += f"    MIX {mix_idx}\n"
            for par_key in self.neutronics_model.parameter_values_dict.keys():
                parameter_value = self.neutronics_model.parameter_values_dict[par_key][mix_idx-1]
                mix_definitions += f"    SET '{par_key}' {parameter_value:5E}\n"
            mix_definitions += "    ENDMIX\n"

        return mix_definitions

class THM:
    """
        Call to the THM: module responsible for providing a simplified thermal-hydraulics solution.
    """
    def __init__(self, inlet_temp, outlet_press, pdrop, dfm, 
                 fuel_radius, gap_radius, clad_radius, acool_profile, dh_profile, pch_profile,
                 kexp_profile=None, kcon_profile=None, rsin_profile=None,
                 wr_holes=None, idelchik_exit=None, idelchik_enter=None,
                 acool_wr=None, dh_wr=None, kexp_wr=None, kcon_wr=None, rsin_wr=None,
                 pch_wr_out=None, rwall_wr=None, split_wr=None):
        self.inlet_temp = inlet_temp
        self.outlet_press = outlet_press
        self.pdrop = pdrop
        self.dfm = dfm
        self.fuel_radius = fuel_radius
        self.gap_radius = gap_radius
        self.clad_radius = clad_radius
        self.acool_profile = acool_profile
        self.dh_profile = dh_profile
        self.pch_profile = pch_profile
        self.kexp_profile = kexp_profile
        self.kcon_profile = kcon_profile
        self.rsin_profile = rsin_profile
        self.wr_holes = wr_holes or []
        self.idelchik_exit = idelchik_exit or []
        self.idelchik_enter = idelchik_enter or []
        self.acool_wr = acool_wr
        self.dh_wr = dh_wr
        self.kexp_wr = kexp_wr
        self.kcon_wr = kcon_wr
        self.rsin_wr = rsin_wr
        self.pch_wr_out = pch_wr_out
        self.rwall_wr = rwall_wr
        self.split_wr = split_wr
        self.inlet_area = acool_profile[0]
        # Use the minimum ACOOL (active zone = most restricted) as mass-flow reference,
        # so plenum nodes at inlet do not artificially inflate the mass flow rate.
        _ref_acool = min(acool_profile)
        self.mass_flow = 8.407E-02 * (_ref_acool / 8.470E-05)

    def write_c2m(self):
        acool_str = format_cle2000_float_array(self.acool_profile)
        hd_str = format_cle2000_float_array(self.dh_profile)
        pch_str = format_cle2000_float_array(self.pch_profile)

        # KEXP-P block: only written when there are non-zero expansion losses
        if self.kexp_profile is not None and any(k != 0.0 for k in self.kexp_profile):
            kexp_p_str = format_cle2000_float_array(self.kexp_profile)
            kexp_p_block = "    KEXP-P\n" + f"{kexp_p_str}\n"
        else:
            kexp_p_block = ""

        # KCON-P and RSIN-P blocks: written together when there are non-zero contraction losses
        if self.kcon_profile is not None and any(k != 0.0 for k in self.kcon_profile):
            kcon_p_str = format_cle2000_float_array(self.kcon_profile)
            rsin_p_str = format_cle2000_float_array(self.rsin_profile)
            kcon_p_block = "    KCON-P\n" + f"{kcon_p_str}\n" + "    RSIN-P\n" + f"{rsin_p_str}\n"
        else:
            kcon_p_block = ""

        # WR-HOLE block: written when holes are defined
        # A_hole = pi*(D_hole/2)^2, D_hole in cm -> A_hole in m^2
        if self.wr_holes:
            hole_z_vals   = " ".join(f"{h['z'] * 1e-2:.5E}" for h in self.wr_holes)
            hole_a_vals   = " ".join(
                f"{np.pi * (h['D_hole'] * 0.5e-2) ** 2:.5E}" for h in self.wr_holes
            )
            wr_hole_block = (
                f"    WR-HOLE {len(self.wr_holes)}\n"
                f"    HOLE-Z {hole_z_vals}\n"
                f"    HOLE-A {hole_a_vals} (*m2*)\n"
            )
            if self.idelchik_exit:
                wr_hole_block += self._idelchik_line('IDELCHIK-EXIT', self.idelchik_exit)
            if self.idelchik_enter:
                wr_hole_block += self._idelchik_line('IDELCHIK-ENTER', self.idelchik_enter)
        else:
            wr_hole_block = ""

        # ACOOL-WR / HD-WR / KSING-WR blocks for water rod interior
        if self.acool_wr is not None:
            acool_wr_str = format_cle2000_float_array(self.acool_wr)
            acool_wr_block = "    ACOOL-WR\n" + f"{acool_wr_str}\n\n"
        else:
            acool_wr_block = ""

        if self.dh_wr is not None:
            dh_wr_str = format_cle2000_float_array(self.dh_wr)
            dh_wr_block = "    HD-WR\n" + f"{dh_wr_str}\n\n"
        else:
            dh_wr_block = ""

        if self.kexp_wr is not None and any(k != 0.0 for k in self.kexp_wr):
            kexp_wr_str = format_cle2000_float_array(self.kexp_wr)
            kexp_wr_block = "    KEXP-WR\n" + f"{kexp_wr_str}\n"
        else:
            kexp_wr_block = ""

        if self.kcon_wr is not None and any(k != 0.0 for k in self.kcon_wr):
            kcon_wr_str = format_cle2000_float_array(self.kcon_wr)
            rsin_wr_str = format_cle2000_float_array(self.rsin_wr)
            kcon_wr_block = "    KCON-WR\n" + f"{kcon_wr_str}\n" + "    RSIN-WR\n" + f"{rsin_wr_str}\n"
        else:
            kcon_wr_block = ""

        if self.pch_wr_out is not None:
            pch_wr_str = format_cle2000_float_array(self.pch_wr_out)
            pch_wr_block = "    PCH-WR-OUT\n" + f"{pch_wr_str}\n\n"
        else:
            pch_wr_block = ""

        if self.rwall_wr is not None:
            rwall_str = format_cle2000_float_array(self.rwall_wr)
            rwall_block = "    RWALL-WR\n" + f"{rwall_str}\n"
        else:
            rwall_block = ""

        # SPLIT-WR block: flow split between active channel and water rod
        #   split_wr='AUTO' -> bisection to equalize outlet pressures, initial guess x=0.05
        #   split_wr=<float> -> fixed fraction of total mass flow to water rod
        if self.split_wr is not None and (self.acool_wr is not None):
            if str(self.split_wr).upper() == 'AUTO':
                split_wr_block = "    SPLIT-WR AUTO 5.0E-02\n"
            else:
                split_wr_block = f"    SPLIT-WR {float(self.split_wr):.5E}\n"
        else:
            split_wr_block = ""

        block = (
            "*--------------------------------------------------------\n"
            "* THM single-stage calculation\n"
            "*--------------------------------------------------------\n"
            "Thm Fmap := THM: Fmap ::\n"
            "    EDIT 1\n"
            "    FLUID H2O\n"
            "    FPUISS 1.0\n"
            "    CRITFL 5.0E7\n"
            f"    INLET {self.outlet_press:.2E} (*Pa*) {self.inlet_temp:.2f} (*K*)\n"
            f"    INLET-Q {self.inlet_area:.5E} (*m2*) {self.mass_flow:.5E} (*kg/s*)\n"
            "    ASSMB 1 0\n"
            f"    RADIUS {self.fuel_radius:.3E} {self.gap_radius:.3E} {self.clad_radius:.3E} 0.000E+00 (*m*)\n"
            "    RODMESH 15 20\n"
            "    HGAP 10000.0\n"
            "    CONDC 0 21.5 KELVIN\n"
            "    CONDF 0 4.18 KELVIN\n"
            "    SAHA\n"
            f"    PDROP {self.pdrop}\n"
            f"    DFM {self.dfm}\n\n"
            "    * --- PROFILS AXIAUX GEOMETRIQUES (CANAL ACTIF) ---\n"
            "    ACOOL-P\n"
            f"{acool_str}\n\n"
            "    HD-P\n"
            f"{hd_str}\n\n"
            "    PCH-P\n"
            f"{pch_str}\n"
            f"{kexp_p_block}"
            f"{kcon_p_block}"
            "    * --- PROFILS AXIAUX GEOMETRIQUES (WATER ROD) ---\n"
            f"{acool_wr_block}"
            f"{dh_wr_block}"
            f"{kexp_wr_block}"
            f"{kcon_wr_block}"
            f"{pch_wr_block}"
            f"{rwall_block}"
            f"{split_wr_block}"
            f"{wr_hole_block}"
            "    * ------------------------------------------------\n"
            ";\n"
        )
        return block

    def _idelchik_line(self, keyword, table, max_len=119):
        """
            Format an IDELCHIK-EXIT/ENTER keyword line, wrapping at max_len chars
        """
        prefix = f"    {keyword} {len(table)} "
        cont   = "        "  # continuation indent (must be shorter than prefix)
        tokens = [f"{pt[0]:.4f} {pt[1]:.4f}" for pt in table]
        lines  = []
        current = prefix
        for tok in tokens:
            candidate = current + tok + " "
            if len(candidate) > max_len and current != prefix:
                lines.append(current.rstrip())
                current = cont + tok + " "
            else:
                current = candidate
        lines.append(current.rstrip())
        return "\n".join(lines) + "\n"

######################################################################################################################################
##############################                         Procedure orchestrators :                           ###########################
######################################################################################################################################


# -----------------------------------------------------------------------------------------------
#   Procedure to initialise a DONJON Case
# -----------------------------------------------------------------------------------------------
class INIT:
    """
    Wrapper class creating a CLE2000 procedure initalizing the Donjon Model
    """

    def __init__(self, scheme, proc_name = "IniDonjon", case_name="donjon_case"):
        """
        scheme : DonjonCalculationScheme
        """
        self.scheme = scheme
        self.geometry_lcm_name = "Geom"
        self.matex_lcm_name = "Matex"
        self.fmap_lcm_name = "Fmap"
        self.proc_name = proc_name
        self.case_name = case_name

        self.geo_block = self.build_GEO_call()
        self.matex_block = self.build_MATEX_call()
        self.resini_block = self.build_RESINI_call()


        

    def build_GEO_call(self):
        """
        Build a call to GEO: initializing the geometry model for the donjon case 
        """
        geo_proc = GEO(lcm_geometry_name=self.geometry_lcm_name, model_initialisation=self.scheme.initialisation_step)
        geom_block = geo_proc.write_c2m()
        return geom_block

    def build_MATEX_call(self):
        """
        Build a call to USPLIT: performing material indexation and producing a Matex lcm object.
        """
        usplit_proc = USPLIT(geometry_name=self.geometry_lcm_name, matex_name=self.matex_lcm_name, model_initialisation=self.scheme.initialisation_step)
        usplit_block = usplit_proc.write_c2m()
        return usplit_block

    def build_RESINI_call(self):
        """
        Build a call to RESINI: 
        """
        resini_proc = RESINI(fuel_map_name=self.fmap_lcm_name, matex_name=self.matex_lcm_name, 
                             model_initialiser = self.scheme.initialisation_step, 
                             )
        resini_block = resini_proc.write_c2m()
        return resini_block

    def build_initialisation_procedure(self):
        """build the body of the IniDonjon.c2m proceudre"""
        body = ""
        if self.geo_block:
            body += (
                f"* Build geometry for {self.case_name}\n"
            )
            body += self.geo_block
            body += "\n"
        if self.matex_block:
            body += (
                f"* Build material indexation for {self.case_name}\n"
            )
            body += self.matex_block
            body += "\n"
        if self.resini_block:
            body += (
                f"* Build fuel map definition for {self.case_name}\n"
            )
            body += self.resini_block
            body += "\n"
        return body

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the complete TRK sub-procedure to a ``.c2m`` file.
        """
        from .CLE2000 import (
            CLE2000_MAX_LINE, CLE2000_MAX_VARNAME,
            validate_varname, wrap_cle2000_line,
        )


        # --- PARAMETER block ---
        param_items = ["FMap", "Matex", "Cpo", "Track", "THData"]
        
        header = (
            f"* PROCEDURE {proc_name}.c2m : problem initialisation\n"
            "* --------------------------------\n"
            "* Procedure generated by starterDD\n"
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
        )

        # PARAMETER declaration
        param_block = "PARAMETER"
        for item in param_items:
            param_block += f" {item}"
        param_block += " ::\n"
        # Linked-list declarations
        param_block += f"::: LINKED_LIST"
        for lcm_obj in param_items:
            param_block += f" {lcm_obj} "
        
        param_block += "; ;\n"

        # MODULE declaration
        mod_block = "MODULE GEO: USPLIT: RESINI: END: ;\n"

        body = self.build_initialisation_procedure()

        footer = "END: ;\nQUIT .\n"

        content = (
            f"{header}{param_block}\n"
            f"{mod_block}\n{body}{footer}"
        )

        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(
            path_to_procs, f"{proc_name}.c2m"
        )
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[INIT] Wrote procedure to {filepath}")
        return filepath



# -----------------------------------------------------------------------------------------------
#    Procedure solving the neutronics problem
# -----------------------------------------------------------------------------------------------

class NEUTRONICS:
    """
    Wrapper class orchestrating the neutronics solution
    """

    def __init__(self, scheme, proc_name = "Neutronics", case_name = "Donjon_case"):
        self.scheme = scheme

    def build_RESINI_call(self):
        """
        Build a call to the RESINI: module to update Fuel Map parameters
        """
        
        resini_proc = RESINI(fuel_map_name="FMap", matex_name="Matex", model_initialiser=self.scheme.initialisation_step, model_updater=self.scheme.neutronics_step)
        resini_block = resini_proc.write_c2m()
        return resini_block

    def build_NCR_call(self):
        """
        Build a call to the NCR: interpolation module
        """
        ncr_proc = NCR(libname="MicroF", fuelmap_name="FMap", initialModel= self.scheme.initialisation_step, neutronicsSolve=self.scheme.neutronics_step)
        ncr_block = ncr_proc.write_c2m()
        return ncr_block


# -----------------------------------------------------------------------------------------------
# Test Case Orchestrator for THM / Donjon
# -----------------------------------------------------------------------------------------------

class DonjonTHM1DProcedure:
    """
    Orchestrator class that instantiates DONJON modules and generates the final .c2m file.

    Optional flags
    --------------
    use_acool_profile : bool (default True)
        If True, use the full axial ACOOL-P / HD-P / PCH-P profiles from the analyser.
        If False, replace them with constant arrays equal to the mid-channel value
        (overrides all other geometry flags).  Useful for fully-uniform channel validation.
    use_ksing : bool (default True)
        If True, include the KSING-P (singular loss) profile from the analyser.
        If False, set KSING-P = 0 everywhere (overrides all other KSING flags).
    use_spacer_grids : bool (default True)
        If True, spacer-grid regions (axial_region == "GRID") are included as-is.
        If False, ACOOL-P / HD-P / PCH-P at GRID nodes are replaced by the value
        of the immediately preceding non-GRID node (forward-fill), and KSING-P is
        zeroed at all boundaries involving a GRID tranche.  The vanished-rod geometry
        (VAN region) is unaffected by this flag.
    use_vanished_rods : bool (default True)
        If True, the vanished-rod region (axial_region == "VAN") is included as-is.
        If False, all nodes at z >= z_van_start (including GRID nodes within the VAN
        territory) have their ACOOL-P / HD-P / PCH-P replaced by the DOM reference
        value (the last stable DOM-region value before the VAN zone), and KSING-P is
        zeroed for all inter-tranche boundaries within that territory.  The spacer-grid
        geometry (GRID nodes in the DOM region) is unaffected by this flag.
    use_wr : bool (default True)
        If True, include all water-rod geometry and coupling blocks.
        If False, skip the WR entirely (equivalent to a standard active channel).
    use_wr_holes : bool (default True)
        If True (and use_wr=True), include the lateral orifice model (Idelchik
        tables + SPLIT-WR) that drives mass-flow exchange between the active
        channel and the WR.
        If False, the WR is still thermally coupled (heat exchange through the
        wall via pch_wr_out / rwall_wr) but there is no lateral mass-flow
        exchange.  Useful for isolating the thermal coupling from the hydraulic
        coupling when validating the WR solver.
        Has no effect when use_wr=False.

    Flag independence
    -----------------
    use_spacer_grids and use_vanished_rods are fully independent:

    +------------------+-------------------+------------------------------------------+
    | use_spacer_grids | use_vanished_rods | Effect                                   |
    +------------------+-------------------+------------------------------------------+
    | True             | True              | full geometry (reference)                |
    | False            | True              | no spacers; VAN region kept as-is        |
    | True             | False             | spacers kept (incl. in VAN territory);   |
    |                  |                   | VAN non-GRID nodes → DOM ref             |
    | False            | False             | all variations removed; uniform channel  |
    +------------------+-------------------+------------------------------------------+

    When use_spacer_grids=True and use_vanished_rods=False, KSING at spacer
    boundaries in the VAN territory is kept (computed from original VAN/GRID area
    ratios — slightly approximate since surrounding area is now DOM-equivalent).

    use_acool_profile=False and use_ksing=False are applied last and override the
    region-based filters above.
    """
    def __init__(self, model_initialiser, analyser, nz, power_kw,
                 axial_pform=None, pdrop=1, dfm=1,
                 inlet_temp=543.15, outlet_press=7.20E+06,
                 use_acool_profile=True, use_ksing=True,
                 use_spacer_grids=True, use_vanished_rods=True,
                 use_wr=True, use_wr_holes=True):
        
        self.nz = nz
        
        # --- 1. Data preparation ---
        _, pitch_cm = analyser.get_x_global_bounds()
        pin_geom = analyser.slices_data[0]['dragon_assembly_model'].pin_geometry_dict
        fuel_radius = pin_geom['fuel_radius'] * 1E-2
        gap_radius = pin_geom['gap_radius'] * 1E-2
        clad_radius = pin_geom['clad_radius'] * 1E-2

        
        z_min, maxh = analyser.get_z_global_bounds()
        z_bounds = np.linspace(z_min, maxh, self.nz + 1).tolist()
        
        dz = (maxh - z_min) / self.nz
        geom_profiles = analyser.execute_profile_z(
            ['cv', [0, 0, pitch_cm, pitch_cm]], 
            dz, dz, z_min, maxh
        )
        
        acool_profile  = [a / 10000.0 for a in geom_profiles[2]] 
        dh_profile     = [dh * 1e-2   for dh in geom_profiles[3]]
        pch_profile    = [pch * 1e-2  for pch in geom_profiles[4]]
        # kexp/kcon/rsin lists are dimensionless, no unit conversion needed
        kexp_profile  = list(geom_profiles[5])
        kcon_profile  = list(geom_profiles[6])
        rsin_profile  = list(geom_profiles[7])

        # --- Apply geometry flags (region-based, then global overrides) ---
        # Order: use_vanished_rods first, then use_spacer_grids.
        # When both are False, GRID_VAN nodes forward-fill from the already-DOM-
        # substituted preceding value, giving a fully uniform channel.

        # Build node-to-region mapping from analyser.slices_data
        z_node_mids = [z_min + (k + 0.5) * dz for k in range(nz)]

        def _region_at(z):
            for t in analyser.slices_data:
                if t['z_start'] - 1e-9 <= z <= t['z_end'] + 1e-9:
                    return t.get('axial_region', '')
            return ''

        def _node_region(z_mid):
            # A node whose range [z_mid-dz/2, z_mid+dz/2] overlaps ANY GRID tranche
            # is classified as GRID, even if the midpoint itself is not inside
            # (handles GRIDs narrower than dz).
            z_lo = z_mid - dz / 2
            z_hi = z_mid + dz / 2
            for t in analyser.slices_data:
                if t.get('axial_region') == 'GRID':
                    if t['z_start'] < z_hi - 1e-9 and t['z_end'] > z_lo + 1e-9:
                        return 'GRID'
            return _region_at(z_mid)

        node_regions = [_node_region(z) for z in z_node_mids]

        # Build inter-tranche boundary list: (z_b, node_k, reg_before, reg_after)
        boundaries = []
        for k_idx in range(len(analyser.slices_data) - 1):
            t_b = analyser.slices_data[k_idx]
            t_a = analyser.slices_data[k_idx + 1]
            z_b = float(t_b['z_end'])
            if not (z_min < z_b < z_min + nz * dz):
                continue
            node_k = min(int((z_b - z_min) / dz), nz - 1)
            boundaries.append((z_b, node_k,
                                t_b.get('axial_region', ''),
                                t_a.get('axial_region', '')))

        # -- use_vanished_rods=False: replace only nodes whose axial_region is exactly
        #    'VAN' with the last DOM reference value.  GRID nodes are intentionally
        #    left unchanged — whether they are in DOM or VAN territory — so that
        #    use_spacer_grids can manage them independently.  This means a GRID_VAN
        #    tranche (spacer straddling the DOM-VAN boundary on the VAN side) is
        #    correctly preserved when use_spacer_grids=True, and forward-filled from
        #    already-substituted DOM values when use_spacer_grids=False (applied next).
        #
        #    KSING is zeroed at boundaries where BOTH adjacent regions are 'VAN' and
        #    at the direct DOM→VAN boundary (if no spacer separates them).
        #    Boundaries involving any GRID tranche are left for use_spacer_grids.
        #    'VAN' in (rb, ra) covers: DOM→VAN, VAN→VAN, and VAN→PLENUM entries;
        #    'GRID' not in (rb, ra) excludes all spacer-adjacent boundaries.
        if not use_vanished_rods:
            dom_ref_a = dom_ref_dh = dom_ref_pch = None
            for k, region in enumerate(node_regions):
                if region == 'DOM':
                    dom_ref_a   = acool_profile[k]
                    dom_ref_dh  = dh_profile[k]
                    dom_ref_pch = pch_profile[k]
                elif region == 'VAN':
                    # Replace with the last DOM value seen so far.
                    if dom_ref_a is not None:
                        acool_profile[k] = dom_ref_a
                        dh_profile[k]    = dom_ref_dh
                        pch_profile[k]   = dom_ref_pch
                # 'GRID' nodes: skip — handled by use_spacer_grids.
            # Zero KEXP/KCON only at non-GRID boundaries that involve the VAN region.
            for z_b, node_k, rb, ra in boundaries:
                if 'VAN' in (rb, ra) and 'GRID' not in rb and 'GRID' not in ra:
                    kexp_profile[node_k] = 0.0
                    kcon_profile[node_k] = 0.0

        # -- use_spacer_grids=False: replace GRID nodes with the preceding non-GRID
        #    value (forward-fill); zero KSING at all boundaries involving a GRID tranche.
        #    Applied after use_vanished_rods so that GRID-VAN nodes forward-fill from
        #    the already-substituted DOM value when both flags are False.
        if not use_spacer_grids:
            last_a, last_dh, last_pch = acool_profile[0], dh_profile[0], pch_profile[0]
            for k, region in enumerate(node_regions):
                if 'GRID' not in region:
                    last_a   = acool_profile[k]
                    last_dh  = dh_profile[k]
                    last_pch = pch_profile[k]
                else:
                    acool_profile[k] = last_a
                    dh_profile[k]    = last_dh
                    pch_profile[k]   = last_pch
            for z_b, node_k, rb, ra in boundaries:
                if 'GRID' in rb or 'GRID' in ra:
                    kexp_profile[node_k] = 0.0
                    kcon_profile[node_k] = 0.0

        # -- use_acool_profile=False: override with fully uniform channel (mid value).
        #    Applied last so it overrides all region-based filters above.
        if not use_acool_profile:
            mid = nz // 2
            acool_profile = [acool_profile[mid]] * nz
            dh_profile    = [dh_profile[mid]]    * nz
            pch_profile   = [pch_profile[mid]]   * nz

        # -- use_ksing=False: zero ALL singular losses.
        #    Applied last so it overrides all region-based filters above.
        if not use_ksing:
            kexp_profile = [0.0] * nz
            kcon_profile = [0.0] * nz

        # --- Water rod interior hydraulic profiles ---
        if use_wr:
            wr_profiles = analyser.execute_profile_z(
                ('wr_tube',),
                dz, dz, z_min, maxh
            )
            acool_wr = [a / 10000.0 for a in wr_profiles[2]]
            dh_wr    = [dh * 1e-2   for dh in wr_profiles[3]]
            kexp_wr  = list(wr_profiles[5])
            kcon_wr  = list(wr_profiles[6])
            rsin_wr  = list(wr_profiles[7])

            # Outer perimeter and wall resistance profiles for WR thermal coupling
            pch_wr_out = []
            rwall_wr   = []
            curr_z = z_min
            while curr_z + dz <= maxh + 1e-10:
                z1, z2 = curr_z, curr_z + dz
                pch_wr_out.append(analyser.get_pch_wr_outer(z1, z2) * 1e-2)
                rwall_wr.append(analyser.get_rwall_wr_tube(z1, z2))
                curr_z += dz

            wr_holes, idelchik_exit, idelchik_enter = analyser.get_wr_hole_data()
            if not use_wr_holes:
                wr_holes = []
                idelchik_exit = []
                idelchik_enter = []
        else:
            acool_wr = dh_wr = kexp_wr = kcon_wr = rsin_wr = pch_wr_out = rwall_wr = None
            wr_holes = idelchik_exit = idelchik_enter = []
        
        power_mw = power_kw / 1000.0
        axial_pform = axial_pform if axial_pform else [1.0] * self.nz

        # --- 2. Instantiation of DONJON modules ---
        self.geo_module = GEO(lcm_geometry_name="Geom", model_initialisation=model_initialiser)
        self.resini_module = RESINI(fuel_map_name="FMap", matex_name="Matex", model_initialiser=model_initialiser)
        self.thm_module = THM(
            inlet_temp=inlet_temp, outlet_press=outlet_press, 
            pdrop=pdrop, dfm=dfm, fuel_radius=fuel_radius, gap_radius=gap_radius, clad_radius=clad_radius,
            acool_profile=acool_profile, dh_profile=dh_profile, pch_profile=pch_profile,
            kexp_profile=kexp_profile, kcon_profile=kcon_profile, rsin_profile=rsin_profile,
            wr_holes=wr_holes, idelchik_exit=idelchik_exit, idelchik_enter=idelchik_enter,
            acool_wr=acool_wr, dh_wr=dh_wr, kexp_wr=kexp_wr, kcon_wr=kcon_wr, rsin_wr=rsin_wr,
            pch_wr_out=pch_wr_out, rwall_wr=rwall_wr,
            split_wr='AUTO' if use_wr else None
        )

    def write_to_c2m(self, path_to_procs, proc_name, pinlet_ref=None, epsout_ref=None):
        header = (
            "*************************************************************************\n"
            f"* Input file : {proc_name} \n"
            "* Generated by starterDD - THM 1D Test Case\n"
            "*************************************************************************\n\n"
            "LINKED_LIST Geom Matex Fmap Thm ;\n"
            "MODULE GEO: RESINI: USPLIT: THM: GREP: UTL: DELETE: ABORT: END: ;\n\n"
        )

        if pinlet_ref is not None and epsout_ref is not None:
            assert_block = (
                f"REAL PINLET EPSOUT DELTA ;\n"
                f"REAL REFPVAL := {pinlet_ref:.6E} ;\n"
                f"REAL REFEPSVAL := {epsout_ref:.6E} ;\n\n"
                "GREP: Thm :: STEP UP 'HISTORY-DATA' STEP UP 'TIMESTEP0000' STEP UP 'CHANNEL'\n"
                "            STEP AT 1 GETVAL 'PINLET' 1 1 1 >>PINLET<< ;\n"
                "EVALUATE DELTA := PINLET REFPVAL - REFPVAL / ABS ;\n"
                "IF DELTA 1.0E-2 < THEN\n"
                f"  PRINT \"TEST SUCCESSFUL ({proc_name} PINLET); DELTA=\" DELTA ;\n"
                "ELSE\n"
                "  PRINT \"TEST FAILURE\" ;\n"
                "  PRINT \"REFERENCE=\" REFPVAL \" CALCULATED=\" PINLET ;\n"
                "  ABORT: ;\n"
                "ENDIF ;\n\n"
                "GREP: Thm :: STEP UP 'HISTORY-DATA' STEP UP 'TIMESTEP0000' STEP UP 'CHANNEL'\n"
                "            STEP AT 1 GETVAL 'EPSOUT' 1 1 1 >>EPSOUT<< ;\n"
                "EVALUATE DELTA := EPSOUT REFEPSVAL - REFEPSVAL / ABS ;\n"
                "IF DELTA 1.0E-2 < THEN\n"
                f"  PRINT \"TEST SUCCESSFUL ({proc_name} EPSOUT); DELTA=\" DELTA ;\n"
                "ELSE\n"
                "  PRINT \"TEST FAILURE\" ;\n"
                "  PRINT \"REFERENCE=\" REFEPSVAL \" CALCULATED=\" EPSOUT ;\n"
                "  ABORT: ;\n"
                "ENDIF ;\n\n"
            )
        else:
            assert_block = ""

        content = (
            f"{header}"
            f"{self.geo_module.write_c2m()}\n"
            f"{self.resini_module.write_c2m()}\n"
            f"{self.thm_module.write_c2m()}\n"
            "UTL: Thm :: DIR DUMP ;\n\n"
            "Thm := UTL: Thm :: STEP UP 'HISTORY-DATA' STEP UP 'TIMESTEP0000'"
            " STEP UP 'CHANNEL' STEP AT 1 DUMP ;\n\n"
            f"{assert_block}"
            "END: ;\nQUIT .\n"
        )

        if not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(path_to_procs, f"{proc_name}.c2m")
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[DONJON_THM] Procédure générée avec succès : {filepath}")
        return filepath