## collection of classes to handle DONJON module calls
# Author: B. Godard
#Date: 16/04/2026 (creation)
# Purpose : Define class structures to handle DONJON module calls for starterDD package
# ----------------------------------------------------------------------------------------------- 

import os
import numpy as np


def format_cle2000_array(vector, max_per_line=10):
    """Formate une liste Python en un bloc de texte CLE-2000."""
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


# -----------------------------------------------------------------------------------------------
# Classes representing DONJON modules
# -----------------------------------------------------------------------------------------------

class GEO:
    """Génère l'appel au module GEO: (et USPLIT:) pour la géométrie."""
    def __init__(self, nz, pitch_cm, z_bounds):
        self.nz = nz
        self.pitch_cm = pitch_cm
        self.z_bounds = z_bounds

    def write_c2m(self):
        z_str = " ".join([f"{z:.4f}" for z in self.z_bounds])
        block = (
            f"Geom := GEO: :: CAR3D 1 1 {self.nz}\n"
            " X- REFL X+ REFL Y- REFL Y+ REFL Z- REFL Z+ REFL\n"
            f" MESHX 0.0 {self.pitch_cm:.5f}\n"
            f" MESHY 0.0 {self.pitch_cm:.5f}\n"
            f" MESHZ {z_str}\n"
            " MIX\n"
        )
        for i in range(1, self.nz + 1):
            if i == 1:
                block += f"  PLANE {i}  1\n"
            else:
                block += f"  PLANE {i} SAME 1\n"
                
        block += " ;\n"
        block += "Geom Matex := USPLIT: Geom :: NGRP 2 MAXR 10000 NFUEL 1 FMIX 1 ;\n"
        return block


class RESINI:
    """Génère l'appel au module RESINI: pour la création de la Fuel Map."""
    def __init__(self, nz, power_mw, axial_pform, pitch_cm, z_bounds):
        self.nz = nz
        self.power_mw = power_mw
        self.axial_pform = axial_pform
        self.pitch_cm = pitch_cm
        self.z_bounds = z_bounds

    def write_c2m(self):
        pform_str = format_cle2000_array(self.axial_pform)
        z_str = " ".join([f"{z:.4f}" for z in self.z_bounds])

        geo_inner = (
            " X- REFL X+ REFL Y- REFL Y+ REFL Z- REFL Z+ REFL\n"
            f" MESHX 0.0 {self.pitch_cm:.5f}\n"
            f" MESHY 0.0 {self.pitch_cm:.5f}\n"
            f" MESHZ {z_str}\n"
            " MIX\n"
        )
        for i in range(1, self.nz + 1):
            if i == 1:
                geo_inner += f"  PLANE {i}  1\n"
            else:
                geo_inner += f"  PLANE {i} SAME 1\n"

        block = (
            "*--------------------------------------------------------\n"
            "* Fuel map definition\n"
            "*--------------------------------------------------------\n"
            "Fmap Matex := RESINI: Matex ::\n"
            f"    ::: GEO: CAR3D 1 1 {self.nz}\n"
            f"{geo_inner}"
            "    NXNAME '01' NYNAME 'A' NCOMB 1 B-ZONE 1\n\n"
            "    ADD-PARAM PNAME 'T-FUEL' PARKEY 'TFuel' GLOBAL\n"
            "    ADD-PARAM PNAME 'T-COOL' PARKEY 'TCool' GLOBAL\n"
            "    ADD-PARAM PNAME 'D-COOL' PARKEY 'DCool' GLOBAL\n"
            "    BTYPE INST-BURN INST-BVAL CHAN 0.0\n"
            f"    REACTOR-POW {self.power_mw:.5E} AXIAL-PFORM\n"
            f"{pform_str}\n"
            "    SET-PARAM 'T-FUEL' 900.0\n"
            "    SET-PARAM 'T-COOL' 543.15\n"
            "    SET-PARAM 'D-COOL' 0.65\n"
            "    FUEL WEIGHT 6.464E-3\n"
            ";\n"
        )
        return block


class THM:
    """Génère l'appel au module THM: pour le calcul thermohydraulique."""
    def __init__(self, inlet_temp, outlet_press, pdrop, dfm, 
                 acool_profile, dh_profile, pch_profile):
        self.inlet_temp = inlet_temp
        self.outlet_press = outlet_press
        self.pdrop = pdrop
        self.dfm = dfm
        self.acool_profile = acool_profile
        self.dh_profile = dh_profile
        self.pch_profile = pch_profile

    def write_c2m(self):
        acool_str = format_cle2000_array(self.acool_profile)
        hd_str = format_cle2000_array(self.dh_profile)
        pch_str = format_cle2000_array(self.pch_profile)

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
            "    INLET-Q 8.470E-05 (*m2*) 8.407E-02 (*kg/s*)\n"
            "    ASSMB 1 0\n"
            "    RADIUS 4.435E-03 4.520E-03 5.140E-03 0.000E+00 (*m*)\n"
            "    RODMESH 15 20\n"
            "    HGAP 10000.0\n"
            "    CONDC 0 21.5 KELVIN\n"
            "    CONDF 0 4.18 KELVIN\n"
            "    SAHA\n"
            f"    PDROP {self.pdrop}\n"
            f"    DFM {self.dfm}\n\n"
            "    * --- NOUVEAUX PROFILS AXIAUX GEOMETRIQUES ---\n"
            "    ACOOL-P\n"
            f"{acool_str}\n\n"
            "    HD-P\n"
            f"{hd_str}\n\n"
            "    PCH-P\n"
            f"{pch_str}\n"
            "    * --------------------------------------------\n"
            ";\n"
        )
        return block


# -----------------------------------------------------------------------------------------------
# Test Case Orchestrator
# -----------------------------------------------------------------------------------------------

class DonjonTHM1DProcedure:
    """
    Orchestrator class that instantiates DONJON modules and generates the final .c2m file.
    """
    def __init__(self, analyser, nz, power_kw, 
                 axial_pform=None, pdrop=1, dfm=1, 
                 inlet_temp=543.15, outlet_press=7.20E+06):
        
        self.nz = nz
        
        # --- 1. Data preparation ---
        ass_geo = analyser.data_ref['ASSEMBLY_GEOMETRY']
        pitch_cm = ass_geo['assembly_pitch'] 
        
        z_min, maxh = analyser.get_z_global_bounds()
        z_bounds = np.linspace(z_min, maxh, self.nz + 1).tolist()
        
        dz = maxh / self.nz
        geom_profiles = analyser.execute_profile_z(
            ['cv', [0, 0, pitch_cm, pitch_cm]], 
            dz, dz, 0.0, maxh
        )
        
        acool_profile = [a / 10000.0 for a in geom_profiles[2]] 
        dh_profile = geom_profiles[3]
        pch_profile = geom_profiles[4]
        
        power_mw = power_kw / 1000.0
        axial_pform = axial_pform if axial_pform else [1.0] * self.nz

        # --- 2. Instantiation of DONJON modules ---
        self.geo_module = GEO(nz=self.nz, pitch_cm=pitch_cm, z_bounds=z_bounds)
        self.resini_module = RESINI(
            nz=self.nz, power_mw=power_mw, 
            axial_pform=axial_pform, pitch_cm=pitch_cm, z_bounds=z_bounds
        )
        self.thm_module = THM(
            inlet_temp=inlet_temp, outlet_press=outlet_press, 
            pdrop=pdrop, dfm=dfm,
            acool_profile=acool_profile, dh_profile=dh_profile, pch_profile=pch_profile
        )

    def write_to_c2m(self, path_to_procs, proc_name):
        header = (
            "*************************************************************************\n"
            f"* Input file : {proc_name} \n"
            "* Generated by starterDD - THM 1D Test Case\n"
            "*************************************************************************\n\n"
            "LINKED_LIST Geom Matex Fmap Thm ;\n"
            "MODULE GEO: RESINI: USPLIT: THM: GREP: UTL: DELETE: ABORT: END: ;\n\n"
        )
        
        content = (
            f"{header}"
            f"{self.geo_module.write_c2m()}\n"
            f"{self.resini_module.write_c2m()}\n"
            f"{self.thm_module.write_c2m()}\n"
            "END: ;\nQUIT .\n"
        )

        if not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(path_to_procs, f"{proc_name}.c2m")
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[DONJON_THM] Procédure générée avec succès : {filepath}")
        return filepath