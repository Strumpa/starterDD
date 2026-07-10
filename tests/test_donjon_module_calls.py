import os
import math
import numpy as np
import pytest

from starterDD.GeometryAnalysis.geometry_analysis import GeometricAnalyser
from starterDD.InterfaceToDD.donjon_module_calls import DonjonTHM1DProcedure
from starterDD.DDModel.DonjonModel import CoreModel

from conftest import GE14_CORE_YAML, OUTPUTS_DIR

def generate_case_name(pdrop, power_kw, profile_type):
    name = "dfm"
    if pdrop == 1: name += "p"
    name += str(int(power_kw))
    if profile_type.lower() in ['cosinus', 'cos', 'c']: name += "c"
    elif profile_type.lower() in ['sinus', 'sin', 's']: name += "s"
    else: name += "u"
    return name

def generate_power_profile(profile_type, nz):
    profile = []
    for i in range(nz):
        z_norm = (i + 0.5) / nz 
        if profile_type == 'cosinus':
            # Quart de cosinus : max en bas (z=0), nul en haut (z=1)
            val = math.cos((math.pi / 2.0) * z_norm) 
        elif profile_type == 'sinus':
            # Demi-sinus : cloche symétrique max au centre
            val = math.sin(math.pi * z_norm)
        else:
            val = 1.0
        profile.append(val)
    
    mean = sum(profile) / nz
    return [v / mean for v in profile]

def test_geometric_analyser():

    core_pos = (1, 1)
    # create a CoreModel from DonjonModel
    minicore_model = CoreModel(name="GE14_4x4_minicore", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml="GEOM_4x4_mini_CORE.yaml")
    minicore_model.createAssemblyModels()
    analyser = GeometricAnalyser(minicore_model, core_i=core_pos[0], core_j=core_pos[1])

    assert analyser.data_ref['ASSEMBLY_GEOMETRY']['assembly_pitch'] == 15.24
    assert analyser.data_ref['WATER_ROD_GEOMETRY']['centers'][0] == (8.920000000000002, 6.320000000000001)
    assert len(analyser.slices_data) == 2
    assert analyser.slices_data[0]["z_start"] == 0.0
    assert analyser.slices_data[0]["z_end"] == 222.0595
    assert analyser.slices_data[1]["z_start"] == 222.0595
    assert analyser.slices_data[1]["z_end"] == 347.1291

def test_DONJON_THM_generator_and_porosity_calculation():

    core_pos = (1, 1)
    # create a CoreModel from DonjonModel
    minicore_model = CoreModel(name="GE14_4x4_minicore", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml="GEOM_4x4_mini_CORE.yaml")
    minicore_model.createAssemblyModels()
    analyser = GeometricAnalyser(minicore_model, core_i=core_pos[0], core_j=core_pos[1])
    # --- TEST CASE MATRIX ---
    power = 20.0
    profiles_to_test = ['cosinus', 'sinus']
    pdrop = 1
    nz = 40
    dfm = 1
    
    # loop over possible combinations
    for profile_type in profiles_to_test:
        
        # 1. Automated case name generation
        case_name = generate_case_name(pdrop, power, profile_type)
        
        # 2. Generate the power profile
        axial_pform = generate_power_profile(profile_type, nz)
        
        # 3. Create a c2m procedure
        procedure = DonjonTHM1DProcedure(
            analyser=analyser,
            nz=nz,
            power_kw=power,
            axial_pform=axial_pform,
            pdrop=pdrop,
            dfm=dfm
        )

        assert len(procedure.thm_module.acool_profile) == nz
        assert len(procedure.thm_module.pch_profile) == nz
        assert len(procedure.thm_module.dh_profile) == nz

        assert procedure.thm_module.acool_profile[0] == 0.009237286
        assert procedure.thm_module.acool_profile[-1] == 0.010403689

        assert procedure.thm_module.acool_profile[-1] - procedure.thm_module.acool_profile[0] == pytest.approx(14*np.pi*(0.515e-2)**2, abs=1e-5)
        print(procedure.thm_module.acool_profile)
        assert procedure.thm_module.acool_profile[25] == 0.009717707
        assert procedure.thm_module.acool_profile[26] == 0.010403689

        # 4. Écriture du fichier
        c2m_path = procedure.write_to_c2m(OUTPUTS_DIR, case_name)