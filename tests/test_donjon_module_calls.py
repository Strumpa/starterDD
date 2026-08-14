import os
import math
import numpy as np
import pytest

from conftest import GE14_CORE_YAML, OUTPUTS_DIR

from starterDD.GeometryAnalysis.cartesian_geometry_analysis import CartesianGeometricAnalyser
from starterDD.InterfaceToDD.donjon_module_calls import DonjonTHM1DProcedure
from starterDD.DDModel.DonjonModel import CoreModel 

def generate_case_name(pdrop, power_kw, profile_type):
    nom = "dfm"
    if pdrop == 1: nom += "p"
    nom += str(int(power_kw))
    if profile_type.lower() in ['cosine', 'cos', 'c']: nom += "c"
    elif profile_type.lower() in ['sine', 'sin', 's']: nom += "s"
    else: nom += "u"
    nom += "_v"
    return nom

def generate_power_profile(profile_type, nz):
    profile = []
    for i in range(nz):
        z_norm = (i + 0.5) / nz 
        if profile_type == 'cosine':
            # Quarter cosine : max at z=0, 0 at z=zmax
            val = math.cos((math.pi / 2.0) * z_norm) 
        elif profile_type == 'sine':
            # Half-sine : 0 at z=0, max at z=zmax/2, 0 at z=zmax 
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
    analyser = CartesianGeometricAnalyser(minicore_model, core_i=core_pos[0], core_j=core_pos[1])

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
    analyser = CartesianGeometricAnalyser(minicore_model, core_i=core_pos[0], core_j=core_pos[1])
    # --- TEST CASE MATRIX ---
    power = 20.0
    profiles_to_test = ['cosine', 'sine']
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
        
def test_4x4_minicore_surfaces_definition():

    # Reference values (PINLET, EPSOUT) for non-regression assertions.
    REFERENCE_VALS = {
        'dfmp10s_v': (7.23205200E+06, 6.43932670E-02),
        'dfmp20s_v': (7.23644300E+06, 4.47797000E-01),
        'dfmp40s_v': (7.25449950E+06, 6.77849114E-01),
        'dfmp10c_v': (7.23193850E+06, 6.07140325E-02),
        'dfmp20c_v': (7.23913300E+06, 4.56152081E-01),
        'dfmp40c_v': (7.26596800E+06, 6.83239937E-01),
        'dfm10c_v':  (7.20000000E+06, 6.04344234E-02),
        'dfm10s_v':  (7.20000000E+06, 6.40267283E-02),
        'dfm20c_v':  (7.20000000E+06, 4.43293720E-01),
        'dfm20s_v':  (7.20000000E+06, 4.34969038E-01),
        'dfm40c_v':  (7.20000000E+06, 6.60635710E-01),
        'dfm40s_v':  (7.20000000E+06, 6.55953765E-01),
    }

    core_pos = (1, 1)
    
    # Core mode creation
    core_model = CoreModel(
        name="GE14_mini_core", 
        path_to_yaml_configs=GE14_CORE_YAML, 
        core_description_yaml="GEOM_4x4_mini_CORE_surfaces.yaml"
    )
    
    # 3. Create assembly models : runs geometry analysis for each slice's geometry
    core_model.createAssemblyModels()
    
    # 4. Analyse the core model to build a 3D porosity map
    analyser = CartesianGeometricAnalyser(core_model=core_model, core_i=core_pos[0], core_j=core_pos[1])
    # -------------------------

    # --- TEST CASE MATRIX ---
    powers_to_test = [20.0, 40.0] #kW
    power_profile = 'sine'
    pdrop = 1
    nz = 40
    dfm = 1
    
    for power in powers_to_test:
        # 1. generate case name
        case_name = generate_case_name(pdrop, power, power_profile)
        
        # 2. generate power profile form
        axial_pform = generate_power_profile(power_profile, nz)
        
        # 3. Create a THM1D procedure
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

        assert procedure.thm_module.acool_profile[0] == pytest.approx(0.009821292, abs=1e-5)
        assert procedure.thm_module.acool_profile[-1] == pytest.approx(0.010987695, abs=1e-5)

        assert procedure.thm_module.acool_profile[-1] - procedure.thm_module.acool_profile[0] == pytest.approx(14*np.pi*(0.515e-2)**2, abs=1e-5)
        print(procedure.thm_module.acool_profile)
        assert procedure.thm_module.acool_profile[25] == 0.009039768
        assert procedure.thm_module.acool_profile[26] == 0.010310048
        # 4. Write the c2m file
        refs = REFERENCE_VALS.get(case_name, (None, None))
        c2m_path = procedure.write_to_c2m(OUTPUTS_DIR, case_name,
                                            pinlet_ref=refs[0], epsout_ref=refs[1])