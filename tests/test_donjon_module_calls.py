import os
import math
import numpy as np
import pytest

from conftest import GE14_CORE_YAML, OUTPUTS_DIR

from starterDD.GeometryAnalysis.cartesian_geometry_analysis import CartesianGeometricAnalyser
from starterDD.InterfaceToDD.donjon_module_calls import DonjonTHM1DProcedure, INIT, NEUTRONICS
from starterDD.DDModel.DonjonModel import CoreModel
from starterDD.DDModel.DonjonCalculationScheme import DonjonCalculationScheme, ModelInitialisation, NeutronicsSolve

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


@pytest.fixture
def single_channel_model_initialiser():
    GE14_core_description_yaml = "GEOM_single_assembly_CORE.yaml"
    core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    single_channel_model_initialiser = ModelInitialisation(core_model=core_model, 
                                        radial_homogenization_strategy="by_assembly",
                                        boundary_conditions_dict={"x": "reflective",
                                                                  "y": "reflective",
                                                                  "z": "void"},
                                        number_of_axial_materials_per_slice=[10, 5],
                                        number_of_energy_groups=2,
                                        interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                 ("TCool", "local", "coolant_temperature"),
                                                                 ("DCool", "local", "coolant_density")])
    
    single_channel_model_initialiser.resolve_mesh()

    single_channel_model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 222.0595], compo_name="CPODOM", dir_name="EDI_2G")
    single_channel_model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[222.0595, 347.1291], compo_name="CPOVAN", dir_name="EDI_2G")
    single_channel_model_initialiser.set_reactor_power(870e6)
    single_channel_model_initialiser.set_initial_axial_power_form("uniform")
    single_channel_model_initialiser.set_initial_parameters({"TFuel": 900, "TCool": 600.0, "DCool": 0.73669})
    single_channel_model_initialiser.set_fuel_mass(total_fuel_mass = 0.5*16)
    
    return single_channel_model_initialiser

@pytest.fixture
def minicore_model_initialiser():
    GE14_core_description_yaml = "GEOM_mini_CORE_reflector.yaml"
    core_model = CoreModel(name="GE14_single_assembly_core", path_to_yaml_configs=GE14_CORE_YAML, core_description_yaml=GE14_core_description_yaml)
    minicore_model_initialiser = ModelInitialisation(core_model=core_model, 
                                        radial_homogenization_strategy="by_assembly",
                                        boundary_conditions_dict={"x": "reflective",
                                                                  "y": "reflective",
                                                                  "z": "void"},
                                        number_of_axial_materials_per_slice=[10, 5],
                                        number_of_energy_groups=2,
                                        interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                 ("TCool", "local", "coolant_temperature"),
                                                                 ("DCool", "local", "coolant_density")])
    
    minicore_model_initialiser.resolve_mesh()

    minicore_model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 222.0595], compo_name="CPODOM", dir_name="EDI_2G")
    minicore_model_initialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[222.0595, 347.1291], compo_name="CPOVAN", dir_name="EDI_2G")
    minicore_model_initialiser.set_reactor_power(870e6)
    minicore_model_initialiser.set_initial_axial_power_form("uniform")
    minicore_model_initialiser.set_initial_parameters({"TFuel": 900, "TCool": 600.0, "DCool": 0.73669})
    minicore_model_initialiser.set_fuel_mass(total_fuel_mass = 0.5*16)
    
    return minicore_model_initialiser

@pytest.fixture
def neutronics_step(single_channel_model_initialiser):
    neutronics_step = NeutronicsSolve(single_channel_model_initialiser, "diffusion", "linear", ["MCFD", 1])
    nz = len(single_channel_model_initialiser.meshz) - 1
    neutronics_step.set_interpolation_parameter("TFuel", [1200.0]*nz)
    neutronics_step.set_interpolation_parameter("TCool", [559.0]*nz)
    neutronics_step.set_interpolation_parameter("DCool", [0.600]*nz)
    neutronics_step.resolve_cpo_dir_to_mix()
    return neutronics_step


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


def test_DONJON_THM_generator_and_porosity_calculation(minicore_model_initialiser):

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
            model_initialiser=minicore_model_initialiser,
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


def test_4x4_minicore_surfaces_definition(minicore_model_initialiser):

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
            model_initialiser=minicore_model_initialiser,
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

        assert procedure.thm_module.acool_profile[0] == pytest.approx(0.009741326, abs=1e-5)
        assert procedure.thm_module.acool_profile[-1] == pytest.approx(0.010907729, abs=1e-5)

        assert procedure.thm_module.acool_profile[-1] - procedure.thm_module.acool_profile[0] == pytest.approx(14*np.pi*(0.515e-2)**2, abs=1e-5)
        assert procedure.thm_module.acool_profile[25] == pytest.approx(0.008959802, abs=1e-5)
        assert procedure.thm_module.acool_profile[26] == pytest.approx(0.010230082, abs=1e-5)
        # 4. Write the c2m file
        refs = REFERENCE_VALS.get(case_name, (None, None))
        c2m_path = procedure.write_to_c2m(OUTPUTS_DIR, case_name,
                                            pinlet_ref=refs[0], epsout_ref=refs[1])


def test_donjon_single_channel_model_INIT_procedure(single_channel_model_initialiser, neutronics_step):

    
    scheme = DonjonCalculationScheme(name="test_scheme")
    scheme.add_initialisation_step(single_channel_model_initialiser)
    scheme.add_neutronics_step(neutronics_step)

    init_proc = INIT(scheme)
    #init_proc.build_GEO_call()
    assert init_proc.scheme.initialisation_step.nfuel == 15
    #init_proc.build_MATEX_call()
    #init_proc.build_RESINI_call()
    assert init_proc.scheme.initialisation_step.nfuel == 15

    neutronics_proc = NEUTRONICS(scheme)
    neutronics_proc.build_RESINI_call()
    neutronics_proc.build_NCR_call()
    init_proc.write_to_c2m("tests/outputs", "IniDonjon_single_channel")


def test_donjon_minicore_model_INIT_procedure(minicore_model_initialiser):

    
    scheme = DonjonCalculationScheme(name="test_scheme")
    scheme.add_initialisation_step(minicore_model_initialiser)

    init_proc = INIT(scheme)
    #init_proc.build_GEO_call()
    assert init_proc.scheme.initialisation_step.nfuel == 16*15
    #init_proc.build_MATEX_call()
    assert init_proc.scheme.initialisation_step.nfuel == 16*15
    #init_proc.build_RESINI_call()
    assert init_proc.scheme.initialisation_step.nfuel == 16*15
    init_proc.write_to_c2m("tests/outputs", "IniDonjon_minicore")
