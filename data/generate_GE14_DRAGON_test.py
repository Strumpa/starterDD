# test for GE-14 controlled assembly
# R.Guasch
# Date : 05/02/2026

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition, XSData
from starterDD.InterfaceToDD.dragon_module_calls import MAC

def generate_GE14_DRAGON_test():
    """
    Docstring for generate_GE14_DRAGON_test
    """
    # Parse "glow_data/tdt_data/GE14_ctrl_simplified.dat" to extract material mixture - index mapping

    

    # define material mixtures with MACRO cross sections : 
    material_mixtures = []
    # Example material mixture for UOX 4.5w% fuel at 900K
    composition_45UOX = Composition(
        material_name="UOX_4.5w%",
        isotopic_composition={}
    )
    material_mixture_45UOX = MaterialMixture(
        material_name="Fuel_U45_Mix1",
        material_mixture_index=1,
        composition=composition_45UOX,
        temperature=900.0,
        isdepletable=True
    )
    # set cross section data for 45UOX : 
    xs_data_45UOX = XSData(
        mixture_index=material_mixture_45UOX.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [4.64771211E-01],
              "scattering": [[3.97388905E-01]],
              "nufission": [9.81869251E-02],
              "chi": [1.0]
             }
    )
    # set cross-sectional data for the material mixture
    material_mixture_45UOX.set_cross_sectional_data(xs_data_45UOX)

    material_mixture_45Gd = MaterialMixture(
        material_name="Gd2O3_4.5w%",
        material_mixture_index=2,
        composition=Composition(
            material_name="Gd2O3_4.5w%",
            isotopic_composition={}
        ),
        temperature=900.0,
        isdepletable=False
    )
    # set cross section data for Gd poisonned pin
    # 45Gd :
    xs_data_45Gd = XSData(
        mixture_index=material_mixture_45Gd.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [5.11014581E-01],
              "scattering": [[3.98968190E-01]],
              "nufission": [3.49357314E-02],
              "chi": [1.0]
             }
    )
    # set cross-sectional data for the material mixture
    material_mixture_45Gd.set_cross_sectional_data(xs_data_45Gd)
    # Gap 
    material_mixture_gap = MaterialMixture(
        material_name="GAP",
        material_mixture_index=3,
        composition=Composition(
            material_name="GAP",
            isotopic_composition={}
        ),
        temperature=900.0,
        isdepletable=False
    )
    xs_data_gap = XSData(
        mixture_index=material_mixture_gap.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [3.00460903E-04],
              "scattering": [[3.00460844E-04]]
             }
    )
    material_mixture_gap.set_cross_sectional_data(xs_data_gap)
    # Clad
    """
    TOTAL 0.4029        SCAT 1 1 0.4000
    """
    material_mixture_clad = MaterialMixture(
        material_name="CLAD",
        material_mixture_index=4,
        composition=Composition(
            material_name="CLAD",
            isotopic_composition={}
        ),
        temperature=900.0,
        isdepletable=False
    )
    xs_data_clad = XSData(
        mixture_index=material_mixture_clad.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [0.4029],
              "scattering": [[0.4000]]
             }
    )
    material_mixture_clad.set_cross_sectional_data(xs_data_clad)
    # Moderator
    """
    TOTAL 0.3683        SCAT 1 1 0.3661
    """
    material_mixture_moderator = MaterialMixture(
        material_name="MODERATOR",
        material_mixture_index=5,
        composition=Composition(
            material_name="MODERATOR",
            isotopic_composition={}
        ),
        temperature=900.0,
        isdepletable=False
    )
    xs_data_moderator = XSData(
        mixture_index=material_mixture_moderator.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [0.3683],
              "scattering": [[0.3661]]
             }
    )
    material_mixture_moderator.set_cross_sectional_data(xs_data_moderator)
    # Coolant 
    """
    TOTAL 0.3683        SCAT 1 1 0.3661
    """ 
    material_mixture_coolant = MaterialMixture(
        material_name="COOLANT",
        material_mixture_index=6,
        composition=Composition(
            material_name="COOLANT",
            isotopic_composition={}
        ),
        temperature=900.0,
        isdepletable=False
    )
    xs_data_coolant = XSData(
        mixture_index=material_mixture_coolant.material_mixture_index,
        cross_section_type="macroscopic",
        data={"total": [0.3683],
              "scattering": [[0.3661]]
             }
    )
    material_mixture_coolant.set_cross_sectional_data(xs_data_coolant)

    ### Need to add SS304 and B4C material mixtures too for control blades
    material_mixtures.extend([
        material_mixture_45UOX,
        material_mixture_45Gd,
        material_mixture_gap,
        material_mixture_clad,
        material_mixture_moderator,
        material_mixture_coolant
    ])


    # create MAC module
    macrolib = MAC(macro_lib_name="GE14_MACRO_LIB", create_new=True)
    macrolib.iprint = 1
    macrolib.ngroup = 1
    macrolib.anisotropy_level = 1
    # add material mixtures to MAC module
    for mat_mix in material_mixtures:
        macrolib.add_material_mixture(mat_mix)
    # create MAC module .c2m file (here we just simulate the writing process)
    path_to_procs = "./outputs/"
    proc_name = "GE14_DRAGON_test_macro_lib"
    macrolib.write_to_c2m(path_to_procs, proc_name)


