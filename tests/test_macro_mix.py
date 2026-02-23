# testing mix handling classes of the newly implemented startDD package
# Date : 05/02/2026
# R.Guasch

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition, XSData
from starterDD.InterfaceToDD.dragon_module_calls import MAC
from starterDD.InterfaceToDD.dragon_module_calls import LIB
from test_material_mix import test_material_mixture_instantiation


def test_macro_library_creation():
    # instatiate a composition from BWRProgressionProblems benchmark GE-14 fuel assembly
    U40_components = {
        "U235": 9.277668758368378E-4,
        "U234": 8.292511730252297E-6,
        "U236": 4.249616794477542E-6,
        "U238": 2.226267248428196E-2,
        "O16": 4.640596297728707E-2
    }
    composition = Composition(material_name="Fuel_U40", isotopic_composition=U40_components)
    material_mixture = MaterialMixture(
        material_name="Fuel_U40_Mix1",
        material_mixture_index=1,
        composition=composition,
        temperature=900.0,
        isdepletable=True
    )
    # material cross section data (dummy data for testing) from DRAGON output of UOX 4.5w% fuel at 900K, 1 energy group from ATRIUM-10 cell
    xs_data_24UOX = XSData(
        mixture_index=material_mixture.material_mixture_index,
        cross_section_type="macroscopic",
        data={
            "total": [4.64771211E-01],
            "scattering": [[3.97388905E-01]],
            "nufission": [9.81869251E-02],
            "chi": [1.0]
        }
    )
    
    # create MAC module
    macrolib = MAC(macro_lib_name="MACROLIB", create_new=True)
    macrolib.iprint = 1
    macrolib.ngroup = 1
    macrolib.anisotropy_level = 1

    # add material mixture to MAC module
    macrolib.add_material_mixture(material_mixture)
    # set cross-sectional data for the material mixture
    material_mixture.set_cross_sectional_data(xs_data_24UOX)

    # create MAC module .c2m file (here we just simulate the writing process)
    path_to_procs = "./outputs/"
    proc_name = "test_macro_lib"
    macrolib.write_to_c2m(path_to_procs, proc_name)

    # Here we would normally write to a .c2m file, but for testing purposes, we will just check the attributes
    assert macrolib.macro_lib_name == "MACROLIB"
    assert macrolib.count_mixtures == 1
    assert macrolib.material_mixtures[0].material_name == "Fuel_U40_Mix1"

if __name__ == "__main__":
    
    ## test MACROscopic library creation
    test_macro_library_creation()
    print("All tests passed.")