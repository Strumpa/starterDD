# testing mix handling classes of the newly implemented startDD package
# Date : 05/02/2026
# R.Guasch

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.InterfaceToDD.dragon_module_calls import MAC
from starterDD.InterfaceToDD.dragon_module_calls import LIB


# instatiate a composition from BWRProgressionProblems benchmark GE-14 fuel assembly

def test_material_mixture_instantiation():
    """
    %  U40
    mat  28  sum rgb  230  56  71
      92235.82c 9.277668758368378E-4
      92234.82c 8.292511730252297E-6
      92236.82c 4.249616794477542E-6
      92238.82c 2.226267248428196E-2
      8016.82c 4.640596297728707E-2
    """

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

    assert material_mixture.material_name == "Fuel_U40_Mix1"
    assert material_mixture.material_mixture_index == 1
    assert material_mixture.composition.material_name == "Fuel_U40"
    assert material_mixture.composition.isotopic_composition["U235"] == 9.277668758368378E-4
    assert material_mixture.temperature == 900.0
    assert material_mixture.isdepletable is True


if __name__ == "__main__":
    test_material_mixture_instantiation()
    print("All tests passed.")