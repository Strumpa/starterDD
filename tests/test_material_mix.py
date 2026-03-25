# testing mix handling classes of the newly implemented startDD package
# Date : 05/02/2026
# R.Guasch

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition


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

def test_natural_elements():
    """
    Test the handling of natural elements in the composition.
    For example, if we have a material that is composed of natural iron, we should be able to specify it as "Fe_NAT" and the code should automatically use the natural isotopic composition of iron.
    """
    Fe_variants_to_test = ["Fe_NAT", "Fe_nat", "Fe"]  # Test different variants of natural element specification
    for variant in Fe_variants_to_test:
        # Define a composition with natural Iron
        natural_iron_composition = Composition(
            material_name="Natural_Iron",
            isotopic_composition={variant: 1.0}  # This indicates that we want to use natural iron
        )

        # Create a material mixture using this composition
        material_mixture = MaterialMixture(
            material_name="Natural_Iron_Mix",
            material_mixture_index=2,
            composition=natural_iron_composition,
            temperature=300.0,
            isdepletable=False
        )

        # Check that the isotopic composition has been correctly set to the natural composition of iron
        assert round(material_mixture.composition.isotopic_composition["Fe54"], 4) == 0.0584  # Natural abundance of Fe-54
        assert round(material_mixture.composition.isotopic_composition["Fe56"], 4) == 0.9175  # Natural abundance of Fe-56
        assert round(material_mixture.composition.isotopic_composition["Fe57"], 5) == 0.02119  # Natural abundance of Fe-57
        assert round(material_mixture.composition.isotopic_composition["Fe58"], 5) == 0.00282  # Natural abundance of Fe-58
        assert variant not in material_mixture.composition.isotopic_composition.keys()  # Check that the original natural element specification is not in the isotopic composition  

    U_variants_to_test = ["U_NAT", "U_nat", "U"]  # Test different variants of natural element specification
    for variant in U_variants_to_test:
        # Define a composition with natural uranium
        natural_uranium_composition = Composition(
            material_name="Natural_Uranium",
            isotopic_composition={variant: 1.0}  # This indicates that we want to use natural uranium
        )

        # Create a material mixture using this composition
        material_mixture = MaterialMixture(
            material_name="Natural_Uranium_Mix",
            material_mixture_index=2,
            composition=natural_uranium_composition,
            temperature=300.0,
            isdepletable=False
        )

        # Check that the isotopic composition has been correctly set to the natural composition of uranium
        assert round(material_mixture.composition.isotopic_composition["U234"], 6) == 0.000054  # Natural abundance of U-234
        assert round(material_mixture.composition.isotopic_composition["U235"], 6) == 0.007204  # Natural abundance of U-235
        assert round(material_mixture.composition.isotopic_composition["U238"], 6) == 0.992742  # Natural abundance of U-238
        assert variant not in material_mixture.composition.isotopic_composition.keys()  # Check that the original natural element specification is not in the isotopic composition

    Zr_variants_to_test = ["Zr_NAT", "Zr_nat", "Zr"]  # Test different variants of natural element specification
    for variant in Zr_variants_to_test:
        # Define a composition with natural zirconium
        natural_zirconium_composition = Composition(
            material_name="Natural_Zirconium",
            isotopic_composition={variant: 1.0}  # This indicates that we want to use natural zirconium
        )
        # Create a material mixture using this composition
        material_mixture = MaterialMixture(
            material_name="Natural_Zirconium_Mix",
            material_mixture_index=3,
            composition=natural_zirconium_composition,
            temperature=300.0,
            isdepletable=False
        )

        # Check that the isotopic composition has been correctly set to the natural composition of zirconium
        assert round(material_mixture.composition.isotopic_composition["Zr90"], 5) == 0.5145  # Natural abundance of Zr-90
        assert round(material_mixture.composition.isotopic_composition["Zr91"], 5) == 0.1122  # Natural abundance of Zr-91
        assert round(material_mixture.composition.isotopic_composition["Zr92"], 5) == 0.1715  # Natural abundance of Zr-92
        assert round(material_mixture.composition.isotopic_composition["Zr94"], 5) == 0.1738  # Natural abundance of Zr-94
        assert round(material_mixture.composition.isotopic_composition["Zr96"], 5) == 0.0280  # Natural abundance of Zr-96
        assert variant not in material_mixture.composition.isotopic_composition.keys()  # Check that the original natural element specification is not in the isotopic composition


    # Check the cladding material from OECD Phase IIIB BWR benchmark, which is composed of natural Zirconium, Iron and Cromium
    cladding_composition = Composition(
        material_name="Cladding",
        isotopic_composition={"Zr_NAT": 4.2982E-02, "Fe_NAT": 1.4838E-04, "Cr_NAT": 7.5891E-05}  # This indicates the isotopic composition of the cladding material
    )
    material_mixture = MaterialMixture(
        material_name="Cladding_Mix",
        material_mixture_index=4,
        composition=cladding_composition,
        temperature=300.0,
        isdepletable=False
    )
    # Check that the isotopic composition has been correctly set to the natural composition of the elements in the cladding
    # ref values from NAT to isotopic composition conversion from OECD Phase IIIB BWR benchmark specification
    reference_values = {
        "Cr50": 3.29746E-06,
        "Cr52": 6.35883E-05,
        "Cr53": 7.21040E-06,
        "Cr54": 1.79482E-06,
        "Fe54": 8.67281E-06,
        "Fe56": 1.36145E-04,
        "Fe57": 3.14417E-06,
        "Fe58": 4.18432E-07,
        "Zr90": 2.21142E-02,
        "Zr91": 4.82258E-03,
        "Zr92": 7.37141E-03,
        "Zr94": 7.47027E-03,
        "Zr96": 1.20350E-03
    }
    for nuclide, expected_value in reference_values.items():
        rel_diff_percent = abs(material_mixture.composition.isotopic_composition[nuclide] - expected_value)*100/expected_value
        assert rel_diff_percent < 1e-3 # Check that the relative difference is less than 0.001%
    assert variant not in material_mixture.composition.isotopic_composition.keys() # Check that the original natural element specification is not in the isotopic composition


if __name__ == "__main__":
    mm = test_material_mixture_instantiation()
    test_natural_elements()
    print("All tests passed.")