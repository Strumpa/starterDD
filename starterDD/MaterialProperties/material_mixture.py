# Class to handle abstract material mixtures
# Author : R. Guasch
# Date : 04/02/2026
# Purpose : Define a material mixture with its index and composition
# uses the composition class which is unique per material, 
# Material mixture allows to have different mixtures with same materials but different properties
# -----------------------------------------------------------------------------------------------
import yaml
HM_isotopes = ["U234", "U235", "U236", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Am243", "Cm244", "Cm245", "Cm246"] 

class Composition:
    def __init__(self, material_name : str, isotopic_composition: dict):
        """
        Composition initialization.
        material name + isotopic composition : dictionnary of isotope name and isotope density in iso/barn*cm
        
        :param material_name (str): material given name eg. "UOX_4.5w%", "B4C" or "moderator" 
        :param isotopic_composition (dict): dictionnary of isotope name and isotope density in iso/barn*cm
        """
        self.material_name = material_name
        self.isotopic_composition = isotopic_composition  # isotopic_composition is a dict of {isotope_name: density}

        # identify if the isotopes are in zaid format or in isotope name format and convert to isotope name format if needed
        if all(zaid.isdigit() for zaid in isotopic_composition.keys()):
            self.zaid_to_isotope()

    def setDepletable(self, depletable: bool):
        """
        Set the depletable flag for the composition.
        
        :param depletable (bool): Flag indicating if the composition is depletable
        """
        self.depletable = depletable

    def zaid_to_isotope(self):
        """
        Convert a ZAID to an isotope name.
        
        :param zaid (str): ZAID string (e.g., "92235" for U-235)
        :return: Isotope name corresponding to the ZAID
        """
        self.isotope_name_composition = {}
        for zaid in self.isotopic_composition.keys():
            #print(f"Converting ZAID {zaid} to isotope name...")
            # Extract Z and A from the ZAID
            Z = int(float(zaid)//1000)  # Z is the part before the last three digits
            A = int(float(zaid) - Z*1000)   # A is the last three digits
            element_symbol = get_element_symbol(Z)
            isotope_name = f"{element_symbol}{A}"
            # update the isotopic composition dict with the isotope name as key instead of ZAID
            self.isotope_name_composition[isotope_name] = self.isotopic_composition[zaid]  # copy the density value to the new key

    def get_isotope_name_composition(self):
        """
        Return the isotope-name-keyed composition dictionary.
        If the original composition used ZAID keys, the converted dict is returned;
        otherwise the original isotopic_composition is returned as-is (assumed to
        already use isotope name keys).
        """
        return getattr(self, 'isotope_name_composition', self.isotopic_composition)

class MaterialMixture:
    def __init__(self, material_name: str, material_mixture_index: int, composition: Composition, temperature:float, isdepletable: bool = False):
        """
        Material mixture initialization.

        In DRAGON5, a material mixture is defined by a material index.
        This class allows for a flexible definition of material mixtures
        with different compositions and properties.

        For self-shielding problems, it is recommended to split fuel materials
        into different material mixtures. This allows to:

        - solve for different self-shielded cross sections per mixture.
        - define a different temperature per mixture.

        :param material_name: material given name e.g. ``"UOX_4.5w%"``, ``"B4C"`` or ``"moderator"``
        :param material_mixture_index: Unique index assigned to the material mixture
        :param composition: Composition object defining the isotopic composition
        :param temperature: Temperature of the material mixture in Kelvin
        :param isdepletable: Flag indicating if the material mixture is depletable
        """
        
        self.material_name = material_name
        self.material_mixture_index = material_mixture_index
        self.composition = composition
        self.temperature = temperature
        self.isdepletable = isdepletable
        self.is_generating = False  # True if this is a "generating" mix (first mix with a given composition)
        self.generating_mix = None  # reference to the generating MaterialMixture if this is a "daughter" mix

    def set_cross_sectional_data(self, xs_data):
        """
        Set the cross-sectional data for the material mixture.
        
        :param xs_data (XSData): XSData object containing cross-sectional data
        """
        self.xs_data = xs_data

    def set_unique_material_mixture_name(self, unique_material_mixture_name: str):
        """
        Set a unique name for the material mixture.
        
        :param unique_material_mixture_name (str): Unique name for the material mixture
        """
        self.unique_material_mixture_name = unique_material_mixture_name

class XSData:
    def __init__(self, mixture_index: int, cross_section_type: str, data: dict):
        """
        Cross-section data initialization.
        
        :param mixture_index (int): Unique index assigned to the material mixture
        :param cross_section_type (str): Type of cross-section data (e.g., "microscopic", "macroscopic")
        :param data (dict): Dictionary containing cross-section data
        """
        self.mixture_index = mixture_index
        self.cross_section_type = cross_section_type
        self.data = data  # data is a dict containing cross-section information

        # minimal is total for each energy group + scattering for self-scattering and transfer matrices
        self.required_keys_MACRO = ["total", "scattering"]

        # check completeness of data
        if not self.check_completeness_MACRO():
            raise ValueError("Cross-section data is incomplete for MACRO type.")


    def check_completeness_MACRO(self):
        """
        Check if the cross-section data contains all required keys.        
        :return (bool): True if all required keys are present, False otherwise
        """
        for key in self.required_keys_MACRO:
            if key not in self.data:
                return False
        return True
    
    def get_cross_sections(self, key: str):
        """
        Get the cross-section data for a specific key.
        
        :param key (str): Key for which to retrieve the cross-section data
        :return: Cross-section data corresponding to the key
        """
        return self.data.get(key, None)
    
    def set_cross_sections(self, key: str, value):
        """
        Set the cross-section data for a specific key.
        
        :param key (str): Key for which to set the cross-section data
        :param value: Cross-section data to be set
        """
        self.data[key] = value

def parse_all_compositions_from_yaml(path_to_yaml_data: str):
    """
    Parse all material compositions from a YAML file and create Composition objects.
    Handles both 'isotopic_composition' and 'composition' + 'density' cases.
    Returns a list of Composition objects.
    """
    with open(path_to_yaml_data, 'r') as file:
        yaml_data = yaml.safe_load(file)

    compositions = []
    mix_list = yaml_data.get('MIX_COMPOSITIONS', [])
    for entry in mix_list:
        name = entry.get('name')
        if 'isotopic_composition' in entry:
            iso_densities = entry['isotopic_composition']
        elif 'composition' in entry and 'density' in entry:
            density = entry['density']
            iso_densities = {iso: prop * density for iso, prop in entry['composition'].items()}
            print(iso_densities)
            # actually update to iso / b*cm
            # assume its water lol
            iso_densities = DensToIsoDens_water(density)
        else:
            raise ValueError(f"Entry for {name} missing isotopic_composition or (composition + density)")
        compositions.append(Composition(name, iso_densities))
    # recover depletable flag for each composition from the yaml file and set it in the composition object
    for comp in compositions:
        for entry in mix_list:
            if entry['name'] == comp.material_name:
                comp.setDepletable(entry.get('depletable', False))
    return compositions


def get_element_symbol(Z: int):
    """
    Get the element symbol for a given atomic number Z.
    
    :param Z (int): Atomic number
    :return: Element symbol corresponding to the atomic number
    """
    # A simple mapping of atomic numbers to element symbols (up to Z=92 for simplicity)
    periodic_table = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
        11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
        19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
        27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
        35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
        43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
        51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce",
        59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
        72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
        81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
        93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr"
    }
    return periodic_table[Z]


def DensToIsoDens_water(density):
    """
    density : density of the material in g/cm3
    isotopic_fractions : dict of isotopic fractions for each isotope in the material (e.g., {"H1": 0.666, "O16": 0.333})
    """
    # Calculation of moderator data
    # AVOGADRO's number
    A = 6.022094E-1 # Normalizing by 10E-24 to obtain isotopic density in # / b*cm

    M_H2O = 15.9994 + 2.0*1.00794
    # compute molar 
    N_MAT = density*A/M_H2O
    N_O = N_MAT
    N_H = 2.0*N_MAT 
    return {"H1": N_H, "O16": N_O}