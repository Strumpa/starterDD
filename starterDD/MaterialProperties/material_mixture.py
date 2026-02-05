# Class to handle abstract material mixtures
# Author : R. Guasch
# Date : 04/02/2026
# Purpose : Define a material mixture with its index and composition
# uses the composition class which is unique per material, 
# Material mixture allows to have different mixtures with same materials but different properties
# -----------------------------------------------------------------------------------------------

HM_isotopes = ["U234", "U235", "U236", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Am243", "Cm244", "Cm245", "Cm246"]


class Composition:
    def __init__(self, material_name : str, isotopic_composition: dict):
        """
        Composition initialization.
        material name + isotopic composition : dictionnary of isotope name and isotope density in iso/barn
        
        :param material_name (str): material given name eg. "UOX_4.5w%", "B4C" or "moderator" 
        :param isotopic_composition (dict): dictionnary of isotope name and isotope density in iso/barn
        """
        self.material_name = material_name
        self.isotopic_composition = isotopic_composition  # isotopic_composition is a dict of {isotope_name: density}

class MaterialMixture:
    def __init__(self, material_name: str, material_mixture_index: int, composition: Composition, temperature:float, isdepletable: bool = False):
        """
        Material mixture initialization.
        In DRAGON5 : a material mixture is defined by a material index. 
        This class allows for a flexible definition of material mixtures with different compositions and properties.
        For self-shielding problems, it is recommended to split fuel materials into different material mixture : 
        This allows to : 
            - solve for different self-shielded cross sections per mixture.
            - define a different temperature per mixture.

        
        :param material_name (str): material given name eg. "UOX_4.5w%", "B4C" or "moderator" 
        :param material_mixture_index (int): Unique index assigned to the material mixture
        :param composition (Composition): Composition object defining the isotopic composition of the material mixture
        :param temperature (float): Temperature of the material mixture in Kelvin
        :param isdepletable (bool): Flag indicating if the material mixture is depletable
        """
        
        self.material_name = material_name
        self.material_mixture_index = material_mixture_index
        self.composition = composition
        self.temperature = temperature
        self.isdepletable = isdepletable

    def set_cross_sectional_data(self, xs_data):
        """
        Set the cross-sectional data for the material mixture.
        
        :param xs_data (XSData): XSData object containing cross-sectional data
        """
        self.xs_data = xs_data

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