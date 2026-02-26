# Class to handle abstract material mixtures
# Author : R. Guasch
# Date : 04/02/2026
# Purpose : Define a material mixture with its index and composition
# uses the composition class which is unique per material, 
# Material mixture allows to have different mixtures with same materials but different properties
# -----------------------------------------------------------------------------------------------
import yaml
import warnings

HM_isotopes = ["U234", "U235", "U236", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Am243", "Cm244", "Cm245", "Cm246"]

# ---------------------------------------------------------------------------
#  Default thermal scattering registry
# ---------------------------------------------------------------------------
# Maps an isotope name to the default thermal scattering treatment applied
# when the YAML entry sets ``therm: true``.
#
# * ``dragon_alias``     – evaluation name in the DRAGLIB (LIB: module).
# * ``serpent2_therm_name`` – identifier used in the Serpent2 ``therm`` and
#                            ``moder`` cards (e.g. ``lwtr``).
# * ``serpent2_zaid``    – ZAID attached to the ``moder`` keyword inside the
#                          ``mat`` card.
#
# Users can extend this dict or override individual aliases via
# ``LIB.set_isotope_alias`` / ``Serpent2Model.add_thermal_scattering``.

DEFAULT_THERMAL_SCATTERING = {
    "H1": {
        "dragon_alias":        "H1_H2O",
        "serpent2_therm_name": "lwtr",
        "serpent2_zaid":       "1001",
    },
    "H2": {
        "dragon_alias":        "D2_D2O",
        "serpent2_therm_name": "hwtr",
        "serpent2_zaid":       "1002",
    },
    # Extensible – add graphite, Be, ZrH, etc. as needed.
} 

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

        # Thermal scattering flag & metadata
        # ----------------------------------
        # ``therm`` is ``True`` when this composition requires a bound
        # thermal scattering law (e.g. H in H₂O).  ``therm_data`` holds
        # per-isotope entries resolved from :data:`DEFAULT_THERMAL_SCATTERING`
        # or from explicit YAML data.  It is a list of dicts, each with keys
        # ``isotope``, ``dragon_alias``, ``serpent2_therm_name``,
        # ``serpent2_zaid``.
        self.therm = False
        self.therm_data = []  # populated by setTherm()

        # identify if the isotopes are in zaid format or in isotope name format and convert to isotope name format if needed
        if all(zaid.isdigit() for zaid in isotopic_composition.keys()):
            self.zaid_to_isotope()

    def setDepletable(self, depletable: bool):
        """
        Set the depletable flag for the composition.
        
        :param depletable (bool): Flag indicating if the composition is depletable
        """
        self.depletable = depletable

    def setTherm(self, therm_value):
        """Configure thermal scattering treatment for this composition.

        *therm_value* can be:

        * ``True``  – auto-detect from :data:`DEFAULT_THERMAL_SCATTERING`
          by inspecting which isotopes in the composition have a default
          entry.
        * A **dict** ``{isotope: {dragon_alias, serpent2_therm_name,
          serpent2_zaid}}`` for full manual control.
        * ``False`` / ``None`` – no thermal scattering.

        :param therm_value: bool, dict, or None
        """
        if therm_value is True:
            self.therm = True
            self.therm_data = []
            iso_names = set(self.get_isotope_name_composition().keys())
            for iso, info in DEFAULT_THERMAL_SCATTERING.items():
                if iso in iso_names:
                    self.therm_data.append({
                        "isotope":              iso,
                        "dragon_alias":         info["dragon_alias"],
                        "serpent2_therm_name":  info["serpent2_therm_name"],
                        "serpent2_zaid":        info["serpent2_zaid"],
                    })
        elif isinstance(therm_value, dict):
            self.therm = True
            self.therm_data = []
            for iso, info in therm_value.items():
                self.therm_data.append({
                    "isotope":              iso,
                    "dragon_alias":         info.get("dragon_alias", iso),
                    "serpent2_therm_name":  info.get("serpent2_therm_name"),
                    "serpent2_zaid":        info.get("serpent2_zaid"),
                })
        else:
            self.therm = False
            self.therm_data = []

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

    Supported YAML entry formats:

    1. **Direct isotopic number densities** (iso/barn·cm)::

        isotopic_composition:
          "92235": 3.758E-4
          "8016":  4.639E-2

    2. **Total atomic density** (atoms/barn·cm) + atomic fractions::

        atomic_density: 7.389E-2
        composition:
          "1001": 0.6667
          "8016": 0.3333

    3. **Mass density** (g/cm³) + mass fractions::

        mass_density: 0.998
        composition:
          "1001": 0.1119
          "8016": 0.8881

    :param path_to_yaml_data: Path to the YAML file containing mix compositions
    :return: List of :class:`Composition` objects
    """
    with open(path_to_yaml_data, 'r') as file:
        yaml_data = yaml.safe_load(file)

    compositions = []
    mix_list = yaml_data.get('MIX_COMPOSITIONS', [])
    for entry in mix_list:
        name = entry.get('name')
        if 'isotopic_composition' in entry:
            # Case 1: direct number densities already in iso/barn·cm
            iso_densities = entry['isotopic_composition']

        elif 'composition' in entry and 'atomic_density' in entry:
            # Case 2: total atomic density + atomic (number) fractions
            iso_densities = fractions_to_iso_densities(
                composition=entry['composition'],
                density_type="atomic_density",
                density_value=entry['atomic_density'],
            )

        elif 'composition' in entry and 'mass_density' in entry:
            # Case 3: mass density (g/cm³) + mass (weight) fractions
            iso_densities = fractions_to_iso_densities(
                composition=entry['composition'],
                density_type="mass_density",
                density_value=entry['mass_density'],
            )

        else:
            raise ValueError(
                f"Entry for '{name}' must provide one of: "
                f"'isotopic_composition', "
                f"'atomic_density' + 'composition', or "
                f"'mass_density' + 'composition'."
            )
        compositions.append(Composition(name, iso_densities))

    # recover depletable flag and therm flag for each composition from the yaml file
    for comp in compositions:
        for entry in mix_list:
            if entry['name'] == comp.material_name:
                comp.setDepletable(entry.get('depletable', False))
                comp.setTherm(entry.get('therm', False))
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


# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------
AVOGADRO = 6.022094e23       # Avogadro's number [1/mol]
CM2_TO_BARN = 1e24           # barn -> cm² conversion factor


def get_isotope_atomic_mass(zaid: str) -> float:
    """Return the atomic mass (in g/mol) for the isotope identified by *zaid*.

    For now the mass number *A* is used as an approximation.
    Replace the body of this function with a lookup table of evaluated
    (AUDI) atomic masses for higher accuracy.

    :param zaid: ZAID string, e.g. ``"92235"`` for U-235 or ``"1001"`` for H-1
    :return: Atomic mass in g/mol (≈ mass number *A*)
    """
    Z = int(float(zaid) // 1000)
    A = int(float(zaid) - Z * 1000)
    if A == 0:
        raise ValueError(
            f"Cannot determine mass number from ZAID '{zaid}'. "
            "Natural-element ZAIDs (A=0) are not supported; "
            "please specify individual isotopes."
        )
    return float(A)


def fractions_to_iso_densities(
    composition: dict,
    density_type: str,
    density_value: float,
    fraction_tolerance: float = 1e-3,
) -> dict:
    """Convert fraction-based composition + bulk density to isotopic number densities.

    Supports two density modes:

    * ``"atomic_density"`` – *density_value* is the total atomic density
      :math:`N_{tot}` in atoms/barn·cm and *composition* values are **atomic
      (number) fractions** :math:`f_i` with :math:`\\sum f_i = 1`.

      .. math:: N_i = f_i \\cdot N_{tot}

    * ``"mass_density"`` – *density_value* is the mass density :math:`\\rho`
      in g/cm³ and *composition* values are **mass (weight) fractions**
      :math:`w_i` with :math:`\\sum w_i = 1`.

      .. math:: N_i = \\frac{\\rho \\, N_A \\, w_i}{A_i \\cdot 10^{24}}

    :param composition: ``{zaid_string: fraction}`` dictionary.
    :param density_type: ``"atomic_density"`` or ``"mass_density"``.
    :param density_value: Bulk density value (see above for units).
    :param fraction_tolerance: Allowed deviation of the fraction sum from 1.0
        (default 1e-3).
    :return: ``{zaid_string: number_density}`` in atoms/barn·cm.
    :raises ValueError: On unknown *density_type* or if fractions do not sum
        to ~1.
    """
    # --- validate fractions ------------------------------------------------
    frac_sum = sum(composition.values())
    if abs(frac_sum - 1.0) > fraction_tolerance:
        raise ValueError(
            f"Composition fractions sum to {frac_sum:.6g}, expected ~1.0 "
            f"(tolerance={fraction_tolerance})."
        )

    iso_densities: dict = {}

    if density_type == "atomic_density":
        # N_i = f_i * N_tot
        N_tot = density_value  # already in atoms/barn·cm
        for zaid, frac in composition.items():
            iso_densities[zaid] = frac * N_tot

    elif density_type == "mass_density":
        # N_i = rho * N_A * w_i / (A_i * 1e24)
        rho = density_value  # g/cm³
        for zaid, w_i in composition.items():
            A_i = get_isotope_atomic_mass(zaid)
            iso_densities[zaid] = (rho * AVOGADRO * w_i) / (A_i * CM2_TO_BARN)

    else:
        raise ValueError(
            f"Unknown density_type '{density_type}'. "
            f"Use 'atomic_density' or 'mass_density'."
        )

    return iso_densities


def DensToIsoDens_water(density):
    """Convert water mass density to isotopic number densities.

    .. deprecated::
        Use :func:`fractions_to_iso_densities` with ``density_type='mass_density'``
        and explicit mass fractions instead.

    :param density: Water mass density in g/cm³.
    :return: ``{"H1": N_H, "O16": N_O}`` in atoms/barn·cm.
    """
    warnings.warn(
        "DensToIsoDens_water is deprecated. "
        "Use fractions_to_iso_densities(composition, 'mass_density', density) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    M_H2O = 15.9994 + 2.0 * 1.00794
    N_MAT = density * AVOGADRO / (M_H2O * CM2_TO_BARN)
    N_O = N_MAT
    N_H = 2.0 * N_MAT
    return {"H1": N_H, "O16": N_O}