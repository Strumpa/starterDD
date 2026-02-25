## Collection of classes to write Serpent2 cards from a starterDD AssemblyModel
# Goal: generate a complete Serpent2 input deck equivalent to a Dragon5 model
#
# Author : R. Guasch
# Date : 02/24/2026

import re
import string
import warnings
from typing import Optional, Union

# ── Supported energy meshes ──

# 2g collapsed at 0.625 eV

ene2g = [1.1E-11, 6.25E-7, 1.9640E+1]


# 26g condensed from SHEM295 : used in 2nd level DRAGON flux calculations
ene26g = [
    1.10002700E-10,  1.04505001E-08,  3.43997590E-08,  7.64968619E-08,  1.37999400E-07,
    2.31192306E-07,  3.52993488E-07,  6.24998689E-07,  8.80024374E-07,  9.63959813E-07,
    1.00903499E-06,  1.10395002E-06,  1.14796901E-06,  1.25093901E-06,  4.00000000E-05,
    5.29895287E-05,  4.19093597E-04,  1.81183301E-03,  9.11880762E-03,  2.49990801E-02,
    6.73793828E-02,  1.95066203E-01,  4.94001812E-01,  1.33694100E+00,  2.23129900E+00,
    4.96584700E+00,  1.96403000E+01
]

#295g SHEM295 energy mesh used in light water reactor calculations in DRAGON5
ene295g = [
    1.10003E-10,  2.49990E-09,  4.55602E-09,  7.14526E-09,  1.04505E-08,  1.48300E-08,  2.00104E-08,  2.49394E-08,  2.92989E-08,  3.43998E-08, 
    4.02999E-08,  4.73019E-08,  5.54981E-08,  6.51994E-08,  7.64969E-08,  8.97968E-08,  1.04298E-07,  1.19995E-07,  1.37999E-07,  1.61895E-07,  
    1.90005E-07,  2.09610E-07,  2.31192E-07,  2.54997E-07,  2.79989E-07,  3.05012E-07,  3.25008E-07,  3.52993E-07,  3.90001E-07,  4.31579E-07,  
    4.75017E-07,  5.20011E-07,  5.54990E-07,  5.94993E-07,  6.24999E-07,  7.19999E-07,  8.20037E-07,  8.80024E-07,  9.19978E-07,  9.44022E-07,  
    9.63960E-07,  9.81959E-07,  9.96500E-07,  1.00903E-06,  1.02101E-06,  1.03499E-06,  1.07799E-06,  1.09198E-06,  1.10395E-06,  1.11605E-06,  
    1.12997E-06,  1.14797E-06,  1.16999E-06,  1.21397E-06,  1.25094E-06,  1.29304E-06,  1.33095E-06,  1.38098E-06,  1.41001E-06,  1.44397E-06,  
    1.51998E-06,  1.58803E-06,  1.66895E-06,  1.77997E-06,  1.90008E-06,  1.98992E-06,  2.07010E-06,  2.15695E-06,  2.21709E-06,  2.27299E-06,  
    2.33006E-06,  2.46994E-06,  2.55000E-06,  2.59009E-06,  2.62005E-06,  2.64004E-06,  2.70011E-06,  2.71990E-06,  2.74092E-06,  2.77512E-06,  
    2.88405E-06,  3.14211E-06,  3.54307E-06,  3.71209E-06,  3.88217E-06,  4.00000E-06,  4.21983E-06,  4.30981E-06,  4.40216E-06,  4.63249E-06,  
    4.95846E-06,  5.31798E-06,  5.61867E-06,  5.91857E-06,  6.22824E-06,  6.50841E-06,  6.82843E-06,  7.21451E-06,  7.62242E-06,  7.96530E-06,  
    8.41566E-06,  8.85599E-06,  9.31005E-06,  9.78738E-06,  1.02379E-05,  1.07091E-05,  1.12020E-05,  1.17411E-05,  1.22325E-05,  1.28854E-05,  
    1.35732E-05,  1.42121E-05,  1.49707E-05,  1.57382E-05,  1.65617E-05,  1.73761E-05,  1.83035E-05,  1.90315E-05,  2.01075E-05,  2.05343E-05,  
    2.12232E-05,  2.20011E-05,  2.26712E-05,  2.46578E-05,  2.78852E-05,  3.16929E-05,  3.30855E-05,  3.45392E-05,  3.56980E-05,  3.60568E-05,  
    3.64191E-05,  3.68588E-05,  3.73038E-05,  3.77919E-05,  3.87874E-05,  3.97295E-05,  4.12270E-05,  4.21441E-05,  4.31246E-05,  4.41721E-05,  
    4.52904E-05,  4.62053E-05,  4.75173E-05,  4.92591E-05,  5.17847E-05,  5.29895E-05,  5.40600E-05,  5.70595E-05,  5.99250E-05,  6.23083E-05,  
    6.36306E-05,  6.45922E-05,  6.50460E-05,  6.55029E-05,  6.58312E-05,  6.61612E-05,  6.64929E-05,  6.68261E-05,  6.90682E-05,  7.18869E-05,  
    7.35595E-05,  7.63322E-05,  7.93679E-05,  8.39393E-05,  8.87740E-05,  9.33256E-05,  9.73287E-05,  1.00594E-04,  1.01098E-04,  1.01605E-04,  
    1.02115E-04,  1.03038E-04,  1.05646E-04,  1.10288E-04,  1.12854E-04,  1.15480E-04,  1.16524E-04,  1.17577E-04,  1.20554E-04,  1.26229E-04,  
    1.32701E-04,  1.39504E-04,  1.46657E-04,  1.54176E-04,  1.63056E-04,  1.67519E-04,  1.75229E-04,  1.83294E-04,  1.84952E-04,  1.86251E-04,  
    1.87559E-04,  1.88877E-04,  1.90204E-04,  1.93078E-04,  1.95996E-04,  2.00958E-04,  2.12108E-04,  2.24325E-04,  2.35590E-04,  2.41796E-04,  
    2.56748E-04,  2.68297E-04,  2.76468E-04,  2.84888E-04,  2.88327E-04,  2.95922E-04,  3.19927E-04,  3.35323E-04,  3.53575E-04,  3.71703E-04,  
    3.90760E-04,  4.19094E-04,  4.53999E-04,  5.01746E-04,  5.39204E-04,  5.77146E-04,  5.92941E-04,  6.00099E-04,  6.12834E-04,  6.46837E-04,  
    6.77286E-04,  7.48517E-04,  8.32218E-04,  9.09681E-04,  9.82494E-04,  1.06432E-03,  1.13467E-03,  1.34358E-03,  1.58620E-03,  1.81183E-03,  
    2.08410E-03,  2.39729E-03,  2.70024E-03,  2.99618E-03,  3.48107E-03,  4.09735E-03,  5.00451E-03,  6.11252E-03,  7.46585E-03,  9.11881E-03,  
    1.11377E-02,  1.36037E-02,  1.48997E-02,  1.62005E-02,  1.85847E-02,  2.26994E-02,  2.49991E-02,  2.61001E-02,  2.73944E-02,  2.92810E-02,  
    3.34596E-02,  3.69786E-02,  4.08677E-02,  4.99159E-02,  5.51656E-02,  6.73794E-02,  8.22974E-02,  9.46645E-02,  1.15623E-01,  1.22773E-01,  
    1.40098E-01,  1.65065E-01,  1.95066E-01,  2.30060E-01,  2.67826E-01,  3.20646E-01,  3.83883E-01,  4.12501E-01,  4.56021E-01,  4.94002E-01,  
    5.78442E-01,  7.06511E-01,  8.60006E-01,  9.51119E-01,  1.05115E+00,  1.16205E+00,  1.28696E+00,  1.33694E+00,  1.40577E+00,  1.63654E+00,  
    1.90139E+00,  2.23130E+00,  2.72531E+00,  3.32871E+00,  4.06569E+00,  4.96585E+00,  6.06530E+00,  6.70319E+00,  7.40817E+00,  8.18730E+00,  
    9.04836E+00,  9.99999E+00,  1.16183E+01,  1.38403E+01,  1.49182E+01,  1.96403E+01 
] 

# ── Supported energy mesh registry ──
# Maps user-facing names to boundary lists so that detectors and energy
# grids can be built by name without manually supplying boundaries.

SUPPORTED_ENERGY_MESHES = {
    "2g":   ene2g,
    "26g":  ene26g,
    "295g": ene295g,
}

# ── Reaction MT number mappings ──

REACTION_TO_MT_NUMBER = {
    'total': 1,
    'elastic scattering': 2,
    'fission': 18,
    'n,2n': 16,
    'n,3n': 17,
    'n,np': 28,
    'n,4n': 37,
    'n,gamma': 102,
    'n,proton': 103,
    'n,deutron': 104,
    'n,triton': 105,
    'n,alpha': 107,
    'n,2alpha': 108,
    'absorption': 27,
    'disappearance': 101,
}

DRAGON_REAC_TO_REACTION_NAME = {
    'NTOT0': 'total',
    'SIGS00': 'neutronic scattering',
    'NFTOT': 'fission',
    'NP': 'n,proton',
    'ND': 'n,deuton',
    'NT': 'n,triton',
    'NA': 'n,alpha',
    'N2A': 'n,2alpha',
    'N2N': 'n,2n',
    'N3N': 'n,3n',
    'N4N': 'n,4n',
    'NG': 'n,gamma',
    'NNP': 'n,np',
}

# ── Element symbol to atomic number mapping ──

ELEMENT_TO_Z = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
    'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
    'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
    'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
    'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
    'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
    'Fm': 100,
}

# Reverse mapping for ZAID to element symbol
Z_TO_ELEMENT = {v: k for k, v in ELEMENT_TO_Z.items()}

# ── Temperature to cross-section suffix mapping ──
# These correspond to common ENDF/B and JEFF library temperature points.
# Format: temperature_K -> suffix string
# The suffix convention is .XXc where XX encodes the temperature.
# This mapping should be adapted to the specific nuclear data library used.

TEMPERATURE_TO_XS_SUFFIX = {
    293.6: '.03c',
    300.0: '.03c',
    600.0: '.05c',
    900.0: '.09c',
    1200.0: '.12c',
    1500.0: '.15c',
}

def get_xs_suffix(temperature: float, suffix_map: dict = None) -> str:
    """Get the closest cross-section suffix for a given temperature.
    
    Args:
        temperature: Temperature in Kelvin.
        suffix_map: Optional custom mapping of temperature -> suffix.
                   If None, uses TEMPERATURE_TO_XS_SUFFIX.
    
    Returns:
        Cross-section suffix string (e.g., '.81c').
    """
    if suffix_map is None:
        suffix_map = TEMPERATURE_TO_XS_SUFFIX
    
    # Exact match
    if temperature in suffix_map:
        return suffix_map[temperature]
    
    # Find closest temperature
    closest_temp = min(suffix_map.keys(), key=lambda t: abs(t - temperature))
    return suffix_map[closest_temp]


def _reaction_name_to_mt(reaction: Union[str, int]) -> int:
    """Convert a reaction name to its MT number.
    
    Args:
        reaction: Either a reaction name string (e.g., 'fission', 'absorption')
                 or an MT number (int). If already an int, returns it unchanged.
    
    Returns:
        MT number (int).
    
    Raises:
        ValueError: If reaction name is not found in REACTION_TO_MT_NUMBER.
    """
    if isinstance(reaction, int):
        return reaction
    
    if reaction in REACTION_TO_MT_NUMBER:
        return REACTION_TO_MT_NUMBER[reaction]
    
    # Try case-insensitive match
    reaction_lower = reaction.lower()
    for name, mt in REACTION_TO_MT_NUMBER.items():
        if name.lower() == reaction_lower:
            return mt
    
    raise ValueError(
        f"Unknown reaction name: '{reaction}'. "
        f"Valid names: {list(REACTION_TO_MT_NUMBER.keys())}"
    )


def parse_isotope_name(name: str) -> tuple:
    """Parse a Dragon-style isotope name to (element_symbol, mass_number, Z, A).
    
    Handles formats like: 'U235', 'O16', 'Gd155', 'Pu239', 'B10', 'H1'
    Also handles ZAID-style names like '92235', '8016'
    
    Args:
        name: Isotope name string.
    
    Returns:
        Tuple of (element_symbol, mass_number, Z, A) or raises ValueError.
    """
    # Try Dragon-style name: ElementMassNumber (e.g., U235, Gd155)
    match = re.match(r'^([A-Z][a-z]?)(\d+)$', name)
    if match:
        element = match.group(1)
        mass_number = int(match.group(2))
        if element in ELEMENT_TO_Z:
            Z = ELEMENT_TO_Z[element]
            return (element, mass_number, Z, mass_number)
    
    # Try pure ZAID format (e.g., 92235 -> U235)
    try:
        zaid = int(name)
        Z = zaid // 1000
        A = zaid % 1000
        if Z in Z_TO_ELEMENT:
            element = Z_TO_ELEMENT[Z]
            return (element, A, Z, A)
    except ValueError:
        pass
    
    raise ValueError(f"Cannot parse isotope name: '{name}'")


def isotope_name_to_zaid_str(name: str) -> str:
    """Convert a Dragon-style isotope name to a ZAID string.
    
    e.g., 'U235' -> '92235', 'O16' -> '8016', 'Gd155' -> '64155'
    
    Args:
        name: Isotope name string (Dragon convention).
    
    Returns:
        ZAID string (e.g., '92235').
    """
    _, _, Z, A = parse_isotope_name(name)
    return f"{Z}{A:03d}" if Z < 100 else f"{Z}{A:03d}"


# ── Helper to safely extract pin_pitch from an assembly model ──

def _get_pin_pitch(assembly_model) -> float:
    """Extract pin pitch from a CartesianAssemblyModel.
    
    The assembly stores the pitch in ``pin_geometry_dict['pin_pitch']``
    rather than as a top-level attribute.  This helper handles both
    conventions so the Serpent2 interface is resilient to future changes.
    
    Args:
        assembly_model: A starterDD CartesianAssemblyModel.
    
    Returns:
        Pin pitch in cm.
    
    Raises:
        RuntimeError: If the pitch cannot be determined.
    """
    # Direct attribute (in case a future version exposes it)
    if hasattr(assembly_model, 'pin_pitch'):
        return assembly_model.pin_pitch
    # Standard location in the geometry dict
    if hasattr(assembly_model, 'pin_geometry_dict') and assembly_model.pin_geometry_dict:
        pitch = assembly_model.pin_geometry_dict.get('pin_pitch', None)
        if pitch is not None:
            return pitch
    raise RuntimeError(
        "Cannot determine pin pitch from the assembly model.  "
        "Ensure 'pin_pitch' is set in the geometry description YAML."
    )


# ═══════════════════════════════════════════════════════════════
#  S2_Material
# ═══════════════════════════════════════════════════════════════

class S2_Material:
    """Represents a single Serpent2 `mat` card.
    
    Can be created from a starterDD MaterialMixture or from raw data.
    """
    
    def __init__(self, name: str, composition: dict, temperature: float,
                 is_burnable: bool = False, density_type: str = 'sum',
                 total_density: float = None, rgb: tuple = None,
                 moder: tuple = None, xs_suffix: str = None):
        """
        Args:
            name: Material name (Serpent2 identifier).
            composition: Dict of {zaid_string: number_density} pairs.
                        zaid_string should be like '92235' (without suffix).
                        number_density in atoms/barn-cm.
            temperature: Temperature in Kelvin.
            is_burnable: Whether this material is depletable.
            density_type: 'sum' for sum of isotopic densities, 'total' for total given density.
            total_density: Total atom density (used when density_type='total').
            rgb: Optional RGB color tuple for plotting.
            moder: Optional tuple (thermal_name, nuclide_zaid) for thermal scattering.
            xs_suffix: Override for cross-section temperature suffix.
        """
        self.name = name
        self.composition = composition  # {zaid_str: number_density}
        self.temperature = temperature
        self.is_burnable = is_burnable
        self.density_type = density_type
        self.total_density = total_density
        self.rgb = rgb
        self.moder = moder
        self._xs_suffix = xs_suffix
    
    @property
    def xs_suffix(self) -> str:
        if self._xs_suffix:
            return self._xs_suffix
        return get_xs_suffix(self.temperature)
    
    @xs_suffix.setter
    def xs_suffix(self, value: str):
        self._xs_suffix = value
    
    @classmethod
    def from_material_mixture(cls, material_mixture, name_override: str = None,
                              xs_suffix: str = None, rgb: tuple = None):
        """Create an S2_Material from a starterDD MaterialMixture object.
        
        Args:
            material_mixture: A starterDD MaterialMixture object with attributes:
                - material_name (str): base material name
                - unique_material_mixture_name (str): full zone+pin name
                - composition (Composition object with get_isotope_name_composition())
                - temperature (float)
                - isdepletable (bool)
            name_override: Explicit Serpent2 material name.  If None, the
                ``unique_material_mixture_name`` is used (preserving the
                Dragon naming scheme).
            xs_suffix: Override for cross-section suffix.
            rgb: Optional RGB color.
        
        Returns:
            S2_Material instance.
        """
        composition = {}
        # Use the accessor that handles both ZAID and isotope-name keyed dicts
        iso_comp = material_mixture.composition.get_isotope_name_composition()
        for isotope_name, density in iso_comp.items():
            try:
                zaid = isotope_name_to_zaid_str(isotope_name)
                composition[zaid] = density
            except ValueError:
                print(f"Warning: Could not convert isotope '{isotope_name}' to ZAID, skipping.")
        
        # Determine the Serpent2 material name
        if name_override is not None:
            s2_name = name_override
        else:
            s2_name = getattr(material_mixture, 'unique_material_mixture_name',
                              material_mixture.material_name)
        
        return cls(
            name=s2_name,
            composition=composition,
            temperature=material_mixture.temperature,
            is_burnable=material_mixture.isdepletable,
            xs_suffix=xs_suffix,
            rgb=rgb,
        )
    
    @classmethod
    def from_raw(cls, name: str, isotope_densities: dict, temperature: float,
                 is_burnable: bool = False, xs_suffix: str = None,
                 rgb: tuple = None, moder: tuple = None):
        """Create from a dict of {isotope_name_or_zaid: number_density}.
        
        Accepts both Dragon-style names ('U235') and ZAID strings ('92235').
        
        Args:
            name: Material name.
            isotope_densities: Dict of isotope densities.
            temperature: Temperature in Kelvin.
            is_burnable: Whether depletable.
            xs_suffix: Override cross-section suffix.
            rgb: RGB color tuple.
            moder: Thermal scattering law tuple (therm_name, nuclide_zaid).
        
        Returns:
            S2_Material instance.
        """
        composition = {}
        for iso, dens in isotope_densities.items():
            # Check if already a ZAID string (digits only)
            if re.match(r'^\d+$', iso):
                composition[iso] = dens
            else:
                try:
                    zaid = isotope_name_to_zaid_str(iso)
                    composition[zaid] = dens
                except ValueError:
                    print(f"Warning: Could not convert isotope '{iso}' to ZAID, skipping.")
        
        return cls(
            name=name,
            composition=composition,
            temperature=temperature,
            is_burnable=is_burnable,
            xs_suffix=xs_suffix,
            rgb=rgb,
            moder=moder,
        )
    
    def format_card(self) -> str:
        """Format the complete Serpent2 `mat` card as a string.
        
        Returns:
            Formatted mat card string.
        """
        parts = [f"mat {self.name}"]
        
        if self.density_type == 'sum':
            parts.append("sum")
        elif self.total_density is not None:
            parts.append(f"{self.total_density:.8E}")
        
        parts.append(f"tmp {self.temperature:.1f}")
        
        if self.is_burnable:
            parts.append("burn 1")
        
        if self.rgb:
            parts.append(f"rgb {self.rgb[0]} {self.rgb[1]} {self.rgb[2]}")
        
        if self.moder:
            parts.append(f"moder {self.moder[0]} {self.moder[1]}")
        
        header = " ".join(parts)
        
        lines = [header]
        suffix = self.xs_suffix
        for zaid, density in self.composition.items():
            lines.append(f"    {zaid}{suffix} {density:.6E}")
        
        return "\n".join(lines)
    
    def __repr__(self):
        return f"S2_Material(name='{self.name}', n_isotopes={len(self.composition)}, T={self.temperature}K)"


# ═══════════════════════════════════════════════════════════════
#  S2_PinUniverse
# ═══════════════════════════════════════════════════════════════

class S2_PinUniverse:
    """Represents a single Serpent2 `pin` universe card.
    
    A pin is defined as a series of concentric cylindrical regions,
    each filled with a material, defined from innermost to outermost.
    The outermost region has no bounding radius.
    """
    
    def __init__(self, universe_name: str, material_names: list, radii: list):
        """
        Args:
            universe_name: Name/identifier for the pin universe.
            material_names: List of material names, length = len(radii) + 1.
                           Last material fills the outermost unbounded region.
            radii: List of outer radii (cm) for each bounded region.
        """
        if len(material_names) != len(radii) + 1:
            raise ValueError(
                f"Pin '{universe_name}': len(materials)={len(material_names)} "
                f"must equal len(radii)+1={len(radii)+1}"
            )
        self.universe_name = universe_name
        self.material_names = material_names
        self.radii = radii
    
    @classmethod
    def from_fuel_pin_model(cls, fuel_pin_model, pin_idx: int,
                            gap_material_name: str = "gap",
                            clad_material_name: str = "clad",
                            coolant_material_name: str = "cool"):
        """Create from a starterDD FuelPinModel.
        
        The FuelPinModel has:
            - radii: list of radii [fuel_zone_radii..., gap_outer_r, clad_outer_r]
            - fuel_material_mixtures: list of MaterialMixture for fuel zones
            - fuel_material_name: base name like 'UOX24'
        
        The pin universe will have regions:
            fuel_zone_1 | fuel_zone_2 | ... | fuel_zone_N | gap | clad | coolant
        
        Material names for fuel zones use the ``unique_material_mixture_name``
        from each ``MaterialMixture`` object, which preserves the Dragon naming
        scheme (e.g. ``UOX24_zone1_pin3``).
        
        Args:
            fuel_pin_model: starterDD FuelPinModel object.
            pin_idx: Unique index for this pin position.
            gap_material_name: Name of the gap material.
            clad_material_name: Name of the cladding material.
            coolant_material_name: Name of the coolant material.
        
        Returns:
            S2_PinUniverse instance.
        """
        pin = fuel_pin_model
        
        # Build material name list — use the unique mixture name that carries
        # zone and pin information, matching the Dragon naming scheme
        material_names = []
        for mix in pin.fuel_material_mixtures:
            mat_name = getattr(mix, 'unique_material_mixture_name',
                               mix.material_name)
            material_names.append(mat_name)
        
        material_names.append(gap_material_name)
        material_names.append(clad_material_name)
        material_names.append(coolant_material_name)
        
        # Build radii list (all bounded regions)
        # pin.radii should contain [fuel_zone_radii..., gap_r, clad_r]
        radii = list(pin.radii)
        
        # Universe name uses the base fuel material name + pin index
        universe_name = f"{pin.fuel_material_name}_{pin_idx}"
        
        return cls(universe_name, material_names, radii)
    
    @classmethod
    def water_rod(cls, universe_name: str, inner_radius: float,
                  outer_radius: float, wall_material: str = "clad",
                  water_material: str = "cool"):
        """Create a water rod pin universe.
        
        Structure: water | wall | water(surrounding)
        
        Args:
            universe_name: Name for the water rod universe.
            inner_radius: Inner radius of the water rod tube.
            outer_radius: Outer radius of the water rod tube.
            wall_material: Material name for the tube wall.
            water_material: Material name for water inside and outside.
        
        Returns:
            S2_PinUniverse instance.
        """
        return cls(
            universe_name=universe_name,
            material_names=[water_material, wall_material, water_material],
            radii=[inner_radius, outer_radius],
        )
    
    @classmethod
    def empty(cls, universe_name: str = "empty",
              fill_material: str = "cool"):
        """Create an empty lattice position (filled with coolant).
        
        Args:
            universe_name: Name for the empty universe.
            fill_material: Material to fill the position with.
        
        Returns:
            S2_PinUniverse instance.
        """
        return cls(
            universe_name=universe_name,
            material_names=[fill_material],
            radii=[],
        )
    
    def get_fuel_material_names(self) -> list:
        """Get the names of fuel materials only (excluding gap, clad, coolant).
        
        Convention: fuel materials are all materials before the last 3
        (gap, clad, coolant). For water rods and empty pins, returns [].
        
        Returns:
            List of fuel material name strings.
        """
        # For a standard fuel pin: materials = [fuel_zones..., gap, clad, cool]
        # Fuel zones are everything except the last 3
        if len(self.material_names) > 3:
            return self.material_names[:-3]
        return []
    
    def format_card(self) -> str:
        """Format the Serpent2 `pin` card.
        
        Returns:
            Formatted pin card string.
        """
        lines = [f"pin {self.universe_name}"]
        for i, radius in enumerate(self.radii):
            lines.append(f"    {self.material_names[i]}  {radius:.6f}")
        # Outermost material (no radius)
        lines.append(f"    {self.material_names[-1]}")
        return "\n".join(lines)
    
    def __repr__(self):
        return f"S2_PinUniverse('{self.universe_name}', n_regions={len(self.material_names)})"


# ═══════════════════════════════════════════════════════════════
#  S2_Lattice
# ═══════════════════════════════════════════════════════════════

class S2_Lattice:
    """Represents a Serpent2 `lat` card for a 2D square lattice.
    
    Type 1: infinite 2D square lattice.
    """
    
    def __init__(self, name: str, center_x: float, center_y: float,
                 nx: int, ny: int, pitch: float,
                 universe_map: list):
        """
        Args:
            name: Lattice universe name/identifier.
            center_x, center_y: Lattice center coordinates (cm).
            nx, ny: Number of elements in x and y.
            pitch: Lattice pitch (cm).
            universe_map: 2D list (ny rows x nx cols) of pin universe name strings.
                         Row ordering: y-decreasing (top row = most positive y).
        """
        if len(universe_map) != ny:
            raise ValueError(f"Lattice '{name}': universe_map has {len(universe_map)} rows, expected {ny}")
        for i, row in enumerate(universe_map):
            if len(row) != nx:
                raise ValueError(f"Lattice '{name}': row {i} has {len(row)} elements, expected {nx}")
        
        self.name = name
        self.center_x = center_x
        self.center_y = center_y
        self.nx = nx
        self.ny = ny
        self.pitch = pitch
        self.universe_map = universe_map
    
    @classmethod
    def from_assembly_model(cls, assembly_model, lattice_name: str = "10",
                            pin_universe_lookup: dict = None,
                            empty_universe_name: str = "empty"):
        """Build a lattice from a starterDD CartesianAssemblyModel.
        
        The assembly model has:
            - lattice: 2D list of pin model objects (FuelPinModel/DummyPinModel)
            - pin_geometry_dict with 'pin_pitch'
            - nx, ny dimensions
        
        DummyPinModel positions (water rod placeholders, vanished rods) are
        mapped to ``empty_universe_name``.  In Serpent2 the empty pin is
        typically filled with coolant; water-rod geometry is handled
        separately via surface/cell definitions that cut into the lattice.
        
        Args:
            assembly_model: starterDD CartesianAssemblyModel.
            lattice_name: Name for the lattice universe.
            pin_universe_lookup: Dict mapping pin model -> universe name string.
                                If None, uses default naming convention.
            empty_universe_name: Name for empty/dummy positions.
        
        Returns:
            S2_Lattice instance.
        """
        assembly = assembly_model
        ny = len(assembly.lattice)
        nx = len(assembly.lattice[0])
        pin_pitch = _get_pin_pitch(assembly)
        
        # Build universe map
        universe_map = []
        for j in range(ny):
            row = []
            for i in range(nx):
                pin = assembly.lattice[j][i]
                if pin_universe_lookup and pin in pin_universe_lookup:
                    row.append(pin_universe_lookup[pin])
                elif hasattr(pin, 'fuel_material_name') and hasattr(pin, 'pin_idx'):
                    # FuelPinModel — use naming convention matching from_fuel_pin_model
                    row.append(f"{pin.fuel_material_name}_{pin.pin_idx}")
                else:
                    # DummyPinModel or unknown → empty (filled with coolant)
                    row.append(empty_universe_name)
            universe_map.append(row)
        
        # Pad the lattice with one row/column of empty universes on each
        # side to account for the coolant strip between the outermost
        # fuel pins and the channel box inner wall.
        padded_nx = nx + 2
        padded_ny = ny + 2
        empty_row = [empty_universe_name] * padded_nx
        padded_map = [list(empty_row)]
        for row in universe_map:
            padded_map.append([empty_universe_name] + row + [empty_universe_name])
        padded_map.append(list(empty_row))
        universe_map = padded_map
        nx = padded_nx
        ny = padded_ny
        
        # Compute lattice center: lower-left assembly corner at (0,0),
        # everything in the +x/+y quadrant.
        assembly_pitch = getattr(assembly, 'assembly_pitch', nx * pin_pitch)
        center_x = assembly_pitch / 2.0
        center_y = assembly_pitch / 2.0
        
        return cls(
            name=lattice_name,
            center_x=center_x,
            center_y=center_y,
            nx=nx,
            ny=ny,
            pitch=pin_pitch,
            universe_map=universe_map,
        )
    
    def format_card(self) -> str:
        """Format the Serpent2 `lat` card.
        
        Returns:
            Formatted lat card string.
        """
        header = (f"lat {self.name}  1  {self.center_x:.6f} {self.center_y:.6f}  "
                  f"{self.nx} {self.ny}  {self.pitch:.4f}")
        
        lines = [header]
        # Find maximum universe name length for alignment
        max_len = max(
            len(name)
            for row in self.universe_map
            for name in row
        )
        
        for row in self.universe_map:
            formatted_row = "  ".join(f"{name:<{max_len}}" for name in row)
            lines.append(formatted_row)
        
        return "\n".join(lines)
    
    def __repr__(self):
        return f"S2_Lattice('{self.name}', {self.nx}x{self.ny}, pitch={self.pitch})"


# ═══════════════════════════════════════════════════════════════
#  S2_ChannelGeometry
# ═══════════════════════════════════════════════════════════════

class S2_ChannelGeometry:
    """Generates surf/cell cards for BWR channel box and outer boundary.
    
    Supports:
        - Rounded-corner square channel box (sqc surface)
        - Rectangular outer assembly boundary (cuboid)
        - Simple reflective boundary conditions
    """
    
    def __init__(self, center_x: float, center_y: float,
                 channel_inner_half_width: float,
                 channel_outer_half_width: float,
                 assembly_half_pitch: float,
                 corner_radius_inner: float = 0.0,
                 corner_radius_outer: float = 0.0,
                 lattice_universe_name: str = "10",
                 channel_box_material: str = "clad",
                 inner_water_material: str = "cool",
                 outer_water_material: str = "cool",
                 surf_id_start: int = 1000,
                 cell_id_start: int = 500,
                 water_rods: list = None):
        """
        Args:
            center_x, center_y: Center of the channel box (same as lattice center).
            channel_inner_half_width: Half-width of channel box inner surface.
            channel_outer_half_width: Half-width of channel box outer surface.
            assembly_half_pitch: Half of the full assembly pitch.
            corner_radius_inner: Corner radius for inner channel box (0 = sharp).
            corner_radius_outer: Corner radius for outer channel box.
            lattice_universe_name: Name of the lattice universe filling the channel.
            channel_box_material: Material for channel box walls.
            inner_water_material: Material for moderator inside the channel (in-channel water).
            outer_water_material: Material for water gap outside channel box.
            surf_id_start: Starting surface ID number.
            cell_id_start: Starting cell ID number.
            water_rods: Optional list of dicts describing water rod geometry:
                        ``[{"type": "circular", "center": (x, y),
                            "inner_radius": r_in, "outer_radius": r_out,
                            "moderator_material": "cool",
                            "wall_material": "clad", "rod_id": "WR_1"}, ...]``
                        or ``{"type": "square", ...}`` for square moderator boxes.
        """
        self.cx = center_x
        self.cy = center_y
        self.inner_hw = channel_inner_half_width
        self.outer_hw = channel_outer_half_width
        self.assembly_hp = assembly_half_pitch
        self.cr_inner = corner_radius_inner
        self.cr_outer = corner_radius_outer
        self.lattice_name = lattice_universe_name
        self.box_mat = channel_box_material
        self.inner_water_mat = inner_water_material
        self.outer_water_mat = outer_water_material
        self.sid = surf_id_start
        self.cid = cell_id_start
        self.water_rods = water_rods or []
    
    @classmethod
    def from_assembly_model(cls, assembly_model, lattice_universe_name: str = "10",
                            channel_box_material: str = "zr4",
                            inner_water_material: str = "moderator",
                            outer_water_material: str = "coolant"):
        """Create channel geometry from assembly model parameters.
        
        Handles:
            - Channel box inner/outer surfaces with optional corner radius
            - Assembly boundary surface
            - Water rod surface/cell overlays (circular or square)
        
        Expects assembly_model to have attributes like:
            - assembly_pitch
            - channel_box_inner_side (derived: assembly_pitch - 2*channel_box_thickness - 2*gap_wide)
            - channel_box_thickness
            - corner_inner_radius_of_curvature
            - water_rods: list of CircularWaterRodModel / SquareWaterRodModel
            - lattice: 2D grid for center computation
        """
        am = assembly_model
        pin_pitch = _get_pin_pitch(am)
        
        # Compute center: lower-left assembly corner at (0,0),
        # everything in the +x/+y quadrant.
        nx = len(am.lattice[0])
        ny = len(am.lattice)
        assembly_pitch = getattr(am, 'assembly_pitch', nx * pin_pitch)
        cx = assembly_pitch / 2.0
        cy = assembly_pitch / 2.0
        
        # ── Channel box dimensions ──────────────────────────────
        # The DragonModel (CartesianAssemblyModel) defines:
        #   - channel_box_thickness       (from YAML)
        #   - corner_inner_radius_of_curvature  (from YAML)
        #   - channel_box_inner_side      (derived: assembly_pitch - 2*cbt - 2*gap_wide)
        # We derive the Serpent2 half-widths from these.
        thickness = getattr(am, 'channel_box_thickness', None)
        cr_inner = getattr(am, 'corner_inner_radius_of_curvature', 0.0) or 0.0
        
        # Inner half-width from channel_box_inner_side
        channel_box_inner_side = getattr(am, 'channel_box_inner_side', None)
        if channel_box_inner_side is not None:
            inner_hw = channel_box_inner_side / 2.0
        else:
            inner_hw = nx * pin_pitch / 2.0 + 0.05  # small gap estimate
        
        # Outer half-width = inner + thickness
        if thickness is not None:
            outer_hw = inner_hw + thickness
        else:
            outer_hw = inner_hw + 0.2  # default 2mm thickness
        
        # Outer corner radius = inner corner radius + wall thickness
        cr_outer = cr_inner + (thickness if thickness is not None else 0.0)
        
        assembly_hp = getattr(am, 'assembly_pitch', None)
        if assembly_hp is None:
            assembly_hp = outer_hw + 0.5  # default gap
        else:
            assembly_hp = assembly_hp / 2.0
        
        # ── Extract water rod geometry ──────────────────────────────
        #
        # Compute the translation offset (gap + cbt + intra-assembly
        # coolant strip), mirroring glow_builder.build_full_assembly_geometry.
        gap = getattr(am, 'gap_wide', 0.0) or 0.0
        cbt = getattr(am, 'channel_box_thickness', 0.0) or 0.0
        channel_box_inner_side = getattr(am, 'channel_box_inner_side', None)
        if channel_box_inner_side is not None:
            intra_coolant = (channel_box_inner_side - nx * pin_pitch) / 2.0
        else:
            intra_coolant = 0.0
        translation = getattr(am, 'translation_offset', None)
        if translation is None:
            translation = gap + cbt + intra_coolant
        
        # Map Dragon-convention material names (set by
        # water_rod_model.set_materials) to the Serpent2 names chosen
        # by the caller.
        dragon_to_s2 = {
            "MODERATOR": inner_water_material,
            "CLAD": channel_box_material,
            "COOLANT": outer_water_material,
        }
        
        wr_list = []
        for wr in getattr(am, 'water_rods', []):
            wr_info = {"rod_id": getattr(wr, 'rod_ID', 'WR')}
            # Water rod centers in the YAML are already in cm in the
            # lattice-local coordinate system.  Apply the same
            # translation used in glow_builder to place them in the
            # assembly coordinate system (lower-left corner at 0,0).
            center = getattr(wr, 'center', None)
            if center is not None:
                wr_info["center"] = (center[0] + translation,
                                     center[1] + translation)
            else:
                wr_info["center"] = (cx, cy)
            
            # Circular water rod
            inner_r = getattr(wr, 'inner_radius', None)
            outer_r = getattr(wr, 'outer_radius', None)
            if inner_r is not None and outer_r is not None:
                wr_info["type"] = "circular"
                wr_info["inner_radius"] = inner_r
                wr_info["outer_radius"] = outer_r
            
            # Square water rod (moderator box)
            inner_side = getattr(wr, 'moderator_box_inner_side', None)
            outer_side = getattr(wr, 'moderator_box_outer_side', None)
            if inner_side is not None and outer_side is not None:
                wr_info["type"] = "square"
                wr_info["inner_half_side"] = inner_side / 2.0
                wr_info["outer_half_side"] = outer_side / 2.0
            
            # Materials — map Dragon names to Serpent2 names
            raw_mod = getattr(wr, 'moderator_material_name', None)
            raw_wall = getattr(wr, 'cladding_material_name', None)
            wr_info["moderator_material"] = dragon_to_s2.get(
                raw_mod, raw_mod if raw_mod else inner_water_material)
            wr_info["wall_material"] = dragon_to_s2.get(
                raw_wall, raw_wall if raw_wall else channel_box_material)
            
            if "type" in wr_info:
                wr_list.append(wr_info)
        
        return cls(
            center_x=cx,
            center_y=cy,
            channel_inner_half_width=inner_hw,
            channel_outer_half_width=outer_hw,
            assembly_half_pitch=assembly_hp,
            corner_radius_inner=cr_inner,
            corner_radius_outer=cr_outer,
            lattice_universe_name=lattice_universe_name,
            channel_box_material=channel_box_material,
            inner_water_material=inner_water_material,
            outer_water_material=outer_water_material,
            water_rods=wr_list,
        )
    
    def format_cards(self) -> str:
        """Format all surface and cell cards for the channel geometry.
        
        Includes:
            - Channel box inner/outer surfaces
            - Assembly boundary surface
            - Water rod surfaces (cylindrical or rectangular)
            - Cells: lattice fill, channel box, water gap, water rods, outside
        
        Returns:
            Complete string with surf and cell definitions.
        """
        lines = []
        lines.append("% --- Channel box and assembly boundary surfaces ---")
        
        s_inner = self.sid          # inner channel box
        s_outer = self.sid + 1      # outer channel box
        s_bndry = self.sid + 2      # assembly boundary
        next_sid = self.sid + 3     # next available surface ID
        
        # Inner channel box surface
        if self.cr_inner > 0:
            lines.append(f"surf {s_inner}  sqc  {self.cx:.6f}  {self.cy:.6f}  "
                        f"{self.inner_hw:.6f}  {self.cr_inner:.6f}")
        else:
            lines.append(f"surf {s_inner}  sqc  {self.cx:.6f}  {self.cy:.6f}  "
                        f"{self.inner_hw:.6f}")
        
        # Outer channel box surface
        if self.cr_outer > 0:
            lines.append(f"surf {s_outer}  sqc  {self.cx:.6f}  {self.cy:.6f}  "
                        f"{self.outer_hw:.6f}  {self.cr_outer:.6f}")
        else:
            lines.append(f"surf {s_outer}  sqc  {self.cx:.6f}  {self.cy:.6f}  "
                        f"{self.outer_hw:.6f}")
        
        # Assembly boundary
        lines.append(f"surf {s_bndry}  sqc  {self.cx:.6f}  {self.cy:.6f}  "
                    f"{self.assembly_hp:.6f}")
        
        # ── Water rod surfaces ──────────────────────────────────
        wr_surf_ids = []  # list of (inner_surf_id, outer_surf_id) per water rod
        for wr in self.water_rods:
            wr_cx, wr_cy = wr["center"]
            rod_id = wr.get("rod_id", "WR")
            
            if wr["type"] == "circular":
                s_wr_inner = next_sid
                s_wr_outer = next_sid + 1
                next_sid += 2
                lines.append(f"")
                lines.append(f"% --- Water rod {rod_id} ---")
                lines.append(f"surf {s_wr_inner}  cyl  {wr_cx:.6f}  {wr_cy:.6f}  "
                            f"{wr['inner_radius']:.6f}")
                lines.append(f"surf {s_wr_outer}  cyl  {wr_cx:.6f}  {wr_cy:.6f}  "
                            f"{wr['outer_radius']:.6f}")
                wr_surf_ids.append({
                    "type": "circular",
                    "inner": s_wr_inner,
                    "outer": s_wr_outer,
                    "moderator_material": wr["moderator_material"],
                    "wall_material": wr["wall_material"],
                    "rod_id": rod_id,
                })
            elif wr["type"] == "square":
                s_wr_inner = next_sid
                s_wr_outer = next_sid + 1
                next_sid += 2
                lines.append(f"")
                lines.append(f"% --- Moderator box {rod_id} ---")
                lines.append(f"surf {s_wr_inner}  sqc  {wr_cx:.6f}  {wr_cy:.6f}  "
                            f"{wr['inner_half_side']:.6f}")
                lines.append(f"surf {s_wr_outer}  sqc  {wr_cx:.6f}  {wr_cy:.6f}  "
                            f"{wr['outer_half_side']:.6f}")
                wr_surf_ids.append({
                    "type": "square",
                    "inner": s_wr_inner,
                    "outer": s_wr_outer,
                    "moderator_material": wr["moderator_material"],
                    "wall_material": wr["wall_material"],
                    "rod_id": rod_id,
                })
        
        lines.append("")
        lines.append("% --- Cells ---")
        
        next_cid = self.cid
        
        # Water rod cells (must be defined BEFORE the lattice fill cell
        # so they take priority over the lattice at the water rod positions)
        for wr_s in wr_surf_ids:
            c_mod = next_cid
            c_wall = next_cid + 1
            next_cid += 2
            lines.append(f"% Water rod / moderator box: {wr_s['rod_id']}")
            # Moderator inside the tube / box
            lines.append(f"cell {c_mod}  0  {wr_s['moderator_material']}  "
                        f"-{wr_s['inner']}")
            # Tube wall / box wall
            lines.append(f"cell {c_wall}  0  {wr_s['wall_material']}  "
                        f"{wr_s['inner']} -{wr_s['outer']}")
        
        # Collect all water-rod outer surface IDs — the lattice fill cell
        # must exclude these regions (pin is inside channel box but outside
        # all water rod outer surfaces)
        wr_outer_surfs = [ws["outer"] for ws in wr_surf_ids]
        
        c_lattice = next_cid
        c_box     = next_cid + 1
        c_gap     = next_cid + 2
        c_outside = next_cid + 3
        
        # Lattice fill cell: inside inner channel box, outside all water rods
        wr_exclusion = " ".join(str(s) for s in wr_outer_surfs)
        if wr_exclusion:
            lines.append(f"cell {c_lattice}  0  fill {self.lattice_name}  "
                        f"-{s_inner} {wr_exclusion}")
        else:
            lines.append(f"cell {c_lattice}  0  fill {self.lattice_name}  "
                        f"-{s_inner}")
        
        # Channel box wall
        lines.append(f"cell {c_box}  0  {self.box_mat}  {s_inner} -{s_outer}")
        # Water gap between channel box and assembly boundary
        lines.append(f"cell {c_gap}  0  {self.outer_water_mat}  "
                    f"{s_outer} -{s_bndry}")
        # Outside world
        lines.append(f"cell {c_outside}  0  outside  {s_bndry}")
        
        return "\n".join(lines)
    
    def __repr__(self):
        return (f"S2_ChannelGeometry(inner_hw={self.inner_hw}, "
                f"outer_hw={self.outer_hw}, assembly_hp={self.assembly_hp})")


# ═══════════════════════════════════════════════════════════════
#  S2_EnergyGrid
# ═══════════════════════════════════════════════════════════════

class S2_EnergyGrid:
    """Represents a Serpent2 `ene` card for detector energy binning.
    
    Type 1: arbitrary energy bin boundaries.
    """
    
    def __init__(self, name: str, boundaries: list):
        """
        Args:
            name: Energy grid name/identifier.
            boundaries: List of energy boundaries in MeV (ascending order).
                       For N energy groups, provide N+1 boundaries.
        
        Raises:
            ValueError: If fewer than 2 boundaries are given (need N+1
                       boundaries for N ≥ 1 energy groups) or if the
                       boundaries are not strictly increasing.
        """
        boundaries = sorted(boundaries)
        if len(boundaries) < 2:
            raise ValueError(
                f"Energy grid '{name}': at least 2 boundaries are required "
                f"(N+1 for N ≥ 1 groups), got {len(boundaries)}."
            )
        for i in range(len(boundaries) - 1):
            if boundaries[i] >= boundaries[i + 1]:
                raise ValueError(
                    f"Energy grid '{name}': boundaries must be strictly "
                    f"increasing, but boundary[{i}]={boundaries[i]} >= "
                    f"boundary[{i+1}]={boundaries[i+1]}."
                )
        self.name = name
        self.boundaries = boundaries
    
    @classmethod
    def from_preset(cls, preset_name: str,
                    name: str = None) -> 'S2_EnergyGrid':
        """Create an energy grid from one of the supported preset meshes.
        
        Supported presets (see ``SUPPORTED_ENERGY_MESHES``):
            - ``"2g"``  : 2-group mesh collapsed at 0.625 eV
            - ``"26g"`` : 26-group mesh condensed from SHEM-295
            - ``"295g"``: 295-group SHEM-295 mesh
        
        Args:
            preset_name: Key in ``SUPPORTED_ENERGY_MESHES``.
            name: Optional override for the grid name.  Defaults to
                  *preset_name* (e.g. ``"26g"``).
        
        Returns:
            S2_EnergyGrid instance.
        
        Raises:
            ValueError: If *preset_name* is not a recognised preset.
        """
        if preset_name not in SUPPORTED_ENERGY_MESHES:
            raise ValueError(
                f"Unknown energy mesh preset: '{preset_name}'. "
                f"Supported presets: {list(SUPPORTED_ENERGY_MESHES.keys())}. "
                f"For a custom mesh, provide a list of N+1 energy "
                f"boundaries (in MeV) instead."
            )
        grid_name = name if name is not None else preset_name
        return cls(grid_name, list(SUPPORTED_ENERGY_MESHES[preset_name]))
    
    @classmethod
    def full_range(cls, name: str = "full",
                   e_min: float = 1.0E-11,
                   e_max: float = 2.0E+1):
        """Create a single-bin energy grid covering the full energy range.
        
        This is the default for reaction rate detectors where energy-resolved
        results are not needed — scores are integrated over all energies.
        
        Args:
            name: Grid name.
            e_min: Minimum energy in MeV (default: 1e-11 MeV = 0.01 meV).
            e_max: Maximum energy in MeV (default: 20 MeV).
        
        Returns:
            S2_EnergyGrid instance with a single energy bin.
        """
        return cls(name, [e_min, e_max])
    
    @classmethod
    def two_group(cls, name: str = "2g",
                  thermal_cutoff: float = None,
                  e_min: float = None,
                  e_max: float = None):
        """Create a standard 2-group energy grid.
        
        When called without arguments the boundaries are taken from the
        ``"2g"`` preset in ``SUPPORTED_ENERGY_MESHES``.  Custom cutoff,
        minimum and maximum values can still be passed to override the
        preset.
        
        Args:
            name: Grid name.
            thermal_cutoff: Thermal/fast boundary in MeV (default from preset).
            e_min: Minimum energy in MeV (default from preset).
            e_max: Maximum energy in MeV (default from preset).
        
        Returns:
            S2_EnergyGrid instance.
        """
        if thermal_cutoff is None and e_min is None and e_max is None:
            return cls.from_preset("2g", name=name)
        # Fall back to explicit values when any override is given
        _e_min = e_min if e_min is not None else ene2g[0]
        _cutoff = thermal_cutoff if thermal_cutoff is not None else ene2g[1]
        _e_max = e_max if e_max is not None else ene2g[2]
        return cls(name, [_e_min, _cutoff, _e_max])
    
    @classmethod
    def from_dragon_energy_mesh(cls, name: str, energy_bounds_eV: list):
        """Create from Dragon5 energy mesh (in eV, descending order).
        
        Dragon5 convention: energy bounds are in eV, ordered from high to low.
        Serpent2 convention: energy bounds in MeV, ordered low to high.
        
        Args:
            name: Grid name.
            energy_bounds_eV: List of energy boundaries in eV (descending).
        
        Returns:
            S2_EnergyGrid instance.
        """
        # Convert eV -> MeV and reverse order
        boundaries_MeV = sorted([e * 1.0E-6 for e in energy_bounds_eV])
        return cls(name, boundaries_MeV)
    
    def format_card(self) -> str:
        """Format the Serpent2 `ene` card.
        
        Returns:
            Formatted ene card string.
        """
        bounds_str = "  ".join(f"{e:.4E}" for e in self.boundaries)
        return f"ene {self.name}  1  {bounds_str}"
    
    @property
    def n_groups(self) -> int:
        return len(self.boundaries) - 1
    
    def __repr__(self):
        return f"S2_EnergyGrid('{self.name}', {self.n_groups} groups)"


# ═══════════════════════════════════════════════════════════════
#  S2_Detector
# ═══════════════════════════════════════════════════════════════

class S2_Detector:
    """Represents a single Serpent2 `det` card.
    
    Supports scoring reaction rates for specified isotopes over specified
    material or universe domains, with energy binning.
    
    In Serpent2, a detector response is defined by:
        dr <MT> <material>  -- where material is a single-isotope "response material"
        dm <material>       -- domain material(s) where the reaction is scored
        du <universe>       -- domain universe where the reaction is scored
        de <energy_grid>    -- energy binning
        dt <flag>           -- detector type flag
    
    Using `du` (detector universe) is preferred when you want to score
    reaction rates integrated over an entire pin universe, as it automatically
    sums over all material regions within that universe.
    
    Using `dm` with `dt -4` provides cumulative scores over materials.
    """
    
    def __init__(self, name: str, energy_grid_name: str = None,
                 domain_materials: list = None,
                 domain_universe: str = None,
                 responses: list = None,
                 detector_type: int = None):
        """
        Args:
            name: Detector name/identifier.
            energy_grid_name: Name of the energy grid (ene card) to use.
            domain_materials: List of material names over which to integrate.
            domain_universe: Universe name for universe-based scoring (du card).
                            If set, scores are integrated over the entire universe.
            responses: List of (MT_number, response_material_name) tuples.
                      MT_number: positive for standard MT, negative for special.
                      response_material_name: single-isotope material for dr card.
                      Use MT=-6 for total fission, MT=-2 for total absorption, etc.
            detector_type: dt flag value (e.g., -4 for cumulative over materials).
        """
        self.name = name
        self.energy_grid_name = energy_grid_name
        self.domain_materials = domain_materials or []
        self.domain_universe = domain_universe
        self.responses = responses or []
        self.detector_type = detector_type
    
    def add_response(self, mt_number: int, response_material_name: str):
        """Add a reaction response to this detector.
        
        Args:
            mt_number: MT reaction number (or Serpent2 special code).
            response_material_name: Name of the single-isotope response material.
        """
        self.responses.append((mt_number, response_material_name))
    
    def add_domain_material(self, material_name: str):
        """Add a material domain for integration.
        
        Args:
            material_name: Name of the material to integrate over.
        """
        self.domain_materials.append(material_name)
    
    def format_card(self) -> str:
        """Format the Serpent2 `det` card.
        
        Returns:
            Formatted det card string.
        """
        lines = [f"det {self.name}"]
        
        if self.energy_grid_name:
            lines.append(f"    de {self.energy_grid_name}")
        
        if self.detector_type is not None:
            lines.append(f"    dt {self.detector_type}")
        
        # Universe-based domain (preferred for pin-wise scoring)
        if self.domain_universe:
            lines.append(f"    du {self.domain_universe}")
        
        # Material-based domains
        if self.domain_materials:
            dm_parts = "  ".join(f"dm {m}" for m in self.domain_materials)
            lines.append(f"    {dm_parts}")
        
        # Responses
        for mt, resp_mat in self.responses:
            lines.append(f"    dr {mt}  {resp_mat}")
        
        return "\n".join(lines)
    
    @classmethod
    def for_pin(cls, pin_universe: 'S2_PinUniverse',
                reaction_isotope_map: dict = None,
                energy_grid_name: str = None,
                detector_type: int = -4,
                suffix: str = "",
                # Deprecated parameters for backward compatibility
                isotopes: list = None,
                reactions: list = None):
        """Create a detector for a specific pin universe.
        
        Creates responses for each (reaction, isotope) combination as specified
        in the reaction_isotope_map. Uses material-based scoring (`dm` cards)
        with the fuel materials of the pin. Combined with `dt -4`, this sums
        scores over all fuel material zones within the pin, matching the Dragon
        _by_pin numbering convention.
        
        Only fuel materials are included in the domain (gap, clad, coolant
        are excluded), ensuring reaction rates are scored in the fuel pellet only.
        
        Args:
            pin_universe: S2_PinUniverse whose fuel materials define the
                         scoring domain.
            reaction_isotope_map: Dict mapping reaction names to lists of isotopes.
                         e.g., {'fission': ['U235', 'U238'], 'absorption': ['U235', 'Gd155']}
                         Reaction names are converted to MT numbers internally.
            energy_grid_name: Energy grid to use.
            detector_type: dt flag (default -4 to sum over all fuel zones).
            suffix: Optional suffix for detector name.
            isotopes: DEPRECATED. Use reaction_isotope_map instead.
            reactions: DEPRECATED. Use reaction_isotope_map instead.
        
        Returns:
            S2_Detector instance.
        """
        det_name = f"det_{pin_universe.universe_name}{suffix}"
        det = cls(
            name=det_name,
            energy_grid_name=energy_grid_name,
            domain_materials=pin_universe.get_fuel_material_names(),
            detector_type=detector_type,
        )
        
        # Handle new API: reaction_isotope_map
        if reaction_isotope_map is not None:
            for reaction_name, isotope_list in reaction_isotope_map.items():
                mt = _reaction_name_to_mt(reaction_name)
                for iso in isotope_list:
                    det.add_response(mt, iso)
        # Handle deprecated API: isotopes + reactions
        elif isotopes is not None and reactions is not None:
            warnings.warn(
                "The 'isotopes' and 'reactions' parameters are deprecated. "
                "Use 'reaction_isotope_map' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            for reaction in reactions:
                mt = _reaction_name_to_mt(reaction) if isinstance(reaction, str) else reaction
                for iso in isotopes:
                    det.add_response(mt, iso)
        
        return det
    
    @classmethod
    def for_assembly(cls, all_fuel_material_names: list,
                     reaction_isotope_map: dict,
                     energy_grid_name: str = None,
                     detector_type: int = -4,
                     name: str = None):
        """Create a single detector integrated over all fuel materials.
        
        Produces one detector whose ``dm`` cards list every fuel material
        in the assembly and whose ``dr`` lines cover every
        (reaction, isotope) pair in *reaction_isotope_map*.  Combined
        with ``dt -4`` this yields spatially-integrated, energy-resolved
        reaction rates for the whole assembly.
        
        Serpent2 outputs a result array of shape
        ``(n_responses × n_energy_bins)`` for this detector.
        
        Args:
            all_fuel_material_names: Flat list of fuel material names
                collected from every pin universe in the model.
            reaction_isotope_map: Dict mapping reaction names to lists of
                isotopes, e.g.
                ``{'fission': ['U238'], 'n,gamma': ['U238']}``.
            energy_grid_name: Name of the ``ene`` card to reference.
            detector_type: ``dt`` flag (default ``-4``, cumulative over
                ``dm`` materials).
            name: Detector name.  Defaults to
                ``"det_assembly_{energy_grid_name}"``.
        
        Returns:
            S2_Detector instance.
        """
        if name is None:
            name = f"det_assembly_{energy_grid_name}" if energy_grid_name else "det_assembly"
        
        det = cls(
            name=name,
            energy_grid_name=energy_grid_name,
            domain_materials=list(all_fuel_material_names),
            detector_type=detector_type,
        )
        
        for reaction_name, isotope_list in reaction_isotope_map.items():
            mt = _reaction_name_to_mt(reaction_name)
            for iso in isotope_list:
                det.add_response(mt, iso)
        
        return det
    
    def __repr__(self):
        n_domains = len(self.domain_materials) + (1 if self.domain_universe else 0)
        return (f"S2_Detector('{self.name}', {len(self.responses)} responses, "
                f"{n_domains} domains)")


# ═══════════════════════════════════════════════════════════════
#  S2_IsotopeResponseMaterial
# ═══════════════════════════════════════════════════════════════

class S2_IsotopeResponseMaterial:
    """Single-isotope material used as a response function in detectors.
    
    In Serpent2, the `dr` card in a detector references a material.
    For isotope-specific reaction rates, we create single-isotope materials
    at unit density to serve as the response function.
    
    Following the AT10/GE14 reference convention, the material name is the
    isotope name itself (e.g., ``U235``, ``Gd155``), not a suffixed variant.
    This allows detectors to use ``dr -6 U235`` directly.
    """
    
    def __init__(self, isotope_name: str, temperature: float = 900.0,
                 xs_suffix: str = None):
        """
        Args:
            isotope_name: Dragon-style isotope name (e.g., 'U235').
            temperature: Temperature for cross-section lookup.
            xs_suffix: Override cross-section suffix.
        """
        self.isotope_name = isotope_name
        self.temperature = temperature
        self._xs_suffix = xs_suffix
        
        # Derive ZAID
        self.zaid = isotope_name_to_zaid_str(isotope_name)
        # Use isotope name directly as the material name (AT10 convention)
        self.material_name = isotope_name
    
    @property
    def xs_suffix(self) -> str:
        if self._xs_suffix:
            return self._xs_suffix
        return get_xs_suffix(self.temperature)
    
    def format_card(self) -> str:
        """Format the single-isotope material card.
        
        Density is set to 1.0 (arbitrary, since it's only used as a response).
        
        Returns:
            Formatted mat card string.
        """
        return (f"mat {self.material_name}  1.0  tmp {self.temperature:.1f}\n"
                f"    {self.zaid}{self.xs_suffix}  1.0")
    
    def __repr__(self):
        return f"S2_IsotopeResponseMaterial('{self.isotope_name}' -> {self.material_name})"


# ═══════════════════════════════════════════════════════════════
#  S2_ThermalScattering
# ═══════════════════════════════════════════════════════════════

class S2_ThermalScattering:
    """Represents a Serpent2 `therm` card for thermal scattering laws."""
    
    def __init__(self, name: str, library_name: str):
        """
        Args:
            name: Thermal scattering law identifier (e.g., 'lwtr').
            library_name: Library file name (e.g., 'lwtr.16t').
        """
        self.name = name
        self.library_name = library_name
    
    def format_card(self) -> str:
        return f"therm {self.name}  {self.library_name}"
    
    def __repr__(self):
        return f"S2_ThermalScattering('{self.name}', '{self.library_name}')"


# ═══════════════════════════════════════════════════════════════
#  S2_Settings
# ═══════════════════════════════════════════════════════════════

class S2_Settings:
    """Global Serpent2 simulation settings.
    
    Collects library paths, neutron population, boundary conditions,
    burnup options, and plotting directives.
    """
    
    def __init__(self):
        """Initialize with default settings."""
        self.title = "starterDD generated Serpent2 model"
        
        # Nuclear data libraries
        self.acelib = None
        self.declib = None
        self.nfylib = None
        
        # Boundary conditions (1=vacuum, 2=reflective, 3=periodic)
        self.bc = 2  # reflective for lattice calculations
        
        # Neutron population
        self.neutrons_per_cycle = 50000
        self.active_cycles = 500
        self.inactive_cycles = 100
        
        # Power density (W/gHM) for burnup calculations
        self.power_density = None
        
        # Burnup steps (list of burnup points in MWd/kgHM)
        self.burnup_steps = None
        
        # Burnup mode (1=TTA, 2=CRAM)
        self.bumode = 2
        
        # Predictor-corrector (0=off, 1=CE/LI, 2=CE/CE)
        self.pcc = 1
        
        # Unresolved resonance probability tables
        self.ures = True
        
        # Group constant generation universe (0 = none, list = universe names)
        self.gcu = None
        
        # Few-group energy boundaries for group constant generation
        self.nfg = None
        
        # Geometry plotting
        self.plot_commands = []
        
        # Mesh plotting
        self.mesh_commands = []
    
    def set_endfb8r1_libraries(self, base_path: str):
        """Set ENDF/B-VIII.0r1 library paths.
        
        Args:
            base_path: Base directory for the nuclear data files.
        """
        self.acelib = f"{base_path}/endfb8r1_pynjoy2012_KERMA.xsdata"
        self.declib = f"{base_path}/endfb8r1.dec"
        self.nfylib = f"{base_path}/endfb8r1.nfy"
    
    def set_jeff311_libraries(self, base_path: str):
        """Set JEFF-3.1.1 library paths.
        
        Args:
            base_path: Base directory for the nuclear data files.
        """
        self.acelib = f"{base_path}/JEFF-311_pynjoy2016.xsdata"
        self.declib = f"{base_path}/sss_jeff311.dec"
        self.nfylib = f"{base_path}/sss_jeff311.nfy"
    
    def add_plot(self, plot_type: int = 3, x_pixels: int = 1000,
                 y_pixels: int = 1000, position: float = 0.5):
        """Add a geometry plot command.
        
        Args:
            plot_type: 1=xz, 2=yz, 3=xy plane.
            x_pixels, y_pixels: Image resolution.
            position: Position along the perpendicular axis.
        """
        self.plot_commands.append(f"plot {plot_type} {x_pixels} {y_pixels} {position}")
    
    def format_cards(self) -> str:
        """Format all settings as Serpent2 input cards.
        
        Returns:
            Complete settings block string.
        """
        lines = []
        
        # Sanitize title to pure ASCII (Serpent2 requires ASCII input files)
        safe_title = self.title.encode('ascii', errors='replace').decode('ascii')
        lines.append(f'set title "{safe_title}"')
        lines.append("")
        
        if self.acelib:
            lines.append(f'% --- Cross section library file path:')
            lines.append(f'set acelib "{self.acelib}"')
            lines.append("")
        
        if self.declib:
            lines.append(f'% --- Decay data library:')
            lines.append(f'set declib "{self.declib}"')
            lines.append("")
        
        if self.nfylib:
            lines.append(f'% --- Neutron-induced fission yield library:')
            lines.append(f'set nfylib "{self.nfylib}"')
            lines.append("")
        
        lines.append(f"% --- Boundary conditions:")
        lines.append(f"set bc {self.bc}")
        lines.append("")
        
        lines.append(f"% --- Neutron population:")
        lines.append(f"set pop {self.neutrons_per_cycle} {self.active_cycles} "
                     f"{self.inactive_cycles}")
        lines.append("")
        
        if self.ures:
            lines.append(f"% --- Unresolved resonance probability tables:")
            lines.append(f"set ures 1")
            lines.append("")
        
        if self.gcu is not None:
            if isinstance(self.gcu, (list, tuple)):
                gcu_str = " ".join(str(g) for g in self.gcu)
                lines.append(f"% --- Group constant generation universes:")
                lines.append(f"set gcu {gcu_str}")
            else:
                lines.append(f"% --- Group constant generation:")
                lines.append(f"set gcu {self.gcu}")
            lines.append("")
        
        if self.nfg is not None:
            if isinstance(self.nfg, (list, tuple)):
                nfg_str = " ".join(f"{e:.6E}" for e in self.nfg)
                lines.append(f"% --- Few-group energy boundaries:")
                lines.append(f"set nfg {len(self.nfg) + 1} {nfg_str}")
            else:
                lines.append(f"set nfg {self.nfg}")
            lines.append("")
        
        if self.power_density is not None and self.burnup_steps:
            lines.append(f"% --- Burnup settings:")
            lines.append(f"set powdens {self.power_density:.6E}")
            lines.append(f"set bumode {self.bumode}")
            lines.append(f"set pcc {self.pcc}")
            bu_str = " ".join(f"{b:.2f}" for b in self.burnup_steps)
            lines.append(f"dep butot {bu_str}")
            lines.append("")
        
        for cmd in self.plot_commands:
            lines.append(cmd)
        
        for cmd in self.mesh_commands:
            lines.append(cmd)
        
        return "\n".join(lines)
    
    def __repr__(self):
        return f"S2_Settings(title='{self.title}')"


# ═══════════════════════════════════════════════════════════════
#  Serpent2Model — Top-level orchestrator
# ═══════════════════════════════════════════════════════════════

class Serpent2Model:
    """Orchestrates the generation of a complete Serpent2 input deck
    from a starterDD AssemblyModel.
    
    Workflow:
        1. Create Serpent2Model with an assembly model
        2. Call build() to generate all cards
        3. Optionally call add_detector_config() for reaction rate tallies
        4. Call write() to output the .serp file
    """
    
    def __init__(self, assembly_model=None, settings: S2_Settings = None):
        """
        Args:
            assembly_model: A starterDD CartesianAssemblyModel (or compatible object).
            settings: S2_Settings instance. Created with defaults if None.
        """
        self.assembly = assembly_model
        self.settings = settings or S2_Settings()
        
        # Generated cards
        self.materials = []                # list of S2_Material
        self.pin_universes = []            # list of S2_PinUniverse
        self.lattice = None                # S2_Lattice
        self.channel_geometry = None       # S2_ChannelGeometry
        self.energy_grids = []             # list of S2_EnergyGrid
        self.detectors = []                # list of S2_Detector
        self.isotope_response_materials = []  # list of S2_IsotopeResponseMaterial
        self.thermal_scattering_laws = []  # list of S2_ThermalScattering
        
        # Internal tracking
        self._pin_universe_map = {}  # pin_idx -> S2_PinUniverse
        self._built = False
    
    def build(self, gap_material_name: str = "gap",
              clad_material_name: str = "clad",
              coolant_material_name: str = "cool",
              outer_water_material_name: str = "cool_outer",
              channel_box_material_name: str = "zr4",
              lattice_name: str = "10",
              empty_universe_name: str = "empty",
              xs_suffix_fuel: str = None,
              xs_suffix_struct: str = None):
        """Build all Serpent2 cards from the assembly model.
        
        Args:
            gap_material_name: Serpent2 name for gap material.
            clad_material_name: Serpent2 name for cladding material.
            coolant_material_name: Serpent2 name for in-channel coolant.
            outer_water_material_name: Serpent2 name for outer water.
            channel_box_material_name: Serpent2 name for channel box.
            lattice_name: Name for the lattice universe.
            empty_universe_name: Name for empty lattice positions.
            xs_suffix_fuel: Override XS suffix for fuel materials.
            xs_suffix_struct: Override XS suffix for structural materials.
        """
        if self.assembly is None:
            raise RuntimeError("No assembly model provided. Cannot build.")
        
        self._build_fuel_materials(xs_suffix=xs_suffix_fuel)
        self._build_pin_universes(
            gap_material_name=gap_material_name,
            clad_material_name=clad_material_name,
            coolant_material_name=coolant_material_name,
        )
        self._build_empty_universe(empty_universe_name, coolant_material_name)
        self._build_water_rod_universes(coolant_material_name, clad_material_name)
        self._build_lattice(lattice_name, empty_universe_name)
        self._build_channel_geometry(
            lattice_name=lattice_name,
            channel_box_material=channel_box_material_name,
            outer_water_material=outer_water_material_name,
        )
        
        self._built = True
    
    def _build_fuel_materials(self, xs_suffix: str = None):
        """Create S2_Material for each fuel MaterialMixture in the assembly."""
        if not hasattr(self.assembly, 'fuel_material_mixtures'):
            print("Warning: assembly model has no 'fuel_material_mixtures' attribute.")
            return
        
        for mix in self.assembly.fuel_material_mixtures:
            mat = S2_Material.from_material_mixture(mix, xs_suffix=xs_suffix)
            self.materials.append(mat)
    
    def add_structural_material(self, name: str, isotope_densities: dict,
                                temperature: float, xs_suffix: str = None,
                                rgb: tuple = None, moder: tuple = None,
                                total_density: float = None):
        """Manually add a structural (non-fuel) material.
        
        Args:
            name: Material name.
            isotope_densities: Dict of {isotope_name_or_zaid: number_density}.
            temperature: Temperature in Kelvin.
            xs_suffix: Cross-section suffix override.
            rgb: RGB color for plotting.
            moder: Thermal scattering law tuple (therm_name, nuclide_zaid).
            total_density: Total atom density (if not using 'sum').
        """
        mat = S2_Material.from_raw(
            name=name,
            isotope_densities=isotope_densities,
            temperature=temperature,
            is_burnable=False,
            xs_suffix=xs_suffix,
            rgb=rgb,
            moder=moder,
        )
        if total_density is not None:
            mat.density_type = 'total'
            mat.total_density = total_density
        
        self.materials.append(mat)
    
    def build_structural_materials_from_assembly(
            self, name_map: dict = None,
            temperature_map: dict = None,
            default_temperature: float = 600.0,
            xs_suffix: str = None):
        """Auto-generate structural (non-fuel) materials from the assembly's
        ``composition_lookup``.
        
        This inspects the assembly's ``composition_lookup`` dictionary and
        creates ``S2_Material`` entries for every composition whose name is
        NOT already present in ``self.materials`` (i.e., not a fuel material).
        
        Args:
            name_map: Optional dict mapping Dragon composition name to the
                      desired Serpent2 material name.
                      E.g. ``{"COOLANT": "cool", "CLAD": "clad"}``.
                      Compositions not in the map keep their original name
                      (lowercased).
            temperature_map: Optional dict mapping composition name to
                             temperature (K).
            default_temperature: Temperature to use when no mapping is given.
            xs_suffix: Cross-section suffix override.
        """
        comp_lookup = getattr(self.assembly, 'composition_lookup', None)
        if comp_lookup is None:
            return
        
        name_map = name_map or {}
        temperature_map = temperature_map or {}
        existing_names = {m.name for m in self.materials}
        
        # Collect base fuel material names so that their parent
        # compositions (e.g. '24UOX', '32UOX', '45Gd') are not emitted
        # as standalone materials — only the per-pin numbered variants
        # (e.g. '24UOX_zone1_pin1') are used in pin universes.
        fuel_base_names = set()
        for row in getattr(self.assembly, 'lattice', []):
            for pin in row:
                fname = getattr(pin, 'fuel_material_name', None)
                if fname:
                    fuel_base_names.add(fname)
                    fuel_base_names.add(fname.upper())
                    fuel_base_names.add(fname.lower())
        
        for comp_name, composition in comp_lookup.items():
            # Skip compositions that already have a material (fuel mixes)
            s2_name = name_map.get(comp_name, comp_name.lower())
            if s2_name in existing_names:
                continue
            
            # Skip base fuel compositions (not assigned to any cell)
            if comp_name in fuel_base_names or s2_name in fuel_base_names:
                continue
            
            temp = temperature_map.get(comp_name, default_temperature)
            iso_dict = composition.get_isotope_name_composition()
            if not iso_dict:
                continue
            
            mat = S2_Material.from_raw(
                name=s2_name,
                isotope_densities=iso_dict,
                temperature=temp,
                is_burnable=False,
                xs_suffix=xs_suffix,
            )
            self.materials.append(mat)
            existing_names.add(s2_name)
    
    def add_thermal_scattering(self, name: str, library_name: str):
        """Add a thermal scattering law.
        
        Args:
            name: Thermal law name (e.g., 'lwtr').
            library_name: Library file (e.g., 'lwtr.16t').
        """
        self.thermal_scattering_laws.append(
            S2_ThermalScattering(name, library_name)
        )
    
    def _build_pin_universes(self, gap_material_name: str,
                             clad_material_name: str,
                             coolant_material_name: str):
        """Create S2_PinUniverse for each unique fuel pin in the assembly."""
        seen_pin_idx = set()
        
        for row in self.assembly.lattice:
            for pin in row:
                # Check if it's a fuel pin
                if not hasattr(pin, 'fuel_material_mixtures'):
                    continue
                
                pin_idx = getattr(pin, 'pin_idx', None)
                if pin_idx is None:
                    continue
                
                if pin_idx in seen_pin_idx:
                    continue
                seen_pin_idx.add(pin_idx)
                
                univ = S2_PinUniverse.from_fuel_pin_model(
                    pin, pin_idx,
                    gap_material_name=gap_material_name,
                    clad_material_name=clad_material_name,
                    coolant_material_name=coolant_material_name,
                )
                self.pin_universes.append(univ)
                self._pin_universe_map[pin_idx] = univ
    
    def _build_empty_universe(self, name: str, coolant_material: str):
        """Create an empty pin universe for non-fuel lattice positions."""
        empty = S2_PinUniverse.empty(universe_name=name,
                                      fill_material=coolant_material)
        self.pin_universes.append(empty)
    
    def _build_water_rod_universes(self, coolant_material: str,
                                    wall_material: str):
        """Handle water rod geometry from the assembly.
        
        Water rods are stored in ``assembly.water_rods`` (a list of
        CircularWaterRodModel or SquareWaterRodModel objects) — they are
        **not** represented in the lattice cells.  Lattice positions
        occupied by water rods are DummyPinModel placeholders that become
        empty (coolant-only) universes in the lattice.
        
        The actual water rod geometry (tubes, moderator boxes) is handled
        via surface and cell definitions appended to the channel geometry
        by ``_build_channel_geometry``.  This follows the AT10/GE14
        reference convention where water rods overlay the lattice rather
        than being embedded in it.
        
        Vanished rod cases (no ``water_rods`` attribute, or empty list)
        are left as empty coolant universes — no extra action needed.
        """
        # Nothing to do here — water rod surface/cell geometry is
        # handled in _build_channel_geometry.  This method is retained
        # for API compatibility and the docstring above.
        pass
    
    def _build_lattice(self, lattice_name: str, empty_universe_name: str):
        """Build the lattice card from the assembly lattice grid.
        
        DummyPinModel positions (water rod placeholders or vanished rods)
        are all filled with the empty/coolant universe.  Actual water rod
        geometry is handled by surface/cell overlays in the channel
        geometry section (AT10/GE14 convention).
        """
        pin_pitch = _get_pin_pitch(self.assembly)
        
        universe_map = []
        for row in self.assembly.lattice:
            univ_row = []
            for pin in row:
                if hasattr(pin, 'fuel_material_mixtures') and hasattr(pin, 'pin_idx'):
                    pin_idx = pin.pin_idx
                    if pin_idx in self._pin_universe_map:
                        univ_row.append(self._pin_universe_map[pin_idx].universe_name)
                    else:
                        univ_row.append(empty_universe_name)
                else:
                    # DummyPinModel or other non-fuel position → coolant
                    univ_row.append(empty_universe_name)
            universe_map.append(univ_row)
        
        ny = len(universe_map)
        nx = len(universe_map[0]) if ny > 0 else 0
        
        # Pad the lattice with one row/column of empty universes on each
        # side to account for the coolant strip between the outermost
        # fuel pins and the channel box inner wall.
        padded_nx = nx + 2
        padded_ny = ny + 2
        empty_row = [empty_universe_name] * padded_nx
        padded_map = [list(empty_row)]
        for row in universe_map:
            padded_map.append([empty_universe_name] + row + [empty_universe_name])
        padded_map.append(list(empty_row))
        universe_map = padded_map
        nx = padded_nx
        ny = padded_ny
        
        # Compute lattice center: lower-left assembly corner at (0,0),
        # everything in the +x/+y quadrant.
        assembly_pitch = getattr(self.assembly, 'assembly_pitch', nx * pin_pitch)
        center_x = assembly_pitch / 2.0
        center_y = assembly_pitch / 2.0
        
        self.lattice = S2_Lattice(
            name=lattice_name,
            center_x=center_x,
            center_y=center_y,
            nx=nx,
            ny=ny,
            pitch=pin_pitch,
            universe_map=universe_map,
        )
    
    def _build_channel_geometry(self, lattice_name: str,
                                channel_box_material: str,
                                outer_water_material: str):
        """Build channel box geometry from assembly parameters."""
        self.channel_geometry = S2_ChannelGeometry.from_assembly_model(
            self.assembly,
            lattice_universe_name=lattice_name,
            channel_box_material=channel_box_material,
            outer_water_material=outer_water_material,
        )
    
    def add_detector_config(self, reaction_isotope_map: dict = None,
                            energy_bounds: list = None,
                            energy_grid_name: str = "full",
                            fuel_temperature: float = 900.0,
                            xs_suffix: str = None,
                            detector_type: int = -4,
                            # Deprecated parameters for backward compatibility
                            isotopes: list = None,
                            reactions: list = None):
        """Configure detectors for reaction rate scoring.
        
        Creates:
            - An energy grid (ene card)
            - Single-isotope response materials for each tracked isotope
            - One detector per unique fuel pin
        
        This is analogous to the Dragon5 EDI/COMPO export functionality.
        
        Args:
            reaction_isotope_map: Dict mapping reaction names to lists of isotopes.
                         e.g., {'fission': ['U235', 'U238'], 'absorption': ['Gd155', 'Gd157']}
                         Reaction names are converted to MT numbers internally.
                         Only isotopes with data for each reaction should be listed.
            energy_bounds: Energy boundaries in MeV for binning (N+1 values
                          for N energy groups).  When *None* the grid is
                          resolved from *energy_grid_name*: if it matches a
                          supported preset (``"2g"``, ``"26g"``, ``"295g"``)
                          the predefined boundaries are used automatically;
                          otherwise a single full-range bin is created.
            energy_grid_name: Name for the energy grid.  When it matches a
                          key in ``SUPPORTED_ENERGY_MESHES`` and no explicit
                          *energy_bounds* are given the preset is used.
                          Default: ``"full"`` (single bin, no preset match).
            fuel_temperature: Temperature for response materials.
            xs_suffix: Override XS suffix for response materials.
            detector_type: dt flag for detectors (default: -4 to sum over
                          all fuel zones in each pin).
            isotopes: DEPRECATED. Use reaction_isotope_map instead.
            reactions: DEPRECATED. Use reaction_isotope_map instead.
        
        Example:
            model.add_detector_config(
                reaction_isotope_map={
                    'fission': ['U234', 'U235', 'U236', 'U238'],
                    'absorption': ['U234', 'U235', 'U236', 'U238', 'Gd155', 'Gd157'],
                },
                energy_grid_name="26g",   # uses SHEM-295 condensed 26-group mesh
            )
        """
        # Handle deprecated API
        if reaction_isotope_map is None and isotopes is not None and reactions is not None:
            warnings.warn(
                "The 'isotopes' and 'reactions' parameters are deprecated. "
                "Use 'reaction_isotope_map' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            # Convert old API to new format (all isotopes for all reactions)
            reaction_isotope_map = {}
            for reaction in reactions:
                reaction_name = reaction if isinstance(reaction, str) else str(reaction)
                reaction_isotope_map[reaction_name] = list(isotopes)
        
        if reaction_isotope_map is None:
            raise ValueError(
                "Must provide 'reaction_isotope_map' argument. "
                "Example: {'fission': ['U235', 'U238'], 'absorption': ['Gd155']}"
            )
        
        # Build energy grid
        # Priority: explicit bounds > preset by name > single full-range bin
        if energy_bounds is not None:
            eg = S2_EnergyGrid(name=energy_grid_name, boundaries=energy_bounds)
        elif energy_grid_name in SUPPORTED_ENERGY_MESHES:
            eg = S2_EnergyGrid.from_preset(energy_grid_name)
        else:
            eg = S2_EnergyGrid.full_range(name=energy_grid_name)
        self.energy_grids.append(eg)
        
        # Collect all unique isotopes from the map
        all_isotopes = set()
        for isotope_list in reaction_isotope_map.values():
            all_isotopes.update(isotope_list)
        
        # Build isotope response materials (one per unique isotope)
        existing_response_names = {rm.isotope_name for rm in self.isotope_response_materials}
        for iso in all_isotopes:
            if iso not in existing_response_names:
                resp_mat = S2_IsotopeResponseMaterial(
                    isotope_name=iso,
                    temperature=fuel_temperature,
                    xs_suffix=xs_suffix,
                )
                self.isotope_response_materials.append(resp_mat)
        
        # Build detectors — one per unique fuel pin
        for univ in self.pin_universes:
            fuel_mats = univ.get_fuel_material_names()
            if not fuel_mats:
                continue  # Skip non-fuel universes
            
            det = S2_Detector.for_pin(
                pin_universe=univ,
                reaction_isotope_map=reaction_isotope_map,
                energy_grid_name=energy_grid_name,
                detector_type=detector_type,
            )
            self.detectors.append(det)
    
    def add_assembly_integrated_detector_config(
            self,
            reaction_isotope_map: dict,
            energy_grid_name: str = "295g",
            energy_bounds: list = None,
            fuel_temperature: float = 900.0,
            xs_suffix: str = None,
            detector_type: int = -4,
            name: str = None):
        """Configure a single detector integrated over all fuel materials.
        
        Creates one ``S2_Detector`` whose ``dm`` cards list **every** fuel
        material across all pin universes in the model.  All (reaction,
        isotope) pairs from *reaction_isotope_map* are added as ``dr``
        responses on that single detector.  Combined with ``dt -4`` this
        yields spatially-integrated, energy-resolved reaction rates for
        the whole assembly.
        
        Serpent2 outputs a result array of shape
        ``(n_responses × n_energy_bins)`` for the detector, giving
        energy-resolved rates for each response without spatial
        pin-by-pin breakdown.
        
        This method complements :meth:`add_detector_config` (which
        creates one detector **per pin**) and
        :meth:`add_flux_detector` (global scalar flux).
        
        Args:
            reaction_isotope_map: Dict mapping reaction names to lists of
                isotopes, e.g.
                ``{'fission': ['U238'], 'n,gamma': ['U238']}``.
            energy_grid_name: Name for the energy grid.  When it matches a
                key in ``SUPPORTED_ENERGY_MESHES`` and no explicit
                *energy_bounds* are given the preset is used.
            energy_bounds: Explicit energy boundaries (N+1 values in MeV).
                Takes priority over presets.
            fuel_temperature: Temperature for isotope response materials.
            xs_suffix: Override XS suffix for response materials.
            detector_type: ``dt`` flag (default ``-4``).
            name: Detector name.  Defaults to
                ``"det_assembly_{energy_grid_name}"``.
        
        Example::
        
            model.add_assembly_integrated_detector_config(
                reaction_isotope_map={
                    'fission':  ['U238'],
                    'n,gamma':  ['U238'],
                    'absorption': ['U238'],
                },
                energy_grid_name="295g",
                fuel_temperature=900.0,
            )
        """
        if reaction_isotope_map is None or len(reaction_isotope_map) == 0:
            raise ValueError(
                "Must provide a non-empty 'reaction_isotope_map' argument."
            )
        
        # ── Energy grid (reuse if already present) ──
        existing_grids = {eg.name for eg in self.energy_grids}
        if energy_grid_name not in existing_grids:
            if energy_bounds is not None:
                eg = S2_EnergyGrid(name=energy_grid_name,
                                   boundaries=energy_bounds)
            elif energy_grid_name in SUPPORTED_ENERGY_MESHES:
                eg = S2_EnergyGrid.from_preset(energy_grid_name)
            else:
                eg = S2_EnergyGrid.full_range(name=energy_grid_name)
            self.energy_grids.append(eg)
        
        # ── Isotope response materials (dedup) ──
        all_isotopes = set()
        for isotope_list in reaction_isotope_map.values():
            all_isotopes.update(isotope_list)
        
        existing_response_names = {
            rm.isotope_name for rm in self.isotope_response_materials
        }
        for iso in all_isotopes:
            if iso not in existing_response_names:
                resp_mat = S2_IsotopeResponseMaterial(
                    isotope_name=iso,
                    temperature=fuel_temperature,
                    xs_suffix=xs_suffix,
                )
                self.isotope_response_materials.append(resp_mat)
        
        # ── Collect all fuel materials across every pin ──
        all_fuel_mats = []
        for univ in self.pin_universes:
            all_fuel_mats.extend(univ.get_fuel_material_names())
        
        if not all_fuel_mats:
            warnings.warn(
                "No fuel materials found in pin universes. "
                "Assembly-integrated detector will have no domain materials.",
                stacklevel=2,
            )
        
        # ── Build a single detector with all responses ──
        det = S2_Detector.for_assembly(
            all_fuel_material_names=all_fuel_mats,
            reaction_isotope_map=reaction_isotope_map,
            energy_grid_name=energy_grid_name,
            detector_type=detector_type,
            name=name,
        )
        self.detectors.append(det)
    
    def add_flux_detector(self, energy_grid_name: str = "2g",
                          name: str = "flux",
                          energy_bounds: list = None):
        """Add a global flux detector (no domain restriction).
        
        Scores the scalar flux integrated over the whole geometry,
        binned by the specified energy grid.
        
        Args:
            energy_grid_name: Energy grid to use.  If it matches a key in
                          ``SUPPORTED_ENERGY_MESHES`` (e.g. ``"2g"``,
                          ``"26g"``, ``"295g"``) and no explicit
                          *energy_bounds* are given, the preset boundaries
                          are used automatically.
            name: Detector name.
            energy_bounds: Optional energy boundaries (N+1 values in MeV).
                          If provided and no matching grid already exists,
                          a new grid is created.
        """
        # Ensure an energy grid exists
        existing_grids = {eg.name for eg in self.energy_grids}
        if energy_grid_name not in existing_grids:
            if energy_bounds is not None:
                eg = S2_EnergyGrid(name=energy_grid_name,
                                   boundaries=energy_bounds)
            elif energy_grid_name in SUPPORTED_ENERGY_MESHES:
                eg = S2_EnergyGrid.from_preset(energy_grid_name)
            else:
                eg = S2_EnergyGrid.two_group(name=energy_grid_name)
            self.energy_grids.append(eg)
        
        det = S2_Detector(name=name, energy_grid_name=energy_grid_name)
        self.detectors.append(det)
    
    def write(self, filepath: str):
        """Write the complete Serpent2 input file.
        
        Args:
            filepath: Output file path (e.g., 'assembly.serp').
        """
        sections = []
        
        # 1. Settings (title, libraries, population, BC)
        sections.append(self.settings.format_cards())
        sections.append("")
        
        # 2. Materials
        sections.append("% " + "=" * 60)
        sections.append("%  MATERIAL DEFINITIONS")
        sections.append("% " + "=" * 60)
        sections.append("")
        
        for mat in self.materials:
            sections.append(mat.format_card())
            sections.append("")
        
        # Thermal scattering laws
        for therm in self.thermal_scattering_laws:
            sections.append(therm.format_card())
        if self.thermal_scattering_laws:
            sections.append("")
        
        # 3. Pin universes
        sections.append("% " + "=" * 60)
        sections.append("%  PIN UNIVERSE DEFINITIONS")
        sections.append("% " + "=" * 60)
        sections.append("")
        
        for pin in self.pin_universes:
            sections.append(pin.format_card())
            sections.append("")
        
        # 4. Lattice
        if self.lattice:
            sections.append("% " + "=" * 60)
            sections.append("%  LATTICE DEFINITION")
            sections.append("% " + "=" * 60)
            sections.append("")
            sections.append(self.lattice.format_card())
            sections.append("")
        
        # 5. Channel geometry
        if self.channel_geometry:
            sections.append("% " + "=" * 60)
            sections.append("%  GEOMETRY (CHANNEL BOX & BOUNDARY)")
            sections.append("% " + "=" * 60)
            sections.append("")
            sections.append(self.channel_geometry.format_cards())
            sections.append("")
        
        # 6. Energy grids
        if self.energy_grids:
            sections.append("% " + "=" * 60)
            sections.append("%  ENERGY GRIDS")
            sections.append("% " + "=" * 60)
            sections.append("")
            for eg in self.energy_grids:
                sections.append(eg.format_card())
            sections.append("")
        
        # 7. Isotope response materials
        if self.isotope_response_materials:
            sections.append("% " + "=" * 60)
            sections.append("%  ISOTOPE RESPONSE MATERIALS (for detectors)")
            sections.append("% " + "=" * 60)
            sections.append("")
            for irm in self.isotope_response_materials:
                sections.append(irm.format_card())
                sections.append("")
        
        # 8. Detectors
        if self.detectors:
            sections.append("% " + "=" * 60)
            sections.append("%  DETECTOR DEFINITIONS")
            sections.append("% " + "=" * 60)
            sections.append("")
            for det in self.detectors:
                sections.append(det.format_card())
                sections.append("")
        
        # Write to file — enforce ASCII encoding for Serpent2 compatibility
        content = "\n".join(sections)
        with open(filepath, 'w', encoding='ascii', errors='replace') as f:
            f.write(content)
        
        print(f"Serpent2 input written to: {filepath}")
    
    def summary(self) -> str:
        """Print a summary of the model.
        
        Returns:
            Summary string.
        """
        lines = [
            "=" * 60,
            "  Serpent2 Model Summary",
            "=" * 60,
            f"  Title:        {self.settings.title}",
            f"  Materials:    {len(self.materials)}",
            f"  Pin universes:{len(self.pin_universes)}",
            f"  Lattice:      {self.lattice}",
            f"  Energy grids: {len(self.energy_grids)}",
            f"  Detectors:    {len(self.detectors)}",
            f"  Response mats:{len(self.isotope_response_materials)}",
            "=" * 60,
        ]
        return "\n".join(lines)
    
    def __repr__(self):
        return (f"Serpent2Model(materials={len(self.materials)}, "
                f"pins={len(self.pin_universes)}, "
                f"detectors={len(self.detectors)})")
