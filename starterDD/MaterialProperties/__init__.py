"""MaterialProperties module - handles material compositions and mixtures.

This module works in two modes:
1. With iapws: Full water property functions (THMSAT, THMTX, THMPH, THMSAP)
2. Without iapws: Only utility functions like compute_water_iso_densities_at_densities
"""

# Import always-available functionality
from .material_mixture import (
    Composition,
    MaterialMixture,
    XSData,
    parse_all_compositions_from_yaml,
    get_element_symbol,
    get_isotope_atomic_mass,
    fractions_to_iso_densities,
    HM_isotopes,
    AVOGADRO,
    CM2_TO_BARN,
)
from .water_properties import compute_water_iso_densities_at_densities

__all__ = [
    "Composition",
    "MaterialMixture",
    "XSData",
    "parse_all_compositions_from_yaml",
    "get_element_symbol",
    "get_isotope_atomic_mass",
    "fractions_to_iso_densities",
    "compute_water_iso_densities_at_densities",
    "HM_isotopes",
    "AVOGADRO",
    "CM2_TO_BARN",
]

# Flag to check if iapws integration is available
IAPWS_AVAILABLE = False

# Try to import iapws-dependent functionality
try:
    from .water_properties import (
        THMSAT,
        THMTX,
        THMPH,
        THMSAP,
    )
    IAPWS_AVAILABLE = True
    __all__.extend([
        "THMSAT",
        "THMTX",
        "THMPH",
        "THMSAP",
        "IAPWS_AVAILABLE",
    ])
except ImportError:
    # iapws not installed - standalone mode
    pass

__all__.append("IAPWS_AVAILABLE")   

