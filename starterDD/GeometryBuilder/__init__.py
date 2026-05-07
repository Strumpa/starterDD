"""GeometryBuilder module - provides geometry building utilities.

This module works in two modes:
1. With glow: Full functionality including glow_builder for SALOME/TDT geometry generation
2. Standalone: Only helpers module for abstract geometry calculations (no glow dependency)
"""

# Import standalone helpers (always available)
from .helpers import computeSantamarinaRadii

__all__ = ["computeSantamarinaRadii"]

# Flag to check if glow integration is available
GLOW_AVAILABLE = False

# Try to import glow-dependent functionality
try:
    from .glow_builder import (
        build_full_assembly_geometry,
    )
    GLOW_AVAILABLE = True
    __all__.extend([
        "build_full_assembly_geometry",
        "GLOW_AVAILABLE",
    ])
except ImportError:
    # glow not installed - standalone mode
    pass

__all__.append("GLOW_AVAILABLE")