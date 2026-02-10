"""starterDD - A package for Donjon&Dragon model building and glow interoperability.

This package provides tools for:
- Material property and composition handling (MaterialProperties)
- Geometry building with optional glow integration (GeometryBuilder)  
- Dragon model definition (DragonModel)
- Geometry analysis from TDT files (GeometryAnalysis)
- Interface to Dragon&Donjon (InterfaceToDD)

The package works in two modes:
1. Standalone: Works without glow, uses abstract classes for Dragon model definition
2. With glow: Full integration with glow for SALOME geometry generation and TDT export
"""

# Import submodules
from starterDD import MaterialProperties
from starterDD import GeometryBuilder
from starterDD import DDModel
from starterDD import GeometryAnalysis
from starterDD import InterfaceToDD

# Re-export the glow availability flag
from starterDD.GeometryBuilder import GLOW_AVAILABLE

__version__ = "0.1.0"
__all__ = [
    "MaterialProperties",
    "GeometryBuilder", 
    "DDModel",
    "GeometryAnalysis",
    "InterfaceToDD",
    "GLOW_AVAILABLE",
]