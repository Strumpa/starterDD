"""GeometryAnalysis module - TDT file parsing and geometry analysis."""

from .tdt_parser import read_material_mixture_indices_from_tdt_file
from .geometry_analysis import build_assembly_geometry, analyse_mesh, GeometricAnalyser

__all__ = ["read_material_mixture_indices_from_tdt_file", "build_assembly_geometry", "analyse_mesh", "GeometricAnalyser"]
