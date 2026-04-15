"""GeometryAnalysis module - TDT file parsing and geometry analysis."""

from .tdt_parser import read_material_mixture_indices_from_tdt_file
from .geometry_analysis import build_global_geometry, analyse_mesh, analyse_3d_volume, GeometricAnalyser

__all__ = ["read_material_mixture_indices_from_tdt_file", "build_global_geometry", "analyse_mesh", "analyse_3d_volume", "GeometricAnalyser"]
