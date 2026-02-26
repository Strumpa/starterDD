"""
Shared pytest configuration for starterDD tests.

All data paths are anchored to ``Path(__file__).parent`` so that tests work
regardless of the working directory from which ``pytest`` is invoked
(e.g. ``starterDD/`` or ``starterDD/tests/``).
"""

from pathlib import Path

import pytest

# ── directory anchors ─────────────────────────────────────────────────────
TESTS_DIR = Path(__file__).resolve().parent          # .../starterDD/tests
print(f"TESTS_DIR: {TESTS_DIR}")  # --- DEBUG ---
DATA_DIR = (TESTS_DIR / ".." / "data").resolve()     # .../starterDD/data
print(f"DATA_DIR: {DATA_DIR}")  # --- DEBUG ---
REF_TDT_DIR = (TESTS_DIR / "reference_tdt_files").resolve()  # .../starterDD/tests/reference_tdt_files
print(f"REF_TDT_DIR: {REF_TDT_DIR}")  # --- DEBUG ---
OUTPUTS_DIR = (TESTS_DIR / "outputs").resolve()      # .../starterDD/tests/outputs
print(f"OUTPUTS_DIR: {OUTPUTS_DIR}")  # --- DEBUG ---

# ── GE14 inputs ───────────────────────────────────────────────────────────
GE14_INPUTS_DIR = DATA_DIR / "BWRProgressionProblems" / "GE14" / "inputs"
GE14_COMPOSITIONS_YAML = str(GE14_INPUTS_DIR / "material_compositions.yaml")
GE14_SIMPLE_GEOMETRY_YAML = str(GE14_INPUTS_DIR / "simplified_geometry.yaml")
GE14_DOM_GEOMETRY_YAML = str(GE14_INPUTS_DIR / "GEOM_GE14_DOM.yaml")
GE14_CALC_SCHEME_YAML = str(GE14_INPUTS_DIR / "CALC_SCHEME_GE14.yaml")
GE14_INPUTS_DIR_STR = str(GE14_INPUTS_DIR)

GE14_TDT_DIR = str(REF_TDT_DIR / "GE14")

# ── ATRIUM-10 inputs ─────────────────────────────────────────────────────
AT10_INPUTS_DIR = DATA_DIR / "ATRIUM10" / "inputs"
AT10_COMPOSITIONS_YAML = str(AT10_INPUTS_DIR / "material_compositions.yaml")
AT10_GEOMETRY_YAML = str(AT10_INPUTS_DIR / "GEOM_ATRIUM10.yaml")
AT10_CALC_SCHEME_YAML = str(AT10_INPUTS_DIR / "CALC_SCHEME_AT10.yaml")
