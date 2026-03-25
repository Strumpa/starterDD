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
DATA_DIR = (TESTS_DIR / ".." / "data").resolve()     # .../starterDD/data
REF_TDT_DIR = (TESTS_DIR / "reference_tdt_files").resolve()  # .../starterDD/tests/reference_tdt_files
OUTPUTS_DIR = (TESTS_DIR / "outputs").resolve()      # .../starterDD/tests/outputs

# ── GE14 inputs ───────────────────────────────────────────────────────────
GE14_INPUTS_DIR = DATA_DIR / "BWRProgressionProblems" / "GE14_inputs"
# CORE TESTS
GE14_CORE_YAML = str(GE14_INPUTS_DIR / "CORE")
# ASSEMBLY TESTS
GE14_COMPOSITIONS_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" / "material_compositions.yaml")
GE14_SIMPLE_GEOMETRY_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" / "simplified_geometry.yaml")
GE14_DOM_GEOMETRY_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" / "GEOM_GE14_DOM.yaml")
GE14_DOM_C_GEOMETRY_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" / "GEOM_GE14_DOM-C.yaml")
GE14_CALC_SCHEME_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" / "CALC_SCHEME_GE14.yaml")
GE14_CALC_SCHEME_1L_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" /  "CALC_SCHEME_1L.yaml")
GE14_CALC_BRANCHES_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" /  "CALC_BRANCHES.yaml")
GE14_CALC_OUTPUTS_YAML = str(GE14_INPUTS_DIR / "ASSEMBLY" /  "CALC_OUTPUTS.yaml")
GE14_INPUTS_DIR_STR = str(GE14_INPUTS_DIR)

GE14_TDT_DIR = str(REF_TDT_DIR / "GE14")

# ── ATRIUM-10 inputs ─────────────────────────────────────────────────────
AT10_INPUTS_DIR = DATA_DIR / "ATRIUM10_inputs"
AT10_COMPOSITIONS_YAML = str(AT10_INPUTS_DIR / "ASSEMBLY" / "material_compositions.yaml")
AT10_GEOMETRY_YAML = str(AT10_INPUTS_DIR / "ASSEMBLY" / "GEOM_ATRIUM10.yaml")
AT10_GEOMETRY_CTRL_YAML = str(AT10_INPUTS_DIR / "ASSEMBLY" / "GEOM_ATRIUM10_CTRL.yaml")
AT10_CALC_SCHEME_YAML = str(AT10_INPUTS_DIR / "ASSEMBLY" / "CALC_SCHEME_AT10.yaml")
AT10_TDT_DIR = str(REF_TDT_DIR / "AT10")

# ── ATRIUM-10 4x4 sub-assembly inputs ────────────────────────────────────
AT10_4x4_INPUTS_DIR = AT10_INPUTS_DIR / "SUB-ASSEMBLY"
AT10_4x4_GEOMETRY_YAML = str(AT10_4x4_INPUTS_DIR / "GEOM_4x4.yaml")
AT10_4x4_COMPOSITIONS_YAML = str(AT10_4x4_INPUTS_DIR / "material_compositions.yaml")
AT10_4x4_CALC_SCHEME_2L_SPLIT_YAML = str(AT10_4x4_INPUTS_DIR / "CALC_SCHEME_2L_split.yaml")
AT10_4x4_CALC_SCHEME_2L_BY_PIN_YAML = str(AT10_4x4_INPUTS_DIR / "CALC_SCHEME_2L_by_pin.yaml")
AT10_4x4_CALC_SCHEME_1L_YAML = str(AT10_4x4_INPUTS_DIR / "CALC_SCHEME_1L.yaml")
