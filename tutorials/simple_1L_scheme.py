# Example of a simple DRAGON 1 level flux calulation scheme with starterDD
# 

# Date : 11/03/2026
# R.Guasch

from pathlib import Path

from starterDD.InterfaceToDD.case_generator import DragonCase

# =====================================================================
# Configuration paths (relative to the project root)
# =====================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
GE14_INPUTS = PROJECT_ROOT / "data" / "BWRProgressionProblems" / "GE14" / "inputs"
TDT_FILES = PROJECT_ROOT / "tests" / "reference_tdt_files" / "GE14"

GE14_DOM_test_case = DragonCase(
        case_name="GE14_DOM",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": str(GE14_INPUTS / "material_compositions.yaml"),
            "GEOM": str(GE14_INPUTS / "GEOM_GE14_DOM.yaml"),
            "CALC_SCHEME": str(GE14_INPUTS / "CALC_SCHEME_1L.yaml"),
        },
        output_path="outputs/GE14_DOM/1L_scheme",
        tdt_path=str(TDT_FILES),
    )

result = GE14_DOM_test_case.generate_cle2000_procedures()

print(result)