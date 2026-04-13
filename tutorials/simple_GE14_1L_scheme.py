# Example of a simple DRAGON 1 level flux calulation scheme with starterDD
# 

# Date : 11/03/2026
# R.Guasch

from pathlib import Path
import os

from starterDD.InterfaceToDD.case_generator import DragonCase

# =====================================================================
# Configuration paths (relative to the project root)
# =====================================================================
try:
    PROJECT_ROOT = Path(__file__).resolve().parent.parent
except NameError:
    # Running inside glow/SALOME — CWD is /home/user/data/
    PROJECT_ROOT = Path("/home/user/data/starterDD")
GE14_INPUTS = PROJECT_ROOT / "data" / "BWRProgressionProblems" / "GE14_inputs" / "ASSEMBLY"
TDT_FILES = PROJECT_ROOT / "tests" / "reference_tdt_files" / "GE14"
DRAGON_EXEC = os.environ.get('dragon_exec', None)
DRAGLIBS_PATH = Path(os.environ.get('DRAGLIB_DIR', "/path/to/draglibs"))

run_dragon=True  # Set to False for a dry run (no Dragon execution)
run_glow = False  # Set to True to call glow for geometry processing (if needed)

GE14_DOM_test_case = DragonCase(
        case_name="GE14_DOM",
        call_glow=run_glow,
        draglib_name_to_alias={
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

# Step 1: Generate CLE2000 procedures (x2m + c2m files)
result = GE14_DOM_test_case.generate_cle2000_procedures()
print("Generated procedures:", result)

# =====================================================================
# Step 2: Execute the case with the Dragon runner
#
# The runner replaces the manual rdragon + .access + .save workflow.
# It handles: executable resolution, input staging (TDT files,
# draglibs), execution in a temp directory, output collection,
# and traceability (manifest + config archival).
#
# Requirements:
#   - Dragon executable: set $dragon_exec or pass explicitly
#   - Draglib files: set $DRAGLIBS directory or pass draglib_paths
#
# Use dry_run=True to stage everything without executing Dragon.
# =====================================================================

# --- Option A: dry run (no Dragon execution) -----------------------
# Useful for verifying the setup before running.
if not run_dragon:
    dry_result = GE14_DOM_test_case.run(
         draglib_paths={
             "draglibendfb8r1SHEM295_v5p1": (DRAGLIBS_PATH / "draglibendfb8r1SHEM295_v5p1"),
         },
         results_root=str(PROJECT_ROOT / "tutorials" / "results"),
         dry_run=True,
    )
    print(f"Dry run directory: {dry_result.run_directory}")
# --- Option B: full execution --------------------------------------
# Requires $dragon_exec and draglib files to be available.
if run_dragon:
     run_result = GE14_DOM_test_case.run(
         dragon_executable=DRAGON_EXEC,  # or None to use $dragon_exec
         draglib_paths={
             "draglibendfb8r1SHEM295_v5p1": (DRAGLIBS_PATH / "draglibendfb8r1SHEM295_v5p1"),
         },
         results_root=str(PROJECT_ROOT / "tutorials" / "results"),
         num_threads=1,
     )
     print(f"Success: {run_result.success}")
     print(f"keff:    {run_result.keff}")
     print(f"Time:    {run_result.wall_time_seconds:.1f}s")
     print(f"Results: {run_result.run_directory}")

# Results directory structure:
#   results/
#   └── GE14_DOM/
#       └── RSE_IC_1L_MOC_5sp/
#           ├── 2026-03-12T14-32-07/
#           │   ├── run_manifest.yaml
#           │   ├── inputs/           (frozen config yamls)
#           │   ├── procedures/       (generated CLE2000 files)
#           │   ├── GE14_DOM.result   (Dragon listing)
#           │   └── _CPO_GE14_DOM     (COMPO output)
#           └── latest -> 2026-03-12T14-32-07/