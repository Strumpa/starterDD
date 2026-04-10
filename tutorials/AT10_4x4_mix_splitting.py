# Example of a simple DRAGON 2-levels flux calculation scheme with starterDD
# Generate a glow geometry with a single AT10 4x4 sub-assembly with 16 pins,
# Run dragon with a RSE+IC self shielding step + first level IC on 295g + second level 26g MOC flux calculation.
# At self shielding step and 1st level flux calculation, each pin is divided in Santamarina radial zones, 
# Material Mixes are numbered by material ie enrichment / fuel type. 
# Cross sections are condensed from 295g to 26g after the first level flux calculation, 
# Mixes are duplicated to the "by pin" numbering scheme and assigned to the respective pins,
# The second level flux calculation is performed on the 26g cross sections with the "by pin" mixes.

# Date : 11/03/2026
# R.Guasch

from pathlib import Path
import os

try:
    from glow.support.types import GeometryType, PropertyType
    from starterDD.starterDD.InterfaceToDD.case_generator import DragonCase
    GLOW_AVAILABLE = True
    
except ImportError:
    
    GLOW_AVAILABLE = False
    from starterDD.InterfaceToDD.case_generator import DragonCase

# =====================================================================
# Configuration paths — anchored to the project root so the script
# works regardless of the working directory it is launched from.
#
# Inside glow/SALOME, __file__ is not defined (the script is exec'd
# in an embedded interpreter), so we fall back to the well-known
# Docker mount point.
# =====================================================================
try:
    PROJECT_ROOT = Path(__file__).resolve().parent.parent
except NameError:
    # Running inside glow/SALOME — CWD is /home/user/data/
    PROJECT_ROOT = Path("/home/user/data/starterDD")

AT10_INPUTS = PROJECT_ROOT / "data" / "ATRIUM10_inputs" / "SUB-ASSEMBLY"
DRAGON_EXEC = os.environ.get('dragon_exec', 'path/to/dragon_executable')
DRAGLIBS_PATH = Path(os.environ.get('DRAGLIB_DIR', "/path/to/draglibs"))

# glow_data sits next to the starterDD project root
GLOW_DATA = PROJECT_ROOT.parent / "glow_data"
AT10_OUTPUT = GLOW_DATA / "starterDD_outputs" / "AT10_4x4" / "2L_scheme"
run_dragon=False  # Set to False for a dry run (no Dragon execution)

AT10_4x4_test_case = DragonCase(
        case_name="AT10_4x4_mix_splitting",
        call_glow=True, # Set to True to call glow for geometry processing.
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295_v5p1": "endfb8r1_295",
        },
        config_yamls={
            "MATS": str(AT10_INPUTS / "material_compositions.yaml"),
            "GEOM": str(AT10_INPUTS / "GEOM_4x4.yaml"),
            "CALC_SCHEME": str(AT10_INPUTS / "CALC_SCHEME_2L_split.yaml"),
        },
        output_path=str(AT10_OUTPUT),
        tdt_path=str(AT10_OUTPUT),
    )

AT10_4x4_test_case.set_fuel_material_temperatures({
    "24UOX": 900.0,
    "32UOX": 900.0,
    "42UOX": 900.0,
    "45UOX": 900.0,
    "48UOX": 900.0,
    "50UOX": 900.0,
    "45Gd" : 900.0,
})
AT10_4x4_test_case.set_non_fuel_temperatures(structural_temperature=600.0, coolant_temperature=600.0, moderator_temperature=600.0, gap_temperature=600.0)


# Step 1: Generate CLE2000 procedures (x2m + c2m files)
result = AT10_4x4_test_case.generate_cle2000_procedures()

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
if not run_dragon:
    # --- Option A: dry run (no Dragon execution) -----------------------
    # Useful for verifying the setup before running.
    #
    dry_result = AT10_4x4_test_case.run(
         draglib_paths={
             "draglibendfb8r1SHEM295_v5p1": (DRAGLIBS_PATH / "draglibendfb8r1SHEM295_v5p1"),
         },
         results_root=str(PROJECT_ROOT / "tutorials" / "results"),
         dry_run=True,
    )
    print(f"Dry run directory: {dry_result.run_directory}")

# --- Option B: full execution --------------------------------------
# Requires $dragon_exec and draglib files to be available.
#
if run_dragon:
    print("Running Dragon... This may take a few moments.")
    print(f"Using Dragon executable: {DRAGON_EXEC}")
    run_result = AT10_4x4_test_case.run(
        dragon_executable=DRAGON_EXEC,  # or None to use $dragon_exec
        draglib_paths={
            "draglibendfb8r1SHEM295_v5p1": (DRAGLIBS_PATH / "draglibendfb8r1SHEM295_v5p1"),
        },
        results_root=str(PROJECT_ROOT / "tutorials" / "results"),
        num_threads=1,
    )
    print(f"Draglibs path used: {DRAGLIBS_PATH / 'draglibendfb8r1SHEM295_v5p1'}")
    print("Dragon run completed.")
    print(f"Success: {run_result.success}")
    print(f"keff:    {run_result.keff}")
    print(f"Time:    {run_result.wall_time_seconds:.1f}s")
    print(f"Results: {run_result.run_directory}")
#
# Results directory structure:
#   results/
#   └── AT10_4x4/
#       └── RSE_IC_2L_MOC_5sp/
#           ├── 2026-03-12T14-32-07/
#           │   ├── run_manifest.yaml
#           │   ├── inputs/           (frozen config yamls)
#           │   ├── procedures/       (generated CLE2000 files)
#           │   ├── AT10_4x4.result   (Dragon listing)
#           │   └── _CPO_AT10_4x4     (COMPO output)
#           └── latest -> 2026-03-12T14-32-07/