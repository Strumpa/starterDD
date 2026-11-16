# Tutorial: a DRAGON5 / Serpent2 code-to-code study with starterDD
#
# A study is the product  assemblies x scheme variants x operating conditions.
# Each (assembly, variant, condition) is one DRAGON case (one multi-statepoint
# COMPO); each statepoint gets one Serpent2 reference, shared by all variants.
# The study keeps a record per case and a study_manifest.yaml, so that runs
# can be split over cluster jobs, skipped when done, and post-processed from
# the manifest alone.
#
#   Part A  runs anywhere (committed data only, ~5 s): a small study built in
#           Python (GE14 DOM, 1L scheme, committed TDT files):
#           instancing, cases and statepoints, dry run, Serpent2 inputs,
#           status / manifest, adding a variant, post-processing calls.
#   Part B  the GE14 2L study of data/BWRProgressionProblems/GE14_inputs/STUDIES/
#           GE14_2L_PT3x3_IC3x3_6x6_HZP_HFP.yaml (DOM, DOM-C, VAN x HZP, HFP),
#           read from its YAML. Needs the TDT files under ${GE14_WORK}
#           (GE14_<A>/tdt_files), set with the other machine paths in the
#           site file (~/.config/starterDD/site.yaml, see docs/site.example.yaml). By default: dry run of the 6 cases and the
#           12 Serpent2 inputs. With --run it runs DRAGON (needs $dragon_exec
#           and the draglib, see the YAML 'dragon:' block).
#
# Results go to tutorials/results/code_to_code_study/ (git-ignored), or --root.
#
# Run:
#   python tutorials/code_to_code_study.py                          # A, then B (dry)
#   python tutorials/code_to_code_study.py --part B --run \
#          --assembly GE14_DOM --condition HZP --threads 10        # one DRAGON case
#
# The same steps from the shell (one case per cluster job):
#   python -m starterDD.Workflow STUDY.yaml status
#   python -m starterDD.Workflow STUDY.yaml run --dry-run
#   python -m starterDD.Workflow STUDY.yaml run --assembly GE14_VAN --condition HFP --threads 10
#   python -m starterDD.Workflow STUDY.yaml serpent
#   python -m starterDD.Workflow STUDY.yaml import-dragon   # 'dragon: import:' block
#   python -m starterDD.Workflow STUDY.yaml post
#   python -m starterDD.Workflow STUDY.yaml figures
#
# Date : 15/11/2026

import argparse
from pathlib import Path

import yaml

from starterDD.Workflow import (CodeToCodeStudy, StudyAssembly, StudyCondition,
                                StudyVariant)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GE14 = PROJECT_ROOT / "data" / "BWRProgressionProblems" / "GE14_inputs"
ASSEMBLY = GE14 / "ASSEMBLY"
STUDY_YAML = GE14 / "STUDIES" / "GE14_2L_PT3x3_IC3x3_6x6_HZP_HFP.yaml"
REF_TDT = PROJECT_ROOT / "tests" / "reference_tdt_files" / "GE14"


def title(text):
    print("\n" + f" {text} ".center(72, "="))


def show_tree(root: Path, depth: int = 6):
    """Print the study folder, without procedures/inputs contents."""
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root)
        if len(rel.parts) > depth or {"procedures", "inputs"} & set(rel.parts[:-1]):
            continue
        print("   " + "  " * (len(rel.parts) - 1) + p.name + ("/" if p.is_dir() else ""))


# =====================================================================
# PART A - a small study built in Python (committed data)
# =====================================================================
def part_a(out: Path):
    # -----------------------------------------------------------------
    # A.1  Instancing: assemblies, conditions, variants, spec
    # -----------------------------------------------------------------
    title("A.1 Instancing a study")
    coolant = {v: ASSEMBLY / f"material_compos_{v}.yaml" for v in ("00", "40", "80")}
    study = CodeToCodeStudy(
        name="GE14_DOM_1L_demo",
        root=out,
        # One lattice: geometry, materials, TDT folder and its base scheme.
        assemblies={"GE14_DOM": StudyAssembly(
            name="GE14_DOM",
            geometry=ASSEMBLY / "GEOM_GE14_DOM.yaml",
            materials=ASSEMBLY / "material_compositions.yaml",
            tdt_path=REF_TDT,
            scheme=ASSEMBLY / "CALC_SCHEME_1L.yaml")},
        # Conditions = branch values. The coolant statepoints carry the exact
        # benchmark isotopics of each file (COOLANT entry).
        conditions={
            "HZP": StudyCondition("HZP", coolant={"00": coolant["00"]},
                                  fuel_temperature=[600.0], coolant_temperature=[600.0]),
            "HFP": StudyCondition("HFP", coolant=coolant,
                                  fuel_temperature=[900.0], coolant_temperature=[600.0]),
        },
        # Variants = step changes applied to every assembly's base scheme.
        variants={"1L_base": StudyVariant("1L_base", description="scheme as in the YAML")},
        # The spec sets the DRAGON outputs of every case and the Serpent detectors.
        spec=dict(preset="pin_rates", isotopes=["U235", "U238", "Gd155", "Gd157"]),
        library="endfb8r1",
        draglib={"draglibendfb8r1SHEM295": "endfb8r1_295"},
        serpent_settings=dict(neutrons_per_cycle=100000, active_cycles=500, inactive_cycles=100),
    )
    for case in study.cases():
        print(f"{case.key:28s} variant_tag={case.scheme.variant_tag}")
        for sp in case.statepoints:
            print(f"    {sp.label:8s} {sp.params}")
        print(f"    outputs: {[o.name for o in case.scheme.outputs]}")
        print(f"    folder : {case.directory.relative_to(out)}")

    # -----------------------------------------------------------------
    # A.2  Dry run: procedures + staged inputs, no DRAGON execution
    # -----------------------------------------------------------------
    title("A.2 DRAGON dry run (full tree, no execution)")
    print(study.run_dragon(dry_run=True))
    # study.run_dragon(num_threads=10) runs DRAGON; completed cases are then
    # skipped by later calls unless force=True (or their inputs changed).

    # -----------------------------------------------------------------
    # A.3  Serpent2 references: one input per (assembly, statepoint)
    # -----------------------------------------------------------------
    title("A.3 Serpent2 inputs (written, not run)")
    refs = study.write_serpent_references()      # only the new ones (existing are kept)
    for r in refs:
        print(f"   {r.input_path.relative_to(out)}")
    ref40 = study.registry.find("GE14_DOM", study.library, "HFP_40")
    print("reference.yaml of HFP_40:")
    print(yaml.safe_dump({k: ref40.meta[k] for k in ("statepoint", "coolant_composition",
                                                      "settings", "spec")}, sort_keys=False))
    print("coverage problems:", study.check_references() or "none")

    show_tree(study.directory)

    # -----------------------------------------------------------------
    # A.4  Status and manifest
    # -----------------------------------------------------------------
    title("A.4 Status and study_manifest.yaml")
    study.print_status()
    man = yaml.safe_load(study.manifest_path.read_text())
    print(yaml.safe_dump(man["cases"][1], sort_keys=False, width=100)[:1200], "...")

    # -----------------------------------------------------------------
    # A.5  Extending the study with a variant: the references are reused
    # -----------------------------------------------------------------
    title("A.5 Adding a variant")
    study2 = CodeToCodeStudy(
        name=study.name, root=out, assemblies=study.assemblies, conditions=study.conditions,
        variants={"1L_base": study.variants["1L_base"],
                  "1L_12az": StudyVariant("1L_12az", steps={"FLUX": {"num_angles_2d": 12}},
                                          description="12 azimuthal angles in the flux step")},
        spec=study.spec, library=study.library, draglib=study.draglib,
        serpent_settings=study.serpent_settings)
    study2.run_dragon(dry_run=True, variants="1L_12az")      # only the new cases
    study2.write_serpent_references()                       # nothing new: same statepoints
    study2.print_status()

    # -----------------------------------------------------------------
    # A.6  Post-processing (once DRAGON and Serpent have run)
    # -----------------------------------------------------------------
    title("A.6 Post-processing")
    rows = study2.post_process()
    print({r["status"] for r in rows})
    print("Once a case is completed and its references have _det0.m/_res.m, the same call\n"
          "writes post/study.h5, summary.csv and summary.json (keff, normalisation, pin\n"
          "rates, spectra), and study2.figures() draws post/figures/ from study.h5 alone.\n"
          "Default normalisation: fission_source (assembly nu-fission, same factor for\n"
          "every isotope set); references without 'void' responses fall back to\n"
          "power_subset, reported in the 'normalisation' row.\n"
          "An existing COMPO is used without re-running DRAGON with\n"
          "study2.import_dragon_run('GE14_DOM', '1L_base', 'HFP', '<run folder or _CPO_ file>').")


# =====================================================================
# PART B - the GE14 2L PT 3x3 / IC 3x3_6x6 study, HZP and HFP
# =====================================================================
def part_b(out: Path, args):
    title("B.1 Reading the study file")
    # Machine paths (GE14_WORK, DRAGLIB_DIR, SERPENT_XS, DRAGON_EXEC) come from
    # the site file of this machine (~/.config/starterDD/site.yaml, template
    # docs/site.example.yaml), the environment, or --set; never from the
    # committed study file. `python -m starterDD.Workflow STUDY.yaml vars`
    # shows each value and where it comes from.
    variables = {"STUDY_ROOT": str(out)}
    try:
        study = CodeToCodeStudy.from_yaml(STUDY_YAML, variables=variables)
    except KeyError as exc:
        print(f"{exc.args[0]}\npart B skipped.")
        return
    missing = [a.name for a in study.assemblies.values() if not a.tdt_path.exists()]
    if missing:
        print(f"TDT folders not found for {missing} (check GE14_WORK with `vars`); part B skipped.")
        return
    sel = dict(assemblies=args.assembly, conditions=args.condition)
    print(f"study folder: {study.directory}")
    study.print_status()

    title("B.2 DRAGON: dry run of the selected cases")
    study.run_dragon(dry_run=True, **sel)

    if args.import_serpent:
        title("B.3a Registering the existing NSE Serpent2 runs")
        study.import_serpent_references()

    title("B.3 Serpent2 inputs")
    from starterDD.PostTreatment.physics_data import UnknownReactionsError
    try:
        study.write_serpent_references(**sel)
    except UnknownReactionsError as exc:     # e.g. U236: add 'serpent: reaction_mts:'
        print(f"Serpent2 inputs not written: {exc}")
    print("coverage problems:", study.check_references(**sel) or "none")

    if args.import_dragon:
        title("B.4a Importing existing DRAGON COMPOs (no re-run)")
        for item in args.import_dragon:
            key, path = item.split("=", 1)
            study.import_dragon_run(*key.split("/"), source=path)

    if args.run:
        title("B.4 DRAGON runs")
        print(study.run_dragon(num_threads=args.threads, **sel))

    if any(c.status == "completed" for c in study.cases(**sel)):
        title("B.5 Post-processing and figures")
        study.post_process(**sel)
        for f in study.figures(items=("keff", "macro_rates", "spectra")):
            print("   ", f.relative_to(study.directory))
    study.print_status()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--part", choices=["A", "B", "AB"], default="AB")
    ap.add_argument("--root", type=Path,
                    default=PROJECT_ROOT / "tutorials" / "results" / "code_to_code_study")
    ap.add_argument("--run", action="store_true", help="part B: run DRAGON")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--assembly", action="append", help="part B selection (repeatable)")
    ap.add_argument("--condition", action="append", help="part B selection (repeatable)")
    ap.add_argument("--import-serpent", action="store_true",
                    help="part B: register the existing NSE Serpent2 runs first")
    ap.add_argument("--import-dragon", action="append", metavar="ASM/VARIANT/COND=PATH",
                    help="part B: use an existing COMPO (or run folder) for a case")
    args = ap.parse_args()
    if "A" in args.part:
        part_a(args.root / "A_python_api")
    if "B" in args.part:
        part_b(args.root / "B_GE14_HZP_HFP", args)
