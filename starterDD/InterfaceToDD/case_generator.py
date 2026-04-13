# Case generator for CLE2000 procedure generation
# A DRAGON5 case is defined by a main procedure and its
# sub-procedures, and the module calls in the main procedure.
#
# R.Guasch — 10/03/2026

from .CLE2000 import (
    main_procedure, sub_procedure,
    validate_varname, wrap_cle2000_line,
)
from .dragon_module_calls import LIB, EDI_COMPO, TRK, EDI_condensation, SPH_correction, MIXEQ
from ..DDModel.DragonCalculationScheme import (
    DragonCalculationScheme,
    BRANCH_TYPE_TO_PARA_NAME,
)
from ..DDModel.helpers import associate_material_to_rod_ID
from ..MaterialProperties.material_mixture import (
    parse_all_compositions_from_yaml,
    compute_water_iso_densities_at_densities,
)
from ..GeometryAnalysis.tdt_parser import (
    read_material_mixture_indices_from_tdt_file,
)

import os

# Optional glow dependency — only needed when call_glow=True.
# glow requires the SALOME platform (typically run via Docker).
try:
    from ..GeometryBuilder.glow_builder import (
        build_full_assembly_geometry,
    )
    from glow.support.types import GeometryType, PropertyType
    GLOW_AVAILABLE = True
    
except ImportError:
    GLOW_AVAILABLE = False


class DragonCase:
    """
    Centralised case generator for DRAGON5/CLE2000 calculations.

    Given three YAML configuration files (materials, geometry,
    calculation scheme) and a library alias mapping the class:

    1. Builds a ``CartesianAssemblyModel`` from the GEOM + MATS
       yamls.
    2. Applies radial subdivision and fuel mixture numbering per
       the calculation scheme.
    3. Generates four CLE2000 files:

       - **<case>.x2m** — main procedure with USS/ASM/FLU calls.
       - **MIX.c2m** — material library (LIB:).
       - **TRK.c2m** — tracking (SALT:/MCCGT:).
       - **EDIR.c2m** — editions/COMPO export.
    4. Generates additional MIXEQ procedures if mix numbering strategy
        transitions are detected between consecutive trackable steps in the calculation scheme.
        These map indices from the current step's strategy to the next step's strategy, ensuring correct material-mixture-index tracking across steps with different numbering schemes.
        e.g. : use material-based numbering for self-shielding and first calculation level steps, use pin-based numbering for MOC flux calculation step. 
    """

    def __init__(self, case_name, call_glow,
                 draglib_name_to_alias, config_yamls,
                 enable_g2s=False, output_path=None,
                 tdt_path=None, tdt_base_name=None):
        """
        Parameters
        ----------
        case_name : str
            Identifier for the DRAGON case.
        call_glow : bool
            If ``True``, call glow to generate TDT geometry
            files.  If ``False`` the TDT files are expected
            to exist already.
        draglib_name_to_alias : dict
            ``{draglib_file_name: cle2000_alias}``.
        config_yamls : dict
            ``{"MATS": path, "GEOM": path,
              "CALC_SCHEME": path}``.
        enable_g2s : bool
            Include G2S: geometry visualization calls in x2m.
        output_path : str or None
            Directory for generated CLE2000 files.  Defaults
            to ``"./cle2000_procs/<case_name>"``.
        tdt_path : str or None
            Directory where TDT ``.dat`` files are located
            (or will be written when ``call_glow=True``).
            Defaults to ``output_path``.
        tdt_base_name : str or None
            Base name prefix for TDT file names.  The full
            TDT file name is built as
            ``<tdt_base_name>_<step>_<method>``.
            Defaults to ``case_name``.
        """
        self.case_name = case_name
        self.call_glow = call_glow
        self.draglib_name_to_alias = draglib_name_to_alias
        self.enable_g2s = enable_g2s
        self.tdt_base_name = tdt_base_name or case_name

        # Temperature dictionaries for material-specific temperatures
        self._fuel_material_temperatures = {}  # {material_name: temperature_K}
        self._non_fuel_temperatures = {}  # {material_type: temperature_K}

        # ----------------------------------------------------------
        # Resolve all paths to absolute so the case is
        # independent of later CWD changes.
        # ----------------------------------------------------------
        from pathlib import Path

        cwd = Path.cwd()

        # Config yamls — resolve each value
        self.config_yamls = {
            key: str(Path(val).resolve())
            if not os.path.isabs(val)
            else val
            for key, val in config_yamls.items()
        }

        # Output path (for generated CLE2000 procedures)
        raw_output = output_path or os.path.join(
            "cle2000_procs", case_name,
        )
        self.output_path = str(
            (cwd / raw_output).resolve()
            if not os.path.isabs(raw_output)
            else Path(raw_output).resolve()
        )

        # TDT path (where .dat files live / will be written)
        if tdt_path is not None:
            self.tdt_path = str(
                (cwd / tdt_path).resolve()
                if not os.path.isabs(tdt_path)
                else Path(tdt_path).resolve()
            )
        else:
            self.tdt_path = self.output_path

        self._validate_config_yamls()

    # -----------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------

    def _validate_config_yamls(self):
        expected = ("MATS", "GEOM", "CALC_SCHEME")
        for key in expected:
            if key not in self.config_yamls:
                raise ValueError(
                    f"Missing key '{key}' in config_yamls. "
                    f"Expected keys: {expected}"
                )
        self.mats_yaml = self.config_yamls["MATS"]
        self.geom_yaml = self.config_yamls["GEOM"]
        self.calc_scheme_yaml = self.config_yamls[
            "CALC_SCHEME"
        ]

    def set_fuel_material_temperatures(self, material_temps: dict):
        """
        Register fuel material temperatures by name.

        Parameters
        ----------
        material_temps : dict
            Mapping of fuel material name → temperature in K.
            Example: {'UOX_4.5': 900.0, 'MOX_8.0': 1100.0}
        """
        self._fuel_material_temperatures.update(material_temps)

    def set_non_fuel_temperatures(self, structural_temperature: float = 600.0,
                                   gap_temperature: float = 600.0,
                                   coolant_temperature: float = 600.0,
                                   moderator_temperature: float = 600.0):
        """
        Register non-fuel material temperatures (single value per type).

        These apply to all non-fuel materials of their type across the assembly.

        Parameters
        ----------
        structural_temperature : float
            Temperature for cladding/absorber tubes in K
        gap_temperature : float
            Temperature for gas gaps (helium, void) in K
        coolant_temperature : float
            Temperature for coolant (single temperature across assembly)
        moderator_temperature : float
            Temperature for moderator (single temperature across assembly)
        """
        self._non_fuel_temperatures.update({
            'structural': structural_temperature,
            'gap': gap_temperature,
            'coolant': coolant_temperature,
            'moderator': moderator_temperature,
        })

    # -----------------------------------------------------------
    # Assembly model construction
    # -----------------------------------------------------------

    def _extract_cle2000_temperatures(self):
        """
        Extract representative temperatures for CLE2000 MIX procedure call.

        Returns dict: {'TFUEL': float, 'TBOX': float, 'TCLAD': float, ...}

        Logic:
        - For fuel: collect temperatures from fuel mixtures
          (all UOX_4.5 instances share one temp, all MOX_8.0 share another, etc.)
        - For non-fuel: all instances of a type share one temp, so take first instance of each type
        """
        if self.assembly is None or self.assembly.compositions is None or len(self.assembly.compositions) == 0:
            # Fallback: use assembly defaults
            return {
                'TFUEL': self.assembly.default_fuel_temperature if self.assembly else 900.0,
                'TBOX': self.assembly.default_structural_temperature if self.assembly else 600.0,
                'TCLAD': self.assembly.default_structural_temperature if self.assembly else 600.0,
                'TCOOL': self.assembly.default_coolant_temperature if self.assembly else 600.0,
                'TMODE': self.assembly.default_moderator_temperature if self.assembly else 600.0,
                'TCTRL': self.assembly.default_structural_temperature if self.assembly else 600.0,
            }

        # Group mixtures by material type
        temps_by_type = {}  # {material_type: temp}
        fuels = []  # Collect fuel temps to extract representative value

        # Collect fuel material mixture temps
        if hasattr(self.assembly, 'fuel_material_mixtures') and self.assembly.fuel_material_mixtures:
            for mix in self.assembly.fuel_material_mixtures:
                if hasattr(mix, 'temperature') and mix.temperature is not None:
                    fuels.append(mix.temperature)

        # Collect non-fuel material mixture temps (first instance of each type)
        if hasattr(self.assembly, 'non_fuel_material_mixtures') and self.assembly.non_fuel_material_mixtures:
            for mix in self.assembly.non_fuel_material_mixtures:
                material_type = getattr(mix, 'material_type_key', None)
                if material_type and material_type not in temps_by_type:
                    if hasattr(mix, 'temperature') and mix.temperature is not None:
                        temps_by_type[material_type] = mix.temperature

        # Map to CLE2000 variables with fallback to assembly defaults
        return {
            'TFUEL': fuels[0] if fuels else self.assembly.default_fuel_temperature,
            'TBOX': temps_by_type.get('structural', self.assembly.default_structural_temperature),
            'TCLAD': temps_by_type.get('structural', self.assembly.default_structural_temperature),
            'TCOOL': temps_by_type.get('coolant', self.assembly.default_coolant_temperature),
            'TMODE': temps_by_type.get('moderator', self.assembly.default_moderator_temperature),
            'TCTRL': temps_by_type.get('structural', self.assembly.default_structural_temperature),
        }

    def _build_assembly_model(self):
        """Build and return a ``CartesianAssemblyModel``."""
        from ..DDModel.DragonModel import CartesianAssemblyModel

        compositions = parse_all_compositions_from_yaml(
            self.mats_yaml
        )
        rod_to_mat = associate_material_to_rod_ID(
            self.mats_yaml, self.geom_yaml,
        )

        assembly = CartesianAssemblyModel(
            name=self.case_name,
            tdt_file="",  # TDT will be generated later
            geometry_description_yaml=self.geom_yaml,
        )
        assembly.set_rod_ID_to_material_mapping(rod_to_mat)

        # Set assembly-level defaults (used as fallback if temperatures not registered via API)
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0,
            gap_temperature=600.0,
            coolant_temperature=600.0,
            moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True)

        # Pass temperature dicts to set_material_compositions
        fuel_temps = self._fuel_material_temperatures
        non_fuel_temps = self._non_fuel_temperatures
        assembly.set_material_compositions(
            compositions,
            fuel_temps=fuel_temps,
            non_fuel_temps=non_fuel_temps
        )
        return assembly, compositions

    def _generate_mixeq_procedures(self, scheme, assembly, draglib_alias):
        """
        Generate MIXEQ procedures for numbering strategy transitions.

        This method detects when consecutive trackable steps use different
        mix numbering strategies and generates MIXEQ procedures to map
        mixture indices between the different numbering schemes.

        Parameters
        ----------
        scheme : DragonCalculationScheme
            Calculation scheme with trackable steps
        assembly : CartesianAssemblyModel
            Assembly with mix state history from multiple steps
        draglib_alias : str
            Alias for the draglib to be used in MIXEQ procedures : used to recover depletion chain information.

        Returns
        -------
        list of tuples
            [(proc_name, file_path, current_step_name, next_step_name), ...] for generated MIXEQ procedures
        """
        from .CLE2000 import validate_varname

        mixeq_procedures = []
        trackable_steps = scheme.get_trackable_steps()

        # Check consecutive steps for strategy transitions
        for i in range(len(trackable_steps) - 1):
            current_step = trackable_steps[i]
            next_step = trackable_steps[i + 1]

            current_strategy = current_step.mix_numbering_strategy
            next_strategy = next_step.mix_numbering_strategy

            if current_strategy != next_strategy:
                print(f"[MIXEQ] Strategy transition detected: {current_strategy} → {next_strategy}")
                print(f"[MIXEQ] Generating MIXEQ for {current_step.name} → {next_step.name}")

                try:

                    input_lib_name = f"LIBRARY2"  # Default input library name for MIXEQ
                    output_lib_name = f"LIBEQ"
                    # Generate MIXEQ procedure
                    mixeq = MIXEQ(assembly, input_lib_name, output_lib_name, current_step.name, next_step.name, draglib_alias)
                    proc_name = f"MIXEQ_{current_step.name}_to_{next_step.name}"

                    # Validate procedure name (CLE-2000 limitation)
                    validate_varname(proc_name[:12])
                    mixeq_path = mixeq.write_to_c2m(self.output_path, proc_name[:12])
                    mixeq_procedures.append((proc_name[:12], mixeq_path, current_step.name, next_step.name))

                    print(f"[MIXEQ] Successfully generated {proc_name[:12]}")

                except Exception as e:
                    print(f"[ERROR] Failed to generate MIXEQ for {current_step.name} → {next_step.name}: {e}")
                    # Continue processing other transitions rather than failing completely
                    continue

        if not mixeq_procedures:
            print("[MIXEQ] No strategy transitions detected - no MIXEQ procedures generated")
        else:
            print(f"[MIXEQ] Generated {len(mixeq_procedures)} MIXEQ procedure(s)")

        return mixeq_procedures

    def _validate_mixeq_generation(self, scheme, assembly):
        """
        Validate that MIXEQ generation prerequisites are met.

        Parameters
        ----------
        scheme : DragonCalculationScheme
            Calculation scheme to validate
        assembly : CartesianAssemblyModel
            Assembly model to validate

        Returns
        -------
        bool
            True if validation passes

        Raises
        ------
        RuntimeError
            If validation fails
        """
        # Check that assembly supports MIXEQ generation
        if not hasattr(assembly, 'build_mixeq_correspondence_table'):
            raise RuntimeError(
                "Assembly model does not support MIXEQ generation. "
                "The build_mixeq_correspondence_table method is missing."
            )

        # Check that assembly has mix state history tracking
        if not hasattr(assembly, 'mix_state_history'):
            raise RuntimeError(
                "Assembly model does not have mix_state_history tracking enabled. "
                "This is required for MIXEQ generation."
            )

        # Validate supported strategy transitions
        trackable_steps = scheme.get_trackable_steps()
        valid_strategies = ["by_material", "by_pin"]

        for step in trackable_steps:
            if step.mix_numbering_strategy not in valid_strategies:
                raise RuntimeError(
                    f"Unsupported mix numbering strategy '{step.mix_numbering_strategy}' "
                    f"in step '{step.name}'. Supported strategies: {valid_strategies}"
                )

        # Check for potential strategy transition issues
        strategy_changes = 0
        for i in range(len(trackable_steps) - 1):
            current_step = trackable_steps[i]
            next_step = trackable_steps[i + 1]

            if current_step.mix_numbering_strategy != next_step.mix_numbering_strategy:
                strategy_changes += 1

                # Validate transition directions (optional warning)
                if (current_step.mix_numbering_strategy == "by_pin" and
                    next_step.mix_numbering_strategy == "by_material"):
                    print(f"[WARNING] Unusual strategy transition {current_step.name}: "
                          f"by_pin → by_material. This may indicate a configuration issue.")

        if strategy_changes > 3:
            print(f"[WARNING] High number of strategy transitions detected ({strategy_changes}). "
                  f"This may result in many MIXEQ procedures being generated.")

        return True

    # -----------------------------------------------------------
    # Case execution
    # -----------------------------------------------------------

    def run(
        self,
        dragon_executable=None,
        draglib_paths=None,
        results_root="./results",
        num_threads=1,
        create_latest_symlink=True,
        dry_run=False,
    ):
        """Generate CLE2000 procedures and execute the Dragon case.

        This is a convenience method that chains procedure
        generation with execution via ``DragonRunner``.

        Parameters
        ----------
        dragon_executable : str or None
            Path to the Dragon binary.  If ``None``, resolved
            from ``$dragon_exec`` or ``$DRAGON_EXEC``.
        draglib_paths : dict or None
            ``{draglib_name: absolute_path}``.  If ``None``,
            resolved from the ``$DRAGLIBS`` directory.
        results_root : str
            Root directory for result storage.
        num_threads : int
            OpenMP thread count.
        create_latest_symlink : bool
            Maintain a ``latest`` symlink to the most recent run.
        dry_run : bool
            If ``True``, stage inputs and write manifest but
            do not execute Dragon.

        Returns
        -------
        RunResult
        """
        from .dragon_runner import DragonRunner

        # Generate procedures if not already done
        if not hasattr(self, "scheme") or self.scheme is None:
            self.generate_cle2000_procedures()
        runner = DragonRunner(
            dragon_case=self,
            dragon_executable=dragon_executable,
            draglib_paths=draglib_paths,
            results_root=results_root,
            num_threads=num_threads,
            create_latest_symlink=create_latest_symlink,
        )
        return runner.run(dry_run=dry_run)

    # -----------------------------------------------------------
    # Main entry point
    # -----------------------------------------------------------

    def generate_cle2000_procedures(self):
        """Generate all CLE2000 files for the case.

        Returns
        -------
        dict
            ``{"x2m": path, "mix": path,
              "trk": path, "edir": path}``
        """
        # Determine the first draglib alias
        draglib_alias = list(
            self.draglib_name_to_alias.values()
        )[0]
        # 1. Load scheme
        scheme = DragonCalculationScheme.from_yaml(
            self.calc_scheme_yaml
        )
        self.scheme = scheme

        # 2. Build assembly model
        assembly, compositions = self._build_assembly_model()
        self.assembly = assembly

        # 3. Apply the SSH step radii & number fuel mixes
        ssh_steps = scheme.get_self_shielding_steps()
        if not ssh_steps:
            raise RuntimeError(
                "Scheme must have at least one "
                "self_shielding step."
            )
        ssh_step = ssh_steps[0]
        ssh_step.apply_radii(assembly)

        # Apply initial mix numbering strategy from SSH step to establish baseline
        # This ensures the assembly has a consistent starting state before per-step processing
        #assembly.apply_mix_numbering_strategy(ssh_step.mix_numbering_strategy)

        # 4. Build geometry (if call_glow) and read TDT files
        #    for each calculation step.  The TDT file provides
        #    the definitive material-mixture-index numbering
        #    (assigned by glow/SALOME) for both fuel and
        #    non-fuel materials.
        for step in scheme.get_trackable_steps():
            # Apply mix numbering strategy for THIS step BEFORE building geometry
            # (Phase 3: hybrid orchestration per-step strategy)
            # This ensures glow/SALOME uses the correct naming scheme when
            # generating geometry for this particular calculation step.
            assembly.apply_mix_numbering_strategy(step.mix_numbering_strategy)
            
            if step.tdt_file_id is not None:
                tdt_file_name = (
                    f"{self.tdt_base_name}_{step.tdt_file_id}"
                    f"_{step.spatial_method}"
            )
            else:
                tdt_file_name = (
                    f"{self.tdt_base_name}_{step.name}"
                    f"_{step.spatial_method}"
                )

            # 4a. Optionally build geometry with glow
            if self.call_glow:
                if not GLOW_AVAILABLE:
                    raise ImportError(
                        "call_glow=True but glow is not "
                        "importable.  glow requires the "
                        "SALOME platform (typically run "
                        "via Docker)."
                    )
                lattice, assembly_box = build_full_assembly_geometry(
                    assembly_model=assembly,
                    calculation_step=step,
                    output_path=self.tdt_path,
                    output_file_name=tdt_file_name,
                )
                lattice.show(
                    geometry_type_to_show=GeometryType.SECTORIZED,
                    property_type_to_show=PropertyType.MATERIAL,
                )
                if step.export_macros:
                    lattice.show(
                        geometry_type_to_show=GeometryType.SECTORIZED,
                        property_type_to_show=PropertyType.MACRO,
                    )
                    lattice.show(
                        geometry_type_to_show=GeometryType.SECTORIZED,
                        property_type_to_show=PropertyType.MATERIAL,
                    )

            # 4b. Read TDT and enforce material indices for THIS step
            #     This records the TDT-assigned indices in the step's state history
            #     for per-step mix tracking (Phase 3: hybrid orchestration).
            #     Note: SSH step defines the LIB numbering used by all subsequent steps.
            tdt_file_path = self.tdt_path
            tdt_full = os.path.join(
                tdt_file_path,
                f"{tdt_file_name}_{step.tracking}"
                + ("_MACRO" if step.export_macros else "")
                + ".dat",
            )
            if not os.path.isfile(tdt_full):
                raise FileNotFoundError(
                    f"TDT file not found: {tdt_full}\n"
                    f"  tdt_path = {tdt_file_path}\n"
                    f"When call_glow=False, the TDT .dat "
                    f"file must already exist.  Pass an "
                    f"explicit tdt_path= to DragonCase "
                    f"pointing at the directory that "
                    f"contains the file."
                )
            tdt_indices = (
                read_material_mixture_indices_from_tdt_file(
                    tdt_file_path=tdt_file_path,
                    tdt_file_name=tdt_file_name,
                    tracking_option=step.tracking,
                    include_macros=step.export_macros,
                    material_names=None,
                )
            )
            # Enforce TDT indices: updates assembly mix indices and records
            # them in mix_state_history for this step
            assembly.enforce_material_mixture_indices_from_tdt(
                tdt_indices, step_name=step.name
            )
            # Identify generating / daughter mixes
            #    (must happen after TDT enforcement so that
            #    indices are final, and before LIB creation)
            assembly.identify_generating_and_daughter_mixes()
            
            if step.step_type == "self_shielding":
                print(f"[SSH Step] Applied TDT indices for step '{step.name}'. ")
                print(f"These indices define the initial LIB numbering.")
                
                # ----- MIX.c2m -----
                mix_proc_name = "MIX"
                print(f"Generating MIX procedure '{mix_proc_name}' with LIB: ...")
                has_density = (
                    scheme.has_branches()
                    and scheme.get_branch("coolant_density") is not None
                )
                lib = LIB(assembly, density_branch=has_density, ssh_calculation_step=ssh_step)
                mix_path = lib.write_to_c2m(
                    self.output_path, mix_proc_name,
                )


        # Build per-step mix mapping
        # This records each trackable step's mix names and TDT indices for later use
        step_mix_mapping = {}
        for step in scheme.get_trackable_steps():
            names, indices = assembly.get_step_mix_names_and_indices(step.name)
            if names is not None:
                state = assembly.get_step_mix_state(step.name)
                step_mix_mapping[step.name] = {
                    "strategy": state["strategy"] if state else None,
                    "mix_names": names,
                    "tdt_indices": indices,
                }

        # 6. Generate MIXEQ procedures for numbering strategy transitions
        try:
            # Validate prerequisites for MIXEQ generation
            self._validate_mixeq_generation(scheme, assembly)
            mixeq_procedures = self._generate_mixeq_procedures(scheme, assembly, draglib_alias)
        except Exception as e:
            print(f"[ERROR] MIXEQ validation failed: {e}")
            print("[MIXEQ] Skipping MIXEQ generation - continuing with standard workflow")
            mixeq_procedures = []

        # ----- TRK.c2m -----
        trk_proc_name = f"TRK"
        trk = TRK(scheme, self.case_name)
        trk_path = trk.write_to_c2m(
            self.output_path, trk_proc_name,
        )

        # ----- EDIR.c2m (EDI/COMPO) -----
        edir_proc_name = f"EDIR"
        if scheme.has_branches():
            edi_compo = EDI_COMPO(
                assembly,
                branches=scheme.branches,
                outputs=scheme.outputs,
            )
        else:
            edi_compo = EDI_COMPO(assembly)
            # Default edition: homogenised condensed
            edi_compo.add_edition(
                name="EDIHOM_COND",
                comment="Condensed, homogenized over fuel",
                isotopes=[
                    "U235", "U238", "U234",
                    "Gd155", "Gd157",
                ],
                spatial_mode="FUEL",
                energy_bounds=[],
            )
        edir_path = edi_compo.write_to_c2m(
            self.output_path, edir_proc_name,
        )

        # ----- main x2m -----
        if scheme.is_two_level_scheme():
            if scheme.has_branches():
                x2m_path = self._build_main_x2m_with_branches_2level(
                    scheme, trk, lib, edi_compo,
                    mix_proc_name, trk_proc_name,
                    edir_proc_name, mixeq_procedures, 
                    draglib_alias,
                    ssh_step,
                )
            else:
                x2m_path = self._build_main_x2m_2level(
                    scheme, trk, lib, edi_compo,
                    mix_proc_name, trk_proc_name,
                    edir_proc_name, mixeq_procedures,  
                    draglib_alias,
                    ssh_step,
                )
        elif scheme.has_branches():
            x2m_path = self._build_main_x2m_with_branches(
                scheme, trk, lib, edi_compo,
                mix_proc_name, trk_proc_name,
                edir_proc_name, mixeq_procedures,  
                draglib_alias,
                ssh_step,
            )
        else:
            x2m_path = self._build_main_x2m(
                scheme, trk, lib, edi_compo,
                mix_proc_name, trk_proc_name,
                edir_proc_name, mixeq_procedures, 
                draglib_alias,
                ssh_step,
            )

        # Build return dictionary with individual MIXEQ procedure entries
        result_dict = {
            "x2m": x2m_path,
            "mix": mix_path,
            "trk": trk_path,
            "edir": edir_path,
        }

        # Add MIXEQ procedures as individual entries
        for i, (proc_name, proc_path, current_step_name, next_step_name) in enumerate(mixeq_procedures):
            result_dict[f"mixeq_{i+1}"] = proc_path

        # Store procedures dictionary for dragon_runner
        self._procedure_files = result_dict

        return result_dict

    # -----------------------------------------------------------
    # Main x2m builder (single statepoint, backward compat)
    # -----------------------------------------------------------

    def _build_main_x2m(self, scheme, trk, lib, edi_compo,
                        mix_proc_name, trk_proc_name,
                        edir_proc_name, mixeq_procedures, 
                        draglib_alias,
                        ssh_step):
        """Assemble and write the main .x2m procedure."""

        proc = main_procedure(self.case_name)
        proc.enable_g2s = self.enable_g2s

        # --- Modules ---
        modules = [
            "LIB", "G2S", "SALT", "MCCGT",
            "USS", "ASM", "FLU", "EDI",
            "SPH", "GREP", "DELETE", "END",
        ]
        for m in modules:
            proc.add_module_call(m)

        # --- Sub-procedures ---
        sub_mix = sub_procedure(mix_proc_name)
        sub_trk = sub_procedure(trk_proc_name)
        sub_edir = sub_procedure(edir_proc_name)
        for sp in (sub_mix, sub_trk, sub_edir):
            proc.add_sub_procedure(sp)
        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                sub_mixeq = sub_procedure(proc_name)
                proc.add_sub_procedure(sub_mixeq)

        # --- Data structures ---
        proc.add_linked_list("LIBRARY")
        proc.add_linked_list("LIBRARY2")
        for trk_ll, _ in trk.get_track_names():
            proc.add_linked_list(trk_ll)
        proc.add_linked_list("SYS")
        proc.add_linked_list("FLUX")
        proc.add_linked_list("COMPO")

        for _, trkfil in trk.get_track_names():
            proc.add_seq_binary(trkfil)

        # --- Variables ---
        aniso = ssh_step.anisotropy_level
        proc.add_variable("o_anis", "INTEGER", aniso)
        proc.add_variable(
            "name_compo", "STRING",
            f"_CPO_{self.case_name}",
        )
        proc.add_variable(
            "Draglib", "STRING", draglib_alias,
        )

        # Temperature variables
        temps = self._extract_cle2000_temperatures()
        for tname, temp_val in temps.items():
            proc.add_variable(tname, "REAL", temp_val)

        # --- TDT file SEQ_ASCII imports ---
        tdt_map = trk.get_tdt_file_mapping()
        for var, fname in tdt_map.items():
            proc.add_seq_ascii(var, fname)

        # --- G2S block ---
        if proc.enable_g2s:
            proc.add_body_block(
                proc.build_g2s_block(tdt_map)
            )

        # --- TRK sub-procedure call ---
        track_names = trk.get_track_names()
        tdt_vars = list(tdt_map.keys())
        lhs = " ".join(
            f"{ll} {bf}" for ll, bf in track_names
        )
        rhs_tdts = " ".join(tdt_vars)

        trk_call = wrap_cle2000_line(
            f"{lhs} := {trk_proc_name} "
            f"{rhs_tdts} ;"
        )
        proc.add_body_line(
            "*" * 50
        )
        proc.add_body_line(
            "* Tracking sub-procedure call"
        )
        proc.add_body_line(
            "*" * 50
        )
        proc.add_body_line(trk_call)
        proc.add_body_line("")

        # --- MIX sub-procedure call ---
        ssh_method = ssh_step.self_shielding_method or "PT"
        if hasattr(ssh_step, 'transport_correction') and ssh_step.transport_correction is not None:
            transport_correction = ssh_step.transport_correction
        else:
            transport_correction = "NONE"
        mix_call_line = (
            f"LIBRARY := {mix_proc_name} :: "
            f"<<Draglib>> '{ssh_method}' <<o_anis>> "
            f"'{transport_correction}' "
            f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
            f"<<TCOOL>> <<TMODE>> <<TCTRL>> ;"
        )
        proc.add_body_line(
            "*" * 50
        )
        proc.add_body_line("* LIBRARY creation")
        proc.add_body_line(
            "*" * 50
        )
        proc.add_body_line(
            wrap_cle2000_line(mix_call_line)
        )
        proc.add_body_line("")

        # --- USS call ---
        ssh_trk, ssh_trkfil = track_names[0]
        proc.add_body_line(
            "*" * 50
        )
        proc.add_body_line("* USS: self-shielding")
        proc.add_body_line(
            "*" * 50
        )
        uss_line = (
            f"LIBRARY2 := USS: LIBRARY "
            f"{ssh_trk} {ssh_trkfil} ::"
        )
        proc.add_body_line(
            wrap_cle2000_line(uss_line)
        )
        proc.add_body_line(" PASS 3")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == ssh_step.name and next_step_name == scheme.get_flux_steps()[0].name:
                    proc.add_body_line(
                        f"* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"LIBRARY2 := {proc_name} LIBRARY2 ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(mixeq_call)
                    )

        asm_keyword = "PIJ" if ssh_step.spatial_method == "CP" else "ARM"
        proc.add_body_line(f"    EDIT 1 {asm_keyword}")
        proc.add_body_line(";")
        proc.add_body_line("")

        # --- ASM + FLU on flux step ---
        flux_steps = scheme.get_flux_steps()
        if flux_steps:
            flux_step = flux_steps[0]
            flux_trk, flux_trkfil = track_names[-1]

            proc.add_body_line(
                "*" * 50
            )
            proc.add_body_line(
                "* ASM: + FLU: flux calculation"
            )
            proc.add_body_line(
                "*" * 50
            )
            asm_line = (
                f"SYS := ASM: LIBRARY2 "
                f"{flux_trk} {flux_trkfil} ::"
            )
            proc.add_body_line(
                wrap_cle2000_line(asm_line)
            )
            asm_keyword = "PIJ" if flux_step.spatial_method == "CP" else "ARM"
            proc.add_body_line(f"    EDIT 1 {asm_keyword}")
            proc.add_body_line(";")
            proc.add_body_line("")

            flu_line = (
                f"FLUX := FLU: SYS LIBRARY2 "
                f"{flux_trk} {flux_trkfil} ::"
            )
            proc.add_body_line(
                wrap_cle2000_line(flu_line)
            )
            proc.add_body_line("    EDIT 1 TYPE K")
            proc.add_body_line(";")

            # GREP keff
            keff_var = "keff"
            validate_varname(keff_var)
            proc.add_variable(keff_var, "REAL")
            proc.add_body_line(
                wrap_cle2000_line(
                    "GREP: FLUX :: GETVAL "
                    "'K-EFFECTIVE' 1 1 1 "
                    f">>keff<< ;"
                )
            )
            proc.add_body_line(
                f'ECHO "{self.case_name} keff=" keff ;'
            )
            proc.add_body_line("")

            # --- EDIR sub-procedure call ---
            proc.add_body_line(
                "*" * 50
            )
            proc.add_body_line(
                "* EDI/COMPO sub-procedure call"
            )
            proc.add_body_line(
                "*" * 50
            )
            edir_call = (
                f"{edir_proc_name} FLUX LIBRARY2 "
                f"{flux_trk} :: <<name_compo>> ;"
            )
            proc.add_body_line(
                wrap_cle2000_line(edir_call)
            )

        return proc.write_to_x2m(self.output_path)

    # -----------------------------------------------------------
    # Main x2m builder with branch loops (multi-statepoint)
    # -----------------------------------------------------------

    def _build_main_x2m_with_branches(
        self, scheme, trk, lib, edi_compo,
        mix_proc_name, trk_proc_name,
        edir_proc_name, mixeq_procedures, 
        draglib_alias,
        ssh_step,
    ):
        """Assemble and write the main .x2m with WHILE loops over branches.

        This generates the AT10-style multi-physics iteration pattern:
        UTL: parameter arrays → COMPO init → nested WHILE loops →
        LIB → USS → ASM → FLU → EDIR → cleanup.
        """
        from .CLE2000 import wrap_cle2000_line

        proc = main_procedure(self.case_name)
        proc.enable_g2s = self.enable_g2s

        # --- Modules ---
        modules = [
            "LIB", "G2S", "SALT", "MCCGT",
            "USS", "ASM", "FLU", "EDI",
            "SPH", "GREP", "UTL", "COMPO",
            "DELETE", "END",
        ]
        for m in modules:
            proc.add_module_call(m)

        # --- Sub-procedures ---
        print(f"Generating main x2m with branches. Sub-procedures: "
                f"{mix_proc_name}, {trk_proc_name}, "
                f"{edir_proc_name}")
        sub_mix = sub_procedure(mix_proc_name)
        sub_trk = sub_procedure(trk_proc_name)
        sub_edir = sub_procedure(edir_proc_name)
        for sp in (sub_mix, sub_trk, sub_edir):
            proc.add_sub_procedure(sp)

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                sub_mixeq = sub_procedure(proc_name)
                proc.add_sub_procedure(sub_mixeq)

        # --- Data structures ---
        proc.add_linked_list("LIBRARY")
        proc.add_linked_list("LIBRARY2")
        for trk_ll, _ in trk.get_track_names():
            proc.add_linked_list(trk_ll)
        proc.add_linked_list("SYS")
        proc.add_linked_list("FLUX")
        proc.add_linked_list("COMPO")
        proc.add_linked_list("PARAMS")

        for _, trkfil in trk.get_track_names():
            proc.add_seq_binary(trkfil)

        # --- Variables ---
        aniso = ssh_step.anisotropy_level
        proc.add_variable("o_anis", "INTEGER", aniso)
        proc.add_variable(
            "name_compo", "STRING",
            f"_CPO_{self.case_name}",
        )
        proc.add_variable(
            "Draglib", "STRING", draglib_alias,
        )

        # Ordered branches for loop generation
        ordered_branches = scheme.get_ordered_branches()

        # Branch count and counter variables
        for branch in ordered_branches:
            proc.add_variable(
                branch.count_var, "INTEGER",
                len(branch.values),
            )
            proc.add_variable(
                branch.counter_var, "INTEGER", 0,
            )

        # Branch parameter REAL variables (current values)
        for branch in ordered_branches:
            proc.add_variable(
                branch.para_name, "REAL",
            )

        # Coolant density isotopic density variables
        density_branch = scheme.get_branch("coolant_density")
        if density_branch:
            proc.add_variable("N_H", "REAL")
            proc.add_variable("N_O", "REAL")

        # keff
        proc.add_variable("keff", "REAL")

        # --- TDT file SEQ_ASCII imports ---
        tdt_map = trk.get_tdt_file_mapping()
        for var, fname in tdt_map.items():
            proc.add_seq_ascii(var, fname)

        # COMPO export file
        proc.add_seq_ascii(
            "_COMPO",
            f"_CPO_{self.case_name}",
        )

        # --- G2S block ---
        if proc.enable_g2s:
            proc.add_body_block(
                proc.build_g2s_block(tdt_map)
            )

        # --- TRK sub-procedure call ---
        track_names = trk.get_track_names()
        tdt_vars = list(tdt_map.keys())
        lhs = " ".join(
            f"{ll} {bf}" for ll, bf in track_names
        )
        rhs_tdts = " ".join(tdt_vars)

        trk_call = wrap_cle2000_line(
            f"{lhs} := {trk_proc_name} "
            f"{rhs_tdts} ;"
        )
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Tracking sub-procedure call")
        proc.add_body_line("*" * 50)
        proc.add_body_line(trk_call)
        proc.add_body_line("")

        # ===================================================
        # UTL: parameter array declaration
        # ===================================================
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            "* Parameter arrays for branch calculations"
        )
        proc.add_body_line("*" * 50)

        utl_block = "PARAMS := UTL: ::\n"
        for branch in ordered_branches:
            vals_str = " ".join(
                f"{v}" for v in branch.values
            )
            utl_block += (
                f"    CREA\n"
                f"        {branch.array_name} "
                f"<<{branch.count_var}>> =\n"
                f"            {vals_str}\n"
            )
        utl_block += ";\n"
        proc.add_body_block(utl_block)

        # Pre-compute N_H and N_O arrays for coolant density branch
        if density_branch:
            iso_densities = compute_water_iso_densities_at_densities(
                density_branch.values
            )
            nh_vals = " ".join(
                f"{d['H1']:.8E}" for d in iso_densities
            )
            no_vals = " ".join(
                f"{d['O16']:.8E}" for d in iso_densities
            )
            n_dens = len(density_branch.values)
            proc.add_variable("n_dens", "INTEGER", n_dens)

            nh_utl = (
                f"PARAMS := UTL: PARAMS ::\n"
                f"    CREA\n"
                f"        N_H_arr <<n_dens>> =\n"
                f"            {nh_vals}\n"
                f"    CREA\n"
                f"        N_O_arr <<n_dens>> =\n"
                f"            {no_vals}\n"
                f";\n"
            )
            proc.add_body_block(nh_utl)

        proc.add_body_line("")

        # ===================================================
        # COMPO initialization (outside the loops)
        # ===================================================
        proc.add_body_line("*" * 50)
        proc.add_body_line("* COMPO initialization")
        proc.add_body_line("*" * 50)
        proc.add_body_block(
            edi_compo.compo.build_compo_init()
        )
        proc.add_body_line("")

        # ===================================================
        # Nested WHILE loops
        # ===================================================
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            "* Main loop over branch statepoints"
        )
        proc.add_body_line("*" * 50)

        # Open loops (outermost first)
        for depth, branch in enumerate(ordered_branches):
            indent = "    " * depth
            proc.add_body_line(
                f"{indent}WHILE {branch.counter_var} "
                f"{branch.count_var} < DO"
            )
            proc.add_body_line(
                f"{indent}    EVALUATE {branch.counter_var} "
                f":= {branch.counter_var} 1 + ;"
            )
            # GREP parameter value
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{indent}    GREP: PARAMS :: "
                    f"GETVAL '{branch.array_name}' "
                    f"<<{branch.counter_var}>> "
                    f">>{branch.para_name}<< ;"
                )
            )
            # If this is the coolant density branch, also GREP N_H, N_O
            if branch.type == "coolant_density":
                proc.add_body_line(
                    wrap_cle2000_line(
                        f"{indent}    GREP: PARAMS :: "
                        f"GETVAL 'N_H_arr' "
                        f"<<{branch.counter_var}>> "
                        f">>N_H<< ;"
                    )
                )
                proc.add_body_line(
                    wrap_cle2000_line(
                        f"{indent}    GREP: PARAMS :: "
                        f"GETVAL 'N_O_arr' "
                        f"<<{branch.counter_var}>> "
                        f">>N_O<< ;"
                    )
                )

        # Inner body (deepest indentation)
        inner_indent = "    " * len(ordered_branches)

        # Echo current statepoint
        echo_parts = " ".join(
            f'"{b.para_name}=" {b.para_name}'
            for b in ordered_branches
        )
        proc.add_body_line(
            wrap_cle2000_line(
                f'{inner_indent}ECHO "State Point:" ; \n'
            )
        )
        proc.add_body_line(
            wrap_cle2000_line(
                f'{inner_indent}ECHO {echo_parts} ; \n'
            )
        )

        # --- Map branch values to MIX temperature variables ---
        # The branch loop GREPs into TFuel / TCool / DCool;
        # the MIX procedure expects TFUEL / TCOOL / TMODE.
        if scheme.get_branch("fuel_temperature"):
            proc.add_body_line(
                f"{inner_indent}EVALUATE TFUEL "
                f":= TFuel ;"
            )
        if scheme.get_branch("coolant_temperature"):
            proc.add_body_line(
                f"{inner_indent}EVALUATE TCOOL "
                f":= TCool ;"
            )
            proc.add_body_line(
                f"{inner_indent}EVALUATE TMODE "
                f":= TCool ;"
            )
        proc.add_body_line("")

        # --- MIX sub-procedure call ---
        ssh_method = ssh_step.self_shielding_method or "PT"
        if hasattr(ssh_step, 'transport_correction') and ssh_step.transport_correction is not None:
            transport_correction = ssh_step.transport_correction
        else:
            transport_correction = "NONE"

        # Build MIX call — pass N_H N_O if density branch exists
        if density_branch:
            mix_call_line = (
                f"LIBRARY := {mix_proc_name} :: "
                f"<<Draglib>> '{ssh_method}' <<o_anis>> "
                f"'{transport_correction}' "
                f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
                f"<<TCOOL>> <<TMODE>> <<TCTRL>> "
                f"<<N_H>> <<N_O>> ;"
            )
        else:
            mix_call_line = (
                f"LIBRARY := {mix_proc_name} :: "
                f"<<Draglib>> '{ssh_method}' <<o_anis>> "
                f"'{transport_correction}' "
                f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
                f"<<TCOOL>> <<TMODE>> <<TCTRL>> ;"
            )
        proc.add_body_line(
            f"{inner_indent}* LIB: creation"
        )
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}{mix_call_line}"
            )
        )
        proc.add_body_line("")

        # --- USS call ---
        ssh_trk, ssh_trkfil = track_names[0]
        proc.add_body_line(
            f"{inner_indent}* USS: self-shielding"
        )
        uss_line = (
            f"LIBRARY2 := USS: LIBRARY "
            f"{ssh_trk} {ssh_trkfil} ::"
        )
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}{uss_line}"
            )
        )
        ssh_asm_keyword = "PIJ" if ssh_step.spatial_method == "CP" else "ARM"
        proc.add_body_line(
            f"{inner_indent}    EDIT 1 {ssh_asm_keyword}"
        )
        proc.add_body_line("    PASS 3")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line("")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == ssh_step.name and next_step_name == scheme.get_flux_steps()[0].name:
                    proc.add_body_line(
                        f"{inner_indent}* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"LIBRARY2 := {proc_name} LIBRARY2 ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(
                            f"{inner_indent}{mixeq_call}"
                        )
                    )
                    proc.add_body_line("")

        # --- ASM + FLU on flux step ---
        flux_steps = scheme.get_flux_steps()
        if flux_steps:
            flux_trk, flux_trkfil = track_names[-1]

            proc.add_body_line(
                f"{inner_indent}* ASM: + FLU:"
            )
            asm_line = (
                f"SYS := ASM: LIBRARY2 "
                f"{flux_trk} {flux_trkfil} ::"
            )
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{inner_indent}{asm_line}"
                )
            )
            asm_keyword = "PIJ" if flux_steps[0].spatial_method == "CP" else "ARM"
            proc.add_body_line(
                f"{inner_indent}    EDIT 1 {asm_keyword}"
            )
            proc.add_body_line(f"{inner_indent};")
            proc.add_body_line("")

            flu_line = (
                f"FLUX := FLU: SYS LIBRARY2 "
                f"{flux_trk} {flux_trkfil} ::"
            )
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{inner_indent}{flu_line}"
                )
            )
            proc.add_body_line(
                f"{inner_indent}    EDIT 1 TYPE K"
            )
            proc.add_body_line(f"{inner_indent};")

            # GREP keff
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{inner_indent}GREP: FLUX :: GETVAL "
                    f"'K-EFFECTIVE' 1 1 1 "
                    f">>keff<< ;"
                )
            )
            proc.add_body_line(
                f'{inner_indent}ECHO '
                f'"{self.case_name} keff=" keff ;'
            )
            proc.add_body_line("")

            # --- EDIR sub-procedure call ---
            proc.add_body_line(
                f"{inner_indent}* EDI/COMPO call"
            )
            # Build the EDIR call with branch params
            # COMPO is passed as I/O data structure
            para_args = " ".join(
                f"<<{b.para_name}>>"
                for b in ordered_branches
            )
            edir_call = (
                f"COMPO := {edir_proc_name} FLUX LIBRARY2 "
                f"{flux_trk} COMPO :: <<name_compo>> "
                f"{para_args} ;"
            )
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{inner_indent}{edir_call}"
                )
            )
            proc.add_body_line("")

            # --- Cleanup ---
            proc.add_body_line(
                f"{inner_indent}* Cleanup"
            )
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{inner_indent}LIBRARY LIBRARY2 := "
                    f"DELETE: LIBRARY LIBRARY2 ;"
                )
            )
            proc.add_body_line(
                f"{inner_indent}FLUX := DELETE: FLUX ;"
            )
            proc.add_body_line(
                f"{inner_indent}SYS := DELETE: SYS ;"
            )

        # Close loops (innermost first → outermost last)
        for depth in range(
            len(ordered_branches) - 1, -1, -1
        ):
            indent = "    " * depth
            branch = ordered_branches[depth]
            proc.add_body_line(
                f"{indent}ENDWHILE ;"
            )
            # Reset counter for next iteration of outer loop
            if depth > 0:
                proc.add_body_line(
                    f"{indent}EVALUATE "
                    f"{branch.counter_var} := 0 ;"
                )

        proc.add_body_line("")

        # ===================================================
        # COMPO export (after all loops)
        # ===================================================
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Export COMPO to file")
        proc.add_body_line("*" * 50)
        proc.add_body_line("_COMPO := COMPO ;")

        # Extract temperatures for non-fuel materials
        temps = self._extract_cle2000_temperatures()

        # Map temperature branch values to assembly temps
        # for the non-loop temperature variables (structural etc.)
        # For branch-driven temps, use branch values as override
        fuel_temp_branch = scheme.get_branch(
            "fuel_temperature"
        )
        cool_temp_branch = scheme.get_branch(
            "coolant_temperature"
        )
        proc.add_variable(
            "TFUEL", "REAL",
            fuel_temp_branch.values[0]
            if fuel_temp_branch else temps.get('TFUEL', 900.0),
        )
        proc.add_variable(
            "TCOOL", "REAL",
            cool_temp_branch.values[0]
            if cool_temp_branch else temps.get('TCOOL', 600.0),
        )
        proc.add_variable(
            "TMODE", "REAL",
            cool_temp_branch.values[0]
            if cool_temp_branch else temps.get('TMODE', 600.0),
        )
        # Non-fuel temperatures from extracted values
        proc.add_variable("TBOX", "REAL", temps.get('TBOX', 600.0))
        proc.add_variable("TCLAD", "REAL", temps.get('TCLAD', 600.0))
        proc.add_variable("TCTRL", "REAL", temps.get('TCTRL', 600.0))

        return proc.write_to_x2m(self.output_path)

    # -----------------------------------------------------------
    # 2-level x2m builders
    # -----------------------------------------------------------

    def _build_main_x2m_2level(self, scheme, trk, lib, edi_compo,
                               mix_proc_name, trk_proc_name,
                               edir_proc_name, mixeq_procedures, 
                               draglib_alias,
                               ssh_step):
        """Build the main .x2m for a 2-level calculation scheme
        (no branch loops, single statepoint).

        Sequence: TRK → MIX → USS → DELETE_SSH →
                  ASM_L1 → FLU_L1 → EDI_cond → [SPH] → DELETE_L1 →
                  ASM_L2 → FLU_L2 → EDIR
        """
        from .CLE2000 import wrap_cle2000_line

        # Retrieve edition step for EDI/SPH helpers
        edi_step = scheme.get_edition_between_levels_steps()[0]
        edi_cond = EDI_condensation(edi_step, lib_name="LIBEQ")
        sph_corr = SPH_correction(edi_step, lib_name="LIBEQ")

        flux_steps = scheme.get_flux_steps()
        l1_step = flux_steps[0]
        l2_step = flux_steps[1]

        proc = main_procedure(self.case_name)
        proc.enable_g2s = self.enable_g2s

        # --- Modules ---
        for m in ["LIB", "G2S", "SALT", "MCCGT",
                  "USS", "ASM", "FLU", "EDI",
                  "SPH", "GREP", "DELETE", "END"]:
            proc.add_module_call(m)

        # --- Sub-procedures ---
        for sp_name in (mix_proc_name, trk_proc_name, edir_proc_name):
            proc.add_sub_procedure(sub_procedure(sp_name))

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                sub_mixeq = sub_procedure(proc_name)
                proc.add_sub_procedure(sub_mixeq)        

        # --- Data structures ---
        proc.add_linked_list("LIBRARY")
        proc.add_linked_list("LIBRARY2")
        for trk_ll, _ in trk.get_track_names():
            proc.add_linked_list(trk_ll)
        proc.add_linked_list("SYS")
        proc.add_linked_list("FLUXL1")
        proc.add_linked_list("FLUXL2")
        proc.add_linked_list("EDITION")
        proc.add_linked_list(edi_cond.lib_name)
        proc.add_linked_list("COMPO")
        for _, trkfil in trk.get_track_names():
            proc.add_seq_binary(trkfil)

        # --- Variables ---
        aniso = ssh_step.anisotropy_level
        proc.add_variable("o_anis", "INTEGER", aniso)
        proc.add_variable("name_compo", "STRING",
                          f"_CPO_{self.case_name}")
        proc.add_variable("Draglib", "STRING", draglib_alias)
        proc.add_variable("keffL1", "REAL")
        proc.add_variable("keffL2", "REAL")

        # Temperature variables
        temps = self._extract_cle2000_temperatures()
        for tname, temp_val in temps.items():
            proc.add_variable(tname, "REAL", temp_val)

        # --- TDT SEQ_ASCII imports ---
        tdt_map = trk.get_tdt_file_mapping()
        for var, fname in tdt_map.items():
            proc.add_seq_ascii(var, fname)

        # --- G2S block ---
        if proc.enable_g2s:
            proc.add_body_block(proc.build_g2s_block(tdt_map))

        # ---- track name aliases ----
        track_names = trk.get_track_names()
        tdt_vars = list(tdt_map.keys())
        ssh_trk,    ssh_trkfil    = track_names[0]
        l1_trk,     l1_trkfil     = track_names[1]
        l2_trk,     l2_trkfil     = track_names[2]

        # --- TRK sub-procedure call ---
        lhs = " ".join(f"{ll} {bf}" for ll, bf in track_names)
        rhs_tdts = " ".join(tdt_vars)
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Tracking sub-procedure call")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"{lhs} := {trk_proc_name} {rhs_tdts} :: <<o_anis>> ;"
            )
        )
        proc.add_body_line("")

        # --- MIX sub-procedure call ---
        ssh_method = ssh_step.self_shielding_method or "PT"
        if hasattr(ssh_step, 'transport_correction') and ssh_step.transport_correction is not None:
            transport_correction = ssh_step.transport_correction
        else:
            transport_correction = "NONE"
        
        proc.add_body_line("*" * 50)
        proc.add_body_line("* LIBRARY creation")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"LIBRARY := {mix_proc_name} :: "
                f"<<Draglib>> '{ssh_method}' <<o_anis>> "
                f"'{transport_correction}' "
                f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
                f"<<TCOOL>> <<TMODE>> <<TCTRL>> ;"
            )
        )
        proc.add_body_line("")

        # --- USS call (SSH tracking) ---
        uss_keyword = "PIJ" if ssh_step.spatial_method == "CP" else "ARM"
        proc.add_body_line("*" * 50)
        proc.add_body_line("* USS: self-shielding")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"LIBRARY2 := USS: LIBRARY {ssh_trk} {ssh_trkfil} ::"
            )
        )
        proc.add_body_line(f"    EDIT 1 {uss_keyword}")
        proc.add_body_line("    PASS 3")
        proc.add_body_line(";")
        proc.add_body_line("")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == ssh_step.name and next_step_name == l1_step.name:
                    proc.add_body_line(
                        f"* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"LIBRARY2 := {proc_name} LIBRARY2 ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(
                            f"{mixeq_call}"
                        )
                    )
                    proc.add_body_line("")

        # --- ASM L1 + FLU L1 ---
        l1_asm_keyword = "PIJ" if l1_step.spatial_method == "CP" else "ARM"
        proc.add_body_line("*" * 50)
        proc.add_body_line("* ASM: + FLU: — first level")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"SYS := ASM: LIBRARY2 {l1_trk} {l1_trkfil} ::"
            )
        )
        proc.add_body_line(f"    EDIT 1 {l1_asm_keyword}")
        proc.add_body_line(";")
        proc.add_body_line("")
        proc.add_body_line(
            wrap_cle2000_line(
                f"FLUXL1 := FLU: SYS LIBRARY2 {l1_trk} {l1_trkfil} ::"
            )
        )
        proc.add_body_line("    EDIT 1 TYPE K")
        proc.add_body_line(";")
        proc.add_body_line(
            wrap_cle2000_line(
                "GREP: FLUXL1 :: GETVAL 'K-EFFECTIVE' 1 1 1 >>keffL1<< ;"
            )
        )
        proc.add_body_line(f'ECHO "{self.case_name} L1 keff=" keffL1 ;')
        proc.add_body_line("")

        # --- EDI condensation ---
        proc.add_body_line("*" * 50)
        proc.add_body_line("* EDI: energy condensation (L1 → coarse)")
        proc.add_body_line("*" * 50)
        proc.add_body_block(
            edi_cond.build_edi_call("FLUXL1", "LIBRARY2", l1_trk)
        )
        proc.add_body_block(edi_cond.build_lib_extract())
        proc.add_body_line("")

        # --- SPH correction (optional) ---
        if edi_step.sph_correction:
            proc.add_body_line("*" * 50)
            proc.add_body_line("* SPH: equivalence correction")
            proc.add_body_line("*" * 50)
            proc.add_body_block(
                sph_corr.build_sph_call(l1_trk, l1_trkfil)
            )
            proc.add_body_line("")

        # --- DELETE L1 objects ---
        delete_l1 = (
            f"FLUXL1 SYS EDITION := "
            f"DELETE: FLUXL1 SYS EDITION ;"
        )
        proc.add_body_line(wrap_cle2000_line(delete_l1))
        proc.add_body_line("")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == l1_step.name and next_step_name == l2_step.name:
                    proc.add_body_line(
                        f"* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"{edi_cond.lib_name} := {proc_name} {edi_cond.lib_name} ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(
                            f"{mixeq_call}"
                        )
                    )
                    proc.add_body_line("")

        # --- ASM L2 + FLU L2 ---
        l2_asm_keyword = "PIJ" if l2_step.spatial_method == "CP" else "ARM"
        proc.add_body_line("*" * 50)
        proc.add_body_line("* ASM: + FLU: — second level")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"SYS := ASM: {edi_cond.lib_name} {l2_trk} {l2_trkfil} ::"
            )
        )
        proc.add_body_line(f"    {l2_asm_keyword} EDIT 1")
        proc.add_body_line(";")
        proc.add_body_line("")
        proc.add_body_line(
            wrap_cle2000_line(
                f"FLUXL2 := FLU: {edi_cond.lib_name} SYS "
                f"{l2_trk} {l2_trkfil} ::"
            )
        )
        proc.add_body_line("    EDIT 1 TYPE K")
        proc.add_body_line(";")
        proc.add_body_line(
            wrap_cle2000_line(
                "GREP: FLUXL2 :: GETVAL 'K-EFFECTIVE' 1 1 1 >>keffL2<< ;"
            )
        )
        proc.add_body_line(f'ECHO "{self.case_name} L2 keff=" keffL2 ;')
        proc.add_body_line("")

        # --- EDIR call (uses coarse library and L2 tracking) ---
        proc.add_body_line("*" * 50)
        proc.add_body_line("* EDI/COMPO sub-procedure call")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"{edir_proc_name} FLUXL2 {edi_cond.lib_name} "
                f"{l2_trk} :: <<name_compo>> ;"
            )
        )

        return proc.write_to_x2m(self.output_path)

    def _build_main_x2m_with_branches_2level(
        self, scheme, trk, lib, edi_compo,
        mix_proc_name, trk_proc_name,
        edir_proc_name, mixeq_procedures,
        draglib_alias,
        ssh_step,
    ):
        """Build the main .x2m for a 2-level scheme with branch WHILE loops.

        Wraps the 2L inner body (USS → ASM_L1/FLU_L1 → EDI/SPH →
        DELETE_L1 → ASM_L2/FLU_L2 → EDIR → cleanup) in the same
        nested WHILE loop pattern used by ``_build_main_x2m_with_branches``.
        """
        from .CLE2000 import wrap_cle2000_line

        edi_step = scheme.get_edition_between_levels_steps()[0]
        edi_cond = EDI_condensation(edi_step, lib_name="LIBEQ")
        sph_corr = SPH_correction(edi_step, lib_name="LIBEQ")

        flux_steps = scheme.get_flux_steps()
        l1_step = flux_steps[0]
        l2_step = flux_steps[1]

        proc = main_procedure(self.case_name)
        proc.enable_g2s = self.enable_g2s

        # --- Modules ---
        for m in ["LIB", "G2S", "SALT", "MCCGT",
                  "USS", "ASM", "FLU", "EDI",
                  "SPH", "GREP", "UTL", "COMPO",
                  "DELETE", "END"]:
            proc.add_module_call(m)

        # --- Sub-procedures ---
        for sp_name in (mix_proc_name, trk_proc_name, edir_proc_name):
            proc.add_sub_procedure(sub_procedure(sp_name))
        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                sub_mixeq = sub_procedure(proc_name)
                proc.add_sub_procedure(sub_mixeq)

        # --- Data structures ---
        proc.add_linked_list("LIBRARY")
        proc.add_linked_list("LIBRARY2")
        for trk_ll, _ in trk.get_track_names():
            proc.add_linked_list(trk_ll)
        proc.add_linked_list("SYS")
        proc.add_linked_list("FLUXL1")
        proc.add_linked_list("FLUXL2")
        proc.add_linked_list("EDITION")
        proc.add_linked_list(edi_cond.lib_name)
        proc.add_linked_list("COMPO")
        proc.add_linked_list("PARAMS")
        for _, trkfil in trk.get_track_names():
            proc.add_seq_binary(trkfil)

        # --- Variables ---
        aniso = ssh_step.anisotropy_level
        proc.add_variable("o_anis", "INTEGER", aniso)
        proc.add_variable("name_compo", "STRING",
                          f"_CPO_{self.case_name}")
        proc.add_variable("Draglib", "STRING", draglib_alias)

        ordered_branches = scheme.get_ordered_branches()
        for branch in ordered_branches:
            proc.add_variable(branch.count_var, "INTEGER",
                              len(branch.values))
            proc.add_variable(branch.counter_var, "INTEGER", 0)
        for branch in ordered_branches:
            proc.add_variable(branch.para_name, "REAL")

        density_branch = scheme.get_branch("coolant_density")
        if density_branch:
            proc.add_variable("N_H", "REAL")
            proc.add_variable("N_O", "REAL")

        proc.add_variable("keffL1", "REAL")
        proc.add_variable("keffL2", "REAL")

        # --- TDT SEQ_ASCII imports ---
        tdt_map = trk.get_tdt_file_mapping()
        for var, fname in tdt_map.items():
            proc.add_seq_ascii(var, fname)
        proc.add_seq_ascii("_COMPO", f"_CPO_{self.case_name}")

        # --- G2S block ---
        if proc.enable_g2s:
            proc.add_body_block(proc.build_g2s_block(tdt_map))

        # ---- track name aliases ----
        track_names = trk.get_track_names()
        tdt_vars = list(tdt_map.keys())
        ssh_trk, ssh_trkfil = track_names[0]
        l1_trk,  l1_trkfil  = track_names[1]
        l2_trk,  l2_trkfil  = track_names[2]

        # --- TRK sub-procedure call (outside the loop) ---
        lhs = " ".join(f"{ll} {bf}" for ll, bf in track_names)
        rhs_tdts = " ".join(tdt_vars)
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Tracking sub-procedure call")
        proc.add_body_line("*" * 50)
        proc.add_body_line(
            wrap_cle2000_line(
                f"{lhs} := {trk_proc_name} {rhs_tdts} :: <<o_anis>> ;"
            )
        )
        proc.add_body_line("")

        # --- UTL: parameter arrays ---
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Parameter arrays for branch calculations")
        proc.add_body_line("*" * 50)
        utl_block = "PARAMS := UTL: ::\n"
        for branch in ordered_branches:
            vals_str = " ".join(f"{v}" for v in branch.values)
            utl_block += (
                f"    CREA\n"
                f"        {branch.array_name} "
                f"<<{branch.count_var}>> =\n"
                f"            {vals_str}\n"
            )
        utl_block += ";\n"
        proc.add_body_block(utl_block)

        if density_branch:
            iso_densities = compute_water_iso_densities_at_densities(
                density_branch.values
            )
            nh_vals = " ".join(
                f"{d['H1']:.8E}" for d in iso_densities
            )
            no_vals = " ".join(
                f"{d['O16']:.8E}" for d in iso_densities
            )
            n_dens = len(density_branch.values)
            proc.add_variable("n_dens", "INTEGER", n_dens)
            proc.add_body_block(
                f"PARAMS := UTL: PARAMS ::\n"
                f"    CREA\n"
                f"        N_H_arr <<n_dens>> =\n"
                f"            {nh_vals}\n"
                f"    CREA\n"
                f"        N_O_arr <<n_dens>> =\n"
                f"            {no_vals}\n"
                f";\n"
            )

        proc.add_body_line("")

        # --- COMPO initialisation ---
        proc.add_body_line("*" * 50)
        proc.add_body_line("* COMPO initialization")
        proc.add_body_line("*" * 50)
        proc.add_body_block(edi_compo.compo.build_compo_init())
        proc.add_body_line("")

        # ===================================================
        # Nested WHILE loops
        # ===================================================
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Main loop over branch statepoints")
        proc.add_body_line("*" * 50)

        for depth, branch in enumerate(ordered_branches):
            indent = "    " * depth
            proc.add_body_line(
                f"{indent}WHILE {branch.counter_var} "
                f"{branch.count_var} < DO"
            )
            proc.add_body_line(
                f"{indent}    EVALUATE {branch.counter_var} "
                f":= {branch.counter_var} 1 + ;"
            )
            proc.add_body_line(
                wrap_cle2000_line(
                    f"{indent}    GREP: PARAMS :: "
                    f"GETVAL '{branch.array_name}' "
                    f"<<{branch.counter_var}>> "
                    f">>{branch.para_name}<< ;"
                )
            )
            if branch.type == "coolant_density":
                proc.add_body_line(
                    wrap_cle2000_line(
                        f"{indent}    GREP: PARAMS :: "
                        f"GETVAL 'N_H_arr' "
                        f"<<{branch.counter_var}>> >>N_H<< ;"
                    )
                )
                proc.add_body_line(
                    wrap_cle2000_line(
                        f"{indent}    GREP: PARAMS :: "
                        f"GETVAL 'N_O_arr' "
                        f"<<{branch.counter_var}>> >>N_O<< ;"
                    )
                )

        inner_indent = "    " * len(ordered_branches)

        echo_parts = " ".join(
            f'"{b.para_name}=" {b.para_name}' for b in ordered_branches
        )
        proc.add_body_line(
            wrap_cle2000_line(f'{inner_indent}ECHO "State Point:" ;\n')
        )
        proc.add_body_line(
            wrap_cle2000_line(f'{inner_indent}ECHO {echo_parts} ;\n')
        )

        if scheme.get_branch("fuel_temperature"):
            proc.add_body_line(
                f"{inner_indent}EVALUATE TFUEL := TFuel ;"
            )
        if scheme.get_branch("coolant_temperature"):
            proc.add_body_line(
                f"{inner_indent}EVALUATE TCOOL := TCool ;"
            )
            proc.add_body_line(
                f"{inner_indent}EVALUATE TMODE := TCool ;"
            )
        proc.add_body_line("")

        # --- MIX ---
        ssh_method = ssh_step.self_shielding_method or "PT"
        if hasattr(ssh_step, 'transport_correction') and ssh_step.transport_correction is not None:
            transport_correction = ssh_step.transport_correction
        else:
            transport_correction = "NONE"

        if density_branch:
            mix_line = (
                f"LIBRARY := {mix_proc_name} :: "
                f"<<Draglib>> '{ssh_method}' <<o_anis>> "
                f"'{transport_correction}' "
                f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
                f"<<TCOOL>> <<TMODE>> <<TCTRL>> "
                f"<<N_H>> <<N_O>> ;"
            )
        else:
            mix_line = (
                f"LIBRARY := {mix_proc_name} :: "
                f"<<Draglib>> '{ssh_method}' <<o_anis>> "
                f"'{transport_correction}' "
                f"<<TFUEL>> <<TBOX>> <<TCLAD>> "
                f"<<TCOOL>> <<TMODE>> <<TCTRL>> ;"
            )
        proc.add_body_line(f"{inner_indent}* LIB: creation")
        proc.add_body_line(
            wrap_cle2000_line(f"{inner_indent}{mix_line}")
        )
        proc.add_body_line("")

        # --- USS ---
        uss_keyword = "PIJ" if ssh_step.spatial_method == "CP" else "ARM"
        proc.add_body_line(f"{inner_indent}* USS: self-shielding")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}LIBRARY2 := USS: LIBRARY "
                f"{ssh_trk} {ssh_trkfil} ::"
            )
        )
        proc.add_body_line(f"{inner_indent}    EDIT 1 {uss_keyword}")
        proc.add_body_line(f"{inner_indent}    PASS 3")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line("")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == ssh_step.name and next_step_name == l1_step.name:
                    proc.add_body_line(
                        f"{inner_indent}* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"{inner_indent}LIBRARY2 := {proc_name} LIBRARY2 ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(
                            f"{mixeq_call}"
                        )
                    )
                    proc.add_body_line("")

        # --- ASM L1 + FLU L1 ---
        l1_asm_keyword = "PIJ" if l1_step.spatial_method == "CP" else "ARM"
        proc.add_body_line(f"{inner_indent}* ASM: + FLU: — first level")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}SYS := ASM: LIBRARY2 "
                f"{l1_trk} {l1_trkfil} ::"
            )
        )
        proc.add_body_line(f"{inner_indent}    EDIT 1 {l1_asm_keyword}")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line("")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}FLUXL1 := FLU: SYS LIBRARY2 "
                f"{l1_trk} {l1_trkfil} ::"
            )
        )
        proc.add_body_line(f"{inner_indent}    EDIT 1 TYPE K")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}GREP: FLUXL1 :: GETVAL "
                f"'K-EFFECTIVE' 1 1 1 >>keffL1<< ;"
            )
        )
        proc.add_body_line(
            f'{inner_indent}ECHO "{self.case_name} L1 keff=" keffL1 ;'
        )
        proc.add_body_line("")

        # --- EDI condensation ---
        proc.add_body_line(
            f"{inner_indent}* EDI: energy condensation (L1 → coarse)"
        )

        def _indent_block(block, indent):
            return "".join(
                f"{indent}{line}\n" if line.strip() else "\n"
                for line in block.rstrip("\n").split("\n")
            )

        proc.add_body_block(
            _indent_block(
                edi_cond.build_edi_call("FLUXL1", "LIBRARY2", l1_trk),
                inner_indent,
            )
        )
        proc.add_body_block(
            _indent_block(edi_cond.build_lib_extract(), inner_indent)
        )
        proc.add_body_line("")

        # --- SPH correction (optional) ---
        if edi_step.sph_correction:
            proc.add_body_line(
                f"{inner_indent}* SPH: equivalence correction"
            )
            proc.add_body_block(
                _indent_block(
                    sph_corr.build_sph_call(l1_trk, l1_trkfil),
                    inner_indent,
                )
            )
            proc.add_body_line("")

        # --- DELETE L1 objects ---
        proc.add_body_line(f"{inner_indent}* Cleanup L1")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}FLUXL1 SYS EDITION := "
                f"DELETE: FLUXL1 SYS EDITION ;"
            )
        )
        proc.add_body_line("")

        if mixeq_procedures:
            for proc_name, _, current_step_name, next_step_name in mixeq_procedures:
                if current_step_name == l1_step.name and next_step_name == l2_step.name:
                    proc.add_body_line(
                        f"{inner_indent}* MIXEQ: {current_step_name} → {next_step_name}"
                    )
                    mixeq_call = (
                        f"{inner_indent}{edi_cond.lib_name} := "
                        f"{proc_name} {edi_cond.lib_name} ;"
                    )
                    proc.add_body_line(
                        wrap_cle2000_line(
                            f"{mixeq_call}"
                        )
                    )
                    proc.add_body_line("")

        # --- ASM L2 + FLU L2 ---
        l2_asm_keyword = "PIJ" if l2_step.spatial_method == "CP" else "ARM"
        proc.add_body_line(f"{inner_indent}* ASM: + FLU: — second level")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}SYS := ASM: {edi_cond.lib_name} "
                f"{l2_trk} {l2_trkfil} ::"
            )
        )
        proc.add_body_line(f"{inner_indent}    {l2_asm_keyword} EDIT 1")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line("")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}FLUXL2 := FLU: {edi_cond.lib_name} SYS "
                f"{l2_trk} {l2_trkfil} ::"
            )
        )
        proc.add_body_line(f"{inner_indent}    EDIT 1 TYPE K")
        proc.add_body_line(f"{inner_indent};")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}GREP: FLUXL2 :: GETVAL "
                f"'K-EFFECTIVE' 1 1 1 >>keffL2<< ;"
            )
        )
        proc.add_body_line(
            f'{inner_indent}ECHO "{self.case_name} L2 keff=" keffL2 ;'
        )
        proc.add_body_line("")

        # --- EDIR call (coarse library + L2 tracking) ---
        proc.add_body_line(f"{inner_indent}* EDI/COMPO call")
        para_args = " ".join(
            f"<<{b.para_name}>>" for b in ordered_branches
        )
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}COMPO := {edir_proc_name} FLUXL2 "
                f"{edi_cond.lib_name} "
                f"{l2_trk} COMPO :: <<name_compo>> "
                f"{para_args} ;"
            )
        )
        proc.add_body_line("")

        # --- Cleanup per statepoint ---
        proc.add_body_line(f"{inner_indent}* Cleanup")
        proc.add_body_line(
            wrap_cle2000_line(
                f"{inner_indent}LIBRARY LIBRARY2 "
                f"{edi_cond.lib_name} := "
                f"DELETE: LIBRARY LIBRARY2 "
                f"{edi_cond.lib_name} ;"
            )
        )
        proc.add_body_line(
            f"{inner_indent}FLUXL2 := DELETE: FLUXL2 ;"
        )
        proc.add_body_line(
            f"{inner_indent}SYS := DELETE: SYS ;"
        )

        # Close loops (innermost first)
        for depth in range(len(ordered_branches) - 1, -1, -1):
            indent = "    " * depth
            branch = ordered_branches[depth]
            proc.add_body_line(f"{indent}ENDWHILE ;")
            if depth > 0:
                proc.add_body_line(
                    f"{indent}EVALUATE {branch.counter_var} := 0 ;"
                )

        proc.add_body_line("")

        # --- COMPO export ---
        proc.add_body_line("*" * 50)
        proc.add_body_line("* Export COMPO to file")
        proc.add_body_line("*" * 50)
        proc.add_body_line("_COMPO := COMPO ;")

        # Temperature variables (from MaterialMixture objects)
        temps = self._extract_cle2000_temperatures()
        for tname, temp_val in temps.items():
            proc.add_variable(
                tname, "REAL", temp_val
            )

        return proc.write_to_x2m(self.output_path)