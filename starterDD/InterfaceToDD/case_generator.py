# Case generator for CLE2000 procedure generation
# A DRAGON5 case is defined by a main procedure and its
# sub-procedures, and the module calls in the main procedure.
#
# R.Guasch — 10/03/2026

from .CLE2000 import (
    main_procedure, sub_procedure,
    validate_varname, wrap_cle2000_line,
)
from .dragon_module_calls import LIB, EDI_COMPO, TRK
from ..DDModel.DragonCalculationScheme import (
    DragonCalculationScheme,
)
from ..DDModel.helpers import associate_material_to_rod_ID
from ..MaterialProperties.material_mixture import (
    parse_all_compositions_from_yaml,
)

import os


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
       - **MIX_<case>.c2m** — material library (LIB:).
       - **TRK_<case>.c2m** — tracking (SALT:/MCCGT:).
       - **EDIR_<case>.c2m** — editions/COMPO export.
    """

    def __init__(self, case_name, call_glow,
                 draglibs_names_to_alias, config_yamls,
                 enable_g2s=False, output_path=None):
        """
        Parameters
        ----------
        case_name : str
            Identifier for the DRAGON case.
        call_glow : bool
            If ``True``, call glow to generate TDT geometry
            files.  If ``False`` the TDT files are expected
            to exist already.
        draglibs_names_to_alias : dict
            ``{draglib_file_name: cle2000_alias}``.
        config_yamls : dict
            ``{"MATS": path, "GEOM": path,
              "CALC_SCHEME": path}``.
        enable_g2s : bool
            Include G2S: geometry visualization calls in x2m.
        output_path : str or None
            Directory for generated CLE2000 files.  Defaults
            to ``"./cle2000_procs/<case_name>"``.
        """
        self.case_name = case_name
        self.call_glow = call_glow
        self.draglibs_names_to_alias = draglibs_names_to_alias
        self.config_yamls = config_yamls
        self.enable_g2s = enable_g2s
        self.output_path = output_path or os.path.join(
            "cle2000_procs", case_name,
        )
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

    # -----------------------------------------------------------
    # Assembly model construction
    # -----------------------------------------------------------

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
        assembly.set_uniform_temperatures(
            fuel_temperature=900.0,
            gap_temperature=600.0,
            coolant_temperature=600.0,
            moderator_temperature=600.0,
            structural_temperature=600.0,
        )
        assembly.analyze_lattice_description(build_pins=True)
        assembly.set_material_compositions(compositions)
        return assembly, compositions

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
        assembly.number_fuel_material_mixtures_by_pin()

        # 4. Identify generating / daughter mixes
        assembly.identify_generating_and_daughter_mixes()

        # Determine the first draglib alias
        draglib_alias = list(
            self.draglibs_names_to_alias.values()
        )[0]

        # ----- MIX.c2m (LIB procedure) -----
        mix_proc_name = f"MIX_{self.case_name}"
        validate_varname(mix_proc_name[:12])
        lib = LIB(assembly)
        mix_path = lib.write_to_c2m(
            self.output_path, mix_proc_name,
        )

        # ----- TRK.c2m -----
        trk_proc_name = f"TRK_{self.case_name}"
        trk = TRK(scheme, self.case_name)
        trk_path = trk.write_to_c2m(
            self.output_path, trk_proc_name,
        )

        # ----- EDIR.c2m (EDI/COMPO) -----
        edir_proc_name = f"EDIR_{self.case_name}"
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
        x2m_path = self._build_main_x2m(
            scheme, trk, lib, edi_compo,
            mix_proc_name, trk_proc_name,
            edir_proc_name, draglib_alias,
            ssh_step,
        )

        return {
            "x2m": x2m_path,
            "mix": mix_path,
            "trk": trk_path,
            "edir": edir_path,
        }

    # -----------------------------------------------------------
    # Main x2m builder
    # -----------------------------------------------------------

    def _build_main_x2m(self, scheme, trk, lib, edi_compo,
                        mix_proc_name, trk_proc_name,
                        edir_proc_name, draglib_alias,
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

        # --- Data structures ---
        proc.add_linked_list("LIBRARY")
        proc.add_linked_list("LIBRARY2")
        for trk_ll, _ in trk.get_track_names():
            proc.add_linked_list(trk_ll)
        proc.add_linked_list("SYS")
        proc.add_linked_list("FLUX")

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
        temp_mapping = {
            "TFUEL": "fuel_temperature",
            "TBOX": "structural_temperature",
            "TCLAD": "gap_temperature",
            "TCOOL": "coolant_temperature",
            "TMODE": "moderator_temperature",
            "TCTRL": "structural_temperature",
        }
        for tname, attr in temp_mapping.items():
            temp_val = getattr(
                self.assembly, attr, 600.0
            )
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
            f"{rhs_tdts} :: <<o_anis>> ;"
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
        ssh_method = ssh_step.self_shielding_method or "RSE"
        mix_call_line = (
            f"LIBRARY := {mix_proc_name} :: "
            f"<<Draglib>> '{ssh_method}' <<o_anis>> "
            f"'NONE' "
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
        proc.add_body_line("    EDIT 1 ARM")
        proc.add_body_line(";")
        proc.add_body_line("")

        # Clean up SSH tracking
        proc.add_body_line(
            wrap_cle2000_line(
                f"{ssh_trk} {ssh_trkfil} := "
                f"DELETE: {ssh_trk} {ssh_trkfil} ;"
            )
        )
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
            proc.add_body_line("    EDIT 1 ARM")
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

    