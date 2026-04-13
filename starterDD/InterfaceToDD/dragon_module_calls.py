## collection of classes to handle DRAGON module calls
# Author : R. Guasch
# Date : 05/02/2026 (creation)
# Purpose : Define class structures to handle DRAGON module calls for starterDD package
# ----------------------------------------------------------------------------------------------- 

# classes representing DRAGON modules 
import os


class LIB:
    """
    Class for generating DRAGON5 ``LIB:`` module procedure files (``.c2m``).

    Given a ``CartesianAssemblyModel`` whose fuel material mixtures have been
    numbered and whose generating / daughter mixes have been identified, this
    class generates:

    * Full isotopic composition definitions for **generating** mixes.
    * ``COMB`` duplication lines for **daughter** mixes.
    * ``NOEV`` definitions for non-fuel (structural) materials.
    * The complete ``LIB:`` module call.
    * A complete ``.c2m`` procedure file.

    Fuel mixture indices are written as **plain integers** (no CLE-2000
    ``INTEGER`` variable declarations), because each ``MaterialMixture``
    object already carries a unique ``material_mixture_index``.  This avoids
    the 12-character variable-name limitation in CLE-2000.

    Usage
    -----
    ::

        assembly.number_fuel_material_mixtures_by_pin()
        assembly.identify_generating_and_daughter_mixes()

        lib = LIB(assembly)
        # Aliases are auto-populated from the YAML ``therm`` flag.
        # Manual overrides remain available as a fallback:
        #   lib.set_isotope_alias("MODERATOR", "H1", "H1_H2O")
        lib.write_to_c2m("./procs", "MIX_LIB")
    """

    SUPPORTED_SELF_SHIELDING_METHODS = ("PT", "RSE", "SUBG")

    # Heavy-metal + Gd isotopes that receive self-shielding treatment (INRS)
    DEFAULT_SELF_SHIELDED_FUEL_ISOTOPES = [
        "U234", "U235", "U236", "U238",
        "Pu239", "Pu240", "Pu241", "Pu242",
        "Gd152", "Gd154", "Gd155", "Gd156", "Gd157", "Gd158", "Gd160",
    ]
    
    DEFAULT_NON_FUEL_SELF_SHIELDED_ISOTOPES = [
        "Zr90", "Zr91", "Zr92", "Zr94", "Zr96",
        "Fe56",
        "Cr52",
        "Hf174", "Hf176", "Hf177", "Hf178", "Hf179", "Hf180",
    ]

    # CLE-2000 temperature variable used for each non-fuel material
    DEFAULT_TEMPERATURE_MAP = {
        "CHANNEL_BOX": "TBOX",
        "CLAD":        "TCLAD",
        "GAP":         "TFUEL",
        "MODERATOR":   "TMODE",
        "COOLANT":     "TCOOL",
        "ABS":         "TCTRL",
        "SHEATH":      "TCTRL",
    }

    def __init__(self, assembly_model, density_branch=False, ssh_calculation_step = None):
        """
        Initialize LIB procedure generator.

        Parameters
        ----------
        assembly_model : CartesianAssemblyModel
            Assembly model with:

            1. Fuel material mixtures numbered
               (``number_fuel_material_mixtures_by_material`` or ``_by_pin``).
            2. Generating / daughter mixes identified
               (``identify_generating_and_daughter_mixes``).
            3. Optionally, TDT indices enforced
               (``enforce_material_mixture_indices_from_tdt``).
        density_branch : bool
            When ``True``, COOLANT and MODERATOR mix H1/O16 densities
            are replaced by ``<<N_H>>`` / ``<<N_O>>`` CLE2000 variables
            and the procedure receives two extra REAL input parameters.
        """
        self.assembly = assembly_model
        self.density_branch = density_branch
        self.ssh_calculation_step = ssh_calculation_step

        # Fuel self-shielding configuration
        self.self_shielded_fuel_isotopes = list(self.DEFAULT_SELF_SHIELDED_FUEL_ISOTOPES)
        self.fuel_inrs = 1  # INRS value for fuel self-shielding group

        # Non-fuel self-shielding: {material_name: {isotope_name: inrs_value}}
        self.default_non_fuel_inrs = {isotope: 2 for isotope in self.DEFAULT_NON_FUEL_SELF_SHIELDED_ISOTOPES}
        self.non_fuel_inrs = {}

        # CLE-2000 temperature variable per non-fuel material
        self.non_fuel_temperature_map = dict(self.DEFAULT_TEMPERATURE_MAP)

        # Resonance correlation isotopes (get ``CORR`` keyword)
        self.correlation_isotopes = []  # e.g. ["U238", "Pu240", "Gd157"]

        # Library evaluation-name overrides: (material_name, isotope) → lib_name
        self.isotope_aliases = {}

        # Manually supplied non-fuel mixes (when TDT enforcement is skipped)
        self._extra_non_fuel_mixes = []

        # --- Auto-populate aliases from the Composition.therm flag --------
        self._auto_populate_therm_aliases()

        # --- Register control cross material temperature variables --------
        self._register_control_cross_temp_vars()

        if ssh_calculation_step is not None:
            # Override self-shielding configuration based on the provided calculation step
            if ssh_calculation_step.fuel_self_shielded_isotopes is not None:
                self.set_self_shielded_fuel_isotopes(ssh_calculation_step.fuel_self_shielded_isotopes)
            if ssh_calculation_step.clad_self_shielded_isotopes is not None:
                self.set_non_fuel_self_shielding(
                    "CLAD",
                    {iso: 2 for iso in ssh_calculation_step.clad_self_shielded_isotopes}
                )

    # ------------------------------------------------------------------
    #  Configuration helpers
    # ------------------------------------------------------------------

    def set_self_shielded_fuel_isotopes(self, isotopes):
        """Override the list of self-shielded fuel isotopes (default: HM + Gd)."""
        print(f"[LIB] Setting self-shielded fuel isotopes to: {isotopes}")
        self.self_shielded_fuel_isotopes = list(isotopes)

    def set_non_fuel_self_shielding(self, material_name, isotope_inrs_dict):
        """
        Set self-shielding INRS options for a non-fuel material.

        Parameters
        ----------
        material_name : str
            E.g. ``"CLAD"``, ``"CHANNEL_BOX"``.
        isotope_inrs_dict : dict
            ``{isotope_name: inrs_value}`` for isotopes that require
            self-shielding in this material, e.g.
            ``{"Zr90": 3, "Zr91": 3, "Fe56": 3}``.
        """
        self.non_fuel_inrs[material_name] = isotope_inrs_dict

    def set_non_fuel_temperature_variable(self, material_name, temp_variable):
        """
        Override the CLE-2000 temperature variable for a non-fuel material.

        Parameters
        ----------
        material_name : str
            E.g. ``"CHANNEL_BOX"``.
        temp_variable : str
            CLE-2000 variable name, e.g. ``"TBOX"``.
        """
        self.non_fuel_temperature_map[material_name] = temp_variable

    def set_correlation_isotopes(self, isotopes):
        """
        Set isotopes that receive the ``CORR`` keyword (resonance correlation
        model for PT / RSE self-shielding methods).

        Parameters
        ----------
        isotopes : list of str
            e.g. ``["U238", "Pu240", "Gd157"]``
        """
        self.correlation_isotopes = list(isotopes)

    def set_isotope_alias(self, material_name, isotope_name, library_name):
        """
        Set a library evaluation-name override for an isotope in a specific
        material (e.g. bound scattering kernels: ``H1 → H1_H2O`` in water).

        This is a **manual fallback** — aliases are auto-populated from the
        ``therm`` flag in the YAML composition file.  Calling this method
        will override the auto-detected alias for the given
        ``(material_name, isotope_name)`` pair.

        Parameters
        ----------
        material_name : str
            E.g. ``"MODERATOR"``, ``"COOLANT"``.
        isotope_name : str
            Isotope identifier in the ``Composition`` (e.g. ``"H1"``).
        library_name : str
            Library evaluation name in the DRAGLIB (e.g. ``"H1_H2O"``).
        """
        self.isotope_aliases[(material_name, isotope_name)] = library_name

    # ------------------------------------------------------------------
    #  Thermal-scattering auto-detection
    # ------------------------------------------------------------------

    def _auto_populate_therm_aliases(self):
        """Inspect the assembly's ``composition_lookup`` and register Dragon
        isotope aliases for every composition whose ``therm`` flag is set.

        Only entries **not** already present in ``self.isotope_aliases`` are
        added, so manual ``set_isotope_alias`` calls always take priority.
        """
        comp_lookup = getattr(self.assembly, 'composition_lookup', None)
        if comp_lookup is None:
            return
        for mat_name, composition in comp_lookup.items():
            if not getattr(composition, 'therm', False):
                continue
            for entry in getattr(composition, 'therm_data', []):
                key = (mat_name, entry["isotope"])
                if key not in self.isotope_aliases:
                    self.isotope_aliases[key] = entry["dragon_alias"]

    def _register_control_cross_temp_vars(self):
        """Ensure the assembly's control cross material names have entries
        in ``non_fuel_temperature_map``, even when the user chose custom
        material names that are not in ``DEFAULT_TEMPERATURE_MAP``.
        """
        ctrl = getattr(self.assembly, 'control_cross', None)
        if ctrl is None:
            return
        for mat_name in (ctrl.absorber_material, ctrl.sheath_material):
            if mat_name not in self.non_fuel_temperature_map:
                self.non_fuel_temperature_map[mat_name] = "TCTRL"

    def add_non_fuel_mix(self, material_name, mix_index, composition,
                         temperature_variable="TCOOL"):
        """
        Manually register a non-fuel material for inclusion in the ``LIB:``
        call (for when ``enforce_material_mixture_indices_from_tdt`` was not
        called and mix indices are set by hand).

        Parameters
        ----------
        material_name : str
        mix_index : int
        composition : Composition
        temperature_variable : str
        """
        self._extra_non_fuel_mixes.append({
            "name": material_name,
            "index": mix_index,
            "composition": composition,
            "temp_var": temperature_variable,
        })

    # ------------------------------------------------------------------
    #  Build blocks
    # ------------------------------------------------------------------

    def build_generating_mix_lines(self):
        """
        Build ``MIX <index> <<TFUEL>>`` lines for generating mixes with full
        isotopic composition and self-shielding (INRS) flags.

        Returns
        -------
        str
        """
        lines = ""
        for mix in self.assembly.generating_mixes:
            idx = mix.material_mixture_index
            lines += f"    MIX {idx} <<TFUEL>> \n"
            iso_comp = mix.composition.get_isotope_name_composition()
            for isotope, density in iso_comp.items():
                lib_name = self.isotope_aliases.get(
                    (mix.material_name, isotope), isotope
                )
                line = f"    {isotope} = {lib_name} {density:.6E}"
                if isotope in self.correlation_isotopes:
                    line += f" CORR {self.fuel_inrs}"
                elif isotope in self.self_shielded_fuel_isotopes:
                    line += f" {self.fuel_inrs}"
                lines += line + "\n"
        return lines

    def build_daughter_mix_lines(self):
        """
        Build ``MIX <index> COMB <generating_index> 1.0`` lines for daughter
        mixes.

        Returns
        -------
        str
        """
        lines = ""
        for mix in self.assembly.daughter_mixes:
            idx = mix.material_mixture_index
            gen_idx = mix.generating_mix.material_mixture_index
            lines += f"    MIX {idx} COMB {gen_idx} 1.0\n"
        return lines

    def build_non_fuel_mix_lines(self):
        """
        Build ``MIX <index> <<T...>> NOEV`` definitions for non-fuel
        materials.

        Returns
        -------
        str
        """
        entries = []  # (name, index, composition, temp_var)

        # From TDT enforcement
        if hasattr(self.assembly, 'non_fuel_material_mixture_indices'):
            for mat_name, idx in self.assembly.non_fuel_material_mixture_indices.items():
                composition = self.assembly.composition_lookup.get(mat_name)
                if composition is None:
                    print(f"[LIB] Warning: no composition found for "
                          f"non-fuel material '{mat_name}', skipping.")
                    continue
                temp_var = self.non_fuel_temperature_map.get(mat_name, "TCOOL")
                entries.append((mat_name, idx, composition, temp_var))

        # Manually added
        for extra in self._extra_non_fuel_mixes:
            entries.append(
                (extra["name"], extra["index"],
                 extra["composition"], extra["temp_var"])
            )

        # Isotopes whose density is replaced by CLE2000 variables
        # when coolant-density branching is active.
        _WATER_ISO_MAPPING = {
            "H1": "N_H", "H2": "N_H",
            "O16": "N_O", "O17": "N_O",
        }
        # Materials affected by coolant-density branching
        _WATER_MATERIALS = {"COOLANT", "MODERATOR"}

        lines = ""
        for mat_name, idx, composition, temp_var in entries:
            iso_comp = composition.get_isotope_name_composition()
            inrs_map = self.non_fuel_inrs.get(mat_name, self.default_non_fuel_inrs)
            depletable = getattr(composition, 'depletable', False)
            if depletable:
                lines += f"    MIX {idx} <<{temp_var}>>\n"
            else:
                lines += f"    MIX {idx} <<{temp_var}>> NOEV\n"

            use_density_vars = (
                self.density_branch
                and mat_name in _WATER_MATERIALS
            )

            for isotope, density in iso_comp.items():
                lib_name = self.isotope_aliases.get(
                    (mat_name, isotope), isotope
                )
                if use_density_vars and isotope in _WATER_ISO_MAPPING:
                    var_name = _WATER_ISO_MAPPING[isotope]
                    line = f"    {isotope} = {lib_name} <<{var_name}>>"
                else:
                    line = f"    {isotope} = {lib_name} {density:.6E}"
                if isotope in inrs_map:
                    line += f" {inrs_map[isotope]}"
                lines += line + "\n"
            lines += "\n"

        return lines

    def build_control_cross_generating_mix_lines(self):
        """
        Build full isotopic ``MIX`` definitions for control cross absorber
        generating mixes (per-tube numbering only).

        Returns an empty string if there are no per-tube absorber mixtures.

        Returns
        -------
        str
        """
        if (not hasattr(self.assembly, 'control_cross_generating_mixes')
                or not self.assembly.control_cross_generating_mixes):
            return ""

        temp_var = self.non_fuel_temperature_map.get(
            self.assembly.control_cross.absorber_material, "TCTRL"
        )
        lines = ""
        for mix in self.assembly.control_cross_generating_mixes:
            idx = mix.material_mixture_index
            lines += f"    MIX {idx} <<{temp_var}>>\n"
            iso_comp = mix.composition.get_isotope_name_composition()
            if self.non_fuel_inrs:
                inrs_map = self.non_fuel_inrs.get(mix.material_name, {})
            else:
                inrs_map = self.default_non_fuel_inrs
            for isotope, density in iso_comp.items():
                lib_name = self.isotope_aliases.get(
                    (mix.material_name, isotope), isotope
                )
                line = f"    {isotope} = {lib_name} {density:.6E}"
                if isotope in inrs_map:
                    line += f" {inrs_map[isotope]}"
                lines += line + "\n"
        return lines

    def build_control_cross_daughter_mix_lines(self):
        """
        Build ``MIX <index> COMB <generating_index> 1.0`` lines for
        control cross absorber daughter mixes (per-tube numbering only).

        Returns an empty string if there are no per-tube absorber mixtures.

        Returns
        -------
        str
        """
        if (not hasattr(self.assembly, 'control_cross_daughter_mixes')
                or not self.assembly.control_cross_daughter_mixes):
            return ""

        lines = ""
        for mix in self.assembly.control_cross_daughter_mixes:
            idx = mix.material_mixture_index
            gen_idx = mix.generating_mix.material_mixture_index
            lines += f"    MIX {idx} COMB {gen_idx} 1.0\n"
        return lines

    def build_mix_index_comment_block(self):
        """
        Build a comment block listing all fuel and non-fuel mix index
        assignments.  This replaces the CLE-2000 ``INTEGER`` variable
        declarations and serves as human-readable documentation in the
        generated ``.c2m`` file.

        Returns
        -------
        str
        """
        lines = "* --------------------------------\n"
        lines += "*    MIX INDEX ASSIGNMENTS\n"
        lines += "* --------------------------------\n"
        lines += "* -- Fuel mixes --\n"
        for mix in self.assembly.fuel_material_mixtures:
            if mix.is_generating:
                tag = " (generating)"
            else:
                tag = f" (daughter of mix {mix.generating_mix.material_mixture_index})"
            lines += (f"*   {mix.material_mixture_index:4d} : "
                      f"{mix.unique_material_mixture_name}{tag}\n")

        if (hasattr(self.assembly, 'non_fuel_material_mixture_indices')
                and self.assembly.non_fuel_material_mixture_indices):
            lines += "* -- Non-fuel mixes --\n"
            for name, idx in self.assembly.non_fuel_material_mixture_indices.items():
                lines += f"*   {idx:4d} : {name}\n"

        # Control cross absorber tube mixes (per-tube numbering)
        if (hasattr(self.assembly, 'control_cross_absorber_mixtures')
                and self.assembly.control_cross_absorber_mixtures):
            lines += "* -- Control cross absorber tubes --\n"
            for mix in self.assembly.control_cross_absorber_mixtures:
                if mix.is_generating:
                    tag = " (generating)"
                else:
                    tag = f" (daughter of mix {mix.generating_mix.material_mixture_index})"
                lines += (f"*   {mix.material_mixture_index:4d} : "
                          f"{mix.unique_material_mixture_name}{tag}\n")

        for extra in self._extra_non_fuel_mixes:
            lines += f"*   {extra['index']:4d} : {extra['name']}\n"

        return lines

    def build_lib_module_call(self):
        """
        Assemble the complete ``LIBRARY := LIB: :: ... ;`` block.

        Returns
        -------
        str
        """
        max_mix = self._get_max_mix_index()
        generating = self.build_generating_mix_lines()
        daughters = self.build_daughter_mix_lines()
        ctrl_generating = self.build_control_cross_generating_mix_lines()
        ctrl_daughters = self.build_control_cross_daughter_mix_lines()
        non_fuel = self.build_non_fuel_mix_lines()

        if hasattr(self, 'ssh_calculation_step') and self.ssh_calculation_step is not None:
            ssh_method = self.ssh_calculation_step.self_shielding_method
            if ssh_method not in self.SUPPORTED_SELF_SHIELDING_METHODS:
                raise ValueError(
                    f"Unsupported self-shielding method '{ssh_method}' in "
                    f"SSH calculation step '{self.ssh_calculation_step.name}'. "
                    f"Supported methods: {self.SUPPORTED_SELF_SHIELDING_METHODS}"
                )
            if ssh_method == "PT":
                ssh_method_str = "PT CALENDF 4"
            elif ssh_method == "RSE":
                ssh_method_str = "RSE"
            elif ssh_method == "SUBG":
                ssh_method_str = "SUBG"
        else:
            ssh_method_str = "<<ssh_method>>"  # receive as input parameter when not fixed

        call = (
            "LIBRARY := LIB: ::\n"
            "EDIT 0\n"
            f"NMIX {max_mix}  ! MAXIMUM OF MATERIAL MIXTURES\n"
            f"{ssh_method_str}\n"
            "ANIS <<anis_level>>\n"
            "CTRA <<tran_correc>>\n"
            "ADED 4 NELAS N2N N3N N4N\n"
            "DEPL LIB: DRAGON FIL: <<Library>>\n"
            "MIXS LIB: DRAGON FIL: <<Library>>\n"
            f"{generating}"
            f"{daughters}"
            f"{ctrl_generating}"
            f"{ctrl_daughters}"
            f"{non_fuel}"
            ";\n"
        )
        return call

    # ------------------------------------------------------------------
    #  Write to file
    # ------------------------------------------------------------------

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the complete ``LIB:`` module procedure to a ``.c2m`` file.

        The generated file is a standalone CLE-2000 procedure that receives
        the DRAGLIB path, self-shielding method, anisotropy level, transport
        correction and material temperatures as input parameters.

        Parameters
        ----------
        path_to_procs : str
            Directory path for the output file.
        proc_name : str
            Procedure name (without ``.c2m`` extension).

        Returns
        -------
        str
            Relative path of the written file.
        """
        # Extra REAL parameters for coolant density branch
        if self.density_branch:
            density_param_block = (
                "REAL N_H N_O ;\n"
                ":: >>N_H<< >>N_O<< ;\n"
            )
        else:
            density_param_block = ""

        header = (
            f"*PROCEDURE {proc_name}.c2m\n"
            "* --------------------------------\n"
            "* Procedure generated by starterDD\n"
            "* Author: R. Guasch\n"
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
            "PARAMETER LIBRARY ::\n"
            "::: LINKED_LIST LIBRARY ; ;\n"
            "\n"
            "! Library name, ssh method\n"
            "STRING Library ssh_method ;\n"
            ":: >>Library<<  >>ssh_method<< ;\n"
            "INTEGER anis_level ;\n"
            ":: >>anis_level<< ; ! Anisotropy level\n"
            "STRING tran_correc ;\n"
            ":: >>tran_correc<< ; ! Transport correction option\n"
            "REAL TFUEL TBOX TCLAD TCOOL TMODE TCTRL ;\n"
            ":: >>TFUEL<< >>TBOX<< >>TCLAD<< >>TCOOL<< >>TMODE<< >>TCTRL<< ;\n"
            f"{density_param_block}"
            "\n"
            "* -------------------------------\n"
            "*    STRUCTURES AND MODULES\n"
            "* -------------------------------\n"
            "MODULE  LIB: UTL: DELETE: END: ABORT: ;\n"
            "\n"
        )

        mix_comment = self.build_mix_index_comment_block()
        lib_call = self.build_lib_module_call()

        footer = (
            "* -----------------------------------------\n"
            "*         END OF LIBRARY DEFINITION\n"
            "* -----------------------------------------\n"
            "END: ;\n"
            "QUIT ."
        )

        content = (
            f"{header}"
            f"{mix_comment}"
            "\n"
            "* --------------------------------\n"
            "*    MIX COMPOSITION DEFINITION\n"
            "*    CALL TO LIB: MODULE\n"
            "* --------------------------------\n"
            f"{lib_call}"
            f"{footer}"
        )

        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(path_to_procs, f"{proc_name}.c2m")
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[LIB] Wrote procedure to {filepath}")
        return filepath

    # ------------------------------------------------------------------
    #  Internal helpers
    # ------------------------------------------------------------------

    def _get_max_mix_index(self):
        """Return the highest material mixture index across all mixes."""
        all_indices = list(self.assembly.fuel_material_mixture_indices)
        if hasattr(self.assembly, 'non_fuel_material_mixture_indices'):
            all_indices.extend(
                self.assembly.non_fuel_material_mixture_indices.values()
            )
        if (hasattr(self.assembly, 'control_cross_absorber_mixtures')
                and self.assembly.control_cross_absorber_mixtures):
            all_indices.extend(
                m.material_mixture_index
                for m in self.assembly.control_cross_absorber_mixtures
            )
        for extra in self._extra_non_fuel_mixes:
            all_indices.append(extra["index"])
        return max(all_indices) if all_indices else 0


class EDI:
    """
    Class for generating DRAGON5 ``EDI:`` module call blocks.

    Given a ``CartesianAssemblyModel`` whose fuel material mixtures have been
    numbered (and optionally enforced from a TDT file), this class builds a
    single ``EDIRATES := EDI: FLUX LIBRARY2 TRACK :: ... ;`` block with the
    user's choice of:

    * Isotope selection (``MICR`` keyword).
    * Spatial homogenization (``MERG MIX`` keyword).
    * Energy condensation (``COND`` keyword, optional).
    * Output directory name (``SAVE ON`` keyword).

    Spatial homogenization modes
    ----------------------------
    ``"FUEL"``
        All fuel mixes → output region 1, non-fuel mixes → 0 (excluded).
    ``"ALL"``
        Every mix (fuel + non-fuel) → output region 1.
    ``"by_pin"``
        Each unique pin position gets its own output region (zones within
        the same pin are merged). Requires ``number_fuel_material_mixtures_by_pin()``
        numbering.
    ``"by_mix"``
        Each individual fuel ``MaterialMixture`` gets its own output region
        (full zone-level detail). Non-fuel → 0.
    ``"by_material"``
        All mixes sharing the same base composition (``material_name``)
        are merged into one output region. Non-fuel → 0.
    ``"custom"``
        User supplies ``{tdt_mix_index: output_region_index}``.

    Usage
    -----
    ::

        edi = EDI("EDIHOM_COND", assembly)
        edi.set_isotopes(["U235", "U238", "Gd155"])
        edi.set_spatial_homogenization("FUEL")
        edi.set_energy_condensation([])           # COND (1 group)
        print(edi.build_edi_call())
    """

    VALID_SPATIAL_MODES = ("FUEL", "ALL", "by_pin", "by_mix", "by_material", "custom")

    def __init__(self, name, assembly_model):
        """
        Parameters
        ----------
        name : str
            Directory name used in ``SAVE ON <name>`` (and later in COMPO).
        assembly_model : CartesianAssemblyModel
            Assembly with numbered material mixtures (and optionally TDT-
            enforced indices).
        """
        self.name = name
        self.assembly = assembly_model
        self.isotopes = []
        self.spatial_mode = "FUEL"
        self.custom_merge_map = None
        self.energy_bounds = None  # None → no COND; [] → COND; [0.625] → COND 0.625
        self.edit_level = 1

    # ------------------------------------------------------------------
    #  Configuration
    # ------------------------------------------------------------------

    def set_isotopes(self, isotopes):
        """
        Set the isotope list for the ``MICR`` keyword.

        Parameters
        ----------
        isotopes : list of str
            e.g. ``["U235", "U238", "U234", "Gd155", "Gd157"]``
        """
        self.isotopes = list(isotopes)

    def set_spatial_homogenization(self, mode, custom_map=None):
        """
        Set the spatial homogenization strategy for ``MERG MIX``.

        Parameters
        ----------
        mode : str
            One of ``"FUEL"``, ``"ALL"``, ``"by_pin"``, ``"by_mix"``,
            ``"by_material"``, ``"custom"``.
        custom_map : dict, optional
            Required when *mode* is ``"custom"``.
            ``{tdt_mix_index: output_region_index}``.
        """
        if mode not in self.VALID_SPATIAL_MODES:
            raise ValueError(
                f"Invalid spatial mode '{mode}'. "
                f"Valid modes: {self.VALID_SPATIAL_MODES}"
            )
        if mode == "custom" and custom_map is None:
            raise ValueError("custom_map is required when mode='custom'.")
        if mode == "by_pin":
            # Validate that by_pin numbering was used
            if self.assembly.fuel_material_mixtures:
                sample_name = self.assembly.fuel_material_mixtures[0].unique_material_mixture_name
                if "_pin" not in sample_name:
                    raise ValueError(
                        "Spatial mode 'by_pin' requires "
                        "number_fuel_material_mixtures_by_pin() numbering "
                        "(mixture names must contain '_pin<N>')."
                    )
        self.spatial_mode = mode
        self.custom_merge_map = custom_map

    def set_energy_condensation(self, bounds):
        """
        Set energy condensation for the ``COND`` keyword.

        Parameters
        ----------
        bounds : None, list of float
            ``None``    → no ``COND`` keyword (keep original energy mesh).
            ``[]``      → ``COND`` (collapse to 1 group).
            ``[0.625]`` → ``COND 0.625`` (2-group split at 0.625 eV).
        """
        self.energy_bounds = bounds

    def set_edit_level(self, level):
        """Set the ``EDIT`` level (default 1)."""
        self.edit_level = level

    # ------------------------------------------------------------------
    #  Builders
    # ------------------------------------------------------------------

    def _get_max_mix_index(self):
        """Return the highest material mixture index across all mixes."""
        all_indices = list(self.assembly.fuel_material_mixture_indices)
        if hasattr(self.assembly, 'non_fuel_material_mixture_indices'):
            all_indices.extend(
                self.assembly.non_fuel_material_mixture_indices.values()
            )
        return max(all_indices) if all_indices else 0

    def build_merg_mix_vector(self):
        """
        Build the integer vector for the ``MERG MIX`` keyword.

        The vector has length ``max_mix_index``.  Position ``i`` (0-indexed)
        corresponds to material mixture number ``i + 1``.  The value at each
        position determines which output region (in the EDI / COMPO result)
        the mixture is assigned to.

        Returns
        -------
        list of int
        """
        max_mix = self._get_max_mix_index()
        vector = [0] * max_mix  # default: excluded

        if self.spatial_mode == "FUEL":
            for mix in self.assembly.fuel_material_mixtures:
                vector[mix.material_mixture_index - 1] = 1

        elif self.spatial_mode == "ALL":
            for i in range(max_mix):
                vector[i] = 1

        elif self.spatial_mode == "by_pin":
            for mix in self.assembly.fuel_material_mixtures:
                # Parse pin_idx from "UOX28_zone_1_pin_3"
                pin_idx = int(
                    mix.unique_material_mixture_name.split("_")[-1]
                )
                vector[mix.material_mixture_index - 1] = pin_idx

        elif self.spatial_mode == "by_mix":
            for rank, mix in enumerate(self.assembly.fuel_material_mixtures, start=1):
                vector[mix.material_mixture_index - 1] = rank

        elif self.spatial_mode == "by_material":
            families = {}
            family_counter = 0
            for mix in self.assembly.fuel_material_mixtures:
                if mix.material_name not in families:
                    family_counter += 1
                    families[mix.material_name] = family_counter
                vector[mix.material_mixture_index - 1] = families[mix.material_name]

        elif self.spatial_mode == "custom":
            for tdt_idx, output_idx in self.custom_merge_map.items():
                if 1 <= tdt_idx <= max_mix:
                    vector[tdt_idx - 1] = output_idx

        return vector

    def _format_merg_mix_vector(self, vector):
        """Format the MERG MIX integer vector as a multi-line string (10 values per line)."""
        lines = ""
        for i, val in enumerate(vector):
            lines += f"{val} "
            if (i + 1) % 10 == 0:
                lines += "\n"
        return lines.strip()

    def build_cond_option(self):
        """Build the ``COND`` keyword string."""
        if self.energy_bounds is None:
            return ""
        elif self.energy_bounds == []:
            return "COND"
        else:
            bounds_str = " ".join(str(b) for b in self.energy_bounds)
            return f"COND {bounds_str}"

    def build_edi_call(self):
        """
        Assemble the complete ``EDIRATES := EDI: ... ;`` block.

        Returns
        -------
        str
        """
        from .CLE2000 import wrap_cle2000_line

        iso_str = " ".join(self.isotopes)
        vector = self.build_merg_mix_vector()
        merg_mix = self._format_merg_mix_vector(vector)
        cond = self.build_cond_option()

        micr_line = f"  MICR {len(self.isotopes)} {iso_str}"

        call = (
            "EDIRATES := EDI: FLUX LIBRARY2 TRACK ::\n"
            f"   EDIT {self.edit_level}\n"
            f"{wrap_cle2000_line(micr_line)}\n"
            "   MERG MIX\n"
            f"  {merg_mix}\n"
        )
        if cond:
            call += f"  {cond}\n"
        call += f"  SAVE ON {self.name}\n"
        call += ";\n"
        return call


class COMPO:
    """
    Class for generating DRAGON5 ``COMPO:`` module initialization and store
    call blocks.

    A ``COMPO`` object holds a list of named directories (volumes) that act
    as keys in the output COMPO heterogeneous list.  Each directory stores
    reaction rates and cross-section data computed by an ``EDI:`` call.

    When *branches* are provided the initialization block includes
    ``MAXCAL``, ``PARA`` keywords and the store block passes the current
    parameter values (AT10 multi-physics pattern).

    Usage
    -----
    ::

        compo = COMPO()
        compo.add_directory("EDIHOM_COND",
                            "Condensed, Homogenized over all fuel cells",
                            ["U235", "U238"])
        print(compo.build_compo_init())
        print(compo.build_compo_store("EDIHOM_COND"))
    """

    def __init__(self, branches=None):
        """
        Parameters
        ----------
        branches : list[CalculationBranch] or None
            When provided, COMPO initialization will include ``PARA``
            keywords and ``MAXCAL``.
        """
        self.directories = []  # list of (name, comment, isotopes)
        self.branches = branches or []

    def add_directory(self, name, comment, isotopes):
        """
        Register a named directory in the COMPO object.

        Parameters
        ----------
        name : str
            Directory key, e.g. ``"EDIHOM_COND"``.
        comment : str
            Human-readable description stored via ``COMM ... ENDC``.
        isotopes : list of str
            Isotopes tracked in this directory (``ISOT`` keyword).
        """
        self.directories.append((name, comment, list(isotopes)))

    def build_compo_init(self):
        """
        Build the ``COMPO := COMPO: :: STEP UP ... INIT ;`` initialization block.

        When branches are configured, includes ``MAXCAL`` and ``PARA``
        keywords for the parameter tree.

        Returns
        -------
        str
        """
        lines = "COMPO := COMPO: ::\n"
        lines += "   EDIT 0\n"
        for name, comment, isotopes in self.directories:
            iso_str = " ".join(isotopes)
            lines += f"   STEP UP '{name}'\n"
            if self.branches:
                # Compute MAXCAL from product of branch sizes
                total = 1
                for b in self.branches:
                    total *= len(b.values)
                lines += f"       MAXCAL {total}\n"
            lines += f"       COMM '{comment}' ENDC\n"
            if self.branches:
                for b in self.branches:
                    lines += f"       PARA '{b.para_name}' {b.para_keyword}\n"
            isot_line = f"      ISOT {len(isotopes)} {iso_str}"
            from .CLE2000 import wrap_cle2000_line as _wrap_isot
            lines += f"{_wrap_isot(isot_line)}\n"
            lines += "   INIT\n"
        lines += ";\n"
        return lines

    def build_compo_store(self, directory_name):
        """
        Build a ``COMPO := COMPO: COMPO EDIRATES LIBRARY2 :: STEP UP <name> ;``
        store block.

        When branches are configured, passes current parameter values
        using CLE2000 variable references ``<<TFuel>>``, etc.

        Parameters
        ----------
        directory_name : str
            Must match a name previously registered via ``add_directory``.

        Returns
        -------
        str
        """
        if self.branches:
            call = (
                "COMPO := COMPO: COMPO EDIRATES LIBRARY2 ::\n"
                f"   EDIT 1\n"
                f"    STEP UP '{directory_name}'\n"
            )
            for b in self.branches:
                if b.para_name == "Burnup":
                    call += f"    SET <<TIME>> DAY \n"
                else:
                    call += f"    '{b.para_name}' <<{b.para_name}>>\n"
            call += ";\n"
        else:
            call = (
                "COMPO := COMPO: COMPO EDIRATES LIBRARY2 ::\n"
                f"   EDIT 1\n"
                f"    STEP UP {directory_name}\n"
                ";\n"
            )
        return call


class EDI_COMPO:
    """
    Orchestrator that combines ``EDI`` and ``COMPO`` calls into a single
    CLE-2000 ``.c2m`` procedure file.

    Each *edition* added via ``add_edition()`` produces:

    1. A ``STEP UP`` directory in the ``COMPO:`` initialization.
    2. An ``EDIRATES := EDI: ...`` call.
    3. A ``COMPO := COMPO: COMPO EDIRATES LIBRARY2 :: STEP UP <name> ;`` store.
    4. An ``EDIRATES := DELETE: EDIRATES ;`` cleanup.

    When branches and outputs are configured, editions are built from
    ``CalculationOutput`` specifications and the COMPO initialization
    includes ``PARA`` / ``MAXCAL`` keywords.  Selective state-point
    outputs generate CLE2000 ``IF`` guards in the procedure body.

    The generated ``.c2m`` file structure matches the layout of a standard
    DRAGON5 EDI/COMPO procedure (see ``EDICPO_R.c2m`` for reference).

    Usage
    -----
    ::

        assembly.number_fuel_material_mixtures_by_pin()
        assembly.enforce_material_mixture_indices_from_tdt(tdt_indices)
        assembly.identify_generating_and_daughter_mixes()

        edi_compo = EDI_COMPO(assembly)

        edi_compo.add_edition(
            name="EDIHOM_COND",
            comment="Condensed, homogenized over all fuel cells",
            isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
            spatial_mode="FUEL",
            energy_bounds=[],
        )

        edi_compo.add_edition(
            name="H_EDI_REGI_2g",
            comment="Condensed to 2 groups, per unique region",
            isotopes=["U235", "U238", "U234", "Gd155", "Gd157"],
            spatial_mode="by_pin",
            energy_bounds=[0.625],
        )

        edi_compo.write_to_c2m("./procs", "EDICPO_GE14")
    """

    def __init__(self, assembly_model, branches=None, outputs=None):
        """
        Parameters
        ----------
        assembly_model : CartesianAssemblyModel
            Assembly model with numbered material mixtures and (optionally)
            TDT-enforced indices.
        branches : list[CalculationBranch] or None
            When provided, COMPO init includes PARA/MAXCAL and
            stores pass parameter values.
        outputs : list[CalculationOutput] or None
            When provided, editions are built from these specs
            via ``add_editions_from_outputs()``.
        """
        self.assembly = assembly_model
        self.branches = branches or []
        self.compo = COMPO(branches=self.branches)
        self.editions = []  # list of EDI objects
        self._output_specs = []  # parallel list of CalculationOutput (or None)

        if outputs:
            self.add_editions_from_outputs(outputs)

    def add_editions_from_outputs(self, outputs):
        """Build EDI objects from ``CalculationOutput`` specifications.

        Parameters
        ----------
        outputs : list[CalculationOutput]
        """
        for out in outputs:
            self.add_edition(
                name=out.name,
                comment=out.comment,
                isotopes=out.isotopes,
                spatial_mode=out.spatial_integration_mode,
                energy_bounds=out.energy_bounds,
                output_spec=out,
            )

    def add_edition(self, name, comment, isotopes, spatial_mode,
                    energy_bounds=None, custom_map=None,
                    output_spec=None):
        """
        Add an edition (EDI + COMPO store pair).

        Parameters
        ----------
        name : str
            Directory name for ``SAVE ON`` and ``STEP UP``.
        comment : str
            Human-readable description for the COMPO directory.
        isotopes : list of str
            Isotopes tracked (``MICR`` / ``ISOT``).
        spatial_mode : str
            ``"FUEL"``, ``"ALL"``, ``"by_pin"``, ``"by_mix"``,
            ``"by_material"``, ``"custom"``.
        energy_bounds : None or list of float
            ``None`` → no COND; ``[]`` → COND; ``[0.625]`` → COND 0.625.
        custom_map : dict, optional
            Required when ``spatial_mode="custom"``.
        output_spec : CalculationOutput or None
            When provided, used for selective state-point filtering.
        """
        edi = EDI(name, self.assembly)
        edi.set_isotopes(isotopes)
        edi.set_spatial_homogenization(spatial_mode, custom_map=custom_map)
        edi.set_energy_condensation(energy_bounds)
        self.editions.append(edi)
        self._output_specs.append(output_spec)
        self.compo.add_directory(name, comment, isotopes)

    def build_procedure_body(self):
        """
        Build the body of the CLE-2000 procedure (everything between header
        and footer).

        When branches are present and an output has selective
        ``state_points``, the EDI+COMPO block is wrapped in a CLE2000
        ``IF`` guard that checks whether the current parameter values
        match the requested subset.

        Returns
        -------
        str
        """
        body = ""

        # COMPO initialization — only when there are NO branches.
        # When branches are present the COMPO init is emitted in
        # the main x2m (before the loops) so we skip it here.
        if not self.branches:
            body += "* --------------------------------\n"
            body += "*   COMPO INITIALIZATION\n"
            body += "* --------------------------------\n"
            body += self.compo.build_compo_init()

        # EDI + COMPO pairs
        for idx, edi in enumerate(self.editions):
            output_spec = (
                self._output_specs[idx]
                if idx < len(self._output_specs)
                else None
            )
            needs_guard = (
                self.branches
                and output_spec is not None
                and output_spec.state_points != "ALL"
            )

            if needs_guard:
                guard = self._build_statepoint_guard(output_spec)
                body += guard
                indent = "    "
            else:
                indent = ""

            body += f"{indent}* --------------------------------\n"
            body += f"{indent}*    EDI: CALL FOR {edi.name}\n"
            body += f"{indent}* --------------------------------\n"
            # Indent the EDI call block
            from .CLE2000 import wrap_cle2000_line as _wrap_body
            for line in edi.build_edi_call().splitlines():
                body += _wrap_body(f"{indent}{line}") + "\n"
            body += f"{indent}* --------------------------------\n"
            body += f"{indent}*    COMPO: CALL FOR {edi.name}\n"
            body += f"{indent}* --------------------------------\n"
            for line in self.compo.build_compo_store(edi.name).splitlines():
                body += _wrap_body(f"{indent}{line}") + "\n"
            body += f"{indent}EDIRATES := DELETE: EDIRATES ;\n"

            if needs_guard:
                body += "ENDIF ;\n"

        return body

    def _build_statepoint_guard(self, output_spec):
        """Build a CLE2000 IF condition for selective state-point outputs.

        Generates conditions like::

            IF DCool 0.73669000 = TCool 600.0 = * TFuel 900.0 = * * THEN

        using CLE2000 reverse-Polish notation.

        Parameters
        ----------
        output_spec : CalculationOutput

        Returns
        -------
        str
            CLE2000 ``IF ... THEN`` line.
        """
        from ..DDModel.DragonCalculationScheme import (
            BRANCH_TYPE_TO_PARA_NAME,
        )

        conditions = []
        sp = output_spec.state_points
        if isinstance(sp, dict):
            for branch_type, values in sp.items():
                para = BRANCH_TYPE_TO_PARA_NAME.get(branch_type)
                if para is None:
                    continue
                # Build OR over the allowed values for this branch type
                # CLE2000 RPN: val1 var = val2 var = + val3 var = + ...
                # (sum > 0 means at least one matched)
                if len(values) == 1:
                    conditions.append(f"{para} {values[0]} =")
                else:
                    parts = []
                    for v in values:
                        parts.append(f"{para} {v} =")
                    # Combine with OR (+ in RPN, then check > 0)
                    combined = " ".join(parts)
                    combined += " " + "+ " * (len(parts) - 1)
                    combined += "0 >"
                    conditions.append(combined)

        if not conditions:
            return ""

        # AND all conditions together
        cond_str = " ".join(conditions)
        if len(conditions) > 1:
            cond_str += " " + "* " * (len(conditions) - 1)

        from .CLE2000 import wrap_cle2000_line as _wrap_guard
        return _wrap_guard(f"IF {cond_str} THEN") + "\n"

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the complete EDI/COMPO procedure to a ``.c2m`` file.

        When branches are present, the procedure receives additional
        REAL input parameters for the current statepoint values
        (e.g. ``TFuel``, ``TCool``, ``DCool``).

        Parameters
        ----------
        path_to_procs : str
            Directory path for the output file.
        proc_name : str
            Procedure name (without ``.c2m`` extension).

        Returns
        -------
        str
            Relative path of the written file.
        """
        # Build parameter variable recovery block for branch params
        branch_var_decl = ""
        if self.branches:
            para_names = [b.para_name for b in self.branches]
            branch_var_decl = (
                "REAL " + " ".join(para_names) + " ;\n"
                ":: " + " ".join(f">>{p}<<" for p in para_names)
                + " ;\n"
            )

        if self.branches:
            # COMPO is an input/output parameter from x2m
            header = (
                f"* PROCEDURE {proc_name}.c2m : calls to EDI: and COMPO: modules\n"
                "* --------------------------------\n"
                "* Procedure generated by starterDD\n"
                "* Author: R. Guasch\n"
                "* --------------------------------\n"
                "*    INPUT & OUTPUT PARAMETERS\n"
                "* --------------------------------\n"
                "PARAMETER COMPO FLUX LIBRARY2 TRACK ::\n"
                "::: LINKED_LIST COMPO ; \n"
                "::: LINKED_LIST FLUX ;\n"
                "::: LINKED_LIST LIBRARY2 ;\n"
                "::: LINKED_LIST TRACK ; ;\n"
                "STRING name_cpo ;\n"
                ":: >>name_cpo<< ;\n"
                f"{branch_var_decl}\n"
                "* --------------------------------\n"
                "*    MODULES DEFINITION\n"
                "* --------------------------------\n"
                "MODULE EDI: COMPO: DELETE: END: ;\n"
                "* --------------------------------\n"
                "*    LOCAL VARIABLES DEFINITION\n"
                "* --------------------------------\n"
                "LINKED_LIST EDIRATES ;\n"
            )
            header += "REAL spec_pow := 38.6 ; ! W/g\n"
            if "Burnup" in [b.para_name for b in self.branches]:
                header += "REAL TIME := Burnup spec_pow / ; \n"
        else:
            header = (
                f"* PROCEDURE {proc_name}.c2m : calls to EDI: and COMPO: modules\n"
                "* --------------------------------\n"
                "* Procedure generated by starterDD\n"
                "* Author: R. Guasch\n"
                "* --------------------------------\n"
                "*    INPUT & OUTPUT PARAMETERS\n"
                "* --------------------------------\n"
                "PARAMETER FLUX LIBRARY2 TRACK ::\n"
                "::: LINKED_LIST FLUX ;\n"
                "::: LINKED_LIST LIBRARY2 ;\n"
                "::: LINKED_LIST TRACK ; ;\n"
                "STRING name_cpo ;\n"
                ":: >>name_cpo<< ;\n"
                "* --------------------------------\n"
                "*    MODULES DEFINITION\n"
                "* --------------------------------\n"
                "MODULE EDI: COMPO: DELETE: END: ;\n"
                "* --------------------------------\n"
                "*    LOCAL VARIABLES DEFINITION\n"
                "* --------------------------------\n"
                "LINKED_LIST EDIRATES COMPO ;\n"
                "* --------------------------------\n"
                "*    COMPO FILE NAME FOR EXPORT\n"
                "* --------------------------------\n"
                "SEQ_ASCII _COMPO :: FILE <<name_cpo>> ;\n"
            )

        body = self.build_procedure_body()

        if self.branches:
            # When branches are present, COMPO init and export
            # are handled in the main x2m procedure.
            footer = "END: ;\n"
        else:
            footer = (
                "* --------------------------------\n"
                "*    EXPORT COMPO TO ASCII FILE\n"
                "* --------------------------------\n"
                "_COMPO := COMPO ;\n"
                "END: ;\n"
            )

        content = f"{header}{body}{footer}"

        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(path_to_procs, f"{proc_name}.c2m")
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[EDI_COMPO] Wrote procedure to {filepath}")
        return filepath


class MAC:
    """Generate DRAGON5 ``MAC:`` module calls for macroscopic cross-section libraries.

    Creates a macrolib from a collection of :class:`MaterialMixture` objects
    whose ``xs_data`` attribute carries the cross-section values (total,
    absorption, nu·fission, scattering matrix, etc.).

    The module supports two modes:

    * ``create_new=True`` (default): emit ``MACLIB := MAC: :: …``
      which creates a new macroscopic library.
    * ``create_new=False``: emit ``MACLIB := MAC: MACLIB :: …`` which
      updates an existing macrolib.

    Parameters
    ----------
    macro_lib_name : str
        CLE-2000 identifier for the macrolib object.
    create_new : bool
        If ``True`` (default), a fresh macrolib is created.  If
        ``False``, an existing one is updated.

    Attributes
    ----------
    material_mixtures : list of MaterialMixture
        Mixtures added via :meth:`add_material_mixture`.
    iprint : int
        DRAGON print level (default 1).
    ngroup : int
        Number of energy groups (default 1).
    anisotropy_level : int
        Legendre expansion order for scattering (default 0).
    """

    def __init__(self, macro_lib_name: str, create_new: bool = True):
        """
        MAC: module initialization.
        { MACLIB := MAC: [ MACLIB ] :: (descmacinp)
        | MICLIB := MAC: MICLIB :: (descmacinp)
        | MACLIB := MAC: [ MACLIB ] [ OLDLIB ] :: (descmacupd)
        | MACLIB := MAC: MACLIB OPTIM
        } 
        ;
        Purpose : create a macrolib based on a given collection of material mixtures
        :param macro_lib_name (str): Name of the macrolib used as identifier in DRAGON input file
        :param create_new (bool): Flag to indicate if a new macrolib should be created, if False, an existing MACROLIB can be edited
        """
        self.macro_lib_name = macro_lib_name
        self.create_new = create_new
        self.material_mixtures = []
        self.iprint = 1 # default print level
        self.ngroup = 1 # default number of energy groups
        self.anisotropy_level = 0 # default anisotropy level
        self.count_mixtures = 0



    def add_material_mixture(self, material_mixture):
        """
        Add a material mixture to the macrolib.
        
        :param material_mixture (MaterialMixture): MaterialMixture object to be added to the macrolib
        """
        self.material_mixtures.append(material_mixture)


    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the MAC module to a .c2m file for DRAGON processing.
        
        :param path_to_procs (str): Path to the directory where the .c2m file will be saved
        :param proc_name (str): Name of the process for which the .c2m file is being created
        """
        self.count_mixtures = len(self.material_mixtures)
        if self.count_mixtures == 0:
            raise ValueError("No material mixtures have been added to the MAC module.")
        filename = f"{path_to_procs}/{proc_name}.c2m"
        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        with open(filename, 'w') as file:
            if self.create_new:
                file.write(f"{self.macro_lib_name} := MAC:  ::\n")
            file.write(f"  EDIT {self.iprint}\n")
            file.write(f"  NGRO {self.ngroup}\n")
            file.write(f"  ANIS {self.anisotropy_level}\n")
            file.write(f"  NMIX {self.count_mixtures}\n")
            for mix in self.material_mixtures:
                file.write(f"! define {mix.material_name}\n")
                file.write(f"  MIX {mix.material_mixture_index}\n")
                for reaction, values in mix.xs_data.data.items():
                    if reaction == "scattering":
                        for i, row in enumerate(values):
                            row_str = ' '.join(f"{val:.8E}" for val in row)
                            j = values.index(row)
                            file.write(f"  SCAT {j+1} {i+1} {row_str}\n")
                    else:
                        values_str = ' '.join(f"{val:.8E}" for val in values)
                        file.write(f"  {reaction.upper()} {values_str}\n")
            
            file.write(";\n")

class USS:
    """
    Generate DRAGON5 ``USS:`` module calls for the self shielding calculation step.
    Support both IC and CP methods (ie ARM or PIJ kewords).
    Need to develop options for self-shielding options from USS:

    """

    def __init__(self, calculation_step, track_name=None, title=None):
        self.step = calculation_step
        # check step type is self-shielding
        if self.step.step_type != "self_shielding":
            raise ValueError(f"USS module only applicable for self-shielding steps, but got step type '{self.step.type}'")
        step_tag = self.step.name.upper()
        self.track_name = track_name or f"TRK{step_tag}"
        self.trkfil_name = f"TRKFIL{step_tag}"
        self.solution_method = "PIJ" if self.step.spatial_method == "CP" else "ARM"
        # Default names for libraries
        self.self_shielded_library_name = "LIBRARY2"
        self.library_name = "LIBRARY"

    def build_uss_call(self):
        """Return the full ``USS:`` call block as a string."""
        s = self.step
        lines = []
        lines.append(
            
            f"{self.self_shielded_library_name} := USS:"
        )
        lines.append(
            f" {self.track_name} {self.trkfil_name} {self.library_name} ::"
        )
        lines.append("    EDIT 1")
        lines.append(f"    {self.solution_method}")
        lines.append(";")
        return "\n".join(lines) + "\n"


class SALT:
    """
    Generate DRAGON5 ``SALT:`` module call blocks for tracking glow
    generated TDT geometries.

    Builds the correct SALT syntax depending on the spatial method
    (CP, IC, MOC) stored in the ``CalculationStep``.

    Parameters
    ----------
    calculation_step : CalculationStep
        Provides spatial_method, tracking, num_angles_2d,
        line_density, anisotropy_level, polar_angles_quadrature.
    tdt_var_name : str
        CLE-2000 variable name for the SEQ_ASCII TDT geometry.
    track_name : str or None
        Name for the LINKED_LIST tracking object.  Derived from
        ``step.name`` when ``None``.
    title : str or None
        Title string for the SALT call.  Auto-generated if ``None``.
    batch : int or None
        BATCH parameter.  Defaults to 200 for IC/CP, 2000 for MOC.
    """

    AVAILABLE_QUADRATURES = (
        "GAUS", "CACA", "CACB", "LCMD", "OPP1", "OGAU",
    )
    AVAILABLE_TSPC_ANGLES = (2, 6, 8, 12, 14, 18, 20, 24, 30)

    def __init__(self, calculation_step, tdt_var_name,
                 track_name=None, title=None, batch=None):
        self.step = calculation_step
        self.tdt_var_name = tdt_var_name

        step_tag = self.step.name.upper()
        self.track_name = track_name or f"TRK{step_tag}"
        self.trkfil_name = f"TRKFIL{step_tag}"
        self.title = title or (
            f"{step_tag} - {self.step.spatial_method}"
        )
        self.batch = calculation_step.batch_size

        if self.step.tracking == "TSPC" and self.step.num_angles_2d not in self.AVAILABLE_TSPC_ANGLES:
            raise ValueError(
                f"TSPC tracking only supports {self.AVAILABLE_TSPC_ANGLES} angles, "
                f"but got {self.step.num_angles_2d}."
            )

    def build_salt_call(self):
        """Return the full ``SALT:`` call block as a string."""
        s = self.step
        lines = []
        lines.append(
            f"{self.track_name} {self.trkfil_name} "
            f":= SALT: {self.tdt_var_name} ::"
        )
        lines.append("    EDIT 2")
        lines.append(f"    TITLE '{self.title}'")
        if self.batch is not None:
            lines.append(f"    BATCH {self.batch}")
        lines.append(f"    ANIS {s.anisotropy_level}")
        
        if s.spatial_method == "MOC":
            #lines.append(f"    {s.polar_angles_quadrature}")    
            lines.append("    ALLG")
            lines.append(
                f"    {s.tracking} "
                f"{s.num_angles_2d} "
                f"{s.line_density:.1f} REND"
            )
            lines.append("    NOIC")
        elif s.spatial_method == "IC":
            lines.append(
                f"    {s.tracking} "
                f"{s.num_angles_2d} "
                f"{s.line_density:.1f}"
            )
            lines.append("    IC EPSJ 1.0E-5")
        else:
            # CP
            lines.append(
                f"    {s.tracking} "
                f"{s.num_angles_2d} "
                f"{s.line_density:.1f}"
            )
            lines.append("    NOIC")

        lines.append(";")
        return "\n".join(lines) + "\n"


class MCCGT:
    """
    Generate DRAGON5 ``MCCGT:`` module call for MOC tracking.

    Only applicable when ``spatial_method == "MOC"``.

    Parameters
    ----------
    calculation_step : CalculationStep
        Must have ``spatial_method == "MOC"``.
    track_name : str or None
        LINKED_LIST tracking object name.  Defaults to
        ``TRK{step.name}``.
    max_inner_iterations : int
        AAC iteration limit (default 100).
    krylov_dim : int
        Krylov subspace dimension (default 10).
    """

    def __init__(self, calculation_step, track_name=None,
                 max_inner_iterations=100, krylov_dim=10):
        if calculation_step.spatial_method != "MOC":
            raise ValueError(
                "MCCGT is only applicable for MOC steps."
            )
        self.step = calculation_step
        step_tag = self.step.name.upper()
        self.track_name = track_name or f"TRK{step_tag}"
        self.trkfil_name = f"TRKFIL{step_tag}"
        self.max_inner_iterations = max_inner_iterations
        self.nmu = 4  # default number of polar angles for quadrature
        if self.step.anisotropy_level > 4:
            self.nmu = 6

    def build_mccgt_call(self):
        """Return the full ``MCCGT:`` call block as a string."""
        s = self.step
        lines = []
        lines.append(
            f"{self.track_name} := MCCGT: "
            f"{self.track_name} {self.trkfil_name} ::"
        )
        lines.append("    EDIT 1")
        lines.append(
            f"    {s.polar_angles_quadrature} {self.nmu} "
            f"AAC 80 TMT EPSI 1E-5"
        )
        lines.append(";")
        return "\n".join(lines) + "\n"


class TRK:
    """
    Generate a ``TRK.c2m`` sub-procedure that contains all SALT:
    (and MCCGT:) calls for every step in a calculation scheme.

    The procedure receives TDT SEQ_ASCII variables and the anisotropy
    level as input parameters and returns the LINKED_LIST tracking
    objects and their SEQ_BINARY companion files.

    Parameters
    ----------
    scheme : DragonCalculationScheme
        Calculation scheme whose steps define the tracking calls.
    case_name : str
        Used to derive TDT file names.
    """

    def __init__(self, scheme, case_name):
        self.scheme = scheme
        self.case_name = case_name
        self._salt_objects = []
        self._mccgt_objects = []
        self._build_tracking_objects()

    def _build_tracking_objects(self):
        """Create SALT (and MCCGT) objects for each trackable step."""
        for step in self.scheme.get_trackable_steps():
            if step.step_type == "self_shielding":
                self.problem_anisotropy_level = step.anisotropy_level
            else:
                step.anisotropy_level = getattr(self, "problem_anisotropy_level", 1)
            tdt_var = self._tdt_var_name(step)
            salt = SALT(step, tdt_var)
            self._salt_objects.append(salt)
            if step.spatial_method == "MOC":
                mccgt = MCCGT(step)
                self._mccgt_objects.append(mccgt)

    def _tdt_var_name(self, step):
        """CLE-2000 variable name for the TDT file import."""
        if len(step.name) > 9:
            tag = step.name[:9].lower()
        else:
            tag = step.name.lower()
        return f"tdt{tag}"

    def _tdt_file_name(self, step):
        """Physical TDT file name as produced by glow."""
        suffix = "_MACRO" if step.export_macros else ""
        return (
            f"{self.case_name}_{step.name}"
            f"_{step.spatial_method}"
            f"_{step.tracking}{suffix}.dat"
        )

    def get_tdt_file_mapping(self):
        """
        ``{cle2000_var_name: physical_file_name}`` for every trackable step.
        Used by the main x2m to declare SEQ_ASCII imports.
        """
        return {
            self._tdt_var_name(s): self._tdt_file_name(s)
            for s in self.scheme.get_trackable_steps()
        }

    def get_track_names(self):
        """
        ``[(track_ll_name, trkfil_name)]`` per trackable step (ordered).
        """
        result = []
        for step in self.scheme.get_trackable_steps():
            if len(step.name) > 9:
                tagTRK = step.name[:9].upper()
            else:
                tagTRK = step.name.upper()
            if len(step.name) > 6:
                tagTRKFIL = step.name[:6].upper()
            else:                
                tagTRKFIL = step.name.upper()
            result.append((f"TRK{tagTRK}", f"TRKFIL{tagTRKFIL}"))
            print(result)
        return result

    def build_procedure_body(self):
        """Build the body of the TRK.c2m procedure."""
        from .CLE2000 import CLE2000_MAX_LINE, wrap_cle2000_line
        body = ""
        for salt in self._salt_objects:
            body += (
                f"* Tracking for step {salt.step.name}\n"
            )
            body += salt.build_salt_call()
            body += "\n"

        for mccgt in self._mccgt_objects:
            body += (
                f"* MOC tracking for step "
                f"{mccgt.step.name}\n"
            )
            body += mccgt.build_mccgt_call()
            body += "\n"

        return body

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the complete TRK sub-procedure to a ``.c2m`` file.
        """
        from .CLE2000 import (
            CLE2000_MAX_LINE, CLE2000_MAX_VARNAME,
            validate_varname, wrap_cle2000_line,
        )

        track_names = self.get_track_names()
        tdt_map = self.get_tdt_file_mapping()

        # --- PARAMETER block ---
        param_items = []
        for trk_ll, trkfil in track_names:
            print(f"Validating track names: {trk_ll}, {trkfil}")
            validate_varname(trk_ll)
            validate_varname(trkfil)
            param_items.append(trk_ll)
            param_items.append(trkfil)
        for tdt_var in tdt_map:
            validate_varname(tdt_var)
            param_items.append(tdt_var)

        header = (
            f"* PROCEDURE {proc_name}.c2m : tracking\n"
            "* --------------------------------\n"
            "* Procedure generated by starterDD\n"
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
        )

        # PARAMETER declaration
        param_block = "PARAMETER"
        for item in param_items:
            param_block += f" {item}"
        param_block += " ::\n"
        # Linked-list declarations
        for trk_ll, trk_bin in track_names:
            param_block += (
                f"::: LINKED_LIST {trk_ll} ;\n"
                f"::: SEQ_BINARY {trk_bin} ;\n"
            )
        # TDT vars are SEQ_ASCII (no linked-list decl)
        param_block += ";\n"

        # MODULE declaration
        mod_block = "MODULE SALT: MCCGT: END: ;\n"

        body = self.build_procedure_body()

        footer = "END: ;\nQUIT .\n"

        content = (
            f"{header}{param_block}\n"
            f"{mod_block}\n{body}{footer}"
        )

        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        filepath = os.path.join(
            path_to_procs, f"{proc_name}.c2m"
        )
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[TRK] Wrote procedure to {filepath}")
        return filepath


class ASM:
    """
    Defintion of the ASM: module calls for Dragon calculations
    Given a calculation step, determine the appropriate ASM: call and generate the corresponding c2m procedure.
    """

    def __init__(self, calculation_step, track_name=None, title=None):
        self.step = calculation_step
        step_tag = self.step.name.upper()
        self.track_name = track_name or f"TRK{step_tag}"
        self.solution_method = "PIJ" if self.step.spatial_method == "CP" else "ARM"
        

    def build_asm_call(self):
        """Return the full ``ASM:`` call block as a string."""
        s = self.step
        lines = []
        lines.append(
            f"ASM: {self.track_name} ::"
        )
        lines.append("    EDIT 1")
        lines.append(f"    {self.solution_method}")
        lines.append(";")
        return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# EDI_condensation – energy-group condensation between flux levels
# ---------------------------------------------------------------------------

class EDI_condensation:
    """
    Generate the ``EDI:`` call that condenses cross sections from the
    fine energy mesh to a coarser one between two flux levels.

    If lib_name is not specified :
    Produces two CLE-2000 lines:

    1. ``EDITION := EDI: <flux> <lib> <trk> :: EDIT 0 MICR ALL MERG MIX
       COND <bounds> SAVE ON COND<N> ;``
    2. ``LIB<N>G := EDITION :: STEP UP COND<N> ;``

    Parameters
    ----------
    edition_step : EditionBetweenLevelsStep
        Carries ``number_of_macro_groups`` and ``energy_groups_bounds``.
    """

    def __init__(self, edition_step, lib_name=None):
        self.step = edition_step
        n = edition_step.number_of_macro_groups
        if lib_name is None:
            self.lib_name = f"LIB{n}G"
        else:
            self.lib_name = lib_name
        self.cond_name = f"COND{n}"

    def build_edi_call(self, flux_ll, lib_ll, trk_ll):
        """Return the ``EDI:`` call block as a string.

        Parameters
        ----------
        flux_ll : str
            CLE-2000 name of the FLUX linked list (e.g. ``"FLUX"``).
        lib_ll : str
            CLE-2000 name of the self-shielded library (e.g. ``"LIBRARY2"``).
        trk_ll : str
            CLE-2000 name of the L1 tracking linked list.
        """
        groups_str = "  ".join(
            str(g) for g in self.step.energy_groups_bounds
        )
        lines = [
            f"EDITION := EDI: {flux_ll} {lib_ll} {trk_ll} ::",
            "    EDIT 0",
            "    MICR ALL",
            "    MERG MIX",
            f"    COND  {groups_str}",
            f"    SAVE ON {self.cond_name}",
            ";",
        ]
        return "\n".join(lines) + "\n"

    def build_lib_extract(self):
        """Return the ``STEP UP`` extraction line as a string."""
        return (
            f"{self.lib_name} := EDITION :: "
            f"STEP UP {self.cond_name} ;\n"
        )


# ---------------------------------------------------------------------------
# SPH_correction – SPH equivalence correction on the condensed library
# ---------------------------------------------------------------------------

class SPH_correction:
    """
    Generate the ``SPH:`` equivalence correction call applied to the
    condensed library after a cross-section condensation step.

    If lib_name is not specified,
    Produces:
      ``LIB<N>G := SPH: LIB<N>G <trk> <trkfil> :: EDIT 1 GRMAX <n> ;``
    else: uses the provided lib_name as a CLE-2000 variable

    Parameters
    ----------
    edition_step : EditionBetweenLevelsStep
        Carries ``number_of_macro_groups`` and ``max_sph_group``.
    """

    def __init__(self, edition_step, lib_name=None):
        self.step = edition_step
        n = edition_step.number_of_macro_groups
        if lib_name is None:
            self.lib_name = f"LIB{n}G"
        else:
            self.lib_name = lib_name

    def build_sph_call(self, trk_ll, trkfil_ll):
        """Return the ``SPH:`` call block as a string.

        Parameters
        ----------
        trk_ll : str
            CLE-2000 name of the L1 tracking linked list.
        trkfil_ll : str
            CLE-2000 name of the L1 binary tracking file (SEQ_BINARY).
        """
        lines = [
            f"{self.lib_name} := SPH: {self.lib_name} "
            f"{trk_ll} {trkfil_ll} ::",
            "    EDIT 1",
        ]
        if self.step.max_sph_group is not None:
            lines.append(f"    GRMAX {self.step.max_sph_group}")
        lines.append(";")
        return "\n".join(lines) + "\n"


class MIXEQ:
    """
    Class for generating DRAGON5 MIXEQ procedure files (.c2m).

    When calculation steps use different mix numbering strategies,
    DRAGON requires a MIXEQ module to map mixtures between steps.
    This class generates the required correspondence using library edit mode.

    The MIXEQ procedure duplicates mixtures from an initial "by material"
    library to a comprehensive "by pin" or other numbering scheme using
    the LIB: module in library edit mode.

    Usage
    -----
    ::
        # Assuming assembly has steps with different strategies recorded
        mixeq = MIXEQ(assembly, lib_name, "SSH", "FLUX_L2")
        mixeq.write_to_c2m("./procs", "MIXEQ_SSH_to_L2")
    """

    def __init__(self, assembly_model, input_lib_name, output_lib_name, from_step, to_step, draglib_alias):
        """
        Initialize MIXEQ procedure generator.

        Parameters
        ----------
        assembly_model : CartesianAssemblyModel
            Assembly with mix state history from multiple steps.
            Must have run different numbering strategies and recorded
            the state history for both from_step and to_step.
        input_lib_name : str
            CLE-2000 variable name for the input library (e.g., "LIBRARY2)
        output_lib_name : str
            CLE-2000 variable name for the output library (e.g., "LIBEQ")
        from_step : str
            Source step name (e.g., "SSH") - typically uses by_material strategy
        to_step : str
            Target step name (e.g., "FLUX_L2") - typically uses by_pin strategy
        draglib_alias : str
            Alias for the draglib to be used in MIXEQ procedures : used to recover depletion chain information.
        """
        self.assembly = assembly_model
        self.input_lib_name = input_lib_name
        self.output_lib_name = output_lib_name
        self.from_step = from_step
        self.to_step = to_step
        self.draglib_alias = draglib_alias

        # Get correspondence table using existing method
        try:
            self.correspondence_table = assembly_model.build_mixeq_correspondence_table(
                from_step, to_step
            )
        except AttributeError:
            raise RuntimeError(
                f"Assembly model does not have build_mixeq_correspondence_table method. "
                f"This indicates the assembly may not support MIXEQ generation."
            )

        # Validate that correspondence table is not empty
        if not self.correspondence_table:
            raise ValueError(
                f"No correspondence found between steps '{from_step}' and '{to_step}'. "
                f"Ensure both steps have been processed with different numbering strategies."
            )

    def _determine_max_mix_count(self):
        """
        Calculate maximum mixture index needed for the library.

        Returns
        -------
        int
            Maximum mixture index from correspondence table
        """
        if not self.correspondence_table:
            return 0

        # Get max index from correspondence table
        max_from = max(entry[2] for entry in self.correspondence_table)  # from_idx
        max_to = max(entry[3] for entry in self.correspondence_table)    # to_idx
        return max(max_from, max_to)

    def _validate_correspondence_completeness(self):
        """
        Validate that correspondence table covers required mixtures.

        Returns
        -------
        bool
            True if validation passes

        Raises
        ------
        ValueError
            If validation fails
        """
        if not self.correspondence_table:
            raise ValueError(f"Empty correspondence table for {self.from_step} → {self.to_step}")

        # Check if step mix states exist in assembly history
        if hasattr(self.assembly, 'mix_state_history'):
            from_state = self.assembly.get_step_mix_state(self.from_step)
            to_state = self.assembly.get_step_mix_state(self.to_step)

            if not from_state:
                raise ValueError(f"No mix state found for step '{self.from_step}'")
            if not to_state:
                raise ValueError(f"No mix state found for step '{self.to_step}'")

        return True

    def build_lib_module_call(self):
        """
        Generate LIB: module call in library edit mode for mixture duplication.

        Returns
        -------
        str
            Complete LIB: module call with COMB statements for mixture duplication
        """
        if not self.correspondence_table:
            raise ValueError(f"No correspondence found between {self.from_step} and {self.to_step}")

        # Validate correspondence before building
        self._validate_correspondence_completeness()

        max_mix = self._determine_max_mix_count()

        # Start LIB: call in library edit mode
        call = (
            f"{self.output_lib_name} := LIB: {self.input_lib_name} ::\n" # Library edit mode - read existing microlib
            "  EDIT 0\n"  
            f"  NMIX {max_mix}\n"
            f"  DEPL LIB: DRAGON FIL: {self.draglib_alias}\n"
            "   CATL\n"
            "\n"
        )

        # Generate MIX statements for mixture duplication and copying into the new library
        # Non-fuel isotopes are now included in the correspondence table via build_mixeq_correspondence_table()

        # Sort by target index to ensure consistent output
        sorted_table = sorted(self.correspondence_table, key=lambda x: x[3])

        for from_name, to_name, from_idx, to_idx in sorted_table:
            call += f"  ! {from_name} → {to_name}\n"
            call += f"  MIX {to_idx} {from_idx}\n"
            call += "\n"

        call += ";"
        return call

    def _validate_max_mix_count(self, expected_max=None):
        """
        Validate max mix count against expected value.

        Parameters
        ----------
        expected_max : int, optional
            Expected maximum mixture count. If provided, raises error if
            calculated max exceeds this value.

        Returns
        -------
        int
            Calculated maximum mixture count
        """
        calculated_max = self._determine_max_mix_count()

        if expected_max and calculated_max > expected_max:
            raise ValueError(
                f"Calculated max mix count ({calculated_max}) exceeds "
                f"expected maximum ({expected_max}) for transition "
                f"{self.from_step} → {self.to_step}"
            )

        return calculated_max

    def build_mix_index_comment_block(self):
        """
        Build comment block showing mixture correspondence mapping.

        Returns
        -------
        str
            Comment block with mixture mapping information
        """
        if not self.correspondence_table:
            return "! No mixture correspondences found\n"

        comment = (
            f"! --------------------------------\n"
            f"! MIXTURE CORRESPONDENCE MAPPING\n"
            f"! {self.from_step} → {self.to_step}\n"
            f"! --------------------------------\n"
        )

        for from_name, to_name, from_idx, to_idx in self.correspondence_table:
            comment += f"! MIX {from_idx:<3} → MIX {to_idx:<3} ({from_name} → {to_name})\n"

        comment += f"! --------------------------------\n"
        return comment

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write complete MIXEQ procedure to .c2m file.

        The generated file is a standalone CLE-2000 procedure that receives
        an input microlib and returns an enhanced microlib with additional
        mixture duplications.

        Parameters
        ----------
        path_to_procs : str
            Directory path for the output file.
        proc_name : str
            Procedure name (without .c2m extension).

        Returns
        -------
        str
            Relative path of the written file.
        """
        from .CLE2000 import validate_varname

        # Validate procedure name (CLE-2000 limitation)
        validate_varname(proc_name[:12])

        # Validate correspondence before writing
        self._validate_correspondence_completeness()

        header = (
            f"*PROCEDURE {proc_name}.c2m\n"
            "* --------------------------------\n"
            "* MIXEQ procedure generated by starterDD\n"
            "* Mixture duplication for numbering strategy transitions\n"
            "* Author: R. Guasch\n"
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
            f"PARAMETER {self.input_lib_name} {self.output_lib_name} ::\n"
            "       EDIT 0\n"
            f"           ::: LINKED_LIST {self.input_lib_name} {self.output_lib_name} ;\n"
            "   ;\n"
            "\n"
            "* --------------------------------\n"
            "*    STRUCTURES AND MODULES\n"
            "* --------------------------------\n"
            "MODULE LIB: END: ;\n"
            "\n"
        )

        mix_comment = self.build_mix_index_comment_block()
        lib_call = self.build_lib_module_call()

        footer = (
            "* -----------------------------------------\n"
            "*         END OF MIXEQ PROCEDURE\n"
            "* -----------------------------------------\n"
            "END: ;\n"
            "QUIT ."
        )

        content = (
            f"{header}"
            f"{mix_comment}"
            "\n"
            "* --------------------------------\n"
            "*    MIXTURE DUPLICATION\n"
            "*    CALL TO LIB: MODULE\n"
            "* --------------------------------\n"
            f"{lib_call}\n"
            f"{footer}"
        )

        # Ensure output directory exists
        if path_to_procs and not os.path.exists(path_to_procs):
            os.makedirs(path_to_procs)

        # Write file
        filepath = os.path.join(path_to_procs, f"{proc_name}.c2m")
        with open(filepath, 'w') as f:
            f.write(content)

        print(f"[MIXEQ] Wrote procedure to {filepath}")
        print(f"[MIXEQ] Generated {len(self.correspondence_table)} mixture correspondences")
        print(f"[MIXEQ] Max mixture index: {self._determine_max_mix_count()}")

        return filepath