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

    # Heavy-metal + Gd isotopes that receive self-shielding treatment (INRS)
    DEFAULT_SELF_SHIELDED_FUEL_ISOTOPES = [
        "U234", "U235", "U236", "U238",
        "Pu239", "Pu240", "Pu241", "Pu242",
        "Gd154", "Gd155", "Gd156", "Gd157", "Gd158", "Gd160",
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

    def __init__(self, assembly_model):
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
        """
        self.assembly = assembly_model

        # Fuel self-shielding configuration
        self.self_shielded_fuel_isotopes = list(self.DEFAULT_SELF_SHIELDED_FUEL_ISOTOPES)
        self.fuel_inrs = 1  # INRS value for fuel self-shielding group

        # Non-fuel self-shielding: {material_name: {isotope_name: inrs_value}}
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

    # ------------------------------------------------------------------
    #  Configuration helpers
    # ------------------------------------------------------------------

    def set_self_shielded_fuel_isotopes(self, isotopes):
        """Override the list of self-shielded fuel isotopes (default: HM + Gd)."""
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
                line = f"    {isotope} = {lib_name} {density:.5E}"
                if isotope in self.correlation_isotopes:
                    line += " CORR 1"
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

        lines = ""
        for mat_name, idx, composition, temp_var in entries:
            iso_comp = composition.get_isotope_name_composition()
            inrs_map = self.non_fuel_inrs.get(mat_name, {})
            depletable = getattr(composition, 'depletable', False)
            if depletable:
                lines += f"    MIX {idx} <<{temp_var}>>\n"
            else:
                lines += f"    MIX {idx} <<{temp_var}>> NOEV\n"
            for isotope, density in iso_comp.items():
                lib_name = self.isotope_aliases.get(
                    (mat_name, isotope), isotope
                )
                line = f"    {isotope} = {lib_name} {density:.5E}"
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
            inrs_map = self.non_fuel_inrs.get(mix.material_name, {})
            for isotope, density in iso_comp.items():
                lib_name = self.isotope_aliases.get(
                    (mix.material_name, isotope), isotope
                )
                line = f"    {isotope} = {lib_name} {density:.5E}"
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

        call = (
            "LIBRARY := LIB: ::\n"
            "EDIT 0\n"
            f"NMIX {max_mix}  ! MAXIMUM OF MATERIAL MIXTURES\n"
            "<<ssh_method>>\n"
            "ANIS <<anis_level>>\n"
            "CTRA <<tran_correc>>\n"
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
            "DOUBLE DTFUEL DTBOX DTCLAD DTCOOL DTMODE ;\n"
            ":: >>DTFUEL<< >>DTBOX<< >>DTCLAD<< >>DTCOOL<< "
            ">>DTMODE<< ; ! Temperatures\n"
            "\n"
            "* --------------------------------------------\n"
            "*  CONVERT DOUBLE TO REALS for TEMPERATURES\n"
            "* --------------------------------------------\n"
            "REAL TFUEL := DTFUEL D_TO_R ;\n"
            "REAL TBOX := DTBOX D_TO_R ;\n"
            "REAL TCLAD := DTCLAD D_TO_R ;\n"
            "REAL TCOOL := DTCOOL D_TO_R ;\n"
            "REAL TMODE := DTMODE D_TO_R ;\n"
            "\n"
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
                # Parse pin_idx from "UOX28_zone1_pin3"
                pin_idx = int(
                    mix.unique_material_mixture_name.rsplit("_pin", 1)[1]
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
        iso_str = " ".join(self.isotopes)
        vector = self.build_merg_mix_vector()
        merg_mix = self._format_merg_mix_vector(vector)
        cond = self.build_cond_option()

        call = (
            "EDIRATES := EDI: FLUX LIBRARY2 TRACK ::\n"
            f"   EDIT {self.edit_level}\n"
            f"  MICR {len(self.isotopes)} {iso_str}\n"
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

    def __init__(self):
        self.directories = []  # list of (name, comment, isotopes)

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

        Returns
        -------
        str
        """
        lines = "COMPO := COMPO: ::\n"
        for name, comment, isotopes in self.directories:
            iso_str = " ".join(isotopes)
            lines += f"   STEP UP '{name}'\n"
            lines += f"       COMM '{comment}' ENDC\n"
            lines += f"      ISOT {len(isotopes)} {iso_str}\n"
            lines += "   INIT\n"
        lines += ";\n"
        return lines

    def build_compo_store(self, directory_name):
        """
        Build a ``COMPO := COMPO: COMPO EDIRATES LIBRARY2 :: STEP UP <name> ;``
        store block.

        Parameters
        ----------
        directory_name : str
            Must match a name previously registered via ``add_directory``.

        Returns
        -------
        str
        """
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

    def __init__(self, assembly_model):
        """
        Parameters
        ----------
        assembly_model : CartesianAssemblyModel
            Assembly model with numbered material mixtures and (optionally)
            TDT-enforced indices.
        """
        self.assembly = assembly_model
        self.compo = COMPO()
        self.editions = []  # list of EDI objects

    def add_edition(self, name, comment, isotopes, spatial_mode,
                    energy_bounds=None, custom_map=None):
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
        """
        edi = EDI(name, self.assembly)
        edi.set_isotopes(isotopes)
        edi.set_spatial_homogenization(spatial_mode, custom_map=custom_map)
        edi.set_energy_condensation(energy_bounds)
        self.editions.append(edi)
        self.compo.add_directory(name, comment, isotopes)

    def build_procedure_body(self):
        """
        Build the body of the CLE-2000 procedure (everything between header
        and footer).

        Returns
        -------
        str
        """
        body = ""

        # COMPO initialization
        body += "* --------------------------------\n"
        body += "*   COMPO INITIALIZATION\n"
        body += "* --------------------------------\n"
        body += self.compo.build_compo_init()

        # EDI + COMPO pairs
        for edi in self.editions:
            body += "* --------------------------------\n"
            body += f"*    EDI: CALL FOR {edi.name}\n"
            body += "* --------------------------------\n"
            body += edi.build_edi_call()
            body += "* --------------------------------\n"
            body += f"*    COMPO: CALL FOR {edi.name}\n"
            body += "* --------------------------------\n"
            body += self.compo.build_compo_store(edi.name)
            body += "EDIRATES := DELETE: EDIRATES ;\n"

        return body

    def write_to_c2m(self, path_to_procs, proc_name):
        """
        Write the complete EDI/COMPO procedure to a ``.c2m`` file.

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
        header = (
            f"* PROCEDURE {proc_name}.c2m : calls to EDI: and COMPO: modules\n"
            "* --------------------------------\n"
            "* Procedure generated by starterDD\n"
            "* Author: R. Guasch\n"
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
            "PARAMETER COMPO FLUX LIBRARY2 TRACK ::\n"
            "::: LINKED_LIST COMPO ;\n"
            "::: LINKED_LIST FLUX ;\n"
            "::: LINKED_LIST LIBRARY2 ;\n"
            "::: LINKED_LIST TRACK ; ;\n"
            "STRING name_cpo save_opt ;\n"
            ":: >>name_cpo<< >>save_opt<< ; "
            "! save option for COMPO: module, e.g. 'SAVE' or 'NOSAVE'\n"
            "* --------------------------------\n"
            "*    MODULES DEFINITION\n"
            "* --------------------------------\n"
            "MODULE EDI: COMPO: DELETE: END: ;\n"
            "* --------------------------------\n"
            "*    LOCAL VARIABLES DEFINITION\n"
            "* --------------------------------\n"
            "LINKED_LIST EDIRATES ;\n"
            "* --------------------------------\n"
            "*    COMPO FILE NAME FOR EXPORT\n"
            "* --------------------------------\n"
            "SEQ_ASCII _COMPO :: FILE <<name_cpo>> ;\n"
        )

        body = self.build_procedure_body()

        footer = (
            "* --------------------------------\n"
            "IF save_opt 'SAVE' = THEN\n"
            "   _COMPO := COMPO ;\n"
            "ENDIF ;\n"
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
