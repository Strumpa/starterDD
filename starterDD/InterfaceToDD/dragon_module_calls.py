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
        lib.set_isotope_alias("MODERATOR", "H1", "H1_H2O")
        lib.set_isotope_alias("COOLANT",   "H1", "H1_H2O")
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
            Absolute path of the written file.
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
        for extra in self._extra_non_fuel_mixes:
            all_indices.append(extra["index"])
        return max(all_indices) if all_indices else 0

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
