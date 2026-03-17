# Collection of functions to orchestrate CLE2000 procedure generation
# CLE2000 language constraints and procedure classes.
# R.Guasch

import os
import textwrap

# ---------------------------------------------------------------
# CLE2000 language constraints
# ---------------------------------------------------------------
CLE2000_MAX_LINE = 72      # max characters per source line
CLE2000_MAX_VARNAME = 12   # max characters for variable names

COMMON_VARIABLE_TYPES = {
    "number_angles": "INTEGER",
    "line_density": "REAL",
    "temperature": "REAL",
    "density": "REAL",
    "draglib_name": "STRING",
}


def validate_varname(name):
    """Raise ``ValueError`` if *name* exceeds the CLE2000 limit."""
    if len(name) > CLE2000_MAX_VARNAME:
        raise ValueError(
            f"CLE2000 variable name '{name}' is "
            f"{len(name)} characters long "
            f"(max {CLE2000_MAX_VARNAME})."
        )


def wrap_cle2000_line(line, indent=2):
    """Wrap a single logical line to respect ``CLE2000_MAX_LINE``.

    Continuation is indicated by placing the remaining tokens on
    the next line(s) indented by *indent* spaces.  Already
    short-enough lines are returned unchanged.
    """
    if len(line) <= CLE2000_MAX_LINE:
        return line
    # Split on whitespace; reassemble keeping under the limit
    tokens = line.split()
    lines = []
    current = tokens[0]
    for tok in tokens[1:]:
        if len(current) + 1 + len(tok) > CLE2000_MAX_LINE:
            lines.append(current)
            current = " " * indent + tok
        else:
            current += " " + tok
    lines.append(current)
    return "\n".join(lines)


def _wrap_write(f, line):
    """Write *line* through ``wrap_cle2000_line`` + newline."""
    f.write(wrap_cle2000_line(line) + "\n")


# ---------------------------------------------------------------
# Base CLE2000 procedure
# ---------------------------------------------------------------
class CLE2000_procedure:
    def __init__(self, procedure_type, procedure_name,
                 parameters_structures=None,
                 parameters_variables=None):
        self.procedure_type = procedure_type
        self.procedure_name = procedure_name
        self.parameters_structures = parameters_structures or []
        self.parameters_variables = parameters_variables or []
        self.module_calls = []
        self.external_imports = []
        self.variables = []        # (name, type, default)
        self.linked_lists = []     # [name, ...]
        self.seq_binaries = []     # [name, ...]
        self.seq_asciis = []       # [(var, filename), ...]
        self.body_lines = []       # arbitrary CLE2000 lines

    # --- variable / structure declarations -------------------------

    def add_variable(self, name, var_type, default=None):
        validate_varname(name)
        if var_type not in ("STRING", "INTEGER", "REAL"):
            raise ValueError(
                f"Unknown CLE2000 type '{var_type}'."
            )
        self.variables.append((name, var_type, default))

    def add_linked_list(self, name):
        validate_varname(name)
        self.linked_lists.append(name)

    def add_seq_binary(self, name):
        validate_varname(name)
        self.seq_binaries.append(name)

    def add_seq_ascii(self, var_name, filename):
        validate_varname(var_name)
        self.seq_asciis.append((var_name, filename))

    def add_body_line(self, line):
        """Append an arbitrary CLE2000 source line."""
        self.body_lines.append(line)

    def add_body_block(self, block):
        """Append a multi-line CLE2000 block (string)."""
        for ln in block.splitlines():
            self.body_lines.append(ln)

    # --- external imports -------------------------------------------

    def add_external_import(self, imported_file_type,
                            cle2000_var_name,
                            imported_file_name):
        validate_varname(cle2000_var_name)
        self.external_imports.append(
            (imported_file_type, cle2000_var_name,
             imported_file_name)
        )

    # --- module calls -----------------------------------------------

    def add_module_call(self, module_name,
                        module_parameters=None):
        self.module_calls.append(
            (module_name, module_parameters)
        )

    # --- build helpers (return strings) ----------------------------

    def build_module_declaration(self):
        if not self.module_calls:
            return ""
        names = " ".join(
            f"{m}:" for m, _ in self.module_calls
        )
        line = f"MODULE {names} ;"
        return wrap_cle2000_line(line) + "\n"

    def build_variable_declarations(self):
        # Group by type for compact output
        by_type = {"STRING": [], "INTEGER": [], "REAL": []}
        for name, vt, default in self.variables:
            by_type[vt].append((name, default))
        out = ""
        for vt in ("STRING", "INTEGER", "REAL"):
            items = by_type[vt]
            if not items:
                continue
            for name, default in items:
                if default is not None:
                    if vt == "STRING":
                        ln = f'{vt} {name} := "{default}" ;'
                    else:
                        ln = f"{vt} {name} := {default} ;"
                else:
                    ln = f"{vt} {name} ;"
                out += wrap_cle2000_line(ln) + "\n"
        return out

    def build_linked_list_declaration(self):
        if not self.linked_lists:
            return ""
        line = (
            "LINKED_LIST "
            + " ".join(self.linked_lists) + " ;"
        )
        return wrap_cle2000_line(line) + "\n"

    def build_seq_binary_declaration(self):
        if not self.seq_binaries:
            return ""
        line = (
            "SEQ_BINARY "
            + " ".join(self.seq_binaries) + " ;"
        )
        return wrap_cle2000_line(line) + "\n"

    def build_seq_ascii_declarations(self):
        out = ""
        for var, fname in self.seq_asciis:
            ln = (
                f"SEQ_ASCII {var}"
                f" :: FILE '{fname}' ;"
            )
            out += wrap_cle2000_line(ln) + "\n"
        return out

    def build_external_imports_declaration(self):
        if not self.external_imports:
            return ""
        out = ""
        for ftype, var, fname in self.external_imports:
            ln = f"{ftype} {var} :: FILE '{fname}' ;"
            out += wrap_cle2000_line(ln) + "\n"
        return out


# ---------------------------------------------------------------
# Main procedure (.x2m)
# ---------------------------------------------------------------
class main_procedure(CLE2000_procedure):
    def __init__(self, procedure_name,
                 parameters_structures=None,
                 parameters_variables=None):
        super().__init__(
            "main", procedure_name,
            parameters_structures, parameters_variables,
        )
        self.sub_procedures = []
        self.enable_g2s = False

    def add_sub_procedure(self, sub_proc):
        self.sub_procedures.append(sub_proc)

    def build_procedure_declaration(self):
        if not self.sub_procedures:
            return ""
        names = " ".join(
            sp.procedure_name for sp in self.sub_procedures
        )
        line = f"PROCEDURE {names} ;"
        return wrap_cle2000_line(line) + "\n"

    def build_g2s_block(self, tdt_map):
        """Optional G2S: visualization calls.

        Parameters
        ----------
        tdt_map : dict
            ``{tdt_cle_var: tdt_filename}`` — one G2S per TDT.
        """
        if not self.enable_g2s:
            return ""
        out = (
            "*" * 55 + "\n"
            "* Geometry visualization with G2S:\n"
            "*" * 55 + "\n"
        )
        for var, fname in tdt_map.items():
            ps_var = f"{var}.ps"
            ps_file = fname.replace(".dat", ".ps")
            validate_varname(ps_var)
            out += wrap_cle2000_line(
                f"SEQ_ASCII {ps_var}"
                f" :: FILE '{ps_file}' ;"
            ) + "\n"
            out += f"{ps_var} := G2S: {var} :: ;\n"
        out += "\n"
        return out

    def write_to_x2m(self, output_path):
        """Write the complete .x2m procedure file."""
        if output_path and not os.path.exists(output_path):
            os.makedirs(output_path)
        filepath = os.path.join(
            output_path, f"{self.procedure_name}.x2m"
        )

        with open(filepath, 'w') as f:
            # Header comment
            f.write(
                f"* CLE2000 procedure {self.procedure_name}\n"
                "* Generated by starterDD\n"
                "* " + "-" * 40 + "\n\n"
            )
            # MODULE
            _wrap_write(f, self.build_module_declaration())
            # PROCEDURE
            _wrap_write(f, self.build_procedure_declaration())
            # LINKED_LIST
            f.write(
                self.build_linked_list_declaration()
            )
            # SEQ_BINARY
            f.write(
                self.build_seq_binary_declaration()
            )
            # Variables
            f.write(
                self.build_variable_declarations()
            )
            f.write("\n")
            # SEQ_ASCII (TDT imports)
            f.write(
                self.build_seq_ascii_declarations()
            )
            f.write("\n")
            # Body lines
            for ln in self.body_lines:
                _wrap_write(f, ln)
            f.write("\n")
            # Footer
            f.write("END: ;\n")
            f.write("QUIT .\n")

        print(f"[x2m] Wrote {filepath}")
        return filepath


# ---------------------------------------------------------------
# Sub-procedure (.c2m)
# ---------------------------------------------------------------
class sub_procedure(CLE2000_procedure):
    def __init__(self, procedure_name,
                 parameters_structures=None,
                 parameters_variables=None):
        super().__init__(
            "sub-procedure", procedure_name,
            parameters_structures, parameters_variables,
        )
        self.main_proc = None

    def set_main_procedure(self, main_proc):
        self.main_proc = main_proc

    def build_parameter_declaration(self):
        out = "PARAMETER"
        for s in self.parameters_structures:
            out += f" {s}"
        out += " ::\n"
        for s in self.parameters_structures:
            out += f"::: LINKED_LIST {s} ;\n"
        out += ";\n"
        return out

    def build_input_variable_recovery(self):
        out = ""
        for name, vtype, _ in self.variables:
            out += f"{vtype} {name} ;\n"
            out += f":: >>{name}<< ;\n"
        return out

