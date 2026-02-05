## collection of classes to handle DRAGON module calls
# Author : R. Guasch
# Date : 05/02/2026
# Purpose : Define functions to handle DRAGON module calls for starterDD package
# ----------------------------------------------------------------------------------------------- 

# classes representing DRAGON modules 
import os

class LIB:
    def __init__(self, module_name: str):
        """
        DRAGON Module initialization.
        
        :param module_name (str): Name of the DRAGON module
        """
        self.module_name = module_name

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
