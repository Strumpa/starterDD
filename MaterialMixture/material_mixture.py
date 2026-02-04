# Class to handle abstract material mixtures
# Author : R. Guasch
# Date : 04/02/2026
# Purpose : Define a material mixture with its index and composition
# uses the composition class which is unique per material, 
# Material mixture allows to have different mixtures with same materials but different properties
# -----------------------------------------------------------------------------------------------

class Composition:
    def __init__(self, components):
        self.components = components  # components is a dict of {material_name: fraction}


class MaterialMixture:
    def __init__(self, material_mixture_index):
        self.material_mixture_index = material_mixture_index

    def set_composition(self, composition):
        self.composition = composition
        