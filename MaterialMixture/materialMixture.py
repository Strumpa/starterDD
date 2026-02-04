# Class to handle abstract material mixtures
class MaterialMixture:
    def __init__(self, material_mixture_index):
        self.material_mixture_index = material_mixture_index

    def set_composition(self, composition):
        self.composition = composition
        