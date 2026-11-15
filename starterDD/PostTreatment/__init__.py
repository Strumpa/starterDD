"""PostTreatment module - Helpers to assist in the post-treatment of Dragon Cases.

- :mod:`.pylcm`: LCM objects (ASCII parser, path access, HDF5 cache).
- :mod:`.compo_reader`: statepoint-aware MULTICOMPO reader.
- :mod:`.serpent_reader`: Serpent2 det0/res/input readers and detector
  manifests (including manifests inferred from legacy inputs).
- :mod:`.physics_data`: reactions, Q-values, energy meshes, condensation.
"""

from .pylcm import Pylcm, LcmChars, LcmH5
from .compo_reader import CompoReader, Statepoint, StatepointError, open_compo
from .serpent_reader import (
    DetectorManifest,
    ManifestEntry,
    PinResolver,
    SerpentRun,
    Tally,
    parse_serpent_input,
    read_det,
    read_res,
)
from . import physics_data
from .lattice_maps import PinLattice
from .comparison import (
    KeffComparison, PinRateComparison, SpectrumComparison,
    compare_keff, compare_spectra, compute_pin_rates,
)

# Figures: ``from starterDD.PostTreatment import plotting`` (needs matplotlib).

__all__ = [
    "PinLattice", "KeffComparison", "PinRateComparison", "SpectrumComparison",
    "compare_keff", "compare_spectra", "compute_pin_rates",
    "Pylcm", "LcmChars", "LcmH5",
    "CompoReader", "Statepoint", "StatepointError", "open_compo",
    "DetectorManifest", "ManifestEntry", "PinResolver", "SerpentRun", "Tally",
    "parse_serpent_input", "read_det", "read_res",
    "physics_data",
]
