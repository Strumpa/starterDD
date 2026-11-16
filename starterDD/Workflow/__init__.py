"""Workflow module - studies that chain starterDD cases.

- :mod:`.study`: :class:`CodeToCodeStudy`, DRAGON5 cases over assemblies ×
  scheme variants × operating conditions, their Serpent2 references, the
  study manifest and the code-to-code post-processing (plan §3.9, P5).
"""

from .study import (
    CodeToCodeStudy,
    StudyAssembly,
    StudyCase,
    StudyCondition,
    StudyStatepoint,
    StudyVariant,
)

__all__ = [
    "CodeToCodeStudy", "StudyAssembly", "StudyCase", "StudyCondition",
    "StudyStatepoint", "StudyVariant",
]
