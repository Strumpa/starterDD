# -*- coding: utf-8 -*-
"""Physics constants and energy-mesh utilities shared by the post-treatment.

This module is the single place for:

- the MT <-> DRAGON reaction-key correspondence used to build absorption and
  fission rates (formerly ``REACTIONS_FISS`` / ``REACTIONS_ABS`` in the GE14
  analysis scripts),
- recoverable energy per reaction (Q-values) per nuclear data library,
- named energy meshes (eV, ascending),
- condensation of fine-group tallies onto a coarser mesh, including the
  propagation of statistical errors.

Conventions
-----------
- Energies are in **eV** and meshes are **ascending** (Serpent2 order).
  DRAGON stores ``ENERGY`` descending; use :func:`to_ascending` when needed.
- Group-wise arrays have the energy axis **last**.
- Condensation indices follow the legacy convention: ``idx[k]`` is the
  (exclusive) end, in fine groups, of coarse group ``k``
  (e.g. 2g from SHEM295 -> ``[34, 295]``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Reactions
# ─────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Reaction:
    """A partial reaction that can be scored in both codes.

    Attributes:
        mt: ENDF MT number (Serpent2 ``dr`` response).
        dragon_key: Record name of the microscopic cross section in a DRAGON
            COMPO / MICROLIB isotope directory.
        absorption_factor: Contribution of one reaction to the neutron
            balance as "absorption" (+1 for capture-like reactions, -1 for
            (n,2n), -2 for (n,3n), 0 if not counted as absorption).
        name: Human-readable name.
    """

    mt: int
    dragon_key: str
    absorption_factor: int
    name: str


FISSION = Reaction(18, "NFTOT", 1, "fission")

#: Reactions that make up the absorption rate (fission excluded, see
#: :data:`FISSION`). Order and factors match the legacy ``REACTIONS_ABS``.
ABSORPTION_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(102, "NG", 1, "n,gamma"),
    Reaction(103, "NP", 1, "n,p"),
    Reaction(104, "ND", 1, "n,d"),
    Reaction(105, "NT", 1, "n,t"),
    Reaction(107, "NA", 1, "n,alpha"),
    Reaction(108, "N2A", 1, "n,2alpha"),
    Reaction(28, "NNP", 1, "n,np"),
    Reaction(16, "N2N", -1, "n,2n"),
    Reaction(17, "N3N", -2, "n,3n"),
)

#: Every reaction known to the comparison, keyed by MT.
REACTIONS_BY_MT: Dict[int, Reaction] = {
    r.mt: r for r in (FISSION, *ABSORPTION_REACTIONS)
}
#: Same, keyed by DRAGON record name.
REACTIONS_BY_DRAGON_KEY: Dict[str, Reaction] = {
    r.dragon_key: r for r in REACTIONS_BY_MT.values()
}

#: Legacy tuple views, kept so refactored code can be diffed against the
#: original analysis functions.
REACTIONS_FISS = [(FISSION.mt, FISSION.dragon_key)]
REACTIONS_ABS = [(r.mt, r.dragon_key, r.absorption_factor)
                 for r in ABSORPTION_REACTIONS]

#: Serpent2 MT numbers that are *not* partial reactions but appear in legacy
#: detector inputs (``dr 101`` = disappearance). They are kept in manifests
#: but ignored when assembling absorption from partial reactions.
SERPENT_SUMMARY_MTS = {101: "disappearance", 27: "absorption", 1: "total",
                       452: "nubar"}

#: MT numbers scored per isotope by the Serpent detectors that starterDD
#: writes, per library: only the partial reactions the evaluation provides
#: (as in the GE14 references; e.g. U235 has no MT 103 in ENDF/B-VIII.1, U234
#: only 16, 17, 18, 37, 102 among the partial reactions). A response on a
#: reaction absent from the ACE data is not valid in Serpent, so an isotope
#: without a list here (or in ``reaction_mts=`` overrides) raises
#: :class:`UnknownReactionsError` when its detectors are written.
SERPENT_DETECTOR_MTS: Dict[str, Dict[str, Tuple[int, ...]]] = {
    "endfb8r1": {
        "U234": (18, 102, 16, 17),
        "U235": (18, 102, 107, 16, 17),
        "U238": (18, 102, 103, 107, 16, 17),
        "Gd155": (102,),
        "Gd157": (102,),
    },
    "jeff311": {
        "U235": (18, 102, 16, 17, 37),
        "U238": (18, 102, 16, 17, 37),
        "Gd155": (102,),
        "Gd157": (102,),
    },
}
DEFAULT_DETECTOR_MTS: Tuple[int, ...] = (FISSION.mt,) + tuple(r.mt for r in ABSORPTION_REACTIONS)


@dataclass(frozen=True)
class MacroReaction:
    """A macroscopic (all-nuclide) rate scored by a Serpent ``dr <mt> void``
    response and rebuilt in DRAGON from the edited isotopes plus the
    ``*MAC*RES`` residual pseudo-isotope (density 1, macroscopic data).

    Serpent ``-2`` (and ``dr 101``) sums the partial reactions with no
    outgoing neutron (ACE ``TY == 0``): capture only,
    without fission (``TY = 19``) or (n,xn). Checked on the GE14 references:
    ``101`` equals (n,γ)+(n,p)+(n,α) to < 4e-7.
    """

    name: str
    serpent_mts: Tuple[int, ...]       # summed void responses
    dragon_keys: Tuple[str, ...]       # summed DRAGON records


#: Macroscopic quantities: ``absorption`` = capture + fission, the usual
#: definition (and that of :func:`compute_pin_rates`, without (n,xn)).
MACRO_REACTIONS: Dict[str, MacroReaction] = {
    "fission": MacroReaction("fission", (-6,), ("NFTOT",)),
    "capture": MacroReaction("capture", (-2,), ("NG", "NP", "ND", "NT", "NA", "N2A")),
    "absorption": MacroReaction("absorption", (-2, -6),
                                ("NG", "NP", "ND", "NT", "NA", "N2A", "NFTOT")),
    "nufission": MacroReaction("nufission", (-7,), ("NUSIGF",)),
    "total": MacroReaction("total", (-1,), ("NTOT0",)),
}
#: Material name of Serpent macroscopic (collision-material) responses.
VOID = "void"


class UnknownReactionsError(KeyError):
    """No list of available reactions for an isotope (see
    :data:`SERPENT_DETECTOR_MTS`)."""


def serpent_detector_mts(library: str, isotope: str,
                         overrides: Optional[Mapping[str, Sequence[int]]] = None
                         ) -> Tuple[int, ...]:
    """MT numbers of the Serpent responses written for *isotope*: from
    *overrides* (``{isotope: [MT, ...]}``), else
    :data:`SERPENT_DETECTOR_MTS`. Raises :class:`UnknownReactionsError` if
    neither has the isotope: writing reactions the ACE data lacks is invalid,
    and a guess would silently change the compared absorption."""
    if overrides and isotope in overrides:
        return tuple(int(mt) for mt in overrides[isotope])
    table = SERPENT_DETECTOR_MTS.get(canonical_library(library), {})
    if isotope not in table:
        raise UnknownReactionsError(unknown_reactions_message(library, [isotope]))
    return table[isotope]


def unknown_reactions_message(library: str, isotopes: Sequence[str]) -> str:
    return (f"No Serpent reaction list for {list(isotopes)} ({canonical_library(library)}). "
            "Give the partial reactions their ACE data has (Serpent prints them per "
            f"nuclide), among {list(DEFAULT_DETECTOR_MTS)}: reaction_mts={{'"
            f"{isotopes[0]}': [18, 102, 16, 17], ...}} (study YAML 'serpent: reaction_mts:'), "
            "or add them to physics_data.SERPENT_DETECTOR_MTS.")

# ─────────────────────────────────────────────────────────────────────────────
# Q-values (MeV), per library
# ─────────────────────────────────────────────────────────────────────────────

#: Recoverable energy per reaction, MeV. ``Q_VALUES[library][isotope][mt]``.
#: ENDF/B-VIII.1 values are those used by the original GE14 analysis.
Q_VALUES: Dict[str, Dict[str, Dict[int, float]]] = {
    "endfb8r1": {
        "U235": {18: 194.077, 102: 6.54443, 107: 0.0, 16: -5.298, 17: -12.141},
        "U238": {18: 197.786, 102: 4.80643, 103: 0.0, 107: 0.0, 16: -6.153,
                 17: -11.279},
        "Gd155": {102: 8.536},
        "Gd157": {102: 7.937},
    },
}

#: Aliases so that DRAGON library tags and Serpent library names resolve to
#: the same Q-value table.
LIBRARY_ALIASES = {
    "endfb8r1": "endfb8r1",
    "endfb8r1295": "endfb8r1",
    "endfb81": "endfb8r1",
    "endf/b-viii.1": "endfb8r1",
}


def canonical_library(library: str) -> str:
    """Return the canonical key of *library* (e.g. ``"endfb8r1295"`` ->
    ``"endfb8r1"``)."""
    key = library.lower().replace("_", "").replace("-", "")
    for alias, canon in LIBRARY_ALIASES.items():
        if key == alias.replace("_", "").replace("-", ""):
            return canon
    return library.lower()


def q_value(library: str, isotope: str, mt: int) -> float:
    """Recoverable energy (MeV) of reaction *mt* on *isotope*.

    Returns 0.0 when the reaction is not tabulated (same behaviour as the
    legacy ``Q_VALUES.get(iso, 0).get(mt, 0)``), but raises ``KeyError`` for an
    unknown library so that a missing table is never silently zero.
    """
    lib = canonical_library(library)
    if lib not in Q_VALUES:
        raise KeyError(
            f"No Q-value table for library '{library}'. "
            f"Known: {sorted(Q_VALUES)}"
        )
    return Q_VALUES[lib].get(isotope, {}).get(mt, 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# Nuclide naming
# ─────────────────────────────────────────────────────────────────────────────

_ELEMENTS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te "
    "I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir "
    "Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm"
).split()
Z_TO_SYMBOL = {z: s for z, s in enumerate(_ELEMENTS)}
SYMBOL_TO_Z = {s: z for z, s in Z_TO_SYMBOL.items()}

_ZAID_RE = re.compile(r"^(\d{1,3})(\d{3})(?:\.\w+)?$")
_NAME_RE = re.compile(r"^([A-Z][a-z]?)-?(\d{1,3})(m\d?)?$")


def nuclide_name(token: str) -> Optional[str]:
    """Normalise a nuclide token to the DRAGON-style name (``U235``).

    Accepts Serpent ZAIDs with or without library suffix (``92235.09c``,
    ``92235``), DRAGON names (``U235``, ``Gd155``) and dashed names
    (``U-235``). Returns ``None`` if *token* is not recognisable as a nuclide
    (e.g. thermal-scattering libraries such as ``lwj3.11t`` or elements).
    """
    token = token.strip()
    m = _ZAID_RE.match(token)
    if m:
        z, a = int(m.group(1)), int(m.group(2))
        if a == 0 or z not in Z_TO_SYMBOL:
            return None
        if a >= 300:
            # Isomeric states use library-dependent A offsets (NJOY/Serpent
            # vs MCNP conventions); do not guess. Callers fall back to the
            # material name.
            return None
        return f"{Z_TO_SYMBOL[z]}{a}"
    m = _NAME_RE.match(token)
    if m and m.group(1) in SYMBOL_TO_Z:
        return f"{m.group(1)}{int(m.group(2))}{'m' if m.group(3) else ''}"
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Energy meshes (eV, ascending)
# ─────────────────────────────────────────────────────────────────────────────

_MEV = 1.0e6

# All presets are taken from the DRAGON5 library mesh: DRAGON is the code of
# interest, so its bounds define the code-to-code comparison meshes and the
# Serpent2 ``ene`` cards written by the exporter.

# 2 groups, thermal cut at 0.625 eV (SHEM295 bounds 0, 34 and 295).
G2_EV = np.array([1.10002700E-10, 6.24998689E-07, 1.96403000E+01]) * _MEV

# 26 groups condensed from SHEM295 (DRAGON 2-level schemes, 2nd level).
G26_EV = np.array([
    1.10002700E-10, 1.04505001E-08, 3.43997590E-08, 7.64968619E-08,
    1.37999400E-07, 2.31192306E-07, 3.52993488E-07, 6.24998689E-07,
    8.80024374E-07, 9.63959813E-07, 1.00903499E-06, 1.10395002E-06,
    1.14796901E-06, 1.25093901E-06, 4.00000000E-06, 5.29895287E-05,
    4.19093597E-04, 1.81183301E-03, 9.11880762E-03, 2.49990801E-02,
    6.73793828E-02, 1.95066203E-01, 4.94001812E-01, 1.33694100E+00,
    2.23129900E+00, 4.96584700E+00, 1.96403000E+01,
]) * _MEV

#: SHEM295 exactly as stored in the DRAGON5 library ``draglibendfb8r1SHEM295_v5p1``
#: (``ENERGY`` record of a 1-level COMPO, full precision). Canonical ``295g``.
G295_EV = np.array([
    1.10002700E-10, 2.49989703E-09, 4.55602119E-09, 7.14526279E-09, 1.04505001E-08, 1.48299597E-08,
    2.00103503E-08, 2.49394197E-08, 2.92988904E-08, 3.43997590E-08, 4.02999297E-08, 4.73018587E-08,
    5.54981492E-08, 6.51993603E-08, 7.64968619E-08, 8.97968337E-08, 1.04297698E-07, 1.19994901E-07,
    1.37999400E-07, 1.61895305E-07, 1.90004900E-07, 2.09610194E-07, 2.31192306E-07, 2.54996508E-07,
    2.79988796E-07, 3.05011511E-07, 3.25007886E-07, 3.52993488E-07, 3.90001088E-07, 4.31578606E-07,
    4.75016505E-07, 5.20010829E-07, 5.54989696E-07, 5.94992995E-07, 6.24998689E-07, 7.19998896E-07,
    8.20037127E-07, 8.80024374E-07, 9.19977903E-07, 9.44022179E-07, 9.63959813E-07, 9.81959105E-07,
    9.96500492E-07, 1.00903499E-06, 1.02101195E-06, 1.03499305E-06, 1.07798600E-06, 1.09198201E-06,
    1.10395002E-06, 1.11604905E-06, 1.12997401E-06, 1.14796901E-06, 1.16998899E-06, 1.21396804E-06,
    1.25093901E-06, 1.29303801E-06, 1.33095205E-06, 1.38098097E-06, 1.41000700E-06, 1.44396698E-06,
    1.51997602E-06, 1.58802998E-06, 1.66894901E-06, 1.77996600E-06, 1.90007699E-06, 1.98992002E-06,
    2.07009506E-06, 2.15694809E-06, 2.21708703E-06, 2.27298594E-06, 2.33006096E-06, 2.46994090E-06,
    2.55000305E-06, 2.59009409E-06, 2.62005305E-06, 2.64004111E-06, 2.70011497E-06, 2.71989799E-06,
    2.74092197E-06, 2.77512097E-06, 2.88404703E-06, 3.14210892E-06, 3.54307294E-06, 3.71208692E-06,
    3.88216996E-06, 4.00000000E-06, 4.21982813E-06, 4.30981207E-06, 4.40215588E-06, 4.63248920E-06,
    4.95845604E-06, 5.31798410E-06, 5.61866713E-06, 5.91856718E-06, 6.22824383E-06, 6.50840521E-06,
    6.82842684E-06, 7.21451092E-06, 7.62242317E-06, 7.96529818E-06, 8.41566086E-06, 8.85599327E-06,
    9.31004906E-06, 9.78738499E-06, 1.02378798E-05, 1.07091103E-05, 1.12020197E-05, 1.17410898E-05,
    1.22324800E-05, 1.28853903E-05, 1.35731602E-05, 1.42121096E-05, 1.49706898E-05, 1.57382498E-05,
    1.65617294E-05, 1.73760796E-05, 1.83035393E-05, 1.90314808E-05, 2.01075306E-05, 2.05342503E-05,
    2.12231903E-05, 2.20011406E-05, 2.26711807E-05, 2.46578293E-05, 2.78851509E-05, 3.16929493E-05,
    3.30854683E-05, 3.45391808E-05, 3.56979904E-05, 3.60567589E-05, 3.64191399E-05, 3.68587990E-05,
    3.73037682E-05, 3.77918816E-05, 3.87873611E-05, 3.97295113E-05, 4.12270393E-05, 4.21440887E-05,
    4.31246300E-05, 4.41721382E-05, 4.52903709E-05, 4.62052917E-05, 4.75173187E-05, 4.92591095E-05,
    5.17846794E-05, 5.29895287E-05, 5.40599899E-05, 5.70594902E-05, 5.99250298E-05, 6.23082809E-05,
    6.36305885E-05, 6.45922470E-05, 6.50459824E-05, 6.55028992E-05, 6.58312302E-05, 6.61612091E-05,
    6.64928513E-05, 6.68261414E-05, 6.90681992E-05, 7.18869171E-05, 7.35594788E-05, 7.63321609E-05,
    7.93679276E-05, 8.39393387E-05, 8.87740479E-05, 9.33255920E-05, 9.73287430E-05, 1.00594200E-04,
    1.01098396E-04, 1.01605202E-04, 1.02114502E-04, 1.03037598E-04, 1.05646103E-04, 1.10287903E-04,
    1.12853897E-04, 1.15479698E-04, 1.16523697E-04, 1.17577103E-04, 1.20553596E-04, 1.26228600E-04,
    1.32700500E-04, 1.39504196E-04, 1.46656693E-04, 1.54175903E-04, 1.63056107E-04, 1.67518600E-04,
    1.75229095E-04, 1.83294495E-04, 1.84951599E-04, 1.86250793E-04, 1.87559204E-04, 1.88876694E-04,
    1.90203506E-04, 1.93078003E-04, 1.95996002E-04, 2.00957703E-04, 2.12107697E-04, 2.24324707E-04,
    2.35590302E-04, 2.41796005E-04, 2.56747803E-04, 2.68296906E-04, 2.76467804E-04, 2.84887512E-04,
    2.88326691E-04, 2.95921509E-04, 3.19927490E-04, 3.35322998E-04, 3.53574585E-04, 3.71702698E-04,
    3.90760315E-04, 4.19093597E-04, 4.53998688E-04, 5.01746185E-04, 5.39204224E-04, 5.77145508E-04,
    5.92940674E-04, 6.00098816E-04, 6.12834229E-04, 6.46836975E-04, 6.77286499E-04, 7.48517273E-04,
    8.32217896E-04, 9.09681274E-04, 9.82494080E-04, 1.06432300E-03, 1.13466699E-03, 1.34358203E-03,
    1.58619702E-03, 1.81183301E-03, 2.08410400E-03, 2.39729004E-03, 2.70023608E-03, 2.99618311E-03,
    3.48106812E-03, 3.93668506E-03, 4.85660205E-03, 5.99148389E-03, 7.39156201E-03, 9.11880762E-03,
    1.11377402E-02, 1.36036602E-02, 1.48996699E-02, 1.62004502E-02, 1.85847109E-02, 2.26994102E-02,
    2.49990801E-02, 2.61000996E-02, 2.73944102E-02, 2.92810098E-02, 3.34596094E-02, 3.69785898E-02,
    4.08676602E-02, 4.99158711E-02, 5.51655703E-02, 6.73793828E-02, 8.22973594E-02, 9.46645000E-02,
    1.15623500E-01, 1.22773203E-01, 1.40097594E-01, 1.65065000E-01, 1.95066203E-01, 2.30059797E-01,
    2.67826406E-01, 3.20646406E-01, 3.83883500E-01, 4.12501188E-01, 4.56021094E-01, 4.94001812E-01,
    5.78442500E-01, 7.06511188E-01, 8.60005812E-01, 9.51118875E-01, 1.05114900E+00, 1.16204800E+00,
    1.28696100E+00, 1.33694100E+00, 1.40576800E+00, 1.63653900E+00, 1.90138700E+00, 2.23129900E+00,
    2.72531400E+00, 3.32870700E+00, 4.06569100E+00, 4.96584700E+00, 6.06529900E+00, 6.70319200E+00,
    7.40817300E+00, 8.18729700E+00, 9.04836300E+00, 9.99998700E+00, 1.16183300E+01, 1.38402900E+01,
    1.49182300E+01, 1.96403000E+01,
]) * _MEV

#: Deprecated alias of :data:`G295_EV` (P1 name, kept for compatibility).
G295_DRAGON_EV = G295_EV

#: The ``295g`` grid written by the starterDD Serpent2 exporter **before
#: 2026-11-14**, used by every existing GE14 reference (NSE paper runs). It
#: differs from :data:`G295_EV` in four bounds (indices 235-238, 3.9-7.5 keV,
#: up to 4 %): it continues the Delta-u = 0.2 ladder below 9.12 keV where the
#: DRAGON library has four Delta-u = 0.21 groups. The outer bounds 3481.07 and
#: 9118.81 eV are shared, so ascending groups 234-238 (DRAGON groups 57-61)
#: merge exactly (``mesh="common"``, 291 groups). Every other bound agrees
#: within 5e-6. Only used to recognise legacy references.
G295_STARTERDD_LEGACY_EV = np.array([
    1.10003E-10, 2.49990E-09, 4.55602E-09, 7.14526E-09, 1.04505E-08, 1.48300E-08, 2.00104E-08, 2.49394E-08, 2.92989E-08, 3.43998E-08,
    4.02999E-08, 4.73019E-08, 5.54981E-08, 6.51994E-08, 7.64969E-08, 8.97968E-08, 1.04298E-07, 1.19995E-07, 1.37999E-07, 1.61895E-07,
    1.90005E-07, 2.09610E-07, 2.31192E-07, 2.54997E-07, 2.79989E-07, 3.05012E-07, 3.25008E-07, 3.52993E-07, 3.90001E-07, 4.31579E-07,
    4.75017E-07, 5.20011E-07, 5.54990E-07, 5.94993E-07, 6.24999E-07, 7.19999E-07, 8.20037E-07, 8.80024E-07, 9.19978E-07, 9.44022E-07,
    9.63960E-07, 9.81959E-07, 9.96500E-07, 1.00903E-06, 1.02101E-06, 1.03499E-06, 1.07799E-06, 1.09198E-06, 1.10395E-06, 1.11605E-06,
    1.12997E-06, 1.14797E-06, 1.16999E-06, 1.21397E-06, 1.25094E-06, 1.29304E-06, 1.33095E-06, 1.38098E-06, 1.41001E-06, 1.44397E-06,
    1.51998E-06, 1.58803E-06, 1.66895E-06, 1.77997E-06, 1.90008E-06, 1.98992E-06, 2.07010E-06, 2.15695E-06, 2.21709E-06, 2.27299E-06,
    2.33006E-06, 2.46994E-06, 2.55000E-06, 2.59009E-06, 2.62005E-06, 2.64004E-06, 2.70011E-06, 2.71990E-06, 2.74092E-06, 2.77512E-06,
    2.88405E-06, 3.14211E-06, 3.54307E-06, 3.71209E-06, 3.88217E-06, 4.00000E-06, 4.21983E-06, 4.30981E-06, 4.40216E-06, 4.63249E-06,
    4.95846E-06, 5.31798E-06, 5.61867E-06, 5.91857E-06, 6.22824E-06, 6.50841E-06, 6.82843E-06, 7.21451E-06, 7.62242E-06, 7.96530E-06,
    8.41566E-06, 8.85599E-06, 9.31005E-06, 9.78738E-06, 1.02379E-05, 1.07091E-05, 1.12020E-05, 1.17411E-05, 1.22325E-05, 1.28854E-05,
    1.35732E-05, 1.42121E-05, 1.49707E-05, 1.57382E-05, 1.65617E-05, 1.73761E-05, 1.83035E-05, 1.90315E-05, 2.01075E-05, 2.05343E-05,
    2.12232E-05, 2.20011E-05, 2.26712E-05, 2.46578E-05, 2.78852E-05, 3.16929E-05, 3.30855E-05, 3.45392E-05, 3.56980E-05, 3.60568E-05,
    3.64191E-05, 3.68588E-05, 3.73038E-05, 3.77919E-05, 3.87874E-05, 3.97295E-05, 4.12270E-05, 4.21441E-05, 4.31246E-05, 4.41721E-05,
    4.52904E-05, 4.62053E-05, 4.75173E-05, 4.92591E-05, 5.17847E-05, 5.29895E-05, 5.40600E-05, 5.70595E-05, 5.99250E-05, 6.23083E-05,
    6.36306E-05, 6.45922E-05, 6.50460E-05, 6.55029E-05, 6.58312E-05, 6.61612E-05, 6.64929E-05, 6.68261E-05, 6.90682E-05, 7.18869E-05,
    7.35595E-05, 7.63322E-05, 7.93679E-05, 8.39393E-05, 8.87740E-05, 9.33256E-05, 9.73287E-05, 1.00594E-04, 1.01098E-04, 1.01605E-04,
    1.02115E-04, 1.03038E-04, 1.05646E-04, 1.10288E-04, 1.12854E-04, 1.15480E-04, 1.16524E-04, 1.17577E-04, 1.20554E-04, 1.26229E-04,
    1.32701E-04, 1.39504E-04, 1.46657E-04, 1.54176E-04, 1.63056E-04, 1.67519E-04, 1.75229E-04, 1.83294E-04, 1.84952E-04, 1.86251E-04,
    1.87559E-04, 1.88877E-04, 1.90204E-04, 1.93078E-04, 1.95996E-04, 2.00958E-04, 2.12108E-04, 2.24325E-04, 2.35590E-04, 2.41796E-04,
    2.56748E-04, 2.68297E-04, 2.76468E-04, 2.84888E-04, 2.88327E-04, 2.95922E-04, 3.19927E-04, 3.35323E-04, 3.53575E-04, 3.71703E-04,
    3.90760E-04, 4.19094E-04, 4.53999E-04, 5.01746E-04, 5.39204E-04, 5.77146E-04, 5.92941E-04, 6.00099E-04, 6.12834E-04, 6.46837E-04,
    6.77286E-04, 7.48517E-04, 8.32218E-04, 9.09681E-04, 9.82494E-04, 1.06432E-03, 1.13467E-03, 1.34358E-03, 1.58620E-03, 1.81183E-03,
    2.08410E-03, 2.39729E-03, 2.70024E-03, 2.99618E-03, 3.48107E-03, 4.09735E-03, 5.00451E-03, 6.11252E-03, 7.46585E-03, 9.11881E-03,
    1.11377E-02, 1.36037E-02, 1.48997E-02, 1.62005E-02, 1.85847E-02, 2.26994E-02, 2.49991E-02, 2.61001E-02, 2.73944E-02, 2.92810E-02,
    3.34596E-02, 3.69786E-02, 4.08677E-02, 4.99159E-02, 5.51656E-02, 6.73794E-02, 8.22974E-02, 9.46645E-02, 1.15623E-01, 1.22773E-01,
    1.40098E-01, 1.65065E-01, 1.95066E-01, 2.30060E-01, 2.67826E-01, 3.20646E-01, 3.83883E-01, 4.12501E-01, 4.56021E-01, 4.94002E-01,
    5.78442E-01, 7.06511E-01, 8.60006E-01, 9.51119E-01, 1.05115E+00, 1.16205E+00, 1.28696E+00, 1.33694E+00, 1.40577E+00, 1.63654E+00,
    1.90139E+00, 2.23130E+00, 2.72531E+00, 3.32871E+00, 4.06569E+00, 4.96585E+00, 6.06530E+00, 6.70319E+00, 7.40817E+00, 8.18730E+00,
    9.04836E+00, 9.99999E+00, 1.16183E+01, 1.38403E+01, 1.49182E+01, 1.96403E+01,
]) * _MEV

#: Named meshes. Keys are the lower-case Serpent grid names of the naming
#: convention (§3.2 of the plan).
ENERGY_MESHES: Dict[str, np.ndarray] = {
    "2g": G2_EV,
    "26g": G26_EV,
    "295g": G295_EV,
}

#: Grids written by earlier starterDD versions that differ from the presets.
#: Recognised by :func:`mesh_code`, never written.
LEGACY_MESHES: Dict[str, np.ndarray] = {
    "295g": G295_STARTERDD_LEGACY_EV,
}

#: Default relative tolerance when matching energy bounds between codes.
#: Legacy ``ene`` cards and det0.m bounds have 5-6 significant digits.
BOUNDS_RTOL = 1.0e-4


def mesh(name: str) -> np.ndarray:
    """Return a copy of the named mesh (eV, ascending)."""
    try:
        return ENERGY_MESHES[name.lower()].copy()
    except KeyError:
        raise KeyError(
            f"Unknown energy mesh '{name}'. Known: {sorted(ENERGY_MESHES)}"
        ) from None


def mesh_code(bounds: Sequence[float]) -> str:
    """Naming-convention mesh code for a set of bounds: ``"2G"``, ``"295G"``,
    or ``"C{n}G"`` for a custom mesh that is not one of the presets.

    The legacy exporter grid (:data:`G295_STARTERDD_LEGACY_EV`) is also
    ``"295G"``, so manifests of existing references keep their names."""
    b = to_ascending(bounds)
    n = len(b) - 1
    for name, ref in ENERGY_MESHES.items():
        if len(ref) == len(b) and bounds_match(ref, b):
            return name.upper()
    for name, ref in LEGACY_MESHES.items():
        if len(ref) == len(b) and bounds_match(ref, b):
            return name.upper()
    return f"C{n}G"


def is_legacy_mesh(bounds: Sequence[float]) -> bool:
    """True if *bounds* is one of :data:`LEGACY_MESHES` (e.g. the pre-fix
    Serpent ``295g`` grid)."""
    b = to_ascending(bounds)
    return any(len(ref) == len(b) and bounds_match(ref, b)
               for ref in LEGACY_MESHES.values())


def to_ascending(bounds: Sequence[float]) -> np.ndarray:
    """Return *bounds* as an ascending float array."""
    b = np.asarray(bounds, dtype=float).ravel()
    if b.size > 1 and b[0] > b[-1]:
        b = b[::-1]
    return b


def bounds_match(a: Sequence[float], b: Sequence[float],
                 rtol: float = BOUNDS_RTOL) -> bool:
    """True if two meshes have the same number of bounds and agree within
    *rtol* (bounds compared after sorting ascending)."""
    a, b = to_ascending(a), to_ascending(b)
    return a.shape == b.shape and bool(np.allclose(a, b, rtol=rtol, atol=0.0))


class MeshMismatchError(ValueError):
    """Raised when a coarse mesh is not a sub-mesh of a fine mesh."""


def condensation_indices(fine: Sequence[float], coarse: Sequence[float],
                         rtol: float = BOUNDS_RTOL) -> np.ndarray:
    """Indices that condense *fine* groups onto *coarse* groups.

    Both meshes are given as bounds (any order, eV). Every coarse bound must
    coincide (within *rtol*) with a fine bound and the outer bounds must
    match, i.e. *fine* must be a refinement of *coarse*.

    Returns:
        ``idx`` of length ``n_coarse`` with ``idx[k]`` the exclusive end, in
        fine groups (ascending energy), of coarse group ``k``. Example: SHEM295
        -> 2g gives ``[34, 295]``.

    Raises:
        MeshMismatchError: listing the coarse bounds without a fine match.
    """
    f, c = to_ascending(fine), to_ascending(coarse)
    pos = []
    missing = []
    for e in c:
        j = int(np.argmin(np.abs(f - e)))
        if abs(f[j] - e) <= rtol * abs(e):
            pos.append(j)
        else:
            missing.append(e)
    if missing:
        raise MeshMismatchError(
            f"The {len(f) - 1}-group mesh is not a refinement of the "
            f"{len(c) - 1}-group mesh: no fine bound within rtol={rtol:g} of "
            f"coarse bound(s) {', '.join(f'{e:.6e}' for e in missing)} eV."
        )
    if pos[0] != 0 or pos[-1] != len(f) - 1:
        raise MeshMismatchError(
            "Outer bounds differ: fine mesh spans "
            f"[{f[0]:.6e}, {f[-1]:.6e}] eV, coarse mesh spans "
            f"[{c[0]:.6e}, {c[-1]:.6e}] eV."
        )
    if np.any(np.diff(pos) <= 0):
        raise MeshMismatchError("Coarse bounds map to non-increasing fine "
                                "bounds (duplicate or unsorted mesh).")
    return np.asarray(pos[1:], dtype=int)


def common_mesh(a: Sequence[float], b: Sequence[float],
                rtol: float = BOUNDS_RTOL) -> np.ndarray:
    """Coarsest mesh whose bounds belong to both *a* and *b* (ascending, values
    taken from *a*). Both meshes can be condensed onto it exactly. The outer
    bounds must match."""
    a, b = to_ascending(a), to_ascending(b)
    keep = [e for e in a if np.min(np.abs(b - e)) <= rtol * abs(e)]
    if not keep or abs(keep[0] - a[0]) > 0 or abs(keep[-1] - a[-1]) > 0:
        raise MeshMismatchError("Meshes do not share their outer bounds.")
    return np.asarray(keep)


def merged_groups(fine: Sequence[float], coarse: Sequence[float]) -> List[Tuple[int, int]]:
    """``[(first, last), ...]`` 0-based ranges of *fine* groups merged into a
    single *coarse* group (ranges of length 1 are omitted)."""
    idx = condensation_indices(fine, coarse)
    starts = np.concatenate(([0], idx[:-1]))
    return [(int(s), int(e) - 1) for s, e in zip(starts, idx) if e - s > 1]


def condense(values: np.ndarray, indices: Sequence[int],
             axis: int = -1) -> np.ndarray:
    """Sum fine-group *values* into coarse groups along *axis*.

    *indices* as returned by :func:`condensation_indices`.
    """
    values = np.asarray(values, dtype=float)
    starts = np.concatenate(([0], np.asarray(indices)[:-1]))
    return np.add.reduceat(values, starts, axis=axis)


def condense_abs_error(abs_err: np.ndarray, indices: Sequence[int],
                       mode: str = "quadrature", axis: int = -1) -> np.ndarray:
    """Condense **absolute** 1-sigma errors.

    Args:
        abs_err: Absolute errors (convert Serpent relative errors with
            ``rel * value`` first).
        mode: ``"quadrature"`` (sigma^2 = sum sigma_i^2, statistically
            correct for independent bins) or ``"linear"`` (sum sigma_i, the
            legacy behaviour kept for the P1 regression).
    """
    abs_err = np.asarray(abs_err, dtype=float)
    if mode == "quadrature":
        return np.sqrt(condense(abs_err ** 2, indices, axis=axis))
    if mode == "linear":
        return condense(abs_err, indices, axis=axis)
    raise ValueError(f"Unknown error condensation mode '{mode}'.")


def lethargy_widths(bounds: Sequence[float]) -> np.ndarray:
    """ln(E_hi / E_lo) for each group of an ascending mesh."""
    b = to_ascending(bounds)
    return np.log(b[1:] / b[:-1])


def equal_lethargy_mesh(n_groups: int, e_min: float, e_max: float) -> np.ndarray:
    """Equal-lethargy mesh (Serpent ``ene`` type 2), ascending."""
    return np.geomspace(e_min, e_max, n_groups + 1)


def uniform_mesh(n_groups: int, e_min: float, e_max: float) -> np.ndarray:
    """Equal-width mesh (linear spacing in energy), ascending."""
    return np.linspace(e_min, e_max, n_groups + 1)


__all__ = [
    "Reaction", "FISSION", "ABSORPTION_REACTIONS", "REACTIONS_BY_MT",
    "REACTIONS_BY_DRAGON_KEY", "REACTIONS_FISS", "REACTIONS_ABS",
    "SERPENT_SUMMARY_MTS", "Q_VALUES", "canonical_library", "q_value",
    "nuclide_name", "Z_TO_SYMBOL", "SYMBOL_TO_Z",
    "G2_EV", "G26_EV", "G295_EV", "G295_DRAGON_EV", "G295_STARTERDD_LEGACY_EV",
    "ENERGY_MESHES", "LEGACY_MESHES", "is_legacy_mesh",
    "BOUNDS_RTOL", "mesh", "mesh_code", "to_ascending", "bounds_match",
    "MeshMismatchError", "condensation_indices", "common_mesh", "merged_groups",
    "condense", "condense_abs_error",
    "lethargy_widths", "equal_lethargy_mesh", "uniform_mesh",
]
