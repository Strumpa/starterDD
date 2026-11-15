"""P0 post-treatment readers: physics_data, Pylcm, CompoReader, serpent_reader.

Test data tiers (see ``tests/data/README.md``):

1. synthetic fixtures written in DRAGON / Serpent2 formats by this module;
2. the committed real GE14 subset ``tests/data/ge14_mini`` (always run);
3. full local data, never committed, marked ``fulldata`` and skipped when
   absent: GE14 samples (``$STARTERDD_TEST_DATA``, default ``temp_examples/``)
   and DRAGON5/Donjon test COMPOs (``$DONJON_DATA_DIR``).
"""

from __future__ import annotations

import itertools
import json
import os
import random
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

from starterDD.PostTreatment import physics_data as phys
from starterDD.PostTreatment.pylcm import LcmChars, LcmH5, Pylcm, h5py_available
from starterDD.PostTreatment.compo_reader import (
    CompoReader, StatepointError, cache_is_valid, decode_tree,
    default_cache_path, open_compo,
)
from starterDD.PostTreatment.serpent_reader import (
    DetectorManifest, ManifestConflictError, PinResolver, SerpentRun,
    default_manifest_path, parse_serpent_input, read_det, read_res,
    serpent_simple_det_reader, serpent_simple_res_reader,
    serpent_simple_serp_reader,
)

from conftest import (
    DONJON_DATA_DIR, FULL_DATA_DIR, GE14_MATERIAL_COMPOS_00_YAML, GE14_MINI_DIR,
)

REPO = Path(__file__).resolve().parents[1]
needs_h5 = pytest.mark.skipif(not h5py_available(), reason="h5py not installed")


# ═════════════════════════════════════════════════════════════════════════════
# physics_data
# ═════════════════════════════════════════════════════════════════════════════

class TestPhysicsData:
    def test_condensation_indices_reproduce_legacy_constants(self):
        legacy_26 = 295 - np.flip(np.array(
            [0, 10, 14, 18, 26, 33, 40, 49, 56, 66, 84, 150, 210, 241, 244,
             247, 252, 255, 258, 261, 268, 273, 277, 281, 286, 291]))
        assert np.array_equal(
            phys.condensation_indices(phys.G295_EV, phys.G26_EV), legacy_26)
        assert np.array_equal(
            phys.condensation_indices(phys.G295_EV, phys.G2_EV), [34, 295])

    def test_condensation_accepts_descending_dragon_bounds(self):
        idx = phys.condensation_indices(phys.G295_EV[::-1], phys.G2_EV[::-1])
        assert np.array_equal(idx, [34, 295])

    def test_non_refinement_raises_with_bounds(self):
        bad = phys.G2_EV.copy()
        bad[1] = 1.0  # 1 eV is not a SHEM295 bound
        with pytest.raises(phys.MeshMismatchError, match="1.000000e\\+00"):
            phys.condensation_indices(phys.G295_EV, bad)

    def test_condense_and_errors(self):
        v = np.arange(1.0, 7.0)                      # 6 fine groups
        idx = np.array([2, 6])
        assert np.allclose(phys.condense(v, idx), [3.0, 18.0])
        err = np.ones(6)
        assert np.allclose(phys.condense_abs_error(err, idx, "linear"), [2, 4])
        assert np.allclose(phys.condense_abs_error(err, idx), [np.sqrt(2), 2])
        two_d = np.vstack([v, 2 * v])
        assert np.allclose(phys.condense(two_d, idx), [[3, 18], [6, 36]])

    def test_meshes_match_serpent_exporter(self):
        from starterDD.InterfaceToDD.Serpent2_exports import SUPPORTED_ENERGY_MESHES
        for name, bounds_mev in SUPPORTED_ENERGY_MESHES.items():
            # Same DRAGON library bounds, to full precision.
            assert np.allclose(np.array(bounds_mev) * 1e6, phys.mesh(name),
                               rtol=1e-12, atol=0), name

    def test_exporter_ene_card_reproduces_dragon_bounds(self):
        from starterDD.InterfaceToDD.Serpent2_exports import S2_EnergyGrid
        card = S2_EnergyGrid.from_preset("295g").format_card()
        written = np.array([float(x) for x in card.split()[3:]]) * 1e6
        assert written.size == 296
        assert np.allclose(written, phys.G295_EV, rtol=1e-8, atol=0)
        assert phys.mesh_code(written) == "295G" and not phys.is_legacy_mesh(written)

    def test_mesh_code(self):
        assert phys.mesh_code(phys.G2_EV) == "2G"
        assert phys.mesh_code(phys.G295_EV[::-1]) == "295G"
        assert phys.mesh_code([1e-5, 1.0, 10.0, 2e7]) == "C3G"

    @pytest.mark.parametrize("tok,name", [
        ("92235.09c", "U235"), ("92235", "U235"), ("64155.06c", "Gd155"),
        ("U235", "U235"), ("Gd-157", "Gd157"), ("1001.06c", "H1"),
        ("lwj3.11t", None), ("water", None), ("95342.09c", None),
    ])
    def test_nuclide_name(self, tok, name):
        assert phys.nuclide_name(tok) == name

    def test_reactions_and_q_values(self):
        assert phys.REACTIONS_ABS[0] == (102, "NG", 1)
        assert phys.REACTIONS_BY_MT[16].absorption_factor == -1
        assert phys.q_value("endfb8r1295", "U235", 18) == pytest.approx(194.077)
        assert phys.q_value("endfb8r1", "Gd155", 18) == 0.0
        with pytest.raises(KeyError):
            phys.q_value("jeff33", "U235", 18)


# ═════════════════════════════════════════════════════════════════════════════
# Synthetic DRAGON ASCII writer (test helper)
# ═════════════════════════════════════════════════════════════════════════════

def _hdr(level, ntype, ctype, size, idx=""):
    return (f"->{level:8d}{ntype:8d}{ctype:8d}{size:8d}" + " " * 33
            + "<-   " + idx + "\n")


class _Ascii:
    """Minimal writer reproducing the layout of DRAGON ASCII exports."""

    def __init__(self):
        self.out = []

    def _key(self, key):
        self.out.append(key.ljust(80) + "\n")

    def table(self, level, key=None, idx=None):
        if idx is None:
            self.out.append(_hdr(level, 12, 0, -1))
            self._key(key)
        else:
            self.out.append(_hdr(level, 0, 0, -1, f"{idx:08d}"))

    def lst(self, level, key, n):
        self.out.append(_hdr(level, 12, 10, n))
        self._key(key)

    def close(self, child_level):
        self.out.append(_hdr(-child_level, 0, 0, 0))

    def ints(self, level, key, vals, ctype=1):
        vals = list(np.atleast_1d(vals))
        self.out.append(_hdr(level, 12, ctype, len(vals)))
        self._key(key)
        for i in range(0, len(vals), 8):
            self.out.append("".join(f"{int(v):10d}" for v in vals[i:i + 8]) + "\n")

    def reals(self, level, key, vals, fmt="{:16.8E}"):
        vals = list(np.atleast_1d(vals))
        self.out.append(_hdr(level, 12, 2, len(vals)))
        self._key(key)
        for i in range(0, len(vals), 5):
            self.out.append("".join(fmt.format(v) for v in vals[i:i + 5]) + "\n")

    def raw_reals(self, level, key, tokens):
        self.out.append(_hdr(level, 12, 2, len(tokens)))
        self._key(key)
        self.out.append("".join(f"{t:>16s}" for t in tokens) + "\n")

    def chars(self, level, key, text):
        nw = -(-len(text) // 4)
        text = text.ljust(4 * nw)
        self.out.append(_hdr(level, 12, 3, nw))
        self._key(key)
        for i in range(0, nw, 8):
            self.out.append("".join(f"{4:10d}" for _ in range(min(8, nw - i))) + "\n")
        for i in range(0, len(text), 80):
            self.out.append(text[i:i + 80] + "\n")

    def undef(self, level, key):
        self.out.append(_hdr(level, 12, 99, 0))
        self._key(key)

    def text(self):
        return "".join(self.out)


def build_tree(leaves):
    """Build (NVP, DEBARB, ARBVAL) from ``{value_index_tuple: calc}``
    following DRAGON's convention (inverse of decode_tree)."""
    n_par = len(next(iter(leaves)))
    levels = [[()]]
    for k in range(1, n_par + 1):
        levels.append(sorted({t[:k] for t in leaves}))
    number = {}
    for lev in levels:
        for p in lev:
            number[p] = len(number) + 1
    nvp = len(number)
    deb = [0] * (nvp + 1)
    arb = [0] * nvp
    for p, i in number.items():
        if p:
            arb[i - 1] = p[-1]
        if len(p) < n_par:
            kids = [number[c] for c in levels[len(p) + 1] if c[:len(p)] == p]
            deb[i - 1] = min(kids)
            deb[i] = max(kids) + 1   # overwritten by the next node if internal
    for p, i in number.items():
        if len(p) == n_par:
            deb[i] = leaves[p]
    return nvp, deb, arb


TFUEL = [600.0, 900.0]
DCOOL = [0.2, 0.45, 0.7366899]
# Deliberately shuffled calculation numbering (insertion order != tree order)
CALC_OF = {(1, 3): 1, (1, 1): 2, (2, 2): 3, (1, 2): 4, (2, 1): 5, (2, 3): 6}
N_MIX = 3
ISOS = ["U235", "U238", "Gd155", "Gd157", "Xe135", "Sm149", "Pu239", "*MAC*RES"]


def keff_of(t, d):
    return 1.0 + t / 1.0e4 + d / 10.0


def write_synthetic_compo(path: Path):
    """Two outputs: OUT_2G (2 parameters, 6 calcs, 3 mixtures) and a
    single-statepoint HOM_1G. Group-wise values encode (mix, calc, group)."""
    w = _Ascii()
    w.chars(1, "SIGNATURE", "L_MULTICOMPO")
    for out, ng in (("OUT_2G", 2), ("HOM_1G", 1)):
        params = (("TFuel", TFUEL), ("DCool", DCOOL)) if out == "OUT_2G" \
            else (("TFuel", [900.0]), ("DCool", [0.7366899]))
        calc_of = CALC_OF if out == "OUT_2G" else {(1, 1): 1}
        nmix = N_MIX if out == "OUT_2G" else 1
        w.table(1, out)
        w.chars(2, "COMMENT", f"synthetic {out}")
        w.table(2, "GLOBAL")
        w.chars(3, "PARKEY", "".join(n.ljust(12) for n, _ in params))
        w.chars(3, "PARTYP", "VALU" * len(params))
        w.chars(3, "PARFMT", "REAL    " * len(params))
        w.ints(3, "NVALUE", [len(v) for _, v in params])
        for i, (_, vals) in enumerate(params, start=1):
            w.reals(3, f"pval{i:08d}", vals)
        w.close(3)
        w.lst(2, "MIXTURES", nmix)
        nvp, deb, arb = build_tree(calc_of)
        for m in range(1, nmix + 1):
            w.table(3, idx=m)
            w.lst(4, "CALCULATIONS", len(calc_of))
            for c in range(1, len(calc_of) + 1):
                tup = next(t for t, cc in calc_of.items() if cc == c)
                t, d = params[0][1][tup[0] - 1], params[1][1][tup[1] - 1]
                w.table(5, idx=c)
                isos = ISOS if m != 2 else ISOS[:2] + ISOS[4:]  # no Gd in mix 2
                w.lst(6, "ISOTOPESLIST", len(isos))
                for k, iso in enumerate(isos, start=1):
                    w.table(7, idx=k)
                    w.chars(8, "ALIAS", iso.ljust(8) + f"{k:04d}")
                    # NWT0 / NG stored DRAGON-style: group 1 = fast.
                    w.reals(8, "NWT0", [100 * m + 10 * c + g for g in range(1, ng + 1)])
                    if iso not in ("*MAC*RES", "Xe135"):
                        w.reals(8, "NG", [k + 0.1 * g for g in range(1, ng + 1)])
                    if iso == "U235":
                        w.reals(8, "NFTOT", [50.0 + g for g in range(1, ng + 1)])
                    w.close(8)
                w.chars(6, "ISOTOPESUSED", "".join(i.ljust(12) for i in isos))
                w.reals(6, "ISOTOPESDENS", [1e-3 * (k + 1) * m for k in range(len(isos))])
                w.reals(6, "ISOTOPESVOL", [2.5 * m] * len(isos))
                w.reals(6, "MIXTURESVOL", 2.5 * m)
                w.reals(6, "K-EFFECTIVE", keff_of(t, d))
                w.ints(6, "FLAGS", [1, 0, 1], ctype=5)   # logical record
                w.undef(6, "UNDEFINED")
                bounds = [1.964e7, 0.625, 1.1e-4] if ng == 2 else [1.964e7, 1.1e-4]
                w.reals(6, "ENERGY", bounds)
                w.close(6)
            w.table(4, "TREE")
            w.ints(5, "NVP", [nvp, nvp + 4])
            w.ints(5, "NCALS", len(calc_of))
            w.ints(5, "DEBARB", deb + [0, 0, 0, 0])   # padded like DRAGON
            w.ints(5, "ARBVAL", arb + [0, 0, 0, 0])
            w.close(5)
            w.close(4)
        w.table(2, "DEPL-CHAIN")
        # Fortran 3-digit exponent without "E" and a "D" exponent
        w.raw_reals(3, "TINY", ["1.0000000-100", "2.5000000D+01"])
        w.close(3)
        w.close(2)
    w.close(1)
    path.write_text(w.text())
    return path


@pytest.fixture
def compo_file(tmp_path):
    return write_synthetic_compo(tmp_path / "_CPO_SYNTH")


# ═════════════════════════════════════════════════════════════════════════════
# Pylcm
# ═════════════════════════════════════════════════════════════════════════════

class TestPylcm:
    def test_parse_structure_and_path_access(self, compo_file):
        lcm = Pylcm.from_ascii(compo_file, name="MULTICOMPO")
        assert lcm.keys() == ["SIGNATURE", "OUT_2G", "HOM_1G"]
        assert lcm["SIGNATURE"] == "L_MULTICOMPO"
        mix = lcm.get("OUT_2G/MIXTURES")
        assert mix.type == "list" and mix.size == N_MIX and mix.keys() == [1, 2, 3]
        c = CALC_OF[(2, 3)]
        k = lcm.get(f"OUT_2G/MIXTURES/1/CALCULATIONS/{c}/K-EFFECTIVE")
        assert k == pytest.approx(keff_of(900.0, 0.7366899))
        # Original chained access still works
        same = (lcm.content["OUT_2G"].content["MIXTURES"].content[1]
                .content["CALCULATIONS"].content[c].content["K-EFFECTIVE"])
        assert same == k
        assert "OUT_2G/GLOBAL/PARKEY" in lcm
        assert "OUT_2G/NOPE" not in lcm
        with pytest.raises(KeyError):
            lcm.get("OUT_2G/NOPE")

    def test_records_types(self, compo_file):
        lcm = Pylcm.from_ascii(compo_file)
        calc = lcm.get("OUT_2G/MIXTURES/1/CALCULATIONS/1")
        used = calc["ISOTOPESUSED"]
        assert isinstance(used, LcmChars)
        assert used.split() == ISOS             # multi-line char record
        assert used.chunks(12) == ISOS
        assert calc["FLAGS"].dtype == bool and list(calc["FLAGS"]) == [True, False, True]
        assert calc["UNDEFINED"] is None and "UNDEFINED" in calc
        assert isinstance(lcm.get("OUT_2G/GLOBAL/NVALUE"), np.ndarray)
        assert isinstance(lcm.get("OUT_2G/MIXTURES/1/TREE/NCALS"), int)
        tiny = lcm.get("OUT_2G/DEPL-CHAIN/TINY")
        assert tiny[0] == pytest.approx(1e-100) and tiny[1] == 25.0

    def test_selective_loading(self, compo_file):
        lcm = Pylcm.from_ascii(compo_file, only=["HOM_1G"])
        assert lcm.keys() == ["HOM_1G"]
        assert lcm.get("HOM_1G/MIXTURES/1/CALCULATIONS/1/K-EFFECTIVE") == \
            pytest.approx(keff_of(900.0, 0.7366899))

    def test_generate_from_ascii_compat(self, compo_file):
        lcm = Pylcm.generate_from_ascii(str(compo_file), name="MULTICOMPO")
        assert lcm.name == "MULTICOMPO" and "OUT_2G" in lcm

    def test_print_tree(self, compo_file, tmp_path):
        Pylcm.from_ascii(compo_file, only=["HOM_1G"]).print_tree(str(tmp_path), "tree")
        assert "ISOTOPESLIST" in (tmp_path / "tree.txt").read_text()

    @needs_h5
    def test_hdf5_round_trip_and_lazy_view(self, compo_file, tmp_path):
        lcm = Pylcm.from_ascii(compo_file)
        h5 = lcm.to_hdf5(tmp_path / "c.h5")
        back = Pylcm.from_hdf5(h5)
        a, b = list(lcm.walk()), list(back.walk())
        assert [p for p, _ in a] == [p for p, _ in b]       # order preserved
        for (p, x), (_, y) in zip(a, b):
            if isinstance(x, np.ndarray):
                assert np.array_equal(x, y) and x.dtype == y.dtype, p
            else:
                assert x == y and type(x) is type(y), p
        with LcmH5(h5) as view:
            assert view.keys() == ["SIGNATURE", "OUT_2G", "HOM_1G"]
            assert view.keys("OUT_2G/MIXTURES") == [1, 2, 3]
            k = view.get("OUT_2G/MIXTURES/2/CALCULATIONS/3/K-EFFECTIVE")
            assert k == lcm.get("OUT_2G/MIXTURES/2/CALCULATIONS/3/K-EFFECTIVE")
            assert view.get("OUT_2G/MIXTURES/1/CALCULATIONS/1/ISOTOPESUSED").chunks(12) == ISOS
            assert view.get("OUT_2G/MIXTURES/1/CALCULATIONS/1/UNDEFINED") is None
        part = Pylcm.from_hdf5(h5, only=["HOM_1G"])
        assert part.keys() == ["HOM_1G"]


# ═════════════════════════════════════════════════════════════════════════════
# CompoReader
# ═════════════════════════════════════════════════════════════════════════════

class TestDecodeTree:
    def test_single_statepoint_ge14_records(self):
        # Records of the GE14 sample (_CPO_GE14_DOM_DIAG, 3 parameters)
        assert decode_tree(4, [2, 3, 4, 5, 1], [0, 1, 1, 1], 3) == [((1, 1, 1), 1)]

    def test_one_parameter(self):
        # root 1, leaves 2..4 -> calcs 3, 1, 2
        assert decode_tree(4, [2, 5, 3, 1, 2], [0, 1, 2, 3], 1) == \
            [((1,), 3), ((2,), 1), ((3,), 2)]

    @pytest.mark.parametrize("seed", range(5))
    def test_random_trees_round_trip(self, seed):
        rng = random.Random(seed)
        n_par = rng.randint(1, 4)
        sizes = [rng.randint(1, 4) for _ in range(n_par)]
        full = list(itertools.product(*[range(1, s + 1) for s in sizes]))
        tuples = rng.sample(full, rng.randint(1, len(full)))   # sparse tree
        calcs = list(range(1, len(tuples) + 1))
        rng.shuffle(calcs)
        leaves = dict(zip(tuples, calcs))
        nvp, deb, arb = build_tree(leaves)
        got = decode_tree(nvp, deb + [0] * 3, arb + [0] * 3, n_par)
        assert dict(got) == leaves
        assert [t for t, _ in got] == sorted(leaves)          # tree order


class TestCompoReader:
    @pytest.fixture
    def reader(self, compo_file):
        return CompoReader(compo_file, cache=False,
                           labels={"00": {"DCool": 0.7366899, "TFuel": 900.0},
                                   "cold": {"TFuel": 600.0, "DCool": 0.2}})

    def test_outputs_and_parameters(self, reader):
        assert reader.outputs() == ["OUT_2G", "HOM_1G"]
        p = reader.parameters("OUT_2G")
        assert list(p) == ["TFuel", "DCool"]
        assert p["TFuel"] == TFUEL and np.allclose(p["DCool"], DCOOL)
        assert reader.n_mixtures("OUT_2G") == N_MIX
        assert reader.comment("OUT_2G") == "synthetic OUT_2G"

    def test_statepoints_decode_shuffled_tree(self, reader):
        sps = reader.statepoints("OUT_2G", mixture=2)
        assert [s.calc for s in sps] == list(range(1, 7))
        for s in sps:
            t, d = s.values["TFuel"], s.values["DCool"]
            # keff encodes the parameters: proves calc <-> tuple mapping
            assert reader.keff("OUT_2G", sp=s.calc) == pytest.approx(keff_of(t, d))
            assert CALC_OF[(s.indices["TFuel"], s.indices["DCool"])] == s.calc
        assert {s.label for s in sps} >= {"00", "cold"}

    def test_calc_index_selection(self, reader):
        assert reader.calc_index("OUT_2G", {"TFuel": 900, "DCool": 0.73669}) == CALC_OF[(2, 3)]
        assert reader.calc_index("OUT_2G", "cold") == CALC_OF[(1, 1)]
        assert reader.calc_index("OUT_2G", 4) == 4
        with pytest.raises(StatepointError, match="2 calculation"):
            reader.calc_index("OUT_2G", {"DCool": 0.2})          # ambiguous
        with pytest.raises(StatepointError, match="No"):
            reader.calc_index("OUT_2G", {"TFuel": 750.0, "DCool": 0.2})
        with pytest.raises(StatepointError):
            reader.keff("OUT_2G")                                 # 6 calcs, sp=None
        with pytest.raises(StatepointError):
            reader.calc_index("OUT_2G", 7)
        with pytest.raises(KeyError, match="unknown parameter"):
            reader.calc_index("OUT_2G", {"Burnup": 0.0})
        assert reader.keff("HOM_1G") == pytest.approx(keff_of(900.0, 0.7366899))

    def test_groupwise_accessors_ascending(self, reader):
        sp = "00"
        c = CALC_OF[(2, 3)]
        assert np.allclose(reader.energy_bounds("OUT_2G", sp), [1.1e-4, 0.625, 1.964e7])
        flux = reader.flux("OUT_2G", sp)
        assert flux.shape == (N_MIX, 2)
        # stored [g1 (fast), g2 (thermal)] = [100m+10c+1, 100m+10c+2] -> flipped
        assert np.allclose(flux[:, 0], [100 * m + 10 * c + 2 for m in (1, 2, 3)])
        ng = reader.micro_xs("OUT_2G", "U238", "NG", sp)
        assert np.allclose(ng, [[2.2, 2.1]] * N_MIX)
        gd = reader.micro_xs("OUT_2G", "Gd155", "NG", sp)
        assert np.allclose(gd[1], 0.0) and np.all(gd[[0, 2]] > 0)   # no Gd in mix 2
        assert reader.micro_xs("OUT_2G", "U238", "N2N", sp) is None
        with pytest.raises(KeyError):
            reader.micro_xs("OUT_2G", "U238", "N2N", sp, missing="raise")
        assert reader.micro_xs("OUT_2G", "U235", "NFTOT", sp, mixtures=3).shape == (2,)
        # ISOTOPESDENS = 1e-3 * (k+1) * m, U238 is k=1
        assert reader.densities("OUT_2G", sp, mixture=3)["U238"] == pytest.approx(6e-3)
        assert reader.isotope_volumes("OUT_2G", sp, mixture=2)["U235"] == pytest.approx(5.0)
        assert reader.mixture_volume("OUT_2G", sp, mixture=3) == pytest.approx(7.5)

    @needs_h5
    def test_cache_build_reuse_and_invalidation(self, compo_file):
        cache = default_cache_path(compo_file)
        with CompoReader(compo_file) as r:          # builds the cache
            assert isinstance(r.root, LcmH5)
            k = r.keff("OUT_2G", {"TFuel": 600, "DCool": 0.45})
        assert cache.exists() and cache_is_valid(compo_file, cache)
        mtime = cache.stat().st_mtime_ns
        with CompoReader(compo_file) as r:          # reuses it
            assert r.keff("OUT_2G", {"TFuel": 600, "DCool": 0.45}) == k
        assert cache.stat().st_mtime_ns == mtime
        # Touching the source without changing it: hash fallback keeps the cache
        os.utime(compo_file, None)
        assert cache_is_valid(compo_file, cache)
        # Changing the source invalidates it
        compo_file.write_text(compo_file.read_text().replace("synthetic OUT_2G",
                                                             "modified OUT_2G"))
        assert not cache_is_valid(compo_file, cache)
        with CompoReader(compo_file) as r:
            assert r.comment("OUT_2G") == "modified OUT_2G"
        assert isinstance(open_compo(cache), LcmH5)


# ═════════════════════════════════════════════════════════════════════════════
# Serpent2: synthetic legacy run (arbitrary names, include, comments)
# ═════════════════════════════════════════════════════════════════════════════

SERP_MAIN = """\
% Synthetic legacy input: detector names follow no convention.
set title "synthetic legacy run, void 40%"
set pop 1000 50 10
mat fuelA_zone_1_pin_1 sum tmp 900.0 burn 1
    92235.09c 1.0E-03  92238.09c 2.0E-02  8016.09c 4.0E-02
mat fuelA_zone_2_pin_1 sum tmp 900.0
    92235.09c 1.0E-03  92238.09c 2.0E-02  8016.09c 4.0E-02
mat gdfuel_zone_1_pin_2 sum tmp 900.0
    92235.09c 1.0E-03 92238.09c 2.0E-02 64155.09c 1.0E-04 8016.09c 4.0E-02
mat water sum tmp 600.0 moder lwtr 1001  1001.06c 0.05  8016.06c 0.025
/* response materials: names unrelated to the nuclide */
mat resp_a 1.0 tmp 900.0  92235.09c 1.0
mat U238 1.0 tmp 900.0  92238.09c 1.0
mat gd_r 1.0 tmp 900.0  64155.09c 1.0
therm lwtr lwj3.11t
ene grid_coarse 1  1.1E-10 6.25E-07 1.964E+01
ene grid_fine   1  1.1E-10 1.0E-07 6.25E-07 1.0E-01 1.964E+01
include "detectors.inc"
"""

SERP_DETS = """\
det weird_A  de grid_fine  dt -4  dm fuelA_zone_1_pin_1  dm fuelA_zone_2_pin_1
    dr 18 resp_a   dr 102 U238   % trailing comment
det weird_B
    de grid_fine
    dt -4
    dm gdfuel_zone_1_pin_2
    dr 18 resp_a
    dr 102 gd_r
    dr 102 U238
det fluxA de grid_fine dm fuelA_zone_1_pin_1 dm fuelA_zone_2_pin_1
det fluxB de grid_fine dm gdfuel_zone_1_pin_2
det everything_fuel de grid_coarse dt -4
    dm fuelA_zone_1_pin_1 dm fuelA_zone_2_pin_1 dm gdfuel_zone_1_pin_2
    dr 18 resp_a
det whole de grid_fine
det in_water dm water
"""

N_FINE = 4


def _det_rows(n_resp, n_ene, base):
    rows = []
    i = 0
    for r in range(1, n_resp + 1):
        for e in range(1, n_ene + 1):
            i += 1
            val = base + 10.0 * r + e
            rows.append(f"{i:5d} {e:4d} 1 1 1 1 {r:4d} 1 1 1  {val:.5E} {0.001 * e:.5f}")
    return rows


def write_synthetic_serpent(d: Path, drop_response_in_output=False):
    (d / "detectors.inc").write_text(SERP_DETS)
    serp = d / "legacy.serp"
    serp.write_text(SERP_MAIN)
    fine = [1.1e-10, 1.0e-07, 6.25e-07, 1.0e-01, 1.964e1]
    coarse = [1.1e-10, 6.25e-07, 1.964e1]
    spec = {"weird_A": (2, fine), "weird_B": (3, fine), "fluxA": (1, fine),
            "fluxB": (1, fine), "everything_fuel": (1, coarse),
            "whole": (1, fine), "in_water": (1, None)}
    if drop_response_in_output:
        spec["weird_B"] = (2, fine)
    out = ["", "% synthetic det0.m"]
    for k, (name, (nr, grid)) in enumerate(spec.items()):
        ne = 1 if grid is None else len(grid) - 1
        out.append(f"DET{name} = [")
        out += _det_rows(nr, ne, 100.0 * k)
        out.append("];")
        if grid is not None:
            out.append(f"DET{name}E = [")
            out += [f"{lo:.5E} {hi:.5E} {(lo + hi) / 2:.5E}" for lo, hi in zip(grid[:-1], grid[1:])]
            out.append("];")
    (d / "legacy.serp_det0.m").write_text("\n".join(out) + "\n")
    (d / "legacy.serp_res.m").write_text(
        "% Serpent results\n"
        "% Version, title and date:\n"
        "VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;\n"
        "TITLE                     (idx, [1: 30])  = 'synthetic legacy run = void 40%' ;\n"
        "POP                       (idx, 1)        = 1000 ;\n"
        "IMP_KEFF                  (idx, [1:   2]) = [  1.01234E+00 5.0E-05 ];\n"
        "ANA_KEFF                  (idx, [1:   6]) = [  1.01200E+00 6.0E-05  1.0E+00 1.0E-05 1.0E-02 1.0E-03 ];\n"
    )
    return serp


@pytest.fixture
def serp_file(tmp_path):
    return write_synthetic_serpent(tmp_path)


class TestSerpentInput:
    def test_parse_cards(self, serp_file):
        inp = parse_serpent_input(serp_file)
        assert inp.title == "synthetic legacy run, void 40%"
        assert set(inp.detectors) == {"weird_A", "weird_B", "fluxA", "fluxB",
                                      "everything_fuel", "whole", "in_water"}
        a = inp.detectors["weird_A"]
        assert a.grid == "grid_fine" and a.dtype == ["-4"]
        assert a.materials == ["fuelA_zone_1_pin_1", "fuelA_zone_2_pin_1"]
        assert a.responses == [(18, "resp_a"), (102, "U238")]
        assert a.file.endswith("detectors.inc")
        w = inp.materials["water"]
        assert w.density is None and w.temperature == 600.0
        assert w.options["moder"] == ["lwtr", "1001"]
        assert w.composition() == {"H1": 0.05, "O16": 0.025}
        assert inp.materials["resp_a"].single_nuclide == "U235"
        assert inp.materials["fuelA_zone_1_pin_1"].options["burn"] == ["1"]
        assert len(inp.grids["grid_fine"].bounds_eV) == N_FINE + 1


class TestManifestInference:
    def test_classification_from_card_content(self, serp_file):
        man = DetectorManifest.infer(serp_file)
        by = {e.detector: e for e in man}
        assert by["weird_A"].key == ("PIN_C4G", "rr", 1)
        assert by["weird_B"].key == ("PIN_C4G", "rr", 2)
        assert by["fluxA"].key == ("PIN_C4G", "flux", 1)
        assert by["everything_fuel"].key == ("FUEL_2G", "rr", None)
        assert by["everything_fuel"].pins == [1, 2]
        assert by["whole"].key == ("ALL_C4G", "flux", None)
        assert by["in_water"].key == ("MAT_1G", "flux", None)
        assert by["in_water"].n_groups == 1 and by["in_water"].bounds_eV is None
        # Response nuclides come from the response materials' mat cards
        assert [(r.mt, r.nuclide) for r in by["weird_B"].responses] == \
            [(18, "U235"), (102, "Gd155"), (102, "U238")]
        assert by["weird_B"].response_index(102, "Gd155") == 1
        assert man.pins("PIN_C4G", "rr") == [1, 2]
        assert "PIN_C4G" in man.summary()

    def test_custom_pin_resolution(self, serp_file):
        resolver = PinResolver(material_patterns=(),
                               material_to_pin={"fuelA_zone_1_pin_1": 7,
                                                "fuelA_zone_2_pin_1": 7,
                                                "gdfuel_zone_1_pin_2": 9,
                                                "water": 99})
        man = DetectorManifest.infer(serp_file, pin_resolver=resolver)
        assert man.by_detector("weird_A").pin_idx == 7
        assert man.by_detector("in_water").key == ("PIN_1G", "flux", 99)

    def test_overrides_exclude_and_conflicts(self, serp_file):
        # Make two detectors collide: unresolvable pins -> both MAT_C4G flux
        resolver = PinResolver(material_patterns=())
        with pytest.raises(ManifestConflictError, match="fluxA"):
            DetectorManifest.infer(serp_file, pin_resolver=resolver)
        man = DetectorManifest.infer(serp_file, pin_resolver=resolver,
                                     exclude=["fluxB"],
                                     overrides={"weird_B": {"observable": "GD_C4G"}})
        assert "fluxB" not in [e.detector for e in man]
        assert man.by_detector("weird_B").observable == "GD_C4G"
        assert man.by_detector("weird_B").source == "override"
        with pytest.raises(KeyError):
            DetectorManifest.infer(serp_file, exclude=["nope"])

    def test_json_round_trip(self, serp_file, tmp_path):
        man = DetectorManifest.infer(serp_file)
        p = man.write(tmp_path / "m.json")
        back = DetectorManifest.read(p)
        assert [e.to_dict() for e in back] == [e.to_dict() for e in man]
        assert back.source["input_sha256"] == man.source["input_sha256"]
        assert json.loads(p.read_text())["format"] == "starterDD-serpent-detectors"


class TestSerpentRun:
    def test_tallies_keff_and_condensation(self, serp_file):
        run = SerpentRun(serp_file)
        assert run.check() == []
        assert run.keff() == (1.01234, 5.0e-05)
        assert run.keff("ANA_KEFF")[0] == pytest.approx(1.012)
        t = run.tally("PIN_C4G", "rr", pin=2)
        assert t.values.shape == (3, N_FINE)
        v, e = t.response(102, "Gd155")
        # weird_B is the 2nd detector (base 100), reaction bin 2
        assert np.allclose(v, [100 + 20 + g for g in range(1, N_FINE + 1)])
        assert np.allclose(e, [0.001 * g for g in range(1, N_FINE + 1)])
        assert t.response(18, "Pu239") is None
        c = t.condense(phys.G2_EV)
        assert c.values.shape == (3, 2)
        assert np.allclose(c.values[1], [v[0] + v[1], v[2] + v[3]])
        abs_q = np.sqrt((v[0] * e[0]) ** 2 + (v[1] * e[1]) ** 2)
        assert c.rel_err[1, 0] == pytest.approx(abs_q / (v[0] + v[1]))
        assert run.tally("MAT_1G", "flux").bounds_eV is None
        with pytest.raises(phys.MeshMismatchError):
            t.condense([1.1e-4, 1.0, 1.964e7])

    def test_frozen_manifest_is_reused(self, serp_file):
        man = DetectorManifest.infer(serp_file, overrides={"whole": {"observable": "ALL_4G"}})
        man.write(default_manifest_path(serp_file))
        run = SerpentRun(serp_file)
        assert run.manifest.by_detector("whole").observable == "ALL_4G"

    def test_changed_input_warns(self, serp_file):
        DetectorManifest.infer(serp_file).write(default_manifest_path(serp_file))
        serp_file.write_text(serp_file.read_text() + "\n% edited\n")
        with pytest.warns(UserWarning, match="changed since"):
            SerpentRun(serp_file)

    def test_output_inconsistent_with_input_is_reported(self, tmp_path):
        serp = write_synthetic_serpent(tmp_path, drop_response_in_output=True)
        run = SerpentRun(serp)
        problems = run.check()
        assert len(problems) == 1 and "weird_B" in problems[0]
        with pytest.raises(ValueError, match="3 responses"):
            run.tally("PIN_C4G", "rr", pin=2)

    def test_det_reader_aux_arrays(self, serp_file):
        dets = read_det(serp_file.with_name("legacy.serp_det0.m"))
        assert "weird_AE" not in dets and "E" in dets["weird_A"].aux
        assert np.allclose(dets["weird_A"].energy_bounds_eV[[0, -1]], [0.11e-3, 1.964e7])
        assert dets["in_water"].energy_bounds_eV is None
        assert dets["weird_B"].values().shape == (4, 1, 1, 1, 1, 3, 1, 1, 1)

    def test_res_reader(self, serp_file):
        res = read_res(serp_file.with_name("legacy.serp_res.m"))
        assert res["TITLE"] == "synthetic legacy run = void 40%"
        assert res["POP"] == 1000 and isinstance(res["POP"], int)

    def test_legacy_wrappers(self, serp_file):
        det = serpent_simple_det_reader(str(serp_file.with_name("legacy.serp_det0.m")))
        assert "weird_AE" in det and det["weird_A"].shape[1] == 12
        info = serpent_simple_serp_reader(str(serp_file))
        assert info["weird_B"] == {"energy": 4, (18, "resp_a"): 1,
                                   (102, "gd_r"): 2, (102, "U238"): 3}
        assert info["in_water"]["energy"] == 1
        res = serpent_simple_res_reader(str(serp_file.with_name("legacy.serp_res.m")))
        assert np.allclose(res["IMP_KEFF"], [1.01234, 5e-5])


# ═════════════════════════════════════════════════════════════════════════════
# Tier 2: committed GE14 subset (tests/data/ge14_mini, built by build_fixture.py)
# ═════════════════════════════════════════════════════════════════════════════

MINI_CPO = GE14_MINI_DIR / "_CPO_GE14_DOM_DIAG_mini.gz"
MINI_SERP = GE14_MINI_DIR / "GE14_DOM_00_endfb8r1.serp.gz"
MINI_PINS = [1, 2, 15, 28]
# Legacy (golden) per-pin factor: 2 on the diagonal, 1 elsewhere. See the
# P1 open point in the plan before relying on this convention.
GE14_DOM_DIAG_AXIS_PINS = {1, 11, 20, 28, 33, 37, 42, 46, 49, 51}
GE14_FUEL_RADIUS = 0.438


class TestGE14Mini:
    @pytest.fixture(scope="class")
    def golden(self):
        return {
            "ref": json.loads((GE14_MINI_DIR / "golden" / "reference_values.json").read_text()),
            "npz": np.load(GE14_MINI_DIR / "golden" / "2G_data.npz"),
        }

    def test_compo_from_gzip(self, golden, tmp_path):
        r = CompoReader(MINI_CPO, cache=tmp_path / "mini.h5",
                        labels={"00": {"DCool": 0.73669}})
        assert r.outputs() == ["EDIHOM_COND", "EDIR_2G"]
        k = r.keff("EDIHOM_COND", sp="00")
        assert k == golden["ref"]["keff_dragon_EDIHOM_COND"] == pytest.approx(0.989089489, abs=1e-12)
        assert r.n_mixtures("EDIR_2G") == 51
        assert np.allclose(r.energy_bounds("EDIR_2G"), [1.100027e-4, 0.624998689, 1.96403e7])
        assert r.micro_xs("EDIR_2G", "U238", "NG", "00").shape == (51, 2)
        assert r.flux("EDIR_2G", "00").shape == (51, 2)
        r.close()

    def test_legacy_serpent_with_frozen_manifest(self, golden):
        run = SerpentRun(MINI_SERP)          # picks up the committed manifest
        assert run.det_path.name.endswith("_det0.m.gz")
        assert run.manifest.source["subset_pins"] == MINI_PINS
        assert run.check() == []
        assert run.manifest.pins("PIN_295G", "rr") == MINI_PINS
        assert run.manifest.pins("PIN_295G", "flux") == MINI_PINS
        assert run.manifest.find("FUEL_295G", "rr").detector == "det_assembly_295g"
        assert run.manifest.find("ALL_295G", "flux").detector == "flux_295g"
        k, dk = run.keff()
        assert (k, dk) == (golden["ref"]["keff_serpent_IMP_KEFF"],
                           golden["ref"]["keff_serpent_IMP_KEFF_rel_err"])
        assert (golden["ref"]["keff_dragon_EDIHOM_COND"] - k) * 1e5 == pytest.approx(-117.1511, abs=1e-3)

    def test_manifest_inference_matches_frozen_subset(self):
        inferred = DetectorManifest.infer(MINI_SERP)
        frozen = SerpentRun(MINI_SERP).manifest
        for e in frozen:
            assert inferred.by_detector(e.detector).to_dict() == e.to_dict()

    def test_reproduces_golden_serpent_pin_rates(self, golden):
        """The fixture is enough to rebuild the legacy 2G_data.npz Serpent
        columns of its pins (fission rates and linear-condensed errors)."""
        import yaml
        comp = {m["name"]: m["isotopic_composition"] for m in
                yaml.safe_load(open(GE14_MATERIAL_COMPOS_00_YAML))["MIX_COMPOSITIONS"]}
        npz = golden["npz"]
        run = SerpentRun(MINI_SERP)
        vol = np.pi * GE14_FUEL_RADIUS ** 2
        for p in MINI_PINS:
            e = run.manifest.find("PIN_295G", "rr", p)
            mat = e.materials[0].split("_zone_")[0]
            t = run.tally_of(e).condense(phys.G2_EV, error_mode="linear")
            sym = 2 if p in GE14_DOM_DIAG_AXIS_PINS else 1
            assert npz["symmetry_factors"][p - 1] == sym
            rate = np.zeros(2)
            err = np.zeros(2)
            for iso, zaid in (("U235", "92235"), ("U238", "92238")):
                v, rel = t.response(18, iso)
                rate += v * comp[mat][zaid] * vol * sym
                err += v * rel * comp[mat][zaid] * vol * sym
            assert np.allclose(rate, npz["serp_groups_fiss_rate_tot"][:, p - 1], rtol=1e-12)
            assert np.allclose(err, npz["serp_groups_fiss_rate_err_tot"][:, p - 1], rtol=1e-12)


# ═════════════════════════════════════════════════════════════════════════════
# Tier 3: full local data (never committed), marked ``fulldata``
# ═════════════════════════════════════════════════════════════════════════════

GE14_DOM = FULL_DATA_DIR / "GE14_DOM"
GE14_CPO = (GE14_DOM / "results" / "coolant_density_1L_PT_3x3_HFP" / "GE14_DOM_DIAG"
            / "PT_IC_1L_MOC" / "2026-07-19T10-46-27" / "_CPO_GE14_DOM_DIAG")
GE14_SERP = GE14_DOM / "Serpent" / "GE14_DOM_00_endfb8r1.serp"
DONJON_DATA = DONJON_DATA_DIR


@pytest.mark.fulldata
@pytest.mark.skipif(not GE14_CPO.exists(), reason="GE14 sample COMPO not present")
def test_sample_ge14_keff_vs_serpent():
    """P0 acceptance: keff 0.989089489 (DOM 1L, 0% void), -117.15 pcm vs
    Serpent, as in sortie_DOM_1L_PT_3x3_HFP.csv."""
    lcm = Pylcm.from_ascii(GE14_CPO, only=["EDIHOM_COND", "EDIR_2G"])
    r = CompoReader(lcm, labels={"00": {"DCool": 0.73669}})
    assert r.parameters("EDIR_2G") == {"TFuel": [900.0], "TCool": [600.0],
                                       "DCool": [pytest.approx(0.73669, rel=1e-6)]}
    k = r.keff("EDIHOM_COND", sp="00")
    assert k == pytest.approx(0.989089489, abs=1e-12)
    assert r.n_mixtures("EDIR_2G") == 51
    assert np.allclose(r.energy_bounds("EDIR_2G"), [1.100027e-4, 0.624998689, 1.96403e7])
    if GE14_SERP.exists():
        ks, _ = SerpentRun(GE14_SERP).keff()
        assert (k - ks) * 1e5 == pytest.approx(-117.1511, abs=1e-3)


@pytest.mark.fulldata
@pytest.mark.skipif(not GE14_SERP.exists(), reason="GE14 Serpent sample not present")
def test_sample_ge14_legacy_serpent_manifest():
    run = SerpentRun(GE14_SERP, manifest=DetectorManifest.infer(GE14_SERP))
    m = run.manifest
    assert len(m.select("PIN_295G", "rr")) == 51
    assert len(m.select("PIN_295G", "flux")) == 51
    assert m.find("FUEL_295G", "rr").detector == "det_assembly_295g"
    assert m.find("ALL_295G", "flux").detector == "flux_295g"
    gd = m.find("PIN_295G", "rr", pin=15)
    assert gd.detector == "det_295g_UOX40Gd8_15" and len(gd.materials) == 6
    assert run.check() == []
    t = run.tally("PIN_295G", "rr", pin=1).condense(phys.G2_EV)
    assert t.values.shape == (17, 2)


@pytest.mark.fulldata
@pytest.mark.skipif(not (DONJON_DATA / "rep900_proc" / "_ACompo_rep900").exists(),
                    reason="Donjon test data not found (set DONJON_DATA_DIR)")
def test_donjon_multiparameter_tree_physics():
    """Tree decoding on a real 24-calculation COMPO: Pu239 must increase
    with the decoded burnup and B10 with the decoded boron concentration."""
    r = CompoReader(Pylcm.from_ascii(DONJON_DATA / "rep900_proc" / "_ACompo_rep900"))
    out = "EDI2B"
    p = r.parameters(out)
    assert [len(p[k]) for k in ("burnup", "ppmBore", "TF", "TCA", "DCA")] == [6, 4, 1, 1, 1]
    assert len(r.statepoints(out)) == 24
    for ppm in p["ppmBore"]:
        sps = sorted(r.find_statepoints(out, ppmBore=ppm), key=lambda s: s.values["burnup"])
        pu = [r.densities(out, sp=s.calc)["Pu239"] for s in sps]
        assert np.all(np.diff(pu) > 0)
    for bu in p["burnup"]:
        sps = sorted(r.find_statepoints(out, burnup=bu), key=lambda s: s.values["ppmBore"])
        b10 = [r.densities(out, sp=s.calc)["B10"] for s in sps]
        assert np.all(np.diff(b10) > 0)


@pytest.mark.fulldata
@pytest.mark.skipif(not (DONJON_DATA / "compo" / "MultiCompo").exists(),
                    reason="Donjon test data not found (set DONJON_DATA_DIR)")
def test_donjon_string_parameters():
    r = CompoReader(Pylcm.from_ascii(DONJON_DATA / "compo" / "MultiCompo"))
    assert r.parameters("MODE") == {"MTYPE": ["CELL18", "CELL20"]}
    assert len(r.statepoints("FUEL")) == 206
    assert len(r.find_statepoints("MODE", MTYPE="CELL20")) == 1
