"""Tests for src/chord_invariants.py.  Run:  .venv/bin/python -m pytest tests/ -q   (or python tests/test_chord_invariants.py)"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from chord_invariants import (BPSProjector, crossing_number, free_compression_moments, free_cumulants,  # noqa: E402
                              nc_partitions, pair_invariants, single_invariants, six_letter_words,
                              triple_invariants)


def test_nc_partition_counts_are_catalan():
    assert [len(nc_partitions(n)) for n in range(1, 7)] == [1, 2, 5, 14, 42, 132]


def test_free_cumulants_semicircle_and_bernoulli():
    k = free_cumulants([0, 1, 0, 2, 0, 5])            # standard semicircle
    assert np.allclose(k, [0, 1, 0, 0, 0, 0])
    k = free_cumulants([0.5] * 4)                     # Bern(1/2): kappa2 = 1/4, kappa4 = -1/16
    assert np.allclose(k[:4], [0.5, 0.25, 0.0, -1 / 16])


def test_free_compression_reproduces_wachter_moments():
    a, b = 0.37, 0.61
    m = free_compression_moments(b, a, 3)
    assert np.isclose(m[0], b)
    assert np.isclose(m[1], b * (a + b - a * b))
    assert np.isclose(m[2], b * (a**2 + b**2 + 3*a*b - 3*a**2*b - 3*a*b**2 + 2*a**2*b**2))


def test_crossing_numbers():
    assert crossing_number("iijjkk") == 0 and crossing_number("ijijkk") == 1
    assert crossing_number("ijikjk") == 2 and crossing_number("ijkijk") == 3
    assert len(six_letter_words()) == 16


def _compare_methods(Np, P, q, haar=False):
    pb = BPSProjector(Np, P, q=q, seed=3, method="basis", haar=False)
    pc = BPSProjector(Np, P, q=q, seed=3, method="complement", haar=False)
    assert pb.d == pc.d
    x1, x2, x3 = pb.slot([Np - 1]), pb.slot([0]), pb.slot([1, 2])
    s1, lam_b = single_invariants(pb, x1)
    s2, lam_c = single_invariants(pc, x1)
    assert np.allclose(lam_b, lam_c, atol=1e-9)
    assert np.isclose(s1["mean"], pb.tau([x1]))
    p1, p2 = pair_invariants(pb, x1, x2), pair_invariants(pc, x1, x2)
    for key in ("c", "F", "X"):
        assert np.isclose(p1[key], p2[key], atol=1e-9), key
    t1, t2 = triple_invariants(pb, x1, x2, x3), triple_invariants(pc, x1, x2, x3)
    assert np.allclose([w["ratio"] for w in t1], [w["ratio"] for w in t2], atol=1e-8)


def test_basis_and_complement_agree_q3():
    _compare_methods(8, 4, 3)


def test_basis_and_complement_agree_q5():
    _compare_methods(11, 5, 5)


def test_eigh_nullspace_matches_svd():
    """null_method='eigh' (H = {Q, Q^dag}) and 'svd' give the same BPS projector and invariants."""
    for Np, q in ((10, 3), (11, 5)):
        ps = BPSProjector(Np, Np // 2, q=q, seed=4, method="basis", null_method="svd")
        pe = BPSProjector(Np, Np // 2, q=q, seed=4, method="basis", null_method="eigh")
        assert ps.d == pe.d
        Ps, Pe = ps.B @ ps.B.conj().T, pe.B @ pe.B.conj().T
        assert np.abs(Ps - Pe).max() < 1e-9
        x1, x2 = ps.slot([Np - 1]), ps.slot([0])
        a, b = pair_invariants(ps, x1, x2), pair_invariants(pe, x1, x2)
        assert all(np.isclose(a[k], b[k], atol=1e-10) for k in ("c", "F", "X"))


def test_tiny_cache_gives_same_words():
    """A byte budget that forces constant eviction must not change any word trace."""
    big = BPSProjector(11, 5, q=5, seed=1, method="complement")
    small = BPSProjector(11, 5, q=5, seed=1, method="complement", cache_gb=1e-6)
    x1, x2, x3 = big.slot([10]), big.slot([0]), big.slot([1])
    t1, t2 = triple_invariants(big, x1, x2, x3), triple_invariants(small, x1, x2, x3)
    assert np.allclose([w["ratio"] for w in t1], [w["ratio"] for w in t2], atol=1e-12)
    assert small._cache.nbytes <= max(v.nbytes for v in small._cache.d.values()) + 1


def test_decoder_matches_decoder_fast():
    from decoder_fast import decoder_spectrum
    lam_ref, d = decoder_spectrum(10, 5, 2)
    pb = BPSProjector(10, 5, q=3, seed=2, method="basis")
    _, lam = single_invariants(pb, pb.slot([9]))
    assert pb.d == d and np.allclose(np.sort(lam), np.sort(lam_ref), atol=1e-9)


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
