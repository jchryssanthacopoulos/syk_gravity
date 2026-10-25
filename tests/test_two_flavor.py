"""Tests for the two-flavor model and its exact pair uplift (src/two_flavor.py, src/two_flavor_uplift.py; D8)."""
import os
import sys
from math import factorial

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from two_flavor import TwoFlavor, random_generic_couplings, random_symmetric_couplings, vacuum_tower  # noqa: E402
from two_flavor_uplift import CHANNELS, SLOTS, PairUplift  # noqa: E402


def pair(N, seed, generic=False):
    rng = np.random.default_rng(seed)
    return PairUplift((random_generic_couplings if generic else random_symmetric_couplings)(N + 1, rng))


def sectors(N):
    return [(a, b) for a in range(N + 1) for b in range(N + 1)]


def close(x, y, rel=0.0, abs_=0.0):
    return abs(x - y) <= max(abs_, rel * abs(y))


# ------------------------------------------------------------------------------------------------ the model
def test_nilpotent():
    rng = np.random.default_rng(1)
    for generic in (False, True):
        m = TwoFlavor((random_generic_couplings if generic else random_symmetric_couplings)(4, rng))
        for a, b in sectors(4):
            if m.dim(a, b) and m.dim(a + 2, b - 1) and m.dim(a + 4, b - 2):
                assert np.abs((m.Q(a + 2, b - 1) @ m.Q(a, b)).toarray()).max() < 1e-12


def test_V_commutes_only_when_symmetrized():
    rng = np.random.default_rng(2)
    for generic, expect in ((False, True), (True, False)):
        m = TwoFlavor((random_generic_couplings if generic else random_symmetric_couplings)(4, rng))
        worst = 0.0
        for a, b in sectors(4):
            if m.dim(a, b) and m.dim(a + 3, b) and m.dim(a + 2, b - 1) and m.dim(a + 1, b + 1):
                comm = m.Q(a + 1, b + 1) @ m.V(a, b) - m.V(a + 2, b - 1) @ m.Q(a, b)
                worst = max(worst, np.abs(comm.toarray()).max())
        assert (worst < 1e-12) == expect


def test_tower_is_bps_with_cv_norms():
    rng = np.random.default_rng(3)
    N = 5
    m = TwoFlavor(random_symmetric_couplings(N, rng))
    for n in range(N + 1):
        t = vacuum_tower(m, n)
        assert abs(np.linalg.norm(t) ** 2 - factorial(n) * factorial(N) / factorial(N - n)) < 1e-8 * np.linalg.norm(t) ** 2
        assert close(np.linalg.norm(m.bps(n, n).conj().T @ t) ** 2, np.linalg.norm(t) ** 2, rel=1e-10)


# ------------------------------------------------------------------------------------------------ block form and identities
def test_block_form_of_new_supercharge():
    """Q_new = Q0 + psi_e^dag E + K chi_f + psi_e^dag Z chi_f as an operator identity on the new Fock space."""
    for generic in (False, True):
        _block_form(generic)


def _block_form(generic):
    P = pair(3, 10, generic)
    N, new = P.N, P.new
    e, f = 2 * N, 2 * N + 1
    # the old modes keep their bit positions in the new model (the new pair sits at bits 2N, 2N+1)
    assert new.psi_pos[:N] == P.old.psi_pos and new.chi_pos[:N] == P.old.chi_pos
    terms = list(new.Q_terms(sites=range(N)))
    terms += [(c, [(e, True)] + ops) for c, ops in P._pair_terms("E")]
    terms += [(c, ops + [(f, False)]) for c, ops in P._pair_terms("K")]
    terms += [(c, [(e, True)] + ops + [(f, False)]) for c, ops in P._pair_terms("Z")]
    for A in range(N + 2):
        for B in range(N + 2):
            if new.dim(A, B) and new.dim(A + 2, B - 1):
                diff = new.op_block(terms, (A, B), (A + 2, B - 1)) - new.Q(A, B)
                assert np.abs(diff.toarray()).max() < 1e-12


def test_pair_operator_identities():
    """[Q0,E] = 0, [Q0,K] = 0, [E,K] = {Q0,Z} (both models); [E,V] = K, [K,V] = 0 (symmetrized only)."""
    for generic in (False, True):
        _pair_identities(generic)


def _pair_identities(generic):
    P = pair(4, 11, generic)
    old = P.old
    worst = dict(QE=0, QK=0, EKZ=0, EV=0, KV=0)
    for a, b in sectors(P.N):
        if not P.old_dim(a, b):
            continue
        if P.old_dim(a + 3, b - 2):
            worst["QE"] = max(worst["QE"], np.abs(P.old_Q(a + 1, b - 1) @ P.E(a, b) - P.E(a + 2, b - 1) @ P.old_Q(a, b)).max())
        if P.old_dim(a + 4, b - 1):
            worst["QK"] = max(worst["QK"], np.abs(P.old_Q(a + 2, b) @ P.K(a, b) - P.K(a + 2, b - 1) @ P.old_Q(a, b)).max())
        if P.old_dim(a + 3, b - 1):                                   # (a,b) -> (a+3,b-1)
            EK = P.E(a + 2, b) @ P.K(a, b) - P.K(a + 1, b - 1) @ P.E(a, b)
            QZ = P.old_Q(a + 1, b) @ P.Z(a, b) + P.Z(a + 2, b - 1) @ P.old_Q(a, b)
            worst["EKZ"] = max(worst["EKZ"], np.abs(EK - QZ).max())
        if P.old_dim(a + 2, b):                                       # [E, V]: (a,b) -> (a+2,b)
            EV = P.E(a + 1, b + 1) @ old.V(a, b).toarray() - old.V(a + 1, b - 1).toarray() @ P.E(a, b) \
                if P.old_dim(a + 1, b + 1) and P.old_dim(a + 1, b - 1) else None
            if EV is not None:
                worst["EV"] = max(worst["EV"], np.abs(EV - P.K(a, b)).max())
        if P.old_dim(a + 3, b + 1):
            KV = P.K(a + 1, b + 1) @ old.V(a, b).toarray() - old.V(a + 2, b).toarray() @ P.K(a, b)
            worst["KV"] = max(worst["KV"], np.abs(KV).max())
    for k in ("QE", "QK", "EKZ", "KV"):
        assert worst[k] < 1e-12, k
    assert (worst["EV"] < 1e-12) == (not generic)


# ------------------------------------------------------------------------------------------------ exact uplift
def test_two_step_dimension_formula():
    for N, generic in ((3, False), (4, False), (4, True)):
        P = pair(N, 20 + N, generic)
        for A in range(N + 2):
            for B in range(N + 2):
                r = P.two_step_dimensions(A, B)
                assert r["M_pred"] == r["M"] and r["new_pred"] == r["new"], (A, B, r)


def test_slot_equations():
    """The eight slot equations (D8 eq. S, operator convention) hold on every new BPS state."""
    for generic in (False, True):
        _slot_equations(generic)


def _slot_equations(generic):
    P = pair(4, 30, generic)
    worst = 0.0
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            X = P.new_B(A, B)
            if X.size == 0 or X.shape[1] == 0:
                continue
            S = P.slot_components(A, B, X, operator_convention=True)
            P00, P10, P01, P11 = S[(0, 0)], S[(1, 0)], S[(0, 1)], S[(1, 1)]

            def op(f, *args):
                return f(*args) if P.old_ok(*args[:2]) else None

            eqs = []
            Q0 = P.old_Q
            # (Q00) Q0 P00 + K P01 = 0                    in old (A+2, B-1)
            if P.old_dim(A + 2, B - 1):
                eqs.append((Q0(A, B) @ P00 if P00.size else 0) + (P.K(A, B - 1) @ P01 if P01.size else 0))
            # (Q10) E P00 - Q0 P10 + Z P01 - K P11 = 0     in old (A+1, B-1)
            if P.old_dim(A + 1, B - 1):
                eqs.append((P.E(A, B) @ P00 if P00.size else 0) - (Q0(A - 1, B) @ P10 if P10.size else 0)
                           + (P.Z(A, B - 1) @ P01 if P01.size else 0) - (P.K(A - 1, B - 1) @ P11 if P11.size else 0))
            # (Q11) Q0 P11 + E P01 = 0                    in old (A+1, B-2)
            if P.old_dim(A + 1, B - 2):
                eqs.append((Q0(A - 1, B - 1) @ P11 if P11.size else 0) + (P.E(A, B - 1) @ P01 if P01.size else 0))
            # (D00) Q0^dag P00 + E^dag P10 = 0           in old (A-2, B+1)
            if P.old_dim(A - 2, B + 1):
                eqs.append((Q0(A - 2, B + 1).conj().T @ P00 if P00.size else 0)
                           + (P.E(A - 2, B + 1).conj().T @ P10 if P10.size else 0))
            # (D01) Q0^dag P01 = K^dag P00 + Z^dag P10 + E^dag P11   in old (A-2, B)
            if P.old_dim(A - 2, B):
                eqs.append((Q0(A - 2, B).conj().T @ P01 if P01.size else 0)
                           - (P.K(A - 2, B).conj().T @ P00 if P00.size else 0)
                           - (P.Z(A - 2, B).conj().T @ P10 if P10.size else 0)
                           - (P.E(A - 2, B).conj().T @ P11 if P11.size else 0))
            # (D11) Q0^dag P11 = K^dag P10                in old (A-3, B)
            if P.old_dim(A - 3, B):
                eqs.append((Q0(A - 3, B).conj().T @ P11 if P11.size else 0)
                           - (P.K(A - 3, B).conj().T @ P10 if P10.size else 0))
            for e in eqs:
                if not np.isscalar(e):
                    worst = max(worst, np.abs(e).max())
    assert worst < 1e-10


def test_four_slot_injectivity():
    for generic in (False, True):
        P = pair(4, 40, generic)
        for A in range(P.N + 2):
            for B in range(P.N + 2):
                assert P.injectivity_deficit(A, B) == 0


def test_slot_complex_E2_bounds_and_d2():
    """dim B_new <= dim E2 (slot complex of old BPS data); equality fails exactly where d2 acts (N = 6: four sectors)."""
    P = pair(6, 83)
    strict = []
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            if P.new_dim(A, B) == 0:
                continue
            e2 = sum(P.slot_complex_E2(A, B))
            assert P.d_new(A, B) <= e2
            if P.d_new(A, B) < e2:
                strict.append((A, B))
    assert sorted(strict) == [(2, 3), (3, 5), (4, 2), (5, 4)], strict


def test_theorem_A_explicit_step1():
    """Theorem A for step 1 (psi_e, creation type), solved explicitly: admissible harmonic sources give a unique BPS state
    of the intermediate model M; obstructed sources give none; B_M splits orthogonally into B_h (+) B_g."""
    from two_flavor_uplift import nullspace, orth_complement, rank
    for generic in (False, True):
        P = pair(4, 90, generic)
        rng = np.random.default_rng(5)
        for A in range(P.N + 2):
            for B in range(P.N + 1):
                idxM = P._M_idx(A, B)
                if len(idxM) == 0:
                    continue
                masks = P.new.sector(A, B)[idxM]
                low = masks & (P.be - 1)
                empty, filled = np.where((masks & P.be) == 0)[0], np.where((masks & P.be) != 0)[0]
                d = len(masks)

                def slot_basis(rows, a, b):
                    Y = np.zeros((d, P.b(a, b) if P.old_ok(a, b) else 0), complex)
                    if Y.shape[1]:
                        Y[rows] = P.old_B(a, b)[np.searchsorted(P.old.sector(a, b), low[rows])]
                    return Y
                H0, H1 = slot_basis(empty, A, B), slot_basis(filled, A - 1, B)
                L = [x for x in (P._QM(A, B), P._QM(A - 2, B + 1).conj().T) if x.size]
                if not L:
                    continue
                L = np.vstack(L)
                Nm = orth_complement(np.hstack([H0, H1]), d)
                # admissible / obstructed sources (Theorem A (b))
                tgt_h, tgt_g = P.old_B(A + 1, B - 1), P.old_B(A - 2, B + 1)
                Mh = tgt_h.conj().T @ P.E(A, B) @ P.old_B(A, B) if (H0.shape[1] and tgt_h.size and tgt_h.shape[1]) else None
                Mg = tgt_g.conj().T @ P.E(A - 2, B + 1).conj().T @ P.old_B(A - 1, B) if (H1.shape[1] and tgt_g.size and tgt_g.shape[1]) else None
                kh = nullspace(Mh, H0.shape[1]) if Mh is not None else np.eye(H0.shape[1])
                kg = nullspace(Mg, H1.shape[1]) if Mg is not None else np.eye(H1.shape[1])
                if Nm.shape[1]:
                    assert rank(L @ Nm) == Nm.shape[1]                       # uniqueness
                for H, k in ((H0, kh), (H1, kg)):
                    if H.shape[1] == 0:
                        continue
                    if k.shape[1]:                                           # admissible: solvable
                        src = H @ (k @ rng.standard_normal(k.shape[1]))
                        c, *_ = np.linalg.lstsq(L @ Nm, -L @ src, rcond=None)
                        assert np.linalg.norm(L @ (src + Nm @ c)) < 1e-10 * max(1, np.linalg.norm(src))
                    if k.shape[1] < H.shape[1]:                              # obstructed: no solution
                        kc = orth_complement(k, H.shape[1])
                        src = H @ (kc @ rng.standard_normal(kc.shape[1]))
                        c, *_ = np.linalg.lstsq(L @ Nm, -L @ src, rcond=None)
                        assert np.linalg.norm(L @ (src + Nm @ c)) > 1e-3 * np.linalg.norm(src)
                # orthogonal splitting of B_M
                X = P.M_B(A, B)
                if X.size and X.shape[1] and H0.shape[1] and H1.shape[1]:
                    Bh = X @ nullspace(H1.conj().T @ X, X.shape[1])
                    Bg = X @ nullspace(H0.conj().T @ X, X.shape[1])
                    assert Bh.shape[1] + Bg.shape[1] == X.shape[1]
                    if Bh.shape[1] and Bg.shape[1]:
                        assert np.abs(Bh.conj().T @ Bg).max() < 1e-10


# ------------------------------------------------------------------------------------------------ the tower
def test_tower_is_bare_with_exact_slot_weights():
    P = pair(5, 50)
    N = P.N
    for n in range(N + 2):
        w, bare = P.tower_slot_weights(n)
        assert close(bare, 1.0, abs_=1e-10)
        assert close(w["00"], (N + 1 - n) / (N + 1), abs_=1e-12)
        assert close(w["11"], n / (N + 1), abs_=1e-12)
        assert w["10"] < 1e-20 and w["01"] < 1e-20


def test_tower_cost_attains_cv_bound_iff_kernel_trivial():
    P = pair(5, 4005)
    for n in range(1, P.N + 1):
        for ch in ("00", "11"):
            r = P.tower_cost(n, ch)
            assert r["kappa"] <= r["bound"] + 1e-10
            assert (abs(r["kappa"] - r["bound"]) < 1e-9 * max(1, r["bound"])) == (r["kerD"] == 0), r


# ------------------------------------------------------------------------------------------------ arithmetic of exact lifts
def test_room_count_predicts_exact_lift_loci():
    """rk D_T = min(d_new, R_T) and dim(B cap im D_T) = max(0, b + rk D_T - R_T) [D8 sec. 6; numerical law]."""
    for N in (3, 4):
        _room_count(N)


def _room_count(N):
    P = pair(N, 60 + N)
    for a, b in sectors(N):
        rec = P.cv_channels(a, b)
        if rec is None:
            continue
        for ch in CHANNELS:
            c = rec[ch]
            assert c["rkD"] == c["rk_pred"], (a, b, ch, c)
            assert c["dimA"] == c["a_gen"], (a, b, ch, c)
            if a == b and ch in ("00", "11"):
                assert close(c["U_tow"], 1.0, abs_=1e-9)


def test_bare_count_is_generic_count_plus_tower():
    """N = 6 -> 7: bare states = generic count everywhere except the central diagonal sectors, where the tower is the
    only bare state (generic count 0) [D8 sec. 5; numerical law]."""
    for generic in (False, True):
        P = pair(6, 11, generic)
        for A in range(P.N + 2):
            for B in range(P.N + 2):
                if P.new_dim(A, B) == 0:
                    continue
                bare, gen = P.bare_count(A, B), P.bare_count_generic(A, B)["gen"]
                expect = gen + (1 if (not generic and A == B and gen == 0 and 1 <= A <= P.N) else 0)
                if A == B == 0 or A == B == P.N + 1:
                    expect = gen
                assert bare == expect, (generic, A, B, bare, gen)


# ------------------------------------------------------------------------------------------------ first order
def test_first_order_single_slot_is_two_point_function():
    P = pair(4, 70)
    rng = np.random.default_rng(0)
    checked = 0
    for A in range(P.N + 2):
        for B in range(P.N + 2):
            for T in SLOTS:
                a, b = A - T[0], B - T[1]
                if not P.old_ok(a, b) or P.b(a, b) == 0 or P.new_dim(A, B) == 0:
                    continue
                h = P.old_B(a, b) @ (rng.standard_normal(P.b(a, b)) + 1j * rng.standard_normal(P.b(a, b)))
                d = P.first_order_dressing(A, B, P.embed_slot(A, B, T, h[:, None]))
                assert close(np.linalg.norm(d) ** 2, P.single_slot_two_point(A, B, T, h), rel=1e-8, abs_=1e-12)
                checked += 1
    assert checked > 20


def test_tower_is_dark_at_first_order():
    P = pair(5, 80)
    for n in range(1, P.N + 1):
        t = P.tower_state(n, "new")
        H, _ = P.harmonic_embedding(n, n)
        s = H @ (H.conj().T @ t)
        assert np.linalg.norm(P.first_order_dressing(n, n, s[:, None])) < 1e-10


# ------------------------------------------------------------------------------------------------ closed forms (D8 sec. 9)
def test_hpl_dimension_formula():
    """dim B_new = dim H(E1, D), D = d1 + d2 with d2 = P_H Zc P_H, in every sector (both models)."""
    for N, generic in ((3, False), (4, False), (5, False), (4, True), (5, True)):
        P = pair(N, 200 + N, generic)
        for A in range(N + 2):
            for B in range(N + 2):
                if P.new_dim(A, B):
                    assert P.hpl_dim(A, B) == P.d_new(A, B), (N, generic, A, B)


def test_closed_form_uplift_reproduces_bps_states():
    """Psi = (1 + Kc)^{-1} h reconstructs every BPS state from its harmonic slot parts; the admissible sources
    (vanishing harmonic parts of the currents) span exactly B_new."""
    from two_flavor_uplift import rank as _rank
    for generic in (False, True):
        P = pair(4, 210, generic)
        for A in range(P.N + 2):
            for B in range(P.N + 2):
                if P.new_dim(A, B) == 0:
                    continue
                Hadm, Psi, R = P.cf_bps(A, B)
                X = P.new_B(A, B)
                d = X.shape[1] if X.size else 0
                assert Hadm.shape[1] == d
                if d == 0:
                    continue
                Y = P.cf_to_slots(A, B, X)
                Hh = P.cf_harmonic(A, B)
                assert np.linalg.norm(R @ (Hh @ (Hh.conj().T @ Y)) - Y) < 1e-10 * np.linalg.norm(Y)
                L = P.bps_constraints(A, B)
                code = P.cf_from_slots(A, B, Psi)
                if L is not None:
                    assert np.linalg.norm(L @ code) < 1e-10 * np.linalg.norm(code)
                assert _rank(np.hstack([X, code])) == d          # same span


def test_reduced_system_structure():
    """A, A' Hermitian PSD and equal to the projected self-energies; L' = -L^dag + C with the factorized leak;
    Herm(M) >= 1 - |C|/2 (so M is invertible when |C| < 2); C = 0 when the middle slots carry no cohomology."""
    for generic in (False, True):
        P = pair(5, 220, generic)
        for A in range(P.N + 2):
            for B in range(P.N + 2):
                if P.new_dim(A, B) == 0:
                    continue
                r = P.reduced_system(A, B)
                for X in (r["A"], r["Ap"]):
                    if X.size:
                        assert np.abs(X - X.conj().T).max() < 1e-10
                        assert np.linalg.eigvalsh((X + X.conj().T) / 2).min() > -1e-10
                if r["L"].size:
                    assert np.abs(r["Lp"] + r["L"].conj().T - r["C"]).max() < 1e-9
                M = r["M"]
                if M.size:
                    cn = np.linalg.norm(r["C"], 2) if r["C"].size else 0.0
                    assert np.linalg.eigvalsh((M + M.conj().T) / 2).min() >= 1 - cn / 2 - 1e-9
                    bmid = (P.b(A, B) if P.old_ok(A, B) else 0) + (P.b(A - 1, B - 1) if P.old_ok(A - 1, B - 1) else 0)
                    if bmid == 0:
                        assert cn < 1e-10


def test_closed_form_decoders_match_exact():
    for generic in (False, True):
        P = pair(4, 230, generic)
        for A in range(P.N + 2):
            for B in range(P.N + 2):
                r = P.cf_decoders(A, B) if P.new_dim(A, B) else None
                if r is None:
                    continue
                assert np.abs(r["mu"] - np.sort(P.bare_fraction_spectrum(A, B))).max() < 1e-9
                X = P.new_B(A, B)
                S = P.slot_components(A, B, X)
                for T in SLOTS:
                    key = "".join(map(str, T))
                    if key in r["lam"]:
                        ex = np.sort(np.linalg.eigvalsh(S[T].conj().T @ S[T]))
                        assert np.abs(r["lam"][key] - ex).max() < 1e-9



def test_dressing_levels_consistency():
    """The closed-form dressing quadratic form equals |Psi|^2 - |h|^2; self-dressing moves the mean bare fraction
    from first order towards the exact value (central sector, N = 5)."""
    P = pair(5, 240)
    r = P.cf_dressing_levels(3, 3)
    assert r["consistency"] < 1e-10
    assert np.abs(r["exact"] - np.sort(P.bare_fraction_spectrum(3, 3))).max() < 1e-9
    mf, ms, me = np.mean(r["first"]), np.mean(r["self"]), np.mean(r["exact"])
    assert mf < ms <= me + 1e-12


if __name__ == "__main__":
    import time
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t0 = time.time()
            fn()
            print(f"ok  {name}  [{time.time() - t0:.1f}s]", flush=True)
