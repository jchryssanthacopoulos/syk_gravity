r"""
two_flavor_uplift.py -- the exact BPS uplift N -> N+1 of the symmetrized two-flavor N=2 SYK model (docs/derivations.md
D8): the pair (psi_e, chi_f) is added as two single-mode steps, each a mapping cone with an obstruction.

Conventions (D8 sec. 1).  The new pair has coupling index N (0-based) and sits at bits 2N (psi_e) and 2N+1 (chi_f) of the
new model, so the old model is the restriction to the lower 2N bits.  The supercharge splits as
    Q_new = Q0 + psi_e^dag E + K chi_f + psi_e^dag Z chi_f ,
    E = sum_{i,k} (C_eik - C_iek) psi_i^dag chi_k     (even, charge (+1,-1)),
    K = sum_{i<j} (C_ijf - C_jif) psi_i^dag psi_j^dag (even, charge (+2, 0)),
    Z = sum_i (C_eif - C_ief) psi_i^dag               (odd,  charge (+1, 0)),   e = f = N (coupling index).
Slots T = (t_psi, t_chi) in the OPERATOR convention
    Psi = Psi00 + psi_e^dag Psi10 + chi_f^dag Psi01 + psi_e^dag chi_f^dag Psi11 .
The code basis |l> (x) new bits = c^dag_{l...} (psi_e^dag)^t_psi (chi_f^dag)^t_chi |0> differs from it by the sign
(-1)^{|l| (t_psi + t_chi)} (|l| = old fermion number), i.e. a fixed sign on slots 10 and 01 of each sector.
`slot_components(..., operator_convention=True)` applies it.  Ranks, fidelities, bare fractions and dressing costs are
insensitive to these signs; the slot equations (D8 eq. S) need them.

Old sectors (a, b) = (N_psi, N_chi); slot T of new sector (A, B) is old sector (A - t_psi, B - t_chi).
Dense linear algebra throughout: fine up to N' = N + 1 = 8 (largest sector 4900, peak ~2.5 GB).
"""
from __future__ import annotations

from math import comb

import numpy as np

from two_flavor import TwoFlavor, vacuum_tower

SLOTS = ((0, 0), (1, 0), (0, 1), (1, 1))
CHANNELS = ("00", "10", "01", "11")
RTOL = 1e-9


# ---------------------------------------------------------------------------------------------------- linear algebra
def colspace(M, rtol=RTOL):
    """Orthonormal basis of the column space of M (relative tolerance on singular values)."""
    if M.size == 0 or M.shape[1] == 0:
        return np.zeros((M.shape[0], 0), complex)
    U, s, _ = np.linalg.svd(M, full_matrices=False)
    if len(s) == 0 or s.max() == 0:
        return np.zeros((M.shape[0], 0), complex)
    return U[:, s > rtol * s.max()]


def rank(M, rtol=RTOL):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > rtol * max(1.0, s.max()))) if len(s) else 0


def nullspace(M, n, rtol=RTOL):
    """Orthonormal basis of ker M (M has n columns; M may have zero rows)."""
    if M is None or M.shape[0] == 0:
        return np.eye(n, dtype=complex)
    _, s, Vh = np.linalg.svd(M, full_matrices=True)
    r = int(np.sum(s > rtol * max(1.0, s.max()))) if len(s) else 0
    return Vh[r:].conj().T


def orth_complement(H, d):
    """Orthonormal basis of span(H)^perp in C^d."""
    if H.size == 0 or H.shape[1] == 0:
        return np.eye(d, dtype=complex)
    U, s, _ = np.linalg.svd(H, full_matrices=True)
    return U[:, int(np.sum(s > 1e-10)):]


def laplacian_kernel(Qout, Qin, d, tol=RTOL):
    """Harmonic subspace ker(Qout^dag Qout + Qin Qin^dag) of a sector of dimension d (Qout/Qin may be None)."""
    if d == 0:
        return np.zeros((0, 0), complex)
    L = np.zeros((d, d), complex)
    if Qout is not None and Qout.size:
        L += Qout.conj().T @ Qout
    if Qin is not None and Qin.size:
        L += Qin @ Qin.conj().T
    w, U = np.linalg.eigh((L + L.conj().T) / 2)
    return U[:, w < tol * max(1.0, w.max())]


# ---------------------------------------------------------------------------------------------------- the pair uplift
class PairUplift:
    """Old model (N pairs) and new model (N+1 pairs, new pair last) for couplings C1 of shape (N+1, N+1, N+1)."""

    def __init__(self, C1):
        C1 = np.asarray(C1)
        self.N = N = C1.shape[0] - 1
        self.C1 = C1
        self.old = TwoFlavor(C1[:N, :N, :N].copy())
        self.new = TwoFlavor(C1, psi_pos=list(range(N)) + [2 * N], chi_pos=list(range(N, 2 * N)) + [2 * N + 1])
        self.be, self.bf = 1 << (2 * N), 1 << (2 * N + 1)
        self._c = {}

    def _cache(self, key, fn):
        if key not in self._c:
            self._c[key] = fn()
        return self._c[key]

    def clear_cache(self, *kinds):
        """Drop cached objects whose key starts with one of `kinds` (e.g. 'MB', 'fo', 'nQ'); all if none given."""
        for k in list(self._c):
            if not kinds or k[0] in kinds:
                del self._c[k]

    # ------------------------------------------------------------------ old model
    def old_ok(self, a, b):
        return 0 <= a <= self.N and 0 <= b <= self.N

    def old_dim(self, a, b):
        return self.old.dim(a, b) if self.old_ok(a, b) else 0

    def old_Q(self, a, b):
        """Q0: (a, b) -> (a+2, b-1), dense."""
        def f():
            d0, d1 = self.old_dim(a, b), self.old_dim(a + 2, b - 1)
            return self.old.Q(a, b).toarray() if d0 and d1 else np.zeros((d1, d0), complex)
        return self._cache(("oQ", a, b), f)

    def old_B(self, a, b):
        """Orthonormal basis of the old harmonic (BPS) space in sector (a, b)."""
        def f():
            d = self.old_dim(a, b)
            if d == 0:
                return np.zeros((0, 0), complex)
            return laplacian_kernel(self.old_Q(a, b), self.old_Q(a - 2, b + 1), d)
        return self._cache(("oB", a, b), f)

    def b(self, a, b):
        B = self.old_B(a, b)
        return B.shape[1] if B.size else 0

    def rank_out(self, a, b):
        """rank of Q0 out of (a, b) = dim of the coexact forms in (a, b) = dim of the exact forms in (a+2, b-1)."""
        return self._cache(("rk", a, b), lambda: rank(self.old_Q(a, b)) if self.old_dim(a, b) else 0)

    def closed_dim(self, a, b):
        return self.old_dim(a, b) - self.rank_out(a, b) if self.old_dim(a, b) else 0

    def coclosed_dim(self, a, b):
        return self.old_dim(a, b) - self.rank_out(a - 2, b + 1) if self.old_dim(a, b) else 0

    def exact_basis(self, a, b):
        return self._cache(("ex", a, b), lambda: colspace(self.old_Q(a - 2, b + 1)))

    def coexact_basis(self, a, b):
        return self._cache(("co", a, b), lambda: colspace(self.old_Q(a, b).conj().T))

    def closed_basis(self, a, b):
        return self._cache(("cl", a, b), lambda: nullspace(self.old_Q(a, b), self.old_dim(a, b)))

    def coclosed_basis(self, a, b):
        def f():
            Qi = self.old_Q(a - 2, b + 1)
            return nullspace(Qi.conj().T, self.old_dim(a, b)) if Qi.size else np.eye(self.old_dim(a, b), dtype=complex)
        return self._cache(("ccl", a, b), f)

    # ------------------------------------------------------------------ pair operators (old Fock space, operator convention)
    def _pair_terms(self, which):
        N, C = self.N, self.C1
        ps, ch = self.old.psi_pos, self.old.chi_pos
        if which == "E":
            return [(C[N, i, k] - C[i, N, k], [(ps[i], True), (ch[k], False)]) for i in range(N) for k in range(N)]
        if which == "K":
            return [(C[i, j, N] - C[j, i, N], [(ps[i], True), (ps[j], True)]) for i in range(N) for j in range(i + 1, N)]
        if which == "Z":
            return [(C[N, i, N] - C[i, N, N], [(ps[i], True)]) for i in range(N)]
        raise ValueError(which)

    _SHIFT = {"E": (1, -1), "K": (2, 0), "Z": (1, 0)}

    def pair_op(self, which, a, b):
        """E: (a,b)->(a+1,b-1); K: (a,b)->(a+2,b); Z: (a,b)->(a+1,b); dense old-sector matrices."""
        def f():
            da, db = self._SHIFT[which]
            d0, d1 = self.old_dim(a, b), self.old_dim(a + da, b + db)
            if not (d0 and d1):
                return np.zeros((d1, d0), complex)
            return self.old.op_block(self._pair_terms(which), (a, b), (a + da, b + db)).toarray()
        return self._cache(("op", which, a, b), f)

    def E(self, a, b):
        return self.pair_op("E", a, b)

    def K(self, a, b):
        return self.pair_op("K", a, b)

    def Z(self, a, b):
        return self.pair_op("Z", a, b)

    # ------------------------------------------------------------------ new model
    def new_ok(self, A, B):
        return 0 <= A <= self.N + 1 and 0 <= B <= self.N + 1

    def new_dim(self, A, B):
        return self.new.dim(A, B) if self.new_ok(A, B) else 0

    def new_Q(self, A, B):
        def f():
            d0, d1 = self.new_dim(A, B), self.new_dim(A + 2, B - 1)
            return self.new.Q(A, B).toarray() if d0 and d1 else np.zeros((d1, d0), complex)
        return self._cache(("nQ", A, B), f)

    def new_B(self, A, B):
        def f():
            d = self.new_dim(A, B)
            if d == 0:
                return np.zeros((0, 0), complex)
            return laplacian_kernel(self.new_Q(A, B), self.new_Q(A - 2, B + 1), d)
        return self._cache(("nB", A, B), f)

    def d_new(self, A, B):
        X = self.new_B(A, B)
        return X.shape[1] if X.size else 0

    def bps_constraints(self, A, B):
        """[Q_out; Q_in^dag] on new sector (A, B) (None if both vanish)."""
        L = []
        if self.new_dim(A + 2, B - 1):
            L.append(self.new_Q(A, B))
        if self.new_dim(A - 2, B + 1):
            L.append(self.new_Q(A - 2, B + 1).conj().T)
        return np.vstack(L) if L else None

    # ------------------------------------------------------------------ slots
    def slot_index(self, A, B):
        """dict T -> integer array idx such that new-sector entries idx[j] <-> old-sector (A-t_psi, B-t_chi) basis j."""
        def f():
            masks = self.new.sector(A, B)
            low = masks & (self.be - 1)
            tp, tc = (masks & self.be) != 0, (masks & self.bf) != 0
            out = {}
            for T in SLOTS:
                a, b = A - T[0], B - T[1]
                if not self.old_ok(a, b):
                    out[T] = np.zeros(0, dtype=int)
                    continue
                sel = np.where((tp == bool(T[0])) & (tc == bool(T[1])))[0]
                om = self.old.sector(a, b)
                idx = np.empty(len(om), dtype=int)
                pos = np.searchsorted(om, low[sel])
                assert np.all(om[pos] == low[sel])
                idx[pos] = sel
                out[T] = idx
            return out
        return self._cache(("sidx", A, B), f)

    def slot_sign(self, A, B, T):
        """Code amplitude -> operator-convention amplitude sign for slot T of new sector (A, B)."""
        nold = (A - T[0]) + (B - T[1])
        return (-1) ** (nold * (T[0] + T[1]))

    def slot_components(self, A, B, X, operator_convention=False):
        """dict T -> old-sector components (rows = old basis of sector (A-t_psi, B-t_chi)) of the columns of X."""
        idx = self.slot_index(A, B)
        X = X.reshape(X.shape[0], -1)
        out = {}
        for T in SLOTS:
            Y = X[idx[T]] if len(idx[T]) else np.zeros((0, X.shape[1]), complex)
            out[T] = self.slot_sign(A, B, T) * Y if operator_convention else Y
        return out

    def embed_slot(self, A, B, T, Y, operator_convention=False):
        """Embed old-sector columns Y as slot-T components of new sector (A, B)."""
        out = np.zeros((self.new_dim(A, B), Y.shape[1]), complex)
        idx = self.slot_index(A, B)[T]
        if len(idx):
            out[idx] = (self.slot_sign(A, B, T) if operator_convention else 1) * Y
        return out

    def harmonic_embedding(self, A, B):
        """(H, labels): orthonormal columns spanning (+)_T old-harmonic(slot T) inside new sector (A, B); labels = slot of
        each column."""
        def f():
            cols, labs = [], []
            for t, T in enumerate(SLOTS):
                a, b = A - T[0], B - T[1]
                if not self.old_ok(a, b) or self.b(a, b) == 0:
                    continue
                cols.append(self.embed_slot(A, B, T, self.old_B(a, b)))
                labs += [t] * self.b(a, b)
            H = np.hstack(cols) if cols else np.zeros((self.new_dim(A, B), 0), complex)
            return H, np.array(labs, dtype=int)
        return self._cache(("hemb", A, B), f)

    def slot_labels(self, A, B):
        """Slot index (0..3 in SLOTS order) of every basis state of new sector (A, B)."""
        masks = self.new.sector(A, B)
        tp, tc = (masks & self.be) != 0, (masks & self.bf) != 0
        return np.array([SLOTS.index((int(p), int(c))) for p, c in zip(tp, tc)], dtype=int)

    # ------------------------------------------------------------------ Theorem A twice: the two-step dimension formula
    def _M_idx(self, A, B):
        """Intermediate model M (chi_f frozen empty): indices of new sector (A, B) with chi_f empty."""
        m = self.new.sector(A, B)
        return np.where((m & self.bf) == 0)[0]

    def _QM(self, A, B):
        """Q_M: M(A, B) -> M(A+2, B-1)."""
        Q = self.new_Q(A, B)
        r, c = self._M_idx(A + 2, B - 1), self._M_idx(A, B)
        return Q[np.ix_(r, c)] if Q.size else np.zeros((len(r), len(c)), complex)

    def M_B(self, A, B):
        return self._cache(("MB", A, B), lambda: laplacian_kernel(self._QM(A, B), self._QM(A - 2, B + 1),
                                                                  len(self._M_idx(A, B))))

    def _Kprime(self, A, B):
        """K' = K + psi_e^dag Z (chi_f-annihilating block of Q_new): M(A, B-1) -> M(A+2, B-1), in M bases."""
        Q = self.new_Q(A, B)
        src, dst = self.new.sector(A, B), self.new.sector(A + 2, B - 1)
        mS = self.new.sector(A, B - 1)[self._M_idx(A, B - 1)] if self.new_dim(A, B - 1) else np.zeros(0, np.int64)
        mT = dst[self._M_idx(A + 2, B - 1)] if len(dst) else np.zeros(0, np.int64)
        Kp = np.zeros((len(mT), len(mS)), complex)
        if Q.size:
            cs = np.where((src & self.bf) != 0)[0]
            rd = np.where((dst & self.bf) == 0)[0]
            if len(cs) and len(rd):
                Kp[np.ix_(np.searchsorted(mT, dst[rd]), np.searchsorted(mS, src[cs] & ~self.bf))] = Q[np.ix_(rd, cs)]
        return Kp

    def two_step_dimensions(self, A, B):
        """Predicted (Theorem A, twice) vs actual dim B_M(A,B) and dim B_new(A,B).
        Step 1 (psi_e, creation type): dim B_M(A,B) = dim ker(P E P | B(A,B)) + dim ker(P E^dag P | B(A-1,B)).
        Step 2 (chi_f, annihilation type): dim B_new(A,B) = dim ker(P K'^dag P | B_M(A,B)) + dim ker(P K' P | B_M(A,B-1))."""
        def kdim(M, n):
            return n - rank(M) if n else 0
        h, g = self.old_B(A, B), self.old_B(A - 1, B)
        k1 = 0
        if h.size and h.shape[1]:
            t = self.old_B(A + 1, B - 1)
            k1 += kdim(t.conj().T @ self.E(A, B) @ h, h.shape[1]) if t.size and t.shape[1] else h.shape[1]
        if g.size and g.shape[1]:
            t = self.old_B(A - 2, B + 1)
            k1 += kdim(t.conj().T @ self.E(A - 2, B + 1).conj().T @ g, g.shape[1]) if t.size and t.shape[1] else g.shape[1]
        hM, gM = self.M_B(A, B), self.M_B(A, B - 1)
        k2 = 0
        if hM.size and hM.shape[1]:
            t = self.M_B(A - 2, B)
            k2 += kdim(t.conj().T @ self._Kprime(A - 2, B + 1).conj().T @ hM, hM.shape[1]) if t.size and t.shape[1] \
                else hM.shape[1]
        if gM.size and gM.shape[1]:
            t = self.M_B(A + 2, B - 1)
            k2 += kdim(t.conj().T @ self._Kprime(A, B) @ gM, gM.shape[1]) if t.size and t.shape[1] else gM.shape[1]
        dM = self.M_B(A, B).shape[1] if self.M_B(A, B).size else 0
        return dict(M_pred=k1, M=dM, new_pred=k2, new=self.d_new(A, B))

    # ------------------------------------------------------------------ the slot complex (E2 page of the w-filtration)
    def _hblock(self, which, a, b, ta, tb):
        """P_H(target) op P_H(source) in harmonic coordinates (old (a,b) -> old (ta,tb))."""
        if not (self.old_ok(a, b) and self.old_ok(ta, tb)) or self.b(a, b) == 0 or self.b(ta, tb) == 0:
            return np.zeros((self.b(ta, tb) if self.old_ok(ta, tb) else 0, self.b(a, b) if self.old_ok(a, b) else 0), complex)
        return self.old_B(ta, tb).conj().T @ self.pair_op(which, a, b) @ self.old_B(a, b)

    def slot_complex_E2(self, A, B):
        """dims (w=-1, w=0, w=+1) of the cohomology of the slot complex at new sector (A, B):
            H(01) --d1--> H(00) (+) H(11) --d1--> H(10),   d1 = ([K], [E]) and (h00, h11) -> [E] h00 - [K] h11,
        built from old BPS states (E1 = old cohomology of each slot; w = t_psi - t_chi).  dim B_new = sum iff the
        spectral sequence degenerates at E2; a nonzero d2 ~ [Z - E Q0^+ K + K Q0^+ E] removes classes (D8 sec. 3)."""
        def bb(a, b):
            return self.b(a, b) if self.old_ok(a, b) else 0

        def rk(blocks, axis):
            blocks = [x for x in blocks if x.size]
            if not blocks:
                return 0
            return rank(np.vstack(blocks) if axis == 0 else np.hstack(blocks))
        n_m1 = bb(A, B - 1)
        e_m1 = n_m1 - (rk([self._hblock("K", A, B - 1, A + 2, B - 1), self._hblock("E", A, B - 1, A + 1, B - 2)], 0)
                       if n_m1 else 0)
        n0 = bb(A, B) + bb(A - 1, B - 1)
        out0 = rk([self._hblock("E", A, B, A + 1, B - 1), -self._hblock("K", A - 1, B - 1, A + 1, B - 1)], 1) \
            if (n0 and bb(A + 1, B - 1)) else 0
        in0 = rk([self._hblock("K", A - 2, B, A, B), self._hblock("E", A - 2, B, A - 1, B - 1)], 0) \
            if (n0 and bb(A - 2, B)) else 0
        n1 = bb(A - 1, B)
        in1 = rk([self._hblock("E", A - 2, B + 1, A - 1, B), -self._hblock("K", A - 3, B, A - 1, B)], 1) \
            if (n1 and bb(A - 2, B + 1) + bb(A - 3, B)) else 0
        return e_m1, n0 - out0 - in0, n1 - in1

    # ------------------------------------------------------------------ four-slot harmonic map
    def injectivity_deficit(self, A, B):
        """dim of {BPS states of new (A,B) with vanishing old-harmonic part in every slot} (0 = injective)."""
        d = self.new_dim(A, B)
        if d == 0:
            return 0
        H, _ = self.harmonic_embedding(A, B)
        Nm = orth_complement(H, d)
        L = self.bps_constraints(A, B)
        if L is None:
            return Nm.shape[1]
        return Nm.shape[1] - rank(L @ Nm)

    def bare_fraction_spectrum(self, A, B, return_vectors=False):
        """Eigenvalues (ascending) of G = sum_T S_T^dag P_H(T) S_T on B_new(A, B): the fraction of each BPS state's norm
        carried by old-harmonic slot components.  Eigenvalue 1 = bare (all slots harmonic)."""
        X = self.new_B(A, B)
        if X.size == 0 or X.shape[1] == 0:
            return (np.zeros(0), X) if return_vectors else np.zeros(0)
        H, _ = self.harmonic_embedding(A, B)
        Cm = H.conj().T @ X
        w, U = np.linalg.eigh(Cm.conj().T @ Cm)
        return (w, X @ U) if return_vectors else w

    # ------------------------------------------------------------------ bare states: the bare map and its generic count
    def bare_count(self, A, B):
        """dim ker of the bare map beta = [Q_out; Q_in^dag] restricted to (+)_T old-harmonic(slot T)."""
        H, _ = self.harmonic_embedding(A, B)
        if H.shape[1] == 0:
            return 0
        L = self.bps_constraints(A, B)
        return H.shape[1] - (rank(L @ H) if L is not None else 0)

    def bare_count_generic(self, A, B, exact_subkernels=False):
        """Generic-position count of bare states from the slot equations (D8 sec. 5):
          h01: K h01 = 0 in closed(A+2,B-1), E h01 = 0 in closed(A+1,B-2)
          h10: E^dag h10 = 0 in coclosed(A-2,B+1), K^dag h10 = 0 in coclosed(A-3,B)
          J10 = E h00 - K h11 + Z h01 = 0 in closed(A+1,B-1);  J01 = K^dag h00 + E^dag h11 + Z^dag h10 = 0 in coclosed(A-2,B)
        bare_gen = max(0, b00 + b11 + n01 + n10 - c(J10) - c(J01)), n01 = max(0, b01 - c1 - c2), n10 likewise.
        exact_subkernels: use the exact dims of ker[K;E] on B(01) and ker[E^dag;K^dag] on B(10) instead of n01, n10."""
        bs = {T: (self.b(A - T[0], B - T[1]) if self.old_ok(A - T[0], B - T[1]) else 0) for T in SLOTS}
        c1, c2 = self.closed_dim(A + 2, B - 1), self.closed_dim(A + 1, B - 2)
        c3, c4 = self.coclosed_dim(A - 2, B + 1), self.coclosed_dim(A - 3, B)
        c5, c6 = self.closed_dim(A + 1, B - 1), self.coclosed_dim(A - 2, B)
        if exact_subkernels:
            n01 = self._subkernel(A, B - 1, "EK")
            n10 = self._subkernel(A - 1, B, "EdKd")
        else:
            n01 = max(0, bs[(0, 1)] - c1 - c2)
            n10 = max(0, bs[(1, 0)] - c3 - c4)
        gen = max(0, bs[(0, 0)] + bs[(1, 1)] + n01 + n10 - c5 - c6)
        return dict(gen=gen, n01=n01, n10=n10, slots={"".join(map(str, T)): v for T, v in bs.items()},
                    rooms=[c1, c2, c3, c4, c5, c6])

    def _subkernel(self, a, b, kind):
        """dim ker [K; E] (kind 'EK') or [E^dag; K^dag] (kind 'EdKd') on the old harmonic space B(a, b)."""
        if not self.old_ok(a, b) or self.b(a, b) == 0:
            return 0
        H = self.old_B(a, b)
        rows = []
        if kind == "EK":
            if self.old_dim(a + 2, b):
                rows.append(self.K(a, b) @ H)
            if self.old_dim(a + 1, b - 1):
                rows.append(self.E(a, b) @ H)
        else:
            if self.old_dim(a - 1, b + 1):
                rows.append(self.E(a - 1, b + 1).conj().T @ H)
            if self.old_dim(a - 2, b):
                rows.append(self.K(a - 2, b).conj().T @ H)
        return H.shape[1] - (rank(np.vstack(rows)) if rows else 0)

    # ------------------------------------------------------------------ CV decoders: rooms, exact-lift loci, fidelities
    def rooms(self, a, b):
        """Constrained rooms R_T in old sector (a, b) containing im D_T (D8 sec. 6):
          01: ker Q0;  10: ker Q0^dag;
          00: harm + Q0^+ P_ex K(closed(a,b-1)) + (Q0^dag)^+ P_co E^dag(coclosed(a-1,b));
          11: harm + Q0^+ P_ex E(closed(a+1,b)) + (Q0^dag)^+ P_co K^dag(coclosed(a,b+1))."""
        bh = self.b(a, b)
        out = {"01": self.closed_dim(a, b), "10": self.coclosed_dim(a, b), "full": self.old_dim(a, b)}

        def part(src_ok, Pbasis, op, srcbasis):
            if not src_ok or Pbasis.shape[1] == 0:
                return 0
            return colspace(Pbasis.conj().T @ op @ srcbasis).shape[1]

        r00 = 0
        if self.old_dim(a, b - 1) and self.old_dim(a + 2, b - 1):
            r00 += part(True, self.exact_basis(a + 2, b - 1), self.K(a, b - 1), self.closed_basis(a, b - 1))
        if self.old_dim(a - 1, b) and self.old_dim(a - 2, b + 1):
            r00 += part(True, self.coexact_basis(a - 2, b + 1), self.E(a - 2, b + 1).conj().T, self.coclosed_basis(a - 1, b))
        r11 = 0
        if self.old_dim(a + 1, b) and self.old_dim(a + 2, b - 1):
            r11 += part(True, self.exact_basis(a + 2, b - 1), self.E(a + 1, b), self.closed_basis(a + 1, b))
        if self.old_dim(a, b + 1) and self.old_dim(a - 2, b + 1):
            r11 += part(True, self.coexact_basis(a - 2, b + 1), self.K(a - 2, b + 1).conj().T, self.coclosed_basis(a, b + 1))
        out["00"], out["11"] = bh + r00, bh + r11
        return out

    def tower_state(self, n, model="old"):
        m = self.old if model == "old" else self.new
        t = vacuum_tower(m, n)
        return t / np.linalg.norm(t)

    def cv_channels(self, a, b, symmetrized=True):
        """For old sector (a, b): per channel T (target new sector (a+t_psi, b+t_chi)) the decoder rank, room, exact-lift
        locus dimension dim(B cap im D_T), its generic-position prediction, and CV's U_T (eq. 89) on the fortuitous
        subspace (B minus the tower when a == b in the symmetrized model) and on the tower."""
        Bo = self.old_B(a, b)
        if Bo.size == 0 or Bo.shape[1] == 0:
            return None
        S, tow = Bo, None
        if symmetrized and a == b:
            tow = self.tower_state(a)
            c = Bo.conj().T @ tow
            S = colspace(Bo @ (np.eye(Bo.shape[1]) - np.outer(c, c.conj())))
        rm = self.rooms(a, b)
        rec = dict(a=a, b=b, b_old=Bo.shape[1], fort=S.shape[1], dim=self.old_dim(a, b), room=rm)
        for T, key in zip(SLOTS, CHANNELS):
            A, Bn = a + T[0], b + T[1]
            X = self.new_B(A, Bn)
            dn = X.shape[1] if X.size else 0
            if dn == 0:
                Im = np.zeros((self.old_dim(a, b), 0), complex)
            else:
                Im = colspace(self.slot_components(A, Bn, X)[T])
            rk = Im.shape[1]
            sv = np.linalg.svd(Bo.conj().T @ Im, compute_uv=False) if rk else np.zeros(0)
            dimA = int(np.sum(sv > 1 - 1e-8))
            Uf = Ut = None
            if S.shape[1]:
                Mx = S.conj().T @ Im
                Uf = float(max(0.0, np.linalg.eigvalsh(Mx @ Mx.conj().T).min())) if rk else 0.0
            if tow is not None:
                Ut = float(np.linalg.norm(Im.conj().T @ tow) ** 2) if rk else 0.0
            rec[key] = dict(d_new=dn, rkD=rk, R=rm[key], dimA=dimA, a_gen=max(0, Bo.shape[1] + rk - rm[key]),
                            rk_pred=min(dn, rm[key]), U_fort=Uf, U_tow=Ut)
        return rec

    # ------------------------------------------------------------------ the tower: CV dressing costs (eqs. 87-88)
    def tower_cost(self, n, channel):
        """kappa_T of the old tower state V^n|0> (T = '00' -> new (n,n); '11' -> new (n+1,n+1)) as an exact min-norm
        uplift, CV's bound, and dim ker D_T (BPS states of the target sector with vanishing slot-T component)."""
        T = (0, 0) if channel == "00" else (1, 1)
        A, B = n + T[0], n + T[1]
        alpha = self.tower_state(n)
        a = self.embed_slot(A, B, T, alpha[:, None])[:, 0]
        rest = np.setdiff1d(np.arange(self.new_dim(A, B)), self.slot_index(A, B)[T])
        L = self.bps_constraints(A, B)
        if L is None:
            Psi, kerD = a, len(rest)
        else:
            y, *_ = np.linalg.lstsq(L[:, rest], -L @ a, rcond=None)
            Psi = a.copy()
            Psi[rest] += y
            kerD = len(rest) - rank(L[:, rest])
            assert np.linalg.norm(L @ Psi) < 1e-8, "tower has no exact uplift?"
        bound = n / (self.N + 1 - n) if channel == "00" else (self.N - n) / (n + 1)
        return dict(n=n, channel=channel, kappa=float(np.linalg.norm(Psi) ** 2 - 1), bound=bound, kerD=int(kerD))

    def tower_slot_weights(self, n):
        """Slot weights of the new tower V_{N+1}^n|0> (normalized) and its bare fraction (should be 1)."""
        t = self.tower_state(n, "new")
        S = self.slot_components(n, n, t[:, None])
        w = {"".join(map(str, T)): float(np.linalg.norm(S[T]) ** 2) for T in SLOTS}
        H, _ = self.harmonic_embedding(n, n)
        return w, float(np.linalg.norm(H.conj().T @ t) ** 2)

    # ------------------------------------------------------------------ first-order dressing
    def _first_order_ops(self, A, B):
        """Slot-diagonal pseudo-inverse and cross-slot parts of Q_out and Q_in^dag on new sector (A, B)."""
        lab = self.slot_labels(A, B)
        ops = []
        if self.new_dim(A + 2, B - 1):
            Q = self.new_Q(A, B)
            same = self.slot_labels(A + 2, B - 1)[:, None] == lab[None, :]
            ops.append((np.linalg.pinv(np.where(same, Q, 0), rcond=1e-10), np.where(same, 0, Q)))
        if self.new_dim(A - 2, B + 1):
            Qd = self.new_Q(A - 2, B + 1).conj().T                       # new (A,B) -> new (A-2,B+1)
            same = self.slot_labels(A - 2, B + 1)[:, None] == lab[None, :]
            ops.append((np.linalg.pinv(np.where(same, Qd, 0), rcond=1e-10), np.where(same, 0, Qd)))
        return ops, lab

    def first_order_dressing(self, A, B, S):
        """delta1(S) = -sum_ops Qd^+ Qc S for harmonic source columns S (new-sector vectors).  Its slot components are the
        Green's-operator images of the slot currents (D8 sec. 7), e.g. z10 = Q0^+ (E h00 - K h11 + Z h01)."""
        ops, _ = self._cache(("fo", A, B), lambda: self._first_order_ops(A, B))
        out = np.zeros_like(S)
        for Pi, Qc in ops:
            out -= Pi @ (Qc @ S)
        return out

    def first_order(self, A, B):
        """Exact vs first-order bare-fraction spectrum on B_new(A, B), and the per-state first-order dressing split into
        the incoherent (single-slot) sum and the interference term, in the eigenbasis of the exact bare-fraction decoder."""
        X = self.new_B(A, B)
        if X.size == 0 or X.shape[1] == 0:
            return None
        H, _ = self.harmonic_embedding(A, B)
        lab = self.slot_labels(A, B)
        Cm = H.conj().T @ X
        mu, U = np.linalg.eigh(Cm.conj().T @ Cm)
        S0 = H @ (Cm @ U)                                   # admissible sources of the eigenvectors
        n2 = np.linalg.norm(S0, axis=0) ** 2
        D = self.first_order_dressing(A, B, S0)
        k1 = np.linalg.norm(D, axis=0) ** 2 / n2
        kinc = np.zeros(len(mu))
        for t in range(4):
            kinc += np.linalg.norm(self.first_order_dressing(A, B, np.where((lab == t)[:, None], S0, 0)), axis=0) ** 2
        kinc /= n2
        # first-order spectrum: generalized eigenvalues of (G0, G0 + delta^dag delta) on the source space
        S = H @ Cm
        D = self.first_order_dressing(A, B, S)
        G0 = S.conj().T @ S
        G1 = G0 + D.conj().T @ D
        w0, U0 = np.linalg.eigh((G0 + G0.conj().T) / 2)
        keep = w0 > 1e-10 * w0.max()
        Tm = U0[:, keep] / np.sqrt(w0[keep])
        Mm = Tm.conj().T @ G1 @ Tm
        mu1 = np.sort(1 / np.linalg.eigvalsh((Mm + Mm.conj().T) / 2))
        return dict(mu_exact=mu, mu_first=mu1, kappa1=k1, kappa_inc=kinc)

    def single_slot_two_point(self, A, B, T, h):
        """First-order dressing norm^2 of a single-slot source h (old-sector vector in slot T) from the time-integrated
        two-point functions of D8 sec. 7 (Green's operator G = H0^+ of the old model):
          00: <E h, G E h> + <K^dag h, G K^dag h>          11: <K h, G K h> + <E^dag h, G E^dag h>
          01: <K h, G K h> + <E h, G E h> + <Z h, G Z h>   10: <E^dag h, G E^dag h> + <K^dag h, G K^dag h> + <Z^dag h, G Z^dag h>
        where each current is projected on its (co)exact part, as the pseudo-inverse does."""
        a, b = A - T[0], B - T[1]
        terms = []
        if T == (0, 0):
            terms = [("E", a, b, False), ("K", a - 2, b, True)]
        elif T == (1, 1):
            terms = [("K", a, b, False), ("E", a - 1, b + 1, True)]
        elif T == (0, 1):
            terms = [("K", a, b, False), ("E", a, b, False), ("Z", a, b, False)]
        elif T == (1, 0):
            terms = [("E", a - 1, b + 1, True), ("K", a - 2, b, True), ("Z", a - 1, b, True)]
        tot = 0.0
        for which, sa, sb, dag in terms:
            if dag:                                         # op^dag maps (a,b) -> source sector (sa, sb) of op
                if not self.old_ok(sa, sb) or self.old_dim(sa, sb) == 0:
                    continue
                j = self.pair_op(which, sa, sb).conj().T @ h
                tgt = (sa, sb)
            else:
                da, db = self._SHIFT[which]
                tgt = (a + da, b + db)
                if not self.old_ok(*tgt) or self.old_dim(*tgt) == 0:
                    continue
                j = self.pair_op(which, a, b) @ h
            tot += self._green_norm(tgt, j)
        return tot

    def _green_norm(self, sector, j):
        """||Q0^+ P_ex j||^2 + ||(Q0^dag)^+ P_co j||^2 = <j, H0^+ j> restricted to the non-harmonic part of sector."""
        a, b = sector
        out = 0.0
        Qin = self.old_Q(a - 2, b + 1)                      # exact part: Q0 from (a-2,b+1)
        if Qin.size:
            x = np.linalg.pinv(Qin, rcond=1e-10) @ j
            out += float(np.linalg.norm(x) ** 2)
        Qout = self.old_Q(a, b)                             # coexact part: Q0^dag from (a+2,b-1)
        if Qout.size:
            y = np.linalg.pinv(Qout.conj().T, rcond=1e-10) @ j
            out += float(np.linalg.norm(y) ** 2)
        return out

    # ================================================================== closed forms (D8 sec. 9)
    # Slot order by w = t_psi - t_chi:  01, 00, 11, 10.  New sector (A,B) -> old sectors (A,B-1), (A,B), (A-1,B-1), (A-1,B).
    # "Slot vectors" are the concatenated operator-convention slot components in this order.
    CF_SLOTS = ((0, 1), (0, 0), (1, 1), (1, 0))

    def _d(self, a, b):
        return self.old_dim(a, b) if self.old_ok(a, b) else 0

    def _Q(self, a, b):
        """Q0: (a,b) -> (a+2,b-1), with correct (possibly empty) shape."""
        return self.old_Q(a, b) if self._d(a, b) else np.zeros((self._d(a + 2, b - 1), 0), complex)

    def Qp(self, a, b):
        """Q0^+ : (a+2,b-1) -> (a,b)."""
        def f():
            Q = self._Q(a, b)
            return np.linalg.pinv(Q, rcond=1e-10) if Q.size else np.zeros((self._d(a, b), self._d(a + 2, b - 1)), complex)
        return self._cache(("Qp", a, b), f)

    def Qdp(self, a, b):
        """(Q0^dag)^+ : (a-2,b+1) -> (a,b)."""
        def f():
            Q = self._Q(a - 2, b + 1)
            return np.linalg.pinv(Q.conj().T, rcond=1e-10) if Q.size else np.zeros((self._d(a, b), self._d(a - 2, b + 1)), complex)
        return self._cache(("Qdp", a, b), f)

    def op(self, which, a, b):
        """E, K, Z from old sector (a,b) (shape-safe)."""
        da, db = self._SHIFT[which]
        if self._d(a, b) and self._d(a + da, b + db):
            return self.pair_op(which, a, b)
        return np.zeros((self._d(a + da, b + db), self._d(a, b)), complex)

    def opd(self, which, a, b):
        """Adjoint of E, K or Z whose TARGET is old sector (a,b): maps (a,b) back to the source sector."""
        da, db = self._SHIFT[which]
        return self.op(which, a - da, b - db).conj().T

    def H(self, a, b):
        B = self.old_B(a, b) if self._d(a, b) else np.zeros((0, 0))
        return B if B.size else np.zeros((self._d(a, b), 0), complex)

    def green(self, a, b):
        """Old Green's operator G = H0^+ on sector (a,b)."""
        def f():
            d = self._d(a, b)
            if d == 0:
                return np.zeros((0, 0), complex)
            L = np.zeros((d, d), complex)
            Qo, Qi = self._Q(a, b), self._Q(a - 2, b + 1)
            if Qo.size:
                L += Qo.conj().T @ Qo
            if Qi.size:
                L += Qi @ Qi.conj().T
            return np.linalg.pinv((L + L.conj().T) / 2, rcond=1e-10, hermitian=True)
        return self._cache(("G", a, b), f)

    def Zc(self, A, B):
        """Zc = Z - E Q0^+ K + K Q0^+ E : old (A,B-1) -> (A+1,B-1)  (the d2 kernel; D8 eq. Zc)."""
        return (self.op("Z", A, B - 1) - self.op("E", A, B) @ self.Qp(A, B) @ self.op("K", A, B - 1)
                + self.op("K", A - 1, B - 1) @ self.Qp(A - 1, B - 1) @ self.op("E", A, B - 1))

    # ------------------------------------------------------------------ (I) cohomology by homological perturbation
    def hpl_D(self, A, B):
        """Transferred differential D = d1 + d2 : E1(A,B) -> E1(A+2,B-1), E1 = (+)_T old BPS(slot T), slot order 01,00,11,10:
            D(eta01, eta00, eta11, eta10) = (0, P K eta01, P E eta01, P(E eta00 - K eta11 + Zc eta01))."""
        src = [self.H(A, B - 1), self.H(A, B), self.H(A - 1, B - 1), self.H(A - 1, B)]
        tgt = [self.H(A + 2, B - 2), self.H(A + 2, B - 1), self.H(A + 1, B - 2), self.H(A + 1, B - 1)]
        ns, nt = [x.shape[1] for x in src], [x.shape[1] for x in tgt]
        os_, ot = np.concatenate([[0], np.cumsum(ns)]), np.concatenate([[0], np.cumsum(nt)])
        M = np.zeros((ot[-1], os_[-1]), complex)

        def put(r, c, X):
            if nt[r] and ns[c]:
                M[ot[r]:ot[r + 1], os_[c]:os_[c + 1]] = tgt[r].conj().T @ X @ src[c]
        put(1, 0, self.op("K", A, B - 1))
        put(2, 0, self.op("E", A, B - 1))
        put(3, 0, self.Zc(A, B))
        put(3, 1, self.op("E", A, B))
        put(3, 2, -self.op("K", A - 1, B - 1))
        return M

    def hpl_dim(self, A, B):
        """dim H(E1, D) at new sector (A,B) = dim B_new(A,B) (D8 Theorem B)."""
        Dout, Din = self.hpl_D(A, B), self.hpl_D(A - 2, B + 1)
        return Dout.shape[1] - (rank(Dout) if Dout.size else 0) - (rank(Din) if Din.size else 0)

    # ------------------------------------------------------------------ (II) explicit two-sided uplift
    def cf_layout(self, A, B):
        secs = [(A - T[0], B - T[1]) for T in self.CF_SLOTS]
        dims = [self._d(*s) for s in secs]
        return secs, dims, np.concatenate([[0], np.cumsum(dims)]).astype(int)

    def cf_K(self, A, B):
        """Kc = Qd^+ Qc + (Qd^dag)^+ Qc^dag on slot vectors (01,00,11,10), operator convention (D8 eq. Kc)."""
        secs, dims, off = self.cf_layout(A, B)
        M = np.zeros((off[-1], off[-1]), complex)

        def put(r, c, X):
            if X.size and dims[r] and dims[c]:
                M[off[r]:off[r + 1], off[c]:off[c + 1]] += X
        put(0, 1, -self.Qdp(A, B - 1) @ self.opd("K", A, B))
        put(0, 2, -self.Qdp(A, B - 1) @ self.opd("E", A - 1, B - 1))
        put(0, 3, -self.Qdp(A, B - 1) @ self.opd("Z", A - 1, B))
        put(1, 0, self.Qp(A, B) @ self.op("K", A, B - 1))
        put(1, 3, self.Qdp(A, B) @ self.opd("E", A - 1, B))
        put(2, 0, self.Qp(A - 1, B - 1) @ self.op("E", A, B - 1))
        put(2, 3, -self.Qdp(A - 1, B - 1) @ self.opd("K", A - 1, B))
        put(3, 0, -self.Qp(A - 1, B) @ self.op("Z", A, B - 1))
        put(3, 1, -self.Qp(A - 1, B) @ self.op("E", A, B))
        put(3, 2, self.Qp(A - 1, B) @ self.op("K", A - 1, B - 1))
        return M

    def cf_harmonic(self, A, B):
        """Columns: the old-harmonic basis of each slot, embedded in slot-vector coordinates."""
        secs, dims, off = self.cf_layout(A, B)
        cols = []
        for i, s in enumerate(secs):
            Hs = self.H(*s)
            if Hs.shape[1]:
                Y = np.zeros((off[-1], Hs.shape[1]), complex)
                Y[off[i]:off[i + 1]] = Hs
                cols.append(Y)
        return np.hstack(cols) if cols else np.zeros((off[-1], 0), complex)

    def cf_to_slots(self, A, B, X):
        S = self.slot_components(A, B, X, operator_convention=True)
        return np.vstack([S[T] if S[T].size else np.zeros((0, X.shape[1]), complex) for T in self.CF_SLOTS])

    def cf_from_slots(self, A, B, Y):
        secs, dims, off = self.cf_layout(A, B)
        out = np.zeros((self.new_dim(A, B), Y.shape[1]), complex)
        idx = self.slot_index(A, B)
        for i, T in enumerate(self.CF_SLOTS):
            if dims[i]:
                out[idx[T]] = self.slot_sign(A, B, T) * Y[off[i]:off[i + 1]]
        return out

    def cf_obstruction(self, A, B):
        """Rows of Psi -> (P_H Qc Psi, P_H Qc^dag Psi): the harmonic parts of the six currents (D8 sec. 9.2)."""
        secs, dims, off = self.cf_layout(A, B)
        rows = []

        def row(parts, tgt):
            Ht = self.H(*tgt)
            if Ht.shape[1] == 0:
                return
            R = np.zeros((Ht.shape[1], off[-1]), complex)
            for i, X in parts:
                if X.size and dims[i]:
                    R[:, off[i]:off[i + 1]] += Ht.conj().T @ X
            rows.append(R)
        row([(0, self.op("K", A, B - 1))], (A + 2, B - 1))
        row([(0, self.op("E", A, B - 1))], (A + 1, B - 2))
        row([(0, self.op("Z", A, B - 1)), (1, self.op("E", A, B)), (2, -self.op("K", A - 1, B - 1))], (A + 1, B - 1))
        row([(1, self.opd("K", A, B)), (2, self.opd("E", A - 1, B - 1)), (3, self.opd("Z", A - 1, B))], (A - 2, B))
        row([(3, self.opd("E", A - 1, B))], (A - 2, B + 1))
        row([(3, -self.opd("K", A - 1, B))], (A - 3, B))
        return np.vstack(rows) if rows else np.zeros((0, off[-1]), complex)

    def cf_bps(self, A, B):
        """Closed-form BPS space: admissible harmonic sources h (columns, slot coordinates) and Psi = (1 + Kc)^{-1} h.
        Returns (Hadm, Psi, R) with R = (1 + Kc)^{-1}."""
        Kc = self.cf_K(A, B)
        R = np.linalg.inv(np.eye(Kc.shape[0]) + Kc)
        Hh = self.cf_harmonic(A, B)
        O = self.cf_obstruction(A, B)
        if Hh.shape[1] == 0:
            return Hh, Hh, R
        adm = nullspace(O @ R @ Hh, Hh.shape[1]) if O.shape[0] else np.eye(Hh.shape[1])
        Hadm = Hh @ adm
        return Hadm, R @ Hadm, R

    # ------------------------------------------------------------------ (III) the reduced two-slot system
    def self_energies(self, A, B):
        """S_G = E G E^dag + K G K^dag on old (A-1,B) [slot 10];  S'_G = K^dag G K + E^dag G E on old (A,B-1) [slot 01]."""
        SG = (self.op("E", A - 2, B + 1) @ self.green(A - 2, B + 1) @ self.op("E", A - 2, B + 1).conj().T
              + self.op("K", A - 3, B) @ self.green(A - 3, B) @ self.op("K", A - 3, B).conj().T)
        SpG = (self.op("K", A, B - 1).conj().T @ self.green(A + 2, B - 1) @ self.op("K", A, B - 1)
               + self.op("E", A, B - 1).conj().T @ self.green(A + 1, B - 2) @ self.op("E", A, B - 1))
        return SG, SpG

    def reduced_system(self, A, B):
        """The 2x2 block system for (z10, x01) in orthonormal bases of Co(10) and Ex(01) (D8 sec. 9.3):
             M = [[1 + A, -L], [-L', 1 + A']],   A = P S_G P,  A' = P S'_G P,  L = Q0^+ Zc,  L' = (Q0^dag)^+ Zc'^dag,
             L' = -L^dag + C,   C = alpha00 gamma00^dag - alpha11 gamma11^dag   (leak through the middle cohomology)."""
        Uz = colspace(self._Q(A - 1, B).conj().T) if (self._d(A - 1, B) and self._d(A + 1, B - 1)) \
            else np.zeros((self._d(A - 1, B), 0), complex)
        Ux = colspace(self._Q(A - 2, B)) if (self._d(A - 2, B) and self._d(A, B - 1)) \
            else np.zeros((self._d(A, B - 1), 0), complex)
        SG, SpG = self.self_energies(A, B)
        Am = Uz.conj().T @ SG @ Uz if Uz.shape[1] else np.zeros((0, 0), complex)
        Apm = Ux.conj().T @ SpG @ Ux if Ux.shape[1] else np.zeros((0, 0), complex)
        L = Uz.conj().T @ self.Qp(A - 1, B) @ self.Zc(A, B) @ Ux
        Lp = Ux.conj().T @ self.Qdp(A, B - 1) @ self.Zc(A - 2, B + 1).conj().T @ Uz
        nz, nx = Uz.shape[1], Ux.shape[1]
        C = np.zeros((nx, nz), complex)
        H00, H11 = self.H(A, B), self.H(A - 1, B - 1)
        if nz and nx and H00.shape[1]:
            a00 = Ux.conj().T @ self.Qdp(A, B - 1) @ self.opd("K", A, B) @ H00
            g00 = Uz.conj().T @ self.Qp(A - 1, B) @ self.op("E", A, B) @ H00
            C += a00 @ g00.conj().T
        if nz and nx and H11.shape[1]:
            a11 = Ux.conj().T @ self.Qdp(A, B - 1) @ self.opd("E", A - 1, B - 1) @ H11
            g11 = Uz.conj().T @ self.Qp(A - 1, B) @ self.op("K", A - 1, B - 1) @ H11
            C -= a11 @ g11.conj().T
        M = np.eye(nz + nx, dtype=complex)
        if nz:
            M[:nz, :nz] += Am
        if nx:
            M[nz:, nz:] += Apm
        if nz and nx:
            M[:nz, nz:] -= L
            M[nz:, :nz] -= Lp
        return dict(A=Am, Ap=Apm, L=L, Lp=Lp, C=C, M=M, nz=nz, nx=nx, Uz=Uz, Ux=Ux)

    def reduced_diagnostics(self, A, B):
        r = self.reduced_system(A, B)
        M = r["M"]
        if M.size == 0:
            return None
        herm = np.linalg.eigvalsh((M + M.conj().T) / 2).min()
        cn = float(np.linalg.norm(r["C"], 2)) if r["C"].size else 0.0
        out = dict(nz=r["nz"], nx=r["nx"], sigma_min_M=float(np.linalg.svd(M, compute_uv=False).min()),
                   herm_min=float(herm), leak_norm=cn, accretive_bound=1 - cn / 2,
                   b_mid=[self.b(A, B) if self.old_ok(A, B) else 0, self.b(A - 1, B - 1) if self.old_ok(A - 1, B - 1) else 0])
        if r["L"].size:
            out["L_norm"] = float(np.linalg.norm(r["L"], 2))
            out["adjoint_resid"] = float(np.abs(r["Lp"] + r["L"].conj().T - r["C"]).max())
        return out

    def cf_dressing_levels(self, A, B):
        """Bare-fraction spectra at three levels of the closed form, on the exact admissible sources (D8 sec. 9.6):
          first : sources at first order (Z instead of Zc, no self-energy terms), M -> 1, middle dressing from h only;
          self  : exact sources, M -> 1 + A (+) A'  (self-energies kept, transfer L, L' and leak dropped);
          exact : M u = sigma(h).
        Also returns the relative mismatch between the exact quadratic form and |Psi|^2 (a consistency check)."""
        Hadm, Psi, _ = self.cf_bps(A, B)
        k = Hadm.shape[1]
        if k == 0:
            return None
        secs, dims, off = self.cf_layout(A, B)
        h01, h00, h11, h10 = (Hadm[off[i]:off[i + 1]] for i in range(4))
        r = self.reduced_system(A, B)
        Uz, Ux, M, nz, nx = r["Uz"], r["Ux"], r["M"], r["nz"], r["nx"]
        SG, SpG = self.self_energies(A, B)
        t10 = np.zeros((self._d(A + 1, B - 1), k), complex)
        if t10.shape[0]:
            if h00.size:
                t10 += self.op("E", A, B) @ h00
            if h11.size:
                t10 -= self.op("K", A - 1, B - 1) @ h11
        t10_1 = t10.copy()
        if t10.shape[0] and h01.size:
            t10 += self.Zc(A, B) @ h01
            t10_1 += self.op("Z", A, B - 1) @ h01
        t01 = np.zeros((self._d(A - 2, B), k), complex)
        if t01.shape[0]:
            if h00.size:
                t01 += self.opd("K", A, B) @ h00
            if h11.size:
                t01 += self.opd("E", A - 1, B - 1) @ h11
        t01_1 = t01.copy()
        if t01.shape[0] and h10.size:
            t01 += self.Zc(A - 2, B + 1).conj().T @ h10
            t01_1 += self.opd("Z", A - 1, B) @ h10
        z = lambda d: np.zeros((d, k), complex)
        sz = Uz.conj().T @ self.Qp(A - 1, B) @ t10 - (Uz.conj().T @ SG @ h10 if h10.size else 0) if nz else z(0)
        sx = Ux.conj().T @ self.Qdp(A, B - 1) @ t01 - (Ux.conj().T @ SpG @ h01 if h01.size else 0) if nx else z(0)
        sz1 = Uz.conj().T @ self.Qp(A - 1, B) @ t10_1 if nz else z(0)
        sx1 = Ux.conj().T @ self.Qdp(A, B - 1) @ t01_1 if nx else z(0)
        sig, sig1 = np.vstack([sz, sx]), np.vstack([sz1, sx1])
        u_ex = np.linalg.solve(M, sig) if M.size else sig
        Dg = np.eye(nz + nx, dtype=complex)
        if nz:
            Dg[:nz, :nz] += r["A"]
        if nx:
            Dg[nz:, nz:] += r["Ap"]
        u_self = np.linalg.solve(Dg, sig) if Dg.size else sig

        def dressing(u, outer):
            P10 = (h10 + (Uz @ u[:nz] if (outer and nz) else 0)) if h10.size else None
            P01 = (h01 + (Ux @ u[nz:] if (outer and nx) else 0)) if h01.size else None
            Dm = u.conj().T @ u
            if P10 is not None:
                Dm = Dm + P10.conj().T @ SG @ P10
            if P01 is not None:
                Dm = Dm + P01.conj().T @ SpG @ P01
            return Dm
        N0 = Hadm.conj().T @ Hadm
        out = {}
        for name, (u, outer) in (("first", (sig1, False)), ("self", (u_self, True)), ("exact", (u_ex, True))):
            Gm = N0 + dressing(u, outer)
            Li = np.linalg.inv(np.linalg.cholesky((Gm + Gm.conj().T) / 2))
            out[name] = np.sort(np.linalg.eigvalsh(Li @ N0 @ Li.conj().T))
        Gt = Psi.conj().T @ Psi
        out["consistency"] = float(np.abs(N0 + dressing(u_ex, True) - Gt).max() / np.abs(Gt).max())
        return out

    # ------------------------------------------------------------------ (IV) closed-form decoders
    def cf_decoders(self, A, B):
        """Bare-fraction spectrum and CV slot-decoder spectra from the closed-form pencil on the admissible sources:
            mu: <h,h> = mu <Psi(h),Psi(h)>;   lambda_T: <Psi_T(h),Psi_T(h)> = lambda <Psi(h),Psi(h)>."""
        Hadm, Psi, _ = self.cf_bps(A, B)
        if Hadm.shape[1] == 0:
            return None
        Gm = Psi.conj().T @ Psi
        Li = np.linalg.inv(np.linalg.cholesky((Gm + Gm.conj().T) / 2))
        mu = np.sort(np.linalg.eigvalsh(Li @ (Hadm.conj().T @ Hadm) @ Li.conj().T))
        secs, dims, off = self.cf_layout(A, B)
        lam = {}
        for i, T in enumerate(self.CF_SLOTS):
            if dims[i]:
                PT = Psi[off[i]:off[i + 1]]
                lam["".join(map(str, T))] = np.sort(np.linalg.eigvalsh(Li @ (PT.conj().T @ PT) @ Li.conj().T))
        return dict(mu=mu, lam=lam)


# ---------------------------------------------------------------------------------------------------- utilities
def memory_estimate_gb(N):
    """Peak memory estimate for the full analysis at N -> N+1 (dense new sectors): ~12 dense copies of the largest new
    sector (Laplacian, eigenvectors, cached Q blocks and BPS bases of all sectors) plus the cached old blocks.
    Measured: 3.7 GB at N = 7 (estimate 4.9 GB)."""
    Np = N + 1
    dmax = max(comb(Np, A) * comb(Np, B) for A in range(Np + 1) for B in range(Np + 1))
    return dict(dmax=dmax, total=12 * dmax ** 2 * 16 / 1e9 + 0.3)
