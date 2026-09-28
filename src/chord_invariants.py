r"""
chord_invariants.py -- "chord data" of BPS-projected occupation operators in one-flavor N=2 SYK.

Purpose (research action item 1, see docs/research_questions.md Q1, Q4)
----------------------------------------------------------------------
A gravitational / chord-diagram calculation does not output a histogram. It outputs moments organized by
(i) two-point weights, (ii) crossing weights between different operators, (iii) higher connected pieces.
This module measures exactly those invariants for the operators

    Pi_i  = 1 - n_i                  (single-mode slot projector; the uplift decoder for i = N'-1)
    Pi_T  = prod_{i in T} (1 - n_i)  (k-mode slot projector, |T| = k; CV's multi-mode decoder)

compressed to the BPS space B = ker Q cap ker Q^dag of the (enlarged) theory with N' modes at charge P.
Note: with i.i.d. couplings all N' modes are statistically equivalent, so the decoder of the "new" mode N'-1 is
statistically identical to the compression of any other single-mode projector. We average over several modes.

Normalized trace tau(.) = Tr_B(.)/d, d = dim B, D = C(N',P), a = d/D.

Invariants (per realization)
----------------------------
single operator X = P Pi P (law on [0,1], d eigenvalues):
    b_UV   = Tr_D(Pi)/D            (UV mean; b = 1 - P/N' for one mode)
    free cumulants kappa_n, n = 2..6 of the law of X  (NC moment-cumulant relation)
    ratio rho_n = kappa_n / kappa_n^{Bern(b_UV)}.  Wachter/Haar (free compression) predicts rho_n = a^{n-1};
    a "Wachter with effective parameter" law predicts rho_n = a_eff^{n-1} for all n.
    r := rho_2 = Var/(b(1-b))       (fraction of UV variance surviving BPS projection)
pairs of *different* operators (disjoint mode sets), centered Ahat = P (Pi - mu) P, mu = tau(P Pi P):
    c   = tau(Ai Aj)/sqrt(tau(Ai^2) tau(Aj^2))       (covariance; nonzero since sum_i n_i = P)
    F   = tau(Ai Ai Aj Aj)/(tau(Ai^2) tau(Aj^2))      (non-crossing factorization; chord/free/q-Gaussian => 1)
    X   = tau(Ai Aj Ai Aj)/(tau(Ai^2) tau(Aj^2))      (crossing weight; free => 0, commuting UV => 1)
triples, all 6-letter "pair words" w (each of i,j,k twice), up to rotation:
    tau(w)/prod tau(A^2)  vs the chord (q-Gaussian) rule X^{cr(w)}, cr = number of chord crossings.

Two representations of P_B (identical results, chosen for cost):
    'basis'      : orthonormal BPS basis B (D x d); ops become d x d matrices.   Good when d is small (q=3).
    'complement' : orthonormal basis R (D x r) of (ker Q cap ker Q^dag)^perp, P = 1 - R R^dag;
                   traces expanded in R.  Good when r << d (q=5, where a ~ 0.9).
'haar=True' replaces P_B by a Haar-random projector of the same rank (null model; free compression).

Everything is exact linear algebra (no stochastic trace estimation). Couplings: i.i.d. unit complex Gaussian,
generated exactly as in src/q_scan.rand_form / src/decoder_fast.couplings (same seeds => same realizations).

Memory (added 2026-09-28 after the q=3 N'=16 / q=5 N'=17 runs were OOM-killed)
-------------------------------------------------------------------------------
* The BPS basis for 'basis' is obtained either from a full SVD of M = [Q; Q^dag] (null_method='svd'; used for all
  runs up to D <= 7000) or from the lowest eigenvectors of the dense Hermitian H = M^dag M = {Q, Q^dag}
  (null_method='eigh'; default for D > 7000).  B = ker H; the rank cut is certified by a spectral gap
  (all returned eigenvalues <= tol while the window extends to 1e-6 ||H||).  Peak ~ D^2 + D d complex numbers
  instead of ~ 2 D^2 + rows^2 + LAPACK workspace.
* Compressed letters (d x d or r x r matrices) are kept in a byte-bounded LRU cache (cache_gb).  Keys are the
  *multiset* of letters in a segment (diagonal operators commute), so a 6-letter triple needs <= 26 entries.
  The old unbounded cache keyed on raw bytes grew to tens of GB for r ~ 3000.
* estimate_memory_gb() gives the expected peak before running; use it and scripts/memwatch.py (CLAUDE.md).
"""
from __future__ import annotations

import hashlib
import itertools
from collections import OrderedDict
from functools import lru_cache
from math import comb

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

from q_scan import form_wedge, rand_form, subset_index


# ----------------------------------------------------------------------------- free cumulants
@lru_cache(maxsize=None)
def nc_partitions(n: int):
    """All non-crossing partitions of {0..n-1} as tuples of blocks (tuples)."""
    out = []

    def gen(rest, blocks):
        if not rest:
            out.append(tuple(tuple(b) for b in blocks))
            return
        first, others = rest[0], rest[1:]
        for r in range(len(others) + 1):
            for comb_ in itertools.combinations(others, r):
                blk = (first,) + comb_
                gen(tuple(x for x in others if x not in comb_), blocks + [blk])

    gen(tuple(range(n)), [])

    def crossing(p):
        for A, B in itertools.combinations(p, 2):
            for a1, a2 in itertools.combinations(sorted(A), 2):
                for b1, b2 in itertools.combinations(sorted(B), 2):
                    if a1 < b1 < a2 < b2 or b1 < a1 < b2 < a2:
                        return True
        return False

    return tuple(p for p in out if not crossing(p))


def free_cumulants(moments):
    """moments = [m_1, ..., m_K] (raw).  Returns [kappa_1, ..., kappa_K] via the NC moment-cumulant formula."""
    K = len(moments)
    kap = [0.0] * (K + 1)
    for n in range(1, K + 1):
        s = 0.0
        for p in nc_partitions(n):
            if len(p) == 1:
                continue
            prod = 1.0
            for blk in p:
                prod *= kap[len(blk)]
            s += prod
        kap[n] = moments[n - 1] - s
    return kap[1:]


def bernoulli_free_cumulants(b, K=6):
    return free_cumulants([b] * K)


def free_compression_moments(b, t, K=6):
    """Moments of the free compression (in the compressed algebra) of Bern(b) by a free projection of trace t:
    kappa_n = t^{n-1} kappa_n(Bern b)  (Nica-Speicher).  t = a gives the Wachter(a,b) law."""
    kB = bernoulli_free_cumulants(b, K)
    kap = [t ** (n - 1) * kB[n - 1] for n in range(1, K + 1)]
    m = []
    for n in range(1, K + 1):
        s = 0.0
        for p in nc_partitions(n):
            prod = 1.0
            for blk in p:
                prod *= kap[len(blk) - 1]
            s += prod
        m.append(s)
    return m


# ----------------------------------------------------------------------------- pair words
def crossing_number(word):
    """word: sequence of labels, each appearing exactly twice. Number of crossing pairs of chords on a circle."""
    pos = {}
    for i, c in enumerate(word):
        pos.setdefault(c, []).append(i)
    chords = [tuple(v) for v in pos.values()]
    cr = 0
    for (a1, a2), (b1, b2) in itertools.combinations(chords, 2):
        if a1 < b1 < a2 < b2 or b1 < a1 < b2 < a2:
            cr += 1
    return cr


def six_letter_words():
    """Representatives (up to rotation) of words using letters 0,1,2 each twice."""
    seen, reps = set(), []
    for w in set(itertools.permutations((0, 0, 1, 1, 2, 2))):
        rots = [w[i:] + w[:i] for i in range(6)]
        key = min(rots)
        if key in seen:
            continue
        seen.add(key)
        reps.append(key)
    return sorted(reps)


# ----------------------------------------------------------------------------- memory
_C16 = 16  # bytes per complex128


def constraint_rows(Np, P, q):
    return (comb(Np, P + q) if P + q <= Np else 0) + (comb(Np, P - q) if P - q >= 0 else 0)


def estimate_memory_gb(Np, P, q, method="basis", null_method="auto", cache_gb=4.0, d=None):
    """Rough peak memory (GB) of one BPSProjector + word evaluations.  d defaults to the generic-rank value
    D - rows (a lower bound on d; the exact d is slightly larger because Q^2 = 0 makes M rank deficient)."""
    D, rows = comb(Np, P), constraint_rows(Np, P, q)
    d = max(D - rows, 1) if d is None else d
    r = D - d
    if null_method == "auto" or method == "complement":
        null_method = "svd" if (D <= 7000 or method == "complement") else "eigh"
    if method == "complement":
        build = rows * D + rows * rows + r * D + 5 * rows * rows / 2    # M, U, Vh (reduced), rwork (in real)
        k = r
        keep = D * r
    elif null_method == "svd":
        m = min(rows, D)
        build = rows * D + rows * rows + D * D + D * d + 5 * m * m / 2 + 2 * m * max(rows, D) / 2
        k = d
        keep = D * d
    else:
        build = 1.5 * D * D + D * d                                     # sparse H (~30% dense) freed, dense H, Z
        k = d
        keep = D * d
    words = 4 * k * k                                                    # product temporaries + spectrum
    peak = max(build, keep + words) * _C16 / 1e9 + min(cache_gb, 26 * k * k * _C16 / 1e9)
    return dict(D=D, rows=rows, d_est=d, r_est=r, null_method=null_method, peak_gb=round(peak + 0.3, 2))


class _LRU:
    """Byte-bounded LRU cache of numpy arrays."""

    def __init__(self, max_bytes):
        self.max_bytes, self.nbytes, self.d = max_bytes, 0, OrderedDict()

    def get(self, key, make):
        if key in self.d:
            self.d.move_to_end(key)
            return self.d[key]
        val = make()
        self.d[key] = val
        self.nbytes += val.nbytes
        while self.nbytes > self.max_bytes and len(self.d) > 1:
            _, old = self.d.popitem(last=False)
            self.nbytes -= old.nbytes
        return val

    def clear(self):
        self.d.clear()
        self.nbytes = 0


def _digest(x):
    return hashlib.blake2b(np.ascontiguousarray(x).view(np.uint8), digest_size=16).digest()


# ----------------------------------------------------------------------------- BPS projector
class BPSProjector:
    """P_B for the q-form wedge model on Lambda^P(C^{N'}) (or a Haar projector of the same rank)."""

    def __init__(self, Np, P, q=3, seed=0, method="auto", haar=False, haar_seed=None, null_method="auto",
                 cache_gb=4.0):
        self.Np, self.P, self.q, self.seed, self.haar = Np, P, q, seed, haar
        self.masks, _ = subset_index(Np, P)
        self.D = len(self.masks)
        C = rand_form(seed, Np, q)
        blocks = [form_wedge(C, Np, P, q)]
        if P - q >= 0:
            blocks.append(form_wedge(C, Np, P - q, q).conj().transpose())
        Ms = sp.vstack(blocks).tocsr()
        # B = null space of M; its complement = row space of M.
        want_full = (method == "basis") or (method == "auto" and self.D <= 7000)
        if null_method == "auto":
            null_method = "svd" if self.D <= 7000 else "eigh"
        self.null_method = null_method if want_full else "svd"
        Vh = Bnull = None
        if Ms.shape[0] == 0:
            rank = 0
            Bnull = np.eye(self.D, dtype=complex)
        elif want_full and null_method == "eigh":
            Bnull = self._null_eigh(Ms)
            rank = self.D - Bnull.shape[1]
        else:
            M = Ms.toarray()
            _, s, Vh = np.linalg.svd(M, full_matrices=want_full)
            del M
            tol = max(Ms.shape) * s[0] * np.finfo(float).eps * 10
            rank = int((s > tol).sum())
        self.r = rank
        self.d = self.D - rank
        self.a = self.d / self.D
        if method == "auto":
            method = "basis" if (want_full and (self.d <= self.r or self.d <= 1500)) else "complement"
        self.method = method
        rng = np.random.default_rng(10_000_019 + 7919 * seed if haar_seed is None else haar_seed)
        self.B = self.R = None
        if method == "complement":
            if haar:
                G = rng.standard_normal((self.D, self.r)) + 1j * rng.standard_normal((self.D, self.r))
                self.R, _ = np.linalg.qr(G)
            else:
                if Vh is None:
                    raise RuntimeError("complement method needs the SVD row space; use null_method='svd'")
                self.R = np.ascontiguousarray(Vh[:rank].conj().T)       # D x r, orthonormal columns
        else:
            if haar:
                G = rng.standard_normal((self.D, self.d)) + 1j * rng.standard_normal((self.D, self.d))
                self.B, _ = np.linalg.qr(G)
            elif Bnull is not None:
                self.B = np.ascontiguousarray(Bnull)                    # D x d, orthonormal columns
            else:
                if Vh.shape[0] != self.D:
                    raise RuntimeError("basis method needs the full SVD; call with method='basis'")
                self.B = np.ascontiguousarray(Vh[rank:].conj().T)       # D x d, orthonormal columns
        del Vh, Bnull
        self._cache = _LRU(int(cache_gb * 1e9))

    def _null_eigh(self, Ms):
        """Orthonormal basis of ker M = ker H, H = M^dag M = {Q, Q^dag} on the charge-P sector.
        Eigenvectors with eigenvalue < thr = 1e-6 ||H|| are returned; they must all be <= tol (numerical zero),
        which certifies a gap of at least thr/tol between BPS and non-BPS states."""
        Hs = (Ms.conj().T @ Ms).tocsr()
        norm_bound = float(np.abs(Hs).sum(axis=1).max())                # Gershgorin bound on ||H||
        H = Hs.toarray()
        del Hs
        w, Z = sla.eigh(H, subset_by_value=(-np.inf, 1e-6 * norm_bound), driver="evr", overwrite_a=True,
                        check_finite=False)
        del H
        tol = self.D * norm_bound * np.finfo(float).eps * 100
        if len(w) and w.max() > tol:
            raise RuntimeError(f"no clean BPS gap: eigenvalue {w.max():.3e} in ({tol:.1e}, {1e-6 * norm_bound:.1e})")
        self.gap_window = (float(w.max()) if len(w) else 0.0, 1e-6 * norm_bound)
        return Z

    def clear_cache(self):
        self._cache.clear()

    # ---- diagonal operators
    def occ(self, i):
        return np.array([float((m >> i) & 1) for m in self.masks])

    def slot(self, T):
        """Pi_T = prod_{i in T} (1 - n_i) as a diagonal vector."""
        x = np.ones(self.D)
        for i in T:
            x *= 1.0 - self.occ(i)
        return x

    # ---- spectrum of P x P on B for a 0/1 diagonal x (a projector)
    def projector_spectrum(self, x):
        """Eigenvalues (d of them) of P_B X P_B restricted to B, X = diag(x) a coordinate projector."""
        if self.method == "basis":
            A = self.B.conj().T @ (x[:, None] * self.B)
            return np.clip(np.linalg.eigvalsh((A + A.conj().T) / 2), 0.0, 1.0)
        # complement: nonzero spectrum of X P X on ran X equals 1 - eig(Rs^dag Rs) (+ ones),
        # Rs = rows of R inside ran X.  Pad with zeros to d eigenvalues.
        idx = np.nonzero(x > 0.5)[0]
        m = len(idx)
        Rs = self.R[idx]
        s2 = np.linalg.eigvalsh(Rs.conj().T @ Rs).real          # r values in [0,1]
        mu = 1.0 - np.clip(s2, 0.0, 1.0)                         # eigenvalues of X P X on ran X (r of them)
        mu = np.concatenate([mu, np.ones(max(m - self.r, 0))])   # remaining ran X directions are in B
        nz = mu[mu > 1e-12]
        lam = np.concatenate([nz, np.zeros(self.d - len(nz))])
        return np.sort(np.clip(lam, 0.0, 1.0))

    # ---- traces of words  Tr_B( P x1 P x2 ... P xk )  for diagonal x's
    def _RYR(self, letters, keys):
        """R^dag diag(prod letters) R, cached by the multiset of letter digests (diagonals commute).
        The product is formed in sorted-key order, so the value does not depend on the word it came from."""
        order = sorted(range(len(keys)), key=lambda t: keys[t])
        ckey = ("R",) + tuple(keys[t] for t in order)

        def make():
            y = letters[order[0]].copy()
            for t in order[1:]:
                y *= letters[t]
            return (self.R.conj().T * y) @ self.R
        return self._cache.get(ckey, make)

    def _A(self, y, key=None):
        key = ("A", _digest(y) if key is None else key)
        return self._cache.get(key, lambda: self.B.conj().T @ (y[:, None] * self.B))

    def trace_word(self, xs):
        """Tr over B of the product (P x1 P)(P x2 P)...(P xk P) = Tr_D(P x1 P x2 ... P xk)."""
        k = len(xs)
        if self.method == "basis":
            M = self._A(xs[0])
            for x in xs[1:]:
                M = M @ self._A(x)
            return float(np.trace(M).real)
        keys = [_digest(x) for x in xs]
        # complement: expand each P = 1 - E, E = R R^dag, and sum over subsets S of positions carrying E.
        # For S = {s_0 < ... < s_{m-1}}:  Tr(E Y_0 E Y_1 ... E Y_{m-1}) = tr(prod_t R^dag Y_t R),
        # Y_t = x_{s_t} x_{s_t+1} ... x_{next-1}  (cyclically; next = s_{t+1}, or s_0 + k for the last segment).
        total = 0.0
        for S_mask in range(1 << k):
            S = [s_ for s_ in range(k) if (S_mask >> s_) & 1]
            sign = -1.0 if len(S) % 2 else 1.0
            if not S:
                y = np.ones(self.D)
                for x in xs:
                    y = y * x
                total += sign * float(y.sum())
                continue
            M = None
            m = len(S)
            for t in range(m):
                start_, nxt = S[t], (S[t + 1] if t + 1 < m else S[0] + k)
                seg = [u % k for u in range(start_, nxt)]
                Y = self._RYR([xs[u] for u in seg], [keys[u] for u in seg])
                M = Y if M is None else M @ Y
            total += sign * float(np.trace(M).real)
        return total

    def tau(self, xs):
        return self.trace_word(xs) / self.d


# ----------------------------------------------------------------------------- measurement
def single_invariants(proj, x, K=6):
    lam = proj.projector_spectrum(x)
    b_uv = float(x.sum() / proj.D)
    mom = [float(np.mean(lam ** k)) for k in range(1, K + 1)]
    kap = free_cumulants(mom)
    kB = bernoulli_free_cumulants(b_uv, K)
    rho = [kap[n] / kB[n] if abs(kB[n]) > 1e-14 else float("nan") for n in range(K)]
    return dict(b_uv=b_uv, mean=mom[0], moments=mom, kappa=kap, kappaB=kB, rho=rho,
                r=rho[1], n0=int(np.sum(lam < 1e-9)), n1=int(np.sum(lam > 1 - 1e-9))), lam


def centered_letters(proj, xs):
    """Return the centered diagonal vectors x - mu (mu = tau(P x P)) for the given diagonals."""
    out = []
    for x in xs:
        mu = proj.tau([x])
        out.append(x - mu)
    return out


def pair_invariants(proj, x1, x2):
    y1, y2 = centered_letters(proj, [x1, x2])
    v1, v2 = proj.tau([y1, y1]), proj.tau([y2, y2])
    c = proj.tau([y1, y2]) / np.sqrt(v1 * v2)
    F = proj.tau([y1, y1, y2, y2]) / (v1 * v2)
    X = proj.tau([y1, y2, y1, y2]) / (v1 * v2)
    return dict(c=float(c), F=float(F), X=float(X), X_over_F=float(X / F))


def triple_invariants(proj, x1, x2, x3):
    ys = centered_letters(proj, [x1, x2, x3])
    v = [proj.tau([y, y]) for y in ys]
    norm = v[0] * v[1] * v[2]
    out = []
    for w in six_letter_words():
        val = proj.tau([ys[c] for c in w]) / norm
        out.append(dict(word="".join("ijk"[c] for c in w), cr=crossing_number(w), ratio=float(val)))
    return out
