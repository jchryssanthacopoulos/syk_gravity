#!/usr/bin/env python3
r"""
Chord (gravitational) picture of the two-flavor pair uplift (docs/derivations.md D9; src/two_flavor_chords.py).

Tasks (one JSONL record per realization and sector / rung):
  window     BPS dimension, chord charge s, double-scaled BPS fraction D(s), window membership, and the lowest nonzero
             eigenvalues of the exact (B) and coexact (A) families vs the chord continuum edges, for every sector of the
             N'-pair model (--Np).  With --dmax only sectors up to that dimension (used for N' = 9, 10).   [per sector]
  blocks     BPS-averaged two-point functions of the uplift matter operators E, K, Z and daggers (unit-variance
             independent couplings, exact average over couplings estimated with random phases) in the old N-pair model:
             BPS -> BPS weight and time-integrated BPS -> continuum function, exact vs chord.            [per sector, op]
  resonance  exact bare-fraction spectrum mu (dressing kappa = 1/mu - 1) of every new BPS sector of the uplift N -> N+1,
             the chord charge s' of the sector, the old-model target sector of each first-order dressing term with its
             chord charge, resonance flag and exact lowest gap; the tower's bare fraction (symmetrized).  [per new sector]
  pairing    slot weights of the fortuitous states (BPS minus tower) of every new sector and of the tower:
             <n_psi_e>, <n_chi_f>, connected correlator c = w11 - <n_psi_e><n_chi_f>.                    [per new sector]
  edgegap    odd old N: the low spectrum of the family that is gapless (in the chord model) in each window-edge
             sector s = +-1/2 (exact family at +1/2, coexact at -1/2), relative to the sector's mean energy.  [per edge sector]
  twisted    coherence (G metric) of the E-emission of B(n,n) with the K-emission of B(n-1,n-1), tower projected out:
             the twisted chord's visibility on black-hole states, symmetrized vs generic; the tower's own cosine. [per n]
  tower      symmetrized only: the tower's individual emissions E Psi00, K Psi11 (and K^dag Psi00, E^dag Psi11) into
             the old target sectors, G-weighted, vs the black-hole (old BPS) averages of the same operators; the
             cancellation residuals; target charge and lowest gap.                                         [per rung n]

Couplings: rng = default_rng(seed); C1 = random_symmetric_couplings(N+1, rng) (or random_generic_couplings with
--model generic), exactly as scripts/run_two_flavor_uplift.py, so "old model N, seed s" is the same realization as in
results/data/two_flavor_2026-10-23.  For --tasks window the N'-pair model is TwoFlavor(C1) with C1 drawn at size N'.

Usage:  .venv/bin/python scripts/run_two_flavor_chords.py --tasks blocks tower --N 3 4 5 6 7 --seeds 1 2 3
        .venv/bin/python scripts/run_two_flavor_chords.py --tasks window --Np 4 5 6 7 8 --seeds 1 2 3
        .venv/bin/python scripts/run_two_flavor_chords.py --tasks resonance pairing --N 7 --seeds 1 --dry-run
Run N = 7 (resonance/pairing) and N' >= 8 (window) under scripts/memwatch.py.
"""
import argparse
import itertools
import json
import os
import sys
import time
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from run_two_flavor_uplift import git_state  # noqa: E402
from two_flavor import TwoFlavor, random_generic_couplings, random_symmetric_couplings, vacuum_tower  # noqa: E402
from two_flavor_chords import (DRESSING_TERMS, UPLIFT_OPS, bps_fraction, dressing_targets, edge_energies,  # noqa: E402
                               kappa_in, lam, s_eff, two_point)
from two_flavor_uplift import SLOTS, PairUplift, memory_estimate_gb  # noqa: E402

OUT = os.path.join(ROOT, "results/data/two_flavor_chords_2026-10-25")


def couplings(n_pairs, seed, model):
    rng = np.random.default_rng(seed)
    return (random_symmetric_couplings if model == "sym" else random_generic_couplings)(n_pairs, rng)


def q(x, ps=(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)):
    return [float(v) for v in np.quantile(x, ps)] if len(x) else []


# ----------------------------------------------------------------------------------------------- exact spectra
class Spectra:
    """Dense Laplacian spectra of the sectors of a TwoFlavor model, with the exact (B) / coexact (A) family split."""

    def __init__(self, model: TwoFlavor):
        self.M = model
        self.c = {}

    def Qd(self, a, b):
        return self.M.Q(a, b).toarray() if self.M.dim(a, b) and self.M.dim(a + 2, b - 1) else None

    @staticmethod
    def _gram_eigs(Q):
        """Nonzero spectrum of Q^dag Q (= that of Q Q^dag): eigenvalues of the smaller Gram matrix."""
        G = Q @ Q.conj().T if Q.shape[0] <= Q.shape[1] else Q.conj().T @ Q
        return np.linalg.eigvalsh((G + G.conj().T) / 2)

    def get(self, a, b, vectors=False):
        key = (a, b, vectors)
        if key in self.c:
            return self.c[key]
        d = self.M.dim(a, b)
        Qo, Qi = self.Qd(a, b), self.Qd(a - 2, b + 1)
        # The Laplacian L = Qo^dag Qo + Qi Qi^dag is the direct sum of the coexact (A) and exact (B) families, whose
        # nonzero spectra are those of the Gram matrices of Qo and Qi; the rest is the BPS kernel.
        wA = self._gram_eigs(Qo) if Qo is not None else np.zeros(0)
        wB = self._gram_eigs(Qi) if Qi is not None else np.zeros(0)
        scale = max([1.0] + [float(x.max()) for x in (wA, wB) if len(x)])
        tol = 1e-9 * scale
        nzA, nzB = wA[wA > tol], wB[wB > tol]
        U = None
        if vectors:
            L = np.zeros((d, d), complex)
            if Qo is not None:
                L += Qo.conj().T @ Qo
            if Qi is not None:
                L += Qi @ Qi.conj().T
            w, U = np.linalg.eigh((L + L.conj().T) / 2)
        else:
            w = np.sort(np.concatenate([np.zeros(d - len(nzA) - len(nzB)), nzA, nzB]))
        rec = dict(w=w, U=U, tol=tol, nbps=int(d - len(nzA) - len(nzB)), m1=float((nzA.sum() + nzB.sum()) / d) if d else 0.0,
                   gapA=float(nzA.min()) if len(nzA) else None, gapB=float(nzB.min()) if len(nzB) else None,
                   lowB=[float(x) for x in np.sort(nzB)[:6]], lowA=[float(x) for x in np.sort(nzA)[:6]])
        self.c[key] = rec
        return rec


# ------------------------------------------------------------------------------------------------------- tasks
def task_window(Np, seed, model, dmax):
    M = TwoFlavor(couplings(Np, seed, model))
    S = Spectra(M)
    out = []
    for a in range(Np + 1):
        for b in range(Np + 1):
            d = M.dim(a, b)
            if d == 0 or (dmax and d > dmax):
                continue
            r = S.get(a, b)
            s = s_eff(a, b, Np)
            eB, eA = edge_energies(a, b, Np, r["m1"]) if r["m1"] > 0 else (None, None)
            out.append(dict(Np=Np, a=a, b=b, dim=d, s=s, lam=lam(Np), nbps=r["nbps"], frac=r["nbps"] / d,
                            D_chord=bps_fraction(a, b, Np), in_window=bool(abs(s) < 0.5 - 1e-9), edge=bool(abs(abs(s) - 0.5) < 1e-9),
                            m1=r["m1"], gapB=r["gapB"], gapA=r["gapA"], lowB=r["lowB"], lowA=r["lowA"],
                            chord_edgeB=eB, chord_edgeA=eA, tower_sector=bool(a == b)))
            S.c.clear()
    return out


def task_twisted(P, symmetrized):
    """Coherence of the E-emission of B(n,n) with the K-emission of B(n-1,n-1) into old (n+1, n-1), in the G metric:
    mean squared cosine between random elements, c2 = |A^dag B|_F^2 / (|A|_F^2 |B|_F^2), A = G^1/2 E W00, B = G^1/2 K W11,
    with the tower projected out of both BPS spaces (symmetrized model).  E and K are independent in the generic model;
    in the symmetrized one their correlation is the twisted chord (D9 sec. 1.2).  Also the random-subspace baseline
    1/rank of the emission space and, for the tower, its own cosine (= 1 by the Ward identity)."""
    So = Spectra(P.old)
    out = []
    for n in range(2, P.N):
        ta, tb = n + 1, n - 1
        if not (P.old_dim(ta, tb) and P.b(n, n) > 1 and P.b(n - 1, n - 1) > 1):
            continue
        r = So.get(ta, tb, vectors=True)
        w, U = r["w"], r["U"]
        nz = w > r["tol"]
        Gh = (U[:, nz] / np.sqrt(w[nz])) @ U[:, nz].conj().T
        W00, W11 = P.old_B(n, n), P.old_B(n - 1, n - 1)
        rec = dict(n=n, target=[ta, tb], s_target=s_eff(ta, tb, P.N), b00=int(W00.shape[1]), b11=int(W11.shape[1]))
        if symmetrized:
            t0, t1 = P.tower_state(n), P.tower_state(n - 1)
            a, b = Gh @ (P.E(n, n) @ t0), Gh @ (P.K(n - 1, n - 1) @ t1)
            rec["tower_cos2"] = float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))

            def perp(W, t):
                W = W - np.outer(t, t.conj() @ W)
                Uq, sv, _ = np.linalg.svd(W, full_matrices=False)
                return Uq[:, sv > 1e-8]
            W00, W11 = perp(W00, t0), perp(W11, t1)
        A, Bm = Gh @ (P.E(n, n) @ W00), Gh @ (P.K(n - 1, n - 1) @ W11)
        nA, nB = np.linalg.norm(A) ** 2, np.linalg.norm(Bm) ** 2
        rec.update(c2_mean=float(np.linalg.norm(A.conj().T @ Bm) ** 2 / (nA * nB)),
                   c2_max=float(np.linalg.svd(A.conj().T @ Bm, compute_uv=False).max() ** 2 /
                                (np.linalg.svd(A, compute_uv=False).max() * np.linalg.svd(Bm, compute_uv=False).max()) ** 2),
                   rank_target=int(np.sum(nz)))
        out.append(rec)
    return out


def task_edgegap(N, seed, model):
    """Old N-pair model (odd N): for every window-edge sector (s = +-1/2) the low spectrum of the family that the chord
    model makes gapless there (exact/B at s = +1/2: Q from (a-2, b+1); coexact/A at s = -1/2: Q to (a+2, b-1)).  Each Q map
    is diagonalized once (its nonzero spectrum is shared by the B family of its target and the A family of its source).
    m1 = (|Q_out|_F^2 + |Q_in|_F^2) / dim (no diagonalization needed)."""
    C1 = couplings(N + 1, seed, model)
    M = TwoFlavor(C1[:N, :N, :N].copy())
    gram = {}

    def spec(src):
        if src not in gram:
            Q = M.Q(*src).toarray()
            gram[src] = (np.sort(Spectra._gram_eigs(Q)), min(Q.shape))
        return gram[src]

    out = []
    for a in range(N + 1):
        for b in range(N + 1):
            d = M.dim(a, b)
            s = s_eff(a, b, N)
            if d == 0 or abs(abs(s) - 0.5) > 1e-9:
                continue
            src = (a - 2, b + 1) if s > 0 else (a, b)
            if not (M.dim(*src) and M.dim(src[0] + 2, src[1] - 1)):
                continue
            w, _ = spec(src)
            tol = 1e-9 * max(1.0, w.max())
            nz = w[w > tol]
            fro = 0.0
            for (x, y) in ((a, b), (a - 2, b + 1)):
                if M.dim(x, y) and M.dim(x + 2, y - 1):
                    fro += float(np.sum(np.abs(M.Q(x, y).data) ** 2))
            m1 = fro / d
            out.append(dict(N=N, a=a, b=b, dim=d, s=s, family="B" if s > 0 else "A", family_dim=int(len(nz)), m1=m1,
                            gap=float(nz.min()) if len(nz) else None, gap_rel=float(nz.min() / m1) if len(nz) else None,
                            low_rel=[float(x / m1) for x in nz[:8]], mean_family=float(nz.mean()) if len(nz) else None))
    return out


def _strings(model, op, sign):
    P, X, N = model.psi_pos, model.chi_pos, model.Ns
    if op == "E":
        s = [[(P[i], True), (X[k], False)] for i in range(N) for k in range(N)]
    elif op == "K":
        s = [[(P[i], True), (P[j], True)] for i, j in itertools.combinations(range(N), 2)]
    else:
        s = [[(P[i], True)] for i in range(N)]
    if sign < 0:                                  # dagger: reverse order, flip daggers
        s = [[(pos, not dag) for pos, dag in reversed(st)] for st in s]
    return s


def task_blocks(N, seed, model, ndraw):
    C1 = couplings(N + 1, seed, model)
    M = TwoFlavor(C1[:N, :N, :N].copy())
    S = Spectra(M)
    rng = np.random.default_rng(10_000 + seed)
    out = []
    for a in range(N + 1):
        for b in range(N + 1):
            if M.dim(a, b) == 0:
                continue
            src = S.get(a, b, vectors=True)
            W = src["U"][:, src["w"] <= src["tol"]]
            d = W.shape[1]
            if d == 0:
                continue
            for op in ("E", "K", "Z"):
                for sign in (+1, -1):
                    (da, db), _ = UPLIFT_OPS[op]
                    ta, tb = a + sign * da, b + sign * db
                    if M.dim(ta, tb) == 0:
                        continue
                    tg = S.get(ta, tb, vectors=True)
                    w, U = tg["w"], tg["U"]
                    nz = w > tg["tol"]
                    strs = _strings(M, op, sign)
                    vals = []
                    for _ in range(ndraw):
                        eta = np.exp(2j * np.pi * rng.random(len(strs)))
                        X = M.op_block([(e, st) for e, st in zip(eta, strs)], (a, b), (ta, tb))
                        Y = U.conj().T @ (X @ W)
                        p2 = np.sum(np.abs(Y) ** 2, axis=1) / d
                        vals.append((p2.sum(), np.sum(p2[nz] / w[nz]), p2[~nz].sum()))
                    vals = np.array(vals)
                    c0, T, pm = vals.mean(0)
                    ch = two_point(op, sign, a, b, N)
                    kap = kappa_in(ta, tb, N, tg["m1"], sign)
                    out.append(dict(N=N, a=a, b=b, d=d, op=op, sign=sign, target=[ta, tb], s_out=ch["s_out"], s_in=ch["s_in"],
                                    delta=ch["delta"], c0=float(c0), T_exact=float(T), dT=float(vals[:, 1].std(ddof=1) / np.sqrt(ndraw)),
                                    pm_exact=float(pm / c0) if c0 > 0 else None, pm_chord=ch["point_mass"],
                                    T_chord=float(c0 * ch["I"] / kap) if (np.isfinite(ch["I"]) and kap > 0) else None,
                                    m1_target=tg["m1"], edge_target=bool(abs(ch["s_in"] - 0.5) < 1e-9)))
    return out


def task_new_sectors(P, symmetrized, want_res=True, want_pair=True):
    """resonance and pairing records, one pass over the new sectors (each new BPS basis is computed once)."""
    So = Spectra(P.old)
    res, pair = [], []
    Np = P.N + 1
    for A in range(Np + 1):
        for B in range(Np + 1):
            if not P.new_ok(A, B) or P.d_new(A, B) == 0:
                continue
            if want_res:
                mu = np.clip(P.bare_fraction_spectrum(A, B), 1e-300, 1.0)
                kap = 1 / mu - 1
                terms = {}
                for name, t in dressing_targets(A, B, Np).items():
                    ta, tb = t["target"]
                    if not (P.old_ok(ta, tb) and P.old_dim(ta, tb) and P.old_ok(*t["source"]) and P.b(*t["source"])):
                        continue
                    r = So.get(ta, tb)
                    gap = r["gapB"] if DRESSING_TERMS[name][2] > 0 else r["gapA"]
                    terms[name] = dict(target=[ta, tb], s_target=t["s_target"], resonant=t["resonant"], gap=gap,
                                       gap_rel=gap / r["m1"] if gap else None)
                rec = dict(A=A, B=B, d_new=len(mu), s_new=s_eff(A, B, Np), kappa_mean=float(kap.mean()), kappa_q=q(kap),
                           kappa_max=float(kap.max()), mu_min=float(mu.min()), n_bare=int(np.sum(mu > 1 - 1e-9)),
                           resonant=[k for k, v in terms.items() if v["resonant"]], terms=terms, kappa=kap.tolist())
                if symmetrized and A == B:
                    rec["tower_bare_fraction"] = float(P.tower_slot_weights(A)[1])
                res.append(rec)
            if want_pair:
                Bn = P.new_B(A, B)
                rec = dict(A=A, B=B, s_new=s_eff(A, B, Np))
                if symmetrized and A == B:
                    t = P.tower_state(A, model="new")
                    comp_t = P.slot_components(A, B, t[:, None])
                    wt = {f"{T[0]}{T[1]}": float(np.sum(np.abs(comp_t[T]) ** 2)) for T in SLOTS}
                    rec["tower_w"] = wt
                    rec["tower_c"] = wt["11"] - (wt["10"] + wt["11"]) * (wt["01"] + wt["11"])
                    Bn = Bn - np.outer(t, t.conj() @ Bn)
                    U, s, _ = np.linalg.svd(Bn, full_matrices=False)
                    Bn = U[:, s > 1e-8]
                d = Bn.shape[1]
                rec["d_fort"] = d
                if d:
                    comp = P.slot_components(A, B, Bn)
                    ws = {f"{T[0]}{T[1]}": np.sum(np.abs(comp[T]) ** 2, axis=0) for T in SLOTS}
                    cs = ws["11"] - (ws["10"] + ws["11"]) * (ws["01"] + ws["11"])
                    w = {k: float(v.mean()) for k, v in ws.items()}
                    rec.update(w=w, n_psi=w["10"] + w["11"], n_chi=w["01"] + w["11"], n_psi_site=A / Np, n_chi_site=B / Np,
                               c_mean=w["11"] - (w["10"] + w["11"]) * (w["01"] + w["11"]),
                               c_state_q=q(cs, (0.05, 0.5, 0.95)), c_state_max=float(cs.max()))
                pair.append(rec)
            P.clear_cache("nQ", "nB", "hemb", "sidx")
    return res, pair


def task_tower(P):
    So = Spectra(P.old)
    out = []
    for n in range(1, P.N + 1):
        psi00 = vacuum_tower(P.old, n)
        psi11 = n * vacuum_tower(P.old, n - 1)
        n00, n11 = np.linalg.norm(psi00), np.linalg.norm(psi11)
        for cur, (ta, tb), ops in (("J10", (n + 1, n - 1), ("E", "K")), ("J01", (n - 2, n), ("K+", "E+"))):
            if not P.old_ok(ta, tb) or P.old_dim(ta, tb) == 0:
                continue
            r = So.get(ta, tb, vectors=True)
            w, U = r["w"], r["U"]
            nz = w > r["tol"]
            if cur == "J10":
                Ox, Oy = P.E(n, n), P.K(n - 1, n - 1)                        # (n,n), (n-1,n-1) -> (n+1,n-1)
                gap = r["gapB"]
            else:
                Ox, Oy = P.K(n - 2, n).conj().T, P.E(n - 2, n).conj().T      # (n,n), (n-1,n-1) -> (n-2,n)
                gap = r["gapA"]
            x, y = Ox @ psi00, Oy @ psi11
            Wx, Wy = P.old_B(n, n), P.old_B(n - 1, n - 1)
            # the currents are J10 = E Psi00 - K Psi11, J01 = K^dag Psi00 + E^dag Psi11 in the operator convention; the code
            # basis can differ by a sign per slot, so report the smaller of |x - y|, |x + y|
            resid = min(np.linalg.norm(x - y), np.linalg.norm(x + y)) / max(np.linalg.norm(x), 1e-300)

            def G2(v):
                yv = U.conj().T @ v
                return float(np.sum(np.abs(yv[nz]) ** 2 / w[nz]))

            def bh(Oop, W):
                if W.size == 0 or W.shape[1] == 0:
                    return None
                Y = U.conj().T @ (Oop @ W)
                return float(np.sum(np.abs(Y[nz]) ** 2 / w[nz, None]) / W.shape[1])

            out.append(dict(n=n, current=cur, target=[ta, tb], s_target=s_eff(ta, tb, P.N), gap=gap,
                            gap_rel=gap / r["m1"] if gap else None, resid=float(resid),
                            tower_G_first=G2(x) / n00 ** 2 if n00 else None, tower_G_second=G2(y) / n11 ** 2 if n11 else None,
                            bh_G_first=bh(Ox, Wx), bh_G_second=bh(Oy, Wy), ops=list(ops)))
    return out


# -------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tasks", nargs="+", required=True, choices=["window", "blocks", "resonance", "pairing", "tower", "edgegap", "twisted"])
    ap.add_argument("--N", type=int, nargs="*", default=[], help="old-model sizes (blocks, resonance, pairing, tower)")
    ap.add_argument("--Np", type=int, nargs="*", default=[], help="model sizes for --tasks window")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    ap.add_argument("--model", choices=["sym", "generic"], default="sym")
    ap.add_argument("--dmax", type=int, default=0, help="window: skip sectors larger than this (0 = all)")
    ap.add_argument("--ndraw", type=int, default=6, help="blocks: random-phase draws per operator")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tag", default="")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    est = {}
    for N in a.N:
        dedge = max([comb(N, x) * comb(N, y) for x in range(N + 1) for y in range(N + 1)
                     if abs(abs(s_eff(x, y, N)) - 0.5) < 1e-9] or [0])
        est[f"N={N}"] = dict(resonance_pairing_tower=memory_estimate_gb(N)["total"],
                             blocks=round((comb(2 * N, N)) ** 2 * 48 / 1e9 + 0.3, 2),
                             edgegap=round(3.5 * dedge ** 2 * 16 / 1e9 + 0.3, 2))
    for Np in a.Np:
        dmax = max(comb(Np, x) * comb(Np, y) for x in range(Np + 1) for y in range(Np + 1))
        if a.dmax:
            dmax = min(dmax, a.dmax)
        est[f"Np={Np}"] = dict(window=round(6 * dmax ** 2 * 16 / 1e9 + 0.3, 2))
    print("memory estimate (GB):", json.dumps(est))
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    meta = dict(args=vars(a), git=git_state(), started=stamp, estimate_gb=est)
    with open(os.path.join(a.out, f"meta_{a.model}_{'-'.join(a.tasks)}{a.tag}_{stamp}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    sym = a.model == "sym"

    def dump(task, size_key, size, recs):
        fn = os.path.join(a.out, f"{task}_{size_key}{size}_{a.model}{a.tag}.jsonl")
        with open(fn, "a") as fh:
            for r in recs:
                fh.write(json.dumps(r) + "\n")

    for Np in a.Np:
        for seed in a.seeds:
            if "window" in a.tasks:
                t0 = time.time()
                recs = [dict(seed=seed, model=a.model, **r) for r in task_window(Np, seed, a.model, a.dmax)]
                dump("window", "Np", Np, recs)
                print(f"window N'={Np} seed={seed}: {len(recs)} sectors, {time.time() - t0:.1f}s", flush=True)
    for N in a.N:
        for seed in a.seeds:
            if "edgegap" in a.tasks:
                t0 = time.time()
                recs = [dict(seed=seed, model=a.model, **r) for r in task_edgegap(N, seed, a.model)]
                dump("edgegap", "N", N, recs)
                print(f"edgegap N={N} seed={seed}: {len(recs)} edge sectors, {time.time() - t0:.1f}s", flush=True)
            if "blocks" in a.tasks:
                t0 = time.time()
                recs = [dict(seed=seed, model=a.model, **r) for r in task_blocks(N, seed, a.model, a.ndraw)]
                dump("blocks", "N", N, recs)
                print(f"blocks N={N} seed={seed}: {len(recs)} records, {time.time() - t0:.1f}s", flush=True)
            if any(t in a.tasks for t in ("resonance", "pairing", "tower", "twisted")):
                P = PairUplift(couplings(N + 1, seed, a.model))
                if "twisted" in a.tasks:
                    t0 = time.time()
                    recs = [dict(N=N, seed=seed, model=a.model, **r) for r in task_twisted(P, sym)]
                    dump("twisted", "N", N, recs)
                    print(f"twisted N={N} seed={seed}: {len(recs)} records, {time.time() - t0:.1f}s", flush=True)
                if "tower" in a.tasks and sym:
                    t0 = time.time()
                    recs = [dict(N=N, seed=seed, model=a.model, **r) for r in task_tower(P)]
                    dump("tower", "N", N, recs)
                    print(f"tower N={N} seed={seed}: {len(recs)} records, {time.time() - t0:.1f}s", flush=True)
                wr, wp = "resonance" in a.tasks, "pairing" in a.tasks
                if wr or wp:
                    t0 = time.time()
                    res, pair = task_new_sectors(P, sym, wr, wp)
                    if wr:
                        dump("resonance", "N", N, [dict(N=N, seed=seed, model=a.model, **r) for r in res])
                    if wp:
                        dump("pairing", "N", N, [dict(N=N, seed=seed, model=a.model, **r) for r in pair])
                    print(f"resonance/pairing N={N} seed={seed}: {len(res)}/{len(pair)} records, {time.time() - t0:.1f}s",
                          flush=True)
                del P


if __name__ == "__main__":
    main()
