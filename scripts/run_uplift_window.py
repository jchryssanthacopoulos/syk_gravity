#!/usr/bin/env python3
r"""
The uplift anywhere in the fortuity window (docs/derivations.md D10): per sector (N', F) of the N'-mode model, the exact
uplift from N = N' - 1 modes, with its obstruction, decoder, dressing and inner spectra.

Task 'sectors'.  For each N', seed and window sector |l'| <= lmax (l' = 2F - N'), using uplift.OldModel(Np, c, rng, F=F):
  - slot dimensions, the singular values and ranks of the connecting maps delta_h = P_H E, delta_g = P_H E^dag, and the
    admissible dimensions d_h, d_g (BPS dimension of the sector = d_h + d_g, Theorem A of D8 / D10 §2);
  - direct check for N' <= --direct-max: BPS basis of the N'-mode model, residual of the uplift, rank, and the
    largest deviation between the closed-form decoder and the directly computed decoder;
  - decoder lambda (weight in n_e = 0) exactly (closed-form pencil), first order (1/(1+s^2), s'^2/(1+s'^2)), self-dressed
    (blocks without h-g mixing), and the D6 per-mode mean field; T2 = <mu^2>, T4, H = <y(1-y)>, walls, atoms;
  - dressing costs: first order (s^2, s'^2) and exact minimum-norm uplift (eigenvalues of G_hh - 1, G_gg - 1);
  - inner spectra (D10 §4): energies s_j^2 of the old supercharge through which each family dresses, the coupling
    weights |e_j|^2 = sum over sources |<u_j|E h>|^2, and r_j = W_jj/s_j^2 = <v_j|R|v_j>.
Couplings: rng = default_rng(seed), c = random_couplings(rng, C(N', p)): for p = 3 the same realization as
results/data/uplift_2026-10-12 (seeds 9000, 9001).  All sectors of one (N', seed) share c.

Task 'forms'.  For each N: the bilinear connecting form B(beta, alpha) = int beta ^ eta ^ alpha on harmonic forms
(uplift.connecting_form).  Even N (odd p): the self-paired level 1 - p, with its symmetry (-1)^k and radical (the
exceptional states at N + 1).  Odd N: the pairing H(-p) x H(2-p) (rank = exceptional count when it is full).

Usage (run under scripts/memwatch.py):
  .venv/bin/python scripts/run_uplift_window.py --task sectors --Np 5 6 7 8 9 10 11 12 13 14 --seeds 9000 9001 9002
  .venv/bin/python scripts/run_uplift_window.py --task forms --N 4 5 6 7 8 9 10 11 12 13 --seeds 9000 9001
  --dry-run prints the memory estimate (dense old-model maps: ~14 Dmax^2 x 16 B; N' = 14: 1.4 GB, N' = 16: 9.6 GB).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from math import comb

import numpy as np
import scipy.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import Sector, random_couplings  # noqa: E402
from uplift import (OldModel, bps_dim_q3, connecting_form, decoder_from_basis, harmonic_basis, level,  # noqa: E402
                    split_couplings, strip_walk_count, uplift_bookkeeping)

OUT = os.path.join(ROOT, "results/data/uplift_window_2026-10-26")


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", "src", "scripts"],
                               capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def mem_estimate_gb(Np, p=3, lmax=None):
    N = Np - 1
    D = max(comb(N, k) for k in range(N + 1))
    return 14 * D * D * 16 / 1e9 + 0.3


def window_sectors(Np, p, lmax):
    return [F for F in range(1, Np) if abs(level(Np, F)) <= lmax]


def moments(mu):
    y = np.asarray(mu, float) ** 2
    if len(y) == 0:
        return {}
    lam = (1 + np.asarray(mu)) / 2
    w = np.minimum(lam, 1 - lam)
    return dict(mean_lam=float(lam.mean()), T2=float(y.mean()), T4=float(np.mean(y * y)), H=float(np.mean(y * (1 - y))),
                walls={str(e): float(np.mean(w < e)) for e in (0.1, 0.03, 0.01, 0.001)},
                atoms=int(np.sum(w < 1e-9)))


def tail(x):
    x = np.sort(np.asarray(x, float))
    if len(x) == 0:
        return {}
    q = lambda f: float(np.quantile(x, f))  # noqa: E731
    return dict(n=len(x), mean=float(x.mean()), q50=q(0.5), q90=q(0.9), q99=q(0.99), max=float(x[-1]),
                frac_gt10=float(np.mean(x > 10)), frac_gt100=float(np.mean(x > 100)))


def mean_field_mu(KH, R):
    """D6 §5 per-mode mean field: along each singular direction of K_H, R -> rho_j = <k_j, R k_j>/|k_j|^2,
    mu_j = [1 - s_j^2 (1 - rho_j)/(1 + rho_j)^2] / [1 + s_j^2/(1 + rho_j)].  (Cloud-dominated limit: lambda -> rho/(1+rho).)"""
    if KH.shape[1] == 0:
        return np.zeros(0)
    U, s, _ = np.linalg.svd(KH, full_matrices=False)
    rho = np.real(np.einsum("ij,ij->j", U.conj(), R @ U))
    s2 = np.concatenate([s ** 2, np.zeros(KH.shape[1] - len(s))])
    rho = np.concatenate([rho, np.zeros(KH.shape[1] - len(s))])
    return (1 - s2 * (1 - rho) / (1 + rho) ** 2) / (1 + s2 / (1 + rho))


def inner_record(om, fam, nlow=12):
    dy = om.dyson(fam)
    s2, e, W = dy["s2"], dy["e"], dy["W"]
    if len(s2) == 0 or e.shape[1] == 0:
        return dict(n=int(len(s2)))
    we = np.sum(np.abs(e) ** 2, axis=1)                     # coupling of inner mode j to all admissible sources
    contrib = we / s2                                        # its share of Tr(first-order dressing)
    o = np.argsort(s2)
    tot = contrib.sum()
    r = np.real(np.diag(W)) / s2
    return dict(n=int(len(s2)), s2_min=float(s2.min()), s2_mean=float(s2.mean()), s2_low=[float(x) for x in s2[o[:nlow]]],
                we_low=[float(x) for x in we[o[:nlow]]], we_mean=float(we.mean()), r_low=[float(x) for x in r[o[:nlow]]],
                r_mean=float(r.mean()),
                share_low1=float(contrib[o[:max(1, len(o) // 100)]].sum() / tot),
                share_low5=float(contrib[o[:max(1, len(o) // 20)]].sum() / tot))


def sector_record(Np, F, c, seed, p, direct_max, rng):
    t0 = time.time()
    N = Np - 1
    bk = uplift_bookkeeping(Np, F, p)
    om = OldModel(Np, c, rng, p=p, F=F)
    ob = om.obstruction
    d_h, d_g = om.VF.shape[1], om.VF1.shape[1]
    d = d_h + d_g
    rec = dict(task="sectors", Np=Np, N=N, F=F, lp=level(Np, F), p=p, seed=seed, bookkeeping=bk,
               slot_dims=dict(a=ob["h"]["dim"], b=ob["g"]["dim"]),
               obstruction=dict(h=dict(rank=ob["h"]["rank"], sv=ob["h"]["sv"][:24]),
                                g=dict(rank=ob["g"]["rank"], sv=ob["g"]["sv"][:24])),
               d_h=d_h, d_g=d_g, d=d)
    try:
        rec["expected_d"] = bps_dim_q3(Np, F) if p == 3 else (strip_walk_count(Np, p, F) if abs(level(Np, F)) < p else None)
    except ValueError:
        rec["expected_d"] = None
    if d:
        Nm, G = om.decoder_closed_form()
        mu = sla.eigvalsh(Nm, G)
        mu_self = np.concatenate([sla.eigvalsh(Nm[:d_h, :d_h], G[:d_h, :d_h]) if d_h else np.zeros(0),
                                  sla.eigvalsh(Nm[d_h:, d_h:], G[d_h:, d_h:]) if d_g else np.zeros(0)])
        sig, sigp = om.first_order_singular_values()
        mu1 = np.concatenate([(1 - sig ** 2) / (1 + sig ** 2), -(1 - sigp ** 2) / (1 + sigp ** 2)])
        KH, KpH, L = om.three_operators()
        R, Rp = L @ L.conj().T, L.conj().T @ L
        mu_mf = np.concatenate([mean_field_mu(KH, R), -mean_field_mu(KpH, Rp)])
        Mh, Mg, C2, overlap = om.decoder_blocks()
        kap_h = np.linalg.eigvalsh(G[:d_h, :d_h]) - 1 if d_h else np.zeros(0)
        kap_g = np.linalg.eigvalsh(G[d_h:, d_h:]) - 1 if d_g else np.zeros(0)
        rec.update(lam=sorted(float(x) for x in (1 + mu) / 2), lam_first=sorted(float(x) for x in (1 + mu1) / 2),
                   lam_self=sorted(float(x) for x in (1 + mu_self) / 2),
                   exact=moments(mu), first=moments(mu1), self=moments(mu_self), mean_field=moments(mu_mf),
                   split_overlap=float(overlap), mixing_T2=float(2 * np.linalg.norm(C2) ** 2 / d) if C2.size else 0.0,
                   dress_first=dict(h=tail(sig ** 2) if d_h else {}, g=tail(sigp ** 2) if d_g else {}),
                   dress_exact=dict(h=tail(kap_h), g=tail(kap_g)),
                   sigma2_h=[float(x) for x in np.sort(sig ** 2)], sigma2_g=[float(x) for x in np.sort(sigp ** 2)],
                   kappa_h=[float(x) for x in np.sort(kap_h)], kappa_g=[float(x) for x in np.sort(kap_g)],
                   inner=dict(h=inner_record(om, "h"), g=inner_record(om, "g")))
    if Np <= direct_max:
        SN = Sector(Np, p, F=F)
        V = harmonic_basis(SN, c, rng)
        rec["d_direct"] = int(V.shape[1])
        if d and V.shape[1]:
            Psi = om.exact_uplift(SN)
            sv = np.linalg.svd(Psi, compute_uv=False)
            rec["uplift_residual"] = float(np.linalg.norm(Psi - V @ (V.conj().T @ Psi)) / np.linalg.norm(Psi))
            rec["uplift_rank"] = int(np.sum(sv > 1e-8 * sv[0]))
            lam_dir = np.sort(decoder_from_basis(SN, V, Np))
            if len(lam_dir) == d:
                rec["decoder_dev"] = float(np.abs(lam_dir - np.array(rec["lam"])).max())
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def form_record(N, p, seed):
    Np = N + 1
    rng = np.random.default_rng(seed)
    c = random_couplings(rng, comb(Np, p))
    cA, eta = split_couplings(Np, p, c)
    rec = dict(task="forms", N=N, p=p, seed=seed)
    if (N - p + 1) % 2 == 0:
        k = (N - p + 1) // 2
        M, V, _ = connecting_form(N, p, k, cA, eta, rng)
        sgn = (-1) ** k
        s = np.linalg.svd(M, compute_uv=False) if M.size else np.zeros(0)
        rad = int(np.sum(s < 1e-9 * s[0])) if len(s) else 0
        rec.update(kind="self-paired", level=1 - p, k=k, dim=int(M.shape[0]), sign=sgn,
                   asym=float(np.abs(M.T - sgn * M).max() / max(np.abs(M).max(), 1e-300)) if M.size else 0.0,
                   radical=rad, sv_low=[float(x) for x in s[-8:]], sv_max=float(s[0]) if len(s) else 0.0)
    else:
        k = (N - p) // 2                                       # level -p
        M, Va, Vb = connecting_form(N, p, k, cA, eta, rng)
        s = np.linalg.svd(M, compute_uv=False) if M.size else np.zeros(0)
        rec.update(kind="pairing", levels=[-p, 2 - p], shape=list(M.shape),
                   rank=int(np.sum(s > 1e-9 * s[0])) if len(s) else 0, sv=[float(x) for x in s[:8]])
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", choices=["sectors", "forms"], default="sectors")
    ap.add_argument("--Np", type=int, nargs="*", default=[], help="N' values (task sectors)")
    ap.add_argument("--N", type=int, nargs="*", default=[], help="old-model sizes (task forms)")
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000, 9001, 9002])
    ap.add_argument("--lmax", type=int, default=None, help="largest |l'| (default p + 1: the window and one beyond)")
    ap.add_argument("--F", type=int, nargs="*", default=None, help="restrict to these charges")
    ap.add_argument("--direct-max", type=int, default=14, help="largest N' with the direct BPS check")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    lmax = a.p + 1 if a.lmax is None else a.lmax
    os.makedirs(OUT, exist_ok=True)
    out = a.out or os.path.join(OUT, f"{a.task}_p{a.p}.jsonl")
    if a.dry_run:
        for Np in a.Np:
            print(f"N'={Np}: estimate {mem_estimate_gb(Np, a.p):.2f} GB; sectors {window_sectors(Np, a.p, lmax)}")
        for N in a.N:
            print(f"N={N}: forms estimate {mem_estimate_gb(N + 1, a.p) / 4:.2f} GB")
        return
    meta = dict(argv=sys.argv, git=git_state(), started=time.strftime("%Y-%m-%dT%H:%M:%S"), out=out)
    with open(os.path.join(OUT, f"meta_{a.task}_p{a.p}_{time.strftime('%Y%m%dT%H%M%S')}_{os.getpid()}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    if a.task == "sectors":
        for Np in a.Np:
            for seed in a.seeds:
                rng = np.random.default_rng(seed)
                c = random_couplings(rng, comb(Np, a.p))
                for F in (a.F if a.F else window_sectors(Np, a.p, lmax)):
                    rec = sector_record(Np, F, c, seed, a.p, a.direct_max, rng)
                    rec["git"] = meta["git"]
                    with open(out, "a") as fh:
                        fh.write(json.dumps(rec) + "\n")
                    ex = rec.get("exact", {})
                    print(f"N'={Np} seed={seed} F={F} l'={rec['lp']:+d} d={rec['d_h']}+{rec['d_g']} (direct {rec.get('d_direct')},"
                          f" expected {rec.get('expected_d')}) ranks {rec['obstruction']['h']['rank']},{rec['obstruction']['g']['rank']}"
                          f" T2={ex.get('T2', float('nan')):.4f} res={rec.get('uplift_residual', float('nan')):.1e}"
                          f" dev={rec.get('decoder_dev', float('nan')):.1e} [{rec['seconds']}s]", flush=True)
    else:
        for N in a.N:
            for seed in a.seeds:
                rec = form_record(N, a.p, seed)
                rec["git"] = meta["git"]
                with open(out, "a") as fh:
                    fh.write(json.dumps(rec) + "\n")
                print(json.dumps({k: v for k, v in rec.items() if k not in ("git",)}), flush=True)


if __name__ == "__main__":
    main()
