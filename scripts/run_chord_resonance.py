#!/usr/bin/env python3
r"""
Layer 1 of the resonant chord problem (docs/derivations.md D10 §4 -> D11): the first-order dressing spectrum of the
resonant family from the chord rules (heavy-tailed free compound Poisson), against finite N.

For each N', seed and family, the D7 error decomposition (src/chord_dressing.analyze_family):
  own      exact first-order spectrum sigma^2 for the realization's own eta (= sigma2_h of the window runs)
  sampled  own + (neta - 1) fresh etas (the eta ensemble)
  kron     separable model X = S^1/2 Z C^1/2 with the exact inner S and C = Dbar/t1
  wishart  exact inner S, C = 1 (planar law, factorized BPS faces)
  chord    chord population of S (D7 jump law at the sector's (s_out, s_in)), C = 1, at the exact (e, b)
Resonant sector: l' = -2 (F = (N'-2)/2), h family, s_in = 1/2, where the jump law ~ 1/phi^2.
Control: the central sector l' = 0 (F = N'/2), both families, which reproduces D7.
Inner spectra come from SVDs (DressingFamily(inner='svd')); the Gram cut of D7 would drop the near-zero edge modes.
Couplings: rng = default_rng(seed), c = random_couplings(rng, C(N', 3)), as in results/data/uplift_window_2026-10-26.
Admissible sources from uplift.OldModel(F=...).

Recorded per (N', seed, l', family): b, e, inner E_min/E_mean, t1 (exact) and t1 (chord; inf if resonant).
For each layer: W = <y(1-y)>, T2 = <y> (y = ((1-s)/(1+s))^2), m1, atoms, quantiles, P(s > K) sqrt(K) at
K = 10, 100, 1000, and Hill exponents.  The own/sampled/wishart/chord spectra go to an .npz per record.

Usage (under scripts/memwatch.py):
  .venv/bin/python scripts/run_chord_resonance.py --Np 8 10 12 14 --seeds 9000 9001 9002 --neta 24
  .venv/bin/python scripts/run_chord_resonance.py --Np 16 --seeds 9000 9001 --neta 8 --no-control
  --dry-run prints a memory estimate (N' = 16: about 6 GB).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from math import comb

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_chaos import random_couplings  # noqa: E402
from chord_dressing import DressingFamily, analyze_family, old_model_from_uplift  # noqa: E402
from uplift import OldModel, level  # noqa: E402

OUT = os.path.join(ROOT, "results/data/chord_resonance_2026-10-26")


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", "src", "scripts"],
                               capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def mem_estimate_gb(Np):
    N = Np - 1
    D = max(comb(N, k) for k in range(N + 1))
    return 16 * D * D * 16 / 1e9 + 0.3            # OldModel dense maps + P, V, G, lift scratch (all ~D^2)


def tail_stats(x):
    x = np.asarray(x, float)
    out = dict(n=len(x), median=float(np.median(x)), q90=float(np.quantile(x, 0.9)), q99=float(np.quantile(x, 0.99)),
               max=float(x.max()))
    for K in (10.0, 100.0, 1000.0):
        out[f"Psqrt_{int(K)}"] = float(np.mean(x > K) * np.sqrt(K))
    for K in (10.0, 100.0):
        t = x[x > K]
        out[f"hill_{int(K)}"] = float(len(t) / np.sum(np.log(t / K))) if len(t) > 2 else float("nan")
        out[f"hill_{int(K)}_n"] = int(len(t))
    return out


def run_one(Np, seed, F, fam, neta, draws, third, rng_seed):
    rng = np.random.default_rng(seed)
    c = random_couplings(rng, comb(Np, 3))
    om = old_model_from_uplift(OldModel(Np, c, rng, F=F))
    fobj = DressingFamily(om, fam, inner="svd")
    rec, spec = analyze_family(fobj, om["eta"], neta, np.random.default_rng(rng_seed), draws=draws, third_moment=third,
                               Np=Np, Fsec=F)
    rec.update(Np=Np, N=Np - 1, F=F, lp=level(Np, F), fam=fam, seed=seed, Emin_rel=float(fobj.we.min() / fobj.we.mean()))
    for k, x in spec.items():
        rec[k].update(tail_stats(x))
    rec.pop("Dbar_eigs", None)
    return rec, spec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000, 9001, 9002])
    ap.add_argument("--neta", type=int, default=24, help="etas per realization (own + fresh)")
    ap.add_argument("--draws", type=int, default=3, help="separable-model draws")
    ap.add_argument("--no-control", action="store_true", help="skip the central control sector l' = 0")
    ap.add_argument("--third-moment", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.dry_run:
        for Np in a.Np:
            print(f"N'={Np}: estimate {mem_estimate_gb(Np):.2f} GB")
        return
    os.makedirs(OUT, exist_ok=True)
    out = a.out or os.path.join(OUT, "chord_resonance.jsonl")
    meta = dict(argv=sys.argv, git=git_state(), started=time.strftime("%Y-%m-%dT%H:%M:%S"), out=out)
    with open(os.path.join(OUT, f"meta_{time.strftime('%Y%m%dT%H%M%S')}_{os.getpid()}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    for Np in a.Np:
        jobs = [((Np - 2) // 2, "h")]                                    # resonant: l' = -2, h family
        if not a.no_control and Np % 2 == 0:
            jobs += [(Np // 2, "h"), (Np // 2, "g")]                     # control: l' = 0 (D7)
        for seed in a.seeds:
            for F, fam in jobs:
                t0 = time.time()
                rec, spec = run_one(Np, seed, F, fam, a.neta, a.draws, a.third_moment, rng_seed=10 * seed + F)
                rec["seconds"], rec["git"] = round(time.time() - t0, 1), meta["git"]
                with open(out, "a") as fh:
                    fh.write(json.dumps(rec) + "\n")
                np.savez_compressed(os.path.join(OUT, f"spec_N{Np}_F{F}_{fam}_s{seed}.npz"),
                                    **{k: spec[k] for k in ("own", "sampled", "wishart", "chord") if k in spec})
                print(f"N'={Np} seed={seed} l'={rec['lp']:+d} {fam}: b={rec['b']} e={rec['e']} Emin/Emean={rec['Emin_rel']:.1e}"
                      f" | W own {rec['own']['W']:.4f} sampled {rec['sampled']['W']:.4f} wishart {rec['wishart']['W']:.4f}"
                      f" chord {rec['chord']['W']:.4f} | T2 own {rec['own']['T2']:.4f} sampled {rec['sampled']['T2']:.4f}"
                      f" chord {rec['chord']['T2']:.4f} | P*sqrtK(100) sampled {rec['sampled']['Psqrt_100']:.3f}"
                      f" chord {rec['chord']['Psqrt_100']:.3f} [{rec['seconds']}s]", flush=True)


if __name__ == "__main__":
    main()
