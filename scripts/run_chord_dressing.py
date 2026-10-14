#!/usr/bin/env python3
r"""
Finite-N test of the chord (double-scaled) description of the first-order dressing spectrum (docs/derivations.md D7).

For each realization (seed s: couplings exactly as results/data/uplift_2026-10-12 / uplift.OldModel) and each family
(h: B^F -> exact forms of degree F+2; g: B^{F-1} -> coexact forms of degree F-3):
  own        first-order spectrum sigma^2 for the realization's own eta (= the 'first' data of run_uplift_walls.py)
  sampled    own + (neta - 1) fresh etas
  kron       X = S^1/2 Z C^1/2, exact S = G^1/2 V G^1/2 and C = Dbar/t1  (planar law, exact second moment)
  wishart    same with C = 1                                              (planar law, factorized BPS faces)
  chord      chord population of S (closed form, D7), C = 1               (gravitational prediction)
  chord_resc chord population rescaled to the exact t1
plus t1 (exact vs chord), Dbar heterogeneity, diagonal fraction of S in the energy basis, the second-moment identity
and the third moment (sampled vs exact non-crossing terms).  W = <y(1-y)>, T2 = <y>, y = ((1-s)/(1+s))^2.

Usage: .venv/bin/python scripts/run_chord_dressing.py --Np 14 --seeds 9000 9001 9002 9003 --neta 24
       .venv/bin/python scripts/run_chord_dressing.py --Np 17 --dry-run          # memory estimate only
Run heavy cases (N' >= 15) under scripts/memwatch.py.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_dressing import DressingFamily, analyze_family, memory_estimate_gb, old_model_lean  # noqa: E402

OUT = os.path.join(ROOT, "results/data/chord_dressing_2026-10-14")


def git_state():
    try:
        g = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
        return g + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[9000])
    ap.add_argument("--neta", type=int, default=24, help="number of eta samples (the first is the realization's own)")
    ap.add_argument("--draws", type=int, default=3, help="Monte Carlo draws of each separable model")
    ap.add_argument("--no-third-moment", action="store_true", help="skip the exact non-crossing third-moment terms")
    ap.add_argument("--families", nargs="+", default=["h", "g"])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tag", default="", help="suffix of the output file name")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    est = memory_estimate_gb(a.Np)
    print(f"N'={a.Np}: estimated peak memory {est['total']:.1f} GB (basis {est['basis']:.1f}, held {est['held']:.1f}, "
          f"analysis {est['analysis']:.1f})", flush=True)
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"decomp_N{a.Np}{a.tag}.jsonl")
    meta = dict(args=vars(a), git=git_state(), memory_estimate_gb=est, started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    with open(os.path.join(a.out, f"meta_decomp_N{a.Np}{a.tag}_{time.strftime('%Y%m%dT%H%M%S')}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    for seed in a.seeds:
        t0 = time.time()
        om = old_model_lean(a.Np, seed)
        rng = np.random.default_rng(seed + 777)
        rec = dict(Np=a.Np, seed=seed)
        for fam in a.families:
            Fm = DressingFamily(om, fam)
            eta_own = om["eta"]
            rec[fam], _ = analyze_family(Fm, eta_own, a.neta, rng, draws=a.draws, third_moment=not a.no_third_moment,
                                         Np=a.Np)
            del Fm
        if set(a.families) == {"h", "g"}:
            bh, bg = rec["h"]["b"], rec["g"]["b"]
            for lev in ("own", "sampled", "kron", "wishart", "chord", "chord_rescaled"):
                rec[lev] = {k: (bh * rec["h"][lev][k] + bg * rec["g"][lev][k]) / (bh + bg) for k in ("W", "T2")}
        rec["secs"] = time.time() - t0
        with open(path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        if "own" in rec:
            print(f"N'={a.Np} seed={seed} [{rec['secs']:.0f}s] W: " + " ".join(f"{k} {rec[k]['W']:.4f}" for k in
                  ("own", "sampled", "kron", "wishart", "chord", "chord_rescaled")), flush=True)
        for fam in a.families:
            r = rec[fam]
            print(f"   {fam}: b={r['b']} e={r['e']} t1={r['t1']:.4f} (chord {r['t1_chord']:.4f}) Dbar rsd {r['Dbar_rel_sd']:.3f} "
                  f"diagS {r['diag_frac_S']:.4f} m2 identity {r['m2_identity']:.4f} sampled {r['m2_sampled']:.4f}", flush=True)


if __name__ == "__main__":
    main()
