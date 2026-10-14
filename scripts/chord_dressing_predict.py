#!/usr/bin/env python3
r"""
Gravitational (double-scaled chord) prediction of the first-order dressing spectrum, docs/derivations.md D7:
free compound Poisson law with rate e/b and the chord jump law, b -> infinity (Monte Carlo at large effective b),
for both families; t1 from the 2phi1 closed form.  No input from finite-N diagonalization.

Usage: .venv/bin/python scripts/chord_dressing_predict.py --Np 10 11 12 13 14 15 16 17 18 20 24 30 40 --json out.json
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from chord_dressing import chord_prediction  # noqa: E402
from run_chord_dressing import OUT, git_state  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, nargs="+", required=True)
    ap.add_argument("--b-eff", type=int, default=1500)
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--json", default=os.path.join(OUT, "chord_predictions.json"))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    res = dict(git=git_state(), args=vars(a), predictions={})
    for Np in a.Np:
        p = chord_prediction(Np, rng, b_eff=a.b_eff, draws=a.draws)
        res["predictions"][str(Np)] = p
        h, g = p["h"], p["g"]
        print(f"N'={Np}: W1 = {p['W']:.4f}  T2_1 = {p['T2']:.4f} | h: r={h['r']:.3g} t1={h['t1_closed']:.4f} "
              f"sqrt(t2)/t1={np.sqrt(h['t2']) / h['t1']:.3f} | g: r={g['r']:.3g} t1={g['t1_closed']:.4f} "
              f"sqrt(t2)/t1={np.sqrt(g['t2']) / g['t1']:.3f}", flush=True)
    os.makedirs(os.path.dirname(a.json), exist_ok=True)
    with open(a.json, "w") as fh:
        json.dump(res, fh, indent=1)
    print("wrote", a.json)


if __name__ == "__main__":
    main()
