#!/usr/bin/env python3
r"""
Test of the conjectured size <-> chord-number dictionary  k ~ (p - 1) n  (docs/derivations.md D3).

Motivation (Lin 2022, arXiv:2208.07032, eqs. 53-59): in double-scaled SYK the chord number of the two-sided state
equals its operator size divided by the size of one Hamiltonian term.  In one-flavor N=2 SYK the supercharge terms
psi_I have p fermions, but {psi_I, psi_J^dag} = 0 unless I and J overlap (odd p), so H = {Q, Q^dag} contains at most
2p - 2 fermions (FGMS 2016, Sec. 2), i.e. H is at most (p-1)-body in the U(N) sense.  If one unit of chord number
carries the size of one H, the U(N) size spectrum of the BPS projector should sit at k ~ (p-1) n:
p = 3 -> k = 2n (the empirical match w_{2n} ~ P_n of research/pdfs/progress_2026-09-29.pdf), p = 5 -> k = 4n.

For each realization this script records
  f_k(H) = ||H^{(k)}||_F^2 / ||H||_F^2      (U(N) size spectrum of the Hamiltonian on Lambda^F)
  w_k(P) = ||P^{(k)}||_F^2 / d               (size spectrum of the BPS projector; w_0 = a exactly)
and the BLY chord-number distribution P_n at lambda = 2 p^2 / N', j = 0 (even N', half filling).

Usage: .venv/bin/python scripts/memwatch.py --limit-gb 4 -- .venv/bin/python scripts/check_size_per_chord.py \
           --q 5 --Np 12 --seeds 0-1
Memory: ~ (F + 4) dense D x D complex matrices (N' = 14: ~2 GB).
Output: results/data/size_per_chord_2026-09-30/size_per_chord_q{q}_N{N'}.jsonl (append) + meta json.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from chord_invariants import BPSProjector  # noqa: E402
from dssyk_n2 import hh_length_distribution  # noqa: E402
from q_scan import form_wedge, rand_form  # noqa: E402
from size_decomposition import HopTable, size_weights  # noqa: E402


def hamiltonian(Np, F, q, seed):
    """H = M^dag M on Lambda^F, with M = [Q restricted to Lambda^F ; Q^dag restricted to Lambda^F] exactly as in
    chord_invariants.BPSProjector (same couplings for the same seed)."""
    C = rand_form(seed, Np, q)
    blocks = [form_wedge(C, Np, F, q)]
    if F - q >= 0:
        blocks.append(form_wedge(C, Np, F - q, q).conj().transpose())
    Ms = sp.vstack(blocks).tocsr()
    return (Ms.conj().T @ Ms).toarray()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=3, help="p, the number of fermions in each supercharge term")
    ap.add_argument("--Np", type=int, required=True)
    ap.add_argument("--seeds", default="0-1")
    ap.add_argument("--out", default=os.path.join(ROOT, "results/data/size_per_chord_2026-09-30"))
    a = ap.parse_args()
    lo, _, hi = a.seeds.partition("-")
    seeds = list(range(int(lo), int(hi or lo) + 1))
    os.makedirs(a.out, exist_ok=True)
    try:
        git = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"], capture_output=True,
                               text=True).stdout.strip()
        git += "-dirty" if dirty else ""
    except OSError:
        git = "unknown"
    F = a.Np // 2
    meta = dict(args=vars(a), git=git, started=time.strftime("%Y%m%dT%H%M%S"), F=F,
                description="U(N) size spectra of H and of the BPS projector; test of k ~ (p-1) n")
    with open(os.path.join(a.out, f"meta_q{a.q}_N{a.Np}_{meta['started']}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    table = HopTable(a.Np, F)
    lam = 2 * a.q ** 2 / a.Np
    Pn = hh_length_distribution(lam, 0.0) if a.Np % 2 == 0 else None
    for seed in seeds:
        t0 = time.time()
        H = hamiltonian(a.Np, F, a.q, seed)
        fH, residH = size_weights(H, table)                   # normalized by Tr H; rescale to ||H||^2
        fH = fH / fH.sum()
        del H
        pr = BPSProjector(a.Np, F, q=a.q, seed=seed, method="basis")
        Pm = pr.B @ pr.B.conj().T
        w, resid = size_weights(Pm, table)
        row = dict(q=a.q, Np=a.Np, F=F, seed=seed, a=pr.a, d=pr.d, lam=lam, fH=fH.tolist(), residH=residH,
                   w=w.tolist(), resid=resid, Pn=(Pn[:F + 2].tolist() if Pn is not None else None),
                   secs=time.time() - t0)
        with open(os.path.join(a.out, f"size_per_chord_q{a.q}_N{a.Np}.jsonl"), "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"q={a.q} N'={a.Np} seed={seed} d={pr.d} a={pr.a:.4f} ({row['secs']:.0f}s)", flush=True)
        print("   size of H  f_k : " + " ".join(f"{x:.4f}" for x in fH) + f"   (resid {residH:.1e})", flush=True)
        print("   size of P  w_k : " + " ".join(f"{x:.4f}" for x in w) + f"   (resid {resid:.1e})", flush=True)
        if Pn is not None:
            print(f"   BLY P_n (lam={lam:.3f}): " + " ".join(f"{x:.4f}" for x in Pn[:F + 1]), flush=True)


if __name__ == "__main__":
    main()
