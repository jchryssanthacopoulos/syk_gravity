#!/usr/bin/env python3
r"""
Selection-rule check: maximal number of rows of the shapes in Sym^n(Lambda^p) (plethysm h_n[e_p]).

The chaos-n (chord-number-n) component of the BPS projector transforms under U(N) inside
Sym^n(Lambda^p) (x) conj(Sym^n(Lambda^p)), and the size-k irrep V_k = [1^k; 1^k] occurs there only if k <= l(mu)
for a shape mu of h_n[e_p].  For odd p the prediction is max l(mu) = n (p - 1) + 1
(research/notes/chord_size_dictionary_2026-10-02.md, Sec. 3).

Usage: .venv/bin/python scripts/check_plethysm_rows.py --cases 2,3 3,3 2,5 --out results/data/chord_chaos_2026-10-01/checks/plethysm_rows.txt
Cost: m! antisymmetrization per shape with m = n p variables (m <= 9: seconds).
"""
import argparse
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from plethysm import max_rows  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", nargs="+", default=["2,3", "3,3", "2,5"], help="n,p pairs (n p <= 9)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    lines = []
    for case in a.cases:
        n, p = (int(x) for x in case.split(","))
        mr, res = max_rows(n, p)
        pred = n * (p - 1) + 1 if p % 2 else n * p
        lines.append(f"h_{n}[e_{p}] = " + " + ".join(f"{c} s{lam}" if c != 1 else f"s{lam}" for lam, c in res.items())
                     + f"   max rows = {mr}   prediction n(p-1)+1 = {pred}   {'OK' if mr == pred else 'MISMATCH'}")
    txt = "\n".join(lines)
    print(txt)
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write(txt + "\n")


if __name__ == "__main__":
    main()
