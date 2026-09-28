"""Exploratory: principal angles between the enlarged BPS space B_{N'}^P and the
LES reference subspace B0 = B_N^P (mode N' empty)  (+)  B_N^{P-1} ^ e_{N'} (occupied).
Also: are small-lambda decoder eigenvectors concentrated in the occupied copy of B_N^{P-1}?"""
import sys; import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np
from decoder_fast import couplings, bps_basis, subset_index
for N in [9, 10, 11]:
    Np, P = N + 1, (N + 1) // 2
    top = 1 << N
    for seed in range(2):
        C = couplings(seed, Np)
        Cold = {k: v for k, v in C.items() if not (k & top)}
        B = bps_basis(C, Np, P)                        # enlarged BPS basis (C(N',P) x d)
        Be = bps_basis(Cold, N, P); Bo = bps_basis(Cold, N, P - 1)
        masks, idx = subset_index(Np, P)
        mP, _ = subset_index(N, P); mPm, _ = subset_index(N, P - 1)
        E0 = np.zeros((len(masks), Be.shape[1] + Bo.shape[1]), complex)
        for j, mk in enumerate(mP):  E0[idx[mk], :Be.shape[1]] = Be[j]
        for j, mk in enumerate(mPm): E0[idx[mk | top], Be.shape[1]:] = Bo[j]
        d, d0 = B.shape[1], E0.shape[1]
        cos = np.linalg.svd(B.conj().T @ E0, compute_uv=False)
        # decoder eigenvectors
        empty = np.array([0. if (mk & top) else 1. for mk in masks])
        O = B.conj().T @ (empty[:, None] * B); lam, V = np.linalg.eigh((O + O.conj().T) / 2)
        Psi = B @ V
        Pocc = E0[:, Be.shape[1]:]; Pemp = E0[:, :Be.shape[1]]
        w_occ = np.sum(np.abs(Pocc.conj().T @ Psi) ** 2, axis=0)   # weight in occupied copy of B_N^{P-1}
        w_emp = np.sum(np.abs(Pemp.conj().T @ Psi) ** 2, axis=0)
        lo, hi = lam < 0.05, lam > 0.95
        print(f"N={N}->{Np} P={P} seed={seed}: d={d} dimB0={d0} | principal cos^2: "
              f"mean={np.mean(cos**2):.3f}, frac>0.9={np.mean(cos**2>0.9):.3f}, frac<0.1={np.mean(cos**2<0.1):.3f} | "
              f"lam<.05: n={lo.sum()}, mean w_occ={w_occ[lo].mean() if lo.any() else np.nan:.3f} | "
              f"lam>.95: n={hi.sum()}, mean w_emp={w_emp[hi].mean() if hi.any() else np.nan:.3f} | "
              f"bulk .2<lam<.8: w_occ+w_emp={np.mean((w_occ+w_emp)[(lam>.2)&(lam<.8)]):.3f}")
