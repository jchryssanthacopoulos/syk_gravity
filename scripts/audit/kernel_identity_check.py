"""Is  Tr(O^k) = d*m_k^W(a,b) + sum_l K_{k,l} T_l  an identity for ANY pair of projectors?
K_{k,l} = sum_i C(k,l+2i) C(l+2i,i) u^{k-l-2i} v^i ,  T_l = Tr[(P-a)(Pi-b)]^l.
Test on deliberately NON-free pairs (random coordinate projectors, nested, nearly aligned)."""
import numpy as np
from math import comb
rng = np.random.default_rng(0)

def wachter_moment(a, b, k):   # per-rank-d normalisation, via Lagrange inversion of S-transform
    num = np.array([a*b, a+b, 1.0]); inv = np.array([(-1.0)**j for j in range(k+1)])
    h = np.convolve(num, inv)[:k+1]; p = np.array([1.0])
    for _ in range(k): p = np.convolve(p, h)[:k+1]
    return p[k-1] / k / a

def K(k, l, u, v):
    return sum(comb(k, l+2*i)*comb(l+2*i, i)*u**(k-l-2*i)*v**i for i in range(0, (k-l)//2+1))

def proj(X):
    Q, _ = np.linalg.qr(X); return Q @ Q.conj().T

D = 60
cases = {}
d, m = 25, 33
cases['haar'] = (proj(rng.normal(size=(D, d)) + 1j*rng.normal(size=(D, d))), np.diag([1.]*m + [0.]*(D-m)))
cases['coord-coord'] = (np.diag(rng.permutation([1.]*d + [0.]*(D-d))), np.diag([1.]*m + [0.]*(D-m)))
A = rng.normal(size=(D, d)); A[:m//2] *= 5          # strongly aligned with the slot
cases['aligned'] = (proj(A), np.diag([1.]*m + [0.]*(D-m)))
for name, (P, Pi) in cases.items():
    a, b = np.trace(P).real/D, np.trace(Pi).real/D
    u, v = a+b-2*a*b, a*b*(1-a)*(1-b)
    M = (P - a*np.eye(D)) @ (Pi - b*np.eye(D))
    O = P @ Pi @ P
    errs = []
    for k in range(1, 8):
        lhs = np.trace(np.linalg.matrix_power(O, k)).real
        T = [np.trace(np.linalg.matrix_power(M, l)).real for l in range(1, k+1)]
        rhs = a*D*wachter_moment(a, b, k) + sum(K(k, l, u, v)*T[l-1] for l in range(1, k+1))
        rhs_paper = D*wachter_moment(a, b, k) + sum(K(k, l, u, v)*T[l-1] for l in range(1, k+1))
        errs.append((lhs - rhs, lhs - rhs_paper))
    print(f"{name:12s} max|identity err| (d m^W) = {max(abs(e[0]) for e in errs):.2e};"
          f"  with literal 'D m^W' of eq.(errdef) = {max(abs(e[1]) for e in errs):.2e}")
