"""Var(lambda) of O_T vs Wachter prediction ab(1-b). r := Var/(b(1-b)) is the fraction of the UV
variance of Pi_T that survives BPS projection. Wachter/Haar predicts r = a exactly (disorder avg)."""
import sys; import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np, time
from math import comb
from decoder_fast import decoder_spectrum
print(" N'  P    d      a     b    r=Var/b(1-b)  (sd)   r/a   dm2=Var-ab(1-b)")
for Np in range(7, 15):
    P = Np // 2
    rs = []; t=time.time()
    for s in range(4 if Np < 14 else 2):
        lam, d = decoder_spectrum(Np, P, s)
        D = comb(Np, P); b = comb(Np-1, P)/D; a = d/D
        # use disorder-average mean b (not per-realization mean)
        rs.append(np.mean((lam-b)**2)/(b*(1-b)))
    rs=np.array(rs)
    print(f"{Np:3d} {P:2d} {d:5d} {a:6.3f} {b:5.3f}   {rs.mean():.4f}   ({rs.std():.4f})  {rs.mean()/a:5.2f}   {(rs.mean()-a)*b*(1-b):+.4f}   [{time.time()-t:.0f}s]", flush=True)
