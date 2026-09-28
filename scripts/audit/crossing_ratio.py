"""Chord crossing weight of BPS-projected occupations: X = tau(AiAjAiAj)/tau(AiAiAjAj), A = P(n-<n>)P, i != j.
X->0: free (non-crossing, Wachter-like); X->1: light probe (crossings unsuppressed). Compare SYK vs Haar P_B."""
import sys; import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np
from decoder_fast import couplings, bps_basis, subset_index
rng=np.random.default_rng(7)
def X_of(B, masks, Np, pairs):
    out=[]
    for i,j in pairs:
        ni=np.array([float((m>>i)&1) for m in masks]); nj=np.array([float((m>>j)&1) for m in masks])
        Ai=B.conj().T@(ni[:,None]*B); Aj=B.conj().T@(nj[:,None]*B); d=B.shape[1]
        Ai-=np.trace(Ai).real/d*np.eye(d); Aj-=np.trace(Aj).real/d*np.eye(d)
        cr=np.trace(Ai@Aj@Ai@Aj).real; un=np.trace(Ai@Ai@Aj@Aj).real
        out.append(cr/un)
    return np.mean(out)
print(" N'  P    d     a    X_SYK   X_Haar(same D,d)")
for Np in range(8,14):
    P=Np//2; masks,_=subset_index(Np,P); D=len(masks)
    pairs=[(Np-1,j) for j in range(0,Np-1,max(1,(Np-1)//4))][:4]
    xs=[];xh=[]
    for s in range(3):
        B=bps_basis(couplings(s,Np),Np,P); d=B.shape[1]
        xs.append(X_of(B,masks,Np,pairs))
        G=rng.standard_normal((D,d))+1j*rng.standard_normal((D,d)); H,_=np.linalg.qr(G)
        xh.append(X_of(H,masks,Np,pairs))
    print(f"{Np:3d} {P:2d} {d:5d} {d/D:.3f}   {np.mean(xs):+.3f}   {np.mean(xh):+.3f}", flush=True)
