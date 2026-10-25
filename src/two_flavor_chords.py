r"""
two_flavor_chords.py -- double-scaled chord description of the two-flavor N=2 SYK model (CCSY arXiv:2412.06902 eq. 3.18;
HTZ arXiv:2206.14900 eq. 3.1) and of its pair uplift N -> N+1 (docs/derivations.md D9).

Chord dictionary (D9 sec. 1; derived by repeating BBNR arXiv:2003.04405 secs. 4.2-4.5 with two flavors).  Let every term
of Q create p_psi psi's and annihilate p_chi chi's (our model: p_psi = 2, p_chi = 1).  Then the normalized
fixed-charge moments of H in the sector (N_psi, N_chi) of N pairs are those of the one-flavor N=2 chord model (BBNR, BLY)
at
    lambda = 2 (p_psi^2 + p_chi^2) / N ,       q = e^{-lambda} ,
    s = [p_psi (N_psi - N/2) - p_chi (N_chi - N/2)] / (p_psi^2 + p_chi^2)            (Q raises s by 1).
A matter operator with independent Gaussian couplings that creates p'_psi psi's and annihilates p'_chi chi's (chiral:
the orientation of Q in both flavors) raises s by delta and crosses H-chords with weight q^delta,
    delta = (p_psi p'_psi + p_chi p'_chi) / (p_psi^2 + p_chi^2) ,
so the uplift operators E (psi^dag chi), K (psi^dag psi^dag), Z (psi^dag) have delta = 3/5, 4/5, 2/5; their daggers
lower s by the same amounts (the mirror s -> -s).

Consequences (one-flavor results of src/dssyk_n2.py with lambda, s, delta as above):
  * BPS fraction D(s) (BLY eq. 3.17), nonzero only in the black-hole window |s| < 1/2, one sector per Q-complex;
  * continuum edges (BBNR eq. 5.3): the exact (B) family at charge s starts at 2 q^{-1/2}[cosh lam(s - 1/2) - 1], the
    coexact (A) family at 2 q^{-1/2}[cosh lam(s + 1/2) - 1]; the window edges s = +1/2 (B) and s = -1/2 (A) are gapless;
  * BPS -> BPS weight (point mass) and time-integrated BPS -> continuum two-point functions of E, K, Z (D7 secs. 1-2).

Uplift bookkeeping (D9 sec. 3).  For a new BPS state in sector (A, B) at N' = N + 1 with chord charge s' (at N'), the
old slot sectors sit at s' + 1/10 (00), s' - 1/10 (11), s' - 3/10 (10), s' + 3/10 (01), and the six dressing terms of
D8 sec. 7 land at s' + 7/10 (J10), s' - 7/10 (J01), s' +- 9/10 (z11 = -Q0^+ E h01, x00 = -(Q0^dag)^+ E^dag h10) and
s' +- 11/10 (z00 = -Q0^+ K h01, x11 = (Q0^dag)^+ K^dag h10).  A raising term is resonant when it lands on s = +1/2 (gapless
B edge), a lowering term when it lands on s = -1/2 (gapless A edge).

Energy scale: as in D7, kappa = m1 / (q^s + q^{-s}) with m1 the exact mean eigenvalue of the old Laplacian in the
relevant sector (the first chord moment is q^s + q^{-s}).  No fitted parameters.
"""
from __future__ import annotations

from math import cosh, exp

import numpy as np

from dssyk_n2 import bps_bps_weight, bps_continuum_weight, bps_fraction_chord, chord_energy

P_PSI, P_CHI = 2, 1
P2 = P_PSI ** 2 + P_CHI ** 2

# matter operators of the pair uplift (D8 sec. 1): name -> ((d N_psi, d N_chi), (p'_psi, p'_chi))
UPLIFT_OPS = {"E": ((1, -1), (1, 1)), "K": ((2, 0), (2, 0)), "Z": ((1, 0), (1, 0))}


def lam(N: int) -> float:
    """lambda = 2 (p_psi^2 + p_chi^2) / N  (= 10/N for the two-flavor model)."""
    return 2 * P2 / N


def s_eff(a: int, b: int, N: int) -> float:
    """Effective chord charge of sector (N_psi, N_chi) = (a, b) of the N-pair model."""
    return (P_PSI * (a - N / 2) - P_CHI * (b - N / 2)) / P2


def delta(pp_psi: int, pp_chi: int) -> float:
    """Charge shift = crossing exponent of a chiral matter chord with pp_psi psi-creators and pp_chi chi-annihilators."""
    return (P_PSI * pp_psi + P_CHI * pp_chi) / P2


DELTA = {k: delta(*v[1]) for k, v in UPLIFT_OPS.items()}           # E: 3/5, K: 4/5, Z: 2/5


def in_window(a: int, b: int, N: int, tol: float = 1e-9) -> bool:
    return abs(s_eff(a, b, N)) < 0.5 - tol


def bps_fraction(a: int, b: int, N: int) -> float:
    """Double-scaled BPS fraction D(s) (BLY 3.17) of sector (a, b); 0 outside the window (BBNR 5.25-5.28)."""
    s = s_eff(a, b, N)
    return bps_fraction_chord(lam(N), s) if abs(s) < 0.5 else 0.0


def edge_energies(a: int, b: int, N: int, m1: float) -> tuple[float, float]:
    """Chord lower edges (E_B, E_A) of the exact and coexact continua of sector (a, b), in the units of the exact mean
    eigenvalue m1 of the Laplacian of that sector (kappa = m1 / (q^s + q^-s))."""
    L, s = lam(N), s_eff(a, b, N)
    q = exp(-L)
    kappa = m1 / (q ** s + q ** -s)
    return (kappa * 2 * q ** -0.5 * (cosh(L * (s - 0.5)) - 1), kappa * 2 * q ** -0.5 * (cosh(L * (s + 0.5)) - 1))


_GL = {}


def _nodes(n: int):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = ((x + 1) * np.pi / 2, w * np.pi / 2)
    return _GL[n]


def two_point(op: str, sign: int, a: int, b: int, N: int, nphi: int = 2001) -> dict:
    """Chord BPS -> (BPS + continuum) two-point function of uplift operator op in {'E','K','Z'} (sign=+1) or its dagger
    (sign=-1), from the BPS states of sector (a, b) of the N-pair model.  Returns the point mass (BPS -> BPS weight,
    i.e. the obstruction weight ||P_H O h||^2 / ||O h||^2), the continuous mass, and I = int f / Lambda over the continuum
    (multiply by c0 / kappa_in for the time-integrated two-point function).  s_in = +-1/2 exactly gives I = inf."""
    L = lam(N)
    d = DELTA[op]
    s_out = sign * s_eff(a, b, N)                 # lowering operators: mirror s -> -s (D7, g family)
    s_in = s_out + d
    phi, w = _nodes(nphi)
    if abs(s_in - 0.5) < 1e-12:
        I = float("inf")
        f = bps_continuum_weight(L, s_out, d, phi)
    else:
        f = bps_continuum_weight(L, s_out, d, phi)
        I = float(np.sum(w * f / chord_energy(L, s_in, phi)))
    return dict(s_out=s_out, s_in=s_in, delta=d, point_mass=bps_bps_weight(L, s_out, d), cont_mass=float(np.sum(w * f)),
                I=I)


def kappa_in(a_t: int, b_t: int, N: int, m1_target: float, sign: int) -> float:
    """Energy scale of the target sector (a_t, b_t) from its exact mean Laplacian eigenvalue (mirror charge if sign=-1)."""
    q = exp(-lam(N))
    s = sign * s_eff(a_t, b_t, N)
    return m1_target / (q ** s + q ** -s)


# ------------------------------------------------------------------------------------------------- uplift bookkeeping
SLOT_SHIFT = {(0, 0): 1 / 10, (1, 1): -1 / 10, (1, 0): -3 / 10, (0, 1): 3 / 10}     # s_old(slot) - s_new, p = (2, 1)

# dressing terms of D8 sec. 7: name -> (source slot, operator, sign); target charge = s_slot + sign * delta
DRESSING_TERMS = {
    "J10[E h00]": ((0, 0), "E", +1), "J10[K h11]": ((1, 1), "K", +1), "J10[Z h01]": ((0, 1), "Z", +1),
    "J01[K+ h00]": ((0, 0), "K", -1), "J01[E+ h11]": ((1, 1), "E", -1), "J01[Z+ h10]": ((1, 0), "Z", -1),
    "z11[E h01]": ((0, 1), "E", +1), "z00[K h01]": ((0, 1), "K", +1),
    "x00[E+ h10]": ((1, 0), "E", -1), "x11[K+ h10]": ((1, 0), "K", -1),
}


def dressing_targets(A: int, B: int, Np: int) -> dict:
    """Old-model target sector and chord charge of every first-order dressing term of a new BPS state in (A, B), N' = Np.
    'resonant' = the term lands on the gapless edge of its continuum (s = +1/2 for raising, -1/2 for lowering)."""
    N = Np - 1
    out = {}
    for name, (T, op, sg) in DRESSING_TERMS.items():
        a, b = A - T[0], B - T[1]
        (da, db), _ = UPLIFT_OPS[op]
        ta, tb = a + sg * da, b + sg * db
        s_t = s_eff(ta, tb, N)
        out[name] = dict(slot=T, source=(a, b), target=(ta, tb), s_target=s_t,
                         resonant=bool(abs(s_t - 0.5 * sg) < 1e-9))
    return out


def resonant_terms(A: int, B: int, Np: int) -> list[str]:
    return [k for k, v in dressing_targets(A, B, Np).items() if v["resonant"]]
