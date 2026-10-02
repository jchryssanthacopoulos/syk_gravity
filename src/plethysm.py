r"""
plethysm.py -- Schur expansion of the plethysm h_n[e_p], i.e. the GL(m) decomposition of Sym^n(Lambda^p C^m).

Used for the chord-number / operator-size selection rule (research/notes/chord_size_dictionary_2026-10-02.md, Sec. 3):
the n-th Wiener chaos of the couplings carries Sym^n(Lambda^p) (x) conj(Sym^n(Lambda^p)); the U(N) size-k irrep
V_k = [1^k; 1^k] occurs only if k <= l(mu) for some mu in h_n[e_p].  For odd p, max l(mu) = n (p - 1) + 1.

Method: the character of Sym^n(Lambda^p C^m) is a sum of monomials x^e over multisets of n p-subsets; the Schur
coefficient of lambda is the coefficient of x^{lambda + delta} in (character) * a_delta, computed by antisymmetrizing
over S_m.  Exact for m >= n p (enough variables for every shape); cost ~ m! per shape, fine for m <= 9.
"""
from __future__ import annotations

import itertools
from collections import Counter


def partitions(n, maxpart, maxlen):
    """Partitions of n with parts <= maxpart and at most maxlen parts (weakly decreasing tuples)."""
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k, maxlen - 1):
            yield (k,) + rest


def _perm_sign(pm):
    m, s, seen = len(pm), 1, [False] * len(pm)
    for i in range(m):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = pm[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def plethysm_h_e(n, p, m):
    """Dict {shape: multiplicity} of h_n[e_p] in m variables (shapes with at most n columns and n p boxes)."""
    subsets = list(itertools.combinations(range(m), p))
    f = Counter()
    for ms in itertools.combinations_with_replacement(range(len(subsets)), n):
        e = [0] * m
        for s in ms:
            for i in subsets[s]:
                e[i] += 1
        f[tuple(e)] += 1
    delta = list(range(m - 1, -1, -1))
    perms = list(itertools.permutations(range(m)))
    signs = [_perm_sign(pm) for pm in perms]
    out = {}
    for lam in partitions(n * p, n, m):
        lam = lam + (0,) * (m - len(lam))
        tgt = [lam[i] + delta[i] for i in range(m)]
        c = 0
        for pm, sg in zip(perms, signs):
            e = tuple(tgt[i] - delta[pm[i]] for i in range(m))
            if min(e) < 0:
                continue
            c += sg * f.get(e, 0)
        if c:
            out[tuple(x for x in lam if x)] = c
    return out


def max_rows(n, p, m=None):
    """Largest number of rows among the shapes of h_n[e_p] (m defaults to n p)."""
    res = plethysm_h_e(n, p, m if m is not None else n * p)
    return max(len(lam) for lam in res), res
