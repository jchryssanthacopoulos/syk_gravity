"""Independent check of the Haar genus expansion claims in gravitational_dual_fragility.tex.

m_k = (1/d) E Tr (U P0 U^dag Pi)^k,  P0 rank d = aD, Pi rank m = bD, U Haar on U(D).
E Tr = sum_{sigma,tau in S_k} Wg(sigma tau^-1; D) m^{c(gamma sigma)} d^{c(tau)}.
Wg obtained exactly by inverting G = sum_s D^{c(s)} s in the centre of C[S_k].
Compare the D^0 and D^-2 coefficients with the paper's Wachter moments and with the
series expansion of the claimed W1(x) = -b(1-a)(1-b) x(x-1)/[(x-l-)(x-l+)]^{5/2}.
"""
import itertools
from collections import Counter, defaultdict
import sympy as sp

a, b, D, x, eps = sp.symbols('a b D x epsilon', positive=True)

def compose(p, q):  # (p o q)(i) = p[q[i]]
    return tuple(p[i] for i in q)

def inverse(p):
    inv = [0] * len(p)
    for i, pi in enumerate(p):
        inv[pi] = i
    return tuple(inv)

def cycle_type(p):
    n, seen, ct = len(p), [False] * len(p), []
    for i in range(n):
        if not seen[i]:
            L, j = 0, i
            while not seen[j]:
                seen[j] = True; j = p[j]; L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))

def weingarten(k):
    perms = list(itertools.permutations(range(k)))
    classes = sorted(set(cycle_type(p) for p in perms))
    rep = {c: next(p for p in perms if cycle_type(p) == c) for c in classes}
    ncyc = {c: len(c) for c in classes}
    # unknown class function w; condition (G * w)(sigma) = delta_{sigma,e}
    wsym = {c: sp.Symbol('w_%s' % '_'.join(map(str, c))) for c in classes}
    eqs = []
    ident = tuple(range(k))
    for c in classes:
        s = rep[c]
        expr = 0
        for t in perms:  # (G*w)(s) = sum_t D^{c(t)} w(t^-1 s)
            expr += D ** len(cycle_type(t)) * wsym[cycle_type(compose(inverse(t), s))]
        eqs.append(sp.Eq(expr, 1 if s == ident else 0))
    sol = sp.solve(eqs, list(wsym.values()), dict=True)[0]
    return {c: sp.factor(sol[wsym[c]]) for c in classes}

def haar_moment(k):
    Wg = weingarten(k)
    perms = list(itertools.permutations(range(k)))
    gamma = tuple((i + 1) % k for i in range(k))
    cnt = Counter()
    for s in perms:
        cgs = len(cycle_type(compose(gamma, s)))
        for t in perms:
            cnt[(cycle_type(compose(s, inverse(t))), cgs, len(cycle_type(t)))] += 1
    m_, d_ = b * D, a * D
    tot = sum(n * Wg[c] * m_ ** i * d_ ** j for (c, i, j), n in cnt.items())
    return sp.simplify(tot / d_)

# claimed objects
lm = (sp.sqrt(a * (1 - b)) - sp.sqrt(b * (1 - a))) ** 2
lp = (sp.sqrt(a * (1 - b)) + sp.sqrt(b * (1 - a))) ** 2
paper_m0 = {1: b, 2: b * (a + b - a * b),
            3: b * (a**2 + b**2 + 3*a*b - 3*a**2*b - 3*a*b**2 + 2*a**2*b**2)}
paper_m1 = {1: 0, 2: -b*(1-a)*(1-b),
            3: b*(1-a)*(1-b)*(10*a*b - 5*a - 5*b + 1),
            4: -5*b*(1-a)*(1-b)*(14*a**2*b**2 - 14*a**2*b + 3*a**2 - 14*a*b**2 + 10*a*b - a + 3*b**2 - b),
            5: 5*b*(1-a)*(1-b)*(84*a**3*b**3 - 126*a**3*b**2 + 56*a**3*b - 7*a**3 - 126*a**2*b**3
                                 + 154*a**2*b**2 - 49*a**2*b + 3*a**2 + 56*a*b**3 - 49*a*b**2 + 8*a*b
                                 - 7*b**3 + 3*b**2)}
# W1 large-x expansion: x -> 1/eps; (x-l-)(x-l+) = x^2 - u x ... use u, v
u = a + b - 2*a*b
v = a*b*(1-a)*(1-b)
sig2 = x**2 - 2*u*x + (u**2 - 4*v)   # (x-l-)(x-l+) since l+- = u +- 2 sqrt(v)
W1 = -b*(1-a)*(1-b) * x*(x-1) / sig2**sp.Rational(5, 2)
W1e = sp.series(W1.subs(x, 1/eps), eps, 0, 8).removeO()
W1_coef = {kk: sp.factor(sp.expand(W1e).coeff(eps, kk + 1)) for kk in range(1, 7)}

if __name__ == '__main__':
    import sys
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for k in range(1, kmax + 1):
        mk = haar_moment(k)
        ser = sp.series(mk.subs(D, 1/eps), eps, 0, 4).removeO()
        c0 = sp.factor(sp.expand(ser).coeff(eps, 0))
        c1 = sp.factor(sp.expand(ser).coeff(eps, 1))
        c2 = sp.factor(sp.expand(ser).coeff(eps, 2))
        line = f"k={k}: odd-order 1/D term = {c1};"
        if k in paper_m0:
            line += f" disk==paper? {sp.simplify(c0 - paper_m0[k]) == 0};"
        if k in paper_m1:
            line += f" genus1==paper? {sp.simplify(c2 - paper_m1[k]) == 0};"
        line += f" genus1==W1 series? {sp.simplify(c2 - W1_coef[k]) == 0}"
        print(line, flush=True)
        if k == 2:
            print("   exact m2 =", sp.factor(mk))
            print("   at a=b=1/2:", sp.series(mk.subs({a: sp.Rational(1, 2), b: sp.Rational(1, 2)}).subs(D, 1/eps), eps, 0, 3))
