#!/usr/bin/env python3
"""End-to-end check of the BHKK narrow sieve for p-packing of q-sets (arXiv:1007.1161, Section 4).
Python 3, standard library only; about 25 seconds.

For a set family F, an (n1, n-n1)-partition and parameters (p0,p1,p2), compute at a random point of GF(2^16):
  sieve : sum_{J subset L} [w0^{2p0} w1^{2p1} w2^{2p2}] S(J) det T(J)      (their (20), (21)-(24))
          with det T evaluated numerically and the w-coefficients recovered by interpolation
  fm    : sum over (p0,p1,p2)-packings A, compatible involutions iota, bijective labelings lambda
          of the monomial (18) with the w's stripped                       (their Lemma 13 / Section 4.6)
and compare with brute force: does F contain a (p0,p1,p2)-packing for this partition?
"""
import itertools, random, sys, time
from math import comb

Q = 65535; POLY = 0x1100B
EXP = [0] * (2 * Q); LOG = [0] * (Q + 1); v = 1
for i in range(Q):
    EXP[i] = v; LOG[v] = i; v <<= 1
    if v & 0x10000: v ^= POLY
assert v == 1 and len(set(EXP[:Q])) == Q
for i in range(Q, 2 * Q): EXP[i] = EXP[i - Q]
def mul(a, b): return EXP[LOG[a] + LOG[b]] if a and b else 0
def inv(a): return EXP[Q - LOG[a]]
def sq(a): return mul(a, a)

def det(M):
    n = len(M)
    if n == 0: return 1
    M = [r[:] for r in M]; d = 1
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c]), None)
        if p is None: return 0
        M[c], M[p] = M[p], M[c]; piv = M[c][c]; d = mul(d, piv); ip = inv(piv); Rc = M[c]
        for i in range(c + 1, n):
            if M[i][c]:
                f = mul(M[i][c], ip); Ri = M[i]
                for j in range(c, n):
                    if Rc[j]: Ri[j] ^= mul(f, Rc[j])
    return d

def interp_coeff(xs, ys, t):
    N = [1]
    for a in xs:
        R = [0] * (len(N) + 1)
        for i, c in enumerate(N): R[i + 1] ^= c; R[i] ^= mul(c, a)
        N = R
    res = 0
    for xi, yi in zip(xs, ys):
        deg = len(N) - 1; Qd = [0] * deg; Qd[deg - 1] = N[deg]
        for k in range(deg - 1, 0, -1): Qd[k - 1] = N[k] ^ mul(xi, Qd[k])
        den = 0; pw = 1
        for c in Qd: den ^= mul(c, pw); pw = mul(pw, xi)
        res ^= mul(yi, mul(Qd[t], inv(den))) if t < len(Qd) else 0
    return res

class Vars:
    def __init__(self, rng): self.rng, self.val = rng, {}
    def __call__(self, *key):
        if key not in self.val: self.val[key] = self.rng.randrange(Q + 1)
        return self.val[key]

def sieve(F, U1, q, p0, p1, p2, var):
    U1 = list(U1); U1s = set(U1); n1 = len(U1); idx = {u: i for i, u in enumerate(U1)}
    Lsize = q * p0 + (q - 1) * p1 + (q - 2) * p2
    total = 0
    alphas = [EXP[i] for i in range(1, n1 + 2)]          # alpha stands for w1^2, degree <= n1
    betas = [EXP[i] for i in range(1, n1 + 2)]           # w2, degree <= n1
    for J in range(1 << Lsize):
        keep = [l for l in range(Lsize) if not J >> l & 1]
        def zsum(u, squared):
            s = 0
            for l in keep: s ^= sq(var('z', u, l)) if squared else var('z', u, l)
            return s
        # S(J): coefficient of w0^{2 p0} = elementary symmetric e_{p0} of a_A over sets missing U1
        e = [1] + [0] * p0
        for A in F:
            if A & U1s: continue
            a = sq(var('x', A))
            for u in A: a = mul(a, zsum(u, True))
            for k in range(p0, 0, -1): e[k] ^= mul(e[k - 1], a)
        sJ = e[p0]
        if not sJ: continue
        # det T(J) as a polynomial in alpha = w1^2 and w2; extract [alpha^{p1} w2^{2 p2}]
        diag_sum = [0] * n1; off = {}
        for A in F:
            I = A & U1s
            if len(I) == 1:
                (u,) = I; t = sq(var('x', A))
                for u2 in A - U1s: t = mul(t, zsum(u2, True))
                diag_sum[idx[u]] ^= t
            elif len(I) == 2:
                u, w = sorted(I); t = var('x', A)
                for u2 in A - U1s: t = mul(t, zsum(u2, False))
                off[(u, w)] = off.get((u, w), 0) ^ t
        cs = []
        for al in alphas:
            ys = []
            for be in betas:
                M = [[0] * n1 for _ in range(n1)]
                for i, u in enumerate(U1):
                    M[i][i] = mul(sq(var('y', u)), 1 ^ mul(al, diag_sum[i]))
                for i in range(n1):
                    for j in range(i + 1, n1):
                        u, w = U1[i], U1[j]
                        val = mul(var('y', u, w), 1 ^ mul(be, off.get((u, w), 0)))
                        M[i][j] = M[j][i] = val
                ys.append(det(M))
            cs.append(interp_coeff(betas, ys, 2 * p2))
        tJ = interp_coeff(alphas, cs, p1)
        total ^= mul(sJ, tJ)
    return total

def involutions(points):
    points = list(points)
    if not points: yield []; return
    a, rest = points[0], points[1:]
    for inv_ in involutions(rest): yield [(a, a)] + inv_
    for i, b in enumerate(rest):
        for inv_ in involutions(rest[:i] + rest[i + 1:]): yield [(a, b)] + inv_

def formula(F, U1, q, p, p0, p1, p2, var):
    U1s = set(U1); Lsize = q * p0 + (q - 1) * p1 + (q - 2) * p2; total = 0
    for A in itertools.combinations(F, p):
        if any(a & b for a, b in itertools.combinations(A, 2)): continue
        k = [len(S & U1s) for S in A]
        if max(k) > 2 or (k.count(0), k.count(1), k.count(2)) != (p0, p1, p2): continue
        covered = set().union(*(S & U1s for S in A)); fixed_pairs = []
        for S in A:
            I = sorted(S & U1s)
            if len(I) == 1: fixed_pairs.append((I[0], I[0]))
            if len(I) == 2: fixed_pairs.append((I[0], I[1]))
        base = 1
        for S in A: base = mul(base, sq(var('x', S)))
        for (a, b) in fixed_pairs: base = mul(base, sq(var('y', a) if a == b else var('y', a, b)))
        ysum = 0
        for iv in involutions(sorted(U1s - covered)):
            t = 1
            for (a, b) in iv: t = mul(t, sq(var('y', a) if a == b else var('y', a, b)))
            ysum ^= t
        dom = [(u, S) for S in A for u in sorted(S - U1s)]
        assert len(dom) == Lsize
        lsum = 0
        for perm in itertools.permutations(range(Lsize)):
            t = 1
            for (u, S), l in zip(dom, perm): t = mul(t, sq(var('z', u, l)))
            lsum ^= t
        total ^= mul(base, mul(ysum, lsum))
    return total

def has_pattern(F, U1, p, p0, p1, p2):
    U1s = set(U1)
    for A in itertools.combinations(F, p):
        if any(a & b for a, b in itertools.combinations(A, 2)): continue
        k = [len(S & U1s) for S in A]
        if max(k) <= 2 and (k.count(0), k.count(1), k.count(2)) == (p0, p1, p2): return True
    return False

def triples(p):
    return [(p - a - b, a, b) for a in range(p + 1) for b in range(p + 1 - a)]

def run(name, n, q, p, F, rng, extra=(0, 1), sample=None):
    stats = dict(cases=0, match_formula=0, mismatch=0, sound_zero=0, detected=0, sz_miss=0, false_pos=0)
    t0 = time.time()
    for (p0, p1, p2) in triples(p):
        for ex in extra:
            n1 = p1 + 2 * p2 + ex
            if n1 > n or n1 - p1 - 2 * p2 > n - p * q: continue
            parts = list(itertools.combinations(range(n), n1))
            if sample and len(parts) > sample: parts = rng.sample(parts, sample)
            for U1 in parts:
                var = Vars(rng)
                s = sieve(F, U1, q, p0, p1, p2, var); f = formula(F, U1, q, p, p0, p1, p2, var)
                truth = has_pattern(F, U1, p, p0, p1, p2)
                stats['cases'] += 1
                stats['match_formula' if s == f else 'mismatch'] += 1
                if truth: stats['detected' if s else 'sz_miss'] += 1
                else: stats['false_pos' if s else 'sound_zero'] += 1
    print(f"{name}: n={n}, q={q}, p={p}, |F|={len(F)}, {time.time()-t0:.0f}s")
    print("   " + ", ".join(f"{k}={v}" for k, v in stats.items()))

def rate_check(n, q, p, F, rng):
    """Unique p-packing: the partitions on which the sieve's polynomial is nonzero must be exactly
    as many as Lemma 12 predicts. A zero at the random point on a partition that has the pattern is
    re-tested at fresh points; a nonzero value there is a Schwartz-Zippel miss, not an error."""
    print(f"Lemma 12 end to end: n={n}, q={q}, p={p}, |F|={len(F)} (unique packing)")
    for (p0, p1, p2) in triples(p):
        n1 = p1 + 2 * p2 + (n - p * q) // 2
        if not (0 <= n1 - p1 - 2 * p2 <= n - p * q): continue
        fired = misses = wrong = 0; parts = list(itertools.combinations(range(n), n1))
        for U1 in parts:
            if sieve(F, U1, q, p0, p1, p2, Vars(rng)): fired += 1
            elif has_pattern(F, U1, p, p0, p1, p2):
                if any(sieve(F, U1, q, p0, p1, p2, Vars(random.Random(s))) for s in range(3)): misses += 1
                else: wrong += 1
        exact = comb(p, p1 + p2) * comb(p1 + p2, p2) * q**p1 * comb(q, 2)**p2 * comb(n - p*q, n1 - p1 - 2*p2)
        note = f" (+{misses} random-point miss, nonzero on retry)" if misses else ""
        flag = "" if fired + misses == exact and not wrong else "   <-- MISMATCH"
        print(f"   (p0,p1,p2)=({p0},{p1},{p2}) n1={n1}: fired on {fired}{note} of {len(parts)}; Lemma 12 predicts {exact}{flag}")

def count_packings(F, p):
    return sum(1 for A in itertools.combinations(F, p) if not any(a & b for a, b in itertools.combinations(A, 2)))

if __name__ == "__main__":
    rng = random.Random(2026)
    def rand_family(n, q, m, planted, seed):
        r = random.Random(seed); F = set(frozenset(s) for s in planted)
        while len(F) < m: F.add(frozenset(r.sample(range(n), q)))
        return sorted(F, key=sorted)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "a"):
        F = rand_family(8, 3, 10, [{0, 1, 2}, {3, 4, 5}], 1)
        run("random 3-sets, planted 2-packing", 8, 3, 2, F, rng)
        F = rand_family(9, 4, 8, [{0, 1, 2, 3}, {4, 5, 6, 7}], 2)
        run("random 4-sets, planted 2-packing", 9, 4, 2, F, rng)
        F = [frozenset({0, a, b}) for a, b in itertools.combinations(range(1, 7), 2)][:10]
        run("star family (every set contains 0): no 2-packing", 7, 3, 2, F, rng)
    if which in ("all", "b"):
        F = rand_family(9, 3, 9, [{0, 1, 2}, {3, 4, 5}, {6, 7, 8}], 3)
        print(f"   exact covers in this family: {count_packings(F, 3)}")
        run("exact cover by 3-sets, n=9", 9, 3, 3, F, rng, extra=(0,), sample=12)
    if which in ("all", "c"):
        r = random.Random(4)
        while True:
            F = rand_family(10, 3, 9, [{0, 1, 2}, {3, 4, 5}, {6, 7, 8}], r.randrange(10**6))
            if count_packings(F, 3) == 1: break
        rate_check(10, 3, 3, F, rng)
