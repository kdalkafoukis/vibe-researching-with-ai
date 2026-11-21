#!/usr/bin/env python3
"""End-to-end check of Björklund's general-case test (arXiv:1008.0541, Lemmas 6, 9 and 11) on small graphs.
Python 3, standard library only; about 4 minutes.

For a graph G, a balanced partition V1|V2 and a cutoff m, compute the labelled cycle cover sum
Lambda(D, V2 ∪ L_m, f) at a random point of GF(2^16) three ways:
  p  : exponential-space route, Lemma 6/7: [r^|L|] sum_Y det(sum_{Z⊆Y} r^|Z| M_f(Z))
  q  : polynomial-space route, Section 4 / Lemma 11: walks B A^l B^T instead of paths
  fm : closed form of Lemma 9 I (sum over oriented Hamiltonian cycles of G with exactly m
       direct arcs, times a permanent over the m extra labels)
and compare with brute force: is there a Hamiltonian cycle with exactly m arcs inside V1?
Soundness requires Lambda == 0 exactly whenever there is none.
"""
import itertools, random, sys, time
from math import comb

# ---------------- GF(2^16) via log tables ----------------
Q = 65535
POLY = 0x1100B                                   # x^16+x^12+x^3+x+1
EXP = [0] * (2 * Q); LOG = [0] * (Q + 1)
v = 1
for i in range(Q):
    EXP[i] = v; LOG[v] = i
    v <<= 1
    if v & 0x10000: v ^= POLY
assert v == 1 and len(set(EXP[:Q])) == Q, "modulus not primitive"
for i in range(Q, 2 * Q): EXP[i] = EXP[i - Q]

def mul(a, b): return EXP[LOG[a] + LOG[b]] if a and b else 0
def inv(a): return EXP[Q - LOG[a]]

def det(M):                       # char 2: det == permanent, row swaps free
    n = len(M); M = [row[:] for row in M]; d = 1
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c]), None)
        if p is None: return 0
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]; d = mul(d, piv); ip = inv(piv); Rc = M[c]
        for i in range(c + 1, n):
            if M[i][c]:
                f = mul(M[i][c], ip); Ri = M[i]
                for j in range(c, n):
                    if Rc[j]: Ri[j] ^= mul(f, Rc[j])
    return d

def interp_coeff(xs, ys, t):      # coefficient of r^t of the interpolating polynomial
    N = [1]
    for a in xs:                  # N(r) = prod (r + a)
        R = [0] * (len(N) + 1)
        for i, c in enumerate(N): R[i + 1] ^= c; R[i] ^= mul(c, a)
        N = R
    res = 0
    for xi, yi in zip(xs, ys):
        deg = len(N) - 1; Qd = [0] * deg; Qd[deg - 1] = N[deg]
        for k in range(deg - 1, 0, -1): Qd[k - 1] = N[k] ^ mul(xi, Qd[k])
        den = 0; pw = 1
        for c in Qd: den ^= mul(c, pw); pw = mul(pw, xi)
        res ^= mul(yi, mul(Qd[t], inv(den)))
    return res

# ---------------- variables (s-oriented mirror) ----------------
class Vars:
    def __init__(self, rng, s): self.rng, self.s, self.val = rng, s, {}
    def _g(self, key):
        if key not in self.val: self.val[key] = self.rng.randrange(Q + 1)
        return self.val[key]
    def x(self, u, w):            # traverse edge u -> w
        return self._g(('x', u, w) if self.s in (u, w) else ('x', min(u, w), max(u, w)))
    def xd(self, u, v, d):        # direct arc u -> v with extra label d
        return self._g(('d', u, v, d) if self.s in (u, v) else ('d', min(u, v), max(u, v), d))

# ---------------- the three evaluations ----------------
def path_table(adj, V1, V2, var):
    idx = {w: i for i, w in enumerate(V2)}; memo = {}
    def h(w, v, mask):            # simple paths w -> ... -> v, interior (incl. w) == mask
        key = (w, v, mask)
        if key in memo: return memo[key]
        rest = mask & ~(1 << idx[w])
        if rest == 0: val = var.x(w, v) if v in adj[w] else 0
        else:
            val = 0
            for z in adj[w]:
                if z in idx and rest >> idx[z] & 1:
                    t = h(z, v, rest)
                    if t: val ^= mul(var.x(w, z), t)
        memo[key] = val; return val
    f = {}
    for u in V1:
        for v in V1:
            if u == v: continue
            for mask in range(1, 1 << len(V2)):
                val = 0
                for w in adj[u]:
                    if w in idx and mask >> idx[w] & 1:
                        t = h(w, v, mask)
                        if t: val ^= mul(var.x(u, w), t)
                f[(u, v, mask)] = val
    return f

def lam_p(adj, V1, V2, m, var):
    d, k = len(V1), len(V2); Ls = k + m
    f = path_table(adj, V1, V2, var)
    xs = [EXP[i] for i in range(1, d * Ls + 2)]; ys = []
    for r in xs:
        rp = [1]
        for _ in range(Ls): rp.append(mul(rp[-1], r))
        S = {}
        for u in V1:
            for v in V1:
                if u == v: continue
                arr = [0] * (1 << k)
                for mask in range(1, 1 << k): arr[mask] = mul(rp[bin(mask).count('1')], f[(u, v, mask)])
                for i in range(k):
                    for mask in range(1 << k):
                        if mask >> i & 1: arr[mask] ^= arr[mask ^ (1 << i)]
                S[(u, v)] = arr
        tot = 0
        for Y in range(1 << Ls):
            Yv, Yd = Y & ((1 << k) - 1), Y >> k
            M = [[0] * d for _ in range(d)]
            for a, u in enumerate(V1):
                for b, v in enumerate(V1):
                    if u == v: continue
                    e = S[(u, v)][Yv]
                    if v in adj[u]:
                        for j in range(m):
                            if Yd >> j & 1: e ^= mul(r, var.xd(u, v, j))
                    M[a][b] = e
            tot ^= det(M)
        ys.append(tot)
    return interp_coeff(xs, ys, Ls)

def lam_q(adj, V1, V2, m, var, n):
    d, k = len(V1), len(V2); Ls = k + m
    xs = [EXP[i] for i in range(1, d * n + 2)]; ys = []
    for r in xs:
        tot = 0
        for Y in range(1 << Ls):
            X1 = [V2[i] for i in range(k) if Y >> i & 1]; Yd = Y >> k
            M = [[0] * d for _ in range(d)]
            for a, u in enumerate(V1):
                vec = {w: mul(r, var.x(u, w)) for w in X1 if w in adj[u]}   # r^1 * Bout
                row = [0] * d
                for l in range(1, n + 1):
                    for b, v in enumerate(V1):
                        if v == u: continue
                        s_ = 0
                        for w, c in vec.items():
                            if v in adj[w]: s_ ^= mul(c, var.x(w, v))
                        row[b] ^= s_
                    if l == n: break
                    nv = {}
                    for w, c in vec.items():
                        for z in X1:
                            if z in adj[w]: nv[z] = nv.get(z, 0) ^ mul(c, mul(r, var.x(w, z)))
                    vec = nv
                for b, v in enumerate(V1):
                    if v != u and v in adj[u]:
                        for j in range(m):
                            if Yd >> j & 1: row[b] ^= mul(r, var.xd(u, v, j))
                M[a] = row
            tot ^= det(M)
        ys.append(tot)
    return interp_coeff(xs, ys, Ls)

def oriented_hcs(n, adj):
    out = []
    for p in itertools.permutations(range(1, n)):
        if p[0] > p[-1]: continue
        cyc = (0,) + p
        if all(cyc[(i + 1) % n] in adj[cyc[i]] for i in range(n)):
            out.append(cyc); out.append((0,) + tuple(reversed(p)))
    return out

def lam_formula(hcs, V1s, m, var):
    tot = 0
    for H in hcs:
        n = len(H); arcs = [(H[i], H[(i + 1) % n]) for i in range(n)]
        U = [(u, v) for u, v in arcs if u in V1s and v in V1s]
        if len(U) != m: continue
        base = 1
        for u, v in arcs:
            if not (u in V1s and v in V1s): base = mul(base, var.x(u, v))
        perm = 1 if m == 0 else det([[var.xd(u, v, j) for j in range(m)] for u, v in U])
        tot ^= mul(base, perm)
    return tot

def direct_counts(hcs, V1s):
    return {sum(1 for i in range(len(H)) if H[i] in V1s and H[(i + 1) % len(H)] in V1s) for H in hcs}

def run(name, n, edges, parts, use_q, seed):
    adj = {u: set() for u in range(n)}
    for u, v in edges: adj[u].add(v); adj[v].add(u)
    hcs = oriented_hcs(n, adj); rng = random.Random(seed)
    stats = dict(cases=0, sound_zero=0, nonzero_ok=0, sz_miss=0, mismatch=0, false_pos=0)
    hist = {}
    t0 = time.time()
    for V1 in parts:
        V1 = sorted(V1); V2 = [x for x in range(n) if x not in V1]; V1s = set(V1)
        have = direct_counts(hcs, V1s)
        for m in range(0, n // 2):
            var = Vars(rng, s=V1[0])
            lp = lam_p(adj, V1, V2, m, var); lf = lam_formula(hcs, V1s, m, var)
            lq = lam_q(adj, V1, V2, m, var, n) if use_q else lp
            stats['cases'] += 1
            if not (lp == lf == lq): stats['mismatch'] += 1
            if m in have:
                if lp: stats['nonzero_ok'] += 1; hist[m] = hist.get(m, 0) + 1
                else: stats['sz_miss'] += 1
            else:
                if lp: stats['false_pos'] += 1
                else: stats['sound_zero'] += 1
    print(f"{name}: n={n}, {len(hcs)//2} undirected Hamiltonian cycles, {len(parts)} partitions, "
          f"{'p+q+formula' if use_q else 'p+formula'}, {time.time()-t0:.0f}s")
    print("   " + ", ".join(f"{k}={v}" for k, v in stats.items()))
    return hist

if __name__ == "__main__":
    def allparts(n): return [set(c) for c in itertools.combinations(range(n), n // 2)]
    def randgraph(n, p, seed):
        r = random.Random(seed)
        return [(u, v) for u in range(n) for v in range(u + 1, n) if r.random() < p]
    for n, p, sd in [(6, 0.6, 1), (6, 1.0, 2), (8, 0.5, 3), (8, 0.7, 4)]:
        run(f"random graph G(n={n}, p={p})", n, randgraph(n, p, sd), allparts(n), use_q=True, seed=sd)
    hand6 = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3), (0, 4)]
    run("hand-made 6-vertex graph", 6, hand6, allparts(6), use_q=True, seed=9)
    hand8 = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (2, 6)]
    run("hand-made 8-vertex graph", 8, hand8, allparts(8), use_q=True, seed=10)
    # A single Hamiltonian cycle: the test must fire on exactly the partitions counted by the
    # paper's exact direct-arc formula (Appendix A).
    for n in (8, 10):
        cyc = [(i, (i + 1) % n) for i in range(n)]
        hist = run(f"cycle C_{n}", n, cyc, allparts(n), use_q=False, seed=n)
        k = n // 2; tot = comb(n, k)
        print("   m  fired on    exact formula (n/r)C(k-1,r-1)^2")
        for m in range(k):
            r_ = k - m
            print(f"   {m}  {hist.get(m, 0):4d}/{tot}    {n * comb(k - 1, r_ - 1) ** 2 // r_:4d}/{tot}")
