"""Checks for the amplified narrow sieve for p-packing of q-sets (Python 3, standard library only).

Model (Björklund-Husfeldt-Kaski-Koivisto, Section 4): a trial with parameters p0+p1+p2 = p costs
2^{q p0 + (q-1) p1 + (q-2) p2} and succeeds with probability P (their Lemma 12 / eq. (16)).
Classical cost ~ cost/P, quantum (amplitude amplification) ~ cost/sqrt(P).

1) closed form vs independent grid optimisation, and vs the classical table printed in BHKK
2) exhaustive check of the admissibility probability (BHKK Lemma 12) on small instances
3) exact finite-p optimum (perfect packing, n = pq) converging to the closed form
4) data for the figure
"""
import itertools, os
from math import comb, log2, sqrt, lgamma, log

H = lambda x: 0.0 if x <= 0 or x >= 1 else -x*log2(x) - (1-x)*log2(1-x)

def closed_form(q, w):
    """log2 of the optimal base per packed set, 2^q (Z^{q-1}/W^q)^w, x the positive root of
    C(q-1,2) x^2 - (c-1)(q-1) x - (c-1) = 0 with c = 2^{1/w}."""
    c = 2 ** (1 / w); a, b, d = comb(q - 1, 2), -(c - 1) * (q - 1), -(c - 1)
    x = (-b + sqrt(b * b - 4 * a * d)) / (2 * a)
    Z = 1 + q * x + comb(q, 2) * x * x
    W = 1 + (q - 1) * x + comb(q - 1, 2) * x * x
    beta = (1 / Z, q * x / Z, comb(q, 2) * x * x / Z)
    return q + w * ((q - 1) * log2(Z) - q * log2(W)), beta     # log2 of the base per packed set

def exponent(q, w, b1, b2):     # log2(base per packed set) for given fractions
    b0 = 1 - b1 - b2
    if min(b0, b1, b2) < 0: return None
    s = b1 + b2; t = b1 + 2 * b2
    lp = H(s) + (s * H(b2 / s) if s > 0 else 0) + b1 * log2(q) + b2 * log2(comb(q, 2)) - q * H(t / q)
    return q - t - w * lp

def grid(q, w, N=400):
    best = None
    for i in range(N + 1):
        for j in range(N + 1 - i):
            e = exponent(q, w, i / N, j / N)
            if e is not None and (best is None or e < best[0]): best = (e, i / N, j / N)
    e, b1, b2 = best; st = 1 / N
    for _ in range(60):
        for d1 in (-st, 0, st):
            for d2 in (-st, 0, st):
                v = exponent(q, w, b1 + d1, b2 + d2)
                if v is not None and v < e: e, b1, b2 = v, b1 + d1, b2 + d2
        st /= 2
    return 2 ** e

print("1) Closed form vs grid search. Base per packed set; per covered element in brackets.")
bhkk = {3: 1.4953, 4: 1.6413, 5: 1.7205, 6: 1.7707, 7: 1.8055, 8: 1.83118,
        10: 1.8663, 20: 1.9345, 50: 1.9741, 100: 1.9871, 500: 1.9975}   # BHKK table, c with f = c^{qp}
print("   q   classical          (BHKK table)  quantum            grid check (cl / qu)")
for q in [3, 4, 5, 6, 7, 8, 10, 20, 50, 100, 500]:
    cl, qu = closed_form(q, 1.0)[0], closed_form(q, 0.5)[0]
    g = f"{grid(q, 1.0):.4f} / {grid(q, 0.5):.4f}" if q <= 8 else "-"
    big = lambda e: f"{2**e:10.4f}" if q <= 10 else f"{'2^' + format(e, '.1f'):>10s}"
    print(f"  {q:3d}  {big(cl)} [{2**(cl/q):.5f}]  [{bhkk[q]}]  {big(qu)} [{2**(qu/q):.5f}]  {g}")
print(f"   q=3 exact: classical 9-4*sqrt(2) = {9-4*sqrt(2):.6f}, quantum 8-3*sqrt(3) = {8-3*sqrt(3):.6f}")
print(f"   q=3 quantum fractions (b0,b1,b2) = {tuple(round(b, 6) for b in closed_form(3, 0.5)[1])}")
print(f"   exact cover by 3-sets (p = n/3): classical {(9-4*sqrt(2))**(1/3):.5f}^n, quantum {(8-3*sqrt(3))**(1/3):.5f}^n")

print("\n2) Admissibility probability (BHKK Lemma 12), exhaustive over all (n1, n-n1)-partitions.")
def lemma12(n, q, p, n1, p1, p2):
    return comb(p, p1 + p2) * comb(p1 + p2, p2) * q**p1 * comb(q, 2)**p2 * comb(n - p*q, n1 - p1 - 2*p2) \
        if 0 <= n1 - p1 - 2*p2 <= n - p*q else 0
for n, q, p in [(10, 3, 3), (11, 4, 2), (12, 3, 4), (13, 5, 2)]:
    sets = [set(range(i*q, (i+1)*q)) for i in range(p)]; ok = True; checked = 0
    for n1 in range(n + 1):
        cnt = {}
        for U1 in itertools.combinations(range(n), n1):
            U1 = set(U1); k = [len(A & U1) for A in sets]
            if max(k) > 2: continue
            key = (k.count(1), k.count(2)); cnt[key] = cnt.get(key, 0) + 1
        for p1 in range(p + 1):
            for p2 in range(p + 1 - p1):
                checked += 1
                if cnt.get((p1, p2), 0) != lemma12(n, q, p, n1, p1, p2): ok = False
    print(f"   n={n}, q={q}, p={p}: {checked} (n1,p1,p2) cases, all counts match: {ok}")

print("\n3) Exact finite-p optimum for perfect packings (n = pq, n1 = p1+2p2), log-gamma binomials.")
def lc(a, b): return (lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)) / log(2)
conv = []
for p in [30, 100, 300, 1000, 3000]:
    q = 3; best = {}
    for w in (1.0, 0.5):
        e = min((q*(p-p1-p2) + (q-1)*p1 + (q-2)*p2
                 - w*(lc(p, p1+p2) + lc(p1+p2, p2) + p1*log2(q) + p2*log2(comb(q, 2)) - lc(p*q, p1+2*p2))) / p
                for p1 in range(p + 1) for p2 in range(p + 1 - p1))
        best[w] = 2 ** e
    conv.append((p, best[1.0], best[0.5]))
    print(f"   p={p:5d}: classical {best[1.0]:.4f}   quantum {best[0.5]:.4f}")
print(f"   limits: {9-4*sqrt(2):.4f} and {8-3*sqrt(3):.4f}")

os.makedirs("figs", exist_ok=True)
with open("figs/bases.dat", "w") as f:
    f.write("q classical quantum\n")
    for q in range(3, 31):
        f.write(f"{q} {2**(closed_form(q, 1.0)[0]/q):.6f} {2**(closed_form(q, 0.5)[0]/q):.6f}\n")
print("\nwrote figs/bases.dat")
