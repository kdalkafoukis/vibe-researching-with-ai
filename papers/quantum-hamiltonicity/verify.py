"""Independent checks of the quantum Björklund exponent.
1) Monte Carlo: distribution of m (number of V1-V1 consecutive pairs on a fixed Hamiltonian cycle)
   under a uniformly random balanced partition, vs the exact formula C(n/2,m)^2/C(n,n/2) (cyclic, up to poly).
2) Exact large-n exponent: (1/n) log2 [ 2^{n/2+m} / sqrt(P(m)) ] minimised over m, using exact binomials,
   converging to 0.5*log2(8/3) = 0.70752.
"""
import random, math
from math import comb, log2
def exact_count(n, m):
    # cyclic binary strings of length n with k=n/2 ones and exactly m '11' adjacencies (cyclically): count via runs
    k = n // 2; r = k - m        # number of runs of ones = number of runs of zeros
    if r <= 0: return 0
    # place r runs of ones (compositions of k into r parts) and r runs of zeros (compositions of n-k into r parts), times n/r rotations
    return comb(k - 1, r - 1) * comb(n - k - 1, r - 1) * n // r
n = 24; rng = random.Random(1); trials = 200000; hist = {}
for _ in range(trials):
    V1 = set(rng.sample(range(n), n // 2))
    m = sum(1 for i in range(n) if i in V1 and (i + 1) % n in V1)
    hist[m] = hist.get(m, 0) + 1
print("n=24: m  MonteCarlo  exact")
tot = comb(n, n // 2)
for m in sorted(hist):
    print(f"  {m:2d}  {hist[m]/trials:.5f}  {exact_count(n, m)/tot:.5f}")
print("\nExact exponent vs n (quantum: 2^{n/2+m}/sqrt(P), classical: 2^{n/2+m}/P):")
for n in [60, 200, 1000, 4000, 20000]:
    lt = log2(comb(n, n // 2))   # logs of exact integers: the ratio underflows for large n
    best_q = min((n/2 + m - 0.5*(log2(exact_count(n, m)) - lt), m) for m in range(1, n//2))
    best_c = min((n/2 + m - (log2(exact_count(n, m)) - lt), m) for m in range(1, n//2))
    print(f"  n={n:6d}  quantum {best_q[0]/n:.5f} (m/n={best_q[1]/n:.4f})   classical {best_c[0]/n:.5f} (m/n={best_c[1]/n:.4f})")
print(f"limit: quantum 0.5*log2(8/3) = {0.5*log2(8/3):.5f} -> base {math.sqrt(8/3):.5f};  classical base {2**0.72845:.4f}")
