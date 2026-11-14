"""Data for the paper's figures (all exact or seeded, reproducible)."""
import random
from math import comb, log2, sqrt
H = lambda x: 0.0 if x <= 0 or x >= 1 else -x*log2(x) - (1-x)*log2(1-x)
phi = lambda mu: 0.5 + mu + 0.5*(1 - H(2*mu))   # quantum exponent
psi = lambda mu: 0.5 + mu + (1 - H(2*mu))       # classical exponent
with open("figs/curves.dat", "w") as f:
    f.write("mu quantum classical\n")
    for i in range(0, 301):
        mu = 0.05 + i*(0.30/300)
        f.write(f"{mu:.5f} {2**phi(mu):.6f} {2**psi(mu):.6f}\n")
def cnt(n, m):
    k = n//2; r = k - m
    return 0 if r <= 0 else comb(k-1, r-1)*comb(n-k-1, r-1)*n//r
with open("figs/convergence.dat", "w") as f:
    f.write("n quantum classical\n")
    for n in [60, 100, 200, 400, 1000, 2000, 4000, 8000, 12000]:
        lt = log2(comb(n, n//2)); lp = lambda m: log2(cnt(n, m)) - lt
        q = min(n/2 + m - 0.5*lp(m) for m in range(1, n//2))
        c = min(n/2 + m - lp(m) for m in range(1, n//2))
        f.write(f"{n} {2**(q/n):.6f} {2**(c/n):.6f}\n")
n = 24; rng = random.Random(1); T = 200000; hist = {}
for _ in range(T):
    V1 = set(rng.sample(range(n), n//2))
    m = sum(1 for i in range(n) if i in V1 and (i+1) % n in V1)
    hist[m] = hist.get(m, 0) + 1
tot = comb(n, n//2); maxdev = 0
with open("figs/montecarlo.dat", "w") as f:
    f.write("m mc exact\n")
    for m in range(0, n//2):
        mc = hist.get(m, 0)/T; ex = cnt(n, m)/tot; maxdev = max(maxdev, abs(mc-ex))
        f.write(f"{m} {mc:.6f} {ex:.6f}\n")
x = sqrt(2) - 1
print("max MC deviation", round(maxdev, 5), "| total exact prob", sum(cnt(n, m) for m in range(n//2))/tot)
print("quantum min", sqrt(8/3), "classical min", 2**(1.5 + x/2 - H(x)), "classical mu*", x/2)
