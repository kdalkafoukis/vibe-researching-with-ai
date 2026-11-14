"""Quantum amplification of the BHKK narrow-sieve k-path algorithm.
Per trial cost 2^{k1+l2}; success probability (BHKK Lemma 7):
  P = 2^{-k} C(k1+1, k-k1-l2) C(k-k1-1, l2),  k1 = g1 k, l2 = g2 k.
Classical exponent per k: g1+g2 - log2(P)/k ; quantum: g1+g2 - log2(P)/(2k)."""
import numpy as np
from math import comb, log2
from scipy.optimize import minimize
H = lambda x: 0.0 if x <= 0 or x >= 1 else -x*np.log2(x) - (1-x)*np.log2(1-x)
def logP(g1, g2):   # asymptotic log2(P)/k
    a, b = (1 - g1 - g2), g2
    if g1 <= 0 or a < 0 or a > g1 or b < 0 or b > 1 - g1: return -np.inf
    return -1 + g1*H(a/g1) + (1-g1)*H(b/(1-g1))
for name, w in [("classical", 1.0), ("quantum", 0.5)]:
    f = lambda p: (p[0] + p[1] - w*logP(*p)) if np.isfinite(logP(*p)) else 9
    best = min((minimize(f, [a, b], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=5000))
                for a in np.linspace(0.3, 0.7, 9) for b in np.linspace(0.05, 0.3, 6)), key=lambda r: r.fun)
    print(f"{name:9s} g1={best.x[0]:.5f} g2={best.x[1]:.5f} exponent={best.fun:.5f} base={2**best.fun:.5f}")
print("sqrt(8/3) =", np.sqrt(8/3))
# exact finite-k check at the quantum optimum
for k in [300, 3000, 30000]:
    best = None
    for k1 in range(int(0.4*k), int(0.6*k)+1, max(1, k//300)):
        for l2 in range(max(0, k-1-2*k1), k-k1):
            if not (k1 + l2 <= k-1 <= 2*k1 + l2): continue
            lp = -k + log2(comb(k1+1, k-k1-l2)) + log2(comb(k-k1-1, l2))
            v = (k1 + l2 - 0.5*lp) / k
            if best is None or v < best[0]: best = (v, k1/k, l2/k)
    print(f"k={k}: quantum exponent {best[0]:.5f} base {2**best[0]:.5f} at k1/k={best[1]:.3f} l2/k={best[2]:.4f}")
