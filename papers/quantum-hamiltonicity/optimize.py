"""Exponent of Björklund's determinant-sum algorithm, classical vs amplitude-amplified.

A run: random partition V = V1 ∪ V2 with |V1| = a n. A fixed Hamiltonian cycle has
m = mu n "direct" arcs (consecutive V1-V1 pairs). The run detects it if m <= m_max,
at cost 2^{|V2| + m_max} (labelled cycle-cover sum over 2^{|L|} subsets, |L| = |V2| + m).
P(a cycle has exactly m direct arcs) = C(k,m) C(n-k, k-m) / C(n,k)   (k = a n, up to poly)
  [cyclic 0/1 strings with k ones forming k-m runs; matches Björklund's Lemma 10 at a=1/2]
Classical: cost * 1/P.   Quantum (amplitude amplification over the random choices): cost * 1/sqrt(P).
"""
import numpy as np
from scipy.optimize import minimize
H = lambda x: 0.0 if x <= 0 or x >= 1 else -x*np.log2(x) - (1-x)*np.log2(1-x)

def logP(a, mu):   # bits per n, <= 0
    if mu < 0 or mu > a or a - mu > 1 - a: return -np.inf
    return a*H(mu/a) + (1-a)*H((a-mu)/(1-a)) - H(a)

def expo(params, w):    # w = 1 classical, 1/2 quantum
    a, mu = params
    lp = logP(a, mu)
    return np.inf if not np.isfinite(lp) else (1 - a) + mu - w*lp

for label, w in [("classical", 1.0), ("quantum (amp. amplification)", 0.5)]:
    # balanced partition, as in Björklund
    mus = np.linspace(1e-4, 0.4999, 200001)
    vals = [expo((0.5, m), w) for m in mus]; i = int(np.argmin(vals))
    print(f"{label:30s} balanced a=1/2: mu*={mus[i]:.4f}  exponent={vals[i]:.5f}  base={2**vals[i]:.4f}  (P=2^-{-logP(0.5,mus[i]):.4f}n)")
    # free partition size
    best = None
    for a0 in np.linspace(0.3, 0.75, 46):
        for m0 in np.linspace(0.01, 0.5, 50):
            if not np.isfinite(expo((a0, m0), w)): continue
            r = minimize(lambda p: expo(p, w), [a0, m0], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=4000))
            if best is None or r.fun < best.fun: best = r
    a, mu = best.x
    print(f"{label:30s} free a:        a*={a:.4f} mu*={mu:.4f}  exponent={best.fun:.5f}  base={2**best.fun:.4f}  (P=2^-{-logP(a,mu):.4f}n)")
