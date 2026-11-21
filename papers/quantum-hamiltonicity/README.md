# Amplifying Determinant Sums

**A Quantum Algorithm for Undirected Hamiltonicity in O\*((8/3)^(n/2)) Time**

Konstantinos Dalkafoukis, with Claude Opus 5.5 (Anthropic) · **[Read the paper (PDF)](paper.pdf)**

> ⚠️ Draft, not peer reviewed.

## Main result

**Theorem.** There is a quantum algorithm that decides whether an undirected graph on n vertices has a Hamiltonian cycle in time

**O\*((8/3)^(n/2)) ⊂ O(1.6330ⁿ)**

It uses poly(n) qubits and no QRAM. It never answers "yes" on a graph without a Hamiltonian cycle, and it answers "yes" with probability at least 2/3 when one exists. The same holds for Hamiltonian path, and a cycle can be found with polynomial overhead.

**Corollary.** The undirected k-path problem can be decided in time O\*(1.6330ᵏ)·poly(n).

| Algorithm | Model | Graphs | Bound |
|---|---|---|---|
| Bellman; Held–Karp (1962) | classical | general | O\*(2ⁿ) |
| Björklund (2010) | classical, Monte Carlo | undirected | O\*(1.657ⁿ) |
| Björklund (2010) | classical, Monte Carlo | bipartite | O\*(1.415ⁿ) |
| Ambainis et al. (2019) | quantum, exponential QRAM | general (TSP) | O\*(1.728ⁿ) |
| Moylett–Linden–Montanaro (2017) | quantum | max degree 3 / 4 / 5–6 | O\*(1.110ⁿ) / O\*(1.301ⁿ) / O\*(1.680ⁿ) |
| **This paper** | **quantum, poly(n) qubits** | **undirected** | **O\*(1.6330ⁿ)** |

## How it works

1. **Björklund's test.** Pick a random half/half partition V₁ ∪ V₂ of the vertices. Call an arc of a Hamiltonian cycle *direct* if both its ends lie in V₁. The test detects cycles with exactly m direct arcs. It evaluates a sum of 2^(n/2+m) determinants over GF(2ᵏ) in polynomial space, and it never gives a false positive.
2. **Classical repetition.** A fixed cycle has m direct arcs with probability P ≈ C(n/2, m)² / 2ⁿ, so the test must be repeated about 1/P times. Balancing the cost 2^(n/2+m) against 1/P gives m ≈ 0.207n and Björklund's 1.657ⁿ.
3. **Quantum amplification.** The quantum algorithm runs the whole test reversibly on a superposition of all its random choices, then applies amplitude amplification. That needs only about 1/√P repetitions. The total cost 2^(n/2+m)/√P has exponent

   φ(μ) = ½ + μ + ½(1 − H(2μ)),  where μ = m/n and H is the binary entropy.

   It is minimised at μ = 1/6, where 2^φ = √(8/3) ≈ 1.63299.

Because rare successes are cheaper quantumly, the best cutoff moves from 0.207n to n/6, and each test gets cheaper.

## What's in the paper

- **§1 Introduction:** the result, the landscape of known bounds, and the relation to ETH (the result is consistent with it).
- **§2 Preliminaries:** the parts of Björklund's framework that are used (his Lemmas 9–11), and amplitude amplification.
- **§3 The algorithm:** two algorithm boxes, the register layout and the gate count.
- **§4 Analysis:** soundness, completeness, and the optimisation lemma with a closed-form optimum.
- **§5 Extensions:** Hamiltonian path, finding the cycle, k-path via the narrow sieve of Björklund–Husfeldt–Kaski–Koivisto, bounded-degree graphs, and weighted TSP in O\*(w·1.633ⁿ).
- **§6 Numerical checks:** exact finite-n convergence to √(8/3), and a Monte Carlo test of the counting formula.
- **§7 Limitations and open problems:** novelty, directed graphs, and other sieves (worked out for set packing in [quantum-set-packing](../quantum-set-packing/)).
- **Appendix:** an exact formula for the direct-arc distribution, and reproducibility notes.

## Status

The bound is **proved analytically** from published lemmas: Björklund's Lemmas 9–11 and Lemma 7 of the narrow-sieve paper. The scripts only cross-check it.

Still open:
- **Novelty.** Amplifying a Monte Carlo algorithm is a standard technique, so experts may already know this.
- **Independent review.**

## Reproduce

| What | Command | Needs |
|---|---|---|
| Build the PDF | `latexmk -pdf paper.tex` | TeX Live with TikZ/pgfplots |
| Classical vs quantum optimum, including unbalanced partitions | `python3 optimize.py` | NumPy, SciPy |
| Monte Carlo check and exact finite-n exponents | `python3 verify.py` | Python 3 |
| k-path exponents | `python3 kpath.py` | NumPy, SciPy |
| Figure data in `figs/` | `python3 make_data.py` | Python 3 |

`kpath.py`'s exact check at k = 30,000 takes several minutes. Set `\arxivtrue` at the top of `paper.tex` for the arXiv-compliant author line, which credits the model in a footnote instead of as an author.
