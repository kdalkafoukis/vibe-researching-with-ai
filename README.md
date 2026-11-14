# Vibe researching with AI

Research done by directing AI models. I choose the problems, set the direction and review the results. The model does the derivations, proofs, numerical checks and writing. Each paper credits the model it was made with.

## Papers

### Quantum speedup of Björklund's algorithm

**Amplifying Determinant Sums: A Quantum Algorithm for Undirected Hamiltonicity in O\*((8/3)^{n/2}) Time**
Konstantinos Dalkafoukis, with Claude Opus 5.5 (Anthropic) · [PDF](papers/quantum-hamiltonicity/paper.pdf)

> ⚠️ Draft, not peer reviewed.

**Result.** This is a quantum algorithm that decides whether an undirected graph on n vertices has a Hamiltonian cycle. It runs in O\*(1.6330ⁿ) time, with poly(n) qubits and no QRAM, and it never gives a false "yes". The same argument gives O\*(1.6330ᵏ) for the k-path problem.

| Algorithm | Model | General undirected graphs |
|---|---|---|
| Björklund (2010) | classical | O\*(1.657ⁿ) |
| Ambainis et al. (2019) | quantum, exponential QRAM | O\*(1.728ⁿ) |
| **This paper** | **quantum, poly(n) qubits** | **O\*(1.633ⁿ)** |

To our knowledge, this is the first quantum algorithm faster than the best known classical one on general undirected graphs. On graphs of maximum degree 3 or 4, quantum backtracking (Moylett–Linden–Montanaro) was already faster than classical.

**Idea.** Björklund's algorithm repeats a random test. Each test costs 2^(n/2+m) and succeeds with an exponentially small probability P. Amplitude amplification needs 1/√P repetitions instead of 1/P. Re-optimising the cutoff then moves m from 0.207n to n/6, which gives (8/3)^(n/2).

**Status.** The bound is proved analytically from Björklund's published lemmas. The scripts only cross-check it.
- *Checked:* the closed-form optimum, Björklund's lemmas and his polynomial-space variant (read against the original paper), and the exact success-probability formula.
- *Open:* whether the result is novel, and independent review. Amplifying a Monte Carlo algorithm is a standard technique, so experts may already know it.

**Files** in `papers/quantum-hamiltonicity/`:
- `paper.tex`: build with `latexmk -pdf paper.tex`. Set `\arxivtrue` for the arXiv author line.
- `optimize.py`, `verify.py`, `kpath.py`: numerical checks (Python 3, NumPy, SciPy).
- `make_data.py`: generates the figure data in `figs/`.
