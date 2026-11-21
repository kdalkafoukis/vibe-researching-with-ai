# Amplifying Narrow Sieves

**Quantum 3-Set Packing in O\*((8 − 3√3)ᵖ) Time**

Konstantinos Dalkafoukis, with Claude Opus 5.5 (Anthropic) · **[Read the paper (PDF)](paper.pdf)**

> ⚠️ Draft, not peer reviewed.

## Main result

**Theorem.** There is a quantum algorithm that decides whether a family of 3-element sets contains p pairwise disjoint members in time

**O\*((8 − 3√3)ᵖ) ⊂ O(2.8039ᵖ)**

It uses poly(n) qubits and no QRAM, and it never answers "yes" when no packing exists. The fastest known randomized classical algorithm is the narrow sieve of Björklund, Husfeldt, Kaski and Koivisto, which runs in O\*(3.3432ᵖ). The paper shows that this classical base is exactly 9 − 4√2.

**Also in the paper:**
- A closed form for the optimal base for sets of any size q, classical and quantum.
- Exact cover by 3-sets in O\*(1.4102ⁿ), compared with 1.4953ⁿ classically.

| q | Classical (narrow sieve) | Quantum (this paper) |
|---|---|---|
| 3 | 3.3431ᵖ | **2.8038ᵖ** |
| 4 | 7.2561ᵖ | **5.9752ᵖ** |
| 5 | 15.0717ᵖ | **12.3154ᵖ** |
| 6 | 30.8178ᵖ | **25.0732ᵖ** |

## How it works

This is the same idea as the [Hamiltonicity paper](../quantum-hamiltonicity/):

1. Each run of the narrow sieve picks a random partition of the universe and evaluates a sum of 2^|L| terms.
2. The run detects a fixed packing only if the partition cuts it in a prescribed pattern, which happens with an exponentially small probability P.
3. Classically about 1/P runs are needed. Amplitude amplification needs about 1/√P.
4. Re-optimising which pattern to look for gives the closed form. The optimum is a Gibbs distribution whose parameter solves a quadratic equation.

## Status

The bound is **proved analytically** from the narrow-sieve paper's Lemmas 12 and 13 (§4). The scripts only cross-check it. `sieve_check.py` runs the sieve itself end to end on small instances. Across 1,602 cases it matched the closed formula every time, never gave a false positive, and on an instance with a unique packing it fired on exactly the fraction of partitions that Lemma 12 predicts.

Still open:
- **Novelty.** Amplifying a Monte Carlo algorithm is a standard technique, so experts may already know this.
- **Independent review.** Tests on small instances are not a proof.

## Reproduce

| What | Command |
|---|---|
| Build the PDF | `latexmk -pdf paper.tex` |
| Analysis checks, plus the figure data in `figs/` | `python3 packing.py` (standard library only, about 15 s) |
| The sieve end to end on small instances | `python3 sieve_check.py` (standard library only, about 25 s) |

`packing.py` checks four things:
- the closed form against a grid search;
- the classical values against the narrow-sieve paper's table;
- Lemma 12 by exhaustive enumeration;
- convergence at finite p.

Set `\arxivtrue` at the top of `paper.tex` for the arXiv-compliant author line.
