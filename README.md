# Vibe researching with AI

Research done by directing AI models. I choose the problems, set the direction and review the results. The model does the derivations, proofs, numerical checks and writing. Each paper credits the model it was made with.

> ⚠️ Everything here is a draft and has not been peer reviewed.

## Papers

| Paper | Result | Best known before |
|---|---|---|
| [Amplifying Determinant Sums](papers/quantum-hamiltonicity/): quantum undirected Hamiltonicity | **O\*(1.633ⁿ)** = O\*((8/3)^(n/2)) | classical O\*(1.657ⁿ), quantum O\*(1.728ⁿ) |
| [Amplifying Narrow Sieves](papers/quantum-set-packing/): quantum 3-set packing | **O\*(2.804ᵖ)** = O\*((8 − 3√3)ᵖ) | classical O\*(3.343ᵖ) |

Both papers are by Konstantinos Dalkafoukis, with Claude Opus 5.5 (Anthropic). Each folder contains the PDF, the LaTeX source and the scripts that check the numbers.

## The common idea

The best classical algorithms for both problems repeat a cheap random test that rarely succeeds: Björklund's determinant sums for Hamiltonicity, and the narrow sieve of Björklund, Husfeldt, Kaski and Koivisto for set packing. A test that succeeds with probability P must be repeated about 1/P times. Amplitude amplification needs only about 1/√P repetitions. Once rare successes are cheaper, it pays to re-tune the test to be cheaper still. The re-tuned optimum has a closed form in both cases.

## Status

Both bounds are proved analytically from published lemmas, and the scripts only cross-check them. Two questions are still open: novelty (amplifying a Monte Carlo algorithm is standard, so experts may already know these results) and independent review.

## License

- **Papers** (LaTeX sources, PDFs, figures, data and text): [CC BY 4.0](LICENSE-papers). Please credit the authors when reusing them.
- **Code** (all `*.py` files): [MIT](LICENSE).
