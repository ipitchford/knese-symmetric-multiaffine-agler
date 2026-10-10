# Stable symmetric multiaffine polynomials are Agler denominators

Anonymous · 1.1.0-candidate · Unrefereed written-proof candidate

For every finite number of variables, every nonzero symmetric multiaffine polynomial with no zeros in the open unit polydisk is an Agler denominator, using its actual multidegree and no extra monomial. Boundary zeros, repeated roots, common factors and complex coefficients are included.

Read [paper.pdf](paper.pdf), [paper.tex](paper.tex), [AI_INDEX.md](AI_INDEX.md), [ASSURANCE.md](ASSURANCE.md) and [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md). The unrestricted result rests on the written proof, not the sampled dimensions in the checkers.

## Replay

Python 3.12 or 3.13 and SymPy 1.14.0:

```sh
python3 -m pip install sympy==1.14.0
python3 code/run_checks.py
python3 -O code/run_checks.py
```

The entry point runs three exact checkers in temporary directories, preserves assertion semantics even under optimization, and rejects a deliberately corrupted inverse-binomial weight. Build the paper with two pdfLaTeX passes. No empirical dataset or proof-assistant formalisation is claimed.
