# Exact checks and limits

`python3 code/run_checks.py` and `python3 -O code/run_checks.py` run:

1. Standard-library rational permutation-orbit identities and non-Hermitian/singular tensor recurrences for n=1 through 8, branching/matching scalar normalization through n=100.
2. SymPy full-matrix Casimir-projector checks through n=6, heterogeneous coefficient and matching checks through n=5, seven exact Clark interpolation cases including degree loss.
3. Reviewer endpoint checks: complex-scale interpolation, formal-degree loss, and a commuting but not doubly commuting 3-by-3 strict-contraction tuple.
4. An incorrect inverse-binomial weight is required to fail. Assertion ASTs are transformed into explicit exception branches and compiled with optimize=0, including when the entry point is invoked with -O.

Execution occurs in fresh temporary directories and does not overwrite archived results. These checks corroborate normalization and difficult finite cases, not the universal proof. SymPy is pinned to 1.14.0. Original historical evidence is preserved outside the public package; source and checker hashes bind the current replay.
