# AI index

## Identity and status
Anonymous; version 1.1.0-candidate; unrefereed written proof. [STATUS.md](STATUS.md), [PROVENANCE.md](PROVENANCE.md), [LICENSES.md](LICENSES.md).

## Exact claim and proof route
[paper.tex](paper.tex), Theorem thm:main: all finite n, nonzero symmetric multiaffine p stable on the open polydisk; actual-multidegree quotient is Schur–Agler, no monomial. Boundary zeros and common factors allowed. Read lem:twirl, prop:equal, lem:clark, eq:Agler and eq:operators in order. The spin-zero and n=1 cases are explicit; the full-union coefficient equations are not inferred from the solved recurrence alone. Corollary cor:matching uses unnormalised vectors; the subsequent D_l congruence is required for eigenvalues.

## Machine checking and trust boundary
[code/run_checks.py](code/run_checks.py), [requirements.txt](requirements.txt), [VERIFICATION.md](VERIFICATION.md). Run normally and under -O. Finite computations supplement rather than prove the unrestricted theorem. The negative control corrupts the weight denominator. [ASSURANCE.md](ASSURANCE.md) states what is not independently or formally verified.

## Review and predecessors
[REVIEW_RESPONSE.md](REVIEW_RESPONSE.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [prior-art.md](prior-art.md). Standard Schur–Weyl branching, scalar polarisation, Clark interpolation and radial approximation are inherited. The proposed contribution is automatic positivity and full coefficient transfer in this class. No absolute priority guarantee.

## Entry points
Human-readable [paper.pdf](paper.pdf); canonical editable [paper.tex](paper.tex); current normal/optimized replay logs and manifest supplied with final archive. Read claims before running code; successful code execution is not a proof audit.
