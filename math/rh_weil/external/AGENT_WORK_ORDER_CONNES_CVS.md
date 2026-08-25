---
status: HISTORICAL
superseded_by: ATLAS-RH-ENG-007
---

> **HISTORICAL / SUPERSEDED.** This is external Connes-CvS oracle: the work order that added the external cross-check layer. The layer exists; see external/PROVENANCE.md and external/CONNES_CVS_MAPPING.md for its live description.
> It is kept as the record of what was asked and when, not as instruction. For
> the live state read [`../AGENT_INSTRUCTIONS.md`](../AGENT_INSTRUCTIONS.md) and
> [`../certificates/work_order_status.json`](../certificates/work_order_status.json).
> No RH proof claim is made here or anywhere in this program.

# Agent Work Order — Connes-CvS Cross-Validation Layer

## Mission

Add an independent external oracle around `connes-cvs` without changing Atlas's
canonical RH/Weil formulas, certificate semantics, or claim boundary.

## Hard rules

1. `connes-cvs` is optional; core Atlas tests must run without it.
2. Pin the reviewed package/version in the cross-check environment.
3. Never compare CvS matrix entries to polynomial Gram entries by index.
4. Never promote an upstream floating result to Atlas E1.
5. Preserve the statement: finite matrices do not prove RH.
6. Every external artifact records package version and, when available, git SHA.

## WO-CVS-01 — install and provenance

- Install `connes-cvs==0.3.1` with `python-flint` in a dedicated environment.
- Record Python, mpmath, python-flint, FLINT/Arb, and package versions.
- Record upstream git commit if source checkout is used.
- Run upstream fast tests before Atlas cross-checks.

## WO-CVS-02 — low-level invariant cross-checks

Implement XC-01 through XC-03 from `CONNES_CVS_MAPPING.md`.

Required grid includes:

- `L = log(3), 1.1059498113, log(4)`;
- `t = 0, 1e-12/L, 1e-8/L, 0.1, 7, 21.218, 84`;
- all prime powers active on the selected cell.

Failure tolerance must scale with requested precision; do not hard-code a
binary64 tolerance for arbitrary-precision runs.

## WO-CVS-03 — scalar finite-Fourier cross-check

At `T=84`, reproduce Atlas scalar anchor values using Atlas's direct formula
while sourcing `h_plus` independently from connes-cvs.  This isolates the
archimedean multiplier implementation.

Output an E3 cross-check report.  If interval wrapping is added on both sides,
a later work order may promote the comparison itself to E1.

## WO-CVS-04 — projection bridge

Construct Fourier coefficients of Atlas basis functions in the CvS basis.
For each basis function `p`:

1. compute coefficients at increasing `N`;
2. measure projection residual in an explicit norm;
3. compare projected CvS quadratic forms with Atlas finite-Fourier values;
4. report convergence in `N` separately from arithmetic precision.

Do not claim equality until the normalization/basis map is symbolically proved.

## WO-CVS-05 — certificate-engineering cross-check

Use `connes_cvs.validation.arb_eigenpair_residual_bound` on controlled Atlas
finite matrices and compare with an Atlas-native residual verifier.

Acceptance requires matching scope statements: finite real-symmetric matrix
only, no infinite-operator conclusion.

## WO-CVS-06 — Atlas/PIR output

Emit a cross-validation certificate under

`math/rh_weil/certificates/external/`

with fields:

- upstream repository;
- upstream version/commit;
- dependency versions;
- Atlas source hash;
- test IDs XC-01..XC-06;
- evidence class per test;
- pass/fail/inconclusive;
- numeric precision and tolerances;
- explicit `rh_proof_claim: false`.

Do not modify existing E0/E1 RH certificates when an external cross-check fails.
Instead append a conflict/investigation artifact through PIR.

## Acceptance gate

The cross-validation layer is accepted when XC-01, XC-02 and XC-03 pass at two
precision levels and the adapter behaves cleanly when the optional dependency is
absent.  XC-04 and beyond are research extensions and should not block core CI.
