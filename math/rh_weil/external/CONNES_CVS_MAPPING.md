# Connes-CvS -> Atlas RH/Weil Technical Mapping v0.1


**No RH proof claim is made** by this mapping or by the oracle it describes.
Atlas's claim scope throughout is `finite_dimensional_weil_compression`.

## Decision

Use `connes-cvs` as an **optional external cross-validation oracle**, not as the
canonical Atlas proof engine.

Reviewed upstream package: `connes-cvs` v0.3.1 (MIT), repository
`akivag613/connes-cvs-`.

The two projects share the Weil explicit-formula ingredients but use different
finite-dimensional coordinate systems:

- **connes-cvs:** trigonometric Connes--van Suijlekom Galerkin basis
  `e_k(t)=exp(2*pi*i*k*t/L)`, matrix size `2N+1`.
- **Atlas RH/Weil:** polynomial, midpoint-parity, and endpoint-neutral bubble
  blocks such as `{1, q1, b}`.

Therefore **do not compare matrix entries by index**.

## File-level mapping

| Upstream | Atlas use | Decision |
|---|---|---|
| `connes_cvs/kernels.py::stable_A` | Independent stable `H0(t;L)` oracle | WRAP |
| `connes_cvs/kernels.py::stable_B` | Low-frequency linear-taper transform; useful cancellation reference | WRAP, do not confuse with quadratic bubble |
| `connes_cvs/operator.py::h_plus` | Archimedean multiplier cross-check | WRAP |
| `connes_cvs/operator.py::prime_powers_up_to` | von Mangoldt support/weight cross-check | WRAP |
| `connes_cvs/operator.py::build_galerkin_matrix` | Independent spectral/Galerkin diagnostic | WRAP with basis warning |
| `connes_cvs/operator.py::compute_ground_state` | Spectral diagnostic / regression oracle | OPTIONAL WRAP |
| `connes_cvs/validation.py::arb_eigenpair_residual_bound` | External finite-matrix residual certificate | WRAP as external certificate only |
| `runner.py`, `sweep.py` | Precision escalation, checkpoint, sweep engineering patterns | STUDY/REPLICATE patterns |
| zero extraction / Aitken extrapolation | Diagnostic research only | DO NOT USE for RH proof claims |

## What Atlas should independently replicate

These remain Atlas-native so the cross-check is genuinely independent:

1. polynomial overlap kernels `C_ij`, `K_ij`;
2. midpoint parity identities;
3. quadratic bubble basis `b=x(L-x)`;
4. compact real-space archimedean formula;
5. endpoint-jet tail bounds;
6. direct Fourier Taylor jets in support length `L`;
7. interval coverage in `L` and determinant/LDL certificates;
8. Atlas certificate/PIR evidence emission.

Do not copy these from upstream even if an equivalent implementation appears.

## Shared mathematical invariants to cross-check

### XC-01 — archimedean multiplier
For selected exact decimal `tau` values, compare

`Atlas h_plus(tau)` vs `connes_cvs.operator.h_plus(tau,dps)`.

Target: agreement far inside the requested working precision.  This is E3
cross-validation until Atlas wraps both values in interval enclosures.

### XC-02 — scalar transform
Compare

`H0(t;L) = integral_0^L exp(i t x) dx`

against `connes_cvs.kernels.stable_A(t,L)`, especially at `|tL| << 1`.

### XC-03 — prime-power ledger
Compare exact integer prime-power support and the weights

`Lambda(n)/sqrt(n)`

for every cell used by Atlas.  Compare `n` and prime base exactly; numeric logs
and weights may be compared at high precision.

### XC-04 — finite Fourier scalar quadratic form
At fixed `(L,T)`, compare the Atlas direct-Fourier scalar entry with an
independently assembled scalar from shared low-level ingredients.  Do not use
CvS matrix index `[N,N]` as a substitute unless equivalence is proved.

### XC-05 — basis-projection experiment
Approximate Atlas test functions `{1,q1,b}` in the CvS trigonometric basis,
construct a coefficient vector `v`, and compare `v^T Q v` with Atlas's finite
Fourier quadratic form.  Record projection/truncation error separately.

This is the first meaningful matrix-level comparison.

### XC-06 — finite eigenpair residual machinery
On Atlas-generated finite symmetric test matrices, feed an approximate
`(lambda,v)` into `arb_eigenpair_residual_bound` and compare with Atlas's own
finite-matrix validation path.  This tests certificate engineering, not RH.

## Evidence policy

- Upstream source identity/version/hash: provenance metadata.
- Upstream numerical value alone: E3 / HEURISTIC diagnostic.
- Upstream Arb residual on its exact supplied finite matrix: external E1, scope
  restricted to that finite matrix.
- Atlas theorem/certificate promotion: only after Atlas-native regeneration.
- Never infer infinite-dimensional positivity or RH from finite CvS spectra.

## Dependency strategy

Preferred: optional pinned dependency in a research/dev environment:

```bash
pip install 'connes-cvs==0.3.1'
```

For publication/reproducibility, also record the reviewed upstream git commit.
Do not make it a mandatory Atlas runtime dependency.

The upstream package is MIT licensed.  If code is ever vendored rather than
wrapped, preserve license/attribution and record provenance.
