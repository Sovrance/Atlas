# Preregistration 004 — B15-POS at Hankel order t = 2 (EFT-hedron, gapped forward limit)

**Status: DRAFT-FOR-FREEZE.** Becomes ACTIVE at the commit Erick pushes after the pre-freeze
items below are answered. The commit hash is then written to `docs/preregistrations/prereg-004.freeze`
(`commit=` / `file=` lines, same format as `prereg-002.freeze`). Every B15-POS2 certificate carries
`prereg_ref = prereg-004@<commit>` and `prereg_sha256`. Rules may not change between ACTIVE and
evaluation (SPEC §6); a needed change is a new preregistration.

Origin: prereg-003 A3 (§9-OI-1). At t = 1 the registered bound `μ1² ≤ μ2 ≤ μ1` is nearly trivial,
so this preregistration is the stronger test of the same verifier against the same published
geometry. It is a new benchmark, not an upgrade of the B15-POS certificate, which stays as
certified under prereg-002/003.

**What is fixed before any code runs:** every wall value, point and expected verdict below was
derived in closed form (Schur complements and quadrature) without calling the repository's
verifier. The verifier has not been run on any t = 2 point.

## Source (transcribed 2026-09-24, arXiv:2012.15849v2)
- eq. (7.17)–(7.19): forward-limit moments `a_{k,0} = Σ_a p'_a x_a^k`, `p'_a > 0`. The Hankel
  matrix of `(1, a_{3,0}/a_{2,0}, a_{4,0}/a_{2,0}, …)` is totally positive.
- eq. (7.34)–(7.37): with a gap (units `M_Gap = 1`), `x_a ∈ [0, 1]`.
- eq. (7.45): example of a gapped (discrete-derivative) Hankel condition,
  `(a3 − a4)(a5 − a6) − (a4 − a5)² > 0`.
- eq. (7.46): with a known gap, the sequence, its first differences, its second differences, …
  each have a totally positive Hankel matrix ("the Hausdorff moment problem").

## Registered target (truncation t = 2)
Moments `μ_j = a_{j+2,0}/a_{2,0}`, `j = 0…4` (so `μ0 = 1`; moments through `a_{6,0}`), of a
positive measure on `[0, 1]`. The registered system is the necessary-and-sufficient condition
for the even truncated Hausdorff problem of degree 4:
- **Hankel wall:** `H = [μ_{i+j}]_{i,j=0..2} ⪰ 0` (3×3);
- **gap wall (localizing, weight x(1−x)):** `B = [μ_{i+j+1} − μ_{i+j+2}]_{i,j=0..1} ⪰ 0` (2×2).
  This is the index-shifted form of eq. (7.45) and is contained in the eq. (7.46) tower.

At fixed `(μ1, μ2, μ3)` the feasible `μ4` is an interval `[L, U]` with rational endpoints:
- `L = v^T A^{-1} v`, where `A = [[1, μ1], [μ1, μ2]]` and `v = (μ2, μ3)` (Schur complement of `H`);
- `U = μ3 − (μ2 − μ3)² / (μ1 − μ2)` (Schur complement of `B`).

The walls are single determinant conditions on a slice but not single inequalities in the
moments. This is the step up from t = 1.

## Registered slices, points and expected verdicts
| Slice | (μ1, μ2, μ3) | Ground-truth measure | L (Hankel wall) | U (gap wall) |
|---|---|---|---|---|
| S1 | (1/2, 1/3, 1/4) | Lebesgue on [0, 1] | **7/36** | **5/24** |
| S2 | (7/12, 7/16, 73/192) | atoms {1/4, 1/2, 1}, weights 1/3 each | **631/1792** | **641/1792** |

| Point | μ4 | Expected | Violated wall |
|---|---|---|---|
| S1-int | 1/5 | PERMITTED (pivots H: 1, 1/12, 1/180; B: 1/6, 1/120) | — |
| S1-lo | 13/72 | REJECTED (H pivot −1/72) | Hankel |
| S1-hi | 11/48 | REJECTED (B pivot −1/48) | gap |
| S2-int | 91/256 | PERMITTED (pivots H: 1, 7/72, 3/896; B: 7/48, 1/448) | — |
| S2-lo | 45/128 | REJECTED (H pivot −1/1792) | Hankel |
| S2-hi | 321/896 | REJECTED (B pivot −1/1792) | gap |

The pivots listed were computed independently with a plain rational LDLᵀ, so they serve as
registered expectations for `SCHUR_PIVOT_EXACT`.

**Independent route for S1 walls (quadrature ground truth).** `L = 7/36` is `μ4` of the
2-point Gauss–Legendre rule on [0, 1] (the lower principal representation; nodes are roots
of `x² − x + 1/6`). `U = 5/24` is `μ4` of Simpson's rule (atoms 0, 1/2, 1; weights 1/6, 2/3,
1/6; the upper principal representation). Both are checked exactly: the Gauss value by
reduction modulo `x² − x + 1/6`.

## Tests (registered)
- **T1 primal:** S1-int and S2-int → PERMITTED; witness = pivot chains of H and B.
- **T2 dual:** each exterior point → REJECTED with (a) the negative pivot and (b) an exact dual
  exclusion functional `y` (`y·μ_ext < 0`) whose nonnegativity on [0, 1] is certified in Gram
  form: `Σ_j y_j x^j = z₂ᵀ Q0 z₂ + x(1−x) z₁ᵀ Q1 z₁`, where `z₂ = (1, x, x²)`, `z₁ = (1, x)`, and
  `Q0, Q1` are rational PSD matrices. The identity is checked by coefficient matching through
  `pir.symbolic.linear.verify_solution`, and `Q0, Q1 ⪰ 0` by `SCHUR_PIVOT_EXACT`. This is
  **`DUAL_EXCLUSION_FUNCTIONAL` rev. 2**: the rev. 1 witness `s0 (v0+v1x)² + s1 x^p(1−x)` cannot
  express degree 4. No Farkas companion is registered at t = 2: the Hankel-wall boundary
  measure has irrational atoms (Gauss nodes), so the equality system is not over ℚ.
- **T3 brackets:** on each slice, exact bisection in μ4 from the interior point to each
  exterior point, to width `≤ 1/2^20`. Each `certified_inner_interval` must contain its
  registered wall (L or U).
- **T4 quadrature route:** the S1 walls equal the Gauss–Legendre and Simpson `μ4` exactly.
- **T5 t = 1 blindness control:** the prereg-002 t = 1 verifier (`b15_surf.pos`) returns
  PERMITTED on the `(μ1, μ2)` projection of every registered point. This is the demonstration
  that t = 2 adds content.
- **T6 negative control:** `μ4 → μ4 − 1/36` at S1-int (gives 31/180 < 7/36) flips the verdict to
  REJECTED.
- **T7 tower consistency:** at both interior points, every contiguous minor of the eq. (7.46)
  Hankel tower that is defined through μ4 (sequence, first and second differences) is
  positive. This is a consistency check of the registered two-matrix reduction against the
  paper's totally-positive statement.
- **T8 standalone re-verify + schema negatives:** as B15-POS T5/T6.

## Registered falsifiers
| ID | Trigger | Consequence |
|---|---|---|
| G1 | a bracket excludes its registered wall | verifier REJECTED for t = 2 external bounds; `asm:B15-POS2-unverified` on B15-POS2 facts |
| G2 | a REJECTED point without a verified Gram dual (identity or Q-pivots fail) | same as G1 |
| G3 | an interior point REJECTED, or an exterior point PERMITTED | same as G1 |
| G4 | T4 quadrature ≠ registered S1 wall | transcription/derivation erratum; stop before certification |
| G5 | T7 finds a negative tower minor at an interior point | the two-matrix reduction is wrong as registered; stop and open an item |

## Machinery and artifacts
- Frozen ops: `SCHUR_PIVOT_EXACT`, `DUAL_EXCLUSION_FUNCTIONAL` (rev. 2 witness form above),
  `pir.symbolic.linear.verify_solution`. `fractions.Fraction` on every certified path; no floats.
- Planned code: `b15_surf/pos/certify_t2.py`, `tests/test_b15_pos_t2.py`,
  `certificates/b15_pos_t2_certificate.json` (schema `schemas/b15_certificate.schema.json`,
  benchmark id `B15-POS2`), with the standalone tool `tools/verify_b15_certificate.py` extended
  for the Gram witness. No randomness; no seeds.
- Truncation at μ4 (`a_{6,0}`) is carried as a located HEURISTIC warning. The certificate
  itself is SOUND/E0 because every arithmetic step is exact.
- On success, claims-table row: *B15-POS2 verifier vs EFT-hedron gapped bound at t = 2 —
  recovered known result — BENCHMARK / E0*, appended only after sign-off.

## Pre-freeze items (to answer before this becomes ACTIVE)
- **P4-OI-1 — sufficiency citation.** The paper states necessity (eq. 7.46). That the two-matrix
  system is also *sufficient* for degree 4 on [0, 1] is the classical truncated Hausdorff
  result (Krein–Nudelman, *The Markov Moment Problem*, 1977; Curto–Fialkow). Confirm the exact
  theorem reference to cite, or accept T7 plus the explicit ground-truth measures as sufficient
  support for the registered walls.
- **P4-OI-2 — slices.** Confirm S1 (uniform) and S2 (three-atom) or name others. S2's band is
  deliberately narrow (width 10/1792).
- **P4-OI-3 — op revision.** Approve `DUAL_EXCLUSION_FUNCTIONAL` rev. 2 (Gram witness) for
  `docs/verifier-ops-v0.1.md`.
- **P4-OI-4 — width.** Keep `1/2^20` (as accepted for t = 1 in prereg-003 A1).
