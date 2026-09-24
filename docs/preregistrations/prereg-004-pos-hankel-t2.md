# Preregistration 004 — B15-POS at Hankel order t = 2 (EFT-hedron, gapped forward limit)

**Status: ACTIVE (frozen).** Pre-freeze items answered 2026-09-24 (see the end of this file). The commit hash is then written to `docs/preregistrations/prereg-004.freeze`
(`commit=` / `file=` lines, same format as `prereg-002.freeze`). Every B15-POS2 certificate carries
`prereg_ref = prereg-004@<commit>` and `prereg_sha256`. Rules may not change between ACTIVE and
evaluation (SPEC §6); a needed change is a new preregistration.

Origin: prereg-003 A3 (§9-OI-1). At t = 1 the registered bound `μ1² ≤ μ2 ≤ μ1` is nearly trivial,
so this preregistration is the stronger test of the same verifier against the same published
geometry. It is a new benchmark, not an upgrade of the B15-POS certificate, which stays as
certified under prereg-002/003.

**What is fixed before any code runs:** every wall value, point and expected verdict below was
derived in closed form (Schur complements and quadrature) without calling the repository's
verifier, and merged to `main` (PR #19) in that state.

**Disclosure — checked once out-of-band before freeze.** After the merge and before freeze,
Erick re-derived every wall, pivot chain and quadrature value with exact arithmetic. That check
included running the pivot routine of `b15_surf.pos` on the six t = 2 points and the t = 1
verifier on their projections. All values agreed with the registered expectations, and T5's
premise held (t = 1 PERMITTED on all six). No registered value was changed as a result. The
registration was fixed before that run.

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
for the even truncated Hausdorff problem of degree 4 (warrant: see P4-OI-1 below — Kreĭn–Nudel′man
Thm III.2.3 as stated in Curto–Fialkow 1991, Remark 4.4):
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
- **T7 consistency with the paper's tower:** at both interior points, every contiguous minor of
  the eq. (7.46) Hankel tower that is defined through μ4 (sequence, first and second
  differences) is positive. This is a consistency check against the paper's totally-positive
  statement, not the warrant for the walls; the warrant is the theorem cited in P4-OI-1.
- **T8 standalone re-verify + schema negatives:** as B15-POS T5/T6.

## Registered falsifiers
| ID | Trigger | Consequence |
|---|---|---|
| G1 | a bracket excludes its registered wall | verifier REJECTED for t = 2 external bounds; `asm:B15-POS2-unverified` on B15-POS2 facts |
| G2 | a REJECTED point without a verified Gram dual (identity or Q-pivots fail) | same as G1 |
| G3 | an interior point REJECTED, or an exterior point PERMITTED | same as G1 |
| G4 | T4 quadrature ≠ registered S1 wall | transcription/derivation erratum; stop before certification |
| G5 | T7 finds a negative tower minor at an interior point | inconsistency between the cited theorem and the transcription of eq. (7.46); stop and open an item |

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

## Pre-freeze items — answered 2026-09-24 (Erick)
- **P4-OI-1 — sufficiency citation.** Cite a theorem; T7 is consistency, not warrant.
  The registered two-matrix system is the classical solvability criterion for the *even*
  truncated Hausdorff problem: for `γ_0…γ_{2k}` on `[a, b]`, a representing measure exists **iff**
  `A(k) ⪰ 0` and `(a+b)B(k−1) ⪰ ab·A(k−1) + C`, where `A(k) = (γ_{i+j})_{i,j=0..k}`,
  `B(k−1) = (γ_{i+j+1})_{i,j=0..k−1}`, and `C = (γ_{i+j+2})_{i,j=0..k−1}`. On `[0, 1]` with
  `k = 2` the second matrix is exactly `B = [μ_{i+j+1} − μ_{i+j+2}]`.
  **References (checked against the paper at freeze):** Kreĭn & Nudel′man, *The Markov Moment
  Problem and Extremal Problems*, AMS Transl. Math. Monographs 50 (1977), **Theorem III.2.3**
  (the explicit criterion, via Markov–Lukács), as stated in R. E. Curto & L. A. Fialkow,
  *Recursiveness, positivity, and truncated moment problems*, Houston J. Math. 17 (1991)
  603–635, **Remark 4.4** (p. 624). Curto–Fialkow's own even-case theorem is **Theorem 4.3**
  (range/extension form, equivalent). Their Theorem 4.1 is the *odd* case (m = 2k+1). The
  theorem number in Remark 4.4 comes from a text extraction of the scan, which reads
  "11.2.3". Remark 4.2 in the same extraction reads "III.2.4" for the odd case, so the
  intended reading is III.2.3. Recorded point: the even truncated Hausdorff case is one of the
  few truncated problems where *positivity alone* is sufficient: there is no rank or flatness
  condition, unlike the Hamburger case B1 handles. That is why the reduction to two matrices
  is exact.
- **P4-OI-2 — slices.** Keep S1 and S2 as registered, with no third slice. S1's walls are the
  two principal representations (lower = 2-point Gauss, upper = Simpson), which gives the
  quadrature route. S2's narrow band (10/1792 ≈ 2^−7.5) stresses the bracket and the
  exterior-point design. A slice with singular `A` (`μ2 = μ1²`), which exercises the
  zero-pivot path, is deferred to a future prereg-005.
- **P4-OI-3 — `DUAL_EXCLUSION_FUNCTIONAL` rev. 2: approved.** Rev. 2 is *complete* at this
  degree, not merely sufficient. By Markov–Lukács, every polynomial of degree `2m` that is
  nonnegative on [0, 1] equals `σ0 + x(1−x)σ1` with `deg σ0 ≤ 2m` and `deg σ1 ≤ 2m−2`; for
  m = 2 that is exactly a 3×3 Gram `Q0` and a 2×2 Gram `Q1`. Rev. 2 supersedes rev. 1 as the
  canonical witness form. Rev. 1 witnesses stay valid certificates and can be re-expressed in
  rev. 2 form (e.g. `1 − x = (1 − x)² + x(1 − x)`). Conditions written into the op: (i) `Q0` and
  `Q1` are rational; this is guaranteed for pivot-derived functionals, since a negative Schur
  pivot of `H` (or `B`) yields the rank-1 witness `Q = vvᵀ` with rational `v`. (ii)
  `verify_solution` checks the coefficient identity and `SCHUR_PIVOT_EXACT` checks
  `Q0, Q1 ⪰ 0`; both must pass, otherwise the rejection is uncertified (G2). There is no
  Farkas companion at t = 2 (the Gauss nodes are irrational); forcing one would be a
  forbidden repair.
- **P4-OI-4 — width.** Keep `1/2^20`. The narrowest band is ~2^−7.5, so the brackets sit about
  13 bits inside it. Every wall has a non-dyadic denominator (36, 24, 1792), so bisection
  never lands on a wall exactly, which makes the width threshold the operative guarantee.
