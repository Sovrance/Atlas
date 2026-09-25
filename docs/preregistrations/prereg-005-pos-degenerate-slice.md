# Preregistration 005 — B15-POS3: degenerate (variance-zero) slices at Hankel order t = 2

**Status: ACTIVE (frozen).** Pre-freeze items answered 2026-09-25 (see the end of this file). The
freeze commit hash is written to `docs/preregistrations/prereg-005.freeze` (same format as
prereg-002/004). Every B15-POS3 certificate carries `prereg_ref = prereg-005@<commit>` and
`prereg_sha256`. Rules may not change between ACTIVE and evaluation (SPEC §6); a needed change is
a new preregistration.

Origin: prereg-004 P4-OI-2 (Erick, 2026-09-24): "a slice where A = [[1, μ1], [μ1, μ2]] is
singular (μ2 = μ1²) … the L-formula breaks and the code path for zero pivots gets exercised".

**What is fixed before any code runs:** every value below was derived in closed form and checked
with a standalone rational elimination script written for this draft. The repository's verifier
(`b1_moment_solver.exact.psd_certificate`, `b15_surf.pos.*`) has **not** been run on any point
registered here. The draft was merged to `main` in that state (PR #21).

**Disclosure — checked once out-of-band before freeze.** After the merge and before freeze,
Erick reproduced every pivot chain, all five dual directions (`vᵀMv = −1/64, −1/64, −1, −1/64,
−1/64`), the single-atom Farkas system (UNIQUE at D-on/E-on, INCONSISTENT with a verified Farkas
vector at all five exterior points) and the two U3 neighbours at D (`65535/2^20` fails H,
`65537/2^20` fails B). That check used the repository's pivot routine. Every value agreed with
the registration; none was changed. The values were fixed and merged before that run.

## Target
The system is the same as prereg-004: `H = [μ_{i+j}]_{i,j=0..2} ⪰ 0` and
`B = [μ_{i+j+1} − μ_{i+j+2}]_{i,j=0..1} ⪰ 0` (degree-4 truncated Hausdorff problem on [0, 1];
warrant [KrN, Theorem II.2.3] as cited in Curto–Fialkow 1991, Remark 4.4 — see prereg-004 errata E1;
the KrN chapter number is **VERIFY-BEFORE-USE** against the book itself).

**Why these slices are degenerate.** If `μ2 = μ1²`, the variance of any representing measure
is `μ2 − μ1² = 0`, so the measure is the Dirac mass `δ_{μ1}`. That forces `μ3 = μ1³` and
`μ4 = μ1⁴`, so the feasible set on the slice is a **single point**. This is the singular case in
which Curto–Fialkow 1991 Remark 4.5 notes the solution is unique. Consequences for the machinery:
- prereg-004's lower-wall formula `L = vᵀA⁻¹v` is undefined (`det A = 0`);
- on slice E (atom at the gap edge, `μ1 = 1`) the upper-wall formula `U` is also undefined
  (`μ1 − μ2 = 0`);
- the pivot chains contain **zero pivots**, both with vanishing rows (PSD) and, at D-off, with a
  non-vanishing row (indefinite via the zero-diagonal rule of `SCHUR_PIVOT_EXACT`).

## Registered slices and points (a = μ1; μ = (1, a, a², μ3, μ4))
| Point | (μ3, μ4) | Expected | H pivots | B pivots | Violated |
|---|---|---|---|---|---|
| D-on | (1/8, 1/16) | PERMITTED (PSD, not PD) | 1, 0, 0 | 1/4, 0 | — |
| D-lo | (1/8, 3/64) | REJECTED | 1, 0, **−1/64** | 1/4, 1/64 (PD) | Hankel (negative pivot after a zero pivot) |
| D-hi | (1/8, 5/64) | REJECTED | 1, 0, 1/64 (PSD) | 1/4, **−1/64** | gap |
| D-off | (9/64, 1/16) | REJECTED | 1, **0** (row not zero; stops) | 1/4, 31/1024 (PD) | Hankel (zero-pivot indefiniteness) |
| E-on | (1, 1) | PERMITTED (PSD, not PD) | 1, 0, 0 | **0, 0** | — |
| E-lo | (1, 63/64) | REJECTED | 1, 0, **−1/64** | 0, 1/64 (PSD) | Hankel |
| E-hi | (1, 65/64) | REJECTED | 1, 0, 1/64 (PSD) | 0, **−1/64** | gap (negative pivot after a zero pivot) |

Slice D: a = 1/2. Slice E: a = 1 (atom at the gap, where B = 0 at E-on). The feasible points are
exactly D-on and E-on (`μ4 = a⁴`). On D, `U = μ3 − (μ2 − μ3)²/(μ1 − μ2) = 1/16 = a⁴`, which is
consistent with the point.

## Registered dual construction (extends `DUAL_EXCLUSION_FUNCTIONAL` rev. 2 to zero pivots)
Let the symmetric elimination of `M` (`H` or `B`) stop at index k, with current Schur complement `S`.
- **Negative pivot** (`S_kk < 0`): `w = e_k`.
- **Zero pivot with a non-vanishing row** (`S_kk = 0`, first j > k with `S_kj ≠ 0`):
  `w = t e_k + e_j` with `t = −(S_jj + 1) / (2 S_kj)`, so `wᵀ S w = −1` exactly.
- **Lift** to `v` in the original coordinates. Set `v_i = 0` at earlier zero pivots (their rows
  vanish). Solve `M[I,I] v_I = −M[I, k:] w_{k:}` over the earlier nonzero-pivot indices `I`. Then
  `vᵀ M v = wᵀ S w`. The Gram witness is `Q = vvᵀ` (rational), checked exactly as in rev. 2
  (`verify_solution` for the identity, `SCHUR_PIVOT_EXACT` for Q ⪰ 0).
- **Generality of the normalisation.** For `w = t e_k + e_j`, `wᵀSw = 2t S_kj + S_jj` (since
  `S_kk = 0`), so `t = −(S_jj + 1)/(2 S_kj)` gives `−1` for every value of `S_jj`, including
  `S_jj = 0`. No further special case exists; this is the complete rule.
- **Recorded lift.** The dual records, per matrix, the stopping index `k`, the kind of stop
  (`negative_pivot` or `zero_pivot_nonzero_row`, with `j` and `t` for the latter), the earlier
  indices set to 0 as zero-pivot rows, and the indices solved from `M[I,I]`. A re-verifier must
  reproduce `v` from these, not only check `vᵀMv`.
- **Value is a certificate, not a margin.** Where the zero-pivot normalisation applies, `y·μ = −1`
  is fixed by construction. No report may read `y·μ` as a distance to the feasible set.

Registered directions and values (computed by the standalone script):
| Point | matrix | v | vᵀMv = y·μ |
|---|---|---|---|
| D-lo | H | (−1/4, 0, 1) | −1/64 |
| D-hi | B | (−1/2, 1) | −1/64 |
| D-off | H | (63/4, −32, 1) | −1 |
| E-lo | H | (−1, 0, 1) | −1/64 |
| E-hi | B | (0, 1) | −1/64 |

**Farkas certificate (complete on these slices).** Every representing measure is `δ_a`, so the
single-atom equality system `[1, a, a², a³, a⁴]ᵀ w = μ` is consistent (w = 1) at D-on and E-on and
inconsistent at every other point. Because `a` is rational, the system is over ℚ, and the Farkas
vector from `pir.symbolic.linear.solve`, checked by `verify_farkas`, is a genuine and here
*complete* certificate of infeasibility **on these slices**. It is not only a companion, unlike at
the prereg-002 flat boundary, and unlike prereg-004, where the atoms are irrational.
- **Recorded premise.** Completeness holds only because the slice imposes `μ2 = μ1²`: variance
  zero ⇒ the unique candidate measure is `δ_a`, so "not `δ_a`" is equivalent to "no measure on
  [0, 1]". Each certificate records this premise as a registered fact of the slice (not an
  `asm:` taint). Off these slices the Farkas system reverts to companion status (prereg-002) or
  is unavailable (prereg-004).
- **Two statements that coincide here.** The Gram dual certifies "no representing measure on
  [0, 1]"; the Farkas vector certifies "not the unique candidate `δ_a`". Both are required on
  every exterior point; neither substitutes for the other in the certificate.

## Tests (registered)
- **U1 primal:** D-on, E-on → PERMITTED with the registered chains (status PSD, not PD; zero
  pivots with vanishing rows accepted).
- **U2 dual:** the five exterior points → REJECTED. Stopping pivots equal the registered ones,
  and the violated wall is as registered. The dual follows the registered construction, with
  `v` and `vᵀMv` exactly as tabulated; the Gram witness is certified (identity + Q ⪰ 0).
- **U3 degenerate brackets:** on each slice, `μ4 = a⁴ ± 1/2^20` are both REJECTED and `μ4 = a⁴` is
  PERMITTED. Since there is no interval, the feasible set on the slice is the point `{a⁴}` at the
  registered resolution (`1/2^20`, as prereg-004). Expected at D: `65535/1048576` fails H and
  `65537/1048576` fails B. **Schema:** U3 is not a bracket, so the certificate must not populate
  `certified_inner_interval` (whose B1/B7 meaning is "inner feasible interval"). It uses a
  distinct field `certified_point`: the permitted value `a⁴` plus the two rejected neighbours,
  each with its violated wall.
- **U4 Farkas:** the single-atom system is consistent at D-on and E-on (w = 1), and at all five
  exterior points it is inconsistent with `verify_farkas` true.
- **U5 formula controls:** `det A = 0` on both slices, so `L` is undefined on D and on E; `U` is
  defined and equals `a⁴` on D; `μ1 − μ2 = 0` on E, so `U` is **also** undefined on E (asserted
  explicitly: on E neither formula exists). The implementation must not call either formula on
  these slices.
- **U6 no silent path:** reading prereg-004's `certify_t2.negative_direction` shows it refuses
  these points. It asserts on D-lo, E-lo and E-hi (singular or zero leading block) and raises `ValueError`
  on D-off.
  This is the expected behaviour **to be confirmed at evaluation**; it was not run for this draft.
  The prereg-005 implementation must certify all seven points through the registered
  construction. Neither version may emit a verdict through an unregistered path.
- **U7 standalone re-verify + schema negatives:** as B15-POS2 T8.

## Registered falsifiers
| ID | Trigger | Consequence |
|---|---|---|
| K1 | D-on or E-on REJECTED, or an exterior point PERMITTED | zero-pivot handling REJECTED; `asm:B15-POS3-unverified` |
| K2 | an exterior point without a certified Gram dual, or with `v` / `vᵀMv` ≠ registered | same as K1 |
| K3 | a `±1/2^20` neighbour of `a⁴` PERMITTED | the feasible set is not a point; the reduction is wrong as registered; stop |
| K4 | the Farkas system inconsistent at D-on/E-on or consistent at an exterior point | derivation erratum; stop |
| K5 | any verdict reached by calling either wall formula on these slices — L (undefined on D and E) or U (undefined on E) | implementation REJECTED |

## Machinery and artifacts (planned)
- Extend `b15_surf/pos/certify_t2.py` with the registered zero-pivot direction. The prereg-004 code
  path stays unchanged for B15-POS2 reproducibility; the new function is used by B15-POS3.
- `tests/test_b15_pos3.py`, `certificates/b15_pos3_certificate.json` (benchmark `B15-POS3`,
  prereg binding `prereg-005@<commit>`), standalone-tool support for zero pivots and the Farkas
  check at t = 2.
- `DUAL_EXCLUSION_FUNCTIONAL`: document the zero-pivot construction as a rev. 2 clarification. The
  witness form is unchanged; only the derivation of the rank-1 `Q` is extended.
- On success, the claims-table row *B15-POS3 zero-pivot / degenerate-slice verification —
  framework capability — BENCHMARK / E0* is appended after sign-off.

## Pre-freeze items — answered 2026-09-25 (Erick)
- **P5-OI-1 — slices.** Confirm D and E; no a = 0 slice. D exercises "zero pivot with a vanishing
  row, then a decisive pivot" in H, and D-off is the only point that reaches the
  zero-diagonal/non-zero-row indefiniteness rule of `SCHUR_PIVOT_EXACT`, a path never certified
  before. E adds the case prereg-004 could not reach: `B ≡ 0` at E-on (both gap pivots zero), so
  the gap wall is itself degenerate. An a = 0 slice would give H the chain of D and B the chain of
  E, so no new path. Written in before freeze: on E the upper formula `U` is undefined too (U5
  asserts it explicitly; K5 fires on either formula).
- **P5-OI-2 — zero-pivot direction: approved as a rev. 2 clarification, not rev. 3.** The witness
  form (rank-1 Gram `Q = vvᵀ`, identity by `verify_solution`, `Q ⪰ 0` by pivots) is unchanged;
  only the derivation of `v` is extended. The normalisation `wᵀSw = −1` is general (see
  "Generality of the normalisation"). Added to the op text: (i) the lift is recorded so a
  re-verifier reproduces `v`; (ii) `y·μ` on these points is a certificate, not a margin.
- **P5-OI-3 — Farkas role.** Complete on the registered slices, with the Gram dual still
  mandatory. The slice condition is a recorded premise of the completeness claim (registered
  fact, not `asm:`). The two certificates certify different statements that coincide here (see
  "Two statements that coincide here").
- **P5-OI-4 — resolution.** Keep `1/2^20`. U3 populates `certified_point`, not
  `certified_inner_interval`.
- **KrN numbering.** Unresolved from here; keep "Theorem II.2.3 as cited in Curto–Fialkow 1991,
  Remark 4.4", marked VERIFY-BEFORE-USE (prereg-004 errata E1).
