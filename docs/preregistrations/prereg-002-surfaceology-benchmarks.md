# Preregistration 002 — Surfaceology / positivity-bootstrap benchmarks (B15-SURF)

**Status: DRAFT-FOR-FREEZE.** Becomes ACTIVE at the commit Erick pushes; the
commit hash of this file is written into every B15 certificate as
`prereg_ref = prereg-002@<commit>` (via `docs/preregistrations/prereg-002.freeze`)
together with `prereg_sha256` (content hash of this file). Rules may not change
between ACTIVE and evaluation (SPEC §6). Thresholds marked **[proposal]** are
carried as `asm:` taint until §9-OI-2 of the work order is answered.

Work order: `archive/sprint-drops/gv-b15-surf.zip` → `AGENT_WORK_ORDER_B15-SURF.md`.
Research report: `docs/notes/surfaceology-hidden-zeros-review.md`.

## Frozen machinery
- Verifier ops (no new op): `SCHUR_PIVOT_EXACT` = `b1_moment_solver.exact.psd_certificate`,
  `hankel`; `RANK_TEST` = `b3_electroweak.rank.matrix_rank`; linear bridge =
  `pir.symbolic.linear.solve / verify_solution / verify_farkas`; fingerprints =
  `pir.fingerprints.KnownGrammarDB / full_fingerprint / specific_fingerprint`
  (with `pole_set`, `zero_loci`, `split_structure` appended to `INVARIANT_KEYS`);
  candidate lattice = `pir.candidates.apply_rules / evaluate / lattice_fact`;
  certificates via `b1_moment_solver.certificate.save_certificate`;
  schema `schemas/b15_certificate.schema.json` (validated by `pir/jsonschema_mini.py`);
  standalone re-verifier `tools/verify_b15_certificate.py`.
- New code: `b15_surf/` (stdlib only; `fractions.Fraction` on every certified path;
  no `*_explore.py` files exist in this drop — no floats anywhere in B15).
- Tests: `tests/test_b15_pos.py`, `tests/test_b15_zero.py`, `tests/test_b15_gid.py`
  (auto-discovered by the unmodified `ci/run_all_certified.py` glob `tests/test_b*.py`).
- Seeds: dataset seed `20260923`; per-test seeds `20260923 + k`; `PIR_CI_SEED` recorded.

## §POS — B15-POS: verifier vs. a published EFT-hedron bound (R26)

**Source (transcribed at 2026-09-23, arXiv:2012.15849v2):**
- eq. (2.16)/(3.1): `M_IR(s,t) = {massless poles} + Σ_{k,q} a_{k,q} s^{k−q} t^q`,
  `s = (p1+p2)²`, `t = (p2+p3)²`, `u = −s−t`; `a_{k,q}` ↔ dimension-`2k+4` operators.
- eq. (7.17): forward limit `a_{k,0} = Σ_a p'_a x_a^k`, `p'_a > 0`, `x_a = 1/m_a²`, `k ≥ 2`.
- eq. (7.18)–(7.19): `(1, a_{3,0}/a_{2,0}, a_{4,0}/a_{2,0}, …)` lies in the convex hull of the
  half moment curve; its Hankel matrix `K[ã_0]` is totally positive.
- eq. (7.34)–(7.35): with a mass gap, `x_a ≤ 1` (units `M_Gap = 1`), so
  `a_{2,0} ≥ M_Gap² a_{3,0} ≥ … ≥ M_Gap^{2(k−2)} a_{k,0} ≥ 0`.

**Registered target bound (§9-OI-1 pending; carried as `asm:B15-OI-1-pending`):**
moments `μ0 = 1, μ1 = a_{3,0}/a_{2,0}, μ2 = a_{4,0}/a_{2,0}` (units `a_{2,0}=1`, `M_Gap=1`),
truncated at Hankel order `t = 1` (HEURISTIC truncation, located warning):
- lower wall (Hankel 2×2 minor, eq. 7.19): `μ0 μ2 − μ1² ≥ 0`;
- upper wall (gap chain, eq. 7.35, as 1×1 shifted-Hankel blocks): `μ0 − μ1 ≥ 0`, `μ1 − μ2 ≥ 0`.
At the registered slice `μ1 = 1/2` this is the two-sided statement **`1/4 ≤ μ2 ≤ 1/2`**;
the published wall values are `1/4` and `1/2`.

**Registered points / rays / thresholds**
- interior point `(μ1, μ2) = (1/2, 3/8)` → expected PERMITTED with pivot witness;
- exterior points `(1/2, 1/5)` (below the Hankel wall) and `(1/2, 3/5)` (above the gap wall)
  → expected REJECTED with (a) a negative Schur pivot and (b) an exact dual exclusion
  functional `y` (`y·μ ≥ 0` on the truncated moment cone, `y·μ_ext < 0`), nonnegativity
  certified by the identity `Σ y_j x^j = s0 (v0+v1 x)² + s1 x^p (1−x)`, `s0,s1 ≥ 0`, checked
  by coefficient matching through `verify_solution`; on the Hankel wall additionally the
  equality-form Farkas certificate on the single-atom system `[1, x*, x*²] w = μ`,
  `x* = μ1/μ0`, checked by `verify_farkas`;
- rays: `μ2` from `3/8` down to `0` (must bracket `1/4`) and up to `1` (must bracket `1/2`);
  bracket width threshold `1/2^20` **[proposal]**; reported as `certified_inner_interval`
  exactly as B1/B7;
- negative control: `μ2 → μ2 − 1/4` at the interior point → verdict must flip to REJECTED.

**Note on `verify_farkas` (recorded before the run):** the linear-equality Farkas lemma
cannot certify a conic (moment-cone) exclusion; the dual functional above is the exact
dual, and `verify_farkas` is used only where it is a genuine certificate (flat boundary).
See `docs/sprint-B15-SURF-reconciliation.md` §3.

## §ZERO — B15-ZERO: exact Tr(φ³) ground truth (R37)

**Conventions (transcribed, arXiv:2312.16282v3 §2):**
`X_{i,j} = (p_i+…+p_{j−1})²` (2.2); `c_{i,j} = −2 p_i·p_j` (2.3);
`c_{i,j} = X_{i,j} + X_{i+1,j+1} − X_{i,j+1} − X_{i+1,j}` (2.4); `X_{i,i+1} = 0`; `X_{i,j} = X_{j,i}`;
`n(n−3)/2` independent `X`'s. Route A: `A_n = Σ_T Π_{chords} 1/X` (e.g. 2405.09608 eqs. 2.3, 2.6);
Route B seeds `A_3 = 1`, `A_4 = 1/X_{1,3} + 1/X_{2,4}` (2401.05483 §I).

**Zero loci (transcribed):** §3.2 "Zeros" (causal diamond anchored at `X_B`); explicit index
ranges eq. (6.14) with `n_{a,b} = 0`: for the diamond anchored at `X_{1,i}`,
`c_{a,b} = 0` for `1 ≤ a ≤ i−2, i ≤ b ≤ n−1`; remaining loci = cyclic images; total
`n(n−3)/2` distinct loci (§6.1). Registered: all loci at `n ∈ {5, 6, 8}` (5, 9, 20 loci),
3 exact points per locus in the tests, 3 per locus per theory in the dataset.

**Negative control:** the same number of dense random rational linear constraints; the
amplitude must be nonzero at **≥ 95 %** of ≥ 60 points per `n` (registered).

**2-splits (transcribed, arXiv:2405.09608v1 §2):** eq. (1.7) `A_S = A_{S1}(x) A_{S2}(y)`;
sub-surfaces on a triangle `τ = (i,j,k)`: `P = (i,…,k)`, `Q = (i,j,k,k+1,…,i−1)`; blocks
`B = (i..j−1), C = (j..k−1), A = (k..i−1)` (2.9); map (2.10) `X_{a,b} → x_{k,b} + y_{a,i}`,
`X_{a,c} → x_{c,k} + y_{j,a}`, chords inside `P`/`Q` map trivially. Transcription guards:
eq. (2.7) [n=5, τ=(3,5,1)] and eq. (2.1) [n=6, τ=(5,1,3)] must be reproduced verbatim by
`b15_surf.zero.splits.split_point`; eq. (2.5) `{c_{1,3}=0, c_{3,5}=0}` must be the exact
vanishing set at n=6. Registered splits: n=5 `(3,5,1),(1,3,4)`; n=6 `(5,1,3),(1,3,5),(2,4,6)`;
n=8 `(1,3,5),(1,4,7),(2,5,7)`; ≥ 50 exact points each. Near-zero factorization
(2312.16282 eq. 3.3) at n=5 with `c_{1,4}=0`: `A_5 → (1/X_{1,3}+1/X_{2,5})(1/X_{1,4}+1/X_{3,5})`.

**δ-shift → NLSM (transcribed, arXiv:2401.05483v3 §I):** `X̃_{i,j} = X_{i,j} + δ_{i,j}`,
`δ_{e,e} = −δ_{o,o} = δ`, `δ_{o,e} = 0` (eq. 1, c-preserving);
**`lim_{δ→∞} δ^{2n−2} A^δ_{2n}(X̃) = A^NLSM_{2n}(X)`** (eq. 2), units `f_π = 1`;
4-point check `δ²(X13+X24)/((X13−δ)(X24+δ)) → −(X13+X24)`. The lower orders `O(δ^{−m})`
must cancel exactly (asserted, not assumed). Cross-check at 6 points: 2312.16282 eq. (7.20).
Independent route: `nlsm_feynman.py`, flavour-ordered Feynman rules derived from
`L_NLSM = (1/8λ²) Tr(∂U†∂U)`, `U = (I+λΦ)(I−λΦ)^{−1}` (2312.16282 eq. 7.4), `λ = 1`.
**Stipulation S1:** the metric signature (mostly-plus) is fixed once so that the 4-point
Feynman amplitude equals the letter's `−(X13+X24)`; the 6- and 8-point equalities then carry
no residual freedom. Registered test: exact equality at `n = 4, 6, 8` on 50/50/20 points.

**Dataset:** `data/b15/zero_dataset_v1.json`, classes `TrPhi3 (n=4,5,6,8)`,
`NLSM (n=4,6,8)`, `TrPhi3_shifted (n=4,6,8; registered hidden δ = 7/3)`; record kinds
`generic / locus / control / pole_probe / split_grid / split_ctrl`; probe values
`t ∈ {1/2,1/3,1/5,1/7,2/3,3/2,5/2,4/3,−1/2,−3/4}`; split grids `3×3`; canonical JSON,
`seed`, `generator_sha256`, `dataset_sha256`. **Two-class + deformation scope:** WP2.8
(YM-scaffolded) is NOT run (§9-OI-3 stretch); GID scope is `{TrPhi3, NLSM, TrPhi3_shifted}`.

## §GID — B15-GID: blind identification

**Input:** label-stripped dataset only (class ids `C0/C1/C2`, δ withheld); the pass never
imports `b15_surf.zero`. **Feature keys:** `pole_set` (exact single-variable model fit,
polynomial degree ≤ 3, pole location a linear unknown; overdetermined UNIQUE = certificate;
otherwise NONIDENTIFIABLE(model order)), `zero_loci` (exact vanishing), `split_structure`
(exact rank of the 3×3 grid; rank 1 = factorization), representation-specific
`pole_location_pattern` (δ-normalised); `redundancy_metric = C_{n−2}/(#atoms)` recorded only.
**Relabelling group:** dihedral `D_n` (cyclic + reflection); canonical form = orbit minimum;
registered test: a relabelled dataset gives byte-identical canonical features.
**δ recovery:** `X_{e,e}` poles at `−δ`, `X_{o,o}` poles at `+δ`; all must agree exactly →
FORCED δ; none present → `NONIDENTIFIABLE(delta-underdetermined)` (registered expectation
for NLSM); disagreement → `NONIDENTIFIABLE(delta-inconsistent)`.
**Menu / DB:** `data/b15/known_grammars_v1.json`; GVAR rules on predicates
`poles:{oe,ee,oo}`, `zeros:all_registered`, `splits:rank1`, `delta:{zero,nonzero,absent}`.
**Correlators (separate similarity & confidence, both `p/q`):**
`b15.gid.jaccard_atoms_v1` (primary: pole, zero, split and δ-class atoms; Jaccard similarity
to best; confidence = fraction of best-vs-runner-up discriminating atoms agreeing with best);
`b15.gid.invariants_only_v1` (δ-class removed); `b15.gid.zeros_only_v1` (poles and δ removed).
**Thresholds:** similarity ≥ 9/10 and confidence ≥ 9/10 **[proposal]**; no other tunable.
**Held-out:** fit `n ∈ {5, 6}`, apply to `n = 8`.
**Outcomes (primary correlator, n = 8):**
- **A (separable):** every class → single compatible family = DB identification, δ recovered
  on the shifted class, thresholds met → `FORCED` on the grammar label; claims-table
  `framework capability / BENCHMARK / E0`.
- **B (not separable):** ≥ 2 compatible families → `OBSERVATIONALLY_EQUIVALENT([…] mod δ-shift)`,
  `selected: null`; memo line on structural limits of fingerprint-only identification.
- **Tie** among ≥ 2 menu grammars with δ recovered → `AMBIGUOUS` (B8 semantics).
- **Anything else** → `NONIDENTIFIABLE(cause)` + open item.
The secondary correlators are evaluated and filed with the same vocabulary.

## Registered falsifiers
| ID | Trigger | Consequence |
|---|---|---|
| F1 | POS bracket excludes a published wall, or REJECTED without verified dual | verifier `REJECTED` for external bounds; `asm:B15-POS-unverified` on all B15 facts; program-level review line |
| F2 | Route A ≠ Route B anywhere | generator `REJECTED`; no dataset shipped |
| F3 | a registered zero does not vanish exactly | transcription erratum in edit-011, or open item if source-faithful; never a numerical patch |
| F4 | GID separable at n ≤ 6 but not at n = 8 | `fingerprint_v1` `REJECTED` (over-fit); not shipped |

## Stop rules
- Any F1–F4 trigger stops certification of the affected benchmark; the remaining WPs run
  with the registered taint. A schema or standalone-verifier failure is a test FAIL
  (Score = ∞), never a warning.
- No threshold, point, seed or index range above may be edited after freeze; a needed
  change is a new preregistration (prereg-003).

## Evaluation protocol
Run `python3 ci/run_all_certified.py` (unmodified); results to
`docs/sprint-B15-SURF-report.md` with PASS/FAIL per test, the GID outcome(s) obtained under
each correlator, redundancy metrics, errata, and the §9 open items with status. Both GID
outcomes are results. Claims-table rows are appended only after §9-OI-5 sign-off.
