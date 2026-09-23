# AGENT WORK ORDER — B15-SURF (Surfaceology / Positivity-Bootstrap Integration)

**Repository:** github.com/Sovrance/Atlas (main)
**Branch:** `feature/b15-surf` — commit in coherent increments with prefix `B15-SURF:`; **do not merge, do not push.** Erick integrates.
**Status:** DRAFT-FOR-EXECUTION — §9 open items need Erick's answer before WP1.2 (target bound) and before any `docs/claims-table.md` row is added.
**Supersedes:** `B15-SURF_engineering_spec.md` (2026-09-23 draft written against the wrong repo layout; retired).
**Reconciled against:** `Atlas-main` tree as uploaded 2026-09-23 (488 files, 37 test files, certificates v0.3, SPEC v1.0, ledger through R35, atlas v0.6).
**Author:** Rosy 🌹 · **Owner:** Erick L. Gonzalez

---

## Mission

Turn the surfaceology / hidden-zeros deep-research report (2026-09-21) into one certified, pre-registered benchmark sprint that:

1. checks the repo's exact-PSD verifier against an **externally published** positivity bound (EFT-hedron, R26) with primal **and** dual certificates — the first time the verifier is scored against a bound the program did not itself construct;
2. builds a **ground-truth amplitude dataset** (Tr φ³ tree amplitudes, hidden zeros, 2-splits, δ-shift to NLSM) in exact `Fraction` arithmetic;
3. runs the existing blind-identification machinery (B8 lineage → `pir/fingerprints.KnownGrammarDB`, M1 `canon`, B12-RGRC rival framing) on that dataset, pre-registered for **both** a separable and a non-separable outcome;
4. produces one design note on constraint-manifest representations (P2 of the report) — no code.

**Without** changing any existing B1–B13 verdict, certificate, or atlas cell, and without expanding `ci/run_all_certified.py` (`.cursor/rules` prohibition; the new `tests/test_b15*.py` are auto-discovered by its existing glob).

## Required reading (in order)

1. `docs/SPECIFICATION.md` — §3 evidence levels, **§4 outcome vocabulary (LAW)**, §5 equivalence, §6 promotion/claims.
2. `docs/adr/ADR-0002-conditionality-as-taint.md` — conditionality is `asm:X` taint, **not** a `CONDITIONAL(X)` verdict.
3. `docs/verifier-ops-v0.1.md` — frozen op set; B15 adds **no new op**, it composes `SCHUR_PIVOT_EXACT` + `RANK_TEST` + the `pir/symbolic/linear.py` Farkas bridge.
4. `README.md` §§B1, B8, B12-RGRC, M1, M2, PIR (Stages 1–3).
5. `docs/pir-specification-v0.1.md` — fact schema (`witness`, `impossibility_certificate` fields already reserved).
6. `docs/preregistrations/prereg-001-gwtc-next-area.md` — prereg house style.
7. `docs/atlas-edits/edit-009-alpha-tension-memo.md` — memo house style; **note edit-010 is reserved** by edit-009's decision rule.
8. `docs/notes/related-programs-review.md` + `docs/references.md` R25–R27 — the bootstrap import already on file; B15 executes what R26 asked for ("add dual/exclusion functionals to the verifier layer").
9. `.cursor/rules/*` — work-zip protocol.
10. The research report *Surfaceology and Hidden Zeros: Discovery Heuristics for Constraint-Certified Physics* (attached as `docs/notes/surfaceology-hidden-zeros-review.md` in this drop).

## Pre-flight checks (do first; record answers in `docs/sprint-B15-SURF-reconciliation.md`)

- [ ] **Sprint ID collision.** Program records mention a merged `B14-KMAP`; this tree has no `b14_*`. Run `ls | grep -i b14; grep -rn "B14\|B15" docs/ README.md`. If B15 is taken on main, rename to the next free B\* and update every path below.
- [ ] **Vocabulary check.** `grep -rn ENTAILMENT_ .` returns nothing in this tree. If main has extended SPEC §4 with `ENTAILMENT_*` verdicts, they may be used where §3.3 says so; otherwise use the §4 set only.
- [ ] **Edit number.** `ls docs/atlas-edits/` → next free number after 009, skipping 010 (reserved). Expect `edit-011`.
- [ ] **Ledger number.** `tail docs/references.md` → next free R-number. Expect R36.
- [ ] `python3 ci/run_all_certified.py` is green on the untouched tree. Record the run-manifest hash. This is the baseline; **nothing in B15 may degrade it.**

## Create

```
b15_surf/
  __init__.py
  exact.py                 Fraction helpers (copy the b7_onsager/exact.py pattern: fmt, q)
  pos/
    __init__.py
    dispersive.py          moment/Hankel construction for the registered EFT-hedron bound (transcribed, cited)
    certify.py             PERMITTED / REJECTED with primal pivots + dual Farkas functional; exact bisection bracket
  zero/
    __init__.py
    kinematics.py          planar X_{i,j}, mesh c_{i,j}, rational point sampler (seeded)
    triangulations.py      Route A: Catalan enumeration, A_n = Σ_T Π 1/X
    recursion.py           Route B: factorization recursion from A_3=1, A_4
    zeros.py               registered zero loci (linear solve, exact) + random-constraint negative control
    splits.py              2-split near zeros (arXiv:2405.09608, transcribed)
    delta_shift.py         X_ee → X_ee+δ, X_oo → X_oo−δ; NLSM limit per arXiv:2401.05483 (transcribed)
    nlsm_feynman.py        independent NLSM flavor-ordered amplitudes (4-, 6-pt) from the paper's Feynman rules
    dataset.py             writes data/b15/zero_dataset_v1.json (exact "p/q" strings, seed, generator sha256)
  gid/
    __init__.py
    fingerprint.py         pole set / zero-locus set / split structure → feature dict for pir.fingerprints
    identify.py            blind identification → candidate forest (pir.candidates lattice) + verdict
tests/
  test_b15_pos.py          T1 primal PERMITTED · T2 dual REJECTED · T3 bracket · T4 negative control · T5 standalone re-verify
  test_b15_zero.py         T1 Catalan counts · T2 Route A=B (≥200 pts/n) · T3 zeros vanish · T4 negative control · T5 2-split · T6 δ-shift=NLSM · T7 dataset hash
  test_b15_gid.py          T1 fingerprint build · T2 M1 quotient · T3 DB identify · T4 δ recovery · T5 held-out n=8 · T6 outcome A/B/NONIDENTIFIABLE adjudication
certificates/
  b15_pos_certificate.json
  b15_zero_certificate.json
  b15_gid_certificate.json
schemas/
  b15_certificate.schema.json   validated with pir/jsonschema_mini.py (b13_cdl pattern)
data/b15/
  zero_dataset_v1.json
  known_grammars_v1.json
tools/
  verify_b15_certificate.py     standalone re-verifier: recomputes pivots + Farkas evaluation from certificate contents only
docs/
  preregistrations/prereg-002-surfaceology-benchmarks.md
  notes/surfaceology-hidden-zeros-review.md            (research report, verbatim, from this drop)
  notes/constraint-manifest-representation-v0.md       (WP4 design note, ≤3 pages)
  atlas-edits/edit-011-b15-surf-memo.md                (memo, no cell change)
  sprint-B15-SURF-reconciliation.md
  sprint-B15-SURF-report.md
```

**Modify (additively only):**
- `docs/references.md` — append block R36–R40 (see §12).
- `docs/claims-table.md` — append three rows **after Erick's sign-off** (§9-OI-5): B15-POS `recovered known result / BENCHMARK / E0`; B15-ZERO `recovered known result / BENCHMARK / E0`; B15-GID `framework capability / BENCHMARK / E0–E3` (E-level fixed by which outcome obtains).
- `README.md` — one B15 section in the existing per-benchmark style (run command, N/N, certificate path).
- `archive/sprint-drops/README.md` — one row for `gv-b15-surf.zip`.
- `pir/fingerprints.py` — **only** if needed: add `"pole_set", "zero_loci", "split_structure"` to `INVARIANT_KEYS`. `_project` skips absent keys, so existing full-fingerprint hashes are unchanged; add a regression test asserting that (`tests/test_pir_candidates.py` or a new `tests/test_b15_gid.py::T0`).

**Do not touch:** `ci/run_all_certified.py`, `docs/constant-atlas-v0.6.md`, any `b1_*`…`b13_*` module, any existing certificate, `math/rh_weil/`.

## Hard constraints

1. **Stdlib only; `fractions.Fraction` on every certified path.** Floats are confined to files named `*_explore.py`, which never write to `certificates/` or `data/`.
2. **Verdict vocabulary = SPEC §4 exactly:** `FORCED | PERMITTED | REJECTED | NONIDENTIFIABLE(cause) | OBSERVATIONALLY_EQUIVALENT(class) | APPARATUS_LIMITED | REPRESENTATION_DEPENDENT | AMBIGUOUS`. Conditionality is emitted as assumption-taint `asm:<id>` on the PIR fact (ADR-0002), never as a verdict string. Non-SPEC verdict strings must fail schema validation (negative test required).
3. **Certificate format = existing v0.3 shape** (`certificate_version`, `certificate_class`, `problem`, `timestamp_utc`, `headline`, `m_layer_stipulations`, `calibration_route`, `results`) **plus** the §5 fields. Written via `b1_moment_solver.certificate.save_certificate`. Content-addressed `certificate_id` as in `b13_cdl/docs/pir-bridge-v0.1.md`.
4. **Every pass tagged `SOUND | HEURISTIC`**; HEURISTIC passes carry located `warnings[]`; a HEURISTIC pass may not assert E0 (existing `pir/passes.py` negative test applies).
5. **`Score = ∞` on any uncertified hard constraint** — a certificate that fails `schemas/b15_certificate.schema.json` or `tools/verify_b15_certificate.py` is a test FAIL.
6. **No atlas cell changes; no new H or P.** B15 files a memo (`edit-011`) and claims-table rows only. B15 facts enter PIR at layer `DOMAIN`, namespace `domain:` (POS, ZERO) / `analyst:` (GID candidate forest); no `global:` promotion.
7. **Three-layer no-circularity:** the sprint report's first paragraph states that reproducing a DOMAIN-layer result (a known QFT identity) is not evidence for the UNIVERSAL-layer conjecture that 𝕽 exists.
8. **Pre-registration before certification:** `prereg-002` is committed first; its commit hash goes into every B15 certificate as `prereg_ref`. Thresholds are proposals until Erick confirms (§9-OI-2); if a calibration route exists, use it and name it.
9. **Transcribe, don't reconstruct.** Every kinematic convention, zero-locus index set, split statement, δ-shift prefactor/limit, and the EFT-hedron bound are copied from the cited arXiv source with equation numbers into `prereg-002`. Re-verify each arXiv ID at transcription time; any discrepancy with §12 is filed in `edit-011`.
10. **Determinism:** every generator takes `seed`; certificates record `seed`, `python_version`, `generator_sha256`. Use the CI gate's pinned hash seed convention.
11. **Package:** `gv-b15-surf.zip` (< 25 MB; split `part1/part2` if larger, live tree in part1), with this work order at the root, per `.cursor/rules` work-zip protocol.

---

## Work packages

### WP0 — Reconciliation (½ day)
Produce `docs/sprint-B15-SURF-reconciliation.md`: pre-flight answers; a table mapping each "Create" path above to the actual path used if any convention forced a change; the public signatures you will call — `b1_moment_solver.exact.psd_certificate`, `hankel`, `hankel_rank`; `pir.symbolic.linear.solve`, `verify_solution`, `verify_farkas`; `pir.fingerprints.KnownGrammarDB`, `full_fingerprint`, `specific_fingerprint`; `pir.candidates.apply_rules`, `evaluate`, `lattice_fact`; `canon.engine.canonicalize`; `b1_moment_solver.certificate.save_certificate`. No new abstractions where these suffice.

### WP1 — B15-POS: verifier vs. a published EFT-hedron bound (gate for the sprint)

*Why:* `docs/verifier-ops-v0.1.md` lists Schur pivots as primal certificates; R26 already asked for dual/exclusion functionals; neither has been scored against a bound the program did not build. If the verifier cannot reproduce a known bound with both certificates, its PERMITTED/REJECTED outputs elsewhere inherit `asm:B15-POS-unverified` taint.

- **WP1.1 Target.** From arXiv:2012.15849 select the lowest-order 2→2 single-scalar forward-limit bound that is **two-sided** (not a bare `g ≥ 0`). Transcribe equation number, normalization, Mandelstam conventions into `prereg-002 §POS`. **Stop and file §9-OI-1 for Erick before proceeding.**
- **WP1.2 `dispersive.py`.** Build the truncated moment/Hankel matrices in `Fraction` from the paper's dispersive representation. Truncation order is a `HEURISTIC` choice → `warnings[]` names it.
- **WP1.3 Primal (T1).** Registered interior point → `psd_certificate` chain, all pivots ≥ 0 → `PERMITTED`, `witness = {pivots}`.
- **WP1.4 Dual (T2).** Registered exterior point → (a) exact negative Schur pivot **and** (b) a rational Farkas vector `y` via `pir/symbolic/linear.py` such that `verify_farkas(A, b, y)` is `True` → `REJECTED`, `witness = {pivots}`, `impossibility_certificate = {y, evaluation}`. Both (a) and (b) required.
- **WP1.5 Bracket (T3).** Exact bisection along one registered ray; report `[inner feasible, outer infeasible]` exactly as B1/B7 report `certified_inner_interval`. Registered width threshold (proposal `1/2^20` in the bound's units).
- **WP1.6 Negative control (T4).** Perturb one Hankel entry by a registered rational; verdict must flip with a witness.
- **WP1.7 Standalone re-verify (T5).** `tools/verify_b15_certificate.py` re-derives pivots and re-evaluates `y` from the certificate JSON alone — no import of `b15_surf`.

**Pass:** bracket contains the published value at ≤ threshold; T1–T5 green; certificate validates. **Fail:** file `REJECTED` in `edit-011`, tag `asm:B15-POS-unverified` on B15-ZERO/GID facts, continue WP2–3, and add a `NEEDS-REVIEW` line to the report — do not silently repair.

### WP2 — B15-ZERO: exact ground-truth amplitude dataset

Conventions (transcribe into prereg §ZERO from arXiv:2312.16282 §2): `X_{i,j} = (p_i+…+p_{j−1})²`, `c_{i,j} = X_{i,j} + X_{i+1,j+1} − X_{i,j+1} − X_{i+1,j}`, and the exact hidden-zero locus statement (the "causal diamond / maximal rectangle" condition) **with its index ranges**.

- **WP2.1 Route A (T1, T2).** Catalan-enumerated triangulations (assert 5 for n=5, 14 for n=6, 132 for n=8); `A_n = Σ_T Π_{chords} 1/X`.
- **WP2.2 Route B (T2).** Factorization recursion from `A_3 = 1`, `A_4 = 1/X_{1,3} + 1/X_{2,4}`. A = B exactly on ≥ 200 registered seeded rational points per n ∈ {5,6,8}. Any mismatch → `REJECTED`, dataset not shipped (falsifier F2).
- **WP2.3 Route C (optional, HEURISTIC).** ABHY canonical form (arXiv:1711.09102) for n=5 only if WP2.1–2.2 finish early; tag HEURISTIC until orientation conventions are pinned.
- **WP2.4 Zeros (T3, T4).** For each registered split, solve the locus exactly, evaluate by A and B, assert `A_n == 0`. Negative control: same count of random rational linear constraints → nonzero at ≥ 95% of points (registered).
- **WP2.5 2-split (T5).** Per arXiv:2405.09608: relax exactly one mesh variable `c_⋆` off the locus; assert exact factorization into the two lower currents the paper specifies, ≥ 50 registered points.
- **WP2.6 δ-shift (T6).** Implement the even/odd shift; take the NLSM prescription **exactly as arXiv:2401.05483 states it** (prefactor + limit transcribed into prereg; this work order deliberately does not restate it). Compare with `nlsm_feynman.py` at 4 and 6 points — exact equality.
- **WP2.7 Dataset (T7).** `data/b15/zero_dataset_v1.json`: per theory ∈ {TrPhi3, NLSM}, per n, records `(point, value, on_zero_locus, split_id|null)` as `"p/q"` strings; plus `seed`, `generator_sha256`, hash-stable canonical JSON (`pir/canonical.py`). YM-scaffolded rows only if WP2.8 completes.
- **WP2.8 (stretch, §9-OI-3).** YM-scaffolded amplitudes per arXiv:2401.00041 at the lowest closed-form multiplicity. If skipped, prereg records two-class scope **before** WP3 runs.

**Evidence:** E0, tier BENCHMARK, category *recovered known result*. Certificate `certificate_class: "EXACT-RATIONAL (same class as B1/B2/B7); perturbative-amplitude identities"`.

### WP3 — B15-GID: blind identification from zero/pole fingerprints

*Why:* this is B8 → `KnownGrammarDB` → B12-RGRC machinery applied to a case whose equivalence-class answer is known in advance. Both registered outcomes are results.

- **WP3.1 Fingerprint (T1).** From the label-stripped dataset only (never import `b15_surf.zero`): exact pole set (which `X_{i,j}` diverge under registered scaling), zero-locus set (which registered splits annihilate), 2-split structure, and `redundancy_metric = (#triangulation terms) / (#independent fingerprint constraints)` (P1 — record, don't act). Pole/zero detection is `SOUND` (exact); any "generically nonzero" tolerance is `HEURISTIC`.
- **WP3.2 Quotient (T2).** `canon.engine.canonicalize` under the registered relabeling group (cyclic + reflection of color order) before hashing — R15 rule as code.
- **WP3.3 Identify (T3).** Register TrPhi3 / NLSM / (YM) in a `KnownGrammarDB` from `data/b15/known_grammars_v1.json` (expected pole sets, zero loci, deformation parameter and range, sources). Two-hash semantics as documented in `pir/fingerprints.py`: full = δ-insensitive invariants; specific = includes δ. Output a **candidate forest** via `pir.candidates` (monotonic, provenance-annotated), `selected: null` unless one candidate survives.
- **WP3.4 δ recovery (T4).** Solve δ exactly from shifted pole positions, or `NONIDENTIFIABLE(delta-underdetermined)`.
- **WP3.5 Held-out (T5).** Fit on n ≤ 6, apply to n = 8.
- **WP3.6 Adjudicate (T6)** — pre-registered:
  - **Outcome A (separable):** classes distinguished and δ recovered at n=8, similarity and confidence each ≥ registered threshold (proposal 9/10, reported as separate fields with `correlator` named). Verdict `FORCED` on the grammar label; claims-table `framework capability / E0`.
  - **Outcome B (not separable, because the theories share zeros by construction):** `OBSERVATIONALLY_EQUIVALENT([TrPhi3, NLSM(, YM)] mod δ-shift)`, `selected: null`; memo line: fingerprint-only identification is structurally limited for deformation-related grammars. **This is the more important epistemic result and is not a classifier failure.**
  - **Tie between ≥ 2 menu grammars** with δ recovered: `AMBIGUOUS` (B8 semantics).
  - **Anything else:** `NONIDENTIFIABLE(cause)` + open item.
  - **Falsifier F4:** separable at n ≤ 6 but not n = 8 → `fingerprint_v1` `REJECTED` as over-fit; not shipped.

### WP4 — Design note (no code)
`docs/notes/constraint-manifest-representation-v0.md`, ≤ 3 pages: for B10 (Gaussian channels) or B12-b (CPTP), which of Pos/Uni/Cau is currently *checked* by a pass; candidate coordinates under which it would hold *identically* (u-variable analogue: symplectic normal form, Choi/Stinespring parameterization, Williamson form, …); what each makes automatic, what it hides; and a pre-registerable Stage-3 threshold ("≥ 1 previously-checked constraint becomes an identity, verified on B10 + B15-ZERO data"). Ends with §9-OI-4 for Erick.

---

## 5. Certificate fields added for B15 (`schemas/b15_certificate.schema.json`)

On top of the v0.3 shape (hard constraint 3):

| Field | Type | Required | Note |
|---|---|---|---|
| `soundness` | `SOUND\|HEURISTIC` | yes | weakest of composed passes |
| `warnings` | `[{location, text}]` | if HEURISTIC | located |
| `ground_truth_route` | string | yes | `"planar-feynman-sum"`, `"published-bound:arXiv:2012.15849 eq.N"` — the honesty marker for the report's ground-truth disanalogy |
| `witness` | object | PERMITTED/REJECTED | pivots / split factorization (PIR fact schema already reserves it) |
| `impossibility_certificate` | object\|null | REJECTED in POS | Farkas `y` + evaluation (field already reserved in `pir/schema/fact.schema.json`) |
| `similarity`, `confidence`, `correlator` | `"p/q"`, `"p/q"`, string | GID | separate, never merged |
| `evidence_level`, `pir_level`, `layer` | `E0..E4`, `L0..L3`, `DOMAIN` | yes | orthogonal axes |
| `assumptions` | `["asm:…"]` | yes (may be empty) | ADR-0002 taint |
| `prereg_ref` | string | yes | `prereg-002@<commit>` |
| `falsifier_direction` | string | yes | what flips this verdict |
| `seed`, `python_version`, `generator_sha256`, `certificate_id` | — | yes | determinism + content address |

Negative tests: non-SPEC verdict rejected; HEURISTIC + E0 rejected; missing `impossibility_certificate` on a POS `REJECTED` rejected.

## 6. PIR emission
Each certificate is also emitted as a `pir.Fact` (b13 bridge pattern): B15-POS and B15-ZERO facts in `domain:` at `L2/E0`; the GID candidate forest as `pir.candidates.lattice_fact` in `analyst:` with `candidate_class` tags only (never as verdicts). Any `asm:` taint from WP1 failure propagates by the existing invalidation traversal — no bespoke logic.

## 7. Pre-registration — `docs/preregistrations/prereg-002-surfaceology-benchmarks.md`
House style of prereg-001. Sections: **Frozen machinery** (paths + commit); **POS** (bound, conventions, interior/exterior points, ray, width threshold, perturbation); **ZERO** (conventions, split list per n, point counts, nonzero-rate threshold, δ prescription, NLSM source); **GID** (feature keys, relabeling group, thresholds, held-out n, outcomes A/B/AMBIGUOUS/NONIDENTIFIABLE); **Registered falsifiers** F1–F4 (§8); **Stop rules**; **Evaluation protocol** (results to `sprint-B15-SURF-report.md`, PASS/FAIL per test, both GID outcomes filed as results). Status `DRAFT-FOR-FREEZE` → ACTIVE at the commit Erick pushes.

## 8. Registered falsifiers
| ID | Trigger | Consequence |
|---|---|---|
| F1 | POS bracket excludes the published value | verifier `REJECTED` for external bounds; `asm:B15-POS-unverified` on all B15 facts; program-level review line in report |
| F2 | Route A ≠ Route B anywhere | generator `REJECTED`; no dataset shipped |
| F3 | a registered zero does not vanish exactly | transcription erratum in `edit-011`, or (if source-faithful) open item — never a numerical patch |
| F4 | GID separable at n≤6, not at n=8 | `fingerprint_v1` `REJECTED` (over-fit) |

## 9. Open items — Erick decides; agent does not resolve
- **OI-1** Confirm the EFT-hedron target bound the agent proposes from arXiv:2012.15849.
- **OI-2** Confirm/replace the two proposed thresholds (bracket width `1/2^20`; GID similarity & confidence `9/10`) or name calibration routes.
- **OI-3** WP2.8 YM-scaffolded: required or stretch.
- **OI-4** Open B16 on the constraint-manifest representation (WP4 output)?
- **OI-5** Approve the three claims-table rows before they are appended.
- **OI-6** Should pre-existing certificates be annotated `ground_truth_route: null` in a follow-up (schema-level, no verdict change)? Recommendation: yes, separate sprint.
- **OI-7** Sprint-ID/vocabulary/edit-number pre-flight results if they differ from this order's assumptions (B15, §4 vocabulary without `ENTAILMENT_*`, edit-011, R36).

## 10. Backlog (recorded, not scheduled)
- ML conjecture generation over `zero_dataset_v1` (numeric→symbolic, PySR-style; every proposal must pass `tools/verify_b15_certificate.py` or be filed `REJECTED`/`NONIDENTIFIABLE`). Precedents R38–R40.
- Loop-level zeros / 2-splits (arXiv:2604.13810); cosmological-wavefunction zeros (arXiv:2503.23579) as a second domain for the ≥2-domains rule.
- B11 composition classifier (roadmap) can absorb GID fingerprint keys.
- MnemesisOS transfer of the candidate-forest pattern — separate session.

## 11. Definition of done
- `python3 ci/run_all_certified.py` green, **unmodified**, with `test_b15_pos.py`, `test_b15_zero.py`, `test_b15_gid.py` picked up by its existing glob; run-manifest hash recorded in the report alongside the WP0 baseline hash; zero degradations.
- All three certificates validate against `schemas/b15_certificate.schema.json` and pass `tools/verify_b15_certificate.py` standalone.
- `prereg-002` committed **before** the first certified run; its hash in every certificate.
- `edit-011` memo filed (no cell change); R36–R40 appended; README section; sprint-drops row.
- `sprint-B15-SURF-report.md`: no-circularity paragraph first; per-test PASS/FAIL; which GID outcome obtained; redundancy metrics; errata; §9 items reproduced with status.
- `gv-b15-surf.zip` < 25 MB with this order at root.

## 12. Ledger entries to append (`docs/references.md`, R36–R40; verify IDs at transcription)
| R | Reference | Role |
|---|---|---|
| R36 | Arkani-Hamed, Frost, Salvatori, Plamondon, Thomas, *All Loop Scattering as a Counting Problem*, arXiv:2309.15913; *All Loop Scattering for All Multiplicity*, arXiv:2311.09284 | METHOD — curve-integral formalism; representation-redundancy pattern (P1/P7) |
| R37 | Arkani-Hamed, Cao, Dong, Figueiredo, He, *Hidden zeros…*, arXiv:2312.16282; *NLSM ⊂ Tr(φ³)*, arXiv:2401.05483; *Scalar-scaffolded gluons*, arXiv:2401.00041; Arkani-Hamed & Figueiredo, *All-order splits*, arXiv:2405.09608 | TARGET (B15-ZERO ground truth) + METHOD (zeros/splits as fingerprints, δ-shift as base-grammar+deformation) |
| R38 | Arkani-Hamed, Bai, He, Yan, *Scattering forms and the positive geometry of kinematics…*, arXiv:1711.09102 | METHOD — ABHY associahedron (Route C) |
| R39 | Cheung, Dersy, Schwartz, arXiv:2408.04720; Moynihan, arXiv:2602.15169 | METHOD — ML/symbolic rediscovery precedents (backlog); rediscovery only |
| R40 | Guevara, Lupsasca, Skinner, Strominger, Weil, arXiv:2602.12176 | METHOD — AI-conjectured, human-proved amplitude result; narrow slice; preprint |

R26 (EFT-hedron) already covers WP1; cite it, do not duplicate.
