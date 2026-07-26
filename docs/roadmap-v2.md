# Roadmap v2 — Proving "domain-independent, infer, same, genuinely shared" (2026-07-14)

Successor to sprint-specs-v1 (S1–S4 DONE; B1–B10 landed). Organizing
question, adopted from the cross-agent assessment: the next phase must
show the vocabulary *selects the right theory over alternatives* and
*predicts something its designers did not already know*.

## 0. The four words → four machines

| Word | Machine | Deliverable | Acceptance test |
|---|---|---|---|
| **domain-independent** | **M1 Canonicalization engine** (`canon/`): every domain artifact → invariant feature vector (positivity-cone data, rank sequence, symmetry group, factorization poset, flow exponents). | Layer-2 of the L0–L4 architecture; REPRESENTATION_DEPENDENT verdicts automatic. | Two chart-presentations of the same model yield identical features (unit-tested on B8/B9 generators); a coordinate change never flips a verdict. |
| **infer** | **M2 Generator engine** (`generator/`): typed expression grammar + dimensional analysis + sparsity/MDL ranking, feeding the EXISTING certified verifiers (B1/B2/B6/B10 gates). Replaces B3's supplied templates. | Candidate laws with complexity budgets and certificates. | B3 and B5 rediscovered from raw tables with NO template hints; false-law rate measured on scrambled controls. |
| **same** | **M3 Equivalence layer**: SPEC §5 implemented — outputs are equivalence classes; `docs/no-go.md` ledger of impossibility results (six entries seeded from the assessment §12). | Engine emits OBSERVATIONALLY_EQUIVALENT(class) and NONIDENTIFIABLE(cause) with machine-readable causes. | Dual presentations (Hamiltonian/Lagrangian toy; gauge pair; 1D-reduced B9) grouped into one class, never resolved. |
| **genuinely shared** | **M4 B12-RGRC** (rival-grammar recovery challenge, renamed from the assessment's "B10"): 8 generator families incl. GPT-beyond-quantum, shared-latent causal model, and spurious-shared-calibration. + **M6 preregistration protocol**. | Blind verdicts: {genuine shared latent | shared mathematics | shared apparatus | shared units | accidental low rank | independent}. | Confusion matrix with registered thresholds; false-positive rate on the spurious-calibration family is the headline number. |

Supporting machines: **M5 coarse-graining engine** (𝖱 interface unifying
B7/B8/B9: semigroup tests, fixed points, restraint degradation) and
**M7 measurement-interface layer** (`measurement_interface/`: apparatus
transfer functions, censoring, calibration priors, injection-recovery —
generalizing B9-T3 and the S4 APPARATUS_LIMITED lesson).

## 1. The prospective prediction (M6, Priority top)

Freeze rules at tag `v0.7-frozen`; preregister in `preregistrations/`
(timestamped, hash-signed) ONE quantitative prediction on data unavailable
during development. Primary candidate: **GWTC-next area-theorem
population** — predicted η_A distribution parameters and per-event CP-bound
pass rates from the frozen B4 pipeline, registered before the next LVK
catalog release. Secondary: predicted α_s from the *next* independent
extraction class (frozen S3 pipeline). Evaluation against registered
thresholds; failure filed as program-level REJECTED per SPEC §6.

## 2. The experiment puzzle-chain (each solves a prerequisite of the next)

Chain rule: the certificate of experiment k becomes a mandatory
`m_layer_stipulations` input of experiment k+1. Solving one unlocks the
next — the puzzle structure is explicit and auditable.

1. **EXP-D — α-consistency (zero hardware, S3-EM).** DONE (5/5):
   `s3em_alpha/` + `tests/test_s3em.py` certify the live 5.5σ Rb–Cs
   tension, localize it to h/m (shared chain ≤7%), refuse single-value
   FORCED, emit NONIDENTIFIABLE(systematics vs new physics). Memo:
   edit-009. **Cells: α Pos/Uni memo (no change).**
2. **EXP-B — FDT audit (tabletop condensate/optomech).** Uses EXP-D's
   combination protocol for its thermometry channels; B9 gates
   (independent calibration route mandatory).
   → *Output:* certified noise-and-dissipation floor. **Cells: k_B
   Thm/Pos (ANALOGUE-REAL).**
3. **EXP-A — CV channel tomography (B10 REAL).** Requires EXP-B's noise
   certificate as its stipulated floor; reconstructs (T, N); Caves
   boundary and CP gate at E2.
   → *Output:* certified channel library. **Cells: ħ Pos/Uni.**
4. **EXP-E — composition classifier (B11 REAL).** Trains/negative-controls
   on EXP-A's library; blind product-vs-joint verdicts with AMBIGUOUS.
   → *Output:* calibrated composition detector. **Cells: ħ Cmp.**
5. **EXP-C — information-Gram across an analogue horizon.** Requires
   EXP-A tomography (Gram entries) + EXP-E discrimination (kinematic vs
   informational structure); measures the S4 signature: kinematic rank
   ∝ system size vs information-Gram rank ℓ-independent ∝ ln(1/ε).
   → *Output:* rank-vs-entropy law in the lab. **Cells: ?₇ typing
   (ANALOGUE-REAL, E4→E2 within the analogue).**
6. **EXP-F — two-site known-hidden-variable benchmark** (assessment §17,
   adopted verbatim): controlled shared parameter, blinded engine, FP/FN
   rates measured — the instrument for "genuinely shared," validated on a
   KNOWN shared variable before any residual search. Requires the full
   upstream chain. **Cells: ?₈ typing methodology.**

## 3. Priorities (aligned)

P0 claims table (SPEC §6) retrofitted to README — one sprint.
P1 SPECIFICATION.md (DONE this sprint) + ontology/ schemas.
P2 B12-RGRC (M4) — DONE core (4/4) + B12-b quantum layer (5/5); GPT state-space families deferred to B12-c.
P3 M1 canonicalization — DONE (5/5) core (`canon/`); broader gauge
   families deferred to M1-b. **M7 measurement layer — DONE (5/5)**
   (`measurement_interface/`; SPEC §2 layer 3 executable — apparatus transfer,
   dual calibration routes, APPARATUS_LIMITED gate, censoring correction;
   generalizes B9-T3 + the S4 ε-rank plateau). The gate before any new
   cosmology/biology/BEC claim is now in place; M7-b (continuous transfers,
   calibration priors, real-posterior recovery, PIR analyzer coupling) deferred.
P4 M2 generator engine — DONE (5/5) monomial core (`generator/`); B3/B5
   rediscovered from raw tables; wider grammars deferred to M2-b.
P5 M6 preregistration: prereg-001 DRAFT-FOR-FREEZE filed; ACTIVE at tag `v0.7-frozen`.
P6 External replication pack (`external_replication/`) + CI hardening.

Benchmark namespace going forward: B11 = composition classifier;
**B12 = RGRC**; B13+ reserved by roadmap only (no ad-hoc numbering).
