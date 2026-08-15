# Constant-Atlas Physics AI Research

**Observational decompilation of physics.** Treat physical law as an unknown rule system; recover its grammar — types, invariants, operators, composition rules — one certified constraint at a time. Constants are conjectured to be invariants of a positive, compositional operator structure 𝕽, not primitive numbers.

Core documents:
- [`docs/constant-atlas-v0.5.md`](docs/constant-atlas-v0.5.md) — constant atlas, restraint matrix, bridge formula, discipline formula, decoding chain (the base algorithm). *(Prior: [v0.4](docs/constant-atlas-v0.4.md), [v0.3](docs/constant-atlas-v0.3.md), [v0.2](docs/constant-atlas-v0.2.md).)*
- [`docs/conjectures-v0.1.md`](docs/conjectures-v0.1.md) — the ten blank-filling conjectures ?₁–?₁₀, each with formalization, test protocol, and falsifier.

## B1 — Truncated moment solver (implemented ✅)

Benchmark B1 of the decoding chain: hidden moments of a finite measure are recovered through positivity + Curto–Fialkow flat extension, with **exact rational certificates**.

Epistemic outcomes (Discipline Formula v0.2):

| Outcome | Meaning | Mechanism |
|---|---|---|
| `FORCED` | uniquely determined | flat extension → exact recurrence, exact recovery |
| `PERMITTED` | constrained to a certified interval | exact-PSD bisection; inner interval certified feasible, outer bracket certified infeasible; unbounded directions reported honestly |
| `REJECTED` | hard constraint violated | exact negative Schur pivot as witness |

Design guarantees:
- All positivity/rank/forcing claims use **exact `Fraction` arithmetic** (the Schur pivots generalize the `G00 > 0`, `S2 = G22 − G02²/G00 > 0` verifier standard).
- Numerics are **quarantined** to atom extraction only, with an audited residual against the exact moments reported in the certificate.
- `Score = ∞` on any violated **or uncertified** hard constraint — uncertified ≠ pass.

Run:
```bash
python3 tests/test_b1.py
```
Output: 4/4 tests pass; JSON certificate written to `certificates/b1_certificate.json`.

```
b1_moment_solver/
  exact.py        exact Hankel matrices, PSD certificates, rank, flat-extension recurrence
  recover.py      forced-tail recovery, permitted intervals, audited atom extraction
  certificate.py  forced/permitted/free/rejected JSON certificate format
tests/test_b1.py  T1 forced · T2 permitted · T3 rejected (negative control) · T4 atom audit
```

## B4 — Area theorem as gravity's composition law (implemented ✅, REAL mode)

Benchmark B4 tests conjecture **?₃**: A_f ≥ A₁ + A₂ under merger, via the dimensionless invariant η_A = (A_f − A₁ − A₂)/A_f. Certificate class: **STATISTICAL** (Monte Carlo + exact Clopper–Pearson binomial lower bound) — deliberately distinct from B1's exact-rational class.

```bash
python3 tests/test_b4.py          # DEMO regression
python3 scripts/run_b4_real.py    # REAL mode (see data/README.md)
```

**REAL results** (official posteriors; certificate `certificates/b4_certificate.json`, mode REAL):
GW250114 → SUPPORTED (P_lower 0.99977, η_A ≈ 0.366) · GW150914 → SUPPORTED (P_lower 0.99862, η_A ≈ 0.363). Final-state provenance in these PE releases is NR-fit (ringdown columns absent); loader still auto-prefers ringdown keys when present. Edit-001 apply conditions are checked; atlas bump to v0.3 remains for the promotion pass.

Demo regression (reconstructed summary-statistic posteriors): GW150914 / GW250114 SUPPORTED; GW190521 / GW190814 INCONCLUSIVE **by design**; falsifier injection → VIOLATION_CANDIDATE.

```
b4_area_pipeline/
  kerr.py       Kerr horizon area, η_A invariant
  remnant.py    validated nonspinning NR fits (cross-check layer only)
  demo_data.py  published summary stats, provenance-flagged
  pipeline.py   MC test, exact binomial bounds, statistical certificate
  loader.py     GWOSC/PESummary HDF5 loader (REAL mode, ringdown-preferred)
```

## B2 — Qubit process completion via Choi positivity (implemented ✅)

Extends the exact-rational certificate class to **Gaussian-rational Hermitian matrices** (`cexact.py`) and demonstrates **restraint stacking** — the atlas mechanism — on a finite quantum process:

| Restraints active | Hidden Choi entry | Outcome |
|---|---|---|
| PSD alone | diagonal of rank-2 mixed channel | PERMITTED — certified interval [1/9, ∞) |
| PSD + trace preservation | same entry | **FORCED** to exact truth (1) |
| PSD + rank-1 (pure process) | complex off-diagonal | **FORCED** to −12i/25 exactly (flat-extension analogue) |

Plus exact CPTP audit (Choi's theorem as the gate) and a negative control (corrupted Choi rejected with exact witness pivot). Run: `python3 tests/test_b2.py` — 5/5 pass, certificate in `certificates/b2_certificate.json`.

## B3 — Electroweak closed system (implemented ✅)

Benchmark B3 runs the decoding chain on tree-level electroweak structure. Pseudo-data for `(e, g, g′, θ_W, M_W, M_Z, v, G_F)` is presented as an unlabeled exact-Fraction table; the engine must *discover* the relations by template forcing (same epistemic class as B1/B2) and report **d_identifiable = 3**.

```bash
python3 tests/test_b3.py
```

Results (5/5):
- Headline relations found: `e = g sinθ_W`, `M_W = gv/2`, `M_Z = (v/2)√(g²+g′²)`, `G_F = 1/(√2 v²)` (√2 certified via `1/G_F² = 2 v⁴`)
- Exact Jacobian rank over ℚ: n=9, rank=6 → **d_identifiable = 3**
- FREE basis recovered: `{g, g′, v}`; forced: `{e, sinθ_W, cosθ_W, M_W, M_Z, G_F_inv_sq}`
- Held-out `M_Z` predicted exactly from the free basis; corrupted `e` breaks R1/R2/R9 (negative control)

```
b3_electroweak/
  exact.py        Fraction helpers, Pythagorean Weinberg angle
  pseudo_data.py  closed-system generator from free basis (g, g′, v)
  relations.py    exact template discovery (R1–R9)
  rank.py         Jacobian rank → d_identifiable
  graph.py        FREE/FORCED dependency labels
  discover.py     end-to-end pipeline
```

## B5 — Cluster / multi-channel factorization (?₅ → H) (implemented ✅)

Benchmark B5 is the B3 extension named in conjecture ?₅: 2→2 and 2→4 toy channels with exact-Fraction amplitudes. The engine must discover that **one coupling** controls all factorization residues (`R_24 = A_22²`) and reject channel-dependent couplings.

```bash
python3 tests/test_b5.py
```

Results (5/5): factorization identities hold; **d_identifiable = 1**; falsifier FAIL as required. Atlas edit-002 promotes gauge Cmp (?₅) → **H** for α, α_s, and the electroweak block.

```
b5_cluster/
  exact.py        Fraction helpers
  pseudo_data.py  local vs channel-dependent counterexample tables
  relations.py    factorization template discovery (F1–F6)
  rank.py         Jacobian rank → d_identifiable
  discover.py     end-to-end pipeline
```

## B6 — 2D QNEC as Schur-pivot positivity (?₂ → P) (implemented ✅)

Benchmark B6 places the strengthened 2D QNEC inside the certified verifier: `2π⟨T⟩ ≥ S'' + (6/c)(S')²` **⇔** PSD of `M = [[c/6, S'], [S', 2πT − S'']]`. Saturation ⇔ rank-1 flat boundary (same structure B1 uses for forcing).

```bash
python3 tests/test_b6.py
```

Results (5/5): 500-case equivalence audit; vacuum saturation (pivot exactly 0); coherent-state strictness; thermal symbolic identity; negative control (witness pivot −1/12). Atlas edit-003 promotes G Pos (?₂) → **P (H-track)**. Epistemic scope: exact reformulation + certified instances — not a new QNEC proof; GNS realization of M remains open for H.

```
b6_qnec/
  qnec.py   matrix, Schur pivot, vacuum / coherent / thermal / falsifier instances
```

## B7 — Onsager transport-matrix completion (implemented ✅)

Benchmark B7 (seeded by R15) exercises restraint stacking on **effective coefficients**: near-equilibrium transport `J = L X` with `L ⪰ 0` (second law) and `L = Lᵀ` (Onsager reciprocity).

```bash
python3 tests/test_b7.py
```

Results (5/5): ground-truth PD+symmetric; PSD alone → PERMITTED; +reciprocity → FORCED; second-law and reciprocity falsifiers rejected. Certificate: `certificates/b7_certificate.json`.

```
b7_onsager/
  exact.py         Fraction helpers
  pseudo_data.py   thermoelectric L, hide/corrupt
  complete.py      PSD audit, reciprocity forcing, permitted intervals
```

## B8 — Blind grammar identification (implemented ✅) & Atlas Engine v0.1 (implemented ✅)

**B8** (6/6): identifies Hamiltonian / gradient-flow / GENERIC / quantum grammars from unlabeled data using similarity-invariant signatures only (R15 canonicalization). Exactly-zero damping is never claimed — conservative verdict is **HAMILTONIAN_WITHIN_BOUND** with a certified damping bound; drift test detects weak damping below the spectral floor. Quantum split is exact via Choi rank (B2).

```bash
python3 tests/test_b8.py
python3 tests/test_atlas_engine.py
```

**Atlas Engine v0.1** (5/5): restraint matrix as a constraint-satisfaction object. Eight theorem-anchored implication rules propagate to fixpoint → **7 candidate upgrades** (Cmp/Top/Λ-Pos), each with derivation chains; engine cannot mint H. Certificates: `certificates/b8_certificate.json`, `certificates/atlas_engine_certificate.json`.

**Sprint work orders:** [`docs/sprint-specs-v1.md`](docs/sprint-specs-v1.md) — S1–S4 DONE; S2 remains PARTIAL on the GNS holdout checklist. Successor roadmap: [`docs/roadmap-v2.md`](docs/roadmap-v2.md).

## B9 — Circuit decompilation: 2025 Nobel methodology as gates (implemented ✅)

Benchmark from the Clarke–Devoret–Martinis paper chain (ledger **R16–R22**). Linearized RCSJ in flux–charge coordinates = GENERIC grammar; B9 (6/6) runs the Nobel recipe as certified gates: M/W recovery + FDT; hidden-bath REJECT; **Gibbs circularity** (independent calibration required); model-order inversion; effective-T prohibition; quantum held-out spectroscopy.

```bash
python3 tests/test_b9.py
```

Note: [`docs/notes/nobel2025-circuit-decompilation.md`](docs/notes/nobel2025-circuit-decompilation.md). Certificate: `certificates/b9_certificate.json`.

```
b9_circuit/
  model.py      RCSJ ground truth + simulation
  estimate.py   equilibrium decompilation (K = −J_b H_V⁻¹)
  quantum.py    charge-basis H + held-out spectroscopy
```

## Atlas edits

- [edit-001](docs/atlas-edits/edit-001-conjecture3-promotion.md) — **?₃ → P (H-track)** (atlas v0.3)
- [edit-002](docs/atlas-edits/edit-002-conjecture5-promotion.md) — **?₅ → H** (atlas v0.4)
- [edit-003](docs/atlas-edits/edit-003-conjecture2-qnec.md) — **?₂ → P (H-track)** (atlas v0.5)
- [edit-004](docs/atlas-edits/edit-004-cluster-cmp-extension.md) — gauge_qft Cmp → P (atlas v0.6)
- [edit-005](docs/atlas-edits/edit-005-weak-conjecture6-top.md) — weak ?₆ Top → P (atlas v0.6)
- [edit-006](docs/atlas-edits/edit-006-weak-conjecture7-and-findings.md) — weak ?₇ + Sym/Uni findings (atlas v0.6)
- [edit-007](docs/atlas-edits/edit-007-gns-partial-draft.md) — **DRAFT** GNS partial for ?₂ (no atlas bump)
- [edit-008](docs/atlas-edits/edit-008-s3-empirical-note.md) — S3 empirical memo (no cell change)
- [edit-009](docs/atlas-edits/edit-009-alpha-tension-memo.md) — α-row tension memo (no cell change)

Current atlas: [`docs/constant-atlas-v0.6.md`](docs/constant-atlas-v0.6.md).

## S2 — GNS program for ?₂ (partial ✅)

Two complementary layers; ?₂ stays **P (H-track)**.

**Operator probe** (6/6): continuum vacuum Gram ray exact; geometric (modular weight, cut stress) pair fails ratio tracking. `certificates/s2_gns_certificate.json`.

**State layer** (5/5): free Dirac Gaussian kernels populate B6 from data — c_fit = 1.0001 blind, vacuum flat-boundary saturation, thermal energy-density identity (1.9%), entanglement first law (4× ε-scaling), positivity gate. Obstruction localized to operator ansatz (missing bilocal modular term). `certificates/s2_certificate.json`.

Notes: [`docs/notes/s2-gns-free-fermion.md`](docs/notes/s2-gns-free-fermion.md), [`docs/notes/s2-gns-status.md`](docs/notes/s2-gns-status.md).

```bash
python3 tests/test_s2_gns.py
python3 tests/test_s2_state_layer.py
python3 tests/test_s2b.py
```

## S3 — one α_s across channels and scales (implemented ✅, 5/5)

`s3_pdg/` + `tests/test_s3.py`: four measurement classes spanning 51× in scale combine to α_s(M_Z) = 0.1186 ± 0.0009; no-running null rejected ~13σ. Certificate: `certificates/s3_certificate.json`. Memo: edit-008.

## S4 — ?₇ rank saturation + information-Gram refinement (implemented ✅, 5/5)

`s4_ds/` + `tests/test_s4.py`: horizon-weighted rank saturates; unweighted control grows; information Gram `D C(1−C) D` has ℓ-independent saturated rank. Certificate: `certificates/s4_certificate.json`. BEC bridge: [`docs/notes/experimental-program-bec.md`](docs/notes/experimental-program-bec.md).

## B10 — CV Gaussian-channel completion (implemented ✅, 6/6) + program v2

**B10** (`b10_cv_channel/`, EXACT + tomography): CP gate `M = N + (i/2)(Ω − TΩTᵀ) ⪰ 0`; Caves amplifier bound as a certified PERMITTED boundary; symplectic forcing; quantum-limited attenuator semigroup exact; tomography + held-out squeezed prediction. ħ Pos/Uni/Cmp at BENCHMARK tier — EXP-A software back-end. Schema: [`schemas/bec_experiment_record.schema.json`](schemas/bec_experiment_record.schema.json).

```bash
python3 tests/test_b10.py
```

**Program v2:** [`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) — formal 𝕽, three-layer separation, evidence levels E0–E4, extended outcomes. [`docs/notes/related-programs-review.md`](docs/notes/related-programs-review.md) — bootstrap/EFT-hedron, GPT, categorical QM, equation discovery (R25–R30). [`docs/roadmap-v2.md`](docs/roadmap-v2.md) — machines M1–M7, B12-RGRC, experiment puzzle-chain EXP-D→B→A→E→C→F. Claims: [`docs/claims-table.md`](docs/claims-table.md).

## B12-RGRC — rival-grammar recovery challenge, core (implemented ✅, 4/4)

Six generator families (genuine shared latent, calibration artifact, sensor common mode, shared-mathematics trap, direct coupling, independence). **18/18** blind classification; **0/24** false `GENUINE_SHARED_LATENT` on look-alikes; passive F1 vs F6 → `NONIDENTIFIABLE`; clamp interventions separate them; held-out intervention prediction <12%. GPT/CPTP process families: see B12-b below. Certificate: `certificates/b12_certificate.json`.

```bash
python3 tests/test_b12.py
```

## B12-b — rival grammars, quantum layer (implemented ✅, 5/5) + prereg-001

Process families G1–G4 (QUANTUM_CPTP / CLASSICAL_EMBEDDABLE / BEYOND_QUANTUM / NOT_A_PROCESS) + CHSH ladder. **12/12** blind; **0/12** false BEYOND_QUANTUM; population-only access → `NONIDENTIFIABLE(apparatus)`; held-out channel prediction trace distance < 0.02. Certificate: `certificates/b12b_certificate.json`.

**M6:** [`docs/preregistrations/prereg-001-gwtc-next-area.md`](docs/preregistrations/prereg-001-gwtc-next-area.md) — DRAFT-FOR-FREEZE prospective registration (P1–P3 GWTC-next + α_s P4). Activates at tag `v0.7-frozen`.

```bash
python3 tests/test_b12b.py
```

## S3-EM — α-row tension, certified (implemented ✅, 5/5, FUNDAMENTAL-REAL)

Puzzle-chain step 1 (EXP-D). Same pipeline that FORCED α_s (S3) here **refuses** a single-value α verdict: Rb 2020 / Cs 2018 / a_e → χ²/dof ≈ 17, **5.5σ** Rb–Cs pull; tension localized to h/m (shared chain ≤7%); CODATA ×2.5 softens but does not resolve. Cause: `NONIDENTIFIABLE(systematics vs new physics)`. Memo: edit-009; ledger R31–R33. Certificate: `certificates/s3em_certificate.json`.

```bash
python3 tests/test_s3em.py
```

## M1 — canonicalization engine (implemented ✅, 5/5)

`canon/`: chart-invariant feature hashes (R15 rule as code). Similarity / unitary / congruence tests across B9 generators, B12-b Choi, S2/S4 Gram; cross-domain spectral class without collapsing factor_poset. Broader gauges deferred to M1-b. Certificate: `certificates/m1_certificate.json`.

```bash
python3 tests/test_m1.py
```

## M2 — typed candidate generator (implemented ✅, 5/5)

`generator/`: monomial grammar + MDL + M1-style canonical dedup. Rediscoveries from raw tables: G_F ~ e² M_W⁻² s_W⁻² (B3) and R_24 ~ A_22² (B5); scrambled controls reject noise-fitting. Engine proposes only — verifiers dispose. Wider grammars deferred to M2-b. Certificate: `certificates/m2_certificate.json`.

```bash
python3 tests/test_m2.py
```

## M7 — measurement-interface layer (implemented ✅, 5/5)

`measurement_interface/`: SPEC §2 **layer 3** (the declared 𝖬 in `y = 𝖬[𝒰(x)] + ε`) as executable code — the gate the roadmap requires before any new cosmology/biology/BEC claim. Generalizes two one-off results: **B9-T3** (an independent spectroscopy calibration route detects a hidden bath that the self-consistent Gibbs route is fooled by and cannot promote — the Voss–Webb ambiguity, same FDT residual threshold 0.15) and the **S4** ε-rank plateau (a `rank ≥ 3` verdict flipping inside the detector's declared resolution window is refused as SPEC-§4 `APPARATUS_LIMITED`, naming the bounding threshold). Also: unbiased injection–recovery with a certified interval, and declared-truncation censoring correction. Honesty invariant (T5): an ample apparatus is never cried apparatus-limited. Certificate carries first-class `calibration_route` + `m_layer_stipulations`. Scope: [`docs/M7-measurement-interface-scope.md`](docs/M7-measurement-interface-scope.md). Certificate: `certificates/m7_certificate.json`. M7-b (continuous/functional transfers, calibration priors, recovery over real B4/B10 posteriors, PIR analyzer coupling) deferred.

```bash
python3 tests/test_m7.py
```

## Truncation-audit probe + GfE review (implemented ✅, 5/5)

Deep-research sprint on Gravity-from-Entropy (Bianconi, R34–R35) with the truncation critique **converted into an executable gate**: on an exact Gaussian ladder, the truncated "emergent constant" differs from the full Schur elimination by an exact rational shift (**−6/3311** in the demo) unless the boundary coupling vanishes — then robustness is *certified* (Fractions-exact), not assumed; the shift scales quadratically in the neglected coupling (**4.002×** under halving). Three verdicts wired to SPEC §§2/4: TRUNCATION_ROBUST / TRUNCATION_ARTIFACT / **TRUNCATION_UNAUDITED — GfE's current status in our classification**, pending their consistency proof. **Mathematics imported:** GfE's metric-pair invariant spectrum + the Burg/Stein divergence (the GQRE's eigenvalue-log core) now live in M1 as `canonicalize(kind="metric_pair")`, congruence-invariance tested. Matrix mapping: Λ×Pos rival typing (vs ?₇/?₈), G×RG rival for ?₄, ?₃-row theory-side note, B12-c family spec. Note: `docs/notes/gfe-bianconi-review.md`.

## PIR — evidence substrate (Stage 1, implemented ✅)

**PIR (Physics Intermediate Representation) is an evidence substrate, not a proof engine.** It records what was measured (L0), what was done (L1), what was certified (L2), and which grammars remain in play (L3) — each fact carrying **orthogonal** representation (L0–L3) and warrant (E0–E4) coordinates, full provenance, and assumption taint. **The atlas engine never reads raw datasets directly; it reads certified PIR facts.** The mathematics that decides verdicts stays in B1–B12 and the frozen [verifier ops](docs/verifier-ops-v0.1.md); PIR is where their inputs, outputs, and dependencies are stored honestly.

Stage 1 delivers the substrate only — **no atlas edits, no changes to any B1–B12 verdict or certificate**. Enforced invariants: append-only facts (invalidation downgrades transitively, deleting nothing); verdict vocabulary locked to SPEC §4/§6 (candidate-class labels are hypothesis taxonomy, never verdicts); SOUND/HEURISTIC honesty (HEURISTIC E3/E4 must carry located warnings); typed cross-namespace transforms (no silent promotion); hash-stable canonical JSON (exact `Fraction`s as `"p/q"`, never floats).

```
pir/            types · models · namespaces · canonical · provenance · passes · schema/
ci/run_all_certified.py   reruns tests/test_b*.py, diffs regenerated vs committed
                          certificates (jitter-tolerant, degradation-strict),
                          signed run manifest; exits nonzero on any degradation
docs/pir-specification-v0.1.md · docs/verifier-ops-v0.1.md · docs/adr/ADR-PIR-0001.md
examples/pir/minimal_circuit.json   B9-style circuit → L0/L1 → exact L2 fact →
                                    OBSERVATIONALLY_EQUIVALENT hypothesis pair
```

Run: `python3 tests/test_pir_models.py` · `test_pir_provenance.py` · `test_pir_namespaces.py` · `test_pir_passes.py` (positive + the six required negative tests) and `python3 tests/test_ci_guard.py` (CI green on the repo, red on a corrupted certificate). Substrate contract: [`docs/pir-specification-v0.1.md`](docs/pir-specification-v0.1.md).

**Stages 2–3 — analysis + forward loop (implemented, additive):** B9 lowering + Circuit Domain Semantics (4 contracts) + structural graph ([`pir/domains/`](pir/domains/)); analyzer runtime + MeasurementProvenance/ObservationalEquivalence/GlobalVariableCandidate analyzers ([`pir/runtime.py`](pir/runtime.py), [`pir/analyzers.py`](pir/analyzers.py)); exact-rational symbolic constraint bridge with SAT witness + Farkas UNSAT certificate ([`pir/symbolic/`](pir/symbolic/)); candidate lattice + GVAR rules and two-hash fingerprints ([`pir/candidates.py`](pir/candidates.py), [`pir/fingerprints.py`](pir/fingerprints.py)); forward recompilation with held-out comparison ([`pir/forward.py`](pir/forward.py)); intervention search / OED ([`pir/intervention_search.py`](pir/intervention_search.py)); and a BEC domain + cross-domain PIR-Diff ([`pir/domains/bec.py`](pir/domains/bec.py), [`pir/diff.py`](pir/diff.py)). Run: `python3 tests/test_pir_b9_lowering.py` · `test_pir_runtime.py` · `test_pir_symbolic.py` · `test_pir_candidates.py` · `test_pir_s3.py`. Conditional results ride as assumption-taint (ADR-0002). No B9/benchmark code, verdict, certificate, or atlas cell is changed; only the P9 workbench UI is deferred. Progress table: [`docs/pir-specification-v0.1.md §9`](docs/pir-specification-v0.1.md).

## B13-CDL — Conditional-Dependency Ledger (implemented ✅, v0.1)

Turns the Conditional-Dependency Survey (v1.0 + v1.1 addendum) into a **machine-readable ledger + certified pipeline suite** (Survey Recommendation Stages 1–5). Self-contained, **stdlib-only**, under [`b13_cdl/`](b13_cdl/): 25 ledger entries (A1–A15, M1–M10) with block×predicate tags, 15 registered falsifiers (liveness computed, not asserted), and 8 pipelines emitting schema-validated certificates. Certified results include: `ew_vacuum` HEURISTIC METASTABLE-side 3.4σ (warning carried, input-conditionality made explicit), `strong_cp` SOUND `|θ̄| < 1.18e-10` (exact rationals), `varying_constants` SOUND Dirac-LNH REJECTED (~2.9 OOM witness), `grh_miller`/`tunnell_bsd` SOUND with the epistemic asymmetry preserved (composite/not-congruent UNCONDITIONAL, prime/congruent CONDITIONAL(GRH/BSD)), `muon_g2`/`h0_tension` SOUND NONIDENTIFIABLE(method) honest refusals, and `coverage` SOUND (Uni-column scarcity confirmed, Stage-4 audit PASS).

Run: `python3 b13_cdl/tests/test_b13.py` (schemas + 8 pipelines + certificates + PIR facts, exit 0). Design: [`b13_cdl/docs/B13-CDL-design.md`](b13_cdl/docs/B13-CDL-design.md); survey addendum: [`b13_cdl/docs/survey-v1.1-addendum.md`](b13_cdl/docs/survey-v1.1-addendum.md).

**Certificate identity + PIR bridge** ([`b13_cdl/docs/pir-bridge-v0.1.md`](b13_cdl/docs/pir-bridge-v0.1.md)): certificates are **content-addressed** (`certificate_id` = hash of the canonical body) with **stable filenames** that reruns overwrite reproducibly; `b13_cdl/tests/test_b13.py` runs **under the CI gate** (`ci/run_all_certified.py` diffs its certificates alongside B1–B12); and each certificate is emitted as a **`pir.Fact`** — modelling `CONDITIONAL(X)` as PIR **assumption-taint** (not a new verdict), so invalidating e.g. `asm:GRH` downgrades exactly the GRH-conditional facts.

**DRAFT items awaiting sign-off — deliberately NOT filed to `docs/claims-table.md`, SPEC, or the atlas:** (1) a `CONDITIONAL(X)` verdict-vocabulary extension to SPEC §6, and (2) the mathematics-entry relabeling B\* → M\* (benchmark-ID collision fix). These live only inside `b13_cdl/` as drafts.

## Research library

[`docs/references.md`](docs/references.md) — ledger through **R33**. **R16–R22** = 2025 Nobel / Berkeley (B9). **R23–R24** = BEC / interferometry. **R25–R30** = related-programs imports. **R31–R33** = α-row tension (S3-EM). **R15** + addendum: B7/B8/engine. **R14** Maudlin: M-layer + arrival-time track.

## Roadmap (proposed next steps)

Successor document: [`docs/roadmap-v2.md`](docs/roadmap-v2.md). Immediate:

1. **S5** — B7 Kelvin-relation real data.
2. **EXP-B** — FDT audit (next puzzle-chain step after S3-EM).
3. **B11** — composition classifier (expands B8 for the lab).
4. **M2-b** — wider grammars; **M1-b** / **M7** measurement layer.
5. Freeze tag `v0.7-frozen` to activate prereg-001.

## Discipline

Every `?` is a **prediction with a falsifier**, never evidence for the synthesis that produced it. Promotion requires: all hard constraints certified · ≥2 domains from one 𝕽 · d_identifiable strictly reduced · survival on held-out **domains** · ≥1 novel blind-tested prediction.
