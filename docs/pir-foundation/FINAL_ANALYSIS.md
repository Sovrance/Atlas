# Reverse Engineering → Decompilation of Physics
## Final Complete Analysis (Post Cross-Agent Adjudication)

**Project:** Global-variables (github.com/ErickLGonzalez/Global-variables)
**Inputs merged:** Rosy 🌹 (Claude) deep-research report + ChatGPT "Reverse-Engineering Architecture" package
**Adjudication:** 20 patterns ACCEPTED, 9 ACCEPT-MODIFIED, 2 REJECTED/RECONSOLIDATED, 2 DEFERRED (full table in ADJUDICATION.md)
**Date:** 2026-07-15 · **Status:** Consolidated design basis for the PIR workstream · **Portable to:** MnemesisOS (memory/state/provenance sections flagged [→M])

---

## 0. One-paragraph thesis

The Global-variables program is already, structurally, a decompiler: the restraint matrix, evidence levels E0–E4, certified pipelines, and NONIDENTIFIABLE / OBSERVATIONALLY_EQUIVALENT verdicts are direct analogues of Ghidra's SLEIGH specs, P-code IR, sound-vs-heuristic analysis passes, and undecidability handling. The transferable object from software reverse engineering is **not** the metaphor that nature is a compiled program — it is the mature **engineering discipline for reconstructing latent structure from lossy observations**: preserve raw evidence immutably; formalize loaders and measurement semantics declaratively; lift heterogeneous observations into a small typed intermediate representation; run modular append-only analyses; retain ambiguity with provenance as a first-class output; combine passive data with declared interventions; validate reconstructions by forward-recompiling them into predictions; and define success as an intervention-tested equivalence class, never a unique source.

## 1. The Physical Intermediate Representation (PIR) — merged design

**Pattern source:** Ghidra P-code (small orthogonal op set over typeless varnodes in explicit address spaces, SSA form) + Binary Ninja BNIL (stacked IL levels) + LLVM (typed IR, verifiable passes).

**PIR has four representation levels (L-axis):**
- **L0 — raw records.** Samples, calibration lineage, apparatus manifests, source spans. Content-addressed, immutable, no theoretical labels.
- **L1 — operational acts.** What was done: PREPARE, INTERVENE, EVOLVE, COUPLE, SPLIT, RECOMBINE, TRANSPORT, PHASE_ACCUMULATE, MEASURE, RECORD, CONDITION, TRACE_OUT, COARSE_GRAIN, SYMMETRY_ACT, COMPOSE, PARALLEL. Ops carry typed ports, units, uncertainty, timing, apparatus identity, calibration route, source spans, assumptions.
- **L2 — structural facts.** Inferred causal edges, symmetry candidates, rank sequences, positivity cones, factorization structure, loop relations, scale-flow features.
- **L3 — candidate grammars.** Competing explanatory families (Hamiltonian / gradient-flow / GENERIC / quantum-CPTP …) with compatibility matrices. Every L3 node retains a derivation path to L0.

**Two op families, one IR (adjudicated merge):**
- **Act-ops** (L1): the experimental vocabulary above — what happened.
- **Verifier-ops** (the certified analysis layer): the project's existing de-facto opcode set — exact-Fraction Schur pivot, rank tests, flat-extension recurrence, Jacobian-rank identifiability, CPTP/symplectic gates. These are the ops B1–B12 already share; freezing and documenting them per-op (P-code-reference style) is Stage 1 work, not new invention.

**Critical orthogonality rule:** PIR **L-levels** are *representation abstraction*; evidence **E-levels (E0–E4)** are *warrant strength*. They are independent axes. An L2 structural fact can be E0 (exact certificate) or E4 (proxy). Every fact carries both coordinates. Conflating them is a spec violation.

**Namespaces as physics address spaces [→M]:** `raw:`, `apparatus:`, `operational:`, `domain:`, `latent:`, `gauge:`, `invariant:`, `global:`, `effective:`, `analyst:`. Every cross-namespace transform is named, typed, and recorded. This mechanically prevents detector voltages from being read as fundamental observables, and effective parameters from being silently promoted to universal constants.

**Typed structure:** State, Effect, Observable, Process, Coupling, SymmetryParameter, GaugeCoordinate, Invariant, NuisanceParameter, ApparatusParameter, BoundarySector, GlobalCandidate, ScaleDependentEffectiveParameter. Type uncertainty stays explicit; invalid compositions are detectable at lift time.

## 2. Declarative domain specification (SLEIGH pattern)

Each physics domain gets a **Domain Semantics Module** with four contracts, cleanly separating spec from engine:
1. **Source syntax mapping** — how the raw artifact's columns/fields are recognized (a flux-bias column in a Josephson dataset).
2. **Apparatus semantics** — calibrated transfer functions, resolution/bandwidth limits (declared numerically, so APPARATUS_LIMITED verdicts can name the exact bounding threshold).
3. **PIR lowering rules** — how recognized syntax becomes L1 act-ops (INTERVENE(loop_flux, value)).
4. **Invariant/gauge declarations** — what the domain claims is invariant vs coordinate (closed-loop phase as candidate invariant).

Likewise, each restraint-matrix predicate (Sym, Pos, Uni, Cau, Cmp, RG, Top, Thm) is a **declarative spec object with a fixed certificate schema**, separate from the engine evaluating it. One engine, many auditable, versionable specs — the reason Ghidra supports dozens of processors with one analysis core.

## 3. The fact store: ambiguity, provenance, honesty [→M]

**3.1 Superset/PGSD candidate retention (both agents, independently).** Following Superset Decompilation (Liu et al. 2026, arXiv:2603.28002): the fact store is **append-only and monotonic**; passes derive facts, never mutate or invalidate each other's conclusions; **ambiguous interpretations persist as parallel candidates with provenance**; resolution is a separate, optional, final selection phase. OBSERVATIONALLY_EQUIVALENT *is* a candidate forest; NONIDENTIFIABLE *is* a certified refusal to select. Rival-generator testing (M4/B12-RGRC) is the "selecting a representative" phase — with non-selection as a legal outcome.

**3.2 Assumption-taint + invalidation traversal (ChatGPT's best contribution) [→M].** Every derived fact carries queryable dependencies on its assumptions (Gaussian noise, detector linearity, Markov approximation, gauge choice, truncation order, calibration fit) using PROV-O's Entity–Activity–Agent core. **Revising or retracting an assumption automatically identifies and downgrades every downstream fact and certificate.** This is the mechanized form of the erratum culture and the no-circularity rule.

**3.3 SOUND/HEURISTIC pass taxonomy (Rosy addition).** Every analyzer pass is tagged **SOUND** (exact, E0–E1: changes only what it can certify — Ghidra's auto-analyzer rule) or **HEURISTIC** (E3–E4). Heuristic passes must emit a located `warnings[]` field (Ghidra's `WARNING:` comment pattern) — the boundary between proven and guessed is surfaced in-line, never buried.

**3.4 Four storage layers [→M].** Immutable **Evidence** → machine **Facts** → **Analyst annotations** (attributed, versioned, reversible — never promoted to facts) → **Certified claims**. Analysis is a layered, diffable overlay on immutable ground truth; hard-won annotation survives data updates via correlator-scored porting (Ghidra Version Tracking pattern: when GWTC-next or a revised α_s lands, atlas markup migrates with scored associations and retained history rather than being redone).

**3.5 Cross-agent merge semantics (Rosy addition) [→M].** From Ghidra Server: non-conflicting atlas edits auto-merge; conflicting same-cell claims require explicit resolution under deterministic per-cell rules; **"keep both, with provenance" is a legal resolution** (the candidate forest again); exclusive locks = one agent owns a benchmark during a sprint.

## 4. Analysis machinery

**4.1 Analyzer runtime.** Small passes over the shared store, dependency-resolved, deterministic where possible: UnitNormalizer, MeasurementProvenanceAnalyzer, DimensionalConsistencyAnalyzer, PositivityAnalyzer, SymmetryAnalyzer, CausalGraphAnalyzer, CompositionAnalyzer, GaugeDependenceAnalyzer, IdentifiabilityAnalyzer, TopologyAnalyzer, ScaleFlowAnalyzer, ObservationalEquivalenceAnalyzer, GlobalVariableCandidateAnalyzer. Passes append; conflicts are stored as conflicts. Diagnostic hooks (adjudicated down from "debugger"): assertions on first positivity violation / first rank transition / first apparatus-resolution crossing; watched quantities (min eigenvalue, determinant, rank, residuals, beta functions); every failed claim yields a backtrace from verdict to raw records and assumptions.

**4.2 Symbolic constraint bridge (SAT/UNSAT certification).** Unknown constants, couplings, latent states, noise variables, boundary sectors, apparatus parameters are symbolic; each restraint-matrix column emits constraints; solvers (Z3/SymPy) return feasible regions, forced values, excluded regions, underdetermined dimensions. **Certificate vocabulary:** SAT model = *witness* (inner interval feasible), UNSAT bracket = *impossibility certificate*, ideally with UNSAT core (B1's exact negative Schur pivot is already one). Target: certificates a third party re-verifies from contents alone, without rerunning the pipeline.

**4.3 Fingerprinting & signature library (two-hash design).** Ghidra Function ID pattern: a **full hash** (similarity-invariant canonical signature — R15 already computes these) for robust recognition, plus a **specific hash** (finer invariants) for disambiguation, over a versioned database of known grammars and motifs (phase loops, avoided crossings, Lindblad decay, tunneling, Gaussian channels, diffusion kernels, Schur-complement positivity, Onsager reciprocity, critical scaling). Matches emit RECOVERED_PATTERN facts with confidence and required distinguishing tests — never ontological identity. B8's blind grammar ID becomes a fingerprint lookup. A hash collision between physically distinct grammars is filed as REPRESENTATION_DEPENDENT, not patched.

**4.4 Cross-domain diffing with dual scores.** PIR-Diff / Invariant-Diff across experiments and domains, reporting shared motifs, differing apparatus assumptions, invariant matches, scale mismatches — and per BinDiff discipline, **similarity and confidence as separate numbers with the correlator identity attached**. High-similarity/low-confidence = a lead, not a claim. Each diff report names the least-cost intervention that would decide whether similarity is structural.

**4.5 GVAR rule language (constrained).** Declarative, versioned, fixture-tested YARA-style rules (GlobalHolonomyCandidate, HiddenCommonCauseCandidate, ApparatusAliasingCandidate, NonMarkovianMemoryCandidate, CrossScaleInvariantCandidate). **Rules emit candidate facts and test obligations only — never verdicts.** Verdict vocabulary is locked to SPEC §6; ChatGPT's proposed labels (GLOBAL_CANDIDATE, TOPOLOGICAL_CANDIDATE, HIDDEN_COMMON_CAUSE_CANDIDATE, REPRESENTATION_ARTIFACT, NOT_DETECTED) survive only as **candidate-class tags inside the hypothesis store**. No discovery is ever emitted automatically.

**4.6 Global-candidate bar.** A global-variable candidate is not a highly connected node. It must show: cross-subsystem influence, resistance to local decomposition, gauge-invariant observable consequences, and intervention-supported discrimination from common-cause and apparatus explanations.

## 5. Validation & forward loop

**5.1 Five gates.** Structural validity → mathematical validity → retrodictive validity → held-out validity → prospective validity. Failures produce structured counterexamples that guide the next pass. **Only prospective success under frozen rules can materially promote a universal claim** (the existing promotion rule, mechanized).

**5.2 Forward recompilation / counterfactual patching.** Decompile candidate grammar H; patch one component (remove the global variable; replace with locals; break a symmetry; alter a topological sector) → H′; forward-compile H′ through realization functor + apparatus model into predicted records; compare residuals on observed and held-out data. A useful reconstruction must generate a discriminating prospective prediction. This is RE's round-trip differential testing, run in physics.

**5.3 Intervention-relative observational equivalence.** Static analysis consumes completed datasets; **dynamic analysis formalizes declared interventions** (apparatus swaps, boundary changes, timing permutations, loop deformations, scale sweeps). Observational equivalence is implemented **relative to the declared intervention set**, not only passive distributions — two grammars OBSERVATIONALLY_EQUIVALENT under passive data may separate under a cheap declared intervention, and the engine should say which one.

**5.4 Experiment search (reframed).** Objective-driven intervention search — maximize candidate disagreement, expected information gain, d_identifiable reduction — with feasibility/cost/safety filters. This is **optimal experimental design**, not fuzzing (the AFL++ analogy supplied the retain-what-separates loop; the statistics is OED and should cite that literature). Whole machine tagged HEURISTIC; scheduled late; an intervention is retained when it eliminates a candidate, reaches an untested structural regime, or improves identifiability. Negative control: a case where no admissible intervention identifies the model (must emit NONIDENTIFIABLE).

**5.5 CI-enforced frozen certificates [→M].** One headless entrypoint reruns every certified pipeline on every commit, regenerates certificates, and **fails the build on any certificate degradation or failure-to-certify**. This is what makes `v0.7-frozen` and prereg-001 enforceable rather than aspirational. Reproducibility: pinned dependencies, deterministic seeds, content hashes, signed run manifests.

**5.6 Failure-mode discipline.** From decompiler pathology: readable ≠ correct; compilation is many-to-one; metrics get gamed when elegance is optimized without correctness gates. **Never equate elegance, compression, or cross-domain resemblance with truth.** Track evidence independence — multiple tools/models sharing assumptions are one source, not several (this applies to the two agent reports themselves: their agreement is correlated via shared primary sources).

## 6. Target architecture (adjudicated v0.2)

```
Raw Artifact → Loader → Domain Semantics Module → PIR (L0/L1)
   → Append-only Fact Store (L2; provenance, assumption-taint, namespaces)
   → Analyzer Runtime (SOUND/HEURISTIC passes; diagnostic hooks)
   → Candidate Grammar Store (L3; compatibility matrix; candidate-class tags)
   → { Symbolic Constraint Bridge · Fingerprint/Signature Library ·
       PIR-Diff (similarity+confidence) · GVAR Rules · Intervention Search }
   → Forward Recompiler (patch → predict → held-out compare)
   → Certification & Atlas Bridge (SPEC-§6 verdicts only; tier/evidence
     rules enforced; draft atlas edits with full derivation chains;
     cannot mint promotions)
   → Headless CLI + CI (frozen-certificate enforcement)     [UI: deferred]
```
Hard rules: the atlas engine never reads raw datasets directly; raw evidence never contains inferred universal claims; every boundary has a schema and validation contract; a paper parser never becomes a theory interpreter.

## 7. Staged engineering plan (compressed from 10 phases to 3 gated stages)

**Stage 1 — Substrate (P0–P1, ~3 sprints).**
ADR + versioning policy; PIR v0.1 schemas (artifact, event, fact, hypothesis, intervention, provenance) with L×E coordinates, namespaces, assumption-taint, SOUND/HEURISTIC tag; Python models; canonical hash-stable serialization; provenance query + invalidation traversal; **verifier-op reference doc** (freeze the existing Schur-pivot opcode set); **CI headless entrypoint** wrapping existing `tests/test_bN.py`. Acceptance: round-trip serialization; invalid cross-namespace casts rejected; invalidation traversal works; existing tests untouched and passing; CI fails on synthetic certificate degradation.

**Stage 2 — First lowering & analysis (P2–P5, gated on Stage 1, ~6 sprints).**
B9 circuit Domain Semantics Module + loader; reproduce the existing B9 verdict from PIR facts (negative controls preserved; no atlas changes). Analyzer runtime + MeasurementProvenance / ObservationalEquivalence / GlobalVariableCandidate analyzers (deterministic reruns; conflicts retained; analyzer failure cannot corrupt the store). Symbolic bridge (reproduce one B1/B2 forcing case + one B3 identifiability case; witnesses stored). Candidate lattice + GVAR rules (tie → OBSERVATIONALLY_EQUIVALENT; insufficient interventions → NONIDENTIFIABLE). Two-hash fingerprint stub over B8 grammars.

**Stage 3 — Forward loop & cross-domain (P6–P8, gated on Stage 2, ~7 sprints).**
Forward recompiler + counterfactual patching (alter one B9 latent component → changed held-out spectroscopy prediction; no held-out reuse during reconstruction). Intervention search with OED objectives + NONIDENTIFIABLE negative control. Second domain: BEC interferometry spec + loader; PIR-Diff vs B9 with dual scores and a proposed discriminator. **Deferred indefinitely:** workbench UI. **Downgraded:** security threat model → single rule (sandbox untrusted parsers; no code execution from papers/data).

**Plan-changing thresholds:** if the monotonic store is too heavy at current scale, fall back to SOUND/HEURISTIC tagging + dual-score diffs as the minimum honesty layer; fingerprint collisions between distinct grammars are filed as findings, not patched; certificate instability across CI runs freezes numerics (seeds, solver versions) before any promotion.

## 8. MnemesisOS transfer notes [→M]

The patterns marked [→M] form a coherent spec for a provenance-aware memory system, independent of physics:

1. **Immutable substrate + layered scored annotation + retained history** (§3.4): memories are an overlay on append-only evidence; interpretation never overwrites experience.
2. **SSA discipline** (§1): every write creates a new value; merges are explicit phi-points → memory states become diffable and analyzable.
3. **PGSD candidate retention** (§3.1): store competing interpretations of state with provenance and defer collapse; ambiguity is a representable memory condition, not an error.
4. **Assumption-taint + invalidation traversal** (§3.2): when a belief/assumption is revised, every dependent memory and conclusion is automatically found and downgraded — mechanized belief revision.
5. **Namespaces + typed casts** (§1): raw-percept vs interpreted vs analyst/self-annotation namespaces with named transforms prevent inference laundering into observation.
6. **Version-Tracking markup porting** (§3.4): when underlying data updates, annotations migrate via correlator-scored associations with history, rather than being lost or blindly copied.
7. **Merge semantics for multi-agent memory** (§3.5): deterministic per-field conflict rules; "keep both + provenance" as a legal resolution.
8. **CI-style frozen-state verification** (§5.5): memory-integrity invariants re-verified on every mutation batch; degradation fails the operation.
9. **Dual similarity/confidence scoring** (§4.4) for memory retrieval/matching: how alike vs how much the matcher is trusted, reported separately.

## 9. Caveats (both agents, merged)

- **The mapping is a design hypothesis**, not a validated result; every "maps onto" must be tested by building it.
- **Analogy risk is the deepest caveat:** decompilation recovers structure that provably existed; physics decompilation assumes a latent 𝕽 whose existence is the conjecture. RE patterns supply *engineering discipline for certified inference under uncertainty* — they cannot supply evidence that 𝕽 exists. The three-layer no-circularity rule is the guard and is never relaxed because an analogy feels tight.
- **Superset Decompilation is a 2026 preprint** (not yet peer-reviewed); its formalism is compelling, its maturity claims unverified.
- **P-code op count (~65) is derived from the manual's table, not an official figure**; BinDiff's largest similarity term (~50%) is flow-graph MD-index, so the equivalence analogy must not over-index on one structural feature.
- **Agent agreement is correlated:** both reports drew on the same primary sources; convergence counts once.
- **Success is redefined, permanently:** recovery of the minimal, compositional, intervention-tested equivalence class of latent structures compatible with declared measurements — disciplined narrowing, explicit nonidentifiability, and decisive-experiment generation. Not nature's unique source code.
