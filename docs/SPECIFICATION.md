# SPECIFICATION.md — Normative definitions (v1.0, 2026-07-14)

This document is the machine-checkable contract of the program. Anything
claimed by a benchmark, sprint, or atlas edit is interpreted against these
definitions. Changes require an edit record.

## 1. The object 𝕽 (working formalization)

𝕽 = (𝒞, ⊗, 𝟙, †, Ω, 𝖬, 𝖤, 𝖢, 𝖱) where:

- **𝒞** — typed process category: objects = system types; morphisms =
  admissible processes; sequential composition ∘.
- **⊗** — parallel composition (monoidal); **𝟙** — trivial system.
- **†** — adjoint/reversal structure on morphisms.
- **Ω** — admissible state family: for each type A, a convex set Ω_A of
  states with positive evaluation against effects.
- **𝖬** — measurement interface: maps (state, apparatus model) →
  distributions over records. *No claim about 𝕽 is evaluable except
  through some declared 𝖬.*
- **𝖤** — realization functor: latent structure → observable statistics.
- **𝖢** — hard-constraint family (the restraint columns): Sym, Pos, Uni,
  Cau, Cmp, RG, Top, Thm — each a predicate on (𝒞, Ω, 𝖤) with a
  certificate format.
- **𝖱** — scale/coarse-graining evolution: a semigroup
  𝖱_{μ₂←μ₁}: 𝕽_{μ₁} → 𝕽_{μ₂} with fixed-point and monotonicity data.

This tuple is a *working* formalization: its own adequacy is a falsifiable
claim (§6). The constants of the atlas are conjectured to be invariants of
𝕽; every atlas cell is a (row = constant block, column = 𝖢-member) claim
with a certificate.

## 2. Mandatory three-layer separation

Every artifact must declare which layer it lives in:

1. **UNIVERSAL STRUCTURE** — properties claimed for all domains
   (positivity, composition, state evaluation, coarse-graining,
   invariance, consistency).
2. **DOMAIN REPRESENTATION** — how a field realizes the structure (Choi
   matrix, moment matrix, S-matrix, transport matrix, modular operator,
   horizon). Success of an *encoding* is evidence at this layer ONLY.
3. **MEASUREMENT INTERFACE** — the declared 𝖬: apparatus transfer
   functions, preprocessing provenance, calibration route, noise family,
   selection/censoring. Fields `calibration_route` and
   `m_layer_stipulations` are mandatory in every certificate (≥ v0.3).
   This layer is executable as **M7** (`measurement_interface/`,
   `docs/M7-measurement-interface-scope.md`): apparatus transfer + injection
   recovery, independent-vs-self-consistent calibration routes, the
   `APPARATUS_LIMITED` gate (§4), and censoring correction.

An encoding success at layer 2 is never, by itself, evidence for layer 1.
(This codifies the standing no-circularity rule.)

## 3. Evidence levels (every certificate and atlas cell carries one)

- **E0** — exact theorem / exact arithmetic (Fraction / Gaussian-rational
  Schur pivots; symbolic identities).
- **E1** — interval-certified: rigorous input intervals propagated to a
  certified output interval.
- **E2** — statistical: explicit likelihood with stated coverage
  (e.g. Clopper–Pearson in B4; χ² combination in S3).
- **E3** — simulation-conditioned: verdict depends on a numerical model
  or generator (most BENCHMARK-tier results).
- **E4** — proxy: the measured object is an indirect proxy for the latent
  claim (all ANALOGUE-tier physics; the S4 free-field entropy).

Tier labels (BENCHMARK / ANALOGUE-REAL / FUNDAMENTAL-REAL) compose with
evidence levels; only FUNDAMENTAL-REAL at E0–E2 can change an atlas cell.

## 4. Outcome vocabulary (extended)

Cell/claim outcomes: **H, P, ?, —** (atlas) and **FORCED, PERMITTED,
REJECTED** (completions), plus the following first-class verdicts, which
the engine and benchmarks must emit *with a stated cause*:

- **NONIDENTIFIABLE(cause)** — the data cannot single out a value/structure
  (cause ∈ {insufficient intervention, gauge freedom, model order, ...}).
- **OBSERVATIONALLY_EQUIVALENT(class)** — output is an equivalence class
  of structures indistinguishable under the declared 𝖬 (see §5).
- **APPARATUS_LIMITED** — the verdict would change within the apparatus's
  declared resolution/bandwidth; no structural conclusion permitted
  (e.g. rank plateaus at detector cutoffs).
- **REPRESENTATION_DEPENDENT** — the quantity is chart-dependent
  (R15 canonicalization rule); only invariants may be promoted.
- **AMBIGUOUS** — classification tie (B8 semantics), distinct from
  NONIDENTIFIABLE: the menu contains ≥2 compatible grammars.

*Conditional results* (a verdict that holds only under a stipulated unknown,
e.g. GRH/BSD/lattice-coefficient) are **not** a separate verdict: they are
recorded as **assumption-taint** on the result (`asm:X`), so withdrawing X
downgrades exactly the dependent results via the PIR invalidation traversal.
This is the provisional resolution of the proposed `CONDITIONAL(X)` extension —
see `docs/adr/ADR-0002-conditionality-as-taint.md` (Status: PROVISIONAL).

## 5. Observational equivalence and the "same" relation

𝕽₁ ~_obs 𝕽₂ ⇔ for all declared measurement interfaces 𝖬 in scope and all
admissible interventions I: P_{𝕽₁}(D | I, 𝖬) = P_{𝕽₂}(D | I, 𝖬).

A decompilation output is an equivalence class [𝕽]_{~obs}, never a favored
representative. Claims of "the same structure across domains" must state
the equivalence used: equality / isomorphism / Morita-type equivalence /
gauge-orbit / feature-level (canonical-invariant) equality. The default
claim level for cross-domain results is **feature-level**: equality of
canonical invariants (positivity-cone data, rank sequences, symmetry
group, factorization poset, flow exponents) — nothing stronger, unless an
intervention-based test upgrades it.

## 6. Promotion, falsification, and the claims table

Promotion to H additionally requires (unchanged rules + new):
certified constraints at stated evidence level; ≥2 domains; reduced
d_identifiable; survival on held-out DOMAINS; ≥1 novel prediction; AND the
three-layer declaration of §2 with tier/evidence labels.

**Claims table (Priority-0 discipline):** every result is filed as exactly
one of {known theorem | recovered known result | empirical consistency |
framework prediction | novel prospective prediction | speculative
interpretation}. README badges must not upgrade a category.

**Global falsifiers of the program itself:** (i) a rival grammar family
(§ roadmap, B12-RGRC) matches all benchmark verdicts with fewer
primitives; (ii) the frozen-rules prospective prediction protocol fails at
registered thresholds; (iii) canonical invariants of two presentations of
the same physics disagree (canonicalization failure). Any of these is
filed as a program-level REJECTED, not silently repaired.
