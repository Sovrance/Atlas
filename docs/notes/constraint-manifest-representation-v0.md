# Constraint-manifest representations — design note v0 (B15-SURF WP4)

**Status:** design note, no code. Output of work-order WP4 (report pattern P2: prefer
representations in which hard constraints are *structural*, not *checked*). Ends with §9-OI-4.
**Scope:** B10 (CV Gaussian channels, `b10_cv_channel/`) as the primary sector, B2/B12-b
(qubit CPTP, `b2_process_solver/`, `b12_rgrc/quantum_*`) as the secondary one.

## 1. What is *checked* today

| Sector | Constraint column | Current representation | How it is enforced | Pass tag |
|---|---|---|---|---|
| B10 | **Pos** (complete positivity) | `(T, N)` real rational matrices | `M_CP = N + (i/2)(Ω − TΩTᵗ) ⪰ 0` **checked** by `SCHUR_PIVOT_EXACT` on a Gaussian-rational Hermitian matrix (`certify_channel`) | SOUND (E0) |
| B10 | **Uni** (unitarity / symplecticity of the noiseless part) | same | `TΩTᵗ = Ω` **checked** as exact polynomial relations; `det T = 1` witness; hidden entries FORCED by solving those relations (`force_symplectic_entry`) | SOUND (E0) |
| B10 | **Cau** | — | not represented: Gaussian channels are static maps; causality (no signalling / input–output ordering) is implicit in the `V → TVTᵗ + N` form and is never a pass | — |
| B2 / B12-b | **Pos** (CP) | Choi matrix `J` (Gaussian-rational 4×4) | `J ⪰ 0` **checked** (`CPTP_GATE` = Choi + `SCHUR_PIVOT_EXACT`); B12-b's `classify_process` bootstraps `min eig` on reconstructed Choi (HEURISTIC, E3) | SOUND / HEURISTIC |
| B2 / B12-b | **Uni** (trace preservation) | same | `Tr_out J = 𝟙` **checked** by `is_trace_preserving` / `block_traces` | SOUND |
| B2 / B12-b | **Cau** | Choi with input/output tensor order | implicit (a Choi matrix *is* a map with a fixed direction); the "no-signalling" analogue would be the block-trace identity, which is checked as TP | — |

Summary: **Pos and Uni are checked predicates in both sectors; Cau is structural by
construction but undeclared.** Every PERMITTED/FORCED verdict in B10 and B2 therefore rests on
a positivity *test* of a matrix whose entries were free to violate it — exactly the "checked,
not manifest" situation P2 warns about, and exactly the pattern that the surfaceology
`u`-variables replace (locality/factorization become identities `u_C + Π u_D^{n(C,D)} = 1`,
R36).

## 2. Candidate coordinates in which a constraint holds identically

### 2.1 B10 — Gaussian channels

| Chart | Parameters | What becomes automatic | What it hides / costs |
|---|---|---|---|
| **(a) Stinespring / dilation chart** — `T = P S`, `N = P S₀ Sᵗ P^ᵗ`-type: any Gaussian channel is a symplectic `S` on system ⊕ environment (`2(n+k)` modes) restricted to the system, with the environment in a Gaussian state of covariance `V_E ⪰ (i/2)Ω_E` | `S ∈ Sp(2(n+k), ℚ)` (rational symplectic via Cayley/exponential-free parametrisations, e.g. products of rational shears and rotations) + `V_E` | **Pos becomes an identity**: `M_CP ⪰ 0` holds for *every* parameter value, no pivot test needed; **Uni is a sub-chart** (`k = 0`) | environment dimension `k` is a model-order choice (NONIDENTIFIABLE(model order) risk); the same `(T,N)` has many dilations (gauge orbit under `Sp(2k)` on the environment) — the R15 quotient is mandatory; rational symplectic parametrisations of *all* of `Sp` are not surjective over ℚ (Cayley transform misses `−1` eigenvalues) |
| **(b) Williamson normal form** — `N`-side: `M_CP = D (⊕ νⱼ 𝟙₂) Dᵗ`, symplectic `D` and symplectic eigenvalues `νⱼ ≥ ½` | `D ∈ Sp`, `νⱼ ∈ ℚ_{≥1/2}` | Pos is manifest as `νⱼ ≥ ½` (a *linear* box, not a PSD cone); the quantum-limited boundary (B10's rank-deficient FORCED case, Caves bound) is the face `νⱼ = ½` — an **extremal-ray coordinate**, exactly R26's "FORCED boundaries as extremal rays" | `D` and `ν` are only rational for special channels (symplectic diagonalisation needs square roots in general) — E0 arithmetic is not guaranteed; hides the `T`-dependence |
| **(c) Symplectic normal form of `T` alone** (Cooper/Holevo canonical forms: `T = S₁ T_c S₂`, `T_c ∈ {𝟙, τ𝟙, diag(1,0), …}`) | `S₁, S₂ ∈ Sp`, class label + one gain `τ` | Uni manifest (`T_c = 𝟙`); channel *class* (attenuator/amplifier/…) is a discrete invariant → M1 canonical feature | Pos still a check on `N`; the normal-form class is itself a `REPRESENTATION_DEPENDENT`-free invariant but the decomposition is not unique |

### 2.2 B2 / B12-b — qubit CPTP

| Chart | Parameters | Automatic | Hides / costs |
|---|---|---|---|
| **Kraus / Stinespring** `J = Σ_k vec(K_k) vec(K_k)†`, `Σ K_k†K_k = 𝟙` | `K_k` (Gaussian-rational) | **CP is an identity** (`J` is a Gram matrix); TP becomes an *equality* constraint on the `K_k` | Kraus rank is a model-order choice; unitary freedom `K_k → Σ u_{kl} K_l` (R15 quotient); TP is still checked unless one uses an isometry parametrisation |
| **Isometry chart** `V: ℂ² → ℂ²⊗ℂ^k`, `V†V = 𝟙` | Stiefel-manifold coordinates | CP *and* TP identities | rational points on the Stiefel manifold are sparse → E0 arithmetic only on a sub-family; hides nothing physical but makes exact benchmarks harder to seed |
| **Choi (current)** | `J` | Cau/direction manifest | Pos, Uni checked |

## 3. What each makes automatic vs. what it hides — the honest ledger

- Making Pos structural (charts a, Kraus) trades a **PSD test** for a **gauge quotient +
  model-order choice**. Nothing is free: the checked inequality becomes an unobservable
  redundancy (`Sp(2k)` on the environment; `U(k)` on Kraus indices) — which is precisely the
  "representation redundancy signals wrong representation" warning (P1) turned around: the
  *right* redundancy is the one that only re-parametrises the environment we already declared
  we cannot see (M-layer).
- Making Uni structural (`k = 0`, isometries) removes the *noise* sector; B10's most valuable
  result (Caves amplifier noise FORCED) lives exactly where Uni fails and Pos saturates — so a
  Uni-manifest chart hides the extremal-ray physics unless combined with (b).
- Cau: in both sectors it is already structural but **undeclared**. The cheapest concrete
  gain is to *declare* it as a chart property in the PIR fact (`constraint_manifest: [Cau]`),
  so that a future non-Markovian or bidirectional representation (where Cau would have to
  be checked) is visibly a different chart.
- Exactness: only chart (a) with rational shear/rotation generators and the Kraus chart keep
  the E0 class for generic parameters; (b)/(c)/isometries generally leave ℚ.

## 4. Pre-registerable Stage-3 threshold (proposal for prereg-003)

**Claim to register:** "≥ 1 previously-checked constraint becomes an identity, verified on
B10 + B15-ZERO data."

- **B10 leg.** Re-express every B10 channel in the test set through the dilation chart (a)
  with `k = 1` environment mode and rational symplectic generators; assert (i) round-trip
  `(T,N) → (S, V_E) → (T,N)` exact, (ii) `M_CP ⪰ 0` is *never* evaluated as a pass on the
  chart-side and holds identically (proof-level: `M_CP` is a principal compression of a PSD
  matrix), (iii) the FORCED Caves boundary reappears as the face `V_E = (i/2)Ω_E`-saturated
  (`ν = ½`). PASS = all three exact on the existing B10 certificate cases + the same
  negative controls (a non-CP `(T,N)` must have **no** rational preimage — the chart's
  "REJECTED" is an *inconsistent linear system*, i.e. a genuine `verify_farkas` certificate,
  which is the equality-form dual B15-POS could only partially use).
- **B15-ZERO leg.** The amplitude analogue: `u`-variables (R36) make factorization/locality
  identities. Registered sub-claim: on the B15-ZERO split grids, the rank-1 factorization
  (checked today by `RANK_TEST`) becomes an identity of the chart — i.e. in
  `(x, y)`-coordinates of the sub-surfaces (2405.09608 eq. 2.10) the split is manifest and
  the *checked* object is instead the linear relation between `X` and `(x, y)`. PASS = the
  grid matrix is rank 1 *by construction* of the coordinates and the reverse map recovers the
  registered vanishing mesh set (`c_{1,3}, c_{3,5}` at n=6) exactly.
- **Threshold:** both legs PASS with zero degradations of existing certificates; a FAIL on
  either is filed as REJECTED for that chart, not repaired. Evidence level of a PASS: E0,
  category *known theorem (re-expression)* — the same category as B6/B10 today; **no atlas
  cell moves** on this (a representation change cannot promote H).

## 5. §9-OI-4 — for Erick

Open B16 on the constraint-manifest representation? If yes, the recommended scope is the two
legs of §4 as a `B16-MANI` sprint (stdlib, E0, `prereg-003`), with chart (a) for B10 and the
`(x,y)` sub-surface chart for B15-ZERO, and an explicit `constraint_manifest: [...]` field
added to the PIR fact schema so that "checked" vs "identity" is machine-readable per fact
(this is the same honesty axis as `ground_truth_route`, §9-OI-6). If no, this note stays as
the recorded design and the `constraint_manifest` field goes to the backlog.
