# `zeta-23-lean` — theorem mapping


**No RH proof claim is made** by Atlas here or anywhere. Every row below is
reference material; none of it warrants an Atlas result.

Which upstream declarations correspond to which Atlas objects, at upstream commit
`cec57f919ccf34e5fa5372b4ba332f7c848bbb6e`. See `PROVENANCE.md` for the pin, the
license chain, and the toolchain incompatibility that keeps every row below at
**reference only**.

The three uses ENG-007 §11 allows are separated here, because they carry very
different weight: (1) architecture, which Atlas has adopted; (2) theorem mapping,
which is what this file is; (3) a future formal dependency, which does not exist.

---

## 1. Architecture — adopted

| upstream | Atlas | note |
|---|---|---|
| `comparator/ChallengeDeps.lean` — definitions from Mathlib alone | `AtlasRH/Definitions.lean` | Layer A |
| `comparator/Challenge.lean` — statements, proofs `sorry`, Mathlib-only imports | `AtlasRH/Statements.lean` + `comparator/TrustedStatements.lean` | Layer B. Atlas differs deliberately: the statements are `def … : Prop` with no placeholder proofs, so **no `sorry` exists anywhere in the Atlas project**, and drift is caught by hashing the statement source rather than by freezing a second copy. |
| `comparator/Solution.lean` — same statements, proved by delegating to the library | `comparator/Solution.lean` | Layer C |
| `comparator/PrintAxioms.lean` — `#print axioms` per statement | `comparator/PrintAxioms.lean` | Layer D. Atlas emits machine-readable tab-separated records and additionally re-derives the `isDefEq` comparison in-process. |
| `AUDIT.md` — reproducible record: sorry count, axiom count, verbatim axiom lines | `manifests/theorem_manifest.json` + `scripts/check_formal_manifest.py` | Atlas makes the audit a gate rather than a document. |
| axiom allow-list `{propext, Classical.choice, Quot.sound}` | `allowed_axioms` in the manifest | identical set, enumerated in both |

Atlas does **not** use the [`leanprover/comparator`](https://github.com/leanprover/comparator)
tool itself. Upstream's Layer D is an external sandboxed replay with an
independent kernel; Atlas's is the Lean elaborator plus a source-hash gate. The
Atlas version is weaker — it does not defend against a compromised toolchain —
and that is stated here rather than glossed.

---

## 2. Theorem mapping — reference only

### Inertia and the positive index

| Atlas | upstream | relationship |
|---|---|---|
| `AtlasRH.PosIndexAtLeast` (`MatrixInertia.lean`) — ∃ a `k`-dimensional subspace on which the form is positive | `Zeta23.RHLinalg.PosDefOn` + `posIndex` (`LinAlg/Sylvester.lean`, `LinAlg/PosIndex.lean`) | same idea, different shape. Upstream defines `posIndex` by counting positive eigenvalues and *proves* the subspace characterization; Atlas takes the subspace form as the definition and never mentions eigenvalues, because its runtime has no eigenvalue solver on any rigorous path. |
| — | `posIndex_eq_max_finrank_posDefOn` — Sylvester's law of inertia, subspace characterization | the theorem that would let Atlas's index-formulation and an eigenvalue-count formulation be identified. Atlas does not need it, and does not have it. |
| `AtlasRH.posIndexAtLeast_congruence_iff` — congruence by an invertible `S` preserves the positive index | `posIndex_conj_le` (`LinAlg/Inertia.lean`) — `n₊(Bᴴ Q B) ≤ n₊(Q)` for arbitrary `B` | upstream's is the more general one-sided bound; applying it twice at `B = S` and `B = S⁻¹` gives Atlas's equality. Atlas proves the equality directly instead, over `ℝ` only. |
| `AtlasRH.negIndexAtLeast_congruence_iff` | `negIndex` + the same `posIndex_conj_le` applied to `-Q` | same route |
| `AtlasRH.rank_congruence` | `posIndex_eq_rank_of_posSemidef`, `rank_hermPosPart` | different decomposition of the same fact |

### Rank–trace

| Atlas | upstream | relationship |
|---|---|---|
| `AtlasRH.rank_trace_zero_Q` (`RankTrace.lean`) — proved: `2 tr P − ‖P‖²_HS ≤ rank P` for a spectrum in `[0,1]` | — | the `Q = 0`, `b = 0` case. This is the case the ENG-006 degree-3 certificate actually uses. |
| `AtlasRH.RankTraceGeneralStatement` — **recorded, not proved**, carried as a `def … : Prop` with no inhabitant | `Zeta23.RHLinalg.rank_trace_ineq_two` (`LinAlg/RankTrace.lean:260`) — **proved upstream**: `2 tr P + 4 tr Q − 4b − ‖P+Q‖²_F ≤ r` under `P ⪰ 0`, `rank P ≤ r`, `Q` Hermitian, `n₊(Q) ≤ b` | **the same theorem.** Upstream has it; Atlas does not. This is the single highest-value row in this table: if Atlas ever wants the general case formally, the target exists, is Apache-2.0, and its hypotheses match the ones `ranktrace/theorem.py` already records. |
| — | `rank_trace_ineq` (`LinAlg/RankTrace.lean:163`) — the general `c > 0` form | the `c = 2` specialization above is the form Atlas's runtime implements |

The `math/rh_weil/ranktrace/` engine currently implements the general inequality
in Python under `THEOREM_ID = "rank_trace_hs_v1"`, with its hypotheses checked at
call time and its warrant recorded as E1 numeric. That warrant is unaffected by
this row: an unimported upstream theorem strengthens nothing.

### Weil explicit formula

| Atlas | upstream | relationship |
|---|---|---|
| the finite Weil compression in `src/finite_weil.py`, `src/weil_entries.py` — Gram entries assembled by rigorous quadrature | `Zeta23.ExplicitFormulaPaper` (`Hypotheses.lean`), `explicitFormulaPaper_of_lit` (`ExplicitFormula/Bridge.lean:191`), `EF_lit` (`ExplicitFormula.lean:81`) | **not a mapping yet.** Upstream formalizes the explicit formula in the analytic form the Alpöge–Furman paper needs, conditional on `EF_lit` and Stirling-type `GammaFacts`. Atlas's compression is a different object: a finite Gram matrix in a fixed basis over a fixed interval, with the archimedean and prime contributions evaluated as intervals. Relating the two is a research task, not a lookup, and ENG-007 §8 explicitly declines to formalize the digamma integral. |

---

## 3. Future formal dependency — does not exist

No Atlas build fetches, imports, or links anything from upstream. Making one a
dependency requires, in order:

1. proving toolchain compatibility (upstream Lean v4.33.0-rc2 / Mathlib
   `51e6992…`; Atlas v4.34.0-rc2 / Mathlib `f1c1e67…` — see `PROVENANCE.md`);
2. building `Zeta23` and recording the **elaborated** statement of the imported
   theorem, not the source-text hash `theorem_manifest.json` currently carries;
3. an axiom report for the imported declaration;
4. the provenance record §11 requires, in full.

Until all four exist, the honest reading of this file is: Atlas knows where the
general rank–trace theorem is proved, and has not used it.
