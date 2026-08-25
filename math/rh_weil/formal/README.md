# AtlasRH — the formal finite theorem boundary

A Lean 4 / Mathlib project holding the finite theorems the Python runtime relies
on, and the machinery that keeps what is proved and what is claimed from drifting
apart.

**No RH proof claim is made here or anywhere in this program.** Every theorem in
this project is finite linear algebra over `ℝ`. None of them is about the
Riemann zeta function, and proving all of them says nothing whatever about the
Riemann Hypothesis. What they do is make one step honest: the step from "Arb
certified these bounds" to "therefore this block is positive definite".

## Pins

| | |
|---|---|
| Lean | `leanprover/lean4:v4.34.0-rc2` (`lean-toolchain`) |
| Mathlib | `f1c1e67f08f57b6d7088b1a98fdceab6da4407ee` (`lakefile.toml`) |
| axioms | `propext`, `Classical.choice`, `Quot.sound` — and nothing else |
| `sorry` | none, anywhere |

The Mathlib revision is an exact commit, not a branch or a range. A floating
dependency would let the meaning of a proved theorem change with no commit in
this repository, which is precisely the drift the rest of this project exists to
prevent.

## Layers

```
AtlasRH/Definitions.lean    Layer A — canonical definitions, no assertions
AtlasRH/Statements.lean     Layer B — the propositions Atlas claims to rely on
AtlasRH/Positivity.lean     Layer C — congruence and the 2x2 criterion
AtlasRH/MatrixInertia.lean  Layer C — Sylvester's law and the 3x3 criterion
AtlasRH/WeilBasis.lean      Layer C — midpoint parity and determinant identities
AtlasRH/CertificateSemantics.lean
                            Layer C — what a rigorous certificate implies
AtlasRH/RankTrace.lean      Layer C — rank-trace, Q = 0 case
comparator/TrustedStatements.lean
                            Layer B — the comparator's view of Layer B
comparator/Solution.lean    Layer C — every theorem, typed against its statement
comparator/PrintAxioms.lean Layer D — the comparison and the axiom report
manifests/theorem_manifest.json
                            the hashed record the Python gate checks
```

`AtlasRH/Statements.lean` imports `AtlasRH.Definitions` and Mathlib and nothing
else — no module containing a proof. A statement written in an implementation
lemma's own vocabulary can be weakened by editing that vocabulary, so the import
list is enforced by `scripts/check_formal_manifest.py`. Each statement is also
spelled out rather than abbreviated: the text is longer, and a reader checking
whether Atlas proved what Atlas says it proved should not have to chase
definitions through a proof library.

## What is proved

Eighteen theorems, each named in the manifest:

* `inertia_congruence_positive`, `inertia_congruence_negative`,
  `inertia_congruence_rank` — Sylvester's law of inertia in the index
  formulation the runtime uses: congruence by an invertible `S` preserves the
  positive index, the negative index and the rank. This is what makes the
  signature the LDL elimination reports a property of the matrix rather than of
  the elimination order.
* `pd_two_by_two`, `pd_three_by_three` — the leading-principal-minor criteria in
  dimensions two and three, both directions. The 3×3 case rests on an exact,
  division-free completed-square identity whose three coefficients are exactly
  the three leading minors.
* `certificate_even2_implies_pd` — honest enclosures plus a positive leading
  entry plus a positive determinant bound at the worst corner implies positive
  definiteness. This is the theorem that gives a rigorous interval run its
  meaning.
* `weil_basis_parity`, `odd_degree3_cross_block`, `odd_degree3_factorization` —
  midpoint parity of `1, b, q1, b3`; cross-parity blocks of a
  reflection-invariant pairing vanish; the parity determinant factorization.
* `pd_three_by_three_certificate`, `diagonal_congruence_preserves_pd`,
  `preconditioned_certificate3`, `diagonal_congruence_preserves_index`,
  `diagonal_congruence_preserves_rank` — the ENG-008 3×3 certificate semantics:
  minor bounds imply definiteness, and the exact dyadic preconditioner changes
  neither the answer nor the signature nor the rank.
* `generalized_rayleigh`, `generalized_pencil_congruence`,
  `preconditioned_gap_certificate3` — the ENG-009 generalized-gap implication:
  shifted positivity of `G − λM` is the Rayleigh bound `λ·vᵀMv ≤ vᵀGv`; the
  certified gap transports across any invertible simultaneous congruence with
  the same `λ`; and the preconditioned 3×3 minor bounds compose all the way to
  the Rayleigh bound for the original pencil.
* `rank_trace_hs` — the rank–trace inequality in the `Q = 0` case, which is the
  case the degree-3 certificate actually uses.

## What is not proved

The general rank–trace inequality, with `Q ≠ 0` and a positive-index bound `b`,
is carried in `AtlasRH/RankTrace.lean` as

```lean
def RankTraceGeneralStatement : Prop := ...
```

a `def` returning `Prop`, with **no inhabitant anywhere in this project**. It is
recorded in the manifest as `EXTERNAL_THEOREM_PENDING_FORMAL_PROOF` with a null
warrant. Carrying it this way is deliberate: it gives a future proof a fixed
target, and it cannot be mistaken for something proved.

The E1 runtime rank–trace result remains valid under its own warrant. Lean does
not retroactively strengthen it.

## The comparator

`comparator/Solution.lean` declares each theorem at the type of a trusted
statement, so Lean's elaborator performs the comparison: a proof of a drifted
proposition does not typecheck under the name it claims.
`comparator/PrintAxioms.lean` re-derives that comparison with `isDefEq` and emits
one tab-separated record per theorem — id, solution theorem, trusted statement,
axiom set, and the elaborated statement pretty-printed by Lean.

`scripts/check_formal_manifest.py` consumes those records and gates two layers:

* **offline**, with no Lean at all: source hashes of every file that defines what
  Atlas claims, no `sorry` under `AtlasRH/` or in `Solution.lean`, no
  project-local `axiom`, the manifest's own content hash, and the pinned
  toolchain and Mathlib commit matching the lakefile. This alone is the drift
  gate — a statement cannot change without changing `Statements.lean`, whose hash
  is recorded;
* **with Lean**: `lake build`, the comparator run, statement-hash equality, and
  an axiom set contained in the three enumerated standard axioms.

Both were exercised against deliberate regressions before being trusted:
weakening `0 < det` to `0 ≤ det` in a statement fails the offline hash check and
fails `Solution.lean` elaboration; an appended `sorry` is caught by the source
scan.

## Commands

```bash
lake build                                    # AtlasRH + comparator
lake env lean comparator/PrintAxioms.lean     # comparator + axiom report
python3 ../scripts/check_formal_manifest.py            # gate (Lean layer optional)
python3 ../scripts/check_formal_manifest.py --with-lean
python3 ../scripts/check_formal_manifest.py --write    # regenerate the manifest
python3 ../scripts/ci_formal.py                        # the rh-formal gate
```

If `lake` is not on `PATH`, set `ATLAS_LEAN_BIN` to a Lean toolchain's `bin`
directory. Without it the manifest gate still runs its offline layer and says so.

## The boundary this project must not cross

A formal theorem strengthens an **exact theorem dependency**. It never converts
interval numerical evidence to FORMAL:

```
Arb interval enclosure                E1     the bounds hold
Lean: positive bounds => PD           FORMAL the bounds suffice
```

Those are separate fields on every PIR fact — `numeric_warrant` and
`logical_implication_warrant` — and `src/formal_evidence.py` is written so the
formal side cannot see, let alone raise, the numeric one. See
[`../external/zeta23/MAPPING.md`](../external/zeta23/MAPPING.md) for the
architecture this borrows from and what it deliberately does differently.
