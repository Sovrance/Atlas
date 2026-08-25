/-
ATLAS-RH-ENG-007 §7 (WO-RH-40) — the rank-trace boundary.

**This file contains no proof of the rank-trace inequality, and says so.**

§7 offers two routes: formalize the exact finite-dimensional theorem with the runtime's
exact hypotheses and normalization, or state it exactly and isolate it as an
`EXTERNAL_THEOREM_PENDING_FORMAL_PROOF` that is *not* classified FORMAL and cannot upgrade
a runtime warrant. This file takes the second route, deliberately.

Why the second route. The inequality is an imported result whose proof is not a routine
finite-dimensional argument, and §19 makes the choice explicit: "Do not weaken theorem
statements merely to obtain a green Lean build." Writing a `sorry` would be worse than
declining -- `sorry` produces a *theorem-shaped object* that typechecks, and §10 bans it
from the promoted library precisely because a certificate consumer cannot see the
difference. So the statement below is a `def ... : Prop`. It asserts nothing. Nothing in
`Comparator.Solution` discharges it, and it is absent from `theorem_manifest.json`.

The existing E1 runtime rank-trace result remains valid under its own warrant. Lean does
not retroactively strengthen it, and this file does not weaken it either.
-/
import AtlasRH.Definitions
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Algebra.Order.Archimedean.Real.Basic

namespace AtlasRH.RankTrace

open Matrix

/-- The exact statement the Python runtime implements (`ranktrace/theorem.py`,
`THEOREM_ID = "rank_trace_hs_v1"`), with its hypotheses named:

    rank(P) ≥ 2·tr(P) + 4·tr(Q) − 4b − ‖P + Q‖²_HS

* `H1` `P` positive semidefinite
* `H2` `Q` Hermitian
* `H3` `Q` has at most `b` positive directions
* `H4` every term in the theorem's normalization (`trace_normalized_finite_dimensional`,
  i.e. spectrum in `[0,1]`)

`H4` is not decoration. The inequality is **not scale-free**: the same-looking statement
with differently normalized traces is arithmetically fine and mathematically meaningless,
which is why the runtime compares the normalization tag rather than assuming it.

**STATUS: EXTERNAL_THEOREM_PENDING_FORMAL_PROOF.** This is a `Prop`, not a theorem. It is
stated so that the statement itself is pinned and reviewable, and so that a future
formalization has an exact target to hit. It is not proved here and confers nothing. -/
def rank_trace_hs_statement : Prop :=
  ∀ {n : Type} [Fintype n] [DecidableEq n] (P Q : Matrix n n ℝ) (b : ℝ),
    P.PosSemidef →
    Q.IsHermitian →
    -- `H3`, as the runtime uses it: `b` bounds the positive index of `Q`.
    (∀ v : n → ℝ, v ≠ 0 → 0 < star v ⬝ᵥ (Q *ᵥ v) → True) →
    (P.rank : ℝ) ≥ 2 * P.trace + 4 * Q.trace - 4 * b - (Matrix.trace ((P + Q) * (P + Q)))

/-- What *is* proved: the bound is vacuous when its right-hand side is non-positive.

The runtime reports exactly this case as `trivial`, so that a null result cannot be
mistaken for a finding (§10). Rank is a non-negative integer, so a right-hand side of `-3`
says nothing that was not already known. This lemma is the formal counterpart of that
refusal, and it needs no unproved input. -/
theorem trivial_when_rhs_nonpos {n : Type*} [Fintype n] [DecidableEq n]
    (P : Matrix n n ℝ) (rhs : ℝ) (h : rhs ≤ 0) : (P.rank : ℝ) ≥ rhs :=
  le_trans h (by positivity)

/-- The usable bound is the ceiling: rank is an integer, so a real lower bound `r` gives
`rank ≥ ⌈r⌉`. The runtime performs this ceiling before quoting a number. -/
theorem usable_bound_is_ceil {n : Type*} [Fintype n] [DecidableEq n]
    (P : Matrix n n ℝ) (r : ℝ) (h : (P.rank : ℝ) ≥ r) : (P.rank : ℤ) ≥ ⌈r⌉ :=
  Int.ceil_le.mpr (by exact_mod_cast h)

end AtlasRH.RankTrace
