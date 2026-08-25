/-
ATLAS-RH-ENG-007 §7 — the rank–trace boundary.

The Python runtime implements

    rank P ≥ 2 tr P + 4 tr Q − 4b − ‖P + Q‖²_HS

for `P ⪰ 0`, `Q` Hermitian with at most `b` positive directions, under a
normalization. §7 offers two acceptable outcomes: prove it, or state it exactly
and isolate the unproved part so it cannot be mistaken for formal.

**What is proved here.** The specialisation `Q = 0`, `b = 0`, which is the case
the ENG-006 degree-3 certificate actually invokes -- and the case that fixes the
normalization. The inequality is not scale free (rank is scale invariant while
`tr P` is degree 1 and the HS term degree 2), so it holds only under a
normalization, and equality at a projection identifies that normalization as
*spectrum in `[0,1]`*. Under that hypothesis the `Q = 0` case is a short
argument: each eigenvalue contributes `λ(2 − λ) ≤ 1` to the right-hand side and
each nonzero eigenvalue contributes `1` to the rank.

**What is not proved here.** The general statement with `Q ≠ 0`. It is stated,
its hypotheses are named, and it is marked `EXTERNAL_THEOREM_PENDING_FORMAL_PROOF`
in the manifest. There is no `sorry` in this file: an unproved statement is
carried as a *definition of a proposition*, never as a theorem with a hole. A
consumer asking whether the general form is formal gets "no", not a build error
and not a silent yes.

This does not weaken the ENG-006 runtime result. That result stands on its own
warrant (E1, hypotheses checked at runtime); Lean does not retroactively
strengthen it, and §7 says so explicitly.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import Mathlib.Data.List.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

namespace AtlasRH

/-- The normalization under which the rank–trace inequality is stated: every
eigenvalue of `P` lies in `[0, 1]`. Equality at a projection is what pins this
down; without it the inequality is false by scaling. -/
def SpectrumInUnitInterval (l : List ℝ) : Prop := ∀ x ∈ l, 0 ≤ x ∧ x ≤ 1

/-- Rank of a diagonalised PSD operator: the number of nonzero eigenvalues. -/
noncomputable def rankOf (l : List ℝ) : ℕ := (l.filter (fun x => !decide (x = 0))).length

/-- `tr P` for a diagonalised operator. -/
def traceOf (l : List ℝ) : ℝ := (l.map id).sum

/-- `‖P‖²_HS` for a diagonalised operator. -/
def hsSqOf (l : List ℝ) : ℝ := (l.map (fun x => x ^ 2)).sum

/-- The right-hand side of the rank–trace bound at `Q = 0`, `b = 0`. -/
def rankTraceRhs (l : List ℝ) : ℝ := 2 * traceOf l - hsSqOf l

/-- Each eigenvalue's contribution `λ(2 − λ)` is at most its contribution to the
rank: at most `1` when `λ ≠ 0`, and exactly `0` when `λ = 0`. This is the whole
content of the `Q = 0` case, isolated as a single-eigenvalue fact. -/
theorem term_le_rank_contrib {x : ℝ} (hx : 0 ≤ x) (hx1 : x ≤ 1) :
    2 * x - x ^ 2 ≤ (if x = 0 then (0:ℝ) else 1) := by
  by_cases h : x = 0
  · subst h; norm_num
  · simp only [h, ite_false]
    nlinarith [sq_nonneg (x - 1)]

/-- The left-hand side is a sum of per-eigenvalue terms. -/
theorem rankTraceRhs_eq_sum (l : List ℝ) :
    rankTraceRhs l = (l.map (fun x => 2 * x - x ^ 2)).sum := by
  unfold rankTraceRhs traceOf hsSqOf
  induction l with
  | nil => simp
  | cons x xs ih => simp only [List.map_cons, List.sum_cons, id_eq] at *; linarith

/-- And so is the rank. -/
theorem rankOf_eq_sum (l : List ℝ) :
    (rankOf l : ℝ) = (l.map (fun x => if x = 0 then (0:ℝ) else 1)).sum := by
  unfold rankOf
  induction l with
  | nil => simp
  | cons x xs ih =>
      by_cases h : x = 0
      · simp [h, ih]
      · simp [h, ih]; linarith

/-- **The rank–trace inequality, `Q = 0` case.**

Under the unit-interval normalization, `rank P ≥ 2 tr P − ‖P‖²_HS`. This is the
form the ENG-006 degree-3 certificate uses (its `Q` is the zero matrix and its
`b` is `0`), so the case that is actually invoked is the case that is proved. -/
theorem rank_trace_zero_Q (l : List ℝ) (h : SpectrumInUnitInterval l) :
    rankTraceRhs l ≤ (rankOf l : ℝ) := by
  rw [rankTraceRhs_eq_sum, rankOf_eq_sum]
  exact List.sum_le_sum (fun i hi => term_le_rank_contrib (h i hi).1 (h i hi).2)

/-- The general statement, as a proposition rather than a theorem.

Carrying it this way is deliberate. It records exactly what the runtime relies
on when `Q ≠ 0`, so a future proof has a fixed target and cannot drift; and it
cannot be mistaken for something proved, because it is a `def` returning `Prop`
and there is no term inhabiting it anywhere in this project. -/
def RankTraceGeneralStatement : Prop :=
  ∀ (lP lQ : List ℝ) (b : ℕ),
    SpectrumInUnitInterval lP →
    (lQ.filter (fun x => decide (0 < x))).length ≤ b →
    2 * traceOf lP + 4 * traceOf lQ - 4 * (b : ℝ) - hsSqOf (lP ++ lQ) ≤ (rankOf lP : ℝ)

end AtlasRH
