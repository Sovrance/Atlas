/-
ATLAS-RH-ENG-007 §6 — exact finite linear algebra.

The theorems here are the ones the Python runtime already relies on when it pivots a Gram
block and reports a sign. Formalizing them does not make the numerics rigorous; it makes
the *implication* from a numerical premise to a mathematical conclusion machine-checked.
-/
import AtlasRH.Definitions

namespace AtlasRH

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-! ### §6.1 Congruence preserves positive definiteness -/

/-- Congruence by an invertible `S` preserves positive definiteness, in both directions.

Mathlib already proves this (`Matrix.IsUnit.posDef_star_left_conjugate_iff`); ENG-007 §6.2
says to prefer an existing theorem over a bespoke one, so this is a named restatement in
Atlas vocabulary rather than a second proof. The `iff` is what matters to the runtime: it
pivots *into* a congruent block and reads the sign off there, which is only sound because
positivity transfers back. -/
theorem congruence_posDef_iff {S A : Matrix n n ℝ} (hS : IsUnit S) :
    (congruence S A).PosDef ↔ A.PosDef := by
  have h : congruence S A = star S * A * S := by
    simp [congruence, star_eq_conjTranspose]
  rw [h]
  exact hS.posDef_star_left_conjugate_iff

/-- The forward direction on its own, which is the shape the pivoting code uses. -/
theorem congruence_posDef {S A : Matrix n n ℝ} (hS : IsUnit S) (hA : A.PosDef) :
    (congruence S A).PosDef :=
  (congruence_posDef_iff hS).mpr hA

/-! ### §6.3 The 2x2 positive-definite criterion

This is the criterion the degree-2 and degree-3 certificates are actually read through:
the runtime reports a lower bound on `G₀₀` and a lower bound on `det G`, and the consumer
concludes positive definiteness. -/

/-- The quadratic form of a real symmetric 2x2 block, on an explicit vector.

`PosDef` in Mathlib is stated over `Finsupp`; this is the plain-vector shape the criterion
is easier to reason about in, bridged by `posDef_iff_dotProduct_mulVec`. -/
lemma sym2_quadratic (a b c : ℝ) (x : Fin 2 → ℝ) :
    star x ⬝ᵥ (sym2 a b c *ᵥ x)
      = a * (x 0) ^ 2 + 2 * b * (x 0) * (x 1) + c * (x 1) ^ 2 := by
  simp [dotProduct, mulVec, Fin.sum_univ_two, sym2, star_trivial]
  ring

/-- Sylvester's criterion for a real symmetric 2x2 matrix, forward direction.

`a > 0` and `ac - b² > 0` imply `!![a, b; b, c]` is positive definite.

The proof is the completed square: `a·Q(x) = (a x₀ + b x₁)² + (ac - b²)·x₁²`. With `a > 0`
the sign of `Q` is decided by two nonnegative terms that cannot both vanish unless `x = 0`. -/
theorem sym2_posDef_of (a b c : ℝ) (ha : 0 < a) (hdet : 0 < a * c - b * b) :
    (sym2 a b c).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos (sym2_isHermitian a b c) (fun x hx => ?_)
  rw [sym2_quadratic]
  have hne : x 0 ≠ 0 ∨ x 1 ≠ 0 := by
    by_contra h
    push_neg at h
    exact hx (by ext i; fin_cases i <;> simp [h.1, h.2])
  rcases eq_or_ne (x 1) 0 with h1 | h1
  · have h0 : x 0 ≠ 0 := by
      rcases hne with h | h
      · exact h
      · exact absurd h1 h
    have hx0 : (0:ℝ) < (x 0) ^ 2 := by positivity
    have hmul : 0 < a * (x 0) ^ 2 := mul_pos ha hx0
    rw [h1]
    simpa using hmul
  · have hx1 : (0:ℝ) < (x 1) ^ 2 := by positivity
    nlinarith [sq_nonneg (a * x 0 + b * x 1), mul_pos hdet hx1]

/-- The converse. ENG-007 §6.3 asks for this direction "if practical"; it is.

The witness is `x = (-b, a)`, on which the form evaluates to `a·(ac - b²)`, so positive
definiteness together with `a > 0` forces the determinant positive. Having both directions
makes the criterion an exact characterisation rather than a sufficient condition that a
certificate consumer might over-read as necessary. -/
theorem sym2_of_posDef {a b c : ℝ} (h : (sym2 a b c).PosDef) :
    0 < a ∧ 0 < a * c - b * b := by
  have ha : 0 < a := by simpa [sym2] using h.diag_pos (i := 0)
  refine ⟨ha, ?_⟩
  have hx : (![-b, a] : Fin 2 → ℝ) ≠ 0 := by
    intro hzero
    have h1 : (![-b, a] : Fin 2 → ℝ) 1 = 0 := by rw [hzero]; rfl
    simp only [Matrix.cons_val_one, Matrix.head_cons] at h1
    exact ha.ne' h1
  have hpos := h.dotProduct_mulVec_pos hx
  rw [sym2_quadratic] at hpos
  have hval : a * (![-b, a] : Fin 2 → ℝ) 0 ^ 2
      + 2 * b * (![-b, a] : Fin 2 → ℝ) 0 * (![-b, a] : Fin 2 → ℝ) 1
      + c * (![-b, a] : Fin 2 → ℝ) 1 ^ 2 = a * (a * c - b * b) := by
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
    ring
  rw [hval] at hpos
  nlinarith [hpos, ha]

/-- The 2x2 criterion as an exact characterisation. -/
theorem sym2_posDef_iff (a b c : ℝ) :
    (sym2 a b c).PosDef ↔ 0 < a ∧ 0 < a * c - b * b :=
  ⟨sym2_of_posDef, fun h => sym2_posDef_of a b c h.1 h.2⟩

end AtlasRH
