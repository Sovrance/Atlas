/-
ATLAS-RH-ENG-007 §6 — finite positivity: congruence and the 2×2 criterion.

These are the two facts the Python runtime leans on hardest. `inertia/ldl.py`
performs symmetric elimination and reads the signature off the pivots; that is
only valid because congruence preserves definiteness. And every 2×2 block in the
program -- the even degree-2 block, the odd degree-3 block -- has its positivity
decided from its `(0,0)` entry and its determinant, never by an eigenvalue
solver.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions

namespace AtlasRH

open Matrix

section Congruence

variable {n : ℕ}

/-- Congruence by an invertible matrix preserves positive definiteness.

This is the direction the runtime uses: it eliminates, and needs to know the
eliminated matrix is definite exactly when the original was. -/
theorem posDef_congruence {S A : SymMatrix n} (hS : IsUnit S.det) (hA : A.PosDef) :
    (congruence S A).PosDef := by
  have hunit : IsUnit S := (Matrix.isUnit_iff_isUnit_det S).2 hS
  have hinj : Function.Injective S.mulVec := (Matrix.mulVec_injective_iff_isUnit).2 hunit
  have h := Matrix.PosDef.conjTranspose_mul_mul_same hA (B := S) hinj
  simpa [congruence, Matrix.conjTranspose_eq_transpose_of_trivial] using h

/-- Congruence by `S` then by `S⁻¹` is the identity. This is what makes the
converse work: the elimination can be undone. -/
theorem congruence_inv_congruence {S A : SymMatrix n} (hS : IsUnit S.det) :
    congruence S⁻¹ (congruence S A) = A := by
  have hmul : S * S⁻¹ = 1 := Matrix.mul_nonsing_inv S hS
  have htr : (S⁻¹)ᵀ * Sᵀ = 1 := by
    rw [← Matrix.transpose_mul, hmul, Matrix.transpose_one]
  calc congruence S⁻¹ (congruence S A)
      = ((S⁻¹)ᵀ * Sᵀ) * A * (S * S⁻¹) := by
        simp [congruence, Matrix.mul_assoc]
    _ = A := by rw [htr, hmul]; simp

/-- The converse: if the congruence is definite then so was the original.
Together with `posDef_congruence` this is the `↔` that makes a pivot signature a
matrix signature rather than an artifact of the elimination order. -/
theorem posDef_of_congruence {S A : SymMatrix n} (hS : IsUnit S.det)
    (h : (congruence S A).PosDef) : A.PosDef := by
  have hSinv : IsUnit (S⁻¹).det := Matrix.isUnit_nonsing_inv_det S hS
  have := posDef_congruence (S := S⁻¹) hSinv h
  rwa [congruence_inv_congruence hS] at this

/-- Congruence-invariance of definiteness, as an iff. -/
theorem posDef_congruence_iff {S A : SymMatrix n} (hS : IsUnit S.det) :
    (congruence S A).PosDef ↔ A.PosDef :=
  ⟨posDef_of_congruence hS, posDef_congruence hS⟩

end Congruence

section TwoByTwo

/-- `sym2` is symmetric. -/
theorem sym2_isHermitian (a b c : ℝ) : (sym2 a b c).IsHermitian := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [sym2]

/-- Determinant of the 2×2 symmetric block. -/
@[simp] theorem det_sym2 (a b c : ℝ) : (sym2 a b c).det = a * c - b * b := by
  simp [sym2, Matrix.det_fin_two]

/-- **The 2×2 positive-definite criterion.**

`a > 0` together with `ac - b² > 0` gives positive definiteness. This is exactly
the test `inertia.congruence.inertia_2x2` applies, and exactly what the degree-2
and degree-3 certificates report: a positive `(0,0)` entry and a positive
determinant. -/
theorem posDef_sym2 {a b c : ℝ} (ha : 0 < a) (hdet : 0 < a * c - b * b) :
    (sym2 a b c).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos (sym2_isHermitian a b c) ?_
  intro x hx
  have hx0 : x 0 ≠ 0 ∨ x 1 ≠ 0 := by
    by_contra hcon
    simp only [not_or, not_not] at hcon
    exact hx (funext fun i => by fin_cases i <;> simp [hcon.1, hcon.2])
  have expand : star x ⬝ᵥ (sym2 a b c).mulVec x
      = a * (x 0) ^ 2 + 2 * b * (x 0) * (x 1) + c * (x 1) ^ 2 := by
    simp [sym2, Matrix.mulVec, dotProduct, Fin.sum_univ_two, star_trivial]
    ring
  rw [expand]
  -- The completed square: a * Q = (a x0 + b x1)^2 + (ac - b^2) x1^2, so a > 0
  -- and a positive determinant force Q > 0 unless x vanishes entirely.
  rcases hx0 with h0 | h1
  · rcases eq_or_ne (x 1) 0 with h1 | h1
    · have hx0sq : 0 < (x 0) ^ 2 := sq_pos_of_ne_zero h0
      have : a * (x 0) ^ 2 + 2 * b * (x 0) * (x 1) + c * (x 1) ^ 2 = a * (x 0) ^ 2 := by
        rw [h1]; ring
      rw [this]; positivity
    · have h1sq : 0 < (x 1) ^ 2 := sq_pos_of_ne_zero h1
      nlinarith [sq_nonneg (a * x 0 + b * x 1), mul_pos hdet h1sq]
  · have h1sq : 0 < (x 1) ^ 2 := sq_pos_of_ne_zero h1
    nlinarith [sq_nonneg (a * x 0 + b * x 1), mul_pos hdet h1sq]

/-- The converse: a positive definite 2×2 has positive `(0,0)` entry and
positive determinant. Both directions matter -- the forward one lets a
certificate conclude definiteness, the reverse says the test is not merely
sufficient, so a failure to clear it is genuine information. -/
theorem sym2_pos_of_posDef {a b c : ℝ} (h : (sym2 a b c).PosDef) :
    0 < a ∧ 0 < a * c - b * b := by
  have ha : 0 < a := by
    have := Matrix.PosDef.diag_pos h (i := 0)
    simpa [sym2] using this
  refine ⟨ha, ?_⟩
  -- Evaluate the form at (-b, a): it equals a * (ac - b²), which is positive,
  -- and a > 0. No determinant lemma is needed.
  have hx : (![-b, a] : Fin 2 → ℝ) ≠ 0 := by
    intro hz
    have : (![-b, a] : Fin 2 → ℝ) 1 = 0 := by rw [hz]; rfl
    simp at this
    exact ha.ne' this
  have hpos := Matrix.PosDef.dotProduct_mulVec_pos h hx
  have expand : star (![-b, a] : Fin 2 → ℝ) ⬝ᵥ (sym2 a b c).mulVec ![-b, a]
      = a * (a * c - b * b) := by
    simp [sym2, Matrix.mulVec, dotProduct, Fin.sum_univ_two, star_trivial]
    ring
  rw [expand] at hpos
  nlinarith [hpos, ha]

/-- The 2×2 criterion as an iff. -/
theorem posDef_sym2_iff {a b c : ℝ} :
    (sym2 a b c).PosDef ↔ 0 < a ∧ 0 < a * c - b * b :=
  ⟨sym2_pos_of_posDef, fun h => posDef_sym2 h.1 h.2⟩

/-- A positive definite 2×2 block has inertia `(2, 0, 0)` -- the signature the
ENG-006 degree-3 certificate reports. -/
theorem definiteInertia_sym2 {a b c : ℝ} (ha : 0 < a) (hdet : 0 < a * c - b * b) :
    HasDefiniteInertia (sym2 a b c) :=
  posDef_sym2 ha hdet

end TwoByTwo

end AtlasRH
