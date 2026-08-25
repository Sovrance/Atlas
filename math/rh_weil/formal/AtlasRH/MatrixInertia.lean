/-
ATLAS-RH-ENG-007 §6 (WO-RH-39) — exact finite linear algebra.

This module formalizes the finite theorems the Python runtime actually leans on
when it reduces a Gram block by congruence and reads a signature off the result:

* congruence by an invertible matrix preserves the quadratic form up to a change
  of variable (`qform_congruence`);
* it therefore preserves the positive and negative indices of the form
  (`posIndexAtLeast_congruence_iff`, `negIndexAtLeast_congruence_iff`) and the
  rank (`rank_congruence`), which together are Sylvester's law of inertia in the
  index formulation the runtime uses;
* the leading-principal-minor criteria in dimensions 2 (in `AtlasRH.Positivity`)
  and 3 (`posDef_sym3_iff` here), which is the criterion the ENG-008 pilot block
  will be tested against.

No RH proof claim is made here or anywhere in this project. Everything below is
finite linear algebra over `ℝ`.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Matrix.ToLin

namespace AtlasRH

open Matrix

variable {n : ℕ}

/-! ## The quadratic form of a symmetric matrix -/

/-- The real quadratic form `x ↦ xᵀ A x`. On `ℝ` this is the form appearing in
Mathlib's `Matrix.PosDef`, because `star` is the identity. -/
def qform (A : SymMatrix n) (x : Fin n → ℝ) : ℝ := x ⬝ᵥ A.mulVec x

theorem qform_eq_star (A : SymMatrix n) (x : Fin n → ℝ) :
    qform A x = star x ⬝ᵥ A.mulVec x := by
  simp [qform, star_trivial]

/-- Positive definiteness restated in terms of `qform`, so the theorems below can
be read without unfolding `star`. -/
theorem posDef_iff_qform_pos {A : SymMatrix n} :
    A.PosDef ↔ A.IsHermitian ∧ ∀ x : Fin n → ℝ, x ≠ 0 → 0 < qform A x := by
  constructor
  · intro h
    exact ⟨h.1, fun x hx => by
      rw [qform_eq_star]; exact h.dotProduct_mulVec_pos hx⟩
  · rintro ⟨h1, h2⟩
    refine Matrix.PosDef.of_dotProduct_mulVec_pos h1 ?_
    intro x hx
    have := h2 x hx
    rwa [qform_eq_star] at this

/-! ## Congruence changes the variable, not the form -/

/-- **Congruence is a change of variable.** `xᵀ (Sᵀ A S) x = (Sx)ᵀ A (Sx)`.

Every step the runtime's LDL elimination performs -- symmetric row/column
elimination and symmetric permutation alike -- is of this shape, so this single
identity is what licenses reading a signature off the reduced matrix. -/
theorem qform_congruence (S A : SymMatrix n) (x : Fin n → ℝ) :
    qform (congruence S A) x = qform A (S.mulVec x) := by
  simp [qform, congruence, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
    Matrix.vecMul_transpose]

/-- The linear automorphism `x ↦ S x` attached to an invertible `S`. -/
noncomputable def congEquiv {S : SymMatrix n} (hS : IsUnit S.det) :
    (Fin n → ℝ) ≃ₗ[ℝ] (Fin n → ℝ) :=
  LinearEquiv.ofLinearMap (Matrix.mulVecLin S) (Matrix.mulVecLin S⁻¹)
    (by rw [← Matrix.mulVecLin_mul, Matrix.mul_nonsing_inv _ hS, Matrix.mulVecLin_one])
    (by rw [← Matrix.mulVecLin_mul, Matrix.nonsing_inv_mul _ hS, Matrix.mulVecLin_one])

@[simp] theorem congEquiv_apply {S : SymMatrix n} (hS : IsUnit S.det) (x : Fin n → ℝ) :
    congEquiv hS x = S.mulVec x := rfl

/-! ## Inertia in the index formulation

The runtime reports a triple `(n₊, n₋, n₀)`. The dimension-theoretic content of
those numbers is: `n₊` is the largest dimension of a subspace on which the form
is positive, `n₋` the same for negative, and `n₀ = n - rank`. Below, each of the
three is shown to be a congruence invariant. -/

/-- The form of `A` is positive on some `k`-dimensional subspace. -/
def PosIndexAtLeast (A : SymMatrix n) (k : ℕ) : Prop :=
  ∃ V : Submodule ℝ (Fin n → ℝ),
    Module.finrank ℝ V = k ∧ ∀ x ∈ V, x ≠ 0 → 0 < qform A x

/-- The form of `A` is negative on some `k`-dimensional subspace. -/
def NegIndexAtLeast (A : SymMatrix n) (k : ℕ) : Prop := PosIndexAtLeast (-A) k

theorem qform_neg (A : SymMatrix n) (x : Fin n → ℝ) : qform (-A) x = -qform A x := by
  simp [qform, Matrix.neg_mulVec]

/-- A positive definite `n × n` matrix has positive index `n`: the whole space
works. This is the bridge from `PosDef` to the signature `(n, 0, 0)`. -/
theorem posIndexAtLeast_of_posDef {A : SymMatrix n} (h : A.PosDef) :
    PosIndexAtLeast A n := by
  refine ⟨⊤, ?_, ?_⟩
  · simp
  · intro x _ hx
    exact (posDef_iff_qform_pos.mp h).2 x hx

/-- The positive index never exceeds the dimension. -/
theorem posIndexAtLeast_le {A : SymMatrix n} {k : ℕ} (h : PosIndexAtLeast A k) :
    k ≤ n := by
  obtain ⟨V, hV, -⟩ := h
  have := Submodule.finrank_le V
  simpa [hV] using this

/-- **Congruence does not decrease the positive index.** -/
theorem posIndexAtLeast_congruence {S A : SymMatrix n} {k : ℕ} (hS : IsUnit S.det)
    (h : PosIndexAtLeast A k) : PosIndexAtLeast (congruence S A) k := by
  obtain ⟨V, hV, hpos⟩ := h
  refine ⟨V.map ((congEquiv hS).symm : (Fin n → ℝ) →ₗ[ℝ] (Fin n → ℝ)), ?_, ?_⟩
  · rw [LinearEquiv.finrank_map_eq]; exact hV
  · rintro x hx hx0
    obtain ⟨v, hv, rfl⟩ := Submodule.mem_map.mp hx
    simp only [LinearEquiv.coe_coe] at hx0 ⊢
    have hv0 : v ≠ 0 := by
      rintro rfl
      exact hx0 (by simp)
    have hq : qform (congruence S A) ((congEquiv hS).symm v) = qform A v := by
      rw [qform_congruence, ← congEquiv_apply hS, LinearEquiv.apply_symm_apply]
    rw [hq]
    exact hpos v hv hv0

/-- Congruence by an invertible `S` is undone by congruence by `S⁻¹`. -/
theorem congruence_symm {S A : SymMatrix n} (hS : IsUnit S.det) :
    congruence S⁻¹ (congruence S A) = A := by
  have hST : IsUnit (Sᵀ).det := by simpa [Matrix.det_transpose] using hS
  simp only [congruence, Matrix.transpose_nonsing_inv]
  have h1 : (Sᵀ)⁻¹ * (Sᵀ * A * S) * S⁻¹ = ((Sᵀ)⁻¹ * Sᵀ) * A * (S * S⁻¹) := by
    simp [Matrix.mul_assoc]
  rw [h1, Matrix.nonsing_inv_mul _ hST, Matrix.mul_nonsing_inv _ hS, Matrix.one_mul,
    Matrix.mul_one]

/-- **Sylvester's law of inertia, positive half.** The positive index is a
congruence invariant. -/
theorem posIndexAtLeast_congruence_iff {S A : SymMatrix n} {k : ℕ} (hS : IsUnit S.det) :
    PosIndexAtLeast (congruence S A) k ↔ PosIndexAtLeast A k := by
  refine ⟨fun h => ?_, posIndexAtLeast_congruence hS⟩
  have hS' : IsUnit (S⁻¹).det := (Matrix.isUnit_nonsing_inv_det S hS)
  have := posIndexAtLeast_congruence (S := S⁻¹) hS' h
  rwa [congruence_symm hS] at this

/-- **Sylvester's law of inertia, negative half.** -/
theorem negIndexAtLeast_congruence_iff {S A : SymMatrix n} {k : ℕ} (hS : IsUnit S.det) :
    NegIndexAtLeast (congruence S A) k ↔ NegIndexAtLeast A k := by
  unfold NegIndexAtLeast
  have hneg : -(congruence S A) = congruence S (-A) := by
    simp [congruence]
  rw [hneg]
  exact posIndexAtLeast_congruence_iff hS

/-- **Sylvester's law of inertia, null half.** Congruence preserves rank, hence
the count of zero eigenvalues `n₀ = n - rank`. -/
theorem rank_congruence {S : SymMatrix n} (A : SymMatrix n) (hS : IsUnit S.det) :
    (congruence S A).rank = A.rank := by
  have hST : IsUnit (Sᵀ).det := by simpa [Matrix.det_transpose] using hS
  simp only [congruence]
  rw [Matrix.rank_mul_eq_left_of_isUnit_det _ _ hS,
    Matrix.rank_mul_eq_right_of_isUnit_det _ _ hST]

/-! ## The 3×3 leading-principal-minor criterion (ENG-007 §6.4)

The ENG-008 pilot block is real symmetric 3×3, and the test the runtime will
apply to it is Sylvester's leading-principal-minor criterion. Both directions
are proved: the forward one lets a certificate conclude definiteness from three
positive minors, the reverse says a failure to clear the test is real
information about the block rather than an artefact of the test. -/

section ThreeByThree

/-- The real symmetric 3×3 matrix with diagonal `a, d, f` and off-diagonal
entries `b = A₀₁`, `c = A₀₂`, `e = A₁₂`. -/
def sym3 (a b c d e f : ℝ) : SymMatrix 3 := !![a, b, c; b, d, e; c, e, f]

/-- Second leading principal minor of `sym3`. -/
def minor2 (a b _c d _e _f : ℝ) : ℝ := a * d - b * b

/-- Third leading principal minor of `sym3`, i.e. its determinant. -/
def minor3 (a b c d e f : ℝ) : ℝ :=
  a * d * f - a * e ^ 2 - b ^ 2 * f + 2 * b * c * e - c ^ 2 * d

theorem sym3_isHermitian (a b c d e f : ℝ) : (sym3 a b c d e f).IsHermitian := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [sym3]

@[simp] theorem det_sym3 (a b c d e f : ℝ) :
    (sym3 a b c d e f).det = minor3 a b c d e f := by
  simp [sym3, minor3, Matrix.det_fin_three]
  ring

/-- The leading 2×2 block of `sym3` is `sym2`. -/
theorem sym3_submatrix_two (a b c d e f : ℝ) :
    (sym3 a b c d e f).submatrix ![0, 1] ![0, 1] = sym2 a b d := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [sym3, sym2]

theorem qform_sym3 (a b c d e f : ℝ) (x : Fin 3 → ℝ) :
    qform (sym3 a b c d e f) x
      = a * x 0 ^ 2 + d * x 1 ^ 2 + f * x 2 ^ 2
        + 2 * b * x 0 * x 1 + 2 * c * x 0 * x 2 + 2 * e * x 1 * x 2 := by
  simp [qform, sym3, Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  ring

/-- **The completed-square identity behind the 3×3 criterion.**

An exact polynomial identity -- no division, no hypotheses -- expressing
`a · D₂ · Q` as a nonnegative combination whose three coefficients are exactly
the three leading principal minors. Everything else in this section is a
case split on which of `x 0, x 1, x 2` is nonzero. -/
theorem sym3_completed_square (a b c d e f : ℝ) (x : Fin 3 → ℝ) :
    a * minor2 a b c d e f * qform (sym3 a b c d e f) x
      = minor2 a b c d e f * (a * x 0 + b * x 1 + c * x 2) ^ 2
        + (minor2 a b c d e f * x 1 + (a * e - b * c) * x 2) ^ 2
        + a * minor3 a b c d e f * x 2 ^ 2 := by
  rw [qform_sym3]
  simp only [minor2, minor3]
  ring

/-- **Sylvester's leading-principal-minor criterion in dimension 3.** -/
theorem posDef_sym3 {a b c d e f : ℝ} (h1 : 0 < a)
    (h2 : 0 < minor2 a b c d e f) (h3 : 0 < minor3 a b c d e f) :
    (sym3 a b c d e f).PosDef := by
  refine posDef_iff_qform_pos.mpr ⟨sym3_isHermitian a b c d e f, ?_⟩
  intro x hx
  have hx0 : x 0 ≠ 0 ∨ x 1 ≠ 0 ∨ x 2 ≠ 0 := by
    by_contra hcon
    simp only [not_or, not_not] at hcon
    exact hx (funext fun i => by fin_cases i <;> simp [hcon.1, hcon.2.1, hcon.2.2])
  have key := sym3_completed_square a b c d e f x
  have hlead : 0 < a * minor2 a b c d e f := mul_pos h1 h2
  -- The right-hand side of `key` is a sum of three nonnegative terms; the case
  -- split exhibits one of them as strictly positive.
  have hrhs : 0 < minor2 a b c d e f * (a * x 0 + b * x 1 + c * x 2) ^ 2
      + (minor2 a b c d e f * x 1 + (a * e - b * c) * x 2) ^ 2
      + a * minor3 a b c d e f * x 2 ^ 2 := by
    rcases eq_or_ne (x 2) 0 with h2z | h2z
    · rcases eq_or_ne (x 1) 0 with h1z | h1z
      · have h0z : x 0 ≠ 0 := by
          rcases hx0 with h | h | h
          · exact h
          · exact absurd h1z h
          · exact absurd h2z h
        rw [h1z, h2z]
        have hne : a * x 0 + b * 0 + c * 0 ≠ 0 := by
          simpa using mul_ne_zero h1.ne' h0z
        nlinarith [mul_pos h2 (sq_pos_of_ne_zero hne)]
      · rw [h2z]
        have hne : minor2 a b c d e f * x 1 + (a * e - b * c) * 0 ≠ 0 := by
          simpa using mul_ne_zero h2.ne' h1z
        have hT1 : 0 ≤ minor2 a b c d e f * (a * x 0 + b * x 1 + c * 0) ^ 2 :=
          mul_nonneg h2.le (sq_nonneg _)
        nlinarith [hT1, sq_pos_of_ne_zero hne]
    · have hT1 : 0 ≤ minor2 a b c d e f * (a * x 0 + b * x 1 + c * x 2) ^ 2 :=
        mul_nonneg h2.le (sq_nonneg _)
      have hT3 : 0 < a * minor3 a b c d e f * x 2 ^ 2 :=
        mul_pos (mul_pos h1 h3) (sq_pos_of_ne_zero h2z)
      nlinarith [hT1, hT3, sq_nonneg (minor2 a b c d e f * x 1 + (a * e - b * c) * x 2)]
  have hprod : 0 < a * minor2 a b c d e f * qform (sym3 a b c d e f) x := by
    rw [key]; exact hrhs
  nlinarith [hprod, hlead]

/-- The converse: a positive definite real symmetric 3×3 has all three leading
principal minors positive. -/
theorem sym3_minors_of_posDef {a b c d e f : ℝ} (h : (sym3 a b c d e f).PosDef) :
    0 < a ∧ 0 < minor2 a b c d e f ∧ 0 < minor3 a b c d e f := by
  refine ⟨?_, ?_, ?_⟩
  · have := Matrix.PosDef.diag_pos h (i := 0)
    simpa [sym3] using this
  · have hinj : Function.Injective (![0, 1] : Fin 2 → Fin 3) := by decide
    have h2 := (h.submatrix hinj)
    rw [sym3_submatrix_two] at h2
    exact (sym2_pos_of_posDef h2).2
  · have := Matrix.PosDef.det_pos h
    simpa using this

/-- The 3×3 criterion as an iff. -/
theorem posDef_sym3_iff {a b c d e f : ℝ} :
    (sym3 a b c d e f).PosDef ↔
      0 < a ∧ 0 < minor2 a b c d e f ∧ 0 < minor3 a b c d e f :=
  ⟨sym3_minors_of_posDef, fun h => posDef_sym3 h.1 h.2.1 h.2.2⟩

/-- A 3×3 block clearing the criterion has inertia `(3, 0, 0)`. -/
theorem definiteInertia_sym3 {a b c d e f : ℝ} (h1 : 0 < a)
    (h2 : 0 < minor2 a b c d e f) (h3 : 0 < minor3 a b c d e f) :
    HasDefiniteInertia (sym3 a b c d e f) :=
  posDef_sym3 h1 h2 h3

end ThreeByThree

end AtlasRH
