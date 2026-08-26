/-
ATLAS-RH-ENG-010 §WO-RH-73 — the 4×4 leading-minor criterion and its
certificate semantics.

The runtime certifies four positive lower bounds on the leading principal
minors of the exactly preconditioned 4×4 even block (and of the shifted pencil
`G - λM`). This file proves what those bounds imply, in the same style as the
ENG-007 3×3 development: one exact completed-square identity carries all the
mathematical content, and everything else is a case split.

The identity is Jacobi's: with `Δ₀ = 1` and `y_k` the leading-minor linear
forms,

    Δ₁Δ₂Δ₃ · xᵀAx = Δ₂Δ₃·y₁² + Δ₃·y₂² + Δ₁·y₃² + Δ₁Δ₂Δ₄·(x 3)².

Every coefficient is a product of leading minors, so positivity of the four
minors makes the right side a positive combination -- Sylvester's criterion in
dimension 4, division-free.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity
import AtlasRH.MatrixInertia
import AtlasRH.CertificateSemantics
import AtlasRH.GeneralizedGap

namespace AtlasRH

open Matrix

/-! ## The symmetric 4×4 and its minors -/

/-- The generic real symmetric 4×4, spelled entry by entry. -/
def sym4 (a b c d e f g h i j : ℝ) : SymMatrix 4 :=
  !![a, b, c, d; b, e, f, g; c, f, h, i; d, g, i, j]

/-- `Δ₂` of `sym4` -- the same quantity `minor2` computes for `sym3`'s leading
block, restated in `sym4`'s entry names. -/
def minor2of4 (a b e : ℝ) : ℝ := a * e - b * b

/-- `Δ₃` of `sym4`: the determinant of its leading 3×3 `[[a,b,c],[b,e,f],[c,f,h]]`. -/
def minor3of4 (a b c e f h : ℝ) : ℝ :=
  a * e * h - a * f ^ 2 - b ^ 2 * h + 2 * b * c * f - c ^ 2 * e

/-- `Δ₄`: the full determinant, expanded once here and nowhere else. -/
def minor4of4 (a b c d e f g h i j : ℝ) : ℝ :=
  a * e * h * j - a * e * i ^ 2 - a * f ^ 2 * j + 2 * a * f * g * i
    - a * g ^ 2 * h - b ^ 2 * h * j + b ^ 2 * i ^ 2 + 2 * b * c * f * j
    - 2 * b * c * g * i - 2 * b * d * f * i + 2 * b * d * g * h
    - c ^ 2 * e * j + c ^ 2 * g ^ 2 + 2 * c * d * e * i - 2 * c * d * f * g
    - d ^ 2 * e * h + d ^ 2 * f ^ 2

theorem sym4_isHermitian (a b c d e f g h i j : ℝ) :
    (sym4 a b c d e f g h i j).IsHermitian := by
  ext r s
  fin_cases r <;> fin_cases s <;> simp [sym4, Matrix.conjTranspose_apply]

theorem qform_sym4 (a b c d e f g h i j : ℝ) (x : Fin 4 → ℝ) :
    qform (sym4 a b c d e f g h i j) x
      = a * x 0 ^ 2 + e * x 1 ^ 2 + h * x 2 ^ 2 + j * x 3 ^ 2
        + 2 * b * x 0 * x 1 + 2 * c * x 0 * x 2 + 2 * d * x 0 * x 3
        + 2 * f * x 1 * x 2 + 2 * g * x 1 * x 3 + 2 * i * x 2 * x 3 := by
  simp [qform, sym4, Matrix.mulVec, dotProduct, Fin.sum_univ_four]
  ring

/-- **The completed-square identity behind the 4×4 criterion** (Jacobi).

An exact polynomial identity -- no division, no hypotheses. The four
coefficients on the right are `Δ₂Δ₃`, `Δ₃`, `Δ₁` and `Δ₁Δ₂Δ₄`, so four
positive minors make the right side a positive combination. The linear forms
are the leading-minor forms: `y₁` is the first row, `y₂` and `y₃` are the
bordered determinants `[rows 1..k, cols 1..k-1, j]`. -/
theorem sym4_completed_square (a b c d e f g h i j : ℝ) (x : Fin 4 → ℝ) :
    a * minor2of4 a b e * minor3of4 a b c e f h
        * qform (sym4 a b c d e f g h i j) x
      = minor2of4 a b e * minor3of4 a b c e f h
            * (a * x 0 + b * x 1 + c * x 2 + d * x 3) ^ 2
        + minor3of4 a b c e f h
            * (minor2of4 a b e * x 1 + (a * f - b * c) * x 2
                + (a * g - b * d) * x 3) ^ 2
        + a * (minor3of4 a b c e f h * x 2
                + (a * e * i - a * f * g - b ^ 2 * i + b * c * g + b * d * f
                    - c * d * e) * x 3) ^ 2
        + a * minor2of4 a b e * minor4of4 a b c d e f g h i j * x 3 ^ 2 := by
  rw [qform_sym4]
  simp only [minor2of4, minor3of4, minor4of4]
  ring

/-- **Sylvester's leading-principal-minor criterion in dimension 4.** -/
theorem posDef_sym4 {a b c d e f g h i j : ℝ} (h1 : 0 < a)
    (h2 : 0 < minor2of4 a b e) (h3 : 0 < minor3of4 a b c e f h)
    (h4 : 0 < minor4of4 a b c d e f g h i j) :
    (sym4 a b c d e f g h i j).PosDef := by
  refine posDef_iff_qform_pos.mpr ⟨sym4_isHermitian a b c d e f g h i j, ?_⟩
  intro x hx
  have hx0 : x 0 ≠ 0 ∨ x 1 ≠ 0 ∨ x 2 ≠ 0 ∨ x 3 ≠ 0 := by
    by_contra hcon
    simp only [not_or, not_not] at hcon
    exact hx (funext fun k => by
      fin_cases k <;> simp [hcon.1, hcon.2.1, hcon.2.2.1, hcon.2.2.2])
  have key := sym4_completed_square a b c d e f g h i j x
  have hlead : 0 < a * minor2of4 a b e * minor3of4 a b c e f h :=
    mul_pos (mul_pos h1 h2) h3
  set m2 := minor2of4 a b e
  set m3 := minor3of4 a b c e f h
  set m4 := minor4of4 a b c d e f g h i j
  set E := a * e * i - a * f * g - b ^ 2 * i + b * c * g + b * d * f - c * d * e
  have hrhs : 0 < m2 * m3 * (a * x 0 + b * x 1 + c * x 2 + d * x 3) ^ 2
      + m3 * (m2 * x 1 + (a * f - b * c) * x 2 + (a * g - b * d) * x 3) ^ 2
      + a * (m3 * x 2 + E * x 3) ^ 2
      + a * m2 * m4 * x 3 ^ 2 := by
    have hT1 : 0 ≤ m2 * m3 * (a * x 0 + b * x 1 + c * x 2 + d * x 3) ^ 2 :=
      mul_nonneg (mul_pos h2 h3).le (sq_nonneg _)
    have hT2 : 0 ≤ m3 * (m2 * x 1 + (a * f - b * c) * x 2
        + (a * g - b * d) * x 3) ^ 2 := mul_nonneg h3.le (sq_nonneg _)
    have hT3 : 0 ≤ a * (m3 * x 2 + E * x 3) ^ 2 :=
      mul_nonneg h1.le (sq_nonneg _)
    have hT4 : 0 ≤ a * m2 * m4 * x 3 ^ 2 :=
      mul_nonneg (mul_pos (mul_pos h1 h2) h4).le (sq_nonneg _)
    rcases eq_or_ne (x 3) 0 with h3z | h3z
    · rcases eq_or_ne (x 2) 0 with h2z | h2z
      · rcases eq_or_ne (x 1) 0 with h1z | h1z
        · have h0z : x 0 ≠ 0 := by
            rcases hx0 with hh | hh | hh | hh
            · exact hh
            · exact absurd h1z hh
            · exact absurd h2z hh
            · exact absurd h3z hh
          have hne : a * x 0 + b * x 1 + c * x 2 + d * x 3 ≠ 0 := by
            rw [h1z, h2z, h3z]
            simpa using mul_ne_zero h1.ne' h0z
          have : 0 < m2 * m3 * (a * x 0 + b * x 1 + c * x 2 + d * x 3) ^ 2 :=
            mul_pos (mul_pos h2 h3) (sq_pos_of_ne_zero hne)
          linarith
        · have hne : m2 * x 1 + (a * f - b * c) * x 2
              + (a * g - b * d) * x 3 ≠ 0 := by
            rw [h2z, h3z]
            simpa using mul_ne_zero h2.ne' h1z
          have : 0 < m3 * (m2 * x 1 + (a * f - b * c) * x 2
              + (a * g - b * d) * x 3) ^ 2 :=
            mul_pos h3 (sq_pos_of_ne_zero hne)
          linarith
      · have hne : m3 * x 2 + E * x 3 ≠ 0 := by
          rw [h3z]
          simpa using mul_ne_zero h3.ne' h2z
        have : 0 < a * (m3 * x 2 + E * x 3) ^ 2 :=
          mul_pos h1 (sq_pos_of_ne_zero hne)
        linarith
    · have : 0 < a * m2 * m4 * x 3 ^ 2 :=
        mul_pos (mul_pos (mul_pos h1 h2) h4) (sq_pos_of_ne_zero h3z)
      linarith
  have := key
  nlinarith [hrhs, hlead, key]

/-! ## The diagonal congruence in dimension 4 -/

def diag4 (w x y z : ℝ) : SymMatrix 4 := !![w, 0, 0, 0; 0, x, 0, 0; 0, 0, y, 0; 0, 0, 0, z]

theorem diag4_eq_diagonal (w x y z : ℝ) :
    diag4 w x y z = Matrix.diagonal ![w, x, y, z] := by
  ext r s
  fin_cases r <;> fin_cases s <;> simp [diag4, Matrix.diagonal]

@[simp] theorem det_diag4 (w x y z : ℝ) : (diag4 w x y z).det = w * x * y * z := by
  rw [diag4_eq_diagonal, Matrix.det_diagonal]
  simp [Fin.prod_univ_four]

theorem isUnit_det_diag4 {w x y z : ℝ} (hw : w ≠ 0) (hx : x ≠ 0) (hy : y ≠ 0)
    (hz : z ≠ 0) : IsUnit (diag4 w x y z).det := by
  rw [det_diag4]
  exact (isUnit_iff_ne_zero).mpr (by
    simpa using mul_ne_zero (mul_ne_zero (mul_ne_zero hw hx) hy) hz)

theorem posDef_of_diagonal_congruence4 {w x y z : ℝ} {A : SymMatrix 4}
    (hw : w ≠ 0) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (h : (congruence (diag4 w x y z) A).PosDef) : A.PosDef :=
  posDef_of_congruence (isUnit_det_diag4 hw hx hy hz) h

/-! ## Certificate semantics -/

/-- The semantic payload of a 4×4 positive-definiteness certificate: a positive
certified lower bound for each leading principal minor (ENG-010 §WO-RH-69). -/
structure PD4Certificate where
  d1Lower : ℝ
  d2Lower : ℝ
  d3Lower : ℝ
  d4Lower : ℝ
  h_d1_pos : 0 < d1Lower
  h_d2_pos : 0 < d2Lower
  h_d3_pos : 0 < d3Lower
  h_d4_pos : 0 < d4Lower

/-- A concrete symmetric 4×4 is *described by* the certificate when each of its
leading principal minors is at least the corresponding certified bound. -/
def PD4Certificate.Describes (cert : PD4Certificate)
    (a b c d e f g h i j : ℝ) : Prop :=
  cert.d1Lower ≤ a ∧
  cert.d2Lower ≤ minor2of4 a b e ∧
  cert.d3Lower ≤ minor3of4 a b c e f h ∧
  cert.d4Lower ≤ minor4of4 a b c d e f g h i j

/-- **The 4×4 certificate implication.** -/
theorem posDef_of_certificate4 (cert : PD4Certificate)
    {a b c d e f g h i j : ℝ}
    (hd : cert.Describes a b c d e f g h i j) :
    (sym4 a b c d e f g h i j).PosDef := by
  obtain ⟨hb1, hb2, hb3, hb4⟩ := hd
  exact posDef_sym4 (lt_of_lt_of_le cert.h_d1_pos hb1)
    (lt_of_lt_of_le cert.h_d2_pos hb2) (lt_of_lt_of_le cert.h_d3_pos hb3)
    (lt_of_lt_of_le cert.h_d4_pos hb4)

/-- **Preconditioned certificate ⟹ the original 4×4 block is positive
definite.** The composition the ENG-010 runtime performs for the block itself. -/
theorem posDef_of_preconditioned_certificate4 (cert : PD4Certificate)
    {w x y z : ℝ} {A : SymMatrix 4} {a b c d e f g h i j : ℝ}
    (hw : w ≠ 0) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (hcong : congruence (diag4 w x y z) A = sym4 a b c d e f g h i j)
    (hd : cert.Describes a b c d e f g h i j) : A.PosDef := by
  refine posDef_of_diagonal_congruence4 hw hx hy hz ?_
  rw [hcong]
  exact posDef_of_certificate4 cert hd

/-- **The composed 4×4 gap implication** (ENG-010 §WO-RH-70): certified minor
bounds on the exactly preconditioned *shifted* pencil `D (G - λM) D` force the
generalized Rayleigh bound for the original pencil at the certified `λ`. -/
theorem gap_of_preconditioned_certificate4 (cert : PD4Certificate)
    {w x y z lam : ℝ} {G M : SymMatrix 4} {a b c d e f g h i j : ℝ}
    (hw : w ≠ 0) (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (hcong : congruence (diag4 w x y z) (G - lam • M)
      = sym4 a b c d e f g h i j)
    (hd : cert.Describes a b c d e f g h i j) (v : Fin 4 → ℝ) :
    lam * qform M v ≤ qform G v :=
  rayleigh_lower_of_shifted_psd
    (posDef_of_preconditioned_certificate4 cert hw hx hy hz hcong hd).posSemidef v

end AtlasRH
