/-
ATLAS-RH-ENG-009 §WO-RH-63 — the generalized-gap implication.

The runtime certifies, by interval covers, that `G - λ M` is positive
(semi)definite for a concrete dyadic `λ` and the exact L² reference metric `M`.
This file proves the two implications that give that computation its meaning:

* **The Rayleigh bound.** If `G - λ M` is positive semidefinite, then
  `λ * xᵀMx ≤ xᵀGx` for every `x`. With `M` positive definite this is exactly
  "every generalized eigenvalue of the pencil `(G, M)` is at least `λ`", stated
  without eigenvalues.

* **Simultaneous-congruence invariance.** `Sᵀ(G - λM)S = SᵀGS - λ(SᵀMS)`, so a
  shifted-positivity certificate transports across any invertible change of
  basis applied to both forms at once. This is the precise sense in which the
  generalized gap is basis-invariant while raw eigenvalues and raw determinants
  are not.

As everywhere in this development, Lean proves the implication and Arb is
responsible for the premises. No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity
import AtlasRH.MatrixInertia
import AtlasRH.CertificateSemantics

namespace AtlasRH

open Matrix

variable {n : ℕ}

/-- The quadratic form is linear in the matrix argument through a shift:
`qform (G - λ • M) x = qform G x - λ * qform M x`. -/
theorem qform_sub_smul (G M : SymMatrix n) (lam : ℝ) (x : Fin n → ℝ) :
    qform (G - lam • M) x = qform G x - lam * qform M x := by
  simp [qform, sub_mulVec, Matrix.smul_mulVec, dotProduct_sub,
        dotProduct_smul, smul_eq_mul]

/-- **The Rayleigh lower bound.** Shifted positive semidefiniteness pushes the
quadratic form of `G` above `λ` times that of `M`, for every vector. This is
the implication the generalized-gap certificate rests on: the interval run
establishes the premise for a concrete `λ`, and this theorem is what that
premise means. -/
theorem rayleigh_lower_of_shifted_psd {G M : SymMatrix n} {lam : ℝ}
    (h : (G - lam • M).PosSemidef) (x : Fin n → ℝ) :
    lam * qform M x ≤ qform G x := by
  have h0 : 0 ≤ qform (G - lam • M) x := by
    simpa [qform_eq_star] using h.dotProduct_mulVec_nonneg x
  have := qform_sub_smul G M lam x
  linarith [h0, this.symm.le, this.le]

/-- Strict version: shifted *definiteness* gives a strict Rayleigh bound on
nonzero vectors. -/
theorem rayleigh_lower_of_shifted_posDef {G M : SymMatrix n} {lam : ℝ}
    (h : (G - lam • M).PosDef) {x : Fin n → ℝ} (hx : x ≠ 0) :
    lam * qform M x < qform G x := by
  have h0 : 0 < qform (G - lam • M) x :=
    (posDef_iff_qform_pos.mp h).2 x hx
  have := qform_sub_smul G M lam x
  linarith [h0, this.le, this.symm.le]

/-- The pencil shift commutes with simultaneous congruence:
`Sᵀ(G - λM)S = SᵀGS - λ(SᵀMS)`. Pure matrix algebra, and the identity that
makes "generalized eigenvalues are basis-invariant" a theorem rather than a
slogan. -/
theorem congruence_sub_smul (S G M : SymMatrix n) (lam : ℝ) :
    congruence S (G - lam • M) = congruence S G - lam • congruence S M := by
  simp [congruence, Matrix.mul_sub, Matrix.sub_mul]

/-- **Invariance of the certified gap.** A shifted-positivity certificate for
`(G, M)` yields one for `(SᵀGS, SᵀMS)` under any invertible `S` — both forms
transformed together. The certified lower bound `λ` itself is unchanged. -/
theorem shifted_posDef_congruence {S G M : SymMatrix n} {lam : ℝ}
    (hS : IsUnit S.det) (h : (G - lam • M).PosDef) :
    (congruence S G - lam • congruence S M).PosDef := by
  rw [← congruence_sub_smul]
  exact posDef_congruence hS h

/-- And back: certifying in any congruent coordinate system certifies the
original pencil. Together with `shifted_posDef_congruence` this says the two
statements are one statement. -/
theorem shifted_posDef_of_congruence {S G M : SymMatrix n} {lam : ℝ}
    (hS : IsUnit S.det)
    (h : (congruence S G - lam • congruence S M).PosDef) :
    (G - lam • M).PosDef := by
  rw [← congruence_sub_smul] at h
  exact posDef_of_congruence hS h

/-- **The composed 3×3 gap implication.** What the ENG-009 runtime actually
establishes for the 3×3 even block: certified positive lower bounds on the
three leading minors of the exactly-preconditioned shifted block. Composing the
ENG-008 certificate semantics with the Rayleigh bound: those minors force the
generalized Rayleigh bound for the *original* pencil, at the certified `λ`. -/
theorem gap_of_preconditioned_certificate3
    (c : PD3Certificate) {x y z lam : ℝ} {G M : SymMatrix 3}
    {a b cc d e f : ℝ}
    (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (hcong : congruence (diag3 x y z) (G - lam • M) = sym3 a b cc d e f)
    (h : c.Describes a b cc d e f) (v : Fin 3 → ℝ) :
    lam * qform M v ≤ qform G v :=
  rayleigh_lower_of_shifted_psd
    (posDef_of_preconditioned_certificate3 c hx hy hz hcong h).posSemidef v

end AtlasRH
