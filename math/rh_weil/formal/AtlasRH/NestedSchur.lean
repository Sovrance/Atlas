/-
ATLAS-RH-ENG-011 §WO-RH-85 — shifted-positivity composition, nested-gap
regression, and the Schur witness theorem.

Three implications the 5×5 runtime uses, each stated at the same trust
boundary as always: Arb establishes the premises, Lean proves they suffice.

* **Shifted positivity transfer.** The n = 5 certificates are LDL sweeps of
  the *shifted* pencil `G − λM`. Two theorems turn one sweep into three
  conclusions: a smaller shift stays positive definite (which is how the
  ENG-010 4×4 gap certificate licenses `G₄ − λ₅M₄ ≻ 0` at the smaller `λ₅`
  with no new cover), and `G` itself is positive definite because
  `G = (G − λM) + λM` with `λ > 0` and `M ≻ 0`.

* **Nested-gap regression.** `M₄` and `G₄` are exactly the leading blocks of
  `M₅` and `G₅`, so a shifted-PSD certificate at n = 5 *restricts*: the same
  `λ` works for every sub-pencil. This is the formal direction of generalized
  Cauchy interlacing the runtime quotes — certified gap lower bounds regress
  upward through nesting.

* **The Schur witness theorem.** If `A ≻ 0`, `A y = c`, and `d − c ⬝ᵥ y > 0`,
  the bordered quadratic form `q(x, t) = xᵀAx + 2t·(c ⬝ᵥ x) + d·t²` — which
  is the quadratic form of `[[A, c], [cᵀ, d]]` — is positive away from zero:
  it equals `(x + t·y)ᵀ A (x + t·y) + (d − c ⬝ᵥ y)·t²` exactly. This is the
  semantic payload of the §WO-RH-79 verified-solve analysis; the runtime's
  interval solve encloses the true `y`, and the hypotheses hold for it.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity
import AtlasRH.MatrixInertia
import AtlasRH.GeneralizedGap

namespace AtlasRH

open Matrix

variable {n : ℕ}

/-- Nonnegative scaling preserves positive semidefiniteness, spelled out for
real matrices so no star-order instance is needed. -/
theorem smul_posSemidef {M : SymMatrix n} (hM : M.PosSemidef) {c : ℝ}
    (hc : 0 ≤ c) : (c • M).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec] at hM ⊢
  refine ⟨?_, fun x => ?_⟩
  · have ht : Mᵀ = M := by
      ext i j
      have := congrFun (congrFun hM.1.eq i) j
      simpa [Matrix.conjTranspose_apply, Matrix.transpose_apply] using this
    ext i j
    simp [Matrix.conjTranspose_apply]
    exact Or.inl (by simpa using congrFun (congrFun ht i) j)
  · have := hM.2 x
    have hval : star x ⬝ᵥ (c • M).mulVec x = c * (star x ⬝ᵥ M.mulVec x) := by
      simp [Matrix.smul_mulVec, dotProduct_smul]
    rw [hval]
    exact mul_nonneg hc this

/-- A positive-definite matrix plus a positive-semidefinite one is positive
definite, on `qform` where everything downstream consumes it. -/
theorem posDef_add_posSemidef {A B : SymMatrix n}
    (hA : A.PosDef) (hB : B.PosSemidef) : (A + B).PosDef := by
  refine posDef_iff_qform_pos.mpr ⟨hA.1.add hB.1, ?_⟩
  intro x hx
  have h1 : 0 < qform A x := (posDef_iff_qform_pos.mp hA).2 x hx
  have h2 : 0 ≤ qform B x := by
    simpa [qform_eq_star] using hB.dotProduct_mulVec_nonneg x
  have hsum : qform (A + B) x = qform A x + qform B x := by
    simp [qform, add_mulVec, dotProduct_add]
  linarith [hsum.le, hsum.symm.le]

/-- **Shifted positivity composes back.** `G = (G − λM) + λM`: a certified
shifted sweep, a positive exact `λ`, and the E0-positive reference metric
force `G ≻ 0` with no additional cover. -/
theorem posDef_of_shifted_posDef_add {G M : SymMatrix n} {lam : ℝ}
    (hlam : 0 < lam) (hM : M.PosDef) (h : (G - lam • M).PosDef) : G.PosDef := by
  have hG : G = (G - lam • M) + lam • M := by abel
  rw [hG]
  exact posDef_add_posSemidef h (smul_posSemidef hM.posSemidef hlam.le)

/-- **A smaller shift stays positive definite.**
`G − bM = (G − aM) + (a − b)M` with `a − b ≥ 0` and `M ⪰ 0`. This is how the
4×4 gap certificate at `λ₄` licenses the leading shifted block at any
`λ₅ ≤ λ₄` — the theorem the n = 5 Schur analysis rests on. -/
theorem shifted_posDef_of_le {G M : SymMatrix n} {a b : ℝ}
    (hba : b ≤ a) (hM : M.PosSemidef) (h : (G - a • M).PosDef) :
    (G - b • M).PosDef := by
  have hG : G - b • M = (G - a • M) + (a - b) • M := by
    rw [sub_smul]; abel
  rw [hG]
  exact posDef_add_posSemidef h (smul_posSemidef hM (by linarith))

/-- **Nested-gap regression.** A shifted-PSD certificate restricts along any
index map: the certified `λ` transfers to every sub-pencil, in particular to
the leading 4×4 pair inside the 5×5. -/
theorem gap_bound_restricts {m : ℕ} {G M : SymMatrix n} {lam : ℝ}
    (h : (G - lam • M).PosSemidef) (e : Fin m → Fin n) :
    (G.submatrix e e - lam • M.submatrix e e).PosSemidef := by
  have := h.submatrix e
  simpa [Matrix.submatrix_sub, Matrix.submatrix_smul] using this

/-- The Rayleigh form of the restriction: certified gap lower bounds regress
upward through nesting — the formal regression bound §WO-RH-84 asked for. -/
theorem rayleigh_lower_restricts {m : ℕ} {G M : SymMatrix n} {lam : ℝ}
    (h : (G - lam • M).PosSemidef) (e : Fin m → Fin n) (x : Fin m → ℝ) :
    lam * qform (M.submatrix e e) x ≤ qform (G.submatrix e e) x :=
  rayleigh_lower_of_shifted_psd (gap_bound_restricts h e) x

/-! ## The Schur witness theorem -/

/-- The bordered quadratic form of `[[A, c], [cᵀ, d]]`, written out. -/
def borderedForm (A : SymMatrix n) (c : Fin n → ℝ) (d : ℝ)
    (x : Fin n → ℝ) (t : ℝ) : ℝ :=
  qform A x + 2 * t * (c ⬝ᵥ x) + d * t ^ 2

/-- Real symmetry moves `A` across the dot product. -/
theorem dotProduct_mulVec_symm {A : SymMatrix n} (hAt : Aᵀ = A)
    (v w : Fin n → ℝ) : v ⬝ᵥ A.mulVec w = w ⬝ᵥ A.mulVec v := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hAt,
      dotProduct_comm]

/-- **The Schur witness identity.** With `A y = c` and `A` symmetric,
`q(x, t) = qform A (x + t·y) + (d − c ⬝ᵥ y)·t²` exactly. -/
theorem borderedForm_eq (A : SymMatrix n) (hAt : Aᵀ = A)
    {c : Fin n → ℝ} (d : ℝ) {y : Fin n → ℝ} (hy : A.mulVec y = c)
    (x : Fin n → ℝ) (t : ℝ) :
    borderedForm A c d x t
      = qform A (x + t • y) + (d - c ⬝ᵥ y) * t ^ 2 := by
  have hxy : x ⬝ᵥ A.mulVec y = c ⬝ᵥ x := by
    rw [hy, dotProduct_comm]
  have hyx : y ⬝ᵥ A.mulVec x = c ⬝ᵥ x := by
    rw [dotProduct_mulVec_symm hAt, hy, dotProduct_comm]
  have hyy : y ⬝ᵥ A.mulVec y = c ⬝ᵥ y := by
    rw [hy, dotProduct_comm]
  have hexpand : qform A (x + t • y)
      = qform A x + t * (x ⬝ᵥ A.mulVec y) + t * (y ⬝ᵥ A.mulVec x)
        + t ^ 2 * (y ⬝ᵥ A.mulVec y) := by
    simp only [qform, Matrix.mulVec_add, Matrix.mulVec_smul, dotProduct_add,
               add_dotProduct, smul_dotProduct, dotProduct_smul, smul_eq_mul]
    ring
  rw [hexpand, hxy, hyx, hyy]
  simp only [borderedForm]
  ring

/-- **The Schur witness theorem.** `A ≻ 0`, an exact witness `A y = c`, and a
positive residual `d − c ⬝ᵥ y` make the bordered form positive away from zero.
This is what the §WO-RH-79 certified solves establish about adding the `b⁴`
direction to a positive block. -/
theorem borderedForm_pos (A : SymMatrix n) (hA : A.PosDef)
    {c : Fin n → ℝ} {d : ℝ} {y : Fin n → ℝ} (hy : A.mulVec y = c)
    (hd : 0 < d - c ⬝ᵥ y)
    (x : Fin n → ℝ) (t : ℝ) (hxt : x ≠ 0 ∨ t ≠ 0) :
    0 < borderedForm A c d x t := by
  have hAt : Aᵀ = A := by
    ext i j
    have := congrFun (congrFun hA.1.eq i) j
    simpa [Matrix.conjTranspose_apply, Matrix.transpose_apply] using this
  rw [borderedForm_eq A hAt d hy x t]
  have hq : 0 ≤ qform A (x + t • y) := by
    rcases eq_or_ne (x + t • y) 0 with hz | hz
    · rw [hz]; simp [qform]
    · exact ((posDef_iff_qform_pos.mp hA).2 _ hz).le
  rcases eq_or_ne t 0 with ht | ht
  · have hx : x ≠ 0 := by
      rcases hxt with hx | hcon
      · exact hx
      · exact absurd ht hcon
    have hz : x + t • y ≠ 0 := by
      rw [ht]; simpa using hx
    have : 0 < qform A (x + t • y) := (posDef_iff_qform_pos.mp hA).2 _ hz
    nlinarith [sq_nonneg t]
  · have : 0 < (d - c ⬝ᵥ y) * t ^ 2 :=
      mul_pos hd (by positivity)
    linarith
end AtlasRH
