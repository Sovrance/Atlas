/-
ATLAS-RH-ENG-007 §5 Layer A — canonical finite-dimensional definitions.

Nothing in this file mentions the RH/Weil program. It fixes the finite linear-algebra
vocabulary that the trusted statements are phrased in, so that a statement can be read
without reading the implementation that discharges it.

Scope: finite-dimensional real/complex matrices only. No analytic content, no zeta
function, no proof claim about the Riemann Hypothesis.
-/
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Rank

namespace AtlasRH

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Congruence by `S`: the map `A ↦ Sᴴ A S`.

This is the operation the runtime performs when it pivots a Gram block. It is *not*
similarity: congruence preserves the inertia of a Hermitian form, while similarity
preserves eigenvalues. Conflating the two is the mistake this definition exists to make
impossible to write down by accident. -/
def congruence (S A : Matrix n n ℝ) : Matrix n n ℝ := Sᴴ * A * S

/-- The inertia of a Hermitian matrix: the counts of positive, zero and negative
eigenvalues, as a triple `(n₊, n₀, n₋)`.

Atlas reports inertia as a triple of naturals, e.g. the degree-3 odd block's `(2,0,0)`.
Stated abstractly here so that a claimed signature is a proposition about the matrix
rather than about whatever the runtime happened to print. -/
structure Inertia where
  pos : ℕ
  zero : ℕ
  neg : ℕ
deriving DecidableEq, Repr

/-- A real symmetric 2x2 matrix given by its three independent entries.

The finite Gram blocks Atlas certifies at degree 2 and degree 3 are 2x2, so this concrete
shape carries the certificate-semantics theorems rather than a general `n`. -/
def sym2 (a b c : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![a, b; b, c]

@[simp] lemma sym2_apply_00 (a b c : ℝ) : sym2 a b c 0 0 = a := rfl
@[simp] lemma sym2_apply_01 (a b c : ℝ) : sym2 a b c 0 1 = b := rfl
@[simp] lemma sym2_apply_10 (a b c : ℝ) : sym2 a b c 1 0 = b := rfl
@[simp] lemma sym2_apply_11 (a b c : ℝ) : sym2 a b c 1 1 = c := rfl

lemma sym2_isHermitian (a b c : ℝ) : (sym2 a b c).IsHermitian := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [sym2, Matrix.conjTranspose_apply]

/-- The determinant of a real symmetric 2x2 block, in the form the certificate reports it. -/
lemma sym2_det (a b c : ℝ) : (sym2 a b c).det = a * c - b * b := by
  simp [sym2, Matrix.det_fin_two_of]

end AtlasRH
