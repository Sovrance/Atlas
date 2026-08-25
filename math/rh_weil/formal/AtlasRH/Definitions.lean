/-
ATLAS-RH-ENG-007 §5 Layer A — canonical finite-dimensional definitions.

Only definitions live here. Nothing in this file asserts anything Atlas relies
on; the statements Atlas relies on are in `comparator/TrustedStatements.lean`,
and their proofs are in the `AtlasRH.*` modules. Keeping the three apart is the
whole point of the exercise: a proof that quietly drifts to a weaker theorem
should not be able to satisfy a consumer of the original one.

No RH proof claim is made anywhere in this project. Everything here is finite
linear algebra over `ℝ`.
-/
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.Matrix.Basic

namespace AtlasRH

open Matrix

/-- A real symmetric `n × n` matrix, the only kind of matrix this project
certifies anything about. -/
abbrev SymMatrix (n : ℕ) := Matrix (Fin n) (Fin n) ℝ

/-- Congruence by `S`: the transformation the runtime's LDL elimination performs
at every step. Symmetric row/column elimination and symmetric permutation are
both of this shape, which is why the pivot signature is the matrix signature. -/
def congruence {n : ℕ} (S A : SymMatrix n) : SymMatrix n := Sᵀ * A * S

/-- The 2×2 real symmetric matrix with diagonal `a, c` and off-diagonal `b`.
The even degree-2 block and the odd degree-3 block are both of this shape. -/
def sym2 (a b c : ℝ) : SymMatrix 2 := !![a, b; b, c]

/-- Signature of a real symmetric matrix, as the runtime reports it: counts of
positive, negative and zero eigenvalues. Carried as data so a certificate can
state one without the consumer having to recompute it. -/
structure Inertia where
  nPos : ℕ
  nNeg : ℕ
  nZero : ℕ
deriving DecidableEq, Repr

/-- The signature of a positive definite `n × n` matrix. -/
def Inertia.definite (n : ℕ) : Inertia := ⟨n, 0, 0⟩

/-- `A` has inertia `(n, 0, 0)` exactly when it is positive definite. This is the
bridge between the runtime's signature vocabulary and Mathlib's `PosDef`. -/
def HasDefiniteInertia {n : ℕ} (A : SymMatrix n) : Prop := A.PosDef

/-! ### Midpoint parity (ENG-006 §7)

The Weil basis splits by parity about `x = L/2`. These are the definitions the
parity identities in `AtlasRH.WeilBasis` are stated over. -/

/-- Reflection about the midpoint of `[0, L]`. -/
def reflect (L x : ℝ) : ℝ := L - x

/-- `f` is even about `L/2`. -/
def EvenAbout (L : ℝ) (f : ℝ → ℝ) : Prop := ∀ x, f (reflect L x) = f x

/-- `f` is odd about `L/2`. -/
def OddAbout (L : ℝ) (f : ℝ → ℝ) : Prop := ∀ x, f (reflect L x) = -f x

/-- The constant basis element `1`. -/
def basisOne : ℝ → ℝ := fun _ => 1

/-- The bubble `b(x) = x(L - x)`. -/
def basisB (L : ℝ) : ℝ → ℝ := fun x => x * (L - x)

/-- The midpoint-odd linear element `q1(x) = x - L/2`. -/
noncomputable def basisQ1 (L : ℝ) : ℝ → ℝ := fun x => x - L / 2

/-- The odd degree-3 element `b3(x) = x(L - x)(x - L/2)`. -/
noncomputable def basisB3 (L : ℝ) : ℝ → ℝ := fun x => x * (L - x) * (x - L / 2)

end AtlasRH
