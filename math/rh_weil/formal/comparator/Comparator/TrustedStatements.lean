/-
ATLAS-RH-ENG-007 §5 Layer B — the statements Atlas claims to rely on.

Phrased in Mathlib vocabulary only. Nothing here imports `AtlasRH`, and that is the whole
point: if a trusted statement were written in terms of Atlas's own definitions, then
redefining `sym2` or `congruence` would silently redefine what Atlas claims to have proved.
A future proof of a subtly different theorem would then still satisfy the old certificate
consumer -- the exact failure mode ENG-007 §5 Layer D is designed to catch.

These are `def`s of `Prop`, not theorems. They assert nothing on their own; `Solution.lean`
must discharge each one, and the Lean kernel performs the statement comparison by
typechecking that assignment.
-/
import Mathlib.LinearAlgebra.Matrix.PosDef

namespace Comparator.TrustedStatements

open Matrix

/-- Congruence by an invertible matrix preserves positive definiteness, both ways.

This is what licenses the runtime to pivot into a congruent block and read a sign there. -/
def congruence_preserves_posDef : Prop :=
  ∀ {n : Type} [Fintype n] [DecidableEq n] (S A : Matrix n n ℝ), IsUnit S →
    ((Sᴴ * A * S).PosDef ↔ A.PosDef)

/-- Sylvester's 2x2 criterion, as an exact characterisation.

This is what licenses a consumer to read `G₀₀ > 0` and `det > 0` off a certificate and
conclude positive definiteness. -/
def pd_two_by_two : Prop :=
  ∀ a b c : ℝ, ((!![a, b; b, c] : Matrix (Fin 2) (Fin 2) ℝ).PosDef ↔ 0 < a ∧ 0 < a * c - b * b)

/-- A positive Schur pivot transcript implies positive definiteness of a 2x2 block. -/
def schur_pivot_implies_posDef : Prop :=
  ∀ a b c : ℝ, 0 < a → 0 < c - b * b / a →
    (!![a, b; b, c] : Matrix (Fin 2) (Fin 2) ℝ).PosDef

/-- The certificate implication: rigorous lower bounds on the leading entry and the
determinant, holding pointwise over a domain, imply positive definiteness over that domain.

Stated without reference to Atlas's `PD2Certificate` structure, so that changing that
structure cannot change what was claimed. -/
def certificate_even2_implies_pd : Prop :=
  ∀ (g00Lower detLower : ℝ), 0 < g00Lower → 0 < detLower →
    ∀ (G : ℝ → Matrix (Fin 2) (Fin 2) ℝ) (D : Set ℝ),
      (∀ L ∈ D, ∃ a b c : ℝ,
        G L = !![a, b; b, c] ∧ g00Lower ≤ a ∧ detLower ≤ a * c - b * b) →
      ∀ L ∈ D, (G L).PosDef

end Comparator.TrustedStatements
