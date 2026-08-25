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

/-- Even/odd cross terms vanish for any reflection-invariant bilinear form.

This is the structural fact that makes the parity basis worth using: it block-diagonalises
the Gram matrix without any analytic input about how the form was assembled. -/
def even_odd_cross_vanishes : Prop :=
  ∀ {V : Type} [AddCommGroup V] [Module ℝ V]
    (B : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (σ : V →ₗ[ℝ] V),
    (∀ u v : V, B (σ u) (σ v) = B u v) →
    ∀ f g : V, σ f = f → σ g = -g → B f g = 0

/-- The parity factorization of the degree-≤2 Gram determinant: `det G = O1 · E2`.

This is the `det_matches_O1_times_E2` invariant the E1 certifier asserts. -/
def det_parity_factorization : Prop :=
  ∀ g00 g02 g22 o1 : ℝ,
    (!![g00, 0, g02; 0, o1, 0; g02, 0, g22] : Matrix (Fin 3) (Fin 3) ℝ).det
      = o1 * (g00 * g22 - g02 * g02)

/-- A congruence by a determinant-one matrix preserves the Gram determinant.

Licenses computing a determinant in whichever basis is convenient. -/
def det_congruence_invariant : Prop :=
  ∀ {n : Type} [Fintype n] [DecidableEq n] (S G : Matrix n n ℝ),
    S.det = 1 → (Sᵀ * G * S).det = G.det

end Comparator.TrustedStatements
