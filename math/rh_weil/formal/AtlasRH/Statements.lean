/-
ATLAS-RH-ENG-007 §5 Layer B — trusted statements.

This file holds the exact propositions Atlas claims to rely on. It is the
comparator's fixed target, so two rules govern it:

1. **It imports `AtlasRH.Definitions` and Mathlib, and nothing else.** In
   particular it imports no module that contains a proof. A statement written
   in terms of an implementation lemma's own vocabulary can be weakened by
   editing that vocabulary, which is exactly the drift this layer exists to
   catch. `scripts/check_formal_manifest.py` enforces the import list.

2. **Each statement is spelled out.** Where an implementation module has a
   convenient abbreviation -- `qform`, `PosIndexAtLeast`, `sym3`, `minor2`,
   `PD2Certificate` -- the statement below expands it. The text is longer and
   that is the point: a reader checking whether Atlas proved what Atlas says it
   proved should not have to chase definitions through the proof library.

`comparator/Solution.lean` proves every one of these, and Lean itself performs
the Layer D comparison: a proof whose proposition drifted would no longer
typecheck against the name it claims to inhabit.

No RH proof claim is made here. Every statement below is finite linear algebra
over `ℝ`.
-/
import AtlasRH.Definitions
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Rank

namespace AtlasRH

open Matrix

/-! ## Sylvester's law of inertia -/

/-- **Congruence preserves the positive index.**

The positive index of a real symmetric form is the largest dimension of a
subspace on which the form is positive. Congruence by an invertible `S` leaves
it unchanged. This is what makes the signature the runtime's LDL elimination
reports a property of the matrix rather than of the elimination order. -/
def InertiaCongruencePositiveStatement : Prop :=
  ∀ (n : ℕ) (S A : SymMatrix n) (k : ℕ), IsUnit S.det →
    ((∃ V : Submodule ℝ (Fin n → ℝ), Module.finrank ℝ V = k ∧
        ∀ x ∈ V, x ≠ 0 → 0 < x ⬝ᵥ (congruence S A).mulVec x)
      ↔ (∃ V : Submodule ℝ (Fin n → ℝ), Module.finrank ℝ V = k ∧
        ∀ x ∈ V, x ≠ 0 → 0 < x ⬝ᵥ A.mulVec x))

/-- **Congruence preserves the negative index.** Same statement for `-A`. -/
def InertiaCongruenceNegativeStatement : Prop :=
  ∀ (n : ℕ) (S A : SymMatrix n) (k : ℕ), IsUnit S.det →
    ((∃ V : Submodule ℝ (Fin n → ℝ), Module.finrank ℝ V = k ∧
        ∀ x ∈ V, x ≠ 0 → 0 < x ⬝ᵥ (congruence S (-A)).mulVec x)
      ↔ (∃ V : Submodule ℝ (Fin n → ℝ), Module.finrank ℝ V = k ∧
        ∀ x ∈ V, x ≠ 0 → 0 < x ⬝ᵥ (-A).mulVec x))

/-- **Congruence preserves rank**, hence the zero count `n₀ = n - rank`. -/
def InertiaCongruenceRankStatement : Prop :=
  ∀ (n : ℕ) (S A : SymMatrix n), IsUnit S.det → (congruence S A).rank = A.rank

/-! ## Leading-principal-minor criteria -/

/-- **The 2×2 criterion**, as an iff. Exactly the test the runtime's
`inertia_2x2` applies, and exactly what the degree-2 and degree-3 certificates
report: a positive leading entry and a positive determinant. -/
def PdTwoByTwoStatement : Prop :=
  ∀ a b c : ℝ,
    (!![a, b; b, c] : SymMatrix 2).PosDef ↔ (0 < a ∧ 0 < a * c - b * b)

/-- **The 3×3 leading-principal-minor criterion**, as an iff. This is the test
the ENG-008 higher-dimensional pilot block is to be judged by. -/
def PdThreeByThreeStatement : Prop :=
  ∀ a b c d e f : ℝ,
    (!![a, b, c; b, d, e; c, e, f] : SymMatrix 3).PosDef ↔
      (0 < a ∧ 0 < a * d - b * b ∧
        0 < a * d * f - a * e ^ 2 - b ^ 2 * f + 2 * b * c * e - c ^ 2 * d)

/-! ## Certificate semantics -/

/-- **The 2×2 certificate implication.**

Given honest enclosures for the three distinct entries of a real symmetric 2×2
matrix, a positive certified lower bound on the leading entry, and a positive
certified lower bound on the determinant taken at the worst corner of those
enclosures, the matrix is positive definite.

This is the implication that gives a rigorous interval run its meaning: Arb
establishes the premises, and this says they suffice. It converts nothing --
the numeric warrant stays E1; only the implication is FORMAL. -/
def CertificateEven2ImpliesPdStatement : Prop :=
  ∀ a b d g00lo g00hi g0blo g0bhi gbblo gbbhi g00Lower detLower : ℝ,
    0 < g00Lower →
    0 < detLower →
    g00Lower ≤ g00lo →
    detLower ≤ g00lo * gbblo - max (g0blo ^ 2) (g0bhi ^ 2) →
    g00lo ≤ a → a ≤ g00hi →
    g0blo ≤ b → b ≤ g0bhi →
    gbblo ≤ d → d ≤ gbbhi →
    (!![a, b; b, d] : SymMatrix 2).PosDef

/-! ## RH/Weil exact identities -/

/-- **Parity of the Weil basis about the midpoint `L/2`.** -/
def WeilBasisParityStatement : Prop :=
  ∀ L : ℝ,
    EvenAbout L basisOne ∧ EvenAbout L (basisB L) ∧
      OddAbout L (basisQ1 L) ∧ OddAbout L (basisB3 L)

/-- **Cross-parity blocks vanish.**

For any pairing that is linear in its second argument and invariant under
simultaneous reflection of both arguments, an even function pairs to zero with
an odd one. Nothing about digamma, primes or integration enters -- which is why
this is provable here and the numerical Gram entries are not. -/
def OddDegree3CrossBlockStatement : Prop :=
  ∀ (L : ℝ) (pair : (ℝ → ℝ) → (ℝ → ℝ) → ℝ),
    (∀ f g, pair f (fun x => -g x) = -pair f g) →
    (∀ f g, pair (fun x => f (reflect L x)) (fun x => g (reflect L x)) = pair f g) →
    ∀ f g : ℝ → ℝ, EvenAbout L f → OddAbout L g → pair f g = 0

/-- **Degree-2 parity factorization.** With the cross entries gone, the
determinant of the parity-split Gram is the product of the two block
determinants. -/
def OddDegree3FactorizationStatement : Prop :=
  ∀ e00 e0b ebb o11 o1b obb : ℝ,
    (!![e00, e0b; e0b, ebb] : SymMatrix 2).det *
        (!![o11, o1b; o1b, obb] : SymMatrix 2).det
      = (e00 * ebb - e0b ^ 2) * (o11 * obb - o1b ^ 2)

/-! ## Rank–trace -/

/-- **The rank–trace inequality, `Q = 0` case.**

For a spectrum in `[0, 1]`, `2·tr P − ‖P‖²_HS ≤ rank P`. This is the case the
ENG-006 degree-3 certificate actually uses: its `Q` is the zero matrix and its
bound is `b = 0`.

Rank, trace and the Hilbert–Schmidt norm are spelled out on the eigenvalue list
so that no implementation abbreviation appears in the statement. -/
def RankTraceZeroQStatement : Prop :=
  ∀ l : List ℝ, (∀ x ∈ l, 0 ≤ x ∧ x ≤ 1) →
    2 * (l.map id).sum - (l.map (fun x => x ^ 2)).sum
      ≤ ((l.filter (fun x => !decide (x = 0))).length : ℝ)

end AtlasRH
