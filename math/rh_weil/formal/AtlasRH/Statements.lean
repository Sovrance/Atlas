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

/-! ## The 3×3 certificate and its preconditioner (ENG-008) -/

/-- **The 3×3 certificate implication.**

Certified positive lower bounds on the three leading principal minors imply the
block is positive definite. This is what the ENG-008 degree-4 result rests on:
the interval run establishes the three bounds, and this says they suffice.

Stated on the bounds rather than on entry enclosures because that is where the
runtime establishes them. At 2×2 the worst corner of a determinant is findable
in closed form, so that certificate can carry enclosures; at 3×3 it is not, and
encoding interval arithmetic here would be the wrong division of labour. -/
def PdThreeByThreeCertificateStatement : Prop :=
  ∀ a b c d e f d1Lower d2Lower d3Lower : ℝ,
    0 < d1Lower →
    0 < d2Lower →
    0 < d3Lower →
    d1Lower ≤ a →
    d2Lower ≤ a * d - b * b →
    d3Lower ≤ a * d * f - a * e ^ 2 - b ^ 2 * f + 2 * b * c * e - c ^ 2 * d →
    (!![a, b, c; b, d, e; c, e, f] : SymMatrix 3).PosDef

/-- **A diagonal preconditioner does not change the answer.**

If the rescaled block is positive definite then so is the block itself. The
ENG-008 runtime rescales by an exact dyadic diagonal before eliminating, because
the raw block spans ten orders of magnitude; this is the step that carries the
conclusion back. Nothing is needed of the preconditioner beyond its diagonal
entries being nonzero. -/
def DiagonalCongruencePreservesPdStatement : Prop :=
  ∀ (x y z : ℝ) (A : SymMatrix 3),
    x ≠ 0 → y ≠ 0 → z ≠ 0 →
    (congruence (!![x, 0, 0; 0, y, 0; 0, 0, z] : SymMatrix 3) A).PosDef →
    A.PosDef

/-- **The composition the runtime actually performs.**

Certified positive lower bounds on the three leading minors of the *rescaled*
block imply the *original* block is positive definite. The two steps above,
composed -- which is worth stating as one theorem because the composition is
what a reader of the certificate has to trust. -/
def PreconditionedCertificate3Statement : Prop :=
  ∀ (x y z : ℝ) (A : SymMatrix 3) (a b c d e f d1Lower d2Lower d3Lower : ℝ),
    x ≠ 0 → y ≠ 0 → z ≠ 0 →
    congruence (!![x, 0, 0; 0, y, 0; 0, 0, z] : SymMatrix 3) A
      = !![a, b, c; b, d, e; c, e, f] →
    0 < d1Lower →
    0 < d2Lower →
    0 < d3Lower →
    d1Lower ≤ a →
    d2Lower ≤ a * d - b * b →
    d3Lower ≤ a * d * f - a * e ^ 2 - b ^ 2 * f + 2 * b * c * e - c ^ 2 * d →
    A.PosDef

/-- **The preconditioner preserves the whole signature, not just definiteness.**

Stated separately because the runtime reads a signature rather than a yes/no: had
the block turned out indefinite, this is the theorem that would have carried
*that* result across instead. -/
def DiagonalCongruencePreservesIndexStatement : Prop :=
  ∀ (x y z : ℝ) (A : SymMatrix 3) (k : ℕ),
    x ≠ 0 → y ≠ 0 → z ≠ 0 →
    ((∃ V : Submodule ℝ (Fin 3 → ℝ), Module.finrank ℝ V = k ∧
        ∀ v ∈ V, v ≠ 0 → 0 < v ⬝ᵥ
          (congruence (!![x, 0, 0; 0, y, 0; 0, 0, z] : SymMatrix 3) A).mulVec v)
      ↔ (∃ V : Submodule ℝ (Fin 3 → ℝ), Module.finrank ℝ V = k ∧
        ∀ v ∈ V, v ≠ 0 → 0 < v ⬝ᵥ A.mulVec v))

/-- **And the rank, hence the zero count.** -/
def DiagonalCongruencePreservesRankStatement : Prop :=
  ∀ (x y z : ℝ) (A : SymMatrix 3),
    x ≠ 0 → y ≠ 0 → z ≠ 0 →
    (congruence (!![x, 0, 0; 0, y, 0; 0, 0, z] : SymMatrix 3) A).rank = A.rank

/-! ## The generalized gap (ENG-009) -/

/-- **Shifted positivity is a generalized Rayleigh bound.**

If `G - λM` is positive semidefinite, then `λ·vᵀMv ≤ vᵀGv` for every vector.
With `M` positive definite this says every generalized eigenvalue of the pencil
`(G, M)` is at least `λ` -- stated without eigenvalues, which is how the
runtime certifies it. This is the implication the ENG-009 generalized-gap
certificate rests on. -/
def GeneralizedRayleighStatement : Prop :=
  ∀ (n : ℕ) (G M : SymMatrix n) (lam : ℝ),
    (G - lam • M).PosSemidef →
    ∀ v : Fin n → ℝ, lam * (v ⬝ᵥ M.mulVec v) ≤ v ⬝ᵥ G.mulVec v

/-- **The certified gap does not depend on the coordinates.**

A shifted-definiteness certificate transports across any invertible change of
basis applied to both forms at once, in either direction, with the same `λ`.
This is the precise sense in which the generalized gap is basis-invariant while
raw eigenvalues and raw determinants are not -- the load-bearing claim of
ENG-009 §WO-RH-57. -/
def GeneralizedPencilCongruenceStatement : Prop :=
  ∀ (n : ℕ) (S G M : SymMatrix n) (lam : ℝ), IsUnit S.det →
    ((congruence S G - lam • congruence S M).PosDef ↔ (G - lam • M).PosDef)

/-- **The composition the ENG-009 runtime actually performs.**

Certified positive lower bounds on the three leading minors of the exactly
preconditioned *shifted* block `D(G - λM)D` imply the generalized Rayleigh
bound for the original pencil at that `λ`. Three prior theorems composed --
the minor criterion, the diagonal congruence, and the Rayleigh bound -- and
stated as one because the composition is what a reader of the gap certificate
has to trust. -/
def PreconditionedGapCertificate3Statement : Prop :=
  ∀ (x y z lam : ℝ) (G M : SymMatrix 3)
    (a b c d e f d1Lower d2Lower d3Lower : ℝ),
    x ≠ 0 → y ≠ 0 → z ≠ 0 →
    congruence (!![x, 0, 0; 0, y, 0; 0, 0, z] : SymMatrix 3) (G - lam • M)
      = !![a, b, c; b, d, e; c, e, f] →
    0 < d1Lower →
    0 < d2Lower →
    0 < d3Lower →
    d1Lower ≤ a →
    d2Lower ≤ a * d - b * b →
    d3Lower ≤ a * d * f - a * e ^ 2 - b ^ 2 * f + 2 * b * c * e - c ^ 2 * d →
    ∀ v : Fin 3 → ℝ, lam * (v ⬝ᵥ M.mulVec v) ≤ v ⬝ᵥ G.mulVec v

/-! ## The 4×4 block (ENG-010) -/

/-- **A rigorous 4×4 minor certificate implies positive definiteness.**

The minors are spelled out: `Δ₂`, `Δ₃`, `Δ₄` of the symmetric matrix
`[[a,b,c,d],[b,e,f,g],[c,f,h,i],[d,g,i,j]]`, each as its polynomial
expansion. -/
def PdFourByFourCertificateStatement : Prop :=
  ∀ (a b c d e f g h i j d1 d2 d3 d4 : ℝ),
    0 < d1 → 0 < d2 → 0 < d3 → 0 < d4 →
    d1 ≤ a →
    d2 ≤ a * e - b * b →
    d3 ≤ a * e * h - a * f ^ 2 - b ^ 2 * h + 2 * b * c * f - c ^ 2 * e →
    d4 ≤ a * e * h * j - a * e * i ^ 2 - a * f ^ 2 * j + 2 * a * f * g * i
        - a * g ^ 2 * h - b ^ 2 * h * j + b ^ 2 * i ^ 2 + 2 * b * c * f * j
        - 2 * b * c * g * i - 2 * b * d * f * i + 2 * b * d * g * h
        - c ^ 2 * e * j + c ^ 2 * g ^ 2 + 2 * c * d * e * i
        - 2 * c * d * f * g - d ^ 2 * e * h + d ^ 2 * f ^ 2 →
    (!![a, b, c, d; b, e, f, g; c, f, h, i; d, g, i, j] : SymMatrix 4).PosDef

/-- **The composition the ENG-010 runtime performs for the block:** certified
minor bounds on the exactly preconditioned block imply the original block is
positive definite. -/
def PreconditionedCertificate4Statement : Prop :=
  ∀ (w x y z : ℝ) (A : SymMatrix 4)
    (a b c d e f g h i j d1 d2 d3 d4 : ℝ),
    w ≠ 0 → x ≠ 0 → y ≠ 0 → z ≠ 0 →
    congruence (!![w, 0, 0, 0; 0, x, 0, 0; 0, 0, y, 0; 0, 0, 0, z] : SymMatrix 4) A
      = !![a, b, c, d; b, e, f, g; c, f, h, i; d, g, i, j] →
    0 < d1 → 0 < d2 → 0 < d3 → 0 < d4 →
    d1 ≤ a →
    d2 ≤ a * e - b * b →
    d3 ≤ a * e * h - a * f ^ 2 - b ^ 2 * h + 2 * b * c * f - c ^ 2 * e →
    d4 ≤ a * e * h * j - a * e * i ^ 2 - a * f ^ 2 * j + 2 * a * f * g * i
        - a * g ^ 2 * h - b ^ 2 * h * j + b ^ 2 * i ^ 2 + 2 * b * c * f * j
        - 2 * b * c * g * i - 2 * b * d * f * i + 2 * b * d * g * h
        - c ^ 2 * e * j + c ^ 2 * g ^ 2 + 2 * c * d * e * i
        - 2 * c * d * f * g - d ^ 2 * e * h + d ^ 2 * f ^ 2 →
    A.PosDef

/-- **The composed 4×4 gap implication** (ENG-010 §WO-RH-70): the same minor
bounds applied to the preconditioned *shifted* pencil `D(G − λM)D` force the
generalized Rayleigh bound for the original pencil at that `λ`. -/
def PreconditionedGapCertificate4Statement : Prop :=
  ∀ (w x y z lam : ℝ) (G M : SymMatrix 4)
    (a b c d e f g h i j d1 d2 d3 d4 : ℝ),
    w ≠ 0 → x ≠ 0 → y ≠ 0 → z ≠ 0 →
    congruence (!![w, 0, 0, 0; 0, x, 0, 0; 0, 0, y, 0; 0, 0, 0, z] : SymMatrix 4)
        (G - lam • M)
      = !![a, b, c, d; b, e, f, g; c, f, h, i; d, g, i, j] →
    0 < d1 → 0 < d2 → 0 < d3 → 0 < d4 →
    d1 ≤ a →
    d2 ≤ a * e - b * b →
    d3 ≤ a * e * h - a * f ^ 2 - b ^ 2 * h + 2 * b * c * f - c ^ 2 * e →
    d4 ≤ a * e * h * j - a * e * i ^ 2 - a * f ^ 2 * j + 2 * a * f * g * i
        - a * g ^ 2 * h - b ^ 2 * h * j + b ^ 2 * i ^ 2 + 2 * b * c * f * j
        - 2 * b * c * g * i - 2 * b * d * f * i + 2 * b * d * g * h
        - c ^ 2 * e * j + c ^ 2 * g ^ 2 + 2 * c * d * e * i
        - 2 * c * d * f * g - d ^ 2 * e * h + d ^ 2 * f ^ 2 →
    ∀ v : Fin 4 → ℝ, lam * (v ⬝ᵥ M.mulVec v) ≤ v ⬝ᵥ G.mulVec v

/-! ## The 5×5 block: shifted composition, nesting, Schur (ENG-011) -/

/-- **Shifted positivity composes back.** A certified shifted sweep, a
positive `λ`, and a positive definite reference metric force `G` itself
positive definite: `G = (G − λM) + λM`. -/
def ShiftedPositivityTransferStatement : Prop :=
  ∀ (n : ℕ) (G M : SymMatrix n) (lam : ℝ),
    0 < lam → M.PosDef → (G - lam • M).PosDef → G.PosDef

/-- **A smaller shift stays positive definite** — how a certified gap at a
larger `λ` licenses every smaller shift with no new computation. -/
def ShiftedShiftMonotoneStatement : Prop :=
  ∀ (n : ℕ) (G M : SymMatrix n) (a b : ℝ),
    b ≤ a → M.PosSemidef → (G - a • M).PosDef → (G - b • M).PosDef

/-- **Certified gap lower bounds regress upward through nesting.** A
shifted-PSD certificate restricts along any index map, so the certified `λ`
bounds every sub-pencil's Rayleigh quotients — the formal direction of
generalized Cauchy interlacing the ENG-011 runtime uses. -/
def NestedGapRegressionStatement : Prop :=
  ∀ (n m : ℕ) (G M : SymMatrix n) (lam : ℝ),
    (G - lam • M).PosSemidef →
    ∀ (e : Fin m → Fin n) (x : Fin m → ℝ),
      lam * (x ⬝ᵥ (M.submatrix e e).mulVec x)
        ≤ x ⬝ᵥ (G.submatrix e e).mulVec x

/-- **The Schur witness theorem.** `A ≻ 0`, an exact witness `A y = c`, and a
positive residual `d − c ⬝ᵥ y` make the bordered quadratic form
`xᵀAx + 2t(c ⬝ᵥ x) + d t²` positive away from zero — the semantic payload of
the ENG-011 verified-solve Schur analysis of the `b⁴` direction. -/
def SchurWitnessBlockStatement : Prop :=
  ∀ (n : ℕ) (A : SymMatrix n) (c : Fin n → ℝ) (d : ℝ) (y : Fin n → ℝ),
    A.PosDef → A.mulVec y = c → 0 < d - c ⬝ᵥ y →
    ∀ (x : Fin n → ℝ) (t : ℝ), x ≠ 0 ∨ t ≠ 0 →
      0 < x ⬝ᵥ A.mulVec x + 2 * t * (c ⬝ᵥ x) + d * t ^ 2

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
