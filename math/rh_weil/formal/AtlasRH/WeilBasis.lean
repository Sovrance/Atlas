/-
ATLAS-RH-ENG-007 §8 (WO-RH-41) — exact RH/Weil basis and parity identities.

Everything here is exact algebra over the polynomial basis Atlas uses on a cell `[0, L]`:

    1                      even
    b  = x(L-x)            even
    q1 = x - L/2           odd
    b3 = x(L-x)(x - L/2)   odd

"Even" and "odd" are with respect to the **midpoint reflection** `σ(x) = L - x`. The whole
point of the basis is that this reflection block-diagonalises the Gram matrix, which is what
turns a 3x3 determinant into `O1 · E2` and lets the runtime read the two parity sectors
independently.

§8 says not to formalize Arb or the digamma integral here, and this file does not: it
formalizes the *algebraic* content -- the parity of the basis, the vanishing of even/odd
cross terms under any reflection-invariant form, and the determinant factorization that
follows. The analytic construction of the form stays outside Lean.
-/
import AtlasRH.Definitions

namespace AtlasRH.WeilBasis

open Matrix

/-! ### Midpoint reflection and basis parity -/

/-- Midpoint reflection on the cell `[0, L]`. -/
def reflect (L x : ℝ) : ℝ := L - x

/-- The bubble basis element `b(x) = x(L-x)`. -/
def bubble (L x : ℝ) : ℝ := x * (L - x)

/-- The midpoint-odd element `q1(x) = x - L/2`. -/
noncomputable def q1 (L x : ℝ) : ℝ := x - L / 2

/-- The odd degree-3 element `b3(x) = x(L-x)(x - L/2)`. -/
noncomputable def b3 (L x : ℝ) : ℝ := x * (L - x) * (x - L / 2)

/-- The constant `1` is reflection-even (trivially). -/
theorem one_even (L x : ℝ) : (1 : ℝ) = 1 := rfl

/-- `b` is reflection-**even**: `b(L - x) = b(x)`. -/
theorem bubble_even (L x : ℝ) : bubble L (reflect L x) = bubble L x := by
  simp only [bubble, reflect]; ring

/-- `q1` is reflection-**odd**: `q1(L - x) = -q1(x)`. -/
theorem q1_odd (L x : ℝ) : q1 L (reflect L x) = -q1 L x := by
  simp only [q1, reflect]; ring

/-- `b3` is reflection-**odd**: `b3(L - x) = -b3(x)`.

This is the parity that puts the degree-3 element in the same sector as `q1`, which is why
ENG-006's degree-3 block is the odd pair `{q1, b3}` and not a mixed block. -/
theorem b3_odd (L x : ℝ) : b3 L (reflect L x) = -b3 L x := by
  simp only [b3, reflect]; ring

/-- Reflection is an involution. -/
theorem reflect_involutive (L x : ℝ) : reflect L (reflect L x) = x := by
  simp only [reflect]; ring

/-! ### Even/odd cross terms vanish

This is the structural fact the basis is chosen for, and it needs no analysis: it holds for
*any* bilinear form invariant under the reflection. -/

/-- **Even/odd orthogonality.**

If a bilinear form `B` is invariant under an involution `σ`, then any `σ`-even vector is
`B`-orthogonal to any `σ`-odd vector.

The proof is one line of algebra and is the entire reason the Gram matrix block-diagonalises:
`B(f,g) = B(σf, σg) = B(f, -g) = -B(f,g)`, so `2·B(f,g) = 0`.

Stated abstractly over a real vector space so that it applies to whichever concrete form the
runtime assembles -- `G = G⁰ - Gᵖ + G^∞` -- without Lean needing to know how that form is
built. What must be checked outside Lean is the *invariance hypothesis*; given it, the
vanishing is automatic. -/
theorem even_odd_orthogonal {V : Type*} [AddCommGroup V] [Module ℝ V]
    (B : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (σ : V →ₗ[ℝ] V)
    (hinv : ∀ u v : V, B (σ u) (σ v) = B u v)
    (f g : V) (hf : σ f = f) (hg : σ g = -g) :
    B f g = 0 := by
  have h1 : B (σ f) (σ g) = B f g := hinv f g
  rw [hf, hg] at h1
  simp only [map_neg] at h1
  linarith [h1]

/-! ### Determinant identities

The two the runtime checks by name are `det_matches_O1_times_E2` and
`D2_matches_E2_plus_L2_G00_O1`. -/

/-- The 2x2 determinant identity, in the form the even sector reports it:
`E2 = G₀₀·G_bb − G₀ᵦ²`. -/
theorem even_block_det (g00 g0b gbb : ℝ) :
    (sym2 g00 g0b gbb).det = g00 * gbb - g0b * g0b := sym2_det g00 g0b gbb

/-- **Parity factorization of the degree-≤2 Gram determinant.**

Order the basis `(1, q1, b)` so that indices `0, 2` are even and index `1` is odd. Even/odd
orthogonality forces the cross entries `G₀₁`, `G₁₂` to vanish, and the determinant then
factors as

    det G = O1 · E2,   O1 = G₁₁,   E2 = G₀₀·G₂₂ − G₀₂².

This is exactly the `det_matches_O1_times_E2` invariant the E1 certifier asserts, and it is
why the odd pivot can be read off independently of the even block. -/
theorem det_parity_factorization (g00 g02 g22 o1 : ℝ) :
    (!![g00, 0, g02; 0, o1, 0; g02, 0, g22] : Matrix (Fin 3) (Fin 3) ℝ).det
      = o1 * (g00 * g22 - g02 * g02) := by
  simp [Matrix.det_fin_three]
  ring

/-- **Basis congruence preserves the Gram determinant.**

The parity basis `(1, q1, b)` and the monomial basis `(1, x, x²)` are related by a
unitriangular change of basis, and a congruence by any determinant-one matrix leaves the
Gram determinant unchanged:

    det (Sᵀ G S) = (det S)² · det G = det G   when   det S = 1.

This is the §8 "basis-congruence relationship used by the finite Gram blocks", and it is
what licenses the runtime to compute a determinant in whichever basis is convenient and
compare it against the parity-split product.

Note what is *not* claimed here. The runtime's `D2 = E2 + L²·G₀₀·O1` is, in
`core.degree2_raw_det`, the definition of `D2` rather than a theorem about it; the content
of the `D2_matches_E2_plus_L2_G00_O1` check is that an independently computed monomial-basis
determinant agrees with it numerically. That cross-check remains an E1 runtime assertion.
Formalizing it would require carrying the explicit basis-change matrix, which §8 does not
ask for and which is not attempted here rather than being faked with a definitional
identity. -/
theorem det_congruence_invariant {n : Type*} [Fintype n] [DecidableEq n]
    (S G : Matrix n n ℝ) (hS : S.det = 1) :
    (Sᵀ * G * S).det = G.det := by
  rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hS]
  ring

end AtlasRH.WeilBasis
