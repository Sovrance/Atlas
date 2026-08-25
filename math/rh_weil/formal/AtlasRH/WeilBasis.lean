/-
ATLAS-RH-ENG-007 §8 — exact RH/Weil basis and parity identities.

These are the algebraic facts the block decomposition rests on. The Gram matrix
is treated as block diagonal by parity about `x = L/2`, which is what lets the
odd degree-3 block be analysed as a standalone 2×2 rather than as a corner of a
4×4. If the parity classification were wrong, that reduction would be wrong, and
every certificate built on it would be internally consistent and meaningless.

§8 is explicit that the digamma integral and Arb arithmetic are out of scope
here. What is in scope is the finite algebra: which basis elements are even,
which are odd, and why an even/odd cross entry of a reflection-invariant form
vanishes.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity

namespace AtlasRH

section Parity

variable (L : ℝ)

/-- Reflection is an involution on `[0, L]`. -/
@[simp] theorem reflect_reflect (x : ℝ) : reflect L (reflect L x) = x := by
  simp [reflect]

/-- The constant is even about the midpoint. -/
theorem basisOne_even : EvenAbout L basisOne := fun _ => rfl

/-- The bubble `x(L - x)` is even about the midpoint. -/
theorem basisB_even : EvenAbout L (basisB L) := by
  intro x; simp only [basisB, reflect]; ring

/-- `q1(x) = x - L/2` is odd about the midpoint. -/
theorem basisQ1_odd : OddAbout L (basisQ1 L) := by
  intro x; simp only [basisQ1, reflect]; ring

/-- `b3(x) = x(L - x)(x - L/2)` is odd about the midpoint.

`b3` is the product of the even bubble and the odd `q1`, so its parity is forced;
this is why `{q1, b3}` is a parity-matched pair and forms the odd block. -/
theorem basisB3_odd : OddAbout L (basisB3 L) := by
  intro x; simp only [basisB3, reflect]; ring

/-- `b3 = b * q1` pointwise: the odd degree-3 element is the even bubble times
the odd linear element. -/
theorem basisB3_eq_mul (x : ℝ) : basisB3 L x = basisB L x * basisQ1 L x := by
  simp only [basisB3, basisB, basisQ1]

/-- An even function times an odd function is odd. This is the general reason
`b3` is odd, stated once rather than rediscovered per basis element. -/
theorem odd_of_even_mul_odd {f g : ℝ → ℝ} (hf : EvenAbout L f) (hg : OddAbout L g) :
    OddAbout L (fun x => f x * g x) := by
  intro x
  simp only [hf x, hg x]
  ring

/-- A product of two odd functions is even. -/
theorem even_of_odd_mul_odd {f g : ℝ → ℝ} (hf : OddAbout L f) (hg : OddAbout L g) :
    EvenAbout L (fun x => f x * g x) := by
  intro x
  simp only [hf x, hg x]
  ring

end Parity

section CrossBlock

variable {L : ℝ}

/-- **Even/odd cross entries vanish.**

A bilinear form built from a reflection-invariant pairing sends an even and an
odd function to something equal to its own negation, hence to zero. This is the
statement that makes the Gram matrix block diagonal by parity.

The hypothesis is exactly what the Weil pairing supplies: invariance of the
pairing under simultaneous reflection of both arguments. Nothing about digamma,
primes or integration enters -- only that one argument flips sign and the other
does not. -/
theorem cross_block_vanishes
    (pair : (ℝ → ℝ) → (ℝ → ℝ) → ℝ)
    (bilin_neg : ∀ f g, pair f (fun x => -g x) = -pair f g)
    (refl_inv : ∀ f g, pair (fun x => f (reflect L x)) (fun x => g (reflect L x)) = pair f g)
    {f g : ℝ → ℝ} (hf : EvenAbout L f) (hg : OddAbout L g) :
    pair f g = 0 := by
  have h1 : pair (fun x => f (reflect L x)) (fun x => g (reflect L x)) = pair f g := refl_inv f g
  have hfe : (fun x => f (reflect L x)) = f := funext hf
  have hgo : (fun x => g (reflect L x)) = (fun x => -g x) := funext hg
  rw [hfe, hgo, bilin_neg] at h1
  linarith

end CrossBlock

section DegreeTwo

/-- The 2×2 determinant identity used throughout: `E2 = G00·Gbb - G0b²`. -/
theorem det_two_identity (g00 g0b gbb : ℝ) :
    (sym2 g00 g0b gbb).det = g00 * gbb - g0b ^ 2 := by
  rw [det_sym2]; ring

/-- **Degree-2 parity factorization.** With the basis ordered `{1, b, q1, b3}`
and the cross entries vanishing by parity, the determinant of the whole
degree-3 Gram is the product of the even and odd block determinants. Stated on
the block data itself, since that is the form the runtime produces. -/
theorem parity_block_det (e00 e0b ebb o11 o1b obb : ℝ) :
    (sym2 e00 e0b ebb).det * (sym2 o11 o1b obb).det
      = (e00 * ebb - e0b ^ 2) * (o11 * obb - o1b ^ 2) := by
  rw [det_sym2, det_sym2]; ring

end DegreeTwo

end AtlasRH
