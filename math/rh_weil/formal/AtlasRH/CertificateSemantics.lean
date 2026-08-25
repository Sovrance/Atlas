/-
ATLAS-RH-ENG-007 §9 (WO-RH-42) — what a rigorous finite certificate *means*.

This is the product ENG-007 exists for. It is not a JSON parser and not a verified
integrator. It is the implication layer, and the trust boundary is deliberately explicit:

  * Python/Arb is responsible for producing TRUE numerical enclosures. Nothing here
    checks that, and nothing here could.
  * Lean is responsible for proving that those enclosures, if true, imply the finite
    mathematical conclusion the certificate advertises.

The failure this prevents is the one ENG-007 §0 names: the code proving one finite theorem
while a certificate consumer believes it proved another. A consumer that reads
`e1_degree2_compact_log3_log4.json` and concludes "the Gram block is positive definite on
the cell" is now relying on a machine-checked implication rather than on someone's memory
of what the number meant.

No RH proof claim. Everything below is finite-dimensional.
-/
import AtlasRH.Positivity

namespace AtlasRH

open Matrix

/-- The semantic payload of a rigorous 2x2 positive-definiteness certificate.

Mirrors what the Atlas degree-2/degree-3 certificates actually report: a rigorous lower
bound on the leading entry and a rigorous lower bound on the determinant, both strictly
positive. The `h_*` fields are the *claims*; producing enclosures that justify them is the
interval engine's job, not Lean's. -/
structure PD2Certificate where
  /-- Rigorous lower bound on `G₀₀`, from the scalar canary. -/
  g00Lower : ℝ
  /-- Rigorous lower bound on `det G = G₀₀G₁₁ - G₀₁²`, e.g. `E2`'s certified bound. -/
  detLower : ℝ
  h_g00 : 0 < g00Lower
  h_det : 0 < detLower

namespace PD2Certificate

variable (cert : PD2Certificate)

/-- What the interval engine is asserting when it releases the certificate over a cell:
at every parameter in the domain the block is real symmetric with entries enclosed so that
the two reported lower bounds hold.

This proposition is the *premise*. Arb discharges it numerically; Lean consumes it. -/
def Covers (cert : PD2Certificate) (G : ℝ → Matrix (Fin 2) (Fin 2) ℝ) (D : Set ℝ) : Prop :=
  ∀ L ∈ D, ∃ a b c : ℝ,
    G L = sym2 a b c ∧ cert.g00Lower ≤ a ∧ cert.detLower ≤ a * c - b * b

/-- **The certificate implication.**

A valid 2x2 positivity certificate over a domain implies the block is positive definite at
every point of that domain.

This is the theorem a consumer of `e1_degree2_compact_log3_log4.json` is implicitly using.
Note what it does *not* say: it does not make the enclosures rigorous, and it does not
upgrade the numerical evidence class. It converts a numerical premise into a mathematical
conclusion, and that conversion is now machine-checked. -/
theorem posDef_of_covers {G : ℝ → Matrix (Fin 2) (Fin 2) ℝ} {D : Set ℝ}
    (h : cert.Covers G D) : ∀ L ∈ D, (G L).PosDef := by
  intro L hL
  obtain ⟨a, b, c, hG, ha, hdet⟩ := h L hL
  rw [hG]
  exact sym2_posDef_of a b c (lt_of_lt_of_le cert.h_g00 ha)
    (lt_of_lt_of_le cert.h_det hdet)

/-- The same conclusion stated as the finite signature the inertia channel reports.

ENG-006 reports the odd degree-3 block as inertia `(2,0,0)`. For a 2x2 real symmetric
block, "positive definite" and "inertia `(2,0,0)`" are the same claim; stating it in the
inertia vocabulary keeps the certificate consumer honest about which channel it read.

The distinction ENG-006 insists on is preserved: an inertia certificate reporting, say,
`(1,0,1)` would satisfy this statement's shape but not its hypotheses, and would never
license a PSD consumer. -/
theorem inertia_two_zero_zero_of_covers {G : ℝ → Matrix (Fin 2) (Fin 2) ℝ} {D : Set ℝ}
    (h : cert.Covers G D) (L : ℝ) (hL : L ∈ D) :
    (G L).PosDef ∧ (Inertia.mk 2 0 0).pos = 2 :=
  ⟨cert.posDef_of_covers h L hL, rfl⟩

end PD2Certificate

/-- A pivot transcript for the 2x2 congruence the runtime performs: the two leading
pivots it observed, both required positive.

Sylvester's law of inertia is what makes reading signs off a congruent block sound, and
`congruence_posDef_iff` is the formal version of that soundness. -/
structure PivotTranscript where
  pivot0 : ℝ
  pivot1 : ℝ
  h_pivot0 : 0 < pivot0
  h_pivot1 : 0 < pivot1

/-- A positive pivot transcript for `!![a,b;b,c]`, where the pivots are `a` and the Schur
complement `c - b²/a`, implies positive definiteness.

This is the exact-arithmetic path the B1 solver and the inertia engine share: the Schur
pivot `S = c - b²/a` generalises the `G₀₀ > 0`, `det > 0` verifier standard, and the
runtime reports its sign. -/
theorem posDef_of_schur_pivots (a b c : ℝ) (ha : 0 < a) (hschur : 0 < c - b * b / a) :
    (sym2 a b c).PosDef := by
  refine sym2_posDef_of a b c ha ?_
  have h : a * c - b * b = a * (c - b * b / a) := by field_simp
  rw [h]
  exact mul_pos ha hschur

end AtlasRH
