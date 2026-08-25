/-
ATLAS-RH-ENG-007 §9 — certificate semantics. The most important product here.

This file answers one question: *given that the Python/Arb side produced honest
numerical enclosures, what finite mathematical conclusion follows?* It does not
model JSON, and it does not model Arb. It models the semantic payload of a
certificate -- the numbers and the inequalities among them -- and proves what
those inequalities imply.

The trust boundary, stated once and enforced by the shape of the definitions:

* **Python/Arb is responsible** for the enclosures being true: that the real
  entry `G i j` really does lie in the reported interval.
* **Lean is responsible** for the implication: that those intervals, if true,
  force the advertised conclusion.

Nothing here claims a verified interval integrator, and nothing here upgrades a
numerical result to a formal one. A certificate whose enclosures are wrong
yields a true implication about false premises, which is exactly the honest
division of labour.

No RH proof claim is made in this file.
-/
import AtlasRH.Definitions
import AtlasRH.Positivity
import AtlasRH.MatrixInertia

namespace AtlasRH

open Matrix

/-- An enclosure of a real number: the payload Arb hands over for one entry. -/
structure Enclosure where
  lo : ℝ
  hi : ℝ
  le : lo ≤ hi

/-- The proposition that an enclosure is honest about `x`. This is the premise
Python is responsible for and Lean never proves. -/
def Enclosure.Encloses (e : Enclosure) (x : ℝ) : Prop := e.lo ≤ x ∧ x ≤ e.hi

/-- The semantic payload of a 2×2 positive-definiteness certificate: enclosures
for the three distinct entries, together with the two positive lower bounds the
runtime reports.

This mirrors exactly what `e1_degree3_odd_positivity_log3_log4.json` carries:
a positive lower bound on the leading entry and a positive lower bound on the
determinant. -/
structure PD2Certificate where
  /-- Enclosure of the `(0,0)` entry. -/
  g00 : Enclosure
  /-- Enclosure of the off-diagonal entry. -/
  g0b : Enclosure
  /-- Enclosure of the `(1,1)` entry. -/
  gbb : Enclosure
  /-- Certified lower bound on the leading entry. -/
  g00Lower : ℝ
  /-- Certified lower bound on the determinant. -/
  detLower : ℝ
  h_g00_pos : 0 < g00Lower
  h_det_pos : 0 < detLower
  /-- The reported bound really is a lower bound for the enclosure. -/
  h_g00_bound : g00Lower ≤ g00.lo
  /-- The determinant bound holds against the worst corner of the entry
  enclosures: `g00` and `gbb` at their lower ends, and the off-diagonal at
  whichever end makes `g0b²` largest. Taking the maximum of the two squares is
  what makes this the worst case rather than an average one. -/
  h_det_bound : detLower ≤ g00.lo * gbb.lo - max (g0b.lo ^ 2) (g0b.hi ^ 2)

/-- A concrete matrix is *described by* a certificate when every entry lies in
its enclosure and the matrix is symmetric with the reported shape. -/
def PD2Certificate.Describes (c : PD2Certificate) (a b d : ℝ) : Prop :=
  c.g00.Encloses a ∧ c.g0b.Encloses b ∧ c.gbb.Encloses d

/-- A square is bounded by the larger of the endpoint squares, whichever side of
zero the interval sits on. This is the only real content in the determinant
step: an interval containing zero has its largest square at an endpoint too. -/
theorem sq_le_max_endpoints {lo hi x : ℝ} (h1 : lo ≤ x) (h2 : x ≤ hi) :
    x ^ 2 ≤ max (lo ^ 2) (hi ^ 2) := by
  rcases le_total 0 x with hx | hx
  · exact le_max_of_le_right (by nlinarith)
  · exact le_max_of_le_left (by nlinarith)

/-- **The certificate implication.**

If a certificate's enclosures are honest about a real symmetric 2×2 matrix, and
the certificate reports a positive lower bound on the leading entry and on the
determinant, then that matrix is positive definite.

This is the theorem that gives the ENG-006 degree-3 result its meaning: the Arb
run establishes the premises, and this establishes that they suffice. -/
theorem posDef_of_certificate (c : PD2Certificate) {a b d : ℝ}
    (h : c.Describes a b d) : (sym2 a b d).PosDef := by
  obtain ⟨⟨ha_lo, ha_hi⟩, ⟨hb_lo, hb_hi⟩, ⟨hd_lo, hd_hi⟩⟩ := h
  have hg00_lo_pos : 0 < c.g00.lo := lt_of_lt_of_le c.h_g00_pos c.h_g00_bound
  have ha : 0 < a := lt_of_lt_of_le hg00_lo_pos ha_lo
  refine posDef_sym2 ha ?_
  -- The off-diagonal square is bounded by the larger endpoint square, and the
  -- diagonal entries are at least their lower ends. Every one of those pushes
  -- the determinant down, so the reported bound survives the substitution.
  have hmax_nonneg : (0:ℝ) ≤ max (c.g0b.lo ^ 2) (c.g0b.hi ^ 2) :=
    le_trans (sq_nonneg _) (le_max_left _ _)
  have hb2 : b ^ 2 ≤ max (c.g0b.lo ^ 2) (c.g0b.hi ^ 2) := sq_le_max_endpoints hb_lo hb_hi
  have hprod_pos : 0 < c.g00.lo * c.gbb.lo := by
    have h1 := c.h_det_bound
    have h2 := c.h_det_pos
    linarith
  have hgbb_lo_pos : 0 < c.gbb.lo := by
    by_contra hcon
    rw [not_lt] at hcon
    nlinarith [hg00_lo_pos, hprod_pos]
  have hprod : c.g00.lo * c.gbb.lo ≤ a * d :=
    mul_le_mul ha_lo hd_lo (le_of_lt hgbb_lo_pos) (le_of_lt ha)
  have hbb : b * b = b ^ 2 := by ring
  rw [hbb]
  linarith [c.h_det_pos, c.h_det_bound, hprod, hb2]

/-- The signature the certificate licenses: a certified 2×2 block has inertia
`(2, 0, 0)`. This is the exact statement the ENG-006 artifact makes. -/
theorem definiteInertia_of_certificate (c : PD2Certificate) {a b d : ℝ}
    (h : c.Describes a b d) : HasDefiniteInertia (sym2 a b d) :=
  posDef_of_certificate c h

/-! ### Inertia transcripts

The runtime's LDL engine reports a *transcript*: a sequence of pivots with
determined signs. For the definite case the transcript's content is exactly that
every pivot was positive, and the congruence theorems say the signature is then
the matrix's own. -/

/-! ## The 3×3 certificate (ENG-008 §WO-RH-53)

The 2×2 certificate above carries entry enclosures and lets the determinant
bound be derived from them, because at 2×2 the worst corner of a determinant is
findable in closed form. At 3×3 it is not, and encoding interval arithmetic in
Lean would be the wrong division of labour anyway: §WO-RH-53 is explicit that
Lean verifies the *implication*, not the Arb arithmetic.

So the 3×3 certificate's semantic payload is stated where the runtime actually
establishes it -- as certified positive lower bounds on the three leading
principal minors. Arb proves those hold; what follows is that they suffice. -/

/-- The semantic payload of a 3×3 positive-definiteness certificate: a positive
certified lower bound for each leading principal minor. -/
structure PD3Certificate where
  /-- Certified lower bound on `Δ₁ = a`. -/
  d1Lower : ℝ
  /-- Certified lower bound on `Δ₂ = ad − b²`. -/
  d2Lower : ℝ
  /-- Certified lower bound on `Δ₃ = det`. -/
  d3Lower : ℝ
  h_d1_pos : 0 < d1Lower
  h_d2_pos : 0 < d2Lower
  h_d3_pos : 0 < d3Lower

/-- A concrete symmetric 3×3 is *described by* the certificate when each of its
leading principal minors is at least the corresponding certified bound. -/
def PD3Certificate.Describes (c : PD3Certificate) (a b cc d e f : ℝ) : Prop :=
  c.d1Lower ≤ a ∧
  c.d2Lower ≤ minor2 a b cc d e f ∧
  c.d3Lower ≤ minor3 a b cc d e f

/-- **The 3×3 certificate implication.**

Certified positive lower bounds on the three leading principal minors imply the
block is positive definite. This is the theorem the ENG-008 degree-4 result
rests on: the interval run establishes the three bounds, and this establishes
that they suffice. -/
theorem posDef_of_certificate3 (c : PD3Certificate) {a b cc d e f : ℝ}
    (h : c.Describes a b cc d e f) : (sym3 a b cc d e f).PosDef := by
  obtain ⟨h1, h2, h3⟩ := h
  exact posDef_sym3 (lt_of_lt_of_le c.h_d1_pos h1)
    (lt_of_lt_of_le c.h_d2_pos h2) (lt_of_lt_of_le c.h_d3_pos h3)

/-- And the signature it licenses. -/
theorem definiteInertia_of_certificate3 (c : PD3Certificate) {a b cc d e f : ℝ}
    (h : c.Describes a b cc d e f) : HasDefiniteInertia (sym3 a b cc d e f) :=
  posDef_of_certificate3 c h

/-! ### The certificate the runtime actually produces

The bounds are certified for the *preconditioned* block `DᵀGD`, not for `G`.
Composing the two steps -- the certificate implication and the diagonal
congruence -- is the whole chain from what Arb proved to what the certificate
claims, and it is worth having as one theorem rather than two, because the
composition is the thing a reader has to trust. -/

/-- **Preconditioned certificate ⟹ the original block is positive definite.**

Given a diagonal preconditioner with nonzero entries, and certified positive
lower bounds on the three leading minors of the rescaled block, the *original*
block is positive definite. -/
theorem posDef_of_preconditioned_certificate3
    (c : PD3Certificate) {x y z : ℝ} {A : SymMatrix 3} {a b cc d e f : ℝ}
    (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (hcong : congruence (diag3 x y z) A = sym3 a b cc d e f)
    (h : c.Describes a b cc d e f) : A.PosDef := by
  refine posDef_of_diagonal_congruence hx hy hz ?_
  rw [hcong]
  exact posDef_of_certificate3 c h

/-- A congruence transcript witnesses that `A` is congruent to `D` by an
invertible `S` -- the record `inertia/ldl.py` produces as it eliminates. -/
structure CongruenceTranscript (n : ℕ) where
  source : SymMatrix n
  target : SymMatrix n
  transform : SymMatrix n
  h_unit : IsUnit transform.det
  h_congr : target = congruence transform source

/-- **A transcript transfers definiteness.** If the eliminated form is positive
definite then so is the original. This is the formal content of "the pivot
signature is the matrix signature": what the engine observes after elimination
is a fact about what it started with. -/
theorem posDef_of_transcript {n : ℕ} (t : CongruenceTranscript n)
    (h : t.target.PosDef) : t.source.PosDef := by
  rw [t.h_congr] at h
  exact posDef_of_congruence t.h_unit h

/-- And conversely, so a transcript can be run in either direction. -/
theorem transcript_posDef {n : ℕ} (t : CongruenceTranscript n)
    (h : t.source.PosDef) : t.target.PosDef := by
  rw [t.h_congr]
  exact posDef_congruence t.h_unit h

end AtlasRH
