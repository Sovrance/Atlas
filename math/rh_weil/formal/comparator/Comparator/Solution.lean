/-
ATLAS-RH-ENG-007 §5 Layer C/D — the proofs, and the statement comparison.

Each declaration below has its type ascribed to a `TrustedStatements` definition and is
discharged by an `AtlasRH` theorem. That ascription *is* the comparator: if an AtlasRH
theorem ever drifts from the statement Atlas claims to rely on, this file stops
typechecking. No separate diffing tool is needed, and none could be more trustworthy than
the kernel.

ENG-007 §5: "The point is not byte-level source equality. The point is to prevent a future
proof of a subtly changed theorem from satisfying the original certificate consumer."
-/
import AtlasRH
import Comparator.TrustedStatements

namespace Comparator.Solution

open Matrix

/-- Discharges `congruence_preserves_posDef`. -/
theorem congruence_preserves_posDef : TrustedStatements.congruence_preserves_posDef := by
  intro n _ _ S A hS
  exact AtlasRH.congruence_posDef_iff (S := S) (A := A) hS

/-- Discharges `pd_two_by_two`. -/
theorem pd_two_by_two : TrustedStatements.pd_two_by_two := by
  intro a b c
  exact AtlasRH.sym2_posDef_iff a b c

/-- Discharges `schur_pivot_implies_posDef`. -/
theorem schur_pivot_implies_posDef : TrustedStatements.schur_pivot_implies_posDef := by
  intro a b c ha hs
  exact AtlasRH.posDef_of_schur_pivots a b c ha hs

/-- Discharges `certificate_even2_implies_pd`.

Built by packaging the raw bounds into `PD2Certificate` and applying the implication
theorem, which is how a real consumer uses it. -/
theorem certificate_even2_implies_pd : TrustedStatements.certificate_even2_implies_pd := by
  intro g00Lower detLower hg hd G D hcov L hL
  exact AtlasRH.PD2Certificate.posDef_of_covers
    { g00Lower := g00Lower, detLower := detLower, h_g00 := hg, h_det := hd }
    hcov L hL

end Comparator.Solution
