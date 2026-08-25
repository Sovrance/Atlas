/-
ATLAS-RH-ENG-007 §10 — axiom audit and statement emission.

Run: `lake env lean comparator/Comparator/PrintAxioms.lean`

§10 requires that allowed standard axioms be "enumerated, not hidden". The point is not an
empty list -- `propext`, `Classical.choice` and `Quot.sound` are unavoidable with Mathlib --
but that the list is visible, so a project-specific `axiom` or a `sorry` (which appears as
`sorryAx`) shows up as a diff instead of hiding behind a green build.

`scripts/check_formal_manifest.py` parses this output, rebuilds the manifest, and fails on
any axiom outside the allowlist or any statement that has drifted from the committed hash.
-/
import Comparator.Solution

-- Statements, as elaborated propositions. These are what get hashed into the manifest.
#print Comparator.TrustedStatements.congruence_preserves_posDef
#print Comparator.TrustedStatements.pd_two_by_two
#print Comparator.TrustedStatements.schur_pivot_implies_posDef
#print Comparator.TrustedStatements.certificate_even2_implies_pd

-- Axiom dependencies of every exported theorem.
#print axioms Comparator.Solution.congruence_preserves_posDef
#print axioms Comparator.Solution.pd_two_by_two
#print axioms Comparator.Solution.schur_pivot_implies_posDef
#print axioms Comparator.Solution.certificate_even2_implies_pd
