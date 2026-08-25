/-
ATLAS-RH-ENG-007 §5 Layer B — the comparator's view of the trusted statements.

The statement text itself lives in `AtlasRH/Statements.lean`, which imports only
`AtlasRH.Definitions` and Mathlib. This module is the comparator-facing alias
layer: it names, in one place, exactly the propositions the manifest tracks and
`comparator/Solution.lean` must inhabit.

Keeping the text in the library rather than here avoids the failure mode this
work order exists to prevent. Two independently maintained copies of a
proposition drift; one copy with a hashed source region does not.
-/
import AtlasRH.Statements

namespace Comparator

/-- Sylvester's law of inertia, positive index. -/
abbrev inertia_congruence_positive_statement : Prop :=
  AtlasRH.InertiaCongruencePositiveStatement

/-- Sylvester's law of inertia, negative index. -/
abbrev inertia_congruence_negative_statement : Prop :=
  AtlasRH.InertiaCongruenceNegativeStatement

/-- Sylvester's law of inertia, rank (hence the zero count). -/
abbrev inertia_congruence_rank_statement : Prop :=
  AtlasRH.InertiaCongruenceRankStatement

/-- The 2×2 leading-principal-minor criterion. -/
abbrev pd_two_by_two_statement : Prop := AtlasRH.PdTwoByTwoStatement

/-- The 3×3 leading-principal-minor criterion. -/
abbrev pd_three_by_three_statement : Prop := AtlasRH.PdThreeByThreeStatement

/-- A rigorous 2×2 enclosure certificate implies positive definiteness. -/
abbrev certificate_even2_implies_pd_statement : Prop :=
  AtlasRH.CertificateEven2ImpliesPdStatement

/-- Parity of the Weil basis about the midpoint. -/
abbrev weil_basis_parity_statement : Prop := AtlasRH.WeilBasisParityStatement

/-- Cross-parity blocks of a reflection-invariant pairing vanish. -/
abbrev odd_degree3_cross_block_statement : Prop :=
  AtlasRH.OddDegree3CrossBlockStatement

/-- Degree-2 parity factorization of the Gram determinant. -/
abbrev odd_degree3_factorization_statement : Prop :=
  AtlasRH.OddDegree3FactorizationStatement

/-- Certified positive leading minors imply a 3×3 block is positive definite. -/
abbrev pd_three_by_three_certificate_statement : Prop :=
  AtlasRH.PdThreeByThreeCertificateStatement

/-- A nonzero diagonal preconditioner preserves positive definiteness. -/
abbrev diagonal_congruence_preserves_pd_statement : Prop :=
  AtlasRH.DiagonalCongruencePreservesPdStatement

/-- The composition: certified minors of the rescaled block imply the original
block is positive definite. -/
abbrev preconditioned_certificate3_statement : Prop :=
  AtlasRH.PreconditionedCertificate3Statement

/-- The preconditioner preserves the positive index. -/
abbrev diagonal_congruence_preserves_index_statement : Prop :=
  AtlasRH.DiagonalCongruencePreservesIndexStatement

/-- The preconditioner preserves the rank. -/
abbrev diagonal_congruence_preserves_rank_statement : Prop :=
  AtlasRH.DiagonalCongruencePreservesRankStatement

/-- The rank–trace inequality in the `Q = 0` case the runtime uses. -/
abbrev rank_trace_hs_statement : Prop := AtlasRH.RankTraceZeroQStatement

end Comparator
