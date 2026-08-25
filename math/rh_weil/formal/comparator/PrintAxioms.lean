/-
ATLAS-RH-ENG-007 §10 (WO-RH-43) — axiom audit and statement comparator report.

Run with:

    lake env lean comparator/PrintAxioms.lean

For every theorem `comparator/Solution.lean` exports, this emits one
tab-separated line:

    ATLAS_FORMAL_THEOREM <id> <solution theorem> <trusted statement> <axioms> <normalized statement>

`scripts/check_formal_manifest.py` consumes those lines. It checks that the
axiom set is a subset of the three standard Lean axioms, that no
project-specific `axiom` declaration appears, and that the normalized
pretty-printed statement still hashes to the value recorded in
`manifests/theorem_manifest.json`.

The hashed statement is the elaborated body of the trusted `def` in
`AtlasRH/Statements.lean`, pretty-printed by Lean rather than read off the
source, so a change of notation that means the same thing does not trip the
gate and a change of meaning does. The link from that body to the theorem is
checked separately, by `isDefEq` against the solution theorem's own type.
-/
import comparator.Solution

open Lean Elab Command Meta

/-- The theorems under audit, paired with the trusted statement each inhabits.
Adding a theorem here without adding it to the manifest fails CI, and vice
versa. -/
def atlasAudited : List (String × Name × Name) :=
  [ ("inertia_congruence_positive",
      ``Comparator.atlas_inertia_congruence_positive,
      ``AtlasRH.InertiaCongruencePositiveStatement)
  , ("inertia_congruence_negative",
      ``Comparator.atlas_inertia_congruence_negative,
      ``AtlasRH.InertiaCongruenceNegativeStatement)
  , ("inertia_congruence_rank",
      ``Comparator.atlas_inertia_congruence_rank,
      ``AtlasRH.InertiaCongruenceRankStatement)
  , ("pd_two_by_two",
      ``Comparator.atlas_pd_two_by_two,
      ``AtlasRH.PdTwoByTwoStatement)
  , ("pd_three_by_three",
      ``Comparator.atlas_pd_three_by_three,
      ``AtlasRH.PdThreeByThreeStatement)
  , ("certificate_even2_implies_pd",
      ``Comparator.atlas_certificate_even2_implies_pd,
      ``AtlasRH.CertificateEven2ImpliesPdStatement)
  , ("weil_basis_parity",
      ``Comparator.atlas_weil_basis_parity,
      ``AtlasRH.WeilBasisParityStatement)
  , ("odd_degree3_cross_block",
      ``Comparator.atlas_odd_degree3_cross_block,
      ``AtlasRH.OddDegree3CrossBlockStatement)
  , ("odd_degree3_factorization",
      ``Comparator.atlas_odd_degree3_factorization,
      ``AtlasRH.OddDegree3FactorizationStatement)
  , ("pd_three_by_three_certificate",
      ``Comparator.atlas_pd_three_by_three_certificate,
      ``AtlasRH.PdThreeByThreeCertificateStatement)
  , ("diagonal_congruence_preserves_pd",
      ``Comparator.atlas_diagonal_congruence_preserves_pd,
      ``AtlasRH.DiagonalCongruencePreservesPdStatement)
  , ("preconditioned_certificate3",
      ``Comparator.atlas_preconditioned_certificate3,
      ``AtlasRH.PreconditionedCertificate3Statement)
  , ("diagonal_congruence_preserves_index",
      ``Comparator.atlas_diagonal_congruence_preserves_index,
      ``AtlasRH.DiagonalCongruencePreservesIndexStatement)
  , ("diagonal_congruence_preserves_rank",
      ``Comparator.atlas_diagonal_congruence_preserves_rank,
      ``AtlasRH.DiagonalCongruencePreservesRankStatement)
  , ("rank_trace_hs",
      ``Comparator.atlas_rank_trace_hs,
      ``AtlasRH.RankTraceZeroQStatement)
  , ("generalized_rayleigh",
      ``Comparator.atlas_generalized_rayleigh,
      ``AtlasRH.GeneralizedRayleighStatement)
  , ("generalized_pencil_congruence",
      ``Comparator.atlas_generalized_pencil_congruence,
      ``AtlasRH.GeneralizedPencilCongruenceStatement)
  , ("preconditioned_gap_certificate3",
      ``Comparator.atlas_preconditioned_gap_certificate3,
      ``AtlasRH.PreconditionedGapCertificate3Statement)
  , ("pd_four_by_four_certificate",
      ``Comparator.atlas_pd_four_by_four_certificate,
      ``AtlasRH.PdFourByFourCertificateStatement)
  , ("preconditioned_certificate4",
      ``Comparator.atlas_preconditioned_certificate4,
      ``AtlasRH.PreconditionedCertificate4Statement)
  , ("preconditioned_gap_certificate4",
      ``Comparator.atlas_preconditioned_gap_certificate4,
      ``AtlasRH.PreconditionedGapCertificate4Statement) ]

run_cmd do
  let env ← getEnv
  for (id, solution, trusted) in atlasAudited do
    let some solInfo := env.find? solution
      | throwError "missing solution theorem {solution}"
    let some trustedInfo := env.find? trusted
      | throwError "missing trusted statement {trusted}"
    -- The Layer D comparison. `Solution` declares the theorem at the aliased
    -- type, so Lean has already checked it; this re-derives the unfolded form
    -- that gets hashed, and fails loudly if the two ever come apart.
    let solType := solInfo.type
    let trustedBody := trustedInfo.value?.getD (mkConst trusted)
    let same ← liftTermElabM <| withoutModifyingEnv <| Meta.isDefEq solType trustedBody
    unless same do
      throwError "statement comparator FAILED for {id}: {solution} does not inhabit {trusted}"
    let axs ← collectAxioms solution
    let axStr := if axs.isEmpty then "none" else
      String.intercalate "," (axs.toList.map toString)
    let stmt ← liftTermElabM <| withOptions (fun o => o.setBool `pp.unicode.fun true) <|
      do return (← Meta.ppExpr trustedBody).pretty (width := 1000000)
    IO.println s!"ATLAS_FORMAL_THEOREM\t{id}\t{solution}\t{trusted}\t{axStr}\t{stmt}"
  IO.println s!"ATLAS_FORMAL_THEOREM_COUNT\t{atlasAudited.length}"
