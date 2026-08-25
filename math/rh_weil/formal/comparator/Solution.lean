/-
ATLAS-RH-ENG-007 §5 Layer C — the proofs.

Every theorem below is typed against a name from `comparator/TrustedStatements`.
That typing is the Layer D comparison, and Lean performs it: if a future edit
weakens one of the propositions in `AtlasRH/Statements.lean`, the theorem here
still typechecks but its recorded statement hash changes and
`scripts/check_formal_manifest.py` fails; if instead someone tries to prove a
different proposition under one of these names, the elaborator rejects it.

No `sorry` appears in this file, and `comparator/PrintAxioms.lean` reports the
axiom dependencies of every theorem named here.

No RH proof claim is made. Everything below is finite linear algebra over `ℝ`.
-/
import AtlasRH
import comparator.TrustedStatements

namespace Comparator

open Matrix AtlasRH

/-! ## Sylvester's law of inertia -/

theorem atlas_inertia_congruence_positive : inertia_congruence_positive_statement := by
  intro n S A k hS
  have h := AtlasRH.posIndexAtLeast_congruence_iff (S := S) (A := A) (k := k) hS
  simpa [AtlasRH.PosIndexAtLeast, AtlasRH.qform] using h

theorem atlas_inertia_congruence_negative : inertia_congruence_negative_statement := by
  intro n S A k hS
  have h := AtlasRH.posIndexAtLeast_congruence_iff (S := S) (A := -A) (k := k) hS
  simpa [AtlasRH.PosIndexAtLeast, AtlasRH.qform] using h

theorem atlas_inertia_congruence_rank : inertia_congruence_rank_statement := by
  intro n S A hS
  exact AtlasRH.rank_congruence A hS

/-! ## Leading-principal-minor criteria -/

theorem atlas_pd_two_by_two : pd_two_by_two_statement := by
  intro a b c
  exact AtlasRH.posDef_sym2_iff (a := a) (b := b) (c := c)

theorem atlas_pd_three_by_three : pd_three_by_three_statement := by
  intro a b c d e f
  have h := AtlasRH.posDef_sym3_iff (a := a) (b := b) (c := c) (d := d) (e := e) (f := f)
  simpa [AtlasRH.sym3, AtlasRH.minor2, AtlasRH.minor3] using h

/-! ## Certificate semantics -/

theorem atlas_certificate_even2_implies_pd : certificate_even2_implies_pd_statement := by
  intro a b d g00lo g00hi g0blo g0bhi gbblo gbbhi g00Lower detLower hg00pos hdetpos
    hg00bound hdetbound ha_lo ha_hi hb_lo hb_hi hd_lo hd_hi
  have hg00le : g00lo ≤ g00hi := le_trans ha_lo ha_hi
  have hg0ble : g0blo ≤ g0bhi := le_trans hb_lo hb_hi
  have hgbble : gbblo ≤ gbbhi := le_trans hd_lo hd_hi
  let c : AtlasRH.PD2Certificate :=
    { g00 := ⟨g00lo, g00hi, hg00le⟩
      g0b := ⟨g0blo, g0bhi, hg0ble⟩
      gbb := ⟨gbblo, gbbhi, hgbble⟩
      g00Lower := g00Lower
      detLower := detLower
      h_g00_pos := hg00pos
      h_det_pos := hdetpos
      h_g00_bound := hg00bound
      h_det_bound := hdetbound }
  have h := AtlasRH.posDef_of_certificate c (a := a) (b := b) (d := d)
    ⟨⟨ha_lo, ha_hi⟩, ⟨hb_lo, hb_hi⟩, ⟨hd_lo, hd_hi⟩⟩
  simpa [AtlasRH.sym2] using h

/-! ## RH/Weil exact identities -/

theorem atlas_weil_basis_parity : weil_basis_parity_statement := by
  intro L
  exact ⟨AtlasRH.basisOne_even L, AtlasRH.basisB_even L, AtlasRH.basisQ1_odd L,
    AtlasRH.basisB3_odd L⟩

theorem atlas_odd_degree3_cross_block : odd_degree3_cross_block_statement := by
  intro L pair bilin_neg refl_inv f g hf hg
  exact AtlasRH.cross_block_vanishes pair bilin_neg refl_inv hf hg

theorem atlas_odd_degree3_factorization : odd_degree3_factorization_statement := by
  intro e00 e0b ebb o11 o1b obb
  have h := AtlasRH.parity_block_det e00 e0b ebb o11 o1b obb
  simpa [AtlasRH.sym2] using h

/-! ## Rank–trace -/

theorem atlas_rank_trace_hs : rank_trace_hs_statement := by
  intro l h
  have h' := AtlasRH.rank_trace_zero_Q l h
  simpa [AtlasRH.rankTraceRhs, AtlasRH.traceOf, AtlasRH.hsSqOf, AtlasRH.rankOf] using h'

end Comparator
