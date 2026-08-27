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

/-! ## The 3×3 certificate and its preconditioner (ENG-008) -/

theorem atlas_pd_three_by_three_certificate : pd_three_by_three_certificate_statement := by
  intro a b c d e f d1 d2 d3 h1 h2 h3 hb1 hb2 hb3
  -- The certificate is inlined rather than bound with `have`: `have` is
  -- proof-irrelevant, so `cert.d1Lower` would stop reducing to `d1` and the
  -- hypotheses would no longer typecheck against the fields.
  have h := AtlasRH.posDef_of_certificate3
    { d1Lower := d1, d2Lower := d2, d3Lower := d3
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3 }
    (a := a) (b := b) (cc := c) (d := d) (e := e) (f := f) ⟨hb1, hb2, hb3⟩
  simpa [AtlasRH.sym3] using h

theorem atlas_diagonal_congruence_preserves_pd :
    diagonal_congruence_preserves_pd_statement := by
  intro x y z A hx hy hz h
  refine AtlasRH.posDef_of_diagonal_congruence hx hy hz ?_
  simpa [AtlasRH.diag3] using h

theorem atlas_preconditioned_certificate3 : preconditioned_certificate3_statement := by
  intro x y z A a b c d e f d1 d2 d3 hx hy hz hcong h1 h2 h3 hb1 hb2 hb3
  refine AtlasRH.posDef_of_preconditioned_certificate3
    { d1Lower := d1, d2Lower := d2, d3Lower := d3
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3 }
    hx hy hz ?_ ⟨hb1, hb2, hb3⟩
  simpa [AtlasRH.diag3, AtlasRH.sym3] using hcong

theorem atlas_diagonal_congruence_preserves_index :
    diagonal_congruence_preserves_index_statement := by
  intro x y z A k hx hy hz
  have h := AtlasRH.posIndexAtLeast_diagonal_congruence_iff
    (x := x) (y := y) (z := z) (A := A) (k := k) hx hy hz
  simpa [AtlasRH.PosIndexAtLeast, AtlasRH.qform, AtlasRH.diag3] using h

theorem atlas_diagonal_congruence_preserves_rank :
    diagonal_congruence_preserves_rank_statement := by
  intro x y z A hx hy hz
  have h := AtlasRH.rank_diagonal_congruence A hx hy hz
  simpa [AtlasRH.diag3] using h

/-! ## The generalized gap (ENG-009) -/

theorem atlas_generalized_rayleigh : generalized_rayleigh_statement := by
  intro n G M lam h v
  simpa [AtlasRH.qform] using
    AtlasRH.rayleigh_lower_of_shifted_psd (G := G) (M := M) (lam := lam) h v

theorem atlas_generalized_pencil_congruence :
    generalized_pencil_congruence_statement := by
  intro n S G M lam hS
  constructor
  · exact AtlasRH.shifted_posDef_of_congruence hS
  · exact AtlasRH.shifted_posDef_congruence hS

theorem atlas_preconditioned_gap_certificate3 :
    preconditioned_gap_certificate3_statement := by
  intro x y z lam G M a b c d e f d1 d2 d3 hx hy hz hcong h1 h2 h3 hb1 hb2 hb3 v
  have h := AtlasRH.gap_of_preconditioned_certificate3
    { d1Lower := d1, d2Lower := d2, d3Lower := d3
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3 }
    (x := x) (y := y) (z := z) (lam := lam) (G := G) (M := M)
    hx hy hz ?_ ⟨hb1, hb2, hb3⟩ v
  · simpa [AtlasRH.qform] using h
  · simpa [AtlasRH.diag3, AtlasRH.sym3] using hcong

/-! ## The 4×4 block (ENG-010) -/

theorem atlas_pd_four_by_four_certificate : pd_four_by_four_certificate_statement := by
  intro a b c d e f g h i j d1 d2 d3 d4 h1 h2 h3 h4 hb1 hb2 hb3 hb4
  have hpd := AtlasRH.posDef_of_certificate4
    { d1Lower := d1, d2Lower := d2, d3Lower := d3, d4Lower := d4
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3, h_d4_pos := h4 }
    (a := a) (b := b) (c := c) (d := d) (e := e) (f := f) (g := g) (h := h)
    (i := i) (j := j)
    ⟨hb1, by simpa [AtlasRH.minor2of4] using hb2,
     by simpa [AtlasRH.minor3of4] using hb3,
     by simpa [AtlasRH.minor4of4] using hb4⟩
  simpa [AtlasRH.sym4] using hpd

theorem atlas_preconditioned_certificate4 : preconditioned_certificate4_statement := by
  intro w x y z A a b c d e f g h i j d1 d2 d3 d4 hw hx hy hz hcong
    h1 h2 h3 h4 hb1 hb2 hb3 hb4
  refine AtlasRH.posDef_of_preconditioned_certificate4
    { d1Lower := d1, d2Lower := d2, d3Lower := d3, d4Lower := d4
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3, h_d4_pos := h4 }
    hw hx hy hz ?_
    ⟨hb1, by simpa [AtlasRH.minor2of4] using hb2,
     by simpa [AtlasRH.minor3of4] using hb3,
     by simpa [AtlasRH.minor4of4] using hb4⟩
  simpa [AtlasRH.diag4, AtlasRH.sym4] using hcong

theorem atlas_preconditioned_gap_certificate4 :
    preconditioned_gap_certificate4_statement := by
  intro w x y z lam G M a b c d e f g h i j d1 d2 d3 d4 hw hx hy hz hcong
    h1 h2 h3 h4 hb1 hb2 hb3 hb4 v
  have hq := AtlasRH.gap_of_preconditioned_certificate4
    { d1Lower := d1, d2Lower := d2, d3Lower := d3, d4Lower := d4
      h_d1_pos := h1, h_d2_pos := h2, h_d3_pos := h3, h_d4_pos := h4 }
    (w := w) (x := x) (y := y) (z := z) (lam := lam) (G := G) (M := M)
    hw hx hy hz ?_
    ⟨hb1, by simpa [AtlasRH.minor2of4] using hb2,
     by simpa [AtlasRH.minor3of4] using hb3,
     by simpa [AtlasRH.minor4of4] using hb4⟩ v
  · simpa [AtlasRH.qform] using hq
  · simpa [AtlasRH.diag4, AtlasRH.sym4] using hcong

/-! ## The 5×5 block (ENG-011) -/

theorem atlas_shifted_positivity_transfer : shifted_positivity_transfer_statement := by
  intro n G M lam hlam hM h
  exact AtlasRH.posDef_of_shifted_posDef_add hlam hM h

theorem atlas_shifted_shift_monotone : shifted_shift_monotone_statement := by
  intro n G M a b hba hM h
  exact AtlasRH.shifted_posDef_of_le hba hM h

theorem atlas_nested_gap_regression : nested_gap_regression_statement := by
  intro n m G M lam h e x
  simpa [AtlasRH.qform] using AtlasRH.rayleigh_lower_restricts h e x

theorem atlas_schur_witness_block : schur_witness_block_statement := by
  intro n A c d y hA hy hd x t hxt
  simpa [AtlasRH.borderedForm, AtlasRH.qform] using
    AtlasRH.borderedForm_pos A hA hy hd x t hxt

/-! ## Rank–trace -/

theorem atlas_rank_trace_hs : rank_trace_hs_statement := by
  intro l h
  have h' := AtlasRH.rank_trace_zero_Q l h
  simpa [AtlasRH.rankTraceRhs, AtlasRH.traceOf, AtlasRH.hsSqOf, AtlasRH.rankOf] using h'

end Comparator
