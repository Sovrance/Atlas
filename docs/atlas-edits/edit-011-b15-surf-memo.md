# Atlas Edit 011 — B15-SURF memo (no cell change)

**Status: MEMO — records three BENCHMARK-tier certificates at layer DOMAIN; changes no cell,
adds no H or P.** (Edit 010 remains reserved by edit-009's decision rule.)

**Certificates:** `certificates/b15_pos_certificate.json` (6/6), `certificates/b15_zero_certificate.json`
(7/7), `certificates/b15_gid_certificate.json` (7/7, T0–T6). Pre-registration `prereg-002`
(`docs/preregistrations/prereg-002-surfaceology-benchmarks.md`, DRAFT-FOR-FREEZE) is bound into
every certificate by commit hash and content hash. Standalone re-verifier:
`tools/verify_b15_certificate.py`.

**No-circularity statement.** Everything below reproduces DOMAIN-layer results — known
perturbative QFT identities and a published EFT bound. None of it is evidence for the
UNIVERSAL-layer conjecture that 𝕽 exists (SPEC §2); it scores the program's *verifier and
identification machinery* on cases with ground truth.

**Findings (filed, not promoted):**
1. **B15-POS — verifier vs. external bound (R26).** For the first time the exact-PSD verifier
   (`SCHUR_PIVOT_EXACT`) was scored against a bound the program did not construct — the
   EFT-hedron forward-limit two-sided statement `μ₁² ≤ μ₂ ≤ μ₁` (arXiv:2012.15849 eqs. 7.19,
   7.35; slice μ₁ = ½). Interior PERMITTED with pivot witness; both exterior points REJECTED
   with (a) a negative pivot and (b) an exact **dual exclusion functional** whose
   nonnegativity on the moment cone is certified by an exact polynomial identity; both walls
   bracketed to ≤ 2⁻²⁰; negative control flips. R26's request ("add dual/exclusion
   functionals to the verifier layer") is now executed, with one recorded caveat: the
   repo's `verify_farkas` (linear-equality Farkas) is *not* the dual of a conic constraint;
   it is used only where genuine (flat boundary, single-atom system). Verdict PERMITTED,
   taint `asm:B15-OI-1-pending` until the target bound is confirmed.
2. **B15-ZERO — exact ground truth (R37).** Route A (triangulation sum, 5/14/132 terms) equals
   Route B (factorization recursion) on 200 points per n ∈ {5,6,8}; all n(n−3)/2 registered
   hidden-zero loci (2312.16282 eq. 6.14 + cyclic images) annihilate Aₙ exactly; the 2-split
   (2405.09608 eq. 1.7) holds exactly on every registered split; and the δ-shift prescription
   `lim δ^{2n−2} A^δ_{2n} = A^{NLSM}_{2n}` (2401.05483 eq. 2) equals an independent
   Feynman-rule computation **exactly at 4, 6 and 8 points**, with the claimed cancellation of
   all lower orders in 1/δ verified. Category: recovered known result, E0.
3. **B15-GID — blind identification.** From the label-stripped dataset alone the pass
   recovers exact pole sets, zero loci, split ranks and the deformation parameter
   (δ = 7/3 on the withheld shifted class; `NONIDENTIFIABLE(delta-underdetermined)` on the
   NLSM class, as registered). **Outcome A** obtains on the primary correlator at the
   held-out n = 8 (all three classes `FORCED`, similarity = confidence = 1; F4 not
   triggered). **Outcome B** obtains on the registered zeros-only correlator:
   `OBSERVATIONALLY_EQUIVALENT([NLSM, TrPhi3, TrPhi3_shifted] mod δ-shift)` — the theories
   share every hidden zero by construction, so *fingerprint-only identification from zeros
   is structurally limited for deformation-related grammars*; what separates them is the
   pole set (NLSM lacks `X_{e,e}`, `X_{o,o}` poles) and the δ-class. The δ-insensitive
   invariants alone leave `[TrPhi3, TrPhi3_shifted]` `AMBIGUOUS`, exactly as the two-hash
   design predicts. Redundancy metric (P1) at n = 8: 132/43 (Tr φ³), 132/31 (NLSM) — recorded,
   not acted on.

**Errata / transcription notes (filed, not repaired):**
- **E1 (arXiv:2405.09608 eq. 2.13).** Our literal reading of the index ranges
  "`c_{a−1,c−1} = 0`, `c_{b−1,c−1} = 0` for all a ∈ A, b ∈ B, c ∈ C" does not reproduce the
  paper's own explicit vanishing sets (eq. 2.5: `{c_{1,3}, c_{3,5}}` at n = 6; the single
  `c_{1,3}` implied by eq. 2.8 at n = 5). The registered split locus is therefore derived
  exactly from the linear map (2.10) and reproduces (2.5) verbatim; the factorization (1.7)
  is unaffected. Most likely a labelling/offset convention we have not matched — flagged for
  a second reading, not asserted as an error in the source.
- **E2 (naming).** arXiv:2405.09608 §2.2 uses `y` for S1 and `x` for S2 while the general
  formula (2.10) is written with `x` on the `(i..k)` surface; prereg-002 fixes `x := P = (i..k)`,
  `y := Q` and guards it with the explicit examples.
- **E3 (sign convention).** The NLSM Feynman route's overall sign is fixed by the metric
  signature; prereg-002 stipulation S1 pins it at 4 points to the letter's `−(X₁₃+X₂₄)`.
- **R39 title.** arXiv:2408.04720 (Cheung, Dersy, Schwartz) is titled *Learning the
  Simplicity of Scattering Amplitudes* (verified 2026-09-23).
- No discrepancy was found in the arXiv IDs of §12 of the work order (all re-verified).

**Registered decision rules (unchanged from prereg-002):** F1–F4. None triggered.

**Ledger:** R36–R40 appended (`docs/references.md`); R26 covers WP1 and is not duplicated.

**Not done (recorded):** WP2.3 Route C (ABHY canonical form, R38) — not attempted; WP2.8
YM-scaffolded rows — skipped pending §9-OI-3 (prereg records the three-class scope).

**Addendum (2026-09-24, §9 decisions — `docs/preregistrations/prereg-003-b15-open-item-amendments.md`).**
The target bound is confirmed and `asm:B15-OI-1-pending` is removed. The B15-GID
identification verdict is now `PERMITTED` (menu-relative identification) rather than `FORCED`,
and its threshold is exact agreement (the result is unchanged; every score was 1). The POS dual
is named as verifier op 7, `DUAL_EXCLUSION_FUNCTIONAL`. The body above is kept as written at
sprint time.
