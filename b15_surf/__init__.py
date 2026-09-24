"""B15-SURF — surfaceology / positivity-bootstrap integration benchmarks.

Three certified, pre-registered sub-benchmarks (docs/preregistrations/
prereg-002-surfaceology-benchmarks.md):

* ``b15_surf.pos``  — B15-POS: the repo's exact-PSD verifier scored against a
  published EFT-hedron forward-limit bound (R26, arXiv:2012.15849) with primal
  (Schur pivots) AND dual (exclusion functional) certificates.
* ``b15_surf.zero`` — B15-ZERO: exact ground-truth Tr(phi^3) tree amplitudes,
  hidden zeros, splits, and the delta-shift to the NLSM (R37), in Fraction.
* ``b15_surf.gid``  — B15-GID: blind grammar identification from pole / zero /
  split fingerprints only (B8 -> pir.fingerprints -> pir.candidates lineage).

Stdlib only; ``fractions.Fraction`` on every certified path (hard constraint 1).
Verdict vocabulary = SPEC §4 exactly (hard constraint 2). Layer = DOMAIN; no
atlas cell is touched (hard constraint 6).
"""
