# Sprint B15-POS2 — report (prereg-004, Hankel order t = 2)

**Scope and warrant.** This is a reproduction of a known mathematical bound, the degree-4 gapped
forward-limit EFT-hedron walls (arXiv:2012.15849 eqs. 7.45/7.46), with the program's own exact
verifier. It measures the machinery (`SCHUR_PIVOT_EXACT`, `DUAL_EXCLUSION_FUNCTIONAL` rev. 2,
bracketing) on a case with ground truth. It is not evidence about the universal layer (SPEC §2).

Preregistration: `docs/preregistrations/prereg-004-pos-hankel-t2.md`, frozen at `1ccea83`
(`prereg-004.freeze`), `prereg_sha256 8d68d570fa13020a…`. The warrant for the two-matrix
system is Kreĭn–Nudel′man Thm III.2.3, as stated in Curto–Fialkow 1991, Remark 4.4.
Certificate: `certificates/b15_pos_t2_certificate.json`, **`b15-pos2-6facf90a799a`**,
verdict **PERMITTED**, SOUND/E0.

## Disclosure — checked once out-of-band before freeze
After prereg-004 was merged (PR #19) and before it was frozen, Erick re-derived every wall, pivot
chain and quadrature value with exact arithmetic. That check included running `b15_surf.pos`'s
pivot routine on the six t = 2 points and the t = 1 verifier on their projections. Everything
agreed with the registered expectations, and no registered value was changed. The walls, points
and expected pivots had been fixed and merged before that run. This is also recorded in
prereg-004 itself.

## Gate
`python3 ci/run_all_certified.py`: **PASS**, 17 suites (16 + `test_b15_pos_t2.py`), 0 failures,
0 degradations. Existing certificates are byte-unchanged. Manifest signature `23beaab6dbf68071…`.
The regenerated B15-POS2 certificate hash is `ba84736f…`. A rerun of the test regenerates the
certificate byte-for-byte (excluding the timestamp).

## Results
| Test | Result |
|---|---|
| T1 primal | S1-int (μ4 = 1/5), S2-int (91/256): PERMITTED. Pivot chains **equal the registered ones** (H: 1, 1/12, 1/180 / B: 1/6, 1/120; H: 1, 7/72, 3/896 / B: 7/48, 1/448) |
| T2 dual | all four exterior points REJECTED. Negative pivots equal the registered values (−1/72, −1/48, −1/1792, −1/1792), and each point violates only its registered wall. Every dual is `DUAL_EXCLUSION_FUNCTIONAL` rev. 2 with identity ✓ and `Q0`, `Q1` PSD-certified. Example (S1-lo): `y = (1/36, −1/3, 4/3, −2, 1)`, i.e. `(x² − x + 1/6)²`, the squared Gauss–Legendre node polynomial, gives −1/72 |
| T3 brackets | S1: `[152917/786432, 1146881/5898240] ∋ 7/36`, `[819199/3932160, 109227/524288] ∋ 5/24`. S2: `[369225/1048576, 184613/524288] ∋ 631/1792`, `[93769/262144, 2625537/7340032] ∋ 641/1792`. All widths ≤ 2⁻²⁰, and no endpoint lands on a wall. G1 not triggered |
| T4 quadrature | Gauss–Legendre 2-point μ4 = 7/36 = L(S1); Simpson μ4 = 5/24 = U(S1). G4 not triggered |
| T5 t = 1 blindness | the prereg-002 verifier PERMITS all six (μ1, μ2) projections, while t = 2 rejects the four exterior points: **t = 2 adds content** |
| T6 negative control | μ4 1/5 → 31/180 at S1-int flips to REJECTED with a certified Gram dual |
| T7 consistency with eq. (7.46) | 19 contiguous tower minors at each interior point, all > 0 (min 1/2160, 1/3072). G5 not triggered |
| T8 | standalone tool VERIFIED. Schema/tool negatives 7/7: non-SPEC verdict, HEURISTIC+E0, REJECTED without an impossibility certificate, wrong prereg binding, tampered pivot, tampered Gram `Q1`, bracket not containing the wall |

The registered falsifiers G1–G5 were not triggered.

## Changes
- `b15_surf/pos/certify_t2.py`: primal on the 3×3 Hankel and 2×2 localizing matrices; pivot-derived
  rank-1 Gram duals (rev. 2) checked with `verify_solution` and `psd_certificate`; μ4 bisection;
  quadrature and tower routes.
- `tests/test_b15_pos_t2.py`: T1–T8, with the registered values transcribed from prereg-004.
- `b15_surf/certificate.py`: the prereg binding is per benchmark (B15-POS2 → `prereg-004@<commit>`,
  others unchanged; unfrozen refs rejected). The POS REJECTED rule now also covers B15-POS2.
- `schemas/b15_certificate.schema.json`: `B15-POS2` added to the benchmark enum.
- `tools/verify_b15_certificate.py`: `verify_pos2` recomputes the pivots, Gram identity and PSD,
  brackets against the registered walls, and the quadrature, without importing `b15_surf`.
- `docs/verifier-ops-v0.1.md`: `DUAL_EXCLUSION_FUNCTIONAL` rev. 2 (complete by Markov–Lukács,
  supersedes rev. 1 as the canonical form, conditions (i) rationality and (ii) both checks).

## Proposed claims-table row (append after sign-off)
| B15-POS2 verifier vs EFT-hedron gapped bound at t = 2 | recovered known result | BENCHMARK / E0 | two slices; Gram-form duals (DUAL_EXCLUSION_FUNCTIONAL rev. 2); brackets ≤ 2⁻²⁰; S1 walls = Gauss–Legendre / Simpson; t = 1 blind to all four exterior points |

## Backlog
prereg-005: a slice with singular `A` (μ2 = μ1²), which exercises the zero-pivot path (P4-OI-2).
