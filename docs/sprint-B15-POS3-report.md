# Sprint B15-POS3 — report (prereg-005, degenerate slices at Hankel order t = 2)

**Scope and warrant.** B15-POS3 tests the machinery, not physics. It checks the exact verifier
on the variance-zero slices of the degree-4 truncated Hausdorff problem (`μ2 = μ1²`). There the
feasible set on each slice is the single point `μ4 = a⁴`, the lower-wall formula `L` is undefined,
and `SCHUR_PIVOT_EXACT` meets zero pivots of every kind. It is a framework-capability benchmark,
not evidence about the universal layer (SPEC §2). System warrant: [KrN, Theorem II.2.3] as cited in
Curto–Fialkow 1991, Remark 4.4 (prereg-004 errata E1). The KrN chapter number is
VERIFY-BEFORE-USE against the book itself.

Preregistration: `docs/preregistrations/prereg-005-pos-degenerate-slice.md`, frozen at `017cf65`
(`prereg-005.freeze`), `prereg_sha256 b133b1d894e6a211…`.
Certificate: `certificates/b15_pos3_certificate.json`, **`b15-pos3-57497b7ab161`**, verdict
**PERMITTED** (the two on-slice points), SOUND/E0.

## Disclosure — checked once out-of-band before freeze
After the draft was merged (PR #21) and before freeze, Erick reproduced every pivot chain, all five
dual directions and values, the single-atom Farkas system and the two U3 neighbours at D, using the
repository's pivot routine. Every value agreed with the registration and none was changed. Before
freezing, I separately checked the E-side U3 neighbours (`1 ± 2⁻²⁰`) with a standalone rational
elimination (not the repository routine). They follow D's pattern (`−` fails H, `+` fails B). This
is also recorded in prereg-005.

## Gate
`python3 ci/run_all_certified.py`: **PASS**, 18 suites (17 + `test_b15_pos3.py`), 0 failures,
0 degradations. Existing certificates are byte-unchanged, including B15-POS2: the prereg-004 code
path in `certify_t2.py` is untouched. Manifest signature `733c5cf423de2a30`. A rerun of the test
regenerates the certificate identically, excluding the timestamp.

## Results
| Test | Result |
|---|---|
| U1 primal | D-on, E-on: PERMITTED, both PSD but not PD. Chains **equal the registered ones** (D: H 1, 0, 0 / B 1/4, 0; E: H 1, 0, 0 / B 0, 0, so the gap wall is itself degenerate) |
| U2 dual | All five exterior points REJECTED. Each stops at its registered pivot and violates only its registered wall. `v` and `vᵀMv` equal the registered values: D-lo `(−1/4, 0, 1)` → −1/64; D-hi `(−1/2, 1)` → −1/64; **D-off `(63/4, −32, 1)` → −1**, the first certified pass through the zero-diagonal/non-zero-row rule; E-lo `(−1, 0, 1)` → −1/64; E-hi `(0, 1)` → −1/64. Every Gram witness passes both the identity and the PSD check. The lift is recorded (stop index/kind, `j`, `t`, zero-pivot rows, solved indices, `w`) |
| U3 certified point | D: `1/16` PERMITTED; `65535/2²⁰` fails H, `65537/2²⁰` fails B. E: `1` PERMITTED; `1 − 2⁻²⁰` fails H, `1 + 2⁻²⁰` fails B. Each neighbour carries a certified zero-pivot dual. Stored as `certified_point`; no `certified_inner_interval` anywhere (enforced by the validator and by the standalone tool) |
| U4 Farkas | Single-atom system UNIQUE (`w = 1`) at D-on and E-on. INCONSISTENT at all five exterior points, each with a verified Farkas vector. Each carries the recorded premise (`μ2 = μ1²` ⇒ unique candidate `δ_a`) and states what it certifies ("not the unique candidate `δ_a`"; the Gram dual: "no representing measure on [0, 1]") |
| U5 formula controls | `det A = 0` on D and E, so L is undefined on both. U(D) = `1/16` = a⁴. `μ1 − μ2 = 0` on E, so U is also undefined there. The verdict route is pivots plus the Gram dual only; the formulas are computed only as controls |
| U6 no silent path | Confirmed at evaluation, as the draft predicted from reading the code: the prereg-004 routine refuses D-lo, E-lo and E-hi (AssertionError, singular leading block) and D-off (ValueError, zero-pivot indefiniteness). B15-POS3 certifies all seven through the registered construction |
| U7 | Standalone tool VERIFIED. Validator negatives 5/5: non-SPEC verdict, HEURISTIC+E0, REJECTED without an impossibility certificate, wrong prereg binding (prereg-004), `certified_inner_interval` on a point-valued slice. Standalone-tool negatives 5/5, each **re-hashed** so the content-hash check passes and the substantive check has to fire: tampered pivot, tampered lift (zero-pivot rows), tampered direction `v`, tampered Farkas vector, `certified_inner_interval` present |

None of the registered falsifiers K1–K5 was triggered.

## Finding outside this sprint (B15-POS2 test evidence, not its certificate)
B15-POS2's T8 tool-level negatives (tampered pivot, tampered Gram `Q1`, bracket excluding its wall)
tamper the certificate without re-computing its id. So the standalone tool's first failure line is
`certificate_id != content hash`, and that line is what the certificate records. When re-hashed in
a scratch run, all three are rejected by their substantive checks: pivot mismatch; Gram identity
fails and Q1 not PSD; stored wall ≠ registered. **The B15-POS2 verdict and certificate are
unaffected.** Only the recorded evidence of T8 is weaker than it reads. It is left as generated,
because regenerating the certificate with changed T8 evidence changes a certified artifact.
B15-POS3 uses the re-hashed form. Suggested follow-up: apply the same change to
`test_b15_pos_t2.py` in a separately reviewed commit, with the certificate regenerated.

## Changes
- `b15_surf/pos/certify_t2.py`: `zero_pivot_direction` (the registered construction with the
  recorded lift), `dual_functional_zp`, `farkas_single_atom`, `certify_point_zp`. The prereg-004
  functions are unchanged.
- `tests/test_b15_pos3.py`: U1–U7, with the registered values transcribed from prereg-005.
- `b15_surf/certificate.py`: B15-POS3 is bound to `prereg-005@<commit>`. The POS REJECTED rule
  covers B15-POS3. B15-POS3 may not carry `certified_inner_interval`.
- `schemas/b15_certificate.schema.json`: `B15-POS3` added to the benchmark enum.
- `tools/verify_b15_certificate.py`: `verify_pos3` reproduces `v` from the recorded lift with its
  own elimination and re-checks the single-atom system with its premise and the `certified_point`.
  It does not import `b15_surf`.
- `docs/verifier-ops-v0.1.md`: `DUAL_EXCLUSION_FUNCTIONAL` rev. 2 zero-pivot clarification (the
  general normalisation, the recorded lift, certificate-not-margin, Farkas role on variance-zero
  slices), plus a coverage-map row.

## Claims-table row (to be appended to `docs/claims-table.md` after sign-off)
| B15-POS3 zero-pivot / degenerate-slice verification | framework capability | BENCHMARK / E0 | prereg-005 (frozen 017cf65); variance-zero slices D (a = 1/2), E (a = 1); zero-pivot Gram duals (DUAL_EXCLUSION_FUNCTIONAL rev. 2 clarification), incl. zero-diagonal indefiniteness at D-off; single-atom Farkas complete on these slices; feasible set = point a⁴ at 2⁻²⁰ |
