# Preregistration 003 — B15-SURF open-item amendments to prereg-002

**Status: ACTIVE on merge.** Bound into the B15 certificates by content hash
(`inputs.prereg_amendment = {ref: prereg-003, sha256}`, via
`b15_surf/certificate.py::amendment_ref`); prereg-002 stays frozen and unedited
(`prereg_ref` / `prereg_sha256` unchanged in meaning).

**Why a new preregistration.** prereg-002 §"Stop rules": *no threshold … may be edited
after freeze; a needed change is a new preregistration (prereg-003)*. Two of the §9
answers below change registered rules (A1 threshold, A2 verdict vocabulary), so they are
filed here rather than patched into prereg-002.

**Disclosure — these amendments are post-observation.** They were decided after the
certified B15 run. Neither changes any outcome: every primary GID score at the held-out
n = 8 was exactly 1 (so "= 1" and "≥ 9/10" select identically), and A2 renames the label
on the same adjudication (the lattice already emitted `PERMITTED`). They are recorded as
honesty/vocabulary corrections, not as a re-run that could have gone differently.

Decisions by Erick (2026-09-24), §9 of the B15-SURF work order.

## Amendments to registered rules

- **A1 (OI-2) — GID threshold = exact agreement.** Replaces prereg-002 §GID
  "similarity ≥ 9/10 and confidence ≥ 9/10 **[proposal]**" with **similarity = 1 and
  confidence = 1**. Every GID detection is exact, so a fractional agreement can arise only
  from a wrong menu entry, never from noise; a 9/10 threshold would pretend to a statistical
  tolerance the pass does not have. The POS bracket width `1/2^20` is **accepted as
  registered** (a pure reporting choice; bisection is exact at any width).
  `asm:B15-OI-2-thresholds-pending` is removed.

- **A2 (OI-9) — a single surviving menu member is `PERMITTED`, not `FORCED`.** Replaces
  prereg-002 §GID Outcome A "`FORCED` on the grammar label" with **`PERMITTED` (menu-relative
  identification)**. SPEC §4 reserves `FORCED` for unique determination by constraints; one
  survivor of a finite menu is identification relative to that menu — B8 semantics and what
  `pir.candidates.evaluate` already returns. The certificate and the PIR lattice fact now use
  one vocabulary (reconciliation D3 closed). δ recovery keeps `FORCED`: δ is uniquely fixed
  by exact pole-location constraints, which is SPEC §4 `FORCED`.

- **A3 (OI-1) — POS target bound confirmed as registered.** The lowest-order two-sided
  forward-limit statement of arXiv:2012.15849 (Hankel wall eq. 7.19 + gap wall eq. 7.35) at
  μ₁ = ½, walls 1/4 and 1/2. `asm:B15-OI-1-pending` is removed. Known weakness, recorded: at
  Hankel order t = 1 the bound is nearly trivial (μ₁² ≤ μ₂ ≤ μ₁). The stronger test (t = 2,
  moments through a_{6,0}, 3×3 Hankel + shifted Hankel, walls no longer single inequalities)
  is a **separate future preregistration (prereg-004)**, not a silent upgrade of this one.

## Rulings with no rule change (recorded here for one place of record)

- **OI-3** — WP2.8 YM-scaffolded rows: separate sprint (needs the stringy integrand,
  α′δ = 1; a different generator with its own transcription burden). B15-ZERO stays E0.
- **OI-4** — open B16, narrowly: the two legs of `docs/notes/constraint-manifest-representation-v0.md`
  §4 (B10 dilation chart k = 1; B15-ZERO (x, y) sub-surface chart), category *known theorem
  (re-expression)*, no cell movement; deliverable is the per-fact schema field
  `constraint_manifest: [...]`.
- **OI-5** — claims-table rows for B15-ZERO and B15-GID approved; B15-POS row held until
  OI-1, which A3 resolves.
- **OI-6** — `ground_truth_route: null` on legacy certificates: yes, separate sprint,
  schema-only; the degradation predicate treats added keys as non-degrading.
- **OI-7** — pre-flight: no decision; the gate's Python dependencies (`mpmath`, `numpy`,
  `scipy`) are recorded in `ci/requirements.txt`.
- **OI-8** — the dual-functional + Farkas-companion form is accepted as the verifier-layer
  dual. It is named as the seventh verifier op `DUAL_EXCLUSION_FUNCTIONAL`
  (`docs/verifier-ops-v0.1.md`, rev. 1); no rational-SDP primitive is required at these
  truncation orders. The B15-POS dual object carries `"op": "DUAL_EXCLUSION_FUNCTIONAL"`.
