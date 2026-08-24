# RH / Weil Positivity Program (Research Notebook V2)

This directory integrates the RH Research Notebook V2 into Atlas as the reproducible source behind the Constant Atlas positivity-verifier example (`G00 > 0`, Schur pivot positivity).

## Scope and claim boundary

This program studies finite-dimensional polynomial compressions of the localized Weil quadratic form. **It does not prove the Riemann Hypothesis.** A finite block certificate is evidence only for that stated block, interval, normalization, and cutoff/representation.

Accepted assembly convention:

`G = G^0 - G^p + G^infinity`.

Current certified/research cell: `L in [log(3), log(4)]`.

## Evidence policy

- E0 / SOUND: exact algebraic identities independently re-derived/tested.
- E1 / SOUND: interval-certified numerical statements with reproducible interval engine and certificate.
- E3 / HEURISTIC: floating scans and localization of candidate minima.
- Spectral/quantum analogies: diagnostics only; never proof warrants.

Imported notebook claims are intentionally marked `IMPORTED_PENDING_REGENERATION` until an Atlas agent reproduces them from the code in this repository. Do not promote them merely because they appeared in the prior notebook transcript.

## Normalization is frozen (WO-RH-17, P0) and now lives in one module (ENG-004)

The pole block is **adjudicated**. The adopted form is derived from the explicit
formula,

`G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+`,  `E_i^± = ∫_0^L h_i(x) e^{±x/2} dx`,

and the legacy even block `(√3/2)(v₊v₊ᵀ+v₋v₋ᵀ)` is **rejected**: it equals the
adopted pole times `(√3/2)cosh(L/2)`, a factor that is 1 only at `L = log 3` — a
calibration fitted at one test point (+8.25 % at `L = log 4`). The odd pivot
`−8A²` was already correct and is unchanged.

`src/pole.py` is the **single** implementation (`laplace_plus`, `laplace_minus`,
`pole_gram_entry`, `pole_gram_matrix`, `pole_gram_entry_dL`). Every production
path routes through it. The rejected candidate has been moved out of production
into `src/rejected_pole.py`, which is archival: `tests/test_production_imports.py`
fails CI if anything under `src/` imports it or spells its scale in executable
code. Two copies of the rejected block existed before this change — the pole
assembly and, separately, its `L`-derivative in the jet module.

See [`docs/NORMALIZATION_ADJUDICATION_v0.1.md`](docs/NORMALIZATION_ADJUDICATION_v0.1.md),
`certificates/normalization_adjudication.json` and
`certificates/normalization_crosscheck.json`. That comparison is a
**three-way internal cross-check** (explicit formula, compact real space, direct
Fourier); the external Connes/CvS provider quantifies no projection or truncation
error, so it reports `NOT_COMPARABLE` and never certifies.

Certificates that depended on the rejected block are
`QUARANTINED_NORMALIZATION_ADJUDICATION` — preserved, not deleted, not relabelled —
and the promotion predicate refuses them. The quarantine is enforced at the point
of write, so re-running a legacy certifier cannot erase it. They must be
**regenerated** under the frozen normalization, never reinterpreted.

## Scalar E1 canary (ATLAS-RH-ENG-004)

The scalar cell entry is the first artifact regenerated under Candidate A and the
only one ENG-004 releases from quarantine. `certificates/e1_scalar_log3_log4.json`
carries a **uniform** rigorous lower bound over `[log 3, log 4]`, computed with
python-flint/Arb — there is no mpmath path that may emit E1.

The bound rests on convexity, which is algebraic rather than numeric. `Gp` is
piecewise linear with breakpoints exactly at the cell endpoints, the pole
contributes `G0'' = 4cosh(L/2)`, and the archimedean term's cosine transform gives
`Ginf'' = −e^{L/2}/sinh L`, so

`G00''(L) = 4cosh(L/2) − e^{L/2}/sinh(L) = 2(r³−r−1)/(√r(r²−1))`,  `r = eᴸ`,

which is exactly the repository's E0 curvature `W00''`, already proved positive on
the cell. That identity is itself independent evidence for the adjudication: the
`4cosh(L/2)` term is produced by Candidate A's pole and by nothing else, so the
rejected block cannot reproduce the certified curvature.

Regression note: the historical notebook minimum `0.0753795566…` — which the
rejected calibration had put out of reach, holding this entry above ≈0.1276 on the
cell — falls **inside** the recovered enclosure. It is recorded as regression
evidence only; the acceptance gate is the positive bound, never a fitted constant.

Degree-1, compact degree-2 and both T=84 E1 artifacts remain quarantined. ENG-005
recovers them.

```bash
python3 scripts/derive_normalization.py           # derivation + adjudication certificate (requires SymPy)
python3 scripts/run_normalization_crosscheck.py   # three-way internal cross-check (--no-arch to skip the slow route)
python3 scripts/quarantine_normalization.py       # idempotent quarantine
python3 scripts/run_rigorous_scalar.py            # the rigorous scalar path, in order
```

### Rigorous dependencies

SymPy and python-flint are **required** for the rigorous/research path, not
optional extras. A missing one now fails the job rather than degrading an artifact
to weaker evidence and dirtying the tree:

```bash
pip install -r requirements-rigorous.txt
```

`scripts/run_rh_weil_suite.py` is the fast path. Passing it does **not** mean the
rigorous certificates are current — it does not re-derive the scalar canary. Read
`scripts/run_rigorous_scalar.py`'s exit code before believing an E1 claim is fresh.

## Current priority

Regenerate the E1 chain under the frozen normalization (WO-RH-19), then the direct
Fourier-side uniform degree-2 certificate at `T=84` (WO-RH-20), then the
midpoint-odd degree-3 block (WO-RH-27, gated on P0/P1).

**Executed in-repo so far (see `certificates/work_order_status.json`):**
- WO-RH-01/03/04 — exact identities, f1 audit, bubble block (E0)
- WO-RH-02 — scalar cell `[log 3, log 4]` algebraic `W00''>0` (E0)
- WO-RH-05 — stable `H0`/`Hb` + L-jets; E3 energy-probe scan only (interval E1 **open**)
- WO-RH-06/07 — regenerated E0 certificates + dedicated runner (does not expand root CI)

Run:

```bash
python math/rh_weil/scripts/run_rh_weil_suite.py
```

See `AGENT_INSTRUCTIONS.md` and `notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`.

Optional external oracle (not required at runtime): `external/` wraps
`connes-cvs` for independent cross-checks of shared Weil ingredients. See
`external/CONNES_CVS_MAPPING.md` and `external/PROVENANCE.md`. Do not compare
Galerkin matrix entries to Atlas polynomial Gram blocks by index.

Regenerate the external cross-validation certificate (research env):

```bash
pip install 'connes-cvs==0.3.1' python-flint mpmath
python math/rh_weil/scripts/run_connes_cvs_crosschecks.py
```

