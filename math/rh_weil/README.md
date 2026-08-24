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

## Normalization is frozen (WO-RH-17, P0)

The pole block is **adjudicated**. The adopted form is derived from the explicit
formula,

`G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+`,  `E_i^± = ∫_0^L h_i(x) e^{±x/2} dx`,

and the legacy even block `(√3/2)(v₊v₊ᵀ+v₋v₋ᵀ)` is **rejected**: it equals the
adopted pole times `(√3/2)cosh(L/2)`, a factor that is 1 only at `L = log 3` — a
calibration fitted at one test point (+8.25 % at `L = log 4`). The odd pivot
`−8A²` was already correct and is unchanged.

See [`docs/NORMALIZATION_ADJUDICATION_v0.1.md`](docs/NORMALIZATION_ADJUDICATION_v0.1.md),
`certificates/normalization_adjudication.json` and
`certificates/normalization_crosscheck.json` (four independent providers).

Certificates that depended on the rejected block are
`QUARANTINED_NORMALIZATION_ADJUDICATION` — preserved, not deleted, not relabelled —
and the PIR bridge refuses to promote them. They must be **regenerated** under the
frozen normalization (WO-RH-19/20), never reinterpreted.

```bash
python3 scripts/derive_normalization.py           # derivation + adjudication certificate
python3 scripts/run_normalization_crosscheck.py   # four-way cross-check (--no-arch to skip the slow route)
python3 scripts/quarantine_normalization.py       # idempotent quarantine
```

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

