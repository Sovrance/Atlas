---
status: HISTORICAL
superseded_by: ../../AGENT_INSTRUCTIONS.md
superseded_at: ATLAS-RH-ENG-007
---

> **HISTORICAL — do not follow these instructions.**
>
> This is the original integration work order (WO-RH-01 … WO-RH-08), preserved
> verbatim as provenance. It is **out of date in ways that matter**: it says
> WO-RH-05 is open and that degree 3 must not start until WO-RH-05 closes. Both
> statements were true when written and are false now — WO-RH-05 was recovered by
> ENG-005 and degree 3 is E1-certified positive definite by ENG-006.
>
> The live entrypoint is [`AGENT_INSTRUCTIONS.md`](../../AGENT_INSTRUCTIONS.md).
> Current machine-readable status lives in
> [`certificates/work_order_status.json`](../../certificates/work_order_status.json).
>
> It is kept because WO-RH-17 forbids deleting superseded evidence: the sequence
> of decisions is itself part of the record.

---

# Agent Work Order — RH Research Notebook V2 -> Atlas

## Mission

Turn the imported RH/Weil notebook state into a reproducible Atlas-native certified program without weakening Atlas's epistemic rules or overwriting existing work.

## Non-negotiable rules

1. Never emit or imply `RH PROVED` from finite polynomial blocks.
2. Keep the sign convention `G = G0 - Gp + Ginf` under an explicit normalization audit.
3. `uncertified != pass`: imported numerical values remain pending until regenerated.
4. Exact identities and interval claims must be separated from floating scans.
5. Every E1 claim must record interval, precision, algorithm/version, source hash, and outward-rounded lower/upper bounds.
6. Spectral/quantum material is diagnostic only.
7. Do not modify unrelated Atlas certificates or benchmark outputs.

## Implementation order

### WO-RH-01 — normalization + exact identities
- Implement and test the polynomial overlap formula
  `C_ij(a,L)=sum_{r=0}^j binom(j,r) a^(j-r)(L-a)^(i+r+1)/(i+r+1)`.
- Verify midpoint-odd and bubble kernels in `src/core.py`.
- Verify parity identities and the degree-2 determinant factorization.
- Gate: stdlib tests pass with no optional numerical dependency.

### WO-RH-02 — scalar interval verifier
- Implement prime-power cell splitting for arbitrary bounded L intervals.
- Use `W00'' = 2(r^3-r-1)/(sqrt(r)(r^2-1))`, r=e^L.
- Verify downward derivative jump `-2 Lambda(q)/sqrt(q)` at q=p^k.
- Certify at most one interior scalar minimizer per cell.
- Reproduce the `[log 3, log 4]` scalar certificate before importing any lower bound as E1.

### WO-RH-03 — f1 audit
- Independently derive `K_q1q1` and its sign threshold.
- Reproduce the midpoint-reflection identity and odd pivot.
- Add a normalization/sign regression test.

### WO-RH-04 — even {0,2} block
- Use bubble basis `b=x(L-x)` rather than raw monomials for numerical work.
- Reproduce `K00`, `K0b`, `Kbb` and individual prime-kernel determinant.
- Reproduce the compact real-space degree-2 certificate.

### WO-RH-05 — direct Fourier cross-check
- Cutoff T=84.
- Implement stable entire low-frequency forms for H0/Hb; never divide naively near t=0.
- Integrate Taylor jets in L directly, not finite differences.
- Target interval statements:
  * `E2,84'' > 0` on `[log 3, 1.20]`;
  * `E2,84' > 0` on `[1.20, log 4]`.
- Combine with one interval point ball near `L=1.1059498113`.
- Only after uniform coverage may the direct-Fourier certificate become E1/SOUND.

### WO-RH-06 — Atlas certificate + PIR bridge
- Emit stable JSON into `math/rh_weil/certificates/` and optionally mirror a headline certificate into root `certificates/` only after review.
- Add content hash and explicit claim scope.
- Lower exact identities as E0/SOUND facts and interval results as E1/SOUND facts.
- Imported transcript values must never masquerade as regenerated Atlas evidence.

### WO-RH-07 — CI
- Add RH tests to a dedicated runner first. Do not silently expand `ci/run_all_certified.py` until certificates are reproducible on the target CI image.
- CI must fail on certificate degradation, missing coverage cells, sign-convention change, or evidence-level promotion without regenerated evidence.

### WO-RH-08 — degree 3
Only after WO-RH-05 closes, start odd block `{q1, b3}` with `b3=x(L-x)(x-L/2)` and endpoint-jet tail bounds.

## Acceptance criteria

A coding agent is finished only when:
- exact tests pass;
- every E1 number is regenerated locally;
- certificate schema validates;
- no imported claim is silently promoted;
- README claim boundary remains intact;
- `git diff` contains no unrelated changes.

## External cross-validation extension — Connes-CvS

After WO-RH-05 is stable, execute the optional external work order at
`external/AGENT_WORK_ORDER_CONNES_CVS.md`.  This does not replace any Atlas-native
certificate and is not a required runtime dependency.
