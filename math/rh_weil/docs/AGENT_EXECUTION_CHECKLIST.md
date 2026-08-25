---
status: HISTORICAL
superseded_by: ATLAS-RH-ENG-007
---

> **HISTORICAL / SUPERSEDED.** This is ATLAS-RH-ENG-002: the execution checklist for that work order. Every box has since been closed, and several of its instructions were later reversed by the WO-RH-17 normalization adjudication.
> It is kept as the record of what was asked and when, not as instruction. For
> the live state read [`../AGENT_INSTRUCTIONS.md`](../AGENT_INSTRUCTIONS.md) and
> [`../certificates/work_order_status.json`](../certificates/work_order_status.json).
> No RH proof claim is made here or anywhere in this program.

# Agent Execution Checklist --- ATLAS-RH-ENG-002

## Baseline

-   [ ] Confirm merge baseline includes PR #5 / 1fda9b4 or descendant.
-   [ ] Run current `math/rh_weil/scripts/run_rh_weil_suite.py`.
-   [ ] Record python, python-flint, Arb/FLINT, mpmath versions.

## E1 regeneration

-   [ ] WO-RH-09 scalar E1.
-   [ ] WO-RH-10 degree-1 E1.
-   [ ] WO-RH-11 compact degree-2 E1.

## True T=84 Fourier

-   [ ] Keep E3 energy probe quarantined.
-   [ ] Implement h_plus.
-   [ ] Implement pole block.
-   [ ] Implement prime block.
-   [ ] Implement rigorous \[0,84\] archimedean integral.
-   [ ] Reproduce three Run-16 point balls.
-   [ ] Implement H0/Hb analytic L-jets.
-   [ ] Certify E2'' \> 0 on \[log3,1.20\].
-   [ ] Certify E2' \> 0 on \[1.20,log4\].
-   [ ] Emit uniform E1 degree-2 Fourier certificate.

## PIR/CI

-   [ ] Emit E0/E1/E3 facts.
-   [ ] Validate PIR schema.
-   [ ] HEURISTIC warning attached to E3 scan.
-   [ ] Rigorous CI job added.
-   [ ] No RH-proof language.

## Final

-   [ ] Update work_order_status.json.
-   [ ] Update README/notebook integration status.
-   [ ] Regenerate SHA256 manifest.
-   [ ] Atlas CI green.
-   [ ] WO-RH-08 unblocked only after uniform T=84 E1.
