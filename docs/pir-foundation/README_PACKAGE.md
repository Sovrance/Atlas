# PIR Foundation Package v2 — Post-Adjudication
**Global-variables project · Rosy 🌹 · 2026-07-15**

Merged deliverable from cross-agent review: Rosy deep-research report + ChatGPT "Reverse-Engineering Architecture" package, adjudicated per project cross-agent review practice.

## Contents

| File | Purpose |
|---|---|
| `FINAL_ANALYSIS.md` | The consolidated, post-adjudication analysis. **Read this first.** Also the document to carry into the MnemesisOS chat (§8 flags the transferable patterns). |
| `ADJUDICATION.md` | Formal accept / accept-modified / reject / defer verdicts on all 26 ChatGPT patterns + package-level items, plus the Rosy additions merged in. |
| `AGENT_WORK_ORDER_v2.md` | Instructions for the coding agent — Stage 1 (PIR substrate) scope, hard constraints, negative tests, definition of done. Supersedes ChatGPT's work order. |
| `architecture.yaml` | Target architecture v0.2 (adjudicated: SPEC-§6 verdict lock, L×E orthogonality, SOUND/HEURISTIC tags, dual diff scores, OED reframing, UI deferred). |
| `engineering_backlog.csv` | Stage-gated backlog (S1 substrate → S2 lowering/analysis → S3 forward loop; P9 UI deferred). |
| `schemas/*.schema.json` | Six starter JSON-schema **reference drafts** (artifact, event, fact, hypothesis, intervention, provenance) encoding the adjudicated constraints. The agent adapts these to repo conventions. |
| `examples/minimal_circuit.json` | End-to-end minimal example: B9-style circuit lowered to L0/L1, one exact L2 fact with Schur-pivot witness, two hypotheses retained as OBSERVATIONALLY_EQUIVALENT, one declared discriminating intervention. |

## Integration workflow (standard)
1. Erick reviews `FINAL_ANALYSIS.md` and `ADJUDICATION.md`.
2. Hand `AGENT_WORK_ORDER_v2.md` + `schemas/` + `examples/` to the coding agent with repo access; branch `feature/pir-foundation-v0.1`; no merge, no atlas edits.
3. Agent returns Stage 1 deliverables zipped (<25 MB per zip) for verification; Erick integrates and confirms with commit hash.

## Non-negotiables carried through every file
- Verdict vocabulary locked to SPEC §6; candidate-class labels are hypothesis-store taxonomy only.
- PIR L-levels ⊥ evidence E-levels (both on every fact).
- Append-only store; assumption-taint with invalidation traversal; SOUND/HEURISTIC tagging with located warnings.
- Similarity and confidence reported separately, correlator named.
- Existing B1–B12 verdicts, certificates, and atlas untouched.
