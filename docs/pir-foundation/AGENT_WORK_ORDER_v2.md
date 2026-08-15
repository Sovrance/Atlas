# AGENT WORK ORDER v2 — PIR Foundation (Stage 1)
**Supersedes:** ChatGPT REMOTE_AGENT_WORK_ORDER.md (adjudicated ACCEPT-MODIFIED; see ADJUDICATION.md)
**Repository:** github.com/ErickLGonzalez/Global-variables
**Branch:** `feature/pir-foundation-v0.1` — commit in coherent increments; **do not merge**.

## Mission
Implement Stage 1 (Substrate) of the adjudicated Physics Reverse-Engineering Architecture: PIR v0.1 schemas and models, provenance with invalidation traversal, the verifier-op reference, and the CI headless entrypoint — **without changing any existing scientific verdict, certificate, or atlas cell.**

## Required reading (in order)
1. `docs/SPECIFICATION.md` (normative; verdict vocabulary §6 is LAW for this work)
2. `README.md`
3. `docs/constant-atlas-v0.6.md`
4. `docs/roadmap-v2.md`
5. `FINAL_ANALYSIS.md` (this package) — esp. §1 (L×E orthogonality), §3 (fact store), §7 (Stage 1)
6. `ADJUDICATION.md` (this package) — esp. verdicts #4, #16, #24

## Create
```
docs/adr/ADR-PIR-0001.md
docs/pir-specification-v0.1.md
docs/verifier-ops-v0.1.md            # NEW vs v1 work order: freeze existing opcode set
pir/__init__.py
pir/models.py
pir/types.py
pir/namespaces.py
pir/provenance.py
pir/canonical.py
pir/passes.py                        # NEW: pass registry with SOUND/HEURISTIC tag + warnings[]
pir/schema/artifact.schema.json
pir/schema/event.schema.json
pir/schema/fact.schema.json
pir/schema/hypothesis.schema.json
pir/schema/intervention.schema.json
pir/schema/provenance.schema.json
tests/test_pir_models.py
tests/test_pir_provenance.py
tests/test_pir_namespaces.py
tests/test_pir_passes.py             # NEW
examples/pir/minimal_circuit.json
ci/run_all_certified.py              # NEW: headless CI entrypoint
```
Starter schema drafts are provided in this package under `schemas/` and a starter example under `examples/` — treat them as reference drafts to adapt to repository conventions, not as final.

## Hard constraints (adjudicated)
1. **Append-only facts.** No pass may mutate or delete another pass's facts. Conflicts are stored as conflicts.
2. **Dual coordinates on every fact:** `pir_level` (L0–L3, representation abstraction) AND `evidence_level` (E0–E4, warrant strength). These are orthogonal; schema must not derive one from the other.
3. **Verdict vocabulary is locked to SPEC §6.** The fields `verdict`/`outcome` accept ONLY the SPEC vocabulary. ChatGPT's candidate labels (GLOBAL_CANDIDATE, TOPOLOGICAL_CANDIDATE, HIDDEN_COMMON_CAUSE_CANDIDATE, REPRESENTATION_ARTIFACT, NOT_DETECTED) go in a separate `candidate_class` tag field on hypothesis objects only. Add a negative test asserting a non-SPEC verdict is rejected.
4. **SOUND/HEURISTIC tag mandatory** on every registered pass; HEURISTIC passes must support a located `warnings[]` output. SOUND passes may only assert facts they can certify.
5. **Every derived fact carries:** analyzer ID + version, dependency fact IDs, assumption IDs (taint), source spans, evidence level, pir level, layer (UNIVERSAL/DOMAIN/MEASUREMENT), namespace, status.
6. **Namespaces:** `raw:, apparatus:, operational:, domain:, latent:, gauge:, invariant:, global:, effective:, analyst:`. Cross-namespace references require an explicit, typed transform record. Illegal promotion (e.g., `effective:` → `global:` without transform) must raise and is a required negative test.
7. **Raw evidence (L0) cannot contain inferred universal claims.** Loader-layer validation enforces this.
8. **Invalidation traversal:** given an assumption ID, return all transitively dependent facts/claims and mark them `DOWNGRADED(reason)` — without deleting anything (append a downgrade record).
9. **Canonical JSON must be hash-stable** (sorted keys, fixed float handling, exact Fractions serialized as strings `"p/q"` — never floats).
10. **Do not change existing B1–B12 code** except import-safe compatibility fixes. **No new atlas promotions. No atlas edits.**
11. **CI entrypoint** (`ci/run_all_certified.py`): discovers and reruns all `tests/test_b*.py`, regenerates certificates to a build directory, diffs against committed certificates, and **exits nonzero on any degradation or failure-to-certify**. Pin seeds and record solver/library versions in a signed run manifest (name, version, content hashes).
12. **Dependencies:** standard library + Pydantic (if already acceptable in repo) + jsonschema. No services, no databases, no UI.
13. **Certificate fields for future symbolic bridge:** fact schema reserves `witness` and `impossibility_certificate` (UNSAT-core-style) optional fields, so Stage 2 needs no schema migration.
14. **Security note (downgraded from threat model):** parsers for external artifacts must not execute code from data; treat all paper/dataset inputs as untrusted strings/bytes.

## Negative tests (required)
- Provenance cycle detection (A depends on B depends on A → reject).
- Missing measurement interface on a DOMAIN-layer fact → reject.
- Illegal namespace promotion without transform record → reject.
- Non-SPEC verdict string → reject.
- HEURISTIC pass asserting `evidence_level: E0` → reject.
- Mutation of an existing fact → reject (append-only).

## Definition of done
- Existing test suite passes, unmodified.
- All new tests (positive + the six negative tests above) pass.
- All six schemas validate `examples/pir/minimal_circuit.json`.
- `ci/run_all_certified.py` runs green on the current repo and demonstrably fails on a synthetically corrupted certificate (include that test).
- `docs/verifier-ops-v0.1.md` documents each existing verifier op (Schur pivot, rank, flat extension, Jacobian-rank identifiability, CPTP/symplectic gates) with signature, certificate format, and the benchmarks that use it.
- README gains a short section: PIR is an evidence substrate, not a proof engine; the atlas engine never reads raw datasets directly.
- ADR-PIR-0001 records rejected alternatives (including: mutable store, single-vocabulary L/E conflation, OpenTelemetry provenance — all rejected in adjudication).
- Machine-readable manifest lists created files with hashes.
- Sprint deliverable zipped for review at `/mnt/user-data/outputs/`, **each zip < 25 MB** (split part1/part2/... if larger), per project workflow. Erick integrates and pushes; agent never pushes to remote.
