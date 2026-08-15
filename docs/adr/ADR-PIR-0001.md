# ADR-PIR-0001 — PIR v0.1 substrate: architecture decisions and rejected alternatives

- **Status:** Accepted (Stage 1 substrate)
- **Date:** 2026-07-16
- **Context:** AGENT_WORK_ORDER_v2 (PIR Foundation, Stage 1), adjudicated from
  the ChatGPT "Decompilation of Physics" package. Provenance package docs live
  under `docs/pir-foundation/` (`FINAL_ANALYSIS.md`, `ADJUDICATION.md`,
  `AGENT_WORK_ORDER_v2.md`, `architecture.yaml`). See
  `docs/pir-specification-v0.1.md` for the normative substrate contract and
  `docs/verifier-ops-v0.1.md` for the frozen verifier-op set.
- **Scope:** substrate only — schemas, models, provenance, pass registry,
  verifier-op reference, CI. No new science, no atlas edits, no changes to
  existing B1–B12 verdicts or certificates.

## Decisions

### D1 — Append-only, content-addressed fact store
Facts are immutable and identified by the hash of their canonical content +
provenance skeleton. Re-adding an identical fact is idempotent; re-adding a
different payload under an existing ID is an `AppendOnlyViolation`. Invalidation
never deletes — it appends a `DowngradeRecord` and swaps the stored object for a
`DOWNGRADED` copy while the prior status survives in the downgrade log.
*Consequence:* every historical claim remains auditable; conflicts are stored as
conflicts (hard constraint #1).

### D2 — Orthogonal L and E axes on every fact
`pir_level` (L0–L3, representation abstraction) and `evidence_level` (E0–E4,
warrant strength) are independent enums. Nothing in the models or schema derives
one from the other; `test_pir_models.py::m5` exercises four L×E combinations
(ADJUDICATION #4).

### D3 — Verdict vocabulary locked to the specification
The `verdict` field accepts only the SPEC §4/§6 vocabulary
(`FORCED, PERMITTED, REJECTED, NONIDENTIFIABLE, OBSERVATIONALLY_EQUIVALENT,
APPARATUS_LIMITED, REPRESENTATION_DEPENDENT, AMBIGUOUS`). ChatGPT's candidate
labels (`GLOBAL_CANDIDATE`, …) are a hypothesis-store taxonomy on
`hypothesis.candidate_class` and raise if placed in a `verdict` field
(ADJUDICATION #24; hard constraint #3).
**Deviation from the starter draft, recorded here:** the package's `fact.schema`
draft listed seven verdicts; SPEC §4 additionally makes `AMBIGUOUS` a first-class
verdict (a classification tie, distinct from `NONIDENTIFIABLE`). Per the
work-order instruction to adapt drafts to repository conventions, the lock
tracks the normative SPEC superset (eight), not the draft (seven).

### D4 — SOUND/HEURISTIC honesty is enforced at construction and at emission
A SOUND pass may assert only E0–E2 facts (what it can certify); a HEURISTIC pass
may not assert E0, and any E3/E4 fact it asserts must carry a non-empty, located
`warnings[]` (the Ghidra WARNING-comment pattern). The rule is enforced both in
`Fact.__post_init__` (so a hand-built fact cannot bypass it) and in
`PassRegistry.check_emission` (defense in depth) — hard constraint #4.

### D5 — Namespace promotion requires a typed transform
Namespaces carry an abstraction rank; a fact depending on a fact in a
higher-warrant namespace (e.g. `effective:` → `global:`) requires an explicit,
typed `NamespaceTransform`, or it raises `IllegalNamespacePromotion`. A transform
whose declared endpoints do not match the edge is also rejected (hard
constraint #6). Stage 1 checks that the transform is *declared and typed*; the
*checkable* type system is Stage 2's symbolic bridge.

### D6 — Measurement interface mandatory on DOMAIN/MEASUREMENT facts
Per SPEC §2, a DOMAIN- or MEASUREMENT-layer fact must name its declared
measurement interface (`measurement_interface`: calibration-route / apparatus
IDs). UNIVERSAL structural facts are exempt. This is the model-level home of the
required "missing measurement interface on a DOMAIN-layer fact → reject" test.
**Deviation recorded:** `measurement_interface` is added to the fact schema
beyond the starter draft, which grounded 𝖬 only implicitly via `source_spans`;
making it an explicit, checkable field is the repository-convention adaptation
(certificates already require `calibration_route` at ≥ v0.3).

### D7 — Canonical JSON is hash-stable; exact rationals are strings
Sorted keys, shortest-round-trip floats, non-finite floats rejected, and
`Fraction` serialized as `"p/q"` in lowest terms — never a float (hard
constraint #9). This is what makes content-addressing reproducible across
machines.

### D8 — CI fails on *degradation*, tolerant to representation jitter
`ci/run_all_certified.py` reruns every `tests/test_b*.py`, captures regenerated
certificates to a build directory, restores the committed copies (leaving the
tree unchanged), and compares with a degradation predicate: certification flags,
status/verdict labels, and result-key presence are strict; timestamps, added
keys, float jitter within tolerance, and numeric-representation differences
(e.g. `-1.8` vs `-1.8+0.0j`) are not degradations. It writes a self-hash–signed
run manifest with pinned seed and library versions and exits nonzero on any
failure or degradation (hard constraint #11).

## Rejected alternatives

- **Mutable fact store** — rejected (D1). A store that overwrites cannot support
  invalidation-without-deletion or an audit trail; conflicts would be lost.
- **Single conflated L/E vocabulary** — rejected (D2). Collapsing representation
  abstraction and warrant strength into one axis is the exact confusion the
  program's three-layer rule exists to prevent.
- **Candidate labels as verdicts** — rejected (D3). Emitting `GLOBAL_CANDIDATE`
  as an outcome would let a hypothesis-store taxonomy masquerade as a certified
  verdict.
- **OpenTelemetry-style provenance (ChatGPT S24)** — rejected. Distributed-trace
  spans model service calls, not epistemic derivation with assumption taint and
  invalidation traversal; we use a PROV-O core (Entity/Activity/Agent) instead.
- **Frida-style dynamic instrumentation (ChatGPT S9)** — rejected as out of
  scope for a data-first substrate; not needed to lower certified benchmarks.
- **Pydantic as the model layer** — rejected for Stage 1. The work order permits
  Pydantic only "if already acceptable in repo"; it is not an existing
  dependency, so the substrate uses stdlib dataclasses and a vendored,
  auditable JSON-Schema-subset validator (`pir/jsonschema_mini.py`). Result:
  tests and CI run on a clean checkout with **no new required dependency**; the
  real `jsonschema` package, when present, is used only as a cross-check.
- **Regenerating committed certificates in place during CI** — rejected (D8).
  CI must never mutate the committed scientific record; it diffs against a
  snapshot and restores it.

## Consequences
The substrate is import-safe alongside B1–B12 (nothing under `pir/`, `ci/`, or
the new tests imports or mutates benchmark code) and adds no required runtime
dependency. Stage 2 (lowering B9, analyzer runtime, symbolic bridge) can build
on these contracts without a schema migration: `fact.witness` and
`fact.impossibility_certificate` are already reserved (hard constraint #13).
