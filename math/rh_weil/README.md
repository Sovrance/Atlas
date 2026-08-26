# RH / Weil finite-compression verifier

Reproducible source behind the Constant Atlas positivity-verifier example
(`G00 > 0`, Schur pivot positivity), grown into a program of its own.

Live agent instructions: [`AGENT_INSTRUCTIONS.md`](AGENT_INSTRUCTIONS.md).
Machine-readable status:
[`certificates/work_order_status.json`](certificates/work_order_status.json) —
when this file disagrees with that one, that one wins and this one is a bug.
[`scripts/check_docs.py`](scripts/check_docs.py) fails CI on exactly that
disagreement.

---

## 1. Scope and claim boundary

This program studies **finite-dimensional polynomial compressions of the
localized Weil quadratic form**. **It does not prove the Riemann Hypothesis.**
A finite block certificate is evidence for that stated block, interval,
normalization and cutoff, and for nothing else. Every artifact carries
`rh_proof_claim: false` and `claim_scope: finite_dimensional_weil_compression`.

Accepted assembly convention: `G = G0 - Gp + Ginf`.
Current research cell: `L ∈ [log 3, log 4]`.

## 2. Active normalization

The pole block is **adjudicated** (WO-RH-17) and frozen. The adopted form comes
from the explicit formula,

`G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+`,
`E_i^± = ∫_0^L h_i(x) e^{±x/2} dx`,

and the legacy even block `(√3/2)(v₊v₊ᵀ + v₋v₋ᵀ)` is **rejected**: it equals the
adopted pole times `(√3/2)cosh(L/2)`, a factor equal to 1 only at `L = log 3` —
a calibration fitted at one test point (+8.25 % at `L = log 4`). The odd pivot
`−8A²` was already correct and is unchanged.

`src/pole.py` is the single implementation and every production path routes
through it. The rejected candidate is archival in `src/rejected_pole.py`;
`tests/test_production_imports.py` fails CI if anything under `src/` imports it
or spells its scale in executable code.

Active normalization id: `norm_sha256_f84b2fae2e13c777b1f829ef2567699c`. Every
rigorous certificate binds it, and the promotion predicate refuses a certificate
whose id has drifted.

Certificates that depended on the rejected block are
`QUARANTINED_NORMALIZATION_ADJUDICATION` — preserved, not deleted, not
relabelled. The quarantine is enforced at the point of write, so re-running a
legacy certifier cannot erase it. They must be **regenerated**, never
reinterpreted. See
[`docs/NORMALIZATION_ADJUDICATION_v0.1.md`](docs/NORMALIZATION_ADJUDICATION_v0.1.md).

## 3. Current certified results

Numbers are read from the certificates by
[`scripts/check_docs.py`](scripts/check_docs.py), not from prose memory; a value
here that disagrees with its certificate fails the gate.

| Object | Domain | Result | Warrant | Certificate |
|---|---|---|---|---|
| scalar `G00` | `[log 3, log 4]` | `≥ 0.06962397439120689` | E1 | `e1_scalar_log3_log4.json` |
| degree-1 odd `O1` | cell | `≥ 0.015026786943232338` | E1 | `e1_degree1_log3_log4.json` |
| compact degree-2 even `E2` | cell | `≥ 2.0652586666890377e-06` | E1 | `e1_degree2_compact_log3_log4.json` |
| `T=84` degree-2 `E2,84`, points | `log 3`, `1.20`, `log 4` | positive at each | E1 | `e1_fourier_T84_points.json` |
| `T=84` interior minimum | cell | `L*` isolated to `[1.1059498108971377, 1.1059498114329873]`, a strict interior minimum | E1 | `e1_fourier_T84_interior_minimum.json` |
| `T=84` degree-2 uniform | cell | `E2,84(L) ≥ 3.4251152511218656e-06 > 0` | E1 | `e1_fourier_T84_uniform_degree2.json` |
| odd degree-3 block | cell | positive definite, inertia `(2,0,0)`, one stratum, no transition regions; `O1 ≥ 0.015331267702040234`, `det ≥ 1.073120529992708e-06` | E1 | `e1_degree3_odd_positivity_log3_log4.json` |
| rank–trace | 3 degree-3 sample points | `rank ≥ 1` — nontrivial but weak, against a true rank of 2 | E1 | `e1_degree3_odd_moments_log3_log4.json` |
| spectral moments `m₁..m₄` | 3 degree-3 sample points | dimension 2; mixed conclusive / insufficient B1 queries | E1 | `e1_degree3_odd_moments_log3_log4.json` |
| **3×3 even block `{1, b, b²}`** | cell | **positive definite, inertia `(3,0,0)`**, one stratum, no transition regions; `Δ1 ≥ 0.07537591825740127`, `Δ2 ≥ 3.4244304067666463e-06`, `Δ3 ≥ 6.451586222238981e-15` | E1 | `e1_degree4_even3_positivity_log3_log4.json` |
| 3×3 even block, inertia route | cell | `(3,0,0)` by interval LDL* congruence, 2838 boxes, max depth 3 | E1 | `e1_degree4_even3_inertia_log3_log4.json` |
| 3×3 even moments `m₁..m₄` | 5 sample points | moments do **not** force the inertia at `n = 3`; rank–trace `≥ 1` against a true rank of 3 | E1 | `e1_degree4_even3_moments_log3_log4.json` |
| finite theorem boundary | — | 21 theorems proved in Lean 4 / Mathlib, no `sorry`, three standard axioms | FORMAL (implication only) | `formal_theorem_certificate.json` |
| 3×3 even independent assembly | 5 points | agrees with the rigorous assembly to `1.5637e-13` | **E3 — regression only, never a warrant** | `e3_degree4_even3_crosscheck.json` |
| **reference metric `M`** | every `L > 0` | the exact `L²` Gram of the basis is positive definite — exact rational Sylvester + diagonal congruence | **E0** | `e0_eng009_reference_metric.json` |
| **generalized gap, scalar** | cell | `λmin(G, M) ∈ [0.052666425704956055, 0.058518285644676]` | E1 | `e1_eng009_generalized_gap_log3_log4.json` |
| generalized gap, odd `{q1}` | cell | `λmin ∈ [0.12558043003082275, 0.1395338743236466]` | E1 | same file |
| generalized gap, even `{1, b}` | cell | `λmin ∈ [0.0005867481231689453, 0.0006520075621795992]` | E1 | same file |
| generalized gap, odd `{q1, b3}` | cell | `λmin ∈ [0.004970192909240723, 0.005522444641828629]` | E1 | same file |
| **generalized gap, even `{1, b, b²}`** | cell | `λmin ∈ [3.606081008911133e-05, 4.007732435284577e-05]` | E1 | same file |
| structural dataset + verdict | 5 blocks | determinant collapse is `MOSTLY_COORDINATE_DRIVEN_BUT_THE_GAP_ALSO_DECAYS` | E1 | `eng009_structural_dataset.json` |
| scaling models | 2 families | exponential vs power-law fits with explicit n = 4 falsifiers | **E3 — plans, never promoted** | `e3_eng009_scaling_models.json` |
| ENG-010 preview, even 4×4 | midpoint | `{1, b, b², b³}` float minors and scouted gap — superseded by the certified row below, and off by 20× at the bottleneck, which is why previews never promote | **E3 — never a warrant** | `e3_eng010_even4_preview.json` |
| **4×4 even block `{1, b, b², b³}`** | cell | **positive definite, inertia `(4,0,0)`**, one stratum, by LDL* (26,854 boxes) *and* Sylvester; `Δ4 ≥ 5.2788569975917046e-24` raw | E1 | `e1_degree6_even4_positivity_log3_log4.json` |
| **generalized gap, even 4×4** | cell | `λmin(G, M) ∈ [1.9073486328125e-06, 2.415977410246834e-06]`, bottleneck at `L → log 4` | E1 | `e1_eng010_even4_generalized_gap_log3_log4.json` |
| scaling-model adjudication | n = 4 | **`NEITHER_FALSIFIED`** — the enclosure lands in the overlap of both preregistered ×5 windows, between the two point predictions | E1 adjudicating E3 | `eng010_scaling_model_adjudication.json` |
| 4×4 independent assembly | 7 points | agrees with the rigorous assembly to `1.243e-13` | **E3 — regression only, never a warrant** | `e3_degree6_even4_crosscheck.json` |

The 3×3 block is the first here where the determinant does not fix the spectrum:
`det > 0` on a 3×3 is consistent with `(3,0,0)` and with `(1,2,0)`, so
determinant-only reporting could not have told a positive definite block from one
with a two-dimensional negative subspace. Its preconditioner is a diagonal matrix
of exact powers of two, frozen for the cell, which is a congruence performed
without rounding — and the ENG-007 congruence theorems are what say the inertia
it reports is the block's own.

**Why raw determinant magnitude is not a basis-invariant distance-to-RH
signal.** A determinant moves by `det(S)²` under any change of basis `S`, so
its size measures the coordinates as much as the matrix. Concretely: the
certified 3×3 raw determinant bound is `6.45e-15`, which looks like a spectrum
about to fail — but the determinant of the exact reference metric `M` (a Hankel
moment matrix, positive definite for *every* `L > 0`, minors exactly `1`,
`1/180`, `1/7938000`) collapses by seven orders of magnitude over the same
ladder while being the healthiest object in the program. Of the thirteen orders
the raw determinant loses from the scalar block to the 3×3, about seven are the
coordinates. The object that does not move under simultaneous congruence is the
generalized eigenvalue of the pencil `(G, M)` — `det(SᵀGS − λSᵀMS) =
det(S)²·det(G − λM)`, checked exactly in the tests and proved in Lean — and
*its* certified decay over the same span is about three orders, slow and so far
unclassified. Margins in this program are therefore quoted against a named
metric or not at all; see
[`docs/GENERALIZED_GAP_ENG009_v0.1.md`](docs/GENERALIZED_GAP_ENG009_v0.1.md).

Two of the older rows deserve their qualifiers. The rank–trace bound is *weak*: it
proves `rank ≥ 1` where the rank is 2, and saying so is the point — a bound that
is true and uninformative is still a result, and pretending otherwise is how a
program starts believing its own machinery. And the B1 moment queries come back
**INSUFFICIENT_INFORMATION** for "do the moments force PSD" even though the block
*is* positive definite, because PSD-ness of a truncated localizing matrix is
necessary and not sufficient; only the refuting direction is available from four
moments.

### Work-order status

| Order | State |
|---|---|
| WO-RH-01…04 | exact identities, f1 audit, even block — E0 |
| WO-RH-05 | **recovered by ENG-005** — cutoff-free uniform E1 |
| WO-RH-08 | **done by ENG-006** — odd degree-3 implemented and certified |
| WO-RH-09 | **recovered by ENG-004** — scalar canary PROMOTED |
| WO-RH-15 | **recovered by ENG-005** — T=84 uniform E1 + interior minimum |
| WO-RH-17/18 | normalization adjudicated, Candidate A adopted, three-way internal cross-check |
| WO-RH-28…33 | inertia, rank–trace, moment engines; degree-3 exact block and E3 scan |
| WO-RH-34 | degree-3 E1 certificate — positive definite |
| WO-RH-36 | positivity-vs-inertia/moment information report |
| WO-RH-38 | pinned Lean project and theorem boundary |
| WO-RH-43 | statement comparator, axiom audit, theorem manifest |
| WO-RH-46 | 3×3 even pilot prepared — E0 identities and an E3 preview only |
| WO-RH-47 | basis `{1, b, b²}` frozen; dyadic preconditioner congruence certified |
| WO-RH-48 | six exact entries; independent assembly agrees to `1.56e-13` |
| WO-RH-49 | derivative provider generalized; prior certified bounds unchanged |
| WO-RH-51 | **3×3 even block certified positive definite, inertia `(3,0,0)`** |
| WO-RH-53 | Lean 3×3 certificate replay |
| WO-RH-55 | cross-block diagnostics prepared for ENG-009 |
| WO-RH-58 | reference metric E0; **generalized gap certified for all five blocks** |
| WO-RH-59 | determinant-collapse verdict: mostly coordinates, gap decays slowly |
| WO-RH-62 | ENG-010 target selected: even 4×4; `bcube` E0-prepared |
| WO-RH-63 | Lean generalized-gap implication |
| WO-RH-69 | **4×4 even block certified positive definite, inertia `(4,0,0)`, two routes** |
| WO-RH-70 | 4×4 generalized gap enclosed; bottleneck found at the right cell edge |
| WO-RH-71 | preregistered models adjudicated `NEITHER_FALSIFIED`, before any refit |
| WO-RH-75 | ENG-011 target selected: even 5×5, refits diverge 11.9× there — 21 Lean theorems, no `sorry` |

The authoritative per-order state, including the pre-quarantine values WO-RH-17
forbids deleting, is
[`certificates/work_order_status.json`](certificates/work_order_status.json).

## 4. Verifier architecture

```
src/pole.py            the one Candidate-A pole primitive
src/weil_entries.py    prime kernels, frequency-space archimedean route, assembly
src/archimedean_realspace.py
                       the exact real-space archimedean form and its L-jets
src/interval_cover.py  the shared adaptive interval branch-and-bound
src/promotion.py       the one promotion predicate; source hashes; normalization binding
inertia/               interval Hermitian LDL congruence, signature, stratification
ranktrace/             the rank-trace / Hilbert-Schmidt theorem with enforced hypotheses
moments/               m1..m4 as traces of powers, fed to the Atlas B1 solver
formal/                the Lean 4 project: definitions, statements, proofs, comparator
src/formal_evidence.py the FORMAL warrant and the boundary it may not cross
src/content_kinds.py   the one registry of content kinds and what each licenses
src/basis_algebra.py   exact overlap kernels, derived from the basis coefficients
src/even3.py           the rigorous 3x3 even block and its dyadic preconditioner
src/independent_even3.py
                       a second assembly that imports none of the above
src/pilot3.py          the ENG-007 pilot, superseded by src/even3.py
```

Four facts about this architecture are load-bearing:

* **No eigenvalue solver is reachable from any rigorous path.** A gate in
  `scripts/ci_inertia.py` proves it by import scan. Signatures come from interval
  LDL congruence, cross-checked against Descartes' rule of signs on the
  characteristic polynomial — exact, because a symmetric matrix is real-rooted.
* **mpmath never certifies.** E1 requires python-flint/Arb. mpmath appears in E3
  previews and in the pilot, both labelled.
* **Derivatives are exact jets, never finite differences** in a rigorous path.
  Finite differences appear only as a check on the jets — and earned it: they
  caught a `d²/dL²(L³/6)` coefficient written as `L/2` instead of `L`.
* **Certificates are build outputs.** Each records the source hashes it depends
  on; editing a certifier's dependency makes its certificate stale by
  construction and the chain refuses it.

## 5. Canonical commands

```bash
python3 math/rh_weil/scripts/run_rh_weil_suite.py        # fast path — does NOT re-derive E1
python3 math/rh_weil/scripts/run_rigorous_chain.py --release
python3 math/rh_weil/scripts/ci_inertia.py --gate fast
python3 math/rh_weil/scripts/ci_inertia.py --gate rigorous
python3 math/rh_weil/scripts/check_docs.py               # rh-docs gate
python3 math/rh_weil/scripts/ci_formal.py                # rh-formal gate
python3 math/rh_weil/scripts/certify_even3.py            # the 3x3 even block
python3 math/rh_weil/scripts/report_even3_information.py # §WO-RH-52/55 reports
```

Passing the fast suite does **not** mean the rigorous certificates are current —
it does not re-derive the scalar canary. Read the rigorous chain's exit code
before believing an E1 claim is fresh.

The formal layer, from `math/rh_weil/formal/`:

```bash
lake build                                    # AtlasRH + comparator
lake env lean comparator/PrintAxioms.lean     # statement comparator + axiom report
python3 ../scripts/check_formal_manifest.py   # manifest gate (Lean layer optional)
```

Rigorous dependencies are required, not optional extras — a missing one fails
the job rather than degrading an artifact to weaker evidence:

```bash
pip install -r math/rh_weil/requirements-rigorous.txt
```

## 6. Evidence and certificate taxonomy

| Class | Meaning | May it warrant a claim? |
|---|---|---|
| **E0 / SOUND** | exact algebraic identity, independently re-derived | yes |
| **E1 / SOUND** | interval-certified, python-flint/Arb, outward rounded, hashes bound | yes |
| **E3 / HEURISTIC** | floating scan, topology preview, conditioning report | never |
| **FORMAL** | a machine-checked *implication*, Lean 4 / Mathlib, pinned toolchain | only the implication |

`FORMAL` is deliberately not a rung of the E-ladder: the ladder grades how
reliable a **number** is, and `FORMAL` grades whether an **implication** was
checked. The two travel in separate fields, `numeric_warrant` and
`logical_implication_warrant`, on every PIR fact. So the degree-3 result reads:

```
Arb interval enclosure                E1     the bounds hold
Lean: positive bounds => PD           FORMAL the bounds suffice
```

Remove either half and the claim is gone. A formal theorem may strengthen an
exact theorem dependency; it **never** converts interval numerical evidence to
FORMAL.

Content kinds are also load-bearing. An **inertia** certificate never satisfies
a consumer requiring PSD, even when its signature is `(2,0,0)` — "I know the
signature" must not be read as "it is positive". The degree-3 artifact answers
such a consumer as a *positivity* certificate carrying certified bounds, with
the inertia object nested inside it still refusing.

Imported notebook claims stay `IMPORTED_PENDING_REGENERATION` until regenerated
in-repo. `uncertified != pass`.

## 7. External cross-validation

**Connes / CvS** (`external/`) — an optional external oracle for shared Weil
ingredients. It quantifies no projection or truncation error, so it reports
`NOT_COMPARABLE` and **never certifies**. Do not compare Galerkin matrix entries
to Atlas polynomial Gram blocks by index. See
[`external/CONNES_CVS_MAPPING.md`](external/CONNES_CVS_MAPPING.md) and
[`external/PROVENANCE.md`](external/PROVENANCE.md).

```bash
pip install 'connes-cvs==0.3.1' python-flint mpmath
python3 math/rh_weil/scripts/run_connes_cvs_crosschecks.py
```

The internal cross-check is separate and is **three-way** — explicit formula,
compact real space, direct Fourier (`certificates/normalization_crosscheck.json`).

**zeta-23-lean** (`external/zeta23/`) — architecture reference for the formal
layer, pinned at commit `cec57f9`, Apache-2.0. Status is `REFERENCE_ONLY`:
nothing is vendored, imported or depended on, and its toolchain does not
currently compose with this project's. Its mapping records one thing worth
knowing: the general rank–trace inequality Atlas carries as *unproved* is proved
upstream as `Zeta23.RHLinalg.rank_trace_ineq_two`. See
[`external/zeta23/MAPPING.md`](external/zeta23/MAPPING.md).

## 8. Current frontier — ATLAS-RH-ENG-010, the first prediction test

ENG-009 made the generalized gap the primary margin and preregistered two
decay models that disagreed by ~8× about dimension 4. **ENG-010 certified the
4×4 even block and adjudicated them**:

> `G[{1, b, b², b³}](L)` is positive definite for every `L ∈ [log 3, log 4]`,
> inertia `(4, 0, 0)` by two independent routes, and its generalized gap
> satisfies `λmin(G, M) ∈ [1.9073486328125e-06, 2.415977410246834e-06]` —
> with the bottleneck at the right cell edge, twenty times below the E3
> midpoint scout.

The adjudication verdict, recorded against the preregistered artifact before
any refit, is **`NEITHER_FALSIFIED`**: the certified enclosure is tight, but
it lands in the overlap of both models' ×5 falsifier windows, between the two
point predictions. The experiment was decisive about the number and
indecisive between the models — expected outcome E of the work order, and a
statement about the preregistered tolerance rather than about the block. See
[`docs/EVEN4_GAP_ADJUDICATION_ENG010_v0.1.md`](docs/EVEN4_GAP_ADJUDICATION_ENG010_v0.1.md).

**Next: ENG-011.** The refitted models (four certified points) diverge by
11.9× at n = 5, so `eng011_target_selection.json` selects the even 5×5
`{1, b, b², b³, b⁴}` — plus the new structural question the certification
surfaced: why does the pencil weaken toward `log 4`? No E1 work on n = 5 is
launched before ENG-010 is interpreted.

## 9. History

Work-order records, each accurate for the work it describes:

* [`docs/NORMALIZATION_ADJUDICATION_v0.1.md`](docs/NORMALIZATION_ADJUDICATION_v0.1.md) — WO-RH-17/18, the pole adjudication (still in force)
* [`docs/SCALAR_E1_CANARY_ENG004_v0.1.md`](docs/SCALAR_E1_CANARY_ENG004_v0.1.md) — ENG-004
* [`docs/CORE_E1_RECOVERY_ENG005_v0.1.md`](docs/CORE_E1_RECOVERY_ENG005_v0.1.md) — ENG-005
* [`docs/INERTIA_RANKTRACE_MOMENTS_ENG006_v0.1.md`](docs/INERTIA_RANKTRACE_MOMENTS_ENG006_v0.1.md) — ENG-006
* [`docs/FORMAL_BOUNDARY_ENG007_v0.1.md`](docs/FORMAL_BOUNDARY_ENG007_v0.1.md) — ENG-007
* [`docs/HIGHER_DIMENSIONAL_BLOCK_ENG008_v0.1.md`](docs/HIGHER_DIMENSIONAL_BLOCK_ENG008_v0.1.md) — ENG-008

Superseded instructions, preserved and labelled rather than deleted:

* [`docs/history/agent-instructions-initial-integration.md`](docs/history/agent-instructions-initial-integration.md) — the original integration work order
* [`docs/ATLAS_RH_ENG_002_Mathematical_Parity_Run18.md`](docs/ATLAS_RH_ENG_002_Mathematical_Parity_Run18.md), [`docs/README.md`](docs/README.md), [`docs/AGENT_EXECUTION_CHECKLIST.md`](docs/AGENT_EXECUTION_CHECKLIST.md) — ENG-002
* `certificates/history/` — superseded certificates, kept as provenance
* [`notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`](notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md) — the imported notebook checkpoint and what became of each of its claims
