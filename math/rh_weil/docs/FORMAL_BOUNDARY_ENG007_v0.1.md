# The Formal Theorem Boundary and the Documentation Truth Pass (ATLAS-RH-ENG-007)

**Baseline:** ENG-006, merged. **Scope:** formalize the stabilized finite
theorem boundary in Lean 4 / Mathlib, add statement-comparator and axiom-audit
discipline, make the documentation mechanically truthful, and prepare the first
genuinely three-dimensional parity block for ENG-008.

**No RH proof claim.** Claim scope is `finite_dimensional_weil_compression`.
Every theorem this work order proves is finite linear algebra over `ℝ`; none is
about the Riemann zeta function.

---

## 1. Why formalize now

Through ENG-006 the program's guarantees were of one kind: a number, produced by
interval arithmetic, with its hypotheses checked at call time. That is a real
guarantee and it is not the only one available. Every certified number reaches a
claim through an *implication* — "the leading entry is bounded below by
`0.0153…` and the determinant by `1.073…e-06`, **therefore** the block is
positive definite" — and until ENG-007 that implication was carried by a Python
function and a comment.

So ENG-007 splits the warrant in two. The numbers keep the warrant they had. The
step from the numbers to the claim gets a separate one, and that step is now
machine-checked against a pinned Lean toolchain.

The distinction is enforced rather than described. `src/formal_evidence.py`
takes the numeric warrant as an argument and returns it unchanged; it has no
path by which a proved implication could raise an E1 to anything else. PIR facts
carry `numeric_warrant` and `logical_implication_warrant` as separate fields, so
the degree-3 fact reads `E1 / FORMAL` and the formal certificate itself reads
`None / FORMAL`.

## 2. What is proved (WO-RH-39, 40, 41, 42)

Ten theorems, in `math/rh_weil/formal/`, with no `sorry` and no axioms beyond
`propext`, `Classical.choice` and `Quot.sound`.

**Sylvester's law of inertia**, in the formulation the runtime uses. The positive
index is defined as the largest dimension of a subspace on which the form is
positive — not as a count of positive eigenvalues, because no rigorous path in
this program may reach an eigenvalue solver. Congruence by an invertible `S`
preserves that index, the corresponding negative index, and the rank. Together
those three are what make the signature the LDL elimination reports a property of
the matrix rather than of the elimination order.

The proof rests on one identity, `qform_congruence`:
`xᵀ(SᵀAS)x = (Sx)ᵀA(Sx)`. Congruence is a change of variable; a witness subspace
for `A` maps to a witness subspace of the same dimension for `SᵀAS` under
`S⁻¹`, and back again.

**The 2×2 and 3×3 leading-principal-minor criteria**, both directions. The 3×3
case is an exact, division-free polynomial identity:

`a·D₂·Q = D₂(au + bv + cw)² + (D₂v + (ae − bc)w)² + a·D₃·w²`

with `D₂ = ad − b²` and `D₃ = det`. The three coefficients on the right are
exactly the three leading principal minors, and everything else in the proof is a
case split on which of `u, v, w` is nonzero. The converse comes from Mathlib's
`PosDef.det_pos` and `PosDef.submatrix`.

**The certificate implication.** Given honest enclosures for the three distinct
entries of a real symmetric 2×2, a positive certified lower bound on the leading
entry, and a positive certified lower bound on the determinant taken at the worst
corner of those enclosures, the matrix is positive definite. The only real
content is that a square is bounded by the larger endpoint square even when the
interval straddles zero, which is why the certificate's determinant bound uses
`max(g0b.lo², g0b.hi²)` rather than an average.

**The Weil parity identities.** Midpoint reflection; parity of `1, b, q1, b3`;
cross-parity blocks of a reflection-invariant pairing vanish; the parity
determinant factorization. The cross-block theorem is stated over an abstract
pairing with two hypotheses — linearity in the second argument and invariance
under simultaneous reflection — so nothing about digamma, primes or integration
enters. That is why it is provable here while the numerical Gram entries are not.

**Rank–trace, `Q = 0`.** `2·tr P − ‖P‖²_HS ≤ rank P` for a spectrum in `[0,1]`,
which is the case the ENG-006 degree-3 certificate uses.

## 3. What is deliberately not proved (WO-RH-40)

The general rank–trace inequality — `Q ≠ 0`, positive-index bound `b` — is
carried as

```lean
def RankTraceGeneralStatement : Prop := ...
```

with no inhabitant anywhere in the project, and recorded in the manifest as
`EXTERNAL_THEOREM_PENDING_FORMAL_PROOF` with a null warrant.

This is the honest shape for it. Recording the statement gives a future proof a
fixed target and stops the theorem being quietly assumed; making it a `def`
rather than a `theorem` means no term inhabits it, so it cannot be mistaken for
proved. The E1 runtime rank–trace result keeps its own warrant either way — Lean
does not retroactively strengthen a numerical certificate, and the absence of a
formal proof does not weaken one.

`external/zeta23/MAPPING.md` records that this exact theorem *is* proved
upstream, as `Zeta23.RHLinalg.rank_trace_ineq_two`, Apache-2.0, at a pinned
commit. Atlas does not import it: the toolchains do not currently compose, and
§11 permits an external formal dependency only after that compatibility is
proven.

## 4. The comparator and the axiom audit (WO-RH-43)

Four layers, following the architecture `anthropics/zeta-23-lean` demonstrates:

* **A — Definitions.** Canonical objects, no assertions.
* **B — Statements.** The propositions Atlas claims to rely on, written over
  Layer A and Mathlib only, each spelled out rather than abbreviated.
* **C — Solution.** The proofs, each declared *at the type of* its trusted
  statement, so the Lean elaborator performs the comparison.
* **D — Comparator.** `PrintAxioms.lean` re-derives the comparison with `isDefEq`
  and emits the axiom set and the elaborated statement, which
  `scripts/check_formal_manifest.py` hashes against the recorded manifest.

The gate has two layers because they need different things to run. The offline
layer needs no Lean: source hashes of every file that defines what Atlas claims,
no `sorry`, no project-local `axiom`, the manifest's own content hash, and the
pinned toolchain matching the lakefile. The Lean layer adds `lake build`, the
comparator run, statement-hash equality and the axiom containment.

Both were exercised against deliberate regressions before being trusted.
Weakening `0 < det` to `0 ≤ det` in a statement fails the offline hash check and
fails `Solution.lean` elaboration under Lean. An appended `sorry` is caught by
the source scan — which strips Lean comments first, because the modules' own
prose about `sorry` and `axiom` would otherwise trip it, and a gate that flags
its own documentation is a gate people learn to ignore.

Atlas's version differs from upstream's in two respects, stated rather than
glossed. It keeps the statement text in the library and hashes the source region
instead of freezing an inlined copy, which avoids two independently maintained
copies of a proposition. And it does not run the `leanprover/comparator` tool:
upstream's Layer D is a sandboxed replay with an independent kernel, Atlas's is
the elaborator plus a hash. Atlas's is weaker — it does not defend against a
compromised toolchain.

## 5. The documentation truth pass (WO-RH-37, WO-RH-45)

Every prior work order in this program was executed by an agent that read the
repository's documentation first. A stale instruction there is not a tidiness
problem: it sends the next run to regenerate a closed result, or to reinterpret a
certificate under a normalization that was rejected. So documentation is a merge
gate.

What was actually stale, and is now fixed:

<!-- docs-check: superseded-quote start -->
* `math/rh_weil/README.md` claimed WO-RH-05's interval E1 was open and listed
  only WO-RH-01 through 07 as executed — six work orders and three engineering
  specifications behind.
* `AGENT_INSTRUCTIONS.md` was the original integration work order, still
  instructing agents not to begin degree 3. It is now the live entry point; the
  original is preserved at `docs/history/agent-instructions-initial-integration.md`.
* `notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md` described the direct-Fourier
  uniform check as a future target, when the certificate exists and is promoted.
* `INTEGRATION_MANIFEST.md` still carried the ENG-002 work-order table with
  WO-RH-08 marked blocked.
* The root `README.md` linked v0.5 of the Constant Atlas in its opening block
  while stating later that v0.6 is current.
* `certificates/work_order_status.json` was nineteen orders stale, because only
  the fast path writes it and the rigorous chain never did.
<!-- docs-check: superseded-quote end -->

Two checksum manifests, `SHA256SUMS.txt` and `docs/SHA256SUMS.txt`, were also
wholesale stale — twelve of twelve and three of four entries failing
verification. They date from the delivery-package era and list build outputs
that have been regenerated many times since, including under a normalization
that was later rejected. A checksum manifest that fails `sha256sum -c` is worse
than none, because a reader concludes the tree is corrupt. Both now carry a
header saying what they are, and point at the mechanism that actually enforces
integrity here: each certificate's own `content_hash`, its
`dependencies.source_hashes`, and the promotion predicate that refuses it the
moment one drifts.

`scripts/check_docs.py` now enforces the invariant rather than the fix: local
link integrity, no live document repeating a status the repository has closed,
the claim boundary present in every live document, the README's status table
agreeing with `work_order_status.json`, and the README's certified numbers read
from the certificates at check time. The current Constant Atlas version is
detected from repository state rather than hard-coded, since pinning a version in
the checker only moves the staleness from the README into the gate.

Superseded documents are preserved and labelled, never deleted — the same rule
WO-RH-17 applied to contrary numerical evidence.

## 6. The 3×3 pilot (WO-RH-46)

ENG-006's "degree-3" odd block is 2×2, and its own information report says what
that cost: on a 2×2 the trace and determinant fix the spectrum, so inertia,
moments and rank–trace all collapsed onto information the determinant already
carried.

Writing `u = x − L/2`, the odd sector through degree 3 is exactly `span{u, u³}`,
and `q1³ = (L²/4)q1 − b3` lies in it — extending the odd sector at degree 3 is
impossible, which is why ENG-006 stopped there. The even sector extends:
`{1, b, b²}` spans `{1, u², u⁴}`. It is also the sector worth extending, since it
contains the constant function and its leading 2×2 determinant is exactly the
`E2` that ENG-004 and ENG-005 certified.

ENG-007 **prepares** it and does not certify it: exact kernel identities (E0) and
a conditioning and topology preview (E3), with no new E1 degree result, per §15.

The preview found inertia `(3,0,0)` across a nine-point grid, all three leading
minors positive, and a condition number up to `1.3e5` dropping to at most `2.4e1`
under the Jacobi rescaling `D M D`. That rescaling is a congruence, so it changes
no part of the inertia — which is `posIndexAtLeast_congruence_iff` and
`rank_congruence`, both proved above. The conditioning fix is free, and the
theorem is what says so.

Two things recorded now rather than rediscovered later. `b2` is quadratic in `L`,
unlike every basis element ENG-005 works with, and `pole.py`'s `_laplace_d2L`
drops the second integral precisely because `d²_L h = 0` for a linear-in-`L`
element — so an `E2''`-style curvature argument on this block needs that
machinery extended first. And the preview's Fourier transform needed a Taylor
branch for small `|t|`: the integration-by-parts recursion divides by `s = it` at
every step, and mpmath's tanh-sinh quadrature clusters nodes exponentially close
to panel endpoints, so `t` reaches `1e-25` and the naive recursion returned Gram
entries around `1e30`.

The pilot lives in its own module with its own basis table rather than extending
`pole.py` or `core.py`. Extending those would change their source hashes, which
every promoted E1 certificate binds — a preparatory pilot must not be able to
invalidate a certified result. A test asserts the separation holds.

## 7. Validation

* Full rigorous chain re-run end to end: OK. Nine rigorous certificates with
  current source hashes; seventeen PIR facts.
* The pilot's independent mpmath assembly reproduces the certified Arb 2×2 at
  `L = log 3, T = 84` to every digit displayed:
  `G00 = 0.107356700415`, `G0b = 0.000461820208772`,
  `Gbb = 3.42537786464e-5`, `E2 = 3.46409474697e-6`.
* `lake build` clean, no warnings, no `sorry`; ten theorems audited by the
  comparator; axiom set exactly the three standard ones.
* `scripts/ci_formal.py` and `scripts/check_docs.py` both green.
