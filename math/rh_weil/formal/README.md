---
status: CURRENT
work_order: ATLAS-RH-ENG-007
---

# AtlasRH — the formal theorem boundary

Machine-checked statements of the finite theorems the RH/Weil runtime relies on.

**Claim boundary: finite-dimensional linear algebra and certificate semantics only. No RH
proof claim, and nothing here can produce one.** `rh_proof_claim` is `false`.

## What problem this solves

ENG-007 §0: the remaining risk is no longer primarily numerical, it is **statement drift** —
the code proves one finite theorem while documentation, certificate consumers, or a future
agent believe it proved another.

Concretely: `e1_degree2_compact_log3_log4.json` reports a certified lower bound on
`E2 = G00·Gbb − G0b²`. A consumer reads that and concludes "the Gram block is positive
definite on the cell". That conclusion is a *theorem*, and until now it lived in people's
heads. It is now proved, and the proof is checked on every build.

## The four layers

| Layer | File | Role |
|---|---|---|
| A Definitions | `AtlasRH/Definitions.lean` | canonical finite-dimensional vocabulary |
| B TrustedStatements | `comparator/Comparator/TrustedStatements.lean` | the propositions Atlas claims to rely on |
| C Solution | `comparator/Comparator/Solution.lean` | the proofs, delegating to `AtlasRH` |
| D Comparator / audit | `Solution` type ascriptions + `PrintAxioms.lean` | statement comparison and axiom audit |

**TrustedStatements does not import AtlasRH.** That is deliberate. A trusted statement
phrased in the implementation's own vocabulary changes silently whenever that vocabulary
changes, which is exactly the drift being defended against.

**The comparator is the Lean kernel.** Each theorem in `Solution` has its type ascribed to a
`TrustedStatements` definition; if a proof ever drifts, the file stops typechecking. No
external diffing tool is involved, and none could be more trustworthy. Verified negatively:
a true-but-weaker one-directional theorem offered in place of the trusted `iff` is a type
error, not a silent pass.

## What is proved

| Theorem | Statement |
|---|---|
| `congruence_preserves_posDef` | `Sᴴ A S ≻ 0 ↔ A ≻ 0` for invertible `S` — licenses pivoting into a congruent block and reading a sign there |
| `pd_two_by_two` | `!![a,b;b,c] ≻ 0 ↔ a > 0 ∧ ac − b² > 0`, both directions |
| `schur_pivot_implies_posDef` | a positive Schur pivot transcript implies positive definiteness |
| `certificate_even2_implies_pd` | rigorous lower bounds on `G₀₀` and `det`, pointwise over a domain, imply `G(L) ≻ 0` on that domain |

## The trust boundary, stated plainly

- **Python/Arb** is responsible for producing true numerical enclosures. Nothing in Lean
  checks that, and nothing here could.
- **Lean** is responsible for proving that those enclosures, if true, imply the advertised
  finite mathematical conclusion.

A formal theorem therefore strengthens an *exact theorem dependency*. It does **not** convert
interval numerical evidence into FORMAL evidence. An E1 certificate stays E1.

This is not a verified interval integrator and does not claim to be.

## Build and audit

```bash
cd math/rh_weil/formal
lake build                                        # AtlasRH + Comparator
lake env lean comparator/Comparator/PrintAxioms.lean   # statements + axiom audit
python3 ../../../scripts/check_formal_manifest.py      # manifest + allowlist gate
```

Pinned: Lean `v4.33.0`, Mathlib `v4.33.0` (`db584cd`), recorded in `lake-manifest.json`.
No floating resolution — a formal layer that cannot be replayed is not evidence.

Allowed axioms are `propext`, `Classical.choice`, `Quot.sound`, enumerated in
`manifests/theorem_manifest.json`. `sorryAx` is deliberately outside the allowlist: an
incomplete proof must not be able to wear a FORMAL label.

## Not yet formal

The rank–trace bound is **not** proved here. It remains an E1 runtime result under its
existing warrant, and ENG-007 §7 is explicit that Lean does not retroactively strengthen it.
See `AtlasRH/RankTrace.lean` for the exact statement and its isolated non-formal status.
