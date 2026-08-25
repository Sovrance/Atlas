# `zeta-23-lean` — provenance


**No RH proof claim is made** by Atlas here or anywhere. Upstream proves
zero-proportion theorems about the zeta function; Atlas takes none of them, and
nothing in this directory carries a warrant.

**Status in Atlas: reference only. Nothing from this project is imported, vendored,
or depended on by any Atlas build, and no Atlas certificate derives any warrant
from it.** This directory records the pin and the mapping so that if that ever
changes, it changes deliberately and with provenance attached.

## The upstream project

| | |
|---|---|
| repository | <https://github.com/anthropics/zeta-23-lean> |
| pinned commit | `cec57f919ccf34e5fa5372b4ba332f7c848bbb6e` (branch `main`) |
| pin resolved | 2026-08-25 (`git ls-remote`, then `git clone --depth 1`) |
| library | `Zeta23` |
| license | Apache License 2.0 (`LICENSE` at the pinned commit) |
| copyright | Copyright 2026 Anthropic, PBC (`NOTICE`) |
| Lean toolchain | `leanprover/lean4:v4.33.0-rc2` |
| Mathlib commit | `51e6992efd06126df61a496bebf8f49482a4e129` |
| subject | a Lean 4 formalization of critical-line zero-proportion theorems (Alpöge–Furman, arXiv:2608.13637) |

Upstream itself derives files under `Zeta23/FromPNTPlus/` from
[PrimeNumberTheoremAnd](https://github.com/AlexKontorovich/PrimeNumberTheoremAnd)
(Apache-2.0), some of which in turn derive from Mathlib (Apache-2.0). That chain
is recorded in upstream's `NOTICE` and is reproduced here so a future importer
does not have to rediscover it.

## Toolchain compatibility

Atlas's formal project pins `leanprover/lean4:v4.34.0-rc2` and Mathlib
`f1c1e67f08f57b6d7088b1a98fdceab6da4407ee`. Upstream pins one Lean release
earlier and a different Mathlib commit. **The two do not currently compose.**
§11 of ENG-007 permits an external formal dependency "only after toolchain
compatibility is proven"; it has not been proven, and until it is, this
directory is documentation.

## What Atlas took, and what it did not

Atlas took the **architecture**: a trusted statement layer written over Mathlib
alone, a solution layer that may import the implementation library, a comparator
that machine-checks the two describe the same proposition, and a `#print axioms`
report enumerating the standard axioms rather than hiding them. That shape is
implemented independently in `math/rh_weil/formal/` — see `AtlasRH/Statements.lean`,
`comparator/Solution.lean` and `comparator/PrintAxioms.lean`. Atlas's version
differs where its needs differ: it keeps the statement text in the library and
hashes the source region, rather than inlining a frozen copy, and its manifest
gate runs offline as well as under Lean.

Atlas took **no Lean source, no theorem, and no proof.** No file in this
repository is copied from upstream. Everything in `math/rh_weil/formal/` was
written for Atlas and is proved against Atlas's own pinned Mathlib.

## If Atlas ever imports an upstream theorem

§11 requires all of the following to be recorded before an import is relied on,
and `theorem_manifest.json` in this directory carries the fields for them:

* upstream commit (pinned above, and per-theorem in the manifest);
* theorem name, exactly as declared upstream;
* statement hash;
* license and provenance (Apache-2.0, chain above);
* axiom report for the imported declaration.

Note what the manifest's per-theorem hashes currently are and are not. They are
**source-text hashes** of the declaration as it reads at the pinned commit,
whitespace-normalized. They are not elaborated-statement hashes, because
elaborating them means building `Zeta23` under its own toolchain, which has not
been done here. Promoting a mapping from reference to dependency requires that
build, and the elaborated hash it produces — the source hash is a pin, not a
statement comparison.

**External theorem reuse must never become provenance-free E0.** An imported
theorem is only ever as good as this record; if the record is missing, the
import does not happen.
