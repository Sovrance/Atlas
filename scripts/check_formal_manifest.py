#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §10 (WO-RH-43) — formal theorem manifest and axiom audit.

    python3 scripts/check_formal_manifest.py            # verify committed manifest
    python3 scripts/check_formal_manifest.py --write    # regenerate it

What this gate is for
---------------------
The Lean kernel already guarantees that every theorem in ``Comparator.Solution`` *is* the
proposition declared in ``Comparator.TrustedStatements`` -- that comparison is the type
ascription, and it cannot be fooled. What the kernel does not do is tell a *certificate
consumer* that the statement it relied on last month is the statement being proved today.

So the manifest hashes the elaborated statements. A drifted statement still builds -- it is
simply a different theorem -- but its hash moves, and this gate fails. That is the ENG-007
§5 failure mode made mechanical: "prevent a future proof of a subtly changed theorem from
satisfying the original certificate consumer".

The axiom audit is the second half. ``sorryAx`` in the list means a proof is incomplete;
anything outside the allowlist means a project-specific axiom was introduced. Both are
merge-blocking, because either would let an unproved claim wear a FORMAL label.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FORMAL = REPO / "math" / "rh_weil" / "formal"
MANIFEST = FORMAL / "manifests" / "theorem_manifest.json"
AUDIT_SRC = Path("comparator") / "Comparator" / "PrintAxioms.lean"

#: Standard Lean/Mathlib axioms. Enumerated rather than hidden (§10). `sorryAx` is
#: deliberately ABSENT: a `sorry` reaching the promoted library must fail this gate.
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}

#: The theorems Atlas exports. Adding one here without adding it to PrintAxioms.lean fails.
EXPORTED = [
    "congruence_preserves_posDef",
    "pd_two_by_two",
    "schur_pivot_implies_posDef",
    "certificate_even2_implies_pd",
    "even_odd_cross_vanishes",
    "det_parity_factorization",
    "det_congruence_invariant",
]


def _lean_env() -> dict:
    env = dict(os.environ)
    elan = Path("/home/user/.elan")
    if elan.is_dir():
        env["ELAN_HOME"] = str(elan)
        env["PATH"] = f"{elan / 'bin'}:{env.get('PATH', '')}"
    return env


def run_audit() -> str:
    proc = subprocess.run(
        ["lake", "env", "lean", str(AUDIT_SRC)],
        cwd=str(FORMAL), env=_lean_env(), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(
            "formal audit did not run (is the project built? `lake build`):\n"
            + proc.stderr[-2000:]
        )
    return proc.stdout


def parse_audit(text: str) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Split the audit output into statements and axiom lists.

    A `#print` block runs until the next `def ` or `'...' depends on axioms:` line, so the
    statements are reassembled by accumulating continuation lines rather than assuming one
    statement per line -- several of them wrap.
    """
    statements: dict[str, str] = {}
    axioms: dict[str, list[str]] = {}

    current: str | None = None
    buf: list[str] = []

    def flush() -> None:
        if current is not None:
            # Normalise whitespace: pretty-printer line breaks are layout, not content.
            statements[current] = " ".join(" ".join(buf).split())

    for line in text.splitlines():
        m_def = re.match(r"def Comparator\.TrustedStatements\.(\w+) : Prop :=(.*)$", line)
        m_ax = re.match(r"'Comparator\.Solution\.(\w+)' depends on axioms: \[(.*)\]\s*$", line)
        if m_def:
            flush()
            current = m_def.group(1)
            buf = [m_def.group(2)]
            continue
        if m_ax:
            flush()
            current = None
            buf = []
            axioms[m_ax.group(1)] = [a.strip() for a in m_ax.group(2).split(",") if a.strip()]
            continue
        if current is not None:
            buf.append(line)
    flush()
    return statements, axioms


def statement_hash(statement: str) -> str:
    return hashlib.sha256(statement.encode("utf-8")).hexdigest()


def toolchain() -> str:
    return (FORMAL / "lean-toolchain").read_text(encoding="utf-8").strip()


def mathlib_rev() -> str:
    manifest = json.loads((FORMAL / "lake-manifest.json").read_text(encoding="utf-8"))
    for pkg in manifest.get("packages", []):
        if pkg.get("name") == "mathlib":
            return pkg.get("rev", "")
    return ""


def build_manifest() -> dict:
    statements, axioms = parse_audit(run_audit())
    missing = [t for t in EXPORTED if t not in statements or t not in axioms]
    if missing:
        raise SystemExit(
            "these theorems are declared exported but were not emitted by "
            f"{AUDIT_SRC}: {missing}"
        )
    return {
        "formal_project": "AtlasRH",
        "work_order": "ATLAS-RH-ENG-007",
        "lean_toolchain": toolchain(),
        "mathlib_rev": mathlib_rev(),
        "mathlib_inputRev": "v4.33.0",
        "rh_proof_claim": False,
        "claim_scope": "finite_dimensional_linear_algebra_and_certificate_semantics",
        "allowed_axioms": sorted(ALLOWED_AXIOMS),
        "theorems": [
            {
                "id": name,
                "trusted_statement": f"Comparator.TrustedStatements.{name}",
                "solution_theorem": f"Comparator.Solution.{name}",
                "statement": statements[name],
                "statement_hash": statement_hash(statements[name]),
                "axioms": sorted(axioms[name]),
            }
            for name in EXPORTED
        ],
    }


def check() -> int:
    fresh = build_manifest()
    failures: list[str] = []

    for entry in fresh["theorems"]:
        extra = sorted(set(entry["axioms"]) - ALLOWED_AXIOMS)
        if extra:
            kind = "INCOMPLETE PROOF (sorry)" if "sorryAx" in extra else "non-standard axiom"
            failures.append(f"{entry['id']}: {kind}: {extra}")

    if not MANIFEST.exists():
        failures.append(f"missing manifest {MANIFEST.relative_to(REPO)}; run with --write")
    else:
        committed = json.loads(MANIFEST.read_text(encoding="utf-8"))
        by_id = {t["id"]: t for t in committed.get("theorems", [])}
        for entry in fresh["theorems"]:
            old = by_id.get(entry["id"])
            if old is None:
                failures.append(f"{entry['id']}: not in committed manifest")
            elif old.get("statement_hash") != entry["statement_hash"]:
                failures.append(
                    f"{entry['id']}: STATEMENT DRIFT\n"
                    f"    committed: {old.get('statement')}\n"
                    f"    current  : {entry['statement']}"
                )
        for stale in sorted(set(by_id) - {t["id"] for t in fresh["theorems"]}):
            failures.append(f"{stale}: in committed manifest but no longer exported")
        for key in ("lean_toolchain", "mathlib_rev"):
            if committed.get(key) != fresh[key]:
                failures.append(
                    f"{key}: pinned {committed.get(key)!r}, built {fresh[key]!r}"
                )

    if failures:
        print("FORMAL MANIFEST GATE: FAIL", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    print(f"FORMAL MANIFEST GATE: OK ({len(fresh['theorems'])} theorems, "
          f"axioms within allowlist, statements unchanged)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="regenerate the manifest")
    args = ap.parse_args()
    if args.write:
        manifest = build_manifest()
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
        print(f"wrote {MANIFEST.relative_to(REPO)} ({len(manifest['theorems'])} theorems)")
        return 0
    return check()


if __name__ == "__main__":
    raise SystemExit(main())
