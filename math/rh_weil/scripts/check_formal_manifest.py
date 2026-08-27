#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §10 (WO-RH-43) -- formal theorem manifest gate.

Two layers, because they need different things to run.

**Offline layer** (always runs, no Lean required)

  * every Lean source the manifest names still hashes to the recorded value;
  * no ``sorry`` anywhere under ``formal/AtlasRH/`` or in
    ``formal/comparator/Solution.lean``;
  * no project-specific ``axiom`` declaration anywhere in the formal project;
  * the manifest's own ``manifest_id`` is the content hash of its body;
  * the pinned toolchain and Mathlib commit in the manifest match
    ``lean-toolchain`` and ``lakefile.toml``.

  This layer alone is the drift gate. A statement cannot change without
  changing ``AtlasRH/Statements.lean``, and that file's hash is recorded.

**Lean layer** (``--with-lean``, or automatically when ``lake`` is on PATH)

  * ``lake build`` succeeds;
  * ``lake env lean comparator/PrintAxioms.lean`` runs, which itself performs
    the Layer D comparison (``isDefEq`` between each solution theorem's type
    and the trusted statement's body) and fails loudly if they ever come apart;
  * the elaborated, pretty-printed statement of every theorem hashes to the
    recorded ``statement_hash``;
  * every theorem's axiom set is a subset of the enumerated allowed axioms.

``--write`` regenerates the manifest from a live Lean run. It refuses to run
without Lean, since a manifest written from stale numbers would be worse than
no manifest.

No RH proof claim is made or checked here. This gate is about whether the
finite theorems Atlas cites are the finite theorems Atlas proved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
FORMAL = ROOT / "formal"
MANIFEST = FORMAL / "manifests" / "theorem_manifest.json"

# The three axioms every Mathlib development already depends on. Enumerated
# rather than hidden, per §10: an axiom outside this set is a finding, not a
# detail.
ALLOWED_AXIOMS = ("Classical.choice", "Quot.sound", "propext")

# Files whose content defines what Atlas claims. Any edit to one of these
# changes the manifest, and therefore has to be a deliberate act.
HASHED_SOURCES = (
    "formal/lean-toolchain",
    "formal/lakefile.toml",
    "formal/AtlasRH.lean",
    "formal/AtlasRH/Definitions.lean",
    "formal/AtlasRH/Statements.lean",
    "formal/AtlasRH/Positivity.lean",
    "formal/AtlasRH/MatrixInertia.lean",
    "formal/AtlasRH/WeilBasis.lean",
    "formal/AtlasRH/RankTrace.lean",
    "formal/AtlasRH/CertificateSemantics.lean",
    "formal/AtlasRH/GeneralizedGap.lean",
    "formal/AtlasRH/Sylvester4.lean",
    "formal/AtlasRH/NestedSchur.lean",
    "formal/comparator/TrustedStatements.lean",
    "formal/comparator/Solution.lean",
    "formal/comparator/PrintAxioms.lean",
)

# §5 Layer B: the statement file must not import a module that contains a
# proof, or a statement could be weakened by editing the vocabulary it is
# written in.
STATEMENTS_ALLOWED_IMPORTS = ("AtlasRH.Definitions",)

# Statements Atlas records but has NOT proved. These are named here so the
# absence of a proof is a fact in the manifest rather than an omission, and so
# a future proof has a fixed target. None of them carries any warrant.
UNPROVED = (
    {
        "id": "rank_trace_general",
        "lean_name": "AtlasRH.RankTraceGeneralStatement",
        "status": "EXTERNAL_THEOREM_PENDING_FORMAL_PROOF",
        "warrant": None,
        "note": (
            "the rank-trace inequality with Q != 0 and a positive-index bound b. "
            "Carried as a `def ... : Prop` with no inhabitant anywhere in the "
            "project, so it cannot be mistaken for proved and cannot upgrade the "
            "warrant of the E1 runtime result that uses the Q = 0 case."
        ),
    },
)


# --------------------------------------------------------------------------- #
# helpers                                                                      #
# --------------------------------------------------------------------------- #
def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def file_sha256(rel: str) -> Optional[str]:
    p = ROOT / rel
    if not p.exists():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


def canonical_id(body: Dict[str, Any]) -> str:
    """Content hash of the manifest body, excluding its own id and timestamp."""
    stripped = {k: v for k, v in body.items() if k not in ("manifest_id", "generated_utc")}
    blob = json.dumps(stripped, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return "formal_sha256_" + sha256_text(blob)


def normalize_statement(s: str) -> str:
    """Collapse whitespace. Lean's pretty printer is deterministic under a
    pinned toolchain, but line wrapping depends on the terminal width it was
    invoked with, and that is not part of the statement's meaning."""
    return re.sub(r"\s+", " ", s).strip()


def read_toolchain() -> str:
    return (FORMAL / "lean-toolchain").read_text(encoding="utf-8").strip()


def read_mathlib_commit() -> Optional[str]:
    text = (FORMAL / "lakefile.toml").read_text(encoding="utf-8")
    block = re.search(
        r'name\s*=\s*"mathlib".*?rev\s*=\s*"([0-9a-f]{40})"', text, re.S
    )
    return block.group(1) if block else None


def lean_env() -> Dict[str, str]:
    env = dict(os.environ)
    extra = os.environ.get("ATLAS_LEAN_BIN")
    if extra:
        env["PATH"] = f"{extra}{os.pathsep}" + env.get("PATH", "")
    return env


def lake_available() -> bool:
    env = lean_env()
    return shutil.which("lake", path=env.get("PATH")) is not None


# --------------------------------------------------------------------------- #
# the Lean run                                                                 #
# --------------------------------------------------------------------------- #
def run_lean() -> Tuple[List[Dict[str, Any]], List[str]]:
    """Build the project and collect the comparator report.

    Returns ``(theorems, notes)``. Raises ``RuntimeError`` on any Lean failure,
    including a failed Layer D comparison -- ``PrintAxioms.lean`` throws in that
    case, so a mismatch surfaces as a nonzero exit rather than as a missing line.
    """
    env = lean_env()
    notes: List[str] = []
    build = subprocess.run(
        ["lake", "build"], cwd=FORMAL, env=env, capture_output=True, text=True
    )
    if build.returncode != 0:
        raise RuntimeError(f"lake build failed:\n{build.stdout}\n{build.stderr}")
    combined = build.stdout + build.stderr
    if "declaration uses 'sorry'" in combined:
        raise RuntimeError("lake build reported a declaration using 'sorry'")

    proc = subprocess.run(
        ["lake", "env", "lean", "comparator/PrintAxioms.lean"],
        cwd=FORMAL,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"comparator/PrintAxioms.lean failed:\n{proc.stdout}\n{proc.stderr}"
        )

    theorems: List[Dict[str, Any]] = []
    count: Optional[int] = None
    for line in proc.stdout.splitlines():
        if line.startswith("ATLAS_FORMAL_THEOREM_COUNT\t"):
            count = int(line.split("\t")[1])
            continue
        if not line.startswith("ATLAS_FORMAL_THEOREM\t"):
            continue
        _, tid, solution, trusted, axioms, statement = line.split("\t", 5)
        axs = [] if axioms == "none" else sorted(axioms.split(","))
        stmt = normalize_statement(statement)
        theorems.append(
            {
                "id": tid,
                "trusted_statement": trusted,
                "solution_theorem": solution,
                "statement_hash": "sha256:" + sha256_text(stmt),
                "axioms": axs,
            }
        )
    if count is None:
        raise RuntimeError("comparator report carried no theorem count")
    if count != len(theorems):
        raise RuntimeError(
            f"comparator reported {count} theorems but emitted {len(theorems)} lines"
        )
    notes.append(f"lean: {len(theorems)} theorems audited")
    return theorems, notes


# --------------------------------------------------------------------------- #
# offline checks                                                               #
# --------------------------------------------------------------------------- #
def strip_lean_comments(text: str) -> str:
    """Blank out Lean comments while preserving line numbering.

    Needed because this file's own subject matter -- the words ``sorry`` and
    ``axiom`` -- appears in the prose of the modules being scanned. A scanner
    that flagged its own documentation would train everyone to ignore it.
    Nested ``/- -/`` blocks are handled, since Lean allows them.
    """
    out: List[str] = []
    i, n, depth = 0, len(text), 0
    while i < n:
        ch = text[i]
        if depth == 0 and text.startswith("/-", i):
            depth = 1
            out.append("  ")
            i += 2
            continue
        if depth > 0:
            if text.startswith("/-", i):
                depth += 1
                out.append("  ")
                i += 2
                continue
            if text.startswith("-/", i):
                depth -= 1
                out.append("  ")
                i += 2
                continue
            out.append("\n" if ch == "\n" else " ")
            i += 1
            continue
        if text.startswith("--", i):
            j = text.find("\n", i)
            if j < 0:
                out.append(" " * (n - i))
                break
            out.append(" " * (j - i))
            i = j
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def check_no_sorry() -> List[str]:
    """A `sorry` in the proof library would silently void every downstream
    claim, so this looks at the source rather than trusting the build log."""
    problems: List[str] = []
    targets = sorted((FORMAL / "AtlasRH").glob("*.lean"))
    targets.append(FORMAL / "comparator" / "Solution.lean")
    pattern = re.compile(r"(?<![A-Za-z_.])sorry(?![A-Za-z_])")
    for p in targets:
        if not p.exists():
            continue
        code = strip_lean_comments(p.read_text(encoding="utf-8"))
        for n, line in enumerate(code.splitlines(), 1):
            if pattern.search(line):
                problems.append(f"{p.relative_to(FORMAL)}:{n}: sorry")
    return problems


def check_no_local_axioms() -> List[str]:
    """`axiom` declarations of our own would put an unaudited assumption under
    the whole library. Mathlib's three are allowed; ours are not."""
    problems: List[str] = []
    for p in sorted(FORMAL.rglob("*.lean")):
        if ".lake" in p.parts:
            continue
        code = strip_lean_comments(p.read_text(encoding="utf-8"))
        for n, raw in enumerate(code.splitlines(), 1):
            if re.match(r"^\s*(private\s+|protected\s+)?axiom\s", raw):
                problems.append(f"{p.relative_to(FORMAL)}:{n}: {raw.strip()}")
    return problems


def check_statements_imports() -> List[str]:
    p = FORMAL / "AtlasRH" / "Statements.lean"
    problems: List[str] = []
    code = strip_lean_comments(p.read_text(encoding="utf-8"))
    for n, raw in enumerate(code.splitlines(), 1):
        m = re.match(r"^import\s+(\S+)", raw)
        if not m:
            continue
        mod = m.group(1)
        if mod.startswith("Mathlib"):
            continue
        if mod not in STATEMENTS_ALLOWED_IMPORTS:
            problems.append(
                f"AtlasRH/Statements.lean:{n}: imports {mod}, "
                "which is not a Layer A definitions module"
            )
    return problems


def check_sources(manifest: Dict[str, Any]) -> List[str]:
    problems: List[str] = []
    recorded = manifest.get("sources") or {}
    for rel in HASHED_SOURCES:
        if rel not in recorded:
            problems.append(f"{rel}: not recorded in the manifest")
    for rel, want in sorted(recorded.items()):
        got = file_sha256(rel)
        if got is None:
            problems.append(f"{rel}: missing")
        elif got != want:
            problems.append(f"{rel}: source hash changed since the manifest was written")
    return problems


# --------------------------------------------------------------------------- #
# build / verify                                                               #
# --------------------------------------------------------------------------- #
def build_manifest(theorems: List[Dict[str, Any]]) -> Dict[str, Any]:
    body: Dict[str, Any] = {
        "certificate_version": "0.1",
        "content_kind": "FORMAL_THEOREM_MANIFEST",
        "formal_project": "AtlasRH",
        "work_order": "ATLAS-RH-ENG-007",
        "program": "RH/Weil finite theorem boundary",
        "rh_proof_claim": False,
        "claim_scope": "finite_dimensional_weil_compression",
        "lean_toolchain": read_toolchain(),
        "mathlib_commit": read_mathlib_commit(),
        "allowed_axioms": list(ALLOWED_AXIOMS),
        "theorems": theorems,
        "unproved_statements": [dict(u) for u in UNPROVED],
        "sources": {rel: file_sha256(rel) for rel in HASHED_SOURCES},
        "note": (
            "Every theorem listed here is finite linear algebra over the reals. "
            "A formal theorem may strengthen an exact theorem dependency; it "
            "never converts interval numerical evidence to FORMAL."
        ),
    }
    body["manifest_id"] = canonical_id(body)
    body["generated_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return body


def verify(manifest: Dict[str, Any], with_lean: bool) -> Tuple[List[str], List[str]]:
    problems: List[str] = []
    notes: List[str] = []

    if manifest.get("rh_proof_claim") is not False:
        problems.append("manifest does not carry rh_proof_claim: false")

    want_id = canonical_id(manifest)
    if manifest.get("manifest_id") != want_id:
        problems.append(
            f"manifest_id {manifest.get('manifest_id')} does not match its body ({want_id})"
        )
    else:
        notes.append(f"manifest_id: {want_id}")

    tc = read_toolchain()
    if manifest.get("lean_toolchain") != tc:
        problems.append(
            f"manifest pins {manifest.get('lean_toolchain')!r} but lean-toolchain says {tc!r}"
        )
    mc = read_mathlib_commit()
    if mc is None:
        problems.append("lakefile.toml does not pin mathlib to a 40-hex commit")
    elif manifest.get("mathlib_commit") != mc:
        problems.append(
            f"manifest pins mathlib {manifest.get('mathlib_commit')} but lakefile says {mc}"
        )
    else:
        notes.append(f"mathlib pinned at {mc}")

    problems += check_no_sorry()
    problems += check_no_local_axioms()
    problems += check_statements_imports()
    problems += check_sources(manifest)

    allowed = set(manifest.get("allowed_axioms") or ())
    if allowed != set(ALLOWED_AXIOMS):
        problems.append(
            "manifest's allowed_axioms differ from the enumerated standard set"
        )
    for t in manifest.get("theorems") or ():
        extra = sorted(set(t.get("axioms") or ()) - allowed)
        if extra:
            problems.append(f"{t.get('id')}: depends on non-standard axioms {extra}")

    if not manifest.get("theorems"):
        problems.append("manifest lists no theorems")

    if with_lean:
        live, lean_notes = run_lean()
        notes += lean_notes
        recorded = {t["id"]: t for t in manifest.get("theorems") or ()}
        seen = {t["id"]: t for t in live}
        for tid in sorted(set(recorded) - set(seen)):
            problems.append(f"{tid}: in the manifest but not exported by the comparator")
        for tid in sorted(set(seen) - set(recorded)):
            problems.append(f"{tid}: exported by the comparator but not in the manifest")
        for tid in sorted(set(recorded) & set(seen)):
            r, s = recorded[tid], seen[tid]
            for field in ("trusted_statement", "solution_theorem", "statement_hash"):
                if r.get(field) != s.get(field):
                    problems.append(
                        f"{tid}: {field} drifted -- manifest {r.get(field)!r}, "
                        f"live {s.get(field)!r}"
                    )
            if sorted(r.get("axioms") or ()) != sorted(s.get("axioms") or ()):
                problems.append(f"{tid}: axiom set drifted")
    else:
        notes.append("lean layer skipped (no lake on PATH; source hashes still checked)")

    return problems, notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="regenerate the manifest (requires Lean)")
    ap.add_argument("--with-lean", action="store_true", help="require the Lean layer to run")
    ap.add_argument("--no-lean", action="store_true", help="offline layer only")
    args = ap.parse_args()

    have_lake = lake_available()
    if args.write:
        if not have_lake:
            print("check_formal_manifest: --write needs `lake` on PATH "
                  "(set ATLAS_LEAN_BIN to the Lean toolchain bin directory)")
            return 2
        theorems, _ = run_lean()
        body = build_manifest(theorems)
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(
            json.dumps(body, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {MANIFEST.relative_to(ROOT)}")
        print(f"  manifest_id: {body['manifest_id']}")
        print(f"  theorems: {len(body['theorems'])}")
        return 0

    if not MANIFEST.exists():
        print(f"check_formal_manifest: missing {MANIFEST.relative_to(ROOT)}")
        return 1
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))

    with_lean = have_lake and not args.no_lean
    if args.with_lean and not have_lake:
        print("check_formal_manifest: --with-lean requested but `lake` is not on PATH")
        return 2

    try:
        problems, notes = verify(manifest, with_lean)
    except RuntimeError as exc:
        print("check_formal_manifest: FAIL")
        print(f"  {exc}")
        return 1

    for n in notes:
        print(f"  {n}")
    if problems:
        print("check_formal_manifest: FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(
        f"check_formal_manifest: PASS ({len(manifest.get('theorems') or ())} theorems, "
        f"{len(manifest.get('unproved_statements') or ())} recorded unproved)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
