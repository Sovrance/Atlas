#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §16 — the RH CI gates, as runnable names.

    python3 scripts/ci_rh.py --gate rh-docs      # documentation truth (fast, stdlib only)
    python3 scripts/ci_rh.py --gate rh-formal    # Lean build + comparator + axiom audit
    python3 scripts/ci_rh.py --gate all

Why a script and not YAML: `math/rh_weil/scripts/ci_inertia.py` already establishes the
convention -- "the gates are exposed the same way rather than as YAML that nothing would
run" -- and this repository has no workflow runner. A gate that exists only as
configuration for a system nobody invokes is not a gate.

§16 is explicit that a docs-only failure blocks merge, and gives the reason: stale
instructions can cause a future agent to regenerate or reinterpret the wrong mathematics.
`rh-docs` is therefore a hard failure, not a warning.

`rh-fast`, `rh-rigorous` and the inertia gates are unchanged and still live in
`math/rh_weil/scripts/`.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FORMAL = REPO / "math" / "rh_weil" / "formal"


def _run(cmd: list[str], cwd: Path, env: dict | None = None) -> int:
    print(f"    $ {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, cwd=str(cwd), env=env).returncode


def _lean_env() -> dict | None:
    """Lean lives outside the default PATH in some environments.

    Returns None when no toolchain can be found, so the caller can report the gate as
    unrunnable rather than silently passing it.
    """
    env = dict(os.environ)
    for home in (Path(env.get("ELAN_HOME", "")), Path.home() / ".elan", Path("/home/user/.elan")):
        if home and (home / "bin" / "lake").exists():
            env["ELAN_HOME"] = str(home)
            env["PATH"] = f"{home / 'bin'}:{env.get('PATH', '')}"
            return env
    from shutil import which
    return env if which("lake") else None


def gate_docs() -> int:
    print("== rh-docs ==")
    return _run([sys.executable, "scripts/check_docs.py"], REPO)


def gate_formal() -> int:
    print("== rh-formal ==")
    env = _lean_env()
    if env is None:
        print("    ERROR: no Lean toolchain found (elan/lake). The formal gate cannot be\n"
              "    satisfied by skipping it: an unrun gate is not a passing gate.",
              file=sys.stderr)
        return 1
    rc = _run(["lake", "build"], FORMAL, env)
    if rc:
        return rc
    return _run([sys.executable, "scripts/check_formal_manifest.py"], REPO, env)


GATES = {"rh-docs": gate_docs, "rh-formal": gate_formal}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate", choices=[*GATES, "all"], default="all")
    args = ap.parse_args()
    names = list(GATES) if args.gate == "all" else [args.gate]

    failures = []
    for name in names:
        if GATES[name]():
            failures.append(name)
        print()

    if failures:
        print(f"RH CI: FAIL ({', '.join(failures)})", file=sys.stderr)
        return 1
    print(f"RH CI: OK ({', '.join(names)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
