#!/usr/bin/env python3
"""Regenerate and verify `SHA256SUMS.txt` for the RH/Weil certificates.

    python3 math/rh_weil/scripts/write_sha256sums.py            # verify
    python3 math/rh_weil/scripts/write_sha256sums.py --write     # regenerate

Why this script exists
----------------------
`SHA256SUMS.txt` was committed but nothing generated it and nothing checked it. Audited
during ENG-007 it was stale in every single entry -- 0 of 12 matched, and it was already
stale at the ENG-006 baseline `28ec698`. A hash manifest that nobody regenerates is worse
than no manifest: it looks like an integrity control and is one only by accident.

Note what this is NOT. The `dependencies.source_hashes` recorded *inside* each certificate
are a different and working mechanism -- those pin the source files a certificate was
produced from, they are written by the certifier, and at the ENG-007 audit all 89 of them
were valid. This file pins the certificate artifacts themselves, which is a coarser and
more fragile thing: every regeneration changes `generated_utc` and therefore every hash.

Because of that, `--write` is expected after any run that regenerates certificates, and the
verify mode is deliberately NOT wired into the fast docs gate: it would fail on every
legitimate rerun and train people to ignore it. It is a release-time check.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

PROGRAM = Path(__file__).resolve().parents[1]
MANIFEST = PROGRAM / "SHA256SUMS.txt"
CERT_MANIFEST = PROGRAM / "certificates" / "SHA256SUMS.txt"


def tracked_files() -> list[Path]:
    """Certificates, in a stable order. Excludes the manifests themselves."""
    out: list[Path] = []
    for p in sorted((PROGRAM / "certificates").rglob("*.json")):
        out.append(p)
    return out


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lines() -> list[str]:
    return [f"{digest(p)}  {p.relative_to(PROGRAM).as_posix()}" for p in tracked_files()]


def write() -> int:
    body = "\n".join(lines()) + "\n"
    MANIFEST.write_text(body, encoding="utf-8")
    # The certificates/ copy uses the same program-relative paths as the committed file did.
    CERT_MANIFEST.write_text(body, encoding="utf-8")
    print(f"wrote {MANIFEST.relative_to(PROGRAM.parents[1])} and "
          f"{CERT_MANIFEST.relative_to(PROGRAM.parents[1])} ({len(lines())} entries)")
    return 0


def verify() -> int:
    if not MANIFEST.exists():
        print("SHA256SUMS.txt missing; run with --write", file=sys.stderr)
        return 1
    want = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        h, _, rel = line.partition("  ")
        want[rel.strip()] = h.strip()

    stale, missing, extra = [], [], []
    have = {p.relative_to(PROGRAM).as_posix(): digest(p) for p in tracked_files()}
    for rel, h in want.items():
        if rel not in have:
            missing.append(rel)
        elif have[rel] != h:
            stale.append(rel)
    for rel in have:
        if rel not in want:
            extra.append(rel)

    if stale or missing or extra:
        print("SHA256SUMS: STALE", file=sys.stderr)
        for rel in stale:
            print(f"  - changed: {rel}", file=sys.stderr)
        for rel in missing:
            print(f"  - listed but absent: {rel}", file=sys.stderr)
        for rel in extra:
            print(f"  - present but unlisted: {rel}", file=sys.stderr)
        print("\nRegenerate with --write after a certificate run. Every regeneration "
              "rewrites generated_utc, so this file goes stale on any rerun by design.",
              file=sys.stderr)
        return 1
    print(f"SHA256SUMS: OK ({len(want)} entries)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    return write() if args.write else verify()


if __name__ == "__main__":
    raise SystemExit(main())
