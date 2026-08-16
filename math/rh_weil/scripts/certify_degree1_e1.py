#!/usr/bin/env python3
"""WO-RH-10 degree-1 certification entrypoint."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from interval_backend import FlintUnavailable  # noqa: E402
from weil_degree1 import certify_degree1_e1  # noqa: E402


def main() -> int:
    out = ROOT / "certificates" / "e1_degree1_log3_log4.json"
    try:
        body = certify_degree1_e1()
    except FlintUnavailable as exc:
        print(exc, file=sys.stderr)
        return 2
    out.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} evidence_class={body['evidence_class']} status={body['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
