#!/usr/bin/env python3
"""WO-RH-11 compact degree-2 certification entrypoint."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from certificate_io import write_certificate  # noqa: E402
from interval_backend import FlintUnavailable  # noqa: E402
from weil_degree2 import certify_degree2_compact_e1  # noqa: E402


def main() -> int:
    out = ROOT / "certificates" / "e1_degree2_compact_log3_log4.json"
    try:
        body = certify_degree2_compact_e1()
    except FlintUnavailable as exc:
        print(exc, file=sys.stderr)
        return 2
    # Written through certificate_io so the WO-RH-17 quarantine guard applies:
    # this file is regenerated from the REJECTED even pole block and must not
    # come back promotable just because the script was re-run.
    write_certificate(out.name, body)
    print(f"wrote {out} evidence_class={body['evidence_class']} status={body['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
