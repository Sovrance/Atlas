#!/usr/bin/env python3
"""Export RH/Weil PIR facts (WO-RH-16)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pir_bridge import export_pir_facts  # noqa: E402


def main() -> int:
    path = export_pir_facts()
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
