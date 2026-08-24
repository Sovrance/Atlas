#!/usr/bin/env python3
"""WO-RH-17 §3.4 — quarantine certificates that depend on the disputed pole block.

Marks the affected certificates non-promotable **without deleting anything and
without rewriting their claimed evidence class**:

    promotion_state = "QUARANTINED_NORMALIZATION_ADJUDICATION"
    hard_constraints_certified = false
    rh_proof_claim = false

The prior values are preserved verbatim under ``quarantine.prior_state`` (and the
prior content hash under ``quarantine.pre_quarantine_content_hash``) so the
historical claim remains auditable. Certificates are *not* relabelled E3: the
record of what was once claimed is evidence too.

Also flips the affected work orders in ``work_order_status.json`` to
``quarantined_pending_WO-RH-17``.

Idempotent. Usage::

    python3 scripts/quarantine_normalization.py [--release]   # --release only after WO-RH-19
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import normalization as N  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

CERT_DIR = ROOT / "certificates"

# Single source of truth lives in src/normalization.py so that the certificate
# writer can enforce the same list at the point of write.
QUARANTINE_STATE = N.QUARANTINE_STATE
AFFECTED = list(N.QUARANTINED_CERTIFICATES)

AFFECTED_ORDERS = ["WO-RH-05"] + [f"WO-RH-{n:02d}" for n in range(9, 16)]

REASON = N.QUARANTINE_REASON


def quarantine_certificate(path: Path, release: bool = False) -> str:
    body = json.loads(path.read_text(encoding="utf-8"))
    if release:
        if body.get("promotion_state") != QUARANTINE_STATE:
            return "not-quarantined"
        prior = body.get("quarantine", {}).get("prior_state", {})
        body.pop("quarantine", None)
        body["promotion_state"] = "RELEASED_PENDING_REGENERATION"
        for k, v in prior.items():
            body[k] = v
        write_certificate(path.name, body, allow_quarantine_change=True)
        return "released"

    if body.get("promotion_state") == QUARANTINE_STATE:
        return "already-quarantined"

    # prior_state (the claim being suspended) is captured by quarantine_block.
    body["quarantine"] = N.quarantine_block(body)
    body["promotion_state"] = QUARANTINE_STATE
    body["hard_constraints_certified"] = False
    body["rh_proof_claim"] = False
    write_certificate(path.name, body)
    return "quarantined"


def update_work_order_status(release: bool = False) -> str:
    path = CERT_DIR / "work_order_status.json"
    body = json.loads(path.read_text(encoding="utf-8"))
    orders = body.setdefault("orders", {})
    prior = body.setdefault("pre_quarantine_orders", {})
    changed = 0
    for wo in AFFECTED_ORDERS:
        if wo not in orders:
            continue
        if release:
            if wo in prior:
                orders[wo] = prior.pop(wo)
                changed += 1
            continue
        if orders[wo] != "quarantined_pending_WO-RH-17":
            prior.setdefault(wo, orders[wo])
            orders[wo] = "quarantined_pending_WO-RH-17"
            changed += 1
    if not prior:
        body.pop("pre_quarantine_orders", None)
    orders.setdefault("WO-RH-17", "done_normalization_adjudicated")
    orders.setdefault("WO-RH-18", "done_four_way_crosscheck")
    notes = body.setdefault("notes", [])
    marker = ("Even pole outer-product (sqrt(3)/2) REJECTED by WO-RH-17: it is the "
              "explicit-formula pole times (sqrt(3)/2)cosh(L/2), a calibration fitted at L=log3")
    if marker not in notes:
        notes.append(marker)
    body["active_normalization_id"] = N.normalization_id()
    write_certificate(path.name, body)
    return f"{changed} order(s) updated"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true",
                    help="lift the quarantine (only legitimate after WO-RH-19/20 regeneration)")
    args = ap.parse_args()

    print(f"{'releasing' if args.release else 'quarantining'} certificates:")
    for name in AFFECTED:
        p = CERT_DIR / name
        if not p.exists():
            print(f"  [skip] {name} (absent)")
            continue
        print(f"  [{quarantine_certificate(p, args.release)}] {name}")
    print(f"work_order_status.json: {update_work_order_status(args.release)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
