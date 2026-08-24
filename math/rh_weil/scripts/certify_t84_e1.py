#!/usr/bin/env python3
"""ENG-005 §6/§7/§8/§10 — fresh Candidate-A T=84 topology, point E1 and uniform E1.

Three artifacts, in dependency order:

  e3_fourier_T84_scan.json          fresh Candidate-A topology scan (§6)
  e1_fourier_T84_points.json        rigorous point balls (§7)
  e1_fourier_T84_uniform_degree2.json   uniform E2_84 lower bound (§10)

The scan comes first because §8 requires the uniform proof topology to be chosen
from it rather than precommitted. The previous Candidate-B monotonicity topology
is not reused; the superseded scan is preserved under
``certificates/history/`` as rejected-normalization provenance.

    python3 scripts/certify_t84_e1.py [--release] [--quick] [--stage scan|points|uniform]
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import e1_t84 as ET  # noqa: E402
import promotion  # noqa: E402
import t84  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from interval_backend import FlintUnavailable, require_flint, set_precision_bits  # noqa: E402

CERT_DIR = ROOT / "certificates"
HISTORY_DIR = CERT_DIR / "history"

SCAN_FILE = "e3_fourier_T84_scan.json"
POINTS_FILE = "e1_fourier_T84_points.json"
UNIFORM_FILE = "e1_fourier_T84_uniform_degree2.json"
INTERIOR_FILE = "e1_fourier_T84_interior_minimum.json"

DEPENDENCIES = (
    "src/pole.py",
    "src/weil_entries.py",
    "src/t84.py",
    "src/e1_t84.py",
    "src/interval_cover.py",
    "src/rigorous_integration.py",
    "src/normalization.py",
    "scripts/certify_t84_e1.py",
)


def preserve_superseded_scan() -> str | None:
    """Keep the Candidate-B scan as provenance rather than overwriting it (§6)."""
    src = CERT_DIR / SCAN_FILE
    if not src.exists():
        return None
    try:
        body = json.loads(src.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if body.get("pole_candidate") == "A":
        return None  # already regenerated; nothing superseded to keep
    HISTORY_DIR.mkdir(parents=True, exist_ok=True)
    dest = HISTORY_DIR / "e3_fourier_T84_scan_candidateB_superseded.json"
    body["superseded_by"] = "ATLAS-RH-ENG-005 §6 fresh Candidate-A scan"
    body["retained_as"] = "rejected-normalization provenance; never a warrant"
    body["pole_candidate"] = "B_REJECTED"
    dest.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
    return str(dest.relative_to(ROOT))


def build_scan(*, n_points: int) -> dict:
    _, arb, acb, _ = require_flint()
    set_precision_bits(ET.DEFAULT_PRECISION_BITS)
    scan = t84.topology_scan(arb, acb, n_points=n_points)
    topo = ET.describe_topology(scan)
    return {
        "certificate_version": "1.0",
        "program": "RH/Weil T=84 topology scan — fresh under Candidate A",
        "work_order": "ATLAS-RH-ENG-005 §6",
        "evidence_class": "E3",
        "status": "SCANNED",
        "hard_constraints_certified": False,
        "rh_proof_claim": False,
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "pole_candidate": "A",
        "jets": "exact support-length jets (§9); no finite differences",
        "scan": scan,
        "topology": topo,
        "note": ("E3 topology evidence. Apparent features are located on a grid; "
                 "nothing here is a warrant. The rejected Candidate-B monotonicity "
                 "topology is not reused."),
    }


def build_points(scan_body: dict, *, quick: bool) -> dict:
    _, arb, acb, _ = require_flint()
    set_precision_bits(ET.DEFAULT_PRECISION_BITS)
    pts = ET.selected_points(scan_body.get("scan"))
    rows = ET.point_rows(pts, arb, acb)
    ok = all(r["E2_definitely_positive"] and r["O1_definitely_positive"] for r in rows)
    body = {
        "certificate_version": "1.0",
        "program": "RH/Weil direct-Fourier T=84 point entries — Candidate A",
        "work_order": "ATLAS-RH-ENG-005 §7",
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if ok else "FAIL",
        "hard_constraints_certified": bool(ok),
        "promotion_state": promotion.PROMOTED_STATE if ok else "REFUSED",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "claim_scope": ET.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "T": t84.T84,
        "pole_candidate": "A",
        "point_scoped": True,
        "claim": ("At each listed L, the true finite Weil T=84 matrix has "
                  "E2 > 0 and O1 > 0, certified by Arb interval balls."),
        "points_selected_from": "the fresh Candidate-A topology scan (§7)",
        "points": rows,
        "quadrature": {"panel_schedule": [list(p) for p in t84.PANELS_T84],
                       "method": "arb_acb_integral_panelled"},
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "precision_bits": ET.DEFAULT_PRECISION_BITS,
        "mpmath_used": False,
        "quick_mode": quick,
    }
    if quick:
        body.update({"promotion_state": "REFUSED", "hard_constraints_certified": False,
                     "status": "QUICK_SMOKE_TEST_NOT_PROMOTABLE"})
    return body


def build_interior(scan_body: dict, *, quick: bool) -> dict:
    """§8 interior-minimum: locate the minimiser, don't merely bound the function."""
    detail = ET.certify_interior_minimum(
        scan=scan_body.get("topology") or scan_body.get("scan"),
        window=4e-5 if quick else ET.STATIONARY_WINDOW,
        progress=lambda m: print(f"  {m}", flush=True),
    )
    ok = (detail["stationary_point"]["sign_change_certified"]
          and float(detail["curvature"]["E2_d2_certified_lower_bound"]) > 0
          and float(detail["basin_bound"]["certified_lower_bound"]) > 0)
    body = {
        "certificate_version": "1.0",
        "program": "RH/Weil T=84 degree-2 interior minimum — Candidate A",
        "work_order": "ATLAS-RH-ENG-005 §8",
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if ok else "FAIL",
        "hard_constraints_certified": bool(ok),
        "promotion_state": promotion.PROMOTED_STATE if ok else "REFUSED",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "claim_scope": ET.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "pole_candidate": "A",
        "claim": ("E2_84 has a unique critical point L* in [log 3, log 4], it is a "
                  "strict minimum, and E2_84 >= the window bound everywhere on the "
                  "closed cell."),
        "topology_source": ("starting bracket taken from the fresh Candidate-A scan "
                            "(E3); every sign and every bound below is certified"),
        "interior_minimum": detail,
        "quadrature": {"panel_schedule": [list(p) for p in t84.PANELS_T84],
                       "method": "arb_acb_integral_panelled"},
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "precision_bits": ET.UNIFORM_PRECISION_BITS,
        "mpmath_used": False,
        "quick_mode": quick,
    }
    if quick:
        body.update({"promotion_state": "REFUSED", "hard_constraints_certified": False,
                     "status": "QUICK_SMOKE_TEST_NOT_PROMOTABLE"})
    return body


def build_uniform(scan_body: dict, *, quick: bool, workers: int = 4) -> dict:
    interior = _interior_summary()
    governed = interior.get("governed_interval")
    exclude = (float(governed[0]), float(governed[1])) if governed else None
    result = ET.certify_uniform_parallel(
        "E2",
        initial_boxes=32 if quick else 896,
        max_depth=6 if quick else 8,
        workers=1 if quick else workers,
        exclude=exclude,
        progress=lambda p: print(
            f"  chunk {p['index']}: bound={p['certified_lower_bound']:.6e} "
            f"boxes={p['boxes_examined']} depth={p['max_depth']}", flush=True),
    )
    ok = result.certified_lower_bound > 0
    topo = ET.describe_topology(scan_body.get("scan", {}))
    body = {
        "certificate_version": "1.0",
        "program": "RH/Weil direct-Fourier T=84 uniform degree-2 — Candidate A",
        "work_order": "ATLAS-RH-ENG-005 §10",
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if ok else "FAIL",
        "hard_constraints_certified": bool(ok),
        "promotion_state": promotion.PROMOTED_STATE if ok else "REFUSED",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "claim_scope": ET.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "T": t84.T84,
        "pole_candidate": "A",
        "domain": {
            "L_interval": list(ET.CELL_LABEL),
            "L_left": repr(ET.CELL[0]),
            "L_right": repr(ET.CELL[1]),
            "closed": True,
        },
        "claim": ("E2_84(L) = G00 Gbb - G0b^2 >= certified_lower_bound > 0 for every "
                  "L in [log 3, log 4], for the truncated T=84 finite Weil matrix."),
        "certified_lower_bound": repr(_combined_bound(result, interior)),
        "bound_provenance": _bound_provenance(result, interior),
        "certified_topology": topo,
        "interval_coverage": result.to_dict(),
        "point_anchors": "see e1_fourier_T84_points.json (§7)",
        "interior_minimum": _interior_summary(),
        "jets": "exact support-length jets (§9); no finite differences",
        "quadrature": {"panel_schedule": [list(p) for p in t84.PANELS_T84],
                       "method": "arb_acb_integral_panelled"},
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "precision_bits": ET.DEFAULT_PRECISION_BITS,
        "mpmath_used": False,
        "quick_mode": quick,
    }
    if quick:
        body.update({"promotion_state": "REFUSED", "hard_constraints_certified": False,
                     "status": "QUICK_SMOKE_TEST_NOT_PROMOTABLE"})
    return body


def _combined_bound(result, interior: dict) -> float:
    """The sharper of the two warrants, where both apply.

    The exhaustive cover alone bounds E2_84 over the whole cell. Where the §8
    interior-minimum argument also applies — on ``governed_interval``, by
    certified derivative and curvature signs — its window bound is sharper,
    because it does not pay the box width of a cover that has to resolve a very
    flat minimum. Outside that interval only the cover applies. So the combined
    bound is the smaller of the interior window bound and the cover's bound over
    the boxes lying wholly outside the governed interval.
    """
    if not interior.get("available") or interior.get("status") != "PASS":
        return result.certified_lower_bound
    outside = result.lower_bound_outside
    if not math.isfinite(outside):
        return result.certified_lower_bound
    combined = min(float(interior["window_bound"]), outside)
    # Never report a bound weaker than the plain cover already proves.
    return max(combined, result.certified_lower_bound)


def _bound_provenance(result, interior: dict) -> dict:
    return {
        "exhaustive_cover_whole_cell": repr(result.certified_lower_bound),
        "exhaustive_cover_outside_governed_interval": (
            repr(result.lower_bound_outside)
            if math.isfinite(result.lower_bound_outside) else None),
        "interior_minimum_window_bound": (interior.get("window_bound")
                                          if interior.get("available") else None),
        "governed_interval": interior.get("governed_interval"),
        "rule": ("min(interior window bound, cover bound outside the governed "
                 "interval), floored at the plain whole-cell cover bound so the "
                 "headline is never weaker than what the cover alone proves"),
    }


def _interior_summary() -> dict:
    """Summarise the §8 interior-minimum certificate, if it has been produced."""
    path = CERT_DIR / INTERIOR_FILE
    if not path.exists():
        return {"available": False,
                "note": f"run --stage interior to produce {INTERIOR_FILE}"}
    body = json.loads(path.read_text(encoding="utf-8"))
    d = body.get("interior_minimum", {})
    return {
        "available": True,
        "certificate": INTERIOR_FILE,
        "status": body.get("status"),
        "minimiser_enclosure": d.get("stationary_point", {}).get("rigorous_interval"),
        "minimiser_approx": d.get("stationary_point", {}).get("approximate_location"),
        "window_bound": d.get("basin_bound", {}).get("certified_lower_bound"),
        "governed_interval": d.get("governed_interval"),
        "E2_d2_lower_bound": d.get("curvature", {}).get("E2_d2_certified_lower_bound"),
        "note": ("independent second warrant: locates the minimiser and proves no "
                 "lower values elsewhere, where the exhaustive cover assumes no "
                 "topology at all. If the two ever disagreed the cover is the warrant."),
    }


def _write(name: str, body: dict, release: bool):
    import normalization as N

    # Only a certificate the WO-RH-17 quarantine actually covers gets a release
    # block. The interior-minimum certificate is new in ENG-005 and was never
    # quarantined, so stamping it "previous_state: QUARANTINED..." would be a
    # false provenance claim about a file that has no such history.
    if release and body.get("status") == "PASS" and N.is_quarantined_certificate(name):
        body["quarantine_released"] = {
            "by": body["work_order"],
            "previous_state": N.QUARANTINE_STATE,
            "authorised_because": "rigorous Candidate-A T=84 regeneration passed",
        }
        return write_certificate(name, body, allow_quarantine_change=True)
    return write_certificate(name, body)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--stage", choices=["scan", "points", "interior", "uniform", "all"],
                    default="all")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--scan-points", type=int, default=17)
    args = ap.parse_args()

    try:
        stages = (["scan", "points", "interior", "uniform"]
                  if args.stage == "all" else [args.stage])

        scan_body = None
        if "scan" in stages:
            kept = preserve_superseded_scan()
            if kept:
                print(f"preserved superseded Candidate-B scan at {kept}")
            scan_body = build_scan(n_points=9 if args.quick else args.scan_points)
            print(f"wrote {write_certificate(SCAN_FILE, scan_body)}")
            print(f"  topology: {scan_body['topology']['classification']} "
                  f"({scan_body['topology']['reason']})")
        else:
            scan_body = json.loads((CERT_DIR / SCAN_FILE).read_text(encoding="utf-8"))

        rc = 0
        if "points" in stages:
            body = build_points(scan_body, quick=args.quick)
            print(f"wrote {_write(POINTS_FILE, body, args.release)}")
            print(f"  points: status={body['status']} promotion_state={body['promotion_state']} "
                  f"n={len(body['points'])}")
            rc = rc or (0 if body["status"] in ("PASS", "QUICK_SMOKE_TEST_NOT_PROMOTABLE") else 1)

        if "interior" in stages:
            body = build_interior(scan_body, quick=args.quick)
            print(f"wrote {_write(INTERIOR_FILE, body, args.release)}")
            d = body["interior_minimum"]
            print(f"  interior: status={body['status']} "
                  f"L*in{d['stationary_point']['rigorous_interval']} "
                  f"bound={d['basin_bound']['certified_lower_bound']}")
            rc = rc or (0 if body["status"] in ("PASS", "QUICK_SMOKE_TEST_NOT_PROMOTABLE") else 1)

        if "uniform" in stages:
            body = build_uniform(scan_body, quick=args.quick, workers=args.workers)
            print(f"wrote {_write(UNIFORM_FILE, body, args.release)}")
            print(f"  uniform: status={body['status']} bound={body['certified_lower_bound']} "
                  f"boxes={body['interval_coverage']['boxes_examined']} "
                  f"depth={body['interval_coverage']['max_subdivision_depth']}")
            rc = rc or (0 if body["status"] in ("PASS", "QUICK_SMOKE_TEST_NOT_PROMOTABLE") else 1)
        return rc
    except FlintUnavailable as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001 - surface the stop condition verbatim
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
