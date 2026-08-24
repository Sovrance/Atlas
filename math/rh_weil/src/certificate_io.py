"""Certificate I/O for RH/Weil (WO-RH-06). Stdlib JSON only."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[1]
CERT_DIR = ROOT / "certificates"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def source_hash(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for p in sorted(paths, key=lambda x: str(x)):
        h.update(p.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def _enforce_quarantine(name: str, body: Dict[str, Any]) -> None:
    """Re-assert the WO-RH-17 quarantine on any affected certificate being written.

    The quarantine marks certificates whose numbers came from the REJECTED even
    pole block. Those files are still produced by the original ``certify_*.py``
    scripts, which know nothing about the adjudication -- so without this guard a
    single re-run would silently restore ``hard_constraints_certified: true`` and
    drop the marker. This is the same failure mode already closed for
    ``work_order_status.json``; the fix belongs at the one write choke point.

    Only ``quarantine_normalization.py --release`` may lift the marker, and it
    does so by passing ``allow_quarantine_change=True``.
    """
    import normalization as N

    if not N.is_quarantined_certificate(name):
        return
    if body.get("promotion_state") == N.QUARANTINE_STATE and "quarantine" in body:
        return  # already carries the marker; leave the recorded prior state alone
    body["quarantine"] = N.quarantine_block(body)
    body["promotion_state"] = N.QUARANTINE_STATE
    body["hard_constraints_certified"] = False
    body["rh_proof_claim"] = False


def write_certificate(
    name: str, body: Dict[str, Any], *, allow_quarantine_change: bool = False
) -> Path:
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    path = CERT_DIR / name
    if not allow_quarantine_change:
        _enforce_quarantine(name, body)
    if "rh_proof_claim" not in body:
        body["rh_proof_claim"] = False
    if "generated_utc" not in body:
        body["generated_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    # Content hash excludes volatile timestamp.
    stable = {k: v for k, v in body.items() if k not in {"generated_utc", "content_hash"}}
    body["content_hash"] = sha256_bytes(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )
    path.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
    return path


def build_e0_exact_certificate() -> Dict[str, Any]:
    import core

    return {
        "certificate_version": "0.1",
        "program": "RH/Weil exact identities",
        "work_order": "WO-RH-01/03/04",
        "evidence_class": "E0",
        "status": "REGENERATED",
        "hard_constraints_certified": True,
        "normalization": core.NORMALIZATION,
        "claim_boundary": core.CLAIM_BOUNDARY,
        "rh_proof_claim": False,
        "checks": [
            "overlap_c / K_ij kernels",
            "midpoint-odd K_q1q1 + sign threshold",
            "bubble even block det",
            "degree-2 parity factorization identities",
        ],
        "note": "Algebraic identities only; not an interval or RH certificate.",
    }


def build_e0_scalar_cell_certificate(report: Dict[str, Any]) -> Dict[str, Any]:
    import core

    return {
        "certificate_version": "0.1",
        "program": "RH/Weil scalar cell [log3, log4]",
        "work_order": "WO-RH-02",
        "evidence_class": "E0",
        "status": "REGENERATED",
        "hard_constraints_certified": bool(report.get("w00_second_positive")),
        "normalization": core.NORMALIZATION,
        "claim_boundary": core.CLAIM_BOUNDARY,
        "rh_proof_claim": False,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "report": report,
        "note": (
            "Algebraic positivity of W00'' on the cell and jump bookkeeping. "
            "Does not promote imported notebook numeric lower bounds."
        ),
    }


def build_e3_fourier_scan_certificate(scan: Dict[str, Any]) -> Dict[str, Any]:
    import core

    return {
        "certificate_version": "0.1",
        "program": "RH/Weil direct-Fourier T=84 probe scan",
        "work_order": "WO-RH-05",
        "evidence_class": "E3",
        "status": "HEURISTIC_SCAN_PENDING_INTERVAL_COVERAGE",
        "hard_constraints_certified": False,
        "normalization": core.NORMALIZATION,
        "claim_boundary": core.CLAIM_BOUNDARY,
        "rh_proof_claim": False,
        "domain": {"fourier_cutoff_T": 84, "L_interval": ["log(3)", "log(4)"]},
        "targets_pending_E1": [
            "E2,84'' > 0 on [log(3), 1.20]",
            "E2,84' > 0 on [1.20, log(4)]",
            "interval point ball near L=1.1059498113",
        ],
        "scan": scan,
        "note": (
            "Stable H0/Hb forms and L-jets are implemented. Uniform interval "
            "coverage of the true Weil E2,84 Gram is NOT closed; this scan is "
            "heuristic energy-probe only and must not be promoted to E1."
        ),
    }


def _active_normalization_id():
    try:
        import normalization

        return normalization.normalization_id()
    except Exception:  # pragma: no cover
        return None


def build_work_order_status() -> Dict[str, Any]:
    return {
        "certificate_version": "0.2",
        "program": "RH/Weil work-order status",
        "rh_proof_claim": False,
        "eng_spec": "ATLAS-RH-ENG-003 / WO-RH-17 normalization adjudication",
        "orders": {
            "WO-RH-01": "done",
            "WO-RH-02": "done_E0_scalar_cell",
            "WO-RH-03": "done",
            "WO-RH-04": "done_algebraic",
            "WO-RH-05": "quarantined_pending_WO-RH-17",
            "WO-RH-06": "done_partial_E0_certs_no_imported_promotion",
            "WO-RH-07": "done_dedicated_runner",
            "WO-RH-08": "unblocked_pending_degree3_implementation",
            "WO-RH-09": "quarantined_pending_WO-RH-17",
            "WO-RH-10": "quarantined_pending_WO-RH-17",
            "WO-RH-11": "quarantined_pending_WO-RH-17",
            "WO-RH-12": "quarantined_pending_WO-RH-17",
            "WO-RH-13": "quarantined_pending_WO-RH-17",
            "WO-RH-14": "quarantined_pending_WO-RH-17",
            "WO-RH-15": "quarantined_pending_WO-RH-17",
            "WO-RH-16": "done_pir_export_partial",
            "WO-RH-17": "done_normalization_adjudicated",
            "WO-RH-18": "done_four_way_crosscheck",
        },
        # Historical values retained verbatim: WO-RH-17 forbids deleting contrary
        # evidence. These are what the tree claimed before the adjudication.
        "pre_quarantine_orders": {
            "WO-RH-05": "done_E1_uniform_true_weil_gram",
            "WO-RH-09": "partial_absolute_G00_regenerated_pending_full_cell_cover",
            "WO-RH-10": "partial_assembly_pending_tail_bound",
            "WO-RH-11": "partial_pole_wired_pending_tight_tail",
            "WO-RH-12": "done_true_weil_gram_with_pole_and_quad_bound",
            "WO-RH-13": "done_E1_T84_points",
            "WO-RH-14": "done_analytic_jets",
            "WO-RH-15": "done_E1_uniform_true_weil_gram",
        },
        "active_normalization_id": _active_normalization_id(),
        "notes": [
            "E3 energy probe quarantined as fourier_energy_probe",
            "No RH proof claim",
            "Even pole outer-product (sqrt(3)/2) REJECTED by WO-RH-17: it equals the "
            "explicit-formula pole times (sqrt(3)/2)cosh(L/2), a calibration fitted at L=log3",
            "Adopted pole: G0_ij = E_i^+E_j^- + E_i^-E_j^+ (see docs/NORMALIZATION_ADJUDICATION_v0.1.md)",
            "WO-RH-05/09..15 quarantined pending regeneration under WO-RH-19/20",
            "Uniform T=84 uses regenerated monotone E2' > 0 (not notebook split)",
            "E1 filename prefixes may hold non-E1 status until gates close",
        ],
    }
