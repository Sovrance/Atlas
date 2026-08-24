"""PIR bridge for RH/Weil certificates (WO-RH-16)."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[3]  # repo root
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import promotion

CERT_DIR = Path(__file__).resolve().parents[1] / "certificates"

try:
    import pir
    from pir import AnalyzerRef, Fact, Warning_
    from pir.canonical import content_id

    _PIR = True
except Exception:  # pragma: no cover
    _PIR = False


def available() -> bool:
    return _PIR


def _load(name: str) -> Dict[str, Any] | None:
    p = CERT_DIR / name
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


# ENG-004 §3: the promotion rules live in one place. ``pir_bridge`` asks that
# module rather than carrying its own copy, so PIR can never promote something
# the runner or the release path would refuse.
QUARANTINE_PREFIX = promotion.QUARANTINE_PREFIX
promotion_refusal = promotion.promotion_refusal


def _active_normalization_id() -> str | None:
    """Active id, read from the adjudication artifact (never from a filename)."""
    return promotion.active_normalization_id()


def refused_promotions() -> List[Dict[str, Any]]:
    """Every certificate the central predicate currently refuses to promote."""
    return promotion.refused_promotions(CERT_DIR)


def certs_to_facts() -> List[Any]:
    if not _PIR:
        raise RuntimeError("pir unavailable")
    facts = []
    mapping = [
        ("e0_exact_identities.json", "E0", "SOUND", None),
        ("e0_scalar_cell_log3_log4.json", "E0", "SOUND", None),
        ("normalization_adjudication.json", "E0", "SOUND", None),
        ("normalization_crosscheck.json", "E2", "HEURISTIC",
         "three-way internal cross-check; only interval_certified rows may support E1"),
        # ENG-005 §11: the recovered Candidate-A numerical chain. Each of these
        # still has to clear the central promotion predicate before it becomes a
        # fact; being listed here only makes it eligible.
        ("e1_scalar_log3_log4.json", "E1", "SOUND",
         "scalar canary: rigorous Arb lower bound on the recovered Candidate-A cell entry"),
        ("e1_degree1_log3_log4.json", "E1", "SOUND",
         "degree-1 odd pivot O1, cutoff-free, uniform on [log3, log4]"),
        ("e1_degree2_compact_log3_log4.json", "E1", "SOUND",
         "compact degree-2 even determinant E2, cutoff-free, uniform on [log3, log4]"),
        ("e1_fourier_T84_points.json", "E1", "SOUND",
         "direct-Fourier T=84 entries at selected L — point-scoped, not a cell claim"),
        ("e1_fourier_T84_interior_minimum.json", "E1", "SOUND",
         "direct-Fourier T=84 interior minimum: minimiser enclosure, curvature and "
         "derivative-sign inequalities on the governed interval"),
        ("e1_fourier_T84_uniform_degree2.json", "E1", "SOUND",
         "direct-Fourier T=84 uniform degree-2 lower bound on [log3, log4]"),
        ("e3_fourier_T84_scan.json", "E3", "HEURISTIC",
         "fresh Candidate-A T=84 topology scan — E3 evidence, never a warrant"),
        ("external/connes_cvs_crossvalidation_v0.1.json", "E3", "HEURISTIC", "external cross-check only"),
    ]
    for fname, ev, tag, warn in mapping:
        cert = _load(fname)
        if cert is None:
            continue
        refusal = promotion_refusal(cert)
        if refusal:  # WO-RH-17 guard: never promote a quarantined/stale certificate
            continue
        warnings = ()
        if warn:
            warnings = (Warning_(location=fname, message=warn),)
        content = {
            "certificate_file": fname,
            "point_scoped": bool(cert.get("point_scoped", False)),
            "program": cert.get("program"),
            "status": cert.get("status"),
            "claim_scope": "finite_dimensional_weil_compression",
            "rh_proof_claim": False,
            "evidence_class_declared": cert.get("evidence_class", ev),
        }
        analyzer = AnalyzerRef(id=f"rh_weil.{Path(fname).stem}", version="0.2.0", tag=tag)
        _nid = _active_normalization_id()
        assumptions = (
            "asm:normalization:G=G0-Gp+Ginf",
            "asm:domain:L_in_[log3,log4]",
        ) + ((f"asm:normalization_id:{_nid}",) if _nid else ())
        fid = content_id(
            "fct",
            {
                "content": content,
                "analyzer": analyzer.to_dict(),
                "assumptions": list(assumptions),
            },
        )
        facts.append(
            Fact(
                fact_id=fid,
                pir_level="L2",
                evidence_level=ev,
                layer="UNIVERSAL",
                namespace="invariant",
                status="SUPPORTED",
                analyzer=analyzer,
                content=content,
                created_at=cert.get("generated_utc", "1970-01-01T00:00:00Z"),
                assumptions=assumptions,
                warnings=warnings,
                measurement_interface=(),
                verdict=None,
            )
        )
    return facts


def export_pir_facts(path: Path | None = None) -> Path:
    path = path or (CERT_DIR / "pir_facts.json")
    if not _PIR:
        payload = {
            "status": "PIR_UNAVAILABLE",
            "rh_proof_claim": False,
            "facts": [],
            "note": "pir package not importable",
        }
    else:
        facts = certs_to_facts()
        payload = {
            "certificate_version": "0.2",
            "program": "RH/Weil PIR export",
            "work_order": "WO-RH-16",
            "rh_proof_claim": False,
            "active_normalization_id": _active_normalization_id(),
            "n_facts": len(facts),
            "refused_promotions": refused_promotions(),
            "facts": [
                {
                    "fact_id": f.fact_id,
                    "pir_level": f.pir_level,
                    "evidence_level": f.evidence_level,
                    "analyzer": f.analyzer.to_dict(),
                    "content": f.content,
                    "assumptions": list(f.assumptions),
                    "warnings": [w.to_dict() for w in f.warnings],
                }
                for f in facts
            ],
        }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path
