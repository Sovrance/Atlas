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


QUARANTINE_PREFIX = "QUARANTINED"


def promotion_refusal(cert: Dict[str, Any]) -> str | None:
    """WO-RH-17 promotion guard: why this certificate may NOT be promoted.

    Returns a reason string, or ``None`` when promotion is permitted. A
    certificate under normalization quarantine is never promoted, whatever its
    historical evidence label says.
    """
    state = str(cert.get("promotion_state", "") or "")
    if state.startswith(QUARANTINE_PREFIX):
        return f"promotion_state={state}"
    declared = str(cert.get("evidence_class", "") or "")
    if declared.startswith("E1") and not cert.get("hard_constraints_certified", False):
        return "declares E1 but hard_constraints_certified is false"
    norm_id = cert.get("active_normalization_id") or cert.get("normalization_id")
    if norm_id is not None and _active_normalization_id() and norm_id != _active_normalization_id():
        return f"stale normalization_id={norm_id}"
    return None


def _active_normalization_id() -> str | None:
    try:
        import normalization as _N

        return _N.normalization_id()
    except Exception:
        return None


def refused_promotions() -> List[Dict[str, Any]]:
    """Every certificate the guard currently refuses to promote."""
    out = []
    for p in sorted(CERT_DIR.glob("*.json")):
        try:
            cert = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        reason = promotion_refusal(cert)
        if reason:
            out.append({"certificate_file": p.name, "reason": reason,
                        "evidence_class_declared": cert.get("evidence_class")})
    return out


def certs_to_facts() -> List[Any]:
    if not _PIR:
        raise RuntimeError("pir unavailable")
    facts = []
    mapping = [
        ("e0_exact_identities.json", "E0", "SOUND", None),
        ("e0_scalar_cell_log3_log4.json", "E0", "SOUND", None),
        ("normalization_adjudication.json", "E0", "SOUND", None),
        ("normalization_crosscheck.json", "E2", "HEURISTIC",
         "numeric four-way cross-check; only interval_certified rows may support E1"),
        ("e3_fourier_T84_scan.json", "E3", "HEURISTIC", "E3 energy probe — not the true Weil Gram"),
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
