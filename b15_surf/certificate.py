"""B15 certificate builder / validator (v0.3 shape + work-order §5 fields).

* Shape: existing v0.3 keys (``certificate_version``, ``certificate_class``,
  ``problem``, ``timestamp_utc``, ``headline``, ``m_layer_stipulations``,
  ``calibration_route``, ``results``) PLUS the §5 fields; written through
  ``b1_moment_solver.certificate.save_certificate`` (hard constraint 3).
* ``certificate_id`` is content-addressed as in b13_cdl/docs/pir-bridge-v0.1.md:
  sha256 of the canonical body with the id and ``timestamp_utc`` excluded.
* ``validate`` = schema (``pir/jsonschema_mini``) + the §5 semantic rules that a
  JSON schema cannot express: HEURISTIC may not assert E0; a POS ``REJECTED``
  must carry an ``impossibility_certificate``; a HEURISTIC certificate must
  carry located warnings; ``asm:`` prefix on every assumption.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional

from b1_moment_solver.certificate import save_certificate
from pir.canonical import canonical_json
from pir.jsonschema_mini import SchemaError, validate as _schema_validate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCHEMA_PATH = os.path.join(ROOT, "schemas", "b15_certificate.schema.json")
PREREG_PATH = os.path.join(ROOT, "docs", "preregistrations",
                           "prereg-002-surfaceology-benchmarks.md")
PREREG_FREEZE = os.path.join(ROOT, "docs", "preregistrations", "prereg-002.freeze")

SPEC_VERDICTS = ("FORCED", "PERMITTED", "REJECTED", "NONIDENTIFIABLE",
                 "OBSERVATIONALLY_EQUIVALENT", "APPARATUS_LIMITED",
                 "REPRESENTATION_DEPENDENT", "AMBIGUOUS")


class B15CertificateError(Exception):
    pass


def load_schema() -> Dict:
    with open(SCHEMA_PATH) as f:
        return json.load(f)


def prereg_ref() -> Dict[str, str]:
    """``prereg-002@<commit>`` from the freeze file + content hash of the prereg."""
    with open(PREREG_PATH, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    commit = "UNFROZEN"
    if os.path.exists(PREREG_FREEZE):
        with open(PREREG_FREEZE) as f:
            for line in f:
                if line.startswith("commit="):
                    commit = line.strip().split("=", 1)[1]
    return {"prereg_ref": f"prereg-002@{commit}", "prereg_sha256": sha}


def content_hash(cert: Dict) -> str:
    body = {k: v for k, v in cert.items() if k not in ("certificate_id", "timestamp_utc")}
    return hashlib.sha256(canonical_json(body).encode("utf-8")).hexdigest()[:12]


def verdict_display(verdict: str, cause: Optional[str] = None,
                    cls: Optional[List[str]] = None) -> str:
    if verdict == "NONIDENTIFIABLE" and cause:
        return f"NONIDENTIFIABLE({cause})"
    if verdict == "OBSERVATIONALLY_EQUIVALENT" and cls:
        return f"OBSERVATIONALLY_EQUIVALENT([{', '.join(cls)}])"
    return verdict


def build(*, benchmark: str, problem: str, headline: str, certificate_class: str,
          results: Dict, soundness: str, warnings: List[Dict[str, str]],
          ground_truth_route: str, evidence_level: str, pir_level: str,
          assumptions: List[str], falsifier_direction: str, seed: int,
          generator_sha256: str, m_layer_stipulations: List[str],
          calibration_route: str, verdict: str, witness: Optional[Dict],
          impossibility_certificate: Optional[Dict], inputs: Optional[Dict] = None,
          verdict_cause: Optional[str] = None, verdict_class: Optional[List[str]] = None,
          similarity: Optional[str] = None, confidence: Optional[str] = None,
          correlator: Optional[str] = None, pir_facts: Optional[List[Dict]] = None) -> Dict:
    hard_ok = all(
        r.get("status") in ("PASS", "FORCED", "PERMITTED", "PD_CERTIFIED", "PSD_CERTIFIED")
        for r in results.values() if isinstance(r, dict) and "status" in r)
    pr = prereg_ref()
    cert: Dict[str, Any] = {
        "certificate_version": "0.3",
        "certificate_class": certificate_class,
        "benchmark": benchmark,
        "problem": problem,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "headline": headline,
        "arithmetic": "exact rational (Fraction) on every certified path; no floats",
        "hard_constraints_certified": hard_ok,
        "score": "finite" if hard_ok else "infinity (violated or uncertified)",
        "m_layer_stipulations": list(m_layer_stipulations),
        "calibration_route": calibration_route,
        "inputs": inputs or {},
        "results": results,
        "soundness": soundness,
        "warnings": list(warnings),
        "ground_truth_route": ground_truth_route,
        "verdict": verdict,
        "verdict_cause": verdict_cause,
        "verdict_class": verdict_class,
        "verdict_display": verdict_display(verdict, verdict_cause, verdict_class),
        "witness": witness,
        "impossibility_certificate": impossibility_certificate,
        "similarity": similarity,
        "confidence": confidence,
        "correlator": correlator,
        "evidence_level": evidence_level,
        "pir_level": pir_level,
        "layer": "DOMAIN",
        "assumptions": list(assumptions),
        "prereg_ref": pr["prereg_ref"],
        "prereg_sha256": pr["prereg_sha256"],
        "falsifier_direction": falsifier_direction,
        "seed": int(seed),
        "ci_seed_env": os.environ.get("PIR_CI_SEED", ""),
        "python_version": sys.version.split()[0],
        "generator_sha256": generator_sha256,
        "pir_facts": pir_facts or [],
    }
    cert["certificate_id"] = f"b15-{benchmark.split('-')[1].lower()}-{content_hash(cert)}"
    validate(cert)
    return cert


def validate(cert: Dict) -> None:
    """Schema + §5 semantic rules. Raises B15CertificateError."""
    try:
        _schema_validate(load_schema(), cert)
    except SchemaError as e:
        raise B15CertificateError(f"schema: {e}") from e
    if cert["verdict"] not in SPEC_VERDICTS:
        raise B15CertificateError("non-SPEC verdict")
    if cert["soundness"] == "HEURISTIC" and cert["evidence_level"] == "E0":
        raise B15CertificateError("a HEURISTIC certificate may not assert E0 (hard constraint 4)")
    if cert["soundness"] == "HEURISTIC" and not cert["warnings"]:
        raise B15CertificateError("HEURISTIC requires located warnings[] (hard constraint 4)")
    if cert.get("benchmark") == "B15-POS" and cert["verdict"] == "REJECTED" \
            and not cert.get("impossibility_certificate"):
        raise B15CertificateError("POS REJECTED requires impossibility_certificate (§5)")
    if cert["verdict"] in ("PERMITTED", "REJECTED", "FORCED") and cert["witness"] is None:
        raise B15CertificateError(f"{cert['verdict']} requires a witness (§5)")
    if cert["verdict"] == "NONIDENTIFIABLE" and not cert.get("verdict_cause"):
        raise B15CertificateError("NONIDENTIFIABLE requires a stated cause (SPEC §4)")
    if (cert["similarity"] is not None or cert["confidence"] is not None) and not cert["correlator"]:
        raise B15CertificateError("similarity/confidence require a named correlator")
    for a in cert["assumptions"]:
        if not a.startswith("asm:"):
            raise B15CertificateError(f"assumption {a!r} must be asm:-prefixed (ADR-0002)")
    if not cert["prereg_ref"].startswith("prereg-002@"):
        raise B15CertificateError("prereg_ref must be prereg-002@<commit>")
    expected = content_hash(cert)
    if not cert["certificate_id"].endswith(expected):
        raise B15CertificateError("certificate_id does not match content hash")


def save(cert: Dict, name: str) -> str:
    validate(cert)
    path = os.path.join(ROOT, "certificates", name)
    save_certificate(cert, path)
    return path
