"""Central promotion predicate for RH/Weil certificates (ATLAS-RH-ENG-004 §3).

One place decides whether a certificate may be promoted. Everything else --
``pir_bridge``, the runner, the release path -- asks this module rather than
re-implementing the rules, so a rule can never hold in one caller and not
another.

A certificate is refused when any of these hold:

* it is quarantined (``promotion_state`` starts with ``QUARANTINED``);
* it is *rigorous* (declares E1 or sets ``rigorous: true``) but carries no
  ``normalization_certificate_id``;
* its ``normalization_certificate_id`` does not match the **active** id;
* a recorded source/dependency hash no longer matches the file on disk;
* it declares E1 without ``hard_constraints_certified``.

The active id is read from the adjudication artifact, never inferred from a
filename (§3). It is additionally cross-checked against the content id computed
from the frozen normalization definition: a disagreement means the normalization
changed underneath the artifact, which is an ENG-004 §14 stop condition rather
than something to paper over.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
CERT_DIR = ROOT / "certificates"
SRC_DIR = ROOT / "src"

ADJUDICATION_CERTIFICATE = "normalization_adjudication.json"
NORMALIZATION_ID_FIELD = "normalization_certificate_id"
QUARANTINE_PREFIX = "QUARANTINED"

#: Promotion states that mean "this certificate is allowed to carry a claim".
PROMOTED_STATE = "PROMOTED"


class NormalizationUnavailable(RuntimeError):
    """The adjudication artifact is missing or does not declare an active id."""


# --------------------------------------------------------------------------- #
# Active normalization id                                                      #
# --------------------------------------------------------------------------- #
def _read_adjudication() -> Optional[Dict[str, Any]]:
    p = CERT_DIR / ADJUDICATION_CERTIFICATE
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):  # pragma: no cover - corrupt artifact
        return None


def active_normalization_id(*, strict: bool = False) -> Optional[str]:
    """The active normalization content id, read from the adjudication artifact.

    Never derived from a filename. When the artifact is present its declared
    ``active_normalization_id`` wins; ``strict=True`` raises instead of
    returning ``None`` when it cannot be read.
    """
    body = _read_adjudication()
    if body is not None:
        declared = body.get("active_normalization_id")
        if declared:
            return str(declared)
    if strict:
        raise NormalizationUnavailable(
            f"{ADJUDICATION_CERTIFICATE} missing or has no active_normalization_id; "
            "run scripts/derive_normalization.py"
        )
    return None


def computed_normalization_id() -> Optional[str]:
    """The content id computed from the frozen normalization definition."""
    try:
        import normalization as _N

        return _N.normalization_id()
    except Exception:  # pragma: no cover - normalization import failure
        return None


def normalization_id_consistent() -> Tuple[bool, str]:
    """Does the artifact's declared id still match the code's content id?

    A mismatch is ENG-004 §14 ("normalization ID changes unexpectedly"): the
    frozen definition moved without the adjudication being re-run, so every
    dependent numerical certificate is stale.
    """
    declared = active_normalization_id()
    computed = computed_normalization_id()
    if declared is None:
        return False, "no active_normalization_id in the adjudication artifact"
    if computed is None:
        return False, "normalization module could not compute a content id"
    if declared != computed:
        return False, f"artifact id {declared} != computed id {computed}"
    return True, f"consistent ({declared})"


# --------------------------------------------------------------------------- #
# Source / dependency hashes                                                   #
# --------------------------------------------------------------------------- #
def file_sha256(relpath: str) -> Optional[str]:
    """sha256 of a repo-relative path under ``math/rh_weil``."""
    p = ROOT / relpath
    if not p.exists():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_hashes(relpaths) -> Dict[str, str]:
    """Record ``{relpath: sha256}`` for the modules a certificate depends on."""
    out: Dict[str, str] = {}
    for rel in sorted(relpaths):
        h = file_sha256(rel)
        if h is None:
            raise FileNotFoundError(f"cannot hash missing dependency {rel!r}")
        out[rel] = h
    return out


def stale_dependencies(cert: Dict[str, Any]) -> List[str]:
    """Recorded dependencies whose file no longer hashes to the recorded value."""
    recorded = (cert.get("dependencies") or {}).get("source_hashes") or {}
    bad = []
    for rel, want in sorted(recorded.items()):
        got = file_sha256(rel)
        if got is None:
            bad.append(f"{rel} (missing)")
        elif got != want:
            bad.append(rel)
    return bad


# --------------------------------------------------------------------------- #
# The predicate                                                                #
# --------------------------------------------------------------------------- #
def is_rigorous(cert: Dict[str, Any]) -> bool:
    """Does this certificate carry a rigorous numerical claim (E1 or better)?"""
    if cert.get("rigorous") is True:
        return True
    return str(cert.get("evidence_class", "") or "").startswith("E1")


def promotion_refusal(cert: Dict[str, Any]) -> Optional[str]:
    """Why this certificate may NOT be promoted, or ``None`` if it may.

    Order matters only for which reason is reported first; every rule is
    independently sufficient to refuse.
    """
    state = str(cert.get("promotion_state", "") or "")
    if state.startswith(QUARANTINE_PREFIX):
        return f"promotion_state={state}"

    declared = str(cert.get("evidence_class", "") or "")
    if declared.startswith("E1") and not cert.get("hard_constraints_certified", False):
        return "declares E1 but hard_constraints_certified is false"

    active = active_normalization_id()
    cert_id = cert.get(NORMALIZATION_ID_FIELD)
    # Legacy artifacts recorded the id under other names; accept them for the
    # match test but never as a substitute for the required field on E1.
    legacy_id = cert.get("active_normalization_id") or cert.get("normalization_id")

    if is_rigorous(cert):
        if not cert_id:
            return f"rigorous certificate is missing {NORMALIZATION_ID_FIELD}"
        if active is None:
            return "no active normalization id available to validate against"
        if cert_id != active:
            return f"stale {NORMALIZATION_ID_FIELD}={cert_id} (active {active})"
        ok, why = normalization_id_consistent()
        if not ok:
            return f"normalization id inconsistent: {why}"
    elif legacy_id and active and legacy_id != active:
        return f"stale normalization_id={legacy_id}"

    bad = stale_dependencies(cert)
    if bad:
        return "stale dependency hashes: " + ", ".join(bad)

    return None


def may_promote(cert: Dict[str, Any]) -> bool:
    return promotion_refusal(cert) is None


def refused_promotions(cert_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Every certificate under ``cert_dir`` the predicate currently refuses."""
    cert_dir = cert_dir or CERT_DIR
    out = []
    for p in sorted(cert_dir.glob("*.json")):
        try:
            cert = json.loads(p.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        reason = promotion_refusal(cert)
        if reason:
            out.append(
                {
                    "certificate_file": p.name,
                    "reason": reason,
                    "evidence_class_declared": cert.get("evidence_class"),
                }
            )
    return out
