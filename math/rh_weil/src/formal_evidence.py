"""Formal theorem evidence for PIR (ATLAS-RH-ENG-007 §12, WO-RH-44).

The distinction this module exists to keep visible
--------------------------------------------------
A Lean theorem and an Arb interval certificate are different kinds of evidence, and the
composite claim needs both::

    Arb interval enclosure                    E1 / RIGOROUS_COMPUTATION
    Lean: positive lower bounds => PD         FORMAL theorem dependency

    composite warrant = rigorous numerical certificate
                      + formally verified implication

§12 is explicit: "A formal certificate may strengthen an **exact theorem dependency**, but
it does not convert interval numerical evidence to FORMAL." So this module never touches a
certificate's ``evidence_class``. It adds two *separate* fields -- ``numeric_warrant`` and
``logical_implication_warrant`` -- so that a consumer can see which half is which.

The failure being designed against is a certificate that reads FORMAL and is taken to mean
the numbers were formally verified. They were not, and no theorem here could make them so:
Lean proves the implication, Arb produces the enclosures, and the trust boundary between
them is exactly where this module sits.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

#: §12's new content kind.
KIND_FORMAL_THEOREM = "FORMAL_THEOREM_CERTIFICATE"

_PROGRAM = Path(__file__).resolve().parents[1]
MANIFEST_PATH = _PROGRAM / "formal" / "manifests" / "theorem_manifest.json"

#: Which formal theorem underwrites which certificate's *conclusion*.
#:
#: Deliberately an explicit mapping and not a heuristic. A certificate is linked to an
#: implication because a work order says the consumer relies on it, never because the names
#: happen to look related -- the same discipline `quarantine_normalization.RECOVERED_ORDERS`
#: applies to quarantine release.
CERTIFICATE_IMPLICATIONS: Dict[str, tuple[str, ...]] = {
    "e1_degree2_compact_log3_log4.json": ("certificate_even2_implies_pd", "pd_two_by_two"),
    "e1_degree1_log3_log4.json": ("pd_two_by_two",),
    "e1_fourier_T84_uniform_degree2.json": ("certificate_even2_implies_pd", "pd_two_by_two"),
    "e1_degree3_odd_positivity_log3_log4.json": (
        "pd_two_by_two", "det_parity_factorization",
    ),
}


def load_manifest(path: Optional[Path] = None) -> Optional[Dict[str, Any]]:
    """Read the formal theorem manifest, or None when the formal layer is absent.

    Absence is not an error: the Python chain must keep working without a Lean toolchain.
    What must never happen is a certificate *claiming* a formal dependency that cannot be
    resolved, which `formal_dependency_for` handles by returning nothing rather than a
    plausible-looking placeholder.
    """
    p = path or MANIFEST_PATH
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return None


def manifest_id(manifest: Optional[Dict[str, Any]]) -> Optional[str]:
    """A stable identifier for the formal layer a fact was built against.

    Uses the pinned toolchain and Mathlib revision rather than a content hash of the file,
    so that the identifier names the *replayable environment*: two manifests built from the
    same Lean and Mathlib prove the same theorems.
    """
    if not manifest:
        return None
    return f"formal:{manifest.get('lean_toolchain', '?')}+mathlib:{manifest.get('mathlib_rev', '?')[:12]}"


def formal_dependency_for(
    certificate_file: str, manifest: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """The formal-evidence annotation for one certificate.

    Returns ``{}`` when the formal layer is unavailable or the certificate has no declared
    implication. An empty annotation is honest; a partially-filled one would let a consumer
    read a formal warrant that nothing backs.
    """
    manifest = manifest if manifest is not None else load_manifest()
    wanted = CERTIFICATE_IMPLICATIONS.get(certificate_file)
    if not manifest or not wanted:
        return {}
    available = {t["id"]: t for t in manifest.get("theorems", [])}
    resolved = [t for t in wanted if t in available]
    if not resolved:
        return {}
    return {
        "formal_manifest_id": manifest_id(manifest),
        "formal_theorem_ids": resolved,
        "formal_statement_hashes": {t: available[t]["statement_hash"] for t in resolved},
        # The two halves, kept apart on purpose.
        "numeric_warrant": "E1",
        "logical_implication_warrant": "FORMAL",
        "formal_scope_note": (
            "Lean proves the implication from the certified enclosures to the finite "
            "conclusion. It does not verify the enclosures and does not upgrade the "
            "numerical evidence class, which remains E1."
        ),
    }


def formal_theorem_facts_content(manifest: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Content dicts for the formal theorems themselves, one per exported theorem.

    These carry ``numeric_warrant: None`` -- a theorem has no numerical content at all --
    which is what stops a consumer from reading a proved implication as a measurement.
    """
    manifest = manifest if manifest is not None else load_manifest()
    if not manifest:
        return []
    mid = manifest_id(manifest)
    out: List[Dict[str, Any]] = []
    for t in manifest.get("theorems", []):
        out.append({
            "content_kind": KIND_FORMAL_THEOREM,
            "formal_manifest_id": mid,
            "theorem_id": t["id"],
            "trusted_statement": t["trusted_statement"],
            "solution_theorem": t["solution_theorem"],
            "statement": t["statement"],
            "statement_hash": t["statement_hash"],
            "axioms": t["axioms"],
            "lean_toolchain": manifest.get("lean_toolchain"),
            "mathlib_rev": manifest.get("mathlib_rev"),
            "numeric_warrant": None,
            "logical_implication_warrant": "FORMAL",
            "claim_scope": "finite_dimensional_linear_algebra_and_certificate_semantics",
            "rh_proof_claim": False,
            # A proved implication is not a PSD measurement. Spelled out because §11's rule
            # -- a consumer requiring PSD must not be satisfied by the wrong content kind --
            # applies to this kind exactly as it does to inertia.
            "satisfies_psd_requirement": False,
        })
    return out
