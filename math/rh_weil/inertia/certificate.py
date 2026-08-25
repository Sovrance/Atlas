"""Inertia certificate construction (§3, §11, §14).

An inertia certificate is a *different content kind* from a positivity
certificate, and the difference is load-bearing. "This block has inertia
(1, 1, 0)" is a complete, useful, certified statement, and it is also a
statement that the block is **not** positive semidefinite. A consumer that
requires PSD must never be satisfied by one. :func:`satisfies_psd_requirement`
is the single place that decision is made, so no downstream reader has to
re-derive it -- and it answers ``True`` only for a signature that is literally
PSD, with the zero count known.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

CLAIM_SCOPE = "finite_dimensional_weil_compression"

#: PIR content kinds introduced by ENG-006 (§11).
KIND_INERTIA = "WEIL_INERTIA_CERTIFICATE"
KIND_STRATIFICATION = "WEIL_INERTIA_STRATIFICATION"


#: Content kinds that are *inertia* artifacts. §11 is categorical about these:
#: they may never satisfy a consumer requiring PSD, whatever their signature.
INERTIA_KINDS = (KIND_INERTIA, KIND_STRATIFICATION)

#: Kinds that can never satisfy a PSD-requiring consumer, whatever their fields say.
#: Inertia kinds know a signature; a formal theorem knows an implication. Neither has
#: measured the block in front of the consumer.
NON_MEASUREMENT_KINDS = frozenset(INERTIA_KINDS) | frozenset({"FORMAL_THEOREM_CERTIFICATE"})


def satisfies_psd_requirement(cert: Dict[str, Any]) -> bool:
    """True only if this certificate *claims* PSD and its signature backs it.

    §11: "An inertia certificate must never satisfy a consumer that explicitly
    requires PSD." That is categorical and it binds on the content kind, not on
    how favourable the signature happens to be -- an inertia artifact with
    signature ``(n, 0, 0)`` is refused exactly like an indefinite one. The point
    of the rule is that "I know the signature" must not be silently read as "it
    is positive"; a consumer wanting positivity should be handed something that
    claims positivity.

    So two independent conditions have to hold. The certificate must not be an
    inertia kind, and it must *say* it is positive via ``psd_claim`` -- an
    explicit declaration by the producer rather than an inference drawn here
    from fields the producer never meant that way. The signature is then checked
    against that claim: zero negative directions, with the zero multiplicity
    known, since an unresolved ``n_zero`` leaves open a negative direction
    hiding in the part that did not resolve.

    ENG-007 §12 added a third refused kind for the same reason, not a new one.
    A ``FORMAL_THEOREM_CERTIFICATE`` records that Lean proved an implication; it
    carries no interval evidence about any particular matrix. It is about
    positive definiteness, which is exactly what makes it dangerous here -- a
    body naming ``pd_two_by_two`` with ``psd_claim: true`` reads as favourable
    and would have satisfied a PSD consumer while measuring nothing. The rule is
    the same one the inertia kinds are refused under: knowing a theorem, like
    knowing a signature, is not the same as having measured this block.

    An earlier version refused only stratifications and inferred the rest from
    the signature. That let a passing ``WEIL_INERTIA_CERTIFICATE`` with
    ``(2, 0, 0)`` satisfy a PSD consumer while its own body said
    ``psd_claim: false`` -- the predicate contradicting the certificate it was
    reading.
    """
    if cert.get("rh_proof_claim") is not False:
        return False
    if cert.get("content_kind") in NON_MEASUREMENT_KINDS:
        return False
    if cert.get("psd_claim") is not True:
        return False
    if cert.get("status") != "PASS":
        return False
    if cert.get("evidence_class") not in ("E0", "E1"):
        return False
    n_neg, n_zero = cert.get("n_negative"), cert.get("n_zero")
    if n_neg is None or n_zero is None:
        return False
    return int(n_neg) == 0


def build_inertia_certificate(
    result,
    *,
    dimension: int,
    program: str,
    work_order: str,
    evidence_class: str,
    normalization_certificate_id: str,
    source_hashes: Optional[Dict[str, str]] = None,
    parameter_domain: Optional[Dict[str, Any]] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """A single-matrix inertia certificate matching ``inertia_certificate.schema.json``."""
    body: Dict[str, Any] = {
        "certificate_version": "1.0",
        "content_kind": KIND_INERTIA,
        "program": program,
        "work_order": work_order,
        "status": result.status if result.status in ("PASS", "FAIL") else "INCONCLUSIVE",
        "claim_scope": CLAIM_SCOPE,
        "dimension": int(dimension),
        "n_positive": result.n_positive,
        "n_negative": result.n_negative,
        "n_zero": result.n_zero,
        "method": result.method,
        "sign_oracle": result.oracle,
        "parameter_domain": parameter_domain,
        "pivot_intervals": [
            {"step": p.step, "kind": p.kind, "indices": list(p.indices),
             "signature": list(p.sign), "value": p.value}
            for p in result.pivots
        ],
        "transition_regions": [],
        "evidence_class": evidence_class,
        "normalization_certificate_id": normalization_certificate_id,
        "rh_proof_claim": False,
        "psd_claim": False,
        "note": (
            "Inertia is not positivity. This certificate states a signature; it "
            "satisfies a PSD requirement only when n_negative == 0 with n_zero "
            "known, which is decided by inertia.certificate.satisfies_psd_requirement."
        ),
    }
    if result.blocker:
        body["blocker"] = result.blocker
    if source_hashes:
        body["source_hashes"] = source_hashes
    if extra:
        body.update(extra)
    return body


def build_stratification_certificate(
    strat,
    *,
    dimension: int,
    program: str,
    work_order: str,
    evidence_class: str,
    normalization_certificate_id: str,
    source_hashes: Optional[Dict[str, str]] = None,
    parameter_domain: Optional[Dict[str, Any]] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """An inertia-stratification certificate over a parameter cell."""
    d = strat.to_dict()
    constant = strat.signature_if_constant()
    body: Dict[str, Any] = {
        "certificate_version": "1.0",
        "content_kind": KIND_STRATIFICATION,
        "program": program,
        "work_order": work_order,
        "status": "PASS" if strat.status.startswith("PASS") else strat.status,
        "claim_scope": CLAIM_SCOPE,
        "dimension": int(dimension),
        # A stratification reports a single signature only when it proved the
        # same one everywhere with nothing left over.
        "n_positive": constant[0] if constant else None,
        "n_negative": constant[1] if constant else None,
        "n_zero": constant[2] if constant else None,
        "constant_inertia": strat.is_constant,
        "method": d["method"],
        "parameter_domain": parameter_domain or {"cell": d["cell"]},
        "strata": d["strata"],
        "transition_regions": d["transition_regions"],
        "coverage": d["coverage"],
        "boxes_examined": d["boxes_examined"],
        "max_subdivision_depth": d["max_subdivision_depth"],
        "subdivision_policy": d["subdivision_policy"],
        "pivot_intervals": [],
        "evidence_class": evidence_class,
        "normalization_certificate_id": normalization_certificate_id,
        "rh_proof_claim": False,
        "psd_claim": False,
        "note": (
            "A stratification never satisfies a PSD requirement, even when every "
            "stratum is PSD: the transition regions between strata are unresolved "
            "by construction, and a negative direction could live in one."
        ),
    }
    if source_hashes:
        body["source_hashes"] = source_hashes
    if extra:
        body.update(extra)
    return body


def validate_against_schema(body: Dict[str, Any], schema: Dict[str, Any]) -> List[str]:
    """Minimal structural validation — required keys, enums, and consts.

    Deliberately dependency-free: the fast CI gate must run without jsonschema
    installed, and the schemas here use only the handful of keywords checked.
    """
    errors: List[str] = []
    for key in schema.get("required", []):
        if key not in body:
            errors.append(f"missing required property {key!r}")
    props = schema.get("properties", {})
    for key, spec in props.items():
        if key not in body:
            continue
        val = body[key]
        if "const" in spec and val != spec["const"]:
            errors.append(f"{key!r} must be {spec['const']!r}, got {val!r}")
        if "enum" in spec and val not in spec["enum"]:
            errors.append(f"{key!r} must be one of {spec['enum']!r}, got {val!r}")
        if "type" in spec:
            types = spec["type"] if isinstance(spec["type"], list) else [spec["type"]]
            ok = False
            for t in types:
                ok = ok or (
                    (t == "object" and isinstance(val, dict))
                    or (t == "array" and isinstance(val, list))
                    or (t == "string" and isinstance(val, str))
                    or (t == "integer" and isinstance(val, int) and not isinstance(val, bool))
                    or (t == "number" and isinstance(val, (int, float)) and not isinstance(val, bool))
                    or (t == "boolean" and isinstance(val, bool))
                    or (t == "null" and val is None)
                )
            if not ok:
                errors.append(f"{key!r} must be of type {spec['type']!r}, got {type(val).__name__}")
        if "minimum" in spec and isinstance(val, (int, float)) and val < spec["minimum"]:
            errors.append(f"{key!r} must be >= {spec['minimum']}, got {val}")
    return errors
