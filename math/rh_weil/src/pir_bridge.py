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

_PROGRAM = Path(__file__).resolve().parents[1]
if str(_PROGRAM) not in sys.path:
    sys.path.insert(0, str(_PROGRAM))

# ENG-006 §11: the new content kinds. A consumer that requires PSD must never be
# satisfied by an inertia certificate, so the fact is *tagged* with what it
# actually establishes and the decision is delegated to the one predicate that
# owns it rather than re-implemented from the fact's fields.
from inertia.certificate import (  # noqa: E402
    KIND_INERTIA,
    KIND_STRATIFICATION,
    satisfies_psd_requirement,
)

# ENG-007 §12: the formal channel. A formal theorem strengthens an exact
# theorem dependency; it never converts interval numerical evidence to FORMAL.
# The two warrants therefore travel in two separate fields on every fact, and
# the formal artifact itself carries no numeric warrant at all.
import formal_evidence  # noqa: E402
from formal_evidence import KIND_FORMAL  # noqa: E402

# ENG-008 §WO-RH-54: the kinds are declared in one place, with an explicit
# answer to "may this ever satisfy a PSD consumer?" attached to each. This list
# used to be maintained here and re-asserted as a frozen literal in
# scripts/ci_inertia.py; the two drifted the moment ENG-007 added a kind, and
# the gate went red on merge. Deriving both from the registry means adding a
# kind is one edit, and the gate proves the licensing decision was made.
from content_kinds import CONTENT_KINDS, REGISTRY, unregistered  # noqa: E402,F401

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
        # ENG-006: the inertia / moment channels. The degree-3 result is emitted
        # under whichever filename the certification actually reached, so both
        # are listed and the absent one is simply skipped.
        ("e1_degree3_odd_positivity_log3_log4.json", "E1", "SOUND",
         "odd degree-3 block: certified inertia over [log3, log4]"),
        ("e1_degree3_odd_inertia_log3_log4.json", "E1", "SOUND",
         "odd degree-3 block: certified inertia stratification over [log3, log4]"),
        ("e1_degree3_odd_moments_log3_log4.json", "E1", "SOUND",
         "odd degree-3 spectral moments m1..m4 and rank-trace, via the Atlas B1 adapter"),
        ("e3_degree3_odd_scan_log3_log4.json", "E3", "HEURISTIC",
         "fresh Candidate-A degree-3 scan — E3 evidence, never a warrant"),
        ("external/connes_cvs_crossvalidation_v0.1.json", "E3", "HEURISTIC", "external cross-check only"),
        # ENG-007: the formal boundary. Emitted at E0 because every theorem it
        # reports is an exact finite algebraic fact, machine-checked under a
        # pinned toolchain -- PIR's ladder has no FORMAL rung, so the FORMAL
        # warrant travels in the content as logical_implication_warrant rather
        # than being smuggled into the evidence level.
        ("formal_theorem_certificate.json", "E0", "SOUND",
         "machine-checked finite theorems only — proves implications, never a numeric bound"),
        # ENG-008: the 3x3 even block -- the first block where inertia,
        # moments and conditioning-by-congruence are exercised beyond what a
        # 2x2 determinant already carries.
        ("e0_degree4_even3_exact_identities.json", "E0", "SOUND", ""),
        ("e1_degree4_even3_inertia_log3_log4.json", "E1", "SOUND",
         "3x3 even block: certified whole-cell inertia over [log3, log4]"),
        ("e1_degree4_even3_positivity_log3_log4.json", "E1", "SOUND",
         "3x3 even block: uniform positive definiteness over [log3, log4]"),
        ("e1_degree4_even3_moments_log3_log4.json", "E1", "SOUND",
         "3x3 even block: spectral moments m1..m4 and rank-trace at sample points"),
        ("e3_degree4_even3_crosscheck.json", "E3", "HEURISTIC",
         "independent mpmath/SymPy assembly — E3 regression evidence, never a warrant"),
        ("e3_pilot3_even_conditioning_log3_log4.json", "E3", "HEURISTIC",
         "ENG-008 preparation: mpmath preview of the 3x3 even block, floating "
         "eigenvalue solver — E3 evidence, never a warrant"),
        # ENG-009: the pencil. The reference metric is E0 (exact rational
        # Sylvester + congruence); the gap enclosures are E1; the structural
        # dataset quotes only certified numbers; the scaling models and the
        # next-block artifacts are plans and say so.
        ("e0_eng009_reference_metric.json", "E0", "SOUND",
         "exact L2 reference metric: PD for every L > 0, exact rational proof"),
        ("e1_eng009_generalized_gap_log3_log4.json", "E1", "SOUND",
         "generalized gap lambda_min(G, M): uniform lower bounds by shifted "
         "positivity, Rayleigh upper bounds at the bottleneck"),
        ("eng009_structural_dataset.json", "E1", "SOUND",
         "cross-block structural dataset regenerated from promoted certificates"),
        ("e3_eng009_scaling_models.json", "E3", "HEURISTIC",
         "exploratory finite scaling models — E3, never promotable to an "
         "infinite claim"),
        ("eng009_next_block_selection.json", "E3", "HEURISTIC",
         "the ENG-010 target selection — a plan, not a numeric fact"),
        ("e3_eng010_even4_preview.json", "E3", "HEURISTIC",
         "ENG-010 preparation: float preview of the 4x4 even block — E3 "
         "evidence, never a warrant"),
        # ENG-010: the 4x4 block itself and the model adjudication.
        ("e0_degree6_even4_exact_identities.json", "E0", "SOUND", ""),
        ("e0_eng010_even4_reference_metric.json", "E0", "SOUND",
         "exact L2 reference metric on {1, b, b^2, b^3}: PD for every L > 0"),
        ("e3_degree6_even4_crosscheck.json", "E3", "HEURISTIC",
         "independent mpmath/SymPy assembly — E3 regression evidence, never a warrant"),
        ("e1_degree6_even4_inertia_log3_log4.json", "E1", "SOUND",
         "4x4 even block: certified whole-cell inertia over [log3, log4]"),
        ("e1_degree6_even4_positivity_log3_log4.json", "E1", "SOUND",
         "4x4 even block: uniform positive definiteness over [log3, log4]"),
        ("e1_eng010_even4_generalized_gap_log3_log4.json", "E1", "SOUND",
         "4x4 even block: generalized gap enclosure against the exact L2 metric"),
        ("e1_degree6_even4_moments_log3_log4.json", "E1", "SOUND",
         "4x4 even block: spectral moments and rank-trace at sample points"),
        ("eng010_scaling_model_adjudication.json", "E1", "SOUND",
         "the preregistered ENG-009 models adjudicated against the certified "
         "n=4 gap, before any refit"),
        ("eng010_information_comparison_report.json", "E1", "SOUND",
         "information-channel comparison at n = 4, certified numbers only"),
        ("e3_eng010_scaling_models_refit.json", "E3", "HEURISTIC",
         "post-adjudication exploratory refit — E3, never promotable"),
        ("eng011_target_selection.json", "E3", "HEURISTIC",
         "the ENG-011 selection — a plan, not a numeric fact"),
        # ENG-011: the 5x5 block, boundary/Schur analyses and adjudication.
        ("e0_degree8_even5_exact_identities.json", "E0", "SOUND", ""),
        ("e0_eng011_even5_reference_metric.json", "E0", "SOUND",
         "exact L2 reference metric on {1, b, b^2, b^3, b^4}: PD for L > 0"),
        ("e3_degree8_even5_crosscheck.json", "E3", "HEURISTIC",
         "independent mpmath/SymPy assembly — E3 regression evidence, never a warrant"),
        ("e1_degree8_even5_inertia_log3_log4.json", "E1", "SOUND",
         "5x5 even block: certified whole-cell inertia over [log3, log4]"),
        ("e1_degree8_even5_positivity_log3_log4.json", "E1", "SOUND",
         "5x5 even block: uniform positive definiteness over [log3, log4]"),
        ("e1_eng011_even5_generalized_gap_log3_log4.json", "E1", "SOUND",
         "5x5 even block: generalized gap enclosure with interior bottleneck"),
        ("e1_degree8_even5_moments_log3_log4.json", "E1", "SOUND",
         "5x5 even block: spectral moments and rank-trace at sample points"),
        ("eng011_boundary_bottleneck_analysis.json", "E1", "SOUND",
         "boundary bottleneck decomposition with certified derivative signs"),
        ("eng011_even5_schur_analysis.json", "E1", "SOUND",
         "Schur complement of the b^4 direction via certified solves"),
        ("eng011_scaling_model_adjudication.json", "E1", "SOUND",
         "the frozen ENG-010 refits adjudicated against the certified n=5 gap"),
        ("eng011_information_comparison_report.json", "E1", "SOUND",
         "information-channel comparison at n = 5, certified numbers only"),
        ("eng012_target_selection.json", "E3", "HEURISTIC",
         "the ENG-012 selection — a plan, not a numeric fact"),
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
        kind = cert.get("content_kind")
        if kind:
            content["content_kind"] = kind
        # §11: what this fact does and does not license. ``satisfies_psd`` is
        # answered by the inertia module's predicate for every certificate, so a
        # PSD-requiring consumer reads one field instead of inferring positivity
        # from a signature it may not understand.
        content["satisfies_psd_requirement"] = bool(satisfies_psd_requirement(cert))
        # §12: which warrant grades the numbers, and which grades the step from
        # those numbers to the claim. They are different questions and a
        # consumer that conflates them would read a proved implication as a
        # certified value.
        if kind == KIND_FORMAL:
            content["numeric_warrant"] = None
            content["logical_implication_warrant"] = cert.get(
                "logical_implication_warrant", formal_evidence.FORMAL_WARRANT
            )
            content["formal_manifest_id"] = cert.get("formal_manifest_id")
            content["formal_theorem_ids"] = list(cert.get("formal_theorem_ids") or ())
            content["lean_toolchain"] = cert.get("lean_toolchain")
            content["mathlib_commit"] = cert.get("mathlib_commit")
            content["axioms"] = list(cert.get("axioms") or ())
            content["unproved_statements"] = [
                {"id": u.get("id"), "status": u.get("status"), "warrant": u.get("warrant")}
                for u in cert.get("unproved_statements") or ()
            ]
        else:
            formal = formal_evidence.formal_block(
                fname, content["evidence_class_declared"]
            )
            if formal:
                content.update(formal)
        if kind in (KIND_INERTIA, KIND_STRATIFICATION):
            content["inertia"] = {
                "n_positive": cert.get("n_positive"),
                "n_negative": cert.get("n_negative"),
                "n_zero": cert.get("n_zero"),
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
            "content_kinds": list(CONTENT_KINDS),
            "formal_manifest_id": formal_evidence.manifest_id(),
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
