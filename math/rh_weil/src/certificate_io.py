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
        "eng_spec": ("ATLAS-RH-ENG-008 first higher-dimensional certified Weil "
                     "block: the 3x3 even block {1, b, b^2} (baseline: ENG-007)"),
        "current_work_order": "ATLAS-RH-ENG-008",
        "latest_completed_work_order": "ATLAS-RH-ENG-007",
        "orders": {
            "WO-RH-01": "done",
            "WO-RH-02": "done_E0_scalar_cell",
            "WO-RH-03": "done",
            "WO-RH-04": "done_algebraic",
            "WO-RH-05": "recovered_ENG-005_cutoff_free_uniform_E1",
            "WO-RH-06": "done_partial_E0_certs_no_imported_promotion",
            "WO-RH-07": "done_dedicated_runner",
            "WO-RH-08": "done_ENG-006_odd_degree3_implemented_and_certified",
            "WO-RH-09": "recovered_ENG-004_scalar_canary_PROMOTED",
            "WO-RH-10": "recovered_ENG-005_assembly_real_space_archimedean",
            "WO-RH-11": "recovered_ENG-005_candidate_a_pole_wired",
            "WO-RH-12": "recovered_ENG-005_true_weil_gram_candidate_a",
            "WO-RH-13": "recovered_ENG-005_E1_T84_points",
            "WO-RH-14": "recovered_ENG-005_exact_support_length_jets",
            "WO-RH-15": "recovered_ENG-005_E1_T84_uniform_plus_interior_minimum",
            "WO-RH-16": "done_pir_export_partial",
            "WO-RH-17": "done_normalization_adjudicated",
            "WO-RH-18": "done_three_way_internal_crosscheck",
            "ENG-004-P1": "done_candidate_a_centralised_scalar_canary_promoted",
            "ENG-004-P2": "done_ENG-005_degree1_degree2_T84_recovered",
            "ENG-005": "done_core_E1_recovery_and_T84_reconstruction",
            "WO-RH-28": "done_generic_interval_hermitian_inertia_engine",
            "WO-RH-29": "done_exact_congruence_sylvester_regression",
            "WO-RH-30": "done_rank_trace_hilbert_schmidt_engine",
            "WO-RH-31": "done_spectral_moments_and_b1_adapter",
            "WO-RH-32": "done_odd_degree3_exact_block",
            "WO-RH-33": "done_fresh_degree3_E3_scan",
            "WO-RH-34": "done_degree3_E1_certificate",
            "WO-RH-35": "done_pir_runner_ci_integration",
            "WO-RH-36": "done_positivity_vs_inertia_moment_report",
            "ENG-006": "done_inertia_ranktrace_moments_and_degree3_pilot",
            "WO-RH-37": "done_ENG-007_documentation_truth_audit",
            "WO-RH-38": "done_ENG-007_pinned_lean_project_and_theorem_boundary",
            "WO-RH-39": "done_ENG-007_congruence_inertia_and_2x2_3x3_criteria",
            "WO-RH-40": "partial_ENG-007_rank_trace_zero_Q_proved_general_case_recorded_unproved",
            "WO-RH-41": "done_ENG-007_weil_parity_and_determinant_identities",
            "WO-RH-42": "done_ENG-007_certificate_semantics_theorems",
            "WO-RH-43": "done_ENG-007_statement_comparator_axiom_audit_and_manifest",
            "WO-RH-44": "done_ENG-007_formal_evidence_in_pir",
            "WO-RH-45": "done_ENG-007_readme_refresh_and_docs_gate",
            "WO-RH-46": "done_ENG-007_3x3_even_pilot_prepared_E0_and_E3_only",
            "ENG-007": "done_formal_theorem_boundary_and_documentation_gate",
            "WO-RH-47": "done_ENG-008_basis_frozen_and_dyadic_preconditioner_certified",
            "WO-RH-48": "done_ENG-008_six_exact_entries_and_independent_crosscheck",
            "WO-RH-49": "done_ENG-008_derivative_provider_generalized_prior_results_identical",
            "WO-RH-50": "done_ENG-008_rigorous_arb_3x3_assembly",
            "WO-RH-51": "done_ENG-008_uniform_E1_inertia_3_0_0_and_positivity",
            "WO-RH-52": "done_ENG-008_moments_ranktrace_information_comparison",
            "WO-RH-53": "done_ENG-008_lean_3x3_certificate_replay",
            "WO-RH-54": "done_ENG-008_pir_runner_ci_and_docs",
            "WO-RH-55": "done_ENG-008_cross_block_diagnostics_for_ENG-009",
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
            "WO-RH-05/10..15 RECOVERED by ENG-005: regenerated from scratch under "
            "Candidate A and explicitly released; the pre-quarantine values above are "
            "retained as the contrary evidence WO-RH-17 forbids deleting",
            "WO-RH-09 (scalar cell) RECOVERED by ENG-004: regenerated under Candidate A "
            "with a uniform rigorous Arb lower bound and explicitly released",
            "ENG-004: src/pole.py is the single pole implementation; the rejected "
            "(sqrt(3)/2) block is archival in src/rejected_pole.py and production may not import it",
            "ENG-004: G00'' = 4cosh(L/2) - e^{L/2}/sinh(L) equals the E0 curvature "
            "2(r^3-r-1)/(sqrt(r)(r^2-1)); Candidate B cannot reproduce it",
            "Cross-check is three-way INTERNAL; Connes/CvS reports NOT_COMPARABLE and never certifies",
            "ENG-005: the T=84 topology was rescanned fresh under Candidate A. The "
            "earlier 'monotone E2' > 0' reading was a Candidate-B artifact and is NOT "
            "reused: E2' changes sign once, at L* ~ 1.10595, and the minimum is interior",
            "ENG-005: two independent warrants at T=84 -- an exhaustive interval cover "
            "assuming no topology, and an interior-minimum argument locating L*; the "
            "headline bound is the sharper, floored at what the cover alone proves",
            "E1 filename prefixes may hold non-E1 status until gates close",
            "ENG-006: WO-RH-08 (degree 3) is no longer blocked -- the odd degree-3 "
            "block is implemented, its kernels are E0-verified, and it carries a "
            "rigorous inertia result",
            "ENG-006: inertia, rank-trace and spectral-moment certificates are "
            "distinct content kinds from positivity; an inertia certificate never "
            "satisfies a consumer that requires PSD",
            "ENG-006: the prime-shift blocks of the odd degree-3 pair are each "
            "indefinite on the cell, so termwise PSD domination is unavailable and "
            "the assembled entry must be bounded as a whole",
            "ENG-007: the Lean project under formal/ pins an exact toolchain and "
            "mathlib commit; ten finite theorems are proved with no sorry and only "
            "the three standard axioms",
            "ENG-007: a formal theorem strengthens an exact theorem dependency and "
            "never converts interval numerical evidence to FORMAL; PIR facts carry "
            "numeric_warrant and logical_implication_warrant as separate fields",
            "ENG-007: WO-RH-40 is deliberately partial. The Q = 0 case the runtime "
            "uses is proved; the general case is recorded as "
            "EXTERNAL_THEOREM_PENDING_FORMAL_PROOF with a null warrant and no "
            "inhabitant anywhere in the project",
            "ENG-007: zeta-23-lean is an architecture reference pinned at "
            "cec57f9, REFERENCE_ONLY -- nothing is vendored, imported or depended on, "
            "and its toolchain does not currently compose with Atlas's",
            "ENG-007: the 3x3 even pilot block {1, b, b2} is prepared, not certified "
            "-- E0 kernel identities and an E3 conditioning preview only, per the "
            "work order's instruction not to promote a new E1 degree result",
            "ENG-008: the 3x3 even block G[{1, b, b2}] is certified positive "
            "definite with inertia (3,0,0) uniformly on [log 3, log 4], by two "
            "independent routes -- interval LDL* congruence stratified over the "
            "cell, and Sylvester's criterion as three separate adaptive covers",
            "ENG-008: the preconditioner is a diagonal matrix of exact powers of "
            "two, frozen for the cell. Scaling an Arb ball by a power of two is "
            "exact, so the congruence adds no width; invertibility is exact; and "
            "the inertia is unchanged by the ENG-007 congruence theorems",
            "ENG-008: the overlap kernels and the L-derivative machinery are now "
            "derived from the basis coefficients rather than hand-tabulated. The "
            "prior certified bounds were required to come back unchanged and did",
            "ENG-008: at n = 3 the moments no longer force the inertia. The "
            "ENG-006 finding that they did was an artefact of n = 2, where the "
            "map from a spectrum to (m1, m2) is injective",
            "ENG-008: rank-trace gives rank >= 1 against a true rank of 3 -- "
            "weaker than at 2x2, and preserved as a weak result",
        ],
    }
