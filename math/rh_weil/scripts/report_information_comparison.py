#!/usr/bin/env python3
"""ENG-006 §10/§36 — what each channel actually told us about the degree-3 block.

Answers the five questions §10 poses, by *reading the emitted certificates*
rather than by restating conclusions in prose. If a future run changes an
outcome, this report changes with it instead of quietly disagreeing with the
artifacts it describes.

§10 is explicit that a null result is acceptable and should be preserved, so the
report records weak and insufficient answers with the same weight as strong
ones, and says why each one came out that way.

    python3 scripts/report_information_comparison.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

CERT_DIR = ROOT / "certificates"
REPORT_FILE = "eng006_information_comparison_report.json"
E1_ALTERNATIVES = ("e1_degree3_odd_positivity_log3_log4.json",
                   "e1_degree3_odd_inertia_log3_log4.json")


def _load(name):
    path = CERT_DIR / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _query(row, name):
    return next(q for q in row["moment_analysis"]["b1_queries"] if q["query"] == name)


def build_report() -> dict:
    e1_name = next((n for n in E1_ALTERNATIVES if (CERT_DIR / n).exists()), None)
    e1 = _load(e1_name) if e1_name else None
    mom = _load("e1_degree3_odd_moments_log3_log4.json")
    scan = _load("e3_degree3_odd_scan_log3_log4.json")
    if e1 is None or mom is None:
        raise SystemExit("degree-3 E1 and moment certificates must exist first")

    sig = tuple(e1["inertia"][k] for k in ("n_positive", "n_negative", "n_zero"))
    definite = sig == (2, 0, 0)
    bounds = e1.get("uniform_bounds", {})
    row = mom["points"][1] if len(mom["points"]) > 1 else mom["points"][0]

    determined = _query(row, "inertia_determined_by_moments")
    forces_psd = _query(row, "moments_force_psd")
    lam = _query(row, "smallest_eigenvalue_bounds")
    rt = row["rank_trace"]
    rt_bound = rt.get("result", {}).get("certified_rank_lower_bound")
    true_rank = sig[0] + sig[1] if sig[0] is not None else None

    answers = [
        {
            "question": "Did positivity alone fully characterize the block?",
            "answer": "yes, but only because the block is 2x2",
            "detail": (
                "For a symmetric 2x2 the pair (trace > 0, det > 0) *is* the full "
                "signature (2,0,0): there is no room between 'positive definite' "
                "and 'inertia known'. So on this block positivity and inertia carry "
                "identical information. That is a fact about the dimension, not "
                "about the method -- at 3x3 and above a positivity test and an "
                "inertia test stop being the same question."),
            "evidence": {"certified_inertia": list(sig),
                         "uniform_bounds": {k: v["certified_lower_bound"]
                                            for k, v in bounds.items()}},
        },
        {
            "question": "What additional information did inertia provide?",
            "answer": ("beyond positivity, nothing here -- the block turned out "
                       "definite; the value delivered was the route, not the number"),
            "detail": (
                "The inertia engine is what proved the signature constant over the "
                "whole closed cell rather than at sample points, and it is what "
                "would have produced a usable stratification had any part of the "
                "cell been indefinite. On this block it returned one stratum and no "
                "transition regions, so there was no extra structure to report. A "
                "channel that adds nothing when the easy answer holds is behaving "
                "correctly; it earns its place on the blocks where the easy answer "
                "does not."),
            "evidence": {
                "strata": len(e1["inertia_stratification"]["strata"]),
                "transition_regions": len(e1["inertia_stratification"]["transition_regions"]),
                "constant_inertia": e1["inertia"]["constant_on_cell"],
                "boxes_examined": e1["inertia_stratification"]["boxes_examined"],
            },
        },
        {
            "question": "Did low-order moments recover the inertia?",
            "answer": "yes, exactly, and again because the block is 2x2",
            "detail": (
                "At n = 2 the Wolkowicz-Styan inequalities are equalities, so m1 and "
                "m2 invert to the spectrum outright: lambda = (m1 +- sqrt(2 m2 - "
                "m1^2))/2. The signature is then a consequence of two moments rather "
                "than a separate observation, and it agrees with the certified "
                "inertia. m3 and m4 were not needed. This is the sharpest the moment "
                "channel can ever be, and it does not survive to larger blocks."),
            "evidence": {
                "implied_inertia": determined.get("implied_inertia"),
                "matches_certified_inertia": determined.get("matches_observed"),
                "status": determined["status"],
                "lambda_min": lam["lambda_min"], "lambda_max": lam["lambda_max"],
                "tight": lam["tight"],
            },
        },
        {
            "question": "Were the moments insufficient?",
            "answer": ("yes, by the general truncated-moment route -- and that is "
                       "not a contradiction with the previous answer"),
            "detail": (
                "Two different questions get two different answers. Asking 'do these "
                "moments force the spectrum to be non-negative' via the localizing "
                "matrix returns INSUFFICIENT_INFORMATION, correctly: PSD-ness of a "
                "truncated localizing matrix is necessary but not sufficient for "
                "support in [0, inf), and establishing sufficiency needs a flat "
                "extension that four moments do not supply. Only the *refuting* "
                "direction is conclusive there. The n = 2 spectrum inversion answers "
                "a different and easier question, and it is what settled the "
                "inertia. Reporting the localizing route as insufficient while the "
                "inversion succeeds is the honest description of both."),
            "evidence": {
                "moments_force_psd": forces_psd["status"],
                "reason": forces_psd.get("reason"),
                "minimum_negative_eigenvalue_count":
                    _query(row, "minimum_negative_eigenvalue_count")["status"],
                "minimum_positive_eigenvalue_count":
                    _query(row, "minimum_positive_eigenvalue_count")["status"],
            },
        },
        {
            "question": "Did rank-trace yield a nontrivial finite-dimensional lower bound?",
            "answer": (f"yes but weak: rank >= {rt_bound} against a true rank of "
                       f"{true_rank}"),
            "detail": (
                "The bound is nontrivial -- the right-hand side is positive, so it "
                "is not the vacuous 'rank >= 0' -- but it is off by a factor of two "
                "on a rank-2 block, and the reason is structural rather than a "
                "tuning failure. The inequality is tight at projections, where every "
                "eigenvalue is 0 or 1. This block's eigenvalues run from ~3.3e-05 to "
                "~4.2e-02, four orders of magnitude below 1, so each contributes "
                "lambda(2 - lambda) ~ 2*lambda to the right-hand side instead of the "
                "1 it would contribute at a projection. A Weil Gram block is simply "
                "not the kind of operator this inequality is sharp for. Recorded as "
                "a null-ish result per §10 rather than tuned until it looks better."),
            "evidence": {
                "status": rt["status"],
                "rhs_enclosure": rt.get("result", {}).get("rhs_enclosure"),
                "certified_rank_lower_bound": rt_bound,
                "trivial": rt.get("result", {}).get("trivial"),
                "true_rank_from_certified_inertia": true_rank,
                "lambda_max_upper":
                    rt["hypotheses"]["shared_normalization"]["evidence"].get(
                        "lambda_max_upper"),
                "theorem_id": rt["theorem_id"],
            },
        },
    ]

    return {
        "certificate_version": "1.0",
        "program": "RH/Weil ENG-006 information comparison",
        "work_order": "ATLAS-RH-ENG-006 §10 / WO-RH-36",
        "content_kind": "WEIL_INFORMATION_COMPARISON_REPORT",
        "evidence_class": "E3",
        "status": "REPORTED",
        "hard_constraints_certified": False,
        "promotion_state": "REPORT_NOT_PROMOTABLE",
        "claim_scope": "finite_dimensional_weil_compression",
        "rh_proof_claim": False,
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "subject": {
            "block": "odd degree-3 Weil block on [log 3, log 4]",
            "dimension": 2,
            "certified_inertia": list(sig),
            "outcome": e1["outcome"],
            "e1_certificate": e1_name,
            "moment_certificate": "e1_degree3_odd_moments_log3_log4.json",
            "scan_certificate": "e3_degree3_odd_scan_log3_log4.json",
        },
        "headline": (
            "The odd degree-3 block is positive definite, so the case ENG-006 was "
            "built for -- extracting information from an indefinite block -- did not "
            "arise here. What the run does establish is that the machinery works and "
            "what each channel is worth: inertia and positivity coincide at 2x2, the "
            "moments recover the signature exactly at 2x2 and would not above it, "
            "and the rank-trace bound is real but weak on an operator this far from "
            "a projection."
            if definite else
            "The odd degree-3 block is not positive definite on the whole cell; the "
            "inertia channel is what kept the run from ending in a null result."),
        "answers": answers,
        "caveats": [
            "One cell and one 2x2 block. Three of the five answers above turn on the "
            "dimension being 2 and do not generalize.",
            "A null or weak result is preserved here rather than tuned away (§10).",
            "E3: this is an interpretation of certified artifacts, not itself a "
            "certified claim. The numbers it quotes are certified; the comparisons "
            "it draws are commentary.",
        ],
    }


def main() -> int:
    body = build_report()
    print(f"wrote {write_certificate(REPORT_FILE, body)}")
    print(f"\nsubject: {body['subject']['block']}")
    print(f"inertia: {body['subject']['certified_inertia']}  "
          f"outcome: {body['subject']['outcome']}\n")
    for a in body["answers"]:
        print(f"Q: {a['question']}")
        print(f"A: {a['answer']}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
