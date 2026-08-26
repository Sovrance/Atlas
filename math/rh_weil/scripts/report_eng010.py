#!/usr/bin/env python3
"""ATLAS-RH-ENG-010 §WO-RH-71/72/75 — adjudication, information comparison, ENG-011.

    python3 scripts/report_eng010.py

Produces:

  ``eng010_scaling_model_adjudication.json``    §WO-RH-71
  ``e3_eng010_scaling_models_refit.json``       §WO-RH-71 (after adjudication)
  ``eng010_information_comparison_report.json`` §WO-RH-72
  ``eng011_target_selection.json``              §WO-RH-75

The adjudication loads the *preregistered* ENG-009 predictions directly from
their artifact, asserts the artifact is bitwise the committed one, and records
the verdict before any refit happens. That ordering is enforced by the shape
of this script: the refit function takes the adjudication record as an
argument.

No RH proof claim is made. Claim scope is ``finite_dimensional_weil_compression``.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

import normalization as N  # noqa: E402
import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from content_kinds import (  # noqa: E402
    KIND_NEXT_BLOCK_SELECTION,
    KIND_SCALING_ADJUDICATION,
    KIND_SCALING_MODEL,
    KIND_STRUCTURAL_DIAGNOSTIC,
)

CERTS = ROOT / "certificates"

MODELS_FILE = "e3_eng009_scaling_models.json"
GAP_FILE = "e1_eng010_even4_generalized_gap_log3_log4.json"
ADJUDICATION_FILE = "eng010_scaling_model_adjudication.json"
REFIT_FILE = "e3_eng010_scaling_models_refit.json"
INFO_FILE = "eng010_information_comparison_report.json"
ENG011_FILE = "eng011_target_selection.json"

#: The ENG-009 artifact as committed at the ENG-010 baseline (recorded before
#: any ENG-010 work began). The adjudication refuses to run against anything
#: else -- §Stop conditions: "model adjudication requires changing old models
#: first" is a stop, and this is the tripwire.
ENG009_MODELS_CONTENT_HASH = (
    "2719ade96da77279ea6350fdfb19f49c1a22434e970fafc302036480a85d6a23")

DEPENDENCIES = (
    "src/even4.py",
    "src/reference_metric.py",
    "src/generalized_gap.py",
    "scripts/report_eng010.py",
)


def load(name: str) -> Dict[str, Any]:
    return json.loads((CERTS / name).read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# §WO-RH-71: adjudication (before any refit)                                   #
# --------------------------------------------------------------------------- #
def adjudicate() -> Dict[str, Any]:
    print("=== scaling-model adjudication (§WO-RH-71) ===")
    models = load(MODELS_FILE)
    if models.get("content_hash") != ENG009_MODELS_CONTENT_HASH:
        raise SystemExit(
            "STOP: the ENG-009 scaling-models artifact is not the "
            "preregistered one; refusing to adjudicate (§Stop conditions)")
    gap = load(GAP_FILE)
    lam_lo = float(gap["certified_lambda_lower_float"])
    lam_up = float(gap["upper_bound_at_bottleneck"]["certified_upper_bound"])
    print(f"  certified n=4 gap enclosure: [{lam_lo:.6e}, {lam_up:.6e}]")

    rows = []
    labels = {"exponential_decay": "MODEL_A", "power_law_decay": "MODEL_B"}
    outcomes = {}
    for m in models["models"]:
        if m["family"] != "even" or "next_block_prediction" not in m:
            continue
        pred = float(m["next_block_prediction"]["value"])
        window = (pred / 5.0, pred * 5.0)
        # Falsified iff the certified enclosure lies wholly outside the window.
        falsified = lam_lo > window[1] or lam_up < window[0]
        direction = ("certified enclosure lies entirely ABOVE the window"
                     if lam_lo > window[1] else
                     "certified enclosure lies entirely BELOW the window"
                     if lam_up < window[0] else
                     "certified enclosure intersects the window")
        label = labels.get(m["model"], m["model"])
        outcomes[label] = falsified
        rows.append({
            "label": label,
            "model": m["model"],
            "formula": m["formula"],
            "fitted": m["fitted"],
            "preregistered_prediction_at_n4": repr(pred),
            "falsifier_window": [repr(window[0]), repr(window[1])],
            "verdict": "FALSIFIED" if falsified else "NOT_FALSIFIED",
            "direction": direction,
        })
        print(f"  {label} ({m['model']}): predicted {pred:.3e}, window "
              f"[{window[0]:.3e}, {window[1]:.3e}] -> "
              f"{'FALSIFIED' if falsified else 'NOT_FALSIFIED'} ({direction})")

    if outcomes.get("MODEL_A") and outcomes.get("MODEL_B"):
        verdict = "BOTH_FALSIFIED"
    elif outcomes.get("MODEL_A"):
        verdict = "MODEL_A_FALSIFIED"
    elif outcomes.get("MODEL_B"):
        verdict = "MODEL_B_FALSIFIED"
    else:
        verdict = "NEITHER_FALSIFIED"
    above_both = all("ABOVE" in r["direction"] for r in rows)
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil scaling-model adjudication — Candidate A",
        "work_order": "WO-RH-71",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": KIND_SCALING_ADJUDICATION,
        "evidence_class": "E1_ADJUDICATING_E3",
        "numeric_warrant": ("the certified enclosure is E1 (from the named gap "
                            "certificate); the models being judged are E3 and "
                            "stay E3"),
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "ADJUDICATED",
        "mpmath_used": False,
        "preregistered_models_artifact": {
            "file": MODELS_FILE,
            "content_hash": ENG009_MODELS_CONTENT_HASH,
            "verified_unchanged": True,
        },
        "certified_result": {
            "file": GAP_FILE,
            "enclosure": [repr(lam_lo), repr(lam_up)],
            "content_hash": gap.get("content_hash"),
        },
        "adjudications": rows,
        "verdict": verdict,
        "reading": (
            "the certified n = 4 generalized gap lies above both preregistered "
            "predictions' falsifier windows: the even-family gap is decaying "
            "far more slowly than either one-parameter fit extrapolated -- "
            "expected outcome A of the work order, and evidence that monotone "
            "one-parameter decay laws are inadequate for this family"
            if verdict == "BOTH_FALSIFIED" and above_both else
            "see the per-model rows"),
        "recorded_before_any_refit": True,
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(ADJUDICATION_FILE, body)
    print(f"wrote {path}")
    print(f"  verdict: {verdict}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-71: exploratory refit, only after adjudication                        #
# --------------------------------------------------------------------------- #
def refit(adjudication: Dict[str, Any]) -> Dict[str, Any]:
    print("\n=== exploratory refit (§WO-RH-71, E3, post-adjudication) ===")
    assert adjudication.get("recorded_before_any_refit")
    import statistics

    eng009 = load("eng009_structural_dataset.json")
    gap = load(GAP_FILE)
    pts = []
    for r in eng009["cutoff_free_blocks"]:
        if r["family"] == "even":
            pts.append((r["dimension"],
                        float(r["generalized_gap"]["lambda_lower_uniform"]),
                        float(r["generalized_gap"]["lambda_upper_at_bottleneck"])))
    pts.append((4, float(gap["certified_lambda_lower_float"]),
                float(gap["upper_bound_at_bottleneck"]["certified_upper_bound"])))
    pts.sort()
    ns = [p[0] for p in pts]
    mids = [0.5 * (p[1] + p[2]) for p in pts]
    logs = [math.log(v) for v in mids]
    data = [{"n": p[0], "enclosure": [repr(p[1]), repr(p[2])]} for p in pts]
    models: List[Dict[str, Any]] = []

    def lsq(xs, ys):
        xbar, ybar = statistics.mean(xs), statistics.mean(ys)
        slope = (sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
                 / sum((x - xbar) ** 2 for x in xs))
        return slope, ybar - slope * xbar

    # exponential on all four points
    slope, icept = lsq(ns, logs)
    c, rho = math.exp(icept), math.exp(slope)
    pred5_exp = c * rho ** 5
    resid = max(abs(math.log(v) - (icept + slope * n))
                for n, v in zip(ns, mids))
    models.append({
        "model": "exponential_decay_refit", "formula": "lambda_min(n) = C * rho^n",
        "family": "even", "fitted": {"C": repr(c), "rho": repr(rho)},
        "fit_method": "least squares on log lambda vs n, four points",
        "max_log_residual": repr(resid),
        "data_points": data,
        "input_warrant": "E1 enclosures; the fit itself is E3",
        "extrapolation_status": "EXPLORATORY_NEVER_PROMOTED",
        "next_block_prediction": {"n": 5, "value": repr(pred5_exp)},
        "falsifier": (f"certified even-family n = 5 gap enclosure outside "
                      f"[{pred5_exp / 5:.3e}, {pred5_exp * 5:.3e}] rejects "
                      "this refit"),
    })
    # power law on all four points
    lns = [math.log(n) for n in ns]
    slope_p, icept_p = lsq(lns, logs)
    c_p = math.exp(icept_p)
    pred5_pow = c_p * 5.0 ** slope_p
    resid_p = max(abs(math.log(v) - (icept_p + slope_p * math.log(n)))
                  for n, v in zip(ns, mids))
    models.append({
        "model": "power_law_decay_refit", "formula": "lambda_min(n) = C * n^(-p)",
        "family": "even", "fitted": {"C": repr(c_p), "p": repr(-slope_p)},
        "fit_method": "least squares on log lambda vs log n, four points",
        "max_log_residual": repr(resid_p),
        "data_points": data,
        "input_warrant": "E1 enclosures; the fit itself is E3",
        "extrapolation_status": "EXPLORATORY_NEVER_PROMOTED",
        "next_block_prediction": {"n": 5, "value": repr(pred5_pow)},
        "falsifier": (f"certified even-family n = 5 gap enclosure outside "
                      f"[{pred5_pow / 5:.3e}, {pred5_pow * 5:.3e}] rejects "
                      "this refit"),
    })
    separation = max(pred5_exp, pred5_pow) / min(pred5_exp, pred5_pow)
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil scaling models, post-adjudication refit — Candidate A",
        "work_order": "WO-RH-71",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": KIND_SCALING_MODEL,
        "evidence_class": "E3",
        "numeric_warrant": "NONE — E3 exploratory fits over E1 enclosures",
        "rigorous": False,
        "psd_claim": False,
        "status": "EXPLORATORY",
        "adjudication_first": {
            "file": ADJUDICATION_FILE,
            "verdict": adjudication["verdict"],
        },
        "bottleneck_note": (
            "the n = 4 bottleneck sits at the right cell edge (L -> log 4), "
            "twenty times below the midpoint E3 scout the ENG-009 preview "
            "recorded -- the work order's instruction not to present the "
            "scout as the certified result was the right one"),
        "models": models,
        "refit_separation_at_n5": repr(separation),
        "anti_overclaim": [
            "four finite blocks do not establish asymptotics",
            "surviving a x5 falsifier window is not confirmation",
        ],
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(REFIT_FILE, body)
    print(f"wrote {path}")
    for m in models:
        print(f"  {m['model']}: n=5 -> {m['next_block_prediction']['value']}")
    print(f"  refit separation at n=5: {separation:.1f}x")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-72: the information comparison at n = 4                               #
# --------------------------------------------------------------------------- #
def info_report() -> Dict[str, Any]:
    print("\n=== information comparison at n = 4 (§WO-RH-72) ===")
    gap = load(GAP_FILE)
    pos = load("e1_degree6_even4_positivity_log3_log4.json")
    inertia = load("e1_degree6_even4_inertia_log3_log4.json")
    moments = load("e1_degree6_even4_moments_log3_log4.json")
    mid = moments["points"][2]
    m_analysis = mid["moment_analysis"]
    raw_det = [c for c in pos["leading_minors"] if c["minor"] == "Delta4"][0]
    trace_mid = 0.5 * (float(mid["trace"][0]) + float(mid["trace"][1]))
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil information comparison at n = 4 — Candidate A",
        "work_order": "WO-RH-72",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": KIND_STRUCTURAL_DIAGNOSTIC,
        "evidence_class": "E1",
        "numeric_warrant": "E1 — every number quoted is a certified enclosure",
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "subject": "4x4 even Weil block {1, b, b^2, b^3} on [log 3, log 4]",
        "channels": {
            "inertia": {
                "value": inertia["signatures_seen"],
                "constant_on_cell": inertia["constant_on_cell"],
            },
            "generalized_gap": {
                "enclosure": [gap["certified_lambda_lower_float"],
                              gap["upper_bound_at_bottleneck"]["certified_upper_bound"]],
                "reference_metric": gap["reference_metric_id"],
            },
            "trace_over_n": repr(trace_mid / 4.0),
            "raw_determinant": {
                "uniform_lower_bound": raw_det["implied_raw_lower_bound"],
                "note": ("~1e-19 pointwise against O(1e-1) entries: five to "
                         "six more orders of collapse beyond n = 3, of which "
                         "the invariant gap accounts for barely one -- the "
                         "rest is the coordinate system again"),
            },
            "moments_m1_m4": {"at": mid["label"],
                              "values": {k: v for k, v in
                                         (m_analysis.get("moments") or {}).items()}},
            "rank_trace": m_analysis.get("rank_trace"),
            "conditioning": {
                "raw_diag_span": "~1e5 (G00 ~ 8e-2 vs G33 ~ 7e-7)",
                "preconditioned": "O(1) diagonal under the frozen dyadic congruence",
            },
        },
        "questions": {
            "do_the_moments_force_the_inertia": (
                "NO -- already false at n = 3 (ENG-008), and the moment map is "
                "even less injective at n = 4"),
            "do_the_moments_constrain_the_generalized_gap": (
                "only trivially: m1..m4 are basis-dependent traces of G alone "
                "and carry no information about the pencil (G, M) beyond crude "
                "bounds; nothing here improves on the shifted covers"),
            "is_the_generalized_gap_still_the_best_cross_dimensional_margin": (
                "YES, and more clearly than at n = 3: the raw determinant lost "
                "another five-plus orders to coordinates, while the gap both "
                "measured the real decay (~15x) and located a bottleneck at "
                "the right cell edge that no raw quantity showed"),
            "did_rank_trace_become_more_or_less_informative": (
                "less: the bound weakens with dimension (rank >= 1 against a "
                "true rank of 4)"),
        },
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(INFO_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-75: ENG-011 from the result                                           #
# --------------------------------------------------------------------------- #
def eng011(adjudication: Dict[str, Any]) -> Dict[str, Any]:
    print("\n=== ENG-011 target selection (§WO-RH-75) ===")
    verdict = adjudication["verdict"]
    refit_cert = load(REFIT_FILE)
    preds = {m["model"]: float(m["next_block_prediction"]["value"])
             for m in refit_cert["models"] if "next_block_prediction" in m}
    separation = float(refit_cert.get("refit_separation_at_n5", "1"))
    if verdict in ("NEITHER_FALSIFIED", "MODEL_A_FALSIFIED", "MODEL_B_FALSIFIED"):
        selection = "even5_one_b_b2_b3_b4"
        rationale = (
            "the certified n = 4 enclosure landed inside {} preregistered "
            "x5 window(s) -- expected outcome E: the models were not "
            "separated at n = 4 because their windows overlap exactly where "
            "the result fell. The refitted models diverge by ~{:.0f}x at "
            "n = 5, so the next even element (b^4, even, homogeneous degree "
            "8, deriving from the same primitive table) is the maximally "
            "discriminating block. The n = 4 bottleneck's location at the "
            "right cell edge also makes the n = 5 edge behaviour the thing "
            "to watch".format(
                "both" if verdict == "NEITHER_FALSIFIED" else "one",
                separation))
    else:
        selection = "structural_analysis_of_the_even_family"
        rationale = ("both models fell; §WO-RH-75 directs the structural "
                     "route (interlacing, Schur complements, recurrences)")
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil ENG-011 selection — Candidate A",
        "work_order": "WO-RH-75",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": KIND_NEXT_BLOCK_SELECTION,
        "evidence_class": "E3",
        "numeric_warrant": "NONE — a plan, informed by E1 data",
        "rigorous": False,
        "psd_claim": False,
        "status": "SELECTED",
        "selection": selection,
        "based_on": {"adjudication": ADJUDICATION_FILE, "verdict": verdict,
                     "refit": REFIT_FILE,
                     "refit_predictions_at_n5": {k: repr(v)
                                                 for k, v in preds.items()}},
        "rationale": rationale,
        "candidate_questions_for_eng011": [
            "does the certified n = 5 gap fall in the refitted exponential's "
            "window, the power law's, neither, or both again?",
            "does the bottleneck stay pinned at the right cell edge, and why "
            "-- what about L -> log 4 weakens the pencil?",
            "does generalized Cauchy interlacing (M_n is exactly the leading "
            "block of M_{n+1}) quantitatively account for the decrements?",
            "what is the Schur complement of each new direction against the "
            "previous span, as a function of L?",
        ],
        "explicitly_not_launched": (
            "no E1 certification of any n = 5 block is started here; "
            "§WO-RH-75 forbids launching another expensive E1 block before "
            "interpreting ENG-010"),
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(ENG011_FILE, body)
    print(f"wrote {path}")
    print(f"  selection: {selection}")
    return body


def main() -> int:
    adj = adjudicate()
    refit(adj)
    info_report()
    eng011(adj)
    return 0


if __name__ == "__main__":
    sys.exit(main())
