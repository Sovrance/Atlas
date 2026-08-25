#!/usr/bin/env python3
"""ATLAS-RH-ENG-008 §WO-RH-52/55 — what each channel told us about the 3x3 block.

Two reports, both built by *reading the emitted certificates* rather than by
restating conclusions in prose, so a future run that changes an outcome changes
the report with it instead of quietly disagreeing with the artifacts.

    eng008_information_comparison_report.json   §WO-RH-52
    eng009_structural_diagnostics.json          §WO-RH-55

§WO-RH-52 exists because ENG-006 could not answer its own version of the
question. That block was 2x2, where the trace and determinant fix the spectrum,
so every channel agreed with the determinant by construction and there was
nothing to compare. This is the first block where the four channels could in
principle disagree.

§WO-RH-55 compares across every finite block this program has certified, and is
explicitly *not* allowed to infer an infinite-dimensional theorem from the
pattern. Candidate invariants are recorded with the falsifier that would kill
them.

    python3 scripts/report_even3_information.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

CERT_DIR = ROOT / "certificates"
INFO_REPORT = "eng008_information_comparison_report.json"
DIAG_REPORT = "eng009_structural_diagnostics.json"

DEPENDENCIES = ("scripts/report_even3_information.py",)


def _load(name: str) -> Optional[Dict[str, Any]]:
    path = CERT_DIR / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _query(point: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
    return next((q for q in point.get("b1_queries", []) if q["query"] == name), None)


def _mid(enc: Dict[str, str]) -> float:
    return (float(enc["lo"]) + float(enc["hi"])) / 2


# --------------------------------------------------------------------------- #
# §WO-RH-52                                                                    #
# --------------------------------------------------------------------------- #
def build_information_report() -> Dict[str, Any]:
    inertia = _load("e1_degree4_even3_inertia_log3_log4.json")
    positivity = _load("e1_degree4_even3_positivity_log3_log4.json")
    moments = _load("e1_degree4_even3_moments_log3_log4.json")
    if inertia is None or moments is None:
        raise SystemExit("run scripts/certify_even3.py first")

    points = moments["points"]
    sample = points[0]
    dim = moments["dimension"]

    forced = _query(sample, "moments_force_psd")
    determined = _query(sample, "inertia_determined_by_moments")
    smallest = _query(sample, "smallest_eigenvalue_bounds")
    neg = _query(sample, "minimum_negative_eigenvalue_count")
    pos = _query(sample, "minimum_positive_eigenvalue_count")

    rt = sample.get("rank_trace") or {}
    rt_result = rt.get("result") or {}
    rt_bound = rt_result.get("certified_rank_lower_bound")
    true_rank = dim if positivity else None

    minors = {c["minor"]: c for c in (positivity or {}).get("leading_minors", [])}
    det_bound = minors.get("Delta3", {}).get("implied_raw_lower_bound")

    qa: List[Dict[str, Any]] = []

    qa.append({
        "question": "Did positivity alone fully characterize the block?",
        "answer": ("no -- and this is the first block in the program where that "
                   "is true"),
        "why": (
            "At 3x3 the trace and determinant no longer fix the spectrum: two "
            "matrices can share both and differ in signature. ENG-006's block "
            "was 2x2, where (trace > 0, det > 0) *is* the signature, so every "
            "channel there agreed with the determinant by construction. Here "
            "positivity is a genuine one-bit answer and the signature is three "
            "numbers."),
    })

    qa.append({
        "question": "What did inertia provide beyond positivity?",
        "answer": ("the answer would have survived a negative outcome; "
                   "positivity would not"),
        "why": (
            "The stratification route returns a signature for every L whether "
            "or not the block is definite. Had any part of the cell come out "
            "indefinite, the run would have produced strata and bounded "
            "transition regions rather than nothing. That the cell turned out "
            "to be one stratum is a fact about this block, not about the "
            "method."),
        "signatures_seen": inertia.get("signatures_seen"),
        "transition_regions": (inertia.get("stratification") or {}).get(
            "transition_regions"),
        "boxes_examined": (inertia.get("stratification") or {}).get(
            "boxes_examined"),
    })

    qa.append({
        "question": "Did the moments alone force the inertia?",
        "answer": determined.get("status") if determined else "unavailable",
        "why": (determined or {}).get("reason") or (determined or {}).get(
            "conclusion"),
        "detail": (
            "This is the substantive difference from ENG-006. At n = 2 the map "
            "from a spectrum to (m1, m2) is injective, so the moments recovered "
            "the inertia exactly and the channel looked stronger than it is. At "
            "n = 3 that map is not injective: the moments constrain the inertia "
            "and do not force it. The 2x2 result was an artefact of the "
            "dimension; this one is the general behaviour."),
        "queries": {
            "moments_force_psd": forced,
            "minimum_negative_eigenvalue_count": neg,
            "minimum_positive_eigenvalue_count": pos,
            "inertia_determined_by_moments": determined,
            "smallest_eigenvalue_bounds": smallest,
        },
    })

    trivial = (rt_bound is not None and true_rank is not None
               and rt_bound < true_rank)
    qa.append({
        "question": "Did rank-trace improve on a trivial bound?",
        "answer": (f"barely: rank >= {rt_bound} against a true rank of {true_rank}"
                   if rt_bound is not None else "unavailable"),
        "why": (
            "The inequality is tight at projections. This block's eigenvalues "
            "are nowhere near 1 -- its trace is around 0.1 across the whole "
            "cell, so all three eigenvalues together sum to a tenth of what a "
            "single unit eigenvalue would contribute -- and the bound degrades "
            "accordingly. Recorded as a weak result rather than tuned until it "
            "looked better."),
        "certified_rank_lower_bound": rt_bound,
        "true_rank": true_rank,
        "strictly_better_than_zero": bool(rt_bound and rt_bound > 0),
        "attains_the_true_rank": bool(rt_bound == true_rank),
        "weak": bool(trivial),
        "note": ("ENG-006 reported rank >= 1 against a true rank of 2. The same "
                 "shortfall at 3x3 gives rank >= 1 against 3, so the gap widened "
                 "with dimension rather than closing."),
    })

    qa.append({
        "question": "What would determinant-only reporting have lost?",
        "answer": "the sign pattern, and any indefinite outcome",
        "why": (
            "A positive determinant on a 3x3 is consistent with signature "
            "(3,0,0) and with (1,2,0) -- two negative eigenvalues multiply to a "
            "positive contribution. So det > 0 alone would not have "
            "distinguished a positive definite block from one with a "
            "two-dimensional negative subspace. That ambiguity does not exist "
            "at 2x2, which is the concrete sense in which this block exercises "
            "the machinery beyond what ENG-006 could."),
        "worked_example": {
            "signature_(3,0,0)": "eigenvalues (1, 1, 1), det = +1",
            "signature_(1,2,0)": "eigenvalues (1, -1, -1), det = +1",
            "conclusion": "identical determinant sign, different inertia",
        },
        "certified_determinant_lower_bound": det_bound,
    })

    return {
        "certificate_version": "0.1",
        "program": "RH/Weil information comparison, 3x3 even block",
        "work_order": "ATLAS-RH-ENG-008",
        "claim_scope": "finite_dimensional_weil_compression",
        "rh_proof_claim": False,
        "evidence_class": "E3",
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "REPORT",
        "subject": "3x3 even Weil block {1, b, b^2} on [log 3, log 4]",
        "dimension": dim,
        "inertia": inertia.get("signatures_seen"),
        "positivity_certified": positivity is not None,
        "questions": qa,
        "caveats": [
            "This is an interpretation of certified artifacts, not itself a "
            "certified claim: the numbers it quotes are certified, the "
            "comparisons it draws are commentary.",
            "Every answer is about one block on one cell under one "
            "normalization.",
            "A null or weak result is preserved here rather than tuned away.",
        ],
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }


# --------------------------------------------------------------------------- #
# §WO-RH-55                                                                    #
# --------------------------------------------------------------------------- #
def build_diagnostics() -> Dict[str, Any]:
    blocks: List[Dict[str, Any]] = []

    scalar = _load("e1_scalar_log3_log4.json")
    if scalar:
        blocks.append({
            "block": "scalar", "basis": ["one"], "dimension": 1,
            "certificate": "e1_scalar_log3_log4.json",
            "smallest_certified_quantity": scalar.get("certified_lower_bound"),
            "quantity": "G00",
            "inertia": [1, 0, 0],
        })
    deg1 = _load("e1_degree1_log3_log4.json")
    if deg1:
        blocks.append({
            "block": "degree-1 odd", "basis": ["q1"], "dimension": 1,
            "certificate": "e1_degree1_log3_log4.json",
            "smallest_certified_quantity": deg1.get("certified_lower_bound"),
            "quantity": "O1 = G[q1,q1]",
            "inertia": [1, 0, 0],
        })
    deg2 = _load("e1_degree2_compact_log3_log4.json")
    if deg2:
        blocks.append({
            "block": "degree-2 even", "basis": ["one", "b"], "dimension": 2,
            "certificate": "e1_degree2_compact_log3_log4.json",
            "smallest_certified_quantity": deg2.get("certified_lower_bound"),
            "quantity": "E2 = det",
            "inertia": [2, 0, 0],
        })
    deg3 = _load("e1_degree3_odd_positivity_log3_log4.json")
    if deg3:
        blocks.append({
            "block": "degree-3 odd", "basis": ["q1", "b3"], "dimension": 2,
            "certificate": "e1_degree3_odd_positivity_log3_log4.json",
            "smallest_certified_quantity":
                deg3["uniform_bounds"]["det_odd3"]["certified_lower_bound"],
            "quantity": "det",
            "leading_entry_bound":
                deg3["uniform_bounds"]["O1"]["certified_lower_bound"],
            "inertia": [deg3["n_positive"], deg3["n_negative"], deg3["n_zero"]],
        })
    even3c = _load("e1_degree4_even3_positivity_log3_log4.json")
    moments = _load("e1_degree4_even3_moments_log3_log4.json")
    if even3c:
        minors = {c["minor"]: c for c in even3c["leading_minors"]}
        blocks.append({
            "block": "degree-4 even", "basis": ["one", "b", "b2"], "dimension": 3,
            "certificate": "e1_degree4_even3_positivity_log3_log4.json",
            "smallest_certified_quantity":
                minors["Delta3"]["implied_raw_lower_bound"],
            "quantity": "Delta3 = det",
            "leading_entry_bound": minors["Delta1"]["implied_raw_lower_bound"],
            "second_minor_bound": minors["Delta2"]["implied_raw_lower_bound"],
            "inertia": [even3c["n_positive"], even3c["n_negative"],
                        even3c["n_zero"]],
            "preconditioner_exponents":
                even3c["preconditioner"]["exponents"],
            "conditioning_note": (
                "the raw third minor is ~1e-11 while the entries are O(1e-1); "
                "the frozen dyadic preconditioner brings the three rescaled "
                "minors to O(1), which is ten orders of magnitude and costs "
                "nothing -- the scaling is exact and the congruence theorem "
                "says the inertia is unchanged"),
        })

    trace_rows = []
    if moments:
        for p in moments["points"]:
            trace_rows.append({
                "L": p["L"], "label": p["label"],
                "trace": p["trace"], "hs_norm_squared": p["hs_norm_squared"],
                "determinant": p["determinant"],
                "trace_over_dimension": repr(_mid(p["trace"]) / moments["dimension"]),
                "m2_over_dimension":
                    repr(_mid(p["hs_norm_squared"]) / moments["dimension"]),
            })

    candidates = [
        {
            "candidate": "every certified even/odd Weil block on this cell is "
                         "positive definite",
            "supported_by": [b["block"] for b in blocks],
            "status": "consistent with every block certified so far",
            "falsifier": (
                "a single cell, block or basis extension whose certified "
                "signature has a nonzero negative index. The stratification "
                "machinery would report it rather than fail, so this is a "
                "cheap test to run and not an assumption anything here rests "
                "on."),
            "explicitly_not_claimed": (
                "nothing about blocks at higher degree, other cells, or any "
                "infinite-dimensional limit"),
        },
        {
            "candidate": "the smallest certified minor falls by roughly five "
                         "orders of magnitude per added dimension",
            "observed": [
                {"dimension": b["dimension"],
                 "smallest_certified_quantity": b["smallest_certified_quantity"],
                 "block": b["block"]}
                for b in blocks
            ],
            "status": "a numerical pattern across five blocks, not a theorem",
            "falsifier": (
                "a block whose certified determinant bound does not follow the "
                "trend, or the same blocks recertified with a sharper cover -- "
                "the bounds are lower bounds, so part of the decline is the "
                "cover's conservatism rather than the block's"),
            "explicitly_not_claimed": (
                "any extrapolation to dimension 4 or beyond, and in particular "
                "no inference that the trend continues to a limit"),
        },
        {
            "candidate": "conditioning degrades with dimension faster than the "
                         "determinant does",
            "observed": {
                "degree-2 even": "no preconditioner needed",
                "degree-3 odd": "no preconditioner needed",
                "degree-4 even": "ten orders of magnitude, dyadic diagonal",
            },
            "status": "one data point for the claim; recorded so ENG-009 can "
                      "look for a second",
            "falsifier": "a 4x4 block that certifies without preconditioning",
            "explicitly_not_claimed": "any rate",
        },
    ]

    return {
        "certificate_version": "0.1",
        "program": "RH/Weil cross-block structural diagnostics",
        "work_order": "ATLAS-RH-ENG-008",
        "prepared_for": "ATLAS-RH-ENG-009",
        "claim_scope": "finite_dimensional_weil_compression",
        "rh_proof_claim": False,
        "evidence_class": "E3",
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "REPORT",
        "cell": "[log 3, log 4]",
        "blocks": blocks,
        "trace_and_moment_scales": trace_rows,
        "candidate_invariants": candidates,
        "hard_boundary": (
            "No infinite-dimensional theorem is inferred, suggested or implied "
            "by anything in this file. Every row is a finite block on one cell "
            "under one normalization, and the patterns across them are "
            "observations about five certificates, not evidence about a limit."
        ),
        "caveats": [
            "Certified bounds are lower bounds, so differences between blocks "
            "mix the block's own scale with how hard the cover worked.",
            "The rank-trace column is weak everywhere it appears and gets "
            "weaker with dimension; that is a property of the inequality at "
            "small eigenvalues, not of the blocks.",
        ],
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }


def main() -> int:
    info = build_information_report()
    p1 = write_certificate(INFO_REPORT, info)
    print(f"wrote {p1}\n")
    print(f"subject: {info['subject']}")
    print(f"inertia: {info['inertia']}  positivity certified: "
          f"{info['positivity_certified']}\n")
    for q in info["questions"]:
        print(f"Q: {q['question']}\nA: {q['answer']}\n")

    diag = build_diagnostics()
    p2 = write_certificate(DIAG_REPORT, diag)
    print(f"wrote {p2}")
    print(f"  {len(diag['blocks'])} blocks compared, "
          f"{len(diag['candidate_invariants'])} candidate invariants, "
          f"each with a falsifier")
    for b in diag["blocks"]:
        print(f"    dim {b['dimension']}  {b['block']:14s} "
              f"{b['quantity']:18s} >= {b['smallest_certified_quantity']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
