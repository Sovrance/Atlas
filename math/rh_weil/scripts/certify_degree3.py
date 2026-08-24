#!/usr/bin/env python3
"""ENG-006 §7/§8/§9/§10 — odd degree-3 block: scan, E1 result, moments.

Three artifacts, in the order the work order requires:

  e3_degree3_odd_scan_log3_log4.json        fresh Candidate-A scan (§8, E3)
  e1_degree3_odd_<positivity|inertia>_...   the rigorous result (§9, E1)
  e1_degree3_odd_moments_log3_log4.json     moments and rank-trace (§10)

The scan comes first because §9 says the E1 strategy is chosen *after* it, not
before. Which of the two E1 filenames is written is decided by what the
certification actually establishes, not by what was hoped for: a positive
definite block yields a positivity certificate, an indefinite one yields an
inertia stratification, and §17 counts both as success.

    python3 scripts/certify_degree3.py [--release] [--quick] [--stage ...]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import degree3 as D3  # noqa: E402
import interval_cover as IC  # noqa: E402
import promotion  # noqa: E402
import weil_entries as WE  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from inertia.certificate import (  # noqa: E402
    build_stratification_certificate,
    satisfies_psd_requirement,
)
from inertia.stratify import certify_inertia_family  # noqa: E402
from interval_backend import (  # noqa: E402
    FlintUnavailable,
    backend_info,
    interval_box,
    require_flint,
    set_precision_bits,
)
from moments.adapter import analyse  # noqa: E402
from ranktrace.theorem import HYPOTHESES, rank_trace_lower_bound  # noqa: E402

CERT_DIR = ROOT / "certificates"
SCAN_FILE = "e3_degree3_odd_scan_log3_log4.json"
POSITIVITY_FILE = "e1_degree3_odd_positivity_log3_log4.json"
INERTIA_FILE = "e1_degree3_odd_inertia_log3_log4.json"
MOMENTS_FILE = "e1_degree3_odd_moments_log3_log4.json"

#: Content kind for an Outcome-A artifact. Deliberately *not* one of §11's four
#: inertia/moment kinds: this one does claim positivity, and is allowed to.
POSITIVITY_KIND = "WEIL_DEGREE3_POSITIVITY_CERTIFICATE"

DEPENDENCIES = (
    "src/pole.py",
    "src/core.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/degree3.py",
    "src/interval_cover.py",
    "src/interval_backend.py",
    "src/normalization.py",
    "inertia/ldl.py",
    "inertia/stratify.py",
    "inertia/certificate.py",
    "moments/spectral_moments.py",
    "moments/adapter.py",
    "ranktrace/theorem.py",
    "scripts/certify_degree3.py",
)

PRECISION_BITS = 160
#: det bottoms out near 1.4e-6 and the interval-L enclosure widens like ~0.15*r,
#: so separation needs r <~ 1e-5. Starting near that spacing means most boxes
#: clear without splitting. This is enough to decide the *signature*, which is
#: all the stratification needs.
INITIAL_CELLS = 8192

#: The magnitude of the bound is a different question from the sign of it. A
#: cover accepts a box the moment its lower end clears the target, so a cover
#: sized just to decide positivity reports whatever the worst box happened to
#: give -- for det that was 9.6e-11 against a true minimum near 1.4e-6, four
#: orders low. A bound that far below the quantity it bounds is a weak
#: certificate even though it is a true one (the same point ENG-005 made about
#: E2). Halving the box radius halves the enclosure width, so a finer cover is
#: used for the reported bounds: at r ~ 4.4e-6 the width is ~6.6e-7 and the
#: bound lands near 1.1e-6, within a factor of 1.3 of the truth.
BOUND_CELLS = {"O1": 8192, "det_odd3": 32768}


def _common(status_ok: bool, work_order: str, quick: bool) -> dict:
    return {
        "certificate_version": "1.0",
        "program": "RH/Weil odd degree-3 block — Candidate A",
        "work_order": work_order,
        "claim_scope": D3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "pole_candidate": "A",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "domain": {"L_interval": list(D3.CELL_LABEL),
                   "L_left": repr(D3.CELL[0]), "L_right": repr(D3.CELL[1]),
                   "closed": True},
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "precision_bits": PRECISION_BITS,
        "mpmath_used": False,
        "quick_mode": quick,
    }


def build_scan(*, n_points: int, quick: bool) -> dict:
    _, arb, acb, _ = require_flint()
    set_precision_bits(PRECISION_BITS)
    scan = D3.topology_scan(arb, acb, n_points=n_points)
    body = _common(True, "ATLAS-RH-ENG-006 §8", quick)
    body.update({
        "evidence_class": "E3",
        "status": "SCANNED",
        "hard_constraints_certified": False,
        "content_kind": "WEIL_DEGREE3_ODD_SCAN",
        "basis": {"q1": "x - L/2", "b3": "x(L-x)(x-L/2)",
                  "parity": "both odd about x = L/2"},
        "scan": scan,
        "parity_identities": D3.parity_identities(),
        "note": ("E3 heuristic only. Eigenvalues here come from the closed 2x2 form "
                 "at grid points; they locate apparent features and can never "
                 "promote an E1 inertia claim (§14.4)."),
    })
    return body


def _evaluator(quantity: str, arb, acb, primes):
    def evaluate(lo: float, hi: float):
        blk = D3.odd_block(interval_box(lo, hi), arb, acb, prime_powers=primes)
        v = blk[quantity]
        return float(v.lower()), float(v.upper())

    return evaluate


def build_e1(scan_body: dict, *, quick: bool) -> tuple:
    """Certify the odd block's inertia, then bound what §9 asks for."""
    _, arb, acb, _ = require_flint()
    set_precision_bits(PRECISION_BITS)
    primes = WE.prime_powers_below(sum(D3.CELL) / 2)
    # Quick mode cannot simply use fewer boxes: separation genuinely needs a
    # radius near 1e-5, so a coarse start just splits its way back down and ends
    # up *slower* than the real run. It exercises the same code on a short
    # sub-interval instead, and is never promotable.
    cell = (1.30, 1.32) if quick else D3.CELL
    cells = 128 if quick else INITIAL_CELLS

    def fam(lo, hi):
        return D3.odd_matrix_at(interval_box(lo, hi), arb, acb, prime_powers=primes)

    strat = certify_inertia_family(
        fam, cell,
        subdivision_policy={"initial_cells": cells, "max_depth": 20})

    constant = strat.signature_if_constant()
    positive_definite = constant == (2, 0, 0)

    body = _common(True, "ATLAS-RH-ENG-006 §9", quick)
    strat_cert = build_stratification_certificate(
        strat, dimension=2, program=body["program"],
        work_order="ATLAS-RH-ENG-006 §9", evidence_class="E1",
        normalization_certificate_id=body[promotion.NORMALIZATION_ID_FIELD])

    bounds = {}
    if positive_definite:
        # §9 Outcome A names these two explicitly.
        for label, quantity in (("O1", "Oqq"), ("det_odd3", "det")):
            boxes = cells if quick else BOUND_CELLS[label]
            cov = IC.adaptive_cover(
                _evaluator(quantity, arb, acb, primes),
                quantity=f"degree3_{label}", cell=cell, target=0.0,
                initial_boxes=boxes, max_depth=20)
            bounds[label] = {"certified_lower_bound": repr(cov.certified_lower_bound),
                             "grid_minimum_for_scale": (
                                 scan_body.get("scan", {}).get("det_min_on_grid", {})
                                 .get("det") if label == "det_odd3" else None),
                             "cover": cov.to_dict()}

    ok = bool(positive_definite) and all(
        float(b["certified_lower_bound"]) > 0 for b in bounds.values())
    outcome = ("A_POSITIVE_DEFINITE" if positive_definite
               else "B_INERTIA_STRATIFICATION" if strat.status.startswith("PASS")
               else "C_INCONCLUSIVE")
    # What kind of artifact this is depends on what was established, and the
    # distinction is not cosmetic. §11 says an *inertia* certificate must never
    # satisfy a consumer requiring PSD -- the point being that "I know the
    # signature" should not be silently read as "it is positive". That rule
    # binds the inertia object, which is nested below and still refuses.
    #
    # Outcome A is a different claim: the block is positive definite everywhere
    # on the cell, with certified positive lower bounds on O1 and det. That is a
    # positivity certificate that happens to have been proved by an inertia
    # computation, so it carries a positivity content kind and answers a PSD
    # consumer honestly. Outcome B stays an inertia stratification and does not.
    kind = (POSITIVITY_KIND if positive_definite else strat_cert["content_kind"])
    body.update({
        "certified_cell": [repr(cell[0]), repr(cell[1])],
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if (ok or outcome == "B_INERTIA_STRATIFICATION") else "INCONCLUSIVE",
        "hard_constraints_certified": bool(ok),
        "promotion_state": promotion.PROMOTED_STATE if ok else "REFUSED",
        "outcome": outcome,
        "content_kind": kind,
        "psd_claim": bool(positive_definite and ok),
        # Top level so the PSD predicate can read the signature directly; also
        # repeated under "inertia" for readers who want it named.
        "n_positive": constant[0] if constant else None,
        "n_negative": constant[1] if constant else None,
        "n_zero": constant[2] if constant else None,
        "inertia": {"n_positive": constant[0] if constant else None,
                    "n_negative": constant[1] if constant else None,
                    "n_zero": constant[2] if constant else None,
                    "constant_on_cell": strat.is_constant},
        "inertia_stratification": strat_cert,
        "uniform_bounds": bounds,
        "claim": (
            "The odd degree-3 Weil block [[G_q1q1, G_q1b3], [G_q1b3, G_b3b3]] is "
            "positive definite for every L in [log 3, log 4]: inertia (2, 0, 0), "
            "with O1 and det(G_odd,3) both uniformly bounded below by the certified "
            "constants above."
            if positive_definite else
            "The odd degree-3 Weil block's inertia over [log 3, log 4], as a "
            "stratification with any unresolved transitions reported."),
        "strategy_chosen_from": ("the fresh Candidate-A E3 scan (§8); the scan is a "
                                 "clue, this cover is the warrant"),
        "scan_said": scan_body.get("scan", {}).get("distinct_floating_inertias"),
        "method": ("interval Hermitian LDL congruence on the assembled block, with "
                   "adaptive subdivision of the L cell; no eigenvalue solver, no "
                   "termwise PSD domination of the prime block"),
        "backend": backend_info(PRECISION_BITS).to_dict(),
    })
    if quick:
        body.update({"promotion_state": "REFUSED", "hard_constraints_certified": False,
                     "status": "QUICK_SMOKE_TEST_NOT_PROMOTABLE"})
    filename = POSITIVITY_FILE if positive_definite else INERTIA_FILE
    return filename, body


def build_moments(e1_body: dict, *, quick: bool) -> dict:
    """§10: moments and rank-trace, once the inertia result is in hand."""
    _, arb, acb, _ = require_flint()
    set_precision_bits(PRECISION_BITS)
    primes = WE.prime_powers_below(sum(D3.CELL) / 2)
    observed = [e1_body["inertia"]["n_positive"], e1_body["inertia"]["n_negative"],
                e1_body["inertia"]["n_zero"]]

    points = [("log3", D3.CELL[0]), ("det_min_1.2377", 1.2376586236864704),
              ("log4", D3.CELL[1])]
    rows = []
    for label, Lv in points:
        blk = D3.odd_block(arb(repr(Lv)), arb, acb, prime_powers=primes)
        G = D3.odd_matrix(blk)
        report = analyse(G, observed_inertia=observed if all(
            v is not None for v in observed) else None)

        m1 = report["moments"]["m1"]
        m2 = report["moments"]["m2"]
        lam = next(q for q in report["b1_queries"]
                   if q["query"] == "smallest_eigenvalue_bounds")
        lam_hi = float(lam["lambda_max"]["hi"])
        # The rank-trace inequality is not scale free: rank is scale invariant,
        # tr(P) is degree 1 and ||P||_HS^2 is degree 2. It therefore holds under a
        # normalization, and the one that makes it sharp is spectrum in [0, 1]
        # (a projection attains equality). That is checkable here, so it is
        # checked rather than declared.
        spectrum_ok = 0.0 <= float(lam["lambda_min"]["lo"]) and lam_hi <= 1.0
        hyp = {
            "P_positive_semidefinite": {
                "verified": observed[1] == 0 and observed[2] == 0,
                "evidence": {"from": e1_body["work_order"],
                             "inertia": observed,
                             "certificate": "this run's E1 inertia result"}},
            "Q_hermitian": {"verified": True,
                            "evidence": {"Q": "zero matrix", "trivially Hermitian": True}},
            "Q_positive_index_at_most_b": {
                "verified": True,
                "evidence": {"Q": "zero matrix", "positive_index": 0, "b": 0}},
            "shared_normalization": {
                "verified": bool(spectrum_ok),
                "evidence": {
                    "requirement": "spectrum of P contained in [0, 1]",
                    "lambda_min_lower": lam["lambda_min"]["lo"],
                    "lambda_max_upper": lam["lambda_max"]["hi"],
                    "why": ("the inequality is not scale free -- rank is scale "
                            "invariant while tr(P) is degree 1 and the HS term "
                            "degree 2 -- so it holds under a normalization, and "
                            "equality at a projection identifies that "
                            "normalization as spectrum in [0, 1]")}},
        }
        rt = rank_trace_lower_bound(
            trace_P=float(m1["lo"]), trace_Q=0.0,
            hs_sq_P_plus_Q=float(m2["hi"]), positive_index_Q_bound=0,
            hypotheses=hyp)
        rows.append({"label": label, "L": repr(Lv),
                     "moment_analysis": report,
                     "rank_trace": rt.to_dict()})

    body = _common(True, "ATLAS-RH-ENG-006 §10", quick)
    body.update({
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS",
        "hard_constraints_certified": True,
        "promotion_state": promotion.PROMOTED_STATE,
        "content_kind": "WEIL_SPECTRAL_MOMENT_CERTIFICATE",
        "dimension": 2,
        "inertia_dependency": e1_body.get("content_kind"),
        "points": rows,
        "moments": {r["label"]: r["moment_analysis"]["moments"] for r in rows},
        "b1_queries": [q for r in rows for q in r["moment_analysis"]["b1_queries"]],
        "note": ("Moments are enclosures of traces of matrix powers. "
                 "INSUFFICIENT_INFORMATION answers are certified results (§6), not "
                 "failures to compute."),
    })
    if quick:
        body.update({"promotion_state": "REFUSED", "hard_constraints_certified": False,
                     "status": "QUICK_SMOKE_TEST_NOT_PROMOTABLE"})
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    # Accepted for symmetry with the other certifiers and ignored on purpose:
    # every artifact here is new in ENG-006, so none of them was ever under the
    # WO-RH-17 quarantine and none has a marker to release.
    ap.add_argument("--release", action="store_true",
                    help="accepted and ignored; these artifacts were never quarantined")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--stage", choices=["scan", "e1", "moments", "all"], default="all")
    ap.add_argument("--scan-points", type=int, default=61)
    args = ap.parse_args()

    try:
        stages = ["scan", "e1", "moments"] if args.stage == "all" else [args.stage]

        if "scan" in stages:
            scan_body = build_scan(n_points=15 if args.quick else args.scan_points,
                                   quick=args.quick)
            print(f"wrote {write_certificate(SCAN_FILE, scan_body)}")
            sc = scan_body["scan"]
            print(f"  scan: inertias {sc['distinct_floating_inertias']} "
                  f"det_min {sc['det_min_on_grid']['det']} "
                  f"crossings {len(sc['det_zero_crossings'])}")
        else:
            scan_body = json.loads((CERT_DIR / SCAN_FILE).read_text(encoding="utf-8"))

        e1_body = None
        if "e1" in stages:
            name, e1_body = build_e1(scan_body, quick=args.quick)
            print(f"wrote {write_certificate(name, e1_body)}")
            print(f"  e1: outcome={e1_body['outcome']} status={e1_body['status']} "
                  f"inertia={tuple(e1_body['inertia'][k] for k in ('n_positive','n_negative','n_zero'))}")
            for k, v in e1_body["uniform_bounds"].items():
                print(f"      {k} >= {v['certified_lower_bound']} "
                      f"({v['cover']['boxes_examined']} boxes, "
                      f"depth {v['cover']['max_subdivision_depth']})")
        else:
            for cand in (POSITIVITY_FILE, INERTIA_FILE):
                if (CERT_DIR / cand).exists():
                    e1_body = json.loads((CERT_DIR / cand).read_text(encoding="utf-8"))
                    break

        if "moments" in stages:
            if e1_body is None:
                print("ERROR: the E1 stage must run before moments (§10)", file=sys.stderr)
                return 1
            body = build_moments(e1_body, quick=args.quick)
            print(f"wrote {write_certificate(MOMENTS_FILE, body)}")
            for r in body["points"]:
                rt = r["rank_trace"]
                got = rt.get("result", {}).get("certified_rank_lower_bound")
                print(f"      {r['label']:16} rank-trace: {rt['status']} bound={got}")
        return 0
    except FlintUnavailable as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
