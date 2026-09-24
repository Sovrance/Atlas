"""B15-POS tests — exact-PSD verifier vs. a published EFT-hedron bound (R26).

Registered target (prereg-002 §POS; §9-OI-1 pending): forward-limit moments
(mu0, mu1, mu2) = (1, a30/a20, a40/a20) with M_Gap = 1, at the slice mu1 = 1/2,
lower wall mu2 >= mu1^2 = 1/4 (Hankel total positivity, arXiv:2012.15849 eq. 7.19)
and upper wall mu2 <= mu1 = 1/2 (gap chain, eq. 7.35).

T1  primal: registered interior point -> all Schur pivots >= 0 -> PERMITTED, witness = pivots
T2  dual:   registered exterior points -> negative pivot AND an exact exclusion functional
            (nonnegativity identity verified via pir.symbolic.linear) -> REJECTED
T3  bracket: exact bisection on both registered rays; the published walls 1/4 and 1/2
            lie inside brackets of width <= 2^-20
T4  negative control: perturb mu2 by the registered rational -> verdict flips with witness
T5  standalone re-verify: tools/verify_b15_certificate.py recomputes everything from the JSON
T6  schema negative tests (§5): non-SPEC verdict, HEURISTIC+E0, POS REJECTED without an
    impossibility certificate are all rejected; a tampered certificate fails the standalone tool
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fractions import Fraction as F

from b15_surf import certificate as C
from b15_surf.exact import fmt
from b15_surf.pir_bridge import domain_fact
from b15_surf.pos import certify as P
from b15_surf.pos.dispersive import (EXTERIOR_POINT_GAP, EXTERIOR_POINT_HANKEL,
                                     INTERIOR_POINT, MU1_SLICE, PERTURBATION,
                                     PUBLISHED_LOWER_WALL, PUBLISHED_UPPER_WALL,
                                     warnings_for_truncation)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT = os.path.join(ROOT, "certificates", "b15_pos_certificate.json")
WIDTH = F(1, 2 ** 20)          # registered bracket width (proposal, §9-OI-2)
SEED = 20260923


def t1_primal():
    r = P.certify_point(*INTERIOR_POINT)
    assert r["verdict"] == "PERMITTED", r
    for chain in r["primal"]["chains"].values():
        assert all(F(p) >= 0 for p in chain["pivots"])
    print(f"T1 primal PERMITTED at {r['point']}: pivots={r['witness']['pivots']}")
    return {"status": "PASS", **r}


def t2_dual():
    out = {}
    for name, pt in (("hankel_wall", EXTERIOR_POINT_HANKEL), ("gap_wall", EXTERIOR_POINT_GAP)):
        r = P.certify_point(*pt)
        assert r["verdict"] == "REJECTED", r
        neg = [F(p) for ch in r["primal"]["chains"].values() for p in ch["pivots"] if F(p) < 0]
        assert neg, "no negative pivot"
        ic = r["impossibility_certificate"]
        assert ic["dual_functional"]["identity_verified"], ic
        assert F(ic["dual_functional"]["evaluation"]) < 0
        if name == "hankel_wall":
            assert ic["farkas_atom"]["verify_farkas"] is True
        out[name] = {"status": "PASS", **r}
        print(f"T2 dual REJECTED at {r['point']} [{ic['dual_functional']['wall']}]: "
              f"negative pivot {fmt(neg[0])}, functional y={ic['dual_functional']['y']} "
              f"-> {ic['dual_functional']['evaluation']}")
    return out


def t3_bracket():
    lo = P.bisect_ray(MU1_SLICE, INTERIOR_POINT[1], F(0), WIDTH)
    hi = P.bisect_ray(MU1_SLICE, INTERIOR_POINT[1], F(1), WIDTH)
    for br, wall in ((lo, PUBLISHED_LOWER_WALL), (hi, PUBLISHED_UPPER_WALL)):
        a, b = (F(x) for x in br["certified_inner_interval"])
        assert a <= wall <= b, (br, wall)
        assert F(br["width"]) <= WIDTH
        br["published_wall"] = fmt(wall)
        br["contains_published"] = True
    print(f"T3 brackets: lower {lo['certified_inner_interval']} ∋ 1/4, "
          f"upper {hi['certified_inner_interval']} ∋ 1/2 (width <= 2^-20)")
    return {"status": "PASS", "lower_wall": lo, "upper_wall": hi, "width_threshold": fmt(WIDTH)}


def t4_negative_control():
    mu1, mu2 = INTERIOR_POINT
    r = P.certify_point(mu1, mu2 + PERTURBATION)
    assert r["verdict"] == "REJECTED", r
    assert r["impossibility_certificate"]["dual_functional"]["identity_verified"]
    print(f"T4 negative control: mu2 {fmt(mu2)} -> {fmt(mu2 + PERTURBATION)} flips to REJECTED "
          f"with witness pivots {r['witness']['pivots']['hankel']}")
    return {"status": "PASS", "perturbation": fmt(PERTURBATION), **r}


def t5_standalone(cert_path):
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"),
                           cert_path], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    line = proc.stdout.strip().splitlines()[-1]
    assert line.startswith("VERIFIED "), line
    print("T5 standalone re-verify: " + line)
    # The intermediate certificate id is not stored (it changes once T5/T6 results are
    # appended and the final id is recomputed); the verification outcome is.
    return {"status": "PASS", "tool": "tools/verify_b15_certificate.py", "verified": True,
            "verdict_line": line.split("(", 1)[1].rstrip(")")}


def t6_schema_negatives(cert):
    import copy
    import json
    import tempfile
    bad = []
    for name, mutate in (
        ("non-SPEC verdict", lambda c: c.update(verdict="CONDITIONAL(X)")),
        ("HEURISTIC asserting E0", lambda c: c.update(soundness="HEURISTIC", warnings=[{"location": "x", "text": "y"}])),
        ("POS REJECTED without impossibility certificate",
         lambda c: c.update(verdict="REJECTED", impossibility_certificate=None)),
        ("assumption without asm: prefix", lambda c: c.update(assumptions=["GRH"])),
    ):
        c2 = copy.deepcopy(cert)
        mutate(c2)
        c2["certificate_id"] = f"b15-pos-{C.content_hash(c2)}"
        try:
            C.validate(c2)
            raise AssertionError(f"negative test did not fire: {name}")
        except C.B15CertificateError as e:
            bad.append({"case": name, "rejected_with": str(e)[:80]})
    # tampered certificate must fail the standalone tool (Score = infinity)
    c3 = copy.deepcopy(cert)
    c3["results"]["T1_primal_permitted"]["primal"]["chains"]["hankel"]["pivots"][1] = "-1/8"
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(c3, f)
        tmp = f.name
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"), tmp],
                          capture_output=True, text=True, cwd=ROOT)
    os.unlink(tmp)
    assert proc.returncode == 1, proc.stdout
    bad.append({"case": "tampered pivot -> standalone tool FAIL", "rejected_with": proc.stdout.strip().splitlines()[0]})
    print(f"T6 schema negatives: {len(bad)}/5 rejected as required")
    return {"status": "PASS", "cases": bad}


if __name__ == "__main__":
    r1 = t1_primal()
    r2 = t2_dual()
    r3 = t3_bracket()
    r4 = t4_negative_control()
    results = {"T1_primal_permitted": r1, "T2_dual_rejected": r2, "T3_bracket": r3,
               "T4_negative_control": r4}
    warnings = warnings_for_truncation()
    fact = domain_fact(
        benchmark="B15-POS",
        content={"subject": "EFT-hedron forward-limit two-sided bound, slice mu1=1/2",
                 "bracket_lower": r3["lower_wall"]["certified_inner_interval"],
                 "bracket_upper": r3["upper_wall"]["certified_inner_interval"],
                 "published_walls": ["1/4", "1/2"]},
        verdict="PERMITTED", evidence_level="E0", soundness="SOUND",
        witness=r1["witness"], impossibility_certificate=None, assumptions=[],
        warnings=[], measurement_interface="published bound arXiv:2012.15849 eqs. 7.19/7.35 (no apparatus)")
    cert = C.build(
        benchmark="B15-POS",
        problem="B15-POS: exact-PSD verifier (SCHUR_PIVOT_EXACT + dual exclusion functional) "
                "scored against the published EFT-hedron forward-limit bound (R26)",
        headline="published two-sided bound reproduced: interior PERMITTED (primal pivots), "
                 "exterior REJECTED (negative pivot + exact dual functional), both walls "
                 "bracketed to <= 2^-20; negative control flips",
        certificate_class="EXACT-RATIONAL (same class as B1/B2/B7); published-bound reproduction",
        results=results, soundness="SOUND", warnings=warnings,
        ground_truth_route="published-bound:arXiv:2012.15849 eq.7.19 (Hankel) + eq.7.35 (gap)",
        evidence_level="E0", pir_level="L2",
        assumptions=[],
        falsifier_direction="F1: a bracket that excludes a published wall, or a REJECTED point "
                            "without a verified dual functional, files the verifier REJECTED "
                            "for external bounds and taints all B15 facts asm:B15-POS-unverified",
        seed=SEED, generator_sha256=C.content_hash({"module": "b15_surf.pos", "registered": {
            "interior": [fmt(x) for x in INTERIOR_POINT],
            "exterior": [[fmt(x) for x in EXTERIOR_POINT_HANKEL], [fmt(x) for x in EXTERIOR_POINT_GAP]],
            "width": fmt(WIDTH), "perturbation": fmt(PERTURBATION)}}),
        m_layer_stipulations=[
            "Bound is a theorem-level statement about EFT coefficients; no apparatus, no data",
            "Moments truncated at a_{4,0} (HEURISTIC truncation, see warnings)",
            "Units a_{2,0} = 1, M_Gap = 1; slice mu1 = 1/2 registered in prereg-002",
            "Target bound confirmed as registered (§9-OI-1, prereg-003 A3); bracket width 2^-20 accepted (§9-OI-2)",
        ],
        calibration_route="none required (exact rational arithmetic; bracket width registered)",
        verdict="PERMITTED", witness=r1["witness"], impossibility_certificate=None,
        inputs={"interior_point": [fmt(x) for x in INTERIOR_POINT],
                "exterior_points": {"hankel_wall": [fmt(x) for x in EXTERIOR_POINT_HANKEL],
                                    "gap_wall": [fmt(x) for x in EXTERIOR_POINT_GAP]},
                "bracket_width_threshold": fmt(WIDTH), "perturbation": fmt(PERTURBATION),
                "prereg_amendment": C.amendment_ref()},
        pir_facts=[fact.to_dict()],
    )
    # T5 is scored on the saved certificate itself (standalone tool).
    C.save(cert, "b15_pos_certificate.json")
    r5 = t5_standalone(CERT)
    cert["results"]["T5_standalone_reverify"] = r5
    cert["results"]["T6_schema_negatives"] = t6_schema_negatives(cert)
    cert["certificate_id"] = f"b15-pos-{C.content_hash(cert)}"
    C.save(cert, "b15_pos_certificate.json")
    print("\nCertificate written: certificates/b15_pos_certificate.json")
    print("ALL B15-POS TESTS PASS")
