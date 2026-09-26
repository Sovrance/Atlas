"""B15-POS3 tests — zero-pivot / degenerate (variance-zero) slices at Hankel order t = 2.

Registered in docs/preregistrations/prereg-005-pos-degenerate-slice.md (ACTIVE, frozen; see
prereg-005.freeze). Same system as prereg-004: H = [mu_{i+j}]_{0..2} >= 0,
B = [mu_{i+j+1} - mu_{i+j+2}]_{0..1} >= 0 (warrant [KrN, Theorem II.2.3] as cited in
Curto-Fialkow 1991 Remark 4.4; KrN chapter number VERIFY-BEFORE-USE). On mu2 = mu1^2 the
feasible set on the slice is the single point mu4 = a^4 (representing measure delta_a).

U1  primal: D-on, E-on -> PERMITTED (PSD, not PD) with the registered chains                  [K1]
U2  dual: five exterior points -> REJECTED; stop pivots, violated wall, v and v^T M v ==
    registered; zero-pivot Gram dual (rev. 2 clarification) certified                          [K1, K2]
U3  degenerate certified point: a^4 PERMITTED, a^4 +- 2^-20 REJECTED (certified_point field)    [K3]
U4  single-atom Farkas: UNIQUE (w = 1) at D-on/E-on, INCONSISTENT + verified at all five        [K4]
U5  formula controls: det A = 0 (L undefined) on D and E; U = a^4 on D; U undefined on E        [K5]
U6  no silent path: the prereg-004 routine refuses D-lo, E-lo, E-hi (assert) and D-off (ValueError)
U7  standalone re-verify + schema negatives
"""

import copy
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fractions import Fraction as F

from b15_surf import certificate as C
from b15_surf.exact import fmt
from b15_surf.pir_bridge import domain_fact
from b15_surf.pos import certify_t2 as P2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT_NAME = "b15_pos3_certificate.json"
CERT = os.path.join(ROOT, "certificates", CERT_NAME)
RESOLUTION = F(1, 2 ** 20)                  # prereg-005 P5-OI-4
SEED = 0                                    # no randomness in B15-POS3 (prereg-005)

# ---- registered values, transcribed from prereg-005 (do not edit after freeze) ----
SLICES = {"D": F(1, 2), "E": F(1)}
POINTS = {   # point: (slice, mu3, mu4, expected, H pivots, B pivots, violated wall)
    "D-on":  ("D", F(1, 8), F(1, 16), "PERMITTED", ["1", "0", "0"], ["1/4", "0"], None),
    "D-lo":  ("D", F(1, 8), F(3, 64), "REJECTED", ["1", "0", "-1/64"], ["1/4", "1/64"], "hankel"),
    "D-hi":  ("D", F(1, 8), F(5, 64), "REJECTED", ["1", "0", "1/64"], ["1/4", "-1/64"], "gap"),
    "D-off": ("D", F(9, 64), F(1, 16), "REJECTED", ["1", "0"], ["1/4", "31/1024"], "hankel"),
    "E-on":  ("E", F(1), F(1), "PERMITTED", ["1", "0", "0"], ["0", "0"], None),
    "E-lo":  ("E", F(1), F(63, 64), "REJECTED", ["1", "0", "-1/64"], ["0", "1/64"], "hankel"),
    "E-hi":  ("E", F(1), F(65, 64), "REJECTED", ["1", "0", "1/64"], ["0", "-1/64"], "gap"),
}
DUALS = {    # point: (v, v^T M v)
    "D-lo":  (["-1/4", "0", "1"], F(-1, 64)),
    "D-hi":  (["-1/2", "1"], F(-1, 64)),
    "D-off": (["63/4", "-32", "1"], F(-1)),
    "E-lo":  (["-1", "0", "1"], F(-1, 64)),
    "E-hi":  (["0", "1"], F(-1, 64)),
}
U3_EXPECTED_WALL = {"-": "hankel", "+": "gap"}   # D registered; E pattern identical (verified at U3)


def _mu(name):
    s, m3, m4 = POINTS[name][:3]
    a = SLICES[s]
    return a, a * a, m3, m4


def _chains(r):
    return {k: [F(p) for p in c["pivots"]] for k, c in r["primal"]["chains"].items()}


def u1_primal():
    out = {}
    for name, (s, m3, m4, exp, hp, bp, _) in POINTS.items():
        if exp != "PERMITTED":
            continue
        r = P2.certify_point_zp(*_mu(name))
        assert r["verdict"] == "PERMITTED", ("K1", name, r["verdict"])                        # K1
        assert _chains(r) == {"hankel": [F(x) for x in hp], "gap": [F(x) for x in bp]}, (name, _chains(r))
        for k in ("hankel", "gap"):
            assert r["primal"]["chains"][k]["status"] in ("PSD_CERTIFIED", "PD_CERTIFIED")
        assert r["primal"]["chains"]["hankel"]["status"] == "PSD_CERTIFIED"                     # not PD
        out[name] = {"status": "PASS", "registered_pivots_match": True, "psd_not_pd": True, **r}
        print(f"U1 primal {name}: PERMITTED (PSD, not PD), H {hp}, B {bp} (== registered)")
    return {"status": "PASS", "points": out}


def u2_dual():
    out = {}
    for name, (s, m3, m4, exp, hp, bp, wall) in POINTS.items():
        if exp != "REJECTED":
            continue
        r = P2.certify_point_zp(*_mu(name))
        assert r["verdict"] == "REJECTED", ("K1", name)                                        # K1
        assert _chains(r) == {"hankel": [F(x) for x in hp], "gap": [F(x) for x in bp]}, (name, _chains(r))
        other = "gap" if wall == "hankel" else "hankel"
        assert r["primal"]["chains"][other]["status"] in ("PSD_CERTIFIED", "PD_CERTIFIED"), \
            (name, "violates more than the registered wall")
        d = r["impossibility_certificate"]["dual_functional"]
        v, val = DUALS[name]
        assert d["op"] == "DUAL_EXCLUSION_FUNCTIONAL" and d["op_rev"] == 2
        assert d["matrix"] == wall, (name, d["matrix"], wall)
        assert [F(x) for x in d["direction"]] == [F(x) for x in v], ("K2", name, d["direction"], v)  # K2
        assert F(d["evaluation"]) == val, ("K2", name, d["evaluation"], val)                  # K2
        assert d["certified"] and d["identity_verified"], ("K2", name, d)                     # K2
        out[name] = {"status": "PASS", "violated_wall": wall, "registered_v": v,
                     "registered_value": fmt(val), **r}
        print(f"U2 dual {name}: REJECTED ({wall}; {d['lift']['stop_kind']} at {d['lift']['stop_index']}, "
              f"zero-pivot rows {d['lift']['zero_pivot_rows']}); v={d['direction']} -> {d['evaluation']} "
              f"(== registered); Gram certified")
    return {"status": "PASS", "points": out}


def u3_certified_point():
    out = {}
    for s, a in SLICES.items():
        sl = (a, a * a, a ** 3)
        on = P2.certify_point_zp(*sl, a ** 4)
        assert on["verdict"] == "PERMITTED", ("K3", s)
        neighbours = []
        for sign, delta in (("-", -RESOLUTION), ("+", RESOLUTION)):
            r = P2.certify_point_zp(*sl, a ** 4 + delta)
            assert r["verdict"] == "REJECTED", ("K3", s, sign)                                  # K3
            wall = r["impossibility_certificate"]["dual_functional"]["matrix"]
            assert wall == U3_EXPECTED_WALL[sign], (s, sign, wall)
            assert r["impossibility_certificate"]["dual_functional"]["certified"]
            neighbours.append({"mu4": fmt(a ** 4 + delta), "violated_wall": wall, "certificate": r})
        out[s] = {"permitted": fmt(a ** 4), "rejected_neighbours": neighbours}
        print(f"U3 {s}: mu4 = {fmt(a ** 4)} PERMITTED; {neighbours[0]['mu4']} fails "
              f"{neighbours[0]['violated_wall']}, {neighbours[1]['mu4']} fails {neighbours[1]['violated_wall']}")
    return {"status": "PASS", "resolution": fmt(RESOLUTION),
            "note": "feasible set on each slice is the point {a^4}; not a bracket (prereg-005 U3)",
            "slices": out}


def u4_farkas():
    rows = {}
    for name, (s, m3, m4, exp, *_ ) in POINTS.items():
        mu = P2.moments(*_mu(name))
        fk = P2.farkas_single_atom(SLICES[s], mu)
        if exp == "PERMITTED":
            assert fk["status"] == "UNIQUE" and fk["w"] == ["1/1"] and fk["verified"], ("K4", name, fk)  # K4
        else:
            assert fk["status"] == "INCONSISTENT" and fk["verified"], ("K4", name, fk)          # K4
        rows[name] = fk
    print("U4 single-atom system: UNIQUE (w = 1) at D-on, E-on; INCONSISTENT with verified Farkas "
          "vector at D-lo, D-hi, D-off, E-lo, E-hi")
    return {"status": "PASS", "role": "complete on the registered slices (premise: mu2 = mu1^2); "
                                     "Gram dual still mandatory (prereg-005 P5-OI-3)", "points": rows}


def u5_formula_controls():
    out = {}
    for s, a in SLICES.items():
        mu1, mu2, mu3, mu4 = a, a * a, a ** 3, a ** 4
        det_A = 1 * mu2 - mu1 * mu1
        assert det_A == 0, s                                   # L = v^T A^-1 v undefined
        gap = mu1 - mu2
        row = {"det_A": fmt(det_A), "L": "undefined (det A = 0)", "mu1_minus_mu2": fmt(gap)}
        if gap != 0:
            U = mu3 - (mu2 - mu3) ** 2 / gap
            assert U == a ** 4, (s, U)
            row["U"] = fmt(U)
        else:
            row["U"] = "undefined (mu1 - mu2 = 0)"
        out[s] = row
    assert out["D"]["U"] == "1/16" and out["E"]["U"].startswith("undefined")
    print(f"U5 formula controls: L undefined on D and E; U(D) = {out['D']['U']} = a^4; U undefined on E")
    return {"status": "PASS", "slices": out,
            "verdict_route": "certify_point_zp: SCHUR_PIVOT_EXACT pivots + zero-pivot Gram dual; "
                             "no L/U wall formula on the verdict path (K5); formulas computed here "
                             "only as controls"}


def u6_no_silent_path():
    rows = {}
    expect = {"D-lo": AssertionError, "E-lo": AssertionError, "E-hi": AssertionError, "D-off": ValueError}
    for name, exc in expect.items():
        try:
            P2.certify_point(*_mu(name))
        except exc as e:
            rows[name] = {"prereg004_routine": f"refuses ({type(e).__name__})"}
            continue
        except Exception as e:  # noqa: BLE001 - any other outcome is a finding
            raise AssertionError((name, "unexpected exception", type(e).__name__)) from e
        raise AssertionError((name, "prereg-004 routine emitted a verdict through an unregistered path"))
    print("U6 no silent path: prereg-004 routine refuses D-lo, E-lo, E-hi (AssertionError) and "
          "D-off (ValueError), as registered; B15-POS3 certifies all seven via the registered construction")
    return {"status": "PASS", "points": rows}


def u7_standalone(cert_path):
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"),
                           cert_path], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    line = proc.stdout.strip().splitlines()[-1]
    assert line.startswith("VERIFIED "), line
    print("U7 standalone re-verify: " + line)
    return {"status": "PASS", "tool": "tools/verify_b15_certificate.py", "verified": True,
            "verdict_line": line.split("(", 1)[1].rstrip(")")}


def u7_schema_negatives(cert):
    bad = []
    for name, mutate in (
        ("non-SPEC verdict", lambda c: c.update(verdict="CONDITIONAL(X)")),
        ("HEURISTIC asserting E0", lambda c: c.update(soundness="HEURISTIC", warnings=[{"location": "x", "text": "y"}])),
        ("POS3 REJECTED without impossibility certificate",
         lambda c: c.update(verdict="REJECTED", impossibility_certificate=None)),
        ("wrong preregistration binding", lambda c: c.update(prereg_ref=c["prereg_ref"].replace("prereg-005", "prereg-004"))),
        ("certified_inner_interval on a point-valued slice",
         lambda c: c["results"]["U3_certified_point"].update(certified_inner_interval=["1/16", "1/16"])),
    ):
        c2 = copy.deepcopy(cert)
        mutate(c2)
        c2["certificate_id"] = f"b15-pos3-{C.content_hash(c2)}"
        try:
            C.validate(c2)
            raise AssertionError(f"negative test did not fire: {name}")
        except C.B15CertificateError as e:
            bad.append({"case": name, "rejected_with": str(e)[:80]})
    dual = lambda c, n: c["results"]["U2_dual"]["points"][n]["impossibility_certificate"]
    for name, mutate in (
        ("tampered pivot", lambda c: c["results"]["U1_primal"]["points"]["D-on"]["primal"]["chains"]["hankel"]["pivots"].__setitem__(2, "-1/64")),
        ("tampered lift (zero-pivot rows)", lambda c: dual(c, "D-lo")["dual_functional"]["lift"].__setitem__("zero_pivot_rows", [])),
        ("tampered direction v", lambda c: dual(c, "D-off")["dual_functional"].__setitem__("direction", ["63/4", "-31", "1"])),
        ("tampered Farkas vector", lambda c: dual(c, "E-hi")["farkas_atom"].__setitem__("farkas_y", ["0", "0", "0", "0", "0"])),
        ("certified_inner_interval on a point-valued slice",
         lambda c: c["results"]["U3_certified_point"].update(certified_inner_interval=["1/16", "1/16"])),
    ):
        c3 = copy.deepcopy(cert)
        mutate(c3)
        # re-hash so the content-hash check passes and the substantive check must fire
        c3["certificate_id"] = f"b15-pos3-{C.content_hash(c3)}"
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(c3, f)
            tmp = f.name
        proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"), tmp],
                              capture_output=True, text=True, cwd=ROOT)
        os.unlink(tmp)
        assert proc.returncode == 1, (name, proc.stdout)
        lines = [x.strip()[2:] for x in proc.stdout.strip().splitlines()[1:]]
        assert not any("content hash" in x for x in lines), (name, lines)
        bad.append({"case": f"{name} -> standalone tool FAIL (re-hashed)", "rejected_with": "; ".join(lines)[:160]})
    print(f"U7 schema negatives: {len(bad)}/10 rejected as required")
    return {"status": "PASS", "cases": bad}


if __name__ == "__main__":
    r1, r2, r3, r4 = u1_primal(), u2_dual(), u3_certified_point(), u4_farkas()
    r5, r6 = u5_formula_controls(), u6_no_silent_path()
    results = {"U1_primal": r1, "U2_dual": r2, "U3_certified_point": r3, "U4_farkas": r4,
               "U5_formula_controls": r5, "U6_no_silent_path": r6}
    witness = {name: r1["points"][name]["witness"] for name in r1["points"]}
    registered = {
        "slices": {s: {"a": fmt(a)} for s, a in SLICES.items()},
        "points": {n: {"slice": p[0], "a": fmt(SLICES[p[0]]), "mu3": fmt(p[1]), "mu4": fmt(p[2]),
                       "expected": p[3], "hankel_pivots": p[4], "gap_pivots": p[5], "violated": p[6],
                       **({"matrix": p[6], "v": DUALS[n][0], "value": fmt(DUALS[n][1])} if n in DUALS else {})}
                   for n, p in POINTS.items()},
        "resolution": fmt(RESOLUTION)}
    warnings = [{"location": "b15_surf/pos/certify_t2.py::moments",
                 "text": "moment sequence truncated at a_{6,0} (Hankel order t=2); higher-order "
                         "EFT-hedron walls are not tested (HEURISTIC choice)"}]
    fact = domain_fact(
        benchmark="B15-POS3",
        content={"subject": "zero-pivot / degenerate-slice verification at Hankel order t=2 "
                            "(variance-zero slices of the degree-4 Hausdorff problem)",
                 "certified_points": {s: r3["slices"][s]["permitted"] for s in SLICES},
                 "resolution": fmt(RESOLUTION)},
        verdict="PERMITTED", evidence_level="E0", soundness="SOUND",
        witness=witness, impossibility_certificate=None, assumptions=[], warnings=[],
        measurement_interface="published bound arXiv:2012.15849 eqs. 7.45/7.46 (no apparatus)")
    cert = C.build(
        benchmark="B15-POS3",
        problem="B15-POS3: exact-PSD verifier on variance-zero slices at Hankel order t = 2 — zero "
                "pivots (vanishing and non-vanishing rows), degenerate gap wall, undefined wall formulas",
        headline="D-on and E-on PERMITTED (PSD, not PD); all five exterior points REJECTED with the "
                 "registered zero-pivot Gram duals (incl. zero-diagonal indefiniteness at D-off, "
                 "v^T H v = -1) and verified single-atom Farkas vectors; feasible set on each slice "
                 "is the point a^4 at resolution 2^-20",
        certificate_class="EXACT-RATIONAL (same class as B1/B2/B7/B15-POS/B15-POS2); framework capability",
        results=results, soundness="SOUND", warnings=warnings,
        ground_truth_route="closed form: variance zero => unique measure delta_a, feasible set {a^4}; "
                           "system warrant [KrN, Theorem II.2.3] as cited in Curto-Fialkow 1991 Rem. 4.4 "
                           "(prereg-004 errata E1; KrN number VERIFY-BEFORE-USE)",
        evidence_level="E0", pir_level="L2", assumptions=[],
        falsifier_direction="K1/K2: a wrong on-slice verdict or an exterior point without the registered "
                            "certified Gram dual files zero-pivot handling REJECTED (asm:B15-POS3-unverified); "
                            "K3/K4 stop certification; K5 rejects any verdict reached via the L/U formulas",
        seed=SEED, generator_sha256=C.content_hash({"module": "b15_surf.pos.certify_t2::certify_point_zp",
                                                    "registered": registered}),
        m_layer_stipulations=[
            "Bound is a theorem-level statement about EFT coefficients; no apparatus, no data",
            "Moments truncated at a_{6,0} (HEURISTIC truncation, see warnings)",
            "Slices D (a = 1/2) and E (a = 1) impose mu2 = mu1^2 (variance zero); registered in prereg-005",
            "Farkas completeness premise: variance zero => unique candidate delta_a (registered fact, not asm:)",
            "Pre-freeze out-of-band check disclosed in prereg-005 (values unchanged)",
        ],
        calibration_route="none required (exact rational arithmetic; resolution registered)",
        verdict="PERMITTED", witness=witness, impossibility_certificate=None,
        inputs={"registered": registered,
                "prereg": "docs/preregistrations/prereg-005-pos-degenerate-slice.md"},
        pir_facts=[fact.to_dict()],
    )
    C.save(cert, CERT_NAME)
    r7 = u7_standalone(CERT)
    cert["results"]["U7_standalone_reverify"] = r7
    cert["results"]["U7_schema_negatives"] = u7_schema_negatives(cert)
    cert["certificate_id"] = f"b15-pos3-{C.content_hash(cert)}"
    C.save(cert, CERT_NAME)
    print(f"\nCertificate written: certificates/{CERT_NAME}")
    print("ALL B15-POS3 TESTS PASS")
