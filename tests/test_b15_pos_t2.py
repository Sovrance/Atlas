"""B15-POS2 tests — exact-PSD verifier vs. the EFT-hedron gapped bound at Hankel order t = 2.

Registered in docs/preregistrations/prereg-004-pos-hankel-t2.md (ACTIVE, frozen; see
prereg-004.freeze). Target: the even truncated Hausdorff system of degree 4 on [0, 1]
(arXiv:2012.15849 eqs. 7.45/7.46; warrant [KrN, Theorem II.2.3] as cited in Curto-Fialkow
1991 Remark 4.4 — see prereg-004 errata E1; certificate strings keep the frozen wording):  H = [mu_{i+j}]_{0..2} >= 0,  B = [mu_{i+j+1} - mu_{i+j+2}]_{0..1} >= 0.

T1  primal: interior points S1-int, S2-int -> PERMITTED; pivot chains == registered
T2  dual: four exterior points -> REJECTED; negative pivot == registered; Gram-form
    DUAL_EXCLUSION_FUNCTIONAL rev. 2 witness certified (identity + Q0, Q1 >= 0)      [G2, G3]
T3  brackets: mu_4 bisection on each slice to width <= 2^-20; each contains its wall  [G1]
T4  quadrature route: S1 walls == Gauss-Legendre 2-pt and Simpson mu_4                 [G4]
T5  t = 1 blindness control: the prereg-002 verifier PERMITS every (mu1, mu2) projection
T6  negative control: mu_4 -> mu_4 - 1/36 at S1-int flips to REJECTED
T7  consistency with the paper's eq. (7.46) tower: all contiguous minors > 0 at interior points [G5]
T8  standalone re-verify + schema negatives
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
from b15_surf.pos import certify as P1
from b15_surf.pos import certify_t2 as P2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT_NAME = "b15_pos_t2_certificate.json"
CERT = os.path.join(ROOT, "certificates", CERT_NAME)
WIDTH = F(1, 2 ** 20)                       # prereg-004 P4-OI-4
NEG_CONTROL = F(-1, 36)                     # prereg-004 T6
SEED = 0                                    # no randomness in B15-POS2 (prereg-004)

# ---- registered values, transcribed from prereg-004 (do not edit after freeze) ----
SLICES = {
    "S1": {"slice": (F(1, 2), F(1, 3), F(1, 4)), "L": F(7, 36), "U": F(5, 24),
           "measure": "Lebesgue on [0, 1]"},
    "S2": {"slice": (F(7, 12), F(7, 16), F(73, 192)), "L": F(631, 1792), "U": F(641, 1792),
           "measure": "atoms {1/4, 1/2, 1}, weights 1/3 each"},
}
INTERIOR = {
    "S1-int": ("S1", F(1, 5), {"hankel": ["1/1", "1/12", "1/180"], "gap": ["1/6", "1/120"]}),
    "S2-int": ("S2", F(91, 256), {"hankel": ["1/1", "7/72", "3/896"], "gap": ["7/48", "1/448"]}),
}
EXTERIOR = {   # point: (slice, mu4, violated chain, registered negative pivot, which wall)
    "S1-lo": ("S1", F(13, 72), "hankel", F(-1, 72), "L"),
    "S1-hi": ("S1", F(11, 48), "gap", F(-1, 48), "U"),
    "S2-lo": ("S2", F(45, 128), "hankel", F(-1, 1792), "L"),
    "S2-hi": ("S2", F(321, 896), "gap", F(-1, 1792), "U"),
}


def _pt(name, table):
    s, m4 = table[name][0], table[name][1]
    return (*SLICES[s]["slice"], m4)


def t1_primal():
    out = {}
    for name, (s, m4, reg) in INTERIOR.items():
        r = P2.certify_point(*SLICES[s]["slice"], m4)
        assert r["verdict"] == "PERMITTED", (name, r)
        got = {k: [fmt(F(p)) for p in c["pivots"]] for k, c in r["primal"]["chains"].items()}
        want = {k: [fmt(F(p)) for p in v] for k, v in reg.items()}
        assert got == want, (name, got, want)
        out[name] = {"status": "PASS", "registered_pivots_match": True, **r}
        print(f"T1 primal {name}: PERMITTED, pivots H {got['hankel']}, B {got['gap']} (== registered)")
    return {"status": "PASS", "points": out}


def t2_dual():
    out = {}
    for name, (s, m4, chain, neg, wall) in EXTERIOR.items():
        r = P2.certify_point(*SLICES[s]["slice"], m4)
        assert r["verdict"] == "REJECTED", (name, r)
        piv = [F(p) for p in r["primal"]["chains"][chain]["pivots"]]
        assert piv[-1] == neg < 0, (name, piv, neg)
        other = "gap" if chain == "hankel" else "hankel"
        assert r["primal"]["chains"][other]["status"] in ("PSD_CERTIFIED", "PD_CERTIFIED"), \
            (name, "violates more than the registered wall")
        d = r["impossibility_certificate"]["dual_functional"]
        assert d["op"] == "DUAL_EXCLUSION_FUNCTIONAL" and d["op_rev"] == 2
        assert d["certified"] and d["identity_verified"], (name, d)             # G2
        assert F(d["evaluation"]) == neg, (name, d["evaluation"], neg)
        out[name] = {"status": "PASS", "violated_wall": wall, "registered_negative_pivot": fmt(neg), **r}
        print(f"T2 dual {name}: REJECTED, {chain} pivot {fmt(neg)} (== registered); "
              f"y={d['y']} -> {d['evaluation']}; Gram Q0 {d['Q0_status']}, Q1 {d['Q1_status']}")
    return {"status": "PASS", "points": out}


def t3_brackets():
    out = {}
    for s, spec in SLICES.items():
        sl = spec["slice"]
        inner = INTERIOR[f"{s}-int"][1]
        lo = P2.bisect_mu4(sl, inner, EXTERIOR[f"{s}-lo"][1], WIDTH)
        hi = P2.bisect_mu4(sl, inner, EXTERIOR[f"{s}-hi"][1], WIDTH)
        for br, wall in ((lo, spec["L"]), (hi, spec["U"])):
            a, b = (F(x) for x in br["certified_inner_interval"])
            assert a <= wall <= b, ("G1", s, br, wall)                            # G1
            assert F(br["width"]) <= WIDTH
            assert wall not in (a, b), "bisection landed on a wall (non-dyadic walls expected)"
            br["registered_wall"] = fmt(wall)
            br["contains_registered_wall"] = True
        out[s] = {"lower_wall": lo, "upper_wall": hi}
        print(f"T3 {s}: lower {lo['certified_inner_interval']} ∋ {fmt(spec['L'])} ({lo['steps']} steps); "
              f"upper {hi['certified_inner_interval']} ∋ {fmt(spec['U'])} ({hi['steps']} steps)")
    return {"status": "PASS", "width_threshold": fmt(WIDTH), "slices": out}


def t4_quadrature():
    g, sp = P2.gauss_legendre2_mu4(), P2.simpson_mu4()
    assert g == SLICES["S1"]["L"] and sp == SLICES["S1"]["U"], ("G4", g, sp)       # G4
    print(f"T4 quadrature: Gauss-Legendre 2-pt mu4 = {fmt(g)} == L(S1); Simpson mu4 = {fmt(sp)} == U(S1)")
    return {"status": "PASS", "gauss_legendre_2pt_mu4": fmt(g), "simpson_mu4": fmt(sp),
            "S1_L": fmt(SLICES["S1"]["L"]), "S1_U": fmt(SLICES["S1"]["U"])}


def t5_t1_blindness():
    rows = {}
    for name, table in [(n, INTERIOR) for n in INTERIOR] + [(n, EXTERIOR) for n in EXTERIOR]:
        mu1, mu2 = SLICES[table[name][0]]["slice"][:2]
        r1 = P1.certify_point(mu1, mu2)
        assert r1["verdict"] == "PERMITTED", (name, r1)
        rows[name] = {"t1_verdict": r1["verdict"], "projection": [fmt(mu1), fmt(mu2)],
                      "t2_verdict": "PERMITTED" if name in INTERIOR else "REJECTED"}
    print("T5 t=1 blindness: t=1 verifier PERMITS all 6 projections; t=2 rejects the 4 exterior points")
    return {"status": "PASS", "points": rows}


def t6_negative_control():
    s, m4, _ = INTERIOR["S1-int"]
    r = P2.certify_point(*SLICES[s]["slice"], m4 + NEG_CONTROL)
    assert r["verdict"] == "REJECTED", r
    assert r["impossibility_certificate"]["dual_functional"]["certified"]
    print(f"T6 negative control: mu4 {fmt(m4)} -> {fmt(m4 + NEG_CONTROL)} flips to REJECTED")
    return {"status": "PASS", "perturbation": fmt(NEG_CONTROL), **r}


def t7_tower():
    out = {}
    for name, (s, m4, _) in INTERIOR.items():
        minors = P2.tower_minors(P2.moments(*SLICES[s]["slice"], m4))
        assert all(F(m["minor"]) > 0 for m in minors), ("G5", name)               # G5
        out[name] = {"n_minors": len(minors), "min_minor": fmt(min(F(m["minor"]) for m in minors)),
                     "minors": minors}
        print(f"T7 tower consistency {name}: {len(minors)} contiguous minors, all > 0 "
              f"(min {out[name]['min_minor']})")
    return {"status": "PASS", "role": "consistency with eq. (7.46); warrant is KrN III.2.3 (prereg-004 P4-OI-1)",
            "points": out}


def t8_standalone(cert_path):
    proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"),
                           cert_path], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    line = proc.stdout.strip().splitlines()[-1]
    assert line.startswith("VERIFIED "), line
    print("T8 standalone re-verify: " + line)
    return {"status": "PASS", "tool": "tools/verify_b15_certificate.py", "verified": True,
            "verdict_line": line.split("(", 1)[1].rstrip(")")}


def t8_schema_negatives(cert):
    bad = []
    for name, mutate in (
        ("non-SPEC verdict", lambda c: c.update(verdict="CONDITIONAL(X)")),
        ("HEURISTIC asserting E0", lambda c: c.update(soundness="HEURISTIC", warnings=[{"location": "x", "text": "y"}])),
        ("POS2 REJECTED without impossibility certificate",
         lambda c: c.update(verdict="REJECTED", impossibility_certificate=None)),
        ("wrong preregistration binding", lambda c: c.update(prereg_ref=c["prereg_ref"].replace("prereg-004", "prereg-002"))),
    ):
        c2 = copy.deepcopy(cert)
        mutate(c2)
        c2["certificate_id"] = f"b15-pos2-{C.content_hash(c2)}"
        try:
            C.validate(c2)
            raise AssertionError(f"negative test did not fire: {name}")
        except C.B15CertificateError as e:
            bad.append({"case": name, "rejected_with": str(e)[:80]})
    for name, mutate in (
        ("tampered pivot", lambda c: c["results"]["T1_primal"]["points"]["S1-int"]["primal"]["chains"]["hankel"]["pivots"].__setitem__(2, "-1/180")),
        ("tampered Gram Q1", lambda c: c["results"]["T2_dual"]["points"]["S1-hi"]["impossibility_certificate"]["dual_functional"]["gram"]["Q1"][0].__setitem__(0, "-1/1")),
        ("bracket excluding wall", lambda c: c["results"]["T3_brackets"]["slices"]["S2"]["upper_wall"].__setitem__("registered_wall", "1/2")),
    ):
        c3 = copy.deepcopy(cert)
        mutate(c3)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(c3, f)
            tmp = f.name
        proc = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_b15_certificate.py"), tmp],
                              capture_output=True, text=True, cwd=ROOT)
        os.unlink(tmp)
        assert proc.returncode == 1, (name, proc.stdout)
        bad.append({"case": f"{name} -> standalone tool FAIL", "rejected_with": proc.stdout.strip().splitlines()[1][:80]})
    print(f"T8 schema negatives: {len(bad)}/7 rejected as required")
    return {"status": "PASS", "cases": bad}


if __name__ == "__main__":
    r1, r2, r3, r4 = t1_primal(), t2_dual(), t3_brackets(), t4_quadrature()
    r5, r6, r7 = t5_t1_blindness(), t6_negative_control(), t7_tower()
    results = {"T1_primal": r1, "T2_dual": r2, "T3_brackets": r3, "T4_quadrature": r4,
               "T5_t1_blindness": r5, "T6_negative_control": r6, "T7_tower_consistency": r7}
    witness = {name: r1["points"][name]["witness"] for name in INTERIOR}
    warnings = [{"location": "b15_surf/pos/certify_t2.py::moments",
                 "text": "moment sequence truncated at a_{6,0} (Hankel order t=2); higher-order "
                         "EFT-hedron walls are not tested (HEURISTIC choice)"}]
    brackets = {s: {k: v["certified_inner_interval"] for k, v in r3["slices"][s].items()} for s in SLICES}
    fact = domain_fact(
        benchmark="B15-POS2",
        content={"subject": "EFT-hedron gapped forward-limit bound at Hankel order t=2 (degree-4 Hausdorff)",
                 "brackets": brackets,
                 "registered_walls": {s: [fmt(v["L"]), fmt(v["U"])] for s, v in SLICES.items()}},
        verdict="PERMITTED", evidence_level="E0", soundness="SOUND",
        witness=witness, impossibility_certificate=None, assumptions=[], warnings=[],
        measurement_interface="published bound arXiv:2012.15849 eqs. 7.45/7.46 (no apparatus)")
    registered = {"slices": {s: {"slice": [fmt(x) for x in v["slice"]], "L": fmt(v["L"]), "U": fmt(v["U"])}
                             for s, v in SLICES.items()},
                  "interior": {n: fmt(v[1]) for n, v in INTERIOR.items()},
                  "exterior": {n: fmt(v[1]) for n, v in EXTERIOR.items()},
                  "width": fmt(WIDTH), "negative_control": fmt(NEG_CONTROL)}
    cert = C.build(
        benchmark="B15-POS2",
        problem="B15-POS2: exact-PSD verifier (SCHUR_PIVOT_EXACT + DUAL_EXCLUSION_FUNCTIONAL rev. 2) "
                "scored against the EFT-hedron gapped forward-limit bound at Hankel order t = 2",
        headline="degree-4 Hausdorff walls reproduced on two slices: interior PERMITTED, all four "
                 "exterior points REJECTED with registered negative pivots and certified Gram duals; "
                 "walls bracketed to <= 2^-20; S1 walls equal Gauss-Legendre / Simpson; t=1 is blind to all four",
        certificate_class="EXACT-RATIONAL (same class as B1/B2/B7/B15-POS); published-bound reproduction",
        results=results, soundness="SOUND", warnings=warnings,
        ground_truth_route="published-bound:arXiv:2012.15849 eqs.7.45/7.46 + KrN Thm III.2.3 "
                           "(Curto-Fialkow 1991 Rem. 4.4); quadrature (Gauss-Legendre, Simpson) for S1",
        evidence_level="E0", pir_level="L2", assumptions=[],
        falsifier_direction="G1-G3: a bracket excluding its wall, a REJECTED point without a certified "
                            "Gram dual, or a wrong interior/exterior verdict files the verifier REJECTED "
                            "for t=2 external bounds (asm:B15-POS2-unverified); G4/G5 stop certification",
        seed=SEED, generator_sha256=C.content_hash({"module": "b15_surf.pos.certify_t2", "registered": registered}),
        m_layer_stipulations=[
            "Bound is a theorem-level statement about EFT coefficients; no apparatus, no data",
            "Moments truncated at a_{6,0} (HEURISTIC truncation, see warnings)",
            "Units a_{2,0} = 1, M_Gap = 1; slices S1, S2 and all points registered in prereg-004",
            "Pre-freeze out-of-band check disclosed in prereg-004 (values unchanged)",
        ],
        calibration_route="none required (exact rational arithmetic; bracket width registered)",
        verdict="PERMITTED", witness=witness, impossibility_certificate=None,
        inputs={"registered": registered, "prereg": "docs/preregistrations/prereg-004-pos-hankel-t2.md"},
        pir_facts=[fact.to_dict()],
    )
    C.save(cert, CERT_NAME)
    r8 = t8_standalone(CERT)
    cert["results"]["T8_standalone_reverify"] = r8
    cert["results"]["T8_schema_negatives"] = t8_schema_negatives(cert)
    cert["certificate_id"] = f"b15-pos2-{C.content_hash(cert)}"
    C.save(cert, CERT_NAME)
    print(f"\nCertificate written: certificates/{CERT_NAME}")
    print("ALL B15-POS2 TESTS PASS")
