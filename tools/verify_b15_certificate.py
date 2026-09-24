#!/usr/bin/env python3
"""Standalone re-verifier for B15 certificates (WP1.7 / hard constraint 5).

Recomputes, from the certificate JSON alone and with NO import of ``b15_surf``
(nor of any benchmark module): Schur pivot chains from the stored moments,
the dual exclusion functional's evaluation and its polynomial nonnegativity
identity, the Farkas companion (y^T A = 0, y^T b != 0), the bracket
containment of the published walls, and the content-addressed certificate id.
For B15-POS2 (prereg-004, Hankel order t = 2) it recomputes the 3x3 Hankel and
2x2 localizing pivot chains, the Gram-form dual (DUAL_EXCLUSION_FUNCTIONAL rev. 2:
coefficient identity and Q0, Q1 >= 0), the brackets against the registered walls,
and the quadrature values of the S1 walls.
For B15-ZERO / B15-GID certificates it checks the schema-level invariants
(SPEC verdict, soundness/E-level honesty, content hash, dataset hash binding).

Exit status 0 = verified; 1 = FAIL (Score = infinity).
"""

import hashlib
import json
import sys
from fractions import Fraction as F

SPEC = {"FORCED", "PERMITTED", "REJECTED", "NONIDENTIFIABLE", "OBSERVATIONALLY_EQUIVALENT",
        "APPARATUS_LIMITED", "REPRESENTATION_DEPENDENT", "AMBIGUOUS"}


def canonical(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def content_hash(cert):
    body = {k: v for k, v in cert.items() if k not in ("certificate_id", "timestamp_utc")}
    return hashlib.sha256(canonical(body).encode("utf-8")).hexdigest()[:12]


def pivots(M):
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    out = []
    for k in range(n):
        d = A[k][k]
        if d < 0:
            out.append(d)
            return out, False
        if d == 0:
            if any(A[k][j] != 0 for j in range(k + 1, n)):
                out.append(d)
                return out, False
            out.append(F(0))
            continue
        out.append(d)
        for i in range(k + 1, n):
            f = A[i][k] / d
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return out, True


def chains_from_mu(mu):
    mu = [F(x) for x in mu]
    H = [[mu[0], mu[1]], [mu[1], mu[2]]]
    out = {"hankel": pivots(H)}
    for j in range(2):
        out[f"gap_{j}"] = pivots([[mu[j] - mu[j + 1]]])
    return out


def check_point(r, fails):
    mu = r["primal"]["mu"]
    ch = chains_from_mu(mu)
    feasible = all(ok for _, ok in ch.values())
    for name, (piv, _) in ch.items():
        stored = [F(x) for x in r["primal"]["chains"][name]["pivots"]]
        if stored != piv:
            fails.append(f"pivot mismatch {name}: {stored} vs {piv}")
    expected = "PERMITTED" if feasible else "REJECTED"
    if r["verdict"] != expected:
        fails.append(f"verdict {r['verdict']} but recomputed {expected}")
    if not feasible:
        ic = r.get("impossibility_certificate") or {}
        d = ic.get("dual_functional")
        if not d:
            fails.append("REJECTED without dual functional")
            return
        y = [F(x) for x in d["y"]]
        val = sum(a * b for a, b in zip(y, [F(x) for x in mu]))
        if val != F(d["evaluation"]) or val >= 0:
            fails.append(f"dual evaluation {val} inconsistent")
        w = d["nonnegativity_witness"]
        v0, v1 = F(w["v"][0]), F(w["v"][1])
        s0, s1 = F(w["s0"]), F(w["s1"])
        p = int(w.get("x_power", 0))
        poly = [F(0)] * 4
        sq = [v0 * v0, 2 * v0 * v1, v1 * v1]
        for k, c in enumerate(sq):
            poly[k] += s0 * c
        poly[p] += s1
        poly[p + 1] -= s1
        yy = y + [F(0)] * (4 - len(y))
        if s0 < 0 or s1 < 0 or poly != yy:
            fails.append("nonnegativity identity fails")
        fa = ic.get("farkas_atom")
        if fa:
            A = [[F(x) for x in row] for row in fa["A"]]
            b = [F(x) for x in fa["b"]]
            fy = [F(x) for x in fa["farkas_y"]]
            yA = [sum(fy[i] * A[i][j] for i in range(len(A))) for j in range(len(A[0]))]
            yb = sum(fy[i] * b[i] for i in range(len(A)))
            if any(t != 0 for t in yA) or yb == 0:
                fails.append("Farkas certificate fails")


def verify_pos(cert, fails):
    res = cert["results"]
    check_point(res["T1_primal_permitted"], fails)
    for k in ("hankel_wall", "gap_wall"):
        check_point(res["T2_dual_rejected"][k], fails)
    check_point(res["T4_negative_control"], fails)
    for k in ("lower_wall", "upper_wall"):
        br = res["T3_bracket"][k]
        a, b = (F(x) for x in br["certified_inner_interval"])
        wall = F(br["published_wall"])
        mu1 = F(br["mu1"])
        if not (a <= wall <= b):
            fails.append(f"bracket {k} excludes published wall")
        if (b - a) > F(res["T3_bracket"]["width_threshold"]):
            fails.append(f"bracket {k} too wide")
        # endpoints: inner feasible, outer infeasible, recomputed
        fin = all(ok for _, ok in chains_from_mu(["1", str(mu1), br["inner_feasible"]]).values())
        fout = all(ok for _, ok in chains_from_mu(["1", str(mu1), br["outer_infeasible"]]).values())
        if not fin or fout:
            fails.append(f"bracket {k} endpoints not certified")


def _poly_gram(Q0, Q1):
    """Coefficients (degree 4) of z2^T Q0 z2 + x(1-x) z1^T Q1 z1."""
    poly = [F(0)] * 5
    for i in range(3):
        for j in range(3):
            poly[i + j] += F(Q0[i][j])
    for i in range(2):
        for j in range(2):
            poly[i + j + 1] += F(Q1[i][j])
            poly[i + j + 2] -= F(Q1[i][j])
    return poly


def chains_t2(mu):
    mu = [F(x) for x in mu]
    H = [[mu[i + j] for j in range(3)] for i in range(3)]
    B = [[mu[i + j + 1] - mu[i + j + 2] for j in range(2)] for i in range(2)]
    return {"hankel": pivots(H), "gap": pivots(B)}


def check_point_t2(r, fails, tag):
    mu = r["primal"]["mu"]
    ch = chains_t2(mu)
    feasible = all(ok for _, ok in ch.values())
    for name, (piv, _) in ch.items():
        stored = [F(x) for x in r["primal"]["chains"][name]["pivots"]]
        if stored != piv:
            fails.append(f"{tag}: pivot mismatch {name}: {stored} vs {piv}")
    expected = "PERMITTED" if feasible else "REJECTED"
    if r["verdict"] != expected:
        fails.append(f"{tag}: verdict {r['verdict']} but recomputed {expected}")
    if feasible:
        return
    d = (r.get("impossibility_certificate") or {}).get("dual_functional")
    if not d:
        fails.append(f"{tag}: REJECTED without dual functional")
        return
    if d.get("op") != "DUAL_EXCLUSION_FUNCTIONAL" or d.get("op_rev") != 2:
        fails.append(f"{tag}: dual is not DUAL_EXCLUSION_FUNCTIONAL rev. 2")
    y = [F(x) for x in d["y"]]
    val = sum(a * b for a, b in zip(y, [F(x) for x in mu]))
    if val != F(d["evaluation"]) or val >= 0:
        fails.append(f"{tag}: dual evaluation {val} inconsistent")
    Q0, Q1 = d["gram"]["Q0"], d["gram"]["Q1"]
    if _poly_gram(Q0, Q1) != y:
        fails.append(f"{tag}: Gram identity fails")
    for nm, Q in (("Q0", Q0), ("Q1", Q1)):
        _, ok = pivots([[F(x) for x in row] for row in Q])
        if not ok:
            fails.append(f"{tag}: Gram {nm} not PSD")


def verify_pos2(cert, fails):
    res = cert["results"]
    reg = cert["inputs"]["registered"]
    for name, r in res["T1_primal"]["points"].items():
        check_point_t2(r, fails, name)
        if r["verdict"] != "PERMITTED":
            fails.append(f"{name}: interior point not PERMITTED")
    for name, r in res["T2_dual"]["points"].items():
        check_point_t2(r, fails, name)
        if r["verdict"] != "REJECTED":
            fails.append(f"{name}: exterior point not REJECTED")
    check_point_t2(res["T6_negative_control"], fails, "negative_control")
    width = F(res["T3_brackets"]["width_threshold"])
    if width != F(reg["width"]):
        fails.append("bracket width threshold != registered")
    for s, br2 in res["T3_brackets"]["slices"].items():
        rs = reg["slices"][s]
        for k, key in (("lower_wall", "L"), ("upper_wall", "U")):
            br = br2[k]
            wall = F(rs[key])
            if F(br["registered_wall"]) != wall:
                fails.append(f"{s} {k}: stored wall != registered")
            a, b = (F(x) for x in br["certified_inner_interval"])
            if not (a <= wall <= b):
                fails.append(f"{s} {k}: bracket excludes registered wall")
            if (b - a) > width:
                fails.append(f"{s} {k}: bracket too wide")
            sl = rs["slice"]
            fin = all(ok for _, ok in chains_t2(["1", *sl, br["inner_feasible"]]).values())
            fout = all(ok for _, ok in chains_t2(["1", *sl, br["outer_infeasible"]]).values())
            if not fin or fout:
                fails.append(f"{s} {k}: bracket endpoints not certified")
    # quadrature route for S1 (Simpson rational; Gauss via x^4 mod x^2 - x + 1/6)
    simpson = F(1, 6) * 0 + F(2, 3) * F(1, 16) + F(1, 6) * 1
    a, b = F(1), F(0)
    for _ in range(3):
        a, b = a + b, -a / 6
    gauss = a / 2 + b
    q = res["T4_quadrature"]
    if F(q["gauss_legendre_2pt_mu4"]) != gauss or F(q["simpson_mu4"]) != simpson:
        fails.append("quadrature values do not recompute")
    if gauss != F(reg["slices"]["S1"]["L"]) or simpson != F(reg["slices"]["S1"]["U"]):
        fails.append("quadrature route disagrees with registered S1 walls")


def verify_common(cert, fails):
    if cert["verdict"] not in SPEC:
        fails.append("non-SPEC verdict")
    if cert["soundness"] == "HEURISTIC" and cert["evidence_level"] == "E0":
        fails.append("HEURISTIC asserting E0")
    if cert["soundness"] == "HEURISTIC" and not cert["warnings"]:
        fails.append("HEURISTIC without warnings")
    if not cert["certificate_id"].endswith(content_hash(cert)):
        fails.append("certificate_id != content hash")
    want = "prereg-004@" if cert.get("benchmark") == "B15-POS2" else "prereg-002@"
    if not cert["prereg_ref"].startswith(want) or "UNFROZEN" in cert["prereg_ref"]:
        fails.append(f"prereg_ref missing or not {want}<commit>")
    for a in cert["assumptions"]:
        if not a.startswith("asm:"):
            fails.append(f"bad assumption {a}")


def verify_zero(cert, fails):
    res = cert["results"]
    for k, v in res.items():
        if isinstance(v, dict) and v.get("status") not in ("PASS", None):
            fails.append(f"{k}: {v.get('status')}")
    ds = cert["inputs"].get("dataset_sha256")
    if ds and cert["inputs"].get("dataset_path"):
        try:
            with open(cert["inputs"]["dataset_path"]) as f:
                data = json.load(f)
            if data.get("dataset_sha256") != ds:
                fails.append("dataset hash does not match certificate")
        except OSError:
            fails.append("dataset file missing")


def main(argv):
    if len(argv) != 2:
        print("usage: verify_b15_certificate.py <certificate.json>")
        return 2
    with open(argv[1]) as f:
        cert = json.load(f)
    fails = []
    verify_common(cert, fails)
    b = cert.get("benchmark")
    if b == "B15-POS":
        verify_pos(cert, fails)
    elif b == "B15-POS2":
        verify_pos2(cert, fails)
    elif b in ("B15-ZERO", "B15-GID"):
        verify_zero(cert, fails)
    else:
        fails.append("unknown benchmark")
    if fails:
        print("FAIL (Score = infinity):")
        for x in fails:
            print("  - " + x)
        return 1
    print(f"VERIFIED {cert['certificate_id']} ({b}, verdict {cert['verdict_display']})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
