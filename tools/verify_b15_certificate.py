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
For B15-POS3 (prereg-005, variance-zero slices) it additionally reproduces each
zero-pivot dual direction v from its recorded lift (stop index/kind, t, zero-pivot
rows, solved indices), re-checks the single-atom Farkas / consistency system, and
the degenerate certified_point (a^4 permitted, a^4 +- resolution rejected).
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


def _solve_small(A, b):
    """Gauss-Jordan over Q for a nonsingular square system; None if singular."""
    n = len(A)
    M = [[F(x) for x in row] + [F(b[i])] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def _reduced_at_stop(M):
    """Symmetric elimination; returns (pivots, stop index or None, reduced trailing block)."""
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    piv = []
    for k in range(n):
        d = A[k][k]
        piv.append(d)
        if d < 0 or (d == 0 and any(A[k][j] != 0 for j in range(k + 1, n))):
            return piv, k, [row[k:] for row in A[k:]]
        if d == 0:
            continue
        for i in range(k + 1, n):
            f = A[i][k] / d
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return piv, None, None


def check_zero_pivot_dual(mu, d, fails, tag):
    """Reproduce v from the recorded lift (prereg-005, P5-OI-2 (i))."""
    mu = [F(x) for x in mu]
    M = ([[mu[i + j] for j in range(3)] for i in range(3)] if d["matrix"] == "hankel"
         else [[mu[i + j + 1] - mu[i + j + 2] for j in range(2)] for i in range(2)])
    piv, k, S = _reduced_at_stop(M)
    lift = d["lift"]
    if k is None or lift["stop_index"] != k:
        fails.append(f"{tag}: recorded stop index does not recompute")
        return
    Z = [i for i in range(k) if piv[i] == 0]
    I = [i for i in range(k) if piv[i] != 0]
    if lift["zero_pivot_rows"] != Z or lift["solved_indices"] != I:
        fails.append(f"{tag}: recorded zero-pivot rows / solved indices do not recompute")
    w = [F(0)] * len(S)
    if piv[k] < 0:
        kind = "negative_pivot"
        w[0] = F(1)
    else:
        kind = "zero_pivot_nonzero_row"
        jj = next(c for c in range(1, len(S)) if S[0][c] != 0)
        t = -(S[jj][jj] + 1) / (2 * S[0][jj])
        if lift["j"] != k + jj or F(lift["t"]) != t:
            fails.append(f"{tag}: recorded j / t do not recompute")
        w[0], w[jj] = t, F(1)
    if lift["stop_kind"] != kind or [F(x) for x in lift["w"]] != w:
        fails.append(f"{tag}: recorded stop kind / w do not recompute")
    n = len(M)
    v = [F(0)] * n
    for ci in range(len(S)):
        v[k + ci] = w[ci]
    if I:
        sol = _solve_small([[M[a][b] for b in I] for a in I],
                           [-sum(M[a][c] * v[c] for c in range(k, n)) for a in I])
        if sol is None:
            fails.append(f"{tag}: lift system singular")
            return
        for m, i in enumerate(I):
            v[i] = sol[m]
    if [F(x) for x in d["direction"]] != v:
        fails.append(f"{tag}: direction v not reproduced from the recorded lift")
    val = sum(v[i] * M[i][j] * v[j] for i in range(n) for j in range(n))
    if val != F(d["evaluation"]) or val >= 0:
        fails.append(f"{tag}: v^T M v != recorded evaluation")
    if kind == "zero_pivot_nonzero_row" and val != -1:
        fails.append(f"{tag}: zero-pivot normalisation w^T S w = -1 violated")


def check_single_atom(sys_, mu, fails, tag, expect_consistent):
    a = F(sys_["a"])
    A = [[a ** i] for i in range(5)]
    b = [F(x) for x in mu]
    if [[F(x) for x in r] for r in sys_["A"]] != A or [F(x) for x in sys_["b"]] != b:
        fails.append(f"{tag}: single-atom system does not match the point")
        return
    if F(mu[2]) != F(mu[1]) ** 2 or not sys_.get("premise"):
        fails.append(f"{tag}: variance-zero premise missing or false")
    if expect_consistent:
        wv = [F(x) for x in sys_.get("w", [])]
        if sys_["status"] != "UNIQUE" or wv != [F(1)] or any(A[i][0] * wv[0] != b[i] for i in range(5)):
            fails.append(f"{tag}: single-atom system not UNIQUE with w = 1")
    else:
        y = [F(x) for x in sys_.get("farkas_y", [])]
        if sys_["status"] != "INCONSISTENT" or len(y) != 5 \
                or sum(y[i] * A[i][0] for i in range(5)) != 0 or sum(y[i] * b[i] for i in range(5)) == 0:
            fails.append(f"{tag}: Farkas vector fails (y^T A = 0, y^T b != 0)")


def _has_key(obj, key):
    if isinstance(obj, dict):
        return key in obj or any(_has_key(v, key) for v in obj.values())
    if isinstance(obj, list):
        return any(_has_key(v, key) for v in obj)
    return False


def verify_pos3(cert, fails):
    res = cert["results"]
    reg = cert["inputs"]["registered"]
    if _has_key(res, "certified_inner_interval"):
        fails.append("certified_inner_interval present (feasible set is a point; prereg-005 U3)")
    for name, r in {**res["U1_primal"]["points"], **res["U2_dual"]["points"]}.items():
        spec = reg["points"][name]
        mu = r["primal"]["mu"]
        if [F(x) for x in mu] != [F(1), F(spec["a"]), F(spec["a"]) ** 2, F(spec["mu3"]), F(spec["mu4"])]:
            fails.append(f"{name}: stored moments != registered point")
        check_point_t2(r, fails, name)
        if r["verdict"] != spec["expected"]:
            fails.append(f"{name}: verdict {r['verdict']} != registered {spec['expected']}")
        if r["verdict"] == "PERMITTED":
            check_single_atom(r["single_atom_system"], mu, fails, name, True)
        else:
            ic = r["impossibility_certificate"]
            d = ic["dual_functional"]
            check_zero_pivot_dual(mu, d, fails, name)
            if [F(x) for x in d["direction"]] != [F(x) for x in spec["v"]] \
                    or F(d["evaluation"]) != F(spec["value"]) or d["matrix"] != spec["matrix"]:
                fails.append(f"{name}: dual direction / value / matrix != registered")
            check_single_atom(ic["farkas_atom"], mu, fails, name, False)
    cp = res["U3_certified_point"]
    if F(cp["resolution"]) != F(reg["resolution"]):
        fails.append("certified_point resolution != registered")
    for s, e in cp["slices"].items():
        a4 = F(reg["slices"][s]["a"]) ** 4
        sl = ["1", reg["slices"][s]["a"], str(F(reg["slices"][s]["a"]) ** 2), str(F(reg["slices"][s]["a"]) ** 3)]
        if F(e["permitted"]) != a4 or not all(ok for _, ok in chains_t2([*sl, e["permitted"]]).values()):
            fails.append(f"{s}: certified point a^4 not recomputed as PERMITTED")
        for nb in e["rejected_neighbours"]:
            ch = chains_t2([*sl, nb["mu4"]])
            if abs(F(nb["mu4"]) - a4) != F(cp["resolution"]) or ch[nb["violated_wall"]][1] \
                    or ch["gap" if nb["violated_wall"] == "hankel" else "hankel"][1] is False:
                fails.append(f"{s}: neighbour {nb['mu4']} not rejected by exactly {nb['violated_wall']}")
            check_point_t2(nb["certificate"], fails, f"{s} neighbour {nb['mu4']}")
            check_zero_pivot_dual(nb["certificate"]["primal"]["mu"],
                                  nb["certificate"]["impossibility_certificate"]["dual_functional"],
                                  fails, f"{s} neighbour {nb['mu4']}")


def verify_common(cert, fails):
    if cert["verdict"] not in SPEC:
        fails.append("non-SPEC verdict")
    if cert["soundness"] == "HEURISTIC" and cert["evidence_level"] == "E0":
        fails.append("HEURISTIC asserting E0")
    if cert["soundness"] == "HEURISTIC" and not cert["warnings"]:
        fails.append("HEURISTIC without warnings")
    if not cert["certificate_id"].endswith(content_hash(cert)):
        fails.append("certificate_id != content hash")
    want = {"B15-POS2": "prereg-004@", "B15-POS3": "prereg-005@"}.get(cert.get("benchmark"), "prereg-002@")
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
    elif b == "B15-POS3":
        verify_pos3(cert, fails)
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
