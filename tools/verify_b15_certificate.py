#!/usr/bin/env python3
"""Standalone re-verifier for B15 certificates (WP1.7 / hard constraint 5).

Recomputes, from the certificate JSON alone and with NO import of ``b15_surf``
(nor of any benchmark module): Schur pivot chains from the stored moments,
the dual exclusion functional's evaluation and its polynomial nonnegativity
identity, the Farkas companion (y^T A = 0, y^T b != 0), the bracket
containment of the published walls, and the content-addressed certificate id.
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


def verify_common(cert, fails):
    if cert["verdict"] not in SPEC:
        fails.append("non-SPEC verdict")
    if cert["soundness"] == "HEURISTIC" and cert["evidence_level"] == "E0":
        fails.append("HEURISTIC asserting E0")
    if cert["soundness"] == "HEURISTIC" and not cert["warnings"]:
        fails.append("HEURISTIC without warnings")
    if not cert["certificate_id"].endswith(content_hash(cert)):
        fails.append("certificate_id != content hash")
    if not cert["prereg_ref"].startswith("prereg-002@"):
        fails.append("prereg_ref missing")
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
