"""B15-GID tests — blind grammar identification from pole / zero / split fingerprints.

T0  regression: adding pole_set / zero_loci / split_structure to pir.fingerprints.INVARIANT_KEYS
    leaves every pre-B15 full-fingerprint hash unchanged (_project skips absent keys)
T1  fingerprint build from the label-stripped dataset only: exact pole sets (single-variable
    model fit, overdetermined UNIQUE), zero-locus set, split ranks (RANK_TEST), redundancy metric
T2  M1-style quotient: features are invariant under the registered relabelling group
    (cyclic + reflection) — the canonical form of a relabelled dataset is byte-identical
T3  KnownGrammarDB identification with two-hash semantics; candidate forest via pir.candidates
T4  delta recovery: exact delta on the shifted class; NONIDENTIFIABLE(delta-underdetermined) on NLSM
T5  held-out: fit n <= 6, apply to n = 8 (falsifier F4: separable at n<=6, not at n=8)
T6  pre-registered adjudication — Outcome A / B / AMBIGUOUS / NONIDENTIFIABLE per prereg-002 §GID,
    filed for the primary correlator AND the registered zeros-only correlator
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fractions import Fraction as F

import pir.fingerprints as FP
from b15_surf import certificate as C
from b15_surf.exact import fmt
from b15_surf.gid import fingerprint as G
from b15_surf.gid import identify as I
from b15_surf.zero import dataset as D
from b15_surf.zero.kinematics import relabel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = D.SEED
CLASS_MAP = {"TrPhi3": "C0", "NLSM": "C1", "TrPhi3_shifted": "C2"}       # hidden from the pass
TRUTH = {v: k for k, v in CLASS_MAP.items()}
THRESHOLD = F(1)                  # exact agreement (prereg-003 A1; §9-OI-2)
FIT_N = (5, 6)
HELD_OUT_N = 8


def t0_invariant_keys_regression():
    for k in ("pole_set", "zero_loci", "split_structure"):
        assert k in FP.INVARIANT_KEYS
    legacy = {"psd_signature": {"dim": 3, "rank": 3, "n_pos": 3, "n_zero": 0, "n_neg": 0},
              "rank_sequence": [1, 2, 3], "spectral_class": "stable", "representation": "chart-A"}
    h_now = FP.full_fingerprint(legacy)
    # recompute with the pre-B15 key set: identical because absent keys are skipped
    old_keys = tuple(k for k in FP.INVARIANT_KEYS if k not in ("pole_set", "zero_loci", "split_structure"))
    h_old = "fpf_" + FP.sha256_hex({"invariant": FP._project(legacy, old_keys)})[:16]
    assert h_now == h_old, (h_now, h_old)
    print(f"T0 INVARIANT_KEYS extension is additive: legacy full hash {h_now} unchanged")
    return {"status": "PASS", "legacy_full_hash": h_now}


def t1_fingerprints(stripped):
    feats = {}
    for cls in sorted(TRUTH):
        for n in (5, 6, 8):
            try:
                f = G.features_for(stripped, cls, n)
            except KeyError:
                continue
            assert not f["warnings"], f["warnings"]
            assert f["n_vanishing_loci"] == f["n_registered_loci"]
            assert all(r == 1 for _, r in f["split_structure"])
            assert all(r > 1 for r in f["split_control_ranks"].values())
            assert F(f["control_nonzero_rate"]) >= F(95, 100)
            feats[(cls, n)] = f
    summary = {f"{cls}:n{n}": {"poles": len(f["raw_pole_set"]), "vanishing_loci": f["n_vanishing_loci"],
                               "split_ranks": [r for _, r in f["split_structure"]],
                               "split_control_ranks": sorted(f["split_control_ranks"].values()),
                               "redundancy_metric": f["redundancy_metric"]}
               for (cls, n), f in feats.items()}
    print("T1 fingerprints (label-stripped): " + "; ".join(
        f"{k}: {v['poles']} poles, {v['vanishing_loci']} zeros, splits {v['split_ranks']}, "
        f"redundancy {v['redundancy_metric']}" for k, v in summary.items()))
    return feats, {"status": "PASS", "summary": summary,
                   "soundness": "pole/zero/rank detection SOUND (exact overdetermined fits, exact rank)"}


def t2_quotient(ds, feats):
    # relabel every point of one class at n=6 by a rotation + reflection and re-extract
    import copy
    n = 6
    perm = {v: ((-(v - 1) + 2) % n) + 1 for v in range(1, n + 1)}       # reflection composed with shift
    ds2 = copy.deepcopy(ds)
    for r in ds2["records"]:
        if r["n"] != n:
            continue
        pt = {tuple(int(x) for x in k.split(",")): F(v) for k, v in r["point"].items()}
        r["point"] = {f"{i},{j}": fmt(v) for (i, j), v in sorted(relabel(pt, n, perm).items())}
        if r.get("probe"):
            i, j = (int(x) for x in r["probe"]["chord"].split(","))
            from b15_surf.zero.kinematics import norm
            a, b = norm(n, perm[i], perm[j])
            r["probe"]["chord"] = f"{a},{b}"
        if r.get("triangle"):
            r["triangle"] = [perm[v] for v in r["triangle"]]
    # loci / splits registered tables are relabelled too (they are kinematic metadata)
    reg = ds2["registered"]
    from b15_surf.zero.kinematics import norm
    reg["loci"][str(n)] = {lid: [",".join(map(str, norm(n, perm[int(p.split(',')[0])], perm[int(p.split(',')[1])])))
                                 for p in pairs] for lid, pairs in reg["loci"][str(n)].items()}
    st2 = D.label_stripped(ds2, CLASS_MAP)
    same = {}
    for cls in sorted(TRUTH):
        f1 = feats[(cls, n)]
        f2 = G.features_for(st2, cls, n)
        h1 = I._hashable(f1)
        h2 = I._hashable(f2)
        assert h1 == h2, (cls, h1, h2)
        assert FP.full_fingerprint(h1) == FP.full_fingerprint(h2)
        assert FP.specific_fingerprint(h1) == FP.specific_fingerprint(h2)
        same[cls] = FP.full_fingerprint(h1)
    print(f"T2 dihedral quotient: relabelled dataset gives byte-identical canonical features at n=6 ({same})")
    return {"status": "PASS", "relabelling": "reflection ∘ rotation (D_6)", "full_hashes": same}


def t3_identify(stripped, feats, known):
    res = {}
    for (cls, n), f in feats.items():
        r = I.identify(f, known, stripped["registered"], THRESHOLD)
        res[(cls, n)] = r
    # two-hash semantics: TrPhi3 and TrPhi3_shifted share the FULL hash, differ in SPECIFIC
    for n in (6, 8):
        a, b = res[("C0", n)], res[("C2", n)]
        assert a["full"] == b["full"] and a["specific"] != b["specific"]
        assert res[("C1", n)]["full"] != a["full"]
    print("T3 DB identification: " + "; ".join(
        f"{cls}:n{n} -> {r['db_lookup'].get('grammar_id', r['db_lookup']['status'])} "
        f"(lattice {r['lattice']['verdict']} {r['lattice']['compatible']})" for (cls, n), r in res.items()))
    return res, {"status": "PASS",
                 "two_hash": "TrPhi3 / TrPhi3_shifted share full (delta-insensitive) hash, differ in specific",
                 "lookups": {f"{cls}:n{n}": r["db_lookup"] for (cls, n), r in res.items()}}


def t4_delta(feats):
    out = {}
    for (cls, n), f in feats.items():
        d = f["delta"]
        truth = TRUTH[cls]
        if truth == "TrPhi3_shifted":
            assert d["verdict"] == "FORCED" and F(d["value"]) == D.DELTA_SHIFTED, d
        elif truth == "TrPhi3":
            assert d["verdict"] == "FORCED" and F(d["value"]) == 0, d
        else:
            assert d["verdict"] == "NONIDENTIFIABLE" and d["cause"] == "delta-underdetermined", d
        out[f"{cls}:n{n}"] = d
    print(f"T4 delta recovery: shifted class -> delta = {fmt(D.DELTA_SHIFTED)} exactly; "
          f"NLSM -> NONIDENTIFIABLE(delta-underdetermined)")
    return {"status": "PASS", "per_class": out}


def t5_held_out(res):
    def separable(n):
        return all(r["verdict"] == "PERMITTED" and r["selected"] == TRUTH[cls]
                   for (cls, m), r in res.items() if m == n)
    fit = {n: separable(n) for n in FIT_N}
    held = separable(HELD_OUT_N)
    f4 = all(fit.values()) and not held
    assert not f4, "F4: fingerprint_v1 over-fit (separable at n<=6, not at n=8) -> REJECTED"
    print(f"T5 held-out: separable at n<=6: {fit}; at n=8: {held}; F4 {'TRIGGERED' if f4 else 'not triggered'}")
    return {"status": "PASS", "separable_fit": {str(k): v for k, v in fit.items()},
            "separable_held_out": held, "falsifier_F4": "not triggered"}


def t6_adjudicate(res):
    n = HELD_OUT_N
    rows = {cls: res[(cls, n)] for cls in sorted(TRUTH) if (cls, n) in res}
    # primary correlator
    all_identified = all(r["verdict"] == "PERMITTED" and r["selected"] == TRUTH[cls] for cls, r in rows.items())
    all_thresh = all(F(r["primary"]["similarity"]) >= THRESHOLD and F(r["primary"]["confidence"]) >= THRESHOLD
                     for r in rows.values())
    deltas_ok = res[("C2", n)]["delta"]["verdict"] == "FORCED"
    if all_identified and all_thresh and deltas_ok:
        primary = {"outcome": "A", "verdict": "PERMITTED",
                   "note": "classes distinguished (menu-relative identification) and delta recovered at n=8"}
    else:
        primary = {"outcome": "other", "verdict": "NONIDENTIFIABLE"}
    # secondary (zeros-only) correlator
    sec_ties = {cls: r["secondary"]["ties"] for cls, r in rows.items()}
    if all(len(t) >= 2 for t in sec_ties.values()):
        secondary = {"outcome": "B", "verdict": "OBSERVATIONALLY_EQUIVALENT",
                     "class": sorted(set().union(*sec_ties.values())),
                     "note": "fingerprint-only (zeros + splits) identification is structurally limited "
                             "for deformation-related grammars: the theories share zeros by construction"}
    else:
        secondary = {"outcome": "A", "verdict": "PERMITTED"}
    inv_ties = {cls: r["invariants_only"]["ties"] for cls, r in rows.items()}
    print(f"T6 adjudication at n={n}: primary -> Outcome {primary['outcome']} ({primary['verdict']}); "
          f"zeros-only -> Outcome {secondary['outcome']} ({secondary['verdict']} {secondary.get('class')}); "
          f"invariants-only ties {inv_ties}")
    return {"status": "PASS", "held_out_n": n, "threshold": fmt(THRESHOLD),
            "primary": {**primary, "correlator": "b15.gid.jaccard_atoms_v1",
                        "per_class": {cls: {"selected": r["selected"], "similarity": r["primary"]["similarity"],
                                            "confidence": r["primary"]["confidence"]} for cls, r in rows.items()}},
            "invariants_only": {"correlator": "b15.gid.invariants_only_v1", "ties": inv_ties,
                                "verdict": "AMBIGUOUS" if any(len(t) >= 2 for t in inv_ties.values()) else "PERMITTED",
                                "note": "delta-insensitive invariants cannot separate a base grammar from its "
                                        "deformation (by construction); NLSM still separated by its pole set"},
            "secondary": {**secondary, "correlator": "b15.gid.zeros_only_v1", "ties": sec_ties}}


if __name__ == "__main__":
    if not os.path.exists(D.DATASET_PATH):
        D.write()
    I.write_known_grammars()
    known = I.load_known_grammars()
    ds = D.load()
    stripped = D.label_stripped(ds, CLASS_MAP)

    r0 = t0_invariant_keys_regression()
    feats, r1 = t1_fingerprints(stripped)
    r2 = t2_quotient(ds, feats)
    res, r3 = t3_identify(stripped, feats, known)
    r4 = t4_delta(feats)
    r5 = t5_held_out(res)
    r6 = t6_adjudicate(res)
    results = {"T0_invariant_keys_regression": r0, "T1_fingerprint_build": r1, "T2_quotient": r2,
               "T3_db_identify": r3, "T4_delta_recovery": r4, "T5_held_out": r5, "T6_adjudication": r6}
    n = HELD_OUT_N
    verdict = r6["primary"]["verdict"]
    witness = {"selected_at_n8": {cls: res[(cls, n)]["selected"] for cls in sorted(TRUTH)},
               "delta_recovered": res[("C2", n)]["delta"],
               "predicates": {cls: res[(cls, n)]["predicates"] for cls in sorted(TRUTH)}}
    forest = [res[(cls, n)]["pir_fact"] for cls in sorted(TRUTH)]
    cert = C.build(
        benchmark="B15-GID",
        problem="B15-GID: blind identification of {TrPhi3, NLSM, TrPhi3_shifted} from pole / zero / "
                "split fingerprints of the label-stripped dataset (B8 -> KnownGrammarDB -> pir.candidates)",
        headline=f"primary correlator: Outcome {r6['primary']['outcome']} at held-out n=8 "
                 f"({verdict}); zeros-only correlator: Outcome {r6['secondary']['outcome']} "
                 f"({r6['secondary']['verdict']}) — zeros are shared by construction",
        certificate_class="EXACT-RATIONAL fingerprints; framework-capability benchmark",
        results=results, soundness="SOUND", warnings=[],
        ground_truth_route="planar-feynman-sum (dataset data/b15/zero_dataset_v1.json, label-stripped)",
        evidence_level="E0", pir_level="L3", assumptions=[],
        falsifier_direction="F4: separable at n<=6 but not at n=8 -> fingerprint_v1 REJECTED (over-fit); "
                            "a wrong selected grammar at any n -> REJECTED",
        seed=SEED, generator_sha256=ds["generator_sha256"],
        m_layer_stipulations=[
            "Input is the exact label-stripped dataset; class ids are opaque; registered delta withheld",
            "Grammar menu = {TrPhi3, NLSM, TrPhi3_shifted}: AMBIGUOUS/NONIDENTIFIABLE are menu-relative",
            "Similarity and confidence are separate exact ratios; correlator named per field",
            "Threshold = exact agreement (similarity = confidence = 1), prereg-003 A1 (§9-OI-2)",
            "Single surviving menu member -> PERMITTED (menu-relative), not FORCED: prereg-003 A2 (§9-OI-9)",
        ],
        calibration_route="threshold = exact agreement (prereg-003 A1); no data-fitted parameter — the pass has "
                          "no tunable numeric tolerance (all detection exact)",
        verdict=verdict, verdict_cause=None if verdict != "NONIDENTIFIABLE" else "threshold/lattice",
        verdict_class=None, witness=witness, impossibility_certificate=None,
        similarity=min((res[(cls, n)]["primary"]["similarity"] for cls in sorted(TRUTH)), key=F),
        confidence=min((res[(cls, n)]["primary"]["confidence"] for cls in sorted(TRUTH)), key=F),
        correlator="b15.gid.jaccard_atoms_v1",
        inputs={"dataset_path": "data/b15/zero_dataset_v1.json", "dataset_sha256": ds["dataset_sha256"],
                "known_grammars": "data/b15/known_grammars_v1.json", "threshold": fmt(THRESHOLD),
                "prereg_amendment": C.amendment_ref(),
                "fit_n": list(FIT_N), "held_out_n": HELD_OUT_N, "poly_degree": G.POLY_DEGREE},
        pir_facts=forest,
    )
    C.save(cert, "b15_gid_certificate.json")
    print("\nCertificate written: certificates/b15_gid_certificate.json")
    print("ALL B15-GID TESTS PASS")
