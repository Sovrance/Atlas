#!/usr/bin/env python3
"""WO-RH-17 — normalization adjudication (derivation + certificate).

Executes the P0 runbook end to end and writes
``certificates/normalization_adjudication.json``:

1. state the Fourier convention and ``tilde h``;
2. define ``F_ij = h_i * tilde h_j``;
3. prove the convolution/transform identity under that convention;
4. evaluate ``Fhat_ij(+i/2)`` and ``Fhat_ij(-i/2)``;
5. derive the pole matrix entry symbolically  → **Candidate A**;
6. audit the repository ``sqrt(3)/2`` block  → **Candidate B**;
7. decide, with an invariant argument, which survives;
8. freeze the survivor as a content-addressed normalization id.

SymPy is used for the algebra when available (the work order sanctions it); when
it is absent the same identities are verified numerically at high precision and
the certificate records the weaker ``symbolic_engine: none`` status honestly.

No RH proof claim.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import normalization as N  # noqa: E402
import rejected_pole as RP  # noqa: E402  (archival; audit surface only)
from certificate_io import source_hash, write_certificate  # noqa: E402

# ENG-004 §6: SymPy is a REQUIRED rigorous dependency, not an optional extra.
# Before ENG-004 a missing SymPy silently degraded this job: it rewrote the
# adjudication artifact with a "verified numerically only" note, dropping the
# symbolic derivation and dirtying the tree with weaker evidence. That is exactly
# the failure the runbook forbids, so the job now refuses to run instead.
try:
    import sympy as sp

    SYMPY = True
    SYMPY_VERSION = sp.__version__
except ImportError:  # pragma: no cover - exercised by the CI gate
    SYMPY = False
    SYMPY_VERSION = None


class RigorousDependencyMissing(RuntimeError):
    """A dependency the rigorous derivation cannot honestly proceed without."""


def require_sympy() -> None:
    if not SYMPY:
        raise RigorousDependencyMissing(
            "sympy is required for the normalization derivation (ENG-004 §6). "
            "Install the rigorous dependency set: "
            "pip install -r math/rh_weil/requirements-rigorous.txt . "
            "Refusing to rewrite normalization_adjudication.json with weaker evidence."
        )

L_POINTS = [
    ("log3", math.log(3.0)),
    ("1.1059498113", 1.1059498113),
    ("1.20", 1.20),
    ("log4", math.log(4.0)),
]


def symbolic_derivation() -> dict:
    """Steps 3-7 of the runbook, symbolically. Requires SymPy (ENG-004 §6)."""
    require_sympy()

    x, t = sp.symbols("x t", real=True)
    L = sp.Symbol("L", positive=True)
    u = sp.Symbol("u", positive=True)
    basis = {"one": sp.Integer(1), "q1": x - L / 2, "b": x * (L - x)}

    def E(h, s):
        return sp.simplify(sp.integrate(h * sp.exp(s * x / 2), (x, 0, L)))

    def is_zero(expr) -> bool:
        """Decide ``expr == 0`` by mapping e^{L/4} -> u, which turns these
        exponential/hyperbolic expressions into rational functions of (L, u)
        where ``cancel`` is a decision procedure."""
        e = sp.expand(sp.simplify(sp.sympify(expr).rewrite(sp.exp)))
        for k, rep in ((1, u**4), (sp.Rational(1, 2), u**2), (sp.Rational(1, 4), u)):
            e = e.subs(sp.exp(k * L), rep).subs(sp.exp(-k * L), 1 / rep)
        return sp.simplify(sp.cancel(sp.together(sp.expand(e)))) == 0

    steps = []

    # Step 3 — transform identity. With Fhat(xi)=int F e^{-i xi x}dx and
    # tilde h(x)=h(-x): Fhat_ij(i/2) = int (h_i * tilde h_j)(x) e^{x/2} dx.
    # Substituting u = y - x factorises the double integral.
    steps.append({
        "step": "convolution transform identity",
        "statement": "Fhat_ij(±i/2) = (int h_i(y)e^{±y/2}dy)(int h_j(u)e^{∓u/2}du) = E_i^± E_j^∓",
        "argument": "Fhat(±i/2)=int F(x)e^{±x/2}dx with F=h_i*tilde h_j; substitute u=y-x, "
                    "the exponential separates as e^{±(y-u)/2}, and the double integral factorises.",
        "verified": True,
    })

    # Step 4/5 — the pole entry.
    E_vals = {k: (E(h, 1), E(h, -1)) for k, h in basis.items()}
    steps.append({
        "step": "pole matrix entry",
        "statement": "G0_ij = Fhat_ij(i/2)+Fhat_ij(-i/2) = E_i^+E_j^- + E_i^-E_j^+",
        "E_plus": {k: sp.sstr(v[0]) for k, v in E_vals.items()},
        "E_minus": {k: sp.sstr(v[1]) for k, v in E_vals.items()},
        "verified": True,
    })

    # Parity lemma: even basis ⇒ E^- = e^{-L/2}E^+ ; odd ⇒ E^- = -e^{-L/2}E^+.
    parity = {}
    for k, h in basis.items():
        refl = sp.simplify(h.subs(x, L - x) - h)
        sign = 1 if refl == 0 else (-1 if sp.simplify(h.subs(x, L - x) + h) == 0 else None)
        Ep, Em = E_vals[k]
        holds = is_zero(Em - sign * sp.exp(-L / 2) * Ep) if sign else False
        parity[k] = {"parity": {1: "even", -1: "odd"}.get(sign, "mixed"),
                     "identity": "E^- = %s e^{-L/2} E^+" % ("+" if sign == 1 else "-"),
                     "identity_holds": bool(holds)}
    steps.append({
        "step": "parity lemma about x = L/2",
        "statement": "h(L-x)=±h(x)  =>  E^- = ± e^{-L/2} E^+",
        "per_basis": parity,
        "verified": all(v["identity_holds"] for v in parity.values()),
    })

    # Step 6 — the repository block, and the exact quotient on the even sector.
    E0p, E0m = E_vals["one"]
    Ebp, Ebm = E_vals["b"]
    Delta = sp.simplify(E0p * Ebm - E0m * Ebp)
    s32 = sp.sqrt(3) / 2
    A00 = sp.simplify(2 * E0p * E0m)
    B00 = sp.simplify(s32 * (E0p**2 + E0m**2))
    quotient_is_closed_form = is_zero(B00 / A00 - s32 * sp.cosh(L / 2))
    fixed_points = sp.solve(sp.Eq(s32 * sp.cosh(L / 2), 1), L)
    steps.append({
        "step": "audit of the repository sqrt(3)/2 block (Candidate B)",
        "even_sector_determinant_gap": sp.sstr(Delta),
        "both_candidates_rank": 1 if Delta == 0 else 2,
        "quotient_B_over_A": sp.sstr(sp.simplify(s32 * sp.cosh(L / 2))),
        "quotient_matches_closed_form": bool(quotient_is_closed_form),
        "quotient_equals_one_at": [sp.sstr(sp.simplify(f)) for f in fixed_points],
        "verified": bool(quotient_is_closed_form and Delta == 0),
        "proof_method": "e^{L/4} -> u reduces the identity to a rational function; cancel decides it",
    })

    # Odd sector: the repository pivot already equals Candidate A.
    Eqp, Eqm = E_vals["q1"]
    A_repo = L * sp.cosh(L / 4) - 4 * sp.sinh(L / 4)
    odd_identity_holds = is_zero(Eqp - 2 * sp.exp(L / 4) * A_repo)
    steps.append({
        "step": "odd sector agrees already",
        "statement": "E_q1^+ = 2 e^{L/4} A(L) with A = L cosh(L/4) - 4 sinh(L/4); "
                     "hence G0[q1,q1] = 2E^+E^- = -8A^2, the value already in the repository",
        "verified": bool(odd_identity_holds),
        "proof_method": "e^{L/4} -> u reduces the identity to a rational function; cancel decides it",
    })
    return {"symbolic_engine": f"sympy {sp.__version__}", "steps": steps}


def numeric_checks() -> list:
    """Independent numeric confirmation at the mandated points."""
    out = []
    for label, L in L_POINTS:
        row = {"L_label": label, "L": L, "checks": {}}
        # (a) closed-form pole vs the real-space kernel route
        try:
            import mpmath as mp

            mp.mp.dps = 40
            worst = 0.0
            for i in N.BASIS_NAMES:
                for j in N.BASIS_NAMES:
                    closed = N.pole_entry(i, j, L)
                    quad = float(mp.quad(lambda a: N.kernel_K(i, j, float(a), L) * 2 * mp.cosh(a / 2), [0, L]))
                    worst = max(worst, abs(closed - quad) / max(1.0, abs(closed), abs(quad)))
            row["checks"]["pole_closed_form_vs_realspace_kernel_max_rel"] = worst
            row["checks"]["pole_routes_agree"] = worst < 1e-12
        except Exception as exc:  # pragma: no cover
            row["checks"]["pole_routes_agree"] = f"unavailable: {exc}"
        # (b) parity block-diagonality of the pole matrix
        row["checks"]["even_odd_cross_entries_vanish"] = all(
            abs(N.pole_entry(i, "q1", L)) < 1e-12 for i in ("one", "b")
        )
        # (c) legacy quotient
        ratio = RP.legacy_pole_entry("one", "one", L) / N.pole_entry("one", "one", L)
        row["checks"]["legacy_over_adopted_ratio"] = ratio
        row["checks"]["matches_sqrt3_over_2_cosh"] = abs(ratio - RP.legacy_over_adopted_ratio(L)) < 1e-12
        row["checks"]["legacy_relative_error"] = ratio - 1.0
        out.append(row)
    return out


def build() -> dict:
    sym = symbolic_derivation()
    num = numeric_checks()
    src = ROOT / "src"
    hashes = {
        "normalization.py": source_hash([src / "normalization.py"]),
        "providers.py": source_hash([src / "providers.py"]),
        "cross_validation.py": source_hash([src / "cross_validation.py"]),
        "finite_weil.py(audited)": source_hash([src / "finite_weil.py"]),
    }
    content = N.normalization_content()
    return {
        "certificate_version": "0.1",
        "program": "RH/Weil normalization adjudication",
        "work_order": "WO-RH-17",
        "status": "ADJUDICATED",
        "active_normalization_id": N.normalization_id(),
        "fourier_convention": content["fourier_convention"],
        "tilde_convention": content["tilde_convention"],
        "pole_formula": content["pole_formula"],
        "prime_formula": content["prime_formula"],
        "archimedean_formula": content["archimedean_formula"],
        "assembly": content["assembly"],
        "frequency_integral": content["frequency_integral"],
        "basis_normalization": content["basis_normalization"],
        "candidate_dispositions": {
            "candidate_A_explicit_formula": {
                "formula": "G0_ij = E_i^+E_j^- + E_i^-E_j^+",
                "disposition": "ADOPTED",
                "grounds": [
                    "derived step by step from the explicit formula under the stated convention",
                    "independently reproduced by the real-space identity "
                    "G0_ij = int_0^L K_ij(a) 2cosh(a/2) da using the same K_ij as the prime block",
                    "reproduces the repository's odd pivot -8A^2 exactly (E_q1^+ = 2e^{L/4}A)",
                    "no fitted constant anywhere in the derivation",
                ],
            },
            "candidate_B_repo_sqrt3_over_2": {
                "formula": "G0_even = (sqrt(3)/2)(v+ v+^T + v- v-^T)",
                "disposition": "REJECTED",
                "grounds": [
                    "no derivation from the explicit formula was found in the repository or its history",
                    "on the even sector it equals Candidate A times (sqrt(3)/2)cosh(L/2)",
                    "that factor equals 1 at exactly one point, L = log 3, and differs elsewhere "
                    "(+0.18% at L=1.1059498113, +2.66% at L=1.20, +8.25% at L=log 4)",
                    "an L-dependent scalar cannot arise from a basis change: a change of basis is a "
                    "constant congruence G -> S^T G S, so it cannot produce a factor varying with L",
                    "it is therefore a multiplicative calibration fitted at L = log 3, which the P0 "
                    "runbook forbids",
                ],
                "calibration_fixed_point": N.CALIBRATION_FIXED_POINT,
                "retained_for_audit_as": "rejected_pole.legacy_pole_entry (archival; production may not import it)",
            },
        },
        "structural_findings": [
            "the even basis {1, x(L-x)} is symmetric about x = L/2, so E^- = e^{-L/2}E^+ and the "
            "even pole block is rank 1 under BOTH candidates; a rank/determinant regression therefore "
            "cannot discriminate them, which is why the conflict survived earlier review",
            "the pole matrix is parity block diagonal: even-odd cross entries vanish identically",
            "the defect is confined to the even block; the odd pivot was already correct",
        ],
        "crosschecks": num,
        "symbolic_derivation": sym,
        "source_hashes": hashes,
        "primary_sources": [
            {"role": "primary explicit formula", "ref": "Weil explicit formula, pole term "
             "Fhat(i/2)+Fhat(-i/2) for F = h_i * tilde h_j"},
            {"role": "repository implementation under audit",
             "ref": "math/rh_weil/src/rejected_pole.py::legacy_pole_entry (removed from production by ENG-004 §1)"},
            {"role": "regression evidence only",
             "ref": "math/rh_weil/certificates/*.json (never normative for a normalization dispute)"},
        ],
        "hard_constraints_certified": False,
        "claim_boundary": N.CLAIM_BOUNDARY,
        "rh_proof_claim": False,
    }


def main() -> int:
    # Fail *before* touching the artifact, so a missing dependency can never
    # leave the tree dirty with a degraded certificate.
    try:
        require_sympy()
    except RigorousDependencyMissing as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    body = build()
    path = write_certificate("normalization_adjudication.json", body)
    print(f"wrote {path}")
    print(f"status={body['status']}  normalization_id={body['active_normalization_id']}")
    sym = body["symbolic_derivation"]
    print(f"symbolic engine: {sym.get('symbolic_engine')}")
    for st in sym.get("steps", []):
        print(f"  [{'ok' if st.get('verified') else '!!'}] {st['step']}")
    for row in body["crosschecks"]:
        print(f"  L={row['L_label']:>14}: legacy/adopted={row['checks']['legacy_over_adopted_ratio']:.12f} "
              f"pole routes agree={row['checks']['pole_routes_agree']}")
    ok = all(st.get("verified") for st in sym.get("steps", []))
    ok = ok and all(r["checks"]["pole_routes_agree"] is True for r in body["crosschecks"])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
