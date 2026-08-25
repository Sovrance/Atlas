"""ATLAS-RH-ENG-007 §12 (WO-RH-44) -- formal evidence and its exact boundary.

The distinction this module exists to keep visible:

    Arb interval enclosure               numeric warrant      E1
    Lean: positive lower bounds => PD    implication warrant  FORMAL

A formal theorem strengthens an *exact theorem dependency*. It never converts
interval numerical evidence into FORMAL evidence. The degree-3 block is
positive definite because Arb certified two positive lower bounds **and**
because a proved theorem says those bounds suffice; remove either half and the
claim is gone. So a certificate that gains formal backing keeps its numeric
warrant exactly where it was and gains a second, separate field saying which
implication is now machine-checked.

Nothing here makes an RH proof claim, and no theorem in the manifest is about
the Riemann hypothesis: they are all finite linear algebra over the reals.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "formal" / "manifests" / "theorem_manifest.json"

#: PIR content kind introduced by ENG-007 §12.
KIND_FORMAL = "FORMAL_THEOREM_CERTIFICATE"

#: The warrant a formal theorem carries. Deliberately not a member of the
#: E0/E1/E2/E3 ladder -- it answers a different question. The ladder grades how
#: reliable a *number* is; this grades whether an *implication* was checked.
FORMAL_WARRANT = "FORMAL"

#: Which proved theorems back which numeric certificate. Each entry says: the
#: numbers in this file are E1, and these theorems are why those numbers imply
#: what the file claims. An entry is only justified when the theorem's
#: hypotheses are literally what the certificate reports.
FORMAL_BACKING: Dict[str, Tuple[str, ...]] = {
    # "E2(L) = G00 Gbb - G0b^2 >= bound > 0": the 2x2 criterion, applied to a
    # block whose parity structure and determinant factorization the
    # certificate's own parity_identities re-check numerically.
    "e1_degree2_compact_log3_log4.json": (
        "pd_two_by_two",
        "certificate_even2_implies_pd",
        "weil_basis_parity",
        "odd_degree3_factorization",
    ),
    # The odd degree-3 block, certified positive definite by interval LDL
    # congruence. The congruence step is what the inertia theorems license; the
    # two positive bounds are what the 2x2 criterion consumes; the block exists
    # as a block because of parity.
    "e1_degree3_odd_positivity_log3_log4.json": (
        "inertia_congruence_positive",
        "inertia_congruence_negative",
        "inertia_congruence_rank",
        "pd_two_by_two",
        "certificate_even2_implies_pd",
        "weil_basis_parity",
        "odd_degree3_cross_block",
        "odd_degree3_factorization",
    ),
    # The stratification variant, if the certification ever lands there instead.
    "e1_degree3_odd_inertia_log3_log4.json": (
        "inertia_congruence_positive",
        "inertia_congruence_negative",
        "inertia_congruence_rank",
        "weil_basis_parity",
        "odd_degree3_cross_block",
    ),
    # ENG-008: the 3x3 even block. Its certified numbers are bounds on the
    # leading minors of the *preconditioned* matrix, so the implication it needs
    # is the composed one -- rescale, read the minors, conclude about the
    # original block -- and that is `preconditioned_certificate3`.
    "e1_degree4_even3_positivity_log3_log4.json": (
        "pd_three_by_three",
        "pd_three_by_three_certificate",
        "preconditioned_certificate3",
        "diagonal_congruence_preserves_pd",
        "weil_basis_parity",
        "odd_degree3_cross_block",
    ),
    # The inertia artifact reads a signature rather than a yes/no, so what
    # licenses it is index and rank invariance under the preconditioner, not
    # definiteness.
    "e1_degree4_even3_inertia_log3_log4.json": (
        "inertia_congruence_positive",
        "inertia_congruence_negative",
        "inertia_congruence_rank",
        "diagonal_congruence_preserves_index",
        "diagonal_congruence_preserves_rank",
        "weil_basis_parity",
    ),
    "e1_degree4_even3_moments_log3_log4.json": ("rank_trace_hs",),
    # Rank-trace: the runtime uses the Q = 0, b = 0 case, which is the case
    # that is proved. The general case is recorded in the manifest as unproved
    # and carries no warrant.
    "e1_degree3_odd_moments_log3_log4.json": ("rank_trace_hs",),
}


# --------------------------------------------------------------------------- #
# manifest access                                                              #
# --------------------------------------------------------------------------- #
def manifest_available() -> bool:
    return MANIFEST_PATH.exists()


def load_manifest() -> Optional[Dict[str, Any]]:
    if not MANIFEST_PATH.exists():
        return None
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def manifest_id() -> Optional[str]:
    m = load_manifest()
    return m.get("manifest_id") if m else None


def proved_theorem_ids() -> Tuple[str, ...]:
    m = load_manifest()
    if not m:
        return ()
    return tuple(t["id"] for t in m.get("theorems") or ())


def unproved_statement_ids() -> Tuple[str, ...]:
    m = load_manifest()
    if not m:
        return ()
    return tuple(u["id"] for u in m.get("unproved_statements") or ())


def hashed_sources() -> Dict[str, str]:
    """The Lean sources the manifest binds, for a certificate's dependency block.

    Recording these on the *numeric* certificate would be wrong -- the Arb run
    does not depend on Lean. They belong on the formal certificate, which does.
    """
    m = load_manifest()
    return dict(m.get("sources") or {}) if m else {}


# --------------------------------------------------------------------------- #
# the annotation a numeric certificate gains                                   #
# --------------------------------------------------------------------------- #
def formal_block(certificate_file: str, numeric_warrant: Optional[str]) -> Optional[Dict[str, Any]]:
    """The §12 fields for a numeric certificate, or ``None`` if it has no backing.

    ``numeric_warrant`` is passed in rather than inferred so that this function
    can never quietly upgrade it. It comes back out unchanged; the only new
    information is which implication is formal.
    """
    ids = FORMAL_BACKING.get(certificate_file)
    if not ids:
        return None
    mid = manifest_id()
    if mid is None:
        return None
    proved = set(proved_theorem_ids())
    missing = sorted(set(ids) - proved)
    if missing:
        # A backing claim naming a theorem the manifest does not carry is a bug,
        # not a degraded warrant. Refuse to emit rather than emit something
        # weaker than it looks.
        raise ValueError(
            f"{certificate_file}: claims formal backing by {missing}, "
            "which the theorem manifest does not list as proved"
        )
    return {
        "formal_manifest_id": mid,
        "formal_theorem_ids": list(ids),
        "numeric_warrant": numeric_warrant,
        "logical_implication_warrant": FORMAL_WARRANT,
    }


def describes_only_an_implication(cert: Dict[str, Any]) -> bool:
    """True for a formal certificate: it proves an implication, not a number.

    A consumer asking "is this block PSD?" must not be satisfied by one. The
    formal certificate says *if* the bounds hold *then* the block is definite;
    only the numeric certificate says the bounds hold.
    """
    return cert.get("content_kind") == KIND_FORMAL


# --------------------------------------------------------------------------- #
# the formal certificate itself                                                #
# --------------------------------------------------------------------------- #
def build_formal_certificate() -> Dict[str, Any]:
    """The FORMAL_THEOREM_CERTIFICATE body.

    Carries no numeric bound and no normalization binding, because it depends
    on neither: every theorem it reports is finite linear algebra over the
    reals, true independently of which pole primitive Atlas adopted. What it
    does carry is the manifest id, the theorem list, the axiom set, and the
    Lean source hashes -- so a change to any statement invalidates it the same
    way a source change invalidates a numeric certificate.
    """
    m = load_manifest()
    if m is None:
        raise FileNotFoundError(f"missing theorem manifest at {MANIFEST_PATH}")
    theorems = m.get("theorems") or []
    axioms = sorted({a for t in theorems for a in (t.get("axioms") or ())})
    body: Dict[str, Any] = {
        "certificate_version": "0.1",
        "program": "RH/Weil finite theorem boundary",
        "work_order": "ATLAS-RH-ENG-007",
        "content_kind": KIND_FORMAL,
        "claim_scope": "finite_dimensional_weil_compression",
        "rh_proof_claim": False,
        "status": "PASS" if theorems else "FAIL",
        # Not a rung of the numeric ladder. `promotion.is_rigorous` reads this
        # and correctly answers False: there is no interval computation here to
        # bind to a normalization.
        "evidence_class": FORMAL_WARRANT,
        "rigorous": False,
        "hard_constraints_certified": False,
        "numeric_warrant": None,
        "logical_implication_warrant": FORMAL_WARRANT,
        "psd_claim": False,
        "formal_manifest_id": m.get("manifest_id"),
        "formal_project": m.get("formal_project"),
        "lean_toolchain": m.get("lean_toolchain"),
        "mathlib_commit": m.get("mathlib_commit"),
        "formal_theorem_ids": [t["id"] for t in theorems],
        "theorems": [
            {
                "id": t["id"],
                "trusted_statement": t.get("trusted_statement"),
                "solution_theorem": t.get("solution_theorem"),
                "statement_hash": t.get("statement_hash"),
            }
            for t in theorems
        ],
        "axioms": axioms,
        "allowed_axioms": list(m.get("allowed_axioms") or ()),
        "unproved_statements": [dict(u) for u in (m.get("unproved_statements") or ())],
        "backs_certificates": {
            name: list(ids) for name, ids in sorted(FORMAL_BACKING.items())
        },
        "note": (
            "A formal theorem strengthens an exact theorem dependency. It does "
            "not convert interval numerical evidence to FORMAL: the numeric "
            "warrant of every certificate listed under backs_certificates is "
            "unchanged by this artifact."
        ),
        "dependencies": {"source_hashes": hashed_sources()},
    }
    return body


def formal_certificate_problems(cert: Dict[str, Any]) -> List[str]:
    """Everything wrong with a formal certificate, for the runner and the tests."""
    problems: List[str] = []
    if cert.get("content_kind") != KIND_FORMAL:
        problems.append("not a FORMAL_THEOREM_CERTIFICATE")
    if cert.get("rh_proof_claim") is not False:
        problems.append("missing rh_proof_claim: false")
    if cert.get("psd_claim") is not False:
        problems.append("a formal certificate must not claim PSD of anything")
    if cert.get("numeric_warrant") is not None:
        problems.append("a formal certificate must not carry a numeric warrant")
    if cert.get("logical_implication_warrant") != FORMAL_WARRANT:
        problems.append("logical_implication_warrant is not FORMAL")
    if cert.get("rigorous") is not False:
        problems.append("a formal certificate is not a rigorous numeric artifact")
    live = manifest_id()
    if live is None:
        problems.append("no theorem manifest available")
    elif cert.get("formal_manifest_id") != live:
        problems.append(
            f"stale formal_manifest_id {cert.get('formal_manifest_id')} (live {live})"
        )
    allowed = set(cert.get("allowed_axioms") or ())
    extra = sorted(set(cert.get("axioms") or ()) - allowed)
    if extra:
        problems.append(f"depends on non-standard axioms {extra}")
    for u in cert.get("unproved_statements") or ():
        if u.get("warrant") is not None:
            problems.append(f"unproved statement {u.get('id')} carries a warrant")
    return problems
