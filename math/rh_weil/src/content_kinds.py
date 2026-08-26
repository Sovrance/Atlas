"""The registry of PIR content kinds and what each one licenses (ENG-008 §WO-RH-54).

Content kinds are load-bearing in this program: "this block has inertia
``(2,0,0)``" and "this block is positive semidefinite" are different statements,
and a consumer that requires PSD must never be satisfied by the first. ENG-006
established that rule; what it did not establish was a place to record it, so the
answer lived in three places -- a tuple in ``inertia.certificate``, a second tuple
in ``pir_bridge``, and a frozen literal in ``scripts/ci_inertia.py``.

That drifted immediately. ENG-007 added ``FORMAL_THEOREM_CERTIFICATE`` and
``WEIL_PILOT_CONDITIONING_PREVIEW`` to ``pir_bridge`` and did not update the
gate's literal, so ``rh-inertia-fast`` went red on merge; and
``WEIL_DEGREE3_POSITIVITY_CERTIFICATE`` was being emitted and promoted while
appearing in no declared list at all.

This module is the one place a kind is declared. Adding a kind means adding an
entry here with an explicit answer to "may this ever satisfy a PSD consumer?",
and the gate fails if a kind reaches PIR without one. The decision itself is
still made by :func:`inertia.certificate.satisfies_psd_requirement`; what is
recorded here is the *declared* answer, so the gate can check the predicate
against the declaration rather than against a hand-maintained list.

Dependency direction is deliberate: this module imports from ``inertia``, never
the other way. The ``inertia`` package depends on nothing but the standard
library, which is what lets its exact tests run in any environment.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

from inertia.certificate import KIND_INERTIA, KIND_STRATIFICATION

#: Positivity certificates. These *claim* PSD and carry the signature to back it.
KIND_DEGREE3_POSITIVITY = "WEIL_DEGREE3_POSITIVITY_CERTIFICATE"
KIND_DEGREE4_POSITIVITY = "WEIL_DEGREE4_POSITIVITY_CERTIFICATE"

#: Derived-quantity certificates. Each says something true and none says PSD.
KIND_RANK_TRACE = "WEIL_RANK_TRACE_CERTIFICATE"
KIND_SPECTRAL_MOMENT = "WEIL_SPECTRAL_MOMENT_CERTIFICATE"

#: The formal channel (ENG-007 §12): an implication, never a number.
KIND_FORMAL = "FORMAL_THEOREM_CERTIFICATE"

#: Heuristic previews (E3). Plan inputs; never a warrant.
KIND_PILOT_PREVIEW = "WEIL_PILOT_CONDITIONING_PREVIEW"
KIND_SCAN_PREVIEW = "WEIL_SCAN_PREVIEW"

#: ENG-009 (§WO-RH-63). The generalized-gap kind is deliberately *not* PSD
#: licensable even though its content implies positivity when the bound is
#: positive: a consumer that needs PSD should consume the positivity
#: certificate that says PSD, not infer it from a gap bound whose purpose is
#: cross-dimension comparison. Default-deny until a work order justifies
#: otherwise.
KIND_GENERALIZED_GAP = "WEIL_GENERALIZED_GAP_CERTIFICATE"

#: ENG-010 (§WO-RH-74). The 4x4 positivity kind is PSD-licensable for the same
#: reason the degree-3/4 kinds are: it *claims* PSD and carries the certified
#: minors to back it. The adjudication kind records a verdict about E3 models
#: and licenses nothing.
KIND_DEGREE6_POSITIVITY = "WEIL_DEGREE6_POSITIVITY_CERTIFICATE"
KIND_SCALING_ADJUDICATION = "WEIL_SCALING_ADJUDICATION"
KIND_STRUCTURAL_DIAGNOSTIC = "WEIL_STRUCTURAL_DIAGNOSTIC"
KIND_SCALING_MODEL = "WEIL_SCALING_MODEL"
KIND_NEXT_BLOCK_SELECTION = "WEIL_NEXT_BLOCK_SELECTION"


@dataclass(frozen=True)
class ContentKind:
    """What a content kind is, and what a consumer may conclude from it."""

    name: str
    work_order: str
    summary: str
    #: May a certificate of this kind, at its best, satisfy a consumer that
    #: explicitly requires positive semidefiniteness? This is a property of the
    #: *kind*, not of any particular signature: an inertia artifact answers
    #: ``False`` even when its signature is ``(n, 0, 0)``.
    psd_licensable: bool
    #: Which question this kind's warrant grades -- the value of a number
    #: ("numeric"), the validity of a step ("implication"), or nothing at all
    #: ("preview").
    warrant_role: str


_KINDS: Tuple[ContentKind, ...] = (
    ContentKind(
        KIND_INERTIA, "WO-RH-28/34",
        "a certified signature (n+, n-, n0) for one block on one cell",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_STRATIFICATION, "WO-RH-34",
        "a partition of the cell into strata of constant signature, with bounded "
        "inconclusive transition regions",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_DEGREE3_POSITIVITY, "WO-RH-34",
        "the odd degree-3 block is positive definite on the cell, with certified "
        "lower bounds on the leading entry and the determinant",
        psd_licensable=True,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_DEGREE4_POSITIVITY, "WO-RH-51",
        "the 3x3 even block {1, b, b^2} is positive definite on the cell, with "
        "certified lower bounds on all three leading principal minors",
        psd_licensable=True,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_RANK_TRACE, "WO-RH-30",
        "a rank lower bound from the trace / Hilbert-Schmidt inequality, under "
        "hypotheses checked at the call site",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_SPECTRAL_MOMENT, "WO-RH-31",
        "enclosures of m1..m4 as traces of matrix powers, with the B1 "
        "truncated-moment queries they support",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_FORMAL, "WO-RH-44",
        "machine-checked finite theorems: which implications Atlas may use, "
        "under which pinned toolchain and axioms",
        psd_licensable=False,
        warrant_role="implication",
    ),
    ContentKind(
        KIND_PILOT_PREVIEW, "WO-RH-46",
        "a floating conditioning and topology preview of a block being prepared "
        "for a later work order",
        psd_licensable=False,
        warrant_role="preview",
    ),
    ContentKind(
        KIND_SCAN_PREVIEW, "WO-RH-48",
        "an independent non-interval assembly of a block, for regression against "
        "the rigorous one",
        psd_licensable=False,
        warrant_role="preview",
    ),
    ContentKind(
        KIND_DEGREE6_POSITIVITY, "WO-RH-69",
        "the 4x4 even block {1, b, b^2, b^3} is positive definite on the cell, "
        "with certified lower bounds on all four leading principal minors",
        psd_licensable=True,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_SCALING_ADJUDICATION, "WO-RH-71",
        "the verdict of a certified result against preregistered E3 scaling "
        "models, recorded before any refit; adjudicates plans, asserts no new "
        "numeric fact",
        psd_licensable=False,
        warrant_role="preview",
    ),
    ContentKind(
        KIND_GENERALIZED_GAP, "WO-RH-58",
        "a certified enclosure of lambda_min(G, M) for one block against the "
        "named exact reference metric, by shifted positivity and a Rayleigh "
        "witness -- basis-invariant by simultaneous congruence",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_STRUCTURAL_DIAGNOSTIC, "WO-RH-56",
        "cross-block structural data regenerated from promoted certificates: "
        "minors, traces, moments, conditioning, inertia, warrants",
        psd_licensable=False,
        warrant_role="numeric",
    ),
    ContentKind(
        KIND_SCALING_MODEL, "WO-RH-60",
        "exploratory finite scaling models fitted to certified per-block data, "
        "each with an explicit next-block falsifier; E3, never promotable",
        psd_licensable=False,
        warrant_role="preview",
    ),
    ContentKind(
        KIND_NEXT_BLOCK_SELECTION, "WO-RH-62",
        "the scored selection of the next block to certify, with the criteria "
        "and the losing candidates recorded",
        psd_licensable=False,
        warrant_role="preview",
    ),
)

REGISTRY: Dict[str, ContentKind] = {k.name: k for k in _KINDS}

#: Every declared kind, in declaration order. ``pir_bridge`` publishes this.
CONTENT_KINDS: Tuple[str, ...] = tuple(k.name for k in _KINDS)

#: Kinds that may never satisfy a PSD-requiring consumer.
NON_PSD_KINDS: Tuple[str, ...] = tuple(
    k.name for k in _KINDS if not k.psd_licensable
)

VALID_WARRANT_ROLES = ("numeric", "implication", "preview")


def is_registered(kind: str) -> bool:
    return kind in REGISTRY


def describe(kind: str) -> ContentKind:
    try:
        return REGISTRY[kind]
    except KeyError:
        raise KeyError(
            f"unregistered content kind {kind!r}; declare it in "
            "src/content_kinds.py with an explicit psd_licensable answer"
        ) from None


def psd_licensable(kind: str) -> bool:
    """May this kind ever satisfy a PSD requirement? Unregistered kinds cannot."""
    entry = REGISTRY.get(kind)
    return bool(entry and entry.psd_licensable)


def unregistered(kinds) -> Tuple[str, ...]:
    """Which of ``kinds`` have no declaration here."""
    return tuple(sorted(k for k in kinds if k not in REGISTRY))
