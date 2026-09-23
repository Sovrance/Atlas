# Surfaceology, Carolina Figueiredo's "Hidden Zeros," and What Their Discovery Path Teaches About Simplifying Physical Law

**Provenance:** deep-research report produced 2026-09-21 (web sources; arXiv IDs verified at that date). Filed as required reading #10 for work order B15-SURF. Ledger entries R36–R40 derive from §"Key Findings" and §"AI/ML" below. Treat as reference material; every formula, convention, and ID is re-transcribed from the cited primary source before use in `prereg-002`.

## TL;DR
- **Surfaceology reformulates known perturbative QFT amplitudes; it does not discover new laws of nature.** It computes scattering amplitudes for Tr(φ³) colored scalars — and, via a single kinematic shift, pions (NLSM) and gluons (Yang–Mills) — as one "curve integral" over a surface, replacing the exponentially growing sum over Feynman diagrams at all loop orders and all orders in the topological (1/N) expansion. Carolina Figueiredo's specific contribution (with Nima Arkani-Hamed, Qu Cao, Jin Dong, Song He) was the discovery that these three physically distinct theories share the same "hidden zeros" and are related by a simple deformation — work recognized with the inaugural 2026 Vera Rubin New Frontiers Prize.
- **The transferable discovery heuristics are real and directly relevant to a constraint-certification program:** compute lots of explicit data and distrust unexplained simplicity; change to variables where hard constraints (positivity, unitarity, locality) are *manifest* rather than checked; use zeros/poles/factorization as basis-independent structural fingerprints; find a minimal base theory from which realistic ones descend by a simple deformation; and treat representation redundancy (many Feynman diagrams / gauge choices ↔ one geometric object) as the signature of a wrong representation. These map onto the 𝕽 program's canonicalization, rival-generator testing, and REPRESENTATION_DEPENDENT / OBSERVATIONALLY_EQUIVALENT verdicts with strong-to-moderate analogy strength.
- **The honest limit — and the sharpest caveat for the program:** surfaceology succeeds inside perturbation theory of *already-known* Lagrangian theories, where Feynman-diagram ground truth exists to check against. It has produced new *relations* and computational compression, not new predictions or constants of nature, and it does not yet cover real QCD with massive fermions, non-planar gravity, or the non-perturbative regime. A program lacking ground truth cannot borrow surfaceology's validation guarantees, only its heuristics — and the mature, rigorously-certified analogue to anchor on is the EFT-hedron / S-matrix / conformal bootstrap, which already carves certified PERMITTED/EXCLUDED regions from Sym+Pos+Uni+Cau using exact/semidefinite methods with primal (witness) and dual (certificate) structure.

## Key Findings

**1. What surfaceology is.** The "curve integral formalism," introduced by Nima Arkani-Hamed with Hadleigh Frost, Giulio Salvatori, Pierre-Guy Plamondon, and Hugh Thomas: *All Loop Scattering As A Counting Problem* (arXiv:2309.15913; JHEP 08 (2025) 194) and *All Loop Scattering For All Multiplicity* (arXiv:2311.09284; JHEP 09 (2025) 033). For Tr(φ³), the loop-integrated amplitude at all loop orders and all orders in the 1/N expansion is a single integral over "curve space," determined by a counting problem attached to curves on a surface (fatgraph), with no trace of the conventional sum over Feynman diagrams.

**2. The machinery.** Thicken a Feynman graph into a surface; draw curves recording left/right turns. Each curve gets a **g-vector** in an E-dimensional curve space whose cones form the "Feynman fan" (each top-dimensional cone = one Feynman diagram, propagators = the g-vectors spanning it), and a piecewise-linear **headlight function** α_C(t) computed from a 2×2 curve matrix as α_C = Trop(M₁₁) + Trop(M₂₂) − Trop(M₁₂) − Trop(M₂₁). The headlight function is the tropicalization of a **u-variable**, u_C = F₋₊F₊₋ / (F₋₋F₊₊); the u's satisfy binary-geometry non-crossing relations (u + ∏u = 1 type) that make locality/factorization structural. The formalism is the tropical limit of stringy (Koba–Nielsen) integrals and connects to the **ABHY associahedron** (arXiv:1711.09102), whose canonical form gives the Tr(φ³) tree amplitude, and to Teichmüller theory / the mapping class group (a "Mirzakhani kernel" removes overcounting). Lineage: the amplituhedron (Arkani-Hamed–Trnka, arXiv:1312.2007, N=4 SYM), but with **no supersymmetry required**.

**3. Figueiredo's discovery.** Carolina Figueiredo (B.Sc. IST Lisbon; Ph.D. Princeton 2026, advisor Arkani-Hamed; Harvard Society of Fellows) found that Tr(φ³), NLSM, and Yang–Mills amplitudes share the *same* hidden zeros — kinematic loci where the amplitude vanishes, invisible in the Feynman-diagram representation but manifest in the associahedron/kinematic-mesh picture. Paper: *Hidden zeros for particle/string amplitudes and the unity of colored scalars, pions and gluons* (Arkani-Hamed, Cao, Dong, Figueiredo, He; arXiv:2312.16282; JHEP 10 (2024) 231). Inaugural 2026 Vera Rubin New Frontiers Prize (Breakthrough Prize Foundation).

**4. How she found it.** Per Quanta (Charlie Wood, 25 Sep 2024): singularities (poles) had driven earlier discoveries; Figueiredo targeted the *numerators*, searching for configurations where the amplitude vanishes, by recasting Tr(φ³) amplitudes as associahedron volumes and looking for configurations that flatten the polytope. She then checked pions (brute-force Feynman diagrams) and found the same zeros; likewise gluons. The zeros are rigid: only one part of the curve-integral expression can be deformed while preserving them; one tuning yields pions, another gluons.

**5. The δ-shift.** Planar variables X_{i,j} = (p_i + … + p_{j−1})²: X_{e,e} → X_{e,e} + δ, X_{o,o} → X_{o,o} − δ (even/odd indices). Shifting Tr(φ³) reproduces NLSM to all loops (*NLSM ⊂ Tr(φ³)*, arXiv:2401.05483, PRD 110 (2024) 065018; see also *Circles and Triangles*, arXiv:2403.04826); with δ = 1/α′ in the stringy integral it produces scalar-scaffolded gluons (arXiv:2401.00041, JHEP 04 (2025) 078). The shift preserves all hidden zeros. **Exact prefactor and limit prescription: transcribe from arXiv:2401.05483, not from this note.**

**6. Zeros, splits, factorization.** Near (not on) a hidden zero, amplitudes exhibit a "2-split": factorization into two lower currents without a residue on a physical pole (Arkani-Hamed–Figueiredo, *All-order splits and multi-soft limits*, arXiv:2405.09608, JHEP 10 (2025) 077). The zero condition is set on a "maximal rectangle / causal diamond" of the kinematic mesh; relaxing a single variable c_⋆ ≠ 0 triggers the split. Generalizations: double-copy theories inherit the zeros (Bartsch et al., arXiv:2403.10594); hidden zeros ⇔ enhanced UV scaling and uniqueness (Rodina, arXiv:2406.04234, PRL 134 (2025) 031601); cosmological wavefunction (arXiv:2503.23579); massive theories (arXiv:2601.16860).

**7. Coverage (2024–2026).** Constructed: Tr(φ³) all loops/multiplicities; NLSM all loops; scalar-scaffolded gluons and the canonical planar non-SUSY YM all-loop integrand (arXiv:2408.11891, PRL 134 (2025) 171601); colored Yukawa toy fermions (arXiv:2406.04411). Open: massive fermions / real QCD; gravity (hints only); non-planar beyond conjecture; anything non-perturbative. 2026 extensions: non-orientable curve integral (JHEP 06 (2026) 249); loop-level zeros/2-splits (arXiv:2604.13810); forward-limit recursions (arXiv:2503.15860, arXiv:2509.25129).

**8. Scientific status.** Read by the community as a powerful reorganization ("exponential compactification in information" — Spradlin; Bourjaily checked the zero conspiracy to 14 points: "I have no idea why that's true"). Critics (Woit) note 17 years without the promised replacement for spacetime; Arkani-Hamed himself names non-perturbative blindness the approach's greatest weakness. Locality and unitarity are "emergent" only in the limited sense that they are outputs of the combinatorics/positivity, demonstrated for scalars and at one loop for Tr(φ³)/NLSM.

## Correcting the framing "she simplified Feynman diagrams"
- The curve-integral formalism predates her central result and was built by Arkani-Hamed, Frost, Salvatori, Plamondon, Thomas.
- Her signature contribution is the hidden zeros and the resulting unity of three theories — a discovery about *relations between theories*.
- The deeper claim is representational: the Feynman expansion is a redundant, gauge- and basis-dependent representation of a simpler object. That is the transferable lesson.

## Transferable discovery patterns
| ID | Pattern | Strength | Failure mode |
|---|---|---|---|
| P1 [ARCH] | Distrust unexplained simplicity; redundancy signals wrong representation | Strong | low-order coincidence without a generative mechanism |
| P2 [ARCH] | Variables in which hard constraints are manifest, not checked (spinor-helicity, momentum twistors, u-variables) | Strong — deepest import | theory-specific; can hide other constraints; resists generalization (planar → non-planar) |
| P3 [ARCH] | Zeros / poles / residues / factorization as basis-independent fingerprints | Strong (classification) / Moderate (discovery) | fingerprints can be observationally inaccessible — theory-internal, low measurement warrant |
| P4 [ARCH] | Minimal base theory + simple deformation family (Tr φ³ → NLSM → YM) | Moderate | δ-shift exists because all three share associahedral kinematics; no general guarantee; over-fitting |
| P5 [ARCH] | Locality/unitarity as outputs of geometry rather than inputs | Moderate | "consistent with" mistaken for "derived from" |
| P6 [ARCH] | Positivity / convex geometry as organizing constraint (EFT-hedron) | Strong | not all constraints convex; relaxations loose or over-tight |
| P7 [ARCH] | Counting / tropical problems replacing sums | Moderate | counting object itself intractable (Feynman-fan sampling open, arXiv:2503.07707) |
| P8 [ARCH] | Quotient by representation redundancy before comparing | Strong (concept) / Moderate (ops) | quotienting too much collapses distinct theories |
| P9 [ARCH] | Unexpected cancellations as shadows of hidden symmetry | Moderate | some cancellations are low-order accidents |
| P10 [TOOL] | Compute massive explicit exact data first, then mine | Strong | low-order data misleads; patterns without mechanism are conjecture |

## AI/ML precedents (all rediscovery or narrow)
- Cheung–Dersy–Schwartz, transformer simplification of spinor-helicity expressions (arXiv:2408.04720) — rediscovery/compression; outputs are conjectures needing numerical validation.
- Guevara, Lupsasca, Skinner, Strominger, Weil, *Single-minus gluon tree amplitudes are nonzero* (arXiv:2602.12176, Feb 2026) — key formula conjectured by an LLM, proved by humans/another model, checked by Berends–Giele; non-generic kinematic slice; preprint.
- Moynihan, *Learning the S-matrix from data* (arXiv:2602.15169, Feb 2026) — PySR rediscovers KK/BCJ/Parke–Taylor/KLT to 5 points, fails at 6.
- AI Feynman; Ramanujan Machine; FunSearch; AlphaEvolve — rediscovery or narrow record-breaks.
- Net: ML accelerates and rediscovers; keep all epistemic weight on the exact-certificate layer; prefer numeric→symbolic (checkable) over free-form transformers.

## Historical precedents
Maxwell → forms; Lagrange/Hamilton/symplectic; Noether; Wilsonian RG; Dirac; gauge bundles; twistors; categorical QM; quantum reconstructions (Hardy; CDP; Masanes–Müller); constructor theory; bootstrap. Discriminator between fertile and merely-elegant reformulations: did it make a new invariance manifest, reveal the right equivalence classes, or expose a new regime? If not, it is compression, not discovery.

## Recommendations (as adopted in work order B15-SURF)
- **Stage 0 → WP1 (B15-POS):** anchor on the EFT-hedron / bootstrap certificate stack (SDPB, primal/dual, exact arithmetic); reproduce a published bound with primal + dual certificates before trusting novel claims.
- **Stage 1 → WP2 (B15-ZERO):** 5-/6-point Tr(φ³) by two routes; hidden zeros; 2-split; δ-shift to NLSM — all exact.
- **Stage 2 → WP3 (B15-GID):** blind grammar identification from zero/factorization fingerprints; pre-registered for both separable and non-separable outcomes; the non-separable outcome is the more important epistemic lesson.
- **Stage 3 → WP4:** constraint-manifest representation design note.
- **Stage 4 → backlog:** ML as conjecture generator only.
- Plan-changing events: surfaceology certified for massive fermions/QCD; a non-perturbative or gravity surface formulation; a proven, generic AI-discovered amplitude result.

## Caveats
- Reformulation ≠ new law; surfaceology derives no constants of nature.
- Ground truth is the crucial disanalogy: the heuristics transfer, the epistemic warrant does not.
- Perturbative and largely planar; spacetime emergence remains aspirational.
- Hidden zeros are effectively unobservable; shared-fingerprint equivalence is theory-internal (E0 mathematical warrant, low measurement warrant) — respect the axis distinction.
- All mappings are analogies; P1, P2, P6, P8 are the genuine imports; P4, P5, P7, P9 are hypotheses under test.
- arXiv-ID note: 2602.12176 = single-minus gluon result; 2602.15169 = Moynihan symbolic regression. Both preprints.

## Sources
arXiv:2309.15913 · 2311.09284 · 2312.16282 · 2401.05483 · 2403.04826 · 2401.00041 · 2405.09608 · 2403.10594 · 2406.04234 · 2503.23579 · 2601.16860 · 2408.11891 · 2406.04411 · 2604.13810 · 2503.15860 · 2509.25129 · 1711.09102 · 1312.2007 · 2012.15849 · 1502.02033 · 1909.09745 · 2310.10729 · 2408.04720 · 2602.12176 · 2602.15169 · 1905.11481 · 2506.13131 · 2308.14789 · 2503.07707; Quanta Magazine, C. Wood, 25 Sep 2024; Breakthrough Prize Foundation, Vera Rubin New Frontiers Prize 2026.
