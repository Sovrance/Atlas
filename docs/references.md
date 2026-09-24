# Research Paper Ledger

Papers researched for the program, with role classification. Roles:
**ANCHOR** (external evidence cited in an edit/certificate) · **METHOD** (shapes the pipeline/discipline) · **TARGET** (source of data or falsifiers) · **CANDIDATE** (queued, not yet assessed).

| Ref | Paper | Role | Where it acts |
|---|---|---|---|
| R01 | Isi, Farr, Giesler, Scheel, Teukolsky, PRL 127, 011103 (2021) — GW150914 area-law test | ANCHOR | ?₃ / edit-001 / B4 |
| R02 | LVK, "GW250114: Testing Hawking's area law and the Kerr nature of black holes," PRL 135, 111403 (2025) | ANCHOR | ?₃ / edit-001 / B4 REAL |
| R03 | arXiv:2509.03480 — GW230814/GW231226 area law ≳5σ (GWTC-4) | ANCHOR | ?₃ / edit-001 |
| R04 | Curto & Fialkow — flat extensions of positive moment matrices | ANCHOR | B1 forcing mechanism |
| R05 | Choi (1975) — completely positive maps; Stinespring dilation | ANCHOR | B2 gate |
| R06 | Weinberg, *QTF* Vol. 1 — cluster decomposition | ANCHOR | ?₅ / edit-002 / B5 |
| R07 | Haag, *Local Quantum Physics* | ANCHOR | ?₅ / edit-002 |
| R08 | Wall (2011) — GSL proof; QNEC for free fields | ANCHOR | ?₂ / edit-003 / B6 |
| R09 | Balakrishnan, Faulkner, Khandker, Wang (2017) — QNEC from modular theory | ANCHOR | ?₂ / edit-003 / B6 |
| R10 | Jacobson (1995) — Einstein equations as equation of state | ANCHOR | ?₁ |
| R11 | Adams et al. — positivity bounds from causality/analyticity/unitarity | METHOD | feasible set 𝔉; future EFT gate |
| R12 | Schmidt & Lipson; AI Feynman; SINDy | METHOD | decoding-chain step 5 (candidate generators only) |
| R13 | Chamseddine & Connes — spectral action | METHOD (hypothesis) | bridge formula; flagged conjectural |
| R14 | **Maudlin, "Actual Physics, Observation, and Quantum Theory," arXiv:2512.22618 (physics.hist-ph, Dec 2025)** | **METHOD** | **𝓜-layer (decoding-chain step 1); certificate format; arrival-time track** |
| R15 | **Chen, Soh, Ooi, Vissol-Gaudin, Yu, Novoselov, Hippalgaonkar, Li — "Constructing custom thermodynamics using deep learning," Nat. Comput. Sci. 4, 66–85 (2024); arXiv:2308.04119; code: Zenodo 10.5281/zenodo.10212239** | **METHOD** | **decoding-chain step 5 (constraint-native generators); closure-as-d_identifiable; B7 + B8 done; Atlas Engine seeded** |

## R14 annotation (added 2026-07-11)

Philosophy-of-physics essay; produces no theorems or certificates and changes no matrix cell. Adopted contributions (see `docs/notes/measurement-interface-maudlin.md`):
1. **M-layer stipulation ledger** — every certificate must list its un-derived measurement-interface stipulations explicitly (Einstein's clocks-and-rods "incurred obligation," systematized).
2. **Scope clause for process-type certificates (B2 class)** — completions certify consistency *given* a device↔POVM association; the association itself is an M-layer stipulation, not a certified object.
3. **Arrival-time track (candidate held-out domain)** — standard QM supplies no unique arrival-time observable (Pauli's theorem: spectral positivity of H forbids a conjugate self-adjoint T — an atlas-native restraint fact); rival completions differ measurably. Rare domain where candidate rule-objects 𝕽 genuinely disagree → high U_prediction value in the Score function.

## R15 annotation (added 2026-07-11)

Open-access research article (pointer received via its paywalled Research Briefing, doi:10.1038/s43588-023-00590-4). S-OnsagerNet: learns closed macroscopic stochastic dynamics Ż = −M∇V(Z) (+ stochastic extension) from microscopic trajectories, with M ⪰ 0 (dissipation) and generalized-potential structure ENFORCED BY ARCHITECTURE, while simultaneously learning the minimal closure coordinates (polymer stretching: ~900 DOF → 3 coordinates; DNA-experiment validation). Adopted contributions (see `docs/notes/onsager-sonsagernet-chen2024.md`):
1. **Constraint-native candidate generation** — restraint bundles embedded in the generator's hypothesis space (structured) vs applied as post-hoc gates (unstructured, SINDy/R12). Firewall: by-construction constraints characterize the hypothesis class, never the learned object; certification classes and held-out gates unchanged.
2. **Closure coordinates = dynamical d_identifiable** — algorithmic pattern for extending rank discovery (B3/B5) from static relations to evolution laws.
3. **Benchmark B7 implemented: Onsager transport-matrix completion** — `b7_onsager/` + certificate; PSD alone → PERMITTED; +reciprocity → FORCED.

**R15 addendum (2026-07-12):** cross-agent review filed (`docs/notes/r15-addendum-cross-agent-review.md`). Adopted: canonicalization rule, distill-before-score rule for L(𝕽), **B8 implemented** (blind grammar identification), certified closure dimension research item. Rejected: single-number "encoding N" endpoint. Erratum: full S-OnsagerNet law is Ż = −[M(Z)+W(Z)]∇V(Z)+σ(Z)Ḃ with W antisymmetric.

---

## R16–R22: The 2025 Nobel block (added 2026-07-13) — verbose provenance

Added as a group because the 1984–85 Berkeley sequence is, to date, the cleanest HISTORICAL EXECUTION of this program's discipline formula: independent restraint-wise calibration → structure frozen → parameter-free prediction on a held-out DOMAIN (quantum regime) → promotion. Benchmark B9 encodes it as executable gates. Verified against Royal Swedish Academy materials and secondary coverage (Physics World, Physics Today, 2025-10-07/10).

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R16 | Royal Swedish Academy of Sciences, *Scientific Background: Nobel Prize in Physics 2025* (nobelprize.org, Oct 2025); prize "for the discovery of macroscopic quantum mechanical tunnelling and energy quantisation in an electric circuit" (Clarke, Devoret, Martinis) | METHOD/TARGET | Names the decisive issue in the pre-1984 literature: low-temperature switching saturation was AMBIGUOUS between quantum tunneling and unmodeled environmental noise; resolution required eliminating excess noise and measuring junction/environment parameters INDEPENDENTLY rather than fitting them. That ambiguity-and-resolution pattern is B9's T2/T3/T5 gates and the origin of the "calibration route" field in certificates. |
| R17 | Devoret, Martinis, Esteve, Clarke, "Resonant Activation from the Zero-Voltage State of a Current-Biased Josephson Junction," PRL 53, 1260 (1984) | METHOD | Resonant activation as IN-SITU SYSTEM IDENTIFICATION: plasma frequency ω_p(I), quality factor Q, hence C and R measured on the working device. Traced into B9 as the "SPECTROSCOPY route" for H_V — the independent calibration channel that T3 proves is the degeneracy-breaking ingredient. Also the historical example of decompiling a hidden electrical network from its resonance response. |
| R18 | Martinis, Devoret, Clarke, "Energy-Level Quantization in the Zero-Voltage State of a Current-Biased Josephson Junction," PRL 55, 1543 (1985) | ANCHOR (ħ row) + METHOD | Discrete levels of a MACROSCOPIC degree of freedom, with resonance positions matching a Schrödinger calculation whose parameters came from the CLASSICAL regime, no refitting. (a) Empirical macroscopic support for the ħ row (Pos/Top entries). (b) The parameter-free held-out protocol = promotion criterion 5, executed 1985; miniaturized in B9-T6 (exact rational charge-basis H → asymptotic prediction untouched by calibration). |
| R19 | Devoret, Martinis, Clarke, "Measurements of Macroscopic Quantum Tunneling out of the Zero-Voltage State of a Current-Biased Josephson Junction," PRL 55, 1908 (1985) | ANCHOR/METHOD | Temperature-independent low-T escape agreeing with dissipative-MQT theory with NO adjustable junction parameters. Two imports: the Γ_q formula's exponent carries Q — dissipation (the M matrix) shifts a rare-event rate by orders of magnitude, so PERMITTED intervals on M propagate EXPONENTIALLY into rate predictions (interval-arithmetic lesson for statistical certificates); and switching histograms become hazard functions via Γ(I) = İ·p(I)/S(I), not free-form fit targets. |
| R20 | Martinis, Devoret, Clarke, PRB 35, 4682 (1987) "Experimental tests for the quantum behavior of a macroscopic degree of freedom"; Cleland, Martinis, Clarke, PRB 37, 5950 (1988); Clarke, Cleland, Devoret, Esteve, Martinis, Science 239, 992 (1988) | ANCHOR | Full synthesis + controlled-dissipation test: a shunted junction suppressed tunneling ~×300 as predicted — the M matrix VALIDATED as a causal dial, not a nuisance parameter. Supports treating dissipative structure as first-class physics in the GENERIC/B8/B9 grammar. |
| R21 | Caldeira & Leggett, PRL 46, 211 (1981); Ann. Phys. 149, 374 (1983); Leggett's late-1970s proposals | METHOD (theory layer) | The theory that made dissipation QUANTITATIVE in quantum tunneling (environment as oscillator bath; friction suppresses tunneling). Establishes the M-column ↔ quantum-rate coupling B9 gates against, and the theoretical basis for "M(ω) or enlarged Markovian state space" when white-noise FDT fails. |
| R22 | Precursor/classical layer: Josephson (1962); Stewart (1968) & McCumber (1968) RCSJ model; Ivanchenko–Zil'berman (1969); Kramers (1940) escape theory; Fulton & Dunkleberger, PRB 9, 4760 (1974); Büttiker–Harris–Landauer, PRB 28, 1268 (1983); precursor observations Voss & Webb, PRL 47, 265 (1981), Jackel et al. (1981), Washburn et al. (1985) | TARGET/METHOD | Supplies (a) the RCSJ → tilted-washboard mapping B9's model.py linearizes, (b) the CLASSICAL NULL HYPOTHESIS (thermal escape) every quantum claim had to defeat, and (c) the cautionary precursors: correct-looking saturation with unresolved noise ambiguity. B9-T2/T3 exist because of this row — the program treats "plausible signal, uncertified environment" as REJECTED-pending, permanently. |

**Cross-agent note:** a second AI analysis of this block was reviewed the same day; convergences, corrections, and adopted/rejected items in `docs/notes/nobel2025-circuit-decompilation.md`.

---

## R23–R24: Experimental-program block (added 2026-07-14) — verbose provenance

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R23 | Howl & Fuentes, "Quantum Frequency Interferometry: with applications ranging from gravitational wave detection to dark matter searches," arXiv:2103.02618 | METHOD/TARGET | The estimation-theory half of the coupling-matrix machine: frequency modes trapped in one volume, in contact at all times, permitting estimation of GLOBAL multi-mode Gaussian channels (mode-mixing + two-mode squeezing generators) with QFIM bounds — i.e., experimental reconstruction of exactly the object stack our B2 (Choi), B7 (susceptibility), B8 (grammar), B9 (M/W/σ) machinery certifies. Direct seed of EXP-A/B10; the QFIM ⪰ 0 layer is the Pos column at the metrology tier. Read directly 2026-07-14. |
| R24 | Zhang, del Aguila, Mazzoni, Poli, Tino, "A trapped atom interferometer with ultracold Sr atoms," arXiv:1609.06092 | METHOD (engineering) | The platform half: Ramsey–Bordé Bragg + Bloch oscillations holding ⁸⁸Sr against gravity, 1 s in-lattice interference, decoherence budget; Sr's zero nuclear spin and >100 s Bloch coherence. Contributions traced: trapped long interrogation without towers; loss-vs-decoherence discrimination; T-scaling diagnostics; the platform class feeding future h/m → α inputs (EXP-D). Read directly 2026-07-14. |

Companion analysis + cross-agent review: `docs/notes/experimental-program-bec.md`. Adopted from the reviewed second-agent document: BENCHMARK/ANALOGUE-REAL/FUNDAMENTAL-REAL promotion policy (verbatim), falsifiability pentad, Phase A/B/C demonstrator, experiment-record schema (+2 mandatory fields: calibration_route, m_layer_stipulations). Rejected: G/Λ rows as primary cells of a BEC classifier (tier violation); finite-rank design lacking the information-Gram refinement (S4).

---

## R25–R30: Related-programs block (added 2026-07-14) — verbose provenance

Filed with `docs/notes/related-programs-review.md` (independent sweep + cross-agent comparison) and `docs/SPECIFICATION.md`.

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R25 | Wan & Zhou, "Matrix moment approach to positivity bounds and UV reconstruction from IR," JHEP 02 (2025) 168, arXiv:2411.11964 | METHOD (convergence witness) | The physics community independently running B1's mathematics (moment matrices, positivity cones, reconstruction from truncated data) inside EFT positivity bounds. Elevates our moment/flat-extension layer from "one benchmark" to the central verifier, with an external anchor. |
| R26 | Arkani-Hamed, Huang & Huang, "The EFT-Hedron," JHEP 05 (2021) 259, arXiv:2012.15849; He & Kruczenski, "S-matrix bootstrap in 3+1D: regularization and dual convex problem," JHEP 08 (2021) 125, arXiv:2103.11484 | METHOD | Convex geometry of allowed theory space; DUAL certificates that exclude. Import: FORCED boundaries as extremal rays (B10's quantum-limited flatness is one); add dual/exclusion functionals to the verifier layer (roadmap M-item). |
| R27 | Guerrieri, Häring & Su, "From data to the analytic S-matrix: a bootstrap fit of the pion scattering amplitude," arXiv:2410.23333 | METHOD/TARGET | Bootstrap meets REAL data — the community's B4/S3 analogue. Confirms REAL-mode certificates as a shared frontier; candidate future cross-check domain. |
| R28 | Hardy, "Reformulating and Reconstructing Quantum Theory," arXiv:1104.2066 (GPT reconstruction lineage) | METHOD | Operational axioms recover quantum structure — our program restricted to the ħ row, done rigorously. Import: neutral operational substrate (SPEC §1) and GPT rival family for B12-RGRC (guards Pos against Hilbert-space circularity). |
| R29 | Categorical quantum mechanics lineage (Abramsky–Coecke); symmetric monoidal process theories | METHOD | The Cmp column's proper mathematics: typed processes, ⊗, †, diagram equivalence. Import: the 𝒞-skeleton of SPECIFICATION §1; process equivalence as the syntactic layer of "same." |
| R30 | Equation/law discovery block: SINDy lineage (R12); differential-invariants discovery arXiv:2505.18798; AI Poincaré 2.0 arXiv:2203.12610; causal discovery review arXiv:2305.13341 | METHOD | Mature generator-side machinery (candidate enumeration, sparsity, train/test, intervention-based identifiability). Import: M2 generator engine feeding OUR certified verifiers; NONIDENTIFIABLE(cause) vocabulary; intervention metadata in B12-RGRC. |

---

## R31–R33: α-row tension block (added 2026-07-15) — verbose provenance

Input to S3-EM (edit-009 memo). All values web-verified 2026-07-15.

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R31 | Morel, Yao, Cladé, Guellati-Khélifa (LKB), "Determination of the fine-structure constant with an accuracy of 81 parts per trillion," *Nature* **588**, 61 (2020) | ANCHOR (α row) | Rb photon-recoil α⁻¹ = 137.035999206(11); identified new beam-profile-class systematics as the likely explanation for the earlier 2.4σ shift in the Rb series. S3-EM T1 pull, T3 jackknife pivot. |
| R32 | Parker, Yu, Zhong, Estey, Müller (Berkeley), "Measurement of the fine-structure constant as a test of the Standard Model," *Science* **360**, 191 (2018) | ANCHOR (α row) | Cs recoil α⁻¹ = 137.035999046(27); the Rb–Cs gap this creates is the >5σ tension analyzed by S3-EM. Provides the second h/m channel that makes localization possible. |
| R33 | Hanneke, Fogwell, Gabrielse, *PRL* **100**, 120801 (2008) [electron g−2] + Aoyama, Kinoshita, Nio, *PRD* **97**, 036001 (2018) [QED-input α⁻¹]. Newer a_e determinations flagged VERIFY-BEFORE-USE, not input | ANCHOR/METHOD | Third, independent channel via g−2 route — the "tie-breaker" the jackknife exploits. |

---

## R34–R35: Gravity-from-Entropy block (added 2026-07-16) — verbose provenance

Review + code: `docs/notes/gfe-bianconi-review.md`, `probes/truncation_audit.py`, M1 `metric_pair` kind. Both papers read directly 2026-07-16.

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R34 | Bianconi, "Gravity from entropy," arXiv:2408.14391 (v7, 2025), Phys. Rev. D lineage | TARGET (rival grammar) + METHOD (imported mathematics) | Entropic action = GQRE between spacetime metric and matter-induced metric; Dirac–Kähler 0⊕1⊕2-form matter; auxiliary G-field → dressed Einstein–Hilbert with emergent positive Λ(G). Imports: metric-pair invariant spectrum (→ M1 `metric_pair`), Burg/Stein divergence core (→ canon). Rival typing of Λ filed against ?₇/?₈; rival grammar for ?₄. Truncation status: UNAUDITED per our gate (probes/truncation_audit.py) pending a consistency proof — the audit demand is now executable, not rhetorical. |
| R35 | Bianconi, "The Thermodynamics of the Gravity from Entropy Theory," arXiv:2510.22545 (v4, 2026) | TARGET (rival grammar, cosmological layer) | GfE FRW thermodynamics: k-temperatures/k-pressures with a first law; dynamical Λ_G dark energy; entropy growth reconciled with GQRE/volume non-increase; Schwarzschild area law from GfE first principles (companion). Marks the future FUNDAMENTAL-REAL comparison surface (dark-energy dynamics) and an R-KMS-adjacent Thm structure for a future FRW toy sprint. |

---

## R36–R40: Surfaceology / hidden-zeros block (added 2026-09-23) — verbose provenance

Input to B15-SURF (`docs/preregistrations/prereg-002-surfaceology-benchmarks.md`, `b15_surf/`, edit-011). All arXiv IDs re-verified 2026-09-23 at transcription; every convention used in code is transcribed with equation numbers in prereg-002. R26 (EFT-hedron) already covers B15-POS and is not duplicated.

| Ref | Source | Role | Original contribution traced into this project |
|---|---|---|---|
| R36 | Arkani-Hamed, Frost, Salvatori, Plamondon, Thomas, *All Loop Scattering as a Counting Problem*, arXiv:2309.15913 (JHEP 08 (2025) 194); *All Loop Scattering for All Multiplicity*, arXiv:2311.09284 (JHEP 09 (2025) 033) | METHOD | Curve-integral formalism: amplitudes as counting problems on surfaces; `u`-variables satisfying `u_C + Π u_D^{n(C,D)} = 1` make locality/factorization structural rather than checked. Imports: the representation-redundancy pattern (P1/P7 — one geometric object vs. ~4ⁿ diagrams → `redundancy_metric` recorded by B15-GID) and the constraint-manifest design question of WP4 (`docs/notes/constraint-manifest-representation-v0.md`). |
| R37 | Arkani-Hamed, Cao, Dong, Figueiredo, He, *Hidden zeros for particle/string amplitudes and the unity of colored scalars, pions and gluons*, arXiv:2312.16282 (JHEP 10 (2024) 231); *NLSM ⊂ Tr(Φ³)*, arXiv:2401.05483 (PRD 110 (2024) 065018); *Scalar-Scaffolded Gluons and the Combinatorial Origins of Yang-Mills Theory*, arXiv:2401.00041; Arkani-Hamed & Figueiredo, *All-order splits and multi-soft limits for particle and string amplitudes*, arXiv:2405.09608 | TARGET (B15-ZERO ground truth) + METHOD (zeros/splits as fingerprints; δ-shift as base grammar + deformation) | Transcribed and executed exactly: planar/mesh conventions (2312.16282 eqs. 2.2–2.4), hidden-zero loci (eq. 6.14 + §3.2, n(n−3)/2 loci), near-zero factorization (eq. 3.3/3.8), δ-shift and NLSM limit (2401.05483 eqs. 1–2; 2312.16282 eq. 7.20), 2-splits (2405.09608 eqs. 1.7, 2.1, 2.7, 2.10). Every statement holds exactly in `Fraction` at n ≤ 8; one index-range reading (2405.09608 eq. 2.13) filed as erratum-candidate E1 in edit-011. YM-scaffolded rows (2401.00041) deferred (§9-OI-3). |
| R38 | Arkani-Hamed, Bai, He, Yan, *Scattering Forms and the Positive Geometry of Kinematics, Color and the Worldsheet*, arXiv:1711.09102 | METHOD | ABHY associahedron whose canonical form is the Tr(φ³) tree amplitude; the geometric origin of the causal-diamond zeros. Reserved for Route C (WP2.3, not attempted this sprint; HEURISTIC until orientation conventions are pinned). |
| R39 | Cheung, Dersy, Schwartz, *Learning the Simplicity of Scattering Amplitudes*, arXiv:2408.04720; Moynihan, *Learning the S-matrix from data: Rediscovering gravity from gauge theory via symbolic regression*, arXiv:2602.15169 (preprint, Feb 2026) | METHOD (backlog) | ML / symbolic rediscovery precedents — transformer simplification of spinor-helicity expressions; PySR rediscovery of KK/BCJ/Parke–Taylor/KLT to 5 points (fails at 6). Import: ML as *conjecture generator only*; every proposal must pass `tools/verify_b15_certificate.py`-class exact checks or be filed REJECTED/NONIDENTIFIABLE (work-order §10). |
| R40 | Guevara, Lupsasca, Skinner, Strominger, Weil, *Single-minus gluon tree amplitudes are nonzero*, arXiv:2602.12176 (preprint, Feb 2026) | METHOD (backlog) | AI-conjectured, human-proved amplitude result on a non-generic ("half-collinear") kinematic slice, checked by Berends–Giele recursion. Narrow slice; preprint; VERIFY-BEFORE-USE. Precedent for the numeric→symbolic→exact-certificate pipeline, not for any claim here. |
