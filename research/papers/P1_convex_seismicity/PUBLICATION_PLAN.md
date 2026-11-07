> **C50 (2026-11-04): this plan is superseded by `../PORTFOLIO.md`.** Changes: the Äspö part moved to P3; the JGR long version is retired; Kim et al. 2025 (Otaniemi control) was added as the closest prior work; references were re-audited.

# P1 publication plan

## Readiness checklist
| Step | Status |
|---|---|
| Theory: convexity, dose invariance, peak-rate LP | done (proof sketches; write full proofs in an appendix) |
| Independent verification (ODE, convexity, KKT, FD, multistart) | done |
| Novelty search: OpenAlex + web, 2018–2026, ~15 queries (lit/) | done, not found. **Before submission:** Google Scholar check of the citations of Heimisson & Segall 2018 and Hill 2024; check the Stefanou/Tzortzopoulos group's 2024–26 papers |
| Field replication: Otaniemi 2018 (public injection and catalogue data, Kwiatek et al. 2019 Sci Adv), or Basel 2006 | **Done (Cycle 3)** for Otaniemi: E1 PASS, E2 WEAK, E3 FAIL, E6 PASS; gain 1.6 % explained by the θ-law (F10, candidate). A second site (Basel 2006) is recommended before submission so F10 can be VERIFIED |
| Engine demonstration: heterogeneous field via ResSim Green's functions, with forward check in the nonlinear simulator | done (C20, F8: −18.0 % events, model error 0.15 % after one correction) |
| Optimum structure + closed form G(Π, θ) (Cycle 4, F12) | done (numerics ≤ 0.001 %); write Appendix C (3-arc derivation, κ₂) |
| Kernel-class dependence (F13) | exploratory; Cycle 5 tested the kernel class on field data: at Äspö the Dieterich kernel is rejected for HF-dominated AE (F14) → use as a domain-of-validity statement |
| Second site | **Äspö HRL done (Cycle 5, F14–F16): a scope/limitation site, not a positive replication.** A positive second site is still needed (Basel 2006, H8: digitise the hourly log) |
| Full proofs + appendix | **C24: Theorem 1 (Hölder proof), Corollaries 1–2, Propositions 2–3 written in full (manuscript §2).** Appendix C (three-arc derivation) still a sketch |
| Draft with all sections + figures from our own app | **C24: v0.4 done** — P1_manuscript_draft.pdf; figures rendered by ResSim's figure engine (tools/paper_figs.cpp, mkfigs.py); Fig. 5–6 = coupled fault run of the ResSim engine |
| Second stimulation, pre-registered (C25, Cycle 6) | **Otaniemi 2020 OTN-2: E4′ (gain ≤ 3 %) REFUTED (9.1 %), E8′ θ-peak CONFIRMED (0.26), E9 transfer PASS, E6 FAIL.** F10 stays CANDIDATE (post-hoc forecast 6.5 %). Same reservoir → not independent; Utah FORGE (public GDR) is the third-site candidate |
| Appendix C (three-arc derivation, J_c, κ₂) | **C25: written in full** (manuscript Appendix C; GRL SI Text S1) |
| Paper pack (C25) | **v0.5 done:** grl/ (AGU LaTeX letter + SI + cover letter + refs.bib + data availability, compiled with pdflatex), jgr/ (long version PDF + markdown), README.md, release zip. Open: third site, reference titles marked [TITLE TO VERIFY], authors, Zenodo DOI |
| Code archive (Zenodo DOI, MIT licence for the scripts) | at submission |

## Target venues (in order)
1. **Geophysical Research Letters**: short, high-impact theory plus the reconciliation. Fits the 12-unit limit: 3 figures, 1 table.
2. **JGR: Solid Earth**: the long version with poroelastic regimes, field case and simulator.
3. Alternatives: *Geothermics*; *IJGGC* (CCS framing); *Geomech. Energy Environ.*

**Preprint:** ESSOAr or EarthArXiv on the day of submission (accepted by both AGU and Elsevier).

## Authorship and ethics
- The user (and any collaborators) are the authors. The AI cannot be an author; its assistance is disclosed in Methods/Acknowledgements as AGU requires.
- The user must understand and be able to defend every step. The appendix proofs are written for that purpose.

## Spin-off
- **F15 comment / data note** (optional, low effort): stage-level reanalysis of the Äspö cyclic-injection claim (SRL "Data Mine" or Solid Earth short communication).
- **P2 (Cycle 6):** optimal injection against simulated rate-and-state slip (G3 + G5, H7). Is convexity preserved with elastic interaction?

## Risks
- **Prior art (structure):** Stefanou group control papers (2024–26, e.g. arXiv 2604.26802) — tracking controllers, no horizon optimum found; re-check before submission.
- **Prior art:** the convexity may be stated implicitly in control or optimisation papers. If so, reframe around the window artefact, the Π law, the regime-dependent taper and the simulator integration.
- **Reviewers will ask about** magnitudes, aftershocks (ETAS) and stochasticity. Plan: a Gutenberg–Richter exceedance appendix, plus a note that ETAS triggering breaks convexity.

## Timeline (suggested)
| Weeks | Task |
|---|---|
| 1–2 | Field replication (Cycle 3) |
| 3 | Engine demo (C20) |
| 4 | Full draft, then user review |
| 5 | Preprint + GRL submission |
