# Publication portfolio: decision record (C50, 2026-11-04)

This file is the single place that decides **what is published, where, in what order, and which finding belongs to which paper**.
It supersedes the scattered "paper plans" in earlier cycle reports. Findings are numbered as in `../FINDINGS_LOG.md` (F1–F40).

## 1. Verdict: three manuscripts, one candidate (P5), one deferred paper (P4) and a technical report (TR1)

| # | Working title | Target (fallback) | Core findings | Status | Honest publishability |
|---|---|---|---|---|---|
| **P1** | Minimum-seismicity injection design is a convex programme | GRL letter + SI (JGR Solid Earth; GJI) | F1–F4, F6, F7, F10, F12, F13, F19 | LaTeX compiled (10 + 3 pp), references Crossref-checked | **Strongest.** Clear general result plus field tests, including one honest refuted prediction. Main risks: a reviewer calls the theorem "elementary" (log of an integral of an exponential of a linear map), and the θ-law rests on two stimulations at one site |
| **P2** | When does fault-permeability feedback break convex injection design? | JGR Solid Earth (GJI; Geomech. Energy Environ.) | F25–F27, F31, F35, F37, F38, F20/F22/F28/F33/F36 (RSF part), F8 | Markdown v1.0, PDF 11 pp, 5 figures, 2 tables | **Moderate.** Pre-registered, quantitative threshold and decision rule. Weaknesses: 2-D synthetic, one simulator, explicit coupling (the lag inflates the water amplitude about 3.5×). Expect a major revision asking for 3-D or an independent code |
| **P3** | Hydraulic energy, not Dieterich pressure, controls AE at Äspö: a reanalysis of the continuous catalogue | GJI (Solid Earth; SRL) | F39, F40 (with F14–F16 as superseded history) | Markdown v1.0, PDF 10 pp, 5 figures, 2 tables | **Moderate to good.** Public data, pre-registered, strong ΔAIC with bootstrap, data defect found. Risk: "reanalysis, incremental". The HF3 result is deliberately not overclaimed (p 0.04–0.31) |
| P4 | Grid-converged leak-off for embedded hydraulic fractures | (Comput. Geosci.; SPE J.) | F11, F21, F24, F30, F32 | **Not ready**: dossier only (`P4_leakoff/DOSSIER.md`) | Novelty weakened: Norbeck & Horne (2015) already show the coarse-grid leak-off error and a converged 1-D leak-off remedy. Needs a PKN/Carter benchmark plus SGL vs 1-D leak-off; likely demoted to a technical note |
| P5 (candidate) | Convex, memory-robust scheduling of cyclic subsurface storage | (GRL; IJGGC / Applied Energy) | F41, F42 | Cycle 13 done, 1 figure (`P5_storage/figs`) | **Promising, not ready.** One 2-D model; needs the gas case, a second geometry, multi-cycle runs and an equilibrated state (NEW_FIELDS) |
| TR1 | Technical report: OPM regression references and well-reference-depth semantics | Zenodo report / OPM issue | F29, F34 | **Written**: `TR1_opm_reference_drift/TECH_REPORT.md` + PDF (2 pp) | Useful to the community, not a journal paper |
| — | Notes, not papers | — | F5 (P1 SI, or drop), F9 and F17 (engineering notes in ENGINES_AND_MATH), F18 (one case: P1 sentence only), F23 (P3 caveat) | — | — |

**Why not more papers.** We considered splitting P2 (CO₂ convexity / permeability feedback / RSF surrogate). Each part alone is a
thin synthetic result, and the split would be salami slicing. The RSF part (F20, F22, F28, F33, F36) stays in P2 as the
"what to do when convexity fails" section. Merging P1 and P2 was also considered and rejected: P1's claim is a theorem with field tests
(letter format); P2 is a simulation study of where the theorem's assumptions break.

## 2. One home per finding (duplicate-publication matrix)

| Finding | Primary home | May be mentioned (one sentence, citing the primary) |
|---|---|---|
| F1–F4, F6, F7, F10, F12, F13, F19 | P1 | P2 §2 (theory recap) |
| F18 (+36 %, one case) | P1 (§4, labelled one test case) | P2 introduction |
| F25–F27, F31, F35, F37, F38, F8 | P2 | P1 §4 (companion study, in preparation) |
| F20, F22, F28, F33, F36 | P2 (§4.5–4.6) | — |
| F39, F40 | P3 | P1 §3 (scope sentence only, **no numbers**) |
| F14–F16 | superseded by F39/F40; history in P3 §2 | — |
| F23 (FORGE) | P3 limitation | — |
| F11, F21, F24, F30, F32 | P4 | — |

**C50 changes made to enforce this:** the Äspö paragraph and Figure S2 were removed from the P1 letter (replaced by a scope
sentence; Figure S3 renumbered S2). The P1 JGR long version v0.6 and `manuscript.md` are marked SUPERSEDED, because they would
duplicate P2 and P3. P1 no longer claims that convexity "holds within a sequential-convex iteration" under two-way coupling, which P2
showed is false above m₁ ≈ 1.7.

## 3. Critical self-review (as a referee would read them)

### P1 (GRL)
* **Novelty.** Convexity of the Dieterich expected count follows from the log-sum-exp structure. The contribution is the
  *consequences*: window-dependent rankings (F2), dose invariance and the Π/θ law (F3, F10), bang–singular–bang optima (F12) and
  identifiability (F19), plus field tests. The cover letter must state this plainly rather than oversell the theorem.
* **Closest prior work (re-checked C50).** Kim, Gutiérrez-Oribio, Stefanou, Acosta & Avouac (2025, Geothermics 132:103396) optimise
  the *same* Otaniemi stimulation with a tracking controller and a scenario search. This paper was missing and is now cited and
  differentiated: no global-optimality certificate, but a consistent conclusion. Also added: Gutiérrez-Oribio & Stefanou 2024,
  Walter et al. 2026 (EarthArXiv).
* **Weak points to pre-empt.** (i) θ-law: two stimulations, one site (CANDIDATE). A Basel or Pohang test would make it VERIFIED, but
  those injection logs are not public. Say so. (ii) The Otaniemi 2018 gain is only 1.6 %; the paper's honesty about small gains is a
  strength, not a weakness. (iii) The pre-registrations are internal files, not registry-stamped (now disclosed in SI Text S3).
* **Verdict.** Submit. GRL desk rejection is a real possibility (GRL wants broad interest); fallback is JGR Solid Earth with the
  same SI promoted to appendices. No new computations are needed.

### P2 (JGR Solid Earth)
* **Strengths.** Two pre-registered tests supported (H10, H10g), one refuted with controls (H10e); a dose–response threshold with a CI
  (5 % violations at m₁ 1.73 [1.53, 1.93]); a practical decision rule; lag sensitivity reported.
* **Weak points.** 2-D section; one fault; one simulator; explicit coupling; overdispersed pooled fit (φ 1.9, partly exploratory);
  geometry term p = 0.06. A referee will ask: "Does this hold in 3-D, and with implicit coupling?" Our answer is in §7 (open items),
  but it remains open.
* **Verdict.** Submit after the P1 preprint is public (P2 builds on P1's theorem). Expect a major revision. If JGR declines, go to
  Geomechanics for Energy and the Environment, where simulator-method papers fit better.

### P3 (GJI)
* **Strengths.** Public data; a data defect found and documented (F39); pre-registered; ΔAIC 545 with the stage bootstrap never below
  70; robust across completeness thresholds; the post-shut-in deficit has z = −27.
* **Weak points.** 31 stages; per-test factors absorb much of the variance; HF4 is a different rock. The rate exponent depends on the
  completeness threshold (0.1–2.3) and must be read as qualitative. The physical interpretation (tip-driven AE) is plausible, not proven.
* **Ethics / courtesy.** Report the `origin_time` defect to the data authors (GFZ) **before** submission. A draft note is in
  `P3_aspo_kernels/DATA_DEFECT_NOTE.md`.
* **Verdict.** The quickest to submit; independent of P1/P2.

## 4. Order and dependencies

1. Notify GFZ of the catalogue defect (P3 courtesy) and wait about 2 weeks for a reply. Meanwhile, do the human steps below.
2. Post P1 and P3 preprints on EarthArXiv together. Submit P1 → GRL and P3 → GJI.
3. Archive code, scripts, pre-registration files and figure CSVs on Zenodo (one DOI per paper or one concept DOI). Fill the
   "[DOI to be added]" placeholders.
4. Submit P2 → JGR Solid Earth, citing the P1 preprint.
5. P4: only after the KGD/Carter benchmark and a fine-grid reference (see dossier).

## 5. Steps only the human research team can do

* **Authorship.** Decide authors, affiliations and ORCIDs. **AI tools cannot be authors**; AI use is disclosed in all three drafts and
  must be adapted to AGU/GJI policy wording.
* Read every manuscript line by line and take responsibility for each claim. The drafts are written to be checkable: every number
  traces to a log or JSON file named in the text.
* Choose the code licence (the ResSim engine uses OPM-derived keyword data under GPL-3.0, so a GPL-compatible licence is needed).
* Upload to Zenodo, submit, and handle correspondence.

## 6. File map

| Paper | Source | Built output | Figures |
|---|---|---|---|
| P1 | `P1_convex_seismicity/grl/P1_GRL.tex`, `P1_GRL_SI.tex`, `refs.bib` | `P1_GRL.pdf` (10 pp), `P1_GRL_SI.pdf` (3 pp) | `grl/figs/Figure1–4, S1–S2` (from `../mkfigs.py`) |
| P2 | `P2_coupled_convexity/manuscript.md` | `P2_manuscript_v1.0.pdf` (11 pp) | `figs/fig1–fig5` (from `scripts/mkfigs.py`, `scripts/fig_dose.py`) |
| P3 | `P3_aspo_kernels/manuscript.md` | `P3_manuscript_v1.0.pdf` (10 pp) | `figs/fig1–fig5` (from `scripts/mkfigs.py`) |
| P4 | `P4_leakoff/DOSSIER.md` (B3 literature check: problem known, Norbeck & Horne 2015) | — | — |
| P5 | `P5_storage/` (candidate; CYCLE13_PLAN results) | — | `figs/fig_h12_storage` (`mkfig_h12.py`) |
| TR1 | `TR1_opm_reference_drift/TECH_REPORT.md` | `TR1_tech_report_v1.0.pdf` | — |

Rebuild: P1 `pdflatex P1_GRL && bibtex P1_GRL && pdflatex P1_GRL && pdflatex P1_GRL` (needs TeX Live); P2/P3
`python3 build_md_pdf.py <md> <pdf> "<footer>"` from this folder.

## 7. Reference audit (C50)

All DOIs in P1 `refs.bib` (20 entries; 19 DOIs plus one arXiv) and in the P2/P3 reference lists were resolved against Crossref/DataCite with title matching.
Errors found and fixed:
* P1 `Zang2018` pointed to the wrong paper: the DOI belonged to Hofmann et al. 2018 (CSS, Geothermal Energy). It is now Zang et al.
  2018, Rock Mech. Rock Eng. 52:475–493.
* The P3 draft had a wrong author list for Niemz et al. 2021.
* "and others" author lists were expanded.
* Zang 2017 / Kwiatek 2018 page ranges were added.

**Rule kept: no reference enters a manuscript without a resolved DOI or an archive identifier.**
