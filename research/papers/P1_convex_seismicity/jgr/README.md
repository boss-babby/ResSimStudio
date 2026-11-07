**STATUS C50: SUPERSEDED, not to be submitted.** The GRL letter + SI (`../grl`) is the P1 submission; the Äspö part is now P3; two-way coupling is now P2. A JGR long version would duplicate P2/P3.

# P1 — JGR Solid Earth long version (v0.6)

* `P1_JGR_long_v0.6.md` — master text (long version; `../manuscript.md` is the shorter GRL text): full theory (Theorem 1, Corollaries, Propositions 2–3),
  synthetic study, Otaniemi 2018 + pre-registered Otaniemi 2020, Äspö scope test, simulator integration, Appendix C (θ-law).
* (PDF deleted in C50: superseded) — was rendered by `../build_pdf.py` (7 figures, all produced by ResSim Studio's figure engine via `../mkfigs.py`).
* Figures: `../figs/fig1 … fig6*.pdf|svg|csv` (vector PDF + data CSV for every panel).
* LaTeX: use `../grl/agujournal2019.cls` with `\journalname{Journal of Geophysical Research: Solid Earth}`; the GRL source
  `../grl/P1_GRL.tex` and `../grl/refs.bib` are the skeleton — sections 2–3 and Appendix C of the markdown go in unchanged.
  JGR has no length limit; the long version keeps all derivations in the main text.
* Status: pre-submission. Before submitting: third independent site, reference check (all refs.bib titles verified via Crossref / arXiv in C47),
  novelty re-check (OpenAlex 2026-10-21: no paper states convexity), author list, Zenodo DOI.
