# P1 paper pack — "Minimum-seismicity injection design is convex" (v0.5, 2026-10-21)

| Path | Content |
|---|---|
| `grl/P1_GRL.tex/.pdf` | GRL Letter, AGU template (agujournal2019.cls), 4 figures, ≈ 11 PU (limit 12), abstract ≤ 150 words, 3 key points, PLS |
| `grl/P1_GRL_SI.tex/.pdf` | Supporting Information: Text S1–S4, Table S1, Figures S1–S3 |
| `grl/COVER_LETTER.md`, `grl/DATA_AVAILABILITY.md`, `grl/refs.bib` | submission material |
| `jgr/` | JGR Solid Earth long version (markdown master + PDF), see jgr/README.md |
| `manuscript.md`, `build_pdf.py` | long-version master and renderer |
| `mkfigs.py`, `fig_data.py`, `figs/` | figure pipeline (ResSim figure engine) + vector PDF/SVG + CSV data |
| `PUBLICATION_PLAN.md` | venue plan and checklist |
| `../../data/`, `../../exp/` | processed data and experiment scripts (reproduction: grl/DATA_AVAILABILITY.md) |

Build GRL: `cd grl && pdflatex P1_GRL && bibtex P1_GRL && pdflatex P1_GRL && pdflatex P1_GRL`.
Honest status: not submission-ready — see jgr/README.md "Status".
