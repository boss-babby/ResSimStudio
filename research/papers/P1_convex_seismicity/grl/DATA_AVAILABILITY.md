# Open Research / data availability (GRL and JGR versions)

| Item | Source | Licence | Used in |
|---|---|---|---|
| Otaniemi 2018 catalogue | GFZ Data Services, doi 10.5880/GFZ.4.2.2021.001 | CC BY 4.0 | Fig. 3 (GRL) / Fig. 3 (JGR) |
| Otaniemi 2020 catalogue | GFZ Data Services, doi 10.5880/GFZ.4.2.2022.001 | CC BY 4.0 | Fig. 2 / Fig. 4 |
| Otaniemi hydraulics OTN-2, OTN-3 | Mendeley Data, doi 10.17632/ycdcw2rdfj.1 | CC BY 4.0 | Figs. 2–3 |
| Äspö HRL HF data | GFZ, doi 10.5880/GFZ.2.6.2023.004 (ESSD 16:295, 2024) | CC BY 4.0 | Fig. S2 / Fig. 5 |
| Processed data | `research/data/` (otaniemi, otaniemi2020) and the figure CSVs in `figs/*.csv` | CC BY 4.0 | all |
| Scripts | `research/exp/*.py` (h1e_otaniemi.py with SITE=otn2, h8_fetch.py, h8_theta.py, h1f_theta_law.py) | MIT/GPL (authors decide) | all |
| Simulator | ResSim source (`engine/src/faultmech.cpp`, `seismicity.cpp`, `rsf.cpp`, `app/src/figure.cpp`, `tools/paper_figs.cpp`), deck `examples/FAULT_INJECTION.DATA` | authors decide | Fig. 4 / Figs. 6–7 |
| Pre-registrations | `research/CYCLE3_PLAN.md`, `CYCLE6_PLAN.md`, `FINDINGS_LOG.md` | — | Text S3 |

Reproduce: `cd research/exp && python3 h8_fetch.py && SITE=otn2 python3 -u h1e_otaniemi.py && python3 h8_theta.py`; figures:
`python3 mkfigs.py <paper_figs>`; simulator figures: `paper_figs faultrun examples/FAULT_INJECTION.DATA figs/fig6_coupled_fault`.
Zenodo DOI: to be minted from the archive `P1_paper_pack.zip` at submission.
