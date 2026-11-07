# P4 dossier: grid-converged leak-off for embedded hydraulic fractures (NOT a manuscript yet)

**Status (C50, 2026-11-04): deferred.** The evidence is internally consistent but not publishable as it stands. This dossier records
what exists, what a referee would demand, and the minimum work needed to change that verdict.

## Claim (candidate)
In EDFM-coupled LEFM fracture growth with leak-off, the grid error of the growth time is dominated by host-cell leak-off
quantisation (F11). Treating the exposed fraction of the tip host cell as its own sub-cell storage ("sub-grid leak-off", SGL) removes
the grid dependence (F21), including on a dipping corner-point grid (F24). For dry gas a compressibility-aware sub-cell (SGL-C,
engine option `SUBGRID LEAKC`) removes the over-correction of plain SGL (F30 → F32).

## Evidence in hand
| Finding | Test | Result | Files |
|---|---|---|---|
| F11 | grid bracket, 2-D EDFM, PKN-type K_I | SG converges from above, column opening from below | `exp/h2_driver.cpp`, `exp/h2_mkdeck.py` |
| F21 | pre-registered H2b | t60 spread 0.85 d (pass < 1 d) vs 6.6 d | `exp/h2b_results.txt` |
| F24 | pre-registered H2c (corner-point, 2° dip) | SUPPORTED; gas variant non-discriminating | `exp/h2c_mkdeck.py`, `exp/h2c_results.txt` |
| F30 | H2d / H2d′ (gas) | inconclusive / over-correction at t90 | `exp/h2d*.log`, `exp/h2d_run.py` |
| F32 | pre-registered H2e | SUPPORTED (t90 spread below SG) | `exp/h2e.log`, `exp/h2e_results.json` |

## Why it is not ready (referee view)
1. **Novelty is contested.** Embedded/EFM fracture models with leak-off and sub-cell treatments exist: Norbeck & Horne (2015), and
   Li et al. (2020, JGR, 10.1029/2020JB020115). Our sub-cell storage may be a variant of existing near-tip leak-off corrections. A
   systematic literature comparison is needed, covering Carter leak-off in P3D/PKN codes and Vinsome–Westerveld-type semi-analytic
   leak-off.
2. **No analytic benchmark.** Grid convergence shows self-consistency, not correctness. Required: KGD or PKN with Carter leak-off
   (Nordgren; Economides & Nolte), in toughness- and leak-off-dominated limits.
3. **No independent reference.** Required: a fine-grid explicit-fracture (DFM) reference or an independent code.
4. **One deck family**, single-phase water plus one dry-gas case; no multiphase leak-off, no 3-D.

## B3 literature check: partly done (C50, 2026-11-04)
* **The problem is known.** Norbeck & Horne (2015, embedded fracture modelling framework, Stanford/Springer) report in a grid-refinement
  study that coarse EFM grids underestimate leak-off and give longer fractures. Their solutions converge toward a *1-D leak-off
  approximation* (a sub-grid, Carter-type analytic leak-off), and they benchmark against KGD and Carter-type semi-analytic lengths
  (−1 to −4.4 %). Li et al. (2020, JGR) also note that coarse matrix grids underestimate leak-off.
* **Our remedy differs** (partial tip-cell storage rather than a 1-D analytic leak-off profile), but nothing shows it is better.
  **The decisive comparison is now SGL vs the Norbeck–Horne 1-D leak-off approximation** on the same deck, in addition to B1.
* Novelty rating lowered: "variant of a known remedy". The most likely outcome is a technical note in ENGINES_AND_MATH, not a paper.

## Minimum work to reach "submittable" (about 1–2 research cycles)
* B1: Carter-leak-off PKN benchmark (analytic) in both limits, at three grids (pre-register the tolerance).
* B2: fine-grid explicit reference on one water deck (dx 0.625 m).
* B3: literature table (started above); add the SGL vs 1-D leak-off comparison (B4).
* Decision rule: if B1 fails or B3 shows SGL ≡ an existing method, demote it to a technical note in ENGINES_AND_MATH
  (`SUBGRID LEAKC` documentation) and stop.

Target if it passes: Computational Geosciences (fallback SPE Journal).
