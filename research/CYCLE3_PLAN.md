# Research Cycle 3: pre-registration (written before any Cycle-3 experiment was run)

Date: 2026-10-20. Order set by the user: research first, then the G-series, plus the user's performance report (4-core DD and GPU not faster than serial).

| ID | Question | Decides |
|---|---|---|
| H1e | Does the convex seismicity-design theory (F1–F4) survive contact with **field data** (Otaniemi 2018)? | P1 publication (field replication was the main TODO) |
| H1d | Magnitude-aware objectives: are Gutenberg–Richter exceedance objectives / chance constraints convex-equivalent? Does a **stress-dependent b-value** break convexity? | P1 appendix + designer feature (hazard cap, max safe volume) |
| H3 | Why does MPI DD (4 ranks) or the GPU give ≈ serial speed on the user's laptop while the CPU shows high usage? | engineering priority for HPC |
| H2 | Does a sub-grid tip position in EDFM HFGROW give grid-converged propagation (10 / 5 / 2.5 m)? | HF quality |

## H1e — Otaniemi 2018 field replication
**Data** (both CC BY 4.0): OTN-3 1-min pumping data (Kukkonen et al. 2023, Mendeley 10.17632/ycdcw2rdfj.1); catalogue (Leonhardt et al. 2021, GFZ 10.5880/GFZ.4.2.2021.001). Processing: `exp/h1e_fetch.py`. 49 days of injection (shut-in 22 July 2018 15:52 UTC), 63 days post-injection monitoring, ~18,000 m³.

**Model class (stated up front):** one aggregate fault population with Dieterich (1994) rate-and-state seismicity. The Coulomb stress proxy is proportional to a reservoir pressure built as a linear convolution of the injection rate. The kernel is a non-negative sum of exponentials fitted to the wellhead pressure; its instantaneous terms (tubing friction, skin) are excluded from the stress. An optional diffusion delay τ_d to the seismogenic volume is allowed. Parameters: r, P* = Aσ̄/α (bar), t_a, τ_d. Poisson likelihood of hourly counts with M_L ≥ M_c.

**Pre-registered criteria**

| # | Test | Pass | Kill / weaken |
|---|---|---|---|
| E1 | Linear rate → pressure kernel: fitted on the first 60 % of the injection, predicts the last 40 % of the WHP | RMSE ≤ 10 % of the WHP range | > 20 % → the pressure response is non-stationary (permeability creation); the design must then be adaptive (SCP) |
| E2 | Dieterich fit, in-sample | correlation of cumulative counts ≥ 0.98; post-shut-in decay reproduced (count in the 63-day tail within ±30 %) | tail off by > 2× → the model misses post-injection physics (e.g. poroelastic or aseismic slip) |
| E3 | Temporal hold-out: fit to stages 1–3 (before 1 July), forecast stages 4–5 + tail | cumulative count at shut-in and at the end within ±30 %; better than the seismogenic-index baseline N = Σ·V at shut-in | worse than the baseline → the field model is not predictive enough to design with |
| E4 | Retrospective optimum (same volume and window, 0 ≤ q ≤ q_max observed), events counted over the **full** monitoring window | report Π, the actual / constant / optimal expected counts and the Π-law prediction of the gain (F3) vs the realised gain | — (descriptive) |
| E5 | Window artefact (F2) on field data: rankings with counting stopped at shut-in vs over the full window | report whether the ranking flips | — |
| E6 | Robustness: the optimum computed with the hold-out parameters, evaluated with the full-data parameters | keeps ≥ 50 % of the full-information gain | < 0 → designs are not robust to calibration error |

## H1d — magnitudes
D1. Proof: with independent Poisson patches, P(M ≥ m in window) = 1 − exp(−Σ_k N_k 10^{−b_k (m − M_c,k)}). This is monotone in a non-negative weighted sum of convex functions, so it is quasi-convex with the same minimiser as the weighted sum. A chance constraint P(M ≥ m*) ≤ α is a convex constraint. "Maximum volume under a hazard cap" is therefore a convex programme. Check numerically (midpoint convexity of the weighted sum, KKT).
D2. Stress-dependent b-value (b = b0 − β ΔS, Scholz 2015 / laboratory): test convexity of the exceedance rate with 2,000 random midpoint checks over β ∈ [0, 0.5]/MPa. Pass = no violation (convexity survives); otherwise report the β threshold and whether multistart optima differ.
D3. Otaniemi: b-value (Aki–Utsu MLE), P(M_L ≥ 2) for the actual / constant / optimal schedules, and the maximum volume that keeps P(M_L ≥ 2) ≤ the actual schedule's value.

## H3 — parallel performance (user report)
P1. Time breakdown per Newton phase (property evaluation, assembly, linear solve, well update), serial, on SPE9LIKE and a larger model.
P2. Strong scaling of OpenMP threads (1, 2) and MPI DD (1, 2 ranks) on this 2-core box. Amdahl extrapolation to 4 cores from the measured serial fraction.
P3. Busy-wait check: CPU use of the ranks that are waiting (MPI polling shows as 100 % CPU in Task Manager).
P4. GPU path: count of host↔device transfers per linear iteration and the measured break-even model size.
Pass = a quantified explanation that predicts the user's observation, plus a fix list ranked by expected speed-up. Kill = an unexplained slowdown → instrument more.

## H2 — sub-grid tip (engine)
T1. Baseline: HFGROW propagation time to reach 100 m at grid sizes 10 / 5 / 2.5 m (piece-quantised).
T2. Sub-grid tip: continuous tip coordinate, partial-cell fracture transmissibility ∝ open fraction, tip criterion evaluated at the tip coordinate.
Pass = the spread of times over the 3 grids drops by ≥ 3× vs baseline and the 2.5 m result is within 10 % of the PKN/KGD semi-analytic reference.
