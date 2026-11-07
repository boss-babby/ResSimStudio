# Research Cycle 3 — report

Date: 2026-10-20 · Pre-registration: `CYCLE3_PLAN.md` (written before any experiment) · Scripts: `exp/h1e_*.py`, `exp/h1d_magnitude.py`,
`exp/h2_*`, `exp/h3_mkdeck.py` · Raw numbers: `exp/h1e_results.json`, `exp/h1d_results.json`.

## Summary

| ID | Question | Verdict | Consequence |
|---|---|---|---|
| H1e | Does the convex seismicity-design theory survive field data (Otaniemi 2018)? | **Partly.** Pressure model PASS, seismicity fit good in-sample, temporal forecast of the post-injection tail FAILS; retrospective optimum gains only 1.6 % vs constant rate, although the Π law predicted 13 % | New candidate finding **F10**: the gain is a function of two numbers, G(Π, θ), θ = kernel time / injection time; it peaks at θ ≈ 0.1–1. Otaniemi has θ ≪ 0.1 for the fast pressure modes |
| H1d | Are magnitude-aware objectives convex-equivalent? Does a stress-dependent b-value break it? | D1 **proven**; D2 **no violation found** (numerical, proof open); D3 tapered G-R at Otaniemi | Designer feature: hazard cap / max safe volume (1.09× volume at the same P(M ≥ 2)) |
| H3 | Why are DD / GPU not faster on the user's laptop, with high CPU use? | **Explained and fixed** | thread default = physical cores, polite MPI wait + PASSIVE OpenMP, hybrid DD, lazy CPR, advisor (ENGINES_AND_MATH §6.13) |
| H2 | Does a sub-grid tip give grid-converged HF growth at 10 / 5 / 2.5 m? | **Refuted** as stated | New finding **F11**: the dominant error is host-cell leak-off quantisation; the discrete and sub-grid schemes bracket the converged answer from both sides |

## H1e — Otaniemi 2018 field replication

Data: OTN-3 1-min pumping (Kukkonen et al. 2023, CC BY 4.0) and the GFZ catalogue (Leonhardt et al. 2021, CC BY 4.0). Injected volume 18,433 m³,
shut-in at 1168 h (22 Jul 2018 15:52 UTC; pumping times are UTC+3 in the source and were shifted), 61,163 catalogue events, M_c = −0.95
(max curvature −1.15 + 0.2), 23,854 events above M_c, b = 1.278 ± 0.008 (Aki–Utsu), M_L max 1.90.

| Test | Pre-registered pass | Result | Verdict |
|---|---|---|---|
| E1 rate → WHP linear kernel, fit on 60 %, predict 40 % | RMSE ≤ 10 % of range | 60.8 bar = 6.1 % of 990 bar | **PASS** |
| E2 Dieterich fit in-sample | corr ≥ 0.98, tail ± 30 % | corr 0.998; tail 2756 vs 1839 observed (+50 %); τ_d 17.1 h | **WEAK** |
| E3 temporal hold-out (fit stages 1–3) | ± 30 % at shut-in and end, beat the SI baseline | shut-in +2.7 % (SI baseline −0.9 %), end +28.7 % | **FAIL** (worse than the SI baseline at shut-in) |
| E4 retrospective optimum (same volume, 0 ≤ q ≤ q_max, full window) | descriptive | actual 23,828 / constant 21,546 / optimal 21,211 expected events; Π = 3.06 → Π law predicts 13.2 % gain, realised 1.6 % vs constant (11 % vs actual) | the operator's schedule cost +10.6 % vs constant |
| E5 window artefact | descriptive | actual and constant keep their order; only the two optima swap (each wins its own window) | no practical flip |
| E6 robustness (hold-out parameters, full-data evaluation) | keep ≥ 50 % of the gain | 70 % | **PASS** |

Exploratory (not pre-registered, so stated as hypotheses):
* **E7** — splitting the kernel: with the fast modes only (T ≤ 30 h) the actual schedule is +30 % vs constant but the optimum gains only 1.8 %; with the slow
  modes only (T ≥ 100 h) the optimum gains 4.0 %.
* **E8** — rescaling the fitted kernel time by a factor gives the gain as a function of θ = T_k / T_inj:

| θ | 0.0009 | 0.0086 | 0.026 | 0.086 | 0.26 | 0.86 | 2.57 |
|---|---|---|---|---|---|---|---|
| optimal vs constant | 0.0 % | −0.6 % | −2.0 % | −5.5 % | −9.2 % | −7.0 % | −3.9 % |
| actual vs constant | +91 % | +51 % | +28 % | +8 % | +0.2 % | −0.3 % | −0.3 % |

Interpretation (**F10**, candidate). The Cycle-1 Π law assumed a response time comparable with the injection time. With θ → 0 the stress follows the
rate instantly: every schedule with the same rate-duration distribution gives the same count, and the constant rate is already optimal (Jensen gap of a
convex function of a constant argument). With θ → ∞ the stress integrates the rate and the schedule does not matter either. The design benefit is
G(Π, θ) with a maximum near θ ≈ 0.1–1; the *penalty for a bad schedule* grows as θ → 0 (+91 %). Practical rule: a fast-responding reservoir needs a
**smooth** schedule (avoid spikes) but gains nothing from optimisation; a reservoir with θ ≈ 0.1–1 gains up to ~10 % from an optimised schedule.
Earlier claim refuted: "long kernel memory explains the small Otaniemi gain" — the cause is θ ≪ 0.1 for the modes that carry the stress.
Literature check (OpenAlex, Oct 2026): no paper states a two-parameter (Π, θ) law; related are the GRL 2026 cycled-injection lab study
(doi 10.1029/2025gl119289: cadence changes M_max, not total moment) and the Basel model comparison (Phil Trans A 2024, doi 10.1098/rsta.2023.0179).

D3 (magnitudes at Otaniemi): G-R with b = 1.28 predicts 4.0 events with M_L ≥ 2 and 17.6 with M_L ≥ 1.5; observed 0 and 17. P(0 | 4.04) = 0.018, so the
tail is tapered above M_L ≈ 1.5–2 (b is stable 1.24–1.30 for M_c −0.95 … 1.0). At the actual schedule's hazard P(M_L ≥ 2) = 0.982 the optimal schedule could
inject 20,036 m³ (1.09×).

## H1d — magnitude objectives

* D1 (proof): with independent Poisson patches P(M ≥ m) = 1 − exp(−Σ_k N_k 10^{−b_k(m − M_c,k)}), a monotone function of a non-negative combination of convex
  N_k → quasi-convex, same minimiser as the weighted count; the chance constraint P ≤ α is the convex constraint Σ w_k N_k ≤ −ln(1 − α). Numerical check:
  2000 random midpoint pairs, 0 violations.
* D2 (stress-dependent b = b0 − βΔS): 24 configurations (t_a 100 / 0.05 yr, β 0 … 0.5 /MPa, Δm 1–3) and 12 synthetic extremes (γΔS_max up to 30), 2000
  pairs each: **0 violations**. The composition proof fails (the exponent 1 + γa exceeds 1), so D2 stays a numerical result.

## H3 — laptop performance (user report: DD and GPU not faster, high CPU)

Measured on this 2-core box; wall times. SPE9LIKE (9k cells): 1 thread 9.33 s, OpenMP 2 6.83 s, old DD 2 ranks 10.58 s, hybrid DD + PASSIVE 7.27 s.
CH15 (54k cells): 1 thread 48.6 s, OpenMP 2 31.4 s (1.55×), OpenMP 2 + level ILU 34.2 s, old DD 48.4 s (1324 vs 949 linear iterations), hybrid DD 36.0 s.
Causes: busy-polling MPI ranks and actively spinning OpenMP threads (Task Manager shows ~100 % without progress; hybrid DD with active wait was 3.5×
slower), logical-CPU thread default (hyper-threads), CPR rebuilt every Newton iteration, DD's block-Jacobi coupling (+40 % iterations), no device AMG on the
GPU path. Fixes in the engine/app (C21). Expected on the user's 4-core laptop: OpenMP ≈ 2.4× (Amdahl from the measured serial fraction), DD below that.
Not a research finding (known HPC practice); recorded for the engineering log.

## H2 — sub-grid HF tip

T1 baseline (discrete column opening) and T2 (continuous tip, leak-off area ∝ open fraction, K_I root solved in the step). Time for the + wing to reach
60 / 90 m (HF_dx decks, 200 m trace, K_IC 1 MPa√m, σ_min 340 bar):

| scheme | dx 10 m | dx 5 m | dx 2.5 m |
|---|---|---|---|
| discrete | 22.0 / 56.0 d | 24.0 / 56.0 d | 25.0 / 58.5 d |
| sub-grid | 29.5 / 70.2 d | 28.0 / 73.6 d | 26.7 / 60.9 d |

Pre-registered pass (spread over the grids drops ≥ 3×) is **not met**. Mechanism (debug trace at the 35 m piece boundary): one step after the tip enters a
fresh host cell (φ = 0.02) the net pressure drops 2.64 → 2.00 bar because that cell is still at the initial pressure (Δp ≈ 42 bar vs ≈ 1.5 bar behind the
tip): leak-off per area is ~30× larger. Growth therefore stalls at each host-cell boundary and then runs. The discrete scheme makes the opposite error
(the whole column opens at once). Earlier attempts recorded as dead ends: K_I from the tip piece's pressure (saw-tooth), from the outermost full column
(jumps), φ-linear fracture conductance (acts as a step because k_f = w²/12 ≫ k_host).

**F11** (candidate): in EDFM-coupled LEFM growth with leak-off, the grid error of the growth time is dominated by the host-cell leak-off quantisation, not by
the tip location; a continuous tip alone does not reduce it, but the discrete (column) and continuous (area-weighted) schemes converge from opposite sides
and bracket the converged growth time (t60 at dx 2.5: 25.0 … 26.7 d). Follow-up (C22-a): a sub-grid leak-off model (Carter-type fresh-rock correction for the
tip host cell); that should make the sub-grid scheme converge at coarse dx.

## Engineering outputs of Cycle 3 (in the app)
* Thread default, polite MPI waiting, hybrid DD, lazy CPR, HPC performance advisor (Run → HPC).
* `HFGROW … SUBGRID` (research option with the bracket note in the UI).
* G1-c flow ↔ contact coupling (Tools → Geomechanics, `--geocouple`), not part of the hypotheses but of the G-plan.

## Next cycle (4) candidates
* H2b: sub-grid leak-off — pre-registered pass: spread of t60 over dx 10/5/2.5 below 1 d.
* H1f: θ-law — closed form of G(Π, θ) for one exponential kernel; check against Basel 2006 (second field case) to move F10 to VERIFIED.
* H5: coupled G1-c in the seismicity designer (slip-based instead of Coulomb proxies).
