# Research Cycle 6 — report (2026-10-21)

Pre-registration: `CYCLE6_PLAN.md` (written before any fit on the new data). Code: `exp/h8_fetch.py` (data),
`exp/h1e_otaniemi.py` with `SITE=otn2` (same pipeline as Cycle 3, plus the new test E9), `exp/h8_theta.py` (post hoc).
Outputs: `exp/h8_results.json`, `exp/h8_curves.json`, `exp/h8_theta.json`, `exp/h8_run.log`. Runtime: 6 s.

## Site and data
Otaniemi 2020, OTN-2 well (St1 Deep Heat). Basel 2006 (H8) was replaced because its hourly log is not public (see
the plan). Hydraulics: Kukkonen et al. 2023 (Mendeley doi 10.17632/ycdcw2rdfj.1, OTN2.xlsx, 1-min, local time UTC+3).
Catalogue: Kwiatek et al. 2022 (GFZ doi 10.5880/GFZ.4.2.2022.001).

| quantity | value |
|---|---|
| injected volume | 2,884 m³ (paper: 2,875 m³) |
| shut-in | 387.3 h after 2020-05-05 00:00 UTC (16.1 d) |
| max rate / max WHP | 1.05 m³/min / 716 bar |
| events after T0 / above M_c | 6,121 / 2,846 (M_c = −1.45, b = 1.39 ± 0.03, M_L max 1.09) |
| operator pauses ≥ 3 h | 8 pauses, incl. one of 116 h (49–165 h) |

## Pre-registered tests

| test | criterion | result | prediction | verdict |
|---|---|---|---|---|
| E1′ kernel | hold-out RMSE ≤ 10 % of WHP range | 81.8 bar = 11.4 % of 716 bar | PASS | **WEAK** (prediction missed narrowly) |
| E2′ Dieterich fit | corr ≥ 0.98, tail ±30 % | corr 0.9984; tail 351 vs 325 (+8 %) | WEAK | **PASS** (better than predicted) |
| E3′ temporal hold-out | shut-in and end ±30 %, beat SI | shut-in −31 %, end −38 %, SI −42 % | FAIL (end) | **FAIL** — the hold-out fit hits the P* bound (0.5 bar): the first 60 % of the time contain 2 days of injection and the 116-h pause, so P* is not identifiable |
| E9 cross-stimulation transfer | 2018 P*, t_a, τ_d fixed; corr ≥ 0.95, shut-in ±30 % | corr 0.985; shut-in −23 % | PASS | **PASS**; tail over-predicted (910 vs 325, 2.8×) — same bias as 2018 (1.5×), larger |
| E4′ optimal vs constant | gain ≤ 3 % | **9.1 %** (2,078 vs 2,286); Π 3.43, Π-law 15.6 % | ≤ 3 % | **REFUTED** |
| E4′ operator vs constant | within ±10 % | **+24 %** (2,834 vs 2,286) | ±10 % | **REFUTED** |
| E8′ θ-scan | gain peaks for θ ∈ [0.05, 1], vanishes at both ends | 0.05 / 0.4 / 2.4 / 7.6 / **14.0** / 11.5 / 5.8 / 3.2 % at θ = 0.0026 … 7.7; peak at θ = 0.26 | peak in [0.05, 1] | **CONFIRMED** |

Verdict rule for F10 (fixed in the plan): VERIFIED only if E4′ and E8′ both meet their predictions → **F10 stays
CANDIDATE**.

## Why E4′ failed: the input to the law was wrong, not (apparently) the law
The prediction "≤ 3 %" assumed θ ≪ 0.1 as in 2018. That assumption was wrong. The 2020 injection lasted 16 days instead of
49, so the same crustal response times give a θ three times larger. The fitted 2020 kernel carries 73 % of its
steady-state pressure in components with T = 100–300 h, i.e. θ = 0.26–0.77, which is the peak of the θ-curve.

**Post-hoc θ-law forecast** (`exp/h8_theta.py`; exploratory, *not* pre-registered): weight the 2018 θ-scan gain
G₂₀₁₈(θ_j) by the 2020 kernel's pressure shares a_j. The result is **6.5 %**, against **9.1 % realised** and **15.6 %**
from the Π-law. The θ-law computed from the 2018 data alone predicts the 2020 gain to within a factor of 1.4. The Π-law
over-predicts by 1.7× (in 2018 it was 8×). One caveat: Π is 3.43 in 2020 against 3.06 in 2018, and no Π correction was
applied. Because the forecast was made after seeing the result, it is reported as consistent, not as a confirmation.

## Fairness of the comparator (exploratory, E4b)
The operator paused 8 times (≥ 3 h), 51 % of the time. Two comparisons when the same pauses are imposed:
- **Constant rate inside the active windows: 3,816 events.** This is *worse* than the operator's schedule (2,834):
  compressing the volume into fewer hours raises Π.
- **Optimum with the pauses: 2,747 events, 3.1 % better than the operator's schedule.**

So the 9.1 % / 26.7 % gains assume the pauses could be avoided. If the pauses are operationally fixed, the 2020 schedule
was within 3 % of optimal. Both numbers go into the paper.

## Robustness (E6) failed
The design made from the hold-out parameters was evaluated with the full parameters and produced **6,723 events**,
2.9× the constant-rate count. The cause is the degenerate hold-out fit (P* at its 0.5-bar bound), which makes the
designer believe the stress response is extremely sharp. Lesson for practice: **design only with an identifiable
P\***. In the app's designer this needs a guard: refuse a design or warn when a parameter sits on a bound. This is now
ROADMAP item C25-r1.

## Window artefact (E5) replicated
The ranking at shut-in (opt_shutin < opt_full < constant < actual) differs from the ranking over the full window
(opt_full < opt_shutin < constant < actual). This is the F2 flip, now seen in a second field data set.

## Novelty re-check (OpenAlex, 2026-10-21, publications since 2025-06)
Four queries: convex injection schedule; optimal injection rate rate-and-state; protocol optimisation Dieterich;
bang-bang injection control. Newly found:
- Hannover thesis 2026 (doi 10.15488/22844): spring-slider RSF models, cyclic and depressurisation protocols.
- EarthArXiv review 2026 (doi 10.31223/x5p75x): forecasting and control.
- JGR 2026 Bayesian stress inversion (doi 10.1029/2026jb033922).

None of these states convexity, a Pontryagin bang–singular–bang structure, a closed-form gain, or a θ-dependence. The
status is unchanged: **not found**. Google Scholar was not reachable from the sandbox, so the Scholar check is still
open for the authors.

## Consequences
- F10: evidence extended (second stimulation, E8′ shape confirmed, post-hoc forecast 6.5 vs 9.1 %); status stays
  CANDIDATE by the pre-registered rule. The paper frames it as "supported by two stimulations, one pre-registered
  prediction missed because θ was mis-estimated."
- New finding F19: the operator-pause comparator and the identifiability guard.
- P1 now has a second field data set. The paper's main claims rest on theory (F1, F2, F4, F12, all VERIFIED) plus field
  consistency, and P1 is ready to be packaged for GRL (letter) and JGR (long version).
