# Research Cycle 4 — report

Date: 2026-10-20. Pre-registration: `CYCLE4_PLAN.md`, written before any experiment. Scripts: `exp/h1f_theta_law.py`
(log `exp/h1f.log`, data `exp/h1f_results.json`, figure `exp/h1f_fig.png`) and `exp/h1g_kernel.py` (exploratory, `exp/h1g_results.json`).
Engine integration: `engine/src/seismicity.cpp` (`seis::scheduleGain`, design potential, last-interval bound) and tests `tests/c22_tests.inc` (`-c22`).

## Summary

1. **The pre-registered optimum structure is wrong: F1 fails as registered.** For one exponential stress kernel the event-count
   optimum is **bang–singular–bang**:
   - a maximum-rate *charging* arc;
   - a constant-stress plateau;
   - a short **maximum-rate terminal burst at shut-in**, which the derivation had missed.

   The costate is not free at the end of injection. The tail integral ∫(e^{Π s(t)}−1) after shut-in depends on s(T_i), and its
   marginal value is smaller than the volume it frees up. So it pays to put the last few percent of the volume in at the highest
   rate, after which the stress decays outside the injection window. With the burst included, the **3-arc closed form matches
   the numerical optimum to ≤ 0.001 %** at every tested point (table F1/F2).
2. **The small-θ law G = θ κ(Π, u) + O(θ²) is refuted with the pre-registered κ (F2 FAIL).** The numerical gain was 1.17–1.71× θκ
   in every case. The revised coefficient κ₂ adds the burst term B*:

   κ₂ = e^Π [Π m(u) + R(Π, u) − Ein(Π) + B*] / (e^Π − 1).

   It is within 3–9 % of the optimum on **fresh** Π values (2, 4.5, 8) that were not used to derive it, at θ = 0.015 and u ∈ {3, 10}.
   The remaining 5–8 % is the O(θ²) term. It always has the same sign (κ₂ overestimates), so θκ₂ is a mild upper estimate.
3. **G(θ) has an interior maximum (F3 PASS).**
   - For u_max = 3, θ* ≈ 0.27–0.38. G_max is 4.2 %, 16.4 %, 39.4 % and 65.3 % for Π = 1, 3, 6 and 10, and the closed form equals the numerics there.
   - **For u_max = 10 and θ ≳ 0.1 the 3-arc family is no longer optimal.** The numerical optimum was 18.0 % against 11.6 % at
     Π = 3, and 43.5 % against 32.8 % at Π = 6. A fourth arc (coasting at u = 0 just before the terminal burst) appears. The 12-interval
     runs in H1g show it directly: u = …, 0.09, 2.56. In that range the closed form is only a lower bound. The engine flags it
     (`GainEstimate::exact`, `DesignOutput::gainTheoryIsBound`) and clamps the value at ≥ 0, because constant rate is always feasible.
4. **Otaniemi consistency (F4 PASS, within the factor 2).**
   - The E8 small-θ slope is 0.70. The single-exponential κ at Π = 3.06 is 0.90–1.21 without the burst (u 2–4.6) and 1.36–1.81 with it.
   - The missing ingredients are **(a) a diffusive kernel instead of an exponential one, (b) coarse control intervals, and (c) no burst in E8**,
     whose schedule had a fixed number of steps.
   - H1g gives (a) its size: at the same Π, θ, u and interval count, a radial (logarithmic) or linear (√t) kernel reaches only 29–31 % (Π = 3)
     or 59–74 % (Π = 10) of the exponential-kernel gain. 0.70/1.21 = 0.58 falls inside this band.
5. **Kernel shape matters more than expected (H1g, exploratory, not pre-registered).**
   - Diffusive kernels keep the charging arc, which gives a front-loaded schedule. Their terminal uptick is weak: 1.2× and 1.6× the mean
     with u_max = 3, against 2.6× for the exponential kernel.
   - The engine designer on a 41×41 single-phase case (radial flow) behaves the same way. The schedule is
     2.08, 1.20, 1.07, …, 0.67, 0.84 × mean, and the designed gain is 20.4 %, against a first-order reference of 83.9 %.
   - **The closed form is therefore a reference for saturating (first-order) responses, not a predictor for diffusive ones.** The app labels it that way.

## F1 — structure (Π ∈ {1, 3, 6}, θ = 0.05, u_max = 10)

| Π | numeric G | 3-arc | 2-arc (pre-reg) | s* | s(T_i) | t_b/T_i | plateau variation | registered verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.176 % | 2.176 % | 1.237 % | 0.922 | 1.576 | 0.0037 | 1.07 % | FAIL (burst) |
| 3 | 7.946 % | 7.946 % | 5.392 % | 0.932 | 1.392 | 0.0026 | 0.58 % | FAIL (burst) |
| 6 | 18.135 % | 18.135 % | 14.058 % | 0.939 | 1.263 | 0.0018 | 0.24 % | FAIL (burst) |

The plateau criterion (< 2 %) holds. The verdict is FAIL only because the structure has a third arc.

**The "early stop" flag at Π = 3 is an artefact.** It was raised by optimiser ringing just before the burst. On the refined grid, u dips to
0.37 over 3·10⁻⁴ T_i (`h1f.log`, inspected step by step), and the cost equals the 3-arc closed form to 10⁻⁵. No early stop exists in the
optimum. Its first-order cost is positive, as derived in the plan.

## F2 — small-θ law (numeric = discretised optimum; ratio = law / numeric)

| set | Π | u | θ | numeric | pre-reg θκ (ratio) | revised θκ₂ (ratio) |
|---|---|---|---|---|---|---|
| reg | 1 | 3 | 0.01 / 0.02 | 0.418 / 0.808 % | 1.71 / 1.65 | 0.97 / 0.94 |
| reg | 3 | 10 | 0.01 / 0.02 | 1.845 / 3.554 % | 1.43 / 1.38 | 0.96 / 0.93 |
| reg | 6 | 10 | 0.01 / 0.02 | 4.300 / 8.233 % | 1.26 / 1.21 | 0.96 / 0.92 |
| fresh | 2 | 3 / 10 | 0.015 | 1.385 / 1.659 % | — | 0.95 / 0.95 |
| fresh | 4.5 | 3 / 10 | 0.015 | 3.771 / 4.468 % | — | 0.95 / 0.94 |
| fresh | 8 | 3 / 10 | 0.015 | 7.487 / 8.792 % | — | 0.94 / 0.93 |

The full set of 18 points is in `h1f.log`.

## Interpretation

- **Physics.** The optimum uses the fact that the event-count objective treats stress *during* injection and stress *after* shut-in
  asymmetrically. Stress stored at shut-in decays freely, and the volume it represents does not have to be injected while the fault is
  loaded. This is the mathematical form of a known operational worry, raised in the other direction: large post-shut-in events
  (Basel 2006, Pohang 2017). **An event-count optimum raises s(T_i) by 26–58 %.** If the shut-in magnitude hazard matters more than the
  count, the burst is a bad trade. The engine therefore offers a last-interval rate bound. A bound of 1 removes the burst, at the cost
  shown as "without a shut-in burst". On the near-field engine case, a bound of 0.5 × mean costs 2.6 points of gain (20.4 % → 17.8 %).
- **What is new**, as far as the OpenAlex checks show (C20, C21, and today's queries on optimal control / bang-bang / shut-in and
  seismicity, from 2024-06):
  - No paper gives the Pontryagin structure, the terminal burst, or a closed-form gain law.
  - The nearest work is a feedback-control approach to seismicity-rate regulation with actuator saturation (Gutiérrez-Oribio &
    Stefanou, arXiv 2604.26802, Groningen). It is a tracking controller, not a horizon-optimal schedule, and it gives no closed form.
  - Other related work: the LP design of Hill et al. 2024 (JGR, doi 10.1029/2024JB028759) and the cycled-injection lab result (GRL 2026,
    doi 10.1029/2025GL119289).
- **Limits.**
  - A single exponential kernel.
  - The Dieterich linear (t ≪ t_a) regime.
  - Event *count* as the objective.
  - Closed form validated for θ ≤ 0.1, or u_max ≤ 3 and θ ≤ 0.4.

## Second field site (Basel 2006)

- The Basel catalogue is openly available: Herrmann et al. 2019, ETH Research Collection, doi 10.3929/ethz-b-000357354.
- An openly licensed hourly injection log was **not found**. The hydraulic data are published only as figures (Häring et al. 2008).
- **F10 stays CANDIDATE.** Next-cycle option: the Äspö HRL hydraulic stimulation data set (ESSD 2024, 16, 295). It is open and has
  hydraulic time series plus seismicity for conventional, cyclic and pulse injection, which tests the schedule-shape question directly
  at mine scale.

## H2b, H5

- **H2b is deferred.** The cycle went to H1f/H1g and the engine integration. The pre-registered criterion is unchanged.
- **H5 is deferred** (as planned).

## Engine and app integration (C22)

- `seis::scheduleGain(Π, θ, u_max)` uses the 3-arc closed form with a golden-section search over the burst length. It returns the gain,
  the no-burst gain, the burst volume fraction, s*, s(T_i) and an `exact` flag. It reproduces `h1f.log` to < 0.01 % (`-c22`).
- Designer additions:
  - `DesignOutput.{Pi, theta, gainTheory, gainTheoryNoBurst, gainTheoryIsBound, farField, endsWithBurst}`.
  - The option `DesignSetup::lastRateFactor`.
  - A UI field "Last-interval rate / mean" and a "Design potential (first-order reference)" line, with a tooltip and a shut-in-burst warning.
