# Research Cycle 8 report (C29, 2026-10-23)

Three hypotheses were tested with criteria fixed before the runs (HYPOTHESES.md rows H7b / H2c / H8b,
CYCLE8_PLAN.md, header of `exp/h7b_rsf_ensemble.cpp`). One engine defect was found while doing it, and it changed an
earlier finding (F21). Everything is reproducible from `exp/`.

| hypothesis | question | verdict | finding |
|---|---|---|---|
| H7b | Is the *ensemble-mean* rate-and-state event count convex in the rate history? | **REFUTED** | F22 |
| H2c | Does sub-grid leak-off stay grid-converged on a dipping corner-point grid and for gas injection? | **SUPPORTED** (CP); GAS non-discriminating | F24 (+ F21 amended) |
| H8b | Utah FORGE 16A(78)-32, April 2022: volume kernel or Dieterich pressure kernel? | **INCONCLUSIVE** | F23 |

## H7b: does averaging over fault heterogeneity restore convexity?

**Why.** F20 (Cycle 7) showed that the event count of one rate-and-state (RSF) fault is not convex in the rate history.
The violations were never larger than one event, so integer quantisation could have caused them. H7b averages over
E = 16 realisations of a ±20 % heterogeneous normal stress, using common random numbers across schedules.

**Result** (`exp/h7b_rsf_ensemble.cpp`, 60 pre-registered pairs, 1192 s):
* Strict midpoint violations of the expected count: **14 / 60 (23 %)**, which is above the pre-registered REFUTED band of
  15 %. 12 of the 14 exceed one event in the ensemble total (0.06–0.75 events per realisation).
* Cumulative moment: 19 / 60 violations.
* Dieterich expected count with the same pressure: 0 / 60. Over all of H7 + H7b this is 0 / 290.
* The heterogeneity really decorrelates the catalogues: only 16 / 60 pairs have counts that are exact multiples of E.

**Meaning for P2.** A rate-schedule optimiser with a global-optimum guarantee must use the Dieterich expected rate
(convex, F1). Neither a single nor an ensemble-averaged RSF catalogue is a convex objective. P2 therefore needs either a
stochastic (non-convex) optimiser or a Dieterich surrogate calibrated to the RSF fault.

## H2c: sub-grid leak-off beyond the H2b decks

**Design** (CYCLE8_PLAN.md). The decks are made from `examples/HFGROW_INJECTION.DATA` by `exp/h2c_mkdeck.py`, at
dx = 10 / 5 / 2.5 m, with the sub-grid tip only (SG) and with sub-grid leak-off (SGL):
* **CP:** COORD / ZCORN grid with a 2° dip, water injection.
* **GAS:** dry gas at 3000 sm³/d into the oil zone.

The H2b water decks were rerun with the C29 engine as the new baseline (**WAT**). The measure is t60 / t90, the times at
which a wing first reaches 60 / 90 m.

| variant | scheme | t60 at 10 / 5 / 2.5 m (d) | spread | t90 at 10 / 5 / 2.5 m (d) | cuts at 10 / 5 / 2.5 m |
|---|---|---|---|---|---|
| CP | SG | 25.84 / 24.53 / 24.00 | 1.84 | 57.16 / 55.40 / 54.66 | 36 / 41 / 105 |
| CP | **SGL** | 23.00 / 22.80 / 23.65 | **0.85** | 53.41 / 52.63 / 53.62 | 16 / 52 / 99 |
| WAT (C29 engine) | SG | 28.67 / 27.42 / 26.84 | 1.83 | 63.41 / 61.47 / 60.85 | 38 / 70 / 171 |
| WAT (C29 engine) | **SGL** | 25.31 / 25.55 / 26.00 | **0.69** | 58.70 / 58.69 / 59.55 | 17 / 65 / 141 |
| GAS | SG | 0.49 / 0.46 / 0.44 | 0.05 | 1.06 / 1.02 / 0.98 | 7 / 6 / 2 |
| GAS | **SGL** | 0.44 / 0.45 / 0.42 | **0.03** | 0.95 / 0.98 / 0.98 | 9 / 4 / 1 |

**Verdicts (pre-registered rule: SGL spread < 1 d and smaller than SG → SUPPORTED).**
* **CP: SUPPORTED.** 0.85 d against 1.84 d. SGL is not monotone in dx (the 5 m grid is fastest), so the remaining error
  is about 3 % noise rather than a converging sequence.
* **GAS: formally SUPPORTED, but not informative.** Gas reaches 60 m in under half a day, and the 1-day threshold was
  chosen for growth that takes about 25 days. Relative spreads are 7 % (SGL) against 11 % (SG). A fair gas test needs a
  lower rate or a higher K_IC, so the growth takes weeks. That is left for a later cycle.

**F21 amended.** In C28 the SG baseline was 33.58 / 28.77 / 26.97 d (spread 6.6 d). With the C29 engine the same decks
give 28.67 / 27.42 / 26.84 d (spread 1.83 d). SGL hardly moved (25.49 → 25.31 d at 10 m). So about two thirds of the
"grid dependence" that F21 removed was caused by the old well start (next section), not by the leak-off quantisation. The
leak-off correction still works: with SGL the spread is 0.69 d against 1.83 d, and the 10 m grid needs half the time-step
cuts. However, the gain is about 2.7×, not the 8× that Cycle 7 reported. The user guide and ENGINES_AND_MATH.md now quote
the new numbers.

## Engine defect found: well start in a fracture cell (well repair v1 → v2 → edge fix)

Gas injected into an EDFM fracture cell flip-flopped between "well dead" and "well alive" at t = 0 and cut the first step
many times. The fix went through three versions:
1. **v1** reset the BHP to a repaired value. It could stall at a fixed point.
2. **v2** bisects between the failed BHP and the static well pressure b0.
3. **Edge fix.** v2 took the static "edge" as the maximum over all connections. For water injectors on the CP grid, that
   made the well start too early and broke three phase-field tests. The edge is now the **minimum over connections for
   injectors and the maximum for producers**: the pressure at which the first connection can flow.

With the edge fix `-c29` passes and the full suite passes (ENGINES_AND_MATH.md §6.22). Time-step paths changed, so
reference values shifted by about 2 % (WAT 10 m SGL t60 25.49 → 25.31 d). The large SG shift (33.58 → 28.67 d) is the
F21 amendment above.

## H8b: Utah FORGE 16A(78)-32 as the third P1 site

**Data.** Pump records come from GDR 1379 (1 s treating pressure and slurry rate, stages 1–3). The catalogue is GES v2
(GDR 1429, status-2 events, M ≥ −1.5). `data/forge/` holds them; they are not shipped in the release.

**Clock problem (found before any fit).** If the pump clock is read as local MDT, every stage has 197–1623 events up to
50 min *before* pumping starts. A clock one hour off (UTC−5) removes them, so the pre-registered primary shift was
−60 min.

**Fit** (5-min bins, 3 stages, 204 bins, Poisson with one factor per stage, `exp/h8b_forge.py`, `h8b_results.json`):

| kernel | AIC at −60 min |
|---|---|
| constant | 4922 |
| volume (1 % floor) | 39357 |
| energy (1 % floor) | 40523 |
| volume with memory (τ ≈ 2 h) | **4091** |
| Dieterich, surface pressure (Aσ ≈ 50 MPa) | 4556 |

* At −60 min the volume type wins by 465.
* At −50 min Dieterich wins by 1694.
* At −70 min the volume type wins by 1899.

The pre-registered robustness condition fails, so the verdict is **INCONCLUSIVE**.

**Robust parts.**
* Dieterich with a realistic Aσ (0.1–1 MPa) is rejected at every shift. This is because the *surface* treating pressure
  jumps about 40 MPa at pump start; that jump is friction, not reservoir stress.
* The event rate stays near 400 per 10 min while the rate ramps from 5 to 35 bpm. This suggests the catalogue saturates
  at high activity.
* Post-shut-in fractions are 0.147 / 0.151 / 0.475.

**Consequence.** FORGE becomes a P1 site only with downhole (or friction-corrected) pressure and a confirmed pump clock.
P1 keeps two field sites (Otaniemi, Äspö). OpenAlex screening (Oct 2026) found FORGE 2022 modelling and rate-forecasting
papers, but no kernel comparison with a clock-alignment check.

## Literature watch
OpenAlex queries (FORGE induced seismicity and stimulation 2022–2026; RSF ensemble convexity) found nothing that
pre-empts F22 / F23 / F24. Titles were screened only; Scholar is still to be checked for F21 / F24 (sub-cell storage in
EDFM growth).

## Next cycle (proposal)
* Gas H2c at a growth time of weeks (rate 300 sm³/d or K_IC 3 MPa√m).
* For P2, a Dieterich surrogate calibrated to the RSF fault, plus a stochastic optimiser comparison.
* FORGE: look for a downhole gauge or a friction-corrected bottom-hole pressure in GDR 1558 / 1674.
