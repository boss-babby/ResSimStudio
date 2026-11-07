# Research Cycle 12 plan (C50, 2026-11-04) — written before any pair run

**Purpose.** Publication target P2 rests on F35. F35's dose-response is *exploratory* and comes from *one geometry*. A reviewer
will ask whether the "≈ 2× transmissibility" threshold is a property of the coupled problem or an accident of the
fault at x = 150 m. H10g is a pre-registered confirmation on a new geometry.

## Why not simply vary permFactor?

Inside the ramp, the transmissibility multiplier is m = P^x with x = CFF/(w σ'n). Hence ln m = k·CFF with
k = ln P/(w σ'n). Changing P at fixed w, or w at fixed P, gives *the same function* as long as the cap x = 1 is not
reached: max CFF ≈ 1.3–2.4 bar ≪ w σ'n ≈ 4.7 bar at w = 0.05. A permFactor sweep would reproduce the ramp sweep
exactly. It is an identity, not a test. An independent test must change the *stress path*, so that the same k reaches
a different m. This is done here by moving the fault.

## H10g (pre-registered)

**Claim under test (from F35).** In the two-way coupled CO₂ / water section, the frequency of midpoint-convexity
violations of the Dieterich expected event count is governed by the transmissibility multiplier m₁ reached at u = 1:
* m₁ ≲ 1.6: violations are rare.
* m₁ ≳ 2.5: violations are frequent.

This is claimed independently of the fault position.

**Geometry G2.** As H10e: CO₂STORE or water, 40×1×10, injector at column 1, layers 8–10, sealing fault MULTFLT 0.02,
FAULTMECH friction 0.6, permFactor 50. The fault is at **column 22 (x = 220 m)** instead of 15 (x = 150 m).

**Calibration (u = 1 runs only; no pair information; fixed rules).**
1. *Pre-load (rule A1 of Cycle 10).* Run u = 1 with the sealing fault and no pre-load, and read max CFF₀. Set the
   pre-load c = 2 − CFF₀ bar, so that the most-stressed face reaches CFF ≈ +2 bar. Check with the switching run that
   3 or more faces reactivate. If not, raise c in 1-bar steps until 3 faces do; this is recorded.
2. *Ramp widths.* For each fluid, choose two ramp widths w by bisection on u = 1 runs. The targets are:
   * LOW: m₁ = 1.40 ± 0.10;
   * HIGH: m₁ = 3.50 ± 0.30.

   If the HIGH target cannot be reached for w ≥ 0.01, use w = 0.01 and record m₁; the prediction still applies if
   m₁ ≥ 2.5.
3. Calibration values are written to `exp/h10g_calib.json` and appended below *before* the pair runs start.

**Cases (4).** G2N-LOW, G2N-HIGH (CO₂); G2W-LOW, G2W-HIGH (water). Same M = 40 pairs, seed 13 (common random numbers,
as H10e / H10f). Same ε rule and the same violation definition as H10.

**Predictions.**
* LOW cases: ≤ 5 % violations (≤ 2 of 40).
* HIGH cases: > 15 % violations (≥ 7 of 40) with median violator δ > 1 %.

**Verdict.**
* SUPPORTED if all 4 predictions hold.
* REFUTED if any LOW case exceeds 15 %, or any HIGH case is ≤ 5 %.
* INCONCLUSIVE otherwise.

The alternative explanation, "violations depend on k (log-slope) alone, not on m₁", is checked descriptively. If k
alone governed the outcome, G2 would need the *same* k as G1 for the same violation rate. If m₁ governs, G2 needs a
*larger* k, because the fault is farther away and sees lower CFF.

**Kill / honesty rules.**
* If the calibration cannot reactivate 3 faces with c ≤ 60 bar, H10g is non-discriminating and is reported as such.
* No re-tuning of w after a pair run has started.
* Anything run after the verdict is labelled exploratory.

## Exploratory (declared now, not part of the verdict)

* E12a: G1 NSR08 / WNSR08 with a second pair seed (29), to show the sampling variability of the violation count at
  M = 40.

## Calibration record (appended before the pair runs)

Recorded 2026-11-04, before any pair run. Source: `exp/h10g_calib.py` → `exp/h10g_calib.json`, about 60 u = 1 runs.

* G1 sanity check: NSR08 m₁ = 2.127 and WNSR08 m₁ = 2.923, the same as F35 (2.13 / 2.92). The binary and deck
  generator reproduce Cycle 11.
* CO₂, G2: CFF₀ = −30.68 bar, so c = 32.68 bar; 7 faces reactivate at w = 0.05, where m₁ = 2.86.
  * LOW: w = 0.27425, m₁ = 1.31 (inside 1.40 ± 0.10).
  * HIGH: w = 0.03761, m₁ = 3.60.
* Water, G2: CFF₀ = −25.06 bar, so c = 27.06 bar; 4 faces reactivate, m₁(w = 0.05) = 3.26.
  * LOW: w = 0.27425, m₁ = 1.35.
  * HIGH: w = 0.04438, m₁ = 3.61.

**Honest note recorded before the pair runs.** Rule A1 normalises the peak CFF to about +2 bar, so G2 reaches almost
the same m₁ as G1 at the same w (2.86 vs 2.92). The descriptive "k alone vs m₁" comparison is therefore
**non-discriminating by construction** and will not be claimed. H10g still tests what P2 needs: whether the
m₁-threshold behaviour transfers to a different fault position and stress path. The LOW/HIGH predictions above are
unchanged.

Pair-run cases: `X:22:0:32.68:0.27425:50`, `X:22:0:32.68:0.03761:50`, `X:22:1:27.06:0.27425:50`,
`X:22:1:27.06:0.04438:50`. M = 40, seed 13.

---

## H6c (pre-registered 2026-11-04, before any stage-level count on the new catalogue was computed)

**Trigger: a data defect found during P3 preparation.** In both GFZ catalogue files (10.5880/GFZ.2.6.2023.004), the
`origin_time` seconds field equals the **day of the month**:
* Kwiatek-match file: 196/196 events, seconds ∈ {03, 04, 09, 11, 12} = the test days.
* Continuous file: 19,756 events, the same pattern.

`event_name` (YYYYMMDDHHMMSS.µs) is monotonic and self-consistent, so it is taken as the true origin time. Cycle 5
(F14–F16) used `origin_time`, so every event was mis-timed by up to ±60 s. This affects stage assignment and the
co-injection / post-shut-in split. The continuous catalogue of Niemz et al. (2020), 19,756 events with AE magnitude
MAE for 4,302, is public in the same archive (`seismic_data/AE_continuous/CATALOG`). Cycle 5 report item R3, which
said "needs a data request", was wrong.

**Analysis.** The Cycle 5 pipeline (`h6_aspo.py`: stage definition, kernels, Poisson fits, training tests HF1+HF2)
is reused **unchanged**, except for the event times and the catalogue:
* (i) the 196-event catalogue with corrected times (Niemz event-name times of the matched events; events whose name
  cannot be parsed are dropped and counted);
* (ii) the continuous catalogue, all events;
* (iii) the continuous catalogue with MAE ≥ 2.5, the completeness magnitude of Niemz et al. 2020. Exploratory
  sensitivity checks: MAE ≥ 2.15 and ≥ 3.0.

**Predictions (applied to (ii) and (iii); (i) decides whether F14–F16 must be corrected).**
* P-a: the Dieterich pressure kernel K1 is worse than the better of energy (K2) and volume (K3) by ΔAIC ≥ 10.
* P-b: K2 and K3 both beat duration (K5) and constant (K0) by ΔAIC ≥ 10.
* P-c (on (iii)): the HF3 rate factor is below the AG-conventional mean (ratio < 1) and is **not** significant at the
  one-sided 5 % level against the inter-test scatter of conventional tests (the F16 method).
* P-d (on (iii)): the observed post-shut-in fraction (HF1+HF2) is below the Dieterich prediction.

Each prediction is SUPPORTED or REFUTED separately. Results contradicting F14–F16 replace them in FINDINGS_LOG
(marked CORRECTED); the old values are kept for the record.

**Robustness (declared now):**
* leave-one-test-out for the H6b factor ratio;
* a stage-bootstrap (1000 resamples of stages within tests) of ΔAIC(K1 − best of K2/K3);
* the energy-vs-volume difference reported with its bootstrap interval. A difference < 2 is not claimed as a
  preference.

---

## RESULTS (appended after all runs, C50, 2026-11-04; the text above was not edited)

### H10g: transfer of the convexity-loss threshold to fault position G2 (x = 220 m): **SUPPORTED**
| Case | Calibrated | m₁ | Violations / 40 | Median violator δ | Prediction | Verdict |
|---|---|---|---|---|---|---|
| CO₂ LOW | c 32.68 bar, w 0.27425 | 1.31 | 0 | — | ≤ 2 | met |
| water LOW | c 27.06 bar, w 0.27425 | 1.35 | 0 | — | ≤ 2 | met |
| CO₂ HIGH | w 0.03761 | 3.60 | 21 | +4.56 % | > 6 and > 1 % | met |
| water HIGH | w 0.04438 | 3.61 | 16 | +5.12 % | > 6 and > 1 % | met |
Files: exp/h10g_calib.py, h10g_calib.json, h10g.log, h10g_results.json. Finding F37.

### Exploratory (not pre-registered)
* **E12 (seed 29, intermediate m₁):** CO₂ m₁ 2.02 → 6/40, 2.60 → 13/40; water 2.01 → 1/40, 2.61 → 11/40; NSR08 5, WNSR08 14.
  Pooled logit over 18 configurations: p = −4.965 + 3.693 ln m₁; 5 % at m₁ 1.73 [1.53, 1.93]; 15 % at 2.40 [2.25, 2.56]; geometry
  p 0.058, fluid p 0.082, φ 1.90.
* **E12c (violators only, dt 2.5 d):** regression-to-the-mean biased; superseded by E12d.
* **E12d (all pairs):** CO₂ 2.5 d 20/40 (median violator 3.6 %); water 2.5 d 13/40 (1.6 %); water 1.25 d 13 of 16 original violators
  persist (1.45 %). The loss is physical; the lag inflates the water amplitude about 3.5×. Finding F38.

### H6c: Äspö kernels on corrected times and the continuous catalogue: **P-a, P-b, P-c, P-d SUPPORTED** (primary M_AE ≥ 2.5)
* Data defect: origin_time seconds = day of month in both GFZ files (F39); times taken from event names.
* M ≥ 2.5 (2,322 events, 31 stages): AIC energy 246.4, volume 299.4, Dieterich 791.4, duration 982.4. Bootstrap ΔAIC(K1 − best)
  [70, 981], all ≥ 10. Post-shut-in 0.143 vs 0.500 (z −27). HF3 0.41, one-sided p 0.04–0.31 (8 analysis variants).
* Rate at fixed volume: β 0.45–0.65 (K6), 0.4–0.9 (K8): **corrects F14**. HF5 is not covered by the continuous catalogue.
Files: exp/h6c_aspo.py, h6c.log, h6c_hf4.log, h6c_results.json. Finding F40.

### Consequences
* Manuscripts: P2 v1.0 and P3 v1.0 written; P1 GRL corrected (Äspö moved to P3; two-way overclaim removed; Kim et al. 2025 cited).
  See papers/PORTFOLIO.md.
* Next cycle: N1 (cyclic storage), see NEW_FIELDS.md; H12a/H12b to be pre-registered in CYCLE13_PLAN.md.
