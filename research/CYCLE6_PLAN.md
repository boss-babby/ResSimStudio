# Research Cycle 6 — plan and pre-registration (written 2026-10-21, BEFORE any model fit on the new data)

## Why the site changed (H8 → H8′)
H8 asked for Basel 2006 as the positive second site. A data search (2026-10-21) found no public hourly
injection log: the Basel flow rate and wellhead pressure exist only as figures (Häring et al. 2008) and on the
EPOS EPISODES platform, which needs a registered account. Digitising a figure would make the main test depend
on hand-made data. **H8′** uses instead the 2020 stimulation of the **OTN-2** well (St1 Deep Heat, Helsinki):
- hydraulics: Kukkonen et al. (2023), Mendeley Data doi 10.17632/ycdcw2rdfj.1, file OTN2.xlsx (CC BY 4.0);
- catalogue: Kwiatek, Martínez-Garzón & Karjalainen (2022), GFZ Data Services doi 10.5880/GFZ.4.2.2022.001
  (6,318 events: 197 before, 5,427 during, 694 after the stimulation).

It is a *different well* (OTN-2, open hole 4.86–5.77 km b.s.l.; 2018 was OTN-3, 5.7–6.1 km), a *different
schedule* (≈ 400 L/min for 84 % of the time, short 1,000 L/min periods, four phases with pauses; 2,875 m³ in
16 days vs 18,433 m³ in 49 days) and a *different year*, but the same crustal block. It is therefore an
out-of-sample test of the model and of the θ-law, **not** an independent tectonic setting. Basel stays open
(H8, needs EPISODES access by the user).

## Pre-registered tests (thresholds fixed now; code = exp/h1e_otaniemi.py with SITE=otn2)
Data handling fixed in advance: local time EEST = UTC+3; 10-min bins; origin 2020-05-05 00:00 UTC; shut-in
= last 10-min bin with rate > 0.01 m³/min; catalogue window to the last catalogue event; M_c by maximum
curvature + 0.2 (as Cycle 3).

- **E1′ kernel (pressure model).** Linear rate → WHP kernel (same 10 time constants, NNLS) fitted on the first
  60 % of the injection time, RMSE on the last 40 % ≤ 10 % of the WHP range → PASS (≤ 20 % WEAK).
  *Prediction:* PASS.
- **E2′ Dieterich fit.** Same likelihood as Cycle 3. PASS if corr(cumulative) ≥ 0.98 and the post-shut-in count
  within ±30 %. *Prediction:* WEAK (Cycle 3 found the tail over-predicted by 50 %).
- **E3′ temporal hold-out.** Fit on the first 60 % of the injection time, forecast shut-in and end counts; PASS if
  both within ±30 % and better than the seismogenic-index baseline at shut-in. *Prediction:* FAIL at the end
  count (stationary Dieterich over-predicts the tail).
- **E9 cross-stimulation transfer (new, the main out-of-sample test).** Keep the 2018 Dieterich parameters
  (P*, t_a, τ_d from h1e_results.json E2), drive them with the 2020 pressure model, refit only the productivity r.
  PASS if corr(cumulative) ≥ 0.95 and the shut-in count within ±30 %. *Prediction:* PASS for the correlation,
  shut-in count uncertain (deeper interval, previous stress relaxation by the M_W 1.2 event before the stimulation).
- **E4′ / F10 θ-law transfer (the P1 claim).** Optimal vs constant rate at the same volume, bounds and window:
  *Prediction:* gain ≤ 3 % (θ ≪ 0.1 as in 2018; the 2020 schedule is close to constant). Operator schedule vs
  constant: *prediction* within ±10 % (few rate excursions) — the sign is not predicted, because the pauses are
  part of the operator schedule and of no constant comparator.
- **E8′ θ-scan (exploratory as in Cycle 3).** Single-exponential kernels at the fitted Π: *prediction* the gain
  peaks for θ between 0.05 and 1 and vanishes at both ends (shape of F10/F12).

Verdict rule for F10: VERIFIED (out-of-sample) if E4′ and E8′ both meet their predictions; otherwise F10 stays
CANDIDATE and the paper says so.

## Also in Cycle 6
- **H7 / P2 (rate-and-state design):** scoped, not run (needs G5); recorded in HYPOTHESES.md.
- **R4 volume-allocation LP:** not run this cycle.
- **Paper Pack P1:** GRL letter (LaTeX, AGU style), Supporting Information, cover letter, data/code archive, JGR long
  version draft; figures rendered by the ResSim figure engine.
