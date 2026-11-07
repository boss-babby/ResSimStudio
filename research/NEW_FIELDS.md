# New research fields: critical screening memo (C50, 2026-11-04)

**Rule:** existing work (P1–P3) is released first. A new field enters only if it (a) reuses what we already have (theorem, engine,
FAULTMECH, data pipelines), (b) has a real prior-art gap after a search, and (c) has a test that can fail.

## Screening table

| ID | Field | Why it may be new (prior art checked) | Reuses | Can fail by | Verdict |
|---|---|---|---|---|---|
| **N1** | **Seismicity-aware cyclic storage scheduling (UGS / UHS / CAES)** | Current storage work uses reactivation thresholds: pressure bandwidths (Frontiers Earth Sci. 2026, 10.3389/feart.2026.1933751), safe envelopes from parametric runs (PNAS 2026, 10.1073/pnas.2617910123), deep-learning surrogate QRA (OSTI 2588094), forward HM studies (PMC11363680; EGU26-18385). None optimises a seismicity-rate objective, and none states the structure of the design problem | P1 theorem; engine gas storage + FAULTMECH; SCP loop (F27) | Convexity fails in the engine because of nonlinear gas PVT, hysteretic relperm or the explicit lag; or Dieterich and Kaiser memory give opposite designs | **TOP**: Cycle 13 |
| N2 | Robust / Bayesian design | sup over a parameter set and expectation over a posterior both preserve convexity → robust design stays a convex programme. Not stated for seismicity design (P1 search) | P1, F19 (identifiability) | Only theory; the test is whether robust designs beat nominal ones out of sample at Otaniemi 2020 | fold into P5 (with N1) |
| N3 | Hazard-constrained design | P(M ≥ M_max) = 1 − exp(−N·10^(−b(M_max − M_c))) is increasing in N, so {P ≤ α} = {N ≤ N_α} is a convex set: maximum injection under a hazard cap is convex | P1 | — (a corollary; testable only through N1/N2 applications) | corollary in P5 |
| N4 | Geothermal flexible dispatch (load-following doublets) with a seismicity budget | Reinjection rate schedules for flexible geothermal output; seismicity treated by traffic lights (Lindner et al. 2026 GJI ETLS) | P1, thermal engine | Thermal stress makes the map nonlinear (as F26: check sign) | later |
| N5 | CO₂ storage with brine extraction (pressure management) | Signed multi-well rates; F25/F26 suggest convexity survives | P2 deck family | Extraction creates a non-monotone ΔS (still linear → convex; the test is the engine nonlinearity) | later; small extension of P2's follow-up |
| N6 | ETAS triggering on a Dieterich background | **Correction of an earlier note:** with a background λ₀ that is convex in q and a stationary branching ratio n < 1, the expected total count is λ₀-integral/(1 − n), still convex. Convexity breaks only if productivity or n depends on pressure or rate | P1 | A pressure-dependent productivity (literature: Otaniemi aftershock productivity) | theory note; check before claiming |

## N1 in detail: why it is the best next field

**Theory sketch (to be verified, not yet a finding).** In storage operation q(t) is signed (injection > 0, withdrawal < 0). Under
linear poroelastic coupling ΔS(x, t) = Σ_w ∫ g_w(x, t − s) q_w(s) ds is still linear in q. The Dieterich count is a convex function of
ΔS (P1, Theorem 1), so the count is convex in the signed schedule. Typical storage constraints are linear (inventory balance, working
gas volume, deliverability bounds, a seasonal demand profile), so the feasible set is a polytope. Minimising expected events subject to
the demand profile is therefore a convex programme.

**The scientific question that can fail.** Field storage seismicity (e.g. Norg, Grijpskerk; Frontiers 2026 discussion) suggests a
*stress-memory* (Kaiser-type) behaviour: events occur when stress exceeds its previous maximum. A Kaiser-type objective
f(max_t ΔS(x, t)) with f convex and increasing is also convex (max of linear functions is convex). So **both memory models give convex
design problems, but they may give different optima.** Dieterich penalises time spent at high stress; Kaiser penalises only new
maxima. Whether the optimal cycle shapes differ, and by how much, is the testable hypothesis H12b.

**Pre-registration draft (Cycle 13, to be finalised in CYCLE13_PLAN.md before any run).**
* **H12a (engine convexity).** Gas-storage deck (dry gas, one fault, FAULTMECH one-way, then two-way with m₁ < 1.5), seasonal
  signed schedules on 4 control intervals. Midpoint test with 40 pairs, verdict rules as H10. Prediction: SUPPORTED for one-way.
* **H12b (memory-model discrimination).** Optimise the same demand-constrained schedule under Dieterich and under the Kaiser
  objective. Prediction: the optimal schedules differ by > 10 % in the L1 norm, and each model's optimum costs > 5 % under the
  other model. If neither holds, the memory question does not matter for design, which is a useful negative.
* **Data feasibility (no hypothesis).** Check public availability of Norg / Grijpskerk pressure histories (NLOG) and KNMI catalogues.
  Without data, N1 stays a synthetic method paper.

**Expected output:** P5, "Convex scheduling of cyclic subsurface storage under seismicity-rate and hazard constraints" (N1 + N2 + N3).
Target: GRL or IJGGC/Applied Energy, depending on whether field data can be used.

## Fields considered and rejected
* **Machine-learning surrogates for seismicity design:** crowded field (OSTI 2588094 and many others), and convexity makes them
  unnecessary in the convex regime. Rejected except as the GP-EI tool already in P2 for non-convex objectives.
* **Full dynamic-rupture optimisation:** outside the engine's scope (quasi-static FAULTMECH); no realistic test.
* **Basel / Pohang replication of the θ-law:** scientifically wanted (P1 weak point), but the injection logs are not public. Revisit if
  data appear.

## Cycle 13 outcome (C50, 2026-11-04)
* **H12a SUPPORTED** (F41): convexity holds in the engine for signed storage cycles, water and depleted gas.
* **H12b SUPPORTED** in both regimes (F42). Exploratory: a peak-first (Kaiser) design with Dieterich tie-break is robust (< 0.35 %
  regret under both models); a Dieterich-only design can lose 40 % at low Π.
* **N1 status: promising, not yet a paper.** Needed for P5:
  * (i) the same test on the gas case (response matrix for SG via SCP, since it is nonlinear);
  * (ii) a second geometry and multi-cycle horizons; Kaiser memory across cycles is where the models should diverge most;
  * (iii) an equilibrated initial state, to remove the 0.68-bar drift artefact;
  * (iv) a field data feasibility check (Norg / Grijpskerk, NLOG + KNMI).
  Estimated 1–2 cycles.
