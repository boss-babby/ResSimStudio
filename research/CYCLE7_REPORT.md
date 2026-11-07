# Research Cycle 7 report (C28, 2026-10-22)

Two hypotheses were tested with criteria fixed before the runs (pre-registered in HYPOTHESES.md):
* **H7**: convexity of the event count under rate-and-state slip.
* **H2b**: sub-grid leak-off for in-run fracture growth.

The literature watch was also run. Everything is reproducible from `exp/`.

## H7: is the event count still convex when the seismicity comes from a rate-and-state fault?

**Why it matters.** P1's main result is F1: the expected Dieterich event count is convex in the injection-rate history. That
is why the designer's optimum is global. The open question for a second paper (P2) is whether this survives when the events
come from a physical fault, with elastic interaction, nucleation and stress drops, instead of a rate model.

**Design.** `exp/h7_rsf_convexity.cpp`:
* *Fault.* The engine's 2-D quasi-dynamic rate-and-state fault (`rsf.h`, G3): 1.5 km, 64 elements, velocity-weakening
  centre 25–75 %. There is no tectonic loading, so every event is injection-induced.
* *Forcing.* A well 20 m from the fault centre, with 2-D Theis superposition (P* 1.5 MPa, c 0.05 m²/s). Each schedule has
  4 intervals of 100 d at the same total volume.
* *Pairs.* Random schedule pairs (A, B) and their midpoint C = (A+B)/2.
* *Violation.* N(C) > (N_A + N_B)/2.
* *Verdict rule.* Supported if ≤ 10 % violations, refuted if > 25 %, otherwise inconclusive.

**Results.**

| set | pairs | count violations | of which > 1 event | moment violations | Dieterich count violations |
|---|---|---|---|---|---|
| pre-registered (seed 7) | 30 | 7 (23 %) | 0 | 12 (40 %) | 0 |
| robustness (seed 11) | 200 | 30 (15 %) | 0 | 82 (41 %) | 0 |

Each run gives 2–4 events. The midpoint moment reaches up to 1.27× the average of A and B.

**Verdict.** INCONCLUSIVE for counts (the pre-registered band is 10–25 %) and clearly non-convex for moment.
* No violation exceeds one event (0 of 230). The counts are therefore convex up to the integer quantisation of a single
  fault's catalogue, and the violations come from rounding, not from a systematic concave response.
* The cumulative moment is not convex: one large rupture more or less dominates it.
* The Dieterich count driven by the same pressure history never violates convexity (F1 confirmed again).

**Finding F20.** The convex design of P1 applies to the *expected* count, which describes a population of faults, or
probabilistically a single fault over many realisations. It does not apply to the deterministic catalogue or the moment of
one fault.

**Consequence for P2.**
* State the objective as an expectation, for example over random initial stress and state, or over a population of faults.
  A test with an ensemble of initial stress fields is the natural next step (H7b).
* A moment or Mw-based objective needs a different, non-convex treatment.

## H2b: sub-grid leak-off for EDFM fracture growth

**Why it matters.** F11 (Cycle 3) found that the grid error of in-run fracture growth comes from host-cell leak-off
quantisation. When the tip enters a fresh host cell, the whole cell is still at the initial pressure, so leak-off per area is
about 30× too large. The column scheme and the area-weighted sub-grid scheme bracket the converged answer from opposite sides.

**Method.** The engine now has `HFGROW … SUBGRID LEAK` (simulator.cpp `setHfTips`; research switch RESSIM_SUBLEAK=1).
* The exposed part φ of the tip host cell is treated as its own storage. Everything the cell absorbed since the tip entered
  sits in that part: p_sub = p_rest + (p_cell − p_rest)/φ.
* The leak-off of the tip piece becomes φT(p_f − p_sub) instead of φT(p_f − p_cell), applied as a lagged multiplier
  m ∈ [0, 1].
* m → 1 continuously as φ → 1, so the correction removes itself and needs no tuning parameter.

**Pre-registered pass.** The spread of t60 over dx 10 / 5 / 2.5 m is below 1 d.

**Results.** HFGROW_INJECTION deck family, 200 m trace, K_IC 1 MPa√m, σ_min 340 bar (`exp/h2b_results.txt`).

| scheme | t60 at 10 / 5 / 2.5 m (d) | spread | t90 (d) | cuts at 10 m | wall at 10 m |
|---|---|---|---|---|---|
| column opening (C21 data) | 22.0 / 24.0 / 25.0 | 3.0 | 56.0 / 56.0 / 58.5 | – | – |
| sub-grid tip, area-weighted leak-off | 33.58 / 28.77 / 26.97 | 6.6 | 71.44 / 69.60 / 62.32 | 272 | 42 s |
| **sub-grid tip + sub-grid leak-off** | **25.49 / 26.34 / 26.00** | **0.85** | 58.96 / 59.78 / 60.21 | 147 | 35 s |

**Verdict: SUPPORTED** (finding F21; F11 moved to VERIFIED).
* At 10 m cells the coarse grid gives the 2.5 m answer within 2 %. It needs fewer time-step cuts and is faster.
* The corrected value lies inside the F11 bracket (25.0 … 26.97 d), as predicted.
* The option is in the app: Hydraulic fractures → growth table → *sub-grid leak-off (0/1)*, and the deck `SUBGRID LEAK`.

**Limits.**
* One deck family: single-phase water, 2-D EDFM, PKN-type K_I.
* The correction ignores lateral matrix flow inside the tip cell and is lagged by one step.
* Next: multiphase, and corner-point grids.

**Novelty.** A quick OpenAlex check found no sub-cell storage correction for EDFM growth. The usual practice is cell-wise
opening. Google Scholar was not reachable, so a proper check is still open. Its value is as a methods note (for example
Computational Geosciences or SPE Journal technical note) and as an engineering improvement in the app.

## Literature watch (OpenAlex, since 2026-06-24)
* 32 new items. Two queries were added (rsf-schedule, hf-leakoff-grid); the rsf-schedule query picks up computing "fault
  injection" noise and should be narrowed next cycle.
* Relevant: Science Advances 2026 on fluid-driven aseismic slip fronts (doi 10.1126/sciadv.aeb0234; context for H7);
  SPE 233898 on Utah FORGE THM + ML forecasting (a possible third-site data source for P1); PFLOTRAN sequential
  thermo-poroelasticity verification (arXiv 2607.01156).
* None of them states the convexity, the schedule structure or the sub-grid leak-off results.
* Digest: `lit/digest_2026-10-22.md`.

## Status of the publication plan
* **P1** (convex induced-seismicity scheduling, GRL + JGR pack v0.5): unchanged. F20 adds one honest scope sentence: the
  result is about expected counts, not single-fault catalogues. A third site (Utah FORGE) is still the main gap before
  submission.
* **P2** (rate-and-state scheduling): F20 redefines it as an ensemble-expectation problem. Next is H7b (an ensemble of
  initial stresses).
* **Methods note** (sub-grid leak-off, F21): ready for a short paper after a multiphase check.

## Next cycle (8) candidates
* H7b: ensemble-averaged RSF counts (random initial τ/σ′ and state), testing whether the expectation is convex
  (pre-register ≤ 5 % violations).
* H2c: sub-grid leak-off for multiphase and gas injection, and on corner-point grids.
* H8: Basel as the positive second site (digitised injection log), or Utah FORGE.
