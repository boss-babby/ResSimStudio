# Research Cycle 5 report: a second site (Äspö HRL), the kernel test, and an honest verdict

Pre-registration: `CYCLE5_PLAN.md`, written before any fit. Script: `exp/h6_aspo.py` (numpy only, about 5 s).
Outputs: `exp/h6.log`, `exp/h6_results.json`, `exp/h6_fig.png`. Data: `data/aspo/` (1 s means of GFZ doi
10.5880/GFZ.2.6.2023.004, CC-BY 4.0, plus the 196-event triggered catalogue of Kwiatek et al. 2018).
Engine integration: `seis::fitStageKernels` (seismicity.h) and the *Calibrate from observed stage events* panel in the
seismicity designer. Regression test: `tests/c23_tests.inc` testKernelFitAspo, run with `-c23`.

## Data and stage table

- **Stages.** I found 33 injection stages in six tests. The two longest tests (HF1, HF2) supply 12 stages with 140 events between them.
- **Event assignment.** 184 of the 196 catalogue events fall into a stage window: HF1 48, HF2 92, HF3 16, HF4 1, HF5 0, HF6 27.
- **Deviation from the pre-registration (documented).** Four intervals had flow above 0.3 L/min while the interval pressure stayed near 0 MPa. These are packer inflation or line filling, not injection into the test interval. A stage now also needs a median interval pressure of at least 1 MPa. This removed the 4 events of HF2.0.0. The rule was fixed before any kernel was compared, but after I looked at the stage table.

## H6a: kernel discrimination (pre-registered)

Training set: HF1 + HF2, Poisson likelihood with one rate factor per test.

| Kernel | log L | AIC |
|---|---|---|
| K2 hydraulic energy ∫p q dt | −38.2 | **80.4** |
| K3 volume | −41.2 | 86.3 |
| K1 Dieterich (S = p, best Aσ = 1.50 MPa) | −65.6 | 137.2 |
| K5 duration (reference) | −77.0 | 158.0 |
| K0 constant per stage (reference) | −89.7 | 183.4 |
| K4 Kaiser (energy above the previous maximum) | −130.9 | 265.7 |

**Pre-registered prediction "K1 is not preferred": CONFIRMED.** Dieterich driven by the interval pressure is 57 AIC units
worse than the energy kernel, and the Kaiser kernel is worse still. Energy is preferred over volume by 5.9 AIC units.
The stage pressures vary only between 8 and 11 MPa, so energy and volume are nearly proportional.

The out-of-sample check on HF6 was uninformative. HF6 has only 4 stages and 27 events, and its deviances lie between 2.3 and 4.8 for all
kernels except K4 (90.9).

**Exploratory, not pre-registered: does the rate matter at fixed volume?**

- On its own, V·q^β gives β = 0.70 (95 % interval 0.35–1.05).
- In this step-rate protocol, however, the rate rises with stage number. Rate and earlier injected volume correlate at ρ = 0.97.
- V·(1+V_cum)^γ fits better (AIC 71.5, γ = 0.45 [0.25, 0.70]).
- In the joint model the rate exponent collapses to β ≈ −0.1 to −0.25, with a 95 % profile interval of **[−1.25, 0.75]**.

So the apparent rate effect is the growing fracture: later stages extend a larger fracture. **I could not resolve a rate effect at fixed volume.**

## H6b: is the "cyclic reduction" distinguishable from scatter between intervals? (pre-registered)

- **Rate factors.** With the energy kernel, the per-test rate factors are HF1 0.21, HF2 0.33, HF3 0.064, HF4 0.026, HF5 0 and HF6 0.11.
- **Scatter.** The conventional tests spread with log-σ = 1.11.
- **HF3 against same-rock tests.** HF3, the progressive scheme, is 0.25× the mean of the same-rock conventional tests. Its t-score is −1.03, and the 5 % predictive bound is 0.04×.

**Prediction "not significant": CONFIRMED.** The published 6× reduction (102 → 16 events) rests on a single progressive test.
One conventional test in another rock type (HF4) had even fewer events per unit energy than HF3.

This does not show that cyclic injection fails. It shows that the Äspö data cannot demonstrate a benefit: n = 1, and the scatter between intervals is about a factor of 3.

## H1h: post-shut-in fraction (pre-registered)

- **Observed.** 35 of 140 training events (0.25) occurred after shut-in.
- **Predicted.** Dieterich with the fitted Aσ predicts 0.48: the interval pressure holds for minutes after shut-in, so the model keeps producing events.

**Prediction "Dieterich under-predicts post-shut-in events": REFUTED.** It over-predicts them. The AE follow the flow (the growing fracture),
not the pressure. For hydraulic-fracture-dominated AE, the terminal-burst hazard argument of F12 has no empirical support at this scale.

## What this means: honest verdict

1. **Does Cycle 5 produce a good output?** Yes, but it is a *scope* result, not a breakthrough.
   - At mine scale, AE from hydraulic fracturing are **volume-driven**: N ∝ V·(1+V_cum)^0.45.
   - They are not Dieterich-pressure-driven, and the rate has no resolvable effect.
   - Under a volume kernel, the count at fixed volume does not depend on the schedule. The P1 scheduling gain is therefore zero there.
   - This is a real and useful constraint, and it is new for the Äspö data. I found no stage-level kernel comparison or confound analysis of Äspö: OpenAlex, 2019–2026, 2 queries, recorded in `lit/`.
2. **Does it go into P1?** Yes, as a **second-site / domain-of-validity section** with one figure (h6_fig).
   - P1 states its domain: shear reactivation of pre-existing faults by pore-pressure diffusion (Otaniemi/Basel type), where Dieterich applies.
   - It excludes tensile-fracture-growth AE, where only volume and site matter.
   - Reviewers will ask exactly this, so it makes P1 stronger. It does **not** replace a second *positive* site; Basel F10 remains the candidate.
3. **Novelty.**
   - **P1 core (F1 convexity, F12 bang–singular–bang + κ₂ law): still unprecedented** in the literature I searched, and still P1-worthy (GRL → JGR SE).
   - **F14 (volume-driven AE, rate confounded by step-rate protocols):** medium novelty. The empirical energy correlation is known (Zang 2017, Niemz 2020). The formal kernel test and the confound are not.
   - **F15 (cyclic reduction not significant):** a reanalysis. It suits a short comment or data note, or a paragraph in P1. It is low impact on its own.
4. **Is the line exhausted?** Not entirely. The schedule-design line pays off only where Dieterich-type shear reactivation dominates. The next scientifically valuable step is G5 (below), not more Dieterich-proxy optimisation.

## Next directions (research)

- **R1, G5 with G3: optimal injection against simulated rate-and-state slip.**
  - Open question: does convexity survive elastic interaction and nucleation physics?
  - The G3 module (rsf.h) now produces the event catalogue needed to test this.
  - It is a genuinely new paper candidate, P2: "optimal injection under rate-and-state fault interaction".
- **R2, Basel 2006 as a positive second site.** Digitise the hourly log from Häring 2008 figures. The catalogue is open at ETH.
- **R3, a better-powered Äspö test.** Use the continuous catalogue of Niemz et al. 2020, about 19,600 events, b-values per stage. This needs the waveform picks or a data request.
- **R4, the seismogenic index under volume control.** Since N ∝ V at mine scale, the design variable is volume allocation between stages and wells (F6: worth about 8× the schedule shape). Allocation under a volume kernel with per-interval site factors is an LP. That is known in principle (Hill 2024), and the app can expose it.

## Product directions (if the research line stalls)

- **P-a: seismicity kernel diagnostics.** Done this round: `fitStageKernels` and the UI panel. An operator pastes past stage data and gets a yes/no on whether rate scheduling is worth trying. This is useful even where the answer is "no".
- **P-b: adaptive traffic-light forecaster.** Fit the seismogenic index and RSF to streaming rate/pressure/catalogue data during a stimulation. It would forecast the next-hour rate and Mmax (ETAS + Σ index), with a stop/continue recommendation. Geothermal operators and regulators pay for this today, using Python scripts.
- **P-c: EGS virtual stimulation studio.**
  - Chain: HFGROW (growth), EFEM contact (G1-c), RSF fault cycles (G3) and the event catalogue, then report export.
  - Result: one workflow from deck to "expected catalogue + permeability gain".
  - No commercial simulator combines all of these in a desktop GUI.
- **P-d: CCS screening and permitting pack.** Fast pressure front, fault-stability and RSF hazard over a licence area, with automatic publication-quality figures and report tables. Permit applications in the EU and US need exactly these figures.
- **P-e: batch HM/UQ service.** Headless CLI plus MPI batch for ensemble history matching, behind a small web front end.
