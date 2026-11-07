# Research Cycle 5 — plan and pre-registration (written 2026-10-20, BEFORE any model fit)

Site: Äspö HRL 2015 mine-scale hydraulic-fracturing experiment (Zang et al., ESSD 16:295, 2024;
data doi 10.5880/GFZ.2.6.2023.004, CC-BY 4.0). Six tests HF1–HF6 at 410 m depth, 196-event
triggered AE catalogue (Kwiatek et al. 2018, Mw −3.5 … −4.2), 1 s hydraulics (flow, interval pressure).
Local copies: `research/data/aspo/` (1 s means; catalogue csv; README).

Schemes (Zang 2017 GJI Tab. 1; Kwiatek 2018 JGR Tab. 1): HF1, HF2 conventional (Ävrö granodiorite, AG);
HF3 progressive cyclic initial frac then constant refracs (AG); HF4 conventional (fine-grained
diorite-gabbro); HF5 pulse/hammer (fgDG); HF6 conventional (fine-grained granite).

Why this site: P1 (convex seismicity scheduling) rests on the Dieterich rate-and-state seismicity
kernel. C4/F13 showed that the gain depends on the kernel class (exponential ≫ log ≈ sqrt). Äspö has
many short stages with very different rates/pressures inside the same test → the stage-to-stage AE
counts can discriminate kernels, which Otaniemi (one long stimulation) could not.

## Pre-registered tests

**Stage definition (fixed now).** Stage = maximal interval with flow > 0.3 L/min, gaps < 5 s merged,
duration ≥ 10 s. Event window of a stage = [start, min(next stage start, end + 600 s)). Co-injection
events = before shut-in, post events = after. Catalogue and hydraulics are both UTC.

**H6a — kernel discrimination (stage level).** Poisson likelihood of the stage counts with one free
rate factor per test plus the kernel parameters. Kernels:
- K1 Dieterich: R/r from the exact Dieterich ODE driven by S(t)=p(t) (interval pressure), t_a = 1e8 s,
  free Aσ;  N = r·∫R dt over the event window.
- K2 hydraulic energy N ∝ ∫ p q dt (the empirical Äspö measure, Zang 2017).
- K3 injected volume N ∝ V.
- K4 Kaiser (stress memory) N ∝ max(0, E_s − max previous E) with E = stage hydraulic energy.
Fit on the conventional tests HF1+HF2 (train). Decision: a kernel is "preferred" if its AIC is ≥ 2
below every other. Out-of-sample: within-test stage proportions of HF6 (conventional, other rock),
Poisson deviance per kernel, same ordering expected.
*Prediction registered:* K1 (exponential in pressure) is NOT preferred, because stage pressures
vary only 8–11 MPa while counts vary 0–50; energy/volume kernels explain more. (If K1 wins, the
exponential-kernel gains of P1 stand at mine scale.)

**H6b — is the cyclic reduction distinguishable from inter-interval variability?** Using the
preferred kernel's per-test rate factors for the conventional tests (HF1, HF2, HF4, HF6), estimate
the log-normal between-test spread σ and the predictive interval of a new conventional test in the
same rock (AG: HF1, HF2). HF3 "reduced" if its rate factor lies below the 5 % predictive quantile.
*Prediction registered:* NOT significant (n = 2 same-rock references; the claimed 6× reduction is
within inter-interval scatter once HF4 is considered).

**H1h — post-shut-in fraction.** Observed post/total fraction per stage vs the K1 Dieterich prediction
with the fitted Aσ. *Prediction registered:* Dieterich under-predicts post-shut-in events (fracture
growth continues after shut-in via pressure redistribution; Dieterich with local pressure decays
instantly when p drops) → the "terminal burst" of the P1 optimum is penalised more in reality than
in the horizon objective, supporting a shut-in-aware objective.

**H1g (follow-up).** The kernel preferred in H6a selects the gain row of the C4 table
(exponential 16.4 % / log 5.1 % / sqrt 4.8 % at Π 3) quoted as the mine-scale expectation.

Statistical caveat stated in advance: 196 events, ~20 stages; all conclusions are at the level of
"consistent / inconsistent", not effect sizes.
