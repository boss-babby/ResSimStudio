# When does a certified-optimal injection schedule survive coupled reservoir physics? Convexity of induced-seismicity design under multiphase CO₂ flow, fault-permeability feedback and rate-and-state faults

*Manuscript draft v1.0, 2026-11-04. Target: JGR Solid Earth (alternative: Geophysical Journal International). Authors: [to be decided by the research team]. Companion paper: P1 (convexity theory, GRL, in preparation). Generative-AI assistance was used for code, analysis and drafting; this must be disclosed per AGU policy. AI tools are not authors.*

> **Status (honest, for the authors).** All experiments are complete and reproducible from `research/exp` (scripts, logs, JSON).
> Figures come from `scripts/mkfigs.py` and `scripts/fig_dose.py`. Two pre-registered tests supported their predictions (H10, H10g);
> one was refuted (H10e, against our expectation); one was non-discriminating (H10f-b). Exploratory results are labelled as such.
> Open items before submission are listed in §7.

## Abstract

For the Dieterich (1994) seismicity-rate model with a Coulomb stress that is linear in the injection rates, the expected number
of induced events is a convex function of the injection schedule (companion paper P1). Minimum-seismicity design is then a convex
programme with a certified global optimum. Real reservoirs violate the linearity assumption. We ask which kinds of coupled physics
preserve convexity and which destroy it, and test this with a reservoir simulator in which fault Coulomb stress and Dieterich
seismicity are updated during the flow run.

In a CO₂-storage section with a monitored fault, 0 of 400 midpoint-convexity tests were violated across ten configurations. These
included a plume that crosses the fault, a fault 60 m from the injector and a threefold rate. An exact decomposition shows that the
nonlinearity from multiphase flow is 1–6 % of the Jensen margin and mostly *helps* convexity. Sequential convex programming
reaches the global optimum in about 11 simulations, 1.25 % better than a 231-point grid.

Convexity breaks once a reactivated fault becomes more permeable (two-way coupling). With an abrupt permeability switch, 9–21 of
40 tests per case were violated (pre-registered, refuted). The violation frequency is governed by the transmissibility multiplier
reached on the fault. A pre-registered transfer test at a second fault position confirmed this: 0/40 and 0/40 violations at
m₁ ≈ 1.3, and 21/40 and 16/40 at m₁ ≈ 3.6. Pooled over 18 configurations, violations reach 5 % at m₁ = 1.7 (95 % interval
1.5–2.0) and 15 % at m₁ = 2.4 (2.2–2.6). The loss of convexity persists when the coupling lag is reduced fourfold; for water
injection its size shrinks from about 5 % to 1.5 %.

On a rate-and-state fault, discrete event counts (15–23 % violations) and seismic moment (32–41 %) are not convex, even after
ensemble averaging. A Gaussian-process optimiser nevertheless designs for moment within 0.2 % of an evolution strategy that uses
7.5 times more simulations.

These results give a decision rule: use certified convex design while the expected fault-permeability enhancement stays below
about 1.7×; above that, use global or Bayesian optimisation and treat the convex optimum as a starting point.

## Key points

1. Multiphase CO₂ flow and a plume crossing the fault preserve the convexity of the expected induced-event count in the injection
   schedule. Simulator-driven sequential convex programming then finds the global optimum in about 11 runs.
2. Fault-permeability feedback breaks convexity once the reached transmissibility multiplier exceeds about 1.7. A pre-registered
   test at a second fault position confirms this threshold.
3. Rate-and-state event counts and seismic moment are not convex, but Bayesian optimisation designs for them efficiently. Together
   this gives a practical decision rule for schedule design.

## Plain language summary

Injecting water or CO₂ underground can trigger earthquakes. Operators can reduce this risk by choosing how the injection rate
varies in time. In a companion paper we showed that, for the most common physics-based forecasting model, this choice is a
"convex" mathematical problem. That means a computer can find the best possible schedule and prove that nothing better exists.
Real reservoirs are more complicated than that model.

Here we use a reservoir simulator to test when the guarantee still holds. It holds when CO₂ and brine flow together and when the
CO₂ reaches the fault. It fails when the fault becomes more permeable as it starts to slip, because the fault then changes the
pressure that drives it. We measure how much permeability change is tolerable (about 1.7 times) and confirm this with a test
whose outcome we wrote down in advance. For models that simulate individual earthquakes, the guarantee does not hold, but a
learning-based optimiser still finds good schedules.

## 1. Introduction

Injection-rate management is the main operational lever against induced seismicity in enhanced geothermal systems, wastewater
disposal and CO₂ storage. Recommended protocols disagree. Proposals include constant rate (Boyet et al., 2024), gradual step-up
(Alghannam & Juanes, 2020), cyclic soft stimulation (Zang et al., 2017) and tapered shut-in (Segall & Lu, 2015).
Optimisation-based approaches have appeared:
* constrained optimisation of a Dieterich-type forecast with stochastic gradients (Candela et al., 2022);
* ensemble-based control of CO₂ injection under geomechanical constraints (Zheng et al., 2023);
* a response-matrix linear programme on Coulomb-stress-rate constraints (Hill et al., 2024).
* robust feedback control that tracks a target seismicity or pressure (Gutiérrez-Oribio & Stefanou, 2024; Kim et al., 2025), applied to Otaniemi by Kim et al. (2025).

None of these studies establishes whether its optimum is global.

Our companion paper (P1) proves the following. If the Coulomb stress on every fault patch is an affine function of the injection
rates (linear poroelasticity, any heterogeneity, any number of wells), the Dieterich expected event count over any finite window
is a log-integral-exp of affine functions of the rate history, and therefore convex. Minimum-seismicity design is then a convex
programme. Its optimum is global, the optimal stress history on each fault is unique, and a solver can certify optimality.

The weak point of this result is the affine stress map. In practice the map is produced by a nonlinear reservoir simulator. The
nonlinearities come from compressibility, relative permeability, CO₂ dissolution and buoyancy. Most importantly, faults change
their own permeability when they are reactivated (Cappa & Rutqvist, 2011; Guglielmi et al., 2015). The forecasting model itself
may also be replaced by rate-and-state fault simulations (McClure & Horne, 2011), whose event counts are discrete. A design method
used by operators must say which of these departures still allow certified design.

This paper answers that question empirically, under pre-registered hypotheses. We use a reservoir simulator (ResSim) in which
fault Coulomb stress, Dieterich seismicity and reactivation-dependent fault transmissibility are updated inside the flow run.
Section 2 summarises the theory and introduces an exact decomposition of convexity gaps. Section 3 describes the simulator, the
test design and the pre-registration. Section 4 reports five results:
* passive faults keep convexity;
* simulator-driven sequential convex programming finds the global optimum;
* fault-permeability feedback breaks convexity, with a threshold that transfers to a second geometry;
* rate-and-state objectives are not convex;
* Bayesian optimisation handles the non-convex cases.

Section 5 turns the results into a decision rule and discusses their limits.

## 2. Theory

### 2.1 Convexity of the Dieterich expected count (P1)

Let u ∈ ℝᴷ be the injection rates over K intervals. Let S_j(t; u) be the Coulomb stress history on fault patch j and r_j the
background rate. For the Dieterich (1994) model the expected number of events in [0, T] is

N(u) = Σ_j r_j ∫₀ᵀ R_j(t; u) dt,

where the state variable γ_j obeys dγ/dt = (1 − γ Ṡ_j)/(Aσ). P1 shows that, if S_j is affine in u, R_j(t; u) is a
log-integral-exp of affine functions of u. N is therefore convex for any window, any patch population and any linear stress
operator. Expected exceedance probabilities, exp(−N) transformations and peak-rate constraints inherit convexity or
quasi-convexity.

### 2.2 Midpoint test and decomposition of the gap

For schedules A and B and their midpoint C = (A + B)/2 we define the relative midpoint gap

δ = [N(C) − ½(N(A) + N(B))] / N(C).

Convexity requires δ ≤ 0. Let D_A, D_B and D_C be the simulated stress histories and D_avg = ½(D_A + D_B), and write R(D) for the
expected count of a stress history. Then exactly

δ·N(C) = [R(D_C) − R(D_avg)] + [R(D_avg) − ½(R(D_A) + R(D_B))] ≡ NL + CV.

* CV ≤ 0 is the Jensen term guaranteed by the convexity of R in the stress history (P1).
* NL measures how far the simulator's rate-to-stress map departs from linearity, in units of events.

A violation requires NL > |CV|. The ratio NL/|CV| is therefore a natural measure of how close a physical system is to losing
convexity.

## 3. Methods

### 3.1 Simulator and fault mechanics

All flow simulations use ResSim, a fully implicit black-oil and compositional simulator developed by the authors. Its CO2STORE
mode is a CO₂–brine model with dissolution, and the same keyword semantics are used as in ECLIPSE/OPM. Faults are defined by
`FAULTS`. A fault-mechanics module evaluates each fault face after every report step:
* face pressure p_f;
* total stress from a gradient stress field with a poroelastic stress path;
* the Coulomb function CFF = τ − μ(σ_n − p_f) − c;
* the Dieterich state of each face, updated exactly for piecewise-linear stress.

The module also reports the expected event count. In **two-way** mode, the transmissibility of the cell connections crossing a
reactivated face is multiplied by

m = P^x,  x = clip(CFF/(w σ'_n), 0, 1),

where P is the permeability factor and w the ramp width as a fraction of the effective normal stress σ'_n (Figure 1b). The
coupling is sequential, with a lag of one report step: the multipliers computed at the end of step n apply in step n+1.

### 3.2 Model

The model is a 40 × 1 × 10 vertical section with cells of 10 m × 100 m × 5 m (Figure 1a):
* reservoir top at 1000 m, porosity 0.25, k_h = 50 mD, k_v = 5 mD;
* a large pore-volume boundary column;
* an injector perforated in layers 8–10.

The CO₂ cases use CO2STORE at 50 °C with salinity 0.7 and a mean rate of 40,000 sm³/d. The water control cases use a
water-filled black-oil model at 150 sm³/d, which gives a similar reservoir volume rate.

A vertical fault crosses the section, at x = 150 m in geometry G1 or x = 220 m in G2. Exploratory passive cases also placed it at
60 m and 300 m. Fault parameters:
* friction 0.6;
* Aσ = 1 bar, t_a = 100 yr;
* background rate 5 events per year per face;
* stress field: S_v gradient 0.226 bar/m, S_h/S_v = 0.7, S_H/S_v = 0.95, azimuth 30°, stress path 0.5.

In the two-way cases the fault is sealing (MULTFLT 0.02) with P = 50. Without help, the default stress state keeps the fault 25–37
bar from failure. A uniform shear pre-load (rule A1, fixed before any test) is therefore applied: it is chosen so that at constant
rate the most-stressed face reaches CFF ≈ +2 bar.

A schedule is K = 4 injection intervals of 90 d with rates u_k ∈ [0.25, 1.75] and Σu_k = 4, followed by 360 d of shut-in. The
expected count N is accumulated over the whole 720 d.

### 3.3 Convexity tests and pre-registration

Each test draws M = 40 random pairs (A, B), with common random numbers across cases (seed 13 unless stated). It runs A, B and C
and records δ. A pair violates convexity if δ > ε. The tolerance ε = 3|N₁₀ − N₅|/N₁₀ is three times the relative change of the
constant-rate count when the report step is halved from 10 to 5 d. It ranges from 0.04 % to 1.4 %.

The verdict rule was stated before each test:
* SUPPORTED (convex) if ≤ 5 % of pairs violate;
* REFUTED if > 15 % violate and the median violator δ > 1 %;
* inconclusive otherwise.

Hypotheses, predictions and amendments were written to dated plan files before the runs (CYCLE9_PLAN to CYCLE12_PLAN in the data
archive). Runs added after a verdict are labelled *exploratory*. These plan files are internal project records: they were not deposited with an external registry, so their dates are not independently time-stamped. They will be archived unchanged with the code so that readers can compare plans and results.

### 3.4 Sequential convex programming (SCP)

Each SCP iterate linearises the stress history around the current schedule. This costs one base run plus four pulse runs. The
iterate then solves the convex Dieterich subproblem (SLSQP) within a trust region, and accepts or rejects the step on the true
simulated count. For comparison, an exhaustive grid over u_k ∈ {0.25, 0.5, …, 1.75} with Σu = 4 has 231 schedules.

### 3.5 Rate-and-state fault model

The discrete-event tests use a 2-D quasi-dynamic rate-and-state fault, 1.5 km long with 64 elements and a velocity-weakening
centre, solved with radiation damping (Rice, 1993). It has no tectonic loading and is forced by 2-D Theis pore pressure from the
same 4-interval schedules at fixed volume. Ensembles use E = 16 realisations of a ±20 % heterogeneous normal stress, with common
random numbers across schedules. For moment design we compare four approaches:
* a Dieterich surrogate fitted to RSF runs;
* an evolution strategy (ES; 60 evaluations × 3 seeds);
* a Gaussian-process expected-improvement optimiser (GP-EI; Jones et al., 1998), using 12 training and 12 acquisition steps;
* constant rate.

## 4. Results

### 4.1 Passive faults: multiphase CO₂ flow preserves convexity (Figure 2, left)

Under the pre-registered test H10, the four base cases had **0 of 40** violations each:
* CO₂ with the fault at 150 m (inside the plume) and at 300 m (pressure only);
* the two water controls.

The median δ was −6.9 %, −0.93 %, −9.2 % and −1.2 %, and the worst δ −0.44 %, −0.04 %, −1.0 % and −0.11 %. The exploratory
extensions also had 0/40 each: a fault at 60 m, and a threefold rate at which the plume crosses the fault by day 180 and the
fault produces 11.5 expected events. The one-way controls of the two-way experiments (sealing fault and pre-load, no permeability
switching) likewise had 0/40.

In total, 0 of 400 passive-fault tests were violated.

**Decomposition.** We decomposed 10 pairs per case (Section 2.2). The history-norm nonlinearity of the CO₂ rate-to-stress map is
0.7–1.1 %, about 10⁴ times that of water. Even so, the median NL/|CV| is only −0.010 to −0.058, and NL is positive in only 2–4
of 10 pairs. The largest harmful value, +0.18, occurs when the plume crosses the fault.

Our own pre-stated expectation was that the CO₂ map would be concave, so that NL > 0 and convexity would erode. That expectation
was refuted. The nonlinearity mainly redistributes stress in time: the signed residual Σ(D_C − D_avg)/Σ|D_C| is only −2 to
−4 × 10⁻⁴. The event functional is insensitive to such redistribution.

### 4.2 Simulator-driven sequential convex programming (Figure 5)

Starting from a deliberately poor schedule [1.75, 1.75, 0.25, 0.25] (N = 1.523), SCP reached N = 0.70729 at iteration 2, after 11
simulations. Its optimum, [1.12, 0.99, 0.89, 1.00], is 1.25 % lower than the best grid point (constant rate, 0.71624). For the
pressure-only fault, SCP gives 0.222538 against 0.222739 for the grid. Starts from the poor schedule and from constant rate agree
to 2 × 10⁻⁷, as convexity implies, and the trust-region ratio was 0.995–1.03.

Two points are practically important. First, the optimum is not constant rate, but constant rate is within 1.3 %. The small
value of design here is consistent with the θ-law of P1, because the stress-response time is short compared with 360 d.
Second, plain step-response superposition is *not* accurate enough for design: in an earlier single-phase test the event-count
error was 18 %. The reason is that the exponential sensitivity exp(ΔS/Aσ) amplifies small stress errors. The affine offset
correction used inside SCP removes this error (0.15 % after one pass) without breaking convexity.

### 4.3 Fault-permeability feedback breaks convexity (Figure 2, right)

When the reactivated fault becomes conductive with a steep switch (w = 0.05, so m₁ = 2.9 in the CO₂ case and 4.5 in water), the
pre-registered test H10e gave 9, 21, 19 and 16 violations of 40 for:
* CO₂ (NS);
* CO₂ at triple rate (NHS);
* water (WNS);
* water at triple rate (WNHS).

The median violator δ was 2.9–6.3 % and the maximum 26 %. All four cases are REFUTED as convex. The same faults without
permeability switching have 0/40 (Section 4.1), so the cause is the feedback itself, not the pre-load or multiphase flow.

The decomposition confirms this. NL/|CV| rises from about −0.03 in the passive cases to +0.73 to +0.89, and NL is positive in 8–10
of 10 pairs. Physically, the stress map now depends on the schedule through the fault's own transmissibility, S = 𝒢(u)u. A
schedule that reactivates the fault earlier spreads pressure along the fault plane earlier. That raises the stress on other
faces, and the midpoint of two schedules can reactivate more fault area than either endpoint.

### 4.4 The threshold is set by the reached transmissibility multiplier and transfers to a second geometry (Figure 4)

In the parameterisation of Section 3.1, the ramp width w and the permeability factor P enter only through the log-slope
k = ln P/(w σ'_n), as long as the cap x = 1 is not reached. A sweep in P is therefore *identical* to a sweep in w and cannot
separate "abruptness" from "amplitude". We varied w and report the multiplier m₁ actually reached on the fault at constant rate:
* w = 0.5 (pre-registered H10f-a; m₁ = 1.17 and 1.26): 0/40 and 0/40;
* exploratory w = 0.15 and 0.08 (m₁ = 1.58–2.92): violations appear only above m₁ ≈ 2.

To test whether this threshold belongs to the coupled problem or to one fault position, we pre-registered H10g. The fault was
moved to x = 220 m (G2). The pre-load and the ramp widths were calibrated by fixed rules on constant-rate runs only, with targets
m₁ = 1.4 (LOW) and m₁ = 3.5 (HIGH). The predictions were ≤ 5 % violations for LOW, and > 15 % with a median violator δ > 1 %
for HIGH.

All four predictions held:
* LOW: CO₂ 0/40 (m₁ = 1.31) and water 0/40 (m₁ = 1.35);
* HIGH: CO₂ 21/40 (m₁ = 3.60, median violator δ +4.6 %) and water 16/40 (m₁ = 3.61, +5.1 %).

H10g is SUPPORTED.

**Pooled dose-response.** We fitted a binomial model, logit p = a + b ln m₁, to all 18 configurations (720 pairs; Table 1). This
includes exploratory intermediate points in G2 and a second pair seed for two G1 cases. The fit gives b = 3.7, and the
violation frequency reaches:

| violation frequency | m₁ (fit) | 95 % interval (binomial bootstrap) | 95 % interval (quasi-binomial, φ = 1.9) |
|---|---|---|---|
| 5 % | 1.73 | 1.53–1.93 | 1.46–2.01 |
| 15 % | 2.40 | 2.25–2.56 | 2.20–2.63 |
| 50 % | 3.84 | 3.56–4.23 | 3.46–4.39 |

*Table 1. Reached transmissibility multiplier at which the convexity-violation frequency reaches 5, 15 and 50 %.* Data:
`figs/fig4_dose_response.csv`.

Three caveats apply:
* Adding a geometry term improves the fit only marginally (likelihood-ratio p = 0.058; G2 slightly more prone, odds ×1.5).
* A fluid term gives p = 0.08.
* The residual deviance (30.4 on 16 dof) indicates overdispersion. We therefore also report quasi-binomial intervals.

The second seed changes the counts by 2 of 40 (3 → 5 and 12 → 14), which is consistent with binomial noise. Because rule A1
fixes the peak CFF at about +2 bar in both geometries, the experiment cannot distinguish "m₁ governs" from "k governs". We state
this limitation, which we recorded before the pair runs.

**Coupling lag.** The two-way coupling is explicit with a lag of one report step (10 d). Could the loss of convexity be an artefact of this lag? We re-ran the H10g HIGH cases with
shorter report steps (exploratory E12c/E12d, after the verdict). For each step length we used all 40 pairs, each with its own
tolerance.

| case | report step | violations / 40 | median violator δ |
|---|---|---|---|
| CO₂ HIGH (G2) | 10 d | 21 | 4.6 % |
| CO₂ HIGH (G2) | 2.5 d | 20 | 3.6 % |
| water HIGH (G2) | 10 d | 16 | 5.1 % |
| water HIGH (G2) | 2.5 d | 13 | 1.6 % |
| water HIGH (G2) | 1.25 d | 13 of the 16 original violators | 1.45 % |

*Table 2. Lag refinement of the H10g HIGH cases.*

For CO₂ the result is insensitive to the lag. For water, the violation frequency stays above 15 %, but the amplitude falls by a
factor of about 3.5 and converges near 1.4 %. That is still above the 1 % criterion. The lag therefore inflates the size of the
effect for water but does not create it. Quantitative amplitudes need implicit coupling (Section 5.3).

### 4.5 Rate-and-state objectives are not convex (Figure 3a, b)

On the single quasi-dynamic RSF fault, the midpoint event count exceeded the endpoint average in 7 of 30 pre-registered pairs and
30 of 200 exploratory pairs (15–23 %). None of these violations was larger than one event. This is consistent with integer
catalogues of 2–4 events, so the single-fault count result is inconclusive under our rule.

The cumulative seismic moment is clearly not convex: 12/30 and 82/200 pairs violate, with a midpoint up to 1.27 times the
endpoint average. The Dieterich expected count computed from the same pressure histories had 0 violations in 230 pairs.

Ensemble averaging over 16 heterogeneous realisations does not restore convexity. The pre-registered H7b gave 14 of 60 count
violations (23 %, REFUTED), 12 of them larger than one event in the ensemble total. Moment violations were 19 of 60. Only 16 of 60
pairs had counts that were exact multiples of E, so the non-convexity is not a quantisation artefact of identical realisations.

### 4.6 Designing for a non-convex objective (Figure 3c)

A Dieterich surrogate fitted to 12 RSF moment evaluations (R² = 0.85) proposes a schedule with 1.83 times the moment found by
the evolution strategy (pre-registered H9b, REFUTED). The convex surrogate is therefore not a safe proxy for moment.

A GP-EI optimiser on log M₀ reached 2.2193 × 10¹⁵ N m in 24 evaluations. That is 1.002 times the ES median best, which used 180
evaluations, and 0.941 times constant rate (pre-registered H11, SUPPORTED). Four other acquisition seeds gave ratios of
1.001–1.006. Four new training sets (pre-registered H11b) gave 1.000–1.060 (median 1.008). One unlucky training set got only
0.4 % below constant rate in 24 evaluations. All optima front-load the injection (u₁ ≈ 1.4–1.5).

## 5. Discussion

### 5.1 A decision rule for schedule design

The experiments sort coupled physics into two classes.
* **Smooth state dependence** keeps convexity (Section 4.1). This covers multiphase flow, dissolution, compressibility and a plume
  crossing the fault. Here the expected count is effectively convex, and simulator-driven SCP gives a certified optimum in about
  ten runs.
* **Threshold feedback** breaks convexity (Sections 4.3–4.4). This covers faults whose transmissibility jumps on reactivation.
  The violation frequency is controlled by the size of the jump: below m₁ ≈ 1.7 violations are rare (≤ 5 %), and above 2.4
  they are common.

Rate-and-state objectives with discrete events or moment (Section 4.5) are non-convex but smooth enough for Bayesian
optimisation (Section 4.6). We propose:

1. Estimate the plausible fault-permeability enhancement for the design range, for example from hydraulic tests,
   injection-pressure history or published fault-valve behaviour.
2. If it is below about 1.7×, design with SCP and report the convex certificate. A posteriori, check δ on a few random pairs.
3. If it is above, or unknown, use GP-EI or another global optimiser. Start it from the SCP solution, budget more than 24
   simulations, and restart from several training sets.
4. When the objective is moment or a discrete catalogue, never use a Dieterich surrogate as the objective. Use it only as a
   prior mean.

### 5.2 Relation to previous work

Hill et al. (2024) optimise well rates with a linear programme on stress-rate constraints derived from a poroelastic response
matrix. Our results support their linear-response assumption for passive faults and multiphase flow. They also identify when it
fails: under fault-permeability feedback the response matrix depends on the schedule. The offset correction of Section 4.2 is a
cheap way to keep a linear-response design accurate where convexity holds.

Candela et al. (2022) optimise a Dieterich/RSF forecast for Oklahoma with gradient methods. Our results imply that their optimum
is global if permeability feedback is weak, and they explain why multi-start strategies are needed if it is not.

Boyet et al. (2024) find, in a hydro-mechanical model of Basel with permeability enhancement, that constant injection outperforms
cyclic and progressive protocols. Our threshold suggests a reason why protocol rankings in such models may not transfer between
sites: whether the comparison takes place in the convex regime depends on the site-specific enhancement.

### 5.3 Limitations

1. **Geometry and dimension.** We used a 2-D section with one fault and two fault positions. The threshold may depend on the
   fault-to-injector distance, the network topology and 3-D pressure spreading. The geometry term was marginal (p = 0.058). A 3-D
   model with several faults is the obvious next test.
2. **Permeability law.** The exponential multiplier with a cap is one of many (Cappa & Rutqvist, 2011). Laws with hysteresis
   (permanent enhancement) are likely to be worse for convexity, because they add memory.
3. **Pre-load.** Reactivation was ensured by a uniform shear pre-load (rule A1). Real faults have heterogeneous stress, which
   smooths the population response and may help convexity.
4. **Explicit coupling.** The lag reduction of Section 4.4 (Table 2) shows that the effect is physical. For water injection,
   however, the amplitude at a 10-d lag is about 3.5 times the converged value. An implicit coupling is needed for quantitative
   amplitudes.
5. **Statistics.** Each configuration uses M = 40 pairs. The pooled fit is overdispersed, and pairs are shared across
   configurations (common random numbers), which makes the counts correlated.
6. **Rate-and-state model.** It is quasi-dynamic, 2-D, unloaded and has 64 elements. Counts per run are small.
7. **Schedules.** K = 4 intervals with fixed total volume. The convexity results hold for any K, but the cost of SCP grows linearly
   with K.

## 6. Conclusions

1. The convex structure of Dieterich-based injection design (P1) survives multiphase CO₂ flow, dissolution and a plume crossing
   the fault (0/400 violations). The simulator nonlinearity is 1–6 % of the Jensen margin, and mostly favourable.
2. Simulator-driven sequential convex programming finds the global optimum in about 11 simulations.
3. Fault-permeability feedback breaks convexity. The violation frequency grows with the reached transmissibility multiplier:
   5 % at about 1.7× and 15 % at about 2.4×. This threshold transferred to a second fault position in a pre-registered test.
4. Rate-and-state counts and moment are not convex, even in ensembles. Gaussian-process optimisation designs for them within
   0.2 % of a much more expensive evolution strategy.
5. The resulting decision rule tells practitioners when a certified convex design is justified and when global or Bayesian
   optimisation is required.

## 7. Open items before submission (for the authors; delete before submission)

- [ ] Author list, affiliations, CRediT statement.
- [ ] Zenodo DOI for the simulator version (ResSim C50 sources and binaries) and for the `research/exp` archive (scripts, logs,
      JSON). Mention both in Open Research.
- [ ] Decide whether a 3-D multi-fault case is needed for JGR (recommended; it is a reviewer-1 question).
- [ ] An implicit two-way coupling run for the water HIGH case (E12d: the amplitude converges near 1.4 % only at 1.25-d steps).
- [ ] P1 must be submitted first, or simultaneously, and cited as "submitted". P2 must not repeat the proof.
- [ ] Convert to the AGU LaTeX template (`agujournal2019.cls`, as in P1/grl).

## Open Research

The simulator inputs (decks generated by `h10_mkdeck.py`), all experiment scripts, logs and result files, the figure scripts and
the pre-registration files will be archived on Zenodo [DOI to be added]. The simulator version used is ResSim C50 [DOI to be
added]. No external data are used.

## References

Alghannam, M., & Juanes, R. (2020). Understanding rate effects in injection-induced earthquakes. *Nature Communications*, 11. https://doi.org/10.1038/s41467-020-16860-y

Boyd, S., & Vandenberghe, L. (2004). *Convex Optimization*. Cambridge University Press. https://doi.org/10.1017/CBO9780511804441

Boyet, A., Vilarrasa, V., Rutqvist, J., & De Simone, S. (2024). Forecasting fluid-injection induced seismicity to choose the best injection strategy for safety and efficiency. *Philosophical Transactions of the Royal Society A*, 382. https://doi.org/10.1098/rsta.2023.0179

Candela, T., Goncalves Machado, C., Leeuwenburgh, O., & ter Heege, J. (2022). A physics-informed optimization workflow to manage injection while constraining induced seismicity: The Oklahoma case. *Frontiers in Earth Science*, 10, 1053951. https://doi.org/10.3389/feart.2022.1053951

Cappa, F., & Rutqvist, J. (2011). Modeling of coupled deformation and permeability evolution during fault reactivation induced by deep underground injection of CO₂. *International Journal of Greenhouse Gas Control*, 5, 336–346. https://doi.org/10.1016/j.ijggc.2010.08.005

Dieterich, J. (1994). A constitutive law for rate of earthquake production and its application to earthquake clustering. *Journal of Geophysical Research*, 99(B2), 2601–2618. https://doi.org/10.1029/93JB02581

Guglielmi, Y., Cappa, F., Avouac, J.-P., Henry, P., & Elsworth, D. (2015). Seismicity triggered by fluid injection–induced aseismic slip. *Science*, 348, 1224–1226. https://doi.org/10.1126/science.aab0476

Gutiérrez-Oribio, D., & Stefanou, I. (2024). Insights of using control theory for minimizing induced seismicity in underground reservoirs. *Geomechanics for Energy and the Environment*, 39, 100570. https://doi.org/10.1016/j.gete.2024.100570

Heimisson, E. R., & Segall, P. (2018). Constitutive law for earthquake production based on rate-and-state friction: Dieterich 1994 revisited. *Journal of Geophysical Research: Solid Earth*, 123, 4141–4156. https://doi.org/10.1029/2018JB015656

Hill, R. G., Weingarten, M., Langenbruch, C., & Fialko, Y. (2024). Mitigation and optimization of induced seismicity using physics-based forecasting. *Journal of Geophysical Research: Solid Earth*, 129. https://doi.org/10.1029/2024JB028759

Jones, D. R., Schonlau, M., & Welch, W. J. (1998). Efficient global optimization of expensive black-box functions. *Journal of Global Optimization*, 13, 455–492. https://doi.org/10.1023/A:1008306431147

Kim, T., Gutiérrez-Oribio, D., Stefanou, I., Acosta, M., & Avouac, J.-P. (2025). Single-well based control and optimization of hydraulic stimulation and induced seismicity: Application to the Otaniemi geothermal project. *Geothermics*, 132, 103396. https://doi.org/10.1016/j.geothermics.2025.103396

McClure, M. W., & Horne, R. N. (2011). Investigation of injection-induced seismicity using a coupled fluid flow and rate/state friction model. *Geophysics*, 76, WC181–WC198. https://doi.org/10.1190/geo2011-0064.1

Rice, J. R. (1993). Spatio-temporal complexity of slip on a fault. *Journal of Geophysical Research*, 98(B6), 9885–9907. https://doi.org/10.1029/93JB00191

Segall, P., & Lu, S. (2015). Injection-induced seismicity: Poroelastic and earthquake nucleation effects. *Journal of Geophysical Research: Solid Earth*, 120, 5082–5103. https://doi.org/10.1002/2015JB012060

Zang, A., et al. (2017). Hydraulic fracture monitoring in hard rock at 410 m depth with an advanced fluid-injection protocol and extensive sensor array. *Geophysical Journal International*, 208, 790–813. https://doi.org/10.1093/gji/ggw430

Zheng, F., Jha, B., & Jafarpour, B. (2023). Controlled CO₂ injection into storage reservoirs to minimize geomechanical risks under geologic uncertainty. *SPE Annual Technical Conference and Exhibition*, SPE-215155-MS. https://doi.org/10.2118/215155-MS

[P1] Companion paper: Minimum-seismicity injection design is a convex problem (GRL, in preparation).

## Figures

![Figure 1](figs/fig1_model.pdf)

**Figure 1.** (a) The 2-D section: injector perforated in layers 8–10, fault position G1 (x = 150 m) or G2 (x = 220 m), and a
pore-volume boundary. (b) The two-way transmissibility law m = P^x with x = CFF/(wσ'_n), for P = 50 and four ramp widths; the
grey line marks m = 2.

![Figure 2](figs/fig2_midpoint_gaps.pdf)

**Figure 2.** Relative midpoint gaps δ of 40 random schedule pairs per configuration. Black bars mark the tolerance ε; the
numbers at the top give violations per pairs. Blue: passive faults (one-way). Green: two-way coupling with a smooth law
(m₁ ≈ 1.2–1.35). Red: two-way coupling with a steep law (m₁ = 2.9–4.5). Configurations are listed in `figs/fig2_midpoint_gaps.csv`.

![Figure 3](figs/fig3_rsf_nonconvex.pdf)

**Figure 3.** Rate-and-state objectives. (a, b) Ratio of the midpoint value to the endpoint average for the event count, the
seismic moment and the Dieterich expected count computed from the same pressure histories: (a) one fault, 200 pairs; (b) a
16-member heterogeneous ensemble, 60 pairs. Ratios above 1 violate convexity. The counts printed are full-precision; the plotted
moments come from 3-significant-figure logs. (c) Designing for the ensemble-mean moment: GP-EI (12 training + 12 acquisition
evaluations) against an evolution strategy (median of three seeds, 60 evaluations each), constant rate and the optimum of a
Dieterich surrogate.

![Figure 4](figs/fig4_dose_response.pdf)

**Figure 4.** Fraction of pairs violating convexity against the transmissibility multiplier m₁ reached at constant rate. Blue:
G1; red: G2. Circles: CO₂; squares: water. Filled symbols are pre-registered, open symbols exploratory. The line is the pooled
logistic fit (Table 1). The grey band marks ≤ 5 % (the convex verdict) and the dashed line 15 % (the non-convex verdict).

![Figure 5](figs/fig5_scp.pdf)

**Figure 5.** Simulator-driven sequential convex programming: expected events against the number of simulations, compared with
the best of an exhaustive 231-point grid. (a) Fault inside the CO₂ plume; (b) pressure-only fault.
