> **SUPERSEDED (C50, 2026-11-04).** The submission version of P1 is `grl/P1_GRL.tex` (+ SI). This draft is kept for its derivations only. Its Äspö section moved to P3 (`../P3_aspo_kernels`), and its Äspö numbers used the defective catalogue times (F39). Do not submit or quote from this file.

# Minimum-seismicity injection design is a convex problem: global optima, the counting-window artefact, and when the schedule matters

*Draft v0.5, 2026-10-21. Authors: [user / team to decide]. AI assistance must be disclosed per journal policy.*

> **Status (honest).** Mathematics (Section 2 and Appendix C), synthetic experiments, two field replications (Otaniemi 2018 and
> the pre-registered Otaniemi 2020 test), the Äspö scope test and the simulator experiments are complete and reproducible. The
> pre-registered 2020 test **refuted** our own quantitative prediction for the gain (≤ 3 % predicted, 9.1 % realised) while
> **confirming** the θ-curve; the θ-law (F10) is therefore reported as a *candidate supported post hoc* (forecast 6.5 % vs 9.1 %),
> not as verified. Before submission: Google Scholar novelty check, author list, Zenodo archive DOI. All figures were rendered by the
> publication figure engine of ResSim Studio (`app/src/figure.cpp`), driven by `tools/paper_figs.cpp`; Figures 6–7 are ResSim runs.

## Abstract
Operators can modify injection rates to limit induced seismicity, yet published studies disagree on which protocol is safest. We show
that within the Dieterich (1994) seismicity-rate framework, with Coulomb stress linear in the injection rates (linear poroelasticity, any
heterogeneity, any number of wells and fault patches), the expected number of events over any finite window is a log-integral-exp of affine
functions of the rate history and therefore convex. Minimum-seismicity design is then a convex programme with a certified global optimum,
the optimal stress history on each fault is unique, and peak-rate minimisation reduces to a linear programme. A stress-dose invariance turns
schedule effects into Jensen gaps governed by Π = ΔS_max/(Aσ̄), and the value of design also depends on θ, the ratio of the stress-response
time to the injection time. In synthetic tests the ranking of protocols reverses with the counting window: a four-step step-up schedule gives
25 % fewer events than constant rate when counting stops at shut-in but 4 % more when the post-injection tail is counted. At the Otaniemi 2018
stimulation (θ ≪ 0.1) the global optimum gains only 1.6 % over constant rate, whereas the operator's irregular schedule cost 11 %. In a
pre-registered test on the 2020 Otaniemi stimulation, which was three times shorter (θ ≈ 0.3), the gain rose to 9.1 % — refuting our
prediction of ≤ 3 % but confirming the predicted peak of the gain at θ ≈ 0.1–1; a post-hoc θ-law forecast from the 2018 data gives 6.5 %. At the
Äspö hard-rock laboratory, acoustic emissions of hydraulic fracturing follow the injected volume, not the Dieterich response to pressure
(ΔAIC = 66), which delimits the theory to shear reactivation by pore-pressure diffusion. Finally, the framework is integrated into a reservoir
simulator, where fault Coulomb stress, Dieterich seismicity and reactivation-dependent fault transmissibility are updated during the flow run;
in a test case the two-way coupling raises the expected event count by 36 % relative to one-way monitoring.

## Key points
1. Under Dieterich seismicity with linear poroelastic coupling, the expected number of induced events is convex in the injection-rate histories of all wells, so optimal injection design has a certified global optimum.
2. Protocol rankings depend on the counting window; step-up schedules can flip from best to worse once post-shut-in events are counted.
3. Whether the schedule matters is set by Π and θ: at Otaniemi 2018 (θ ≪ 0.1) optimisation gains 1.6 %, at the shorter 2020 stimulation (θ ≈ 0.3) 9.1 %. Hydraulic-fracture AE (Äspö) lie outside the theory.

## 1. Introduction
Induced seismicity constrains enhanced geothermal systems, CO₂ storage and wastewater disposal. Recommended protocols conflict: constant
rate (Boyet et al. 2024), gradual step-up (Alghannam & Juanes 2020), cyclic soft stimulation (Zang et al. 2018; Hofmann et al. 2019) and
shut-in / taper strategies (Segall & Lu 2015). Existing optimisation approaches use constrained optimisation of a Dieterich-rate forecast (Candela et al. 2022),
stochastic gradients, linear programmes on the seismogenic index (Hill et al. 2024) or feedback tracking control (Gutiérrez-Oribio & Stefanou 2024, 2026; Walter et al. 2026). For Otaniemi itself, Kim et al. (2025) combine a robust tracking controller with a scenario search and find stimulations that outperform the 2018 schedule; their result is consistent with our finding that the operator's schedule cost 11 %, but it carries no global-optimality certificate. None establishes convexity or global
optimality, and none explains why the recommendations disagree. This paper (i) proves convexity and derives its consequences
(Section 2), (ii) tests the predictions in synthetic experiments, two field data sets and a reservoir simulator (Section 3), and
(iii) states where the theory applies and where it does not (Sections 4–6).

## 2. Mathematical framework

### 2.1 Seismicity model
For a fault patch with constant effective normal stress σ̄, direct-effect parameter A and background stressing rate τ̇_r, the
Dieterich (1994) seismicity rate R satisfies, in the state form of Heimisson & Segall (2018),

  dγ/dt = (1 − γ Ṡ)/a,  R = r/(γ τ̇_r),  a = Aσ̄,  t_a = a/τ̇_r,      (1)

with Coulomb stress S(t) = τ̇_r t + ΔS(t) and background rate r. For any S(t) the solution integrates in closed form,

  N(T)/r = t_a ln( 1 + t_a⁻¹ ∫₀ᵀ e^{S(s)/a} ds ),      (2)

and the excess over the background is E(T) = N(T)/r − T. Linear poroelasticity (any heterogeneity, any number of wells w) gives

  ΔS(x, t) = Σ_w ∫₀ᵗ g_w(x, t − s) q_w(s) ds = (𝒢q)(x, t),      (3)

where g_w is the Coulomb stress response at x to a unit rate step at well w (computed analytically or by a reservoir simulator).

### 2.2 Theorem 1 (convexity)
*Let ΔS = 𝒢q with 𝒢 linear. Then for every T > 0 and patch x, E(T; q) is a convex function of the rate history q. A positive-weighted
sum over patches is convex. With linear constraints (volume, rate bounds, pressure caps that are linear in q), minimum-seismicity design is a
convex programme.*

*Proof.* Write the argument of the logarithm in (2) as an integral over the measure μ = δ_* + t_a⁻¹ ds on {*} ∪ [0, T], where the atom *
carries the integrand e⁰ = 1:

  N/r = t_a ln ∫ e^{h_s(q)} dμ(s),  h_*(q) = 0,  h_s(q) = (τ̇_r s + (𝒢q)(x, s))/a.

Each h_s is affine in q. For q₁, q₂ and λ ∈ (0, 1), Hölder's inequality with exponents 1/λ and 1/(1 − λ) gives

  ∫ e^{λ h_s(q₁) + (1−λ) h_s(q₂)} dμ ≤ (∫ e^{h_s(q₁)} dμ)^λ (∫ e^{h_s(q₂)} dμ)^{1−λ}.

Taking logarithms, ln ∫ e^{h(λq₁+(1−λ)q₂)} dμ ≤ λ ln ∫ e^{h(q₁)} dμ + (1 − λ) ln ∫ e^{h(q₂)} dμ, i.e. N/r is convex in q, and E = N/r − T
differs by a constant. Non-negative sums preserve convexity. ∎

*Corollary 1 (unique stress history).* Because the atom fixes h_* = 0, ln ∫ e^{h} dμ is strictly convex in the function h on the support of
μ (equality in Hölder requires h(q₁) − h(q₂) to be constant μ-almost everywhere, and the atom forces the constant to be 0). Hence the optimal
Coulomb stress history 𝒢q* is unique, while q* itself is unique only if 𝒢 is injective; diffusion filters high frequencies, so in practice many
rate histories share the optimal stress history.

*Corollary 2 (optimality conditions).* With the volume constraint Σ ∫ q dt = V and bounds 0 ≤ q ≤ q_max, the KKT conditions state that the
marginal number of events per unit injected volume, ∂E/∂q_w(t) = (t_a/a) ∫ g_w(x, s − t) ω(s) dμ(s) with the softmax weights
ω(s) = e^{h_s}/∫ e^{h} dμ, is equal on the set where the bounds are inactive (equal marginal events per unit volume), which yields bang–singular–bang structures under rate bounds.

### 2.3 Proposition 2 (stress-dose invariance and Jensen gaps)
*If ∫₀^∞ g(x, t) dt = S_ss(x) is finite (the steady-state stress per unit rate), then ∫₀^∞ ΔS(x, t) dt = V S_ss(x) for every schedule with
volume V.*

*Proof.* By Fubini, ∫₀^∞ ∫₀ᵗ g(t − s) q(s) ds dt = ∫₀^∞ q(s) ∫₀^∞ g(u) du ds = V S_ss. ∎

For times short compared with t_a, (2) gives E ≈ ∫ (e^{ΔS/a} − 1) dt. With the time integral of ΔS fixed by Proposition 2, Jensen's
inequality shows that the flattest admissible stress history minimises E, and the gap between schedules scales with the variance of ΔS/a,
i.e. with Π = ΔS_max/a. Schedules matter only when Π is of order one or larger.

### 2.4 Proposition 3 (peak rate is a linear programme)
*For t ≪ t_a, max_t R ≈ r exp(max_t ΔS(t)/a). Minimising the peak rate is equivalent to minimising max_t (𝒢q)(t), a linear programme in
(q, z): min z subject to (𝒢q)(t_i) ≤ z for all evaluation times t_i and the linear constraints on q.* The ε-constraint problem
min E(q) s.t. max ΔS ≤ z is convex for each z, so the count / peak-rate Pareto front is convex.

### 2.5 Structure of the optimum and the θ-law
For a single first-order kernel T_k ṡ = u − s and the post-injection tail included in the objective, Pontryagin's principle gives a
bang–singular–bang optimum: a maximum-rate charging arc, a constant-stress plateau (singular arc, ∂H/∂u ≡ 0) and a maximum-rate terminal burst
at shut-in (Fig. 1a). The exact three-arc cost yields the gain of the optimum over constant rate, G(Π, θ, u_max), with θ = T_k/T_inj; for θ → 0,
G = θ κ₂(Π, u) + O(θ²) (Appendix C). The gain vanishes for θ → 0 (the stress follows the rate, nothing to shape) and for θ → ∞ (the
stress integrates the volume), and peaks near θ ≈ 0.1–1.

### 2.6 Coupling to a reservoir simulator
*Design.* Step responses g_w are computed by one simulator run per well; the convex programme is solved by FISTA; the design is re-simulated
in the nonlinear simulator and the remaining model error is removed by a convex offset correction (the additive mismatch δ(t) of the stress
history enters (2) as h_s → h_s + δ(s)/a, which preserves convexity).

*Monitoring inside the flow run (new, C24).* Within each report step of a simulator, the Coulomb stress of each fault face changes linearly,
so Ṡ is constant and (1) has the exact solution

  γ(t + h) = (γ(t) − 1/Ṡ) e^{−Ṡh/a} + 1/Ṡ,      (4)

and N(t) = r k [t + τ ln((1 + (k − 1)e^{−t/τ})/k)] for k = Ṡ/τ̇_r, τ = a/Ṡ from a background steady state. The simulator evaluates the
total stress on each fault face from an in-situ stress field with a poroelastic stress path, CFF = τ − μ(σ_n − p) − c, advances (4) per face,
samples a stochastic catalogue (Poisson counts, truncated Gutenberg–Richter magnitudes) and — for two-way coupling — multiplies the
transmissibility of the connections crossing reactivated faces (CFF ≥ 0) before the next report step.

## 3. Experimental study

### 3.1 Design
Each cycle was pre-registered (hypotheses, tests, pass/fail criteria) before the data were examined (research/CYCLE*_PLAN.md). Findings are
logged with their verdicts, including failures (research/FINDINGS_LOG.md).

### 3.2 Synthetic experiments
3-D point-source pore-pressure diffusion with exact erfc step responses; permeability 10⁻¹⁴ m², viscosity 10⁻³ Pa s, diffusivity
0.1 m²/s, friction 0.6, Aσ̄ = 0.05 MPa, t_a = 100 yr, one patch at 300 m, injected volume 8640 m³ in T_inj = 10 days, rate bound three times
the mean rate, 100 piecewise-constant control intervals, counting window up to 5 T_inj. Schedules: constant, four-step step-up, taper of the
last 30 %, and the global optimum (SLSQP with the analytic adjoint gradient; the softmax weights of the log-integral-exp).

### 3.3 Verification protocol
ODE (1) against the integral (2) to 10⁻⁵; 2000 random midpoint convexity tests; gradient against finite differences; multistart from random
schedules; KKT residual; dose invariance to 0.1 %; in the simulator: the in-run Dieterich update (4) against the closed form for a constant
stressing rate and the resolved fault stresses against a hand calculation.

### 3.4 Field data I: Otaniemi 2018 (OTN-3)
Public pumping data (Kukkonen 2023) and the GFZ catalogue (Leonhardt et al. 2021), both CC BY 4.0: 18,433 m³ in 49 days, 23,854 events
above M_c = −0.95, b = 1.28. Pre-registered tests: E1 a linear rate → wellhead pressure kernel on a 40 % hold-out; E2 Dieterich fit and the 63-day
post-injection tail; E3 temporal hold-out forecast; E6 design from the hold-out calibration.

### 3.5 Field data II: Otaniemi 2020 (OTN-2), pre-registered
Same pipeline, fixed before any fit (research/CYCLE6_PLAN.md): pumping data OTN2.xlsx (Kukkonen 2023), catalogue of
Kwiatek et al. (2022; doi 10.5880/GFZ.4.2.2022.001). 2,884 m³ in 16.1 days, 2,846 events above M_c = −1.45, b = 1.39. Predictions
registered in advance: E1′ kernel hold-out error ≤ 10 %; E2′ Dieterich fit WEAK; E3′ FAIL; E9 (new) transfer of the 2018 parameters;
**E4′ gain of the optimum ≤ 3 %** (assuming θ ≪ 0.1 as in 2018); **E8′ θ-scan peak in [0.05, 1]**. F10 would be called verified only
if E4′ and E8′ both met their predictions.

### 3.6 Field data III: Äspö HRL 2015 (scope test)
Six hydraulic-fracturing tests (HF1–HF6) with 1 s flow and interval-pressure records (Zang et al. 2017; ESSD 2024 data set) and the 196-event
AE catalogue (Kwiatek et al. 2018). Stage-level counting models (volume, energy, Dieterich response to the interval pressure, duration,
constant, and volume with rate / cumulative-volume factors) compared by AIC; pressure filter p_med ≥ 1 MPa.

### 3.7 Simulator experiments
(a) Seismicity-aware design in ResSim on a 41 × 41 single-phase water case. (b) The coupled fault case `examples/FAULT_INJECTION.DATA`:
30 × 20 × 5 cells of 50 × 50 × 10 m at 1975 m depth, fault F1 (100 faces, MULTFLT 0.05) 150 m east of a water injector (800 sm³/d,
280 bar cap, shut in at 500 days), producer east of the fault; Sv gradient 0.24 bar/m, Shmin/Sv 0.65, SHmax/Sv 1.10 at 60°, stress path 0.5,
μ = 0.6, t_a = 100 yr, r₀ = 5 events km⁻² yr⁻¹ (M ≥ −1), b = 1; two-way transmissibility factor 20 versus one-way monitoring; a
rate-and-state trace (a = 0.010, b = 0.015, D_c = 1 mm) on layer 3.

## 4. Results

![Figure 1](figs/fig1_window_artefact.pdf)

**Figure 1.** Synthetic test. (a) Equal-volume schedules; the optimum has the bang–singular–bang structure. (b) Coulomb stress on the patch.
(c) Cumulative excess events; the dashed line marks shut-in. (d) Events relative to constant rate when counting stops at shut-in (light) and
when the post-injection tail is counted (dark): 1 constant, 2 step-up, 3 taper, 4 global optimum.

**Counting-window artefact.** Relative to constant rate (16.36 / 38.83 days of background-equivalent events at shut-in / total), the step-up
schedule gives 12.21 (−25.3 %) at shut-in but 40.49 (+4.3 %) in total; the taper gives 20.81 (+27.2 %) and 40.16 (+3.4 %); the global optimum
gives 14.14 (−13.5 %) and 37.01 (−4.7 %), the lowest total of all admissible schedules. The "best" protocol therefore depends on when counting
stops. In 26 % of randomly sampled sites (Cycle 1 Monte Carlo) step-up flips from best to worse once post-shut-in events are counted.

![Figure 2](figs/fig2_convexity_theta.pdf)

**Figure 2.** (a) E along straight segments between random feasible schedules lies on or below the chord (Theorem 1); 2000 random midpoint
tests gave no violation. (b) Otaniemi: change of the expected count relative to constant rate for the global optimum and for the operator's
schedule as a function of θ (kernel time rescaled).

**Convexity.** No midpoint violation in 2000 tests; multistart optima agree to solver tolerance; KKT residuals vanish; the dose invariance holds to
0.1 %.

![Figure 3](figs/fig3_otaniemi.pdf)

**Figure 3.** Otaniemi 2018. (a) Injection rate and the optimal schedule with the same volume. (b) Wellhead pressure and the linear kernel (all
data and first 60 %). (c) Hourly event rate, Dieterich fit and hold-out forecast. (d) Cumulative events. Vertical lines: end of the training
window and shut-in.

**Otaniemi.** E1 PASS (hold-out WHP error 6.1 % of the range). E2 WEAK: the in-sample fit is excellent (correlation of cumulative counts 0.998)
but the 63-day tail is over-predicted by 50 % (2756 vs 1839 events). E3 FAIL: the hold-out forecast misses the end count by +29 %. At the same
volume and rate bounds the expected counts are 23,828 (actual schedule), 21,546 (constant rate) and 21,211 (optimum), Π = 3.06; the optimum gains
1.6 %, an order of magnitude below the Π-law estimate (13 %). Rescaling the kernel time (Fig. 2b) shows that the gain peaks at θ ≈ 0.1–1
(−9.2 % at θ = 0.26) and that the Otaniemi stress-carrying modes have θ ≪ 0.1, where shaping gains little but irregular schedules cost much
(+28 % at θ = 0.026, +10.6 % for the fitted kernel). E6 PASS: designs from the hold-out calibration keep 70 % of the gain.

![Figure 4](figs/fig4_otaniemi2020_theta.pdf)

**Figure 4.** Otaniemi 2020 (OTN-2) and the θ-law at two sites. (a) Injection rate and the optimal schedule with the same volume.
(b) Cumulative events and the Dieterich fit. (c) Gain of the global optimum over constant rate versus θ = T_k/T_inj (kernel time
rescaled) for both stimulations. (d) Extra events of the operator's schedule relative to constant rate.

**Otaniemi 2020 (pre-registered).** E1′ WEAK (hold-out error 11.4 % vs the 10 % criterion); E2′ PASS (correlation 0.998, tail +8 %,
P* = 72.9 bar, τ_d = 10 h); E3′ FAIL as predicted (the hold-out window is degenerate, P* at its bound); E9 PASS (the 2018 parameters
transfer with correlation 0.985). **E4′ REFUTED:** 2,834 (actual), 2,286 (constant), 2,078 (optimum) expected events — the optimum gains
**9.1 %**, three times the predicted bound, and the operator schedule costs +24 %. **E8′ CONFIRMED:** the θ-scan gain is 2.4 / 7.6 / 14.0 /
11.5 / 5.8 % at θ = 0.026 / 0.078 / 0.26 / 0.77 / 2.58, peaking inside the predicted interval (Fig. 4c). The prediction failed because its
input was wrong: the 2020 injection was three times shorter, so the same crustal response times give θ ≈ 0.26–0.77 (73 % of the fitted
kernel's steady-state pressure), the peak of the θ-curve. Weighting the 2018 θ-curve by the 2020 kernel shares gives a *post-hoc*
forecast of 6.5 % against 9.1 % realised (Π-law: 15.6 %); we report it as consistent, not as a confirmation. With the operator's eight
pauses imposed, constant rate inside the active windows gives 3,816 events (worse than the operator) and the optimum 2,747 (−3.1 %).
E6 FAIL (2.9× constant): a design from a degenerate calibration is dangerous — design only with an identifiable P*.

![Figure 5](figs/fig5_aspo.pdf)

**Figure 5.** Äspö HRL. (a) AE events per injection stage versus injected volume for HF1–HF6. (b) AIC of the stage-count models relative to
the best: 1 V(1+V_cum)^γ, 2 V q^β, 3 V q^β (1+V_cum)^γ, 4 hydraulic energy, 5 volume, 6 Dieterich response to pressure, 7 duration, 8 constant.

**Äspö.** The best model is N ∝ V (1 + V_cum)^0.45 (AIC 71.5); hydraulic energy (80.4) and volume (86.3) follow; the Dieterich response to the
interval pressure is rejected (137.2, ΔAIC = 66 against the best). At fixed volume the rate exponent is unresolved (β = −0.10, 95 % interval
[−1.25, 0.75]). The post-shut-in fraction is half the Dieterich prediction (0.25 vs 0.48). The reported reduction of AE by the cyclic protocol
of HF3 (0.25 × the conventional count) lies within the scatter between intervals.

**Simulator design (Table 1).** On the heterogeneous 41 × 41 case the designed schedule reduces the simulated expected events by 18.0 %; the
uncorrected Green's-function model is biased by 17.8 % because small stress errors are amplified exponentially (Π ≈ 9), and one convex offset
correction reduces the model error to 0.15 %.

| Table 1 | Constant rate | Designed |
|---|---|---|
| Expected events, design model | 127.9 | 100.2 |
| Expected events, simulated | 143.8 | 117.8 (−18.0 %) |
| Peak rate / background, model | 1.23e4 | 6.11e3 |
| Model error | 17.8 % (first pass) | 0.15 % (after one correction) |

![Figure 6](figs/fig6_coupled_fault.pdf)

**Figure 6.** Fault mechanics inside the reservoir simulator (ResSim, `examples/FAULT_INJECTION.DATA`), two-way coupled (reactivated faces
conduct) versus one-way monitoring. (a) Maximum Coulomb stress change on F1. (b) Reactivated fault area. (c) Expected and sampled events.
(d) Pressure jump across F1 in the injector row. Dashed vertical line: shut-in.

![Figure 7](figs/fig6_coupled_fault_mw.pdf)

**Figure 7.** Stochastic catalogue of the two-way run: (a) magnitudes versus time, (b) magnitude–frequency distribution with the input b = 1.

**Coupled fault run.** The resolved initial stresses reproduce the hand calculation to 10⁻³ bar (CFF = −4.8 bar at 1980 m) and the in-run
Dieterich count matches the closed form for a constant stressing rate (1.32734 vs 1.32734 events on a test face). In the injection case the
maximum ΔCFS reaches 14.6 bar with two-way coupling and 11.2 bar one-way; 43 % of F1 is reactivated after 100 days and all of it after 200 days
in both runs. Two-way coupling turns the sealing fault into a conduit (the pressure jump across F1 falls from 43.7 bar to ≈ 5 bar by shut-in,
while it stays at 47 bar one-way) and raises the expected count from 352 to 479 events (+36 %), because the conducting fault spreads pressure
over a larger part of the fault plane. 0.3 % of the events occur after shut-in: the Coulomb stress drops below its initial value within about
100 days. The sampled catalogue (483 events, largest Mw 2.27, cap 2.69) follows the input Gutenberg–Richter law (Fig. 7). The rate-and-state
trace on layer 3 (h* = 920 m, L_b/Δx = 25.7) produced no dynamic rupture: the 1 km trace is shorter than the nucleation length, consistent
with aseismic slip under rate-and-state friction.

**Multiphase CO₂ injection with a coupled fault (Cycle 9, F25–F27).**
* *Setup.* ResSim CO2STORE (40 × 1 × 10 cells, dx 10 m) with a FAULTMECH fault 60–300 m from the injector; four 90-d
  periods (Σu = 4, 0.25 ≤ u ≤ 1.75) and a 360-d shut-in.
* *Convexity.* The Dieterich count stayed midpoint-convex in all 320 tested pairs (8 cases: CO₂ and water, near and far fault,
  and a plume crossing the fault). The worst δ was −0.037 %.
* *Why it survives.* The gap splits exactly into δ = NL + CV. CV ≤ 0 is the Jensen term of Theorem 1, and NL is the cost of
  the nonlinear rate → stress map. NL was mostly negative (the CO₂ map leans slightly convex). Where it was positive, it was
  at most 18 % of abs(CV) (median abs(NL)/abs(CV) 1–6 %).
* *Design.* Sequential convex programming on the pulse-linearised history (5 runs per iterate) reached the optimum in
  11 simulations. The result was 1.25 % below the best point of a 231-run exhaustive grid, and both starting points agreed
  to 2·10⁻⁷. This is the behaviour Section 2.6 predicts when the linearisation error is small compared with the convexity
  margin.

## 5. Discussion
*When does design pay?* Compute Π and θ. Fast-responding reservoirs (θ ≪ 0.1) need a smooth schedule — irregularity is expensive — but gain
little from optimisation; θ ≈ 0.1–1 gives up to about 10 %. This reconciles small field gains with large synthetic ones.

*Why do protocol studies disagree?* Rankings flip with the counting window (Fig. 1d). Studies that stop counting at shut-in favour step-up
schedules; studies that count the tail favour constant or optimised schedules. Reports should state the window and the post-shut-in hazard.

*Kernel class.* The closed form is exact for saturating responses; diffusive responses reach 30–75 % of it, which is what the Otaniemi slope
shows (0.70 vs 1.21). Report Π, θ and the kernel class.

*Domain of validity.* Äspö shows that tensile hydraulic-fracture AE follow the injected volume (and fracture growth), not the Dieterich response
to pressure. The theory applies to shear reactivation of pre-existing faults by pore-pressure diffusion (Otaniemi / Basel type), not to
tensile-growth AE.

*Coupled simulation.* Monitoring faults one-way assumes that reactivation does not change the flow. Fig. 6 shows that a reactivated sealing fault
changes the pressure field enough to change the expected count by a third. Because transmissibility changes make 𝒢 depend on q, Theorem 1 then
holds only within a sequential-convex iteration (re-linearise 𝒢 around the current design), which the simulator supports through the
re-simulation and offset correction of Section 2.6.

*Practical recipe.* Compute Π and θ first; optimise allocation between wells before temporal shape (allocation outweighs shaping about
eightfold in Cycle 2); use an ε-constraint to trade count against peak rate; for stimulation with a flowback budget, flow back immediately
at maximum rate; check the terminal burst against post-shut-in hazard.

## 6. Limitations and remaining work
- Dieterich population assumptions, constant σ̄ per patch, expected counts (the catalogue is a sample; no ETAS triggering, which breaks convexity).
- Permeability evolution and multiphase flow break the linearity of (3). Sequential convex programming still applies, and global
  optimality is then local to the linearisation. In the Cycle 9 CO₂ cases the nonlinearity was 1–6 % of the convexity margin and
  mostly favourable. Two-way (slip-dependent permeability) coupling has not been tested yet.
- The simulator coupling uses an analytical stress path without elastic stress transfer between fault faces, and a one-report-step lag.
- The θ-law (F10) is supported post hoc at a second stimulation of the same reservoir, not by an independent pre-registered forecast at a
  different site; a third site with a public hourly log (e.g. Utah FORGE) is the next test.
- Expected counts are compared, not realised catalogues; Π differs between the two Otaniemi stimulations (3.06 vs 3.43) and no Π correction
  was applied in the θ-law forecast.

## 7. Conclusions
Minimum-seismicity injection design under Dieterich seismicity with linear coupling is convex: its optimum is global and cheap to compute, the
optimal stress history is unique, and peak-rate design is a linear programme. Protocol rankings depend on the counting window. The value of
design is governed by Π and θ; the Otaniemi stimulation is consistent with the theory (small optimisation gain, large penalty for irregular
schedules); the shorter 2020 stimulation (θ ≈ 0.3) shows the predicted rise of the gain to 9.1 %, which refuted our pre-registered bound
but follows the θ-curve; both show the limits of a stationary Dieterich model for post-injection forecasts. Hydraulic-fracture AE at Äspö fall outside
the theory. Integrated into a reservoir simulator, the same model monitors fault stress and seismicity during the flow run and shows that
reactivation-dependent fault permeability can change the expected count substantially.

## Data and code availability
Experiment scripts: `research/exp` (Python). Figure data: `research/papers/P1_convex_seismicity/fig_data.py`, `mkfigs.py`; rendering:
`tools/paper_figs.cpp` with ResSim Studio's figure engine; Figures 6–7: `paper_figs faultrun examples/FAULT_INJECTION.DATA`. Field data:
Otaniemi 2018 (doi 10.5880/GFZ.4.2.2021.001), Otaniemi 2020 (doi 10.5880/GFZ.4.2.2022.001), hydraulics of both (doi 10.17632/ycdcw2rdfj.1), Äspö (doi 10.5880/GFZ.2.6.2023.004). To be archived on Zenodo before submission.

## References (verify every entry before submission)
- Alghannam, M., & Juanes, R. (2020). Understanding rate effects in injection-induced earthquakes. Nature Communications, 11.
- Boyet, A., Vilarrasa, V., Rutqvist, J., & De Simone, S. (2024). Forecasting fluid-injection induced seismicity to choose the best injection strategy for safety and efficiency. Phil. Trans. R. Soc. A (doi 10.1098/rsta.2023.0179).
- Dieterich, J. (1994). A constitutive law for rate of earthquake production. JGR, 99, 2601–2618.
- Gutiérrez-Oribio, D., & Stefanou, I. (2026). A Control Framework for Induced Seismicity Mitigation in Groningen Gas Reservoir. arXiv 2604.26802.
- Gutiérrez-Oribio, D., & Stefanou, I. (2024). Insights of using control theory for minimizing induced seismicity in underground reservoirs. Geomech. Energy Environ. 39, 100570. doi 10.1016/j.gete.2024.100570.
- Kim, T., Gutiérrez-Oribio, D., Stefanou, I., Acosta, M., & Avouac, J.-P. (2025). Single-well based control and optimization of hydraulic stimulation and induced seismicity: Application to the Otaniemi geothermal project. Geothermics 132, 103396. doi 10.1016/j.geothermics.2025.103396.
- Walter, L., Boyet, A., Stefanou, I., & Vilarrasa, V. (2026). Real-time forecasting and operational control of perceivable induced seismicity in geo-energies. EarthArXiv, doi 10.31223/X5P75X.
- Heimisson, E. R., & Segall, P. (2018). JGR Solid Earth, 123 (doi 10.1029/2018JB015656).
- Candela, T., Goncalves Machado, C., Leeuwenburgh, O., & Ter Heege, J. (2022). A physics-informed optimization workflow to manage injection while constraining induced seismicity: The Oklahoma case. Front. Earth Sci. 10:1053951 (doi 10.3389/feart.2022.1053951).
- Hill, R. G., Weingarten, M., Langenbruch, C., & Fialko, Y. (2024). Mitigation and Optimization of Induced Seismicity Using Physics-Based Forecasting. JGR Solid Earth (doi 10.1029/2024JB028759).
- Hofmann, H., et al. (2019). GJI.
- Kukkonen, I. (2023). St1 Deep Heat Project stimulation data. Mendeley Data (doi 10.17632/ycdcw2rdfj.1).
- Kwiatek, G., et al. (2022). Otaniemi 2020 stimulation, JGR Solid Earth (doi 10.1029/2022JB024354); catalogue doi 10.5880/GFZ.4.2.2022.001.
- Kwiatek, G., et al. (2018). JGR Solid Earth (doi 10.1029/2017JB014715); (2019) Science Advances (doi 10.1126/sciadv.aav7224).
- Leonhardt, M., et al. (2021). Otaniemi catalogue (doi 10.5880/GFZ.4.2.2021.001).
- Segall, P., & Lu, S. (2015). JGR Solid Earth (doi 10.1002/2015JB012060).
- Zang, A., et al. (2017). GJI (doi 10.1093/gji/ggw430); (2018) Rock Mech. Rock Eng. 52:475-493 (doi 10.1007/s00603-018-1467-4).
- Äspö HRL data set: ESSD 16:295 (2024), doi 10.5880/GFZ.2.6.2023.004.

## Appendix C. Three-arc optimum and the θ-law (single first-order kernel)
Normalise time by T_inj and the volume to one. Stress s obeys θ ṡ = u − s, s(0) = 0, with rate u ∈ [0, u_max], ∫₀¹ u dt = 1, and the
post-injection tail decays freely. Per unit background rate the expected excess events are J = ∫₀^∞ (e^{Π s/s_c} − 1) dt with
s_c = 1 − e^{−1/θ} (the stress reached by constant rate at shut-in, so that Π keeps its meaning).

*Arcs.* The Hamiltonian H = e^{Πs/s_c} − 1 + λ(u − s)/θ + ν u is linear in u, so u is bang–bang except on singular arcs where
λ/θ + ν ≡ 0. Differentiating ∂H/∂u ≡ 0 along the arc gives λ̇ = 0, hence (Π/s_c) e^{Πs/s_c} = λ/θ = const: the stress is constant,
s ≡ s*, and the rate equals it, u = s*. Optimal structure: (1) charging at u_max for t_r = θ ln(u_max/(u_max − s*)); (2) plateau
u = s = s*; (3) terminal burst at u_max of length t_b before shut-in (the tail weights a final stress pulse less than a plateau because
it decays at rate 1/θ after shut-in). The volume constraint fixes s* given (t_r, t_b), and the transversality condition at shut-in
fixes t_b; for u_max = 10 and θ ≳ 0.1 a coasting arc (u = 0) can precede the burst.

*Constant-rate reference.* With u ≡ 1, s = 1 − e^{−t/θ} up to shut-in and s = s_c e^{−(t−1)/θ} after. Substituting x = Π s/s_c,
J_c = (e^Π − 1) − θ e^Π Ein(Π) + θ[Ei(Π) − γ − ln Π]  (s_c → 1 for θ ≪ 1), where Ein(x) = ∫₀^x (1 − e^{−y})/y dy, Ei is the exponential
integral and γ Euler's constant (the first two terms come from the injection arc with its start-up layer, the bracket from the tail).

*Small-θ expansion.* For θ → 0 the stress follows the rate except in boundary layers of width θ at start-up and shut-in. Expanding the
three-arc cost to first order gives G = 1 − J*/J_c = θ κ₂(Π, u_max) + O(θ²) with
κ = e^Π[Π m + R − Ein(Π)]/(e^Π − 1), m = (u_max − 1) ln(u_max/(u_max − 1)),
where R collects the start-up layer. G → 0 for θ → 0 (nothing to shape) and for θ → ∞ (the stress integrates the volume).

*Verification.* `research/exp/h1f_theta_law.py` compares the closed form with a direct optimal-control solution (discretised, convex
solver): agreement to ~1e-5 for the plateau value and the cost; the first-order term θκ₂ is accurate to +5–8 % up to θ = 0.1.
The θ-scans of Figs. 2b and 4c are computed with the fitted multi-mode kernels, not with this single-mode idealisation.

