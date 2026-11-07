# Research Cycle 1: optimal injection under rate-and-state seismicity (H1)

Scripts are in exp/: `h1_seis_opt.py` (main), `h1_verify.py`, `h1_dose.py`, `h1_fig.py`. Results: `h1_results.json` and `h1_fig.png`.
The full cycle runs in about 35 s on 2 cores.

## 1. Theory

The model uses Dieterich (1994) seismicity at constant effective normal stress, in the Heimisson & Segall (2018) / Heimisson (2019) form:

R/r = exp((S − τ̇_r N/r)/a),  with a = Aσ̄,  t_a = a/τ̇_r,  S = τ̇_r t + ΔS(t).

This integrates exactly to

N(T)/r = t_a ln(1 + t_a⁻¹ ∫₀ᵀ e^{S/a} ds).

The **excess events** E(T) = N(T)/r − T are then

E = t_a · ln( e^{−T/t_a} + Σ_i (w_i/t_a) e^{(s_i−T)/t_a} e^{ΔS(s_i)/a} ).

**Theorem (C1).** Suppose ΔS = 𝒢q is linear in the injection-rate history q. This holds for single-phase linear poroelasticity in heterogeneous media, with any number of wells. Then E is a log-sum-exp of affine functions of q, so **E is convex**.

The same holds for any positive-weighted sum over fault patches (a sum of log-convex functions is log-convex). With linear constraints (fixed volume, rate bounds, pressure limits), *minimum-seismicity injection design is a convex programme*:
- every local optimum is global;
- the optimal fault-stress history 𝒢q* is unique (exp is strictly convex);
- q* itself is generally **not** unique, because diffusion is a low-pass filter. 𝒢 has near-null high-frequency directions, which is also why cyclic ≈ constant at distant faults.

**KKT condition.** Write φ = softmax weights. Then 𝒢ᵀφ = μ on the interior of the support: *the marginal expected events per m³ is equal at every active injection time.* With bound constraints this gives bang–singular–bang structure, as Pontryagin predicts for a Hamiltonian that is linear in the control.

**Corollary (C2, dose invariance).** For any linear diffusion/poroelastic system whose response decays, ∫₀^∞ ΔS dt = V · S_ss(x), where S_ss is the steady response to a unit rate. This is **independent of the schedule**.

So when t ≪ t_a, E ≈ ∫(e^{ΔS/a} − 1) dt with a fixed linear dose. Schedules differ only through the Jensen gap, which is controlled by **Π = ΔS_max/(Aσ̄)**:
- Π ≪ 1: the schedule is irrelevant (second-order effect);
- Π ≳ 1: the schedule matters strongly.

## 2. Pre-registered criteria (written before the runs)

| Criterion | Pass if | Result |
|---|---|---|
| Integral form matches the original γ-ODE | rel < 1e-4 | **1e-5 to 2e-5** ✔ |
| No convexity violations on random pairs | 0 of 2000 | **0** ✔ |
| Adjoint gradient = FD | 4 digits | **7 digits** ✔ |
| Multistart reaches the same optimum | ΔE < 1e-4 rel | **3e-6** ✔ (q differs by up to 46%: non-unique q, as predicted) |
| KKT equal-marginal condition | spread < 1e-3 | **5e-4** ✔ |
| Dose invariance (C2) | spread < 0.5% | **0.1%** ✔ (87.0% inside the horizon for every schedule; the rest is the analytic tail) |
| Window artefact explains the literature | the ranking flips between co-injection and full-tail windows | **yes** ✔ |
| Kill criterion | convexity violated, or the optimum depends on the start | not triggered |

## 3. Results

Setup: 3D point source; k = 10 mD; c = 0.1 m²/s; fault at 300 m; V = 8640 m³ over 10 d; Aσ̄ = 0.05 MPa; t_a = 100 yr.

| Evaluation window | step-up | ramp-up | step-down | taper | cyclic | **optimal** |
|---|---|---|---|---|---|---|
| During injection only | **−25%** | **−32%** | +37% | +27% | +7% | −63% |
| + post-shut-in tail (2–20 Ti) | +3 to +4% | +6 to +7% | +3 to +6% | +3 to +5% | +0.2 to +0.9% | −4 to −6% |

All percentages are induced events relative to constant-rate injection.

**Finding 1 (reconciliation).** Step-up "reduces seismicity" only when the count stops at shut-in. It defers events into the post-shut-in tail. Once the tail is counted, the ranking reverses.

This matches the observation that peak rates often occur after shut-in (Segall & Lu 2015). It offers a single explanation for why Alghannam & Juanes (during-injection metrics) and Boyet (cumulative) disagree. That explanation is a model-based hypothesis, still to be checked against their exact set-ups.

**Finding 2 (sensitivity law).**

| Π | 0.19 | 0.38 | 0.76 | 1.5 | 2.6 | 3.8 | 7.6 |
|---|---|---|---|---|---|---|---|
| Optimal gain vs constant | 0.6% | 1.0% | 2.0% | 4.7% | 9.8% | 18% | 47% |
| Worst standard schedule | +0.7% | +1.0% | +2.1% | +6.9% | +18% | +42% | +176% |

**Finding 3.** With several patches (150/300/600 m), the optimum gains 19%. Step-up costs +105% and cyclic +21%: near-well faults "see" the cycling. With strong population depletion (t_a = 0.05 yr), front-loading becomes slightly better (clock-advance regime).

**Finding 4.** The global optimum takes about 1 s (SLSQP, 100 controls). Published approaches need stochastic gradients (100 perturbations per iteration) or LP surrogates.

## 4. Rethink: limitations and threats to validity

- **Physics:** the Dieterich population model assumes non-interacting sources near failure and no magnitudes; it counts expected events. Linearity fails with permeability enhancement, multiphase or thermal nonlinearity, and slip-induced stress transfer. Constant-σ̄ Coulomb approximation: Heimisson & Segall show small differences when normal stress varies.
- **The objective sets the answer:** event count, peak rate and Mmax probability give different optima. Max-rate is not proven convex (Cycle 2).
- **Poroelastic stress kernels are not yet included** (they are linear, so C1 still holds). Segall–Lu tapering benefits come from the poroelastic shut-in jump, so tapering is not judged here.
- **Novelty:** the convexity statement and the window-artefact explanation were not found in a targeted search (GAP_MAP.md). A deeper search is needed before claiming priority.

## 5. Decision

**PROMOTE.**
- **Cycle 2 (H1b):** add poroelastic kernels, multiple wells, max-rate / Mmax objectives, and a replication of the Alghannam–Juanes and Segall–Lu set-ups.
- **C20 (H1c):** a ResSim "Seismicity-aware injection designer". Engine Green's functions at fault cells (unit pulses or adjoint) feed the convex optimiser, with a UI to compare schedules. The result can be exported as WCONINJE/WCONHIST.
