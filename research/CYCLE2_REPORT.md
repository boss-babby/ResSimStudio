# Research Cycle 2 (H1b): poroelastic kernel, peak-rate, flowback, robustness, two wells

Scripts: `exp/h1b_poro.py` and `exp/h1b_extra.py`. Results: `h1b_results.json`, `h1b_extra.json` and `h1b_fig.png`. Run time is about 2 minutes on 2 cores.

## Pre-registered questions and criteria
| # | Question | Decision rule |
|---|---|---|
| T1 | Does convexity hold, and does taper help, once poroelastic stress is included? | Multistart optimum is the same; compare taper with constant |
| T2 | Is min-peak-rate tractable? | LP surrogate (min max ΔS) within 3% of the exact smoothed optimum |
| T3 | Flowback: is it useful (a) at fixed net volume (storage) or (b) at fixed injected volume (stimulation)? | Report optimum structure |
| T4 | Is the window flip (F2) generic? | Monte Carlo over r ∈ [100, 1000] m, c ∈ [0.01, 1] m²/s, Aσ ∈ [0.01, 0.3] MPa, T_i ∈ [1, 100] d |
| T5 | Two wells: does allocation matter more than shape? | Compare the gain with the single-well shape gain (~5%) |
| T6 | Event count vs peak rate | Pareto front via ε-constraint |

## Poroelastic kernel (exact for an irrotational infinite medium)
For a radially symmetric p(r,t), u = ∇φ with ∇²φ = κp, where κ = α(1−2ν)/(2G(1−ν)). This gives:
- σ_ij = 2G(φ_,ij − δ_ij∇²φ);
- σ_rr = −4η M/r³ and σ_θθ = 2η(M/r³ − p), with η = α(1−2ν)/(2(1−ν)) and M = ∫₀ʳ s²p ds.

The Coulomb change on a plane with normal n_f and slip s is
ΔCFS = (σ_rr−σ_θθ)(n·n_f)(n·s) + f[σ_θθ + (σ_rr−σ_θθ)(n·n_f)² + p].

All terms are linear in p and therefore in q, so Theorem C1 still applies.

## Results
**T1: convexity holds and the optimum is unique.** With poroelastic stress, the optimum gains −0.9% (stabilising geometries) and −4.7% (destabilising 45°) against constant rate. The "taper-last-30%" schedule costs **+1.7 to +4.1%** more events in every geometry. It lowers post-shut-in events (15.7 → 12.9 d) but raises co-injection events by more.
→ *Taper moves events in time; it does not remove them.* We could not reproduce a count benefit from tapering in this model. Segall & Lu's benefit was for a stabilising-poroelastic fault with an abrupt shut-in jump; our destabilising case shows peak rate at shut-in but still no net benefit.

**T2: peak rate.** In the t ≪ t_a limit, R ≈ r·e^{ΔS/a}, so min max R ⇔ min max ΔS, which is an **LP**.
- Exact smoothed optimisation, 4 multistarts: R/r = 1.7338–1.7343 (consistent).
- LP surrogate: **1.7025**, slightly better than the smoothed exact attempt.
- Constant rate: 2.058. The LP gives −17% peak rate.

**T3:**
- **(a) Fixed net volume (storage):** the convex optimum uses **zero flowback** for both caps tested (10% and 25% of q_max).
- **(b) Fixed injection, then a flowback budget B:**

| B | 5% | 10% | 20% | 40% |
|---|---|---|---|---|
| Event reduction | −6.0% | −11.0% | −19.1% | −31.7% |

The optimal flowback schedule = **immediate flowback at maximum rate** (identical to 3 significant figures). This gives an optimal-control derivation of the "immediate extraction" recommendation of the Pohang study (EPSL 2025, 10.1016/j.epsl.2025.119579).

**T4: robustness of the window artefact** (150 random sites):

| Statement | Share of sites |
|---|---|
| Step-up beats constant when counting during injection only | 61% |
| Step-up beats constant when the post-shut-in tail is counted | 36% |
| Ranking flips from "better" to "worse" | **26%** |

→ The artefact is common but not universal. Any comparison of schedules must report the counting window.

**T5: two wells, two fault patches.** The global optimum gives **−39%** against an equal constant split, by moving 97% of the volume to the well that is far from the critical patch. Allocation (*where*) is worth roughly 8× the schedule shape (*when*). This agrees with the qualitative conclusion of Sci Rep 2026 (10.1038/s41598-026-74920-7), and here it comes with a certified global optimum.

**T6: Pareto front (count vs peak).**

| Point on the front | Peak R/r | Events |
|---|---|---|
| Count-optimal end | 4.55 | 37.01 |
| Peak-optimal end | 3.40 (**−25%**) | 38.18 (**+3.2%**) |

The front is convex (both objectives are convex in the t ≪ t_a limit). → A large peak-hazard reduction costs almost nothing in count.

## Rethink
- **Strong:** convexity (theorem plus checks), the LP for peak rate, and the immediate-flowback optimality are certified within the model class.
- **Weak:** absolute gains depend on Π. The model has no magnitudes and no permeability change. Real sites need the engine (H1c), because heterogeneity and multiphase flow are captured only by its Green's functions and are only approximately linear.
- **Next:**
  - H1c: put the method in the engine, then forward-verify the optimum in the nonlinear simulator.
  - H1d: magnitude-aware objective (Mmax via Gutenberg–Richter: P(M ≥ m) = 1 − exp(−N·10^{−b(m−Mc)}), which is monotone in N, so it is still convex-equivalent for one patch).
  - H1e: replicate a field case (Otaniemi 2018 public injection data).

## Addendum: T7, strongly coupled poroelastic regime (exp/h1b_regime.py)
Settings: α = 1, ν = 0.15, η = 0.41. This test was prompted by ARMA 2026 (10.56952/arma-2026-0542), found by the literature watch.

| Geometry (early-time kernel sign) | Taper effect on count | Count-optimum, last 20% of rate | Peak-optimum, last 20% of rate | Peak reduction |
|---|---|---|---|---|
| Stabilising, nf∥n (−, later +) | worse | 1.14× (end-loaded) | **0.39× (taper)** | −9.3% |
| Stabilising, 60° (− throughout, stress shadow) | worse | 1.18× | **0.43× (taper)** | −5.8% |
| Pressure-dominated (+) | worse | 1.27× | 0.74× (front-loaded) | −23.7% |

**F7 revised:** taper is objective- and regime-dependent. The convex/LP framework chooses the right one automatically, with no rule imposed. This is consistent with the ARMA 2026 regime criterion, which we reach here as a consequence of optimality.
