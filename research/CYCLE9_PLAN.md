# Research Cycle 9 plan (C47, 2026-11-02), pre-registered before any run

Cycle 9 has two leftover hypotheses from Cycle 8 and one new research direction.

**Why a new direction.** P1 (F1–F19) proved convexity for Dieterich seismicity when the Coulomb stress is *linear* in
the injection rates. P2 (F20, F22) showed that a rate-and-state (RSF) fault breaks convexity. The industrially important
case sits between the two and was never tested: **CO2 storage**. There the seismicity law is still Dieterich, but the map
from rate to fault pressure is nonlinear. Three effects make it nonlinear:
* gas compressibility;
* relative permeability, since the CO2 plume is more mobile than brine;
* dissolution.

## Literature scan (OpenAlex, 10 queries, 2022–2026; `lit/q47.py`, `lit/q47.json`)

| Paper | What it does | What it does not do |
|---|---|---|
| Hill et al., JGR 2024 (10.1029/2024jb028759) | LP on the seismogenic index, Raton Basin, single phase | No multiphase flow, no convexity statement |
| Zheng, Jha & Jafarpour, SPE 2023 (10.2118/215155-ms) | Multiphase flow + geomechanics, CO2 rate optimisation (caprock fracturing objective, ensemble) | No seismicity-rate objective, no convexity |
| Cui et al., Sci Rep 2026 (10.1038/s41598-026-74920-7) | CO2–brine, 2-D: rate history reorganises the fault-yield stress paths. Calls for joint optimisation of schedules | No optimisation, no Dieterich |
| Acosta et al., Flow2Quake, IJGGC 2025 (10.1016/j.ijggc.2025.104388) | Multiphase VFE flow → stress → seismicity forecasting | Forecasting only, no schedule design |
| Walter, Boyet, Stefanou & Vilarrasa, EarthArXiv 2026 (10.31223/x5p75x) | Review: real-time control | No convexity / CO2 rate-map analysis |
| Gutiérrez-Oribio & Stefanou, GETE 2024 | Tracking control, single-phase diffusion | — |

**Not found:** any test of whether the expected induced-event count stays convex in the CO2 injection-rate history when
pressure comes from multiphase flow. Also not found: a sequential convex scheme driven by the simulator for this
objective.

---

## H10 (new): convexity of the Dieterich count under multiphase CO2 injection

**Question.** The simulated quantity is N(u): the expected Dieterich events on a fault over a fixed window, with the stress
from the full nonlinear ResSim CO2STORE run. Is N(u) midpoint-convex in the rate history u?

**Design (`exp/h10_co2_convexity.py`, decks from `exp/h10_mkdeck.py`).**
* **Model.** 2-D vertical section, 40 × 1 × 10 cells of 10 m × 100 m × 5 m. The aquifer is at 1000 m, k = 50 mD and
  φ = 0.25. A large pore volume at the far boundary gives open-aquifer support.
* **Fault.** Fault F1 (X faces) at x = 150 m (case **N**, near: the plume reaches the fault) or x = 300 m (case **F**, far:
  pressure only).
* **Fault mechanics.** One-way poroelastic STRESSFIELD (no FE, no two-way coupling, no fracture growth), so the
  forward map is smooth. FAULTMECH: μ 0.6, Aσ 0.1 MPa, t_a 100 y.
* **Injection.** Injector at the bottom of column 1. K = 4 intervals of 90 d, then 360 d of shut-in. The window
  (720 d) includes the post-shut-in tail, so the F2 window artefact is excluded. Mean rate q̄ (CO2: 40 000 sm³/d;
  surface volume fixed), u_k ∈ [0.25, 1.75] q̄ with mean(u) = 1.
* **Control (case W).** The same grid and fault in a black-oil deck with water injection at the same reservoir
  volume rate. It is single phase, so the pressure is nearly linear in u and N should be convex (F1).
* **Test.** M = 40 random pairs (A, B) per case, midpoint C = (A + B)/2. A violation is
  N(C) > (N(A) + N(B))/2 + ε. The numerical noise ε is measured *before* the pairs: N(constant u) under
  time-step perturbations (TSTEP halved), ε = 3 × the observed spread. The relative violation is
  δ = (N(C) − (N(A) + N(B))/2) / N(C).

**Pre-registered verdicts (per case N, F).**
* **SUPPORTED** (convexity survives multiphase flow): ≤ 5 % of the pairs violate.
* **REFUTED:** > 15 % of the pairs violate **and** the median δ of the violators is > 1 %.
* Otherwise **INCONCLUSIVE**.
* **Validity check:** the control W must show ≤ 5 % violations. If it does not, the test is invalid and is
  redesigned before any interpretation.
* **Secondary (descriptive, not a verdict):** the nonlinearity index of u → mean fault ΔP, i.e.
  |ΔP(C) − (ΔP(A) + ΔP(B))/2| / ΔP(C), reported per case.

**What each result would mean.**
* SUPPORTED: the P1 convex programme carries over to CCS schedule design (strong practical claim).
* REFUTED: CCS needs a sequential approach, which motivates H10b.

## H10c (exploratory, added after the H10 verdicts and before any H10c run)

**Why it was added.**
* **The plume barely reached the fault in case N.** The plume check (`exp/h10_plume.cpp`) shows the front reaching the
  fault-adjacent column only after shut-in, with Sg ≤ 0.04.
* **The descriptive "max dCFS nonlinearity" metric of H10 is flawed.** max over time is itself nonlinear in u, so the
  linear water control also scores about 3 %. It is replaced by the history-norm index
  ‖dCFS(C) − (dCFS(A) + dCFS(B))/2‖₂ / ‖dCFS(C)‖₂ over faces × report steps (`exp/h10_nonlin.py`).

**Strong two-phase cases.**
* **NN:** fault at x = 60 m (the plume crosses it during injection).
* **NH:** fault at x = 150 m at 3× the rate (120 000 sm³/d).
* **Controls:** water controls WNN and WNH, same geometry and reservoir volume rate.

**Rules.** Same M = 40 pairs (seed 11), ε rule and verdict thresholds as H10. Labelled exploratory: it cannot change the
H10 verdict, it only maps where the result holds.

## H10b (conditional on H10 running): sequential convex programming driven by the simulator

**Method (SCP).**
* At iterate u_j, linearise the fault pressure response by finite-difference unit pulses from the simulator
  (K + 1 runs).
* Solve the convex Dieterich subproblem in a trust region ‖u − u_j‖∞ ≤ Δ.
* Accept or reject on the true simulated N, and update Δ.

**Reference.** The best N over a 4-D grid search (step 0.25, about 300 feasible points, run once).

**Pre-registered.**
* **SUPPORTED:** SCP reaches ≤ 1 % above the grid-search best within ≤ 6 iterations (≤ 30 simulator runs) in both cases
  N and F.
* **REFUTED:** > 5 % above the best after 6 iterations in either case.
* **Start points (addendum, fixed before any SCP run).** Primary start: front-loaded u0 = (1.75, 1.75, 0.25, 0.25).
  Secondary start: constant u0 = 1. Pulse step h = 0.1, initial trust radius Δ = 0.5. A step is accepted if
  ρ = actual / predicted decrease > 0.1. Δ doubles if ρ > 0.75 (cap 0.75) and halves if ρ < 0.25.
  Each iteration costs 4 pulse runs + 1 evaluation run. The verdict uses the primary start.
* **Grid caveat.** If the grid search is too expensive (> 45 min), the reference becomes the best of 200 random
  feasible schedules. This substitution is reported as such.

## H9 (leftover P2): a Dieterich surrogate calibrated to the RSF fault versus stochastic search

**Setting.** Same fault, forcing and heterogeneity ensemble (E = 16, common random numbers) as H7b/F22, so the
objective is the ensemble-mean RSF count.

**Surrogate.**
* Dieterich at the fault-centre pressure. Fit Aσ, t_a and r0 by Poisson least squares to the RSF counts of
  **T = 12 training schedules** (random, mean 1).
* Optimise the surrogate globally (convex, F1; SLSQP from 5 starts).
* Evaluate the optimum in RSF: **12 + 1 = 13 RSF ensemble evaluations**.

**Competitor.** Random search plus (1+1)-ES with a **budget of 60 RSF ensemble evaluations**, 3 seeds; the result is the
median over seeds. The constant schedule u = 1 is reported as a baseline.

**Pre-registered.**
* **SUPPORTED:** the surrogate optimum's RSF count is ≤ 1.05 × the stochastic median best, at ≤ 1/4 of its RSF
  evaluations.
* **REFUTED:** > 1.15 × the stochastic median best.
* Otherwise **INCONCLUSIVE**.

## H9b (added after the H9 verdict, before any H9b run)

**Why it was added.** H9 met its criterion (ratio 1.000), but the test turned out to be non-discriminating:
* the ensemble count is quantized (system-wide events, multiples of 16);
* the constant schedule already reaches the floor of 32;
* the surrogate, the ES and the constant schedule all score 32.

**Change.** H9b repeats H9 with the continuous objective **ensemble-mean seismic moment M0**. M0 broke convexity in 32 %
of the H7b pairs. The surrogate is fitted to M0 by least squares on log M0 (scale + Aσ + τ̇).

**Rules.** Same budgets (12 + 1 vs 60 × 3 seeds) and the same thresholds (≤ 1.05 SUPPORTED, > 1.15 REFUTED). The constant
schedule is reported. If the constant schedule is within 1 % of the best found by any method, H9b is reported as
non-discriminating as well.

## H10d (exploratory, added after the H10 / H10c / history-norm results, before any H10d run)

**Why.** In case F the history-norm nonlinearity of dCFS (median 0.66 %) is as large as the convexity margin (median
δ = −0.72 %), yet no pair violated convexity. A norm does not show the *sign* of the nonlinearity or how much it costs in
events. **Question:** does multiphase flow push the midpoint toward or away from violation, and by how much?

**Method (exact decomposition).** Let R be the Dieterich functional of the per-face Coulomb history. It is replicated from
`--fault-csv` to 1e-7 (h10b_scp.py). For a pair A, B with midpoint C, let D̄ = (D_A + D_B)/2. Then

  δ = R(D_C) − ½[R(D_A) + R(D_B)] = **NL** + **CV**,  NL = R(D_C) − R(D̄),  CV = R(D̄) − ½[R(D_A) + R(D_B)] ≤ 0

CV ≤ 0 holds by Jensen, since R is convex in the history (F1). NL is the price of the rate → stress map being nonlinear.
We also report the signed history residual s = Σ(D_C − D̄) / Σ|D_C| (over faces and steps). s > 0 means the map is
concave, i.e. sub-linear pressure build-up.

**Setup.** Cases N, F, NN, NH, W; M = 10 pairs each. Seed 9, same sampler as H10.

**Pre-stated expectation (not a pass / fail test).**
* CO2 cases: NL > 0 (concave map, from the mobility gain of the CO2 plume) in most pairs, and |NL| < |CV|.
* Water case W: |NL| ≈ 0.

We report the median of NL/|CV|. A case with median NL/|CV| > 0.5 is flagged as **near the edge**: a stronger
nonlinearity could break convexity there.

## H2d (leftover F24 gas variant): gas growth slow enough to discriminate

**Change from H2c GAS.** Rate 300 sm³/d instead of 3000 (the growth now takes weeks). The criterion is unchanged from
H2c: SGL t60 spread < 1 d over dx 10 / 5 / 2.5 m and smaller than SG → SUPPORTED. A pre-check is added: if SG t60 at
10 m is < 10 d, the variant is still non-discriminating and the rate is lowered to 100 sm³/d once.

**H2d outcome and addendum H2d′ (written after the H2d runs, before any H2d′ run).** At 300 sm³/d the gas-driven fracture
stalls at 35 m (K_I ≈ K_IC by day 84), so t60 is undefined on every grid. The pre-check only covered growth that is too
fast; it wrongly fell back to 100 sm³/d, where there is no growth either. **H2d is INCONCLUSIVE (design failure).**
H2d′ is exploratory, not a confirmatory test: rate 600 sm³/d, where t60 at 10 m (SG) is 18.8 d. It uses the same three
grids, the same SG / SGL pair and the same 1-d spread criterion.

## Also

* Scholar / web check of F21 / F24 (sub-cell storage in EDFM growth): pending since Cycle 8.
* Literature watch run (`lit/watch.py`).
