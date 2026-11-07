# Research Cycle 9 report (C47, 2026-11-02)

The plan and the pre-registrations are in `CYCLE9_PLAN.md`; the findings are F25–F29 in `FINDINGS_LOG.md`. Every run
uses the ResSim C47 engine. Scripts and logs are in `exp/`.

**The new direction in this cycle:** move the P1 convexity theory from analytic pressure into a **multiphase reservoir
simulator coupled to a fault**, and ask three questions:
1. Does convexity survive multiphase flow?
2. *Why* does it survive?
3. Can it be used for design?

| Test | Status | Result |
|---|---|---|
| H10 CO₂ / water convexity (4 cases × 40 pairs) | pre-registered | **SUPPORTED** 0/160 violations |
| H10c plume across the fault, near fault (4 × 40) | exploratory | **SUPPORTED** 0/160 |
| H10d exact gap decomposition (5 × 10) | exploratory, expectation stated | sign expectation **REFUTED**; abs(NL) ≪ abs(CV) |
| H10b simulator-driven SCP (N, F) | pre-registered | **SUPPORTED** (≤ 1 % at iteration 2, 11 runs) |
| H9 Dieterich surrogate vs ES, event count | pre-registered | SUPPORTED but non-discriminating |
| H9b same, seismic moment | pre-registered | **REFUTED** (ratio 1.83) |
| H2d gas fracture growth at 300 sm³/d | pre-registered | **INCONCLUSIVE** (fracture stalls at 35 m) |
| H2d′ same at 600 sm³/d | exploratory addendum | t60 SUPPORTED (0.60 vs 1.25 d); t90 not (2.19 vs 2.22 d) — F30 |

## H10 / H10c: convexity under multiphase CO₂ injection (F25)

**Setup.**
* Deck: CO2STORE section of 40×1×10 cells, dx 10 m. The injector is at column 1, layers 8–10.
* Fault: FAULTMECH fault (friction 0.6) at column 15 (N), 30 (F) or 6 (NN).
* Schedule: four 90-d periods with Σu = 4 and 0.25 ≤ u ≤ 1.75, then 360 d shut-in.
* Objective: the Dieterich expected event count of the fault.
* Water controls (W, WF, WNN, WNH) inject 150 sm³/d of water.

**Results.**
* There are no midpoint violations in any of the 8 cases (320 pairs).
* The weakest margin is case F: median δ = −0.93 % and worst −0.037 %.
* In NH, CO₂ crosses the fault by day 180 (11.5 expected events), and convexity still holds (median δ −5.4 %).

**Plume caveat.** In N and F the fault feels the *pressure* front, not the CO₂ itself:
* N: Sg ≤ 0.04 at column 15, mostly after shut-in;
* F: the plume stays within 150 m of a fault at 300 m.

NH and NN were therefore added (H10c) to put the CO₂ at the fault.

## H10d: why convexity survives (F26, new)

The "nonlinearity" metric (history norm of the dCFS midpoint residual) was 0.66 % in F, about the size of the convexity
margin, yet F had no violations. H10d decomposes the gap exactly:

  δ = R(D_C) − ½[R(D_A) + R(D_B)] = NL + CV,
  CV = R(D̄) − ½[R(D_A) + R(D_B)] ≤ 0 (Jensen; R is convex in the history),
  NL = R(D_C) − R(D̄).

R is the Dieterich functional, replicated from `--fault-csv` to 2·10⁻⁷.

| case | median NL/abs(CV) | min / max | NL > 0 | median s |
|---|---|---|---|---|
| N | −0.018 | −0.94 / +0.052 | 3/10 | −2.7·10⁻⁴ |
| F | −0.058 | −2.21 / +0.051 | 2/10 | −3.0·10⁻⁴ |
| NN | −0.014 | −0.73 / +0.064 | 4/10 | −4.1·10⁻⁴ |
| NH | −0.010 | −0.081 / +0.181 | 4/10 | −2.6·10⁻⁴ |
| W | +0.000 | 0 / 0 | 10/10 | +4.7·10⁻⁶ |

**What this shows.**
* The expectation stated beforehand was a concave map (sub-linear pressure build-up), so NL > 0. That expectation is
  **refuted**: the CO₂ map leans slightly convex (s < 0), so NL mostly *helps*.
* Where NL hurts, it is at most 18 % of the Jensen margin. The largest value occurs with the plume across the fault (NH).
* The 0.66 % history-norm nonlinearity is a temporal *redistribution* of stress (s is only ~3·10⁻⁴). The event functional is
  almost insensitive to it.
* **Practical meaning:** a norm-based nonlinearity check overstates the risk to convexity. The event-weighted decomposition
  is the right diagnostic, and it costs no extra runs.

## H10b: design by sequential convex programming (F27)

* **Method.** Per iterate: 1 history run plus 4 pulse runs give a linearised Coulomb history. The Dieterich subproblem on
  it is convex and is solved with SLSQP. A trust region is checked against the true simulated count.
* **N.** Starting at [1.75, 1.75, 0.25, 0.25], SCP reaches 0.70729 at iteration 2 (11 runs). That is 1.25 % *below* the best
  of the 231-point exhaustive grid (constant schedule, 0.71624).
* **F.** −0.09 % against the grid.
* **Agreement.** The primary and constant starts converge to the same point to 2·10⁻⁷. That is the global optimum that
  convexity predicts.
* **Shape.** The optimum front-loads slightly ([1.12, 0.99, 0.89, 1.00] for N), but constant is within 1.3 %. This agrees
  with P1: in this regime (θ small) the schedule matters little.

## H9 / H9b: surrogate design on the rate-and-state fault (F28)

* **H9 (event count).** The count is quantised: system-wide events of the 16-member ensemble. Constant, surrogate and ES all
  reach 32 events, so the test is non-discriminating.
* **H9b (continuous objective: ensemble-mean seismic moment).** The surrogate fits the 12 training schedules with R² 0.85.
  Its optimum [0.89, 0, 0.61, 2.50], however, gives 4.05·10¹⁵ N·m against the ES median best of 2.215·10¹⁵ (ratio 1.83).
  H9b is therefore **REFUTED**.
* **Constant schedule.** It is 6.4 % above the ES best, so this test does discriminate.
* **Reading.** This is consistent with F20 (moment is not convex). A Dieterich surrogate is *not* a valid design proxy for
  the moment of an RSF fault. Moment-aware design needs the RSF model in the loop, or a moment surrogate.

## H2d: gas fracture growth (leftover from F24)

* At 300 sm³/d the gas-driven fracture stalls at 35 m (K_I ≈ K_IC by day 84), so t60 is undefined on every grid.
* The pre-registered pre-check only covered growth that is *too fast*. It fell back to 100 sm³/d, where there is no growth
  either. **H2d is INCONCLUSIVE (design failure).**
* The script's automatic verdict ("REFUTED", from NaN comparisons) was corrected in the log, and the script was fixed.
* H2d′ (exploratory, declared before running): 600 sm³/d, where t60 = 18.8 d at 10 m.

| dx (m) | SG t60 | SGL t60 | SG t90 | SGL t90 |
|---|---|---|---|---|
| 10 | 18.80 | 16.60 | 40.78 | 36.21 |
| 5 | 17.80 | 16.60 | 39.32 | 37.30 |
| 2.5 | 17.55 | 17.20 | 38.56 | 38.40 |
| spread | 1.25 | **0.60** | 2.22 | 2.19 |

* The t60 criterion is met (0.60 d < 1 d, and below SG).
* At t90, however, the SGL spread is no smaller than the SG spread, and the two converge toward ≈ 38 d from opposite
  sides. In gas, SGL over-corrects on coarse grids (5 % fast at 10 m), which water (F21) and corner point (F24) did not
  show (**F30**).
* A likely cause, still to be tested: the lagged sub-cell storage keeps the leak-off reduction while compressible gas fills
  the fracture. A gas-compressibility-aware sub-cell storage is the next step (H2e).

## Engineering found by the research track

* **Low-pressure gas-table padding** (OPM PR 3779) and the **reference-drift study** (F29) are documented in
  ENGINES_AND_MATH §9.t.
* **`--fault-csv`**, the per-face fault history, was added for H10b / H10d and is now a user feature.

## Literature watch

* **Citation gap closed: Candela et al. 2022** (Front. Earth Sci. 10:1053951). It optimises injection under a Dieterich-rate
  constraint for Oklahoma, with no statement about the structure of the problem (convexity). GAP_MAP had listed it since
  Cycle 1 ("stochastic gradients"), but P1 never cited it by name. It is now cited in the manuscript, the GRL tex and the
  JGR long version.
* All eight "[TITLE TO VERIFY]" references of P1 were resolved through Crossref / arXiv.
* Nothing was found on EDFM sub-cell leak-off (F21 / F24) or on a decomposition like H10d.

## Next cycle (proposal)

1. **P1 addendum.** Add a short section with F25–F27: convexity in a coupled multiphase simulator, the decomposition and
   SCP. The GRL letter stays as it is; this goes into the JGR long version.
2. **H10e.** Two-way coupling (pore pressure ↔ fault permeability with FAULTMECH slip-dependent transmissibility). This is
   where NL could turn harmful.
3. **H2e.** Sub-cell storage with gas compressibility (F30) on the H2d′ deck.
4. **H11.** Moment-aware design on the RSF fault (F28): a GP surrogate on M0 against ES, with the same budgets.
