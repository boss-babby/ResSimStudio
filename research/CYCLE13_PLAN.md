# Research Cycle 13 plan: seismicity-aware cyclic storage (NEW_FIELDS N1), pre-registered 2026-11-04 (C50)

Written **before** any pair run or design optimisation. The only runs made before writing this were u-base smoke tests
(constant schedule, listed under "Pre-checks"); they carry no pair information. The text below is not edited after the runs;
results are appended at the end.

## Question
In cyclic storage, q(t) is signed (inject > 0, withdraw < 0). Under linear coupling the Dieterich expected count is convex in the
signed schedule (P1 theorem; NEW_FIELDS N1). (a) Does this survive in the simulator, including the strongly nonlinear pressure
response of a depleted gas field? (b) Do Dieterich and Kaiser-memory objectives give different optimal storage cycles?

## Model (exp/h12_mkdeck.py)
2-D section 40×10 cells (as H10); a large-pore-volume boundary column; fault F1 at x = 150 m, not sealing; FAULTMECH one-way,
μ 0.6, cohesion 0, Aσ 1 bar, t_a 100 yr, r0 5, b 1; STRESSFIELD 0.226 0.7 0.95 30 0.5. Schedule: 8 intervals of 45 d
(1 to 4 injection, 5 to 8 withdrawal), then 180 d rest; report step 5 d. Separate injector INJ and producer PRD in column 1.
* **SW**: water storage, initially water-filled; Q_W = 150 sm³/d.
* **SG**: depleted dry-gas field (Sw 0.2, Sg 0.8, p 110 bar), PVDG with Bg 0.025 to 0.003 rm³/sm³.

## Pre-checks (u-base = [1,1,1,1,−1,−1,−1,−1], done before writing)
* SW: max ΔCFS 3.98 bar, N = 0.2535; injected = withdrawn = 27,000 sm³.
* SG at 12,000 sm³/d: max ΔCFS 0.136 bar. That is near-linear at Aσ = 1 bar, so a convexity test would be non-discriminating.
  Max ΔCFS is nonlinear in rate: 0.94 / 2.66 / 4.33 bar at 120k / 360k / 600k sm³/d.
* **Rule G (calibration, base schedule only):** Q_G is set by bisection on [360k, 600k] so that max ΔCFS(u-base) = 3.98 bar ± 5 %.
  This gives a realistic large UGS pressure swing.

## H12a: convexity of signed storage schedules in the simulator
* 40 random pairs (seed 31). Injection multipliers u1..4 drawn as in H10 (normalised to mean 1, each in [0.25, 1.75]);
  withdrawal u5..8 = −(an independent draw of the same kind). Every schedule is inventory-neutral; midpoints are feasible.
* Tolerance ε = 3|N(dt 5) − N(dt 2.5)|/N(dt 5) at u-base, per case.
* Verdict as H10: **SUPPORTED** if ≤ 2/40 violate; **REFUTED** if > 6/40 violate and median violator δ > 1 %; inconclusive otherwise.
* **Predictions:** SW SUPPORTED (near-linear control). SG SUPPORTED: F26 found that rate-to-stress nonlinearity helped convexity for
  CO₂. A REFUTED SG would mean gas PVT nonlinearity under withdrawal can break the convex design. That is a publishable
  negative, and it is what makes the test able to fail.

## H12b: Dieterich vs Kaiser-memory optimal cycles (SW only; response matrix)
* Response matrix: 8 engine runs with u = e_k (unit pulse in interval k) give per-face ΔCFS_f(t) = Σ_k u_k G_fk(t). Linearity
  check: the Python Dieterich count (same Simpson / sub-step scheme as faultmech.cpp) evaluated on the matrix prediction must
  match the engine count for 5 random H12a schedules within 2 % (median). **If not, H12b is declared NON-DISCRIMINATING and not
  interpreted.**
* Feasible set: u1..4 ∈ [0, 2] with mean 1; u5..8 ∈ [−2, 0] with mean −1 (same demand and inventory as u-base).
* Objectives:
  * N_D = Dieterich expected count (Aσ as stated);
  * N_K = Σ_f max(0, max_t ΔCFS_f(t)) (Kaiser: only new stress maxima relative to the storage start count, linear productivity).
* Optima: u_D (SLSQP, 8 starts, best kept) and u_K (LP epigraph form, scipy linprog).
* Metrics: L1 = Σ|u_D − u_K| / Σ|u_K|; cross-costs X_D = N_D(u_K)/N_D(u_D) − 1 and X_K = N_K(u_D)/N_K(u_K) − 1.
* Two Dieterich levels (one response matrix, no extra engine runs): **Aσ = 1 bar** (Π = max ΔCFS/Aσ ≈ 4) and **Aσ = 10 bar** (Π ≈ 0.4).
* **Predictions (from theory, stated in both directions):**
  * **H12b-1 (Π ≈ 4):** large Π makes Dieterich a soft maximum, so the two optima nearly coincide: X_D < 5 % and X_K < 5 %.
  * **H12b-2 (Π ≈ 0.4):** Dieterich becomes near-linear in the time-integrated stress, while Kaiser still sees only the peak, so
    the optima differ: L1 > 10 % and max(X_D, X_K) > 5 %.
  * If H12b-1 fails, the memory model matters even at large Π. If H12b-2 fails, it does not matter in either regime. Both outcomes
    are reportable.
* Engine verification: u_D and u_K (Aσ 1) are rerun in the engine. The Dieterich ordering N_D(u_D) ≤ N_D(u_K) must hold in the
  engine; this is a check of the matrix design, not a hypothesis.

## Honest limits (stated in advance)
2-D, one passive fault, one cycle, one-way coupling. The Kaiser objective is a simple linear-productivity form. The depleted-field
history (stress before storage) is ignored: memory is relative to the storage start. For a depleted field the pre-production
maximum would make Kaiser predict little seismicity at all, which is a separate question.

## Amendment A1 (written after the first H12b run, before any interpretation; H12a unaffected)
The pre-registered linearity check failed (median 8.2 %, NON-DISCRIMINATING; first run kept in h12.log). Diagnosis:
* the Python Dieterich copy reproduces the engine exactly on the engine's own ΔCFS history (0.25345858 vs 0.25345858);
* a zero-rate run has max |ΔCFS| 0.68 bar, from gravity re-equilibration of the uniform initial pressure. The response is affine.

The matrix ΔCFS = Σ u_k G_k counts this drift Σu_k times (0 for balanced schedules) instead of once. Correction: D(u) = G₀ + Σ_k u_k (G_k − G₀),
with G₀ the zero-rate run (convexity is unaffected by an affine offset). Everything else is unchanged: the 2 % check, objectives,
bounds, both Aσ levels, predictions and verdict rules. The drift is a deck artefact (non-equilibrated initial state). It acts
identically on both objectives and is listed as a limitation; the deck is not changed, so that H12a stays valid.

---

## RESULTS (appended after the runs, 2026-11-04)

### Calibration (rule G)
Q_G = 547,500 sm³/d gives max ΔCFS 3.974 bar (target 3.978).

### H12a: **SUPPORTED for both cases**
| Case | ε | Violations / 40 | Median δ | Max δ | Verdict |
|---|---|---|---|---|---|
| SW water storage | 0.11 % | 0 | −17.4 % | −0.27 % | SUPPORTED |
| SG depleted gas | 0.56 % | 0 | −24.5 % | −0.29 % | SUPPORTED |
Files: exp/h12_storage.py, h12a.log, h12_results.json.

### H12b (with amendment A1): **H12b-1 SUPPORTED, H12b-2 SUPPORTED**
* The first run (pure-linear matrix) failed the linearity check (8.2 %) and was not interpreted (exp/h12b_run1_invalid.log). With the
  affine matrix the check passed: median error 0.049 %.
* Aσ 1 bar (Π 3.98): X_D +4.0 %, X_K +4.8 % (both < 5 %); L1 0.59.
* Aσ 10 bar (Π 0.40): X_D +1.1 %, **X_K +40.1 %**, L1 0.68.
* Engine check: N(base) 0.25346, N(u_D) 0.25192, N(u_K) 0.26200; ordering OK.
* **Caveat found after the run:** the Kaiser LP optimum is not unique, because withdrawal timing barely changes new maxima. L1 is
  therefore not meaningful, and X_D depends on which optimum the LP returns. The H12b-2 verdict rests on X_K = 40 %, which is unique.

### Exploratory E13 (not pre-registered; exp/h12c_explore.py, h12c.log, h12c_results.json)
* Lexicographic Kaiser optimum (best Dieterich value within the Kaiser-optimal set): X_D 0.34 % (Π 4) and 0.15 % (Π 0.4).
* Robust min-max schedule: within 0.29 % (Π 4) and 0.15 % (Π 0.4) of **both** model optima.
* Constant schedule: +0.6 / +0.2 % under Dieterich, +8.0 % under Kaiser.
* Pattern: the models disagree only in the injection phase (Dieterich ramps up, Kaiser tapers). The optimal withdrawal is the same
  under both: front-loaded, declining.
* Practical rule (candidate): minimise new stress peaks first, then break ties with Dieterich. This costs < 0.35 % under either
  memory model, whereas a Dieterich-only design can cost 40 % if the fault has Kaiser memory, exactly in the low-Π regime where
  Dieterich says the design hardly matters.

Figure: papers/P5_storage/figs/fig_h12_storage.pdf (script papers/P5_storage/mkfig_h12.py).
