# Research Cycle 11 plan and pre-registration (C49, 2026-11-03)

Written before any Cycle 11 run. Engine: ResSim C49 (adds FAULTMECH item 10, the ramp width of the two-way permeability
factor; default 0.05 reproduces C48 exactly). Scripts and logs go to `exp/`. Post-hoc runs will be labelled exploratory.

Cycle 10 left one main open question (F31): convexity of the expected event count is destroyed when a fault face changes
its hydraulic state on reactivation. Cycle 11 asks whether that can be repaired, and checks the robustness of F33 (H11).

| Test | Question |
|---|---|
| H10f-a | Does a 10× wider (smoother) permeability ramp restore convexity? |
| H10f-b | Does a reactivation constraint (design only where no face reaches CFF = 0) give a convex feasible set on which the objective is convex? |
| H11b | Is the GP-EI result of H11 robust to the choice of the training set? |

## H10f-a: smoothed permeability ramp

**Model.** H10e cases NS (CO₂) and WNS (water), unchanged: sealing fault MULTFLT 0.02, permFactor 50, A1 pre-load
(32 / 27 bar). The only change is FAULTMECH item 10 = 0.5: the factor is permFactor^x, x = clamp(CFF / (0.5 σ'n), 0, 1),
instead of 0.05 σ'n (H10e). The factor stays irreversible (maximum over time) and explicit (one report-step lag).

**Pre-check.** At u = 1 the transmissibility multiplier of at least one face must exceed 1 (reactivation happens). If not,
the case is non-discriminating and is reported as such (no re-tuning of the pre-load).

**Test.** M = 40 pairs, seed 13 (the H10e pairs), ε from the 10 d / 5 d report-step difference, as H10e.

**Criteria (as H10 / H10e).** ≤ 5 % violations → SUPPORTED (convexity restored). > 15 % violations with median violator
δ > 1 % → REFUTED. Otherwise INCONCLUSIVE. Each case is judged separately.

**Stated expectation (both directions).** Fewer violations than H10e (9/40 NS, 19/40 WNS), because the response to CFF is
less abrupt. But the switch remains irreversible and history dependent, so I expect it to stay above 5 % (INCONCLUSIVE or
REFUTED). If it reaches ≤ 5 %, smoothing alone is a usable repair; if violations do not fall at all, the ramp width is not
the cause and the irreversibility is.

## H10f-b: reactivation constraint

**Model.** NS and WNS with permFactor 50, ramp 0.05 (H10e switching), but a pre-load 3 bar lower than A1 (29 / 24 bar), so
that the u = 1 design stays below failure. A design is *feasible* when no fault face reaches CFF ≥ 0 at any report time
(from `--fault-csv`, column `cff`). A feasible design therefore never switches the permeability.

**Pre-check.** At u = 1, max CFF < 0. At least 10 of the 40 pairs must have both A and B feasible; otherwise the case is
non-discriminating.

**Test.** The same 40 pairs (seed 13). For pairs with A and B feasible:
1. *Set convexity:* is C = (A + B)/2 feasible?
2. *Objective convexity on the feasible set:* for pairs with A, B and C feasible, a violation as in H10e.

**Criteria.** Set convexity: ≥ 95 % of C feasible → SUPPORTED; < 85 % → REFUTED; otherwise INCONCLUSIVE. Objective:
≤ 5 % violations → SUPPORTED; > 15 % with median violator δ > 1 % → REFUTED; otherwise INCONCLUSIVE. The hypothesis
H10f-b holds when both parts are SUPPORTED.

**Stated expectation.** SUPPORTED. With no switching, max CFF over faces and times is a maximum of nearly linear
functions of u (F26: small, favourable nonlinearity), hence nearly convex, so its sublevel set is nearly convex. The
objective is then the H10e no-switch control (0/40). A failure of the set test would mean CFF is noticeably concave in u.
**Design cost (descriptive, not a criterion).** The fraction of the 80 sampled A/B designs that are feasible, and the
lowest feasible N against the lowest N over all sampled designs.

## H11b: GP-EI robustness to the training set

**Model.** As H11 (h9_rsf_surrogate.cpp, mode m0gp, 12 training + 12 EI = 24 RSF evaluations; objective: ensemble-mean
seismic moment). The training designs are drawn with four new training seeds (environment `H11_TRAIN_SEED` = 101, 102,
103, 104; H11 used the fixed generator seed 9). The EI candidate sampler is unchanged.

**Criteria.** Ratio = best GP-EI moment / ES median best of H9b (2.21503·10¹⁵ N m). Median ratio over the four seeds
≤ 1.05 and every ratio ≤ 1.10 → SUPPORTED. Median > 1.15 → REFUTED. Otherwise INCONCLUSIVE.

**Stated expectation.** SUPPORTED (H11: 1.002; EI-seed robustness 1.001–1.006). A training set that misses the low-moment
region could, however, leave EI exploring with too few evaluations, which would show as ratios above 1.05.

## Order and resources

H11b first (≈ 3 min per seed), then H10f-a and H10f-b (≈ 120 runs per case each). At most two experiment processes at a
time, OMP_NUM_THREADS = 1, never together with the full test suite or the MinGW build.
