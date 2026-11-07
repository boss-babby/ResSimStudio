# Research Cycle 10 report (C48, 2026-11-03)

The plan and pre-registrations are in `CYCLE10_PLAN.md`; Amendment A1 was written before any pair run. The findings are
F31–F33 in `FINDINGS_LOG.md`, plus F34, an engineering finding from the OPM comparison. Every run uses the ResSim C48
engine. Scripts and logs are in `exp/`.

This cycle answers the three open questions that Cycle 9 left:
1. Does a *threshold* nonlinearity break convexity?
2. Can the dry-gas sub-grid leak-off over-correction be fixed?
3. Can a cheap surrogate design for seismic moment?

| Test | Status | Result |
|---|---|---|
| H10e convexity with fault-permeability switching (4 × 40 pairs) | pre-registered (+ A1) | **REFUTED** in 4/4 (9–21 of 40 violate) |
| H10e decomposition (4 × 10) | pre-registered | NL/abs(CV) median +0.73 … +0.89 (F26: −0.01 … −0.06) |
| NSC / WNSC controls without switching (2 × 40) | exploratory | 0/40 and 0/40 |
| H2e gas-aware sub-cell storage (SUBGRID LEAKC) | pre-registered | **SUPPORTED** (t90 spread 0.81 d vs 2.19 / 2.22) |
| H11 GP-EI surrogate for seismic moment | pre-registered | **SUPPORTED** (1.002 × ES median best in 24 evaluations) |
| H11 seed robustness (4 seeds) | exploratory | 1.001 – 1.006 × |

## H10e: convexity with two-way fault permeability (F31)

**Setup.** This is the H10 section (CO2STORE or water, 40×1×10 cells) with a sealing fault at column 15: MULTFLT 0.02,
and FAULTMECH permFactor 50, so a reactivated face reopens to T × 1.0. The four cases are:
* NS: CO₂;
* NHS: CO₂ at 3× rate;
* WNS: water;
* WNHS: water at 3× rate.

**Amendment A1.** With the default stress no face reactivated (max CFF −11 to −30 bar), so all four cases were
non-discriminating. A uniform shear pre-load was therefore chosen from the u = 1 run alone so that the most-stressed face
reaches CFF ≈ +2 bar: c = 32 / 23 / 27 / 14 bar. FAULTSLIP stays off, which is why entering the pre-load as cohesion is
valid. After A1, 3–8 faces reactivate at u = 1.

**Result.** Seed 13, M = 40:

| Case | Violations | Median δ (violators) | Max δ | Median NL/abs(CV) (10 pairs) | NL > 0 |
|---|---|---|---|---|---|
| NS | 9/40 (22.5 %) | +2.94e-2 | +8.59e-2 | +0.758 | 10/10 |
| NHS | 21/40 | +4.49e-2 | +1.28e-1 | +0.725 | 8/10 |
| WNS | 19/40 | +6.32e-2 | +2.59e-1 | +0.892 | 9/10 |
| WNHS | 16/40 | +4.94e-2 | +1.75e-1 | +0.806 | 8/10 |
| NSC (control, no switching) | 0/40 | — (median δ all −0.151) | −3.9e-3 | — | — |
| WNSC (control) | 0/40 | — (median δ all −0.177) | −5.4e-3 | — | — |

**Reading.**
* The switch alone breaks convexity, in water as clearly as in CO₂. The sealing fault and the pre-load alone do not.
* The nonlinearity term NL, which in F26 was small and mostly helpful, now reaches 70–90 % of the Jensen margin and is
  positive in most pairs. This happens because A, B and the midpoint C reactivate faces at different times.
* P1 claims convexity for a hydraulically passive fault. **A fault that changes its hydraulic state during injection is
  outside the guarantee.** In that case a convex programme (for example the SCP of F27) is a local method only and needs
  global checks such as multi-start or a reactivation-time constraint.

**Limits.**
* 2-D section, 4 periods.
* The pre-load was tuned per case. The shear-stress level controls how often reactivation occurs, so the violation
  fraction is not a universal number.

## H2e: gas-aware sub-grid leak-off storage (F32)

**Code reading first, as the plan required.** The sub-cell had storage, but it used only the pore-volume fraction. The
leaked gas was therefore stored as if it had the cell's mean compressibility, which over-states the leak-off reduction once
the fracture is full of compressible gas.

**Change.** The new option is SUBGRID LEAKC; the default is off.
* Sub-cell pressure: psub = p_res + (p_c − p_res)·(c_t,cell/c_t,sub)/φ.
* Sub-cell gas saturation: Sg,sub = Sg0 + ΔSg/φ.

**Result.** Times in days, at dx 10 / 5 / 2.5 m:

| Scheme | t60 | t90 | t90 spread |
|---|---|---|---|
| SGL-C | 17.41 / 16.92 / 17.00 | 38.27 / 37.47 / 38.25 | **0.81** |
| SGL (F30) | 16.60 / 16.60 / 17.20 | 36.21 / 37.30 / 38.40 | 2.19 |
| SG (F30) | 18.80 / 17.80 / 17.55 | 40.78 / 39.32 / 38.56 | 2.22 |

* The water guard is bit-identical to SGL (25.35 / 58.50 d), as expected: c_t,sub = c_t,cell.
* **SUPPORTED.** SUBGRID LEAKC is recommended for gas or CO₂ fracturing with sub-grid leak-off.

## H11: Bayesian surrogate for seismic moment (F33)

**Setup.**
* Squared-exponential GP on log M0 with ML hyper-parameters.
* The same 12 training schedules as H9b, then 12 EI steps.

**Result.**
* Best 2.2193·10¹⁵ N·m at u = [1.45, 0.55, 0.92, 1.08], reached at evaluation 14.
* That is 1.002 × the ES median best (ES used 60 evaluations × 3) and 0.941 × constant: **SUPPORTED**.

**Exploratory seed robustness.** Seeds 1–4 for the EI phase gave 2.2223 / 2.2167 / 2.2287 / 2.2163·10¹⁵ (1.001–1.006 ×).
Every optimum front-loads: u₁ 1.40–1.48 and u₂ 0.55–0.69.

**Reading.**
* F28 showed that a *physics* surrogate (Dieterich) with the wrong objective structure fails.
* A *non-parametric* surrogate with an acquisition rule succeeds, using 2.5 × fewer evaluations than one ES run.
* For P2: use the GP route for moment objectives, and the convex Dieterich route (F1/F27) for expected counts on passive
  faults.

**Limits.** The training set was not resampled; 4-D only.

## Engineering found along the way (F34)

* **BHP reference depth.** A defaulted WELSPECS reference depth now uses the well's first defined connections, SHUT ones
  included, as OPM's Well::updateRefDepth does.
  * ACTIONX_WEFAC: 9.46 % → 1.8 % (PASS).
  * ACTIONX_GCONINJE: 10.3 % → 1.04 % (PASS).
* **Group UDA.** GCONPROD and GCONINJE rate items can now be given as an FU_ or GU_ UDQ name (group UDA).
  * ACTIONX_GCONPROD: exact (0.00 %) up to 1066 d.
  * The rest comes from a deck that is sensitive to step size (a compounding UDA re-evaluated per time step).
* **Closed as unresolved (time-boxed).** KRNUM-03 directional relperm: the logic matches OPM.
* **Documented as an OPM treatment difference.** 0A1 GRCTRL STW: ResSim STW is within 0.9–2.9 % of OPM MSW, while OPM's own
  STW and MSW differ by 17 %.

## Next cycle (proposal)

* **H10f:** a reactivation-time constraint (or a smoothed permFactor ramp) that restores convexity, and how much of the design
  gain it costs.
* **H11b:** GP-EI with resampled training sets and 8 periods; compare with the SCP of F27 on passive faults.
* **H2f:** SUBGRID LEAKC on corner point with CO₂ (compositional) fracturing.
* Literature: OpenAlex re-query for "fault reactivation permeability switching convex optimisation injection".
