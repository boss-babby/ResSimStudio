# Research Cycle 10 plan (C48, 2026-11-03) — pre-registered before any run

Builds on Cycle 9 (F25–F30). The three questions left open there:
* **H10e:** does a *threshold* nonlinearity, from two-way fault permeability, break convexity?
* **H2e:** does a gas-aware sub-cell storage remove the t90 over-correction of F30?
* **H11:** can a cheap surrogate design for seismic moment where the Dieterich surrogate failed (F28)?

## H10e: convexity with two-way fault-permeability switching

**Why.** F26 showed that smooth multiphase nonlinearity is small and mostly favourable. Two-way coupling is different: a
fault face whose CFF reaches 0 multiplies its transmissibility (FAULTMECH item 9). That is a discontinuous map
u → pressure, and it is the case most likely to break convexity.

**Cases (h10_mkdeck.py):** each sealing case uses MULTFLT F1 0.02 and permFactor 50, so a reactivated face reopens to
T × 1.0.
* NS: CO₂, fault column 15.
* NHS: CO₂, column 15, 3× rate.
* WNS: water, column 15.
* WNHS: water, 3× rate.

**Pre-check.** In a single u = 1 run, reactivated faces > 0 at some time (from `--fault-csv`, cff ≥ 0). A case without
reactivation is reported as non-discriminating and dropped.

**Test.** M = 40 midpoint pairs per case (seed 13, same sampler and tolerance ε as H10). Exact decomposition δ = NL + CV
(H10d method) for the first 10 pairs.

**Criteria (as H10).** ≤ 5 % violations → SUPPORTED; > 15 % with median δ of the violators > 1 % → REFUTED; otherwise
INCONCLUSIVE.

**Stated expectation.** Violations appear where the switch occurs at different times for A, B and C (NL > 0 can exceed
abs(CV)).

**Amendment A1 (C48, after the u = 1 pre-check, before any pair run).** With the default stress state no case reactivates:
the most-stressed face stays 30–37 bar below failure (max CFF −29.97 NS, −20.74 NHS, −24.35 WNS, −11.47 WNHS), so all four
cases were non-discriminating. Re-specification: a uniform shear pre-load c (entered as FAULTMECH cohesion −c; slip is off, so
cohesion only shifts CFF) chosen from the u = 1 run alone, so that the most-stressed face would reach CFF ≈ +2 bar:
c = 32 / 23 / 27 / 14 bar (NS / NHS / WNS / WNHS). Re-check at u = 1: 8 / 7 / 6 / 3 reactivated faces, max CFF +1.26 / +0.92 /
+1.61 / +0.90 bar, and max dCFS fell (NS 6.02 → 5.25 bar) because the reopened fault passes pressure — two-way coupling is
active in both the compositional (CO2STORE) and black-oil paths. Criteria, M, seed and ε are unchanged.

## H2e: gas-aware sub-cell storage for sub-grid leak-off

**Why.** In F30 (gas, 600 sm³/d), SGL is 5 % fast at dx 10 m, and t90 converges from below while SG converges from above.

**Change (engine option, default off).** Read the HF sub-grid leak-off code first. The intended change: the storage of the
exposed sub-cell includes fluid compressibility (φ·c_t·V_sub, with c_t from the cell's phase compressibilities), not only the
pore volume fraction.

**Test.** The H2d′ deck at 600 sm³/d on three grids (10 / 5 / 2.5 m), for SGL-c (new) against SGL and SG.

**Criteria.**
* SUPPORTED: the SGL-c spread is < 1 d at both t60 and t90, and the t90 spread is below the SG spread.
* REFUTED: the SGL-c t90 spread ≥ the SGL t90 spread (2.19 d).
* Otherwise INCONCLUSIVE.

**Guard.** The water case (F21 deck, dx 10 / 5 / 2.5) must stay within 0.3 d of the F21 SGL values.

*If reading the code shows that the over-correction cannot come from storage (for example, the sub-cell has no storage
term), H2e is re-specified here before any run.*

## H11: moment-aware schedule design on the rate-and-state fault

**Why.** F28: the Dieterich surrogate optimum is 1.83× the ES best for the ensemble-mean M0.

**Method.**
* Fit a Gaussian-process surrogate (squared-exponential kernel, ML hyper-parameters) to log M0 on the same 12 training
  schedules as H9b.
* Then run expected-improvement steps up to a total budget of 24 RSF evaluations. H9b gave the surrogate 13 and ES 60 × 3.

**Criteria (ES median best from H9b = 2.215·10¹⁵ N·m).**
* SUPPORTED: best ≤ 1.05 × ES median best within 24 evaluations.
* REFUTED: > 1.15 ×.
* Otherwise INCONCLUSIVE.

The constant schedule (2.358·10¹⁵, 1.064 ×) is the reference that any method must beat by more than 1 % to count as useful.
