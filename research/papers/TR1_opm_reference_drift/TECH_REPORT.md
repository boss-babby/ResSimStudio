# Technical report TR1: Reproducing OPM Flow regression references (moving reference histories and a well reference-depth subtlety)

*Draft v1.0, 2026-11-04. Intended for Zenodo (technical report) and as issues/notes to the OPM project. Not a journal paper.
Authors: [to be decided]. Generative-AI assistance was used for code, analysis and drafting; AI tools are not authors.*

## Summary
Independent simulators are often validated against the reference summary files distributed with the OPM Flow regression suite
(opm-tests). We compared an independent black-oil engine (ResSim) with those references over 209 cases. Two observations matter to
anyone doing the same:
1. Several reference histories are **moving targets**: OPM's own references for the same deck differ between repository versions
   by up to 53 %.
2. The default **bottom-hole reference depth** of a well in OPM is taken from its first connection, *including connections defined
   SHUT*. A simulator that uses the shallowest open connection instead will mismatch BHP-controlled wells by several per cent.

## 1. Method
* Cases: all 209 decks fetched by `tools/fetch_opm_tests.py` (opm-tests default folders).
* Metric: worst relative RMS over the compared summary vectors, after time alignment; pass at 5 %. A "material" variant ignores
  vectors below 1 % of their peak. The metric was never loosened during the work.
* Reference drift: the reference SMSPEC/UNSMRY of each case was taken at four opm-tests commits (2021-01 `66a7fbe991`,
  2023-04 `d5a01f9c5f`, 2024-09 `67db106e68`, 2025-05 `da782d48cd`) and compared with each other using the same metric.
* Final status of the engine (C49): 102 PASS / 77 FAIL / 30 SKIP. SKIP means unsupported features (e.g. binary IMPORT grids,
  gas-water systems without oil). Full list: `docs/regression_209_full.txt`.

## 2. Moving reference histories (F29)
| Group | Cases | Spread between OPM's own reference versions |
|---|---|---|
| stable | 6_ENDSCALE, 9_1A, 0A1, 0B_WELPI, ACTIONX_M1, ACTIONX_WEFAC, 7_HYSTERESIS, WSEGSICD, WSEGAICD | ≤ 2.5 % |
| moderate | 9_3A | 12–16 % |
| unstable | NETWORK-01-WTEST, T1L, UDQ_WCONPROD, UDQ_ACTIONX | 35–53 % |

**Bisection of two cases.**
* T1L changed with opm-simulators PR 6311.
* 9_3A changed with opm-common PR 3779, which pads low-pressure gas PVT tables (`fvf_lim = max(10·B0, 1)`). The change adds table
  rows only below 3 bar, yet the reference FWCT differs by 45 % at t = 151 d (31 vs 26 ministeps). This points to transient
  instability of the case rather than a physics change. ResSim reproduces the 2024-09 reference to 5.4 %.

**Recommendation.** Validation studies should state the opm-tests commit of their references. Cases in the "unstable" group should
not be used as pass/fail criteria without an inter-version tolerance.

## 3. Well reference depth includes SHUT connections (F34)
* **OPM behaviour** (opm-common `Well::updateRefDepth`): when WELSPECS item 5 is defaulted, the reference depth is taken once from
  connection 0, whatever its state. (We checked OPM; we did not check the ECLIPSE manual wording for this edge case.)
* **Symptom** in ACTIONX_WEFAC / GCONINJE: well OP02 has K3–4 defined SHUT, and K7–10 are opened later by WELOPEN on COMPLUMP 2.
  Using the shallowest *open* connection made WBHP 3 bar (1.2 %) high, so the 200-bar limit was reached 180 d late. Late-time
  WOPR/WGPR were 15–20 % off.
* **After matching OPM's rule:** WEFAC 9.46 % → 1.8 % (PASS); GCONINJE 10.3 % → 1.04 % (PASS).
* **Related:** GCONPROD/GCONINJE rate items may be UDQ names (group UDA). Supporting this brought ACTIONX_GCONPROD FGPR to 0.00 %
  up to 1066 d. The remaining 19 % after an ACT-03 control jump comes from the deck itself: a compounding UDA (FU_FGPR = 0.9·FGPR)
  depends on time-step size.

## 4. Reproducibility
Scripts: `tools/fetch_opm_tests.py`, `tools/dev/refdrift*`, `tools/dev/smrydiff.py`, `tools/dev/wdump.cpp`, `tools/dev/cmpw.py`.
Engine notes: ENGINES_AND_MATH.md §9.t. The opm-tests data are GPL-3.0 / ODbL as distributed by OPM. Only commit hashes and
derived statistics are published here.

## 5. Next step (human)
Open a short issue or discussion in OPM/opm-tests on the inter-version spread of the four unstable cases. Link this report and offer
the comparison script.
