# ResSim Studio vs industrial simulators — developer's gap assessment (C30, October 2026)

Written from the point of view of the developer of the app. Compared against ECLIPSE 100/300 + Petrel/ResInsight (SLB),
IMEX / GEM / STARS + Builder / Results / CMOST (CMG), tNavigator (RFD), OPM Flow + ResInsight, TOUGH3/4 + iTOUGH2 +
TOUGHREACT (LBNL), COMSOL / MOOSE (multiphysics). Earlier, narrower comparison: `COMPARISON_CMG_TOUGH3.md`.

**Short answer.** The main difference is not one missing feature. ResSim has **breadth**: black oil, compositional,
thermal, CO2, EGS, geomechanics, fracture growth, phase field, adjoint history matching and a GUI, all in one
executable, which no single commercial product has. What it lacks is **depth and evidence**:

* scale (million-cell field models in minutes);
* independent verification on many real fields;
* well / facility detail (facilities; multisegment wells since C32, ACTIONX / UDQ since C31 are evaluated by the own engine);
* a QA and support process.

In the six areas the user listed, the ranking is:

| Area | Where we stand | Severity of the gap |
|---|---|---|
| 1. UI / workflow | Comparable for single-model work. Geomodelling since C32 (ResGeo) with a 3-D seismic / geomodel scene, out-of-core SEG-Y, MPS and an uncertainty loop since C36; still weaker for project / realisation data management | low-medium |
| 2. Simulation physics | Broad, often wider than any one vendor. Thinner in production engineering (facilities; MSW done in C32, UDQ / ACTIONX in C31) | medium |
| 3. Results / post-processing | Good after C30; report generator (Word / PowerPoint / HTML) added in C33. Missing: cross-run dashboards, well-log / RFT / PLT overlays | low-medium |
| 4. Engine maturity (numerics, scale) | Robust to ~10⁵ cells. Not tuned for 10⁶–10⁷ cells or real GPU clusters | **high** |
| 5. Verification | ~750 automated checks, analytic and benchmark comparisons, Norne vs ECLIPSE. No external QA, no blind field portfolio | **high** |
| 6. Customisability | Strong: open source, decks round-trip, RESSIM directives, CLI, Python package + C API with step-wise schedule steering (C33). No plug-in kernels, no PYACTION yet | low |

---

## 1. UI and workflow

**What we have.**
* Model / Deck Editor / Results / Studies workspaces, with schema-driven keyword forms (all OPM keywords) and an
  Excel-like editor for every table.
* Corner-point 3-D viewer with real well trajectories, ECLIPSE / OPM binary reader and writer, and streamlines.
* Since C30: 15 colour maps with reverse, plus *Fields & sections*: plan / row / column sections of any array,
  fracture-face maps, opening profiles, phase-field mesh views and post-processed geomechanics.
* Studies with sensitivity, uncertainty and history matching (manual, adjoint, ES-MDA, GP proxy), plus vector PDF
  export of every plot.

**What industrial tools have that we do not.**
* **Geomodelling** (Petrel, RMS, Builder): **C32 added the foundation (ResGeo)**: well / surface / fault-stick / SEG-Y
  data, seeded horizon tracking + depth conversion, faulted pillar grid with well tie, facies (SIS / TGS / channels),
  SGS with seismic co-simulation, volumetrics with P90/P50/P10, flow-based upscaling, direct simulation case.
  **C36 closed most of the remaining gaps:**
  * a 3-D geomodel / seismic scene (slices, volume rendering, attributes incl. coherence);
  * out-of-core SEG-Y;
  * flattening;
  * horizon brush editing and fault-stick interpretation in 3-D;
  * fault-parallel sloped pillars;
  * Direct Sampling MPS;
  * a structural + petrophysical uncertainty loop with tornado and exact P90/P50/P10 rebuilds.

  Still missing against Petrel / RMS:
  * truly stratigraphic grids with faults as continuous pillar surfaces through Y-faults;
  * seismic inversion / AVO;
  * facies from seismic probability cubes in MPS;
  * a project database with many realisations.
* **Project / data management**: many realisations, cases and wells, with lineage, in a database. We save single
  projects plus result files.
* **Large-model rendering**: 10⁶–10⁷ cells with filters, intersections along well paths and well-log tracks in 3-D
  (ResInsight, Petrel). Our viewer is fine to a few 10⁵ cells; a well-section curtain with log tracks exists since C31 (2-D, not yet in the 3-D scene).
* **Production-data workflows**: import of production / pressure / RFT / PLT data with QC and allocation, and
  automatic WCONHIST generation. Since C38: production / injection / pressure tables (CSV / TSV / paste) with QC,
  volume → rate with uptime (WEFAC), WCONHIST / WCONINJH merged at DATES, observed CSV for history matching (UI + CLI).
  Still open: RFT / PLT import and allocation of commingled production to layers.
* **Reporting**: a templated report (PowerPoint / Word) with fixed plots per well / group.

* **Platform and product UX (C39 audit, details in `docs/UI_REVIEW.md`).** Commercial tools run on Windows and Linux
  (Petrel: Windows; ECLIPSE / tNavigator / CMG / ResInsight: Windows + Linux) and all have undo, autosave and themes.
  **C39 closed:** Linux and macOS builds with OpenGL back ends (CI on three platforms), undo / redo of every model edit,
  autosave with crash recovery, dark / light / high-contrast themes with accent colours.
  **C40 closed:** command palette for every screen and action (G3), object tree (G4), new-model wizard and empty states
  (G7), linked report-time cursor (G9), case manager with variants / diffs / compare-with-base (G10), Inspector filter (G11),
  status icons (G12), light 3-D backgrounds; Petrel parity items log calculator / cross-plots, variogram fitting, well tie
  and geomodel parameters (variogram range, realisation seed) as Study variables.
  **Still open:** localisation (G13), the rest of the design tokens (G5), a fault-network modeller and a workflow recorder.

**Plan.** ~~Well-path intersection view and well-log tracks (C31)~~ done in C31 (Results → Fields & sections → Well section). A project tree of many cases with lineage. A report
generator (docx / pptx from saved plot selections). ~~Production-data import with WCONHIST generation~~ done in C38.

## 2. Simulation physics

**Where we are at least at par.**
* Black-oil E100 core: Norne matches ECLIPSE within 0.3 % FPR.
* EOS compositional with CO2 storage and all trapping mechanisms.
* Thermal p-N-U with IAPWS-IF97 (steam, CSS, SAGD); EGS; THM(C).
* FE geomechanics on corner-point grids with plasticity.
* DDM fault slip and seismicity.
* In-run hydraulic-fracture growth: LEFM, sub-grid tip, DDM width and phase field (one-way and two-way).
* Salt / EWASG, tracers with decay chains, polymer, dual porosity / MINC, LGR.
* Adjoint gradients for black oil, compositional and thermal.

Several of these (phase field in the field run, adjoint for thermal, DDM seismicity in the same run) are not in any
one commercial package.

**Gaps that matter for daily industrial use.**

| Missing / partial | Industry reference | Impact |
|---|---|---|
| ~~Multisegment wells in the own engine~~ | ECLIPSE MSW, CMG FlexWell | **Done in C32** (WELSEGS / COMPSEGS / WSEGVALV / WSEGSICD / WSEGAICD; open: acceleration term, drift-flux slip, segment storage, WSEGLABY) |
| ~~ACTIONX / UDQ evaluated by the own engine~~ | ECLIPSE, OPM | **Done in C31** (UDQ / UDA / ACTIONX); connections added by an action **done in C32**; open: wells created by an action, PYACTION |
| Facilities: compressors, chokes, separators in the network, gas lift optimisation | ECLIPSE network, IPSM, Petex GAP | Medium: gas-lift optimisation done in C35 (LIFTOPT / WLIFTOPT / GLIFTOPT, matches the commercial simulator on opm-tests GASLIFT-04); automatic chokes, satellites (GSATPROD), fuel / import (GCONSUMP), sales (GCONSALE) and standard / re-routed networks match OPM within 0.6 % on opm-tests network decks (C35-4b); compressors (NETCOMPA), separators and NETBALAN open |
| Miscible solvent (Todd–Longstaff MISC / TLMIXPAR) and gas-water solvent in the own engine (immiscible SOLVENT, surfactant, alkaline, foam, polymer done in C35 / C2x) | ECLIPSE, OPM | Medium (miscible decks run through OPM Flow) |
| True Killough hysteresis, Pc hysteresis | ECLIPSE | Low-medium |
| Transient wellbore (T2Well / FlexWell full), wellbore heat transfer in black oil | CMG, TOUGH | Medium for CO2 / geothermal wells |
| Reactive transport with full speciation | TOUGHREACT, GEM-GHG | Medium for CO2 / H2 storage studies |
| Dynamic gridding (DynaGrid) | CMG GEM | Low (static LGR exists; design note C29-7) |

## 3. Results and post-processing

**After C30.**
* Any 3-D array of any engine, including imported ECLIPSE / OPM runs, can be shown as 3-D, plan or vertical sections,
  with fixed / per-step / symmetric ranges and 15 colour maps.
* Hydraulic fractures have:
  * width (HFWIDTH);
  * net pressure;
  * stress shadow (HFDSNN, HFDSXX…);
  * displacement and volumetric strain around them;
  * fracture-face maps and the opening profile of every report step.
* Phase-field fractures show damage, stress change, displacement and crack-opening profiles on their own mesh.
* FE geomechanics can be post-processed for **any finished run**, including ECLIPSE / OPM results: uplift,
  displacement, effective stresses, σ_min, √J2, Mohr–Coulomb failure index, plastic strain.

**Still missing.**
* Cross-run dashboards (KPI matrix of 100 cases, ensemble fan charts of maps; we have them for summary vectors only).
* Well-section / RFT / PLT plots against observed data.
* Volumetrics reports by region and zone.
* An animation export (mp4 / gif) of the 3-D views.

## 4. Engine maturity (numerics and scale)

This is our largest real gap.
* **Scale.** Commercial engines run 10⁶–10⁸ cells on clusters or GPUs (tNavigator, ECHELON, INTERSECT).
  * Our CPR-AMG + BILU0 FIM engine with MPI domain decomposition is verified on 2–4 ranks.
  * The device CPR-AMG is verified on the D3D11 path and in float emulation, not on a production GPU cluster.
  * At 10⁶ cells memory and assembly cost are not yet engineered: no matrix-free well terms, no reordering for cache,
    no out-of-core I/O.
* **Robustness on hostile decks.** We fixed many classes (THP starts, gas injection in EDFM, hysteresis, network limit
  cycles in C30). ECLIPSE has 40 years of time-step heuristics for every pathological deck. Expect more edge cases on
  unseen field models.
* **Not developed or deliberately deferred.**
  * Fully implicit coupling of every physics at once in one Jacobian: geomechanics two-way is sequential / iterative,
    phase field is staggered.
  * Unstructured / PEBI grids (not planned; corner-point is the industry standard).
  * DynaGrid.
  * Discrete GPU backend through CUDA / HIP.

**Plan.**
1. A 1–2 M cell benchmark (SPE10 full plus a refined Norne) with a profiler-driven optimisation round: assembly
   threading, ordering, well Schur complement, I/O.
2. MSW in the own engine (it also removes the biggest physics gap).
3. A CUDA / HIP path behind the existing device interface.

## 5. Verification

**What we have.**
* ~750 automated checks:
  * analytic solutions (Buckley–Leverett, Theis / Gringarten, Terzaghi / Mandel / Cryer, Lauwerier, Sneddon /
    Griffith, PKN / KGD, Okada / TDE);
  * SPE1 / 3 / 4 / 5 / 9 / 10 / 11A;
  * OPM Flow on every example (relRMS 0.1–2 %);
  * Norne against ECLIPSE;
  * adjoint gradients against finite differences;
  * C30: post-processing = in-run coupling, traction and opening checks of the fracture field.
* MPI and MinGW builds tested every round.

**What industry has.**
* Decades of client fields, internal QA suites with thousands of decks, and regression on every release.
* Benchmark participation (SPE CSP, Code Comparison projects); for TOUGH, a published verification report per module.

**Honest gap.**
* Our tests are written by the developer of the code, so the same mind writes both and errors can be shared.
* There is no blind comparison and no external user base.

**Plan.**
1. Join an open code-comparison (SPE11 B/C full, GeoMech benchmark), with the results in the research track.
2. Run the complete opm-tests deck set nightly against OPM Flow, with an automatic difference report.
3. Publish a verification report: ENGINES_AND_MATH §14 turned into a document with one figure per test.

## 6. Customisability

**Our strengths.**
* Full source; ECLIPSE decks round-trip unchanged.
* `-- RESSIM:` directives for every non-ECLIPSE feature.
* CLI and MPI batch; studies on any deck parameter.
* Every plot exportable as vector graphics; restart files readable by ResInsight.

**What we lack.**
* A **Python API** (Petrel Ocean / CMG Python / OPM opm-common Python bindings): user scripts that read results, drive
  studies and add custom responses.
* **User-defined kernels**: MOOSE- / COMSOL-style user equations, or STARS user-defined reactions without recompiling.
* **UDQ-type user expressions** inside the run.

**Plan.** A Python module (pybind11) over Case / Simulator / results (C33), UDQ evaluation in the own engine (done in C31), and
a user-reaction table in thermal / comp.

---

## What C30 closed (this round)

| Item | Status |
|---|---|
| Network oscillation (user's screenshot) | Root cause found and fixed: zero-rate GOR/WCT limits, implicit Illinois balance (`-c30net`) |
| Visualisation of HF / EDFM / PF / geomech for all engines | Fields & sections tab (sections, fracture face, phase field, geomechanics post); HFFIELD arrays; PF arrays; 3-D fracture panels coloured by any array |
| Geomechanics for any run, incl. ECLIPSE / OPM results | `geoseq::postProcess`, button in Results and in the result viewer (`-c30post`) |
| Colour maps | 15 maps + reverse in 3-D, sections and panels (`-c30cmap`) |

## Order of the next rounds (developer's priority)

1. **C31 (done).** UDQ, UDA and ACTIONX in the own engine (ACTIONX moved forward from C32), run-time logic editor, UDQ
   results tab, well-section view with log tracks.
2. **C32 (done).** Multisegment wells in the own engine; connections added by ACTIONX at run time; geomodelling
   foundation (Geomodel workspace + standalone ResGeo, SEG-Y).
3. **C33 (done).** Python API (ctypes over a C API instead of pybind11, so no compiler is needed). Report generator. opm-tests regression: 64 / 202 pass against OPM Flow at 5 % relRMS, 30 skipped. The 1 M cell performance round moved to C34.
4. **C34 (done).** opm-tests compatibility round. Covers guide rates and group hierarchy, VREP coupling, numerical / table
   aquifers, KRNUMX/Y/Z, WELPI, rock compaction tables, EDITNNC, OPERATER and WECON workovers. 102 / 307 cases pass,
   53 are skipped.
5. **C35 (done).** surfactant, alkaline and foam in the own engine (ECLIPSE SURFACT / ALKALINE / FOAM,
   mass-conservative, with an editor), SOLVENT (immiscible, matches OPM on SPE1CASE2_SOLVENT_FOAM) and gas-lift
   optimisation (LIFTOPT / WLIFTOPT / GLIFTOPT), surface network facilities (automatic chokes, GSATPROD, GCONSUMP,
   GCONSALE, standard networks / re-routing, OPM network semantics).
6. **C36 (done).** Geomodel round 2: a 3-D seismic / geomodel scene, out-of-core SEG-Y, attributes, flattening, 3-D
   interpretation and horizon editing, sloped pillars, MPS and an uncertainty loop. MSW: labyrinth ICDs, acceleration and
   drift flux. Network: NETBALAN, automatic compressors and COMPOFF. Fixes: gas-lift NONE caps and the MSW reference depth.
   Regression Cancel.
7. **C37 (done).** PYACTION in the own engine (persistent Python child, opm_embedded shim, editor), MSW tubing flash,
   oil-water slip, segment storage, WSEGFLIM / WSEGPULL, the 1 M-cell performance round (memory -11 %, scaling to
   200 k measured, 1 M extrapolated), WCYCLE / WTEST P / GRIDUNIT (opm-tests wcycle 9/9).
8. **C38 (done).** Level-parallel float32 BILU0 (100 k cells 12.6 % faster, factor memory halved), production-data
   import with WCONHIST generation, opm-tests NONNC / UDT / EDIT TRANX / SWATINIT < SWL, OPM Newton minimum and TSFMAX 3.
   Not done (moved to C39): dropping the scalar Jacobian copy (1 M cells < 3 GB not yet reached), 9_* group control,
   NETWORK-01, WPWE, CO2STORE, KRNUM, CECON, ROCKTAB, GRCTRL clusters.

   **What C38 closed**

   | Gap | C38 result |
   |---|---|
   | No production-data workflow | table import + QC + WCONHIST / WCONINJH / WEFAC merge + observed CSV (UI and `--prod-import`); ORAT history reproduces the measured oil exactly |
   | float32 preconditioner not faster | level-parallel float BILU0: 100 k cells 58.0 → 50.8 s, identical iterations and FOPT |
   | UDT (user-defined tables) unsupported | UDT / UDTDIMS, LC / NV / LL, union operators, OPM UDQ evaluation order, UDT editor; UDT-1D-01 4 digits to t = 33 |
   | NONNC ignored | fault-throw NNCs dropped (opm-tests NONNC 13.7 → 3.05 %, PASS) |
   | EDIT TRANX applied as MULTX | neighbour connections only, exported as EDIT MULTIPLY (8_MULTIPLY 5.94 → 5.03 %) |
   | SWATINIT < SWL kept below connate | raised to SWL as OPM (5_SWATINIT 4.61 → 1.65 %) |
   | converged steps without an update, growth 2 under TUNING | at least one Newton update per step (OPM), TSFMAX default 3 |

   **What C37 closed**

   | Gap | C37 result |
   |---|---|
   | PYACTION only via OPM Flow | runs in the ResSim engine (SPE1 ACTIONX = PYACTION = legacy run(), restart exact) |
   | MSW without flash / slip / storage | tubing flash, Shi oil-water slip, segment storage, WSEGFLIM, WSEGPULL |
   | Large-model memory / profile unknown | 4 kB/cell, linear per-iteration cost, -11 % memory, 1 M = ~4 GB / ~14 min (2 cores) |
   | WCYCLE ignored | OPM semantics, 9/9 opm-tests wcycle pass |
   | GRIDUNIT ignored | CM / FEET converted (geometry exact) |
   | Injectors without cross-flow | C41: injector and producer cross-flow (OPM StandardWell); ROCK region decks fixed |
   | ROCKTABH ignored, CECON ignored | C41: ROCKCOMP HYSTER elastic curves (OPM interpolation), CECON CON / +CON / WELL |
   | No workflow automation, no fault-network QC, English-only UI | C41: workflow recorder (.rswf), fault junctions / truncation / throw QC, Spanish UI catalogue (menus) |
9. **Continuous.** Research track and code-comparison participation.
