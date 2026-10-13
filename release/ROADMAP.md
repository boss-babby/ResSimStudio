# ResSim Studio → OPM Flow / ECLIPSE 100 class: Roadmap

## 0. Where we stand (October 2026)

| Area | ResSim today | OPM Flow / E100 |
|---|---|---|
| Physics | 3-phase black oil, live oil (Rs), dead oil, water/gas, rock compressibility | + wet gas (Rv / VAPOIL), DRSDT/DRVDT, hysteresis, end-point scaling, polymer, solvent, thermal, CO2STORE, … |
| Grid | Cartesian, corner-point (COORD/ZCORN), ACTNUM, MULT*, faults/MULTFLT, NNC, MINPV, MULTPV, MULTREGT | + PINCH, EDITNNC, dual porosity, LGR, numerical aquifers |
| Regions | SATNUM, PVTNUM, EQLNUM, FIPNUM (multi-region tables) | + IMBNUM, ROCKNUM, MISCNUM, FIP sets, THPRES |
| Init | EQUIL (multi-region), RSVD/PBVD, enumerated PRESSURE/SWAT/SGAS/RS | + SWATINIT Pc scaling, RVVD/PDVD, threshold pressure |
| Wells | standard wells, WCONPROD/INJE/HIST/INJH, WELOPEN, rate→BHP switching, WAG phase per control | + groups (GRUPTREE/GCONPROD/GCONINJE), VFP/THP, multisegment wells, WTEST, WECON, ACTIONX/UDQ |
| Aquifers | Fetkovich, Carter-Tracy (built-in table) | + AQUTAB tables, numerical aquifers (AQUNUM) |
| Solver | fully implicit Newton, CPR + AMG (unsmoothed aggregation) + ILU0, OpenMP | + DUNE-ISTL/AMGCL, GPU (CUDA/OpenCL), MPI domain decomposition |
| I/O | 2-way .DATA (39/40 opm-tests decks import), .rsim projects | + binary restart/summary (UNRST/SMSPEC/UNSMRY/EGRID/INIT), restart from file |
| Workflow | GUI, 3D viewer, HM (LM + ES-MDA), sensitivity/MC, streamlines | ResInsight, ERT, everest (external tools) |

**Deck coverage measured on 39 opm-tests decks.** These are the keywords that were ignored or only partly honoured, ranked by how many decks use them:

VAPOIL/PVTG-Rv 17 · GRUPTREE 13 · late COMPDAT 12 · THPRES 10 · FILLEPS 8 · SWATINIT 7 · WELSEGS/COMPSEGS 5 · VFPPROD 5 · RVVD 5 · GCONINJE 5 · PINCH 4 · AQUANCON coeffs 4 · STONE1 3 · DRSDT 3 · UDQ 2 · TRACER 2 · EDITNNC 2 · AQUTAB 2 · GCONPROD 2.

These numbers set the priority order below. The rule is: implement whatever most decks need first.

---

## Status after the "core" round (October 2026)

Done in this round (all covered by the test suite, 61/61 PASS, and compared with OPM Flow 2024.10):

| Item | What was built | Check against OPM Flow |
|---|---|---|
| OPM Flow as an engine | Toolbar/Run menu engine switch (ResSim / OPM Flow). flow.exe is found automatically (Program Files\OPM Flow\bin, PATH) or picked by the user; MS-MPI check; deck export, live progress from the flow log, cancel, results imported back into the normal Results views | SPE1 via the app runner: 120 report steps, FPR/FOPT identical to a manual flow run |
| Result Viewer | EGRID/INIT/UNRST/SMSPEC/UNSMRY, unified and formatted (.F*), several cases overlaid, differences, 3D property per report step | the reference files of all 12 examples open |
| OPM compare tool | `compare_opm` (CLI) and compare in the Result Viewer | 12 examples: worst field relRMS 0.2–1.3 % (BL1D 5.8 % = OPM step artefact) |
| A2 schedule timing | late COMPDAT (opening at its report step), COMPDAT SHUT later, wells without open connections shut, ECLIPSE rule "control on an unperforated well leaves it SHUT", re-issued control on export, From/Until columns in the Wells editor | late-perforated injector + re-perforated producer: 0.14 % |
| A3 SWATINIT | per-cell Pcow scaling so that equilibrium reproduces SWATINIT | Pc + SWATINIT + THPRES model with wells: worst 1.0 % |
| A5 THPRES | EQLOPTS THPRES, explicit and defaulted (from the initial state) threshold pressures between EQLNUM regions | same model as above |
| D23 True-IMPES CPR | accumulation-based decoupling weights (new default; quasi-IMPES selectable) | linear iterations −21 to −26 % on GASCAP/REGIONS/XSEC/BL1D, +12–14 % on CONING/FAULTED; identical solutions |

Still open from Phase A: A1 wet gas (Rv, the biggest one), A4 end-point scaling, A6 PINCH/EDITNNC, A7–A10.
Next robustness items: Appleyard chopping and dS/dP-based time-step control (D24).

---

## Status after the wet-gas round (October 2026)

**A1 wet gas is done.** 72/72 tests pass (11 new). All 13 examples are compared with OPM Flow.

| Item | What was built | Check against OPM Flow |
|---|---|---|
| Rv primary variable | 3-state switch per cell (Rs-undersaturated oil / Sg three-phase / Rv-undersaturated gas), consistent wells, history RESV split with Rs and Rv, in-place split into liquid and vaporised oil | SPE3CASE1 (opm-tests, gas cycling + blow-down): FOPR 0.09 %, FPR 0.04 %, FGOR 0.11 % relRMS |
| PVTG | full table with undersaturated branches; records without a branch are extended from the next higher-pressure record (OPM rule); "right-extreme" interpolation that keeps the distance below the dew-point Rv; Bg and mu interpolated separately | same deck (a plain bilinear interpolation gave 2.2 %) |
| Equilibration | RVVD, constant Rv, or saturated Rv at the contact (ECLIPSE default); **gas directly on water** (GOC = OWC) now uses the gas-water Pc = Pcow + Pcog for Sw, as OPM does | initial FOIP 0.05 %, FWIP was 6 % off before the gas-water rule |
| Deck I/O | VAPOIL, PVTG, RVVD, EQUIL item 8, RV array: import + export + validation messages | CONDENSATE example exported and run by OPM: 0.56 % |
| Wells | a rate target of zero stops the well (no residual cross-flow), as in ECLIPSE | FGIR after injection stop: exactly 0 |
| UI | VAPOIL switch, PVTG editor (records + undersaturated branch + plots), initial Rv (saturated / constant), RSVD and RVVD table editors, vaporised oil and total oil in the in-place reports, RV in the 3D view (own engine and OPM restart files) | |
| Example | "Gas condensate cycling (SPE3)" in the example library | |

## Status after the keyword pass-through round (October 2026)

**Any OPM Flow deck now survives import → edit → export → run unchanged.** 78/78 tests pass (6 new).

| Item | What was built | Check against OPM Flow |
|---|---|---|
| Raw keyword pass-through | Keywords the ResSim model does not hold (groups, VFP tables, MSW, ACTIONX…ENDACTIO, UDQ, OPM-only options (tracers, polymer, ...), extra SUMMARY requests, RUNSPEC dimensions, unknown keywords) are kept as `ExtraKeyword {name, section, report step, body}` in deck order. Partly modelled keywords keep the unmodelled records (e.g. WELOPEN on completion numbers). Saved in the project file. | opm-tests decks exported and run by OPM, compared with the original run (all summary vectors): **SPE9_CP_GROUP 3.8e-8 %**, **3D_MSW 3.1e-7 %**, **ACTIONX_GCONPROD 0 %** (bit-identical) |
| Step ordering | Each schedule extra remembers whether it came before or after the well controls of its report step, so VFP tables are defined before WCONPROD refers to them and ACTIONX/UDQ come after | the "VFP table 5 not defined" failure is gone |
| Wells | WCONPROD/WCONINJE GRUP, THP, NONE, CRAT, RESV modes; THP limit; VFP table; ALQ are kept for OPM (the own engine uses the best explicit target). **Connection-level WELOPEN (I J K, 0 = all)** is modelled natively. SHUT connections added later by COMPDAT are kept. WELSPECS items 7+ and COMPDAT sat-table, D-factor and r0 are kept. The "OPEN request ignored" rule now applies only to wells that have no connections at all (as in OPM). | ACTIONX_GCONPROD: FOPT 1.72718e6 in both runs (it was 0 before the WELOPEN fix) |
| UI | **Advanced keywords** page: filter, add (by section), reorder, delete, edit the body in a monospace editor with live syntax check, report-step slider, "before/after the well controls" order, support level and the keyword manual. Schedule editor: "OPM Flow control" block (mode sent to OPM, THP, VFP table, ALQ). Wells editor: WELSPECS items 7+, completion Sat-table / D-factor / r0 columns. The own engine lists the OPM-only keywords in its run log. | |

## Status after the editor / results / EPS round (October 2026)

85/85 tests pass (7 new). Norne no longer crashes on the Wells page.

| Item | What was built | Check |
|---|---|---|
| Norne Wells crash | ImPlot colormap lookup received NaN (wells outside the active grid / defaulted values): NaN guard in implot `LerpTable`, NaN cells drawn grey in the wells map | Norne (36 wells, 44k active cells) Wells page |
| Workspaces | Toolbar switch **Model / Deck Editor / Results** (Ctrl+1/2/3, View menu). The docking layout of the Model workspace is kept alive while another workspace fills the window. | |
| Deck Editor (writing mode) | Full-window editor (ImGuiColorTextEdit, MIT): ECLIPSE/OPM highlighting, keyword / value / well / group autocomplete, hover docs with the active parameter (0.45 s rest delay), live lint (unknown keyword, wrong section, too many items, missing `/`, column 1), F7 full parser check, outline + keyword browser, docs and **Form** tabs, Align columns, record templates (auto on a new keyword), Ctrl+click INCLUDE, multi-tab, Apply to model / Run. "Model deck" writes the current case with large arrays as INCLUDE files. Keyword data: opm-flow-editor-support index (3161 keywords). | |
| Schema-driven keyword forms | `deckintel::Form`: records as tables with parameter names, units, defaults and option lists from the keyword index; shared by the editor's Form tab and the **Advanced keywords** page (Form / Text toggle) | |
| Group editor | Wells page: group tree FIELD → groups → wells (WELSPECS + GRUPTREE), parent changes and new groups write GRUPTREE, well move, list of group controls (GCONPROD, GCONINJE, GCONSUMP, ...) with Edit (schema form) and templates | run by OPM Flow |
| Results workspace | Run dashboard + Result Viewer in a full window. **Distinct curve colours** (palettes: default, Okabe-Ito, Tol bright/muted, Tableau, grayscale; golden-angle hues after the palette), colour by case / vector / curve. **Plots fit their space** (grid chosen from the available size, scroll only below a minimum size), SI tick labels. | |
| Publication figure export | Every plot set records what it draws; Export window with journal presets (85 / 120 / 174 mm, thesis, slide, poster), font family/size, line width, palette, legend position, grid, minor ticks, editable titles / axis labels / legend entries / colours / dashes, live WYSIWYG preview. **Vector PDF** (base-14 fonts) and **SVG**, **PNG** 150-1200 dpi (Direct2D/WIC), CSV of the data. | PDF/SVG rendered and checked with MuPDF |
| End-point scaling (own engine) | ENDSCALE + SWL, SWCR, SWU, SGL, SGCR, SGU, SOWCR, SOGCR: ECLIPSE two-point kr scaling, Pc scaled between SWL-SWU / SGL-SGU, 3-phase default kro with cell SWL, equilibration with the scaled Pc (and SWATINIT). Import (BOX/EQUALS/COPY work), export (defaulted cells get the table end point), SCAL page tab (constant / per layer / clear, scaled-curve preview). | identity scaling reproduces the unscaled run exactly; Norne: all 8 arrays read and applied. Comparison with OPM Flow on an EPS deck is still to do (needs flow.exe) |

---

## Status after the EOR / Phase-A round (October 2026)

| Item | Result |
|---|---|
| Polymer flooding (own engine) | Sequential-implicit polymer transport after OPM blackoilpolymermodules: Todd-Longstaff mixing, irreversible/reversible adsorption, RRF, dead pore volume, PLMIXNUM regions, PLYSHLOG shear thinning (cells + well injectivity). POLYMER-02/03 and the new POLYMER example match OPM Flow (0.7 % worst rate). |
| Polymer UI | EOR section (viscosity, adsorption & rock, mixing, shear, injection slug helper, results), POLYMER/PADS in the 3D viewer, FCIR/FCPR/FCIT/FCPT/FCIP/FCAD/WCIR/WCPR in Results, validateCase checks. |
| Equilibration + wellbore | OPM oil-pressure rule in the equilibration; per-connection segregated well heads (OPM StandardWellConnections). POLYMER-03 cumulatives within 1 % of OPM. |
| STONE1 / STONE2 | OPM EclStone1/2Material formulas, EPS-aware; GASCAP variants within 1 % of OPM on FOPT. |
| DRSDT | lastRs + DRSDT·dt limit with ALL/FREE, primary-variable switching respects the cap; SPE1 + DRSDT 0 within 0.04 % of OPM. |
| Tests | 105 native checks pass. |

## Status after the polish / studies round (October 2026)

| Item | Result |
|---|---|
| SPE10 Model 1 (100×1×20, gas injection) | Runs in both engines. FOPT 75 661 (own) vs 75 660.9 (OPM Flow). The deck parser now accepts section headers with trailing `=====`. OPM Flow needs `--shut-unsolvable-wells=false` together with `--min-time-step-before-shutting-problematic-wells-in-days=1e-12`; the OPM engine settings have a "Keep problematic wells open" switch that passes both. The deck is in `examples/spe10/`. |
| One 3D viewer for Model and Results | Shared viewer_ui: the same property/report/filter/display/slice/opacity/streamline options in both workspaces. Real deviated/horizontal well trajectories. Clickable well heads and name tags focus the camera on the well. CMG-style well-focused streamlines (injector→producer, coloured by time of flight). |
| Deck editor | Keyword docs embedded in the exe (3 161 OPM keywords, RCDATA). Right-click menu: docs, form, record template, align columns, select block, cut/copy/paste, comment, fold/unfold all, find. Snippet tab-stops only for newly inserted snippets. Keyword- and section-level folding. Icon toolbar. |
| Icons + modern UI | Embedded Lucide icon font (ISC licence): activity bar (Model / Deck / Results / Studies / Engine / About), toolbar, menus, section list, status bar. App icon (.ico). |
| Studies workspace (Ctrl+4) | Dedicated screen for sensitivity + history matching. Methods: tornado, Monte Carlo, experimental design (LHS / 2-level factorial / CCD / Plackett-Burman / random) + quadratic proxy with R²/Q² and a proxy optimum, Levenberg-Marquardt, ES-MDA, and manual sliders. 24 parameter kinds and an observed-data editor with objective grouping. Views: profiles (archive, P10-P90 band, observed error bars), cross-plot, Pareto front, tornado, MC histograms/fan, convergence, proxy, log. The sortable run archive has CSV export and "apply to case". |
| ECLIPSE as a third engine | Run > Engine: ResSim / OPM Flow / ECLIPSE. Auto-detects a licensed installation: path or folder in the settings, then %ECLIPSE_EXE% / %ECLRUN% / %ECLPATH%, PATH, C:\\ecl … F:\\ecl (macros\\eclrun.exe, or the newest <version>\\bin\\pc_x86_64\\eclipse.exe), Program Files\\Schlumberger and \\SLB. Browse is the fallback. Uses `eclrun [-v VERSION] eclipse/e300 CASE`, with the version list from `eclrun --report-versions`. A direct eclipse.exe gets the case name on stdin. Live progress comes from the STEP/TIME lines; @--ERROR/@--BUG messages and licence failures are reported (with a .PRT fallback). Results are read back like OPM runs. A Windows job object makes Stop end the whole process tree. Tested with mock launchers only, because no ECLIPSE licence is available here. |
| Tests | 115 native checks pass (also under Wine). |

## Status after the robustness / dual-porosity round (October 2026)

| Item | Result |
|---|---|
| External runs for any deck | **Pass-through mode** (Run > Engine > deck): while the imported case is unedited, OPM Flow / ECLIPSE run the *original* deck file in its own folder, so INCLUDEs, keywords ResSim does not model and section order stay exactly as written (SPE9, Norne, SPE10). After an edit, the generated deck is used. Generated decks are also more robust now. Lines are wrapped at 132 characters (the ECLIPSE limit). WCONPROD `NONE` is written as `1*`. *DIMS keywords are written as max(source, model need). PVDG is repaired to strictly decreasing Bg (PVDO to non-increasing Bo), with a comment. Keywords found in a section where they are invalid (e.g. summary vectors with no SUMMARY header) are moved to their own section. ECLIPSE `.PRT` errors (`@--  ERROR`) are parsed and shown. |
| OPM vs ECLIPSE on Norne | Run through pass-through. relRMS over about 430 d: FOPR 2.95 %, FPR 1.86 %, FOPT 1.56 %, FGOR 6.1 %, FWPR 9.3 %. Worst well vector: WWPR:B-4H at 22.6 %. |
| Dual porosity / dual permeability | DUALPORO / DUALPERM in the own engine. NZ is doubled (matrix layers first, then fracture). Matrix-fracture transfer is T = CDARCY·K·V·σ·MULTMF, using SIGMA/SIGMAV or LTOSIGMA with LX/LY/LZ. Fracture permeability is multiplied by fracture porosity unless NODPPM is set. DUALPORO has no matrix-matrix flow; wells connect to fracture cells. Grid editor UI and an example "dualporo" are included. OPM Flow does not support DUALPORO, so the app blocks that combination with a clear message. ECLIPSE runs it. |
| PINCH | Thin cells (below the threshold) are deactivated. Vertical pinch-out NNCs connect across runs of inactive or thin cells, with GAP/NOGAP, a maximum gap and MULTZ TOP or ALL. The connections never cross the matrix/fracture boundary. The keyword round-trips through the deck, and the Norne re-export runs in OPM. Not yet supported: PINCHREG, PINCH item 4 = ALL (warned on), and EDITNNC. |
| ECLIPSE streamlines | FrontSim (ECLIPSE's streamline simulator) is the third ECLIPSE program option, next to E100 and E300. In the 3D viewer, streamlines on external results use the simulator's own inter-cell fluxes (restart FLR{OIL,WAT,GAS}{I,J,K}+). These are present when the deck has `RPTRST FLOWS` (OPM Flow also needs `FLORES`). Without them, the model fluxes are rebuilt from PRESSURE/SWAT/SGAS. Against OPM fluxes the mapping gives corr 1.00000, rel L2 1.1e-5. The streamline popup shows which source is used. |
| History matching UI | The Studies workspace has a Sensitivity / History matching mode toggle. Observed data can be imported straight from a summary file (history vectors WOPRH/WWPRH/… or simulated vectors, with well/quantity filters and time thinning). New views: **Observed** (data vs best/initial run), **Match quality** (per-well/quantity normalised misfit table + heat map, before/after), and **Parameters** (prior range bars with best value, apply best, ES-MDA posterior histograms, LM path). |
| UI | Icons are now an embedded icon font with the Segoe PUA excluded, so they no longer show as broken glyphs on real Windows. The activity bar is removed (workspaces via the toolbar/View menu, Ctrl+1..4). Settings > Interface size: 60–160 % slider with 80/90/100 % presets (default auto-scales on smaller screens). |
| Threading | Measured on SPE9 (9 000 cells): 17.3 s on 1 thread vs 12.5 s on 2 threads. Small models do **not** run better serially, so the planned "serial below 20k cells" rule (D25) was dropped. |
| Tests | 138 native checks pass (also under Wine). |

**What is still open in Phases A–D (honest list):** groups with guide-rate allocation (B12), VFP/THP (B13), MSW (B14), WECON/WTEST/ACTIONX/UDQ (B15/B16) and tracers (B17) are run by OPM/ECLIPSE through pass-through, but not by the own engine. Also open: the binary output writer (C18), RESTART (C19), EDITNNC/PINCHREG, AQUTAB/numerical aquifers, hysteresis, block-CSR/AVX2 assembly, AMG reuse, and the GPU backend. These are multi-week items each and were not finished in this round.

## Industry-grade gap audit (living list, updated every round)

Legend: ✅ done · 🟡 partial · ❌ missing. "Own" is the ResSim engine. "OPM GUI" means the case can be built and edited in the app and run by OPM Flow.

| Area | Own engine | OPM GUI (build in the app) | Next action |
|---|---|---|---|
| Black-oil PVT (PVTO/PVDO/PVTG/PVDG/PVTW, DISGAS, VAPOIL) | ✅ | ✅ | PVCDO, PVZG, RSCONST |
| Multiple PVT/SAT/EQL/FIP regions | ✅ | 🟡 region tables are editable only for region 1 | region selector in the PVT/SCAL/Init editors |
| SCAL (SWOF/SGOF), 3-phase kro: default, STONE1 (+STONE1EX), STONE2 | ✅ (Stone within 1 % of OPM) | ✅ model combo in Rel. perm. | SWFN/SGFN/SOF3 family, SLGOF |
| End-point scaling (SWL, SWCR, …, ENDSCALE) | ✅ two-point (kr + Pc) | ✅ | 3-point (SCALECRS), KRW/KRO/KRG/PCW value scaling, directional; verify against OPM Flow |
| Hysteresis (EHYSTR), capillary-pressure scaling (PCW) | ❌ | ❌ | after A4 |
| Grid: Cartesian, corner-point, faults, MULTFLT, NNC, MINPV, NEWTRAN, PINCH | ✅ | ✅ | EDITNNC, PINCHREG (A6), MULTREGT, LGR |
| Dual porosity / permeability (DUALPORO, DUALPERM, SIGMA(V), LTOSIGMA, MULTMF, NODPPM) | ✅ | ✅ (ECLIPSE; OPM Flow has no DUALPORO) | gravity drainage (GRAVDR), DPKRMOD, verify against ECLIPSE |
| Equilibration: EQUIL, RSVD, RVVD, SWATINIT, THPRES, enumerated | ✅ | ✅ | PBVD/PDVD, THPRES IRREVERS |
| Aquifers: Fetkovich, Carter-Tracy | ✅ | ✅ | AQUTAB, numerical aquifers |
| Wells: WCONPROD/INJE/HIST/INJH, COMPDAT (timed), WELOPEN (well + connection level), WPIMULT | ✅ | ✅ | COMPLUMP numbers (kept raw), per-connection sat-table in the own engine |
| Groups (GRUPTREE, GCONPROD, GCONINJE) | ❌ wells run on their own targets | ✅ group tree editor + schema forms, run by OPM | B12: native group rate allocation (guide rates, actions), the next big physics item |
| VFP tables, THP control | ❌ | 🟡 VFP tables via schema form / Deck Editor; THP/VFP/ALQ in the schedule editor | B13: VFP table viewer (curves) + THP control in the own engine |
| Multisegment wells | ❌ | 🟡 WELSEGS/COMPSEGS/WSEGVALV as text, run by OPM | B14 |
| Sensitivity / history matching | ✅ Studies workspace (tornado, MC, DoE + proxy, LM, ES-MDA); HM mode with summary-file observed data, match-quality table/heat map, parameter posteriors | ❌ studies run on the own engine only | OPM Flow as the study engine (parallel batch runs), MCMC / adaptive proxy, multi-objective GA |
| WECON, WTEST, WELTARG, ACTIONX, UDQ | 🟡 WELTARG only | 🟡 as text, run by OPM (UDQ expressions kept verbatim) | B15/B16 |
| Rs/Rv rate limits (DRSDT ALL/FREE; DRVDT) | ✅ DRSDT (SPE1 + DRSDT 0 within 0.04 % of OPM); ❌ DRVDT | ✅ DRSDT in the PVT section | DRVDT, DRSDTCON |
| Polymer (POLYMER, PLYVISC/PLYADS/PLYROCK/PLYMAX/PLMIXPAR/PLMIXNUM/PLYSHLOG/WPOLYMER/SPOLY) | ✅ Todd-Longstaff, adsorption, RRF, DPV, shear thinning, well injectivity (within 0.4–2 % of OPM) | ✅ EOR section, POLYMER/PADS maps, FC*/WC* vectors, example POLYMER | SHRATE form, PLYATEMP, brine coupling |
| Tracers, solvent, thermal, CO2STORE, foam, brine, surfactant | ❌ | 🟡 keywords kept as text and run by OPM; schema forms | next EOR modules after groups |
| Unknown / unsupported keywords | ignored by the own engine (listed in the log) | ✅ kept, editable, written back in order; unedited imports run the original deck (pass-through) | — |
| Streamlines | ✅ own TOF tracer (Pollock) | ✅ ECLIPSE FrontSim engine option; viewer uses simulator FLR* fluxes | streamline-based saturation transport in the own engine |
| Schema-driven keyword forms (all OPM keywords) | n/a | ✅ forms for record keywords and tables (Deck Editor, Advanced page) | forms for array keywords (box editor) |
| Binary output (EGRID/INIT/UNRST/SMSPEC) | 🟡 reads them; writes its own .rsim | ✅ via OPM | writer for the own engine (C18) |
| Restart (RESTART, SAVE) | ❌ | ❌ | C19 |
| Time-step control (dS/dP targets), Appleyard chopping | 🟡 Appleyard-style dS/dp chopping per Newton update, Newton-count-based step growth/cut; no dS/dP target control yet | n/a | target-change step control, well Schur complement (D24) |
| Linear solver (CPR true-IMPES + AMG, ILU0) | ✅ | n/a | AMG reuse, block-CSR (D22/D23) |
| Million-cell performance | 🟡 100k cells works, slowly | ✅ via OPM | Phase D |

**Gridding decision:** corner-point + NNC (what ECLIPSE/OPM/Petrel use) stays the standard. Still to add: PINCH/MINPV pinch-out NNCs and EDITNNC (A6), then CARFIN/LGR (Phase E). PEBI/unstructured grids are not needed for ECLIPSE-class work and are not planned.

**Plan for the next rounds, in order:**
1. ~~Raw keyword pass-through~~ done.
2. ~~Schema-driven keyword forms, Deck Editor, group editor, end-point scaling~~ done (see above).
3. Own engine: group rate allocation (GCONPROD/GCONINJE with guide rates), then VFP tables + THP control, PINCH/EDITNNC, dS/dP time-step control, MSW, ACTIONX.
4. EPS check against OPM Flow (Norne and an SPE1 variant), hysteresis (EHYSTR).
5. Cross-platform port (GLFW/SDL + OpenGL/Vulkan, portable dialogs and process launch).

---

## Phase A: Correctness parity with E100 black oil (≈ 3–4 months)

The goal is that SPE1/3/5/9, model1–6, Norne and Drogon match OPM Flow within 1–2 % on FOPR, FWCT, FGOR and FPR.

1. ✅ **Wet gas (Rv, VAPOIL, PVTG, RVVD)** (done, SPE3 within 0.11 % of OPM Flow). This was the biggest gap (17 decks).
   - Add an Rv primary variable and a 4-state variable switch (Sg / Rs / Rv / undersaturated gas).
   - Add the PVTG table with undersaturated branches.
   - Apply RVVD in equilibration.
   - Adds a gas-condensate example to the menu.
2. ✅ **Schedule timing** (late COMPDAT done; still open: WELOPEN on single connections, COMPLUMP). COMPDAT/WELSPECS added after t=0 must open at their own report step (12 decks today treat them as present from time zero). Also implement WELOPEN on single connections, COMPLUMP and the time-dependent well list.
3. ✅ **SWATINIT** (7 decks, done): scale Pc per cell so that equilibrium reproduces SWATINIT, as E100 does (PPCWMAX).
4. **End-point scaling** (FILLEPS, SWL/SWCR/SWU/SGCR/SOWCR/KRW/KRO arrays, ENDSCALE): two-point and three-point scaling.
5. ✅ **THPRES / threshold pressure** (done; IRREVERS is not distinguished) between EQLNUM regions (10 decks).
6. 🟡 **PINCH and EDITNNC** (PINCH done; EDITNNC/PINCHREG open). Generate pinch-out NNCs across thin or inactive layers (PINCH, MINPV, PINCHREG), and apply EDITNNC multipliers.
7. ✅ **Rs/Rv dynamics:** DRSDT done (ALL/FREE); DRVDT still open.
8. ✅ **Three-phase relative permeability:** STONE1/STONE2 done, next to the ECLIPSE default model.
9. **Aquifers:** AQUTAB influence tables, AQUANCON influx coefficients and multipliers, and numerical aquifers (AQUNUM/AQUCON).
10. **Units:** the LAB unit system (PINCH05 is the last failing deck).
11. **Benchmark harness.** Build a nightly run of every opm-tests deck, compared against reference OPM Flow summary files (opm-tests ships `*.SMSPEC/UNSMRY` for many cases). Report the RMS deviation per vector and fail on regressions.

## Phase B: Production wells and field management (≈ 3 months)

12. **Groups:** GRUPTREE, GCONPROD and GCONINJE (rate targets, guide rates, voidage replacement), GEFAC/WEFAC, and group-level summary vectors (GOPR, GWCT, …).
13. **VFP tables and THP control:** VFPPROD/VFPINJ, WCONPROD THP limit, and a 5-D interpolation with ALQ/GFR/WFR axes.
14. **Multisegment wells:** WELSEGS/COMPSEGS, a segment pressure-drop model (hydrostatic + friction), WSEGVALV and ICD/AICD.
15. **Economic and testing controls:** WECON, WTEST, CECON, WELTARG and WLIST.
16. **ACTIONX + UDQ.** Add an expression evaluator over summary vectors, plus ACTIONX blocks that re-enter the schedule (2 decks today; essential for real field decks).
17. **Tracers:** TRACER/WTRACER passive tracers, solved explicitly after each step.

## Phase C: Industrial I/O (≈ 1.5 months)

18. **Binary ECLIPSE output** (via `resdata`-compatible writers):
    - EGRID, INIT, UNRST, SMSPEC and UNSMRY.
    - Results open directly in ResInsight and existing post-processing tools.
19. **Restart:** RESTART from UNRST (SOLUTION/RESTART keyword) and SAVE/LOAD of the run state.
20. **Complete SUMMARY section:** honour the requested vectors (field/group/well/connection/region/block), RPTSCHED/RPTRST frequency and the PRT report.
21. **Deck editor in the app:**
    - Syntax-highlighted .DATA editing with INCLUDE navigation.
    - Live parse errors and keyword hover docs (the keyword reference window already exists).
    - "Jump to line" from the import report.

## Phase D: Performance at million-cell scale (≈ 3 months, partly in parallel)

Today Drogon (46×73×31, about 100k cells) takes 47.6 s for step 1. Profiling shows spmv at 36 %, ILU0 at 27.5 % and about 27 linear iterations per Newton.

22. **Assembly:**
    - Block-CSR (3×3) storage with SIMD (AVX2) block kernels.
    - A colour-partitioned parallel Jacobian with no atomics.
    - Automatic differentiation with fixed-size duals (already AD-based; remove the heap allocations).
23. **Linear solver:**
    - ✅ True-IMPES / quasi-IMPES CPR weights (done).
    - AMG setup reuse across Newton iterations and time steps.
    - ILU0 with level scheduling, or a block-Jacobi + ILU0 smoother per thread.
    - Target: under 10 linear iterations per Newton.
24. **Nonlinear:**
    - Appleyard chopping and saturation-dependent damping.
    - Time-step control from the target saturation and pressure change (TUNING honoured safely).
    - Well-equation Schur complement.
25. **Threading heuristics:** ~~serial below 20k cells~~ dropped: SPE9 (9k cells) runs 1.4× faster on 2 threads. Add NUMA-aware first touch for large models.
26. **GPU (optional backend):**
    - Port the CPR stages (spmv, ILU0 apply, AMG V-cycle) to a D3D11/D3D12 compute or CUDA path behind the existing solver interface.
    - Expected gain: 3–6× on more than 1M cells.
27. **Distributed memory:** METIS/Zoltan partitioning, halo exchange and an MPI build. This is a long-term item and only matters past about 5M cells.

Performance targets (8-core desktop): SPE9 under 2 s, Norne full history under 3 min, Drogon full history under 10 min, and 1M cells × 10 years under 1 h.

## Phase E: Beyond E100 (later)

- Hysteresis (Killough, Carlson) and IMBNUM.
- Dual porosity / dual permeability (DUALPORO, SIGMA).
- ✅ Polymer (done). Solvent (4-component), brine, surfactant, foam: next EOR modules.
- Thermal (E300-like energy equation).
- CO2STORE / H2STORE for CCS.
- Local grid refinement (LGR, CARFIN).
- Compositional (E300 class: EOS flash, K-values). This would be a new simulator core.

---

## Verification strategy (applies to every phase)

1. **Unit tests** per feature (`tests/run_tests.cpp`): material balance below 1e-6 and analytical solutions (Buckley-Leverett, Fetkovich, Theis).
2. **Round trip:** every Case → .DATA → Case must be text-identical. This is already enforced for all 12 examples.
3. **Deck regression:**
   - All opm-tests decks run nightly.
   - Summary vectors are compared with OPM Flow reference results.
   - Tolerances: 1 % (black oil) and 3 % (when features are approximated).
4. **Performance regression:** wall-clock time per deck, alerting on a slowdown of more than 10 %.

## Suggested order and milestones

| Milestone | Contents | Exit criterion |
|---|---|---|
| M1 | A1–A3, A11 | SPE3 (wet gas), Drogon and Norne within 2 % of OPM |
| M2 | A4–A10 | all 40 opm-tests decks import with no "ignored" physics warnings except the Phase B/E items |
| M3 | B12–B13, C18 | Norne runs with groups and VFP; results open in ResInsight |
| M4 | D22–D25 | 5× speed-up on Drogon; 1M-cell synthetic case under 1 h |
| M5 | B14–B17, C19–C21 | full Norne/Drogon schedules incl. MSW, ACTIONX, restart |
| M6 | D26, E items | GPU backend, hysteresis, dual porosity |
