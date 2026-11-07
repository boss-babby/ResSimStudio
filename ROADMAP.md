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

## Status after the solver-framework / lenient-import round (October 2026)

| Item | Result |
|---|---|
| Why Norne stalled in the own engine | The Newton iterations were fine. The linear solves failed: scalar ILU(0) diverged (residual 2e13) on cells where a phase vanishes (gas cap, water zone), because their scalar diagonal is ~0. CPR's second stage is that ILU(0), so CPR-GMRES also failed (500 iterations). Every failed solve forced a time-step cut. |
| New linear-algebra framework (`bsr.h`) | **BSR** block matrices (nb x nb blocks, the native layout for black-oil, polymer, thermal and compositional systems, and for GPU kernels). **BILU0** inverts each diagonal block exactly. **WellSchur** eliminates the well BHP unknowns exactly (S = Acc - Acw D^-1 Awc) and recovers the BHPs afterwards. **BlockCPR** does quasi/true-IMPES decoupling, then an AMG V-cycle, then BILU0. **AMG reuse**: `AMG::refresh` keeps the aggregates and recomputes the Galerkin operators (full setup at the start of each time step, refresh inside Newton). The Krylov solvers (BiCGSTAB, flexible GMRES) now take any `LinearOperator`. On the Norne system that failed: BILU0 alone converges in 29 iterations, block CPR in 13 iterations (0.49 s), and an AMG refresh costs 0.05 s. The scalar ladder stays as a fallback. Tuning flag: `blockSolver` (default on). |
| Nonlinear toolkit (`nonlinear.h`) | `Dual<N>` (forward AD) with abs/min/max/smoothMax/pow/tanh. ForwardDiff-style `jacobian()` / `derivative()`. `localNewton()`: exact AD Jacobian, pivoted LU, Armijo line search and bounds, for flash, well and thermal sub-problems. `OscillationDetector`: OPM-style relaxation of oscillating or stagnating Newton updates, now active in the engine (tuning `relaxOscillations`). `PidStepController` for change-based time-step control (not yet wired into the engine). |
| Norne in the own engine | **Runs to the end:** 3312 d, 247 report steps, 2823 Newton iterations, 15 cuts, 23 min on 2 threads. Before, it stalled in the first Newton solve. Zero linear-solver failures. Second fix: the producer rate-control equation counted connections with blocked back-flow, so targets looked unreachable (switch to BHP) and then over-met (switch back). That caused 1208 control flip-flops on D-2H. The equation now counts only flowing connections. A well moved to BHP by the stall rule stays there until the next report step. vs ECLIPSE at 3312 d: FOPT −23 %, FGPT +28 %, FWPT −41 %, FPR +8 % (too much gas, too little oil). Norne uses physics the engine does not have yet: **EHYSTR hysteresis, SCALECRS 3-point scaling, VAPPARS, GCONINJE groups and WCONPROD NONE wells (which run on BHP)**. These are the next parity items. |
| Lenient deck import | Decks with features the own engine cannot model (gas-water systems, PVZG, RSCONST, lab/radial grids, keywords that fail to parse, ...) are no longer rejected. They import with placeholders, so the grid, wells, schedule and deck text can be viewed and edited. The reasons are stored in `Case::engineUnsupported` and shown in the toolbar ("OPM Flow / ECLIPSE deck" tooltip). OPM Flow and ECLIPSE always run the original file; ResSim gives a clear message instead. When a case lists no wells, well names come from the summary file. |
| Editor | Deck checker: a section header with trailing `====` no longer breaks the next keyword. Tab close crosses were invisible below 100 % UI size, because ImGui draws them with thickness `(int)_MainScale` = 0. They are now always drawn. |
| Tests | 152 native checks pass (14 new numerics checks), and 152 under Wine. |

**Still open (honest list):** Norne parity physics (EHYSTR, SCALECRS, VAPPARS), groups/guide rates (B12), VFP/THP (B13), MSW (B14), WECON/WTEST/ACTIONX/UDQ (B15/B16), RESTART (C19), the binary output writer (C18), EDITNNC/PINCHREG, hysteresis, a parallel/multicolour BILU0 and the GPU backend. The BSR/LinearOperator layer is the basis for the GPU backend: the same block kernels map onto CUDA/OpenCL/D3D compute.

## Round: Norne parity, groups, VFP/THP, real ECLIPSE example decks

| Item | What was done | Verified by |
|---|---|---|
| Real ECLIPSE example decks | BRILLIG (9239 lines, old fixed format), ALKALINE_ASP and FOAM import and check cleanly. Parser: a column-1 keyword line with trailing text is accepted only when every record is closed; `,` separates items; padded quoted names are trimmed; a list may end at the next keyword. Physics options the engine does not have (ASP, FOAM, ...) are reported and kept as text, so pass-through runs are unchanged. | Deck checker: 0 errors / 0 warnings on all three decks and on Norne. BRILLIG export → re-import gives the same results; ENDSCALE table counts (NTENDP/NSENDP) are now kept for ENPTVD. |
| **EQUALS / COPY records that start with an array name** | Inside an open `EQUALS` list, the line `SWL 0.05 1 46 1 112 2 2 /` was read as a new keyword `SWL`. Only the first record was applied. On Norne this lost SWL in 21 of 22 layers, and all of SWCR/SGU (COPY), SGCR, SOWCR, SOGCR and SWU. The parser and the editor checker now keep such lines as records while a "records until an empty /" list is open. | All 8 end-point arrays and ISGCR/ISOGCR now match the ECLIPSE INIT file cell by cell (max difference 3e-8). New regression test. |
| **PINCH across ACTNUM = 0 layers** | With `PINCH … GAP`, the grid bridged thick layers switched off by ACTNUM (Norne layer 4 and others). That created 2263 vertical connections that ECLIPSE does not have (Σ T ≈ 2·10⁵, about 13× all of Norne's fault NNCs) and opened the shale barriers. Now only cells removed by MINPV / zero pore volume (with GAP), or cells thinner than the threshold, are bridged. | TRANX/TRANY sums within 0.2 % of ECLIPSE and TRANZ exact. All 13455 ECLIPSE NNCs are present (Σ T +1.3 %). MINPV bridges match ECLIPSE TRANZ. Test updated. |
| Hysteresis | SATOPTS HYSTER, EHYSTR (model, curvature), NOHYST, IMBNUM, ISWL…ISOGCR. Carlson for the non-wetting phases (gas; oil in the oil-water system); Killough is approximated by Carlson; wetting phases follow drainage. | Gascap test: less gas produced with trapping. Deck round trip. |
| SCALECRS | Three-point scaling (krw at 1 − SOWCR − SGL, krg at 1 − SOGCR − SWL, kro at the displacing phase's critical saturation). | Test: the kr at the scaled residual point equals the table value. |
| VAPPARS | Rs/Rv limited by (So/So,max)^v. | Condensate test, deck round trip. |
| Groups | GRUPTREE; GCONPROD ORAT/WRAT/GRAT/LRAT/RESV targets and limits (GCONPRI is read as limits); WCONPROD GRUP; GCONINJE RATE/RESV/REIN/VREP. Shares are proportional to the well potentials. Water filling caps each well at its own target and potential and gives the rest to the other wells. Exceed checks use potentials, not the previous step's rates, so wells do not flip on and off. | SPE9LIKE FIELD ORAT 22500: FOPR 22500 (uncapped 37500). BRILLIG GCONPRI 100000 stb/d: year-1 average 99970. |
| VFP / THP | VFPPROD/VFPINJ with 5-D multilinear interpolation and extrapolation. The THP limit is solved as the intersection of the IPR and the VFP curve by bisection on BHP, with a head correction to the table datum. The earlier lagged fixed point locked wells at impossible BHPs. | Unit tests (exact for linear tables; THP producer). BRILLIG runs with VFPPROD 1/2 and VFPINJ 1. |
| Tests | 165 native checks pass, and the same under Wine. New: VFP interpolation, SCALECRS, hysteresis round trip and run, VAPPARS, group target, THP, EQUALS records, PINCH/ACTNUM. | |

**Norne (own engine) vs ECLIPSE 2014.2:** The full history runs in 28 min on 2 threads (3070 Newton iterations, 23 cuts). At 3312 d: **FOPT 6.611e7 vs 6.732e7 (−1.8 %, was −23.7 %)**, FGPT 1.556e10 vs 1.464e10 (+6.3 %, was +21.6 %), FWPT 1.619e7 vs 2.320e7 (−30 %), FPR 303.4 vs 271.6 bar. Early time is close: at 201 d FOPT −0.7 % and FGPT +1 %; at 817 d FOPT −3.6 % and FGPT +2.5 %. Most of the gain came from the two input bugs above (the EQUALS end points and the pinch-outs across ACTNUM barriers), not from new physics. **Remaining gap:** the field pressure is 6–12 bar high from about 300 d. All well BHPs are higher than ECLIPSE (producers by 25–45 bar, injectors by 30–95 bar), and water production is low. After 2700 d the pressure climbs (279.6 → 303.4 bar) while ECLIPSE stays at about 273. Next checks: well-block pressure (WPAVE/WBP), injector connection mobility, the RESV conversion at region-average pressure, and true Killough hysteresis.

**Still open (honest list):** true Killough hysteresis and Pc hysteresis; guide-rate formulas, priority allocation, GCONSALE/GCONSUMP, networks and gas lift; VFP curve viewer; MSW (B14); WECON/WTEST/ACTIONX/UDQ (B15/B16); RESTART (C19); the binary output writer (C18); EDITNNC/PINCHREG; AQUTAB/numerical aquifers; and the GPU backend.

## Round: Norne physics parity II (MULTREGT, re-opened connections, RESV conversion, three-point krw)

Method: every ECLIPSE output file in the Norne folder was used as a reference, not only the summary. That means INIT TRAN*/TRANNNC per FIPNUM region pair, the PRT well-connection report (BLOCK rows), and the RFT file (cell P/Sw/Sg). Each gap was traced to a single cause before anything was changed. New tools in `tools/`: `norne_diag` (region RPR/ROIP + well table at several times), `norne_conn` (per-connection rates vs the PRT report), `rft_cmp` (cell state vs RFT), `norne_regtran` (FIPNUM-pair transmissibility vs INIT/EGRID), `norne_tran`, `norne_resv` (RESV conversion vs FVPR), `compare_csv`, `rft_dump`.

| Item | What was wrong / what was done | Verified by |
|---|---|---|
| **MULTREGT on fault NNCs** | MULTREGT was converted to MULTX/Y/Z, so it never reached non-neighbour connections. 3158 units of fault-NNC transmissibility across `m = 0` region pairs stayed open, and regions 5–15 were about 20 bar too high by 450 d. MULTREGT is now applied in the grid to every connection: last record wins; region order is symmetric; defaulted region = any; directions XYZ; NNC/NONNC/ALL; region set M/F/O (MULTNUM/FLUXNUM/OPERNUM). Connections with T = 0 are removed. The exporter writes MULTREGT plus the region arrays. | ECLIPSE TRANNNC across m = 0 pairs: 0 = ours 0. Region RPRs at 452 d are within 1 bar (except region 13). 5 new tests (ALL, NNC with reversed order, NONNC = 0, direction filter, deck round trip). |
| **COMPDAT / WELOPEN shut → re-open** | Each connection had only one [from, until) window, so a connection that was shut and later re-opened was open the whole time (15 Norne wells). A re-open now adds a new open period. Later COMPDAT/WELOPEN records act on the latest period. New connection data (CF, Kh, skin, …) for an open connection starts a new period instead of overwriting the past retroactively (C-2H CF 21.685 → 22.834 was applied from 441 d instead of 616 d). | Test: a producer shut at step 10 and re-opened at 20 produces 0 at step 15 and > 0 at 25; deck round trip is stable. |
| **WCONHIST RESV conversion (ECLIPSE/OPM)** | We converted the observed rates at the mean pressure of the connected cells and controlled the in-situ voidage. ECLIPSE/OPM (RateConverter) convert both the target and the well's own surface rates with **field hydrocarbon-pore-volume-averaged p, Rs, Rv** (Rs ≤ qg/qo, Rv ≤ qo/qg), updated every time step. Implemented: `updateResvConversion()` gives coefficients resv = cW·qw + cO·qo + cG·qg, used in the control equation, the potentials and the reported RESV. | `norne_resv`: ECLIPSE FVPR = history rates converted at FPR to within 1 %. Norne 8 d: FOPT 35048 vs 35047. |
| **SCALECRS three-point krw with SOWCR_t = 0** | The Norne SWOF has kro = 0 only at Sw = 1, so the middle point equals the maximum (SR_t = SWU_t). Our code then fell back to two-point scaling, which made water mobility up to 30 % too low. It now matches OPM `EclEpsTwoPhaseLaw::scaledToUnscaledSatThreePoint_` exactly: [SWCR, SR] maps onto [SWCR_t, SR_t], the result is clamped, and there is no fallback. | D-1H at 8 d: water 4.51 vs 4.507 sm3/d, BHP 219.6/219.6. Per-connection water at 424 d had been 15–40 % low on every connection. New unit test. |
| WCONHIST / WCONINJH observed BHP | Item 10 (WCONHIST) / item 5 (WCONINJH) is history data, but it was used as a BHP limit. It is now kept in `bhpHist`; the limit stays 1 atm / 100000 psi unless WELTARG sets one, and BHP mode uses it as the target. A history-well BHP limit is exported as WELTARG. | Test with round trip. |
| Importer messages | No more "re-opened … treated as open throughout"; a SHUT `WCONPROD 'B-4H' 'SHUT' 15*` no longer produces the misleading "control mode NONE" warning. | `rundeck -steps 0` |
| API | `Simulator::connectionRates()`: per-connection surface rates, block and connection pressure, Sw/Sg (the same data as the ECLIPSE PRT well-connection report). It can feed a connection-level results table in the app. | `norne_conn` |
| Tests | **174** native checks pass, and 174 under Wine. All 15 examples give identical FOPT (no RESV, no MULTREGT). BRILLIG is unchanged. | |

**Cell-level agreement (RFT, own engine vs ECLIPSE):** B-4H 143 d, D-4H 201 d, C-1H 229 d, E-3H 297 d: every connection cell has P within 0.1 bar and Sw/Sg identical to 3 decimals. C-3H 544 d: P +1.4 bar (layers 5–18) and +4.3 bar in the water layers 19–20 of FIPNUM 13 (C-2H water injection). F-1H 560 d: +0.6 bar.

**Norne field (own engine, 40 min, 2 threads, 3323 Newton, 35 cuts) vs ECLIPSE:**

| t (d) | FPR | FOPT | FGPT | FWPT |
|---|---|---|---|---|
| 414 | 241.6 / 241.5 | +0.06 % | +0.3 % | +2.0 % |
| 694 | 242.8 / 239.7 | −0.5 % | +0.8 % | −0.1 % |
| 1039 | 252.8 / 244.7 | −0.6 % | +0.2 % | +4.6 % |
| 1913 | 270.1 / 257.5 | −3.0 % | +1.1 % | +18 % |
| 3312 | 308.2 / 271.6 | **−4.4 %** (was −1.8 % for the wrong reasons) | **+0.4 %** (was +6.3 %) | **+9.5 %** (was −30 %) |

**Remaining gap (diagnosed, not yet fixed):** from about 550 d the connected main field (FIPNUM 5–15) becomes uniformly 8–10 bar high, while regions 1–4 stay exact. ROIP per region matches (region 13: 4.3507e6 vs 4.352e6), and FIPNUM-pair transmissibility matches within about 3 %. The excess starts in the FIPNUM 13 water layers around the C-2H/C-1H/C-4H injectors (RFT C-3H +4.3 bar at 544 d). It then grows by feedback: a higher FPR gives a smaller RESV voidage target. After 3000 d ECLIPSE injectors C-1H/F-4H reach the 600 bar limit (ECL FWIT is 1.8 % below history) while ours keep injecting. Next: injector connection-level comparison (no injection report in the PRT, so use WBHP/WWIR per injector), the WPAVE-defined WBP, the pore-volume pressure function (ROCK, ECLIPSE quadratic), and THPRES on the region-13 water layers.

## Round: LGR, full dual porosity / permeability, compositional plan, thermodynamics library (C1) + compositional FIM core (C2), Norne parity III

**Gridding: local grid refinement (done)**
- CARFIN / REFINE / ENDFIN boxes with NXFIN/NYFIN/NZFIN and HXFIN/HYFIN/HZFIN (non-uniform splits), RUNSPEC LGR dims, and
  property arrays inside CARFIN (PORO/PERM*/NTG/ACTNUM/MULT*). Import and export round-trip; multiple LGRs.
- Cartesian and corner-point hosts. Local cells inherit the host geometry (sub-cells of the host hexahedron, so faults,
  MULTFLT and NNCs are kept). Fine–coarse connections are the redistributed host connection: never a direct TPFA between
  local cells of different hosts, which would lose fault throw / MULTFLT.
- Wells: WELSPECL / COMPDATL (local ijk), aquifers on refined hosts, streamlines with multi-link faces.
- UI: LGR editor (boxes, splits, per-LGR arrays), well-editor grid combo, viewer glyphs/probes, wells-map overlay.
  Example `examples/LGR.DATA`.
- Verified: a trivial (1×1×1) LGR reproduces the parent run exactly (T 73659.30, FOPT 7305000.6). A full 2×2×2 refinement
  matches a natively fine grid (FOPT 7304999.67 vs 7305000.60). A faulted corner-point grid keeps its connection count
  and pore volume.
- WELOPENL (local-ijk WELOPEN for LGR wells) ✅ added in the C2 round.
- Gaps: LGR + dual porosity together, ECL flux mapping for LGR cells,
  amalgamation (AMALGAM), auto-refinement.

**Dual porosity / permeability (completed)**
- New: GRAVDR gravity drainage (Gilman–Kazemi matrix–fracture potential with the height DZMTRX / DZMTRXV, oil–water and
  gas–oil). Verified against the hydrostatic equilibrium in a single block: Sw matrix 0.3196 vs fracture 0.3216 against
  the analytical 0.32, water balance 1.8e-7.
- New: DPNUM (single-porosity regions inside a dual-porosity model). The matrix cell is folded into the fracture cell,
  with no M–F connection in those columns. Active cells 1350→1125, oil balance 1e-8.
- Already present: DUALPORO / DUALPERM, SIGMA / SIGMAV / LTOSIGMA, LX/LY/LZ shape factors, NODPPM, fracture
  permeability × porosity.
- Not supported (kept for OPM/ECLIPSE, with a warning): GRAVDRM (alternative gravity drainage), GRAVDRB, the DPGRID
  option. The fracture always takes the matrix geometry, which is what DPGRID would give anyway.

**Compositional / thermal / CO2STORE: plan and order (decision)**

Compositional comes first. CO2STORE is a special case of compositional (CO2–H2O–NaCl with an aqueous-phase solubility
model), so it is the C4 milestone. Thermal (T1) adds an energy equation to the compositional core; it reuses the EOS
residual enthalpy that C1 already provides.

| Milestone | Content | Verification |
|---|---|---|
| **C1 ✅ (this round)** thermodynamics library `eos.h/.cpp`, `eos_deck.cpp` | PR76 / PR78 (PRCORR) / SRK, BIC, Peneloux SSHIFT, per-component OMEGAA/B. Michelsen–Mollerup reduced-Helmholtz formulation with **analytic** ∂lnφ/∂n, ∂P, ∂T, partial molar volumes, residual enthalpy. Min-Gibbs root. Michelsen stability (Wilson vapour/liquid trials, SS, then second-order in α = 2√W with Levenberg damping). Rachford–Rice (Whitson–Michelsen window, safeguarded Newton). SS + GDEM, then Newton flash on vapour moles. K warm start for the simulator. Saturation pressure (Newton in ln P, with a stability-bisection fallback for near-critical fluids). Li phase label. LBC viscosity (LBCCOEF, VCRITVIS/ZCRITVIS). E300 keywords NCOMPS/COMPS, EOS, PRCORR, CNAMES, TCRIT, PCRIT, VCRIT, ZCRIT, ACF, MW, BIC, SSHIFT, PARACHOR, VCRITVIS, ZCRITVIS, LBCCOEF, OMEGAA/B, RTEMP/RTEMPA, TEMPVD, ZMFVD, ZI, STCOND, CO2STORE. The import report shows the fluid and its saturation pressure; decks still run in E300/OPM | SPE5 oil Pb **2302.89 psia vs 2302.9 published**. Fugacity equality 3e-13, mass balance 4e-18. FD derivative errors ≤ 2e-7, Gibbs–Duhem 3e-14. C3 Psat(300 K) 9.977 bar (NIST 9.97). CO2 Psat(280 K) 41.59 bar (41.6). CH4 μ 0.0112 cP. Retrograde dew point consistent with flash. SRK |
| **C2 🟡 core done (C2 round)** `compsim.h/.cpp` | E300-style natural-composition FIM: primary variables p, N_1..N_nc (moles per cell), equations = nc component balances + the volume balance Σ V_phase = PV(p) (Acs–Doleschall–Farkas / Coats). One flash per cell per Newton iteration with the previous-K warm start and stability only when the cell was single-phase or RR leaves (0,1). All phase properties (ξ, x/y, vm, ρ, μ_LBC, kr) carry exact derivatives in a dual-number type built from `flashSensitivities` (implicit-function theorem on the fugacity equations), ∂vm/∂x, and the new analytic LBC derivative overload. TPFA with gravity, phase-potential upwinding, Corey kr on the HC pore space with immobile Swc, rock compressibility. Mobile immiscible water as an optional third phase (block size nc+2): PVTW-style Bw(p), Corey krw / krg / krow / krog and the ECLIPSE default three-phase oil model kro = (Sg·krog + (Sw−Swc)·krow)/(Sg + Sw − Swc). Wells: molar-rate injector with any injection composition, water injector; producers on BHP, hydrocarbon molar rate or reservoir-volume (RESV) rate with a BHP limit. Under rate control the BHP is eliminated analytically (exact Schur complement of the well equation, so perforated cells are coupled in the Jacobian), with an in-Newton switch to BHP and no crossflow. Linear solver: existing BSR (kMaxBlock raised to 16, so nc ≤ 15) + BlockCPR + GMRES; damped Newton (|Δp| ≤ 20 bar, N ≥ 0) and Δt cuts. Convergence: per-cell 1e-7 plus a global material-balance test 1e-12 per component (E300 TOLMB style). Next: well hydrostatic head between connections, injector BHP limits, capillary pressure, near-critical kr (MISCIBLE), deck hook-up (C3) | Analytic Jacobian vs central differences on a 3-cell vertical column with single- + two-phase cells, gravity and wells: **max scaled error 2.7e-9**. Slim tube (40 cells, SPE5 oil, SPE5 injection gas, 1.2 HCPV): component balance **2e-12**; C10+ recovery **97.7 % at 4500 psia** (above the 3875 psia FCM pressure) vs **83.1 % at 2600 psia**; ~3.7 Newton iterations per step. Three-phase Jacobian (mobile water, water + gas injectors) vs FD: **4.7e-9**. Waterflood of the SPE5 oil (100 cells): breakthrough **0.515 PV vs 0.543 PV Buckley–Leverett** (Welge tangent, μo 0.23 cP from LBC); the early arrival is the usual first-order upwind dispersion. Two-connection rate-controlled producer Jacobian vs FD: **5.3e-9**. A RESV voidage-replacement producer holds the mean pressure at 3500.3 psia (initial 3500) |
| C2 (remaining) ✅ mostly done (C3 round) | ✅ SWFN/SGFN/SOF3 tables, ✅ Stone II, ✅ capillary pressure (Pcow, Pcgo), ✅ implicit injectors with BHP limit, ✅ hydrostatic heads. Still open: near-critical kr interpolation (E300 MISCIBLE / Coats), Baker kro, natural-variable option | 1-D slim tube vs analytic MOC solution (Orr) for a 3-component system |
| C3 wells and surface 🟡 deck hook-up done (C3 round) | ✅ Well model with mixture mobility, ✅ FIELDSEP separator train to stock tank, ✅ WCONPROD/WCONINJE with WELLSTRE/WINJGAS injected composition, ✅ field/well summary, ✅ EQUIL + ZMFVD initialisation, ✅ app engine. Still open: SEPCOND/WSEPCOND per-well separators, group controls, secondary rate limits, XMF/YMF/FMWPR vectors, SPE3 deck | **SPE5** (WAG, 6-component) ✅ runs 20 years; **SPE3** (gas cycling, 9-component) next; E300 reference results still needed |
| C4 CO2STORE | CO2–H2O–NaCl. Spycher–Pruess (2005/2010) mutual solubilities, Duan–Sun / Ezrokhi brine density, Garcia apparent CO2 volume, PR gas with Peneloux, Fenghour CO2 viscosity, CO2SOL, SALINITY / SALT. Also the H2STORE variant | OPM CO2STORE test decks; SPE11 type-A geometry comparison |
| T1 thermal | Energy equation (fluids via EOS departure + ideal-gas Cp polynomials, rock HEATCR/HEATCRT, THCONR conduction, HEATVAP), heat loss to over/underburden (Vinsome–Westerveld), temperature as primary variable. Steam tables / K-value thermal (E300 THERMAL) after | Analytical conduction (Lauwerier) and the SPE4 thermal problem |

**Norne parity III (this round's Norne work)**
- Injector check (ours vs ECLIPSE 2014.2, cumulative WWIT / WGIT at our report times): C-1H, C-4H and all other injectors
  match essentially exactly through 3092 d. The earlier per-well "mismatch" was a time-alignment artefact.
- **Fix 1: DRSDT / VAPPARS precedence.** Norne enters `DRSDT 0` and then, in the same step, `VAPPARS 0.5 0`. In OPM
  (OilVaporizationProperties: the last keyword sets the type, and drsdtActive() requires type DRDT) and ECLIPSE the two
  exclude each other, so VAPPARS switches the DRSDT limit off. We applied both, so gas could never re-dissolve. That left
  too much free gas and pressure ran high from ~450 d. The importer, engine and exporter now keep the order (last
  keyword wins). FPR at 737 d: 244.53 → **240.30 (ECL 240.34)**. At 2266 d FPR is 266.1 / 266.5, FOPT −0.2 %, FGPT −0.2 %.
- **Fix 2: group injection with stopped injectors.** From 2947 d `GCONINJE FIELD WATER RATE 42000` is active. In 2006
  the C-injectors are stopped (WELOPEN STOP + WCONINJE RATE 0). Our group allocation took an explicit target of 0 as
  "no limit", so a stopped injector received the group's share: C-2H took 21 900 sm3/d instead of 0, an extra
  ~4.5e6 sm3 of water. An explicit RATE 0 now caps the well at 0 for the group sum and for the water-filling share.
  Regression test: before the fix the zero-target injector took all 900 and the real injector 0.
- Before these fixes (last round), ours/ECL at 3312 d: FPR 308.2 / 271.6, FOPT −4.4 %, FWPT +9.5 %. After fix 1 only:
  FOPT +0.18 %, FGPT +0.6 %, FWPT +3.0 %, FPR within 0.1–0.5 bar to 2266 d and 2–5 bar to 3092 d. After both fixes, ours/ECL:

| t (d) | FPR (bar) | FOPT | FGPT | FWPT |
|---|---|---|---|---|
| 737 | 240.30 / 240.34 | +0.7 % | | |
| 2266 | 266.08 / 266.52 | −0.2 % | −0.2 % | +2.6 % |
| 3092 | 277.37 / 283.24 | +0.15 % | +0.9 % | +2.3 % |
| 3312 | **265.46 / 271.59** (was 308.2) | **+0.06 %** (was −4.4 %) | +0.9 % | **+2.5 %** (was +9.5 %) |

  Field in place at 3312 d: FOIP −0.04 %, FWIP −0.15 %, FGIP −0.6 %. WWIT / WGIT of every injector match ECL exactly,
  except C-1H water (−1.5 %), which reaches its 600 bar BHP limit (593 bar at 3092 d) a little earlier than in ECL.
- Remaining Norne gap: from ~2500 d FPR runs 2–6 bar below ECL, with gas production +0.9 % and water +2.5 %. Next
  suspects: C-1H injectivity near the BHP limit and gas-cap / GOC movement in the late gas-injection period.

**Gap audit additions this round**: ~~connectionFluxes (streamline flux) lacks the GRAVDR term~~ ✅ fixed in the C2 round (verified: flux = accumulation). DUALPERM wells on DPNUM-folded
matrix cells are untested. The DPNUM importer ignores BOX.

## Round: compositional C2 completion + C3 deck hook-up (E300 decks run on the own engine)

**Core (`compsim.h/.cpp`, rewritten).**
- **General geometry.** The core takes either the Cartesian box or a cell list with PV, depth and a connection list {a, b, T}. The deck path passes the processed black-oil `Grid` in SI units, so corner point, faults, NNC, MULTREGT, PINCH, LGR and dual porosity all carry over.
- **Saturation functions and capillary pressure.**
  - SWFN / SGFN / SOF3 tables (piecewise linear, exact slopes), or Corey when no tables are given.
  - Three-phase oil kr: ECLIPSE default model or Stone II.
  - Pcow and Pcgo enter the phase potentials: p_w = p − Pcow(Sw), p_g = p + Pcgo(Sg), with oil pressure as the primary unknown. Pc also enters the well connection potentials.
- **Unified implicit well model (producers and injectors).**
  - Every flowing term has outflow B_i·(Φ_t − p_bh). Φ_t = p_k + Pc_α − H_k, where H_k is the wellbore head ρ_well·g·(D_k − D_ref), lagged one step (mobility-weighted mixture density for producers, injected-fluid density for injectors).
  - **Injectors are implicit:** total-mobility injectivity WI·Σkr/μ·ξ_inj(p). They have a BHP limit with an in-Newton switch, and no crossflow.
  - **Rate control:** BHP is eliminated analytically, with dBHP/du_j = Σ_{t∈j} (dΛ_t·(Φ_t − BHP) + Λ_t·dΦ_t)/ΣΛ. Potentials are taken relative to the first connection, which removes the cancellation at large WI.
  - **Producer targets:** BHP, HC molar, RESV, ORAT, GRAT, WRAT, LRAT.
  - **Surface rates** use lagged per-component separator yields f_o,i and f_g,i (E300 practice). The rate weight Σ_i x_iα·f_o,i is linear in the phase compositions, so composition and mobility changes within a step are implicit. ORAT overshoot at the bubble point is 0.13% (2.8% with a single mixed-stream factor).
- **Separator train:** stages (T, p), with the last liquid flashed to STCOND.
- **Reporting:** surface rates come from flashing the produced stream; per-well reports and cumulatives; FPR (HCPV-weighted); FOIP / FGIP / FWIP through the separator.

**Deck hook-up (`comp_deck.h/.cpp`, `tools/runcomp.cpp`).**
- **Sources:** E300 PROPS via `eos::fluidFromDeck`; grid, wells and schedule from the black-oil importer.
- **Mapping:**
  - SWFN/SGFN/SOF3 are rebuilt exactly from the importer's SWOF/SGOF (union of oil-saturation nodes).
  - WCONPROD / WCONINJE / WELOPEN controls per report step; COMPDAT intervals with Peaceman CF (same formula as the black-oil engine).
  - WELLSTRE + WINJGAS STREAM compositions are tracked per report step.
  - FIELDSEP stages; ROCK, PVTW, DENSITY.
- **EQUIL + ZMFVD/ZI initialisation:**
  - Oil pressure: RK2 integration of the EOS mixture density of z(depth) from the datum.
  - Water pressure: integrated from the OWC.
  - Sw from inverting Pcow; water-zone cells take p_o = p_w + Pcow(Sw).
- **Notes:** unsupported items (group controls, secondary limits, SEPCOND, hysteresis, aquifers, multiple EOS regions) are listed in the run notes instead of being dropped silently.
- **Two-way deck:** E300 keywords are kept verbatim; export → re-import rebuilds an identical model.

**App and importer.**
- The ResSim engine now runs COMPS decks on the compositional core: original deck for the fluid, (edited) case for grid, wells and schedule. Summary, well rates and BHP, pressure / Sw / Sg / So snapshots go to the Results / 3D views.
- The import report says COMPS decks run on the compositional core. E300 PROPS keywords show as supported in the editor checker. The black-oil "PVDO missing" noise is gone for COMPS decks.

**Verification (`-comp`).**
- **Jacobian vs central differences (max scaled error):**

| Case | Error |
|---|---|
| Immobile water | 2.7e-9 |
| Mobile water | 4.7e-9 |
| Rate producer | 5.3e-9 |
| General connections (NNC) + SWFN/SGFN/SOF3 + Pc + implicit two-connection injector + ORAT producer with heads | 5.9e-8 |
| Stone II + BHP-limited injector + LRAT producer | 8.0e-8 |

- **Run test:** the injector switches to its 3300 psia limit (|ΔBHP| < 1 Pa); the ORAT producer holds 2.0000 sm3/d through a 2-stage separator.
- **Equilibrium:** 24-cell column with a Pc transition zone across the OWC (EQUIL + ZMFVD, METRIC). Datum pressure 280.0000 bar; after 2 years without wells, max |Δp| = 1.5e-5 bar and max |ΔSw| = 1.3e-5.
- **SPE5 compositional WAG deck** (`decks/spe5/SPE5_E300.DATA`, authored: SPE5 grid/SCAL/wells + the 6-component PR78 fluid, WELLSTRE injection gas):
  - Separator GOR above Pb is 0.556 Mscf/stb (black-oil table Rs 0.573).
  - ORAT 12000 is held until ~150 d; FWIR / FGIR are exactly 12000.
  - The full 20-year run takes 11 s (733 steps, 5 cuts).
  - Water balance is exact: ΔFWIP = FWIT − FWPT = 22.45 MMSTB.
  - FOPT 24.28 MMSTB of 35.42 MMSTB OOIP.
  - Comparing against E300 / published SPE5 curves is still open.
- **Waterflood / slim tube:**
  - Waterflood breakthrough 0.525 PV vs 0.543 PV Buckley–Leverett.
  - Slim tube C10+ recovery 98.1% at 4500 psia vs 83.5% at 2600 psia.
  - Balance 3e-11.

## Norne parity IV: well-control fixes (C3 round)

The late-time FPR gap (−5 to −6 bar after 2500 d) came from **well control**, not reservoir physics. Two bugs hid each other.

1. **Stale 1 atm BHP and a wrong stall lock.**
   - A producer at its 1 atm BHP limit, stopped by `WCONHIST … ORAT 0 0 0` and then given a small RESV target (Norne B-1H, Aug–Oct 2004), restarted Newton from BHP = 1 atm.
   - The BHP update limit 0.3·max(BHP, 10 atm) needs about 15 iterations to climb to about 220 bar. Before that, the "rate target cannot be met" stall rule locked the well at 1 atm for the whole report step. Result: 10–40× the history rate (936 stb/d instead of 74).
   - Fixes:
     - The BHP limiter is now relative to max(BHP, cell pressure).
     - The stall rule never locks a well that is producing *more* than its target (the target is deliverable).
     - A well stopped by a zero target restarts from 0.9·p_cell, like a newly opened well.
2. **WCONHIST observed rates were treated as limits.**
   - When switching back from BHP control, every observed ORAT/WRAT/GRAT/LRAT was a candidate target. Norne D-1H at 1091 d went to WRAT 270 stb/d instead of producing 4478 stb/d at 1 atm, as ECLIPSE does.
   - Previously this bug was hidden by bug 1's lock. Now history wells return only to their WCONHIST control mode (ECLIPSE / OPM semantics).

| t (d) | FPR ours / ECL | FOPT | FWPT | FGPT |
|---|---|---|---|---|
| 737 | 240.30 / 240.41 | +0.2 % | −1.5 % | +0.1 % |
| 2266 | 266.15 / 266.52 | −0.18 % | +2.6 % | −0.19 % |
| 2917 | 280.10 / 280.96 (was 275.57) | −0.13 % | +0.8 % | −0.12 % |
| 3092 | 282.32 / 283.24 (was 277.37) | −0.08 % | +0.5 % | −0.10 % |
| 3312 | **270.86 / 271.59** (was 265.46) | **−0.08 %** | **+0.37 %** (was +2.5 %) | **−0.08 %** (was +0.9 %) |

- **Run:** 986 s, 3006 Newton, 26 cuts (was 1100 s, 3588, 38).
- **Tests:** "a producer reopened after a zero target honours its new rate target" and "WCONHIST: observed rates are not limits".
- **Remaining Norne items:**
  - FWPT +2.6 % around 2266 d, which recovers later.
  - C-1H injectivity is about 4 % low while on its 600 bar limit.
  - WTEST (D-3BH, K-3H) is not modelled.

## Round C4: dedicated CO2STORE module + full compositional / CO2 UI

**Engine: `co2brine` (new module, engine/src/co2brine.cpp)**
- CO2–brine phase equilibrium: Spycher–Pruess (2005) Redlich–Kwong mutual solubility (b = 27.80 / 18.18 cm³/mol) with Duan–Sun salting-out (NaCl molality), Batzle–Wang brine density / viscosity, Garcia (2001) apparent molar volume of dissolved CO2 (dissolution makes the brine *denser*, which drives convective mixing), Fenghour (1998) CO2 viscosity.
- Checks at 50 °C, pure water: x_CO2 at 200 bar 0.0229 (experiment 0.0223–0.0239); ρ_CO2 764 kg/m³ at 200 bar (NIST 784), 864 at 300 bar (NIST 871); 1 m / 4 m NaCl lowers m_CO2 from 1.30 to 1.05 / 0.62.
- Compositional FIM core (`CompSimulator`) runs CO2STORE cases with two components (CO2, H2O), using the co2brine flash in place of the cubic EOS. The inventory is split into mobile (Sg > Sgc), residually trapped (Sg ≤ Sgc) and dissolved CO2 (`CompSummary::co2Mob/co2Immob/co2Diss`; FGMMO / FGMTR / FGMDS / FGMIP in the app).
- Deck: OPM CO2STORE format in both forms. `OIL GAS CO2STORE` (OIL = brine, SGOF) and `GAS WATER CO2STORE` (SWOF/SGOF or SWFN/SGFN). Also SALINITY (mol/kg), RTEMPVD/RTEMP, EQUIL, RSVD. Import → export → re-run gives identical results (test).
- `validateCompositional()`: compositional cases are validated on their own (EOS data, streams, separator, ZI/ZMFVD); the black-oil PVT checks are skipped.

**App: full UI**
- **Fluid section** has a model selector: Black oil / Compositional (EOS) / CO2 storage (CO2STORE). Switching converts the case (CO2STORE fluid at the current temperature / salinity, or the SPE5 fluid as the compositional starting point).
- **Compositional editor** tabs:
  - Components: table editing Tc/Pc/ω/MW/Zc/Vc/shift/parachor, a component library, the SPE5 preset, and a plus-fraction pseudo-component.
  - Binary interaction: matrix editor plus Chueh–Prausnitz generation.
  - Initial composition: ZI or ZMFVD versus depth.
  - Injection & separator: WELLSTRE streams, WINJGAS per well and step, FIELDSEP stages.
  - Phase behaviour: P–T phase envelope at the current composition, an interactive flash (phase split, x/y, K-value plot), and LBC viscosities.
  - Deck keywords: live preview of the E300 keywords that will be written.
- **CO2 storage editor** tabs:
  - Fluid & initial state: temperature, salinity in mol/kg, ppm or wt %, brine phase OIL ↔ WATER (with SGOF ↔ SWOF conversion), contact and datum, and a property table.
  - Property plots (2×2): solubility, CO2 density, viscosities, brine density versus pressure.
  - Deck keywords.
- **Results:** new "CO2 storage" tab with the trapping inventory (stacked mobile / trapped / dissolved) versus time, and an XCO2 cell property: brine CO2 mole fraction for CO2STORE, overall z₁ for EOS cases. XCO2 is also in the Results 3D viewer property list. FGMMO/FGMTR/FGMDS/FGMIP are plot vectors.
- **Initial state preview** (Initialization → compute) uses the compositional / CO2–brine equilibration for these cases.
- **Run:** compositional / CO2STORE cases built in the UI or imported run on the own compositional core. The status bar shows "Compositional / CO2 storage model valid".
- **Examples:**
  - `co2`: 40×1×20 dipping aquifer, 300 t/d for 10 y, then 40 y of migration. 11 s run; in-place CO2 equals injected CO2 (1.09 Mt).
  - `comp`: SPE5 fluid, lean-gas WINJGAS after one year of depletion. 3 s.
  - Both are exported as examples/CO2.DATA and COMP.DATA, alongside the OPM-format examples/CO2STORE.DATA.
- **Tests:** 242 checks pass. New checks: CO2 thermodynamics, CO2STORE deck run and two-way export, the CO2 example's mass balance and > 10 % dissolution, and the compositional example's rate.

**Known limitations (to do)**
- RK CO2 density error is ±3–7 % between 80 and 200 bar (near-critical). Next: Span–Wagner or a tabulated EOS.
- Spycher–Pruess (2005) is validated only up to 100 °C. Above that, the 2010 Spycher–Pruess high-temperature model is needed.
- No gas-saturation hysteresis, so residual trapping only counts Sg ≤ Sgc cells. Next: Land / Carlson trapping (EHYSTR) in the compositional core.
- Isothermal only: RTEMPVD takes the mean temperature. No thermal or Joule–Thomson effects at the injector.
- No CO2 in the oil phase (CO2-EOR) in CO2STORE mode; use the EOS compositional mode for that.
- No H2STORE / gas-water (GASWAT) mode or diffusion (DIFFC) yet.
- Studies (sensitivity / HM) are not yet available for compositional cases; they wrap the black-oil engine.
- The .rsim project format has no compositional support; the deck is the save format.

## Round C5: CO2 trapping mechanisms + their interactions, compositional completion

**Engine (compsim / co2brine / comp_deck): every trapping mechanism, coupled in one fully implicit step**
- **Structural / capillary trapping:** SATNUM regions in the compositional core (caprock tables with their own entry pressure); the inventory reports *structural* CO2 = mobile gas at rest under a seal (outflow below a tolerance, RESSIM `STRUCTTOL`).
- **Residual (hysteresis) trapping:** gas-phase hysteresis in both CO2STORE and EOS modes. EHYSTR model 0/1 → Carlson, 2/3/4 → Killough scanning curves. Imbibition from an IMBNUM table or from Land's coefficient (RESSIM `SGTMAX`). Sg,max is tracked per cell (lagged to the step end). Residual CO2 = Land trapped saturation of the current scanning curve (OPM FGMTR-consistent). Without hysteresis, the cut-off is SGCR.
- **Solubility trapping + diffusion:** DIFFUSE with DIFFCWAT / DIFFCGAS (CO2STORE) and DIFFC (EOS, per component in oil and gas), using D_eff = φ·S·D on the harmonic φA/d transmissibility, fully implicit. Check: 1-D diffusion within 0.4 % of the analytical solution. Dissolution also makes the brine denser, which drives convective mixing.
- **Salt precipitation (dry-out):** PRECSALT / VAPWAT / SALTSOL with PERMFACT or Verma–Pruess (RESSIM `VERMAPRUESS`). Dissolved salt moves with the converged brine fluxes; solid NaCl lowers φ and k (lagged one step), which feeds back on injectivity. Test: 384 kg solid, φ/φ0 0.967, k/k0 0.719, mass balance 2e-16.
- **Mineral trapping:** kinetic transition-state rate law with Palandri–Kharaka k25 / Ea (presets for calcite, anorthite, albite, kaolinite), lumped CO2(aq) + silicate → carbonate, an implicit source term, the reactant inventory and the porosity change (RESSIM `MINERAL` lines).
- **Interactions** are coupled by construction:
  - Dissolution removes mobile gas.
  - Hysteresis traps the gas left behind by the imbibing brine (the chase-brine test doubles residual CO2: 931 t → 1848 t, mobile → 0).
  - Dry-out lowers injectivity.
  - Minerals consume the dissolved CO2, so more gas dissolves.
- **Mass balance:** in-place / injected = 1.000002 with all mechanisms on.
- **Six-class inventory:** mobile, structural, residual, krg=0, dissolved, mineral, plus solid salt. Exposed through CompSummary / SummaryPoint (FGMMO, FGMST, FGMTR, FGMDS, FGMMN …) and per-cell maps SGMAX, SGTRAP, SSALT, MINCO2, PORFAC, PERMFAC.
- **Compositional completion:**
  - Secondary producer limits (ORAT/GRAT/WRAT/LRAT/RESV, the most restrictive wins inside Newton).
  - SEPCOND / WSEPCOND separator trains per well (wildcard well names). Reported surface rates use each well's own train.
  - ZMF / XMF / YMF per-component maps.
  - **Group control:** GCONPROD (ORAT/GRAT/WRAT/LRAT/RESV) and GCONINJE RATE targets for GRUP wells in the group subtree. The target is split by guide rates (productivity × drawdown to the BHP limit × last yield of the target phase). Members that hit their BHP limit are fixed and the remainder is redistributed by re-solving the step. Test: a FIELD ORAT 4000 target is held at 4000.0 by two producers; a GAS RATE 5000 target gives FGIR 5000.
- **Deck (two-way):** DIFFUSE, DIFFCWAT, DIFFCGAS, DIFFC, PRECSALT, VAPWAT, SALTSOL, PERMFACT, SEPCOND, WSEPCOND, EHYSTR/IMBNUM, plus `-- RESSIM:` lines for MINERAL / SGTMAX / STRUCTTOL / VERMAPRUESS.

**App (workflow)**
- CO2 storage editor, new **Trapping & interactions** tab: hysteresis model and Land Sgt,max, diffusion coefficients, salt precipitation (SALTSOL, PERMFACT table or Verma–Pruess), a mineral reaction table with presets, the structural tolerance, and a live scanning-curve plot.
- Compositional editor:
  - **Diffusion & hysteresis** tab with DIFFC per component.
  - SEPCOND separator-train editor and per-well WSEPCOND.
- SCAL editor: new **Hysteresis** tab.
- Results → CO2 storage: five-class stacked inventory versus time, trapping fractions, and the new maps (also in the 3D viewer). ZMF/XMF/YMF maps come with a component selector.
- New example **`co2trap`**: the CO2 aquifer plus a two-layer caprock, Carlson hysteresis (Sgt,max 0.25), 2 m NaCl, diffusion, salt precipitation and anorthite carbonation. 50 y in ~35–70 s (1233 steps, 490 cuts, caused by hysteresis switching at the plume edge; step control is a next item). At year 50: mobile 188 kt (171 kt of it structural, under the caprock), residual 501 kt, dissolved 276 kt, mineral 127 kt (A = 20 m²/m³). Total in place equals the 1.09 Mt injected.
- **Tests:** 250 checks pass. New checks: Land / Carlson / Killough curves, the chase-brine residual trapping, diffusion versus analytical, salt dry-out, mineral carbonation, secondary limits, SEPCOND round trip, the all-mechanisms deck, and compositional group control.

**Known limitations (to do)**
- CO2SOL (CO2 in the aqueous phase in EOS mode) and the SPE3 gas-condensate comparison are deferred.
- No group guide-rate formulas (GUIDERAT) or REIN / VREP / RESV group injection in the compositional core; these are run by OPM.
- Pc hysteresis is not modelled. Salt in the coarse examples does not dry out (cells are too large); a near-well refinement (LGR) would be needed.
- Mineral kinetics use one lumped reaction per mineral; there is no full aqueous speciation (PHREEQC-class geochemistry).

## Round C6: compositional speed, ECL/OPM result detail, CO2 as an EOS component (CO2SOL), CO2 EOR workflow

**1. Compositional engine speed (advanced numerics, same physics)**
- Flash: in-place LU with partial pivoting (luFactor / luSolveF) instead of per-call allocations. The previous step's K-values warm-start a direct Newton (the stability test is skipped while the cell stays two-phase).
- **Absent components (z_i = 0):** found and fixed this round. These are not unknowns of the flash Newton (identity row, zero column). Before the fix, Newton never converged when a component was missing from a cell (CO2 before it arrives, for example) and the flash fell back to successive substitution. Their K now comes from φ_L/φ_V, a trace partition that also gives the warm start when the component arrives. flashSensitivities uses the trace partition dv_i = K_iV/(L+K_iV)·dz_i, and the finite-difference check passes at 5e-7.
- flashSensitivities reuses the phase derivatives of the flash (no re-evaluation; T derivatives only when needed).
- evalCell is OpenMP-parallel. The CPR keeps its AMG hierarchy across Newton iterations (BlockCPR setup reuse). The bsr kernels are templated on the block size. Inexact Newton uses a linear tolerance of 1e-5. Steps are cut early when divergence is detected. Time-step control follows Aziz–Settari (target ΔS / Δz / Δp). Relaxed CNV applies only on the last iterations. Per-phase timers are available (`RESSIM_COMP_DEBUG`).
- **Timings (2 cores):**

  | Run | Before | After |
  |---|---|---|
  | SPE5 (20 y, 665 steps) | 10.5 s | 6.2 s |
  | `-co2sol` test group | 110 s | 17 s |
  | Full test suite | 2.5–5 min | 128 s |
  | co2eor (10 y) | – | 6.3 s |

- The SPE5 results move by less than 0.05 % (FOPT 2.4284e7 vs 2.4277e7) because the flashes are now fully converged.

**2. ECLIPSE / OPM results in full detail**
- The restart reader exposes every cell array (generic NAME1..k per-component arrays included) through a 60 MB snapshot budget.
- Results/3D viewer: **AMF / XMF / YMF / ZMF** per component, plus all OPM CO2STORE arrays and any other restart array, under "Engine output arrays".
- RPTRST is written per case type. Compositional runs request XMF/YMF/ZMF/DEN/VISC/PSAT; CO2SOL adds AMF/RSW.

**3. Separator semantics (the GOR / FOIP offset in the user's ResSim-vs-ECLIPSE plots)**
- E300 FIELDSEP: the **last stage is the stock tank**. The old comp deck had a single 80 °F / 315 psia stage. ECLIPSE took that stage as the stock tank (FOIP 44.0e6, GOR 0.44), while ResSim added a standard-conditions stage (FOIP 42.4e6, GOR 0.54).
- ResSim now uses the E300 rule. Examples and exported decks carry an explicit 2-stage train (separator + 60 °F / 14.7 psia stock tank), so both simulators see the same process. The compositional editor shows the rule and has an "Add stock-tank stage" button.

**4. CO2 = an ordinary EOS component (no separate CO2 module for EOR) + CO2SOL**
- CO2 is component 7 of the PR78 model (`spe5Co2CompFluid`: Tc 304.13 K, Pc 73.77 bar, ω 0.2239; BIC 0.105–0.125 with the hydrocarbons). Injection uses WELLSTRE / WINJGAS streams.
- **CO2SOL (E300):** the CO2 component also dissolves in the brine. Aqueous CO2 comes from fugacity equality, f_CO2(HC phases, EOS) = x_aq·H(p, T, salinity). H is calibrated so pure CO2 reproduces Spycher–Pruess / Duan–Sun. Aqueous moles are solved per cell by safeguarded Newton.
- Jacobian derivatives follow the implicit-function theorem. The finite-difference check passes at 1.4e-8.
- The aqueous phase volume uses the Garcia apparent molar volume; viscosity uses Islam–Carlson.
- The initial brine is CO2-saturated (the E300 default without RSWVD).
- Deck: CO2SOL + SALINITY import/export; SOLUBILI / SOLUBILT / RSWVD / AMFVD / AMF / H2SSOL are known keywords. Without SOLUBILI tables, ECLIPSE uses its Chang–Coats–Nolen defaults.
- Validation: 150 bar, 50 °C, 1 m NaCl gives x_aq 0.01758 vs Spycher–Pruess 0.01766. The co2eor balance error is 5e-12.

**5. CO2 EOR workflow (Fluid → Compositional → CO2 EOR tools)**
- **Swelling test / P-X diagram:** bubble/dew pressure versus mol % injection gas, swelling factor, gas added (scf/bbl), saturated density. The **first-contact miscibility pressure** is the refined cricondenbar of the pseudo-binary. A stability-test bisection handles near-critical mixtures where the Psat Newton finds the trivial solution.
- **Multi-contact MMP, mixing-cell method (Ahmadi & Johns 2011):** combined vaporising/condensing drive. The minimum key tie-line length is extrapolated to zero. Tie-line length and per-cell profile plots are shown, with Yellig–Metcalfe and Cronquist correlation cross-checks.
  - SPE5 oil at 160 °F: CO2 MMP is 2790 psia (Cronquist 2653). The SPE5 lean solvent MMP is 3386 psia. FCMP with CO2 is 3174 psia. Pb is 2302.89 psia.
- **WAG schedule generator:** injector, gas stream, start, gas/water half-cycle lengths, cycles, rates, BHP limit, tapered gas slugs, water chase, shut-in. It splits report steps at the half-cycle boundaries and replaces the injector's controls. The exported deck carries alternating WCONINJE plus WINJGAS STREAM.
- Results: a "CO2 (EOR / CO2SOL)" tab splits CO2 between the hydrocarbon phases and the brine. 3D maps: AMF (aqueous CO2), XCO2AQ, total CO2 mole fraction, ZMF/XMF/YMF.
- Compositional editor: **SPE5 + CO2** preset, CO2SOL toggle with salinity.
- **New example `co2eor`:** SPE5 + CO2 with CO2SOL, 7×7×3. One year of depletion, then 1-year CO2/water WAG. 10 y run: 205 steps, 6.3 s. CO2 at the end: 5.2e5 t in the HC phases, 5.2e4 t dissolved in the brine.
- **ECLIPSE testability:** `release/examples/CO2EOR.DATA` is a complete E300 deck (COMPS, CO2SOL, SALINITY, EOS PR + PRCORR, CNAMES …, FIELDSEP 2 stages, WELLSTRE, WINJGAS, RPTRST with AMF/RSW). OPM Flow has no CO2SOL, so run this deck with E300.

**Tests:** 259 checks pass, all new this round:
- CO2SOL equilibrium and Jacobian
- co2eor balance
- warm-start vs cold flash (9e-13)
- absent-component sensitivities
- swelling / FCMP
- mixing-cell MMP
- WAG generator
- CO2EOR deck round trip

**Known limitations (to do)**
- No SOLUBILI/SOLUBILT tables in the own engine; it uses the Spycher–Pruess-calibrated Henry function. Water vaporisation into the HC phases is ignored in EOS mode.
- MMP is from the mixing-cell method only; a 1-D slim-tube simulation option is still to do.
- SPE3 gas-condensate comparison and GUIDERAT for compositional groups.

## Round C7: compositional completion + STARS-like thermal simulator (E300 THERMAL)

**Status:** 271 tests pass, both native and on Windows (MinGW build under wine). The release includes ResSimStudio.exe, the 3 new thermal example decks, and `ENGINES_AND_MATH.md`, the new reference for every engine's mathematics and solution methods.

**1. Compositional items completed**

| Item | What was built | Check |
|---|---|---|
| SOLUBILI | E300 table (P, Rs_w, Bw, cw, μw) for CO2SOL: the saturated aqueous CO2 fraction comes from the table instead of the Spycher–Pruess calibration. Bw, cw and μw are read but not used. | CO2 dissolved after the test run: 3378 t vs 1713 t with the default model (higher table solubility) |
| Slim-tube MMP | 1-D compositional displacements on the FIM engine (100 cells, 12 m, 1.2 HCPV) at several pressures. MMP from the break-over of recovery versus pressure, plus the pressure where recovery reaches 90 %. UI shows the recovery and profile plots. | SPE5 oil with CO2: 2743 psia (mixing cell 2790, Cronquist 2653) |
| GUIDERAT for compositional groups | GUIDERAT formula (OIL/LIQ/GAS/WAT/RES phases, A–F coefficients) for allocating group targets in the compositional core. Import, export and a UI button. | test: guide-rate split 3239 / 761 |
| Water vaporisation | VAPWAT / `-- RESSIM: WATERVAP`: water vapour in the gas phase (free-water Raoult law, y_w = Psat(T)/p from IAPWS-IF97), dry-out, exact derivatives, water flows in the gas phase, `vaporisedWater()` report, UI checkbox | 400 K comp example: 17,700 kmol of water in the gas phase, water balance 6e-14, deck round trip |
| SPE3 comparison | **Deferred.** The SPE3 EOS parameters are needed from the user's E300 SPE3 deck. | – |

**2. Thermal simulator (STARS-like, unified with the compositional component model)**
- **Core (`thermsim`):** FIM in (p, N_1..N_nc, N_w, U), with mass, energy and volume equations. Uses the IAPWS-IF97 steam
  tables and IAPWS-2008 viscosity, KVCR K-values with a free-water Rachford–Rice flash, a (p, N, U) energy flash with a
  univariant / near-univariant energy blend, and Vinsome–Westerveld heat loss. Rock heat capacity is HEATCR/HEATCRT and
  conductivity THCONR/THCONSF. See `ENGINES_AND_MATH.md` §6.
- **E300 THERMAL deck, two-way:**
  - Section keywords: THERMAL, DEADOIL, COMPS, CNAMES, MW, TCRIT, PCRIT, CVTYPE, KVCR, KVWI, DREF, CREF, THERMEX1,
    SPECHA/B/G/H, HEATVAP(S/E), ZFACTOR, OILVISCT, GASVISCT/GASVISCF, WATVISCT, STCOND.
  - Rock and initial state: HEATCR, HEATCRT, THCONR, THCONSF, TEMPI/TEMPVD, ZI, XMF/YMF.
  - Wells: WINJTEMP.
  - Heat loss: `-- RESSIM: HEATLOSS`.
  - Export to ECLIPSE always writes THSTT97.
  - Round trip: reimport gives the same model, with summaries identical across both passes.
- **Examples (`examples.cpp`, METRIC heavy oil HEAVY + LIGHT):**

  | Example | Setup | Run time (2 cores) |
  |---|---|---|
  | `steamflood` | quarter five-spot, 5 years | ≈ 80 s |
  | `css` | 3 cycles, 2 m near-well refinement | ≈ 100 s |
  | `sagd` | 2-D well pair, LRAT producer | ≈ 110 s |

- **Robustness:**
  - Path-independent energy flash.
  - Appleyard chopping in (p, N, U).
  - Oscillation relaxation.
  - Relaxed CNV with MB-checked convergence.
  - Well-connection hysteresis.
  - Small restart steps after a well control change.
- **App:**
  - **Fluid model:** a "Thermal (STARS-like)" option opens the thermal editor. Tabs:
    - Components: KVCR table, K(T) plot, ZI, heavy-oil preset.
    - Viscosity: OILVISCT / GASVISCT tables with a μ(T) plot including the IAPWS water curve.
    - Rock & heat loss.
    - Initial temperature: TEMPI or TEMPVD.
    - Steam injection: WINJTEMP table with saturated state and enthalpy, plus a steam-table calculator.
  - **Run, init and validation:** runs on the thermal core (`SimRunner::startThermal`), initial-state preview, "Thermal
    model valid" status.
  - **ECLIPSE engine:** runs e300 automatically for THERMAL and COMPS decks.
  - **Results:**
    - New Thermal tab: average and maximum T, steam PV, energy balance (net injected enthalpy, heat lost), cSOR and oil
      rate.
    - Summary vectors FAVGT, FMAXT, FSTPV, FHLT, FHIN, FCSOR.
    - TEMP, SSTEAM and ZMFk maps in the 3D viewer.
- **Tests (`-thermal`):**
  - IF97 / IAPWS reference values.
  - Univariant flash.
  - Conduction box.
  - Vinsome–Westerveld within 5 % of the analytic solution.
  - 1-D steam vs cold flood energy balance below 1e-5.
  - Deck round trip of the three examples.
  - Short runs of the three examples with an energy balance check.
- **Fix:** `ReportStep{len, {}}` aggregate initialisation had been silently writing into `guiderat` since GUIDERAT was
  added (app demo case and "Insert step after").

**Still open (thermal):**
- Capillary pressure.
- Injector wellbore heat / pressure solve, and discretised wellbore (STARS FlexWell).
- Steam-trap (subcool) producer control.
- Combustion / chemical reactions (STARS REACTION).
- Rock dilation / geomechanics for CSS.
- Analytic AD derivatives of the energy flash.
- EOS-based thermal phase behaviour: plugging the cubic EOS into the thermal core; the component model is already
  shared.
- ECLIPSE ROCKCON / ROCKPROP heat-loss keywords: the item formats are unknown, so they are kept as `-- RESSIM` lines.

## Round C8: full thermal — reactions, in-situ combustion, UCG, geomechanics (THM / THMC)

**Status:** 276 tests (5 new in the `-thm` group) pass natively and on Windows (MinGW build under wine). Release with ResSimStudio.exe and the new example decks COMBUSTION / UCG / SAGD_THM.DATA. Mathematics: `ENGINES_AND_MATH.md` §6.1–6.3.

| Area | What was built | Check |
|---|---|---|
| Component types | K-value fluid / non-condensable gas / solid (coke, coal, char); solids occupy pore space and carry energy | deck `-- RESSIM: CTYPE`, CVTYPE DEAD for solids |
| Reactions | STARS-form Arrhenius kinetics, any stoichiometry (nc + 1 incl. water), reaction orders, O2 partial-pressure factor, minimum temperature, reversible (K_eq(T)) reactions, exothermic / endothermic heat, fully implicit | energy balance with reaction heat 4e-8 |
| In-situ combustion | combustion-tube example: cracking → coke, heavy / light / coke burning; air injection via WELLSTRE / WINJGAS gas injectors | front propagates 0.1 m/h, peak 630 °C, steam plateau and oil bank |
| UCG / coal | pyrolysis, char combustion, partial oxidation, steam gasification, Boudouard, water-gas shift (reversible), CO / H2 / CH4 burning; cavity growth through solid consumption; Kozeny-Carman / exponential permeability | syngas CO + H2 + CH4, φ 0.05 → 0.6, k × 3700 |
| Steam tables | IAPWS-IF97 region 5 (to 2273 K) | v, h at 1500 K / 0.5 MPa to 1e-8 |
| Geomechanics | hex8 FEM thermo-poroelasticity, padding over / under / side-burden, rollers / free top, PCG-SGS; fixed-stress iterative or one-way coupling; initial stress, Mohr-Coulomb shear and tensile failure indices, heave | uniaxial column exact (4e-10); SAGD THM heave 3.6 mm surface / 19 mm reservoir top at 91 d |
| Flow physics review | harmonic conduction, P_c in thermal, PV law (c_r, c_T, solids, mechanics), permeability alteration on transmissibilities and WI | all C7 thermal tests unchanged (energy balance < 1e-5) |
| Solver review | relaxation recovery scaled with ω, Aziz-Settari step control with post-cut growth limit, per-cell damping, rate/BHP switch limit | SAGD 91 d: cuts 66 → 30, 25 → 19 s |
| Deck | `-- RESSIM: CTYPE / REACTION / SOLIDCONC / PERMLAW / GEOMECH`; WELLSTRE / WINJGAS in thermal decks; OPM MECH / YMODULE / PRATIO / BIOTCOEF / THERMEXR import | export → import round trip of all THMC data (test) |
| UI | Thermal editor: component Type column, **Reactions** tab (stoichiometry grid, mass-balance check, Arrhenius plot, solids, permeability law), **gas streams** + WINJGAS assignment, **Geomechanics** tab (moduli, padding, stress, failure, coupling); Results: Reactions and Geomechanics panels, maps UPLIFT / STRAINV / SIGEFF* / FAILURE / POROF / PERMMULT / SOLID_* | MinGW build |

**Not yet (C9 candidates):** plasticity (return mapping) and fracture / dilation permeability; AMG preconditioner for elasticity; calibrated kinetics libraries; ECLIPSE/OPM export of reactions (no equivalent keywords: kept as `-- RESSIM:` lines with a warning); BCMECH boundary loads; mechanics on true corner-point geometry; SAGD steam-trap (subcool) producer control; thermal runtime target (SAGD 2 y ≈ 90 s).

## Round C9: adjoint history matching, multi-component product analysis, steam trap, mechanical boundary loads

**Status:** 291 tests pass natively and on Windows (MinGW build under wine). New test groups: `-adj` and `-mc`, plus new checks in `-thermal` and `-thm`. Mathematics: `ENGINES_AND_MATH.md` §12.1, §12.2, §6 (steam trap) and §6.2 (boundary conditions).

| Area | What was built | Check |
|---|---|---|
| Discrete adjoint | Exact misfit gradient of the black-oil model: recorded forward run, replay Jacobians, reverse recursion with the explicit well heads as an augmented lagged state (multipliers μ), and rate observations differentiated through the connection heads. ∂R/∂θ comes from residual differences at the recorded state, so every study parameter type works. Parameter contractions run in parallel. | adjoint = forward mode to 9e-9; vs FD of full runs 6e-7 on a horizon without switching, 1.4e-3 across control switches |
| Forward sensitivities | Tangent mode with head sensitivities. All dy/dθ go to the LM Jacobian (no perturbed runs); parameter columns run in parallel. | LM analytic recovers the truth in 17 runs |
| Adjoint L-BFGS | Projected L-BFGS + Armijo in the normalised box. One forward + one backward pass per iteration, independent of the parameter count. | truth recovered, misfit 14.7 → 0.022 in 10 runs (was 21 before the head coupling) |
| Studies UI | *Optimiser / gradient* selector (LM-FD, LM-analytic, adjoint L-BFGS) with its options; **Check gradient** table (adjoint / forward / FD, relative errors, timings); HM views label the method | MinGW build |
| Multi-component analysis | For every multi-component engine (EOS compositional, CO2STORE, thermal / combustion / UCG): per-well and field component molar rates (report-step averages from cumulative moles), cumulatives, species recognition, HHV/LHV, stream analysis (dry / gas composition, H2, syngas CO + H2, H2/CO, LHV/HHV power and MJ/sm³, CO2 / O2 slip) | Σ rate·Δt = produced moles (exact) for comp and UCG; UCG product gas H2/CO 0.31, LHV 10.4 MJ/sm³ |
| Results UI | **Components** tab: stream selector (field / well), rate basis (kmol/d, t/d, sm³/d), component table at the time cursor, per-well breakdown, and plots of rates, dry composition, cumulatives, hydrogen / syngas and heating value (vector export) | MinGW build |
| Steam trap | STARS OPERATE STEAMTRAP: implicit dynamic BHP floor P_sat(T_hot + ΔT_sub) combined with the rate / BHP controls. `-- RESSIM: STEAMTRAP well dT` (SCHEDULE); thermal editor table; SAGD example uses 10 °C | subcool held to 1e-5 K, producer-cell Jacobian 8e-4 |
| Mechanical BCs | BCMECH-like per-face roller / free / fixed with normal stress increments (consistent nodal loads), rigid-body guard; `-- RESSIM: MECHBC`; Geomechanics editor table | oedometric surface load: settlement, Δσ'_v, Δσ'_h exact to 1e-8; deck round trip |

**Not yet (C10 candidates):**
- Adjoint for the compositional / thermal cores.
- Per-cell multiplier-field gradients (face-flux adjoint) for pilot-point / zonal matching.
- Gradient-enhanced ES-MDA.
- Group / THP targets differentiated (currently frozen).
- OPM BCMECH import (needs BCCON face regions).
- Plasticity (return mapping) and dilation permeability.
- AMG for elasticity.
- Mechanics on true corner-point geometry.
- E300 / OPM per-component summary vectors in the Components tab.
- Thermal runtime target (SAGD 2 y ≈ 90 s).

## Round C10: unified multi-component framework — solubility, PR gas and trapping in the thermal core, per-component 3D maps, UCG + CO2

**Status:** 299 tests pass natively (291 + 8 new in `-sol`) and on Windows (MinGW build under wine). Mathematics and the
research notes: `ENGINES_AND_MATH.md` §6.4 and §6.5.

**Goal (user request).** Thermal runs show dynamic per-component 3D maps like compositional runs. CO2 tracking and
trapping work in thermal exactly as in compositional, with no separate module (for example UCG syngas + CO2 storage in
one run). The design is universal, so users build cases from reactions upward. Research came first.

| Area | What was built | Check |
|---|---|---|
| Gas / solute library | `gasdata.h/.cpp`: 14 species with PR critical constants, IAPWS G7-04 Henry coefficients, v∞, Sechenov k_s, c_p and formulas; Wagner–Pruss p_sat; formula parser (parentheses, fractional counts); element-balance and carbonate-CO2 helpers | G7-04 check values to 5e-5 |
| Three-phase convex flash | `flashGeneral`: aqueous / oil / gas from fugacity constants (Henry + Krichevsky–Kasarnovsky + salting-out; KVCR; φ p). Minimises Michelsen's convex Q(β) with a projected Newton active-set method; outer φ loop for the PR gas. Used automatically when a component is soluble | ≡ legacy flash to 1e-13 (6 states incl. supercritical water) |
| PR gas | `GASEOS PR`: PR sub-mixture (library / CRIT constants) for φ, volume, residual enthalpy and LBC viscosity | CO2 solubility vs Spycher–Pruess within 5.4 % (8–24 MPa, 50 °C) |
| Aqueous phase with solutes | Dilute partial molar volumes / enthalpies (with dissolution enthalpy from d ln k_H/dT), Islam–Carlson CO2 viscosity; dissolved gas reported as surface gas; aqueous reaction concentrations (`PHASE 0`) | closed cell: exact mass / energy, Jacobian 2.6e-7 |
| Gas hysteresis | `HYSTGAS`: Land trapping + Carlson k_rg from the historical S_g,max (same as the compositional core); SGMAX / SGTRAP maps | 1D CO2 + water chase: residual 129 t (Land) vs 62 t (S_gc) |
| Storage inventory | `TRACK` (CO2 by default) with dissolved-water / dissolved-oil / residual / structural / mobile / mineral classes; mineral = Σ ΔN · trackEq · M of solids with `FORMULA ... trackEq` | inventory closes to round-off |
| Per-component 3D maps | Thermal runs write ZMF / XMF / YMF / AMF for every component (+ water) at each report step, plus SGTRAP / SGMAX; the viewer's component selector uses the run's component list (it was empty for thermal) | MinGW build |
| Results UI | The "CO2 storage" tab becomes "<tracked> storage" for thermal runs (in-oil class added) | MinGW build |
| Thermal editor | New **Solubility & storage** tab: PR gas, salinity, hysteresis, tracked component, structural tolerance; per-component Henry model / coefficients / v∞ / k_s / Tc / Pc / ω / formula / trackEq with a **Library** button; solubility-vs-T preview. Reactions now show an **element-balance** check (closed / per-element imbalance / missing formulas) | MinGW build |
| Deck | `-- RESSIM:` AQSOL / GASEOS / CRIT / SALINITY / HYSTGAS / TRACK / FORMULA; long RESSIM lines now continue on `-- RESSIM+` lines (they were truncated at 128 columns, which broke long REACTION lines on re-import, including the C8 UCG example) | round trip of `ucg_co2` |
| Thermal summary | Gas injection rate / total (FGIR / FGIT: air, O2, CO2 and N2 streams at standard conditions) now reported for thermal runs | `-sol` UCG run |
| Example `ucg_co2` | 30 d UCG gasification → 30 d water quench → 60 d CO2 injection into the cavity → shut-in; dissolved / residual / mobile / mineral (ash CaO carbonation by aqueous CO2) classes all develop; Boudouard conversion of CO2 on hot char is visible | 20-report test run, energy balance 4e-7 |

**Research outcome and decisions (details in §6.5):**
- EGS needs EDFM fracture cells and hydro-shear (Barton–Bandis / Coulomb). These are the next gridding / physics
  additions. Radial and LGR refinement around wells is also more valuable than unstructured PEBI.
- UCG CO2 storage literature points to char/coal adsorption: a competitive extended-Langmuir sorbed state is planned.
- Nanoparticle transport / retention / wettability is planned on the polymer-adsorption pattern.
- H2 colours (brown / black / blue) are covered by gasification chemistry + product analysis + CO2 storage in one run.

**Open (C10):**
- continuous supercritical water;
- extended-Langmuir sorption;
- EDFM / hydro-shear;
- nanoparticles;
- adjoint for the comp / thermal cores (carried over).

## Round C10.1: table and text-cursor fixes (UI), verification status, adjoint vs JutulDarcy

| Problem (user report) | Cause | Fix |
|---|---|---|
| Table columns could not be resized | Most data tables were created without `ImGuiTableFlags_Resizable` | All bordered / striped tables go through `ui::beginTable` / `ui::tableColumn` / `ui::endTable` (`app/src/ui_table.cpp`): always resizable (double-click a border = auto-fit). Property grids (label / value) are resizable too |
| Values cut off; column width followed the header text | Input cells sized with `-FLT_MIN` report no content width, so fixed-fit columns took the header width, and stretch columns got weights from the header text | Every data column is at least 5.5 em wide. When the columns fit, they stretch with weights max(header, minimum). Otherwise fixed widths + horizontal scrolling (with hysteresis, so the mode does not flicker) |
| Reaction tables got thinner with every reaction and disappeared | Scrolling tables without a height take "the rest of the window", which goes to ~0 further down a scrolling page | Scrolling tables without an explicit height use their content height from the previous frame (+ scrollbar), capped at max(visible space, 22 lines). Panel tables (deck report, problems, observations) explicitly fill their panel |
| No visible text cursor anywhere | `ImGuiStyle::ScaleAllSizes` truncates `InputTextCursorSize` (1 px × 0.9 → 0 px), and re-applying the UI scale compounded the scaling of every metric not reset by the theme | The style is rebuilt from unscaled defaults before scaling. Caret thickness = max(1, round(1.5 × scale)) px with a bright caret colour (input fields and the deck editor) |

Checked visually in the MinGW build under wine (virtual display): UCG + CO2 example, all 10 reactions expanded and their
stoichiometry tables fully visible, column drag-resize, and the caret in a numeric field.

**Verification / benchmark plan (next rounds):**
- SPE11 (CO2 storage, 2023 comparative project) against published results;
- SPE4 thermal comparative solution (steam injection);
- Egg model and Norne history matching with the adjoint against JutulDarcy;
- a STARS / E300-THERMAL cross-check when the user can supply runs.

**Adjoint parity with JutulDarcy (gaps):**
- per-cell / per-face parameters by analytic dR/dθ (perm, poro, transmissibility, WI), instead of residual differences
  whose cost grows with the parameter count;
- block-CPR solves of the transposed system (now ILU(0)-GMRES, which will not scale to Norne);
- implicit well heads instead of a frozen-coefficient treatment of lagged data;
- adjoint for the compositional and thermal cores;
- checkpointing;
- measured forward speed parity.

## Round C11: supercritical water, injector stall fix, analytic verifications, adjoint plan

**C10 open items, status:**

| Item | Status | What was done |
|---|---|---|
| Supercritical water (closed hot cavities above p_c) | ✅ | IF97 region 3 in `steam.cpp`, with the root bracketed by the Wagner–Pruss auxiliary saturated densities (the earlier ×1.06 density walk could skip the near-T_c loop); C1 extension of K_w / H_wA through T_c; enthalpy and internal energy of water in a gas mixture capped at the 0.75 ρ_c isochore (`steam::vapourMixCap`), so U(T) stays monotone (scan: negative slopes 1.4e4 → ≤ 150 J/K) |
| Injector stall at the BHP limit (ucg_co2 at 800 sm3/d, 50 bar, failed at day 62.6) | ✅ | Root cause: backflow of injection-only components, clamped at N = 0 (the BHP-kink hypothesis was wrong). Fix: the same C1 no-backflow ramp in both control modes, scaled per well. Now 412 steps, 31 cuts, runs to the end |
| ucg_co2 at 150 sm3/d, 120 bar | 🟡 | Still fails at day 244.3: one cell at 22.86 MPa sits on the 646–648 K fold of the CO2–water mixture energy (ideal mixing of real gases near water's critical point). Next: a real-mixture EOS for water + gases, or an equal-area construction across the fold in the energy flash |
| Langmuir sorption, EDFM, nanoparticles | ❌ | next rounds (research in §6.5 of ENGINES_AND_MATH) |
| Adjoint parity with JutulDarcy | ❌ planned | `adjoint.cpp` was read: ∂R/∂θ comes from perturbed-Case replays, so the cost grows with the parameter count. Plan: analytic per-connection dR/dT (flux is linear in T; harmonic average ∂T/∂m_a = T²/(m_a² t_a)), WI ∝ k, plus a transposed CPR solve. This makes HM gradients fast for any number of parameters, but **does not** speed up forward runs (forward speed is a separate item: linear solver / GPU) |

**New verification group `-scw` (7 checks, all PASS):**
- region-3 check values (6.9e-8) and isobar continuity;
- flash continuity at T_c (7.8e-9);
- supercritical box heated from below (p_inj −152 Pa, 0 cuts);
- injector Jacobians (rate 1.0e-3, BHP 9.4e-4);
- 1D conduction against erf (0.015 K of 20 K);
- Theis against E1 (worst 1.16 %; implicit Euler needs report steps ≤ 0.025 d, since 0.25 d gives 6 %).

Suite: 307 PASS natively and under wine (299 + 8). Lauwerier and Terzaghi / Mandel are next on the analytic list. The big runs (SPE11, SPE4,
Norne HM against JutulDarcy) are left for later rounds because of their time cost on this machine.

## Round C12: remaining verifications, sorption, EDFM, field adjoint, speed profile

Scope requested: every remaining verification except the Norne history match; start every partial or not-started item;
the HM adjoint upgrade if feasible; forward-speed work.

| Item | Status | What was done |
|---|---|---|
| Lauwerier, Terzaghi (analytic) | ✅ | `-ver`: Lauwerier \|ΔT_D\| < 0.05; Terzaghi error 0.001 |
| Mandel consolidation | ❌ planned | needs a 2-D plane-strain load case; left for the next round |
| SPE4 (steam CSP) | ✅ problem 1A | `-spe4`: OOIP / OWIP within the participants' range; heat loss after cycle 1 616 GJ (491–703); injected heat 3,430 GJ vs 3,594 (BHP-limited start-up, ECLIPSE / CMG injectivity); 18.3 MSTB oil after 3 cycles. Problems 2A / 3 (nine-spot, grid orientation) not run |
| SPE11 (CO2 CSP) | 🟡 SPE11A-like | `-spe11`: SPE11A fluids, facies 1/3/4/5 data, wells and schedule on a layered slice. Solubility within 2 % of Henry; exact mass balance; seal holds free gas; full dissolution at 5 d. The real facies geometry and SPE11B/C are still open |
| Egg model | ✅ runs | EclGym deck, 3600 d: FOPT 558,594 m³, 0 cuts, 28.1 s on 2 threads (29.1 s before the CPR change) (JutulDarcy reports 8.6 s on its own hardware); profile in ENGINES_AND_MATH §7 |
| Langmuir / ECBM sorption | ✅ | `SORB` engine + deck + UI (Thermal → Sorption); example `ecbm`; `-sorb` |
| EDFM + Barton–Bandis + hydro-shear | ✅ | `FRACTURE` / `FRACSTRESS`, well–fracture WI, 25 % aperture cap per step, FRACAP / FRACP / FRACT maps, UI (Thermal → Fractures), example `egs_edfm`; `-edfm` (cubic law 0.07 %, influx 2.1 %, deck round trip) |
| Near-critical fold (ucg_co2 150 sm3/d) | ✅ | two-state fold construction on a fixed sample grid (§6.6); runs through day 270 (failed at 244.3 in C11) |
| HM adjoint upgrade | ✅ | analytic per-cell / per-face / per-connection dR/dθ (`-adjf`: 4.3e-6 vs FD; gradient 0.25 s vs forward 0.70 s); per-cell field HM in Studies |
| Forward speed | 🟡 | Newton-loop profile (`rundeck` PROFILE, `RunStats` timers); Jacobian dump + benchmark; cheaper CPR second-stage residual. The Krylov phase is bound by memory bandwidth. Next: mixed-precision CPR, multicolour smoothers |
| Compositional / thermal adjoint, plasticity, AMG for elasticity | ❌ not started | listed in ENGINES_AND_MATH §15 |

**Fixes found by the new cases (self-review):**
- **Thermal Newton iteration 0.** Converged-at-start was accepted with the loose absolute tolerances, so slowly diffusing
  cells skipped steps: a bias that depended on Δt (Lauwerier / Theis). Iteration 0 now requires a 1e4× tighter state.
- **Thermal temperature chop.** |ΔT| ≤ 80 K per update, needed for EDFM fracture cells, which have no rock heat capacity.
- **gasViscT / oilViscT tables** in `exEgsEdfm` and in the E300 exporter had nc + 2 columns instead of nc + 1. Re-import
  then misread the table. Both are fixed, and the EDFM deck round-trip test now covers it.
- **Injector drawdown-ramp FD check.** It was a test artifact: the state sat within the FD step of the BHP limit. The check
  now runs at day 10 on a copy of the state.
- **SPE11A set-up lesson:** a producer alone is not a constant-pressure boundary. Pair it with a BHP water injector.

Suite: **331 PASS natively and under wine** (307 in C11 + 24 new: `-ver`, `-sorb`, `-edfm`, `-spe4`, `-spe11`, `-adjf` and the EDFM deck round trip).

## Round C13: dedicated adjoint HM method, EGS designer + fracture visualisation, industry comparison

Scope requested: make the gradient / adjoint HM a dedicated method with full result views (iterations, final parameters,
per-response check of the best case); make EDFM / EGS usable from the UI, draw the fractures in the grid, and provide
detailed EGS analysis; keep running bug analysis; do the feasible open items first; compare with CMG / ECLIPSE / OPM Flow.

| Item | Status | What was done |
|---|---|---|
| Adjoint HM as its own method | ✅ | Studies → History matching → **"Adjoint gradient (L-BFGS / per-cell field)"** (method 6). Unknowns: the scalar study parameters, or per-cell PERM / PORO multiplier fields. The old hidden per-cell block under LM is gone. LM keeps its own FD / analytic Jacobian switch |
| Adjoint HM centre view | ✅ | new **Adjoint HM** tab: headline (initial / best misfit, reduction, runs, forward / adjoint seconds); iteration slider; **Iterations** (misfit + \|g\| on log axes, click to select, parameter or multiplier-range evolution, iteration table); **Final parameters** (prior / iterate / best, normalised position, gradient; apply best or any iterate); **Field maps** (ln k / ln φ multipliers and the initial / final sensitivity maps per layer, symmetric colour scale, histogram, *Show in the 3D viewer* as HMMULTK / HMMULTPHI / HMSENSK); **Responses** (normalised RMS per observation group for prior / iterate / best, per-response fit plots, grid mode); **Sensitivity** (gradient bars at iteration 0 vs best, layer-wise \|∂J/∂ln k\|). The field result also appears in the Match quality table |
| Engine records | ✅ | L-BFGS keeps every iterate (values, gradient, run, cumulative runs and time); the field HM stores ln k / ln φ snapshots per iterate and the initial / final sensitivity maps; `-adj` checks both |
| EDFM editor + EGS designer | ✅ | Thermal → **Fractures & EGS**: interactive plan view (draw / drag fracture traces over the grid and the wells), vertical section of the selected fracture, fracture table with Barton–Bandis / Coulomb parameters, in-situ stress, EDFM pre-check (fracture cells, well–fracture connections, effective normal stress and slip tendency per fracture). The **multistage designer** (stages, spacing, well separation, half-length, height, aperture, stresses, rates, temperatures, near-fracture refinement) builds a complete case; new example `egs_multistage` |
| Fractures in the grid | ✅ | fracture planes drawn in the 3D viewport of the Model and Results workspaces and in the result viewer, coloured by the fracture value of the current property (aperture, temperature, pressure) |
| EGS analysis | ✅ | Results → **Geothermal**: *Performance* (production / injection temperature, normalised thermal drawdown with a lifetime threshold, thermal / electric power, cumulative heat and heat recovery factor, injectivity / productivity / impedance, per-fracture aperture band, temperature, conductivity and cooled rock volume); *Fractures* (fracture-face maps and along-fracture profiles per time step); *Cooling map* (plan-view cold plume); *Connections & summary* (well–fracture links, lifetime and energy figures) |
| Gringarten verification | ✅ | `-grin`: EDFM single fracture vs Gringarten (1975) T_D; near-fracture refinement gives max \|ΔT_D\| 0.031, a uniform 25 m grid 0.242. The designer therefore refines by default (2 m first cell, ratio 1.8) |
| MINC, mixed-precision BILU0, Mandel, SPE4 2A | ❌ next round | listed in the gap audit and ENGINES_AND_MATH §15 |

**Fixes found in this round (self-review):**
- The coarse EGS grid over-predicted early cooling by ~35 °C (Gringarten test). Fixed by the near-fracture refinement in the designer.
- `propIndexByName` was defined inside an anonymous namespace while declared globally (would not link); moved.
- A vexing-parse `std::vector<double> kk(size_t(nz))` in the sensitivity view; fixed.

### Gaps and advantages against CMG, ECLIPSE and OPM Flow (honest status, C13)

Reference products: CMG IMEX / GEM / STARS + Builder / Results / CMOST; SLB ECLIPSE 100 / 300 (+ Petrel, INTERSECT);
OPM Flow 2025.10 (open source). ✅ = comparable · 🟡 = partial · ❌ = missing in ResSim.

| Area | CMG | ECLIPSE | OPM Flow | ResSim (own engine + GUI) |
|---|---|---|---|---|
| Black oil | IMEX | E100 | Flow | ✅ E100-class FIM; Norne FPR within 0.3 % of ECLIPSE; Egg, SPE1/9 |
| EOS compositional | GEM | E300 | 🟡 young compositional module | ✅ FIM natural variables, SPE5 within a few % of E300; ❌ SPE3 reference comparison |
| Thermal / steam / combustion | STARS (industry standard) | E300 THERMAL | 🟡 energy equation, CO2STORE thermal | ✅ STARS-like p-N-U core, IF97 incl. supercritical, reactions, ISC, UCG, SPE4 1A in range; ❌ FlexWell-class wellbore, ❌ SPE4 2A/3 |
| CO2 storage | GEM (CCS, geochemistry) | E300 CO2STORE | ✅ CO2STORE, SPE11 participant | ✅ all trapping mechanisms in one unified framework, SPE11A-like; ❌ real SPE11 geometry, ❌ full geochemical speciation |
| Geomechanics | coupled module (IMEX/GEM/STARS) | Visage (external) | 🟡 coupled, under development | ✅ THM / THMC fixed-stress, Terzaghi; ❌ plasticity, ❌ Mandel check |
| Fractures | DUALPOR/DUALPERM, DFN, hydraulic fractures | DUALPORO/PERM | 🟡 no DUALPORO | ✅ dual porosity / permeability, EDFM + Barton–Bandis + hydro-shear, EGS designer; ✅ MINC in black-oil, compositional and thermal (C14/C17), ✅ PKN / KGD propagation + multi-cluster stress shadowing (C14/C17), ✅ HFRAC on corner-point hosts (C15); ❌ inclined EDFM, ❌ growth inside the coupled run |
| Wells | segmented wells, FlexWell | MSW, VFP, groups, ACTIONX | MSW, ACTIONX, UDQ, groups | 🟡 groups + VFP/THP ✅; ❌ multisegment wells, ❌ ACTIONX/UDQ, WECON/WTEST in the own engine |
| Performance | GPU (IMEX-GPU, GEM), MPI, DynaGrid | MPI; INTERSECT for giants | MPI, CPR-AMG, GPU in progress | 🟡 OpenMP CPR-AMG; Egg 28 s vs JutulDarcy 8.6 s; ❌ MPI, ❌ GPU solver in production, ❌ million-cell scale |
| History matching / optimisation | CMOST (proxy, Bayesian, AI-assisted, cloud) | Petrel / MEPO (external) | ❌ (external tools such as ERT / Everest) | ✅ **built in**: tornado, MC, DoE + proxy, LM, ES-MDA, **adjoint gradients (scalar and per-cell fields)**; ❌ production optimisation, ❌ multi-objective GA, ❌ parallel batch runs on OPM |
| Pre / post-processing | Builder, Results | Petrel, Office | ResInsight (separate) | ✅ one application: deck editor (two-way ECLIPSE/OPM), all editors, 3D viewer for own / ECL / OPM output, streamlines, vector PDF/PNG export |
| Restart, binary output | ✅ | ✅ | ✅ | ❌ restart; 🟡 reads EGRID/UNRST, writes its own format |
| Validation and QA | decades, large user base | decades | open benchmark suite | 🟡 331+ automated tests with analytic / SPE references; no field track record |

**Where ResSim is ahead or unusual:**
1. **Adjoint history matching inside the GUI.** CMOST, Petrel and ERT are gradient-free (proxy / ensemble). ResSim has exact discrete-adjoint gradients for scalar parameters and per-cell fields, with an iteration / field-map / response viewer.
2. **One unified multi-component framework.** Thermal, compositional and CO2 trapping share one mechanism (for example UCG + CO2 trapping in the same run). CMG splits this across GEM and STARS.
3. **Integrated ECLIPSE / OPM workflow.** The same application edits decks two-way, runs ResSim, OPM Flow or ECLIPSE, and views all their results.
4. **Openness.** Equations and algorithms are documented (ENGINES_AND_MATH.md) and every verification is reproducible from the test suite.

**Most important gaps, in priority order:**
1. Multisegment wells and ACTIONX / UDQ / WECON in the own engine.
2. Restart plus binary EGRID / UNRST writing.
3. Speed: mixed-precision / multicolour CPR, then MPI or GPU, then million-cell models.
4. MINC and fracture propagation.
5. Production optimisation (adjoint-based and proxy-based) next to HM.
6. Plasticity, the Mandel check, and a geochemical speciation library.
7. Validation against field cases (Norne history match) and more SPE problems (SPE3, SPE4 2A/3, SPE10 full).

## Round C14: restart + result files, custom plots, single-plot export, MINC, fracture propagation

Scope requested: own binary restart plus result files; user-defined saved plot selections viewable after the run;
single-plot export on every plot; MINC, fracture propagation and other unsupported types; MPI (MS-MPI) and a GPU solver
via sokol_gfx; SPE3, SPE11 true geometry and an EGS paper verification; Cartesian + corner-point everywhere;
geomechanics enhancement (continuum damage / XFEM - assistant's choice).

| Item | Status | What was done |
|---|---|---|
| Binary restart (own format) | ✅ | black oil, compositional and thermal cores write the complete dynamic state (`stateIO`: primary unknowns, lagged well / hysteresis / sorption / fracture / geomechanics state, cumulatives, time-step controller) at chosen report steps; *Run → Restart from…* window. Restart = schedule replay + state load. `-rst`: every example restarts **bitwise identical** to the continuous run (black oil, comp wells, thermal incl. SAGD + THM geomechanics) |
| ECLIPSE-compatible output | ✅ | UNRST / UNSMRY / SMSPEC / INIT / EGRID writer (`ecl_write.cpp`) for the own engine; the viewer and OPM / ResInsight read them |
| Result files | ✅ | `.rsres` (complete run results incl. summary, maps, wells, components) save / open from the Results workspace |
| Custom saved plots | ✅ | Results → **Custom plots**: build plots from any vector (field / well / group / component / geothermal / block / derived), name and save them into a library (`custom_plots.txt` in the settings folder, kept across sessions), per-plot style; plots marked ★ get their own tab after *Convergence* and open for any run |
| Single-plot export on every plot | ✅ | right-click on any plot → **Export this plot (PDF / SVG / PNG / CSV)…** - that plot alone, through the figure exporter; installed globally through an ImPlot capture hook so every plot of the app (Results, Studies, Geothermal, components, EGS, adjoint HM, stimulation) has it |
| MINC | ✅ | `-- RESSIM: MINC <shells> [dims] [fractions]` (DUALPORO decks) and the dual-porosity editor (shell count, fracture sets, volume fractions, cells added). Nested matrix shells with Pruess–Narasimhan transmissibilities; default fractions ∝ r^j, r = 1 + 3/n (thin outer shell). `-minc`: slab drainage vs the analytic series - MINC10 max error 0.013, MINC4 0.026, Warren–Root 0.373. The deck checker explains the RESSIM-only status (ECLIPSE / OPM run Warren–Root) |
| Hydraulic fracture propagation | ✅ phase 1 | new module `hydrofrac` - PKN and KGD (displacement-discontinuity elasticity), lubrication flow, Carter leak-off, toughness (calibrated tip-element K_I), implicit Newton on a fixed grid. `-hf`: PKN vs Nordgren −1.7 % (L) / +1.7 % (w0); KGD M-vertex +4.0 % converging; KGD K-vertex −0.4 %; Carter limit +2.8 %; exact mass balance. **EGS designer → Stimulation sizing**: each stage grown from the injector (half-length + retained aperture), preview with L(t) and width profile, warning when a stage does not reach the producer |
| Restart self-review fix | ✅ | the thermal `loadState` refreshed the derived (T, p, S) arrays; a continuous run keeps those of its last Newton evaluation, which seed the next flash and the mechanics → SAGD-THM restarts drifted (1.4e-2). Now restored exactly |
| Tests | ✅ | native suite **360 PASS** (was 335 in C13); Windows build under wine 360 PASS + 3 file-deletion failures of the test itself (fixed: the reader is closed before the temporary files are removed; `-rst` 16/16 PASS under wine); MinGW build without warnings |
| Fixes found in this round | ✅ | thermal restart derived-state (above); tab marker ★ was not in the UI font (now the icon-font star); the custom-plot help mentioned a non-existent corner button |
| Fracture propagation phase 2 (in-simulator growth coupled to the thermal / EDFM run, stress shadowing, P3D / planar 3D) | ❌ next | |
| Geomechanics: Mohr–Coulomb / Drucker–Prager return mapping + continuum damage with permeability coupling (chosen over XFEM: works on the existing hexahedral FE grid, corner-point hosts, no enrichment bookkeeping); Mandel + elastoplastic wellbore checks | ❌ next | design fixed (see ENGINES_AND_MATH §15) |
| MPI phase 1 (MS-MPI / MPICH ABI loaded at run time, partition + halo, `-n 1/2/4` checks) | ❌ next | |
| GPU phase 1 (sokol_gfx compute, float32 + defect correction) | ❌ next | |
| SPE3, SPE11 true geometry, EGS paper layout (Pruess 2006 / GTO-CCS) | ❌ next | |
| Mixed-precision BILU0, Mandel, SPE4 2A | ❌ next | |
| Corner-point audit of the C12–C14 features | 🟡 | black-oil MINC / restart / result files work on corner-point grids; EDFM and the EGS designer are Cartesian-host only |

## Round C15: plot export icon, Excel-style tables, hydraulic fracturing for every reservoir, fault & fracture stability, user guide

| Item | Status | What was done |
|---|---|---|
| Export icon on every plot | ✅ | the icon in the plot corner opens the single-plot exporter (PDF / SVG / PNG / CSV), same as right-click. The press also activated ImPlot's own pan item, which hid the click; the button now reacts to the press itself |
| Excel-style tables | ✅ | one table component (`gridTable`) for layers, schedule, completions, faults, NNC, components / BIC and all editable property tables: click / Shift / drag selection, row, column and whole-table selection from the headers, Ctrl+C (Ctrl+Shift+C with headers) / X / V / A / D (fill down), Del, F2 / Enter / typing to edit, arrow / Tab navigation, context menu (insert / delete rows, fill, clear, replace all). Paste from Excel / LibreOffice / another table: tab, comma, semicolon or space separated, header row skipped, exact cell-to-cell paste from the selected cell, one value fills the whole selection, rows appended when needed |
| Hydraulic fracturing in black-oil / unconventional / compositional runs | ✅ | `HFRAC` planes (EDFM cells in any host, Cartesian or corner-point), multistage designer on the well (spacing, half-length from PKN / KGD, proppant width), LANGMUIR adsorbed gas, PERMSTRESS stress-dependent permeability (matrix + per-fracture). Thermal keeps its FRACTURE / EDFM path. Example **SHALE_HF** |
| Geomechanics: fault & fracture stability | ✅ | Tools → Geomechanics: stress state from gradients, Mohr circles + Coulomb envelope, slip / dilation tendency, CFF and critical pressure per plane (faults, hydraulic fractures, user planes), evolution along the run's pressure history with the poroelastic stress path, stereonet, per-cell map |
| Realistic fracture movement | ✅ | rate-and-state friction slip (quasi-dynamic, radiation damping), stick-slip vs stable creep (k vs k_c), event slip / M_0 / M_w, Barton–Bandis closure + shear dilation aperture and permeability multiplier vs time |
| User guide | ✅ | `docs/USER_GUIDE.md` → `USER_GUIDE.html` (tools/mkguide.py, self-contained), Help → User guide / F1; every feature with its example file; change log updated every round |
| Fixes found in this round | ✅ | the 3-D static view and grid statistics indexed the input arrays with appended (HF / LGR / MINC) cells (out of range); duplicate stereonet labels; the geomechanics window now opens as a tab next to the 3-D viewer (more room for its plots); figure preview showed "?" for the typographic minus |
| Tests | ✅ | native suite **377 PASS** (was 360 in C14) incl. the new `-stab` and `-hfbo`; Windows build under wine **377 PASS**; MinGW build without warnings |
| Still open | ❌ next | MPI (MS-MPI), GPU solver via sokol_gfx, FEM plasticity / damage, SPE3 / SPE11 true geometry / EGS paper layout, mixed-precision BILU0, SPE4 2A, fracture growth coupled inside the flow run (stress shadowing), aquifer-connection and thermal-component tables still use the older table widget |

## Round C16: MPI, GPU, plasticity, all tables Excel-style, remaining verifications

| Item | Status | What was done |
|---|---|---|
| MPI | ✅ batch | MS-MPI / MPICH loaded at run time; `ResSimCLI.exe --mpi-batch` runs Studies case-parallel with progress / cancel in the app (local fallback); *Test MPI* button. ❌ domain decomposition of a single run |
| GPU solver | ✅ phase 1 | D3D11 compute (own device): BSR SpMV, block-Jacobi BiCGSTAB in float, refined in double on the CPU, automatic fallback; *Test GPU*. Untested on a physical GPU (sandbox has none); sokol_gfx kept for the rendering migration |
| Mixed precision | ✅ | float BILU0 (half the memory), double outer solve, double rescue |
| Geomechanics detail | ✅ | Drucker–Prager / Mohr–Coulomb matched plasticity, dilation, softening / hardening to residual cohesion, damage permeability, rigid-plate BCs, EPSPLAS / EPSPV / KDAMAGE outputs, UI + deck lines |
| Tables | ✅ | `RecordTable<T>` (one line per column) — aquifers, thermal components / solids / sorption / water injection / traps, EDFM, minerals, separator, plates now Excel-style; new tables need only a column list |
| Verifications | ✅ | SPE4 2A (near 370 d, far 1,645 d — both inside the participants' range), SPE11A official geometry (coarse), Pruess CO2-EGS (heat ratio 1.46, mass ratio 3.46 vs 1.5 / 3.7), Mandel + Cryer, plastic platen test |
| Fixes found | ✅ | thermal CO2 injection enthalpy without the PR residual (cold CO2 entered warm), LBC dense-CO2 viscosity ×2 → Fenghour correction, produced-steam rate reported per producer |
| Tests | ✅ | native suite **390 PASS** (377 in C15) + flag suites `-mp`, `-gmp`, `-spe4b`, `-spe11g`, `-egs` all pass; Windows build under wine **389 PASS** (the two MPI-run checks are replaced by the "no MPI installed" check — wine has no MS-MPI); MinGW build without warnings; GUI smoke test (thermal THM deck, aquifer RecordTable) |
| Still open | ❌ next | single-run MPI domain decomposition; GPU timing on real hardware and device ILU; SPE3 compositional (EOS data); SPE11 dispersion + fine-grid run, SPE11B/C; implicit steam-rate constraint; MINC inside the thermal deck path (EGS used hand-built MINC); coupled fracture growth with stress shadowing; figure export log axes |

## Round C17: implicit steam limit, MINC in every engine, dispersion, multi-cluster fracturing

| Item | Status | What was done |
|---|---|---|
| Implicit steam-rate constraint | ✅ | producer steam limit as a Newton well control (steam-measure Jacobian, hysteresis band); SPE4 2A holds the limit in 100 % of the reports (C16: 67 %); `STEAMTRAP` item 3 + *Max steam* column |
| MINC in thermal and compositional decks | ✅ | host-cell map for MINC shells / LGR cells in both deck builders; thermal slab (`-tminc`) and compositional pressure-drainage slab (`-cminc`): MINC 8 within 0.006 of the fine model, Warren–Root 0.14 |
| Conduction / diffusion geometry | ✅ | `Connection::geo` (A/d without k, NTG); layered conduction exact (18.18 m); dual-porosity conduction halving fixed; compositional diffusion uses the same geometry; matrix–fracture diffusion added |
| Mechanical dispersion | ✅ | OPM `DISPERC` (per-cell dispersivity) in the compositional / CO2 engine, full Jacobian, deck round trip, UI; `-cdisp` recovers α within 4–7 %; SPE11A now runs with the specified 1e-2 m |
| Multi-cluster fracturing | ✅ | simultaneous PKN clusters sharing the stage rate: Sneddon stress shadow, per-cluster stress offsets, limited-entry perforation friction; designer panel with cluster table and plots; `-hf` checks (Nordgren limits, shadow, limited entry) |
| Figure export | ✅ | log axes in PDF/SVG export, y-label fit |
| Fixes found by code review | ✅ | compositional / thermal deck builders now expand the grid arrays like the black-oil setup — well indices of LGR completions read PERMX out of range, LGR cells used porosity 0.2 / no diffusion / SATNUM 1, DPNUM folding was ignored; DP 'M' conduction halved; conduction T/k_h wrong for unequal half-cells with a permeability contrast; cluster-plot log axis with t = 0 |
| Tests | ✅ | native suite **407 PASS, 0 FAIL** (390 in C16) incl. the new `-cminc`, `-cdisp`, `-dlgr`, multi-cluster `-hf`, `-tminc` checks; `-spe4b` and `-spe11g` (with dispersion) pass; Windows build (MinGW, no warnings) + wine suite **406 PASS, 0 FAIL** (the MPI-run check is replaced by "no MPI installed" under wine) |
| Still open | ❌ next | single-run MPI domain decomposition; GPU device ILU + real-GPU timing; SPE3 compositional (EOS data); SPE11 full dispersion tensor, 1 cm runs and SPE11B/C; fracture growth inside the coupled flow run (EDFM planes growing with K_I), P3D height growth; adjoint for the compositional / thermal cores |

## Round C18: compositional + thermal adjoint, in-run fracture growth, single-run MPI DD, GPU ILU, Intel MPI fix

| Item | Status | What was done |
|---|---|---|
| Intel MPI "error code 5 Access is denied" (user screenshot) | ✅ code / ⏳ user check | libfabric bin dirs prepended to PATH, job object allows breakaway, `RESSIM_MPI_LIB` = launcher's DLL, `ResSimCLI --mpi-info` and an error-5 diagnosis in *Test MPI* (Unblock, antivirus, Intel < 2021.6 → `-localonly` / MS-MPI) |
| Compositional adjoint | ✅ | generic step-adjoint driver (`stepadjoint`), AD-exact residual; porosity gradient within 7.5e-6 of central FD (`-cadj`) |
| Thermal adjoint | ✅ | same driver; FD Jacobian blocks of the steam tables limit it to 1e-4 – 4e-4 relative error (`-cadj`) |
| Fracture growth inside the flow run | ✅ | `HFGROW` (K_IC, σ_min, dσ/dz, initial half-length, closed multiplier); LEFM tip test K_I = p_net √(π min(L, H/2)) after every step opens EDFM pieces; restart-safe; `HFOPEN` results array; Excel-style growth table in *Hydraulic fractures*; example `HFGROW_INJECTION.DATA`; `-hfg` (8 checks) |
| Single-run MPI domain decomposition | ✅ phase 1 | distributed BiCGSTAB with CPR (global AMG pressure stage + block-Jacobi BILU0), halo Sendrecv; `ResSimCLI --mpi-dd`; app *Settings → Parallel → MPI → Single runs: domain decomposition*; SPE9LIKE 2 ranks = serial results (1176 vs 849 linear its) |
| GPU ILU + real timing | ✅ | multicolour BILU0 triangular solves on the device (one dispatch per colour), device-busy timestamps vs wall, readback count, *Benchmark* table and `--gpu-bench`; zero-residual guard |
| Results file v2 | ✅ | size-prefixed RunStats in `.rsres` and in the black-oil restart state; C17 (v1) files still load |
| XFEM question | ✅ answered | see below |
| Tests | ✅ | native suite **433 PASS, 0 FAIL** (407 in C17) incl. `-cadj`, `-mp`, `-mpidd`, `-hfg`; `mpiexec -n 2 … -mpidd` passes; Windows build (MinGW, no warnings) + wine suite **432 PASS, 0 FAIL** (MPI run skipped under wine); app smoke test loads the new example (inputs valid). Review fix: thermal-adjoint FD check retries a forced-step run with a looser Newton tolerance (Windows rounding) |
| Still open | ❌ next | distributed assembly / AMG and DD for compositional + thermal; adjoint through growing fractures; P3D variable-height growth; GPU numbers on a real discrete GPU; SPE3 (EOS data); SPE11 α_L ≠ α_T, 1 cm, B/C |

**C18.1 fix (user report):** the 3D result colours were rescaled at every report step (auto range per step), so colours could not be compared over time. Auto range now defaults to the min / max over **all** report steps of the run (cached per run and property; ECLIPSE / OPM restart steps read once), with *This step only* as an option; the setting is saved.

**Would XFEM on the current grids give better geomechanical property analysis?** Yes for *fracture-scale* quantities, no for the
*bulk* ones. The present hexahedral FE (poro-elastic + Drucker–Prager / Mohr–Coulomb plasticity with damage, fixed-stress coupling)
already gives the stress field, compaction, plastic zones and fault stability (Coulomb / slip tendency) well; XFEM would not change
those. What XFEM adds is a displacement *jump* across the HFRAC/EDFM planes without remeshing: aperture from mechanics (instead of
the Barton–Bandis / exponential laws), shear slip and dilation with friction on the plane, and the tip stress-intensity factor
from the interaction integral (instead of the PKN/KGD estimate used by HFGROW). That is a real gain for stimulation, EGS shear
stimulation and fault reactivation studies. Cost: cut-cell (sub-tetrahedra) integration, tip enrichment, ill-conditioning for
small cut volumes (needs a stable preconditioner), contact (Lagrange/penalty) for closed fractures, and corner-point hexes cut by
arbitrary planes. Plan: (1) Heaviside enrichment of hexes cut by fixed planes, aperture → EDFM conductivity (cubic law) in the
fixed-stress loop, Sneddon / Griffith crack checks; (2) frictional contact (penalty + augmented Lagrangian), slip-induced dilation,
Coulomb checks vs the analytical sliding crack; (3) tip SIF via the interaction integral feeding HFGROW (replacing p_net √(πc));
(4) curving growth along the max-hoop-stress direction (needs EDFM re-cutting at run time). Phase 1 is a round of work; phases 2–4
one round each.

## Round C20: research track, seismicity designer, comp/thermal DD, HFGROW adjoint, P3D, deep geomechanics G1

| Item | Status | What was done |
|---|---|---|
| Research as a standing track | ✅ | `research/`: Cycle 2 report, FINDINGS_LOG F1–F9 (publication-grade rows), literature watch via OpenAlex (`lit/watch.py`, digest), paper P1 (convex seismicity design) manuscript + publication plan (GRL, then JGR Solid Earth) |
| C20-1 Seismicity-aware injection designer | ✅ | `seismicity.cpp` (Dieterich rate, convex optimiser, engine Green's functions, sequential convex correction F8), *Tools → Geomechanics: seismicity-aware injection design*; `-seis` |
| C20-2 Research Cycle 2 | ✅ | H1b poroelastic regimes (F7, consistent with ARMA 2026), peak-rate LP (F4), flowback (F5), allocation (F6) |
| Compositional / thermal MPI DD | ✅ (CLI) | DD hooks in compsim / thermsim; `ResSimCLI --mpi-dd` for any deck; COMP / CO2STORE identical, STEAMFLOOD +0.025 %. App runs comp/thermal locally (CLI results are summary-only) |
| C20-4 HFGROW adjoint | ✅ | frozen-topology discrete adjoint through fracture growth (`AdjRec::hfo`); `-hfadj` 8e-5 vs FD |
| C20-3 P3D | ✅ | stress-contrast height growth, analytic section slopes, height cap + runaway warning, profile plot in the designer; `-p3d` |
| **Deep geomechanics G1** | ✅ 2-D | `efem.cpp`: embedded fractures with continuous Heaviside-enriched jumps, penalty contact, Coulomb friction (return mapping), Newton; Sneddon 0°/30° within 2 %, inclined frictional slip within 5 %, locked case; classic constant-jump EFEM kept for comparison (F9); *Tools → Geomechanics: fracture contact & slip (2-D)*; `-efem` |
| GUI smoke test (Wine) | ✅ | fracture-contact window (solve, plots, stress map → continuous colormap fix) and seismicity designer run in the GUI. **Bug found and fixed:** on HFGROW_INJECTION a correction run hit fracture growth + time-step collapse and aborted the whole design → now keeps the best verified schedule, falls back to constant rate, never returns a schedule worse than constant, warns when the linear-model error > 20 % (`testSeismicityRobust`) |
| Tests | ✅ | native full suite 463 PASS / 0 FAIL (+2 robustness checks in `-seis`, PASS); subsets `-seis -hfadj -p3d -efem` PASS; Windows build (MinGW) OK, ressim_tests.exe rebuilt; Wine full suite not rerun this round |
| Research note | 🔬 | Growing fractures make the rate → stress map strongly nonlinear (46 % linear-model error): the convex design (F1) holds for a fixed fracture state only. Candidate H1f: piecewise-linear / growth-aware Green's functions or adjoint-based SCP with trust region |

### Deep-geomechanics plan (G-series, started C20)

| Step | Content | Status |
|---|---|---|
| G1 | 2-D embedded fracture contact mechanics (EFEM/XFEM-type jumps, friction, pressure loading), UI | ✅ C20 |
| G1-c | Couple to the flow run: EDFM fracture pressure → contact problem, aperture / slip → fracture transmissibility (Barton–Bandis + dilation), staggered then monolithic | ✅ C21 staggered (one-step lag, PKN height correction, UI + `--geocouple`); iterative (saveState/loadState) and poroelastic matrix coupling next |
| G1b | 3-D on hexahedra (planar fractures cutting the Q8 geomech mesh), shared fracture object with EDFM | **Deferred to C24** (design unchanged: trilinear hex + nodal Heaviside enrichment as G1, planar polygon cut by marching tetrahedra, same penalty/return-mapping kernel). Reason: a 3-D enriched contact solve needs a 3-D AMG-preconditioned Newton that we cannot verify on this 2-core / 2 GB machine at useful sizes (≥ 40³ hexes ≈ 200k dofs); doing it half-verified would break the "industry grade" rule. G1-c (2-D slice + PKN height) remains the coupled path |
| G2 | DDM / BEM hydraulic fracture (planar 3-D, tip asymptote / ILSA front) as the high-fidelity HF model | **Deferred to C24–C25.** Reason: a credible ILSA implementation (Peirce & Detournay 2008) is a multi-round project (dense hypersingular kernel + H-matrix, level-set front, multiscale tip asymptote, toughness/viscosity regimes) and must be verified against the radial/KGD/PKN self-similar solutions; G4 (phase field) now covers the "no prescribed path" near-wellbore case, P3D covers the field case |
| G3 | Rate-and-state friction on the embedded faults (quasi-dynamic, radiation damping) → event catalogues from the coupled run | ✅ C23: `rsf.cpp` aging law + radiation damping, plane-strain mode-II dislocation kernel, Dormand–Prince RK45, injection (Theis line source) and G1-c pressure tables, initial stress from the EFEM contact solve; spring-slider k/kc, finite-fault cycles (h*, Lb), injection clock advance; UI tab + catalogue; `-c23` |
| G4 | Phase-field fracture sub-model (near-wellbore / lab) | ✅ experimental C23: `phasefield.cpp` AT2 + Amor split (hybrid), pressure work (Bourdin–Chukwudozie–Yoshioka), Q4 + Wilson incompatible modes, AMG; Sneddon volume and Griffith p_c converge like O(ℓ/a), l→0 extrapolation 0.96 / 0.975 (`-c23`). Limits: 2-D, horizontal initial crack, pressure control (unstable growth not followed), no fluid flow inside the crack. Next: volume control (injection), arbitrary initial crack, coupling to flow |
| G5 | Adjoint injection design with contact (extends the seismicity designer from Coulomb proxies to simulated slip) | **Moved to the research track as H7 / paper P2** (research/HYPOTHESES.md). Reason: with G3 in place this is no longer an engineering gap but an open research question (is the rate-and-state design problem still convex / what replaces F1?), so it is run as a pre-registered cycle (C24-r), not a feature |

## Master plan C19+ (research cycles interleaved with engineering)

The research programme lives in `research/` (CHARTER, GAP_MAP, HYPOTHESES, CYCLE*_REPORT, exp/).

**Research Cycle 1 (done):** H1, the convex theory of minimum-seismicity injection under Dieterich rate-and-state with linear poroelastic coupling.
- Result: SUPPORTED. 7 of 7 pre-registered checks passed.
- Finding: the literature contradictions are an artefact of the evaluation window (step-up −25% during injection vs +4% once the post-shut-in tail is counted).
- Sensitivity law: Π = ΔS_max/(Aσ). See research/CYCLE1_REPORT.md.

| Order | Item | Kind | Decision and reason |
|---|---|---|---|
| C20-1 ✅ | **Seismicity-aware injection designer** (H1c): fault patches in the UI; Green's functions at fault cells from the engine (unit pulse or adjoint); convex optimiser; compare and export schedules | new capability | **Do.** Not found in any commercial simulator; builds on the existing adjoint and geomech |
| C20-2 ✅ | Research Cycle 2: H1b (poroelastic kernels, multiple wells, peak-rate / Mmax objectives, replication of Segall–Lu and Alghannam–Juanes) | research | Do |
| C20-3 ✅ | **P3D** height growth for HF (layered σ_min, Mack–Warpinski style) | engineering | **Do.** Biggest HF quality gap; prerequisite for the DDM HF round |
| C20-4 ✅ | **HFGROW adjoint** (treat growth events as frozen topology; FD on KIC/σ for the growth time) | engineering | Do after P3D. Enables HF design and inversion |
| C21-1 ✅ | Research Cycle 3: H1e Otaniemi field replication, H1d magnitude objective, H2 sub-grid tip (refuted: leak-off quantisation dominates; two-sided bracket), H3 laptop performance | research | Done, research/CYCLE3_REPORT.md |
| C21-4 ✅ | Laptop performance: physical-core thread default, polite MPI wait + PASSIVE OpenMP, hybrid DD, lazy CPR, HPC advisor | HPC | Done (user report: no DD/GPU speed-up, high CPU) |
| C21-5 ✅ | G1-c flow ↔ contact coupling | geomech | Done |
| C22-1 ✅ | Research Cycle 4: H1f closed form of the optimal schedule (bang–singular–bang with terminal burst, F12), kernel-class dependence (H1g, F13), Basel open-data check (no open injection log) | research | Done, research/CYCLE4_REPORT.md |
| C22-2 ✅ | Seismicity designer: design potential (Π, θ, first-order reference), shut-in-burst control (last-interval bound), far-field flag | engineering | Done (-c22) |
| C22-3 ✅ | Exact (WYSIWYG) plot export for every plot type: grid, heatmaps, fracture markers, annotations, colour scale; vector + PNG; print colours | UI | Done |
| C23-r ✅ | Research Cycle 5: H6 Äspö HRL open data (second site), H1h shut-in-aware objective | research | Done, research/CYCLE5_REPORT.md (F14–F16; H1h refuted; Dieterich post-shut-in prediction refuted on Äspö) |
| C23-1 ✅ | Native publication export (figure layers: heatmaps, shapes, colour bars, fracture markers, grids drawn natively, not screenshotted) + Exact mode | UI | Done |
| C23-2 ✅ | Stage-kernel calibration wizard (`seis::fitStageKernels`, Seismicity designer) | research → app | Done (-c23) |
| C23-3 ✅ | G3 rate-and-state faults, G4 phase-field fracture (experimental) | geomech | Done (-c23) |
| C24-1 ✅ | Integration: fault mechanics (Coulomb, Dieterich, catalogue, two-way transmissibility, RSF traces) and G1-c contact coupling run inside the main Run; fault panels + events in the main 3D view; Results tab *Faults & seismicity*; deck keywords STRESSFIELD / FAULTMECH / FAULTRSF / GEOCOUPLE | geomech/UI | Done (-c24) |
| C24-2 ✅ | P1 paper draft v0.4 (math with proofs, experimental study, results, discussion; figures rendered by the app figure engine, tools/paper_figs.cpp) | research | Draft; not submission-ready until a positive second site (H8 Basel) |
| C25-1 ✅ | Non-XFEM geomechanics path: isoparametric hexahedra on corner-point grids, sequential FE coupling (pore volume, stress-dependent permeability) for black oil and compositional/CO2, THM for thermal; HFGROW FESTRESS (growth against the FE σ_min), STRESSFIELD FE / TH (fault stress from the FE stress, thermo-elastic term); UI: Grid & rock → Reservoir geomechanics, stress source, FE-stress column | geomech/UI | Done (-c25) |
| C25-2 ✅ | Compositional Newton robustness for tiny EDFM cells (2-cycle damping, relaxed PV-fraction acceptance) | comp | Done |
| C25-t ✅ | Tests: 589 native checks pass (34 new since C24, incl. -c25 BO one/two-way, analytic stress, CP uplift, CO2 growth + fault, EGS thermo-elastic, thermal geo round trip) | tests | Done |
| C25-3 ✅ | Research Cycle 6 (Otaniemi 2020 OTN-2 replication; F10 stays CANDIDATE, θ-law E8′ confirmed) + Paper Pack v0.5 (GRL letter in AGU LaTeX + SI + cover letter; JGR long PDF/markdown; data availability; zip) | research | Done as pre-submission pack; open before submission: independent third site (e.g. Utah FORGE), reference/novelty check by the authors, author list, Zenodo DOI |
| C26-U ✅ | Usability: start page, command palette (Ctrl+Shift+P / Ctrl+K), recent files, workflow checklist, notification history, run history with *what changed*; ResSimCLI `--check --set --json --csv --version` + coupled runs (fault mechanics, FE geomechanics, DDM) from the command line | UI/CLI | Done. Open: `--list-examples`; studies (sensitivity / HM) still run the flow model only (no fault mechanics / geomechanics inside study runs) |
| C26-1 ✅ (C24-a) | DDM with triangular dislocation elements (Nikkhoo & Walter, cutde, half/full space): FAULTSLIP quasi-static slip + stress transfer between fault faces, slip weakening cascades, moment / Mw, slip-corrected Dieterich | geomech | Done (-c26). C28: thermo-elastic stress in FAULTMECH and FAULTSLIP in the thermal (EGS) and compositional (CO2) engines done (headless coupled runner, CLI, studies, MPI KPIs; -c28) |
| C26-2 ✅ (C25-4b) | DDM apertures of in-run hydraulic fractures (HFGROW … DDMWIDTH) → conductance | geomech | Done (-c26), black oil |
| C26-3 ✅ (C25-4a) | Iterative two-way FE coupling (GEOMECH coupling 2): report step repeated from saved state to a pore-volume fixed point; loadState keeps history on rewind | geomech | Done for black oil (GUI + CLI). 150 d report steps: error 0.24 → 0.095 bar; < 0.01 bar difference at 5–50 d steps. Compositional keeps explicit |
| C26-4 ✅ (C24-b) | Phase field: volume control (stable growth), any crack angle, anisotropic in-situ stress (T = pI − S), PCG + reused AMG (3.7× faster) | geomech | Done as the 2-D tool (-c26). **C28: coupled into the field run** (HFGROW … PHASEFIELD: persistent volume-controlled model per fracture, advanced each report step by η·V_inj/H, opens the EDFM pieces; L_PF 32.1 m vs Griffith 35.5 m; -c28pf). **C29: two-way coupling** (TWOWAY: V = C·p_net,flow, mass cap, sub-stepping). Open: fluid flow inside the phase-field crack, 3-D |
| C26-t ✅ | Tests: 623 native checks pass (34 new in -c26: DDM kernels, penny / rectangle / shear self-terms, fault slip run, DDM apertures, deck round trip, iterative coupling, phase field volume control / orientation / stress); MinGW build 0 warnings | tests | Done |
| C26-5 ✅ | HPC: CPR-AMG hierarchy reuse across the time steps of a report step, quality monitor (TUNING field 7, default 10): Egg −3 % (−9 % across report steps, but that broke bitwise restarts), identical results | HPC | Done — small gain, stated honestly |
| C27-1 ✅ | Comparison with CMG (IMEX / GEM / STARS / CMOST / Builder / CoFlow) and TOUGH3 (EOS modules, iTOUGH2): COMPARISON_CMG_TOUGH3.md | audit | Done |
| C27-2 ✅ | Passive tracers in the black-oil engine (TRACER / WTRACER / TBLKF verbatim + RESSIM TRACERPROP decay / adsorption), implicit upwind transport, Tracers results tab, TR_<name> 3-D arrays, restart state | engine/UI | Done (-c27). Open: partitioning tracers, parent → daughter chains, dispersion |
| C27-3 ✅ | WECON / WTEST in the ResSim engine (restart state v4) | engine | Done (-c27) |
| C27-4 ✅ | Target-change time-step control (dsTarget / dpTarget) | engine | Done (-c27): within 2 % of the uncontrolled production, bounded dS per step |
| C27-5 ✅ | iTOUGH2-style inversion statistics after HM (covariance, correlation, 95 % CI, CSS, eigen-analysis, residual statistics, studentised residuals, leverage, dataset contributions) | studies/UI | Done (-c27) |
| C27-6 ✅ | Coupled fault mechanics / FE geomechanics inside study runs (closes the C26-U open item); 5 new responses (expected seismic events, max Mw, max dCFS, max uplift, max failure index) | studies | Done (-c27); the MPI batch runner still runs flow only |
| C27-7 ✅ | Van Genuchten–Mualem SWOF / SGOF generator (SCAL editor); ResSimCLI --list-examples / --example | UI/CLI | Done |
| C27-t ✅ | Tests: 659 native checks pass (36 new in -c27: time-step targets, WECON / WTEST incl. restart + round trip, van Genuchten, HM statistics, coupled study, examples, tracers: balance, maximum principle, adsorption / decay, exact decay, restart, round trip); MinGW build 0 warnings | tests | Done |
| C27-g 🟡 | From the comparison: **done in C28** FOSM / data worth, partitioning tracers + decay chains + dispersion, SGS / kriging, GP (ML) proxy, decline-curve + nodal analysis + VFPPROD generation; **done in C29** EWASG-type salt in black oil, FlexWell-lite wellbore + ICD, surface networks; **deferred** DynaGrid | backlog | Ordered by value in COMPARISON_CMG_TOUGH3.md §4 |
| C28-1 ✅ | Comparison items: FOSM prediction uncertainty + data worth, tracers (partitioning / PITT, dispersion, Bateman chains), kriging / SGS, GP proxy, production analysis (Arps, IPR / VLP, VFPPROD) — all with UI | studies/engine/UI | Done (-c28) |
| C28-t ✅ | Tests: 709 native checks pass (50 new in -c28: FOSM / data worth, tracers PITT / dispersion / chains, kriging / SGS, GP proxy, production analysis, FAULTSLIP thermal + comp, distributed assembly in-process + MPI, contiguous BILU0, phase field in the run, SUBGRID LEAK round trip); MPI -mpidd / -mpida pass; MinGW build | tests | Done |
| C29-1 ✅ | Salt precipitation / dry-out in black oil (EWASG type: NaCl tracer, evaporation into gas, halite, Verma–Pruess k), restart v5, SSALT / KSALT, UI panel, example SALT_DRYOUT | engine/UI | Done (-c29salt) |
| C29-2 ✅ | FlexWell-lite: wellbore friction (Darcy–Weisbach / Haaland, heel → toe) + ICD Δp = Kρ/ρ_cal q² solved implicitly per connection; WDPWB; UI panel; example WELLBORE_ICD | engine/UI | Done (-c29wb) |
| C29-3 ✅ | Extended surface network: NODEPROP / BRANPROP balanced per step with branch VFP tables, node pressure → producer THP limit (well VFP or head); WTHPN; UI panel (+ Create network); example NETWORK | engine/UI | Done (-c29net). Open: implicit network Newton, chokes, compressors |
| C29-4 ✅ | Two-way phase-field coupling (HFGROW … TWOWAY) + UI column; example HFGROW_PF_TWOWAY | geomech | Done (-c29pf) |
| C29-5 ✅ | CPR-AMG on the GPU (host hierarchy, device V-cycle + BILU0 stage; 23 vs 43 its) | HPC | Done (-c29gpu) |
| C29-6 ✅ | Gas-injection well start in EDFM fracture cells: bisection repair toward the static BHP, edge-aware (injector min / producer max over connections) | engine | Done (-c29, full suite) |
| C29-7 ⏸ | **DynaGrid (deferred) — design note.** Coarsen / refine blocks of cells by front indicators (ΔS_g, Δz_CO2, ΔT) every N report steps. Needs: (1) a parent–child map like CARFIN LGR, (2) conservative remap of every state array (p, S, R_s, tracers, salt, geomech displacements: volume-weighted for extensive, pore-volume-weighted for saturations), (3) re-computed transmissibilities incl. NNC / EDFM, (4) well-connection remap, (5) restart state of the current grid. Cost: weeks; benefit mainly for large CO2 / thermal fronts. Static LGR covers the current examples | engine | Deferred |
| C29-t ✅ | Tests: 731 native checks pass (22 new in -c29: salt, wellbore / ICD, network, device CPR-AMG, PF two-way; the C21 sub-grid K_I check now exempts host-cell edges where K_I jumps and the root finder stops on its 0.05 m bracket — off-edge max 3.4 %); MPI -mpidd pass; MinGW build 0 warnings | tests | Done |
| C29-R ✅ | Research Cycle 8: H7b (F22 refuted), H2c (F24), H8b Utah FORGE (F23 inconclusive) — research/CYCLE8_REPORT.md | research | Done |
| C30-1 ✅ | Network oscillation (user report): zero-rate GOR / WCT limits in updateThp + implicit node-pressure balance (Illinois regula falsi, step re-solved, ≤ 12 trials); NETWORK.DATA smooth, FOPT 1.054e7 | engine | Done (-c30net) |
| C30-2 ✅ | HFFIELD: DDM opening + stress / displacement / strain field around in-run fractures (HFWIDTH, HFPNET, HFDS*, HFDSNN, HFDISP*, HFSTRV), deck + UI | geomech | Done (-c30hff) |
| C30-3 ✅ | Phase-field field output (onField): PFDAMAGE / PFDS* / PFDISP / PFWIDTH on the grid, per-step fields kept for the viewer | geomech | Done (-c30pf) |
| C30-4 ✅ | geoseq::postProcess: FE geomechanics of any finished run / imported ECLIPSE-OPM case (worker thread, Results tab + result-viewer button) | geomech/UI | Done (-c30post) |
| C30-5 ✅ | Results → Fields & sections (section maps, fracture face + profiles, phase-field mesh, geomechanics post); fracture panels coloured by any array; 15 colour maps + reverse | UI | Done (-c30cmap) |
| C30-6 ✅ | Industrial gap assessment INDUSTRY_GAP.md (UI, physics, results, engine, verification, customisability) and next-round order | docs | Done |
| C30-t ✅ | Tests: 755 native checks pass (24 new in -c30: network, HF field traction / opening, PF output vs Sneddon, post-processing = in-run coupling + ECLIPSE round trip, colour maps; -c29net branch check now at the implicit balance tolerance); MPI -mpidd pass; MinGW build | tests | Done |
| C31-1 ✅ | UDQ in the ResSim engine (udq.cpp): DEFINE / ASSIGN / UPDATE / UNITS, F / W / G quantities with patterns, set functions, SORTA / SORTD, NaN semantics; ASSIGN-only at report start (restart-safe); FVPR / FVIR (ResSim UDQ field quantities); UDQ vectors in the ECLIPSE summary output | engine | Done (-c31expr) |
| C31-2 ✅ | UDA targets (WCONPROD / WCONINJE / WELTARG items naming a UDQ; import / export / validate; next-step lag) | engine | Done (-c31deck) |
| C31-3 ✅ | ACTIONX in the ResSim engine (moved forward from C32): AND / OR / parentheses, well sets + '?', runs / interval, body through applyScheduleText, state v6 restart; import messages info (PYACTION warning) | engine | Done (-c31deck) |
| C31-4 ✅ | UI: Model → Schedule → Run-time logic editor (UDQ / ACTIONX forms, templates, live checks, deck text, UDA attach); Results → UDQ & actions (vectors, firing markers, log, deck logic) | UI | Done |
| C31-5 ✅ | Well section view (curtain along the trajectory, any array / step / colour map, connections open / shut, hover) + log tracks vs TVD / MD | UI | Done |
| C31-6 ✅ | Example ACTIONX_UDQ (voidage UDA, water-cut shut-ins, pressure-triggered target change) | examples | Done |
| C31-t ✅ | Tests: 784 native checks pass (29 new in -c31: 17 expression / ACTIONX parser checks in -c31expr, 12 deck checks in -c31deck incl. UDA lag, run limits, restart after fired actions, summary output, import messages); MPI -mpidd pass | tests | Done |
| C32-1 ✅ | Multisegment wells in the own engine: WELSEGS (INC / ABS, branches) / COMPSEGS (nearest node rule) / WSEGVALV / WSEGSICD / WSEGAICD, Haaland friction, hydrostatic mixture heads, device drops balanced with the subtree inflow (conductance factor), SHUT subtrees, restart; WSEGAICD items 13 / 14 fixed | engine | Done (-c32msw) |
| C32-2 ✅ | ACTIONX COMPDAT adds connections at run time (connection set + Jacobian pattern rebuilt, well state kept, restart replay) | engine | Done (-c32actx; = schedule COMPDAT to 5e-12) |
| C32-3 ✅ | Geomodelling engine (geomodel.h): LAS / heads / surveys (minimum curvature) / tops, IRAP / ZMAP / XYZ, fault sticks, SEG-Y rev1 IBM / IEEE, seeded auto-tracking + v0 + kz depth conversion + attributes, faulted pillar grid (fault blocks, split ZCORN nodes, FAULTS), zones, well tie, blocking, SIS / TGS / channels, SGS + collocated co-simulation, poro-perm, volumetrics + saturation-height + percentiles, flow-based upscaling (PV conserved), GRDECL, toCase, demo project | geomodel | Done (-c32geo, 28 checks) |
| C32-4 ✅ | Geomodel workspace (Ctrl+5) + standalone ResGeo (ResGeo.exe / --geo): 10-step workflow, map / seismic section / well correlation / grid cross-section views, background jobs, send to the Model workspace; ResSimCLI --geo-demo | UI | Done |
| C32-5 ✅ | MSW UI: Model → Wells → Multisegment well generator (tubing + SICD / AICD / valve per connection), Results → Multisegment wells (profiles, heel / toe history, segment table); example MSW_HORIZONTAL | UI | Done |
| C32-t ✅ | Tests: 831 native checks pass (47 new: 15 -c32msw, 28 -c32geo, 4 -c32actx); MPI -mpidd pass; MinGW build 0 warnings | tests | Done |
| C33-1 ✅ | Scripting API: api::Session (open / set / run, step-wise begin / step / applySchedule incl. new wells), C API (capi.h; libressim.so / ressim.dll), Python package (ctypes, no compiler; numpy arrays), examples quickstart / sweep / compare_opm / steering | API | Done (test_api.py 19/19, -c33) |
| C33-2 ✅ | Report generator: DOCX with native charts, PPTX, self-contained HTML, Markdown (own ZIP / DEFLATE / PNG writers); File → Generate report dialog; ResSimCLI --report | QA/UI | Done (-c33) |
| C33-3 ✅ | Deck regression harness (manifest, relRMS per vector, time history, HTML / DOCX / JUnit), fetch_opm_tests.py (202 opm-tests cases with OPM Flow references), nightly scripts, Run → Deck regression window, ResSimCLI --regress | QA/UI | Done |
| C33-4 ✅ | opm-tests fixes: WEFAC (rates / cumulatives / polymer / tracers), run-time wells (WELSPECS in ACTIONX / scripts) and producer ↔ injector switching, WLIST + full glob, stopped-well BHP, MULTX- / MULTY- / MULTZ- | engine | Done (see USER_GUIDE §7g for the pass count) |
| C33-5 ✅ | Linear-solver threading: OMP_NUM_THREADS honoured, BILU0 level-parallel gate, vector-kernel threshold (5–16× on 2 cores) | HPC | Done |
| C33-t ✅ | Tests: 848 native checks pass (17 new in -c33; 2 older checks adapted to the OMP_NUM_THREADS rule and the solver tolerance); MPI -mpidd pass; MinGW 0 warnings; test_api.py 19/19 | tests | Done |
| C34-1 ✅ | Group control: GUIDERAT potential formula (OPM GuideRateModel), hierarchical subgroup guide rates (GCONPROD 9/10 value / FORM / POTN), most-violated group constraint, VREP / REIN coupled to the current step (NUPCOL-like re-solve) | engine | Done (GRCTRL MSW pass; STW 15–22 %) |
| C34-2 ✅ | Aquifers: AQUNUM / AQUCON numerical aquifers, AQUTAB influence tables, AQUANCON coefficient | engine | Done (CT / NUM / FETKOVICH pass) |
| C34-3 ✅ | KRNUMX/Y/Z directional relative permeability (SATOPTS DIRECT), WELPI, rock compaction (ROCKCOMP / ROCKOPTS / ROCKTAB / ROCK2D(TR) / ROCKWNOD / OVERBURD / ROCKNUM), end-point arrays without ENDSCALE (OPM), EDITNNC / EDITNNCR, OPERATE / OPERATER, WECON CON / +CON / PLUG + COMPLUMP, POTN, secondary WCT; restart state v7 | engine | Done |
| C34-t ✅ | Tests: 855 native checks pass (7 new in -c34), MPI -mpidd pass, MinGW 0 warnings, test_api pass; regression 307 cases: 102 PASS / 147 FAIL / 5 ERROR / 53 SKIP | tests | Done |
| C35-1/2 ✅ | Chemical EOR in the own engine: SURFACT (viscosity, IFT -> capillary number -> SURFCAPD miscible kr blend with SURFNUM, Pc scaling, adsorption), ALKALINE (ALSURFST / ALSURFAD / ALPOLADS / ALKADS), FOAM (FOAMMOB x FOAMMOBP x FOAMMOBS, adsorption, FOAMDCYW / FOAMDCYO decay); EOR editor; FT*/FS*/WT* vectors; restart v8; SPE1FOAM matches OPM (0.8 % relRMS); -c35 tests (15) | engine/UI | Done |
| C35-3 ✅ | SOLVENT in the own engine (four-component black oil, immiscible SSFN split, flash of SgT into Sg / Ss, sequential-implicit composition with outer coupling, exact solvent balance, WSOLVENT, SSOL, FOAMOPTS SOLVENT, FN*/WN* vectors, restart v9); SPE1CASE2_SOLVENT_FOAM matches OPM (FGOR 0.4 %, BHP 0.9 % relRMS); chemSolve line search; EOR editor Solvent tab; -solv tests (8). Miscible model (MISC / TLMIXPAR) and gas-water solvent deferred | engine/UI | Done |
| C35-4 ✅ | Gas-lift optimisation in the own engine (LIFTOPT / WLIFTOPT / GLIFTOPT): lift curves per THP producer from the VFP ALQ axis at the frozen state, greedy concave-hull allocation on the OPM economic gradient (weights, incremental-gas weight, min / max / fixed lift gas, dead wells as one block), GLIFTOPT lift-gas supply / total gas and GCONPROD caps along the group chain, FGLIR / FGLIT / WGLIR / WGLIT, restart v10 (also fixes the THP-derived BHP limit after a restart: THP cases now restart bit for bit), Wells-page editor + allocation plot, custom-plot vectors (gas lift + solvent); opm-tests GASLIFT-01 FOPT -5.3 % and GASLIFT-04 FOPT / FGPT -0.3 / -0.1 % vs the commercial simulator; -glift tests (11). GASLIFT-05..15 (MSW / network gas lift, GLIFTLIM) remain OPM-only | engine/UI | Done |
| C35-4b ✅ | Surface network & facilities with OPM semantics: BRANPROP 0 = drop / 9999 = no loss, GRUPNET standard networks (re-routing via GRUPTREE), NEFAC, leaf inflow = efficiency-weighted subtree rate (+ lift gas, satellites), GSATPROD, GCONSUMP (node / NO-NODES / default leaf), GCONSALE (sales target caps REIN gas), automatic chokes (AS_CHOKE: node pressure bisection to the GCONPROD target), GCONPROD items 11-13 + procedure NONE, robust THP operating point (stable root + golden-section fallback), FGCR / FGIMR / FGSR, restart v10; UI facilities row; opm-tests NETWORK-01 / _STANDARD 0.14 %, REROUTE(_STD) 0.3 %, GSATPROD-GCONSUMP 0.6 %, 9_3C FGSR 0.4 % vs OPM; -facil tests (18: 11 always + 7 against opm-tests, skipped unless fetched). Open: NETBALAN, NETCOMPA compressors, separators, GCONSALE procedures | engine/UI | Done |
| C35 (plan) | Moved from C34 (not started): 1 M-cell performance round; geomodel round 2 (stratigraphic pillars, horizon / fault editing, flattening, uncertainty loop, MPS, out-of-core SEG-Y, 3-D scene); facilities + gas-lift optimisation; solvent / foam / surfactant; transient wellbore; PYACTION; MSW acceleration / drift-flux / WSEGLABY; network re-routing; remaining opm-tests failures (GRCTRL STW, ROCKTAB, AQUFLUX, KRNUM-03, 9_*, WELPI, 3D_WECON, UDT) | engine/UI/HPC/geomodel | Planned |
| C36-1/2 ✅ | Out-of-core SEG-Y (header index, least-squares plan geometry, on-demand slices / sub-cubes, proxy), complex-trace + coherence / variance attributes (slices, cubes, streamed volumes), flattening; 3-D scene renderer (D3D11 slices / meshes / solids / lines, ray-marched volume, time / depth via v0 + kz, DirectWrite labels, PNG) | geomodel/UI | Done (-c36seis 11/11) |
| C36-3 ✅ | Geomodel 3D view: seismic slices / volume / attributes, horizons, faults, wells, grid with filters; seed and fault-stick picking on sections and in 3-D; headless UI harness | UI | Done (harness 30/30) |
| C36-4 ✅ | Horizon brush (map + 3-D) and surface tools with undo; fault-parallel sloped pillars (COORD from fault dip, fixed-point node depths and cellOf); Direct Sampling MPS with servosystem and GSLIB TIs; uncertainty loop (LHS, spectral structural field conditioned to well tops, triangular inputs, tornado, exact realisation rebuild, MPS facies) with UI | geomodel/UI | Done (-c36geo 11/11, harness 35/35) |
| C36-5 ✅ | MSW: WSEGLABY labyrinth ICD (loss-coefficient model from the items; OPM does not simulate it), acceleration (HFA, momentum flux balance), drift flux (DF, Shi et al. 2005 + Hasan-Kabir, WSEGDFPA), holdup / acceleration in Results and the generator; fix: MSW BHP reference depth = WELSEGS top segment (OPM updateRefDepth) | engine/UI | Done (-c36msw 9/9, -c32msw unchanged) |
| C36-6 ✅ | Network: NETBALAN (interval, node tolerance, max iterations, max balance error -> step cut), NETCOMPA automatic compressors (group-target switching by sequence, ON-OFF / MULTI levels, VFP / ALQ, fuel gas in FGCR), COMPOFF; restart state v11; fix: gas-lift group caps ignore limits with exceed procedure NONE (OPM); regression dialog Cancel | engine/UI | Done (-c36net 7/7: compressor +10.3 % FOPT on NETWORK.DATA, MULTI fuel 300/600/900, NETBALAN 365 d holds node pressures, restart exact; -glift, -c33, -c30net, -facil pass) |
| C36-7 | Deferred to C37 (with reason): PYACTION needs an embedded Python interpreter inside the engine (the C API / ctypes direction exists, the reverse does not); the 1 M-cell round and the remaining opm-tests failures need longer measurement runs than this round allowed; transient wellbore (segment storage) is planned on top of the C36 MSW model; WSEGPULL / WSEGFLIM, oil-water slip and tubing flash are still open | engine | Next |
| C37-1 | MSW tubing flash (surface-rate re-split at segment pressure, Rs/Rv), oil-water slip (WSEGDFPA 7-11, Shi 2005), segment storage (swept-volume mixing, restart state), MSWOPTS RESSIM line, generator + Results-tab additions (`-c37msw`) | engine + UI | Done |
| C37-2 | WSEGFLIM (smooth conductance limiter, lagged potential with growth extrapolation) and WSEGPULL (downhole pump boost) - assumed semantics documented (ENGINES 6.36); not simulated by OPM | engine + UI | Done |
| C37-3 | PYACTION in the ResSim engine: persistent Python child over line-JSON pipes, opm_embedded shim (summary on-demand queries, schedule open/shut/insert_keywords, legacy run() + actionx_callback), SINGLE / UNLIMITED / FIRST_TRUE, restart v12, PYACTION editor with templates and syntax check (`-c37py`) | engine + UI | Done |
| C37-4 | 1 M-cell performance round: gen_bigcase.py benchmark, memory probe + gprof, BILU0 contiguous mode without the LU value copy, WellSchur pattern reuse (hash), grid connection trim (100 k: 454 -> 405 MB, precond 6.7 -> 5.0 s, identical results), 9 k - 200 k scaling, 1 M extrapolation (~4 GB, ~14 min on 2 cores); float32 CPR BILU0 measured (not faster) | engine | Done |
| C37-5 | opm-tests triage: WCYCLE (OPM semantics, step landing, start-up ramp, WELOPEN / WTEST re-open, restart v13), WTEST 'P' physical closures, GRIDUNIT (CM / FEET); wcycle 9/9 pass (was 0/9), gridunit geometry exact (remaining 8 % = OPM time-step discretisation after the BHP switch) (`-c37opm`) | engine | Done |
| C38-1 | Level-parallel float32 BILU0 in the CPR second stage (`RESSIM TUNING ... 1`, > ~20 k cells): factor memory halved, 100 k cells 58.0 -> 50.8 s with identical iterations / FOPT (`-c38ilu`) | engine + UI | Done |
| C38-2 | opm-tests: NONNC (PASS 3.05 %), UDT / UDTDIMS with LC / NV / LL and union operators + OPM UDQ evaluation order + UDT editor (`-c38udt`), EDIT-section TRANX / TRANY / TRANZ on neighbour connections only, SWATINIT < SWL raised to SWL, OPM Newton minimum (one update per step) and TSFMAX 3 (`-c38nonnc`, `-c38tran`); model2 + editnnc + udt manifest 13 / 59 pass (C34: 12), no case worse | engine + UI | Done |
| C38-3 | Production-data import: CSV / TSV / paste, column roles, QC, volumes -> rates with uptime (WEFAC), inclusive end dates, unit plausibility (GOR), WCONHIST / WCONINJH merge at DATES, history include, observed CSV; Schedule screen + `ResSimCLI --prod-import` (`-c38prod`) | engine + UI | Done |
| C38-4 | Deferred to C39 (with reason): dropping the scalar Jacobian copy for 1 M cells < 3 GB (needs an assembly refactor larger than this round), 9_* group-control cluster (~30 cases, common ~4 % base offset plus GCONPROD dynamics), NETWORK-01, WPWE, CO2STORE, KRNUM, CECON, ROCKTAB, GRCTRL | engine | Next |
| C39-1 | Cross-platform: Linux / macOS builds from the same sources. Windows keeps D3D11 / D2D / DWrite; Linux / macOS use GLFW + OpenGL 3.3 core back ends for Viewport3D and GeoScene3D (shared CPU geometry, one overlay painter for D2D and ImGui draw lists), platform layer (dialogs via tinyfiledialogs, config dir, process id, message boxes), embedded resources (.incbin), plutovg PNG export, GPU-solver stub, CMake for all platforms, build_linux.sh, CI (Ubuntu 24.04 / macOS 14 / MinGW), desktop entry, `--smoke-test` (verified: Mesa llvmpipe, Model 3-D view + Geomodel seismic scene + light theme) | platform | Done |
| C39-2 | UI / workflow review (`docs/UI_REVIEW.md`: verdict "functionally yes, as a product not yet", gaps G1-G13, modernisation plan, Petrel-parity plan) and first fixes: G1 undo / redo of every model edit (whole-case snapshots, memory-budgeted depth), G2 autosave + crash recovery, G6 dark / light / high-contrast themes with 5 accents, theme-aware labels / status colours / primary buttons | UI | Done |
| C39-3 | Next (UI_REVIEW order): G3 feature registry + command palette for every action, G4 object tree (wells / groups / regions / runs), geomodel log calculator + cross-plots, interactive variogram fitting, G7 wizards / empty states, G9 linked time cursor, well tie, geomodel parameters in Studies, G10 case manager. Engine carry-over: scalar Jacobian copy (1 M cells < 3 GB), 9_* group control, NETWORK-01, WPWE, CO2STORE, KRNUM, CECON, ROCKTAB, GRCTRL | UI / geomodel / engine | Done in C40 (see C40 rows) |
| C40-1 | UI gaps from `docs/UI_REVIEW.md`: **G3** navigation layer (`nav.h`: every tab is deep-linkable and indexed; the command palette lists all tabs, Geomodel steps / views, wells; recents), **G4** object tree in the Case Explorer (field › groups › wells with rename / duplicate / delete / focus in 3-D, faults, LGRs, aquifers, regions, runs, observed data), **G7** new-model wizard (reservoir → fluid / rock-fluid → well pattern → schedule, live pattern preview, runnable case) + empty states, **G9** linked report-time cursor on every time plot (drag = report step of all views; implot PreEnd hook; not exported) | UI | Done |
| C40-2 | Geomodel **Data analysis** view (`geo_analysis.h/.cpp`, `ui_geoanalysis.cpp`): log calculator (expression parser with null propagation, 10 petrophysical presets), cross-plots (colour by well / curve, zone filter, linear / log-linear regression → poro-perm transform of the Petrophysics step), histograms with P10/P50/P90, **interactive variogram fitting** (draggable nugget / sill / range, Cressie-weighted least squares, stratigraphic depth, apply to Facies / Petrophysics; same fit button in the Model geostatistics panel), **well tie** (sonic + density → AI → RC * Ricker, bulk-shift correlation search, time–depth plot, calibrated v0 + k z velocity model). Demo LAS files gain DT / RHOB consistent with the demo SEG-Y. Tests: `-c40geo` (7 checks; demo tie correlation 0.95) | geomodel | Done |
| C40-3 | **Geomodel parameters in Studies**: `Case::propModels` (`PropModelInput`: property, kriging / SGS, variogram, seed, hard data, poro-perm transform) written by Geomodel *Send to Model* and the Model geostatistics panel, persisted as `-- RESSIM: PROPMODEL / PROPDATA`; new `ParamKind::GeoRange` (range multiplier) and `GeoSeed` (realisation) re-simulate the model(s) before the other parameters act (`geostat::applyPropModel`). Usable in sensitivity, ensembles and history matching (seed as a discrete variable). Facies are not re-simulated (single-population porosity recipe). Tests: `-c40geo` | studies / geomodel | Done |
| C40-4 | UI gaps continued: **G10 case manager** (`ui_cases.cpp`: variants of the model with base flag, notes, diffs vs base / current model, run linking, load & run, compare in the Result Viewer; stored in `<project>.cases/`), **G11 Inspector filter** (Ctrl+F; all `ui::prop*` rows), **G12** status icons, **View › Appearance** (theme / accent / 3-D background incl. light and white backgrounds for the Model / Results 3-D views) | UI | Done |
| C40-V | C40 verification: native test suite **1004 PASS / 0 FAIL** (`-c40geo` 10, `-c32geo` 28, `-c36geo` 11); MinGW release build 0 warnings; Linux smoke screenshots of the wizard, variogram fit, well tie, time cursor, case manager, Inspector filter and light 3-D background (`docs/img/c40_*.jpg`) | QA | Done |
| C40-5 | Still open after C40 (next round): fault-network modeller, workflow recorder, G5 rest, G13 i18n; deferred engine items (scalar Jacobian copy for 1M cells < 3 GB, 9_* group-control cluster, NETWORK-01, WPWE, CO2STORE, KRNUM, CECON, ROCKTAB, GRCTRL) | mixed | Done in C41 except the items listed in C41-5 |
| C41-1 | **Workflow recorder** (Petrel 6): record palette commands, wait / note / for-each-variant steps, replay, `.rswf` files (`workflow.h`, `ui_workflow.cpp`); **Inspector "show only modified"** (baseline per row, right-click = new baseline); **G5 design tokens** (surface levels, spacing scale, card headers, primary / secondary / danger buttons); **G13 i18n** (`TR` / `TRN`, `lang/<code>.txt` catalogues, View › Language, Spanish, `tools/i18n_extract.py`; main menus, toolbar and workspace switch translated) | UI | Done |
| C41-2 | **Fault-network modeller** (Petrel 5, `geomodel_faultnet.cpp`): X / Y junction detection with snap tolerance, truncation of minor faults on the major (9 depth levels), throw profiles with reverse / tip-taper QC, revert; **facies by truncated Gaussian simulation** (`geostat::applyFaciesModel`, PROPMODEL 0..3 + `-- RESSIM: PROPFAC`), re-simulated by the Studies GeoRange / GeoSeed parameters | geomodel / studies | Done |
| C41-3 | Engine: **producer and injector cross-flow** (WELSPECS item 10; OPM StandardWell mixture volume ratio for producers, back-production at cell mobilities for injectors with the net injected phase in the rate control) – fixes the ROCK / ROCKTAB / ROCK_REGIONS_DEFAULT opm-tests (WWPR:PROD-3 42.5 % → ~5 %); **ROCKCOMP HYSTER + ROCKTABH** (deflation curve + OPM LeftExtreme guided interpolation of the elastic curves; ROCKTABH FPR 283.93 vs 283.91 bar); **CECON** connection economic limits (WCT / GOR / WGR per completion, CON / +CON / WELL; CECON-02 / -03 rates and shut times match OPM); guide rate unscaled when the current phase rate is zero (OPM GuideRate::get); regression report **material worst** metric | engine | Done |
| C41-V | C41 verification: native suite 1011 tests – 1009 PASS in the full run, the 2 failures found there (gas-lift restart: cross-flow mixture not in the restart state → state version 14; HFA acceleration 6.1 psi: well-level mixture in multisegment wells → cross-flow off for MSW wells) fixed and re-run PASS (`-glift`, `-c36msw`, `-c35`, `-c41` 8 / 8); MinGW release build 0 warnings | QA | Done |
| C41-5 | Still open: segment-wise cross-flow mixture for multisegment wells (cross-flow is off for MSW wells); scalar Jacobian copy (1 M cells < 3 GB); 9_* cluster (strict metric kept: 9_1A fails only on the tiny PROD4 water rate), 9_3E WWPR:PROD2; GUIDERAT_*; NETWORK-01; KRNUM-02Z / 03*; 6_ENDSCALE 9.3 %; WPWE summary vectors (WPWE0-7); CO2STORE decks; GRCTRL STW 15–22 %; CECON PLUG and minimum-rate items (not in OPM either); i18n coverage beyond menus / toolbar; native smoke screenshots of the C41 UI | mixed | C42: CECON PLUG + minimum rates done (ResSim convention, not verified against ECLIPSE); smoke screenshots of the geomodel UI done (`docs/img/c42_*.jpg`); the rest moves to C42-4 |
| C42-1 | Geomodel import with column assignment: sniffing, synonym + positional role guess, manual roles and unit choice per column, several files at once, multi-well surveys (MD/INC/AZI, X/Y/TVDSS, DX/DY), Petrel tops, name matching by column / header / file name | geomodel | Done (18 tests) |
| C42-2 | Well manager: manual wells, per-well trajectory editor (MD/INC/AZI or MD/X/Y/TVDSS, paste, file into one well, vertical to TD), tops editor; LAS batch | geomodel | Done |
| C42-3 | 2-D seismic: several SEG-Y lines, 2-D / 3-D detection, byte locations, section view, tracking / manual picks / erase, intersections and mis-ties, gridding to depth, map + 3-D curtains | geomodel | Done (mis-tie *correction* open) |
| C42-P | Petrel parity (`docs/PETREL_PARITY.md`): surface calculator, property calculator (new / overwrite), 3-D value filter, statistics + histogram; compact (i) help with the *Help text* setting | geomodel / UI | Done; plan items P1–P10 open (object tree with settings, templates, intersection window, polygons, well section window, mis-tie correction, shared filters, layer-cake velocity, seismic calculator, Python) |
| C42-4 | Still open: MSW segment cross-flow mixture; scalar Jacobian copy (1 M cells < 3 GB); 9_* cluster, 9_3E; GUIDERAT_*; NETWORK-01; KRNUM-02Z / 03*; 6_ENDSCALE 9.3 %; WPWE vectors; CO2STORE decks; GRCTRL STW; i18n beyond menus / toolbar; Metal back end for macOS | mixed | C43: WPWE0–7 done (C43-E); i18n extended to the Geomodel workspace; the rest moves to C43-4 |
| C43-G | Petrel parity P1–P10 (`docs/PETREL_PARITY.md` §5): object tree, templates, shared filters, intersection view, polygons (clip, volumetrics limit, contacts per block), well section (tracks, flatten, top picking), mis-tie correction, layer-cake velocity model + checkshots, seismic calculator / crop / SEG-Y export, Python / C geomodel API (`rs_geo_*`, `ressim.geo`) | geomodel | Done (`-c43`) |
| C43-E | WPWE0–7 well performance events (OPM PR 7395 definitions): snapshot tracker per accepted step, '+CON' sweep excluded, closed-to-bottom rule, simulator status events; SUMMARY / UDQ access; WTEST 'E' re-opens economically closed connections | engine | Done (`-c43`); not compared with opm-tests `wpwe/` reference output |
| C43-U | Geomodel workspace translated (steps, views, toolbars, Data panel, object tree; 314 catalogue entries), `--lang <code>` start option; well-section auto-fit and track widths | UI | Done |
| C43-4 | Still open: MSW segment cross-flow mixture; scalar Jacobian copy (1 M cells < 3 GB); 9_* cluster, 9_3E; GUIDERAT_*; NETWORK-01; KRNUM-02Z / 03*; 6_ENDSCALE 9.3 %; CO2STORE decks; GRCTRL STW; WPWE reference comparison (opm-tests `wpwe/`), WELSPECS item 9 STOP for economic closures; translation of help texts and of the model / results panels; faults in the intersection view; 2-D fault sticks on lines; Metal back end for macOS (cannot be built or tested here) | mixed | C44: partly done (rows C44-*) |
| C44-E | WPWE against the OPM reference output of opm-tests `wpwe/` (5 decks identical, `tools/wpwe_compare.cpp`); WELSPECS item 9 STOP (WPWE4); WELOPEN on one layer of a multi-layer COMPDAT record | engine | Done |
| C44-W | Injector cross-flow recirculation (back-produced fluid leaves through the injecting connections, not reported as production); first guide rate from the potential for non-flowing wells; WPIMULT per connection and step (cumulative, reset by COMPDAT, export per layer, completion table column); FBHPDEF; WELPI after the step controls (injector PI of the injected phase); no SWCR/SGCR >= SWL/SGL clamp; EHYSTR 1/3 wetting imbibition; restart state v15 (economic connection closures, WPWE baseline) | engine | Done (ENGINES_AND_MATH §9.w) |
| C44-G | Geomodel: faults cut on intersection paths and 2-D lines, well section width / template / order | geomodel | Done |
| C44-4 | Still open: MSW segment cross-flow mixture; scalar Jacobian copy (1 M cells < 3 GB); KRNUM-03Z (vertical gas); NETWORK-01-WTEST WBHP; 3_A_MPI_MULTFLT_SCHED; 6_ENDSCALE 4.7 % material; CO2STORE; 9_* late FOPR drift; translation of help texts and model / results panels; object-tree folders; Metal back end (cannot be built here) | mixed | Open |
| C45-M | Memory of large models: block assembly without the scalar Jacobian (bit-identical), lean BILU0 (float32 factorised in place, auto ≥ 500k cells, double rescue), AMG takes the fine matrix, accD_ release, solvent arrays only with SOLVENT, Connection 56 B; 100k 416 → 322 MB (lean 294), 200k lean 584 MB → 1 M ≈ 2.75 GiB; UI *Low-memory linear solver* | HPC | Done (1 M cells < 3 GB by extrapolation; the 1 M run needs > 2 GB RAM) |
| C45-X | MSW segment cross-flow: per-segment mixture (OPM MultisegmentWell), storage swept by the injected volume, top-mixture first fill, restart v16 | engine | Experimental, **off by default** (MSWOPTS third value): explicit heads flip-flop on MSW-2D-VERT-02 (FGOR 3.6 → 23.8 %); needs implicit segment coupling |
| C45-U | Translation of the Model / Results panels, widgets and help tooltips (2 011 keys, es 1 660 + 37 kept English), `TRS` / `TRL` | UI | Done |
| C45-4 | Moved to C46-4 (moved from C44-4): MSW cross-flow with implicit segment coupling (default on only after that); KRNUM-03Z; NETWORK-01-WTEST WBHP; 3_A_MPI_MULTFLT_SCHED; 6_ENDSCALE 4.7 % material; CO2STORE decks; 9_* late FOPR drift; GRCTRL STW; 0B_WELPI remainder; object-tree folders; translation of plot legends / long keyword docs; Metal back end (cannot be built or tested here) | mixed | Open |
| C46-G | Geomodel project save / open (`.rgeo`, SEG-Y re-read), object-tree folders (move / rename / remove / drag and drop), in-app round-trip smoke test `gt-proj` | geomodel | Done |
| C46-W | WELPI with Rs / Rv cross terms, WPI* vectors, WELOPEN / WPIMULT completion numbers + COMPLUMP, WPIMULT on SHUT connections, WELPI reset by COMPDAT, ACTIONX WPIMULT (54.1 → 8.0 %), 0B_WELPI 11.7 → 8.1 % (remainder: a transient spike in the reference only) | engine | Done |
| C46-F | MULTFLT in SCHEDULE time dependent and cumulative (3_A_MPI_MULTFLT_SCHED now follows the reference) | engine | Done |
| C46-4 | Still open: MSW implicit segment coupling; 6_ENDSCALE krow front (finding in ENGINES_AND_MATH §9.u); KRNUM-03Z; NETWORK-01-WTEST; CO2STORE decks; 9_* drift; GRCTRL STW; plot legend / keyword-doc translation; Metal back end; ACTIONX_GCONINJE / ACTIONX_WCONINJE time-step failures; full 209-case OPM rerun (the C46 rerun covered model2 + actionx, 67 cases: 27 pass, no regressions against C44 - `docs/regression_C46_model2_actionx.html`; full suite 1096 PASS) | mixed | Open |
| C47-P | OPM low-pressure padding of PVDG / PVTG (opm-common PR 3779), `gasPaddingRows()`; neutral on 9_3A / 9_4B / 9_1A / 6_ENDSCALE / 0A1 | engine | Done |
| C47-D | OPM reference drift study: 18 cases' reference histories, bisection (T1L → opm-simulators PR 6311, 9_3A → opm-common PR 3779), 9_3A reference-unstable (45 % FWCT between OPM versions), `tools/dev/refdrift*`, `smrydiff.py`; ENGINES_AND_MATH §9.t | regression | Done |
| C47-A | ACTIONX: defaulted reference depth of wells completed at run time (ACTIONX_GCONINJE / ACTIONX_WCONINJE: ERROR → FAIL 10.3 / 14.3 %, remainder = reference spike + OP02 family offset), month JLY in conditions | engine | Done |
| C47-R | Research Cycle 9: H10 / H10c SUPPORTED (F25), H10d gap decomposition (F26), H10b SCP SUPPORTED (F27), H9b REFUTED (F28), reference drift (F29), H2d inconclusive + H2d′ (F30); Candela 2022 cited in P1, all P1 references verified | research | Done, research/CYCLE9_REPORT.md |
| C47-4 | Still open: MSW implicit segment coupling; 6_ENDSCALE krow front; KRNUM-03Z; GRCTRL STW (0A1 9.8 %); WGPR:OP02 offset of the ACTIONX model family; plot legend / keyword-doc translation; Metal back end; full 209-case OPM rerun; comp / thermal CLI snapshots (runDeckEngine); research: H10e two-way coupling, H2e gas sub-cell storage, H11 moment-aware design | mixed | Open |
| C48-R | Research Cycle 10: H10e REFUTED 4/4 (fault-permeability switching breaks convexity; no-switch controls 0/40), H2e SUPPORTED (SUBGRID LEAKC), H11 SUPPORTED (GP-EI 1.002×, seeds 1.001–1.006×); F31–F34, CYCLE10_REPORT.md; P1 abstract / results / limitations updated | research | Done |
| C48-W | Defaulted BHP reference depth = first defined connections incl. SHUT (OPM Well::updateRefDepth): ACTIONX_WEFAC and ACTIONX_GCONINJE pass | engine | Done |
| C48-U | Group UDA: GCONPROD / GCONINJE rate items given by FU_ / GU_ UDQ names (ACTIONX_GCONPROD 69.8 % → 17.1 %; rest is step-size dependence of a compounding UDA) | engine | Done |
| C48-S | comp / thermal cell snapshots shared by the GUI and `ResSimCLI -o` (`cell_snapshot.h`); CLI `--csv` per-well columns | engine / CLI | Done |
| C48-L | SUBGRID LEAKC gas-aware sub-cell storage (parse, export, engine) | engine | Done |
| C48-I | i18n: plot axis labels, legends and titles through an ImPlot translation hook (IDs unchanged), Case Explorer tree, section help; 335 new Spanish entries | app | Done |
| C48-T | `-c48` tests (19): LEAKC round trip, SHUT reference depth, group UDA, comp / thermal snapshots | tests | Done |
| C48-4 | Still open: rest of the 209-case OPM rerun (115 done: 48 / 52 / 15 skip, `docs/regression_C48_full_partial.txt`); MSW implicit segment coupling (not started); 6_ENDSCALE krow front; KRNUM-03* and GRCTRL STW (investigated, unresolved: logic matches OPM); ACTIONX_GCONPROD step-size dependence; window / popup title translation (needs ID-stable titles) and long keyword docs (English); Metal back end; Cycle 11: H10f reactivation constraint, H11b, H2f | mixed | Open |
| C49-H | Converged multisegment connection heads (per-connection Illinois-bracketed Picard after each converged step, ≤ 8 passes, 0.02 atm, "located" bracket criterion), default on; MSWOPTS 4th value, RESSIM_MSW_HEADS, UI combo. TEST_WSEGSICD 93.4 → 34.0 %, no OPM MSW case worse; MSW cross-flow FGOR 23.8 → 1.87 % (stays opt-in: neutral on OPM cases) | engine | Done |
| C49-F | FAULTMECH item 10 ramp width; `--fault-csv` columns sne / mult; UI column | engine / app | Done |
| C49-W | STOP-status wells report the zero-inflow BHP as WBHP (OPM; SPE1CASE1_SHUT_CONTROLS 99 → < 0.1 % on the stopped wells) | engine | Done |
| C49-C | Compositional / CO2STORE / thermal runs through api::runCase keep per-well rates, BHP, FWCT, FGOR (regression and studies had WBHP = 0): CO2STORE, CO2STORE_DIFFUSIVE now pass | engine | Done |
| C49-I | i18n: window titles, dock tabs and popup titles via a display-only ImGui hook (IDs unchanged); deck-editor toolbar and empty state; MSW / ramp strings | app | Done |
| C49-R | Research Cycle 11: H10f-a SUPPORTED but confounded → exploratory dose-response (violations only above a ≈ 2× transmissibility enhancement, F35); H10f-b non-discriminating; H11b SUPPORTED (F36); P1 JGR draft v0.6 | research | Done, research/CYCLE11_REPORT.md |
| C49-G | Full 209-case OPM rerun completed: 102 PASS / 77 FAIL / 30 SKIP (C48: 31 of 73 in the short set; `docs/regression_209_full.txt`) | tests | Done |
| C49-T | `-c49` tests (26): MSWOPTS heads round trip, converged heads + restart, FAULTMECH ramp, STOP-well WBHP, comp per-well summary; full suite 1159 PASS / 0 FAIL | tests | Done |
| C49-4 | Still open: BCCON / BCPROP boundary conditions (CO2STORE_DRSDTCON, CO2STORE_GW_DIRICHLET; warned "not modelled"); CO2STORE_GW gas rate 6.8 % (VAPWAT); MSW fully implicit segment equations (the Picard coupling leaves 40/58 VERT-02 steps at the pass limit); 6_ENDSCALE krow front; KRNUM-03*; GRCTRL STW; ACTIONX_GCONPROD Newton cut at 1097 d (documented); SOLVENT decks without PVDS (CO2SOL / H2SOL); long keyword docs in English; Metal back end; Cycle 12 (H10g, H10f-b re-designed, H2f) | mixed | Open |
| C47-F | `--fault-csv` / RESSIM_FAULT_CSV per-face fault history; CO2STORE false warning removed; tests `-c47` | engine / CLI | Done |
| C28-r ✅ | Research Cycle 7: H7 rate-and-state convexity (counts convex up to one event, moment not; F20), H2b sub-grid leak-off (SUPPORTED, F21), literature watch | research | Done, research/CYCLE7_REPORT.md |
| C28-2 ✅ | FAULTSLIP thermal + compositional, MPI distributed assembly, contiguous BILU0, phase field in the run | geomech/HPC | Done (rows C26-1, C21-3, C23-b, C26-4) |
| C24-r | Research Cycle 6: H8 Basel (needs digitised injection log), H7 rate-and-state design (P2) | research | Next |
| C23-a ✅ | Sub-grid leak-off (Carter-type fresh-rock correction of the tip host cell) for HFGROW | research/HF | **C28 done** (HFGROW … SUBGRID LEAK; H2b SUPPORTED: t60 spread 0.85 d over 10/5/2.5 m vs 6.6 d; F21) |
| C23-b ✅ | Level-ordered contiguous BILU0 (makes level scheduling pay) ; device AMG for the GPU path | HPC | **C28: contiguous BILU0 done** (topological renumbering by forward levels, bitwise = serial, default on for n ≥ 20k; 54k cells / 2 threads: 2.09 vs 2.18 ms serial vs 2.39 ms level-scheduled, 19 its each; -c28ilu). Device AMG still open. Earlier: **Deferred (C26, with reason):** on 2 cores the level-scheduled BILU0 was slower (C21), and block-Jacobi ILU per thread/rank costs +40 % iterations (C21 DD measurement); can only be judged on ≥ 4 physical cores. AMG reuse (C26-5) was done instead. Appleyard chopping and the well Schur complement (Phase D 24) already exist |
| C21-2 | SPE11 α_L≠α_T (tensor dispersion), then SPE11B on a coarse grid | verification | Do. 11C (3D) only on the user's machine: too large for 2 cores and 2 GB |
| C21-3 ✅ | DD distributed assembly (each rank assembles its own rows; the master stops gathering) | HPC | **C28 done** (`--mpi-da`, HPC Assembly combo; bitwise identical over ranks = serial; SPE9 2 ranks 5.82 s vs 7.32 s rank-0 assembly vs 7.47 s serial; -c28da, -mpida). Earlier: **Defer.** On 2 ranks the solve is already ~1.07× serial; it only pays off at ≥ 4 ranks on the user's PC. Planned once MPI is confirmed working there. Still deferred in C26 (same reason) |
| C22 | Compositional/thermal DD | HPC | **Done early in C20 (CLI, solve distributed)**: the hook was risk-free (same distributed solver, block size ≤ 16) |
| — | **SPE3 compositional** | verification | **Blocked by design.** Kenyon & Behie gave lab PVT only; each participant fitted their own EOS. Option: fit a PR EOS ourselves to the published CCE/CVD tables, labelled "ResSim fit". The black-oil SPE3CASE1 (decks/spe3) can run now |
| — | Solver profiling (H3) | engineering | Run between rounds; generic tricks are not novel |

## Industry-grade gap audit (living list, updated every round)

Legend: ✅ done · 🟡 partial · ❌ missing. "Own" is the ResSim engine. "OPM GUI" means the case can be built and edited in the app and run by OPM Flow.

| Area | Own engine | OPM GUI (build in the app) | Next action |
|---|---|---|---|
| Black-oil PVT (PVTO/PVDO/PVTG/PVDG/PVTW, DISGAS, VAPOIL) | ✅ | ✅ | PVCDO, PVZG, RSCONST |
| Multiple PVT/SAT/EQL/FIP regions | ✅ | 🟡 region tables are editable only for region 1 | region selector in the PVT/SCAL/Init editors |
| SCAL (SWOF/SGOF), 3-phase kro: default, STONE1 (+STONE1EX), STONE2 | ✅ (Stone within 1 % of OPM) | ✅ model combo in Rel. perm. | SWFN/SGFN/SOF3 family, SLGOF |
| End-point scaling (SWL, SWCR, …, ENDSCALE) | ✅ two-point (kr + Pc), ✅ three-point (SCALECRS) | ✅ | KRW/KRO/KRG/PCW value scaling, directional; all 8 Norne arrays now match the ECLIPSE INIT file cell by cell |
| Hysteresis (SATOPTS HYSTER, EHYSTR, IMBNUM, I* end points), VAPPARS | ✅ Carlson for gas/oil (non-wetting), Killough approximated by Carlson, wetting phases on drainage; VAPPARS (So/So,max)^v on Rs/Rv | ✅ round trip | true Killough, Pc hysteresis, PCW scaling |
| Grid: Cartesian, corner-point, faults, MULTFLT, NNC, MINPV, NEWTRAN, PINCH, MULTREGT (incl. NNCs, M/F/O sets) | ✅ | ✅ | EDITNNC, PINCHREG (A6), SCHEDULE-section MULTREGT, LGR |
| Dual porosity / permeability (DUALPORO, DUALPERM, SIGMA(V), LTOSIGMA, MULTMF, NODPPM) | ✅ | ✅ (ECLIPSE; OPM Flow has no DUALPORO) | gravity drainage (GRAVDR), DPKRMOD, verify against ECLIPSE |
| Equilibration: EQUIL, RSVD, RVVD, SWATINIT, THPRES, enumerated | ✅ | ✅ | PBVD/PDVD, THPRES IRREVERS |
| Aquifers: Fetkovich, Carter-Tracy | ✅ | ✅ | AQUTAB, numerical aquifers |
| Wells: WCONPROD/INJE/HIST/INJH, COMPDAT (timed, multiple open periods), WELOPEN (well + connection level), WPIMULT, RESV at field-average conditions | ✅ | ✅ | COMPLUMP numbers (kept raw), per-connection sat-table in the own engine |
| Groups (GRUPTREE, GCONPROD, GCONPRI, GCONINJE) | ✅ group tree, production targets/limits ORAT/WRAT/GRAT/LRAT/RESV (GCONPRI read as limits), WCONPROD GRUP wells, injection RATE/RESV/REIN/VREP; shares ∝ well potentials with water filling (each well capped by its own target and potential) | ✅ group tree editor + schema forms, run by OPM | guide-rate formulas (GUIDERAT), priority allocation, sales/fuel (GCONSALE/GCONSUMP), networks, GLIFT |
| VFP tables, THP control | ✅ VFPPROD/VFPINJ 5-D multilinear interpolation (FLO OIL/LIQ/GAS, WFR WOR/WCT/WGR, GFR GOR/GLR/OGR, ALQ); THP limit solved as the IPR–VFP intersection (bisection, OPM computeBhpAtThpLimit style) with datum-depth head correction | 🟡 VFP tables via schema form / Deck Editor; THP/VFP/ALQ in the schedule editor | VFP curve viewer, ALQ optimisation, unstable-branch detection |
| Multisegment wells | ❌ | 🟡 WELSEGS/COMPSEGS/WSEGVALV as text, run by OPM | B14 |
| Sensitivity / history matching | ✅ Studies workspace (tornado, MC, DoE + proxy, LM, ES-MDA); **C9: discrete adjoint gradient + L-BFGS, LM with exact forward-sensitivity Jacobian, gradient check**; **C12: analytic per-cell PERM / PORO / per-face T / WI adjoint gradients and per-cell field HM**; **C13: dedicated adjoint method + Adjoint HM view (iterations, final parameters / field maps, per-response check of the best case)**; HM mode with summary-file observed data, match-quality table/heat map, parameter posteriors | ❌ studies run on the own engine only | OPM Flow as the study engine (parallel batch runs), MCMC / adaptive proxy, multi-objective GA |
| WECON, WTEST, WELTARG, ACTIONX, UDQ | 🟡 WELTARG; ✅ WECON (rate / ratio limits, WELL / NONE, end run) + WTEST 'E' re-opening (C27); ACTIONX / UDQ not evaluated | 🟡 as text, run by OPM (UDQ expressions kept verbatim) | B15/B16 |
| Rs/Rv rate limits (DRSDT ALL/FREE; DRVDT) | ✅ DRSDT (SPE1 + DRSDT 0 within 0.04 % of OPM); ❌ DRVDT | ✅ DRSDT in the PVT section | DRVDT, DRSDTCON |
| Polymer (POLYMER, PLYVISC/PLYADS/PLYROCK/PLYMAX/PLMIXPAR/PLMIXNUM/PLYSHLOG/WPOLYMER/SPOLY) | ✅ Todd-Longstaff, adsorption, RRF, DPV, shear thinning, well injectivity (within 0.4–2 % of OPM) | ✅ EOR section, POLYMER/PADS maps, FC*/WC* vectors, example POLYMER | SHRATE form, PLYATEMP, brine coupling |
| Compositional EOS (PR/PR78/SRK, COMPS/CNAMES/TCRIT/…/BIC, ZI/ZMFVD, WELLSTRE/WINJGAS, FIELDSEP, CO2SOL) | ✅ FIM natural-variable core, SPE5 within a few % of E300, CO2 as an EOS component with CO2SOL brine solubility, E300 FIELDSEP semantics | ✅ Fluid → Compositional editor (components, BIC, composition, streams, separator, phase envelope, flash, CO2SOL, CO2 EOR tools: swelling / FCMP, mixing-cell MMP, WAG generator) | SPE3 comparison (needs the user's E300 SPE3 deck); ✅ SOLUBILI, slim-tube MMP, GUIDERAT, water vaporisation done in C7 |
| CO2STORE (CO2 + brine, SALINITY, RTEMPVD) | ✅ dedicated co2brine module (Spycher–Pruess / Duan–Sun); all trapping mechanisms: structural, residual (Carlson/Killough/Land hysteresis), solubility + diffusion, salt dry-out, mineral kinetics; six-class inventory | ✅ Fluid → CO2 storage editor incl. Trapping & interactions, CO2 results tab, SGTRAP/SSALT/MINCO2/PORFAC maps, examples `co2`, `co2trap` | Span–Wagner density, T > 100 °C, thermal, H2STORE, geochemical speciation, real SPE11A geometry / SPE11B-C (✅ C12 SPE11A-like slice) |
| Thermal (E300 THERMAL, STARS-like: steamflood / CSS / SAGD) | ✅ FIM p-N-U core, IAPWS-IF97, KVCR, heat loss, WINJTEMP; examples steamflood / css / sagd | ✅ thermal editor + Thermal results tab; ECLIPSE runs e300 | Pc, wellbore heat, EOS thermal flash (✅ steam trap C9, ✅ reactions C8, ✅ supercritical water C11, ✅ C12 near-critical fold construction, ✅ C12 Langmuir sorption / ECBM, ✅ C12 EDFM + Barton–Bandis / hydro-shear, ✅ C12 SPE4 1A, ✅ C13 EGS designer + fracture visualisation + Geothermal analysis, ✅ C13 Gringarten check); real-mixture EOS for water + gases, implicit fracture aperture, nanoparticles |
| Tracers, solvent, foam, brine, surfactant | 🟡 tracers ✅ (C27: TRACER / WTRACER / TBLKF, decay + adsorption via RESSIM TRACERPROP; no partitioning / dispersion); solvent, foam, brine, surfactant ❌ | 🟡 keywords kept as text and run by OPM; schema forms | next EOR modules |
| Unknown / unsupported keywords | ignored by the own engine (listed in the log) | ✅ kept, editable, written back in order; unedited imports run the original deck (pass-through) | — |
| Streamlines | ✅ own TOF tracer (Pollock) | ✅ ECLIPSE FrontSim engine option; viewer uses simulator FLR* fluxes | streamline-based saturation transport in the own engine |
| Schema-driven keyword forms (all OPM keywords) | n/a | ✅ forms for record keywords and tables (Deck Editor, Advanced page) | forms for array keywords (box editor) |
| Binary output (EGRID/INIT/UNRST/SMSPEC) | 🟡 reads them; writes its own .rsim | ✅ via OPM | writer for the own engine (C18) |
| Restart (RESTART, SAVE) | ❌ | ❌ | C19 |
| Time-step control (dS/dP targets), Appleyard chopping | ✅ Appleyard-style chopping, Newton-count growth / cut, and (C27) dS / dP target-change control (RESSIM TUNING items 8 / 9, Aziz–Settari NORM) | n/a | target-change step control, well Schur complement (D24) |
| Linear solver (well Schur + BSR, block CPR true/quasi-IMPES + AMG with reuse, BILU0 / float BILU0; GPU block-Jacobi) | ✅ | n/a | multicolour parallel BILU0, device ILU, domain decomposition |
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
    - ✅ AMG setup reuse across Newton iterations and time steps (C26, default 10 steps, quality monitor).
    - ILU0 with level scheduling, or a block-Jacobi + ILU0 smoother per thread.
    - Target: under 10 linear iterations per Newton.
24. **Nonlinear:**
    - ✅ Appleyard chopping and saturation-dependent damping (present).
    - Time-step control from the target saturation and pressure change (TUNING honoured safely).
    - ✅ Well-equation Schur complement (present).
25. **Threading heuristics:** ~~serial below 20k cells~~ dropped: SPE9 (9k cells) runs 1.4× faster on 2 threads. Add NUMA-aware first touch for large models.
26. **GPU (optional backend):**
    - Port the CPR stages (spmv, ILU0 apply, AMG V-cycle) to a D3D11/D3D12 compute or CUDA path behind the existing solver interface.
    - Expected gain: 3–6× on more than 1M cells.
27. **Distributed memory:** METIS/Zoltan partitioning, halo exchange and an MPI build. This is a long-term item and only matters past about 5M cells.

Performance targets (8-core desktop): SPE9 under 2 s, Norne full history under 3 min, Drogon full history under 10 min, and 1M cells × 10 years under 1 h.

## Phase E: Beyond E100 (later)

- Hysteresis (Killough, Carlson) and IMBNUM.
- ✅ Dual porosity / dual permeability (DUALPORO, DUALPERM, SIGMA, GRAVDR, DPNUM).
- ✅ Polymer (done). Solvent (4-component), brine, surfactant, foam: next EOR modules.
- Thermal (E300-like energy equation).
- CO2STORE / H2STORE for CCS.
- ✅ Local grid refinement (LGR, CARFIN, REFINE, WELSPECL/COMPDATL).
- Compositional (E300 class): 🟡 C1 thermodynamics library, the C2 FIM core (water, tables, Pc, implicit wells) and the C3 deck hook-up are done, so E300 decks run on the own engine (SPE5 WAG). Next: C3 remainder (SEPCOND, groups, SPE3), then C4 CO2STORE, then T1 thermal.

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
