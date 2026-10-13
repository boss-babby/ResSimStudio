# ResSim Studio — Phases 1–3

ResSim Studio is a desktop black-oil reservoir simulator with a modern dockable UI. The engine is portable C++17. The front-end is Win32 with Direct3D 11 for the 3D grid viewer, Direct2D/DirectWrite for viewport text and legends, and Dear ImGui (docking) plus ImPlot for the panels and charts.

You set up every input through forms, tables, generators and live plots. No deck editing is needed.

---

## 1. Build (Visual Studio 2022)

Requirements:
- Windows 10/11.
- VS 2022 with the "Desktop development with C++" workload. It includes CMake and the Windows SDK.

**Option A: open the folder.** Choose *File → Open → Folder…* and select `ressim/`. Pick the **x64-Release** configuration, then set **ResSimStudio.exe** as the startup item and press F5.

**Option B: command line.** Use the "x64 Native Tools Command Prompt for VS 2022":
```bat
cd ressim
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
build\Release\ressim_tests.exe      :: engine regression tests
build\Release\ResSimStudio.exe
```

Options:
- `-DRESSIM_AVX2=OFF` for older CPUs.
- `-DRESSIM_BUILD_TESTS=OFF` to skip the tests.
- OpenMP (`/openmp`) is enabled automatically. Set the thread count in *Numerics*, where 0 means all cores.

**Option C: MinGW-w64 cross-build from Linux.** This produces a single static exe:
```bash
sudo apt install g++-mingw-w64-x86-64-posix
./build_mingw.sh        # -> build/mingw/ResSimStudio.exe and ressim_tests.exe
```
A ready-built `release/ResSimStudio.exe` is included. It is fully static and needs no extra DLLs on Windows 10/11.

The engine and tests also build on Linux:
```bash
g++ -std=c++17 -O2 -fopenmp -Iengine/include engine/src/*.cpp tests/run_tests.cpp -o ressim_tests
```

## 2. Using the application

| Area | What it does |
|---|---|
| **Case Explorer** (left) | Lists the input sections with a live status dot (green = OK, red = problem; hover for the message). It also shows results status and a model summary. |
| **Inspector** (right) | Editor for the selected section. See the table below. |
| **3D Viewer** (centre) | Shows any static, initial or dynamic property.<br>• Controls: orbit (LMB), pan (RMB/MMB), zoom (wheel), view presets, Z-exaggeration.<br>• Filtering: I/J/K slicing and a value filter.<br>• Display: colour maps, log scale, edges, wells with perforations, legend and axis triad.<br>• Hover a cell to read its value; animate through the report steps. |
| **Plots** (bottom) | Tabs for Field rates and cumulatives, per-well rates/BHP/water cut/GOR with history overlay, a custom vector plotter, and solver convergence. The time cursor is synced with the 3D viewer; double-click a plot to jump there. |
| **Run Monitor** | Shows progress, KPIs (simulated time, average pressure, oil recovery factor, wall-clock time), current rates, cumulative and in-place volumes, and solver statistics. |
| **Log** | Engine log with colour coding, a filter and copy. |

What each Inspector section contains:

| Section | Contents |
|---|---|
| Run specification | Title, units (FIELD/METRIC), active phases, dissolved gas. |
| Grid & rock | Geometry, dip, a layer table (DZ, k, φ, kv/kh, NTG) and a geostatistical heterogeneity generator (σ ln k, correlation lengths, seed). Also shows statistics and a permeability histogram. |
| Fluid PVT | PVTO/PVDO, PVDG, PVTW, ROCK and DENSITY tables with plots. A correlation generator covers Standing, Beggs–Robinson, Vasquez–Beggs, Sutton, Hall–Yarborough and Lee–Gonzalez–Eakin. |
| Rel. perm & Pc | Corey generator, SWOF/SGOF tables, kr curves, and the Buckley–Leverett fractional-flow curve with the Welge tangent and shock-front saturation. |
| Initialization | EQUIL contacts and datum. *Compute initial state* shows STOIIP and GIIP plus pressure- and saturation-vs-depth plots with the contacts drawn in. |
| Wells | Producers and water/gas injectors. Place a well by clicking on the top-view permeability map. Completions use Peaceman well indices, or you can enter explicit CF values. |
| Schedule | Report steps (uniform builder or per-step) and well controls: ORAT, WRAT, GRAT, LRAT, RESV, BHP and RATE. A WCONHIST history flag is included. |
| Numerics | Time step control, Newton/CNV/material-balance tolerances, linear solver settings and CPU threads. |

Shortcuts:

| Key | Action |
|---|---|
| F5 | Run |
| Shift+F5 | Stop |
| Ctrl+S | Save |
| Ctrl+Shift+S | Save as |
| Ctrl+O | Open |
| Ctrl+N | New blank model |

You can also drag and drop a `.rsim` file onto the window. *File → Export ECLIPSE deck* writes an equivalent `.DATA` file.

Built-in cases:
- **Demo:** a 3-phase waterflood on a 20×20×4 grid using SPE1 PVT. It runs in about 4 s.
- **1D Buckley–Leverett:** validation against the analytic solution (front within 4 %).
- **Blank model:** generated from correlations and Corey curves.

## 3. Engine (what is simulated)

**Formulation:**
- Fully-implicit black-oil model: water, oil and gas components with dissolved gas (Rs).
- Variable switching between Sg and Rs when gas appears or disappears.
- Gravity, capillary pressure (Pcow, Pcog), rock and water compressibility.

**Saturation functions:** ECLIPSE default 3-phase kro (Stone-type segregated model).

**Wells:**
- Fully coupled, with BHP as an extra unknown.
- Rate/BHP control with automatic switching and hydrostatic wellbore head.
- Peaceman well index.

**Initialization:** hydrostatic equilibrium by phase-pressure integration, with NACC sub-cell averaging.

**Numerics:**
- Newton–Raphson with forward-mode automatic differentiation for exact Jacobians.
- Appleyard-type chopping and adaptive time stepping with cuts.
- Convergence on CNV plus material balance per component.

**Linear solver:**
- Block-sparse CSR matrix.
- ILU(0)-preconditioned BiCGSTAB, with GMRES fallback.
- OpenMP-parallel assembly, SpMV and vector kernels.

**I/O:** native `.rsim` project files (text, run-length compressed arrays) and ECLIPSE deck export.

Source layout:
```
engine/include/ressim/  case.h grid.h fluid.h ad.h linsolve.h simulator.h tools.h io.h
engine/src/             implementation
app/src/                main.cpp (Win32/D3D11/D2D), app.cpp (frame, menus, actions),
                        ui_inputs.cpp (editors), ui_results.cpp (viewer/plots/monitor),
                        viewport3d.* (D3D11 renderer + D2D overlay), sim_runner.* (worker thread),
                        theme.*, widgets.h
tests/run_tests.cpp     IO round-trip, Buckley-Leverett, correlations, demo material balance
third_party/            Dear ImGui (docking), ImPlot
```

**Threading:**
- The simulation runs on a worker thread, and results stream to the UI under a mutex.
- The UI stays responsive, and you can stop a run at any time.
- When idle, the app throttles its frame rate to save CPU and GPU.

## 4. Phase 2 and 3 features

### Grids, faults and aquifers
- **Corner-point geometry (COORD/ZCORN):** imported from decks or converted from a Cartesian grid ("Convert to corner-point").
  - Transmissibilities come from face-overlap geometry, so fault throws create non-neighbour connections (NNCs) automatically.
- **Faults:** the *Grid & rock* section has a fault-throw tool (orientation, position, throw and MULTFLT). Fault multipliers are editable in a table.
- **Explicit NNCs:** edited in a table, with 1-based indices.
- **Analytic aquifers:** set up in their own *Aquifers* section.
  - **Fetkovich:** pseudo-steady state (AQUFETP).
  - **Carter–Tracy:** infinite acting (AQUCT).
  - Both connect to cell faces (AQUANCON). An influx plot is shown after a run.

### ECLIPSE deck, two-way
- **Import:** File → Import, or Open a `.DATA` file. Warnings are listed in the Log.
- **Export:** File → Export ECLIPSE deck.
- **Round trip:** export → import → export gives an identical deck, verified by the test suite.
- **Keywords supported:**

  | Section | Keywords |
  |---|---|
  | RUNSPEC | DIMENS, phases, units, START, TITLE |
  | GRID | DX/DY/DZ/TOPS, DXV/DYV/DZV, SPECGRID/COORD/ZCORN, PORO/PERM*/NTG/ACTNUM/MULT* |
  | GRID editing | BOX/ENDBOX, EQUALS/MULTIPLY/ADD/COPY/MINVALUE/MAXVALUE |
  | GRID faults | FAULTS/MULTFLT/NNC, AQUANCON |
  | PROPS | PVTO/PVDO/PVCDO/PVDG/PVTW/ROCK/DENSITY/SWOF/SGOF |
  | SOLUTION | EQUIL/RSVD/AQUFETP/AQUCT |
  | SCHEDULE | WELSPECS/COMPDAT/WCONPROD/WCONINJE/WCONHIST/WCONINJH/WELOPEN/WPIMULT/TSTEP/DATES/TUNING |
  | Includes | INCLUDE |

- **ResSim extensions:** things ECLIPSE has no keyword for are carried in `-- RESSIM:` comment lines, so other simulators ignore them. These are the solver choice, observed data (`OBS`) and study parameters (`PARAM`).

### Linear solver: CPR-AMG
- **Method:** two-stage CPR. A quasi-IMPES pressure system is solved with smoothed-aggregation AMG (V-cycle, Gauss–Seidel smoothing), followed by ILU(0) on the full system.
- **Setting:** *Numerics → Preconditioner*. Automatic uses CPR for grids of 500 cells or more.
- **Benchmark** (28 800 cells, heterogeneous; identical results):

  | Solver | Linear iterations / Newton | Run time |
  |---|---|---|
  | ILU(0) | 76 | 71.6 s |
  | CPR-AMG | 7.4 | 23.7 s |

### Studies window (sensitivity and history matching)
- **Parameters:**
  - Box-limited perm X/Y, perm Z and porosity multipliers.
  - MULTFLT and aquifer strength.
  - krw end point and curvatures, water viscosity, rock compressibility.
  - OWC, datum pressure, well skin and PI.
  - Each parameter has low/base/high values and an optional log scale.
- **Observed data:**
  - From WCONHIST history, CSV import (QTY, WELL, TIME, VALUE[, SIGMA]), or a synthetic "twin experiment" from the last run with noise.
  - Misfit = mean of ((sim − obs)/σ)².
- **Sensitivity:**
  - Tornado (one at a time, 2N+1 runs).
  - Latin-hypercube Monte Carlo with P90/P50/P10, a histogram, a FOPT fan chart and Spearman rank correlations.
  - Responses: FOPT, RF, FWPT, FGPT, FWIT, FPR, FWCT, NPV (with an economics model) and misfit.
- **History matching:**
  - **Manual:** sliders with a live misfit (debounced background run) and observed-vs-simulated plots.
  - **Levenberg–Marquardt:** finite-difference Jacobian computed with parallel runs, plus parallel λ trials.
  - **ES-MDA:** ensemble smoother with multiple data assimilation, giving a posterior mean, P10–P90 and the ensemble envelope in the fit plot.
  - "Apply best match" writes the result into the case.
- **Performance:** every study uses a parallel batch runner, with one simulation per hardware thread.

### Streamlines
- **Method:** Pollock semi-analytical tracing in each cell's reference unit cube, using total reservoir face fluxes. It works unchanged across faults and NNCs.
- **Display:** in the 3D viewer (Streamlines checkbox), coloured by time of flight, injector or producer. "Hide cells" and "X-ray" options keep the lines visible.
- **Cell properties:** TOF (time of flight), DRAINAGE (producer region) and SWEEP (injector region).
- **Well allocation:** a table of injector → producer allocation factors.
- **Any report step:** the simulator restores that step's state before tracing.

## 5. Roadmap
- GPU linear solver (D3D11 compute) and GPU-side property evaluation.
- Multi-colour parallel ILU, LGR, group controls and economic limits.
