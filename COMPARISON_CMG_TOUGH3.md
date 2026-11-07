# ResSim Studio vs CMG and TOUGH3: feature comparison, C27 ideas, open gaps

*Round C27 (October 2026). Sources: the CMG product pages (IMEX, GEM, STARS, CMOST, Builder, CoFlow, WinProp, 2026), the
TOUGH3 page and User's Guide (LBNL-2001093; Jung et al. 2017/2018), the iTOUGH2 command reference and theory
(LBNL-1005838, Finsterle). This is our own reading of the public documentation. It is not a vendor benchmark, and the
"ResSim" column only marks features that have a test in `ressim_tests`.*

Legend: ✅ available · 🟡 partial / simpler model · ❌ not available · **C27** = added in this round.

## 1. What the two products are

| | CMG | TOUGH3 (LBNL) |
|---|---|---|
| Scope | Commercial oil and gas, CCS and geothermal suite: IMEX (black oil), GEM (EOS compositional + geochemistry), STARS (thermal / chemical / geomechanics), CMOST (AI-assisted HM / optimisation / UQ), Builder (pre-processing), Results, WinProp (PVT), CoFlow (integrated production / networks) | Research code for non-isothermal multiphase, multicomponent flow in porous and fractured media. One core plus exchangeable EOS modules. Companion tools: iTOUGH2 (inverse modelling), TOUGHREACT (reactive transport), TOUGH-FLAC (geomechanics), PyTOUGH |
| Numerics | Adaptive-implicit / fully implicit FV, massively parallel, IMEX-GPU | Integral finite differences (IFD), fully implicit Newton, PETSc solvers, MPI |
| Strengths | Production engineering, field development, EOR and thermal recovery workflows, polished pre/post-processing, AI-assisted HM | Geothermal, vadose zone, CO2 storage, nuclear waste: wide EOS family (water / brine / air / CO2 / H2 / VOC / radionuclides), MINC, rigorous parameter estimation and uncertainty statistics (iTOUGH2) |

## 2. Feature matrix

### 2.1 Flow physics

| Feature | CMG | TOUGH3 | ResSim Studio |
|---|---|---|---|
| Black oil (live oil, wet gas, ECLIPSE keywords) | ✅ IMEX | 🟡 EOS8 (water, dead oil, NCG) | ✅ own FIM engine + OPM Flow / ECLIPSE runners |
| EOS compositional (PR / SRK, CO2, FIELDSEP) | ✅ GEM | 🟡 ECO2N / ECO2M / EOS7C / TMVOC (fixed component sets) | ✅ natural-variable FIM, any component list |
| CO2 storage (brine solubility, trapping inventory, salt) | ✅ GEM | ✅ ECO2N (V2), ECO2M, EWASG (salt precipitation, porosity / permeability change) | ✅ CO2STORE module, all trapping classes, salt dry-out, mineral kinetics |
| Thermal (steam, SAGD, CSS, combustion) | ✅ STARS | ✅ EOS1 / EOS3 / EOS4 (non-isothermal water, air, vapour-pressure lowering) | ✅ p-N-U thermal core, IAPWS-IF97, reactions, STARS-like inputs |
| Geothermal / EGS | ✅ STARS (EGS, closed loop) | ✅ core strength | ✅ EGS designer, EDFM, Barton–Bandis, THM |
| Richards / unsaturated zone | ❌ | ✅ EOS9 | 🟡 van Genuchten–Mualem tables (**C27** generator), no dedicated single-phase Richards equation |
| Passive tracers | ✅ IMEX "advanced tracer module" | ✅ EOS1 (water + tracer), EOS7C / EOS7R | ✅ **C27** TRACER / WTRACER / TBLKF + decay and linear adsorption (TRACERPROP) in the black-oil engine |
| Radionuclide decay chains (parent → daughter) | ❌ | ✅ EOS7R | 🟡 **C27** first-order decay of independent tracers; no parent–daughter chain |
| Dual porosity / MINC | ✅ | ✅ MINC | ✅ DUALPORO / DUALPERM / MINC in every engine |
| Polymer, surfactant, ASP, foam | ✅ STARS / GEM / IMEX polymer | ❌ | 🟡 polymer native, other chemical EOR via OPM Flow |
| Geochemistry (aqueous speciation, minerals) | ✅ GEM | ✅ TOUGHREACT | 🟡 CO2 mineral kinetics only |
| Geomechanics | ✅ 1/2-way FE incl. thermal | ✅ TOUGH-FLAC | ✅ FE (one-way, two-way, **C26** iterative), DDM fault slip, LEFM growth, phase field (2-D tool) |
| Induced seismicity (Coulomb, Dieterich, Mw catalogue) | ❌ | 🟡 through TOUGH-FLAC research workflows | ✅ in the main run (C24–C26), seismicity designer, **C27** seismicity responses in studies |

### 2.2 Wells and field management

| Feature | CMG | TOUGH3 | ResSim Studio |
|---|---|---|---|
| Rate / BHP control, groups, VFP / THP | ✅ | 🟡 deliverability / rate sources | ✅ |
| Economic limits (WECON-type) and well re-tests | ✅ | ❌ | ✅ **C27** WECON (min oil / gas / liquid, max WCT / GOR / WGR, WELL / NONE, end run) + WTEST 'E' re-opening |
| Discretised transient wellbore, ICD / AICD | ✅ FlexWell | 🟡 T2Well (separate code) | ❌ (standard well model; MSW keywords run in OPM Flow) |
| Surface networks, nodal analysis, decline curves | ✅ CoFlow | ❌ | ❌ |

### 2.3 Numerics and HPC

| Feature | CMG | TOUGH3 | ResSim Studio |
|---|---|---|---|
| Fully implicit Newton, CPR-AMG | ✅ | ✅ (PETSc) | ✅ own CPR-AMG, BSR, mixed precision |
| Time-step control on target changes (ds / dp per step) | ✅ (DTMAX, NORM) | ✅ (automatic time-step adjustment) | ✅ **C27** TUNING dsTarget / dpTarget (Aziz–Settari NORM rule) |
| Dynamic grid refinement / coarsening | ✅ DynaGrid (GEM) | ❌ | ❌ (static LGR only) |
| MPI domain decomposition | ✅ | ✅ | ✅ MS-MPI / Intel MPI; the GPU path is experimental |

### 2.4 History matching, sensitivity, uncertainty

| Feature | CMG CMOST | iTOUGH2 | ResSim Studio |
|---|---|---|---|
| Tornado / DoE / Monte Carlo / proxy | ✅ (incl. ML proxies) | ✅ | ✅ |
| Gradient HM (LM, adjoint) | 🟡 | ✅ LM | ✅ LM, ES-MDA, discrete adjoint |
| Parameter covariance, correlation, 95 % intervals | 🟡 | ✅ | ✅ **C27** HM statistics panel |
| Residual analysis (normality, studentised residuals, leverage / influence) | ❌ | ✅ | ✅ **C27** |
| Contribution of each dataset to the objective | 🟡 | ✅ | ✅ **C27** |
| Identifiability: eigen-analysis, composite scaled sensitivity | ❌ | ✅ | ✅ **C27** (CSS per parameter, eigenvalues, condition number) |
| Data worth / FOSM prediction uncertainty (J C_pp Jᵀ) | ❌ | ✅ | ❌ (next) |
| Coupled-physics KPIs inside studies (seismicity, uplift) | ❌ | 🟡 | ✅ **C27** (fault mechanics / geomechanics run in every study case) |

### 2.5 Pre-processing and usability

| Feature | CMG Builder / Results | TOUGH3 (+ PyTOUGH, PetraSim) | ResSim Studio |
|---|---|---|---|
| Rel-perm / Pc correlations | ✅ Corey / LET / Stone | ✅ built-in functions (van Genuchten, Brooks–Corey, …) | ✅ Corey + **C27** van Genuchten–Mualem SWOF / SGOF generator |
| Example library from the command line | 🟡 | ✅ sample problems | ✅ **C27** `ResSimCLI --list-examples / --example` (31 built-in) |
| Geostatistics (SGS, kriging) | ✅ | ❌ | ❌ |
| Scripting | ✅ Python (CMG) | ✅ PyTOUGH | 🟡 CLI with `--set` / `--json` / `--csv`; no embedded Python |

## 3. Implemented in C27 (ideas taken from this comparison)

1. **Passive tracers** (IMEX tracer module, TOUGH3 EOS1 / EOS7R): ECLIPSE TRACER / WTRACER / TBLKF are simulated by the
   ResSim black-oil engine, with first-order decay and linear equilibrium adsorption added through
   `-- RESSIM: TRACERPROP`. The transport is upwind and fully implicit on the converged fluxes, so it is an M-matrix with
   a maximum principle. Results: Tracers tab (breakthrough curve, produced concentration, recovery of the injected
   tracer), custom-plot vectors, 3-D `TR_<name>`, restart and deck round trip.
2. **Economic limits and well tests** (CMG / ECLIPSE WECON + WTEST).
3. **Target-change time-step control** (TUNING dsTarget / dpTarget).
4. **iTOUGH2-style inversion statistics**: covariance with the s0² scaling, correlation, 95 % intervals, CSS,
   eigen-analysis, residual statistics, studentised residuals, leverage and dataset contributions.
5. **Coupled-physics study runs**: seismicity and uplift responses in tornado, Monte Carlo and HM; 14 study responses.
6. **Van Genuchten–Mualem** SWOF / SGOF generator (TOUGH relative-permeability / capillary functions as ECLIPSE tables).
7. **Example library on the CLI.**

## 4. Gaps after C28 (status, ordered by value)

| Gap | Reference | Status |
|---|---|---|
| Data-worth / FOSM prediction uncertainty | iTOUGH2 | **Done in C28** (Studies → HM → Inversion statistics → Prediction uncertainty and data worth) |
| Parent → daughter decay chains, partitioning tracers (SOR tests), tracer dispersion | TOUGH3 EOS7R, IMEX | **Done in C28** (TRACERPROP K / parent / alphaL / Dm; PITT S_or 0.2013 for 0.2) |
| Geostatistics (SGS / kriging) | CMG Builder | **Done in C28** (Grid & rock → Kriging / SGS conditioned to wells) |
| ML proxy / AI-assisted HM | CMG CMOST | **Done in C28** as a Gaussian-process proxy (ARD, ML hyper-parameters, prediction band); no neural proxy |
| Nodal analysis, decline curves | CMG CoFlow | **Done in C28** (Results → Production analysis: Arps, IPR / VLP, VFPPROD); surface networks and tank models still open |
| Phase field inside the field run, C23-b level-ordered BILU0, C21-3 distributed assembly | own roadmap | **Done in C28** (HFGROW … PHASEFIELD; contiguous BILU0; `--mpi-da`) |
| Salt precipitation coupled to black-oil brine (EWASG-type) | TOUGH3 EWASG | **Done in C29** (`-- RESSIM: SALT`: dry-out by evaporation into gas, halite precipitation / redissolution, Verma–Pruess k(φ); SSALT / KSALT; example SALT_DRYOUT) |
| Discretised wellbore (FlexWell / T2Well), ICD / AICD | CMG STARS, TOUGH | **Partly done in C29** (FlexWell-lite `-- RESSIM: WELLBORE`: heel-to-toe friction + per-connection ICD, implicit per connection; example WELLBORE_ICD). Open: transient wellbore, slip / hold-up, annulus, AICD; MSW decks run in OPM Flow |
| Dynamic grid refinement (DynaGrid) | CMG GEM | **Deferred (C29 design note in ROADMAP)**: needs re-gridding of every state array, wells and restart; static LGR exists |
| Surface network solver | CMG CoFlow / IPSM | **Done in C29, implicit since C30** (NODEPROP / BRANPROP with branch VFP tables → producer THP limits, node pressures iterated to agree with the well rates; example NETWORK). Open: chokes / compressors / separators |

C30: the broader industrial assessment (ECLIPSE, CMG, tNavigator, TOUGH, COMSOL) is in `INDUSTRY_GAP.md`.
