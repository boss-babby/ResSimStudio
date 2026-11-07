# ResSim Studio — UI / workflow review and modernisation plan (C39)

Honest, code-based audit of the desktop application (app/src, ~28 k lines of UI code over 40 files), written before
changing the UI. It answers three questions: (1) is the UI industry standard / product grade, (2) what are the gaps and
how are they closed, (3) how does the Geomodel workspace get to Petrel-class. Status markers: **[done C39]** = implemented
in this round, **[next]** = planned order of work.

---

## 0. Cross-platform foundation (done first, C39)

| Layer | Windows (unchanged) | Linux / macOS (new) |
|---|---|---|
| Window / input | Win32 | GLFW 3.3+ |
| 3-D views (Model, Results, 3D scene) | Direct3D 11 (HLSL) | OpenGL 3.3 core (GLSL), same geometry code |
| View overlays (legend, axes, labels) | Direct2D + DirectWrite | Dear ImGui draw lists into the same FBO — one shared layout (`paintOverlay`) |
| Figure PNG export | Direct2D + WIC | plutovg (+ DPI written as pHYs) |
| File dialogs, message boxes | common dialogs | tinyfiledialogs (zenity / kdialog / osascript) |
| Embedded data (fonts, icons, keyword manual) | RCDATA resources | `.incbin` (still a single self-contained executable) |
| GPU linear solver | sokol_gfx D3D11 compute | CPU solvers (reported as unavailable) |
| Build | build_mingw.sh / CMake+MSVC | CMake (`build_linux.sh`), GitHub Actions CI for all three |

Design rule: the `*_gl.cpp` files only contain GPU calls; camera, picking, geometry generation, colour mapping and the
overlay layout are compiled once for both back ends, so features cannot drift apart between platforms.

---

## 1. Verdict: is the UI industry standard?

**Short answer: functionally yes, as a product not yet.** The breadth of *what* can be edited, run and viewed is at or
beyond commercial pre-/post-processors (ECLIPSE/OPM two-way decks, compositional/thermal/geomechanics editors, HM and
sensitivity, a geomodelling workspace, 3-D seismic). *How* it is presented is still an expert tool built panel by
panel, round by round. Measured against Petrel, CMG Builder/Results, tNavigator and ResInsight:

| Criterion | Rating | Evidence in the code |
|---|---|---|
| Feature coverage / depth | ★★★★★ | 11 input sections, 5 workspaces, ~60 specialised editors, every keyword has UI |
| Correctness feedback | ★★★★☆ | live validation dot per section, problem box in the Inspector, toolbar status, deck checker |
| Discoverability | ★★☆☆☆ | features grew as collapsing headers inside sections (e.g. 7 in Inputs, EGS/HF/thermal panels), several Tools-menu windows; command palette exists but is the only search |
| Consistency | ★★★☆☆ | widgets.h gives shared props/headers, but buttons, tables and plot toolbars differ between older and newer panels |
| Object model | ★★☆☆☆ | the Case Explorer lists *sections*, not *objects* (wells, faults, regions, surfaces, cases, results); Petrel/tNavigator users expect an object tree with right-click actions |
| Undo / redo | ★☆☆☆☆ | only the deck editor (text) and horizon editing; no undo for model edits |
| Safety | ★★☆☆☆ | dirty flag + save prompt; no project autosave / crash recovery (only results autosave) |
| Visual design | ★★★☆☆ | coherent dark theme, Inter font, Lucide icons, toasts; but dense, few hierarchy levels, no light theme, no accent choice |
| Multi-case work | ★★★☆☆ | Result Viewer compares runs, Studies has ensembles; no case manager with side-by-side model diffs |
| Performance | ★★★★★ | GPU views, idle throttling, background runs, out-of-core seismic |
| Accessibility / i18n | ★☆☆☆☆ | fixed English strings, keyboard nav only via ImGui defaults, colour-only status dots |
| Platform | ★★★★☆ | Windows + (new) Linux/macOS; HiDPI via UI-size setting |

Conclusion: to be sold as a product the gaps are **workflow structure, object model, undo/autosave, visual hierarchy and
onboarding** — not missing features. None of these require changing the engine.

---

## 2. Gap list and how each is removed

Ordered by user impact per unit of work. (P1 = blocks "product grade", P2 = expected by professionals, P3 = polish.)

| # | Gap | Pri | Fix |
|---|---|---|---|
| G1 | No undo/redo for model edits | P1 | Whole-case snapshots (`shared_ptr<const Case>`, no serialisation, restored through the open-project path) taken at the end of a frame when no widget is active, so one drag / typed value = one step; Ctrl+Z / Ctrl+Y / Ctrl+Shift+Z, Edit menu with step count, toolbar arrows; depth 4-100 by a 640 MB budget; covers every Model-workspace editor at once (the text editors and the Geomodel keep their own undo) **[done C39]** |
| G2 | No project autosave / crash recovery | P1 | Background autosave (default 2 min, 15 s-1 h or off) to the per-user settings folder via temp file + rename; removed on save / clean exit; on start a *Recover unsaved work* dialog offers Restore / Discard **[done C39]** |
| G3 | Features hidden in collapsing headers and Tools windows | P1 | Feature registry: every editor gets a name, keywords, section and an "enabled in this case" predicate. The command palette, a searchable "Add physics…" menu and the Case Explorer all read it, so nothing is unreachable **[done C40]** |
| G4 | Case Explorer shows sections, not objects | P1 | Object tree: Grid ▸ properties, faults, regions · Fluids ▸ PVT tables / EOS · Wells ▸ each well (status icon, right-click: edit, duplicate, delete, show in 3-D, plot) · Schedule ▸ events · Results ▸ runs · Studies. Selection drives the Inspector **[done C40]** |
| G5 | Visual hierarchy too flat | P2 | Design tokens (spacing 4/8/12/16, radius 6, 3 text levels, elevation via surface tints), card-style groups in the Inspector, 2-column property grids everywhere (ui::beginProps), consistent primary/secondary/danger buttons **[done C39 + C41: themed surfaces; C41 tokens `themeSurface` / `tok::sp` / radius 6, card headers, primary / secondary / danger buttons]** |
| G6 | Dark theme only | P2 | Dark / Light / High-contrast themes and 5 accent colours (Settings > Interface), applied to ImGui and ImPlot; 3-D views keep their own background setting **[done C39; C40: View › Appearance submenu (theme, accent, 3-D background Dark / Follow theme / Light grey / White / Black)]** |
| G7 | No onboarding for general users | P2 | Start page already lists examples; add guided "new case" wizard (template → grid → fluids → wells → run) with progress breadcrumbs, and contextual empty-states ("No wells yet — Add producer / Import WELSPECS") **[done C40]** |
| G8 | Workspace switching | P2 | Already present: segmented workspace switcher in the toolbar (Model · Deck · Results · Studies · Geomodel), Ctrl+1…5. (No left activity bar — user decision.) C39 adds Undo/Redo next to Save and an Edit menu **[done]** |
| G9 | Plots: good export, weak linking | P2 | Linked cursors (time cursor shared by 3-D, plots and well table), brushing a time range in any plot sets the 3-D report step **[done C40]** |
| G10 | Case management | P2 | Case manager: list of cases/variants with parameter diffs, run status, compare-to-base in Results **[done C40: View › Case manager – named variants saved in `<project>.cases/`, base flag, input diff vs base / current model, runs linked automatically, load / load & run / update / duplicate / export, *Compare results with the base* in the Result Viewer]** |
| G11 | Long property lists without filtering | P3 | Filter box on Inspector (filters property rows by label), "show only modified" toggle **[done C40: filter box on top of the Inspector (Ctrl+F, Esc clears) filtering all standard property rows of the open groups; "show only modified" not done – no per-field defaults exist yet]** |
| G12 | Status by colour only | P3 | Icon + colour (✓ / ⚠ / ✕) for validation, high-contrast option **[done C40: Case Explorer sections and the workflow checklist use ✓ / ⚠ icons plus colour; high-contrast theme since C39]** |
| G13 | Strings not localisable | P3 | `_("...")` wrapper + catalogue loader; Bengali / Spanish / Chinese / Arabic later (fonts already merge CJK on Linux) **[partly done C41: `TR` / `TRN`, `lang/` catalogues, View › Language, Spanish for menus / toolbar / workspace switch; scripts needing shaping (Bengali, Arabic) not possible with the ImGui font renderer]** |

---

## 3. From "traditional" to modern — target design

Principles (what Petrel 2024, tNavigator 24, ResInsight 2024, CMG Launcher 2025 and modern IDEs converge on):

1. **One object model, many views.** Selecting a well in the tree, the 3-D view, a plot legend or a table selects it
   everywhere; the Inspector always edits the current selection.
2. **Task-oriented workspaces** (already present) with a visible switcher and per-workspace layouts.
3. **Progressive disclosure.** Basic properties first; "Advanced" expanders; physics modules added explicitly
   ("Add physics: thermal / compositional / geomechanics / chemical …") instead of always-visible headers.
4. **Command-first.** Ctrl+K palette searches commands, objects, keywords, settings and docs (exists; becomes the
   feature registry front end).
5. **Never lose work.** Undo everywhere, autosave, recovery, non-blocking long operations (already non-blocking).
6. **Calm visuals.** Neutral surfaces, one accent colour, 8-pt spacing grid, 3 type sizes, icons with labels in
   toolbars, subtle elevation; light and dark themes.
7. **Live feedback.** Validation as you type (exists), mini-previews (relperm curve next to the table — exists in parts),
   run progress with ETA and convergence sparkline in the status bar.

Layout target (Model workspace):

```
┌ menu ─────────────────────────────────────────────────────────────────────────────────┐
│ [New][Open][Save] [Undo][Redo]   ( Model | Deck | Results | Studies | Geomodel )   [▶ Run][■] engine ▾ │
├ Explorer (objects) ┬──────────── 3-D / map / section (tabs) ─────────────┬ Inspector ────┤
│ ▸ Grid            │                                                     │ selection     │
│ ▸ Fluids          │                                                     │ cards, filter │
│ ▾ Wells (6)       │                                                     │ advanced ▸    │
│   ● P1  ● I1 …    ├──────────── plots · run monitor · log · studies ────┤               │
│ ▸ Schedule        │                                                     │               │
├───────────────────┴─────────────────────────────────────────────────────┴───────────────┤
│ status: ready · file · units · engine · cells · notifications · threads                 │
└─────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Geomodel workspace → Petrel parity

What exists (C32–C36): data import (IRAP/ZMAP/XYZ, LAS, deviation, tops, fault sticks, SEG-Y incl. out-of-core),
velocity model + depth conversion, seeded auto-tracking, structural editing (brush, smooth, fill, isochore),
fault-parallel pillar gridding with zones/layering, well blocking, facies (SIS, TGS, objects, MPS), petrophysics (SGS /
kriging / co-simulation, poro-perm transform), volumetrics with P90/P50/P10, upscaling, GRDECL export / send to case,
3-D scene (slices, volume rendering, horizons, faults, wells, grid), seismic attributes.

Gap to Petrel, by Petrel process group:

| Petrel process | Status | What is missing | Plan |
|---|---|---|---|
| Seismic interpretation | ● partial | well tie **[done C40]** (synthetic, bulk shift, calibrated velocity); still missing: stretch / squeeze, extracted (statistical) wavelets, multi-horizon / fault interpretation sessions, guided fault auto-tracking | P2: fault likelihood attribute + auto fault extraction |
| Velocity modelling | ● basic | interval velocities per zone, well-calibrated V0 maps, time–depth from checkshots | P1: layer-cake velocity model with V0 maps from tops |
| Structural framework | ● partial | true fault network (fault–fault truncation, Y/λ junctions), horizon–fault contacts with throw consistency | P1: fault-network modeller (implicit fault surfaces + truncation rules) |
| Pillar / stair-step gridding | ● good | stair-step option for complex networks, local grid refinement in geomodel | P2 |
| Make horizons / zones | ● good | isochore from well data with trend, conformance rules (erosional / base-lap) | P2 |
| Well correlation | ● view only | interactive top picking with flattening, multi-well panels with cross-plots | P1 |
| Petrophysics / log calculator | ● good **[done C40]** | log calculator with presets (Vsh linear / Larionov, φ density / sonic, Archie Sw, Timur k, net flag, AI), cross-plots coloured by well / curve with zone filter and regression → poro-perm transform; still missing: Simandoux / dual-water Sw, multi-mineral solver | P3 |
| Data analysis | ● good **[done C40]** | interactive variogram fitting (draggable model, weighted auto fit, apply to steps) and histograms are in; still missing: declustering, vertical proportion curves, variogram maps | P2 |
| Facies / property modelling | ● good | trend maps and seismic-guided probability cubes per facies, filter by zone | P2 |
| Uncertainty & optimisation | ● partial (in Studies) | geomodel parameters (contacts, variogram ranges, seeds) as Study variables → structural uncertainty workflow | P1: expose geomodel parameters to Studies **[done C40: variogram range multiplier + realisation seed on the case's property models; contacts were already Study variables]** |
| Volumetrics | ● good | per-segment contacts, report export | P3 |
| Workflow editor | ✕ | record the 10 steps as a re-playable workflow with parameters | P2 (makes structural uncertainty possible) |
| Simulation link | ● good | — (Send to Model, GRDECL, MSW wells from trajectories) | — |

Recommended order: (1) log calculator + cross-plots, (2) interactive variogram analysis, (3) synthetic seismogram /
well tie, (4) geomodel parameters in Studies (structural & property uncertainty → HM), (5) fault-network modeller,
(6) workflow recorder.

---

## 5. Order of work after C39

1. G3 feature registry + palette + "Add physics" menu.
2. G4 object tree in the Case Explorer (wells first, then faults/regions, results).
3. Geomodel: log calculator + cross-plots, interactive variogram fitting.
4. G7 new-case wizard and empty states; G9 linked time cursor.
5. Geomodel: well tie, Studies integration; G10 case manager.

Status after C40: items 1–4 done; in 5 the well tie and the Studies integration are done. G10 case manager, G11 Inspector
filter, G12 status icons, the Appearance menu and the light 3-D background were also done in C40. C41 completed G5, the fault-network modeller (Petrel 5) and the workflow recorder (Petrel 6) and started G13 (Spanish
catalogue for menus, toolbar and workspace switch). Open: G13 coverage of the remaining panels.
C42: the geomodel import was rebuilt around column assignment (wizard, several files, well manager, 2-D lines). Surface /
property calculators, the 3-D value filter and statistics were added. Long help texts moved behind an (i) icon (*Appearance ›
Help text* restores the inline text). The deeper Petrel comparison and the next plan are in `docs/PETREL_PARITY.md`.
Screenshots: `docs/img/c42_gw-import.jpg`, `c42_gw-wells.jpg`, `c42_geo-scalc.jpg`, `c42_geo-stats.jpg`.

## 6. C39 verification screenshots (Linux, OpenGL back end, Mesa software renderer)

| | |
|---|---|
| Model workspace, dark theme | `docs/img/c39_linux_dark_model.jpg` |
| Model workspace, light theme + teal accent | `docs/img/c39_linux_light_model.jpg` |
| Geomodel 3-D scene (seismic slices, horizons, fault, wells, tops) | `docs/img/c39_linux_geomodel3d.jpg` |
| Crash-recovery dialog | `docs/img/c39_recovery_dialog.jpg` |

Produced by `ResSimStudio --smoke-test <frames> <png> <workspace> [--geo-demo]` under `xvfb-run` with `LIBGL_ALWAYS_SOFTWARE=1`.

## 7. C40 verification screenshots (Linux, OpenGL back end, Mesa software renderer, `--smoke-action`)

| | |
|---|---|
| New-model wizard over the Model workspace; object tree in the Case Explorer (G7, G4) | `docs/img/c40_wizard_objects.jpg` |
| Geomodel › Data analysis › Variogram: vertical PHIE variogram of the reservoir zone, auto fit (nugget 0.21, range 7.4 m) | `docs/img/c40_variogram_fit.jpg` |
| Geomodel › Data analysis › Well tie: W1, correlation 0.94 → 0.97, +2 ms bulk shift, calibrated v(z) | `docs/img/c40_well_tie.jpg` |
| Results after a run: linked report-time cursor on every time plot (G9) | `docs/img/c40_time_cursor.jpg` |
| Case manager: base + 2 variants, line diff of *High permeability* against the base (G10) | `docs/img/c40_case_manager.jpg` |
| Inspector filter "perm" on Grid & rock; ✓ / ⚠ section icons in the Case Explorer (G11, G12) | `docs/img/c40_inspector_filter.jpg` |
| 3-D background *Light grey* (View › Appearance) | `docs/img/c40_light_3d.jpg` |

