# Petrel parity: deep gap analysis and plan (C42)

The question: why does the Geomodel workspace feel "far behind Petrel" for **customisation, modification and
analysis**, and what does it take to close the gap? The C39 review in `UI_REVIEW.md` §4 compared the
*processes* (gridding, facies, petrophysics ...). This document goes one level deeper. It looks at the **interaction
model** that makes Petrel feel customisable, and it scores every Petrel domain on four axes:

* **I**mport: getting data in.
* **E**dit: modifying data by hand.
* **A**nalyse: statistics, filters, calculators, plots.
* **C**ustomise: styles, templates, settings per object, scripting.

Scores run from 0 (missing) to 3 (Petrel-class). An `→` marks a C42 change.

## 1. Why Petrel feels customisable: six mechanisms

| # | Petrel mechanism | What the user experiences | ResSim before C42 |
|---|---|---|---|
| M1 | **Every item is an object** in the Input / Models / Results trees. Each has a *Settings* dialog (Info, Style, Colours, Statistics, Histogram, Operations, History) and a right-click menu | "I can always find the thing and change how it looks or what it is" | Geomodel data lived in step panels. The only per-object actions were deleting a horizon or a well |
| M2 | **Templates**: a property / log template (unit, colour table, range, log/linear) is shared by every object of that kind | change once, everything updates; consistent maps | a fixed colour map per view |
| M3 | **Calculators everywhere**: surface calculator, property calculator, log calculator, seismic calculator, point-attribute calculator, each with IF, math and filters | arbitrary modifications without waiting for a feature | log calculator only (C40) |
| M4 | **Filters**: zone / segment / facies / value filters on grids, well filters, seismic cropping. They apply to all windows | "show me cells with φ > 0.2 in zone 2" | I/J/K ranges in the 3-D scene |
| M5 | **Many synchronised windows**: Map, Intersection (arbitrary polyline through seismic, grid, wells and surfaces), Well section (templated tracks), Histogram, Function, Cross-plot, 3-D | analysis from every angle with the same objects | one view at a time; well correlation was view-only |
| M6 | **Interactive editing tools**: surface / point / polygon / fault editing, make polygons, boundaries, well top picking in the well section, manual well and trajectory entry | "modify the model by hand" | horizon brush (C36), fault-stick digitising. No manual wells, no polygons |

Plus the plumbing: CRS and unit systems, the workflow editor and uncertainty, and Python / Ocean plug-ins.

## 2. Domain scorecard (I / E / A / C, 0–3)

| Domain | I | E | A | C | Gap (what Petrel has that we do not) |
|---|---|---|---|---|---|
| Well data (heads, surveys, tops, logs) | 1 → **3** | 0 → **2** | 2 | 1 | Before C42: fixed columns, one file at a time, no manual wells or surveys. **C42**: column-assignment import, several files at once, well manager, per-well survey and tops editors, paste from Excel. Still missing: well attributes (symbol, status, UWI fields), checkshots, well filters (C46: well folders in the object tree) |
| Surfaces & points | 2 | 2 | 1 → **2** | 1 | **C42**: surface calculator and statistics. Still missing: polygons and boundaries as objects, point attributes, surface operations menu (smooth / eliminate / clip by polygon as one-click) |
| Seismic 3-D | 2 | 1 | 2 | 1 | byte-location dialog → **C42**. Still missing: cropping / realisation of sub-volumes, seismic calculator, colour templates |
| Seismic 2-D | 0 → **2** | 0 → **2** | 0 → **2** | 1 | **C42**: several SEG-Y lines, section view, seeded tracking / manual picks / erasing, intersections with mis-ties, gridding to depth, 3-D curtains. Still missing: mis-tie *correction* (bulk shift / phase per line), 2-D fault interpretation, line-to-3-D ties |
| Interpretation | 1 | 2 | 1 | 1 | multi-horizon sessions, interpretation folders, fault interpretation with auto-extraction |
| Velocity / depth conversion | 1 | 1 | 1 | 0 | layer-cake model with per-zone V0 / k maps, checkshot import, residual correction at the wells |
| Structural framework & gridding | — | 2 | 1 | 1 | fault network done (C41). Still missing: stair-step option, zone conformance rules, LGR in the geomodel |
| Property modelling | — | 1 → **2** | 2 | 1 | **C42**: property calculator (new property or overwrite φ / k / facies). Still missing: trend maps per facies, VPC, declustering |
| Data analysis | — | — | 2 | 1 | variogram maps, VPC, transforms pipeline (normal score, trends) per zone |
| Grid filters / views | — | — | 1 → **2** | 1 | **C42**: value filter in the 3-D scene. Still missing: zone / segment filters shared by all windows, arbitrary intersection window |
| Volumetrics | — | — | 2 | 1 | per-segment / per-zone contacts table, report export with templates |
| Workflow & uncertainty | — | — | 2 | 2 | done C40–C41 (recorder, Studies). Still missing: workflow variables bound to every object parameter |
| Customisation (M1 / M2) | — | — | — | 1 | object settings dialogs and templates are the largest structural gap (see §4) |

## 3. Implemented in C42 (this round)

1. **Column-assignment import** (`geo_import.h`, `ui_geoimport.cpp`). Delimiter / header / units-row / decimal-comma /
   null / Petrel `BEGIN HEADER` sniffing. A role + unit factor per column, guessed from about 150 synonyms and the column
   order, always user-overridable and remembered per layout. Typed builders for heads, surveys (MD/INC/AZI, X/Y/TVD(SS)/Z,
   DX/DY offsets, several wells per file), tops (MD or depth), point sets and fault sticks.
2. **Several files at once** for every importer (Windows multi-select, tinyfd on Linux / macOS). Name matching by
   well column, file header (`# WELL NAME:`) or file name (`A1_survey.txt` → A1, never A10).
3. **Well manager**: add / edit / delete wells. Each well has a trajectory editor (MD/INC/AZI or MD/X/Y/TVDSS, insert /
   delete stations, paste from Excel, import a file into this well, vertical to TD, copy), plan and section previews,
   a tops editor and a logs list.
4. **2-D seismic**: several SEG-Y lines with 2-D / 3-D auto-detection and configurable trace-header bytes (shared with
   the 3-D readers). Section view with seeded tracking / manual picks / erase, intersections and mis-tie table, gridding
   of the picks to a depth surface, curtains and picks in the 3-D scene, lines on the map.
5. **Surface calculator**: expressions over all horizons resampled to one lattice, with `if()`, min / max, X / Y.
6. **Property calculator**: expressions over the grid properties (PORO, PERMX, NTG, FACIES, ZONE, DEPTH, I, J, K). The
   result goes into a new user property or overwrites porosity / permeability / facies.
7. **Value filter** in the 3-D scene (hide cells outside a range of the displayed property) and **statistics** (n, min,
   max, mean, std, P10 / P50 / P90) for horizons and grid properties.
8. **Compact help**: explanatory text behind an (i) icon (hover / click to pin). Setting: Appearance → Help text.

## 4. Plan for the remaining gap (priority order)

| P | Item | Mechanism | Size | Notes |
|---|---|---|---|---|
| 1 | **Geomodel object tree** (Input: wells / surfaces / points / polygons / faults / seismic / lines; Models: grid, properties; Results: volumetrics). Each object has a right-click menu (settings, rename, delete, export, statistics, calculator) | M1 | L | reuse the G4 Case Explorer tree code; the step panels stay as the guided workflow |
| 2 | **Templates** (property kind → unit, colour map, range, log scale) used by map, sections, 3-D and correlation | M2 | M | one table in the project; colour maps already shared (colormaps.h) |
| 3 | **Intersection window**: arbitrary polyline (digitised on the map, along a well or a 2-D line) through seismic, grid, horizons, faults and wells | M5 | L | the 2-D line section code and the grid section are the building blocks |
| 4 | **Polygons / boundaries** as objects: digitise on the map, import / export, clip surfaces, limit volumetrics and gridding | M6 | M | |
| 5 | **Well section window**: templated tracks (logs, tops, facies, blocked logs), top picking by click, flatten on a top | M5 / M6 | L | the correlation view is the start |
| 6 | **Mis-tie correction** for 2-D lines (bulk shift per line by least squares over the intersections) and 2-D fault sticks | — | S | intersections already computed |
| 7 | Zone / segment / facies filters shared by all windows | M4 | M | |
| 8 | Layer-cake velocity model with checkshots and V0 maps | — | M | |
| 9 | Seismic calculator and sub-volume cropping | M3 | S | the attribute engine exists |
| 10 | Python scripting of the geomodel (pybind of the geo API, like PYACTION) | C | L | Ocean-equivalent extensibility |

Status after C43: all ten items of this plan are implemented (see §5). Items still open are listed at the end of §5.

## 5. Implemented in C43

| P | Status | Where | What exists now / what is still missing |
|---|---|---|---|
| 1 | **done** | `ui_geotools.cpp` geoObjectTree | Input / Models / Results tree in the left column, right-click menus (edit, show / hide in the well section, map, statistics, IRAP export, template settings, intersection, limit volumetrics, 3-D / line section, calculator, SEG-Y export, delete), rename; C46: user folders per object kind (Move to folder, rename / remove folder, drag and drop), saved in the `.rgeo` project (Open / Save project on the Data step) |
| 2 | **done** | GeoState templates | colour map, direction, fixed range and log scale per property kind, used by map, grid section, intersection and 3-D scene. Missing: templates in the well-section tracks (they use their own log flag) and template files shared between projects |
| 3 | **done** | geoIntersectionView | digitised path / polygon as section line; depth mode (grid property, horizons, wells within a distance) and seismic TWT mode. Missing: faults in the section, paths along a well or a 2-D line |
| 4 | **done** | `geo_tools.cpp` polygons | digitise, import / export (text), from the grid outline, clip a surface, limit the volumetrics, contacts per fault block |
| 5 | **done** | geoWellSection | tracks (curve, log scale, order), discrete curves as blocks, tops, flatten on a top, click top picking, hover values, show / hide wells. Missing: track width, blocked logs, well order |
| 6 | **done** (shifts) | `geo_tools.cpp` mis-tie | least-squares bulk time shift per line, one reference per connected group, applied to traces and picks. Missing: phase rotation, 2-D fault sticks |
| 7 | **done** | GeoState cf | zone / fault-block / facies filter shared by 3-D, grid section and (optionally) volumetrics |
| 8 | **done** | `geo_tools.cpp` layer cake | interval velocities from tops, mean / IDW map / user value per interface, residual correction, V0 + kZ checkshot fit. Missing: checkshot import per well from a file format other than the table |
| 9 | **done** | seismicCalculator / crop | expression over AMP, T, I, J and the attributes; crop / decimate; SEG-Y export |
| 10 | **done** | `capi_geo.cpp`, `python/ressim/geo.py` | C API `rs_geo_*` and Python `ressim.geo.Geomodel` (surfaces, wells, logs, polygons, layer cake, seismic calculator / crop, IRAP / SEG-Y output); example `python/examples/geo_scripting.py` |

The Geomodel workspace is also translated (View › Language, `--lang`).

Scorecard changes (§2): surfaces & points A 2 → 3; seismic 2-D A 2 → 3 (mis-tie correction); velocity I / E / A 1 → 2;
grid filters A 2 → 3; volumetrics A 2 → 3; customisation C 1 → 2 (object tree + templates).

Still open, in priority order:
1. faults and well / 2-D-line paths in the intersection view;
2. ~~rename / folders in the object tree~~ (done C46);
3. templates in the well-section tracks;
4. phase mis-tie;
5. 2-D fault sticks;
6. workflow variables bound to object parameters;
7. VPC / declustering / variogram maps.
