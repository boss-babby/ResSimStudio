// C42 geomodel import UI: column-assignment import wizard (files / pasted text, several files at once), well manager
// (manual heads, per-well trajectory and tops editors), 2-D seismic lines (loading, section view, picking, mis-ties,
// gridding). Implemented in ui_geoimport.cpp.
#pragma once
#include "geo_state.h"
#include "imgui.h"
#include <string>
#include <vector>

// file dialog (several files) -> import wizard for 'kind'; targetWell: import into this well only (deviation / tops)
void geoOpenImport(App::GeoState& S, ressim::geo::ImportKind kind, const std::string& targetWell = std::string());
// pasted / clipboard text -> import wizard
void geoOpenImportText(App::GeoState& S, ressim::geo::ImportKind kind, const std::string& text, const std::string& targetWell = std::string());
void geoImportWindow(App::GeoState& S);           // draws the wizard while open (call every frame)
void geoWellManager(App::GeoState& S);            // draws the well manager window while open
// LAS files (several): well name from the ~Well section or the file name; returns the number of logs attached
int geoImportLas(App::GeoState& S, const std::vector<std::string>& paths);
// 2-D lines (already read) appended to the project
void geoAddLines(App::GeoState& S, std::vector<ressim::geo::SeismicLine> lines);
void geoLinesPanel(App::GeoState& S);             // Seismic step: lines, horizons, picking, mis-ties, gridding
void geoLineSection(App::GeoState& S, ImVec2 size);   // centre view: the current 2-D line
void geoSegyBytesUi(App::GeoState& S);            // trace-header byte locations (2-D / 3-D readers)
// C42 Analysis view tabs: surface calculator, property calculator, statistics (ui_geocalc.cpp)
void geoCalcTabs(App::GeoState& S, int forceTab = -1);   // forceTab 0..2 selects a tab once (smoke tests)
// C43 geomodel tools (ui_geotools.cpp)
void geoMapOverlay(App::GeoState& S);             // inside the map plot: polygons / paths, digitising (pickMode 4 / 5)
void geoPolygonsPanel(App::GeoState& S);          // Data step
// C46 project file (.rgeo; ui_geoproject.cpp)
ressim::geo::Project geoToProject(const App::GeoState& S);
void geoFromProject(App::GeoState& S, ressim::geo::Project&& P, const std::string& path);
bool geoSaveProject(App::GeoState& S, bool saveAs);
bool geoOpenProject(App::GeoState& S, const std::string& path);
void geoProjectBar(App::GeoState& S);             // Data step: open / save / save as
void geoVolumeLimitsUi(App::GeoState& S);         // Volumetrics step: polygon limit, contacts per fault block
void geoPrepareVolSpec(const App::GeoState& S, ressim::geo::VolSpec& v);
void geoMistiePanel(App::GeoState& S);            // Seismic step: 2-D mis-tie correction
void geoLayerCakePanel(App::GeoState& S);         // Seismic step: layer-cake velocity model, checkshots
void geoSeisCalcPanel(App::GeoState& S);          // Seismic step: seismic calculator, crop / decimate
void geoIntersectionView(App::GeoState& S, ImVec2 size);   // centre view 6: arbitrary intersection along a path
// C43 P2 templates / P7 shared cell filter
const App::GeoState::PropTemplate* geoTemplate(App::GeoState& S, const std::string& propName);
void geoTemplateRange(const App::GeoState::PropTemplate* t, double& lo, double& hi);   // data range -> displayed (fixed / log10)
double geoTemplateValue(const App::GeoState::PropTemplate* t, double v);              // log10 when the template is logarithmic
int geoPlotColormap(int cmap, bool reverse);                                          // ImPlotColormap of a colormaps.h map
ImU32 geoTemplateColor(const App::GeoState::PropTemplate* t, double v, double lo, double hi);   // v, lo, hi displayed values
bool geoCellVisible(const App::GeoState& S, size_t cell);                             // natural cell index
void geoFilterButton(App::GeoState& S);                                               // view toolbar: filter popup + templates
void geoTemplatesWindow(App::GeoState& S);
void geoObjectTree(App::GeoState& S);   // C43 P1: input / model / result objects with context menus (left column)
void geoWellSection(App::GeoState& S, ImVec2 size);   // C43 P5: centre view 2 (well section)
