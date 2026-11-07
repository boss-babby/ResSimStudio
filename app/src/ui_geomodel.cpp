// C32 Geomodel workspace / ResGeo (geomodel.h): Petrel / RMS-style workflow from data to a simulation case.
// Left: workflow steps; centre: map / seismic section / well correlation / grid cross-section; right: step settings.
// Long computations (tracking, grid building, facies, petrophysics, realisations, upscaling) run on a worker thread.
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "platform.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "colormaps.h"
#include "plotkit.h"
#include "ressim/geomodel.h"
#include "ressim/examples.h"
#include <utility>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <future>
#include <map>


namespace {
const char* kSteps[] = {TRN(ICON_FOLDER_OPEN " Data"), TRN(ICON_ACTIVITY " Seismic"), TRN(ICON_MOUNTAIN " Structure"), TRN(ICON_LAYERS " Grid"), TRN(ICON_ROUTE " Blocking"),
                        TRN(ICON_GRID_3X3 " Facies"), TRN(ICON_DROPLET " Petrophysics"), TRN(ICON_SIGMA " Volumetrics"), TRN(ICON_SQUARE_STACK " Upscaling"), TRN(ICON_LOG_OUT " Export")};
const int kNSteps = 10;
ImPlotSpec lineSpec(ImVec4 c, float w) { ImPlotSpec s; s.LineColor = c; s.LineWeight = w; return s; }
ImPlotSpec markSpec(ImPlotMarker m, float sz, ImVec4 c) { ImPlotSpec s; s.Marker = m; s.MarkerSize = sz; s.MarkerFillColor = c; s.MarkerLineColor = c; s.LineColor = c; return s; }
std::string baseName(const std::string& p) { const size_t s = p.find_last_of("/\\"); std::string b = p.substr(s == std::string::npos ? 0 : s + 1); const size_t d = b.find('.'); return d == std::string::npos ? b : b.substr(0, d); }
// grid property per cell for display / export
std::vector<double> gridProperty(const App::GeoState& S, int which) {
    const GeoGrid& g = S.g;
    std::vector<double> v(g.bulk.size(), std::nan(""));
    const int nxy = g.grid.nx * g.grid.ny;
    for (size_t c = 0; c < v.size(); ++c) {
        if (!g.grid.actnum.empty() && !g.grid.actnum[c]) continue;
        switch (which) {
            case 0: v[c] = g.top[c]; break;
            case 1: if (S.hasFacies && c < S.fr.facies.size()) v[c] = S.fr.facies[c]; break;
            case 2: if (S.hasPetro && c < S.phi.v.size()) v[c] = S.phi.v[c]; break;
            case 3: if (S.hasPetro && c < S.kx.size() && S.kx[c] > 0) v[c] = std::log10(S.kx[c]); break;
            case 4: v[c] = g.block.empty() ? 0 : g.block[c % size_t(nxy)]; break;
            case 5: { const size_t k = c / size_t(nxy); v[c] = k < g.layerZone.size() ? g.layerZone[k] : 0; } break;
            case 6: v[c] = g.bot[c] - g.top[c]; break;
            default: {
                const size_t u = size_t(which - 7);
                if (which >= 7 && u < S.userProps.size() && c < S.userProps[u].second.size()) v[c] = S.userProps[u].second[c];
            } break;
        }
    }
    return v;
}
}  // namespace
std::vector<double> geoGridProperty(const App::GeoState& S, int which) { return gridProperty(S, which); }
const char* kGeoGridProps[7] = {"Top depth", "Facies", "Porosity", "log10 permeability X", "Fault block", "Zone", "Cell thickness"};
void geoLaunch(App::GeoState& S, const std::string& name, std::function<std::function<void(App::GeoState&)>()> f) {
    if (S.job.valid()) { S.jobQueue.emplace_back(name, std::move(f)); S.msg(name + " (queued)"); return; }
    S.jobName = name; S.jobStart = std::chrono::steady_clock::now();
    S.job = std::async(std::launch::async, std::move(f));
    S.msg(name + " ...");
}
int geoGridPropCount(const App::GeoState& S) { return 7 + int(S.userProps.size()); }
std::string geoGridPropName(const App::GeoState& S, int which) {
    if (which >= 0 && which < 7) return kGeoGridProps[which];
    const size_t u = size_t(which - 7);
    return which >= 7 && u < S.userProps.size() ? S.userProps[u].first : std::string("?");
}
bool geoGridPropCombo(const char* id, int* which, const App::GeoState& S) {
    const int n = geoGridPropCount(S);
    if (*which < 0 || *which >= n) *which = std::clamp(*which, 0, n - 1);
    bool ch = false;
    if (ImGui::BeginCombo(id, geoGridPropName(S, *which).c_str())) {
        for (int i = 0; i < n; ++i) {
            if (i == 7) ImGui::SeparatorText(TR("Calculated"));
            if (ImGui::Selectable(geoGridPropName(S, i).c_str(), i == *which)) { *which = i; ch = true; }
        }
        ImGui::EndCombo();
    }
    return ch;
}

void App::setGeoMode(bool on) { geoMode_ = on; if (on) workspace_ = Workspace::Geomodel; }

void App::geomodelWorkspace() {
    if (!geo_) {
        geo_ = std::make_shared<GeoState>();
        geo_->fs.vg.rangeMax = geo_->fs.vg.rangeMin = 800; geo_->fs.vg.rangeVert = 6;
        geo_->ps.vg.rangeMax = geo_->ps.vg.rangeMin = 1000; geo_->ps.vg.rangeVert = 8; geo_->ps.maxV = 0.4;
        geo_->vs.contacts.owc = 2080; geo_->vs.satHeight = false;
        geo_->us.fz = 3;
    }
    GeoState& S = *geo_;
    if (geoStepReq_ >= 0) { S.step = geoStepReq_; geoStepReq_ = -1; }   // C40 deep links from the command palette
    if (geoViewReq_ >= 0) { S.view = geoViewReq_; geoViewReq_ = -1; }
    // C42 smoke-test states (after the demo project has loaded): well manager / column-assignment wizard
    if ((smokeAction_ == "gw-wells" || smokeAction_ == "gw-import") && !S.wells.empty() && !S.job.valid()) {
        if (smokeAction_ == "gw-wells") { S.we.open = true; S.we.sel = std::min<int>(1, int(S.wells.size()) - 1); }
        else geoOpenImportText(S, ImportKind::Deviation,
                               "Well;Measured depth;Inclination;Azimuth\nB-7;0;0;0\nB-7;800;12,5;45\nB-7;1600;30;47\nB-8;0;0;0\nB-8;1200;8;310\n");
        smokeAction_.clear();
    }
    // C43 smoke-test states: polygons on the map, intersection view along a path, seismic calculator
    if (smokeAction_.rfind("gt-", 0) == 0 && !S.wells.empty() && S.hasCube && !S.job.valid()) {
        if (smokeAction_ == "gt-poly" || smokeAction_ == "gt-isec") {
            ressim::geo::Polygon P; P.name = "Lease A";
            P.pts = {{600, 500}, {2600, 400}, {3000, 1800}, {1400, 2500}, {500, 1600}};
            S.polys.push_back(P);
            ressim::geo::Polygon Q; Q.name = "Path 1"; Q.closed = false;
            Q.pts = {{300, 300}, {1200, 1500}, {2400, 1700}, {3700, 2700}};
            S.polys.push_back(Q);
            S.curPoly = 0; S.volPoly = 0; S.isecPath = 1;
            S.step = 0; S.openTool = 0; S.view = smokeAction_ == "gt-poly" ? 0 : 6;
        } else if (smokeAction_ == "gt-isecw") {   // C44: intersection along the deviated well's trajectory
            for (int w = 0; w < int(S.wells.size()); ++w) if (!S.wells[size_t(w)].traj.empty()) S.isecPath = int(S.polys.size()) + w;
            S.step = 0; S.view = 6;
        } else if (smokeAction_ == "gt-wsec2") {   // C44: well order, track widths
            S.view = 2; S.wsOrder = {"W6", "W3"}; S.wsHidden = {"W4", "W5"};
            App::GeoState::WsTrack a; a.curve = "GR"; a.width = 2.0f; App::GeoState::WsTrack b; b.curve = "PHIE"; App::GeoState::WsTrack c; c.curve = "PERM"; c.width = 0.6f;
            S.wsTracks = {a, b, c};
        } else if (smokeAction_ == "gt-wsec") {
            S.view = 2; S.wsFlat = "TopRes"; S.wsPick = true; S.wsPickTop = "MidRes";
        } else if (smokeAction_ == "gt-proj") {   // C46: folders, project save -> open round trip in the app
            using ressim::geo::ObjKind;
            S.folders.assign(ObjKind::Well, S.wells[0].name, "North block");
            if (S.wells.size() > 2) S.folders.assign(ObjKind::Well, S.wells[2].name, "North block");
            if (!S.hz.empty()) S.folders.assign(ObjKind::Surface, S.hz[0].name, "Interpreted");
            const std::string p = (std::filesystem::temp_directory_path() / "resgeo_smoke.rgeo").string();
            const size_t nw = S.wells.size(), nh = S.hz.size(), nf = S.faults.size(); const int ni = S.cube.ni;
            S.projectPath = p;
            const bool saved = geoSaveProject(S, false);
            const bool opened = saved && geoOpenProject(S, p);
            std::fprintf(stderr, "gt-proj: saved %d opened %d: wells %zu->%zu surfaces %zu->%zu faults %zu->%zu cube ni %d->%d (%s); folder of %s = '%s'\n",
                         int(saved), int(opened), nw, S.wells.size(), nh, S.hz.size(), nf, S.faults.size(), ni, S.cube.ni, S.cubePath.c_str(),
                         S.wells.empty() ? "-" : S.wells[0].name.c_str(), S.wells.empty() ? "" : S.folders.folderOf(ObjKind::Well, S.wells[0].name).c_str());
            S.treeExpand = true; S.view = 0;
        } else if (smokeAction_ == "gt-scalc") {
            S.step = 1; S.openTool = 2; S.view = 1; S.secAttr = 0;
            geoLaunch(S, "Seismic calculator", [c = S.cube]() -> std::function<void(GeoState&)> {
                SeismicCube o; std::string err;
                const bool ok = seismicCalculator(c, "ENVELOPE * (T > 1700)", 6, o, err);
                return [ok, o, err](GeoState& S) { if (ok) { S.cube = o; ++S.cubeRev; S.msg("seismic calculator: ENVELOPE * (T > 1700)"); } else S.msg("Error: " + err); };
            });
        }
        smokeAction_.clear();
    }
    // ---- finished background job?
    if (S.job.valid() && S.job.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        try { auto apply = S.job.get(); if (apply) apply(S); }
        catch (const std::exception& e) { S.msg(std::string("Error: ") + e.what()); notify(std::string("Geomodel: ") + e.what()); }
    }
    if (!S.job.valid() && !S.jobQueue.empty()) {   // C43: next queued job
        auto q = std::move(S.jobQueue.front()); S.jobQueue.pop_front();
        geoLaunch(S, q.first, std::move(q.second));
    }
    const bool busy = S.job.valid();
    auto launch = [&](const std::string& name, std::function<std::function<void(GeoState&)>()> f) {
        if (busy) return;   // buttons are disabled while busy; several launches in one frame are queued by geoLaunch
        geoLaunch(S, name, std::move(f));
    };

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
    ImGui::Begin("##geomodel", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings |
                                            ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoMove);
    ImGui::PopStyleVar(2);
    const float fs = ImGui::GetFontSize();
    const float leftW = fs * 12.5f, rightW = fs * 23.0f;

    // ================================================================ left: workflow steps
    ImGui::BeginChild("##geosteps", ImVec2(leftW, 0), ImGuiChildFlags_Borders);
    ImGui::TextDisabled("%s", geoMode_ ? TR("RESGEO WORKFLOW") : TR("GEOMODEL WORKFLOW"));
    ImGui::Separator();
    const bool done[kNSteps] = {!S.hz.empty() || !S.wells.empty(), S.hasTh, S.hz.size() >= 2, S.hasGrid, !S.blocked.empty(), S.hasFacies, S.hasPetro, S.hasVol, S.hasUp, false};
    for (int i = 0; i < kNSteps; ++i) {
        ImGui::PushID(i);
        char b[96]; std::snprintf(b, sizeof b, "%d  %s", i + 1, TR(kSteps[i]));
        if (ImGui::Selectable(b, S.step == i, 0, ImVec2(0, fs * 1.7f))) S.step = i;
        if (done[i]) { ImGui::SameLine(leftW - fs * 2.2f); ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(Pal::Ok), ICON_CHECK); }
        ImGui::PopID();
    }
    ImGui::Separator();
    geoObjectTree(S);   // C43 P1: project objects with right-click menus (replaces the C32 summary lines)
    ImGui::TextDisabled("%s: %s", TR("units"), S.field ? "FIELD (ft)" : "METRIC (m)");
    if (busy) {
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - S.jobStart).count();
        ImGui::Separator();
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(Pal::Warn), "%s  %.1f s", S.jobName.c_str(), t);
    }
    if (geoMode_) {
        ImGui::Separator();
        if (ImGui::Button(TR(ICON_BOX " Open ResSim Studio"))) { geoMode_ = false; workspace_ = Workspace::Model; }
        ui::hint("Leave the ResGeo-only mode: model, run and results workspaces");
    }
    ImGui::EndChild();
    ImGui::SameLine();

    // ================================================================ centre: views
    const float centreW = ImGui::GetContentRegionAvail().x - rightW - ImGui::GetStyle().ItemSpacing.x;
    ImGui::BeginChild("##geoview", ImVec2(centreW, 0), ImGuiChildFlags_Borders);
    static const char* views[] = {TRN(ICON_IMAGE " Map"), TRN(ICON_ACTIVITY " Seismic section"), TRN(ICON_CHART_LINE " Well section"), TRN(ICON_LAYERS " Grid cross-section"), TRN(ICON_BOX " 3D scene"), TRN(ICON_CHART_SCATTER " Data analysis"), TRN(ICON_SCAN_LINE " Intersection")};
    for (int i = 0; i < 7; ++i) {
        if (i) ImGui::SameLine();
        const bool on = S.view == i;
        if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::ColorConvertU32ToFloat4(Pal::Accent));
        if (ImGui::Button(TR(views[i]))) S.view = i;
        if (on) ImGui::PopStyleColor();
    }
    ImGui::SameLine(0, fs * 1.5f); geoFilterButton(S);   // C43 P7 / P2: shared filter, templates
    geoTemplatesWindow(S);
    ImGui::Separator();
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const ImVec2 plotSize(avail.x, avail.y - fs * 7.5f);
    if (S.view == 0) {   // ---------------- map
        std::vector<std::string> names; for (auto& h : S.hz) names.push_back(h.name);
        if (S.hasGrid) names.push_back("Grid layer property");
        if (!names.empty()) {
            S.mapSurface = std::clamp(S.mapSurface, 0, int(names.size()) - 1);
            ImGui::SetNextItemWidth(fs * 12); ImGui::Combo("##ms", &S.mapSurface, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &names, int(names.size()));
            if (S.hasGrid && S.mapSurface == int(S.hz.size())) {
                ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 10); geoGridPropCombo("##gp", &S.gridProp, S);
                ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8); ImGui::SliderInt("layer", &S.gridLayer, 0, S.g.grid.nz - 1);
            }
        }
        if (ImPlot::BeginPlot("##map", plotSize, ImPlotFlags_Equal)) {
            ImPlot::SetupAxes("x (east)", "y (north)");
            double lo = 0, hi = 1;
            if (S.mapSurface < int(S.hz.size()) && S.mapSurface >= 0) {
                const Surface& h = S.hz[size_t(S.mapSurface)];
                h.range(lo, hi);
                const auto* tpl = geoTemplate(S, h.name);   // C43 P2
                geoTemplateRange(tpl, lo, hi);
                std::vector<float> img(size_t(h.nx) * size_t(h.ny));
                for (int j = 0; j < h.ny; ++j) for (int i = 0; i < h.nx; ++i) { const double v = geoTemplateValue(tpl, h.node(i, j)); img[size_t((h.ny - 1 - j) * h.nx + i)] = std::isfinite(v) ? float(v) : float(hi); }
                ImPlot::PushColormap(geoPlotColormap(tpl ? tpl->cmap : 1, tpl ? tpl->reverse : false));
                // rotation is not drawn (heatmap in the lattice frame at the origin)
                ImPlot::PlotHeatmap("##hz", img.data(), h.ny, h.nx, lo, hi, nullptr, ImPlotPoint(h.x0, h.y0), ImPlotPoint(h.x0 + (h.nx - 1) * h.dx, h.y0 + (h.ny - 1) * h.dy));
                ImPlot::PopColormap();
            } else if (S.hasGrid) {
                const auto v = gridProperty(S, S.gridProp);
                const int nx = S.g.grid.nx, ny = S.g.grid.ny, k = std::clamp(S.gridLayer, 0, S.g.grid.nz - 1);
                std::vector<float> img(size_t(nx) * size_t(ny));
                lo = 1e300; hi = -1e300;
                for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { const double x = v[size_t(i + nx * (j + ny * k))]; if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); } }
                if (lo > hi) { lo = 0; hi = 1; }
                const auto* tpl = geoTemplate(S, geoGridPropName(S, S.gridProp));   // C43 P2 template, P7 filter
                geoTemplateRange(tpl, lo, hi);
                for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
                    const size_t c = size_t(i + nx * (j + ny * k));
                    const double x = geoCellVisible(S, c) ? geoTemplateValue(tpl, v[c]) : NAN;
                    img[size_t((ny - 1 - j) * nx + i)] = float(std::isfinite(x) ? x : lo);
                }
                ImPlot::PushColormap(geoPlotColormap(tpl ? tpl->cmap : 1, tpl ? tpl->reverse : false));
                const GridSpec& G = S.g.spec;   // the spec the grid was built with (rotation not drawn)
                ImPlot::PlotHeatmap("##gl", img.data(), ny, nx, lo, hi, nullptr, ImPlotPoint(G.x0, G.y0), ImPlotPoint(G.x0 + nx * G.dx, G.y0 + ny * G.dy));
                ImPlot::PopColormap();
            }
            // fault traces at mid-depth, grid outline, wells
            const double zm = S.zMid();
            for (const auto& f : S.faults) {
                const auto tr = f.traceAt(zm);
                std::vector<double> xs, ys; for (auto& p : tr) { xs.push_back(p[0]); ys.push_back(p[1]); }
                ImPlot::PlotLine(("fault " + f.name).c_str(), xs.data(), ys.data(), int(xs.size()), lineSpec(ImVec4(1, 0.25f, 0.2f, 1), 2.5f));
            }
            if (!S.fnJ.empty()) {   // C41 fault-network junctions
                std::vector<double> jx, jy;
                for (const auto& J : S.fnJ) { jx.push_back(J.x); jy.push_back(J.y); }
                ImPlot::PlotScatter("fault junctions", jx.data(), jy.data(), int(jx.size()), markSpec(ImPlotMarker_Diamond, 6, ImVec4(1, 0.85f, 0.2f, 1)));
            }
            if (S.hasGrid || S.gs.nx > 0) {
                GeoGrid tmp; tmp.spec = S.gs;
                double xs[5], ys[5];
                const int ci[5] = {0, S.gs.nx, S.gs.nx, 0, 0}, cj[5] = {0, 0, S.gs.ny, S.gs.ny, 0};
                for (int q = 0; q < 5; ++q) tmp.cornerXY(ci[q], cj[q], xs[q], ys[q]);
                ImPlot::PlotLine("grid outline", xs, ys, 5, lineSpec(ImVec4(1, 1, 1, 0.7f), 1.0f));
            }
            for (const auto& w : S.wells) {
                std::vector<double> xs, ys;
                if (w.traj.empty()) { xs.push_back(w.x); ys.push_back(w.y); }
                else for (auto& s : w.traj) { xs.push_back(s.x); ys.push_back(s.y); }
                if (xs.size() > 1) ImPlot::PlotLine(("##t" + w.name).c_str(), xs.data(), ys.data(), int(xs.size()), lineSpec(ImVec4(1, 1, 1, 1), 2.0f));
                ImPlot::PlotScatter(("##w" + w.name).c_str(), &w.x, &w.y, 1, markSpec(ImPlotMarker_Circle, 5, ImVec4(1, 1, 1, 1)));
                ImPlot::PlotText(w.name.c_str(), w.x, w.y, ImVec2(0, -12));
            }
            for (size_t l = 0; l < S.lines.size(); ++l) {   // C42 2-D lines in plan
                if (l < S.lineOn.size() && !S.lineOn[l]) continue;
                const SeismicLine& L = S.lines[l];
                if (L.ntr() < 2) continue;
                ImPlot::PlotLine(("##l2" + L.name).c_str(), L.x.data(), L.y.data(), L.ntr(), lineSpec(int(l) == S.curLine ? ImVec4(1, 0.85f, 0.2f, 1) : ImVec4(0.55f, 0.8f, 1, 0.9f), int(l) == S.curLine ? 3.0f : 1.8f));
                ImPlot::PlotText(L.name.c_str(), L.x[0], L.y[0], ImVec2(0, -10));
            }
            // C36-4 horizon brush on the displayed horizon (the heatmap is drawn in the unrotated lattice frame)
            if (S.pickMode == 3 && S.mapSurface >= 0 && S.mapSurface < int(S.hz.size()) && ImPlot::IsPlotHovered()) {
                S.editHz = S.mapSurface;
                const Surface& h = S.hz[size_t(S.editHz)];
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                double cx[49], cy[49];
                for (int q = 0; q <= 48; ++q) { const double a = 2 * 3.14159265358979 * q / 48; cx[q] = mp.x + S.brushR * std::cos(a); cy[q] = mp.y + S.brushR * std::sin(a); }
                ImPlot::PlotLine("##brush", cx, cy, 49, lineSpec(ImVec4(1, 1, 1, 0.9f), 2.0f));
                double x, y; geoLatticeXY(h, (mp.x - h.x0) / h.dx, (mp.y - h.y0) / h.dy, x, y);
                const double z = h.at(x, y);
                if (std::isfinite(z)) ImGui::SetTooltip("%s depth %.2f\nclick: %s (R %.0f)", h.name.c_str(), z, S.brushMode == 0 ? "raise / lower" : S.brushMode == 1 ? "smooth" : "flatten", S.brushR);
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && ImGui::GetIO().MouseDragMaxDistanceSqr[0] < 16.0f) geoApplyBrush(S, x, y);
            }
            geoMapOverlay(S);   // C43 polygons / intersection paths, digitising (pickMode 4 / 5)
            ImPlot::EndPlot();
        }
    } else if (S.view == 1) {   // ---------------- seismic section (C36: inline / crossline, attributes, flattening, picking)
        geoSeismicSection(S, plotSize);
    } else if (S.view == 4) {   // ---------------- 3-D scene (C36)
        geomodel3D(S, plotSize);
    } else if (S.view == 5) {   // ---------------- C40 data analysis (log calculator, cross-plot, histogram, variogram fit, well tie)
        geoAnalysis(S, ImVec2(avail.x, avail.y - fs * 0.5f));
    } else if (S.view == 6) {   // ---------------- C43 arbitrary intersection along a path
        geoIntersectionView(S, plotSize);
    } else if (S.view == 2) {   // ---------------- C43 P5 well section (tracks, flattening, top picking)
        geoWellSection(S, plotSize);
    } else {   // ---------------- grid cross-section along a row
        if (!S.hasGrid) ImGui::TextDisabled("Build the grid first (Grid step).");
        else {
            const GeoGrid& g = S.g;
            const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
            S.sectionRow = std::clamp(S.sectionRow, 0, ny - 1);
            ImGui::SetNextItemWidth(fs * 12); ImGui::SliderInt("row J", &S.sectionRow, 0, ny - 1);
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 10); geoGridPropCombo("##gp2", &S.gridProp, S);
            const auto v = gridProperty(S, S.gridProp);
            double lo = 1e300, hi = -1e300, zlo = 1e300, zhi = -1e300;
            for (int k = 0; k < nz; ++k) for (int i = 0; i < nx; ++i) {
                const size_t c = size_t(i + nx * (S.sectionRow + ny * k));
                if (std::isfinite(v[c])) { lo = std::min(lo, v[c]); hi = std::max(hi, v[c]); }
                zlo = std::min(zlo, g.top[c]); zhi = std::max(zhi, g.bot[c]);
            }
            if (lo > hi) { lo = 0; hi = 1; }
            const auto* tpl = geoTemplate(S, geoGridPropName(S, S.gridProp));   // C43 P2
            geoTemplateRange(tpl, lo, hi);
            if (ImPlot::BeginPlot("##xs", plotSize)) {
                ImPlot::SetupAxes("distance along the row", "depth", 0, ImPlotAxisFlags_Invert);
                ImPlot::SetupAxesLimits(0, nx * g.spec.dx, zlo - 5, zhi + 5, ImGuiCond_Always);
                ImPlot::PushPlotClipRect();
                ImDrawList* dl = ImPlot::GetPlotDrawList();
                const int j = S.sectionRow;
                for (int k = 0; k < nz; ++k)
                    for (int i = 0; i < nx; ++i) {
                        const size_t c = size_t(i + nx * (j + ny * k));
                        if (!std::isfinite(v[c]) || !geoCellVisible(S, c)) continue;
                        const double tv = geoTemplateValue(tpl, v[c]);
                        if (!std::isfinite(tv)) continue;
                        auto zc = [&](int f, int a) { return 0.5 * (g.grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j) * 2 * nx + 2 * i + a)] + g.grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j + 1) * 2 * nx + 2 * i + a)]); };
                        const double x0 = i * g.spec.dx, x1 = (i + 1) * g.spec.dx;
                        const ImVec2 p[4] = {ImPlot::PlotToPixels(x0, zc(0, 0)), ImPlot::PlotToPixels(x1, zc(0, 1)), ImPlot::PlotToPixels(x1, zc(1, 1)), ImPlot::PlotToPixels(x0, zc(1, 0))};
                        dl->AddConvexPolyFilled(p, 4, geoTemplateColor(tpl, tv, lo, hi));
                    }
                ImPlot::PopPlotClipRect();
                ImPlot::EndPlot();
            }
            ImGui::Text("%s: %.4g .. %.4g", geoGridPropName(S, S.gridProp).c_str(), lo, hi);
        }
    }
    // log
    ImGui::BeginChild("##geolog", ImVec2(0, 0), ImGuiChildFlags_Borders);
    for (const auto& l : S.log) ImGui::TextUnformatted(l.c_str());
    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4) ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
    ImGui::EndChild();
    ImGui::SameLine();

    // ================================================================ right: step settings
    ImGui::BeginChild("##geoprops", ImVec2(0, 0), ImGuiChildFlags_Borders);
    ImGui::BeginDisabled(busy);
    ui::header(TR(kSteps[S.step]));
    if (S.step == 0) {   // ---------------- data
        const bool smokeDemo = std::exchange(smokeGeoDemo_, false);
        if (smokeDemo) S.view = smokeAction_.rfind("geo-", 0) == 0 ? 5 : 4;   // C40: --smoke-action geo-* opens Data analysis
        if (ui::accentButton(ICON_PLAY " Load the demo project") || smokeDemo) {
            const bool fieldU = S.field;
            launch("Writing and loading the demo project", [fieldU]() -> std::function<void(GeoState&)> {
                const std::string dir = (std::filesystem::temp_directory_path() / (fieldU ? "ResGeo_demo_field" : "ResGeo_demo")).string();
                std::string e;
                const auto files = writeDemoProject(dir, fieldU, e);
                std::vector<Surface> hz; std::vector<Fault> fl; std::vector<GeoWell> ws; SeismicCube cube;
                for (const char* h : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(dir + "/" + h + ".irap", e); s.name = h; hz.push_back(s); }
                fl = readFaultSticks(dir + "/F1_sticks.txt", e);
                ws = readWellHeads(dir + "/well_heads.txt", e);
                for (auto& w : ws) {
                    if (std::filesystem::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
                    w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
                }
                for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : ws) if (w.name == t.first) w.tops.push_back(t.second);
                cube = readSegy(dir + "/seismic.segy", e);
                return [=](GeoState& S) {
                    S.dir = dir; S.hz = hz; S.faults = fl; S.wells = ws; S.cube = cube; ++S.cubeRev; S.hasCube = cube.ni > 0; S.cubePath = dir + "/seismic.segy"; S.projectPath.clear(); S.folders = {}; S.segy.reset(); S.stick.clear(); S.draft = Fault();
                    S.hasGrid = S.hasFacies = S.hasPetro = S.hasVol = S.hasUp = S.hasTh = false; S.blocked.clear(); S.seeds.clear();
                    const double L = fieldU ? 3.28084 : 1.0;
                    S.vm.v0 = 1800 * L; S.vm.k = 0.6; S.vm.datum = 0;
                    S.gs = GridSpec(); S.gs.dx = S.gs.dy = 100 * L; S.gs.nx = 40; S.gs.ny = 30;
                    Zone z1; z1.name = "Upper"; z1.top = 0; z1.base = 1; z1.nLayers = 8;
                    Zone z2; z2.name = "Lower"; z2.top = 1; z2.base = 2; z2.nLayers = 8;
                    S.gs.zones = {z1, z2}; S.gs.tieRange = 1500 * L;
                    S.vs.contacts.owc = 2080 * L; S.vs.field = fieldU;
                    S.fs.vg.rangeMax = S.fs.vg.rangeMin = 800 * L; S.fs.vg.rangeVert = 6 * L;
                    S.ps.vg.rangeMax = S.ps.vg.rangeMin = 1000 * L; S.ps.vg.rangeVert = 8 * L;
                    S.msg("Demo project in " + dir + ": " + std::to_string(files.size()) + " files (" + (e.empty() ? std::string("ok") : e) + ")");
                };
            });
        }
        ui::hint("Synthetic faulted anticline: 3 IRAP horizons, fault sticks, 6 wells (LAS, tops, one horizontal survey), a SEG-Y volume");
        geoProjectBar(S);   // C46 project file
        ImGui::Checkbox(TR("FIELD units (ft) for new projects"), &S.field);
        ImGui::Separator();
        std::string p;
        // C42: every importer takes several files; tables go through the column-assignment dialog
        std::vector<std::string> files;
        if (ImGui::Button(TR(ICON_IMAGE " Horizons (IRAP / ZMAP)")) && plat::fileDialogMulti("Surfaces (*.irap;*.grd;*.zmap;*.dat)\0*.irap;*.grd;*.zmap;*.dat;*.txt\0All files (*.*)\0*.*\0", files)) {
            for (const auto& p : files) {
                std::string e; Surface s = readIrap(p, e);
                if (!s.valid()) { e.clear(); s = readZmap(p, e); }
                if (s.valid()) { s.name = baseName(p); S.hz.push_back(s); ++S.hzRev; S.msg("horizon " + s.name + ": " + std::to_string(s.nx) + " x " + std::to_string(s.ny)); }
                else S.msg("Error: " + e);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_IMAGE " Points (XYZ)..."))) geoOpenImport(S, ImportKind::Points);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Scattered points (one or several sets per file) -> kriged surfaces; columns assigned in a dialog");
        if (ImGui::Button(TR(ICON_PIN " Well heads..."))) geoOpenImport(S, ImportKind::WellHeads);
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_ROUTE " Deviation..."))) geoOpenImport(S, ImportKind::Deviation);
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_LOCATE_FIXED " Tops..."))) geoOpenImport(S, ImportKind::Tops);
        if (ImGui::Button(TR(ICON_CHART_LINE " LAS logs")) && plat::fileDialogMulti("LAS (*.las)\0*.las;*.LAS\0All files (*.*)\0*.*\0", files)) geoImportLas(S, files);
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_SPLIT " Fault sticks..."))) geoOpenImport(S, ImportKind::FaultSticks);
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_TABLE " Well manager"))) { S.we.open = true; }
        ui::hint("Select several files at once. Text tables open the column-assignment dialog: pick the role of each column (well, X, Y, KB, MD, INC, AZI, TVDSS ...), "
                 "units and the file format; the choice is remembered for files with the same columns. The well manager adds wells by hand and edits each well's survey and tops.");
        auto openOutOfCore = [&](const std::string& p) {
            const SegyByteMap bm = S.segyMap;
            launch("Indexing SEG-Y (out-of-core)", [p, bm]() -> std::function<void(GeoState&)> {
                auto v = std::make_shared<SegyVolume>();
                v->byteIl = bm.inlineByte; v->byteXl = bm.xlineByte; v->byteX = bm.xByte; v->byteY = bm.yByte; v->byteScalar = bm.scalarByte;
                std::string e;
                if (!v->open(p, e)) return [e](GeoState& S) { S.msg("Error: " + e); };
                SeismicCube proxy = v->proxy(160, 512);
                return [v, proxy, p](GeoState& S) {
                    S.segy = v; S.cube = proxy; S.cubePath = p; ++S.cubeRev; S.hasCube = true; S.hasTh = false; S.seeds.clear(); S.il = v->ni / 2; S.xl = v->nj / 2; S.tSlice = -1;
                    char b[400];
                    std::snprintf(b, sizeof b, "SEG-Y %s (out-of-core): %d x %d traces x %d samples, %.1f MB, format %d%s, %lld missing, geometry fit %.3g; proxy %d x %d x %d in memory for tracking / volume",
                                  v->name.c_str(), v->ni, v->nj, v->ns, v->fileBytes / 1048576.0, v->format, v->littleEndian ? " LE" : "", v->missing, v->fitError, proxy.ni, proxy.nj, proxy.ns);
                    S.msg(b);
                };
            });
        };
        bool go3d = false;
        if (ImGui::Button(TR(ICON_ACTIVITY " SEG-Y (2-D lines / 3-D)")) && plat::fileDialogMulti("SEG-Y (*.sgy;*.segy)\0*.sgy;*.segy;*.SGY;*.SEGY\0All files (*.*)\0*.*\0", files)) {
            // C42: several files; each is detected as a 2-D line or a 3-D volume from the trace headers
            std::vector<std::string> l2, l3;
            for (const auto& f : files) {
                std::string info, e; const int k = segyGeometryKind(f, S.segyMap, info, e);
                if (k == 2) l2.push_back(f); else if (k == 3) l3.push_back(f); else S.msg("Error: " + e);
                if (k) S.msg(baseName(f) + ": " + (k == 2 ? "2-D line" : "3-D volume") + " - " + info);
            }
            if (!l2.empty()) {
                const SegyByteMap bm = S.segyMap;
                launch("Reading " + std::to_string(l2.size()) + " 2-D line(s)", [l2, bm]() -> std::function<void(GeoState&)> {
                    std::vector<SeismicLine> out; std::vector<std::string> errs;
                    for (const auto& f : l2) { std::string e; SeismicLine L = readSegy2D(f, bm, e); if (!e.empty()) errs.push_back(e); if (L.ntr() > 0) out.push_back(std::move(L)); }
                    return [out, errs](GeoState& S) mutable {
                        for (auto& e : errs) S.msg("Warning: " + e);
                        const size_t n = out.size(); geoAddLines(S, std::move(out));
                        S.msg(std::to_string(n) + " 2-D line(s) loaded - Seismic step / Seismic section view / 3-D scene");
                    };
                });
            }
            if (l3.size() > 1) S.msg("Note: " + std::to_string(l3.size()) + " 3-D volumes selected - only " + baseName(l3[0]) + " is loaded (one 3-D volume per project)");
            if (!l3.empty()) { p = l3[0]; go3d = true; }   // C43: queued behind the 2-D lines when both are selected
        }
        if (go3d) {
            std::error_code ec; const auto bytes = std::filesystem::file_size(p, ec);
            if (!ec && bytes > 400ull * 1048576ull) { S.msg("SEG-Y larger than 400 MB - opened out-of-core (full-resolution sections from the file, proxy in memory)"); openOutOfCore(p); }
            else launch("Reading SEG-Y", [p, bm = S.segyMap]() -> std::function<void(GeoState&)> {
                std::string e; SeismicCube c = readSegy(p, bm, e);
                return [c, e, p](GeoState& S) { if (!e.empty()) { S.msg("Error: " + e); return; } S.cube = c; S.cubePath = p; ++S.cubeRev; S.hasCube = true; S.segy.reset(); S.hasTh = false; S.seeds.clear(); S.msg("SEG-Y " + c.name + ": " + std::to_string(c.ni) + " x " + std::to_string(c.nj) + " traces, " + std::to_string(c.ns) + " samples"); };
            });
        }
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_ACTIVITY " SEG-Y (large, out-of-core)")) && fileDialog(false, "SEG-Y (*.sgy;*.segy)\0*.sgy;*.segy;*.SGY;*.SEGY\0All files (*.*)\0*.*\0", "sgy", p)) openOutOfCore(p);
        geoSegyBytesUi(S);
        ui::hint("Large volumes: only a trace index is built; sections and slices are read from disk on demand at full resolution (LRU trace cache). "
                 "Tracking, attributes cubes and the volume rendering use a decimated in-memory proxy");
        ImGui::Separator();
        if (ui::beginTable("##wells", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders, ImVec2(0, fs * 10))) {
            ui::tableColumn("well"); ui::tableColumn("x"); ui::tableColumn("y"); ui::tableColumn("logs / tops"); ui::tableColumn("");
            ImGui::TableHeadersRow();
            for (size_t i = 0; i < S.wells.size(); ++i) {
                auto& w = S.wells[i];
                ImGui::TableNextRow(); ImGui::PushID(int(i));
                ImGui::TableNextColumn(); ImGui::TextUnformatted(w.name.c_str());
                ImGui::TableNextColumn(); ImGui::Text("%.0f", w.x);
                ImGui::TableNextColumn(); ImGui::Text("%.0f", w.y);
                ImGui::TableNextColumn(); ImGui::Text("%zu / %zu%s", w.logs.size(), w.tops.size(), w.traj.empty() ? "" : " dev");
                ImGui::TableNextColumn();
                if (ImGui::SmallButton(ICON_PENCIL)) { S.we.open = true; S.we.sel = int(i); }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("edit the head, survey and tops (well manager)");
                ImGui::SameLine(); if (ImGui::SmallButton(ICON_X)) { S.wells.erase(S.wells.begin() + long(i)); ImGui::PopID(); break; }
                ImGui::PopID();
            }
            ui::endTable();
        }
        for (size_t i = 0; i < S.hz.size(); ++i) {
            ImGui::PushID(int(i) + 1000);
            ImGui::BulletText("%s  %dx%d  (%.0f x %.0f)", S.hz[i].name.c_str(), S.hz[i].nx, S.hz[i].ny, S.hz[i].dx, S.hz[i].dy);
            ImGui::SameLine(); if (ImGui::SmallButton(ICON_X)) { S.hz.erase(S.hz.begin() + long(i)); ImGui::PopID(); break; }
            ImGui::PopID();
        }
        geoPolygonsPanel(S);   // C43
    } else if (S.step == 1) {   // ---------------- seismic
        if (!S.lines.empty()) { geoLinesPanel(S); geoMistiePanel(S); }
        if (!S.hasCube) { if (S.lines.empty()) ImGui::TextDisabled("Import SEG-Y files (2-D lines or a 3-D volume) in the Data step."); }
        else {
            if (!S.lines.empty()) ImGui::SeparatorText(TR(ICON_BOX " 3-D volume"));
            ImGui::Text("%s: %d inlines x %d crosslines x %d samples (dt %.1f ms)", S.cube.name.c_str(), S.cube.ni, S.cube.nj, S.cube.ns, S.cube.dt);
            if (ui::beginProps("##vel")) {
                ui::propD("velocity v0", &S.vm.v0, "%.1f", "interval velocity at the datum (length unit / s)");
                ui::propD("gradient k", &S.vm.k, "%.3f", "v(z) = v0 + k (z - datum), 1/s");
                ui::propD("datum depth", &S.vm.datum, "%.1f");
                static const char* modes[] = {"peak", "trough", "zero crossing + to -"};
                ui::propCombo("track event", &S.trackMode, modes, 3);
                ui::propD("search window (ms)", &S.window, "%.1f");
                ui::propD("max jump between traces (ms)", &S.maxJump, "%.1f");
                ui::propString("horizon name", &S.thName);
                ui::endProps();
            }
            ImGui::Text("%zu seed(s) - Pick: tracking seeds, then click on the seismic section or the 3-D slices; one per fault block", S.seeds.size());
            if (ImGui::Button(TR(ICON_X " Clear seeds"))) S.seeds.clear();
            ImGui::SameLine();
            if (!S.hz.empty() && ImGui::Button(TR("Seed from the map at the wells"))) {
                // a seed at every well head, at the time of the horizon named like the target (or the first)
                const Surface* h = &S.hz[0]; for (auto& x : S.hz) if (x.name == S.thName) h = &x;
                for (const auto& w : S.wells) { double i, j; if (S.cube.ij(w.x, w.y, i, j) && std::isfinite(h->at(w.x, w.y))) S.seeds.push_back({i, j, depthToTime(S.vm, h->at(w.x, w.y))}); }
            }
            if (ui::accentButton(ICON_PLAY " Auto-track") && !S.seeds.empty()) {
                const SeismicCube c = S.cube; const auto seeds = S.seeds; const int m = S.trackMode; const double w = S.window, mj = S.maxJump;
                launch("Auto-tracking", [c, seeds, m, w, mj]() -> std::function<void(GeoState&)> {
                    TimeHorizon th = autoTrack(c, seeds, m, w, mj);
                    int n = 0; for (double t : th.t) n += std::isfinite(t);
                    return [th, n](GeoState& S) { S.th = th; S.hasTh = true; S.msg("tracked " + std::to_string(n) + " / " + std::to_string(th.ni * th.nj) + " traces"); };
                });
            }
            if (S.hasTh && ImGui::Button(TR(ICON_ARROW_DOWN_TO_LINE " Depth-convert -> horizon"))) {
                Surface d = horizonToDepth(S.cube, S.th, S.vm, S.thName);
                bool rep = false; for (auto& h : S.hz) if (h.name == S.thName) { h = d; rep = true; }
                if (!rep) S.hz.push_back(d);
                S.msg("horizon " + S.thName + (rep ? " replaced" : " added") + " from the tracked time horizon (undefined traces are in-painted by the grid builder)");
            }
            ui::hint("Velocity model: linear v0 + k z; time-depth conversion t(z) = (2000 / k) ln(1 + k (z - datum) / v0) ms TWT");
            ImGui::Separator();
            ImGui::TextDisabled("Attributes (sections, 3-D slices, volume)");
            if (ui::beginProps("##attr")) {
                static const char* an[] = {"Amplitude", "Envelope", "Inst. phase", "Inst. frequency", "Cosine of phase", "Sweetness", "RMS amplitude", "Relative impedance", "Coherence (semblance)", "Variance"};
                ui::propCombo("attribute", &S.secAttr, an, int(SeisAttr::Count), "Hilbert-transform attributes (envelope, phase, frequency from the phase difference), RMS / relative impedance in a window, 3x3 semblance coherence and variance");
                ui::propI("window (samples)", &S.attrWin, 1, 64, "half window of the RMS / coherence / variance attributes");
                ui::propSliderD("clip percentile", &S.clipPct, 0.8, 1.0, "%.3f", "colour range = this percentile of |value|");
                ui::endProps();
            }
            if (S.secAttr != 0 && ImGui::Button(TR(ICON_SAVE " Write the attribute volume (SEG-Y)..."))) {
                std::string p = "attribute.sgy";
                if (fileDialog(true, "SEG-Y (*.sgy)\0*.sgy;*.segy\0", "sgy", p)) {
                    const auto vol = S.segy; const SeismicCube c = S.segy ? SeismicCube() : S.cube; const SeisAttr a = SeisAttr(S.secAttr); const int w = S.attrWin;
                    launch("Writing the attribute volume", [vol, c, a, w, p]() -> std::function<void(GeoState&)> {
                        std::string e; bool ok;
                        if (vol) ok = volumeAttribute(*vol, a, w, p, e);   // streamed inline by inline at full resolution
                        else { ok = writeSegy(p, cubeAttribute(c, a, w)); if (!ok) e = "cannot write " + p; }
                        return [ok, e, p, a](GeoState& S) { S.msg(ok ? std::string(seisAttrName(a)) + " volume written: " + p : "Error: " + e); };
                    });
                }
            }
        }
        geoSeisCalcPanel(S);    // C43
        geoLayerCakePanel(S);   // C43
    } else if (S.step == 2) {   // ---------------- structure
        ImGui::TextDisabled("Horizons, top to bottom (order defines the zones)");
        for (size_t i = 0; i < S.hz.size(); ++i) {
            ImGui::PushID(int(i));
            double lo, hi; S.hz[i].range(lo, hi);
            ImGui::Text("%zu  %s   %.1f .. %.1f", i, S.hz[i].name.c_str(), lo, hi);
            ImGui::SameLine(); if (i > 0 && ImGui::SmallButton(TR("up"))) std::swap(S.hz[i], S.hz[i - 1]);
            ImGui::PopID();
        }
        ImGui::Separator();
        ImGui::TextDisabled("Faults (sticks -> stair-stepped grid faults)");
        for (size_t f = 0; f < S.faults.size(); ++f) {
            ImGui::PushID(int(f) + 500);
            ImGui::Text("%s: %zu sticks", S.faults[f].name.c_str(), S.faults[f].sticks.size());
            if (ui::beginProps("##f")) {
                ui::propD("MULTFLT", &S.faults[f].mult, "%.3g", "transmissibility multiplier of the generated fault faces");
                ui::propD("throw hint", &S.faults[f].throwHint, "%.1f", "added to the horizons on the hanging-wall side when the maps are smooth across the fault (0: the maps carry the throw)");
                ui::endProps();
            }
            if (ImGui::SmallButton(TR(ICON_X " remove"))) { S.faults.erase(S.faults.begin() + long(f)); ImGui::PopID(); break; }
            ImGui::PopID();
        }
        // ---- C41 fault network: junctions, truncation rules, throw QC
        ImGui::Separator();
        if (ImGui::TreeNodeEx("Fault network (junctions, truncation, throw QC)", S.faults.size() >= 2 ? ImGuiTreeNodeFlags_DefaultOpen : 0)) {
            const double tol = S.fnTol > 0 ? S.fnTol : 2.0 * std::max(S.gs.dx, S.gs.dy);
            ImGui::SetNextItemWidth(fs * 7);
            ImGui::InputDouble("snap tolerance (0 = 2 cells)", &S.fnTol, 0, 0, "%.4g");
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("A minor fault whose tip ends within this plan distance of another fault is snapped onto it (Y junction).");
            ImGui::BeginDisabled(S.faults.size() < 2);
            if (ImGui::Button(TR(ICON_SEARCH " Detect junctions"))) {
                S.fnJ = detectJunctions(S.faults, S.zMid(), tol);
                int nx_ = 0; for (const auto& j : S.fnJ) nx_ += j.type == 'X';
                S.fnMsg = std::to_string(S.fnJ.size()) + " junction(s) at " + std::to_string(int(S.zMid())) + " depth (" + std::to_string(nx_) + " crossing)";
            }
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::BeginDisabled(S.fnJ.empty());
            if (ImGui::Button(TR(ICON_SCISSORS " Apply truncation"))) {
                FaultNetReport rep;
                if (S.fnBackup.empty()) S.fnBackup = S.faults;
                S.faults = truncateFaults(S.faults, S.fnJ, tol, &rep);
                S.fnMsg = std::to_string(rep.truncated) + " fault(s) truncated, " + std::to_string(rep.sticksRemoved) + " stick(s) removed";
                for (const auto& n : rep.notes) S.msg("Fault network: " + n);
                S.fnJ = detectJunctions(S.faults, S.zMid(), tol);
                if (S.hasGrid) S.msg("Fault network changed: rebuild the 3-D grid to use it");
            }
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::BeginDisabled(S.fnBackup.empty());
            if (ImGui::Button(TR(ICON_UNDO_2 " Revert"))) { S.faults = S.fnBackup; S.fnBackup.clear(); S.fnJ.clear(); S.fnMsg = "faults restored"; }
            ImGui::EndDisabled();
            if (!S.fnMsg.empty()) ImGui::TextDisabled("%s", S.fnMsg.c_str());
            if (!S.fnJ.empty() && ImGui::BeginTable("##fnj", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupColumn("type"); ImGui::TableSetupColumn("major"); ImGui::TableSetupColumn("minor");
                ImGui::TableSetupColumn("gap"); ImGui::TableSetupColumn("keep side"); ImGui::TableSetupColumn("truncate");
                ImGui::TableHeadersRow();
                for (size_t q = 0; q < S.fnJ.size(); ++q) {
                    FaultJunction& J = S.fnJ[q];
                    if (J.major < 0 || J.minor < 0 || size_t(J.major) >= S.faults.size() || size_t(J.minor) >= S.faults.size()) continue;
                    ImGui::PushID(int(q));
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%c", J.type);
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip(J.type == 'X' ? "Crossing: both faults continue. Tick 'truncate' to cut the minor one at the major." : "Y / lambda: the minor fault ends at the major.");
                    ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(S.faults[size_t(J.major)].name.c_str());
                    ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(S.faults[size_t(J.minor)].name.c_str());
                    ImGui::SameLine();
                    if (ImGui::SmallButton(ICON_SHUFFLE)) std::swap(J.major, J.minor);
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Swap major and minor");
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%.3g", J.gap);
                    ImGui::TableSetColumnIndex(4);
                    static const char* kS[] = {"auto", "left", "right"};
                    int ks = J.keepSide == 0 ? 0 : J.keepSide > 0 ? 1 : 2;
                    ImGui::SetNextItemWidth(fs * 4.5f);
                    if (ImGui::Combo("##ks", &ks, kS, 3)) J.keepSide = ks == 0 ? 0 : ks == 1 ? 1 : -1;
                    ImGui::TableSetColumnIndex(5); ImGui::Checkbox("##tr", &J.truncate);
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            // horizon-fault throw QC
            if (!S.faults.empty() && !S.hz.empty()) {
                S.fnThrowF = std::clamp(S.fnThrowF, 0, int(S.faults.size()) - 1); S.fnThrowH = std::clamp(S.fnThrowH, 0, int(S.hz.size()) - 1);
                std::vector<std::string> fn, hn; for (const auto& f : S.faults) fn.push_back(f.name); for (const auto& h : S.hz) hn.push_back(h.name);
                ImGui::SetNextItemWidth(fs * 7);
                if (ImGui::BeginCombo("##tf", fn[size_t(S.fnThrowF)].c_str())) { for (int q = 0; q < int(fn.size()); ++q) if (ImGui::Selectable(fn[size_t(q)].c_str(), q == S.fnThrowF)) { S.fnThrowF = q; S.fnThrowOk = false; } ImGui::EndCombo(); }
                ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
                if (ImGui::BeginCombo("##th", hn[size_t(S.fnThrowH)].c_str())) { for (int q = 0; q < int(hn.size()); ++q) if (ImGui::Selectable(hn[size_t(q)].c_str(), q == S.fnThrowH)) { S.fnThrowH = q; S.fnThrowOk = false; } ImGui::EndCombo(); }
                ImGui::SameLine();
                if (ImGui::Button(TR("Throw profile"))) {
                    S.fnThrow = faultThrowProfile(S.faults[size_t(S.fnThrowF)], S.hz[size_t(S.fnThrowH)], 1.5 * std::max(S.gs.dx, S.gs.dy), 31);
                    S.fnThrowOk = true;
                }
                if (S.fnThrowOk) {
                    const ThrowProfile& P = S.fnThrow;
                    if (P.s.empty()) ImGui::TextDisabled("The horizon is not defined on both sides of the fault.");
                    else {
                        ImGui::TextDisabled("max |throw| %.4g%s%s", P.maxThrow, P.reverses ? "  - sign reversal (scissor fault or mis-picked side)" : "",
                                            P.tipsTaper ? "" : "  - throw does not taper at the tips (fault should end at another fault or the model edge)");
                        if (ImPlot::BeginPlot("##throw", ImVec2(-1, fs * 9))) {
                            ImPlot::SetupAxes("distance along the trace", "throw (HW - FW)");
                            ImPlot::PlotLine("throw", P.s.data(), P.throwV.data(), int(P.s.size()));
                            ImPlot::EndPlot();
                        }
                    }
                }
            }
            ImGui::TreePop();
        }
        // ---- C36-4 horizon editing
        ImGui::Separator();
        ImGui::TextDisabled("Horizon editing");
        if (S.hz.empty()) ImGui::TextDisabled("Import or track horizons first.");
        else {
            S.editHz = std::clamp(S.editHz, 0, int(S.hz.size()) - 1);
            std::vector<std::string> hn; for (const auto& h : S.hz) hn.push_back(h.name);
            static const char* bm[] = {"raise / lower", "smooth", "flatten"};
            if (ui::beginProps("##hed")) {
                if (ui::propCombo("horizon", &S.editHz, hn)) S.flatDepth = 0;
                if (ui::propCombo("brush", &S.brushMode, bm, 3)) S.flatDepth = 0;
                ui::propD("brush radius", &S.brushR, "%.0f", "Gaussian taper exp(-3 r^2 / R^2) inside R");
                if (S.brushMode == 0) ui::propD("dz per click", &S.brushAmt, "%+.2f", "> 0 deepens, < 0 raises the horizon at the brush centre");
                else if (S.brushMode == 1) ui::propD("strength (0..1)", &S.brushStrength, "%.2f");
                else ui::propD("flatten to depth", &S.flatDepth, "%.2f", "0: the depth under the first click");
                ui::propI("smoothing passes", &S.smoothPasses, 1, 50);
                ui::propD("minimum zone thickness", &S.minThick, "%.2f", "horizon h is pushed to at least h-1 + this (top-down)");
                ui::propD("isochore thickness", &S.isoThick, "%.2f", "new conformable horizon below the selected one");
                ui::endProps();
                S.brushR = std::max(S.brushR, 1e-3); S.brushStrength = std::clamp(S.brushStrength, 0.0, 1.0);
            }
            bool on = S.pickMode == 3;
            if (ImGui::Checkbox(TR("Brush: click on the map / 3-D scene"), &on)) { S.pickMode = on ? 3 : 0; if (on) { if (S.view != 4) S.view = 0; S.mapSurface = S.editHz; } }
            Surface& H = S.hz[size_t(S.editHz)];
            auto edited = [&](const std::string& m) { ++S.hzRev; S.msg(m + (S.hasGrid ? " - rebuild the grid to apply" : "")); };
            if (ImGui::Button(TR("Smooth"))) { S.pushUndo(); smoothSurface(H, S.smoothPasses); edited(H.name + ": smoothed (" + std::to_string(S.smoothPasses) + " passes)"); }
            ImGui::SameLine();
            if (ImGui::Button(TR("Fill holes"))) { S.pushUndo(); const int n = fillSurfaceHoles(H); if (n) edited(H.name + ": " + std::to_string(n) + " undefined nodes in-painted"); else { S.hzUndo.pop_back(); S.msg(H.name + ": no holes"); } }
            ImGui::SameLine();
            if (ImGui::Button(TR("Min. thickness"))) { S.pushUndo(); const int n = enforceThickness(S.hz, S.minThick); if (n) edited("minimum thickness: " + std::to_string(n) + " nodes moved"); else { S.hzUndo.pop_back(); S.msg("no crossing / thin zones"); } }
            ui::hint("Fill holes: Laplace in-painting of undefined nodes. Min. thickness: removes crossing horizons top-down (all horizons, same lattice)");
            if (ImGui::Button(TR(ICON_PLUS " Conformable horizon below"))) {
                S.pushUndo();
                Surface n = isochoreSurface(H, S.isoThick, H.name + "_iso" + std::to_string(int(std::lround(S.isoThick))));
                S.hz.insert(S.hz.begin() + S.editHz + 1, n); S.gs.zones.clear();
                edited("new horizon " + n.name + " (" + std::to_string(S.isoThick) + " below " + H.name + "); zones reset");
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(S.hzUndo.empty());
            if (ImGui::Button(TR(ICON_UNDO_2 " Undo"))) { S.hz = S.hzUndo.back(); S.hzUndo.pop_back(); ++S.hzRev; S.msg("horizon edit undone"); }
            ImGui::EndDisabled();
            ImGui::SameLine(); ImGui::TextDisabled("%zu", S.hzUndo.size());
            if (ImGui::Button(TR(ICON_SAVE " Save horizon (IRAP)..."))) {
                std::string p = H.name + ".irap";
                if (fileDialog(true, "IRAP classic (*.irap)\0*.irap;*.gri\0", "irap", p)) S.msg(writeIrap(p, H) ? "written " + p : "Error: cannot write " + p);
            }
        }
    } else if (S.step == 3) {   // ---------------- grid
        if (ui::beginProps("##gs")) {
            ui::propD("origin x", &S.gs.x0, "%.1f"); ui::propD("origin y", &S.gs.y0, "%.1f"); ui::propD("rotation (deg)", &S.gs.rotDeg, "%.2f");
            ui::propD("dx", &S.gs.dx, "%.2f"); ui::propD("dy", &S.gs.dy, "%.2f");
            ui::propI("nx", &S.gs.nx, 1, 2000); ui::propI("ny", &S.gs.ny, 1, 2000);
            ui::propD("min. cell thickness", &S.gs.minThickness, "%.3f", "thinner cells are inactive (pinch-out)");
            ui::propBool("tie horizons to the well tops", &S.gs.tieToWells);
            ui::propD("tie correction range", &S.gs.tieRange, "%.0f", "Gaussian variogram range of the kriged residual (0 = 10 cells)");
            static const char* pl[] = {"vertical", "fault-parallel (sloped near faults)"};
            ui::propCombo("pillars", &S.gs.pillars, pl, 2, "fault-parallel: the COORD lines follow the dip of the nearest fault plane at the fault and become vertical "
                                                          "over the taper distance - fault faces are planar instead of stair-stepped in depth");
            if (S.gs.pillars == 1) ui::propD("pillar taper distance", &S.gs.pillarRange, "%.0f", "Gaussian taper of the pillar slope away from the fault (0 = 4 cells)");
            ui::endProps();
        }
        if (!S.hz.empty() && ImGui::Button(TR("Fit the lattice to the first horizon"))) {
            const Surface& h = S.hz[0];
            S.gs.x0 = h.x0; S.gs.y0 = h.y0; S.gs.rotDeg = h.rotDeg;
            S.gs.nx = std::max(1, int(std::floor((h.nx - 1) * h.dx / S.gs.dx))); S.gs.ny = std::max(1, int(std::floor((h.ny - 1) * h.dy / S.gs.dy)));
        }
        ImGui::TextDisabled("Zones");
        if (S.gs.zones.empty() && S.hz.size() >= 2) for (size_t h = 0; h + 1 < S.hz.size(); ++h) { Zone z; z.name = S.hz[h].name; z.top = int(h); z.base = int(h + 1); S.gs.zones.push_back(z); }
        static const char* lm[] = {"proportional", "parallel to top", "parallel to base"};
        for (size_t z = 0; z < S.gs.zones.size(); ++z) {
            ImGui::PushID(int(z) + 200);
            Zone& Z = S.gs.zones[z];
            ImGui::SetNextItemWidth(fs * 6); ui::inputString("##zn", &Z.name);
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4); ImGui::InputInt("##nl", &Z.nLayers, 0); Z.nLayers = std::clamp(Z.nLayers, 1, 500);
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); ImGui::Combo("##lm", &Z.mode, lm, 3);
            ImGui::SameLine(); ImGui::Text("%d-%d", Z.top, Z.base);
            ImGui::PopID();
        }
        if (ui::accentButton(ICON_PLAY " Build the 3-D grid") && S.hz.size() >= 2) {
            const GridSpec gs = S.gs; const auto hz = S.hz; const auto fl = S.faults; const auto ws = S.wells;
            launch("Building the structural grid", [gs, hz, fl, ws]() -> std::function<void(GeoState&)> {
                std::string e;
                GeoGrid g = buildGrid(gs, hz, fl, ws, e);
                return [g, e](GeoState& S) {
                    if (!e.empty()) { S.msg("Error: " + e); return; }
                    S.g = g; S.hasGrid = true; S.hasFacies = S.hasPetro = S.hasVol = S.hasUp = false; S.blocked.clear();
                    size_t nf = 0; for (auto& f : g.grid.faults) nf += f.faces.size();
                    int act = 0; for (int a : g.grid.actnum) act += a != 0;
                    S.msg("grid " + std::to_string(g.grid.nx) + " x " + std::to_string(g.grid.ny) + " x " + std::to_string(g.grid.nz) + ", " + std::to_string(act) + " active, " +
                          std::to_string(g.nBlocks) + " fault blocks, " + std::to_string(nf) + " fault face records");
                    for (const auto& l : g.tieReport) S.msg("  tie " + l);
                    if (g.spec.pillars == 1) { char b[96]; std::snprintf(b, sizeof b, "  fault-parallel pillars: max slope %.3f (dx/dz)", g.maxPillarSlope()); S.msg(b); }
                };
            });
        }
    } else if (S.step == 4) {   // ---------------- blocking
        if (!S.hasGrid) ImGui::TextDisabled("Build the grid first.");
        else {
            ImGui::TextWrapped("Logs are averaged over the trajectory length inside each cell: PERM* geometric, FACIES / LITH* most-of, others arithmetic.");
            if (ui::accentButton(ICON_PLAY " Block the wells")) {
                S.blocked = blockWells(S.g, S.wells);
                std::map<std::string, int> per; for (auto& b : S.blocked) per[b.well]++;
                std::string m = "blocked " + std::to_string(S.blocked.size()) + " cells:";
                for (auto& kv : per) m += " " + kv.first + " " + std::to_string(kv.second);
                S.msg(m);
            }
            ImGui::Text("%zu blocked cells", S.blocked.size());
        }
    } else if (S.step == 5) {   // ---------------- facies
        static const char* fm[] = {"Sequential indicator simulation", "Truncated Gaussian simulation", "Object-based channels", "Multiple-point statistics (Direct Sampling)"};
        if (ui::beginProps("##fs")) {
            ui::propCombo("method", &S.fs.method, fm, 4);
            ui::propString("facies log", &S.fs.curve);
            if (S.fs.method != 3) {
                ui::propD("range major", &S.fs.vg.rangeMax, "%.0f"); ui::propD("range minor", &S.fs.vg.rangeMin, "%.0f"); ui::propD("range vertical", &S.fs.vg.rangeVert, "%.1f");
                ui::propD("azimuth (deg)", &S.fs.vg.azimuthDeg, "%.1f");
            }
            if (S.fs.method == 2) {
                ui::propD("channel width", &S.fs.chWidth, "%.0f"); ui::propD("channel thickness", &S.fs.chThick, "%.1f");
                ui::propD("flow azimuth (deg)", &S.fs.chAzimuth, "%.0f"); ui::propD("sinuosity amplitude", &S.fs.chAmplitude, "%.0f"); ui::propD("wavelength", &S.fs.chWavelength, "%.0f");
            }
            int sd = int(S.fs.seed); ui::propI("seed", &sd, 1, 1 << 30); S.fs.seed = unsigned(sd);
            ui::endProps();
        }
        if (S.fs.method == 3) {   // ---- C36-4 training image + Direct Sampling parameters
            ImGui::TextDisabled("Training image %s", S.ti ? ("(" + S.tiName + ", " + std::to_string(S.ti->nx) + " x " + std::to_string(S.ti->ny) + " x " + std::to_string(S.ti->nz) + ")").c_str() : "(none)");
            if (ui::beginProps("##ti")) {
                ui::propI("TI size x (px)", &S.tiNx, 20, 2000); ui::propI("TI size y (px)", &S.tiNy, 20, 2000);
                ui::propD("channel width (px)", &S.tiWidth, "%.1f"); ui::propD("meander wavelength (px)", &S.tiWave, "%.0f"); ui::propD("meander amplitude (px)", &S.tiAmp, "%.1f");
                ui::propD("channel fraction", &S.tiFrac, "%.2f");
                int sd = int(S.tiSeed); ui::propI("TI seed", &sd, 1, 1 << 30); S.tiSeed = unsigned(sd);
                ui::endProps();
            }
            if (ImGui::Button(TR(ICON_WAND_SPARKLES " Generate channel TI"))) {
                S.ti = std::make_shared<TrainingImage>(channelTrainingImage(S.tiNx, S.tiNy, 1, S.tiWidth, S.tiWave, S.tiAmp, std::clamp(S.tiFrac, 0.02, 0.9), S.tiSeed));
                S.tiName = "channels"; S.msg("training image: sinuous channels " + std::to_string(S.tiNx) + " x " + std::to_string(S.tiNy));
            }
            ImGui::SameLine();
            { std::string p;
              if (ImGui::Button(TR(ICON_FOLDER_OPEN " Load (GSLIB)...")) && fileDialog(false, "GSLIB / SGeMS training image (*.gslib;*.dat;*.sgems)\0*.gslib;*.dat;*.sgems;*.txt\0All\0*.*\0", "gslib", p)) {
                  std::string e; TrainingImage t = readTrainingImage(p, e);
                  if (!e.empty() || t.v.empty()) S.msg("Error: " + (e.empty() ? std::string("empty training image") : e));
                  else { S.ti = std::make_shared<TrainingImage>(t); S.tiName = baseName(p); S.msg("training image " + p + ": " + std::to_string(t.ncat()) + " categories"); }
              } }
            ImGui::SameLine();
            ImGui::BeginDisabled(!S.ti);
            { std::string p = "training_image.gslib";
              if (ImGui::Button(TR(ICON_SAVE " Save...")) && fileDialog(true, "GSLIB (*.gslib)\0*.gslib\0", "gslib", p)) S.msg(writeTrainingImage(p, *S.ti) ? "written " + p : "Error: cannot write " + p); }
            ImGui::EndDisabled();
            if (ui::beginProps("##mps")) {
                ui::propI("conditioning neighbours", &S.mps.maxNbr, 4, 200);
                ui::propD("distance threshold", &S.mps.threshold, "%.3f", "accepted fraction of mismatching neighbours (0 = exact pattern)");
                ui::propD("max. scan fraction", &S.mps.maxScan, "%.2f", "fraction of the training image scanned per node");
                ui::propI("max. scanned positions", &S.mps.maxScanAbs, 100, 10000000);
                ui::propD("TI azimuth (deg)", &S.mps.azimuthDeg, "%.1f", "map direction of the TI x axis (clockwise from north)");
                ui::propD("TI pixel size", &S.mps.tiCell, "%.1f", "plan size of a TI pixel (0 = grid dx)");
                ui::propD("search radius (TI px)", &S.mps.radius, "%.1f");
                ui::propD("TI layers per grid layer", &S.mps.tiLayer, "%.2f", "3-D training images only");
                ui::propD("servosystem", &S.mps.servo, "%.2f", "0 = plain Direct Sampling; > 0 steers the facies proportions to the TI marginal");
                int sd = int(S.mps.seed); ui::propI("seed##mps", &sd, 1, 1 << 30); S.mps.seed = unsigned(sd);
                ui::endProps();
            }
            if (S.ti && ImPlot::BeginPlot("##tiprev", ImVec2(-1, fs * 12), ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxes("TI x (px)", "TI y (px)");
                const TrainingImage& T = *S.ti;
                std::vector<float> img(size_t(T.nx) * size_t(T.ny));
                for (int j = 0; j < T.ny; ++j) for (int i = 0; i < T.nx; ++i) img[size_t((T.ny - 1 - j) * T.nx + i)] = float(T.at(i, j, 0));
                ImPlot::PushColormap(ImPlotColormap_Viridis);
                ImPlot::PlotHeatmap("##ti", img.data(), T.ny, T.nx, 0, std::max(1, T.ncat() - 1), nullptr, ImPlotPoint(0, 0), ImPlotPoint(T.nx, T.ny));
                ImPlot::PopColormap();
                ImPlot::EndPlot();
            }
        }
        if (S.blocked.empty()) ImGui::TextDisabled("Block the wells first.");
        else if (S.fs.method == 3 && !S.ti) ImGui::TextDisabled("Generate or load a training image first.");
        else if (ui::accentButton(ICON_PLAY " Simulate facies")) {
            const GeoGrid g = S.g; const auto B = S.blocked; FaciesSpec f = S.fs;
            if (f.method == 2 && f.codes.empty()) { f.codes = {0, 1}; }
            const std::shared_ptr<const TrainingImage> ti = S.ti; MpsSpec mp = S.mps; mp.seed = f.method == 3 ? S.mps.seed : mp.seed;
            launch(f.method == 3 ? "MPS facies (Direct Sampling)" : "Facies modelling", [g, B, f, ti, mp]() -> std::function<void(GeoState&)> {
                FaciesResult r = (f.method == 3 && ti) ? modelFaciesMps(g, B, f, *ti, mp) : modelFacies(g, B, f);
                return [r](GeoState& S) {
                    S.fr = r; S.hasFacies = true; S.hasPetro = S.hasVol = S.hasUp = false;
                    std::string m = "facies:"; char b[64];
                    for (size_t k = 0; k < r.codes.size(); ++k) { std::snprintf(b, sizeof b, " code %d %.1f %%", r.codes[k], 100 * r.proportions[k]); m += b; }
                    m += ", honoured " + std::to_string(r.conditioned) + " cells";
                    S.msg(m); S.gridProp = 1;
                };
            });
        }
    } else if (S.step == 6) {   // ---------------- petrophysics
        static const char* pm[] = {"Kriging", "Sequential Gaussian simulation"};
        if (ui::beginProps("##ps")) {
            ui::propString("porosity log", &S.ps.curve);
            ui::propCombo("method", &S.ps.method, pm, 2);
            ui::propD("range major", &S.ps.vg.rangeMax, "%.0f"); ui::propD("range minor", &S.ps.vg.rangeMin, "%.0f"); ui::propD("range vertical", &S.ps.vg.rangeVert, "%.1f");
            ui::propBool("per facies", &S.ps.perFacies);
            ui::propBool("seismic secondary (co-simulation)", &S.useSeis, "collocated co-simulation with the RMS amplitude at the top of each column");
            if (S.useSeis) ui::propD("correlation", &S.corr, "%.2f");
            ui::propBool("NTG from facies (code 0 = non-net)", &S.ntgFromFacies);
            ui::propD("kv / kh", &S.kvkh, "%.3f");
            ui::propD("perm scatter (log10 sd)", &S.scatter, "%.2f");
            int sd = int(S.ps.seed); ui::propI("seed", &sd, 1, 1 << 30); S.ps.seed = unsigned(sd);
            ui::endProps();
        }
        if (S.userFit.n > 0) {   // C40
            char lb[160]; std::snprintf(lb, sizeof lb, "Use the cross-plot transform: log10 k = %.3f + %.3f phi (r2 %.2f)", S.userFit.a, S.userFit.b, S.userFit.r2);
            ImGui::Checkbox(lb, &S.useUserFit);
            ui::hint("Set on Data analysis > Cross-plot. Off: a least-squares fit per facies of the blocked PERM against the modelled curve.");
        }
        if (S.blocked.empty()) ImGui::TextDisabled("Block the wells first.");
        else if (ui::accentButton(ICON_PLAY " Model porosity + permeability")) {
            const GeoGrid g = S.g; const auto B = S.blocked; PetroSpec p = S.ps; const bool hasF = S.hasFacies; const FaciesResult fr = S.fr;
            const bool sec = S.useSeis && S.hasCube; const SeismicCube cube = sec ? S.cube : SeismicCube(); const VelocityModel vm = S.vm; const double corr = S.corr;
            const bool ntgF = S.ntgFromFacies; const double kvkh = S.kvkh, sc = S.scatter;
            const bool uf = S.useUserFit; const PoroPerm ufit = S.userFit;
            launch("Petrophysical modelling", [=]() mutable -> std::function<void(GeoState&)> {
                if (sec) {
                    const int nx = g.grid.nx, ny = g.grid.ny;
                    std::vector<P3> colTop;
                    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { const size_t c = size_t(i + nx * j); colTop.push_back({g.center[c][0], g.center[c][1], g.top[c]}); }
                    const auto am = extractAttribute(cube, vm, colTop, 1, 12);
                    p.secondary.assign(g.bulk.size(), 0.0);
                    for (size_t c = 0; c < g.bulk.size(); ++c) p.secondary[c] = am[c % size_t(nx * ny)];
                    p.corrSecondary = corr;
                }
                PetroResult r = modelProperty(g, B, hasF ? &fr.facies : nullptr, p);
                auto fits = fitPoroPerm(B, p.curve, "PERM", hasF ? &fr.facies : nullptr, fr.codes);
                if (uf) for (auto& f : fits) f = ufit;   // C40 transform chosen on the cross-plot
                auto kx = permFromPoro(r.v, fits, hasF ? &fr.facies : nullptr, fr.codes, sc, p.seed + 1);
                std::vector<double> kz(kx.size()), ntg(kx.size(), 1.0);
                for (size_t c = 0; c < kx.size(); ++c) { kz[c] = kvkh * kx[c]; if (ntgF && hasF && fr.facies[c] == 0) ntg[c] = 0.0; }
                return [r, fits, kx, kz, ntg, fr](GeoState& S) {
                    S.phi = r; S.fits = fits; S.kx = kx; S.kz = kz; S.ntg = ntg; S.hasPetro = true; S.hasVol = S.hasUp = false; S.gridProp = 2;
                    char b[200]; std::snprintf(b, sizeof b, "porosity mean %.4f sd %.4f (%d data) - %s", r.mean, r.sd, r.data, r.note.c_str()); S.msg(b);
                    for (size_t q = 0; q < fits.size(); ++q) { std::snprintf(b, sizeof b, "  poro-perm %s: log10 k = %.3f + %.3f phi (r2 %.2f, n %d)", q < fr.codes.size() ? ("facies " + std::to_string(fr.codes[q])).c_str() : "all", fits[q].a, fits[q].b, fits[q].r2, fits[q].n); S.msg(b); }
                };
            });
        }
        if (S.hasPetro) for (size_t q = 0; q < S.fits.size(); ++q) ImGui::BulletText("log10 k = %.2f + %.2f phi  (r2 %.2f)", S.fits[q].a, S.fits[q].b, S.fits[q].r2);
    } else if (S.step == 7) {   // ---------------- volumetrics
        if (ui::beginProps("##vs")) {
            ui::propD("oil-water contact", &S.vs.contacts.owc, "%.1f");
            ui::propD("gas-oil contact", &S.vs.contacts.goc, "%.1f", "-1e30: no gas cap");
            ui::propBool("saturation-height function", &S.vs.satHeight, "Brooks-Corey: Sw = Swirr + (1 - Swirr)(Pd / Pc)^lambda, Pc = drho g h above the contact");
            if (S.vs.satHeight) { ui::propD("Swirr", &S.vs.sh.swirr, "%.3f"); ui::propD("entry pressure Pd", &S.vs.sh.pd, "%.3f"); ui::propD("lambda", &S.vs.sh.lambda, "%.2f"); ui::propD("drho g (pressure / length)", &S.vs.sh.drhoOW, "%.4f"); }
            else ui::propD("Sw (constant)", &S.vs.swConst, "%.3f");
            ui::propD("Bo", &S.vs.bo, "%.3f"); ui::propD("Bg", &S.vs.bg, "%.5f");
            ui::propI("porosity realisations", &S.nReal, 1, 200);
            ui::endProps();
        }
        geoVolumeLimitsUi(S);   // C43 polygon limit, contacts per fault block
        if (!S.hasPetro) ImGui::TextDisabled("Model the petrophysics first.");
        else if (ui::accentButton(ICON_PLAY " Compute volumes")) {
            const GeoGrid g = S.g; const auto B = S.blocked; const PetroSpec p = S.ps; const bool hasF = S.hasFacies; const FaciesResult fr = S.fr; const auto phi = S.phi.v, ntg = S.ntg; VolSpec v = S.vs; v.field = S.field; const int nr = S.nReal; geoPrepareVolSpec(S, v);
            launch("Volumetrics", [=]() -> std::function<void(GeoState&)> {
                VolResult V = volumetrics(g, phi, ntg, hasF ? &fr.facies : nullptr, v);
                std::vector<double> rv;
                for (int r = 0; r < nr; ++r) { PetroSpec q = p; q.seed = p.seed + 1000 + unsigned(r); rv.push_back(volumetrics(g, modelProperty(g, B, hasF ? &fr.facies : nullptr, q).v, ntg, nullptr, v).total.stoiip); }
                Percentiles P = percentiles(rv);
                return [V, P](GeoState& S) { S.vol = V; S.pct = P; S.hasVol = true; char b[160]; std::snprintf(b, sizeof b, "STOIIP %.4g, GIIP %.4g; realisations P90 %.4g P50 %.4g P10 %.4g", V.total.stoiip, V.total.giip, P.p90, P.p50, P.p10); S.msg(b); };
            });
        }
        if (S.hasVol) {
            const char* su = S.field ? "stb" : "sm3";
            auto table = [&](const char* id, const std::vector<VolRow>& rows) {
                if (rows.empty()) return;
                if (ui::beginTable(id, 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
                    ui::tableColumn(id + 2); ui::tableColumn("bulk"); ui::tableColumn("pore"); ui::tableColumn(su);
                    ImGui::TableHeadersRow();
                    for (const auto& r : rows) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.key.c_str());
                        ImGui::TableNextColumn(); ImGui::Text("%.4g", r.bulk);
                        ImGui::TableNextColumn(); ImGui::Text("%.4g", r.pore);
                        ImGui::TableNextColumn(); ImGui::Text("%.4g", r.stoiip);
                    }
                    ui::endTable();
                }
            };
            table("##zone", S.vol.byZone); table("##block", S.vol.byBlock); table("##facies", S.vol.byFacies);
            ImGui::Text("Total STOIIP %.4g %s", S.vol.total.stoiip, su);
            ImGui::Text("P90 %.4g   P50 %.4g   P10 %.4g", S.pct.p90, S.pct.p50, S.pct.p10);
        }
        // ---- C36-4 uncertainty loop: structure -> grid -> facies -> petrophysics -> volumes per realisation
        ImGui::Separator();
        ImGui::TextDisabled("Uncertainty loop (structure, contacts, petrophysics; tornado)");
        if (ui::beginProps("##unc")) {
            ui::propI("realisations", &S.unc.nReal, 2, 2000);
            int sd = int(S.unc.seed); ui::propI("seed##unc", &sd, 1, 1 << 30); S.unc.seed = unsigned(sd);
            ui::propBool("Latin hypercube", &S.unc.latinHypercube, "stratified sampling of the input distributions (lower noise for few realisations)");
            ui::propBool("new facies / porosity seeds", &S.unc.stochastic, "geostatistical uncertainty: every realisation re-simulates facies and porosity");
            ui::propD("structural sd (depth)", &S.unc.zSd, "%.2f", "Gaussian random field added to all horizons, zero at the well tops (0 = fixed structure; > 0 rebuilds the grid per realisation)");
            if (S.unc.zSd > 0) ui::propD("structural range", &S.unc.zRange, "%.0f", "practical range of the Gaussian covariance");
            ui::propI("tornado: stochastic runs", &S.unc.tornadoStochastic, 0, 100);
            ui::endProps();
        }
        {   // triangular inputs (fixed = not varied)
            const double swB = S.vs.satHeight ? S.vs.sh.swirr : S.vs.swConst;
            struct Row { const char* name; Tri* t; double base, dlo, dhi; const char* fmt; };
            Row rows[] = {{"OWC depth", &S.unc.owc, S.vs.contacts.owc, S.field ? -60.0 : -20.0, S.field ? 60.0 : 20.0, "%.1f"},
                          {"GOC depth", &S.unc.goc, S.vs.contacts.goc, S.field ? -60.0 : -20.0, S.field ? 60.0 : 20.0, "%.1f"},
                          {"porosity multiplier", &S.unc.poroMult, 1.0, -0.1, 0.1, "%.3f"},
                          {"NTG multiplier", &S.unc.ntgMult, 1.0, -0.1, 0.1, "%.3f"},
                          {S.vs.satHeight ? "Swirr" : "Sw", &S.unc.sw, swB, -0.05, 0.05, "%.3f"},
                          {"Bo", &S.unc.bo, S.vs.bo, -0.05 * S.vs.bo, 0.05 * S.vs.bo, "%.4f"}};
            if (ui::beginTable("##tri", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable)) {
                ui::tableColumn("input"); ui::tableColumn("vary"); ui::tableColumn("low"); ui::tableColumn("most likely"); ui::tableColumn("high");
                ImGui::TableHeadersRow();
                for (auto& r : rows) {
                    ImGui::PushID(r.name);
                    const bool gocOff = r.t == &S.unc.goc && S.vs.contacts.goc < -1e29;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::TextUnformatted(r.name);
                    ImGui::TableNextColumn();
                    bool on = !r.t->fixed();
                    ImGui::BeginDisabled(gocOff);
                    if (ImGui::Checkbox("##v", &on)) *r.t = on ? Tri{r.base + r.dlo, r.base, r.base + r.dhi} : Tri{};
                    ImGui::EndDisabled();
                    ImGui::BeginDisabled(!on);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##lo", &r.t->lo, 0, 0, r.fmt);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##mo", &r.t->mode, 0, 0, r.fmt);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##hi", &r.t->hi, 0, 0, r.fmt);
                    ImGui::EndDisabled();
                    if (on && r.t->hi < r.t->lo) std::swap(r.t->lo, r.t->hi);
                    if (on) r.t->mode = std::clamp(r.t->mode, r.t->lo, r.t->hi);
                    ImGui::PopID();
                }
                ui::endTable();
            }
        }
        if (!S.hasPetro) ImGui::TextDisabled("Run the workflow up to the petrophysics first (its settings are re-used).");
        else if (S.fs.method == 3 && S.hasFacies && !S.ti) ImGui::TextDisabled("MPS facies need a training image (Facies step).");
        else if (ui::accentButton(ICON_SHUFFLE " Run the uncertainty loop")) {
            UncInputs in;
            in.gs = S.g.spec; in.hz = S.hz; in.faults = S.faults; in.wells = S.wells;
            in.facies = S.hasFacies; in.fs = S.fs; in.ps = S.ps; in.vs = S.vs; in.vs.field = S.field;
            if (in.fs.method == 2 && in.fs.codes.empty()) in.fs.codes = {0, 1};
            in.ti = S.ti; in.mps = S.mps; in.ntgFromFacies = S.ntgFromFacies && S.hasFacies;
            if (!S.ntgFromFacies && S.ntg.size() == S.g.bulk.size() && S.unc.zSd <= 0) in.ntgConst = S.ntg;
            const UncSpec u = S.unc;
            S.uncIn = in;
            launch("Uncertainty loop (" + std::to_string(u.nReal) + " realisations)", [in, u]() -> std::function<void(GeoState&)> {
                UncResult R = uncertaintyLoop(in, u);
                return [R](GeoState& S) {
                    if (!R.err.empty()) { S.msg("Error: uncertainty loop: " + R.err); return; }
                    S.ur = R; S.hasUnc = true;
                    std::string l; for (char c : R.report) { if (c == '\n') { S.msg(l); l.clear(); } else l += c; }
                    if (!l.empty()) S.msg(l);
                };
            });
        }
        if (S.hasUnc && !S.ur.real.empty()) {
            const char* su = S.field ? "stb" : "sm3";
            ImGui::Text("STOIIP  P90 %.4g   P50 %.4g   P10 %.4g %s   (base %.4g, mean %.4g)", S.ur.stoiip.p90, S.ur.stoiip.p50, S.ur.stoiip.p10, su, S.ur.base, S.ur.stoiip.mean);
            if (S.ur.giip.p10 > 0) ImGui::Text("GIIP    P90 %.4g   P50 %.4g   P10 %.4g", S.ur.giip.p90, S.ur.giip.p50, S.ur.giip.p10);
            // distribution: exceedance curve (industry P90 / P50 / P10 convention)
            std::vector<double> v; for (const auto& r : S.ur.real) v.push_back(r.stoiip);
            std::sort(v.begin(), v.end());
            std::vector<double> ex(v.size()); for (size_t q = 0; q < v.size(); ++q) ex[q] = 100.0 * (1.0 - (q + 0.5) / double(v.size()));
            if (ImPlot::BeginPlot("STOIIP distribution", ImVec2(-1, fs * 13))) {
                ImPlot::SetupAxes(su, "probability of exceedance (%)"); setupSI(ImAxis_X1);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 100, ImPlotCond_Always);
                ImPlot::PlotStairs("realisations", v.data(), ex.data(), int(v.size()), lineSpec(ImVec4(0.30f, 0.65f, 1.0f, 1), 2.0f));
                ImPlot::PlotScatter("##pts", v.data(), ex.data(), int(v.size()), markSpec(ImPlotMarker_Circle, 3, ImVec4(0.30f, 0.65f, 1.0f, 1)));
                const double px[3] = {S.ur.stoiip.p90, S.ur.stoiip.p50, S.ur.stoiip.p10}, py[3] = {90, 50, 10};
                ImPlot::PlotScatter("P90 / P50 / P10", px, py, 3, markSpec(ImPlotMarker_Diamond, 7, ImVec4(1, 0.75f, 0.2f, 1)));
                ImPlot::PlotInfLines("base case", &S.ur.base, 1, lineSpec(ImVec4(1, 1, 1, 0.7f), 1.5f));
                ImPlot::EndPlot();
            }
            // tornado (swing vs the base case, sorted by swing)
            const int nb = int(S.ur.tornado.size());
            if (nb > 0 && S.ur.base > 0 && ImPlot::BeginPlot("Tornado", ImVec2(-1, fs * (4.0f + 1.6f * nb)))) {
                const size_t nbs = static_cast<size_t>(nb); std::vector<double> lo(nbs), hi(nbs), pos(nbs);
                std::vector<std::string> nm(nbs); std::vector<const char*> lab(nbs);
                for (int q = 0; q < nb; ++q) {
                    const auto& t = S.ur.tornado[size_t(q)];
                    pos[size_t(q)] = nb - 1 - q;   // largest swing on top
                    lo[size_t(q)] = 100 * (t.low / S.ur.base - 1); hi[size_t(q)] = 100 * (t.high / S.ur.base - 1);
                    nm[size_t(q)] = t.name.substr(0, t.name.find(" (")); lab[size_t(q)] = nm[size_t(q)].c_str();
                }
                ImPlot::SetupAxes("STOIIP change vs base case (%)", nullptr, 0, ImPlotAxisFlags_NoGridLines);
                ImPlot::SetupAxisTicks(ImAxis_Y1, pos.data(), nb, lab.data());
                ImPlot::SetupAxisLimits(ImAxis_Y1, -0.7, nb - 0.3, ImPlotCond_Always);
                ImPlotSpec bl; bl.FillColor = ImVec4(0.85f, 0.35f, 0.30f, 1); bl.LineColor = bl.FillColor; bl.Flags = ImPlotBarsFlags_Horizontal;
                ImPlotSpec bh = bl; bh.FillColor = ImVec4(0.30f, 0.70f, 0.40f, 1); bh.LineColor = bh.FillColor;
                ImPlot::PlotBars("low", lo.data(), pos.data(), nb, 0.6, bl);
                ImPlot::PlotBars("high", hi.data(), pos.data(), nb, 0.6, bh);
                const double z0 = 0; ImPlot::PlotInfLines("##base", &z0, 1, lineSpec(ImVec4(1, 1, 1, 0.8f), 1.5f));
                ImPlot::EndPlot();
            }
            // rebuild a realisation into the workspace
            auto loadReal = [&](int idx, const char* tag) {
                if (idx < 0 || idx >= int(S.ur.real.size())) return;
                const UncInputs in = S.uncIn; const UncSample smp = S.ur.real[size_t(idx)].s;
                const auto fits = S.fits; const double sc = S.scatter, kvkh = S.kvkh; const std::string name = std::string(tag) + " (#" + std::to_string(idx + 1) + ")";
                launch("Rebuilding realisation " + name, [in, smp, fits, sc, kvkh, name]() -> std::function<void(GeoState&)> {
                    std::string e;
                    UncModel M = realise(in, smp, e);
                    std::vector<double> kx, kz;
                    std::vector<Blocked> B;
                    if (e.empty()) {
                        kx = permFromPoro(M.poro, fits, in.facies ? &M.fr.facies : nullptr, M.fr.codes, sc, smp.petroSeed + 1);
                        kz.resize(kx.size()); for (size_t c = 0; c < kx.size(); ++c) kz[c] = kvkh * kx[c];
                        B = blockWells(M.g, in.wells);
                    }
                    return [M, kx, kz, B, e, name, in](GeoState& S) {
                        if (!e.empty()) { S.msg("Error: " + e); return; }
                        S.g = M.g; S.hasGrid = true; S.blocked = B;
                        S.fr = M.fr; S.hasFacies = in.facies;
                        S.phi.v = M.poro; double m = 0; int n = 0; for (double p : M.poro) if (p > 0) { m += p; ++n; } S.phi.mean = n ? m / n : 0; S.phi.note = "realisation " + name;
                        S.kx = kx; S.kz = kz; S.ntg = M.ntg; S.hasPetro = true;
                        S.vol = M.vol; S.hasVol = true; S.hasUp = false;
                        char b[200]; std::snprintf(b, sizeof b, "realisation %s loaded: STOIIP %.4g (grid, facies, porosity, permeability; horizons unchanged)", name.c_str(), M.vol.total.stoiip); S.msg(b);
                    };
                });
            };
            if (ImGui::Button(TR("Load P90"))) loadReal(S.ur.iP90, "P90");
            ImGui::SameLine(); if (ImGui::Button(TR("Load P50"))) loadReal(S.ur.iP50, "P50");
            ImGui::SameLine(); if (ImGui::Button(TR("Load P10"))) loadReal(S.ur.iP10, "P10");
            ui::hint("Rebuilds that realisation exactly (structural field, contacts, seeds) as the current geomodel - then upscale / export / send it to simulation");
            ImGui::SameLine();
            { std::string p = "uncertainty.csv";
              if (ImGui::Button(TR(ICON_SAVE " CSV...")) && fileDialog(true, "CSV (*.csv)\0*.csv\0", "csv", p)) {
                  FILE* f = std::fopen(p.c_str(), "w");
                  if (f) {
                      std::fprintf(f, "realisation,stoiip,giip,pore,bulk,owc,goc,poro_mult,ntg_mult,sw,bo,field_seed,facies_seed,petro_seed\n");
                      for (size_t q = 0; q < S.ur.real.size(); ++q) { const auto& r = S.ur.real[q];
                          std::fprintf(f, "%zu,%.8g,%.8g,%.8g,%.8g,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%u,%u,%u\n", q + 1, r.stoiip, r.giip, r.pore, r.bulk, r.s.owc, r.s.goc, r.s.poroMult, r.s.ntgMult, r.s.sw, r.s.bo, r.s.fieldSeed, r.s.faciesSeed, r.s.petroSeed); }
                      std::fprintf(f, "\ntornado,low,high\n");
                      for (const auto& t : S.ur.tornado) std::fprintf(f, "\"%s\",%.8g,%.8g\n", t.name.c_str(), t.low, t.high);
                      std::fclose(f); S.msg("written " + p);
                  } else S.msg("Error: cannot write " + p);
              } }
        }
    } else if (S.step == 8) {   // ---------------- upscaling
        static const char* um[] = {"arithmetic", "harmonic", "geometric", "flow-based (TPFA)"};
        if (ui::beginProps("##us")) {
            ui::propI("factor I", &S.us.fx, 1, 50); ui::propI("factor J", &S.us.fy, 1, 50); ui::propI("factor K (inside each zone)", &S.us.fz, 1, 50);
            ui::propCombo("permeability", &S.us.permMethod, um, 4);
            ui::endProps();
        }
        if (!S.hasPetro) ImGui::TextDisabled("Model the petrophysics first.");
        else if (ui::accentButton(ICON_PLAY " Upscale")) {
            const GeoGrid g = S.g; const auto phi = S.phi.v, ntg = S.ntg, kx = S.kx, kz = S.kz; const UpscaleSpec us = S.us;
            launch("Upscaling", [=]() -> std::function<void(GeoState&)> {
                Upscaled U = upscale(g, phi, ntg, kx, kx, kz, us);
                return [U](GeoState& S) { S.up = U; S.hasUp = true; S.msg("upscaled to " + std::to_string(U.g.grid.nx) + " x " + std::to_string(U.g.grid.ny) + " x " + std::to_string(U.g.grid.nz) + " (pore volume conserved)"); };
            });
        }
        if (S.hasUp) ImGui::Text("coarse grid %d x %d x %d", S.up.g.grid.nx, S.up.g.grid.ny, S.up.g.grid.nz);
    } else {   // ---------------- export
        static int which = 1;
        ImGui::RadioButton(TR("fine geomodel"), &which, 0); ImGui::SameLine(); ImGui::RadioButton(TR("upscaled grid"), &which, 1);
        if (which == 1 && !S.hasUp) ImGui::TextDisabled("(no upscaled grid yet - the fine grid is used)");
        const bool up = which == 1 && S.hasUp;
        if (!S.hasPetro) ImGui::TextDisabled("Model the petrophysics first.");
        else {
            const GeoGrid& G = up ? S.up.g : S.g;
            const std::vector<double>& po = up ? S.up.poro : S.phi.v; const std::vector<double>& nt = up ? S.up.ntg : S.ntg;
            const std::vector<double>& kx = up ? S.up.kx : S.kx; const std::vector<double>& kz = up ? S.up.kz : S.kz;
            std::string p;
            if (ImGui::Button(TR(ICON_SAVE " Write GRDECL")) && fileDialog(true, "GRDECL (*.GRDECL)\0*.GRDECL;*.grdecl\0", "GRDECL", p)) {
                const bool ok = writeGrdecl(p, G, {{"PORO", &po}, {"NTG", &nt}, {"PERMX", &kx}, {"PERMY", &kx}, {"PERMZ", &kz}});
                S.msg(ok ? "written " + p : "Error: cannot write " + p);
            }
            ui::hint("SPECGRID / COORD / ZCORN / ACTNUM / PORO / NTG / PERMX-Z / FAULTS / MULTFLT - INCLUDE it in the GRID section of any ECLIPSE / OPM deck");
            if (ui::accentButton(ICON_LOG_OUT " Send to the Model workspace (simulation case)")) {
                Case c = case_;
                // the fluid / rock data come from the current case when its unit system matches, else from the
                // channelised example in the project units (three-phase black oil)
                if ((c.units == UnitSystem::Field) != S.field) {
                    c = makeExample("channels");
                    if ((c.units == UnitSystem::Field) != S.field) S.msg("Warning: the fluid template uses other units than the geomodel - check the PVT data");
                }
                c.wells.clear();
                toCase(c, G, po, nt, kx, kx, kz, S.wells);
                c.propModels.clear();
                if (S.hasPetro && !S.blocked.empty()) {   // C40: porosity recipe for Studies (variogram range / realisation parameters)
                    PropModelInput m;
                    m.property = 0; m.method = S.ps.method == 0 ? 0 : 1; m.model = int(S.ps.vg.model);
                    m.nugget = S.ps.vg.nugget; m.sill = S.ps.vg.sill; m.rangeMax = S.ps.vg.rangeMax; m.rangeMin = S.ps.vg.rangeMin;
                    m.rangeVert = S.ps.vg.rangeVert; m.azimuthDeg = S.ps.vg.azimuthDeg; m.seed = S.ps.seed;
                    m.minV = std::max(0.0, S.ps.minV); m.maxV = std::min(1.0, S.ps.maxV);
                    for (const auto& b : S.blocked) {
                        const auto it = b.v.find(S.ps.curve);
                        if (it == b.v.end() || !std::isfinite(it->second) || b.cell < 0 || size_t(b.cell) >= S.g.center.size()) continue;
                        const P3& q = S.g.center[size_t(b.cell)];
                        m.data.push_back({q[0], q[1], q[2], it->second});
                    }
                    const auto ff = fitPoroPerm(S.blocked, S.ps.curve, "PERM", nullptr, {});
                    const PoroPerm pp = S.useUserFit ? S.userFit : (ff.empty() ? PoroPerm() : ff[0]);
                    m.permFromPoro = pp.n > 1; m.ppA = pp.a; m.ppB = pp.b; m.kvkh = S.kvkh;
                    m.source = "Geomodel";
                    if (!m.data.empty()) c.propModels.push_back(m);
                }
                if (S.hasFacies && !S.blocked.empty() && S.fr.codes.size() >= 2 && S.fr.proportions.size() == S.fr.codes.size()) {
                    // C41: facies recipe (truncated Gaussian re-simulation) - per-facies mean porosity and geometric-mean
                    // permeability of the blocked logs. It replaces the porosity recipe: a later SGS porosity model would
                    // overwrite the facies contrast (the base case keeps the full geomodel heterogeneity either way)
                    PropModelInput f;
                    f.property = 3; f.method = 1; f.model = int(S.fs.vg.model);
                    f.nugget = S.fs.vg.nugget; f.sill = S.fs.vg.sill; f.rangeMax = S.fs.vg.rangeMax; f.rangeMin = S.fs.vg.rangeMin;
                    f.rangeVert = S.fs.vg.rangeVert; f.azimuthDeg = S.fs.vg.azimuthDeg; f.seed = S.fs.seed; f.kvkh = S.kvkh;
                    const size_t K = S.fr.codes.size();
                    f.facProp = S.fr.proportions;
                    std::vector<double> sp(K, 0.0), sk(K, 0.0); std::vector<int> np(K, 0), nk(K, 0);
                    for (const auto& b : S.blocked) {
                        const auto fc = b.v.find(S.fs.curve);
                        if (fc == b.v.end() || !std::isfinite(fc->second) || b.cell < 0 || size_t(b.cell) >= S.g.center.size()) continue;
                        const auto it = std::find(S.fr.codes.begin(), S.fr.codes.end(), int(std::lround(fc->second)));
                        if (it == S.fr.codes.end()) continue;
                        const size_t k = size_t(it - S.fr.codes.begin());
                        const P3& q = S.g.center[size_t(b.cell)];
                        f.data.push_back({q[0], q[1], q[2], double(k + 1)});
                        const auto ph = b.v.find(S.ps.curve);
                        if (ph != b.v.end() && std::isfinite(ph->second)) { sp[k] += ph->second; ++np[k]; }
                        const auto pk = b.v.find("PERM");
                        if (pk != b.v.end() && std::isfinite(pk->second) && pk->second > 0) { sk[k] += std::log(pk->second); ++nk[k]; }
                    }
                    for (size_t k = 0; k < K; ++k) {
                        f.facPoro.push_back(np[k] ? sp[k] / np[k] : -1.0);
                        f.facPerm.push_back(nk[k] ? std::exp(sk[k] / nk[k]) : -1.0);
                    }
                    f.source = "Geomodel facies";
                    if (!f.data.empty()) {
                        c.propModels.erase(std::remove_if(c.propModels.begin(), c.propModels.end(), [](const PropModelInput& q) { return q.property == 0; }), c.propModels.end());
                        c.propModels.push_back(f);
                    }
                }
                double zt = 1e300; for (size_t q = 0; q < G.top.size(); ++q) if (G.grid.actnum.empty() || G.grid.actnum[q]) zt = std::min(zt, G.top[q]);
                c.equil.datumDepth = S.vs.contacts.owc < 1e29 ? S.vs.contacts.owc : zt; c.equil.owc = S.vs.contacts.owc < 1e29 ? S.vs.contacts.owc : zt + 1000;
                c.equil.goc = S.vs.contacts.goc > -1e29 ? S.vs.contacts.goc : zt - (S.field ? 330 : 100); c.equil.rsvdDepth.clear(); c.equil.rsvdRs.clear();
                // hydrostatic datum pressure (normally pressured reservoir; edit in Initialization)
                c.equil.datumPressure = S.field ? 14.7 + 0.433 * c.equil.datumDepth : 1.01325 + 0.0981 * c.equil.datumDepth;
                if (c.schedule.empty()) c.schedule.resize(12);   // 12 x 30 d
                for (auto& st : c.schedule) { st.controls.clear(); st.gprod.clear(); st.ginj.clear(); }
                for (const auto& w : c.wells) { WellControl wc; wc.well = w.name; wc.mode = Control::LRAT; wc.lrat = S.field ? 1000 : 150; wc.bhp = S.field ? 1500 : 100; c.schedule[0].controls.push_back(wc); }
                c.title = "ResGeo geomodel " + std::to_string(G.grid.nx) + "x" + std::to_string(G.grid.ny) + "x" + std::to_string(G.grid.nz);
                case_ = c; markDirty(true);
                geoMode_ = false; workspace_ = Workspace::Model;
                notify("Geomodel sent to the Model workspace: corner-point grid, properties, faults and " + std::to_string(c.wells.size()) + " wells (LRAT producers - edit them in Wells / Schedule; EQUIL datum = OWC)");
            }
            ui::hint("Keeps the current case's PVT / SCAL; replaces the grid, wells and well controls; EQUIL from the volumetrics contacts");
        }
    }
    ImGui::EndDisabled();
    ImGui::EndChild();
    geoImportWindow(S);   // C42 column-assignment dialog (modal)
    geoWellManager(S);    // C42 well manager window
    ImGui::End();
}
