// C43 Geomodel tools UI (docs/PETREL_PARITY.md P4 / P6 / P8 / P9):
//   polygons (digitise on the map, import / export, clip horizons, limit volumetrics), per-segment contacts,
//   2-D mis-tie correction, layer-cake velocity model with checkshots, seismic calculator and cropping
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "platform.h"
#include "ressim/geo_tools.h"
#include "ressim/geo_analysis.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <set>
#include <string>
#include <vector>

namespace {
using namespace ressim;
using namespace ressim::geo;
using GPoly = ressim::geo::Polygon;   // windows.h declares a GDI Polygon()

ImPlotSpec lspec(ImVec4 c, float w) { ImPlotSpec s; s.LineColor = c; s.LineWeight = w; return s; }
ImPlotSpec mspec(ImVec4 c, float r) { ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = r; s.MarkerFillColor = c; s.MarkerLineColor = c; s.LineWeight = 0; return s; }
std::string base(const std::string& p) {
    const size_t s = p.find_last_of("/\\");
    std::string b = p.substr(s == std::string::npos ? 0 : s + 1);
    const size_t d = b.rfind('.');
    return d == std::string::npos ? b : b.substr(0, d);
}
std::vector<std::string> topNames(const App::GeoState& S) {
    std::set<std::string> u;
    for (const auto& w : S.wells) for (const auto& t : w.tops) u.insert(t.surface);
    return {u.begin(), u.end()};
}
bool strCombo(const char* id, std::string* v, const std::vector<std::string>& items) {
    bool ch = false;
    if (ImGui::BeginCombo(id, v->empty() ? "(none)" : v->c_str())) {
        if (ImGui::Selectable(TR("(none)"), v->empty())) { v->clear(); ch = true; }
        for (const auto& it : items) if (ImGui::Selectable(it.c_str(), it == *v)) { *v = it; ch = true; }
        ImGui::EndCombo();
    }
    return ch;
}
// index combo over names; with `none` an entry for -1 is offered first
bool vcombo(const char* id, int* v, const std::vector<std::string>& items, const char* none = nullptr) {
    bool ch = false;
    const char* cur = (*v >= 0 && *v < int(items.size())) ? items[size_t(*v)].c_str() : (none ? none : "");
    if (ImGui::BeginCombo(id, cur)) {
        if (none && ImGui::Selectable(none, *v < 0)) { *v = -1; ch = true; }
        for (int k = 0; k < int(items.size()); ++k) {
            ImGui::PushID(k);
            if (ImGui::Selectable(items[size_t(k)].c_str(), k == *v)) { *v = k; ch = true; }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    return ch;
}
}  // namespace

// ================================================================ map overlay: polygons + digitising
void geoMapOverlay(App::GeoState& S) {
    for (size_t p = 0; p < S.polys.size(); ++p) {
        const GPoly& P = S.polys[p];
        if (P.pts.size() < 2) continue;
        std::vector<double> xs, ys;
        for (const auto& q : P.pts) { xs.push_back(q[0]); ys.push_back(q[1]); }
        if (P.closed) { xs.push_back(P.pts[0][0]); ys.push_back(P.pts[0][1]); }
        const bool cur = int(p) == S.curPoly;
        const ImVec4 col = P.closed ? ImVec4(0.3f, 1.0f, 0.55f, cur ? 1.0f : 0.8f) : ImVec4(1.0f, 0.6f, 0.2f, cur ? 1.0f : 0.8f);
        ImPlot::PlotLine(("##pg" + std::to_string(p)).c_str(), xs.data(), ys.data(), int(xs.size()), lspec(col, cur ? 3.0f : 1.8f));
        ImPlot::PlotText(P.name.c_str(), P.pts[0][0], P.pts[0][1], ImVec2(0, 10));
    }
    if (S.pickMode != 4 && S.pickMode != 5) return;
    const bool closed = S.pickMode == 4;
    GPoly& D = S.polyDraft;
    std::vector<double> xs, ys;
    for (const auto& q : D.pts) { xs.push_back(q[0]); ys.push_back(q[1]); }
    if (ImPlot::IsPlotHovered()) {
        const ImPlotPoint mp = ImPlot::GetPlotMousePos();
        xs.push_back(mp.x); ys.push_back(mp.y);
        ImGui::SetTooltip("%s: click adds a point (%zu so far)\nright-click or Enter finishes, Backspace removes the last point, Esc cancels",
                          closed ? "Polygon" : "Intersection path", D.pts.size());
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && ImGui::GetIO().MouseDragMaxDistanceSqr[0] < 16.0f) D.pts.push_back({mp.x, mp.y});
    }
    if (!xs.empty()) {
        ImPlot::PlotLine("##draft", xs.data(), ys.data(), int(xs.size()), lspec(ImVec4(1, 1, 0.3f, 1), 2.0f));
        ImPlot::PlotScatter("##draftp", xs.data(), ys.data(), int(D.pts.size()), mspec(ImVec4(1, 1, 0.3f, 1), 4));
    }
    const bool finish = (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) || ImGui::IsKeyPressed(ImGuiKey_Enter);
    if (ImGui::IsKeyPressed(ImGuiKey_Backspace) && !D.pts.empty()) D.pts.pop_back();
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) { D = GPoly{}; S.pickMode = 0; }
    if (finish) {
        if (D.pts.size() >= (closed ? 3u : 2u)) {
            D.closed = closed;
            if (D.name.empty()) D.name = (closed ? "Polygon " : "Path ") + std::to_string(S.polys.size() + 1);
            S.polys.push_back(D);
            S.curPoly = int(S.polys.size()) - 1;
            char b[160];
            std::snprintf(b, sizeof b, "%s %s: %zu points, %s %.4g", closed ? "polygon" : "path", D.name.c_str(), D.pts.size(), closed ? "area" : "length",
                          closed ? polygonArea(D) : polylineLength(D));
            S.msg(b);
            if (!closed) { S.isecPath = S.curPoly; S.view = 6; }
        } else S.msg("need at least " + std::string(closed ? "3" : "2") + " points");
        D = GPoly{};
        S.pickMode = 0;
    }
}

// ================================================================ polygons panel (Data step)
void geoPolygonsPanel(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    if (S.openTool == 0) { ImGui::SetNextItemOpen(true); S.openTool = -1; }
    if (!ImGui::CollapsingHeader(TR(ICON_ROUTE " Polygons and intersection paths"))) return;
    if (ImGui::Button(S.pickMode == 4 ? ICON_X " Stop digitising" : ICON_PENCIL " Digitise polygon")) {
        S.pickMode = S.pickMode == 4 ? 0 : 4; S.polyDraft = GPoly{}; if (S.pickMode == 4) S.view = 0;
    }
    ImGui::SameLine();
    if (ImGui::Button(S.pickMode == 5 ? ICON_X " Stop" : ICON_SCAN_LINE " Digitise intersection path")) {
        S.pickMode = S.pickMode == 5 ? 0 : 5; S.polyDraft = GPoly{}; if (S.pickMode == 5) S.view = 0;
    }
    ui::hint("Click points on the Map view; right-click or Enter finishes. Polygons (closed) limit volumetrics and clip horizons; "
             "paths (open) define arbitrary intersections (Intersection view).");
    std::vector<std::string> files;
    if (ImGui::Button(TR(ICON_FOLDER_OPEN " Import...")) && plat::fileDialogMulti("Polygons (*.txt;*.pol;*.xyz;*.dat)\0*.txt;*.pol;*.xyz;*.dat\0All files (*.*)\0*.*\0", files)) {
        int n = 0;
        for (const auto& f : files) {
            std::string e;
            auto v = readPolygons(f, e);
            if (!e.empty()) S.msg("Error: " + e);
            for (auto& p : v) { if (v.size() == 1) p.name = base(f); S.polys.push_back(p); ++n; }
        }
        S.msg(std::to_string(n) + " polygon(s) imported");
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(S.polys.empty());
    if (ImGui::Button(TR(ICON_SAVE " Export..."))) {
        std::string p = "polygons.txt";
        if (plat::fileDialog(true, "Polygons (*.txt)\0*.txt\0", "txt", p)) S.msg(writePolygons(p, S.polys) ? "polygons written: " + p : "Error: cannot write " + p);
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if ((S.hasGrid || S.gs.nx > 0) && ImGui::Button(TR("From grid outline"))) {
        GeoGrid tmp; tmp.spec = S.hasGrid ? S.g.spec : S.gs;
        GPoly P; P.name = "Grid outline";
        const int ci[4] = {0, tmp.spec.nx, tmp.spec.nx, 0}, cj[4] = {0, 0, tmp.spec.ny, tmp.spec.ny};
        for (int q = 0; q < 4; ++q) { double x, y; tmp.cornerXY(ci[q], cj[q], x, y); P.pts.push_back({x, y}); }
        S.polys.push_back(P);
    }
    if (S.polys.empty()) { ImGui::TextDisabled("no polygons"); return; }
    if (ImGui::BeginTable("##polys", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("name", ImGuiTableColumnFlags_WidthStretch, 3);
        ImGui::TableSetupColumn("type", ImGuiTableColumnFlags_WidthStretch, 1.4f);
        ImGui::TableSetupColumn("points", ImGuiTableColumnFlags_WidthStretch, 1);
        ImGui::TableSetupColumn("area / length", ImGuiTableColumnFlags_WidthStretch, 2);
        ImGui::TableSetupColumn("##x", ImGuiTableColumnFlags_WidthFixed, fs * 2.4f);
        ImGui::TableHeadersRow();
        for (int p = 0; p < int(S.polys.size()); ++p) {
            GPoly& P = S.polys[size_t(p)];
            ImGui::PushID(p);
            ImGui::TableNextRow(); ImGui::TableNextColumn();
            char nb[96]; std::snprintf(nb, sizeof nb, "%s", P.name.c_str());
            ImGui::SetNextItemWidth(-1);
            if (ImGui::InputText("##n", nb, sizeof nb)) P.name = nb;
            if (ImGui::IsItemActivated()) S.curPoly = p;
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(P.closed ? "polygon" : "path")) P.closed = !P.closed && P.pts.size() >= 3;
            ImGui::TableNextColumn(); ImGui::Text("%zu", P.pts.size());
            ImGui::TableNextColumn(); ImGui::Text("%.4g", P.closed ? polygonArea(P) : polylineLength(P));
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_TRASH_2)) {
                S.polys.erase(S.polys.begin() + p);
                if (S.volPoly == p) S.volPoly = -1; else if (S.volPoly > p) --S.volPoly;
                if (S.isecPath == p) S.isecPath = -1; else if (S.isecPath > p) --S.isecPath;
                S.curPoly = -1;
                ImGui::PopID();
                break;
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    // clip a horizon by the current polygon
    if (!S.hz.empty()) {
        static int clipHz = 0, clipPoly = 0; static bool keepIn = true;
        std::vector<std::string> hn, pn;
        for (const auto& h : S.hz) hn.push_back(h.name);
        std::vector<int> closedIdx;
        for (int p = 0; p < int(S.polys.size()); ++p) if (S.polys[size_t(p)].closed) { closedIdx.push_back(p); pn.push_back(S.polys[size_t(p)].name); }
        if (!closedIdx.empty()) {
            clipHz = std::clamp(clipHz, 0, int(hn.size()) - 1); clipPoly = std::clamp(clipPoly, 0, int(pn.size()) - 1);
            ImGui::SetNextItemWidth(fs * 8); vcombo("##chz", &clipHz, hn); ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 8); vcombo("##cpo", &clipPoly, pn); ImGui::SameLine();
            ImGui::Checkbox(TR("keep inside"), &keepIn); ImGui::SameLine();
            if (ImGui::Button(TR(ICON_SCISSORS " Clip horizon"))) {
                S.pushUndo();
                S.hz[size_t(clipHz)] = clipSurface(S.hz[size_t(clipHz)], S.polys[size_t(closedIdx[size_t(clipPoly)])], keepIn);
                ++S.hzRev;
                S.msg("horizon " + S.hz[size_t(clipHz)].name + " clipped by " + pn[size_t(clipPoly)] + (S.hasGrid ? " - rebuild the grid to apply" : ""));
            }
        }
    }
}

// ================================================================ volumetrics limits (Volumetrics step)
void geoVolumeLimitsUi(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    std::vector<std::string> pn; std::vector<int> idx;
    for (int p = 0; p < int(S.polys.size()); ++p) if (S.polys[size_t(p)].closed) { idx.push_back(p); pn.push_back(S.polys[size_t(p)].name); }
    int sel = -1;
    for (size_t k = 0; k < idx.size(); ++k) if (idx[k] == S.volPoly) sel = int(k);
    ImGui::SetNextItemWidth(fs * 12);
    if (vcombo("Limit to polygon", &sel, pn, "(whole grid)")) S.volPoly = sel >= 0 ? idx[size_t(sel)] : -1;
    ui::hint("Only cells whose centre lies inside the polygon are counted (licence / lease / prospect boundaries). Digitise or "
             "import polygons on the Data step.");
    if (S.cf.on) { ImGui::Checkbox(TR("only cells passing the view filter (zones / blocks / facies)"), &S.volUseFilter); }
    if (S.hasGrid && S.g.nBlocks > 1) {
        ImGui::Checkbox(TR("contacts per fault block"), &S.volPerBlock);
        if (S.volPerBlock) {
            if (int(S.blockContacts.size()) != S.g.nBlocks) S.blockContacts.assign(size_t(S.g.nBlocks), S.vs.contacts);
            if (ImGui::BeginTable("##bc", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH)) {
                ImGui::TableSetupColumn("block"); ImGui::TableSetupColumn("OWC"); ImGui::TableSetupColumn("GOC");
                ImGui::TableHeadersRow();
                for (int b = 0; b < S.g.nBlocks; ++b) {
                    ImGui::PushID(b);
                    ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::Text("block %d", b + 1);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##o", &S.blockContacts[size_t(b)].owc, 0, 0, "%.1f");
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##g", &S.blockContacts[size_t(b)].goc, 0, 0, "%.1f");
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
        }
    }
}
void geoPrepareVolSpec(const App::GeoState& S, VolSpec& v) {
    v.mask.clear(); v.blockContacts.clear();
    if (S.hasGrid && S.volPoly >= 0 && S.volPoly < int(S.polys.size()) && S.polys[size_t(S.volPoly)].closed) v.mask = polygonMask(S.g, S.polys[size_t(S.volPoly)]);
    if (S.hasGrid && S.cf.on && S.volUseFilter) {
        const size_t n = S.g.bulk.size();
        if (v.mask.size() != n) v.mask.assign(n, 1);
        for (size_t c = 0; c < n; ++c) if (!geoCellVisible(S, c)) v.mask[c] = 0;
    }
    if (S.volPerBlock && S.hasGrid && int(S.blockContacts.size()) == S.g.nBlocks) v.blockContacts = S.blockContacts;
}

// ================================================================ 2-D mis-tie correction (Seismic step)
void geoMistiePanel(App::GeoState& S) {
    if (S.lines.size() < 2 || S.lineHz.empty()) return;
    const float fs = ImGui::GetFontSize();
    if (S.openTool == 1) { ImGui::SetNextItemOpen(true); S.openTool = -1; }
    if (!ImGui::CollapsingHeader(TR(ICON_GIT_COMPARE " Mis-tie correction"))) return;
    std::vector<std::string> ln;
    for (const auto& l : S.lines) ln.push_back(l.name);
    S.mistieRef = std::clamp(S.mistieRef, 0, int(ln.size()) - 1);
    ImGui::SetNextItemWidth(fs * 10); vcombo("reference line", &S.mistieRef, ln);
    ui::hint("Bulk time shift per line by least squares over the intersections of every horizon picked on both lines; the "
             "reference line keeps its times (one reference per group of connected lines). Apply shifts the line data and the picks.");
    if (ImGui::Button(TR(ICON_SIGMA " Compute"))) {
        S.mistie = fitMisties(S.lines, lineIntersections(S.lines), S.lineHz, S.mistieRef);
        S.mistieHas = true;
        char b[160]; std::snprintf(b, sizeof b, "mis-ties: %d pair(s), rms %.3f -> %.3f ms after the shifts", S.mistie.pairs, S.mistie.rmsBefore, S.mistie.rmsAfter);
        S.msg(b);
    }
    if (!S.mistieHas || S.mistie.shift.size() != S.lines.size()) return;
    ImGui::SameLine();
    if (ui::accentButton(ICON_CHECK " Apply shifts")) {
        applyLineShifts(S.lines, S.lineHz, S.mistie.shift);
        ++S.lineRev; S.mistieHas = false;
        S.msg("line shifts applied");
    }
    ImGui::Text("pairs %d   rms %.3f ms -> %.3f ms", S.mistie.pairs, S.mistie.rmsBefore, S.mistie.rmsAfter);
    if (ImGui::BeginTable("##mt", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH)) {
        ImGui::TableSetupColumn("line"); ImGui::TableSetupColumn("shift (ms)"); ImGui::TableSetupColumn("group");
        ImGui::TableHeadersRow();
        for (size_t l = 0; l < S.lines.size(); ++l) {
            ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(S.lines[l].name.c_str());
            ImGui::TableNextColumn(); ImGui::Text("%+.3f", S.mistie.shift[l]);
            ImGui::TableNextColumn();
            if (S.mistie.component[l] < 0) ImGui::TextDisabled("no ties"); else ImGui::Text("%d%s", S.mistie.component[l] + 1, int(l) == S.mistieRef ? " (ref)" : "");
        }
        ImGui::EndTable();
    }
}

// ================================================================ layer-cake velocity model (Seismic step)
void geoLayerCakePanel(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    if (S.openTool == 3) { ImGui::SetNextItemOpen(true); S.openTool = -1; }
    if (!ImGui::CollapsingHeader(TR(ICON_LAYERS " Layer-cake velocity model / checkshots"))) return;
    ui::hint("Interfaces top down, each a time surface (tracked 3-D horizon, a 2-D line horizon gridded in time, or an imported "
             "time grid) tied to a well top. Interval velocity per layer = 2000 dz / dt at the wells: the mean, an inverse-distance "
             "map, or a value you give. Depth surfaces are added to the horizons (replacing the same names); the residual "
             "correction makes them honour the tops.");
    std::vector<std::string> src;
    std::vector<int> srcId;
    if (S.hasTh) { src.push_back("Tracked 3-D horizon (" + S.thName + ")"); srcId.push_back(0); }
    for (size_t h = 0; h < S.lineHz.size(); ++h) { src.push_back("2-D lines: " + S.lineHz[h].name); srcId.push_back(1 + int(h)); }
    for (size_t h = 0; h < S.hz.size(); ++h) { src.push_back("Time grid: " + S.hz[h].name); srcId.push_back(1000 + int(h)); }
    if (src.empty()) { ImGui::TextDisabled("needs time interpretation (tracked 3-D horizon, 2-D line horizons) or imported time grids"); }
    const auto tops = topNames(S);
    static const char* modes[] = {"mean at wells", "map (inverse distance)", "user value"};
    if (ImGui::BeginTable("##lc", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("time surface", ImGuiTableColumnFlags_WidthStretch, 3);
        ImGui::TableSetupColumn("well top", ImGuiTableColumnFlags_WidthStretch, 2);
        ImGui::TableSetupColumn("interval velocity", ImGuiTableColumnFlags_WidthStretch, 2);
        ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch, 1.3f);
        ImGui::TableSetupColumn("##x", ImGuiTableColumnFlags_WidthFixed, fs * 2.2f);
        ImGui::TableHeadersRow();
        for (int r = 0; r < int(S.lc.size()); ++r) {
            auto& R = S.lc[size_t(r)];
            ImGui::PushID(r);
            ImGui::TableNextRow(); ImGui::TableNextColumn();
            int sel = -1; for (size_t k = 0; k < srcId.size(); ++k) if (srcId[k] == R.src) sel = int(k);
            ImGui::SetNextItemWidth(-1); if (vcombo("##s", &sel, src, "(choose)")) R.src = sel >= 0 ? srcId[size_t(sel)] : -1;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); strCombo("##t", &R.top, tops);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::Combo("##m", &R.mode, modes, 3);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::InputDouble("##v", &R.v, 0, 0, "%.0f");
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_X)) { S.lc.erase(S.lc.begin() + r); ImGui::PopID(); break; }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (ImGui::Button(TR(ICON_PLUS " Interface"))) { App::GeoState::LcRow R; R.src = srcId.empty() ? -1 : srcId[0]; if (tops.size() > S.lc.size()) R.top = tops[S.lc.size()]; S.lc.push_back(R); }
    static bool residual = true;
    ImGui::SameLine(); ImGui::Checkbox(TR("residual correction"), &residual);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); ImGui::InputDouble("datum depth", &S.vm.datum, 0, 0, "%.1f");
    ImGui::BeginDisabled(S.lc.empty() || S.job.valid());
    if (ui::accentButton(ICON_PLAY " Build depth surfaces")) {
        // time surfaces on one lattice (first horizon, or the bounding lattice of the line picks)
        LayerCakeInput in; in.datumZ = S.vm.datum; in.residual = residual;
        std::string err;
        Surface tmpl;
        if (!S.hz.empty()) tmpl = S.hz[0];
        for (const auto& R : S.lc) {
            Surface t;
            if (R.src == 0 && S.hasTh && S.hasCube) t = timeSurface(S.cube, S.th, S.thName);
            else if (R.src >= 1 && R.src < 1000 && size_t(R.src - 1) < S.lineHz.size()) {
                const auto pts = lineHorizonTimePoints(S.lines, S.lineHz[size_t(R.src - 1)], 2);
                if (pts.size() < 3) { err = "line horizon " + S.lineHz[size_t(R.src - 1)].name + " has too few picks"; break; }
                if (tmpl.nx == 0) tmpl = latticeFor(pts, 81);
                double x0 = 1e300, x1 = -1e300; for (auto& q : pts) { x0 = std::min(x0, q[0]); x1 = std::max(x1, q[0]); }
                geostat::Variogram vg; vg.model = geostat::VarioModel::Gaussian; vg.rangeMax = vg.rangeMin = std::max(1.0, 0.35 * (x1 - x0));
                t = gridPoints(pts, tmpl, vg, 16); t.name = S.lineHz[size_t(R.src - 1)].name;
            } else if (R.src >= 1000 && size_t(R.src - 1000) < S.hz.size()) t = S.hz[size_t(R.src - 1000)];
            else { err = "choose a time surface for every interface"; break; }
            in.twt.push_back(t); in.topName.push_back(R.top); in.mode.push_back(R.mode); in.vUser.push_back(R.v);
        }
        if (!err.empty()) S.msg("Error: " + err);
        else {
            const auto wells = S.wells;
            geoLaunch(S, "Layer-cake depth conversion", [in, wells]() -> std::function<void(App::GeoState&)> {
                LayerCakeResult R = layerCakeDepth(in, wells);
                return [R](App::GeoState& S) {
                    for (const auto& w : R.warnings) S.msg("Warning: " + w);
                    if (R.depth.empty()) return;
                    S.pushUndo();
                    for (const auto& d : R.depth) {
                        bool rep = false;
                        for (auto& h : S.hz) if (h.name == d.name) { h = d; rep = true; }
                        if (!rep) S.hz.push_back(d);
                    }
                    ++S.hzRev;
                    S.lcRes = R; S.lcHas = true;
                    char b[160]; std::snprintf(b, sizeof b, "layer cake: %zu depth surface(s); misfit at the tops %.2f -> %.2f (rms)", R.depth.size(), R.rmsBefore, R.rmsAfter);
                    S.msg(b);
                };
            });
        }
    }
    ImGui::EndDisabled();
    if (S.lcHas) {
        const auto& R = S.lcRes;
        for (size_t f = 0; f < R.vMean.size(); ++f)
            ImGui::Text("layer %zu: mean interval velocity %.1f from %d well(s)", f + 1, R.vMean[f], f < R.nWells.size() ? R.nWells[f] : 0);
        ImGui::Text("misfit at the tops: rms %.3f before, %.3f after the residual correction", R.rmsBefore, R.rmsAfter);
        if (!R.ties.empty() && ImGui::TreeNode(TR("ties at the wells"))) {
            if (ImGui::BeginTable("##lct", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH)) {
                for (const char* h : {"well", "interface", "top", "model", "corrected"}) ImGui::TableSetupColumn(h);
                ImGui::TableHeadersRow();
                for (const auto& t : R.ties) {
                    ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(t.well.c_str());
                    ImGui::TableNextColumn(); ImGui::Text("%d", t.iface + 1);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f", t.zTop);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f", t.zBefore);
                    ImGui::TableNextColumn(); ImGui::Text("%.2f", t.zAfter);
                }
                ImGui::EndTable();
            }
            ImGui::TreePop();
        }
    }
    // checkshots -> linear velocity function
    ImGui::SeparatorText(TR("Checkshots"));
    static char cs[2048] = "";
    ImGui::InputTextMultiline("##cs", cs, sizeof cs, ImVec2(-1, fs * 4), ImGuiInputTextFlags_AllowTabInput);
    ImGui::TextDisabled("depth  TWT(ms) per line (paste from Excel); fits v = v0 + k (z - datum) and sets the velocity model");
    if (ImGui::Button(TR(ICON_TRENDING_UP " Fit checkshots"))) {
        std::vector<double> z, t;
        std::string all(cs);
        for (char& c : all) if (c == ',' || c == ';' || c == '\t') c = ' ';
        size_t p = 0;
        while (p < all.size()) {
            size_t e = all.find('\n', p); if (e == std::string::npos) e = all.size();
            double a, b;
            if (std::sscanf(all.substr(p, e - p).c_str(), "%lf %lf", &a, &b) == 2) { z.push_back(a); t.push_back(b); }
            p = e + 1;
        }
        if (z.size() < 2) S.msg("Error: need at least two depth / TWT pairs");
        else {
            double rms = 0;
            const VelocityModel vm = fitCheckshots(z, t, S.vm.datum, &rms);
            S.vm.v0 = vm.v0; S.vm.k = vm.k;
            char b[160]; std::snprintf(b, sizeof b, "checkshots (%zu): v0 %.1f, k %.4f, rms %.2f ms - velocity model updated", z.size(), vm.v0, vm.k, rms);
            S.msg(b);
        }
    }
}

// ================================================================ seismic calculator / cropping (Seismic step)
void geoSeisCalcPanel(App::GeoState& S) {
    if (!S.hasCube) return;
    const float fs = ImGui::GetFontSize();
    if (S.openTool == 2) { ImGui::SetNextItemOpen(true); S.openTool = -1; }
    if (!ImGui::CollapsingHeader(TR(ICON_SIGMA " Seismic calculator / crop"))) return;
    if (S.segy) { ImGui::TextDisabled("out-of-core volume: the calculator and cropping work on the in-memory proxy only"); }
    static char expr[256] = "";
    if (!expr[0] && S.cube.name.find("calc") != std::string::npos) std::snprintf(expr, sizeof expr, "%s", "ENVELOPE * (T > 1700)");
    ImGui::SetNextItemWidth(-1);
    ImGui::InputTextWithHint("##sx", "e.g.  ENVELOPE * (T > 1500)      or      AMP / (RMS + 1e-6)", expr, sizeof expr);
    for (const auto& n : seismicCalculatorNames()) {
        if (ImGui::SmallButton(n.c_str())) { const size_t L = std::strlen(expr); std::snprintf(expr + L, sizeof expr - L, "%s%s", L ? " " : "", n.c_str()); }
        ImGui::SameLine();
        if (ImGui::GetContentRegionAvail().x < fs * 6) ImGui::NewLine();
    }
    ImGui::NewLine();
    ui::hint("Attributes used in the expression are computed on demand with the attribute window of the Attributes section. "
             "The result replaces the volume in the project (write it with File > SEG-Y export or the attribute writer); "
             "re-load the original file to go back.");
    ImGui::BeginDisabled(!expr[0] || S.job.valid());
    if (ImGui::Button(TR(ICON_CHECK " Apply to the volume"))) {
        const SeismicCube c = S.cube; const std::string e = expr; const int w = S.attrWin;
        geoLaunch(S, "Seismic calculator", [c, e, w]() -> std::function<void(App::GeoState&)> {
            SeismicCube o; std::string err;
            const bool ok = seismicCalculator(c, e, w, o, err);
            return [ok, o, err, e](App::GeoState& S) {
                if (!ok) { S.msg("Error: " + err); return; }
                S.cube = o; S.hasTh = false; S.msg("seismic calculator: " + e); ++S.cubeRev;
            };
        });
    }
    ImGui::EndDisabled();
    ImGui::SeparatorText(TR("Crop / decimate"));
    static int ir[2] = {0, 0}, jr[2] = {0, 0}, sr[2] = {0, 0}, dec[3] = {1, 1, 1};
    static std::string cubeKey;
    const std::string key = S.cube.name + std::to_string(S.cube.ni) + "x" + std::to_string(S.cube.nj) + "x" + std::to_string(S.cube.ns);
    if (key != cubeKey) { cubeKey = key; ir[0] = jr[0] = sr[0] = 0; ir[1] = S.cube.ni - 1; jr[1] = S.cube.nj - 1; sr[1] = S.cube.ns - 1; dec[0] = dec[1] = dec[2] = 1; }
    ImGui::SetNextItemWidth(fs * 12); ImGui::DragIntRange2("inline index", &ir[0], &ir[1], 0.3f, 0, S.cube.ni - 1);
    ImGui::SetNextItemWidth(fs * 12); ImGui::DragIntRange2("crossline index", &jr[0], &jr[1], 0.3f, 0, S.cube.nj - 1);
    ImGui::SetNextItemWidth(fs * 12); ImGui::DragIntRange2("sample", &sr[0], &sr[1], 0.5f, 0, S.cube.ns - 1);
    ImGui::SetNextItemWidth(fs * 12); ImGui::InputInt3("decimation (i, j, s)", dec);
    for (int& d : dec) d = std::clamp(d, 1, 16);
    const long long n = (long long)((ir[1] - ir[0]) / dec[0] + 1) * ((jr[1] - jr[0]) / dec[1] + 1) * ((sr[1] - sr[0]) / dec[2] + 1);
    ImGui::TextDisabled("result: %lld samples (%.1f MB), times %.1f .. %.1f ms", n, double(n) * 4 / 1048576.0, S.cube.t0 + sr[0] * S.cube.dt, S.cube.t0 + sr[1] * S.cube.dt);
    if (ImGui::Button(TR(ICON_SCISSORS " Crop the volume"))) {
        S.cube = cropCube(S.cube, ir[0], ir[1], jr[0], jr[1], sr[0], sr[1], dec[0], dec[1], dec[2]); ++S.cubeRev;
        S.hasTh = false; S.seeds.clear();
        S.msg("volume cropped to " + std::to_string(S.cube.ni) + " x " + std::to_string(S.cube.nj) + " x " + std::to_string(S.cube.ns));
    }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_SAVE " Save volume as SEG-Y..."))) {
        std::string p = S.cube.name + ".sgy";
        if (plat::fileDialog(true, "SEG-Y (*.sgy)\0*.sgy;*.segy\0", "sgy", p)) S.msg(writeSegy(p, S.cube) ? "written " + p : "Error: cannot write " + p);
    }
}

// ================================================================ P3 arbitrary intersection along a path (view 6)
namespace {
struct PathSample { double x, y, s; };
std::vector<PathSample> samplePath(const GPoly& P, int n) {
    std::vector<PathSample> out;
    const double L = polylineLength(P);
    if (P.pts.size() < 2 || L <= 0) return out;
    size_t seg = 0; double acc = 0;
    for (int q = 0; q < n; ++q) {
        const double s = L * q / double(n - 1);
        while (seg + 2 < P.pts.size()) {
            const double sl = std::hypot(P.pts[seg + 1][0] - P.pts[seg][0], P.pts[seg + 1][1] - P.pts[seg][1]);
            if (acc + sl >= s) break;
            acc += sl; ++seg;
        }
        const auto& a = P.pts[seg]; const auto& b = P.pts[seg + 1];
        const double sl = std::max(1e-12, std::hypot(b[0] - a[0], b[1] - a[1]));
        const double f = std::clamp((s - acc) / sl, 0.0, 1.0);
        out.push_back({a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), s});
    }
    return out;
}
// distance along the path of the closest point to (x, y) and the plan offset
void projectOnPath(const GPoly& P, double x, double y, double& s, double& d) {
    s = 0; d = 1e300; double acc = 0;
    for (size_t k = 0; k + 1 < P.pts.size(); ++k) {
        const auto& a = P.pts[k]; const auto& b = P.pts[k + 1];
        const double ux = b[0] - a[0], uy = b[1] - a[1], l2 = std::max(1e-12, ux * ux + uy * uy);
        const double f = std::clamp(((x - a[0]) * ux + (y - a[1]) * uy) / l2, 0.0, 1.0);
        const double dd = std::hypot(x - a[0] - f * ux, y - a[1] - f * uy);
        if (dd < d) { d = dd; s = acc + f * std::sqrt(l2); }
        acc += std::sqrt(l2);
    }
}
}  // namespace

void geoIntersectionView(App::GeoState& S, ImVec2 size) {
    const float fs = ImGui::GetFontSize();
    std::vector<std::string> pn; std::vector<int> idx;
    for (int p = 0; p < int(S.polys.size()); ++p) { idx.push_back(p); pn.push_back(S.polys[size_t(p)].name + (S.polys[size_t(p)].closed ? " (polygon)" : "")); }
    // C44: also along a deviated well (plan trajectory) or a 2-D seismic line (trace positions)
    const int np = int(S.polys.size()), nw = int(S.wells.size());
    for (int w = 0; w < nw; ++w) { idx.push_back(np + w); pn.push_back("well " + S.wells[size_t(w)].name); }
    for (int l = 0; l < int(S.lines.size()); ++l) { idx.push_back(np + nw + l); pn.push_back("2-D line " + S.lines[size_t(l)].name); }
    if (idx.empty()) {
        ImGui::TextDisabled("No intersection path yet.");
        if (ImGui::Button(TR(ICON_SCAN_LINE " Digitise a path on the map"))) { S.pickMode = 5; S.polyDraft = GPoly{}; S.view = 0; }
        ImGui::SameLine(); ImGui::TextDisabled("(or import polylines on the Data step > Polygons)");
        return;
    }
    if (S.isecPath < 0 || S.isecPath >= int(idx.size())) S.isecPath = idx[0];
    ImGui::SetNextItemWidth(fs * 12); vcombo("path", &S.isecPath, pn);
    static int mode = 0; static float wellTol = 0;
    static const char* modes[] = {"depth: grid + horizons + wells", "seismic (TWT)"};
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 13); ImGui::Combo("##imode", &mode, modes, S.hasCube ? 2 : 1);
    if (mode == 0 && S.hasGrid) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9); geoGridPropCombo("##igp", &S.gridProp, S); }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); ImGui::DragFloat("well distance", &wellTol, 5.0f, 0.0f, 1e6f, wellTol <= 0 ? "auto" : "%.0f");
    GPoly P;
    if (S.isecPath < np) P = S.polys[size_t(S.isecPath)];
    else if (S.isecPath < np + nw) {   // plan trajectory of the well (stations every ~1 % of MD), extended 15 % at both ends
        const auto& w = S.wells[size_t(S.isecPath - np)];
        P.name = w.name; P.closed = false;
        const double md1 = std::max(w.totalMd(), w.td);
        for (int q = 0; q <= 100; ++q) {
            const P3 a = w.at(md1 * q / 100);
            if (P.pts.empty() || std::hypot(a[0] - P.pts.back()[0], a[1] - P.pts.back()[1]) > 1e-6 * (1 + md1)) P.pts.push_back({a[0], a[1]});
        }
        if (P.pts.size() >= 2) {
            auto ext = [](const std::array<double, 2>& A, const std::array<double, 2>& B, double f) { return std::array<double, 2>{A[0] + (A[0] - B[0]) * f, A[1] + (A[1] - B[1]) * f}; };
            const double Lw = std::hypot(P.pts.back()[0] - P.pts.front()[0], P.pts.back()[1] - P.pts.front()[1]);
            const std::array<double, 2> A = P.pts.front(), B = P.pts[1], C = P.pts.back(), D = P.pts[P.pts.size() - 2];   // copies
            const double sab = std::hypot(A[0] - B[0], A[1] - B[1]), scd = std::hypot(C[0] - D[0], C[1] - D[1]);
            if (sab > 0 && scd > 0) { P.pts.insert(P.pts.begin(), ext(A, B, 0.15 * Lw / sab)); P.pts.push_back(ext(C, D, 0.15 * Lw / scd)); }
        } else { ImGui::TextDisabled("vertical well: no plan extent - digitise a path through it on the map"); return; }
    } else {
        const auto& ln = S.lines[size_t(S.isecPath - np - nw)];
        P.name = ln.name; P.closed = false;
        const int step = std::max(1, ln.ntr() / 200);
        for (int t = 0; t < ln.ntr(); t += step) P.pts.push_back({ln.x[size_t(t)], ln.y[size_t(t)]});
        if (ln.ntr() > 0 && (ln.ntr() - 1) % step) P.pts.push_back({ln.x.back(), ln.y.back()});
    }
    if (P.closed && P.pts.size() >= 3) P.pts.push_back(P.pts[0]);
    const double L = polylineLength(P);
    const int ns = std::clamp(int(size.x), 80, 1200);
    const auto smp = samplePath(P, ns);
    if (smp.empty()) { ImGui::TextDisabled("path has zero length"); return; }
    const double tol = wellTol > 0 ? wellTol : 0.03 * L;
    const ImVec2 psz(size.x, size.y - fs * 1.8f);
    if (mode == 1 && S.hasCube) {
        const SeismicCube& c = S.cube;
        std::vector<float> img(size_t(c.ns) * size_t(ns), std::nanf(""));
        float amax = 0;
        for (int q = 0; q < ns; ++q) {
            double fi, fj;
            if (!c.ij(smp[size_t(q)].x, smp[size_t(q)].y, fi, fj)) continue;
            const int i = int(std::lround(fi)), j = int(std::lround(fj));
            if (i < 0 || j < 0 || i >= c.ni || j >= c.nj) continue;
            for (int s = 0; s < c.ns; ++s) { const float a = c.at(i, j, s); img[size_t(s) * size_t(ns) + size_t(q)] = a; amax = std::max(amax, std::fabs(a)); }
        }
        if (amax <= 0) amax = 1;
        for (auto& v : img) if (!std::isfinite(v)) v = 0;
        if (ImPlot::BeginPlot("##isec_s", psz)) {
            ImPlot::SetupAxes("distance along the path", "TWT (ms)", 0, ImPlotAxisFlags_Invert);
            ImPlot::SetupAxesLimits(0, L, c.t0, c.t0 + (c.ns - 1) * c.dt, ImGuiCond_Always);
            ImPlot::PushColormap(ImPlotColormap_RdBu);
            ImPlot::PlotHeatmap("##amp", img.data(), c.ns, ns, -amax, amax, nullptr, ImPlotPoint(0, c.t0 + (c.ns - 1) * c.dt), ImPlotPoint(L, c.t0));
            ImPlot::PopColormap();
            if (S.hasTh && S.th.ni == c.ni && S.th.nj == c.nj) {
                std::vector<double> xs, ys;
                for (const auto& p : smp) {
                    double fi, fj; double t = NAN;
                    if (c.ij(p.x, p.y, fi, fj)) {
                        const int i = int(std::lround(fi)), j = int(std::lround(fj));
                        if (i >= 0 && j >= 0 && i < c.ni && j < c.nj) t = S.th.t[size_t(i) * size_t(c.nj) + size_t(j)];
                    }
                    xs.push_back(p.s); ys.push_back(t);
                }
                ImPlot::PlotLine(S.thName.c_str(), xs.data(), ys.data(), int(xs.size()), lspec(ImVec4(1, 0.9f, 0.2f, 1), 2.0f));
            }
            ImPlot::EndPlot();
        }
        ImGui::TextDisabled("path length %.0f, %d traces sampled (nearest trace), amplitude +-%.3g", L, ns, amax);
        return;
    }
    // depth domain
    std::vector<double> gp;
    const App::GeoState::PropTemplate* tpl = S.hasGrid ? geoTemplate(S, geoGridPropName(S, S.gridProp)) : nullptr;
    double lo = 1e300, hi = -1e300, zlo = 1e300, zhi = -1e300;
    if (S.hasGrid) {
        gp = geoGridProperty(S, S.gridProp);
        for (double v : gp) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
        if (lo > hi) { lo = 0; hi = 1; }
        geoTemplateRange(tpl, lo, hi);
    }
    std::vector<std::vector<double>> hzz(S.hz.size());
    for (size_t h = 0; h < S.hz.size(); ++h)
        for (const auto& p : smp) { const double z = S.hz[h].at(p.x, p.y); hzz[h].push_back(z); if (std::isfinite(z)) { zlo = std::min(zlo, z); zhi = std::max(zhi, z); } }
    struct Col { int i, j; double s0, s1; };
    std::vector<Col> cols;
    if (S.hasGrid) {
        const GeoGrid& g = S.g;
        const int nx = g.grid.nx, ny = g.grid.ny;
        for (int q = 0; q < ns; ++q) {
            double fi, fj; g.columnOf(smp[size_t(q)].x, smp[size_t(q)].y, fi, fj);
            const int i = int(std::floor(fi)), j = int(std::floor(fj));
            const double s0 = q ? 0.5 * (smp[size_t(q - 1)].s + smp[size_t(q)].s) : 0.0;
            const double s1 = q + 1 < ns ? 0.5 * (smp[size_t(q)].s + smp[size_t(q + 1)].s) : L;
            if (i < 0 || j < 0 || i >= nx || j >= ny) continue;
            if (!cols.empty() && cols.back().i == i && cols.back().j == j && std::fabs(cols.back().s1 - s0) < 1e-9) { cols.back().s1 = s1; continue; }
            cols.push_back({i, j, s0, s1});
            for (int k = 0; k < g.grid.nz; ++k) {
                const size_t c = size_t(i + nx * (j + ny * k));
                if (std::isfinite(g.top[c])) zlo = std::min(zlo, g.top[c]);
                if (std::isfinite(g.bot[c])) zhi = std::max(zhi, g.bot[c]);
            }
        }
    }
    // wells isNear the path
    struct WP { std::string name; std::vector<double> s, z; std::vector<std::pair<double, double>> tops; std::vector<std::string> topN; };
    std::vector<WP> wps;
    for (const auto& w : S.wells) {
        WP wp; wp.name = w.name;
        const double md1 = std::max(w.totalMd(), w.td);
        const int nst = 60;
        bool isNear = false;
        for (int q = 0; q <= nst; ++q) {
            const P3 p = w.at(md1 * q / nst);
            double s, d; projectOnPath(P, p[0], p[1], s, d);
            if (d <= tol) isNear = true;
            wp.s.push_back(s); wp.z.push_back(p[2]);
        }
        if (!isNear) continue;
        for (const auto& t : w.tops) { const P3 p = w.at(t.md); double s, d; projectOnPath(P, p[0], p[1], s, d); wp.tops.push_back({s, p[2]}); wp.topN.push_back(t.surface); }
        wps.push_back(std::move(wp));
    }
    if (zlo > zhi)   // no horizons / grid: the well paths set the depth range (otherwise the interval of interest is framed)
        for (const auto& wp : wps) for (double z : wp.z) { zlo = std::min(zlo, z); zhi = std::max(zhi, z); }
    if (zlo > zhi) { zlo = 0; zhi = 1000; }
    const double pad = 0.25 * (zhi - zlo) + 1;
    if (ImPlot::BeginPlot("##isec", psz)) {
        ImPlot::SetupAxes("distance along the path", "depth", 0, ImPlotAxisFlags_Invert);
        ImPlot::SetupAxesLimits(0, L, zlo - pad, zhi + pad, ImGuiCond_Appearing);
        if (S.hasGrid) {
            const GeoGrid& g = S.g;
            const int nx = g.grid.nx, ny = g.grid.ny;
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            for (const auto& cc : cols)
                for (int k = 0; k < g.grid.nz; ++k) {
                    const size_t c = size_t(cc.i + nx * (cc.j + ny * k));
                    if (c >= gp.size() || !std::isfinite(gp[c]) || !std::isfinite(g.top[c]) || !std::isfinite(g.bot[c]) || !geoCellVisible(S, c)) continue;
                    const double tv = geoTemplateValue(tpl, gp[c]);
                    if (!std::isfinite(tv)) continue;
                    dl->AddRectFilled(ImPlot::PlotToPixels(cc.s0, g.top[c]), ImPlot::PlotToPixels(cc.s1, g.bot[c]), geoTemplateColor(tpl, tv, lo, hi));
                }
            ImPlot::PopPlotClipRect();
        }
        static const ImVec4 hc[] = {{1, 0.85f, 0.2f, 1}, {1, 0.4f, 0.3f, 1}, {0.4f, 0.9f, 1, 1}, {0.9f, 0.5f, 1, 1}, {0.6f, 1, 0.5f, 1}, {1, 1, 1, 1}};
        std::vector<double> ss; for (const auto& p : smp) ss.push_back(p.s);
        for (size_t h = 0; h < S.hz.size(); ++h) ImPlot::PlotLine(S.hz[h].name.c_str(), ss.data(), hzz[h].data(), ns, lspec(hc[h % 6], 2.0f));
        for (const auto& f : S.faults) {   // C44: fault cuts (stick interpolation per depth level)
            const auto br = faultCutOnPath(f, P, zlo - pad, zhi + pad, 60);
            for (size_t b = 0; b < br.size(); ++b) {
                std::vector<double> xs, ys; for (const auto& q : br[b]) { xs.push_back(q[0]); ys.push_back(q[1]); }
                ImPlot::PlotLine(b ? ("##f" + f.name + std::to_string(b)).c_str() : f.name.c_str(), xs.data(), ys.data(), int(xs.size()), lspec(ImVec4(1, 0.25f, 0.25f, 1), 2.5f));
                if (!xs.empty()) { ImPlotSpec ts; ts.Flags = ImPlotItemFlags_NoFit; ImPlot::PlotText(f.name.c_str(), xs.front(), ys.front(), ImVec2(0, -fs * 0.7f), ts); }
            }
        }
        for (const auto& wp : wps) {
            ImPlot::PlotLine(("##w" + wp.name).c_str(), wp.s.data(), wp.z.data(), int(wp.s.size()), lspec(ImVec4(0.95f, 0.95f, 0.95f, 1), 2.0f));
            ImPlot::PlotText(wp.name.c_str(), wp.s.front(), wp.z.front(), ImVec2(0, -10));
            for (size_t t = 0; t < wp.tops.size(); ++t) {
                const double xs = wp.tops[t].first, ys = wp.tops[t].second;
                ImPlot::PlotScatter(("##t" + wp.name + std::to_string(t)).c_str(), &xs, &ys, 1, mspec(ImVec4(1, 0.3f, 0.3f, 1), 4));
                ImPlot::PlotText(wp.topN[t].c_str(), xs, ys, ImVec2(fs * 2, 0));
            }
        }
        ImPlot::EndPlot();
    }
    if (S.hasGrid) ImGui::TextDisabled("path length %.0f; %s %.4g .. %.4g%s (template %s); %zu well(s) within %.0f", L, geoGridPropName(S, S.gridProp).c_str(),
                                       tpl && tpl->logScale ? std::pow(10.0, lo) : lo, tpl && tpl->logScale ? std::pow(10.0, hi) : hi, tpl && tpl->logScale ? " log" : "",
                                       tpl ? tpl->name.c_str() : "-", wps.size(), tol);
    else ImGui::TextDisabled("path length %.0f; build the grid to see properties; %zu well(s) within %.0f", L, wps.size(), tol);
}

// ================================================================ P2 property templates
#include "colormaps.h"
namespace {
std::string upper(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
void defaultTemplates(App::GeoState& S) {
    using T = App::GeoState::PropTemplate;
    auto mk = [](const char* n, const char* m, const char* u, Colormap c, bool rev, bool fixed, double lo, double hi, bool lg) {
        T t; t.name = n; t.match = m; t.unit = u; t.cmap = int(c); t.reverse = rev; t.fixed = fixed; t.lo = lo; t.hi = hi; t.logScale = lg; return t;
    };
    S.templates = {
        mk("Porosity", "PORO|PHI", "v/v", Colormap::Viridis, false, true, 0.0, 0.35, false),
        mk("Permeability", "PERM|KX|KY|KZ", "mD", Colormap::Turbo, false, false, 0.1, 5000, true),
        mk("Water saturation", "SW|SAT", "v/v", Colormap::RdBu, true, true, 0.0, 1.0, false),
        mk("Net-to-gross", "NTG|NET", "v/v", Colormap::Cividis, false, true, 0.0, 1.0, false),
        mk("Facies / zones", "FACIES|ZONE|BLOCK|SEGMENT", "code", Colormap::Spectral, false, false, 0, 1, false),
        mk("Thickness", "ISOCH|THICK|DZ", "length", Colormap::Plasma, false, false, 0, 1, false),
        mk("Velocity", "VINT|VEL", "length/s", Colormap::Jet, false, false, 0, 1, false),
        mk("Depth / time (default)", "*", "length", Colormap::Viridis, true, false, 0, 1, false),
    };
    S.templatesInit = true;
}
}  // namespace

const App::GeoState::PropTemplate* geoTemplate(App::GeoState& S, const std::string& propName) {
    if (!S.templatesInit) defaultTemplates(S);
    const std::string u = upper(propName);
    for (const auto& t : S.templates) {
        if (t.match == "*") return &t;
        size_t p = 0;
        const std::string m = upper(t.match);
        while (p <= m.size()) {
            size_t e = m.find('|', p); if (e == std::string::npos) e = m.size();
            const std::string tok = m.substr(p, e - p);
            if (!tok.empty() && u.find(tok) != std::string::npos) return &t;
            p = e + 1;
        }
    }
    return S.templates.empty() ? nullptr : &S.templates.back();
}
void geoTemplateRange(const App::GeoState::PropTemplate* t, double& lo, double& hi) {
    if (t && t->fixed) { lo = t->lo; hi = t->hi; }
    if (t && t->logScale) {
        if (!t->fixed) { lo = std::max(lo, 1e-6); hi = std::max(hi, lo * 10); }
        lo = std::log10(std::max(lo, 1e-12)); hi = std::log10(std::max(hi, 1e-12));
    }
    if (!(hi > lo)) hi = lo + 1;
}
double geoTemplateValue(const App::GeoState::PropTemplate* t, double v) {
    return (t && t->logScale) ? (v > 0 ? std::log10(v) : NAN) : v;
}
int geoPlotColormap(int cmap, bool reverse) {
    cmap = std::clamp(cmap, 0, int(Colormap::Count) - 1);
    const std::string nm = std::string("rs_") + colormapName(Colormap(cmap)) + (reverse ? "_r" : "");
    int id = ImPlot::GetColormapIndex(nm.c_str());
    if (id < 0) {
        ImVec4 c[64];
        for (int q = 0; q < 64; ++q) { float rgb[3]; colormapRGB(Colormap(cmap), q / 63.0f, rgb, reverse); c[q] = ImVec4(rgb[0], rgb[1], rgb[2], 1); }
        id = ImPlot::AddColormap(nm.c_str(), c, 64, false);
    }
    return id;
}
ImU32 geoTemplateColor(const App::GeoState::PropTemplate* t, double v, double lo, double hi) {
    float rgb[3];
    const float f = float(std::clamp((v - lo) / std::max(1e-300, hi - lo), 0.0, 1.0));
    colormapRGB(t ? Colormap(std::clamp(t->cmap, 0, int(Colormap::Count) - 1)) : Colormap::Viridis, f, rgb, t && t->reverse);
    return IM_COL32(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255), 255);
}

void geoTemplatesWindow(App::GeoState& S) {
    if (!S.showTemplates) return;
    if (!S.templatesInit) defaultTemplates(S);
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextWindowSize(ImVec2(fs * 46, fs * 20), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin(ICON_PALETTE " Property templates", &S.showTemplates)) { ImGui::End(); return; }
    ImGui::TextWrapped("A property (grid property, horizon, calculated surface) uses the first template whose match list "
                       "(| separated, case-insensitive) occurs in its name; * matches everything. Map, grid cross-section, "
                       "intersection and 3-D scene share these colour maps and ranges.");
    std::vector<std::string> cm;
    for (int c = 0; c < int(Colormap::Count); ++c) cm.push_back(colormapName(Colormap(c)));
    if (ImGui::BeginTable("##tpl", 9, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp)) {
        for (const char* h : {"name", "match", "unit", "colour map", "rev", "fixed", "min", "max", "log"}) ImGui::TableSetupColumn(h);
        ImGui::TableHeadersRow();
        int del = -1;
        for (int r = 0; r < int(S.templates.size()); ++r) {
            auto& t = S.templates[size_t(r)];
            ImGui::PushID(r);
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); { char b[64]; std::snprintf(b, sizeof b, "%s", t.name.c_str()); if (ImGui::InputText("##n", b, sizeof b)) t.name = b; }
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); { char b[96]; std::snprintf(b, sizeof b, "%s", t.match.c_str()); if (ImGui::InputText("##m", b, sizeof b)) t.match = b; }
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); { char b[32]; std::snprintf(b, sizeof b, "%s", t.unit.c_str()); if (ImGui::InputText("##u", b, sizeof b)) t.unit = b; }
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); vcombo("##c", &t.cmap, cm);
            ImGui::TableNextColumn(); ImGui::Checkbox("##r", &t.reverse);
            ImGui::TableNextColumn(); ImGui::Checkbox("##f", &t.fixed);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::BeginDisabled(!t.fixed); ImGui::InputDouble("##lo", &t.lo, 0, 0, "%.4g"); ImGui::EndDisabled();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::BeginDisabled(!t.fixed); ImGui::InputDouble("##hi", &t.hi, 0, 0, "%.4g"); ImGui::EndDisabled();
            ImGui::TableNextColumn(); ImGui::Checkbox("##l", &t.logScale);
            ImGui::SameLine(); if (t.match != "*" && ImGui::SmallButton(ICON_X)) del = r;
            ImGui::PopID();
        }
        ImGui::EndTable();
        if (del >= 0) S.templates.erase(S.templates.begin() + del);
    }
    if (ImGui::Button(TR(ICON_PLUS " Template"))) { App::GeoState::PropTemplate t; t.name = "New"; t.match = "NAME"; S.templates.insert(S.templates.end() - 1, t); }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_ROTATE_CCW " Defaults"))) defaultTemplates(S);
    ImGui::End();
}

// ================================================================ P7 shared cell filter
bool geoCellVisible(const App::GeoState& S, size_t c) {
    if (!S.cf.on || !S.hasGrid) return true;
    const GeoGrid& g = S.g;
    const size_t nxy = size_t(g.grid.nx) * size_t(g.grid.ny);
    if (nxy == 0) return true;
    const size_t k = c / nxy, col = c % nxy;
    if (k < g.layerZone.size()) { const int z = g.layerZone[k]; if (z >= 0 && size_t(z) < S.cf.zoneOff.size() && S.cf.zoneOff[size_t(z)]) return false; }
    if (col < g.block.size()) { const int b = g.block[col]; if (b >= 0 && size_t(b) < S.cf.blockOff.size() && S.cf.blockOff[size_t(b)]) return false; }
    if (S.hasFacies && c < S.fr.facies.size() && !S.cf.faciesOff.empty()) {
        const int code = S.fr.facies[c];
        for (size_t q = 0; q < S.fr.codes.size() && q < S.cf.faciesOff.size(); ++q) if (S.fr.codes[q] == code) return !S.cf.faciesOff[q];
    }
    return true;
}
void geoFilterButton(App::GeoState& S) {
    const bool act = S.cf.on;
    if (act) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::ColorConvertU32ToFloat4(Pal::Accent));
    if (ImGui::Button(TR(ICON_FILTER " Filter"))) ImGui::OpenPopup("##geofilter");
    if (act) ImGui::PopStyleColor();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Zones / fault blocks / facies shown in every window (map, cross-sections, intersection, 3-D) and used by the volumetrics");
    ImGui::SameLine();
    if (ImGui::Button(ICON_PALETTE)) S.showTemplates = !S.showTemplates;
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Property templates (colour map, range, log scale per property kind)");
    if (ImGui::BeginPopup("##geofilter")) {
        ImGui::Checkbox(TR("filter on"), &S.cf.on);
        if (!S.hasGrid) ImGui::TextDisabled("build the grid first");
        else {
            ImGui::BeginDisabled(!S.cf.on);
            const int nz = int(S.gs.zones.size());
            if (nz > 0) {
                ImGui::SeparatorText(TR("Zones"));
                S.cf.zoneOff.resize(size_t(nz), 0);
                for (int z = 0; z < nz; ++z) { bool on = !S.cf.zoneOff[size_t(z)]; if (ImGui::Checkbox((S.gs.zones[size_t(z)].name + "##z" + std::to_string(z)).c_str(), &on)) S.cf.zoneOff[size_t(z)] = !on; }
            }
            if (S.g.nBlocks > 1) {
                ImGui::SeparatorText(TR("Fault blocks"));
                S.cf.blockOff.resize(size_t(S.g.nBlocks), 0);
                for (int b = 0; b < S.g.nBlocks; ++b) { bool on = !S.cf.blockOff[size_t(b)]; if (ImGui::Checkbox(("block " + std::to_string(b + 1)).c_str(), &on)) S.cf.blockOff[size_t(b)] = !on; }
            }
            if (S.hasFacies && !S.fr.codes.empty()) {
                ImGui::SeparatorText(TR("Facies"));
                S.cf.faciesOff.resize(S.fr.codes.size(), 0);
                for (size_t q = 0; q < S.fr.codes.size(); ++q) { bool on = !S.cf.faciesOff[q]; if (ImGui::Checkbox(("facies " + std::to_string(S.fr.codes[q])).c_str(), &on)) S.cf.faciesOff[q] = !on; }
            }
            if (ImGui::SmallButton(TR("show all"))) { std::fill(S.cf.zoneOff.begin(), S.cf.zoneOff.end(), 0); std::fill(S.cf.blockOff.begin(), S.cf.blockOff.end(), 0); std::fill(S.cf.faciesOff.begin(), S.cf.faciesOff.end(), 0); }
            ImGui::EndDisabled();
        }
        ImGui::EndPopup();
    }
}

// ================================================================ P5 well section (templated tracks, flattening, top picking)
void geoWellSection(App::GeoState& S, ImVec2 size) {
    const float fs = ImGui::GetFontSize();
    std::vector<std::string> curves;
    for (const auto& w : S.wells) for (const auto& l : w.logs) for (const auto& c : l.curves) if (std::find(curves.begin(), curves.end(), c.mnem) == curves.end()) curves.push_back(c.mnem);
    if (S.wells.empty()) { ImGui::TextDisabled("No wells loaded."); return; }
    if (S.wsTracks.empty()) {   // default template: up to three curves, permeability on a log scale
        for (const auto& c : curves) {
            if (S.wsTracks.size() >= 3) break;
            App::GeoState::WsTrack t; t.curve = c; t.logX = upper(c).find("PERM") != std::string::npos || upper(c) == "K"; S.wsTracks.push_back(t);
        }
    }
    const auto tops = topNames(S);
    // toolbar: tracks, wells, flatten, pick
    if (ImGui::Button(TR(ICON_LIST " Tracks"))) ImGui::OpenPopup("##wstracks");
    if (ImGui::BeginPopup("##wstracks")) {
        int del = -1;
        for (int t = 0; t < int(S.wsTracks.size()); ++t) {
            ImGui::PushID(t);
            int ci = -1; for (size_t q = 0; q < curves.size(); ++q) if (curves[q] == S.wsTracks[size_t(t)].curve) ci = int(q);
            ImGui::SetNextItemWidth(fs * 8); if (vcombo("##c", &ci, curves) && ci >= 0) S.wsTracks[size_t(t)].curve = curves[size_t(ci)];
            ImGui::SameLine(); ImGui::Checkbox(TR("log"), &S.wsTracks[size_t(t)].logX);
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4); ImGui::DragFloat("##w", &S.wsTracks[size_t(t)].width, 0.02f, 0.4f, 4.0f, "width %.1f");
            ImGui::SameLine(); ImGui::Checkbox(TR("template"), &S.wsTracks[size_t(t)].useTpl);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("use the property template of this curve (log scale, fixed range) when one matches");
            ImGui::SameLine(); if (t > 0 && ImGui::SmallButton(TR("<"))) std::swap(S.wsTracks[size_t(t)], S.wsTracks[size_t(t - 1)]);
            ImGui::SameLine(); if (ImGui::SmallButton(ICON_X)) del = t;
            ImGui::PopID();
        }
        if (del >= 0) S.wsTracks.erase(S.wsTracks.begin() + del);
        if (!curves.empty() && ImGui::SmallButton(TR(ICON_PLUS " track"))) { App::GeoState::WsTrack t; t.curve = curves[0]; S.wsTracks.push_back(t); }
        ImGui::TextDisabled("discrete curves (FACIES ...) are drawn as coloured blocks");
        ImGui::EndPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_EYE " Wells"))) ImGui::OpenPopup("##wswells");
    // C44 well order: names in wsOrder first, the rest in project order
    std::vector<int> order;
    for (const auto& n : S.wsOrder) for (int w = 0; w < int(S.wells.size()); ++w) if (S.wells[size_t(w)].name == n && std::find(order.begin(), order.end(), w) == order.end()) order.push_back(w);
    for (int w = 0; w < int(S.wells.size()); ++w) if (std::find(order.begin(), order.end(), w) == order.end()) order.push_back(w);
    if (ImGui::BeginPopup("##wswells")) {
        int up = -1;
        for (int q = 0; q < int(order.size()); ++q) {
            const auto& w = S.wells[size_t(order[size_t(q)])];
            ImGui::PushID(q);
            const auto it = std::find(S.wsHidden.begin(), S.wsHidden.end(), w.name);
            bool on = it == S.wsHidden.end();
            if (ImGui::Checkbox(w.name.c_str(), &on)) { if (on) S.wsHidden.erase(it); else S.wsHidden.push_back(w.name); }
            if (q > 0) { ImGui::SameLine(fs * 8); if (ImGui::SmallButton(TR("move up"))) up = q; }
            ImGui::PopID();
        }
        if (up > 0) { std::swap(order[size_t(up)], order[size_t(up - 1)]); S.wsOrder.clear(); for (int w : order) S.wsOrder.push_back(S.wells[size_t(w)].name); }
        if (ImGui::SmallButton(TR("project order"))) S.wsOrder.clear();
        ImGui::EndPopup();
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8); strCombo(TR("flatten on"), &S.wsFlat, tops);
    ImGui::SameLine();
    if (ImGui::Checkbox(TR("pick top"), &S.wsPick) && S.wsPick && S.wsPickTop.empty() && !tops.empty()) S.wsPickTop = tops[0];
    if (S.wsPick) {
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8);
        char b[64]; std::snprintf(b, sizeof b, "%s", S.wsPickTop.c_str());
        if (ImGui::InputText("##ptop", b, sizeof b)) S.wsPickTop = b;
        if (ImGui::BeginPopupContextItem("##ptopc")) { for (const auto& t : tops) if (ImGui::Selectable(t.c_str())) S.wsPickTop = t; ImGui::EndPopup(); }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("top name (right-click: existing tops); click in a well's track sets / moves this top, right-click a top line deletes it");
    }
    std::vector<int> wi;
    for (int w : order) if (std::find(S.wsHidden.begin(), S.wsHidden.end(), S.wells[size_t(w)].name) == S.wsHidden.end()) wi.push_back(w);
    const int nt = std::max(1, int(S.wsTracks.size()));
    if (wi.empty()) { ImGui::TextDisabled("all wells hidden"); return; }
    static const ImVec4 pal[] = {{0.95f, 0.75f, 0.2f, 1}, {0.35f, 0.75f, 1, 1}, {0.4f, 0.9f, 0.45f, 1}, {1, 0.45f, 0.4f, 1}, {0.8f, 0.55f, 1, 1}, {0.95f, 0.95f, 0.95f, 1}, {0.3f, 0.95f, 0.9f, 1}, {1, 0.6f, 0.85f, 1}};
    auto topColor = [&](const std::string& n) { for (size_t q = 0; q < tops.size(); ++q) if (tops[q] == n) return pal[q % 8]; return pal[5]; };
    // depth window: the logs / tops of the shown wells (relative to the flattening top); refitted when the selection changes
    double zlo = 1e300, zhi = -1e300;
    for (int w : wi) {
        const GeoWell& W = S.wells[size_t(w)];
        double zf = 0; bool ok = S.wsFlat.empty();
        for (const auto& t : W.tops) if (!S.wsFlat.empty() && t.surface == S.wsFlat) { zf = W.at(t.md)[2]; ok = true; }
        if (!ok) continue;
        for (const auto& l : W.logs) if (!l.md.empty()) { zlo = std::min(zlo, W.at(l.md.front())[2] - zf); zhi = std::max(zhi, W.at(l.md.back())[2] - zf); }
        for (const auto& t : W.tops) { const double z = W.at(t.md)[2] - zf; zlo = std::min(zlo, z); zhi = std::max(zhi, z); }
    }
    if (zlo > zhi) { zlo = 0; zhi = 1000; }
    static std::string fitKey;
    std::string key = S.wsFlat + "|" + std::to_string(S.wells.size());
    for (int w : wi) key += "," + std::to_string(w);
    const bool refit = key != fitKey; fitKey = key;
    const double zpad = 0.03 * (zhi - zlo) + 1;
    const ImVec2 psz(size.x, size.y - fs * 0.5f);
    const bool narrow = psz.x / float(std::max<size_t>(1, wi.size() * size_t(nt))) < fs * 7;   // x tick labels would overlap
    int pickWell = -1; double pickZ = NAN; int delWell = -1; std::string delTop;
    std::vector<float> colRatio;   // the first track of each well also carries the depth axis
    float wsum = 0; for (int t = 0; t < nt; ++t) wsum += t < int(S.wsTracks.size()) ? S.wsTracks[size_t(t)].width : 1.0f;
    for (size_t w = 0; w < wi.size(); ++w)
        for (int t = 0; t < nt; ++t) {
            const float tw = t < int(S.wsTracks.size()) ? S.wsTracks[size_t(t)].width : 1.0f;
            colRatio.push_back(tw + (t == 0 ? fs * 3.2f * wsum / nt / std::max(1.0f, psz.x / float(wi.size() * size_t(nt))) : 0.0f));
        }
    if (ImPlot::BeginSubplots("##wsec", 1, int(wi.size()) * nt, psz, ImPlotSubplotFlags_LinkRows | ImPlotSubplotFlags_NoResize, nullptr, colRatio.data())) {
        for (int w : wi) {
            const GeoWell& W = S.wells[size_t(w)];
            double zf = 0;
            if (!S.wsFlat.empty()) {
                bool found = false;
                for (const auto& t : W.tops) if (t.surface == S.wsFlat) { zf = W.at(t.md)[2]; found = true; }
                if (!found) zf = NAN;
            }
            for (int t = 0; t < nt; ++t) {
                const std::string title = (t == 0 ? W.name : std::string()) + "##" + W.name + std::to_string(t);
                if (!ImPlot::BeginPlot(title.c_str(), ImVec2(-1, -1), ImPlotFlags_NoLegend)) continue;
                const App::GeoState::WsTrack* T = t < int(S.wsTracks.size()) ? &S.wsTracks[size_t(t)] : nullptr;
                // C44: the curve's property template (not the catch-all default) sets log scale and a fixed range
                const App::GeoState::PropTemplate* tp = T && T->useTpl ? geoTemplate(S, T->curve) : nullptr;
                if (tp && tp->match == "*") tp = nullptr;
                const bool logX = T && (T->logX || (tp && tp->logScale));
                const bool fixedX = tp && tp->fixed && tp->hi > tp->lo && (!logX || tp->lo > 0);
                ImPlot::SetupAxes(T ? T->curve.c_str() : "", t == 0 ? (S.wsFlat.empty() ? "TVD SS" : ("TVD - " + S.wsFlat).c_str()) : nullptr,
                                  (fixedX ? 0 : ImPlotAxisFlags_AutoFit) | (narrow ? ImPlotAxisFlags_NoTickLabels : 0), ImPlotAxisFlags_Invert | (t ? ImPlotAxisFlags_NoTickLabels : 0));
                if (logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                if (fixedX) ImPlot::SetupAxisLimits(ImAxis_X1, tp->lo, tp->hi, ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_Y1, zlo - zpad, zhi + zpad, refit ? ImPlotCond_Always : ImPlotCond_Once);
                ImPlot::SetupFinish();
                if (!std::isfinite(zf)) { ImPlot::PlotText(("no " + S.wsFlat).c_str(), 0.5, 0.5); ImPlot::EndPlot(); continue; }
                if (T)
                    for (const auto& l : W.logs) {
                        const LogCurve* c = l.curve(T->curve);
                        if (!c) continue;
                        const bool disc = c->discrete || upper(c->mnem).find("FACIES") != std::string::npos;
                        if (disc) {
                            ImPlot::PushPlotClipRect();
                            ImDrawList* dl = ImPlot::GetPlotDrawList();
                            const ImPlotRect lim = ImPlot::GetPlotLimits();
                            for (size_t k = 0; k + 1 < l.md.size() && k + 1 < c->v.size(); ++k) {
                                const double v = c->v[k];
                                if (!std::isfinite(v) || v == l.nullValue) continue;
                                const double z0 = W.at(l.md[k])[2] - zf, z1 = W.at(l.md[k + 1])[2] - zf;
                                const ImVec4 col = pal[size_t(std::abs(int(std::lround(v)))) % 8];
                                dl->AddRectFilled(ImPlot::PlotToPixels(lim.X.Min, z0), ImPlot::PlotToPixels(lim.X.Max, z1), ImGui::ColorConvertFloat4ToU32(col));
                            }
                            ImPlot::PopPlotClipRect();
                        } else {
                            std::vector<double> xs, ys;
                            for (size_t k = 0; k < l.md.size() && k < c->v.size(); ++k) {
                                const double v = c->v[k];
                                if (!std::isfinite(v) || v == l.nullValue || (logX && v <= 0)) continue;
                                xs.push_back(v); ys.push_back(W.at(l.md[k])[2] - zf);
                            }
                            ImPlot::PlotLine("##c", xs.data(), ys.data(), int(xs.size()), lspec(pal[size_t(t) % 8], 1.3f));
                        }
                    }
                for (const auto& tp : W.tops) {
                    const double z = W.at(tp.md)[2] - zf;
                    ImPlot::PlotInfLines(("##t" + tp.surface).c_str(), &z, 1, [&] { ImPlotSpec s = lspec(topColor(tp.surface), 1.8f); s.Flags = ImPlotInfLinesFlags_Horizontal; return s; }());
                    if (t == 0) {   // label at the left edge; NoFit: must not feed the x auto-fit
                        const ImPlotRect lim = ImPlot::GetPlotLimits(); ImPlotSpec ts; ts.Flags = ImPlotItemFlags_NoFit;
                        ImPlot::PlotText(tp.surface.c_str(), lim.X.Min, z, ImVec2(fs * 2.2f, -fs * 0.6f), ts);
                    }
                }
                if (ImPlot::IsPlotHovered()) {
                    const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                    const double tvd = mp.y + zf;
                    const double mdh = W.mdAtTvd(tvd);
                    std::string vals;
                    for (const auto& tr : S.wsTracks) for (const auto& l : W.logs) if (const LogCurve* c = l.curve(tr.curve)) { char vb[64]; std::snprintf(vb, sizeof vb, "\n%s %.4g", tr.curve.c_str(), l.valueAt(*c, mdh)); vals += vb; break; }
                    ImGui::SetTooltip("%s  TVD %.2f  MD %.2f%s%s", W.name.c_str(), tvd, mdh, vals.c_str(), S.wsPick ? ("\nclick: set " + S.wsPickTop).c_str() : "");
                    if (S.wsPick && ImGui::IsMouseReleased(ImGuiMouseButton_Left) && ImGui::GetIO().MouseDragMaxDistanceSqr[0] < 16.0f) { pickWell = w; pickZ = tvd; }
                    if (S.wsPick && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                        double best = 1e300; std::string bn;
                        for (const auto& tp : W.tops) { const double d = std::fabs(W.at(tp.md)[2] - tvd); if (d < best) { best = d; bn = tp.surface; } }
                        const ImPlotRect lim = ImPlot::GetPlotLimits();
                        if (!bn.empty() && best < 0.02 * (lim.Y.Max - lim.Y.Min)) { delWell = w; delTop = bn; }
                    }
                }
                ImPlot::EndPlot();
            }
        }
        ImPlot::EndSubplots();
    }
    if (pickWell >= 0 && !S.wsPickTop.empty()) {
        GeoWell& W = S.wells[size_t(pickWell)];
        const double md = W.mdAtTvd(pickZ);
        if (std::isfinite(md)) {
            bool rep = false;
            for (auto& tp : W.tops) if (tp.surface == S.wsPickTop) { tp.md = md; rep = true; }
            if (!rep) { Top tp; tp.surface = S.wsPickTop; tp.md = md; W.tops.push_back(tp); std::sort(W.tops.begin(), W.tops.end(), [](const Top& a, const Top& b) { return a.md < b.md; }); }
            char b[160]; std::snprintf(b, sizeof b, "%s: top %s %s at MD %.2f (TVD %.2f)%s", W.name.c_str(), S.wsPickTop.c_str(), rep ? "moved" : "set", md, pickZ, S.hasGrid ? " - rebuild the grid to tie" : "");
            S.msg(b);
        }
    }
    if (delWell >= 0) {
        auto& T = S.wells[size_t(delWell)].tops;
        T.erase(std::remove_if(T.begin(), T.end(), [&](const Top& t) { return t.surface == delTop; }), T.end());
        S.msg(S.wells[size_t(delWell)].name + ": top " + delTop + " deleted");
    }
}

// ================================================================ P1 object tree (left column)
namespace {
void rangeOf(const std::vector<double>& v, double& lo, double& hi, double& mean, int& n) {
    lo = 1e300; hi = -1e300; mean = 0; n = 0;
    for (double x : v) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); mean += x; ++n; }
    if (n) mean /= n;
}
bool renameField(std::string& name) {
    char b[96]; std::snprintf(b, sizeof b, "%s", name.c_str());
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
    if (ImGui::InputText("name", b, sizeof b, ImGuiInputTextFlags_EnterReturnsTrue) && b[0]) { name = b; return true; }
    return false;
}

// C46 object-tree folders: drag source, drop targets, "Move to folder", folder nodes per kind
using ressim::geo::ObjKind;
struct DragObj { int kind, index; };
std::string objNameOf(const App::GeoState& S, int kind, int i) {
    switch (ObjKind(kind)) {
    case ObjKind::Well: return i >= 0 && i < int(S.wells.size()) ? S.wells[size_t(i)].name : std::string();
    case ObjKind::Surface: return i >= 0 && i < int(S.hz.size()) ? S.hz[size_t(i)].name : std::string();
    case ObjKind::Polygon: return i >= 0 && i < int(S.polys.size()) ? S.polys[size_t(i)].name : std::string();
    case ObjKind::Fault: return i >= 0 && i < int(S.faults.size()) ? S.faults[size_t(i)].name : std::string();
    case ObjKind::Line: return i >= 0 && i < int(S.lines.size()) ? S.lines[size_t(i)].name : std::string();
    }
    return std::string();
}
void dragSource(ObjKind k, int i, const std::string& name) {
    if (ImGui::BeginDragDropSource()) {
        const DragObj d{int(k), i};
        ImGui::SetDragDropPayload("RSGEO_OBJ", &d, sizeof d);
        ImGui::Text("%s %s", ICON_FOLDER, name.c_str());
        ImGui::EndDragDropSource();
    }
}
void dropTarget(App::GeoState& S, ObjKind k, const std::string& folder) {
    if (!ImGui::BeginDragDropTarget()) return;
    if (const ImGuiPayload* pl = ImGui::AcceptDragDropPayload("RSGEO_OBJ")) {
        DragObj d; std::memcpy(&d, pl->Data, sizeof d);
        const std::string n = objNameOf(S, d.kind, d.index);
        if (d.kind == int(k) && !n.empty()) S.folders.assign(k, n, folder);
    }
    ImGui::EndDragDropTarget();
}
bool newFolderField(const char* id, std::string& out) {
    static char nb[64] = "";
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
    if (ImGui::InputTextWithHint(id, TR("new folder + Enter"), nb, sizeof nb, ImGuiInputTextFlags_EnterReturnsTrue) && nb[0]) { out = nb; nb[0] = 0; return true; }
    return false;
}
void moveToMenu(App::GeoState& S, ObjKind k, const std::string& name) {
    if (!ImGui::BeginMenu(TR(ICON_FOLDER " Move to folder"))) return;
    const std::string cur = S.folders.folderOf(k, name);
    if (ImGui::MenuItem(TR("(top level)"), nullptr, cur.empty())) S.folders.assign(k, name, "");
    for (const auto& f : S.folders.foldersOf(k)) if (ImGui::MenuItem(f.c_str(), nullptr, cur == f)) S.folders.assign(k, name, f);
    ImGui::Separator();
    std::string nf;
    if (newFolderField("##mvnf", nf)) { S.folders.assign(k, name, nf); ImGui::CloseCurrentPopup(); }
    ImGui::EndMenu();
}
// the kind node (Wells, Surfaces ...): drop target = back to the top level, right-click = new folder
void kindNodeExtras(App::GeoState& S, ObjKind k) {
    dropTarget(S, k, "");
    if (ImGui::BeginPopupContextItem("##kindc")) {
        std::string nf;
        ImGui::TextDisabled("%s", TR("Folder"));
        if (newFolderField("##knf", nf)) { S.folders.addFolder(k, nf); ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}
template <class NameOf, class Item>
void kindTree(App::GeoState& S, ObjKind k, int n, NameOf nameOf, Item item) {
    const std::vector<std::string> fl = S.folders.foldersOf(k);
    std::string renFrom, renTo, delF;
    for (const auto& f : fl) {
        int cnt = 0; for (int i = 0; i < n; ++i) if (S.folders.folderOf(k, nameOf(i)) == f) ++cnt;
        char b[200]; std::snprintf(b, sizeof b, "%s %s (%d)###fd_%s", ICON_FOLDER, f.c_str(), cnt, f.c_str());
        if (S.treeExpand) ImGui::SetNextItemOpen(true);
        const bool open = ImGui::TreeNodeEx(b, ImGuiTreeNodeFlags_SpanAvailWidth);
        dropTarget(S, k, f);
        ImGui::PushID(f.c_str());
        if (ImGui::BeginPopupContextItem("##fdc")) {
            char rb[64]; std::snprintf(rb, sizeof rb, "%s", f.c_str());
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
            if (ImGui::InputText(TR("name"), rb, sizeof rb, ImGuiInputTextFlags_EnterReturnsTrue) && rb[0]) { renFrom = f; renTo = rb; }
            if (k == ObjKind::Well) {
                std::vector<std::string> in; for (int i = 0; i < n; ++i) if (S.folders.folderOf(k, nameOf(i)) == f) in.push_back(nameOf(i));
                if (ImGui::MenuItem(TR(ICON_CHART_LINE " Well section: only this folder"))) {
                    S.wsHidden.clear();
                    for (int i = 0; i < n; ++i) if (std::find(in.begin(), in.end(), nameOf(i)) == in.end()) S.wsHidden.push_back(nameOf(i));
                    S.view = 2;
                }
                if (ImGui::MenuItem(TR(ICON_EYE_OFF " Well section: hide this folder")))
                    for (const auto& w : in) if (std::find(S.wsHidden.begin(), S.wsHidden.end(), w) == S.wsHidden.end()) S.wsHidden.push_back(w);
            }
            ImGui::Separator();
            if (ImGui::MenuItem(TR(ICON_TRASH_2 " Remove folder (objects stay)"))) delF = f;
            ImGui::EndPopup();
        }
        ImGui::PopID();
        if (open) {
            for (int i = 0; i < n; ++i) if (S.folders.folderOf(k, nameOf(i)) == f) item(i);
            if (cnt == 0) ImGui::TextDisabled("  %s", TR("drag objects here"));
            ImGui::TreePop();
        }
    }
    for (int i = 0; i < n; ++i) {
        const std::string f = S.folders.folderOf(k, nameOf(i));
        if (f.empty() || std::find(fl.begin(), fl.end(), f) == fl.end()) item(i);
    }
    if (!renFrom.empty()) S.folders.renameFolder(k, renFrom, renTo);
    if (!delF.empty()) S.folders.removeFolder(k, delF);
}
}  // namespace

void geoObjectTree(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    ImGui::TextDisabled("%s", TR("PROJECT OBJECTS"));
    ImGui::BeginChild("##objtree", ImVec2(0, std::max(fs * 6, ImGui::GetContentRegionAvail().y - fs * 4.5f)));
    const ImGuiTreeNodeFlags leaf = ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanAvailWidth;
    char b[160];
    // ---- input data (C46: user folders per kind, drag and drop onto a folder / the kind node, "Move to folder")
    if (ImGui::TreeNodeEx(TR("Input###gtInput"), ImGuiTreeNodeFlags_DefaultOpen)) {
        std::snprintf(b, sizeof b, ICON_LOCATE_FIXED " %s (%zu)###gtWells", TR("Wells"), S.wells.size());
        {
            if (S.treeExpand) ImGui::SetNextItemOpen(true);
            const bool open = ImGui::TreeNode(b);
            kindNodeExtras(S, ObjKind::Well);
            if (open) {
                int del = -1;
                kindTree(S, ObjKind::Well, int(S.wells.size()), [&](int w) { return S.wells[size_t(w)].name; }, [&](int w) {
                    GeoWell& W = S.wells[size_t(w)];
                    ImGui::PushID(w);
                    size_t nc = 0; for (const auto& l : W.logs) nc += l.curves.size();
                    std::snprintf(b, sizeof b, "%s  (%zu logs, %zu tops%s)", W.name.c_str(), nc, W.tops.size(), W.traj.empty() ? "" : ", dev");
                    ImGui::TreeNodeEx(b, leaf);
                    dragSource(ObjKind::Well, w, W.name);
                    if (ImGui::IsItemClicked()) { S.view = 2; S.wsHidden.erase(std::remove(S.wsHidden.begin(), S.wsHidden.end(), W.name), S.wsHidden.end()); }
                    if (ImGui::BeginPopupContextItem("##wc")) {
                        const std::string old = W.name;
                        if (renameField(W.name)) S.folders.renameObject(ObjKind::Well, old, W.name);
                        if (ImGui::MenuItem(TR(ICON_PENCIL " Edit (well manager)"))) { S.we.open = true; S.we.sel = w; }
                        if (ImGui::MenuItem(TR(ICON_CHART_LINE " Show in the well section"))) { S.view = 2; S.wsHidden.erase(std::remove(S.wsHidden.begin(), S.wsHidden.end(), W.name), S.wsHidden.end()); }
                        if (ImGui::MenuItem(TR(ICON_EYE_OFF " Hide in the well section"))) S.wsHidden.push_back(W.name);
                        if (!W.traj.empty() && ImGui::MenuItem(TR(ICON_SCAN_LINE " Intersection along the trajectory"))) { S.isecPath = int(S.polys.size()) + w; S.view = 6; }
                        moveToMenu(S, ObjKind::Well, W.name);
                        ImGui::Separator();
                        if (ImGui::MenuItem(TR(ICON_TRASH_2 " Delete"))) del = w;
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                });
                if (del >= 0) { S.msg("well " + S.wells[size_t(del)].name + " deleted"); S.folders.forget(ObjKind::Well, S.wells[size_t(del)].name); S.wells.erase(S.wells.begin() + del); }
                ImGui::TreePop();
            }
        }
        std::snprintf(b, sizeof b, ICON_MOUNTAIN " %s (%zu)###gtSurf", TR("Surfaces"), S.hz.size());
        {
            if (S.treeExpand) ImGui::SetNextItemOpen(true);
            const bool open = ImGui::TreeNode(b);
            kindNodeExtras(S, ObjKind::Surface);
            if (open) {
                int del = -1;
                kindTree(S, ObjKind::Surface, int(S.hz.size()), [&](int h) { return S.hz[size_t(h)].name; }, [&](int h) {
                    Surface& H = S.hz[size_t(h)];
                    ImGui::PushID(h);
                    const auto* tp = geoTemplate(S, H.name);
                    std::snprintf(b, sizeof b, "%s  %dx%d", H.name.c_str(), H.nx, H.ny);
                    ImGui::TreeNodeEx(b, leaf | (S.view == 0 && S.mapSurface == h ? ImGuiTreeNodeFlags_Selected : 0));
                    dragSource(ObjKind::Surface, h, H.name);
                    if (ImGui::IsItemClicked()) { S.mapSurface = h; S.view = 0; }
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("template: %s", tp ? tp->name.c_str() : "-");
                    if (ImGui::BeginPopupContextItem("##hc")) {
                        const std::string old = H.name;
                        if (renameField(H.name)) { ++S.hzRev; S.folders.renameObject(ObjKind::Surface, old, H.name); }
                        if (ImGui::MenuItem(TR(ICON_IMAGE " Show on the map"))) { S.mapSurface = h; S.view = 0; }
                        if (ImGui::MenuItem(TR(ICON_SIGMA " Statistics"))) {
                            double lo, hi, mean; int n; rangeOf(H.z, lo, hi, mean, n);
                            std::snprintf(b, sizeof b, "%s: %d defined nodes of %d, min %.4g max %.4g mean %.4g", H.name.c_str(), n, H.nx * H.ny, lo, hi, mean);
                            S.msg(b);
                        }
                        if (ImGui::MenuItem(TR(ICON_SAVE " Export IRAP..."))) {
                            std::string p = H.name + ".irap";
                            if (plat::fileDialog(true, "IRAP classic (*.irap)\0*.irap;*.grd\0", "irap", p)) S.msg(writeIrap(p, H) ? "written " + p : "Error: cannot write " + p);
                        }
                        if (ImGui::MenuItem(TR(ICON_PALETTE " Template settings..."))) S.showTemplates = true;
                        moveToMenu(S, ObjKind::Surface, H.name);
                        ImGui::Separator();
                        if (ImGui::MenuItem(TR(ICON_TRASH_2 " Delete"))) del = h;
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                });
                if (del >= 0) { S.pushUndo(); S.folders.forget(ObjKind::Surface, S.hz[size_t(del)].name); S.hz.erase(S.hz.begin() + del); ++S.hzRev; S.mapSurface = 0; }
                ImGui::TreePop();
            }
        }
        std::snprintf(b, sizeof b, ICON_ROUTE " %s (%zu)###gtPoly", TR("Polygons"), S.polys.size());
        {
            if (S.treeExpand) ImGui::SetNextItemOpen(true);
            const bool open = ImGui::TreeNode(b);
            kindNodeExtras(S, ObjKind::Polygon);
            if (open) {
                int del = -1;
                kindTree(S, ObjKind::Polygon, int(S.polys.size()), [&](int p) { return S.polys[size_t(p)].name; }, [&](int p) {
                    GPoly& P = S.polys[size_t(p)];
                    ImGui::PushID(p);
                    std::snprintf(b, sizeof b, "%s  (%s)", P.name.c_str(), P.closed ? "polygon" : "path");
                    ImGui::TreeNodeEx(b, leaf | (S.curPoly == p ? ImGuiTreeNodeFlags_Selected : 0));
                    dragSource(ObjKind::Polygon, p, P.name);
                    if (ImGui::IsItemClicked()) { S.curPoly = p; S.view = 0; }
                    if (ImGui::BeginPopupContextItem("##pc")) {
                        const std::string old = P.name;
                        if (renameField(P.name)) S.folders.renameObject(ObjKind::Polygon, old, P.name);
                        if (ImGui::MenuItem(TR(ICON_SCAN_LINE " Intersection along it"))) { S.isecPath = p; S.view = 6; }
                        if (P.closed && ImGui::MenuItem(TR(ICON_SIGMA " Limit the volumetrics"), nullptr, S.volPoly == p)) S.volPoly = S.volPoly == p ? -1 : p;
                        if (ImGui::MenuItem(TR(ICON_SAVE " Export..."))) {
                            std::string f = P.name + ".txt";
                            if (plat::fileDialog(true, "Polygons (*.txt)\0*.txt\0", "txt", f)) S.msg(writePolygons(f, {P}) ? "written " + f : "Error: cannot write " + f);
                        }
                        moveToMenu(S, ObjKind::Polygon, P.name);
                        ImGui::Separator();
                        if (ImGui::MenuItem(TR(ICON_TRASH_2 " Delete"))) del = p;
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                });
                if (del >= 0) {
                    S.folders.forget(ObjKind::Polygon, S.polys[size_t(del)].name);
                    S.polys.erase(S.polys.begin() + del);
                    if (S.volPoly == del) S.volPoly = -1; else if (S.volPoly > del) --S.volPoly;
                    if (S.isecPath == del) S.isecPath = -1; else if (S.isecPath > del) --S.isecPath;
                    S.curPoly = -1;
                }
                ImGui::TreePop();
            }
        }
        std::snprintf(b, sizeof b, ICON_SPLIT " %s (%zu)###gtFaults", TR("Faults"), S.faults.size());
        {
            if (S.treeExpand) ImGui::SetNextItemOpen(true);
            const bool open = ImGui::TreeNode(b);
            kindNodeExtras(S, ObjKind::Fault);
            if (open) {
                int del = -1;
                kindTree(S, ObjKind::Fault, int(S.faults.size()), [&](int f) { return S.faults[size_t(f)].name; }, [&](int f) {
                    ImGui::PushID(f);
                    ImGui::TreeNodeEx(S.faults[size_t(f)].name.c_str(), leaf);
                    dragSource(ObjKind::Fault, f, S.faults[size_t(f)].name);
                    if (ImGui::BeginPopupContextItem("##fc")) {
                        const std::string old = S.faults[size_t(f)].name;
                        if (renameField(S.faults[size_t(f)].name)) S.folders.renameObject(ObjKind::Fault, old, S.faults[size_t(f)].name);
                        moveToMenu(S, ObjKind::Fault, S.faults[size_t(f)].name);
                        if (ImGui::MenuItem(TR(ICON_TRASH_2 " Delete"))) del = f;
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                });
                if (del >= 0) { S.folders.forget(ObjKind::Fault, S.faults[size_t(del)].name); S.faults.erase(S.faults.begin() + del); }
                ImGui::TreePop();
            }
        }
        std::snprintf(b, sizeof b, ICON_ACTIVITY " %s (%s%zu %s)###gtSeis", TR("Seismic"), S.hasCube ? "3-D, " : "", S.lines.size(), TR("lines"));
        {
            if (S.treeExpand) ImGui::SetNextItemOpen(true);
            const bool open = ImGui::TreeNode(b);
            kindNodeExtras(S, ObjKind::Line);
            if (open) {
                if (S.hasCube) {
                    std::snprintf(b, sizeof b, "%s  %dx%dx%d", S.cube.name.c_str(), S.cube.ni, S.cube.nj, S.cube.ns);
                    ImGui::TreeNodeEx(b, leaf | (S.curLine < 0 && S.view == 1 ? ImGuiTreeNodeFlags_Selected : 0));
                    if (ImGui::IsItemClicked()) { S.curLine = -1; S.view = 1; }
                    if (ImGui::BeginPopupContextItem("##cc")) {
                        if (ImGui::MenuItem(TR(ICON_ACTIVITY " Seismic section"))) { S.curLine = -1; S.view = 1; }
                        if (ImGui::MenuItem(TR(ICON_BOX " 3-D scene"))) S.view = 4;
                        if (ImGui::MenuItem(TR(ICON_SIGMA " Calculator / crop"))) { S.step = 1; S.openTool = 2; }
                        if (ImGui::MenuItem(TR(ICON_SAVE " Export SEG-Y..."))) {
                            std::string p = "volume.sgy";
                            if (plat::fileDialog(true, "SEG-Y (*.sgy)\0*.sgy;*.segy\0", "sgy", p)) S.msg(writeSegy(p, S.cube) ? "written " + p : "Error: cannot write " + p);
                        }
                        ImGui::EndPopup();
                    }
                }
                kindTree(S, ObjKind::Line, int(S.lines.size()), [&](int l) { return S.lines[size_t(l)].name; }, [&](int l) {
                    ImGui::PushID(l);
                    ImGui::TreeNodeEx(S.lines[size_t(l)].name.c_str(), leaf | (S.curLine == l ? ImGuiTreeNodeFlags_Selected : 0));
                    dragSource(ObjKind::Line, l, S.lines[size_t(l)].name);
                    if (ImGui::IsItemClicked()) { S.curLine = l; S.view = 1; }
                    if (ImGui::BeginPopupContextItem("##lc")) {
                        const std::string old = S.lines[size_t(l)].name;
                        if (renameField(S.lines[size_t(l)].name)) S.folders.renameObject(ObjKind::Line, old, S.lines[size_t(l)].name);
                        if (ImGui::MenuItem(TR(ICON_ACTIVITY " Line section"))) { S.curLine = l; S.view = 1; }
                        if (S.lines.size() > 1 && ImGui::MenuItem(TR(ICON_GIT_COMPARE " Mis-tie correction..."))) { S.step = 1; S.openTool = 1; S.mistieRef = l; }
                        moveToMenu(S, ObjKind::Line, S.lines[size_t(l)].name);
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                });
                ImGui::TreePop();
            }
        }
        ImGui::TreePop();
    }
    // ---- models
    if (ImGui::TreeNodeEx(TR("Models###gtModels"), ImGuiTreeNodeFlags_DefaultOpen)) {
        if (!S.hasGrid) ImGui::TextDisabled("  %s", TR("no grid (Grid step)"));
        else {
            std::snprintf(b, sizeof b, ICON_GRID_3X3 " Grid %dx%dx%d", S.g.grid.nx, S.g.grid.ny, S.g.grid.nz);
            if (ImGui::TreeNode(b)) {
                for (int p = 0; p < geoGridPropCount(S); ++p) {
                    ImGui::PushID(p);
                    const std::string nm = geoGridPropName(S, p);
                    ImGui::TreeNodeEx(nm.c_str(), leaf | (S.gridProp == p ? ImGuiTreeNodeFlags_Selected : 0));
                    if (ImGui::IsItemClicked()) { S.gridProp = p; if (S.view != 3 && S.view != 4 && S.view != 6) S.view = 3; }
                    if (ImGui::BeginPopupContextItem("##gpc")) {
                        if (ImGui::MenuItem(TR(ICON_LAYERS " Grid cross-section"))) { S.gridProp = p; S.view = 3; }
                        if (ImGui::MenuItem(TR(ICON_BOX " 3-D scene"))) { S.gridProp = p; S.view = 4; }
                        if (ImGui::MenuItem(TR(ICON_SCAN_LINE " Intersection"))) { S.gridProp = p; S.view = 6; }
                        if (ImGui::MenuItem(TR(ICON_SIGMA " Statistics"))) {
                            double lo, hi, mean; int n; rangeOf(geoGridProperty(S, p), lo, hi, mean, n);
                            std::snprintf(b, sizeof b, "%s: %d cells, min %.4g max %.4g mean %.4g", nm.c_str(), n, lo, hi, mean);
                            S.msg(b);
                        }
                        if (ImGui::MenuItem(TR(ICON_PALETTE " Template settings..."))) S.showTemplates = true;
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                }
                ImGui::TreePop();
            }
        }
        ImGui::TreePop();
    }
    // ---- results
    if (ImGui::TreeNodeEx(TR("Results###gtResults"), ImGuiTreeNodeFlags_DefaultOpen)) {
        if (S.hasVol) { ImGui::TreeNodeEx(ICON_SIGMA " Volumetrics", leaf); if (ImGui::IsItemClicked()) S.step = 7; }
        if (S.lcHas) { ImGui::TreeNodeEx(ICON_LAYERS " Layer-cake velocity model", leaf); if (ImGui::IsItemClicked()) { S.step = 1; S.openTool = 3; } }
        if (S.mistieHas) { ImGui::TreeNodeEx(ICON_GIT_COMPARE " Mis-tie fit", leaf); if (ImGui::IsItemClicked()) { S.step = 1; S.openTool = 1; } }
        if (!S.hasVol && !S.lcHas && !S.mistieHas) ImGui::TextDisabled("  %s", TR("none yet"));
        ImGui::TreePop();
    }
    ImGui::EndChild();
    S.treeExpand = false;
}
