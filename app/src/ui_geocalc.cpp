// C42 Geomodel > Analysis: Petrel-style calculators and object statistics
//   Surface calculator  - expression over every horizon (resampled to one lattice) -> new or replaced surface
//   Property calculator - expression over the grid properties -> new user property, or overwrite porosity /
//                         permeability / NTG / facies (the edit Petrel users do with the property calculator)
//   Statistics          - n / min / max / mean / sd / P10 / P50 / P90 of every horizon and grid property + histogram
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "widgets.h"
#include "nav.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "plotkit.h"
#include "ressim/geo_analysis.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {
using namespace ressim;
using namespace ressim::geo;

struct Calc {
    char sExpr[512] = "", sName[64] = "Calc surface";
    int sOnto = 0; bool sReplace = false;
    std::string sKey, sErr; Stats sPrev;
    char pExpr[512] = "", pName[64] = "PROP1";
    int pTarget = 0;   // 0 new property, 1 porosity, 2 permeability X, 3 NTG, 4 facies, 5.. user property (overwrite)
    std::string pKey, pErr; Stats pPrev;
    int stSel = 0; int bins = 30;
};
Calc g_calc;

ImVec4 col4(ImU32 c) { return ImGui::ColorConvertU32ToFloat4(c); }
void append(char* buf, size_t cap, const std::string& s) {
    const size_t L = std::strlen(buf);
    std::snprintf(buf + L, cap - L, "%s%s", L && buf[L - 1] != ' ' && buf[L - 1] != '(' ? " " : "", s.c_str());
}
void statsLine(const Stats& s) {
    if (!s.n) { ImGui::TextDisabled("no defined values"); return; }
    ImGui::Text("n %zu  (undefined %zu)   min %.5g   max %.5g   mean %.5g   sd %.4g   P10 %.5g   P50 %.5g   P90 %.5g",
                s.n, s.nNull, s.min, s.max, s.mean, s.sd, s.p10, s.p50, s.p90);
}

// grid property variables for the calculator (inactive cells NaN)
struct GridVars {
    std::vector<std::pair<std::string, std::vector<double>>> v;
    std::map<std::string, const std::vector<double>*> map;
};
void gridVars(const App::GeoState& S, GridVars& G) {
    const GeoGrid& g = S.g;
    const size_t n = g.bulk.size();
    const int nx = g.grid.nx, ny = g.grid.ny;
    auto add = [&](const std::string& name, std::vector<double> x) { G.v.emplace_back(name, std::move(x)); };
    const double nan = std::nan("");
    std::vector<double> poro(n, nan), kx(n, nan), kz(n, nan), ntg(n, nan), fac(n, nan), zone(n, nan), blk(n, nan), dep(n, nan),
        top(n, nan), th(n, nan), I(n), J(n), K(n), X(n, nan), Y(n, nan), bulk(n, nan);
    for (size_t c = 0; c < n; ++c) {
        I[c] = double(c % size_t(nx) + 1); J[c] = double((c / size_t(nx)) % size_t(ny) + 1); K[c] = double(c / size_t(nx * ny) + 1);
        if (!g.grid.actnum.empty() && !g.grid.actnum[c]) continue;
        if (S.hasPetro && c < S.phi.v.size()) poro[c] = S.phi.v[c];
        if (S.hasPetro && c < S.kx.size()) kx[c] = S.kx[c];
        if (S.hasPetro && c < S.kz.size()) kz[c] = S.kz[c];
        if (S.hasPetro && c < S.ntg.size()) ntg[c] = S.ntg[c];
        if (S.hasFacies && c < S.fr.facies.size()) fac[c] = S.fr.facies[c];
        const size_t k = c / size_t(nx * ny);
        zone[c] = k < g.layerZone.size() ? g.layerZone[k] : 0;
        blk[c] = g.block.empty() ? 0 : g.block[c % size_t(nx * ny)];
        if (c < g.center.size()) { X[c] = g.center[c][0]; Y[c] = g.center[c][1]; dep[c] = g.center[c][2]; }
        top[c] = g.top[c]; th[c] = g.bot[c] - g.top[c]; bulk[c] = g.bulk[c];
    }
    add("PORO", std::move(poro)); add("PERMX", std::move(kx)); add("PERMZ", std::move(kz)); add("NTG", std::move(ntg));
    add("FACIES", std::move(fac)); add("ZONE", std::move(zone)); add("BLOCK", std::move(blk)); add("DEPTH", std::move(dep));
    add("TOP", std::move(top)); add("THICK", std::move(th)); add("BULK", std::move(bulk));
    add("I", std::move(I)); add("J", std::move(J)); add("K", std::move(K)); add("X", std::move(X)); add("Y", std::move(Y));
    for (const auto& u : S.userProps)
        if (u.second.size() == n) add(exprName(u.first), u.second);
    for (const auto& p : G.v) G.map[p.first] = &p.second;
}

void surfaceTab(App::GeoState& S) {
    Calc& C = g_calc;
    const float fs = ImGui::GetFontSize();
    std::vector<int> valid;
    for (int i = 0; i < int(S.hz.size()); ++i) if (S.hz[size_t(i)].valid()) valid.push_back(i);
    if (valid.empty()) {
        ui::emptyState(ICON_LAYERS, "No horizons", "Load or build horizons on the Data / Horizons steps; the surface calculator combines them "
                                                   "(isochores, shifted or clipped surfaces, conditional edits).");
        return;
    }
    if (C.sOnto < 0 || C.sOnto >= int(S.hz.size()) || !S.hz[size_t(C.sOnto)].valid()) C.sOnto = valid[0];
    ImGui::SetNextItemWidth(fs * 12);
    if (ImGui::BeginCombo("Lattice of", S.hz[size_t(C.sOnto)].name.c_str())) {
        for (int i : valid) if (ImGui::Selectable(S.hz[size_t(i)].name.c_str(), i == C.sOnto)) C.sOnto = i;
        ImGui::EndCombo();
    }
    ui::hint("All horizons are resampled (bilinear) onto the nodes of this surface; outside a surface's extent its value is "
             "undefined (nan) and so is every result that uses it.");
    ImGui::SameLine(); ImGui::Checkbox(TR("replace it"), &C.sReplace);
    if (!C.sReplace) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 10); ImGui::InputText("Result name", C.sName, sizeof C.sName); }
    ImGui::SetNextItemWidth(-1); ImGui::InputTextWithHint("##sexpr", "e.g.  BASERES - TOPRES      or      if(X > 452000, TOPRES + 10, TOPRES)", C.sExpr, sizeof C.sExpr);
    ImGui::TextDisabled("Surfaces:");
    for (int i : valid) {
        ImGui::SameLine();
        const std::string n = exprName(S.hz[size_t(i)].name);
        if (ImGui::SmallButton((n + "##sv" + std::to_string(i)).c_str())) append(C.sExpr, sizeof C.sExpr, n);
    }
    for (const char* n : {"X", "Y"}) { ImGui::SameLine(); if (ImGui::SmallButton(n)) append(C.sExpr, sizeof C.sExpr, n); }
    ImGui::TextDisabled("Functions: abs sqrt ln log10 exp pow min max clamp if(c,a,b) isnan  |  + - * / ^ < <= > >= == != && || !  |  nan, pi");
    const std::string key = std::string(C.sExpr) + "|" + std::to_string(C.sOnto) + "|" + std::to_string(S.hzRev) + "|" + std::to_string(S.hz.size());
    Surface res;
    bool ok = false;
    if (key != C.sKey || !C.sErr.empty()) {
        C.sKey = key; C.sErr.clear(); C.sPrev = Stats{};
        if (C.sExpr[0]) {
            ok = surfaceCalculator(S.hz, C.sOnto, C.sExpr, res, C.sErr);
            if (ok) C.sPrev = describe(res.z);
        }
    }
    if (!C.sErr.empty()) ImGui::TextColored(col4(Pal::Err), ICON_TRIANGLE_ALERT " %s", C.sErr.c_str());
    else if (C.sExpr[0]) { ImGui::TextDisabled("Preview:"); ImGui::SameLine(); statsLine(C.sPrev); }
    ImGui::BeginDisabled(!C.sExpr[0] || !C.sErr.empty() || (!C.sReplace && !C.sName[0]));
    if (ImGui::Button(C.sReplace ? ICON_CHECK " Replace surface" : ICON_PLUS " Create surface")) {
        std::string err;
        if (!ok) ok = surfaceCalculator(S.hz, C.sOnto, C.sExpr, res, err);
        if (ok) {
            if (C.sReplace) {
                S.pushUndo();
                res.name = S.hz[size_t(C.sOnto)].name;
                S.hz[size_t(C.sOnto)] = res;
                S.msg("surface " + res.name + " replaced by " + C.sExpr + (S.hasGrid ? " - rebuild the grid to apply" : ""));
            } else {
                res.name = C.sName;
                bool rep = false;
                for (auto& h : S.hz) if (h.name == res.name) { h = res; rep = true; }
                if (!rep) S.hz.push_back(res);
                S.msg("surface " + res.name + " = " + C.sExpr);
            }
            ++S.hzRev;
        } else S.msg("surface calculator: " + err);
    }
    ImGui::EndDisabled();
    if (C.sReplace) { ImGui::SameLine(); ImGui::TextDisabled("(undo with the horizon-edit undo)"); }
}

void propertyTab(App::GeoState& S) {
    Calc& C = g_calc;
    const float fs = ImGui::GetFontSize();
    if (!S.hasGrid) {
        ui::emptyState(ICON_GRID_3X3, "No grid", "Build the grid (Grid step) first; the property calculator works on its cells.");
        return;
    }
    std::vector<std::string> targets = {"New property", "Porosity (overwrite)", "Permeability X, mD (overwrite)", "NTG (overwrite)", "Facies (overwrite)"};
    for (const auto& u : S.userProps) targets.push_back(u.first + " (overwrite)");
    C.pTarget = std::clamp(C.pTarget, 0, int(targets.size()) - 1);
    ImGui::SetNextItemWidth(fs * 14);
    if (ImGui::BeginCombo("Result", targets[size_t(C.pTarget)].c_str())) {
        for (int i = 0; i < int(targets.size()); ++i) if (ImGui::Selectable(targets[size_t(i)].c_str(), i == C.pTarget)) C.pTarget = i;
        ImGui::EndCombo();
    }
    if (C.pTarget == 0) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9); ImGui::InputText("Name", C.pName, sizeof C.pName); }
    const bool needPetro = C.pTarget >= 1 && C.pTarget <= 3, needFac = C.pTarget == 4;
    ui::hint("Undefined results (nan) leave the cell unchanged when overwriting. Overwriting permeability X scales permeability Z "
             "by the same factor (kv/kh kept). Facies results are rounded to integer codes. Volumetrics and upscaling must be "
             "re-run after an edit.");
    ImGui::SetNextItemWidth(-1);
    ImGui::InputTextWithHint("##pexpr", "e.g.  if(FACIES == 1 && ZONE == 2, PORO * 0.9, PORO)      or      PORO > 0.2 && PERMX > 50", C.pExpr, sizeof C.pExpr);
    GridVars G;
    gridVars(S, G);
    ImGui::TextDisabled("Properties:");
    for (const auto& p : G.v) {
        ImGui::SameLine();
        if (ImGui::GetContentRegionAvail().x < fs * 4) ImGui::NewLine();
        if (ImGui::SmallButton(p.first.c_str())) append(C.pExpr, sizeof C.pExpr, p.first);
    }
    ImGui::TextDisabled("I J K are 1-based; DEPTH = cell centre depth, TOP = cell top, THICK = cell thickness, BULK = bulk volume");
    const std::string key = std::string(C.pExpr) + "|" + std::to_string(S.propRev) + "|" + std::to_string(size_t(S.phi.v.data())) + "|" +
                            std::to_string(size_t(S.fr.facies.data())) + "|" + std::to_string(S.userProps.size()) + "|" + std::to_string(S.g.bulk.size());
    std::vector<double> res;
    bool ok = false;
    if (key != C.pKey) {
        C.pKey = key; C.pErr.clear(); C.pPrev = Stats{};
        if (C.pExpr[0]) {
            ok = evalExpression(C.pExpr, G.map, S.g.bulk.size(), res, C.pErr);
            if (ok) C.pPrev = describe(res);
        }
    }
    if (!C.pErr.empty()) ImGui::TextColored(col4(Pal::Err), ICON_TRIANGLE_ALERT " %s", C.pErr.c_str());
    else if (C.pExpr[0]) { ImGui::TextDisabled("Preview:"); ImGui::SameLine(); statsLine(C.pPrev); }
    const bool blocked = (needPetro && !S.hasPetro) || (needFac && !S.hasFacies);
    if (blocked) ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT " run the %s step first", needFac ? "Facies" : "Petrophysics");
    ImGui::BeginDisabled(!C.pExpr[0] || !C.pErr.empty() || blocked || (C.pTarget == 0 && !C.pName[0]));
    if (ImGui::Button(TR(ICON_CHECK " Apply"))) {
        std::string err;
        if (!ok) ok = evalExpression(C.pExpr, G.map, S.g.bulk.size(), res, err);
        if (!ok) S.msg("property calculator: " + err);
        else {
            size_t changed = 0;
            auto over = [&](std::vector<double>& dst) {
                for (size_t c = 0; c < res.size() && c < dst.size(); ++c) if (std::isfinite(res[c])) { dst[c] = res[c]; ++changed; }
            };
            std::string what;
            if (C.pTarget == 0) {
                what = C.pName;
                bool rep = false;
                for (auto& u : S.userProps) if (u.first == what) { u.second = res; rep = true; }
                if (!rep) S.userProps.emplace_back(what, res);
                for (double x : res) changed += std::isfinite(x) ? 1 : 0;
                for (size_t u = 0; u < S.userProps.size(); ++u) if (S.userProps[u].first == what) S.gridProp = 7 + int(u);
            } else if (C.pTarget == 1) { what = "porosity"; over(S.phi.v); S.gridProp = 2; }
            else if (C.pTarget == 2) {
                what = "permeability X";
                for (size_t c = 0; c < res.size() && c < S.kx.size(); ++c) {
                    if (!std::isfinite(res[c]) || res[c] < 0) continue;
                    if (c < S.kz.size()) S.kz[c] = S.kx[c] > 0 ? S.kz[c] * res[c] / S.kx[c] : S.kvkh * res[c];
                    S.kx[c] = res[c]; ++changed;
                }
                S.gridProp = 3;
            } else if (C.pTarget == 3) { what = "NTG"; over(S.ntg); }
            else if (C.pTarget == 4) {
                what = "facies";
                for (size_t c = 0; c < res.size() && c < S.fr.facies.size(); ++c)
                    if (std::isfinite(res[c])) { S.fr.facies[c] = int(std::lround(res[c])); ++changed; }
                S.gridProp = 1;
            } else {
                const size_t u = size_t(C.pTarget - 5);
                if (u < S.userProps.size()) { what = S.userProps[u].first; over(S.userProps[u].second); S.gridProp = 7 + int(u); }
            }
            if (C.pTarget != 0) S.hasVol = S.hasUp = false;
            ++S.propRev;
            C.pKey.clear();
            S.msg(what + " = " + C.pExpr + ": " + std::to_string(changed) + " cells");
        }
    }
    ImGui::EndDisabled();
    if (!S.userProps.empty()) {
        ImGui::SeparatorText(TR("Calculated properties"));
        for (size_t u = 0; u < S.userProps.size(); ++u) {
            ImGui::PushID(int(u));
            if (ImGui::SmallButton(TR("X"))) { S.userProps.erase(S.userProps.begin() + long(u)); S.gridProp = 2; ++S.propRev; ImGui::PopID(); break; }
            ImGui::SameLine(); ImGui::TextUnformatted(S.userProps[u].first.c_str());
            ImGui::SameLine(fs * 10); statsLine(describe(S.userProps[u].second));
            ImGui::PopID();
        }
        ImGui::TextDisabled("Calculated properties appear in the property lists of the map, sections and the 3-D scene.");
    }
}

void statisticsTab(App::GeoState& S) {
    Calc& C = g_calc;
    const float fs = ImGui::GetFontSize();
    struct Row { std::string kind, name; Stats s; std::vector<double> v; };
    std::vector<Row> rows;
    for (const auto& h : S.hz) if (h.valid()) rows.push_back({"Surface", h.name, describe(h.z), h.z});
    if (S.hasGrid)
        for (int p = 0; p < geoGridPropCount(S); ++p) {
            auto v = geoGridProperty(S, p);
            const Stats st = describe(v);
            if (st.n) rows.push_back({"Grid", geoGridPropName(S, p), st, std::move(v)});
        }
    if (rows.empty()) { ui::emptyState(ICON_TABLE, "Nothing to describe", "Load horizons or build a grid first."); return; }
    C.stSel = std::clamp(C.stSel, 0, int(rows.size()) - 1);
    const float th = std::min(fs * 1.6f * float(rows.size() + 1), ImGui::GetContentRegionAvail().y * 0.45f);
    if (ImGui::BeginTable("##gst", 10, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable, ImVec2(-1, th))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        for (const char* h : {"Object", "Kind", "n", "min", "max", "mean", "sd", "P10", "P50", "P90"}) ImGui::TableSetupColumn(h);
        ImGui::TableHeadersRow();
        for (int r = 0; r < int(rows.size()); ++r) {
            const Row& R = rows[size_t(r)];
            ImGui::TableNextRow(); ImGui::TableNextColumn();
            if (ImGui::Selectable((R.name + "##st" + std::to_string(r)).c_str(), r == C.stSel, ImGuiSelectableFlags_SpanAllColumns)) C.stSel = r;
            ImGui::TableNextColumn(); ImGui::TextUnformatted(R.kind.c_str());
            ImGui::TableNextColumn(); ImGui::Text("%zu", R.s.n);
            for (double x : {R.s.min, R.s.max, R.s.mean, R.s.sd, R.s.p10, R.s.p50, R.s.p90}) { ImGui::TableNextColumn(); ImGui::Text("%.5g", x); }
        }
        ImGui::EndTable();
    }
    if (ImGui::SmallButton(TR(ICON_COPY " Copy table"))) {
        std::string t = "Object\tKind\tn\tmin\tmax\tmean\tsd\tP10\tP50\tP90\n";
        char b[256];
        for (const auto& R : rows) {
            std::snprintf(b, sizeof b, "\t%zu\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n", R.s.n, R.s.min, R.s.max, R.s.mean, R.s.sd, R.s.p10, R.s.p50, R.s.p90);
            t += R.name + "\t" + R.kind + b;
        }
        ImGui::SetClipboardText(t.c_str());
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8); ImGui::SliderInt("bins", &C.bins, 5, 100);
    const Row& R = rows[size_t(C.stSel)];
    std::vector<double> f;
    for (double x : R.v) if (std::isfinite(x)) f.push_back(x);
    if (f.empty()) return;
    const double lo = R.s.min, hi = R.s.max > R.s.min ? R.s.max : R.s.min + 1;
    if (ImPlot::BeginPlot("##gsth", ImVec2(-1, std::max(fs * 10, ImGui::GetContentRegionAvail().y)))) {
        recPanel(("Histogram " + R.name).c_str(), R.name.c_str(), "count", "");
        ImPlot::SetupAxes(R.name.c_str(), "count", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        const ImU32 c = seriesColor(0);
        ImPlotSpec s; s.FillColor = col4(c); s.LineColor = col4(c); s.FillAlpha = 0.6f;
        ImPlot::PlotHistogram(R.name.c_str(), f.data(), int(f.size()), C.bins, 1.0, ImPlotRange(lo, hi), s);
        std::vector<double> cnt(size_t(C.bins), 0.0), bx, by;
        const double bw = (hi - lo) / C.bins;
        for (double x : f) cnt[size_t(std::clamp(int((x - lo) / bw), 0, C.bins - 1))] += 1;
        for (int b = 0; b < C.bins; ++b) { bx.push_back(lo + b * bw); by.push_back(cnt[size_t(b)]); bx.push_back(lo + (b + 1) * bw); by.push_back(cnt[size_t(b)]); }
        recSeries(R.name.c_str(), bx.data(), by.data(), int(bx.size()), c, 1.5f);
        ImPlot::EndPlot();
    }
}
}  // namespace

void geoCalcTabs(App::GeoState& S, int forceTab) {
    auto fl = [&](int i) { return forceTab == i ? ImGuiTabItemFlags_SetSelected : 0; };
    if (forceTab == 0 && !g_calc.sExpr[0]) std::snprintf(g_calc.sExpr, sizeof g_calc.sExpr, "BASERES - TOPRES");
    if (ui::tabItem(ICON_LAYERS " Surface calculator", nullptr, fl(0))) { surfaceTab(S); ImGui::EndTabItem(); }
    if (ui::tabItem(ICON_GRID_3X3 " Property calculator", nullptr, fl(1))) { propertyTab(S); ImGui::EndTabItem(); }
    if (ui::tabItem(ICON_TABLE " Statistics", nullptr, fl(2))) { statisticsTab(S); ImGui::EndTabItem(); }
}
