// C40 Geomodel > Data analysis view (Petrel "Data analysis" / "Well tie" level):
//   Log calculator - expressions over the curves of every well (presets: Vsh, density / sonic porosity, Archie, Timur...)
//   Cross-plot     - curve vs curve over all wells, coloured by well or by a third curve, regression (linear / log-linear)
//                    and "use as the poro-perm transform" of the Petrophysics step
//   Histogram      - per well or all, statistics (P10 / P50 / P90)
//   Variogram      - experimental variogram from the logs (vertical / horizontal / azimuth, stratigraphic depth),
//                    draggable nugget / sill / range handles, weighted least-squares auto fit, apply to Facies / Petrophysics
//   Well tie       - synthetic seismogram (sonic + density -> AI -> reflectivity * Ricker), bulk-shift search against the
//                    seismic trace at the well, time-depth plot and a calibrated velocity model
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "colormaps.h"
#include "plotkit.h"
#include "ressim/geomodel.h"
#include "ressim/geo_analysis.h"
#include "ressim/geostat.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

struct GA {
    int tab = 0;
    int rev = 0;                                  // bumped when curves are added / removed (cache keys)
    // log calculator
    int calcWell = -1;                            // -1 all wells
    char expr[512] = "clamp((GR - 30) / (120 - 30), 0, 1)";
    char mnem[32] = "VSH", unit[16] = "V/V";
    int preset = 0;
    std::string calcMsg, prevKey, prevErr;
    std::vector<double> prevX, prevZ;
    // zone filter (shared by cross-plot, histogram, variogram)
    int zoneTop = -1, zoneBase = -1;
    // cross-plot
    int cx = -1, cy = -1, cz = -1, colorMode = 1, reg = 2;
    bool logX = false, logY = true;
    // histogram
    int hc = -1, bins = 30;
    bool hLog = false, hByWell = false;
    // variogram
    int vc = -1, vTrans = 0, vCode = 1, vDir = 1, vLags = 20, vModel = 0;
    double vAz = 0, vLag = 0;
    bool vStd = true, vStrat = true, fixNug = false;
    geostat::ExpVario ev;
    double evVar = 0, evLagW = 0;
    int evN = 0;
    std::string evKey, evMsg;
    double nug = 0, sill = 1, range = 10, rmse = 0;
    bool haveModel = false, smokeFit = false;
    int applyTo = 1;
    // well tie
    int tw = 0, tdt = -1, trho = -1;
    geo::TieInput ti;
    geo::TieResult tr;
    std::string tKey;
};
GA g_ga;

bool isNull(double v, double nullValue) { return !std::isfinite(v) || std::fabs(v - nullValue) < 1e-6 * std::max(1.0, std::fabs(nullValue)); }

std::vector<std::string> curveNames(const App::GeoState& S) {
    std::vector<std::string> v;
    for (const auto& w : S.wells) for (const auto& l : w.logs) for (const auto& c : l.curves)
        if (std::find(v.begin(), v.end(), c.mnem) == v.end()) v.push_back(c.mnem);
    return v;
}
std::vector<std::string> topNames(const App::GeoState& S) {
    std::vector<std::string> v;
    for (const auto& w : S.wells) for (const auto& t : w.tops) if (std::find(v.begin(), v.end(), t.surface) == v.end()) v.push_back(t.surface);
    return v;
}
int findName(const std::vector<std::string>& v, const char* a, int fallback) {
    for (size_t i = 0; i < v.size(); ++i) if (v[i] == a) return int(i);
    return fallback;
}
bool combo(const char* label, int* idx, const std::vector<std::string>& items, const char* noneLabel = nullptr) {
    const char* cur = *idx >= 0 && *idx < int(items.size()) ? items[size_t(*idx)].c_str() : (noneLabel ? noneLabel : "-");
    bool ch = false;
    if (ImGui::BeginCombo(label, cur)) {
        if (noneLabel && ImGui::Selectable(noneLabel, *idx < 0)) { *idx = -1; ch = true; }
        for (int i = 0; i < int(items.size()); ++i) if (ImGui::Selectable(items[size_t(i)].c_str(), *idx == i)) { *idx = i; ch = true; }
        ImGui::EndCombo();
    }
    return ch;
}
ImPlotSpec lineSpec(ImVec4 c, float w) { ImPlotSpec s; s.LineColor = c; s.LineWeight = w; return s; }
ImPlotSpec dotSpec(ImVec4 c, float r) { ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = r; s.MarkerFillColor = c; s.MarkerLineColor = c; s.LineWeight = 0; return s; }
ImVec4 c4f(ImU32 c) { return ImGui::ColorConvertU32ToFloat4(c); }
bool sliderD(const char* l, double* v, double lo, double hi, const char* fmt) { return ImGui::SliderScalar(l, ImGuiDataType_Double, v, &lo, &hi, fmt); }

// MD interval of the zone filter in one well; false when the well lacks a requested top
bool zoneMd(const geo::GeoWell& w, const std::vector<std::string>& tops, double& md0, double& md1) {
    md0 = -1e300; md1 = 1e300;
    auto find = [&](int i, double& md) {
        if (i < 0 || i >= int(tops.size())) return true;
        for (const auto& t : w.tops) if (t.surface == tops[size_t(i)]) { md = t.md; return true; }
        return false;
    };
    return find(g_ga.zoneTop, md0) && find(g_ga.zoneBase, md1);
}

struct Sample { int well; double x, y, z, zs; double v[3]; };
// every log sample of every well inside the zone with the values of up to three curves (NaN where missing / null);
// zs: depth below the zone top of that well (stratigraphic depth; z when no top filter)
std::vector<Sample> gather(const App::GeoState& S, const std::vector<std::string>& names, const int* idx, int nc, const std::vector<std::string>& tops) {
    std::vector<Sample> out;
    for (int wi = 0; wi < int(S.wells.size()); ++wi) {
        const auto& w = S.wells[size_t(wi)];
        double md0, md1;
        if (!zoneMd(w, tops, md0, md1)) continue;
        const double zTop = md0 > -1e299 ? w.at(md0)[2] : 0.0;
        for (const auto& l : w.logs) {
            const geo::LogCurve* c[3] = {nullptr, nullptr, nullptr};
            bool any = false;
            for (int k = 0; k < nc; ++k) if (idx[k] >= 0 && idx[k] < int(names.size())) { c[k] = l.curve(names[size_t(idx[k])]); any = any || c[k]; }
            if (!any) continue;
            for (size_t s = 0; s < l.md.size(); ++s) {
                if (l.md[s] < md0 || l.md[s] >= md1) continue;
                const geo::P3 p = w.at(l.md[s]);
                Sample q{wi, p[0], p[1], p[2], p[2] - zTop, {NAN, NAN, NAN}};
                for (int k = 0; k < nc; ++k) if (c[k] && s < c[k]->v.size() && !isNull(c[k]->v[s], l.nullValue)) q.v[k] = c[k]->v[s];
                out.push_back(q);
            }
        }
    }
    return out;
}

void zoneFilterUi(const std::vector<std::string>& tops, float fs) {
    ImGui::TextUnformatted("Zone"); ImGui::SameLine();
    ImGui::SetNextItemWidth(fs * 9); combo("##ztop", &g_ga.zoneTop, tops, "(log top)");
    ImGui::SameLine(); ImGui::TextUnformatted("to"); ImGui::SameLine();
    ImGui::SetNextItemWidth(fs * 9); combo("##zbase", &g_ga.zoneBase, tops, "(log base)");
}

void gradientBar(Colormap m, float w, float h) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const int n = 32;
    for (int i = 0; i < n; ++i) {
        float a[3], b[3];
        colormapRGB(m, float(i) / n, a); colormapRGB(m, float(i + 1) / n, b);
        const ImU32 ca = IM_COL32(int(a[0] * 255), int(a[1] * 255), int(a[2] * 255), 255), cb = IM_COL32(int(b[0] * 255), int(b[1] * 255), int(b[2] * 255), 255);
        dl->AddRectFilledMultiColor(ImVec2(p.x + w * i / n, p.y), ImVec2(p.x + w * (i + 1) / n, p.y + h), ca, cb, cb, ca);
    }
    ImGui::Dummy(ImVec2(w, h));
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const double r = q * double(v.size() - 1);
    const size_t i = size_t(r);
    return i + 1 < v.size() ? v[i] + (r - double(i)) * (v[i + 1] - v[i]) : v.back();
}

}  // namespace

void App::geoAnalysis(GeoState& S, ImVec2 size) {
    GA& G = g_ga;
    const float fs = ImGui::GetFontSize();
    (void)size;   // the tabs size their plots from the remaining content region
    if (S.wells.empty()) {   // C42: the calculators / statistics work without wells
        if (ImGui::BeginTabBar("##ga0")) {
            if (ui::tabItem(ICON_CHART_SCATTER " Well logs")) {
                ui::emptyState(ICON_CHART_SCATTER, "No wells loaded",
                               "Log analysis works on the well logs: load LAS files (or the demo project) on the Data step, then compute "
                               "curves, cross-plot them, fit variograms and tie the wells to the seismic here.", {});
                ImGui::EndTabItem();
            }
            geoCalcTabs(S);
            ImGui::EndTabBar();
        }
        return;
    }
    const auto names = curveNames(S);
    const auto tops = topNames(S);
    std::vector<std::string> wellNames;
    for (const auto& w : S.wells) wellNames.push_back(w.name);
    const bool busy = S.job.valid();

    int forceTab = -1;
    if (!smokeAction_.empty() && smokeAction_.rfind("geo-", 0) == 0) {   // C40 smoke tests: open a tab once
        static const char* kSm[] = {"geo-calc", "geo-xplot", "geo-hist", "geo-vario", "geo-tie", "geo-scalc", "geo-pcalc", "geo-stats"};
        for (int i = 0; i < 8; ++i) if (smokeAction_ == kSm[i]) forceTab = i;
        if (forceTab == 3) {   // reservoir zone, vertical PHIE variogram, auto fit once the experimental variogram exists
            G.vc = findName(names, "PHIE", G.vc); G.vDir = 1;
            G.zoneTop = findName(tops, "TopRes", -1); G.zoneBase = findName(tops, "BaseRes", -1);
            G.smokeFit = true;
        }
        smokeAction_.clear();
    }
    auto tabFlags = [&](int i) { return forceTab == i ? ImGuiTabItemFlags_SetSelected : 0; };
    if (!ImGui::BeginTabBar("##ga")) return;

    // ================================================================ log calculator
    if (ui::tabItem(ICON_SIGMA " Log calculator", nullptr, tabFlags(0))) {
        const auto presets = geo::calcPresets();
        ImGui::SetNextItemWidth(fs * 9); combo("Well##cw", &G.calcWell, wellNames, "All wells");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 18);
        if (ImGui::BeginCombo("Preset", "Insert a preset...")) {
            for (int i = 0; i < int(presets.size()); ++i)
                if (ImGui::Selectable(presets[size_t(i)].name.c_str())) {
                    std::snprintf(G.expr, sizeof G.expr, "%s", presets[size_t(i)].expr.c_str());
                    std::snprintf(G.mnem, sizeof G.mnem, "%s", presets[size_t(i)].mnem.c_str());
                    std::snprintf(G.unit, sizeof G.unit, "%s", presets[size_t(i)].unit.c_str());
                } else if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s\n%s", presets[size_t(i)].expr.c_str(), presets[size_t(i)].help.c_str());
            ImGui::EndCombo();
        }
        ImGui::SetNextItemWidth(fs * 6); ImGui::InputText("New curve", G.mnem, sizeof G.mnem, ImGuiInputTextFlags_CharsUppercase);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputText("Unit", G.unit, sizeof G.unit);
        ImGui::SameLine(); ImGui::TextUnformatted("=");
        ImGui::SameLine(); ImGui::SetNextItemWidth(-1); ImGui::InputText("##expr", G.expr, sizeof G.expr);
        // clickable names / functions append to the expression
        ImGui::TextDisabled("Curves:"); 
        const int pw = std::max(0, G.calcWell);
        std::vector<std::string> avail = S.wells[size_t(pw)].logs.empty() ? std::vector<std::string>{"MD", "TVD"} : geo::expressionNames(S.wells[size_t(pw)].logs[0]);
        for (const auto& n : avail) {
            ImGui::SameLine();
            if (ImGui::SmallButton(n.c_str())) { const size_t L = std::strlen(G.expr); std::snprintf(G.expr + L, sizeof G.expr - L, "%s%s", L && G.expr[L - 1] != ' ' && G.expr[L - 1] != '(' ? " " : "", n.c_str()); }
        }
        ImGui::TextDisabled("Functions: abs sqrt ln log10 exp pow(a,b) min max clamp(x,lo,hi) if(c,a,b) isnan  |  operators + - * / ^ < <= > >= == != && || !");
        // live preview on one well
        const std::string key = std::string(G.expr) + "|" + std::to_string(pw) + "|" + std::to_string(G.rev);
        if (key != G.prevKey) {
            G.prevKey = key; G.prevX.clear(); G.prevZ.clear(); G.prevErr.clear();
            const auto& w = S.wells[size_t(pw)];
            if (!w.logs.empty()) {
                geo::WellLog tmp = w.logs[0];
                std::string err;
                if (geo::addCalculatedCurve(tmp, &w, "__PREVIEW", "", G.expr, err)) {
                    const geo::LogCurve* c = tmp.curve("__PREVIEW");
                    for (size_t s = 0; c && s < tmp.md.size(); ++s) if (!isNull(c->v[s], tmp.nullValue)) { G.prevX.push_back(c->v[s]); G.prevZ.push_back(w.at(tmp.md[s])[2]); }
                } else G.prevErr = err;
            }
        }
        if (!G.prevErr.empty()) ImGui::TextColored(c4f(Pal::Err), ICON_TRIANGLE_ALERT " %s", G.prevErr.c_str());
        else {
            double mn = 1e300, mx = -1e300, sum = 0;
            for (double v : G.prevX) { mn = std::min(mn, v); mx = std::max(mx, v); sum += v; }
            if (G.prevX.empty()) ImGui::TextDisabled("Preview on %s: no valid samples", wellNames[size_t(pw)].c_str());
            else ImGui::TextDisabled("Preview on %s: %zu samples, min %.4g, mean %.4g, max %.4g", wellNames[size_t(pw)].c_str(), G.prevX.size(), mn, sum / G.prevX.size(), mx);
        }
        ImGui::BeginDisabled(busy || !G.prevErr.empty() || !G.mnem[0]);
        if (ui::accentButton(ICON_CHECK " Compute curve")) {
            int ok = 0; std::string fails;
            for (int wi = 0; wi < int(S.wells.size()); ++wi) {
                if (G.calcWell >= 0 && wi != G.calcWell) continue;
                auto& w = S.wells[size_t(wi)];
                for (auto& l : w.logs) {
                    std::string err;
                    if (geo::addCalculatedCurve(l, &w, G.mnem, G.unit, G.expr, err)) ++ok;
                    else fails += "  " + w.name + ": " + err + "\n";
                }
            }
            ++G.rev;
            const std::string b = "Log calculator: " + std::string(G.mnem) + " = " + G.expr + " computed in " + std::to_string(ok) + " log(s)";
            S.msg(b);
            G.calcMsg = b + (fails.empty() ? std::string() : "\nNot computed:\n" + fails);
        }
        ImGui::EndDisabled();
        ImGui::SameLine(); ImGui::TextDisabled("Results are new curves of the logs: they appear in the well correlation, blocking, cross-plots and variograms.");
        if (!G.calcMsg.empty()) ImGui::TextWrapped("%s", G.calcMsg.c_str());
        // calculated curves of the preview well (removable)
        if (!S.wells[size_t(pw)].logs.empty()) {
            std::string del;
            for (const auto& c : S.wells[size_t(pw)].logs[0].curves) {
                if (c.descr.rfind("= ", 0) != 0) continue;
                ImGui::PushID(c.mnem.c_str());
                if (ImGui::SmallButton(ICON_X)) del = c.mnem;
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Remove %s from every well", c.mnem.c_str());
                ImGui::SameLine(); ImGui::Text("%s [%s] %s", c.mnem.c_str(), c.unit.c_str(), c.descr.c_str());
                ImGui::PopID();
            }
            if (!del.empty() && !busy) {
                for (auto& w : S.wells) for (auto& l : w.logs)
                    l.curves.erase(std::remove_if(l.curves.begin(), l.curves.end(), [&](const geo::LogCurve& c) { return c.mnem == del && c.descr.rfind("= ", 0) == 0; }), l.curves.end());
                ++G.rev; S.msg("Log calculator: removed " + del);
            }
        }
        const float h = std::max(fs * 10, ImGui::GetContentRegionAvail().y);
        if (!G.prevX.empty() && ImPlot::BeginPlot("##calcprev", ImVec2(-1, h))) {
            const std::string t = std::string(G.mnem) + " (preview) - " + wellNames[size_t(pw)];
            recPanel(t.c_str(), G.mnem, "TVD", "");
            ImPlot::SetupAxes(G.mnem, "TVD", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
            ImPlot::PlotLine(G.mnem, G.prevX.data(), G.prevZ.data(), int(G.prevX.size()), lineSpec(c4f(seriesColor(0)), 1.4f));
            recSeries(G.mnem, G.prevX.data(), G.prevZ.data(), int(G.prevX.size()), seriesColor(0), 1.4f);
            for (const auto& tp : S.wells[size_t(pw)].tops) {
                const double z = S.wells[size_t(pw)].at(tp.md)[2];
                ImPlot::TagY(z, c4f(seriesColor(3)), "%s", tp.surface.c_str());
                double zz[1] = {z};
                ImPlotSpec hs = lineSpec(c4f(seriesColor(3)), 1.0f); hs.Flags = ImPlotInfLinesFlags_Horizontal;
                ImPlot::PlotInfLines("##top", zz, 1, hs);
            }
            recEnd();
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ================================================================ cross-plot
    if (ui::tabItem(ICON_CHART_SCATTER " Cross-plot", nullptr, tabFlags(1))) {
        if (G.cx < 0) G.cx = findName(names, "PHIE", 0);
        if (G.cy < 0) G.cy = findName(names, "PERM", std::min(1, int(names.size()) - 1));
        ImGui::SetNextItemWidth(fs * 7); combo("X##xx", &G.cx, names);
        ImGui::SameLine(); ImGui::Checkbox(TR("log##lx"), &G.logX);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); combo("Y##yy", &G.cy, names);
        ImGui::SameLine(); ImGui::Checkbox(TR("log##ly"), &G.logY);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
        static const char* kCol[] = {"Single colour", "By well", "By curve"};
        ImGui::Combo("Colour", &G.colorMode, kCol, 3);
        if (G.colorMode == 2) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); if (G.cz < 0) G.cz = findName(names, "GR", 0); combo("##zz", &G.cz, names); }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9);
        static const char* kReg[] = {"No regression", "Linear y = a + b x", "Log-linear log10 y = a + b x"};
        ImGui::Combo("##reg", &G.reg, kReg, 3);
        zoneFilterUi(tops, fs);
        const int idx[3] = {G.cx, G.cy, G.colorMode == 2 ? G.cz : -1};
        const auto smp = gather(S, names, idx, 3, tops);
        std::vector<double> X, Y, Z; std::vector<int> W;
        for (const auto& q : smp) {
            if (std::isnan(q.v[0]) || std::isnan(q.v[1])) continue;
            if ((G.logX && q.v[0] <= 0) || (G.logY && q.v[1] <= 0)) continue;
            if (G.colorMode == 2 && std::isnan(q.v[2])) continue;
            X.push_back(q.v[0]); Y.push_back(q.v[1]); Z.push_back(q.v[2]); W.push_back(q.well);
        }
        geo::Regression R;
        if (G.reg > 0) R = geo::regress(X, Y, G.reg == 2);
        if (G.reg > 0 && R.n >= 2) {
            ImGui::Text("%s  a = %.4g, b = %.4g, R\xc2\xb2 = %.3f, n = %d", G.reg == 2 ? "log10 y = a + b x:" : "y = a + b x:", R.a, R.b, R.r2, R.n);
            if (G.reg == 2) {
                ImGui::SameLine();
                ImGui::BeginDisabled(busy);
                if (ImGui::SmallButton(TR(ICON_CHECK " Use as the poro-perm transform"))) {
                    S.userFit = geo::PoroPerm{R.a, R.b, R.r2, R.n}; S.useUserFit = true;
                    char b[200]; std::snprintf(b, sizeof b, "Poro-perm transform from the cross-plot: log10 k = %.4f + %.4f phi (R2 %.3f, n %d) - used by the Petrophysics step", R.a, R.b, R.r2, R.n);
                    S.msg(b);
                }
                ImGui::EndDisabled();
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("X must be the porosity curve the Petrophysics step models (fraction) and Y the permeability (mD).");
            }
        }
        double zlo = 1e300, zhi = -1e300;
        if (G.colorMode == 2) for (double z : Z) { zlo = std::min(zlo, z); zhi = std::max(zhi, z); }
        if (G.colorMode == 2 && zhi > zlo) {
            ImGui::SameLine(); ImGui::Text("   %s %.3g", names[size_t(G.cz)].c_str(), zlo); ImGui::SameLine(); gradientBar(Colormap::Viridis, fs * 7, fs * 0.8f); ImGui::SameLine(); ImGui::Text("%.3g", zhi);
        }
        const float h = std::max(fs * 10, ImGui::GetContentRegionAvail().y);
        if (!X.empty() && ImPlot::BeginPlot("##xplot", ImVec2(-1, h))) {
            const std::string xl = names[size_t(G.cx)], yl = names[size_t(G.cy)];
            recPanel(("Cross-plot " + yl + " vs " + xl).c_str(), xl.c_str(), yl.c_str(), "");
            ImPlot::SetupAxes(xl.c_str(), yl.c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (G.logX) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            if (G.logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (G.colorMode == 1) {
                for (int wi = 0; wi < int(S.wells.size()); ++wi) {
                    std::vector<double> xs, ys;
                    for (size_t k = 0; k < X.size(); ++k) if (W[k] == wi) { xs.push_back(X[k]); ys.push_back(Y[k]); }
                    if (xs.empty()) continue;
                    ImPlot::PlotScatter(wellNames[size_t(wi)].c_str(), xs.data(), ys.data(), int(xs.size()), dotSpec(c4f(seriesColor(wi)), 2.2f));
                    recSeries(wellNames[size_t(wi)].c_str(), xs.data(), ys.data(), int(xs.size()), seriesColor(wi), 1.0f, true);
                }
            } else if (G.colorMode == 2 && zhi > zlo) {
                std::vector<ImU32> cols(X.size());
                for (size_t k = 0; k < X.size(); ++k) { float c[3]; colormapRGB(Colormap::Viridis, float((Z[k] - zlo) / (zhi - zlo)), c); cols[k] = IM_COL32(int(c[0] * 255), int(c[1] * 255), int(c[2] * 255), 255); }
                ImPlotSpec s = dotSpec(c4f(seriesColor(0)), 2.2f); s.MarkerFillColors = cols.data(); s.MarkerLineColors = cols.data();
                ImPlot::PlotScatter("##pts", X.data(), Y.data(), int(X.size()), s);
                recSeries("samples", X.data(), Y.data(), int(X.size()), seriesColor(0), 1.0f, true);
            } else {
                ImPlot::PlotScatter("samples", X.data(), Y.data(), int(X.size()), dotSpec(c4f(seriesColor(0)), 2.2f));
                recSeries("samples", X.data(), Y.data(), int(X.size()), seriesColor(0), 1.0f, true);
            }
            if (G.reg > 0 && R.n >= 2) {
                const double x0 = *std::min_element(X.begin(), X.end()), x1 = *std::max_element(X.begin(), X.end());
                std::vector<double> lx(60), ly(60);
                for (int k = 0; k < 60; ++k) {
                    lx[size_t(k)] = G.logX ? x0 * std::pow(x1 / x0, k / 59.0) : x0 + (x1 - x0) * k / 59.0;
                    const double yy = R.a + R.b * lx[size_t(k)];
                    ly[size_t(k)] = G.reg == 2 ? std::pow(10.0, yy) : yy;
                }
                ImPlot::PlotLine("regression", lx.data(), ly.data(), 60, lineSpec(c4f(Pal::Err), 2.0f));
                recSeries("regression", lx.data(), ly.data(), 60, Pal::Err, 2.0f);
            }
            recEnd();
            ImPlot::EndPlot();
        } else if (X.empty()) ImGui::TextDisabled("No samples with both curves in the zone.");
        ImGui::EndTabItem();
    }

    // ================================================================ histogram
    if (ui::tabItem(ICON_CHART_COLUMN " Histogram", nullptr, tabFlags(2))) {
        if (G.hc < 0) G.hc = findName(names, "PHIE", 0);
        ImGui::SetNextItemWidth(fs * 7); combo("Curve##hc", &G.hc, names);
        ImGui::SameLine(); ImGui::Checkbox(TR("log10"), &G.hLog);
        ImGui::SameLine(); ImGui::Checkbox(TR("per well"), &G.hByWell);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); ImGui::SliderInt("bins", &G.bins, 5, 100);
        zoneFilterUi(tops, fs);
        const int idx[1] = {G.hc};
        const auto smp = gather(S, names, idx, 1, tops);
        std::vector<double> V; std::vector<int> W;
        for (const auto& q : smp) { if (std::isnan(q.v[0]) || (G.hLog && q.v[0] <= 0)) continue; V.push_back(G.hLog ? std::log10(q.v[0]) : q.v[0]); W.push_back(q.well); }
        if (!V.empty()) {
            double m = 0, sd = 0;
            for (double v : V) m += v;
            m /= double(V.size());
            for (double v : V) sd += (v - m) * (v - m);
            sd = std::sqrt(sd / double(V.size()));
            ImGui::Text("n %zu   mean %.4g   sd %.4g   P10 %.4g   P50 %.4g   P90 %.4g   min %.4g   max %.4g", V.size(), m, sd, quantile(V, 0.1), quantile(V, 0.5), quantile(V, 0.9),
                        *std::min_element(V.begin(), V.end()), *std::max_element(V.begin(), V.end()));
            const double lo = *std::min_element(V.begin(), V.end()), hi = *std::max_element(V.begin(), V.end());
            const float h = std::max(fs * 10, ImGui::GetContentRegionAvail().y);
            if (ImPlot::BeginPlot("##hist", ImVec2(-1, h))) {
                const std::string xl = (G.hLog ? "log10 " : "") + names[size_t(G.hc)];
                recPanel(("Histogram " + xl).c_str(), xl.c_str(), "count", "");
                ImPlot::SetupAxes(xl.c_str(), "count", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                const ImPlotRange rg(lo, hi > lo ? hi : lo + 1);
                auto draw = [&](const char* lab, const std::vector<double>& v, ImU32 c, float alpha) {
                    ImPlotSpec s; s.FillColor = c4f(c); s.LineColor = c4f(c); s.FillAlpha = alpha;
                    ImPlot::PlotHistogram(lab, v.data(), int(v.size()), G.bins, 1.0, rg, s);
                    // export: bar outline as a step polyline
                    std::vector<double> cnt(size_t(G.bins), 0.0), bx, by;
                    const double bw = (rg.Max - rg.Min) / G.bins;
                    for (double x : v) cnt[size_t(std::clamp(int((x - rg.Min) / bw), 0, G.bins - 1))] += 1;
                    for (int b = 0; b < G.bins; ++b) { bx.push_back(rg.Min + b * bw); by.push_back(cnt[size_t(b)]); bx.push_back(rg.Min + (b + 1) * bw); by.push_back(cnt[size_t(b)]); }
                    recSeries(lab, bx.data(), by.data(), int(bx.size()), c, 1.5f);
                };
                if (G.hByWell) {
                    for (int wi = 0; wi < int(S.wells.size()); ++wi) {
                        std::vector<double> v;
                        for (size_t k = 0; k < V.size(); ++k) if (W[k] == wi) v.push_back(V[k]);
                        if (!v.empty()) draw(wellNames[size_t(wi)].c_str(), v, seriesColor(wi), 0.45f);
                    }
                } else draw(names[size_t(G.hc)].c_str(), V, seriesColor(0), 0.75f);
                recEnd();
                ImPlot::EndPlot();
            }
        } else ImGui::TextDisabled("No samples of this curve in the zone.");
        ImGui::EndTabItem();
    }

    // ================================================================ variogram
    if (ui::tabItem(ICON_ACTIVITY " Variogram", nullptr, tabFlags(3))) {
        if (G.vc < 0) G.vc = findName(names, "PHIE", 0);
        ImGui::SetNextItemWidth(fs * 7); combo("Curve##vc", &G.vc, names);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8);
        static const char* kTr[] = {"Values", "log10", "Indicator ="};
        ImGui::Combo("##vtr", &G.vTrans, kTr, 3);
        if (G.vTrans == 2) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4); ImGui::InputInt("##vcode", &G.vCode, 0); }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9);
        static const char* kDir[] = {"Horizontal (omni)", "Vertical", "Along azimuth"};
        int dirUi = G.vDir == 0 ? 0 : G.vDir == 1 ? 1 : 2;
        if (ImGui::Combo("##vdir", &dirUi, kDir, 3)) G.vDir = dirUi;
        if (G.vDir == 2) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("az##vaz", &G.vAz, 0, 0, "%.0f"); }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 4); ImGui::InputInt("lags", &G.vLags, 0); G.vLags = std::clamp(G.vLags, 4, 80);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("lag (0 auto)", &G.vLag, 0, 0, "%.3g");
        zoneFilterUi(tops, fs);
        ImGui::SameLine(); ImGui::Checkbox(TR("Standardise (sill 1)"), &G.vStd);
        ImGui::SameLine(); ImGui::BeginDisabled(G.zoneTop < 0); ImGui::Checkbox(TR("Stratigraphic depth"), &G.vStrat); ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Depth below the zone top of each well (flattens dipping layers for horizontal pairs; needs a zone top).");
        // experimental variogram (recomputed when an input changes)
        char kb[256];
        std::snprintf(kb, sizeof kb, "%d|%d|%d|%d|%.6g|%d|%.6g|%d|%d|%d|%d|%d|%zu", G.vc, G.vTrans, G.vCode, G.vDir, G.vAz, G.vLags, G.vLag, G.zoneTop, G.zoneBase,
                      int(G.vStd), int(G.vStrat), G.rev, S.wells.size());
        if (kb != G.evKey) {
            G.evKey = kb; G.ev = geostat::ExpVario(); G.evMsg.clear();
            const int idx[1] = {G.vc};
            const auto smp = gather(S, names, idx, 1, tops);
            std::vector<geostat::Point> pts;
            for (const auto& q : smp) {
                double v = q.v[0];
                if (std::isnan(v)) continue;
                if (G.vTrans == 1) { if (v <= 0) continue; v = std::log10(v); }
                else if (G.vTrans == 2) v = std::lround(v) == G.vCode ? 1.0 : 0.0;
                pts.push_back({q.x, q.y, G.vStrat && G.zoneTop >= 0 ? q.zs : q.z, v});
            }
            G.evN = int(pts.size());
            if (pts.size() >= 10) {
                double m = 0, var = 0;
                for (const auto& p : pts) m += p.v;
                m /= double(pts.size());
                for (const auto& p : pts) var += (p.v - m) * (p.v - m);
                var /= double(pts.size());
                double lw = G.vLag;
                if (lw <= 0) {
                    if (G.vDir == 1) {   // vertical: zone thickness / 2 / lags (median over wells of the sampled thickness)
                        std::map<int, std::pair<double, double>> zr;
                        for (const auto& q : smp) { auto& r = zr.emplace(q.well, std::make_pair(1e300, -1e300)).first->second; r.first = std::min(r.first, q.z); r.second = std::max(r.second, q.z); }
                        std::vector<double> th; for (const auto& kv : zr) th.push_back(kv.second.second - kv.second.first);
                        lw = std::max(1e-3, quantile(th, 0.5) / (2.0 * G.vLags));
                    } else {             // horizontal: largest well-to-well distance / 2 / lags
                        double dmax = 0;
                        for (size_t a = 0; a < S.wells.size(); ++a) for (size_t b = a + 1; b < S.wells.size(); ++b) dmax = std::max(dmax, std::hypot(S.wells[a].x - S.wells[b].x, S.wells[a].y - S.wells[b].y));
                        lw = std::max(1e-3, dmax / (1.5 * G.vLags));
                    }
                }
                G.ev = geostat::experimental(pts, lw, G.vLags, G.vDir, G.vAz);
                G.evVar = var; G.evLagW = lw;
                if (G.vStd && var > 0) for (double& g : G.ev.gamma) g /= var;
                if (!G.haveModel) {   // initial model: total sill = variance, range = half the lag span
                    G.sill = G.vStd ? 1.0 : var; G.nug = 0; G.range = lw * G.vLags * 0.5; G.haveModel = true;
                }
            } else G.evMsg = "fewer than 10 samples in the zone";
        }
        ImGui::SetNextItemWidth(fs * 8);
        static const char* kMod[] = {"Spherical", "Exponential", "Gaussian"};
        ImGui::Combo("Model", &G.vModel, kMod, 3);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("nugget", &G.nug, 0, 0, "%.3g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("sill", &G.sill, 0, 0, "%.3g");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::InputDouble("range", &G.range, 0, 0, "%.4g");
        G.nug = std::max(0.0, G.nug); G.sill = std::max(0.0, G.sill); G.range = std::max(1e-6, G.range);
        ImGui::SameLine(); ImGui::Checkbox(TR("fix nugget"), &G.fixNug);
        ImGui::SameLine();
        if ((ImGui::Button(TR(ICON_WAND_SPARKLES " Auto fit")) || G.smokeFit) && !G.ev.lag.empty()) {
            G.smokeFit = false;
            const auto f = geo::fitVariogram(G.ev, geostat::VarioModel(G.vModel), G.fixNug ? G.nug : -1.0);
            if (f.ok) { G.nug = f.nugget; G.sill = f.sill; G.range = f.range; G.rmse = f.rmse; }
            else G.evMsg = "the fit failed (flat or too few lags)";
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Weighted least squares (Cressie weights N(h) / gamma(h)^2) over the range, nugget and sill >= 0.");
        if (G.rmse > 0) { ImGui::SameLine(); ImGui::TextDisabled("rms %.3g", G.rmse); }
        ImGui::SetNextItemWidth(fs * 9);
        static const char* kTo[] = {"Facies variogram", "Petrophysics variogram"};
        ImGui::Combo("##vto", &G.applyTo, kTo, 2);
        ImGui::SameLine();
        ImGui::BeginDisabled(busy || G.ev.lag.empty());
        if (ImGui::Button(TR(ICON_CHECK " Apply"))) {
            geostat::Variogram& vg = G.applyTo == 0 ? S.fs.vg : S.ps.vg;
            const double tot = std::max(1e-30, G.nug + G.sill);
            vg.model = geostat::VarioModel(G.vModel); vg.nugget = G.nug / tot; vg.sill = G.sill / tot;   // the simulations use a unit total sill
            if (G.vDir == 1) vg.rangeVert = G.range;
            else if (G.vDir == 0) vg.rangeMax = vg.rangeMin = G.range;
            else { vg.rangeMax = G.range; vg.azimuthDeg = G.vAz; }
            char b[220];
            std::snprintf(b, sizeof b, "Variogram fit applied to the %s: %s, nugget %.3f, %s range %.4g", G.applyTo == 0 ? "Facies step" : "Petrophysics step", kMod[G.vModel], vg.nugget,
                          G.vDir == 1 ? "vertical" : G.vDir == 0 ? "horizontal" : "major", G.range);
            S.msg(b);
            G.evMsg = b;
        }
        ImGui::EndDisabled();
        ImGui::SameLine(); ImGui::TextDisabled("%d samples, lag %.3g%s%s", G.evN, G.evLagW, G.evMsg.empty() ? "" : "  |  ", G.evMsg.c_str());
        const float h = std::max(fs * 10, ImGui::GetContentRegionAvail().y);
        if (!G.ev.lag.empty() && ImPlot::BeginPlot("##vario", ImVec2(-1, h))) {
            const char* yl = G.vStd ? "gamma / variance" : "gamma";
            recPanel(("Variogram " + names[size_t(std::clamp(G.vc, 0, int(names.size()) - 1))] + " (" + kDir[dirUi] + ")").c_str(), "lag distance", yl, "");
            ImPlot::SetupAxes("lag distance", yl, 0, 0);
            double gmax = 0, hmax = 0;
            for (size_t k = 0; k < G.ev.lag.size(); ++k) if (G.ev.pairs[k] > 0) { gmax = std::max(gmax, G.ev.gamma[k]); hmax = std::max(hmax, G.ev.lag[k]); }
            ImPlot::SetupAxesLimits(0, std::max(hmax, G.range) * 1.1, 0, std::max(gmax, G.nug + G.sill) * 1.25, ImPlotCond_Once);
            std::vector<double> ex, ey;
            for (size_t k = 0; k < G.ev.lag.size(); ++k) if (G.ev.pairs[k] > 0) { ex.push_back(G.ev.lag[k]); ey.push_back(G.ev.gamma[k]); }
            ImPlot::PlotScatter("experimental", ex.data(), ey.data(), int(ex.size()), dotSpec(c4f(seriesColor(1)), 4.0f));
            recSeries("experimental", ex.data(), ey.data(), int(ex.size()), seriesColor(1), 1.0f, true, true);
            if (ImPlot::IsPlotHovered()) {   // pairs per lag on hover
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                for (size_t k = 0; k < G.ev.lag.size(); ++k)
                    if (G.ev.pairs[k] > 0 && std::fabs(G.ev.lag[k] - mp.x) < 0.5 * G.evLagW) { ImGui::SetTooltip("lag %.4g\ngamma %.4g\n%ld pairs", G.ev.lag[k], G.ev.gamma[k], G.ev.pairs[k]); break; }
            }
            if (G.vStd) { double one[1] = {1.0}; ImPlotSpec rs = lineSpec(c4f(ImGui::GetColorU32(ImGuiCol_TextDisabled)), 1.0f); rs.Flags = ImPlotInfLinesFlags_Horizontal; ImPlot::PlotInfLines("##sill1", one, 1, rs); }
            const double hEnd = std::max(hmax, G.range) * 1.1;
            std::vector<double> mx(121), my(121);
            for (int k = 0; k <= 120; ++k) { mx[size_t(k)] = hEnd * k / 120.0; my[size_t(k)] = k == 0 ? 0.0 : G.nug + G.sill * geo::structure(geostat::VarioModel(G.vModel), mx[size_t(k)] / G.range); }
            ImPlot::PlotLine("model", mx.data(), my.data(), 121, lineSpec(c4f(Pal::Accent), 2.0f));
            recSeries("model", mx.data(), my.data(), 121, Pal::Accent, 2.0f);
            {   // currently applied variogram of the target step (dashed, for comparison)
                const geostat::Variogram& vg = G.applyTo == 0 ? S.fs.vg : S.ps.vg;
                const double r = G.vDir == 1 ? vg.rangeVert : vg.rangeMax, sc = G.vStd ? 1.0 : G.evVar;
                std::vector<double> cy(121);
                for (int k = 0; k <= 120; ++k) cy[size_t(k)] = k == 0 ? 0.0 : sc * (vg.nugget + vg.sill * geo::structure(vg.model, mx[size_t(k)] / std::max(1e-9, r)));
                ImPlotSpec ds = lineSpec(c4f(seriesColor(2)), 1.2f); ds.Flags = ImPlotLineFlags_Segments;
                std::vector<double> sx, sy;   // dashes: every other segment
                for (int k = 0; k + 1 <= 120; k += 2) { sx.push_back(mx[size_t(k)]); sy.push_back(cy[size_t(k)]); sx.push_back(mx[size_t(k) + 1]); sy.push_back(cy[size_t(k) + 1]); }
                ImPlot::PlotLine(G.applyTo == 0 ? "current facies" : "current petrophysics", sx.data(), sy.data(), int(sx.size()), ds);
            }
            // draggable handles: nugget at h = 0, (range, nugget + sill)
            double hx = G.range, hy = G.nug + G.sill;
            if (ImPlot::DragPoint(0x7C41, &hx, &hy, c4f(Pal::Accent), 7)) { G.range = std::max(1e-6, hx); G.sill = std::max(0.0, hy - G.nug); G.rmse = 0; }
            double nx = 0, ny = G.nug;
            if (ImPlot::DragPoint(0x7C42, &nx, &ny, c4f(Pal::Warn), 6)) { const double tot = G.nug + G.sill; G.nug = std::clamp(ny, 0.0, tot); G.sill = tot - G.nug; G.rmse = 0; }
            ImPlot::TagX(G.range, c4f(Pal::Accent), "range %.4g", G.range);
            recEnd();
            ImPlot::EndPlot();
        } else if (G.ev.lag.empty()) ImGui::TextDisabled("No experimental variogram: %s", G.evMsg.empty() ? "choose a curve with samples" : G.evMsg.c_str());
        ImGui::EndTabItem();
    }

    // ================================================================ well tie
    if (ui::tabItem(ICON_WAVES " Well tie", nullptr, tabFlags(4))) {
        G.tw = std::clamp(G.tw, 0, int(S.wells.size()) - 1);
        if (G.tdt < 0) G.tdt = findName(names, "DT", findName(names, "DTC", -1));
        if (G.trho < 0) G.trho = findName(names, "RHOB", -1);
        ImGui::SetNextItemWidth(fs * 7); combo("Well##tw", &G.tw, wellNames);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); combo("Sonic##tdt", &G.tdt, names, "(none)");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); combo("Density##trho", &G.trho, names, "(constant)");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6);
        static const char* kU[] = {"auto", "us/ft", "us/m"};
        ImGui::Combo("unit##tdu", &G.ti.dtUnit, kU, 3);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); sliderD("Ricker Hz", &G.ti.freq, 5, 80, "%.0f");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); sliderD("max shift ms", &G.ti.maxShift, 0, 200, "%.0f");
        G.ti.dt = G.tdt >= 0 ? names[size_t(G.tdt)] : std::string();
        G.ti.rhob = G.trho >= 0 ? names[size_t(G.trho)] : std::string();
        G.ti.field = S.field;
        if (!S.hasCube) ImGui::TextDisabled("Load a seismic volume (SEG-Y) on the Data / Seismic step to tie the wells.");
        else if (G.tdt < 0) ImGui::TextDisabled("No sonic curve: load DT (or compute it with the log calculator, e.g. DT = 304800 / VP for VP in m/s).");
        else {
            char kb[256];
            std::snprintf(kb, sizeof kb, "%d|%s|%s|%d|%.3f|%.3f|%.4f|%.4f|%.4f|%d|%d", G.tw, G.ti.dt.c_str(), G.ti.rhob.c_str(), G.ti.dtUnit, G.ti.freq, G.ti.maxShift,
                          S.vm.v0, S.vm.k, S.vm.datum, G.rev, S.cube.ns);
            if (kb != G.tKey) {
                G.tKey = kb;
                const auto& w = S.wells[size_t(G.tw)];
                const geo::WellLog* lg = nullptr;
                for (const auto& l : w.logs) if (l.curve(G.ti.dt)) { lg = &l; break; }
                if (lg) G.tr = geo::wellTie(w, *lg, S.cube, S.vm, G.ti);
                else { G.tr = geo::TieResult(); G.tr.msg = w.name + " has no " + G.ti.dt + " curve"; }
            }
            const geo::TieResult& T = G.tr;
            if (!T.ok) ImGui::TextColored(c4f(Pal::Warn), ICON_TRIANGLE_ALERT " %s", T.msg.c_str());
            else {
                ImGui::TextColored(c4f(T.corr >= 0.6 ? Pal::Ok : T.corr >= 0.3 ? Pal::Warn : Pal::Err), "%s %s", T.corr >= 0.6 ? ICON_CHECK : ICON_TRIANGLE_ALERT, T.msg.c_str());
                ImGui::BeginDisabled(busy);
                if (ImGui::Button(TR(ICON_CHECK " Apply the calibrated velocity model"))) {
                    S.vm = T.vmFit;
                    char b[200]; std::snprintf(b, sizeof b, "Well tie %s: velocity model set to v(z) = %.1f + %.4f (z - %.1f); re-run the depth conversion on the Seismic step",
                                               wellNames[size_t(G.tw)].c_str(), S.vm.v0, S.vm.k, S.vm.datum);
                    S.msg(b);
                }
                ImGui::EndDisabled();
                ImGui::SameLine(); ImGui::TextDisabled("current v(z) = %.1f + %.4f (z - %.1f)", S.vm.v0, S.vm.k, S.vm.datum);
                // tops in time from the tie's time-depth pairs
                std::vector<std::pair<std::string, double>> topsT;
                const auto& w = S.wells[size_t(G.tw)];
                for (const auto& tp : w.tops) {
                    const double z = w.at(tp.md)[2];
                    double t = NAN;
                    for (size_t k = 0; k + 1 < T.tdZ.size(); ++k)
                        if (z >= T.tdZ[k] && z <= T.tdZ[k + 1]) { t = T.tdT[k] + (z - T.tdZ[k]) / std::max(1e-12, T.tdZ[k + 1] - T.tdZ[k]) * (T.tdT[k + 1] - T.tdT[k]); break; }
                    if (std::isfinite(t)) topsT.push_back({tp.surface, t});
                }
                // normalise the synthetic to the RMS of the trace inside the logged window
                double rs = 0, rt = 0; int nn = 0;
                for (size_t s = 0; s < T.traceT.size(); ++s) {
                    if (T.traceT[s] < T.twt.front() || T.traceT[s] > T.twt.back()) continue;
                    rt += T.trace[s] * T.trace[s]; ++nn;
                }
                for (double v : T.synth) rs += v * v;
                const double scale = nn && rs > 0 ? std::sqrt(rt / nn) / std::sqrt(rs / double(T.synth.size())) : 1.0;
                std::vector<double> synS(T.synth.size());
                for (size_t k = 0; k < synS.size(); ++k) synS[k] = T.synth[k] * scale;
                const float h = std::max(fs * 12, ImGui::GetContentRegionAvail().y);
                const double tLo = T.twt.front() - 30, tHi = T.twt.back() + 30;
                if (ImPlot::BeginSubplots("##tie", 1, 4, ImVec2(-1, h), ImPlotSubplotFlags_LinkRows)) {
                    auto tracks = [&](const char* title, const char* xl) {
                        if (!ImPlot::BeginPlot(title)) return false;
                        recPanel(title, xl, "TWT (ms)", "");
                        ImPlot::SetupAxes(xl, "TWT (ms)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_Invert);
                        ImPlot::SetupAxisLimits(ImAxis_Y1, tLo, tHi, ImPlotCond_Once);
                        return true;
                    };
                    auto topLines = [&]() {
                        for (const auto& tt : topsT) {
                            double v[1] = {tt.second};
                            ImPlotSpec hs = lineSpec(c4f(seriesColor(3)), 1.0f); hs.Flags = ImPlotInfLinesFlags_Horizontal;
                            ImPlot::PlotInfLines("##tops", v, 1, hs);
                            ImPlot::TagY(tt.second, c4f(seriesColor(3)), "%s", tt.first.c_str());
                        }
                    };
                    const int n = int(T.twt.size());
                    if (tracks("Acoustic impedance", "AI (g/cc m/s)")) {
                        ImPlot::PlotLine("AI", T.ai.data(), T.twt.data(), n, lineSpec(c4f(seriesColor(0)), 1.3f));
                        recSeries("AI", T.ai.data(), T.twt.data(), n, seriesColor(0), 1.3f);
                        topLines(); recEnd(); ImPlot::EndPlot();
                    }
                    if (tracks("Reflectivity", "RC")) {
                        std::vector<double> sx, sy;   // stems
                        for (int k = 0; k < n; ++k) if (std::fabs(T.rc[size_t(k)]) > 1e-4) { sx.push_back(0); sy.push_back(T.twt[size_t(k)]); sx.push_back(T.rc[size_t(k)]); sy.push_back(T.twt[size_t(k)]); }
                        ImPlotSpec ss = lineSpec(c4f(seriesColor(1)), 1.2f); ss.Flags = ImPlotLineFlags_Segments;
                        ImPlot::PlotLine("RC", sx.data(), sy.data(), int(sx.size()), ss);
                        recSeries("RC", T.rc.data(), T.twt.data(), n, seriesColor(1), 1.0f);
                        topLines(); recEnd(); ImPlot::EndPlot();
                    }
                    if (tracks("Synthetic vs seismic", "amplitude")) {
                        ImPlot::PlotLine("seismic", T.trace.data(), T.traceT.data(), int(T.trace.size()), lineSpec(c4f(Pal::Err), 1.6f));
                        ImPlot::PlotLine("synthetic", synS.data(), T.twt.data(), n, lineSpec(c4f(seriesColor(0)), 1.6f));
                        recSeries("seismic", T.trace.data(), T.traceT.data(), int(T.trace.size()), Pal::Err, 1.6f);
                        recSeries("synthetic", synS.data(), T.twt.data(), n, seriesColor(0), 1.6f);
                        topLines(); recEnd(); ImPlot::EndPlot();
                    }
                    if (tracks("Time-depth", "TVD")) {
                        ImPlot::PlotLine("sonic (tied)", T.tdZ.data(), T.tdT.data(), int(T.tdZ.size()), lineSpec(c4f(seriesColor(0)), 2.0f));
                        recSeries("sonic (tied)", T.tdZ.data(), T.tdT.data(), int(T.tdZ.size()), seriesColor(0), 2.0f);
                        std::vector<double> zc, tc, tf;
                        if (!T.tdZ.empty()) {
                            const double z0 = T.tdZ.front(), z1 = T.tdZ.back();
                            for (int k = 0; k <= 40; ++k) { const double z = z0 + (z1 - z0) * k / 40.0; zc.push_back(z); tc.push_back(geo::depthToTime(S.vm, z)); tf.push_back(geo::depthToTime(T.vmFit, z)); }
                            ImPlot::PlotLine("current v(z)", zc.data(), tc.data(), int(zc.size()), lineSpec(c4f(Pal::Warn), 1.3f));
                            ImPlot::PlotLine("calibrated v(z)", zc.data(), tf.data(), int(zc.size()), lineSpec(c4f(Pal::Accent), 1.3f));
                            recSeries("current v(z)", zc.data(), tc.data(), int(zc.size()), Pal::Warn, 1.3f);
                            recSeries("calibrated v(z)", zc.data(), tf.data(), int(zc.size()), Pal::Accent, 1.3f);
                        }
                        recEnd(); ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
            }
        }
        ImGui::EndTabItem();
    }
    geoCalcTabs(S, forceTab >= 5 ? forceTab - 5 : -1);
    ImGui::EndTabBar();
}
