// EOR section: polymer flooding (ECLIPSE / OPM Flow POLYMER option) - viscosity, adsorption / rock, Todd-Longstaff
// mixing regions, shear thinning, injection schedule (WPOLYMER), initial concentration and run results.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include <array>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <string>

using namespace ressim;

namespace {
using Row2 = std::array<double, 2>;
using Row3 = std::array<double, 3>;
std::vector<Row2> toRows(const std::vector<double>& f) {
    std::vector<Row2> r;
    for (size_t i = 0; i + 1 < f.size(); i += 2) r.push_back({f[i], f[i + 1]});
    return r;
}
std::vector<double> fromRows(const std::vector<Row2>& r) {
    std::vector<double> f;
    for (auto& x : r) { f.push_back(x[0]); f.push_back(x[1]); }
    return f;
}
const char* concUnit(UnitSystem u) { return u == UnitSystem::Field ? "lb/stb" : "kg/sm3"; }
const char* velUnit(UnitSystem u) { return u == UnitSystem::Field ? "ft/day" : "m/day"; }

// typical data (opm-tests POLYMER-02), converted to FIELD when needed
void seedPolymer(Case& c) {
    PolymerInput& P = c.polymer;
    const bool field = c.units == UnitSystem::Field;
    const double cs = field ? 0.350507 : 1.0;          // kg/sm3 -> lb/stb
    const int npvt = std::max(1, c.props.numPvtRegions()), nsat = std::max(1, c.props.numSatRegions());
    if (P.plyvisc.empty()) {
        std::vector<double> t;
        const double v[][2] = {{0, 1}, {0.5, 3}, {1, 6}, {1.5, 12}, {2, 24}, {3, 48}};
        for (auto& r : v) { t.push_back(r[0] * cs); t.push_back(r[1]); }
        P.plyvisc.assign(size_t(npvt), t);
    }
    if (P.plyads.empty()) {
        std::vector<double> t;
        const double v[][2] = {{0, 0}, {0.25, 1.2e-5}, {0.5, 1.6e-5}, {1, 2e-5}, {1.5, 2.3e-5}, {1.75, 2.5e-5}, {3, 2.5e-5}};
        for (auto& r : v) { t.push_back(r[0] * cs); t.push_back(r[1]); }
        P.plyads.assign(size_t(nsat), t);
    }
    if (P.plyrock.empty()) {
        PolymerRock R;
        R.deadPoreVolume = 0.05; R.rrf = 1.2; R.rockDensity = field ? 2500 * 0.350507 : 2500; R.adsIndex = 2; R.maxAdsorption = 2.5e-5;
        P.plyrock.assign(size_t(nsat), R);
    }
    if (P.cmax.empty() || (P.cmax.size() == 1 && P.cmax[0] <= 0)) { P.cmax = {3 * cs}; P.cmaxSalt = {0}; P.mixing = {0.3}; }
}
template <class T> void fitRegions(std::vector<T>& v, int n) {
    if (v.empty() || n <= 0) return;
    if (int(v.size()) < n) v.resize(size_t(n), v.back());
}
void regionCombo(const char* id, const char* what, int* sel, int n) {
    if (n <= 1) { *sel = 0; return; }
    *sel = std::clamp(*sel, 0, n - 1);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
    char cur[48]; std::snprintf(cur, sizeof(cur), "%s %d", what, *sel + 1);
    if (ImGui::BeginCombo(id, cur)) {
        for (int r = 0; r < n; ++r) {
            char b[48]; std::snprintf(b, sizeof(b), "%s %d", what, r + 1);
            if (ImGui::Selectable(b, r == *sel)) *sel = r;
        }
        ImGui::EndCombo();
    }
}
} // namespace

void App::editEOR() {
    PolymerInput& P = case_.polymer;
    const UnitSystem us = case_.units;
    bool ch = false;
    ui::header("Polymer flooding (POLYMER)");
    ImGui::BeginDisabled(!case_.phases.water);
    bool on = case_.phases.polymer;
    if (ImGui::Checkbox("Polymer component in the water phase (POLYMER)", &on)) {
        case_.phases.polymer = on;
        if (on) seedPolymer(case_);
        ch = true;
    }
    ImGui::EndDisabled();
    if (!case_.phases.water) ui::hint("Polymer is carried by the water phase - activate WATER in Run specification first.");
    ui::hint("ECLIPSE / OPM Flow polymer model: Todd-Longstaff viscosity mixing (PLYVISC, PLYMAX, PLMIXPAR per PLMIXNUM region), "
             "adsorption on the rock (PLYADS, PLYROCK; reversible or irreversible), permeability reduction (residual resistance "
             "factor), inaccessible pore volume, shear thinning / thickening in faces and well connections (PLYSHLOG) and injector "
             "full mixing. Transport is solved sequentially-implicitly with outer coupling to the flow equations each time step; "
             "validated against OPM Flow (opm-tests POLYMER-02 / 03).");
    if (!case_.phases.polymer) {
        if (ch) markDirty();
        editChemEOR();
        return;
    }
    const int npvt = std::max(1, case_.props.numPvtRegions()), nsat = std::max(1, case_.props.numSatRegions());
    fitRegions(P.plyvisc, npvt); fitRegions(P.plyads, nsat); fitRegions(P.plyrock, nsat);
    char cu[64]; std::snprintf(cu, sizeof(cu), "c (%s)", concUnit(us));

    if (ImGui::BeginTabBar("eor")) {
        // ---------------------------------------------------------------- viscosity
        if (ui::tabItem("Viscosity (PLYVISC)")) {
            regionCombo("##pvtr", "PVT region", &eorPvtSel_, npvt);
            auto rows = toRows(P.plyvisc[size_t(eorPvtSel_)]);
            if (ui::editTable<Row2>("plyvisc", rows, {cu, "water viscosity multiplier"}, [](Row2& r, int c) -> double* { return &r[size_t(c)]; })) {
                P.plyvisc[size_t(eorPvtSel_)] = fromRows(rows); ch = true;
            }
            ui::hint("Multiplier of the water viscosity for a fully mixed polymer solution at concentration c. The first row must be "
                     "c = 0 with multiplier 1; concentrations increase down the table.");
            if (ImPlot::BeginPlot("##plyvisc", ImVec2(-1, ImGui::GetFontSize() * 13))) {
                ImPlot::SetupAxes(cu, "viscosity multiplier", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                std::vector<double> x, y; for (auto& r : rows) { x.push_back(r[0]); y.push_back(r[1]); }
                ImPlotSpec a; a.LineColor = col4(Pal::Water); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                if (!x.empty()) ImPlot::PlotLine("PLYVISC", x.data(), y.data(), int(x.size()), a);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        // ---------------------------------------------------------------- adsorption and rock
        if (ui::tabItem("Adsorption & rock (PLYADS, PLYROCK)")) {
            regionCombo("##satr", "Sat. region", &eorSatSel_, nsat);
            const size_t r = size_t(eorSatSel_);
            ui::header("Adsorption isotherm (PLYADS)");
            auto rows = toRows(P.plyads[r]);
            if (ui::editTable<Row2>("plyads", rows, {cu, "adsorbed (mass / rock mass)"}, [](Row2& x, int c) -> double* { return &x[size_t(c)]; },
                                    0, nullptr, "%.6g")) {
                P.plyads[r] = fromRows(rows); ch = true;
            }
            ui::header("Polymer-rock properties (PLYROCK)");
            PolymerRock& R = P.plyrock[r];
            if (ui::beginProps("plyrock")) {
                ch |= ui::propD("Inaccessible pore volume", &R.deadPoreVolume, "%.4g", "Fraction of the pore volume polymer cannot enter (dead pore volume).");
                ch |= ui::propD("Residual resistance factor", &R.rrf, "%.4g", "Permeability reduction to water at the maximum adsorption (>= 1).");
                char rd[64]; std::snprintf(rd, sizeof(rd), "Rock density (%s)", us == UnitSystem::Field ? "lb/rb" : "kg/rm3");
                ch |= ui::propD(rd, &R.rockDensity, "%.6g");
                static const char* adsT[] = {"1 - desorption allowed", "2 - irreversible (no desorption)"};
                int ai = R.adsIndex == 2 ? 1 : 0;
                if (ui::propCombo("Adsorption index", &ai, adsT, 2)) { R.adsIndex = ai ? 2 : 1; ch = true; }
                ch |= ui::propD("Maximum adsorption", &R.maxAdsorption, "%.6g", "Adsorbed polymer (mass / rock mass) at which the full RRF applies.");
                ui::endProps();
            }
            if (R.rrf < 1) { R.rrf = 1; ch = true; }
            if (ImPlot::BeginPlot("##plyads", ImVec2(-1, ImGui::GetFontSize() * 11))) {
                ImPlot::SetupAxes(cu, "adsorbed", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                std::vector<double> x, y; for (auto& q : rows) { x.push_back(q[0]); y.push_back(q[1]); }
                ImPlotSpec a; a.LineColor = col4(Pal::Warn); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                if (!x.empty()) ImPlot::PlotLine("PLYADS", x.data(), y.data(), int(x.size()), a);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        // ---------------------------------------------------------------- mixing regions
        if (ui::tabItem("Mixing (PLYMAX, PLMIXPAR, PLMIXNUM)")) {
            const int nm = std::max(1, P.numMixRegions());
            P.cmax.resize(size_t(nm), P.cmax.empty() ? 1.0 : P.cmax.back());
            P.mixing.resize(size_t(nm), P.mixing.empty() ? 1.0 : P.mixing.back());
            P.cmaxSalt.resize(size_t(nm), 0.0);
            std::vector<Row3> mr;
            for (int k = 0; k < nm; ++k) mr.push_back({P.cmax[size_t(k)], P.mixing[size_t(k)], P.cmaxSalt[size_t(k)]});
            char cm[64]; std::snprintf(cm, sizeof(cm), "c max (%s)", concUnit(us));
            if (ui::editTable<Row3>("plmix", mr, {cm, "mixing \xCF\x89", "salt (export)"}, [](Row3& x, int c) -> double* { return &x[size_t(c)]; })) {
                if (mr.empty()) mr.push_back({1.0, 1.0, 0.0});
                P.cmax.clear(); P.mixing.clear(); P.cmaxSalt.clear();
                for (auto& x : mr) { P.cmax.push_back(x[0]); P.mixing.push_back(std::clamp(x[1], 0.0, 1.0)); P.cmaxSalt.push_back(x[2]); }
                ch = true;
            }
            ui::hint("One row per PLMIXNUM region. Todd-Longstaff: \xCE\xBC_eff = \xCE\xBC_m^\xCF\x89 \xC2\xB7 \xCE\xBC_p^(1-\xCF\x89); "
                     "\xCF\x89 = 1 full mixing, \xCF\x89 = 0 fully segregated. c max is the concentration of the injected solution "
                     "used for the segregated viscosity.");
            ui::header("Region map (PLMIXNUM)");
            const int n = case_.grid.size();
            auto& mx = case_.grid.plmixnum;
            std::vector<int> cnt(size_t(nm) + 1, 0);
            if (mx.empty()) cnt[1] = n; else for (int v : mx) if (v >= 1 && v <= nm) ++cnt[size_t(v)]; else ++cnt[0];
            for (int k = 1; k <= nm; ++k) { ImGui::Text("Region %d: %d cells", k, cnt[size_t(k)]); if (k < nm) ImGui::SameLine(0, 18); }
            if (cnt[0]) ImGui::TextColored(col4(Pal::Err), "%d cells refer to a region without a row above", cnt[0]);
            static int box[7] = {1, 1, 1, 1, 1, 1, 2};
            const GridInput& g = case_.grid;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("I1, I2", &box[0]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("J1, J2", &box[2]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9); ImGui::InputInt2("K1, K2", &box[4]);
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputInt("region", &box[6], 0);
            box[6] = std::clamp(box[6], 1, nm);
            if (ImGui::Button("Assign box")) {
                if (mx.size() != size_t(n)) mx.assign(size_t(n), 1);
                const int i1 = std::clamp(std::min(box[0], box[1]), 1, g.nx), i2 = std::clamp(std::max(box[0], box[1]), 1, g.nx);
                const int j1 = std::clamp(std::min(box[2], box[3]), 1, g.ny), j2 = std::clamp(std::max(box[2], box[3]), 1, g.ny);
                const int k1 = std::clamp(std::min(box[4], box[5]), 1, g.nz), k2 = std::clamp(std::max(box[4], box[5]), 1, g.nz);
                for (int k = k1; k <= k2; ++k) for (int j = j1; j <= j2; ++j) for (int i = i1; i <= i2; ++i)
                    mx[size_t((i - 1) + g.nx * ((j - 1) + g.ny * (k - 1)))] = box[6];
                ch = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("All cells region 1")) { mx.clear(); ch = true; }
            ImGui::EndTabItem();
        }
        // ---------------------------------------------------------------- shear
        if (ui::tabItem("Shear (PLYSHLOG)")) {
            bool shear = !P.plyshlog.empty();
            if (ImGui::Checkbox("Shear thinning / thickening (PLYSHLOG)", &shear)) {
                if (shear) {
                    PolymerInput::Shear s; s.refConc = P.cmax.empty() ? 1.0 : 0.5 * P.cmax[0];
                    const double vs = us == UnitSystem::Field ? 3.28084 : 1.0;
                    const double t[][2] = {{1e-7, 1.0}, {1e-6, 1.1}, {1e-4, 1.3}, {1e-3, 1.47}, {1e-2, 1.67}, {0.1, 2.0}, {1, 2.2}, {10, 2.3}, {100, 2.4}, {1000, 2.4}};
                    for (auto& r : t) { s.rows.push_back(r[0] * vs); s.rows.push_back(r[1]); }
                    P.plyshlog.assign(size_t(npvt), s);
                } else P.plyshlog.clear();
                ch = true;
            }
            if (!P.plyshlog.empty()) {
                fitRegions(P.plyshlog, npvt);
                regionCombo("##shpvt", "PVT region", &eorPvtSel_, npvt);
                PolymerInput::Shear& S = P.plyshlog[size_t(eorPvtSel_)];
                if (ui::beginProps("shl")) {
                    char rc[64]; std::snprintf(rc, sizeof(rc), "Reference concentration (%s)", concUnit(us));
                    ch |= ui::propD(rc, &S.refConc, "%.5g", "Concentration at which the multiplier table was measured.");
                    ui::endProps();
                }
                char vc[64]; std::snprintf(vc, sizeof(vc), "water velocity (%s)", velUnit(us));
                auto rows = toRows(S.rows);
                if (ui::editTable<Row2>("plyshlog", rows, {vc, "viscosity multiplier"}, [](Row2& x, int c) -> double* { return &x[size_t(c)]; })) {
                    S.rows = fromRows(rows); ch = true;
                }
                ui::hint("Multiplier of the polymer solution viscosity at the reference concentration versus the Darcy water velocity "
                         "(interpolated in log velocity, ECLIPSE / OPM Flow). Applied in the reservoir faces and at the well "
                         "connections (well-bore area 2\xCF\x80\xE2\x88\x9A(r_e r_w) h).");
                if (ImPlot::BeginPlot("##shl", ImVec2(-1, ImGui::GetFontSize() * 12))) {
                    ImPlot::SetupAxes(vc, "multiplier", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                    std::vector<double> x, y; for (auto& q : rows) if (q[0] > 0) { x.push_back(q[0]); y.push_back(q[1]); }
                    ImPlotSpec a; a.LineColor = col4(Pal::Accent); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                    if (!x.empty()) ImPlot::PlotLine("PLYSHLOG", x.data(), y.data(), int(x.size()), a);
                    ImPlot::EndPlot();
                }
            }
            ImGui::EndTabItem();
        }
        // ---------------------------------------------------------------- injection schedule and initial concentration
        if (ui::tabItem("Injection (WPOLYMER, SPOLY)")) {
            auto& S = case_.schedule;
            std::vector<std::string> inj;
            for (auto& w : case_.wells) if (w.type == WellType::WaterInjector) inj.push_back(w.name);
            ui::header("Polymer slug");
            static int slugWell = 0, slugFrom = 0, slugTo = 0; static double slugC = 1.0;
            if (inj.empty()) ui::hint("No water injector in the case - add one in Wells.");
            else if (!S.empty()) {
                slugWell = std::clamp(slugWell, 0, int(inj.size()) - 1);
                slugFrom = std::clamp(slugFrom, 0, int(S.size()) - 1); slugTo = std::clamp(slugTo, slugFrom, int(S.size()) - 1);
                std::vector<double> tStart(S.size(), 0.0);
                for (size_t k = 1; k < S.size(); ++k) tStart[k] = tStart[k - 1] + S[k - 1].length;
                if (ui::beginProps("slug")) {
                    ui::label("Injector");
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::BeginCombo("##sw", inj[size_t(slugWell)].c_str())) {
                        for (int k = 0; k < int(inj.size()); ++k) if (ImGui::Selectable(inj[size_t(k)].c_str(), k == slugWell)) slugWell = k;
                        ImGui::EndCombo();
                    }
                    char b[64]; std::snprintf(b, sizeof(b), "First step (t = %.0f d)", tStart[size_t(slugFrom)]);
                    ui::propI(b, &slugFrom, 0, int(S.size()) - 1);
                    std::snprintf(b, sizeof(b), "Last step (ends t = %.0f d)", tStart[size_t(slugTo)] + S[size_t(slugTo)].length);
                    ui::propI(b, &slugTo, slugFrom, int(S.size()) - 1);
                    char cc[64]; std::snprintf(cc, sizeof(cc), "Concentration (%s)", concUnit(us));
                    ui::propD(cc, &slugC, "%.5g");
                    ui::endProps();
                }
                if (ui::accentButton("Apply slug")) {
                    const std::string& wn = inj[size_t(slugWell)];
                    auto setAt = [&](size_t k, double c) {
                        auto& v = S[k].wpolymer;
                        v.erase(std::remove_if(v.begin(), v.end(), [&](const WellPolymer& p) { return p.well == wn; }), v.end());
                        v.push_back({wn, c});
                    };
                    for (size_t k = size_t(slugFrom) + 1; k <= size_t(slugTo) && k < S.size(); ++k) {
                        auto& v = S[k].wpolymer;
                        v.erase(std::remove_if(v.begin(), v.end(), [&](const WellPolymer& p) { return p.well == wn; }), v.end());
                    }
                    setAt(size_t(slugFrom), slugC);
                    if (size_t(slugTo) + 1 < S.size()) setAt(size_t(slugTo) + 1, 0.0);
                    ch = true;
                }
                ImGui::SameLine();
                ui::hint("Sets WPOLYMER at the first step and resets it to 0 after the last step (chase water).");
            }
            ui::header("WPOLYMER records by report step");
            if (!S.empty()) {
                eorStepSel_ = std::clamp(eorStepSel_, 0, int(S.size()) - 1);
                double t0 = 0; for (int k = 0; k < eorStepSel_; ++k) t0 += S[size_t(k)].length;
                char b[96]; std::snprintf(b, sizeof(b), "Step %d (t = %.1f d)", eorStepSel_ + 1, t0);
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
                if (ImGui::SliderInt("##stp", &eorStepSel_, 0, int(S.size()) - 1, b)) {}
                auto& v = S[size_t(eorStepSel_)].wpolymer;
                int del = -1;
                if (ui::beginTable("wpol", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame)) {
                    ui::tableColumn("Well (name or pattern*)");
                    char cc[64]; std::snprintf(cc, sizeof(cc), "Concentration (%s)", concUnit(us));
                    ui::tableColumn(cc);
                    ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
                    ImGui::TableHeadersRow();
                    for (int r = 0; r < int(v.size()); ++r) {
                        ImGui::PushID(r);
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::SetNextItemWidth(-FLT_MIN);
                        if (ui::inputString("##w", &v[size_t(r)].well)) ch = true;
                        ImGui::TableSetColumnIndex(1); ImGui::SetNextItemWidth(-FLT_MIN);
                        if (ImGui::InputDouble("##c", &v[size_t(r)].conc, 0, 0, "%.5g")) { v[size_t(r)].conc = std::max(0.0, v[size_t(r)].conc); ch = true; }
                        ImGui::TableSetColumnIndex(2); if (ImGui::SmallButton("x")) del = r;
                        ImGui::PopID();
                    }
                    ui::endTable();
                }
                if (del >= 0) { v.erase(v.begin() + del); ch = true; }
                if (ImGui::SmallButton("+ record")) { v.push_back({inj.empty() ? std::string("INJ*") : inj.front(), 0.0}); ch = true; }
                int nrec = 0; for (auto& s : S) nrec += int(s.wpolymer.size());
                ImGui::SameLine(); ImGui::TextDisabled("%d WPOLYMER record(s) in the schedule; a concentration stays in force until changed", nrec);
            }
            ui::header("Initial polymer concentration (SPOLY)");
            static double spoly = 0;
            if (ui::beginProps("spoly")) {
                char cc[64]; std::snprintf(cc, sizeof(cc), "Uniform value (%s)", concUnit(us));
                ui::propD(cc, &spoly, "%.5g");
                ui::endProps();
            }
            if (ImGui::Button("Apply to all cells")) { case_.grid.spoly.assign(size_t(case_.grid.size()), std::max(0.0, spoly)); ch = true; }
            ImGui::SameLine();
            if (ImGui::Button("Clear (no initial polymer)")) { case_.grid.spoly.clear(); ch = true; }
            ImGui::SameLine();
            ImGui::TextDisabled(case_.grid.spoly.empty() ? "SPOLY not set (0)" : "SPOLY set per cell");
            ImGui::EndTabItem();
        }
        // ---------------------------------------------------------------- results of the last run
        if (ui::tabItem("Results")) {
            std::vector<double> t, cir, cpr, cit, cpt, cip, cad;
            bool have = false;
            {
                std::lock_guard<std::mutex> l(runner_.mutex());
                const RunResults& R = runner_.results();
                if (R.valid && R.summary.size() > 1) {
                    have = true;
                    for (auto& sp : R.summary) {
                        t.push_back(sp.time / 365.25); cir.push_back(sp.fcir); cpr.push_back(sp.fcpr);
                        cit.push_back(sp.fcit); cpt.push_back(sp.fcpt); cip.push_back(sp.fcip); cad.push_back(sp.fcad);
                    }
                }
            }
            if (!have) ui::hint("Run the case to see polymer injection, production, retention and the mass balance.");
            else {
                const char* mu = us == UnitSystem::Field ? "lb" : "kg";
                char ya[64]; std::snprintf(ya, sizeof(ya), "rate (%s/day)", mu);
                if (ImPlot::BeginPlot("Polymer rates", ImVec2(-1, ImGui::GetFontSize() * 12))) {
                    ImPlot::SetupAxes("years", ya, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec a; a.LineColor = col4(Pal::Water); a.LineWeight = 2;
                    ImPlotSpec b; b.LineColor = col4(Pal::Oil); b.LineWeight = 2;
                    ImPlot::PlotLine("FCIR injected", t.data(), cir.data(), int(t.size()), a);
                    ImPlot::PlotLine("FCPR produced", t.data(), cpr.data(), int(t.size()), b);
                    ImPlot::EndPlot();
                }
                char yb[64]; std::snprintf(yb, sizeof(yb), "mass (%s)", mu);
                if (ImPlot::BeginPlot("Polymer mass", ImVec2(-1, ImGui::GetFontSize() * 13))) {
                    ImPlot::SetupAxes("years", yb, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec a; a.LineColor = col4(Pal::Water); a.LineWeight = 2;
                    ImPlotSpec b; b.LineColor = col4(Pal::Oil); b.LineWeight = 2;
                    ImPlotSpec c; c.LineColor = col4(Pal::Accent); c.LineWeight = 2;
                    ImPlotSpec d; d.LineColor = col4(Pal::Warn); d.LineWeight = 2;
                    ImPlot::PlotLine("FCIT injected", t.data(), cit.data(), int(t.size()), a);
                    ImPlot::PlotLine("FCPT produced", t.data(), cpt.data(), int(t.size()), b);
                    ImPlot::PlotLine("FCIP in solution", t.data(), cip.data(), int(t.size()), c);
                    ImPlot::PlotLine("FCAD adsorbed", t.data(), cad.data(), int(t.size()), d);
                    ImPlot::EndPlot();
                }
                const double bal = cit.back() - cpt.back() - cip.back() - cad.back();
                ImGui::Text("Polymer mass balance at the end: injected %.6g = produced %.6g + in solution %.6g + adsorbed %.6g  (error %.2e)",
                            cit.back(), cpt.back(), cip.back(), cad.back(), std::fabs(bal) / std::max(1e-30, std::fabs(cit.back())));
                ui::hint("3-D maps: choose POLYMER (concentration) or PADS (adsorbed) in the viewer / Results property list.");
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    if (ch) markDirty();
    editChemEOR();
}
