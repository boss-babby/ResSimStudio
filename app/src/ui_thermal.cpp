// Thermal (E300 THERMAL / STARS-like) model editor: components and K-values, oil / gas viscosity tables, rock heat
// capacity and conductivity, heat loss, initial temperature and steam injection (WINJTEMP), with steam-table and
// K-value / viscosity previews. All values are in deck units of the case (METRIC / FIELD), as in the exported deck.
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "implot.h"
#include "ressim/steam.h"
#include "ressim/examples.h"
#include "ressim/geomech.h"
#include "ressim/gasdata.h"
#include <cctype>
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ressim;

namespace {
struct TU {   // thermal deck units
    bool f;
    const char* t() const { return f ? "\xC2\xB0" "F" : "\xC2\xB0" "C"; }
    const char* p() const { return f ? "psia" : "bar"; }
    const char* dens() const { return f ? "lb/ft3" : "kg/m3"; }
    const char* cp() const { return f ? "Btu/lb/\xC2\xB0R" : "kJ/kg/K"; }
    const char* hv() const { return f ? "Btu/lb" : "kJ/kg"; }
    const char* hcr() const { return f ? "Btu/ft3/\xC2\xB0R" : "kJ/m3/K"; }
    const char* thc() const { return f ? "Btu/ft/d/\xC2\xB0R" : "kJ/m/d/K"; }
    double toK(double t) const { return f ? (t + 459.67) / 1.8 : t + 273.15; }
    double fromK(double K) const { return f ? K * 1.8 - 459.67 : K - 273.15; }
    double toPa(double p) const { return f ? p * 6894.757293168 : p * 1e5; }
    double fromPa(double Pa) const { return f ? Pa / 6894.757293168 : Pa / 1e5; }
    double absT(double t) const { return f ? t + 459.67 : t + 273.15; }   // deck T -> absolute (K | R)
};

void reshapeThermal(ThermalInput& T, int at, bool ins) {
    auto fixRow = [&](std::vector<double>& r, double def) {   // row: T, mu_1..mu_nc
        const size_t k = size_t(at) + 1;
        if (ins) r.insert(r.begin() + long(std::min(k, r.size())), def);
        else if (k < r.size()) r.erase(r.begin() + long(k));
    };
    for (auto& r : T.oilViscT) fixRow(r, r.size() > 1 ? r[1] : 1.0);
    for (auto& r : T.gasViscT) fixRow(r, r.size() > 1 ? r[1] : 0.01);
    if (ins) T.zi.insert(T.zi.begin() + long(std::min<size_t>(size_t(at), T.zi.size())), 0.0);
    else if (at < int(T.zi.size())) T.zi.erase(T.zi.begin() + at);
    T.xmf.clear(); T.ymf.clear();   // per-cell compositions no longer match: ZI is used
    auto fixV = [&](std::vector<double>& v, double def) {   // per-component vectors (nc or nc + 1 long)
        if (v.empty()) return;
        if (ins) v.insert(v.begin() + long(std::min<size_t>(size_t(at), v.size())), def);
        else if (at < int(v.size())) v.erase(v.begin() + at);
    };
    for (auto& r : T.reactions) {
        fixV(r.sto, 0.0); fixV(r.prod, 0.0); fixV(r.order, 0.0);
        if (!r.phase.empty()) { if (ins) r.phase.insert(r.phase.begin() + long(std::min<size_t>(size_t(at), r.phase.size())), -1); else if (at < int(r.phase.size())) r.phase.erase(r.phase.begin() + at); }
        if (r.o2pp == at && !ins) r.o2pp = -1; else if (r.o2pp >= at) r.o2pp += ins ? 1 : -1;
    }
    fixV(T.solidConc, 0.0);
    for (auto& st : T.streams) fixV(st.z, 0.0);
}
const char* kCompTypes[] = {"Fluid (K-value)", "Gas (non-condensable)", "Solid"};
}  // namespace

void App::editThermal() {
    ThermalInput& T = case_.thermal;
    const TU U{case_.units == UnitSystem::Field};
    auto changed = [&]() { T.structured = true; markDirty(); };
    if (T.comps.empty()) { thermalPreset(case_, 25); changed(); }
    const int nc = T.nc();
    if (int(T.zi.size()) != nc) { T.zi.assign(size_t(nc), 0.0); if (nc) T.zi[0] = 1; }
    for (auto& r : T.oilViscT) if (int(r.size()) != nc + 1) r.resize(size_t(nc) + 1, 1.0);
    for (auto& r : T.gasViscT) if (int(r.size()) != nc + 1) r.resize(size_t(nc) + 1, 0.01);
    char lab[128];
    if (!ImGui::BeginTabBar("therm")) return;

    // ------------------------------------------------------------------ components
    if (ui::tabItem("Components")) {
        ui::hint("STARS-like thermal core, fully implicit in pressure, component moles and internal energy (p, N, U). Water / steam: "
                 "IAPWS-IF97 + IAPWS-2008 viscosity. Hydrocarbons: dead (non-volatile) or volatile with K = (A + B/p + C p) exp(-D/(T - E)) "
                 "(E300 KVCR, Crookston form). Free-water Rachford-Rice flash; steam is the water component in the gas phase.");
        if (ui::beginProps("thg")) {
            if (ui::propBool("Dead oil (DEADOIL)", &T.deadOil, "RUNSPEC DEADOIL: all hydrocarbons non-volatile")) changed();
            double st = T.stdT < -1e8 ? (U.f ? 60 : 15.56) : T.stdT, sp = T.stdP <= 0 ? (U.f ? 14.696 : 1.01325) : T.stdP;
            std::snprintf(lab, sizeof lab, "Standard temperature (%s)", U.t());
            if (ui::propD(lab, &st, "%.2f", "STCOND")) { T.stdT = st; changed(); }
            std::snprintf(lab, sizeof lab, "Standard pressure (%s)", U.p());
            if (ui::propD(lab, &sp, "%.4f", "STCOND") && sp > 0) { T.stdP = sp; changed(); }
            ui::endProps();
        }
        char hD[40], hCp[48], hHv[48], hB[40], hDk[40];
        std::snprintf(hD, sizeof hD, "DREF (%s)", U.dens()); std::snprintf(hCp, sizeof hCp, "SPECHA (%s)", U.cp());
        std::snprintf(hHv, sizeof hHv, "HEATVAP (%s)", U.hv()); std::snprintf(hB, sizeof hB, "KV B (%s)", U.p());
        std::snprintf(hDk, sizeof hDk, "KV D (%s)", U.f ? "\xC2\xB0R" : "K");
        {   // C16: Excel-style component table (copy / paste whole blocks of KVCR / thermal data from spreadsheets)
            static int compSel = 0;
            ThermalComponent* base = T.comps.data();
            ui::RecordTable<ThermalComponent> t;
            t.text("Name", [](ThermalComponent& q) -> std::string& { return q.name; })
             .choice("Type", [](ThermalComponent& q) -> int& { return q.type; }, {kCompTypes[0], kCompTypes[1], kCompTypes[2]})
             .num("ZI", [&T, base](ThermalComponent& q) -> double& { return T.zi[size_t(&q - base)]; })
             .num("MW", [](ThermalComponent& q) -> double& { return q.mw; })
             .flag("Dead", [](ThermalComponent& q) -> bool& { return q.dead; })
             .num("KV A", [](ThermalComponent& q) -> double& { return q.kvcr[0]; })
             .num(hB, [](ThermalComponent& q) -> double& { return q.kvcr[1]; })
             .num("KV C", [](ThermalComponent& q) -> double& { return q.kvcr[2]; })
             .num(hDk, [](ThermalComponent& q) -> double& { return q.kvcr[3]; })
             .num("KV E", [](ThermalComponent& q) -> double& { return q.kvcr[4]; })
             .num(hD, [](ThermalComponent& q) -> double& { return q.dref; })
             .num("CREF", [](ThermalComponent& q) -> double& { return q.cref; })
             .num("THERMEX1", [](ThermalComponent& q) -> double& { return q.thermex1; })
             .num(hCp, [](ThermalComponent& q) -> double& { return q.specha; })
             .num(hHv, [](ThermalComponent& q) -> double& { return q.heatvap; });
            t.fixedRows();   // add / delete with the buttons below (all per-component tables are reshaped together)
            t.onChange = [&T, nc](ThermalComponent& q) {
                q.type = std::clamp(q.type, 0, 2);
                if (q.type == 2) { q.dead = true; if (int(T.solidConc.size()) != nc) T.solidConc.resize(size_t(nc), 0.0); }
                if (q.type == 1) q.dead = false;
                if (q.dead) q.tcrit = std::max(q.tcrit, 1e5);
                if (q.name.empty()) q.name = "C";
            };
            const float h = std::min(ImGui::GetFrameHeightWithSpacing() * (nc + 1.8f), ImGui::GetFontSize() * 18);
            if (t.draw("tcomps", T.comps, h, &compSel, "%.6g", false)) changed();
            compSel = std::clamp(compSel, 0, std::max(0, nc - 1));
            ImGui::TextDisabled("Type: Fluid = K-value (KVCR) or dead oil | Gas = always in the gas phase (O2, N2, CO2, CO, H2, CH4) | Solid = immobile "
                                "(coke, coal, char). K = (A + B/p + C p) exp(-D / (T - E)).");
            if (nc > 1 && ImGui::Button(ICON_TRASH_2 " Delete selected component")) {
                const int del = compSel;
                T.comps.erase(T.comps.begin() + del); reshapeThermal(T, del, false); changed();
                compSel = std::max(0, del - 1);
                ImGui::EndTabItem(); ImGui::EndTabBar(); return;   // nc changed: redraw next frame
            }
            ImGui::SameLine();
        }
        if (nc < 13 && ImGui::Button(ICON_PLUS " Add component")) {
            ThermalComponent q; q.name = "C" + std::to_string(nc + 1); q.dead = false; q.tcrit = 600; q.pcrit = 20; q.mw = 150;
            q.kvcr[1] = U.f ? 83740 * 14.5037738 : 83740; q.kvcr[3] = U.f ? 5537 * 1.8 : 5537; q.dref = U.f ? 51.2 : 820;
            q.specha = U.f ? 0.5 : 2.1; q.heatvap = U.f ? 125 : 290;
            T.comps.push_back(q); reshapeThermal(T, nc, true); changed();
        }
        ImGui::SameLine();
        double s = 0; for (double z : T.zi) s += z;
        if (std::fabs(s - 1) > 1e-6) {
            ImGui::TextColored(col4(Pal::Warn), "ZI sums to %.4f", s); ImGui::SameLine();
            if (ImGui::SmallButton("Normalise")) { for (double& z : T.zi) z = s > 0 ? z / s : 1.0 / nc; changed(); }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset to heavy-oil preset")) { thermalPreset(case_, 25); changed(); }
        // K-value preview at a chosen pressure
        static double pK = 0;
        if (pK <= 0) pK = U.f ? 500 : 35;
        ui::header("K-values (KVCR) vs temperature");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        std::snprintf(lab, sizeof lab, "Pressure (%s)##pk", U.p());
        ImGui::InputDouble(lab, &pK, 0, 0, "%.2f");
        if (ImPlot::BeginPlot("##kv", ImVec2(-1, ImGui::GetFontSize() * 16))) {
            std::snprintf(lab, sizeof lab, "Temperature (%s)", U.t());
            ImPlot::SetupAxes(lab, "K", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < nc; ++i) {
                const ThermalComponent& q = T.comps[size_t(i)];
                if (q.dead) continue;
                std::vector<double> x, y;
                for (int k = 0; k <= 60; ++k) {
                    const double t = (U.f ? 50 : 10) + k * (U.f ? 10 : 5.5);
                    const double Ta = U.absT(t), den = Ta - q.kvcr[4];
                    if (den <= 0) continue;
                    const double K = (q.kvcr[0] + q.kvcr[1] / pK + q.kvcr[2] * pK) * std::exp(-q.kvcr[3] / den);
                    if (K > 0 && std::isfinite(K)) { x.push_back(t); y.push_back(K); }
                }
                ImPlot::PlotLine(q.name.c_str(), x.data(), y.data(), int(x.size()));
            }
            {   // water: K = psat(T) / p
                std::vector<double> x, y;
                for (int k = 0; k <= 60; ++k) { const double t = (U.f ? 50 : 10) + k * (U.f ? 10 : 5.5); x.push_back(t); y.push_back(steam::psat(U.toK(t)) / U.toPa(pK)); }
                ImPlot::PlotLine("H2O (psat/p)", x.data(), y.data(), int(x.size()));
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ C10: solubility, PR gas, trapping (unified with the CO2 cores)
    if (ui::tabItem("Solubility & storage")) {
        ui::hint("Unified multi-component phase behaviour: gases dissolve in the water phase by Henry's law (IAPWS G7-04 constants on the "
                 "saturation curve, Krichevsky-Kasarnovsky pressure correction with the partial molar volume, Sechenov salting-out), the gas "
                 "phase can use Peng-Robinson fugacities / volumes / residual enthalpy / LBC viscosity, and a three-phase (aqueous / oil / gas) "
                 "convex Michelsen-type flash replaces the free-water flash when any component is soluble. The tracked component (CO2 by "
                 "default) gets the same storage inventory as the compositional / CO2STORE cores: structural, residual (Land / Carlson gas "
                 "hysteresis), dissolved (water and oil) and mineral (solids with a FORMULA trackEq). Exported as '-- RESSIM:' lines.");
        if (ui::beginProps("tsol")) {
            if (ui::propBool("Peng-Robinson gas phase (GASEOS PR)", &T.gasPR, "components with critical constants (library by name or CRIT) form a PR mixture; others keep the Z-factor")) changed();
            if (ui::propD("Salinity (mol/kg NaCl)", &T.salinity, "%.3f", "salting-out of dissolved gases: ln(x0 / x) = ksalt m")) { T.salinity = std::max(0.0, T.salinity); changed(); }
            if (ui::propBool("Gas hysteresis (residual trapping)", &T.hystGas, "Land trapping + Carlson free-gas krg on imbibition (as the compositional core)")) changed();
            if (T.hystGas && ui::propD("Maximum trapped gas saturation Sgt,max", &T.sgtMax, "%.3f")) { T.sgtMax = std::clamp(T.sgtMax, 0.01, 0.9); changed(); }
            std::vector<std::string> opts = {"(auto: CO2 if present)"};
            for (int i = 0; i < nc; ++i) opts.push_back(T.comps[size_t(i)].name);
            int sel = 0;
            for (int i = 0; i < nc; ++i) if (!T.track.empty() && T.comps[size_t(i)].name == T.track) sel = i + 1;
            if (ui::propCombo("Tracked component (storage inventory)", &sel, opts, "inventory classes as the CO2 cores; CO2 by default")) {
                T.track = sel == 0 ? std::string() : opts[size_t(sel)]; changed();
            }
            if (ui::propD("Structural tolerance (fraction / year)", &T.structTol, "%.3g", "mobile gas counts as structurally held when its outflow < tol x volume per year")) changed();
            ui::endProps();
        }
        ui::header("Per-component solubility, critical constants and formula");
        {   // C16: Excel-style solubility / critical-constant table
            static int solSel = 0;
            ui::RecordTable<ThermalComponent> t;
            t.text("Name", [](ThermalComponent& q) -> std::string& { return q.name; })
             .choice("Aqueous model", [](ThermalComponent& q) -> int& { return q.aqModel; }, {"Insoluble", "IAPWS G7-04 (A, B, C)", "van't Hoff (kH298, dT)"})
             .num("A | kH298 (Pa)", [](ThermalComponent& q) -> double& { return q.aqModel == 2 ? q.kh298 : q.hA; })
             .num("B | dT (K)", [](ThermalComponent& q) -> double& { return q.aqModel == 2 ? q.hdT : q.hB; })
             .num("C", [](ThermalComponent& q) -> double& { return q.hC; })
             .num("vInf (cm3/mol, <0 Garcia)", [](ThermalComponent& q) -> double& { return q.vInf; }, 1e6)
             .num("ksalt (kg/mol)", [](ThermalComponent& q) -> double& { return q.ksalt; })
             .num("Tc (K, 0 = lib)", [](ThermalComponent& q) -> double& { return q.tcG; })
             .num("Pc (bar, 0 = lib)", [](ThermalComponent& q) -> double& { return q.pcG; }, 1e-5)
             .num("omega", [](ThermalComponent& q) -> double& { return q.omegaG; })
             .text("Formula", [](ThermalComponent& q) -> std::string& { return q.formula; })
             .num("trackEq (solids)", [](ThermalComponent& q) -> double& { return q.trackEq; });
            t.fixedRows();
            t.onChange = [](ThermalComponent& q) { q.aqModel = std::clamp(q.aqModel, 0, 2); if (q.type == 2) q.aqModel = 0; };
            const float h = std::min(ImGui::GetFrameHeightWithSpacing() * (nc + 1.8f), ImGui::GetFontSize() * 18);
            if (t.draw("tsolt", T.comps, h, &solSel, "%.6g", false)) changed();
            solSel = std::clamp(solSel, 0, std::max(0, nc - 1));
            auto fillLib = [&](ThermalComponent& q) {
                const gasdata::GasInfo* gi = gasdata::find(q.name);
                if (!gi) return false;
                if (gi->hA != 0 && q.type != 2) { q.aqModel = 1; q.hA = gi->hA; q.hB = gi->hB; q.hC = gi->hC; q.vInf = gi->vInf; q.ksalt = gi->ksalt; }
                q.tcG = gi->tc; q.pcG = gi->pc; q.omegaG = gi->omega; q.formula = gi->formula;
                if (q.type == 2) q.trackEq = gasdata::carbonateCO2(q.formula);
                return true;
            };
            if (nc > 0) {
                ThermalComponent& q = T.comps[size_t(solSel)];
                const bool have = gasdata::find(q.name) != nullptr;
                ImGui::BeginDisabled(!have);
                char bl[96]; std::snprintf(bl, sizeof bl, "Library data for %s", q.name.c_str());
                if (ImGui::Button(bl) && fillLib(q)) changed();
                ImGui::EndDisabled();
                ImGui::SameLine();
                if (ImGui::Button("Library data for all known components")) { bool any = false; for (auto& c : T.comps) any |= fillLib(c); if (any) changed(); }
            }
        }
        {
            std::string lib;
            for (const auto& g : gasdata::all()) { if (!lib.empty()) lib += ", "; lib += g.name; }
            ImGui::TextDisabled("Library: %s (aliases C1, C2, C3, METHANE, ...)", lib.c_str());
        }
        // solubility preview: x_i in pure water with the pure gas at pressure p (ideal gas fugacity), vs temperature
        static double pS = 0;
        if (pS <= 0) pS = U.f ? 1500 : 100;
        ui::header("Solubility preview (pure gas over water, ideal-gas fugacity)");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        std::snprintf(lab, sizeof lab, "Pressure (%s)##ps", U.p());
        ImGui::InputDouble(lab, &pS, 0, 0, "%.2f");
        if (ImPlot::BeginPlot("##sol", ImVec2(-1, ImGui::GetFontSize() * 15))) {
            std::snprintf(lab, sizeof lab, "Temperature (%s)", U.t());
            ImPlot::SetupAxes(lab, "mole fraction in water", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            const double p = U.toPa(pS);
            for (int i = 0; i < nc; ++i) {
                const ThermalComponent& q = T.comps[size_t(i)];
                if (q.aqModel == 0 || q.type == 2) continue;
                std::vector<double> x, y;
                for (int k = 0; k <= 60; ++k) {
                    const double TK = 278.15 + k * 5.0;
                    if (TK > 620) break;
                    const double kH = q.aqModel == 1 ? gasdata::henryIapws(q.hA, q.hB, q.hC, TK) : q.kh298 * std::exp(-q.hdT * (1.0 / TK - 1.0 / 298.15));
                    const double ps = steam::psat(TK);
                    if (!(kH > 0) || p <= ps) continue;
                    const double vi = q.vInf < 0 ? gasdata::vInfCO2(TK) : q.vInf;
                    const double H = kH * std::exp(vi * (p - ps) / (8.314462618 * TK) + q.ksalt * T.salinity);
                    const double xi = (p - ps) / H;
                    if (xi > 0 && std::isfinite(xi)) { x.push_back(U.fromK(TK)); y.push_back(xi); }
                }
                ImPlot::PlotLine(q.name.c_str(), x.data(), y.data(), int(x.size()));
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ viscosity
    if (ui::tabItem("Viscosity")) {
        ui::hint("OILVISCT: liquid viscosity of each component vs temperature (cP); the oil-phase viscosity is the log-linear mole-fraction "
                 "mix (E300 default). GASVISCT likewise for the gas phase. Water: WATVISCT or IAPWS-2008 (written as WATVISCT on export).");
        for (int tab = 0; tab < 2; ++tab) {
            auto& V = tab == 0 ? T.oilViscT : T.gasViscT;
            ui::header(tab == 0 ? "OILVISCT (cP)" : "GASVISCT (cP)");
            ImGui::PushID(tab);
            int delRow = -1;
            if (ui::beginTable("vt", nc + 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame)) {
                std::snprintf(lab, sizeof lab, "T (%s)", U.t());
                ui::tableColumn(lab);
                for (int i = 0; i < nc; ++i) ui::tableColumn(T.comps[size_t(i)].name.c_str());
                ui::tableColumn("");
                ImGui::TableHeadersRow();
                for (int r = 0; r < int(V.size()); ++r) {
                    ImGui::PushID(r);
                    ImGui::TableNextRow();
                    for (int k = 0; k <= nc; ++k) {
                        ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                        char id[16]; std::snprintf(id, sizeof id, "##%d", k);
                        if (ImGui::InputDouble(id, &V[size_t(r)][size_t(k)], 0, 0, k == 0 ? "%.2f" : "%.5g")) changed();
                    }
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton(ICON_TRASH_2)) delRow = r;
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (delRow >= 0) { V.erase(V.begin() + delRow); changed(); }
            if (ImGui::SmallButton(ICON_PLUS " Row")) {
                std::vector<double> r = V.empty() ? std::vector<double>(size_t(nc) + 1, tab ? 0.01 : 1.0) : V.back();
                r[0] += U.f ? 50 : 25; V.push_back(r); changed();
            }
            ImGui::SameLine();
            if (ImGui::SmallButton("Sort by T")) { std::sort(V.begin(), V.end(), [](const auto& a, const auto& b) { return a[0] < b[0]; }); changed(); }
            for (size_t r = 1; r < V.size(); ++r)
                if (V[r][0] <= V[r - 1][0]) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "temperatures must increase"); break; }
            ImGui::PopID();
        }
        if (ImPlot::BeginPlot("##visc", ImVec2(-1, ImGui::GetFontSize() * 18))) {
            std::snprintf(lab, sizeof lab, "Temperature (%s)", U.t());
            ImPlot::SetupAxes(lab, "Viscosity (cP)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < nc; ++i) {
                std::vector<double> x, y;
                for (const auto& r : T.oilViscT) if (r[size_t(i) + 1] > 0) { x.push_back(r[0]); y.push_back(r[size_t(i) + 1]); }
                ImPlot::PlotLine(T.comps[size_t(i)].name.c_str(), x.data(), y.data(), int(x.size()));
            }
            {   // live oil (ZI mix)
                std::vector<double> x, y;
                for (const auto& r : T.oilViscT) {
                    double lm = 0, zs = 0; bool ok = true;
                    for (int i = 0; i < nc; ++i) { if (r[size_t(i) + 1] <= 0) { ok = false; break; } lm += T.zi[size_t(i)] * std::log(r[size_t(i) + 1]); zs += T.zi[size_t(i)]; }
                    if (ok && zs > 0) { x.push_back(r[0]); y.push_back(std::exp(lm / zs)); }
                }
                ImPlot::PlotLine("Oil (ZI mix)", x.data(), y.data(), int(x.size()));
            }
            {   // water (IAPWS-2008 on the saturated liquid)
                std::vector<double> x, y;
                for (int k = 0; k <= 40; ++k) {
                    const double t = (U.f ? 50 : 10) + k * (U.f ? 13 : 7.5), K = U.toK(t);
                    const steam::Props L = steam::liquid(K, std::max(steam::psat(K), 1e5));
                    x.push_back(t); y.push_back(1e3 * steam::viscosity(K, 1.0 / L.v));
                }
                ImPlot::PlotLine("Water (IAPWS)", x.data(), y.data(), int(x.size()));
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ reactions, solids and permeability alteration
    if (ui::tabItem("Reactions")) {
        ui::hint("Chemical reactions (STARS-like): rate r = A exp(-Ea / RT) prod C_k^order_k (mol/m3 bulk/s, concentrations in mol per m3 of "
                 "bulk volume); optional O2 partial-pressure factor (kPa) instead of the O2 concentration; reversible reactions use "
                 "K_eq = exp(A_eq + B_eq / T). dH > 0 is exothermic (J/mol of reaction). Solids (coke, coal, char) occupy pore space; "
                 "their consumption opens porosity and permeability (Kozeny-Carman or exponential). Exported as '-- RESSIM:' lines "
                 "(ECLIPSE E300 has no reaction kinetics).");
        const int nk = nc + 1;
        auto cname = [&](int k) { return k < nc ? T.comps[size_t(k)].name : std::string("H2O"); };
        int delR = -1;
        for (int r = 0; r < int(T.reactions.size()); ++r) {
            ThermalReaction& R = T.reactions[size_t(r)];
            if (int(R.sto.size()) != nk) R.sto.resize(size_t(nk), 0.0);
            if (int(R.prod.size()) != nk) R.prod.resize(size_t(nk), 0.0);
            if (!R.order.empty() && int(R.order.size()) != nk) R.order.resize(size_t(nk), 0.0);
            ImGui::PushID(r);
            std::string eq;
            for (int k = 0; k < nk; ++k) if (R.sto[size_t(k)] > 0) { if (!eq.empty()) eq += " + "; char b2[48]; std::snprintf(b2, sizeof b2, "%g %s", R.sto[size_t(k)], cname(k).c_str()); eq += b2; }
            eq += R.reversible ? "  <=>  " : "  ->  ";
            bool firstP = true;
            for (int k = 0; k < nk; ++k) if (R.prod[size_t(k)] > 0) { if (!firstP) eq += " + "; firstP = false; char b2[48]; std::snprintf(b2, sizeof b2, "%g %s", R.prod[size_t(k)], cname(k).c_str()); eq += b2; }
            std::snprintf(lab, sizeof lab, "%s   (%s)###rx", R.name.c_str(), eq.c_str());
            const bool open = ImGui::CollapsingHeader(lab);
            if (open) {
                ImGui::Indent();
                {   // element balance (formulas: FORMULA / library by name; water H2O)
                    std::vector<std::string> fo;
                    for (int k = 0; k < nk; ++k) {
                        if (k == nc) { fo.push_back("H2O"); continue; }
                        const ThermalComponent& q = T.comps[size_t(k)];
                        const gasdata::GasInfo* gi = gasdata::find(q.name);
                        fo.push_back(!q.formula.empty() ? q.formula : gi ? std::string(gi->formula) : std::string());
                    }
                    std::map<std::string, double> per; std::vector<int> unk;
                    const double imb = gasdata::elementImbalance(fo, R.sto, R.prod, &per, &unk);
                    std::string u;
                    for (int k : unk) if (R.sto[size_t(k)] != 0 || R.prod[size_t(k)] != 0) { if (!u.empty()) u += ", "; u += cname(k); }
                    if (imb < 1e-6 && u.empty()) ImGui::TextColored(col4(Pal::Ok), ICON_CHECK " Element balance closed");
                    else {
                        std::string e;
                        for (const auto& kv : per) if (std::fabs(kv.second) > 1e-6) { char b3[48]; std::snprintf(b3, sizeof b3, "%s %+.4g  ", kv.first.c_str(), kv.second); e += b3; }
                        if (!e.empty()) ImGui::TextColored(col4(Pal::Warn), "Element imbalance (products - reactants, mol per mol of reaction): %s", e.c_str());
                        if (!u.empty()) ImGui::TextDisabled("No formula for: %s (set it in Solubility & storage)", u.c_str());
                    }
                }
                if (ui::beginProps("rxp")) {
                    std::string nm = R.name;
                    if (ui::propString("Name", &nm) && !nm.empty()) { R.name = nm; changed(); }
                    if (ui::propD("Frequency factor A", &R.A, "%.6g", "units follow the reaction order: (mol/m3)^(1 - sum order) / s [x 1/kPa for the O2 factor]")) changed();
                    double ea = R.Ea * 1e-3, dh = R.dH * 1e-3;
                    if (ui::propD("Activation energy Ea (kJ/mol)", &ea, "%.4g")) { R.Ea = ea * 1e3; changed(); }
                    if (ui::propD("Reaction enthalpy dH (kJ/mol, + exothermic)", &dh, "%.5g")) { R.dH = dh * 1e3; changed(); }
                    double tm = R.tMin > 0 ? U.fromK(R.tMin) : 0;
                    std::snprintf(lab, sizeof lab, "Minimum temperature (%s, 0 = none)", U.t());
                    if (ui::propD(lab, &tm, "%.1f", "the reaction is switched off below this temperature")) { R.tMin = tm != 0 ? U.toK(tm) : 0; changed(); }
                    std::vector<std::string> o2n = {"none (concentration)"};
                    for (int k = 0; k < nc; ++k) o2n.push_back(cname(k));
                    int o2 = R.o2pp + 1;
                    if (ui::propCombo("Partial-pressure factor", &o2, o2n, "STARS O2PP: this component enters as its gas partial pressure (kPa)")) { R.o2pp = o2 - 1; changed(); }
                    if (ui::propBool("Reversible (equilibrium)", &R.reversible, "net rate r (1 - Q / K_eq), ln K_eq = A_eq + B_eq / T")) changed();
                    if (R.reversible) {
                        if (ui::propD("A_eq (-)", &R.keqA, "%.5g")) changed();
                        if (ui::propD("B_eq (K)", &R.keqB, "%.5g")) changed();
                    }
                    ui::endProps();
                }
                // stoichiometry grid: rows reactant / product / order, columns components + water
                bool customOrder = !R.order.empty();
                if (ImGui::Checkbox("Custom reaction orders (default: reactant coefficients)", &customOrder)) { if (customOrder) R.order = R.sto; else R.order.clear(); changed(); }
                if (ui::beginTable("sto", nk + 1, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollX | ImGuiTableFlags_SizingFixedFit)) {
                    ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 6);
                    for (int k = 0; k < nk; ++k) ui::tableColumn(cname(k).c_str(), ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 5);
                    ImGui::TableHeadersRow();
                    for (int row = 0; row < (customOrder ? 3 : 2); ++row) {
                        std::vector<double>& v = row == 0 ? R.sto : row == 1 ? R.prod : R.order;
                        ImGui::TableNextRow(); ImGui::TableNextColumn();
                        ImGui::TextUnformatted(row == 0 ? "Reactants" : row == 1 ? "Products" : "Order");
                        for (int k = 0; k < nk; ++k) {
                            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                            char id[16]; std::snprintf(id, sizeof id, "##s%d_%d", row, k);
                            if (ImGui::InputDouble(id, &v[size_t(k)], 0, 0, "%.4g")) { v[size_t(k)] = std::max(0.0, v[size_t(k)]); changed(); }
                        }
                    }
                    ui::endTable();
                }
                {   // mass balance check (kg per mol of reaction)
                    double mr = 0, mp = 0;
                    for (int k = 0; k < nk; ++k) {
                        const double mw = k < nc ? T.comps[size_t(k)].mw : 18.01528;
                        mr += R.sto[size_t(k)] * mw; mp += R.prod[size_t(k)] * mw;
                    }
                    if (std::fabs(mr - mp) > 1e-3 * std::max(mr, 1.0)) ImGui::TextColored(col4(Pal::Warn), "Mass not conserved: reactants %.4g, products %.4g g/mol", mr, mp);
                    else ImGui::TextDisabled("Mass balanced: %.4g g per mol of reaction", mr);
                }
                if (ImGui::SmallButton(ICON_TRASH_2 " Delete reaction")) delR = r;
                ImGui::Unindent();
            }
            ImGui::PopID();
        }
        if (delR >= 0) { T.reactions.erase(T.reactions.begin() + delR); changed(); }
        if (ImGui::Button(ICON_PLUS " Add reaction")) {
            ThermalReaction R; R.name = "R" + std::to_string(T.reactions.size() + 1);
            R.sto.assign(size_t(nk), 0.0); R.prod.assign(size_t(nk), 0.0); R.A = 1; R.Ea = 8e4;
            T.reactions.push_back(R); changed();
        }
        ImGui::SameLine();
        if (ImGui::Button("Load in-situ combustion example")) { ImGui::EndTabItem(); ImGui::EndTabBar(); newCase(3, "combustion"); return; }
        ImGui::SameLine();
        if (ImGui::Button("Load UCG example")) { ImGui::EndTabItem(); ImGui::EndTabBar(); newCase(3, "ucg"); return; }
        // Arrhenius rate constants
        if (!T.reactions.empty() && ImPlot::BeginPlot("##arr", ImVec2(-1, ImGui::GetFontSize() * 14))) {
            std::snprintf(lab, sizeof lab, "Temperature (%s)", U.t());
            ImPlot::SetupAxes(lab, "k = A exp(-Ea/RT)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (const auto& R : T.reactions) {
                std::vector<double> x, y;
                for (int k = 0; k <= 80; ++k) {
                    const double K = 300 + k * 15.0;
                    const double kk = R.A * std::exp(-R.Ea / (8.314462618 * K));
                    if (kk > 1e-30 && std::isfinite(kk)) { x.push_back(U.fromK(K)); y.push_back(kk); }
                }
                ImPlot::PlotLine(R.name.c_str(), x.data(), y.data(), int(x.size()));
            }
            ImPlot::EndPlot();
        }
        // solids and permeability alteration
        ui::header("Solids and permeability alteration");
        bool anySolid = false;
        for (const auto& q : T.comps) anySolid = anySolid || q.type == 2;
        if (anySolid) {
            if (int(T.solidConc.size()) != nc) T.solidConc.resize(size_t(nc), 0.0);
            if (ui::beginProps("sol")) {
                for (int k = 0; k < nc; ++k) {
                    if (T.comps[size_t(k)].type != 2) continue;
                    std::snprintf(lab, sizeof lab, "Initial %s (mol / m3 pore volume)", cname(k).c_str());
                    if (ui::propD(lab, &T.solidConc[size_t(k)], "%.6g", "e.g. coal: 0.95 x 1350 kg/m3 / 12.5 g/mol per m3 of reactive pore space")) changed();
                    std::snprintf(lab, sizeof lab, "%s solid density (%s)", cname(k).c_str(), U.dens());
                    if (ui::propD(lab, &T.comps[size_t(k)].dref, "%.2f", "molar volume = MW / density (pore space taken by the solid)")) changed();
                }
                ui::endProps();
            }
        } else ImGui::TextDisabled("No solid components (set a component's Type to Solid in Components).");
        if (ui::beginProps("perm")) {
            const std::vector<std::string> laws = {"Off", "Kozeny-Carman (phi/phi0)^3 ((1-phi0)/(1-phi))^2", "Exponential exp(C (phi - phi0))"};
            if (ui::propCombo("Permeability law", &T.permLaw, laws, "fluid porosity changes from solid deposition / consumption and geomechanics")) changed();
            if (T.permLaw == 2 && ui::propD("Exponent C (-)", &T.permC, "%.4g")) changed();
            if (T.permLaw != 0 && ui::propD("Multiplier limit (-)", &T.permMultMax, "%.4g", "clamped to [1 / limit, limit]")) changed();
            ui::endProps();
        }
        // C12: Langmuir sorption (gas adsorption on coal / char / shale, solute retention: nanoparticles, polymer)
        ui::header("Sorption (Langmuir, competitive on shared sites)");
        ui::hint("C + Site <=> C_ads with r = A [c_C Site - C_ads / b(T)], b(T) = b_ref exp(dH/R (1/T - 1/T_ref)): at equilibrium the "
                 "extended Langmuir isotherm q_i = capacity b_i c_i / (1 + sum b_j c_j). c = gas partial pressure (kPa) or aqueous moles "
                 "per m3 bulk. Sorbing components on the same site name compete. The site and sorbed species are added as solids "
                 "(SOLID_<name>_ADS maps); a sorbed density > 0 makes the deposit occupy pore space (permeability law above).");
        {
            std::vector<std::string> fluidNames;
            for (const auto& q : T.comps) if (q.type != 2) fluidNames.push_back(q.name);
            if (!T.sorb.empty()) {   // C16: Excel-style
                ui::RecordTable<ThermalSorb> t;
                t.choiceStr("Component", [](ThermalSorb& S) -> std::string& { return S.comp; }, fluidNames)
                 .num("Capacity (mol/m3 bulk)", [](ThermalSorb& S) -> double& { return S.capacity; })
                 .num("b_ref (1/kPa | m3/mol)", [](ThermalSorb& S) -> double& { return S.b; })
                 .num("T_ref (K)", [](ThermalSorb& S) -> double& { return S.tRef; })
                 .num("dH ads (kJ/mol)", [](ThermalSorb& S) -> double& { return S.dH; }, 1e-3)
                 .num("Time const. (h)", [](ThermalSorb& S) -> double& { return S.tau; }, 1.0 / 3600.0)
                 .choiceFn("Phase", [](ThermalSorb& S) { return S.phase == 2 ? 1 : S.phase == 0 ? 2 : 0; },
                           [](ThermalSorb& S, int k) { S.phase = k == 1 ? 2 : k == 2 ? 0 : -1; }, {"Auto", "Gas (partial pressure)", "Aqueous"})
                 .num("Sorbed density (kg/m3)", [](ThermalSorb& S) -> double& { return S.adsRho; })
                 .text("Site", [](ThermalSorb& S) -> std::string& { return S.site; })
                 .flag("Initial equilibrium", [](ThermalSorb& S) -> bool& { return S.initEq; });
                t.onChange = [](ThermalSorb& S) { S.capacity = std::max(0.0, S.capacity); S.b = std::max(0.0, S.b); S.tau = std::max(1.0, S.tau); };
                const float h = std::min(ImGui::GetFrameHeightWithSpacing() * (float(T.sorb.size()) + 1.8f), ImGui::GetFontSize() * 18);
                if (t.draw("sorb", T.sorb, h, nullptr, "%.5g", true)) changed();
            }
            if (ImGui::Button(ICON_PLUS " Add sorbing component") && !fluidNames.empty()) {
                ThermalSorb S; S.comp = fluidNames.back(); S.capacity = 500; S.b = 1.0 / 6690; S.tRef = 308.15; S.dH = 15e3; S.tau = 86400;
                T.sorb.push_back(S); changed();
            }
            ImGui::SameLine();
            if (ImGui::Button("Load CO2-ECBM example")) { ImGui::EndTabItem(); ImGui::EndTabBar(); newCase(3, "ecbm"); return; }
            // single-component isotherms of the gas-phase entries at T_ref
            bool anyGas = false;
            for (const auto& S : T.sorb) anyGas = anyGas || S.phase != 0;
            if (anyGas && ImPlot::BeginPlot("##iso", ImVec2(-1, ImGui::GetFontSize() * 13))) {
                ImPlot::SetupAxes("Partial pressure (MPa)", "Sorbed (mol / m3 bulk)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                for (const auto& S : T.sorb) {
                    if (S.phase == 0 || !(S.b > 0) || !(S.capacity > 0)) continue;
                    std::vector<double> x, y;
                    for (int k = 0; k <= 60; ++k) { const double pk = k * 250.0; x.push_back(pk * 1e-3); y.push_back(S.capacity * S.b * pk / (1 + S.b * pk)); }
                    ImPlot::PlotLine(S.comp.c_str(), x.data(), y.data(), int(x.size()));
                }
                ImPlot::EndPlot();
            }
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ C13: fractures & EGS (EDFM)
    if (ui::tabItem("Fractures & EGS###fracegs")) {
        if (!egsEditor()) { ImGui::EndTabItem(); ImGui::EndTabBar(); return; }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ rock and heat loss
    if (ui::tabItem("Rock & heat loss")) {
        ui::hint("Rock enthalpy H_r = HEATCR dT + HEATCRT dT^2/2 per rock volume; conductivity THCONR x (1 - THCONSF S_g), harmonic "
                 "averaged between cells. Heat loss to the over / underburden: Vinsome-Westerveld semi-analytic profile per boundary face.");
        auto uni = [&](std::vector<double>& v, double def) { return v.empty() ? def : v[0]; };
        if (ui::beginProps("rock")) {
            double hc = uni(T.heatcr, U.f ? 35 : 2347), hct = uni(T.heatcrt, 0), tk = uni(T.thconr, U.f ? 24 : 150), sf = uni(T.thconsf, 0);
            std::snprintf(lab, sizeof lab, "Rock heat capacity HEATCR (%s)", U.hcr());
            if (ui::propD(lab, &hc, "%.4g", T.heatcr.size() > 1 ? "Per-cell array in the case: editing sets a uniform value" : "per unit rock volume")) { T.heatcr = {hc}; changed(); }
            std::snprintf(lab, sizeof lab, "Its T coefficient HEATCRT (%s/K)", U.hcr());
            if (ui::propD(lab, &hct, "%.4g")) { T.heatcrt = {hct}; changed(); }
            std::snprintf(lab, sizeof lab, "Conductivity THCONR (%s)", U.thc());
            if (ui::propD(lab, &tk, "%.4g", "Saturated-rock conductivity")) { T.thconr = {tk}; changed(); }
            if (ui::propD("Gas-saturation factor THCONSF (-)", &sf, "%.3f")) { T.thconsf = {std::clamp(sf, 0.0, 1.0)}; changed(); }
            if (ui::propBool("Heat loss to over / underburden", &T.heatLoss, "Vinsome-Westerveld; exported as a '-- RESSIM: HEATLOSS' line")) changed();
            if (T.heatLoss) {
                if (T.obCv <= 0) T.obCv = hc;
                if (T.obK <= 0) T.obK = tk;
                std::snprintf(lab, sizeof lab, "Cap / base rock heat capacity (%s)", U.hcr());
                if (ui::propD(lab, &T.obCv, "%.4g")) changed();
                std::snprintf(lab, sizeof lab, "Cap / base rock conductivity (%s)", U.thc());
                if (ui::propD(lab, &T.obK, "%.4g")) changed();
            }
            ui::endProps();
        }
        if (T.heatcr.size() > 1 || T.thconr.size() > 1) ui::hint("The deck defines per-cell HEATCR / THCONR arrays: view them in the 3D viewer (static properties).");
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ initial temperature
    if (ui::tabItem("Initial temperature")) {
        int mode = !T.tempvd.empty() ? 1 : 0;
        ImGui::RadioButton("Uniform / per cell (TEMPI)", &mode, 0); ImGui::SameLine();
        ImGui::RadioButton("Depth table (TEMPVD)", &mode, 1);
        if (mode == 1 && T.tempvd.empty()) { const double t0 = T.tempi.empty() ? (U.f ? 77 : 25) : T.tempi[0]; T.tempvd = {{0, t0}, {U.f ? 10000 : 3000, t0}}; changed(); }
        if (mode == 0 && !T.tempvd.empty()) { T.tempi = {T.tempvd.front().second}; T.tempvd.clear(); changed(); }
        if (mode == 0) {
            double t0 = T.tempi.empty() ? (U.f ? 77 : 25) : T.tempi[0];
            if (ui::beginProps("ti")) {
                std::snprintf(lab, sizeof lab, "Initial temperature (%s)", U.t());
                if (ui::propD(lab, &t0, "%.2f", T.tempi.size() > 1 ? "Per-cell TEMPI in the case: editing sets a uniform value" : nullptr)) { T.tempi = {t0}; changed(); }
                ui::endProps();
            }
        } else {
            int delRow = -1;
            if (ui::beginTable("tvd", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                std::snprintf(lab, sizeof lab, "Depth (%s)", U.f ? "ft" : "m"); ui::tableColumn(lab);
                std::snprintf(lab, sizeof lab, "T (%s)", U.t()); ui::tableColumn(lab);
                ui::tableColumn(""); ImGui::TableHeadersRow();
                for (int r = 0; r < int(T.tempvd.size()); ++r) {
                    ImGui::PushID(r); ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##d", &T.tempvd[size_t(r)].first, 0, 0, "%.2f")) changed();
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##t", &T.tempvd[size_t(r)].second, 0, 0, "%.2f")) changed();
                    ImGui::TableNextColumn(); if (T.tempvd.size() > 1 && ImGui::SmallButton(ICON_TRASH_2)) delRow = r;
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (delRow >= 0) { T.tempvd.erase(T.tempvd.begin() + delRow); changed(); }
            if (ImGui::SmallButton(ICON_PLUS " Row")) { auto r = T.tempvd.back(); r.first += U.f ? 100 : 30; T.tempvd.push_back(r); changed(); }
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ steam injection
    if (ui::tabItem("Steam injection")) {
        ui::hint("WINJTEMP: injected water / steam state per well from a report step on. Quality >= 0 with a temperature (or pressure) gives "
                 "saturated steam of that quality; the injection rate of the well is cold-water equivalent (CWE).");
        std::vector<std::string> inj;
        for (const auto& w : case_.wells) if (w.type != WellType::Producer) inj.push_back(w.name);
        if (!T.winj.empty()) {   // C16: Excel-style; T and P stay on the saturation curve when quality >= 0
            static std::vector<ThermalWinj> prev;
            ui::RecordTable<ThermalWinj> t;
            char hT[32], hP[32], hPs[40], hH[40];
            std::snprintf(hT, sizeof hT, "T (%s)", U.t()); std::snprintf(hP, sizeof hP, "P (%s)", U.p());
            std::snprintf(hPs, sizeof hPs, "psat (%s)", U.p()); std::snprintf(hH, sizeof hH, "h (%s)", U.hv());
            const int nsteps = std::max(0, int(case_.schedule.size()) - 1);
            t.choiceStr("Well", [](ThermalWinj& w) -> std::string& { return w.well; }, inj)
             .integer("From step", [](ThermalWinj& w) -> int& { return w.step; })
             .num("Quality (-, <0 water)", [](ThermalWinj& w) -> double& { return w.quality; })
             .optNum(hT, [](ThermalWinj& w) -> double& { return w.temp; }, [](ThermalWinj& w) { return w.temp < -1e8; },
                     [](ThermalWinj& w) { w.temp = -1e9; }, "Tsat at BHP")
             .num(hP, [](ThermalWinj& w) -> double& { return w.pres; })
             .computed(hPs, [U](const ThermalWinj& w) { return w.temp > -1e8 ? U.fromPa(steam::psat(U.toK(w.temp))) : 0.0; })
             .computed(hH, [U](const ThermalWinj& w) {
                 if (w.temp < -1e8) return 0.0;
                 const double K = U.toK(w.temp), ps = steam::psat(K);
                 const double q = std::clamp(w.quality, 0.0, 1.0), h = (1 - q) * steam::liquid(K, ps).h + q * steam::vapour(K, ps).h;
                 return U.f ? h / 2326.0 : h / 1000.0; });
            // keep T / P consistent on the saturation curve: whichever changed drives the other
            if (prev.size() != T.winj.size()) prev = T.winj;
            t.onChange = [nsteps](ThermalWinj& w) { w.step = std::clamp(w.step, 0, nsteps); w.quality = std::min(w.quality, 1.0); };
            if (t.draw("winj", T.winj, ImGui::GetFrameHeightWithSpacing() * (float(T.winj.size()) + 1.8f), nullptr, "%.4g", true)) {
                if (prev.size() == T.winj.size())
                    for (size_t r = 0; r < T.winj.size(); ++r) {
                        ThermalWinj& w = T.winj[r];
                        if (w.quality < 0) continue;
                        if (w.temp != prev[r].temp && w.temp > -1e8) w.pres = U.fromPa(steam::psat(U.toK(w.temp)));
                        else if (w.pres != prev[r].pres && w.pres > 0) w.temp = U.fromK(steam::tsat(U.toPa(w.pres)));
                    }
                changed();
            }
            prev = T.winj;
        }
        ImGui::BeginDisabled(inj.empty());
        if (ImGui::Button(ICON_PLUS " Add WINJTEMP")) {
            ThermalWinj w; w.well = inj.empty() ? "" : inj.front(); w.quality = 0.8; w.temp = U.f ? 464 : 240;
            w.pres = U.fromPa(steam::psat(U.toK(w.temp))); T.winj.push_back(w); changed();
        }
        ImGui::EndDisabled();
        if (inj.empty()) { ImGui::SameLine(); ImGui::TextDisabled("add a water injector in Wells first"); }
        // producer steam-trap (subcool) constraint
        ui::header("Producer steam trap (subcool) and steam-rate limit");
        ImGui::TextDisabled("SAGD / steamflood producers: the produced fluid is kept the given temperature difference below saturation - a dynamic "
                            "BHP floor Psat(T_hot + subcool), solved implicitly with the rate / BHP controls (STARS OPERATE STEAMTRAP). "
                            "Max steam caps the produced steam (cold-water equivalent, STARS OPERATE MAX STEAM) - also implicit. "
                            "Written as '-- RESSIM: STEAMTRAP' lines in the SCHEDULE section.");
        {
            std::vector<std::string> prod{"*"};
            for (const auto& w : case_.wells) if (w.type == WellType::Producer) prod.push_back(w.name);
            if (!T.traps.empty()) {   // C16: Excel-style
                ui::RecordTable<ThermalTrap> t;
                std::snprintf(lab, sizeof lab, "Subcool (%s, 0 = off)", U.f ? "degF" : "degC");
                const int nsteps = std::max(0, int(case_.schedule.size()) - 1);
                t.choiceStr("Producer", [](ThermalTrap& tp) -> std::string& { return tp.well; }, prod)
                 .integer("From step", [](ThermalTrap& tp) -> int& { return tp.step; })
                 .num(lab, [](ThermalTrap& tp) -> double& { return tp.subcool; })
                 .num(U.f ? "Max steam (STB/d CWE, 0 = off)" : "Max steam (sm3/d CWE, 0 = off)", [](ThermalTrap& tp) -> double& { return tp.steamMax; });
                t.onChange = [nsteps](ThermalTrap& tp) { tp.step = std::clamp(tp.step, 0, nsteps); tp.subcool = std::max(0.0, tp.subcool); tp.steamMax = std::max(0.0, tp.steamMax); };
                if (t.draw("traps", T.traps, ImGui::GetFrameHeightWithSpacing() * (float(T.traps.size()) + 1.8f), nullptr, "%.4g", true)) changed();
            }
            ImGui::BeginDisabled(prod.size() < 2);
            if (ImGui::Button(ICON_PLUS " Add steam trap")) { ThermalTrap tp; tp.well = prod.size() > 1 ? prod[1] : "*"; tp.subcool = U.f ? 18 : 10; T.traps.push_back(tp); changed(); }
            ImGui::EndDisabled();
        }
        // gas injection streams (air / enriched air / O2 / CO2 ...): WELLSTRE + WINJGAS
        ui::header("Gas injection streams (WELLSTRE / WINJGAS)");
        ImGui::TextDisabled("Gas injectors (air for combustion, oxygen / air for UCG, CO2 ...) inject the stream composition; rate in surface gas volume.");
        {
            int delS = -1;
            if (!T.streams.empty() && ui::beginTable("streams", nc + 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollX | ImGuiTableFlags_SizingFixedFit)) {
                ui::tableColumn("Stream", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 7);
                for (int k = 0; k < nc; ++k) ui::tableColumn(T.comps[size_t(k)].name.c_str(), ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 4.5f);
                ui::tableColumn("Sum", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 4);
                ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 2);
                ImGui::TableHeadersRow();
                for (int r = 0; r < int(T.streams.size()); ++r) {
                    CompStream& st = T.streams[size_t(r)];
                    if (int(st.z.size()) != nc) st.z.resize(size_t(nc), 0.0);
                    ImGui::PushID(r); ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    std::string nm = st.name;
                    if (ui::inputString("##sn", &nm) && !nm.empty()) {
                        for (auto& rs : case_.schedule) for (auto& ws : rs.injStreams) if (ws.second == st.name) ws.second = nm;
                        st.name = nm; changed();
                    }
                    double sum = 0;
                    for (int k = 0; k < nc; ++k) {
                        ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                        char id[16]; std::snprintf(id, sizeof id, "##z%d", k);
                        ImGui::BeginDisabled(T.comps[size_t(k)].type == 2);
                        if (ImGui::InputDouble(id, &st.z[size_t(k)], 0, 0, "%.4f")) { st.z[size_t(k)] = std::clamp(st.z[size_t(k)], 0.0, 1.0); changed(); }
                        ImGui::EndDisabled();
                        sum += st.z[size_t(k)];
                    }
                    ImGui::TableNextColumn();
                    if (std::fabs(sum - 1) > 1e-6) { if (ImGui::SmallButton("Norm") && sum > 0) { for (double& z : st.z) z /= sum; changed(); } }
                    else ImGui::TextDisabled("1");
                    ImGui::TableNextColumn(); if (ImGui::SmallButton(ICON_TRASH_2)) delS = r;
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (delS >= 0) { T.streams.erase(T.streams.begin() + delS); changed(); }
            if (ImGui::Button(ICON_PLUS " Add stream")) {
                CompStream st; st.name = "AIR" + std::to_string(T.streams.size() + 1); st.z.assign(size_t(nc), 0.0);
                for (int k = 0; k < nc; ++k) {   // air guess from component names
                    std::string n = T.comps[size_t(k)].name; for (auto& ch : n) ch = char(std::toupper((unsigned char)ch));
                    if (n == "O2") st.z[size_t(k)] = 0.21; else if (n == "N2") st.z[size_t(k)] = 0.79;
                }
                T.streams.push_back(st); changed();
            }
            // WINJGAS: stream per gas injector at report step 0 (later steps: Schedule / deck editor)
            std::vector<std::string> ginj;
            for (const auto& w : case_.wells) if (w.type == WellType::GasInjector) ginj.push_back(w.name);
            if (!ginj.empty() && !T.streams.empty() && !case_.schedule.empty()) {
                std::vector<std::string> sn = {"(none: shut)"};
                for (const auto& st : T.streams) sn.push_back(st.name);
                if (ui::beginProps("winjgas")) {
                    for (const auto& wn : ginj) {
                        auto& lst = case_.schedule[0].injStreams;
                        int cur = 0; size_t at = lst.size();
                        for (size_t i = 0; i < lst.size(); ++i) if (lst[i].first == wn) { at = i; for (size_t j = 0; j < T.streams.size(); ++j) if (T.streams[j].name == lst[i].second) cur = int(j) + 1; }
                        std::snprintf(lab, sizeof lab, "%s stream (WINJGAS, step 0)", wn.c_str());
                        if (ui::propCombo(lab, &cur, sn)) {
                            if (cur == 0) { if (at < lst.size()) lst.erase(lst.begin() + long(at)); }
                            else if (at < lst.size()) lst[at].second = T.streams[size_t(cur - 1)].name;
                            else lst.push_back({wn, T.streams[size_t(cur - 1)].name});
                            changed();
                        }
                    }
                    ui::endProps();
                }
            } else if (ginj.empty()) ImGui::TextDisabled("Add a gas injector in Wells to assign a stream.");
        }
        // steam table calculator
        ui::header("Steam tables (IAPWS-IF97)");
        static double tc = 0;
        if (tc == 0) tc = U.f ? 450 : 230;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        std::snprintf(lab, sizeof lab, "Temperature (%s)##st", U.t());
        ImGui::InputDouble(lab, &tc, 0, 0, "%.2f");
        {
            const double K = U.toK(tc), ps = steam::psat(K);
            const steam::Props L = steam::liquid(K, ps), V = steam::vapour(K, ps);
            const double hs = U.f ? 2326.0 : 1000.0;
            ImGui::Text("psat %.3f %s | liquid h %.1f, vapour h %.1f, latent %.1f %s | rho_l %.1f, rho_v %.3f kg/m3 | mu_l %.4f, mu_v %.5f cP",
                        U.fromPa(ps), U.p(), L.h / hs, V.h / hs, (V.h - L.h) / hs, U.hv(), 1 / L.v, 1 / V.v,
                        1e3 * steam::viscosity(K, 1 / L.v), 1e3 * steam::viscosity(K, 1 / V.v));
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ geomechanics (THM / THMC)
    if (ui::tabItem("Geomechanics")) {
        geomech::Input& G = T.geo;
        ui::hint("Linear thermo-poroelastic finite elements (8-node hexahedra) on the flow grid plus over / under / side-burden padding; "
                 "rollers on the sides and base, free or fixed top. Coupling: fixed-stress split (Kim, Tchelepi & Juanes), iterated per "
                 "time step - unconditionally stable. Reports displacement / surface heave, effective stresses, volumetric strain, "
                 "Mohr-Coulomb shear and tensile failure; porosity (and with a permeability law, permeability) follow the mechanics. "
                 "Exported as '-- RESSIM: GEOMECH' (OPM: MECH / YMODULE / PRATIO / BIOTCOEF / THERMEXR are read on import).");
        if (editGeomechInput(G, true)) changed();
        if (G.enabled && case_.grid.cornerPoint())
            ImGui::TextDisabled("Corner-point grid: the mechanics mesh is built from the corner-point geometry (isoparametric hexahedra on the "
                                "averaged cell corners, fault throws kept), padding follows the edge columns.");
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// C25: the finite-element geomechanics property sheet - shared by the thermal (THM) model and the sequential geomechanics
// of black-oil / compositional runs (Grid & rock > Reservoir geomechanics). Returns true when something changed.
bool App::editGeomechInput(geomech::Input& G, bool thermal) {
    bool ch = false;
        if (ui::beginProps("geo")) {
            if (ui::propBool(thermal ? "Geomechanics (THM)" : "Geomechanics (finite elements)", &G.enabled)) (ch = true);
            ImGui::BeginDisabled(!G.enabled);
            const std::vector<std::string> cpl = thermal ? std::vector<std::string>{"One-way (after each step)", "Iterative fixed-stress (two-way)"}
                                                   : std::vector<std::string>{"One-way (stress output, FE stress for faults / fractures)", "Two-way explicit (pore volume + permeability after each report step)",
                                                                              "Two-way iterative (report step repeated until flow and mechanics agree)"};
            if (ui::propCombo("Coupling", &G.coupling, cpl)) (ch = true);
            auto uniProp = [&](const char* name, std::vector<double>& v, double scale, const char* fmt, const char* tip) {
                double x = (v.empty() ? 0.0 : v[0]) * scale;
                if (ui::propD(name, &x, fmt, v.size() > 1 ? "Per-cell values in the case: editing sets a uniform value" : tip)) { v = {x / scale}; (ch = true); }
            };
            uniProp("Young's modulus E (GPa)", G.E, 1e-9, "%.4g", "drained; OPM YMODULE");
            uniProp("Poisson's ratio (-)", G.nu, 1.0, "%.3f", "OPM PRATIO");
            uniProp("Biot coefficient (-)", G.biot, 1.0, "%.3f", "OPM BIOTCOEF");
            uniProp("Linear thermal expansion (1/K)", G.alphaT, 1.0, "%.3g", "OPM THERMEXR");
            double ks = G.Ks * 1e-9;
            if (ui::propD("Grain bulk modulus Ks (GPa, 0 = incompressible)", &ks, "%.4g", "1/N = (alpha - phi0) / Ks")) { G.Ks = ks * 1e9; (ch = true); }
            if (ui::propI("Overburden layers", &G.overLayers, 0, 50)) (ch = true);
            if (ui::propD("Overburden thickness (m)", &G.overThick, "%.1f")) (ch = true);
            if (ui::propI("Underburden layers", &G.underLayers, 0, 50)) (ch = true);
            if (ui::propD("Underburden thickness (m)", &G.underThick, "%.1f")) (ch = true);
            if (ui::propI("Side-burden columns", &G.sidePad, 0, 20)) (ch = true);
            if (ui::propD("Side-burden width (m)", &G.sideThick, "%.1f")) (ch = true);
            double pe = G.padE * 1e-9;
            if (ui::propD("Padding rock E (GPa)", &pe, "%.4g")) { G.padE = pe * 1e9; (ch = true); }
            if (ui::propD("Padding rock Poisson's ratio", &G.padNu, "%.3f")) (ch = true);
            if (ui::propD("Padding rock thermal expansion (1/K)", &G.padAlphaT, "%.3g")) (ch = true);
            if (ui::propBool("Free top surface", &G.topFree, "off: vertical displacement fixed at the top")) (ch = true);
            double sg = G.sigmaVGrad * 1e-3, c0 = G.cohesion * 1e-6, t0 = G.tensile * 1e-6;
            if (ui::propD("Overburden stress gradient (kPa/m)", &sg, "%.2f")) { G.sigmaVGrad = sg * 1e3; (ch = true); }
            if (ui::propD("k0 = sigma'h / sigma'v (-)", &G.k0, "%.3f")) (ch = true);
            if (ui::propD("Cohesion (MPa)", &c0, "%.3f", "Mohr-Coulomb")) { G.cohesion = c0 * 1e6; (ch = true); }
            if (ui::propD("Friction angle (deg)", &G.frictionDeg, "%.1f")) (ch = true);
            if (ui::propD("Tensile strength (MPa)", &t0, "%.3f")) { G.tensile = t0 * 1e6; (ch = true); }
            if (ui::propI("Fixed-stress iterations per step", &G.maxIter, 1, 20)) (ch = true);
            if (ui::propD("Porosity tolerance (-)", &G.tol, "%.2g")) (ch = true);
            double ps = G.permStress * 1e5;
            if (ui::propD("Stress-dependent permeability (1/bar)", &ps, "%.3g",
                          "k / k0 = exp(-c (sigma'_m - sigma'_m0)) from the FE mean effective stress change; 0 = off. "
                          "Typical 1e-3 - 1e-2 1/bar (fractured / weak rock)")) { G.permStress = std::max(0.0, ps * 1e-5); (ch = true); }
            ImGui::EndDisabled();
            ui::endProps();
        }
        if (G.enabled) {   // C16 elasto-plasticity (Drucker-Prager matched to Mohr-Coulomb), softening, damage permeability
            ui::header("Plasticity and damage");
            if (ui::beginProps("geoplast")) {
                static const char* kPl[] = {"Elastic (failure is only reported)", "Elasto-plastic (Drucker-Prager)"};
                if (ui::propCombo("Constitutive model", &G.plasticity, kPl, 2,
                                  "Elasto-plastic: stresses are returned to the yield surface at every Gauss point (implicit return mapping, "
                                  "committed state per accepted time step); plastic strain, softening and damage permeability follow.")) (ch = true);
                ImGui::BeginDisabled(G.plasticity != 1);
                static const char* kMatch[] = {"Outer cone (compression corners)", "Inner cone (extension corners)", "Plane strain (exact 2D collapse)"};
                if (ui::propCombo("Mohr-Coulomb match", &G.dpMatch, kMatch, 3,
                                  "How the smooth Drucker-Prager cone is fitted to the Mohr-Coulomb cohesion and friction angle above. Plane "
                                  "strain reproduces the Mohr-Coulomb collapse load exactly in plane-strain problems.")) (ch = true);
                if (ui::propD("Dilation angle (deg)", &G.dilationDeg, "%.1f", "= friction angle: associated flow; 0: no plastic volume change (non-associated)")) (ch = true);
                double hM = G.hardening * 1e-6, cr = G.residualCohesion * 1e-6;
                if (ui::propD("Hardening modulus H (MPa)", &hM, "%.4g", "cohesion change per unit equivalent plastic strain; negative = softening")) { G.hardening = hM * 1e6; (ch = true); }
                if (ui::propD("Residual cohesion (MPa)", &cr, "%.3f", "lower limit of the cohesion under softening")) { G.residualCohesion = cr * 1e6; (ch = true); }
                if (ui::propD("Damage permeability coefficient (-)", &G.permPlastic, "%.4g",
                              "k / k0 = exp(coefficient x equivalent plastic strain), capped by the limit below (shear-enhanced damage permeability); 0 = off")) (ch = true);
                if (ui::propD("Damage permeability limit (-)", &G.permPlasticMax, "%.4g")) (ch = true);
                if (ui::propI("Max equilibrium iterations", &G.plasticIter, 5, 1000)) (ch = true);
                if (ui::propD("Equilibrium tolerance (-)", &G.plasticTol, "%.2g")) (ch = true);
                ImGui::EndDisabled();
                ui::endProps();
            }
            if (G.plasticity == 1) {
                const double pi = 3.14159265358979323846, sp = std::sin(G.frictionDeg * pi / 180), Np = (1 + sp) / std::max(1e-9, 1 - sp);
                ImGui::TextDisabled("Mohr-Coulomb: UCS = 2c sqrt(N_phi) = %.3g MPa, N_phi = %.3g. Results: plastic strain, plastic volumetric "
                                    "strain and damage permeability multiplier maps (Results -> Geomechanics).", 2 * G.cohesion * std::sqrt(Np) * 1e-6, Np);
            }
        }
        if (G.enabled) {   // boundary conditions per face (BCMECH-like): type + normal load increment
            ui::header("Mechanical boundary conditions");
            ImGui::TextDisabled("Per face of the mechanics grid (incl. padding). Roller: normal displacement fixed; free: traction-free or loaded by the "
                                "normal stress increment (compression +, e.g. surface load or tectonic stress change); fixed: all components; "
                                "rigid plate: one common normal displacement (prescribed, into the body +) or one net force (Mandel / triaxial tests).");
            static const char* kFace[] = {"West (x-)", "East (x+)", "South (y-)", "North (y+)", "Top surface", "Base"};
            static const char* kType[] = {"Default", "Roller", "Free / loaded", "Fixed", "Rigid plate"};
            if (ui::beginTable("mechbc", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                ui::tableColumn("Face"); ui::tableColumn("Condition"); ui::tableColumn("Normal stress increment (MPa)");
                ui::tableColumn("Plate control"); ui::tableColumn("Plate displacement (m)");
                ImGui::TableHeadersRow();
                for (int f = 0; f < 6; ++f) {
                    ImGui::PushID(f); ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::TextUnformatted(kFace[f]);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    int t = G.bcType[f] + 1;
                    if (ImGui::Combo("##t", &t, kType, 5)) { G.bcType[f] = t - 1; (ch = true); }
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    const bool plate = G.bcType[f] == 3;
                    const bool loadable = G.bcType[f] == 1 || (G.bcType[f] < 0 && f == 4 && G.topFree) || (plate && !G.plateDispControl[f]);
                    ImGui::BeginDisabled(!loadable);
                    double sMPa = G.bcStress[f] * 1e-6;
                    if (ImGui::InputDouble("##s", &sMPa, 0, 0, "%.4g")) { G.bcStress[f] = sMPa * 1e6; (ch = true); }
                    ImGui::EndDisabled();
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::BeginDisabled(!plate);
                    int pc = G.plateDispControl[f] ? 1 : 0;
                    static const char* kPc[] = {"Force (net stress)", "Displacement"};
                    if (ImGui::Combo("##pc", &pc, kPc, 2)) {
                        // at most one force-controlled plate: switching this one to force makes the others displacement-controlled
                        G.plateDispControl[f] = pc == 1;
                        if (pc == 0) for (int o = 0; o < 6; ++o) if (o != f && G.bcType[o] == 3) G.plateDispControl[o] = true;
                        (ch = true);
                    }
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::BeginDisabled(!G.plateDispControl[f]);
                    if (ImGui::InputDouble("##u", &G.bcDisp[f], 0, 0, "%.4g")) (ch = true);
                    ImGui::EndDisabled();
                    ImGui::EndDisabled();
                    ImGui::PopID();
                }
                ui::endTable();
            }
        }
        if (G.enabled && !G.E.empty() && !G.nu.empty()) {
            const double E = G.E[0], nu = G.nu[0], al = G.biot.empty() ? 1.0 : G.biot[0];
            const double K = geomech::bulkModulus(E, nu), Gs = geomech::shearModulus(E, nu), M = K + 4 * Gs / 3;
            ImGui::TextDisabled("K_dr %.3g GPa, G %.3g GPa, oedometric M %.3g GPa, fixed-stress c_p = alpha^2/K_dr = %.3g 1/bar, "
                                "uniaxial heave per 10 bar over 10 m: %.2f mm", K * 1e-9, Gs * 1e-9, M * 1e-9, al * al / K * 1e5, al * 1e6 * 10 / M * 1e3);
        }
    return ch;
}
