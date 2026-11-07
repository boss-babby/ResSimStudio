// Fluid model selector (black oil / compositional EOS / CO2 storage) and the editors of the compositional and the
// dedicated CO2STORE models (Fluid section). Every edit marks the compositional data "structured", so the exporter
// writes COMPS / EOS / CNAMES / TCRIT / ... / WELLSTRE / WINJGAS (or CO2STORE / SALINITY / RTEMPVD) from the case.
#include "app.h"
#include "ressim/examples.h"
#include "ressim/steam.h"
#include "widgets.h"
#include "implot.h"
#include "ressim/co2brine.h"
#include "ressim/eor.h"
#include "ressim/io.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <future>
#include <string>
#include <vector>

using namespace ressim;

namespace {
const char* kLib[] = {"C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6", "N2", "CO2", "H2S", "H2O", "H2"};
constexpr int kNumLib = int(sizeof(kLib) / sizeof(kLib[0]));

struct DispUnits {   // display units of the case unit system
    bool field;
    double tOut(double K) const { return field ? K * 1.8 - 459.67 : K - 273.15; }
    double tIn(double v) const { return field ? (v + 459.67) / 1.8 : v + 273.15; }
    double tcOut(double K) const { return field ? K * 1.8 : K; }
    double tcIn(double v) const { return field ? v / 1.8 : v; }
    double pOut(double Pa) const { return field ? Pa / 6894.757293168 : Pa / 1e5; }
    double pIn(double v) const { return field ? v * 6894.757293168 : v * 1e5; }
    double vOut(double m3mol) const { return field ? m3mol / (0.028316846592 / 453.59237) : m3mol * 1e3; }
    double vIn(double v) const { return field ? v * (0.028316846592 / 453.59237) : v * 1e-3; }
    double dOut(double m) const { return field ? m / 0.3048 : m; }
    double dIn(double v) const { return field ? v * 0.3048 : v; }
    const char* t() const { return field ? "\xC2\xB0" "F" : "\xC2\xB0" "C"; }
    const char* tc() const { return field ? "\xC2\xB0R" : "K"; }
    const char* p() const { return field ? "psia" : "bar"; }
    const char* v() const { return field ? "ft3/lb-mol" : "m3/kg-mol"; }
    const char* d() const { return field ? "ft" : "m"; }
};

eos::CompFluid spe5Preset() { return spe5CompFluid(); }

void normalize(std::vector<double>& z) {
    double s = 0; for (double v : z) s += std::max(0.0, v);
    if (s > 0) for (double& v : z) v = std::max(0.0, v) / s;
}

// keeps all per-component data consistent after inserting (at = index, ins = true) or removing a component
void reshapeComponents(CompositionalInput& C, int at, bool ins) {
    auto fix = [&](std::vector<double>& z) {
        if (ins) z.insert(z.begin() + std::min<long>(long(z.size()), at), 0.0);
        else if (at < int(z.size())) z.erase(z.begin() + at);
    };
    fix(C.fluid.zi);
    for (auto& t : C.fluid.zmfvd) for (auto& z : t.z) fix(z);
    for (auto& st : C.streams) fix(st.z);
    for (auto& m : C.fluid.eosRegions) {
        if (m.kij.empty()) continue;
        const int n = m.nc(), nOld = ins ? n - 1 : n + 1;
        std::vector<double> k(size_t(n) * size_t(n), 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                int oi = i, oj = j;
                if (ins) { if (i == at || j == at) continue; oi = i > at ? i - 1 : i; oj = j > at ? j - 1 : j; }
                else { oi = i >= at ? i + 1 : i; oj = j >= at ? j + 1 : j; }
                if (oi < nOld && oj < nOld) k[size_t(i) * size_t(n) + size_t(j)] = m.kij[size_t(oi) * size_t(nOld) + size_t(oj)];
            }
        m.kij = k;
    }
}

struct Envelope { std::vector<double> Tb, Pb, Td, Pd; double psat = 0; bool vapour = false; std::string comp; };
}  // namespace

// ============================================================================================ model selector
bool App::editFluidModel() {
    CompositionalInput& C = case_.comp;
    int model = case_.thermal.enabled ? 3 : !C.enabled ? 0 : C.co2store ? 2 : 1;
    ui::header("Fluid model");
    const int old = model;
    ImGui::RadioButton("Black oil (PVT tables)", &model, 0); ImGui::SameLine();
    ImGui::RadioButton("Compositional (EOS)", &model, 1); ImGui::SameLine();
    ImGui::RadioButton("CO2 storage (CO2STORE)", &model, 2);
    const bool co2Hover = ImGui::IsItemHovered();
    ImGui::SameLine();
    ImGui::RadioButton("Thermal (STARS-like, E300 THERMAL)", &model, 3);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Steam / hot-water processes (steamflood, CSS, SAGD): energy equation, IAPWS-IF97 steam tables, K-value\n"
                          "components (dead / distillable), T-dependent viscosity, rock heat capacity / conductivity, heat loss.");
    if (co2Hover)
        ImGui::SetTooltip("Dedicated CO2 - brine model (OPM Flow / ECLIPSE CO2STORE): CO2 and H2O are the components; mutual solubility\n"
                          "Spycher-Pruess (2003/2005) with Duan-Sun salting-out, CO2 density Redlich-Kwong, viscosity Fenghour,\n"
                          "brine density / viscosity Batzle-Wang, apparent molar volume of dissolved CO2 Garcia (2001).");
    if (model != old) {
        if (old == 3 && model != 3) case_.thermal.enabled = false;
        if (model == 3) {
            C.enabled = false;
            if (case_.thermal.comps.empty()) thermalPreset(case_, 25);
            case_.thermal.enabled = true; case_.thermal.structured = true;
            case_.phases.oil = case_.phases.water = case_.phases.gas = true;
            case_.phases.disgas = case_.phases.vapoil = false;
            notify("Thermal model: heavy-oil preset (HEAVY + LIGHT, OILVISCT, rock heat, heat loss) - edit below; set WINJTEMP for steam injectors");
            markDirty();
            editThermal();
            return false;
        }
        if (model == 0) {
            C.enabled = false;
            if (!case_.phases.oil) { case_.phases.oil = true; if (case_.props.swof.size() >= 2) for (auto& r : case_.props.swof) if (r.kro <= 0 && r.s < case_.props.swof.back().s) r.kro = std::pow(1 - r.s, 2); }
            notify("Black-oil fluid: PVT tables (the compositional data is kept and can be switched back on)");
        } else if (model == 1) {
            if (C.co2store || C.fluid.eosRegions.empty()) { C.fluid = spe5Preset(); C.streams = {{"GAS", {0.77, 0.20, 0.03, 0, 0, 0}}}; C.separator.clear(); }
            C.enabled = true; C.co2store = false; C.structured = true;
            case_.phases.oil = true; case_.phases.gas = true;
            notify("Compositional fluid: SPE5 six-component Peng-Robinson preset - edit the components, BIC and compositions below");
        } else {
            const double T = (C.fluid.resT > 285 && C.fluid.resT < 373.2) ? C.fluid.resT : 323.15;
            const double m = C.co2store ? C.fluid.salinity : 0.5;
            C = CompositionalInput();
            C.enabled = true; C.co2store = true; C.structured = true;
            C.fluid = co2StoreFluid(T, m);
            Phases& ph = case_.phases;
            ph.gas = true; ph.disgas = false; ph.vapoil = false; ph.polymer = false;
            if (ph.water) ph.oil = false; else ph.oil = true;   // brine = WATER (GAS + WATER) or OIL (OPM convention)
            PropsInput& p = case_.props;
            const double swc = p.swof.empty() ? 0.2 : p.swof.front().s;
            if (p.sgof.size() < 2) coreySGOF(p.sgof, swc, 0.05, 0.0, 2.0, 2.0, 0.8);
            if (ph.water && p.swof.size() < 2) coreySWOF(p.swof, swc, 0.0, 3.0, 2.0, 1.0, 1.0);
            notify(std::string("CO2 storage model: CO2 + brine (") + (ph.oil ? "OIL + GAS, OPM convention" : "WATER + GAS") + ") - gas injectors inject CO2, water injectors brine");
        }
        markDirty();
    }
    if (model == 0) return true;
    if (model == 3) { editThermal(); return false; }
    if (model == 1) editCompositional(); else editCo2Store();
    return false;
}

// ============================================================================================ compositional EOS
void App::editCompositional() {
    CompositionalInput& C = case_.comp;
    const DispUnits D{case_.units == UnitSystem::Field};
    if (C.fluid.eosRegions.empty()) { C.fluid = spe5Preset(); markDirty(); }
    eos::CompFluid& F = C.fluid;
    eos::Model& M = F.eosRegions[0];
    auto changed = [&]() { C.structured = true; markDirty(); };
    const int nc = M.nc();
    if (int(F.names.size()) != nc) { F.names.clear(); for (auto& q : M.comp) F.names.push_back(q.name); }
    if (!ImGui::BeginTabBar("comp")) return;

    // ---------------------------------------------------------------- components
    if (ui::tabItem("Components")) {
        ui::hint("Cubic equation of state (Peng-Robinson / Soave-Redlich-Kwong) with Peneloux volume shift and Lohrenz-Bray-Clark viscosity. "
                 "The ResSim compositional core is fully implicit in pressure and the component moles; phase split by stability test + flash.");
        if (ui::beginProps("eosp")) {
            int et = M.type == eos::EosType::SRK ? 1 : 0;
            static const char* types[] = {"Peng-Robinson", "Soave-Redlich-Kwong"};
            if (ui::propCombo("Equation of state (EOS)", &et, types, 2)) { for (auto& m : F.eosRegions) m.type = et ? eos::EosType::SRK : eos::EosType::PR; changed(); }
            if (M.type == eos::EosType::PR && ui::propBool("PR 1978 correction (PRCORR)", &F.pr78, "m(omega) of the 1978 Peng-Robinson for omega > 0.49")) { for (auto& m : F.eosRegions) m.pr78 = F.pr78; changed(); }
            double t = D.tOut(F.resT);
            char lab[64]; std::snprintf(lab, sizeof lab, "Reservoir temperature (%s)", D.t());
            if (ui::propD(lab, &t, "%.2f", "RTEMP: the compositional core is isothermal") && D.tIn(t) > 0) { F.resT = D.tIn(t); F.tempvd.clear(); changed(); }
            double st = D.tOut(F.stdT), sp = D.pOut(F.stdP);
            std::snprintf(lab, sizeof lab, "Standard temperature (%s)", D.t());
            if (ui::propD(lab, &st, "%.2f", "STCOND")) { F.stdT = D.tIn(st); changed(); }
            std::snprintf(lab, sizeof lab, "Standard pressure (%s)", D.p());
            if (ui::propD(lab, &sp, "%.4f", "STCOND") && sp > 0) { F.stdP = D.pIn(sp); changed(); }
            ui::endProps();
        }
        if (F.eosRegions.size() > 1) ui::hint("The deck has several EOS regions (NTEOS): the table edits region 1; the core uses region 1 everywhere.");
        const bool zi = F.zmfvd.empty();
        if (zi && int(F.zi.size()) != nc) F.zi.resize(size_t(nc), 0.0);
        char hTc[32], hPc[32], hVc[40];
        std::snprintf(hTc, sizeof hTc, "Tc (%s)", D.tc()); std::snprintf(hPc, sizeof hPc, "Pc (%s)", D.p()); std::snprintf(hVc, sizeof hVc, "Vc (%s)", D.v());
        int del = -1;
        {   // spreadsheet grid (paste a component table from a PVT report / Excel: name, z, MW, Tc, Pc, ...)
            ui::GridAccess ga;
            ga.count = [&] { return nc; };
            ga.cols = 10;
            ga.getText = [&](int i, int c) -> std::string {
                if (c == 0) return M.comp[size_t(i)].name;
                if (c == 1 && !zi) { char b[32]; std::snprintf(b, sizeof b, "%.5f", F.zmfvd[0].z.empty() ? 0.0 : F.zmfvd[0].z[0][size_t(i)]); return b; }
                return std::string();
            };
            ga.setText = [&](int i, int c, const std::string& nm) {
                if (c != 0 || nm.empty()) return;
                M.comp[size_t(i)].name = nm; F.names[size_t(i)] = nm;
                for (auto& m : F.eosRegions) if (i < m.nc()) m.comp[size_t(i)].name = nm;
            };
            ga.numGet = [&](int i, int c, double& v) -> bool {
                const eos::Component& q = M.comp[size_t(i)];
                switch (c) {
                    case 1: if (!zi) return false; v = F.zi[size_t(i)]; return true;
                    case 2: v = q.mw * 1e3; return true;
                    case 3: v = D.tcOut(q.tc); return true;
                    case 4: v = D.pOut(q.pc); return true;
                    case 5: v = D.vOut(q.vc); return true;
                    case 6: v = q.omega; return true;
                    case 7: v = q.sshift; return true;
                    case 8: v = q.parachor; return true;
                    case 9: v = D.vOut(q.vcVisc); return true;
                    default: return false;
                }
            };
            ga.numSet = [&](int i, int c, double x) {
                eos::Component& q = M.comp[size_t(i)];
                switch (c) {
                    case 1: F.zi[size_t(i)] = std::max(0.0, x); break;
                    case 2: if (x > 0) q.mw = x * 1e-3; break;
                    case 3: if (x > 0) q.tc = D.tcIn(x); break;
                    case 4: if (x > 0) q.pc = D.pIn(x); break;
                    case 5: q.vc = D.vIn(std::max(0.0, x)); if (q.vc > 0) q.zc = q.pc * q.vc / (eos::Rgas * q.tc); break;
                    case 6: q.omega = x; break;
                    case 7: q.sshift = x; break;
                    case 8: q.parachor = std::max(0.0, x); break;
                    case 9: q.vcVisc = D.vIn(std::max(0.0, x)); break;
                }
            };
            // rows are added / removed with the library buttons below (compositions, streams and BIC are reshaped)
            const float h = std::min(ImGui::GetFrameHeight() * (nc + 1.6f) + 4, ImGui::GetFontSize() * 22);
            if (ui::gridTable("comps", ga, {"Name", zi ? "z (ZI)" : "z (ZMFVD)", "MW (g/mol)", hTc, hPc, hVc, "\xCF\x89", "Shift", "Parachor", "Vc visc"}, h, nullptr, "%.5g", false))
                changed();
            if (nc > 1) {
                static int rm = 0;
                rm = std::max(0, std::min(rm, nc - 1));
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
                if (ImGui::BeginCombo("##rmc", M.comp[size_t(rm)].name.c_str())) {
                    for (int i = 0; i < nc; ++i) if (ImGui::Selectable(M.comp[size_t(i)].name.c_str(), i == rm)) rm = i;
                    ImGui::EndCombo();
                }
                ImGui::SameLine();
                if (ImGui::SmallButton("Remove component")) del = rm;
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Remove the component (compositions, streams and BIC are reshaped)");
            }
        }
        if (del >= 0) {
            for (auto& m : F.eosRegions) if (del < m.nc()) m.comp.erase(m.comp.begin() + del);
            F.names.erase(F.names.begin() + del);
            reshapeComponents(C, del, false);
            normalize(F.zi);
            changed();
        }
        if (zi) {
            double s = 0; for (double v : F.zi) s += v;
            ImGui::TextDisabled("Sum of z = %.6f", s);
            if (std::fabs(s - 1) > 1e-9) { ImGui::SameLine(); if (ImGui::SmallButton("Normalise")) { normalize(F.zi); changed(); } }
        }
        static int lib = 0;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
        ImGui::Combo("##lib", &lib, kLib, kNumLib);
        ImGui::SameLine();
        auto addComp = [&](const eos::Component& q) {
            for (auto& m : F.eosRegions) m.comp.push_back(q);
            F.names.push_back(q.name);
            reshapeComponents(C, int(M.comp.size()) - 1, true);
            changed();
        };
        if (ImGui::Button("Add from library")) {
            eos::Component q = eos::standardComponent(kLib[lib]);
            if (q.vc <= 0 && q.zc > 0) q.vc = q.zc * eos::Rgas * q.tc / q.pc;
            addComp(q);
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Critical properties from Poling-Prausnitz-O'Connell / NIST");
        ImGui::SameLine();
        if (ImGui::Button("Add pseudo-component")) {
            eos::Component q; q.name = "PS" + std::to_string(nc + 1); q.tc = 700; q.pc = 2.0e6; q.omega = 0.55; q.mw = 0.2; q.zc = 0.25; q.vc = q.zc * eos::Rgas * q.tc / q.pc;
            addComp(q);
        }
        ImGui::SameLine();
        bool reset = false;   // the fluid was replaced: M / nc are stale - restart next frame
        if (ImGui::Button("SPE5 preset")) { C.fluid = spe5Preset(); C.streams = {{"GAS", {0.77, 0.20, 0.03, 0, 0, 0}}}; changed(); reset = true; }
        ImGui::SameLine();
        if (ImGui::Button("SPE5 + CO2 (CO2 EOR)")) {
            reset = true;
            C.fluid = spe5Co2CompFluid(C.fluid.salinity > 0 ? C.fluid.salinity : 1.0);
            C.streams = {{"CO2", {0, 0, 0, 0, 0, 0, 1.0}}, {"GAS", {0.77, 0.20, 0.03, 0, 0, 0, 0}}};
            changed();
            notify("Compositional fluid: SPE5 + CO2 (7 components, CO2 BIC 0.105 - 0.125), CO2SOL on, injection stream CO2 = pure CO2");
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("SPE5 six-component fluid plus CO2 as a seventh EOS component (CO2 EOR / WAG). "
                                                      "Inject the 'CO2' stream with a gas injector (WINJGAS STREAM CO2).");
        // CO2SOL: CO2 is an ordinary EOS component that, in addition, dissolves in the aqueous phase (E300 CO2SOL)
        if (!reset) {
            int ico2 = -1;
            for (int i = 0; i < nc; ++i) { std::string u = M.comp[size_t(i)].name; for (auto& ch : u) ch = char(std::toupper((unsigned char)ch)); if (u == "CO2") ico2 = i; }
            ImGui::SeparatorText("CO2 solubility in water (CO2SOL)");
            ImGui::BeginDisabled(ico2 < 0 || !case_.phases.water);
            bool sol = C.fluid.co2sol;
            if (ImGui::Checkbox("CO2SOL - CO2 component also dissolves in the brine", &sol)) { C.fluid.co2sol = sol; changed(); }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("E300 CO2SOL: aqueous CO2 from fugacity equality f_CO2(HC, EOS) = x_aq H(p, T, salinity); H calibrated to "
                                  "Spycher-Pruess / Duan-Sun. Needs a component named CO2 and the WATER phase.");
            if (C.fluid.co2sol) {
                double m = C.fluid.salinity;
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
                if (ImGui::InputDouble("Salinity (mol NaCl / kg H2O, SALINITY)", &m, 0.1, 0.5, "%.3f") && m >= 0) { C.fluid.salinity = std::min(m, 6.0); changed(); }
                ImGui::TextDisabled("%.0f ppm NaCl. Initial brine: CO2-saturated at the local hydrocarbon composition (E300 default without RSWVD).",
                                    co2::ppmFromMolality(C.fluid.salinity));
            }
            ImGui::EndDisabled();
            if (ico2 < 0) ImGui::TextDisabled("Add CO2 from the library (or use 'SPE5 + CO2') to enable CO2SOL.");
            else if (!case_.phases.water) ImGui::TextDisabled("Enable the WATER phase to use CO2SOL.");
            ImGui::SeparatorText("Water vaporisation");
            ImGui::BeginDisabled(!case_.phases.water);
            bool wv = C.waterVap;
            if (ImGui::Checkbox("Water vapour in the gas phase (VAPWAT)", &wv)) { C.waterVap = wv; changed(); }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("Free-water Raoult law: the gas phase carries y_w = Psat(T)/p of water vapour (IAPWS-IF97) while free water\n"
                                  "remains; dry-out when it is used up. Matters for hot reservoirs, gas cycling and CO2 / gas injection\n"
                                  "(near-well drying). Exported as a '-- RESSIM: WATERVAP' line (E300 EOS mode has no VAPWAT).");
            if (C.waterVap) {
                const double T = C.fluid.resT, ps = T > 273.16 && T < 647 ? steam::psat(T) : 0.0;
                ImGui::TextDisabled("Psat(%.1f %s) = %.4g %s: y_w = %.3g %% at %.0f %s", D.tOut(T), D.t(), D.pOut(ps), D.p(),
                                    100.0 * ps / std::max(1.0, D.pIn(case_.equil.datumPressure)), case_.equil.datumPressure, D.p());
            }
            ImGui::EndDisabled();
        }
        ImGui::EndTabItem();
        if (reset) { ImGui::EndTabBar(); return; }
    }

    // ---------------------------------------------------------------- BIC
    if (ui::tabItem("Binary interaction")) {
        ui::hint("Binary interaction coefficients k_ij (BIC, symmetric). Edit the lower triangle; the diagonal is zero.");
        {   // full symmetric matrix as a spreadsheet grid: paste the whole k_ij matrix (or its lower triangle) at once
            ui::GridAccess ga;
            ga.count = [&] { return nc; };
            ga.cols = nc;
            ga.numGet = [&](int i, int j, double& v) -> bool { if (i == j) return false; v = M.k(i, j); return true; };
            ga.numSet = [&](int i, int j, double v) { if (i != j) for (auto& m : F.eosRegions) m.setK(i, j, v); };
            ga.getText = [](int, int) { return std::string("0"); };
            ga.rowLabel = [&](int i) { return M.comp[size_t(i)].name; };
            std::vector<std::string> names; for (int j = 0; j < nc; ++j) names.push_back(M.comp[size_t(j)].name);
            std::vector<const char*> cn; for (auto& x : names) cn.push_back(x.c_str());
            if (ui::gridTable("bic", ga, cn, std::min(ImGui::GetFrameHeight() * (nc + 1.6f) + 4, ImGui::GetFontSize() * 24), nullptr, "%.4f", false)) changed();
            ui::hint("Symmetric: editing k_ij also sets k_ji. A pasted lower triangle (blank upper part) is accepted.");
        }
        static double A = 0.18, theta = 1.2;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputDouble("A", &A, 0, 0, "%.3f"); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputDouble("\xCE\xB8", &theta, 0, 0, "%.2f"); ImGui::SameLine();
        if (ImGui::Button("Chueh-Prausnitz (hydrocarbon pairs)")) {
            auto vc = [](const eos::Component& q) { return q.vc > 0 ? q.vc : q.zc * eos::Rgas * q.tc / q.pc; };
            auto isHC = [](const std::string& n) { return n != "N2" && n != "CO2" && n != "H2S" && n != "H2O" && n != "H2"; };
            for (auto& m : F.eosRegions)
                for (int i = 1; i < nc; ++i)
                    for (int j = 0; j < i; ++j) {
                        if (!isHC(m.comp[size_t(i)].name) || !isHC(m.comp[size_t(j)].name)) continue;
                        const double a = std::cbrt(vc(m.comp[size_t(i)])), b = std::cbrt(vc(m.comp[size_t(j)]));
                        if (!(a > 0 && b > 0)) continue;
                        m.setK(i, j, A * (1 - std::pow(2 * std::sqrt(a * b) / (a + b), theta)));
                    }
            changed();
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("k_ij = A [1 - (2 (Vci Vcj)^(1/6) / (Vci^(1/3) + Vcj^(1/3)))^theta] for hydrocarbon pairs (non-hydrocarbons unchanged)");
        ImGui::SameLine();
        if (ImGui::Button("Zero all")) { for (auto& m : F.eosRegions) m.kij.clear(); changed(); }
        ImGui::EndTabItem();
    }

    // ---------------------------------------------------------------- initial composition
    if (ui::tabItem("Initial composition")) {
        int mode = F.zmfvd.empty() ? 0 : 1;
        if (ImGui::RadioButton("Uniform (ZI)", &mode, 0) && !F.zmfvd.empty()) {
            if (!F.zmfvd[0].z.empty()) F.zi = F.zmfvd[0].z[0];
            F.zmfvd.clear(); changed();
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Composition vs depth (ZMFVD)", &mode, 1) && F.zmfvd.empty()) {
            eos::ZmfTable t;
            std::vector<double> z = F.zi; z.resize(size_t(nc), 0.0); if (z == std::vector<double>(size_t(nc), 0.0)) z[0] = 1;
            const double d0 = case_.equil.datumDepth * (case_.units == UnitSystem::Field ? 0.3048 : 1.0);
            t.depth = {d0 - 100, d0 + 100}; t.z = {z, z};
            F.zmfvd = {t}; changed();
        }
        if (mode == 0) ui::hint("Uniform composition: edit z in the Components tab. Equilibration (EQUIL) integrates the EOS density from the datum.");
        else {
            ui::hint("ZMFVD: overall composition as a function of depth (linear interpolation, constant outside), per equilibration region. "
                     "Gas caps form where the composition is above its dew point.");
            static int reg = 0; reg = std::clamp(reg, 0, int(F.zmfvd.size()) - 1);
            eos::ZmfTable& T = F.zmfvd[size_t(reg)];
            int del = -1;
            if (ui::beginTable("zmf", nc + 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_ScrollX)) {
                char hd[32]; std::snprintf(hd, sizeof hd, "Depth (%s)", D.d());
                ui::tableColumn(hd);
                for (int j = 0; j < nc; ++j) ui::tableColumn(M.comp[size_t(j)].name.c_str());
                ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
                ImGui::TableHeadersRow();
                for (size_t r = 0; r < T.depth.size(); ++r) {
                    ImGui::PushID(int(r));
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    double d = D.dOut(T.depth[r]);
                    if (ImGui::InputDouble("##d", &d, 0, 0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue)) { T.depth[r] = D.dIn(d); changed(); }
                    T.z[r].resize(size_t(nc), 0.0);
                    for (int j = 0; j < nc; ++j) {
                        ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(j);
                        if (ImGui::InputDouble("##z", &T.z[r][size_t(j)], 0, 0, "%.5f", ImGuiInputTextFlags_EnterReturnsTrue)) changed();
                        ImGui::PopID();
                    }
                    ImGui::TableNextColumn();
                    if (T.depth.size() > 1 && ImGui::SmallButton("\xC3\x97")) del = int(r);
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (del >= 0) { T.depth.erase(T.depth.begin() + del); T.z.erase(T.z.begin() + del); changed(); }
            if (ImGui::Button("Add row")) { T.depth.push_back(T.depth.empty() ? 0.0 : T.depth.back() + 50); T.z.push_back(T.z.empty() ? std::vector<double>(size_t(nc), 1.0 / nc) : T.z.back()); changed(); }
            ImGui::SameLine();
            if (ImGui::Button("Normalise rows and sort by depth")) {
                for (auto& z : T.z) normalize(z);
                std::vector<size_t> ix(T.depth.size()); for (size_t i = 0; i < ix.size(); ++i) ix[i] = i;
                std::sort(ix.begin(), ix.end(), [&](size_t a, size_t b) { return T.depth[a] < T.depth[b]; });
                eos::ZmfTable s; for (size_t i : ix) { s.depth.push_back(T.depth[i]); s.z.push_back(T.z[i]); }
                T = s; changed();
            }
        }
        ImGui::EndTabItem();
    }

    // ---------------------------------------------------------------- injection streams and separator
    if (ui::tabItem("Injection & separator")) {
        ui::header("Injection streams (WELLSTRE)");
        int delS = -1;
        if (ui::beginTable("streams", nc + 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_ScrollX)) {
            ui::tableColumn("Stream");
            for (int j = 0; j < nc; ++j) ui::tableColumn(M.comp[size_t(j)].name.c_str());
            ui::tableColumn("Sum"); ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
            ImGui::TableHeadersRow();
            for (size_t s = 0; s < C.streams.size(); ++s) {
                CompStream& st = C.streams[s];
                st.z.resize(size_t(nc), 0.0);
                ImGui::PushID(int(s));
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                if (ui::inputString("##sn", &st.name)) changed();
                double sum = 0;
                for (int j = 0; j < nc; ++j) {
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(j);
                    if (ImGui::InputDouble("##z", &st.z[size_t(j)], 0, 0, "%.4f", ImGuiInputTextFlags_EnterReturnsTrue)) changed();
                    sum += st.z[size_t(j)];
                    ImGui::PopID();
                }
                ImGui::TableNextColumn();
                if (std::fabs(sum - 1) > 1e-9) { if (ImGui::SmallButton("norm")) { normalize(st.z); changed(); } } else ImGui::TextDisabled("1");
                ImGui::TableNextColumn(); if (ImGui::SmallButton("\xC3\x97")) delS = int(s);
                ImGui::PopID();
            }
            ui::endTable();
        }
        if (delS >= 0) { C.streams.erase(C.streams.begin() + delS); changed(); }
        if (ImGui::Button("Add stream")) { CompStream st; st.name = "INJ" + std::to_string(C.streams.size() + 1); st.z.assign(size_t(nc), 0.0); st.z[0] = 1; C.streams.push_back(st); changed(); }

        ui::header("Injection stream per well and report step (WINJGAS)");
        ui::hint("Gas injectors inject the assigned stream from the given report step on (until a later assignment). Without an assignment "
                 "the first component is injected.");
        std::vector<std::string> inj;
        for (auto& w : case_.wells) if (w.type == WellType::GasInjector) inj.push_back(w.name);
        int delStep = -1, delIdx = -1;
        if (ui::beginTable("winj", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("From report step"); ui::tableColumn("Well"); ui::tableColumn("Stream");
            ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
            ImGui::TableHeadersRow();
            double t = 0;
            for (size_t s = 0; s < case_.schedule.size(); ++s) {
                auto& A = case_.schedule[s].injStreams;
                for (size_t a = 0; a < A.size(); ++a) {
                    ImGui::PushID(int(s * 1000 + a));
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::Text("%zu (day %.0f)", s + 1, t);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::BeginCombo("##w", A[a].first.c_str())) {
                        for (auto& n : inj) if (ImGui::Selectable(n.c_str(), n == A[a].first)) { A[a].first = n; changed(); }
                        ImGui::EndCombo();
                    }
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::BeginCombo("##s", A[a].second.c_str())) {
                        for (auto& st : C.streams) if (ImGui::Selectable(st.name.c_str(), st.name == A[a].second)) { A[a].second = st.name; changed(); }
                        ImGui::EndCombo();
                    }
                    ImGui::TableNextColumn(); if (ImGui::SmallButton("\xC3\x97")) { delStep = int(s); delIdx = int(a); }
                    ImGui::PopID();
                }
                t += case_.schedule[s].length;
            }
            ui::endTable();
        }
        if (delStep >= 0) { auto& A = case_.schedule[size_t(delStep)].injStreams; A.erase(A.begin() + delIdx); changed(); }
        static int addStep = 1, addWell = 0, addStream = 0;
        if (!inj.empty() && !C.streams.empty() && !case_.schedule.empty()) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputInt("step", &addStep); addStep = std::clamp(addStep, 1, int(case_.schedule.size()));
            ImGui::SameLine(); addWell = std::clamp(addWell, 0, int(inj.size()) - 1); addStream = std::clamp(addStream, 0, int(C.streams.size()) - 1);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
            if (ImGui::BeginCombo("well", inj[size_t(addWell)].c_str())) { for (int i = 0; i < int(inj.size()); ++i) if (ImGui::Selectable(inj[size_t(i)].c_str(), i == addWell)) addWell = i; ImGui::EndCombo(); }
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
            if (ImGui::BeginCombo("stream", C.streams[size_t(addStream)].name.c_str())) { for (int i = 0; i < int(C.streams.size()); ++i) if (ImGui::Selectable(C.streams[size_t(i)].name.c_str(), i == addStream)) addStream = i; ImGui::EndCombo(); }
            ImGui::SameLine();
            if (ui::accentButton("Assign")) { case_.schedule[size_t(addStep - 1)].injStreams.push_back({inj[size_t(addWell)], C.streams[size_t(addStream)].name}); changed(); }
        } else ui::hint("Needs a gas injector (Wells), a stream and a schedule.");

        ui::header("Field separator (FIELDSEP)");
        ui::hint("Stages in series: the liquid of each stage feeds the next. ECLIPSE 300 semantics: the LAST stage is the stock tank "
                 "(its liquid is the stock-tank oil) - end the train with a stage at standard conditions for a separator + stock-tank "
                 "train. Empty = single flash to standard conditions.");
        if (!C.separator.empty()) {   // C16: Excel-style over a display-unit proxy (temperature offsets are not a pure scale)
            struct Stage { double t = 0, p = 0; };
            std::vector<Stage> rows;
            for (const auto& st : C.separator) rows.push_back({D.tOut(st.first), D.pOut(st.second)});
            char ht[32], hp[32]; std::snprintf(ht, sizeof ht, "T (%s)", D.t()); std::snprintf(hp, sizeof hp, "P (%s)", D.p());
            ui::RecordTable<Stage> t;
            t.num(ht, [](Stage& x) -> double& { return x.t; }).num(hp, [](Stage& x) -> double& { return x.p; });
            if (t.draw("sep", rows, ImGui::GetFrameHeightWithSpacing() * (float(rows.size()) + 1.8f), nullptr, "%.4g", true)) {
                C.separator.clear();
                for (const auto& x : rows) if (x.p > 0) C.separator.push_back({D.tIn(x.t), D.pIn(x.p)});
                changed();
            }
        }
        if (ImGui::Button("Add stage")) { C.separator.push_back(C.separator.empty() ? std::make_pair(D.tIn(case_.units == UnitSystem::Field ? 80.0 : 27.0), D.pIn(case_.units == UnitSystem::Field ? 300.0 : 20.0)) : std::make_pair(F.stdT, F.stdP)); changed(); }
        {
            const bool lastStd = !C.separator.empty() && std::fabs(C.separator.back().second - F.stdP) < 1e3 && std::fabs(C.separator.back().first - F.stdT) < 3.0;
            if (!C.separator.empty() && !lastStd) {
                ImGui::SameLine();
                if (ImGui::Button("Add stock-tank stage")) { C.separator.push_back({F.stdT, F.stdP}); changed(); }
                ui::hint("The last stage is not at standard conditions: its liquid at that pressure is counted as stock-tank oil (higher oil volume, lower GOR).");
            }
        }

        ui::header("Separator trains per well (SEPCOND / WSEPCOND)");
        ui::hint("Named multi-stage trains; wells assigned to a train report their surface rates (and are controlled) through its stages, "
                 "all other wells use the field separator.");
        int delTr = -1;
        for (size_t k = 0; k < C.sepTrains.size(); ++k) {
            auto& tr = C.sepTrains[k];
            ImGui::PushID(int(k) + 5000);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            if (ui::inputString("##trn", &tr.name)) changed();
            ImGui::SameLine(); ImGui::TextDisabled("%zu stage(s)", tr.stages.size());
            ImGui::SameLine(); if (ImGui::SmallButton("Delete train")) delTr = int(k);
            int delSt = -1;
            for (size_t i = 0; i < tr.stages.size(); ++i) {
                ImGui::PushID(int(i));
                double tt = D.tOut(tr.stages[i].first), pp = D.pOut(tr.stages[i].second);
                ImGui::Text("  stage %zu", i + 1); ImGui::SameLine();
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
                if (ImGui::InputDouble(D.t(), &tt, 0, 0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue)) { tr.stages[i].first = D.tIn(tt); changed(); }
                ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
                if (ImGui::InputDouble(D.p(), &pp, 0, 0, "%.3f", ImGuiInputTextFlags_EnterReturnsTrue) && pp > 0) { tr.stages[i].second = D.pIn(pp); changed(); }
                ImGui::SameLine(); if (ImGui::SmallButton("\xC3\x97")) delSt = int(i);
                ImGui::PopID();
            }
            if (delSt >= 0) { tr.stages.erase(tr.stages.begin() + delSt); changed(); }
            if (ImGui::SmallButton("Add stage##tr")) { tr.stages.push_back(tr.stages.empty() ? std::make_pair(D.tIn(case_.units == UnitSystem::Field ? 80.0 : 27.0), D.pIn(case_.units == UnitSystem::Field ? 300.0 : 20.0)) : std::make_pair(F.stdT, F.stdP)); changed(); }
            ImGui::PopID();
        }
        if (delTr >= 0) {
            const std::string nm = C.sepTrains[size_t(delTr)].name;
            C.sepTrains.erase(C.sepTrains.begin() + delTr);
            for (auto it = C.wellSep.begin(); it != C.wellSep.end();) it = it->second == nm ? C.wellSep.erase(it) : std::next(it);
            changed();
        }
        if (ImGui::Button("Add separator train")) {
            CompositionalInput::SepTrain tr; tr.name = "SEP" + std::to_string(C.sepTrains.size() + 1);
            tr.stages = {{D.tIn(case_.units == UnitSystem::Field ? 100.0 : 40.0), D.pIn(case_.units == UnitSystem::Field ? 500.0 : 35.0)}, {F.stdT, F.stdP}};
            C.sepTrains.push_back(tr); changed();
        }
        if (!C.sepTrains.empty()) {
            for (auto& w : case_.wells) {
                if (w.type != WellType::Producer) continue;
                ImGui::PushID(w.name.c_str());
                auto it = C.wellSep.find(w.name);
                const std::string cur = it == C.wellSep.end() ? std::string("field separator") : it->second;
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                if (ImGui::BeginCombo(w.name.c_str(), cur.c_str())) {
                    if (ImGui::Selectable("field separator", it == C.wellSep.end())) { C.wellSep.erase(w.name); changed(); }
                    for (auto& tr : C.sepTrains) if (ImGui::Selectable(tr.name.c_str(), cur == tr.name)) { C.wellSep[w.name] = tr.name; changed(); }
                    ImGui::EndCombo();
                }
                ImGui::PopID();
            }
        }
        ImGui::EndTabItem();
    }

    // ---------------------------------------------------------------- phase behaviour: envelope + flash calculator
    if (ui::tabItem("Phase behaviour")) {
        static Envelope env;
        static int src = 0;
        std::vector<std::string> srcNames = {F.zmfvd.empty() ? "Initial composition (ZI)" : "Initial composition (ZMFVD, first row)"};
        for (auto& st : C.streams) srcNames.push_back("Stream " + st.name);
        src = std::clamp(src, 0, int(srcNames.size()) - 1);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16);
        if (ImGui::BeginCombo("Composition", srcNames[size_t(src)].c_str())) {
            for (int i = 0; i < int(srcNames.size()); ++i) if (ImGui::Selectable(srcNames[size_t(i)].c_str(), i == src)) src = i;
            ImGui::EndCombo();
        }
        std::vector<double> z = src == 0 ? (F.zmfvd.empty() ? F.zi : (F.zmfvd[0].z.empty() ? std::vector<double>() : F.zmfvd[0].z[0])) : C.streams[size_t(src - 1)].z;
        z.resize(size_t(nc), 0.0); normalize(z);
        bool okData = true; for (auto& q : M.comp) okData = okData && q.tc > 0 && q.pc > 0 && q.mw > 0;
        ImGui::SameLine();
        if (ui::accentButton("Compute phase envelope") && okData) {
            env = Envelope(); env.comp = srcNames[size_t(src)];
            double tmax = 0; for (auto& q : M.comp) tmax = std::max(tmax, q.tc);
            const int NT = 70;
            const double t0 = 150.0, t1 = std::min(1100.0, tmax * 1.15);
            double gb = 0, gd = 0;
            for (int k = 0; k <= NT; ++k) {
                const double T = t0 + (t1 - t0) * k / NT;
                try {
                    const double pb = eos::saturationPressure(M, T, z, true, gb, 1e4, 1.5e8);
                    if (pb > 1.1e4) { env.Tb.push_back(D.tOut(T)); env.Pb.push_back(D.pOut(pb)); gb = pb; }
                    const double pd = eos::saturationPressure(M, T, z, false, gd, 1e4, 1.5e8);
                    if (pd > 1.1e4) { env.Td.push_back(D.tOut(T)); env.Pd.push_back(D.pOut(pd)); gd = pd; }
                } catch (...) {}
            }
            try {
                env.vapour = eos::isVapourLike(M, z, F.resT);
                env.psat = eos::saturationPressure(M, F.resT, z, !env.vapour, 0);
            } catch (...) { env.psat = 0; }
        }
        if (!okData) ui::hint("Every component needs Tc, Pc and MW.");
        if (!env.comp.empty()) {
            if (env.psat > 0) ImGui::Text("%s at reservoir temperature: %s pressure %.2f %s", env.comp.c_str(), env.vapour ? "dew-point" : "bubble-point", D.pOut(env.psat), D.p());
            else ImGui::Text("%s: single phase at reservoir temperature for all pressures checked", env.comp.c_str());
        }
        const float hPlot = ImGui::GetFontSize() * 17;
        if (ImPlot::BeginSubplots("##pb", 1, 2, ImVec2(-1, hPlot))) {
            if (ImPlot::BeginPlot("Phase envelope (P-T)")) {
                ImPlot::SetupAxes(D.t(), D.p(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlotSpec sb; sb.LineColor = ImVec4(0.30f, 0.65f, 1.0f, 1); sb.LineWeight = 2;
                ImPlotSpec sd; sd.LineColor = ImVec4(1.0f, 0.45f, 0.35f, 1); sd.LineWeight = 2;
                if (!env.Tb.empty()) ImPlot::PlotLine("Bubble point", env.Tb.data(), env.Pb.data(), int(env.Tb.size()), sb);
                if (!env.Td.empty()) ImPlot::PlotLine("Dew point", env.Td.data(), env.Pd.data(), int(env.Td.size()), sd);
                const double tr = D.tOut(F.resT);
                ImPlotSpec sr; sr.LineColor = ImVec4(1, 1, 1, 0.4f); sr.Flags = ImPlotItemFlags_NoLegend;
                ImPlot::PlotInfLines("##tres", &tr, 1, sr);
                if (env.psat > 0) { const double x = tr, y = D.pOut(env.psat); ImPlotSpec sm; sm.Marker = ImPlotMarker_Circle; sm.MarkerSize = 5; ImPlot::PlotScatter("Psat (T res)", &x, &y, 1, sm); }
                ImPlot::EndPlot();
            }
            // flash calculator
            static double fp = 0, ft = 0;
            if (fp <= 0) fp = D.pOut(case_.equil.datumPressure * (case_.units == UnitSystem::Field ? 6894.757293168 : 1e5));
            if (ft == 0) ft = D.tOut(F.resT);
            static eos::FlashResult fr; static bool frOk = false; static std::string frErr;
            if (ImPlot::BeginPlot("Flash: K-values")) {
                ImPlot::SetupAxes("component", "K = y / x", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                if (frOk && fr.phases == 2 && int(fr.K.size()) == nc) {
                    std::vector<double> xs(static_cast<size_t>(nc)); for (int i = 0; i < nc; ++i) xs[size_t(i)] = i + 1;
                    ImPlotSpec sk; sk.Marker = ImPlotMarker_Circle; sk.MarkerSize = 5; sk.LineWeight = 2;
                    ImPlot::PlotLine("K", xs.data(), fr.K.data(), nc, sk);
                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            ui::header("Flash calculator (stability test + two-phase flash)");
            char lp[32], lt[32]; std::snprintf(lp, sizeof lp, "P (%s)", D.p()); std::snprintf(lt, sizeof lt, "T (%s)", D.t());
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7); ImGui::InputDouble(lp, &fp, 0, 0, "%.3f"); ImGui::SameLine();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7); ImGui::InputDouble(lt, &ft, 0, 0, "%.2f"); ImGui::SameLine();
            if (ui::accentButton("Flash") && okData && fp > 0) {
                try {
                    const eos::TParams tp = eos::atTemperature(M, D.tIn(ft));
                    fr = eos::flash(M, tp, D.pIn(fp), z);
                    frOk = true; frErr.clear();
                } catch (const std::exception& e) { frOk = false; frErr = e.what(); }
            }
            if (!frErr.empty()) ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", frErr.c_str());
            if (frOk) {
                const double T = D.tIn(ft);
                auto mu = [&](const eos::PhaseEval& P, const std::vector<double>& x) { return P.vm > 0 ? eos::lbcViscosity(M, x, T, 1.0 / P.vm) * 1e3 : 0.0; };
                if (fr.phases == 2)
                    ImGui::Text("Two phases: vapour mole fraction %.4f  |  liquid %.1f kg/m3, %.4f cP  |  vapour %.1f kg/m3, %.4f cP  |  Z_L %.4f  Z_V %.4f",
                                fr.beta, fr.L.rho, mu(fr.L, fr.x), fr.V.rho, mu(fr.V, fr.y), fr.L.Z, fr.V.Z);
                else {
                    const eos::PhaseEval& P = fr.vapourLabel ? fr.V : fr.L;
                    ImGui::Text("Single phase (%s-like): %.1f kg/m3, %.4f cP, Z %.4f, MW %.2f g/mol", fr.vapourLabel ? "vapour" : "liquid", P.rho, mu(P, z), P.Z, P.mw * 1e3);
                }
                if (ui::beginTable("fr", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame)) {
                    for (const char* c : {"Component", "z", "x (liquid)", "y (vapour)", "K"}) ui::tableColumn(c);
                    ImGui::TableHeadersRow();
                    for (int i = 0; i < nc; ++i) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(M.comp[size_t(i)].name.c_str());
                        ImGui::TableNextColumn(); ImGui::Text("%.5f", z[size_t(i)]);
                        const bool two = fr.phases == 2 && int(fr.x.size()) == nc;
                        ImGui::TableNextColumn(); if (two) ImGui::Text("%.5f", fr.x[size_t(i)]); else ImGui::TextDisabled("-");
                        ImGui::TableNextColumn(); if (two) ImGui::Text("%.5f", fr.y[size_t(i)]); else ImGui::TextDisabled("-");
                        ImGui::TableNextColumn(); if (two) ImGui::Text("%.4g", fr.K[size_t(i)]); else ImGui::TextDisabled("-");
                    }
                    ui::endTable();
                }
            }
        }
        ImGui::EndTabItem();
    }

    // ---------------------------------------------------------------- CO2 / gas EOR tools
    if (ui::tabItem("CO2 EOR tools")) {
        ui::hint("Gas / CO2 EOR design on the EOS model: swelling test and P-X diagram (first-contact miscibility), multi-contact MMP "
                 "by the mixing-cell method (Ahmadi & Johns 2011, combined vaporising / condensing drive) with correlation cross-checks, "
                 "and a WAG schedule generator. CO2 is an ordinary EOS component - inject it as a WELLSTRE stream.");
        bool okData = true; for (auto& q : M.comp) okData = okData && q.tc > 0 && q.pc > 0 && q.mw > 0;
        // compositions
        static int oilSrc = 0, gasSrc = 0;
        std::vector<std::string> oilNames = {"Initial composition (ZI)"};
        for (size_t t = 0; t < F.zmfvd.size(); ++t) oilNames.push_back("ZMFVD table " + std::to_string(t + 1) + ", first row");
        int ico2 = -1;
        for (int i = 0; i < nc; ++i) { std::string u = M.comp[size_t(i)].name; for (auto& ch : u) ch = char(std::toupper((unsigned char)ch)); if (u == "CO2") ico2 = i; }
        std::vector<std::string> gasNames;
        for (auto& st : C.streams) gasNames.push_back("Stream " + st.name);
        if (ico2 >= 0) gasNames.push_back("Pure CO2");
        oilSrc = std::clamp(oilSrc, 0, int(oilNames.size()) - 1);
        auto combo = [](const char* lab, int& v, const std::vector<std::string>& names) {
            if (names.empty()) { ImGui::TextDisabled("%s: none", lab); return; }
            v = std::clamp(v, 0, int(names.size()) - 1);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
            if (ImGui::BeginCombo(lab, names[size_t(v)].c_str())) {
                for (int i = 0; i < int(names.size()); ++i) if (ImGui::Selectable(names[size_t(i)].c_str(), i == v)) v = i;
                ImGui::EndCombo();
            }
        };
        combo("Oil", oilSrc, oilNames); ImGui::SameLine();
        combo("Injection gas", gasSrc, gasNames); ImGui::SameLine();
        static double tEor = 0;
        if (tEor == 0) tEor = D.tOut(F.resT > 0 ? F.resT : 344.26);
        char ltl[40]; std::snprintf(ltl, sizeof ltl, "T (%s)##eor", D.t());
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(ltl, &tEor, 0, 0, "%.2f");
        std::vector<double> zo = oilSrc == 0 || F.zmfvd.empty() ? F.zi : (F.zmfvd[size_t(oilSrc - 1)].z.empty() ? std::vector<double>() : F.zmfvd[size_t(oilSrc - 1)].z[0]);
        zo.resize(size_t(nc), 0.0); normalize(zo);
        std::vector<double> zg(size_t(nc), 0.0);
        if (!gasNames.empty()) {
            gasSrc = std::clamp(gasSrc, 0, int(gasNames.size()) - 1);
            if (gasSrc < int(C.streams.size())) zg = C.streams[size_t(gasSrc)].z; else if (ico2 >= 0) zg[size_t(ico2)] = 1.0;
            zg.resize(size_t(nc), 0.0); normalize(zg);
        }
        const double TK = D.tIn(tEor);
        const bool canRun = okData && !gasNames.empty() && TK > 200;
        if (!okData) ui::hint("Every component needs Tc, Pc and MW.");
        if (gasNames.empty()) ui::hint("Define an injection stream (Injection & separator tab) or add CO2 to the components.");

        // ---- swelling test / P-X
        static eor::SwellingResult sw; static std::string swFor;
        ui::header("Swelling test and P-X diagram (first-contact miscibility)");
        ImGui::BeginDisabled(!canRun);
        if (ui::accentButton("Run swelling test")) {
            std::vector<double> al; for (int i = 0; i <= 38; ++i) al.push_back(i * 0.025);
            try { sw = eor::swellingTest(M, TK, zo, zg, al, F.stdT, F.stdP); } catch (const std::exception& e) { sw = eor::SwellingResult(); sw.message = e.what(); }
            swFor = oilNames[size_t(oilSrc)] + " + " + (gasNames.empty() ? std::string() : gasNames[size_t(gasSrc)]);
        }
        ImGui::EndDisabled();
        const bool field = case_.units == UnitSystem::Field;
        const double gasConv = field ? 5.614583 : 1.0;   // sm3/m3 -> scf/bbl
        const char* gasUnit = field ? "scf/bbl" : "sm3/m3";
        if (!sw.message.empty() && !sw.ok) ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", sw.message.c_str());
        if (sw.ok) {
            ImGui::Text("%s: saturation pressure of the oil %.1f %s; first-contact miscibility pressure (cricondenbar of the P-X diagram) %.1f %s at %.1f mol%% gas",
                        swFor.c_str(), D.pOut(sw.psatOil), D.p(), D.pOut(sw.fcmp), D.p(), 100 * sw.fcmpAlpha);
            std::vector<double> ab, pb, ad, pd, gx, sf;
            for (auto& r : sw.rows) {
                if (r.psat > 0) { (r.dew ? ad : ab).push_back(100 * r.alpha); (r.dew ? pd : pb).push_back(D.pOut(r.psat)); }
                if (r.psat > 0 && !r.dew) { gx.push_back(r.gasAddedStd * gasConv); sf.push_back(r.swelling); }
            }
            if (ImPlot::BeginSubplots("##swp", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 15))) {
                if (ImPlot::BeginPlot("P-X diagram")) {
                    ImPlot::SetupAxes("injection gas (mol %)", D.p(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec a; a.LineColor = ImVec4(0.30f, 0.65f, 1.0f, 1); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 3;
                    ImPlotSpec b; b.LineColor = ImVec4(1.0f, 0.45f, 0.35f, 1); b.LineWeight = 2; b.Marker = ImPlotMarker_Square; b.MarkerSize = 3;
                    if (!ab.empty()) ImPlot::PlotLine("Bubble point", ab.data(), pb.data(), int(ab.size()), a);
                    if (!ad.empty()) ImPlot::PlotLine("Dew point", ad.data(), pd.data(), int(ad.size()), b);
                    if (sw.fcmp > 0) { const double y = D.pOut(sw.fcmp); ImPlotSpec h; h.LineColor = ImVec4(1, 1, 1, 0.45f); h.Flags = ImPlotInfLinesFlags_Horizontal; ImPlot::PlotInfLines("FCMP", &y, 1, h); }
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Swelling factor")) {
                    char xl[48]; std::snprintf(xl, sizeof xl, "gas added (%s of original oil)", gasUnit);
                    ImPlot::SetupAxes(xl, "V_sat / V_sat,oil", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec a; a.LineColor = ImVec4(0.45f, 0.85f, 0.45f, 1); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 3;
                    if (!gx.empty()) ImPlot::PlotLine("Swelling factor", gx.data(), sf.data(), int(gx.size()), a);
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            if (ImGui::TreeNode("Swelling test table")) {
                if (ui::beginTable("swt", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_ScrollY,
                                      ImVec2(0, ImGui::GetFontSize() * 14))) {
                    char h1[40], h2[40]; std::snprintf(h1, sizeof h1, "Psat (%s)", D.p()); std::snprintf(h2, sizeof h2, "Gas added (%s)", gasUnit);
                    for (const char* c : {"Gas (mol %)", (const char*)h2, (const char*)h1, "Type", "Swelling factor", "Density (kg/m3)"}) ui::tableColumn(c);
                    ImGui::TableSetupScrollFreeze(0, 1); ImGui::TableHeadersRow();
                    for (auto& r : sw.rows) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::Text("%.1f", 100 * r.alpha);
                        ImGui::TableNextColumn(); ImGui::Text("%.1f", r.gasAddedStd * gasConv);
                        ImGui::TableNextColumn(); if (r.psat > 0) ImGui::Text("%.1f", D.pOut(r.psat)); else ImGui::TextDisabled("-");
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.psat > 0 ? (r.dew ? "dew" : "bubble") : "-");
                        ImGui::TableNextColumn(); if (r.psat > 0) ImGui::Text("%.4f", r.swelling); else ImGui::TextDisabled("-");
                        ImGui::TableNextColumn(); if (r.psat > 0) ImGui::Text("%.1f", r.rhoSat); else ImGui::TextDisabled("-");
                    }
                    ui::endTable();
                }
                ImGui::TreePop();
            }
        }

        // ---- MMP
        ui::header("Minimum miscibility pressure (multi-contact, mixing-cell method)");
        static double mLo = 0, mHi = 0; static int mContacts = 60;
        if (mLo <= 0) { mLo = D.pOut(eos::psiaToPa(1000)); mHi = D.pOut(eos::psiaToPa(7000)); }
        char l1[40], l2[40]; std::snprintf(l1, sizeof l1, "from (%s)##mmp", D.p()); std::snprintf(l2, sizeof l2, "to (%s)##mmp", D.p());
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(l1, &mLo, 0, 0, "%.0f"); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(l2, &mHi, 0, 0, "%.0f"); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputInt("contacts", &mContacts, 0); mContacts = std::clamp(mContacts, 10, 200);
        ImGui::SameLine();
        static eor::MmpResult mr; static double mrYM = 0, mrCQ = 0, mrT = 0; static bool mrCo2 = false;
        ImGui::BeginDisabled(!canRun || !(mHi > mLo && mLo > 0));
        if (ui::accentButton("Compute MMP")) {
            eor::MmpOptions o; o.pLo = D.pIn(mLo); o.pHi = D.pIn(mHi); o.contacts = mContacts; o.nP = 18;
            try { mr = eor::mixingCellMmp(M, TK, zo, zg, o); } catch (const std::exception& e) { mr = eor::MmpResult(); mr.message = e.what(); }
            mrT = TK; mrYM = eor::mmpYelligMetcalfe(TK);
            double c5 = 0, mw5 = 0, light = 0;
            for (int i = 0; i < nc; ++i) {
                const auto& q = M.comp[size_t(i)]; std::string u = q.name; for (auto& ch : u) ch = char(std::toupper((unsigned char)ch));
                if (u == "C1" || u == "N2" || u == "CH4") light += zo[size_t(i)];
                else if (q.mw >= 0.0715 && u != "CO2" && u != "H2S") { c5 += zo[size_t(i)]; mw5 += zo[size_t(i)] * q.mw * 1e3; }
            }
            mrCQ = c5 > 0 ? eor::mmpCronquist(TK, mw5 / c5, light) : 0;
            mrCo2 = ico2 >= 0 && zg[size_t(ico2)] > 0.9;
        }
        ImGui::EndDisabled();
        if (!mr.message.empty()) {
            if (mr.ok) ImGui::Text("MMP %.1f %s at %.1f %s (mixing-cell, %d flashes)", D.pOut(mr.mmp), D.p(), D.tOut(mrT), D.t(), mr.flashes);
            else ImGui::TextColored(ImVec4(1, 0.6f, 0.3f, 1), "%s", mr.message.c_str());
            if (mrCo2) ImGui::TextDisabled("Pure-CO2 correlations: Yellig-Metcalfe %.0f %s, Cronquist %.0f %s (cross-check only)", D.pOut(mrYM), D.p(), D.pOut(mrCQ), D.p());
            if (sw.ok && sw.fcmp > 0) ImGui::TextDisabled("First-contact miscibility pressure (swelling test) %.0f %s - the multi-contact MMP must lie below it.", D.pOut(sw.fcmp), D.p());
            if (!mr.p.empty() && ImPlot::BeginSubplots("##mmpp", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
                if (ImPlot::BeginPlot("Minimum tie-line length")) {
                    ImPlot::SetupAxes(D.p(), "TL_min", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    std::vector<double> px; for (double v : mr.p) px.push_back(D.pOut(v));
                    ImPlotSpec a; a.LineColor = ImVec4(0.30f, 0.65f, 1.0f, 1); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 3;
                    ImPlot::PlotLine("TL_min", px.data(), mr.tlMin.data(), int(px.size()), a);
                    if (mr.ok) { const double x = D.pOut(mr.mmp); ImPlotSpec v; v.LineColor = ImVec4(1, 0.8f, 0.3f, 0.8f); ImPlot::PlotInfLines("MMP", &x, 1, v); }
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Tie-line profile (last immiscible pressure)")) {
                    ImPlot::SetupAxes("cell (gas -> oil)", "tie-line length", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    if (!mr.profileTL.empty()) {
                        std::vector<double> cx; for (size_t i = 0; i < mr.profileTL.size(); ++i) cx.push_back(double(i));
                        ImPlotSpec a; a.LineColor = ImVec4(0.45f, 0.85f, 0.45f, 1); a.LineWeight = 1.5f;
                        ImPlot::PlotLine("TL", cx.data(), mr.profileTL.data(), int(cx.size()), a);
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
        }

        // ---- slim tube (1-D compositional displacement on the full simulator)
        ui::header("Slim-tube MMP (1-D compositional displacement)");
        ImGui::TextDisabled("Runs the compositional core on a 1-D tube (100 cells, 12 m, no residual saturations) at each pressure and "
                            "records the oil recovery after 1.2 HCPV injected. MMP = break-over of the recovery curve.");
        static double sLo = 0, sHi = 0; static int sN = 7, sCells = 100; static double sPv = 1.2;
        if (sLo <= 0) { sLo = D.pOut(eos::psiaToPa(1800)); sHi = D.pOut(eos::psiaToPa(4000)); }
        char s1[40], s2[40]; std::snprintf(s1, sizeof s1, "from (%s)##st", D.p()); std::snprintf(s2, sizeof s2, "to (%s)##st", D.p());
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(s1, &sLo, 0, 0, "%.0f"); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(s2, &sHi, 0, 0, "%.0f"); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4); ImGui::InputInt("pressures##st", &sN, 0); sN = std::clamp(sN, 3, 20); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4); ImGui::InputInt("cells##st", &sCells, 0); sCells = std::clamp(sCells, 20, 400); ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4); ImGui::InputDouble("HCPV##st", &sPv, 0, 0, "%.2f"); sPv = std::clamp(sPv, 0.3, 3.0);
        static std::future<eor::SlimTubeResult> stFut; static eor::SlimTubeResult st; static bool stRun = false;
        static std::chrono::steady_clock::time_point stT0;
        if (stRun && stFut.valid() && stFut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            try { st = stFut.get(); } catch (const std::exception& e) { st = eor::SlimTubeResult(); st.message = e.what(); }
            stRun = false;
        }
        ImGui::BeginDisabled(!canRun || stRun || !(sHi > sLo && sLo > 0));
        if (ui::accentButton("Run slim tube")) {
            std::vector<double> ps;
            for (int k = 0; k < sN; ++k) ps.push_back(D.pIn(sLo + (sHi - sLo) * k / double(sN - 1)));
            eor::SlimTubeOptions o; o.cells = sCells; o.pvInjected = sPv;
            const eos::Model Mc = M; const std::vector<double> zoc = zo, zgc = zg; const double Tc = TK;
            stFut = std::async(std::launch::async, [Mc, zoc, zgc, Tc, ps, o] { return eor::slimTubeMmp(Mc, Tc, zoc, zgc, ps, o); });
            stRun = true; stT0 = std::chrono::steady_clock::now();
        }
        ImGui::EndDisabled();
        if (stRun) {
            ImGui::SameLine();
            const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - stT0).count();
            ImGui::TextDisabled("running %d displacements ... %.0f s", sN, el);
        } else if (!st.message.empty() || st.ok) {
            if (st.ok) {
                if (st.mmpRf > 0) ImGui::Text("Slim-tube MMP: break-over %.0f %s; 90 %% recovery reached at %.0f %s", D.pOut(st.mmpBreak), D.p(), D.pOut(st.mmpRf), D.p());
                else ImGui::Text("Slim-tube MMP: break-over %.0f %s", D.pOut(st.mmpBreak), D.p());
            }
            else ImGui::TextColored(ImVec4(1, 0.6f, 0.3f, 1), "%s", st.message.c_str());
            if (st.ok && mr.ok) ImGui::TextDisabled("Mixing-cell MMP %.0f %s (difference %.1f %%)", D.pOut(mr.mmp), D.p(), 100.0 * (st.mmpBreak / mr.mmp - 1));
            if (!st.p.empty() && ImPlot::BeginSubplots("##stp", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
                if (ImPlot::BeginPlot("Recovery vs pressure")) {
                    ImPlot::SetupAxes(D.p(), "oil recovery (fraction)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    std::vector<double> px; for (double v : st.p) px.push_back(D.pOut(v));
                    ImPlotSpec a; a.LineColor = ImVec4(0.30f, 0.65f, 1.0f, 1); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 4;
                    ImPlot::PlotLine("RF", px.data(), st.rf.data(), int(px.size()), a);
                    if (st.mmpBreak > 0) { const double x = D.pOut(st.mmpBreak); ImPlotSpec v; v.LineColor = ImVec4(1, 0.8f, 0.3f, 0.8f); ImPlot::PlotInfLines("break-over MMP", &x, 1, v); }
                    if (mr.ok) { const double x = D.pOut(mr.mmp); ImPlotSpec v; v.LineColor = ImVec4(0.6f, 0.9f, 0.5f, 0.6f); ImPlot::PlotInfLines("mixing-cell MMP", &x, 1, v); }
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Recovery vs injected HCPV")) {
                    ImPlot::SetupAxes("HCPV injected", "oil recovery (fraction)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    for (size_t k = 0; k < st.profiles.size() && k < st.p.size(); ++k) {
                        const auto& pr = st.profiles[k]; if (pr.empty()) continue;
                        std::vector<double> x; for (size_t i = 0; i < pr.size(); ++i) x.push_back(sPv * double(i + 1) / double(pr.size()));
                        char nm[48]; std::snprintf(nm, sizeof nm, "%.0f %s", D.pOut(st.p[k]), D.p());
                        ImPlotSpec a; a.LineColor = ImPlot::GetColormapColor(int(k)); a.LineWeight = 1.5f;
                        ImPlot::PlotLine(nm, x.data(), pr.data(), int(x.size()), a);
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
        }

        // ---- WAG generator
        ui::header("WAG schedule generator");
        std::vector<std::string> inj;
        for (auto& w : case_.wells) if (w.type != WellType::Producer) inj.push_back(w.name);
        static int wInj = 0, wStream = 0;
        static eor::WagSpec ws;
        static bool wsInit = false;
        if (!wsInit) { ws.gasRate = field ? 5000 : 1.5e5; ws.waterRate = field ? 5000 : 800; ws.bhpMax = field ? 6000 : 400; wsInit = true; }
        if (inj.empty()) ImGui::TextDisabled("No injector in the case - add one in the Wells section.");
        else {
            std::vector<std::string> strNames; for (auto& st : C.streams) strNames.push_back(st.name);
            combo("Injector", wInj, inj); ImGui::SameLine();
            combo("Gas stream", wStream, strNames);
            auto num = [&](const char* lab, double* v, const char* fmt) { ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6); ImGui::InputDouble(lab, v, 0, 0, fmt); };
            num("start (day)", &ws.start, "%.1f"); ImGui::SameLine(); num("gas half-cycle (d)", &ws.gasDays, "%.1f"); ImGui::SameLine();
            num("water half-cycle (d)", &ws.waterDays, "%.1f"); ImGui::SameLine();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5); ImGui::InputInt("cycles", &ws.cycles, 0); ws.cycles = std::clamp(ws.cycles, 1, 200);
            char gr[40], wr[40], bh[40];
            std::snprintf(gr, sizeof gr, "gas rate (%s)", field ? "Mscf/d" : "sm3/d"); std::snprintf(wr, sizeof wr, "water rate (%s)", field ? "stb/d" : "sm3/d");
            std::snprintf(bh, sizeof bh, "max BHP (%s)", field ? "psia" : "bar");
            num(gr, &ws.gasRate, "%.1f"); ImGui::SameLine(); num(wr, &ws.waterRate, "%.1f"); ImGui::SameLine(); num(bh, &ws.bhpMax, "%.1f");
            num("taper (gas length x per cycle)", &ws.taper, "%.3f"); ImGui::SameLine(); num("water chase (d)", &ws.chaseWaterDays, "%.1f"); ImGui::SameLine();
            ImGui::Checkbox("gas first", &ws.gasFirst); ImGui::SameLine(); ImGui::Checkbox("shut after", &ws.shutAfter);
            double gasTot = 0, gl = ws.gasDays; for (int k = 0; k < ws.cycles; ++k) { gasTot += gl; gl *= ws.taper > 0 ? ws.taper : 1.0; }
            ImGui::TextDisabled("WAG ratio (water : gas time) %.2f : 1, WAG period ends on day %.1f; replaces the injector's controls from the start day.",
                                ws.waterDays * ws.cycles / std::max(1e-9, gasTot), ws.start + gasTot + ws.waterDays * ws.cycles + ws.chaseWaterDays);
            if (ui::accentButton("Generate WAG schedule")) {
                ws.well = inj[size_t(std::clamp(wInj, 0, int(inj.size()) - 1))];
                ws.stream = strNames.empty() ? std::string() : strNames[size_t(std::clamp(wStream, 0, int(strNames.size()) - 1))];
                std::string msg;
                const int n = eor::applyWagSchedule(case_, ws, &msg);
                if (n > 0) markDirty();
                notify(msg);
            }
        }
        ImGui::EndTabItem();
    }

    // ---------------------------------------------------------------- deck preview
    if (ui::tabItem("Diffusion & hysteresis")) {
        editTrapping(false);
        ImGui::EndTabItem();
    }
    if (ui::tabItem("Deck keywords")) {
        ui::hint(C.structured ? "Written from this data on export (E300 / OPM Flow compositional keywords):"
                              : "Not edited since the import: the original compositional keywords are exported verbatim. Preview of the structured form:");
        std::string txt;
        for (int s = 0; s < 4; ++s) txt += compositionalKeywordsText(case_, s);
        ImGui::InputTextMultiline("##ck", txt.data(), txt.size() + 1, ImVec2(-1, -1), ImGuiInputTextFlags_ReadOnly);
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ============================================================================================ CO2 storage
void App::editCo2Store() {
    CompositionalInput& C = case_.comp;
    eos::CompFluid& F = C.fluid;
    const DispUnits D{case_.units == UnitSystem::Field};
    auto changed = [&]() { C.structured = true; markDirty(); };
    ui::hint("Dedicated CO2 storage model (CO2STORE): two components (CO2, H2O) in a CO2-rich gas phase and an aqueous brine phase. "
             "Mutual solubility Spycher-Pruess (2003, 2005) with Duan-Sun (2003) salting-out; CO2 density Redlich-Kwong (Spycher-Pruess "
             "parameters), viscosity Fenghour et al. (1998); brine density / viscosity Batzle-Wang (1992); dissolved CO2 volume Garcia (2001). "
             "Validity: 12 - 100 \xC2\xB0" "C, up to 600 bar, NaCl up to 6 mol/kg.");
    if (!ImGui::BeginTabBar("co2")) return;
    if (ui::tabItem("Fluid & initial state")) {
        if (ui::beginProps("co2p")) {
            double t = D.tOut(F.resT);
            char lab[80]; std::snprintf(lab, sizeof lab, "Temperature (%s)", D.t());
            if (ui::propD(lab, &t, "%.2f", "Isothermal: RTEMPVD is written with this value") && D.tIn(t) > 273.15) { F.resT = D.tIn(t); F.tempvd.clear(); changed(); }
            static int su = 0;
            static const char* sus[] = {"mol NaCl / kg H2O", "ppm NaCl (mass)", "wt % NaCl"};
            ui::propCombo("Salinity unit", &su, sus, 3);
            double v = su == 0 ? F.salinity : su == 1 ? co2::ppmFromMolality(F.salinity) : co2::ppmFromMolality(F.salinity) / 1e4;
            if (ui::propD("Salinity (SALINITY)", &v, su == 0 ? "%.4f" : "%.1f", "NaCl-equivalent salinity of the brine") && v >= 0) {
                F.salinity = su == 0 ? v : su == 1 ? co2::molalityFromPpm(v) : co2::molalityFromPpm(v * 1e4);
                F.salinity = std::min(F.salinity, 6.0);
                changed();
            }
            int bp = case_.phases.oil ? 0 : 1;
            static const char* bps[] = {"OIL (OPM Flow CO2STORE convention)", "WATER (GAS + WATER)"};
            if (ui::propCombo("Brine phase keyword", &bp, bps, 2, "OPM Flow accepts both; OIL + GAS decks use SGOF, GAS + WATER decks SWFN + SGFN")) {
                PropsInput& p = case_.props;
                if (bp == 1 && case_.phases.oil) {   // SGOF krog -> SWOF krw (Sw = 1 - Sg)
                    std::vector<SatRow> sw;
                    for (size_t i = p.sgof.size(); i-- > 0;) sw.push_back({1 - p.sgof[i].s, p.sgof[i].kro, 0, 0});
                    p.swof = sw; case_.phases.oil = false; case_.phases.water = true;
                    case_.equil.owc = case_.equil.goc;
                } else if (bp == 0 && !case_.phases.oil) {   // SWOF krw -> SGOF krog
                    auto krw = [&](double s) {
                        if (p.swof.empty()) return 1 - s;
                        if (s <= p.swof.front().s) return p.swof.front().krw;
                        for (size_t i = 1; i < p.swof.size(); ++i) if (s <= p.swof[i].s) { const double a = (s - p.swof[i - 1].s) / (p.swof[i].s - p.swof[i - 1].s); return p.swof[i - 1].krw + a * (p.swof[i].krw - p.swof[i - 1].krw); }
                        return p.swof.back().krw;
                    };
                    for (auto& r : p.sgof) r.kro = krw(1 - r.s);
                    case_.phases.oil = true; case_.phases.water = false;
                    case_.equil.goc = case_.equil.owc;
                }
                changed();
            }
            const double ul = case_.units == UnitSystem::Field ? 0.3048 : 1.0;
            (void)ul;
            double& contact = case_.phases.oil ? case_.equil.goc : case_.equil.owc;
            char lc[80]; std::snprintf(lc, sizeof lc, "CO2 - brine contact depth (%s)", Units::of(case_.units).len);
            if (ui::propD(lc, &contact, "%.2f", "Initial CO2 above this depth (EQUIL GOC for OIL + GAS, item 3 for GAS + WATER); below: CO2-free brine")) changed();
            char ld[80]; std::snprintf(ld, sizeof ld, "Datum depth (%s)", Units::of(case_.units).len);
            if (ui::propD(ld, &case_.equil.datumDepth, "%.2f")) changed();
            char lp[80]; std::snprintf(lp, sizeof lp, "Datum pressure (%s)", Units::of(case_.units).p);
            if (ui::propD(lp, &case_.equil.datumPressure, "%.3f")) changed();
            ui::endProps();
        }
        ui::hint("Wells: gas injectors inject CO2, water injectors inject brine; producer ORAT / WRAT / LRAT limit the brine rate, GRAT the CO2 rate. "
                 "Trapping & interactions: hysteresis (residual), diffusion, salt precipitation, mineral reactions, structural inventory. "
                 "Results: five-class CO2 inventory (mobile / structural / residual / dissolved / mineral) and the XCO2, SGTRAP, SGMAX, SSALT, MINCO2, PORFAC maps.");
        // properties at the datum pressure
        const co2::Fluid cf = co2::makeFluid(F.resT, F.salinity);
        const double pd = case_.equil.datumPressure * (case_.units == UnitSystem::Field ? 6894.757293168 : 1e5);
        if (pd > 1e5) {
            double x, y; co2::solubility(cf, pd, x, y);
            ui::header("Properties at the datum pressure");
            if (ui::beginTable("co2tab", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                auto row = [&](const char* k, const char* fmt, double v) { ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(k); ImGui::TableNextColumn(); ImGui::Text(fmt, v); };
                row("CO2 solubility (mol/kg H2O)", "%.4f", co2::co2Molality(cf, pd));
                row("CO2 mole fraction in brine x_CO2", "%.5f", x);
                row("H2O mole fraction in the CO2 phase y_H2O", "%.5f", y);
                const double rc = co2::co2Density(cf, pd);
                row("CO2 density (kg/m3)", "%.1f", rc);
                row("CO2 viscosity (cP)", "%.5f", co2::co2Viscosity(F.resT, rc) * 1e3);
                row("Brine density, CO2-free (kg/m3)", "%.1f", co2::brineDensity(cf, pd));
                row("Brine density, CO2-saturated (kg/m3)", "%.1f", co2::aqueousDensity(cf, pd, x));
                row("Brine viscosity (cP)", "%.4f", co2::brineViscosity(cf) * 1e3);
                row("Salinity (ppm NaCl)", "%.0f", co2::ppmFromMolality(F.salinity));
                ui::endTable();
            }
        }
        ImGui::EndTabItem();
    }
    if (ui::tabItem("Property plots")) {
        const co2::Fluid cf = co2::makeFluid(F.resT, F.salinity), cw = co2::makeFluid(F.resT, 0.0);
        const int N = 80;
        const double pmax = 4e7;
        std::vector<double> P, m, m0, yw, rc, rb, rbs, mc;
        for (int i = 1; i <= N; ++i) {
            const double p = pmax * i / N;
            double x, y; co2::solubility(cf, p, x, y);
            P.push_back(D.pOut(p)); m.push_back(co2::co2Molality(cf, p)); m0.push_back(co2::co2Molality(cw, p)); yw.push_back(y * 1e3);
            const double r = co2::co2Density(cf, p);
            rc.push_back(r); rb.push_back(co2::brineDensity(cf, p)); rbs.push_back(co2::aqueousDensity(cf, p, x)); mc.push_back(co2::co2Viscosity(F.resT, r) * 1e3);
        }
        std::vector<double> S, ms;
        const double pd = std::max(1e6, case_.equil.datumPressure * (case_.units == UnitSystem::Field ? 6894.757293168 : 1e5));
        for (int i = 0; i <= 30; ++i) { const double s = 6.0 * i / 30; S.push_back(s); ms.push_back(co2::co2Molality(co2::makeFluid(F.resT, s), pd)); }
        char tl[64]; std::snprintf(tl, sizeof tl, "T = %.1f %s, salinity %.3g mol/kg", D.tOut(F.resT), D.t(), F.salinity);
        ImGui::TextDisabled("%s", tl);
        if (ImPlot::BeginSubplots("##co2plots", 2, 2, ImVec2(-1, std::max(ImGui::GetContentRegionAvail().y, ImGui::GetFontSize() * 26)))) {
            ImPlotSpec a; a.LineWeight = 2; a.LineColor = ImVec4(0.30f, 0.65f, 1.0f, 1);
            ImPlotSpec b; b.LineWeight = 2; b.LineColor = ImVec4(0.6f, 0.6f, 0.6f, 1);
            ImPlotSpec c; c.LineWeight = 2; c.LineColor = ImVec4(1.0f, 0.45f, 0.35f, 1);
            ImPlotSpec d; d.LineWeight = 2; d.LineColor = ImVec4(0.35f, 0.85f, 0.55f, 1);
            if (ImPlot::BeginPlot("CO2 solubility in brine")) {
                ImPlot::SetupAxes(D.p(), "mol CO2 / kg H2O", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y2, "y_H2O (x 1000)", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                ImPlot::PlotLine("This brine", P.data(), m.data(), N, a);
                ImPlot::PlotLine("Pure water", P.data(), m0.data(), N, b);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                ImPlot::PlotLine("H2O in CO2 (y x 1000)", P.data(), yw.data(), N, c);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Densities")) {
                ImPlot::SetupAxes(D.p(), "kg/m3", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotLine("CO2", P.data(), rc.data(), N, c);
                ImPlot::PlotLine("Brine", P.data(), rb.data(), N, a);
                ImPlot::PlotLine("CO2-saturated brine", P.data(), rbs.data(), N, d);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("CO2 viscosity")) {
                ImPlot::SetupAxes(D.p(), "cP", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotLine("CO2 (Fenghour)", P.data(), mc.data(), N, c);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Salting-out (datum pressure)")) {
                ImPlot::SetupAxes("NaCl (mol/kg)", "mol CO2 / kg H2O", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::PlotLine("Solubility", S.data(), ms.data(), int(S.size()), a);
                const double sx = F.salinity, sy = co2::co2Molality(cf, pd);
                ImPlotSpec sm; sm.Marker = ImPlotMarker_Circle; sm.MarkerSize = 5;
                ImPlot::PlotScatter("This brine", &sx, &sy, 1, sm);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::EndTabItem();
    }
    if (ui::tabItem("Trapping & interactions")) {
        editTrapping(true);
        ImGui::EndTabItem();
    }
    if (ui::tabItem("Deck keywords")) {
        ui::hint(C.structured ? "Written from this data on export (OPM Flow CO2STORE keywords):"
                              : "Not edited since the import: the original keywords are exported verbatim. Preview of the structured form:");
        std::string txt;
        for (int s = 0; s < 4; ++s) txt += compositionalKeywordsText(case_, s);
        ImGui::InputTextMultiline("##ck", txt.data(), txt.size() + 1, ImVec2(-1, -1), ImGuiInputTextFlags_ReadOnly);
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}
