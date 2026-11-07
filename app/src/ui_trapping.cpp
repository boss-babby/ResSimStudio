// Trapping mechanisms and their interactions (round C5):
//   * residual trapping - gas relative-permeability hysteresis (SATOPTS HYSTER + EHYSTR, IMBNUM or Land Sgt,max)
//   * solubility trapping - CO2 dissolution (always on in CO2STORE) + molecular diffusion (DIFFUSE / DIFFCWAT / DIFFCGAS / DIFFC)
//   * salt precipitation (dry-out) - PRECSALT / VAPWAT / SALTSOL / PERMFACT or Verma-Pruess
//   * mineral trapping - kinetic carbonation reactions ('-- RESSIM: MINERAL')
//   * structural trapping - inventory split of the free gas held below seals ('-- RESSIM: STRUCTTOL')
// The hysteresis editor is shared by the black-oil SCAL section and the compositional / CO2STORE trapping tabs.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "ressim/co2brine.h"
#include "ressim/compsim.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace ressim;

namespace {
// SGOF-style table (Sg, krg) of SATNUM region r -> compositional saturation tables (gas part only)
comp::SatTables gasTable(const Case& c, int r) {
    comp::SatTables t;
    const SatRegion R = c.props.satRegion(r);
    for (const SatRow& row : R.sgof) { t.sg.push_back(row.s); t.krg.push_back(row.krw); }
    return t;
}
ImPlotSpec lineSpec(ImU32 col, float w = 2.0f) { ImPlotSpec s; s.LineColor = ImGui::ColorConvertU32ToFloat4(col); s.LineWeight = w; return s; }
}  // namespace

// ------------------------------------------------------------------------------------------------ residual trapping
void App::editHysteresis(bool compCore) {
    PropsInput& p = case_.props;
    CompositionalInput& C = case_.comp;
    auto changedComp = [&]() { C.structured = true; markDirty(); };
    ImGui::TextWrapped("Residual (capillary) trapping: when brine re-imbibes a CO2 / gas plume, part of the gas is snapped off and stays immobile. "
                       "The gas relative permeability follows a scanning curve from the historical maximum gas saturation Sg,max down to the "
                       "trapped saturation Sgt given by Land's relation 1/Sgt* - 1/Sg,max* = C.");
    if (ui::beginProps("hyst", 0.5f)) {
        if (ui::propBool("Hysteresis (SATOPTS HYSTER)", &p.hyst.on, "Switches relative-permeability hysteresis on (EHYSTR is written with it)")) markDirty();
        int model = p.hyst.model >= 2 ? 1 : 0;
        static const char* models[] = {"Carlson (EHYSTR model 0)", "Killough (EHYSTR model 2)"};
        if (ui::propCombo("Non-wetting phase model", &model, models, 2,
                          "Carlson: the scanning curve is the drainage curve shifted by the trapped gas (parallel).\n"
                          "Killough: the imbibition curve is mapped between Sgt and Sg,max and scaled to meet the drainage curve.")) {
            p.hyst.model = model == 1 ? 2 : 0; markDirty();
        }
        if (ui::propD("Curvature parameter (EHYSTR item 1)", &p.hyst.curv, "%.3g", "Capillary-pressure scanning curvature (written for ECLIPSE / OPM; Pc hysteresis is not modelled)")) markDirty();
        int src = C.trap.sgtMax > 0 ? 1 : 0;
        static const char* srcs[] = {"Imbibition table (IMBNUM)", "Land model: maximum trapped saturation"};
        if (compCore && ui::propCombo("Trapped saturation from", &src, srcs, 2, "IMBNUM: critical gas saturation of the imbibition table. Land: Sgt at Sg,max = end point.")) {
            C.trap.sgtMax = src == 1 ? 0.25 : -1; changedComp();
        }
        if (compCore && src == 1) {
            if (ui::propSliderD("Sgt,max (Land, at the maximum Sg)", &C.trap.sgtMax, 0.01, 0.6, "%.3f", "Written as '-- RESSIM: SGTMAX'")) changedComp();
        }
        ui::endProps();
    }
    // imbibition table selection (IMBNUM for all cells)
    const int nsr = p.numSatRegions();
    GridInput& g = case_.grid;
    int imbTab = -1;
    if (!g.imbnum.empty()) imbTab = std::clamp(g.imbnum[0] - 1, 0, nsr - 1);
    {
        ImGui::AlignTextToFramePadding(); ImGui::TextUnformatted("IMBNUM (imbibition table):");
        ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
        std::string cur = imbTab < 0 ? std::string("none (drainage)") : "table " + std::to_string(imbTab + 1);
        if (ImGui::BeginCombo("##imb", cur.c_str())) {
            if (ImGui::Selectable("none (drainage)", imbTab < 0)) { g.imbnum.clear(); markDirty(); }
            for (int r = 0; r < nsr; ++r) {
                const std::string l = "table " + std::to_string(r + 1) + (r == 0 ? " (= SATNUM 1)" : "");
                if (ImGui::Selectable(l.c_str(), imbTab == r)) { g.imbnum.assign(size_t(g.nx) * g.ny * g.nz, r + 1); markDirty(); }
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        if (ImGui::Button("Add imbibition table")) {   // copy of table 1 with a raised critical gas saturation
            SatRegion R = p.satRegion(0);
            const double sgcr = 0.25;
            if (R.sgof.size() >= 2) {
                // rescale the krg curve so that it starts at sgcr and ends at the drainage end point
                const SatRegion D = p.satRegion(0);
                const double s0 = D.sgof.front().s, s1 = D.sgof.back().s;
                for (size_t i = 0; i < R.sgof.size(); ++i) {
                    const double s = R.sgof[i].s;
                    const double sd = s <= sgcr ? s0 : s0 + (s - sgcr) * (s1 - s0) / std::max(1e-9, s1 - sgcr);
                    double kr = 0;
                    for (size_t j = 1; j < D.sgof.size(); ++j)
                        if (sd <= D.sgof[j].s) { const double a = (sd - D.sgof[j - 1].s) / std::max(1e-12, D.sgof[j].s - D.sgof[j - 1].s); kr = D.sgof[j - 1].krw + a * (D.sgof[j].krw - D.sgof[j - 1].krw); break; }
                    R.sgof[i].krw = s <= sgcr ? 0.0 : kr;
                }
            }
            p.satExtra.push_back(R);
            g.imbnum.assign(size_t(g.nx) * g.ny * g.nz, p.numSatRegions());
            if (g.satnum.empty()) g.satnum.assign(size_t(g.nx) * g.ny * g.nz, 1);
            markDirty();
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Adds a saturation table (copy of table 1, krg = 0 below Sg = 0.25, rescaled), sets SATNUM = 1 and IMBNUM = new table for all cells.\nEdit the table in the SCAL section (region tables) or in the deck.");
    }
    if (!p.hyst.on) ui::hint("Hysteresis is off: the drainage curve is used for imbibition too (only gas below the critical saturation is trapped).");

    // ---- preview: drainage curve, scanning curves, Land trapping curve
    comp::CompCase cc;
    cc.sat = gasTable(case_, 0);
    for (int r = 1; r < nsr; ++r) cc.satRegions.push_back(gasTable(case_, r));
    cc.hyst.on = true; cc.hyst.model = p.hyst.model >= 2 ? 1 : 0;
    cc.trap.sgtMax = compCore ? C.trap.sgtMax : -1;
    if (cc.sat.sg.size() < 2) { ui::hint("Gas saturation table (SGOF) needed for the preview."); return; }
    const double smax = cc.sat.sg.back();
    const int N = 120;
    std::vector<double> S(N + 1), kd(N + 1);
    for (int i = 0; i <= N; ++i) { S[size_t(i)] = smax * i / N; kd[size_t(i)] = comp::gasRelPerm(cc, 0, imbTab, S[size_t(i)], 0.0); }
    const float h = std::max(ImGui::GetContentRegionAvail().y, ImGui::GetFontSize() * 16);
    if (ImPlot::BeginSubplots("##hystplots", 1, 2, ImVec2(-1, h))) {
        if (ImPlot::BeginPlot("Gas relative permeability: drainage and scanning curves")) {
            ImPlot::SetupAxes("Sg", "krg", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("Drainage", S.data(), kd.data(), N + 1, lineSpec(IM_COL32(240, 90, 80, 255), 2.5f));
            const ImU32 cols[4] = {IM_COL32(80, 170, 255, 255), IM_COL32(120, 200, 120, 255), IM_COL32(250, 190, 60, 255), IM_COL32(190, 130, 230, 255)};
            for (int k = 0; k < 4; ++k) {
                const double sh = smax * (k + 1) / 4.0;
                std::vector<double> x, y;
                for (int i = 0; i <= N; ++i) { const double s = sh * i / N; x.push_back(s); y.push_back(comp::gasRelPerm(cc, 0, imbTab, s, sh)); }
                char l[48]; std::snprintf(l, sizeof l, "Scanning from Sg,max = %.2f", sh);
                ImPlot::PlotLine(l, x.data(), y.data(), int(x.size()), lineSpec(cols[k], 1.6f));
            }
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Trapped gas saturation (Land)")) {
            ImPlot::SetupAxes("Sg,max (turning point)", "Sgt", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            std::vector<double> x, y, d;
            for (int i = 0; i <= N; ++i) { double t = 0; const double sh = smax * i / N; comp::gasRelPerm(cc, 0, imbTab, sh, sh, &t); x.push_back(sh); y.push_back(t); d.push_back(sh); }
            ImPlot::PlotLine("Sgt (residually trapped)", x.data(), y.data(), int(x.size()), lineSpec(IM_COL32(250, 190, 60, 255), 2.5f));
            ImPlot::PlotLine("Sgt = Sg,max (all trapped)", x.data(), d.data(), int(x.size()), lineSpec(IM_COL32(150, 150, 150, 255), 1.0f));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------------------------------------------------------------------------ all mechanisms
void App::editTrapping(bool co2) {
    CompositionalInput& C = case_.comp;
    co2::TrappingOptions& T = C.trap;
    auto changed = [&]() { C.structured = true; markDirty(); };
    const float fs = ImGui::GetFontSize();
    if (!ImGui::BeginTabBar("trapmech")) return;

    if (ui::tabItem("Residual (hysteresis)")) {
        editHysteresis(true);
        ImGui::EndTabItem();
    }

    if (ui::tabItem(co2 ? "Solubility & diffusion" : "Diffusion")) {
        if (co2) {
            ImGui::TextWrapped("Solubility trapping: CO2 dissolves into the brine (Spycher-Pruess / Duan-Sun phase equilibrium in every cell, always on). "
                               "CO2-saturated brine is denser than the formation brine; with gravity on, it sinks (convective mixing on fine grids). "
                               "Molecular diffusion transports dissolved CO2 and water vapour across concentration gradients (Fick's law, mole-fraction form).");
            if (ui::beginProps("diff", 0.55f)) {
                if (ui::propBool("Molecular diffusion (DIFFUSE)", &T.diffusion)) changed();
                auto coef = [&](const char* l, double* v) {
                    double cm = *v * 1e4;   // cm2/s
                    if (ui::propD(l, &cm, "%.4g", "cm2/s (deck: DIFFCWAT / DIFFCGAS in m2/day or ft2/day)") && cm >= 0) { *v = cm * 1e-4; changed(); }
                };
                if (T.diffusion) {
                    coef("CO2 in brine D (cm2/s)", &T.dCo2Water);
                    coef("H2O in brine D (cm2/s)", &T.dH2oWater);
                    coef("CO2 in the CO2 phase D (cm2/s)", &T.dCo2Gas);
                    coef("H2O in the CO2 phase D (cm2/s)", &T.dH2oGas);
                }
                static const char* tort[] = {"phi S (OPM Flow / E300)", "Millington-Quirk phi^4/3 S^10/3"};
                if (T.diffusion && ui::propCombo("Porous-medium factor", &T.tortuosity, tort, 2, "Written as '-- RESSIM: TORTUOSITY' when not the default")) changed();
                ui::endProps();
            }
            // dissolution capacity at the datum
            const eos::CompFluid& F = C.fluid;
            const double pd = case_.equil.datumPressure * (case_.units == UnitSystem::Field ? 6894.757293168 : 1e5);
            if (pd > 1e5) {
                const co2::Fluid cf = co2::makeFluid(F.resT, F.salinity);
                const double m = co2::co2Molality(cf, pd), rho = co2::brineDensity(cf, pd);
                const double kgPerM3 = m * co2::MW_CO2 * rho / (1.0 + F.salinity * 0.05844);
                ImGui::Text("Dissolution capacity at the datum: %.3f mol/kg H2O  =  %.1f kg CO2 per m3 of brine", m, kgPerM3);
                if (T.diffusion) {
                    const double L = std::sqrt(T.dCo2Water * 365.25 * 86400.0 * 100.0);
                    ImGui::TextDisabled("Diffusion length sqrt(D t) after 100 years: %.2f m (in free brine; times sqrt(phi S) in the rock)", L);
                }
            }
        } else {
            ImGui::TextWrapped("Molecular diffusion of the hydrocarbon components in the oil and gas phases (DIFFUSE + DIFFC). "
                               "Per-component coefficients other than a single value per phase are written as '-- RESSIM: DIFFCOMP' lines.");
            const int nc = C.nc();
            bool on = !C.diffOil.empty() || !C.diffGas.empty();
            if (ImGui::Checkbox("Molecular diffusion (DIFFUSE)", &on)) {
                if (on) { C.diffOil.assign(size_t(nc), 1e-9); C.diffGas.assign(size_t(nc), 1e-7); }
                else { C.diffOil.clear(); C.diffGas.clear(); }
                changed();
            }
            if (on) {
                C.diffOil.resize(size_t(nc), 0.0); C.diffGas.resize(size_t(nc), 0.0);
                if (ui::beginTable("dc", nc + 1, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_ScrollX)) {
                    ui::tableColumn("D (cm2/s)");
                    for (int i = 0; i < nc; ++i) ui::tableColumn(i < int(C.fluid.names.size()) ? C.fluid.names[size_t(i)].c_str() : "?");
                    ImGui::TableHeadersRow();
                    for (int ph = 0; ph < 2; ++ph) {
                        std::vector<double>& D = ph ? C.diffGas : C.diffOil;
                        ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(ph ? "in gas" : "in oil");
                        for (int i = 0; i < nc; ++i) {
                            ImGui::TableNextColumn(); ImGui::PushID(ph * 100 + i); ImGui::SetNextItemWidth(-FLT_MIN);
                            double cm = D[size_t(i)] * 1e4;
                            if (ImGui::InputDouble("##d", &cm, 0, 0, "%.3g", ImGuiInputTextFlags_EnterReturnsTrue) && cm >= 0) { D[size_t(i)] = cm * 1e-4; changed(); }
                            ImGui::PopID();
                        }
                    }
                    ui::endTable();
                }
                static const char* tort[] = {"phi S (E300)", "Millington-Quirk phi^4/3 S^10/3"};
                ImGui::SetNextItemWidth(fs * 16);
                if (ImGui::Combo("Porous-medium factor", &T.tortuosity, tort, 2)) changed();
            }
        }
        // C17 mechanical dispersion (OPM DISPERC): flux alpha |v| xi A/d dx per phase, isotropic alpha_L = alpha_T
        ImGui::SeparatorText("Mechanical dispersion (DISPERC)");
        {
            GridInput& gd = case_.grid;
            const char* lu = case_.units == UnitSystem::Field ? "ft" : "m";
            double lo = 1e300, hi = -1e300;
            for (double x : gd.disperc) { lo = std::min(lo, x); hi = std::max(hi, x); }
            bool on = !gd.disperc.empty();
            if (ImGui::Checkbox("Dispersion (velocity-dependent mixing, J = alpha |v| xi A / d dx)", &on)) {
                if (on) gd.disperc.assign(size_t(gd.size()), case_.units == UnitSystem::Field ? 0.0328 : 0.01);
                else gd.disperc.clear();
                changed();
            }
            if (on) {
                static double uni = 0.01;
                if (lo == hi) uni = lo;
                ImGui::SetNextItemWidth(fs * 10);
                if (ImGui::InputDouble("Dispersivity (all cells)", &uni, 0, 0, "%.4g", ImGuiInputTextFlags_EnterReturnsTrue) && uni >= 0) {
                    gd.disperc.assign(size_t(gd.size()), uni); changed();
                }
                ImGui::SameLine(); ImGui::TextUnformatted(lu);
                if (lo != hi) ImGui::TextDisabled("Per-cell values from the deck / property editor: %.4g - %.4g %s", lo, hi, lu);
                ImGui::TextDisabled("Written as the OPM Flow GRID array DISPERC; ECLIPSE has no equivalent. Compositional / CO2 engine; "
                                    "SPE11A uses 1e-2 m, SPE11B/C 10 m.");
            }
        }
        ImGui::EndTabItem();
    }

    if (co2 && ui::tabItem("Salt precipitation")) {
        ImGui::TextWrapped("Dry-out near CO2 injectors: water evaporates into the dry CO2 stream, the remaining brine concentrates and NaCl precipitates "
                           "once the solubility limit is reached. The solid salt reduces porosity and permeability (and so injectivity); it re-dissolves "
                           "when fresh brine flows back. Interaction: pore-volume loss raises the local pressure build-up, the permeability factor "
                           "acts on the inter-cell transmissibilities and on the well indices.");
        if (ui::beginProps("salt", 0.55f)) {
            if (ui::propBool("Salt precipitation (PRECSALT + VAPWAT)", &T.saltPrecip)) changed();
            if (T.saltPrecip) {
                double kg = T.saltSolubility * 58.44;
                if (ui::propD("NaCl solubility limit (mol / kg H2O)", &T.saltSolubility, "%.3f", "SALTSOL (kg/sm3 = mol/kg x 58.44); halite: 6.15 at 25 C")) { T.saltSolubility = std::max(0.1, T.saltSolubility); changed(); }
                ui::propText("  = SALTSOL (kg/sm3)", (std::to_string(int(std::lround(kg)))).c_str());
                if (ui::propD("Halite density (kg/m3)", &T.rhoSalt, "%.1f")) { T.rhoSalt = std::max(500.0, T.rhoSalt); changed(); }
                int pm = T.permFact.empty() ? 0 : 1;
                static const char* pms[] = {"Verma-Pruess tubes-in-series", "Table (PERMFACT)"};
                if (ui::propCombo("Permeability model", &pm, pms, 2)) {
                    if (pm == 1 && T.permFact.empty()) for (int i = 0; i <= 10; ++i) { const double f = 0.5 + 0.05 * i; T.permFact.push_back({f, co2::permFactor(co2::TrappingOptions{}, f)}); }
                    if (pm == 0) T.permFact.clear();
                    changed();
                }
                if (pm == 0) {
                    if (ui::propSliderD("phi_r (k = 0 at phi/phi0 =)", &T.vpPhiR, 0.0, 0.95, "%.3f", "Verma & Pruess (1988): fraction of the original porosity at which the permeability vanishes")) changed();
                    if (ui::propSliderD("Gamma (tube length fraction)", &T.vpGamma, 0.05, 0.95, "%.3f")) changed();
                }
                ui::endProps();
            } else ui::endProps();
        }
        if (T.saltPrecip) {
            if (!T.permFact.empty()) {
                struct PF { double a, b; };
                std::vector<PF> rows; for (auto& r : T.permFact) rows.push_back({r.first, r.second});
                if (ui::editTable<PF>("pf", rows, {"phi / phi0", "k / k0"}, [](PF& r, int c) -> double* { return c == 0 ? &r.a : &r.b; }, fs * 8)) {
                    T.permFact.clear(); for (auto& r : rows) T.permFact.push_back({r.a, r.b});
                    std::sort(T.permFact.begin(), T.permFact.end());
                    changed();
                }
            }
            std::vector<double> x, y;
            for (int i = 0; i <= 100; ++i) { const double f = 0.3 + 0.7 * i / 100.0; x.push_back(f); y.push_back(co2::permFactor(T, f)); }
            if (ImPlot::BeginPlot("Permeability reduction", ImVec2(-1, std::max(ImGui::GetContentRegionAvail().y, fs * 12)))) {
                ImPlot::SetupAxes("phi / phi0", "k / k0", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1.05, ImPlotCond_Always);
                ImPlot::PlotLine("k / k0", x.data(), y.data(), int(x.size()), lineSpec(IM_COL32(80, 170, 255, 255), 2.5f));
                ImPlot::EndPlot();
            }
        }
        ImGui::EndTabItem();
    }

    if (co2 && ui::tabItem("Mineral trapping")) {
        ImGui::TextWrapped("Kinetic carbonation: dissolved CO2 reacts with reactive silicates to carbonate minerals (e.g. anorthite + CO2 + 2 H2O -> "
                           "kaolinite + calcite). Rate r = A k(T) (1 - Omega) per m3 of pore volume with k(T) = k25 exp(-Ea/R (1/T - 1/298.15)) and the "
                           "saturation ratio Omega = m_CO2,eq / m_CO2. CO2 and water are removed from the brine (implicit sink); the net solid volume "
                           "change alters porosity / permeability (same coupling as salt). The reactive mineral is depleted. "
                           "Rate constants: Palandri & Kharaka (2004), neutral mechanism.");
        if (!T.minerals.empty()) {   // C16: Excel-style
            using MR = co2::MineralReaction;
            ui::RecordTable<MR> t;
            t.text("Reaction", [](MR& m) -> std::string& { return m.name; })
             .num("log10 k25 (mol/m2/s)", [](MR& m) -> double& { return m.logk25; })
             .num("Ea (kJ/mol)", [](MR& m) -> double& { return m.Ea; }, 1e-3)
             .num("Area (m2/m3 PV)", [](MR& m) -> double& { return m.area; })
             .num("CO2 / rxn", [](MR& m) -> double& { return m.nuCO2; })
             .num("H2O / rxn", [](MR& m) -> double& { return m.nuH2O; })
             .num("dV solid (cm3/mol)", [](MR& m) -> double& { return m.dVsolid; }, 1e6)
             .num("Reactant (PV fraction)", [](MR& m) -> double& { return m.reactant0; })
             .num("Vm reactant (cm3/mol)", [](MR& m) -> double& { return m.vmReactant; }, 1e6)
             .num("m_CO2,eq (mol/kg)", [](MR& m) -> double& { return m.mEq; });
            t.onChange = [](MR& m) { m.area = std::max(0.0, m.area); m.reactant0 = std::clamp(m.reactant0, 0.0, 1.0); m.vmReactant = std::max(1e-9, m.vmReactant); };
            if (t.draw("min", T.minerals, std::min(ImGui::GetFrameHeightWithSpacing() * (float(T.minerals.size()) + 1.8f), fs * 14), nullptr, "%.5g", true)) changed();
        }
        static int preset = 0;
        const std::vector<co2::MineralReaction> P = co2::mineralPresets();
        if (!P.empty()) {
            preset = std::clamp(preset, 0, int(P.size()) - 1);
            ImGui::SetNextItemWidth(fs * 18);
            if (ImGui::BeginCombo("##preset", P[size_t(preset)].name.c_str())) {
                for (int i = 0; i < int(P.size()); ++i) if (ImGui::Selectable(P[size_t(i)].name.c_str(), i == preset)) preset = i;
                ImGui::EndCombo();
            }
            ImGui::SameLine();
            if (ui::accentButton("Add reaction")) { T.minerals.push_back(P[size_t(preset)]); changed(); }
        }
        if (!T.minerals.empty()) {   // rate at the reservoir temperature and the time to consume the reactant
            const double Tk = C.fluid.resT > 0 ? C.fluid.resT : 323.15;
            for (const auto& m : T.minerals) {
                const double k = std::pow(10.0, m.logk25) * std::exp(-m.Ea / 8.314462618 * (1.0 / Tk - 1.0 / 298.15));
                const double r = k * m.area;   // mol / (m3 PV s), far from equilibrium
                const double inv = m.reactant0 / std::max(1e-12, m.vmReactant);   // mol / m3 PV
                const double yrs = r > 0 ? inv / r / (365.25 * 86400.0) : INFINITY;
                ImGui::TextDisabled("%s: k(T) = %.3g mol/m2/s, max rate %.3g mol/m3PV/yr (%.3g kg CO2/m3PV/yr), reactant lasts %.3g years",
                                    m.name.c_str(), k, r * 365.25 * 86400.0, r * 365.25 * 86400.0 * m.nuCO2 * co2::MW_CO2, yrs);
            }
        }
        ImGui::EndTabItem();
    }

    if (co2 && ui::tabItem("Structural & inventory")) {
        ImGui::TextWrapped("The CO2 inventory reported every step splits the CO2 in place by mechanism - every molecule is counted once:\n"
                           "  mineral  -  CO2 fixed in carbonates\n"
                           "  dissolved  -  in the brine (solubility trapping)\n"
                           "  residual  -  free gas below the trapped saturation of its scanning curve (Land) / below the critical saturation\n"
                           "  structural  -  remaining free gas held in place below a seal (its cell's gas outflow is negligible)\n"
                           "  mobile  -  the rest of the free gas\n"
                           "Gas in cells with krg = 0 is reported in addition (it overlaps the residual class).");
        if (ui::beginProps("st", 0.55f)) {
            if (ui::propD("Structural tolerance (outflow / gas volume per year)", &T.structTol, "%.3g",
                          "Free gas counts as structurally trapped when the cell's gas outflow over one year is below this fraction of its mobile gas volume")) {
                T.structTol = std::max(0.0, T.structTol); changed();
            }
            ui::endProps();
        }
        ui::header("Interactions between the mechanisms (all coupled in the fully implicit solve)");
        ImGui::TextWrapped(
            "- Hysteresis x dissolution: trapped gas keeps a large brine contact area and dissolves over time; the trapped saturation is "
            "re-evaluated from the historical maximum Sg after every step.\n"
            "- Dissolution x convection: the density of CO2-saturated brine (Garcia) drives downward flow, which brings fresh brine to the plume.\n"
            "- Diffusion x residual trapping: diffusion moves dissolved CO2 out of trapped zones and spreads the dissolved plume.\n"
            "- Dry-out x salt x injectivity: water vaporisation precipitates salt, porosity and permeability fall (Verma-Pruess / PERMFACT), "
            "transmissibilities and well indices are updated at every step.\n"
            "- Minerals x dissolution x rock: carbonation consumes dissolved CO2 and water (shifting the phase equilibrium) and changes the "
            "pore volume.\n"
            "- Salinity x solubility: the Duan-Sun salting-out lowers the CO2 solubility of saltier brine.");
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}
