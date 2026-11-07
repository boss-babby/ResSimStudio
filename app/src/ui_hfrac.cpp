// C15 user interface
//  * editHydraulicFractures(): Grid & rock > "Hydraulic fractures & unconventional" - spreadsheet table of the HFRAC
//    planes (EDFM fractures of the black-oil / compositional grid), the multistage designer (stages along a well,
//    geometry given or sized from the pump schedule with the PKN / KGD propagation model and the pumped proppant),
//    Langmuir adsorbed gas and stress-dependent permeability, and a pre-check (fracture cells of the grid preview).
//  * stabilityWindow(): Tools > Geomechanics - fault & fracture stability: stress state, planes (from the fractures or
//    typed), Mohr diagram with the Coulomb envelope, slip / dilation tendency, critical pore pressure, time evolution
//    driven by a run's pressure history with the poroelastic stress path, a lower-hemisphere slip-tendency stereonet,
//    cells above the Coulomb line, and fracture movement (rate-and-state slip + Barton-Bandis aperture, events / Mw).
//  * openUserGuide(): Help > User guide (USER_GUIDE.html next to the executable, docs/ in a source tree).
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "implot.h"
#include "plotkit.h"
#include "imgui_internal.h"
#include "ressim/hydrofrac.h"
#include "ressim/faultstab.h"
#include "ressim/grid.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include "platform.h"

using namespace ressim;
namespace fs = std::filesystem;

namespace {
// ---- designer state (session)
struct HfDesign {
    int well = 0, stages = 10;
    double spacing = 0, toeOffset = 0, strike = -1;
    bool perpendicular = true, fromPump = false;
    double halfLength = 100, height = 30, width = 0.003, perm = 50000, poro = 0.35, gamma = 0;
    // pump schedule per stage (SI / field-friendly units shown in the UI)
    int model = 0;                    // 0 PKN, 1 KGD, 2 P3D (C20)
    double dsigmaMPa = 5, maxHeightM = 0;   // P3D: stress contrast of the bounding layers, height cap (0 = 20 h)
    std::vector<double> profX, profH, profW; double hWell = 0; bool runaway = false;   // P3D: last height / width profile
    double rateM3min = 8, timeMin = 90, viscCp = 50, eGPa = 30, nu = 0.25, clE5 = 2, kicMPa = 1.0, sigmaMinMPa = 40;
    double heightM = 30;              // pumped fracture height (m)
    double propTonnes = 80, propRho = 2650, packPhi = 0.35, proppedFrac = 0.8, packPermD = 50;
    bool sized = false;
    double Lhyd = 0, wMean = 0, eff = 0, pnet = 0;
    std::string msg;
    // C17 multi-cluster stage with stress shadowing
    int clusters = 5, perfs = 6;
    double clSpacing = 10, perfMm = 10, cd = 0.6, rhoFluid = 1000;
    bool shadow = true;
    struct ClusterRow { double dsMPa = 0; double L = 0, w0 = 0, share = 0, shadowMPa = 0; };
    std::vector<ClusterRow> cl;
    bool clDone = false;
    hf::ClusterResult clRes;
    std::string clMsg;
};
HfDesign g_hd;


const double kPi = 3.14159265358979323846;

// ---- stability state
struct StabState {
    double depth = 2500, svG = 23, sHG = 18, shG = 15, az = 0, ppG = 10.0;   // MPa/km
    double biot = 1.0, nu = 0.25, mu = 0.6, coh = 0.0;
    double dpUser = 0;                // MPa added to the hydrostatic pressure (no run)
    int pSource = 0;                  // 0 hydrostatic + dp, 1 run: field average, 2 run: max cell, 3 run: cell
    int ci = 1, cj = 1, ck = 1;
    bool stressPathOn = true;
    double rsfTauEnd = 0, rsfSnEnd = 0;
    struct PlaneRow { std::string name = "P1"; double strike = 0, dip = 60, radius = 100; };
    std::vector<PlaneRow> planes;
    int sel = 0;
    int snap = -1;
    // rate and state + aperture
    stab::RsfInput rsf;
    stab::RsfResult res;
    bool haveRes = false;
    double rampMPa = 10, rampDays = 365;   // synthetic pressure ramp when there is no run
    std::string msg;
    bool init = false;
    StabState() { rsf.G = 20; rsf.bb.w0 = 0.05; rsf.bb.Vm = 0.4; rsf.bb.Kn0 = 20; }   // UI units: GPa, mm, mm, GPa/m
};
StabState g_st;

double pToPa(const Case& c) { return c.units == UnitSystem::Field ? 6894.757 : 1e5; }
double lenToM(const Case& c) { return c.units == UnitSystem::Field ? 0.3048 : 1.0; }

void plotLine(const char* label, const std::vector<double>& x, const std::vector<double>& y, ImU32 col = 0, float w = 2.0f) {
    if (x.empty() || x.size() != y.size()) return;
    ImPlotSpec sp;
    if (col) sp.LineColor = col4(col);
    sp.LineWeight = w;
    ImPlot::PlotLine(label, x.data(), y.data(), int(x.size()), sp);
}
// C17: one stage of several perforation clusters growing together (PKN), stress shadow + limited entry
bool clusterPanel(double lm) {   // true: the stage half-length was taken from the result
    bool applied = false;
    HfDesign& d = g_hd;
    if (!ImGui::TreeNodeEx("Multi-cluster stage: stress shadowing & limited entry (C17)")) return false;
    ImGui::TextDisabled("Clusters of one stage share the pump rate through one wellbore pressure; every cluster opens against\n"
                        "S_hmin + its stress offset + the Sneddon stress shadow of the others + its perforation friction.");
    if (ui::beginProps("hfcl", 0.5f)) {
        ui::propI("Clusters per stage", &d.clusters, 1, 50);
        ui::propD("Cluster spacing (m)", &d.clSpacing, "%.4g");
        ui::propI("Perforations per cluster", &d.perfs, 1, 100);
        ui::propD("Perforation diameter (mm)", &d.perfMm, "%.4g");
        ui::propD("Discharge coefficient C_d", &d.cd, "%.3f", "0.6 sharp-edged, 0.9 eroded");
        ui::propD("Fluid density (kg/m3)", &d.rhoFluid, "%.5g");
        ui::propBool("Stress shadowing", &d.shadow, "Off: clusters only interact through the shared rate");
        char b[96];
        const double K = hf::perfFrictionCoef(d.rhoFluid, d.perfs, d.perfMm * 1e-3, d.cd), qc = d.rateM3min / 60.0 / std::max(1, d.clusters);
        std::snprintf(b, sizeof b, "%.3g MPa at the even rate %.3g m3/min", K * qc * qc * 1e-6, qc * 60);
        ui::propText("Perforation friction", b);
        ui::endProps();
    }
    d.clusters = std::clamp(d.clusters, 1, 50);
    d.cl.resize(size_t(d.clusters));
    ui::RecordTable<HfDesign::ClusterRow> t;
    t.num("Stress offset (MPa)", [](HfDesign::ClusterRow& r) -> double& { return r.dsMPa; })
     .computed("Half-length (m)", [](const HfDesign::ClusterRow& r) { return r.L; })
     .computed("Width at well (mm)", [](const HfDesign::ClusterRow& r) { return r.w0; })
     .computed("Rate share (%)", [](const HfDesign::ClusterRow& r) { return r.share; })
     .computed("Shadow stress (MPa)", [](const HfDesign::ClusterRow& r) { return r.shadowMPa; })
     .fixedRows();
    t.draw("clusters", d.cl, ImGui::GetFrameHeightWithSpacing() * (float(std::min(d.clusters, 10)) + 1.8f), nullptr, "%.4g", false);
    if (ImGui::Button(ICON_PLAY "  Simulate the stage")) {
        hf::ClusterInput ci;
        ci.base.model = 0; ci.base.E = d.eGPa * 1e9; ci.base.nu = d.nu; ci.base.height = d.heightM; ci.base.mu = d.viscCp * 1e-3;
        ci.base.rate = d.rateM3min / 60.0; ci.base.cl = d.clE5 * 1e-5; ci.base.tEnd = d.timeMin * 60; ci.base.outputs = 30;
        ci.clusters = d.clusters; ci.spacing = std::max(0.1, d.clSpacing); ci.shadow = d.shadow;
        ci.perfFriction = hf::perfFrictionCoef(d.rhoFluid, d.perfs, d.perfMm * 1e-3, d.cd);
        for (const auto& r : d.cl) ci.stressOffset.push_back(r.dsMPa * 1e6);
        // domain: the single-fracture Nordgren length of the whole rate bounds every cluster
        hf::Input one = ci.base;
        const double Lmax = std::max(20.0, 1.3 * hf::nordgrenLength(one, one.tEnd));
        ci.base.maxLength = Lmax; ci.base.cells = std::clamp(int(Lmax / 1.0), 100, 600);   // UI thread: keep the run to about a second
        d.clDone = hf::simulateClusters(ci, d.clRes);
        if (d.clDone) {
            double Q = 0; for (double q : d.clRes.qf) Q += q;
            for (size_t k = 0; k < d.cl.size(); ++k) {
                d.cl[k].L = d.clRes.Lf[k]; d.cl[k].w0 = d.clRes.w0f[k] * 1e3; d.cl[k].share = Q > 0 ? 100 * d.clRes.qf[k] / Q : 0;
                d.cl[k].shadowMPa = d.clRes.shadowf[k] * 1e-6;
            }
            char b[200];
            std::snprintf(b, sizeof b, "wellbore pressure %.2f MPa above S_hmin, efficiency %.2f, %d steps",
                          d.clRes.pw.empty() ? 0.0 : d.clRes.pw.back() * 1e-6, d.clRes.volInjected > 0 ? d.clRes.volFracture / d.clRes.volInjected : 0.0, d.clRes.steps);
            d.clMsg = b;
        } else d.clMsg = "failed: " + d.clRes.error;
    }
    if (!d.clMsg.empty()) { ImGui::SameLine(); ImGui::TextWrapped("%s", d.clMsg.c_str()); }
    if (d.clDone && !d.clRes.t.empty()) {
        const float h = ImGui::GetFontSize() * 13;
        if (ImGui::BeginTable("##clplots", 2, ImGuiTableFlags_SizingStretchSame)) {
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Cluster rate share##clq", ImVec2(-1, h))) {
                ImPlot::SetupAxes("time (min)", "share of the stage rate (%)");
                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                std::vector<double> tm; for (double t : d.clRes.t) if (t > 0) tm.push_back(t / 60);   // log axis: t > 0 only
                const size_t skip = d.clRes.t.size() - tm.size();
                const int M = int(d.clRes.qf.size());
                for (int k = 0; k < M; ++k) {
                    std::vector<double> y;
                    for (size_t o = skip; o < d.clRes.q.size(); ++o) {
                        const auto& q = d.clRes.q[o];
                        double Q = 0; for (double x : q) Q += x;
                        y.push_back(Q > 0 && size_t(k) < q.size() ? 100 * q[size_t(k)] / Q : 0);
                    }
                    char n[32]; std::snprintf(n, sizeof n, "cluster %d", k + 1);
                    plotLine(n, tm, y, seriesColor(k));
                }
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Half-length per cluster##cll", ImVec2(-1, h))) {
                ImPlot::SetupAxes("cluster", "half-length (m)");
                std::vector<double> x, y;
                for (size_t k = 0; k < d.clRes.Lf.size(); ++k) { x.push_back(double(k + 1)); y.push_back(d.clRes.Lf[k]); }
                ImPlot::PlotBars("half-length", x.data(), y.data(), int(x.size()), 0.6);
                ImPlot::EndPlot();
            }
            ImGui::EndTable();
        }
        if (ImGui::Button("Use the mean cluster half-length for the stages")) {
            double s = 0; for (double L : d.clRes.Lf) s += L;
            const double Lm = s / std::max<size_t>(1, d.clRes.Lf.size());
            d.halfLength = std::max(1.0, d.proppedFrac * Lm) / lm;
            applied = true;
        }
    }
    ImGui::TreePop();
    return applied;
}
}  // namespace

// ============================================================================ Grid & rock: hydraulic fractures
void App::editHydraulicFractures() {
    GridInput& g = case_.grid;
    const char* L = case_.units == UnitSystem::Field ? "ft" : "m";
    const char* P = case_.units == UnitSystem::Field ? "psi" : "bar";
    char hdr[96];
    std::snprintf(hdr, sizeof hdr, "Hydraulic fractures & unconventional (%d fractures)###hfsec", int(g.hfracs.size()));
    if (!ImGui::CollapsingHeader(hdr)) return;
    bool ch = false;
    ui::hint("Planar vertical fractures embedded in the grid (EDFM): each plane is cut into fracture cells connected to the host "
             "cells (T = 2 A k_n / <d>), to each other and to every well completed in a host cell. Works for black-oil and "
             "compositional runs, Cartesian and corner-point grids. Coordinates are the grid x / y of the 3D view in the "
             "case length units.");
    // ---- table
    {
        ui::GridAccess ga;
        ga.count = [&] { return int(g.hfracs.size()); };
        ga.cols = 12;
        ga.getText = [&](int r, int c) -> std::string { return c == 0 ? g.hfracs[size_t(r)].name : std::string(); };
        ga.setText = [&](int r, int c, const std::string& v) { if (c == 0 && !v.empty()) g.hfracs[size_t(r)].name = v; };
        ga.numGet = [&](int r, int c, double& v) -> bool {
            const GridInput::HFrac& F = g.hfracs[size_t(r)];
            switch (c) {
                case 1: v = F.x0; return true; case 2: v = F.y0; return true; case 3: v = F.x1; return true; case 4: v = F.y1; return true;
                case 5: v = F.zTop; return true; case 6: v = F.zBot; return true;
                case 7: v = F.width * (case_.units == UnitSystem::Field ? 12.0 : 1000.0); return true;   // in / mm
                case 8: v = F.perm; return true; case 9: v = F.poro; return true; case 10: v = F.gamma; return true;
                case 11: v = F.satnum; return true;
                default: return false;
            }
        };
        ga.numSet = [&](int r, int c, double v) {
            GridInput::HFrac& F = g.hfracs[size_t(r)];
            switch (c) {
                case 1: F.x0 = v; break; case 2: F.y0 = v; break; case 3: F.x1 = v; break; case 4: F.y1 = v; break;
                case 5: F.zTop = v; break; case 6: F.zBot = v; break;
                case 7: if (v > 0) F.width = v / (case_.units == UnitSystem::Field ? 12.0 : 1000.0); break;
                case 8: F.perm = std::max(0.0, v); break; case 9: F.poro = std::clamp(v, 1e-3, 1.0); break;
                case 10: F.gamma = std::max(0.0, v); break; case 11: F.satnum = std::max(0, int(std::lround(v))); break;
            }
        };
        ga.insert = [&](int at, int n) {
            at = std::clamp(at, 0, int(g.hfracs.size()));
            GridInput::HFrac F;
            if (!g.hfracs.empty()) F = g.hfracs[size_t(std::min(at, int(g.hfracs.size()) - 1))];
            else if (preview_) { F.x0 = preview_->bbMin.x; F.x1 = preview_->bbMax.x; F.y0 = F.y1 = 0.5f * (preview_->bbMin.y + preview_->bbMax.y); }
            for (int k = 0; k < n; ++k) { GridInput::HFrac f = F; f.name = "HF" + std::to_string(g.hfracs.size() + k + 1); g.hfracs.insert(g.hfracs.begin() + at + k, f); }
        };
        ga.erase = [&](int r0, int r1) { r1 = std::min(r1, int(g.hfracs.size()) - 1); if (r0 <= r1) g.hfracs.erase(g.hfracs.begin() + r0, g.hfracs.begin() + r1 + 1); };
        char cx0[24], cy0[24], cx1[24], cy1[24], czt[24], czb[24], cw[24], cg[24];
        std::snprintf(cx0, sizeof cx0, "x0 (%s)", L); std::snprintf(cy0, sizeof cy0, "y0 (%s)", L); std::snprintf(cx1, sizeof cx1, "x1 (%s)", L);
        std::snprintf(cy1, sizeof cy1, "y1 (%s)", L); std::snprintf(czt, sizeof czt, "z top (%s)", L); std::snprintf(czb, sizeof czb, "z bottom (%s)", L);
        std::snprintf(cw, sizeof cw, "width (%s)", case_.units == UnitSystem::Field ? "in" : "mm"); std::snprintf(cg, sizeof cg, "gamma (1/%s)", P);
        const float h = std::min(ImGui::GetFontSize() * 16, ImGui::GetFrameHeight() * (float(g.hfracs.size()) + 2.6f) + 4);
        if (ui::gridTable("hfracs", ga, {"Name", cx0, cy0, cx1, cy1, czt, czb, cw, "perm (mD)", "porosity", cg, "SATNUM"}, h, nullptr, "%.5g", true)) ch = true;
        ui::hint("z top / bottom < 0: the whole grid thickness. perm 0: cubic law w^2/12. gamma: k = k0 exp(-gamma (p_i - p)) "
                 "(proppant embedment / closure). SATNUM 0: the host's table (use a straight-line table for open fractures).");
    }
    // ---- C18 in-run growth (LEFM) - one row per fracture, same order as the table above
    if (!g.hfracs.empty() && ImGui::TreeNodeEx("Fracture growth during the run (LEFM)###hfgrow")) {
        ui::hint("K_IC > 0 lets the fracture grow inside the flow run: the trace above is the MAXIMUM extent. At the start only the "
                 "pieces within the initial half-length of the injection point (pieces in the cells of completed wells, else the "
                 "trace centre) are open; closed pieces conduct like the host rock. After every converged time step each wing's tip "
                 "is tested with K_I = (p_frac - sigma_min(z)) sqrt(pi min(L, H/2)) (KGD while short, PKN once L > H/2) and the next "
                 "piece opens when K_I >= K_IC. Black-oil engine; adjoint gradients are disabled for growing fractures.");
        ui::GridAccess ga;
        ga.count = [&] { return int(g.hfracs.size()); };
        ga.cols = 18;
        ga.getText = [&](int r, int c) -> std::string { return c == 0 ? g.hfracs[size_t(r)].name : std::string(); };
        ga.setText = [&](int, int, const std::string&) {};
        ga.readOnly = [](int, int c) { return c == 0; };
        ga.numGet = [&](int r, int c, double& v) -> bool {
            const GridInput::HFrac& F = g.hfracs[size_t(r)];
            switch (c) {
                case 1: v = F.kic; return true; case 2: v = F.sigma; return true; case 3: v = F.sigmaGrad; return true;
                case 4: v = F.initLen; return true; case 5: v = F.closedMult; return true;
                case 6: v = F.subgrid; return true;
                case 7: v = F.feStress; return true;
                case 8: v = F.ddmWidth; return true; case 9: v = F.ddmE; return true; case 10: v = F.ddmNu; return true;
                case 11: v = F.phaseField; return true; case 12: v = F.pfEta; return true; case 13: v = F.pfH; return true;
                case 14: v = F.pfE; return true; case 15: v = F.pfNu; return true;
                case 16: v = F.subLeak; return true;
                case 17: v = F.pfTwoWay; return true;
                default: return false;
            }
        };
        ga.numSet = [&](int r, int c, double v) {
            GridInput::HFrac& F = g.hfracs[size_t(r)];
            switch (c) {
                case 1: F.kic = std::max(0.0, v); break; case 2: F.sigma = v; break; case 3: F.sigmaGrad = v; break;
                case 4: F.initLen = std::max(0.0, v); break; case 5: F.closedMult = std::clamp(v, 1e-12, 1e6); break;
                case 6: F.subgrid = v >= 0.5 ? 1 : 0; break;
                case 7: F.feStress = v >= 0.5 ? 1 : 0; break;
                case 8: F.ddmWidth = v >= 0.5 ? 1 : 0; break;
                case 9: F.ddmE = std::max(0.0, v); break;
                case 10: F.ddmNu = std::clamp(v, 0.0, 0.49); break;
                case 11: F.phaseField = v >= 0.5 ? 1 : 0; if (F.phaseField) F.subgrid = 0; break;
                case 12: F.pfEta = std::clamp(v, 1e-6, 1.0); break;
                case 13: F.pfH = std::max(0.0, v); break;
                case 14: F.pfE = std::max(0.0, v); break;
                case 15: F.pfNu = std::clamp(v, 0.0, 0.49); break;
                case 16: F.subLeak = v >= 0.5 ? 1 : 0; break;
                case 17: F.pfTwoWay = v >= 0.5 ? 1 : 0; break;
            }
        };
        char cs[32], cgr[40], cl[32];
        std::snprintf(cs, sizeof cs, "sigma_min (%s)", P); std::snprintf(cgr, sizeof cgr, "d sigma / dz (%s/%s)", P, L);
        std::snprintf(cl, sizeof cl, "initial half-length (%s)", L);
        const float h = std::min(ImGui::GetFontSize() * 12, ImGui::GetFrameHeight() * (float(g.hfracs.size()) + 2.6f) + 4);
        if (ui::gridTable("hfgrow", ga, {"Fracture", "K_IC (MPa m^0.5)", cs, cgr, cl, "closed multiplier", "sub-grid tip (0/1)", "FE stress (0/1)", "DDM aperture (0/1)", "DDM E (GPa, 0 = 20)", "DDM nu", "phase field (0/1)", "PF efficiency eta", "PF element h (m, 0 = auto)", "PF E (GPa, 0 = DDM/GEOMECH/20)", "PF nu", "sub-grid leak-off (0/1)", "PF two-way (0/1)"}, h, nullptr, "%.5g", false)) ch = true;
        ui::hint("K_IC 0 = fixed fracture. sigma_min is taken at the fracture mid-depth, the gradient adds (z - z_mid). Typical K_IC: "
                 "shale 0.5 - 2, sandstone 0.5 - 1.5, granite 1 - 3 MPa m^0.5. Results: the 'HFOPEN' cell array (1 = open piece) "
                 "and the run log lines 'fracture X wing +/- grows to ...'. Sub-grid tip = 1 (research option, C21): the wing length is "
                 "continuous - the tip piece opens by the fraction it covers (its leak-off area) and the length is solved inside each "
                 "time step so that K_I = K_IC (HFOPEN shows the open fraction 0..1). 10-30x slower and, on coarse grids, not more "
                 "accurate: the default reaches a length too early, the sub-grid tip too late (leak-off into fresh host cells), both "
                 "converge to the same answer with refinement - run both to bracket the grid-converged growth time. FE stress = 1 (C25, FESTRESS): the closure stress is the smallest principal "
                 "total stress of the finite-element geomechanics at the tip piece's host cell (Grid & rock > Reservoir geomechanics must be "
                 "on) - it follows the poro-elastic stress change and the layering of the grid (Cartesian or corner-point) instead of the "
                 "linear sigma_min(z); black-oil and compositional (CO2) runs. DDM aperture = 1 (C26, DDMWIDTH): after every report step the open "
                 "pieces are opened by a 3-D displacement-discontinuity solve (half space, all fractures interact - stress shadow) under the net "
                 "pressure p - closure stress; the conductance of each piece follows the cubic law (w / w_input)^3 (floor 0.1 w_input). Results: "
                 "Faults & seismicity -> Fracture apertures (DDM). Black-oil engine. Phase field = 1 (C28, PHASEFIELD): the wing length comes "
                 "from a persistent 2-D plane-strain phase-field model of the fracture plane (AT2, volume control, irreversible damage) that "
                 "is advanced at the end of every report step to the crack volume V = eta x (water injected by the wells in the fracture's "
                 "cells) / fracture height; the pieces covered by the phase-field crack open. eta = fluid efficiency (fraction of the "
                 "injected volume stored in the crack, 0.05-0.6 for treatments, much smaller for slow injection). One-way volume coupling: "
                 "the phase-field pressure is reported next to the flow-model fracture pressure (Results below). Not saved in restarts. "
                 "Sub-grid leak-off = 1 (C28 research H2b, with sub-grid tip = 1; deck SUBGRID LEAK): the exposed part of the tip host "
                 "cell is its own storage, so a fresh host cell no longer over-leaks - growth time grid-converged at 10 m cells "
                 "(t60 25.5 / 26.3 / 26.0 d at 10 / 5 / 2.5 m vs 33.6 / 28.8 / 27.0 d without). PF two-way = 1 (C29, deck PHASEFIELD ... TWOWAY): "
                 "after the first phase-field step the crack volume follows the flow model, V = compliance x (flow-model fracture pressure - "
                 "closure), capped by the fluid injected in the step (sub-stepped), and the conductance of every piece follows the cubic law "
                 "of the Sneddon aperture of the phase-field crack - eta is then only the start value. Results: the 'note' column says two-way.");
        // C30 field output of the fractures (HFFIELD)
        ImGui::SeparatorText("Stress / opening field output (3-D maps)");
        {
            Case::HfFieldInput& Hf = case_.hfField;
            static const char* modes[] = {"Auto (growing / DDM fractures)", "On (also phase-field / fixed with sigma_min)", "Off"};
            int m = Hf.on < 0 ? 0 : Hf.on == 1 ? 1 : 2;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 18);
            if (ImGui::Combo("Field output##hff", &m, modes, 3)) { Hf.on = m == 0 ? -1 : m == 1 ? 1 : 0; ch = true; }
            if (ui::beginProps("hffield")) {
                ch |= ui::propD("E (GPa, 0 = DDM / GEOMECH / 20)", &Hf.E, "%.4g");
                ch |= ui::propD("Poisson's ratio", &Hf.nu, "%.3g");
                ch |= ui::propD("Region radius (x fracture extent)", &Hf.radius, "%.3g", "matrix cells within radius x the largest fracture half-length get the stress / displacement maps");
                ch |= ui::propD("Budget (million evaluations / report)", &Hf.budget, "%.3g", "caps cells x open fracture triangles per report; the nearest cells are kept");
                ui::endProps();
            }
            Hf.nu = std::clamp(Hf.nu, 0.0, 0.49); Hf.radius = std::clamp(Hf.radius, 0.05, 20.0); Hf.budget = std::clamp(Hf.budget, 0.01, 1000.0); Hf.E = std::max(0.0, Hf.E);
            ui::hint("After every report step the open pieces are opened by a 3-D displacement-discontinuity solve under the net pressure and "
                     "the stress / displacement / strain change of the surrounding matrix is evaluated (display only, no feedback). 3-D "
                     "arrays: HFWIDTH (mm), HFPNET, HFCOND on the fracture and host cells; HFDSXX/YY/ZZ/XY, HFDSNN (stress shadow on planes "
                     "parallel to the fracture), HFDSM, HFDTAU, HFDISPX/Y/Z, HFDISP (mm), HFSTRV (microstrain) in the matrix. Phase-field "
                     "fractures write PFDAMAGE, PFDSXX/YY/XY, PFDSNN, PFDISP and PFWIDTH from their own model. Results -> Fields & sections "
                     "shows the fracture face, the opening profiles of every step, plan / vertical sections and the phase-field mesh. "
                     "Deck: '-- RESSIM: HFFIELD ON|OFF [E nu] [radius] [budget]'.");
        }
        ImGui::TreePop();
    }
    // ---- designer
    if (ImGui::TreeNodeEx("Multistage designer", ImGuiTreeNodeFlags_DefaultOpen)) {
        std::vector<std::string> wells;
        for (const WellInput& w : case_.wells) wells.push_back(w.name);
        if (wells.empty()) ui::hint("Add a (horizontal) well with completions first - the stages are placed along its completed length.");
        else {
            g_hd.well = std::clamp(g_hd.well, 0, int(wells.size()) - 1);
            if (ui::beginProps("hfd", 0.5f)) {
                ui::propCombo("Well", &g_hd.well, wells);
                ui::propI("Stages", &g_hd.stages, 1, 200);
                char lab[64];
                std::snprintf(lab, sizeof lab, "Stage spacing (%s, 0 = even)", L); ui::propD(lab, &g_hd.spacing, "%.4g");
                if (g_hd.spacing > 0) { std::snprintf(lab, sizeof lab, "First stage from the toe (%s)", L); ui::propD(lab, &g_hd.toeOffset, "%.4g"); }
                ui::propBool("Perpendicular to the well", &g_hd.perpendicular, "Transverse fractures (open against S_hmin when the lateral is drilled along S_hmin)");
                if (!g_hd.perpendicular) { if (g_hd.strike < 0) g_hd.strike = 0; ui::propD("Fracture strike (deg from +y, clockwise)", &g_hd.strike, "%.1f", "= S_Hmax azimuth"); }
                ui::propBool("Size from the pump schedule (PKN / KGD)", &g_hd.fromPump);
                ui::endProps();
            }
            if (g_hd.fromPump) {
                if (ui::beginProps("hfp", 0.5f)) {
                    static const char* models[] = {"PKN (height-contained, L >> H)", "KGD (plane strain, H >> L)", "P3D (height growth through stress barriers)"};
                    ui::propCombo("Propagation model", &g_hd.model, models, 3,
                                  "P3D: PKN-type lateral flow; every cross-section grows in height until K_I = K_IC at both vertical tips "
                                  "(equilibrium height, symmetric bounding layers with closure stress S_hmin + contrast)");
                    ui::propD(g_hd.model == 2 ? "Pay-zone thickness h (m)" : "Fracture height (m)", &g_hd.heightM, "%.4g");
                    if (g_hd.model == 2) {
                        ui::propD("Stress contrast of the bounding layers (MPa)", &g_hd.dsigmaMPa, "%.4g", "S_hmin(barrier) - S_hmin(pay). Height runs away when p_net approaches the contrast");
                        ui::propD("Height cap (m, 0 = 20 h)", &g_hd.maxHeightM, "%.4g");
                    }
                    ui::propD("Pump rate per stage (m3/min)", &g_hd.rateM3min, "%.4g");
                    ui::propD("Pumping time (min)", &g_hd.timeMin, "%.4g");
                    ui::propD("Fluid viscosity (cP)", &g_hd.viscCp, "%.4g");
                    ui::propD("Young's modulus (GPa)", &g_hd.eGPa, "%.4g");
                    ui::propD("Poisson's ratio", &g_hd.nu, "%.3f");
                    ui::propD("Carter leak-off C_L (1e-5 m/s^0.5)", &g_hd.clE5, "%.4g");
                    ui::propD("Fracture toughness K_IC (MPa m^0.5)", &g_hd.kicMPa, "%.4g");
                    ui::propD("Closure stress S_hmin (MPa)", &g_hd.sigmaMinMPa, "%.4g");
                    ui::propD("Proppant per stage (t)", &g_hd.propTonnes, "%.4g");
                    ui::propD("Proppant grain density (kg/m3)", &g_hd.propRho, "%.5g");
                    ui::propD("Pack porosity", &g_hd.packPhi, "%.3f");
                    ui::propD("Propped / hydraulic length", &g_hd.proppedFrac, "%.3f");
                    ui::propD("Pack permeability at closure (D)", &g_hd.packPermD, "%.4g");
                    ui::endProps();
                }
                if (ImGui::Button(ICON_PLAY "  Size the fracture")) {
                    hf::Input in;
                    in.model = g_hd.model; in.E = g_hd.eGPa * 1e9; in.nu = g_hd.nu; in.height = g_hd.heightM; in.mu = g_hd.viscCp * 1e-3;
                    in.rate = g_hd.rateM3min / 60.0; in.cl = g_hd.clE5 * 1e-5; in.kic = g_hd.kicMPa * 1e6; in.sigmaMin = g_hd.sigmaMinMPa * 1e6;
                    in.tEnd = g_hd.timeMin * 60; in.cells = 0;
                    in.stressContrast = g_hd.dsigmaMPa * 1e6; in.maxHeight = g_hd.maxHeightM;
                    hf::Result r;
                    if (hf::simulateAdaptive(in, r) && !r.snaps.empty()) {
                        g_hd.sized = true; g_hd.Lhyd = r.last().L; g_hd.wMean = r.meanWidth(); g_hd.eff = r.last().eff; g_hd.pnet = r.last().pnet0;
                        const double Hf = g_hd.model == 2 ? r.height : g_hd.heightM;   // P3D: length-averaged height
                        const double Lp = std::max(1.0, g_hd.proppedFrac * g_hd.Lhyd);
                        const double wp = hf::proppedWidth(g_hd.propTonnes * 1000, g_hd.propRho, g_hd.packPhi, Lp, Hf);
                        const double lm = lenToM(case_);
                        g_hd.halfLength = Lp / lm; g_hd.height = Hf / lm; g_hd.width = wp / lm;
                        g_hd.perm = g_hd.packPermD * 1000; g_hd.poro = g_hd.packPhi;
                        g_hd.profX = r.last().x; g_hd.profH = r.last().H; g_hd.profW.clear();
                        for (double w : r.last().w) g_hd.profW.push_back(w * 1e3);
                        g_hd.hWell = r.heightWell; g_hd.runaway = r.heightRunaway;
                        char b[320];
                        std::snprintf(b, sizeof b, "hydraulic half-length %.1f m, mean width %.2f mm, efficiency %.2f, net pressure %.2f MPa -> propped %.1f m x %.2f mm",
                                      g_hd.Lhyd, g_hd.wMean * 1e3, g_hd.eff, g_hd.pnet / 1e6, Lp, wp * 1e3);
                        g_hd.msg = b;
                        if (g_hd.model == 2) {
                            std::snprintf(b, sizeof b, "; height %.1f m at the well, %.1f m mean%s", r.heightWell, r.height,
                                          r.heightRunaway ? " - HEIGHT RUNAWAY: the barriers do not contain the fracture (cap reached)" : "");
                            g_hd.msg += b;
                        }
                    } else { g_hd.sized = false; g_hd.msg = "sizing failed: " + r.error; }
                }
                if (!g_hd.msg.empty()) { ImGui::SameLine(); ImGui::TextWrapped("%s", g_hd.msg.c_str()); }
                if (g_hd.model == 2 && g_hd.sized && !g_hd.profX.empty() && g_hd.profH.size() == g_hd.profX.size()) {
                    // P3D side view: fracture outline (+-H/2 about the pay-zone centre) and the pay zone
                    if (ImPlot::BeginPlot("P3D fracture profile##p3d", ImVec2(-1, ImGui::GetFontSize() * 12))) {
                        ImPlot::SetupAxes("distance from the well (m)", "height (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        ImPlot::SetupAxis(ImAxis_Y2, "width (mm)", ImPlotAxisFlags_AuxDefault | ImPlotAxisFlags_AutoFit);
                        std::vector<double> up, dn, pu, pd;
                        for (double H : g_hd.profH) { up.push_back(0.5 * H); dn.push_back(-0.5 * H); pu.push_back(0.5 * g_hd.heightM); pd.push_back(-0.5 * g_hd.heightM); }
                        ImPlotSpec sf; sf.FillColor = col4(Pal::Accent); sf.FillAlpha = 0.25f;
                        ImPlot::PlotShaded("fracture", g_hd.profX.data(), dn.data(), up.data(), int(up.size()), sf);
                        plotLine("top tip", g_hd.profX, up, Pal::Accent);
                        plotLine("bottom tip", g_hd.profX, dn, Pal::Accent);
                        plotLine("pay zone", g_hd.profX, pu, Pal::Muted, 1.5f);
                        plotLine("##pay2", g_hd.profX, pd, Pal::Muted, 1.5f);
                        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                        plotLine("width", g_hd.profX, g_hd.profW, Pal::Warn);
                        ImPlot::EndPlot();
                    }
                }
                if (clusterPanel(lenToM(case_))) notify("Stage half-length set from the multi-cluster simulation (propped fraction applied)");
            }
            if (ui::beginProps("hfg", 0.5f)) {
                char lab[64];
                std::snprintf(lab, sizeof lab, "Half-length (%s)", L); ui::propD(lab, &g_hd.halfLength, "%.4g");
                std::snprintf(lab, sizeof lab, "Height (%s)", L); ui::propD(lab, &g_hd.height, "%.4g");
                double wd = g_hd.width * (case_.units == UnitSystem::Field ? 12.0 : 1000.0);
                std::snprintf(lab, sizeof lab, "Propped width (%s)", case_.units == UnitSystem::Field ? "in" : "mm");
                if (ui::propD(lab, &wd, "%.4g") && wd > 0) g_hd.width = wd / (case_.units == UnitSystem::Field ? 12.0 : 1000.0);
                ui::propD("Fracture permeability (mD, 0 = cubic law)", &g_hd.perm, "%.5g");
                ui::propD("Fracture porosity", &g_hd.poro, "%.3f");
                std::snprintf(lab, sizeof lab, "Conductivity decline gamma (1/%s)", P); ui::propD(lab, &g_hd.gamma, "%.4g");
                char b[64]; std::snprintf(b, sizeof b, "%.4g mD.%s", g_hd.perm * g_hd.width, L); ui::propText("Conductivity k_f w", b);
                ui::endProps();
            }
            auto make = [&](bool append) {
                try {
                    hf::StageLayout S;
                    S.well = wells[size_t(g_hd.well)]; S.stages = g_hd.stages; S.spacing = g_hd.spacing; S.toeOffset = g_hd.toeOffset;
                    S.halfLength = g_hd.halfLength; S.height = g_hd.height; S.width = g_hd.width; S.perm = g_hd.perm; S.poro = g_hd.poro;
                    S.gamma = g_hd.gamma; S.strikeDeg = g_hd.perpendicular ? -1 : g_hd.strike; S.prefix = S.well + "_S";
                    Case c0 = case_;
                    auto F = hf::layoutStages(c0, S);
                    if (!append) g.hfracs.clear();
                    g.hfracs.insert(g.hfracs.end(), F.begin(), F.end());
                    ch = true;
                    notify(std::to_string(F.size()) + " hydraulic fractures placed along " + S.well);
                } catch (const std::exception& e) { g_hd.msg = e.what(); notify(std::string("Stage design: ") + e.what()); }
            };
            if (ui::accentButton("Replace the fractures with this design")) make(false);
            ImGui::SameLine();
            if (ImGui::Button("Append")) make(true);
        }
        ImGui::TreePop();
    }
    // ---- unconventional rock
    if (ImGui::TreeNodeEx("Adsorbed gas & stress-dependent permeability", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginProps("unc", 0.5f)) {
            const bool fld = case_.units == UnitSystem::Field;
            char lab[64];
            std::snprintf(lab, sizeof lab, "Langmuir volume V_L (%s)", fld ? "scf/ton" : "sm3/t");
            ch |= ui::propD(lab, &g.langVL, "%.4g", "0 = no adsorbed gas. Typical shale 1-4 sm3/t (30-130 scf/ton), coal 10-25 sm3/t");
            std::snprintf(lab, sizeof lab, "Langmuir pressure P_L (%s)", P);
            ch |= ui::propD(lab, &g.langPL, "%.4g");
            ch |= ui::propD("Grain density (g/cm3)", &g.langRho, "%.4g");
            std::snprintf(lab, sizeof lab, "Matrix permeability decline gamma (1/%s)", P);
            ch |= ui::propD(lab, &g.permGamma, "%.4g", "k = k0 exp(-gamma (p_i - p)) on every matrix connection");
            ui::endProps();
        }
        if (g.langVL > 0 && g.langPL > 0) {
            // isotherm preview
            std::vector<double> pp, vv;
            for (int k = 0; k <= 60; ++k) { const double p = 6 * g.langPL * k / 60.0; pp.push_back(p); vv.push_back(g.langVL * p / (g.langPL + p)); }
            if (ImPlot::BeginPlot("Langmuir isotherm##iso", ImVec2(-1, ImGui::GetFontSize() * 11))) {
                ImPlot::SetupAxes(P, case_.units == UnitSystem::Field ? "scf/ton" : "sm3/t");
                plotLine("adsorbed", pp, vv, Pal::Gas);
                ImPlot::EndPlot();
            }
        }
        ImGui::TreePop();
    }
    if (preview_ && !g.hfracs.empty()) {
        char b[160];
        std::snprintf(b, sizeof b, "Grid preview: %d fracture cells from %d planes", preview_->numHfCells(), int(g.hfracs.size()));
        ImGui::TextColored(col4(preview_->numHfCells() ? Pal::Ok : Pal::Warn), "%s", b);
        if (!preview_->numHfCells()) ui::hint("No fracture cell yet: rebuild the grid (auto-apply) or check that the traces cross active cells.");
    }
    if (!g.hfracs.empty()) ui::hint("OPM Flow / ECLIPSE ignore HFRAC (comment lines in the exported deck) - run fractured cases with the ResSim engine.");
    if (ch) markDirty(true);
}

// ============================================================================ Tools: geomechanical stability
void App::stabilityWindow() {
    if (!showStability_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 70, ImGui::GetFontSize() * 46), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
    if (stabDock_) {   // opened from the menu: a tab of the central dock node (next to the 3D viewer) in the Model workspace
        stabDock_ = false;
        if (workspace_ == Workspace::Model)
            if (ImGuiDockNode* cn = ImGui::DockBuilderGetCentralNode(ImGui::GetID("MainDock"))) ImGui::SetNextWindowDockID(cn->ID, ImGuiCond_Always);
        ImGui::SetNextWindowFocus();
    }
    if (!ImGui::Begin(ICON_ACTIVITY "  Geomechanics: fault & fracture stability###stabwin", &showStability_)) { ImGui::End(); return; }
    StabState& S = g_st;
    if (!S.init) {
        S.init = true;
        S.planes.clear();
        const double dipOpt = 90 - 0.5 * std::atan(1 / S.mu) * 180 / kPi;
        S.planes.push_back({"Optimal normal fault", S.az, dipOpt, 200});
        auto addTrace = [&](const std::string& n, double x0, double y0, double x1, double y1) {
            double st = std::atan2(x1 - x0, y1 - y0) * 180 / kPi; if (st < 0) st += 360;
            S.planes.push_back({n, st, 90, 100});
        };
        for (const auto& F : case_.grid.hfracs) { if (S.planes.size() > 12) break; addTrace(F.name, F.x0, F.y0, F.x1, F.y1); }
        for (const auto& F : case_.thermal.fractures) { if (S.planes.size() > 12) break; addTrace(F.name, F.x0, F.y0, F.x1, F.y1); }
        if (case_.thermal.enabled && case_.thermal.fracStress.sHmaxGrad > 0) { S.sHG = case_.thermal.fracStress.sHmaxGrad / 1e3; S.shG = case_.thermal.fracStress.sHminGrad / 1e3; S.az = case_.thermal.fracStress.azimuthDeg; }
        if (preview_) S.depth = 0.5 * (preview_->bbMin.z + preview_->bbMax.z) * lenToM(case_);
    }
    // ---- run data (pressure history)
    std::vector<double> tRun, pRun;      // days, Pa
    std::vector<float> lastP; std::vector<double> cellDepth;
    double pPa0 = 0; bool haveRun = false;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        if (R.valid && R.grid && !R.snaps.empty() && R.engine.rfind("ResSim", 0) == 0) {
            haveRun = true;
            const double cv = pToPa(R.input);
            const Grid& G = *R.grid;
            int qSel = -1;
            if (S.pSource == 3) {
                const int i = S.ci - 1, j = S.cj - 1, k = S.ck - 1;
                if (i >= 0 && i < G.nx && j >= 0 && j < G.ny && k >= 0 && k < G.nz) qSel = G.activeOf[size_t(i + G.nx * (j + G.ny * k))];
            }
            for (const Snapshot& sn : R.snaps) {
                double v = 0;
                if (S.pSource == 2) { for (float x : sn.pressure) v = std::max(v, double(x)); }
                else if (S.pSource == 3 && qSel >= 0 && size_t(qSel) < sn.pressure.size()) v = sn.pressure[size_t(qSel)];
                else { double pv = 0; for (size_t a = 0; a < sn.pressure.size() && a < G.poreVolRef.size(); ++a) { v += sn.pressure[a] * G.poreVolRef[a]; pv += G.poreVolRef[a]; } v = pv > 0 ? v / pv : 0; }
                tRun.push_back(sn.time); pRun.push_back(v * cv);
            }
            const int sI = S.snap < 0 ? int(R.snaps.size()) - 1 : std::min(S.snap, int(R.snaps.size()) - 1);
            lastP = R.snaps[size_t(sI)].pressure;
            for (float& x : lastP) x = float(x * cv);
            cellDepth.assign(G.depth.begin(), G.depth.end());
            for (double& d : cellDepth) d *= lenToM(R.input);
            pPa0 = pRun.front();
        }
    }
    const bool useRun = haveRun && S.pSource > 0;
    if (ImGui::BeginTable("stablay", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("in", ImGuiTableColumnFlags_WidthStretch, 0.36f);
        ImGui::TableSetupColumn("out", ImGuiTableColumnFlags_WidthStretch, 0.64f);
        ImGui::TableNextColumn();
        ImGui::BeginChild("stabin", ImVec2(0, 0));
        ui::header("Stress state");
        if (ui::beginProps("ss", 0.55f)) {
            ui::propD("Depth (m)", &S.depth, "%.5g");
            ui::propD("S_v gradient (MPa/km)", &S.svG, "%.4g");
            ui::propD("S_Hmax gradient (MPa/km)", &S.sHG, "%.4g");
            ui::propD("S_hmin gradient (MPa/km)", &S.shG, "%.4g");
            ui::propD("S_Hmax azimuth (deg from N)", &S.az, "%.1f");
            ui::propD("Pore pressure gradient (MPa/km)", &S.ppG, "%.4g");
            ui::propD("Biot coefficient", &S.biot, "%.3f");
            ui::propD("Poisson's ratio (stress path)", &S.nu, "%.3f");
            ui::propD("Friction coefficient", &S.mu, "%.3f");
            ui::propD("Cohesion (MPa)", &S.coh, "%.4g");
            ui::endProps();
        }
        const char* regime = (S.svG >= S.sHG && S.sHG >= S.shG) ? "normal faulting (S_v > S_H > S_h)" : (S.sHG >= S.svG && S.svG >= S.shG) ? "strike-slip (S_H > S_v > S_h)" : (S.sHG >= S.shG && S.shG >= S.svG) ? "reverse faulting (S_H > S_h > S_v)" : "check the stress order";
        ImGui::TextDisabled("Regime: %s", regime);
        ui::header("Pore pressure");
        {
            static const char* src[] = {"Hydrostatic + change", "Run: field average", "Run: maximum cell", "Run: cell (I, J, K)"};
            if (ui::beginProps("pp", 0.55f)) {
                ui::propCombo("Source", &S.pSource, src, 4);
                if (S.pSource == 0 || !haveRun) ui::propD("Pressure change (MPa)", &S.dpUser, "%.4g");
                if (S.pSource == 3) { ui::propI("I", &S.ci, 1, 100000); ui::propI("J", &S.cj, 1, 100000); ui::propI("K", &S.ck, 1, 100000); }
                ui::propBool("Poroelastic stress path", &S.stressPathOn, "d S_h = d S_H = alpha (1 - 2 nu) / (1 - nu) dp inside a laterally extensive reservoir");
                ui::endProps();
            }
            if (S.pSource > 0 && !haveRun) ui::hint("No ResSim run with snapshots yet - using the hydrostatic pressure.");
        }
        ui::header("Planes");
        {
            ui::GridAccess ga;
            ga.count = [&] { return int(S.planes.size()); };
            ga.cols = 4;
            ga.getText = [&](int r, int c) { return c == 0 ? S.planes[size_t(r)].name : std::string(); };
            ga.setText = [&](int r, int c, const std::string& v) { if (c == 0) S.planes[size_t(r)].name = v; };
            ga.numGet = [&](int r, int c, double& v) -> bool { const auto& p = S.planes[size_t(r)]; if (c == 1) v = p.strike; else if (c == 2) v = p.dip; else if (c == 3) v = p.radius; else return false; return true; };
            ga.numSet = [&](int r, int c, double v) { auto& p = S.planes[size_t(r)]; if (c == 1) p.strike = v; else if (c == 2) p.dip = std::clamp(v, 0.0, 90.0); else if (c == 3) p.radius = std::max(1.0, v); };
            ga.insert = [&](int at, int n) { at = std::clamp(at, 0, int(S.planes.size())); for (int k = 0; k < n; ++k) S.planes.insert(S.planes.begin() + at, StabState::PlaneRow{"P" + std::to_string(S.planes.size() + 1), 0, 60, 100}); };
            ga.erase = [&](int r0, int r1) { r1 = std::min(r1, int(S.planes.size()) - 1); if (r0 <= r1) S.planes.erase(S.planes.begin() + r0, S.planes.begin() + r1 + 1); };
            ui::gridTable("stabplanes", ga, {"Name", "Strike (deg)", "Dip (deg)", "Patch radius (m)"}, ImGui::GetFontSize() * 11, &S.sel, "%.4g", true);
            if (ImGui::SmallButton("Reset from the fractures")) S.init = false;
        }
        ImGui::EndChild();
        // ---- outputs
        ImGui::TableNextColumn();
        ImGui::BeginChild("stabout", ImVec2(0, 0));
        const double pHyd = S.ppG * 1e3 * S.depth;
        const double p0 = useRun ? pPa0 : pHyd;
        const double pNow = useRun ? (S.snap < 0 || S.snap >= int(pRun.size()) ? pRun.back() : pRun[size_t(S.snap)]) : pHyd + S.dpUser * 1e6;
        auto stateAt = [&](double p) {
            stab::StressState s = stab::fromGradients(S.depth, S.svG * 1e3, S.sHG * 1e3, S.shG * 1e3, S.az, p0, S.biot);
            if (S.stressPathOn) stab::stressPath(s, p - p0, S.nu); else s.p = p;
            return s;
        };
        const stab::Strength str{S.mu, S.coh * 1e6};
        const stab::StressState s0 = stateAt(p0), s1 = stateAt(pNow);
        if (useRun) {
            int sn = S.snap < 0 ? int(tRun.size()) - 1 : S.snap;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 18);
            if (ImGui::SliderInt("Report step", &sn, 0, int(tRun.size()) - 1, "%d")) S.snap = sn;
            ImGui::SameLine(); ImGui::TextDisabled("t = %.1f d, p = %.2f MPa (initial %.2f)", tRun[size_t(sn)], pNow / 1e6, p0 / 1e6);
        } else ImGui::TextDisabled("p = %.2f MPa (hydrostatic %.2f)", pNow / 1e6, pHyd / 1e6);
        if (ui::beginTable("stabres", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, ImGui::GetFontSize() * 9))) {
            for (const char* c : {"Plane", "sigma_n' (MPa)", "tau (MPa)", "T_s", "T_s initial", "T_d", "CFF (MPa)", "dP_crit (MPa)", "State"}) ui::tableColumn(c);
            ImGui::TableHeadersRow();
            for (size_t k = 0; k < S.planes.size(); ++k) {
                const stab::Plane pl{S.planes[k].strike, S.planes[k].dip};
                const stab::Resolved a = stab::resolve(s0, pl, str), b = stab::resolve(s1, pl, str);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                if (ImGui::Selectable(S.planes[k].name.c_str(), S.sel == int(k), ImGuiSelectableFlags_SpanAllColumns)) S.sel = int(k);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", b.snEff / 1e6);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", b.tau / 1e6);
                ImGui::TableNextColumn(); ImGui::Text("%.3f", b.slipTend);
                ImGui::TableNextColumn(); ImGui::Text("%.3f", a.slipTend);
                ImGui::TableNextColumn(); ImGui::Text("%.3f", b.dilTend);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", b.cff / 1e6);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", b.dpCrit / 1e6);
                ImGui::TableNextColumn();
                if (b.cff >= 0) ImGui::TextColored(col4(Pal::Err), "slipping");
                else if (b.dpCrit < 2e6) ImGui::TextColored(col4(Pal::Warn), "critical");
                else ImGui::TextColored(col4(Pal::Ok), "stable");
            }
            ui::endTable();
        }
        if (ImGui::BeginTabBar("stabtabs")) {
            const float ph = std::max(ImGui::GetFontSize() * 16, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 1.5f);
            if (ui::tabItem("Mohr diagram")) {
                if (ImPlot::BeginPlot("Mohr diagram (effective stresses)##mohr", ImVec2(-1, ph), ImPlotFlags_Equal)) {
                    ImPlot::SetupAxes("effective normal stress (MPa)", "shear stress (MPa)");
                    auto circles = [&](const stab::StressState& s, ImU32 col, const char* lab) {
                        double pr[3] = {s.Sv, s.SH, s.Sh};
                        std::sort(pr, pr + 3);
                        for (double& v : pr) v = (v - s.biot * s.p) / 1e6;
                        const int pairs[3][2] = {{0, 2}, {0, 1}, {1, 2}};
                        for (int c = 0; c < 3; ++c) {
                            const double a = pr[pairs[c][0]], b = pr[pairs[c][1]], cx = 0.5 * (a + b), r = 0.5 * (b - a);
                            std::vector<double> x, y;
                            for (int k = 0; k <= 90; ++k) { const double t = kPi * k / 90; x.push_back(cx + r * std::cos(t)); y.push_back(r * std::sin(t)); }
                            plotLine(lab, x, y, col, c == 0 ? 2.2f : 1.2f);
                        }
                    };
                    circles(s0, Pal::Muted, "initial");
                    circles(s1, Pal::Accent, "current");
                    double smax = std::max({s0.Sv, s0.SH, s0.Sh, s1.Sv, s1.SH, s1.Sh}) / 1e6;
                    std::vector<double> ex{0, smax}, ey{S.coh, S.coh + S.mu * smax};
                    plotLine("Coulomb envelope", ex, ey, Pal::Err, 2.0f);
                    for (size_t k = 0; k < S.planes.size(); ++k) {
                        const stab::Resolved b = stab::resolve(s1, {S.planes[k].strike, S.planes[k].dip}, str);
                        const double x = b.snEff / 1e6, y = b.tau / 1e6;
                        ImPlotSpec sp; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = int(k) == S.sel ? 7.0f : 4.5f;
                        ImPlot::PlotScatter(S.planes[k].name.c_str(), &x, &y, 1, sp);
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Evolution")) {
                std::vector<double> tt, pp;
                if (useRun) { tt = tRun; pp = pRun; }
                else { for (int k = 0; k <= 50; ++k) { tt.push_back(S.rampDays * k / 50.0); pp.push_back(pHyd + S.dpUser * 1e6 * k / 50.0); } }
                if (ImPlot::BeginPlot("Slip tendency and Coulomb stress change##evo", ImVec2(-1, ph * 0.5f))) {
                    ImPlot::SetupAxes("time (days)", "slip tendency T_s");
                    for (size_t k = 0; k < S.planes.size(); ++k) {
                        std::vector<double> y;
                        for (double p : pp) y.push_back(stab::resolve(stateAt(p), {S.planes[k].strike, S.planes[k].dip}, str).slipTend);
                        plotLine(S.planes[k].name.c_str(), tt, y);
                    }
                    std::vector<double> mx{tt.front(), tt.back()}, my{S.mu, S.mu};
                    plotLine("friction mu", mx, my, Pal::Err, 1.5f);
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Margin to slip##marg", ImVec2(-1, ph * 0.48f))) {
                    ImPlot::SetupAxes("time (days)", "dP_crit (MPa)");
                    for (size_t k = 0; k < S.planes.size(); ++k) {
                        std::vector<double> y;
                        for (double p : pp) y.push_back(stab::resolve(stateAt(p), {S.planes[k].strike, S.planes[k].dip}, str).dpCrit / 1e6);
                        plotLine(S.planes[k].name.c_str(), tt, y);
                    }
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Stereonet")) {
                // lower-hemisphere equal-area (Schmidt) projection of the poles, coloured by slip tendency.
                // C23: drawn as an ImPlot plot in unit-circle coordinates (north up), so it has the export button and
                // exports natively (cells, circle, poles, labels, colour bar) in the publication layout.
                const ImVec2 avail = ImGui::GetContentRegionAvail();
                const float side = std::max(ImGui::GetFontSize() * 14, std::min(avail.x, avail.y - ImGui::GetFontSize() * 2));
                double tsMax = 0;
                const int nA = 72, nD = 18;
                std::vector<double> ts(size_t(nA) * nD);
                for (int i = 0; i < nA; ++i)
                    for (int j = 0; j < nD; ++j) {
                        const double st = (i + 0.5) * 360.0 / nA, dp = (j + 0.5) * 90.0 / nD;
                        ts[size_t(i) * nD + size_t(j)] = std::min(2.0, stab::resolve(s1, {st, dp}, str).slipTend);
                        tsMax = std::max(tsMax, ts[size_t(i) * nD + size_t(j)]);
                    }
                auto polePos = [&](double strike, double dip) {   // pole trend = strike - 90, plunge = 90 - dip; unit radius
                    const double trend = (strike - 90) * kPi / 180, plunge = (90 - dip) * kPi / 180;
                    const double r = std::sqrt(2.0) * std::sin(0.5 * (kPi / 2 - plunge));
                    return ImPlotPoint(r * std::sin(trend), r * std::cos(trend));
                };
                const double vTop = std::max(tsMax, S.mu);
                auto colour = [&](double v) {
                    const float t = float(std::clamp(v / std::max(1e-9, vTop), 0.0, 1.0));
                    return ImGui::ColorConvertFloat4ToU32(ImVec4(0.15f + 0.85f * t, 0.35f + 0.5f * (1 - std::fabs(2 * t - 1)), 0.9f * (1 - t) + 0.1f, 1.0f));
                };
                if (ImPlot::BeginPlot("Stereonet (poles, slip tendency)##stereo", ImVec2(side * 1.25f, side), ImPlotFlags_Equal | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus)) {
                    ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoDecorations | ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoDecorations | ImPlotAxisFlags_NoGridLines);
                    ImPlot::SetupAxesLimits(-1.12, 1.62, -1.12, 1.16, ImPlotCond_Always);
                    ImPlot::PushPlotClipRect();
                    for (int i = 0; i < nA; ++i)
                        for (int j = 0; j < nD; ++j) {
                            const double s0d = i * 360.0 / nA, s1d = (i + 1) * 360.0 / nA, d0 = j * 90.0 / nD, d1 = (j + 1) * 90.0 / nD;
                            const ImPlotPoint q[4] = {polePos(s0d, d0), polePos(s1d, d0), polePos(s1d, d1), polePos(s0d, d1)};
                            const double qx[4] = {q[0].x, q[1].x, q[2].x, q[3].x}, qy[4] = {q[0].y, q[1].y, q[2].y, q[3].y};
                            pk::quad(qx, qy, colour(ts[size_t(i) * nD + size_t(j)]));
                        }
                    {
                        double cx[97], cy[97];
                        for (int k = 0; k <= 96; ++k) { cx[k] = std::sin(2 * kPi * k / 96); cy[k] = std::cos(2 * kPi * k / 96); }
                        pk::polyline(cx, cy, 97, IM_COL32(220, 220, 220, 255), 1.5f);
                    }
                    pk::text(0, 1.0, "N", IM_COL32(230, 230, 230, 255), ImVec2(0, -ImGui::GetFontSize() * 0.75f), 1);
                    std::vector<ImPlotPoint> labelled;
                    for (size_t k = 0; k < S.planes.size(); ++k) {
                        const ImPlotPoint p = polePos(S.planes[k].strike, S.planes[k].dip);
                        pk::marker(p.x, p.y, int(k) == S.sel ? 6.0f : 4.0f, IM_COL32(255, 255, 255, 255), IM_COL32(20, 20, 20, 255), 0, 1.0f);
                        bool dup = false;   // planes with the same orientation (e.g. all stages): label the first one only
                        for (const ImPlotPoint& q : labelled) dup = dup || (std::fabs(q.x - p.x) < 0.04 && std::fabs(q.y - p.y) < 0.04);
                        if (dup) continue;
                        labelled.push_back(p);
                        pk::text(p.x, p.y, S.planes[k].name.c_str(), IM_COL32(255, 255, 255, 255), ImVec2(7, -7), 0);
                    }
                    // colour bar (plot coordinates, so it exports with the net)
                    const double bx0 = 1.25, bx1 = 1.35;
                    for (int k = 0; k < 50; ++k) {
                        const double v = vTop * (k + 0.5) / 50.0;
                        pk::rect(bx0, -1 + 2.0 * k / 50, bx1, -1 + 2.0 * (k + 1) / 50, colour(v));
                    }
                    pk::rect(bx0, -1, bx1, 1, 0, IM_COL32(200, 200, 200, 255), 1.0f);
                    char lb[48];
                    std::snprintf(lb, sizeof lb, "T_s %.2f", vTop);
                    pk::text(bx1, 1, lb, IM_COL32(230, 230, 230, 255), ImVec2(5, 0), 0);
                    pk::text(bx1, -1, "0", IM_COL32(230, 230, 230, 255), ImVec2(5, 0), 0);
                    ImPlot::PopPlotClipRect();
                    ImPlot::EndPlot();
                }
                ImGui::TextDisabled("Poles to planes (lower hemisphere, equal area), coloured by slip tendency at the current pressure.");
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Reservoir cells")) {
                if (lastP.empty()) ui::hint("Run the case with the ResSim engine: every cell's pressure at the selected report step is "
                                            "checked against the Coulomb line of the selected plane (stress at the cell depth).");
                else if (S.sel >= 0 && S.sel < int(S.planes.size())) {
                    const stab::Plane pl{S.planes[size_t(S.sel)].strike, S.planes[size_t(S.sel)].dip};
                    int nSlip = 0, nCrit = 0; double tsMax = 0;
                    std::vector<double> hist(20, 0.0);
                    for (size_t a = 0; a < lastP.size() && a < cellDepth.size(); ++a) {
                        const double z = cellDepth[a];
                        stab::StressState s = stab::fromGradients(z, S.svG * 1e3, S.sHG * 1e3, S.shG * 1e3, S.az, S.ppG * 1e3 * z, S.biot);
                        if (S.stressPathOn) stab::stressPath(s, lastP[a] - S.ppG * 1e3 * z, S.nu); else s.p = lastP[a];
                        const stab::Resolved r = stab::resolve(s, pl, str);
                        if (r.cff >= 0) ++nSlip; else if (r.dpCrit < 2e6) ++nCrit;
                        tsMax = std::max(tsMax, r.slipTend);
                        hist[size_t(std::clamp(int(r.slipTend / std::max(1e-9, 1.5 * S.mu) * 20), 0, 19))] += 1;
                    }
                    ImGui::Text("Plane %s: %d cells above the Coulomb line, %d within 2 MPa, max T_s %.3f (of %d active cells)",
                                S.planes[size_t(S.sel)].name.c_str(), nSlip, nCrit, tsMax, int(lastP.size()));
                    ui::hint("Initial pressure = hydrostatic gradient at the cell depth; the change is the simulated cell pressure minus it.");
                    std::vector<double> xs(20), ys(20);
                    for (int k = 0; k < 20; ++k) { xs[size_t(k)] = (k + 0.5) * 1.5 * S.mu / 20; ys[size_t(k)] = hist[size_t(k)]; }
                    if (ImPlot::BeginPlot("Slip tendency of the cells##hist", ImVec2(-1, ph * 0.8f))) {
                        ImPlot::SetupAxes("slip tendency", "cells");
                        ImPlotSpec sp; sp.FillColor = col4(Pal::Accent);
                        ImPlot::PlotBars("cells", xs.data(), ys.data(), 20, 1.5 * S.mu / 22, sp);
                        std::vector<double> mx{S.mu, S.mu}, my{0, *std::max_element(ys.begin(), ys.end()) + 1};
                        plotLine("friction", mx, my, Pal::Err, 2.0f);
                        ImPlot::EndPlot();
                    }
                }
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Fracture movement")) {
                if (S.sel < 0 || S.sel >= int(S.planes.size())) ui::hint("Select a plane.");
                else {
                    stab::RsfInput& in = S.rsf;
                    ImGui::TextWrapped("Quasi-dynamic slip of %s (patch radius %.0f m) with rate-and-state friction driven by the pore-pressure "
                                       "history; the aperture follows Barton-Bandis closure plus shear dilation.", S.planes[size_t(S.sel)].name.c_str(), S.planes[size_t(S.sel)].radius);
                    if (ui::beginProps("rsf", 0.55f)) {
                        ui::propD("a (direct effect)", &in.a, "%.4f");
                        ui::propD("b (evolution effect)", &in.b, "%.4f", "a < b: velocity weakening (potentially seismic); a > b: stable creep");
                        ui::propD("D_c (m)", &in.Dc, "%.3g");
                        ui::propD("f_0 at V_0", &in.f0, "%.3f");
                        ui::propD("Shear modulus G (GPa)", &in.G, "%.4g");
                        ui::propD("Residual aperture (mm)", &in.bb.w0, "%.4g");
                        ui::propD("Max. closure V_m (mm)", &in.bb.Vm, "%.4g");
                        ui::propD("Normal stiffness K_n0 (GPa/m)", &in.bb.Kn0, "%.4g");
                        ui::propD("Dilation angle (deg)", &in.bb.psiDeg, "%.3g");
                        if (!useRun) { ui::propD("Ramp: pressure increase (MPa)", &S.rampMPa, "%.4g"); ui::propD("Ramp duration (days)", &S.rampDays, "%.4g"); }
                        ui::endProps();
                    }
                    ImGui::TextDisabled("(G, K_n0 in GPa; apertures in mm - converted on run)");
                    if (ui::accentButton(ICON_PLAY "  Simulate slip")) {
                        stab::RsfInput r = in;
                        r.G = in.G * 1e9; r.bb.w0 = in.bb.w0 * 1e-3; r.bb.Vm = in.bb.Vm * 1e-3; r.bb.Kn0 = in.bb.Kn0 * 1e9; r.bb.dilMax = 2e-3;
                        r.radius = S.planes[size_t(S.sel)].radius;
                        const stab::Plane pl{S.planes[size_t(S.sel)].strike, S.planes[size_t(S.sel)].dip};
                        std::vector<double> tt, pp;
                        if (useRun) { tt = tRun; pp = pRun; }
                        else { tt = {0, S.rampDays}; pp = {pHyd, pHyd + S.rampMPa * 1e6}; }
                        const stab::Resolved a0 = stab::resolve(stateAt(pp.front()), pl, str);
                        r.sn0 = a0.snEff; r.tau0 = a0.tau;
                        r.t.clear(); r.dTau.clear(); r.dSn.clear();
                        for (size_t k = 0; k < tt.size(); ++k) {
                            const stab::Resolved ak = stab::resolve(stateAt(pp[k]), pl, str);
                            r.t.push_back(tt[k] * 86400); r.dTau.push_back(ak.tau - a0.tau); r.dSn.push_back(ak.snEff - a0.snEff);
                        }
                        if (r.tau0 >= r.f0 * r.sn0) S.msg = "the plane is already above the steady-state friction f_0 at the start - lower f_0 or check the stresses";
                        else S.msg.clear();
                        S.rsfTauEnd = r.tau0 + (r.dTau.empty() ? 0.0 : r.dTau.back()); S.rsfSnEnd = r.sn0 + (r.dSn.empty() ? 0.0 : r.dSn.back());
                        S.haveRes = stab::simulateRsf(r, S.res);
                        if (!S.res.error.empty()) S.msg = S.res.error;
                    }
                    if (!S.msg.empty()) ImGui::TextColored(col4(Pal::Warn), "%s", S.msg.c_str());
                    if (S.haveRes && !S.res.samples.empty()) {
                        const stab::RsfResult& o = S.res;
                        ImGui::Text("Patch stiffness %.3g Pa/m vs critical %.3g Pa/m (%s); %d steps; %zu seismic event(s)", o.kPatch, o.kCrit,
                                    o.kPatch < o.kCrit ? "unstable: stick-slip possible" : "stable sliding", o.steps, o.events.size());
                        if (o.events.empty() && !S.res.samples.empty()) {
                            const double fEnd = (in.f0 > 0 && S.rsfTauEnd > 0 && S.rsfSnEnd > 0) ? S.rsfTauEnd / S.rsfSnEnd : 0;
                            ImGui::TextWrapped("No seismic event: at the end of the loading tau/sigma_n' = %.3f (f_0 = %.3f)%s.", fEnd, in.f0,
                                               S.stressPathOn ? " - with the poroelastic stress path a pressure increase also shrinks the deviatoric stress; "
                                                                "switch the stress path off (localised injection near the fault) or raise the ramp" : "");
                        }
                        for (const auto& e : o.events)
                            ImGui::BulletText("t %.2f d: slip %.4f m, V_max %.3g m/s, M_0 %.3g N.m, M_w %.2f", e.tStart / 86400, e.slip, e.Vmax, e.M0, e.Mw);
                        std::vector<double> t, sl, v, ap, km;
                        const double w0 = o.samples.front().aperture;
                        for (const auto& s : o.samples) { t.push_back(s.t / 86400); sl.push_back(s.slip * 1e3); v.push_back(std::max(1e-20, s.V)); ap.push_back(s.aperture * 1e3); km.push_back(std::pow(s.aperture / w0, 3)); }
                        if (ImPlot::BeginPlot("Slip and slip rate##rsf1", ImVec2(-1, std::max(ph * 0.5f, ImGui::GetFontSize() * 15)))) {
                            ImPlot::SetupAxes("time (days)", "slip (mm)");
                            ImPlot::SetupAxis(ImAxis_Y2, "slip rate (m/s)", ImPlotAxisFlags_AuxDefault);
                            ImPlot::SetupAxisScale(ImAxis_Y2, ImPlotScale_Log10);
                            plotLine("slip", t, sl, Pal::Accent);
                            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                            plotLine("slip rate", t, v, Pal::Err, 1.5f);
                            ImPlot::EndPlot();
                        }
                        if (ImPlot::BeginPlot("Aperture and permeability##rsf2", ImVec2(-1, std::max(ph * 0.48f, ImGui::GetFontSize() * 15)))) {
                            ImPlot::SetupAxes("time (days)", "aperture (mm)");
                            ImPlot::SetupAxis(ImAxis_Y2, "k / k_initial", ImPlotAxisFlags_AuxDefault);
                            plotLine("aperture", t, ap, Pal::Water);
                            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                            plotLine("permeability multiplier (w/w0)^3", t, km, Pal::Oil, 1.5f);
                            ImPlot::EndPlot();
                        }
                    }
                }
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::EndChild();
        ImGui::EndTable();
    }
    ImGui::End();
}

// ============================================================================ Help: user guide
void App::openUserGuide() {
    const fs::path exe = fs::u8path(plat::exePath());
    const fs::path cands[] = {exe.parent_path() / "USER_GUIDE.html", exe.parent_path() / "docs" / "USER_GUIDE.html",
                              exe.parent_path().parent_path() / "release" / "USER_GUIDE.html", exe.parent_path().parent_path() / "docs" / "USER_GUIDE.html",
                              exe.parent_path().parent_path() / "share" / "ressim" / "USER_GUIDE.html"};
    for (const fs::path& p : cands) {
        std::error_code ec;
        if (fs::exists(p, ec)) { plat::openPath(p.u8string()); return; }
    }
    notify("USER_GUIDE.html was not found next to ResSimStudio.exe");
}
