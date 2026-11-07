// C24: fault mechanics and induced seismicity as part of the main model
//  * editFaultMechanics(): Grid editor section below FAULTS / MULTFLT - in-situ stress field, per-fault Coulomb / Dieterich /
//    Gutenberg-Richter / two-way coupling / rate-and-state parameters (STRESSFIELD / FAULTMECH / FAULTRSF) and the G1-c
//    fracture contact coupling of the hydraulic fractures (GEOCOUPLE). Everything is stored in the deck and runs inside
//    the normal Run.
//  * faultPanels(): every FAULTS face as a panel of the main 3D view, coloured by the chosen mechanical quantity of the
//    current report step (fixed colour range over the run), plus the event catalogue up to the current time.
//  * faultsTab(): Results tab "Faults & seismicity" - Coulomb stress, slip tendency, reactivated area, expected vs
//    realised events, magnitudes, Gutenberg-Richter, fault-face map, catalogue, rate-and-state traces and G1-c history.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "plotkit.h"
#include "implot.h"
#include "ressim/faultmech.h"
#include "ressim/grid.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <map>

using namespace ressim;

namespace {
const char* kFaultProps[] = {"Coulomb failure function CFF", "Coulomb stress change dCFS", "Slip tendency tau / sigma'n", "Effective normal stress sigma'n",
                             "Seismicity rate R / r (log10)", "Expected events per km2", "Transmissibility multiplier",
                             "Cumulative slip (mm, FAULTSLIP)", "Shear stress change from slip (FAULTSLIP)"};
constexpr int kNumFaultProps = 9;

const std::vector<float>* faultField(const faultmech::StepOut& s, int prop) {
    switch (prop) {
        case 0: return &s.cff; case 1: return &s.dcfs; case 2: return &s.st; case 3: return &s.sne;
        case 4: return &s.rate; case 5: return &s.cum; case 7: return &s.slip; case 8: return &s.dtauSlip; default: return &s.mult;
    }
}
double faultValue(const faultmech::Output& O, const faultmech::StepOut& s, size_t i, int prop) {
    const std::vector<float>* v = faultField(s, prop);
    if (i >= v->size()) return NAN;
    double x = (*v)[i];
    if (prop == 4) x = std::log10(std::max(1e-6, x));
    if (prop == 7) x *= 1e3;
    if (prop == 5) { const double a = O.faces[i].area * O.mPerUnit * O.mPerUnit * 1e-6; x = a > 0 ? x / a : 0; }
    return x;
}
bool divergingProp(int prop) { return prop == 0 || prop == 1 || prop == 8; }

void plotXY(const char* label, const std::vector<double>& X, const std::vector<double>& Y, ImU32 col, float w, bool scatter = false) {
    const int n = int(std::min(X.size(), Y.size()));
    if (n <= 0) return;
    if (fig::recorder().recording()) recSeries(label, X.data(), Y.data(), n, col, w, scatter, scatter);
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = w * g_plotLineScale;
    if (scatter) { sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3.0f; sp.MarkerFillColor = col4(col); sp.MarkerLineColor = col4(col); ImPlot::PlotScatter(label, X.data(), Y.data(), n, sp); }
    else ImPlot::PlotLine(label, X.data(), Y.data(), n, sp);
}
void hline(const char* label, double y, ImU32 col) {
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = 1.0f; sp.Flags = ImPlotItemFlags_NoLegend | ImPlotInfLinesFlags_Horizontal;
    ImPlot::PlotInfLines(label, &y, 1, sp);
}
// step of the fault output closest to a report-step time
size_t stepAt(const faultmech::Output& O, double t) {
    size_t b = 0;
    for (size_t i = 0; i < O.steps.size(); ++i) if (std::fabs(O.steps[i].t - t) < std::fabs(O.steps[b].t - t)) b = i;
    return b;
}
}  // namespace

// ============================================================================ inputs
// C25: sequential finite-element geomechanics of black-oil and compositional (CO2) runs on Cartesian and corner-point grids
void App::editReservoirGeomech() {
    if (case_.thermal.enabled) return;   // thermal decks: Thermal > Geomechanics (THM, fixed-stress)
    geomech::Input& G = case_.geoSeq;
    ui::header("Reservoir geomechanics (finite elements, coupled with the run)");
    ui::hint("Linear poro-elastic (optionally elasto-plastic) 8-node hexahedra on the flow grid itself: Cartesian cells directly, "
             "corner-point cells as isoparametric hexahedra on the averaged cell corners (dip, layering and fault throws kept), plus "
             "over / under / side-burden padding. After every report step the simulated pore pressure loads the mesh; the results are "
             "displacement / surface uplift, the full stress tensor, sigma_min, volumetric strain and failure (3D view and Results). "
             "Two-way: pore volume x phi_FE / phi_flow and transmissibility x the stress-dependent permeability factor before the next step. "
             "The FE stress also drives fracture growth (HFGROW FE stress) and fault stresses (stress source FE) - this replaces XFEM: "
             "fractures propagate along grid pieces with an LEFM tip criterion against the FE closure stress. "
             "Black oil and compositional / CO2. Stored as '-- RESSIM: GEOMECH' (ignored by ECLIPSE / OPM).");
    if (G.enabled && G.E.empty()) { G.E = {10e9}; G.nu = {0.25}; G.biot = {1.0}; G.alphaT = {0.0}; }
    if (editGeomechInput(G, false)) markDirty(true);
    if (G.enabled) {
        int nFe = 0; for (const auto& h : case_.grid.hfracs) nFe += h.feStress;
        ImGui::TextDisabled("%s grid | %d of %zu hydraulic fractures use the FE closure stress | fault stress source: %s",
                            case_.grid.cornerPoint() ? "Corner-point (isoparametric)" : "Cartesian", nFe, case_.grid.hfracs.size(),
                            case_.faultMech.on() ? (case_.faultMech.source == 1 ? "FE" : "analytic") : "off");
    }
}

void App::editFaultMechanics() {
    GridInput& gi = case_.grid;
    FaultMechInput& M = case_.faultMech;
    const bool fld = case_.units == UnitSystem::Field;
    const char* pU = fld ? "psi" : "bar";
    ui::header("Fault mechanics & induced seismicity (coupled with the run)");
    bool on = M.svGrad > 0;
    if (ImGui::Checkbox("Compute fault stresses and seismicity during the run", &on)) {
        if (on) {
            M.svGrad = fld ? 1.0 : 0.2262;   // ~ 2300 kg/m3 overburden
            if (M.faults.empty()) for (const FaultInput& F : gi.faults) { FaultMechInput::Row r; r.name = F.name; if (fld) r.coh = 0; M.faults.push_back(r); }
        } else M.svGrad = 0;
        markDirty(true);
    }
    ImGui::SameLine(); ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Every FAULTS face listed below becomes a mechanical object of the model: after each report step the simulated\n"
                          "pore pressure on both sides gives the Coulomb stress, slip tendency and the Dieterich seismicity rate; reactivated\n"
                          "faces can raise their transmissibility (two-way coupling); FAULTRSF adds a quasi-dynamic rate-and-state rupture\n"
                          "model on one layer driven by the simulated pressure history. Results: 3D view (fault panels + events) and the\n"
                          "Results tab 'Faults & seismicity'. Stored as '-- RESSIM:' lines (ignored by ECLIPSE / OPM).");
    if (!on) {
        if (gi.faults.empty()) ui::hint("Define faults first (FAULTS above), then enable the mechanics.");
    } else {
        const bool feOn = case_.thermal.enabled ? case_.thermal.geo.enabled : case_.geoSeq.enabled;
        if (M.source == 1 && !feOn)
            ImGui::TextColored(col4(Pal::Warn), "Stress source FE needs the finite-element geomechanics (%s) - the analytic field is used until then.",
                               case_.thermal.enabled ? "Thermal > Geomechanics" : "Reservoir geomechanics below");
        bool ch = false;
        if (ui::beginProps("fmstress")) {
            char lab[96];
            std::snprintf(lab, sizeof lab, "Vertical stress gradient (%s/%s)", pU, fld ? "ft" : "m");
            ch |= ui::propD(lab, &M.svGrad, "%.4g", "Total overburden stress gradient, Sv = Sv0 + gradient x depth");
            std::snprintf(lab, sizeof lab, "Sv at depth 0 (%s)", pU);
            ch |= ui::propD(lab, &M.sv0, "%.4g");
            ch |= ui::propD("Shmin / Sv", &M.shRatio, "%.3f", "Minimum horizontal total stress ratio");
            ch |= ui::propD("SHmax / Sv", &M.sHRatio, "%.3f", "Maximum horizontal total stress ratio (> 1: strike-slip regime when Shmin < Sv)");
            ch |= ui::propD("SHmax azimuth (deg from +x to +y)", &M.azimuth, "%.1f");
            ch |= ui::propD("Stress path dSh / dp", &M.stressPath, "%.3f", "Poroelastic horizontal total-stress change per unit pore-pressure change, typically (1-2nu)/(1-nu) x alpha ~ 0.5-0.7");
            static const char* kSrc[] = {"Analytic field + stress path", "Finite-element geomechanics (FE)"};
            ch |= ui::propCombo("Stress change source", &M.source, kSrc, 2,
                                "FE: the face stress = in-situ field above + the change of the full total-stress tensor of the finite-element "
                                "geomechanics since the start (black oil / CO2: Reservoir geomechanics below; thermal: the THM mechanics). "
                                "Captures stress arching, poro-elastic stressing of the far side of the fault and thermal stresses, "
                                "which the scalar stress path cannot. Needs the geomechanics to be enabled.");
            if (case_.thermal.enabled) {
                ch |= ui::propD("Thermal stress: E (GPa, 0 = off)", &M.thE, "%.4g", "Analytic field only: dS_h = dS_H = -E alpha / (1 - nu) dT (cooling unloads the faces)");
                ch |= ui::propD("Thermal stress: alpha (1/K)", &M.thAlpha, "%.3g");
                ch |= ui::propD("Thermal stress: nu", &M.thNu, "%.3f");
                if (M.thE > 0 && M.thAlpha > 0)
                    ImGui::TextDisabled("dS_h per K of cooling: %.3g %s", M.thE * 1e9 * M.thAlpha / std::max(1e-6, 1 - M.thNu) * (fld ? 1.450377e-4 : 1e-5), pU);
            }
            ui::endProps();
        }
        std::vector<std::string> names; for (const FaultInput& F : gi.faults) names.push_back(F.name);
        ui::RecordTable<FaultMechInput::Row> t;
        t.choiceStr("Fault", [](FaultMechInput::Row& r) -> std::string& { return r.name; }, names)
         .num("Friction mu", [](FaultMechInput::Row& r) -> double& { return r.mu; })
         .num(fld ? "Cohesion (psi)" : "Cohesion (bar)", [](FaultMechInput::Row& r) -> double& { return r.coh; })
         .optNum(fld ? "A sigma (psi)" : "A sigma (bar)", [](FaultMechInput::Row& r) -> double& { return r.aSig; }, [](FaultMechInput::Row& r) { return r.aSig <= 0; },
                 [](FaultMechInput::Row& r) { r.aSig = 0; }, "0.003 sigma'n")
         .num("t_a (yr)", [](FaultMechInput::Row& r) -> double& { return r.ta; })
         .num("r0 (/km2/yr)", [](FaultMechInput::Row& r) -> double& { return r.r0; })
         .num("b-value", [](FaultMechInput::Row& r) -> double& { return r.bValue; })
         .num("M min", [](FaultMechInput::Row& r) -> double& { return r.mMin; })
         .num("T factor when slipping", [](FaultMechInput::Row& r) -> double& { return r.permFactor; })
         .num("ramp (x sigma'n)", [](FaultMechInput::Row& r) -> double& { return r.permRamp; })
         .integer("RSF layer (0 off)", [](FaultMechInput::Row& r) -> int& { return r.rsfLayer; })
         .num("RSF a", [](FaultMechInput::Row& r) -> double& { return r.rsfA; })
         .num("RSF b", [](FaultMechInput::Row& r) -> double& { return r.rsfB; })
         .num("RSF Dc (mm)", [](FaultMechInput::Row& r) -> double& { return r.rsfDcMm; });
        t.onChange = [&](FaultMechInput::Row& r) {
            r.mu = std::clamp(r.mu, 0.0, 2.0); r.ta = std::max(1e-3, r.ta); r.r0 = std::max(0.0, r.r0); r.bValue = std::clamp(r.bValue, 0.3, 3.0);
            r.permFactor = std::max(1e-3, r.permFactor); r.permRamp = std::clamp(r.permRamp, 1e-4, 10.0); r.rsfLayer = std::clamp(r.rsfLayer, 0, std::max(1, gi.nz)); r.rsfDcMm = std::max(1e-4, r.rsfDcMm);
        };
        t.newRow = [&](int) { FaultMechInput::Row r; for (const FaultInput& F : gi.faults) { bool used = false; for (auto& q : M.faults) used |= q.name == F.name; if (!used) { r.name = F.name; break; } } return r; };
        const float h = std::min(ImGui::GetFontSize() * 14, ImGui::GetFrameHeight() * (float(M.faults.size()) + 2.6f) + 4);
        ch |= t.draw("fmrows", M.faults, h, nullptr, "%.4g");
        if (ImGui::Button("+ all faults")) {
            for (const FaultInput& F : gi.faults) { bool used = false; for (auto& q : M.faults) used |= q.name == F.name; if (!used) { FaultMechInput::Row r; r.name = F.name; M.faults.push_back(r); ch = true; } }
        }
        ImGui::SameLine();
        ui::hint("Two-way coupling: 'T factor when slipping' > 1 multiplies the transmissibility of the cross-fault connections of every face with CFF >= 0 "
                 "(smooth in the excess stress), applied before the next report step. 1 = one-way (monitoring only).");
        // initial state pre-check on the hydrostatic gradient of the equilibration
        static std::string preMsg; static unsigned preVer = ~0u;
        if (preVer != caseVersion_) { preVer = caseVersion_; preMsg.clear(); }
        if (ImGui::Button("Pre-check initial fault state")) {
            try {
                Grid g; g.build(case_);
                faultmech::Monitor mon; std::string err;
                if (!mon.setup(case_, g, err)) preMsg = "Not usable: " + err;
                else {
                    const EquilInput& E = case_.equil;
                    const double grad = fld ? 0.433 : 0.0981;
                    std::vector<float> p(size_t(g.numActive()));
                    for (int a = 0; a < g.numActive(); ++a) p[size_t(a)] = float(E.datumPressure + grad * (g.center[size_t(g.globalOf[size_t(a)])].z - E.datumDepth));
                    mon.step(0, p);
                    const auto& O = mon.out(); const auto& S = O.steps[0];
                    char b[512]; preMsg.clear();
                    for (size_t f = 0; f < O.names.size(); ++f) {
                        double cmax = -1e300, stmax = 0; int nf = 0, crit = 0;
                        for (size_t i = 0; i < O.faces.size(); ++i) if (O.faces[i].fault == int(f)) { ++nf; cmax = std::max(cmax, double(S.cff[i])); stmax = std::max(stmax, double(S.st[i])); crit += S.cff[i] >= 0; }
                        std::snprintf(b, sizeof b, "%s: %d faces, max CFF %.3g %s (closest to failure), max slip tendency %.3f, %d critical at t0; Mw cap %.2f\n",
                                      O.names[f].c_str(), nf, cmax, pU, stmax, crit, f < O.mMax.size() ? O.mMax[f] : 0.0);
                        preMsg += b;
                    }
                    preMsg += "(hydrostatic water gradient from the datum of EQUIL region 1)";
                }
            } catch (const std::exception& e) { preMsg = e.what(); }
        }
        if (!preMsg.empty()) ImGui::TextWrapped("%s", preMsg.c_str());
        // ---- C26 FAULTSLIP: quasi-static slip + elastic stress transfer between the faces (3-D DDM, triangular dislocations)
        ch |= ImGui::Checkbox("Quasi-static fault slip with stress transfer (FAULTSLIP)", &M.slipOn);
        ImGui::SameLine(); ImGui::TextDisabled("(?)");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("After every report step the faces above their Coulomb strength slip until they are back at the strength\n"
                              "(projected Gauss-Seidel on the 3-D displacement-discontinuity matrix of all faces, Nikkhoo & Walter 2015\n"
                              "triangular dislocations, half space). Slip unloads the slipping face and loads its neighbours (stress transfer),\n"
                              "so reactivation can propagate along the fault. Outputs: cumulative slip per face, seismic moment M0 = G sum(A s),\n"
                              "DDM slip events in the catalogue. Friction weakening dMu > 0 makes slip run away (cascade -> seismic event).\n"
                              "Two-way: the transmissibility factor ramps with slip / slip_ref. Stored as '-- RESSIM: FAULTSLIP'.");
        if (M.slipOn && ui::beginProps("fmslip")) {
            ch |= ui::propD("Young's modulus (GPa)", &M.slipE, "%.4g", "Elastic host rock of the DDM (shear modulus G = E / 2(1+nu))");
            ch |= ui::propD("Poisson's ratio", &M.slipNu, "%.3f");
            ch |= ui::propD("Friction drop when slipping dMu", &M.slipDMu, "%.3f", "Slip weakening: a face that slips loses dMu of its friction (0 = stable sliding, no weakening)");
            ch |= ui::propD("Reference slip for the T factor (mm)", &M.slipRefMm, "%.4g", "Two-way coupling: the transmissibility factor reaches 'T factor when slipping' at this slip");
            ch |= ui::propD("Interaction cut-off (element sizes, 0 = all)", &M.slipCutoff, "%.4g", "Pairs further apart than this many mean element sizes are dropped (sparse matrix)");
            static const char* kHs[] = {"Full space", "Half space (free surface at depth 0)"};
            ch |= ui::propCombo("Elastic domain", &M.slipHalfSpace, kHs, 2, "Half space: traction-free ground surface at depth 0 (depths must be positive)");
            ui::endProps();
            M.slipNu = std::clamp(M.slipNu, 0.0, 0.49); M.slipE = std::max(0.01, M.slipE); M.slipDMu = std::clamp(M.slipDMu, 0.0, 1.0);
        }
        if (ch) markDirty(true);
    }

    // ---- G1-c contact coupling of the hydraulic fractures in the run
    if (!gi.hfracs.empty()) {
        ui::header("Fracture contact coupling in the run (G1-c, GEOCOUPLE)");
        GeoCoupleSpec& G = case_.geoCouple;
        bool ch = ImGui::Checkbox("Update hydraulic-fracture conductance from an EFEM contact solve after every report step", &G.on);
        if (G.on && ui::beginProps("gcin")) {
            char lab[64];
            std::snprintf(lab, sizeof lab, "Far-field Sxx (%s, compression +)", pU); ch |= ui::propD(lab, &G.sxx, "%.4g");
            std::snprintf(lab, sizeof lab, "Far-field Syy (%s, compression +)", pU); ch |= ui::propD(lab, &G.syy, "%.4g");
            ch |= ui::propD("Young's modulus (GPa)", &G.eGPa, "%.4g");
            ch |= ui::propD("Poisson's ratio", &G.nu, "%.3f");
            ch |= ui::propD("Fracture friction", &G.friction, "%.3f");
            ch |= ui::propI("Contact slice layer", &G.layer, 1, std::max(1, gi.nz));
            ui::endProps();
        }
        ui::hint("Total-stress EFEM slice of one layer, one report step lag, black-oil engine; results in Faults & seismicity -> Fracture contact.");
        if (ch) markDirty(true);
    }
}

// ============================================================================ 3D panels
std::vector<PanelQuad> App::faultPanels(const Case& c, const Grid* g, const RunResults* R, int snap) {
    std::vector<PanelQuad> out;
    const faultmech::Output* O = (R && R->faultMech) ? R->faultMech.get() : nullptr;
    faultmech::Output pre;
    if (!O) {   // model mode: the fault faces as scene objects (no values)
        if (!c.faultMech.on() || !g) return out;
        faultmech::Monitor mon; std::string err;
        if (!mon.setup(c, *g, err)) return out;
        pre = mon.out();
        std::vector<bool> named(pre.names.size(), false);
        for (const auto& f : pre.faces) {
            PanelQuad q;
            for (int k = 0; k < 4; ++k) q.p[k] = f.quad[size_t(k)];
            q.rgba[0] = 0.95f; q.rgba[1] = 0.45f; q.rgba[2] = 0.20f; q.rgba[3] = 0.55f;
            if (!named[size_t(f.fault)]) { named[size_t(f.fault)] = true; q.label = pre.names[size_t(f.fault)]; }
            out.push_back(std::move(q));
        }
        return out;
    }
    if (O->steps.empty()) return out;
    const double tNow = (R && !R->snaps.empty()) ? R->snaps[size_t(std::clamp(snap, 0, int(R->snaps.size()) - 1))].time : O->steps.back().t;
    const size_t si = stepAt(*O, tNow);
    const auto& S = O->steps[si];
    const int prop = std::clamp(faultProp_, 0, kNumFaultProps - 1);
    // colour range fixed over the whole run (does not jump between report steps)
    double lo = 1e300, hi = -1e300;
    for (const auto& st : O->steps) for (size_t i = 0; i < O->faces.size(); ++i) { const double v = faultValue(*O, st, i, prop); if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); } }
    if (!(hi > lo)) { hi = lo + 1; }
    if (divergingProp(prop)) { const double m = std::max(std::fabs(lo), std::fabs(hi)); lo = -m; hi = m; }
    faultRange_[0] = lo; faultRange_[1] = hi;
    std::vector<bool> named(O->names.size(), false);
    for (size_t i = 0; i < O->faces.size(); ++i) {
        const auto& f = O->faces[i];
        PanelQuad q;
        for (int k = 0; k < 4; ++k) q.p[k] = f.quad[size_t(k)];
        const double v = faultValue(*O, S, i, prop);
        float rgb[3] = {0.5f, 0.5f, 0.5f};
        if (std::isfinite(v)) colormapRGB(divergingProp(prop) ? Colormap::CoolWarm : Colormap::Turbo, float(std::clamp((v - lo) / (hi - lo), 0.0, 1.0)), rgb);
        q.rgba[0] = rgb[0]; q.rgba[1] = rgb[1]; q.rgba[2] = rgb[2]; q.rgba[3] = 0.92f;
        if (f.fault >= 0 && size_t(f.fault) < named.size() && !named[size_t(f.fault)]) {
            named[size_t(f.fault)] = true;
            const auto& Se = O->series[size_t(f.fault)];
            const size_t k = std::min(si, Se.t.empty() ? 0 : Se.t.size() - 1);
            char b[160];
            if (!Se.t.empty()) std::snprintf(b, sizeof b, "%s  %.0f%% slipping, %.3g events", O->names[size_t(f.fault)].c_str(), 100 * Se.reactFrac[k], Se.expected[k]);
            else std::snprintf(b, sizeof b, "%s", O->names[size_t(f.fault)].c_str());
            q.label = b;
        }
        out.push_back(std::move(q));
    }
    if (faultEvents_ && !O->events.empty()) {
        // events up to the current time: two crossed squares, size ~ rupture radius (Eshelby, 3 MPa stress drop) x scale,
        // floored at 15 % of a face so that small events stay visible; the last report interval highlighted
        double faceL = 0; for (const auto& f : O->faces) faceL += std::sqrt(std::max(0.0, f.area)); faceL /= double(std::max<size_t>(1, O->faces.size()));
        const double tPrev = si > 0 ? O->steps[si - 1].t : -1e300;
        size_t shown = 0;
        for (const auto& e : O->events) {
            if (e.t > tNow + 1e-9) continue;
            if (++shown > 20000) break;
            const double M0 = std::pow(10.0, 1.5 * e.mw + 9.1);
            const double r = std::cbrt(7.0 * M0 / (16.0 * 3e6)) / O->mPerUnit;
            const float s = float(std::max(0.15 * faceL, r) * 0.5 * faultEventScale_);
            const bool recent = e.t > tPrev;
            float col[4] = {1.0f, 0.85f, 0.15f, 0.95f};
            if (e.kind == 1) { col[0] = 1.0f; col[1] = 0.2f; col[2] = 0.9f; }
            else if (!recent) { col[0] = 0.75f; col[1] = 0.75f; col[2] = 0.75f; col[3] = 0.6f; }
            for (int pl = 0; pl < 2; ++pl) {
                PanelQuad q;
                const float dx = pl == 0 ? s : 0, dy = pl == 0 ? 0 : s;
                q.p[0] = {e.x - dx, e.y - dy, e.z - s}; q.p[1] = {e.x + dx, e.y + dy, e.z - s};
                q.p[2] = {e.x + dx, e.y + dy, e.z + s}; q.p[3] = {e.x - dx, e.y - dy, e.z + s};
                for (int k = 0; k < 4; ++k) q.rgba[k] = col[k];
                out.push_back(std::move(q));
            }
        }
    }
    return out;
}

// ============================================================================ results tab
void App::faultsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor) {
    if (!R.faultMech && !R.geoCouple && !R.hfDdm && !R.hfPf) { ui::hint("No coupled fault mechanics in this run."); return; }
    const float fs = ImGui::GetFontSize();
    const bool fld = R.input.units == UnitSystem::Field;
    const char* pU = fld ? "psi" : "bar";
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
    const faultmech::Output* O = R.faultMech.get();
    // 3D controls (main view)
    if (O) {
        ImGui::SetNextItemWidth(fs * 16);
        ImGui::Combo("3D fault colour", &faultProp_, kFaultProps, kNumFaultProps);
        ImGui::SameLine(); ImGui::Checkbox("Events in 3D", &faultEvents_);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("Event size", &faultEventScale_, 0.2f, 5.0f, "%.1f");
        ImGui::SameLine(); ImGui::TextDisabled("colour range %.3g .. %.3g (fixed over the run)", faultRange_[0], faultRange_[1]);
        // KPI strip
        const double tNow = tSnap >= 0 ? tSnap / xs : (O->steps.empty() ? 0 : O->steps.back().t);
        for (size_t f = 0; f < O->series.size(); ++f) {
            const auto& S = O->series[f];
            if (S.t.empty()) continue;
            size_t k = 0; for (size_t i = 0; i < S.t.size(); ++i) if (std::fabs(S.t[i] - tNow) < std::fabs(S.t[k] - tNow)) k = i;
            ImGui::Text("%s @ %.4g d:  max dCFS %.3g %s | max CFF %.3g %s | max slip tendency %.3f | %.0f %% slipping | expected events %.4g, catalogue %.0f, Mw max %.2f (cap %.2f)",
                        O->names[f].c_str(), S.t[k], S.maxDcfs[k], pU, S.maxCff[k], pU, S.maxSt[k], 100 * S.reactFrac[k], S.expected[k], S.events[k], S.maxMw[k],
                        f < O->mMax.size() ? O->mMax[f] : 0.0);
        }
        if (O->twoWayUpdates > 0) ImGui::TextDisabled("two-way coupling: %ld face multiplier updates applied to the flow model", O->twoWayUpdates);
    }
    if (!ImGui::BeginTabBar("faultsub")) return;
    if (O && !O->series.empty() && ui::tabItem("Time series")) {
        const ImVec2 av = ImGui::GetContentRegionAvail();
        bool anySlip = false;
        for (const auto& S : O->series) for (double v : S.maxSlip) anySlip |= v > 0;
        anySlip |= O->slip;
        const PlotGrid G = fitGrid(anySlip ? 8 : 6, av, fs * 22, fs * 12, 1.6f);
        if (G.scroll) ImGui::BeginChild("##fmscroll", av);
        fig::recorder().beginFigure("run.faults", G.rows, G.cols);
        if (ImPlot::BeginSubplots("##fm", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
            auto tx = [&](const std::vector<double>& t) { std::vector<double> x(t); for (double& v : x) v *= xs; return x; };
            char yl[64];
            std::snprintf(yl, sizeof yl, "Stress (%s)", pU);
            if (ImPlot::BeginPlot("Coulomb stress on the faults")) {
                recPanel("Coulomb stress on the faults", xl, yl); ImPlot::SetupAxes(xl, yl, af, af);
                for (size_t f = 0; f < O->series.size(); ++f) {
                    const auto& S = O->series[f]; char l[96];
                    std::snprintf(l, sizeof l, "%s max dCFS", O->names[f].c_str()); plotXY(l, tx(S.t), S.maxDcfs, seriesColor(int(2 * f)), 2.0f);
                    std::snprintf(l, sizeof l, "%s max CFF", O->names[f].c_str()); plotXY(l, tx(S.t), S.maxCff, seriesColor(int(2 * f + 1)), 1.6f);
                }
                hline("##cff0", 0, IM_COL32(255, 255, 255, 90));
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Slip tendency")) {
                recPanel("Slip tendency", xl, "max tau / sigma'n"); ImPlot::SetupAxes(xl, "max tau / sigma'n", af, af);
                for (size_t f = 0; f < O->series.size(); ++f) plotXY(O->names[f].c_str(), tx(O->series[f].t), O->series[f].maxSt, seriesColor(int(f)), 2.0f);
                if (!R.input.faultMech.faults.empty()) hline("##mu", R.input.faultMech.faults[0].mu, IM_COL32(255, 120, 80, 140));
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Reactivated fault area")) {
                recPanel("Reactivated fault area", xl, "area with CFF >= 0 (%)"); ImPlot::SetupAxes(xl, "area with CFF >= 0 (%)", af, af);
                for (size_t f = 0; f < O->series.size(); ++f) { std::vector<double> y(O->series[f].reactFrac); for (double& v : y) v *= 100; plotXY(O->names[f].c_str(), tx(O->series[f].t), y, seriesColor(int(f)), 2.0f); }
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Cumulative events")) {
                recPanel("Cumulative events", xl, "events >= M min"); ImPlot::SetupAxes(xl, "events >= M min", af, af);
                for (size_t f = 0; f < O->series.size(); ++f) {
                    const auto& S = O->series[f]; char l[96];
                    std::snprintf(l, sizeof l, "%s expected (Dieterich)", O->names[f].c_str()); plotXY(l, tx(S.t), S.expected, seriesColor(int(2 * f)), 2.2f);
                    std::snprintf(l, sizeof l, "%s catalogue", O->names[f].c_str()); plotXY(l, tx(S.t), S.events, seriesColor(int(2 * f + 1)), 1.4f);
                }
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Event magnitudes")) {
                recPanel("Event magnitudes", xl, "Mw"); ImPlot::SetupAxes(xl, "Mw", af, af);
                std::vector<double> t0, m0, t1, m1, t2, m2;
                for (const auto& e : O->events) {
                    if (e.kind == 0) { t0.push_back(e.t * xs); m0.push_back(e.mw); }
                    else if (e.kind == 1) { t1.push_back(e.t * xs); m1.push_back(e.mw); }
                    else { t2.push_back(e.t * xs); m2.push_back(e.mw); }
                }
                plotXY("Dieterich + Gutenberg-Richter", t0, m0, IM_COL32(255, 200, 40, 255), 1.0f, true);
                if (!t1.empty()) plotXY("rate-and-state ruptures", t1, m1, IM_COL32(255, 60, 220, 255), 1.0f, true);
                if (!t2.empty()) plotXY("DDM slip events (FAULTSLIP)", t2, m2, IM_COL32(80, 230, 255, 255), 1.0f, true);
                for (size_t f = 0; f < O->mMax.size(); ++f) hline("##mcap", O->mMax[f], IM_COL32(255, 90, 90, 120));
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Magnitude - frequency")) {
                recPanel("Magnitude - frequency", "Mw", "N (>= Mw)");
                ImPlot::SetupAxes("Mw", "N (>= Mw)", af, af);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                std::vector<double> m; for (const auto& e : O->events) if (e.kind == 0) m.push_back(e.mw);
                std::sort(m.begin(), m.end());
                std::vector<double> x, y; for (size_t i = 0; i < m.size(); ++i) { x.push_back(m[i]); y.push_back(double(m.size() - i)); }
                plotXY("catalogue", x, y, IM_COL32(255, 200, 40, 255), 1.0f, true);
                if (!m.empty() && !R.input.faultMech.faults.empty()) {   // input b-value through N(M min)
                    const double b = R.input.faultMech.faults[0].bValue, mm = R.input.faultMech.faults[0].mMin;
                    std::vector<double> gx{mm, m.back()}, gy{double(m.size()), double(m.size()) * std::pow(10.0, -b * (m.back() - mm))};
                    plotXY("G-R, input b", gx, gy, IM_COL32(120, 200, 255, 255), 1.6f);
                }
                recEnd(); ImPlot::EndPlot();
            }
            if (anySlip) {
                if (ImPlot::BeginPlot("Fault slip (DDM)")) {
                    recPanel("Fault slip (DDM)", xl, "max cumulative slip (mm)"); ImPlot::SetupAxes(xl, "max cumulative slip (mm)", af, af);
                    for (size_t f = 0; f < O->series.size(); ++f) { std::vector<double> y(O->series[f].maxSlip); for (double& v : y) v *= 1e3; plotXY(O->names[f].c_str(), tx(O->series[f].t), y, seriesColor(int(f)), 2.0f); }
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Seismic moment")) {
                    recPanel("Seismic moment", xl, "M0 (N m)"); ImPlot::SetupAxes(xl, "M0 (N m)", af, af);
                    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                    for (size_t f = 0; f < O->series.size(); ++f) {
                        std::vector<double> t, y;
                        for (size_t k = 0; k < O->series[f].moment.size() && k < O->series[f].t.size(); ++k) if (O->series[f].moment[k] > 0) { t.push_back(O->series[f].t[k] * xs); y.push_back(O->series[f].moment[k]); }
                        plotXY(O->names[f].c_str(), t, y, seriesColor(int(f)), 2.0f);
                    }
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        if (G.scroll) ImGui::EndChild();
        ImGui::EndTabItem();
    }
    if (O && !O->faces.empty() && ui::tabItem("Fault face map")) {
        static int fsel = 0;
        fsel = std::clamp(fsel, 0, int(O->names.size()) - 1);
        ImGui::SetNextItemWidth(fs * 10);
        if (ImGui::BeginCombo("Fault", O->names[size_t(fsel)].c_str())) {
            for (int f = 0; f < int(O->names.size()); ++f) if (ImGui::Selectable(O->names[size_t(f)].c_str(), f == fsel)) fsel = f;
            ImGui::EndCombo();
        }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 16); ImGui::Combo("Quantity", &faultProp_, kFaultProps, kNumFaultProps);
        // along-strike x layer map of the selected fault at the current report step
        std::vector<size_t> idx; double nx = 0, ny = 0;
        for (size_t i = 0; i < O->faces.size(); ++i) if (O->faces[i].fault == fsel) { idx.push_back(i); nx += O->faces[i].n[0]; ny += O->faces[i].n[1]; }
        if (idx.empty() || O->steps.empty()) ui::hint("No faces.");
        else {
            const double nl = std::hypot(nx, ny), tx = nl > 0 ? -ny / nl : 1, ty = nl > 0 ? nx / nl : 0;
            const auto& c0 = O->faces[idx[0]].c;
            std::vector<double> s(idx.size()); int kmin = 1 << 30, kmax = -1;
            for (size_t j = 0; j < idx.size(); ++j) { const auto& c = O->faces[idx[j]].c; s[j] = (c[0] - c0[0]) * tx + (c[1] - c0[1]) * ty; kmin = std::min(kmin, O->faces[idx[j]].k); kmax = std::max(kmax, O->faces[idx[j]].k); }
            std::vector<double> us(s); std::sort(us.begin(), us.end());
            const double tol = 0.1 * std::sqrt(std::max(1e-12, O->faces[idx[0]].area)) + 1e-6;
            std::vector<double> cols; for (double v : us) if (cols.empty() || v - cols.back() > tol) cols.push_back(v);
            const int nc = int(cols.size()), nr = kmax - kmin + 1;
            const double tNow = tSnap >= 0 ? tSnap / xs : O->steps.back().t;
            const auto& S = O->steps[stepAt(*O, tNow)];
            std::vector<double> V(size_t(nr) * size_t(nc), NAN);
            double lo = faultRange_[0], hi = faultRange_[1];
            for (size_t j = 0; j < idx.size(); ++j) {
                const int ci = int(std::lower_bound(cols.begin(), cols.end(), s[j] - tol) - cols.begin());
                const int ri = O->faces[idx[j]].k - kmin;
                if (ci >= 0 && ci < nc) V[size_t(ri) * size_t(nc) + size_t(ci)] = faultValue(*O, S, idx[j], faultProp_);
            }
            for (double& v : V) if (!std::isfinite(v)) v = lo;
            const double sL = (cols.back() - cols.front()) * O->mPerUnit;
            ImPlot::PushColormap(divergingProp(faultProp_) ? ImPlotColormap_RdBu : ImPlotColormap_Viridis);
            const ImVec2 av = ImGui::GetContentRegionAvail();
            if (ImPlot::BeginPlot("##fmap", ImVec2(av.x - fs * 6, av.y))) {
                ImPlot::SetupAxes("along strike (m)", "layer", 0, ImPlotAxisFlags_Invert);
                ImPlot::PlotHeatmap("##h", V.data(), nr, nc, lo, hi, nullptr, ImPlotPoint(0, kmin + 0.5), ImPlotPoint(std::max(1.0, sL), kmax + 1.5));
                ImPlot::EndPlot();
            }
            ImGui::SameLine();
            ImPlot::ColormapScale("##fscale", lo, hi, ImVec2(fs * 5, av.y));
            ImPlot::PopColormap();
        }
        ImGui::EndTabItem();
    }
    if (O && ui::tabItem("Catalogue")) {
        ImGui::Text("%zu events (%d critical faces at t0, %zu fault faces)", O->events.size(), O->criticalAtStart, O->faces.size());
        for (const auto& r : O->rsf) {
            if (!r.error.empty()) ImGui::TextColored(ImVec4(1, 0.5f, 0.3f, 1), "rate-and-state, fault %s layer %d: %s", O->names[size_t(r.fault)].c_str(), r.layer, r.error.c_str());
            else ImGui::Text("rate-and-state, fault %s layer %d: %d faces, h* %.4g m, Lb %.4g m, dx %.4g m (Lb/dx %.1f), %d ruptures, %.2f s",
                             O->names[size_t(r.fault)].c_str(), r.layer, r.faces, r.hStar, r.Lb, r.dx, r.dx > 0 ? r.Lb / r.dx : 0.0, r.events, r.seconds);
        }
        if (ImGui::BeginTable("##cat", 7, ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV, ImVec2(0, -FLT_MIN))) {
            ImGui::TableSetupScrollFreeze(0, 1);
            for (const char* h : {"Time (d)", "Fault", "X", "Y", "Depth", "Mw", "Model"}) ImGui::TableSetupColumn(h);
            ImGui::TableHeadersRow();
            ImGuiListClipper cl; cl.Begin(int(O->events.size()));
            while (cl.Step()) for (int i = cl.DisplayStart; i < cl.DisplayEnd; ++i) {
                const auto& e = O->events[size_t(i)];
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::Text("%.4g", e.t);
                ImGui::TableNextColumn(); ImGui::TextUnformatted(e.fault >= 0 && size_t(e.fault) < O->names.size() ? O->names[size_t(e.fault)].c_str() : "?");
                ImGui::TableNextColumn(); ImGui::Text("%.1f", e.x);
                ImGui::TableNextColumn(); ImGui::Text("%.1f", e.y);
                ImGui::TableNextColumn(); ImGui::Text("%.1f", e.z);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", e.mw);
                ImGui::TableNextColumn(); ImGui::TextUnformatted(e.kind == 0 ? "Dieterich" : e.kind == 1 ? "rate-and-state" : "DDM slip");
            }
            ImGui::EndTable();
        }
        ImGui::EndTabItem();
    }
    if (R.hfPf && !R.hfPf->empty() && ui::tabItem("Fracture growth (phase field)")) {
        // C28 in-run phase field: half-length, crack volume and pressures per report step (one series per fracture)
        const auto& V = *R.hfPf;
        int nf = 0; for (const auto& r : V) nf = std::max(nf, r.frac + 1);
        const double toPa = R.input.units == UnitSystem::Field ? 6894.757 : 1e5;   // case pressure unit -> Pa
        const auto& last = V.back();
        ImGui::Text("@ %.4g d: half-length %.4g m, crack volume %.4g m2/m, phase-field pressure %.4g MPa (closure %.4g, flow model %.4g MPa)%s%s",
                    last.t, last.halfLength, last.volume, last.pressure * 1e-6, last.sigma * toPa * 1e-6, last.pFlow * toPa * 1e-6,
                    last.active ? "" : " - stopped: ", last.active ? "" : last.note.c_str());
        const ImVec2 av = ImGui::GetContentRegionAvail();
        const PlotGrid G = fitGrid(3, av, fs * 22, fs * 12, 1.6f);
        fig::recorder().beginFigure("run.hfpf", G.rows, G.cols);
        if (ImPlot::BeginSubplots("##pf", G.rows, G.cols, G.size, 0)) {
            const char* titles[3] = {"Phase-field half-length", "Crack volume per height", "Fracture pressure"};
            const char* units[3] = {"m", "m2/m", "MPa"};
            for (int pl = 0; pl < 3; ++pl) {
                if (!ImPlot::BeginPlot(titles[pl])) continue;
                recPanel(titles[pl], xl, units[pl]); ImPlot::SetupAxes(xl, units[pl], af, af);
                for (int f = 0; f < nf; ++f) {
                    std::vector<double> t, y, yf, ys;
                    for (const auto& r : V) if (r.frac == f && r.active && r.volume > 0) {
                        t.push_back(r.t * xs);
                        y.push_back(pl == 0 ? r.halfLength : pl == 1 ? r.volume : r.pressure * 1e-6);
                        yf.push_back(r.pFlow * toPa * 1e-6); ys.push_back(r.sigma * toPa * 1e-6);
                    }
                    const std::string nm = f < int(R.input.grid.hfracs.size()) ? R.input.grid.hfracs[size_t(f)].name : "F" + std::to_string(f + 1);
                    plotXY((pl == 2 ? nm + " phase field" : nm).c_str(), t, y, seriesColor(3 * f), 2.0f);
                    if (pl == 2) { plotXY((nm + " flow model").c_str(), t, yf, seriesColor(3 * f + 1), 1.6f); plotXY((nm + " closure").c_str(), t, ys, seriesColor(3 * f + 2), 1.2f); }
                }
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    if (R.hfDdm && !R.hfDdm->steps.empty() && ui::tabItem("Fracture apertures (DDM)")) {
        const auto& D = *R.hfDdm;
        std::vector<double> t, wm, vol, lo, hi, op, pn;
        for (const auto& s : D.steps) { t.push_back(s.t * xs); wm.push_back(s.maxWidth * 1e3); vol.push_back(s.volume); lo.push_back(s.minMult); hi.push_back(s.maxMult); op.push_back(s.open); pn.push_back(s.maxPnet); }
        // along-fracture aperture profile at the current report step (pieces in input order)
        const double tNow = tSnap >= 0 ? tSnap / xs : D.steps.back().t;
        size_t k = 0; for (size_t i = 0; i < D.steps.size(); ++i) if (std::fabs(D.steps[i].t - tNow) < std::fabs(D.steps[k].t - tNow)) k = i;
        std::vector<double> px, pw;
        for (size_t i = 0; i < D.steps[k].width.size(); ++i) if (D.steps[k].width[i] > 0) { px.push_back(double(i)); pw.push_back(D.steps[k].width[i] * 1e3); }
        ImGui::Text("@ %.4g d: %d / %d pieces open, max aperture %.3g mm, opening volume %.4g m3 (%d DDM updates)", D.steps[k].t, D.steps[k].open, D.elements, D.steps[k].maxWidth * 1e3, D.steps[k].volume, D.updates);
        const ImVec2 av = ImGui::GetContentRegionAvail();
        const PlotGrid G = fitGrid(4, av, fs * 22, fs * 12, 1.6f);
        fig::recorder().beginFigure("run.hfddm", G.rows, G.cols);
        if (ImPlot::BeginSubplots("##hd", G.rows, G.cols, G.size, 0)) {
            if (ImPlot::BeginPlot("Maximum aperture")) {
                recPanel("Maximum aperture", xl, "mm"); ImPlot::SetupAxes(xl, "mm", af, af);
                plotXY("max DDM aperture", t, wm, seriesColor(0), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Opening volume")) {
                recPanel("Opening volume", xl, "m3"); ImPlot::SetupAxes(xl, "m3", af, af);
                plotXY("DDM opening volume", t, vol, seriesColor(1), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Conductance multiplier (DDM)")) {
                recPanel("Conductance multiplier (DDM)", xl, "x"); ImPlot::SetupAxes(xl, "x", af, af);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                plotXY("min", t, lo, seriesColor(2), 2.0f); plotXY("max", t, hi, seriesColor(3), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Aperture per piece (current step)")) {
                recPanel("Aperture per piece (current step)", "piece index", "mm"); ImPlot::SetupAxes("piece index", "mm", af, af);
                plotXY("aperture", px, pw, seriesColor(4), 1.6f, true);
                recEnd(); ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    if (R.geoCouple && !R.geoCouple->steps.empty() && ui::tabItem("Fracture contact (G1-c)")) {
        const auto& Gc = *R.geoCouple;
        std::vector<double> t, wn, wt, lo, hi, op, sl;
        for (const auto& s : Gc.steps) { t.push_back(s.t * xs); wn.push_back(s.maxOpening * 1e3); wt.push_back(s.maxSlip * 1e3); lo.push_back(s.minMult); hi.push_back(s.maxMult); op.push_back(s.openPieces); sl.push_back(s.slipPieces); }
        const ImVec2 av = ImGui::GetContentRegionAvail();
        const PlotGrid G = fitGrid(3, av, fs * 22, fs * 12, 1.6f);
        fig::recorder().beginFigure("run.geocouple", G.rows, G.cols);
        if (ImPlot::BeginSubplots("##gc", G.rows, G.cols, G.size, ImPlotSubplotFlags_LinkAllX)) {
            if (ImPlot::BeginPlot("Fracture opening and slip")) {
                recPanel("Fracture opening and slip", xl, "mm"); ImPlot::SetupAxes(xl, "mm", af, af);
                plotXY("max opening", t, wn, seriesColor(0), 2.0f); plotXY("max slip", t, wt, seriesColor(1), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Conductance multiplier")) {
                recPanel("Conductance multiplier", xl, "x"); ImPlot::SetupAxes(xl, "x", af, af);
                plotXY("min", t, lo, seriesColor(2), 2.0f); plotXY("max", t, hi, seriesColor(3), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Contact state")) {
                recPanel("Contact state", xl, "pieces"); ImPlot::SetupAxes(xl, "pieces", af, af);
                plotXY("open", t, op, seriesColor(4), 2.0f); plotXY("slipping", t, sl, seriesColor(5), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}
