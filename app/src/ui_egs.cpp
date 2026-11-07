// C13: geothermal / EDFM user interface
//  * egsEditor(): "Fractures & EGS" tab of the thermal model editor - interactive plan view (draw / drag fracture traces
//    over the grid columns and the wells), vertical section of the selected fracture (drag its depth range), the
//    fracture table with Barton-Bandis / Coulomb parameters, the in-situ stress, an EDFM pre-check (fracture cells,
//    well - fracture connections, initial effective normal stress and slip tendency per fracture) and the multi-stage
//    EGS designer (makeEgsCase).
//  * fracturePanels(): fracture planes / segments for the 3D viewers (coloured by the fracture value of the current
//    property when it has one).
//  * geothermalTab(): EGS performance analysis of a thermal run - production temperature and normalised thermal
//    drawdown with the lifetime threshold, thermal / electric power, cumulative heat and heat recovery, injectivity /
//    productivity / impedance, fracture aperture / temperature / conductivity per fracture, cooled rock volume,
//    fracture-face maps, along-fracture profiles, plan-view cold-plume map and the well - fracture connections.
#include "app.h"
#include "platform.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "plotkit.h"
#include "implot.h"
#include "ressim/examples.h"
#include "ressim/therm_deck.h"
#include "ressim/thermsim.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <map>

using namespace ressim;

extern Colormap g_fvMap;   // ui_fieldviz.cpp
extern bool g_fvRev;
namespace {
EgsDesign g_egs;   // designer inputs of this session

struct FracCheck {
    bool valid = false;
    unsigned version = 0;
    std::string error;
    int fracCells = 0, ffConns = 0, wellConns = 0;
    struct Row { std::string name; int segs = 0; double area = 0, snEff = 0, tau = 0, slip = 0, p = 0; std::vector<std::string> wells; };
    std::vector<Row> rows;
};
FracCheck g_check;
int g_selFrac = 0;
bool g_drawMode = false, g_drawing = false, g_snap = true;
double g_drawX0 = 0, g_drawY0 = 0;

void plotXYe(const char* label, const std::vector<double>& X, const std::vector<double>& Y, ImU32 col, float w, bool dashed = false) {
    const int n = int(std::min(X.size(), Y.size()));
    if (n <= 0) return;
    if (fig::recorder().recording()) recSeries(label, X.data(), Y.data(), n, col, w);
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = w * g_plotLineScale;
    if (dashed) { sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 2.5f; }
    ImPlot::PlotLine(label, X.data(), Y.data(), n, sp);
}
void hline(const char* label, double y, ImU32 col) {
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = 1.0f; sp.Flags = ImPlotItemFlags_NoLegend;
    ImPlot::PlotInfLines(label, &y, 1, sp);
    (void)sp;
}
ImU32 cmapU32(float t, float a = 1.0f) {   // C30: the colour map chosen in Results -> Fields & sections
    float rgb[3]; colormapRGB(g_fvMap, t, rgb, g_fvRev);
    return ImGui::ColorConvertFloat4ToU32(ImVec4(rgb[0], rgb[1], rgb[2], a));
}
// horizontal colour bar with min / max labels (Turbo, as the 3D viewer default)
void colorBar(const char* label, double lo, double hi, float width) {
    const float fs = ImGui::GetFontSize();
    ImGui::TextDisabled("%s", label);
    ImGui::SameLine();
    char b[32]; std::snprintf(b, sizeof(b), "%.4g", lo);
    ImGui::TextUnformatted(b); ImGui::SameLine();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float h = fs * 0.8f;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const int n = 48;
    for (int i = 0; i < n; ++i) {
        const float x0 = p.x + width * float(i) / n, x1 = p.x + width * float(i + 1) / n;
        dl->AddRectFilled(ImVec2(x0, p.y + fs * 0.15f), ImVec2(x1, p.y + fs * 0.15f + h), cmapU32((float(i) + 0.5f) / n));
    }
    ImGui::Dummy(ImVec2(width, fs));
    ImGui::SameLine();
    std::snprintf(b, sizeof(b), "%.4g", hi);
    ImGui::TextUnformatted(b);
}
// grid column edges (x along I, y along J) of the top layer in metres
void columnEdges(const Grid& g, double lenU, std::vector<double>& xe, std::vector<double>& ye) {
    xe.clear(); ye.clear();
    if (g.corners.empty() || g.nx <= 0 || g.ny <= 0) return;
    for (int i = 0; i < g.nx; ++i) xe.push_back(g.corners[size_t(i)][0].x * lenU);
    xe.push_back(g.corners[size_t(g.nx - 1)][1].x * lenU);
    for (int j = 0; j < g.ny; ++j) ye.push_back(g.corners[size_t(j) * size_t(g.nx)][0].y * lenU);
    ye.push_back(g.corners[size_t(g.ny - 1) * size_t(g.nx)][2].y * lenU);
}
double distPointSeg(double px, double py, double ax, double ay, double bx, double by) {
    const double dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy;
    double t = L2 > 0 ? ((px - ax) * dx + (py - ay) * dy) / L2 : 0;
    t = std::max(0.0, std::min(1.0, t));
    const double qx = ax + t * dx - px, qy = ay + t * dy - py;
    return std::sqrt(qx * qx + qy * qy);
}
} // namespace

// ============================================================================ editor: Fractures & EGS
bool App::egsEditor() {
    ThermalInput& T = case_.thermal;
    const bool fld = case_.units == UnitSystem::Field;
    const double lenU = fld ? 0.3048 : 1.0;   // case length unit -> m
    const float fs = ImGui::GetFontSize();
    auto changed = [&]() { T.structured = true; markDirty(); g_check.valid = false; };

    // ---------------------------------------------------------------- how-to
    if (ImGui::CollapsingHeader(ICON_INFO " How EDFM / EGS works in ResSim", T.fractures.empty() ? ImGuiTreeNodeFlags_DefaultOpen : 0)) {
        ui::hint("1. Start from a thermal water model (Components: water only or water + a gas tracer) on a Cartesian grid,\n"
                 "   or let the EGS designer below build a complete multi-stage doublet for you.\n"
                 "2. Add fractures: 'Draw' in the plan view (click-drag a trace), '+ Add fracture', or type coordinates in the\n"
                 "   table. Drag the end points of the selected trace; drag its top / bottom in the section view.\n"
                 "3. Each fracture is cut by the matrix columns and layers into fracture cells (EDFM, Moinfar et al. 2013):\n"
                 "   f-m transmissibility 2 A k_n / <d>, f-f along the trace and at intersections, cubic-law w^2/12.\n"
                 "4. Wells whose completions lie in a cell cut by a fracture are connected to it automatically\n"
                 "   (WI_f = 2 pi k_f w / ln(r_e / r_w)). 'Check fractures' lists the connections before you run.\n"
                 "5. Optional mechanics: Barton-Bandis closure under the effective normal stress from the stress field,\n"
                 "   thermo-elastic unloading on cooling, Coulomb hydro-shear dilation (irreversible).\n"
                 "6. Run (F5). Results: Plots > Geothermal (temperatures, drawdown, MW, impedance, fracture apertures, maps)\n"
                 "   and the 3D viewers draw the fracture planes coloured by TEMP / PRESSURE / FRACAP / DTEMP.");
    }

    // ---------------------------------------------------------------- EGS designer
    if (ImGui::CollapsingHeader(ICON_SPARKLES " EGS designer (multi-stage horizontal doublet)")) {
        ui::hint("Builds a complete thermal case: injector and producer laterals along x (S_hmin direction), transverse\n"
                 "hydraulic fractures normal to the wells, layered Cartesian grid with geometric refinement at the planes\n"
                 "(the coarse-grid early-cooling error is removed: see the Gringarten verification), in-situ stress and\n"
                 "Barton-Bandis / Coulomb fracture mechanics, quarterly reports. Replaces the current model.");
        EgsDesign& d = g_egs;
        if (ui::beginProps("egsd", 0.5f)) {
            ui::label("Layout");
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ui::propI("Stages (fractures)", &d.stages, 1, 40);
            ui::propD("Stage spacing (m)", &d.spacing, "%.4g");
            ui::propD("Injector - producer separation (m)", &d.wellSep, "%.4g");
            ui::propD("Fracture extent beyond each well (m)", &d.fracExt, "%.4g");
            ui::propD("Fracture height (m)", &d.fracHeight, "%.4g");
            double apmm = d.aperture * 1e3;
            if (ui::propD("Initial aperture (mm)", &apmm, "%.4g")) d.aperture = apmm * 1e-3;
            ui::propBool("Fracture mechanics (Barton-Bandis, hydro-shear)", &d.mech);
            ui::label("Reservoir");
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ui::propD("Depth of the laterals (m)", &d.depth, "%.5g");
            ui::propD("Initial rock temperature (\xC2\xB0""C)", &d.tRes, "%.4g");
            ui::propD("Matrix permeability (mD)", &d.kMatrix, "%.4g");
            ui::propD("Matrix porosity (-)", &d.poro, "%.4g");
            ui::propD("S_Hmax gradient (Pa/m)", &d.sHmaxGrad, "%.5g");
            ui::propD("S_hmin gradient (Pa/m)", &d.sHminGrad, "%.5g");
            ui::label("Operation");
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ui::propD("Injection rate (sm3/d water)", &d.rate, "%.5g", "about 2500 sm3/d = 29 kg/s");
            ui::propD("Injection temperature (\xC2\xB0""C)", &d.tInj, "%.4g");
            ui::propD("Producer drawdown below hydrostatic (bar)", &d.drawdown, "%.4g");
            ui::propD("Simulated years", &d.years, "%.4g");
            ui::label("Stimulation sizing (fracture propagation)");
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ui::propBool("Size the stages from a stimulation simulation", &d.stimulate,
                         "each stage = hydraulic fracture grown from the injector (PKN / KGD, Carter leak-off, toughness);\n"
                         "half-length and initial aperture replace 'extent beyond each well' / 'initial aperture'");
            ImGui::BeginDisabled(!d.stimulate);
            {
                static const char* const models[] = {"PKN (height-contained, L >> H)", "KGD (plane strain, H >> L)"};
                ui::propCombo("Propagation model", &d.stimModel, models, 2);
                ui::propD("Pump rate per stage (m3/s)", &d.stimRate, "%.4g", "0.1 m3/s ~ 38 bbl/min");
                double tmin = d.stimTime / 60;
                if (ui::propD("Pumping time per stage (min)", &tmin, "%.4g")) d.stimTime = std::max(1.0, tmin * 60);
                double mucp = d.stimMu * 1e3;
                if (ui::propD("Fluid viscosity (cP)", &mucp, "%.4g", "slickwater ~1 cP, linear gel 10-50 cP")) d.stimMu = std::max(1e-6, mucp * 1e-3);
                ui::propD("Carter leak-off coefficient (m/s^0.5)", &d.stimCl, "%.3g", "0 = no leak-off; crystalline rock ~1e-6 - 1e-5");
                double kic = d.stimKic * 1e-6;
                if (ui::propD("Fracture toughness (MPa m^0.5)", &kic, "%.3g", "granite ~1 - 2")) d.stimKic = std::max(0.0, kic * 1e6);
                ui::propD("Retained aperture / mean width (-)", &d.stimRetain, "%.3g", "residual hydraulic aperture after shut-in (shear dilation / self-propping)");
                double eg = d.rockE * 1e-9;
                if (ui::propD("Rock Young's modulus (GPa)", &eg, "%.4g")) d.rockE = std::max(0.1, eg) * 1e9;
                ui::propD("Rock Poisson ratio (-)", &d.rockNu, "%.3f");
            }
            ImGui::EndDisabled();
            ui::label("Grid");
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ui::propD("Cell size away from the fractures (m)", &d.cell, "%.4g");
            ui::propI("Layers (cap + fractured zone + base)", &d.nz, 3, 30);
            ui::propBool("Geometric refinement at the planes", &d.refine, "recommended: removes the coarse-grid early-cooling error");
            ImGui::BeginDisabled(!d.refine);
            ui::propD("Width of the cell holding a plane (m)", &d.refineFirst, "%.3g");
            ui::propD("Growth ratio (-)", &d.refineRatio, "%.3g");
            ImGui::EndDisabled();
            ui::endProps();
        }
        // C14: stimulation preview (the same simulation makeEgsCase runs)
        if (d.stimulate) {
            static hf::Result s_stim; static bool s_stimOk = false; static std::string s_stimErr;
            if (ImGui::Button(ICON_PLAY " Simulate stimulation")) {
                s_stimOk = egsStimulation(d, s_stim); s_stimErr = s_stimOk ? std::string() : s_stim.error;
            }
            ImGui::SameLine(); ui::hint("PKN: Nordgren (1972) / KGD: Adachi-Detournay (2002) verified, see -hf tests");
            if (!s_stimErr.empty()) ImGui::TextColored(ImVec4(1, 0.4f, 0.3f, 1), "%s", s_stimErr.c_str());
            if (s_stimOk && !s_stim.snaps.empty()) {
                const hf::Snapshot& S = s_stim.last();
                const double ap = d.stimRetain * s_stim.meanWidth();
                ImGui::Text("Half-length %.1f m  |  wellbore width %.2f mm  |  mean width %.2f mm  |  net pressure %.2f MPa  |  BHP-equivalent %.1f MPa  |  efficiency %.2f",
                            S.L, S.w0 * 1e3, s_stim.meanWidth() * 1e3, S.pnet0 * 1e-6, (d.sHminGrad * d.depth + S.pnet0) * 1e-6, S.eff);
                ImGui::Text("Resulting plane: injector +- %.1f m, initial aperture %.3f mm  ->  %s", S.L, ap * 1e3,
                            S.L > d.wellSep ? "reaches the producer" : "does NOT reach the producer (pump longer / faster or reduce leak-off)");
                const float hgt = 220 * ImGui::GetIO().FontGlobalScale;
                if (ImGui::BeginTable("##stimplots", 2, ImGuiTableFlags_SizingStretchSame)) {
                    ImGui::TableNextColumn();
                    std::vector<double> tt, LL; for (const auto& q : s_stim.snaps) { tt.push_back(q.t / 60); LL.push_back(q.L); }
                    if (ImPlot::BeginPlot("Half-length vs time##stimL", ImVec2(-1, hgt))) {
                        recPanel("Stimulation: half-length vs time", "Pumping time (min)", "Half-length (m)", "");
                        ImPlot::SetupAxes("Pumping time (min)", "Half-length (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        ImPlot::PlotLine("L(t)", tt.data(), LL.data(), int(tt.size()));
                        ImPlot::EndPlot();
                    }
                    ImGui::TableNextColumn();
                    if (ImPlot::BeginPlot("Width profile at shut-in##stimW", ImVec2(-1, hgt))) {
                        recPanel("Stimulation: width profile at shut-in", "Distance from the well (m)", "Width (mm)", "");
                        ImPlot::SetupAxes("Distance from the well (m)", "Width (mm)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        std::vector<double> wm; for (double v : S.w) wm.push_back(v * 1e3);
                        if (!S.x.empty()) ImPlot::PlotLine("w(x)", S.x.data(), wm.data(), int(S.x.size()));
                        ImPlot::EndPlot();
                    }
                    ImGui::EndTable();
                }
            }
        }
        // a quick estimate of the grid size
        {
            const double xLen = (d.stages - 1) * d.spacing + 2 * std::max(100.0, 2 * d.cell);
            int ncol = int(std::ceil(xLen / std::max(5.0, d.cell)));
            if (d.refine) {
                int perHalf = 0; double w = d.refineFirst, acc = 0.5 * d.refineFirst;
                while (acc < 0.5 * d.spacing && perHalf < 30) { w = std::min(w * d.refineRatio, d.cell); acc += w; ++perHalf; }
                ncol = int(std::ceil(2 * std::max(100.0, 2 * d.cell) / d.cell)) + d.stages * (1 + 2 * perHalf);
            }
            const int nyc = int(std::ceil((d.wellSep + 2 * d.fracExt + 2 * std::max(100.0, 2 * d.cell)) / std::max(5.0, d.cell)));
            ImGui::TextDisabled("Approximate grid: %d x %d x %d = %d cells + about %d fracture cells", ncol, nyc, d.nz, ncol * nyc * d.nz,
                                d.stages * std::max(1, int((d.wellSep + 2 * d.fracExt) / d.cell)) * std::max(1, d.nz - 2));
        }
        if (ui::accentButton(ICON_SPARKLES " Create EGS case")) {
            if (maybeSaveChanges()) {
                try {
                    pendingCase_ = std::make_unique<Case>(makeEgsCase(g_egs));
                    newCase(4);
                    g_check = FracCheck{};
                    notify("EGS case created - press Run (F5); results: Plots > Geothermal");
                    return false;
                } catch (const std::exception& e) { plat::msgBox(e.what(), "EGS designer", MB_OK | MB_ICONERROR); }
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset defaults")) g_egs = EgsDesign{};
        ImGui::SameLine();
        if (ImGui::Button("Load single-fracture EDFM example")) { if (maybeSaveChanges()) { newCase(3, "egs_edfm"); return false; } }
    }

    // ---------------------------------------------------------------- plan view
    ui::header("Plan view (x along I, y along J; metres from the grid corner)");
    std::vector<double> xe, ye;
    if (preview_ && preview_->nx == case_.grid.nx && preview_->ny == case_.grid.ny) columnEdges(*preview_, lenU, xe, ye);
    const bool cart = !case_.grid.cornerPoint();
    if (!cart) ui::hint("EDFM fractures need a Cartesian (block-centred) grid; this case is corner-point.");
    ImGui::Checkbox(ICON_PENCIL " Draw mode", &g_drawMode);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("click-drag in the plan view to draw a new vertical fracture trace");
    ImGui::SameLine(); ImGui::Checkbox("Snap to 1 m", &g_snap);
    ImGui::SameLine();
    if (ImGui::Button(ICON_PLUS " Add fracture")) {
        ThermalFracture F; F.name = "F" + std::to_string(T.fractures.size() + 1);
        const double lx = xe.empty() ? 100 : xe.back(), ly = ye.empty() ? 100 : ye.back();
        F.x0 = 0.2 * lx; F.x1 = 0.8 * lx; F.y0 = F.y1 = 0.5 * ly + 1e-3 * ly;
        T.fractures.push_back(F); g_selFrac = int(T.fractures.size()) - 1; changed();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(g_selFrac < 0 || g_selFrac >= int(T.fractures.size()));
    if (ImGui::Button(ICON_COPY " Duplicate") && g_selFrac < int(T.fractures.size())) {
        ThermalFracture F = T.fractures[size_t(g_selFrac)];
        const double off = xe.size() > 2 ? 2.0 * (xe[1] - xe[0]) : 20.0;
        F.x0 += off; F.x1 += off; F.name += "b";
        T.fractures.push_back(F); g_selFrac = int(T.fractures.size()) - 1; changed();
    }
    ImGui::SameLine();
    if (ui::dangerButton(ICON_TRASH_2 " Delete") && g_selFrac < int(T.fractures.size())) {
        T.fractures.erase(T.fractures.begin() + g_selFrac); g_selFrac = std::min(g_selFrac, int(T.fractures.size()) - 1); changed();
    }
    ImGui::EndDisabled();

    const float planH = std::max(fs * 16, std::min(fs * 28, ImGui::GetContentRegionAvail().x * 0.45f));
    const float avW = ImGui::GetContentRegionAvail().x;
    const bool side = avW > fs * 50;
    if (ImPlot::BeginPlot("##egsplan", ImVec2(side ? avW * 0.6f : -1, planH), ImPlotFlags_Equal | ImPlotFlags_NoMenus)) {
        ImPlot::SetupAxes("x (m)", "y (m)", 0, ImPlotAxisFlags_Invert);
        if (!xe.empty() && !ye.empty()) ImPlot::SetupAxesLimits(xe.front(), xe.back(), ye.front(), ye.back(), ImPlotCond_Once);
        // grid columns
        if (!xe.empty() && !ye.empty()) {
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            const ImU32 gc = IM_COL32(150, 160, 175, 70);
            (void)dl;
            for (double x : xe) pk::line(x, ye.front(), x, ye.back(), gc);   // C23: grid columns also in the publication export
            for (double y : ye) pk::line(xe.front(), y, xe.back(), y, gc);
            ImPlot::PopPlotClipRect();
            const double bx[2] = {xe.front(), xe.back()}, by[2] = {ye.front(), ye.back()};
            ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.Marker = ImPlotMarker_None; sp.LineColor = ImVec4(0, 0, 0, 0);
            ImPlot::PlotScatter("##bb", bx, by, 2, sp);
        }
        // wells (completion cell centres, case units -> m)
        if (preview_) {
            for (const WellInput& w : case_.wells) {
                std::vector<double> wx, wy;
                for (const Completion& c : w.completions) {
                    if (c.i < 0 || c.j < 0 || c.i >= preview_->nx || c.j >= preview_->ny) continue;
                    const int g = c.i + preview_->nx * (c.j + preview_->ny * std::max(0, std::min(c.k1, preview_->nz - 1)));
                    if (size_t(g) >= preview_->center.size()) continue;
                    wx.push_back(preview_->center[size_t(g)].x * lenU); wy.push_back(preview_->center[size_t(g)].y * lenU);
                }
                if (wx.empty()) continue;
                const bool inj = w.type != WellType::Producer;
                ImPlotSpec sp; sp.LineColor = col4(inj ? IM_COL32(70, 140, 255, 255) : IM_COL32(60, 200, 110, 255)); sp.LineWeight = 3.0f;
                sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3.5f; sp.Flags = ImPlotItemFlags_NoLegend;
                ImPlot::PlotLine(("##w" + w.name).c_str(), wx.data(), wy.data(), int(wx.size()), sp);
                ImPlot::PlotText(w.name.c_str(), wx.front(), wy.front(), ImVec2(0, -fs * 0.9f));
            }
        }
        // fractures
        for (int f = 0; f < int(T.fractures.size()); ++f) {
            ThermalFracture& F = T.fractures[size_t(f)];
            const double fx[2] = {F.x0, F.x1}, fy[2] = {F.y0, F.y1};
            const bool sel = f == g_selFrac;
            ImPlotSpec sp; sp.LineColor = col4(sel ? IM_COL32(255, 170, 40, 255) : IM_COL32(235, 110, 40, 220)); sp.LineWeight = sel ? 4.0f : 2.5f;
            ImPlot::PlotLine(F.name.c_str(), fx, fy, 2, sp);
            if (sel) {
                bool moved = false;
                double ax = F.x0, ay = F.y0, bx = F.x1, by = F.y1;
                moved |= ImPlot::DragPoint(1, &ax, &ay, ImVec4(1, 0.85f, 0.3f, 1), 6);
                moved |= ImPlot::DragPoint(2, &bx, &by, ImVec4(1, 0.85f, 0.3f, 1), 6);
                if (moved) {
                    if (g_snap) { ax = std::round(ax); ay = std::round(ay); bx = std::round(bx); by = std::round(by); }
                    F.x0 = ax; F.y0 = ay; F.x1 = bx; F.y1 = by; changed();
                }
            }
        }
        // selection / drawing
        if (ImPlot::IsPlotHovered()) {
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            if (g_drawMode) {
                if (ImGui::IsMouseClicked(0)) { g_drawing = true; g_drawX0 = mp.x; g_drawY0 = mp.y; }
            } else if (ImGui::IsMouseClicked(0)) {
                int best = -1; double bd = 1e300;
                const ImVec2 m = ImGui::GetMousePos();
                for (int f = 0; f < int(T.fractures.size()); ++f) {
                    const auto& F = T.fractures[size_t(f)];
                    const ImVec2 a = ImPlot::PlotToPixels(F.x0, F.y0), b = ImPlot::PlotToPixels(F.x1, F.y1);
                    const double d = distPointSeg(m.x, m.y, a.x, a.y, b.x, b.y);
                    if (d < bd) { bd = d; best = f; }
                }
                if (best >= 0 && bd < fs * 0.8) g_selFrac = best;
            }
        }
        if (g_drawing) {
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            const double lx[2] = {g_drawX0, mp.x}, ly[2] = {g_drawY0, mp.y};
            ImPlotSpec sp; sp.LineColor = ImVec4(1, 1, 1, 0.8f); sp.LineWeight = 2.0f; sp.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotLine("##new", lx, ly, 2, sp);
            if (ImGui::IsMouseReleased(0)) {
                g_drawing = false;
                double ax = g_drawX0, ay = g_drawY0, bx = mp.x, by = mp.y;
                if (g_snap) { ax = std::round(ax); ay = std::round(ay); bx = std::round(bx); by = std::round(by); }
                // keep the trace off the cell faces (EDFM: a plane on a face has an undefined host column)
                auto offFace = [&](double& v, const std::vector<double>& e) { for (double q : e) if (std::fabs(v - q) < 1e-6) v += 0.01; };
                if (std::fabs(ax - bx) < 1e-9) { offFace(ax, xe); bx = ax; }
                if (std::fabs(ay - by) < 1e-9) { offFace(ay, ye); by = ay; }
                if (std::hypot(bx - ax, by - ay) > 1.0) {
                    ThermalFracture F; F.name = "F" + std::to_string(T.fractures.size() + 1);
                    F.x0 = ax; F.y0 = ay; F.x1 = bx; F.y1 = by;
                    if (!T.fractures.empty()) { const auto& P = T.fractures.back(); F.aperture = P.aperture; F.mech = P.mech; F.zTop = P.zTop; F.zBot = P.zBot; }
                    T.fractures.push_back(F); g_selFrac = int(T.fractures.size()) - 1; changed();
                }
            }
        }
        ImPlot::EndPlot();
    }
    // ---------------------------------------------------------------- vertical section of the selected fracture
    if (side) ImGui::SameLine();
    if (g_selFrac >= int(T.fractures.size())) g_selFrac = int(T.fractures.size()) - 1;
    if (ImPlot::BeginPlot("##egssec", ImVec2(-1, planH), ImPlotFlags_NoMenus)) {
        ImPlot::SetupAxes("distance along the trace (m)", "depth (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_Invert | ImPlotAxisFlags_AutoFit);
        if (g_selFrac >= 0 && preview_ && !preview_->corners.empty()) {
            ThermalFracture& F = T.fractures[size_t(g_selFrac)];
            const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
            const double top = preview_->bbMin.z * lenU, bot = preview_->bbMax.z * lenU;
            // layer boundaries at the column of the trace midpoint
            int ic = 0, jc = 0;
            const double mx = 0.5 * (F.x0 + F.x1), my = 0.5 * (F.y0 + F.y1);
            for (int i = 0; i + 1 < int(xe.size()); ++i) if (mx >= xe[size_t(i)]) ic = i;
            for (int j = 0; j + 1 < int(ye.size()); ++j) if (my >= ye[size_t(j)]) jc = j;
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            for (int k = 0; k < preview_->nz; ++k) {
                const int g = ic + preview_->nx * (jc + preview_->ny * k);
                if (size_t(g) >= preview_->corners.size()) continue;
                const double z = preview_->corners[size_t(g)][0].z * lenU;
                pk::line(0, z, L, z, IM_COL32(150, 160, 175, 90));
            }
            (void)dl;
            const double zt = F.zTop < 0 ? top : F.zTop, zb = F.zBot < 0 ? bot : F.zBot;
            pk::rect(0, zt, L, zb, IM_COL32(235, 110, 40, 110), IM_COL32(255, 170, 40, 255), 2.0f, F.name.c_str());
            ImPlot::PopPlotClipRect();
            const double bx[2] = {0, L}, by[2] = {top, bot};
            ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
            ImPlot::PlotScatter("##bb", bx, by, 2, sp);
            // wells crossing the trace: perforated depth range at the projected distance
            for (const WellInput& w : case_.wells)
                for (const Completion& c : w.completions) {
                    if (c.i < 0 || c.j < 0 || c.i >= preview_->nx || c.j >= preview_->ny) continue;
                    const int g1 = c.i + preview_->nx * (c.j + preview_->ny * std::max(0, std::min(c.k1, preview_->nz - 1)));
                    const int g2 = c.i + preview_->nx * (c.j + preview_->ny * std::max(0, std::min(c.k2, preview_->nz - 1)));
                    if (size_t(g2) >= preview_->corners.size()) continue;
                    const double cx = preview_->center[size_t(g1)].x * lenU, cy = preview_->center[size_t(g1)].y * lenU;
                    const double hx = 0.5 * preview_->cellDx[size_t(g1)] * lenU, hy = 0.5 * preview_->cellDy[size_t(g1)] * lenU;
                    if (distPointSeg(cx, cy, F.x0, F.y0, F.x1, F.y1) > std::max(hx, hy) * 1.01) continue;
                    const double s = L > 0 ? ((cx - F.x0) * (F.x1 - F.x0) + (cy - F.y0) * (F.y1 - F.y0)) / L : 0;
                    const double wx[2] = {s, s}, wz[2] = {preview_->corners[size_t(g1)][0].z * lenU, preview_->corners[size_t(g2)][4].z * lenU};
                    ImPlotSpec ws; ws.LineColor = col4(w.type != WellType::Producer ? IM_COL32(70, 140, 255, 255) : IM_COL32(60, 200, 110, 255));
                    ws.LineWeight = 4.0f; ws.Flags = ImPlotItemFlags_NoLegend;
                    ImPlot::PlotLine(("##s" + w.name + std::to_string(c.k1)).c_str(), wx, wz, 2, ws);
                    ImPlot::PlotText(w.name.c_str(), s, wz[0], ImVec2(0, -fs * 0.8f));
                }
            double ztd = zt, zbd = zb;
            if (ImPlot::DragLineY(11, &ztd, ImVec4(1, 0.85f, 0.3f, 1), 2)) { F.zTop = std::max(top, std::min(ztd, zb - 0.5)); changed(); }
            if (ImPlot::DragLineY(12, &zbd, ImVec4(1, 0.85f, 0.3f, 1), 2)) { F.zBot = std::min(bot, std::max(zbd, F.zTop < 0 ? top + 0.5 : F.zTop + 0.5)); changed(); }
            ImPlot::Annotation(0.5 * L, zt, ImVec4(0.55f, 0.3f, 0.05f, 0.9f), ImVec2(0, -4), true, "%s  L = %.0f m  H = %.0f m", F.name.c_str(), L, zb - zt);
        }
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("Plan: click a trace to select it, drag its yellow end points; Draw mode: click-drag to add. Section: drag the yellow lines to set the depth range.");

    // ---------------------------------------------------------------- fracture table
    ui::header("Fractures (EDFM, Barton-Bandis, hydro-shear)");
    {
        if (!T.fractures.empty()) {   // C16: Excel-style (paste fracture lists from spreadsheets / DFN generators)
            ui::RecordTable<ThermalFracture> t;
            t.text("Name", [](ThermalFracture& F) -> std::string& { return F.name; })
             .num("x0 (m)", [](ThermalFracture& F) -> double& { return F.x0; }).num("y0 (m)", [](ThermalFracture& F) -> double& { return F.y0; })
             .num("x1 (m)", [](ThermalFracture& F) -> double& { return F.x1; }).num("y1 (m)", [](ThermalFracture& F) -> double& { return F.y1; })
             .num("zTop (m, <0 full)", [](ThermalFracture& F) -> double& { return F.zTop; })
             .num("zBot (m, <0 full)", [](ThermalFracture& F) -> double& { return F.zBot; })
             .computed("Length (m)", [](const ThermalFracture& F) { return std::hypot(F.x1 - F.x0, F.y1 - F.y0); })
             .num("Aperture (mm)", [](ThermalFracture& F) -> double& { return F.aperture; }, 1e3)
             .num("Perm (D, 0 = w2/12)", [](ThermalFracture& F) -> double& { return F.perm; }, 1.0 / 9.869233e-13)
             .flag("Mechanics", [](ThermalFracture& F) -> bool& { return F.mech; })
             .num("Kn0 (GPa/mm)", [](ThermalFracture& F) -> double& { return F.kn0; }, 1e-12)
             .num("Vm (mm)", [](ThermalFracture& F) -> double& { return F.vm; }, 1e3)
             .num("Friction mu", [](ThermalFracture& F) -> double& { return F.mu; })
             .num("Ks (GPa/m)", [](ThermalFracture& F) -> double& { return F.ks; }, 1e-9)
             .num("Dilation (deg)", [](ThermalFracture& F) -> double& { return F.dilationDeg; });
            t.onChange = [](ThermalFracture& F) { F.aperture = std::max(1e-7, F.aperture); F.perm = std::max(0.0, F.perm); F.mu = std::max(0.0, F.mu); };
            t.newRow = [&T](int at) {
                ThermalFracture F = T.fractures.empty() ? ThermalFracture{} : T.fractures[size_t(std::min(at, int(T.fractures.size()) - 1))];
                F.name = "F" + std::to_string(T.fractures.size() + 1);
                return F;
            };
            const float h = std::min(ImGui::GetFrameHeightWithSpacing() * (float(T.fractures.size()) + 1.8f), fs * 18);
            if (t.draw("edfm", T.fractures, h, &g_selFrac, "%.5g", true)) changed();
            ImGui::TextDisabled("Mechanics = stress-dependent aperture (Barton-Bandis) and Coulomb shear dilation; Kn0 ... Dilation apply when it is 'yes'.");
        }
        if (T.fractures.empty()) ui::hint("No fractures yet: draw one in the plan view, use '+ Add fracture', or the EGS designer above.");
    }
    if (!T.fractures.empty()) {
        ui::header("In-situ stress and thermo-elasticity (fracture mechanics)");
        ThermalFracStress& S = T.fracStress;
        auto propS = [&](const char* lbl, double* v, const char* fmt, const char* tip, double scale) {
            double x = *v * scale;
            if (ui::propD(lbl, &x, fmt, tip)) { *v = x / scale; changed(); }
        };
        if (ui::beginProps("fstress")) {
            propS("S_Hmax gradient (kPa/m)", &S.sHmaxGrad, "%.4g", "total maximum horizontal stress per metre depth", 1e-3);
            propS("S_hmin gradient (kPa/m)", &S.sHminGrad, "%.4g", "total minimum horizontal stress per metre depth", 1e-3);
            if (ui::propD("S_Hmax azimuth (deg from +x)", &S.azimuthDeg, "%.2f")) changed();
            propS("Young's modulus (GPa)", &S.E, "%.4g", "thermo-elastic unloading E alpha dT / (1 - nu)", 1e-9);
            if (ui::propD("Poisson ratio (-)", &S.nu, "%.3f")) changed();
            if (ui::propD("Thermal expansion (1/K)", &S.alphaT, "%.3g", "linear; 0 = no thermal unloading")) changed();
            ui::endProps();
        }

        // ------------------------------------------------------------ EDFM pre-check
        ui::header("Check fractures (EDFM pre-processing)");
        if (ImGui::Button(ICON_SEARCH " Check fractures")) {
            g_check = FracCheck{};
            g_check.version = caseVersion_;
            thermal::DeckModel m; std::string err;
            try {
                if (!thermal::deckModel(case_, m, &err)) g_check.error = err.empty() ? "the thermal model could not be built" : err;
            } catch (const std::exception& e) { g_check.error = e.what(); }
            if (g_check.error.empty()) {
                const auto& FC = m.base.fracCells;
                g_check.fracCells = int(FC.size()); g_check.ffConns = int(m.base.fracFF.size());
                g_check.rows.resize(T.fractures.size());
                for (size_t f = 0; f < T.fractures.size(); ++f) g_check.rows[f].name = T.fractures[f].name;
                std::vector<double> wsum(T.fractures.size(), 0);
                for (const auto& c : FC) {
                    if (c.frac < 0 || size_t(c.frac) >= g_check.rows.size()) continue;
                    auto& r = g_check.rows[size_t(c.frac)];
                    const double a = c.segLen * c.segH;
                    r.segs++; r.area += a;
                    r.snEff += a * (c.sn0 - c.p0); r.tau += a * c.tau0; r.p += a * c.p0;
                    wsum[size_t(c.frac)] += a;
                }
                for (size_t f = 0; f < g_check.rows.size(); ++f) {
                    auto& r = g_check.rows[f];
                    if (wsum[f] > 0) { r.snEff /= wsum[f]; r.tau /= wsum[f]; r.p /= wsum[f]; }
                    r.slip = r.snEff > 0 ? r.tau / r.snEff : 0;
                }
                if (!m.stepWells.empty())
                    for (size_t w = 0; w < m.stepWells[0].size() && w < m.wellNames.size(); ++w)
                        for (int q : m.stepWells[0][w].cells)
                            for (const auto& c : FC)
                                if (c.cell == q && c.frac >= 0 && size_t(c.frac) < g_check.rows.size()) {
                                    auto& L = g_check.rows[size_t(c.frac)].wells;
                                    if (std::find(L.begin(), L.end(), m.wellNames[w]) == L.end()) { L.push_back(m.wellNames[w]); g_check.wellConns++; }
                                }
            }
            g_check.valid = true;
        }
        if (g_check.valid) {
            if (!g_check.error.empty()) ImGui::TextColored(col4(Pal::Err), "%s", g_check.error.c_str());
            else {
                ImGui::Text("%d fracture cells, %d fracture-fracture connections, %d well-fracture links", g_check.fracCells, g_check.ffConns, g_check.wellConns);
                if (g_check.version != caseVersion_) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "(inputs changed - check again)"); }
                const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp;
                if (ui::beginTable("fcheck", 7, tf)) {
                    for (const char* c : {"Fracture", "Cells", "Area (m2)", "Mean p0 (MPa)", "Eff. normal stress (MPa)", "Slip tendency tau/s'", "Connected wells"})
                        ui::tableColumn(c);
                    ImGui::TableHeadersRow();
                    for (const auto& r : g_check.rows) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(r.name.c_str());
                        ImGui::TableNextColumn(); r.segs ? ImGui::Text("%d", r.segs) : ImGui::TextColored(col4(Pal::Err), "0 (outside the grid?)");
                        ImGui::TableNextColumn(); ImGui::Text("%.4g", r.area);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", r.p * 1e-6);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", r.snEff * 1e-6);
                        ImGui::TableNextColumn();
                        if (r.slip > 0.6) ImGui::TextColored(col4(Pal::Warn), "%.3f (critically stressed)", r.slip); else ImGui::Text("%.3f", r.slip);
                        ImGui::TableNextColumn();
                        std::string ws; for (const auto& w : r.wells) ws += (ws.empty() ? "" : ", ") + w;
                        if (ws.empty()) ImGui::TextColored(col4(Pal::Warn), "none"); else ImGui::TextUnformatted(ws.c_str());
                    }
                    ui::endTable();
                }
                ImGui::TextDisabled("Slip tendency > friction (about 0.6): hydro-shear will dilate the plane once the pressure rises.");
            }
        }
    }
    return true;
}

// ============================================================================ 3D fracture panels
std::vector<PanelQuad> App::fracturePanels(const Case& c, const Grid* g, const RunResults* R, int snap, const std::string& prop) {
    std::vector<PanelQuad> out;
    if (!c.thermal.enabled && !c.grid.hfracs.empty()) {
        // C15 hydraulic fractures of black-oil / compositional cases (case length units): one panel per fracture cell of
        // the built grid (coloured by the cell value of a run), else the input planes
        std::string P = prop; for (char& ch : P) ch = char(std::toupper(static_cast<unsigned char>(ch)));
        const Snapshot* sn = nullptr;
        if (R && R->valid && !R->snaps.empty() && R->grid.get() == g) sn = &R->snaps[size_t(std::clamp(snap, 0, int(R->snaps.size()) - 1))];
        const std::vector<float>* arr = !sn ? nullptr : P == "PRESSURE" ? &sn->pressure : P == "SWAT" ? &sn->sw : P == "SGAS" ? &sn->sg : P == "SOIL" ? &sn->so : nullptr;
        if (sn && !arr) for (const auto& e : sn->extra) if (e.first == prop || e.first == P) { arr = &e.second; break; }   // C30: any engine array
        // C30: matrix-only arrays (stress shadow, phase-field damage, geomechanics ...) have no fracture-cell values - draw the
        // fracture as a translucent outline so that the field around it stays visible
        bool fracValues = false;
        if (arr && g)
            for (const Grid::HfCell& H : g->hfCells) {
                const int a = size_t(H.global) < g->activeOf.size() ? g->activeOf[size_t(H.global)] : -1;
                if (a >= 0 && size_t(a) < arr->size() && (*arr)[size_t(a)] != 0.0f) { fracValues = true; break; }
            }
        if (!fracValues) arr = nullptr;
        if (g && g->numHfCells() > 0) {
            std::vector<bool> named(c.grid.hfracs.size(), false);
            for (const Grid::HfCell& H : g->hfCells) {
                if (H.frac < 0 || size_t(H.frac) >= c.grid.hfracs.size() || size_t(H.global) >= g->center.size()) continue;
                const GridInput::HFrac& F = c.grid.hfracs[size_t(H.frac)];
                double ux = F.x1 - F.x0, uy = F.y1 - F.y0; const double ul = std::hypot(ux, uy);
                if (!(ul > 0)) continue;
                ux /= ul; uy /= ul;
                const Vec3 m = g->center[size_t(H.global)];
                const float hx = float(0.5 * H.len * ux), hy = float(0.5 * H.len * uy), hz = float(0.5 * H.h);
                PanelQuad pq;
                pq.p[0] = {m.x - hx, m.y - hy, m.z - hz}; pq.p[1] = {m.x + hx, m.y + hy, m.z - hz};
                pq.p[2] = {m.x + hx, m.y + hy, m.z + hz}; pq.p[3] = {m.x - hx, m.y - hy, m.z + hz};
                if (arr) { const int a = g->activeOf[size_t(H.global)]; if (a >= 0 && size_t(a) < arr->size()) pq.value = (*arr)[size_t(a)]; }
                pq.rgba[3] = (sn && !arr) ? 0.25f : 0.9f;
                if (!named[size_t(H.frac)]) { named[size_t(H.frac)] = true; pq.label = F.name; }
                out.push_back(std::move(pq));
            }
            if (!out.empty()) return out;
        }
        const float top = g ? float(g->bbMin.z) : 0.0f, bot = g ? float(g->bbMax.z) : 1.0f;
        for (const GridInput::HFrac& F : c.grid.hfracs) {
            const float zt = F.zTop < 0 ? top : float(F.zTop), zb = F.zBot < 0 ? bot : float(F.zBot);
            PanelQuad pq;
            pq.p[0] = {float(F.x0), float(F.y0), zt}; pq.p[1] = {float(F.x1), float(F.y1), zt};
            pq.p[2] = {float(F.x1), float(F.y1), zb}; pq.p[3] = {float(F.x0), float(F.y0), zb};
            pq.label = F.name;
            out.push_back(std::move(pq));
        }
        return out;
    }
    if (!c.thermal.enabled || c.thermal.fractures.empty()) return out;
    const bool fld = c.units == UnitSystem::Field;
    const double lenU = fld ? 0.3048 : 1.0;
    std::string P = prop; for (char& ch : P) ch = char(std::toupper(static_cast<unsigned char>(ch)));
    if (R && !R->fracSegs.empty() && !R->fracSnap.empty()) {
        const int s = std::max(0, std::min(snap, int(R->fracSnap.size()) - 1));
        const RunResults::FracState& st = R->fracSnap[size_t(s)];
        const RunResults::FracState& s0 = R->fracSnap[0];
        std::vector<bool> named(R->fracNames.size() + 1, false);
        for (size_t i = 0; i < R->fracSegs.size(); ++i) {
            const RunResults::FracSeg& q = R->fracSegs[i];
            PanelQuad pq;
            pq.p[0] = {q.xa, q.ya, q.zt}; pq.p[1] = {q.xb, q.yb, q.zt}; pq.p[2] = {q.xb, q.yb, q.zb}; pq.p[3] = {q.xa, q.ya, q.zb};
            auto at = [&](const std::vector<float>& v) { return i < v.size() ? v[i] : NAN; };
            if (P == "FRACAP") pq.value = at(st.ap);
            else if (P == "TEMP" || P == "FRACT") pq.value = at(st.T);
            else if (P == "PRESSURE" || P == "FRACP") pq.value = at(st.p);
            else if (P == "FRACK" || P == "FRACCOND") pq.value = at(st.cond);
            else if (P == "DTEMP") pq.value = at(st.T) - at(s0.T);
            pq.rgba[3] = 0.9f;
            if (q.frac >= 0 && size_t(q.frac) < R->fracNames.size() && !named[size_t(q.frac)]) { named[size_t(q.frac)] = true; pq.label = R->fracNames[size_t(q.frac)]; }
            out.push_back(std::move(pq));
        }
        return out;
    }
    // input planes (model coordinates in case units)
    const float top = g ? float(g->bbMin.z) : 0.0f, bot = g ? float(g->bbMax.z) : 1.0f;
    for (const ThermalFracture& F : c.thermal.fractures) {
        const float zt = F.zTop < 0 ? top : float(F.zTop / lenU), zb = F.zBot < 0 ? bot : float(F.zBot / lenU);
        const float x0 = float(F.x0 / lenU), y0 = float(F.y0 / lenU), x1 = float(F.x1 / lenU), y1 = float(F.y1 / lenU);
        PanelQuad pq;
        pq.p[0] = {x0, y0, zt}; pq.p[1] = {x1, y1, zt}; pq.p[2] = {x1, y1, zb}; pq.p[3] = {x0, y0, zb};
        pq.label = F.name;
        out.push_back(std::move(pq));
    }
    return out;
}

// ============================================================================ results: Geothermal
void App::geothermalTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor) {
    const float fs = ImGui::GetFontSize();
    const bool fld = R.input.units == UnitSystem::Field;
    const char* tU = fld ? "\xC2\xB0""F" : "\xC2\xB0""C";
    const char* pU = fld ? "psia" : "bar";
    const size_t N = std::min({R.summary.size(), R.wellT.size(), R.wellMW.size()});
    if (N < 2) { ui::hint("The geothermal analysis appears after the first report steps of a thermal run."); return; }
    const size_t nw = R.wellNames.size();
    std::vector<double> t(N);
    for (size_t k = 0; k < N; ++k) t[k] = R.summary[k].time * xs;
    const double T0 = R.summary[0].favgT;          // initial (pore-volume average) reservoir temperature
    const double Tin = R.tInjRef;
    const double massPerRate = fld ? 0.158987294928 * 1000.0 / 86400.0 : 1000.0 / 86400.0;   // surface water rate -> kg/s
    static float eff = 0.12f, lifeDrop = 0.10f;
    static int sub = 0;
    auto isProd = [&](size_t w) { return w < R.wellTypes.size() && R.wellTypes[w] == WellType::Producer; };

    // ---- derived time series
    std::vector<std::vector<double>> prodT(nw), TD(nw);
    std::vector<double> mwP(N, 0), mwI(N, 0), mwNet(N, 0), mwE(N, 0), cumGWh(N, 0), rec(N, 0), imp(N, 0), cool(N, 0), qprod(N, 0);
    std::vector<std::vector<double>> II(nw), PI(nw);
    for (size_t w = 0; w < nw; ++w) { prodT[w].assign(N, NAN); TD[w].assign(N, NAN); II[w].assign(N, 0); PI[w].assign(N, 0); }
    double lifeT = -1;
    for (size_t k = 0; k < N; ++k) {
        double bhpI = 0, nI = 0, bhpP = 0, nP = 0, qP = 0;
        const auto& SW = R.summary[k].wells;
        for (size_t w = 0; w < nw; ++w) {
            const double Tw = w < R.wellT[k].size() ? R.wellT[k][w] : 0.0;
            const double H = w < R.wellMW[k].size() ? std::fabs(R.wellMW[k][w]) : 0.0;
            const double fpr = R.summary[k].fprp > 0 ? R.summary[k].fprp : R.summary[k].fpr;
            if (isProd(w)) {
                mwP[k] += H;
                if (w < SW.size() && SW[w].water > 0 && Tw != 0) {
                    prodT[w][k] = Tw;
                    if (std::fabs(T0 - Tin) > 1e-9) TD[w][k] = (T0 - Tw) / (T0 - Tin);
                    const double q = SW[w].water * massPerRate;
                    qP += q; bhpP += SW[w].bhp; nP += 1;
                    PI[w][k] = fpr - SW[w].bhp > 1e-9 ? q / (fpr - SW[w].bhp) : 0;
                }
            } else {
                mwI[k] += H;
                if (w < SW.size() && SW[w].waterInj > 0) {
                    const double q = SW[w].waterInj * massPerRate;
                    bhpI += SW[w].bhp; nI += 1;
                    II[w][k] = SW[w].bhp - fpr > 1e-9 ? q / (SW[w].bhp - fpr) : 0;
                }
            }
        }
        mwNet[k] = mwP[k] - mwI[k];
        mwE[k] = std::max(0.0, mwNet[k]) * eff;
        qprod[k] = qP;
        imp[k] = (nI > 0 && nP > 0 && qP > 0) ? (bhpI / nI - bhpP / nP) / qP : 0;
        if (k > 0) cumGWh[k] = cumGWh[k - 1] + 0.5 * (mwNet[k] + mwNet[k - 1]) * (R.summary[k].time - R.summary[k - 1].time) * 24.0 / 1000.0;
        rec[k] = R.heatInPlace > 0 ? cumGWh[k] * 3.6e12 / R.heatInPlace * 100.0 : 0;
        cool[k] = k < R.cooledFrac.size() ? R.cooledFrac[k] * 100.0 : 0;
        if (lifeT < 0) for (size_t w = 0; w < nw; ++w) if (std::isfinite(TD[w][k]) && TD[w][k] >= lifeDrop && k > 0) { lifeT = R.summary[k].time; break; }
    }
    const int ks = std::max(0, std::min(snap_, int(N) - 1));

    // ---- KPI strip
    {
        double tp = 0, ntp = 0, td = 0;
        for (size_t w = 0; w < nw; ++w) if (std::isfinite(prodT[w][size_t(ks)])) { tp += prodT[w][size_t(ks)]; td = std::max(td, TD[w][size_t(ks)]); ntp += 1; }
        if (ntp > 0) tp /= ntp;
        ImGui::SetNextItemWidth(fs * 6); ImGui::SliderFloat("Conversion efficiency", &eff, 0.02f, 0.25f, "%.2f");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("net electric / net thermal power (binary ORC at 150-200 \xC2\xB0""C: about 0.10-0.14)");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6); ImGui::SliderFloat("Lifetime threshold T_D", &lifeDrop, 0.02f, 0.5f, "%.2f");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("end of life when the normalised drawdown (T0 - Tp)/(T0 - Tinj) of a producer reaches this value");
        ImGui::SameLine(); ImGui::TextDisabled("  t = %.1f d (%.2f yr)", R.summary[size_t(ks)].time, R.summary[size_t(ks)].time / 365.25);
        const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchSame;
        if (ImGui::BeginTable("egskpi", 8, tf)) {
            char b[64];
            auto kpi = [&](const char* name, const char* val, ImU32 col) {
                ImGui::TableNextColumn();
                ImGui::TextDisabled("%s", name);
                ImGui::PushStyleColor(ImGuiCol_Text, col); ImGui::TextUnformatted(val); ImGui::PopStyleColor();
            };
            std::snprintf(b, sizeof(b), ntp > 0 ? "%.1f %s" : "-", tp, tU); kpi("Production T", b, IM_COL32(255, 140, 80, 255));
            std::snprintf(b, sizeof(b), "%.3f", td); kpi("Drawdown T_D", b, td >= lifeDrop ? IM_COL32(240, 90, 80, 255) : IM_COL32(200, 210, 220, 255));
            std::snprintf(b, sizeof(b), "%.2f MWth", mwNet[size_t(ks)]); kpi("Net thermal power", b, IM_COL32(255, 200, 90, 255));
            std::snprintf(b, sizeof(b), "%.2f MWe", mwE[size_t(ks)]); kpi("Electric (estimate)", b, IM_COL32(120, 220, 140, 255));
            std::snprintf(b, sizeof(b), "%.1f GWh", cumGWh[size_t(ks)]); kpi("Cumulative heat", b, IM_COL32(200, 210, 220, 255));
            std::snprintf(b, sizeof(b), "%.2f %%", rec[size_t(ks)]); kpi("Heat recovery", b, IM_COL32(200, 210, 220, 255));
            std::snprintf(b, sizeof(b), imp[size_t(ks)] > 0 ? "%.3g %s/(kg/s)" : "-", imp[size_t(ks)], pU); kpi("Impedance", b, IM_COL32(150, 190, 255, 255));
            if (lifeT > 0) std::snprintf(b, sizeof(b), "%.1f yr", lifeT / 365.25); else std::snprintf(b, sizeof(b), "> %.1f yr", R.summary[N - 1].time / 365.25);
            kpi("Thermal lifetime", b, lifeT > 0 ? IM_COL32(240, 90, 80, 255) : IM_COL32(120, 220, 140, 255));
            ImGui::EndTable();
        }
    }

    if (!ImGui::BeginTabBar("geosub")) return;
    // ------------------------------------------------------------------ performance plots
    if (ui::tabItem("Performance")) {
        sub = 0;
        const bool hasFr = !R.fracSegs.empty() && R.fracSnap.size() >= N;
        const int nPlots = hasFr ? 8 : 5;
        const ImVec2 av = ImGui::GetContentRegionAvail();
        const PlotGrid G = fitGrid(nPlots, av, fs * 22, fs * 12, 1.6f);
        if (G.scroll) ImGui::BeginChild("##geoscroll", av);
        fig::recorder().beginFigure("run.geo", G.rows, G.cols);
        const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
        if (ImPlot::BeginSubplots("##geo", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
            char yl[48];
            std::snprintf(yl, sizeof(yl), "Temperature (%s)", tU);
            if (ImPlot::BeginPlot("Production / injection temperature")) {
                recPanel("Production / injection temperature", xl, yl, "");
                ImPlot::SetupAxes(xl, yl, af, af);
                int ci = 0;
                for (size_t w = 0; w < nw; ++w) {
                    std::vector<double> y(N, NAN);
                    for (size_t k = 0; k < N; ++k) {
                        const double Tw = w < R.wellT[k].size() ? R.wellT[k][w] : 0;
                        if (Tw != 0) y[k] = Tw;
                    }
                    std::vector<double> xx, yy;
                    for (size_t k = 0; k < N; ++k) if (std::isfinite(y[k])) { xx.push_back(t[k]); yy.push_back(y[k]); }
                    plotXYe(R.wellNames[w].c_str(), xx, yy, isProd(w) ? seriesColor(ci) : IM_COL32(90, 150, 255, 255), isProd(w) ? 2.2f : 1.4f);
                    ++ci;
                }
                hline("##T0", T0, IM_COL32(255, 255, 255, 90));
                ImPlot::TagY(T0, ImVec4(0.5f, 0.5f, 0.5f, 1), "T0");
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Thermal drawdown (normalised)")) {
                recPanel("Thermal drawdown (normalised)", xl, "T_D = (T0 - Tp)/(T0 - Tinj)", "");
                ImPlot::SetupAxes(xl, "T_D (-)", af, af);
                int ci = 0;
                for (size_t w = 0; w < nw; ++w) {
                    if (!isProd(w)) continue;
                    std::vector<double> xx, yy;
                    for (size_t k = 0; k < N; ++k) if (std::isfinite(TD[w][k])) { xx.push_back(t[k]); yy.push_back(TD[w][k]); }
                    plotXYe(R.wellNames[w].c_str(), xx, yy, seriesColor(ci++), 2.2f);
                }
                hline("##life", lifeDrop, IM_COL32(240, 90, 80, 200));
                ImPlot::TagY(lifeDrop, ImVec4(0.9f, 0.35f, 0.3f, 1), "limit");
                if (lifeT > 0) { const double lx = lifeT * xs; ImPlotSpec sp; sp.LineColor = ImVec4(0.9f, 0.35f, 0.3f, 0.8f); sp.Flags = ImPlotItemFlags_NoLegend; ImPlot::PlotInfLines("##lt", &lx, 1, sp); }
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Thermal / electric power")) {
                recPanel("Thermal / electric power", xl, "MW", "");
                ImPlot::SetupAxes(xl, "MW", af, af);
                plotXYe("Produced (MWth)", t, mwP, IM_COL32(255, 140, 80, 255), 1.6f);
                plotXYe("Injected (MWth)", t, mwI, IM_COL32(90, 150, 255, 255), 1.4f);
                plotXYe("Net thermal (MWth)", t, mwNet, IM_COL32(255, 200, 90, 255), 2.6f);
                plotXYe("Electric estimate (MWe)", t, mwE, IM_COL32(120, 220, 140, 255), 2.0f);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Cumulative heat and heat recovery")) {
                recPanel("Cumulative heat and heat recovery", xl, "GWh (thermal)", "Recovery (%)");
                ImPlot::SetupAxes(xl, "GWh_th", af, af);
                ImPlot::SetupAxis(ImAxis_Y2, "Recovery of the heat in place (%)", af | ImPlotAxisFlags_AuxDefault);
                plotXYe("Cumulative net heat", t, cumGWh, IM_COL32(255, 200, 90, 255), 2.4f);
                useY2(true); plotXYe("Heat recovery factor", t, rec, IM_COL32(200, 120, 255, 255), 1.6f); useY2(false);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Injectivity / productivity / impedance")) {
                char l1[48], l2[48];
                std::snprintf(l1, sizeof(l1), "kg/s per %s", pU); std::snprintf(l2, sizeof(l2), "%s per kg/s", pU);
                recPanel("Injectivity / productivity / impedance", xl, l1, l2);
                ImPlot::SetupAxes(xl, l1, af, af);
                ImPlot::SetupAxis(ImAxis_Y2, l2, af | ImPlotAxisFlags_AuxDefault);
                int ci = 0;
                for (size_t w = 0; w < nw; ++w) {
                    std::vector<double> xx, yy;
                    const auto& V = isProd(w) ? PI[w] : II[w];
                    for (size_t k = 1; k < N; ++k) if (V[k] > 0) { xx.push_back(t[k]); yy.push_back(V[k]); }
                    plotXYe(((isProd(w) ? "PI " : "II ") + R.wellNames[w]).c_str(), xx, yy, seriesColor(ci++), 1.6f);
                }
                useY2(true);
                { std::vector<double> xx, yy; for (size_t k = 1; k < N; ++k) if (imp[k] > 0) { xx.push_back(t[k]); yy.push_back(imp[k]); }
                  plotXYe("Impedance (dBHP / q_prod)", xx, yy, IM_COL32(255, 255, 255, 220), 2.2f); }
                useY2(false);
                cursor(); recEnd(); ImPlot::EndPlot();
            }
            if (hasFr) {
                const size_t nf = R.fracNames.size();
                // per fracture mean / min / max over its segments
                auto perFrac = [&](int which, std::vector<std::vector<double>>& mean, std::vector<std::vector<double>>& lo, std::vector<std::vector<double>>& hi) {
                    mean.assign(nf, std::vector<double>(N, 0)); lo.assign(nf, std::vector<double>(N, 1e300)); hi.assign(nf, std::vector<double>(N, -1e300));
                    std::vector<double> cnt(nf, 0);
                    for (size_t k = 0; k < N; ++k) {
                        std::fill(cnt.begin(), cnt.end(), 0.0);
                        const auto& st = R.fracSnap[k];
                        const std::vector<float>& v = which == 0 ? st.ap : which == 1 ? st.T : st.cond;
                        for (size_t i = 0; i < R.fracSegs.size() && i < v.size(); ++i) {
                            const int f = R.fracSegs[i].frac;
                            if (f < 0 || size_t(f) >= nf) continue;
                            mean[size_t(f)][k] += v[i]; cnt[size_t(f)] += 1;
                            lo[size_t(f)][k] = std::min(lo[size_t(f)][k], double(v[i])); hi[size_t(f)][k] = std::max(hi[size_t(f)][k], double(v[i]));
                        }
                        for (size_t f = 0; f < nf; ++f) { if (cnt[f] > 0) mean[f][k] /= cnt[f]; else lo[f][k] = hi[f][k] = 0; }
                    }
                };
                std::vector<std::vector<double>> m, lo, hi;
                if (ImPlot::BeginPlot("Fracture aperture (mean, min-max band)")) {
                    recPanel("Fracture aperture", xl, "Aperture (mm)", "");
                    ImPlot::SetupAxes(xl, "Aperture (mm)", af, af);
                    perFrac(0, m, lo, hi);
                    for (size_t f = 0; f < nf; ++f) {
                        ImPlotSpec sp; sp.FillColor = col4(seriesColor(int(f))); sp.FillAlpha = 0.18f; sp.Flags = ImPlotItemFlags_NoLegend;
                        ImPlot::PlotShaded(("##band" + std::to_string(f)).c_str(), t.data(), lo[f].data(), hi[f].data(), int(N), sp);
                        plotXYe(R.fracNames[f].c_str(), t, m[f], seriesColor(int(f)), 2.0f);
                    }
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Fracture temperature (mean)")) {
                    recPanel("Fracture temperature", xl, yl, "");
                    ImPlot::SetupAxes(xl, yl, af, af);
                    perFrac(1, m, lo, hi);
                    for (size_t f = 0; f < nf; ++f) {
                        ImPlotSpec sp; sp.FillColor = col4(seriesColor(int(f))); sp.FillAlpha = 0.15f; sp.Flags = ImPlotItemFlags_NoLegend;
                        ImPlot::PlotShaded(("##tb" + std::to_string(f)).c_str(), t.data(), lo[f].data(), hi[f].data(), int(N), sp);
                        plotXYe(R.fracNames[f].c_str(), t, m[f], seriesColor(int(f)), 2.0f);
                    }
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Fracture conductivity and cooled rock")) {
                    recPanel("Fracture conductivity and cooled rock", xl, "Conductivity k_f w (mD m)", "Cooled bulk volume (%)");
                    ImPlot::SetupAxes(xl, "k_f w (mD m)", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, "Matrix cooled > 10 K (%)", af | ImPlotAxisFlags_AuxDefault);
                    perFrac(2, m, lo, hi);
                    for (size_t f = 0; f < nf; ++f) plotXYe(R.fracNames[f].c_str(), t, m[f], seriesColor(int(f)), 1.8f);
                    useY2(true); plotXYe("Cooled rock volume", t, cool, IM_COL32(90, 180, 255, 255), 2.4f); useY2(false);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        if (G.scroll) ImGui::EndChild();
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ fracture faces / profiles
    if (!R.fracSegs.empty() && R.fracSnap.size() > size_t(ks) && ui::tabItem("Fractures")) {
        sub = 1;
        const size_t nf = R.fracNames.size();
        static int fsel = 0, qty = 1;
        fsel = std::max(0, std::min(fsel, int(nf) - 1));
        ImGui::SetNextItemWidth(fs * 10);
        if (ImGui::BeginCombo("Fracture", nf ? R.fracNames[size_t(fsel)].c_str() : "-")) {
            for (int f = 0; f < int(nf); ++f) if (ImGui::Selectable(R.fracNames[size_t(f)].c_str(), f == fsel)) fsel = f;
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        static const char* qn[] = {"Aperture (mm)", "Temperature", "Pressure", "Conductivity (mD m)", "Cooling T - T0"};
        ImGui::SetNextItemWidth(fs * 12); ImGui::Combo("Face map", &qty, qn, 5);
        ImGui::SameLine(); ImGui::TextDisabled("time follows the 3D time cursor (double-click a performance plot or use the time slider)");
        const auto& st = R.fracSnap[size_t(ks)];
        const auto& s0 = R.fracSnap[0];
        // segments of the fracture: distance along the trace from the first end point
        const double lenU = fld ? 0.3048 : 1.0;
        double ox = 0, oy = 0, ux = 1, uy = 0;
        if (size_t(fsel) < R.input.thermal.fractures.size()) {
            const auto& F = R.input.thermal.fractures[size_t(fsel)];
            ox = F.x0 / lenU; oy = F.y0 / lenU;
            const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
            if (L > 0) { ux = (F.x1 - F.x0) / L; uy = (F.y1 - F.y0) / L; }
        }
        struct Seg { double sa, sb, zt, zb; float v, ap, T, p, c; };
        std::vector<Seg> segs;
        double vlo = 1e300, vhi = -1e300;
        for (size_t i = 0; i < R.fracSegs.size(); ++i) {
            const auto& q = R.fracSegs[i];
            if (q.frac != fsel) continue;
            Seg g;
            g.sa = (q.xa - ox) * ux + (q.ya - oy) * uy; g.sb = (q.xb - ox) * ux + (q.yb - oy) * uy;
            if (g.sa > g.sb) std::swap(g.sa, g.sb);
            g.zt = q.zt; g.zb = q.zb;
            g.ap = i < st.ap.size() ? st.ap[i] : 0; g.T = i < st.T.size() ? st.T[i] : 0; g.p = i < st.p.size() ? st.p[i] : 0; g.c = i < st.cond.size() ? st.cond[i] : 0;
            const float T0s = i < s0.T.size() ? s0.T[i] : g.T;
            g.v = qty == 0 ? g.ap : qty == 1 ? g.T : qty == 2 ? g.p : qty == 3 ? g.c : g.T - T0s;
            vlo = std::min(vlo, double(g.v)); vhi = std::max(vhi, double(g.v));
            segs.push_back(g);
        }
        if (segs.empty()) ui::hint("This fracture has no cells inside the grid.");
        else {
            if (vhi - vlo < 1e-9) { vhi = vlo + 1; }
            const ImVec2 av = ImGui::GetContentRegionAvail();
            const float h = std::max(fs * 14, av.y * 0.55f - fs * 2);
            colorBar(qn[qty], vlo, vhi, fs * 14);
            if (ImPlot::BeginPlot("##face", ImVec2(av.x * 0.55f, h), ImPlotFlags_NoMenus)) {
                ImPlot::SetupAxes(fld ? "distance along the trace (ft)" : "distance along the trace (m)", fld ? "depth (ft)" : "depth (m)", ImPlotAxisFlags_AutoFit,
                                  ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                double smin = 1e300, smax = -1e300, zmin = 1e300, zmax = -1e300;
                for (const auto& g : segs) { smin = std::min(smin, g.sa); smax = std::max(smax, g.sb); zmin = std::min(zmin, g.zt); zmax = std::max(zmax, g.zb); }
                const double bx[2] = {smin, smax}, by[2] = {zmin, zmax};
                ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
                ImPlot::PlotScatter("##bb", bx, by, 2, sp);
                ImPlot::PushPlotClipRect();
                ImDrawList* dl = ImPlot::GetPlotDrawList();
                (void)dl;
                for (const auto& g : segs) pk::rect(g.sa, g.zt, g.sb, g.zb, cmapU32(float((g.v - vlo) / (vhi - vlo))), IM_COL32(0, 0, 0, 60), 1.0f);
                if (pk::recording()) {   // C23: colour scale of the face map in the export
                    ImU32 lut[64];
                    for (int k = 0; k < 64; ++k) lut[k] = cmapU32(k / 63.0f);
                    pk::colorScale(qn[qty], vlo, vhi, lut, 64);
                }
                ImPlot::PopPlotClipRect();
                if (ImPlot::IsPlotHovered()) {
                    const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                    for (const auto& g : segs)
                        if (mp.x >= g.sa && mp.x <= g.sb && mp.y >= g.zt && mp.y <= g.zb) {
                            ImGui::SetTooltip("s = %.1f - %.1f, depth %.1f - %.1f\naperture %.4g mm\nT %.2f %s\np %.3f %s\nk_f w %.4g mD m", g.sa, g.sb, g.zt, g.zb,
                                              g.ap, g.T, tU, g.p, pU, g.c);
                            break;
                        }
                }
                ImPlot::EndPlot();
            }
            ImGui::SameLine();
            // along-trace profile (height-weighted average over the layers)
            std::map<long long, std::array<double, 5>> col;   // key: rounded s-centre
            for (const auto& g : segs) {
                const long long key = (long long)std::llround(0.5 * (g.sa + g.sb) * 100.0);
                auto& a = col[key];
                const double wgt = g.zb - g.zt;
                a[0] += wgt; a[1] += wgt * g.ap; a[2] += wgt * g.T; a[3] += wgt * g.p; a[4] += wgt * g.c;
            }
            std::vector<double> s, ap, Tt, pp;
            for (const auto& kv : col) { const auto& a = kv.second; if (a[0] <= 0) continue; s.push_back(kv.first / 100.0); ap.push_back(a[1] / a[0]); Tt.push_back(a[2] / a[0]); pp.push_back(a[3] / a[0]); }
            if (ImPlot::BeginPlot("Along-fracture profile", ImVec2(-1, h + fs * 1.4f), ImPlotFlags_NoMenus)) {
                char yl[48]; std::snprintf(yl, sizeof(yl), "Temperature (%s)", tU);
                ImPlot::SetupAxes(fld ? "distance (ft)" : "distance (m)", yl, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y2, "Aperture (mm)", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                plotXYe("Temperature", s, Tt, IM_COL32(255, 140, 80, 255), 2.2f, true);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                plotXYe("Aperture", s, ap, IM_COL32(120, 220, 140, 255), 2.0f, true);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
                ImPlot::EndPlot();
            }
        }
        // per-fracture summary table at the current time
        const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY;
        if (ui::beginTable("fracsum", 8, tf, ImVec2(0, std::max(fs * 6, ImGui::GetContentRegionAvail().y)))) {
            ImGui::TableSetupScrollFreeze(1, 1);
            char h1[48], h2[48];
            std::snprintf(h1, sizeof(h1), "Mean T (%s)", tU); std::snprintf(h2, sizeof(h2), "Mean p (%s)", pU);
            const char* hdr[] = {"Fracture", "Cells", "Area", "Aperture t=0 (mm)", "Aperture now (mm)", h1, h2, "k_f w now (mD m)"};
            for (const char* c : hdr) ui::tableColumn(c);
            ImGui::TableHeadersRow();
            for (size_t f = 0; f < nf; ++f) {
                int n = 0; double area = 0, a0 = 0, a1 = 0, Tm = 0, pm = 0, cm = 0;
                for (size_t i = 0; i < R.fracSegs.size(); ++i) {
                    const auto& q = R.fracSegs[i];
                    if (q.frac != int(f)) continue;
                    const double A = std::hypot(double(q.xb - q.xa), double(q.yb - q.ya)) * double(q.zb - q.zt);
                    ++n; area += A;
                    a0 += A * (i < s0.ap.size() ? s0.ap[i] : 0); a1 += A * (i < st.ap.size() ? st.ap[i] : 0);
                    Tm += A * (i < st.T.size() ? st.T[i] : 0); pm += A * (i < st.p.size() ? st.p[i] : 0); cm += A * (i < st.cond.size() ? st.cond[i] : 0);
                }
                if (area > 0) { a0 /= area; a1 /= area; Tm /= area; pm /= area; cm /= area; }
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(R.fracNames[f].c_str());
                ImGui::TableNextColumn(); ImGui::Text("%d", n);
                ImGui::TableNextColumn(); ImGui::Text("%.4g %s", area, fld ? "ft2" : "m2");
                ImGui::TableNextColumn(); ImGui::Text("%.4g", a0);
                ImGui::TableNextColumn(); ImGui::Text("%.4g", a1);
                ImGui::TableNextColumn(); ImGui::Text("%.2f", Tm);
                ImGui::TableNextColumn(); ImGui::Text("%.3f", pm);
                ImGui::TableNextColumn(); ImGui::Text("%.4g", cm);
            }
            ui::endTable();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ plan-view cold plume map
    if (R.grid && ks < int(R.snaps.size()) && ui::tabItem("Cooling map")) {
        sub = 2;
        const Grid& g = *R.grid;
        static int layer = -1, mapQ = 0;
        if (layer < 0 || layer >= g.nz) layer = g.nz / 2;
        ImGui::SetNextItemWidth(fs * 10); ImGui::SliderInt("Layer (K)", &layer, 0, std::max(0, g.nz - 1), "%d");
        ImGui::SameLine();
        static const char* mq[] = {"Cooling T - T0", "Temperature", "Pressure"};
        ImGui::SetNextItemWidth(fs * 10); ImGui::Combo("Map", &mapQ, mq, 3);
        const Snapshot& sn = R.snaps[size_t(ks)];
        std::vector<float> v;
        if (mapQ == 2) v = sn.pressure; else v = sn.field(mapQ == 0 ? "DTEMP" : "TEMP");
        double lo = 1e300, hi = -1e300;
        for (int i = 0; i < g.nx; ++i)
            for (int j = 0; j < g.ny; ++j) {
                const int q = i + g.nx * (j + g.ny * layer);
                const int a = size_t(q) < g.activeOf.size() ? g.activeOf[size_t(q)] : -1;
                if (a < 0 || size_t(a) >= v.size() || !std::isfinite(v[size_t(a)])) continue;
                lo = std::min(lo, double(v[size_t(a)])); hi = std::max(hi, double(v[size_t(a)]));
            }
        if (hi < lo) { ui::hint("No values on this layer."); ImGui::EndTabItem(); ImGui::EndTabBar(); return; }
        if (hi - lo < 1e-9) hi = lo + 1;
        colorBar(mq[mapQ], lo, hi, fs * 14);
        ImGui::SameLine(); ImGui::TextDisabled("  t = %.1f d; white: fracture traces in this layer; dots: well completions", sn.time);
        if (ImPlot::BeginPlot("##plume", ImVec2(-1, -1), ImPlotFlags_Equal | ImPlotFlags_NoMenus)) {
            ImPlot::SetupAxes(fld ? "x (ft)" : "x (m)", fld ? "y (ft)" : "y (m)", 0, ImPlotAxisFlags_Invert);
            ImPlot::SetupAxesLimits(g.bbMin.x, g.bbMax.x, g.bbMin.y, g.bbMax.y, ImPlotCond_Once);
            const double bx[2] = {g.bbMin.x, g.bbMax.x}, by[2] = {g.bbMin.y, g.bbMax.y};
            ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
            ImPlot::PlotScatter("##bb", bx, by, 2, sp);
            ImPlot::PushPlotClipRect();
            
            int hovered = -1;
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            for (int i = 0; i < g.nx; ++i)
                for (int j = 0; j < g.ny; ++j) {
                    const int q = i + g.nx * (j + g.ny * layer);
                    if (size_t(q) >= g.corners.size()) continue;
                    const auto& C = g.corners[size_t(q)];
                    const int a = size_t(q) < g.activeOf.size() ? g.activeOf[size_t(q)] : -1;
                    const ImVec2 p0 = ImPlot::PlotToPixels(C[0].x, C[0].y), p1 = ImPlot::PlotToPixels(C[1].x, C[1].y);
                    const ImVec2 p3 = ImPlot::PlotToPixels(C[3].x, C[3].y), p2 = ImPlot::PlotToPixels(C[2].x, C[2].y);
                    ImU32 colr = IM_COL32(90, 90, 95, 255);
                    if (a >= 0 && size_t(a) < v.size() && std::isfinite(v[size_t(a)])) colr = cmapU32(float((v[size_t(a)] - lo) / (hi - lo)));
                    const double qx[4] = {C[0].x, C[1].x, C[3].x, C[2].x}, qy[4] = {C[0].y, C[1].y, C[3].y, C[2].y};
                    pk::quad(qx, qy, colr);
                    (void)p0; (void)p1; (void)p2; (void)p3;
                    if (mp.x >= C[0].x && mp.x <= C[1].x && mp.y >= C[0].y && mp.y <= C[2].y) hovered = q;
                }
            if (pk::recording()) {   // C23: colour scale of the map in the export
                ImU32 lut[64];
                for (int k = 0; k < 64; ++k) lut[k] = cmapU32(k / 63.0f);
                pk::colorScale(mq[mapQ], lo, hi, lut, 64);
            }
            for (const auto& s : R.fracSegs) {
                const int k = s.host >= 0 ? s.host / (g.nx * g.ny) : -1;
                if (k != layer) continue;
                pk::line(s.xa, s.ya, s.xb, s.yb, IM_COL32(255, 255, 255, 230), 2.5f, "fracture");
            }
            for (const WellInput& w : R.input.wells)
                for (const Completion& c : w.completions) {
                    if (c.i < 0 || c.j < 0 || c.i >= g.nx || c.j >= g.ny) continue;
                    const int q = c.i + g.nx * (c.j + g.ny * std::max(0, std::min(c.k1, g.nz - 1)));
                    if (size_t(q) >= g.center.size()) continue;
                    pk::marker(g.center[size_t(q)].x, g.center[size_t(q)].y, fs * 0.3f,
                               w.type == WellType::Producer ? IM_COL32(60, 200, 110, 255) : IM_COL32(70, 140, 255, 255), 0, 0, 1.0f,
                               w.type == WellType::Producer ? "producer" : "injector");
                }
            ImPlot::PopPlotClipRect();
            if (ImPlot::IsPlotHovered() && hovered >= 0) {
                const int a = g.activeOf[size_t(hovered)];
                const int i = hovered % g.nx, j = (hovered / g.nx) % g.ny;
                if (a >= 0 && size_t(a) < v.size()) ImGui::SetTooltip("cell (%d, %d, %d)\n%s = %.4g", i + 1, j + 1, layer + 1, mq[mapQ], v[size_t(a)]);
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ connections / summary
    if (ui::tabItem("Connections & summary")) {
        sub = 3;
        ui::header("Well - fracture connections (EDFM)");
        if (R.wellFracLinks.empty()) ui::hint("No well is connected to a fracture (no completion in a cell cut by a fracture).");
        else if (ui::beginTable("wfl", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("Well"); ui::tableColumn("Fracture");
            ImGui::TableHeadersRow();
            for (const auto& L : R.wellFracLinks) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(L.first.c_str());
                ImGui::TableNextColumn(); ImGui::TextUnformatted(L.second.c_str());
            }
            ui::endTable();
        }
        ui::header("Run summary");
        if (ui::beginProps("geosum")) {
            char b[96];
            std::snprintf(b, sizeof(b), "%.2f %s", T0, tU); ui::propText("Initial reservoir temperature (PV average)", b);
            std::snprintf(b, sizeof(b), "%.2f %s", Tin, tU); ui::propText("Injection (reference) temperature", b);
            std::snprintf(b, sizeof(b), "%.4g PJ (%.4g GWh)", R.heatInPlace * 1e-15, R.heatInPlace / 3.6e12); ui::propText("Heat in place above T_inj (model)", b);
            std::snprintf(b, sizeof(b), "%.4g GWh_th", cumGWh[N - 1]); ui::propText("Cumulative net heat produced", b);
            std::snprintf(b, sizeof(b), "%.3f %%", rec[N - 1]); ui::propText("Heat recovery factor", b);
            double avgMW = 0; if (R.summary[N - 1].time > 0) avgMW = cumGWh[N - 1] * 1000.0 / 24.0 / R.summary[N - 1].time;
            std::snprintf(b, sizeof(b), "%.3f MWth  (%.3f MWe at %.2f)", avgMW, avgMW * eff, eff); ui::propText("Time-averaged net power", b);
            std::snprintf(b, sizeof(b), "%.3g kg/s", qprod[N - 1]); ui::propText("Production mass rate (end)", b);
            if (lifeT > 0) std::snprintf(b, sizeof(b), "%.2f years (T_D = %.2f)", lifeT / 365.25, lifeDrop); else std::snprintf(b, sizeof(b), "beyond the run (T_D < %.2f)", lifeDrop);
            ui::propText("Thermal lifetime", b);
            std::snprintf(b, sizeof(b), "%.2f %%", cool[N - 1]); ui::propText("Matrix volume cooled by > 10 K (end)", b);
            std::snprintf(b, sizeof(b), "%zu fractures, %zu fracture cells", R.fracNames.size(), R.fracSegs.size()); ui::propText("EDFM", b);
            ui::endProps();
        }
        ui::hint("Net thermal power = produced - injected enthalpy rate (reference-state independent for a balanced doublet).\n"
                 "Heat recovery = cumulative net heat / heat stored above the injection temperature in the matrix (rock + water).\n"
                 "Impedance = (mean injector BHP - mean producer BHP) / production mass rate. Electric power = efficiency x net thermal.");
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
    (void)tSnap; (void)sub;
}
