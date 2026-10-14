// Inspector editors for every section of the case (no deck editing needed).
#include "app.h"
#include <functional>
#include <map>
#include <set>
#include "widgets.h"
#include "implot.h"
#include "ressim/fluid.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <thread>
#include <random>

using namespace ressim;

static void clampWellsToGrid(Case& c) {
    const GridInput& g = c.grid;
    for (auto& w : c.wells) {
        w.i = std::max(0, std::min(g.nx - 1, w.i)); w.j = std::max(0, std::min(g.ny - 1, w.j));
        for (auto& cp : w.completions) {
            cp.i = std::max(0, std::min(g.nx - 1, cp.i)); cp.j = std::max(0, std::min(g.ny - 1, cp.j));
            cp.k1 = std::max(0, std::min(g.nz - 1, cp.k1)); cp.k2 = std::max(cp.k1, std::min(g.nz - 1, cp.k2));
        }
    }
}

// ============================================================================ RUNSPEC
void App::editRunspec() {
    ui::header("Case");
    if (ui::beginProps("rs")) {
        ui::label("Title");
        if (ui::inputString("##title", &case_.title)) markDirty();
        static const char* units[] = {"METRIC (bar, m, sm3)", "FIELD (psia, ft, stb, Mscf)"};
        int u = case_.units == UnitSystem::Field ? 1 : 0;
        if (ui::propCombo("Unit system", &u, units, 2, "Values are NOT converted when switching; re-enter data in the new units.")) {
            case_.units = u ? UnitSystem::Field : UnitSystem::Metric; markDirty(true);
        }
        ui::endProps();
    }
    ui::header("Active phases");
    bool ch = false;
    ch |= ImGui::Checkbox("Oil", &case_.phases.oil); ImGui::SameLine();
    ch |= ImGui::Checkbox("Water", &case_.phases.water); ImGui::SameLine();
    ch |= ImGui::Checkbox("Gas", &case_.phases.gas); ImGui::SameLine();
    ImGui::BeginDisabled(!(case_.phases.gas && case_.phases.oil));
    ch |= ImGui::Checkbox("Dissolved gas (DISGAS)", &case_.phases.disgas); ImGui::SameLine();
    const bool vapWas = case_.phases.vapoil;
    ch |= ImGui::Checkbox("Vaporised oil (VAPOIL)", &case_.phases.vapoil);
    ImGui::EndDisabled();
    if (case_.phases.vapoil && !(case_.phases.gas && case_.phases.oil)) case_.phases.vapoil = false;
    if (case_.phases.vapoil && !vapWas && case_.props.pvtg.size() < 2) {
        // seed PVTG from PVDG: Rv = 0 everywhere (a dry gas expressed as wet gas) - the user edits Rv afterwards
        auto& P = case_.props; P.pvtg.clear();
        for (auto& r : P.pvdg) { PVTGRecord g; g.pg = r.p; g.rows.push_back({0.0, r.B, r.mu}); P.pvtg.push_back(g); }
        if (!P.pvtg.empty()) { auto& L = P.pvtg.back(); L.rows.front().rv = 1e-4; L.rows.push_back({0.0, L.rows.front().B * 0.99, L.rows.front().mu * 0.98}); }
    }
    if (ch) markDirty();
    ui::hint("Black-oil formulation: water, oil and gas components; gas may dissolve in oil (live oil, PVTO) and oil may vaporise "
             "in gas (wet gas / gas condensate, PVTG). With DISGAS off, oil uses PVDO (dead oil); with VAPOIL off, gas uses PVDG "
             "(dry gas). Supported systems: oil-water, gas-oil, 3-phase, single-phase water or oil.");
    ui::header("Model summary");
    if (ui::beginProps("sum")) {
        char b[128];
        std::snprintf(b, sizeof(b), "%d x %d x %d = %d cells", case_.grid.nx, case_.grid.ny, case_.grid.nz, case_.grid.size());
        ui::propText("Grid", b);
        if (preview_) { std::snprintf(b, sizeof(b), "%d", preview_->numActive()); ui::propText("Active cells", b); }
        std::snprintf(b, sizeof(b), "%d", int(case_.wells.size())); ui::propText("Wells", b);
        double T = 0; for (auto& s : case_.schedule) T += s.length;
        std::snprintf(b, sizeof(b), "%d report steps, %.0f days (%.2f years)", int(case_.schedule.size()), T, T / 365.25);
        ui::propText("Schedule", b);
        ui::endProps();
    }
}

// ============================================================================ GRID
void App::editGrid() {
    bool geomCh = false, rockCh = false;
    const bool cp = case_.grid.cornerPoint();
    if (cp) {
        ui::header("Geometry (corner-point COORD / ZCORN)");
        if (ui::beginProps("cpgeo")) {
            char b[128];
            std::snprintf(b, sizeof(b), "%d x %d x %d  (%d cells)", case_.grid.nx, case_.grid.ny, case_.grid.nz, case_.grid.size());
            ui::propText("Dimensions", b);
            std::snprintf(b, sizeof(b), "%d pillars, %d corner depths", (case_.grid.nx + 1) * (case_.grid.ny + 1), int(case_.grid.zcorn.size()));
            ui::propText("Geometry", b);
            if (preview_) { std::snprintf(b, sizeof(b), "%d (faults / explicit)", preview_->numNNC); ui::propText("Non-neighbour connections", b); }
            ui::endProps();
        }
        ui::hint("Corner-point geometry (imported deck, converted grid or fault throws). Layer DZ and dips come from ZCORN; "
                 "rock properties per layer can still be regenerated below.");
        if (ImGui::Button("Discard corner-point geometry (rebuild Cartesian)")) {
            try { generateGeometry(geom_, case_.grid); case_.grid.coord.clear(); case_.grid.zcorn.clear(); clampWellsToGrid(case_); markDirty(true); }
            catch (const std::exception& e) { gridError_ = e.what(); }
        }
    } else {
    ui::header("Geometry (Cartesian block-centred)");
    if (ui::beginProps("geo")) {
        geomCh |= ui::propI("NX (cells in I)", &geom_.nx, 1, 2000);
        geomCh |= ui::propI("NY (cells in J)", &geom_.ny, 1, 2000);
        geomCh |= ui::propI("NZ (layers)", &geom_.nz, 1, 500);
        geomCh |= ui::propD("DX (cell length I)", &geom_.dx, "%.4g");
        geomCh |= ui::propD("DY (cell length J)", &geom_.dy, "%.4g");
        geomCh |= ui::propD("Top depth", &geom_.top, "%.6g", "Depth of the top face of layer 1 at the origin corner (TOPS).");
        geomCh |= ui::propSliderD("Dip along I (deg)", &geom_.dipX, -30, 30, "%.1f");
        geomCh |= ui::propSliderD("Dip along J (deg)", &geom_.dipY, -30, 30, "%.1f");
        if (ui::propBool("Geometric transmissibility", &case_.grid.newtran,
                         "NEWTRAN: neighbouring columns at different depths connect through their actual face overlap "
                         "(non-neighbour connections where layers are offset). This is what OPM Flow always does.\n"
                         "Off = OLDTRAN (ECLIPSE default for block-centred grids): cell (i,j,k) connects to (i+1,j,k) "
                         "whatever the depth offset. Matters for dipping grids when the column offset approaches the layer thickness."))
            markDirty(true);
        ui::endProps();
    }
    if (ImGui::SmallButton("Convert to corner-point")) { convertToCornerPoint(case_.grid); markDirty(true); notify("Converted to corner-point geometry (COORD/ZCORN)"); }
    ImGui::SameLine(); ImGui::TextDisabled("enables fault throws below");
    }
    // ---- dual porosity / dual permeability (ECLIPSE DUALPORO / DUALPERM)
    ui::header("Dual porosity (naturally fractured reservoirs)");
    {
        GridInput& g = case_.grid;
        static double fPoro = 0.005, fPerm = 500, bL = 10, bH = 10; static bool fRegion = true;
        if (!g.dual) {
            ui::hint("Single porosity. Converting doubles NZ: layers NZ+1..2NZ become the fracture system of the same blocks "
                     "(ECLIPSE layout). Each matrix block couples to its fracture cell with T = CDARCY V sigma K.");
            if (ui::beginProps("dpnew")) {
                ui::propD("Fracture porosity", &fPoro, "%.4g", "Fracture pore volume per bulk volume (typically 0.1 - 1 %).");
                ui::propD("Fracture permeability (mD)", &fPerm, "%.5g", "Intrinsic (bulk) fracture-system permeability; NODPPM is set.");
                ui::propD("Matrix block size L (lateral)", &bL, "%.4g", "Fracture spacing in I and J (length units of the case).");
                ui::propD("Matrix block height H", &bH, "%.4g", "Fracture spacing in the vertical direction.");
                ui::propBool("Straight-line kr in the fractures", &fRegion, "Adds a SATNUM region for the fracture cells: krw = Sw, kro = So, Pc = 0.");
                ui::endProps();
            }
            if (ImGui::Button("Convert to dual porosity (DUALPORO)")) {
                makeDualPorosity(case_, 1, fPoro, fPerm, std::max(1e-3, bL), std::max(1e-3, bH), fRegion);
                syncBuildersFromCase(); markDirty(true);
                notify("Dual porosity: NZ doubled, fracture layers added (SIGMA from the block size)");
            }
            ImGui::SameLine();
            if (ImGui::Button("... dual permeability (DUALPERM)")) {
                makeDualPorosity(case_, 2, fPoro, fPerm, std::max(1e-3, bL), std::max(1e-3, bH), fRegion);
                syncBuildersFromCase(); markDirty(true);
                notify("Dual porosity / dual permeability model created");
            }
        } else {
            bool ch = false;
            if (ui::beginProps("dp")) {
                int mode = g.dual - 1;
                static const char* modes[] = {"Dual porosity (DUALPORO)", "Dual porosity + dual permeability (DUALPERM)"};
                if (ui::propCombo("Model", &mode, modes, 2, "DUALPORO: the matrix blocks only exchange fluid with their fracture cell, wells "
                                  "connect to the fractures.\nDUALPERM: matrix blocks also flow to each other and wells connect to both media.")) {
                    g.dual = mode + 1; ch = true;
                }
                char b[96]; std::snprintf(b, sizeof b, "%d matrix + %d fracture", g.nz / 2, g.nz / 2);
                ui::propText("Layers", b);
                int sm = g.ltosigma ? 1 : 0;
                static const char* sms[] = {"SIGMA (shape factor)", "From block sizes (LTOSIGMA, LX/LY/LZ)"};
                if (ui::propCombo("Coupling", &sm, sms, 2)) { g.ltosigma = sm == 1; ch = true; }
                if (!g.ltosigma) {
                    ch |= ui::propD("SIGMA (1/length^2)", &g.sigma, "%.5g",
                                    "Matrix-fracture interface area per unit volume. Kazemi: sigma = 4 (1/lx^2 + 1/ly^2 + 1/lz^2); "
                                    "Warren-Root / Gilman-Kazemi variants differ by a constant factor.");
                    if (!g.sigmav.empty()) ui::propText("SIGMAV", "per-cell values from the deck (override SIGMA)");
                } else {
                    double lx = g.lx.empty() ? 10.0 : g.lx[0], ly = g.ly.empty() ? 10.0 : g.ly[0], lz = g.lz.empty() ? 10.0 : g.lz[0];
                    bool lch = ui::propD("LX (block length I)", &lx, "%.4g");
                    lch |= ui::propD("LY (block length J)", &ly, "%.4g");
                    lch |= ui::propD("LZ (block height)", &lz, "%.4g");
                    if (lch) { g.lx.assign(size_t(g.size()), lx); g.ly.assign(size_t(g.size()), ly); g.lz.assign(size_t(g.size()), lz); ch = true; }
                    ch |= ui::propD("Factor fx", &g.ltoF[0], "%.3g"); ch |= ui::propD("Factor fy", &g.ltoF[1], "%.3g"); ch |= ui::propD("Factor fz", &g.ltoF[2], "%.3g");
                }
                ch |= ui::propBool("Intrinsic fracture permeability (NODPPM)", &g.nodppm,
                                   "Off (ECLIPSE default): the fracture PERMX/Y/Z are multiplied by the fracture porosity.");
                if (preview_) { std::snprintf(b, sizeof b, "%d", preview_->numMF); ui::propText("Matrix-fracture connections", b); }
                if (g.gravdr) ui::propText("Gravity drainage", "GRAVDR in the deck: not modelled by ResSim (kept for ECLIPSE)");
                ui::endProps();
            }
            if (ch) markDirty(true);
            ui::hint("Rock properties of the fracture system are the lower half of the layer table. OPM Flow does not support "
                     "dual porosity - run it with ResSim or ECLIPSE.");
            if (ImGui::SmallButton("Back to single porosity (drop the fracture layers)")) {
                makeSinglePorosity(case_); syncBuildersFromCase(); markDirty(true);
            }
        }
    }

    const int nz = geom_.nz;
    auto fit = [nz](std::vector<double>& v, double def) { if (int(v.size()) != nz) { const double last = v.empty() ? def : v.back(); v.resize(nz, last); } };
    fit(geom_.dz, 10); fit(rock_.layerPerm, 100); fit(rock_.layerPoro, 0.2); fit(rock_.layerKvKh, 0.1); fit(rock_.layerNtg, 1.0);

    ui::header("Layers");
    if (ImGui::BeginTable("layers", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                          ImVec2(0, std::min(260.0f, ImGui::GetFrameHeightWithSpacing() * (nz + 1.5f))))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("K", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.8f);
        ImGui::TableSetupColumn("DZ"); ImGui::TableSetupColumn("Perm (mD)"); ImGui::TableSetupColumn("Porosity");
        ImGui::TableSetupColumn("Kv/Kh"); ImGui::TableSetupColumn("NTG");
        ImGui::TableHeadersRow();
        for (int k = 0; k < nz; ++k) {
            ImGui::PushID(k);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::Text("%d", k + 1);
            double* cols[5] = {&geom_.dz[k], &rock_.layerPerm[k], &rock_.layerPoro[k], &rock_.layerKvKh[k], &rock_.layerNtg[k]};
            for (int c = 0; c < 5; ++c) {
                ImGui::TableSetColumnIndex(c + 1);
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::PushID(c);
                ImGui::BeginDisabled(cp && c == 0);
                if (ImGui::InputDouble("##v", cols[c], 0, 0, "%.4g")) { if (c == 0) geomCh = true; else rockCh = true; }
                ImGui::EndDisabled();
                ImGui::PopID();
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    ui::header("Heterogeneity (geostatistics)");
    if (ui::beginProps("het", 0.56f)) {
        rockCh |= ui::propSliderD("Perm variability  \xCF\x83(ln k)", &rock_.sigmaLn, 0.0, 2.5, "%.2f",
                                  "Standard deviation of ln(k). 0 = homogeneous layers; 1 = strongly heterogeneous. Mean permeability per layer is preserved.");
        rockCh |= ui::propSliderD("Correlation length I (cells)", &rock_.corrX, 0.0, 30.0, "%.1f");
        rockCh |= ui::propSliderD("Correlation length J (cells)", &rock_.corrY, 0.0, 30.0, "%.1f");
        rockCh |= ui::propSliderD("Poro-perm exponent", &rock_.poroPermExp, 0.0, 0.5, "%.2f", "poro = layer poro * (k / k_layer)^exponent");
        ui::label("Random seed");
        int seed = int(rock_.seed);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - ImGui::GetFontSize() * 5);
        if (ImGui::InputInt("##seed", &seed)) { rock_.seed = uint32_t(seed); rockCh = true; }
        ImGui::SameLine();
        if (ImGui::Button("New seed")) { rock_.seed = std::random_device{}(); rockCh = true; }
        ui::endProps();
    }
    ImGui::Spacing();
    ImGui::Checkbox("Auto-apply changes", &autoBuild_);
    ImGui::SameLine();
    const bool apply = ui::accentButton(cp ? "Regenerate rock properties" : "Build grid", ImVec2(-FLT_MIN, 0));
    if (apply || (autoBuild_ && (geomCh || rockCh))) {
        try {
            if (!cp) generateGeometry(geom_, case_.grid);
            generateRock(rock_, case_.grid);
            clampWellsToGrid(case_);
            markDirty(true);
        } catch (const std::exception& e) { gridError_ = e.what(); }
    }
    if (!gridError_.empty()) ImGui::TextColored(col4(Pal::Err), "%s", gridError_.c_str());

    if (preview_ && preview_->numCells() == case_.grid.size() && case_.grid.permx.size() == size_t(case_.grid.size())) {
        ui::header("Statistics");
        const Grid& g = *preview_;
        double pv = 0, bulk = 0, kavg = 0, kmin = 1e300, kmax = 0;
        const GridInput& gi = case_.grid;
        std::vector<double> logk;
        logk.reserve(g.numActive());
        for (int a = 0; a < g.numActive(); ++a) {
            const int q = g.globalOf[a];
            pv += g.poreVolRef[a];
            bulk += q < int(g.bulkVol.size()) ? g.bulkVol[q] : 0.0;
            kavg += gi.permx[q]; kmin = std::min(kmin, gi.permx[q]); kmax = std::max(kmax, gi.permx[q]);
            logk.push_back(std::log10(std::max(gi.permx[q], 1e-6)));
        }
        if (g.numActive()) kavg /= g.numActive();
        const Units U = Units::of(case_.units);
        if (ui::beginProps("stats")) {
            char b[96];
            std::snprintf(b, sizeof(b), "%d of %d", g.numActive(), g.numCells()); ui::propText("Active cells", b);
            std::snprintf(b, sizeof(b), "%.4g %s3", bulk, U.len); ui::propText("Bulk volume", b);
            std::snprintf(b, sizeof(b), "%.4g %s", pv, case_.units == UnitSystem::Field ? "rb" : "rm3"); ui::propText("Pore volume", b);
            std::snprintf(b, sizeof(b), "%.1f  (min %.2f, max %.1f) mD", kavg, kmin, kmax); ui::propText("Permeability X", b);
            std::snprintf(b, sizeof(b), "%d", int(g.conns.size())); ui::propText("Connections", b);
            std::snprintf(b, sizeof(b), "%.4g to %.4g %s", g.bbMin.z, g.bbMax.z, U.len); ui::propText("Depth range", b);
            ui::endProps();
        }
        if (ImPlot::BeginPlot("Permeability distribution", ImVec2(-1, ImGui::GetFontSize() * 11))) {
            ImPlot::SetupAxes("log10(PERMX, mD)", "cells", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlotSpec sp; sp.FillColor = col4(Pal::Accent); sp.FillAlpha = 0.75f; sp.LineColor = col4(Pal::Accent);
            if (!logk.empty()) ImPlot::PlotHistogram("##h", logk.data(), int(logk.size()), 30, 1.0, ImPlotRange(), sp);
            ImPlot::EndPlot();
        }
    }

    // ---- faults (FAULTS + MULTFLT) and fault-throw tool
    GridInput& gi = case_.grid;
    ui::header("Faults (FAULTS / MULTFLT)");
    {
        static int pos = 0; static double thr = 20; static double mult = 1.0; static bool alongJ = false; static char name[32] = "F1";
        if (pos <= 0) pos = std::max(1, gi.nx / 2);
        if (ImGui::TreeNodeEx("Add a fault with throw", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ui::beginProps("fthrow")) {
                static const char* dirs[] = {"Constant I (fault plane along J)", "Constant J (fault plane along I)"};
                int d = alongJ ? 1 : 0; if (ui::propCombo("Orientation", &d, dirs, 2)) alongJ = d == 1;
                ui::propI(alongJ ? "Between rows J and J+1" : "Between columns I and I+1", &pos, 1, (alongJ ? gi.ny : gi.nx) - 1);
                ui::propD("Throw (down on the far side)", &thr, "%.4g", "Vertical displacement of all cells beyond the fault plane. 0 = vertical barrier without offset.");
                ui::propSliderD("Transmissibility multiplier", &mult, 0.0, 1.0, "%.3f", "MULTFLT: 0 = sealing, 1 = open");
                ui::label("Name"); ImGui::InputText("##fn", name, sizeof(name));
                ui::endProps();
            }
            const int nmax = alongJ ? gi.ny : gi.nx;
            ImGui::BeginDisabled(nmax < 2);
            if (ui::accentButton("Add fault")) {
                try {
                    addFaultThrow(gi, std::min(pos, nmax - 1), thr, mult, name, alongJ);
                    syncBuildersFromCase();
                    markDirty(true);
                    notify(std::string("Fault ") + name + " added (corner-point grid)");
                    std::snprintf(name, sizeof(name), "F%d", int(gi.faults.size()) + 1);
                } catch (const std::exception& e) { gridError_ = e.what(); }
            }
            ImGui::EndDisabled();
            ImGui::TreePop();
        }
    }
    if (gi.faults.empty()) ui::hint("No faults defined.");
    else if (ImGui::BeginTable("faults", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Name", 0, 1.0f); ImGui::TableSetupColumn("MULTFLT", 0, 1.6f); ImGui::TableSetupColumn("Faces", 0, 0.8f);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
        ImGui::TableHeadersRow();
        int del = -1;
        for (int f = 0; f < int(gi.faults.size()); ++f) {
            FaultInput& F = gi.faults[f];
            ImGui::PushID(f);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::SetNextItemWidth(-FLT_MIN); if (ui::inputString("##n", &F.name)) markDirty(true);
            ImGui::TableSetColumnIndex(1); ImGui::SetNextItemWidth(-FLT_MIN);
            float m = float(F.mult);
            if (ImGui::SliderFloat("##m", &m, 0.0f, 1.0f, "%.3f")) { F.mult = m; markDirty(true); }
            ImGui::TableSetColumnIndex(2);
            int cells = 0; for (auto& fc : F.faces) cells += (fc.i2 - fc.i1 + 1) * (fc.j2 - fc.j1 + 1) * (fc.k2 - fc.k1 + 1);
            ImGui::Text("%d (%d)", int(F.faces.size()), cells);
            ImGui::TableSetColumnIndex(3); if (ImGui::SmallButton("x")) del = f;
            ImGui::PopID();
        }
        ImGui::EndTable();
        if (del >= 0) { gi.faults.erase(gi.faults.begin() + del); markDirty(true); }
        ui::hint("Removing a fault keeps the geometric throw (ZCORN) but removes its MULTFLT barrier.");
    }

    // ---- explicit non-neighbour connections
    ui::header("Explicit non-neighbour connections (NNC)");
    {
        int del = -1;
        if (!gi.nnc.empty() && ImGui::BeginTable("nnc", 8, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame)) {
            for (const char* h : {"I1", "J1", "K1", "I2", "J2", "K2", "Trans", ""}) ImGui::TableSetupColumn(h);
            ImGui::TableHeadersRow();
            for (int n = 0; n < int(gi.nnc.size()); ++n) {
                NNCInput& c = gi.nnc[n];
                ImGui::PushID(n);
                ImGui::TableNextRow();
                int* idx[6] = {&c.i1, &c.j1, &c.k1, &c.i2, &c.j2, &c.k2};
                const int lim[6] = {gi.nx, gi.ny, gi.nz, gi.nx, gi.ny, gi.nz};
                for (int q = 0; q < 6; ++q) {
                    ImGui::TableSetColumnIndex(q); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(q);
                    int v = *idx[q] + 1;
                    if (ImGui::DragInt("##i", &v, 0.1f, 1, lim[q])) { *idx[q] = std::clamp(v, 1, lim[q]) - 1; markDirty(true); }
                    ImGui::PopID();
                }
                ImGui::TableSetColumnIndex(6); ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::InputDouble("##t", &c.trans, 0, 0, "%.4g")) markDirty(true);
                ImGui::TableSetColumnIndex(7); if (ImGui::SmallButton("x")) del = n;
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (del >= 0) { gi.nnc.erase(gi.nnc.begin() + del); markDirty(true); }
        if (ImGui::Button("+ NNC")) { NNCInput c; c.i2 = std::min(1, gi.nx - 1); c.trans = 10; gi.nnc.push_back(c); markDirty(true); }
        ImGui::SameLine(); ImGui::TextDisabled("transmissibility in simulator units (cDarcy included), 1-based indices");
    }
}

// ============================================================================ AQUIFERS
void App::editAquifers() {
    auto& aq = case_.aquifers;
    const GridInput& g = case_.grid;
    bool ch = false;
    if (ImGui::Button("+ Fetkovich aquifer")) {
        AquiferInput a; a.id = int(aq.size()) + 1; a.type = AquiferType::Fetkovich;
        a.depth = case_.equil.owc; a.volume = case_.units == UnitSystem::Field ? 2e9 : 3e8; a.pi = case_.units == UnitSystem::Field ? 50 : 5;
        a.ct = case_.units == UnitSystem::Field ? 1e-5 : 1.5e-4;
        AquiferConnection c; c.i1 = 0; c.i2 = g.nx - 1; c.j1 = c.j2 = 0; c.k1 = 0; c.k2 = g.nz - 1; c.face = 'Y'; c.minus = true;
        a.conns.push_back(c); aq.push_back(a); ch = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("+ Carter-Tracy aquifer")) {
        AquiferInput a; a.id = int(aq.size()) + 1; a.type = AquiferType::CarterTracy;
        a.depth = case_.equil.owc; a.perm = 200; a.poro = 0.22; a.thickness = case_.units == UnitSystem::Field ? 100 : 30;
        a.rInner = case_.units == UnitSystem::Field ? 3000 : 900; a.angle = 90;
        a.ct = case_.units == UnitSystem::Field ? 1e-5 : 1.5e-4;
        AquiferConnection c; c.i1 = 0; c.i2 = g.nx - 1; c.j1 = c.j2 = g.ny - 1; c.k1 = 0; c.k2 = g.nz - 1; c.face = 'Y'; c.minus = false;
        a.conns.push_back(c); aq.push_back(a); ch = true;
    }
    if (aq.empty()) {
        ui::hint("No aquifers. A Fetkovich aquifer is a finite tank with a productivity index (pseudo-steady state); Carter-Tracy "
                 "approximates the van Everdingen-Hurst infinite-acting radial aquifer. Both connect to cell faces (AQUANCON).");
    }
    int del = -1;
    static const char* faces[] = {"I- (west)", "I+ (east)", "J- (south)", "J+ (north)", "K- (top)", "K+ (bottom)"};
    for (int n = 0; n < int(aq.size()); ++n) {
        AquiferInput& a = aq[n];
        ImGui::PushID(n);
        char hdr[96]; std::snprintf(hdr, sizeof(hdr), "Aquifer %d  -  %s###aq", a.id, a.type == AquiferType::Fetkovich ? "Fetkovich (AQUFETP)" : "Carter-Tracy (AQUCT)");
        if (ImGui::CollapsingHeader(hdr, ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ui::beginProps("aqp")) {
                static const char* types[] = {"Fetkovich", "Carter-Tracy"};
                int t = int(a.type); if (ui::propCombo("Model", &t, types, 2)) { a.type = AquiferType(t); ch = true; }
                ch |= ui::propI("ID", &a.id, 1, 99);
                ch |= ui::propD("Datum depth", &a.depth, "%.6g");
                bool eq = a.p0 < 0;
                ui::label("Initial pressure", "Unchecked: specify; checked: in equilibrium with the reservoir at the aquifer datum.");
                if (ImGui::Checkbox("Equilibrium", &eq)) { a.p0 = eq ? -1 : case_.equil.datumPressure; ch = true; }
                if (!eq) { ImGui::SameLine(); ImGui::SetNextItemWidth(-FLT_MIN); ch |= ImGui::InputDouble("##p0", &a.p0, 0, 0, "%.6g"); }
                ch |= ui::propD("Total compressibility", &a.ct, "%.4g");
                if (a.type == AquiferType::Fetkovich) {
                    ch |= ui::propD("Initial water volume", &a.volume, "%.4g", "Reservoir volume units (rb or rm3)");
                    ch |= ui::propD("Productivity index", &a.pi, "%.4g", "Reservoir volume / day / pressure");
                } else {
                    ch |= ui::propD("Permeability (mD)", &a.perm, "%.4g");
                    ch |= ui::propD("Porosity", &a.poro, "%.3g");
                    ch |= ui::propD("Inner radius", &a.rInner, "%.4g");
                    ch |= ui::propD("Thickness", &a.thickness, "%.4g");
                    ch |= ui::propSliderD("Angle of influence (deg)", &a.angle, 1, 360, "%.0f");
                    bool pv = a.visc < 0;
                    ui::label("Water viscosity");
                    if (ImGui::Checkbox("From PVTW", &pv)) { a.visc = pv ? -1 : case_.props.pvtw.visc; ch = true; }
                    if (!pv) { ImGui::SameLine(); ImGui::SetNextItemWidth(-FLT_MIN); ch |= ImGui::InputDouble("##mu", &a.visc, 0, 0, "%.4g"); }
                }
                ui::endProps();
            }
            ImGui::TextDisabled("Connections (AQUANCON), 1-based cell ranges");
            int dc = -1;
            if (ImGui::BeginTable("aqc", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                for (const char* h : {"I1", "I2", "J1", "J2", "K1", "K2", "Face", "Mult", ""}) ImGui::TableSetupColumn(h, 0, std::string(h) == "Face" ? 2.2f : 1.0f);
                ImGui::TableHeadersRow();
                for (int c = 0; c < int(a.conns.size()); ++c) {
                    AquiferConnection& C = a.conns[c];
                    ImGui::PushID(c);
                    ImGui::TableNextRow();
                    int* v[6] = {&C.i1, &C.i2, &C.j1, &C.j2, &C.k1, &C.k2};
                    const int lim[6] = {g.nx, g.nx, g.ny, g.ny, g.nz, g.nz};
                    for (int q = 0; q < 6; ++q) {
                        ImGui::TableSetColumnIndex(q); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(q);
                        int x = *v[q] + 1;
                        if (ImGui::DragInt("##v", &x, 0.1f, 1, lim[q])) { *v[q] = std::clamp(x, 1, lim[q]) - 1; ch = true; }
                        ImGui::PopID();
                    }
                    ImGui::TableSetColumnIndex(6); ImGui::SetNextItemWidth(-FLT_MIN);
                    int f = (C.face == 'X' ? 0 : C.face == 'Y' ? 2 : 4) + (C.minus ? 0 : 1);
                    if (ImGui::Combo("##f", &f, faces, 6)) { C.face = "XYZ"[f / 2]; C.minus = (f % 2) == 0; ch = true; }
                    ImGui::TableSetColumnIndex(7); ImGui::SetNextItemWidth(-FLT_MIN);
                    ch |= ImGui::InputDouble("##m", &C.mult, 0, 0, "%.3g");
                    ImGui::TableSetColumnIndex(8); if (ImGui::SmallButton("x")) dc = c;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (dc >= 0) { a.conns.erase(a.conns.begin() + dc); ch = true; }
            if (ImGui::SmallButton("+ connection")) { AquiferConnection c; c.i2 = g.nx - 1; c.k2 = g.nz - 1; c.face = 'Y'; c.minus = true; a.conns.push_back(c); ch = true; }
            ImGui::SameLine();
            if (ui::dangerButton("Remove aquifer")) del = n;
        }
        ImGui::PopID();
    }
    if (del >= 0) { aq.erase(aq.begin() + del); ch = true; }
    if (ch) markDirty();
    // influx results from the last run
    std::lock_guard<std::mutex> l(runner_.mutex());
    const RunResults& R = runner_.results();
    if (R.valid && R.summary.size() > 1 && !aq.empty()) {
        ui::header("Aquifer influx (last run)");
        std::vector<double> t, we, qa;
        const double xs = 1.0 / 365.25;
        for (auto& sp : R.summary) { t.push_back(sp.time * xs); we.push_back(sp.faqt); qa.push_back(sp.faqr); }
        if (ImPlot::BeginPlot("##aqinf", ImVec2(-1, ImGui::GetFontSize() * 13))) {
            ImPlot::SetupAxes("years", "cumulative influx", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "rate / day", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
            ImPlotSpec a; a.LineColor = col4(Pal::Water); a.LineWeight = 2;
            ImPlot::PlotLine("We (cum)", t.data(), we.data(), int(t.size()), a);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlotSpec b; b.LineColor = col4(Pal::Accent); b.LineWeight = 1.5f;
            ImPlot::PlotLine("influx rate", t.data(), qa.data(), int(t.size()), b);
            ImPlot::EndPlot();
        }
    }
}

// ============================================================================ PVT
void App::editPVT() {
    const Units U = Units::of(case_.units);
    PropsInput& p = case_.props;
    const bool live = case_.phases.oil && case_.phases.gas && case_.phases.disgas;
    if (live && !case_.schedule.empty()) {   // DRSDT from the start of the run (later changes: deck editor)
        ReportStep& s0 = case_.schedule[0];
        bool on = s0.drsdt >= 0;
        if (ImGui::Checkbox("Limit gas re-dissolution (DRSDT)", &on)) { s0.drsdt = on ? 0.0 : -2; markDirty(); }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Maximum rate at which Rs may increase in a cell (%s per day).\n"
                              "DRSDT 0: free gas cannot re-dissolve, the GOR of undersaturated oil cannot rise.\n"
                              "Common in field models to suppress unphysical re-dissolution in coarse cells.", U.rs);
        if (on) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            double v = s0.drsdt;
            if (ImGui::InputDouble("per day##drsdt", &v, 0, 0, "%.4g") && v >= 0) { s0.drsdt = v; markDirty(); }
            ImGui::SameLine();
            if (ImGui::Checkbox("only cells with free gas (FREE)", &s0.drsdtFree)) markDirty();
        }
    }
    if (!ImGui::BeginTabBar("pvt")) return;
    if (case_.phases.oil && ImGui::BeginTabItem(live ? "Oil (PVTO)" : "Oil (PVDO)")) {
        static int selRec = 0;
        if (live) {
            ui::hint("Saturated live-oil table: each row is a bubble point (Rs, Pb). Select a row to edit its undersaturated branch (p > Pb).");
            char c0[48], c1[48];
            std::snprintf(c0, sizeof(c0), "Rs (%s)", U.rs); std::snprintf(c1, sizeof(c1), "Pb (%s)", U.p);
            if (ui::editTable<PVTORecord>("pvto", p.pvto, {c0, c1, "Bo (rb/stb)", "\xCE\xBCo (cP)"},
                                          [](PVTORecord& r, int c) -> double* { return c == 0 ? &r.rs : c == 1 ? &r.pb : c == 2 ? &r.bo : &r.mu; },
                                          ImGui::GetFontSize() * 12, &selRec))
                markDirty();
            if (selRec >= 0 && selRec < int(p.pvto.size())) {
                ImGui::TextDisabled("Undersaturated branch of record %d (Rs = %.4g)", selRec + 1, p.pvto[selRec].rs);
                if (ui::editTable<PVDRow>("under", p.pvto[selRec].under, {"p", "Bo", "\xCE\xBCo"},
                                          [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; },
                                          ImGui::GetFontSize() * 6))
                    markDirty();
            }
        } else {
            if (ui::editTable<PVDRow>("pvdo", p.pvdo, {"p", "Bo", "\xCE\xBCo (cP)"},
                                      [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; }))
                markDirty();
        }
        // plots
        if (ImPlot::BeginSubplots("##oilplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
            if (ImPlot::BeginPlot("Bo")) {
                ImPlot::SetupAxes(U.p, "Bo", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (live) {
                    std::vector<double> x, y, rs;
                    for (auto& r : p.pvto) { x.push_back(r.pb); y.push_back(r.bo); rs.push_back(r.rs); }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Oil); sp.LineWeight = 2; sp.Marker = ImPlotMarker_Circle;
                    ImPlot::PlotLine("Saturated", x.data(), y.data(), int(x.size()), sp);
                    for (auto& r : p.pvto) if (!r.under.empty()) {
                        std::vector<double> ux{r.pb}, uy{r.bo};
                        for (auto& u : r.under) { ux.push_back(u.p); uy.push_back(u.B); }
                        ImPlotSpec s2; s2.LineColor = ImVec4(0.6f, 0.9f, 0.6f, 0.7f); s2.Flags = ImPlotItemFlags_NoLegend;
                        ImPlot::PlotLine("##u", ux.data(), uy.data(), int(ux.size()), s2);
                    }
                } else {
                    std::vector<double> x, y; for (auto& r : p.pvdo) { x.push_back(r.p); y.push_back(r.B); }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Oil); sp.LineWeight = 2;
                    ImPlot::PlotLine("Bo", x.data(), y.data(), int(x.size()), sp);
                }
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot(live ? "Rs and viscosity" : "Viscosity")) {
                ImPlot::SetupAxes(U.p, live ? U.rs : "cP", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (live) ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                std::vector<double> x, y, m;
                if (live) for (auto& r : p.pvto) { x.push_back(r.pb); y.push_back(r.rs); m.push_back(r.mu); }
                else for (auto& r : p.pvdo) { x.push_back(r.p); m.push_back(r.mu); }
                ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2;
                ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
                if (live) { ImPlot::PlotLine("Rs", x.data(), y.data(), int(x.size()), a); ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2); }
                ImPlot::PlotLine("\xCE\xBCo", x.data(), m.data(), int(x.size()), b);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::EndTabItem();
    }
    const bool vapo = case_.phases.gas && case_.phases.oil && case_.phases.vapoil;
    if (vapo && ImGui::BeginTabItem("Gas (PVTG)")) {
        static int selG = 0;
        char bc[48]; std::snprintf(bc, sizeof(bc), "Bg (%s)", case_.units == UnitSystem::Field ? "rb/Mscf" : "rm3/sm3");
        char rvc[48]; std::snprintf(rvc, sizeof(rvc), "Rv sat (%s)", U.rv);
        ui::hint("Wet-gas table: each row is a gas pressure with its saturated (dew-point) Rv, Bg and \xCE\xBCg. Select a row to edit its "
                 "undersaturated branch (Rv below saturation, down to 0). Records without a branch borrow the shape of the next "
                 "higher-pressure record that has one (as OPM Flow / ECLIPSE).");
        auto fix = [](PVTGRecord& r) { if (r.rows.empty()) r.rows.push_back({0.0, 1.0, 0.02}); };
        for (auto& r : p.pvtg) fix(r);
        bool chg = ui::editTable<PVTGRecord>("pvtg", p.pvtg, {U.p, rvc, bc, "\xCE\xBCg (cP)"},
            [&](PVTGRecord& r, int c) -> double* { fix(r); return c == 0 ? &r.pg : c == 1 ? &r.rows.front().rv : c == 2 ? &r.rows.front().B : &r.rows.front().mu; },
            ImGui::GetFontSize() * 12, &selG);
        if (selG >= 0 && selG < int(p.pvtg.size())) {
            PVTGRecord& R = p.pvtg[size_t(selG)];
            ImGui::TextDisabled("Undersaturated branch of record %d (p = %.5g, Rv sat = %.4g)", selG + 1, R.pg, R.rows.front().rv);
            std::vector<PVTGRow> und(R.rows.begin() + 1, R.rows.end());
            char rvu[48]; std::snprintf(rvu, sizeof(rvu), "Rv (%s)", U.rv);
            if (ui::editTable<PVTGRow>("pvtgu", und, {rvu, "Bg", "\xCE\xBCg"},
                                       [](PVTGRow& r, int c) -> double* { return c == 0 ? &r.rv : c == 1 ? &r.B : &r.mu; },
                                       ImGui::GetFontSize() * 6)) {
                R.rows.resize(1); R.rows.insert(R.rows.end(), und.begin(), und.end()); chg = true;
            }
        }
        if (chg) {
            p.pvdg.clear(); for (auto& g : p.pvtg) if (!g.rows.empty()) p.pvdg.push_back({g.pg, g.rows.front().B, g.rows.front().mu});
            markDirty();
        }
        if (ImPlot::BeginSubplots("##gasplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 14))) {
            std::vector<double> x, bg, rv, mu;
            for (auto& g : p.pvtg) if (!g.rows.empty()) { x.push_back(g.pg); bg.push_back(g.rows.front().B); rv.push_back(g.rows.front().rv); mu.push_back(g.rows.front().mu); }
            if (ImPlot::BeginPlot("Bg (saturated)")) {
                ImPlot::SetupAxes(U.p, bc, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                ImPlot::PlotLine("Bg", x.data(), bg.data(), int(x.size()), a);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Dew-point Rv and viscosity")) {
                ImPlot::SetupAxes(U.p, U.rv, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                ImPlotSpec a; a.LineColor = col4(Pal::Oil); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle;
                ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
                ImPlot::PlotLine("Rv sat", x.data(), rv.data(), int(x.size()), a);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                ImPlot::PlotLine("\xCE\xBCg", x.data(), mu.data(), int(x.size()), b);
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::EndTabItem();
    }
    if (case_.phases.gas && !vapo && ImGui::BeginTabItem("Gas (PVDG)")) {
        char bc[48]; std::snprintf(bc, sizeof(bc), "Bg (%s)", case_.units == UnitSystem::Field ? "rb/Mscf" : "rm3/sm3");
        if (ui::editTable<PVDRow>("pvdg", p.pvdg, {U.p, bc, "\xCE\xBCg (cP)"},
                                  [](PVDRow& r, int c) -> double* { return c == 0 ? &r.p : c == 1 ? &r.B : &r.mu; }))
            markDirty();
        if (ImPlot::BeginPlot("Dry gas", ImVec2(-1, ImGui::GetFontSize() * 14))) {
            ImPlot::SetupAxes(U.p, bc, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "cP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            std::vector<double> x, y, m; for (auto& r : p.pvdg) { x.push_back(r.p); y.push_back(r.B); m.push_back(r.mu); }
            ImPlotSpec a; a.LineColor = col4(Pal::Gas); a.LineWeight = 2;
            ImPlotSpec b; b.LineColor = col4(Pal::Warn); b.LineWeight = 2;
            ImPlot::PlotLine("Bg", x.data(), y.data(), int(x.size()), a);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlot::PlotLine("\xCE\xBCg", x.data(), m.data(), int(x.size()), b);
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Water, rock, densities")) {
        bool ch = false;
        ui::header("Water (PVTW)");
        if (ui::beginProps("pvtw")) {
            ch |= ui::propD("Reference pressure", &p.pvtw.pref);
            ch |= ui::propD("Bw at reference", &p.pvtw.bw);
            ch |= ui::propD("Compressibility cw (1/p)", &p.pvtw.cw, "%.4g");
            ch |= ui::propD("Viscosity \xCE\xBCw (cP)", &p.pvtw.visc);
            ch |= ui::propD("Viscosibility (1/p)", &p.pvtw.cv, "%.4g");
            ui::endProps();
        }
        ui::header("Rock (ROCK)");
        if (ui::beginProps("rock")) {
            ch |= ui::propD("Reference pressure", &p.rock.pref);
            ch |= ui::propD("Pore compressibility cr (1/p)", &p.rock.cr, "%.4g");
            ui::endProps();
        }
        ui::header("Surface densities (DENSITY)");
        if (ui::beginProps("dens")) {
            ch |= ui::propD("Oil", &p.density.oil); ch |= ui::propD("Water", &p.density.water); ch |= ui::propD("Gas", &p.density.gas, "%.5g");
            char b[64];
            const double rhoW = case_.units == UnitSystem::Field ? 62.37 : 999.0;
            const double sg = p.density.oil / rhoW;
            std::snprintf(b, sizeof(b), "%.1f \xC2\xB0" "API", 141.5 / sg - 131.5); ui::propText("Oil gravity", b);
            std::snprintf(b, sizeof(b), "%.3f (air = 1)", p.density.gas / (case_.units == UnitSystem::Field ? 0.0764 : 1.225)); ui::propText("Gas gravity", b);
            ui::endProps();
        }
        if (ch) markDirty();
        ImGui::EndTabItem();
    }
    if (case_.phases.oil && ImGui::BeginTabItem("Correlations")) {
        ui::hint("Generate a complete black-oil PVT description from a few field measurements. "
                 "Standing (Rs, Bo), Beggs-Robinson + Vasquez-Beggs (oil viscosity), Vasquez-Beggs (undersaturated co), "
                 "Sutton pseudo-criticals + Hall-Yarborough (z-factor, Bg), Lee-Gonzalez-Eakin (gas viscosity).");
        if (ui::beginProps("corr")) {
            ui::propD("Oil gravity (\xC2\xB0" "API)", &corr_.api, "%.1f");
            ui::propD("Gas specific gravity (air=1)", &corr_.gasGravity, "%.3f");
            ui::propD("Reservoir temperature (\xC2\xB0" "F)", &corr_.temperatureF, "%.1f");
            ui::propD("Max bubble point (psia)", &corr_.pbMax, "%.0f");
            ui::propD("Max table pressure (psia)", &corr_.pMax, "%.0f");
            ui::propI("Table rows", &corr_.nRows, 4, 40);
            ui::endProps();
        }
        if (ui::accentButton("Generate PVT tables", ImVec2(-FLT_MIN, 0))) {
            try {
                makeCorrelationPVT(corr_, case_.units, case_.props, !live);
                markDirty();
                notify("PVT tables generated from correlations");
            } catch (const std::exception& e) { notify(std::string("Correlation error: ") + e.what()); }
        }
        ImGui::TextDisabled("Pressures are entered in psia; tables are converted to the case units.");
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ============================================================================ SCAL
void App::editSCAL() {
    PropsInput& p = case_.props;
    const bool wo = case_.phases.water && case_.phases.oil, go = case_.phases.gas && case_.phases.oil;
    if (!wo && !go) { ui::hint("Relative permeability is only needed with two or more phases."); return; }
    if (ImGui::CollapsingHeader("Corey generator", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginProps("corey", 0.56f)) {
            ui::propSliderD("Swc (connate water)", &corey_.swc, 0.0, 0.5, "%.3f");
            if (wo) {
                ui::propSliderD("Sorw (residual oil to water)", &corey_.sorw, 0.0, 0.5, "%.3f");
                ui::propSliderD("nw (water exponent)", &corey_.nw, 1.0, 6.0, "%.2f");
                ui::propSliderD("now (oil exponent, o/w)", &corey_.now, 1.0, 6.0, "%.2f");
                ui::propSliderD("krw at Sorw", &corey_.krwMax, 0.05, 1.0, "%.2f");
            }
            ui::propSliderD("kro at Swc", &corey_.kroMax, 0.05, 1.0, "%.2f");
            if (go) {
                ui::propSliderD("Sgc (critical gas)", &corey_.sgc, 0.0, 0.3, "%.3f");
                ui::propSliderD("Sorg (residual oil to gas)", &corey_.sorg, 0.0, 0.5, "%.3f");
                ui::propSliderD("ng (gas exponent)", &corey_.ng, 1.0, 6.0, "%.2f");
                ui::propSliderD("nog (oil exponent, g/o)", &corey_.nog, 1.0, 6.0, "%.2f");
                ui::propSliderD("krg at max Sg", &corey_.krgMax, 0.05, 1.0, "%.2f");
            }
            if (wo) {
                ui::propD("Pcow entry pressure", &corey_.pcEntry, "%.4g", "Capillary pressure at Swc (0 = no capillary pressure).");
                ui::propSliderD("Pcow exponent", &corey_.pcExponent, 0.5, 6.0, "%.2f");
            }
            ui::endProps();
        }
        if (ui::accentButton("Generate SWOF / SGOF", ImVec2(-FLT_MIN, 0))) {
            try {
                if (wo) makeCoreySWOF(corey_, p.swof);
                if (go) makeCoreySGOF(corey_, p.sgof);
                markDirty();
            } catch (const std::exception& e) { notify(e.what()); }
        }
    }
    if (wo && go) {
        static const char* kModels[] = {"ECLIPSE default (saturation-weighted)", "Stone I (STONE1)", "Stone II (STONE2)"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16);
        if (ImGui::Combo("Three-phase oil kr", &p.kroModel, kModels, 3)) markDirty();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("How kro is built from krow (SWOF) and krog (SGOF) when water, oil and gas flow together.\n"
                              "Default: kro = (Sg krog + (Sw - Swc) krow) / (Sg + Sw - Swc).\n"
                              "Stone II: kro = krocw ((krow/krocw + krw)(krog/krocw + krg) - krw - krg).\n"
                              "Stone I: kro = beta krow krog / krocw (normalised saturations, exponent STONE1EX).");
        if (p.kroModel == 1) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            double e = p.stone1Exp;
            if (ImGui::InputDouble("STONE1EX", &e, 0, 0, "%.3g") && e > 0) { p.stone1Exp = e; markDirty(); }
        }
    }
    if (ImGui::BeginTabBar("scal")) {
        if (wo && ImGui::BeginTabItem("Water-oil (SWOF)")) {
            if (ui::editTable<SatRow>("swof", p.swof, {"Sw", "krw", "krow", "Pcow"},
                                      [](SatRow& r, int c) -> double* { return c == 0 ? &r.s : c == 1 ? &r.krw : c == 2 ? &r.kro : &r.pc; },
                                      ImGui::GetFontSize() * 10))
                markDirty();
            if (ImPlot::BeginSubplots("##wo", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 15))) {
                std::vector<double> s, a, b, pc;
                for (auto& r : p.swof) { s.push_back(r.s); a.push_back(r.krw); b.push_back(r.kro); pc.push_back(r.pc); }
                if (ImPlot::BeginPlot("Relative permeability")) {
                    ImPlot::SetupAxes("Sw", "kr", 0, 0);
                    ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                    ImPlotSpec sw; sw.LineColor = col4(Pal::Water); sw.LineWeight = 2.5f;
                    ImPlotSpec so; so.LineColor = col4(Pal::Oil); so.LineWeight = 2.5f;
                    ImPlot::PlotLine("krw", s.data(), a.data(), int(s.size()), sw);
                    ImPlot::PlotLine("krow", s.data(), b.data(), int(s.size()), so);
                    ImPlot::EndPlot();
                }
                // fractional flow + Welge tangent (Buckley-Leverett)
                if (ImPlot::BeginPlot("Fractional flow (Buckley-Leverett)")) {
                    ImPlot::SetupAxes("Sw", "fw", 0, 0);
                    ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                    double muw = p.pvtw.visc, muo = 1.0;
                    if (case_.phases.gas && case_.phases.disgas && !p.pvto.empty()) muo = p.pvto.back().mu;
                    else if (!p.pvdo.empty()) muo = p.pvdo[p.pvdo.size() / 2].mu;
                    std::vector<double> fs, fw;
                    double swc = s.empty() ? 0 : s.front();
                    double bestSlope = 0, sf = swc, fwf = 0;
                    for (int i = 0; i <= 200 && !s.empty(); ++i) {
                        const double sv = swc + (s.back() - swc) * i / 200.0;
                        auto interp = [&](const std::vector<double>& y) {
                            size_t k = std::upper_bound(s.begin(), s.end(), sv) - s.begin();
                            if (k == 0) return y.front();
                            if (k >= s.size()) return y.back();
                            const double t = (sv - s[k - 1]) / (s[k] - s[k - 1]); return y[k - 1] + t * (y[k] - y[k - 1]);
                        };
                        const double lw = interp(a) / muw, lo = interp(b) / muo;
                        const double f = (lw + lo) > 0 ? lw / (lw + lo) : 1.0;
                        fs.push_back(sv); fw.push_back(f);
                        if (sv > swc + 1e-6) { const double sl = f / (sv - swc); if (sl > bestSlope) { bestSlope = sl; sf = sv; fwf = f; } }
                    }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Water); sp.LineWeight = 2.5f;
                    ImPlot::PlotLine("fw", fs.data(), fw.data(), int(fs.size()), sp);
                    if (bestSlope > 0) {
                        const double tx[2] = {swc, std::min(1.0, swc + 1.0 / bestSlope)}, ty[2] = {0, std::min(1.0, bestSlope * (tx[1] - swc))};
                        ImPlotSpec ts; ts.LineColor = col4(Pal::Warn); ts.LineWeight = 1.5f;
                        ImPlot::PlotLine("Welge tangent", tx, ty, 2, ts);
                        ImPlotSpec ms; ms.Marker = ImPlotMarker_Circle; ms.MarkerSize = 5; ms.MarkerFillColor = col4(Pal::Warn);
                        ImPlot::PlotScatter("Shock front", &sf, &fwf, 1, ms);
                        char b2[96]; std::snprintf(b2, sizeof(b2), "Swf = %.3f, fwf = %.3f", sf, fwf);
                        ImPlot::PlotText(b2, 0.3, 0.08);
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            ImGui::EndTabItem();
        }
        if (go && ImGui::BeginTabItem("Gas-oil (SGOF)")) {
            if (ui::editTable<SatRow>("sgof", p.sgof, {"Sg", "krg", "krog", "Pcog"},
                                      [](SatRow& r, int c) -> double* { return c == 0 ? &r.s : c == 1 ? &r.krw : c == 2 ? &r.kro : &r.pc; },
                                      ImGui::GetFontSize() * 10))
                markDirty();
            if (ImPlot::BeginPlot("Gas-oil relative permeability", ImVec2(-1, ImGui::GetFontSize() * 15))) {
                ImPlot::SetupAxes("Sg", "kr", 0, 0);
                ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
                std::vector<double> s, a, b;
                for (auto& r : p.sgof) { s.push_back(r.s); a.push_back(r.krw); b.push_back(r.kro); }
                ImPlotSpec sg; sg.LineColor = col4(Pal::Gas); sg.LineWeight = 2.5f;
                ImPlotSpec so; so.LineColor = col4(Pal::Oil); so.LineWeight = 2.5f;
                ImPlot::PlotLine("krg", s.data(), a.data(), int(s.size()), sg);
                ImPlot::PlotLine("krog", s.data(), b.data(), int(s.size()), so);
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("End-point scaling")) {
            editEndPointScaling();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
}

// End-point scaling: ENDSCALE + per-cell saturation end points (constant or per layer), live preview of the scaled
// curves against the tables of SATNUM region 1.
void App::editEndPointScaling() {
    GridInput& g = case_.grid;
    const bool wo = case_.phases.water && case_.phases.oil, go = case_.phases.gas && case_.phases.oil;
    const float fs = ImGui::GetFontSize();
    ImGui::TextWrapped("Scales the relative permeability and capillary pressure curves per cell to its own saturation end points "
                       "(ECLIPSE two-point scaling: kr between the critical and maximum saturation, Pc between SWL and SWU). "
                       "Cells without a value use the table end points of their SATNUM region. Written as ENDSCALE + SWL ... in the deck.");
    if (ImGui::Checkbox("Enable end-point scaling (ENDSCALE)", &g.endscale)) markDirty();
    SatFunctions sf;
    bool sfOk = true;
    try { sf.setup(case_, 0); } catch (...) { sfOk = false; }
    const SatFunctions::EndPoints T = sfOk ? sf.tableEndPoints() : SatFunctions::EndPoints{};
    struct Arr { const char* name; const char* meaning; std::vector<double>* v; double T; bool show; };
    Arr arrs[] = {
        {"SWL", "connate water (Pc and kro end)", &g.swl, T.swl, wo}, {"SWCR", "critical water (krw = 0 below)", &g.swcr, T.swcr, wo},
        {"SWU", "maximum water", &g.swu, T.swu, wo}, {"SOWCR", "residual oil to water", &g.sowcr, T.sowcr, wo},
        {"SGL", "connate gas", &g.sgl, T.sgl, go}, {"SGCR", "critical gas (krg = 0 below)", &g.sgcr, T.sgcr, go},
        {"SGU", "maximum gas", &g.sgu, T.sgu, go}, {"SOGCR", "residual oil to gas", &g.sogcr, T.sogcr, go},
    };
    const int n = g.size(), nz = std::max(1, g.nz), perLayer = std::max(1, n / nz);
    static int sel = 0;
    static double constVal = 0.2;
    if (ImGui::BeginTable("##eps", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("Array", ImGuiTableColumnFlags_WidthFixed, fs * 4);
        ImGui::TableSetupColumn("Meaning", ImGuiTableColumnFlags_WidthFixed, fs * 12);
        ImGui::TableSetupColumn("Table", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
        ImGui::TableSetupColumn("Cells", ImGuiTableColumnFlags_WidthFixed, fs * 9);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();
        for (int i = 0; i < 8; ++i) {
            Arr& a = arrs[i];
            if (!a.show) continue;
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            if (ImGui::Selectable(a.name, sel == i, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) sel = i;
            ImGui::TableSetColumnIndex(1); ImGui::TextDisabled("%s", a.meaning);
            ImGui::TableSetColumnIndex(2); ImGui::Text("%.3f", a.T);
            ImGui::TableSetColumnIndex(3);
            if (a.v->empty()) ImGui::TextDisabled("table");
            else {
                double lo = 1e9, hi = -1e9; int set = 0;
                for (double x : *a.v) if (x >= 0) { lo = std::min(lo, x); hi = std::max(hi, x); ++set; }
                if (set == 0) ImGui::TextDisabled("table");
                else if (hi - lo < 1e-12) ImGui::Text("%.3f%s", lo, set < n ? " (part)" : "");
                else ImGui::Text("%.3f - %.3f", lo, hi);
            }
            ImGui::TableSetColumnIndex(4);
            if (ImGui::SmallButton("Constant")) { a.v->assign(size_t(n), std::clamp(constVal, 0.0, 1.0)); g.endscale = true; markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Set every cell to the value below");
            ImGui::SameLine();
            if (ImGui::SmallButton("Per layer")) { sel = i; if (a.v->empty()) a.v->assign(size_t(n), a.T); g.endscale = true; markDirty(); }
            ImGui::SameLine();
            ImGui::BeginDisabled(a.v->empty());
            if (ImGui::SmallButton("Clear")) { a.v->clear(); markDirty(); }
            ImGui::EndDisabled();
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    ImGui::SetNextItemWidth(fs * 8);
    ImGui::InputDouble("Value for 'Constant'", &constVal, 0.01, 0.05, "%.3f");
    Arr* Ap = &arrs[std::clamp(sel, 0, 7)];
    if (!Ap->show) for (auto& a : arrs) if (a.show) { Ap = &a; break; }
    Arr& A = *Ap;
    // per-layer editor of the selected array
    if (!A.v->empty() && int(A.v->size()) == n && ImGui::TreeNodeEx("##layers", ImGuiTreeNodeFlags_DefaultOpen, "%s per layer", A.name)) {
        if (ImGui::BeginTable("##epsl", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit,
                              ImVec2(fs * 14, std::min(fs * 12, fs * 1.6f * (nz + 1))))) {
            ImGui::TableSetupColumn("Layer", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
            ImGui::TableSetupColumn(A.name, ImGuiTableColumnFlags_WidthFixed, fs * 9);
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableHeadersRow();
            for (int k = 0; k < nz; ++k) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("%d", k + 1);
                ImGui::TableSetColumnIndex(1);
                double v = (*A.v)[size_t(k) * size_t(perLayer)];
                ImGui::PushID(k);
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::InputDouble("##v", &v, 0, 0, v < 0 ? "table" : "%.4f")) {
                    for (int q = 0; q < perLayer; ++q) (*A.v)[size_t(k) * size_t(perLayer) + size_t(q)] = std::min(1.0, v);
                    markDirty();
                }
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        ImGui::TextDisabled("Values < 0 = table end point. Cell-by-cell values come from imported decks (SWL ... / EQUALS / BOX).");
        ImGui::TreePop();
    }
    // preview: tables (dashed) against the scaled curves of the average cell end points
    if (!sfOk) { ui::hint("Saturation tables are incomplete - no preview."); return; }
    auto mean = [&](const std::vector<double>& v, double def) {
        double s = 0; int c = 0;
        for (double x : v) if (x >= 0) { s += x; ++c; }
        return c ? s / c : def;
    };
    SatFunctions::EndPoints E;
    E.swl = mean(g.swl, T.swl); E.swcr = std::max(E.swl, mean(g.swcr, T.swcr)); E.swu = std::max(E.swcr + 1e-6, mean(g.swu, T.swu));
    E.sgl = mean(g.sgl, T.sgl); E.sgcr = std::max(E.sgl, mean(g.sgcr, T.sgcr)); E.sgu = std::max(E.sgcr + 1e-6, mean(g.sgu, T.sgu));
    E.sowcr = mean(g.sowcr, T.sowcr); E.sogcr = mean(g.sogcr, T.sogcr);
    const int N = 120;
    if (ImPlot::BeginPlot("Table (dashed) vs scaled with the average cell end points", ImVec2(-1, std::max(fs * 14, ImGui::GetContentRegionAvail().y)))) {
        ImPlot::SetupAxes("Saturation", "kr", 0, 0);
        ImPlot::SetupAxesLimits(0, 1, 0, 1, ImPlotCond_Always);
        auto line = [&](const char* lab, ImU32 col, bool dashed, auto f) {
            std::vector<double> x(N + 1), y(N + 1);
            for (int i = 0; i <= N; ++i) { x[size_t(i)] = double(i) / N; y[size_t(i)] = f(x[size_t(i)]); }
            ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = dashed ? 1.2f : 2.5f;
            if (dashed) sp.LineColor.w = 0.55f;
            ImPlot::PlotLine(lab, x.data(), y.data(), N + 1, sp);
        };
        if (wo) {
            line("krw table", Pal::Water, true, [&](double s) { return sf.krw(AD3(s), nullptr).v; });
            line("krow table", Pal::Oil, true, [&](double s) { return sf.kro(AD3(s), AD3(0.0), nullptr).v; });
            line("krw scaled", Pal::Water, false, [&](double s) { return sf.krw(AD3(s), &E).v; });
            line("krow scaled", Pal::Oil, false, [&](double s) { return sf.kro(AD3(s), AD3(0.0), &E).v; });
        }
        if (go) {
            line("krg table (vs Sg)", Pal::Gas, true, [&](double s) { return sf.krg(AD3(s), nullptr).v; });
            line("krg scaled (vs Sg)", Pal::Gas, false, [&](double s) { return sf.krg(AD3(s), &E).v; });
        }
        ImPlot::EndPlot();
    }
}

// ============================================================================ INITIALISATION
void App::editInit() {
    EquilInput& e = case_.equil;
    const Units U = Units::of(case_.units);
    bool ch = false;
    ui::header("Hydrostatic equilibrium (EQUIL)");
    if (ui::beginProps("eq")) {
        ch |= ui::propD("Datum depth", &e.datumDepth, "%.6g");
        ch |= ui::propD("Pressure at datum", &e.datumPressure, "%.6g");
        if (case_.phases.water && case_.phases.oil) {
            ch |= ui::propD("Oil-water contact (OWC)", &e.owc, "%.6g");
            ch |= ui::propD("Pcow at OWC", &e.pcowOWC, "%.4g");
        }
        if (case_.phases.gas && case_.phases.oil) {
            ch |= ui::propD("Gas-oil contact (GOC)", &e.goc, "%.6g");
            ch |= ui::propD("Pcog at GOC", &e.pcogGOC, "%.4g");
        }
        ch |= ui::propI("Accuracy (NACC)", &e.nacc, 0, 20, "0: saturations at cell centre. N > 0: averaged over 2N horizontal slices per cell (smoother transition zones).");
        if (case_.phases.gas && case_.phases.oil && case_.phases.disgas) {
            bool sat = e.initialRs < 0;
            ui::label("Initial Rs", "Saturated at the GOC pressure is the ECLIPSE default; otherwise use a constant Rs (RSVD).");
            if (ImGui::Checkbox("Saturated at GOC", &sat)) { e.initialRs = sat ? -1 : 1.0; ch = true; }
            if (!sat) { ui::label("Constant Rs"); ch |= ImGui::InputDouble("##rs", &e.initialRs, 0, 0, "%.5g"); }
        }
        if (case_.phases.gas && case_.phases.oil && case_.phases.vapoil) {
            bool sat = e.initialRv < 0;
            ui::label("Initial Rv", "Saturated at the GOC pressure (capped by the local dew point) is the ECLIPSE default; "
                                    "otherwise a constant Rv. An RVVD table below overrides both.");
            if (ImGui::Checkbox("Saturated at GOC##rv", &sat)) { e.initialRv = sat ? -1 : 0.0; ch = true; }
            if (!sat) { ui::label("Constant Rv"); ch |= ImGui::InputDouble("##rvc", &e.initialRv, 0, 0, "%.5g"); }
        }
        ui::endProps();
    }
    // depth tables (RSVD / RVVD) - optional, override the constant / saturated choice
    struct DV { double d, v; };
    auto depthTable = [&](const char* id, const char* title, const char* vcol, std::vector<double>& dep, std::vector<double>& val) {
        if (!ImGui::TreeNode(id, "%s (%d rows)", title, int(dep.size()))) return;
        std::vector<DV> t; for (size_t i = 0; i < dep.size() && i < val.size(); ++i) t.push_back({dep[i], val[i]});
        if (ui::editTable<DV>(id, t, {U.len, vcol}, [](DV& r, int c) -> double* { return c == 0 ? &r.d : &r.v; })) {
            dep.clear(); val.clear(); for (auto& r : t) { dep.push_back(r.d); val.push_back(r.v); } ch = true;
        }
        if (!dep.empty() && ImGui::SmallButton("Clear table")) { dep.clear(); val.clear(); ch = true; }
        ImGui::TreePop();
    };
    if (case_.phases.gas && case_.phases.oil && case_.phases.disgas) {
        char c[48]; std::snprintf(c, sizeof(c), "Rs (%s)", U.rs);
        depthTable("rsvd", "Rs vs depth (RSVD)", c, e.rsvdDepth, e.rsvdRs);
    }
    if (case_.phases.gas && case_.phases.oil && case_.phases.vapoil) {
        char c[48]; std::snprintf(c, sizeof(c), "Rv (%s)", U.rv);
        depthTable("rvvd", "Rv vs depth (RVVD)", c, e.rvvdDepth, e.rvvdRv);
    }
    if (ch) markDirty();
    ui::hint("Phase pressures are integrated from the datum with dP/dz = \xCF\x81(P) g; saturations follow from the capillary pressure "
             "curves. Contacts above the reservoir top give undersaturated oil (no free gas).");
    if (ui::accentButton("Compute initial state", ImVec2(-FLT_MIN, 0))) computeInit();
    if (!init_.error.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextColored(col4(Pal::Err), "%s", init_.error.c_str()); ImGui::PopTextWrapPos(); }
    if (init_.valid) {
        if (init_.caseVersion != caseVersion_) ImGui::TextColored(col4(Pal::Warn), "Inputs changed since this initialisation.");
        ui::header("Fluids in place");
        if (ui::beginProps("fip")) {
            char b[96];
            auto big = [&](double v, const char* u) { if (std::fabs(v) >= 1e6) std::snprintf(b, sizeof(b), "%.4f MM %s", v / 1e6, u); else std::snprintf(b, sizeof(b), "%.6g %s", v, u); return b; };
            if (case_.phases.oil) ui::propText("Oil (STOIIP)", big(init_.ip.oil, U.liq));
            if (case_.phases.water) ui::propText("Water", big(init_.ip.water, U.liq));
            if (case_.phases.gas) ui::propText("Free gas", big(init_.ip.gasFree, U.gas));
            if (case_.phases.disgas && case_.phases.gas) ui::propText("Dissolved gas", big(init_.ip.gasDissolved, U.gas));
            if (case_.phases.vapoil && case_.phases.gas) {
                ui::propText("Vaporised oil (in gas)", big(init_.ip.oilVaporized, U.liq));
                ui::propText("Total oil (liquid + vapour)", big(init_.ip.oil + init_.ip.oilVaporized, U.liq));
            }
            ui::endProps();
        }
        if (ImPlot::BeginSubplots("##initplots", 1, 2, ImVec2(-1, ImGui::GetFontSize() * 18), ImPlotSubplotFlags_LinkAllY)) {
            const int n = int(init_.depth.size());
            if (ImPlot::BeginPlot("Pressure vs depth")) {
                ImPlot::SetupAxes(U.p, U.len, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                std::vector<double> pp(init_.p.begin(), init_.p.end());
                ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 2.5f; s.MarkerFillColor = col4(Pal::Pressure); s.MarkerLineColor = col4(Pal::Pressure);
                ImPlot::PlotScatter("cells", pp.data(), init_.depth.data(), n, s);
                ImPlotSpec h; h.Flags = ImPlotInfLinesFlags_Horizontal; h.LineColor = col4(Pal::Water);
                if (case_.phases.water) ImPlot::PlotInfLines("OWC", &e.owc, 1, h);
                h.LineColor = col4(Pal::Gas);
                if (case_.phases.gas) ImPlot::PlotInfLines("GOC", &e.goc, 1, h);
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Saturation vs depth")) {
                ImPlot::SetupAxes("S", U.len, 0, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0, 1, ImPlotCond_Always);
                std::vector<double> sw(init_.sw.begin(), init_.sw.end()), sg(init_.sg.begin(), init_.sg.end());
                ImPlotSpec a; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 2.5f; a.MarkerFillColor = col4(Pal::Water); a.MarkerLineColor = col4(Pal::Water);
                ImPlot::PlotScatter("Sw", sw.data(), init_.depth.data(), n, a);
                if (case_.phases.gas) {
                    ImPlotSpec b2 = a; b2.MarkerFillColor = col4(Pal::Gas); b2.MarkerLineColor = col4(Pal::Gas);
                    ImPlot::PlotScatter("Sg", sg.data(), init_.depth.data(), n, b2);
                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
    }
}

// ============================================================================ WELLS
// ============================================================================ GROUPS
// Group tree (WELSPECS groups + GRUPTREE) and group controls (GCONPROD / GCONINJE / ... kept as deck keywords for
// OPM Flow; edited with schema forms on the Advanced page).
void App::editGroups() {
    namespace di = ressim::deckintel;
    auto& X = case_.extraKeywords;
    // GRUPTREE (first one in the schedule is the editable tree)
    int treeIdx = -1;
    for (int i = 0; i < int(X.size()); ++i) if (X[size_t(i)].name == "GRUPTREE") { treeIdx = i; break; }
    std::map<std::string, std::string> parent;
    di::Form tree;
    if (treeIdx >= 0) {
        tree = di::parseForm("GRUPTREE", di::splitLines(X[size_t(treeIdx)].body));
        for (auto& r : tree.records) if (r.size() >= 2) parent[di::unquote(r[0])] = di::unquote(r[1]) == "1*" ? "FIELD" : di::unquote(r[1]);
    }
    std::set<std::string> groups;
    for (auto& w : case_.wells) if (!w.group.empty()) groups.insert(w.group);
    for (auto& [c, p] : parent) { groups.insert(c); if (p != "FIELD") groups.insert(p); }
    groups.erase("FIELD");
    auto parentOf = [&](const std::string& g) { auto it = parent.find(g); return it == parent.end() ? std::string("FIELD") : it->second; };
    auto isDescendant = [&](std::string a, const std::string& of) {   // a below 'of'?
        for (int d = 0; d < 64 && a != "FIELD"; ++d) { a = parentOf(a); if (a == of) return true; }
        return false;
    };
    const ImGuiTreeNodeFlags tf = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_SpanAvailWidth;
    std::function<void(const std::string&, int)> node = [&](const std::string& g, int depth) {
        if (depth > 32) return;
        ImGuiTreeNodeFlags f = tf | (selGroup_ == g ? ImGuiTreeNodeFlags_Selected : 0);
        const bool open = ImGui::TreeNodeEx(g.c_str(), f, "%s", g.c_str());
        if (ImGui::IsItemClicked() && !ImGui::IsItemToggledOpen()) selGroup_ = g;
        if (!open) return;
        for (auto& c : groups) if (parentOf(c) == g) node(c, depth + 1);
        for (int i = 0; i < int(case_.wells.size()); ++i) {
            const auto& w = case_.wells[size_t(i)];
            if (w.group != g) continue;
            ImGui::PushID(i);
            const ImU32 col = w.type == WellType::Producer ? Pal::Oil : w.type == WellType::WaterInjector ? Pal::Water : Pal::Gas;
            ImGui::TextColored(col4(col), "%s", w.type == WellType::Producer ? "P " : "I ");
            ImGui::SameLine(0, 0);
            if (ImGui::Selectable(w.name.c_str(), selWell_ == i)) selWell_ = i;
            ImGui::PopID();
        }
        ImGui::TreePop();
    };
    const float fs = ImGui::GetFontSize();
    ImGui::BeginChild("##gtree", ImVec2(fs * 13, fs * 14), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    node("FIELD", 0);
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("##gedit", ImVec2(0, fs * 14), ImGuiChildFlags_Borders);
    if (selGroup_.empty() || (selGroup_ != "FIELD" && !groups.count(selGroup_))) selGroup_ = "FIELD";
    ImGui::Text("Group %s", selGroup_.c_str());
    if (selGroup_ != "FIELD") {
        ImGui::SameLine(0, fs * 2);
        ImGui::SetNextItemWidth(fs * 9);
        const std::string cur = parentOf(selGroup_);
        if (ImGui::BeginCombo("Parent", cur.c_str())) {
            std::vector<std::string> cands = {"FIELD"};
            for (auto& g : groups) if (g != selGroup_ && !isDescendant(g, selGroup_)) cands.push_back(g);
            for (auto& c : cands) {
                if (!ImGui::Selectable(c.c_str(), c == cur) || c == cur) continue;
                // update / create GRUPTREE
                if (treeIdx < 0) {
                    ExtraKeyword e; e.name = "GRUPTREE"; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group tree (ResSim Studio)";
                    X.insert(X.begin(), e); treeIdx = 0;
                    tree = di::parseForm("GRUPTREE", {});
                }
                bool found = false;
                for (auto& r : tree.records) if (r.size() >= 2 && di::unquote(r[0]) == selGroup_) { r[1] = "'" + c + "'"; found = true; }
                if (!found) tree.records.push_back({"'" + selGroup_ + "'", "'" + c + "'"});
                X[size_t(treeIdx)].body = di::joinLines(di::writeForm(tree)) + "\n";
                markDirty();
            }
            ImGui::EndCombo();
        }
    }
    // wells of the group: assignment
    int nw = 0;
    for (auto& w : case_.wells) nw += w.group == selGroup_;
    ImGui::TextDisabled("%d well(s). Move the selected well into this group:", nw);
    ImGui::SameLine();
    const bool canMove = selWell_ >= 0 && selWell_ < int(case_.wells.size()) && case_.wells[size_t(selWell_)].group != selGroup_ && selGroup_ != "FIELD";
    ImGui::BeginDisabled(!canMove);
    if (ImGui::SmallButton(canMove ? ("Move " + case_.wells[size_t(selWell_)].name).c_str() : "Move well")) { case_.wells[size_t(selWell_)].group = selGroup_; markDirty(); }
    ImGui::EndDisabled();
    static char newGroup[32] = "";
    ImGui::SetNextItemWidth(fs * 7);
    ImGui::InputTextWithHint("##ng", "new group", newGroup, sizeof(newGroup), ImGuiInputTextFlags_CharsUppercase);
    ImGui::SameLine();
    if (ImGui::SmallButton("+ Group under this") && newGroup[0]) {
        const std::string ng = newGroup;
        if (!groups.count(ng) && ng != "FIELD") {
            if (treeIdx < 0) {
                ExtraKeyword e; e.name = "GRUPTREE"; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group tree (ResSim Studio)";
                X.insert(X.begin(), e); treeIdx = 0;
                tree = di::parseForm("GRUPTREE", {});
            }
            tree.records.push_back({"'" + ng + "'", "'" + selGroup_ + "'"});
            X[size_t(treeIdx)].body = di::joinLines(di::writeForm(tree)) + "\n";
            selGroup_ = ng;
            newGroup[0] = 0;
            markDirty();
        }
    }
    // group controls of this group
    ImGui::SeparatorText("Group controls");
    static const char* kGroupKw[] = {"GCONPROD", "GCONINJE", "GCONSUMP", "GCONSALE", "GEFAC", "GRUPTARG", "GCONPRI", "GUIDERAT"};
    int shown = 0;
    for (int i = 0; i < int(X.size()); ++i) {
        const ExtraKeyword& e = X[size_t(i)];
        bool isG = false;
        for (auto* k : kGroupKw) if (e.name == k) isG = true;
        if (!isG) continue;
        const di::Form f = di::parseForm(e.name, di::splitLines(e.body));
        for (auto& r : f.records) {
            if (r.empty()) continue;
            const std::string g = di::unquote(r[0]);
            if (g != selGroup_ && !(g.find('*') != std::string::npos && selGroup_.rfind(g.substr(0, g.find('*')), 0) == 0)) continue;
            std::string txt;
            for (size_t k = 1; k < r.size() && k < 8; ++k) txt += r[k] + " ";
            ImGui::PushID(i * 1000 + shown);
            ImGui::Text("%s  step %d:", e.name.c_str(), e.step);
            ImGui::SameLine();
            ImGui::TextDisabled("%s", txt.c_str());
            ImGui::SameLine();
            if (ImGui::SmallButton("Edit")) { advSel_ = i; advFormView_ = true; section_ = Section::Advanced; }
            ImGui::PopID();
            ++shown;
        }
    }
    if (!shown) ImGui::TextDisabled("none");
    auto addCtl = [&](const char* kw, const std::string& rec) {
        ExtraKeyword e; e.name = kw; e.section = DeckSection::Schedule; e.step = 0; e.beforeControls = true; e.note = "group control (ResSim Studio)";
        e.body = rec + "\n/\n";
        X.push_back(e);
        advSel_ = int(X.size()) - 1; advFormView_ = true; section_ = Section::Advanced;
        markDirty();
    };
    const bool fieldU = case_.units == UnitSystem::Field;
    if (ImGui::SmallButton("+ Production target (GCONPROD)")) addCtl("GCONPROD", " '" + selGroup_ + "' ORAT " + (fieldU ? "5000" : "800") + " 1* 1* 1* RATE /");
    ImGui::SameLine();
    if (ImGui::SmallButton("+ Water injection (GCONINJE)")) addCtl("GCONINJE", " '" + selGroup_ + "' WATER VREP 1* 1* 1* 1.0 /");
    if (case_.phases.gas) { ImGui::SameLine(); if (ImGui::SmallButton("+ Gas injection")) addCtl("GCONINJE", " '" + selGroup_ + "' GAS REIN 1* 1* 0.8 /"); }
    ImGui::TextColored(col4(Pal::Warn), "Group controls are simulated by OPM Flow. The ResSim engine runs every well on its own targets.");
    ImGui::EndChild();
}

void App::editWells() {
    const GridInput& g = case_.grid;
    if (ImGui::CollapsingHeader("Groups & group controls")) editGroups();
    auto addWell = [&](WellType t) {
        WellInput w;
        int n = 1;
        const char* prefix = t == WellType::Producer ? "PROD" : t == WellType::WaterInjector ? "WINJ" : "GINJ";
        for (;;) {
            std::string nm = prefix + std::to_string(n);
            bool used = false; for (auto& x : case_.wells) if (x.name == nm) used = true;
            if (!used) { w.name = nm; break; }
            ++n;
        }
        w.type = t; w.i = g.nx / 2; w.j = g.ny / 2;
        Completion cp; cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = g.nz - 1;
        w.completions.push_back(cp);
        case_.wells.push_back(w);
        if (case_.schedule.empty()) case_.schedule.push_back(ReportStep{});
        WellControl wc; wc.well = w.name;
        if (t == WellType::Producer) { wc.mode = Control::ORAT; wc.orat = case_.units == UnitSystem::Field ? 1000 : 150; wc.bhp = case_.units == UnitSystem::Field ? 1000 : 70; }
        else { wc.mode = Control::RATE; wc.rate = case_.units == UnitSystem::Field ? 1000 : 150; wc.bhp = case_.units == UnitSystem::Field ? 6000 : 400; }
        case_.schedule[0].controls.push_back(wc);
        selWell_ = int(case_.wells.size()) - 1;
        markDirty();
    };
    if (ImGui::Button("+ Producer")) addWell(WellType::Producer);
    ImGui::SameLine();
    if (ImGui::Button("+ Water injector")) addWell(WellType::WaterInjector);
    ImGui::SameLine();
    ImGui::BeginDisabled(!case_.phases.gas);
    if (ImGui::Button("+ Gas injector")) addWell(WellType::GasInjector);
    ImGui::EndDisabled();
    if (case_.wells.empty()) { ui::hint("No wells yet. Add a producer and an injector, then set their controls in the Schedule."); return; }
    selWell_ = std::max(0, std::min(selWell_, int(case_.wells.size()) - 1));

    // ---- areal map with wells (click to place)
    static bool placeMode = false;
    ui::header("Well locations (top view)");
    ImGui::Checkbox("Click on the map to move the selected well", &placeMode);
    if (preview_ && preview_->numCells() == g.size() && ImPlot::BeginPlot("##map", ImVec2(-1, ImGui::GetFontSize() * 16), ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("I", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
        std::vector<double> m(size_t(g.nx) * g.ny);
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                double s = 0; int cnt = 0;
                for (int k = 0; k < g.nz; ++k) { const int q = g.index(i, j, k); if (preview_->activeOf[q] >= 0) { s += std::log10(std::max(g.permx[q], 1e-3)); ++cnt; } }
                m[size_t(j) * g.nx + i] = cnt ? s / cnt : NAN;
            }
        // inactive columns: NaN must never reach the colormap (it produced a negative table index -> assertion);
        // they get the lowest value here and are painted grey on top
        double lo = 1e300, hi = -1e300;
        for (double v : m) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
        if (lo > hi) { lo = 0; hi = 1; } else if (hi - lo < 1e-9) hi = lo + 1;
        std::vector<char> inactive(m.size(), 0);
        for (size_t q = 0; q < m.size(); ++q) if (!std::isfinite(m[q])) { m[q] = lo; inactive[q] = 1; }
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        ImPlot::PlotHeatmap("logk", m.data(), g.ny, g.nx, lo, hi, nullptr, ImPlotPoint(0.5, g.ny + 0.5), ImPlotPoint(g.nx + 0.5, 0.5));
        ImPlot::PopColormap();
        {
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            ImPlot::PushPlotClipRect();
            for (int j = 0; j < g.ny; ++j)
                for (int i = 0; i < g.nx; ++i)
                    if (inactive[size_t(j) * g.nx + i]) {
                        const ImVec2 a = ImPlot::PlotToPixels(i + 0.5, j + 0.5), b = ImPlot::PlotToPixels(i + 1.5, j + 1.5);
                        dl->AddRectFilled(ImVec2(std::min(a.x, b.x), std::min(a.y, b.y)), ImVec2(std::max(a.x, b.x), std::max(a.y, b.y)), IM_COL32(48, 50, 56, 255));
                    }
            ImPlot::PopPlotClipRect();
        }
        for (int w = 0; w < int(case_.wells.size()); ++w) {
            const auto& wl = case_.wells[w];
            double x = wl.i + 1, y = wl.j + 1;
            ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = w == selWell_ ? 7.0f : 5.0f;
            s.MarkerFillColor = col4(wl.type == WellType::Producer ? Pal::Oil : wl.type == WellType::WaterInjector ? Pal::Water : Pal::Gas);
            s.MarkerLineColor = w == selWell_ ? ImVec4(1, 1, 1, 1) : ImVec4(0, 0, 0, 1);
            s.LineWeight = 2.0f;
            ImPlot::PlotScatter(wl.name.c_str(), &x, &y, 1, s);
            ImPlot::PlotText(wl.name.c_str(), x, y, ImVec2(0, -14));
        }
        if (placeMode && ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            const int i = std::max(0, std::min(g.nx - 1, int(std::lround(mp.x)) - 1));
            const int j = std::max(0, std::min(g.ny - 1, int(std::lround(mp.y)) - 1));
            WellInput& w = case_.wells[selWell_];
            w.i = i; w.j = j;
            for (auto& cp : w.completions) { cp.i = i; cp.j = j; }
            markDirty();
        }
        ImPlot::EndPlot();
    }

    // ---- list + details
    ImGui::BeginChild("wlist", ImVec2(ImGui::GetFontSize() * 8, ImGui::GetFontSize() * 11.5f), ImGuiChildFlags_Borders);
    for (int w = 0; w < int(case_.wells.size()); ++w) {
        const auto& wl = case_.wells[w];
        const ImU32 c = wl.type == WellType::Producer ? Pal::Oil : wl.type == WellType::WaterInjector ? Pal::Water : Pal::Gas;
        ImGui::PushStyleColor(ImGuiCol_Text, col4(c));
        ImGui::TextUnformatted("\xE2\x97\x8F"); ImGui::PopStyleColor(); ImGui::SameLine();
        if (ImGui::Selectable(wl.name.c_str(), selWell_ == w)) selWell_ = w;
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("wdet", ImVec2(0, ImGui::GetFontSize() * 11.5f));
    WellInput& w = case_.wells[selWell_];
    if (ui::beginProps("well")) {
        ui::label("Name");
        std::string nm = w.name;
        if (ui::inputString("##name", &nm, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsNoBlank)) {
            bool clash = false; for (auto& x : case_.wells) if (&x != &w && x.name == nm) clash = true;
            if (!nm.empty() && !clash) {
                for (auto& s : case_.schedule) for (auto& c : s.controls) if (c.well == w.name) c.well = nm;
                w.name = nm; markDirty();
            }
        }
        static const char* types[] = {"Producer", "Water injector", "Gas injector"};
        int t = int(w.type);
        if (ui::propCombo("Type", &t, types, 3)) {
            w.type = WellType(t);
            for (auto& s : case_.schedule) for (auto& c : s.controls) if (c.well == w.name) {
                if (w.type == WellType::Producer && c.mode == Control::RATE) { c.mode = Control::ORAT; c.orat = std::max(c.rate, 0.0); }
                if (w.type != WellType::Producer && c.mode != Control::RATE && c.mode != Control::BHP) { c.mode = Control::RATE; c.rate = std::max(c.orat, 0.0); }
            }
            markDirty();
        }
        int ij[2] = {w.i + 1, w.j + 1};
        ui::label("Location I, J");
        if (ImGui::InputInt2("##ij", ij)) {
            w.i = std::max(0, std::min(g.nx - 1, ij[0] - 1)); w.j = std::max(0, std::min(g.ny - 1, ij[1] - 1));
            for (auto& cp : w.completions) { cp.i = w.i; cp.j = w.j; }
            markDirty();
        }
        if (ui::propD("PI multiplier (WPIMULT)", &w.piMult, "%.4g", "Multiplies all connection factors of the well (stimulation / damage, history matching).")) {
            w.piMult = std::max(0.0, w.piMult); markDirty();
        }
        bool autoRef = w.refDepth <= 0;
        ui::label("BHP reference depth");
        if (ImGui::Checkbox("Top perforation", &autoRef)) { w.refDepth = autoRef ? -1 : case_.equil.datumDepth; markDirty(); }
        if (!autoRef) { ImGui::SameLine(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##ref", &w.refDepth, 0, 0, "%.6g")) markDirty(); }
        ui::label("WELSPECS items 7+ (OPM)", "Written after the preferred phase in WELSPECS, for OPM Flow:\n"
                  "drainage radius, inflow equation, auto shut-in (SHUT/STOP), crossflow (YES/NO), PVT table, density calc (SEG/AVG), ...\n"
                  "e.g.  0.0 STD SHUT NO 1* SEG   - not used by the ResSim engine.");
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ui::inputString("##wtail", w.specsTail)) markDirty();
        ui::endProps();
    }
    ImGui::EndChild();
    ui::header("Completions");
    ImGui::PushTextWrapPos(0);
    ImGui::TextDisabled("Layers K1..K2 (1-based). CF <= 0: Peaceman well index from diameter, skin and cell properties. "
                        "From / Until: report step (1-based) at which the connection is perforated / shut (late COMPDAT); "
                        "Until 0 = never.");
    ImGui::PopTextWrapPos();
    int del = -1;
    const int nsteps = std::max(1, int(case_.schedule.size()));
    static bool opmCols = false;
    bool anyOpm = false;
    for (auto& q : w.completions) if (q.satTable > 0 || q.dfactor >= 0 || q.r0 > 0) anyOpm = true;
    ImGui::Checkbox("OPM Flow connection items (sat. table, D-factor, r0)", &opmCols);
    if (anyOpm && !opmCols) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "(set for this well)"); }
    const int ncol = opmCols ? 13 : 10;
    const float fsz = ImGui::GetFontSize();
    if (ImGui::BeginTable("comp", ncol, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollX,
                          ImVec2(0, std::min(fsz * 16, ImGui::GetFrameHeightWithSpacing() * (float(w.completions.size()) + 1.6f) + fsz)))) {
        const float cw[] = {2.6f, 2.6f, 3.2f, 4.2f, 4.2f, 3.2f, 2.6f, 2.8f, 2.8f, 3.0f, 3.8f, 3.2f, 1.6f};
        const char* hd[] = {"K1", "K2", "Skin", "Diameter", "CF", "Dir", "Open", "From", "Until", "Sat tbl", "D-factor", "r0", ""};
        for (int h = 0; h < 13; ++h) {
            if (!opmCols && h >= 9 && h <= 11) continue;
            ImGui::TableSetupColumn(hd[h], ImGuiTableColumnFlags_WidthFixed, fsz * cw[h]);
        }
        ImGui::TableSetupScrollFreeze(2, 1);
        ImGui::TableHeadersRow();
        for (int c = 0; c < int(w.completions.size()); ++c) {
            Completion& cp = w.completions[c];
            ImGui::PushID(c);
            ImGui::TableNextRow();
            int k1 = cp.k1 + 1, k2 = cp.k2 + 1;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt("##k1", &k1, 0)) { cp.k1 = std::max(0, std::min(g.nz - 1, k1 - 1)); cp.k2 = std::max(cp.k1, cp.k2); markDirty(); }
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt("##k2", &k2, 0)) { cp.k2 = std::max(cp.k1, std::min(g.nz - 1, k2 - 1)); markDirty(); }
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##s", &cp.skin, 0, 0, "%.3g")) markDirty();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##d", &cp.diameter, 0, 0, "%.3g")) markDirty();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputDouble("##cf", &cp.cf, 0, 0, "%.4g")) markDirty();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            int dir = cp.dir == 'X' ? 0 : cp.dir == 'Y' ? 1 : 2;
            static const char* dirs[] = {"X", "Y", "Z"};
            if (ImGui::Combo("##dir", &dir, dirs, 3)) { cp.dir = "XYZ"[dir]; markDirty(); }
            ImGui::TableNextColumn(); if (ImGui::Checkbox("##open", &cp.open)) markDirty();
            int from = cp.fromStep + 1, until = cp.untilStep < 0 ? 0 : cp.untilStep + 1;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt("##from", &from, 0)) { cp.fromStep = std::max(0, std::min(nsteps - 1, from - 1)); markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Report step at which this connection opens (1 = from the start)");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt("##until", &until, 0)) { cp.untilStep = until <= 0 ? -1 : std::max(cp.fromStep + 1, std::min(nsteps, until - 1)); markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Report step at which this connection is shut (0 = stays open)");
            if (opmCols) {
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputInt("##sat", &cp.satTable, 0)) { cp.satTable = std::max(0, cp.satTable); markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("COMPDAT item 6: saturation table for the connection's relative permeability (0 = the cell's SATNUM). OPM Flow only.");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputDouble("##dfac", &cp.dfactor, 0, 0, "%.4g")) markDirty();
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("COMPDAT item 11: non-Darcy D-factor (negative = default / none). OPM Flow only.");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputDouble("##r0", &cp.r0, 0, 0, "%.4g")) markDirty();
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("COMPDAT item 13: pressure-equivalent radius r0 (negative = Peaceman). OPM Flow only.");
            }
            ImGui::TableNextColumn(); if (ImGui::SmallButton("\xC3\x97")) del = c;
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (del >= 0) { w.completions.erase(w.completions.begin() + del); markDirty(); }
    if (ImGui::SmallButton("+ Add completion")) {
        Completion cp; cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = g.nz - 1; w.completions.push_back(cp); markDirty();
    }
    ImGui::SameLine();
    if (ui::dangerButton("Delete well")) {
        const std::string nm = w.name;
        for (auto& s : case_.schedule)
            s.controls.erase(std::remove_if(s.controls.begin(), s.controls.end(), [&](const WellControl& c) { return c.well == nm; }), s.controls.end());
        case_.wells.erase(case_.wells.begin() + selWell_);
        selWell_ = std::max(0, selWell_ - 1);
        markDirty();
    }
}

// ============================================================================ SCHEDULE
void App::editSchedule() {
    static double totalYears = 10, stepDays = 91.3125;
    ui::header("Report steps");
    double T = 0; for (auto& s : case_.schedule) T += s.length;
    ImGui::Text("%d report steps, total %.0f days (%.2f years)", int(case_.schedule.size()), T, T / 365.25);
    if (ImGui::TreeNode("Rebuild uniform steps")) {
        if (ui::beginProps("rb")) {
            ui::propD("Total time (years)", &totalYears, "%.4g");
            ui::propD("Step length (days)", &stepDays, "%.4g");
            ui::endProps();
        }
        if (ImGui::Button("Rebuild (keeps well controls at their times)")) {
            const int n = std::max(1, int(std::ceil(totalYears * 365.25 / std::max(1e-3, stepDays) - 1e-9)));
            std::vector<std::pair<double, WellControl>> ctrls;
            double t = 0;
            for (auto& s : case_.schedule) { for (auto& c : s.controls) ctrls.push_back({t, c}); t += s.length; }
            std::vector<ReportStep> ns(n);
            for (auto& s : ns) s.length = stepDays;
            for (auto& [tc, c] : ctrls) { int k = std::min(n - 1, int(std::floor(tc / stepDays + 1e-6))); ns[k].controls.push_back(c); }
            case_.schedule = ns;
            selStep_ = 0;
            markDirty();
        }
        ImGui::TreePop();
    }
    if (case_.schedule.empty()) case_.schedule.push_back(ReportStep{});
    selStep_ = std::max(0, std::min(selStep_, int(case_.schedule.size()) - 1));

    // timeline table
    if (ImGui::BeginTable("steps", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                          ImVec2(0, ImGui::GetFontSize() * 11))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 2.2f);
        ImGui::TableSetupColumn("Start (days)"); ImGui::TableSetupColumn("Length (days)"); ImGui::TableSetupColumn("Control changes");
        ImGui::TableHeadersRow();
        double t = 0;
        for (int s = 0; s < int(case_.schedule.size()); ++s) {
            const ReportStep& st = case_.schedule[s];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            char b[16]; std::snprintf(b, sizeof(b), "%d", s + 1);
            if (ImGui::Selectable(b, selStep_ == s, ImGuiSelectableFlags_SpanAllColumns)) selStep_ = s;
            ImGui::TableNextColumn(); ImGui::Text("%.2f", t);
            ImGui::TableNextColumn(); ImGui::Text("%.4g", st.length);
            ImGui::TableNextColumn();
            if (!st.controls.empty()) {
                std::string l;
                for (auto& c : st.controls) { if (!l.empty()) l += ", "; l += c.well; }
                ImGui::TextColored(col4(Pal::Accent), "%s", l.c_str());
            }
            t += st.length;
        }
        ImGui::EndTable();
    }
    ReportStep& st = case_.schedule[selStep_];
    ui::header(("Step " + std::to_string(selStep_ + 1)).c_str());
    if (ui::beginProps("st")) {
        if (ui::propD("Length (days)", &st.length, "%.6g")) { st.length = std::max(1e-3, st.length); markDirty(); }
        ui::endProps();
    }
    ImGui::SameLine();
    ImGui::TextDisabled("Controls set here stay in force until changed in a later step.");
    if (ImGui::SmallButton("Insert step after")) { case_.schedule.insert(case_.schedule.begin() + selStep_ + 1, ReportStep{st.length, {}}); markDirty(); }
    ImGui::SameLine();
    if (case_.schedule.size() > 1 && ImGui::SmallButton("Delete step")) {
        auto ctrls = st.controls;
        case_.schedule.erase(case_.schedule.begin() + selStep_);
        selStep_ = std::max(0, selStep_ - 1);
        for (auto& c : ctrls) case_.schedule[selStep_].controls.push_back(c);
        markDirty();
        return;
    }

    ui::header("Well controls in this step");
    int del = -1;
    for (int ci = 0; ci < int(st.controls.size()); ++ci) {
        WellControl& c = st.controls[ci];
        ImGui::PushID(ci);
        const WellInput* w = nullptr; for (auto& x : case_.wells) if (x.name == c.well) w = &x;
        const bool prod = w ? w->type == WellType::Producer : c.mode != Control::RATE;
        ImGui::BeginChild("ctl", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
        if (ImGui::BeginCombo("##well", c.well.c_str())) {
            for (auto& x : case_.wells) if (ImGui::Selectable(x.name.c_str(), x.name == c.well)) { c.well = x.name; markDirty(); }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        static const char* stats[] = {"OPEN", "SHUT", "STOP"};
        int sidx = int(c.status);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
        if (ImGui::Combo("##status", &sidx, stats, 3)) { c.status = WellStatus(sidx); markDirty(); }
        ImGui::SameLine();
        static const char* pm[] = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"};
        static const char* im[] = {"RATE", "BHP"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
        if (prod) {
            int m = std::min(int(c.mode), 5);
            if (ImGui::Combo("##mode", &m, pm, 6)) { c.mode = Control(m); markDirty(); }
        } else {
            int m = c.mode == Control::BHP ? 1 : 0;
            if (ImGui::Combo("##mode", &m, im, 2)) { c.mode = m ? Control::BHP : Control::RATE; markDirty(); }
        }
        ImGui::SameLine();
        if (prod) { if (ImGui::Checkbox("Hist", &c.history)) markDirty(); if (ImGui::IsItemHovered()) ImGui::SetTooltip("History mode (WCONHIST): rates are observed data"); ImGui::SameLine(); }
        if (ImGui::SmallButton("\xC3\x97")) del = ci;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Remove this control");
        // targets
        const Units U = Units::of(case_.units);
        char lb[64];
        // targets/limits in a 2 x (label, value) grid that fits the panel width
        auto fieldD = [&](const char* label, const char* id, double* v, bool active) {
            ImGui::TableNextColumn();
            if (active) ImGui::TextColored(col4(Pal::Accent), "%s", label); else ImGui::TextDisabled("%s", label);
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (active) ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.16f, 0.27f, 0.45f, 1));
            const bool chg = ImGui::InputDouble(id, v, 0, 0, "%.6g");
            if (active) ImGui::PopStyleColor();
            if (chg) markDirty();
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(active ? "Target of the active control mode" : "Limit (negative = none)");
        };
        if (ImGui::BeginTable("targets", 4, ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("l1", ImGuiTableColumnFlags_WidthStretch, 0.8f);
            ImGui::TableSetupColumn("v1", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("l2", ImGuiTableColumnFlags_WidthStretch, 0.8f);
            ImGui::TableSetupColumn("v2", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            char lb2[64];
            if (prod) {
                std::snprintf(lb, sizeof(lb), "Oil (%s/d)", U.liq);   fieldD(lb, "##o", &c.orat, c.mode == Control::ORAT);
                std::snprintf(lb2, sizeof(lb2), "Water (%s/d)", U.liq); fieldD(lb2, "##w", &c.wrat, c.mode == Control::WRAT);
                std::snprintf(lb, sizeof(lb), "Gas (%s/d)", U.gas);   fieldD(lb, "##g", &c.grat, c.mode == Control::GRAT);
                std::snprintf(lb2, sizeof(lb2), "Liquid (%s/d)", U.liq); fieldD(lb2, "##l", &c.lrat, c.mode == Control::LRAT);
                fieldD("Voidage (RESV)", "##r", &c.resv, c.mode == Control::RESV);
                std::snprintf(lb, sizeof(lb), "BHP (%s)", U.p);       fieldD(lb, "##b", &c.bhp, c.mode == Control::BHP);
            } else {
                std::snprintf(lb, sizeof(lb), "Rate (%s/d)", w && w->type == WellType::GasInjector ? U.gas : U.liq);
                fieldD(lb, "##rate", &c.rate, c.mode == Control::RATE);
                std::snprintf(lb2, sizeof(lb2), "BHP (%s)", U.p);    fieldD(lb2, "##b", &c.bhp, c.mode == Control::BHP);
            }
            ImGui::EndTable();
        }
        // ---- OPM Flow-only control: group / THP control, VFP table, artificial lift
        {
            const bool opmSet = !c.deckMode.empty() || c.thp > 0 || c.vfp > 0 || c.alq > 0;
            if (opmSet) ImGui::SetNextItemOpen(true, ImGuiCond_Once);
            if (ImGui::TreeNodeEx("OPM Flow control (group / THP / VFP)", ImGuiTreeNodeFlags_SpanAvailWidth)) {
                static const char* prodModes[] = {"(same as above)", "GRUP", "THP", "NONE", "CRAT"};
                static const char* injModes[]  = {"(same as above)", "GRUP", "THP", "RESV"};
                const char** modes = prod ? prodModes : injModes;
                const int nModes = prod ? 5 : 4;
                const std::string dm = c.deckMode == "RESV*" ? std::string("RESV") : c.deckMode;
                int cur = 0; for (int m = 1; m < nModes; ++m) if (dm == modes[m]) cur = m;
                ImGui::SetNextItemWidth(170);
                if (ImGui::Combo("Mode sent to OPM Flow", &cur, modes, nModes)) {
                    c.deckMode = cur == 0 ? std::string() : std::string(modes[cur]);
                    if (c.deckMode == "RESV" && c.rate <= 0) c.deckMode = "RESV*";
                    markDirty();
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("OPM Flow runs this control mode (GRUP = under its group's GCONPROD/GCONINJE, THP = tubing-head\n"
                                      "pressure through the VFP table). The ResSim engine does not model groups or VFP and runs\n"
                                      "the mode chosen above, using the targets/limits you give.");
                if (ImGui::BeginTable("opmctl", 4, ImGuiTableFlags_SizingStretchProp)) {
                    ImGui::TableSetupColumn("l1", ImGuiTableColumnFlags_WidthStretch, 0.8f);
                    ImGui::TableSetupColumn("v1", ImGuiTableColumnFlags_WidthStretch, 1.0f);
                    ImGui::TableSetupColumn("l2", ImGuiTableColumnFlags_WidthStretch, 0.8f);
                    ImGui::TableSetupColumn("v2", ImGuiTableColumnFlags_WidthStretch, 1.0f);
                    char l3[64]; std::snprintf(l3, sizeof(l3), "THP (%s)", U.p);
                    fieldD(l3, "##thp", &c.thp, dm == "THP");
                    ImGui::TableNextColumn(); ImGui::TextDisabled("VFP table");
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::InputInt("##vfp", &c.vfp, 0)) { c.vfp = std::max(0, c.vfp); markDirty(); }
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("VFPPROD / VFPINJ table number (0 = none). Define the tables on the Advanced keywords page.");
                    if (prod) fieldD("Artificial lift (ALQ)", "##alq", &c.alq, false);
                    ImGui::EndTable();
                }
                if ((dm == "THP" || c.thp > 0) && c.vfp <= 0)
                    ImGui::TextColored(col4(Pal::Warn), "A THP control/limit needs a VFP table number.");
                ImGui::TreePop();
            }
        }
        ImGui::EndChild();
        ImGui::PopID();
    }
    if (del >= 0) { st.controls.erase(st.controls.begin() + del); markDirty(); }
    ImGui::BeginDisabled(case_.wells.empty());
    if (ImGui::Button("+ Add well control")) {
        WellControl c; c.well = case_.wells[std::min(selWell_, int(case_.wells.size()) - 1)].name;
        const WellInput* w = nullptr; for (auto& x : case_.wells) if (x.name == c.well) w = &x;
        if (w && w->type != WellType::Producer) { c.mode = Control::RATE; c.rate = 0; } else { c.mode = Control::ORAT; c.orat = 0; }
        st.controls.push_back(c);
        markDirty();
    }
    ImGui::EndDisabled();
}

// ============================================================================ NUMERICS
void App::editNumerics() {
    Tuning& t = case_.tuning;
    bool ch = false;
    ui::header("Time stepping");
    if (ui::beginProps("ts")) {
        ch |= ui::propD("Initial step (days)", &t.initialDt, "%.4g");
        ch |= ui::propD("Maximum step (days)", &t.maxDt, "%.4g");
        ch |= ui::propD("Minimum step (days)", &t.minDt, "%.3g");
        ch |= ui::propD("Growth factor", &t.growth, "%.3g");
        ch |= ui::propD("Cut factor after failure", &t.cutFactor, "%.3g");
        ch |= ui::propD("Max saturation change / iter", &t.dsMax, "%.3g", "Appleyard-type chop of Newton updates");
        ch |= ui::propD("Max relative pressure change / iter", &t.dpMaxRel, "%.3g");
        ui::endProps();
    }
    ui::header("Newton and linear solver");
    if (ui::beginProps("nl")) {
        ch |= ui::propI("Max Newton iterations", &t.maxNewton, 2, 50);
        ch |= ui::propD("CNV tolerance (max local)", &t.cnvTol, "%.3g");
        ch |= ui::propD("Material balance tolerance", &t.mbTol, "%.3g");
        ch |= ui::propI("Max linear iterations", &t.maxLinear, 10, 5000);
        ch |= ui::propD("Linear relative tolerance", &t.linTol, "%.3g");
            static const char* ls[] = {"Block ILU(0) - BiCGStab", "CPR-AMG + block ILU(0) - FGMRES",
                                "Automatic (block ILU, CPR when needed)", "GPU (D3D11 compute, mixed precision)"};
    ch |= ui::propCombo("Linear solver", &t.linearSolver, ls, 4,
                        "All options use 3x3 block storage with the wells folded in by a Schur complement.\n"
                        "CPR: constrained pressure residual - true-IMPES pressure equation solved by aggregation AMG (over-corrected "
                        "V-cycles inside BiCGStab), then block ILU(0) on the full system, inside flexible GMRES.\n"
                        "Automatic: block ILU(0) first (cheapest), CPR for the next solves when it fails.\n"
                        "GPU: block-Jacobi BiCGStab in single precision on the graphics card as the preconditioner of a double-"
                        "precision FGMRES on the CPU; falls back to CPR when the device is unavailable or the solve fails.");
        static const char* cw[] = {"Quasi-IMPES", "True-IMPES"};
        ch |= ui::propCombo("CPR pressure weights", &t.cprWeights, cw, 2,
                            "How CPR forms the pressure equation. True-IMPES (default, as in OPM Flow) combines the mass balances with "
                            "weights from the accumulation terms only - more robust for strongly coupled or gravity-dominated flow. "
                            "Quasi-IMPES uses the full diagonal block of the Jacobian.");
        const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
        ch |= ui::propI("CPU threads (0 = all)", &t.threads, 0, hw);
        ui::endProps();
    }
    if (ch) markDirty();
    if (ImGui::Button("Restore defaults")) { t = Tuning{}; markDirty(); }
    ui::hint("Fully-implicit black-oil: Newton-Raphson with automatic-differentiation Jacobians, BiCGSTAB/GMRES with CPR-AMG or ILU(0), "
             "OpenMP-parallel assembly and sparse kernels. Convergence: CNV (max scaled residual) and material balance per component.");
}
