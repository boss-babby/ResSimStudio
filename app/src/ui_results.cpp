// 3D viewer panel, result plots and run monitor.
#include "app.h"
#include "icons.h"
#include "widgets.h"
#include "implot.h"
#include "plotkit.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>

using namespace ressim;

// ----------------------------------------------------------------------------- property catalogue
namespace {
enum PropKind { Static, Initial, Dynamic };
struct PropDef { const char* name; const char* label; PropKind kind; };
const PropDef kProps[] = {
    {"PORO", "Porosity", Static},          {"PERMX", "Permeability X", Static}, {"PERMZ", "Permeability Z", Static},
    {"NTG", "Net-to-gross", Static},       {"DEPTH", "Cell depth", Static},     {"PORV", "Pore volume", Static},
    {"PRESSURE_0", "Initial pressure", Initial}, {"SWAT_0", "Initial water saturation", Initial},
    {"SGAS_0", "Initial gas saturation", Initial},
    {"PRESSURE", "Pressure", Dynamic},     {"SWAT", "Water saturation", Dynamic}, {"SOIL", "Oil saturation", Dynamic},
    {"SGAS", "Gas saturation", Dynamic},   {"RS", "Solution GOR", Dynamic},
    {"RV", "Vaporised oil-gas ratio", Dynamic},
    {"POLYMER", "Polymer concentration", Dynamic}, {"PADS", "Adsorbed polymer", Dynamic},
    {"TOF", "Time of flight (streamlines)", Dynamic}, {"DRAINAGE", "Drainage region (producer)", Dynamic},
    {"SWEEP", "Sweep region (injector)", Dynamic},
};
bool isStreamlineProp(const char* n) { const std::string s = n; return s == "TOF" || s == "DRAINAGE" || s == "SWEEP"; }
const int kNumProps = int(sizeof(kProps) / sizeof(kProps[0]));

const char* unitOf(const char* name, UnitSystem us) {
    const Units U = Units::of(us);
    const std::string n = name;
    if (n.rfind("PERM", 0) == 0) return "mD";
    if (n == "DEPTH") return U.len;
    if (n == "PORV") return us == UnitSystem::Field ? "rb" : "rm3";
    if (n.rfind("PRESSURE", 0) == 0) return U.p;
    if (n == "RS") return U.rs;
    if (n == "RV") return U.rv;
    if (n == "POLYMER") return us == UnitSystem::Field ? "lb/stb" : "kg/sm3";
    if (n == "PADS") return us == UnitSystem::Field ? "lb/lb" : "kg/kg";
    if (n == "TOF") return "days";
    return "";
}
bool isSaturation(const char* name) { const std::string n = name; return n.find("SWAT") == 0 || n.find("SGAS") == 0 || n == "SOIL"; }

// summary vectors
enum Vec { FOPR, FWPR, FGPR, FWIR, FGIR, FOPT, FWPT, FGPT, FWIT, FGIT, FPR, FWCT, FGOR, FOIP, FWIP, FGIP,
           FCIR, FCPR, FCIT, FCPT, FCIP, FCAD, DT, NEWTON, LINEAR,
           WOPR, WWPR, WGPR, WWIR, WGIR, WBHP, WWCT, WGOR, WCIR, WCPR, WOPRH, WWPRH, WGPRH, NumVec };
const char* kVecNames[NumVec] = {"FOPR", "FWPR", "FGPR", "FWIR", "FGIR", "FOPT", "FWPT", "FGPT", "FWIT", "FGIT", "FPR", "FWCT",
                                 "FGOR", "FOIP", "FWIP", "FGIP", "FCIR", "FCPR", "FCIT", "FCPT", "FCIP", "FCAD",
                                 "DT", "NEWTON", "LINEAR", "WOPR", "WWPR", "WGPR", "WWIR",
                                 "WGIR", "WBHP", "WWCT", "WGOR", "WCIR", "WCPR", "WOPRH", "WWPRH", "WGPRH"};
double vecValue(const SummaryPoint& s, int v, int well) {
    const WellRates* w = (well >= 0 && well < int(s.wells.size())) ? &s.wells[well] : nullptr;
    switch (v) {
        case FOPR: return s.fopr; case FWPR: return s.fwpr; case FGPR: return s.fgpr; case FWIR: return s.fwir; case FGIR: return s.fgir;
        case FOPT: return s.fopt; case FWPT: return s.fwpt; case FGPT: return s.fgpt; case FWIT: return s.fwit; case FGIT: return s.fgit;
        case FPR: return s.fpr; case FWCT: return s.fwct; case FGOR: return s.fgor;
        case FOIP: return s.oip; case FWIP: return s.wip; case FGIP: return s.gip;
        case FCIR: return s.fcir; case FCPR: return s.fcpr; case FCIT: return s.fcit; case FCPT: return s.fcpt;
        case FCIP: return s.fcip; case FCAD: return s.fcad;
        case DT: return s.dt; case NEWTON: return s.newton; case LINEAR: return s.linear;
        default: break;
    }
    if (!w) return NAN;
    switch (v) {
        case WOPR: return w->oil; case WWPR: return w->water; case WGPR: return w->gas; case WWIR: return w->waterInj;
        case WGIR: return w->gasInj; case WBHP: return w->open ? w->bhp : NAN;
        case WWCT: { const double l = w->oil + w->water; return l > 0 ? w->water / l : 0; }
        case WGOR: return w->oil > 0 ? w->gas / w->oil : 0;
        case WCIR: return w->polyInj; case WCPR: return w->polyProd;
        case WOPRH: return w->oilH >= 0 ? w->oilH : NAN;
        case WWPRH: return w->waterH >= 0 ? w->waterH : NAN;
        case WGPRH: return w->gasH >= 0 ? w->gasH : NAN;
    }
    return NAN;
}
struct Series { const std::vector<SummaryPoint>* s; int vec; int well; double xs; int skipFirst; };
ImPlotPoint seriesGetter(int idx, void* data) {
    const Series* d = static_cast<const Series*>(data);
    const SummaryPoint& p = (*d->s)[idx + d->skipFirst];
    return ImPlotPoint(p.time * d->xs, vecValue(p, d->vec, d->well));
}
void plotVec(const char* label, const std::vector<SummaryPoint>& s, int vec, int well, double xs, ImU32 col, float weight = 2.0f,
             bool scatter = false) {
    // rates are step averages: skip the t=0 point for rate/ratio vectors
    const bool isRate = (vec <= FGIR) || vec == FCIR || vec == FCPR || vec == FWCT || vec == FGOR || (vec >= WOPR && vec != WBHP) || vec == DT || vec == NEWTON || vec == LINEAR;
    Series d{&s, vec, well, xs, (isRate && s.size() > 1) ? 1 : 0};
    const int n = int(s.size()) - d.skipFirst;
    if (n <= 0) return;
    ImPlotSpec spec;
    spec.LineColor = col4(col);
    spec.LineWeight = weight * g_plotLineScale;
    if (fig::recorder().recording()) {
        std::vector<double> X(static_cast<size_t>(n)), Y(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) { const ImPlotPoint q = seriesGetter(i, &d); X[size_t(i)] = q.x; Y[size_t(i)] = q.y; }
        recSeries(label, X.data(), Y.data(), n, col, weight, scatter);
    }
    if (scatter) {
        spec.Marker = ImPlotMarker_Circle; spec.MarkerSize = 3.0f;
        spec.MarkerFillColor = col4(col); spec.MarkerLineColor = col4(col);
        ImPlot::PlotScatterG(label, seriesGetter, &d, n, spec);
    } else ImPlot::PlotLineG(label, seriesGetter, &d, n, spec);
}
ImU32 wellColor(int i) { return seriesColor(i); }
}

std::vector<LineVertex> streamlineVertices(const StreamlineResult& S, int colorBy, Colormap cm, const std::string& focus, bool focusOnly);

// ----------------------------------------------------------------------------- view data
void App::updateViewData() {
    if (!viewOk_) return;
    const PropDef& pd = kProps[std::min(prop_, kNumProps - 1)];
    // choose grid source
    std::shared_ptr<const Grid> results;
    std::vector<float> vals;
    const Grid* grid = nullptr;
    const Case* src = &case_;
    unsigned runVer = runner_.version();
    std::string caption;
    {
        std::lock_guard<std::mutex> lock(runner_.mutex());
        const RunResults& R = runner_.results();
        if (pd.kind == Dynamic && R.valid && !R.snaps.empty()) {
            results = R.grid;
            grid = results.get();
            src = nullptr;
        }
        ViewKey k;
        k.prop = prop_; k.snap = snap_; k.version = caseVersion_; k.runVer = (pd.kind == Dynamic) ? runVer : 0;
        k.grid = grid ? (const void*)grid : (pd.kind == Initial && init_.valid ? (const void*)init_.grid.get() : (const void*)preview_.get());
        k.autoR = autoRange_; k.a = rmin_; k.b = rmax_; k.log = view_.settings().logScale;
        if (k.prop == key_.prop && k.snap == key_.snap && k.version == key_.version && k.runVer == key_.runVer && k.grid == key_.grid &&
            k.autoR == key_.autoR && k.a == key_.a && k.b == key_.b && k.log == key_.log)
            return;
        key_ = k;
        if (pd.kind == Dynamic) {
            if (!grid) { vals.clear(); }
            else {
                snap_ = std::max(0, std::min(snap_, int(R.snaps.size()) - 1));
                const Snapshot& s = R.snaps[snap_];
                const std::string n = pd.name;
                if (isStreamlineProp(pd.name)) {
                    vals.assign(s.pressure.size(), NAN);
                    if (sl_.res && sl_.snap == snap_) {
                        const StreamlineResult& S = *sl_.res;
                        for (size_t a = 0; a < vals.size() && a < S.cellTof.size(); ++a) {
                            if (n == "TOF") vals[a] = S.cellTof[a] >= 0 ? S.cellTof[a] : NAN;
                            else if (n == "DRAINAGE") vals[a] = S.cellSink[a] >= 0 ? float(S.cellSink[a]) : NAN;
                            else vals[a] = S.cellSource[a] >= 0 ? float(S.cellSource[a]) : NAN;
                        }
                    }
                } else
                vals = n == "PRESSURE" ? s.pressure : n == "SWAT" ? s.sw : n == "SOIL" ? s.so : n == "SGAS" ? s.sg : n == "RV" ? s.rv : n == "POLYMER" ? s.poly : n == "PADS" ? s.polyAds : s.rs;
                if (vals.size() != size_t(grid->numActive())) vals.assign(size_t(grid->numActive()), NAN);
                char b[128];
                std::snprintf(b, sizeof(b), "t = %.1f days (%.2f years)  |  report %d of %d", s.time, s.time / 365.25, snap_,
                              int(R.snaps.size()) - 1);
                caption = b;
                if (resultsCaseVersion_ != caseVersion_) caption += "  |  inputs changed since this run";
                if (grid != viewGrid_) {
                    view_.setGrid(grid);
                    viewHold_ = results;
                    viewGrid_ = grid;
                }
                view_.setWells(wellGlyphs(*grid, R.input, std::max(0, snap_ - 1)));
            }
        }
    }
    if (pd.kind == Initial) {
        if (init_.valid) {
            grid = init_.grid.get();
            const std::string n = pd.name;
            vals = n == "PRESSURE_0" ? init_.p : n == "SWAT_0" ? init_.sw : init_.sg;
            caption = "Initial state (equilibrium)";
            if (init_.caseVersion != caseVersion_) caption += "  |  inputs changed - recompute";
            if (grid != viewGrid_) {
                view_.setGrid(grid);
                viewHold_ = init_.grid;
                viewGrid_ = grid;
            }
            view_.setWells(wellGlyphs(*grid, case_));
        }
    }
    if (!grid) {    // static or fallback: model preview
        grid = preview_.get();
        if (grid && pd.kind == Static && grid->numCells() == (src ? src->grid : case_.grid).size()) {
            const GridInput& gi = src ? src->grid : case_.grid;
            const int na = grid->numActive();
            vals.resize(na);
            const std::string n = pd.name;
            for (int a = 0; a < na; ++a) {
                const int q = grid->globalOf[a];
                double v = 0;
                if (n == "PORO") v = gi.poro[q];
                else if (n == "PERMX") v = gi.permx[q];
                else if (n == "PERMZ") v = gi.permz.empty() ? gi.permx[q] : gi.permz[q];
                else if (n == "NTG") v = gi.ntg.empty() ? 1.0 : gi.ntg[q];
                else if (n == "DEPTH") v = grid->depth[a];
                else if (n == "PORV") v = grid->poreVolRef[a];
                vals[a] = float(v);
            }
            caption = "Model input (" + std::to_string(grid->nx) + " x " + std::to_string(grid->ny) + " x " + std::to_string(grid->nz) + ")";
        } else if (pd.kind != Static && vals.empty()) {
            caption = pd.kind == Dynamic ? "No results yet - press Run (F5)" : "Press 'Compute initial state' (Initialization)";
        }
        if (grid && grid != viewGrid_) {
            view_.setGrid(grid);
            viewHold_ = preview_;
            view_.setWells(wellGlyphs(*grid, case_));
            viewGrid_ = grid;
        } else if (grid) view_.setWells(wellGlyphs(*grid, case_));
    }
    // range
    float lo = FLT_MAX, hi = -FLT_MAX;
    for (float v : vals) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
    if (lo > hi) { lo = 0; hi = 1; }
    if (isSaturation(pd.name)) { lo = 0; hi = 1; }
    if (hi - lo < 1e-9f * std::max(1.0f, std::fabs(hi))) { hi = lo + std::max(1e-6f, 0.01f * std::fabs(lo)); }
    if (autoRange_) { rmin_ = lo; rmax_ = hi; key_.a = rmin_; key_.b = rmax_; }
    view_.setValues(vals, rmin_, rmax_, std::string(pd.name) + " - " + pd.label, unitOf(pd.name, case_.units));
    view_.setCaption(caption);
}

// ----------------------------------------------------------------------------- streamlines
void App::updateStreamlines() {
    ViewSettings& vs = view_.settings();
    const bool needProp = isStreamlineProp(kProps[prop_].name);
    vs.showCells = !(sl_.show && sl_.hideCells);
    vs.linesXray = sl_.xray;
    auto clearLines = [&]() { if (sl_.builtLines != -1) { view_.setLines({}); sl_.builtLines = -1; } };
    if (!sl_.show && !needProp) { clearLines(); return; }
    if (runner_.running()) return;
    std::shared_ptr<Simulator> sim; std::shared_ptr<const Grid> grid; int nsnap = 0;
    const unsigned ver = runner_.version();
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        sim = R.sim; grid = R.grid; nsnap = int(R.snaps.size());
    }
    if (!sim || nsnap == 0) { clearLines(); sl_.res.reset(); return; }
    const int s = std::clamp(snap_, 0, nsnap - 1);
    bool rebuilt = false;
    static int tracedFor = -1;
    if (!sl_.res || sl_.snap != s || sl_.runVer != ver || tracedFor != sl_.lines) {
        const auto t0 = std::chrono::steady_clock::now();
        try {
            sim->restoreReport(s);
            StreamlineOptions o; o.lines = sl_.lines;
            sl_.res = std::make_shared<StreamlineResult>(traceStreamlines(*sim, o));
        } catch (const std::exception& e) { runner_.addLog(std::string("Streamline tracing failed: ") + e.what()); sl_.res = std::make_shared<StreamlineResult>(); }
        sl_.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        sl_.snap = s; sl_.runVer = ver; sl_.grid = grid; tracedFor = sl_.lines;
        key_.prop = -1;   // refresh cell colours (TOF / regions)
        rebuilt = true;
    }
    if (!sl_.show) { clearLines(); return; }
    const auto& W = view_.wells();
    const std::string focus = vs.focusWell >= 0 && vs.focusWell < int(W.size()) ? W[size_t(vs.focusWell)].name : std::string();
    const std::string fkey = focus + (sl_.focusOnly ? "|1" : "|0");
    if (!rebuilt && sl_.builtLines == int(sl_.res->lines.size()) && sl_.builtColor == sl_.colorBy + 10 * int(vs.cmap) && sl_.builtFocus == fkey) return;
    const StreamlineResult& S = *sl_.res;
    std::vector<LineVertex> v = streamlineVertices(S, sl_.colorBy, vs.cmap, focus, sl_.focusOnly);
    sl_.builtFocus = fkey;
    view_.setLines(std::move(v));
    sl_.builtLines = int(S.lines.size());
    sl_.builtColor = sl_.colorBy + 10 * int(vs.cmap);
}

// ----------------------------------------------------------------------------- 3D viewer panel
void App::panelViewer() {
    ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("3D Viewer")) { ImGui::End(); return; }
    if (!viewOk_) { ImGui::TextColored(col4(Pal::Err), "3D renderer failed to initialise (Direct3D 11 / shader compiler)."); ImGui::End(); return; }

    // ---- toolbar row
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 13);
    if (ImGui::BeginCombo("##prop", kProps[prop_].label)) {
        const char* groups[3] = {"Static (input)", "Initial state", "Dynamic (results)"};
        for (int gk = 0; gk < 3; ++gk) {
            ImGui::SeparatorText(groups[gk]);
            for (int i = 0; i < kNumProps; ++i) {
                if (kProps[i].kind != gk) continue;
                if (std::string(kProps[i].name) == "RV" && !case_.phases.vapoil) continue;
                if ((std::string(kProps[i].name) == "POLYMER" || std::string(kProps[i].name) == "PADS") && !case_.phases.polymer) continue;
                char b[96]; std::snprintf(b, sizeof(b), "%-11s %s", kProps[i].name, kProps[i].label);
                if (ImGui::Selectable(b, prop_ == i)) { prop_ = i; autoRange_ = true; }
            }
        }
        ImGui::EndCombo();
    }
    int nsnap = 0;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        nsnap = int(runner_.results().snaps.size());
    }
    if (kProps[prop_].kind == Dynamic && nsnap > 0) {
        ImGui::SameLine();
        if (ImGui::Button(play_ ? ICON_PAUSE "##play" : ICON_PLAY "##play")) { play_ = !play_; if (play_ && snap_ >= nsnap - 1) snap_ = 0; }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(play_ ? "Pause animation" : "Play animation through report steps");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        ImGui::SliderInt("##snap", &snap_, 0, std::max(0, nsnap - 1), "report %d");
        if (play_) {
            playAcc_ += ImGui::GetIO().DeltaTime;
            if (playAcc_ > 0.12) { playAcc_ = 0; if (++snap_ >= nsnap) { snap_ = nsnap - 1; play_ = false; } }
        }
    }
    ViewerBinding vb;
    vb.view = &view_; vb.autoRange = &autoRange_; vb.rmin = &rmin_; vb.rmax = &rmax_; vb.sl = &sl_;
    vb.id = "model"; vb.propName = kProps[prop_].name; vb.units = unitOf(kProps[prop_].name, case_.units);
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        vb.slAvailable = R.sim && !R.snaps.empty() && !runner_.running();
        vb.slHint = runner_.running() ? "Simulation running..." : R.valid && R.engine != "ResSim" ? "External engine results: open them in the Results workspace (streamlines are traced there from the restart file)."
                                                                : "Run the case with the ResSim engine first (F5).";
    }
    vb.grid = static_cast<const Grid*>(viewGrid_);
    if (!vb.grid) vb.grid = preview_.get();
    vb.pngName = case_.title.empty() ? std::string("model_view") : case_.title + "_" + kProps[prop_].name;
    if (kProps[prop_].kind == Dynamic && !isStreamlineProp(kProps[prop_].name)) {
        vb.cellHistory = [this](int cell, std::vector<double>& t, std::vector<double>& v) {
            const std::string n = kProps[prop_].name;
            std::lock_guard<std::mutex> l(runner_.mutex());
            for (const Snapshot& s : runner_.results().snaps) {
                const std::vector<float>& a = n == "PRESSURE" ? s.pressure : n == "SWAT" ? s.sw : n == "SOIL" ? s.so : n == "SGAS" ? s.sg : n == "RV" ? s.rv : n == "POLYMER" ? s.poly : n == "PADS" ? s.polyAds : s.rs;
                if (size_t(cell) >= a.size()) continue;
                t.push_back(s.time); v.push_back(a[size_t(cell)]);
            }
            return !t.empty();
        };
    }
    viewerToolbar(vb, mvUi_);

    updateStreamlines();
    updateViewData();
    vb.grid = static_cast<const Grid*>(viewGrid_);
    viewerCanvas(vb, mvUi_, [this](int cell) {
        const Grid* g = static_cast<const Grid*>(viewGrid_);
        int i, j, k; g->ijk(g->globalOf[size_t(cell)], i, j, k);
        ImGui::Text("Cell (%d, %d, %d)", i + 1, j + 1, k + 1);
        ImGui::Separator();
        ImGui::Text("%s = %.5g %s", kProps[prop_].name, view_.valueOf(cell), unitOf(kProps[prop_].name, case_.units));
        ImGui::TextDisabled("depth %.1f %s, PV %.4g", g->depth[size_t(cell)], Units::of(case_.units).len, g->poreVolRef[size_t(cell)]);
    });
    ImGui::End();
}

// ----------------------------------------------------------------------------- plots
void App::panelPlots() {
    ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Plots")) { ImGui::End(); return; }
    plotsBody();
    ImGui::End();
}

void App::plotsBody() {
    static bool years = false;
    static int customSel[NumVec] = {};
    std::lock_guard<std::mutex> lock(runner_.mutex());
    const RunResults& R = runner_.results();
    if (!R.valid || R.summary.empty()) {
        ImGui::Dummy(ImVec2(0, 40));
        ui::hint("No simulation results yet. Build your model in the Inspector, then press  \xE2\x96\xB6 Run  (F5).");
        return;
    }
    const auto& S = R.summary;
    const double xs = years ? 1.0 / 365.25 : 1.0;
    const Units U = Units::of(R.input.units);
    char liq[32], gas[32], pr[32];
    std::snprintf(liq, sizeof(liq), "%s/d", U.liq); std::snprintf(gas, sizeof(gas), "%s/d", U.gas); std::snprintf(pr, sizeof(pr), "%s", U.p);
    const bool hasGas = R.input.phases.gas, hasWater = R.input.phases.water;
    double tSnap = -1;
    if (snap_ >= 0 && snap_ < int(R.snaps.size())) tSnap = R.snaps[snap_].time * xs;

    ImGui::Checkbox("Time in years", &years);
    ImGui::SameLine();
    ImGui::TextDisabled("  Double-click a plot to move the 3D time cursor there.");
    static int tabNow = 0;
    static const char* tabFig[] = {"run.field", "run.well", "run.custom", "run.conv"};
    ImGui::SameLine();
    exportButton(tabFig[tabNow]);
    auto cursor = [&]() {
        if (tSnap >= 0) {
            ImPlotSpec sp; sp.LineColor = ImVec4(1, 1, 1, 0.35f); sp.LineWeight = 1.0f; sp.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##cursor", &tSnap, 1, sp);
            recVLine(tSnap);
        }
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(0) && !R.snaps.empty()) {
            const double t = ImPlot::GetPlotMousePos().x / xs;
            int best = 0;
            for (int i = 0; i < int(R.snaps.size()); ++i) if (std::fabs(R.snaps[i].time - t) < std::fabs(R.snaps[best].time - t)) best = i;
            snap_ = best;
            if (kProps[prop_].kind != Dynamic) prop_ = 9;   // PRESSURE
        }
    };
    const char* xl = years ? "Time (years)" : "Time (days)";
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;

    if (ImGui::BeginTabBar("plots")) {
        if (ImGui::BeginTabItem("Field")) {
            tabNow = 0;
            const ImVec2 av = ImGui::GetContentRegionAvail();
            const float fs = ImGui::GetFontSize();
            const PlotGrid G = fitGrid(4, av, fs * 22, fs * 11, 1.6f);
            if (G.scroll) ImGui::BeginChild("##fieldscroll", av);
            fig::recorder().beginFigure("run.field", G.rows, G.cols);
            if (ImPlot::BeginSubplots("##field", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                if (ImPlot::BeginPlot("Production / injection rates")) {
                    recPanel("Production / injection rates", xl, liq, hasGas ? gas : "");
                    ImPlot::SetupAxes(xl, liq, af, af); setupSI(ImAxis_Y1);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    plotVec("Oil", S, FOPR, -1, xs, Pal::Oil);
                    if (hasWater) { plotVec("Water", S, FWPR, -1, xs, Pal::Water); plotVec("Water inj", S, FWIR, -1, xs, IM_COL32(120, 190, 255, 255), 1.5f); }
                    if (hasGas) { useY2(true); plotVec("Gas", S, FGPR, -1, xs, Pal::Gas); if (R.input.wells.size()) plotVec("Gas inj", S, FGIR, -1, xs, IM_COL32(255, 150, 140, 255), 1.5f); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Cumulative volumes")) {
                    recPanel("Cumulative volumes", xl, U.liq, hasGas ? U.gas : "");
                    ImPlot::SetupAxes(xl, U.liq, af, af); setupSI(ImAxis_Y1);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    plotVec("FOPT", S, FOPT, -1, xs, Pal::Oil);
                    if (hasWater) { plotVec("FWPT", S, FWPT, -1, xs, Pal::Water); plotVec("FWIT", S, FWIT, -1, xs, IM_COL32(120, 190, 255, 255), 1.5f); }
                    if (hasGas) { useY2(true); plotVec("FGPT", S, FGPT, -1, xs, Pal::Gas); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Pressure")) {
                    recPanel("Pressure", xl, pr, "");
                    ImPlot::SetupAxes(xl, pr, af, af); setupSI(ImAxis_Y1);
                    plotVec("Average (FPR)", S, FPR, -1, xs, Pal::Pressure, 2.5f);
                    for (int w = 0; w < int(R.wellNames.size()); ++w) {
                        const std::string l = "BHP " + R.wellNames[w];
                        plotVec(l.c_str(), S, WBHP, w, xs, wellColor(w), 1.5f);
                    }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Water cut / GOR")) {
                    recPanel("Water cut / GOR", xl, "Water cut (-)", hasGas ? U.rs : "");
                    ImPlot::SetupAxes(xl, "Water cut (-)", af, 0);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.rs, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    if (hasWater) plotVec("FWCT", S, FWCT, -1, xs, Pal::Water);
                    if (hasGas) { useY2(true); plotVec("FGOR", S, FGOR, -1, xs, Pal::Gas); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            if (G.scroll) ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Wells")) {
            tabNow = 1;
            const int nw = int(R.wellNames.size());
            if (nw == 0) ui::hint("The run has no wells.");
            else {
                resultsWell_ = std::min(resultsWell_, nw - 1);
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                if (ImGui::BeginCombo("Well", R.wellNames[resultsWell_].c_str())) {
                    for (int w = 0; w < nw; ++w) if (ImGui::Selectable(R.wellNames[w].c_str(), w == resultsWell_)) resultsWell_ = w;
                    ImGui::EndCombo();
                }
                const int w = resultsWell_;
                const bool prod = R.wellTypes[w] == WellType::Producer;
                ImGui::SameLine();
                ImGui::TextDisabled("%s", toString(R.wellTypes[w]));
                if (!S.empty() && w < int(S.back().wells.size())) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("  |  control now: %s", toString(S.back().wells[w].mode));
                }
                const PlotGrid G = fitGrid(3, ImGui::GetContentRegionAvail(), ImGui::GetFontSize() * 22, ImGui::GetFontSize() * 11, 1.5f);
                if (G.scroll) ImGui::BeginChild("##wellscroll", ImGui::GetContentRegionAvail());
                fig::recorder().beginFigure("run.well", G.rows, G.cols);
                if (ImPlot::BeginSubplots("##well", G.rows, G.cols, ImVec2(G.scroll ? -ImGui::GetStyle().ScrollbarSize : -1, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                    if (ImPlot::BeginPlot(prod ? "Rates (lines = simulated, dots = observed)" : "Injection rate")) {
                        recPanel(prod ? "Rates" : "Injection rate", xl, liq, prod && hasGas ? gas : "");
                        ImPlot::SetupAxes(xl, liq, af, af); setupSI(ImAxis_Y1);
                        if (prod) {
                            if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                            plotVec("WOPR", S, WOPR, w, xs, Pal::Oil);
                            plotVec("WOPRH", S, WOPRH, w, xs, Pal::Oil, 1, true);
                            if (hasWater) { plotVec("WWPR", S, WWPR, w, xs, Pal::Water); plotVec("WWPRH", S, WWPRH, w, xs, Pal::Water, 1, true); }
                            if (hasGas) { useY2(true); plotVec("WGPR", S, WGPR, w, xs, Pal::Gas); plotVec("WGPRH", S, WGPRH, w, xs, Pal::Gas, 1, true); useY2(false); }
                        } else {
                            if (R.wellTypes[w] == WellType::WaterInjector) plotVec("WWIR", S, WWIR, w, xs, Pal::Water);
                            else plotVec("WGIR", S, WGIR, w, xs, Pal::Gas);
                        }
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot("Bottom-hole pressure")) {
                        recPanel("Bottom-hole pressure", xl, pr, "");
                        ImPlot::SetupAxes(xl, pr, af, af); setupSI(ImAxis_Y1);
                        plotVec("WBHP", S, WBHP, w, xs, wellColor(w));
                        plotVec("FPR", S, FPR, -1, xs, Pal::Pressure, 1.0f);
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot("Water cut / GOR")) {
                        recPanel("Water cut / GOR", xl, "WWCT", hasGas ? U.rs : "");
                        ImPlot::SetupAxes(xl, "WWCT", af, 0);
                        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                        if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.rs, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                        if (prod && hasWater) plotVec("WWCT", S, WWCT, w, xs, Pal::Water);
                        if (prod && hasGas) { useY2(true); plotVec("WGOR", S, WGOR, w, xs, Pal::Gas); useY2(false); }
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
                fig::recorder().endFigure();
                if (G.scroll) ImGui::EndChild();
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Custom")) {
            tabNow = 2;
            ImGui::BeginChild("vecs", ImVec2(ImGui::GetFontSize() * 13, 0), ImGuiChildFlags_Borders);
            ImGui::TextDisabled("Field vectors");
            for (int v = 0; v <= FCAD; ++v) { if (v >= FCIR && !case_.phases.polymer) continue; bool b = customSel[v] != 0; if (ImGui::Checkbox(kVecNames[v], &b)) customSel[v] = b; }
            ImGui::TextDisabled("Solver");
            for (int v = DT; v <= LINEAR; ++v) { bool b = customSel[v] != 0; if (ImGui::Checkbox(kVecNames[v], &b)) customSel[v] = b; }
            ImGui::TextDisabled("Well vectors (well: %s)", R.wellNames.empty() ? "-" : R.wellNames[std::min(resultsWell_, int(R.wellNames.size()) - 1)].c_str());
            for (int v = WOPR; v < NumVec; ++v) { if ((v == WCIR || v == WCPR) && !case_.phases.polymer) continue; bool b = customSel[v] != 0; if (ImGui::Checkbox(kVecNames[v], &b)) customSel[v] = b; }
            ImGui::EndChild();
            ImGui::SameLine();
            fig::recorder().beginFigure("run.custom", 1, 1);
            if (ImPlot::BeginPlot("##custom", ImVec2(-1, -1))) {
                recPanel("", xl, "", "");
                ImPlot::SetupAxes(xl, nullptr, af, af); setupSI(ImAxis_Y1);
                int ci = 0;
                for (int v = 0; v < NumVec; ++v)
                    if (customSel[v]) plotVec(kVecNames[v], S, v, R.wellNames.empty() ? -1 : std::min(resultsWell_, int(R.wellNames.size()) - 1), xs, wellColor(ci++), 2.0f, v >= WOPRH);
                cursor();
                recEnd();
                ImPlot::EndPlot();
            }
            fig::recorder().endFigure();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Convergence")) {
            tabNow = 3;
            fig::recorder().beginFigure("run.conv", 1, 2);
            if (ImPlot::BeginSubplots("##conv", 1, 2, ImVec2(-1, std::max(ImGui::GetContentRegionAvail().y, ImGui::GetFontSize() * 13.0f)), ImPlotSubplotFlags_LinkAllX)) {
                if (ImPlot::BeginPlot("Time-step size")) {
                    recPanel("Time-step size", xl, "dt (days)", "");
                    ImPlot::SetupAxes(xl, "dt (days)", af, af); setupSI(ImAxis_Y1);
                    plotVec("dt", S, DT, -1, xs, Pal::Accent, 1.5f, true);
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("Iterations per time step")) {
                    recPanel("Iterations per time step", xl, "Newton", "Linear");
                    ImPlot::SetupAxes(xl, "Newton", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, "Linear", af | ImPlotAxisFlags_AuxDefault);
                    plotVec("Newton", S, NEWTON, -1, xs, Pal::Warn, 1.5f, true);
                    useY2(true);
                    plotVec("Linear", S, LINEAR, -1, xs, Pal::Muted, 1.0f);
                    useY2(false);
                    recEnd();
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
}

// ----------------------------------------------------------------------------- run monitor
static std::string bigNum(double v) {
    char b[64];
    const double a = std::fabs(v);
    if (a >= 1e9) std::snprintf(b, sizeof(b), "%.3f B", v / 1e9);
    else if (a >= 1e6) std::snprintf(b, sizeof(b), "%.3f MM", v / 1e6);
    else if (a >= 1e3) std::snprintf(b, sizeof(b), "%.2f M", v / 1e3);
    else std::snprintf(b, sizeof(b), "%.4g", v);
    return b;
}

void App::panelRun() {
    if (!ImGui::Begin("Run Monitor")) { ImGui::End(); return; }
    const bool running = runner_.running();
    const double prog = runner_.progress();
    char ov[96];
    std::snprintf(ov, sizeof(ov), running ? "Running  %.1f%%" : (prog >= 0.999 ? "Finished" : "Idle  %.0f%%"), prog * 100);
    ImGui::PushStyleColor(ImGuiCol_PlotHistogram, running ? col4(Pal::Accent) : col4(Pal::Ok));
    ImGui::ProgressBar(float(prog), ImVec2(-1, 0), ov);
    ImGui::PopStyleColor();
    const std::string err = runner_.error();
    if (!err.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextColored(col4(Pal::Err), "%s", err.c_str()); ImGui::PopTextWrapPos(); }

    std::lock_guard<std::mutex> lock(runner_.mutex());
    const RunResults& R = runner_.results();
    if (!R.valid) { ui::hint("Run statistics, current rates and recovery appear here during a simulation."); ImGui::End(); return; }
    const Units U = Units::of(R.input.units);
    const SummaryPoint& s = R.summary.back();
    if (ImGui::BeginTable("kpi", 4, ImGuiTableFlags_SizingStretchSame)) {
        auto kpi = [&](const char* label, const std::string& value, ImU32 col) {
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label);
            ImGui::PushFont(g_fonts.title ? g_fonts.title : g_fonts.ui, ImGui::GetStyle().FontSizeBase * 1.25f);
            ImGui::TextColored(col4(col), "%s", value.c_str());
            ImGui::PopFont();
        };
        char b[64];
        std::snprintf(b, sizeof(b), "%.0f / %.0f d", s.time, R.totalTime); kpi("Simulated time", b, IM_COL32(230, 232, 240, 255));
        std::snprintf(b, sizeof(b), "%.1f", s.fpr); kpi((std::string("Avg pressure (") + U.p + ")").c_str(), b, Pal::Pressure);
        std::snprintf(b, sizeof(b), "%.1f %%", R.initial.oil > 0 ? 100.0 * s.fopt / R.initial.oil : 0.0);
        kpi("Oil recovery factor", b, Pal::Oil);
        std::snprintf(b, sizeof(b), "%.2f s", runner_.wallSeconds()); kpi("Wall-clock time", b, IM_COL32(230, 232, 240, 255));
        ImGui::EndTable();
    }
    ImGui::Separator();
    if (ImGui::BeginTable("rates", 2, ImGuiTableFlags_SizingStretchProp)) {
        auto row = [&](const char* k, const std::string& v) { ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextDisabled("%s", k); ImGui::TableNextColumn(); ImGui::TextUnformatted(v.c_str()); };
        char b[96];
        std::snprintf(b, sizeof(b), "%.1f %s/d", s.fopr, U.liq); row("Oil rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d   (WCT %.3f)", s.fwpr, U.liq, s.fwct); row("Water rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d   (GOR %.4g)", s.fgpr, U.gas, s.fgor); row("Gas rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d", s.fwir, U.liq); row("Water injection", b);
        row("Cumulative oil", bigNum(s.fopt) + " " + U.liq);
        row("STOIIP (initial)", bigNum(R.initial.oil + R.initial.oilVaporized) + " " + U.liq);
        if (R.initial.oilVaporized > 0) row("  of which vaporised in gas", bigNum(R.initial.oilVaporized) + " " + U.liq);
        row("Gas in place (initial)", bigNum(R.initial.gasFree + R.initial.gasDissolved) + " " + U.gas);
        std::snprintf(b, sizeof(b), "%d steps, %d cuts, %d Newton, %d linear", R.stats.timesteps, R.stats.cuts, R.stats.newton, R.stats.linear);
        row("Solver", b);
        ImGui::EndTable();
    }
    ImGui::End();
}
