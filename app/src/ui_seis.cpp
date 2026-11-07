// C20: Tools > Geomechanics: seismicity-aware injection design (rate-and-state seismicity + convex optimisation).
//  * fault patches (Excel-style table): cell, friction, stress-path factor, A sigma, t_a, background rate, weight
//  * designed injectors, volumes, control intervals, horizon, objective (expected events / peak rate / blend)
//  * background run: base + one step-response run per well (Green's functions) -> global convex optimum -> nonlinear
//    verification runs -> model-corrected re-optimisation (sequential convex programming)
//  * results: schedules, Coulomb stress, seismicity rate, cumulative induced events, model-correction history;
//    open the designed case, copy the SCHEDULE text (WCONINJE / TSTEP)
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "implot.h"
#include "plotkit.h"
#include "imgui_internal.h"
#include "ressim/seismicity.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <memory>
#include <mutex>
#include <thread>

using namespace ressim;

namespace {
struct SeisState {
    bool init = false;
    std::vector<seis::Patch> patches;
    std::vector<std::string> injectors;          // of the current case
    std::vector<char> use;                       // designed?
    std::vector<double> wellVol;                 // per well (case surface volume)
    bool perWell = false;
    double totalVolume = 0;
    int nCtrl = 10, subSteps = 2, postSteps = 16, corrIters = 3;
    double injDays = 30, postDays = 60, maxRateFactor = 3, blend = 0.5;
    double lastRateFactor = 0;                   // C22 (H1f): 0 = no extra bound on the last interval
    int objective = 0;
    int patchSel = -1;
    // background job
    std::shared_ptr<std::mutex> mx = std::make_shared<std::mutex>();
    std::shared_ptr<std::atomic<bool>> busy = std::make_shared<std::atomic<bool>>(false);
    std::shared_ptr<std::atomic<bool>> cancel = std::make_shared<std::atomic<bool>>(false);
    std::shared_ptr<std::vector<std::string>> log = std::make_shared<std::vector<std::string>>();
    std::shared_ptr<seis::DesignOutput> result;
    std::shared_ptr<seis::DesignSetup> resultSetup;
    std::shared_ptr<Case> resultBase;
    unsigned caseStamp = 0;
    // C23 (research Cycle 5 / F14): kernel calibration from observed stage-wise event counts
    struct StageRow { int test = 1; double V = 0, E = 0, q = 0, dur = 0; int n = 0; };
    std::vector<StageRow> stages;
    std::shared_ptr<seis::KernelFit> kfit;
    bool showKernel = false;
};
SeisState g_seis;

// Aspo HRL 2015, tests HF1 / HF2 (conventional, Avro granodiorite): stage volume (L), hydraulic energy (kJ), mean rate (L/min),
// duration (s), triggered AE count incl. 600 s after shut-in (Zang et al. 2024 ESSD data, research/exp/h6_aspo.py)
void loadAspoStages(std::vector<SeisState::StageRow>& v) {
    const double r[][6] = {{1, 0.77, 7.9, 0.96, 47, 0}, {1, 2.66, 24.3, 1.29, 123, 3}, {1, 4.40, 35.8, 1.58, 166, 10}, {1, 5.14, 42.2, 2.41, 127, 12},
                           {1, 4.67, 41.1, 4.12, 67, 3}, {1, 9.25, 81.5, 4.44, 124, 20}, {2, 2.96, 26.7, 1.67, 105, 4}, {2, 5.03, 44.4, 1.63, 184, 6},
                           {2, 1.22, 10.1, 2.36, 30, 0}, {2, 5.19, 47.8, 2.45, 126, 17}, {2, 4.84, 46.8, 4.61, 62, 15}, {2, 9.42, 101.5, 4.75, 118, 50}};
    v.clear();
    for (const auto& q : r) { SeisState::StageRow s; s.test = int(q[0]); s.V = q[1]; s.E = q[2]; s.q = q[3]; s.dur = q[4]; s.n = int(q[5]); v.push_back(s); }
}

void line(const char* label, const std::vector<double>& x, const std::vector<double>& y, ImU32 col, float w = 2.0f, bool dashed = false) {
    if (x.empty() || x.size() != y.size()) return;
    ImPlotSpec sp;
    if (col) sp.LineColor = col4(col);
    sp.LineWeight = w * g_plotLineScale;
    if (dashed) sp.Flags = ImPlotLineFlags_Segments;
    ImPlot::PlotLine(label, x.data(), y.data(), int(x.size()), sp);
}
void stairs(const char* label, const std::vector<double>& x, const std::vector<double>& y, ImU32 col, float w = 2.0f) {
    if (x.size() < 2 || y.empty()) return;
    std::vector<double> xs(x), ys(y);
    ys.push_back(y.back());            // ImPlot stairs: value held until the next x
    xs.resize(ys.size());
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = w * g_plotLineScale;
    ImPlot::PlotStairs(label, xs.data(), ys.data(), int(xs.size()), sp);
}
} // namespace

void App::seismicityWindow() {
    if (!showSeis_) return;
    SeisState& S = g_seis;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 74, ImGui::GetFontSize() * 48), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
    if (seisDock_) {
        seisDock_ = false;
        if (workspace_ == Workspace::Model)
            if (ImGuiDockNode* cn = ImGui::DockBuilderGetCentralNode(ImGui::GetID("MainDock"))) ImGui::SetNextWindowDockID(cn->ID, ImGuiCond_Always);
        ImGui::SetNextWindowFocus();
    }
    if (!ImGui::Begin(ICON_ACTIVITY "  Geomechanics: seismicity-aware injection design###seiswin", &showSeis_)) { ImGui::End(); return; }
    const bool field = case_.units == UnitSystem::Field;
    const char* rateU = field ? "STB/d" : "sm3/d";
    const char* volU = field ? "STB" : "sm3";
    // ---- injectors of the current case (refresh when the well list changes)
    {
        std::vector<std::string> inj;
        for (const auto& w : case_.wells) if (w.type != WellType::Producer) inj.push_back(w.name);
        if (inj != S.injectors || !S.init) {
            std::vector<std::string> old = S.injectors; std::vector<char> oldUse = S.use; std::vector<double> oldVol = S.wellVol;
            S.injectors = inj; S.use.assign(inj.size(), 1); S.wellVol.assign(inj.size(), 0.0);
            for (size_t i = 0; i < inj.size(); ++i)
                for (size_t k = 0; k < old.size(); ++k) if (old[k] == inj[i]) { S.use[i] = oldUse[k]; S.wellVol[i] = oldVol[k]; }
            // default volumes: first rate control of each injector x injection window
            for (size_t i = 0; i < inj.size(); ++i) {
                if (S.wellVol[i] > 0) continue;
                for (const auto& st : case_.schedule) {
                    bool found = false;
                    for (const auto& wc : st.controls) if (wc.well == inj[i] && wc.rate > 0) { S.wellVol[i] = wc.rate * S.injDays; found = true; break; }
                    if (found) break;
                }
            }
            if (!(S.totalVolume > 0)) { S.totalVolume = 0; for (size_t i = 0; i < inj.size(); ++i) S.totalVolume += S.wellVol[i]; }
        }
        if (!S.init) {
            S.init = true;
            if (S.patches.empty()) {   // default patch: a few cells from the first injector
                seis::Patch p; p.name = "F1";
                for (const auto& w : case_.wells)
                    if (w.type != WellType::Producer) {
                        p.i = std::min(w.i + 5, std::max(0, case_.grid.nx - 1)); p.j = w.j;
                        p.k = w.completions.empty() ? 0 : w.completions.front().k1;
                        break;
                    }
                S.patches.push_back(p);
            }
        }
    }
    const bool busy = S.busy->load();
    if (ImGui::BeginTable("seislay", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("in", ImGuiTableColumnFlags_WidthStretch, 0.38f);
        ImGui::TableSetupColumn("out", ImGuiTableColumnFlags_WidthStretch, 0.62f);
        ImGui::TableNextColumn();
        ImGui::BeginChild("seisin", ImVec2(0, 0));
        ImGui::BeginDisabled(busy);
        ui::header("Fault patches");
        ui::hint("Each patch is one grid cell on a fault. dCFS = friction x stress-path factor x pore-pressure change. "
                 "A sigma and t_a are the rate-and-state constitutive parameters (Dieterich 1994); background rate in events/year.");
        {
            ui::RecordTable<seis::Patch> t;
            t.text("Name", [](seis::Patch& p) -> std::string& { return p.name; })
             .integer("I", [](seis::Patch& p) -> int& { return p.i; }, 1)
             .integer("J", [](seis::Patch& p) -> int& { return p.j; }, 1)
             .integer("K", [](seis::Patch& p) -> int& { return p.k; }, 1)
             .num("Friction", [](seis::Patch& p) -> double& { return p.friction; })
             .num("Stress-path factor", [](seis::Patch& p) -> double& { return p.cfsFactor; })
             .num("A sigma (MPa)", [](seis::Patch& p) -> double& { return p.aSigma; })
             .num("t_a (years)", [](seis::Patch& p) -> double& { return p.taYears; })
             .num("Background (1/yr)", [](seis::Patch& p) -> double& { return p.rate; })
             .num("Weight", [](seis::Patch& p) -> double& { return p.weight; });
            t.draw("seispatch", S.patches, ImGui::GetFontSize() * 9);
        }
        ui::header("Designed injectors");
        if (S.injectors.empty()) ui::hint("The case has no injector. Add a water or gas injector in Wells.");
        if (ui::beginProps("seisw", 0.55f)) {
            ui::propBool("Volume per well", &S.perWell, "Off: one total volume shared optimally between the selected wells (allocation is optimised too)");
            if (!S.perWell) ui::propD((std::string("Total volume (") + volU + ")").c_str(), &S.totalVolume, "%.6g");
            ui::endProps();
        }
        for (size_t i = 0; i < S.injectors.size(); ++i) {
            ImGui::PushID(int(i));
            bool u = S.use[i] != 0;
            if (ImGui::Checkbox(S.injectors[i].c_str(), &u)) S.use[i] = u;
            if (S.perWell && u) {
                ImGui::SameLine(ImGui::GetFontSize() * 9);
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9);
                ImGui::InputDouble((std::string("##v") + volU).c_str(), &S.wellVol[i], 0, 0, "%.6g");
                ImGui::SameLine(); ImGui::TextDisabled("%s", volU);
            }
            ImGui::PopID();
        }
        ui::header("Schedule and objective");
        if (ui::beginProps("seiss", 0.55f)) {
            ui::propD("Injection window (days)", &S.injDays, "%.5g");
            ui::propD("Post shut-in horizon (days)", &S.postDays, "%.5g", "Events after shut-in are counted: rankings change when they are ignored");
            ui::propI("Control intervals", &S.nCtrl, 1, 200);
            ui::propI("Report steps per interval", &S.subSteps, 1, 20);
            ui::propI("Report steps after shut-in", &S.postSteps, 1, 200);
            ui::propD("Max rate / mean rate", &S.maxRateFactor, "%.3g");
            ui::propD("Last-interval rate / mean (0 = off)", &S.lastRateFactor, "%.3g",
                      "Research finding F12 (H1f): the event-count optimum always ends with a burst at the maximum rate - it lowers the total\n"
                      "count but raises the stress at shut-in, where large post-injection events are a known risk (Basel, Pohang).\n"
                      "A value <= 1 removes the burst (the price is shown as 'without a shut-in burst' in the design potential).");
            static const char* objs[] = {"Expected induced events", "Peak seismicity rate", "Blend"};
            ui::propCombo("Objective", &S.objective, objs, 3, "Expected events: exact convex objective. Peak: smooth maximum of dCFS/(A sigma) (= peak rate when t << t_a).");
            if (S.objective == 2) ui::propD("Peak weight (0..1)", &S.blend, "%.3f");
            ui::propI("Model-correction passes", &S.corrIters, 0, 10, "Sequential convex programming: each pass adds the simulated-minus-predicted stress of the last verified schedule and re-optimises (1 run)");
            ui::endProps();
        }
        ImGui::EndDisabled();
        int nUse = 0; for (char u : S.use) nUse += u;
        const std::string valid = validation();
        if (!busy) {
            ImGui::BeginDisabled(nUse == 0 || S.patches.empty() || !valid.empty() || case_.comp.enabled || case_.thermal.enabled);
            if (ui::accentButton(ICON_PLAY " Design injection")) {
                auto setup = std::make_shared<seis::DesignSetup>();
                for (size_t i = 0; i < S.injectors.size(); ++i) if (S.use[i]) { setup->wells.push_back(S.injectors[i]); if (S.perWell) setup->wellVolume.push_back(S.wellVol[i]); }
                setup->totalVolume = S.totalVolume;
                setup->patches = S.patches;
                setup->nCtrl = S.nCtrl; setup->injDays = S.injDays; setup->postDays = S.postDays; setup->subSteps = S.subSteps;
                setup->postSteps = S.postSteps; setup->maxRateFactor = S.maxRateFactor; setup->correctionIters = S.corrIters;
                setup->lastRateFactor = std::max(0.0, S.lastRateFactor);
                setup->opt.obj = S.objective == 0 ? seis::Objective::Events : S.objective == 1 ? seis::Objective::Peak : seis::Objective::Blend;
                setup->opt.blend = S.blend;
                auto base = std::make_shared<Case>(case_);
                S.cancel->store(false); S.busy->store(true);
                { std::lock_guard<std::mutex> l(*S.mx); S.log->clear(); S.result.reset(); }
                auto mx = S.mx; auto lg = S.log; auto bz = S.busy; auto cn = S.cancel;
                std::thread([mx, lg, bz, cn, setup, base]() {
                    auto out = std::make_shared<seis::DesignOutput>(seis::designInjection(*base, *setup, cn.get(),
                        [mx, lg](const std::string& m) { std::lock_guard<std::mutex> l(*mx); lg->push_back(m); }));
                    { std::lock_guard<std::mutex> l(*mx); g_seis.result = out; g_seis.resultSetup = setup; g_seis.resultBase = base; }
                    bz->store(false);
                }).detach();
            }
            ImGui::EndDisabled();
            if (case_.comp.enabled || case_.thermal.enabled) ui::hint("The designer uses the black-oil engine (water / gas injection); compositional and thermal cases are not supported yet.");
            else if (!valid.empty()) ui::hint(("Case not valid: " + valid).c_str());
        } else {
            if (ImGui::Button(ICON_SQUARE " Cancel")) S.cancel->store(true);
            ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "running ...");
        }
        // ---- C23: which kernel do the observed events follow? (stage-wise counts -> Poisson / AIC)
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("Calibrate from observed stage events (does the rate matter?)", S.showKernel ? ImGuiTreeNodeFlags_DefaultOpen : 0)) {
            S.showKernel = true;
            ui::hint("One row per injection stage of past stimulations: volume, hydraulic energy (int p q dt), mean rate, duration and the "
                     "number of located events. Each test / interval gets its own rate factor; the fit ranks constant, duration, volume, energy "
                     "and volume x rate^beta x (1 + earlier volume)^gamma kernels by AIC and tells whether a lower rate reduces events at fixed "
                     "volume - the premise of the schedule design above. Research Cycle 5: at Aspo HRL it does not (finding F14).");
            ui::RecordTable<SeisState::StageRow> t;
            t.integer("Test", [](SeisState::StageRow& r) -> int& { return r.test; })
             .num("Volume", [](SeisState::StageRow& r) -> double& { return r.V; })
             .num("Energy", [](SeisState::StageRow& r) -> double& { return r.E; })
             .num("Mean rate", [](SeisState::StageRow& r) -> double& { return r.q; })
             .num("Duration (s)", [](SeisState::StageRow& r) -> double& { return r.dur; })
             .integer("Events", [](SeisState::StageRow& r) -> int& { return r.n; });
            t.draw("seisstages", S.stages, ImGui::GetFontSize() * 9);
            if (ImGui::Button("Aspo HF1 + HF2 example")) { loadAspoStages(S.stages); S.kfit.reset(); }
            ImGui::SameLine();
            ImGui::BeginDisabled(S.stages.size() < 3);
            if (ui::accentButton("Fit kernels")) {
                std::vector<seis::StageObs> obs;
                for (const auto& r : S.stages) { seis::StageObs o; o.test = r.test; o.volume = r.V; o.energy = r.E; o.rate = r.q; o.duration = r.dur; o.count = std::max(0, r.n); obs.push_back(o); }
                S.kfit = std::make_shared<seis::KernelFit>(seis::fitStageKernels(obs));
            }
            ImGui::EndDisabled();
            if (S.kfit) {
                const seis::KernelFit& F = *S.kfit;
                if (ui::beginTable("seiskf", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    for (const char* c : {"Kernel", "dAIC", "Exponents"}) ui::tableColumn(c);
                    ImGui::TableHeadersRow();
                    for (const auto& m : F.models) {
                        ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(m.name.c_str());
                        ImGui::TableNextColumn(); ImGui::Text("%.1f", m.aic - F.models.front().aic);
                        ImGui::TableNextColumn();
                        if (m.exponent != 0 || m.exponent2 != 0) ImGui::Text("%.2f  %.2f", m.exponent, m.exponent2);
                    }
                    ui::endTable();
                }
                ImGui::TextColored(col4(F.rateLo > 0 ? Pal::Ok : Pal::Warn), "%s", "Verdict:");
                ImGui::SameLine(); ImGui::TextWrapped("%s", F.verdict.c_str());
            }
        } else S.showKernel = false;
        {
            std::lock_guard<std::mutex> l(*S.mx);
            if (!S.log->empty()) {
                ui::header("Log");
                ImGui::BeginChild("seislog", ImVec2(0, ImGui::GetFontSize() * 9), ImGuiChildFlags_Borders);
                for (const auto& s : *S.log) ImGui::TextWrapped("%s", s.c_str());
                if (busy) ImGui::SetScrollHereY(1.0f);
                ImGui::EndChild();
            }
        }
        ImGui::EndChild();
        // ---- results
        ImGui::TableNextColumn();
        ImGui::BeginChild("seisout", ImVec2(0, 0));
        std::shared_ptr<seis::DesignOutput> R; std::shared_ptr<seis::DesignSetup> RS; std::shared_ptr<Case> RB;
        { std::lock_guard<std::mutex> l(*S.mx); R = S.result; RS = S.resultSetup; RB = S.resultBase; }
        if (!R) {
            ui::hint("Define fault patches and injectors, then press Design injection.\n\n"
                     "Method (research/papers/P1): with rate-and-state seismicity and pore-pressure Coulomb stress, the expected "
                     "number of induced events is a convex function of the injection rates, so the optimum is global. "
                     "Green's functions come from step-response runs of this model; the optimum is verified by a full forward "
                     "run and the model is corrected until prediction and simulation agree.");
        } else if (!R->ok) {
            ImGui::TextColored(col4(Pal::Err), "Design failed: %s", R->error.c_str());
        } else {
            if (ui::beginTable("seissum", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                for (const char* c : {"", "Constant rate", "Designed", "Change"}) ui::tableColumn(c);
                ImGui::TableHeadersRow();
                auto row = [&](const char* n, double a, double b) {
                    ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(n);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", a);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", b);
                    ImGui::TableNextColumn();
                    const double ch = std::fabs(a) > 0 ? 100 * (b / a - 1) : 0;
                    ImGui::TextColored(col4(ch <= 0 ? Pal::Ok : Pal::Warn), "%+.1f %%", ch);
                };
                row("Induced events, simulated", R->simEventsConst, R->simEventsOpt);
                row("Induced events, design model", R->eventsConst, R->eventsOpt);
                row("Peak rate / background", R->peakConst, R->peakOpt);
                ui::endTable();
            }
            ImGui::TextDisabled("Model error: first pass %.2f %%, returned schedule %.2f %% after %d correction pass(es); %d simulations, %d optimiser iterations, %.1f s",
                                100 * R->linearityError, 100 * R->finalModelError, R->corrections, R->runs, R->iters, R->seconds);
            if (R->Pi > 0) {   // C22: research H1f / F12 - how much can a schedule gain here at all?
                ImGui::TextColored(col4(Pal::Accent), "Design potential (first-order reference):");
                ImGui::SameLine();
                ImGui::Text("Pi = %.3g, theta = %.3g  ->  %s%.1f %% fewer events (%.1f %% without a shut-in burst)%s", R->Pi, R->theta,
                            R->gainTheoryIsBound ? "at least " : "up to ",
                            100 * R->gainTheory, 100 * R->gainTheoryNoBurst, R->farField ? "  [far-field patch: not applicable]" : "");
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
                    ImGui::SetTooltip("Pi = max Coulomb change / (A sigma) of the constant schedule; theta = stress response time / injection window.\n"
                                      "theta << 0.1: the stress follows the rate - constant rate is already near optimal, but avoid rate spikes.\n"
                                      "theta ~ 0.1-1: the largest design benefit.\nThis is the exact optimum for a first-order (saturating) stress response. Diffusive responses (radial/linear flow)\n"
                                      "reach clearly less (research H1g: about 1/4 to 3/4 of it) - the designed schedule above is the real number.");
                if (R->endsWithBurst)
                    ui::hint("The designed schedule ends at its rate bound (shut-in burst). It minimises the expected count but raises the stress at "
                             "shut-in; set 'Last-interval rate / mean' to 1 to avoid it.");
            }
            if (ImGui::Button(ICON_SPARKLES " Open designed case") && RB && RS) {
                if (maybeSaveChanges()) {
                    try { pendingCase_ = std::make_unique<Case>(seis::designCase(*RB, *RS, R->qOpt)); newCase(4); notify("Designed case opened - press Run (F5)"); }
                    catch (const std::exception& e) { notify(std::string("Could not build the case: ") + e.what()); }
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Copy SCHEDULE text") && RB && RS) { ImGui::SetClipboardText(seis::designScheduleDeck(*RB, *RS, R->qOpt).c_str()); notify("WCONINJE / TSTEP schedule copied to the clipboard"); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("ECLIPSE / OPM Flow SCHEDULE text: one WCONINJE + TSTEP per control interval, then shut-in");
            if (ImGui::BeginTabBar("seistabs")) {
                const float ph = std::max(ImGui::GetFontSize() * 15, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 1.4f);
                const size_t np = R->simOpt.size();
                if (ui::tabItem("Injection schedule")) {
                    if (ImPlot::BeginPlot("Injection rate##seisq", ImVec2(-1, ph))) {
                        ImPlot::SetupAxes("time (days)", (std::string("rate (") + rateU + ")").c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        for (size_t w = 0; w < R->qOpt.size(); ++w) {
                            const std::string nm = RS ? RS->wells[w] : std::to_string(w + 1);
                            stairs((nm + " constant").c_str(), R->ctrlStart, R->qConst[w], Pal::Muted, 1.5f);
                            stairs((nm + " designed").c_str(), R->ctrlStart, R->qOpt[w], seriesColor(int(w)), 2.5f);
                        }
                        ImPlot::EndPlot();
                    }
                    ImGui::EndTabItem();
                }
                auto perPatch = [&](const char* tab, const char* ylab, int what) {
                    if (!ui::tabItem(tab)) return;
                    if (ImPlot::BeginPlot((std::string(tab) + "##seisp").c_str(), ImVec2(-1, ph))) {
                        ImPlot::SetupAxes("time (days)", ylab, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        if (what == 1) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                        for (size_t p = 0; p < np; ++p) {
                            const std::string nm = RS ? RS->patches[p].name : std::to_string(p + 1);
                            for (int k = 0; k < 2; ++k) {
                                const seis::Forecast& f = k == 0 ? R->simConst[p] : R->simOpt[p];
                                std::vector<double> y;
                                if (what == 0) y = f.dS;
                                else if (what == 1) y = f.R;
                                else { y = f.N; const double r = (RS ? RS->patches[p].rate : 1) / 365.25; for (size_t i = 0; i < y.size(); ++i) y[i] -= r * (f.t[i] - f.t[0]); }
                                line((nm + (k == 0 ? " constant" : " designed")).c_str(), f.t, y, k == 0 ? Pal::Muted : seriesColor(int(p)), k == 0 ? 1.5f : 2.5f);
                            }
                        }
                        if (RS) { const double x = RS->injDays; ImPlot::PlotInfLines("shut-in", &x, 1); }
                        ImPlot::EndPlot();
                    }
                    ImGui::EndTabItem();
                };
                perPatch("Coulomb stress change", "dCFS (MPa)", 0);
                perPatch("Seismicity rate", "R (events / year)", 1);
                perPatch("Induced events", "cumulative excess events", 2);
                if (ui::tabItem("Model correction")) {
                    if (ImPlot::BeginPlot("Design model vs nonlinear simulation##seisc", ImVec2(-1, ph))) {
                        ImPlot::SetupAxes("pass", "induced events", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        std::vector<double> x;
                        for (size_t i = 0; i < R->corrPred.size(); ++i) x.push_back(double(i));
                        ImPlotSpec a; a.Marker = ImPlotMarker_Circle; a.LineColor = col4(Pal::Accent); a.MarkerFillColor = a.LineColor;
                        ImPlot::PlotLine("design model", x.data(), R->corrPred.data(), int(x.size()), a);
                        ImPlotSpec b2; b2.Marker = ImPlotMarker_Square; b2.LineColor = col4(Pal::Warn); b2.MarkerFillColor = b2.LineColor;
                        ImPlot::PlotLine("simulated", x.data(), R->corrSim.data(), int(x.size()), b2);
                        ImPlot::EndPlot();
                    }
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }
        }
        ImGui::EndChild();
        ImGui::EndTable();
    }
    ImGui::End();
}
