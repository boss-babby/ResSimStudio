// C14: results files (.rsres), ECLIPSE result export and "restart from a report step"
//  * File > Open ResSim results / Save results: a whole run (summary, every report-step map, the geothermal / EDFM
//    analysis, the case that was run and the engine restart states) in one binary file - reopens without re-running.
//  * File > Export results as ECLIPSE files: EGRID / INIT / SMSPEC / UNSMRY / UNRST for ResInsight, OPM tools, resdata.
//  * Run > Restart from a report step: continues the current results (or a loaded .rsres) from a saved engine state
//    with the current model - typically with a changed schedule after that step (what-if forecasts without re-running
//    the history). The continuation is exact: a restart without changes reproduces the continuous run.
#include "app.h"
#include "platform.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "results_io.h"
#include "opm_engine.h"
#include "ressim/io.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;

namespace {
std::string defaultResultsPath(const std::string& projectPath, const std::string& title) {
    if (!projectPath.empty()) { fs::path p = fs::u8path(projectPath); return (p.parent_path() / (p.stem().u8string() + ".rsres")).u8string(); }
    std::string s;
    for (char ch : title) s += std::isalnum(static_cast<unsigned char>(ch)) ? ch : '_';
    if (s.empty()) s = "results";
    if (s.size() > 40) s.resize(40);
    return s + ".rsres";
}
const char* kindName(const std::string& k) { return k == "BO" ? "black oil" : k == "COMP" ? "compositional" : k == "THERMAL" ? "thermal" : k.c_str(); }
} // namespace

void App::openResultsFile(const std::string& path) {
    if (runner_.running()) { notify("Stop the running simulation first"); return; }
    if (!maybeSaveChanges()) return;
    try {
        RunResults r = loadRunResults(path);
        r.sourceFile = path;
        const std::string v = validateCase(r.input);
        runner_.cancel(); runner_.join();
        case_ = r.input;
        path_.clear();                       // the stored case is not a project file: Save asks for a name
        dirty_ = false;
        init_ = InitPreview{};
        if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
        sl_.res.reset();
        syncBuildersFromCase();
        ++caseVersion_;
        undoReset();
        gridDirty_ = true;
        rebuildPreview();
        selWell_ = selStep_ = resultsWell_ = 0;
        const int ns = int(r.snaps.size());
        const size_t nr = r.restarts.size();
        const std::string eng = r.engine;
        runner_.setResults(std::move(r));
        resultsCaseVersion_ = caseVersion_;
        prop_ = 9; snap_ = std::max(0, ns - 1); play_ = false; viewGrid_ = nullptr;
        if (viewOk_) view_.resetCamera();
        workspace_ = Workspace::Results;
        restartBase_.reset(); restartBaseFile_.clear();
        runner_.addLog("=== Opened results " + path + " (" + eng + ", " + std::to_string(ns) + " report steps, " + std::to_string(nr) + " restart states)");
        notify("Opened results " + fs::u8path(path).filename().u8string() + (v.empty() ? std::string() : " - model warnings: " + v));
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "Open results failed", MB_OK | MB_ICONERROR);
    }
}

void App::saveResultsDialog() {
    std::string p;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        p = runner_.results().sourceFile;
    }
    if (p.empty()) p = defaultResultsPath(path_, case_.title);
    if (!fileDialog(true, "ResSim results (*.rsres)\0*.rsres\0All files (*.*)\0*.*\0", "rsres", p)) return;
    try {
        std::lock_guard<std::mutex> l(runner_.mutex());
        saveRunResults(runner_.results(), p);
        runner_.results().sourceFile = p;
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "Save results failed", MB_OK | MB_ICONERROR);
        return;
    }
    std::error_code ec;
    const auto sz = fs::file_size(fs::u8path(p), ec);
    char b[64]; std::snprintf(b, sizeof b, " (%.1f MB)", ec ? 0.0 : double(sz) / 1048576.0);
    notify("Saved results " + fs::u8path(p).filename().u8string() + b);
}

void App::exportResultsEclipseDialog() {
    std::string p = runDirFor(path_, case_.title, "ECLOUT");
    {
        std::error_code ec; fs::create_directories(fs::u8path(p), ec);
        p = (fs::u8path(p) / (caseBaseName(case_) + ".EGRID")).u8string();
    }
    if (!fileDialog(true, "ECLIPSE grid + results (*.EGRID)\0*.EGRID\0All files (*.*)\0*.*\0", "EGRID", p)) return;
    fs::path base = fs::u8path(p); base.replace_extension();
    std::vector<std::string> files, warn;
    std::string err;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        err = exportRunEclipse(runner_.results(), base.u8string(), &files, &warn);
    }
    if (!err.empty()) { plat::msgBox(err.c_str(), "ECLIPSE export failed", MB_OK | MB_ICONERROR); return; }
    for (const auto& w : warn) runner_.addLog("NOTE (ECLIPSE export): " + w);
    runner_.addLog("=== ECLIPSE result files written: " + base.u8string() + ".EGRID / .INIT / .SMSPEC / .UNSMRY / .UNRST");
    notify("Wrote " + std::to_string(files.size()) + " ECLIPSE result files (" + base.filename().u8string() + ")");
}

void App::afterRunOutputs() {
    if (!settings_.autoSaveResults && !settings_.eclOutput) return;
    std::lock_guard<std::mutex> l(runner_.mutex());
    RunResults& R = runner_.results();
    if (!R.valid || R.engine.rfind("ResSim", 0) != 0) return;
    if (settings_.autoSaveResults) {
        const std::string p = path_.empty() ? (fs::u8path(runDirFor(path_, case_.title, "RESULTS")) / (caseBaseName(case_) + ".rsres")).u8string()
                                            : defaultResultsPath(path_, case_.title);
        try {
            std::error_code ec; fs::create_directories(fs::u8path(p).parent_path(), ec);
            saveRunResults(R, p);
            R.sourceFile = p;
            runner_.log().push_back("=== Results saved: " + p);
        } catch (const std::exception& e) { runner_.log().push_back(std::string("WARN: results not saved: ") + e.what()); }
    }
    if (settings_.eclOutput) {
        const std::string dir = runDirFor(path_, case_.title, "ECLOUT");
        std::error_code ec; fs::create_directories(fs::u8path(dir), ec);
        const std::string base = (fs::u8path(dir) / caseBaseName(case_)).u8string();
        std::vector<std::string> files, warn;
        const std::string err = exportRunEclipse(R, base, &files, &warn);
        if (!err.empty()) runner_.log().push_back("WARN: ECLIPSE result files not written: " + err);
        else {
            for (const auto& w : warn) runner_.log().push_back("NOTE (ECLIPSE output): " + w);
            runner_.log().push_back("=== ECLIPSE result files: " + base + ".EGRID / .INIT / .SMSPEC / .UNSMRY / .UNRST");
        }
    }
}

void App::restartSettingsUi(bool& ch) {
    ui::header("Results and restart (built-in engine)");
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
    if (ImGui::InputInt("Save the restart state every N report steps (0 = never)", &settings_.rstEvery)) { settings_.rstEvery = std::clamp(settings_.rstEvery, 0, 100000); ch = true; }
    ch |= ImGui::Checkbox("Save the restart state at the end of every run", &settings_.rstLast);
    ch |= ImGui::Checkbox("Save the results file (.rsres) after every run", &settings_.autoSaveResults);
    ch |= ImGui::Checkbox("Also write ECLIPSE result files (EGRID / INIT / UNRST / SMSPEC / UNSMRY)", &settings_.eclOutput);
    ui::hint("Restart states make Run > Restart from a report step possible: the run continues exactly from that step, e.g. "
             "with a new schedule. A black-oil state is about 30 bytes per cell. Results files keep everything of a run "
             "(also the restart states); ECLIPSE files are written next to the project (<project>_ECLOUT) for ResInsight / OPM tools.");
}

void App::restartWindow() {
    if (!showRestart_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 46, ImGui::GetFontSize() * 30), ImGuiCond_Appearing);
    if (!ImGui::Begin(ICON_PLAY "  Restart from a report step###restartwin", &showRestart_, ImGuiWindowFlags_NoDocking)) { ImGui::End(); return; }
    // base run: a loaded results file or the current results
    std::vector<ressim::ResultsFile::Restart> list;
    std::string baseEngine; int baseSteps = 0; double baseEnd = 0;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& B = restartBase_ ? *restartBase_ : runner_.results();
        for (const auto& R : B.restarts) { ressim::ResultsFile::Restart x = R; x.state.clear(); list.push_back(std::move(x)); }
        baseEngine = B.engine; baseSteps = int(B.snaps.size()) - 1; baseEnd = B.snaps.empty() ? 0.0 : B.snaps.back().time;
    }
    ImGui::TextWrapped("Continue a run from a saved engine state with the current model. Keep the grid, wells and fluids; the "
                       "schedule after the restart step may be changed (new rates, new controls, shut-ins ...).");
    ImGui::Separator();
    ImGui::Text("Base run: %s", restartBase_ ? restartBaseFile_.c_str() : "current results");
    ImGui::SameLine();
    if (ImGui::SmallButton("Load base from .rsres...")) {
        std::string p;
        if (fileDialog(false, "ResSim results (*.rsres)\0*.rsres\0All files (*.*)\0*.*\0", "rsres", p)) {
            try { restartBase_ = std::make_shared<RunResults>(loadRunResults(p)); restartBaseFile_ = fs::u8path(p).filename().u8string(); restartSel_ = -1; }
            catch (const std::exception& e) { notify(std::string("Cannot load: ") + e.what()); }
        }
    }
    if (restartBase_) { ImGui::SameLine(); if (ImGui::SmallButton("Use current results")) { restartBase_.reset(); restartBaseFile_.clear(); restartSel_ = -1; } }
    ImGui::TextDisabled("%s, %d report steps to t = %.6g, %zu restart states", baseEngine.c_str(), std::max(0, baseSteps), baseEnd, list.size());
    if (list.empty()) {
        ImGui::Spacing();
        ImGui::TextColored(col4(Pal::Warn), "No restart states in this run.");
        ui::hint("Restart states are saved by the built-in engine every N report steps and at the end of a run (Settings > "
                 "Results and restart). Run the model once, or load a results file that contains states.");
    } else if (ImGui::BeginTable("##rst", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable,
                                 ImVec2(0, ImGui::GetFontSize() * 12))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Report step"); ImGui::TableSetupColumn("Time"); ImGui::TableSetupColumn("Engine"); ImGui::TableSetupColumn("Date");
        ImGui::TableHeadersRow();
        for (int i = 0; i < int(list.size()); ++i) {
            const auto& R = list[size_t(i)];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            char b[48]; std::snprintf(b, sizeof b, "%d##rs%d", R.report, i);
            if (ImGui::Selectable(b, restartSel_ == i, ImGuiSelectableFlags_SpanAllColumns)) restartSel_ = i;
            ImGui::TableNextColumn(); ImGui::Text("%.6g d", R.time);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(kindName(R.kind));
            ImGui::TableNextColumn();
            if (R.report >= 0 && R.report <= int(case_.schedule.size())) {
                double d = 0; for (int s = 0; s < R.report && s < int(case_.schedule.size()); ++s) d += case_.schedule[size_t(s)].length;
                ImGui::TextDisabled("model day %.6g", d);
            }
        }
        ImGui::EndTable();
    }
    if (restartSel_ >= int(list.size())) restartSel_ = -1;
    const bool canRun = restartSel_ >= 0 && !runner_.running();
    if (canRun) {
        const auto& R = list[size_t(restartSel_)];
        const int left = int(case_.schedule.size()) - R.report;
        if (left <= 0) ImGui::TextColored(col4(Pal::Warn), "The current schedule ends at or before this step - extend it (Schedule) to continue.");
        else ImGui::Text("Continues with %d report step(s) of the current schedule.", left);
    }
    ImGui::Spacing();
    ImGui::BeginDisabled(!canRun);
    if (ImGui::Button(ICON_PLAY "  Restart run", ImVec2(ImGui::GetFontSize() * 10, 0))) {
        std::shared_ptr<RestartRequest> rq = std::make_shared<RestartRequest>();
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& B = restartBase_ ? *restartBase_ : runner_.results();
            const auto& R = B.restarts[size_t(restartSel_)];
            rq->step = R.report; rq->state = R.state; rq->kind = R.kind;
            auto copy = std::make_shared<RunResults>();
            copy->summary = B.summary; copy->snaps = B.snaps; copy->initial = B.initial; copy->components = B.components; copy->trackName = B.trackName;
            copy->fracNames = B.fracNames; copy->fracSegs = B.fracSegs; copy->wellFracLinks = B.wellFracLinks; copy->heatInPlace = B.heatInPlace;
            copy->tInjRef = B.tInjRef; copy->fracSnap = B.fracSnap; copy->wellT = B.wellT; copy->wellMW = B.wellMW; copy->cooledFrac = B.cooledFrac;
            copy->restarts = B.restarts;
            rq->base = copy;
        }
        if (rq->step >= int(case_.schedule.size())) notify("Nothing to run: the schedule ends at the restart step");
        else {
            runner_.requestRestart(rq);
            runSimulation();
            if (runner_.running()) { showRestart_ = false; workspace_ = Workspace::Results; }
            else runner_.requestRestart(nullptr);
        }
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::TextDisabled("(%s engine: the state must match the engine the model runs on)", "built-in");
    ImGui::End();
}
