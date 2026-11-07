// Result Viewer (ECLIPSE / OPM Flow binary + formatted output and ResSim runs: summary overlay, case comparison,
// 3D INIT / restart properties), simulation-engine selection, OPM Flow settings, persistent settings file.
#include "app.h"
#include "i18n.h"
#include "ecl_engine.h"
#include "widgets.h"
#include "implot.h"
#include "plotkit.h"
#include "icons.h"
#include "theme.h"
#include <chrono>
#include <cfloat>
#include "platform.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

using namespace ressim;
std::vector<LineVertex> streamlineVertices(const StreamlineResult& S, int colorBy, Colormap cm, const std::string& focus, bool focusOnly);
namespace fs = std::filesystem;

namespace {
std::string fileNameOf(const std::string& p) { const size_t s = p.find_last_of("\\/"); return s == std::string::npos ? p : p.substr(s + 1); }
std::string lower(std::string s) { for (char& c : s) c = char(std::tolower((unsigned char)c)); return s; }
const char* kCats[] = {"All", "Field", "Wells", "Groups", "Regions", "Blocks", "Connections", "Aquifers", "Other"};
bool catMatch(int cat, char c) {
    switch (cat) {
        case 0: return true;
        case 1: return c == 'F'; case 2: return c == 'W'; case 3: return c == 'G'; case 4: return c == 'R';
        case 5: return c == 'B'; case 6: return c == 'C'; case 7: return c == 'A';
        default: return !std::strchr("FWGRBCA", c);
    }
}
void openFolder(const std::string& dir) {
    plat::openPath(dir);
}
}

// ============================================================================ settings file
void App::loadSettings() {
    settingsPath_ = (fs::u8path(plat::configDir()) / "ressim_settings.ini").u8string();
    std::ifstream f(fs::u8path(settingsPath_));
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t eq = line.find('=');
        if (eq == std::string::npos || line[0] == '#' || line[0] == ';') continue;
        const std::string k = line.substr(0, eq), v = line.substr(eq + 1);
        try {
            if (k == "ui_scale") settings_.uiScale = std::clamp(std::stof(v), 0.6f, 1.6f);
            else if (k == "engine") settings_.engine = std::clamp(std::stoi(v), 0, 2);
            else if (k == "theme") settings_.theme = std::clamp(std::stoi(v), 0, 2);
            else if (k == "accent") settings_.accent = std::clamp(std::stoi(v), 0, 4);
            else if (k == "help_style") { settings_.helpStyle = std::clamp(std::stoi(v), 0, 1); ui::hintMode() = settings_.helpStyle; }
            else if (k == "autosave") settings_.autosave = v != "0";
            else if (k == "autosave_sec") settings_.autosaveSec = std::clamp(std::stoi(v), 15, 3600);
            else if (k == "time_cursor") settings_.timeCursor = v != "0";
            else if (k == "bg3d") { settings_.bg3d = std::clamp(std::stoi(v), 0, 4); set3dBackground(settings_.bg3d); }
            else if (k == "language") { settings_.language = v; i18n::setLanguage(v); }
            else if (k == "ecl_path") settings_.eclPath = v;
            else if (k == "ecl_version") settings_.eclVersion = v;
            else if (k == "ecl_program") settings_.eclProgram = std::clamp(std::stoi(v), 0, 2);
            else if (k == "ecl_args") settings_.eclArgs = v;
            else if (k == "deck_run_mode") settings_.deckRunMode = std::stoi(v) == 1 ? 1 : 0;
            else if (k == "flow_path") settings_.flowPath = v;
            else if (k == "opm_threads") settings_.opmThreads = std::max(0, std::stoi(v));
            else if (k == "opm_max_dt") settings_.opmMaxDt = std::max(0.0f, std::stof(v));
            else if (k == "opm_tuning") settings_.opmTuning = v == "1";
            else if (k == "opm_auto_view") settings_.opmAutoView = v != "0";
            else if (k == "opm_args") settings_.opmArgs = v;
            else if (k == "python_exe") settings_.pythonExe = v;
            else if (k == "opm_keep_wells") settings_.opmKeepWells = v == "1";
            else if (k == "rst_every") settings_.rstEvery = std::clamp(std::stoi(v), 0, 100000);
            else if (k == "rst_last") settings_.rstLast = v != "0";
            else if (k == "auto_save_results") settings_.autoSaveResults = v == "1";
            else if (k == "ecl_output") settings_.eclOutput = v == "1";
            else if (k == "mpi_enabled") settings_.mpiEnabled = v == "1";
            else if (k == "mpi_exec") settings_.mpiExec = v;
            else if (k == "mpi_procs") settings_.mpiProcs = std::clamp(std::stoi(v), 0, 4096);
            else if (k == "mpi_args") settings_.mpiArgs = v;
            else if (k == "mpi_jobdir") settings_.mpiJobDir = v;
            else if (k == "gpu_enabled") settings_.gpuEnabled = v == "1";
            else if (k == "gpu_prec") settings_.gpuPrec = v == "0" ? 0 : v == "2" ? 2 : 1;
            else if (k == "mpi_dd") settings_.mpiDD = v == "1";
            else if (k == "mpi_dd_procs") settings_.mpiDDProcs = std::clamp(std::stoi(v), 1, 4096);
            else if (k == "mpi_da") settings_.mpiDistAsm = v == "1";
            else if (k == "range_all_steps") settings_.rangeAllSteps = v == "1";
            else if (k == "show_start") settings_.showStartPage = v != "0";
            else if (k == "run_history") settings_.runHistory = v != "0";
            else if (k == "run_history_max") settings_.runHistoryMax = std::clamp(std::stoi(v), 1, 50);
            else if (k.rfind("recent", 0) == 0 && !v.empty() && settings_.recent.size() < 12) settings_.recent.push_back(v);
        } catch (...) {}
    }
}

void App::saveSettings() {
    if (settingsPath_.empty()) return;
    std::ofstream f(fs::u8path(settingsPath_), std::ios::binary);
    if (!f) return;
    f << "# ResSim Studio settings\n";
    f << "ui_scale=" << settings_.uiScale << "\n";
    f << "engine=" << settings_.engine << "\n";
    f << "theme=" << settings_.theme << "\naccent=" << settings_.accent << "\n";
    f << "help_style=" << settings_.helpStyle << "\n";
    f << "autosave=" << (settings_.autosave ? 1 : 0) << "\nautosave_sec=" << settings_.autosaveSec << "\n";
    f << "time_cursor=" << (settings_.timeCursor ? 1 : 0) << "\n";
    f << "bg3d=" << settings_.bg3d << "\n";
    f << "language=" << settings_.language << "\n";
    f << "flow_path=" << settings_.flowPath << "\n";
    f << "opm_threads=" << settings_.opmThreads << "\n";
    f << "opm_max_dt=" << settings_.opmMaxDt << "\n";
    f << "opm_tuning=" << (settings_.opmTuning ? 1 : 0) << "\n";
    f << "opm_auto_view=" << (settings_.opmAutoView ? 1 : 0) << "\n";
    f << "opm_args=" << settings_.opmArgs << "\n";
    f << "python_exe=" << settings_.pythonExe << "\n";
    f << "opm_keep_wells=" << (settings_.opmKeepWells ? 1 : 0) << "\n";
    f << "ecl_path=" << settings_.eclPath << "\n";
    f << "ecl_version=" << settings_.eclVersion << "\n";
    f << "ecl_program=" << settings_.eclProgram << "\n";
    f << "ecl_args=" << settings_.eclArgs << "\n";
    f << "deck_run_mode=" << settings_.deckRunMode << "\n";
    f << "rst_every=" << settings_.rstEvery << "\n";
    f << "rst_last=" << (settings_.rstLast ? 1 : 0) << "\n";
    f << "auto_save_results=" << (settings_.autoSaveResults ? 1 : 0) << "\n";
    f << "ecl_output=" << (settings_.eclOutput ? 1 : 0) << "\n";
    f << "mpi_enabled=" << (settings_.mpiEnabled ? 1 : 0) << "\n";
    f << "mpi_exec=" << settings_.mpiExec << "\n";
    f << "mpi_procs=" << settings_.mpiProcs << "\n";
    f << "mpi_args=" << settings_.mpiArgs << "\n";
    f << "mpi_jobdir=" << settings_.mpiJobDir << "\n";
    f << "gpu_enabled=" << (settings_.gpuEnabled ? 1 : 0) << "\n";
    f << "gpu_prec=" << settings_.gpuPrec << "\n";
    f << "mpi_dd=" << (settings_.mpiDD ? 1 : 0) << "\n";
    f << "mpi_dd_procs=" << settings_.mpiDDProcs << "\n";
    f << "mpi_da=" << (settings_.mpiDistAsm ? 1 : 0) << "\n";
    f << "range_all_steps=" << (settings_.rangeAllSteps ? 1 : 0) << "\n";
    f << "show_start=" << (settings_.showStartPage ? 1 : 0) << "\n";
    f << "run_history=" << (settings_.runHistory ? 1 : 0) << "\n";
    f << "run_history_max=" << settings_.runHistoryMax << "\n";
    for (size_t i = 0; i < settings_.recent.size(); ++i) f << "recent" << i << "=" << settings_.recent[i] << "\n";
}

void App::detectEclipse(bool force) {
    if (eclDetectDone_ && !force) return;
    eclDetectDone_ = true;
    eclDetected_ = detectEclipseExecutable(settings_.eclPath);
    eclVersions_ = eclDetected_.empty() ? std::vector<std::string>() : eclipseVersions(eclDetected_, settings_.eclProgram == 1 ? "e300" : settings_.eclProgram == 2 ? "frontsim" : "eclipse");
}

void App::detectFlow(bool force) {
    if (flowDetectDone_ && !force) return;
    flowDetectDone_ = true;
    flowDetected_ = detectFlowExecutable(settings_.flowPath);
    flowVersion_ = flowDetected_.empty() ? std::string() : flowVersion(flowDetected_);
}

// ============================================================================ OPM Flow settings window
void App::opmSettingsWindow() {
    if (!showOpmSettings_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 44, 0), ImGuiCond_Appearing);
    if (!ImGui::Begin("Settings", &showOpmSettings_, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_AlwaysAutoResize)) { ImGui::End(); return; }
    detectFlow(false);
    if (settings_.engine == 2) detectEclipse(false);
    bool ch = false;
    ui::header("Interface");
    {
        int pct = int(std::lround(settings_.uiScale * 100));
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        if (ImGui::SliderInt("Interface size", &pct, 60, 160, "%d %%", ImGuiSliderFlags_AlwaysClamp)) {}
        if (ImGui::IsItemDeactivatedAfterEdit()) { settings_.uiScale = pct / 100.0f; setUiScale(settings_.uiScale); ch = true; }
        ImGui::SameLine();
        for (int v : {80, 90, 100}) {
            char b[16]; std::snprintf(b, sizeof(b), "%d%%##sz", v);
            if (ImGui::SmallButton(b)) { settings_.uiScale = v / 100.0f; setUiScale(settings_.uiScale); ch = true; }
            ImGui::SameLine();
        }
        ImGui::NewLine();
        ui::hint("Smaller sizes fit more panels and wider layouts on the screen. Shortcuts: Ctrl + / Ctrl - / Ctrl 0.");
        {   // C39 appearance
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9);
            if (ImGui::Combo("Theme", &settings_.theme, themeNames(), 3)) { applyColors(settings_.theme, settings_.accent); ch = true; }
            ImGui::SameLine(0, ImGui::GetFontSize() * 1.5f);
            ImGui::TextUnformatted("Accent");
            for (int a = 0; a < 5; ++a) {
                ImGui::SameLine();
                const ImVec4 col = accentColor(a);
                ImGui::PushID(a);
                if (ImGui::ColorButton(accentNames()[a], col, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoDragDrop,
                                       ImVec2(ImGui::GetFrameHeight(), ImGui::GetFrameHeight()))) { settings_.accent = a; applyColors(settings_.theme, a); ch = true; }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", accentNames()[a]);
                if (a == settings_.accent) ImGui::GetWindowDrawList()->AddRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), ImGui::GetColorU32(ImGuiCol_Text), 4.0f, 0, 2.0f);
                ImGui::PopID();
            }
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9);
            if (ImGui::Combo("Help text", &settings_.helpStyle, "Compact (i) icons\0Inline text\0")) { ui::hintMode() = settings_.helpStyle; ch = true; }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Compact: explanations sit behind an (i) icon - hover to read, click to pin. Inline: the explanations are printed under the controls."));
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9);
            if (ImGui::Combo("3-D background", &settings_.bg3d, background3dNames(), 5)) { set3dBackground(settings_.bg3d); ch = true; }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Model and Results 3-D views (the Geomodel 3-D scene has its own colour in its View panel)"));
            ch |= ImGui::Checkbox(TRL("Report-time cursor on all time plots (drag it to change the report step)"), &settings_.timeCursor);
            ch |= ImGui::Checkbox(TRL("Autosave a recovery copy of unsaved changes"), &settings_.autosave);
            if (settings_.autosave) {
                ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
                if (ImGui::InputInt("s##autosave", &settings_.autosaveSec, 15, 60)) { settings_.autosaveSec = std::clamp(settings_.autosaveSec, 15, 3600); ch = true; }
            }
            ui::hint("Undo / redo of every model edit: Ctrl+Z / Ctrl+Y (Edit menu). The recovery copy is offered at the next start "
                     "if the program did not close normally.");
        }
        ch |= ImGui::Checkbox(TRL("Start page at start-up"), &settings_.showStartPage);
        ch |= ImGui::Checkbox(TRL("Keep every finished run for comparison (run history)"), &settings_.runHistory);
        if (settings_.runHistory) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            if (ImGui::InputInt("runs kept", &settings_.runHistoryMax)) { settings_.runHistoryMax = std::clamp(settings_.runHistoryMax, 1, 50); ch = true; }
        }
        ui::hint("Run history: Results > Cases, vectors & compare lists the last runs with what changed in the inputs (hover a run). "
                 "Command palette: Ctrl+Shift+P or Ctrl+K.");
    }
    restartSettingsUi(ch);
    hpcSettingsUi(ch);
    ui::header("Engine used by Run (F5)");
    ch |= ImGui::RadioButton(TRL("ResSim (built-in, fully implicit black oil)"), &settings_.engine, 0);
    ch |= ImGui::RadioButton(TRL("OPM Flow (external, open-source flow.exe)"), &settings_.engine, 1);
    if (ImGui::RadioButton(TRL("ECLIPSE (external, licensed Schlumberger installation)"), &settings_.engine, 2)) { ch = true; detectEclipse(false); }
    ui::hint("With an external engine, ResSim exports the case as an ECLIPSE deck, runs the simulator with live progress and "
             "reads the EGRID / INIT / UNRST / SMSPEC results back into the plots, the 3D viewer and the Result Viewer.");
    if (settings_.engine != 0) {
        ui::header("Imported decks");
        ch |= ImGui::RadioButton(TRL("Run the original deck file while it is not edited in ResSim (recommended)"), &settings_.deckRunMode, 0);
        ch |= ImGui::RadioButton(TRL("Always run the deck generated from the model"), &settings_.deckRunMode, 1);
        ui::hint("The original deck (with all INCLUDE files, groups, VFP, MSW, ACTIONX, UDQ ...) runs exactly as written, so every "
                 "feature of OPM Flow / ECLIPSE is available. After an edit in ResSim the model is exported as a new deck.");
        if (!srcDeck_.empty())
            ImGui::TextDisabled("Imported: %s  (%s)", srcDeck_.c_str(), caseVersion_ == srcDeckVersion_ ? "unchanged - original deck is run" : "edited - generated deck is run");
    }

    if (settings_.engine == 2) {
        ui::header("ECLIPSE executable");
        ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 11);
        if (ui::inputString("##eclpath", &settings_.eclPath)) { ch = true; eclDetectDone_ = false; }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", TR("eclrun.exe (recommended: picks the version and the licence settings), eclipse.exe / e300.exe, or the \n"
                              "installation folder (e.g. C:\\ecl). Empty = automatic detection."));
        ImGui::SameLine();
        if (ImGui::Button(TRL("Browse...##ecl"))) {
            std::string p = settings_.eclPath;
            if (fileDialog(false, "ECLIPSE (eclrun.exe, eclipse.exe, e300.exe)\0eclrun*.exe;eclipse*.exe;e300*.exe\0Programs (*.exe)\0*.exe\0All files (*.*)\0*.*\0",
                           "exe", p)) {
                settings_.eclPath = p; ch = true; detectEclipse(true);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button(TRL("Detect##ecl"))) detectEclipse(true);
        if (eclDetected_.empty()) {
            ImGui::TextColored(col4(Pal::Warn), "ECLIPSE not found.");
            ui::hint("ResSim does not include ECLIPSE: a licensed installation is needed. Searched: the path above, %ECLIPSE_EXE%, "
                     "%ECLRUN%, %ECLPATH%, PATH, C:\\ecl ... F:\\ecl (macros\\eclrun.exe or <version>\\bin\\pc_x86_64\\eclipse.exe), "
                     "Program Files\\Schlumberger and Program Files\\SLB. Press Browse to choose eclrun.exe or eclipse.exe.");
        } else {
            const bool launcher = isEclrunLauncher(eclDetected_);
            ImGui::TextColored(col4(Pal::Ok), "\xE2\x97\x8F %s", launcher ? "eclrun launcher found" : "ECLIPSE simulator found");
            ImGui::SameLine(); ImGui::TextDisabled("%s", eclDetected_.c_str());
            if (!launcher)
                ui::hint("The simulator is started directly and receives the case name on its input. Its licence must be reachable "
                         "through the environment (SLBSLS_LICENSE_FILE / LM_LICENSE_FILE). eclrun.exe sets this up itself.");
        }
        ui::header("ECLIPSE run options");
        if (ui::beginProps("eclopt", 0.42f)) {
            ui::label("Simulator", "E100: black oil (the ResSim deck format). E300: compositional / advanced black oil - the deck must suit E300.\n"
                      "FrontSim: SLB's streamline simulator (IMPES pressure + 1D streamline transport, same deck format). Its results "
                      "open in the Results workspace; streamlines are traced from the inter-block flows (add RPTRST FLOWS) or from the "
                      "pressure field.");
            static const char* progs[] = {"E100 (eclipse)", "E300 (e300)", "FrontSim (frontsim, streamlines)"};
            if (ImGui::Combo("##eclprog", &settings_.eclProgram, progs, 3)) { ch = true; detectEclipse(true); }
            ui::label("Version", "eclrun -v VERSION. Latest = the newest installed version.");
            const bool launcher = isEclrunLauncher(eclDetected_);
            ImGui::BeginDisabled(!launcher);
            if (ImGui::BeginCombo("##eclver", settings_.eclVersion.empty() ? "latest installed" : settings_.eclVersion.c_str())) {
                if (ImGui::Selectable(TRL("latest installed"), settings_.eclVersion.empty())) { settings_.eclVersion.clear(); ch = true; }
                for (auto& v : eclVersions_)
                    if (ImGui::Selectable(v.c_str(), settings_.eclVersion == v)) { settings_.eclVersion = v; ch = true; }
                ImGui::EndCombo();
            }
            ImGui::EndDisabled();
            if (!launcher && !eclVersions_.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", eclVersions_.front().c_str()); }
            ch |= ui::propBool("Open results in Result Viewer", &settings_.opmAutoView);
            ui::label("Extra arguments", launcher ? "eclrun options placed before the program name, e.g. --np 4 (parallel; the deck needs PARALLEL)"
                                                  : "Appended to the simulator command line");
            ch |= ui::inputString("##eclargs", &settings_.eclArgs);
            ui::endProps();
        }
        const std::string edir = runDirFor(path_, case_.title, "ECL");
        ImGui::TextDisabled("Run folder: %s", edir.c_str());
        if (ImGui::SmallButton(TRL("Open run folder##ecl"))) { std::error_code ec; fs::create_directories(fs::u8path(edir), ec); openFolder(edir); }
        if (ch) saveSettings();
        ImGui::End();
        return;
    }

    ui::header("OPM Flow executable");
    ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 11);
    if (ui::inputString("##flowpath", &settings_.flowPath)) { ch = true; flowDetectDone_ = false; }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("flow.exe / flow_blackoil.exe, or its folder. Empty = automatic detection."));
    ImGui::SameLine();
    if (ImGui::Button(TRL("Browse..."))) {
        std::string p = settings_.flowPath;
        if (fileDialog(false, "OPM Flow (flow*.exe)\0flow*.exe\0Programs (*.exe)\0*.exe\0All files (*.*)\0*.*\0", "exe", p)) {
            settings_.flowPath = p; ch = true; detectFlow(true);
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(TRL("Detect"))) detectFlow(true);
    if (flowDetected_.empty()) {
        ImGui::TextColored(col4(Pal::Warn), "OPM Flow not found.");
        ui::hint("Install OPM Flow for Windows (installer: Program Files\\OPM Flow\\bin\\flow.exe, or unzip the portable package) "
                 "from github.com/GitPaean/opm_flow_windows, then press Detect - or Browse to flow.exe. Searched: the path above, "
                 "%FLOW_EXE%, PATH, Program Files, %LOCALAPPDATA%\\Programs, your profile / Downloads / Desktop and C:\\.");
    } else {
        ImGui::TextColored(col4(Pal::Ok), "\xE2\x97\x8F %s", flowVersion_.empty() ? "found (version unknown)" : flowVersion_.c_str());
        ImGui::SameLine(); ImGui::TextDisabled("%s", flowDetected_.c_str());
        const std::string pre = flowPrerequisiteProblem(flowDetected_);
        if (!pre.empty()) { ImGui::PushTextWrapPos(ImGui::GetFontSize() * 42); ImGui::TextColored(col4(Pal::Warn), "%s", pre.c_str()); ImGui::PopTextWrapPos(); }
        if (flowVersion_.empty()) ui::hint("flow --version did not answer: the runtime (MS-MPI / VC++ redistributable) may be missing.");
    }

    ui::header("OPM Flow run options");
    if (ui::beginProps("opmopt", 0.42f)) {
        ui::label("Threads per process", "--threads-per-process. 0 = automatic (all hardware threads, at most 8).");
        ch |= ImGui::SliderInt("##thr", &settings_.opmThreads, 0, 32, settings_.opmThreads == 0 ? "auto" : "%d");
        ui::label("Max time step (days)", "--solver-max-time-step-in-days. 0 = OPM default (up to 365 days). Smaller steps reduce "
                                          "numerical dispersion and make a comparison with ResSim fairer.");
        ch |= ImGui::DragFloat("##mdt", &settings_.opmMaxDt, 0.5f, 0.0f, 365.0f, settings_.opmMaxDt == 0 ? "OPM default" : "%.3g d");
        ch |= ui::propBool("Use TUNING", &settings_.opmTuning, "--enable-tuning=true: apply the case's TUNING (time-step) settings.");
        ch |= ui::propBool("Keep problematic wells open", &settings_.opmKeepWells,
                           "OPM Flow shuts a well whose equations fail to converge (e.g. an injector that hits its BHP limit at "
                           "start-up, SPE10 model 1). The run then finishes but without that well. On: --shut-unsolvable-wells=false "
                           "and no minimum time step before shutting - OPM cuts the time step instead, as ResSim does.");
        ch |= ui::propBool("Open results in Result Viewer", &settings_.opmAutoView);
        ui::label("Extra arguments", "Appended to the command line, e.g. --linear-solver=cprw");
        ch |= ui::inputString("##args", &settings_.opmArgs);
        ui::endProps();
    }
    const std::string dir = runDirFor(path_, case_.title);
    ImGui::TextDisabled("Run folder: %s", dir.c_str());
    if (ImGui::SmallButton(TRL("Open run folder"))) { std::error_code ec; fs::create_directories(fs::u8path(dir), ec); openFolder(dir); }
    const auto notes = opmCompatibilityNotes(case_);
    if (!notes.empty()) {
        ui::header("Notes for the current case");
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 42);
        for (auto& n : notes) ImGui::TextColored(col4(Pal::Warn), "\xE2\x80\xA2 %s", n.c_str());
        ImGui::PopTextWrapPos();
        if (!case_.grid.cornerPoint() && !case_.grid.newtran && ImGui::Button(TRL("Use geometric transmissibility (NEWTRAN) for this case"))) {
            case_.grid.newtran = true; markDirty(true);
        }
    }
    if (ch) saveSettings();
    ImGui::End();
}

// ============================================================================ Result Viewer: cases
void App::rvOpen(const std::string& path, bool quiet) {
    auto rc = std::make_shared<ecl::ResultCase>();
    std::string err;
    if (!ecl::openResultCase(path, *rc, &err)) {
        if (!quiet) notify("Cannot open results: " + err);
        runner_.addLog("ERROR: Result Viewer: " + err);
        return;
    }
    RVCase c;
    c.rc = rc;
    c.path = rc->files.base;
    c.name = rc->name.empty() ? fileNameOf(path) : rc->name;
    // re-opening the same case replaces it (e.g. after a new OPM run)
    for (auto& o : rv_.cases) {
        if (o.internal || o.path != c.path) continue;
        c.color = o.color; c.show = o.show;
        o = std::move(c);
        rv_.diffs.clear(); rv_.key3d.clear(); rv_.grid3d = nullptr; rv_.cmpRef = rv_.cmpCmp = -1;
        if (!quiet) notify("Reloaded " + o.name);
        return;
    }
    c.color = seriesColor(rv_.colorIdx++);
    rv_.cases.push_back(std::move(c));
    const RVCase& added = rv_.cases.back();
    for (auto& w : added.rc->warnings) runner_.addLog("WARNING: " + added.name + ": " + w);
    if (rv_.selected.empty()) {
        for (const char* k : {"FOPR", "FPR", "FWCT", "FGOR"})
            if (added.sum().find(k) >= 0) rv_.selected.push_back(k);
    }
    if (rv_.cases.size() == 2) { rv_.ref = 0; rv_.cmp = 1; }
    runner_.addLog("=== Result Viewer: opened " + added.name + " (" + std::to_string(added.sum().vars.size()) + " vectors, " +
                   std::to_string(added.rc->hasRestart ? added.rc->restart.steps().size() : 0) + " restart steps)");
    if (!quiet) notify("Opened " + added.name);
    showResultViewer_ = true;
}

void App::rvAddRun() {
    // Copy what we need while holding the runner's mutex, then release it BEFORE notify()/rvOpen():
    // both log through runner_.addLog(), which locks the same (non-recursive) mutex - calling them
    // under the lock deadlocked the UI thread.
    RVCase c;
    std::string title, opmDeck;
    bool ok = false, opm = false;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        ok = R.valid && R.summary.size() >= 2;
        opm = ok && R.engine != "ResSim" && !R.opmDeck.empty();
        if (opm) opmDeck = R.opmDeck;
        else if (ok) {
            c.internal = true;
            c.summary = ecl::summaryFromRun(R.summary, R.wellNames, R.input.units);
            c.grid = R.grid;
            c.snaps = R.snaps;
            c.input = std::make_shared<ressim::Case>(R.input);
            c.sim = R.sim;
            c.egs.fracNames = R.fracNames; c.egs.fracSegs = R.fracSegs; c.egs.fracSnap = R.fracSnap;
            title = R.input.title.substr(0, 40);
        }
    }
    if (!ok) { notify("No simulation results yet - run the case first"); return; }
    if (opm) { rvDeferred_ = [this, p = opmDeck] { rvOpen(p); }; return; }
    char b[96];
    std::snprintf(b, sizeof(b), "ResSim run %d", int(std::count_if(rv_.cases.begin(), rv_.cases.end(), [](const RVCase& x) { return x.internal; })) + 1);
    c.name = std::string(b) + ": " + title;
    c.color = seriesColor(rv_.colorIdx++);
    rv_.cases.push_back(std::move(c));
    if (rv_.selected.empty()) rv_.selected = {"FOPR", "FPR", "FWCT", "FGOR"};
    if (rv_.cases.size() == 2) { rv_.ref = 0; rv_.cmp = 1; }
    notify("Added the current ResSim run to the Result Viewer");
    showResultViewer_ = true;
}

// ============================================================================ Result Viewer window
void App::panelResultViewer() {
    if (!showResultViewer_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 80, ImGui::GetFontSize() * 46), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Result Viewer", &showResultViewer_)) { ImGui::End(); return; }
    rvBody();
    ImGui::End();
}

void App::rvBody() {
    const float fs = ImGui::GetFontSize();
    // ---------------- left: cases + vector tree
    ImGui::BeginChild("rvleft", ImVec2(fs * 19, 0), ImGuiChildFlags_ResizeX | ImGuiChildFlags_Borders);
    if (ImGui::Button(TRL("Open results..."))) {
        std::string p;
        if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.FEGRID;*.UNRST;*.FUNRST;*.INIT;*.FINIT;*.UNSMRY;*.FUNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p))
            rvOpen(p);
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Any file of a case (CASE.DATA, .SMSPEC, .EGRID, .UNRST, formatted .F* files): the\nother files of the same case are found automatically."));
    ImGui::SameLine();
    if (ImGui::Button(TRL("+ Current run"))) rvAddRun();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Add the results of the last Run (ResSim or OPM Flow) as a case"));
    ui::header("Cases");
    int remove = -1;
    for (int i = 0; i < int(rv_.cases.size()); ++i) {
        RVCase& c = rv_.cases[size_t(i)];
        ImGui::PushID(i);
        ImGui::Checkbox(TRL("##show"), &c.show);
        ImGui::SameLine();
        ImVec4 col = col4(c.color);
        ImGui::ColorEdit4("##col", &col.x, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_NoAlpha);
        c.color = ImGui::ColorConvertFloat4ToU32(col);
        ImGui::SameLine();
        ImGui::TextUnformatted(c.name.c_str());
        if (ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            ImGui::Text("%s", c.internal ? "ResSim run (in memory)" : c.path.c_str());
            ImGui::Text("%zu vectors, %zu time steps", c.sum().vars.size(), c.sum().time.size());
            if (!c.changes.empty()) {
                ImGui::Separator();
                ImGui::TextDisabled("%s", TR("Changes against the previous run:"));
                if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
                ImGui::TextUnformatted(c.changes.c_str());
                if (g_fonts.mono) ImGui::PopFont();
            }
            ImGui::EndTooltip();
        }
        if (ImGui::BeginPopupContextItem("ctx")) {
            if (!c.internal && ImGui::MenuItem(TRL("Reload"))) {
                const auto& F = c.rc->files;
                std::string p = !F.smspec.empty() ? F.smspec : !F.egrid.empty() ? F.egrid : !F.restart.empty() ? F.restart[0] : F.init;
                rvDeferred_ = [this, p] { rvOpen(p); };
            }
            if (!c.internal && ImGui::MenuItem(TRL("Open folder"))) openFolder(fs::u8path(c.path).parent_path().u8string());
            if (ImGui::MenuItem(TRL("Remove"))) remove = i;
            ImGui::EndPopup();
        }
        ImGui::PopID();
    }
    if (remove >= 0) { rv_.cases.erase(rv_.cases.begin() + remove); rv_.diffs.clear(); rv_.key3d.clear(); rv_.grid3d = nullptr; rv_.cmpRef = rv_.cmpCmp = -1; }
    if (rv_.cases.empty()) ui::hint("Open ECLIPSE / OPM Flow output (binary or formatted), or add the current run. Run the same case "
                                    "with both engines to compare them (Run > Engine).");

    ui::header("Vectors");
    ImGui::SetNextItemWidth(fs * 7);
    ImGui::Combo("##cat", &rv_.category, kCats, IM_ARRAYSIZE(kCats));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##flt", "filter (e.g. WBHP, PROD)", rv_.filter, sizeof(rv_.filter));
    // union of vectors of the shown cases, grouped by keyword
    std::map<std::string, std::vector<std::string>> groups;   // keyword -> keys
    std::map<std::string, std::string> units;
    {
        const std::string f = lower(rv_.filter);
        std::set<std::string> seen;
        for (auto& c : rv_.cases) {
            if (!c.show) continue;
            for (auto& v : c.sum().vars) {
                if (v.category == 'T' || !catMatch(rv_.category, v.category)) continue;
                if (!f.empty() && lower(v.key).find(f) == std::string::npos) continue;
                if (!seen.insert(v.key).second) continue;
                groups[v.keyword].push_back(v.key);
                units[v.key] = v.unit;
            }
        }
    }
    ImGui::BeginChild("vectree", ImVec2(0, 0));
    auto toggle = [&](const std::string& k) {
        auto it = std::find(rv_.selected.begin(), rv_.selected.end(), k);
        if (it == rv_.selected.end()) rv_.selected.push_back(k); else rv_.selected.erase(it);
    };
    auto isSel = [&](const std::string& k) { return std::find(rv_.selected.begin(), rv_.selected.end(), k) != rv_.selected.end(); };
    for (auto& [kw, keys] : groups) {
        if (keys.size() == 1 && keys[0] == kw) {
            if (ImGui::Selectable(kw.c_str(), isSel(kw))) toggle(kw);
            if (!units[kw].empty()) { ImGui::SameLine(fs * 7); ImGui::TextDisabled("%s", units[kw].c_str()); }
            continue;
        }
        const bool open = ImGui::TreeNodeEx(kw.c_str(), ImGuiTreeNodeFlags_SpanAvailWidth | (rv_.filter[0] ? ImGuiTreeNodeFlags_DefaultOpen : 0), "%s (%zu)", kw.c_str(), keys.size());
        if (open) {
            for (auto& k : keys) {
                const size_t c = k.find(':');
                if (ImGui::Selectable(c == std::string::npos ? k.c_str() : k.c_str() + c + 1, isSel(k))) toggle(k);
            }
            ImGui::TreePop();
        }
    }
    ImGui::EndChild();
    ImGui::EndChild();
    ImGui::SameLine();

    // ---------------- right: tabs
    ImGui::BeginChild("rvright", ImVec2(0, 0));
    if (ImGui::BeginTabBar("rvtabs")) {
        auto tab = [&](const char* name, int id) {
            ImGuiTabItemFlags fl = rv_.wantTab == id ? ImGuiTabItemFlags_SetSelected : 0;
            const bool open = ui::tabItem(name, nullptr, fl);
            if (open) rv_.tab = id;
            return open;
        };
        if (tab("Plots", 0)) { rvPlots(); ImGui::EndTabItem(); }
        if (tab("Compare", 1)) { rvCompare(); ImGui::EndTabItem(); }
        if (tab("3D", 2)) { rv3D(); ImGui::EndTabItem(); }
        if (tab("Info", 3)) { rvInfo(); ImGui::EndTabItem(); }
        rv_.wantTab = -1;
        ImGui::EndTabBar();
    }
    ImGui::EndChild();
}

// ---------------------------------------------------------------------------- plots
void App::rvPlots() {
    const float fs = ImGui::GetFontSize();
    ImGui::Checkbox(TRL("Years"), &rv_.years);
    ImGui::SameLine(); ImGui::Checkbox(TRL("One plot per vector"), &rv_.perVector);
    ImGui::SameLine();
    static const char* colorModes[] = {"Colours: auto", "Colours: by case", "Colours: by vector", "Colours: every curve"};
    ImGui::SetNextItemWidth(fs * 9.5f);
    ImGui::Combo("##colmode", &rv_.colorMode, colorModes, 4);
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
        ImGui::SetTooltip("%s", TR("auto: one plot per vector -> one colour per case (one case: per vector); all in one plot -> every curve\n"
                          "gets its own colour. Palette: Results toolbar."));
    if (rv_.perVector) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 6.5f);
        static const char* colsItems[] = {"Grid: auto", "1 column", "2 columns", "3 columns", "4 columns"};
        ImGui::Combo("##cols", &rv_.gridCols, colsItems, 5);
    }
    ImGui::SameLine();
    if (ImGui::SmallButton(TRL("Clear selection"))) rv_.selected.clear();
    ImGui::SameLine();
    ImGui::TextDisabled("%zu vector(s)", rv_.selected.size());
    ImGui::SameLine();
    exportButton("rv.plots");
    if (rv_.selected.empty() || rv_.cases.empty()) { ui::hint("Select vectors in the tree on the left (FOPR, WBHP:..., BPR:...)."); return; }
    const double xs = rv_.years ? 1.0 / 365.25 : 1.0;
    const char* xl = rv_.years ? "Time (years)" : "Time (days)";
    int shownCases = 0;
    for (auto& c : rv_.cases) if (c.show) ++shownCases;
    // colour of a curve: 0 by case, 1 by vector, 2 every curve
    int mode = rv_.colorMode == 1 ? 0 : rv_.colorMode == 2 ? 1 : rv_.colorMode == 3 ? 2 : -1;
    if (mode < 0) mode = rv_.perVector ? (shownCases > 1 ? 0 : 1) : 2;
    int curve = 0;
    auto drawKey = [&](const std::string& key, int keyIndex, bool labelWithKey) {
        int ci = 0;
        for (auto& c : rv_.cases) {
            if (!c.show) continue;
            const int caseIdx = ci++;
            const ecl::Summary& S = c.sum();
            const int v = S.find(key);
            if (v < 0) continue;
            const auto& y = S.values[size_t(v)];
            const size_t n = std::min(y.size(), S.time.size());
            if (n == 0) continue;
            std::vector<double> X(n), Y(n);
            for (size_t k = 0; k < n; ++k) { X[k] = S.time[k] * xs; Y[k] = y[k]; }
            const ImU32 col = mode == 0 ? c.color : mode == 1 ? seriesColor(keyIndex) : seriesColor(curve);
            ++curve;
            (void)caseIdx;
            ImPlotSpec spec;
            spec.LineColor = col4(col);
            spec.LineWeight = (c.internal ? 2.0f : 1.6f) * g_plotLineScale;
            const bool markers = n < 40;
            if (markers) { spec.Marker = ImPlotMarker_Circle; spec.MarkerSize = 2.5f; spec.MarkerFillColor = spec.LineColor; spec.MarkerLineColor = spec.LineColor; }
            std::string lab;
            if (labelWithKey) lab = shownCases > 1 ? c.name + ": " + key : key;
            else lab = shownCases > 1 || mode == 0 ? c.name : key;
            ImPlot::PlotLine(lab.c_str(), X.data(), Y.data(), int(n), spec);
            recSeries(lab.c_str(), X.data(), Y.data(), int(n), col, c.internal ? 2.0f : 1.6f, false, markers);
        }
    };
    auto unitOfKey = [&](const std::string& key) -> std::string {
        // shown cases first (a hidden case with other units - e.g. a METRIC ECLIPSE run next to a FIELD case - must
        // not label the axis), then any case that has the vector
        for (int pass = 0; pass < 2; ++pass)
            for (auto& c : rv_.cases) {
                if (pass == 0 && !c.show) continue;
                const int v = c.sum().find(key);
                if (v >= 0 && !c.sum().vars[size_t(v)].unit.empty()) return c.sum().vars[size_t(v)].unit;
            }
        return {};
    };
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
    const ImVec2 av = ImGui::GetContentRegionAvail();
    if (rv_.perVector) {
        const int n = int(rv_.selected.size());
        const PlotGrid G = fitGrid(n, av, fs * 20, fs * 10, 1.7f, rv_.gridCols);
        if (G.scroll) ImGui::BeginChild("##rvscroll", av);
        fig::recorder().beginFigure("rv.plots", G.rows, G.cols);
        if (ImPlot::BeginSubplots("##rvsub", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
            int ki = 0;
            for (auto& key : rv_.selected) {
                if (ImPlot::BeginPlot(key.c_str())) {
                    const std::string u = unitOfKey(key);
                    ImPlot::SetupAxes(xl, u.c_str(), af, af);
                    setupSI(ImAxis_Y1);
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast);
                    recPanel(key.c_str(), xl, u.empty() ? key.c_str() : (key + " (" + u + ")").c_str());
                    drawKey(key, ki, false);
                    recEnd();
                    ImPlot::EndPlot();
                }
                ++ki;
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        if (G.scroll) ImGui::EndChild();
    } else {
        fig::recorder().beginFigure("rv.plots", 1, 1);
        if (ImPlot::BeginPlot("##rvall", ImVec2(-1, av.y))) {
            // up to two unit groups: first on the left axis, second on the right
            std::vector<std::string> unitList;
            for (auto& k : rv_.selected) { std::string u = unitOfKey(k); if (std::find(unitList.begin(), unitList.end(), u) == unitList.end()) unitList.push_back(u); }
            ImPlot::SetupAxes(xl, unitList.empty() ? "" : unitList[0].c_str(), af, af);
            setupSI(ImAxis_Y1);
            if (unitList.size() > 1) { ImPlot::SetupAxis(ImAxis_Y2, unitList[1].c_str(), af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
            recPanel("", xl, unitList.empty() ? "" : unitList[0].c_str(), unitList.size() > 1 ? unitList[1].c_str() : "");
            int ki = 0;
            for (auto& k : rv_.selected) {
                const std::string u = unitOfKey(k);
                useY2(unitList.size() > 1 && u == unitList[1]);
                drawKey(k, ki++, true);
            }
            recEnd();
            ImPlot::EndPlot();
        }
        fig::recorder().endFigure();
    }
}

// ---------------------------------------------------------------------------- compare
void App::rvCompare() {
    const float fs = ImGui::GetFontSize();
    if (rv_.cases.size() < 2) { ui::hint("Open at least two cases (e.g. the same deck run with ResSim and with OPM Flow) to compare them."); return; }
    const int nc = int(rv_.cases.size());
    rv_.ref = std::min(rv_.ref, nc - 1); rv_.cmp = std::min(rv_.cmp, nc - 1);
    auto combo = [&](const char* id, int* v) {
        ImGui::SetNextItemWidth(fs * 18);
        if (ImGui::BeginCombo(id, rv_.cases[size_t(*v)].name.c_str())) {
            for (int i = 0; i < nc; ++i) if (ImGui::Selectable(rv_.cases[size_t(i)].name.c_str(), *v == i)) *v = i;
            ImGui::EndCombo();
        }
    };
    ImGui::AlignTextToFramePadding(); ImGui::TextUnformatted("Reference"); ImGui::SameLine(); combo("##ref", &rv_.ref);
    ImGui::SameLine(); ImGui::TextUnformatted("vs"); ImGui::SameLine(); combo("##cmp", &rv_.cmp);
    ImGui::SameLine();
    bool redo = rv_.cmpRef != rv_.ref || rv_.cmpCmp != rv_.cmp;
    if (ImGui::Button(TRL("Swap"))) { std::swap(rv_.ref, rv_.cmp); redo = true; }
    ImGui::SameLine(); redo |= ImGui::Checkbox(TRL("Field only"), &rv_.cmpFieldOnly);
    if (redo && rv_.ref != rv_.cmp) {
        rv_.diffs = ecl::compareSummaries(rv_.cases[size_t(rv_.ref)].sum(), rv_.cases[size_t(rv_.cmp)].sum());
        if (rv_.cmpFieldOnly) rv_.diffs.erase(std::remove_if(rv_.diffs.begin(), rv_.diffs.end(), [](const ecl::VectorDiff& d) { return d.key.find(':') != std::string::npos; }), rv_.diffs.end());
        rv_.cmpRef = rv_.ref; rv_.cmpCmp = rv_.cmp;
    }
    if (rv_.ref == rv_.cmp) { ui::hint("Choose two different cases."); return; }
    // summary line: worst field rate / pressure vector
    double worst = 0; std::string worstKey;
    for (auto& d : rv_.diffs)
        if (d.key.find(':') == std::string::npos && d.key.size() == 4 && d.key[0] == 'F' && (d.key[3] == 'R' || d.key == "FPR") && d.relRms > worst) { worst = d.relRms; worstKey = d.key; }
    ImGui::SameLine();
    if (!worstKey.empty())
        ImGui::TextColored(col4(worst < 0.02 ? Pal::Ok : worst < 0.05 ? Pal::Warn : Pal::Err), "worst field rate/pressure: %s %.2f %%", worstKey.c_str(), worst * 100);
    ImGui::SameLine();
    if (ImGui::SmallButton(TRL("Copy CSV"))) {
        std::ostringstream o;
        o << "key,reference_final,compared_final,max_abs_diff,rms_diff,rel_rms_percent,samples\n";
        for (auto& d : rv_.diffs) o << d.key << ',' << d.refFinal << ',' << d.cmpFinal << ',' << d.maxAbs << ',' << d.rms << ',' << d.relRms * 100 << ',' << d.samples << "\n";
        ImGui::SetClipboardText(o.str().c_str());
        notify("Comparison copied to the clipboard (CSV)");
    }
    ImGui::TextDisabled("%zu common vectors, sampled at the reference report times. relRMS = RMS(difference) / max |reference|. Click a row to plot it.", rv_.diffs.size());
    const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Sortable |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp;
    if (ui::beginTable("cmp", 6, tf)) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ui::tableColumn("Vector", ImGuiTableColumnFlags_DefaultSort, 1.6f);
        ui::tableColumn("Reference final", 0, 1.0f);
        ui::tableColumn("Compared final", 0, 1.0f);
        ui::tableColumn("Max |diff|", 0, 1.0f);
        ui::tableColumn("relRMS %", ImGuiTableColumnFlags_PreferSortDescending, 0.8f);
        ui::tableColumn("Samples", 0, 0.6f);
        ImGui::TableHeadersRow();
        if (ImGuiTableSortSpecs* ss = ImGui::TableGetSortSpecs()) {
            if (ss->SpecsDirty && ss->SpecsCount > 0) {
                const int col = ss->Specs[0].ColumnIndex; const bool asc = ss->Specs[0].SortDirection == ImGuiSortDirection_Ascending;
                std::stable_sort(rv_.diffs.begin(), rv_.diffs.end(), [&](const ecl::VectorDiff& a, const ecl::VectorDiff& b) {
                    double x = 0, y = 0;
                    switch (col) {
                        case 0: return asc ? a.key < b.key : a.key > b.key;
                        case 1: x = a.refFinal; y = b.refFinal; break;
                        case 2: x = a.cmpFinal; y = b.cmpFinal; break;
                        case 3: x = a.maxAbs; y = b.maxAbs; break;
                        case 4: x = a.relRms; y = b.relRms; break;
                        default: x = a.samples; y = b.samples; break;
                    }
                    return asc ? x < y : x > y;
                });
                ss->SpecsDirty = false;
            }
        }
        ImGuiListClipper clip;
        clip.Begin(int(rv_.diffs.size()));
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                const ecl::VectorDiff& d = rv_.diffs[size_t(r)];
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::PushID(r);
                if (ImGui::Selectable(d.key.c_str(), false, ImGuiSelectableFlags_SpanAllColumns)) {
                    if (std::find(rv_.selected.begin(), rv_.selected.end(), d.key) == rv_.selected.end()) rv_.selected.push_back(d.key);
                    rv_.wantTab = 0;
                }
                ImGui::PopID();
                ImGui::TableSetColumnIndex(1); ImGui::Text("%.6g", d.refFinal);
                ImGui::TableSetColumnIndex(2); ImGui::Text("%.6g", d.cmpFinal);
                ImGui::TableSetColumnIndex(3); ImGui::Text("%.4g", d.maxAbs);
                ImGui::TableSetColumnIndex(4);
                ImGui::TextColored(col4(d.relRms < 0.02 ? Pal::Ok : d.relRms < 0.05 ? Pal::Warn : Pal::Err), "%.2f", d.relRms * 100);
                ImGui::TableSetColumnIndex(5); ImGui::Text("%d", d.samples);
            }
        ui::endTable();
    }
}

// ---------------------------------------------------------------------------- 3D: model input / streamlines for result cases
bool App::rvEnsureInput(RVCase& c, bool wantSim) {
    if (!c.input && !c.deckTried && c.rc) {
        c.deckTried = true;
        const std::string base = c.rc->files.base;
        bool found = false;
        for (const char* ext : {".DATA", ".data", ".Data"}) {
            const std::string p = base + ext;
            std::ifstream f(p);
            if (!f) continue;
            found = true;
            try {
                std::vector<std::string> warn;
                c.input = std::make_shared<Case>(importEclipseDeck(p, &warn));
                c.simNote.clear();
            } catch (const std::exception& e) { c.simNote = "Could not import the deck " + fileNameOf(p) + ": " + e.what(); }
            break;
        }
        if (!found) c.simNote = "No deck (" + fileNameOf(base) + ".DATA) next to the result files - wells and streamlines need the model input.";
    }
    if (wantSim && !c.sim && c.input && c.simNote.empty()) {
        try {
            auto s = std::make_shared<Simulator>();
            s->setup(*c.input);
            s->initialize();
            const Grid* g = c.grid3d();
            if (g && s->grid().numActive() != g->numActive()) {
                char b[200];
                std::snprintf(b, sizeof(b), "The deck gives %d active cells but the results have %d (MINPV / PINCH treatment differs) - streamlines unavailable.",
                              s->grid().numActive(), g->numActive());
                c.simNote = b;
            } else c.sim = s;
        } catch (const std::exception& e) { c.simNote = std::string("Could not set up the model for streamlines: ") + e.what(); }
    }
    return c.input != nullptr;
}

void App::rvUpdateStreamlines(RVCase& c, int step) {
    SLState& sl = rv_.sl;
    ViewSettings& vs = rv_.view.settings();
    vs.showCells = !(sl.show && sl.hideCells);
    vs.linesXray = sl.xray;
    const bool needProp = rv_.prop == "TOF" || rv_.prop == "DRAINAGE" || rv_.prop == "SWEEP";
    auto clearLines = [&]() { if (sl.builtLines != -1) { rv_.view.setLines({}); sl.builtLines = -1; } };
    static const void* tracedCase = nullptr;
    static int tracedLines = -1;
    static std::vector<double> extFlux;     // external inter-block flows mapped to the grid connections
    static std::string fluxNote;
    if (tracedCase != &c) { sl.res.reset(); tracedCase = &c; }
    if ((!sl.show && !needProp) || !c.sim) { clearLines(); return; }
    if (c.internal && runner_.running()) return;   // the run's simulator may be busy
    bool rebuilt = false;
    if (!sl.res || sl.snap != step || tracedLines != sl.lines) {
        const auto t0 = std::chrono::steady_clock::now();
        try {
            if (c.rc) {
                const auto& st = c.rc->restart.steps();
                if (step < 0 || step >= int(st.size())) throw std::runtime_error("no restart step");
                Snapshot sn;
                sn.time = st[size_t(step)].days;
                sn.pressure = c.rc->restart.load(step, "PRESSURE");
                const size_t na = sn.pressure.size();
                auto get = [&](const char* n, float def) { std::vector<float> v = c.rc->restart.load(step, n); if (v.size() != na) v.assign(na, def); return v; };
                sn.sw = get("SWAT", 0.0f); sn.sg = get("SGAS", 0.0f); sn.rs = get("RS", 0.0f); sn.rv = get("RV", 0.0f);
                sn.so.resize(na);
                for (size_t a = 0; a < na; ++a) sn.so[a] = std::max(0.0f, 1.0f - sn.sw[a] - sn.sg[a]);
                { std::vector<float> p = c.rc->restart.load(step, "POLYMER"); if (p.size() == na) sn.poly = p; }
                // well BHPs at this time from the summary; wells without any rate are treated as shut
                const ecl::Summary& S = c.rc->summary;
                std::vector<double> bhp;
                for (const auto& wn : c.sim->wellNames()) {
                    auto val = [&](const char* kw) { const int v = S.find(std::string(kw) + ":" + wn); const double x = v >= 0 ? S.at(v, sn.time) : NAN; return std::isfinite(x) ? x : 0.0; };
                    const double rate = std::fabs(val("WOPR")) + std::fabs(val("WWPR")) + std::fabs(val("WGPR")) + std::fabs(val("WWIR")) + std::fabs(val("WGIR"));
                    const double b = val("WBHP");
                    bhp.push_back(rate > 0 && b > 0 ? b : -1.0);
                }
                c.sim->restoreExternal(sn, std::max(0, st[size_t(step)].report - 1), bhp);
                // inter-block flows of the simulator itself (RPTRST FLOWS): exact fluxes for the neighbour connections
                extFlux.clear(); fluxNote.clear();
                const Grid& G = c.sim->grid();
                std::vector<float> F[3];
                int phases = 0;
                for (int d = 0; d < 3; ++d) {
                    static const char* ph[] = {"OIL", "WAT", "GAS"};
                    for (const char* p : ph) {
                        std::vector<float> v = c.rc->restart.load(step, std::string("FLR") + p + "IJK"[d] + "+");
                        if (v.size() != na || int(na) != G.numActive()) continue;
                        if (F[d].empty()) F[d].assign(na, 0.0f);
                        for (size_t a = 0; a < na; ++a) F[d][a] += v[a];
                        if (d == 0) ++phases;
                    }
                }
                if (!F[0].empty() || !F[1].empty() || !F[2].empty()) {
                    extFlux.assign(G.conns.size(), NAN);
                    int used = 0;
                    for (size_t n = 0; n < G.conns.size(); ++n) {
                        const auto& cn = G.conns[n];
                        int ia, ja, ka, ib, jb, kb;
                        G.ijk(G.globalOf[cn.a], ia, ja, ka); G.ijk(G.globalOf[cn.b], ib, jb, kb);
                        const int di = ib - ia, dj = jb - ja, dk = kb - ka;
                        const int axis = (std::abs(di) == 1 && !dj && !dk) ? 0 : (!di && std::abs(dj) == 1 && !dk) ? 1 : (!di && !dj && std::abs(dk) == 1) ? 2 : -1;
                        if (axis < 0 || F[axis].empty() || cn.dir == 'M' || cn.dir == 'N') continue;
                        const bool aLow = (axis == 0 ? di : axis == 1 ? dj : dk) > 0;
                        const double q = F[axis][size_t(aLow ? cn.a : cn.b)];     // flow from the low cell to its + neighbour
                        extFlux[n] = aLow ? q : -q;
                        ++used;
                    }
                    fluxNote = "streamlines from the simulator's inter-block flows (" + std::to_string(phases) + " phase(s), " +
                               std::to_string(used) + " connections; NNCs from the pressure field)";
                } else fluxNote = "streamlines from the restart pressure field (add RPTRST FLOWS to the deck for the simulator's own fluxes)";
            } else c.sim->restoreReport(step);
            StreamlineOptions o; o.lines = sl.lines;
            if (c.rc && !extFlux.empty()) o.connFlux = &extFlux;
            sl.source = c.rc ? fluxNote : "streamlines from the ResSim connection fluxes";
            sl.res = std::make_shared<StreamlineResult>(traceStreamlines(*c.sim, o));
        } catch (const std::exception& e) {
            c.simNote = std::string("Streamline tracing failed: ") + e.what();
            sl.res = std::make_shared<StreamlineResult>();
        }
        sl.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        sl.snap = step; tracedLines = sl.lines;
        rv_.key3d.clear();   // refresh TOF / region colours
        rebuilt = true;
    }
    if (!sl.show) { clearLines(); return; }
    const auto& W = rv_.view.wells();
    const std::string focus = vs.focusWell >= 0 && vs.focusWell < int(W.size()) ? W[size_t(vs.focusWell)].name : std::string();
    const std::string fkey = focus + (sl.focusOnly ? "|1" : "|0");
    if (!rebuilt && sl.builtLines == int(sl.res->lines.size()) && sl.builtColor == sl.colorBy + 10 * int(vs.cmap) && sl.builtFocus == fkey) return;
    rv_.view.setLines(streamlineVertices(*sl.res, sl.colorBy, vs.cmap, focus, sl.focusOnly));
    sl.builtLines = int(sl.res->lines.size());
    sl.builtColor = sl.colorBy + 10 * int(vs.cmap);
    sl.builtFocus = fkey;
}

// C30: a restart array of an external case, or its post-processed geomechanics array
std::vector<float> App::rvLoad(RVCase& c, int step, const std::string& name) {
    if (step >= 0 && size_t(step) < c.post.size())
        for (const auto& e : c.post[size_t(step)].extra) if (e.first == name) return e.second;
    return c.rc ? c.rc->restart.load(step, name) : std::vector<float>{};
}

// ---------------------------------------------------------------------------- 3D
void App::rv3D() {
    const float fs = ImGui::GetFontSize();
    if (!rv_.viewInit) { rv_.viewInit = true; rv_.viewOk = rv_.view.init(); }
    if (!rv_.viewOk) { ImGui::TextColored(col4(Pal::Err), "3D renderer failed to initialise."); return; }
    std::vector<int> with;
    for (int i = 0; i < int(rv_.cases.size()); ++i) if (rv_.cases[size_t(i)].grid3d()) with.push_back(i);
    if (with.empty()) { ui::hint("No case with a grid: open a case that has an EGRID file (or add a ResSim run)."); return; }
    if (std::find(with.begin(), with.end(), rv_.c3d) == with.end()) rv_.c3d = with[0];
    RVCase& C = rv_.cases[size_t(rv_.c3d)];
    ImGui::SetNextItemWidth(fs * 16);
    if (ImGui::BeginCombo("##c3d", C.name.c_str())) {
        for (int i : with) if (ImGui::Selectable(rv_.cases[size_t(i)].name.c_str(), rv_.c3d == i)) { rv_.c3d = i; rv_.step = 0; rv_.view.settings().focusWell = -1; }
        ImGui::EndCombo();
    }
    RVCase& c = rv_.cases[size_t(rv_.c3d)];
    rvEnsureInput(c, rv_.sl.show || rv_.prop == "TOF" || rv_.prop == "DRAINAGE" || rv_.prop == "SWEEP");
    // property lists: static (INIT) and dynamic (restart)
    std::vector<std::string> statics, dyn;
    int nsteps = 0;
    if (c.rc) {
        statics = c.rc->initArrays;
        if (c.rc->hasRestart) { nsteps = int(c.rc->restart.steps().size()); if (nsteps) dyn = c.rc->restart.arrayNames(std::min(rv_.step, nsteps - 1)); }
        if (nsteps && !c.post.empty()) for (const auto& e : c.post[size_t(std::min(int(c.post.size()) - 1, std::min(rv_.step, nsteps - 1)))].extra) dyn.push_back(e.first);   // C30
    } else {
        nsteps = int(c.snaps.size());
        if (nsteps) { dyn = {"PRESSURE", "SWAT", "SOIL", "SGAS", "RS"}; if (!c.snaps.front().rv.empty()) dyn.push_back("RV"); if (!c.snaps.front().poly.empty()) { dyn.push_back("POLYMER"); dyn.push_back("PADS"); } if (!c.snaps.front().xco2.empty()) dyn.push_back("XCO2"); for (const auto& e : c.snaps.back().extra) dyn.push_back(e.first); }
        statics = {"DEPTH", "PORV"};
    }
    if (c.sim && nsteps > 0) { dyn.push_back("TOF"); dyn.push_back("DRAINAGE"); dyn.push_back("SWEEP"); }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(fs * 11);
    if (ImGui::BeginCombo("##p3d", rv_.prop.c_str(), ImGuiComboFlags_HeightLarge)) {
        if (!statics.empty()) ImGui::SeparatorText(TR("Static (INIT)"));
        for (auto& s : statics) if (ImGui::Selectable(("INIT " + s).c_str(), rv_.prop == "INIT " + s)) { rv_.prop = "INIT " + s; rv_.autoRange = true; }
        if (!dyn.empty()) ImGui::SeparatorText(TR("Dynamic (restart)"));
        for (auto& s : dyn) {
            if (s == "TOF") ImGui::SeparatorText(TR("Streamline (traced)"));
            if (ImGui::Selectable(s.c_str(), rv_.prop == s)) { rv_.prop = s; rv_.autoRange = true; }
        }
        ImGui::EndCombo();
    }
    if (nsteps > 0) {   // C30: post-processed FE geomechanics of this case (ECLIPSE / OPM restart or ResSim run)
        ImGui::SameLine();
        if (ImGui::Button(TRL(ICON_LAYERS " Geomechanics##rvgeo"))) ImGui::OpenPopup("rvgeo");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Compute stress / strain / displacement / failure maps of this case from its pressure (and\n"
                                                      "temperature) changes - one-way FE post-processing, any engine"));
        if (ImGui::BeginPopup("rvgeo")) {
            std::shared_ptr<const Grid> gg = c.rc ? std::shared_ptr<const Grid>(c.rc->grid) : c.grid;
            const bool fu = c.rc ? c.rc->unitSystem == 2 : (c.input && c.input->units == UnitSystem::Field);
            ImGui::PushItemWidth(fs * 10);
            postGeoControls(1, gg, fu, nullptr);
            ImGui::PopItemWidth();
            ImGui::EndPopup();
        }
    }
    const bool isStatic = rv_.prop.rfind("INIT ", 0) == 0;
    const bool slProp = rv_.prop == "TOF" || rv_.prop == "DRAINAGE" || rv_.prop == "SWEEP";
    static bool play = false; static double playAcc = 0;
    if (nsteps > 0) {
        ImGui::SameLine();
        if (ImGui::Button(play ? ICON_PAUSE "##rvp" : ICON_PLAY "##rvp")) { play = !play; if (play && rv_.step >= nsteps - 1) rv_.step = 0; }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(play ? "Pause" : "Play through the report steps");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 14);
        rv_.step = std::max(0, std::min(rv_.step, nsteps - 1));
        const double days = c.rc ? c.rc->restart.steps()[size_t(rv_.step)].days : c.snaps[size_t(rv_.step)].time;
        char fmt[64]; std::snprintf(fmt, sizeof(fmt), "step %%d  (%.0f d)", days);
        ImGui::SliderInt("##st3d", &rv_.step, 0, nsteps - 1, fmt);
        if (play) {
            playAcc += ImGui::GetIO().DeltaTime;
            if (playAcc > (rv_.sl.show ? 0.4 : 0.15)) { playAcc = 0; if (++rv_.step >= nsteps) { rv_.step = nsteps - 1; play = false; } }
        }
    } else play = false;
    const Grid* g = c.grid3d();
    // wells from the model input (ResSim run / in-app OPM run / deck next to the results), status at this step
    {
        int schedStep = 0;
        if (nsteps > 0) schedStep = c.rc ? std::max(0, c.rc->restart.steps()[size_t(rv_.step)].report - 1) : std::max(0, rv_.step - 1);
        char wk[128]; std::snprintf(wk, sizeof(wk), "%p|%p|%d", (const void*)g, (const void*)c.input.get(), schedStep);
        static std::string wellKey;
        if (wellKey != wk) {
            wellKey = wk;
            if (c.input && g && c.input->grid.nx == g->nx && c.input->grid.ny == g->ny && c.input->grid.nz == g->nz) rv_.view.setWells(wellGlyphs(*g, *c.input, schedStep));
            else rv_.view.setWells({});
            if (rv_.view.settings().focusWell >= int(rv_.view.wells().size())) rv_.view.settings().focusWell = -1;
        }
    }
    {   // C13: EDFM fracture planes (run segments coloured by the fracture value of the property, else input planes)
        const std::string pn = isStatic ? std::string() : rv_.prop;
        char fk[256];
        std::snprintf(fk, sizeof(fk), "%p|%p|%d|%s|%g|%g|%d", (const void*)g, (const void*)c.input.get(), rv_.step, pn.c_str(), double(rv_.rmin),
                      double(rv_.rmax), rv_.autoRange ? 1 : 0);
        if (fracKeyRv_ != fk) {
            fracKeyRv_ = fk;
            const bool seg = !c.rc && !c.egs.fracSegs.empty() && c.egs.fracSnap.size() == c.snaps.size();
            if (c.input && g) rv_.view.setPanels(fracturePanels(*c.input, g, seg ? &c.egs : nullptr, rv_.step, pn));
            else rv_.view.setPanels({});
        }
    }
    ViewerBinding vb;
    vb.view = &rv_.view; vb.autoRange = &rv_.autoRange; vb.rmin = &rv_.rmin; vb.rmax = &rv_.rmax; vb.sl = &rv_.sl;
    if (!isStatic && nsteps > 1 && !slProp) vb.rangeAllSteps = &settings_.rangeAllSteps;
    vb.id = "results"; vb.propName = isStatic ? rv_.prop.substr(5) : rv_.prop; vb.grid = g;
    vb.pngName = c.name + "_" + vb.propName;
    vb.slAvailable = c.sim != nullptr && nsteps > 0;
    vb.slHint = c.sim ? (nsteps > 0 ? "" : "No report steps in these results.")
                      : !c.simNote.empty() ? c.simNote
                      : c.rc ? "Streamlines need the deck (CASE.DATA) next to the result files and the restart arrays PRESSURE / SWAT / SGAS."
                             : "This run has no simulator attached.";
    if (!c.sim && c.input && c.simNote.empty()) vb.slPrepare = [this, &c]() { rvEnsureInput(c, true); };
    if (!isStatic && !slProp && nsteps > 0) {
        vb.cellHistory = [this, &c, nsteps](int cell, std::vector<double>& t, std::vector<double>& v) {
            for (int s = 0; s < nsteps; ++s) {
                if (c.rc) {
                    const std::vector<float> a = rvLoad(c, s, rv_.prop);
                    if (size_t(cell) >= a.size()) continue;
                    t.push_back(c.rc->restart.steps()[size_t(s)].days); v.push_back(a[size_t(cell)]);
                } else {
                    const Snapshot& sn = c.snaps[size_t(s)];
                    const std::string& p = rv_.prop;
                    const std::vector<float>& a = sn.field(p);
                    if (size_t(cell) >= a.size()) continue;
                    t.push_back(sn.time); v.push_back(a[size_t(cell)]);
                }
            }
            return !t.empty();
        };
    }
    viewerToolbar(vb, rvUi_);
    rvUpdateStreamlines(c, rv_.step);

    // (re)load values when case / property / step / range change
    char kb[256];
    std::snprintf(kb, sizeof(kb), "%p|%s|%d|%d|%g|%g|%d|%d", (const void*)g, rv_.prop.c_str(), isStatic ? -1 : rv_.step, rv_.autoRange ? 1 : 0, rv_.rmin, rv_.rmax, rv_.view.settings().logScale ? 1 : 0, settings_.rangeAllSteps ? 1 : 0);
    if (rv_.key3d != kb) {
        if (g != rv_.grid3d) { rv_.view.setGrid(g); rv_.grid3d = g; rv_.view.resetCamera(); }
        std::vector<float> v;
        std::string units;
        if (slProp) {
            if (rv_.sl.res && rv_.sl.snap == rv_.step) {
                const StreamlineResult& S = *rv_.sl.res;
                v.assign(S.cellTof.size(), NAN);
                for (size_t a = 0; a < v.size(); ++a) {
                    if (rv_.prop == "TOF") v[a] = S.cellTof[a] >= 0 ? S.cellTof[a] : NAN;
                    else if (rv_.prop == "DRAINAGE") v[a] = a < S.cellSink.size() && S.cellSink[a] >= 0 ? float(S.cellSink[a]) : NAN;
                    else v[a] = a < S.cellSource.size() && S.cellSource[a] >= 0 ? float(S.cellSource[a]) : NAN;
                }
            } else if (!rv_.sl.res || rv_.sl.snap != rv_.step) {
                // trace on demand (streamlines need not be displayed for the TOF / region properties)
                rvUpdateStreamlines(c, rv_.step);
            }
            units = rv_.prop == "TOF" ? "days" : "well index";
        } else if (c.rc) {
            if (isStatic) v = c.rc->initArray(rv_.prop.substr(5));
            else if (nsteps > 0) v = rvLoad(c, rv_.step, rv_.prop);
        } else {
            if (isStatic) {
                const std::string n = rv_.prop.substr(5);
                v.resize(size_t(g->numActive()));
                for (int a = 0; a < g->numActive(); ++a) v[size_t(a)] = float(n == "DEPTH" ? g->depth[size_t(a)] : g->poreVolRef[size_t(a)]);
            } else if (nsteps > 0) {
                const Snapshot& s = c.snaps[size_t(rv_.step)];
                v = s.field(rv_.prop);
            }
        }
        if (v.size() != size_t(g->numActive())) v.assign(size_t(g->numActive()), NAN);
        if (rv_.autoRange) {
            float lo = FLT_MAX, hi = -FLT_MAX;
            for (float x : v) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); }
            // C18.1: one colour scale over every report step (cached per case + property; ECLIPSE / OPM restart steps are read once)
            if (!isStatic && !slProp && nsteps > 1 && settings_.rangeAllSteps) {
                char gk[200];
                std::snprintf(gk, sizeof(gk), "%p|%s|%d", (const void*)&c, rv_.prop.c_str(), nsteps);
                if (rv_.grKey != gk) {
                    float glo = FLT_MAX, ghi = -FLT_MAX;
                    auto acc = [&](const std::vector<float>& a) { for (float x : a) if (std::isfinite(x)) { glo = std::min(glo, x); ghi = std::max(ghi, x); } };
                    for (int s = 0; s < nsteps; ++s) {
                        if (c.rc) { try { acc(rvLoad(c, s, rv_.prop)); } catch (...) {} }
                        else acc(c.snaps[size_t(s)].field(rv_.prop));
                    }
                    rv_.grKey = gk; rv_.grLo = glo; rv_.grHi = ghi; rv_.grOk = glo <= ghi;
                }
                if (rv_.grOk) { lo = rv_.grLo; hi = rv_.grHi; }
            }
            if (lo > hi) { lo = 0; hi = 1; }
            if (rv_.prop == "SWAT" || rv_.prop == "SOIL" || rv_.prop == "SGAS") { lo = 0; hi = 1; }
            if (hi - lo < 1e-12f) hi = lo + 1e-6f;
            rv_.rmin = lo; rv_.rmax = hi;
        }
        std::snprintf(kb, sizeof(kb), "%p|%s|%d|%d|%g|%g|%d|%d", (const void*)g, rv_.prop.c_str(), isStatic ? -1 : rv_.step, rv_.autoRange ? 1 : 0, rv_.rmin, rv_.rmax, rv_.view.settings().logScale ? 1 : 0, settings_.rangeAllSteps ? 1 : 0);
        rv_.key3d = (slProp && !(rv_.sl.res && rv_.sl.snap == rv_.step)) ? std::string() : std::string(kb);
        rv_.vals = v;
        rv_.view.setValues(v, rv_.rmin, rv_.rmax, isStatic ? rv_.prop.substr(5) : rv_.prop, units);
        std::string cap = c.name;
        if (!isStatic && nsteps > 0) {
            double days = c.rc ? c.rc->restart.steps()[size_t(rv_.step)].days : c.snaps[size_t(rv_.step)].time;
            char b[96]; std::snprintf(b, sizeof(b), "  |  t = %.1f days (%.2f years)", days, days / 365.25); cap += b;
        }
        if (!c.input && !c.simNote.empty()) cap += "  |  no wells (deck not found)";
        rv_.view.setCaption(cap);
    }
    viewerCanvas(vb, rvUi_, [this, g](int cell) {
        ImGui::Text("Cell %s  active #%d", g->cellLabel(g->globalOf[size_t(cell)]).c_str(), cell + 1);
        ImGui::Separator();
        ImGui::Text("%s = %.6g", rv_.prop.c_str(), rv_.view.valueOf(cell));
        ImGui::TextDisabled("depth %.2f", g->depth[size_t(cell)]);
    });
}

// ---------------------------------------------------------------------------- info
void App::rvInfo() {
    for (auto& c : rv_.cases) {
        if (!ImGui::CollapsingHeader(c.name.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) continue;
        ImGui::PushID(&c);
        const ecl::Summary& S = c.sum();
        if (ui::beginProps("inf", 0.3f)) {
            if (c.rc) {
                const auto& F = c.rc->files;
                ui::propText("Case", F.base.c_str());
                ui::propText("Grid (EGRID)", F.egrid.empty() ? "-" : fileNameOf(F.egrid).c_str());
                ui::propText("INIT", F.init.empty() ? "-" : fileNameOf(F.init).c_str());
                ui::propText("Summary", F.smspec.empty() ? "-" : fileNameOf(F.smspec).c_str());
                std::string rs = F.restart.empty() ? "-" : fileNameOf(F.restart[0]) + (F.restart.size() > 1 ? " ... (" + std::to_string(F.restart.size()) + " files)" : "");
                ui::propText("Restart", rs.c_str());
                char b[160];
                if (c.rc->hasGrid) { std::snprintf(b, sizeof(b), "%d x %d x %d, %d active cells", c.rc->egrid.nx, c.rc->egrid.ny, c.rc->egrid.nz, c.rc->egrid.numActive()); ui::propText("Dimensions", b); }
                static const char* us[] = {"unknown", "METRIC", "FIELD", "LAB", "PVT-M"};
                ui::propText("Units", us[std::max(0, std::min(4, c.rc->unitSystem))]);
                std::snprintf(b, sizeof(b), "%zu", c.rc->hasRestart ? c.rc->restart.steps().size() : size_t(0)); ui::propText("Restart steps", b);
                std::snprintf(b, sizeof(b), "%zu INIT arrays", c.rc->initArrays.size()); ui::propText("Static data", b);
            } else ui::propText("Source", "ResSim run (in memory)");
            char b[160];
            std::snprintf(b, sizeof(b), "%zu vectors, %zu time steps, %.1f days", S.vars.size(), S.time.size(), S.time.empty() ? 0.0 : S.time.back());
            ui::propText("Summary", b);
            std::snprintf(b, sizeof(b), "%02d.%02d.%04d", S.startDay, S.startMonth, S.startYear); ui::propText("Start", b);
            ui::endProps();
        }
        if (c.rc) for (auto& w : c.rc->warnings) ImGui::TextColored(col4(Pal::Warn), "%s", w.c_str());
        ImGui::PopID();
    }
    if (rv_.cases.empty()) ui::hint("No cases open.");
}
