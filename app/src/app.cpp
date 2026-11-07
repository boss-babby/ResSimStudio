// App: frame orchestration, docking layout, menus, toolbar, status bar, case explorer,
// inspector dispatch, log, file actions, run control and initialisation preview.
#include "app.h"
#include "i18n.h"
#include "nav.h"
#include <filesystem>
#include "plotkit.h"
#include "ressim/therm_deck.h"
#include "icons.h"
#include "widgets.h"
#include "theme.h"
#include "gfx.h"
#include "imgui_internal.h"
#include "implot.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/examples.h"
#include "ressim/comp_deck.h"
#include "platform.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <thread>
#include <cstring>
#include <vector>
#include <cctype>
#include <stdexcept>

using namespace ressim;

static const char* kSectionNames[] = {TRN("Run specification"), TRN("Grid & rock"), TRN("Fluid (PVT / EOS / CO2)"), TRN("Rel. perm & Pc"), TRN("Initialization"), TRN("Aquifers"), TRN("EOR (polymer)"), TRN("Wells"), TRN("Schedule"), TRN("Numerics"), TRN("Advanced keywords")};
static const char* kSectionIcons[] = {ICON_SETTINGS, ICON_GRID_3X3, ICON_DROPLET, ICON_CHART_SPLINE, ICON_ARROW_DOWN_TO_LINE, ICON_WAVES, ICON_FLASK_CONICAL, ICON_CIRCLE_DOT, ICON_HISTORY, ICON_CPU, ICON_LIST};
static const char* kSectionHelp[] = {
    "Title, unit system and active phases (black-oil formulation).",
    "Cartesian or corner-point (COORD/ZCORN) grid, layering, heterogeneous permeability/porosity, faults (MULTFLT) and NNCs.",
    "Formation volume factors, solution gas and viscosities vs pressure; or generate from correlations.",
    "Saturation functions (SWOF/SGOF): relative permeability and capillary pressure.",
    "Hydrostatic equilibrium: datum pressure and fluid contacts. Preview the initial state and fluids in place.",
    "Analytic aquifers (Fetkovich pseudo-steady state, Carter-Tracy infinite acting) connected to grid faces.",
    "Enhanced oil recovery: polymer flooding (viscosity, adsorption, mixing regions, shear, injection schedule).",
    "Well locations and completions. Well index from Peaceman's model.",
    "Report steps and time-dependent well controls (rates, BHP, history).",
    "Time-step control, Newton and linear solver tolerances and CPU threads.",
    "Any other ECLIPSE / OPM Flow keyword (groups, VFP tables, multisegment wells, ACTIONX/UDQ, tracers, summary requests, ...). "
    "Kept verbatim, edited here and written to the deck in place - used by OPM Flow; the ResSim engine ignores them."};

static std::string baseName(const std::string& p) {
    const size_t s = p.find_last_of("\\/");
    return s == std::string::npos ? p : p.substr(s + 1);
}
static bool endsWithI(const std::string& s, const char* suf) {
    const size_t n = std::strlen(suf);
    if (s.size() < n) return false;
    for (size_t i = 0; i < n; ++i) if (std::tolower((unsigned char)s[s.size() - n + i]) != std::tolower((unsigned char)suf[i])) return false;
    return true;
}
// Windows file names may not contain \\/:*?"<>| or control characters, nor end in a dot or space.
std::string safeFileName(const std::string& name) {
    std::string r;
    for (unsigned char ch : name) {
        if (ch < 32 || std::strchr("\\/:*?\"<>|", ch)) { if (r.empty() || r.back() != '_') r += '_'; }
        else r += char(ch);
    }
    while (!r.empty() && (r.back() == '.' || r.back() == ' ')) r.pop_back();
    while (!r.empty() && r.front() == ' ') r.erase(r.begin());
    if (r.size() > 120) r.resize(120);
    static const char* reserved[] = {"CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4", "LPT1", "LPT2", "LPT3"};
    std::string up = r.substr(0, r.find('.'));
    for (char& ch : up) ch = char(std::toupper((unsigned char)ch));
    for (const char* x : reserved) if (up == x) { r = "_" + r; break; }
    return r.empty() ? std::string("untitled") : r;
}

// ============================================================================ lifecycle
void App::init() {
    viewOk_ = view_.init();
    loadSettings();
    applyHpcSettings();
    loadCustomPlots();
    if (settings_.uiScale <= 0) {
        // first start: laptop screens (1366 / 1440 / 1536 logical px) get a compact interface
        const float logical = float(plat::screenWidth() > 0 ? plat::screenWidth() : 1920) / std::max(0.5f, g_gfx.dpiScale);
        settings_.uiScale = logical < 1500 ? 0.85f : logical < 1800 ? 0.92f : 1.0f;
    }
    setUiScale(settings_.uiScale);
    // C14: every plot can be exported on its own (right-click menu / corner button) through the figure export window
    installPlotCapture([this](const fig::Figure& f) {
        figx_.fig = f;
        figx_.plotId = lastCapturedPlot();
        figx_.snapWait = false;
        figx_.id.clear();                 // captured single plot (not a recorded dashboard figure)
        figx_.open = true;
        figx_.cols = 0;
        figDeferredFocus_ = true;
    });
    newCase(0);
    showStart_ = settings_.showStartPage;
    applyColors(settings_.theme, settings_.accent);   // C39 colour theme
    recoveryCheck();                                  // C39 autosave left by a crash?
    runner_.addLog("ResSim Studio started. " + std::to_string(std::max(1u, std::thread::hardware_concurrency())) + " hardware threads available.");
    if (!viewOk_) runner_.addLog("ERROR: 3D renderer failed to initialise.");
}

void App::shutdown() {
    runner_.cancel();
    runner_.join();
    study_.cancel = true;
    if (study_.th.joinable()) study_.th.join();
    if (study_.manualTh.joinable()) study_.manualTh.join();
    view_.shutdown();
    if (rv_.viewInit) rv_.view.shutdown();
    saveSettings();
    if (nav::indexChanged()) nav::saveIndex(plat::configDir() + "/nav_index.txt");
    recoveryClear();   // C39: closed normally (changes were saved or discarded in confirmClose)
}

bool App::wantsContinuousFrames() const {
    return runner_.running() || study_.running || study_.manualPending || play_ || ImGui::GetTime() < toastUntil_ || nav::pending();
}

void App::frame() {
    plotCaptureNewFrame();   // C22: hand over single-plot captures / exact snapshots of the previous frame
    nav::beginFrame(int(workspace_), int(section_));   // C40 deep links + tab index
    if (!smokeAction_.empty() && ++smokeFrame_ == 3) {   // C40 scripted smoke-test states
        if (smokeAction_ == "wizard") showWizard_ = true;
        else if (smokeAction_ == "palette") { paletteOpen_ = true; paletteFocus_ = true; }
        else if (smokeAction_ == "run" || smokeAction_ == "results") runSimulation();
        else if (smokeAction_ == "model") workspace_ = Workspace::Model;   // C45 (translated panels)
        else if (smokeAction_ == "cases" || smokeAction_ == "filter") {   // C40: case manager with 3 variants / Inspector filter
            if (smokeAction_ == "cases") {
                casesAdd("Base case", case_);
                Case v = case_; for (auto& x : v.grid.permx) x *= 2.0; v.title += " (k x2)"; casesAdd("High permeability", v);
                variants_.back().note = "PERMX doubled";
                Case w = case_; w.equil.owc += 30; casesAdd("Deeper OWC", w);
                casesSel_ = 1; showCases_ = true;
            } else { section_ = Section::Grid; inspFilter_ = "perm"; }
        }
    }
    if (smokeAction_ == "results" && smokeFrame_ > 3 && !runner_.running()) workspace_ = Workspace::Results;   // C45
    if ((smokeAction_ == "run" || smokeAction_ == "results") && smokeFrame_ > 3 && smokeFrame_ % 20 == 0 && !runner_.running()) {   // finished: report step to mid-run
        std::lock_guard<std::mutex> l(runner_.mutex());
        const int n = int(runner_.results().snaps.size());
        if (runner_.results().valid && n >= 2) snap_ = n / 2;
    }
    {   // C40 linked time cursor: the current report time on every time plot
        double t = NAN, epoch = 0;
        if (settings_.timeCursor && !runner_.running()) {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& R = runner_.results();
            if (R.valid && R.snaps.size() >= 2) {
                t = R.snaps[size_t(std::clamp(snap_, 0, int(R.snaps.size()) - 1))].time;
                // days from 1970-01-01 to START (civil-from-days inverse, proleptic Gregorian)
                int y = R.input.startYear, m = R.input.startMonth, d = R.input.startDay;
                y -= m <= 2;
                const int era = (y >= 0 ? y : y - 399) / 400, yoe = y - era * 400;
                const int doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1, doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
                epoch = (double(era) * 146097.0 + doe - 719468.0) * 86400.0;
            }
        }
        timeCursorSet(t, epoch, settings_.timeCursor);
    }
    {
        static bool navLoaded = false;
        static double navSaveAt = 0;
        const std::string navFile = plat::configDir() + "/nav_index.txt";
        if (!navLoaded) { navLoaded = true; nav::loadIndex(navFile); }
        if (nav::indexChanged() && ImGui::GetTime() > navSaveAt) { nav::saveIndex(navFile); navSaveAt = ImGui::GetTime() + 5.0; }
    }
    if (gridDirty_) rebuildPreview();

    // run completion
    const bool running = runner_.running();
    if (wasRunning_ && !running) {
        runner_.join();
        const std::string err = runner_.error();
        std::string opmDeck, extName;
        bool opmOk = false;
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& R = runner_.results();
            snap_ = std::max(0, int(R.snaps.size()) - 1);
            if (R.engine != "ResSim") { opmDeck = R.opmDeck; opmOk = R.valid; extName = R.engine; }
        }
        if (err.empty() && opmDeck.empty()) { afterRunOutputs(); if (settings_.runHistory) rvAddRunHistory(); }
        if (err.empty()) notify(opmDeck.empty() ? std::string(TR("Simulation finished")) : extName + TR(" run finished"));
        else notify("Run stopped: " + err);
        // decks ResSim cannot represent: their 3-D results only map onto the engine's own grid, so always open them
        // in the result viewer (the Model grid holds placeholders)
        if (opmOk && !opmDeck.empty() && (settings_.opmAutoView || !case_.engineUnsupported.empty())) rvOpen(opmDeck, true);
    }
    wasRunning_ = running;
    pollPostGeo();   // C30 post-processed geomechanics finished?

    // global shortcuts (the Deck Editor handles its own file shortcuts and F5)
    const ImGuiInputFlags G = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_1, G)) workspace_ = Workspace::Model;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_2, G)) workspace_ = Workspace::DeckEditor;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_3, G)) workspace_ = Workspace::Results;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_4, G)) workspace_ = Workspace::Studies;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_5, G)) workspace_ = Workspace::Geomodel;
    if (geoMode_) workspace_ = Workspace::Geomodel;   // ResGeo: the geomodelling front end only
    if (ImGui::Shortcut(ImGuiMod_Shift | ImGuiKey_F5, G)) runner_.cancel();
    if (workspace_ != Workspace::DeckEditor) {
        if (ImGui::Shortcut(ImGuiKey_F5, G)) runSimulation();
        if (ImGui::Shortcut(ImGuiKey_F1, G)) openUserGuide();
        {   // interface size: Ctrl + / Ctrl - / Ctrl 0
            float k = settings_.uiScale;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Equal, G) || ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_KeypadAdd, G)) k += 0.05f;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Minus, G) || ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_KeypadSubtract, G)) k -= 0.05f;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_0, G)) k = 1.0f;
            k = std::clamp(k, 0.6f, 1.6f);
            if (std::fabs(k - settings_.uiScale) > 1e-4f) { settings_.uiScale = k; setUiScale(k); saveSettings(); }
        }
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_S, G)) save(false);
        // C39 undo / redo of model edits (text fields keep their own Ctrl+Z while typing; the Geomodel has its own history)
        if (workspace_ != Workspace::Geomodel && !ImGui::GetIO().WantTextInput) {
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Z, G | ImGuiInputFlags_Repeat)) undoApply(-1);
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Y, G | ImGuiInputFlags_Repeat) ||
                ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_Z, G | ImGuiInputFlags_Repeat)) undoApply(+1);
        }
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_S, G)) save(true);
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_O, G)) openDialog();
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_N, G)) { if (maybeSaveChanges()) newCase(2); }
    }

    menuBar();
    toolbar();
    statusBar();
    commandPalette();   // C26: Ctrl+Shift+P / Ctrl+K
    wfTick();           // C41 workflow replay (one step per frame)
    if (showWorkflow_) panelWorkflow();
    if (showCases_) panelCases();   // C40 (G10)
    startPage();
    if (workspace_ == Workspace::Model) {
        dockSpace();
        panelCase();
        panelEditor();
        panelViewer();
        panelPlots();
        panelRun();
        panelLog();
        panelStudies();
        panelDeckTools();
        panelResultViewer();
    } else {
        // keep the docking layout (and the docked windows) alive while another workspace is shown
        ImGui::DockSpace(ImGui::GetID("MainDock"), ImVec2(0, 0), ImGuiDockNodeFlags_KeepAliveOnly);
        if (workspace_ == Workspace::DeckEditor) deckEditorWorkspace();
        else if (workspace_ == Workspace::Studies) studiesWorkspace();
        else if (workspace_ == Workspace::Geomodel) geomodelWorkspace();
        else resultsWorkspace();
    }
    figureExportWindow();
    opmSettingsWindow();
    restartWindow();
    reportWindow();
    regressionWindow();
    stabilityWindow();
    seismicityWindow();
    efemWindow();
    if (rvDeferred_) { auto f = std::move(rvDeferred_); rvDeferred_ = nullptr; f(); }
    aboutWindow();
    newModelWizard();   // C40
    recoveryWindow();
    if (showImGuiDemo_) ImGui::ShowDemoWindow(&showImGuiDemo_);
    undoTrack();      // C39: snapshot after a finished edit gesture
    autosaveTick();   // C39: background recovery copy
    {   // C40: the time cursor was dragged on a plot -> nearest report step for every view
        double d = 0;
        if (timeCursorDragged(&d)) {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& R = runner_.results();
            int best = snap_; double bd = 1e300;
            for (int k = 0; k < int(R.snaps.size()); ++k) { const double e = std::fabs(R.snaps[size_t(k)].time - d); if (e < bd) { bd = e; best = k; } }
            snap_ = best; play_ = false;
        }
    }

    // toast notification (bottom-right, fades out)
    const double now = ImGui::GetTime();
    if (now < toastUntil_ && !toast_.empty()) {
        const ImGuiViewport* vp = ImGui::GetMainViewport();
        const float alpha = float(std::min(1.0, (toastUntil_ - now) / 0.4));
        ImDrawList* dl = ImGui::GetForegroundDrawList();
        const ImVec2 ts = ImGui::CalcTextSize(toast_.c_str());
        const float pad = ImGui::GetFontSize() * 0.8f;
        const ImVec2 p1(vp->WorkPos.x + vp->WorkSize.x - pad * 2, vp->WorkPos.y + vp->WorkSize.y - pad * 2);
        const ImVec2 p0(p1.x - ts.x - pad * 2 - 6, p1.y - ts.y - pad * 2);
        dl->AddRectFilled(ImVec2(p0.x + 3, p0.y + 4), ImVec2(p1.x + 3, p1.y + 4), IM_COL32(0, 0, 0, int(90 * alpha)), 8);
        { ImVec4 cb = themeCardBg(); cb.w *= alpha; dl->AddRectFilled(p0, p1, ImGui::GetColorU32(cb), 8); }
        if (themeIsLight()) dl->AddRect(p0, p1, ImGui::GetColorU32(ImGuiCol_Border, alpha), 8);
        dl->AddRectFilled(p0, ImVec2(p0.x + 5, p1.y), (Pal::Accent & 0x00FFFFFF) | (ImU32(255 * alpha) << 24), 8, ImDrawFlags_RoundCornersLeft);
        dl->AddText(ImVec2(p0.x + pad + 6, p0.y + pad), ImGui::GetColorU32(ImGuiCol_Text, alpha), toast_.c_str());
    }

    // window title
    std::string title = "ResSim Studio - " + (path_.empty() ? case_.title : baseName(path_)) + (dirty_ ? " *" : "");
    if (std::string(titleBuf_.begin(), titleBuf_.end()) != title) {
        titleBuf_.assign(title.begin(), title.end());
        plat::setWindowTitle(title);
    }
}

// ============================================================================ layout
void App::dockSpace() {
    ImGuiViewport* vp = ImGui::GetMainViewport();
    const ImGuiID id = ImGui::GetID("MainDock");
    if (resetLayout_ || ImGui::DockBuilderGetNode(id) == nullptr) {
        resetLayout_ = false;
        ImGui::DockBuilderRemoveNode(id);
        ImGui::DockBuilderAddNode(id, ImGuiDockNodeFlags_DockSpace);
        ImGui::DockBuilderSetNodeSize(id, vp->WorkSize);
        ImGuiID center = id;
        const ImGuiID left = ImGui::DockBuilderSplitNode(center, ImGuiDir_Left, 0.15f, nullptr, &center);
        const ImGuiID right = ImGui::DockBuilderSplitNode(center, ImGuiDir_Right, 0.30f, nullptr, &center);
        const ImGuiID bottom = ImGui::DockBuilderSplitNode(center, ImGuiDir_Down, 0.44f, nullptr, &center);
        ImGui::DockBuilderDockWindow("Case Explorer", left);
        ImGui::DockBuilderDockWindow("Inspector", right);
        ImGui::DockBuilderDockWindow("3D Viewer", center);
        ImGui::DockBuilderDockWindow("Log", bottom);
        ImGui::DockBuilderDockWindow("Run Monitor", bottom);
        ImGui::DockBuilderDockWindow("Plots", bottom);
        ImGui::DockBuilderDockWindow("Studies", bottom);
        ImGui::DockBuilderFinish(id);
    }
    ImGui::DockSpaceOverViewport(id, vp, ImGuiDockNodeFlags_None);
}

void App::menuBar() {
    if (!ImGui::BeginMainMenuBar()) return;
    if (ImGui::BeginMenu(TR("File"))) {
        if (ImGui::BeginMenu(TR(ICON_FILE_PLUS "  New"))) {
            if (ImGui::MenuItem(TR("Demo: 3-phase waterflood (SPE1-like PVT)"))) { if (maybeSaveChanges()) newCase(0); }
            if (ImGui::MenuItem(TR("Validation: 1D Buckley-Leverett"))) { if (maybeSaveChanges()) newCase(1); }
            if (ImGui::MenuItem(TR(ICON_WAND_SPARKLES "  Model wizard..."))) showWizard_ = true;   // C40
            if (ImGui::MenuItem(TR("Blank model"), "Ctrl+N")) { if (maybeSaveChanges()) newCase(2); }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu(TR(ICON_BOOK_OPEN "  Examples"))) { examplesMenu(false); ImGui::EndMenu(); }
        if (ImGui::BeginMenu(TR(ICON_BOOK_OPEN "  Example decks (examples, decks folders)"))) { deckFilesMenu(); ImGui::EndMenu(); }
        if (ImGui::MenuItem(TR(ICON_FOLDER_OPEN "  Open..."), "Ctrl+O")) openDialog();
        if (ImGui::BeginMenu(TR(ICON_HISTORY "  Recent"))) { recentMenu(); ImGui::EndMenu(); }
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_SAVE "  Save"), "Ctrl+S")) save(false);
        if (ImGui::MenuItem(TR(ICON_SAVE "  Save As..."), "Ctrl+Shift+S")) save(true);
        if (ImGui::MenuItem(TR(ICON_FILE_CODE "  Edit deck file in Deck Editor..."))) {
            std::string p;
            if (fileDialog(false, "ECLIPSE deck / include (*.DATA;*.INC;*.GRDECL)\0*.DATA;*.data;*.INC;*.inc;*.GRDECL;*.grdecl;*.SCH;*.sch;*.PROP;*.prop\0All files (*.*)\0*.*\0", "DATA", p)) {
                edOpenFile(p);
                workspace_ = Workspace::DeckEditor;
            }
        }
        if (ImGui::MenuItem(TR("Open current model as deck text (Deck Editor)"))) { edOpenModelDeck(); workspace_ = Workspace::DeckEditor; }
        if (ImGui::MenuItem(TR(ICON_UPLOAD "  Import ECLIPSE deck (.DATA)..."))) {
            std::string p;
            if (maybeSaveChanges() && fileDialog(false, "ECLIPSE deck (*.DATA;*.data)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "DATA", p)) { dirty_ = false; importDeck(p); }
        }
        if (ImGui::MenuItem(TR(ICON_DOWNLOAD "  Export ECLIPSE deck (.DATA)..."))) exportDeck();
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_FOLDER_OPEN "  Open ResSim results (.rsres)..."))) {
            std::string p;
            if (fileDialog(false, "ResSim results (*.rsres)\0*.rsres\0All files (*.*)\0*.*\0", "rsres", p)) openResultsFile(p);
        }
        {
            bool has = false;
            { std::lock_guard<std::mutex> l(runner_.mutex()); has = runner_.results().valid && !runner_.results().summary.empty(); }
            if (ImGui::MenuItem(TR(ICON_SAVE "  Save results (.rsres)..."), nullptr, false, has && !runner_.running())) saveResultsDialog();
            if (ImGui::MenuItem(TR(ICON_DOWNLOAD "  Export results as ECLIPSE files (EGRID / UNRST / UNSMRY)..."), nullptr, false, has && !runner_.running()))
                exportResultsEclipseDialog();
            if (ImGui::MenuItem(TR(ICON_FILE_TEXT "  Generate report (Word / PowerPoint / HTML)..."), nullptr, showReport_, has)) showReport_ = true;
        }
        if (ImGui::MenuItem(TR(ICON_CHART_LINE "  Open results (ECLIPSE / OPM Flow)..."))) {
            std::string p;
            if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.FEGRID;*.UNRST;*.FUNRST;*.INIT;*.FINIT;*.UNSMRY;*.FUNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p))
                rvOpen(p);
        }
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_LOG_OUT "  Exit"), "Alt+F4")) quit_ = true;
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu(TR("Edit"))) {   // C39
        char u[64], r[64];
        std::snprintf(u, sizeof(u), ICON_UNDO_2 "  Undo%s", canUndo() ? "" : " (nothing to undo)");
        std::snprintf(r, sizeof(r), ICON_REDO_2 "  Redo%s", canRedo() ? "" : " (nothing to redo)");
        if (ImGui::MenuItem(u, "Ctrl+Z", false, canUndo())) undoApply(-1);
        if (ImGui::MenuItem(r, "Ctrl+Y", false, canRedo())) undoApply(+1);
        ImGui::Separator();
        ImGui::TextDisabled("History: %d edit step(s) kept (up to %d for this model size)", std::max(0, int(undo_.stack.size()) - 1), int(undo_.maxLevels) - 1);
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_SEARCH "  Command palette"), "Ctrl+Shift+P")) { paletteOpen_ = true; paletteBuf_[0] = 0; paletteSel_ = 0; paletteFocus_ = true; }
        if (ImGui::MenuItem(TR(ICON_SETTINGS "  Settings..."))) showOpmSettings_ = true;
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu(TR("Run"))) {
        if (ImGui::MenuItem(TR(ICON_PLAY "  Run simulation"), "F5", false, !runner_.running())) runSimulation();
        if (ImGui::MenuItem(TR(ICON_SQUARE "  Stop"), "Shift+F5", false, runner_.running())) runner_.cancel();
        if (ImGui::MenuItem(TR(ICON_PLAY "  Restart from a report step..."), nullptr, showRestart_, !runner_.running())) showRestart_ = true;
        if (ImGui::MenuItem(TR(ICON_CHECK "  Deck regression (OPM / ECLIPSE references)..."), nullptr, showRegress_)) showRegress_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem(TR("Engine: ResSim (built-in)"), nullptr, settings_.engine == 0)) { settings_.engine = 0; saveSettings(); }
        if (ImGui::MenuItem(TR("Engine: OPM Flow (external)"), nullptr, settings_.engine == 1)) { settings_.engine = 1; saveSettings(); detectFlow(false); }
        if (ImGui::MenuItem(TR("Engine: ECLIPSE (external, licensed)"), nullptr, settings_.engine == 2)) { settings_.engine = 2; saveSettings(); detectEclipse(false); }
        if (ImGui::MenuItem(TR(ICON_SETTINGS "  Simulation engine settings (OPM Flow / ECLIPSE)..."))) showOpmSettings_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_ARROW_DOWN_TO_LINE "  Compute initial state"))) { computeInit(); section_ = Section::Init; }
        if (ImGui::MenuItem(TR(ICON_CHECK "  Validate inputs"))) {
            const std::string& v = validation();
            notify(v.empty() ? "All inputs are valid" : v);
        }
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu(TR("View"))) {
        if (ImGui::MenuItem(TR(ICON_BOX "  Model workspace"), "Ctrl+1", workspace_ == Workspace::Model)) workspace_ = Workspace::Model;
        if (ImGui::MenuItem(TR(ICON_FILE_CODE "  Deck Editor (writing mode)"), "Ctrl+2", workspace_ == Workspace::DeckEditor)) workspace_ = Workspace::DeckEditor;
        if (ImGui::MenuItem(TR(ICON_CHART_LINE "  Results workspace"), "Ctrl+3", workspace_ == Workspace::Results)) workspace_ = Workspace::Results;
        if (ImGui::MenuItem(TR(ICON_FLASK_CONICAL "  Studies workspace (sensitivity / HM)"), "Ctrl+4", workspace_ == Workspace::Studies)) workspace_ = Workspace::Studies;
        if (ImGui::MenuItem(TR(ICON_MOUNTAIN "  Geomodel workspace (ResGeo)"), "Ctrl+5", workspace_ == Workspace::Geomodel)) workspace_ = Workspace::Geomodel;
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_LAYOUT_DASHBOARD "  Reset window layout"))) resetLayout_ = true;
        if (ImGui::MenuItem(TR(ICON_MAXIMIZE "  Reset 3D camera")) && viewOk_) view_.resetCamera();
        if (ImGui::BeginMenu(TR(ICON_ZOOM_IN "  Interface size"))) {
            static const float sizes[] = {0.70f, 0.80f, 0.90f, 1.00f, 1.15f, 1.30f, 1.50f};
            static const char* names[] = {"70 %  (very compact)", "80 %  (compact)", "90 %", "100 %  (default)", "115 %", "130 %", "150 %"};
            for (int i = 0; i < 7; ++i)
                if (ImGui::MenuItem(names[i], nullptr, std::fabs(settings_.uiScale - sizes[i]) < 0.01f)) { settings_.uiScale = sizes[i]; setUiScale(sizes[i]); saveSettings(); }
            ImGui::Separator();
            ImGui::TextDisabled("%s", TR("Ctrl + / Ctrl - / Ctrl 0"));
            ImGui::EndMenu();
        }
        if (ImGui::MenuItem(TR(ICON_LAYERS "  Case manager (variants)"), nullptr, showCases_)) showCases_ = !showCases_;
        if (ImGui::MenuItem(TR(ICON_LIST_CHECKS "  Workflow editor (record / replay)"), nullptr, showWorkflow_)) showWorkflow_ = !showWorkflow_;
        if (ImGui::BeginMenu(TR(ICON_PALETTE "  Appearance"))) {   // C40: theme / accent / 3-D background without opening Settings
            ImGui::TextDisabled("%s", TR("Theme"));
            for (int t = 0; t < 3; ++t)
                if (ImGui::MenuItem(themeNames()[t], nullptr, settings_.theme == t)) { settings_.theme = t; applyColors(t, settings_.accent); saveSettings(); }
            ImGui::Separator(); ImGui::TextDisabled("%s", TR("Accent colour"));
            for (int a = 0; a < 5; ++a) {
                const ImVec4 c = accentColor(a);
                ImGui::PushStyleColor(ImGuiCol_Text, c); ImGui::Bullet(); ImGui::PopStyleColor(); ImGui::SameLine();
                if (ImGui::MenuItem(accentNames()[a], nullptr, settings_.accent == a)) { settings_.accent = a; applyColors(settings_.theme, a); saveSettings(); }
            }
            ImGui::Separator(); ImGui::TextDisabled("%s", TR("3-D background (Model / Results)"));
            for (int b = 0; b < 5; ++b)
                if (ImGui::MenuItem(background3dNames()[b], nullptr, settings_.bg3d == b)) { settings_.bg3d = b; set3dBackground(b); saveSettings(); }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu(TR("Language"))) {   // C41 (G13)
            static std::vector<i18n::Language> langs;
            if (ImGui::IsWindowAppearing() || langs.empty()) langs = i18n::available();
            for (const auto& l : langs)
                if (ImGui::MenuItem(l.name.c_str(), l.code.c_str(), settings_.language == l.code)) {
                    settings_.language = l.code; i18n::setLanguage(l.code); saveSettings();
                }
            ImGui::Separator();
            ImGui::TextDisabled("%s", TR("Catalogues: lang/<code>.txt next to the program."));
            ImGui::TextDisabled("%s", TR("Partial translation: menus, toolbar, workspaces."));
            ImGui::TextDisabled("%s", TR("Scripts that need shaping (Bengali, Arabic, ...) are not supported yet."));
            ImGui::EndMenu();
        }
        if (ImGui::MenuItem(TR(ICON_SETTINGS "  Settings..."))) showOpmSettings_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem(TR("Result Viewer (ECLIPSE / OPM output, compare)"), nullptr, showResultViewer_)) showResultViewer_ = !showResultViewer_;
        if (ImGui::MenuItem(TR("Deck import report"), nullptr, showDeckTools_ && deckTab_ == 0)) { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; }
        if (ImGui::MenuItem(TR("ECLIPSE keyword reference"), nullptr, showDeckTools_ && deckTab_ == 1)) { showDeckTools_ = true; deckTab_ = 1; deckTabSet_ = true; }
        ImGui::Separator();
        ImGui::MenuItem(TR("UI demo windows (developer)"), nullptr, &showImGuiDemo_);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu(TR("Tools"))) {
        if (ImGui::MenuItem(TR(ICON_FLASK_CONICAL "  Studies workspace: sensitivity & history matching"), "Ctrl+4")) workspace_ = Workspace::Studies;
        if (ImGui::MenuItem(TR(ICON_PANEL_BOTTOM "  Studies panel (Model workspace)"), nullptr, showStudies_)) { showStudies_ = true; ImGui::SetWindowFocus("Studies"); }
        if (ImGui::MenuItem(TR(ICON_ROUTE "  Streamlines in 3D viewer"), nullptr, sl_.show)) sl_.show = !sl_.show;
        if (ImGui::MenuItem(TR(ICON_ACTIVITY "  Geomechanics: fault & fracture stability"), nullptr, showStability_)) { showStability_ = !showStability_; stabDock_ = showStability_; }
        if (ImGui::MenuItem(TR(ICON_ACTIVITY "  Geomechanics: seismicity-aware injection design"), nullptr, showSeis_)) { showSeis_ = !showSeis_; seisDock_ = showSeis_; }
        if (ImGui::MenuItem(TR(ICON_ACTIVITY "  Geomechanics: fracture contact & slip (2-D)"), nullptr, showEfem_)) { showEfem_ = !showEfem_; efemDock_ = showEfem_; }
        ImGui::Separator();
        if (ImGui::MenuItem(TR("Convert grid to corner-point (COORD/ZCORN)"), nullptr, false, !case_.grid.cornerPoint())) {
            convertToCornerPoint(case_.grid); markDirty(true); notify("Grid converted to corner-point geometry");
        }
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu(TR("Help"))) {
        if (ImGui::MenuItem(TR(ICON_BOOK_OPEN "  User guide (all features & examples)"), "F1")) openUserGuide();
        if (ImGui::MenuItem(TR(ICON_SEARCH "  Command palette"), "Ctrl+Shift+P")) { paletteOpen_ = true; paletteBuf_[0] = 0; paletteSel_ = 0; paletteFocus_ = true; }
        if (ImGui::MenuItem(TR(ICON_HOUSE "  Start page"))) showStart_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem(TR(ICON_INFO "  About ResSim Studio"))) showAbout_ = true;
        ImGui::EndMenu();
    }
    // right-aligned case title
    const std::string t = case_.title + (dirty_ ? "  \xE2\x97\x8F" : "");
    const float w = ImGui::CalcTextSize(t.c_str()).x;
    ImGui::SameLine(ImGui::GetWindowWidth() - w - ImGui::GetStyle().ItemSpacing.x * 3);
    ImGui::TextDisabled("%s", t.c_str());
    ImGui::EndMainMenuBar();
}

void App::toolbar() {
    const float h = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y * 1.2f;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(ImGui::GetStyle().WindowPadding.x, 4));
    if (ImGui::BeginViewportSideBar("##toolbar", ImGui::GetMainViewport(), ImGuiDir_Up, h,
                                    ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings)) {
        auto tip = [](const char* t) { if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", t); };
        if (ImGui::Button(TR(ICON_FILE_PLUS " New"))) ImGui::OpenPopup("newpop");
        tip(TR("New case"));
        if (ImGui::BeginPopup("newpop")) {
            if (ImGui::Selectable(TR("Demo: 3-phase waterflood"))) { if (maybeSaveChanges()) newCase(0); }
            if (ImGui::Selectable(TR("Validation: 1D Buckley-Leverett"))) { if (maybeSaveChanges()) newCase(1); }
            if (ImGui::Selectable(TR(ICON_WAND_SPARKLES "  Model wizard..."))) showWizard_ = true;   // C40
            if (ImGui::Selectable(TR("Blank model"))) { if (maybeSaveChanges()) newCase(2); }
            ImGui::Separator();
            if (ImGui::BeginMenu(TR("Examples"))) { examplesMenu(true); ImGui::EndMenu(); }
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_FOLDER_OPEN " Open"))) openDialog();
        tip(TR("Open project (Ctrl+O)"));
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_SAVE " Save"))) save(false);
        tip(TR("Save project (Ctrl+S)"));
        ImGui::SameLine();
        ImGui::BeginDisabled(!canUndo());
        if (ImGui::Button(TR(ICON_UNDO_2 "##undo"))) undoApply(-1);
        ImGui::EndDisabled();
        tip(TR("Undo the last model edit (Ctrl+Z)"));
        ImGui::SameLine(0, 2);
        ImGui::BeginDisabled(!canRedo());
        if (ImGui::Button(TR(ICON_REDO_2 "##redo"))) undoApply(+1);
        ImGui::EndDisabled();
        tip(TR("Redo (Ctrl+Y)"));
        ImGui::SameLine();
        ImGui::TextDisabled("%s", TR("|"));
        ImGui::SameLine();
        workspaceSwitch();
        ImGui::SameLine();
        ImGui::TextDisabled("%s", TR("|"));
        ImGui::SameLine();
        const bool running = runner_.running();
        ImGui::BeginDisabled(running);
        if (ui::accentButton(ICON_PLAY "  Run")) runSimulation();
        ImGui::EndDisabled();
        tip(TR("Run the full simulation (F5)"));
        ImGui::SameLine();
        ImGui::BeginDisabled(!running);
        if (ui::dangerButton(ICON_SQUARE "  Stop")) runner_.cancel();
        ImGui::EndDisabled();
        tip(TR("Stop the running simulation (Shift+F5)"));
        ImGui::SameLine();
        {
            static const char* engines[] = {"ResSim", "OPM Flow", "ECLIPSE"};
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6.5f);
            ImGui::BeginDisabled(running);
            if (ImGui::Combo("##engine", &settings_.engine, engines, 3)) {
                saveSettings();
                if (settings_.engine == 1) detectFlow(false);
                if (settings_.engine == 2) detectEclipse(false);
            }
            ImGui::EndDisabled();
            tip(TR("Simulation engine used by Run: ResSim (built-in), OPM Flow (external flow.exe) or ECLIPSE (licensed eclrun / eclipse.exe). Settings: Run menu."));
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_SETTINGS "##eng"))) showOpmSettings_ = true;
            tip(TR("Simulation engine settings (OPM Flow / ECLIPSE executable, auto-detect, options)"));
            ImGui::SameLine();
        }
        if (ImGui::Button(TR(ICON_ARROW_DOWN_TO_LINE " Initialise"))) { computeInit(); section_ = Section::Init; }
        tip(TR("Compute the equilibrium initial state and fluids in place"));
        ImGui::SameLine();
        if (ImGui::Button(TR(ICON_CHECK " Validate"))) {
            const std::string& v = validation();
            notify(v.empty() ? "All inputs are valid" : v);
        }
        tip(TR("Check all inputs for consistency"));
        ImGui::SameLine();
        // right side: validation indicator
        const std::string& v = validation();
        const bool ext = !case_.engineUnsupported.empty();
        const char* st = ext ? TR(ICON_INFO " OPM Flow / ECLIPSE deck") : !v.empty() ? TR(ICON_TRIANGLE_ALERT " Input problem")
                         : case_.thermal.enabled ? TR(ICON_CHECK " Thermal model valid")
                         : case_.comp.enabled ? (case_.comp.co2store ? TR(ICON_CHECK " CO2 storage model valid") : TR(ICON_CHECK " Compositional model valid"))
                         : TR(ICON_CHECK " Inputs valid");
        const float w = ImGui::CalcTextSize(st).x;
        // the progress bar shares the row with the validation indicator: shrink it (or drop the indicator) on narrow windows
        const float rowEnd = ImGui::GetWindowWidth() - ImGui::GetStyle().WindowPadding.x * 2;
        bool showStatus = true;
        if (running) {
            ImGui::SameLine();
            char b[64]; std::snprintf(b, sizeof(b), "%.0f%%  (%.1f s)", runner_.progress() * 100.0, runner_.wallSeconds());
            const float fsz = ImGui::GetFontSize(), room = rowEnd - ImGui::GetCursorPosX() - w - fsz;
            showStatus = room >= fsz * 6;
            const float bw = std::max(fsz * 6, std::min(fsz * 14, showStatus ? room : rowEnd - ImGui::GetCursorPosX()));
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col4(Pal::Accent));
            ImGui::ProgressBar(float(runner_.progress()), ImVec2(bw, 0), b);
            ImGui::PopStyleColor();
        }
        if (showStatus) {
        ImGui::SameLine(ImGui::GetWindowWidth() - w - ImGui::GetStyle().WindowPadding.x * 2);
        ImGui::TextColored(col4(ext ? Pal::Accent : v.empty() ? Pal::Ok : Pal::Warn), "%s", st);
        if (ext && ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted("This deck uses features the ResSim engine cannot model. It is shown with placeholders\n"
                                   "where needed and runs unchanged (original file) with OPM Flow or ECLIPSE:");
            for (auto& r : case_.engineUnsupported) ImGui::BulletText("%s", r.c_str());
            ImGui::EndTooltip();
        } else if (!v.empty() && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", v.c_str());
        }
    }
    ImGui::End();
    ImGui::PopStyleVar();
}

void App::statusBar() {
    const float h = ImGui::GetFrameHeight();
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(ImGui::GetStyle().WindowPadding.x, 2));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, themeStatusBg());
    if (ImGui::BeginViewportSideBar("##status", ImGui::GetMainViewport(), ImGuiDir_Down, h, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings)) {
        if (runner_.running()) ImGui::TextColored(col4(Pal::Accent), ICON_ACTIVITY " Running  %.0f%%", runner_.progress() * 100.0);
        else ImGui::TextColored(col4(Pal::Ok), "%s", TR(ICON_CIRCLE_DOT " Ready"));
        if (study_.running) { ImGui::SameLine(0, 16); ImGui::TextColored(col4(Pal::Accent), ICON_FLASK_CONICAL " study %d / %d", int(study_.done), int(study_.total)); }
        ImGui::SameLine(0, 24);
        ImGui::TextDisabled(ICON_FILE " %s%s", path_.empty() ? "(not saved)" : path_.c_str(), dirty_ ? "  - modified" : "");
        ImGui::SameLine(0, 24);
        ImGui::TextDisabled("%s units", case_.units == UnitSystem::Field ? "FIELD" : "METRIC");
        ImGui::SameLine(0, 24);
        ImGui::TextDisabled(ICON_SERVER " %s", engineName(settings_.engine));
        ImGui::SameLine(0, 24);
        if (preview_) ImGui::TextDisabled(ICON_GRID_3X3 " %d active cells", preview_->numActive());
        else ImGui::TextColored(col4(Pal::Err), "grid error");
        ImGui::SameLine(0, 24);
        notificationButton();
        char r[96];
        const int thr = case_.tuning.threads > 0 ? case_.tuning.threads : ressim::physicalCores();
        std::snprintf(r, sizeof(r), ICON_CPU " %d CPU threads  |  %.0f fps", thr, ImGui::GetIO().Framerate);
        ImGui::SameLine(ImGui::GetWindowWidth() - ImGui::CalcTextSize(r).x - ImGui::GetStyle().WindowPadding.x * 2);
        ImGui::TextDisabled("%s", r);
    }
    ImGui::End();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();
}

// ============================================================================ panels
void App::panelCase() {
    if (!ImGui::Begin("Case Explorer")) { ImGui::End(); return; }
    if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
    ImGui::TextWrapped("%s", case_.title.c_str());
    if (g_fonts.bold) ImGui::PopFont();
    ImGui::TextDisabled("%s", path_.empty() ? "unsaved project" : baseName(path_).c_str());
    ImGui::Spacing();
    if (ImGui::CollapsingHeader(ICON_LIST_CHECKS "  Workflow", ImGuiTreeNodeFlags_DefaultOpen)) workflowChecklist();
    if (ImGui::CollapsingHeader(ICON_FOLDER_TREE "  Objects", ImGuiTreeNodeFlags_DefaultOpen)) objectTree();   // C40
    ui::header("Model inputs");
    for (int s = 0; s < int(Section::Count); ++s) {
        const std::string st = sectionStatus(Section(s));
        ImGui::PushID(s);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        char lb[96]; std::snprintf(lb, sizeof(lb), "   %s  %s", kSectionIcons[s], TR(kSectionNames[s]));
        if (ImGui::Selectable(lb, section_ == Section(s), 0, ImVec2(0, ImGui::GetFrameHeight()))) section_ = Section(s);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s%s%s %s", TR(kSectionHelp[s]), st.empty() ? "" : "\n\n", st.empty() ? "" : TR("Problem:"), st.c_str());
        {   // C40 (G12): status by icon and colour (not colour alone)
            const float fs = ImGui::GetFontSize() * 0.8f;
            ImGui::GetWindowDrawList()->AddText(ImGui::GetFont(), fs, ImVec2(p.x + 1, p.y + (ImGui::GetFrameHeight() - fs) * 0.5f), st.empty() ? Pal::Ok : Pal::Err,
                                                st.empty() ? ICON_CHECK : ICON_TRIANGLE_ALERT);
        }
        ImGui::PopID();
    }
    ImGui::Spacing();
    ui::header("Results");
    {
        int nsnap = 0, nsum = 0; bool valid = false;
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            nsnap = int(runner_.results().snaps.size()); nsum = int(runner_.results().summary.size()); valid = runner_.results().valid;
        }
        if (runner_.running()) ImGui::TextColored(col4(Pal::Accent), "Running... %.0f%%", runner_.progress() * 100.0);
        else if (valid) {
            ImGui::TextColored(col4(Pal::Ok), "\xE2\x9C\x93 %d time steps, %d reports", nsum, std::max(0, nsnap - 1));
            if (resultsCaseVersion_ != caseVersion_) ImGui::TextColored(col4(Pal::Warn), "Inputs changed since the run");
        } else ImGui::TextDisabled("No results yet");
        if (init_.valid) ImGui::TextColored(col4(Pal::Ok), "\xE2\x9C\x93 Initial state computed");
    }
    ImGui::Spacing();
    ui::header("Summary");
    if (ui::beginProps("summary", 0.42f)) {
        char b[96];
        std::snprintf(b, sizeof(b), "%d\xC3\x97%d\xC3\x97%d", case_.grid.nx, case_.grid.ny, case_.grid.nz); ui::propText("Grid", b);
        std::string ph;
        if (case_.phases.water) ph += "W ";
        if (case_.phases.oil) ph += "O ";
        if (case_.phases.gas) ph += case_.phases.disgas ? (case_.phases.vapoil ? "G (live, wet)" : "G (live)") : (case_.phases.vapoil ? "G (wet)" : "G");
        ui::propText("Phases", ph.c_str());
        std::snprintf(b, sizeof(b), "%d", int(case_.wells.size())); ui::propText("Wells", b);
        double T = 0; for (auto& s : case_.schedule) T += s.length;
        std::snprintf(b, sizeof(b), "%.1f years", T / 365.25); ui::propText("Duration", b);
        if (init_.valid && case_.phases.oil) {
            const Units U = Units::of(case_.units);
            std::snprintf(b, sizeof(b), "%.3f MM %s", init_.ip.oil / 1e6, U.liq); ui::propText("STOIIP", b);
        }
        ui::endProps();
    }
    ImGui::End();
}

void App::panelEditor() {
    if (!ImGui::Begin("Inspector")) { ImGui::End(); return; }
    const int s = int(section_);
    if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
    ImGui::Text("%s  %s", kSectionIcons[s], TR(kSectionNames[s]));
    if (g_fonts.title) ImGui::PopFont();
    ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
    ImGui::TextWrapped("%s", TR(kSectionHelp[s]));   // C48: translated section help
    ImGui::PopStyleColor();
    const std::string st = sectionStatus(section_);
    if (!st.empty()) {
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.35f, 0.12f, 0.12f, 0.45f));
        ImGui::BeginChild("prob", ImVec2(0, 0), ImGuiChildFlags_AutoResizeY | ImGuiChildFlags_AlwaysUseWindowPadding);
        ImGui::PushTextWrapPos(0);
        ImGui::TextColored(col4(Pal::Err), "\xE2\x9A\xA0 %s", st.c_str());
        ImGui::PopTextWrapPos();
        ImGui::EndChild();
        ImGui::PopStyleColor();
    }
    // C40 (G11): filter the property rows of the section by label (Ctrl+F in the Inspector, Esc clears)
    std::string& flt = inspFilter_;
    {
        ImGui::SetNextItemWidth(-ImGui::GetFrameHeight() * 2 - ImGui::GetStyle().ItemSpacing.x * 2);
        if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_F, false)) ImGui::SetKeyboardFocusHere();
        ui::inputStringHint("##inspFilter", ICON_SEARCH "  Filter properties (Ctrl+F)", &flt);
        if (ImGui::IsItemActive() && ImGui::IsKeyPressed(ImGuiKey_Escape, false)) flt.clear();
        ImGui::SameLine();
        ImGui::BeginDisabled(flt.empty());
        if (ImGui::Button(ICON_X "##fltclr", ImVec2(ImGui::GetFrameHeight(), 0))) flt.clear();
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Clear the filter");
        // C41: show only the properties modified since the baseline (project open / new / import, or "set baseline")
        ui::ModTrack& mt = ui::modTrack();
        ImGui::SameLine();
        const bool on = mt.onlyModified;
        if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
        if (ImGui::Button(ICON_PENCIL "##onlymod", ImVec2(ImGui::GetFrameHeight(), 0))) mt.onlyModified = !mt.onlyModified;
        if (on) ImGui::PopStyleColor();
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s\nRight-click: set the current values as the new baseline",
                                                      on ? "Showing only modified properties (click: show all)" : "Show only properties modified since the project was opened");
        if (ImGui::IsItemClicked(ImGuiMouseButton_Right)) { ui::modBaselineReset(); notify("Inspector baseline set to the current values"); }
    }
    ImGui::Separator();
    ImGui::BeginChild("editor", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_None);
    ui::PropFilter& pf = ui::propFilter();
    pf.active = !flt.empty(); pf.text = flt; pf.shown = pf.hidden = 0;
    ui::ModTrack& mt = ui::modTrack();
    mt.inInspector = true; mt.modified = mt.hidden = 0;
    ImGui::PushID(s);
    switch (section_) {
        case Section::Runspec: editRunspec(); break;
        case Section::Grid: editGrid(); break;
        case Section::PVT: editPVT(); break;
        case Section::SCAL: editSCAL(); break;
        case Section::Init: editInit(); break;
        case Section::Aquifers: editAquifers(); break;
        case Section::EOR: editEOR(); break;
        case Section::Wells: editWells(); break;
        case Section::Schedule: editSchedule(); break;
        case Section::Numerics: editNumerics(); break;
        case Section::Advanced: editAdvanced(); break;
        default: break;
    }
    ImGui::PopID();
    if (pf.active) {
        ImGui::Spacing();
        if (pf.shown == 0) ImGui::TextDisabled(ICON_FILTER "  No property of the open groups matches \"%s\" (closed groups and other tabs are not searched)", flt.c_str());
        else ImGui::TextDisabled(ICON_FILTER "  %d matching, %d hidden by the filter", pf.shown, pf.hidden);
    }
    mt.inInspector = false;
    if (mt.onlyModified) {
        ImGui::Spacing();
        if (mt.modified == 0) ImGui::TextDisabled(ICON_PENCIL "  No modified property in the open groups of this section (closed groups and other tabs are not searched)");
        else ImGui::TextDisabled(ICON_PENCIL "  %d modified, %d unmodified hidden", mt.modified, mt.hidden);
    }
    pf.active = false;
    ImGui::EndChild();
    ImGui::End();
}

void App::panelLog() {
    if (!ImGui::Begin("Log")) { ImGui::End(); return; }
    bool clear = ImGui::SmallButton("Clear");
    ImGui::SameLine();
    const bool copy = ImGui::SmallButton("Copy all");
    ImGui::SameLine();
    ImGui::Checkbox("Auto-scroll", &logAutoScroll_);
    ImGui::SameLine();
    static ImGuiTextFilter filter;
    filter.Draw("Filter", ImGui::GetFontSize() * 12);
    ImGui::Separator();
    ImGui::BeginChild("logtext", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_HorizontalScrollbar);
    if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        auto& lines = runner_.log();
        if (clear) lines.clear();
        if (copy) {
            std::string all;
            for (auto& s : lines) { all += s; all += '\n'; }
            ImGui::SetClipboardText(all.c_str());
        }
        auto colorOf = [](const std::string& s) -> ImU32 {
            auto has = [&](const char* k) { return s.find(k) != std::string::npos; };
            if (has("ERROR") || has("error") || has("fail") || has("Fail")) return Pal::Err;
            if (has("===")) return Pal::Accent;
            if (has("WARN") || has("warn") || has("cut") || has("switch")) return Pal::Warn;
            return 0;
        };
        auto show = [&](const std::string& s) {
            const ImU32 c = colorOf(s);
            if (c) ImGui::PushStyleColor(ImGuiCol_Text, col4(c));
            ImGui::TextUnformatted(s.c_str());
            if (c) ImGui::PopStyleColor();
        };
        if (filter.IsActive()) {
            for (auto& s : lines) if (filter.PassFilter(s.c_str())) show(s);
        } else {
            ImGuiListClipper clip;
            clip.Begin(int(lines.size()));
            while (clip.Step())
                for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) show(lines[i]);
        }
    }
    if (g_fonts.mono) ImGui::PopFont();
    if (logAutoScroll_ && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4) ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
    ImGui::End();
}

void App::aboutWindow() {
    if (!showAbout_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 34, 0), ImGuiCond_Appearing);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (ImGui::Begin("About ResSim Studio", &showAbout_, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_AlwaysAutoResize)) {
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::TextUnformatted("ResSim Studio");
        if (g_fonts.title) ImGui::PopFont();
        ImGui::TextDisabled("Black-oil reservoir simulation studio");
        ImGui::Separator();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 32);
        ImGui::TextUnformatted(
            "Fully-implicit three-phase black-oil simulator with dissolved gas (variable switching), "
            "fully coupled wells (rate/BHP control with automatic switching, WCONHIST), hydrostatic equilibrium, "
            "rock and water compressibility, capillary pressure and gravity.\n\n"
            "Grids: Cartesian or corner-point (COORD/ZCORN) with fault throws, MULTFLT and NNCs; Fetkovich and "
            "Carter-Tracy aquifers. ECLIPSE deck import and export (two-way).\n\n"
            "Numerics: Newton-Raphson with forward-mode automatic differentiation, CPR (AMG pressure) or ILU(0) "
            "preconditioned BiCGSTAB/GMRES, OpenMP-parallel assembly and sparse kernels.\n\n"
            "Studies: parallel batch runner, tornado and Latin-hypercube Monte Carlo sensitivity, manual history "
            "matching, Levenberg-Marquardt (finite-difference or exact forward-sensitivity Jacobian), discrete-adjoint "
            "L-BFGS and ES-MDA assisted history matching. Multi-component product analysis (H2, syngas, heating "
            "values) for compositional, CO2 and thermal / UCG runs. Unified multi-component framework: Henry-law "
            "solubility (IAPWS G7-04), Peng-Robinson gas and convex three-phase flash in the thermal core, with the same "
            "CO2 trapping inventory (structural, residual, dissolved, mineral) as the compositional / CO2STORE cores and "
            "per-component 3D maps over time. Pollock streamlines with "
            "time-of-flight, drainage regions and well allocation.\n\n"
            "Graphics: Direct3D 11 (3D grid viewer with MSAA, picking, filters and colour maps), "
            "Direct2D/DirectWrite (viewport text and legends), Dear ImGui docking UI and ImPlot charts.");
        ImGui::PopTextWrapPos();
        ImGui::Separator();
        ImGui::TextDisabled("Roadmap: GPU linear solver.");
    }
    ImGui::End();
}

// ============================================================================ actions
void App::notify(const std::string& msg) {
    pushNote(msg);
    toast_ = msg;
    toastUntil_ = ImGui::GetTime() + 3.5;
    runner_.addLog("[ui] " + msg);
}

void App::markDirty(bool gridChanged) {
    dirty_ = true;
    ++caseVersion_;
    // rebuild immediately so every panel later in this frame sees a grid consistent with case_.grid
    if (gridChanged) rebuildPreview();
}

const std::string& App::validation() const {
    if (validatedVersion_ != caseVersion_) {
        validation_ = validateCase(case_);
        validatedVersion_ = caseVersion_;
    }
    return validation_;
}

std::string App::sectionStatus(Section s) const {
    if (s == Section::Grid && !gridError_.empty()) return gridError_;
    const std::string& v = validation();
    if (v.empty()) return {};
    auto starts = [&](const char* p) { return v.rfind(p, 0) == 0; };
    auto has = [&](const char* k) { return v.find(k) != std::string::npos; };
    Section owner = Section::Runspec;
    if (starts("GRID")) owner = Section::Grid;
    else if (starts("Aquifer") || has("AQU")) owner = Section::Aquifers;
    else if (starts("POLYMER") || has("PLY") || has("WPOLYMER") || has("PLMIX")) owner = Section::EOR;
    else if (starts("RUNSPEC")) owner = Section::Runspec;
    else if (starts("PROPS")) owner = (has("SWOF") || has("SGOF")) ? Section::SCAL : Section::PVT;
    else if (starts("TUNING")) owner = Section::Numerics;
    else if (starts("SCHEDULE")) owner = (has("well name") || has("no name")) ? Section::Wells : Section::Schedule;
    else if (starts("Well ")) owner = (has("control") || has("target") || has("RATE")) ? Section::Schedule : Section::Wells;
    return owner == s ? v : std::string();
}

void App::rebuildPreview() {
    gridDirty_ = false;
    try {
        auto g = std::make_shared<Grid>();
        g->build(case_);
        preview_ = g;
        gridError_.clear();
    } catch (const std::exception& e) {
        preview_.reset();
        gridError_ = e.what();
    }
}

void App::computeInit() {
    init_ = InitPreview{};
    try {
        const std::string v = validation();
        if (!v.empty()) throw std::runtime_error(v);
        if (case_.thermal.enabled) {   // thermal: initial state of the thermal core (p, T, saturations)
            thermal::DeckModel m; std::string err;
            if (!thermal::deckModel(case_, m, &err)) throw std::runtime_error(err);
            const double pu = case_.units == UnitSystem::Field ? 6894.757293168 : 1e5;
            bool got = false;
            std::vector<thermal::Summary> out;
            thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator& sim) {
                const size_t n = sim.pressure().size();
                init_.p.resize(n); init_.sw.resize(n); init_.sg.resize(n); init_.so.resize(n); init_.rs.assign(n, 0.0f);
                for (size_t q = 0; q < n; ++q) {
                    init_.p[q] = float(sim.pressure()[q] / pu); init_.sw[q] = float(sim.sw()[q]);
                    init_.sg[q] = float(sim.sg()[q]); init_.so[q] = float(sim.so()[q]);
                }
                init_.ip.oil = S.foip;
                got = true;
                return false;
            }, &err, 0);
            if (!got) throw std::runtime_error(err.empty() ? "initialisation produced no state" : err);
            init_.grid = std::make_shared<Grid>(); init_.grid->build(case_);
            if (int(init_.p.size()) != init_.grid->numActive()) throw std::runtime_error("thermal cell count differs from the grid");
            init_.depth = init_.grid->depth;
            init_.valid = true;
            init_.caseVersion = caseVersion_;
            prop_ = 6;
            notify("Initial state computed (thermal core)");
            return;
        }
        if (case_.comp.enabled) {   // compositional / CO2STORE: equilibrate with the EOS (or CO2-brine) model
            comp::DeckModel m; std::string err;
            if (!comp::deckModel(case_, m, &err)) throw std::runtime_error(err);
            const double pu = case_.units == UnitSystem::Field ? 6894.757293168 : 1e5;
            const bool co2 = case_.comp.co2store;
            bool got = false;
            std::vector<comp::CompSummary> out;
            comp::runDeckModel(m, out, [&](const comp::CompSummary& S, const comp::CompSimulator& sim) {
                const auto& P = sim.pressure();
                const std::vector<double> sg = sim.gasSaturation(), sw = sim.waterSaturation();
                const size_t n = P.size();
                init_.p.resize(n); init_.sw.resize(n); init_.sg.resize(n); init_.so.resize(n); init_.rs.assign(n, 0.0f);
                for (size_t q = 0; q < n; ++q) {
                    init_.p[q] = float(P[q] / pu); init_.sg[q] = float(sg[q]);
                    init_.sw[q] = float(co2 ? 1.0 - sg[q] : sw[q]);
                    init_.so[q] = float(co2 ? 0.0 : std::max(0.0, 1.0 - sw[q] - sg[q]));
                }
                init_.ip.oil = S.foip; init_.ip.water = S.fwip; init_.ip.gasFree = S.fgip;
                got = true;
                return false;   // initial state only
            }, &err, 0);
            if (!got) throw std::runtime_error(err.empty() ? "initialisation produced no state" : err);
            init_.grid = std::make_shared<Grid>(); init_.grid->build(case_);
            if (int(init_.p.size()) != init_.grid->numActive()) throw std::runtime_error("compositional cell count differs from the grid");
            init_.depth = init_.grid->depth;
            init_.valid = true;
            init_.caseVersion = caseVersion_;
            prop_ = 6;
            notify(co2 ? "Initial state computed (CO2-brine equilibrium)" : "Initial state computed (EOS equilibration)");
            return;
        }
        Simulator sim;
        sim.setup(case_);
        sim.initialize();
        if (sim.snapshots().empty()) throw std::runtime_error("initialisation produced no state");
        const Snapshot& s = sim.snapshots().front();
        init_.grid = std::make_shared<Grid>(sim.grid());
        init_.p = s.pressure; init_.sw = s.sw; init_.sg = s.sg; init_.so = s.so; init_.rs = s.rs;
        init_.ip = sim.inPlace();
        init_.depth = init_.grid->depth;
        init_.valid = true;
        init_.caseVersion = caseVersion_;
        prop_ = 6;                        // PRESSURE_0
        notify("Initial state computed");
    } catch (const std::exception& e) {
        init_.error = e.what();
        notify(std::string("Initialisation failed: ") + e.what());
    }
}

void App::runSimulation() {
    if (runner_.running()) return;
    runner_.setRestartSaving(settings_.rstEvery, settings_.rstLast);
    // decks with features the ResSim engine cannot model are run by OPM Flow / ECLIPSE on the original file
    const bool extOnly = !case_.engineUnsupported.empty() && !srcDeck_.empty();
    // E300 compositional decks (COMPS): the built-in engine runs them on its EOS compositional core
    if (case_.thermal.enabled && settings_.engine == 0) {   // E300 THERMAL case: STARS-like thermal core
        const std::string tv = validateThermal(case_);
        if (!tv.empty()) { notify("Cannot run: " + tv); return; }
        runner_.startThermal(case_);
        resultsCaseVersion_ = caseVersion_;
        snap_ = 0;
        if (prop_ < 9) prop_ = 9;
        wasRunning_ = true;
        return;
    }
    bool compDeck = case_.comp.enabled;   // compositional / CO2STORE case (imported or built in the UI)
    for (const auto& r : case_.engineUnsupported) if (r.rfind("compositional model", 0) == 0 && !srcDeck_.empty()) compDeck = true;
    if (compDeck && settings_.engine == 0) {
        runner_.startComp(case_, srcDeck_);
        resultsCaseVersion_ = caseVersion_;
        snap_ = 0;
        if (prop_ < 9) prop_ = 9;
        wasRunning_ = true;
        return;
    }
    if (extOnly && settings_.engine == 0) {
        showOpmSettings_ = true;
        notify("This deck can only be run by OPM Flow or ECLIPSE (" + case_.engineUnsupported.front() +
               ") - choose the engine in Run > Engine");
        return;
    }
    const std::string& v = validation();
    if (!v.empty() && !extOnly) {
        for (int s = 0; s < int(Section::Count); ++s) if (!sectionStatus(Section(s)).empty()) { section_ = Section(s); break; }
        notify("Cannot run: " + v);
        return;
    }
    if (case_.wells.empty()) runner_.addLog("WARN: no wells - the model will only show pressure equilibration.");
    if (settings_.engine == 2) {
        detectEclipse(false);
        if (eclDetected_.empty()) {
            showOpmSettings_ = true;
            notify("ECLIPSE not found - choose eclrun.exe or eclipse.exe in the engine settings");
            return;
        }
        OpmOptions o;
        o.kind = 1;
        o.flowPath = eclDetected_;
        o.eclVersion = settings_.eclVersion;
        o.eclProgram = settings_.eclProgram == 1 ? "e300" : settings_.eclProgram == 2 ? "frontsim" : "eclipse";
        if ((case_.thermal.enabled || case_.comp.enabled) && o.eclProgram != "e300") {   // THERMAL / COMPS decks need E300
            o.eclProgram = "e300";
            runner_.addLog("NOTE: thermal / compositional deck - running ECLIPSE 300 (e300)");
        }
        o.runDir = runDirFor(path_, case_.title, "ECL");
        o.extraArgs = settings_.eclArgs;
        applySourceDeck(o);
        runner_.startOpm(case_, o);
    } else if (settings_.engine == 1) {
        detectFlow(false);
        if (flowDetected_.empty()) {
            showOpmSettings_ = true;
            notify("OPM Flow not found - choose flow.exe in the engine settings");
            return;
        }
        OpmOptions o;
        o.flowPath = flowDetected_;
        o.runDir = runDirFor(path_, case_.title);
        o.threads = settings_.opmThreads;
        o.maxDtDays = settings_.opmMaxDt;
        o.enableTuning = settings_.opmTuning;
        o.keepWellsOpen = settings_.opmKeepWells;
        o.extraArgs = settings_.opmArgs;
        if (case_.grid.dual)
            notify("OPM Flow does not support dual porosity (DUALPORO/DUALPERM) and will stop - run this model with ResSim or ECLIPSE");
        applySourceDeck(o);
        runner_.startOpm(case_, o);
    } else {
        // C37 PYACTION: scripts are relative to the deck directory; the interpreter from the settings
        if (case_.deckDir.empty() && !path_.empty()) case_.deckDir = std::filesystem::u8path(path_).parent_path().u8string();
        if (case_.pythonExe.empty()) case_.pythonExe = settings_.pythonExe;
        if (!startMpiDD()) runner_.start(case_);
    }
    resultsCaseVersion_ = caseVersion_;
    snap_ = 0;
    if (prop_ < 9) prop_ = 9;             // PRESSURE (dynamic)
    wasRunning_ = true;
}

// Imported deck not edited since the import: the external engine runs the original file (with its INCLUDEs) unchanged.
void App::applySourceDeck(OpmOptions& o) {
    if (!case_.engineUnsupported.empty() && !srcDeck_.empty()) {
        o.sourceDeck = srcDeck_;
        if (caseVersion_ != srcDeckVersion_)
            runner_.addLog("NOTE: this deck uses features ResSim cannot represent - the ORIGINAL deck is run; model edits are not "
                           "included (edit the deck text in the Deck Editor instead).");
        return;
    }
    if (srcDeck_.empty() || settings_.deckRunMode != 0) return;
    if (caseVersion_ == srcDeckVersion_) { o.sourceDeck = srcDeck_; return; }
    runner_.addLog("NOTE: the imported deck was edited in ResSim - running the deck generated from the model. The original is "
                   "still at " + srcDeck_ + " (Deck Editor: open it to edit the deck text directly).");
}

void App::syncBuildersFromCase() {
    const GridInput& g = case_.grid;
    geom_.nx = g.nx; geom_.ny = g.ny; geom_.nz = g.nz;
    if (!g.dx.empty()) geom_.dx = g.dx[0];
    if (!g.dy.empty()) geom_.dy = g.dy[0];
    geom_.dz.assign(g.nz, 10.0);
    if (g.dz.size() == size_t(g.size())) for (int k = 0; k < g.nz; ++k) geom_.dz[k] = g.dz[g.index(0, 0, k)];
    if (!g.tops.empty()) {
        geom_.top = g.tops[0];
        const double rad = 180.0 / 3.14159265358979;
        geom_.dipX = (g.nx > 1 && g.tops.size() >= size_t(g.nx)) ? std::atan((g.tops[g.nx - 1] - g.tops[0]) / ((g.nx - 1) * geom_.dx)) * rad : 0.0;
        geom_.dipY = (g.ny > 1 && g.tops.size() >= size_t(g.nx * g.ny)) ? std::atan((g.tops[g.index(0, g.ny - 1, 0)] - g.tops[0]) / ((g.ny - 1) * geom_.dy)) * rad : 0.0;
        if (std::fabs(geom_.dipX) < 1e-6) geom_.dipX = 0;
        if (std::fabs(geom_.dipY) < 1e-6) geom_.dipY = 0;
    }
    rock_.layerPerm.assign(g.nz, 100); rock_.layerPoro.assign(g.nz, 0.2); rock_.layerKvKh.assign(g.nz, 0.1); rock_.layerNtg.assign(g.nz, 1.0);
    const int nl = g.nx * g.ny;
    if (g.permx.size() == size_t(g.size()) && g.poro.size() == size_t(g.size())) {
        for (int k = 0; k < g.nz; ++k) {
            double lk = 0, kz = 0, ph = 0, ntg = 0;
            for (int q = k * nl; q < (k + 1) * nl; ++q) {
                lk += std::log(std::max(g.permx[q], 1e-9));
                kz += g.permz.size() == size_t(g.size()) ? g.permz[q] / std::max(g.permx[q], 1e-12) : 0.1;
                ph += g.poro[q];
                ntg += g.ntg.size() == size_t(g.size()) ? g.ntg[q] : 1.0;
            }
            rock_.layerPerm[k] = std::exp(lk / nl);
            rock_.layerKvKh[k] = kz / nl;
            rock_.layerPoro[k] = ph / nl;
            rock_.layerNtg[k] = ntg / nl;
        }
    }
    rock_.sigmaLn = 0;
    const auto& sw = case_.props.swof;
    if (sw.size() >= 2) {
        corey_.swc = sw.front().s;
        double sKro0 = sw.back().s;                       // first Sw where krow reaches zero
        for (auto& r : sw) if (r.kro <= 1e-12) { sKro0 = r.s; break; }
        corey_.sorw = std::max(0.0, 1.0 - sKro0);
        corey_.krwMax = sw.back().krw; corey_.kroMax = sw.front().kro;
        corey_.pcEntry = sw.front().pc;
    }
    const auto& sg = case_.props.sgof;
    if (sg.size() >= 2) {
        corey_.krgMax = sg.back().krw;
        for (auto& r : sg) if (r.krw <= 0) corey_.sgc = r.s;
        double sKrog0 = sg.back().s;                      // first Sg where krog reaches zero
        for (auto& r : sg) if (r.kro <= 1e-12) { sKrog0 = r.s; break; }
        corey_.sorg = std::max(0.0, 1.0 - corey_.swc - sKrog0);
    }
}

void App::newCase(int which, const std::string& example) {
    Case ex;
    if (which == 4) {   // C13: a case built in the app (EGS designer); pendingCase_ must be set
        if (!pendingCase_) return;
        ex = std::move(*pendingCase_); pendingCase_.reset(); which = 3;
    } else if (which == 3) {
        try { ex = makeExample(example); }
        catch (const std::exception& e) { plat::msgBox(e.what(), "Example", MB_OK | MB_ICONERROR); return; }
    }
    runner_.cancel();
    if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
    sl_.res.reset();
    runner_.join();
    runner_.clear();
    if (which == 3) case_ = std::move(ex);
    else if (which == 0) case_ = makeDemoCase();
    else if (which == 1) case_ = makeBuckleyLeverettCase(100);
    else {
        Case c;
        c.title = "New model";
        c.units = UnitSystem::Field;
        geom_ = GeomInput{};
        geom_.nx = 15; geom_.ny = 15; geom_.nz = 3; geom_.dx = geom_.dy = 200; geom_.dz = {30, 30, 30}; geom_.top = 8000;
        rock_ = RockGenInput{};
        rock_.layerPerm = {150, 80, 250}; rock_.layerPoro = {0.22, 0.19, 0.24}; rock_.layerKvKh = {0.1, 0.1, 0.1}; rock_.layerNtg = {1, 1, 1};
        rock_.sigmaLn = 0.6;
        generateGeometry(geom_, c.grid);
        generateRock(rock_, c.grid);
        corr_ = CorrelationInput{};
        makeCorrelationPVT(corr_, c.units, c.props, false);
        corey_ = CoreyInput{};
        makeCoreySWOF(corey_, c.props.swof);
        makeCoreySGOF(corey_, c.props.sgof);
        c.equil.datumDepth = 8045; c.equil.datumPressure = 4600;
        c.equil.owc = 8200; c.equil.goc = 7900; c.equil.initialRs = -1;
        { ReportStep rs; rs.length = 91.3125; c.schedule.assign(20, rs); }
        case_ = c;
    }
    path_.clear();
    variants_.clear(); casesSel_ = -1;   // C40 case manager
    dirty_ = false;
    init_ = InitPreview{};
    if (which != 2) syncBuildersFromCase();
    ++caseVersion_;
    undoReset();
    resultsCaseVersion_ = 0;
    gridDirty_ = true;
    gridError_.clear();
    rebuildPreview();
    selWell_ = selStep_ = resultsWell_ = 0;
    prop_ = 1;                            // PERMX
    snap_ = 0;
    play_ = false;
    autoRange_ = true;
    viewGrid_ = nullptr;
    if (viewOk_) view_.resetCamera();
    runner_.addLog("=== New case: " + case_.title);
}

// Built-in example library, grouped Basic / Intermediate / Advanced, with a description tooltip per case.
void App::examplesMenu(bool popup) {
    (void)popup;
    static const std::vector<ExampleInfo> list = listExamples();
    const ExampleLevel levels[3] = {ExampleLevel::Basic, ExampleLevel::Intermediate, ExampleLevel::Advanced};
    for (ExampleLevel lv : levels) {
        ImGui::SeparatorText(toString(lv));
        for (const ExampleInfo& e : list) {
            if (e.level != lv) continue;
            if (ImGui::MenuItem(e.title.c_str(), e.grid.c_str())) { if (maybeSaveChanges()) newCase(3, e.id); }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort) && ImGui::BeginTooltip()) {
                ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
                ImGui::TextUnformatted(e.title.c_str());
                ImGui::Separator();
                ImGui::TextUnformatted(e.summary.c_str());
                ImGui::Spacing();
                ImGui::TextDisabled("%s", e.details.c_str());
                ImGui::Spacing();
                ImGui::TextDisabled("Grid: %s   |   typical run time: %s", e.grid.c_str(),
                                    e.runSeconds < 2 ? "about a second" : e.runSeconds < 20 ? "a few seconds" : "under a minute");
                ImGui::PopTextWrapPos();
                ImGui::EndTooltip();
            }
        }
    }
}

bool App::maybeSaveChanges() {
    if (!dirty_) return true;
    const std::string msg = "Save changes to \"" + (path_.empty() ? case_.title : baseName(path_)) + "\"?";
    const int r = plat::msgBox(msg.c_str(), "ResSim Studio", MB_YESNOCANCEL | MB_ICONQUESTION);
    if (r == IDCANCEL) return false;
    if (r == IDYES) { save(false); return !dirty_; }
    return true;
}

bool App::confirmClose() {
    if (runner_.running()) {
        if (plat::msgBox("A simulation is running. Stop it and exit?", "ResSim Studio", MB_OKCANCEL | MB_ICONWARNING) != IDOK) {
            quit_ = false;
            return false;
        }
        runner_.cancel();
        runner_.join();
    }
    if (!maybeSaveChanges()) { quit_ = false; return false; }
    if (!edConfirmClose()) { quit_ = false; return false; }
    return true;
}

void App::save(bool saveAs) {
    std::string p = path_;
    if (saveAs || p.empty()) {
        if (p.empty()) p = safeFileName(case_.title) + ".rsim";
        if (!fileDialog(true, "ResSim project (*.rsim)\0*.rsim\0All files (*.*)\0*.*\0", "rsim", p)) return;
    }
    try {
        saveProject(case_, p);
        casesWrite(p);   // C40 case manager variants (<project>.cases)
        path_ = p;
        dirty_ = false;
        recoveryClear();
        addRecent(p);
        notify("Saved " + baseName(p));
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "Save failed", MB_OK | MB_ICONERROR);
    }
}

void App::exportDeck() {
    std::string p;
    if (path_.empty()) p = safeFileName(case_.title) + ".DATA";
    else {
        const size_t dot = path_.find_last_of('.'), sl = path_.find_last_of("\\/");
        p = (dot != std::string::npos && (sl == std::string::npos || dot > sl) ? path_.substr(0, dot) : path_) + ".DATA";
    }
    if (!fileDialog(true, "ECLIPSE deck (*.DATA)\0*.DATA\0All files (*.*)\0*.*\0", "DATA", p)) return;
    try {
        exportEclipseDeck(case_, p);
        notify("Exported ECLIPSE deck " + baseName(p));
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "Export failed", MB_OK | MB_ICONERROR);
    }
}

void App::openDialog() {
    if (!maybeSaveChanges()) return;
    std::string p;
    if (!fileDialog(false, "ResSim project, results or ECLIPSE deck\0*.rsim;*.rsres;*.DATA;*.data\0ResSim project (*.rsim)\0*.rsim\0ResSim results (*.rsres)\0*.rsres\0ECLIPSE deck (*.DATA)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "rsim", p)) return;
    dirty_ = false;                       // already confirmed
    openFile(p);
}

void App::openFile(const std::string& path) {
    showStart_ = false;
    if (endsWithI(path, ".rsres")) { openResultsFile(path); addRecent(path); return; }
    for (const char* ext : {".SMSPEC", ".FSMSPEC", ".EGRID", ".FEGRID", ".UNRST", ".FUNRST", ".INIT", ".FINIT", ".UNSMRY", ".FUNSMRY"})
        if (endsWithI(path, ext)) { rvOpen(path); return; }
    if (endsWithI(path, ".DATA")) { if (maybeSaveChanges()) { dirty_ = false; importDeck(path); } return; }
    if (!maybeSaveChanges()) return;
    try {
        Case c = loadProject(path);
        const std::string v = validateCase(c);
        runner_.cancel(); runner_.join(); runner_.clear();
        case_ = std::move(c);
        path_ = path;
        casesRead(path);   // C40
        dirty_ = false;
        init_ = InitPreview{};
        if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
        sl_.res.reset();
        syncBuildersFromCase();
        ++caseVersion_;
        undoReset();
        resultsCaseVersion_ = 0;
        gridDirty_ = true;
        rebuildPreview();
        selWell_ = selStep_ = resultsWell_ = 0;
        prop_ = 1; snap_ = 0; play_ = false; viewGrid_ = nullptr;
        if (viewOk_) view_.resetCamera();
        runner_.addLog("=== Opened " + path);
        addRecent(path);
        notify(v.empty() ? "Opened " + baseName(path) : "Opened with problems: " + v);
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "Open failed", MB_OK | MB_ICONERROR);
    }
}

void App::importDeck(const std::string& path) {
    try {
        std::vector<std::string> warnings;
        std::vector<DeckMessage> report;
        Case c;
        try { c = importEclipseDeck(path, &warnings, &report); }
        catch (...) {
            deckReport_ = report; deckReportPath_ = path; deckSel_ = -1;
            if (!report.empty()) { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; }
            throw;
        }
        deckReport_ = std::move(report); deckReportPath_ = path; deckSel_ = -1;
        const std::string v = validateCase(c);
        runner_.cancel(); runner_.join(); runner_.clear();
        case_ = std::move(c);
        path_.clear();                    // a deck is imported, saving creates a .rsim project
        variants_.clear(); casesSel_ = -1;
        dirty_ = true;
        init_ = InitPreview{};
        if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
        sl_.res.reset();
        syncBuildersFromCase();
        ++caseVersion_;
        undoReset();
        resultsCaseVersion_ = 0;
        gridDirty_ = true;
        rebuildPreview();
        selWell_ = selStep_ = resultsWell_ = 0;
        prop_ = 1; snap_ = 0; play_ = false; viewGrid_ = nullptr;
        sl_.res.reset();
        if (viewOk_) view_.resetCamera();
        srcDeck_ = path; srcDeckVersion_ = caseVersion_;
        addRecent(path);
        runner_.addLog("=== Imported ECLIPSE deck " + path);
        for (auto& w : warnings) runner_.addLog("  deck warning: " + w);
        char b[160];
        std::snprintf(b, sizeof(b), "Imported %s: %dx%dx%d%s, %d wells, %d report steps, %d warnings", baseName(path).c_str(), case_.grid.nx, case_.grid.ny,
                      case_.grid.nz, case_.grid.cornerPoint() ? " corner-point" : "", int(case_.wells.size()), int(case_.schedule.size()), int(warnings.size()));
        if (!case_.engineUnsupported.empty())
            notify(std::string(b) + " - for OPM Flow / ECLIPSE only: " + case_.engineUnsupported.front() +
                   (case_.engineUnsupported.size() > 1 ? " (+" + std::to_string(case_.engineUnsupported.size() - 1) + " more, see the import report)" : std::string()));
        else notify(v.empty() ? std::string(b) : std::string(b) + " - problems: " + v);
        if (!warnings.empty()) { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; }
    } catch (const std::exception& e) {
        plat::msgBox(e.what(), "ECLIPSE deck import failed", MB_OK | MB_ICONERROR);
    }
}
