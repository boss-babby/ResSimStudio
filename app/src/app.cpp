// App: frame orchestration, docking layout, menus, toolbar, status bar, case explorer,
// inspector dispatch, log, file actions, run control and initialisation preview.
#include "app.h"
#include "icons.h"
#include "widgets.h"
#include "theme.h"
#include "gfx.h"
#include "imgui_internal.h"
#include "implot.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/examples.h"
#include <windows.h>
#include <commdlg.h>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <thread>
#include <cstring>
#include <vector>
#include <cctype>
#include <stdexcept>

using namespace ressim;

static const char* kSectionNames[] = {"Run specification", "Grid & rock", "Fluid PVT", "Rel. perm & Pc", "Initialization", "Aquifers", "EOR (polymer)", "Wells", "Schedule", "Numerics", "Advanced keywords"};
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

// Common open/save dialog. `path` may be a bare file name, a full path, or empty. Only the file-name part is
// sanitised (a drive letter's ':' stays). If the dialog refuses the initial name (GetSaveFileName fails silently
// with an invalid lpstrFile) it is reopened without one.
bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path) {
    std::string dir, name = path;
    const size_t slash = path.find_last_of("\\/");
    if (slash != std::string::npos) { dir = path.substr(0, slash); name = path.substr(slash + 1); }
    if (!name.empty()) name = safeFileName(name);
    if (!dir.empty()) {
        const DWORD attr = GetFileAttributesA(dir.c_str());
        if (attr == INVALID_FILE_ATTRIBUTES || !(attr & FILE_ATTRIBUTE_DIRECTORY)) dir.clear();
    }
    for (int attempt = 0; attempt < 2; ++attempt) {
        std::vector<char> buf(32768, 0);
        if (attempt == 0 && !name.empty() && name.size() < buf.size()) std::memcpy(buf.data(), name.data(), name.size());
        OPENFILENAMEA ofn{};
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = g_gfx.hwnd;
        ofn.lpstrFilter = filter;
        ofn.nFilterIndex = 1;
        ofn.lpstrFile = buf.data();
        ofn.nMaxFile = DWORD(buf.size());
        ofn.lpstrInitialDir = dir.empty() ? nullptr : dir.c_str();
        ofn.lpstrDefExt = defExt;
        ofn.Flags = OFN_EXPLORER | OFN_NOCHANGEDIR | OFN_HIDEREADONLY | OFN_PATHMUSTEXIST |
                    (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
        const BOOL ok = save ? GetSaveFileNameA(&ofn) : GetOpenFileNameA(&ofn);
        if (ok) { path = buf.data(); return true; }
        const DWORD err = CommDlgExtendedError();
        if (err == 0) return false;                            // the user cancelled
        if (attempt == 0) { dir.clear(); continue; }           // e.g. FNERR_INVALIDFILENAME (0x3002): retry with no name
        char msg[160];
        std::snprintf(msg, sizeof(msg), "The file dialog could not be opened (error 0x%04lX).", (unsigned long)err);
        MessageBoxA(g_gfx.hwnd, msg, "ResSim Studio", MB_OK | MB_ICONWARNING);
    }
    return false;
}

// ============================================================================ lifecycle
void App::init() {
    viewOk_ = view_.init();
    loadSettings();
    if (settings_.uiScale <= 0) {
        // first start: laptop screens (1366 / 1440 / 1536 logical px) get a compact interface
        const float logical = float(GetSystemMetrics(SM_CXSCREEN)) / std::max(0.5f, g_gfx.dpiScale);
        settings_.uiScale = logical < 1500 ? 0.85f : logical < 1800 ? 0.92f : 1.0f;
    }
    setUiScale(settings_.uiScale);
    newCase(0);
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
}

bool App::wantsContinuousFrames() const {
    return runner_.running() || study_.running || study_.manualPending || play_ || ImGui::GetTime() < toastUntil_;
}

void App::frame() {
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
        if (err.empty()) notify(opmDeck.empty() ? std::string("Simulation finished") : extName + " run finished");
        else notify("Run stopped: " + err);
        if (opmOk && settings_.opmAutoView && !opmDeck.empty()) rvOpen(opmDeck, true);
    }
    wasRunning_ = running;

    // global shortcuts (the Deck Editor handles its own file shortcuts and F5)
    const ImGuiInputFlags G = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_1, G)) workspace_ = Workspace::Model;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_2, G)) workspace_ = Workspace::DeckEditor;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_3, G)) workspace_ = Workspace::Results;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_4, G)) workspace_ = Workspace::Studies;
    if (ImGui::Shortcut(ImGuiMod_Shift | ImGuiKey_F5, G)) runner_.cancel();
    if (workspace_ != Workspace::DeckEditor) {
        if (ImGui::Shortcut(ImGuiKey_F5, G)) runSimulation();
        {   // interface size: Ctrl + / Ctrl - / Ctrl 0
            float k = settings_.uiScale;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Equal, G) || ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_KeypadAdd, G)) k += 0.05f;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Minus, G) || ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_KeypadSubtract, G)) k -= 0.05f;
            if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_0, G)) k = 1.0f;
            k = std::clamp(k, 0.6f, 1.6f);
            if (std::fabs(k - settings_.uiScale) > 1e-4f) { settings_.uiScale = k; setUiScale(k); saveSettings(); }
        }
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_S, G)) save(false);
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_S, G)) save(true);
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_O, G)) openDialog();
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_N, G)) { if (maybeSaveChanges()) newCase(2); }
    }

    menuBar();
    toolbar();
    statusBar();
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
        else resultsWorkspace();
    }
    figureExportWindow();
    opmSettingsWindow();
    if (rvDeferred_) { auto f = std::move(rvDeferred_); rvDeferred_ = nullptr; f(); }
    aboutWindow();
    if (showImGuiDemo_) ImGui::ShowDemoWindow(&showImGuiDemo_);

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
        dl->AddRectFilled(p0, p1, ImGui::GetColorU32(ImVec4(0.13f, 0.14f, 0.17f, 0.97f * alpha)), 8);
        dl->AddRectFilled(p0, ImVec2(p0.x + 5, p1.y), (Pal::Accent & 0x00FFFFFF) | (ImU32(255 * alpha) << 24), 8, ImDrawFlags_RoundCornersLeft);
        dl->AddText(ImVec2(p0.x + pad + 6, p0.y + pad), ImGui::GetColorU32(ImVec4(1, 1, 1, alpha)), toast_.c_str());
    }

    // window title
    std::string title = "ResSim Studio - " + (path_.empty() ? case_.title : baseName(path_)) + (dirty_ ? " *" : "");
    if (std::string(titleBuf_.begin(), titleBuf_.end()) != title) {
        titleBuf_.assign(title.begin(), title.end());
        SetWindowTextA(g_gfx.hwnd, title.c_str());
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
    if (ImGui::BeginMenu("File")) {
        if (ImGui::BeginMenu(ICON_FILE_PLUS "  New")) {
            if (ImGui::MenuItem("Demo: 3-phase waterflood (SPE1-like PVT)")) { if (maybeSaveChanges()) newCase(0); }
            if (ImGui::MenuItem("Validation: 1D Buckley-Leverett")) { if (maybeSaveChanges()) newCase(1); }
            if (ImGui::MenuItem("Blank model", "Ctrl+N")) { if (maybeSaveChanges()) newCase(2); }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu(ICON_BOOK_OPEN "  Examples")) { examplesMenu(false); ImGui::EndMenu(); }
        if (ImGui::MenuItem(ICON_FOLDER_OPEN "  Open...", "Ctrl+O")) openDialog();
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_SAVE "  Save", "Ctrl+S")) save(false);
        if (ImGui::MenuItem(ICON_SAVE "  Save As...", "Ctrl+Shift+S")) save(true);
        if (ImGui::MenuItem(ICON_FILE_CODE "  Edit deck file in Deck Editor...")) {
            std::string p;
            if (fileDialog(false, "ECLIPSE deck / include (*.DATA;*.INC;*.GRDECL)\0*.DATA;*.data;*.INC;*.inc;*.GRDECL;*.grdecl;*.SCH;*.sch;*.PROP;*.prop\0All files (*.*)\0*.*\0", "DATA", p)) {
                edOpenFile(p);
                workspace_ = Workspace::DeckEditor;
            }
        }
        if (ImGui::MenuItem("Open current model as deck text (Deck Editor)")) { edOpenModelDeck(); workspace_ = Workspace::DeckEditor; }
        if (ImGui::MenuItem(ICON_UPLOAD "  Import ECLIPSE deck (.DATA)...")) {
            std::string p;
            if (maybeSaveChanges() && fileDialog(false, "ECLIPSE deck (*.DATA;*.data)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "DATA", p)) { dirty_ = false; importDeck(p); }
        }
        if (ImGui::MenuItem(ICON_DOWNLOAD "  Export ECLIPSE deck (.DATA)...")) exportDeck();
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_CHART_LINE "  Open results (ECLIPSE / OPM Flow)...")) {
            std::string p;
            if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.FEGRID;*.UNRST;*.FUNRST;*.INIT;*.FINIT;*.UNSMRY;*.FUNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p))
                rvOpen(p);
        }
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_LOG_OUT "  Exit", "Alt+F4")) quit_ = true;
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Run")) {
        if (ImGui::MenuItem(ICON_PLAY "  Run simulation", "F5", false, !runner_.running())) runSimulation();
        if (ImGui::MenuItem(ICON_SQUARE "  Stop", "Shift+F5", false, runner_.running())) runner_.cancel();
        ImGui::Separator();
        if (ImGui::MenuItem("Engine: ResSim (built-in)", nullptr, settings_.engine == 0)) { settings_.engine = 0; saveSettings(); }
        if (ImGui::MenuItem("Engine: OPM Flow (external)", nullptr, settings_.engine == 1)) { settings_.engine = 1; saveSettings(); detectFlow(false); }
        if (ImGui::MenuItem("Engine: ECLIPSE (external, licensed)", nullptr, settings_.engine == 2)) { settings_.engine = 2; saveSettings(); detectEclipse(false); }
        if (ImGui::MenuItem(ICON_SETTINGS "  Simulation engine settings (OPM Flow / ECLIPSE)...")) showOpmSettings_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_ARROW_DOWN_TO_LINE "  Compute initial state")) { computeInit(); section_ = Section::Init; }
        if (ImGui::MenuItem(ICON_CHECK "  Validate inputs")) {
            const std::string& v = validation();
            notify(v.empty() ? "All inputs are valid" : v);
        }
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("View")) {
        if (ImGui::MenuItem(ICON_BOX "  Model workspace", "Ctrl+1", workspace_ == Workspace::Model)) workspace_ = Workspace::Model;
        if (ImGui::MenuItem(ICON_FILE_CODE "  Deck Editor (writing mode)", "Ctrl+2", workspace_ == Workspace::DeckEditor)) workspace_ = Workspace::DeckEditor;
        if (ImGui::MenuItem(ICON_CHART_LINE "  Results workspace", "Ctrl+3", workspace_ == Workspace::Results)) workspace_ = Workspace::Results;
        if (ImGui::MenuItem(ICON_FLASK_CONICAL "  Studies workspace (sensitivity / HM)", "Ctrl+4", workspace_ == Workspace::Studies)) workspace_ = Workspace::Studies;
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_LAYOUT_DASHBOARD "  Reset window layout")) resetLayout_ = true;
        if (ImGui::MenuItem(ICON_MAXIMIZE "  Reset 3D camera") && viewOk_) view_.resetCamera();
        if (ImGui::BeginMenu(ICON_ZOOM_IN "  Interface size")) {
            static const float sizes[] = {0.70f, 0.80f, 0.90f, 1.00f, 1.15f, 1.30f, 1.50f};
            static const char* names[] = {"70 %  (very compact)", "80 %  (compact)", "90 %", "100 %  (default)", "115 %", "130 %", "150 %"};
            for (int i = 0; i < 7; ++i)
                if (ImGui::MenuItem(names[i], nullptr, std::fabs(settings_.uiScale - sizes[i]) < 0.01f)) { settings_.uiScale = sizes[i]; setUiScale(sizes[i]); saveSettings(); }
            ImGui::Separator();
            ImGui::TextDisabled("Ctrl + / Ctrl - / Ctrl 0");
            ImGui::EndMenu();
        }
        if (ImGui::MenuItem(ICON_SETTINGS "  Settings...")) showOpmSettings_ = true;
        ImGui::Separator();
        if (ImGui::MenuItem("Result Viewer (ECLIPSE / OPM output, compare)", nullptr, showResultViewer_)) showResultViewer_ = !showResultViewer_;
        if (ImGui::MenuItem("Deck import report", nullptr, showDeckTools_ && deckTab_ == 0)) { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; }
        if (ImGui::MenuItem("ECLIPSE keyword reference", nullptr, showDeckTools_ && deckTab_ == 1)) { showDeckTools_ = true; deckTab_ = 1; deckTabSet_ = true; }
        ImGui::Separator();
        ImGui::MenuItem("UI demo windows (developer)", nullptr, &showImGuiDemo_);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Tools")) {
        if (ImGui::MenuItem(ICON_FLASK_CONICAL "  Studies workspace: sensitivity & history matching", "Ctrl+4")) workspace_ = Workspace::Studies;
        if (ImGui::MenuItem(ICON_PANEL_BOTTOM "  Studies panel (Model workspace)", nullptr, showStudies_)) { showStudies_ = true; ImGui::SetWindowFocus("Studies"); }
        if (ImGui::MenuItem(ICON_ROUTE "  Streamlines in 3D viewer", nullptr, sl_.show)) sl_.show = !sl_.show;
        ImGui::Separator();
        if (ImGui::MenuItem("Convert grid to corner-point (COORD/ZCORN)", nullptr, false, !case_.grid.cornerPoint())) {
            convertToCornerPoint(case_.grid); markDirty(true); notify("Grid converted to corner-point geometry");
        }
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Help")) {
        if (ImGui::MenuItem(ICON_INFO "  About ResSim Studio")) showAbout_ = true;
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
        if (ImGui::Button(ICON_FILE_PLUS " New")) ImGui::OpenPopup("newpop");
        tip("New case");
        if (ImGui::BeginPopup("newpop")) {
            if (ImGui::Selectable("Demo: 3-phase waterflood")) { if (maybeSaveChanges()) newCase(0); }
            if (ImGui::Selectable("Validation: 1D Buckley-Leverett")) { if (maybeSaveChanges()) newCase(1); }
            if (ImGui::Selectable("Blank model")) { if (maybeSaveChanges()) newCase(2); }
            ImGui::Separator();
            if (ImGui::BeginMenu("Examples")) { examplesMenu(true); ImGui::EndMenu(); }
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button(ICON_FOLDER_OPEN " Open")) openDialog();
        tip("Open project (Ctrl+O)");
        ImGui::SameLine();
        if (ImGui::Button(ICON_SAVE " Save")) save(false);
        tip("Save project (Ctrl+S)");
        ImGui::SameLine();
        ImGui::TextDisabled("|");
        ImGui::SameLine();
        workspaceSwitch();
        ImGui::SameLine();
        ImGui::TextDisabled("|");
        ImGui::SameLine();
        const bool running = runner_.running();
        ImGui::BeginDisabled(running);
        if (ui::accentButton(ICON_PLAY "  Run")) runSimulation();
        ImGui::EndDisabled();
        tip("Run the full simulation (F5)");
        ImGui::SameLine();
        ImGui::BeginDisabled(!running);
        if (ui::dangerButton(ICON_SQUARE "  Stop")) runner_.cancel();
        ImGui::EndDisabled();
        tip("Stop the running simulation (Shift+F5)");
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
            tip("Simulation engine used by Run: ResSim (built-in), OPM Flow (external flow.exe) or ECLIPSE (licensed eclrun / eclipse.exe). Settings: Run menu.");
            ImGui::SameLine();
            if (ImGui::Button(ICON_SETTINGS "##eng")) showOpmSettings_ = true;
            tip("Simulation engine settings (OPM Flow / ECLIPSE executable, auto-detect, options)");
            ImGui::SameLine();
        }
        if (ImGui::Button(ICON_ARROW_DOWN_TO_LINE " Initialise")) { computeInit(); section_ = Section::Init; }
        tip("Compute the equilibrium initial state and fluids in place");
        ImGui::SameLine();
        if (ImGui::Button(ICON_CHECK " Validate")) {
            const std::string& v = validation();
            notify(v.empty() ? "All inputs are valid" : v);
        }
        tip("Check all inputs for consistency");
        ImGui::SameLine();
        if (running) {
            ImGui::SameLine();
            char b[64]; std::snprintf(b, sizeof(b), "%.0f%%  (%.1f s)", runner_.progress() * 100.0, runner_.wallSeconds());
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col4(Pal::Accent));
            ImGui::ProgressBar(float(runner_.progress()), ImVec2(ImGui::GetFontSize() * 14, 0), b);
            ImGui::PopStyleColor();
        }
        // right side: validation indicator
        const std::string& v = validation();
        const char* st = v.empty() ? ICON_CHECK " Inputs valid" : ICON_TRIANGLE_ALERT " Input problem";
        const float w = ImGui::CalcTextSize(st).x;
        ImGui::SameLine(ImGui::GetWindowWidth() - w - ImGui::GetStyle().WindowPadding.x * 2);
        ImGui::TextColored(col4(v.empty() ? Pal::Ok : Pal::Warn), "%s", st);
        if (!v.empty() && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", v.c_str());
    }
    ImGui::End();
    ImGui::PopStyleVar();
}

void App::statusBar() {
    const float h = ImGui::GetFrameHeight();
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(ImGui::GetStyle().WindowPadding.x, 2));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.055f, 0.06f, 0.07f, 1));
    if (ImGui::BeginViewportSideBar("##status", ImGui::GetMainViewport(), ImGuiDir_Down, h, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings)) {
        if (runner_.running()) ImGui::TextColored(col4(Pal::Accent), ICON_ACTIVITY " Running  %.0f%%", runner_.progress() * 100.0);
        else ImGui::TextColored(col4(Pal::Ok), ICON_CIRCLE_DOT " Ready");
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
        char r[96];
        const int thr = case_.tuning.threads > 0 ? case_.tuning.threads : int(std::max(1u, std::thread::hardware_concurrency()));
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
    ui::header("Model inputs");
    for (int s = 0; s < int(Section::Count); ++s) {
        const std::string st = sectionStatus(Section(s));
        ImGui::PushID(s);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        char lb[96]; std::snprintf(lb, sizeof(lb), "   %s  %s", kSectionIcons[s], kSectionNames[s]);
        if (ImGui::Selectable(lb, section_ == Section(s), 0, ImVec2(0, ImGui::GetFrameHeight()))) section_ = Section(s);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s%s%s", kSectionHelp[s], st.empty() ? "" : "\n\nProblem: ", st.c_str());
        const float r = ImGui::GetFontSize() * 0.22f;
        ImGui::GetWindowDrawList()->AddCircleFilled(ImVec2(p.x + r + 2, p.y + ImGui::GetFrameHeight() * 0.5f), r, st.empty() ? Pal::Ok : Pal::Err);
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
    ImGui::Text("%s  %s", kSectionIcons[s], kSectionNames[s]);
    if (g_fonts.title) ImGui::PopFont();
    ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
    ImGui::TextWrapped("%s", kSectionHelp[s]);
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
    ImGui::Separator();
    ImGui::BeginChild("editor", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_None);
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
            "matching, Levenberg-Marquardt and ES-MDA assisted history matching. Pollock streamlines with "
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
    const std::string& v = validation();
    if (!v.empty()) {
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
    } else
    runner_.start(case_);
    resultsCaseVersion_ = caseVersion_;
    snap_ = 0;
    if (prop_ < 9) prop_ = 9;             // PRESSURE (dynamic)
    wasRunning_ = true;
}

// Imported deck not edited since the import: the external engine runs the original file (with its INCLUDEs) unchanged.
void App::applySourceDeck(OpmOptions& o) {
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
    if (which == 3) {
        try { ex = makeExample(example); }
        catch (const std::exception& e) { MessageBoxA(g_gfx.hwnd, e.what(), "Example", MB_OK | MB_ICONERROR); return; }
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
        c.schedule.assign(20, ReportStep{91.3125, {}});
        case_ = c;
    }
    path_.clear();
    dirty_ = false;
    init_ = InitPreview{};
    if (which != 2) syncBuildersFromCase();
    ++caseVersion_;
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
    const int r = MessageBoxA(g_gfx.hwnd, msg.c_str(), "ResSim Studio", MB_YESNOCANCEL | MB_ICONQUESTION);
    if (r == IDCANCEL) return false;
    if (r == IDYES) { save(false); return !dirty_; }
    return true;
}

bool App::confirmClose() {
    if (runner_.running()) {
        if (MessageBoxA(g_gfx.hwnd, "A simulation is running. Stop it and exit?", "ResSim Studio", MB_OKCANCEL | MB_ICONWARNING) != IDOK) {
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
        path_ = p;
        dirty_ = false;
        notify("Saved " + baseName(p));
    } catch (const std::exception& e) {
        MessageBoxA(g_gfx.hwnd, e.what(), "Save failed", MB_OK | MB_ICONERROR);
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
        MessageBoxA(g_gfx.hwnd, e.what(), "Export failed", MB_OK | MB_ICONERROR);
    }
}

void App::openDialog() {
    if (!maybeSaveChanges()) return;
    std::string p;
    if (!fileDialog(false, "ResSim project or ECLIPSE deck\0*.rsim;*.DATA;*.data\0ResSim project (*.rsim)\0*.rsim\0ECLIPSE deck (*.DATA)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "rsim", p)) return;
    dirty_ = false;                       // already confirmed
    openFile(p);
}

void App::openFile(const std::string& path) {
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
        dirty_ = false;
        init_ = InitPreview{};
        if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
        sl_.res.reset();
        syncBuildersFromCase();
        ++caseVersion_;
        resultsCaseVersion_ = 0;
        gridDirty_ = true;
        rebuildPreview();
        selWell_ = selStep_ = resultsWell_ = 0;
        prop_ = 1; snap_ = 0; play_ = false; viewGrid_ = nullptr;
        if (viewOk_) view_.resetCamera();
        runner_.addLog("=== Opened " + path);
        notify(v.empty() ? "Opened " + baseName(path) : "Opened with problems: " + v);
    } catch (const std::exception& e) {
        MessageBoxA(g_gfx.hwnd, e.what(), "Open failed", MB_OK | MB_ICONERROR);
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
        dirty_ = true;
        init_ = InitPreview{};
        if (!study_.running) { study_.hasTornado = study_.hasMC = study_.hasHM = study_.hasManual = false; }
        sl_.res.reset();
        syncBuildersFromCase();
        ++caseVersion_;
        resultsCaseVersion_ = 0;
        gridDirty_ = true;
        rebuildPreview();
        selWell_ = selStep_ = resultsWell_ = 0;
        prop_ = 1; snap_ = 0; play_ = false; viewGrid_ = nullptr;
        sl_.res.reset();
        if (viewOk_) view_.resetCamera();
        srcDeck_ = path; srcDeckVersion_ = caseVersion_;
        runner_.addLog("=== Imported ECLIPSE deck " + path);
        for (auto& w : warnings) runner_.addLog("  deck warning: " + w);
        char b[160];
        std::snprintf(b, sizeof(b), "Imported %s: %dx%dx%d%s, %d wells, %d report steps, %d warnings", baseName(path).c_str(), case_.grid.nx, case_.grid.ny,
                      case_.grid.nz, case_.grid.cornerPoint() ? " corner-point" : "", int(case_.wells.size()), int(case_.schedule.size()), int(warnings.size()));
        notify(v.empty() ? std::string(b) : std::string(b) + " - problems: " + v);
        if (!warnings.empty()) { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; }
    } catch (const std::exception& e) {
        MessageBoxA(g_gfx.hwnd, e.what(), "ECLIPSE deck import failed", MB_OK | MB_ICONERROR);
    }
}
