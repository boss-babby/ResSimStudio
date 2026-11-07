// C26 usability layer (general users, engineers, researchers):
//  * Command palette (Ctrl+Shift+P, also Ctrl+K): every menu action, model section, example, example deck and recent
//    file by name - fuzzy search, arrow keys + Enter. Nothing has to be found in the menus any more.
//  * Start page (Help > Start page, shown at start-up unless switched off): new case, open / import, recent files, the
//    built-in example library and the example DECKS shipped next to the executable (examples\, decks\) with their
//    description (first comment lines) - one click opens a deck.
//  * Recent files (projects, decks, results): File > Recent, start page, palette; kept in the settings file.
//  * Workflow checklist (Case Explorer): grid -> fluids -> rel. perm -> initial state -> wells -> schedule -> validate ->
//    run -> results -> save, with the state of each step (done / needs attention / stale results) and a click to go there.
//  * Notification history (status bar): every toast with its time, errors marked; toasts no longer get lost.
//  * Run history: every finished ResSim run is kept in the Result Viewer compare list (summary vectors, report maps for
//    small grids) with WHAT CHANGED against the previous run (input lines removed / added) - quick what-if comparisons
//    without saving result files by hand. Settings: View > Settings > Interface.
#include "app.h"
#include "i18n.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "ressim/examples.h"
#include "ressim/io.h"
#include "ressim/ecl.h"
#include "plotkit.h"
#include "platform.h"
#include <cfloat>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>

using namespace ressim;
namespace fs = std::filesystem;
std::vector<std::pair<std::string, std::string>> keywordIndexList();   // ui_deck.cpp

namespace {
std::string exeFolder() { return fs::u8path(plat::exePath()).parent_path().u8string(); }
std::string lower(std::string s) { for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c))); return s; }
std::string fileOnly(const std::string& p) { const size_t s = p.find_last_of("\\/"); return s == std::string::npos ? p : p.substr(s + 1); }

// fuzzy subsequence score: < 0 = no match; consecutive and word-start matches score higher, shorter labels win ties
int fuzzyScore(const std::string& label, const std::string& q) {
    if (q.empty()) return 1;
    const std::string L = lower(label), Q = lower(q);
    if (L.find(Q) != std::string::npos) return 1000 - int(L.find(Q)) * 2 - int(L.size()) / 8;   // substring: best
    int score = 0; size_t li = 0; int run = 0;
    for (char c : Q) {
        if (c == ' ') continue;
        bool found = false;
        while (li < L.size()) {
            if (L[li] == c) {
                const bool wordStart = li == 0 || !std::isalnum(static_cast<unsigned char>(L[li - 1]));
                score += 10 + (wordStart ? 15 : 0) + run * 5;
                ++run; ++li; found = true; break;
            }
            run = 0; ++li;
        }
        if (!found) return -1;
    }
    return score - int(L.size()) / 8;
}

std::string timeStr(std::time_t t) { char b[16]; const std::tm* tmv = std::localtime(&t); if (!tmv) return "--:--:--"; std::strftime(b, sizeof b, "%H:%M:%S", tmv); return b; }

// description of an example deck: the leading "--" comment lines (without the dashes), at most ~400 characters
std::string deckDescription(const std::string& path) {
    std::ifstream f(fs::u8path(path));
    std::string line, d;
    int n = 0;
    while (std::getline(f, line) && n < 12) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("--", 0) != 0) { if (line.find_first_not_of(" \t") == std::string::npos) continue; break; }
        std::string t = line.substr(2);
        while (!t.empty() && (t[0] == ' ' || t[0] == '-')) t.erase(0, 1);
        if (t.empty()) continue;
        if (!d.empty()) d += ' ';
        d += t; ++n;
        if (d.size() > 420) { d.resize(420); d += "..."; break; }
    }
    return d;
}

// "what changed" between two project texts: lines removed / added (multiset difference, order kept)
std::string textDiff(const std::string& a, const std::string& b, int maxLines = 14) {
    auto split = [](const std::string& s) { std::vector<std::string> v; size_t p = 0; while (p <= s.size()) { size_t e = s.find('\n', p); if (e == std::string::npos) e = s.size(); std::string l = s.substr(p, e - p); if (!l.empty() && l.back() == '\r') l.pop_back(); if (l.find_first_not_of(" \t") != std::string::npos) v.push_back(l); p = e + 1; } return v; };
    const auto A = split(a), B = split(b);
    std::map<std::string, int> ca, cb;
    for (const auto& l : A) ++ca[l];
    for (const auto& l : B) ++cb[l];
    std::vector<std::string> out;
    int removed = 0, added = 0;
    std::map<std::string, int> used;
    for (const auto& l : A) { if (used[l] < ca[l] - cb[l]) { ++used[l]; ++removed; if (int(out.size()) < maxLines) out.push_back("- " + l.substr(0, 110)); } }
    used.clear();
    for (const auto& l : B) { if (used[l] < cb[l] - ca[l]) { ++used[l]; ++added; if (int(out.size()) < maxLines) out.push_back("+ " + l.substr(0, 110)); } }
    if (removed + added == 0) return "no input change";
    std::string s;
    for (const auto& l : out) s += l + "\n";
    if (removed + added > int(out.size())) s += "... (" + std::to_string(removed) + " lines removed, " + std::to_string(added) + " added in total)\n";
    return s;
}
} // namespace

std::string projectTextDiff(const std::string& a, const std::string& b, int maxLines) { return textDiff(a, b, maxLines); }   // C40 case manager

// ============================================================================ recent files
void App::addRecent(const std::string& path) {
    if (path.empty()) return;
    auto& R = settings_.recent;
    R.erase(std::remove_if(R.begin(), R.end(), [&](const std::string& p) { return lower(p) == lower(path); }), R.end());
    R.insert(R.begin(), path);
    if (R.size() > 12) R.resize(12);
    saveSettings();
}

void App::recentMenu() {
    if (settings_.recent.empty()) { ImGui::TextDisabled("(no recent files)"); return; }
    std::string open;
    for (const std::string& p : settings_.recent) {
        std::error_code ec;
        const bool ok = fs::exists(fs::u8path(p), ec);
        if (ImGui::MenuItem(fileOnly(p).c_str(), nullptr, false, ok)) open = p;
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("%s%s", p.c_str(), ok ? "" : "\n(file not found)");
    }
    ImGui::Separator();
    if (ImGui::MenuItem("Clear the list")) { settings_.recent.clear(); saveSettings(); }
    if (!open.empty()) openFile(open);
}

// ============================================================================ example decks next to the executable
void App::scanDeckFiles() {
    if (deckFilesScanned_) return;
    deckFilesScanned_ = true;
    deckFiles_.clear();
    const std::string root = exeFolder();
    for (const char* sub : {"examples", "decks"}) {
        std::error_code ec;
        fs::path d = fs::u8path(root) / sub;
        // C44: the release folder ships without copies of examples/ and decks/ - they are found one level up (ressim/)
        if (!fs::is_directory(d, ec)) d = fs::u8path(root).parent_path() / sub;
        if (!fs::is_directory(d, ec)) d = fs::u8path(root).parent_path().parent_path() / sub;
        if (!fs::is_directory(d, ec)) continue;
        for (fs::recursive_directory_iterator it(d, ec), end; it != end && !ec; it.increment(ec)) {
            if (!it->is_regular_file(ec)) continue;
            const std::string ext = lower(it->path().extension().u8string());
            if (ext != ".data") continue;
            DeckFile f;
            f.path = it->path().u8string();
            f.name = it->path().stem().u8string();
            f.folder = sub;
            const fs::path rel = fs::relative(it->path().parent_path(), d, ec);
            if (!ec && !rel.empty() && rel != ".") f.folder += "/" + rel.generic_u8string();
            f.desc = deckDescription(f.path);
            deckFiles_.push_back(std::move(f));
        }
    }
    std::sort(deckFiles_.begin(), deckFiles_.end(), [](const DeckFile& a, const DeckFile& b) { return a.folder != b.folder ? a.folder < b.folder : a.name < b.name; });
}

void App::deckFilesMenu() {
    scanDeckFiles();
    if (deckFiles_.empty()) { ImGui::TextDisabled("(no examples\\ or decks\\ folder next to ResSimStudio.exe)"); return; }
    std::string folder, open;
    for (const DeckFile& f : deckFiles_) {
        if (f.folder != folder) { folder = f.folder; ImGui::SeparatorText(folder.c_str()); }
        if (ImGui::MenuItem(f.name.c_str())) open = f.path;
        if (!f.desc.empty() && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
            ImGui::BeginTooltip(); ImGui::PushTextWrapPos(ImGui::GetFontSize() * 32); ImGui::TextUnformatted(f.desc.c_str()); ImGui::PopTextWrapPos(); ImGui::EndTooltip();
        }
    }
    if (!open.empty()) openFile(open);
}

// ============================================================================ notifications
void App::pushNote(const std::string& msg) {
    Note n; n.when = (long long)std::time(nullptr); n.text = msg;
    const std::string l = lower(msg);
    n.error = l.find("fail") != std::string::npos || l.find("error") != std::string::npos || l.find("stopped") != std::string::npos || l.find("problem") != std::string::npos;
    notes_.push_back(std::move(n));
    if (notes_.size() > 200) notes_.erase(notes_.begin());
    ++notesUnread_;
}

void App::notificationButton() {
    char b[48];
    std::snprintf(b, sizeof b, notesUnread_ > 0 ? ICON_HISTORY " %d##notes" : ICON_HISTORY "##notes", notesUnread_);
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
    if (ImGui::SmallButton(b) || notesPopupReq_) { ImGui::OpenPopup("##notespop"); notesUnread_ = 0; notesPopupReq_ = false; }
    ImGui::PopStyleColor();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Notification history (%zu)", notes_.size());
    ImGui::SetNextWindowSizeConstraints(ImVec2(ImGui::GetFontSize() * 24, 0), ImVec2(ImGui::GetFontSize() * 48, ImGui::GetFontSize() * 30));
    if (ImGui::BeginPopup("##notespop")) {
        ImGui::TextDisabled("Notifications (newest first)");
        ImGui::SameLine();
        if (ImGui::SmallButton("Clear")) notes_.clear();
        ImGui::Separator();
        if (notes_.empty()) ImGui::TextDisabled("(none yet)");
        for (int i = int(notes_.size()) - 1; i >= 0; --i) {
            const Note& n = notes_[size_t(i)];
            ImGui::TextDisabled("%s", timeStr(std::time_t(n.when)).c_str());
            ImGui::SameLine();
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 46);
            if (n.error) ImGui::TextColored(col4(Pal::Err), "%s", n.text.c_str());
            else ImGui::TextUnformatted(n.text.c_str());
            ImGui::PopTextWrapPos();
        }
        ImGui::EndPopup();
    }
}

// ============================================================================ workflow checklist
void App::workflowChecklist() {
    struct Step { const char* label; int state; std::string why; std::function<void()> go; };   // state 0 todo, 1 done, 2 attention
    std::vector<Step> S;
    const std::string grid = sectionStatus(Section::Grid), pvt = sectionStatus(Section::PVT), scal = sectionStatus(Section::SCAL);
    const std::string ini = sectionStatus(Section::Init), wel = sectionStatus(Section::Wells), sch = sectionStatus(Section::Schedule);
    auto sec = [this](Section s) { return [this, s] { workspace_ = Workspace::Model; section_ = s; }; };
    S.push_back({TRN("Grid & rock"), grid.empty() && preview_ && preview_->numActive() > 0 ? 1 : 2, grid.empty() ? (preview_ ? std::to_string(preview_->numActive()) + " active cells" : "grid not built") : grid, sec(Section::Grid)});
    S.push_back({TRN("Fluid properties"), pvt.empty() ? 1 : 2, pvt, sec(Section::PVT)});
    S.push_back({TRN("Rel. perm & capillary pressure"), scal.empty() ? 1 : 2, scal, sec(Section::SCAL)});
    S.push_back({TRN("Initial state"), ini.empty() ? 1 : 2, ini.empty() ? "equilibration inputs complete (Initialise to preview fluids in place)" : ini, sec(Section::Init)});
    S.push_back({TRN("Wells"), !wel.empty() ? 2 : case_.wells.empty() ? 0 : 1, !wel.empty() ? wel : std::to_string(case_.wells.size()) + " wells", sec(Section::Wells)});
    S.push_back({TRN("Schedule"), !sch.empty() ? 2 : case_.schedule.empty() ? 0 : 1, !sch.empty() ? sch : std::to_string(case_.schedule.size()) + " report steps", sec(Section::Schedule)});
    const std::string& v = validation();
    S.push_back({TRN("Validate inputs"), v.empty() ? 1 : 2, v.empty() ? "all inputs consistent" : v, [this] { const std::string& vv = validation(); notify(vv.empty() ? "All inputs are valid" : vv); }});
    bool valid = false; int nsum = 0;
    { std::lock_guard<std::mutex> l(runner_.mutex()); valid = runner_.results().valid; nsum = int(runner_.results().summary.size()); }
    const bool running = runner_.running();
    const bool fresh = valid && resultsCaseVersion_ == caseVersion_;
    S.push_back({TRN("Run the simulation"), running ? 0 : fresh ? 1 : valid ? 2 : 0, running ? "running..." : fresh ? std::to_string(nsum) + " time steps, results up to date" : valid ? "results are from an older version of the inputs - run again (F5)" : "not run yet (F5)",
                 [this] { if (!runner_.running()) runSimulation(); }});
    S.push_back({TRN("Review results"), valid ? 1 : 0, valid ? "Results workspace (Ctrl+3): dashboard, compare runs, export figures" : "after the run", [this] { workspace_ = Workspace::Results; }});
    S.push_back({TRN("Save the project"), !path_.empty() && !dirty_ ? 1 : 0, path_.empty() ? "not saved yet (Ctrl+S)" : dirty_ ? "unsaved changes (Ctrl+S)" : fileOnly(path_), [this] { save(false); }});
    int done = 0; for (const auto& s : S) done += s.state == 1;
    ImGui::ProgressBar(float(done) / float(S.size()), ImVec2(-FLT_MIN, ImGui::GetFontSize() * 0.45f), "");
    for (size_t i = 0; i < S.size(); ++i) {
        const Step& s = S[i];
        ImGui::PushID(int(i));
        const char* ic = s.state == 1 ? ICON_CHECK : s.state == 2 ? ICON_TRIANGLE_ALERT : ICON_CIRCLE_DOT;
        const ImU32 col = s.state == 1 ? Pal::Ok : s.state == 2 ? Pal::Warn : Pal::Muted;
        ImGui::PushStyleColor(ImGuiCol_Text, col4(col));
        ImGui::TextUnformatted(ic);
        ImGui::PopStyleColor();
        ImGui::SameLine();
        char lb[96]; std::snprintf(lb, sizeof lb, "%zu. %s", i + 1, TR(s.label));
        if (ImGui::Selectable(lb, false)) { if (s.go) s.go(); }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort) && !s.why.empty()) {
            ImGui::BeginTooltip(); ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28); ImGui::TextUnformatted(s.why.c_str()); ImGui::PopTextWrapPos(); ImGui::EndTooltip();
        }
        ImGui::PopID();
    }
}

// ============================================================================ run history (Result Viewer compare list)
void App::rvAddRunHistory() {
    RVCase c;
    std::string text;
    bool ok = false;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        ok = R.valid && R.summary.size() >= 2 && R.engine == "ResSim";
        if (ok) {
            c.internal = true;
            c.summary = ecl::summaryFromRun(R.summary, R.wellNames, R.input.units);
            c.grid = R.grid;
            // report maps only for modest sizes (memory): cells x report steps <= 4e6
            if (R.grid && double(R.grid->numActive()) * double(R.snaps.size()) <= 4e6) c.snaps = R.snaps;
            c.input = std::make_shared<ressim::Case>(R.input);
            text = projectText(R.input);
        }
    }
    if (!ok) return;
    ++runCounter_;
    c.autoHist = true;
    c.changes = lastRunText_.empty() ? std::string("first run of this session") : textDiff(lastRunText_, text);
    lastRunText_ = std::move(text);
    char b[160];
    std::snprintf(b, sizeof b, "Run #%d  %s  %s", runCounter_, timeStr(std::time(nullptr)).c_str(), c.input ? c.input->title.substr(0, 32).c_str() : "");
    c.name = b;
    if (!variants_.empty() && c.input) casesLinkRun(projectText(*c.input), c.name);   // C40: link the run to an identical variant
    c.color = seriesColor(rv_.colorIdx++);
    rv_.cases.push_back(std::move(c));
    // keep the last N automatic entries (manually added runs / files stay)
    int nAuto = 0; for (const auto& x : rv_.cases) nAuto += x.autoHist;
    while (nAuto > std::max(1, settings_.runHistoryMax)) {
        for (size_t i = 0; i < rv_.cases.size(); ++i) if (rv_.cases[i].autoHist) { rv_.cases.erase(rv_.cases.begin() + long(i)); break; }
        --nAuto;
    }
    const int nc = int(rv_.cases.size());
    if (rv_.selected.empty()) rv_.selected = {"FOPR", "FPR", "FWCT", "FGOR"};
    if (nc >= 2) { rv_.ref = nc - 2; rv_.cmp = nc - 1; }
    rv_.ref = std::clamp(rv_.ref, 0, std::max(0, nc - 1)); rv_.cmp = std::clamp(rv_.cmp, 0, std::max(0, nc - 1));
    rv_.c3d = std::clamp(rv_.c3d, 0, std::max(0, nc - 1));
}

// ============================================================================ command palette
void App::commandPalette() {
    const ImGuiInputFlags G = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_P, G) || (workspace_ != Workspace::DeckEditor && ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_K, G))) {
        paletteOpen_ = !paletteOpen_; paletteBuf_[0] = 0; paletteSel_ = 0; paletteFocus_ = true;
    }
    if (!paletteOpen_) return;
    std::vector<PaletteCmd> C = buildCommands(std::strlen(paletteBuf_) >= 2);
    commandPaletteWindow(C);
}

// C41: the command registry, shared by the palette and the workflow recorder (ui_workflow.cpp)
std::vector<App::PaletteCmd> App::buildCommands(bool withKeywords) {
    std::vector<PaletteCmd> C;
    auto add = [&](std::string l, std::string h, std::function<void()> f, bool en = true) { C.push_back({std::move(l), std::move(h), std::move(f), en}); };
    const bool running = runner_.running();
    add("Run simulation", "F5", [this] { runSimulation(); }, !running);
    add("Stop the simulation", "Shift+F5", [this] { runner_.cancel(); }, running);
    add("Restart from a report step...", "Run", [this] { showRestart_ = true; }, !running);
    add("Generate report (Word / PowerPoint / HTML)...", "File", [this] { showReport_ = true; });
    add("Deck regression against OPM / ECLIPSE reference runs...", "Run", [this] { showRegress_ = true; });
    add("Compute initial state (Initialise)", "Run", [this] { computeInit(); section_ = Section::Init; });
    add("Validate inputs", "Run", [this] { const std::string& v = validation(); notify(v.empty() ? "All inputs are valid" : v); });
    add("Workspace: Model", "Ctrl+1", [this] { workspace_ = Workspace::Model; });
    add("Workspace: Deck Editor (writing mode)", "Ctrl+2", [this] { workspace_ = Workspace::DeckEditor; });
    add("Workspace: Results", "Ctrl+3", [this] { workspace_ = Workspace::Results; });
    add("Workspace: Studies (sensitivity / history matching)", "Ctrl+4", [this] { workspace_ = Workspace::Studies; });
    add("Workspace: Geomodel (ResGeo)", "Ctrl+5", [this] { workspace_ = Workspace::Geomodel; });
    add("New: demo 3-phase waterflood", "File", [this] { if (maybeSaveChanges()) newCase(0); });
    add("New: 1D Buckley-Leverett validation", "File", [this] { if (maybeSaveChanges()) newCase(1); });
    add("New: blank model", "Ctrl+N", [this] { if (maybeSaveChanges()) newCase(2); });
    add("New: model wizard (grid, fluid, well pattern, schedule)", "File", [this] { showWizard_ = true; });
    add("Open project / deck / results...", "Ctrl+O", [this] { openDialog(); });
    add("Save project", "Ctrl+S", [this] { save(false); });
    add("Save project as...", "Ctrl+Shift+S", [this] { save(true); });
    add("Import ECLIPSE deck (.DATA)...", "File", [this] { std::string p; if (maybeSaveChanges() && fileDialog(false, "ECLIPSE deck (*.DATA;*.data)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "DATA", p)) { dirty_ = false; importDeck(p); } });
    add("Export ECLIPSE deck (.DATA)...", "File", [this] { exportDeck(); });
    add("Open current model as deck text", "Deck Editor", [this] { edOpenModelDeck(); workspace_ = Workspace::DeckEditor; });
    add("Open ResSim results (.rsres)...", "File", [this] { std::string p; if (fileDialog(false, "ResSim results (*.rsres)\0*.rsres\0All files (*.*)\0*.*\0", "rsres", p)) openResultsFile(p); });
    add("Save results (.rsres)...", "File", [this] { saveResultsDialog(); }, !running);
    add("Export results as ECLIPSE files (EGRID / UNRST / UNSMRY)...", "File", [this] { exportResultsEclipseDialog(); }, !running);
    add("Open ECLIPSE / OPM Flow results...", "File", [this] { std::string p; if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.FEGRID;*.UNRST;*.FUNRST;*.UNSMRY;*.FUNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p)) rvOpen(p); });
    add("Compare: add the current run to the Result Viewer", "Results", [this] { rvAddRun(); });
    add("Result Viewer (ECLIPSE / OPM output, compare runs)", "View", [this] { showResultViewer_ = true; });
    add("Start page (recent files, examples)", "Help", [this] { showStart_ = true; });
    add("Notification history", "status bar", [this] { notesPopupReq_ = true; });
    add("User guide", "F1", [this] { openUserGuide(); });
    add("About ResSim Studio", "Help", [this] { showAbout_ = true; });
    add("Engine: ResSim (built-in)", "Run", [this] { settings_.engine = 0; saveSettings(); });
    add("Engine: OPM Flow (external)", "Run", [this] { settings_.engine = 1; saveSettings(); detectFlow(false); });
    add("Engine: ECLIPSE (external, licensed)", "Run", [this] { settings_.engine = 2; saveSettings(); detectEclipse(false); });
    add("Settings (engines, results, HPC, interface)...", "View", [this] { showOpmSettings_ = true; });
    add("Reset window layout", "View", [this] { resetLayout_ = true; });
    add("Reset 3D camera", "View", [this] { if (viewOk_) view_.resetCamera(); });
    for (float k : {0.8f, 0.9f, 1.0f, 1.15f, 1.3f}) { char l[48]; std::snprintf(l, sizeof l, "Interface size %.0f %%", k * 100); add(l, "View", [this, k] { settings_.uiScale = k; setUiScale(k); saveSettings(); }); }
    add("Deck import report", "View", [this] { showDeckTools_ = true; deckTab_ = 0; deckTabSet_ = true; });
    add("ECLIPSE keyword reference", "View", [this] { showDeckTools_ = true; deckTab_ = 1; deckTabSet_ = true; });
    add("Streamlines in the 3D viewer (toggle)", "Tools", [this] { sl_.show = !sl_.show; });
    add("Geomechanics: fault & fracture stability", "Tools", [this] { showStability_ = true; stabDock_ = true; });
    add("Geomechanics: seismicity-aware injection design", "Tools", [this] { showSeis_ = true; seisDock_ = true; });
    add("Geomechanics: fracture contact & slip (2-D)", "Tools", [this] { showEfem_ = true; efemDock_ = true; });
    add("Convert grid to corner-point (COORD/ZCORN)", "Tools", [this] { convertToCornerPoint(case_.grid); markDirty(true); notify("Grid converted to corner-point geometry"); }, !case_.grid.cornerPoint());
    static const char* kSec[] = {"Run specification", "Grid & rock", "Fluid (PVT / EOS / CO2)", "Rel. perm & Pc", "Initialization", "Aquifers", "EOR (polymer)", "Wells", "Schedule", "Numerics", "Advanced keywords"};
    for (int s = 0; s < int(Section::Count) && s < 11; ++s) add(std::string("Go to: ") + kSec[s], "Model", [this, s] { workspace_ = Workspace::Model; section_ = Section(s); });
    static const std::vector<ExampleInfo> ex = listExamples();
    for (const ExampleInfo& e : ex) { const std::string id = e.id; add("Example: " + e.title, toString(e.level), [this, id] { if (maybeSaveChanges()) newCase(3, id); }); }
    scanDeckFiles();
    for (const DeckFile& f : deckFiles_) { const std::string p = f.path; add("Example deck: " + f.name, f.folder, [this, p] { openFile(p); }); }
    for (const std::string& p : settings_.recent) add("Recent: " + fileOnly(p), p, [this, p] { openFile(p); });

    // ---- C40 feature registry: edit, themes, geomodel, wells, deep links into every tab, keywords
    add("Undo", "Ctrl+Z", [this] { undoApply(-1); }, canUndo());
    add("Redo", "Ctrl+Y", [this] { undoApply(+1); }, canRedo());
    {
        const char* const* tn = themeNames();
        for (int t = 0; t < 3; ++t) add(std::string("Theme: ") + tn[t], "View", [this, t] { settings_.theme = t; applyColors(settings_.theme, settings_.accent); saveSettings(); });
        const char* const* an = accentNames();
        for (int a = 0; a < 5; ++a) add(std::string("Accent colour: ") + an[a], "View", [this, a] { settings_.accent = a; applyColors(settings_.theme, settings_.accent); saveSettings(); });
        for (int b = 0; b < 5; ++b) add(std::string("3-D background: ") + background3dNames()[b], "View", [this, b] { settings_.bg3d = b; set3dBackground(b); saveSettings(); });
        add("Case manager (model variants, compare with base)", "View", [this] { showCases_ = true; ImGui::SetWindowFocus("###casemgr"); });
        add("Save the current model as a variant", "Model", [this] { casesAdd("Variant " + std::to_string(variants_.size() + 1), case_); showCases_ = true; notify("Saved variant " + variants_.back().name); });
    }
    add("Autosave: " + std::string(settings_.autosave ? "switch off" : "switch on"), "Settings", [this] { settings_.autosave = !settings_.autosave; saveSettings(); });
    static const char* kGeoSteps[] = {"Data", "Seismic", "Structure", "Grid", "Blocking", "Facies", "Petrophysics", "Volumetrics", "Upscaling", "Export"};
    for (int i = 0; i < 10; ++i) add(std::string("Geomodel: ") + kGeoSteps[i] + " step", "Geomodel", [this, i] { workspace_ = Workspace::Geomodel; geoStepReq_ = i; });
    static const char* kGeoViews[] = {"Map", "Seismic section", "Well correlation", "Grid cross-section", "3D scene", "Data analysis"};
    static const char* kGeoTabs[] = {"Log calculator", "Cross-plot", "Histogram", "Variogram", "Well tie"};   // C40 data-analysis tabs
    for (int i = 0; i < 5; ++i)
        add(std::string("Geomodel data analysis: ") + kGeoTabs[i], "Geomodel", [this, i] { workspace_ = Workspace::Geomodel; geoViewReq_ = 5; nav::request({kGeoTabs[i]}); });
    for (int i = 0; i < 6; ++i) add(std::string("Geomodel view: ") + kGeoViews[i], "Geomodel", [this, i] { workspace_ = Workspace::Geomodel; geoViewReq_ = i; });
    for (int w = 0; w < int(case_.wells.size()); ++w)
        add("Well: " + case_.wells[size_t(w)].name, "Model > Wells", [this, w] { workspace_ = Workspace::Model; section_ = Section::Wells; selWell_ = w; });
    {
        static const char* kWs[] = {"Model", "Deck Editor", "Results", "Studies", "Geomodel"};
        struct Link { int ws; int sec; std::vector<std::string> chain; };
        // destinations that exist before they have ever been drawn (the index below adds every tab seen so far)
        static const std::vector<Link> kLinks = {
            {2, -1, {"Field"}}, {2, -1, {"Wells"}}, {2, -1, {"Fields & sections"}}, {2, -1, {"Production analysis"}}, {2, -1, {"Custom plots"}},
            {2, -1, {"Convergence"}}, {2, -1, {"Multisegment wells"}}, {2, -1, {"Surface network"}}, {2, -1, {"UDQ & actions"}},
            {2, -1, {"Faults & seismicity"}}, {2, -1, {"Thermal"}}, {2, -1, {"Components"}},
            {3, -1, {"Method"}}, {3, -1, {"Profiles"}}, {3, -1, {"Cross-plot"}}, {3, -1, {"Tornado"}}, {3, -1, {"Monte Carlo"}}, {3, -1, {"Proxy"}},
            {3, -1, {"Pareto"}}, {3, -1, {"Match quality"}}, {3, -1, {"Adjoint HM"}}};
        auto chainText = [](const std::vector<std::string>& ch) { std::string s; for (size_t i = 0; i < ch.size(); ++i) s += (i ? " > " : "") + ch[i]; return s; };
        std::vector<std::string> seen;
        for (const nav::TabEntry& e : nav::index()) {
            std::vector<std::string> full = e.chain; full.push_back(e.label);
            const std::string t = chainText(full);
            if (std::find(seen.begin(), seen.end(), t) != seen.end()) continue;
            seen.push_back(t);
            const int ws = std::clamp(e.workspace, 0, 4), sec = e.section;
            std::string where = kWs[ws];
            if (ws == 0 && sec >= 0 && sec < 11) where += std::string(" > ") + kSec[sec];
            if (!e.window.empty() && e.window != "Inspector" && e.window.size() < 30) where += " (" + e.window + ")";
            add("Go to: " + t, where, [this, ws, sec, full] {
                workspace_ = Workspace(ws);
                if (ws == 0 && sec >= 0 && sec < int(Section::Count)) section_ = Section(sec);
                nav::request(full);
            });
        }
        for (const Link& l : kLinks) {
            const std::string t = chainText(l.chain);
            if (std::find(seen.begin(), seen.end(), t) != seen.end()) continue;
            const int ws = l.ws; const std::vector<std::string> full = l.chain;
            add("Go to: " + t, kWs[ws], [this, ws, full] { workspace_ = Workspace(ws); nav::request(full); });
        }
    }
    add("Workflow editor (record / replay command sequences)", "Tools", [this] { showWorkflow_ = true; ImGui::SetWindowFocus("###workflow"); });
    if (withKeywords) {
        static const std::vector<std::pair<std::string, std::string>> kw = keywordIndexList();
        for (const auto& k : kw) { const std::string name = k.first; add("Keyword: " + name, k.second, [this, name] { openKeywordReference(name); }); }
    }
    return C;
}

void App::commandPaletteWindow(std::vector<PaletteCmd>& C) {

    // filter + rank (C40: with an empty query the recently used commands come first)
    static std::vector<std::string> recentCmds;
    static bool recentLoaded = false;
    const std::string recentFile = plat::configDir() + "/palette_recent.txt";
    if (!recentLoaded) {
        recentLoaded = true;
        std::ifstream f(fs::u8path(recentFile));
        for (std::string l; std::getline(f, l);) if (!l.empty()) recentCmds.push_back(l);
    }
    std::vector<std::pair<int, int>> hits;   // score, index
    for (int i = 0; i < int(C.size()); ++i) {
        int s = fuzzyScore(C[size_t(i)].label, paletteBuf_);
        if (s < 0) continue;
        const auto r = std::find(recentCmds.begin(), recentCmds.end(), C[size_t(i)].label);
        if (r != recentCmds.end()) s += paletteBuf_[0] ? 6 : 100000 - 10 * int(r - recentCmds.begin());
        hits.push_back({s, i});
    }
    std::stable_sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    paletteSel_ = std::clamp(paletteSel_, 0, std::max(0, int(hits.size()) - 1));

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + fs * 3), ImGuiCond_Always, ImVec2(0.5f, 0));
    ImGui::SetNextWindowSize(ImVec2(std::min(vp->WorkSize.x * 0.9f, fs * 44), 0));
    ImGui::SetNextWindowFocus();
    std::function<void()> run;
    std::string runLabel;
    if (ImGui::Begin("##palette", &paletteOpen_, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize)) {
        if (paletteFocus_) { ImGui::SetKeyboardFocusHere(); paletteFocus_ = false; }
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::InputTextWithHint("##q", ICON_SEARCH "  Type a command, section, example or file...", paletteBuf_, sizeof paletteBuf_)) paletteSel_ = 0;
        if (ImGui::IsKeyPressed(ImGuiKey_DownArrow)) paletteSel_ = std::min(paletteSel_ + 1, std::max(0, int(hits.size()) - 1));
        if (ImGui::IsKeyPressed(ImGuiKey_UpArrow)) paletteSel_ = std::max(0, paletteSel_ - 1);
        if (ImGui::IsKeyPressed(ImGuiKey_Escape)) paletteOpen_ = false;
        const bool enter = ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter);
        ImGui::BeginChild("##pl", ImVec2(0, std::min(float(hits.size()), 14.0f) * ImGui::GetTextLineHeightWithSpacing() + 4));
        for (int r = 0; r < int(hits.size()); ++r) {
            const PaletteCmd& c = C[size_t(hits[size_t(r)].second)];
            ImGui::PushID(r);
            ImGui::BeginDisabled(!c.enabled);
            const bool sel = r == paletteSel_;
            if (ImGui::Selectable(c.label.c_str(), sel, ImGuiSelectableFlags_AllowOverlap) && c.enabled) { run = c.fn; runLabel = c.label; }
            if (sel && ImGui::IsWindowAppearing()) ImGui::SetScrollHereY();
            if (sel && (ImGui::IsKeyPressed(ImGuiKey_DownArrow) || ImGui::IsKeyPressed(ImGuiKey_UpArrow))) ImGui::SetScrollHereY();
            if (!c.hint.empty()) {
                const float w = ImGui::CalcTextSize(c.hint.c_str()).x;
                ImGui::SameLine(ImGui::GetContentRegionMax().x - std::min(w, fs * 16));
                ImGui::TextDisabled("%s", c.hint.size() > 40 ? ("..." + c.hint.substr(c.hint.size() - 37)).c_str() : c.hint.c_str());
            }
            ImGui::EndDisabled();
            ImGui::PopID();
        }
        if (hits.empty()) ImGui::TextDisabled("no match");
        ImGui::EndChild();
        if (enter && !hits.empty() && C[size_t(hits[size_t(paletteSel_)].second)].enabled) { run = C[size_t(hits[size_t(paletteSel_)].second)].fn; runLabel = C[size_t(hits[size_t(paletteSel_)].second)].label; }
        ImGui::TextDisabled("Enter run  |  Up / Down select  |  Esc close  |  Ctrl+Shift+P or Ctrl+K");
        if (!ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && !ImGui::IsWindowAppearing()) paletteOpen_ = false;
    }
    ImGui::End();
    if (run) {
        paletteOpen_ = false;
        if (runLabel.size()) {
            recentCmds.erase(std::remove(recentCmds.begin(), recentCmds.end(), runLabel), recentCmds.end());
            recentCmds.insert(recentCmds.begin(), runLabel);
            if (recentCmds.size() > 8) recentCmds.resize(8);
            std::ofstream f(fs::u8path(recentFile));
            for (const std::string& l : recentCmds) f << l << "\n";
        }
        wfRecord(runLabel);   // C41 workflow recorder
        run();
    }
}

// ============================================================================ start page
void App::startPage() {
    if (!showStart_) return;
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + vp->WorkSize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(std::min(vp->WorkSize.x * 0.92f, fs * 70), std::min(vp->WorkSize.y * 0.9f, fs * 40)), ImGuiCond_Appearing);
    if (!ImGui::Begin(ICON_HOUSE "  Start###startpage", &showStart_, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    std::function<void()> act;
    // ---- left column: get started + recent
    ImGui::BeginChild("##stl", ImVec2(fs * 22, -ImGui::GetFrameHeightWithSpacing()), ImGuiChildFlags_None);
    if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
    ImGui::TextUnformatted("ResSim Studio");
    if (g_fonts.bold) ImGui::PopFont();
    ImGui::TextDisabled("reservoir simulation, geomechanics, studies");
    ImGui::SeparatorText("Get started");
    const ImVec2 bw(-FLT_MIN, 0);
    if (ui::accentButton(ICON_WAND_SPARKLES "  New model (wizard)", bw)) act = [this] { showWizard_ = true; };   // C40
    if (ImGui::Button(ICON_FILE_PLUS "  New demo waterflood", bw)) act = [this] { if (maybeSaveChanges()) newCase(0); };
    if (ImGui::Button(ICON_FILE "  Blank model", bw)) act = [this] { if (maybeSaveChanges()) newCase(2); };
    if (ImGui::Button(ICON_FOLDER_OPEN "  Open project / deck / results...", bw)) act = [this] { openDialog(); };
    if (ImGui::Button(ICON_UPLOAD "  Import ECLIPSE / OPM deck...", bw))
        act = [this] { std::string p; if (maybeSaveChanges() && fileDialog(false, "ECLIPSE deck (*.DATA;*.data)\0*.DATA;*.data\0All files (*.*)\0*.*\0", "DATA", p)) { dirty_ = false; importDeck(p); } };
    if (ImGui::Button(ICON_BOOK_OPEN "  User guide (F1)", bw)) act = [this] { openUserGuide(); };
    ImGui::SeparatorText("Recent");
    if (settings_.recent.empty()) ImGui::TextDisabled("(nothing opened yet)");
    for (const std::string& p : settings_.recent) {
        std::error_code ec;
        const bool ok = fs::exists(fs::u8path(p), ec);
        ImGui::BeginDisabled(!ok);
        if (ImGui::Selectable(fileOnly(p).c_str())) { const std::string q = p; act = [this, q] { openFile(q); }; }
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("%s%s", p.c_str(), ok ? "" : "\n(not found)");
    }
    ImGui::SeparatorText("Tips");
    ImGui::PushTextWrapPos(0);
    ImGui::TextDisabled("Ctrl+Shift+P (or Ctrl+K): command palette - type what you want to do. Ctrl+1..4 switch Model / Deck Editor / Results / Studies. "
                        "F5 runs, Shift+F5 stops. Every finished run is kept for comparison in Results > Cases, vectors & compare with what changed.");
    ImGui::PopTextWrapPos();
    ImGui::EndChild();
    ImGui::SameLine();
    // ---- right column: examples
    ImGui::BeginChild("##str", ImVec2(0, -ImGui::GetFrameHeightWithSpacing()), ImGuiChildFlags_Borders);
    ImGui::SetNextItemWidth(fs * 16);
    ImGui::InputTextWithHint("##exf", ICON_SEARCH "  filter examples", startFilter_, sizeof startFilter_);
    const std::string q = startFilter_;
    if (ImGui::BeginTabBar("##sttabs")) {
        if (ui::tabItem("Built-in examples")) {
            static const std::vector<ExampleInfo> ex = listExamples();
            const ExampleLevel levels[3] = {ExampleLevel::Basic, ExampleLevel::Intermediate, ExampleLevel::Advanced};
            ImGui::BeginChild("##exl");
            for (ExampleLevel lv : levels) {
                ImGui::SeparatorText(toString(lv));
                for (const ExampleInfo& e : ex) {
                    if (e.level != lv || (fuzzyScore(e.title + " " + e.summary, q) < 0)) continue;
                    ImGui::PushID(e.id.c_str());
                    if (ImGui::Selectable(e.title.c_str(), false, 0, ImVec2(0, 0))) { const std::string id = e.id; act = [this, id] { if (maybeSaveChanges()) newCase(3, id); }; }
                    ImGui::Indent(fs);
                    ImGui::PushTextWrapPos(0);
                    ImGui::TextDisabled("%s  |  %s", e.summary.c_str(), e.grid.c_str());
                    ImGui::PopTextWrapPos();
                    ImGui::Unindent(fs);
                    ImGui::PopID();
                }
            }
            ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Example decks (examples\\, decks\\)")) {
            scanDeckFiles();
            ImGui::BeginChild("##dkl");
            if (deckFiles_.empty()) ImGui::TextDisabled("No examples\\ or decks\\ folder next to ResSimStudio.exe.");
            std::string folder;
            for (const DeckFile& f : deckFiles_) {
                if (fuzzyScore(f.name + " " + f.desc, q) < 0) continue;
                if (f.folder != folder) { folder = f.folder; ImGui::SeparatorText(folder.c_str()); }
                ImGui::PushID(f.path.c_str());
                if (ImGui::Selectable(f.name.c_str())) { const std::string p = f.path; act = [this, p] { openFile(p); }; }
                if (!f.desc.empty()) {
                    ImGui::Indent(fs);
                    ImGui::PushTextWrapPos(0);
                    const std::string d = f.desc.size() > 220 ? f.desc.substr(0, 217) + "..." : f.desc;
                    ImGui::TextDisabled("%s", d.c_str());
                    ImGui::PopTextWrapPos();
                    ImGui::Unindent(fs);
                    if (f.desc.size() > 220 && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                        ImGui::BeginTooltip(); ImGui::PushTextWrapPos(fs * 32); ImGui::TextUnformatted(f.desc.c_str()); ImGui::PopTextWrapPos(); ImGui::EndTooltip();
                    }
                }
                ImGui::PopID();
            }
            ImGui::EndChild();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::EndChild();
    if (ImGui::Checkbox("Show this page at start-up", &settings_.showStartPage)) saveSettings();
    ImGui::SameLine(ImGui::GetContentRegionMax().x - fs * 6);
    if (ImGui::Button("Close", ImVec2(fs * 6, 0))) showStart_ = false;
    ImGui::End();
    if (act) { showStart_ = false; act(); }
}
