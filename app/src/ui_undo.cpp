// C39 undo / redo of model edits and autosave / crash recovery.
//
// Undo: every editor already reports a change through markDirty() (caseVersion_ bump). At the end of each frame, once no
// widget is active (a slider drag or a typed number counts as ONE edit gesture), the whole case is snapshotted. Restoring
// a snapshot goes through the same path as opening a project, so every panel, the grid preview and the builders stay
// consistent - one mechanism covers all ~60 editors. The history depth adapts to the model size (memory budget).
//
// Autosave: while there are unsaved changes, a copy of the case is written in the background to the settings folder
// every settings_.autosaveSec seconds. A normal save or exit removes it; after a crash the next start offers to restore.
#include "app.h"
#include "icons.h"
#include "platform.h"
#include "theme.h"
#include "widgets.h"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace ressim;

namespace {
// rough memory of one case snapshot: cell arrays dominate (grid geometry, properties, regions, initial state)
size_t snapshotBytes(const Case& c) {
    const double cells = double(std::max(1, c.grid.nx)) * std::max(1, c.grid.ny) * std::max(1, c.grid.nz);
    return size_t(cells * 8.0 * 14.0) + (size_t(2) << 20);
}
}  // namespace

// ============================================================================ undo / redo
void App::undoReset() {
    ui::modBaselineReset();   // C41 Inspector "show only modified" baseline = the freshly loaded model
    undo_.stack.clear();
    undo_.stack.push_back(std::make_shared<const Case>(case_));
    undo_.cur = 0;
    undo_.version = caseVersion_;
    const size_t budget = size_t(640) << 20;   // ~640 MB of history at most
    undo_.maxLevels = std::clamp<size_t>(budget / snapshotBytes(case_), 4, 100);
}

void App::undoTrack() {
    if (undo_.cur < 0) { undoReset(); return; }
    if (caseVersion_ == undo_.version) return;
    if (ImGui::IsAnyItemActive()) return;   // gesture still in progress: one step when it ends
    undo_.stack.resize(size_t(undo_.cur) + 1);   // a new edit drops the redo branch
    undo_.stack.push_back(std::make_shared<const Case>(case_));
    if (undo_.stack.size() > undo_.maxLevels) undo_.stack.erase(undo_.stack.begin(), undo_.stack.begin() + long(undo_.stack.size() - undo_.maxLevels));
    undo_.cur = int(undo_.stack.size()) - 1;
    undo_.version = caseVersion_;
}

void App::undoApply(int dir) {
    const int to = undo_.cur + dir;
    if (to < 0 || to >= int(undo_.stack.size())) return;
    if (runner_.running() || study_.running) { notify("Undo is paused while a simulation or study is running"); return; }
    undo_.cur = to;
    case_ = *undo_.stack[size_t(to)];
    syncBuildersFromCase();
    dirty_ = true;
    ++caseVersion_;
    undo_.version = caseVersion_;
    init_ = InitPreview{};
    gridDirty_ = true;
    rebuildPreview();
    const int nw = int(case_.wells.size()), ns = int(case_.schedule.size());
    selWell_ = std::clamp(selWell_, 0, std::max(0, nw - 1));
    resultsWell_ = std::clamp(resultsWell_, 0, std::max(0, nw - 1));
    selStep_ = std::clamp(selStep_, 0, std::max(0, ns - 1));
    char b[96];
    std::snprintf(b, sizeof(b), "%s  (step %d of %d)", dir < 0 ? "Undo" : "Redo", to, int(undo_.stack.size()) - 1);
    notify(b);
}

// ============================================================================ autosave / recovery
std::string App::recoveryPath() const { return (fs::u8path(plat::configDir()) / "recovery.rsim").u8string(); }

void App::autosaveTick() {
    if (!settings_.autosave || !dirty_ || autosaveBusy_) return;
    const double now = ImGui::GetTime();
    if (autosavedVersion_ == caseVersion_) return;
    if (autosaveAt_ == 0) { autosaveAt_ = now + std::max(15, settings_.autosaveSec); return; }   // first change: start the clock
    if (now < autosaveAt_ || ImGui::IsAnyItemActive()) return;
    autosaveAt_ = now + std::max(15, settings_.autosaveSec);
    autosavedVersion_ = caseVersion_;
    if (autosaveThr_.joinable()) autosaveThr_.join();
    auto snap = (undo_.cur >= 0 && undo_.version == caseVersion_) ? undo_.stack[size_t(undo_.cur)] : std::make_shared<const Case>(case_);
    const std::string file = recoveryPath(), orig = path_, title = case_.title;
    autosaveBusy_ = true;
    autosaveThr_ = std::thread([this, snap, file, orig, title]() {
        try {
            const std::string tmp = file + ".tmp";
            saveProject(*snap, tmp);
            std::error_code ec;
            fs::rename(fs::u8path(tmp), fs::u8path(file), ec);
            if (ec) { fs::copy_file(fs::u8path(tmp), fs::u8path(file), fs::copy_options::overwrite_existing, ec); fs::remove(fs::u8path(tmp), ec); }
            std::ofstream m(fs::u8path(file + ".txt"), std::ios::binary);
            const std::time_t t = std::time(nullptr);
            char when[64];
            std::strftime(when, sizeof(when), "%Y-%m-%d %H:%M", std::localtime(&t));
            m << "original=" << orig << "\ntitle=" << title << "\ntime=" << when << "\n";
        } catch (...) {}
        autosaveBusy_ = false;
    });
}

void App::recoveryClear() {
    if (autosaveThr_.joinable()) autosaveThr_.join();
    std::error_code ec;
    fs::remove(fs::u8path(recoveryPath()), ec);
    fs::remove(fs::u8path(recoveryPath() + ".txt"), ec);
    autosavedVersion_ = caseVersion_;
    autosaveAt_ = 0;
}

void App::recoveryCheck() {
    std::error_code ec;
    if (!fs::exists(fs::u8path(recoveryPath()), ec)) return;
    recoveryOrig_.clear();
    std::string title, when;
    std::ifstream m(fs::u8path(recoveryPath() + ".txt"), std::ios::binary);
    for (std::string line; std::getline(m, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("original=", 0) == 0) recoveryOrig_ = line.substr(9);
        else if (line.rfind("title=", 0) == 0) title = line.substr(6);
        else if (line.rfind("time=", 0) == 0) when = line.substr(5);
    }
    recoveryInfo_ = "\"" + (title.empty() ? std::string("untitled case") : title) + "\"";
    if (!when.empty()) recoveryInfo_ += ", autosaved " + when;
    if (!recoveryOrig_.empty()) recoveryInfo_ += "\nOriginal file: " + recoveryOrig_;
    recoveryPrompt_ = true;
    showStart_ = false;   // one dialog at a time: the start page returns on Discard
}

void App::recoveryWindow() {
    if (!recoveryPrompt_) return;
    if (!ImGui::IsPopupOpen("Recover unsaved work")) ImGui::OpenPopup("Recover unsaved work");
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 32, 0), ImGuiCond_Appearing);
    if (ImGui::BeginPopupModal("Recover unsaved work", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings)) {
        ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT "  ResSim Studio did not close normally last time.");
        ImGui::Spacing();
        ImGui::TextWrapped("An autosaved copy of the case you were editing is available:");
        ImGui::Indent();
        ImGui::TextUnformatted(recoveryInfo_.c_str());
        ImGui::Unindent();
        ImGui::Spacing();
        if (ui::accentButton(ICON_ROTATE_CCW "  Restore")) {
            try {
                Case c = loadProject(recoveryPath());
                runner_.cancel(); runner_.join(); runner_.clear();
                case_ = std::move(c);
                const bool origProject = recoveryOrig_.size() > 5 && recoveryOrig_.substr(recoveryOrig_.size() - 5) == ".rsim";
                path_ = origProject ? recoveryOrig_ : std::string();
                init_ = InitPreview{};
                syncBuildersFromCase();
                ++caseVersion_;
                resultsCaseVersion_ = 0;
                gridDirty_ = true;
                rebuildPreview();
                selWell_ = selStep_ = resultsWell_ = 0;
                if (viewOk_) view_.resetCamera();
                undoReset();
                dirty_ = true;                  // not saved yet: Save writes it back
                showStart_ = false;
                notify("Recovered the autosaved case - save it to keep it");
            } catch (const std::exception& e) {
                plat::msgBox(e.what(), "Recovery failed", MB_OK | MB_ICONERROR);
            }
            recoveryPrompt_ = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Discard")) {
            showStart_ = settings_.showStartPage;
            recoveryClear();
            recoveryPrompt_ = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        ui::hint("Autosave: Settings > Interface");
        ImGui::EndPopup();
    }
}
