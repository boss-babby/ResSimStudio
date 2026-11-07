// C41 (Petrel parity 6): workflow recorder / editor.
//
// Records the commands run from the command palette (Ctrl+Shift+P) into an editable step list, lets the user add
// steps from the same registry, "wait for the simulation", notes and a "for each model variant" loop, saves /
// loads .rswf files and replays them one step per frame (engine/include/ressim/workflow.h holds the format and
// the runner state machine; tests/c41_tests.inc covers them).
#include "app.h"
#include "icons.h"
#include "theme.h"
#include "widgets.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;
using ressim::wf::Kind;
using ressim::wf::Step;

namespace {
bool containsCI(const std::string& hay, const char* needle) {
    std::string h = hay, n = needle;
    for (char& c : h) c = char(std::tolower(static_cast<unsigned char>(c)));
    for (char& c : n) c = char(std::tolower(static_cast<unsigned char>(c)));
    return h.find(n) != std::string::npos;
}
}  // namespace

void App::wfRecord(const std::string& label) {
    if (!wfRecording_ || label.empty() || wfRun_.active()) return;
    if (label.rfind("Workflow editor", 0) == 0) return;
    wfSteps_.push_back({Kind::Command, label});
    if (label == "Run simulation") wfSteps_.push_back({Kind::Wait, ""});   // replay waits for the run it started
}

void App::wfTick() {
    if (!wfRun_.active()) return;
    auto exec = [this](const std::string& label) {
        std::vector<PaletteCmd> C = buildCommands(label.rfind("Keyword: ", 0) == 0);
        for (PaletteCmd& c : C)
            if (c.label == label) {
                if (!c.enabled) return false;
                c.fn();
                return true;
            }
        return false;
    };
    auto busy = [this] { return runner_.running(); };
    auto nVar = [this] { return int(variants_.size()); };
    auto loadVar = [this](int i) {
        if (i < 0 || i >= int(variants_.size()) || runner_.running()) return false;
        casesLoad(i, false);
        return true;
    };
    const bool was = wfRun_.active();
    wfRun_.tick(exec, busy, nVar, loadVar);
    if (was && !wfRun_.active()) {
        if (wfRun_.error.empty()) notify("Workflow finished (" + std::to_string(wfRun_.executed) + " command(s))");
        else notify("Workflow stopped: " + wfRun_.error);
    }
}

void App::panelWorkflow() {
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 34, ImGui::GetFontSize() * 30), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin(ICON_LIST_CHECKS "  Workflow editor###workflow", &showWorkflow_)) { ImGui::End(); return; }
    const bool running = wfRun_.active();

    // ---- toolbar
    if (!wfRecording_) {
        if (ImGui::Button(ICON_CIRCLE_DOT "  Record")) { wfRecording_ = true; notify("Workflow recording: commands run from the palette (Ctrl+Shift+P) are added"); }
    } else {
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.75f, 0.2f, 0.2f, 1));
        if (ImGui::Button(ICON_CIRCLE_STOP "  Stop recording")) wfRecording_ = false;
        ImGui::PopStyleColor();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(running || wfSteps_.empty());
    if (ui::accentButton(ICON_PLAY "  Run")) {
        wfRecording_ = false;
        if (!wfRun_.start(wfSteps_)) notify("Workflow not valid: " + wfRun_.error);
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(!running);
    if (ImGui::Button(ICON_X "  Abort")) wfRun_.stop("aborted by the user");
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(running);
    if (ui::secondaryButton(ICON_FOLDER_OPEN "  Open...")) {
        std::string p;
        if (fileDialog(false, "ResSim workflow (*.rswf)\0*.rswf\0All files (*.*)\0*.*\0", "rswf", p)) {
            std::ifstream f(fs::u8path(p), std::ios::binary);
            std::stringstream ss; ss << f.rdbuf();
            std::vector<Step> st; std::string err;
            if (ressim::wf::parse(ss.str(), st, err)) { wfSteps_ = st; wfPath_ = p; notify("Workflow loaded: " + std::to_string(st.size()) + " step(s)"); }
            else notify("Could not read the workflow: " + err);
        }
    }
    ImGui::SameLine();
    if (ui::secondaryButton(ICON_SAVE "  Save...")) {
        std::string p = wfPath_;
        if (fileDialog(true, "ResSim workflow (*.rswf)\0*.rswf\0All files (*.*)\0*.*\0", "rswf", p)) {
            std::ofstream f(fs::u8path(p), std::ios::binary);
            f << ressim::wf::serialize(wfSteps_);
            if (f) { wfPath_ = p; notify("Workflow saved: " + p); } else notify("Could not write " + p);
        }
    }
    ImGui::SameLine();
    if (ui::dangerButton(ICON_TRASH_2 "  Clear")) wfSteps_.clear();
    ImGui::EndDisabled();
    if (wfRecording_) { ImGui::SameLine(); ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.35f, 1), ICON_CIRCLE_DOT " REC"); }

    const std::string verr = ressim::wf::validate(wfSteps_);
    if (!verr.empty()) ImGui::TextColored(ImVec4(0.95f, 0.6f, 0.2f, 1), ICON_TRIANGLE_ALERT "  %s", verr.c_str());
    if (running) ImGui::TextColored(ImVec4(0.4f, 0.8f, 0.4f, 1), ICON_PLAY "  running step %d of %d%s", wfRun_.pc + 1, int(wfRun_.steps.size()),
                                    wfRun_.loopStart >= 0 ? ("  (variant " + std::to_string(wfRun_.variant + 1) + ")").c_str() : "");
    else if (!wfRun_.error.empty()) ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.35f, 1), "Last run: %s", wfRun_.error.c_str());
    ImGui::Separator();

    // ---- step list
    int del = -1, up = -1, down = -1;
    const float bottom = ImGui::GetFrameHeightWithSpacing() * 3.4f;
    if (ImGui::BeginChild("##wfsteps", ImVec2(0, -bottom), ImGuiChildFlags_Borders)) {
        if (wfSteps_.empty()) ImGui::TextDisabled("No steps. Press Record and run commands from the palette (Ctrl+Shift+P), or add steps below.");
        int indent = 0;
        for (int i = 0; i < int(wfSteps_.size()); ++i) {
            Step& s = wfSteps_[size_t(i)];
            ImGui::PushID(i);
            if (s.kind == Kind::End) indent = 0;
            const bool cur = running && wfRun_.pc == i;
            ImGui::TextColored(cur ? ImVec4(0.4f, 0.85f, 0.4f, 1) : ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled), "%s%2d", cur ? ">" : " ", i + 1);
            ImGui::SameLine();
            if (indent) { ImGui::Dummy(ImVec2(ImGui::GetFontSize(), 1)); ImGui::SameLine(); }
            switch (s.kind) {
            case Kind::Command: ImGui::TextUnformatted(ICON_TERMINAL); ImGui::SameLine(); ImGui::TextUnformatted(s.arg.c_str()); break;
            case Kind::Wait: ImGui::TextDisabled(ICON_HISTORY "  wait until the simulation has finished"); break;
            case Kind::Note: {
                char b[256]; std::snprintf(b, sizeof b, "%s", s.arg.c_str());
                ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 6);
                ImGui::BeginDisabled(running);
                if (ImGui::InputTextWithHint("##note", "note", b, sizeof b)) s.arg = b;
                ImGui::EndDisabled();
                break;
            }
            case Kind::ForVariants: ImGui::TextColored(ImVec4(0.55f, 0.7f, 1, 1), ICON_LAYERS "  for each model variant (case manager: %d)", int(variants_.size())); break;
            case Kind::End: ImGui::TextColored(ImVec4(0.55f, 0.7f, 1, 1), "end of loop"); break;
            }
            if (s.kind == Kind::ForVariants) indent = 1;
            if (!running) {
                ImGui::SameLine(ImGui::GetContentRegionMax().x - ImGui::GetFontSize() * 4.6f);
                if (ImGui::SmallButton("^")) up = i;
                ImGui::SameLine();
                if (ImGui::SmallButton("v")) down = i;
                ImGui::SameLine();
                if (ImGui::SmallButton(ICON_X)) del = i;
            }
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
    if (del >= 0) wfSteps_.erase(wfSteps_.begin() + del);
    if (up > 0) std::swap(wfSteps_[size_t(up)], wfSteps_[size_t(up - 1)]);
    if (down >= 0 && down + 1 < int(wfSteps_.size())) std::swap(wfSteps_[size_t(down)], wfSteps_[size_t(down + 1)]);

    // ---- add steps
    ImGui::BeginDisabled(running);
    static char filt[96] = {0};
    static int pick = -1;
    static std::vector<std::string> labels;
    if (labels.empty() || ImGui::IsWindowAppearing()) {
        labels.clear();
        for (const PaletteCmd& c : buildCommands(false))
            if (c.label.rfind("Workflow editor", 0) != 0 && c.label.rfind("Recent: ", 0) != 0) labels.push_back(c.label);
    }
    ImGui::SetNextItemWidth(-ImGui::GetFontSize() * 7);
    const char* prev = pick >= 0 && pick < int(labels.size()) ? labels[size_t(pick)].c_str() : "choose a command...";
    if (ImGui::BeginCombo("##wfcmd", prev, ImGuiComboFlags_HeightLarge)) {
        if (ImGui::IsWindowAppearing()) ImGui::SetKeyboardFocusHere();
        ImGui::InputTextWithHint("##wff", ICON_SEARCH " filter", filt, sizeof filt);
        for (int i = 0; i < int(labels.size()); ++i) {
            if (filt[0] && !containsCI(labels[size_t(i)], filt)) continue;
            if (ImGui::Selectable(labels[size_t(i)].c_str(), pick == i)) pick = i;
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(pick < 0 || pick >= int(labels.size()));
    if (ImGui::Button(ICON_PLUS "  Add")) {
        wfSteps_.push_back({Kind::Command, labels[size_t(pick)]});
        if (labels[size_t(pick)] == "Run simulation") wfSteps_.push_back({Kind::Wait, ""});
    }
    ImGui::EndDisabled();
    if (ImGui::Button(ICON_HISTORY "  Wait for run")) wfSteps_.push_back({Kind::Wait, ""});
    ImGui::SameLine();
    if (ImGui::Button(ICON_PENCIL "  Note")) wfSteps_.push_back({Kind::Note, ""});
    ImGui::SameLine();
    if (ImGui::Button(ICON_LAYERS "  Loop over variants")) { wfSteps_.push_back({Kind::ForVariants, ""}); wfSteps_.push_back({Kind::End, ""}); }
    ImGui::EndDisabled();
    ImGui::TextDisabled("Steps replay palette commands by name; one step per frame. .rswf files are plain text.");
    ImGui::End();
}
