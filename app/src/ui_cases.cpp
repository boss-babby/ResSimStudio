// C40 (G10) Case manager: named variants of the model kept with the project.
//   - save the current model as a variant, mark one variant as the base
//   - per variant: input differences against the base or the current model, run status (a run whose input equals the
//     variant is linked automatically), load / load & run / update / duplicate / delete / export
//   - compare the results of a variant with the base in the Result Viewer (reference = base run)
//   - stored next to the project: <project>.cases/index.txt + one .rsim per variant (written on Save, read on Open)
#include "app.h"
#include "icons.h"
#include "theme.h"
#include "widgets.h"
#include "ressim/io.h"
#include <algorithm>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace ressim;
namespace fs = std::filesystem;

std::string projectTextDiff(const std::string& a, const std::string& b, int maxLines);   // ui_usability.cpp

namespace {
std::string stamp(long long t) {
    if (t <= 0) return "-";
    std::time_t tt = std::time_t(t);
    char b[40];
    std::strftime(b, sizeof b, "%Y-%m-%d %H:%M", std::localtime(&tt));
    return b;
}
fs::path casesDir(const std::string& project) {
    fs::path p = fs::u8path(project);
    return p.parent_path() / fs::u8path(p.stem().u8string() + ".cases");
}
std::string clean(std::string s) {   // one-line text for the index file
    for (char& c : s) if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    return s;
}
int changedLines(const std::string& d) {   // "+ " / "- " lines of projectTextDiff, or the total it reports
    if (d == "no input change") return 0;
    const size_t p = d.find("... (");
    if (p != std::string::npos) { int r = 0, a = 0; if (std::sscanf(d.c_str() + p, "... (%d lines removed, %d added", &r, &a) == 2) return r + a; }
    int n = 0; std::istringstream is(d); std::string l; while (std::getline(is, l)) n += (l.size() > 1 && (l[0] == '+' || l[0] == '-') && l[1] == ' ');
    return n;
}
}  // namespace

int App::casesFindRun(const std::string& runName) const {
    if (runName.empty()) return -1;
    for (int i = 0; i < int(rv_.cases.size()); ++i) if (rv_.cases[size_t(i)].name == runName) return i;
    return -1;
}

void App::casesAdd(const std::string& name, const Case& c) {
    CaseVariant v;
    v.name = name;
    for (int k = 2; std::any_of(variants_.begin(), variants_.end(), [&](const CaseVariant& x) { return x.name == v.name; }); ++k) v.name = name + " (" + std::to_string(k) + ")";
    v.c = std::make_shared<const Case>(c);
    v.text = projectText(c);
    v.created = (long long)std::time(nullptr);
    v.base = variants_.empty();
    variants_.push_back(std::move(v));
    casesSel_ = int(variants_.size()) - 1;
    dirty_ = true;   // the variants are saved with the project
}

void App::casesLinkRun(const std::string& inputText, std::string& runName) {
    std::vector<CaseVariant*> hit;
    for (auto& v : variants_) if (v.text == inputText) { hit.push_back(&v); runName += "  [" + v.name + "]"; }
    for (CaseVariant* v : hit) v->run = runName;
}

void App::casesWrite(const std::string& project) {
    std::error_code ec;
    const fs::path d = casesDir(project);
    if (variants_.empty()) { if (fs::exists(d / "index.txt", ec)) fs::remove_all(d, ec); return; }
    fs::create_directories(d, ec);
    // remove old variant files, then write the current set
    for (const auto& e : fs::directory_iterator(d, ec)) if (e.path().extension() == ".rsim") fs::remove(e.path(), ec);
    std::ofstream f(d / "index.txt", std::ios::binary);
    f << "# ResSim case manager index (C40): name<TAB>base<TAB>created<TAB>file<TAB>note\n";
    for (size_t i = 0; i < variants_.size(); ++i) {
        const CaseVariant& v = variants_[i];
        const std::string file = "variant_" + std::to_string(i + 1) + ".rsim";
        try { saveProject(*v.c, (d / file).u8string()); } catch (const std::exception& e) { notify(std::string("Case manager: ") + e.what()); continue; }
        f << clean(v.name) << '\t' << int(v.base) << '\t' << v.created << '\t' << file << '\t' << clean(v.note) << '\n';
    }
}

void App::casesRead(const std::string& project) {
    variants_.clear(); casesSel_ = -1;
    std::error_code ec;
    const fs::path d = casesDir(project);
    std::ifstream f(d / "index.txt", std::ios::binary);
    if (!f) return;
    std::string line;
    int bad = 0;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> t; std::string x; std::istringstream is(line);
        while (std::getline(is, x, '\t')) t.push_back(x);
        if (t.size() < 4) { ++bad; continue; }
        try {
            CaseVariant v;
            v.name = t[0]; v.base = t[1] == "1"; v.created = std::atoll(t[2].c_str()); v.note = t.size() > 4 ? t[4] : "";
            v.c = std::make_shared<const Case>(loadProject((d / fs::u8path(t[3])).u8string()));
            v.text = projectText(*v.c);
            variants_.push_back(std::move(v));
        } catch (...) { ++bad; }
    }
    if (!variants_.empty() && std::none_of(variants_.begin(), variants_.end(), [](const CaseVariant& v) { return v.base; })) variants_[0].base = true;
    if (!variants_.empty()) runner_.addLog("Case manager: " + std::to_string(variants_.size()) + " variant(s) loaded from " + d.u8string());
    if (bad) notify("Case manager: " + std::to_string(bad) + " variant(s) could not be read");
}

void App::panelCases() {
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 52, ImGui::GetFontSize() * 30), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin(ICON_LAYERS "  Case manager###casemgr", &showCases_)) { ImGui::End(); return; }
    const float fs = ImGui::GetFontSize();
    static char nameBuf[96] = "";
    if (!nameBuf[0]) std::snprintf(nameBuf, sizeof nameBuf, "Variant %zu", variants_.size() + 1);
    ImGui::SetNextItemWidth(fs * 14);
    ImGui::InputText("##vname", nameBuf, sizeof nameBuf);
    ImGui::SameLine();
    if (ui::accentButton(ICON_PLUS "  Save the current model as a variant")) {
        casesAdd(nameBuf[0] ? nameBuf : "Variant", case_);
        nameBuf[0] = 0;
        notify("Saved variant " + variants_.back().name + (variants_.back().base ? " (base)" : ""));
    }
    ui::hint("Variants are complete copies of the model (grid, fluids, wells, schedule, study parameters). Load one to edit or run it; "
             "runs whose input equals a variant are linked to it. The variants are saved with the project (<project>.cases folder).");
    if (variants_.empty()) {
        ui::emptyState(ICON_LAYERS, "No variants yet", "Save the current model as the base case, change it, and save the changes as further variants.");
        ImGui::End();
        return;
    }
    int baseIdx = 0;
    for (int i = 0; i < int(variants_.size()); ++i) if (variants_[size_t(i)].base) baseIdx = i;
    const CaseVariant& B = variants_[size_t(baseIdx)];
    static std::string curText;   // the current model's text, rebuilt only when the model changes
    static unsigned curVer = ~0u;
    if (curVer != unsigned(caseVersion_) || curText.empty()) { curText = projectText(case_); curVer = unsigned(caseVersion_); }
    casesSel_ = std::clamp(casesSel_, 0, int(variants_.size()) - 1);

    std::function<void()> deferred;
    const ImGuiTableFlags tf = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp;
    if (ImGui::BeginTable("##variants", 5, tf, ImVec2(0, std::min(fs * 1.9f * float(variants_.size() + 1) + 4, fs * 12)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Variant", ImGuiTableColumnFlags_WidthStretch, 2.2f);
        ImGui::TableSetupColumn("Changes vs base", ImGuiTableColumnFlags_WidthStretch, 1.1f);
        ImGui::TableSetupColumn("Current model", ImGuiTableColumnFlags_WidthStretch, 1.0f);
        ImGui::TableSetupColumn("Run", ImGuiTableColumnFlags_WidthStretch, 1.6f);
        ImGui::TableSetupColumn("Created", ImGuiTableColumnFlags_WidthStretch, 1.1f);
        ImGui::TableHeadersRow();
        for (int i = 0; i < int(variants_.size()); ++i) {
            const CaseVariant& v = variants_[size_t(i)];
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            std::string lab = std::string(v.base ? ICON_STAR "  " : "     ") + v.name;
            if (ImGui::Selectable(lab.c_str(), casesSel_ == i, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowDoubleClick)) {
                casesSel_ = i;
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) deferred = [this, i] { casesLoad(i, false); };
            }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) ImGui::SetTooltip("%s%s", v.note.empty() ? "" : (v.note + "\n").c_str(), "double-click: load into the model");
            ImGui::TableSetColumnIndex(1);
            if (v.base) ImGui::TextDisabled("base");
            else {
                CaseVariant& w = variants_[size_t(i)];
                if (w.diffBase != B.c.get() || w.diffSelf != w.c.get()) { w.nDiff = changedLines(projectTextDiff(B.text, w.text, 4)); w.diffBase = B.c.get(); w.diffSelf = w.c.get(); }
                ImGui::Text(w.nDiff ? "%d lines" : "identical", w.nDiff);
            }
            ImGui::TableSetColumnIndex(2);
            if (v.text == curText) { ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Ok)); ImGui::TextUnformatted(ICON_CHECK " loaded"); ImGui::PopStyleColor(); }
            else ImGui::TextDisabled("-");
            ImGui::TableSetColumnIndex(3);
            const int r = casesFindRun(v.run);
            if (r >= 0) { ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Ok)); ImGui::TextUnformatted(ICON_CHECK); ImGui::PopStyleColor(); ImGui::SameLine(); ImGui::TextUnformatted(rv_.cases[size_t(r)].name.substr(0, 28).c_str()); }
            else if (!v.run.empty()) ImGui::TextDisabled(ICON_HISTORY " no longer in the run history");
            else ImGui::TextDisabled("not run");
            ImGui::TableSetColumnIndex(4);
            ImGui::TextDisabled("%s", stamp(v.created).c_str());
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    CaseVariant& S = variants_[size_t(casesSel_)];
    ImGui::Separator();
    ImGui::SetNextItemWidth(fs * 14);
    if (ui::inputString("##rename", &S.name)) dirty_ = true;
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ui::inputStringHint("##note", "note (what this variant tests)", &S.note)) dirty_ = true;
    const int i = casesSel_;
    if (ImGui::Button(ICON_FOLDER_OPEN "  Load into the model")) deferred = [this, i] { casesLoad(i, false); };
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Replaces the current model (Ctrl+Z restores it)");
    ImGui::SameLine();
    ImGui::BeginDisabled(runner_.running());
    if (ImGui::Button(ICON_PLAY "  Load & run")) deferred = [this, i] { casesLoad(i, true); };
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button(ICON_SAVE "  Update from the model")) { S.c = std::make_shared<const Case>(case_); S.text = curText; S.run.clear(); dirty_ = true; notify("Variant " + S.name + " updated from the current model"); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Overwrite this variant with the current model");
    ImGui::SameLine();
    if (ImGui::Button(ICON_COPY "  Duplicate")) deferred = [this, i] { const Case c = *variants_[size_t(i)].c; casesAdd(variants_[size_t(i)].name + " copy", c); variants_.back().base = false; };
    ImGui::SameLine();
    ImGui::BeginDisabled(S.base);
    if (ImGui::Button(ICON_STAR "  Set as base")) { for (auto& v : variants_) v.base = false; S.base = true; dirty_ = true; }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button(ICON_DOWNLOAD "  Export...")) {
        std::string p = S.name + ".DATA";
        for (char& ch : p) if (std::strchr("\\/:*?\"<>|", ch)) ch = '_';
        if (fileDialog(true, "ECLIPSE deck (*.DATA)\0*.DATA\0ResSim project (*.rsim)\0*.rsim\0", "DATA", p)) {
            try {
                if (p.size() > 5 && (p.substr(p.size() - 5) == ".rsim" || p.substr(p.size() - 5) == ".RSIM")) saveProject(*S.c, p); else exportEclipseDeck(*S.c, p);
                notify("Exported " + p);
            } catch (const std::exception& e) { notify(std::string("Export failed: ") + e.what()); }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(ICON_TRASH_2 "  Delete")) deferred = [this, i] {
        const bool wasBase = variants_[size_t(i)].base;
        variants_.erase(variants_.begin() + i);
        if (wasBase && !variants_.empty()) variants_[0].base = true;
        casesSel_ = std::min(i, int(variants_.size()) - 1); dirty_ = true;
    };
    {   // results comparison with the base
        const int rb = casesFindRun(B.run), rs = casesFindRun(S.run);
        ImGui::BeginDisabled(rb < 0 || rs < 0 || rb == rs);
        if (ImGui::Button(ICON_CHART_LINE "  Compare results with the base")) {
            rv_.ref = rb; rv_.cmp = rs; rv_.c3d = rs;
            for (auto& c : rv_.cases) c.show = false;
            rv_.cases[size_t(rb)].show = rv_.cases[size_t(rs)].show = true;
            workspace_ = Workspace::Results; showResultViewer_ = true;
        }
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
            ImGui::SetTooltip(rb < 0 || rs < 0 ? "Run the base and this variant first (Load & run)" : "Result Viewer: base = reference, this variant = compared case");
    }
    static int diffAgainst = 0;
    ImGui::Spacing();
    ImGui::TextUnformatted("Input differences of");
    ImGui::SameLine(); ImGui::TextColored(col4(Pal::Accent), "%s", S.name.c_str());
    ImGui::SameLine(); ImGui::TextUnformatted("against");
    ImGui::SameLine(); ImGui::RadioButton("the base", &diffAgainst, 0);
    ImGui::SameLine(); ImGui::RadioButton("the current model", &diffAgainst, 1);
    static std::string d;
    static const void* dk[3] = {nullptr, nullptr, nullptr};
    static int dkMode = -1; static unsigned dkVer = ~0u;
    if (dk[0] != S.c.get() || dk[1] != B.c.get() || dkMode != diffAgainst || (diffAgainst == 1 && dkVer != curVer)) {
        d = diffAgainst == 0 ? (S.base ? std::string("(this is the base)") : projectTextDiff(B.text, S.text, 60)) : projectTextDiff(S.text, curText, 60);
        dk[0] = S.c.get(); dk[1] = B.c.get(); dkMode = diffAgainst; dkVer = curVer;
    }
    ImGui::BeginChild("##diff", ImVec2(0, 0), ImGuiChildFlags_Borders);
    ImGui::PushFont(g_fonts.mono ? g_fonts.mono : ImGui::GetFont(), 0.0f);
    std::istringstream is(d); std::string l;
    while (std::getline(is, l)) {
        const bool add = l.rfind("+ ", 0) == 0, rem = l.rfind("- ", 0) == 0;
        if (add || rem) ImGui::PushStyleColor(ImGuiCol_Text, col4(add ? Pal::Ok : Pal::Err));
        ImGui::TextUnformatted(l.c_str());
        if (add || rem) ImGui::PopStyleColor();
    }
    ImGui::PopFont();
    ImGui::EndChild();
    ImGui::End();
    if (deferred) deferred();
}

void App::casesLoad(int i, bool run) {
    if (i < 0 || i >= int(variants_.size())) return;
    if (runner_.running()) { notify("A run is in progress - stop it first"); return; }
    case_ = *variants_[size_t(i)].c;
    syncBuildersFromCase();
    markDirty(true);
    workspace_ = Workspace::Model;
    notify("Loaded variant " + variants_[size_t(i)].name + " (Ctrl+Z restores the previous model)");
    if (run) runSimulation();
}
