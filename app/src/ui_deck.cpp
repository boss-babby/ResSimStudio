// Deck tools window: ECLIPSE deck import report + keyword reference (OPM Flow manual index).
#include "app.h"
#include "resources.h"
#include "widgets.h"
#include "imgui_internal.h"
#include "ressim/io.h"
#include "ressim/deck_parser.h"
#include <windows.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace ressim;

namespace {

struct KwParam { std::string index, name, manual, def, units, options, desc; };
struct KwDoc {
    std::string sections, alias, summary, example;
    std::vector<KwParam> params;
};

std::vector<std::string> splitTabs(const std::string& s) {
    std::vector<std::string> out; size_t a = 0;
    for (;;) { size_t b = s.find('\t', a); out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a)); if (b == std::string::npos) break; a = b + 1; }
    return out;
}
std::string unescapeNl(const std::string& s) {
    std::string o; o.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) { if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == 'n') { o += '\n'; ++i; } else o += s[i]; }
    return o;
}

// keyword_docs.txt lives next to the executable; loaded once on first use.
const std::map<std::string, KwDoc>& keywordDocs(bool* found = nullptr) {
    static std::map<std::string, KwDoc> docs;
    static bool loaded = false, ok = false;
    if (!loaded) {
        loaded = true;
        char buf[MAX_PATH] = {0};
        GetModuleFileNameA(nullptr, buf, MAX_PATH);
        std::string dir(buf); const size_t sl = dir.find_last_of("\\/"); dir = sl == std::string::npos ? std::string() : dir.substr(0, sl + 1);
        std::ifstream ff(dir + "keyword_docs.txt");     // optional override; the manual is embedded in the exe
        std::istringstream fe;
        if (!ff) fe.str(std::string(appResource(RES_KEYWORD_DOCS)));
        std::istream& f = ff ? static_cast<std::istream&>(ff) : static_cast<std::istream&>(fe);
        ok = bool(ff) || !fe.str().empty();
        std::string line; KwDoc* cur = nullptr;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.size() < 2 || line[0] == '#' || line[1] != '\t') continue;
            const auto t = splitTabs(line);
            if (line[0] == 'K' && t.size() >= 2) {
                KwDoc& d = docs[t[1]];
                if (t.size() > 2) d.sections = t[2];
                if (t.size() > 3) d.alias = t[3];
                if (t.size() > 4) d.summary = t[4];
                cur = &d;
            } else if (line[0] == 'P' && cur && t.size() >= 3) {
                KwParam p; p.index = t[1]; p.name = t[2];
                if (t.size() > 3) p.manual = t[3];
                if (t.size() > 4) p.def = t[4];
                if (t.size() > 5) p.units = t[5] == "|" ? std::string() : t[5];
                if (t.size() > 6) p.options = t[6];
                if (t.size() > 7) p.desc = t[7];
                cur->params.push_back(p);
            } else if (line[0] == 'E' && cur && t.size() >= 2) {
                if (!cur->example.empty()) cur->example += "\n\n";
                cur->example += unescapeNl(t[1]);
            }
        }
    }
    if (found) *found = ok;
    return docs;
}

ImVec4 supportColor(deck::Support s) {
    switch (s) {
    case deck::Support::Full: return ImVec4(0.45f, 0.85f, 0.50f, 1);
    case deck::Support::Partial: return ImVec4(0.95f, 0.80f, 0.35f, 1);
    case deck::Support::Ignored: return ImVec4(0.60f, 0.62f, 0.68f, 1);
    case deck::Support::Unsupported: return ImVec4(0.95f, 0.45f, 0.40f, 1);
    default: return ImVec4(0.75f, 0.55f, 0.95f, 1);
    }
}

std::string sectionsOf(uint16_t bits) {
    std::string s;
    for (int b = 0; b < 8; ++b) if (bits & (1u << b)) { if (!s.empty()) s += ' '; s += deck::sectionName(b); }
    return s;
}

} // namespace

void App::openKeywordReference(const std::string& kw) {
    showDeckTools_ = true;
    deckTab_ = 1; deckTabSet_ = true;
    kwSel_ = kw;
    std::snprintf(kwSearch_, sizeof(kwSearch_), "%s", kw.c_str());
    ImGui::SetWindowFocus("Deck tools");
}

void App::keywordDocView(const std::string& kw) {
    std::string note;
    const deck::Support sup = deck::ressimSupport(kw, &note);
    const deck::KeywordSpec* spec = deck::findKeyword(kw);
    ImGui::PushFont(nullptr);
    ImGui::TextColored(ImVec4(0.55f, 0.75f, 1.0f, 1), "%s", kw.c_str());
    ImGui::PopFont();
    ImGui::SameLine();
    ImGui::TextColored(supportColor(sup), "[ResSim: %s]", deck::toString(sup));
    if (!note.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", note.c_str()); }
    if (spec) ImGui::TextDisabled("Sections: %s", sectionsOf(spec->sections).c_str());
    bool haveDocs = false;
    const auto& docs = keywordDocs(&haveDocs);
    auto it = docs.find(kw);
    if (it == docs.end()) {
        ImGui::Separator();
        ImGui::TextWrapped(haveDocs ? "No manual entry for this keyword." : "The embedded keyword manual could not be loaded - documentation is unavailable.");
        return;
    }
    const KwDoc& d = it->second;
    if (!d.alias.empty() && d.alias != kw) ImGui::TextDisabled("Family: %s", d.alias.c_str());
    ImGui::Separator();
    ImGui::TextWrapped("%s", d.summary.c_str());
    if (!d.params.empty() && ImGui::CollapsingHeader("Items", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::BeginTable("kwp", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 2);
            ImGui::TableSetupColumn("Item", ImGuiTableColumnFlags_WidthStretch, 1.2f);
            ImGui::TableSetupColumn("Default / units", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Description", ImGuiTableColumnFlags_WidthStretch, 3.5f);
            ImGui::TableHeadersRow();
            for (const auto& p : d.params) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(p.index.c_str());
                ImGui::TableNextColumn(); ImGui::TextWrapped("%s", (p.manual.empty() ? p.name : p.manual).c_str());
                ImGui::TableNextColumn();
                std::string du = p.def; if (!p.units.empty()) du += (du.empty() ? "" : "  ") + p.units;
                ImGui::TextWrapped("%s", du.c_str());
                ImGui::TableNextColumn();
                std::string desc = p.desc; if (!p.options.empty()) desc += "  Options: " + p.options;
                ImGui::TextWrapped("%s", desc.c_str());
            }
            ImGui::EndTable();
        }
    }
    if (!d.example.empty() && ImGui::CollapsingHeader("Example")) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.80f, 0.85f, 0.75f, 1));
        ImGui::TextWrapped("%s", d.example.c_str());
        ImGui::PopStyleColor();
    }
}

void App::panelDeckTools() {
    if (!showDeckTools_) return;
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 60, ImGui::GetFontSize() * 34), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Deck tools", &showDeckTools_)) { ImGui::End(); return; }
    if (ImGui::BeginTabBar("decktabs")) {
        // ------------------------------------------------ import report
        if (ImGui::BeginTabItem("Import report", nullptr, deckTab_ == 0 && deckTabSet_ ? ImGuiTabItemFlags_SetSelected : 0)) {
            if (deckTab_ == 0) deckTabSet_ = false;
            if (deckReportPath_.empty()) {
                ImGui::TextDisabled("No deck imported yet. Use File > Import ECLIPSE deck...");
            } else {
                int nInfo = 0, nWarn = 0, nErr = 0;
                for (const auto& m : deckReport_) (m.level >= 2 ? nErr : m.level == 1 ? nWarn : nInfo)++;
                ImGui::Text("%s", deckReportPath_.c_str());
                ImGui::TextColored(ImVec4(0.95f, 0.45f, 0.40f, 1), "%d errors", nErr); ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.95f, 0.80f, 0.35f, 1), "%d warnings", nWarn); ImGui::SameLine();
                ImGui::TextDisabled("%d notes", nInfo);
                ImGui::Checkbox("Notes", &deckShowInfo_); ImGui::SameLine();
                ImGui::Checkbox("Warnings", &deckShowWarn_); ImGui::SameLine();
                static ImGuiTextFilter filter;
                filter.Draw("Filter", ImGui::GetFontSize() * 14); ImGui::SameLine();
                if (ImGui::SmallButton("Copy all")) {
                    std::ostringstream os;
                    for (const auto& m : deckReport_)
                        os << (m.level >= 2 ? "ERROR" : m.level == 1 ? "WARN " : "NOTE ") << "  " << m.keyword << "  " << m.file << ':' << m.line << "  " << m.text << "\n";
                    ImGui::SetClipboardText(os.str().c_str());
                }
                ImGui::TextDisabled("Double-click a row to open the keyword reference.");
                if (ImGui::BeginTable("deckrep", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp)) {
                    ImGui::TableSetupScrollFreeze(0, 1);
                    ImGui::TableSetupColumn("Level", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 4.5f);
                    ImGui::TableSetupColumn("Keyword", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 7);
                    ImGui::TableSetupColumn("Location", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 11);
                    ImGui::TableSetupColumn("Message", ImGuiTableColumnFlags_WidthStretch);
                    ImGui::TableHeadersRow();
                    for (size_t i = 0; i < deckReport_.size(); ++i) {
                        const DeckMessage& m = deckReport_[i];
                        if (m.level == 0 && !deckShowInfo_) continue;
                        if (m.level == 1 && !deckShowWarn_) continue;
                        const std::string loc = m.file.empty() ? std::string() : m.file + (m.line > 0 ? ":" + std::to_string(m.line) : std::string());
                        if (filter.IsActive() && !filter.PassFilter((m.keyword + " " + loc + " " + m.text).c_str())) continue;
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::PushID(int(i));
                        const ImVec4 col = m.level >= 2 ? ImVec4(0.95f, 0.45f, 0.40f, 1) : m.level == 1 ? ImVec4(0.95f, 0.80f, 0.35f, 1) : ImVec4(0.60f, 0.62f, 0.68f, 1);
                        ImGui::PushStyleColor(ImGuiCol_Text, col);
                        if (ImGui::Selectable(m.level >= 2 ? "error" : m.level == 1 ? "warning" : "note", deckSel_ == int(i),
                                              ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowDoubleClick)) {
                            deckSel_ = int(i);
                            if (ImGui::IsMouseDoubleClicked(0) && !m.keyword.empty()) openKeywordReference(m.keyword);
                        }
                        ImGui::PopStyleColor();
                        ImGui::PopID();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(m.keyword.c_str());
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%s", loc.c_str());
                        ImGui::TableNextColumn(); ImGui::TextWrapped("%s", m.text.c_str());
                    }
                    ImGui::EndTable();
                }
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ keyword reference
        if (ImGui::BeginTabItem("Keyword reference", nullptr, deckTab_ == 1 && deckTabSet_ ? ImGuiTabItemFlags_SetSelected : 0)) {
            if (deckTab_ == 1) deckTabSet_ = false;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
            ImGui::InputTextWithHint("##kws", "search keyword", kwSearch_, sizeof(kwSearch_), ImGuiInputTextFlags_CharsUppercase);
            ImGui::SameLine();
            ImGui::Checkbox("Supported only", &kwSupportedOnly_);
            const float listW = ImGui::GetFontSize() * 14;
            ImGui::BeginChild("kwlist", ImVec2(listW, 0), ImGuiChildFlags_Borders);
            const std::string q(kwSearch_);
            int shown = 0;
            for (const auto& k : deck::keywordTable()) {
                const std::string name(k.name);
                if (!q.empty() && name.find(q) == std::string::npos) continue;
                const deck::Support s = deck::ressimSupport(name);
                if (kwSupportedOnly_ && (s == deck::Support::Unsupported || s == deck::Support::Unknown)) continue;
                if (++shown > 1500) { ImGui::TextDisabled("... refine the search"); break; }
                ImGui::PushStyleColor(ImGuiCol_Text, supportColor(s));
                if (ImGui::Selectable(k.name, kwSel_ == name)) kwSel_ = name;
                ImGui::PopStyleColor();
            }
            ImGui::EndChild();
            ImGui::SameLine();
            ImGui::BeginChild("kwdoc", ImVec2(0, 0), ImGuiChildFlags_Borders);
            if (kwSel_.empty()) {
                ImGui::TextDisabled("Select a keyword. Colours: ");
                ImGui::TextColored(supportColor(deck::Support::Full), "Full"); ImGui::SameLine();
                ImGui::TextColored(supportColor(deck::Support::Partial), "Partial"); ImGui::SameLine();
                ImGui::TextColored(supportColor(deck::Support::Ignored), "Ignored (harmless)"); ImGui::SameLine();
                ImGui::TextColored(supportColor(deck::Support::Unsupported), "Unsupported");
            } else keywordDocView(kwSel_);
            ImGui::EndChild();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}
