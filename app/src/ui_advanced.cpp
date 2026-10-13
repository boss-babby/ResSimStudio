// Advanced keywords editor: every ECLIPSE / OPM Flow keyword the ResSim engine does not model natively (groups,
// VFP tables, multisegment wells, ACTIONX / UDQ, tracers, extra summary requests, ...) is kept verbatim in
// Case::extraKeywords. This page lists, edits, adds, re-orders and checks them; the deck exporter writes them in
// place, so OPM Flow runs the full case.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "ressim/deck_parser.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>

using namespace ressim;

namespace {

const char* kSecNames[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};

std::string upperTrim(const std::string& s) {
    std::string r;
    for (char ch : s) if (!std::isspace((unsigned char)ch)) r += char(std::toupper((unsigned char)ch));
    return r;
}

// Parses "SECTION / NAME / body" with the generic deck reader and returns its messages (empty = clean)
std::vector<std::string> checkKeyword(const ExtraKeyword& e) {
    std::vector<std::string> out;
    if (e.name.empty()) { out.push_back("Error: keyword name is empty"); return out; }
    if (!deck::findKeyword(e.name)) out.push_back("Warning: " + e.name + " is not in the OPM Flow / ECLIPSE keyword index (typing error?)");
    std::string text = std::string(kSecNames[int(e.section)]) + "\n" + e.name + "\n" + e.body;
    if (!text.empty() && text.back() != '\n') text += '\n';
    try {
        const deck::RawDeck D = deck::parseDeckText(text, "", "keyword");
        for (auto& m : D.messages) {
            if (m.level == deck::Level::Info) continue;
            out.push_back(std::string(m.level == deck::Level::Error ? "Error" : "Warning") + " (line " + std::to_string(std::max(0, m.line - 2)) + "): " + m.text);
        }
        int n = 0; for (auto& k : D.keywords) if (k.name == e.name) ++n;
        if (n == 0) out.push_back("Error: the keyword was not recognised");
        else if (D.keywords.size() > 2) out.push_back("Warning: the text contains more than one keyword - each one is written as is");
    } catch (const std::exception& ex) { out.push_back(std::string("Error: ") + ex.what()); }
    return out;
}

}  // namespace

void App::editAdvanced() {
    auto& X = case_.extraKeywords;
    const float fs = ImGui::GetFontSize();
    ui::hint("Run with OPM Flow to simulate them. SCHEDULE keywords are written at their report step.");

    // ---- toolbar: filter + add
    ImGui::SetNextItemWidth(fs * 10);
    ui::inputString("##flt", advFilter_, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase);
    if (advFilter_.empty()) { ImGui::SameLine(fs * 0.6f); ImGui::TextDisabled("filter"); }
    ImGui::SameLine(fs * 11.5f);
    ImGui::SetNextItemWidth(fs * 7);
    ui::inputString("##newkw", advNewName_, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
    ImGui::Combo("##newsec", &advNewSection_, kSecNames, 8);
    ImGui::SameLine();
    const std::string nn = upperTrim(advNewName_);
    ImGui::BeginDisabled(nn.empty());
    if (ImGui::Button("+ Add keyword")) {
        ExtraKeyword e; e.name = nn; e.section = DeckSection(advNewSection_); e.note = "added in ResSim Studio";
        const deck::KeywordSpec* sp = deck::findKeyword(nn);
        if (sp && sp->kind == deck::KwSize::Slash) e.body = " /\n/\n";
        else if (sp && sp->kind == deck::KwSize::None) e.body = "";
        else e.body = " /\n";
        if (e.section == DeckSection::Schedule) e.step = 0;
        // insert after the last keyword of the same section (keeps the deck order readable)
        int at = int(X.size());
        for (int i = int(X.size()) - 1; i >= 0; --i) if (int(X[size_t(i)].section) <= advNewSection_) { at = i + 1; break; }
        if (X.empty() || int(X.front().section) > advNewSection_) at = 0;
        X.insert(X.begin() + at, e);
        advSel_ = at; advNewName_.clear(); markDirty();
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Type a keyword name (e.g. GRUPTREE, GCONPROD, VFPPROD, WELSEGS, ACTIONX, WECON, TRACER),\n"
                                                  "pick its section and add it. The keyword manual is shown below the editor.");

    // ---- list
    double day = 0; std::vector<double> stepDay(case_.schedule.size() + 1, 0.0);
    for (size_t s = 0; s < case_.schedule.size(); ++s) { stepDay[s] = day; day += case_.schedule[s].length; }
    stepDay.back() = day;
    // list height: the rows it needs, at most ~40 % of the panel (the editor goes below)
    const float listMax = std::max(fs * 10, ImGui::GetContentRegionAvail().y * 0.4f);
    const float listH = std::min(listMax, ImGui::GetTextLineHeightWithSpacing() * 1.25f * (float(X.size()) + 1.5f) + fs);
    const bool wide = ImGui::GetContentRegionAvail().x > fs * 30;   // room for the Note column
    int moveFrom = -1, moveTo = -1, del = -1;
    if (ImGui::BeginTable("xkw", wide ? 6 : 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit,
                          ImVec2(0, std::max(listH, fs * 5)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("#", 0, fs * 2);
        ImGui::TableSetupColumn("Section", 0, fs * 5);
        ImGui::TableSetupColumn("Keyword", 0, fs * 5.5f);
        ImGui::TableSetupColumn("When", wide ? 0 : ImGuiTableColumnFlags_WidthStretch, fs * 6.5f);
        if (wide) ImGui::TableSetupColumn("Note", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("", 0, fs * 4.2f);
        ImGui::TableHeadersRow();
        for (int i = 0; i < int(X.size()); ++i) {
            const ExtraKeyword& e = X[size_t(i)];
            if (!advFilter_.empty() && e.name.find(upperTrim(advFilter_)) == std::string::npos &&
                std::string(kSecNames[int(e.section)]).find(upperTrim(advFilter_)) == std::string::npos) continue;
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            char b[16]; std::snprintf(b, sizeof(b), "%d", i + 1);
            if (ImGui::Selectable(b, advSel_ == i, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) advSel_ = i;
            if (!wide && !e.note.empty() && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", e.note.c_str());
            ImGui::TableNextColumn(); ImGui::TextDisabled("%s", kSecNames[int(e.section)]);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(e.name.c_str());
            ImGui::TableNextColumn();
            if (e.section == DeckSection::Schedule) {
                const int s = std::clamp(e.step, 0, int(case_.schedule.size()));
                ImGui::Text("step %d, day %.0f", s, stepDay[size_t(s)]);
            }
            if (wide) { ImGui::TableNextColumn(); ImGui::TextDisabled("%s", e.note.c_str()); }
            ImGui::TableNextColumn();
            if (ImGui::SmallButton("\xE2\x86\x91") && i > 0) { moveFrom = i; moveTo = i - 1; }
            ImGui::SameLine(0, 2);
            if (ImGui::SmallButton("\xE2\x86\x93") && i + 1 < int(X.size())) { moveFrom = i; moveTo = i + 1; }
            ImGui::SameLine(0, 2);
            if (ImGui::SmallButton("\xC3\x97")) del = i;
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (moveFrom >= 0) { std::swap(X[size_t(moveFrom)], X[size_t(moveTo)]); if (advSel_ == moveFrom) advSel_ = moveTo; else if (advSel_ == moveTo) advSel_ = moveFrom; markDirty(); }
    if (del >= 0) { X.erase(X.begin() + del); if (advSel_ >= int(X.size())) advSel_ = int(X.size()) - 1; markDirty(); }
    if (X.empty()) ImGui::TextDisabled("No extra keywords. Imported decks keep every keyword ResSim does not model here; or add one above.");

    // ---- editor of the selected keyword
    if (advSel_ < 0 || advSel_ >= int(X.size())) return;
    ExtraKeyword& e = X[size_t(advSel_)];
    ui::header(("Keyword " + std::to_string(advSel_ + 1) + ": " + e.name).c_str());
    bool ch = false;
    if (ui::beginProps("xkwp")) {
        ui::label("Name");
        if (ui::inputString("##nm", e.name, ImVec2(0, 0), false, ImGuiInputTextFlags_CharsUppercase)) { e.name = upperTrim(e.name); ch = true; }
        ui::label("Section");
        int sec = int(e.section);
        if (ImGui::Combo("##sec", &sec, kSecNames, 8)) { e.section = DeckSection(sec); ch = true; }
        if (e.section == DeckSection::Schedule) {
            ui::label("Report step", "The keyword is written after the well controls of this report step, before its TSTEP. "
                                     "Step 0 = start of the schedule; the last value = after the final step.");
            int st = e.step;
            if (ImGui::SliderInt("##step", &st, 0, int(case_.schedule.size()), "%d")) { e.step = st; ch = true; }
            const int s = std::clamp(e.step, 0, int(case_.schedule.size()));
            ui::label("Time"); ImGui::Text("day %.1f", stepDay[size_t(s)]);
            ui::label("Order in the step", "Before: written ahead of the step's WCONPROD / WCONINJE records (needed for VFP tables, "
                                           "group trees, UDQs they refer to). After: written behind them (e.g. records that modify them).");
            int bc = e.beforeControls ? 0 : 1;
            if (ImGui::Combo("##bc", &bc, "Before the well controls\0After the well controls\0")) { e.beforeControls = bc == 0; ch = true; }
        }
        std::string sn; const deck::Support sup = deck::ressimSupport(e.name, &sn);
        ui::label("ResSim engine");
        ImGui::TextWrapped("%s%s%s", deck::toString(sup), sn.empty() ? "" : " - ", sn.c_str());
        if (sup == deck::Support::Full)
            ImGui::TextColored(col4(Pal::Warn), "ResSim edits this keyword natively - the native data is also written (possible duplicate).");
        ui::endProps();
    }
    const float bodyH = std::max(fs * 10, ImGui::GetContentRegionAvail().y - fs * 5);
    // schema form view (parameter names, units, defaults, option lists from the OPM keyword index)
    namespace di = ressim::deckintel;
    const bool hugeBody = e.body.size() > 256 * 1024;
    if (!hugeBody) {
        if (ImGui::RadioButton("Form", advFormView_)) advFormView_ = true;
        ImGui::SameLine();
        if (ImGui::RadioButton("Text", !advFormView_)) advFormView_ = false;
        ImGui::SameLine();
    }
    bool formShown = false;
    if (advFormView_ && !hugeBody) {
        if (advFormKey_ != advSel_ || advFormSrc_ != e.body || advForm_.keyword != e.name) {
            advForm_ = di::parseForm(e.name, di::splitLines(e.body));
            advFormKey_ = advSel_;
            advFormSrc_ = e.body;
        }
        if (advForm_.table || advForm_.editable()) {
            formShown = true;
            ImGui::TextDisabled("Records of %s - hover a column for its meaning, units and default", e.name.c_str());
            ImGui::BeginChild("##advform", ImVec2(0, bodyH), ImGuiChildFlags_Borders, ImGuiWindowFlags_HorizontalScrollbar);
            if (deckFormRecords(advForm_, case_.units == UnitSystem::Metric ? "METRIC" : "FIELD")) {
                e.body = di::joinLines(di::writeForm(advForm_)) + "\n";
                advFormSrc_ = e.body;
                ch = true;
            }
            ImGui::EndChild();
        } else {
            ImGui::TextDisabled("(no record schema for a form - text view)");
        }
    }
    if (!formShown) ImGui::TextDisabled("Records (ECLIPSE syntax: items separated by spaces, 1* = default, each record ends with /)");
    constexpr size_t kBig = 256 * 1024;      // larger bodies (big VFP tables) are shown read-only unless asked
    if (formShown) {
    } else if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
    if (formShown) {
    } else if (e.body.size() > kBig && advBigEdit_ != advSel_) {
        size_t lines = 0; for (char ch2 : e.body) lines += ch2 == '\n';
        size_t cut = 0; for (int n = 0; n < 300 && cut != std::string::npos; ++n) cut = e.body.find('\n', cut + 1);
        std::string head = e.body.substr(0, cut == std::string::npos ? e.body.size() : cut);
        head += "\n...";
        ImGui::InputTextMultiline("##bodyro", head.data(), head.size() + 1, ImVec2(-FLT_MIN, bodyH), ImGuiInputTextFlags_ReadOnly);
        if (g_fonts.mono) ImGui::PopFont();
        ImGui::TextColored(col4(Pal::Warn), "Large keyword: %.1f MB, %zu lines - the first 300 lines are shown.", double(e.body.size()) / 1048576.0, lines);
        ImGui::SameLine();
        if (ImGui::SmallButton("Edit all")) advBigEdit_ = advSel_;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Edit the full text (can be slow for very large tables)");
    } else {
        if (ui::inputString("##body", e.body, ImVec2(-FLT_MIN, bodyH), true, ImGuiInputTextFlags_AllowTabInput)) ch = true;
        if (g_fonts.mono) ImGui::PopFont();
    }
    if (ch) markDirty();
    // syntax check (re-parsed only when the keyword changes)
    if (advCheckSel_ != advSel_ || advCheckBody_ != e.body || advCheckName_ != e.name + char('0' + int(e.section))) {
        advCheckSel_ = advSel_; advCheckBody_ = e.body; advCheckName_ = e.name + char('0' + int(e.section));
        advCheckMsgs_ = checkKeyword(e);
    }
    if (advCheckMsgs_.empty()) ImGui::TextColored(col4(Pal::Ok), "\xE2\x9C\x93 Syntax OK (record structure as OPM Flow reads it)");
    for (auto& m : advCheckMsgs_) ImGui::TextColored(col4(m.rfind("Error", 0) == 0 ? Pal::Err : Pal::Warn), "%s", m.c_str());
    if (ImGui::CollapsingHeader("Keyword manual")) keywordDocView(e.name);
}
