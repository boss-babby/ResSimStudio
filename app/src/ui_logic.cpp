// C31 Model-mode editor of the run-time logic of a deck (Schedule section):
//   UDQ blocks   - DEFINE / ASSIGN / UPDATE / UNITS records as a form (live expression check) or as deck text
//   ACTIONX      - header (name, max runs, min interval), condition rows (lhs op rhs AND|OR), body keywords, templates
//   UDA          - a UDQ as the target of a well control (WCONPROD / WCONINJE / WELTARG item)
// Everything is stored as Case::extraKeywords (UDQ / ACTIONX, verbatim deck text at a report step) and WellControl::uda,
// so the deck written back is the ECLIPSE / OPM text and the ResSim engine evaluates the same logic.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "ressim/udq.h"
#include <algorithm>
#include <cstdio>
#include <set>
#include <sstream>

using namespace ressim;

namespace {
std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::string upper(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::string stripComment(const std::string& l) { const size_t p = l.find("--"); return p == std::string::npos ? l : l.substr(0, p); }
std::vector<std::string> tokens(const std::string& s) {   // whitespace split, quoted strings kept together (with quotes)
    std::vector<std::string> t; std::string cur; bool q = false;
    for (char c : s) {
        if (c == '\'') { q = !q; cur += c; continue; }
        if (!q && (c == ' ' || c == '\t')) { if (!cur.empty()) t.push_back(cur), cur.clear(); continue; }
        cur += c;
    }
    if (!cur.empty()) t.push_back(cur);
    return t;
}

// ---------------------------------------------------------------- UDQ form
struct UdqRow { int kind = 0; std::string name, text, sel; };   // kind 0 DEFINE 1 ASSIGN 2 UPDATE 3 UNITS
const char* kUdqKinds[] = {"DEFINE", "ASSIGN", "UPDATE", "UNITS"};
bool hasComments(const std::string& body) { std::istringstream is(body); std::string l; while (std::getline(is, l)) if (l.find("--") != std::string::npos) return true; return false; }
std::vector<UdqRow> parseUdqRows(const std::string& body) {
    std::vector<UdqRow> rows; std::string rec;
    std::istringstream is(body); std::string l;
    while (std::getline(is, l)) {
        l = trim(stripComment(l));
        if (l.empty()) continue;
        if (rec.empty() && l == "/") break;
        rec += (rec.empty() ? "" : " ") + l;
        if (rec.back() != '/') continue;
        rec = trim(rec.substr(0, rec.size() - 1));
        std::vector<std::string> t = tokens(rec);
        rec.clear();
        if (t.size() < 2) continue;
        UdqRow r; const std::string k = upper(t[0]);
        r.kind = k == "DEFINE" ? 0 : k == "ASSIGN" ? 1 : k == "UPDATE" ? 2 : k == "UNITS" ? 3 : -1;
        if (r.kind < 0) continue;
        r.name = upper(t[1]);
        std::string rest; for (size_t i = 2; i < t.size(); ++i) rest += (i > 2 ? " " : "") + t[i];
        if (r.kind == 1 && t.size() >= 4) { r.sel = t[2]; r.text = t.back(); }
        else if (r.kind == 3) { r.text = rest; if (r.text.size() >= 2 && r.text.front() == '\'') r.text = r.text.substr(1, r.text.size() - 2); }
        else r.text = rest;
        rows.push_back(r);
    }
    return rows;
}
std::string writeUdqRows(const std::vector<UdqRow>& rows) {
    std::string b;
    for (const auto& r : rows) {
        if (r.name.empty()) continue;
        b += std::string(" ") + kUdqKinds[r.kind] + " " + r.name;
        if (r.kind == 1 && !trim(r.sel).empty()) b += " " + trim(r.sel);
        if (r.kind == 3) b += " '" + r.text + "'"; else b += " " + trim(r.text);
        b += " /\n";
    }
    return b + "/\n";
}

// ---------------------------------------------------------------- ACTIONX form
struct CondRow { std::string lhs, op = ">", rhs; int conn = 0; };   // conn 0 none 1 AND 2 OR
const char* kOps[] = {">", "<", ">=", "<=", "=", "!="};
const char* kConn[] = {"(end)", "AND", "OR"};
struct ActForm { std::string name = "ACT1", body; int runs = 1; double wait = 0; std::vector<CondRow> conds; bool structured = true; };
ActForm parseActForm(const std::string& raw) {
    ActForm f; std::istringstream is(raw); std::string l; int stage = 0;
    while (std::getline(is, l)) {
        const std::string c = trim(stripComment(l));
        if (stage == 2) { if (upper(c).rfind("ENDACTIO", 0) == 0) break; f.body += l + "\n"; continue; }
        if (c.empty()) continue;
        if (stage == 0) {
            std::vector<std::string> t = tokens(c.back() == '/' ? c.substr(0, c.size() - 1) : c);
            if (!t.empty()) f.name = t[0];
            if (t.size() > 1 && t[1] != "1*") f.runs = std::atoi(t[1].c_str());
            if (t.size() > 2 && t[2] != "1*") f.wait = std::atof(t[2].c_str());
            stage = 1; continue;
        }
        if (c == "/") { stage = 2; continue; }
        std::string s = c; if (s.back() == '/') s = trim(s.substr(0, s.size() - 1));
        if (s.find('(') != std::string::npos || s.find(')') != std::string::npos) f.structured = false;
        std::vector<std::string> t = tokens(s);
        CondRow r;
        if (!t.empty() && (upper(t.back()) == "AND" || upper(t.back()) == "OR")) { r.conn = upper(t.back()) == "AND" ? 1 : 2; t.pop_back(); }
        size_t op = t.size();
        for (size_t i = 0; i < t.size(); ++i) for (const char* o : kOps) if (t[i] == o) { op = i; break; }
        if (op == t.size()) { f.structured = false; continue; }
        r.op = t[op];
        for (size_t i = 0; i < op; ++i) r.lhs += (i ? " " : "") + t[i];
        for (size_t i = op + 1; i < t.size(); ++i) r.rhs += (i > op + 1 ? " " : "") + t[i];
        f.conds.push_back(r);
    }
    while (!f.body.empty() && (f.body.back() == '\n' || f.body.back() == ' ')) f.body.pop_back();
    return f;
}
std::string writeActForm(const ActForm& f) {
    char h[160]; std::snprintf(h, sizeof h, " %s %d %g /\n", f.name.c_str(), std::max(1, f.runs), std::max(0.0, f.wait));
    std::string b = h;
    for (size_t i = 0; i < f.conds.size(); ++i) {
        const CondRow& r = f.conds[i];
        b += "  " + trim(r.lhs) + " " + r.op + " " + trim(r.rhs);
        if (i + 1 < f.conds.size()) b += r.conn == 2 ? " OR" : " AND";
        b += " /\n";
    }
    b += "/\n" + f.body + (f.body.empty() || f.body.back() == '\n' ? "" : "\n") + "ENDACTIO\n";
    return b;
}

struct Template { const char* label; const char* kw; const char* body; };
const Template kTemplates[] = {
    {"UDQ: field water cut + voidage replacement ratio", "UDQ",
     " DEFINE FUWCT FWPR / (FWPR + FOPR) /\n DEFINE FUVRR FVIR / FVPR /\n UNITS FUWCT 'frac' /\n UNITS FUVRR 'ratio' /\n/\n"},
    {"UDQ: well liquid rate and the average producer BHP", "UDQ",
     " DEFINE WULPR WOPR + WWPR /\n DEFINE FUAVBHP AVEA(WBHP 'P*') /\n/\n"},
    {"UDQ: injection target = 1.0 x produced reservoir volume (use as UDA 'FUINJ')", "UDQ",
     " DEFINE FUINJ 1.0 * FVPR /\n UNITS FUINJ 'RB/DAY' /\n/\n"},
    {"ACTIONX: shut producers above 90 % water cut", "ACTIONX",
     " WCUT_SHUT 10 /\n  WWCT 'P*' > 0.9 /\n/\nWELOPEN\n '?' SHUT /\n/\nENDACTIO\n"},
    {"ACTIONX: shut producers with a high GOR", "ACTIONX",
     " GOR_SHUT 10 /\n  WGOR 'P*' > 5 /\n/\nWELOPEN\n '?' SHUT /\n/\nENDACTIO\n"},
    {"ACTIONX: cut producer rates when the field pressure falls", "ACTIONX",
     " PLOW 1 /\n  FPR < 2000 /\n/\nWELTARG\n 'P*' ORAT 500 /\n/\nENDACTIO\n"},
    {"ACTIONX: change controls from a calendar date", "ACTIONX",
     " DATE2030 1 /\n  YEAR >= 2030 /\n/\nWELTARG\n 'P*' BHP 1000 /\n/\nENDACTIO\n"},
};
}  // namespace

void App::editRunLogic() {
    auto& X = case_.extraKeywords;
    const float fs = ImGui::GetFontSize();
    const int nSteps = std::max(1, int(case_.schedule.size()));
    static int sel = -1, tab = 0;
    static std::string formKey;           // body the form was built from
    static std::vector<UdqRow> urows;
    static ActForm aform;
    std::vector<int> items;
    for (int i = 0; i < int(X.size()); ++i) if (X[size_t(i)].name == "UDQ" || X[size_t(i)].name == "ACTIONX") items.push_back(i);
    if (sel >= int(X.size()) || (sel >= 0 && X[size_t(sel)].name != "UDQ" && X[size_t(sel)].name != "ACTIONX")) sel = -1;
    ui::hint("UDQ: user-defined quantities (FU field, WU per well, GU per group) evaluated after every time step. A UDQ name can "
             "be the target of a well control (UDA, below). ACTIONX: conditions checked at the end of every report step; when "
             "true the keywords of the body apply from the next step ('?' = the wells that satisfied a well condition). "
             "ResSim, OPM Flow and ECLIPSE run the same deck text.");
    // ---- list
    auto add = [&](const char* kw, const std::string& body, int step) {
        ExtraKeyword e; e.name = kw; e.section = DeckSection::Schedule; e.step = std::clamp(step, 0, nSteps - 1); e.body = body;
        e.note = "run-time logic (C31 editor)"; e.beforeControls = std::string(kw) == "UDQ";
        X.push_back(e); sel = int(X.size()) - 1; formKey.clear(); markDirty();
    };
    if (ui::beginTable("##logic", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY, ImVec2(-1, fs * std::clamp(1.6f + 1.35f * float(items.size()), 3.0f, 10.0f)))) {
        ui::tableColumn("Keyword", ImGuiTableColumnFlags_WidthFixed, fs * 6);
        ui::tableColumn("From step", ImGuiTableColumnFlags_WidthFixed, fs * 5);
        ui::tableColumn("Contents");
        ui::tableColumn("Check", ImGuiTableColumnFlags_WidthFixed, fs * 12);
        ImGui::TableHeadersRow();
        for (int i : items) {
            const ExtraKeyword& e = X[size_t(i)];
            ImGui::TableNextRow(); ImGui::PushID(i);
            ImGui::TableNextColumn();
            if (ImGui::Selectable(e.name.c_str(), sel == i, ImGuiSelectableFlags_SpanAllColumns)) { sel = i; formKey.clear(); }
            ImGui::TableNextColumn(); ImGui::Text("%d", e.step + 1);
            ImGui::TableNextColumn();
            std::string err;
            if (e.name == "UDQ") {
                udq::Set s; s.add(e.body, e.step, err);
                std::string nm; for (const auto& d : s.definitions()) { if (!nm.empty()) nm += ", "; nm += d.name; if (nm.size() > 80) { nm += " ..."; break; } }
                ImGui::TextUnformatted(nm.empty() ? "(no DEFINE / ASSIGN)" : nm.c_str());
            } else {
                udq::Action a; udq::parseAction(e.body, e.step, a, err);
                std::string c = a.conditionText; std::replace(c.begin(), c.end(), '\n', ' ');
                ImGui::Text("%s (max %d): %s", a.name.c_str(), a.maxRuns, c.c_str());
            }
            ImGui::TableNextColumn();
            if (err.empty()) ImGui::TextColored(col4(Pal::Ok), "OK");
            else { ImGui::TextColored(col4(Pal::Warn), "%s", err.c_str()); if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", err.c_str()); }
            ImGui::PopID();
        }
        ui::endTable();
    }
    static int newStep = 0;
    ImGui::SetNextItemWidth(fs * 6); ImGui::InputInt("at report step##lg", &newStep); newStep = std::clamp(newStep, 1, nSteps) ;
    ImGui::SameLine(); if (ImGui::Button("New UDQ block")) add("UDQ", " DEFINE FUNEW FOPR /\n/\n", newStep - 1);
    ImGui::SameLine(); if (ImGui::Button("New ACTIONX")) add("ACTIONX", " ACT" + std::to_string(items.size() + 1) + " 1 /\n  FPR < 2000 /\n/\nWELOPEN\n '?' SHUT /\n/\nENDACTIO\n", newStep - 1);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 12);
    if (ImGui::BeginCombo("##tpl", "From a template...", ImGuiComboFlags_HeightLarge)) {
        for (const auto& t : kTemplates) if (ImGui::Selectable(t.label)) add(t.kw, t.body, newStep - 1);
        ImGui::EndCombo();
    }
    if (sel >= 0) {
        ImGui::SameLine();
        if (ui::dangerButton("Delete")) { X.erase(X.begin() + sel); sel = -1; formKey.clear(); markDirty(); return; }
    }
    if (sel < 0) { if (!items.empty()) ImGui::TextDisabled("Select a row to edit it."); }
    else {
        ExtraKeyword& e = X[size_t(sel)];
        ImGui::Separator();
        int st = e.step + 1;
        ImGui::SetNextItemWidth(fs * 6);
        if (ImGui::InputInt("From report step##lgs", &st)) { e.step = std::clamp(st, 1, nSteps) - 1; markDirty(); }
        ImGui::SameLine();
        if (ImGui::RadioButton("Form", tab == 0)) tab = 0;
        ImGui::SameLine(); if (ImGui::RadioButton("Deck text", tab == 1)) tab = 1;
        if (formKey != e.name + "\x1f" + e.body) {
            formKey = e.name + "\x1f" + e.body;
            if (e.name == "UDQ") urows = parseUdqRows(e.body); else aform = parseActForm(e.body);
        }
        if (tab == 1) {
            ImGui::TextDisabled("%s", e.name.c_str());
            if (ui::inputString("##lgraw", e.body, ImVec2(-1, fs * 14), true)) markDirty();
            if (e.name == "ACTIONX") ImGui::TextDisabled("(the text must end with ENDACTIO)");
        } else if (e.name == "UDQ") {
            if (hasComments(e.body)) ui::hint("Editing in the form rewrites the block without its comments (use 'Deck text' to keep them).");
            bool ch = false; int del = -1;
            if (ui::beginTable("##udqrows", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
                ui::tableColumn("Record", ImGuiTableColumnFlags_WidthFixed, fs * 6);
                ui::tableColumn("Name", ImGuiTableColumnFlags_WidthFixed, fs * 7);
                ui::tableColumn("Expression / value / mode / units");
                ui::tableColumn("Check", ImGuiTableColumnFlags_WidthFixed, fs * 9);
                ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 2);
                ImGui::TableHeadersRow();
                for (int r = 0; r < int(urows.size()); ++r) {
                    UdqRow& u = urows[size_t(r)];
                    ImGui::TableNextRow(); ImGui::PushID(r);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ch |= ImGui::Combo("##k", &u.kind, kUdqKinds, 4);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1);
                    if (ui::inputString("##n", &u.name, ImGuiInputTextFlags_CharsUppercase)) ch = true;
                    ImGui::TableNextColumn();
                    std::string chk = "OK";
                    if (u.kind == 1) {
                        ImGui::SetNextItemWidth(fs * 6); ch |= ui::inputString("##sel", &u.sel); ImGui::SameLine();
                        if (ImGui::IsItemHovered()) ImGui::SetTooltip("well / group pattern (WU / GU only; empty = all)");
                    }
                    if (u.kind == 2) {
                        static const char* md[] = {"ON", "OFF", "NEXT"};
                        int m = upper(u.text) == "OFF" ? 1 : upper(u.text) == "NEXT" ? 2 : 0;
                        ImGui::SetNextItemWidth(-1); if (ImGui::Combo("##m", &m, md, 3)) { u.text = md[m]; ch = true; }
                    } else { ImGui::SetNextItemWidth(-1); ch |= ui::inputString("##t", &u.text); }
                    const std::string nm = upper(u.name);
                    if (nm.size() < 3 || (nm.rfind("FU", 0) && nm.rfind("WU", 0) && nm.rfind("GU", 0) && nm.rfind("SU", 0) && nm.rfind("CU", 0) && nm.rfind("RU", 0))) chk = "name must start FU / WU / GU";
                    else if (u.kind == 0) { udq::Expr ex; std::string err; if (!udq::parseExpr(u.text, ex, err)) chk = err; }
                    else if (u.kind == 1) { char* end = nullptr; std::strtod(u.text.c_str(), &end); if (u.text.empty() || (end && *end)) chk = "value must be a number"; }
                    ImGui::TableNextColumn();
                    if (chk == "OK") ImGui::TextColored(col4(Pal::Ok), "OK"); else { ImGui::TextColored(col4(Pal::Warn), "%s", chk.c_str()); if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", chk.c_str()); }
                    ImGui::TableNextColumn(); if (ImGui::SmallButton("x")) del = r;
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (del >= 0) { urows.erase(urows.begin() + del); ch = true; }
            if (ImGui::Button("+ DEFINE")) { UdqRow u; u.name = "FUNEW"; u.text = "FOPR"; urows.push_back(u); ch = true; }
            ImGui::SameLine(); if (ImGui::Button("+ ASSIGN")) { UdqRow u; u.kind = 1; u.name = "FUNEW"; u.text = "0"; urows.push_back(u); ch = true; }
            ImGui::SameLine(); if (ImGui::Button("+ UPDATE")) { UdqRow u; u.kind = 2; u.name = urows.empty() ? "FUNEW" : urows.back().name; u.text = "ON"; urows.push_back(u); ch = true; }
            ImGui::SameLine(); if (ImGui::Button("+ UNITS")) { UdqRow u; u.kind = 3; u.name = urows.empty() ? "FUNEW" : urows.back().name; urows.push_back(u); ch = true; }
            ui::help("Expressions: numbers, F* / W* ['pattern'] / G* ['group'] summary quantities (FOPR, FWPR, FVPR, FPR, WWCT, WBHP, "
                     "WOPT ...), other UDQs, TIME DAYS DAY MNTH YEAR, + - * / ^, == != < > <= >= (1 / 0), ABS EXP LN LOG NINT SQRT "
                     "DEF UNDEF IDV SORTA SORTD and the set functions SUM AVEA AVEG AVEH MAX MIN PROD NORM1 NORM2 NORMI.");
            if (ch) { e.body = writeUdqRows(urows); formKey = e.name + "\x1f" + e.body; markDirty(); }
        } else {   // ACTIONX
            bool ch = false;
            if (!aform.structured) ui::hint("This action uses parentheses or a form the condition rows cannot show - edit it as deck text.");
            if (ui::beginProps("##acth")) {
                ch |= ui::propString("Name", &aform.name);
                ch |= ui::propI("Max. runs", &aform.runs, 1, 100000, "how many times the action may fire");
                ch |= ui::propD("Min. interval (days)", &aform.wait, "%.4g", "minimum simulated time between two runs");
                ui::endProps();
            }
            ImGui::BeginDisabled(!aform.structured);
            int del = -1;
            if (ui::beginTable("##conds", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
                ui::tableColumn("Left side (quantity ['pattern'] / UDQ / DAY ...)");
                ui::tableColumn("Op", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
                ui::tableColumn("Right side");
                ui::tableColumn("Then", ImGuiTableColumnFlags_WidthFixed, fs * 5);
                ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 2);
                ImGui::TableHeadersRow();
                for (int r = 0; r < int(aform.conds.size()); ++r) {
                    CondRow& c = aform.conds[size_t(r)];
                    ImGui::TableNextRow(); ImGui::PushID(r);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ch |= ui::inputString("##l", &c.lhs);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1);
                    int o = 0; for (int k = 0; k < 6; ++k) if (c.op == kOps[k]) o = k;
                    if (ImGui::Combo("##o", &o, kOps, 6)) { c.op = kOps[o]; ch = true; }
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ch |= ui::inputString("##r", &c.rhs);
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1);
                    if (r + 1 < int(aform.conds.size())) { int cn = std::max(1, c.conn); if (ImGui::Combo("##c", &cn, kConn + 0, 3) && cn > 0) { c.conn = cn; ch = true; } }
                    else ImGui::TextDisabled("-");
                    ImGui::TableNextColumn(); if (ImGui::SmallButton("x")) del = r;
                    ImGui::PopID();
                }
                ui::endTable();
            }
            if (del >= 0) { aform.conds.erase(aform.conds.begin() + del); ch = true; }
            if (ImGui::Button("+ condition")) { CondRow c; c.lhs = "FPR"; c.rhs = "2000"; if (!aform.conds.empty()) aform.conds.back().conn = 1; aform.conds.push_back(c); ch = true; }
            ImGui::EndDisabled();
            ImGui::TextDisabled("Body (keywords applied when the conditions hold; '?' = the triggering wells)");
            if (ui::inputString("##abody", aform.body, ImVec2(-1, fs * 8), true)) ch = true;
            if (ch && aform.structured) { e.body = writeActForm(aform); formKey = e.name + "\x1f" + e.body; markDirty(); }
            udq::Action a; std::string err;
            if (!udq::parseAction(e.body, e.step, a, err)) ImGui::TextColored(col4(Pal::Warn), "%s", err.c_str());
            else ImGui::TextColored(col4(Pal::Ok), "Parsed: %s, max %d run(s), condition: %s", a.name.c_str(), a.maxRuns, trim(a.conditionText).c_str());
        }
    }
    // ---------------------------------------------------------------- C38 UDT (user-defined tables for UDQ expressions)
    ImGui::Separator();
    ui::header("UDT - user-defined tables (TU_NAME[x] in UDQ expressions)");
    ui::hint("1-D lookup tables used inside UDQ DEFINE expressions, e.g. DEFINE WU_BHP TU_FBHP[WOPR] UMIN WBHP. Interpolation: "
             "LC linear with constant extrapolation, LL linear with linear extrapolation, NV nearest value. Evaluated by the ResSim "
             "engine and written verbatim for OPM Flow / ECLIPSE (UDTDIMS in RUNSPEC).");
    {
        static int udtSel = -1;
        std::vector<int> idx; for (int i = 0; i < int(X.size()); ++i) if (X[size_t(i)].name == "UDT") idx.push_back(i);
        auto writeUdt = [](const std::vector<udq::Table>& T) {
            std::ostringstream o;
            for (const auto& t : T) {
                o << " '" << t.name << "' 1 /\n '" << t.interp << "'";
                for (double v : t.x) o << ' ' << v;
                o << " /\n ";
                for (double v : t.y) o << ' ' << v;
                o << " /\n/\n";
            }
            o << "/\n";
            return o.str();
        };
        if (ImGui::Button("+ UDT table")) {
            ExtraKeyword e; e.name = "UDT"; e.section = DeckSection::Schedule; e.step = std::clamp(selStep_, 0, nSteps - 1);
            udq::Table t; t.name = "TU_NEW"; t.interp = "LC"; t.x = {0, 1000}; t.y = {100, 50};
            e.body = writeUdt({t});
            X.push_back(e); udtSel = int(X.size()) - 1; markDirty();
            idx.push_back(udtSel);
        }
        for (int i : idx) {
            ImGui::SameLine();
            std::vector<udq::Table> T; std::string err; udq::parseTables(X[size_t(i)].body, 0, T, err);
            const std::string lab = (T.empty() ? std::string("UDT") : T[0].name) + "##udt" + std::to_string(i);
            if (ImGui::Selectable(lab.c_str(), udtSel == i, 0, ImVec2(ImGui::CalcTextSize(lab.c_str()).x, 0))) udtSel = i;
        }
        if (udtSel >= 0 && udtSel < int(X.size()) && X[size_t(udtSel)].name == "UDT") {
            ExtraKeyword& e = X[size_t(udtSel)];
            std::vector<udq::Table> T; std::string err;
            const bool ok = udq::parseTables(e.body, 0, T, err);
            if (!ok) ImGui::TextColored(col4(Pal::Warn), "%s", err.c_str());
            bool ch = false;
            int st = e.step + 1;
            if (ImGui::InputInt("From report step##udt", &st)) { e.step = std::clamp(st - 1, 0, nSteps - 1); markDirty(); }
            for (size_t ti = 0; ti < T.size(); ++ti) {
                udq::Table& t = T[ti];
                ImGui::PushID(int(ti));
                char nb[32]; std::snprintf(nb, sizeof(nb), "%s", t.name.c_str());
                ImGui::SetNextItemWidth(fs * 8);
                if (ImGui::InputText("name", nb, sizeof(nb), ImGuiInputTextFlags_CharsUppercase)) { t.name = nb; ch = true; }
                ImGui::SameLine();
                static const char* kI[] = {"LC", "LL", "NV"};
                int ii = t.interp == "LL" ? 1 : t.interp == "NV" ? 2 : 0;
                ImGui::SetNextItemWidth(fs * 5);
                if (ImGui::Combo("interpolation", &ii, kI, 3)) { t.interp = kI[ii]; ch = true; }
                if (t.name.empty() || t.name[0] != 'T') { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "names start with T (TU_...)"); }
                int del = -1;
                if (ui::beginTable("##udtrows", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
                    ui::tableColumn("x (argument)"); ui::tableColumn("y (value)"); ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 2);
                    ImGui::TableHeadersRow();
                    for (size_t r = 0; r < t.x.size() && r < t.y.size(); ++r) {
                        ImGui::PushID(int(r)); ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); if (ImGui::InputDouble("##x", &t.x[r], 0, 0, "%.6g")) ch = true;
                        ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); if (ImGui::InputDouble("##y", &t.y[r], 0, 0, "%.6g")) ch = true;
                        ImGui::TableNextColumn(); if (ImGui::SmallButton("x")) del = int(r);
                        ImGui::PopID();
                    }
                    ui::endTable();
                }
                if (del >= 0 && t.x.size() > 1) { t.x.erase(t.x.begin() + del); t.y.erase(t.y.begin() + del); ch = true; }
                if (ImGui::SmallButton("+ row")) { const double xl = t.x.empty() ? 0 : t.x.back(); t.x.push_back(xl + 1); t.y.push_back(t.y.empty() ? 0 : t.y.back()); ch = true; }
                // preview of the lookup curve
                if (t.x.size() >= 1) {
                    const double x0 = t.x.front(), x1 = t.x.back(), pad = std::max(1e-9, 0.15 * (x1 - x0 + 1e-9));
                    std::string pv = "y(" ;
                    char b[96]; std::snprintf(b, sizeof(b), "%.4g) = %.4g, y(%.4g) = %.4g (outside the table)", x0 - pad, t.eval(x0 - pad), x1 + pad, t.eval(x1 + pad));
                    ImGui::TextDisabled("%s%s", pv.c_str(), b);
                }
                ImGui::PopID();
            }
            if (ch) { e.body = writeUdt(T); markDirty(); }
            if (ui::dangerButton("Delete UDT keyword")) { X.erase(X.begin() + udtSel); udtSel = -1; markDirty(); }
            bool dims = false;
            for (const auto& k : X) if (k.name == "UDTDIMS") dims = true;
            if (!dims) ImGui::TextDisabled("UDTDIMS is added to RUNSPEC by the deck writer when a UDT is present.");
        }
    }
    // ---------------------------------------------------------------- C37 PYACTION
    ImGui::Separator();
    editPyActions();
    // ---------------------------------------------------------------- UDA
    ImGui::Separator();
    ui::header("UDA - a UDQ as a well-control target");
    std::vector<std::string> udqNames;
    for (const auto& e : X) if (e.name == "UDQ") {
        udq::Set s; std::string err; s.add(e.body, e.step, err);
        for (const auto& d : s.definitions()) if ((d.var == 'F' || d.var == 'W') && std::find(udqNames.begin(), udqNames.end(), d.name) == udqNames.end()) udqNames.push_back(d.name);
    }
    int nUda = 0;
    for (int s = 0; s < int(case_.schedule.size()); ++s)
        for (auto& c : case_.schedule[size_t(s)].controls)
            for (size_t u = 0; u < c.uda.size(); ++u) {
                ++nUda; ImGui::PushID(int(s * 100000 + nUda));
                ImGui::BulletText("step %d  %s  %s = '%s'", s + 1, c.well.c_str(), c.uda[u].first.c_str(), c.uda[u].second.c_str());
                ImGui::SameLine(); if (ImGui::SmallButton("remove")) { c.uda.erase(c.uda.begin() + long(u)); markDirty(); ImGui::PopID(); return; }
                ImGui::PopID();
            }
    if (!nUda) ImGui::TextDisabled("No control uses a UDQ target yet.");
    static int uStep = 1, uCtl = 0, uTgt = 0, uName = 0;
    uStep = std::clamp(uStep, 1, nSteps);
    const auto* ctrls = case_.schedule.empty() ? nullptr : &case_.schedule[size_t(uStep - 1)].controls;
    ImGui::SetNextItemWidth(fs * 5); ImGui::InputInt("step##uda", &uStep); uStep = std::clamp(uStep, 1, nSteps);
    if (!ctrls || ctrls->empty() || udqNames.empty()) {
        ImGui::TextDisabled(udqNames.empty() ? "Define a FU / WU quantity first." : "No well control (WCONPROD / WCONINJE) at this step.");
        return;
    }
    uCtl = std::clamp(uCtl, 0, int(ctrls->size()) - 1);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8);
    if (ImGui::BeginCombo("control##uda", (*ctrls)[size_t(uCtl)].well.c_str())) {
        for (int i = 0; i < int(ctrls->size()); ++i) if (ImGui::Selectable((*ctrls)[size_t(i)].well.c_str(), i == uCtl)) uCtl = i;
        ImGui::EndCombo();
    }
    bool inj = false;
    for (const auto& w : case_.wells) if (w.name == (*ctrls)[size_t(uCtl)].well) inj = w.type != WellType::Producer;
    static const char* prodT[] = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "THP"};
    static const char* injT[] = {"RATE", "RESV", "BHP", "THP"};
    const int nT = inj ? 4 : 7;
    uTgt = std::clamp(uTgt, 0, nT - 1); uName = std::clamp(uName, 0, int(udqNames.size()) - 1);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::Combo("target##uda", &uTgt, inj ? injT : prodT, nT);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7);
    if (ImGui::BeginCombo("UDQ##uda", udqNames[size_t(uName)].c_str())) {
        for (int i = 0; i < int(udqNames.size()); ++i) if (ImGui::Selectable(udqNames[size_t(i)].c_str(), i == uName)) uName = i;
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    if (ImGui::Button("Attach")) {
        WellControl& c = case_.schedule[size_t(uStep - 1)].controls[size_t(uCtl)];
        const std::string t = inj ? injT[uTgt] : prodT[uTgt];
        bool done = false;
        for (auto& u : c.uda) if (u.first == t) { u.second = udqNames[size_t(uName)]; done = true; }
        if (!done) c.uda.push_back({t, udqNames[size_t(uName)]});
        markDirty();
        notify("UDA attached: " + c.well + " " + t + " = '" + udqNames[size_t(uName)] + "'");
    }
    ui::hint("The target follows the UDQ from the next time step on (a WU quantity: the well's own value); an undefined value "
             "keeps the previous target. The control mode itself is set in the well control.");
}
