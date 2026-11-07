// C38 production-data import (Model mode, Schedule): measured well data (CSV / TSV / whitespace table, e.g. a
// production-database export) -> column mapping -> quality control -> rates / uptime -> history keywords
// (WCONHIST / WCONINJH / WEFAC at DATES) merged into the deck (a new .DATA that is then opened), or saved as an include
// file, plus an observed-data CSV for history matching. Engine: engine/src/proddata.cpp.
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "implot.h"
#include "plotkit.h"
#include "ressim/io.h"
#include "ressim/proddata.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

using namespace ressim;
namespace fs = std::filesystem;

namespace {
struct PdState {
    std::string text, source, status;
    pdata::Options opt;
    pdata::Table table;
    bool parsed = false, ok = false;
    int ctrl = 4;                 // RESV
    int showSev = 1;              // minimum severity shown
    int well = 0, period = 0, srcMode = 0;
    std::vector<pdata::Issue> mergeIssues;
};
PdState g_pd;
const char* kCtrl[] = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"};
const char* kSev[] = {"info", "warning", "error"};

std::string readAll(const std::string& p) {
    std::ifstream f(fs::u8path(p), std::ios::binary);
    std::ostringstream s; s << f.rdbuf(); return s.str();
}
bool writeAll(const std::string& p, const std::string& s) {
    std::ofstream f(fs::u8path(p), std::ios::binary);
    f << s; return bool(f);
}
const char* kExample =
    "WELL,DATE,OIL,WATER,GAS,WINJ,BHP,DAYS\n"
    "PROD,2015-01-01,62000,0,8800000,0,3800,31\n"
    "PROD,2015-02-01,54000,120,7900000,0,3700,28\n"
    "PROD,2015-03-01,55000,560,8300000,0,3610,29\n"
    "PROD,2015-04-01,48000,1500,7700000,0,3520,30\n"
    "INJ,2015-02-01,0,0,0,30000,,28\n"
    "INJ,2015-03-01,0,0,0,42000,,31\n"
    "INJ,2015-04-01,0,0,0,40000,,30\n";
ImU32 sevCol(int s) { return s == 2 ? Pal::Err : s == 1 ? Pal::Warn : Pal::Muted; }
void pline(const char* label, const std::vector<double>& x, const std::vector<double>& y, ImU32 col, bool step) {
    const int n = int(std::min(x.size(), y.size()));
    if (n <= 0) return;
    if (fig::recorder().recording()) recSeries(label, x.data(), y.data(), n, col, 2.0f, false);
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = 2.0f * g_plotLineScale;
    if (step) ImPlot::PlotStairs(label, x.data(), y.data(), n, sp); else ImPlot::PlotLine(label, x.data(), y.data(), n, sp);
}
}  // namespace

void App::editProductionData() {
    PdState& S = g_pd;
    const float fs = ImGui::GetFontSize();
    if (!ImGui::CollapsingHeader("Production data import -> history (WCONHIST / WCONINJH / WEFAC)")) return;
    ImGui::PushID("pdimport");
    ui::hint("Measured well data (monthly volumes or daily rates, uptime, BHP / THP) from a CSV / TSV / whitespace table with a header "
             "row - one row per well and date. The table is checked (dates, negative values, duplicates, uptime, water-cut / GOR "
             "jumps, wells missing in the deck), turned into history rates and merged into the deck's SCHEDULE as WCONHIST / "
             "WCONINJH / WEFAC at the period dates. The result is a new .DATA that runs in ResSim, OPM Flow and ECLIPSE.");
    bool reparse = false;
    if (ImGui::Button("Load table...")) {
        std::string p;
        if (fileDialog(false, "Production tables (*.csv;*.tsv;*.txt;*.prn)\0*.csv;*.CSV;*.tsv;*.txt;*.TXT;*.prn\0All files (*.*)\0*.*\0", "csv", p)) {
            S.text = readAll(p); S.source = p; S.opt.roleOverride.clear(); reparse = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Paste")) { if (const char* c = ImGui::GetClipboardText()) { S.text = c; S.source = "(clipboard)"; S.opt.roleOverride.clear(); reparse = true; } }
    ImGui::SameLine();
    if (ImGui::Button("Example")) { S.text = kExample; S.source = "(example: SPE1-style PROD / INJ, field units)"; S.opt.roleOverride.clear(); reparse = true; }
    ImGui::SameLine();
    ImGui::TextDisabled("%s", S.source.empty() ? "no table loaded" : S.source.c_str());
    if (ImGui::TreeNode("Table text (editable)")) {
        if (ui::inputString("##pdtext", S.text, ImVec2(-1, fs * 8), true)) reparse = true;
        ImGui::TreePop();
    }
    // ---- options
    if (ui::beginProps("##pdopt")) {
        int vol = S.opt.volumes ? 0 : 1;
        static const char* kVol[] = {"volumes per period", "rates (per day)"};
        if (ui::propCombo("Values are", &vol, kVol, 2, "Volumes (monthly production) are divided by the period length or by the days on stream")) { S.opt.volumes = vol == 0; reparse = true; }
        int ends = S.opt.datesAreEnds ? 1 : 0;
        static const char* kEnd[] = {"start of the period", "end of the period"};
        if (ui::propCombo("The date marks the", &ends, kEnd, 2, "Start: the row's rates apply from its date to the next date. End: the row covers the days after the previous date up to and including its date (31 JAN = all of January).")) { S.opt.datesAreEnds = ends == 1; reparse = true; }
        int mf = S.opt.monthFirst ? 1 : 0;
        static const char* kMf[] = {"day / month / year", "month / day / year (US)"};
        if (ui::propCombo("Numeric dates", &mf, kMf, 2, "Only for 01/02/2020-style dates; ISO (2020-02-01) and 1-FEB-2020 are unambiguous")) { S.opt.monthFirst = mf == 1; reparse = true; }
        if (ui::propCombo("WCONHIST control mode", &S.ctrl, kCtrl, 6, "RESV (reservoir voidage of the measured rates) is the usual history-matching choice; ORAT honours the oil rate exactly")) { S.opt.prodControl = kCtrl[S.ctrl]; reparse = true; }
        reparse |= ui::propBool("Uptime -> WEFAC", &S.opt.useWefac, "With a DAYS / HOURS column: flowing rate = volume / days on, WEFAC = days on / period (cumulatives honour the downtime)");
        reparse |= ui::propBool("Replace the deck controls of imported wells", &S.opt.replaceControls, "Removes their WCONPROD / WCONHIST / WCONINJE / WCONINJH records so that the history alone controls them");
        reparse |= ui::propBool("Shut wells in periods without data", &S.opt.shutMissing);
        reparse |= ui::propBool("Write measured BHP / THP", &S.opt.writeBhp, "Observed BHP / THP in WCONHIST / WCONINJH (reported as WBHPH / WTHPH)");
        S.opt.writeThp = S.opt.writeBhp;
        reparse |= ui::propD("Liquid unit factor", &S.opt.liqFactor, "%.6g", "Multiplies oil / water / water-injection values (e.g. 1000 for thousand sm3)");
        reparse |= ui::propD("Gas unit factor", &S.opt.gasFactor, "%.6g", "Multiplies gas values (e.g. 1000 for Mscf -> scf in a deck that uses scf, 1000 for 1000 sm3)");
        reparse |= ui::propD("Outlier jump factor (WOR / GOR)", &S.opt.outlierFactor, "%.3g");
        ui::endProps();
    }
    if (reparse || (!S.parsed && !S.text.empty())) {
        S.opt.startYear = case_.startYear; S.opt.startMonth = case_.startMonth; S.opt.startDay = case_.startDay;
        S.ok = pdata::parseTable(S.text, S.opt, S.table);
        S.parsed = true; S.mergeIssues.clear();
        S.well = std::min(S.well, std::max(0, int(S.table.wells.size()) - 1));
    }
    if (S.text.empty()) { ImGui::PopID(); return; }
    const auto& T = S.table;
    // ---- column mapping
    if (!T.columns.empty() && ImGui::TreeNodeEx("Columns", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (S.opt.roleOverride.size() != T.columns.size()) S.opt.roleOverride.assign(T.columns.size(), -2);
        std::vector<std::string> roles{"(automatic)", "(ignore)"};
        for (int r = 0; r < pdata::C_NCOL; ++r) roles.push_back(pdata::colName(r));
        if (ui::beginTable("##pdcols", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
            ui::tableColumn("column"); ui::tableColumn("detected as"); ui::tableColumn("use as");
            ImGui::TableHeadersRow();
            for (size_t c = 0; c < T.columns.size(); ++c) {
                ImGui::PushID(int(c));
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(T.columns[c].c_str());
                ImGui::TableNextColumn(); ImGui::TextUnformatted(c < T.roles.size() ? pdata::colName(T.roles[c]) : "-");
                ImGui::TableNextColumn();
                int sel = S.opt.roleOverride[c] == -2 ? 0 : S.opt.roleOverride[c] + 2;
                ImGui::SetNextItemWidth(-1);
                if (ImGui::BeginCombo("##role", roles[size_t(sel)].c_str())) {
                    for (int k = 0; k < int(roles.size()); ++k) if (ImGui::Selectable(roles[size_t(k)].c_str(), k == sel)) { S.opt.roleOverride[c] = k == 0 ? -2 : k - 2; S.parsed = false; }
                    ImGui::EndCombo();
                }
                ImGui::PopID();
            }
            ui::endTable();
        }
        ImGui::TreePop();
    }
    // ---- QC
    const int ne = T.nErrors(), nw = T.nWarnings();
    ImGui::TextColored(col4(ne ? Pal::Err : nw ? Pal::Warn : Pal::Ok), "Quality control: %d error(s), %d warning(s), %d note(s)", ne, nw, int(T.issues.size()) - ne - nw);
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); ImGui::Combo("show from##sev", &S.showSev, kSev, 3);
    {
        int shown = 0; for (const auto& i : T.issues) shown += i.severity >= S.showSev;
        if (shown && ui::beginTable("##pdqc", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable,
                                    ImVec2(-1, fs * std::clamp(1.6f + 1.35f * float(shown), 3.0f, 10.0f)))) {
            ui::tableColumn("level", ImGuiTableColumnFlags_WidthFixed, fs * 4.5f); ui::tableColumn("line", ImGuiTableColumnFlags_WidthFixed, fs * 3);
            ui::tableColumn("well", ImGuiTableColumnFlags_WidthFixed, fs * 5); ui::tableColumn("date", ImGuiTableColumnFlags_WidthFixed, fs * 6);
            ui::tableColumn("message", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupScrollFreeze(0, 1); ImGui::TableHeadersRow();
            for (const auto& i : T.issues) {
                if (i.severity < S.showSev) continue;
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextColored(col4(sevCol(i.severity)), "%s", kSev[std::clamp(i.severity, 0, 2)]);
                ImGui::TableNextColumn(); if (i.row) ImGui::Text("%d", i.row);
                ImGui::TableNextColumn(); ImGui::TextUnformatted(i.well.c_str());
                ImGui::TableNextColumn(); ImGui::TextUnformatted(i.date.c_str());
                ImGui::TableNextColumn(); ImGui::TextWrapped("%s", i.msg.c_str());
            }
            ui::endTable();
        }
    }
    if (!S.ok) { ImGui::PopID(); return; }
    // ---- wells
    std::map<std::string, int> deckWell; for (const auto& w : case_.wells) deckWell[w.name] = 1;
    if (ImGui::TreeNodeEx("Wells", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginTable("##pdwells", 8, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable,
                           ImVec2(-1, fs * std::clamp(1.6f + 1.35f * float(T.wells.size()), 3.0f, 9.0f)))) {
            for (const char* h : {"well", "in deck", "periods", "first", "cum oil", "cum water", "cum gas", "cum injected"}) ui::tableColumn(h);
            ImGui::TableSetupScrollFreeze(0, 1); ImGui::TableHeadersRow();
            for (int k = 0; k < int(T.wells.size()); ++k) {
                const std::string& w = T.wells[size_t(k)];
                double co = 0, cw = 0, cg = 0, ci = 0; int np = 0, p0 = -1;
                for (const auto& x : T.wp) {
                    if (x.well != w) continue;
                    if (p0 < 0) p0 = x.period;
                    if (x.shut) continue;
                    ++np;
                    const double L = double(T.periods[size_t(x.period)].d1 - T.periods[size_t(x.period)].d0) * x.wefac;
                    co += x.orat * L; cw += x.wrat * L; cg += x.grat * L; ci += (x.winj + x.ginj) * L;
                }
                ImGui::TableNextRow(); ImGui::TableNextColumn();
                if (ImGui::Selectable(w.c_str(), S.well == k, ImGuiSelectableFlags_SpanAllColumns)) S.well = k;
                ImGui::TableNextColumn();
                if (deckWell.count(w)) ImGui::TextColored(col4(Pal::Ok), "yes"); else ImGui::TextColored(col4(Pal::Err), "NO");
                ImGui::TableNextColumn(); ImGui::Text("%d", np);
                ImGui::TableNextColumn(); if (p0 >= 0) { int y, m, d; pdata::civilFromDays(T.periods[size_t(p0)].d0, y, m, d); ImGui::Text("%04d-%02d-%02d", y, m, d); }
                ImGui::TableNextColumn(); ImGui::Text("%.5g", co);
                ImGui::TableNextColumn(); ImGui::Text("%.5g", cw);
                ImGui::TableNextColumn(); ImGui::Text("%.5g", cg);
                ImGui::TableNextColumn(); ImGui::Text("%.5g", ci);
            }
            ui::endTable();
        }
        ImGui::TreePop();
    }
    // ---- plot of the selected well (history rates as stairs)
    if (S.well < int(T.wells.size())) {
        const std::string& w = T.wells[size_t(S.well)];
        std::vector<double> t, qo, qw, qg, qi, pb;
        const double day0 = double(pdata::daysFromCivil(case_.startYear > 0 ? case_.startYear : 2000, case_.startMonth, case_.startDay));
        for (const auto& x : T.wp) {
            if (x.well != w) continue;
            const double f = x.shut ? 0 : x.wefac;
            t.push_back(double(T.periods[size_t(x.period)].d0) - day0);
            qo.push_back(x.orat * f); qw.push_back(x.wrat * f); qg.push_back(x.grat * f); qi.push_back((x.winj + x.ginj) * f);
            pb.push_back(x.bhp > 0 ? x.bhp : NAN);
        }
        if (!t.empty()) {   // close the last period
            const auto& lp = T.periods.back();
            t.push_back(double(lp.d1) - day0); qo.push_back(qo.back()); qw.push_back(qw.back()); qg.push_back(qg.back()); qi.push_back(qi.back()); pb.push_back(pb.back());
        }
        const std::string title = "History of " + w + " (calendar-average rates)";
        if (ImPlot::BeginPlot(title.c_str(), ImVec2(-1, fs * 16))) {
            recPanel(title.c_str(), "days from START", "liquid rate", "gas rate / BHP");
            ImPlot::SetupAxes("days from START", "liquid rate", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "gas rate / BHP", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Opposite);
            pline("oil", t, qo, Pal::Oil, true);
            pline("water", t, qw, Pal::Water, true);
            pline("injected", t, qi, Pal::Accent, true);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            pline("gas", t, qg, Pal::Gas, true);
            pline("BHP", t, pb, Pal::Warn, false);
            recEnd();
            ImPlot::EndPlot();
        }
    }
    // ---- keyword preview
    if (ImGui::TreeNode("Generated keywords (preview)")) {
        S.period = std::clamp(S.period, 0, std::max(0, int(T.periods.size()) - 1));
        ImGui::SetNextItemWidth(fs * 12);
        ImGui::SliderInt("period", &S.period, 0, std::max(0, int(T.periods.size()) - 1));
        ImGui::SameLine(); ImGui::TextDisabled("starts %s", pdata::dateString(T.periods[size_t(S.period)].d0).c_str());
        std::string k = pdata::periodKeywords(T, S.opt, S.period);
        ImGui::InputTextMultiline("##pdkw", k.data(), k.size() + 1, ImVec2(-1, fs * 9), ImGuiInputTextFlags_ReadOnly);
        ImGui::TreePop();
    }
    // ---- actions
    static const char* kSrc[] = {"current model (exported deck)", "original deck file (keeps its layout and INCLUDEs)"};
    const bool haveSrc = !srcDeck_.empty() && fs::exists(fs::u8path(srcDeck_));
    if (!haveSrc) S.srcMode = 0;
    ImGui::BeginDisabled(!haveSrc);
    ImGui::SetNextItemWidth(fs * 22); ImGui::Combo("merge into", &S.srcMode, kSrc, 2);
    ImGui::EndDisabled();
    if (haveSrc && S.srcMode == 1 && srcDeckVersion_ != caseVersion_)
        ImGui::TextColored(col4(Pal::Warn), "The model was edited after the import - the original file does not contain those edits.");
    if (ui::accentButton("Merge into deck and open...")) {
        std::string p = haveSrc ? (fs::u8path(srcDeck_).parent_path() / (fs::u8path(srcDeck_).stem().u8string() + "_HIST.DATA")).u8string() : std::string("HISTORY.DATA");
        if (fileDialog(true, "ECLIPSE deck (*.DATA)\0*.DATA;*.data\0", "DATA", p)) {
            const std::string deck = S.srcMode == 1 && haveSrc ? readAll(srcDeck_) : eclipseDeckToString(case_);
            std::vector<std::string> names; for (const auto& w : case_.wells) names.push_back(w.name);
            S.mergeIssues.clear();
            const std::string merged = pdata::mergeIntoDeck(deck, T, S.opt, S.mergeIssues, names);
            bool fatal = false; for (const auto& i : S.mergeIssues) fatal |= i.severity == 2;
            if (fatal) S.status = "Merge failed - see the messages below";
            else if (!writeAll(p, merged)) S.status = "Cannot write " + p;
            else if (maybeSaveChanges()) { dirty_ = false; S.status = "Merged into " + p; importDeck(p); }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Save history include...")) {
        std::string p = "HISTORY.INC";
        if (fileDialog(true, "Include file (*.INC;*.SCH)\0*.INC;*.inc;*.SCH;*.sch\0", "INC", p))
            S.status = writeAll(p, pdata::historyInclude(T, S.opt)) ? "Saved " + p + " (INCLUDE it in SCHEDULE after the existing report dates)" : "Cannot write " + p;
    }
    ImGui::SameLine();
    if (ImGui::Button("Save observed CSV...")) {
        std::string p = "OBSERVED.csv";
        if (fileDialog(true, "CSV (*.csv)\0*.csv\0", "csv", p))
            S.status = writeAll(p, pdata::observedCsv(T, S.opt)) ? "Saved " + p + " (observed rates for history matching)" : "Cannot write " + p;
    }
    if (!S.status.empty()) ImGui::TextColored(col4(S.status.rfind("Merged", 0) == 0 || S.status.rfind("Saved", 0) == 0 ? Pal::Ok : Pal::Warn), "%s", S.status.c_str());
    for (const auto& i : S.mergeIssues) ImGui::TextColored(col4(sevCol(i.severity)), "%s%s%s", i.well.c_str(), i.well.empty() ? "" : ": ", i.msg.c_str());
    ImGui::PopID();
}
