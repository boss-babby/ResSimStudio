// C42 geomodel import UI (see ui_geoimport.h):
//   * import wizard - Petrel-style "column assignment": delimiter / header / skip / decimal / null settings with a live
//     preview, a role + unit factor per column (auto-guessed from the names, user-overridable), several files with
//     the same layout in one go, the assignment remembered for files with the same columns
//   * well manager - add / edit / delete wells by hand, per-well trajectory editor (MD/INC/AZI or MD/X/Y/TVDSS,
//     paste from Excel, import a file into the selected well, vertical to TD), plan + section preview, tops editor
//   * 2-D seismic - several SEG-Y lines, section view with seeded tracking / manual picking / erasing, line
//     intersections with mis-ties, gridding of the line picks to a depth surface, curtains in the 3-D scene
#include "i18n.h"
#include "ui_geoimport.h"
#include "ressim/geo_tools.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "platform.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <set>

namespace {
using namespace ressim::geo;
std::string stemOf(const std::string& p) {
    const size_t sl = p.find_last_of("/\\");
    std::string n = p.substr(sl == std::string::npos ? 0 : sl + 1);
    const size_t d = n.rfind('.');
    return d == std::string::npos ? n : n.substr(0, d);
}
const char* kFilter[] = {
    "Well heads (*.txt;*.csv;*.dat;*.prn)\0*.txt;*.csv;*.dat;*.prn;*.tsv\0All files (*.*)\0*.*\0",
    "Deviation surveys (*.txt;*.dev;*.csv;*.dat)\0*.txt;*.dev;*.csv;*.dat;*.wellpath;*.tsv\0All files (*.*)\0*.*\0",
    "Well tops (*.txt;*.csv;*.dat)\0*.txt;*.csv;*.dat;*.tsv\0All files (*.*)\0*.*\0",
    "Points (*.xyz;*.txt;*.csv;*.dat)\0*.xyz;*.txt;*.csv;*.dat;*.pts\0All files (*.*)\0*.*\0",
    "Fault sticks (*.txt;*.dat;*.csv)\0*.txt;*.dat;*.csv;*.fsticks\0All files (*.*)\0*.*\0",
};
// remembered assignment per kind for a header signature
struct Remembered { std::vector<std::string> header; ColumnMap map; std::vector<int> unitSel; };
std::map<int, Remembered>& remembered() { static std::map<int, Remembered> m; return m; }
void parseNulls(GeoImportWiz& W) {
    W.fmt.nulls.clear(); W.fmt.nulls.push_back("");
    std::string cur;
    for (char c : W.nullsText + ",") { if (c == ',' || c == ';') { while (!cur.empty() && cur.back() == ' ') cur.pop_back(); size_t a = cur.find_first_not_of(' '); if (a != std::string::npos) W.fmt.nulls.push_back(cur.substr(a)); cur.clear(); } else cur += c; }
}
TextTable readSource(GeoImportWiz& W, int file, int maxRows, std::string& err) {
    TableFormat f = W.fmt; f.maxRows = maxRows;
    if (W.fromText) return parseTextTable(W.text, f, err);
    if (file < 0 || file >= int(W.files.size())) { err = "no file"; return TextTable(); }
    return readTextTable(W.files[size_t(file)], f, err);
}
void loadPreview(GeoImportWiz& W) {
    W.dirty = false; W.err.clear();
    parseNulls(W);
    W.preview = readSource(W, W.previewFile, 200, W.err);
    const int nc = W.preview.ncol;
    const auto& R = remembered();
    auto it = R.find(int(W.kind));
    if (it != R.end() && it->second.header == W.preview.header && int(it->second.map.role.size()) == nc && !W.mapTouched) {
        W.map = it->second.map; W.unitSel = it->second.unitSel; W.mapTouched = true;
    } else if (!W.mapTouched || int(W.map.role.size()) != nc) {
        W.map = guessColumns(W.preview, W.fmt, W.kind); W.unitSel.assign(size_t(nc), 0);
    }
    W.unitSel.resize(size_t(nc), 0);
}
bool numericRole(ColRole r) { return r != ColRole::Ignore && r != ColRole::WellName && r != ColRole::Surface && r != ColRole::FaultName && r != ColRole::StickId; }
ImPlotSpec lspec(ImVec4 c, float w) { ImPlotSpec s; s.LineColor = c; s.LineWeight = w; return s; }
ImPlotSpec mspec(ImPlotMarker m, float sz, ImVec4 c) { ImPlotSpec s; s.Marker = m; s.MarkerSize = sz; s.MarkerFillColor = s.MarkerLineColor = c; s.LineWeight = 0; return s; }
ImVec4 hzColour(int i) { static const ImVec4 c[] = {{1, 0.9f, 0.1f, 1}, {0.2f, 1, 0.4f, 1}, {0.3f, 0.8f, 1, 1}, {1, 0.45f, 0.9f, 1}, {1, 0.6f, 0.2f, 1}, {0.7f, 0.6f, 1, 1}}; return c[i % 6]; }
void ensureLineHz(App::GeoState& S) {
    if (S.lineHz.empty()) { LineHorizon h; h.name = "H1"; S.lineHz.push_back(h); }
    for (auto& h : S.lineHz) {
        h.t.resize(S.lines.size());
        for (size_t l = 0; l < S.lines.size(); ++l) h.t[l].resize(size_t(S.lines[l].ntr()), std::nan(""));
    }
    S.lineOn.resize(S.lines.size(), 1);
    S.curLineHz = std::clamp(S.curLineHz, 0, int(S.lineHz.size()) - 1);
}
// head edited in the manager / table: MD-INC-AZI trajectories follow, X/Y surveys shift
void headEdited(GeoWell& w, double ox, double oy) {
    if (w.traj.empty()) return;
    if (w.trajSrc == 0) minimumCurvature(w);
    else for (auto& s : w.traj) { s.x += w.x - ox; s.y += w.y - oy; }
}

// ---------------------------------------------------------------- apply the wizard to every source
void applyWizard(App::GeoState& S) {
    GeoImportWiz& W = S.wiz;
    parseNulls(W);
    const int nSrc = W.fromText ? 1 : int(W.files.size());
    int okFiles = 0;
    for (int f = 0; f < nSrc; ++f) {
        const std::string src = W.fromText ? std::string("pasted text") : W.files[size_t(f)];
        std::string err;
        TextTable t = readSource(W, f, 0, err);
        if (!err.empty() && t.rows.empty()) { S.msg("Error: " + err); continue; }
        ColumnMap m = W.map;
        for (size_t c = 0; c < m.scale.size() && c < W.unitSel.size(); ++c) if (W.unitSel[c] == 1) m.scale[c] = 0.3048; else if (W.unitSel[c] == 2) m.scale[c] = 1 / 0.3048; else if (W.unitSel[c] == 0) m.scale[c] = 1;
        if (t.header != W.preview.header || int(m.role.size()) != t.ncol) {
            m = guessColumns(t, W.fmt, W.kind);
            const std::string bad = validateMap(m, W.kind);
            if (!bad.empty()) { S.msg("Skipped " + src + ": different columns and no automatic assignment (" + bad + ") - import it on its own"); continue; }
            S.msg(src + ": different columns - assignment re-guessed");
        }
        ImportOptions o = W.opt;
        if (!W.targetWell.empty()) o.defaultWell = W.targetWell;
        else if (W.wellFrom == 2 && W.fixedWell[0]) o.defaultWell = W.fixedWell;
        else if (W.wellFrom == 1 && (t.meta.count("WELL NAME") || t.meta.count("WELL"))) o.defaultWell = t.meta.count("WELL NAME") ? t.meta["WELL NAME"] : t.meta["WELL"];
        else o.defaultWell = W.fromText ? std::string("WELL") : stemOf(src);
        o.defaultSet = W.fromText ? std::string("pasted") : stemOf(src);
        ImportResult r = buildImport(t, W.fmt, m, W.kind, o);
        if (!W.targetWell.empty()) {
            if (W.kind == ImportKind::Deviation && r.wells.size() > 1) { S.msg("note: " + std::to_string(r.wells.size()) + " wells in the survey - only the first goes into " + W.targetWell); r.wells.resize(1); }
            for (auto& w : r.wells) w.name = W.targetWell;
            for (auto& tp : r.tops) tp.well = W.targetWell;
        }
        MergeStats st;
        std::string what;
        switch (W.kind) {
            case ImportKind::WellHeads: st = mergeHeads(S.wells, r); what = std::to_string(r.wells.size()) + " heads (" + std::to_string(st.added) + " new, " + std::to_string(st.updated) + " updated)"; break;
            case ImportKind::Deviation: st = mergeSurveys(S.wells, r); what = std::to_string(r.wells.size()) + " survey(s), " + std::to_string(r.rowsUsed) + " stations" + (r.trajMode == 0 ? " (MD/INC/AZI)" : r.trajMode == 1 ? " (X/Y/depth)" : " (offsets)"); break;
            case ImportKind::Tops: st = mergeTops(S.wells, r.tops, W.createMissing); what = std::to_string(st.updated) + " tops attached"; break;
            case ImportKind::Points: {
                int n = 0;
                for (auto& ps : r.points) {
                    if (ps.second.size() < 3) { S.msg(ps.first + ": fewer than 3 points - not gridded"); continue; }
                    Surface lat = (!W.gridPoints || S.hz.empty()) ? latticeFor(ps.second, 120) : S.hz[0];
                    double x0 = 1e300, x1 = -1e300; for (auto& p : ps.second) { x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]); }
                    geostat::Variogram vg; vg.model = geostat::VarioModel::Gaussian; vg.rangeMax = vg.rangeMin = std::max(1.0, 0.3 * (x1 - x0));
                    Surface s = gridPoints(ps.second, lat, vg); s.name = ps.first;
                    bool rep = false; for (auto& h : S.hz) if (h.name == s.name) { S.pushUndo(); h = s; rep = true; }
                    if (!rep) S.hz.push_back(s);
                    ++S.hzRev; ++n;
                }
                what = std::to_string(r.points.size()) + " point set(s), " + std::to_string(r.rowsUsed) + " points, " + std::to_string(n) + " surface(s) gridded";
                break;
            }
            case ImportKind::FaultSticks: for (auto& x : r.faults) S.faults.push_back(x); what = std::to_string(r.faults.size()) + " fault(s)"; break;
            default: break;
        }
        S.msg(stemOf(src) + ": " + what + (r.rowsSkipped ? ", " + std::to_string(r.rowsSkipped) + " rows skipped" : std::string()));
        int nw = 0;
        for (const auto& w : r.warnings) if (nw++ < 5) S.msg("  warning: " + w);
        for (const auto& w : st.warnings) if (nw++ < 10) S.msg("  " + w);
        if (r.rowsUsed > 0) ++okFiles;
    }
    // remember the assignment for files with the same columns
    Remembered& R = remembered()[int(W.kind)]; R.header = W.preview.header; R.map = W.map; R.unitSel = W.unitSel;
    S.msg(std::string(importKindName(W.kind)) + ": " + std::to_string(okFiles) + " / " + std::to_string(nSrc) + " source(s) imported");
}
}  // namespace

// ================================================================ wizard
void geoOpenImport(App::GeoState& S, ImportKind kind, const std::string& targetWell) {
    std::vector<std::string> files;
    if (!plat::fileDialogMulti(kFilter[int(kind)], files) || files.empty()) return;
    GeoImportWiz& W = S.wiz;
    W.kind = kind; W.files = files; W.fromText = false; W.text.clear(); W.previewFile = 0; W.targetWell = targetWell;
    W.fmt.delim = 0; W.fmt.headerRow = -2; W.fmt.skipLines = 0; W.fmt.decimalComma = -1;
    W.mapTouched = false; W.dirty = true; W.requestOpen = true;
}
void geoOpenImportText(App::GeoState& S, ImportKind kind, const std::string& text, const std::string& targetWell) {
    GeoImportWiz& W = S.wiz;
    W.kind = kind; W.files.clear(); W.fromText = true; W.text = text; W.previewFile = 0; W.targetWell = targetWell;
    W.fmt.delim = 0; W.fmt.headerRow = -2; W.fmt.skipLines = 0; W.fmt.decimalComma = -1;
    W.mapTouched = false; W.dirty = true; W.requestOpen = true;
}

void geoImportWindow(App::GeoState& S) {
    GeoImportWiz& W = S.wiz;
    const char* title = "Import with column assignment##geowiz";
    if (W.requestOpen) { ImGui::OpenPopup(title); W.requestOpen = false; W.open = true; }
    if (!W.open) return;
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowSize(ImVec2(vp->WorkSize.x * 0.82f, vp->WorkSize.y * 0.82f), ImGuiCond_Appearing);
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + vp->WorkSize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal(title, &W.open)) return;
    if (W.dirty) loadPreview(W);
    const float fs = ImGui::GetFontSize();
    // ---- what / from where
    int kind = int(W.kind);
    ImGui::SetNextItemWidth(fs * 11);
    if (ImGui::Combo("data type", &kind, "Well heads\0Deviation surveys\0Well tops\0Points (XYZ)\0Fault sticks\0")) { W.kind = ImportKind(kind); W.mapTouched = false; W.dirty = true; }
    ImGui::SameLine();
    if (W.fromText) ImGui::TextDisabled("source: pasted text (%zu characters)", W.text.size());
    else {
        ImGui::TextDisabled("%zu file(s)", W.files.size());
        if (W.files.size() > 1) {
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 16);
            std::vector<std::string> names; for (auto& f : W.files) names.push_back(stemOf(f));
            if (ImGui::Combo("preview##pf", &W.previewFile, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &names, int(names.size()))) W.dirty = true;
        } else if (!W.files.empty()) { ImGui::SameLine(); ImGui::TextUnformatted(W.files[0].c_str()); }
    }
    if (!W.targetWell.empty()) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Accent), "-> into well %s", W.targetWell.c_str()); }
    // ---- format
    ImGui::SeparatorText(TR("File format"));
    static const char delims[] = {0, ',', ';', '\t', ' ', '|'};
    int di = 0; for (int i = 0; i < 6; ++i) if (W.fmt.delim == delims[i]) di = i;
    ImGui::SetNextItemWidth(fs * 8);
    if (ImGui::Combo("delimiter", &di, "Auto\0Comma ,\0Semicolon ;\0Tab\0Blanks\0Pipe |\0")) { W.fmt.delim = delims[di]; W.dirty = true; }
    if (W.fmt.delim == 0) { ImGui::SameLine(); ImGui::TextDisabled("(%s)", W.preview.delimUsed == ',' ? "comma" : W.preview.delimUsed == ';' ? "semicolon" : W.preview.delimUsed == '\t' ? "tab" : W.preview.delimUsed == '|' ? "pipe" : "blanks"); }
    ImGui::SameLine();
    int hm = W.fmt.headerRow == -2 ? 0 : W.fmt.headerRow == -1 ? 1 : 2;
    ImGui::SetNextItemWidth(fs * 7);
    if (ImGui::Combo("header", &hm, "Auto\0None\0Row\0")) { W.fmt.headerRow = hm == 0 ? -2 : hm == 1 ? -1 : 0; W.dirty = true; W.mapTouched = false; }
    if (hm == 2) { ImGui::SameLine(); int r = W.fmt.headerRow + 1; ImGui::SetNextItemWidth(fs * 5); if (ImGui::InputInt("##hr", &r)) { W.fmt.headerRow = std::max(0, r - 1); W.dirty = true; W.mapTouched = false; } }
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5);
    if (ImGui::InputInt("skip lines", &W.fmt.skipLines)) { W.fmt.skipLines = std::max(0, W.fmt.skipLines); W.dirty = true; }
    ImGui::SameLine();
    int dc = W.fmt.decimalComma + 1; ImGui::SetNextItemWidth(fs * 7);
    if (ImGui::Combo("decimal", &dc, "Auto\0Point .\0Comma ,\0")) { W.fmt.decimalComma = dc - 1; W.dirty = true; }
    char nb[256]; std::snprintf(nb, sizeof nb, "%s", W.nullsText.c_str());
    ImGui::SetNextItemWidth(fs * 16);
    if (ImGui::InputText("null values", nb, sizeof nb)) { W.nullsText = nb; W.dirty = true; }
    ImGui::SameLine();
    char cb[16]; std::snprintf(cb, sizeof cb, "%s", W.fmt.comment.c_str()); ImGui::SetNextItemWidth(fs * 4);
    if (ImGui::InputText("comment chars", cb, sizeof cb)) { W.fmt.comment = cb; W.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_WAND_SPARKLES " Guess columns again"))) { W.mapTouched = false; W.dirty = true; }
    ui::help("Roles are guessed from the column names (synonyms such as Easting / X / XCOORD, KB / RKB / Elevation, INCL / DEVI, AZIM ...), "
             "then from the column order. Change any role in the first row of the table; the unit row converts feet <-> metres per column. "
             "The assignment is remembered for files with the same columns.");
    if (!W.err.empty()) ImGui::TextColored(col4(Pal::Err), "%s", W.err.c_str());
    if (ImGui::TreeNode(TR("Raw file (first lines)"))) {
        ImGui::BeginChild("##raw", ImVec2(0, fs * 7), ImGuiChildFlags_Borders, ImGuiWindowFlags_HorizontalScrollbar);
        for (size_t k = 0; k < W.preview.rawHead.size() && k < 25; ++k) ImGui::Text("%3zu  %s", k + 1, W.preview.rawHead[k].c_str());
        ImGui::EndChild();
        ImGui::TreePop();
    }
    // ---- column assignment table
    ImGui::SeparatorText(TR("Column assignment"));
    const TextTable& T = W.preview;
    const auto roles = rolesFor(W.kind);
    std::vector<std::string> roleNames; for (auto r : roles) roleNames.push_back(colRoleName(r));
    const float optH = fs * 7.5f;
    const float tableH = std::max(fs * 8, ImGui::GetContentRegionAvail().y - optH);
    if (T.ncol > 0 && ImGui::BeginTable("##map", T.ncol + 1, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable,
                                        ImVec2(0, tableH))) {
        ImGui::TableSetupScrollFreeze(1, 3);
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, fs * 2.5f);
        for (int c = 0; c < T.ncol; ++c) ImGui::TableSetupColumn(T.header[size_t(c)].c_str(), ImGuiTableColumnFlags_WidthFixed, fs * 9.5f);
        // row: roles
        ImGui::TableNextRow();
        ImGui::TableNextColumn(); ImGui::TextDisabled("role");
        for (int c = 0; c < T.ncol; ++c) {
            ImGui::TableNextColumn(); ImGui::PushID(c);
            int sel = 0; for (size_t r = 0; r < roles.size(); ++r) if (roles[r] == W.map.role[size_t(c)]) sel = int(r);
            const bool mapped = sel > 0;
            if (mapped) ImGui::PushStyleColor(ImGuiCol_FrameBg, ImGui::ColorConvertU32ToFloat4((Pal::Accent & 0x00FFFFFF) | 0x90000000));
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::Combo("##role", &sel, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &roleNames, int(roleNames.size()))) {
                const ColRole nr = roles[size_t(sel)];
                if (nr != ColRole::Ignore) for (auto& r : W.map.role) if (r == nr) r = ColRole::Ignore;   // a role belongs to one column
                W.map.role[size_t(c)] = nr; W.mapTouched = true;
            }
            if (mapped) ImGui::PopStyleColor();
            ImGui::PopID();
        }
        // row: units
        ImGui::TableNextRow();
        ImGui::TableNextColumn(); ImGui::TextDisabled("unit");
        for (int c = 0; c < T.ncol; ++c) {
            ImGui::TableNextColumn(); ImGui::PushID(1000 + c);
            if (numericRole(W.map.role[size_t(c)])) {
                int& u = W.unitSel[size_t(c)];
                ImGui::SetNextItemWidth(u == 3 ? fs * 4.5f : -FLT_MIN);
                if (ImGui::Combo("##unit", &u, "as is\0ft -> m\0m -> ft\0x factor\0")) W.mapTouched = true;
                if (u == 3) { ImGui::SameLine(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##sc", &W.map.scale[size_t(c)], 0, 0, "%.6g"); }
            } else ImGui::TextDisabled("-");
            ImGui::PopID();
        }
        // row: names
        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        ImGui::TableNextColumn(); ImGui::TextDisabled(T.hasHeader ? "name" : "(none)");
        for (int c = 0; c < T.ncol; ++c) { ImGui::TableNextColumn(); ImGui::TextUnformatted(T.header[size_t(c)].c_str()); }
        // data
        ImGuiListClipper clip; clip.Begin(int(T.rows.size()));
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextDisabled("%d", r + 1);
                for (int c = 0; c < T.ncol; ++c) {
                    ImGui::TableNextColumn();
                    const std::string& v = T.rows[size_t(r)][size_t(c)];
                    if (W.map.role[size_t(c)] == ColRole::Ignore) ImGui::TextDisabled("%s", v.c_str());
                    else {
                        double d; const bool isNum = tableNumber(T, W.fmt, v, d);
                        if (numericRole(W.map.role[size_t(c)]) && !isNum) ImGui::TextColored(col4(Pal::Warn), "%s", v.empty() ? "(null)" : v.c_str());
                        else ImGui::TextUnformatted(v.c_str());
                    }
                }
            }
        ImGui::EndTable();
    }
    // ---- options + dry run
    ImGui::BeginChild("##wizopt", ImVec2(0, 0));
    if (W.kind == ImportKind::Deviation || W.kind == ImportKind::Tops) {
        if (W.targetWell.empty() && W.map.find(ColRole::WellName) < 0) {
            ImGui::TextUnformatted("No well column - well name from:"); ImGui::SameLine();
            ImGui::RadioButton(TR("file name"), &W.wellFrom, 0); ImGui::SameLine();
            ImGui::RadioButton(TR("file header (# WELL NAME:)"), &W.wellFrom, 1); ImGui::SameLine();
            ImGui::RadioButton(TR("fixed"), &W.wellFrom, 2);
            if (W.wellFrom == 2) { ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 8); ImGui::InputText("##fw", W.fixedWell, sizeof W.fixedWell); }
        }
        if (W.kind == ImportKind::Tops) { ImGui::Checkbox(TR("create wells that are not in the project"), &W.createMissing); }
    }
    if (W.kind == ImportKind::Points || W.kind == ImportKind::FaultSticks || W.kind == ImportKind::Tops) {
        ImGui::TextUnformatted("Z column:"); ImGui::SameLine();
        ImGui::RadioButton(TR("auto"), &W.opt.zSign, 0); ImGui::SameLine(); ImGui::RadioButton(TR("elevation (negative down)"), &W.opt.zSign, 1); ImGui::SameLine(); ImGui::RadioButton(TR("depth (positive down)"), &W.opt.zSign, 2);
    }
    if (W.kind == ImportKind::Points) { ImGui::Checkbox(TR("grid on the lattice of the first horizon (else an automatic 120-node lattice)"), &W.gridPoints); }
    // dry run on the preview rows
    const std::string bad = validateMap(W.map, W.kind);
    if (!bad.empty()) ImGui::TextColored(col4(Pal::Err), ICON_TRIANGLE_ALERT " %s", bad.c_str());
    else {
        ColumnMap m = W.map;
        for (size_t c = 0; c < m.scale.size() && c < W.unitSel.size(); ++c) if (W.unitSel[c] == 1) m.scale[c] = 0.3048; else if (W.unitSel[c] == 2) m.scale[c] = 1 / 0.3048; else if (W.unitSel[c] == 0) m.scale[c] = 1;
        ImportOptions o = W.opt; o.defaultWell = W.targetWell.empty() ? (W.fromText || W.files.empty() ? std::string("WELL") : stemOf(W.files[size_t(std::clamp(W.previewFile, 0, int(W.files.size()) - 1))])) : W.targetWell;
        const ImportResult r = buildImport(T, W.fmt, m, W.kind, o);
        char b[300];
        if (W.kind == ImportKind::WellHeads) std::snprintf(b, sizeof b, "%zu well(s)", r.wells.size());
        else if (W.kind == ImportKind::Deviation) std::snprintf(b, sizeof b, "%zu survey(s), %s", r.wells.size(), r.trajMode == 0 ? "MD / INC / AZI (minimum curvature)" : r.trajMode == 1 ? "X / Y / depth stations" : "DX / DY offsets + depth");
        else if (W.kind == ImportKind::Tops) std::snprintf(b, sizeof b, "%zu top(s)", r.tops.size());
        else if (W.kind == ImportKind::Points) std::snprintf(b, sizeof b, "%zu point set(s)", r.points.size());
        else std::snprintf(b, sizeof b, "%zu fault(s)", r.faults.size());
        ImGui::TextColored(col4(r.rowsUsed ? Pal::Ok : Pal::Warn), ICON_CHECK " preview: %d rows used, %d skipped -> %s   (%lld data rows in %s)", r.rowsUsed, r.rowsSkipped, b, T.totalRows, W.fromText ? "the text" : "this file");
        for (size_t k = 0; k < r.warnings.size() && k < 2; ++k) ImGui::TextColored(col4(Pal::Warn), "%s", r.warnings[k].c_str());
    }
    ImGui::Separator();
    const int nSrc = W.fromText ? 1 : int(W.files.size());
    char lb[96]; std::snprintf(lb, sizeof lb, ICON_DOWNLOAD " Import %d source%s", nSrc, nSrc == 1 ? "" : "s");
    ImGui::BeginDisabled(!bad.empty());
    if (ui::accentButton(lb)) { applyWizard(S); W.open = false; ImGui::CloseCurrentPopup(); }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button(TR("Cancel"))) { W.open = false; ImGui::CloseCurrentPopup(); }
    if (nSrc > 1) { ImGui::SameLine(); ImGui::TextDisabled("all files use this format and assignment (files with other columns are re-guessed)"); }
    ImGui::EndChild();
    ImGui::EndPopup();
}

// ================================================================ LAS (several files)
int geoImportLas(App::GeoState& S, const std::vector<std::string>& paths) {
    int n = 0;
    for (const auto& p : paths) {
        std::string e; WellLog l = readLas(p, e);
        if (!e.empty()) { S.msg("Error: " + e); continue; }
        const std::string nm = l.well.empty() ? stemOf(p) : l.well;
        int i = matchWell(S.wells, nm);
        if (i < 0 && !l.well.empty()) i = matchWell(S.wells, stemOf(p));
        if (i < 0) { GeoWell g; g.name = nm; S.wells.push_back(g); i = int(S.wells.size()) - 1; S.msg("well " + nm + " created (no head: at 0,0 - import or type the head)"); }
        S.wells[size_t(i)].logs.push_back(l);
        S.msg("LAS " + stemOf(p) + " -> " + S.wells[size_t(i)].name + ": " + std::to_string(l.curves.size()) + " curves, " + std::to_string(l.md.size()) + " samples");
        ++n;
    }
    return n;
}

// ================================================================ well manager
void geoWellManager(App::GeoState& S) {
    GeoWellEdit& E = S.we;
    if (!E.open) return;
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextWindowSize(ImVec2(fs * 70, fs * 40), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin(ICON_PIN " Well manager##geowells", &E.open, ImGuiWindowFlags_NoDocking)) { ImGui::End(); return; }
    E.sel = std::clamp(E.sel, -1, int(S.wells.size()) - 1);
    if (E.sel < 0 && !S.wells.empty()) E.sel = 0;
    // ---------------- left: wells
    ImGui::BeginChild("##wl", ImVec2(ImGui::GetContentRegionAvail().x * 0.46f, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    if (ImGui::Button(TR(ICON_PLUS " Add well"))) {
        GeoWell w; int k = int(S.wells.size()) + 1;
        auto exists = [&](const std::string& n) { for (auto& x : S.wells) if (wellKey(x.name) == wellKey(n)) return true; return false; };
        do { w.name = "W-" + std::to_string(k++); } while (exists(w.name));
        if (!S.wells.empty()) { w.x = S.wells.back().x + 100; w.y = S.wells.back().y; w.kb = S.wells.back().kb; }
        else if (!S.hz.empty() && S.hz[0].valid()) S.hz[0].nodeXY(S.hz[0].nx / 2, S.hz[0].ny / 2, w.x, w.y);
        S.wells.push_back(w); E.sel = int(S.wells.size()) - 1;
    }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_CLIPBOARD_PASTE " Paste heads"))) { const char* c = ImGui::GetClipboardText(); if (c && *c) geoOpenImportText(S, ImportKind::WellHeads, c); }
    ui::help("Paste rows copied from Excel / a text file (name, x, y, kb, td in any order) - the column-assignment dialog opens");
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_FILE " Import heads..."))) geoOpenImport(S, ImportKind::WellHeads);
    if (ImGui::BeginTable("##wt", 8, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("well", ImGuiTableColumnFlags_WidthStretch, 1.4f);
        ImGui::TableSetupColumn("X", ImGuiTableColumnFlags_WidthStretch, 1.2f);
        ImGui::TableSetupColumn("Y", ImGuiTableColumnFlags_WidthStretch, 1.2f);
        ImGui::TableSetupColumn("KB", ImGuiTableColumnFlags_WidthStretch, 0.7f);
        ImGui::TableSetupColumn("TD (MD)", ImGuiTableColumnFlags_WidthStretch, 0.8f);
        ImGui::TableSetupColumn("survey", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
        ImGui::TableSetupColumn("tops", ImGuiTableColumnFlags_WidthFixed, fs * 2.5f);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.6f);
        ImGui::TableHeadersRow();
        int del = -1;
        for (int i = 0; i < int(S.wells.size()); ++i) {
            GeoWell& w = S.wells[size_t(i)];
            ImGui::TableNextRow(); ImGui::PushID(i);
            if (i == E.sel) ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg1, (Pal::Accent & 0x00FFFFFF) | 0x50000000);
            ImGui::TableNextColumn();
            char nb[96]; std::snprintf(nb, sizeof nb, "%s", w.name.c_str());
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputText("##n", nb, sizeof nb)) w.name = nb;
            if (ImGui::IsItemActivated()) E.sel = i;
            const double ox = w.x, oy = w.y;
            bool moved = false;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##x", &w.x, 0, 0, "%.2f"); moved |= ImGui::IsItemDeactivatedAfterEdit(); if (ImGui::IsItemActivated()) E.sel = i;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##y", &w.y, 0, 0, "%.2f"); moved |= ImGui::IsItemDeactivatedAfterEdit(); if (ImGui::IsItemActivated()) E.sel = i;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##kb", &w.kb, 0, 0, "%.2f"); moved |= ImGui::IsItemDeactivatedAfterEdit(); if (ImGui::IsItemActivated()) E.sel = i;
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputDouble("##td", &w.td, 0, 0, "%.1f") && w.traj.size() <= 2 && w.trajSrc == 0) {}
            if (ImGui::IsItemDeactivatedAfterEdit() && w.td > 0 && w.traj.size() <= 2) setVertical(w, w.td);   // vertical wells follow TD
            if (ImGui::IsItemActivated()) E.sel = i;
            if (moved) headEdited(w, ox, oy);
            ImGui::TableNextColumn();
            if (w.traj.size() > 2) ImGui::Text("%zu st.", w.traj.size()); else if (w.traj.size() == 2) ImGui::TextDisabled("vert."); else ImGui::TextDisabled("-");
            ImGui::TableNextColumn(); ImGui::Text("%zu", w.tops.size());
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_TRASH_2)) del = i;
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("delete %s (with its survey, logs and tops)", w.name.c_str());
            ImGui::PopID();
        }
        ImGui::EndTable();
        if (del >= 0) { S.wells.erase(S.wells.begin() + del); if (E.sel >= int(S.wells.size())) E.sel = int(S.wells.size()) - 1; }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    // ---------------- right: selected well
    ImGui::BeginChild("##wr", ImVec2(0, 0), ImGuiChildFlags_Borders);
    if (E.sel < 0 || E.sel >= int(S.wells.size())) { ImGui::TextDisabled("Add or select a well."); ImGui::EndChild(); ImGui::End(); return; }
    GeoWell& w = S.wells[size_t(E.sel)];
    ImGui::Text(ICON_ROUTE " %s", w.name.c_str()); ImGui::SameLine();
    ImGui::TextDisabled("head %.1f, %.1f  KB %.1f  |  %s", w.x, w.y, w.kb, w.traj.empty() ? "no survey (vertical)" : w.trajSrc == 0 ? "MD/INC/AZI survey" : "X/Y/TVDSS survey");
    if (ImGui::BeginTabBar("##wtabs")) {
        if (ImGui::BeginTabItem(TR("Trajectory"))) {
            ImGui::RadioButton(TR("MD / INC / AZI"), &E.mode, 0); ImGui::SameLine(); ImGui::RadioButton(TR("MD / X / Y / TVDSS"), &E.mode, 1);
            ui::help("Type or paste the survey. MD / INC / AZI: positions by minimum curvature from the head (they follow head edits). "
                     "X / Y / TVDSS: absolute stations; MD left at -1 is computed from the path length. Changes apply when a cell loses focus.");
            if (ImGui::Button(TR(ICON_CLIPBOARD_PASTE " Paste survey"))) {
                const char* c = ImGui::GetClipboardText();
                if (c && *c) {
                    TableFormat f; std::string e;
                    TextTable t = parseTextTable(c, f, e);
                    ColumnMap m = guessColumns(t, f, ImportKind::Deviation);
                    if (validateMap(m, ImportKind::Deviation).empty() && m.find(ColRole::WellName) < 0) {
                        ImportOptions o; o.defaultWell = w.name;
                        ImportResult r = buildImport(t, f, m, ImportKind::Deviation, o);
                        if (!r.wells.empty()) { r.wells.resize(1); r.wells[0].name = w.name; MergeStats st = mergeSurveys(S.wells, r); S.msg(w.name + ": survey pasted, " + std::to_string(r.rowsUsed) + " stations"); for (auto& x : st.warnings) S.msg("  " + x); }
                    } else geoOpenImportText(S, ImportKind::Deviation, c, w.name);   // let the user assign the columns
                }
            }
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_FILE " Import survey file..."))) geoOpenImport(S, ImportKind::Deviation, w.name);
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_COPY " Copy"))) ImGui::SetClipboardText(writeDeviationTable(w).c_str());
            ImGui::SameLine();
            static double vtd = 2000;
            if (ImGui::Button(TR("Vertical to TD"))) setVertical(w, w.td > 0 ? w.td : vtd);
            ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); if (w.td > 0) ImGui::TextDisabled("TD %.0f", w.td); else ImGui::InputDouble("##vtd", &vtd, 0, 0, "%.0f");
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_X " Clear"))) { w.traj.clear(); w.trajSrc = 0; }
            // table
            const bool xyz = E.mode == 1;
            const float tabH = ImGui::GetContentRegionAvail().y * 0.52f;
            bool commit = false; int delRow = -1, insRow = -1;
            if (ImGui::BeginTable("##st", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, tabH))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                const char* h0[] = {"MD", "INC", "AZI", "TVDSS", "X", "Y", ""};
                const char* h1[] = {"MD (-1 auto)", "X", "Y", "TVDSS", "INC", "AZI", ""};
                for (int c = 0; c < 7; ++c) ImGui::TableSetupColumn((xyz ? h1 : h0)[c], c == 6 ? ImGuiTableColumnFlags_WidthFixed : ImGuiTableColumnFlags_WidthStretch, c == 6 ? fs * 3.2f : 1.0f);
                ImGui::TableHeadersRow();
                for (int k = 0; k < int(w.traj.size()); ++k) {
                    Station& s = w.traj[size_t(k)];
                    ImGui::TableNextRow(); ImGui::PushID(k);
                    double* ed[3] = {&s.inc, &s.azi, nullptr};
                    double* ex[3] = {&s.x, &s.y, &s.tvd};
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##md", &s.md, 0, 0, "%.2f"); commit |= ImGui::IsItemDeactivatedAfterEdit();
                    if (!xyz) {
                        for (int c = 0; c < 2; ++c) { ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(c); ImGui::InputDouble("##a", ed[c], 0, 0, "%.3f"); commit |= ImGui::IsItemDeactivatedAfterEdit(); ImGui::PopID(); }
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", s.tvd);
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", s.x);
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", s.y);
                    } else {
                        for (int c = 0; c < 3; ++c) { ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(c); ImGui::InputDouble("##a", ex[c], 0, 0, "%.2f"); commit |= ImGui::IsItemDeactivatedAfterEdit(); ImGui::PopID(); }
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", s.inc);
                        ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", s.azi);
                    }
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton(ICON_PLUS)) insRow = k;
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("insert a station below");
                    ImGui::SameLine(); if (ImGui::SmallButton(ICON_X)) delRow = k;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (ImGui::Button(TR(ICON_PLUS " Add station"))) {
                Station s; if (!w.traj.empty()) { s = w.traj.back(); s.md += 100; if (xyz) { s.tvd += 100; s.md = -1; } }
                w.traj.push_back(s); commit = true;
            }
            if (insRow >= 0) {
                Station s = w.traj[size_t(insRow)];
                if (insRow + 1 < int(w.traj.size())) { const Station& b = w.traj[size_t(insRow + 1)]; s.md = 0.5 * (s.md + b.md); s.inc = 0.5 * (s.inc + b.inc); s.x = 0.5 * (s.x + b.x); s.y = 0.5 * (s.y + b.y); s.tvd = 0.5 * (s.tvd + b.tvd); }
                else s.md += 100;
                w.traj.insert(w.traj.begin() + insRow + 1, s); commit = true;
            }
            if (delRow >= 0) { w.traj.erase(w.traj.begin() + delRow); commit = true; }
            if (commit && !w.traj.empty()) {
                E.err.clear();
                for (auto& s : w.traj) { if (!xyz) { s.inc = std::clamp(s.inc, 0.0, 180.0); s.azi = std::fmod(std::fmod(s.azi, 360.0) + 360.0, 360.0); } }
                if (!xyz) { w.trajSrc = 0; minimumCurvature(w); }
                else { std::vector<Station> st = w.traj; if (!trajFromXyz(w, st, E.err)) w.traj = st; }
                if (!w.traj.empty()) w.td = std::max(w.td, w.traj.back().md);
            }
            if (!E.err.empty()) ImGui::TextColored(col4(Pal::Err), "%s", E.err.c_str());
            // plan + section preview
            const float ph = ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y;
            const float pw = ImGui::GetContentRegionAvail().x * 0.5f - ImGui::GetStyle().ItemSpacing.x;
            std::vector<double> xs, ys, hd, zs;
            if (w.traj.empty()) { xs = {w.x}; ys = {w.y}; hd = {0, 0}; zs = {-w.kb, (w.td > 0 ? w.td : 1000) - w.kb}; }
            else for (auto& s : w.traj) { xs.push_back(s.x); ys.push_back(s.y); hd.push_back(std::hypot(s.x - w.x, s.y - w.y)); zs.push_back(s.tvd); }
            if (ph > fs * 5 && ImPlot::BeginPlot("plan##wp", ImVec2(pw, ph), ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxes("x (east)", "y (north)");
                for (const auto& o : S.wells) if (&o != &w) ImPlot::PlotScatter("##o", &o.x, &o.y, 1, mspec(ImPlotMarker_Circle, 3, ImVec4(0.6f, 0.6f, 0.6f, 1)));
                ImPlot::PlotLine("##p", xs.data(), ys.data(), int(xs.size()), lspec(col4(Pal::Accent), 2.5f));
                ImPlot::PlotScatter("##h", &w.x, &w.y, 1, mspec(ImPlotMarker_Circle, 6, ImVec4(1, 1, 1, 1)));
                ImPlot::EndPlot();
            }
            ImGui::SameLine();
            if (ph > fs * 5 && ImPlot::BeginPlot("section##ws", ImVec2(-1, ph), ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxes("horizontal departure", "TVDSS", 0, ImPlotAxisFlags_Invert);
                ImPlot::PlotLine("##s", hd.data(), zs.data(), int(std::min(hd.size(), zs.size())), lspec(col4(Pal::Accent), 2.5f));
                for (const auto& t : w.tops) { const P3 p = w.at(t.md); const double d = std::hypot(p[0] - w.x, p[1] - w.y), z = p[2]; ImPlot::PlotScatter("##t", &d, &z, 1, mspec(ImPlotMarker_Diamond, 6, ImVec4(1, 0.8f, 0.2f, 1))); ImPlot::PlotText(t.surface.c_str(), d, z, ImVec2(fs * 2, 0)); }
                ImPlot::EndPlot();
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(TR("Tops"))) {
            if (ImGui::Button(TR(ICON_PLUS " Add top"))) { Top t; t.surface = "Top" + std::to_string(w.tops.size() + 1); t.md = w.tops.empty() ? 1000 : w.tops.back().md + 50; w.tops.push_back(t); }
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_CLIPBOARD_PASTE " Paste tops"))) { const char* c = ImGui::GetClipboardText(); if (c && *c) geoOpenImportText(S, ImportKind::Tops, c, w.name); }
            ImGui::SameLine();
            if (ImGui::Button(TR(ICON_FILE " Import tops file..."))) geoOpenImport(S, ImportKind::Tops, w.name);
            int del = -1; bool resort = false;
            if (ImGui::BeginTable("##tt", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY)) {
                ImGui::TableSetupColumn("surface"); ImGui::TableSetupColumn("MD"); ImGui::TableSetupColumn("TVDSS"); ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.6f);
                ImGui::TableHeadersRow();
                for (int k = 0; k < int(w.tops.size()); ++k) {
                    Top& t = w.tops[size_t(k)];
                    ImGui::TableNextRow(); ImGui::PushID(k);
                    ImGui::TableNextColumn(); char nb[96]; std::snprintf(nb, sizeof nb, "%s", t.surface.c_str()); ImGui::SetNextItemWidth(-FLT_MIN); if (ImGui::InputText("##s", nb, sizeof nb)) t.surface = nb;
                    ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::InputDouble("##m", &t.md, 0, 0, "%.2f"); resort |= ImGui::IsItemDeactivatedAfterEdit();
                    ImGui::TableNextColumn(); ImGui::TextDisabled("%.2f", w.at(t.md)[2]);
                    ImGui::TableNextColumn(); if (ImGui::SmallButton(ICON_X)) del = k;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (del >= 0) w.tops.erase(w.tops.begin() + del);
            if (resort) std::stable_sort(w.tops.begin(), w.tops.end(), [](const Top& a, const Top& b) { return a.md < b.md; });
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(TR("Logs"))) {
            if (w.logs.empty()) ImGui::TextDisabled("No logs - import LAS files in the Data step (several at once; matched by the ~Well name or the file name).");
            for (size_t k = 0; k < w.logs.size(); ++k) {
                const WellLog& l = w.logs[k];
                std::string cs; for (auto& c : l.curves) cs += c.mnem + " ";
                ImGui::BulletText("%s: %zu samples, MD %.1f .. %.1f; %s", l.file.empty() ? "log" : stemOf(l.file).c_str(), l.md.size(), l.md.empty() ? 0.0 : l.md.front(), l.md.empty() ? 0.0 : l.md.back(), cs.c_str());
                ImGui::SameLine(); ImGui::PushID(int(k));
                if (ImGui::SmallButton(ICON_X)) { w.logs.erase(w.logs.begin() + long(k)); ImGui::PopID(); break; }
                ImGui::PopID();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::EndChild();
    ImGui::End();
}

// ================================================================ 2-D seismic
void geoAddLines(App::GeoState& S, std::vector<SeismicLine> lines) {
    for (auto& l : lines) {
        std::string nm = l.name; int k = 2;
        auto taken = [&](const std::string& n) { for (auto& x : S.lines) if (x.name == n) return true; return false; };
        while (taken(nm)) nm = l.name + " (" + std::to_string(k++) + ")";
        l.name = nm;
        S.lines.push_back(std::move(l));
    }
    ensureLineHz(S);
    if (S.curLine < 0 && !S.hasCube && !S.lines.empty()) S.curLine = 0;
    ++S.lineRev;
}
void geoSegyBytesUi(App::GeoState& S) {
    if (!ImGui::TreeNode(TR("SEG-Y trace-header byte locations"))) return;
    SegyByteMap& m = S.segyMap;
    if (ui::beginProps("##segyb")) {
        ui::propI("inline", &m.inlineByte, 1, 237, "3-D inline number (rev 1: 189; also common: 9, 17, 221)");
        ui::propI("crossline", &m.xlineByte, 1, 237, "3-D crossline number (rev 1: 193; also common: 13, 21, 225)");
        ui::propI("CDP / trace number", &m.cdpByte, 1, 237, "2-D trace (CDP) number (rev 1: 21)");
        ui::propI("CDP X", &m.xByte, 1, 237, "X coordinate (rev 1: 181; source X: 73; group X: 81)");
        ui::propI("CDP Y", &m.yByte, 1, 237, "Y coordinate (rev 1: 185; source Y: 77; group Y: 85)");
        ui::propI("coordinate scalar", &m.scalarByte, 1, 237, "rev 1: 71 (negative = divide)");
        ui::propBool("source X/Y when zero", &m.fallbackXY);
        ui::endProps();
    }
    if (ImGui::SmallButton(TR("SEG-Y rev 1 defaults"))) m = SegyByteMap();
    ui::hint("Used by the 2-D line reader, the in-memory and the out-of-core 3-D readers. Files are detected as 2-D or 3-D from the inline / crossline values.");
    ImGui::TreePop();
}

void geoLinesPanel(App::GeoState& S) {
    const float fs = ImGui::GetFontSize();
    ensureLineHz(S);
    ImGui::SeparatorText(TR(ICON_ACTIVITY " 2-D lines"));
    if (S.lines.empty()) { ImGui::TextDisabled("No 2-D lines - Data step: SEG-Y (several files at once)."); return; }
    // lines
    int delLine = -1;
    for (int l = 0; l < int(S.lines.size()); ++l) {
        ImGui::PushID(l);
        bool on = S.lineOn[size_t(l)] != 0;
        if (ImGui::Checkbox("##on", &on)) { S.lineOn[size_t(l)] = on; ++S.lineRev; }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("show in the 3-D scene and on the map");
        ImGui::SameLine();
        const SeismicLine& L = S.lines[size_t(l)];
        char b[200]; std::snprintf(b, sizeof b, "%s  (%d tr, %.1f km, %d smp)", L.name.c_str(), L.ntr(), L.length() / 1000.0, L.ns);
        if (ImGui::Selectable(b, S.curLine == l, 0, ImVec2(ImGui::GetContentRegionAvail().x - fs * 1.8f, 0))) { S.curLine = l; S.view = 1; }
        ImGui::SameLine(); if (ImGui::SmallButton(ICON_X)) delLine = l;
        ImGui::PopID();
    }
    if (delLine >= 0) {
        S.lines.erase(S.lines.begin() + delLine); S.lineOn.erase(S.lineOn.begin() + delLine);
        for (auto& h : S.lineHz) if (delLine < int(h.t.size())) h.t.erase(h.t.begin() + delLine);
        if (S.curLine >= int(S.lines.size())) S.curLine = int(S.lines.size()) - 1;
        ++S.lineRev;
        if (S.lines.empty()) return;
    }
    if (S.hasCube) { if (ImGui::RadioButton(TR("section: 3-D volume"), S.curLine < 0)) S.curLine = -1; }
    // horizons
    ImGui::SeparatorText(TR("Line horizons"));
    for (int h = 0; h < int(S.lineHz.size()); ++h) {
        ImGui::PushID(h);
        int n = 0, tot = 0; for (auto& v : S.lineHz[size_t(h)].t) for (double t : v) { ++tot; n += std::isfinite(t); }
        ImGui::ColorButton("##c", hzColour(h), ImGuiColorEditFlags_NoTooltip, ImVec2(fs * 0.8f, fs * 0.8f)); ImGui::SameLine();
        if (ImGui::RadioButton("##cur", S.curLineHz == h)) S.curLineHz = h;
        ImGui::SameLine();
        char nb[64]; std::snprintf(nb, sizeof nb, "%s", S.lineHz[size_t(h)].name.c_str());
        ImGui::SetNextItemWidth(fs * 7); if (ImGui::InputText("##nm", nb, sizeof nb)) S.lineHz[size_t(h)].name = nb;
        ImGui::SameLine(); ImGui::TextDisabled("%d / %d", n, tot);
        ImGui::SameLine();
        if (S.lineHz.size() > 1 && ImGui::SmallButton(ICON_TRASH_2)) { S.lineHz.erase(S.lineHz.begin() + h); S.curLineHz = 0; ++S.lineRev; ImGui::PopID(); break; }
        ImGui::PopID();
    }
    if (ImGui::SmallButton(TR(ICON_PLUS " horizon"))) { LineHorizon h; h.name = "H" + std::to_string(S.lineHz.size() + 1); S.lineHz.push_back(h); ensureLineHz(S); S.curLineHz = int(S.lineHz.size()) - 1; }
    ImGui::SameLine();
    if (ImGui::SmallButton(TR("clear picks on this line")) && S.curLine >= 0) { auto& v = S.lineHz[size_t(S.curLineHz)].t[size_t(S.curLine)]; std::fill(v.begin(), v.end(), std::nan("")); ++S.lineRev; }
    // picking
    ImGui::TextUnformatted("click on the section:");
    ImGui::RadioButton(TR("off"), &S.linePick, 0); ImGui::SameLine();
    ImGui::RadioButton(TR("seed + track"), &S.linePick, 1); ImGui::SameLine();
    ImGui::RadioButton(TR("manual"), &S.linePick, 2); ImGui::SameLine();
    ImGui::RadioButton(TR("erase"), &S.linePick, 3);
    ui::help("Seed + track: a click snaps to the event (track event / window below) and follows it both ways along the line until the jump limit or the event fades.\n"
             "Manual: click points, then right-click (or 'finish') to interpolate between them.\nErase: click / drag removes picks within 5 traces.");
    if (S.linePick == 2 && !S.manualPts.empty()) {
        ImGui::SameLine();
        if (ImGui::SmallButton(TR("finish")) && S.curLine >= 0) {
            const auto t = interpolatePicks(S.lines[size_t(S.curLine)].ntr(), S.manualPts);
            auto& v = S.lineHz[size_t(S.curLineHz)].t[size_t(S.curLine)];
            for (size_t k = 0; k < v.size() && k < t.size(); ++k) if (std::isfinite(t[k])) v[k] = t[k];
            S.manualPts.clear(); ++S.lineRev;
        }
    }
    if (ui::beginProps("##l2p")) {
        static const char* modes[] = {"peak", "trough", "zero crossing + to -"};
        ui::propCombo("track event", &S.trackMode, modes, 3);
        ui::propD("search window (ms)", &S.window, "%.1f");
        ui::propD("max jump between traces (ms)", &S.maxJump, "%.1f");
        ui::propSliderD("clip percentile", &S.lineClip, 0.8, 1.0, "%.3f");
        if (!S.hasCube) {
            ui::propD("velocity v0", &S.vm.v0, "%.1f", "interval velocity at the datum (length unit / s)");
            ui::propD("gradient k", &S.vm.k, "%.3f", "v(z) = v0 + k (z - datum), 1/s");
            ui::propD("datum depth", &S.vm.datum, "%.1f");
        }
        ui::endProps();
    }
    // mis-ties
    const auto ties = lineIntersections(S.lines);
    if (!ties.empty() && ImGui::TreeNode(TR("Intersections / mis-ties"))) {
        const LineHorizon& H = S.lineHz[size_t(S.curLineHz)];
        double sum = 0; int n = 0;
        if (ImGui::BeginTable("##ties", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("line A / trace"); ImGui::TableSetupColumn("line B / trace"); ImGui::TableSetupColumn((H.name + " A - B (ms)").c_str()); ImGui::TableSetupColumn("");
            ImGui::TableHeadersRow();
            for (size_t k = 0; k < ties.size(); ++k) {
                const LineTie& T = ties[k];
                const double ta = H.t[size_t(T.a)][size_t(T.ta)], tb = H.t[size_t(T.b)][size_t(T.tb)];
                ImGui::TableNextRow(); ImGui::PushID(int(k));
                ImGui::TableNextColumn(); ImGui::Text("%s / %d", S.lines[size_t(T.a)].name.c_str(), S.lines[size_t(T.a)].cdp[size_t(T.ta)]);
                ImGui::TableNextColumn(); ImGui::Text("%s / %d", S.lines[size_t(T.b)].name.c_str(), S.lines[size_t(T.b)].cdp[size_t(T.tb)]);
                ImGui::TableNextColumn();
                if (std::isfinite(ta) && std::isfinite(tb)) { const double d = ta - tb; sum += std::fabs(d); ++n; ImGui::TextColored(col4(std::fabs(d) > 2 * S.lines[size_t(T.a)].dt ? Pal::Warn : Pal::Ok), "%+.2f", d); }
                else ImGui::TextDisabled("-");
                ImGui::TableNextColumn(); if (ImGui::SmallButton(TR("go"))) { S.curLine = T.a; S.view = 1; }
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (n) ImGui::Text("mean |mis-tie| %.2f ms over %d intersection(s)", sum / n, n);
        ImGui::TreePop();
    }
    // gridding
    if (ui::accentButton(ICON_ARROW_DOWN_TO_LINE " Grid the line horizon -> depth surface")) {
        const LineHorizon& H = S.lineHz[size_t(S.curLineHz)];
        const auto pts = lineHorizonPoints(S.lines, H, S.vm, 2);
        if (pts.size() < 3) S.msg("Error: pick the horizon on the lines first (fewer than 3 picked traces)");
        else {
            Surface lat = S.hz.empty() ? latticeFor(pts, 120) : S.hz[0];
            double x0 = 1e300, x1 = -1e300; for (auto& p : pts) { x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]); }
            geostat::Variogram vg; vg.model = geostat::VarioModel::Gaussian; vg.rangeMax = vg.rangeMin = std::max(1.0, 0.35 * (x1 - x0));
            Surface s = gridPoints(pts, lat, vg); s.name = H.name;
            bool rep = false; for (auto& h : S.hz) if (h.name == s.name) { S.pushUndo(); h = s; rep = true; }
            if (!rep) S.hz.push_back(s);
            ++S.hzRev;
            S.msg("surface " + s.name + (rep ? " replaced" : " added") + " from " + std::to_string(pts.size()) + " line picks (depth by the velocity model, kriged on " + (S.hz.size() > 1 || rep ? "the first horizon's" : "an automatic") + " lattice)");
        }
    }
    ui::hint("Picks are converted to depth with the velocity model above and kriged; the surface joins the horizon list (Structure step, 3-D scene).");
}

void geoLineSection(App::GeoState& S, ImVec2 size) {
    ensureLineHz(S);
    if (S.curLine < 0 || S.curLine >= int(S.lines.size())) { ImGui::TextDisabled("No 2-D line selected."); return; }
    const float fs = ImGui::GetFontSize();
    const SeismicLine& L = S.lines[size_t(S.curLine)];
    std::vector<std::string> names; for (auto& l : S.lines) names.push_back(l.name);
    int cl = S.curLine;
    ImGui::SetNextItemWidth(fs * 12);
    if (ImGui::Combo("##line", &cl, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &names, int(names.size()))) { S.curLine = cl; S.manualPts.clear(); }
    ImGui::SameLine(); if (ImGui::ArrowButton("##prev", ImGuiDir_Left)) { S.curLine = (S.curLine + int(S.lines.size()) - 1) % int(S.lines.size()); S.manualPts.clear(); }
    ImGui::SameLine(); if (ImGui::ArrowButton("##next", ImGuiDir_Right)) { S.curLine = (S.curLine + 1) % int(S.lines.size()); S.manualPts.clear(); }
    if (S.hasCube) { ImGui::SameLine(); if (ImGui::SmallButton(TR("3-D volume"))) { S.curLine = -1; return; } }
    ImGui::SameLine();
    std::vector<std::string> hn; for (auto& h : S.lineHz) hn.push_back(h.name);
    ImGui::SetNextItemWidth(fs * 7); ImGui::Combo("horizon##lh", &S.curLineHz, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(size_t(i)).c_str(); }, &hn, int(hn.size()));
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9); ImGui::Combo("##lpm", &S.linePick, "Navigate\0Seed + track\0Manual picks\0Erase\0Fault sticks\0");
    if (S.linePick == 4) {   // C44: fault sticks on 2-D lines (same draft as the 3-D sections)
        ImGui::SameLine(); ImGui::Text("| stick %zu pts, fault %zu sticks", S.stick.size(), S.draft.sticks.size());
        auto finish = [&] {
            if (S.stick.size() >= 2) {
                std::sort(S.stick.begin(), S.stick.end(), [](const P3& a, const P3& b) { return a[2] < b[2]; });
                S.draft.sticks.push_back(S.stick);
                S.msg("fault stick " + std::to_string(S.draft.sticks.size()) + ": " + std::to_string(S.stick.size()) + " points");
            }
            S.stick.clear();
        };
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Finish stick"))) finish();
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Undo"))) { if (!S.stick.empty()) S.stick.pop_back(); else if (!S.draft.sticks.empty()) S.draft.sticks.pop_back(); }
        ImGui::SameLine();
        if (ImGui::SmallButton(TR("Save fault"))) {
            finish();
            if (S.draft.sticks.size() >= 2) {
                Fault f = S.draft;
                if (f.name.empty()) f.name = "F" + std::to_string(S.faults.size() + 1);
                S.faults.push_back(f); S.msg("fault " + f.name + " saved: " + std::to_string(f.sticks.size()) + " sticks (rebuild the grid to use it)");
                S.draft = Fault();
            } else S.msg("a fault needs at least two sticks (digitise them on two or more lines)");
        }
        ImGui::SameLine(); if (ImGui::SmallButton(TR("Discard"))) { S.stick.clear(); S.draft = Fault(); }
        if (!ImGui::GetIO().WantTextInput) {
            if (ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter)) finish();
            if (ImGui::IsKeyPressed(ImGuiKey_Backspace) && !S.stick.empty()) S.stick.pop_back();
        }
    }
    ImGui::SameLine(); ImGui::TextDisabled("%d traces, CDP %d..%d, %.2f km", L.ntr(), L.cdp.empty() ? 0 : L.cdp.front(), L.cdp.empty() ? 0 : L.cdp.back(), L.length() / 1000);
    if (L.ntr() < 1 || L.ns < 2) return;
    // amplitude range (cached per line / clip)
    static std::string ck; static float lo = -1, hi = 1;
    char kb[128]; std::snprintf(kb, sizeof kb, "%p|%zu|%.4f", (const void*)L.a.data(), L.a.size(), S.lineClip);
    if (ck != kb) {
        ck = kb;
        std::vector<float> a; a.reserve(L.a.size() / 7 + 1);
        for (size_t k = 0; k < L.a.size(); k += 7) a.push_back(std::fabs(L.a[k]));
        if (!a.empty()) { const size_t q = std::min(a.size() - 1, size_t(S.lineClip * double(a.size() - 1))); std::nth_element(a.begin(), a.begin() + long(q), a.end()); hi = std::max(1e-12f, a[q]); lo = -hi; }
    }
    std::vector<float> img(size_t(L.ns) * size_t(L.ntr()));
    for (int r = 0; r < L.ns; ++r) for (int t = 0; t < L.ntr(); ++t) img[size_t(r) * size_t(L.ntr()) + size_t(t)] = L.at(t, L.ns - 1 - r);
    const double tEnd = L.t0 + (L.ns - 1) * L.dt;
    const auto ties = lineIntersections(S.lines);
    if (ImPlot::BeginPlot("##l2d", ImVec2(size.x, std::max(100.0f, size.y - ImGui::GetFrameHeightWithSpacing() * 1.5f)))) {
        ImPlot::SetupAxes("trace", "TWT (ms)", 0, ImPlotAxisFlags_Invert);
        ImPlot::PushColormap(ImPlotColormap_RdBu);
        ImPlot::PlotHeatmap("##amp", img.data(), L.ns, L.ntr(), lo, hi, nullptr, ImPlotPoint(0, L.t0), ImPlotPoint(L.ntr() - 1, tEnd));
        ImPlot::PopColormap();
        // horizons
        for (int h = 0; h < int(S.lineHz.size()); ++h) {
            const auto& v = S.lineHz[size_t(h)].t[size_t(S.curLine)];
            std::vector<double> xs, ys;
            auto flush = [&] { if (xs.size() > 0) ImPlot::PlotLine(("##hz" + std::to_string(h) + "_" + std::to_string(xs.front())).c_str(), xs.data(), ys.data(), int(xs.size()), lspec(hzColour(h), h == S.curLineHz ? 2.5f : 1.5f)); xs.clear(); ys.clear(); };
            for (int t = 0; t < int(v.size()); ++t) { if (std::isfinite(v[size_t(t)])) { xs.push_back(t); ys.push_back(v[size_t(t)]); } else flush(); }
            flush();
        }
        // intersections: vertical marker + the crossing line's pick (mis-tie at a glance)
        for (const auto& T : ties) {
            if (T.a != S.curLine && T.b != S.curLine) continue;
            const bool isA = T.a == S.curLine;
            const int me = isA ? T.ta : T.tb, other = isA ? T.b : T.a, ot = isA ? T.tb : T.ta;
            const double x = me;
            ImPlotSpec sp = lspec(ImVec4(1, 1, 1, 0.55f), 1.0f);
            ImPlot::PlotInfLines(("##tie" + std::to_string(other)).c_str(), &x, 1, sp);
            ImPlot::PlotText(S.lines[size_t(other)].name.c_str(), x, L.t0 + 0.03 * (tEnd - L.t0), ImVec2(0, 0));
            for (int h = 0; h < int(S.lineHz.size()); ++h) {
                const double to = S.lineHz[size_t(h)].t[size_t(other)][size_t(ot)];
                if (std::isfinite(to)) ImPlot::PlotScatter("##tiep", &x, &to, 1, mspec(ImPlotMarker_Square, 7, hzColour(h)));
            }
        }
        // C44 faults: cut of each saved fault (stick interpolation), sticks digitised near this line, the draft
        {
            ressim::geo::Polygon path; path.closed = false;
            const auto dist = L.distance();
            for (int t = 0; t < L.ntr(); ++t) path.pts.push_back({L.x[size_t(t)], L.y[size_t(t)]});
            auto trOf = [&](double s) {   // distance along the line -> fractional trace
                const auto it = std::lower_bound(dist.begin(), dist.end(), s);
                if (it == dist.begin()) return 0.0;
                if (it == dist.end()) return double(L.ntr() - 1);
                const size_t k = size_t(it - dist.begin());
                const double d0 = dist[k - 1], d1 = dist[k];
                return double(k - 1) + (d1 > d0 ? (s - d0) / (d1 - d0) : 0.0);
            };
            const double zt0 = timeToDepth(S.vm, L.t0), zt1 = timeToDepth(S.vm, tEnd);
            for (const auto& f : S.faults) {
                const auto br = ressim::geo::faultCutOnPath(f, path, std::min(zt0, zt1), std::max(zt0, zt1), 50);
                for (size_t b = 0; b < br.size(); ++b) {
                    std::vector<double> xs, ys;
                    for (const auto& q : br[b]) { xs.push_back(trOf(q[0])); ys.push_back(depthToTime(S.vm, q[1])); }
                    ImPlot::PlotLine(("##fc" + f.name + std::to_string(b)).c_str(), xs.data(), ys.data(), int(xs.size()), lspec(ImVec4(1, 0.25f, 0.25f, 1), 2.5f));
                    if (!xs.empty()) ImPlot::PlotText(f.name.c_str(), xs.front(), ys.front(), ImVec2(0, -fs * 0.7f));
                }
            }
            const double spacing = L.ntr() > 1 ? L.length() / (L.ntr() - 1) : 1.0;
            auto drawStick = [&](const std::vector<P3>& st, ImVec4 col, const char* id) {   // points within 1.5 trace spacings
                std::vector<double> xs, ys;
                for (const auto& p : st) {
                    int best = -1; double bd = 1e300;
                    for (int t = 0; t < L.ntr(); ++t) { const double d = std::hypot(L.x[size_t(t)] - p[0], L.y[size_t(t)] - p[1]); if (d < bd) { bd = d; best = t; } }
                    if (best < 0 || bd > 1.5 * spacing) continue;
                    xs.push_back(best); ys.push_back(depthToTime(S.vm, p[2]));
                }
                ImPlotSpec sp = lspec(col, 2.5f); sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 4; sp.MarkerFillColor = sp.MarkerLineColor = col;
                if (!xs.empty()) ImPlot::PlotLine(id, xs.data(), ys.data(), int(xs.size()), sp);
            };
            for (const auto& st : S.draft.sticks) drawStick(st, ImVec4(1, 0.75f, 0.1f, 1), "##draft");
            if (!S.stick.empty()) drawStick(S.stick, ImVec4(1, 0.95f, 0.2f, 1), "##stick");
        }
        // manual points being drawn
        if (!S.manualPts.empty()) {
            std::vector<double> xs, ys; for (auto& p : S.manualPts) { xs.push_back(p.first); ys.push_back(p.second); }
            ImPlotSpec sp = lspec(ImVec4(1, 1, 1, 0.9f), 1.5f); sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 4;
            ImPlot::PlotLine("##man", xs.data(), ys.data(), int(xs.size()), sp);
        }
        if (ImPlot::IsPlotHovered()) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            const int tr = std::clamp(int(std::lround(p.x)), 0, L.ntr() - 1);
            const int s = std::clamp(int(std::lround((p.y - L.t0) / L.dt)), 0, L.ns - 1);
            ImGui::SetTooltip("%s  trace %d  CDP %d\nx %.1f  y %.1f\nTWT %.1f ms  depth %.1f\namplitude %.4g", L.name.c_str(), tr, L.cdp[size_t(tr)], L.x[size_t(tr)], L.y[size_t(tr)], p.y, timeToDepth(S.vm, p.y), L.at(tr, s));
            auto& v = S.lineHz[size_t(S.curLineHz)].t[size_t(S.curLine)];
            const bool click = ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::GetIO().KeyCtrl;
            if (S.linePick == 1 && click) {
                const auto t = trackLine(L, {{tr, p.y}}, S.trackMode, S.window, S.maxJump);
                int n = 0; for (size_t k = 0; k < v.size() && k < t.size(); ++k) if (std::isfinite(t[k])) { v[k] = t[k]; ++n; }
                S.msg(S.lineHz[size_t(S.curLineHz)].name + " on " + L.name + ": " + std::to_string(n) + " traces tracked from trace " + std::to_string(tr));
                ++S.lineRev;
            } else if (S.linePick == 2) {
                if (click) S.manualPts.push_back({tr, p.y});
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Right) && !S.manualPts.empty()) {
                    const auto t = interpolatePicks(L.ntr(), S.manualPts);
                    for (size_t k = 0; k < v.size() && k < t.size(); ++k) if (std::isfinite(t[k])) v[k] = t[k];
                    S.manualPts.clear(); ++S.lineRev;
                }
            } else if (S.linePick == 4) {   // C44 fault stick point at the clicked trace / time
                if (click) S.stick.push_back({L.x[size_t(tr)], L.y[size_t(tr)], timeToDepth(S.vm, p.y)});
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Right) && S.stick.size() >= 2) {
                    std::sort(S.stick.begin(), S.stick.end(), [](const P3& a, const P3& b) { return a[2] < b[2]; });
                    S.draft.sticks.push_back(S.stick); S.stick.clear();
                    S.msg("fault stick " + std::to_string(S.draft.sticks.size()) + " on " + L.name);
                }
            } else if (S.linePick == 3 && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                for (int k = std::max(0, tr - 5); k <= std::min(L.ntr() - 1, tr + 5); ++k) v[size_t(k)] = std::nan("");
                ++S.lineRev;
            }
        }
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("amplitude clip +-%.4g (%.1f %%) | %s", hi, S.lineClip * 100, S.linePick == 2 ? "manual: left-click points, right-click to finish" : S.linePick == 1 ? "click a reflector to track it" : S.linePick == 3 ? "drag to erase" : "navigate (wheel zoom, drag pan)");
}
