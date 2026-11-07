// Data-table policy for the whole app (ui::beginTable / ui::tableColumn / ui::endTable).
//
// Every bordered / striped table goes through here so that all of them behave the same:
//  * columns can always be resized by dragging the header borders (double-click a border = auto-fit);
//  * a column is never narrower than a readable numeric field (input cells used to collapse to the width of the
//    header text, so values were cut off);
//  * when all columns fit the available width they stretch to fill it (weights from max(header, minimum width),
//    not from the header text alone); otherwise they keep fixed widths and the table scrolls horizontally;
//  * scrolling tables without an explicit height take the height of their content (capped). Previously they took
//    "the rest of the window", which shrinks to nothing further down a scrolling page: with several reactions the
//    stoichiometry tables became thinner and thinner and finally disappeared.
// Layout-only tables (no borders, no row background) are passed straight to Dear ImGui.
#include "widgets.h"
#include "imgui_internal.h"
#include <vector>
#include <cmath>

namespace ui {
namespace {

struct TableState {
    ImGuiTable* table = nullptr;
    bool data = false;          // data table: policy applies
    bool fixedMode = false;     // fixed column widths + horizontal scroll (else stretch to fill)
    bool autoHeight = false;    // height follows the content (stored per table id in the window storage)
    bool scrollX = false;
    float minW = 0;
    ImGuiID key = 0;
    std::vector<char> bump;     // per column: enforce the minimum width on the auto-fit width
};
std::vector<TableState> g_tables;

bool isDataTable(ImGuiTableFlags f) {
    return (f & (ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) != 0;
}

}  // namespace

float tableMinColumnWidth() { return ImGui::GetFontSize() * 5.5f; }

bool beginTable(const char* id, int columns, ImGuiTableFlags flags, ImVec2 size) {
    TableState st;
    st.data = isDataTable(flags);
    if (!st.data) {
        if (!ImGui::BeginTable(id, columns, flags, size)) return false;
        st.table = ImGui::GetCurrentTable();
        g_tables.push_back(std::move(st));
        return true;
    }
    const ImGuiStyle& sty = ImGui::GetStyle();
    const float fs = ImGui::GetFontSize();
    st.minW = tableMinColumnWidth();
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float availW = size.x > 0 ? size.x : avail.x + (size.x < 0 ? size.x : 0.0f);

    // sizing mode with hysteresis (a parent scrollbar appearing must not flip the mode every frame)
    ImGuiStorage* store = ImGui::GetStateStorage();
    st.key = ImGui::GetID(id);
    const ImGuiID modeKey = st.key + 0x51u, hKey = st.key + 0x52u;
    const float need = float(columns) * (st.minW + 2.0f * sty.CellPadding.x + 1.0f);
    bool fixed = store->GetBool(modeKey, need > availW);
    if (fixed && need < availW - 2.0f * fs) fixed = false;
    else if (!fixed && need > availW) fixed = true;
    store->SetBool(modeKey, fixed);
    st.fixedMode = fixed;

    flags |= ImGuiTableFlags_Resizable;
    flags &= ~ImGuiTableFlags_SizingMask_;
    if (fixed) { flags |= ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollX; flags &= ~ImGuiTableFlags_NoHostExtendX; }
    else { flags |= ImGuiTableFlags_SizingStretchProp; flags &= ~ImGuiTableFlags_ScrollX; }
    st.scrollX = (flags & ImGuiTableFlags_ScrollX) != 0;

    const bool scrolling = (flags & (ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY)) != 0;
    if (scrolling && size.y == 0.0f) {
        // content height of the previous frame (+ horizontal scrollbar), capped at the visible space but never tiny
        st.autoHeight = true;
        const float content = store->GetFloat(hKey, fs * 8.0f);
        const float cap = std::max(avail.y, fs * 22.0f);
        size.y = std::min(content, cap);
        size.y = std::max(size.y, ImGui::GetFrameHeightWithSpacing() * 2.0f);
    }
    if (!ImGui::BeginTable(id, columns, flags, size)) return false;
    st.table = ImGui::GetCurrentTable();
    st.bump.assign(size_t(std::max(columns, 0)), 0);
    g_tables.push_back(std::move(st));
    return true;
}

void tableColumn(const char* label, ImGuiTableColumnFlags flags, float w, ImGuiID userId) {
    label = TRS(label);   // C45 translated header (columns are identified by index, not by label)
    if (g_tables.empty() || g_tables.back().table != ImGui::GetCurrentTable() || !g_tables.back().data) {
        ImGui::TableSetupColumn(label, flags, w, userId);
        return;
    }
    TableState& st = g_tables.back();
    const ImGuiStyle& sty = ImGui::GetStyle();
    const float fs = ImGui::GetFontSize();
    const char* end = label ? ImGui::FindRenderedTextEnd(label) : nullptr;
    const bool emptyHeader = !label || end == label;
    const float hw = emptyHeader ? 0.0f : ImGui::CalcTextSize(label, end).x + fs * 0.9f;   // + sort arrow / slack
    const int col = ImGui::GetCurrentTable()->DeclColumnsCount;   // index of the column being declared
    const ImGuiTableColumnFlags sizing = flags & ImGuiTableColumnFlags_WidthMask_;
    flags &= ~ImGuiTableColumnFlags_WidthMask_;

    if (sizing == ImGuiTableColumnFlags_WidthFixed) {
        // explicit fixed column (row labels, buttons, check boxes): keep the requested width, at least the header
        ImGui::TableSetupColumn(label, flags | ImGuiTableColumnFlags_WidthFixed, w > 0 ? std::max(w, hw) : 0.0f, userId);
        return;
    }
    if (emptyHeader && sizing != ImGuiTableColumnFlags_WidthStretch) {
        // unlabeled column: delete buttons, row selectors - fit the content
        ImGui::TableSetupColumn(label, flags | ImGuiTableColumnFlags_WidthFixed, 0.0f, userId);
        return;
    }
    // explicit stretch weights are relative (0.5 .. 4): scale them to widths; a weightless stretch column is a text column
    float want = std::max(st.minW, hw);
    if (sizing == ImGuiTableColumnFlags_WidthStretch) want = std::max(want, w > 0 ? w * fs * 6.0f : fs * 14.0f);
    if (st.fixedMode) {
        // auto-fit to the content on first use, but never below the minimum width (enforced in endTable)
        if (col >= 0 && col < int(st.bump.size())) st.bump[size_t(col)] = 1;
        ImGui::TableSetupColumn(label, flags | ImGuiTableColumnFlags_WidthFixed,
                                sizing == ImGuiTableColumnFlags_WidthStretch ? want : 0.0f, userId);
        (void)sty;
    } else {
        ImGui::TableSetupColumn(label, flags | ImGuiTableColumnFlags_WidthStretch, want, userId);
    }
}

void endTable() {
    if (g_tables.empty() || g_tables.back().table != ImGui::GetCurrentTable()) { ImGui::EndTable(); return; }
    TableState st = std::move(g_tables.back());
    g_tables.pop_back();
    ImGuiTable* t = st.table;
    if (st.data && st.fixedMode && t) {
        // minimum auto-fit width: input cells sized with -FLT_MIN report no content width of their own
        const ImGuiStyle& sty = ImGui::GetStyle();
        for (int c = 0; c < t->ColumnsCount && c < int(st.bump.size()); ++c) {
            if (!st.bump[size_t(c)]) continue;
            ImGuiTableColumn& col = t->Columns[c];
            const float hdr = col.ContentMaxXHeadersIdeal - col.WorkMinX;
            const float target = col.WorkMinX + std::max(st.minW, hdr);
            col.ContentMaxXUnfrozen = std::max(col.ContentMaxXUnfrozen, target);
        }
        (void)sty;
    }
    ImGui::EndTable();
    if (st.autoHeight && t && t->InnerWindow && t->InnerWindow != t->OuterWindow) {
        const ImGuiWindow* iw = t->InnerWindow;
        float content = iw->DC.CursorMaxPos.y - iw->DC.CursorStartPos.y;
        content += 2.0f * ImGui::GetStyle().CellPadding.y + 3.0f;
        if (st.scrollX) content += ImGui::GetStyle().ScrollbarSize;
        if (content > 1.0f && std::isfinite(content)) ImGui::GetStateStorage()->SetFloat(st.key + 0x52u, content);
    }
}

}  // namespace ui

// ============================================================================ spreadsheet-style data grid
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

namespace ui {
namespace {

struct GridState {
    int ar = -1, ac = 0;        // anchor of the selection
    int cr = -1, cc = 0;        // cursor (other corner)
    bool rowMode = false;       // whole rows selected from the # column
    bool dragging = false;
    int er = -1, ec = -1;       // cell being edited
    char buf[128] = {0};
    bool focusEdit = false, cursorToEnd = false;
    std::string msg; double msgUntil = 0;
};
std::map<ImGuiID, GridState> g_grid;
ImGuiID g_activeGrid = 0;

bool selValid(const GridState& s) { return s.ar >= 0 && s.cr >= 0; }
int sr0(const GridState& s) { return std::min(s.ar, s.cr); }
int sr1(const GridState& s) { return std::max(s.ar, s.cr); }
int sc0(const GridState& s) { return std::min(s.ac, s.cc); }
int sc1(const GridState& s) { return std::max(s.ac, s.cc); }
bool inSel(const GridState& s, int r, int c) { return selValid(s) && r >= sr0(s) && r <= sr1(s) && c >= sc0(s) && c <= sc1(s); }

std::string trim(const std::string& t) {
    size_t a = 0, b = t.size();
    while (a < b && (t[a] == ' ' || t[a] == '\t' || t[a] == '"' || t[a] == '\r')) ++a;
    while (b > a && (t[b - 1] == ' ' || t[b - 1] == '\t' || t[b - 1] == '"' || t[b - 1] == '\r')) --b;
    return t.substr(a, b - a);
}
// number from a pasted token: "1.5", "1,5" (decimal comma), "1e-3", "1 234.5" not supported; false if not numeric
bool toNumber(std::string t, bool commaDecimal, double& v) {
    t = trim(t);
    if (t.empty()) return false;
    if (commaDecimal && t.find(',') != std::string::npos && t.find('.') == std::string::npos) for (char& ch : t) if (ch == ',') ch = '.';
    char* e = nullptr;
    v = std::strtod(t.c_str(), &e);
    return e && e != t.c_str() && *e == 0 && std::isfinite(v);
}

void setMsg(GridState& s, const std::string& m) { s.msg = m; s.msgUntil = ImGui::GetTime() + 4.0; }

bool getNum(GridAccess& g, int r, int c, double& v) {
    if (double* d = g.num ? g.num(r, c) : nullptr) { v = *d; return true; }
    return g.numGet && g.numGet(r, c, v);
}
bool isNumCell(GridAccess& g, int r, int c) { double v; return getNum(g, r, c, v); }
void putNum(GridAccess& g, int r, int c, double v) {
    if (double* d = g.num ? g.num(r, c) : nullptr) { *d = v; return; }
    if (g.numSet) g.numSet(r, c, v);
}
std::string cellText(GridAccess& g, int r, int c, const char* fmt) {
    double v;
    if (getNum(g, r, c, v)) { char b[64]; std::snprintf(b, sizeof b, fmt, v); return b; }
    return g.getText ? g.getText(r, c) : std::string();
}
// returns true when the cell changed
bool setCell(GridAccess& g, int r, int c, const std::string& tok, bool commaDecimal, int& skipped) {
    if (g.readOnly && g.readOnly(r, c)) { ++skipped; return false; }
    double old;
    if (getNum(g, r, c, old)) {
        double v;
        if (!toNumber(tok, commaDecimal, v)) { ++skipped; return false; }
        if (old == v) return false;
        putNum(g, r, c, v); return true;
    }
    if (g.setText) { g.setText(r, c, trim(tok)); return true; }
    ++skipped; return false;
}

std::string copyText(GridAccess& g, const GridState& s, const std::vector<const char*>* headers) {
    std::string out;
    if (headers) {
        for (int c = sc0(s); c <= sc1(s); ++c) { if (c > sc0(s)) out += '\t'; out += (*headers)[size_t(c)]; }
        out += "\r\n";
    }
    for (int r = sr0(s); r <= sr1(s); ++r) {
        for (int c = sc0(s); c <= sc1(s); ++c) { if (c > sc0(s)) out += '\t'; out += cellText(g, r, c, "%.10g"); }
        out += "\r\n";
    }
    return out;
}

}  // namespace

std::vector<std::vector<std::string>> parseTableText(const std::string& text) {
    std::vector<std::string> lines;
    { std::string cur; for (char ch : text) { if (ch == '\n') { lines.push_back(cur); cur.clear(); } else if (ch != '\r') cur += ch; } if (!cur.empty()) lines.push_back(cur); }
    while (!lines.empty() && trim(lines.back()).empty()) lines.pop_back();
    bool tab = false, semi = false, comma = false;
    for (const auto& l : lines) { tab = tab || l.find('\t') != std::string::npos; semi = semi || l.find(';') != std::string::npos; comma = comma || l.find(',') != std::string::npos; }
    std::vector<std::vector<std::string>> rows;
    for (const auto& l : lines) {
        std::vector<std::string> cells;
        if (tab || semi || comma) {
            const char sep = tab ? '\t' : semi ? ';' : ',';
            std::string cur;
            for (char ch : l) { if (ch == sep) { cells.push_back(cur); cur.clear(); } else cur += ch; }
            cells.push_back(cur);
        } else {   // whitespace separated (Notepad, PDF tables, ECLIPSE-style rows)
            std::string cur;
            for (char ch : l) { if (ch == ' ' || ch == '\t') { if (!cur.empty()) cells.push_back(cur); cur.clear(); } else cur += ch; }
            if (!cur.empty()) cells.push_back(cur);
            if (!cells.empty() && cells.back() == "/") cells.pop_back();   // ECLIPSE record terminator
        }
        rows.push_back(std::move(cells));
    }
    return rows;
}

bool gridTable(const char* id, GridAccess& g, const std::vector<const char*>& colNames, float height, int* selected,
               const char* fmt, bool deleteButtons) {
    bool changed = false;
    const int nc = g.cols;
    ImGui::PushID(id);
    const ImGuiID gid = ImGui::GetID("##grid");
    GridState& s = g_grid[gid];
    ImGuiIO& io = ImGui::GetIO();
    const float fs = ImGui::GetFontSize();
    int n = g.count();
    const bool canInsert = bool(g.insert), canErase = bool(g.erase);
    // keep the selection inside the table
    if (selValid(s)) {
        if (n == 0) s.ar = s.cr = -1;
        else { s.ar = std::min(s.ar, n - 1); s.cr = std::min(s.cr, n - 1); s.ac = std::min(s.ac, nc - 1); s.cc = std::min(s.cc, nc - 1); }
    }
    if (s.er >= n || s.ec >= nc) s.er = s.ec = -1;
    const bool commaDec = true;

    // ---- operations (shared by keyboard and context menu)
    auto doCopy = [&](bool headers) {
        if (!selValid(s)) return;
        ImGui::SetClipboardText(copyText(g, s, headers ? &colNames : nullptr).c_str());
        char b[96]; std::snprintf(b, sizeof b, "copied %d x %d cells", sr1(s) - sr0(s) + 1, sc1(s) - sc0(s) + 1); setMsg(s, b);
    };
    auto doClear = [&] {
        if (!selValid(s)) return;
        for (int r = sr0(s); r <= sr1(s); ++r)
            for (int c = sc0(s); c <= sc1(s); ++c) {
                double v;
                if (getNum(g, r, c, v)) { if (v != 0) { putNum(g, r, c, 0.0); changed = true; } }
                else if (g.setText) { g.setText(r, c, ""); changed = true; }
            }
    };
    auto doEraseRows = [&] {
        if (!selValid(s) || !canErase) return;
        const int a = sr0(s), b = sr1(s);
        g.erase(a, b); changed = true; n = g.count();
        s.ar = s.cr = n > 0 ? std::min(a, n - 1) : -1; s.rowMode = false;
        if (selected && *selected >= n) *selected = n - 1;
    };
    auto doInsert = [&](bool below) {
        if (!canInsert) return;
        const int k = selValid(s) ? sr1(s) - sr0(s) + 1 : 1;
        const int at = selValid(s) ? (below ? sr1(s) + 1 : sr0(s)) : n;
        g.insert(at, k); changed = true; n = g.count();
        s.ar = at; s.cr = at + k - 1; s.rowMode = true; s.ac = 0; s.cc = nc - 1;
    };
    auto doPaste = [&](bool replaceAll) {
        const char* clip = ImGui::GetClipboardText();
        if (!clip || !*clip) { setMsg(s, "clipboard is empty"); return; }
        auto rows = parseTableText(clip);
        // skip a header row (no numeric token while the table has numeric columns)
        if (!rows.empty() && (g.num || g.numGet)) {
            bool anyNum = false; double v;
            for (const auto& t : rows[0]) anyNum = anyNum || toNumber(t, commaDec, v);
            if (!anyNum && rows.size() > 1) rows.erase(rows.begin());
        }
        if (rows.empty()) { setMsg(s, "nothing to paste"); return; }
        int R = int(rows.size()), C = 0;
        for (const auto& r : rows) C = std::max(C, int(r.size()));
        int top = 0, left = 0;
        if (replaceAll) {
            if (!canInsert && R != n) { char b[96]; std::snprintf(b, sizeof b, "this table has a fixed size of %d rows", n); setMsg(s, b); R = std::min(R, n); }
            if (canInsert && R > n) { g.insert(n, R - n); n = g.count(); }
            if (canErase && R < n) { g.erase(R, n - 1); n = g.count(); }
            changed = true;
        } else {
            if (!selValid(s)) { s.ar = s.cr = 0; s.ac = s.cc = 0; }
            top = sr0(s); left = sc0(s);
            if (R == 1 && C == 1 && (sr1(s) > sr0(s) || sc1(s) > sc0(s))) {   // one value into the whole selection
                int skipped = 0;
                for (int r = sr0(s); r <= sr1(s); ++r) for (int c = sc0(s); c <= sc1(s); ++c) changed |= setCell(g, r, c, rows[0][0], commaDec, skipped);
                setMsg(s, "filled the selection");
                return;
            }
            if (top + R > n) {
                if (canInsert) { g.insert(n, top + R - n); n = g.count(); changed = true; }
                else R = n - top;
            }
        }
        int skipped = 0, cut = 0;
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < int(rows[size_t(i)].size()); ++j) {
                const int r = top + i, c = left + j;
                if (c >= nc) { ++cut; continue; }
                changed |= setCell(g, r, c, rows[size_t(i)][size_t(j)], commaDec, skipped);
            }
        s.ar = top; s.ac = left; s.cr = top + std::max(0, R - 1); s.cc = std::min(nc - 1, left + std::max(0, C - 1)); s.rowMode = false;
        char b[160];
        std::snprintf(b, sizeof b, "pasted %d x %d%s%s", R, std::min(C, nc - left), skipped ? " - non-numeric values skipped" : "",
                      cut ? " - extra columns ignored" : "");
        setMsg(s, b);
    };
    auto doFillDown = [&] {
        if (!selValid(s) || sr1(s) == sr0(s)) return;
        for (int c = sc0(s); c <= sc1(s); ++c) {
            const std::string v = cellText(g, sr0(s), c, "%.17g");
            int sk = 0;
            for (int r = sr0(s) + 1; r <= sr1(s); ++r) changed |= setCell(g, r, c, v, false, sk);
        }
    };
    auto beginEdit = [&](int r, int c, const char* initial) {
        if (g.readOnly && g.readOnly(r, c)) { setMsg(s, "computed column - read only"); return; }
        if (g.choices && g.choices(c)) initial = nullptr;   // list column: drop-down instead of typing
        s.er = r; s.ec = c; s.focusEdit = true;
        if (initial) { std::snprintf(s.buf, sizeof s.buf, "%s", initial); s.cursorToEnd = true; }
        else { std::snprintf(s.buf, sizeof s.buf, "%s", cellText(g, r, c, "%.10g").c_str()); s.cursorToEnd = false; }
    };
    auto commitEdit = [&] {
        if (s.er < 0) return;
        int sk = 0;
        if (s.er < n && s.ec < nc) changed |= setCell(g, s.er, s.ec, s.buf, commaDec, sk);
        if (sk) setMsg(s, "not a number - value kept");
        s.er = s.ec = -1;
    };
    auto select = [&](int r, int c, bool extend) {
        r = std::max(0, std::min(r, n - 1)); c = std::max(0, std::min(c, nc - 1));
        if (!extend || !selValid(s)) { s.ar = r; s.ac = c; }
        s.cr = r; s.cc = c; s.rowMode = false;
        g_activeGrid = gid;
    };

    // ---- keyboard (this grid has the focus, no cell editor active)
    const bool focusedWin = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);
    if (g_activeGrid == gid && focusedWin && s.er < 0 && !ImGui::IsAnyItemActive() && n > 0) {
        const bool ctrl = io.KeyCtrl, shift = io.KeyShift;
        if (ctrl && ImGui::IsKeyPressed(ImGuiKey_C, false)) doCopy(shift);
        else if (ctrl && ImGui::IsKeyPressed(ImGuiKey_X, false)) { doCopy(false); doClear(); }
        else if (ctrl && ImGui::IsKeyPressed(ImGuiKey_V, false)) doPaste(false);
        else if (ctrl && ImGui::IsKeyPressed(ImGuiKey_A, false)) { s.ar = 0; s.ac = 0; s.cr = n - 1; s.cc = nc - 1; s.rowMode = true; }
        else if (ctrl && ImGui::IsKeyPressed(ImGuiKey_D, false)) doFillDown();
        else if (ImGui::IsKeyPressed(ImGuiKey_Delete, false)) { if (s.rowMode && canErase) doEraseRows(); else doClear(); }
        else if (selValid(s)) {
            if (ImGui::IsKeyPressed(ImGuiKey_UpArrow)) select(s.cr - 1, s.cc, shift);
            else if (ImGui::IsKeyPressed(ImGuiKey_DownArrow)) select(s.cr + 1, s.cc, shift);
            else if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow)) select(s.cr, s.cc - 1, shift);
            else if (ImGui::IsKeyPressed(ImGuiKey_RightArrow) || ImGui::IsKeyPressed(ImGuiKey_Tab)) select(s.cr, s.cc + 1, shift);
            else if (ImGui::IsKeyPressed(ImGuiKey_F2) || ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter)) beginEdit(s.cr, s.cc, nullptr);
            else if (!ctrl && io.InputQueueCharacters.Size > 0) {
                const ImWchar ch = io.InputQueueCharacters[0];
                if (ch > 32 && ch < 127) { char t[2] = {char(ch), 0}; beginEdit(s.cr, s.cc, t); io.InputQueueCharacters.resize(0); }
            }
        }
    }

    // ---- table
    const ImGuiTableFlags fl = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchSame;
    const float rowH = ImGui::GetFrameHeight();
    const float h = height > 0 ? height : std::min(fs * 20.0f, rowH * (float(n) + 1.6f) + 4.0f);
    int delRow = -1;
    bool openCtx = false;
    ImRect outer;
    const bool showDel = deleteButtons && canErase;
    if (ui::beginTable("t", nc + (showDel ? 2 : 1), fl, ImVec2(0, h))) {
        outer = ImGui::GetCurrentTable()->OuterRect;
        ImGui::TableSetupScrollFreeze(1, 1);
        ui::tableColumn("#", ImGuiTableColumnFlags_WidthFixed, fs * 2.2f);
        for (int c = 0; c < nc; ++c) ui::tableColumn(colNames[size_t(c)]);
        if (showDel) ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.6f);
        {   // header row: click a column title to select the whole column (Shift: extend), "#" selects all
            ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
            const int ncols = nc + (showDel ? 2 : 1);
            for (int c = 0; c < ncols; ++c) {
                if (!ImGui::TableSetColumnIndex(c)) continue;
                ImGui::PushID(c);
                ImGui::TableHeader(ImGui::TableGetColumnName(c));
                if (ImGui::IsItemClicked(ImGuiMouseButton_Left) && n > 0) {
                    if (c == 0) { s.ar = 0; s.ac = 0; s.cr = n - 1; s.cc = nc - 1; s.rowMode = true; }
                    else if (c <= nc) {
                        const bool ext = ImGui::GetIO().KeyShift && selValid(s);
                        if (!ext) s.ac = c - 1;
                        s.ar = 0; s.cr = n - 1; s.cc = c - 1; s.rowMode = false;
                    }
                    g_activeGrid = gid;
                }
                if (c >= 1 && c <= nc && ImGui::IsItemHovered()) ImGui::SetTooltip("Click: select the column (Shift: extend), then Ctrl+C / Ctrl+V / Del / Ctrl+D");
                ImGui::PopID();
            }
        }
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImU32 selBg = ImGui::GetColorU32(ImVec4(0.26f, 0.45f, 0.80f, 0.35f));
        const ImU32 curBd = ImGui::GetColorU32(ImVec4(0.40f, 0.65f, 1.00f, 1.0f));
        ImGuiListClipper clip;
        clip.Begin(n, rowH);
        if (selValid(s)) clip.IncludeItemByIndex(s.cr);
        if (s.er >= 0) clip.IncludeItemByIndex(s.er);
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                ImGui::PushID(r);
                ImGui::TableNextRow(0, rowH);
                // row header: selects whole rows, record selection for nested tables
                ImGui::TableSetColumnIndex(0);
                {
                    const std::string lab = g.rowLabel ? g.rowLabel(r) : std::to_string(r + 1);
                    const bool rowSel = s.rowMode && r >= sr0(s) && r <= sr1(s);
                    ImGui::Selectable((lab + "##rh").c_str(), rowSel || (selected && *selected == r), ImGuiSelectableFlags_AllowOverlap, ImVec2(0, rowH - 2));
                    if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(0)) {
                        if (!(io.KeyShift && s.rowMode && selValid(s))) s.ar = r;
                        s.cr = r; s.ac = 0; s.cc = nc - 1; s.rowMode = true; s.dragging = true; g_activeGrid = gid;
                        if (selected) *selected = r;
                        commitEdit();
                    }
                    if (s.dragging && s.rowMode && ImGui::IsMouseDown(0) && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem)) s.cr = r;
                    if (ImGui::IsItemHovered() && ImGui::IsMouseClicked(1)) {
                        if (!(s.rowMode && r >= sr0(s) && r <= sr1(s))) { s.ar = s.cr = r; s.ac = 0; s.cc = nc - 1; s.rowMode = true; }
                        g_activeGrid = gid; openCtx = true;
                    }
                }
                for (int c = 0; c < nc; ++c) {
                    ImGui::TableSetColumnIndex(c + 1);
                    ImGui::PushID(c);
                    const std::vector<std::string>* opts = (s.er == r && s.ec == c && g.choices) ? g.choices(c) : nullptr;
                    if (opts) {
                        // list column: drop-down with the allowed values
                        if (s.focusEdit) { ImGui::OpenPopup("##choice"); s.focusEdit = false; }
                        const std::string cur = cellText(g, r, c, fmt);
                        ImGui::TextUnformatted(cur.c_str());
                        if (ImGui::BeginPopup("##choice")) {
                            for (size_t k = 0; k < opts->size(); ++k)
                                if (ImGui::Selectable((*opts)[k].c_str(), (*opts)[k] == cur)) {
                                    if ((*opts)[k] != cur && g.setText) { g.setText(r, c, (*opts)[k]); changed = true; }
                                    s.er = s.ec = -1;
                                }
                            if (ImGui::IsKeyPressed(ImGuiKey_Escape)) { s.er = s.ec = -1; ImGui::CloseCurrentPopup(); }
                            ImGui::EndPopup();
                        } else s.er = s.ec = -1;
                    } else if (s.er == r && s.ec == c) {
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        if (s.focusEdit) { ImGui::SetKeyboardFocusHere(); s.focusEdit = false; }
                        auto cb = [](ImGuiInputTextCallbackData* d) -> int {
                            GridState* gs = static_cast<GridState*>(d->UserData);
                            if (gs->cursorToEnd) { d->CursorPos = d->SelectionStart = d->SelectionEnd = d->BufTextLen; gs->cursorToEnd = false; }
                            return 0;
                        };
                        const bool enter = ImGui::InputText("##ed", s.buf, sizeof s.buf, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CallbackAlways, cb, &s);
                        const bool esc = ImGui::IsKeyPressed(ImGuiKey_Escape);
                        if (enter) { commitEdit(); select(std::min(r + 1, n - 1), c, false); }
                        else if (esc) { s.er = s.ec = -1; }
                        else if (ImGui::IsItemDeactivated()) {
                            const bool tab = ImGui::IsKeyPressed(ImGuiKey_Tab);
                            commitEdit();
                            if (tab) select(r, std::min(c + 1, nc - 1), false);
                        }
                    } else {
                        const std::string txt = cellText(g, r, c, fmt);
                        const bool sel = inSel(s, r, c);
                        const ImVec2 p0 = ImGui::GetCursorScreenPos();
                        const float cw = ImGui::GetContentRegionAvail().x;
                        ImGui::InvisibleButton("##cell", ImVec2(std::max(cw, 1.0f), rowH - 2));
                        const bool hov = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
                        const ImVec2 p1(p0.x + cw, p0.y + rowH - 2);
                        if (sel) dl->AddRectFilled(ImVec2(p0.x - 2, p0.y - 1), ImVec2(p1.x + 2, p1.y + 1), selBg);
                        if (selValid(s) && s.cr == r && s.cc == c && g_activeGrid == gid) dl->AddRect(ImVec2(p0.x - 2, p0.y - 1), ImVec2(p1.x + 2, p1.y + 1), curBd, 0, 0, 2.0f);
                        // right-aligned numbers, left-aligned text
                        const bool isNum = isNumCell(g, r, c);
                        const float tw = ImGui::CalcTextSize(txt.c_str()).x;
                        const float ty = p0.y + (rowH - 2 - fs) * 0.5f;
                        dl->PushClipRect(p0, p1, true);
                        dl->AddText(ImVec2(isNum ? std::max(p0.x, p1.x - tw - 3) : p0.x + 3, ty), ImGui::GetColorU32(ImGuiCol_Text), txt.c_str());
                        dl->PopClipRect();
                        if (hov && ImGui::IsMouseClicked(0)) {
                            commitEdit();
                            select(r, c, io.KeyShift);
                            s.dragging = true;
                        }
                        if (hov && ImGui::IsMouseDoubleClicked(0)) beginEdit(r, c, nullptr);
                        if (s.dragging && !s.rowMode && ImGui::IsMouseDown(0) && hov) { s.cr = r; s.cc = c; }
                        if (hov && ImGui::IsMouseClicked(1)) {
                            if (!inSel(s, r, c)) select(r, c, false);
                            g_activeGrid = gid; openCtx = true;
                        }
                        if (hov && !ImGui::IsMouseDown(0) && ImGui::GetIO().MouseDelta.x == 0 && txt.size() > 10)
                            ImGui::SetItemTooltip("%s", txt.c_str());
                    }
                    ImGui::PopID();
                }
                if (showDel) {
                    ImGui::TableSetColumnIndex(nc + 1);
                    if (ImGui::SmallButton("\xC3\x97")) delRow = r;
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Delete row");
                }
                ImGui::PopID();
            }
        ui::endTable();
    }
    if (!ImGui::IsMouseDown(0)) s.dragging = false;
    // a click elsewhere ends the grid focus
    if (g_activeGrid == gid && ImGui::IsMouseClicked(0) && !outer.Contains(io.MousePos) && !ImGui::IsPopupOpen("##gridctx")) g_activeGrid = 0;
    if (delRow >= 0 && canErase) {
        g.erase(delRow, delRow); changed = true; n = g.count();
        if (selected && *selected >= n) *selected = n - 1;
        if (selValid(s) && (s.ar >= n || s.cr >= n)) s.ar = s.cr = -1;
    }
    // ---- context menu
    if (openCtx) ImGui::OpenPopup("##gridctx");
    if (ImGui::BeginPopup("##gridctx")) {
        const bool any = selValid(s);
        if (ImGui::MenuItem("Copy", "Ctrl+C", false, any)) doCopy(false);
        if (ImGui::MenuItem("Copy with headers", "Ctrl+Shift+C", false, any)) doCopy(true);
        if (ImGui::MenuItem("Cut", "Ctrl+X", false, any)) { doCopy(false); doClear(); }
        if (ImGui::MenuItem("Paste", "Ctrl+V")) doPaste(false);
        if (ImGui::MenuItem("Replace the whole table with the clipboard")) doPaste(true);
        ImGui::Separator();
        if (ImGui::MenuItem("Fill down", "Ctrl+D", false, any && sr1(s) > sr0(s))) doFillDown();
        if (ImGui::MenuItem("Clear cells", "Del", false, any)) doClear();
        ImGui::Separator();
        if (ImGui::MenuItem("Insert rows above", nullptr, false, canInsert)) doInsert(false);
        if (ImGui::MenuItem("Insert rows below", nullptr, false, canInsert)) doInsert(true);
        if (ImGui::MenuItem("Delete selected rows", s.rowMode ? "Del" : nullptr, false, any && canErase)) doEraseRows();
        ImGui::Separator();
        if (ImGui::MenuItem("Select all", "Ctrl+A", false, n > 0)) { s.ar = 0; s.ac = 0; s.cr = n - 1; s.cc = nc - 1; s.rowMode = true; }
        ImGui::EndPopup();
    }
    // ---- footer: add row, clipboard buttons, status
    if (canInsert) {
        if (ImGui::SmallButton("+ Add row")) { g.insert(n, 1); changed = true; n = g.count(); }
        ImGui::SameLine();
    }
    if (ImGui::SmallButton("Paste")) { g_activeGrid = gid; doPaste(false); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Paste the clipboard at the selected cell (rows are added as needed).\nCopy cells in Excel / LibreOffice / a text editor first.");
    ImGui::SameLine();
    if (ImGui::SmallButton("Copy all")) { s.ar = 0; s.ac = 0; s.cr = n - 1; s.cc = nc - 1; s.rowMode = true; doCopy(true); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Copy the whole table with its headers (tab-separated: pastes into Excel as columns)");
    ImGui::SameLine();
    if (!s.msg.empty() && ImGui::GetTime() < s.msgUntil) ImGui::TextColored(ImVec4(0.55f, 0.80f, 1.0f, 1), "%s", s.msg.c_str());
    else ImGui::TextDisabled("click / drag to select, type to edit, Ctrl+C / Ctrl+V, right-click for more");
    ImGui::PopID();
    return changed;
}

}  // namespace ui
