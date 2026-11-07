// Small UI helpers: property grids, editable tables, std::string input.
#pragma once
#include "icons.h"
#include "imgui.h"
#include "theme.h"
#include "nav.h"
#include "i18n.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <functional>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <cstring>

namespace ui {

// ---- data tables (ui_table.cpp): resizable columns with a readable minimum width, stretch-to-fit or horizontal
// scrolling, and content-following height for scrolling tables. Drop-in replacements for ImGui::BeginTable /
// TableSetupColumn / EndTable; layout tables (no borders / row background) are passed through unchanged.
bool beginTable(const char* id, int columns, ImGuiTableFlags flags = 0, ImVec2 size = ImVec2(0, 0));
void tableColumn(const char* label, ImGuiTableColumnFlags flags = 0, float widthOrWeight = 0.0f, ImGuiID userId = 0);
void endTable();
float tableMinColumnWidth();

// std::string text input (single or multi-line) that grows the string as needed
inline int stringResizeCb(ImGuiInputTextCallbackData* d) {
    if (d->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* s = static_cast<std::string*>(d->UserData);
        s->resize(size_t(d->BufTextLen)); d->Buf = s->data();
    }
    return 0;
}
inline bool inputString(const char* id, std::string& s, ImVec2 size = ImVec2(0, 0), bool multiline = false, ImGuiInputTextFlags fl = 0) {
    fl |= ImGuiInputTextFlags_CallbackResize;
    if (s.capacity() < 16) s.reserve(16);
    return multiline ? ImGui::InputTextMultiline(id, s.data(), s.capacity() + 1, size, fl, stringResizeCb, &s)
                     : ImGui::InputText(id, s.data(), s.capacity() + 1, fl, stringResizeCb, &s);
}

// C41 (G5): section header drawn as a card title - raised surface (elevation 1), accent stripe, bold label
inline void header(const char* text) {
    ImGui::Dummy(ImVec2(0, tok::sp(1)));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float w = std::max(ImGui::GetContentRegionAvail().x, 1.0f);
    if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
    const float th = ImGui::GetTextLineHeight(), h = th + 2.0f * tok::sp(1) + 1.0f;
    const float r = tok::kRadius * uiScale();
    dl->AddRectFilled(p, ImVec2(p.x + w, p.y + h), ImGui::GetColorU32(themeSurface(2)), r);
    dl->AddRectFilled(p, ImVec2(p.x + std::max(3.0f, 3.0f * uiScale()), p.y + h), ImGui::GetColorU32(accentColor(accentId())), r, ImDrawFlags_RoundCornersLeft);
    ImGui::SetCursorScreenPos(ImVec2(p.x + tok::sp(3), p.y + (h - th) * 0.5f));
    ImGui::TextUnformatted(TRS(text));   // C45 translated
    if (g_fonts.bold) ImGui::PopFont();
    ImGui::SetCursorScreenPos(ImVec2(p.x, p.y + h));
    ImGui::Dummy(ImVec2(w, tok::sp(1)));
}

inline void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(TRS(text));
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

// C42 help text style: 0 = compact (an (i) icon after the previous item; hover = tooltip, click = pinned popover with
// Copy), 1 = inline wrapped text (the pre-C42 look). Settings -> Appearance -> "Help text".
inline int& hintMode() { static int m = 0; return m; }
inline void hint(const char* text) {
    if (!text || !*text) return;
    const char* const key = text;   // C45: the ImGui id stays the English text
    text = TRS(text);
    if (hintMode() == 1) {
        ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::PopStyleColor();
        return;
    }
    // same line as the previous item when it fits, else at the start of a new line
    const float iconW = ImGui::GetFontSize() * 1.3f;
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    const ImVec2 last = ImGui::GetItemRectMax();
    if (last.x > ImGui::GetWindowPos().x && last.x + ImGui::GetStyle().ItemSpacing.x + iconW < right) ImGui::SameLine();
    ImGui::PushID(key);
    ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
    ImGui::TextUnformatted(ICON_INFO);
    ImGui::PopStyleColor();
    const bool hov = ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort);
    if (hov) {
        ImGui::GetWindowDrawList()->AddText(ImGui::GetItemRectMin(), ImGui::GetColorU32(col4(Pal::Accent)), ICON_INFO);
        if (ImGui::BeginTooltip()) {
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
            ImGui::TextUnformatted(text);
            ImGui::PopTextWrapPos();
            ImGui::TextDisabled("click to pin");
            ImGui::EndTooltip();
        }
    }
    if (ImGui::IsItemClicked()) ImGui::OpenPopup("##hintpin");
    if (ImGui::BeginPopup("##hintpin")) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 34.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        if (ImGui::SmallButton("Copy")) ImGui::SetClipboardText(text);
        ImGui::EndPopup();
    }
    ImGui::PopID();
}

// ---- C40 (G11) Inspector filter: while active, the standard property rows (prop*) whose label does not contain the
// filter text (case-insensitive) are skipped; they return false (unchanged). Custom rows drawn with ui::label stay visible.
struct PropFilter { bool active = false; std::string text; int shown = 0, hidden = 0; };
inline PropFilter& propFilter() { static PropFilter f; return f; }
inline bool propHidden(const char* text) {
    PropFilter& f = propFilter();
    if (!f.active || f.text.empty() || !text) return false;
    std::string l(text);
    const size_t hh = l.find("##"); if (hh != std::string::npos) l.resize(hh);
    auto low = [](std::string x) { for (char& ch : x) ch = char(std::tolower((unsigned char)ch)); return x; };
    bool hit = low(l).find(low(f.text)) != std::string::npos;
    if (!hit) {   // C45: the filter also matches the translated label
        std::string t(TRS(l.c_str()));
        hit = low(t).find(low(f.text)) != std::string::npos;
    }
    ++(hit ? f.shown : f.hidden);
    return !hit;
}

// ---- C41 "show only modified": lazily captured baseline of every property row. The first time a row is drawn
// after a baseline reset (project open / new / import, or the Inspector's "set baseline" button) its value is
// remembered under (ImGui ID, value address); later frames compare against it. Rows show a modified marker and,
// with onlyModified, unmodified rows are hidden. Limitation: a value changed outside the Inspector (deck editor,
// wizard, script) before its row was first drawn becomes part of the baseline.
struct ModTrack {
    std::unordered_map<unsigned long long, std::string> base;
    bool onlyModified = false, inInspector = false, pending = false;
    int modified = 0, hidden = 0;
};
inline ModTrack& modTrack() { static ModTrack m; return m; }
inline void modBaselineReset() { modTrack().base.clear(); }
inline std::string modVal(double v) { char b[40]; std::snprintf(b, sizeof b, "%.17g", v); return b; }
inline std::string modVal(float v) { return modVal(double(v)); }
inline std::string modVal(int v) { return std::to_string(v); }
inline std::string modVal(bool v) { return v ? "1" : "0"; }
inline std::string modVal(const std::string& v) { return v; }
// true = hide the row; sets ModTrack::pending so label() draws the marker
template <class T> inline bool propModSkip(const char* text, const T* v) {
    ModTrack& m = modTrack();
    m.pending = false;
    if (!m.inInspector || !v) return false;
    const unsigned long long key = (static_cast<unsigned long long>(ImGui::GetID(text)) << 32) ^ static_cast<unsigned long long>(reinterpret_cast<uintptr_t>(v));
    const std::string cur = modVal(*v);
    auto it = m.base.find(key);
    if (it == m.base.end()) { m.base.emplace(key, cur); it = m.base.find(key); }
    const bool mod = it->second != cur;
    if (mod) { ++m.modified; m.pending = true; return false; }
    if (m.onlyModified) { ++m.hidden; return true; }
    return false;
}

// ---- two-column property grid
inline bool beginProps(const char* id, float labelFrac = 0.45f) {
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_Resizable)) return false;   // drag the label / value boundary
    ImGui::TableSetupColumn("k", ImGuiTableColumnFlags_WidthStretch, labelFrac);
    ImGui::TableSetupColumn("v", ImGuiTableColumnFlags_WidthStretch, 1.0f - labelFrac);
    return true;
}
inline void endProps() { ImGui::EndTable(); }
inline void label(const char* text, const char* tip = nullptr) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    const bool modMark = modTrack().pending; modTrack().pending = false;
    ImGui::PushStyleColor(ImGuiCol_Text, modMark ? ImGui::GetStyleColorVec4(ImGuiCol_CheckMark) : themeLabelText());
    if (modMark) { ImGui::TextUnformatted("\xE2\x97\x8F"); ImGui::SameLine(0, ImGui::GetStyle().ItemInnerSpacing.x); }   // filled circle
    {   // C45 translated label (the part before an ImGui "##" suffix only)
        const char* t = TRS(text);
        const char* e = std::strstr(t, "##");
        ImGui::TextUnformatted(t, e);
    }
    ImGui::PopStyleColor();
    if (modMark && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort) && !tip) ImGui::SetTooltip("%s", TR("Modified since the baseline (project open / set baseline)"));
    if (tip && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", TRS(tip));
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
}
inline bool propD(const char* text, double* v, const char* fmt = "%.6g", const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::InputDouble("##v", v, 0, 0, fmt);
    ImGui::PopID();
    return c;
}
inline bool propI(const char* text, int* v, int lo, int hi, const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    bool c = ImGui::InputInt("##v", v, 1, 10);
    if (*v < lo) { *v = lo; } if (*v > hi) { *v = hi; }
    ImGui::PopID();
    return c;
}
inline bool propF(const char* text, float* v, float lo, float hi, const char* fmt = "%.3g", const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::SliderFloat("##v", v, lo, hi, fmt);
    ImGui::PopID();
    return c;
}
inline bool propSliderD(const char* text, double* v, double lo, double hi, const char* fmt = "%.3g", const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    float f = float(*v);
    const bool c = ImGui::SliderFloat("##v", &f, float(lo), float(hi), fmt);
    if (c) *v = f;
    ImGui::PopID();
    return c;
}
inline bool propBool(const char* text, bool* v, const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::Checkbox("##v", v);
    ImGui::PopID();
    return c;
}
inline bool propCombo(const char* text, int* v, const char* const items[], int n, const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::Combo("##v", v, items, n);
    ImGui::PopID();
    return c;
}
inline bool propCombo(const char* text, int* v, const std::vector<std::string>& items, const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    bool c = false;
    const char* cur = *v >= 0 && *v < int(items.size()) ? items[size_t(*v)].c_str() : "";
    if (ImGui::BeginCombo("##v", cur)) {
        for (int i = 0; i < int(items.size()); ++i) if (ImGui::Selectable(items[size_t(i)].c_str(), i == *v)) { *v = i; c = true; }
        ImGui::EndCombo();
    }
    ImGui::PopID();
    return c;
}
inline void propText(const char* text, const char* value) {
    if (propHidden(text)) return;
    label(text);
    ImGui::TextUnformatted(value);
}

// ---- std::string input
inline int strResizeCb(ImGuiInputTextCallbackData* d) {
    if (d->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* s = static_cast<std::string*>(d->UserData);
        s->resize(size_t(d->BufTextLen));
        d->Buf = &(*s)[0];
    }
    return 0;
}
inline bool inputString(const char* id, std::string* s, ImGuiInputTextFlags flags = 0) {
    flags |= ImGuiInputTextFlags_CallbackResize;
    return ImGui::InputText(id, &(*s)[0], s->capacity() + 1, flags, strResizeCb, s);
}
inline bool inputStringHint(const char* id, const char* hint, std::string* s, ImGuiInputTextFlags flags = 0) {
    flags |= ImGuiInputTextFlags_CallbackResize;
    return ImGui::InputTextWithHint(id, hint, &(*s)[0], s->capacity() + 1, flags, strResizeCb, s);
}
inline bool propString(const char* text, std::string* v, const char* tip = nullptr) {
    if (propHidden(text) || propModSkip(text, v)) return false;
    label(text, tip);
    ImGui::PushID(text);
    const bool c = inputString("##v", v);
    ImGui::PopID();
    return c;
}

inline bool accentButton(const char* text, const ImVec2& size);
// ---- C40 empty state: large icon, title and explanation centred in the available area; returns the index of the
// clicked action button (-1 = none). The first action is drawn as the primary (accent) button.
inline int emptyState(const char* icon, const char* title, const char* body, const std::vector<const char*>& actions = {}) {
    title = TRS(title); body = TRS(body);   // C45 translated
    const float fs = ImGui::GetFontSize();
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float w = std::min(avail.x - fs, fs * 30);
    const float x0 = ImGui::GetCursorPosX() + std::max(0.0f, (avail.x - w) * 0.5f);
    ImGui::SetCursorPosY(ImGui::GetCursorPosY() + std::max(fs, avail.y * 0.18f));
    ImFont* f = ImGui::GetFont();
    const ImVec2 isz = f->CalcTextSizeA(fs * 3.0f, FLT_MAX, 0.0f, icon);
    ImGui::SetCursorPosX(x0 + (w - isz.x) * 0.5f);
    const ImVec2 ip = ImGui::GetCursorScreenPos();
    ImGui::GetWindowDrawList()->AddText(f, fs * 3.0f, ip, ImGui::GetColorU32(ImGuiCol_TextDisabled), icon);
    ImGui::Dummy(ImVec2(isz.x, isz.y + fs * 0.4f));
    const ImVec2 tsz = ImGui::CalcTextSize(title);
    ImGui::SetCursorPosX(x0 + std::max(0.0f, (w - tsz.x * 1.15f) * 0.5f));
    if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, fs * 1.15f);
    ImGui::TextUnformatted(title);
    if (g_fonts.bold) ImGui::PopFont();
    ImGui::Spacing();
    ImGui::SetCursorPosX(x0);
    ImGui::PushTextWrapPos(x0 + w);
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s", body);
    ImGui::PopStyleColor();
    ImGui::PopTextWrapPos();
    int clicked = -1;
    if (!actions.empty()) {
        ImGui::Spacing(); ImGui::Spacing();
        float tw = 0;
        for (const char* a : actions) tw += ImGui::CalcTextSize(a).x + ImGui::GetStyle().FramePadding.x * 2 + ImGui::GetStyle().ItemSpacing.x;
        ImGui::SetCursorPosX(x0 + std::max(0.0f, (w - tw) * 0.5f));
        for (int i = 0; i < int(actions.size()); ++i) {
            if (i) ImGui::SameLine();
            ImGui::PushID(i);
            const bool hit = i == 0 ? accentButton(actions[size_t(i)], ImVec2(0, 0)) : ImGui::Button(actions[size_t(i)]);
            if (hit) clicked = i;
            ImGui::PopID();
        }
    }
    return clicked;
}

// ---- accent button
inline bool accentButton(const char* text, const ImVec2& size = ImVec2(0, 0)) {
    // C39: follows the accent colour; white label on all five accents
    const ImVec4 a = accentColor(accentId());
    auto mix = [&](float t, float to) { return ImVec4(a.x + (to - a.x) * t, a.y + (to - a.y) * t, a.z + (to - a.z) * t, 1); };
    ImGui::PushStyleColor(ImGuiCol_Button, mix(0.10f, 0.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, mix(0.12f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, mix(0.28f, 0.0f));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
    const bool r = ImGui::Button(TRL(text), size)   /* C45 translated, English id */;
    ImGui::PopStyleColor(4);
    return r;
}
inline bool dangerButton(const char* text, const ImVec2& size = ImVec2(0, 0)) {
    // C41 (G5): white label on red in every theme (was dark-theme only)
    const bool light = themeIsLight();
    ImGui::PushStyleColor(ImGuiCol_Button, light ? ImVec4(0.80f, 0.20f, 0.22f, 1) : ImVec4(0.55f, 0.16f, 0.18f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, light ? ImVec4(0.88f, 0.28f, 0.30f, 1) : ImVec4(0.70f, 0.22f, 0.24f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, light ? ImVec4(0.68f, 0.14f, 0.16f, 1) : ImVec4(0.45f, 0.12f, 0.14f, 1));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
    const bool r = ImGui::Button(TRL(text), size)   /* C45 translated, English id */;
    ImGui::PopStyleColor(4);
    return r;
}
// C41 (G5): secondary action - transparent with an accent outline (primary = accentButton, destructive = dangerButton)
inline bool secondaryButton(const char* text, const ImVec2& size = ImVec2(0, 0)) {
    const ImVec4 a = accentColor(accentId());
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(a.x, a.y, a.z, 0.18f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(a.x, a.y, a.z, 0.30f));
    ImGui::PushStyleColor(ImGuiCol_Border, a);
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
    const bool r = ImGui::Button(TRL(text), size)   /* C45 translated, English id */;
    ImGui::PopStyleVar();
    ImGui::PopStyleColor(4);
    return r;
}

// ---- spreadsheet-style data grid (ui_table.cpp)
// Excel-like editing for any rectangular table: click / drag / Shift+click / Ctrl+A selection (cells, or whole rows
// from the # column), arrow keys, type or double-click / F2 / Enter to edit (Enter moves down, Tab right, Esc cancels),
// Ctrl+C / Ctrl+X / Ctrl+V (tab-, comma-, semicolon- or space-separated text from Excel, LibreOffice, Notepad, PDF
// tables or another grid; a header row is skipped; rows are appended when the pasted block is longer), Delete
// (clears cells or deletes selected rows), Ctrl+D fill down, right-click menu with insert / delete / duplicate rows
// and "replace the whole table with the clipboard" (the table takes the pasted size).
struct GridAccess {
    std::function<int()> count;                         // number of rows
    int cols = 0;
    std::function<double*(int r, int c)> num;           // numeric cell (nullptr = not numeric)
    std::function<bool(int r, int c, double& v)> numGet;   // alternative: converted numeric value (false = not numeric)
    std::function<void(int r, int c, double v)> numSet;    //   and its setter (validation / unit conversion inside)
    std::function<std::string(int r, int c)> getText;   // optional text cells (num returns nullptr)
    std::function<void(int r, int c, const std::string&)> setText;
    std::function<void(int at, int n)> insert;          // insert n rows before 'at' (empty: fixed row count)
    std::function<void(int r0, int r1)> erase;          // erase rows r0..r1 (empty: rows cannot be deleted)
    std::function<std::string(int r)> rowLabel;         // optional row header text (default 1, 2, ...)
    std::function<const std::vector<std::string>*(int c)> choices;   // C16: list column (edit = drop-down; paste by name / number)
    std::function<bool(int r, int c)> readOnly;         // C16: computed cells (shown, copied, never edited / pasted)
};
bool gridTable(const char* id, GridAccess& g, const std::vector<const char*>& colNames, float height = 0,
               int* selected = nullptr, const char* fmt = "%.6g", bool deleteButtons = true);
// ---- C16: declarative Excel-style table over a vector of records. A new table = a list of column declarations:
//   ui::RecordTable<Aquifer> t;
//   t.text("Name", [](Aquifer& a) -> std::string& { return a.name; })
//    .choice("Type", [](Aquifer& a) -> int& { return a.type; }, {"Fetkovich", "Carter-Tracy"})
//    .num("Volume (m3)", [](Aquifer& a) -> double& { return a.volume; })
//    .num("Depth (ft)", [](Aquifer& a) -> double& { return a.depth; }, 1 / 0.3048)     // display scale (units)
//    .integer("I1", [](Aquifer& a) -> int& { return a.i1; }, 1)                         // display offset (0 -> 1 based)
//    .computed("PV", [](const Aquifer& a) { return a.volume * a.poro; });
//   changed |= t.draw("aq", rows, height);
// All cells get the gridTable behaviour (selection, copy / paste incl. from Excel, fill down, insert / delete rows,
// keyboard navigation); list columns open a drop-down and accept pasted names or numbers.
template <class T>
struct RecordTable {
    struct Col {
        std::string title;
        int kind = 0;                                   // 0 number, 1 integer, 2 text, 3 list, 4 flag, 5 computed, 6 optional number
        std::function<double&(T&)> d;
        std::function<int&(T&)> i;
        std::function<bool&(T&)> b;
        std::function<std::string&(T&)> s;
        std::function<double(const T&)> calc;
        std::function<int(T&)> gi;                      // list columns: index getter / setter (any field type)
        std::function<void(T&, int)> si;
        double scale = 1, offset = 0;
        std::vector<std::string> opts;
        int optBase = 0;
        std::function<bool(T&)> isUnset;                // optional numbers: defaulted state shown as unsetText
        std::function<void(T&)> unset;
        std::string unsetText;
    };
    std::vector<Col> cols;
    bool canInsert = true, canErase = true;
    std::function<T(int at)> newRow;                    // default: copy of the neighbouring row (or T{})
    std::function<void(T&)> onChange;                   // validation / clamping after an edit
    RecordTable& num(const char* t, std::function<double&(T&)> f, double scale = 1) { Col c; c.title = t; c.kind = 0; c.d = f; c.scale = scale; cols.push_back(c); return *this; }
    RecordTable& integer(const char* t, std::function<int&(T&)> f, int offset = 0) { Col c; c.title = t; c.kind = 1; c.i = f; c.offset = offset; cols.push_back(c); return *this; }
    RecordTable& text(const char* t, std::function<std::string&(T&)> f) { Col c; c.title = t; c.kind = 2; c.s = f; cols.push_back(c); return *this; }
    RecordTable& choice(const char* t, std::function<int&(T&)> f, std::vector<std::string> o, int base = 0) {
        Col c; c.title = t; c.kind = 3; c.opts = std::move(o); c.optBase = base;
        c.gi = [f, base](T& x) { return f(x) - base; }; c.si = [f, base](T& x, int k) { f(x) = k + base; };
        cols.push_back(c); return *this;
    }
    // list column over any representation (enums, chars, flag pairs): get returns the option index, set receives it
    RecordTable& choiceFn(const char* t, std::function<int(T&)> get, std::function<void(T&, int)> set, std::vector<std::string> o) {
        Col c; c.title = t; c.kind = 3; c.opts = std::move(o); c.gi = get; c.si = set; cols.push_back(c); return *this;
    }
    RecordTable& flag(const char* t, std::function<bool&(T&)> f) { Col c; c.title = t; c.kind = 4; c.b = f; c.opts = {"no", "yes"}; cols.push_back(c); return *this; }
    RecordTable& computed(const char* t, std::function<double(const T&)> f, double scale = 1) { Col c; c.title = t; c.kind = 5; c.calc = f; c.scale = scale; cols.push_back(c); return *this; }
    // optional number (defaulted / "from elsewhere" state): blank, '*', '1*' or the unset text restores the default
    RecordTable& optNum(const char* t, std::function<double&(T&)> f, std::function<bool(T&)> isUnset, std::function<void(T&)> unset,
                        const char* unsetText = "default", double scale = 1) {
        Col c; c.title = t; c.kind = 6; c.d = f; c.scale = scale; c.isUnset = isUnset; c.unset = unset; c.unsetText = unsetText; cols.push_back(c); return *this;
    }
    // list column over a string field (well / component names)
    RecordTable& choiceStr(const char* t, std::function<std::string&(T&)> f, std::vector<std::string> o) {
        Col c; c.title = t; c.kind = 3; c.opts = o;
        c.gi = [f, o](T& x) { for (size_t k = 0; k < o.size(); ++k) if (o[k] == f(x)) return int(k); return -1; };
        c.si = [f, o](T& x, int k) { if (k >= 0 && k < int(o.size())) f(x) = o[size_t(k)]; };
        c.s = f; cols.push_back(c); return *this;
    }
    RecordTable& fixedRows() { canInsert = canErase = false; return *this; }

    static std::string lower(std::string x) { for (char& ch : x) ch = char(std::tolower(static_cast<unsigned char>(ch))); return x; }
    static int matchOpt(const Col& c, const std::string& t) {
        const std::string lt = lower(t);
        for (size_t k = 0; k < c.opts.size(); ++k) if (lower(c.opts[k]) == lt) return int(k);
        for (size_t k = 0; k < c.opts.size(); ++k) if (!lt.empty() && lower(c.opts[k]).rfind(lt, 0) == 0) return int(k);
        char* e = nullptr; const long v = std::strtol(t.c_str(), &e, 10);
        if (e && e != t.c_str() && *e == 0) { const long k = v - c.optBase; if (k >= 0 && k < long(c.opts.size())) return int(k); }
        if (c.kind == 4) { if (lt == "1" || lt == "true" || lt == "on" || lt == "y") return 1; if (lt == "0" || lt == "false" || lt == "off" || lt == "n") return 0; }
        return -1;
    }
    bool draw(const char* id, std::vector<T>& rows, float height = 0, int* selected = nullptr, const char* fmt = "%.6g", bool deleteButtons = true) {
        GridAccess g;
        const char* fmtNum = fmt;
        g.count = [&] { return int(rows.size()); };
        g.cols = int(cols.size());
        g.numGet = [&](int r, int c, double& v) -> bool {
            const Col& C = cols[size_t(c)]; T& x = rows[size_t(r)];
            switch (C.kind) {
                case 0: v = C.d(x) * C.scale; return true;
                case 1: v = double(C.i(x)) + C.offset; return true;
                case 5: v = C.calc(x) * C.scale; return true;
                default: return false;
            }
        };
        g.numSet = [&](int r, int c, double v) {
            const Col& C = cols[size_t(c)]; T& x = rows[size_t(r)];
            if (C.kind == 0) C.d(x) = v / C.scale;
            else if (C.kind == 1) C.i(x) = int(std::lround(v - C.offset));
            else return;
            if (onChange) onChange(x);
        };
        g.getText = [&](int r, int c) -> std::string {
            const Col& C = cols[size_t(c)]; T& x = rows[size_t(r)];
            if (C.kind == 2) return C.s(x);
            if (C.kind == 3) { const int k = C.gi(x); return k >= 0 && k < int(C.opts.size()) ? C.opts[size_t(k)] : (C.s ? C.s(x) : std::to_string(k + C.optBase)); }
            if (C.kind == 6) { if (C.isUnset(x)) return C.unsetText; char b[48]; std::snprintf(b, sizeof b, fmtNum, C.d(x) * C.scale); return b; }
            if (C.kind == 4) return C.b(x) ? "yes" : "no";
            return std::string();
        };
        g.setText = [&](int r, int c, const std::string& t) {
            const Col& C = cols[size_t(c)]; T& x = rows[size_t(r)];
            if (C.kind == 2) C.s(x) = t;
            else if (C.kind == 3) { const int k = matchOpt(C, t); if (k >= 0) C.si(x, k); }
            else if (C.kind == 4) { const int k = matchOpt(C, t); if (k >= 0) C.b(x) = k == 1; }
            else if (C.kind == 6) {
                if (t.empty() || t == "*" || t == "1*" || lower(t) == lower(C.unsetText)) C.unset(x);
                else { char* e = nullptr; const double v = std::strtod(t.c_str(), &e); if (e && e != t.c_str() && *e == 0 && std::isfinite(v)) C.d(x) = v / C.scale; }
            }
            if (onChange) onChange(x);
        };
        g.choices = [&](int c) -> const std::vector<std::string>* { const Col& C = cols[size_t(c)]; return (C.kind == 3 || C.kind == 4) ? &C.opts : nullptr; };
        g.readOnly = [&](int, int c) { return cols[size_t(c)].kind == 5; };
        if (canInsert)
            g.insert = [&](int at, int n) {
                at = std::max(0, std::min(at, int(rows.size())));
                const T v = newRow ? newRow(at) : (rows.empty() ? T{} : rows[size_t(std::min(at, int(rows.size()) - 1))]);
                rows.insert(rows.begin() + at, size_t(std::max(0, n)), v);
            };
        if (canErase)
            g.erase = [&](int r0, int r1) {
                r0 = std::max(0, r0); r1 = std::min(r1, int(rows.size()) - 1);
                if (r0 <= r1) rows.erase(rows.begin() + r0, rows.begin() + r1 + 1);
            };
        std::vector<const char*> names;
        for (const Col& C : cols) names.push_back(C.title.c_str());
        return gridTable(id, g, names, height, selected, fmt, deleteButtons && canErase);
    }
};

// parse clipboard text into rows of cells (also used by tests / other importers)
std::vector<std::vector<std::string>> parseTableText(const std::string& text);

// ---- editable numeric table over a vector of rows
// cols: column titles; get(row, col) returns pointer to the double for that cell.
template <class T>
bool editTable(const char* id, std::vector<T>& rows, const std::vector<const char*>& cols,
               const std::function<double*(T&, int)>& get, float height = 0, int* selected = nullptr,
               const char* fmt = "%.6g") {
    GridAccess g;
    g.count = [&] { return int(rows.size()); };
    g.cols = int(cols.size());
    g.num = [&](int r, int c) { return get(rows[size_t(r)], c); };
    g.insert = [&](int at, int n) {
        at = std::max(0, std::min(at, int(rows.size())));
        const T v = rows.empty() ? T{} : rows[size_t(std::min(at, int(rows.size()) - 1))];
        rows.insert(rows.begin() + at, size_t(n), v);
    };
    g.erase = [&](int r0, int r1) {
        r0 = std::max(0, r0); r1 = std::min(r1, int(rows.size()) - 1);
        if (r1 >= r0) rows.erase(rows.begin() + r0, rows.begin() + r1 + 1);
    };
    return gridTable(id, g, cols, height, selected, fmt);
}
} // namespace ui
