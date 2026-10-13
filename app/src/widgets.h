// Small UI helpers: property grids, editable tables, std::string input.
#pragma once
#include "imgui.h"
#include "theme.h"
#include <string>
#include <vector>
#include <functional>
#include <cfloat>
#include <algorithm>

namespace ui {

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

inline void header(const char* text) {
    ImGui::Spacing();
    if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
    ImGui::SeparatorText(text);
    if (g_fonts.bold) ImGui::PopFont();
}

inline void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

inline void hint(const char* text) {
    ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(text);
    ImGui::PopTextWrapPos();
    ImGui::PopStyleColor();
}

// ---- two-column property grid
inline bool beginProps(const char* id, float labelFrac = 0.45f) {
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingStretchProp)) return false;
    ImGui::TableSetupColumn("k", ImGuiTableColumnFlags_WidthStretch, labelFrac);
    ImGui::TableSetupColumn("v", ImGuiTableColumnFlags_WidthStretch, 1.0f - labelFrac);
    return true;
}
inline void endProps() { ImGui::EndTable(); }
inline void label(const char* text, const char* tip = nullptr) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.72f, 0.75f, 0.80f, 1));
    ImGui::TextUnformatted(text);
    ImGui::PopStyleColor();
    if (tip && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", tip);
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
}
inline bool propD(const char* text, double* v, const char* fmt = "%.6g", const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::InputDouble("##v", v, 0, 0, fmt);
    ImGui::PopID();
    return c;
}
inline bool propI(const char* text, int* v, int lo, int hi, const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    bool c = ImGui::InputInt("##v", v, 1, 10);
    if (*v < lo) { *v = lo; } if (*v > hi) { *v = hi; }
    ImGui::PopID();
    return c;
}
inline bool propF(const char* text, float* v, float lo, float hi, const char* fmt = "%.3g", const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::SliderFloat("##v", v, lo, hi, fmt);
    ImGui::PopID();
    return c;
}
inline bool propSliderD(const char* text, double* v, double lo, double hi, const char* fmt = "%.3g", const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    float f = float(*v);
    const bool c = ImGui::SliderFloat("##v", &f, float(lo), float(hi), fmt);
    if (c) *v = f;
    ImGui::PopID();
    return c;
}
inline bool propBool(const char* text, bool* v, const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::Checkbox("##v", v);
    ImGui::PopID();
    return c;
}
inline bool propCombo(const char* text, int* v, const char* const items[], int n, const char* tip = nullptr) {
    label(text, tip);
    ImGui::PushID(text);
    const bool c = ImGui::Combo("##v", v, items, n);
    ImGui::PopID();
    return c;
}
inline void propText(const char* text, const char* value) {
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

// ---- accent button
inline bool accentButton(const char* text, const ImVec2& size = ImVec2(0, 0)) {
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.20f, 0.45f, 0.90f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.28f, 0.53f, 0.98f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.16f, 0.38f, 0.80f, 1));
    const bool r = ImGui::Button(text, size);
    ImGui::PopStyleColor(3);
    return r;
}
inline bool dangerButton(const char* text, const ImVec2& size = ImVec2(0, 0)) {
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.55f, 0.16f, 0.18f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.70f, 0.22f, 0.24f, 1));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.45f, 0.12f, 0.14f, 1));
    const bool r = ImGui::Button(text, size);
    ImGui::PopStyleColor(3);
    return r;
}

// ---- editable numeric table over a vector of rows
// cols: column titles; get(row, col) returns pointer to the double for that cell.
template <class T>
bool editTable(const char* id, std::vector<T>& rows, const std::vector<const char*>& cols,
               const std::function<double*(T&, int)>& get, float height = 0, int* selected = nullptr,
               const char* fmt = "%.6g") {
    bool changed = false;
    const int nc = int(cols.size());
    ImGui::PushID(id);
    const ImGuiTableFlags fl = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY |
                               ImGuiTableFlags_SizingStretchSame;
    const float h = height > 0 ? height : std::min(320.0f, ImGui::GetFrameHeightWithSpacing() * (float(rows.size()) + 1.6f));
    int del = -1;
    if (ImGui::BeginTable("t", nc + 2, fl, ImVec2(0, h))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.8f);
        for (int c = 0; c < nc; ++c) ImGui::TableSetupColumn(cols[c]);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
        ImGui::TableHeadersRow();
        ImGuiListClipper clip;
        clip.Begin(int(rows.size()));
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                ImGui::PushID(r);
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (selected) {
                    char b[16]; std::snprintf(b, sizeof(b), "%d", r + 1);
                    if (ImGui::Selectable(b, *selected == r, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) *selected = r;
                } else ImGui::Text("%d", r + 1);
                for (int c = 0; c < nc; ++c) {
                    ImGui::TableSetColumnIndex(c + 1);
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::PushID(c);
                    if (ImGui::InputDouble("##c", get(rows[r], c), 0, 0, fmt)) changed = true;
                    ImGui::PopID();
                }
                ImGui::TableSetColumnIndex(nc + 1);
                if (ImGui::SmallButton("\xC3\x97")) del = r;     // multiplication sign
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Delete row");
                ImGui::PopID();
            }
        ImGui::EndTable();
    }
    if (del >= 0) { rows.erase(rows.begin() + del); changed = true; if (selected && *selected >= int(rows.size())) *selected = int(rows.size()) - 1; }
    if (ImGui::SmallButton("+ Add row")) {
        if (!rows.empty()) rows.push_back(rows.back()); else rows.push_back(T{});
        changed = true;
    }
    ImGui::PopID();
    return changed;
}

} // namespace ui
