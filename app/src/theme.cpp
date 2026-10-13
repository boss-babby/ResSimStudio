#include "theme.h"
#include "implot.h"
#include <cstdio>
#include "icons.h"
#include "resources.h"

Fonts g_fonts;
static float g_dpi = 1.0f, g_uiScale = 1.0f;

static bool fileExists(const char* p) { FILE* f = std::fopen(p, "rb"); if (f) { std::fclose(f); return true; } return false; }

static void styleMetrics(ImGuiStyle& s) {
    s.WindowPadding = ImVec2(10, 10);
    s.FramePadding = ImVec2(8, 5);
    s.ItemSpacing = ImVec2(8, 6);
    s.ItemInnerSpacing = ImVec2(6, 4);
    s.CellPadding = ImVec2(6, 3);
    s.IndentSpacing = 18;
    s.ScrollbarSize = 13;
    s.GrabMinSize = 10;
    s.WindowBorderSize = 1; s.ChildBorderSize = 1; s.PopupBorderSize = 1; s.FrameBorderSize = 0; s.TabBorderSize = 0;
    s.WindowRounding = 8; s.ChildRounding = 6; s.FrameRounding = 5; s.PopupRounding = 6; s.ScrollbarRounding = 8;
    s.GrabRounding = 4; s.TabRounding = 5;
    s.WindowTitleAlign = ImVec2(0.0f, 0.5f);
    s.SeparatorTextBorderSize = 1;
    s.DockingSeparatorSize = 3;
}

void applyTheme(float dpi) {
    ImGuiIO& io = ImGui::GetIO();
    // ---- fonts (Segoe UI family ships with every Windows 10/11)
    const char* ui = "C:\\Windows\\Fonts\\segoeui.ttf";
    const char* uib = "C:\\Windows\\Fonts\\seguisb.ttf";
    const char* sym = "C:\\Windows\\Fonts\\seguisym.ttf";
    const char* mono = "C:\\Windows\\Fonts\\consola.ttf";
    const float base = 17.0f;
    // The Lucide icons live in the Private Use Area (U+E038-U+E612). Segoe UI Symbol (and some Segoe UI builds) carry their
    // own PUA glyphs there; ImGui takes a glyph from the FIRST source that has it, so every non-icon font must exclude the
    // icon range - otherwise Windows shows random Segoe symbols instead of the icons (Wine has no Segoe, which hid this).
    static const ImWchar iconRange[] = {ICON_MIN, ICON_MAX, 0};
    static const ImWchar puaExclude[] = {0xE000, 0xF8FF, 0};
    ImFontConfig cfg; cfg.OversampleH = 2; cfg.OversampleV = 2; cfg.GlyphExcludeRanges = puaExclude;
    auto withSymbols = [&](float size) {
        if (!fileExists(sym)) return;
        ImFontConfig m; m.MergeMode = true; m.GlyphMinAdvanceX = size * 0.9f; m.GlyphExcludeRanges = puaExclude;
        io.Fonts->AddFontFromFileTTF(sym, size, &m);
    };
    // embedded Lucide icon font, merged into every UI font so ICON_* glyphs can be mixed with text
    const std::string_view iconTtf = appResource(RES_ICON_FONT);
    auto withIcons = [&](float size) {
        if (iconTtf.empty()) return;
        ImFontConfig m; m.MergeMode = true; m.FontDataOwnedByAtlas = false; m.PixelSnapH = true;
        m.GlyphMinAdvanceX = size * 1.05f; m.GlyphOffset = ImVec2(0.0f, size * 0.12f);
        m.OversampleH = 2; m.OversampleV = 2;
        io.Fonts->AddFontFromMemoryTTF(const_cast<char*>(iconTtf.data()), int(iconTtf.size()), size * 0.95f, &m, iconRange);
    };
    if (fileExists(ui)) { g_fonts.ui = io.Fonts->AddFontFromFileTTF(ui, base, &cfg); withIcons(base); withSymbols(base); }
    else {
        // no system UI font (e.g. Wine): built-in font, still with the symbol / icon glyphs merged in
        ImFontConfig d; d.SizePixels = base; d.GlyphExcludeRanges = puaExclude;
        g_fonts.ui = io.Fonts->AddFontDefault(&d);
        withIcons(base);
    }
    g_fonts.bold = fileExists(uib) ? io.Fonts->AddFontFromFileTTF(uib, base, &cfg) : g_fonts.ui;
    if (g_fonts.bold != g_fonts.ui) { withIcons(base); withSymbols(base); }
    g_fonts.title = fileExists(uib) ? io.Fonts->AddFontFromFileTTF(uib, base * 1.45f, &cfg) : g_fonts.ui;
    if (g_fonts.title != g_fonts.ui) withIcons(base * 1.45f);
    g_fonts.mono = fileExists(mono) ? io.Fonts->AddFontFromFileTTF(mono, base - 2.0f, &cfg) : g_fonts.ui;
    if (g_fonts.mono != g_fonts.ui) withIcons(base - 2.0f);
    io.FontDefault = g_fonts.ui;

    // ---- style
    ImGuiStyle& s = ImGui::GetStyle();
    styleMetrics(s);

    ImVec4* c = s.Colors;
    const ImVec4 bg0(0.075f, 0.078f, 0.090f, 1.00f);   // app background
    const ImVec4 bg1(0.105f, 0.110f, 0.125f, 1.00f);   // panels
    const ImVec4 bg2(0.145f, 0.152f, 0.172f, 1.00f);   // frames
    const ImVec4 bg3(0.190f, 0.200f, 0.228f, 1.00f);   // hovered frames
    const ImVec4 acc(0.251f, 0.518f, 0.965f, 1.00f);   // accent blue
    const ImVec4 accH(0.345f, 0.600f, 1.000f, 1.00f);
    const ImVec4 txt(0.905f, 0.915f, 0.935f, 1.00f);
    c[ImGuiCol_Text] = txt;
    c[ImGuiCol_TextDisabled] = ImVec4(0.50f, 0.53f, 0.58f, 1.00f);
    c[ImGuiCol_WindowBg] = bg1;
    c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_PopupBg] = ImVec4(0.12f, 0.125f, 0.142f, 0.98f);
    c[ImGuiCol_Border] = ImVec4(0.22f, 0.23f, 0.26f, 0.70f);
    c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_FrameBg] = bg2;
    c[ImGuiCol_FrameBgHovered] = bg3;
    c[ImGuiCol_FrameBgActive] = ImVec4(0.22f, 0.24f, 0.29f, 1.00f);
    c[ImGuiCol_TitleBg] = bg0;
    c[ImGuiCol_TitleBgActive] = bg0;
    c[ImGuiCol_TitleBgCollapsed] = bg0;
    c[ImGuiCol_MenuBarBg] = bg0;
    c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ScrollbarGrab] = ImVec4(0.26f, 0.28f, 0.32f, 1.00f);
    c[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.33f, 0.35f, 0.40f, 1.00f);
    c[ImGuiCol_ScrollbarGrabActive] = acc;
    c[ImGuiCol_CheckMark] = accH;
    c[ImGuiCol_SliderGrab] = acc;
    c[ImGuiCol_SliderGrabActive] = accH;
    c[ImGuiCol_Button] = ImVec4(0.17f, 0.18f, 0.205f, 1.00f);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.22f, 0.24f, 0.28f, 1.00f);
    c[ImGuiCol_ButtonActive] = acc;
    c[ImGuiCol_Header] = ImVec4(0.251f, 0.518f, 0.965f, 0.22f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.251f, 0.518f, 0.965f, 0.32f);
    c[ImGuiCol_HeaderActive] = ImVec4(0.251f, 0.518f, 0.965f, 0.45f);
    c[ImGuiCol_Separator] = ImVec4(0.22f, 0.23f, 0.26f, 1.00f);
    c[ImGuiCol_SeparatorHovered] = acc;
    c[ImGuiCol_SeparatorActive] = accH;
    c[ImGuiCol_ResizeGrip] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ResizeGripHovered] = acc;
    c[ImGuiCol_ResizeGripActive] = accH;
    c[ImGuiCol_Tab] = bg0;
    c[ImGuiCol_TabHovered] = bg3;
    c[ImGuiCol_TabSelected] = bg1;
    c[ImGuiCol_TabSelectedOverline] = acc;
    c[ImGuiCol_TabDimmed] = bg0;
    c[ImGuiCol_TabDimmedSelected] = bg1;
    c[ImGuiCol_TabDimmedSelectedOverline] = ImVec4(0.35f, 0.38f, 0.45f, 1.0f);
    c[ImGuiCol_DockingPreview] = ImVec4(0.251f, 0.518f, 0.965f, 0.55f);
    c[ImGuiCol_DockingEmptyBg] = bg0;
    c[ImGuiCol_PlotLines] = acc;
    c[ImGuiCol_PlotHistogram] = acc;
    c[ImGuiCol_TableHeaderBg] = ImVec4(0.13f, 0.137f, 0.155f, 1.00f);
    c[ImGuiCol_TableBorderStrong] = ImVec4(0.22f, 0.23f, 0.26f, 1.00f);
    c[ImGuiCol_TableBorderLight] = ImVec4(0.17f, 0.18f, 0.20f, 1.00f);
    c[ImGuiCol_TableRowBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_TableRowBgAlt] = ImVec4(1, 1, 1, 0.025f);
    c[ImGuiCol_TextSelectedBg] = ImVec4(0.251f, 0.518f, 0.965f, 0.35f);
    c[ImGuiCol_NavCursor] = acc;
    c[ImGuiCol_ModalWindowDimBg] = ImVec4(0, 0, 0, 0.55f);

    s.ScaleAllSizes(dpi * g_uiScale);
    s.FontScaleDpi = dpi;
    s.FontScaleMain = g_uiScale;
    g_dpi = dpi;

    // ---- ImPlot style
    ImPlotStyle& ps = ImPlot::GetStyle();
    ImPlot::StyleColorsDark();
    ps.PlotBorderSize = 0;
    ps.PlotPadding = ImVec2(10, 10);
    ps.LegendPadding = ImVec2(8, 6);
    ps.MinorAlpha = 0.12f;
    ps.Colors[ImPlotCol_PlotBg] = ImVec4(0.085f, 0.09f, 0.10f, 1.0f);
    ps.Colors[ImPlotCol_FrameBg] = ImVec4(0, 0, 0, 0);
    ps.Colors[ImPlotCol_LegendBg] = ImVec4(0.11f, 0.115f, 0.13f, 0.92f);
    ps.Colors[ImPlotCol_AxisGrid] = ImVec4(1, 1, 1, 0.08f);
}

void setUiScale(float uiScale) {
    uiScale = uiScale < 0.6f ? 0.6f : uiScale > 1.6f ? 1.6f : uiScale;
    g_uiScale = uiScale;
    if (!ImGui::GetCurrentContext()) return;
    ImGuiStyle& s = ImGui::GetStyle();
    styleMetrics(s);
    s.ScaleAllSizes(g_dpi * uiScale);
    s.FontScaleDpi = g_dpi;
    s.FontScaleMain = uiScale;
}
float uiScale() { return g_uiScale; }
