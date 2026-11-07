#include "theme.h"
#include "implot.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "icons.h"
#include "resources.h"

Fonts g_fonts;
static float g_dpi = 1.0f, g_uiScale = 1.0f;
static int g_theme = 0, g_accent = 0;

static bool fileExists(const char* p) { FILE* f = std::fopen(p, "rb"); if (f) { std::fclose(f); return true; } return false; }

static void styleMetrics(ImGuiStyle& s) {
    // Start from unscaled defaults every time: ScaleAllSizes() multiplies (and truncates) every metric, so calling it
    // again on an already scaled style compounded the scale of the metrics not set below - the caret thickness
    // InputTextCursorSize became trunc(1 * 0.9) = 0 px at UI scales below 100 %: no visible text cursor anywhere
    // (input fields and the deck editor, which uses the same thickness).
    {
        ImVec4 cols[ImGuiCol_COUNT];
        for (int i = 0; i < ImGuiCol_COUNT; ++i) cols[i] = s.Colors[i];
        const float fsb = s.FontSizeBase;
        s = ImGuiStyle();
        for (int i = 0; i < ImGuiCol_COUNT; ++i) s.Colors[i] = cols[i];
        s.FontSizeBase = fsb;
    }
    s.WindowPadding = ImVec2(10, 10);
    s.FramePadding = ImVec2(8, 5);
    s.ItemSpacing = ImVec2(8, 6);
    s.ItemInnerSpacing = ImVec2(6, 4);
    s.CellPadding = ImVec2(6, 3);
    s.IndentSpacing = 18;
    s.ScrollbarSize = 13;
    s.GrabMinSize = 10;
    s.WindowBorderSize = 1; s.ChildBorderSize = 1; s.PopupBorderSize = 1; s.FrameBorderSize = g_theme == 2 ? 1.0f : 0.0f; s.TabBorderSize = 0;
    // C41 (G5) design tokens: one corner radius (6) for frames / children / popups / tabs, windows one step larger
    s.WindowRounding = 8; s.ChildRounding = 6; s.FrameRounding = 6; s.PopupRounding = 6; s.ScrollbarRounding = 8;
    s.GrabRounding = 4; s.TabRounding = 6;
    s.TabCloseButtonMinWidthSelected = -1.0f; s.TabCloseButtonMinWidthUnselected = -1.0f;   // close crosses always visible
    s.WindowTitleAlign = ImVec2(0.0f, 0.5f);
    s.SeparatorTextBorderSize = 1;
    s.DockingSeparatorSize = 3;
}

void applyTheme(float dpi) {
    ImGuiIO& io = ImGui::GetIO();
    // ---- fonts: Segoe UI family on Windows 10/11; elsewhere (Linux, macOS, Wine) the embedded Inter + JetBrains Mono, with a
    // system font merged in for symbols outside the latin subset (arrows, maths, Greek)
#ifdef _WIN32
    const char* ui = "C:\\Windows\\Fonts\\segoeui.ttf";
    const char* uib = "C:\\Windows\\Fonts\\seguisb.ttf";
    const char* mono = "C:\\Windows\\Fonts\\consola.ttf";
    static const char* symCands[] = {"C:\\Windows\\Fonts\\seguisym.ttf"};
#else
    const char* ui = "";
    const char* uib = "";
    const char* mono = "";
    static const char* symCands[] = {"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/TTF/DejaVuSans.ttf",
                                     "/usr/share/fonts/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
                                     "/usr/share/fonts/noto/NotoSans-Regular.ttf", "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
                                     "/Library/Fonts/Arial Unicode.ttf"};
#endif
    const char* sym = "";
    for (const char* c : symCands) if (fileExists(c)) { sym = c; break; }
    const bool segoe = *ui && fileExists(ui);
    const std::string_view interR = appResource(RES_FONT_UI), interB = appResource(RES_FONT_UI_BOLD), jbm = appResource(RES_FONT_MONO);
    const float base = segoe ? 17.0f : 16.0f;          // Inter has a larger x-height than Segoe UI
    // The Lucide icons live in the Private Use Area (U+E038-U+E612). Segoe UI Symbol (and some Segoe UI builds) carry their
    // own PUA glyphs there; ImGui takes a glyph from the FIRST source that has it, so every non-icon font must exclude the
    // icon range - otherwise Windows shows random Segoe symbols instead of the icons (Wine has no Segoe, which hid this).
    static const ImWchar iconRange[] = {ICON_MIN, ICON_MAX, 0};
    static const ImWchar puaExclude[] = {0xE000, 0xF8FF, 0};
    ImFontConfig cfg; cfg.OversampleH = 2; cfg.OversampleV = 2; cfg.GlyphExcludeRanges = puaExclude;
    auto withSymbols = [&](float size) {
        if (!*sym) return;
        ImFontConfig m; m.MergeMode = true; m.GlyphMinAdvanceX = segoe ? size * 0.9f : 0.0f; m.GlyphExcludeRanges = puaExclude;
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
    auto embedded = [&](std::string_view ttf, float size) -> ImFont* {
        if (ttf.empty()) return nullptr;
        ImFontConfig m = cfg; m.FontDataOwnedByAtlas = false;
        return io.Fonts->AddFontFromMemoryTTF(const_cast<char*>(ttf.data()), int(ttf.size()), size, &m);
    };
    auto uiFont = [&](bool bold, float size) -> ImFont* {
        ImFont* f = nullptr;
        if (segoe) f = io.Fonts->AddFontFromFileTTF(bold && fileExists(uib) ? uib : ui, size, &cfg);
        if (!f) f = embedded(bold ? interB : interR, size);
        if (!f) { ImFontConfig d; d.SizePixels = size; d.GlyphExcludeRanges = puaExclude; f = io.Fonts->AddFontDefault(&d); }
        withIcons(size); withSymbols(size);
        return f;
    };
    g_fonts.ui = uiFont(false, base);
    g_fonts.bold = uiFont(true, base);
    g_fonts.title = uiFont(true, base * 1.45f);
    {
        ImFont* f = (*mono && fileExists(mono)) ? io.Fonts->AddFontFromFileTTF(mono, base - 2.0f, &cfg) : embedded(jbm, base - 2.5f);
        if (f) { withIcons(base - 2.0f); withSymbols(base - 2.0f); }
        g_fonts.mono = f ? f : g_fonts.ui;
    }
    io.FontDefault = g_fonts.ui;

    // ---- style
    ImGuiStyle& s = ImGui::GetStyle();
    styleMetrics(s);

    s.ScaleAllSizes(dpi * g_uiScale);
    s.InputTextCursorSize = std::max(1.0f, std::round(1.5f * (dpi * g_uiScale)));   // caret: never 0 px, 2 px at 100 %+
    // ImGui draws 1-px details (tab close crosses, link underlines, colour-button borders) with (int)_MainScale:
    // below 100 % that is 0 px and they vanish
    s._MainScale = dpi * g_uiScale < 1.0f ? 1.0f : dpi * g_uiScale;   // (ScaleAllSizes multiplies it; set it absolutely)
    s.FontScaleDpi = dpi;
    s.FontScaleMain = g_uiScale;
    g_dpi = dpi;

    applyColors(g_theme, g_accent);   // ImPlot colours (ImPlot context exists now)
}

void setUiScale(float uiScale) {
    uiScale = uiScale < 0.6f ? 0.6f : uiScale > 1.6f ? 1.6f : uiScale;
    g_uiScale = uiScale;
    if (!ImGui::GetCurrentContext()) return;
    ImGuiStyle& s = ImGui::GetStyle();
    styleMetrics(s);
    s.ScaleAllSizes(g_dpi * uiScale);
    s.InputTextCursorSize = std::max(1.0f, std::round(1.5f * (g_dpi * uiScale)));   // caret: never 0 px, 2 px at 100 %+
    s._MainScale = g_dpi * uiScale < 1.0f ? 1.0f : g_dpi * uiScale;   // see applyTheme
    s.FontScaleDpi = g_dpi;
    s.FontScaleMain = uiScale;
}
float uiScale() { return g_uiScale; }
namespace tok { float sp(int step) { return 4.0f * float(step < 0 ? 0 : step) * g_uiScale; } }

// ============================================================================ C39 colour themes
namespace {
struct Rgb { float r, g, b; };
const Rgb kAccents[] = {{0.251f, 0.518f, 0.965f}, {0.078f, 0.722f, 0.651f}, {0.545f, 0.361f, 0.965f}, {0.976f, 0.451f, 0.086f}, {0.133f, 0.773f, 0.369f}};
const char* kThemeNames[] = {"Dark", "Light", "High contrast"};
const char* kAccentNames[] = {"Blue", "Teal", "Violet", "Orange", "Green"};
ImVec4 mixc(const ImVec4& a, const ImVec4& b, float t) { return ImVec4(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t); }
ImVec4 alpha(const ImVec4& c, float a) { return ImVec4(c.x, c.y, c.z, a); }
}  // namespace

const char* const* themeNames() { return kThemeNames; }
const char* const* accentNames() { return kAccentNames; }

static int g_bg3d = 0;
void set3dBackground(int mode) { g_bg3d = mode < 0 || mode > 4 ? 0 : mode; }
int background3dMode() { return g_bg3d; }
const float* background3d() {
    static const float kDark[3] = {0.085f, 0.09f, 0.105f}, kLight[3] = {0.90f, 0.915f, 0.935f}, kWhite[3] = {1.0f, 1.0f, 1.0f},
                       kBlack[3] = {0.0f, 0.0f, 0.0f};
    switch (g_bg3d) {
        case 1: return themeId() == 1 ? kLight : kDark;
        case 2: return kLight;
        case 3: return kWhite;
        case 4: return kBlack;
        default: return kDark;
    }
}
const char* const* background3dNames() { static const char* n[] = {"Dark", "Follow the theme", "Light grey", "White", "Black"}; return n; }
ImVec4 accentColor(int a) { const Rgb& c = kAccents[a < 0 || a > 4 ? 0 : a]; return ImVec4(c.r, c.g, c.b, 1.0f); }
int themeId() { return g_theme; }
int accentId() { return g_accent; }
bool themeIsLight() { return g_theme == 1; }
// C41 (G5): elevation by surface tint - level 0 = window background, each level mixes 4 % (dark) / 3.5 % (light)
// toward the text colour; high contrast uses stronger steps
ImVec4 themeSurface(int level) {
    const ImGuiStyle& st = ImGui::GetStyle();
    const ImVec4 b = st.Colors[ImGuiCol_WindowBg], t = st.Colors[ImGuiCol_Text];
    const float k = (g_theme == 1 ? 0.035f : g_theme == 2 ? 0.08f : 0.04f) * float(level < 0 ? 0 : level > 4 ? 4 : level);
    return ImVec4(b.x + (t.x - b.x) * k, b.y + (t.y - b.y) * k, b.z + (t.z - b.z) * k, 1.0f);
}
ImVec4 themeTextLevel(int level) {
    if (level <= 0) return ImGui::GetStyle().Colors[ImGuiCol_Text];
    if (level == 1) return themeLabelText();
    return col4(Pal::Muted);
}
ImVec4 themeStatusBg() { return g_theme == 1 ? ImVec4(0.90f, 0.91f, 0.93f, 1) : g_theme == 2 ? ImVec4(0, 0, 0, 1) : ImVec4(0.055f, 0.06f, 0.07f, 1); }
ImVec4 themeCardBg() { return g_theme == 1 ? ImVec4(1, 1, 1, 0.98f) : g_theme == 2 ? ImVec4(0.05f, 0.05f, 0.06f, 0.98f) : ImVec4(0.13f, 0.14f, 0.17f, 0.97f); }

void applyColors(int theme, int accent) {
    g_theme = theme < 0 || theme > 2 ? 0 : theme;
    g_accent = accent < 0 || accent > 4 ? 0 : accent;
    const Rgb A = kAccents[g_accent];
    const ImVec4 acc(A.r, A.g, A.b, 1.0f);
    const bool light = g_theme == 1, hc = g_theme == 2;
    const ImVec4 accH = light ? mixc(acc, ImVec4(0, 0, 0, 1), 0.15f) : mixc(acc, ImVec4(1, 1, 1, 1), 0.18f);
    Pal::Accent = ImGui::ColorConvertFloat4ToU32(acc);
    Pal::Muted = light ? IM_COL32(96, 102, 116, 255) : hc ? IM_COL32(190, 196, 208, 255) : IM_COL32(140, 146, 160, 255);
    Pal::Ok = light ? IM_COL32(16, 132, 82, 255) : IM_COL32(72, 199, 142, 255);
    Pal::Warn = light ? IM_COL32(176, 100, 0, 255) : IM_COL32(255, 179, 64, 255);
    Pal::Err = light ? IM_COL32(200, 36, 36, 255) : IM_COL32(255, 92, 92, 255);
    // surfaces: bg0 app background / title bars, bg1 panels, bg2 frames, bg3 hovered frames
    ImVec4 bg0, bg1, bg2, bg3, txt, txtDis, border, btn, btnH, popup, tableHdr, tableLight, scroll, scrollH;
    if (light) {
        bg0 = ImVec4(0.905f, 0.912f, 0.928f, 1); bg1 = ImVec4(0.972f, 0.975f, 0.982f, 1); bg2 = ImVec4(0.900f, 0.908f, 0.925f, 1);
        bg3 = ImVec4(0.855f, 0.866f, 0.890f, 1); txt = ImVec4(0.105f, 0.115f, 0.140f, 1); txtDis = ImVec4(0.43f, 0.45f, 0.50f, 1);
        border = ImVec4(0.78f, 0.79f, 0.82f, 1); btn = ImVec4(0.885f, 0.893f, 0.912f, 1); btnH = ImVec4(0.835f, 0.848f, 0.875f, 1);
        popup = ImVec4(0.99f, 0.99f, 1.0f, 0.99f); tableHdr = ImVec4(0.915f, 0.922f, 0.938f, 1); tableLight = ImVec4(0.86f, 0.87f, 0.89f, 1);
        scroll = ImVec4(0.75f, 0.77f, 0.80f, 1); scrollH = ImVec4(0.65f, 0.67f, 0.71f, 1);
    } else if (hc) {
        bg0 = ImVec4(0, 0, 0, 1); bg1 = ImVec4(0.03f, 0.03f, 0.035f, 1); bg2 = ImVec4(0.11f, 0.11f, 0.13f, 1); bg3 = ImVec4(0.20f, 0.20f, 0.24f, 1);
        txt = ImVec4(1, 1, 1, 1); txtDis = ImVec4(0.72f, 0.74f, 0.78f, 1); border = ImVec4(0.62f, 0.64f, 0.70f, 1);
        btn = ImVec4(0.13f, 0.13f, 0.16f, 1); btnH = ImVec4(0.24f, 0.25f, 0.30f, 1); popup = ImVec4(0.02f, 0.02f, 0.03f, 1);
        tableHdr = ImVec4(0.10f, 0.10f, 0.12f, 1); tableLight = ImVec4(0.35f, 0.36f, 0.40f, 1); scroll = ImVec4(0.45f, 0.47f, 0.52f, 1);
        scrollH = ImVec4(0.60f, 0.62f, 0.68f, 1);
    } else {
        bg0 = ImVec4(0.075f, 0.078f, 0.090f, 1); bg1 = ImVec4(0.105f, 0.110f, 0.125f, 1); bg2 = ImVec4(0.145f, 0.152f, 0.172f, 1);
        bg3 = ImVec4(0.190f, 0.200f, 0.228f, 1); txt = ImVec4(0.905f, 0.915f, 0.935f, 1); txtDis = ImVec4(0.50f, 0.53f, 0.58f, 1);
        border = ImVec4(0.22f, 0.23f, 0.26f, 0.70f); btn = ImVec4(0.17f, 0.18f, 0.205f, 1); btnH = ImVec4(0.22f, 0.24f, 0.28f, 1);
        popup = ImVec4(0.12f, 0.125f, 0.142f, 0.98f); tableHdr = ImVec4(0.13f, 0.137f, 0.155f, 1); tableLight = ImVec4(0.17f, 0.18f, 0.20f, 1);
        scroll = ImVec4(0.26f, 0.28f, 0.32f, 1); scrollH = ImVec4(0.33f, 0.35f, 0.40f, 1);
    }
    ImVec4* c = ImGui::GetStyle().Colors;
    c[ImGuiCol_Text] = txt;
    c[ImGuiCol_TextDisabled] = txtDis;
    c[ImGuiCol_WindowBg] = bg1;
    c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_PopupBg] = popup;
    c[ImGuiCol_Border] = border;
    c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_FrameBg] = bg2;
    c[ImGuiCol_FrameBgHovered] = bg3;
    c[ImGuiCol_FrameBgActive] = mixc(bg3, acc, 0.12f);
    c[ImGuiCol_TitleBg] = bg0;
    c[ImGuiCol_TitleBgActive] = bg0;
    c[ImGuiCol_TitleBgCollapsed] = bg0;
    c[ImGuiCol_MenuBarBg] = bg0;
    c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ScrollbarGrab] = scroll;
    c[ImGuiCol_ScrollbarGrabHovered] = scrollH;
    c[ImGuiCol_ScrollbarGrabActive] = acc;
    c[ImGuiCol_CheckMark] = accH;
    c[ImGuiCol_SliderGrab] = acc;
    c[ImGuiCol_SliderGrabActive] = accH;
    c[ImGuiCol_Button] = btn;
    c[ImGuiCol_ButtonHovered] = btnH;
    c[ImGuiCol_ButtonActive] = acc;
    c[ImGuiCol_Header] = alpha(acc, light ? 0.16f : 0.22f);
    c[ImGuiCol_HeaderHovered] = alpha(acc, light ? 0.24f : 0.32f);
    c[ImGuiCol_HeaderActive] = alpha(acc, light ? 0.34f : 0.45f);
    c[ImGuiCol_Separator] = alpha(border, 1.0f);
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
    c[ImGuiCol_TabDimmedSelectedOverline] = mixc(bg3, txtDis, 0.5f);
    c[ImGuiCol_DockingPreview] = alpha(acc, 0.55f);
    c[ImGuiCol_DockingEmptyBg] = bg0;
    c[ImGuiCol_PlotLines] = acc;
    c[ImGuiCol_PlotHistogram] = acc;
    c[ImGuiCol_TableHeaderBg] = tableHdr;
    c[ImGuiCol_TableBorderStrong] = alpha(border, 1.0f);
    c[ImGuiCol_TableBorderLight] = tableLight;
    c[ImGuiCol_TableRowBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_TableRowBgAlt] = light ? ImVec4(0, 0, 0, 0.03f) : ImVec4(1, 1, 1, 0.025f);
    c[ImGuiCol_TextSelectedBg] = alpha(acc, 0.35f);
    c[ImGuiCol_InputTextCursor] = light ? ImVec4(0.05f, 0.08f, 0.20f, 1) : ImVec4(0.80f, 0.88f, 1.00f, 1.00f);   // visible caret
    c[ImGuiCol_NavCursor] = acc;
    c[ImGuiCol_ModalWindowDimBg] = light ? ImVec4(0.2f, 0.22f, 0.26f, 0.35f) : ImVec4(0, 0, 0, 0.55f);
    ImGui::GetStyle().FrameBorderSize = hc ? 1.0f : 0.0f;   // high contrast: outlined fields (styleMetrics keeps it)

    if (ImPlot::GetCurrentContext()) {
        ImPlotStyle& ps = ImPlot::GetStyle();
        if (light) ImPlot::StyleColorsLight(); else ImPlot::StyleColorsDark();
        ps.PlotBorderSize = 0;
        ps.PlotPadding = ImVec2(10, 10);
        ps.LegendPadding = ImVec2(8, 6);
        ps.MinorAlpha = light ? 0.18f : 0.12f;
        ps.Colors[ImPlotCol_PlotBg] = light ? ImVec4(1, 1, 1, 1) : hc ? ImVec4(0, 0, 0, 1) : ImVec4(0.085f, 0.09f, 0.10f, 1.0f);
        ps.Colors[ImPlotCol_FrameBg] = ImVec4(0, 0, 0, 0);
        ps.Colors[ImPlotCol_LegendBg] = light ? ImVec4(1, 1, 1, 0.92f) : ImVec4(0.11f, 0.115f, 0.13f, 0.92f);
        ps.Colors[ImPlotCol_AxisGrid] = light ? ImVec4(0, 0, 0, 0.10f) : ImVec4(1, 1, 1, hc ? 0.16f : 0.08f);
    }
}

ImVec4 themeLabelText() {
    return g_theme == 1 ? ImVec4(0.27f, 0.29f, 0.34f, 1) : g_theme == 2 ? ImVec4(0.92f, 0.93f, 0.95f, 1) : ImVec4(0.72f, 0.75f, 0.80f, 1);
}
ImVec4 themeSelectedField() {
    const ImVec4 a = col4(Pal::Accent);
    const ImVec4 base = g_theme == 1 ? ImVec4(1, 1, 1, 1) : ImVec4(0.10f, 0.11f, 0.13f, 1);
    const float t = g_theme == 1 ? 0.22f : 0.38f;
    return ImVec4(base.x + (a.x - base.x) * t, base.y + (a.y - base.y) * t, base.z + (a.z - base.z) * t, 1);
}
