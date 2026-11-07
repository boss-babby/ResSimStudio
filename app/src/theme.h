#pragma once
#include "imgui.h"

void applyTheme(float dpiScale);
// C39 colour themes: 0 dark (default), 1 light, 2 high contrast (dark); accent 0 blue, 1 teal, 2 violet, 3 orange, 4 green.
// Safe to call at any time (only colours change).
void applyColors(int theme, int accent);
int themeId();
int accentId();
const char* const* themeNames();    // 3 entries
const char* const* accentNames();   // 5 entries
ImVec4 accentColor(int accent);
// C40 background of the reservoir 3-D views: 0 dark (default), 1 follow the theme, 2 light grey, 3 white, 4 black
void set3dBackground(int mode);
int background3dMode();
const float* background3d();                 // rgb for the current mode / theme
const char* const* background3dNames();      // 5 entries
// C39 design tokens (theme-aware)
ImVec4 themeLabelText();      // property-grid labels: one step below body text
ImVec4 themeSelectedField();  // background of a highlighted / active input field
// surfaces used by custom-drawn chrome (status bar, toasts, cards)
ImVec4 themeStatusBg();
ImVec4 themeCardBg();
bool themeIsLight();
// C41 (G5) design tokens: spacing scale 4 / 8 / 12 / 16 (x interface size), one corner radius, surface elevation
// levels (tints of the window background) and three text levels (0 body, 1 labels, 2 hints / muted)
ImVec4 themeSurface(int level);
ImVec4 themeTextLevel(int level);
namespace tok {
float sp(int step);                 // 1..4 -> 4, 8, 12, 16 px x uiScale
inline constexpr float kRadius = 6.0f;
}
// Interface size (0.6-1.6, 1 = default): fonts and all paddings/spacings; safe to call every frame-change
void setUiScale(float uiScale);
float uiScale();

// Fonts loaded by applyTheme (may be null -> default font)
struct Fonts { ImFont* ui = nullptr; ImFont* bold = nullptr; ImFont* mono = nullptr; ImFont* title = nullptr; };
extern Fonts g_fonts;

// Accent palette shared by plots and widgets
namespace Pal {
constexpr ImU32 Oil = IM_COL32(46, 204, 113, 255);
constexpr ImU32 Water = IM_COL32(52, 152, 219, 255);
constexpr ImU32 Gas = IM_COL32(231, 76, 60, 255);
constexpr ImU32 Pressure = IM_COL32(241, 196, 15, 255);
inline ImU32 Accent = IM_COL32(64, 132, 246, 255);   // C39: follows the accent setting (applyColors)
inline ImU32 Muted = IM_COL32(140, 146, 160, 255);   // C39: darker in the light theme
inline ImU32 Ok = IM_COL32(72, 199, 142, 255);       // C39: status colours are darkened in the light theme (contrast >= 3:1)
inline ImU32 Warn = IM_COL32(255, 179, 64, 255);
inline ImU32 Err = IM_COL32(255, 92, 92, 255);
}
inline ImVec4 col4(ImU32 c) { return ImGui::ColorConvertU32ToFloat4(c); }
