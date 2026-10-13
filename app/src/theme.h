#pragma once
#include "imgui.h"

void applyTheme(float dpiScale);
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
constexpr ImU32 Accent = IM_COL32(64, 132, 246, 255);
constexpr ImU32 Muted = IM_COL32(140, 146, 160, 255);
constexpr ImU32 Ok = IM_COL32(72, 199, 142, 255);
constexpr ImU32 Warn = IM_COL32(255, 179, 64, 255);
constexpr ImU32 Err = IM_COL32(255, 92, 92, 255);
}
inline ImVec4 col4(ImU32 c) { return ImGui::ColorConvertU32ToFloat4(c); }
