// C40 navigation layer: deep links into any tab of the application and an index of every tab the user has seen.
//
// All tab items go through ui::tabItem (a drop-in replacement for ImGui::BeginTabItem). It
//   * selects a tab when a navigation request names it (requests name the whole chain of nested tabs, outermost first,
//     and stay alive for a few frames so nested tab bars can follow one another), and
//   * records every tab it draws together with the workspace / section and the chain of parent tabs, so the command
//     palette can jump straight to it later (the index is persisted between sessions).
#pragma once
#include "imgui.h"
#include <string>
#include <vector>

namespace nav {

struct TabEntry {
    std::string label;                 // visible label (icons and ##ids stripped)
    std::vector<std::string> chain;    // parent tabs, outermost first (normalised labels)
    int workspace = 0, section = -1;   // where it was drawn
    std::string window;                // ImGui window the tab bar lives in (for display only)
};

// label without icon glyphs, "##id" / "###id" suffixes and surrounding blanks
std::string clean(const char* label);

// called by the App at the start of each frame body, before any panel is drawn
void beginFrame(int workspace, int section);

// ask for a chain of tabs (outermost first) to be selected during the next frames
void request(const std::vector<std::string>& chain);
bool pending();

// tab index: all tabs seen in this and earlier sessions
const std::vector<TabEntry>& index();
void loadIndex(const std::string& file);
void saveIndex(const std::string& file);
bool indexChanged();

}  // namespace nav

namespace ui {
// drop-in replacement for ImGui::BeginTabItem (same arguments)
bool tabItem(const char* label, bool* p_open = nullptr, ImGuiTabItemFlags flags = 0);
}  // namespace ui
