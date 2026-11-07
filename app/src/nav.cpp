// C40 navigation layer (see nav.h)
#include "nav.h"
#include "i18n.h"
#include "imgui_internal.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace nav {
namespace {
std::vector<TabEntry> g_index;
bool g_changed = false;
int g_ws = 0, g_sec = -1;
std::vector<std::string> g_open;        // label of the selected tab at each tab-bar depth in this frame
std::vector<std::string> g_req;         // requested chain
int g_reqFrames = 0;                    // frames left for the request
size_t g_reqDone = 0;                   // chain elements already confirmed

std::string lower(std::string s) { for (char& c : s) c = char(std::tolower((unsigned char)c)); return s; }

void remember(const std::string& lab, const std::vector<std::string>& chain) {
    if (lab.empty() || lab.size() > 80) return;
    for (TabEntry& e : g_index)
        if (e.label == lab && e.chain == chain) {
            if (e.workspace != g_ws || e.section != g_sec) { e.workspace = g_ws; e.section = g_sec; g_changed = true; }
            return;
        }
    TabEntry e;
    e.label = lab; e.chain = chain; e.workspace = g_ws; e.section = g_sec;
    if (ImGuiWindow* w = ImGui::GetCurrentWindow()) {
        std::string n = w->Name ? w->Name : "";
        const size_t k = n.find("###"); if (k != std::string::npos) n.resize(k);
        const size_t s = n.find('/'); if (s != std::string::npos) n.resize(s);   // child windows: "Parent/child_id"
        e.window = clean(n.c_str());
    }
    g_index.push_back(std::move(e));
    g_changed = true;
}
}  // namespace

std::string clean(const char* label) {
    std::string s = label ? label : "";
    const size_t h = s.find("##"); if (h != std::string::npos) s.resize(h);
    // icon glyphs are private-use code points (3-byte UTF-8 sequences starting with 0xEE / 0xEF); drop them
    std::string o; o.reserve(s.size());
    for (size_t i = 0; i < s.size();) {
        const unsigned char c = (unsigned char)s[i];
        if ((c == 0xEE || c == 0xEF) && i + 2 < s.size()) { i += 3; continue; }
        o += s[i]; ++i;
    }
    size_t a = 0, b = o.size();
    while (a < b && std::isspace((unsigned char)o[a])) ++a;
    while (b > a && std::isspace((unsigned char)o[b - 1])) --b;
    return o.substr(a, b - a);
}

void beginFrame(int workspace, int section) {
    g_ws = workspace; g_sec = section;
    g_open.clear();
    if (g_reqFrames > 0 && --g_reqFrames == 0) { g_req.clear(); g_reqDone = 0; }
}

void request(const std::vector<std::string>& chain) {
    g_req.clear();
    for (const std::string& c : chain) g_req.push_back(lower(clean(c.c_str())));
    g_reqDone = 0;
    g_reqFrames = g_req.empty() ? 0 : 12 + 4 * int(g_req.size());
}
bool pending() { return g_reqFrames > 0; }

const std::vector<TabEntry>& index() { return g_index; }
bool indexChanged() { return g_changed; }

// file format: one tab per line   workspace \t section \t window \t chain (joined by " > ") \t label
void loadIndex(const std::string& file) {
    std::ifstream f(file, std::ios::binary);
    for (std::string line; std::getline(f, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::vector<std::string> p; std::stringstream ss(line);
        for (std::string t; std::getline(ss, t, '\t');) p.push_back(t);
        if (p.size() != 5 || p[4].empty()) continue;
        TabEntry e;
        e.workspace = std::atoi(p[0].c_str()); e.section = std::atoi(p[1].c_str()); e.window = p[2]; e.label = p[4];
        for (size_t a = 0; !p[3].empty();) {
            const size_t k = p[3].find(" > ", a);
            e.chain.push_back(p[3].substr(a, k == std::string::npos ? std::string::npos : k - a));
            if (k == std::string::npos) break;
            a = k + 3;
        }
        bool dup = false;
        for (const TabEntry& x : g_index) if (x.label == e.label && x.chain == e.chain) { dup = true; break; }
        if (!dup) g_index.push_back(std::move(e));
    }
    g_changed = false;
}

void saveIndex(const std::string& file) {
    std::ofstream f(file, std::ios::binary);
    for (const TabEntry& e : g_index) {
        std::string ch;
        for (size_t i = 0; i < e.chain.size(); ++i) ch += (i ? " > " : "") + e.chain[i];
        f << e.workspace << '\t' << e.section << '\t' << e.window << '\t' << ch << '\t' << e.label << '\n';
    }
    g_changed = false;
}

}  // namespace nav

namespace ui {
bool tabItem(const char* label, bool* p_open, ImGuiTabItemFlags flags) {
    using namespace nav;
    ImGuiContext& g = *GImGui;
    const int depth = std::max(0, g.CurrentTabBarStack.Size - 1);
    const std::string lab = clean(label);
    if (g_reqFrames > 0 && size_t(depth) < g_req.size()) {
        // a chain element matches when the label equals it and the parents drawn this frame match the earlier elements
        bool parentsOk = true;
        for (int d = 0; d < depth && parentsOk; ++d) parentsOk = d < int(g_open.size()) && lower(g_open[size_t(d)]) == g_req[size_t(d)];
        if (parentsOk && lower(lab) == g_req[size_t(depth)]) flags |= ImGuiTabItemFlags_SetSelected;
    }
    const std::string shown = i18n::trLabel(label);   // C45 translated, ImGui id stays the English label
    const bool open = ImGui::BeginTabItem(shown.c_str(), p_open, flags);
    if (open) {
        if (int(g_open.size()) <= depth) g_open.resize(size_t(depth) + 1);
        g_open[size_t(depth)] = lab;
        std::vector<std::string> chain(g_open.begin(), g_open.begin() + depth);
        remember(lab, chain);
        if (g_reqFrames > 0 && size_t(depth) + 1 == g_req.size() && lower(lab) == g_req.back()) { g_req.clear(); g_reqFrames = 0; g_reqDone = 0; }
    } else if (!(ImGui::GetCurrentContext()->CurrentTabBar == nullptr)) {
        // closed tabs are indexed too (with the chain of the currently open parents) so the palette knows them
        std::vector<std::string> chain(g_open.begin(), g_open.begin() + std::min<size_t>(size_t(depth), g_open.size()));
        if (chain.size() == size_t(depth)) remember(lab, chain);
    }
    return open;
}
}  // namespace ui
