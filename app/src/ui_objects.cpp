// C40 object tree (Case Explorer): the model as objects instead of input sections.
//
//   Field > groups (GRUPTREE) > wells        click: select + open the Wells editor; double-click: focus in 3-D
//   Faults, local grids, aquifers, regions    click: open the editor that owns them
//   Runs (Result Viewer / run history)        click: open in Results; context menu: reference / remove
//   Observed data                             click: Studies
// Context menus on wells: edit, focus in 3-D, rename, duplicate, delete. Every change goes through markDirty(), so it is
// undoable (Ctrl+Z) like any other edit; destructive actions therefore do not ask for confirmation.
#include "app.h"
#include "icons.h"
#include "theme.h"
#include "widgets.h"
#include "i18n.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <functional>
#include <map>
#include <set>

using namespace ressim;

namespace {
std::string lowerStr(std::string s) { for (char& c : s) c = char(std::tolower((unsigned char)c)); return s; }
bool matches(const std::string& text, const char* filter) {
    if (!filter || !*filter) return true;
    return lowerStr(text).find(lowerStr(filter)) != std::string::npos;
}
ImU32 wellColor(WellType t) {
    return t == WellType::Producer ? Pal::Oil : t == WellType::WaterInjector ? Pal::Water : Pal::Gas;
}
std::string uniqueWellName(const Case& c, const std::string& base) {
    auto taken = [&](const std::string& n) { for (const auto& w : c.wells) if (w.name == n) return true; return false; };
    if (!taken(base)) return base;
    for (int k = 2; k < 10000; ++k) { const std::string n = base + "_" + std::to_string(k); if (!taken(n)) return n; }
    return base + "_X";
}
// every place the case refers to a well by name (schedule controls, WELPI, observations)
void renameWellRefs(Case& c, const std::string& from, const std::string& to) {
    for (ReportStep& s : c.schedule) {
        for (WellControl& wc : s.controls) if (wc.well == from) wc.well = to;
        for (auto& p : s.welpi) if (p.first == from) p.first = to;
        for (auto& p : s.wpolymer) if (p.well == from) p.well = to;
    }
    for (Observation& o : c.observations) if (o.well == from) o.well = to;
}
}  // namespace

void App::objectTree() {
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##objfilter", ICON_SEARCH "  filter objects (wells, groups, faults ...)", objFilter_, sizeof objFilter_);
    const bool filtering = objFilter_[0] != 0;
    const ImGuiTreeNodeFlags branch = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick | ImGuiTreeNodeFlags_SpanAvailWidth;
    const ImGuiTreeNodeFlags leaf = ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanAvailWidth;
    auto openIf = [&](bool def) { ImGui::SetNextItemOpen(def || filtering, filtering ? ImGuiCond_Always : ImGuiCond_Once); };
    auto count = [](const char* name, size_t n) { return std::string(TR(name)) + "  (" + std::to_string(n) + ")"; };   // C48: translated
    std::function<void()> deferred;   // structural edits run after the tree is drawn

    // ------------------------------------------------------------------ wells by group
    {
        // group hierarchy: GRUPTREE pairs plus the groups the wells name; everything hangs below FIELD
        std::map<std::string, std::vector<std::string>> kids;
        std::map<std::string, std::vector<int>> wellsOf;
        std::set<std::string> known{"FIELD"};
        for (const auto& ct : case_.groupTree) { kids[ct.second].push_back(ct.first); known.insert(ct.first); known.insert(ct.second); }
        for (int w = 0; w < int(case_.wells.size()); ++w) {
            std::string g = case_.wells[size_t(w)].group.empty() ? std::string("FIELD") : case_.wells[size_t(w)].group;
            wellsOf[g].push_back(w);
            if (!known.count(g)) { known.insert(g); kids["FIELD"].push_back(g); }
        }
        for (const std::string& g : known) {   // orphan groups (parent never defined) attach to FIELD
            if (g == "FIELD") continue;
            bool hasParent = false;
            for (const auto& ct : case_.groupTree) if (ct.first == g) { hasParent = true; break; }
            if (!hasParent && std::find(kids["FIELD"].begin(), kids["FIELD"].end(), g) == kids["FIELD"].end()) kids["FIELD"].push_back(g);
        }
        // does a group subtree contain a matching well / group?
        std::function<bool(const std::string&, int)> anyMatch = [&](const std::string& g, int depth) -> bool {
            if (depth > 32) return false;
            if (matches(g, objFilter_)) return true;
            for (int w : wellsOf[g]) if (matches(case_.wells[size_t(w)].name, objFilter_)) return true;
            for (const std::string& k : kids[g]) if (anyMatch(k, depth + 1)) return true;
            return false;
        };
        auto wellRow = [&](int w) {
            const WellInput& W = case_.wells[size_t(w)];
            if (!matches(W.name, objFilter_) && !matches(W.group, objFilter_)) return;
            ImGui::PushID(w);
            const bool sel = section_ == Section::Wells && selWell_ == w;
            const ImVec2 p = ImGui::GetCursorScreenPos();
            char lb[160];
            std::snprintf(lb, sizeof lb, "    %s   (%d, %d)  %d conn.", W.name.c_str(), W.i + 1, W.j + 1, int(W.completions.size()));
            ImGui::TreeNodeEx(TRL("##w"), leaf | (sel ? ImGuiTreeNodeFlags_Selected : 0), "%s", lb);
            const float r = fs * 0.24f;
            ImGui::GetWindowDrawList()->AddCircleFilled(ImVec2(p.x + ImGui::GetTreeNodeToLabelSpacing() + r, p.y + ImGui::GetTextLineHeight() * 0.5f + 1), r, wellColor(W.type));
            if (ImGui::IsItemClicked(ImGuiMouseButton_Left)) { workspace_ = Workspace::Model; section_ = Section::Wells; selWell_ = w; }
            if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && viewOk_) { view_.settings().focusWell = w; view_.focusWellCamera(w); }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal))
                ImGui::SetTooltip("%s  -  %s, %s %s\n%s", W.name.c_str(),
                                  TR(W.type == WellType::Producer ? "producer" : W.type == WellType::WaterInjector ? "water injector" : "gas injector"),
                                  TR("group"), W.group.c_str(), TR("click: edit   double-click: focus in 3-D   right-click: more"));
            if (ImGui::BeginPopupContextItem("##wctx")) {
                if (ImGui::MenuItem(TRL(ICON_PENCIL "  Edit"))) { workspace_ = Workspace::Model; section_ = Section::Wells; selWell_ = w; }
                if (ImGui::MenuItem(TRL(ICON_CROSSHAIR "  Focus in 3-D"), nullptr, false, viewOk_)) { view_.settings().focusWell = w; view_.focusWellCamera(w); }
                if (ImGui::MenuItem(TRL(ICON_EYE "  Show all wells"), nullptr, false, viewOk_)) view_.settings().focusWell = -1;
                ImGui::Separator();
                if (ImGui::BeginMenu(TRL(ICON_PENCIL "  Rename"))) {
                    static char nb[64];
                    if (ImGui::IsWindowAppearing()) std::snprintf(nb, sizeof nb, "%s", W.name.c_str());
                    ImGui::SetNextItemWidth(fs * 10);
                    const bool ok = ImGui::InputText("##rn", nb, sizeof nb, ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsUppercase);
                    ImGui::SameLine();
                    if ((ImGui::Button(TRL("OK")) || ok) && nb[0]) {
                        const std::string from = W.name, to = uniqueWellName(case_, nb);
                        deferred = [this, w, from, to] {
                            case_.wells[size_t(w)].name = to;
                            renameWellRefs(case_, from, to);
                            markDirty();
                            notify("Well " + from + " renamed to " + to + " (schedule, observations updated; Ctrl+Z undoes)");
                        };
                        ImGui::CloseCurrentPopup();
                    }
                    ImGui::TextDisabled("%s", TR("Updates the schedule and observed data. Verbatim keywords\n(ACTIONX, UDQ ...) keep the old name."));
                    ImGui::EndMenu();
                }
                if (ImGui::MenuItem(TRL(ICON_COPY "  Duplicate (with its schedule)"))) {
                    deferred = [this, w] {
                        WellInput c = case_.wells[size_t(w)];
                        const std::string from = c.name;
                        c.name = uniqueWellName(case_, from);
                        c.i = std::min(c.i + 1, std::max(0, case_.grid.nx - 1));
                        for (Completion& k : c.completions) k.i = std::min(k.i + 1, std::max(0, case_.grid.nx - 1));
                        for (ReportStep& s : case_.schedule) {
                            const size_t n = s.controls.size();
                            for (size_t q = 0; q < n; ++q) if (s.controls[q].well == from) { WellControl x = s.controls[q]; x.well = c.name; s.controls.push_back(x); }
                        }
                        const std::string nn = c.name;
                        case_.wells.push_back(std::move(c));
                        selWell_ = int(case_.wells.size()) - 1;
                        section_ = Section::Wells;
                        markDirty();
                        notify("Well " + nn + " added one cell east of " + from + " with the same controls");
                    };
                }
                if (ImGui::MenuItem(TRL(ICON_TRASH_2 "  Delete"))) {
                    deferred = [this, w] {
                        const std::string name = case_.wells[size_t(w)].name;
                        case_.wells.erase(case_.wells.begin() + w);
                        for (ReportStep& s : case_.schedule) {
                            s.controls.erase(std::remove_if(s.controls.begin(), s.controls.end(), [&](const WellControl& c) { return c.well == name; }), s.controls.end());
                            s.welpi.erase(std::remove_if(s.welpi.begin(), s.welpi.end(), [&](const auto& p) { return p.first == name; }), s.welpi.end());
                            s.wpolymer.erase(std::remove_if(s.wpolymer.begin(), s.wpolymer.end(), [&](const auto& p) { return p.well == name; }), s.wpolymer.end());
                        }
                        case_.observations.erase(std::remove_if(case_.observations.begin(), case_.observations.end(), [&](const Observation& o) { return o.well == name; }), case_.observations.end());
                        selWell_ = std::clamp(selWell_, 0, std::max(0, int(case_.wells.size()) - 1));
                        if (viewOk_) view_.settings().focusWell = -1;
                        markDirty();
                        notify("Well " + name + " deleted - Ctrl+Z restores it");
                    };
                }
                ImGui::EndPopup();
            }
            ImGui::PopID();
        };
        std::function<void(const std::string&, int)> groupNode = [&](const std::string& g, int depth) {
            if (depth > 32 || !anyMatch(g, 0)) return;
            size_t nw = wellsOf[g].size();
            openIf(depth <= 1);
            ImGui::PushID(g.c_str());
            const bool open = ImGui::TreeNodeEx(TRL("##g"), branch, "%s  %s%s", depth == 0 ? ICON_FOLDER_TREE : ICON_FOLDER, g.c_str(),
                                                nw ? ("   " + std::to_string(nw) + (nw == 1 ? " well" : " wells")).c_str() : "");
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) ImGui::SetTooltip("group %s (GRUPTREE / WELSPECS)\nGroup controls: Model > Schedule", g.c_str());
            if (open) {
                for (const std::string& k : kids[g]) if (k != g) groupNode(k, depth + 1);
                for (int w : wellsOf[g]) wellRow(w);
                ImGui::TreePop();
            }
            ImGui::PopID();
        };
        openIf(true);
        if (ImGui::TreeNodeEx(TRL("##wells"), branch, "%s", count(ICON_TARGET "  Wells", case_.wells.size()).c_str())) {
            if (case_.wells.empty()) {
                ImGui::TextDisabled("%s", TR("No wells yet."));
                if (ImGui::SmallButton(TRL(ICON_PLUS " Add wells"))) { workspace_ = Workspace::Model; section_ = Section::Wells; }
            } else groupNode("FIELD", 0);
            ImGui::TreePop();
        }
    }

    // ------------------------------------------------------------------ grid objects
    auto simpleList = [&](const char* id, const std::string& title, size_t n, Section sec, const std::function<std::string(size_t)>& label, const char* empty) {
        bool any = !filtering;
        for (size_t i = 0; i < n && !any; ++i) any = matches(label(i), objFilter_);
        if (!any) return;
        openIf(false);
        ImGui::PushID(id);
        if (ImGui::TreeNodeEx(TRL("##n"), branch, "%s", title.c_str())) {
            if (n == 0) ImGui::TextDisabled("%s", TR(empty));
            for (size_t i = 0; i < n; ++i) {
                const std::string l = label(i);
                if (!matches(l, objFilter_)) continue;
                ImGui::PushID(int(i));
                ImGui::TreeNodeEx(TRL("##i"), leaf, "%s", l.c_str());
                if (ImGui::IsItemClicked()) { workspace_ = Workspace::Model; section_ = sec; }
                ImGui::PopID();
            }
            ImGui::TreePop();
        }
        ImGui::PopID();
    };
    const GridInput& G = case_.grid;
    simpleList("faults", count(ICON_SPLIT "  Faults", G.faults.size()), G.faults.size(), Section::Grid,
               [&](size_t i) { char b[96]; std::snprintf(b, sizeof b, "%s   MULTFLT %.3g   %d faces", G.faults[i].name.c_str(), G.faults[i].mult, int(G.faults[i].faces.size())); return std::string(b); },
               "none (Grid & rock: FAULTS / MULTFLT, or convert to corner-point)");
    simpleList("lgrs", count(ICON_GRID_3X3 "  Local grids", G.lgrs.size()), G.lgrs.size(), Section::Grid,
               [&](size_t i) { return G.lgrs[i].name; }, "none (Grid & rock: Local grid refinement)");
    simpleList("aqu", count(ICON_WAVES "  Aquifers", case_.aquifers.size()), case_.aquifers.size(), Section::Aquifers,
               [&](size_t) { return std::string("aquifer"); }, "none (Aquifers)");
    {
        struct Reg { const char* name; const std::vector<int>* v; Section sec; };
        const Reg regs[] = {{"SATNUM", &G.satnum, Section::SCAL}, {"PVTNUM", &G.pvtnum, Section::PVT}, {"EQLNUM", &G.eqlnum, Section::Init}, {"FIPNUM", &G.fipnum, Section::Init}};
        std::vector<std::string> rows; std::vector<Section> secs;
        for (const Reg& r : regs) {
            std::set<int> ids(r.v->begin(), r.v->end());
            rows.push_back(std::string(r.name) + "   " + (r.v->empty() ? std::string("1 region (default)") : std::to_string(ids.size()) + " regions"));
            secs.push_back(r.sec);
        }
        if (!filtering || std::any_of(rows.begin(), rows.end(), [&](const std::string& s) { return matches(s, objFilter_); })) {
            openIf(false);
            if (ImGui::TreeNodeEx(TRL("##regions"), branch, "%s", TR(ICON_LAYERS "  Regions"))) {
                for (size_t i = 0; i < rows.size(); ++i) {
                    if (!matches(rows[i], objFilter_)) continue;
                    ImGui::PushID(int(i));
                    ImGui::TreeNodeEx(TRL("##r"), leaf, "%s", rows[i].c_str());
                    if (ImGui::IsItemClicked()) { workspace_ = Workspace::Model; section_ = secs[i]; }
                    ImGui::PopID();
                }
                ImGui::TreePop();
            }
        }
    }

    // ------------------------------------------------------------------ runs and data
    {
        bool any = !filtering;
        for (const RVCase& c : rv_.cases) any = any || matches(c.name, objFilter_);
        if (any) {
            openIf(!rv_.cases.empty());
            if (ImGui::TreeNodeEx(TRL("##runs"), branch, "%s", count(ICON_HISTORY "  Runs", rv_.cases.size()).c_str())) {
                if (rv_.cases.empty()) ImGui::TextDisabled("%s", TR("Runs appear here after F5 (run history) or File > Open results."));
                for (int i = 0; i < int(rv_.cases.size()); ++i) {
                    const RVCase& c = rv_.cases[size_t(i)];
                    if (!matches(c.name, objFilter_)) continue;
                    ImGui::PushID(i);
                    const ImVec2 p = ImGui::GetCursorScreenPos();
                    std::string lab = "    " + c.name + (i == rv_.ref ? "   [reference]" : "");
                    ImGui::TreeNodeEx(TRL("##run"), leaf | (rv_.c3d == i ? ImGuiTreeNodeFlags_Selected : 0), "%s", lab.c_str());
                    ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x + ImGui::GetTreeNodeToLabelSpacing(), p.y + 3),
                                                              ImVec2(p.x + ImGui::GetTreeNodeToLabelSpacing() + fs * 0.5f, p.y + ImGui::GetTextLineHeight() - 1), c.color, 2.0f);
                    if (ImGui::IsItemClicked()) { rv_.c3d = i; rv_.step = 0; workspace_ = Workspace::Results; showResultViewer_ = true; }
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal) && !c.changes.empty()) {
                        ImGui::BeginTooltip(); ImGui::PushTextWrapPos(fs * 30);
                        ImGui::TextUnformatted((std::string(TR("Changes against the previous run:")) + "\n" + c.changes.substr(0, 600)).c_str());
                        ImGui::PopTextWrapPos(); ImGui::EndTooltip();
                    }
                    if (ImGui::BeginPopupContextItem("##rctx")) {
                        if (ImGui::MenuItem(TRL(ICON_CHART_LINE "  Open in the Result Viewer"))) { rv_.c3d = i; showResultViewer_ = true; }
                        if (ImGui::MenuItem(TRL(ICON_TARGET "  Use as reference for comparisons"))) rv_.ref = i;
                        if (ImGui::MenuItem(TRL(ICON_TRASH_2 "  Remove from the list"))) deferred = [this, i] {
                            if (i < int(rv_.cases.size())) rv_.cases.erase(rv_.cases.begin() + i);
                            const int nc = int(rv_.cases.size());
                            rv_.ref = std::clamp(rv_.ref, 0, std::max(0, nc - 1)); rv_.cmp = std::clamp(rv_.cmp, 0, std::max(0, nc - 1));
                            rv_.c3d = std::clamp(rv_.c3d, 0, std::max(0, nc - 1));
                        };
                        ImGui::EndPopup();
                    }
                    ImGui::PopID();
                }
                ImGui::TreePop();
            }
        }
    }
    if (!filtering || matches("observed data history", objFilter_)) {
        ImGui::TreeNodeEx(TRL("##obs"), leaf, "%s", count(ICON_CHART_LINE "  Observed data", case_.observations.size()).c_str());
        if (ImGui::IsItemClicked()) { workspace_ = Workspace::Studies; nav::request({"Observed"}); }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) ImGui::SetTooltip("%s", TR("History-matching targets (Studies > Observed data)"));
    }
    if (deferred) deferred();
}
