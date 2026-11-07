// C31 Results tab "UDQ & actions": the user-defined quantities of the run (FU* / WU*:well / GU*:group vectors from the
// UDQ keyword) against time with the ACTIONX firings marked, the action log, and the deck logic (UDQ definitions, action
// conditions and bodies) as the ResSim engine understood it.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "plotkit.h"
#include "theme.h"
#include "ressim/udq.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>

using namespace ressim;

void App::udqTab(const RunResults& R, double xs, const char* xl) {
    const float fs = ImGui::GetFontSize();
    const auto& S = R.summary;
    // all UDQ vector names of the run (first appearance order)
    std::vector<std::string> keys;
    {
        std::set<std::string> seen;
        for (const auto& p : S) for (const auto& u : p.udq) if (seen.insert(u.first).second) keys.push_back(u.first);
    }
    static std::set<std::string> sel;
    static char filter[48] = "";
    static bool markActions = true, logY = false;
    if (sel.empty() && !keys.empty()) for (const auto& k : keys) if (k.rfind("FU", 0) == 0 && sel.size() < 4) sel.insert(k);
    // ---------------------------------------------------------------- left: vector list + logic
    ImGui::BeginChild("udqlist", ImVec2(fs * 17, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    ImGui::TextDisabled("UDQ vectors (%d)", int(keys.size()));
    ImGui::SetNextItemWidth(-1); ImGui::InputTextWithHint("##udqf", "filter (e.g. WU*, :P1)", filter, sizeof filter);
    if (ImGui::SmallButton("All FU")) { sel.clear(); for (const auto& k : keys) if (k.rfind("FU", 0) == 0) sel.insert(k); }
    ImGui::SameLine(); if (ImGui::SmallButton("None")) sel.clear();
    for (const auto& k : keys) {
        if (filter[0]) {
            std::string f = filter, kk = k;
            for (char& c : f) c = char(std::toupper(static_cast<unsigned char>(c)));
            if (f.find('*') != std::string::npos ? !udq::wildMatch(f, kk) : kk.find(f) == std::string::npos) continue;
        }
        bool on = sel.count(k) > 0;
        if (ImGui::Checkbox(k.c_str(), &on)) { if (on) sel.insert(k); else sel.erase(k); }
    }
    if (keys.empty()) ui::hint("No UDQ values in this run (UDQ keyword without DEFINE / ASSIGN, or a run by another engine).");
    ImGui::Separator();
    ImGui::Checkbox("Mark actions", &markActions);
    ImGui::Checkbox("Log scale", &logY);
    ImGui::Separator();
    ImGui::TextDisabled("Deck logic");
    int ia = 0;
    for (const auto& e : R.input.extraKeywords) {
        if (e.name == "UDQ") {
            udq::Set set; std::string err;
            set.add(e.body, e.step, err);
            if (ImGui::TreeNode(&e, "UDQ (report step %d): %d records", e.step + 1, int(set.definitions().size()))) {
                for (const auto& d : set.definitions()) {
                    const std::string u = set.units(d.name);
                    if (d.kind == 0) ImGui::TextWrapped("DEFINE %s = %s%s", d.name.c_str(), d.expr.text.c_str(), u.empty() ? "" : ("  [" + u + "]").c_str());
                    else ImGui::TextWrapped("ASSIGN %s%s%s = %g", d.name.c_str(), d.selector.empty() ? "" : " ", d.selector.c_str(), d.value);
                }
                if (!err.empty()) ImGui::TextColored(col4(Pal::Warn), "%s", err.c_str());
                ImGui::TreePop();
            }
        } else if (e.name == "ACTIONX") {
            udq::Action a; std::string err;
            const bool ok = udq::parseAction(e.body, e.step, a, err);
            int fired = 0;
            for (const auto& ev : R.actions) if (ev.name == a.name) ++fired;
            ImGui::PushID(ia++);
            if (ImGui::TreeNode("act", "ACTIONX %s (from step %d) - fired %d / %d", a.name.c_str(), e.step + 1, fired, a.maxRuns)) {
                if (!ok) ImGui::TextColored(col4(Pal::Warn), "%s", err.c_str());
                ImGui::TextDisabled("when");
                ImGui::TextWrapped("%s", a.conditionText.c_str());
                if (a.minWaitDays > 0) ImGui::TextDisabled("min. interval %g days", a.minWaitDays);
                ImGui::TextDisabled("then");
                ImGui::TextWrapped("%s", a.body.c_str());
                ImGui::TreePop();
            }
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    // ---------------------------------------------------------------- right: plot + action log
    ImGui::BeginChild("udqplot", ImVec2(0, 0));
    const float logH = R.actions.empty() ? 0.0f : std::min(fs * 9.0f, fs * (2.5f + 1.3f * float(R.actions.size())));
    const ImVec2 av = ImGui::GetContentRegionAvail();
    fig::recorder().beginFigure("run.udq", 1, 1);
    if (ImPlot::BeginPlot("User-defined quantities", ImVec2(-1, av.y - logH))) {
        recPanel("User-defined quantities", xl, "value", "");
        ImPlot::SetupAxes(xl, "value", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        if (logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        int ci = 0;
        for (const auto& k : keys) {
            if (!sel.count(k)) continue;
            std::vector<double> x, y;
            for (const auto& p : S)
                for (const auto& u : p.udq)
                    if (u.first == k) { if (std::isfinite(u.second) && (!logY || u.second > 0)) { x.push_back(p.time * xs); y.push_back(u.second); } break; }
            const ImU32 col = seriesColor(ci++);
            if (x.empty()) continue;
            if (fig::recorder().recording()) recSeries(k.c_str(), x.data(), y.data(), int(x.size()), col, 2.0f);
            ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = 2.0f * g_plotLineScale;
            ImPlot::PlotLine(k.c_str(), x.data(), y.data(), int(x.size()), sp);
        }
        if (markActions && !R.actions.empty()) {
            ImPlot::PushPlotClipRect();
            const ImPlotRect lim = ImPlot::GetPlotLimits();
            for (const auto& ev : R.actions) {
                const double t = ev.t * xs;
                pk::line(t, lim.Y.Min, t, lim.Y.Max, IM_COL32(220, 60, 60, 200), 1.5f, "action", 1);
                pk::text(t, lim.Y.Max, ev.name.c_str(), IM_COL32(220, 60, 60, 255), ImVec2(3, 2), 0);
            }
            ImPlot::PopPlotClipRect();
        }
        recEnd();
        ImPlot::EndPlot();
    }
    fig::recorder().endFigure();
    if (!R.actions.empty() && ImGui::BeginTable("##actlog", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY, ImVec2(-1, logH - fs * 0.5f))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Time (days)"); ImGui::TableSetupColumn("End of report step"); ImGui::TableSetupColumn("Action"); ImGui::TableSetupColumn("Wells ('?')");
        ImGui::TableHeadersRow();
        for (const auto& ev : R.actions) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("%.4g", ev.t);
            ImGui::TableNextColumn(); ImGui::Text("%d", ev.step + 1);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(ev.name.c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(ev.wells.empty() ? "-" : ev.wells.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::EndChild();
}
