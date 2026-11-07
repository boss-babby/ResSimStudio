// C32 multisegment wells in the UI:
//   * Results -> Multisegment wells: node pressure / segment flow / friction and device drops along the tubing at any
//     accepted time step, heel-toe pressure history, segment table (ResSim runs with WELSEGS)
//   * Model -> Wells -> Multisegment well: WELSEGS / COMPSEGS (+ WSEGSICD / WSEGAICD / WSEGVALV) generator from the
//     completions of the selected well (one tubing segment per connection, optional device segment per connection)
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include "plotkit.h"
#include "figure.h"
#include "ressim/msw.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

using namespace ressim;

namespace {
ImPlotSpec mswLine(int i, float w = 2.0f) { ImPlotSpec s; s.LineColor = col4(seriesColor(i)); s.LineWeight = w; return s; }
// cumulative node length along the tubing per segment position (top node = topLength)
std::vector<double> nodeLength(const msw::Well& W) {   // the parser stores absolute node lengths (INC converted)
    std::vector<double> L(W.segs.size(), 0.0);
    for (size_t p = 0; p < W.segs.size(); ++p) L[p] = W.segs[p].outlet <= 0 ? W.topLength : W.segs[p].length;
    return L;
}
std::string firstToken(const std::string& line) {
    std::istringstream is(line); std::string t; is >> t;
    if (t.size() >= 2 && (t.front() == '\'' || t.front() == '"')) t = t.substr(1, t.size() - 2);
    return t;
}
}  // namespace

void App::mswTab(const RunResults& R, double xs, const char* xl) {
    const float fs = ImGui::GetFontSize();
    if (R.msw.empty() || !R.mswRec || R.mswRec->empty()) { ui::hint("No multisegment-well results in this run (WELSEGS / COMPSEGS)."); return; }
    const bool field = R.input.units == UnitSystem::Field;
    const char* pu = field ? "psi" : "bar";
    const char* qu = field ? "rb/d" : "rm3/d";
    const char* lu = field ? "ft" : "m";
    static int sel = 0, rec = -1, xAxis = 0;
    static bool showFirst = true;
    sel = std::clamp(sel, 0, int(R.msw.size()) - 1);
    const msw::Well& W = R.msw[size_t(sel)];
    std::vector<int> idx;   // records of this well
    for (size_t k = 0; k < R.mswRec->size(); ++k) if ((*R.mswRec)[k].msw == sel) idx.push_back(int(k));
    // ---------------------------------------------------------------- left: selection
    ImGui::BeginChild("mswleft", ImVec2(fs * 16, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    ImGui::TextDisabled("Multisegment wells");
    for (int m = 0; m < int(R.msw.size()); ++m) if (ImGui::Selectable(R.msw[size_t(m)].name.c_str(), sel == m)) { sel = m; rec = -1; }
    ImGui::Separator();
    if (idx.empty()) { ImGui::TextDisabled("no records"); ImGui::EndChild(); return; }
    if (rec < 0 || rec >= int(idx.size())) rec = int(idx.size()) - 1;
    const auto& Rr = (*R.mswRec)[size_t(idx[size_t(rec)])];
    ImGui::Text("time %.4g %s", Rr.t * xs, xl);
    ImGui::SetNextItemWidth(-1); ImGui::SliderInt("##rec", &rec, 0, int(idx.size()) - 1, "step %d");
    static const char* xa[] = {"length along tubing", "segment depth", "segment id"};
    ImGui::SetNextItemWidth(-1); ImGui::Combo("##xa", &xAxis, xa, 3);
    ImGui::Checkbox("overlay the first record", &showFirst);
    ImGui::Separator();
    int nDev = 0, nShut = 0, nLim = 0, nPump = 0;
    for (const auto& s : W.segs) { nDev += s.device > 0; nShut += s.shut; nLim += (s.flimPh[0] > 0 || s.flimPh[1] > 0); nPump += s.pump; }
    ImGui::Text("%zu segments, %zu connections", W.segs.size(), W.conns.size());
    ImGui::Text("%d device segment(s), %d shut", nDev, nShut);
    if (nLim || nPump) ImGui::Text("%d flow limiter(s) (WSEGFLIM), %d pump(s) (WSEGPULL)", nLim, nPump);
    ImGui::Text("pressure drops: %s", W.accelOn ? "HFA (friction + acceleration)" : W.frictionOn ? "HF- (friction)" : "H-- (hydrostatic only)");
    ImGui::Text("flow model: %s%s", W.driftFlux ? "DF drift flux (gas slip)" : "HO homogeneous", W.driftFlux && W.df.owFactor > 0 ? " + oil-water slip" : "");
    ImGui::Text("tubing flash %s, segment storage %s", R.input.mswOpt.flash ? "on" : "off", R.input.mswOpt.storage ? "on" : "off");
    if (!W.warn.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextColored(col4(Pal::Warn), "%s", W.warn.c_str()); ImGui::PopTextWrapPos(); }
    ui::hint("Pressures are node pressures; drops are per segment (p(segment) - p(outlet) parts). Device segments: the drop "
             "across the ICD / valve between the reservoir branch and the tubing.");
    ImGui::EndChild();
    ImGui::SameLine();
    // ---------------------------------------------------------------- right: profiles + history + table
    ImGui::BeginChild("mswright", ImVec2(0, 0));
    const size_t ns = W.segs.size();
    const auto L = nodeLength(W);
    std::vector<double> x(ns), tub;   // tubing (branch 1) positions first
    for (size_t p = 0; p < ns; ++p) x[p] = xAxis == 0 ? L[p] : xAxis == 1 ? W.segs[p].depth : W.segs[p].id;
    std::vector<size_t> ord;          // tubing segments ordered by x
    for (size_t p = 0; p < ns; ++p) if (W.segs[p].branch == 1) ord.push_back(p);
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
    std::vector<size_t> dev;
    for (size_t p = 0; p < ns; ++p) if (W.segs[p].device > 0 || W.segs[p].branch != 1) dev.push_back(p);
    std::sort(dev.begin(), dev.end(), [&](size_t a, size_t b) { const int oa = W.indexOf(W.segs[a].outlet), ob = W.indexOf(W.segs[b].outlet); return (oa >= 0 ? x[size_t(oa)] : x[a]) < (ob >= 0 ? x[size_t(ob)] : x[b]); });
    auto series = [&](const std::vector<double>& v, const std::vector<size_t>& o, bool atOutlet, std::vector<double>& xo, std::vector<double>& yo) {
        xo.clear(); yo.clear();
        for (size_t p : o) if (p < v.size()) { const int ou = W.indexOf(W.segs[p].outlet); xo.push_back(atOutlet && ou >= 0 ? x[size_t(ou)] : x[p]); yo.push_back(v[p]); }
    };
    const float ph = std::max(fs * 12, (ImGui::GetContentRegionAvail().y - fs * 12) * 0.5f);
    // draw + record (tab export button) helpers
    auto L_ = [](const char* lab, const std::vector<double>& a, const std::vector<double>& b, int ci, float w) {
        if (fig::recorder().recording()) recSeries(lab, a.data(), b.data(), int(a.size()), seriesColor(ci), w);
        ImPlot::PlotLine(lab, a.data(), b.data(), int(a.size()), mswLine(ci, w));
    };
    auto S_ = [](const char* lab, const std::vector<double>& a, const std::vector<double>& b, int ci, bool stems) {
        if (fig::recorder().recording()) recSeries(lab, a.data(), b.data(), int(a.size()), seriesColor(ci), 1.5f, true, true);
        ImPlotSpec sp = mswLine(ci, stems ? 1.5f : 0.0f); sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3;
        if (stems) ImPlot::PlotStems(lab, a.data(), b.data(), int(a.size()), 0.0, sp); else ImPlot::PlotScatter(lab, a.data(), b.data(), int(a.size()), sp);
    };
    fig::recorder().beginFigure("run.msw", 2, 2);
    if (ImPlot::BeginSubplots("##mswprof", 2, 2, ImVec2(-1, 2 * ph))) {
        std::vector<double> xo, yo;
        if (ImPlot::BeginPlot("Node pressure")) {
            recPanel("Node pressure", xa[xAxis], pu);
            ImPlot::SetupAxes(xa[xAxis], pu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            series(Rr.p, ord, false, xo, yo); L_("tubing", xo, yo, 0, 2.0f);
            if (!dev.empty()) { series(Rr.p, dev, true, xo, yo); S_("device / branch nodes", xo, yo, 1, false); }
            if (showFirst && idx.size() > 1) { const auto& R0 = (*R.mswRec)[size_t(idx[0])]; series(R0.p, ord, false, xo, yo); L_("tubing (first record)", xo, yo, 2, 1.2f); }
            recEnd();
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Segment flow (+ toward the heel)")) {
            recPanel("Segment flow", xa[xAxis], qu);
            ImPlot::SetupAxes(xa[xAxis], qu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            series(Rr.q, ord, false, xo, yo); L_("tubing", xo, yo, 0, 2.0f);
            if (!dev.empty()) { series(Rr.q, dev, true, xo, yo); S_("inflow per device", xo, yo, 1, true); }
            recEnd();
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Friction / device / acceleration drop per segment")) {
            recPanel("Friction / device drop", xa[xAxis], pu, W.driftFlux ? "gas fraction" : "");
            ImPlot::SetupAxes(xa[xAxis], pu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            series(Rr.dFric, ord, false, xo, yo); L_("tubing friction", xo, yo, 3, 2.0f);
            if (!dev.empty()) { series(Rr.dDev, dev, true, xo, yo); S_("device", xo, yo, 4, true); }
            if (W.accelOn && Rr.dAcc.size() == ns) { series(Rr.dAcc, ord, false, xo, yo); L_("acceleration", xo, yo, 5, 1.5f); }
            bool anyPump = false;
            for (double d : Rr.dPump) anyPump = anyPump || d != 0;
            if (anyPump && Rr.dPump.size() == ns) { series(Rr.dPump, ord, false, xo, yo); S_("pump boost (WSEGPULL)", xo, yo, 8, true); }
            if (W.driftFlux && Rr.holdup.size() == ns) {   // C36: in-situ gas holdup vs no-slip fraction on the right axis
                ImPlot::SetupAxis(ImAxis_Y2, "gas fraction", ImPlotAxisFlags_AuxDefault | ImPlotAxisFlags_AutoFit);
                useY2(true);
                series(Rr.holdup, ord, false, xo, yo); L_("gas holdup (in situ)", xo, yo, 6, 1.5f);
                series(Rr.noSlip, ord, false, xo, yo); L_("no-slip gas fraction", xo, yo, 7, 1.2f);
                if (W.df.owFactor > 0 && Rr.oilHold.size() == ns) { series(Rr.oilHold, ord, false, xo, yo); L_("oil fraction of the liquid (in situ)", xo, yo, 9, 1.5f); }
                useY2(false);
            }
            recEnd();
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Heel / toe pressure history")) {
            recPanel("Heel / toe pressure", xl, pu, "toe - heel");
            ImPlot::SetupAxes(xl, pu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "toe - heel", ImPlotAxisFlags_AuxDefault | ImPlotAxisFlags_AutoFit);
            std::vector<double> t, heel, toe, d;
            const size_t ph0 = ord.empty() ? 0 : ord.front(), pt = ord.empty() ? 0 : ord.back();
            for (int k : idx) { const auto& r = (*R.mswRec)[size_t(k)]; if (r.p.size() < ns) continue; t.push_back(r.t * xs); heel.push_back(r.p[ph0]); toe.push_back(r.p[pt]); d.push_back(r.p[pt] - r.p[ph0]); }
            L_("heel (top node)", t, heel, 0, 2.0f);
            L_("toe", t, toe, 1, 2.0f);
            useY2(true);
            L_("toe - heel", t, d, 2, 1.5f);
            useY2(false);
            const double tc = Rr.t * xs;
            ImPlotSpec v; v.LineColor = ImVec4(1, 1, 1, 0.5f); ImPlot::PlotInfLines("##cur", &tc, 1, v);
            recVLine(tc);
            recEnd();
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
    fig::recorder().endFigure();
    if (ui::beginTable("##mswtab", 17, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable, ImVec2(0, 0))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        for (const char* h : {"seg", "branch", "outlet", "length", "depth", "ID", "device", "p", "q (local)", "d friction", "d device", "d accel", "gas holdup",
                              "oil / liquid", "head density", "pump boost", "limiter factor"}) ui::tableColumn(h);
        ImGui::TableHeadersRow();
        static const char* dn[] = {"-", "VALV", "SICD", "AICD", "LABY"};
        for (size_t p = 0; p < ns; ++p) {
            const auto& s = W.segs[p];
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("%d", s.id);
            ImGui::TableNextColumn(); ImGui::Text("%d", s.branch);
            ImGui::TableNextColumn(); ImGui::Text("%d", s.outlet);
            ImGui::TableNextColumn(); ImGui::Text("%.1f %s", L[p], lu);
            ImGui::TableNextColumn(); ImGui::Text("%.1f", s.depth);
            ImGui::TableNextColumn(); ImGui::Text("%.4g", s.diam);
            ImGui::TableNextColumn(); ImGui::Text("%s%s%s%s", dn[std::clamp(s.device, 0, 4)], (s.flimPh[0] || s.flimPh[1]) ? " +FLIM" : "", s.pump ? " +PUMP" : "", s.shut ? " SHUT" : "");
            ImGui::TableNextColumn(); ImGui::Text("%.2f", p < Rr.p.size() ? Rr.p[p] : 0.0);
            ImGui::TableNextColumn(); ImGui::Text("%.2f", p < Rr.q.size() ? Rr.q[p] : 0.0);
            ImGui::TableNextColumn(); ImGui::Text("%.4g", p < Rr.dFric.size() ? Rr.dFric[p] : 0.0);
            ImGui::TableNextColumn(); ImGui::Text("%.4g", p < Rr.dDev.size() ? Rr.dDev[p] : 0.0);
            ImGui::TableNextColumn(); ImGui::Text("%.4g", p < Rr.dAcc.size() ? Rr.dAcc[p] : 0.0);
            ImGui::TableNextColumn();
            if (p < Rr.holdup.size() && p < Rr.noSlip.size() && Rr.noSlip[p] > 0) ImGui::Text("%.3f (no-slip %.3f)", Rr.holdup[p], Rr.noSlip[p]); else ImGui::TextDisabled("-");
            ImGui::TableNextColumn(); if (p < Rr.oilHold.size() && Rr.oilHold[p] > 0) ImGui::Text("%.3f", Rr.oilHold[p]); else ImGui::TextDisabled("-");
            ImGui::TableNextColumn(); if (p < Rr.rhoHead.size()) ImGui::Text("%.4g", Rr.rhoHead[p]); else ImGui::TextDisabled("-");
            ImGui::TableNextColumn(); if (p < Rr.dPump.size() && Rr.dPump[p] != 0) ImGui::Text("%.4g", -Rr.dPump[p]); else ImGui::TextDisabled("-");
            ImGui::TableNextColumn(); if (p < Rr.flimFac.size() && Rr.flimFac[p] < 1) ImGui::Text("%.4f", Rr.flimFac[p]); else ImGui::TextDisabled("-");
        }
        ui::endTable();
    }
    ImGui::EndChild();
}

// ============================================================================ Model: WELSEGS generator
void App::mswGenerator(int wellIndex) {
    if (wellIndex < 0 || wellIndex >= int(case_.wells.size())) return;
    WellInput& w = case_.wells[size_t(wellIndex)];
    const bool field = case_.units == UnitSystem::Field;
    static double diam = -1, rough = -1, devStrength = 0.0004, rhoCal = 62.4, muCal = 1.45, aicdX = 3.0, aicdY = 0.7;
    static int device = 0;   // 0 none, 1 SICD, 2 AICD, 3 SHUT-able valve, 4 labyrinth (C36)
    static int pcomp = 1, flowModel = 0;   // C36: 0 HFA, 1 HF-, 2 H--; 0 HO, 1 DF
    static int labN = 1, labCh = 2;
    static double labA = -1, labD = -1;
    if (diam < 0) { diam = field ? 0.2034 : 0.062; rough = field ? 1e-4 : 3e-5; rhoCal = field ? 62.4 : 1000; }
    if (labA < 0) { labA = field ? 6e-5 / (0.3048 * 0.3048) : 6e-5; labD = field ? 0.006316 / 0.3048 : 0.006316; }
    // C37: oil-water slip, flow-limiting valve at the heel, downhole pump
    static double owFac = 0.0, owA = 1.2, owB1 = 0.4, owB2 = 0.7, owN = 2.0;
    static int limPh = 0, pumpOn = 0;   // limPh 0 none, 1 OIL, 2 WAT, 3 GAS, 4 LIQ
    static double limQ = 5000, limA = 8, pumpDp = 300, pumpP = 1.0, pumpQo = 0, pumpG = 0;
    // existing definition of this well
    int have = 0;
    for (const auto& e : case_.extraKeywords) if ((e.name == "WELSEGS" || e.name == "COMPSEGS") && firstToken(e.body) == w.name) ++have;
    ImGui::TextDisabled(have ? "%s is a multisegment well (%d WELSEGS / COMPSEGS keyword(s) in the schedule)." : "%s is a standard well.", w.name.c_str(), have);
    if (ui::beginProps("##mswgen")) {
        ui::propD(field ? "tubing ID (ft)" : "tubing ID (m)", &diam, "%.4g");
        ui::propD("roughness", &rough, "%.3g");
        static const char* pc[] = {"HFA: hydrostatic + friction + acceleration", "HF-: hydrostatic + friction", "H--: hydrostatic only"};
        ui::propCombo("pressure drop components", &pcomp, pc, 3);
        static const char* fm[] = {"HO: homogeneous (no slip)", "DF: drift flux (gas slip, Shi et al. 2005)"};
        ui::propCombo("flow model", &flowModel, fm, 2);
        static const char* dv[] = {"none (tubing only)", "spiral ICD per connection (WSEGSICD)", "autonomous ICD per connection (WSEGAICD)", "valve per connection (WSEGVALV)",
                                   "labyrinth ICD per connection (WSEGLABY)"};
        ui::propCombo("inflow devices", &device, dv, 5);
        // C37 tubing physics (ResSim options, written as '-- RESSIM: MSWOPTS' comment lines - OPM / ECLIPSE always flash)
        bool fo = case_.mswOpt.flash, so = case_.mswOpt.storage;
        if (ui::propBool("tubing flash (PVT at segment pressure)", &fo)) { case_.mswOpt.flash = fo; markDirty(); }
        if (ui::propBool("segment storage (transient wellbore)", &so)) { case_.mswOpt.storage = so; markDirty(); }
        {   // C45 cross-flow, C49 converged connection heads (MSWOPTS tokens 3 and 4)
            bool xo = case_.mswOpt.xflow;
            if (ui::propBool("producer cross-flow in segments", &xo)) { case_.mswOpt.xflow = xo; markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("A connection whose cell pressure is below the tubing pressure takes the fluid of its segment (OPM MultisegmentWell)."));
            static const char* hm[] = {"default (converged)", "explicit (start of step)", "converged (re-solved per step)"};
            int hsel = case_.mswOpt.heads < 0 ? 0 : case_.mswOpt.heads == 0 ? 1 : 2;
            if (ui::propCombo("connection heads", &hsel, hm, 3)) { case_.mswOpt.heads = hsel == 0 ? -1 : hsel == 1 ? 0 : 1; markDirty(); }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Converged: the segment network is re-solved from the end-of-step state and the step repeated until the connection heads stop moving (bracketed per connection, at most 8 passes). A Picard coupling, not a full Jacobian."));
        }
        if (flowModel == 1) {
            ui::propD("oil-water slip factor (0 = off)", &owFac, "%.3g", "WSEGDFPA item 7: multiplier of the Shi et al. (2005) oil-water drift velocity");
            if (owFac > 0) { ui::propD("oil-water profile A", &owA, "%.3g"); ui::propD("B1", &owB1, "%.3g"); ui::propD("B2", &owB2, "%.3g"); ui::propD("exponent N", &owN, "%.3g"); }
        }
        static const char* lp[] = {"none", "OIL", "WAT", "GAS", "LIQ"};
        ui::propCombo("flow-limiting valve at the heel (WSEGFLIM)", &limPh, lp, 5);
        if (limPh > 0) { ui::propD(field ? "limit (stb/d | Mscf/d)" : "limit (sm3/d)", &limQ, "%.5g"); ui::propD("sharpness a", &limA, "%.3g", "factor (1 + (Q0/Qlim)^a)^(-1/a); large a = hard limit"); }
        static const char* pp[] = {"none", "pump at the heel segment (WSEGPULL)"};
        ui::propCombo("downhole pump", &pumpOn, pp, 2);
        if (pumpOn) {
            ui::propD(field ? "pump boost DP (psi)" : "pump boost DP (bar)", &pumpDp, "%.4g");
            ui::propD("scale P", &pumpP, "%.3g");
            ui::propD("oil-rate limit", &pumpQo, "%.5g", "boost falls by G per unit oil rate above this limit");
            ui::propD("gradient G", &pumpG, "%.4g");
        }
        if (device == 1 || device == 2) {
            ui::propD("strength", &devStrength, "%.4g", "device strength (pressure / (rate^2) at the calibration fluid)");
            ui::propD("calibration density", &rhoCal, "%.4g");
            ui::propD("calibration viscosity (cP)", &muCal, "%.4g");
            if (device == 2) { ui::propD("flow exponent x", &aicdX, "%.3g"); ui::propD("viscosity exponent y", &aicdY, "%.3g"); }
        }
        if (device == 3) ui::propD("valve Cv", &devStrength, "%.4g");
        if (device == 4) {
            ui::propI("labyrinth units in parallel", &labN, 1, 20);
            ui::propI("channels per unit", &labCh, 1, 6);
            ui::propD(field ? "channel area (ft2)" : "channel area (m2)", &labA, "%.4g");
            ui::propD(field ? "channel hydraulic diameter (ft)" : "channel hydraulic diameter (m)", &labD, "%.4g");
            labN = std::max(1, labN); labCh = std::max(1, labCh);
        }
        ui::endProps();
    }
    if (!preview_) { ui::hint("Build the grid preview first (the segment lengths come from the cell centres)."); return; }
    const GridInput& g = case_.grid;
    // connection cells in order: completions as listed, layers top -> bottom inside each
    struct C { int i, j, k; double x, y, z; };
    std::vector<C> cs;
    for (const auto& cp : w.completions) {
        if (cp.lgr >= 0 || !cp.open) continue;
        for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
            if (cp.i < 0 || cp.i >= g.nx || cp.j < 0 || cp.j >= g.ny || k < 0 || k >= g.nz) continue;
            const int q = g.index(cp.i, cp.j, k);
            if (size_t(q) >= preview_->center.size()) continue;
            const auto& v = preview_->center[size_t(q)];
            cs.push_back({cp.i, cp.j, k, double(v.x), double(v.y), double(v.z)});
        }
    }
    ImGui::Text("%zu connection(s) -> %zu tubing segment(s)%s", cs.size(), cs.size(), device ? " + one device segment each" : "");
    if (cs.empty()) return;
    if (ui::accentButton(have ? ICON_REFRESH_CW " Regenerate the segment definition" : ICON_PLUS " Make it a multisegment well")) {
        // remove the old records of this well
        for (auto it = case_.extraKeywords.begin(); it != case_.extraKeywords.end();) {
            if ((it->name == "WELSEGS" || it->name == "COMPSEGS") && firstToken(it->body) == w.name) { it = case_.extraKeywords.erase(it); continue; }
            if (it->name == "WSEGSICD" || it->name == "WSEGAICD" || it->name == "WSEGVALV" || it->name == "WSEGLABY" || it->name == "WSEGFLIM" || it->name == "WSEGPULL") {
                std::istringstream is(it->body); std::string line, keep; int recs = 0;
                while (std::getline(is, line)) { if (firstToken(line) == w.name) continue; if (!line.empty() && line.find_first_not_of(" \t/") != std::string::npos) ++recs; keep += line + "\n"; }
                if (!recs) { it = case_.extraKeywords.erase(it); continue; }
                it->body = keep;
            }
            ++it;
        }
        // tubing nodes: top node half a spacing before the first connection, one node per connection
        const size_t n = cs.size();
        auto dist = [&](size_t a, size_t b) { return std::sqrt((cs[a].x - cs[b].x) * (cs[a].x - cs[b].x) + (cs[a].y - cs[b].y) * (cs[a].y - cs[b].y) + (cs[a].z - cs[b].z) * (cs[a].z - cs[b].z)); };
        const double first = n > 1 ? std::max(1e-3, 0.5 * dist(0, 1)) : std::max(1e-3, 0.5 * preview_->cellDz[size_t(g.index(cs[0].i, cs[0].j, cs[0].k))]);
        const double dz0 = n > 1 ? 0.5 * (cs[1].z - cs[0].z) : 0.5 * preview_->cellDz[size_t(g.index(cs[0].i, cs[0].j, cs[0].k))];
        const double topDepth = cs[0].z - dz0;
        char b[256];
        std::string ws;
        std::snprintf(b, sizeof b, " '%s' %.6g 0 1* INC '%s' '%s' /\n", w.name.c_str(), topDepth, pcomp == 0 ? "HFA" : pcomp == 1 ? "HF-" : "H--", flowModel ? "DF" : "HO"); ws += b;
        std::vector<double> segLen(n);
        for (size_t c = 0; c < n; ++c) {
            const double len = c == 0 ? first : std::max(1e-3, dist(c - 1, c));
            const double dz = c == 0 ? dz0 : cs[c].z - cs[c - 1].z;
            segLen[c] = len;
            std::snprintf(b, sizeof b, " %d %d 1 %d %.6g %.6g %.6g %.4g /\n", int(c) + 2, int(c) + 2, int(c) + 1, len, dz, diam, rough); ws += b;
        }
        std::string cg = " '" + w.name + "' /\n";
        if (device) {
            for (size_t c = 0; c < n; ++c) {
                std::snprintf(b, sizeof b, " %d %d %d %d 0.1 0 %.6g %.4g /\n", int(n + 2 + c), int(n + 2 + c), int(c) + 2, int(c) + 2, diam, rough); ws += b;
                std::snprintf(b, sizeof b, " %d %d %d %d 0.0 0.1 /\n", cs[c].i + 1, cs[c].j + 1, cs[c].k + 1, int(c) + 2); cg += b;
            }
        } else {
            double s0 = 0;
            for (size_t c = 0; c < n; ++c) {
                const double a = s0 + segLen[c] - 0.5 * (c + 1 < n ? segLen[c + 1] : segLen[c]), e = s0 + segLen[c] + 0.5 * (c + 1 < n ? segLen[c + 1] : segLen[c]);
                std::snprintf(b, sizeof b, " %d %d %d 1 %.6g %.6g /\n", cs[c].i + 1, cs[c].j + 1, cs[c].k + 1, std::max(0.0, a), e); cg += b;
                s0 += segLen[c];
            }
        }
        int step = 0;
        for (const auto& cp : w.completions) step = std::max(step, cp.fromStep);
        step = 0;   // the definition is in force from the start (connections opened later keep their windows)
        auto add = [&](const char* name, const std::string& body) {
            ExtraKeyword e; e.name = name; e.section = DeckSection::Schedule; e.step = step; e.body = body + "/\n"; e.note = "C32 Model -> Wells: multisegment well generator";
            case_.extraKeywords.push_back(e);
        };
        add("WELSEGS", ws);
        add("COMPSEGS", cg);
        if (device == 1) { std::snprintf(b, sizeof b, " '%s' %d %d %.6g 1* %.6g %.6g /\n", w.name.c_str(), int(n + 2), int(2 * n + 1), devStrength, rhoCal, muCal); add("WSEGSICD", b); }
        if (device == 2) { std::snprintf(b, sizeof b, " '%s' %d %d %.6g 1* %.6g %.6g 5* %.4g %.4g /\n", w.name.c_str(), int(n + 2), int(2 * n + 1), devStrength, rhoCal, muCal, aicdX, aicdY); add("WSEGAICD", b); }
        if (device == 3) {
            std::string vb;
            for (size_t c = 0; c < n; ++c) { std::snprintf(b, sizeof b, " '%s' %d %.4g %.6g /\n", w.name.c_str(), int(n + 2 + c), devStrength, 0.25 * M_PI * diam * diam * 0.1); vb += b; }
            add("WSEGVALV", vb);
        }
        if (device == 4) { std::snprintf(b, sizeof b, " '%s' %d %d %d %d %.6g 2* %.6g /\n", w.name.c_str(), int(n + 2), int(2 * n + 1), labN, labCh, labA, labD); add("WSEGLABY", b); }
        // C37: limiter / pump on segment 2 (the heel tubing segment carries the whole well), oil-water slip parameters
        static const char* lpn[] = {"", "OIL", "WAT", "GAS", "LIQ"};
        if (limPh > 0 && limQ > 0) { std::snprintf(b, sizeof b, " '%s' 2 2 %s %.6g %.4g /\n", w.name.c_str(), lpn[std::clamp(limPh, 1, 4)], limQ, std::max(1.0, limA)); add("WSEGFLIM", b); }
        if (pumpOn && pumpDp > 0) { std::snprintf(b, sizeof b, " '%s' 2 2 %.4g %.6g %.6g %.6g /\n", w.name.c_str(), std::max(1e-3, pumpP), pumpDp, std::max(0.0, pumpQo), std::max(0.0, pumpG)); add("WSEGPULL", b); }
        if (flowModel == 1 && owFac > 0) {
            for (auto it = case_.extraKeywords.begin(); it != case_.extraKeywords.end();) { if (it->name == "WSEGDFPA") it = case_.extraKeywords.erase(it); else ++it; }
            ExtraKeyword e; e.name = "WSEGDFPA"; e.section = DeckSection::Schedule; e.step = step; e.note = "C37 oil-water slip";
            std::snprintf(b, sizeof b, " 1* 1* 1* 1* 1* 1* %.4g %.4g %.4g %.4g %.4g /\n", owFac, owA, owB1, owB2, owN); e.body = b;
            case_.extraKeywords.push_back(e);
        }
        markDirty();
        notify(w.name + ": WELSEGS / COMPSEGS written (" + std::to_string(n) + " tubing segments" + (device ? ", devices" : "") + ") - see Deck editor; the ResSim engine and OPM Flow / ECLIPSE use them");
    }
    if (have) {
        ImGui::SameLine();
        if (ui::dangerButton("Make it a standard well")) {
            for (auto it = case_.extraKeywords.begin(); it != case_.extraKeywords.end();) {
                if ((it->name == "WELSEGS" || it->name == "COMPSEGS") && firstToken(it->body) == w.name) it = case_.extraKeywords.erase(it);
                else ++it;
            }
            markDirty();
        }
    }
    ui::hint("INC mode: segment lengths / depth increments between the connection cell centres; device segments (branch = connection "
             "number + 1, length 0.1) sit between each connection and its tubing node, as in examples/MSW_HORIZONTAL.DATA");
}

// ============================================================================ Results: surface network (C35-4b)
void App::networkTab(const RunResults& R, double xs, const char* xl) {
    if (R.netNodes.empty() || !R.netRec || R.netRec->empty()) { ui::hint("No surface network in this run (NODEPROP / BRANPROP or GRUPNET)."); return; }
    const bool field = R.input.units == UnitSystem::Field;
    const char* pu = field ? "psi" : "bar";
    const auto& N = R.netNodes;
    static std::vector<char> show;
    if (show.size() != N.size()) show.assign(N.size(), 1);
    ui::hint("Node pressures per accepted time step (fixed / terminal nodes are flat). A node on an automatic choke sits above its "
             "uptree node by the choke drop. The table shows the network at the end of the run (rates are the branch flows into "
             "each node's parent, after NEFAC and gas consumption).");
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    for (size_t i = 0; i < N.size(); ++i) {
        const float w = ImGui::CalcTextSize(N[i].name.c_str()).x + ImGui::GetFrameHeight() + 3 * ImGui::GetStyle().ItemSpacing.x;
        if (i) { ImGui::SameLine(); if (ImGui::GetCursorScreenPos().x + w > right) ImGui::NewLine(); }
        bool b = show[i] != 0;
        ImGui::PushID(int(i));
        if (ImGui::Checkbox(N[i].name.c_str(), &b)) show[i] = b ? 1 : 0;
        ImGui::PopID();
    }
    const float ph = std::max(220.0f, ImGui::GetContentRegionAvail().y * 0.55f);
    fig::recorder().beginFigure("run.network", 1, 1);
    if (ImPlot::BeginPlot("Network node pressure", ImVec2(-1, ph))) {
        recPanel("Network node pressure", xl, pu);
        ImPlot::SetupAxes(xl, pu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        std::vector<double> t, p;
        for (size_t i = 0; i < N.size(); ++i) {
            if (!show[i]) continue;
            t.clear(); p.clear();
            for (const auto& r : *R.netRec) if (i < r.p.size() && r.p[i] > 0) { t.push_back(r.t * xs); p.push_back(r.p[i]); }
            if (t.empty()) continue;
            if (fig::recorder().recording()) recSeries(N[i].name.c_str(), t.data(), p.data(), int(t.size()), seriesColor(int(i)), 2.0f, false, false);
            ImPlot::PlotLine(N[i].name.c_str(), t.data(), p.data(), int(t.size()), mswLine(int(i)));
        }
        recEnd();
        ImPlot::EndPlot();
    }
    fig::recorder().endFigure();
    if (ui::beginTable("##nettab", 10, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable, ImVec2(0, 0))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        for (const char* h : {"node", "parent", "branch", "ALQ", "NEFAC", "p", "oil", "water", "gas", "choke drop"}) ui::tableColumn(h);
        ImGui::TableHeadersRow();
        for (const auto& n : N) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(n.name.c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(n.parent.empty() ? (n.fixed ? "(terminal)" : "-") : n.parent.c_str());
            ImGui::TableNextColumn();
            if (n.parent.empty()) ImGui::TextUnformatted(n.fixed ? "fixed p" : "-");
            else if (n.noLoss) ImGui::TextUnformatted(n.choke ? "9999 + auto choke" : "9999 (no loss)");
            else ImGui::Text("VFP %d%s%s", n.vfp, n.addLift ? " + lift gas" : "", n.choke ? " + auto choke" : "");
            ImGui::TableNextColumn(); if (n.alqEq) ImGui::TextUnformatted(n.alqEq == 1 ? "oil dens." : "gas dens."); else ImGui::Text("%.4g", n.alq);
            ImGui::TableNextColumn(); ImGui::Text("%.3g", n.nefac);
            ImGui::TableNextColumn(); ImGui::Text("%.3f %s", n.p, pu);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", n.qo);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", n.qw);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", n.qg);
            ImGui::TableNextColumn(); if (n.choke) ImGui::Text("%.3f %s%s", n.dpChoke, pu, n.chokeGroup.empty() ? "" : (" (" + n.chokeGroup + ")").c_str()); else ImGui::TextUnformatted("-");
        }
        ui::endTable();
    }
}
