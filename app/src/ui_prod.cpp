// C28 production analysis tab of the Results mode (CMG CoFlow / IPSM-style engineering tools on the run):
//   decline-curve analysis (Arps fit of a well / field rate history, forecast, EUR and remaining reserves)
//   nodal analysis (composite Vogel IPR from the run, VLP of a tubing on the case PVT, operating point) and
//   generation of a VFPPROD table from the same tubing model into the case being edited.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "plotkit.h"
#include "ressim/prodanalysis.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ressim;

namespace {
void pline(const char* label, const std::vector<double>& x, const std::vector<double>& y, ImU32 col, float w = 2.0f, bool scatter = false) {
    const int n = int(std::min(x.size(), y.size()));
    if (n <= 0) return;
    if (fig::recorder().recording()) recSeries(label, x.data(), y.data(), n, col, w, scatter);
    ImPlotSpec sp;
    sp.LineColor = col4(col); sp.LineWeight = w * g_plotLineScale;
    if (scatter) {
        sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3.0f; sp.MarkerFillColor = col4(col); sp.MarkerLineColor = col4(col);
        ImPlot::PlotScatter(label, x.data(), y.data(), n, sp);
    } else ImPlot::PlotLine(label, x.data(), y.data(), n, sp);
}
struct DcaState { int stream = 0, phase = 0; double tStart = 0, bMax = 2.0, qEcon = 1.0, years = 20; bool fixB = false; double b = 0.5;
                  bool logY = true; prod::ArpsFit fit; double cumNow = 0, tNow = 0; std::string msg; };
struct NodalState { int well = -1; prod::Ipr ipr; prod::Tubing tb; double thp = 0, wct = 0, gor = 0; bool init = false;
                    std::string note, vfpMsg; prod::OperatingPoint op; int vfpId = 0; };
DcaState g_dca;
NodalState g_nod;
}  // namespace

void App::prodAnalysisTab(const RunResults& R) {
    const auto& S = R.summary;
    const float fs = ImGui::GetFontSize();
    const Case& C = R.input;
    const bool field = C.units == UnitSystem::Field;
    const char* qU = field ? "stb/d" : "sm3/d";
    const char* pU = field ? "psi" : "bar";
    const char* lU = field ? "ft" : "m";
    ui::hint("Engineering checks on the run (CMG CoFlow / IPSM-style): Arps decline curves with forecast and EUR, and nodal "
             "analysis - composite Vogel IPR of a producer from the run, vertical lift performance (VLP) of a tubing on the case "
             "PVT (no-slip homogeneous multiphase flow + Haaland friction; no hold-up / flow-pattern map, so gassy wells at low rate "
             "read optimistic), the operating point and a VFPPROD table for THP control in the simulator.");
    if (S.size() < 3) { ImGui::TextDisabled("The run has too few report steps."); return; }
    std::vector<int> prods;
    for (int w = 0; w < int(R.wellNames.size()); ++w) if (w < int(R.wellTypes.size()) && R.wellTypes[size_t(w)] == WellType::Producer) prods.push_back(w);

    if (ui::beginTable("##prodsplit", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableNextRow();
        // ------------------------------------------------------------------ decline curves
        ImGui::TableSetColumnIndex(0);
        ImGui::PushID("dca");
        ui::header("Decline-curve analysis (Arps)");
        DcaState& D = g_dca;
        std::vector<std::string> streams{"Field"};
        for (int w : prods) streams.push_back("Well " + R.wellNames[size_t(w)]);
        if (D.stream >= int(streams.size())) D.stream = 0;
        bool refit = false;
        if (ui::beginProps("dcap", 0.5f)) {
            refit |= ui::propCombo("Stream", &D.stream, streams);
            static const char* kPh[] = {"Oil rate", "Gas rate", "Water rate", "Liquid rate"};
            refit |= ui::propCombo("Rate", &D.phase, kPh, 4);
            refit |= ui::propD("Fit from (days)", &D.tStart, "%.6g", "Start of the decline period; the fit also starts at the rate peak after it");
            refit |= ui::propD("Max b", &D.bMax, "%.3g", "Upper bound of the hyperbolic exponent (b > 1 only for transient / multi-layer flow)");
            refit |= ui::propBool("Fix b", &D.fixB);
            if (D.fixB) refit |= ui::propD("b", &D.b, "%.3g", "0 exponential, 1 harmonic");
            ui::propD("Economic limit", &D.qEcon, "%.6g", "Rate at which the forecast stops (rate units of the stream)");
            ui::propD("Forecast horizon (years)", &D.years, "%.4g");
            ui::propBool("Log rate axis", &D.logY);
            ui::endProps();
        }
        auto rateOf = [&](const SummaryPoint& p) -> double {
            if (D.stream == 0) return D.phase == 0 ? p.fopr : D.phase == 1 ? p.fgpr : D.phase == 2 ? p.fwpr : p.fopr + p.fwpr;
            const size_t w = size_t(prods[size_t(D.stream - 1)]);
            if (w >= p.wells.size()) return 0.0;
            const WellRates& r = p.wells[w];
            return D.phase == 0 ? r.oil : D.phase == 1 ? r.gas : D.phase == 2 ? r.water : r.oil + r.water;
        };
        std::vector<double> t, q, tAll, qAll;
        double cum = 0;
        for (size_t i = 1; i < S.size(); ++i) {
            const double qi = rateOf(S[i]);
            cum += qi * (S[i].time - S[i - 1].time);
            tAll.push_back(S[i].time); qAll.push_back(qi);
            if (S[i].time >= D.tStart && qi > 0) { t.push_back(S[i].time); q.push_back(qi); }
        }
        D.cumNow = cum; D.tNow = S.back().time;
        if (refit || ImGui::Button("Fit decline") || (!D.fit.ok && D.msg.empty())) {
            D.fit = prod::fitArps(t, q, std::max(0.0, D.bMax), D.fixB ? std::max(0.0, D.b) : -1.0);
            D.msg = D.fit.ok ? "" : D.fit.error;
        }
        if (!D.msg.empty()) ImGui::TextColored(col4(Pal::Warn), "%s", D.msg.c_str());
        double eur = 0, tEnd = 0;
        if (D.fit.ok) {
            const prod::ArpsFit& F = D.fit;
            const double tMax = D.years * 365.25;
            eur = prod::arpsEUR(F, D.qEcon, tMax, &tEnd);
            const double cumToNow = prod::arpsCum(F, D.tNow), cumBefore = cum;   // produced so far (simulated)
            double cumAtT0 = 0;
            for (size_t i = 1; i < S.size() && S[i].time <= F.t0 + 1e-9; ++i) cumAtT0 += rateOf(S[i]) * (S[i].time - S[i - 1].time);
            const double remaining = std::max(0.0, eur - cumToNow);
            const char* kind = F.b < 0.05 ? "exponential" : std::fabs(F.b - 1) < 0.05 ? "harmonic" : "hyperbolic";
            const double effDecl = F.qi > 0 ? 100.0 * (1.0 - prod::arpsRate(F, F.t0 + 365.25) / F.qi) : 0.0;   // first-year effective decline
            ImGui::Text("%s: qi %.5g, Di %.4g /yr (%.3g %% first-year effective), b %.3f", kind, F.qi, F.Di * 365.25, effDecl, F.b);
            const ImU32 qc = F.r2 > 0.95 ? Pal::Ok : F.r2 > 0.8 ? Pal::Warn : Pal::Err;
            ImGui::TextColored(col4(qc), "R2 %.4f on %d points from t0 = %.4g d", F.r2, F.points, F.t0);
            ImGui::Text("Produced so far %.5g; EUR (t0 + %.4g d to q = %.4g) %.5g + %.5g before t0", cumBefore, tEnd, D.qEcon, eur, cumAtT0);
            ImGui::Text("Remaining reserves after the run end: %.5g", remaining);
        }
        const float ph = std::max(fs * 14, ImGui::GetContentRegionAvail().y - fs * 0.5f);
        if (ImPlot::BeginPlot("Rate history and Arps forecast", ImVec2(-1, ph))) {
            recPanel("Rate history and Arps forecast", "time (days)", "rate", "");
            ImPlot::SetupAxes("time (days)", "rate", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (D.logY) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            std::vector<double> qa;
            for (double v : qAll) qa.push_back(D.logY ? std::max(v, 1e-6) : v);
            pline("simulated", tAll, qa, Pal::Oil, 1.5f, true);
            if (D.fit.ok) {
                std::vector<double> tf, qf, tc, qc;
                const double t1 = D.fit.t0 + std::max(tEnd, D.tNow - D.fit.t0);
                for (int i = 0; i <= 200; ++i) {
                    const double tt = D.fit.t0 + (t1 - D.fit.t0) * i / 200.0;
                    (tt <= D.tNow ? tf : tc).push_back(tt);
                    (tt <= D.tNow ? qf : qc).push_back(std::max(prod::arpsRate(D.fit, tt), 1e-6));
                }
                pline("Arps fit", tf, qf, Pal::Accent, 2.0f);
                pline("forecast", tc, qc, Pal::Warn, 2.0f);
                const double ec = std::max(D.qEcon, 1e-6);
                std::vector<double> ex{D.fit.t0, t1}, ey{ec, ec};
                pline("economic limit", ex, ey, Pal::Err, 1.0f);
            }
            recEnd();
            ImPlot::EndPlot();
        }
        ImGui::PopID();

        // ------------------------------------------------------------------ nodal analysis
        ImGui::TableSetColumnIndex(1);
        ImGui::PushID("nodal");
        ui::header("Nodal analysis (IPR / VLP)");
        NodalState& N = g_nod;
        if (prods.empty()) ImGui::TextDisabled("The run has no producer.");
        else {
            if (N.well < 0 || std::find(prods.begin(), prods.end(), N.well) == prods.end()) { N.well = prods[0]; N.init = false; }
            std::vector<std::string> pn;
            int sel = 0;
            for (size_t i = 0; i < prods.size(); ++i) { pn.push_back(R.wellNames[size_t(prods[i])]); if (prods[i] == N.well) sel = int(i); }
            auto fromRun = [&]() {
                N.ipr = prod::iprFromResults(C, S, N.well, &N.note);
                const WellRates* W = nullptr;
                for (size_t k = S.size(); k-- > 0;)
                    if (size_t(N.well) < S[k].wells.size() && S[k].wells[size_t(N.well)].oil > 0) { W = &S[k].wells[size_t(N.well)]; break; }
                if (W) {
                    const double liq = W->oil + W->water;
                    N.wct = liq > 0 ? W->water / liq : 0.0;
                    N.gor = W->oil > 0 ? W->gas / W->oil : 0.0;
                }
                double depth = 0;
                for (const auto& ws : C.wells) if (ws.name == R.wellNames[size_t(N.well)]) depth = ws.refDepth;
                if (depth > 0 && !N.init) { N.tb.depthTVD = depth; N.tb.lengthMD = depth; }   // BHP reference depth of the well
                if (!N.init) {
                    N.tb.diameter = field ? 0.25 : 0.0762;           // 3 in
                    N.tb.roughness = field ? 5e-5 : 1.5e-5;
                    N.thp = field ? 200.0 : 14.0;
                    N.init = true;
                }
                N.op = prod::operatingPoint(C, N.ipr, N.tb, N.thp, N.wct, N.gor);
            };
            bool ch = false;
            if (ui::beginProps("ndp", 0.5f)) {
                if (ui::propCombo("Producer", &sel, pn)) { N.well = prods[size_t(sel)]; N.init = false; fromRun(); }
                ch |= ui::propD((std::string("Reservoir pressure (") + pU + ")").c_str(), &N.ipr.pr);
                ch |= ui::propD((std::string("Bubble point (") + pU + ")").c_str(), &N.ipr.pb, "%.6g", "0 = undersaturated straight-line IPR");
                ch |= ui::propD((std::string("PI J (") + qU + "/" + pU + ")").c_str(), &N.ipr.J, "%.6g", "Productivity index above the bubble point");
                ch |= ui::propD((std::string("Tubing TVD (") + lU + ")").c_str(), &N.tb.depthTVD);
                ch |= ui::propD((std::string("Tubing MD (") + lU + ")").c_str(), &N.tb.lengthMD);
                ch |= ui::propD((std::string("Tubing ID (") + lU + ")").c_str(), &N.tb.diameter, "%.5g");
                ch |= ui::propD((std::string("Roughness (") + lU + ")").c_str(), &N.tb.roughness, "%.3g");
                ch |= ui::propD((std::string("THP (") + pU + ")").c_str(), &N.thp);
                ch |= ui::propD("Water cut", &N.wct, "%.4g");
                ch |= ui::propD("Producing GOR", &N.gor, "%.5g", "Case Rs units");
                ui::endProps();
            }
            if (!N.init) fromRun();
            if (ImGui::Button("IPR + fluids from the run")) fromRun();
            N.tb.lengthMD = std::max(N.tb.lengthMD, N.tb.depthTVD);
            N.wct = std::clamp(N.wct, 0.0, 0.999); N.gor = std::max(0.0, N.gor);
            N.tb.diameter = std::max(N.tb.diameter, 1e-4);
            if (ch) N.op = prod::operatingPoint(C, N.ipr, N.tb, N.thp, N.wct, N.gor);
            if (!N.note.empty()) ImGui::TextDisabled("IPR %s", N.note.c_str());
            if (N.op.ok) ImGui::TextColored(col4(Pal::Ok), "Operating point: %.5g %s oil at pwf %.5g %s (AOF %.5g)", N.op.q, qU, N.op.pwf, pU, prod::iprAof(N.ipr));
            else ImGui::TextColored(col4(Pal::Warn), "No natural-flow operating point%s%s", N.op.note.empty() ? "" : ": ", N.op.note.c_str());
            // VFP generation
            ImGui::SetNextItemWidth(fs * 5);
            if (N.vfpId <= 0) { int mx = 0; for (const auto& v : case_.vfp) mx = std::max(mx, v.id); N.vfpId = mx + 1; }
            ImGui::InputInt("VFP id", &N.vfpId); N.vfpId = std::max(1, N.vfpId);
            ImGui::SameLine();
            if (ImGui::Button("Generate VFPPROD into the case")) {
                const double aof = std::max(prod::iprAof(N.ipr), 1.0);
                std::vector<double> rates, thps, wcts{0.0, 0.3, 0.6, 0.9}, gors;
                for (double f : {0.02, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5}) rates.push_back(f * aof);
                for (double f : {0.5, 1.0, 2.0, 4.0}) thps.push_back(std::max(1.0, f * std::max(N.thp, 1.0)));
                const double g0 = std::max(N.gor, field ? 0.3 : 50.0);
                for (double f : {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0}) gors.push_back(f * g0);   // dense GOR axis: BHP(GOR) is strongly non-linear
                VfpTable T = prod::makeVfpProd(C, N.tb, N.vfpId, rates, thps, wcts, gors);
                bool replaced = false;
                for (auto& v : case_.vfp) if (v.id == T.id && !v.inj) { v = T; replaced = true; }
                if (!replaced) case_.vfp.push_back(T);
                markDirty();
                char b[200];
                std::snprintf(b, sizeof b, "VFPPROD %d %s: %zu rates x %zu THP x %zu WCT x %zu GOR, datum %.5g %s. Use it with WCONPROD THP control (VFP table %d).",
                              T.id, replaced ? "replaced" : "added", rates.size(), thps.size(), wcts.size(), gors.size(), T.datum, lU, T.id);
                N.vfpMsg = b;
            }
            if (!N.vfpMsg.empty()) ImGui::TextDisabled("%s", N.vfpMsg.c_str());
            const float ph2 = std::max(fs * 14, ImGui::GetContentRegionAvail().y - fs * 0.5f);
            if (ImPlot::BeginPlot("IPR and VLP", ImVec2(-1, ph2))) {
                recPanel("IPR and VLP", qU, pU, "");
                ImPlot::SetupAxes((std::string("oil rate (") + qU + ")").c_str(), (std::string("bottom-hole pressure (") + pU + ")").c_str(),
                                  ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                const double aof = prod::iprAof(N.ipr);
                if (aof > 0) {
                    std::vector<double> qi, pi;
                    for (int i = 0; i <= 60; ++i) { const double pw = N.ipr.pr * (1.0 - i / 60.0); pi.push_back(pw); qi.push_back(prod::iprRate(N.ipr, pw)); }
                    pline("IPR", qi, pi, Pal::Oil, 2.0f);
                    const double fTh[3] = {0.5, 1.0, 2.0};
                    const ImU32 cTh[3] = {Pal::Muted, Pal::Water, Pal::Gas};
                    for (int k = 0; k < 3; ++k) {
                        std::vector<double> qv, pv;
                        for (int i = 1; i <= 30; ++i) {
                            const double qq = 1.2 * aof * i / 30.0;
                            const double pb = prod::vlpBhp(C, N.tb, qq, N.thp * fTh[k], N.wct, N.gor);
                            if (std::isfinite(pb)) { qv.push_back(qq); pv.push_back(pb); }
                        }
                        char lab[64]; std::snprintf(lab, sizeof lab, "VLP THP %.4g", N.thp * fTh[k]);
                        pline(lab, qv, pv, cTh[k], k == 1 ? 2.0f : 1.2f);
                    }
                    if (N.op.ok) {
                        std::vector<double> ox{N.op.q}, oy{N.op.pwf};
                        pline("operating point", ox, oy, Pal::Err, 3.0f, true);
                    }
                }
                recEnd();
                ImPlot::EndPlot();
            }
        }
        ImGui::PopID();
        ui::endTable();
    }
}
