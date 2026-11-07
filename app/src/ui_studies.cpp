// Studies window: uncertain parameters, observed data, sensitivity analysis (tornado, Monte Carlo) and
// history matching (manual sliders with live misfit, Levenberg-Marquardt, ES-MDA). All runs execute on a
// background thread through the parallel batch runner of the engine.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include "icons.h"
#include "theme.h"
#include <thread>
#include <chrono>
#include <cfloat>
#include "ressim/io.h"
#include "ressim/ecl.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <numeric>

using namespace ressim;

namespace {
const char* kParamKinds[] = {"Perm X/Y multiplier", "Perm Z multiplier", "Porosity multiplier", "Fault transmissibility (MULTFLT)",
                             "Aquifer strength multiplier", "krw end-point multiplier", "krw curvature exponent",
                             "krow curvature exponent", "Water viscosity multiplier", "Rock compressibility multiplier",
                             "OWC depth", "Datum pressure", "Well skin", "Well PI multiplier", "Oil viscosity multiplier", "GOC depth",
                             "Net-to-gross multiplier (box)", "Vertical transmissibility MULTZ (box)", "krg curvature exponent",
                             "krog curvature exponent", "krow end-point multiplier", "Capillary pressure multiplier",
                             "Injection rate multiplier", "Production rate multiplier",
                             "Geomodel variogram range multiplier", "Geomodel realisation (seed)"};
const int kNumKinds = int(sizeof(kParamKinds) / sizeof(kParamKinds[0]));
static_assert(sizeof(kParamKinds) / sizeof(kParamKinds[0]) == size_t(ressim::kNumParamKinds), "parameter kind names");
bool isGeoKind(ParamKind k) { return k == ParamKind::GeoRange || k == ParamKind::GeoSeed; }
const char* kResp[] = {"FOPT - cumulative oil", "Oil recovery factor", "FWPT - cumulative water", "FGPT - cumulative gas",
                       "FWIT - cumulative water injected", "FPR - final average pressure", "FWCT - final water cut",
                       "NPV (economics)", "History-match misfit",
                       // C27: coupled mechanics (need FAULTMECH / GEOMECH in the case - studies now run the coupled physics)
                       "Expected induced events (faults)", "Largest magnitude Mw (faults)", "Max Coulomb stress change (faults)",
                       "Surface uplift [m] (geomechanics)", "Max failure index (geomechanics)"};
static_assert(sizeof(kResp) / sizeof(kResp[0]) == size_t(ressim::kNumResponses), "response names");
const char* kObsNames[] = {"WOPR", "WWPR", "WGPR", "WWIR", "WBHP", "WWCT", "WGOR", "FOPR", "FWPR", "FGPR", "FPR", "FWCT", "FGOR"};
const int kNumObs = 13;

bool needsTarget(ParamKind k) {
    return k == ParamKind::FaultMult || k == ParamKind::AquiferStrength || k == ParamKind::WellSkin || k == ParamKind::WellPI ||
           k == ParamKind::InjRate || k == ParamKind::ProdRate || isGeoKind(k);
}

void setDefaults(StudyParam& p, const Case& c) {
    p.label.clear();
    p.i1 = p.i2 = p.j1 = p.j2 = p.k1 = p.k2 = -1;
    p.base = 1; p.low = 0.5; p.high = 2; p.log = true;
    const bool field = c.units == UnitSystem::Field;
    switch (p.kind) {
        case ParamKind::KrwMax: p.low = 0.5; p.high = 1.5; break;
        case ParamKind::KrwExp: case ParamKind::KroExp: p.low = 0.6; p.high = 1.6; break;
        case ParamKind::RockComp: p.low = 0.5; p.high = 3; break;
        case ParamKind::FaultMult:
            p.base = (p.target >= 0 && p.target < int(c.grid.faults.size())) ? std::max(1e-3, c.grid.faults[p.target].mult) : 1.0;
            p.low = 1e-3; p.high = 1; break;
        case ParamKind::AquiferStrength: p.low = 0.2; p.high = 5; break;
        case ParamKind::OWC: p.base = c.equil.owc; p.low = c.equil.owc - (field ? 50 : 15); p.high = c.equil.owc + (field ? 50 : 15); p.log = false; break;
        case ParamKind::DatumPressure:
            p.base = c.equil.datumPressure; p.low = 0.9 * p.base; p.high = 1.1 * p.base; p.log = false; break;
        case ParamKind::WellSkin: p.base = 0; p.low = -3; p.high = 10; p.log = false; break;
        case ParamKind::GOC: p.base = c.equil.goc; p.low = c.equil.goc - (field ? 50 : 15); p.high = c.equil.goc + (field ? 50 : 15); p.log = false; break;
        case ParamKind::KrgExp: case ParamKind::KrogExp: p.low = 0.6; p.high = 1.6; break;
        case ParamKind::KrowMax: p.low = 0.6; p.high = 1.2; break;
        case ParamKind::Ntg: p.low = 0.6; p.high = 1.2; break;
        case ParamKind::MultZ: p.low = 0.01; p.high = 1; break;
        case ParamKind::PcMult: p.low = 0.2; p.high = 3; break;
        case ParamKind::InjRate: case ParamKind::ProdRate: p.target = -1; p.low = 0.7; p.high = 1.3; p.log = false; break;
        case ParamKind::GeoRange: p.target = -1; p.low = 0.5; p.high = 2; p.log = true; break;
        case ParamKind::GeoSeed:
            p.target = -1; p.base = c.propModels.empty() ? 1.0 : double(c.propModels[0].seed); p.low = 1; p.high = 100; p.log = false; break;
        default: break;
    }
}

std::string fmtv(double v) { char b[32]; std::snprintf(b, sizeof(b), "%.4g", v); return b; }

struct ObsGroup { ObsQty qty; std::string well; std::vector<int> idx; };
std::vector<ObsGroup> obsGroups(const Case& c) {
    std::vector<ObsGroup> g;
    for (int i = 0; i < int(c.observations.size()); ++i) {
        const Observation& o = c.observations[i];
        auto it = std::find_if(g.begin(), g.end(), [&](const ObsGroup& x) { return x.qty == o.qty && x.well == o.well; });
        if (it == g.end()) { g.push_back({o.qty, o.well, {}}); it = g.end() - 1; }
        it->idx.push_back(i);
    }
    return g;
}
std::string groupName(const ObsGroup& g) { return std::string(toString(g.qty)) + (g.well.empty() ? "" : " " + g.well); }

ImU32 seriesCol(int i) {
    static const ImU32 c[] = {IM_COL32(255, 179, 64, 255), IM_COL32(186, 104, 200, 255), IM_COL32(77, 208, 225, 255),
                              IM_COL32(255, 138, 101, 255), IM_COL32(174, 213, 129, 255), IM_COL32(240, 98, 146, 255)};
    return c[i % 6];
}

// horizontal bar chart around a base value (tornado / Spearman) drawn with the ImGui draw list
void barChart(const char* id, const std::vector<std::string>& labels, const std::vector<double>& lo, const std::vector<double>& hi, double base,
              const char* loName, const char* hiName, bool centred) {
    const int n = int(labels.size());
    if (n == 0) return;
    const float fs = ImGui::GetFontSize();
    const float rowH = fs * 1.7f;
    float labW = 0;
    for (auto& l : labels) labW = std::max(labW, ImGui::CalcTextSize(l.c_str()).x);
    labW += fs;
    const float W = ImGui::GetContentRegionAvail().x;
    const float H = rowH * n + fs * 2.6f;
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton(id, ImVec2(W, H));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    double vmin = base, vmax = base;
    for (int i = 0; i < n; ++i) {
        if (std::isfinite(lo[i])) { vmin = std::min(vmin, lo[i]); vmax = std::max(vmax, lo[i]); }
        if (std::isfinite(hi[i])) { vmin = std::min(vmin, hi[i]); vmax = std::max(vmax, hi[i]); }
    }
    if (centred) { const double m = std::max(std::fabs(vmin - base), std::fabs(vmax - base)); vmin = base - m; vmax = base + m; }
    if (vmax - vmin < 1e-12 * std::max(1.0, std::fabs(base))) { vmin -= 1; vmax += 1; }
    const double pad = 0.08 * (vmax - vmin);
    vmin -= pad; vmax += pad;
    // reserve room for the value labels on both sides so they never sit on top of a bar
    float valW = 0;
    for (int i = 0; i < n; ++i) {
        if (std::isfinite(lo[i])) valW = std::max(valW, ImGui::CalcTextSize(fmtv(lo[i]).c_str()).x);
        if (std::isfinite(hi[i])) valW = std::max(valW, ImGui::CalcTextSize(fmtv(hi[i]).c_str()).x);
    }
    valW = std::min(valW + 6.0f, (W - labW) * 0.3f);
    const float x0 = p0.x + labW + valW, x1 = p0.x + W - fs * 0.5f - valW;
    auto X = [&](double v) { return float(x0 + (v - vmin) / (vmax - vmin) * (x1 - x0)); };
    const ImU32 cLo = IM_COL32(64, 132, 246, 230), cHi = IM_COL32(255, 159, 64, 230);
    for (int i = 0; i < n; ++i) {
        const float y = p0.y + i * rowH;
        if (i % 2 == 0) dl->AddRectFilled(ImVec2(p0.x, y), ImVec2(p0.x + W, y + rowH), IM_COL32(255, 255, 255, 8));
        dl->AddText(ImVec2(p0.x + 4, y + (rowH - fs) * 0.5f), IM_COL32(205, 210, 220, 255), labels[i].c_str());
        const float yb0 = y + rowH * 0.18f, yb1 = y + rowH * 0.82f;
        if (std::isfinite(lo[i])) dl->AddRectFilled(ImVec2(std::min(X(base), X(lo[i])), yb0), ImVec2(std::max(X(base), X(lo[i])) + 1, yb1), cLo, 3);
        if (std::isfinite(hi[i])) dl->AddRectFilled(ImVec2(std::min(X(base), X(hi[i])), yb0), ImVec2(std::max(X(base), X(hi[i])) + 1, yb1), cHi, 3);
        // value labels at the bar ends
        // a value left of the base gets its label left of the bar end, right of base -> right of it.
        // When both bars lie on the same side the nearer label falls inside the longer bar: draw it dark.
        auto lab = [&](double v, bool inside) {
            if (!std::isfinite(v)) return;
            const std::string t = fmtv(v);
            const float tw = ImGui::CalcTextSize(t.c_str()).x;
            const bool left = v < base;
            float xx = left ? X(v) - tw - 4 : X(v) + 4;
            if (inside) xx = left ? X(base) - tw - 4 : X(base) + 4;   // nearer value: inside the longer bar, next to the base line
            xx = std::clamp(xx, x0 - valW, x1 + valW - tw);
            dl->AddText(ImVec2(xx, y + (rowH - fs) * 0.5f), inside ? IM_COL32(20, 22, 28, 255) : IM_COL32(230, 232, 238, 220), t.c_str());
        };
        const bool fl = std::isfinite(lo[i]), fh = std::isfinite(hi[i]);
        if (fl && fh && ((lo[i] < base) == (hi[i] < base))) {
            const bool loFar = std::fabs(lo[i] - base) >= std::fabs(hi[i] - base);
            lab(loFar ? lo[i] : hi[i], false); lab(loFar ? hi[i] : lo[i], true);
        } else { if (fl) lab(lo[i], false); if (fh) lab(hi[i], false); }
    }
    const float yA = p0.y + n * rowH;
    dl->AddLine(ImVec2(X(base), p0.y), ImVec2(X(base), yA), IM_COL32(255, 255, 255, 200), 1.5f);
    dl->AddLine(ImVec2(x0, yA), ImVec2(x1, yA), IM_COL32(150, 156, 170, 160));
    const std::string bt = (centred ? "0" : "base " + fmtv(base));
    const float btw = ImGui::CalcTextSize(bt.c_str()).x;
    const float bx = std::clamp(X(base) - btw * 0.5f, x0, x1 - btw);
    dl->AddText(ImVec2(bx, yA + 3), IM_COL32(200, 204, 214, 255), bt.c_str());
    const std::string a = fmtv(vmin + pad), b = fmtv(vmax - pad);
    const float aw = ImGui::CalcTextSize(a.c_str()).x, bw = ImGui::CalcTextSize(b.c_str()).x;
    const float ax = X(vmin + pad) - aw * 0.5f, bxx = X(vmax - pad) - bw * 0.5f;   // centred under the extreme ticks
    if (ax + aw + fs * 0.5f < bx) dl->AddText(ImVec2(ax, yA + 3), IM_COL32(150, 156, 170, 255), a.c_str());
    if (bxx > bx + btw + fs * 0.5f) dl->AddText(ImVec2(bxx, yA + 3), IM_COL32(150, 156, 170, 255), b.c_str());
    for (double tv : {vmin + pad, vmax - pad}) dl->AddLine(ImVec2(X(tv), yA), ImVec2(X(tv), yA + 3), IM_COL32(150, 156, 170, 160));
    // legend
    const float ly = yA + fs * 1.35f;
    float lx = p0.x + 4;
    dl->AddRectFilled(ImVec2(lx, ly + 2), ImVec2(lx + fs, ly + fs - 2), cLo, 2); lx += fs * 1.3f;
    dl->AddText(ImVec2(lx, ly), IM_COL32(200, 204, 214, 255), loName); lx += ImGui::CalcTextSize(loName).x + fs;
    dl->AddRectFilled(ImVec2(lx, ly + 2), ImVec2(lx + fs, ly + fs - 2), cHi, 2); lx += fs * 1.3f;
    dl->AddText(ImVec2(lx, ly), IM_COL32(200, 204, 214, 255), hiName);
}

bool paramSlider(const char* id, const StudyParam& p, double* v) {
    float u = float(paramToUnit(p, *v));
    char fmt[48]; std::snprintf(fmt, sizeof(fmt), "%s", fmtv(*v).c_str());
    if (ImGui::SliderFloat(id, &u, 0.0f, 1.0f, fmt)) { *v = paramFromUnit(p, std::clamp(double(u), 0.0, 1.0)); return true; }
    return false;
}
} // namespace

// ============================================================================ background worker
const char* App::hmName(bool sh) const {
    if (study_.hmMethod == 1) return "ES-MDA";
    if (study_.hmAlgo == 2) return sh ? "Adjoint L-BFGS" : "Adjoint gradient + L-BFGS";
    if (study_.hmAlgo == 1) return sh ? "LM (analytic)" : "Levenberg-Marquardt (analytic Jacobian)";
    return sh ? "LM" : "Levenberg-Marquardt";
}

// misfit gradient at the base parameters by the discrete adjoint, forward sensitivities and central finite differences
// of full runs - verifies the analytic gradient on the user's own case before an adjoint / analytic-LM match
void App::gradientCheckPanel(int act) {
    ui::header("Gradient check");
    ImGui::TextDisabled("d(misfit) / d(normalised parameter) at the base values");
    bool ma = study_.gradMask & 1, mf = study_.gradMask & 2, mfd = study_.gradMask & 4;
    ImGui::Checkbox("Adjoint", &ma); ImGui::SameLine();
    ImGui::Checkbox("Forward sens.", &mf); ImGui::SameLine();
    ImGui::Checkbox("Finite diff.", &mfd);
    study_.gradMask = (ma ? 1 : 0) | (mf ? 2 : 0) | (mfd ? 4 : 0);
    if (mfd && ui::beginProps("gcopt", 0.45f)) {
        ui::propD("FD step (normalised)", &study_.gradFdStep, "%.3g", "Central differences of full runs: 2 runs per parameter. "
                  "Control / phase switching events make the misfit piecewise smooth, so FD over long schedules can differ "
                  "from the exact (branch) gradient by ~1e-3 relative.");
        ui::endProps();
    }
    ImGui::BeginDisabled(study_.running || act == 0 || study_.gradMask == 0 || case_.observations.empty());
    if (ImGui::Button(ICON_ACTIVITY " Check gradient", ImVec2(-FLT_MIN, 0))) {
        const Case base = case_; const int mask = study_.gradMask; const double h = study_.gradFdStep; const double rs = study_.adj.resStep;
        study_.deferred = [this, base, mask, h, rs]() {
            startStudy("Gradient check", [base, mask, h, rs](StudyState& S, const StudyCallbacks& cb) {
                GradientCheck g = misfitGradient(base, base.studyParams, mask, h, cb, rs);
                if (cb.log) {
                    char b[220];
                    if (g.ok) std::snprintf(b, sizeof(b), "gradient check: misfit %.6g, adjoint %.2fs, forward %.2fs, FD %.2fs (%d runs), replay deviation %.1e",
                                            g.misfit, g.secAdjoint, g.secForward, g.secFd, g.fdRuns, g.replayDev);
                    else std::snprintf(b, sizeof(b), "gradient check failed: %s", g.error.c_str());
                    cb.log(b);
                }
                std::lock_guard<std::mutex> l(S.m); S.gradCheck = std::move(g); S.gradParams = base.studyParams; S.hasGradCheck = true;
            });
        };
    }
    ImGui::EndDisabled();
    std::lock_guard<std::mutex> l(study_.m);
    if (!study_.hasGradCheck) return;
    const GradientCheck& G = study_.gradCheck;
    if (!G.ok) { ImGui::TextColored(col4(Pal::Err), "%s", G.error.c_str()); return; }
    ImGui::Text("misfit %.6g", G.misfit);
    if (!G.adjoint.empty()) { ImGui::SameLine(); ImGui::TextDisabled(" adjoint %.2fs", G.secAdjoint); }
    if (!G.forward.empty()) { ImGui::SameLine(); ImGui::TextDisabled(" forward %.2fs", G.secForward); }
    if (!G.fd.empty()) { ImGui::SameLine(); ImGui::TextDisabled(" FD %.2fs (%d runs)", G.secFd, G.fdRuns); }
    const std::vector<double>& ref = !G.adjoint.empty() ? G.adjoint : G.forward;
    double gmax = 0;
    for (double v : ref) gmax = std::max(gmax, std::fabs(v));
    for (double v : G.fd) gmax = std::max(gmax, std::fabs(v));
    const bool rel = !ref.empty() && !G.fd.empty();
    if (ui::beginTable("gck", 4 + (rel ? 1 : 0), ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ui::tableColumn("Parameter"); ui::tableColumn("Adjoint"); ui::tableColumn("Forward");
        ui::tableColumn("Finite diff."); if (rel) ui::tableColumn("|analytic - FD| / |g|max");
        ImGui::TableHeadersRow();
        for (size_t k = 0; k < G.act.size(); ++k) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            const int pi = G.act[k];
            ImGui::TextUnformatted(pi >= 0 && size_t(pi) < study_.gradParams.size() ? paramLabel(study_.gradParams[size_t(pi)], case_).c_str() : "?");
            auto cell = [&](const std::vector<double>& v) { ImGui::TableNextColumn(); if (k < v.size()) ImGui::Text("%+.6e", v[k]); else ImGui::TextDisabled("-"); };
            cell(G.adjoint); cell(G.forward); cell(G.fd);
            if (rel) {
                ImGui::TableNextColumn();
                const double e = gmax > 0 ? std::fabs(ref[k] - G.fd[k]) / gmax : 0;
                ImGui::TextColored(col4(e < 1e-4 ? Pal::Ok : e < 1e-2 ? Pal::Warn : Pal::Err), "%.2e", e);
            }
        }
        ui::endTable();
    }
    ImGui::TextDisabled("replay deviation %.1e, linear residual %.1e", G.replayDev, G.maxLinRes);
}

void App::startStudy(const std::string& task, std::function<void(App::StudyState&, const StudyCallbacks&)> fn) {
    if (study_.running) return;
    if (case_.comp.enabled) {   // the study drivers wrap the black-oil engine
        notify("Sensitivity / history matching runs the black-oil engine - not yet available for compositional or CO2STORE cases (see ROADMAP)");
        return;
    }
    if (study_.th.joinable()) study_.th.join();
    study_.cancel = false; study_.running = true; study_.done = 0; study_.total = 0;
    study_.t0 = ImGui::GetTime();
    // run archive context: every simulation of the study is stored with its misfit, objective terms and response
    auto ctx = std::make_shared<ArchiveCtx>();
    ctx->base = case_;
    ctx->sigma = autoSigma(case_.observations);
    ctx->terms = objectiveTerms(case_.observations, ObjectiveGrouping(study_.grouping));
    ctx->resp = Response(study_.response);
    ctx->eco = study_.eco;
    {
        std::lock_guard<std::mutex> l(study_.m);
        bool same = study_.archParams.size() == case_.studyParams.size();
        for (size_t i = 0; same && i < case_.studyParams.size(); ++i)
            same = study_.archParams[i].kind == case_.studyParams[i].kind && study_.archParams[i].target == case_.studyParams[i].target;
        if (!same) { study_.archive.clear(); study_.selRun = -1; }
        study_.archParams = case_.studyParams;
        study_.termNames.clear();
        for (const auto& t : ctx->terms) study_.termNames.push_back(t.name);
        study_.studyFirstId = study_.nextId;
        study_.task = task;
        study_.error.clear();
        study_.log.push_back("=== " + task);
    }
    study_.th = std::thread([this, fn, task, ctx]() {
        StudyCallbacks cb;
        cb.cancel = &study_.cancel;
        cb.workers = study_.workers;
        cb.log = [this](const std::string& s) { std::lock_guard<std::mutex> l(study_.m); study_.log.push_back(s); };
        cb.progress = [this](int d, int t) { study_.done = d; study_.total = t; };
        cb.onRun = [this, ctx](const char* tag, const std::vector<double>& values, const RunResult& run) { archiveRun(*ctx, tag, values, run); };
        const auto t0 = std::chrono::steady_clock::now();
        try { fn(study_, cb); }
        catch (const std::exception& e) {
            std::lock_guard<std::mutex> l(study_.m);
            study_.error = e.what();
            study_.log.push_back(std::string("ERROR: ") + e.what());
        }
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        {
            std::lock_guard<std::mutex> l(study_.m);
            study_.seconds = sec;
            char b[160]; std::snprintf(b, sizeof(b), "=== %s %s in %.1f s", task.c_str(), study_.cancel ? "cancelled" : "finished", sec);
            study_.log.push_back(b);
        }
        study_.running = false;
    });
}

void App::archiveRun(const ArchiveCtx& ctx, const char* tag, const std::vector<double>& values, const RunResult& run) {
    StudyState::ArchRun a;
    a.tag = tag ? tag : "";
    a.values = values;
    a.run = run;
    a.ok = run.ok;
    if (run.ok) {
        if (!ctx.base.observations.empty()) {
            a.misfit = misfit(run, ctx.base, ctx.base.observations, ctx.sigma);
            a.terms = misfitTerms(run, ctx.base, ctx.base.observations, ctx.sigma, ctx.terms);
        }
        try { a.response = responseValue(run, ctx.base, ctx.resp, ctx.eco, ctx.base.observations, ctx.sigma); } catch (...) {}
    }
    std::lock_guard<std::mutex> l(study_.m);
    a.id = study_.nextId++;
    if (study_.archive.size() >= 4000) study_.archive.erase(study_.archive.begin());   // bounded memory
    study_.archive.push_back(std::move(a));
    ++study_.archiveVersion;
}

void App::runManual() {
    if (study_.manualRunning) { study_.manualPending = true; return; }
    if (study_.manualTh.joinable()) study_.manualTh.join();
    study_.manualPending = false;
    study_.manualRunning = true;
    const Case base = case_;
    const std::vector<StudyParam> params = case_.studyParams;
    const std::vector<double> vals = study_.manualVals;
    study_.manualTh = std::thread([this, base, params, vals]() {
        Case c = applyParams(base, params, vals);
        RunResult r = runSingle(c, nullptr);
        const std::vector<double> sig = autoSigma(base.observations);
        const double m = r.ok ? misfit(r, base, base.observations, sig) : NAN;
        {
            ArchiveCtx ctx; ctx.base = base; ctx.sigma = sig; ctx.terms = objectiveTerms(base.observations, ObjectiveGrouping(study_.grouping));
            ctx.resp = Response(study_.response); ctx.eco = study_.eco;
            archiveRun(ctx, "manual", vals, r);
        }
        {
            std::lock_guard<std::mutex> l(study_.m);
            study_.manual = std::move(r);
            study_.manualMisfit = m;
            study_.hasManual = true;
            ++study_.manualVersion;
        }
        study_.manualRunning = false;
    });
}

// ============================================================================ window
void App::panelStudies() {
    if (!showStudies_) return;
    if (!ImGui::Begin("Studies", &showStudies_)) { ImGui::End(); return; }
    // keep manual values in sync with the parameter list
    if (study_.manualVals.size() != case_.studyParams.size()) {
        study_.manualVals.resize(case_.studyParams.size());
        for (size_t i = 0; i < case_.studyParams.size(); ++i) study_.manualVals[i] = case_.studyParams[i].base;
    }
    if (study_.manualPending && !study_.manualRunning && ImGui::GetTime() - study_.manualEditTime > 0.25) runManual();
    if (ui::accentButton(ICON_FLASK_CONICAL "  Open the Studies workspace (Ctrl+4)")) workspace_ = Workspace::Studies;
    ImGui::SameLine();
    ImGui::TextDisabled("full screen: run archive, profiles overlay, cross-plots, Pareto, proxy models");

    // ---- header: progress of the running study
    if (study_.running) {
        const int d = study_.done, t = study_.total;
        char b[128];
        std::string task; { std::lock_guard<std::mutex> l(study_.m); task = study_.task; }
        std::snprintf(b, sizeof(b), "%s  -  %d / %d runs  (%.0f s)", task.c_str(), d, t, ImGui::GetTime() - study_.t0);
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col4(Pal::Accent));
        ImGui::ProgressBar(t > 0 ? float(d) / t : 0.0f, ImVec2(ImGui::GetContentRegionAvail().x - ImGui::GetFontSize() * 6, 0), b);
        ImGui::PopStyleColor();
        ImGui::SameLine();
        if (ui::dangerButton("Cancel", ImVec2(-FLT_MIN, 0))) study_.cancel = true;
    }

    if (ImGui::BeginTabBar("studytabs")) {
        std::unique_lock<std::mutex> lk(study_.m);
        if (ui::tabItem("Parameters")) { studyParams(); ImGui::EndTabItem(); }
        if (ui::tabItem("Observed data")) { studyObservations(); ImGui::EndTabItem(); }
        if (ui::tabItem("Sensitivity")) { studySensitivity(); ImGui::EndTabItem(); }
        if (ui::tabItem("History matching")) { studyHistoryMatch(); ImGui::EndTabItem(); }
        if (ui::tabItem("Study log")) {
            if (ImGui::SmallButton("Clear")) study_.log.clear();
            ImGui::BeginChild("slog", ImVec2(0, 0), ImGuiChildFlags_Borders);
            if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
            ImGuiListClipper clip; clip.Begin(int(study_.log.size()));
            while (clip.Step())
                for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                    const std::string& s = study_.log[i];
                    if (s.rfind("ERROR", 0) == 0) ImGui::TextColored(col4(Pal::Err), "%s", s.c_str());
                    else if (s.rfind("===", 0) == 0) ImGui::TextColored(col4(Pal::Accent), "%s", s.c_str());
                    else ImGui::TextUnformatted(s.c_str());
                }
            if (g_fonts.mono) ImGui::PopFont();
            if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4) ImGui::SetScrollHereY(1.0f);
            ImGui::EndChild();
            ImGui::EndTabItem();
        }
        lk.unlock();
        ImGui::EndTabBar();
    }
    ImGui::End();
    if (study_.deferred) { auto f = std::move(study_.deferred); study_.deferred = nullptr; f(); }
}

// ============================================================================ parameters
void App::studyParams() {
    auto& P = case_.studyParams;
    const Case& c = case_;
    ui::hint("Uncertain parameters used by the sensitivity and history-matching tools. Multipliers act on the current case "
             "(1 = unchanged); a box limits a property multiplier to a sub-volume (1-based, 0 = all).");
    for (size_t m = 0; m < c.propModels.size(); ++m) {   // C40: what the geomodel parameters re-simulate
        const PropModelInput& q = c.propModels[m];
        static const char* kPM[] = {"PORO", "PERMX", "NTG", "FACIES"};
        ImGui::TextDisabled("Property model #%zu: %s by %s from %zu hard data%s (%s), ranges %.4g / %.4g / %.4g, seed %u%s", m + 1, kPM[std::clamp(q.property, 0, 3)],
                            q.property == 3 ? "truncated Gaussian" : q.method == 0 ? "kriging" : "SGS", q.data.size(), q.fromWells ? " + well cells" : "", q.source.c_str(), q.rangeMax, q.rangeMin, q.rangeVert, q.seed,
                            q.permFromPoro ? ", PERM from the poro-perm transform" : "");
    }
    auto add = [&](ParamKind k, int target = 0) { StudyParam p; p.kind = k; p.target = target; setDefaults(p, c); P.push_back(p); markStudyDirty(); };
    if (ImGui::Button("+ Parameter")) ImGui::OpenPopup("addparam");
    if (ImGui::BeginPopup("addparam")) {
        for (int k = 0; k < kNumKinds; ++k) {
            const ParamKind pk = ParamKind(k);
            const bool en = !(pk == ParamKind::FaultMult && c.grid.faults.empty()) && !(pk == ParamKind::AquiferStrength && c.aquifers.empty()) &&
                            !((pk == ParamKind::WellSkin || pk == ParamKind::WellPI || pk == ParamKind::InjRate || pk == ParamKind::ProdRate) && c.wells.empty()) &&
                            !(isGeoKind(pk) && c.propModels.empty());
            if (ImGui::MenuItem(kParamKinds[k], nullptr, false, en)) add(pk);
            if (isGeoKind(pk) && c.propModels.empty() && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("Needs a property model in the case: Geomodel > Export > Send to Model, or Model > Grid > Geostatistics > Apply.");
        }
        ImGui::EndPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Typical set")) {
        P.clear();
        for (int k = 0; k < c.grid.nz && k < 6; ++k) { StudyParam p; p.kind = ParamKind::PermXY; setDefaults(p, c); p.k1 = p.k2 = k; p.low = 0.33; p.high = 3; P.push_back(p); }
        add(ParamKind::KrwMax);
        add(ParamKind::WaterVisc);
        for (int f = 0; f < int(c.grid.faults.size()); ++f) add(ParamKind::FaultMult, f);
        for (int a = 0; a < int(c.aquifers.size()); ++a) add(ParamKind::AquiferStrength, a);
        markStudyDirty();
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Layer permeability multipliers, krw end point, water viscosity, faults and aquifers");
    ImGui::SameLine();
    if (!P.empty() && ui::dangerButton("Clear all")) { P.clear(); markStudyDirty(); }
    if (P.empty()) { ui::hint("No parameters yet."); return; }

    int del = -1;
    if (ui::beginTable("params", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY,
                          ImVec2(0, std::min(ImGui::GetContentRegionAvail().y - ImGui::GetFontSize() * 4, ImGui::GetFrameHeightWithSpacing() * (P.size() + 1.6f))))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ui::tableColumn("On", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
        ui::tableColumn("Parameter", 0, 2.6f);
        ui::tableColumn("Target / box", 0, 2.0f);
        ui::tableColumn("Base", 0, 1.0f); ui::tableColumn("Low", 0, 1.0f); ui::tableColumn("High", 0, 1.0f);
        ui::tableColumn("Log", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.6f);
        ui::tableColumn("Label", 0, 1.8f);
        ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.4f);
        ImGui::TableHeadersRow();
        for (int i = 0; i < int(P.size()); ++i) {
            StudyParam& p = P[i];
            ImGui::PushID(i);
            ImGui::TableNextRow();
            bool ch = false;
            ImGui::TableSetColumnIndex(0); ch |= ImGui::Checkbox("##a", &p.active);
            ImGui::TableSetColumnIndex(1); ImGui::SetNextItemWidth(-FLT_MIN);
            int k = int(p.kind);
            if (ImGui::Combo("##k", &k, kParamKinds, kNumKinds)) { p.kind = ParamKind(k); p.target = 0; setDefaults(p, c); ch = true; }
            ImGui::TableSetColumnIndex(2); ImGui::SetNextItemWidth(-FLT_MIN);
            if (needsTarget(p.kind)) {
                std::vector<std::string> names;
                static const char* kPM[] = {"PORO", "PERMX", "NTG", "FACIES"};
                if (p.kind == ParamKind::FaultMult) for (auto& f : c.grid.faults) names.push_back(f.name);
                else if (p.kind == ParamKind::AquiferStrength) for (auto& a : c.aquifers) names.push_back("Aquifer " + std::to_string(a.id));
                else if (isGeoKind(p.kind))
                    for (size_t m = 0; m < c.propModels.size(); ++m)
                        names.push_back(std::string(kPM[std::clamp(c.propModels[m].property, 0, 3)]) + " #" + std::to_string(m + 1) +
                                        (c.propModels[m].source.empty() ? std::string() : " (" + c.propModels[m].source + ")"));
                else for (auto& w : c.wells) names.push_back(w.name);
                const bool allOk = p.kind == ParamKind::InjRate || p.kind == ParamKind::ProdRate || isGeoKind(p.kind);
                const char* allName = p.kind == ParamKind::InjRate ? "All injectors" : isGeoKind(p.kind) ? "All property models" : "All producers";
                const char* cur = (p.target >= 0 && p.target < int(names.size())) ? names[p.target].c_str() : allOk ? allName : "?";
                if (ImGui::BeginCombo("##t", cur)) {
                    if (allOk && ImGui::Selectable(allName, p.target < 0)) { p.target = -1; ch = true; }
                    for (int t = 0; t < int(names.size()); ++t) if (ImGui::Selectable(names[t].c_str(), t == p.target)) { p.target = t; if (p.kind == ParamKind::FaultMult) setDefaults(p, c); ch = true; }
                    ImGui::EndCombo();
                }
            } else if (paramUsesBox(p.kind)) {
                char b[64];
                auto r = [](int a, int z) { return a < 0 ? std::string("*") : (a == z ? std::to_string(a + 1) : std::to_string(a + 1) + "-" + std::to_string(z + 1)); };
                std::snprintf(b, sizeof(b), "I %s  J %s  K %s##box", r(p.i1, p.i2).c_str(), r(p.j1, p.j2).c_str(), r(p.k1, p.k2).c_str());
                if (ImGui::Button(b, ImVec2(-FLT_MIN, 0))) ImGui::OpenPopup("box");
                if (ImGui::BeginPopup("box")) {
                    ImGui::TextDisabled("Box (1-based, 0 = whole range)");
                    int v[3][2] = {{p.i1 + 1, p.i2 + 1}, {p.j1 + 1, p.j2 + 1}, {p.k1 + 1, p.k2 + 1}};
                    const int lim[3] = {c.grid.nx, c.grid.ny, c.grid.nz};
                    const char* nm[3] = {"I range", "J range", "K range"};
                    for (int d = 0; d < 3; ++d) {
                        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                        if (ImGui::DragIntRange2(nm[d], &v[d][0], &v[d][1], 0.1f, 0, lim[d])) {
                            int a = std::clamp(v[d][0], 0, lim[d]), z = std::clamp(v[d][1], 0, lim[d]);
                            if (a == 0 || z == 0) { a = 0; z = 0; }
                            if (z < a) z = a;
                            int* dst[3][2] = {{&p.i1, &p.i2}, {&p.j1, &p.j2}, {&p.k1, &p.k2}};
                            *dst[d][0] = a - 1; *dst[d][1] = z - 1;
                            ch = true;
                        }
                    }
                    ImGui::EndPopup();
                }
            } else ImGui::TextDisabled("global");
            double* vals[3] = {&p.base, &p.low, &p.high};
            for (int q = 0; q < 3; ++q) {
                ImGui::TableSetColumnIndex(3 + q); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(q);
                ch |= ImGui::InputDouble("##v", vals[q], 0, 0, "%.4g");
                ImGui::PopID();
            }
            ImGui::TableSetColumnIndex(6); ch |= ImGui::Checkbox("##l", &p.log);
            ImGui::TableSetColumnIndex(7); ImGui::SetNextItemWidth(-FLT_MIN);
            if (p.label.empty()) {
                const std::string auto_ = paramLabel(p, c);
                char buf[64]; std::snprintf(buf, sizeof(buf), "%s", auto_.c_str());
                ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Muted));
                if (ImGui::InputText("##lab", buf, sizeof(buf))) { p.label = buf; ch = true; }
                ImGui::PopStyleColor();
            } else ch |= ui::inputString("##lab", &p.label);
            ImGui::TableSetColumnIndex(8); if (ImGui::SmallButton("x")) del = i;
            if (p.log && (p.low <= 0 || p.high <= 0)) p.log = false;
            if (ch) markStudyDirty();
            ImGui::PopID();
        }
        ui::endTable();
    }
    if (del >= 0) { P.erase(P.begin() + del); markStudyDirty(); }
    // validation hint
    int bad = 0, act = 0;
    for (auto& p : P) { if (p.active) ++act; if (!(p.low < p.high) || p.base < p.low || p.base > p.high) ++bad; }
    if (bad) ImGui::TextColored(col4(Pal::Warn), "%d parameter(s) need low < base < high", bad);
    else ImGui::TextDisabled("%d active of %d parameters  |  saved with the project and in exported decks (-- RESSIM: PARAM)", act, int(P.size()));
}

// ============================================================================ observed data
void App::plotObservationFit(const char* id, const Case& c, const std::vector<const RunResult*>& runs, const std::vector<const char*>& labels,
                             int obsGroup, float height) {
    const auto groups = obsGroups(c);
    if (groups.empty()) return;
    const ObsGroup& G = groups[std::clamp(obsGroup, 0, int(groups.size()) - 1)];
    const std::vector<double> sig = autoSigma(c.observations);
    std::vector<double> t, v, e;
    for (int i : G.idx) { t.push_back(c.observations[i].time / 365.25); v.push_back(c.observations[i].value); e.push_back(sig[i]); }
    if (ImPlot::BeginPlot(id, ImVec2(-1, height))) {
        ImPlot::SetupAxes("Time (years)", toString(G.qty), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        for (size_t r = 0; r < runs.size(); ++r) {
            const RunResult* R = runs[r];
            if (!R || !R->ok || R->summary.size() < 2) continue;
            std::vector<double> x, y;
            for (size_t k = 1; k < R->summary.size(); ++k) {
                const double tt = R->summary[k].time;
                x.push_back(tt / 365.25);
                y.push_back(simulatedValue(*R, c, G.qty, G.well, tt));
            }
            ImPlotSpec sp;
            const bool ens = labels[r] == nullptr;
            sp.LineColor = ens ? ImVec4(0.6f, 0.65f, 0.75f, 0.25f) : col4(r == 0 ? Pal::Muted : seriesCol(int(r) - 1));
            sp.LineWeight = ens ? 1.0f : (r == 0 ? 1.5f : 2.5f);
            if (ens) sp.Flags = ImPlotItemFlags_NoLegend;
            char lb[32]; std::snprintf(lb, sizeof(lb), "##ens%d", int(r));
            ImPlot::PlotLine(ens ? lb : labels[r], x.data(), y.data(), int(x.size()), sp);
        }
        ImPlotSpec es; es.LineColor = ImVec4(1, 1, 1, 0.55f); es.Size = 3;
        ImPlot::PlotErrorBars("observed", t.data(), v.data(), e.data(), int(t.size()), es);
        ImPlotSpec ms; ms.Marker = ImPlotMarker_Circle; ms.MarkerSize = 3.5f; ms.MarkerFillColor = ImVec4(1, 1, 1, 0.95f);
        ms.MarkerLineColor = ImVec4(0, 0, 0, 1); ms.LineWeight = 0;
        ms.LineColor = ImVec4(1, 1, 1, 0.95f);   // legend swatch
        ImPlot::PlotScatter("observed", t.data(), v.data(), int(t.size()), ms);
        ImPlot::EndPlot();
    }
}

void App::studyObservations() {
    auto& O = case_.observations;
    // last simulation (for synthetic data / comparison)
    RunResult last; Case lastCase; bool haveLast = false;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        if (R.valid && R.summary.size() > 1 && !runner_.running()) { last.ok = true; last.summary = R.summary; lastCase = R.input; haveLast = true; }
    }
    ui::hint("Observed data drive the history-match misfit: sum of ((simulated - observed) / sigma)^2 / N. "
             "sigma <= 0 means automatic (5% of the series range, at least a small floor).");
    if (ImGui::Button("From well history (WCONHIST)")) {
        auto h = observationsFromHistory(case_);
        if (h.empty()) notify("No WCONHIST history rates in the schedule");
        else { O = h; markStudyDirty(); notify("Loaded " + std::to_string(h.size()) + " observations from the history controls"); }
    }
    ImGui::SameLine();
    if (ImGui::Button("Synthetic from last run...")) ImGui::OpenPopup("synth");
    if (ImGui::BeginPopup("synth")) {
        ImGui::TextDisabled("Create observed data from the last simulation (+ Gaussian noise).");
        ImGui::TextDisabled("Classic twin experiment: change parameters afterwards and recover them by history matching.");
        for (int q = 0; q < kNumObs; ++q) {
            bool on = (study_.synthMask >> q) & 1;
            if (q % 7) ImGui::SameLine();
            if (ImGui::Checkbox(kObsNames[q], &on)) study_.synthMask = on ? (study_.synthMask | (1 << q)) : (study_.synthMask & ~(1 << q));
        }
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
        ImGui::SliderFloat("Relative noise", &study_.synthNoise, 0.0f, 0.2f, "%.3f");
        ImGui::BeginDisabled(!haveLast);
        if (ui::accentButton("Create")) {
            std::vector<ObsQty> q;
            for (int i = 0; i < kNumObs; ++i) if ((study_.synthMask >> i) & 1) q.push_back(ObsQty(i));
            O = syntheticObservations(last, lastCase, q, study_.synthNoise, unsigned(study_.seed));
            markStudyDirty();
            notify("Created " + std::to_string(O.size()) + " synthetic observations");
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndDisabled();
        if (!haveLast) { ImGui::SameLine(); ImGui::TextColored(col4(Pal::Warn), "run the simulation first"); }
        ImGui::EndPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("From summary file (SMSPEC)...")) ImGui::OpenPopup("obssum");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("ECLIPSE / OPM summary: history vectors WOPRH, WWPRH, WGPRH, WWIRH, WBHPH, WWCTH, WGORH, FOPRH ...\nor the simulated vectors of another run");
    obsImportSummaryPopup();
    ImGui::SameLine();
    if (ImGui::Button("Import CSV...")) {
        std::string p;
        if (fileDialog(false, "CSV (*.csv;*.txt)\0*.csv;*.txt\0All files (*.*)\0*.*\0", "csv", p)) {
            std::ifstream f(p);
            std::string line; int n = 0, bad = 0;
            std::vector<Observation> in;
            while (std::getline(f, line)) {
                for (char& ch : line) if (ch == ',' || ch == ';' || ch == '\t') ch = ' ';
                std::istringstream ss(line);
                std::string q, w; double t, v, s = -1;
                if (!(ss >> q)) continue;
                if (q[0] == '#' || q == "qty" || q == "QTY") continue;
                for (char& ch : q) ch = char(std::toupper((unsigned char)ch));
                int qi = -1;
                for (int i = 0; i < kNumObs; ++i) if (q == kObsNames[i]) qi = i;
                if (qi < 0) { ++bad; continue; }
                if (!(ss >> w >> t >> v)) { ++bad; continue; }
                ss >> s;
                Observation o; o.qty = ObsQty(qi); o.well = (w == "-" || w == "FIELD") ? "" : w; o.time = t; o.value = v; o.sigma = ss.fail() ? -1 : s;
                in.push_back(o); ++n;
            }
            if (n) { O = in; markStudyDirty(); }
            notify("Imported " + std::to_string(n) + " observations" + (bad ? " (" + std::to_string(bad) + " lines skipped)" : ""));
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Columns: QTY, WELL (- for field), TIME (days), VALUE[, SIGMA]\ne.g.  WBHP, PROD, 365.25, 3120, 25");
    ImGui::SameLine();
    ImGui::BeginDisabled(O.empty());
    if (ImGui::Button("Export CSV...")) {
        std::string p = "observations.csv";
        if (fileDialog(true, "CSV (*.csv)\0*.csv\0", "csv", p)) {
            std::ofstream f(p);
            f << "QTY,WELL,TIME_DAYS,VALUE,SIGMA\n";
            for (auto& o : O) f << toString(o.qty) << ',' << (o.well.empty() ? "-" : o.well) << ',' << o.time << ',' << o.value << ',' << o.sigma << '\n';
            notify("Exported " + std::to_string(O.size()) + " observations");
        }
    }
    ImGui::SameLine();
    if (ui::dangerButton("Clear")) { O.clear(); markStudyDirty(); }
    ImGui::EndDisabled();
    if (O.empty()) { ui::hint("No observed data yet."); return; }

    const auto groups = obsGroups(case_);
    study_.plotObs = std::clamp(study_.plotObs, 0, int(groups.size()) - 1);
    const float leftW = ImGui::GetContentRegionAvail().x * 0.42f;
    ImGui::BeginChild("obsL", ImVec2(leftW, 0), ImGuiChildFlags_ResizeX);
    ImGui::TextDisabled("%d observations in %d series", int(O.size()), int(groups.size()));
    for (int gI = 0; gI < int(groups.size()); ++gI) {
        char b[96]; std::snprintf(b, sizeof(b), "%-14s %4d points", groupName(groups[gI]).c_str(), int(groups[gI].idx.size()));
        if (ImGui::Selectable(b, gI == study_.plotObs)) study_.plotObs = gI;
    }
    ImGui::Separator();
    int del = -1;
    {
        // series weight (relative importance in the objective; 0 = ignore the series)
        const auto& gi = groups[study_.plotObs].idx;
        float w = gi.empty() ? 1.0f : float(O[size_t(gi[0])].weight);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        if (ImGui::DragFloat("Series weight", &w, 0.02f, 0.0f, 100.0f, "%.2f")) { for (int q : gi) O[size_t(q)].weight = std::max(0.0f, w); markStudyDirty(); }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Weight of every observation of this series in the misfit (sigma / sqrt(weight)); 0 ignores the series");
    }
    if (ui::beginTable("obs", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchSame, ImVec2(0, -FLT_MIN))) {   // fills the observation panel
        ImGui::TableSetupScrollFreeze(0, 1);
        for (const char* h : {"Time (d)", "Value", "Sigma", ""}) (void)h;
        ui::tableColumn("Series"); ui::tableColumn("Time (d)"); ui::tableColumn("Value"); ui::tableColumn("Sigma"); ui::tableColumn("Weight");
        ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 1.4f);
        ImGui::TableHeadersRow();
        const auto& idx = groups[study_.plotObs].idx;
        ImGuiListClipper clip; clip.Begin(int(idx.size()));
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                Observation& o = O[idx[r]];
                ImGui::PushID(idx[r]);
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(groupName(groups[study_.plotObs]).c_str());
                double* f[4] = {&o.time, &o.value, &o.sigma, &o.weight};
                for (int q = 0; q < 4; ++q) {
                    ImGui::TableSetColumnIndex(1 + q); ImGui::SetNextItemWidth(-FLT_MIN); ImGui::PushID(q);
                    if (ImGui::InputDouble("##v", f[q], 0, 0, "%.5g")) markStudyDirty();
                    ImGui::PopID();
                }
                ImGui::TableSetColumnIndex(5); if (ImGui::SmallButton("x")) del = idx[r];
                ImGui::PopID();
            }
        ui::endTable();
    }
    if (del >= 0) { O.erase(O.begin() + del); markStudyDirty(); }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("obsR", ImVec2(0, 0));
    std::vector<const RunResult*> runs; std::vector<const char*> labels;
    if (haveLast) { runs.push_back(&last); labels.push_back("last run"); }
    plotObservationFit("Observed vs simulated", case_, runs, labels, study_.plotObs, std::max(ImGui::GetContentRegionAvail().y - 4, ImGui::GetFontSize() * 12));
    ImGui::EndChild();
}

// ============================================================================ sensitivity
void App::studySensitivity() {
    const auto& P = case_.studyParams;
    int act = 0; for (auto& p : P) if (p.active) ++act;
    const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
    if (ui::beginProps("sensopt", 0.36f)) {
        ui::propCombo("Response", &study_.response, kResp, ressim::kNumResponses, "Quantity whose sensitivity / uncertainty is analysed");
        ui::propI("Parallel runs (0 = auto)", &study_.workers, 0, hw * 2);
        if (Response(study_.response) == Response::NPV) {
            ui::propD("Oil price / stb", &study_.eco.oilPrice, "%.4g");
            ui::propD("Gas price / Mscf", &study_.eco.gasPrice, "%.4g");
            ui::propD("Water handling / stb", &study_.eco.waterCost, "%.4g");
            ui::propD("Injection cost / stb", &study_.eco.injCost, "%.4g");
            ui::propD("Discount rate / year", &study_.eco.discount, "%.3g");
        }
        ui::endProps();
    }
    if (act == 0) { ui::hint("Define and activate parameters in the Parameters tab first."); return; }
    if (Response(study_.response) == Response::Misfit && case_.observations.empty()) ImGui::TextColored(col4(Pal::Warn), "Misfit needs observed data.");
    const bool busy = study_.running;
    ImGui::BeginDisabled(busy);
    char b[96];
    std::snprintf(b, sizeof(b), "Run tornado  (%d runs)", 2 * act + 1);
    if (ui::accentButton(b)) {
        const Case base = case_; const Response r = Response(study_.response); const Economics eco = study_.eco;
        study_.deferred = [this, base, r, eco]() {
            startStudy("Tornado sensitivity", [base, r, eco](StudyState& S, const StudyCallbacks& cb) {
                TornadoResult t = runTornado(base, base.studyParams, r, eco, cb);
                std::lock_guard<std::mutex> l(S.m);
                S.tornado = std::move(t); S.tornadoResp = r; S.hasTornado = true;
            });
        };
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
    ImGui::DragInt("##mcn", &study_.mcSamples, 1, 8, 2000, "%d samples");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
    ImGui::DragInt("##seed", &study_.seed, 1, 1, 1 << 30, "seed %d");
    ImGui::SameLine();
    if (ui::accentButton("Run Monte Carlo (Latin hypercube)")) {
        const Case base = case_; const Response r = Response(study_.response); const Economics eco = study_.eco;
        const int n = std::max(4, study_.mcSamples); const unsigned seed = unsigned(study_.seed);
        study_.deferred = [this, base, r, eco, n, seed]() {
            startStudy("Monte Carlo uncertainty", [base, r, eco, n, seed](StudyState& S, const StudyCallbacks& cb) {
                MonteCarloResult m = runMonteCarlo(base, base.studyParams, n, seed, r, eco, cb);
                std::lock_guard<std::mutex> l(S.m);
                S.mc = std::move(m); S.mcResp = r; S.hasMC = true;
            });
        };
    }
    ImGui::EndDisabled();

    const float fs = ImGui::GetFontSize();
    if (ui::beginTable("sensres", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableNextRow();
        // ---- tornado
        ImGui::TableSetColumnIndex(0);
        ui::header("Tornado chart");
        if (study_.hasTornado) {
            const TornadoResult& T = study_.tornado;
            std::vector<std::string> lab; std::vector<double> lo, hi;
            for (auto& bar : T.bars) {
                const StudyParam& p = (bar.param < int(P.size())) ? P[bar.param] : StudyParam{};
                lab.push_back(bar.param < int(P.size()) ? paramLabel(p, case_) + "  [" + fmtv(p.low) + " .. " + fmtv(p.high) + "]" : "?");
                lo.push_back(bar.okLow ? bar.low : NAN); hi.push_back(bar.okHigh ? bar.high : NAN);
            }
            ImGui::TextDisabled("%s  -  base %.5g", kResp[int(study_.tornadoResp)], T.base);
            barChart("##tornado", lab, lo, hi, T.base, "parameter at low", "parameter at high", false);
            ui::hint("Each parameter is moved to its low and high value with all others at base (one-at-a-time). Bars are sorted by swing.");
        } else ImGui::TextDisabled("Run the tornado to rank the parameters.");
        // ---- Monte Carlo
        ImGui::TableSetColumnIndex(1);
        ui::header("Monte Carlo (Latin hypercube)");
        if (study_.hasMC) {
            const MonteCarloResult& M = study_.mc;
            std::vector<double> ok;
            for (double v : M.response) if (std::isfinite(v)) ok.push_back(v);
            ImGui::Text("%s:  P90 %.5g   P50 %.5g   P10 %.5g   mean %.5g", kResp[int(study_.mcResp)], M.p90, M.p50, M.p10, M.mean);
            ImGui::TextDisabled("%d runs, %d failed  |  P90 = 90%% probability of exceeding (low case)", int(M.response.size()), M.failed);
            const float h = fs * 12;
            if (ImPlot::BeginSubplots("##mc", 1, 2, ImVec2(-1, h))) {
                if (ImPlot::BeginPlot("Distribution")) {
                    ImPlot::SetupAxes(kResp[int(study_.mcResp)], "runs", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec sp; sp.FillColor = col4(Pal::Accent); sp.FillAlpha = 0.7f; sp.LineColor = col4(Pal::Accent);
                    if (!ok.empty()) ImPlot::PlotHistogram("##h", ok.data(), int(ok.size()), std::max(6, int(std::sqrt(double(ok.size())))), 1.0, ImPlotRange(), sp);
                    double xs[3] = {M.p90, M.p50, M.p10};
                    ImPlotSpec ls; ls.LineColor = ImVec4(1, 1, 1, 0.7f); ls.Flags = ImPlotItemFlags_NoLegend;
                    ImPlot::PlotInfLines("##p", xs, 3, ls);
                    ImPlot::PlotText("P90", M.p90, 0, ImVec2(0, -10)); ImPlot::PlotText("P50", M.p50, 0, ImVec2(0, -10)); ImPlot::PlotText("P10", M.p10, 0, ImVec2(0, -10));
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot("FOPT fan (P10-P50-P90)")) {
                    ImPlot::SetupAxes("years", "FOPT", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    if (!M.times.empty()) {
                        std::vector<double> ty(M.times.size());
                        for (size_t i = 0; i < ty.size(); ++i) ty[i] = M.times[i] / 365.25;
                        ImPlotSpec sh; sh.FillColor = col4(Pal::Oil); sh.FillAlpha = 0.25f;
                        ImPlot::PlotShaded("P90-P10", ty.data(), M.fanP90.data(), M.fanP10.data(), int(ty.size()), sh);
                        ImPlotSpec l; l.LineColor = col4(Pal::Oil); l.LineWeight = 2.5f;
                        ImPlot::PlotLine("P50", ty.data(), M.fanP50.data(), int(ty.size()), l);
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            ui::header("Rank correlation (Spearman) with the response");
            std::vector<std::string> lab; std::vector<double> neg, pos;
            std::vector<int> order(M.spearman.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return std::fabs(M.spearman[a]) > std::fabs(M.spearman[b]); });
            std::vector<int> activeIdx;
            for (int i = 0; i < int(P.size()); ++i) if (P[i].active) activeIdx.push_back(i);
            for (int o : order) {
                const int pi = (o < int(activeIdx.size())) ? activeIdx[o] : o;
                lab.push_back(pi < int(P.size()) ? paramLabel(P[pi], case_) : "?");
                const double s = M.spearman[o];
                neg.push_back(s < 0 ? s : NAN); pos.push_back(s >= 0 ? s : NAN);
            }
            barChart("##spear", lab, neg, pos, 0.0, "negative", "positive", true);
        } else ImGui::TextDisabled("Run Monte Carlo for P10/P50/P90, the FOPT fan chart and rank correlations.");
        ui::endTable();
    }
}

// ============================================================================ history matching
void App::studyHistoryMatch() {
    auto& P = case_.studyParams;
    int act = 0; for (auto& p : P) if (p.active) ++act;
    if (case_.observations.empty()) { ui::hint("Add observed data first (Observed data tab: WCONHIST history, CSV import or a synthetic twin experiment)."); return; }
    if (act == 0) { ui::hint("Define and activate parameters in the Parameters tab first."); return; }
    const float fs = ImGui::GetFontSize();
    static int method = 0;
    ImGui::RadioButton("Manual (sliders, live misfit)", &method, 0); ImGui::SameLine();
    ImGui::RadioButton("Levenberg-Marquardt", &method, 1); ImGui::SameLine();
    ImGui::RadioButton("ES-MDA ensemble", &method, 2);
    const auto groups = obsGroups(case_);
    study_.plotObs = std::clamp(study_.plotObs, 0, std::max(0, int(groups.size()) - 1));

    if (ui::beginTable("hm", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ui::tableColumn("a", ImGuiTableColumnFlags_WidthStretch, 0.42f);
        ui::tableColumn("b", ImGuiTableColumnFlags_WidthStretch, 0.58f);
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        const RunResult* bestRun = nullptr; const RunResult* initRun = nullptr;
        std::vector<const RunResult*> ens;
        if (method == 0) {
            ui::header("Parameters");
            bool ch = false;
            if (ui::beginProps("man", 0.45f)) {
                for (int i = 0; i < int(P.size()); ++i) {
                    if (!P[i].active) continue;
                    ui::label(paramLabel(P[i], case_).c_str());
                    ImGui::PushID(i);
                    ch |= paramSlider("##s", P[i], &study_.manualVals[i]);
                    ImGui::PopID();
                }
                ui::endProps();
            }
            if (ch) { study_.manualEditTime = ImGui::GetTime(); if (study_.manualAuto) study_.manualPending = true; }
            ImGui::Checkbox("Live update", &study_.manualAuto);
            ImGui::SameLine();
            if (ImGui::Button("Run now")) study_.deferred = [this]() { runManual(); };
            ImGui::SameLine();
            if (ImGui::Button("Reset to base")) { for (size_t i = 0; i < P.size(); ++i) study_.manualVals[i] = P[i].base; study_.manualPending = true; }
            ImGui::Spacing();
            if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
            if (study_.manualRunning) ImGui::TextColored(col4(Pal::Accent), "Misfit: running...");
            else if (study_.hasManual) ImGui::TextColored(col4(std::isfinite(study_.manualMisfit) ? Pal::Ok : Pal::Err), "Misfit: %.4g", study_.manualMisfit);
            else ImGui::TextDisabled("Misfit: -");
            if (g_fonts.title) ImGui::PopFont();
            ImGui::TextDisabled("1 = matched within the data uncertainty on average");
            if (study_.hasManual && ImGui::Button("Apply these values to the case")) {
                const std::vector<double> v = study_.manualVals;
                study_.deferred = [this, v]() {
                    case_ = applyParams(case_, case_.studyParams, v);
                    for (size_t i = 0; i < case_.studyParams.size(); ++i) {
                        StudyParam& p = case_.studyParams[i];
                        if (paramIsMultiplier(p.kind)) { p.base = 1; } else p.base = v[i];
                        p.low = std::min(p.low, p.base); p.high = std::max(p.high, p.base);
                    }
                    study_.manualVals.clear();
                    markDirty(true);
                    notify("Matched parameters applied to the case - re-run the simulation (F5)");
                };
            }
            if (study_.hasManual) bestRun = &study_.manual;
        } else {
            ui::header(method == 1 ? "Levenberg-Marquardt (gradient based)" : "ES-MDA (ensemble smoother, multiple data assimilation)");
            if (ui::beginProps("hmopt", 0.55f)) {
                if (method == 1) {
                    ui::propI("Max iterations", &study_.lm.maxIter, 1, 100);
                    ui::propD("Finite-difference step (unit space)", &study_.lm.fdStep, "%.3g");
                    ui::propD("Initial damping lambda", &study_.lm.lambda0, "%.3g");
                    ui::propD("Relative misfit tolerance", &study_.lm.tol, "%.3g");
                } else {
                    ui::propI("Ensemble size", &study_.es.ensemble, 4, 400);
                    ui::propI("Assimilations", &study_.es.assimilations, 1, 12);
                    ui::propD("Prior spread (logit space)", &study_.es.priorStd, "%.3g");
                    int sd = int(study_.es.seed); if (ui::propI("Seed", &sd, 1, 1 << 30)) study_.es.seed = unsigned(sd);
                }
                const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
                ui::propI("Parallel runs (0 = auto)", &study_.workers, 0, hw * 2);
                ui::endProps();
            }
            char b[96];
            if (method == 1) std::snprintf(b, sizeof(b), "Run history match  (~%d runs / iteration)", act + 4);
            else std::snprintf(b, sizeof(b), "Run history match  (%d runs)", study_.es.ensemble * (study_.es.assimilations + 1));
            ImGui::BeginDisabled(study_.running);
            if (ui::accentButton(b, ImVec2(-FLT_MIN, 0))) {
                const Case base = case_; const int m = method; const LMOptions lo = study_.lm; const ESMDAOptions eo = study_.es;
                study_.deferred = [this, base, m, lo, eo]() {
                    startStudy(m == 1 ? "Levenberg-Marquardt history match" : "ES-MDA history match", [base, m, lo, eo](StudyState& S, const StudyCallbacks& cb) {
                        HMResult r = m == 1 ? runLevenbergMarquardt(base, base.studyParams, lo, cb) : runESMDA(base, base.studyParams, eo, cb);
                        std::lock_guard<std::mutex> l(S.m);
                        S.hm = std::move(r); S.hmMethod = m - 1; S.hmParams = base.studyParams; S.hasHM = true;
                    });
                };
            }
            ImGui::EndDisabled();
            if (study_.hasHM) {
                const HMResult& H = study_.hm;
                const auto& HP = study_.hmParams;
                ImGui::Spacing();
                if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
                ImGui::TextColored(col4(Pal::Ok), "%s: misfit %.4g -> %.4g  (%d runs)", hmName(true), H.initialMisfit, H.bestMisfit, H.runs);
                if (g_fonts.bold) ImGui::PopFont();
                // parameter table
                if (ui::beginTable("hmp", study_.hmMethod == 1 ? 5 : 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                    ui::tableColumn("Parameter", 0, 2.2f); ui::tableColumn("Initial"); ui::tableColumn("Best");
                    if (study_.hmMethod == 1) { ui::tableColumn("Posterior mean"); ui::tableColumn("P90 - P10"); }
                    ImGui::TableHeadersRow();
                    for (int i = 0; i < int(HP.size()) && i < int(H.best.size()); ++i) {
                        if (!HP[i].active) continue;
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(paramLabel(HP[i], case_).c_str());
                        ImGui::TableSetColumnIndex(1); ImGui::Text("%.4g", HP[i].base);
                        ImGui::TableSetColumnIndex(2); ImGui::TextColored(col4(Pal::Accent), "%.4g", H.best[i]);
                        if (study_.hmMethod == 1) {
                            std::vector<double> v;
                            for (size_t j = 0; j < H.ensemble.size(); ++j)
                                if (j < H.ensembleMisfit.size() && std::isfinite(H.ensembleMisfit[j]) && i < int(H.ensemble[j].size())) v.push_back(H.ensemble[j][i]);
                            if (!v.empty()) {
                                std::sort(v.begin(), v.end());
                                const double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
                                ImGui::TableSetColumnIndex(3); ImGui::Text("%.4g", mean);
                                ImGui::TableSetColumnIndex(4); ImGui::Text("%.4g - %.4g", v[size_t(0.1 * (v.size() - 1))], v[size_t(0.9 * (v.size() - 1))]);
                            }
                        }
                    }
                    ui::endTable();
                }
                if (ImGui::Button("Apply best match to the case")) {
                    const std::vector<double> v = H.best; const std::vector<StudyParam> hp = HP;
                    study_.deferred = [this, v, hp]() {
                        case_ = applyParams(case_, hp, v);
                        case_.studyParams = hp;
                        for (size_t i = 0; i < hp.size() && i < v.size(); ++i) {
                            StudyParam& p = case_.studyParams[i];
                            p.base = paramIsMultiplier(p.kind) ? 1.0 : v[i];
                            p.low = std::min(p.low, p.base); p.high = std::max(p.high, p.base);
                        }
                        study_.manualVals.clear();
                        markDirty(true);
                        notify("Best match applied to the case - re-run the simulation (F5)");
                    };
                }
                ImGui::SameLine();
                if (ImGui::Button("Load into manual sliders")) {
                    for (size_t i = 0; i < P.size() && i < H.best.size(); ++i) study_.manualVals[i] = H.best[i];
                    method = 0; study_.manualPending = true;
                }
                // convergence plot
                if (!H.history.empty() && ImPlot::BeginPlot("Misfit convergence", ImVec2(-1, fs * 10))) {
                    ImPlot::SetupAxes(study_.hmMethod == 0 ? "iteration" : "assimilation", "misfit", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                    std::vector<double> x, y;
                    for (auto& it : H.history) { x.push_back(it.iter); y.push_back(std::max(1e-12, it.misfit)); }
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Accent); sp.LineWeight = 2.5f; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 4;
                    sp.MarkerFillColor = col4(Pal::Accent);
                    ImPlot::PlotLine(study_.hmMethod == 0 ? "misfit" : "ensemble mean", x.data(), y.data(), int(x.size()), sp);
                    ImPlot::EndPlot();
                }
                initRun = &H.initialRun; bestRun = &H.bestRun;
                if (study_.hmMethod == 1) for (auto& r : H.ensembleRuns) ens.push_back(&r);
            }
        }
        // ---- right column: data fit
        ImGui::TableSetColumnIndex(1);
        ui::header("Match quality");
        ImGui::SetNextItemWidth(fs * 12);
        if (!groups.empty() && ImGui::BeginCombo("Series", groupName(groups[study_.plotObs]).c_str())) {
            for (int g = 0; g < int(groups.size()); ++g) if (ImGui::Selectable(groupName(groups[g]).c_str(), g == study_.plotObs)) study_.plotObs = g;
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("##prev", ImGuiDir_Left)) study_.plotObs = (study_.plotObs + int(groups.size()) - 1) % std::max(1, int(groups.size()));
        ImGui::SameLine();
        if (ImGui::ArrowButton("##next", ImGuiDir_Right)) study_.plotObs = (study_.plotObs + 1) % std::max(1, int(groups.size()));
        std::vector<const RunResult*> runs; std::vector<const char*> labels;
        // run 0 is drawn muted (initial), then ensemble (no label), then best
        RunResult baseRun;
        if (!initRun) {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& R = runner_.results();
            if (R.valid && R.summary.size() > 1 && !runner_.running()) { baseRun.ok = true; baseRun.summary = R.summary; initRun = &baseRun; }
        }
        runs.push_back(initRun); labels.push_back(initRun == &baseRun ? "last simulation" : "initial");
        for (auto* e : ens) { runs.push_back(e); labels.push_back(nullptr); }
        if (bestRun) { runs.push_back(bestRun); labels.push_back(method == 0 ? "manual" : (study_.hmMethod == 1 ? "best member" : "best match")); }
        plotObservationFit("##fit", case_, runs, labels, study_.plotObs, std::max(fs * 14, ImGui::GetContentRegionAvail().y - fs * 9));
        // per-series misfit table
        if (bestRun && bestRun->ok) {
            const std::vector<double> sig = autoSigma(case_.observations);
            const std::vector<double> rb = residuals(*bestRun, case_, case_.observations, sig);
            std::vector<double> ri;
            if (initRun && initRun->ok) ri = residuals(*initRun, case_, case_.observations, sig);
            if (ui::beginTable("mis", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, fs * 8))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ui::tableColumn("Series"); ui::tableColumn(initRun ? "Misfit before" : "-"); ui::tableColumn("Misfit now");
                ImGui::TableHeadersRow();
                for (int g = 0; g < int(groups.size()); ++g) {
                    double a = 0, b = 0;
                    for (int i : groups[g].idx) { if (i < int(rb.size())) b += rb[i] * rb[i]; if (i < int(ri.size())) a += ri[i] * ri[i]; }
                    a /= groups[g].idx.size(); b /= groups[g].idx.size();
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    if (ImGui::Selectable(groupName(groups[g]).c_str(), g == study_.plotObs, ImGuiSelectableFlags_SpanAllColumns)) study_.plotObs = g;
                    ImGui::TableSetColumnIndex(1); if (!ri.empty()) ImGui::Text("%.4g", a);
                    ImGui::TableSetColumnIndex(2); ImGui::TextColored(col4(b <= 1.5 ? Pal::Ok : b <= 10 ? Pal::Warn : Pal::Err), "%.4g", b);
                }
                ui::endTable();
            }
        }
        ui::endTable();
    }
}

// ============================================================================ Studies workspace
// Full-screen sensitivity / history-matching screen: setup (left), result views (centre), run archive (right).
namespace {
const char* kMethods[] = {"Tornado (one at a time)", "Monte Carlo (Latin hypercube)", "Experimental design + proxy",
                          "Levenberg-Marquardt", "ES-MDA ensemble", "Manual (sliders)", "Adjoint gradient (L-BFGS / per-cell field)"};
const char* kMethodIcons[] = {ICON_CHART_COLUMN, ICON_SHUFFLE, ICON_GRID_3X3, ICON_TRENDING_UP, ICON_SPLIT, ICON_SLIDERS_HORIZONTAL, ICON_ACTIVITY};
const char* kDesigns[] = {"Latin hypercube", "Two-level full factorial", "Central composite", "Plackett-Burman", "Random"};
const char* kGroupings[] = {"Quantity (WOPR, WBHP, ...)", "Well", "Quantity + well", "Phase class (oil / water / gas / pressure)"};
const char* kOvVec[] = {"FOPR", "FWPR", "FGPR", "FWIR", "FOPT", "FWPT", "FGPT", "FWIT", "FPR", "FWCT", "FGOR",
                        "WOPR", "WWPR", "WGPR", "WWIR", "WBHP", "WWCT", "WGOR"};
const int kNumOv = 18, kOvFirstWell = 11;

int ovQty(int v) {   // overlay vector -> ObsQty index (-1 = none)
    static const int m[] = {int(ObsQty::FOPR), int(ObsQty::FWPR), int(ObsQty::FGPR), -1, -1, -1, -1, -1, int(ObsQty::FPR), int(ObsQty::FWCT), int(ObsQty::FGOR),
                            int(ObsQty::WOPR), int(ObsQty::WWPR), int(ObsQty::WGPR), int(ObsQty::WWIR), int(ObsQty::WBHP), int(ObsQty::WWCT), int(ObsQty::WGOR)};
    return v >= 0 && v < kNumOv ? m[v] : -1;
}
bool ovIsRate(int v) { return v <= 3 || (v >= kOvFirstWell && v != 15) || v == 9 || v == 10; }
double ovValue(const RunResult& r, const Case& c, int v, const std::string& well, const SummaryPoint& s) {
    switch (v) {
        case 0: return s.fopr; case 1: return s.fwpr; case 2: return s.fgpr; case 3: return s.fwir;
        case 4: return s.fopt; case 5: return s.fwpt; case 6: return s.fgpt; case 7: return s.fwit;
        case 8: return s.fpr; case 9: return s.fwct; case 10: return s.fgor;
        default: return simulatedValue(r, c, ObsQty(ovQty(v)), well, s.time);
    }
}
ImU32 tagColor(const std::string& t) {
    if (t == "tornado") return IM_COL32(100, 160, 255, 255);
    if (t == "MC") return IM_COL32(160, 120, 240, 255);
    if (t == "LM") return IM_COL32(255, 170, 60, 255);
    if (t == "ES-MDA") return IM_COL32(80, 210, 200, 255);
    if (t == "design") return IM_COL32(120, 200, 100, 255);
    if (t == "proxy") return IM_COL32(255, 220, 80, 255);
    if (t == "manual") return IM_COL32(240, 110, 160, 255);
    return IM_COL32(190, 190, 200, 255);
}
double pct(std::vector<double> v, double q) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const double x = q * (v.size() - 1);
    const size_t a = size_t(x), b = std::min(v.size() - 1, a + 1);
    return v[a] + (v[b] - v[a]) * (x - a);
}
} // namespace

void App::applyStudyValues(const std::vector<StudyParam>& hp, const std::vector<double>& v, const std::string& what) {
    study_.deferred = [this, v, hp, what]() {
        case_ = applyParams(case_, hp, v);
        case_.studyParams = hp;
        for (size_t i = 0; i < hp.size() && i < v.size(); ++i) {
            StudyParam& p = case_.studyParams[i];
            p.base = paramIsMultiplier(p.kind) ? 1.0 : v[i];
            p.low = std::min(p.low, p.base); p.high = std::max(p.high, p.base);
        }
        study_.manualVals.clear();
        markDirty(true);
        notify(what + " applied to the case - re-run the simulation (F5)");
    };
}

void App::studiesWorkspace() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0);
    ImGui::Begin("##studiesws", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoDocking |
                                             ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoSavedSettings);
    ImGui::PopStyleVar();
    const float fs = ImGui::GetFontSize();
    if (study_.manualVals.size() != case_.studyParams.size()) {
        study_.manualVals.resize(case_.studyParams.size());
        for (size_t i = 0; i < case_.studyParams.size(); ++i) study_.manualVals[i] = case_.studyParams[i].base;
    }
    if (study_.manualPending && !study_.manualRunning && ImGui::GetTime() - study_.manualEditTime > 0.25) runManual();

    // ---- header
    if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
    ImGui::TextUnformatted(ICON_FLASK_CONICAL "  Studies");
    if (g_fonts.title) ImGui::PopFont();
    ImGui::SameLine(0, fs * 1.5f);
    {   // mode switch (segmented control)
        static const char* modes[] = {ICON_CHART_COLUMN "  Sensitivity & uncertainty", ICON_TARGET "  History matching"};
        for (int m = 0; m < 2; ++m) {
            if (m) ImGui::SameLine(0, 2);
            const bool on = study_.wsMode == m;
            if (on) ImGui::PushStyleColor(ImGuiCol_Button, col4(Pal::Accent));
            if (ImGui::Button(modes[m])) {
                if (!on) {
                    study_.wsMode = m;
                    study_.wsMethod = m == 1 ? 3 : 0;
                    study_.centerTab = m == 1 ? (case_.observations.empty() ? 8 : 9) : 0;
                }
            }
            if (on) ImGui::PopStyleColor();
        }
    }
    ImGui::SameLine(0, fs * 2);
    {
        int act = 0; for (auto& p : case_.studyParams) if (p.active) ++act;
        size_t nrun; { std::lock_guard<std::mutex> l(study_.m); nrun = study_.archive.size(); }
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("%d active parameters  |  %d observations  |  %zu archived runs", act, int(case_.observations.size()), nrun);
    }
    if (study_.running) {
        ImGui::SameLine(0, fs * 2);
        const int d = study_.done, t = study_.total;
        std::string task; { std::lock_guard<std::mutex> l(study_.m); task = study_.task; }
        char b[160]; std::snprintf(b, sizeof(b), "%s  -  %d / %d runs  (%.0f s)", task.c_str(), d, t, ImGui::GetTime() - study_.t0);
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col4(Pal::Accent));
        ImGui::ProgressBar(t > 0 ? float(d) / t : 0.0f, ImVec2(std::max(fs * 10, ImGui::GetContentRegionAvail().x - fs * 7), 0), b);
        ImGui::PopStyleColor();
        ImGui::SameLine();
        if (ui::dangerButton(ICON_CIRCLE_STOP " Cancel", ImVec2(-FLT_MIN, 0))) study_.cancel = true;
    }
    ImGui::Separator();

    // ---- three panes
    const float H = ImGui::GetContentRegionAvail().y;
    ImGui::BeginChild("##stleft", ImVec2(fs * 25, H), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    studiesSetup();
    ImGui::EndChild();
    ImGui::SameLine();
    const float rightW = fs * 23;
    ImGui::BeginChild("##stcenter", ImVec2(std::max(fs * 20, ImGui::GetContentRegionAvail().x - rightW - ImGui::GetStyle().ItemSpacing.x), H), ImGuiChildFlags_Borders);
    {
        std::unique_lock<std::mutex> lk(study_.m);
        studiesCenter();
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("##stright", ImVec2(0, H), ImGuiChildFlags_Borders);
    {
        std::unique_lock<std::mutex> lk(study_.m);
        studiesRuns();
    }
    ImGui::EndChild();
    ImGui::End();
    if (study_.deferred) { auto f = std::move(study_.deferred); study_.deferred = nullptr; f(); }
}

// ---------------------------------------------------------------------------- setup (left)
void App::studiesSetup() {
    const float fs = ImGui::GetFontSize();
    if (!ImGui::BeginTabBar("##setup")) return;
    if (ui::tabItem(ICON_PLAY " Method")) {
        ImGui::BeginChild("##meth");
        const auto& P = case_.studyParams;
        int act = 0; for (auto& p : P) if (p.active) ++act;
        ui::header("Method");
        if (study_.wsMode == 1 && study_.wsMethod < 3) study_.wsMethod = 3;
        if (study_.wsMode == 0 && study_.wsMethod >= 3) study_.wsMethod = 0;
        static const int kOrderSA[] = {0, 1, 2}, kOrderHM[] = {3, 6, 4, 5};
        const int* order = study_.wsMode == 1 ? kOrderHM : kOrderSA;
        for (int oi = 0; oi < (study_.wsMode == 1 ? 4 : 3); ++oi) {
            const int m = order[oi];
            char b[96]; std::snprintf(b, sizeof(b), "%s  %s", kMethodIcons[m], kMethods[m]);
            if (ImGui::Selectable(b, study_.wsMethod == m)) study_.wsMethod = m;
        }
        const int M = study_.wsMethod;
        ui::header("Options");
        const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
        if (ui::beginProps("wsopt", 0.45f)) {
            if (M <= 2 || M == 5) ui::propCombo("Response", &study_.response, kResp, ressim::kNumResponses, "Quantity analysed by sensitivity / optimised by the proxy");
            if (Response(study_.response) == Response::NPV && (M <= 2)) {
                ui::propD("Oil price / stb", &study_.eco.oilPrice, "%.4g");
                ui::propD("Gas price / Mscf", &study_.eco.gasPrice, "%.4g");
                ui::propD("Water handling / stb", &study_.eco.waterCost, "%.4g");
                ui::propD("Injection cost / stb", &study_.eco.injCost, "%.4g");
                ui::propD("Discount rate / year", &study_.eco.discount, "%.3g");
            }
            if (M == 1) { ui::propI("Samples", &study_.mcSamples, 4, 4000); ui::propI("Seed", &study_.seed, 1, 1 << 30); }
            if (M == 2) {
                ui::propCombo("Design", &study_.designKind, kDesigns, 5, "Experimental design of the training runs");
                if (study_.designKind == 0 || study_.designKind == 4) ui::propI("Runs", &study_.designN, 4, 2000);
                ui::propI("Seed", &study_.seed, 1, 1 << 30);
                {
                    static const char* kPk[] = {"Polynomial response surface", "Gaussian process (kriging / ML)"};
                    ui::propCombo("Proxy model", &study_.proxyKind, kPk, 2, "C28: Gaussian process = CMOST-style ML proxy - interpolates the training runs, "
                                  "ARD length scales by maximum likelihood, gives a prediction uncertainty; best for non-polynomial responses");
                }
                if (study_.proxyKind == 0) ui::propBool("Quadratic proxy", &study_.proxyQuad, "Quadratic response surface (else linear)");
                ui::propBool("Maximise response", &study_.proxyMax, "Optimum of the proxy: maximum (else minimum)");
            }
            if (M == 3) {
                static const char* kOpt[] = {"Finite-difference Jacobian",
                                             "Analytic Jacobian (forward sensitivities)"};
                study_.hmOpt = std::min(study_.hmOpt, 1);
                ui::propCombo("Jacobian", &study_.hmOpt, kOpt, 2,
                              "Finite differences: one perturbed run per parameter per iteration.\n"
                              "Analytic Jacobian: exact discrete forward sensitivities of every observation (one recorded "
                              "run + linearised replays, no perturbed runs).\n"
                              "The adjoint-gradient methods (scalar L-BFGS and per-cell field) have their own entry: "
                              "'Adjoint gradient' in the method list.");
                {
                    ui::propI("Max iterations", &study_.lm.maxIter, 1, 100);
                    if (study_.hmOpt == 0) ui::propD("Finite-difference step", &study_.lm.fdStep, "%.3g", "Normalised parameter step of the Jacobian columns");
                    else ui::propD("Residual perturbation", &study_.lm.resStep, "%.3g", "Relative parameter perturbation used to differentiate the discrete residual (central differences of the residual only, not of runs)");
                    ui::propD("Initial damping", &study_.lm.lambda0, "%.3g");
                    ui::propD("Tolerance", &study_.lm.tol, "%.3g");
                    study_.lm.gradient = study_.hmOpt == 1 ? GradientMethod::Analytic : GradientMethod::FiniteDifference;
                }
            }
            if (M == 6) {
                static const char* kAdjModes[] = {"Scalar parameters (Parameters tab) - projected L-BFGS",
                                                  "Per-cell field (every active cell: ln k / ln phi multipliers)"};
                ui::propCombo("Unknowns", &study_.adjMode, kAdjModes, 2,
                              "Both use the exact discrete adjoint: one forward + one backward pass per iteration gives the "
                              "misfit gradient for all unknowns, so the cost does not grow with their number.\n"
                              "Scalar: the active study parameters (multipliers, contacts, kr end points, ...).\n"
                              "Per-cell: a log-multiplier of the permeability (all directions) and / or porosity of every "
                              "active cell, regularised by a smoothness term.");
                if (study_.adjMode == 0) {
                    ui::propI("Max iterations", &study_.adj.maxIter, 1, 200);
                    ui::propI("L-BFGS memory", &study_.adj.memory, 1, 30, "Number of curvature pairs kept");
                    ui::propD("Max step (normalised)", &study_.adj.maxStep, "%.3g", "Largest parameter move per iteration in the normalised [0, 1] box");
                    ui::propD("Residual perturbation", &study_.adj.resStep, "%.3g");
                    ui::propD("Tolerance", &study_.adj.tol, "%.3g", "Stop when the projected gradient / relative misfit change falls below this");
                } else {
                    FieldHMOptions& fo = study_.fieldOpt;
                    ui::propBool("Match permeability (kx, ky, kz)", &fo.perm);
                    ui::propBool("Match porosity", &fo.poro);
                    ui::propI("Iterations", &fo.maxIter, 1, 200);
                    ui::propD("Smoothness weight", &fo.regWeight, "%.3g", "weight of the mean squared ln-multiplier jump between neighbouring cells");
                    ui::propD("Lowest multiplier", &fo.lo, "%.3g");
                    ui::propD("Highest multiplier", &fo.hi, "%.3g");
                    ui::propD("Max ln step / iteration", &fo.maxLogStep, "%.3g");
                    fo.lo = std::max(1e-4, std::min(fo.lo, 1.0)); fo.hi = std::max(1.0, fo.hi);
                }
            }
            if (M == 4) {
                ui::propI("Ensemble size", &study_.es.ensemble, 4, 400);
                ui::propI("Assimilations", &study_.es.assimilations, 1, 12);
                ui::propD("Prior spread (logit space)", &study_.es.priorStd, "%.3g");
                int sd = int(study_.es.seed); if (ui::propI("Seed", &sd, 1, 1 << 30)) study_.es.seed = unsigned(sd);
            }
            if (M == 3 || M == 4 || M == 5 || (M == 6 && study_.adjMode == 0)) ui::propCombo("Objective grouping", &study_.grouping, kGroupings, 4, "How observations are grouped into objective terms (cross-plots, Pareto front)");
            if (M != 5) ui::propI("Parallel runs (0 = auto)", &study_.workers, 0, hw * 2);
            ui::endProps();
        }
        if (act == 0 && !(M == 6 && study_.adjMode == 1)) ui::hint("Define and activate parameters in the Parameters tab first.");
        if ((M >= 3 || Response(study_.response) == Response::Misfit) && case_.observations.empty())
            ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT " History matching needs observed data (Observed tab).");
        if (M == 5) {
            ui::header("Parameters");
            bool ch = false;
            for (size_t i = 0; i < P.size(); ++i) {
                if (!P[i].active) continue;
                ImGui::PushID(int(i));
                ImGui::SetNextItemWidth(-fs * 9);
                if (paramSlider("##s", P[i], &study_.manualVals[i])) ch = true;
                ImGui::SameLine(); ImGui::TextUnformatted(paramLabel(P[i], case_).c_str());
                ImGui::PopID();
            }
            if (ch) { study_.manualPending = true; study_.manualEditTime = ImGui::GetTime(); }
            ImGui::Checkbox("Run automatically after a change", &study_.manualAuto);
            if (!study_.manualAuto) study_.manualPending = false;
            if (ImGui::Button(ICON_PLAY " Run now")) study_.deferred = [this]() { runManual(); };
            ImGui::SameLine();
            if (ImGui::Button(ICON_ROTATE_CCW " Reset to base")) { for (size_t i = 0; i < P.size(); ++i) study_.manualVals[i] = P[i].base; study_.manualPending = study_.manualAuto; }
            {
                std::lock_guard<std::mutex> l(study_.m);
                if (study_.hasManual) {
                    if (study_.manualRunning) ImGui::TextDisabled("running...");
                    else ImGui::TextColored(col4(Pal::Accent), "misfit %.4g", study_.manualMisfit);
                }
            }
            if (ImGui::Button(ICON_CHECK " Apply slider values to the case")) applyStudyValues(P, study_.manualVals, "Slider values");
        } else {
            ImGui::Spacing();
            char b[128];
            int runs = 0;
            if (M == 0) runs = 2 * act + 1;
            else if (M == 1) runs = std::max(4, study_.mcSamples);
            else if (M == 2) runs = int(designPoints(P, DesignKind(study_.designKind), study_.designN, unsigned(study_.seed)).size());
            else if (M == 4) runs = study_.es.ensemble * (study_.es.assimilations + 1);
            if (M == 3 && study_.hmOpt == 0) std::snprintf(b, sizeof(b), ICON_PLAY "  Run %s  (~%d runs / iteration)", kMethods[M], act + 4);
            else if (M == 3) std::snprintf(b, sizeof(b), ICON_PLAY "  Run %s, analytic Jacobian  (~4 runs / iteration)", kMethods[M]);
            else if (M == 6 && study_.adjMode == 0) std::snprintf(b, sizeof(b), ICON_PLAY "  Run adjoint L-BFGS  (~1-2 runs / iteration)");
            else if (M == 6) std::snprintf(b, sizeof(b), ICON_PLAY "  Run per-cell field match  (1 forward + 1 adjoint / iteration)");
            else std::snprintf(b, sizeof(b), ICON_PLAY "  Run %s  (%d runs)", kMethods[M], runs);
            const bool fieldMode = M == 6 && study_.adjMode == 1;
            ImGui::BeginDisabled(study_.running || (fieldMode ? (case_.observations.empty() || (!study_.fieldOpt.perm && !study_.fieldOpt.poro)) : act == 0));
            if (fieldMode && ui::accentButton(b, ImVec2(-FLT_MIN, fs * 2.2f))) {
                const Case base = case_; const FieldHMOptions o = study_.fieldOpt;
                study_.deferred = [this, base, o]() {
                    startStudy("Per-cell field history match (adjoint)", [base, o](StudyState& S, const StudyCallbacks& cb) {
                        FieldHMResult r = runAdjointFieldHM(base, o, cb);
                        if (!r.ok && cb.log) cb.log("Field HM failed: " + r.error);
                        std::lock_guard<std::mutex> l(S.m); S.field = std::move(r); S.hasField = true; S.adjLast = 1; S.adjIter = -1; S.centerTab = 11;
                    });
                };
            }
            if (!fieldMode)
            if (ui::accentButton(b, ImVec2(-FLT_MIN, fs * 2.2f))) {
                const Case base = case_; const Response r = Response(study_.response); const Economics eco = study_.eco;
                const int n = std::max(4, study_.mcSamples); const unsigned seed = unsigned(study_.seed);
                const LMOptions lo = study_.lm; const ESMDAOptions eo = study_.es; const AdjointHMOptions ao = study_.adj; const int hopt = M == 6 ? 2 : study_.hmOpt;
                const DesignKind dk = DesignKind(study_.designKind); const int dn = study_.designN; const bool quad = study_.proxyQuad, mx = study_.proxyMax; const int pkind = study_.proxyKind;
                study_.deferred = [this, M, base, r, eco, n, seed, lo, eo, ao, hopt, dk, dn, quad, mx, pkind]() {
                    if (M == 0) startStudy("Tornado sensitivity", [base, r, eco](StudyState& S, const StudyCallbacks& cb) {
                        TornadoResult t = runTornado(base, base.studyParams, r, eco, cb);
                        std::lock_guard<std::mutex> l(S.m); S.tornado = std::move(t); S.tornadoResp = r; S.hasTornado = true; S.centerTab = 3;
                    });
                    else if (M == 1) startStudy("Monte Carlo uncertainty", [base, r, eco, n, seed](StudyState& S, const StudyCallbacks& cb) {
                        MonteCarloResult m = runMonteCarlo(base, base.studyParams, n, seed, r, eco, cb);
                        std::lock_guard<std::mutex> l(S.m); S.mc = std::move(m); S.mcResp = r; S.hasMC = true; S.centerTab = 4;
                    });
                    else if (M == 2) startStudy("Experimental design + proxy", [base, r, eco, seed, dk, dn, quad, mx, pkind](StudyState& S, const StudyCallbacks& cb) {
                        const auto X = designPoints(base.studyParams, dk, dn, seed);
                        const auto runs = runDesign(base, base.studyParams, X, "design", cb);
                        const auto sig = autoSigma(base.observations);
                        std::vector<std::vector<double>> xs; std::vector<double> y;
                        for (size_t i = 0; i < runs.size() && i < X.size(); ++i) {
                            if (!runs[i].ok) continue;
                            const double v = responseValue(runs[i], base, r, eco, base.observations, sig);
                            if (!std::isfinite(v)) continue;
                            xs.push_back(X[i]); y.push_back(v);
                        }
                        ProxyModel pm = pkind == 1 ? fitProxyGP(base.studyParams, xs, y) : fitProxy(base.studyParams, xs, y, quad);
                        double pred = NAN;
                        std::vector<double> opt = pm.ok ? proxyOptimum(pm, base.studyParams, mx, seed, &pred) : std::vector<double>();
                        if (cb.log) { char b[160]; std::snprintf(b, sizeof(b), "proxy: %d samples, R2 %.4f, Q2 %.4f, RMSE %.4g", pm.samples, pm.r2, pm.q2, pm.rmse); cb.log(b); }
                        std::lock_guard<std::mutex> l(S.m);
                        S.proxy = std::move(pm); S.designX = std::move(xs); S.designY = std::move(y); S.proxyOpt = std::move(opt); S.proxyPred = pred;
                        S.proxyActual = NAN; S.proxyParams = base.studyParams; S.proxyResp = int(r); S.hasProxy = true; S.centerTab = 6;
                    });
                    else startStudy(M == 4 ? "ES-MDA history match" : hopt == 2 ? "Adjoint L-BFGS history match" : "Levenberg-Marquardt history match",
                                    [base, M, lo, eo, ao, hopt](StudyState& S, const StudyCallbacks& cb) {
                        HMResult h = M == 4 ? runESMDA(base, base.studyParams, eo, cb)
                                   : hopt == 2 ? runAdjointLBFGS(base, base.studyParams, ao, cb) : runLevenbergMarquardt(base, base.studyParams, lo, cb);
                        std::lock_guard<std::mutex> l(S.m); S.hm = std::move(h); S.hmMethod = M == 4 ? 1 : 0; S.hmAlgo = M == 4 ? 0 : hopt;
                        S.hmParams = base.studyParams; S.hasHM = true; S.centerTab = hopt == 2 ? 11 : 9;
                        if (hopt == 2) { S.adjLast = 0; S.adjIter = -1; }
                    });
                };
            }
            ImGui::EndDisabled();
            ImGui::TextDisabled("Every run is archived (right pane) and shown in the views.");
            if (M == 6 && study_.adjMode == 0 && study_.hasHM && study_.hmAlgo == 2 && ImGui::Button(ICON_CHART_LINE " Open the Adjoint HM view")) study_.centerTab = 11;
            if (M == 3 || (M == 6 && study_.adjMode == 0)) gradientCheckPanel(act);
            if (fieldMode) {
                ui::hint("Matches a log-multiplier of the permeability (all directions) and / or porosity of EVERY active cell "
                         "against the observed data. One forward run + one adjoint pass per iteration gives all gradients "
                         "(analytic dR/dtheta), so the cost does not grow with the number of cells. A smoothness term keeps "
                         "the multiplier field geologically plausible.");
                std::lock_guard<std::mutex> l(study_.m);
                if (study_.hasField) {
                    const FieldHMResult& F = study_.field;
                    if (!F.ok) ImGui::TextColored(col4(Pal::Err), "%s", F.error.c_str());
                    else {
                        ImGui::Text("misfit %.4g -> %.4g  (%d iterations, %d runs; forward %.1fs, adjoint %.1fs)", F.initialMisfit, F.bestMisfit,
                                    F.iterations, F.runs, F.secForward, F.secAdjoint);
                        if (ImGui::Button(ICON_CHART_LINE " Open the Adjoint HM view")) study_.centerTab = 11;
                        if (F.history.size() > 1 && ImPlot::BeginPlot("##fhm", ImVec2(-1, fs * 9))) {
                            ImPlot::SetupAxes("iteration", "misfit", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                            ImPlot::PlotLine("misfit", F.history.data(), int(F.history.size()));
                            ImPlot::EndPlot();
                        }
                        double kmin = 1e300, kmax = 0;
                        for (double v : F.multK) { kmin = std::min(kmin, v); kmax = std::max(kmax, v); }
                        if (!F.multK.empty()) ImGui::TextDisabled("permeability multipliers %.3g .. %.3g", kmin, kmax);
                        if (ImGui::Button(ICON_CHECK " Apply the matched field to the case")) {
                            const Case best = F.best;
                            study_.deferred = [this, best]() {
                                case_.grid.permx = best.grid.permx; case_.grid.permy = best.grid.permy;
                                case_.grid.permz = best.grid.permz; case_.grid.poro = best.grid.poro;
                                markDirty(true);
                                notify("Matched per-cell permeability / porosity applied to the case - re-run the simulation (F5)");
                            };
                        }
                    }
                }
            }
        }
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    if (ui::tabItem(ICON_SLIDERS_HORIZONTAL " Parameters")) {
        ImGui::BeginChild("##par", ImVec2(0, 0), 0, ImGuiWindowFlags_HorizontalScrollbar);
        studyParams();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    if (study_.wsMode == 0 && ui::tabItem(ICON_TARGET " Observed")) {
        ImGui::BeginChild("##obs");
        std::unique_lock<std::mutex> lk(study_.m);
        studyObservations();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ---------------------------------------------------------------------------- centre: views
void App::studiesCenter() {
    const float fs = ImGui::GetFontSize();
    auto& A = study_.archive;
    const auto& AP = study_.archParams;
    // run filter
    std::vector<std::string> tags;
    for (auto& a : A) if (std::find(tags.begin(), tags.end(), a.tag) == tags.end()) tags.push_back(a.tag);
    auto inFilter = [&](const StudyState::ArchRun& a) {
        if (study_.ovFilter == 0) return true;
        if (study_.ovFilter == 1) return a.id >= study_.studyFirstId;
        const int t = study_.ovFilter - 2;
        return t < int(tags.size()) && a.tag == tags[size_t(t)];
    };
    int best = -1;
    for (int i = 0; i < int(A.size()); ++i)
        if (A[size_t(i)].ok && inFilter(A[size_t(i)]) && std::isfinite(A[size_t(i)].misfit) && (best < 0 || A[size_t(i)].misfit < A[size_t(best)].misfit)) best = i;
    int sel = -1;
    for (int i = 0; i < int(A.size()); ++i) if (A[size_t(i)].id == study_.selRun) sel = i;

    const int want = study_.centerTab; study_.centerTab = -1;
    // a study result asks for a view of the other mode -> switch the mode
    if (want == 3 || want == 4 || want == 6) study_.wsMode = 0;
    if (want == 2 || want == 5 || want >= 8) study_.wsMode = 1;
    const bool hmMode = study_.wsMode == 1;
    auto tabFlags = [&](int t) { return want == t ? ImGuiTabItemFlags_SetSelected : 0; };
    if (!ImGui::BeginTabBar("##views", ImGuiTabBarFlags_FittingPolicyScroll)) return;

    auto filterCombo = [&]() {
        std::vector<std::string> names = {"All runs", "Last study"};
        for (auto& t : tags) names.push_back("Tag: " + t);
        ImGui::SetNextItemWidth(fs * 9);
        study_.ovFilter = std::clamp(study_.ovFilter, 0, int(names.size()) - 1);
        if (ImGui::BeginCombo("##filter", names[size_t(study_.ovFilter)].c_str())) {
            for (int i = 0; i < int(names.size()); ++i) if (ImGui::Selectable(names[size_t(i)].c_str(), study_.ovFilter == i)) study_.ovFilter = i;
            ImGui::EndCombo();
        }
    };

    // ---- history matching: observed data manager, match quality, parameters
    if (hmMode && ui::tabItem(ICON_TARGET " Observed data", nullptr, tabFlags(8))) {
        ImGui::BeginChild("##hmobs");
        studyObservations();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    if (hmMode && ui::tabItem(ICON_GRID_3X3 " Match quality", nullptr, tabFlags(9))) {
        ImGui::BeginChild("##hmmq");
        hmMatchQuality();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    if (hmMode && ui::tabItem(ICON_SLIDERS_HORIZONTAL " Parameters", nullptr, tabFlags(10))) {
        ImGui::BeginChild("##hmpar");
        hmParameters();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }

    if (hmMode && ui::tabItem(ICON_ACTIVITY " Adjoint HM", nullptr, tabFlags(11))) {
        ImGui::BeginChild("##hmadj");
        adjointHMView();
        ImGui::EndChild();
        ImGui::EndTabItem();
    }

    // ---- profiles overlay
    if (ui::tabItem(ICON_CHART_LINE " Profiles", nullptr, tabFlags(0))) {
        ImGui::SetNextItemWidth(fs * 6);
        ImGui::Combo("##ovv", &study_.ovVec, kOvVec, kNumOv, 20);
        const bool wellVec = study_.ovVec >= kOvFirstWell;
        std::string well;
        if (wellVec) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 8);
            study_.ovWell = std::clamp(study_.ovWell, 0, std::max(0, int(case_.wells.size()) - 1));
            if (ImGui::BeginCombo("##ovw", case_.wells.empty() ? "-" : case_.wells[size_t(study_.ovWell)].name.c_str())) {
                for (int w = 0; w < int(case_.wells.size()); ++w) if (ImGui::Selectable(case_.wells[size_t(w)].name.c_str(), study_.ovWell == w)) study_.ovWell = w;
                ImGui::EndCombo();
            }
            if (!case_.wells.empty()) well = case_.wells[size_t(study_.ovWell)].name;
        }
        auto wrap = [&](const char* lab) {
            const float w = ImGui::CalcTextSize(lab).x + ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x * 2;
            ImGui::SameLine();
            if (ImGui::GetContentRegionAvail().x < w) ImGui::NewLine();
        };
        wrap("filter combo  xxxxxxx"); filterCombo();
        wrap("All runs"); ImGui::Checkbox("All runs", &study_.ovArchive);
        wrap("P10-P90 band"); ImGui::Checkbox("P10-P90 band", &study_.ovBand);
        wrap("Observed"); ImGui::Checkbox("Observed", &study_.ovObs);
        // cached curves
        struct Cache { unsigned ver = ~0u; int vec = -1, filter = -1; std::string well; std::vector<std::vector<double>> x, y; std::vector<int> run;
                       std::vector<double> bt, b10, b50, b90; };
        static Cache C;
        if (C.ver != study_.archiveVersion || C.vec != study_.ovVec || C.well != well || C.filter != study_.ovFilter) {
            C = Cache(); C.ver = study_.archiveVersion; C.vec = study_.ovVec; C.well = well; C.filter = study_.ovFilter;
            double tmax = 0;
            for (int i = 0; i < int(A.size()); ++i) {
                const auto& a = A[size_t(i)];
                if (!a.ok || !inFilter(a) || a.run.summary.size() < 2) continue;
                std::vector<double> x, y;
                for (size_t k = ovIsRate(study_.ovVec) ? 1 : 0; k < a.run.summary.size(); ++k) {
                    const double v = ovValue(a.run, case_, study_.ovVec, well, a.run.summary[k]);
                    if (!std::isfinite(v)) continue;
                    x.push_back(a.run.summary[k].time); y.push_back(v);
                }
                if (x.size() < 2) continue;
                tmax = std::max(tmax, x.back());
                C.x.push_back(std::move(x)); C.y.push_back(std::move(y)); C.run.push_back(i);
            }
            if (C.x.size() >= 3 && tmax > 0) {
                for (int k = 0; k <= 100; ++k) {
                    const double t = tmax * k / 100.0;
                    std::vector<double> v;
                    for (size_t r = 0; r < C.x.size(); ++r) {
                        const auto& X = C.x[r]; const auto& Y = C.y[r];
                        if (t < X.front() || t > X.back()) continue;
                        const size_t j = size_t(std::lower_bound(X.begin(), X.end(), t) - X.begin());
                        if (j == 0) { v.push_back(Y[0]); continue; }
                        const double f = (t - X[j - 1]) / std::max(1e-12, X[j] - X[j - 1]);
                        v.push_back(Y[j - 1] + f * (Y[j] - Y[j - 1]));
                    }
                    if (v.size() < 3) continue;
                    C.bt.push_back(t); C.b10.push_back(pct(v, 0.1)); C.b50.push_back(pct(v, 0.5)); C.b90.push_back(pct(v, 0.9));
                }
            }
        }
        if (ImPlot::BeginPlot("##profiles", ImVec2(-1, -1))) {
            ImPlot::SetupAxes("Time (days)", kOvVec[study_.ovVec], ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupLegend(ImPlotLocation_NorthWest);
            if (study_.ovBand && !C.bt.empty()) {
                ImPlotSpec sb; sb.FillColor = ImVec4(0.39f, 0.6f, 1.0f, 0.22f);
                ImPlot::PlotShaded("P10-P90", C.bt.data(), C.b10.data(), C.b90.data(), int(C.bt.size()), sb);
                ImPlotSpec s5; s5.LineColor = ImVec4(0.45f, 0.65f, 1.0f, 0.9f); s5.LineWeight = 1.5f;
                ImPlot::PlotLine("P50", C.bt.data(), C.b50.data(), int(C.bt.size()), s5);
            }
            if (study_.ovArchive) {
                ImPlotSpec sg; sg.LineColor = ImVec4(0.62f, 0.64f, 0.7f, std::clamp(6.0f / std::max<size_t>(1, C.x.size()), 0.12f, 0.55f)); sg.LineWeight = 1.0f;
                for (size_t r = 0; r < C.x.size(); ++r) {
                    if (C.run[r] == best || C.run[r] == sel) continue;
                    ImPlot::PlotLine("archived runs", C.x[r].data(), C.y[r].data(), int(C.x[r].size()), sg);
                }
            }
            // base case: last regular simulation of the model (dashed)
            {
                std::lock_guard<std::mutex> l(runner_.mutex());
                const RunResults& R = runner_.results();
                if (R.valid && R.summary.size() > 1 && !runner_.running()) {
                    RunResult rr; rr.ok = true; rr.summary = R.summary;
                    std::vector<double> x, y;
                    for (size_t k = ovIsRate(study_.ovVec) ? 1 : 0; k < R.summary.size(); ++k) {
                        const double v = ovValue(rr, R.input, study_.ovVec, well, R.summary[k]);
                        if (std::isfinite(v)) { x.push_back(R.summary[k].time); y.push_back(v); }
                    }
                    ImPlotSpec sd; sd.LineColor = ImVec4(1, 1, 1, 0.85f); sd.LineWeight = 1.6f; sd.Flags = ImPlotLineFlags_Segments;
                    // dashed look: draw every other segment
                    std::vector<double> dx, dy;
                    for (size_t k = 1; k < x.size(); ++k) if (k % 2) { dx.push_back(x[k - 1]); dy.push_back(y[k - 1]); dx.push_back(x[k]); dy.push_back(y[k]); }
                    if (x.size() < 12) { sd.Flags = 0; dx = x; dy = y; }
                    if (!dx.empty()) ImPlot::PlotLine("base case (last run)", dx.data(), dy.data(), int(dx.size()), sd);
                }
            }
            auto drawRun = [&](int i, const char* lab, ImVec4 c, float w) {
                for (size_t r = 0; r < C.x.size(); ++r) if (C.run[r] == i) {
                    ImPlotSpec s; s.LineColor = c; s.LineWeight = w;
                    ImPlot::PlotLine(lab, C.x[r].data(), C.y[r].data(), int(C.x[r].size()), s);
                }
            };
            if (best >= 0) { char b[64]; std::snprintf(b, sizeof(b), "best match #%d", A[size_t(best)].id); drawRun(best, b, col4(Pal::Accent), 3.0f); }
            if (sel >= 0 && sel != best) { char b[64]; std::snprintf(b, sizeof(b), "selected #%d", A[size_t(sel)].id); drawRun(sel, b, ImVec4(1.0f, 0.85f, 0.25f, 1), 2.5f); }
            // observed data
            const int q = ovQty(study_.ovVec);
            if (study_.ovObs && q >= 0) {
                std::vector<double> ox, oy, oe;
                const auto sig = autoSigma(case_.observations);
                for (size_t k = 0; k < case_.observations.size(); ++k) {
                    const Observation& o = case_.observations[k];
                    if (int(o.qty) != q || (wellVec && o.well != well)) continue;
                    ox.push_back(o.time); oy.push_back(o.value); oe.push_back(k < sig.size() && sig[k] < 1e20 ? sig[k] : 0.0);
                }
                if (!ox.empty()) {
                    ImPlotSpec se; se.LineColor = ImVec4(1, 0.45f, 0.4f, 0.6f);
                    ImPlot::PlotErrorBars("observed", ox.data(), oy.data(), oe.data(), int(ox.size()), se);
                    ImPlotSpec so; so.Marker = ImPlotMarker_Circle; so.MarkerSize = 3.5f; so.MarkerFillColor = ImVec4(1, 0.4f, 0.35f, 1); so.MarkerLineColor = ImVec4(1, 0.4f, 0.35f, 1);
                    ImPlot::PlotScatter("observed", ox.data(), oy.data(), int(ox.size()), so);
                }
            }
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ---- cross-plot: parameter vs misfit / response / objective term
    if (ui::tabItem(ICON_CHART_SCATTER " Cross-plot", nullptr, tabFlags(1))) {
        std::vector<int> pidx;
        for (int i = 0; i < int(AP.size()); ++i) if (AP[size_t(i)].active) pidx.push_back(i);
        if (A.empty() || pidx.empty()) ui::hint("Run a study first: every run becomes a point (parameter value vs misfit / response / objective term).");
        else {
            std::vector<std::string> xs, ys = {"Misfit (total)", "Response"};
            for (int i : pidx) xs.push_back(paramLabel(AP[size_t(i)], case_));
            for (auto& t : study_.termNames) ys.push_back("Term: " + t);
            study_.cpX = std::clamp(study_.cpX, 0, int(xs.size()) - 1);
            study_.cpY = std::clamp(study_.cpY, 0, int(ys.size()) - 1);
            ImGui::SetNextItemWidth(fs * 12);
            if (ImGui::BeginCombo("X", xs[size_t(study_.cpX)].c_str())) { for (int i = 0; i < int(xs.size()); ++i) if (ImGui::Selectable(xs[size_t(i)].c_str(), i == study_.cpX)) study_.cpX = i; ImGui::EndCombo(); }
            ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 12);
            if (ImGui::BeginCombo("Y", ys[size_t(study_.cpY)].c_str())) { for (int i = 0; i < int(ys.size()); ++i) if (ImGui::Selectable(ys[size_t(i)].c_str(), i == study_.cpY)) study_.cpY = i; ImGui::EndCombo(); }
            ImGui::SameLine(); filterCombo();
            const int pi = pidx[size_t(study_.cpX)];
            auto yOf = [&](const StudyState::ArchRun& a) {
                if (study_.cpY == 0) return a.misfit;
                if (study_.cpY == 1) return a.response;
                const size_t t = size_t(study_.cpY - 2);
                return t < a.terms.size() ? a.terms[t] : NAN;
            };
            if (ImPlot::BeginPlot("##cross", ImVec2(-1, -1))) {
                ImPlot::SetupAxes(xs[size_t(study_.cpX)].c_str(), ys[size_t(study_.cpY)].c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (AP[size_t(pi)].log) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                if (study_.cpY != 1) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                for (auto& t : tags) {
                    std::vector<double> x, y;
                    for (auto& a : A) if (a.tag == t && a.ok && inFilter(a) && size_t(pi) < a.values.size()) {
                        const double v = yOf(a);
                        if (std::isfinite(v) && (study_.cpY == 1 || v > 0)) { x.push_back(a.values[size_t(pi)]); y.push_back(v); }
                    }
                    if (x.empty()) continue;
                    ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 4; s.MarkerFillColor = col4(tagColor(t)); s.MarkerLineColor = col4(tagColor(t));
                    ImPlot::PlotScatter(t.c_str(), x.data(), y.data(), int(x.size()), s);
                }
                auto mark = [&](int i, const char* lab, ImVec4 c) {
                    if (i < 0 || size_t(pi) >= A[size_t(i)].values.size()) return;
                    double x = A[size_t(i)].values[size_t(pi)], y = yOf(A[size_t(i)]);
                    if (!std::isfinite(y)) return;
                    ImPlotSpec s; s.Marker = ImPlotMarker_Diamond; s.MarkerSize = 9; s.MarkerFillColor = c; s.MarkerLineColor = ImVec4(1, 1, 1, 1);
                    ImPlot::PlotScatter(lab, &x, &y, 1, s);
                };
                mark(best, "best match", col4(Pal::Accent));
                mark(sel, "selected", ImVec4(1.0f, 0.85f, 0.25f, 1));
                // click: select the nearest run
                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
                    const ImVec2 m = ImGui::GetMousePos();
                    double bd = 1e30; int bi = -1;
                    for (int i = 0; i < int(A.size()); ++i) {
                        const auto& a = A[size_t(i)];
                        if (!a.ok || !inFilter(a) || size_t(pi) >= a.values.size()) continue;
                        const double y = yOf(a);
                        if (!std::isfinite(y)) continue;
                        const ImVec2 p = ImPlot::PlotToPixels(a.values[size_t(pi)], y);
                        const double d = (p.x - m.x) * (p.x - m.x) + (p.y - m.y) * (p.y - m.y);
                        if (d < bd) { bd = d; bi = i; }
                    }
                    if (bi >= 0 && bd < 400) study_.selRun = A[size_t(bi)].id;
                }
                ImPlot::EndPlot();
            }
        }
        ImGui::EndTabItem();
    }

    // ---- Pareto front of two objective terms
    if (hmMode && ui::tabItem(ICON_TROPHY " Pareto", nullptr, tabFlags(2))) {
        const auto& T = study_.termNames;
        if (T.size() < 2) ui::hint("The Pareto front needs at least two objective terms: add observations of several quantities / wells and pick an objective grouping (Method tab), then run a study.");
        else {
            study_.parA = std::clamp(study_.parA, 0, int(T.size()) - 1);
            study_.parB = std::clamp(study_.parB, 0, int(T.size()) - 1);
            ImGui::SetNextItemWidth(fs * 10);
            if (ImGui::BeginCombo("Objective A", T[size_t(study_.parA)].c_str())) { for (int i = 0; i < int(T.size()); ++i) if (ImGui::Selectable(T[size_t(i)].c_str(), i == study_.parA)) study_.parA = i; ImGui::EndCombo(); }
            ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 10);
            if (ImGui::BeginCombo("Objective B", T[size_t(study_.parB)].c_str())) { for (int i = 0; i < int(T.size()); ++i) if (ImGui::Selectable(T[size_t(i)].c_str(), i == study_.parB)) study_.parB = i; ImGui::EndCombo(); }
            ImGui::SameLine(); filterCombo();
            std::vector<std::vector<double>> F; std::vector<int> who;
            for (int i = 0; i < int(A.size()); ++i) {
                const auto& a = A[size_t(i)];
                if (!a.ok || !inFilter(a) || a.terms.size() <= size_t(std::max(study_.parA, study_.parB))) continue;
                F.push_back({a.terms[size_t(study_.parA)], a.terms[size_t(study_.parB)]}); who.push_back(i);
            }
            std::vector<int> front = paretoFront(F);
            std::sort(front.begin(), front.end(), [&](int x, int y) { return F[size_t(x)][0] < F[size_t(y)][0]; });
            if (ImPlot::BeginPlot("##pareto", ImVec2(-1, -fs * 6))) {
                ImPlot::SetupAxes(T[size_t(study_.parA)].c_str(), T[size_t(study_.parB)].c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10); ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                std::vector<double> x, y;
                for (auto& f : F) { x.push_back(std::max(1e-12, f[0])); y.push_back(std::max(1e-12, f[1])); }
                ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 3.5f; s.MarkerFillColor = ImVec4(0.6f, 0.62f, 0.7f, 0.7f); s.MarkerLineColor = s.MarkerFillColor;
                if (!x.empty()) ImPlot::PlotScatter("runs", x.data(), y.data(), int(x.size()), s);
                std::vector<double> fx, fy;
                for (int k : front) { fx.push_back(x[size_t(k)]); fy.push_back(y[size_t(k)]); }
                ImPlotSpec sf; sf.LineColor = col4(Pal::Accent); sf.LineWeight = 2; sf.Marker = ImPlotMarker_Diamond; sf.MarkerSize = 6; sf.MarkerFillColor = col4(Pal::Accent);
                if (!fx.empty()) ImPlot::PlotLine("Pareto front", fx.data(), fy.data(), int(fx.size()), sf);
                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) {
                    const ImVec2 m = ImGui::GetMousePos(); double bd = 1e30; int bi = -1;
                    for (size_t k = 0; k < x.size(); ++k) { const ImVec2 p = ImPlot::PlotToPixels(x[k], y[k]); const double d = (p.x - m.x) * (p.x - m.x) + (p.y - m.y) * (p.y - m.y); if (d < bd) { bd = d; bi = int(k); } }
                    if (bi >= 0 && bd < 400) study_.selRun = A[size_t(who[size_t(bi)])].id;
                }
                ImPlot::EndPlot();
            }
            ImGui::TextDisabled("%d non-dominated runs of %d (click a point to select the run):", int(front.size()), int(F.size()));
            std::string ids;
            for (int k : front) { ids += "#" + std::to_string(A[size_t(who[size_t(k)])].id) + " "; if (ids.size() > 200) { ids += "..."; break; } }
            ImGui::TextWrapped("%s", ids.c_str());
        }
        ImGui::EndTabItem();
    }

    // ---- tornado
    if (!hmMode && ui::tabItem(ICON_CHART_COLUMN " Tornado", nullptr, tabFlags(3))) {
        if (!study_.hasTornado) ui::hint("Run the tornado method to rank the parameters by their one-at-a-time swing.");
        else {
            const TornadoResult& T = study_.tornado;
            std::vector<std::string> lab; std::vector<double> lo, hi;
            for (auto& bar : T.bars) {
                const bool okp = bar.param < int(AP.size());
                lab.push_back(okp ? paramLabel(AP[size_t(bar.param)], case_) + "  [" + fmtv(AP[size_t(bar.param)].low) + " .. " + fmtv(AP[size_t(bar.param)].high) + "]" : "?");
                lo.push_back(bar.okLow ? bar.low : NAN); hi.push_back(bar.okHigh ? bar.high : NAN);
            }
            ImGui::TextDisabled("%s  -  base %.5g", kResp[int(study_.tornadoResp)], T.base);
            barChart("##wstornado", lab, lo, hi, T.base, "parameter at low", "parameter at high", false);
        }
        ImGui::EndTabItem();
    }

    // ---- Monte Carlo
    if (!hmMode && ui::tabItem(ICON_SHUFFLE " Monte Carlo", nullptr, tabFlags(4))) {
        if (!study_.hasMC) ui::hint("Run the Monte Carlo method: Latin-hypercube samples of all active parameters, P10 / P50 / P90 and rank correlations.");
        else {
            const MonteCarloResult& M = study_.mc;
            ImGui::Text("%s:  P10 %.5g   P50 %.5g   P90 %.5g   mean %.5g   (%d samples, %d failed)", kResp[int(study_.mcResp)], M.p10, M.p50, M.p90, M.mean,
                        int(M.response.size()), M.failed);
            const float h = (ImGui::GetContentRegionAvail().y - fs) * 0.5f;
            if (ui::beginTable("##mc", 2, ImGuiTableFlags_Resizable)) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                std::vector<double> v; for (double x : M.response) if (std::isfinite(x)) v.push_back(x);
                if (!v.empty() && ImPlot::BeginPlot("Distribution", ImVec2(-1, h))) {
                    ImPlot::SetupAxes(kResp[int(study_.mcResp)], "count", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec s; s.FillColor = ImVec4(0.4f, 0.55f, 0.95f, 0.75f);
                    ImPlot::PlotHistogram("samples", v.data(), int(v.size()), std::clamp(int(std::sqrt(double(v.size()))) + 2, 6, 40), 1.0, ImPlotRange(), s);
                    const double px[3] = {M.p10, M.p50, M.p90};
                    ImPlotSpec sl; sl.LineColor = ImVec4(1, 0.8f, 0.3f, 1); sl.LineWeight = 2;
                    ImPlot::PlotInfLines("P10 / P50 / P90", px, 3, sl);
                    ImPlot::EndPlot();
                }
                ImGui::TableSetColumnIndex(1);
                if (!M.times.empty() && ImPlot::BeginPlot("FOPT fan chart", ImVec2(-1, h))) {
                    ImPlot::SetupAxes("Time (days)", "FOPT", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec sb; sb.FillColor = ImVec4(0.39f, 0.6f, 1.0f, 0.25f);
                    ImPlot::PlotShaded("P10-P90", M.times.data(), M.fanP10.data(), M.fanP90.data(), int(M.times.size()), sb);
                    ImPlotSpec s5; s5.LineColor = col4(Pal::Accent); s5.LineWeight = 2;
                    ImPlot::PlotLine("P50", M.times.data(), M.fanP50.data(), int(M.times.size()), s5);
                    ImPlot::EndPlot();
                }
                ui::endTable();
            }
            std::vector<std::string> lab; std::vector<double> lo, hi;
            std::vector<int> order;
            for (int i = 0; i < int(M.spearman.size()); ++i) if (i < int(AP.size()) && AP[size_t(i)].active && std::isfinite(M.spearman[size_t(i)])) order.push_back(i);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return std::fabs(M.spearman[size_t(a)]) > std::fabs(M.spearman[size_t(b)]); });
            for (int i : order) { lab.push_back(paramLabel(AP[size_t(i)], case_)); const double r = M.spearman[size_t(i)]; lo.push_back(r < 0 ? r : NAN); hi.push_back(r >= 0 ? r : NAN); }
            ui::header("Spearman rank correlation with the response");
            barChart("##wsspear", lab, lo, hi, 0.0, "negative", "positive", true);
        }
        ImGui::EndTabItem();
    }

    // ---- convergence
    if (hmMode && ui::tabItem(ICON_ACTIVITY " Convergence", nullptr, tabFlags(5))) {
        std::vector<double> id, mf, bestSoFar;
        double b = 1e300;
        for (auto& a : A) if (a.ok && std::isfinite(a.misfit) && inFilter(a)) { id.push_back(a.id); mf.push_back(std::max(1e-12, a.misfit)); b = std::min(b, a.misfit); bestSoFar.push_back(std::max(1e-12, b)); }
        filterCombo();
        const float h = study_.hasHM && !study_.hm.history.empty() ? (ImGui::GetContentRegionAvail().y - fs) * 0.5f : ImGui::GetContentRegionAvail().y;
        if (id.empty()) ui::hint("No runs with a misfit yet (needs observed data).");
        else if (ImPlot::BeginPlot("Misfit of every run", ImVec2(-1, h))) {
            ImPlot::SetupAxes("run #", "misfit", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 3; s.MarkerFillColor = ImVec4(0.6f, 0.62f, 0.7f, 0.8f); s.MarkerLineColor = s.MarkerFillColor;
            ImPlot::PlotScatter("runs", id.data(), mf.data(), int(id.size()), s);
            ImPlotSpec sb; sb.LineColor = col4(Pal::Accent); sb.LineWeight = 2.5f; sb.Flags = ImPlotLineFlags_None;
            ImPlot::PlotStairs("best so far", id.data(), bestSoFar.data(), int(id.size()), sb);
            ImPlot::EndPlot();
        }
        if (study_.hasHM && !study_.hm.history.empty() && ImPlot::BeginPlot(study_.hmMethod == 0 ? (study_.hmAlgo == 2 ? "Adjoint L-BFGS iterations" : "Levenberg-Marquardt iterations") : "ES-MDA assimilations", ImVec2(-1, -1))) {
            ImPlot::SetupAxes(study_.hmMethod == 0 ? "iteration" : "assimilation", "misfit", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            std::vector<double> x, y, sp;
            for (auto& it : study_.hm.history) { x.push_back(it.iter); y.push_back(std::max(1e-12, it.misfit)); sp.push_back(it.spread); }
            ImPlotSpec s; s.LineColor = col4(Pal::Accent); s.LineWeight = 2.5f; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 4; s.MarkerFillColor = col4(Pal::Accent);
            ImPlot::PlotLine(study_.hmMethod == 0 ? "misfit" : "ensemble mean misfit", x.data(), y.data(), int(x.size()), s);
            ImPlot::EndPlot();
        }
        ImGui::EndTabItem();
    }

    // ---- proxy
    if (!hmMode && ui::tabItem(ICON_SPARKLES " Proxy", nullptr, tabFlags(6))) {
        if (!study_.hasProxy) ui::hint("Run 'Experimental design + proxy': training runs from a design, a polynomial response surface, its quality (R2, leave-one-out Q2) and the optimum.");
        else {
            const ProxyModel& PM = study_.proxy;
            const auto& PP = study_.proxyParams;
            if (!PM.ok) ImGui::TextColored(col4(Pal::Err), "The proxy could not be fitted (too few successful runs for the number of parameters).");
            else {
                const ImU32 q = PM.q2 > 0.9 ? Pal::Ok : PM.q2 > 0.6 ? Pal::Warn : Pal::Err;
                ImGui::Text("%s proxy of %s from %d runs:", PM.kind == 1 ? "Gaussian-process" : PM.quadratic ? "Quadratic" : "Linear", kResp[study_.proxyResp], PM.samples);
                ImGui::SameLine(); ImGui::TextColored(col4(q), "R2 %.4f   Q2 %.4f", PM.r2, PM.q2);
                ImGui::SameLine(); ImGui::TextDisabled("RMSE %.4g", PM.rmse);
                const float h = ImGui::GetContentRegionAvail().y * 0.55f;
                if (ui::beginTable("##px", 2, ImGuiTableFlags_Resizable)) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    if (ImPlot::BeginPlot("Predicted vs simulated", ImVec2(-1, h), ImPlotFlags_Equal)) {
                        ImPlot::SetupAxes("simulated", "proxy", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        std::vector<double> x, y;
                        double lo = 1e300, hi = -1e300;
                        for (size_t i = 0; i < study_.designX.size() && i < study_.designY.size(); ++i) {
                            x.push_back(study_.designY[i]); y.push_back(proxyPredict(PM, PP, study_.designX[i]));
                            lo = std::min({lo, x.back(), y.back()}); hi = std::max({hi, x.back(), y.back()});
                        }
                        if (!x.empty()) {
                            const double lx[2] = {lo, hi};
                            ImPlotSpec sl; sl.LineColor = ImVec4(1, 1, 1, 0.5f);
                            ImPlot::PlotLine("1:1", lx, lx, 2, sl);
                            ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 4; s.MarkerFillColor = col4(tagColor("design")); s.MarkerLineColor = s.MarkerFillColor;
                            ImPlot::PlotScatter("training runs", x.data(), y.data(), int(x.size()), s);
                        }
                        if (std::isfinite(study_.proxyActual) && std::isfinite(study_.proxyPred)) {
                            double ax = study_.proxyActual, ay = study_.proxyPred;
                            ImPlotSpec s; s.Marker = ImPlotMarker_Diamond; s.MarkerSize = 9; s.MarkerFillColor = col4(Pal::Accent); s.MarkerLineColor = ImVec4(1, 1, 1, 1);
                            ImPlot::PlotScatter("optimum (validated)", &ax, &ay, 1, s);
                        }
                        ImPlot::EndPlot();
                    }
                    ImGui::TableSetColumnIndex(1);
                    ui::header("Parameter importance");
                    std::vector<std::string> lab; std::vector<double> lo, hi;
                    std::vector<int> ord(PM.act.size()); std::iota(ord.begin(), ord.end(), 0);
                    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return PM.importance[size_t(a)] > PM.importance[size_t(b)]; });
                    for (int k : ord) { lab.push_back(paramLabel(PP[size_t(PM.act[size_t(k)])], case_)); lo.push_back(NAN); hi.push_back(PM.importance[size_t(k)]); }
                    barChart("##imp", lab, lo, hi, 0.0, "", "effect over the range", false);
                    ui::endTable();
                }
                ui::header(study_.proxyMax ? "Proxy optimum (maximum)" : "Proxy optimum (minimum)");
                if (ui::beginTable("##opt", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                    ui::tableColumn("Parameter", 0, 2.2f); ui::tableColumn("Base"); ui::tableColumn("Optimum");
                    ImGui::TableHeadersRow();
                    for (size_t i = 0; i < PP.size() && i < study_.proxyOpt.size(); ++i) {
                        if (!PP[i].active) continue;
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(paramLabel(PP[i], case_).c_str());
                        ImGui::TableSetColumnIndex(1); ImGui::Text("%.4g", PP[i].base);
                        ImGui::TableSetColumnIndex(2); ImGui::TextColored(col4(Pal::Accent), "%.4g", study_.proxyOpt[i]);
                    }
                    ui::endTable();
                }
                ImGui::Text("Predicted %s: %.5g", kResp[study_.proxyResp], study_.proxyPred);
                if (PM.kind == 1 && !study_.proxyOpt.empty()) {
                    const double sd = proxyStd(PM, PP, study_.proxyOpt);
                    ImGui::SameLine(); ImGui::TextDisabled("+- %.4g (GP 1 sd: add training runs near the optimum if large)", sd);
                }
                if (std::isfinite(study_.proxyActual)) {
                    ImGui::SameLine(); ImGui::TextColored(col4(Pal::Ok), "  simulated: %.5g  (error %.2f%%)", study_.proxyActual,
                                                          100.0 * (study_.proxyPred - study_.proxyActual) / std::max(1e-12, std::fabs(study_.proxyActual)));
                }
                ImGui::BeginDisabled(study_.running || study_.proxyOpt.empty());
                if (ui::accentButton(ICON_PLAY " Validate optimum (simulate it)")) {
                    const Case base = case_; const auto opt = study_.proxyOpt; const auto pp = PP; const Response r = Response(study_.proxyResp); const Economics eco = study_.eco;
                    study_.deferred = [this, base, opt, pp, r, eco]() {
                        startStudy("Proxy optimum validation", [base, opt, pp, r, eco](StudyState& S, const StudyCallbacks& cb) {
                            RunResult rr = runSingle(applyParams(base, pp, opt), cb.cancel);
                            if (cb.onRun) cb.onRun("proxy", opt, rr);
                            const double v = rr.ok ? responseValue(rr, base, r, eco, base.observations, autoSigma(base.observations)) : NAN;
                            std::lock_guard<std::mutex> l(S.m); S.proxyActual = v; S.centerTab = 6;
                        });
                    };
                }
                ImGui::EndDisabled();
                ImGui::SameLine();
                if (ImGui::Button(ICON_CHECK " Apply optimum to the case") && !study_.proxyOpt.empty()) applyStudyValues(PP, study_.proxyOpt, "Proxy optimum");
            }
        }
        ImGui::EndTabItem();
    }

    // ---- log
    if (ui::tabItem(ICON_TERMINAL " Log", nullptr, tabFlags(7))) {
        if (ImGui::SmallButton("Clear")) study_.log.clear();
        ImGui::BeginChild("##wslog", ImVec2(0, 0), ImGuiChildFlags_Borders);
        if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
        ImGuiListClipper clip; clip.Begin(int(study_.log.size()));
        while (clip.Step())
            for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                const std::string& s = study_.log[size_t(i)];
                if (s.rfind("ERROR", 0) == 0) ImGui::TextColored(col4(Pal::Err), "%s", s.c_str());
                else if (s.rfind("===", 0) == 0) ImGui::TextColored(col4(Pal::Accent), "%s", s.c_str());
                else ImGui::TextUnformatted(s.c_str());
            }
        if (g_fonts.mono) ImGui::PopFont();
        if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4) ImGui::SetScrollHereY(1.0f);
        ImGui::EndChild();
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ---------------------------------------------------------------------------- right: run archive
void App::studiesRuns() {
    auto& A = study_.archive;
    const auto& AP = study_.archParams;
    ui::header(ICON_HISTORY " Run archive");
    if (A.empty()) { ui::hint("Every simulation of every study is stored here with its parameter values, misfit and response."); return; }
    std::vector<int> pcols;
    for (int i = 0; i < int(AP.size()) && int(pcols.size()) < 6; ++i) if (AP[size_t(i)].active) pcols.push_back(i);
    const int ncol = 4 + int(pcols.size());
    const float tableH = ImGui::GetContentRegionAvail().y * 0.55f;
    static std::vector<int> order;
    static unsigned orderVer = ~0u; static int sortCol = 2; static bool sortAsc = true;
    if (ui::beginTable("##runs", ncol, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_ScrollX |
                                              ImGuiTableFlags_Sortable | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Resizable, ImVec2(0, tableH))) {
        ImGui::TableSetupScrollFreeze(1, 1);
        ui::tableColumn("#", ImGuiTableColumnFlags_DefaultSort);
        ui::tableColumn("Tag");
        ui::tableColumn("Misfit", ImGuiTableColumnFlags_PreferSortAscending);
        ui::tableColumn("Response");
        for (int i : pcols) ui::tableColumn(paramLabel(AP[size_t(i)], case_).c_str());
        ImGui::TableHeadersRow();
        if (ImGuiTableSortSpecs* ss = ImGui::TableGetSortSpecs()) {
            if (ss->SpecsDirty && ss->SpecsCount > 0) { sortCol = ss->Specs[0].ColumnIndex; sortAsc = ss->Specs[0].SortDirection == ImGuiSortDirection_Ascending; orderVer = ~0u; ss->SpecsDirty = false; }
        }
        if (orderVer != study_.archiveVersion || order.size() != A.size()) {
            order.resize(A.size()); std::iota(order.begin(), order.end(), 0);
            auto key = [&](int i) -> double {
                const auto& a = A[size_t(i)];
                if (sortCol == 0) return a.id;
                if (sortCol == 1) return double(std::hash<std::string>()(a.tag) % 100000);
                if (sortCol == 2) return std::isfinite(a.misfit) ? a.misfit : 1e300;
                if (sortCol == 3) return std::isfinite(a.response) ? a.response : -1e300;
                const int p = pcols[size_t(sortCol - 4)];
                return size_t(p) < a.values.size() ? a.values[size_t(p)] : 0.0;
            };
            std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return sortAsc ? key(x) < key(y) : key(x) > key(y); });
            orderVer = study_.archiveVersion;
        }
        ImGuiListClipper clip; clip.Begin(int(order.size()));
        while (clip.Step())
            for (int r = clip.DisplayStart; r < clip.DisplayEnd; ++r) {
                const auto& a = A[size_t(order[size_t(r)])];
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                char b[32]; std::snprintf(b, sizeof(b), "%d", a.id);
                if (ImGui::Selectable(b, study_.selRun == a.id, ImGuiSelectableFlags_SpanAllColumns)) study_.selRun = a.id;
                ImGui::TableSetColumnIndex(1); ImGui::TextColored(col4(tagColor(a.tag)), "%s", a.tag.c_str());
                ImGui::TableSetColumnIndex(2);
                if (!a.ok) ImGui::TextColored(col4(Pal::Err), "failed");
                else if (std::isfinite(a.misfit)) ImGui::Text("%.4g", a.misfit); else ImGui::TextDisabled("-");
                ImGui::TableSetColumnIndex(3); if (std::isfinite(a.response)) ImGui::Text("%.5g", a.response); else ImGui::TextDisabled("-");
                for (size_t c = 0; c < pcols.size(); ++c) {
                    ImGui::TableSetColumnIndex(4 + int(c));
                    if (size_t(pcols[c]) < a.values.size()) ImGui::Text("%.4g", a.values[size_t(pcols[c])]);
                }
            }
        ui::endTable();
    }
    if (ImGui::SmallButton(ICON_DOWNLOAD " Export CSV")) {
        std::string path = "run_archive.csv";
        if (fileDialog(true, "CSV (*.csv)\0*.csv\0", "csv", path)) {
            std::ofstream f(path);
            f << "id,tag,ok,misfit,response";
            for (auto& p : AP) f << ",\"" << paramLabel(p, case_) << "\"";
            for (auto& t : study_.termNames) f << ",\"" << t << "\"";
            f << "\n";
            for (auto& a : A) {
                f << a.id << "," << a.tag << "," << (a.ok ? 1 : 0) << "," << a.misfit << "," << a.response;
                for (size_t i = 0; i < AP.size(); ++i) f << "," << (i < a.values.size() ? a.values[i] : NAN);
                for (size_t i = 0; i < study_.termNames.size(); ++i) f << "," << (i < a.terms.size() ? a.terms[i] : NAN);
                f << "\n";
            }
            study_.deferred = [this, path]() { notify("Run archive exported: " + path); };
        }
    }
    ImGui::SameLine();
    if (ImGui::SmallButton(ICON_TRASH_2 " Clear archive")) { A.clear(); study_.selRun = -1; ++study_.archiveVersion; }

    // ---- details of the selected run
    const StudyState::ArchRun* s = nullptr;
    for (auto& a : A) if (a.id == study_.selRun) s = &a;
    ui::header(ICON_INFO " Selected run");
    if (!s) { ui::hint("Select a run in the table, the cross-plot or the Pareto chart."); return; }
    ImGui::BeginChild("##rundet");
    ImGui::Text("Run #%d  (%s)", s->id, s->tag.c_str());
    if (!s->ok) ImGui::TextColored(col4(Pal::Err), "Run failed: %s", s->run.error.c_str());
    else {
        ImGui::Text("Misfit %.5g    Response %.6g", s->misfit, s->response);
        ImGui::TextDisabled("%d time steps, %d Newton, %.2f s", s->run.stats.timesteps, s->run.stats.newton, s->run.stats.seconds);
    }
    if (ui::beginTable("##rv", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingStretchProp)) {
        for (size_t i = 0; i < AP.size() && i < s->values.size(); ++i) {
            if (!AP[i].active) continue;
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(paramLabel(AP[i], case_).c_str());
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.5g", s->values[i]);
        }
        for (size_t i = 0; i < study_.termNames.size() && i < s->terms.size(); ++i) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("term %s", study_.termNames[i].c_str());
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.4g", s->terms[i]);
        }
        ui::endTable();
    }
    if (s->ok && ui::accentButton(ICON_CHECK " Apply to case")) applyStudyValues(AP, s->values, "Run #" + std::to_string(s->id));
    ImGui::SameLine();
    if (ImGui::Button(ICON_SLIDERS_HORIZONTAL " To sliders")) {
        for (size_t i = 0; i < study_.manualVals.size() && i < s->values.size(); ++i) study_.manualVals[i] = s->values[i];
        study_.wsMethod = 5;
    }
    ImGui::EndChild();
}

// ============================================================================ history matching workspace
namespace {
ImU32 misfitColour(double v) {   // normalised RMS: < 1 within the measurement error (green) ... > 3 poor (red)
    if (!std::isfinite(v)) return IM_COL32(90, 90, 98, 160);
    const double t = std::clamp((v - 0.5) / 2.5, 0.0, 1.0);
    const ImVec4 g(0.22f, 0.62f, 0.33f, 0.85f), y(0.80f, 0.66f, 0.18f, 0.85f), r(0.78f, 0.24f, 0.22f, 0.85f);
    const float u = float(t < 0.5 ? t * 2 : (t - 0.5) * 2);
    const ImVec4 a = t < 0.5 ? g : y, b = t < 0.5 ? y : r;
    return ImGui::ColorConvertFloat4ToU32(ImVec4(a.x + (b.x - a.x) * u, a.y + (b.y - a.y) * u, a.z + (b.z - a.z) * u, a.w));
}
}  // namespace

// observed data from an ECLIPSE / OPM summary file (history vectors WOPRH ... or simulated vectors)
void App::obsImportSummaryPopup() {
    if (!ImGui::BeginPopup("obssum")) return;
    ImGui::TextDisabled("Observed data from an ECLIPSE / OPM summary file (.SMSPEC + .UNSMRY)");
    ImGui::RadioButton("History vectors (WOPRH, WWPRH, WBHPH, ...)", &study_.sumHistory, 1);
    ImGui::RadioButton("Simulated vectors (another run as the truth)", &study_.sumHistory, 0);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
    ImGui::DragFloat("Minimum spacing (days)", &study_.sumSpacing, 1.0f, 0.0f, 3650.0f, "%.0f");
    ImGui::TextDisabled("Quantities");
    for (int q = 0; q < kNumObs; ++q) {
        bool on = (study_.sumMask >> q) & 1u;
        if (q % 7) ImGui::SameLine();
        if (ImGui::Checkbox(kObsNames[q], &on)) study_.sumMask = on ? (study_.sumMask | (1u << q)) : (study_.sumMask & ~(1u << q));
    }
    ImGui::Checkbox("Append to the existing observations", &study_.sumAppend);
    if (ui::accentButton(ICON_FOLDER_OPEN " Choose summary file...")) {
        std::string p;
        if (fileDialog(false, "Summary (*.SMSPEC;*.FSMSPEC;*.DATA)\0*.SMSPEC;*.FSMSPEC;*.DATA\0All files (*.*)\0*.*\0", "SMSPEC", p)) {
            ecl::Summary s; std::string err;
            const ecl::CaseFiles f = ecl::locateCase(p);
            auto fileNameOf = [](const std::string& x) { const size_t k = x.find_last_of("/\\"); return k == std::string::npos ? x : x.substr(k + 1); };
            if (f.smspec.empty() || !ecl::readSummary(f.smspec, s, &err)) notify("Cannot read the summary: " + (err.empty() ? p : err));
            else {
                SummaryObsOptions o; o.history = study_.sumHistory != 0; o.minSpacingDays = study_.sumSpacing; o.qtyMask = study_.sumMask;
                std::string note;
                auto obs = observationsFromSummary(s, case_, o, &note);
                if (obs.empty()) notify("No usable vectors in " + fileNameOf(f.smspec) + (o.history ? " (no history vectors - try 'Simulated vectors')" : ""));
                else {
                    if (study_.sumAppend) case_.observations.insert(case_.observations.end(), obs.begin(), obs.end());
                    else case_.observations = obs;
                    markStudyDirty();
                    notify(note);
                    study_.log.push_back("Observed data: " + fileNameOf(f.smspec) + " - " + note);
                }
            }
            ImGui::CloseCurrentPopup();
        }
    }
    ImGui::EndPopup();
}

// runs compared in the match-quality views: initial / best of the HM study, manual, selected archive run, last run
void App::hmColumns(std::vector<std::string>& names, std::vector<const RunResult*>& runs, RunResult& lastHolder) {
    names.clear(); runs.clear();
    if (study_.hasHM && study_.hm.initialRun.ok) { names.push_back("Initial"); runs.push_back(&study_.hm.initialRun); }
    if (study_.hasHM && study_.hm.bestRun.ok) { names.push_back(std::string("Best (") + hmName(true) + ")"); runs.push_back(&study_.hm.bestRun); }
    if (study_.hasField && study_.field.ok && study_.field.bestRun.ok) {   // C13: per-cell adjoint field match
        if (!(study_.hasHM && study_.hm.initialRun.ok) && study_.field.initialRun.ok) { names.push_back("Initial"); runs.push_back(&study_.field.initialRun); }
        names.push_back("Best (field adjoint)"); runs.push_back(&study_.field.bestRun);
    }
    if (study_.hasManual && study_.manual.ok) { names.push_back("Manual"); runs.push_back(&study_.manual); }
    for (auto& a : study_.archive)
        if (a.id == study_.selRun && a.ok) { names.push_back("Run #" + std::to_string(a.id)); runs.push_back(&a.run); }
    if (runs.empty()) {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        if (R.valid && R.summary.size() > 1 && !runner_.running()) {
            lastHolder.ok = true; lastHolder.summary = R.summary;
            names.push_back("Last run"); runs.push_back(&lastHolder);
        }
    }
}

void App::hmMatchQuality() {
    const float fs = ImGui::GetFontSize();
    const auto& O = case_.observations;
    if (O.empty()) { ui::hint("No observed data. Load it in the 'Observed data' tab (summary file history vectors, CSV, WCONHIST or a synthetic truth)."); return; }
    static RunResult lastHolder;
    std::vector<std::string> names; std::vector<const RunResult*> runs;
    hmColumns(names, runs, lastHolder);
    if (runs.empty()) { ui::hint("No simulation to compare yet: run the model, a manual trial or a history-matching study."); return; }
    const auto groups = obsGroups(case_);
    // cache: residuals are evaluated for every observation of every compared run
    double chk = double(O.size()); for (auto& o : O) chk += o.value * 1e-3 + o.time * 1e-7 + o.weight + o.sigma * 1e-5;
    const unsigned key = study_.archiveVersion * 2654435761u ^ study_.manualVersion * 40503u ^ unsigned(study_.selRun * 977) ^
                         unsigned(study_.hasHM ? 7 : 0) ^ unsigned(runs.size() * 131);
    auto& C = mq_;
    if (C.key != key || C.chk != chk || C.cols != names || C.runs != runs) {
        C.key = key; C.chk = chk; C.cols = names; C.runs = runs;
        const auto sig = autoSigma(O);
        C.rms.assign(runs.size(), std::vector<double>(groups.size(), NAN));
        C.total.assign(runs.size(), NAN);
        for (size_t r = 0; r < runs.size(); ++r) {
            const auto res = residuals(*runs[r], case_, O, sig);
            double tot = 0; int nt = 0;
            for (size_t g = 0; g < groups.size(); ++g) {
                double s = 0; int n = 0;
                for (int i : groups[g].idx) if (O[size_t(i)].weight > 0 && std::isfinite(res[size_t(i)])) { s += res[size_t(i)] * res[size_t(i)]; ++n; }
                if (n) { C.rms[r][g] = std::sqrt(s / n); tot += s; nt += n; }
            }
            if (nt) C.total[r] = tot / nt;
        }
    }
    // ---- headline
    for (size_t r = 0; r < runs.size(); ++r) {
        int good = 0, all = 0;
        for (double v : C.rms[r]) if (std::isfinite(v)) { ++all; if (v <= 1.0) ++good; }
        if (r) ImGui::SameLine(0, fs * 1.5f);
        ImGui::BeginGroup();
        ImGui::TextDisabled("%s", names[r].c_str());
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(misfitColour(std::sqrt(std::max(0.0, C.total[r])))), "%.4g", C.total[r]);
        if (g_fonts.title) ImGui::PopFont();
        ImGui::TextDisabled("%d / %d series within error", good, all);
        ImGui::EndGroup();
    }
    if (runs.size() >= 2 && std::isfinite(C.total[0]) && C.total[0] > 0) {
        ImGui::SameLine(0, fs * 2);
        ImGui::BeginGroup();
        ImGui::TextDisabled("Misfit reduction");
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::TextColored(col4(Pal::Ok), "%.1f %%", 100.0 * (1.0 - C.total[1] / C.total[0]));
        if (g_fonts.title) ImGui::PopFont();
        ImGui::TextDisabled("%s vs %s", names[1].c_str(), names[0].c_str());
        ImGui::EndGroup();
    }
    ImGui::Separator();
    ImGui::TextDisabled("View"); ImGui::SameLine();
    ImGui::RadioButton("Series x runs", &study_.mqView, 0); ImGui::SameLine();
    ImGui::RadioButton("Wells x quantities", &study_.mqView, 1);
    if (study_.mqView == 1) {
        ImGui::SameLine(0, fs);
        study_.mqRun = std::clamp(study_.mqRun, 0, int(runs.size()) - 1);
        ImGui::SetNextItemWidth(fs * 9);
        if (ImGui::BeginCombo("##mqrun", names[size_t(study_.mqRun)].c_str())) {
            for (int r = 0; r < int(runs.size()); ++r) if (ImGui::Selectable(names[size_t(r)].c_str(), r == study_.mqRun)) study_.mqRun = r;
            ImGui::EndCombo();
        }
    }
    ImGui::SameLine(0, fs);
    ImGui::TextDisabled("cell = RMS of (sim - obs) / sigma:  <1 within error,  >3 poor");
    const float tableH = std::max(fs * 8, (ImGui::GetContentRegionAvail().y) * 0.48f);
    if (study_.mqView == 0) {
        const int nc = 2 + int(runs.size()) + (runs.size() >= 2 ? 1 : 0);
        if (ui::beginTable("mq", nc, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Sortable |
                                            ImGuiTableFlags_SizingStretchProp, ImVec2(0, tableH))) {
            ImGui::TableSetupScrollFreeze(1, 1);
            ui::tableColumn("Series", ImGuiTableColumnFlags_WidthStretch, 2.0f);
            ui::tableColumn("Points", ImGuiTableColumnFlags_WidthStretch, 0.7f);
            for (auto& n : names) ui::tableColumn(n.c_str(), ImGuiTableColumnFlags_WidthStretch, 1.0f);
            if (runs.size() >= 2) ui::tableColumn("Change", ImGuiTableColumnFlags_WidthStretch, 0.9f);
            ImGui::TableHeadersRow();
            // order: sort by the selected column (default: worst series of the last compared run first)
            std::vector<int> order(groups.size()); std::iota(order.begin(), order.end(), 0);
            int sc = int(runs.size()) - 1 + 2; bool asc = false;
            if (ImGuiTableSortSpecs* ss = ImGui::TableGetSortSpecs()) if (ss->SpecsCount > 0) { sc = ss->Specs[0].ColumnIndex; asc = ss->Specs[0].SortDirection == ImGuiSortDirection_Ascending; }
            auto val = [&](int g) -> double {
                if (sc == 0) return 0;
                if (sc == 1) return double(groups[size_t(g)].idx.size());
                if (sc - 2 < int(runs.size())) return C.rms[size_t(sc - 2)][size_t(g)];
                const double a = C.rms[0][size_t(g)], b = C.rms[1][size_t(g)]; return a > 0 ? b / a : NAN;
            };
            if (sc == 0) std::sort(order.begin(), order.end(), [&](int a, int b) { return asc ? groupName(groups[size_t(a)]) < groupName(groups[size_t(b)]) : groupName(groups[size_t(a)]) > groupName(groups[size_t(b)]); });
            else std::stable_sort(order.begin(), order.end(), [&](int a, int b) { double x = val(a), y = val(b); if (!std::isfinite(x)) x = -1; if (!std::isfinite(y)) y = -1; return asc ? x < y : x > y; });
            ImGuiListClipper clip; clip.Begin(int(order.size()));
            while (clip.Step())
                for (int row = clip.DisplayStart; row < clip.DisplayEnd; ++row) {
                    const int g = order[size_t(row)];
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    if (ImGui::Selectable(groupName(groups[size_t(g)]).c_str(), study_.plotObs == g, ImGuiSelectableFlags_SpanAllColumns)) study_.plotObs = g;
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%d", int(groups[size_t(g)].idx.size()));
                    for (size_t r = 0; r < runs.size(); ++r) {
                        ImGui::TableSetColumnIndex(2 + int(r));
                        const double v = C.rms[r][size_t(g)];
                        ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, misfitColour(v));
                        if (std::isfinite(v)) ImGui::Text("%.3g", v); else ImGui::TextDisabled("-");
                    }
                    if (runs.size() >= 2) {
                        ImGui::TableSetColumnIndex(2 + int(runs.size()));
                        const double a = C.rms[0][size_t(g)], b = C.rms[1][size_t(g)];
                        if (std::isfinite(a) && std::isfinite(b) && a > 0) {
                            const double ch = 100.0 * (b / a - 1.0);
                            ImGui::TextColored(col4(ch <= 0 ? Pal::Ok : Pal::Err), "%+.0f %%", ch);
                        } else ImGui::TextDisabled("-");
                    }
                }
            ui::endTable();
        }
    } else {
        // wells x quantities heat map of one run
        std::vector<std::string> wells; std::vector<int> qtys;
        for (auto& g : groups) {
            const std::string w = g.well.empty() ? "FIELD" : g.well;
            if (std::find(wells.begin(), wells.end(), w) == wells.end()) wells.push_back(w);
            if (std::find(qtys.begin(), qtys.end(), int(g.qty)) == qtys.end()) qtys.push_back(int(g.qty));
        }
        std::sort(qtys.begin(), qtys.end());
        if (ui::beginTable("mqh", 1 + int(qtys.size()) + 1, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_ScrollX,
                              ImVec2(0, tableH))) {
            ImGui::TableSetupScrollFreeze(1, 1);
            ui::tableColumn("Well", ImGuiTableColumnFlags_WidthFixed, fs * 7);
            for (int q : qtys) ui::tableColumn(kObsNames[q], ImGuiTableColumnFlags_WidthFixed, fs * 4.2f);
            ui::tableColumn("Well RMS", ImGuiTableColumnFlags_WidthFixed, fs * 5);
            ImGui::TableHeadersRow();
            const auto& R = C.rms[size_t(study_.mqRun)];
            for (auto& w : wells) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(w.c_str());
                double s = 0; int n = 0;
                for (size_t qi = 0; qi < qtys.size(); ++qi) {
                    ImGui::TableSetColumnIndex(1 + int(qi));
                    int g = -1;
                    for (int k = 0; k < int(groups.size()); ++k)
                        if (int(groups[size_t(k)].qty) == qtys[qi] && (groups[size_t(k)].well.empty() ? "FIELD" : groups[size_t(k)].well) == w) g = k;
                    if (g < 0) { ImGui::TextDisabled(" "); continue; }
                    const double v = R[size_t(g)];
                    ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, misfitColour(v));
                    char b[32];
                    if (std::isfinite(v)) std::snprintf(b, sizeof b, "%.2g##%d", v, g); else std::snprintf(b, sizeof b, "-##%d", g);
                    if (ImGui::Selectable(b, study_.plotObs == g)) study_.plotObs = g;
                    if (std::isfinite(v)) { s += v * v; ++n; }
                }
                ImGui::TableSetColumnIndex(1 + int(qtys.size()));
                if (n) { const double v = std::sqrt(s / n); ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, misfitColour(v)); ImGui::Text("%.2g", v); }
            }
            ui::endTable();
        }
    }
    // ---- profile of the selected series for every compared run
    std::vector<const char*> labels; for (auto& n : names) labels.push_back(n.c_str());
    char id[160]; std::snprintf(id, sizeof id, "%s##mqplot", groups.empty() ? "" : groupName(groups[size_t(std::clamp(study_.plotObs, 0, int(groups.size()) - 1))]).c_str());
    plotObservationFit(id, case_, runs, labels, study_.plotObs, std::max(fs * 10, ImGui::GetContentRegionAvail().y - 4));
}

void App::hmParameters() {
    const float fs = ImGui::GetFontSize();
    if (!study_.hasHM) {
        ui::hint("Run Levenberg-Marquardt or ES-MDA (Method tab): this view then shows the matched parameter values against their "
                 "prior ranges, the LM search path and the ES-MDA posterior distributions (P10 / P50 / P90).");
        return;
    }
    const auto& H = study_.hm; const auto& P = study_.hmParams;
    const bool ens = study_.hmMethod == 1 && !H.ensemble.empty();
    std::vector<int> act; for (int i = 0; i < int(P.size()); ++i) if (P[size_t(i)].active) act.push_back(i);
    ImGui::Text("%s: misfit %.4g -> %.4g  (%d simulations)", hmName(false), H.initialMisfit, H.bestMisfit, H.runs);
    ImGui::SameLine(0, fs * 2);
    if (ui::accentButton(ICON_CHECK " Apply best values to the model")) {
        auto hp = P; auto v = H.best;
        study_.deferred = [this, hp, v] { applyStudyValues(hp, v, "history-match result"); };
    }
    const int nc = ens ? 8 : 5;
    if (ui::beginTable("hmp", nc, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ui::tableColumn("Parameter", 0, 2.2f); ui::tableColumn("Range", 0, 1.4f); ui::tableColumn("Base", 0, 0.9f);
        ui::tableColumn(ens ? "Best member" : "Matched", 0, 0.9f); ui::tableColumn("Position in range", 0, 2.6f);
        if (ens) { ui::tableColumn("P10", 0, 0.8f); ui::tableColumn("P50", 0, 0.8f); ui::tableColumn("P90", 0, 0.8f); }
        ImGui::TableHeadersRow();
        for (int i : act) {
            const StudyParam& p = P[size_t(i)];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(paramLabel(p, case_).c_str());
            ImGui::TableSetColumnIndex(1); ImGui::Text("%s - %s%s", fmtv(p.low).c_str(), fmtv(p.high).c_str(), p.log ? " (log)" : "");
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(fmtv(p.base).c_str());
            const double best = size_t(i) < H.best.size() ? H.best[size_t(i)] : p.base;
            ImGui::TableSetColumnIndex(3); ImGui::TextColored(col4(Pal::Accent), "%s", fmtv(best).c_str());
            // range bar: prior range, P10-P90 band (ensemble), base (grey) and matched (accent) markers
            ImGui::TableSetColumnIndex(4);
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImVec2 a = ImGui::GetCursorScreenPos(); const float w = ImGui::GetContentRegionAvail().x, h = fs * 1.1f;
            const float y = a.y + h * 0.5f;
            dl->AddRectFilled(ImVec2(a.x, y - 2), ImVec2(a.x + w, y + 2), IM_COL32(120, 125, 140, 110), 2);
            auto X = [&](double v) { return a.x + w * float(std::clamp(paramToUnit(p, v), 0.0, 1.0)); };
            double p10 = NAN, p50 = NAN, p90 = NAN;
            if (ens) {
                std::vector<double> v; for (auto& m : H.ensemble) if (size_t(i) < m.size()) v.push_back(m[size_t(i)]);
                std::sort(v.begin(), v.end());
                auto q = [&](double f) { return v.empty() ? NAN : v[size_t(std::clamp(f * (v.size() - 1), 0.0, double(v.size() - 1)))]; };
                p10 = q(0.1); p50 = q(0.5); p90 = q(0.9);
                if (std::isfinite(p10)) dl->AddRectFilled(ImVec2(X(p10), y - h * 0.32f), ImVec2(X(p90), y + h * 0.32f), IM_COL32(64, 132, 246, 90), 3);
            }
            dl->AddCircleFilled(ImVec2(X(p.base), y), h * 0.26f, IM_COL32(170, 175, 190, 255));
            dl->AddTriangleFilled(ImVec2(X(best), y - h * 0.45f), ImVec2(X(best) - h * 0.3f, y + h * 0.4f), ImVec2(X(best) + h * 0.3f, y + h * 0.4f), Pal::Accent);
            ImGui::Dummy(ImVec2(w, h));
            if (ens) {
                ImGui::TableSetColumnIndex(5); ImGui::TextUnformatted(fmtv(p10).c_str());
                ImGui::TableSetColumnIndex(6); ImGui::TextUnformatted(fmtv(p50).c_str());
                ImGui::TableSetColumnIndex(7); ImGui::TextUnformatted(fmtv(p90).c_str());
            }
        }
        ui::endTable();
    }
    ImGui::TextDisabled("grey dot = base value, blue triangle = matched value%s", ens ? ", band = posterior P10 - P90" : "");
    hmStatsPanel();
    const float h = std::max(fs * 10, ImGui::GetContentRegionAvail().y - 4);
    if (ens) {
        // posterior histograms, prior range as the x axis
        const int cols = std::max(1, std::min(3, int(act.size())));
        const int rows = (int(act.size()) + cols - 1) / cols;
        if (rows > 0 && ImPlot::BeginSubplots("Posterior ensemble (ES-MDA)", rows, cols, ImVec2(-1, h))) {
            for (int i : act) {
                const StudyParam& p = P[size_t(i)];
                std::vector<double> v; for (auto& m : H.ensemble) if (size_t(i) < m.size() && std::isfinite(m[size_t(i)])) v.push_back(m[size_t(i)]);
                if (ImPlot::BeginPlot(paramLabel(p, case_).c_str())) {
                    ImPlot::SetupAxes(nullptr, "members", 0, ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxisLimits(ImAxis_X1, p.low, p.high, ImPlotCond_Always);
                    if (p.log && p.low > 0) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                    ImPlotSpec s; s.FillColor = col4(Pal::Accent); s.FillAlpha = 0.6f;
                    if (!v.empty()) ImPlot::PlotHistogram("posterior", v.data(), int(v.size()), 12, 1.0, ImPlotRange(p.low, p.high), s);
                    const double base[1] = {p.base}, best[1] = {size_t(i) < H.best.size() ? H.best[size_t(i)] : p.base};
                    ImPlotSpec sb; sb.LineColor = col4(Pal::Muted); sb.LineWeight = 2;
                    ImPlot::PlotInfLines("base", base, 1, sb);
                    ImPlotSpec sm; sm.LineColor = col4(Pal::Ok); sm.LineWeight = 2;
                    ImPlot::PlotInfLines("best", best, 1, sm);
                    ImPlot::EndPlot();
                }
            }
            ImPlot::EndSubplots();
        }
    } else if (!H.history.empty() && ImPlot::BeginPlot(study_.hmAlgo == 2 ? "Adjoint L-BFGS search path (normalised parameter position)" : "Levenberg-Marquardt search path (normalised parameter position)", ImVec2(-1, h))) {
        ImPlot::SetupAxes("iteration", "position in range (0 = low, 1 = high)", ImPlotAxisFlags_AutoFit, 0);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -0.02, 1.02, ImPlotCond_Always);
        int ci = 0;
        for (int i : act) {
            std::vector<double> x, y;
            for (auto& it : H.history) if (size_t(i) < it.values.size()) { x.push_back(it.iter); y.push_back(paramToUnit(P[size_t(i)], it.values[size_t(i)])); }
            ImPlotSpec s; s.LineColor = col4(seriesCol(ci++)); s.LineWeight = 2; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 3;
            if (!x.empty()) ImPlot::PlotLine(paramLabel(P[size_t(i)], case_).c_str(), x.data(), y.data(), int(x.size()), s);
        }
        ImPlot::EndPlot();
    }
}

// ============================================================================ C13: Adjoint HM view (centre tab 11)
// Dedicated result view of the adjoint-gradient history matching: scalar parameters (projected L-BFGS) or the per-cell
// field (ln k / ln phi multiplier of every active cell). Sub-views: iterations (misfit, gradient norm, runs, time per
// accepted iterate), final parameters / field maps per layer at any iteration (with a 3D-viewer export), per-response
// check of initial / selected iterate / best (normalised RMS per observation series + overlay plots) and the gradient
// (sensitivity) of the misfit.
void App::adjointHMView() {
    const float fs = ImGui::GetFontSize();
    const bool hasS = study_.hasHM && study_.hmMethod == 0 && study_.hmAlgo == 2 && study_.hm.ok;
    const bool hasF = study_.hasField && study_.field.ok;
    if (!hasS && !hasF) {
        ui::hint("No adjoint history match yet. Method tab > 'Adjoint gradient (L-BFGS / per-cell field)':\n"
                 " - Scalar parameters: the active study parameters, projected L-BFGS on the exact adjoint gradient.\n"
                 " - Per-cell field: ln k / ln phi multipliers of every active cell with a smoothness term.\n"
                 "This view then shows the iterations, the final parameters / field maps, every response of the best case\n"
                 "against the observed data and the misfit sensitivities.");
        if (study_.hasField && !study_.field.ok) ImGui::TextColored(col4(Pal::Err), "Last field match failed: %s", study_.field.error.c_str());
        return;
    }
    // which result
    if (hasS && hasF) {
        ImGui::TextDisabled("Result:"); ImGui::SameLine();
        if (ImGui::RadioButton("Scalar parameters (L-BFGS)", study_.adjLast == 0)) { study_.adjLast = 0; study_.adjIter = -1; }
        ImGui::SameLine();
        if (ImGui::RadioButton("Per-cell field", study_.adjLast == 1)) { study_.adjLast = 1; study_.adjIter = -1; }
    } else study_.adjLast = hasF ? 1 : 0;
    const bool field = study_.adjLast == 1;
    const HMResult& H = study_.hm;
    const FieldHMResult& F = study_.field;
    const int nIt = field ? int(F.iters.size()) : int(H.history.size());
    if (study_.adjIter < 0 || study_.adjIter >= nIt) study_.adjIter = std::max(0, nIt - 1);
    const int K = study_.adjIter;
    // ---- headline
    {
        const double m0 = field ? F.initialMisfit : H.initialMisfit, mb = field ? F.bestMisfit : H.bestMisfit;
        const int runs = field ? F.runs : H.runs;
        ImGui::BeginGroup(); ImGui::TextDisabled("Initial misfit");
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::Text("%.4g", m0);
        if (g_fonts.title) ImGui::PopFont();
        ImGui::EndGroup(); ImGui::SameLine(0, fs * 2);
        ImGui::BeginGroup(); ImGui::TextDisabled("Best misfit");
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::TextColored(col4(Pal::Accent), "%.4g", mb);
        if (g_fonts.title) ImGui::PopFont();
        ImGui::EndGroup(); ImGui::SameLine(0, fs * 2);
        ImGui::BeginGroup(); ImGui::TextDisabled("Reduction");
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::Text("%.1f %%", m0 > 0 ? 100.0 * (1.0 - mb / m0) : 0.0);
        if (g_fonts.title) ImGui::PopFont();
        ImGui::EndGroup(); ImGui::SameLine(0, fs * 2);
        ImGui::BeginGroup(); ImGui::TextDisabled("Cost");
        ImGui::Text("%d iterations, %d simulations", std::max(0, nIt - 1), runs);
        if (field) ImGui::TextDisabled("forward %.1f s, adjoint %.1f s", F.secForward, F.secAdjoint);
        else if (!H.history.empty()) ImGui::TextDisabled("wall %.1f s", H.history.back().seconds);
        ImGui::EndGroup();
        if (field) { ImGui::SameLine(0, fs * 2); ImGui::BeginGroup(); ImGui::TextDisabled("Unknowns"); ImGui::Text("%d x %d x %d cells", F.nx, F.ny, F.nz); ImGui::EndGroup(); }
    }
    // iteration selector (shared by every sub-view)
    if (nIt > 1) {
        ImGui::SetNextItemWidth(fs * 16);
        ImGui::SliderInt("Iteration", &study_.adjIter, 0, nIt - 1, study_.adjIter == 0 ? "0 (prior)" : "%d");
        ImGui::SameLine(); ImGui::TextDisabled("selected iterate is compared with the prior and the best match");
    }
    const RunResult* r0 = field ? &F.initialRun : &H.initialRun;
    const RunResult* rb = field ? &F.bestRun : &H.bestRun;
    const RunResult* rk = nullptr;
    if (field && K < int(F.iters.size()) && F.iters[size_t(K)].run.ok) rk = &F.iters[size_t(K)].run;
    if (!field && K < int(H.history.size()) && H.history[size_t(K)].run.ok) rk = &H.history[size_t(K)].run;

    if (!ImGui::BeginTabBar("##adjsub")) return;
    // ------------------------------------------------------------------ iterations
    if (ui::tabItem(ICON_TRENDING_UP " Iterations")) {
        const float h = std::max(fs * 12, ImGui::GetContentRegionAvail().y * 0.45f);
        if (ImPlot::BeginSubplots("##adjit", 1, 2, ImVec2(-1, h))) {
            std::vector<double> it, mis, gn, obj, reg;
            for (int k = 0; k < nIt; ++k) {
                it.push_back(k);
                if (field) { const auto& I = F.iters[size_t(k)]; mis.push_back(I.misfit); gn.push_back(I.gradNorm); obj.push_back(I.objective); reg.push_back(I.reg); }
                else { const auto& I = H.history[size_t(k)]; mis.push_back(I.misfit); gn.push_back(I.gradNorm); }
            }
            if (ImPlot::BeginPlot("Misfit and gradient norm")) {
                ImPlot::SetupAxes("iteration", "misfit", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlot::SetupAxis(ImAxis_Y2, "|gradient|", ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_AuxDefault);
                ImPlot::SetupAxisScale(ImAxis_Y2, ImPlotScale_Log10);
                auto pos = [](std::vector<double> v) { for (double& x : v) x = std::max(x, 1e-30); return v; };
                const auto mp = pos(mis), gp = pos(gn);
                ImPlotSpec sp; sp.LineColor = col4(Pal::Accent); sp.LineWeight = 2.5f; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3.5f;
                ImPlot::PlotLine("data misfit", it.data(), mp.data(), nIt, sp);
                if (field) {
                    const auto op = pos(obj);
                    ImPlotSpec so; so.LineColor = ImVec4(1, 1, 1, 0.6f); so.LineWeight = 1.2f;
                    ImPlot::PlotLine("objective (misfit + smoothness)", it.data(), op.data(), nIt, so);
                }
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                ImPlotSpec sg; sg.LineColor = ImVec4(0.55f, 0.75f, 1.0f, 1); sg.LineWeight = 1.6f; sg.Marker = ImPlotMarker_Square; sg.MarkerSize = 2.5f;
                ImPlot::PlotLine("gradient norm", it.data(), gp.data(), nIt, sg);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
                const double kx = K;
                ImPlotSpec sc; sc.LineColor = ImVec4(1, 1, 1, 0.35f); sc.Flags = ImPlotItemFlags_NoLegend;
                ImPlot::PlotInfLines("##sel", &kx, 1, sc);
                if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0)) study_.adjIter = std::clamp(int(std::lround(ImPlot::GetPlotMousePos().x)), 0, nIt - 1);
                ImPlot::EndPlot();
            }
            if (field) {
                if (ImPlot::BeginPlot("Multiplier range per iteration")) {
                    ImPlot::SetupAxes("iteration", "multiplier", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                    std::vector<double> kmin, kmax, kmean, pmin, pmax, pmean;
                    for (const auto& I : F.iters) { kmin.push_back(I.kMin); kmax.push_back(I.kMax); kmean.push_back(I.kMean); pmin.push_back(I.phiMin); pmax.push_back(I.phiMax); pmean.push_back(I.phiMean); }
                    ImPlotSpec sk; sk.FillColor = ImVec4(1.0f, 0.6f, 0.2f, 1); sk.FillAlpha = 0.25f;
                    ImPlot::PlotShaded("k range", it.data(), kmin.data(), kmax.data(), nIt, sk);
                    ImPlotSpec lk; lk.LineColor = ImVec4(1.0f, 0.6f, 0.2f, 1); lk.LineWeight = 2.0f;
                    ImPlot::PlotLine("k mean", it.data(), kmean.data(), nIt, lk);
                    if (!F.multPhi.empty()) {
                        ImPlotSpec sp2; sp2.FillColor = ImVec4(0.4f, 0.8f, 0.5f, 1); sp2.FillAlpha = 0.25f;
                        ImPlot::PlotShaded("phi range", it.data(), pmin.data(), pmax.data(), nIt, sp2);
                        ImPlotSpec lp; lp.LineColor = ImVec4(0.4f, 0.8f, 0.5f, 1); lp.LineWeight = 2.0f;
                        ImPlot::PlotLine("phi mean", it.data(), pmean.data(), nIt, lp);
                    }
                    ImPlot::EndPlot();
                }
            } else if (ImPlot::BeginPlot("Normalised parameters per iteration")) {
                ImPlot::SetupAxes("iteration", "u = position in [low, high]", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxisLimits(ImAxis_Y1, -0.02, 1.02, ImPlotCond_Always);
                const auto& P = study_.hmParams;
                int ci = 0;
                for (size_t i = 0; i < P.size(); ++i) {
                    if (!P[i].active) continue;
                    std::vector<double> u;
                    for (const auto& I : H.history) u.push_back(i < I.values.size() ? paramToUnit(P[i], I.values[i]) : NAN);
                    ImPlotSpec s; s.LineColor = col4(seriesCol(ci++)); s.LineWeight = 2.0f; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 2.5f;
                    ImPlot::PlotLine(paramLabel(P[i], case_).c_str(), it.data(), u.data(), nIt, s);
                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        // table
        const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp;
        const int nc = field ? 11 : 6;
        if (ui::beginTable("##adjtab", nc, tf, ImVec2(0, std::max(fs * 8, ImGui::GetContentRegionAvail().y)))) {
            ImGui::TableSetupScrollFreeze(0, 1);
            if (field) for (const char* c : {"Iter", "Misfit", "Objective", "Smoothness", "|grad|", "ln step", "k min", "k max", "k mean", "Runs", "Time (s)"}) ui::tableColumn(c);
            else for (const char* c : {"Iter", "Misfit", "|grad|", "Runs", "Time (s)", "Reduction vs prior"}) ui::tableColumn(c);
            ImGui::TableHeadersRow();
            for (int k = 0; k < nIt; ++k) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                char b[24]; std::snprintf(b, sizeof b, "%d", k);
                if (ImGui::Selectable(b, k == K, ImGuiSelectableFlags_SpanAllColumns)) study_.adjIter = k;
                if (field) {
                    const auto& I = F.iters[size_t(k)];
                    ImGui::TableNextColumn(); ImGui::Text("%.5g", I.misfit);
                    ImGui::TableNextColumn(); ImGui::Text("%.5g", I.objective);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", I.reg);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.gradNorm);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.step);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.kMin);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.kMax);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.kMean);
                    ImGui::TableNextColumn(); ImGui::Text("%d", I.runs);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", I.seconds);
                } else {
                    const auto& I = H.history[size_t(k)];
                    ImGui::TableNextColumn(); ImGui::Text("%.5g", I.misfit);
                    ImGui::TableNextColumn(); ImGui::Text("%.3g", I.gradNorm);
                    ImGui::TableNextColumn(); ImGui::Text("%d", I.runs);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", I.seconds);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f %%", H.initialMisfit > 0 ? 100.0 * (1 - I.misfit / H.initialMisfit) : 0.0);
                }
            }
            ui::endTable();
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ final parameters / field maps
    if (ui::tabItem(field ? ICON_GRID_3X3 " Field maps" : ICON_SLIDERS_HORIZONTAL " Final parameters")) {
        if (!field) {
            const auto& P = study_.hmParams;
            const auto& vk = K < int(H.history.size()) ? H.history[size_t(K)].values : H.best;
            const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp;
            if (ui::beginTable("##adjpar", 8, tf)) {
                for (const char* c : {"Parameter", "Low", "High", "Prior (base)", "Iterate", "Best", "u best", "d misfit / du (iterate)"}) ui::tableColumn(c);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < P.size(); ++i) {
                    if (!P[i].active) continue;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::TextUnformatted(paramLabel(P[i], case_).c_str());
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", P[i].low);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", P[i].high);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", P[i].base);
                    ImGui::TableNextColumn(); ImGui::Text("%.4g", i < vk.size() ? vk[i] : NAN);
                    ImGui::TableNextColumn(); ImGui::TextColored(col4(Pal::Accent), "%.4g", i < H.best.size() ? H.best[i] : NAN);
                    const double u = i < H.best.size() ? paramToUnit(P[i], H.best[i]) : NAN;
                    ImGui::TableNextColumn();
                    ImGui::ProgressBar(float(std::clamp(std::isfinite(u) ? u : 0.0, 0.0, 1.0)), ImVec2(-FLT_MIN, 0), fmtv(u).c_str());
                    ImGui::TableNextColumn();
                    const auto& g = K < int(H.history.size()) ? H.history[size_t(K)].grad : std::vector<double>();
                    if (i < g.size() && std::isfinite(g[i])) ImGui::Text("%.4g", g[i]); else ImGui::TextDisabled("-");
                }
                ui::endTable();
            }
            if (ImGui::Button(ICON_CHECK " Apply the best parameters to the case")) {
                const auto Pc = P; const auto best = H.best;
                study_.deferred = [this, Pc, best]() { applyStudyValues(Pc, best, "Adjoint L-BFGS best match"); };
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(K >= int(H.history.size()));
            if (ImGui::Button("Apply the selected iterate")) {
                const auto Pc = P; const auto v = H.history[size_t(K)].values;
                study_.deferred = [this, Pc, v]() { applyStudyValues(Pc, v, "Adjoint L-BFGS iterate"); };
            }
            ImGui::EndDisabled();
            // prior vs best positions in the normalised box
            if (ImPlot::BeginPlot("##adjbox", ImVec2(-1, std::max(fs * 10, ImGui::GetContentRegionAvail().y)))) {
                std::vector<double> u0, ub, uk; std::vector<std::string> names;
                for (size_t i = 0; i < P.size(); ++i) {
                    if (!P[i].active) continue;
                    names.push_back(paramLabel(P[i], case_));
                    u0.push_back(paramToUnit(P[i], P[i].base));
                    uk.push_back(i < vk.size() ? paramToUnit(P[i], vk[i]) : 0);
                    ub.push_back(i < H.best.size() ? paramToUnit(P[i], H.best[i]) : 0);
                }
                std::vector<const char*> lab; for (auto& n : names) lab.push_back(n.c_str());
                std::vector<double> pos(names.size()); for (size_t i = 0; i < pos.size(); ++i) pos[i] = double(i);
                ImPlot::SetupAxes(nullptr, "u = position in [low, high]", ImPlotAxisFlags_AutoFit, 0);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Always);
                if (!pos.empty()) ImPlot::SetupAxisTicks(ImAxis_X1, pos.data(), int(pos.size()), lab.data());
                const int n = int(names.size());
                ImPlotSpec s0; s0.FillColor = ImVec4(0.6f, 0.6f, 0.65f, 1);
                ImPlot::PlotBars("prior", u0.data(), n, 0.25, -0.27, s0);
                ImPlotSpec sk; sk.FillColor = ImVec4(0.55f, 0.75f, 1.0f, 1);
                ImPlot::PlotBars("iterate", uk.data(), n, 0.25, 0.0, sk);
                ImPlotSpec sb; sb.FillColor = col4(Pal::Accent);
                ImPlot::PlotBars("best", ub.data(), n, 0.25, 0.27, sb);
                ImPlot::EndPlot();
            }
        } else {
            static const char* qn[] = {"ln k multiplier", "ln phi multiplier", "d misfit / d ln k  (prior)", "d misfit / d ln k  (final)",
                                       "d misfit / d ln phi  (prior)", "d misfit / d ln phi  (final)"};
            ImGui::SetNextItemWidth(fs * 14); ImGui::Combo("Map", &study_.adjFieldQty, qn, 6);
            ImGui::SameLine();
            const int nz = std::max(1, F.nz), nx = std::max(1, F.nx), ny = std::max(1, F.ny);
            study_.adjLayer = std::clamp(study_.adjLayer, 0, nz - 1);
            ImGui::SetNextItemWidth(fs * 10); ImGui::SliderInt("Layer (K)", &study_.adjLayer, 0, nz - 1, "%d");
            const int q = study_.adjFieldQty;
            const std::vector<float>* src = nullptr;
            if (q == 0 && K < int(F.iters.size())) src = &F.iters[size_t(K)].lnK;
            else if (q == 1 && K < int(F.iters.size())) src = &F.iters[size_t(K)].lnPhi;
            else if (q == 2) src = &F.gradK0; else if (q == 3) src = &F.gradK; else if (q == 4) src = &F.gradPhi0; else if (q == 5) src = &F.gradPhi;
            const size_t ncell = size_t(nx) * size_t(ny) * size_t(nz);
            if (!src || src->size() < ncell) ui::hint("This map is not available (quantity not matched in this run).");
            else {
                std::vector<double> v(size_t(nx) * size_t(ny), 0.0);
                double m = 0;
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i) {
                        const float x = (*src)[size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(study_.adjLayer))];
                        const double d = std::isfinite(x) ? double(x) : 0.0;
                        v[size_t(j) * size_t(nx) + size_t(i)] = d; m = std::max(m, std::fabs(d));
                    }
                // symmetric range over all layers of this quantity (comparable layers)
                double mAll = 0; for (size_t c = 0; c < ncell; ++c) if (std::isfinite((*src)[c])) mAll = std::max(mAll, double(std::fabs((*src)[c])));
                if (mAll <= 0) mAll = 1;
                const float avW = ImGui::GetContentRegionAvail().x;
                const float hh = std::max(fs * 14, ImGui::GetContentRegionAvail().y - fs * 3);
                ImPlot::PushColormap(ImPlotColormap_RdBu);
                if (ImPlot::BeginPlot("##adjmap", ImVec2(avW * 0.62f, hh), ImPlotFlags_Equal | ImPlotFlags_NoLegend)) {
                    ImPlot::SetupAxes("I", "J (row 1 at the top)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::PlotHeatmap("##hm", v.data(), ny, nx, -mAll, mAll, nullptr, ImPlotPoint(0, 0), ImPlotPoint(nx, ny));
                    if (ImPlot::IsPlotHovered()) {
                        const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                        const int i = int(std::floor(mp.x)), jr = int(std::floor(mp.y));
                        const int j = ny - 1 - jr;
                        if (i >= 0 && i < nx && j >= 0 && j < ny) {
                            const double d = v[size_t(j) * size_t(nx) + size_t(i)];
                            if (q <= 1) ImGui::SetTooltip("cell (%d, %d, %d)\nln m = %.4g   multiplier = %.4g", i + 1, j + 1, study_.adjLayer + 1, d, std::exp(d));
                            else ImGui::SetTooltip("cell (%d, %d, %d)\n%s = %.4g", i + 1, j + 1, study_.adjLayer + 1, qn[q], d);
                        }
                    }
                    ImPlot::EndPlot();
                }
                ImGui::SameLine();
                ImPlot::ColormapScale(q <= 1 ? "ln multiplier" : "gradient", -mAll, mAll, ImVec2(fs * 5, hh), "%.3g");
                ImPlot::PopColormap();
                ImGui::SameLine();
                ImGui::BeginGroup();
                // histogram of the active-cell values
                std::vector<double> all;
                for (size_t c = 0; c < ncell; ++c) if (std::isfinite((*src)[c]) && (*src)[c] != 0.0f) all.push_back((*src)[c]);
                if (!all.empty() && ImPlot::BeginPlot("Distribution (all layers)", ImVec2(-1, hh * 0.6f), ImPlotFlags_NoLegend)) {
                    ImPlot::SetupAxes(q <= 1 ? "ln multiplier" : "gradient", "cells", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlotSpec hs; hs.FillColor = col4(Pal::Accent);
                    ImPlot::PlotHistogram("##h", all.data(), int(all.size()), 30, 1.0, ImPlotRange(), hs);
                    ImPlot::EndPlot();
                }
                ImGui::TextDisabled("layer max |value| %.3g, all layers %.3g", m, mAll);
                ImGui::TextDisabled("Blue = increase (multiplier > 1 /\npositive gradient), red = decrease");
                if (ImGui::Button(ICON_ROTATE_3D " Show in the 3D viewer")) {
                    std::vector<float> mk(ncell, NAN), mp(ncell, NAN), sk(ncell, NAN);
                    const auto& I = F.iters[size_t(std::min(K, int(F.iters.size()) - 1))];
                    for (size_t c = 0; c < ncell; ++c) {
                        if (c < I.lnK.size() && std::isfinite(I.lnK[c])) mk[c] = std::exp(I.lnK[c]);
                        if (c < I.lnPhi.size() && std::isfinite(I.lnPhi[c])) mp[c] = std::exp(I.lnPhi[c]);
                        const auto& G = F.gradK.empty() ? F.gradK0 : F.gradK;
                        if (c < G.size()) sk[c] = G[c];
                    }
                    if (I.lnK.empty()) mk.clear();
                    if (I.lnPhi.empty()) mp.clear();
                    if (F.gradK.empty() && F.gradK0.empty()) sk.clear();
                    const bool wantPhi = q == 1 || q >= 4, wantSens = q >= 2;
                    study_.deferred = [this, mk, mp, sk, wantPhi, wantSens]() {
                        hmMultK_ = mk; hmMultPhi_ = mp; hmSensK_ = sk; ++hmViewVer_;
                        const char* want = wantSens ? "HMSENSK" : wantPhi ? "HMMULTPHI" : "HMMULTK";
                        { const int pi = propIndexByName(want); if (pi >= 0) prop_ = pi; }
                        autoRange_ = true;
                        workspace_ = Workspace::Model;
                        notify("Adjoint HM map shown in the Model 3D viewer (property HMMULTK / HMMULTPHI / HMSENSK)");
                    };
                }
                if (ImGui::Button(ICON_CHECK " Apply the matched field to the case")) {
                    const Case best = F.best;
                    study_.deferred = [this, best]() {
                        case_.grid.permx = best.grid.permx; case_.grid.permy = best.grid.permy;
                        case_.grid.permz = best.grid.permz; case_.grid.poro = best.grid.poro;
                        markDirty(true);
                        notify("Matched per-cell permeability / porosity applied to the case - re-run the simulation (F5)");
                    };
                }
                ImGui::EndGroup();
            }
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ responses of the best case
    if (ui::tabItem(ICON_TARGET " Responses")) {
        const auto& O = case_.observations;
        const auto groups = obsGroups(case_);
        if (O.empty() || groups.empty()) ui::hint("No observed data in the case.");
        else {
            std::vector<const RunResult*> runs = {r0, rk, rb};
            std::vector<std::string> names = {"Prior", "Iteration " + std::to_string(K), "Best"};
            // normalised RMS per series (cached per selection)
            struct Cache { const void* a = nullptr; const void* b = nullptr; const void* c = nullptr; size_t nobs = 0; std::vector<std::vector<double>> rms; };
            static Cache C;
            if (C.a != r0 || C.b != rk || C.c != rb || C.nobs != O.size()) {
                C = Cache(); C.a = r0; C.b = rk; C.c = rb; C.nobs = O.size();
                const auto sig = autoSigma(O);
                C.rms.assign(3, std::vector<double>(groups.size(), NAN));
                for (size_t r = 0; r < 3; ++r) {
                    if (!runs[r] || !runs[r]->ok) continue;
                    const auto res = residuals(*runs[r], case_, O, sig);
                    for (size_t g = 0; g < groups.size(); ++g) {
                        double s = 0; int n = 0;
                        for (int i : groups[g].idx) if (O[size_t(i)].weight > 0 && size_t(i) < res.size() && std::isfinite(res[size_t(i)])) { s += res[size_t(i)] * res[size_t(i)]; ++n; }
                        if (n) C.rms[r][g] = std::sqrt(s / n);
                    }
                }
            }
            ImGui::Checkbox("All responses (grid)", &study_.adjAllResp);
            ImGui::SameLine(); ImGui::TextDisabled("normalised RMS = sqrt(mean((sim - obs) / sigma)^2): < 1 within the measurement error");
            if (!study_.adjAllResp) {
                const float tw = std::min(fs * 30, ImGui::GetContentRegionAvail().x * 0.42f);
                ImGui::BeginChild("##adjrl", ImVec2(tw, 0));
                const ImGuiTableFlags tf = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp;
                if (ui::beginTable("##adjresp", 4, tf, ImVec2(0, -FLT_MIN))) {
                    ImGui::TableSetupScrollFreeze(1, 1);
                    ui::tableColumn("Response"); ui::tableColumn("Prior"); ui::tableColumn("Iter"); ui::tableColumn("Best");
                    ImGui::TableHeadersRow();
                    for (size_t g = 0; g < groups.size(); ++g) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        if (ImGui::Selectable((groupName(groups[g]) + "##rg" + std::to_string(g)).c_str(), study_.adjResp == int(g), ImGuiSelectableFlags_SpanAllColumns)) study_.adjResp = int(g);
                        for (int r = 0; r < 3; ++r) {
                            ImGui::TableNextColumn();
                            const double v = C.rms[size_t(r)][g];
                            ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, misfitColour(v));
                            if (std::isfinite(v)) ImGui::Text("%.3g", v); else ImGui::TextDisabled("-");
                        }
                    }
                    ui::endTable();
                }
                ImGui::EndChild();
                ImGui::SameLine();
                ImGui::BeginChild("##adjrp");
                study_.adjResp = std::clamp(study_.adjResp, 0, int(groups.size()) - 1);
                std::vector<const RunResult*> pr; std::vector<const char*> lab;
                for (size_t r = 0; r < 3; ++r) if (runs[r] && runs[r]->ok) { pr.push_back(runs[r]); lab.push_back(names[r].c_str()); }
                char id[160]; std::snprintf(id, sizeof id, "%s##adjfit", groupName(groups[size_t(study_.adjResp)]).c_str());
                plotObservationFit(id, case_, pr, lab, study_.adjResp, std::max(fs * 12, ImGui::GetContentRegionAvail().y - 4));
                ImGui::EndChild();
            } else {
                const int ng = int(groups.size());
                const int cols = std::max(1, std::min(4, int(ImGui::GetContentRegionAvail().x / (fs * 22))));
                ImGui::BeginChild("##adjgrid");
                if (ImGui::BeginTable("##adjg", cols, ImGuiTableFlags_SizingStretchSame)) {
                    std::vector<const RunResult*> pr; std::vector<const char*> lab;
                    for (size_t r = 0; r < 3; ++r) if (runs[r] && runs[r]->ok) { pr.push_back(runs[r]); lab.push_back(names[r].c_str()); }
                    for (int g = 0; g < ng; ++g) {
                        ImGui::TableNextColumn();
                        char id[160]; std::snprintf(id, sizeof id, "%s  (RMS %.2g -> %.2g)##adjg%d", groupName(groups[size_t(g)]).c_str(), C.rms[0][size_t(g)], C.rms[2][size_t(g)], g);
                        plotObservationFit(id, case_, pr, lab, g, fs * 13);
                    }
                    ImGui::EndTable();
                }
                ImGui::EndChild();
            }
        }
        ImGui::EndTabItem();
    }

    // ------------------------------------------------------------------ sensitivity
    if (ui::tabItem(ICON_CHART_COLUMN " Sensitivity")) {
        if (!field) {
            const auto& P = study_.hmParams;
            std::vector<std::string> names; std::vector<double> g0, gk;
            for (size_t i = 0; i < P.size(); ++i) {
                if (!P[i].active) continue;
                names.push_back(paramLabel(P[i], case_));
                const auto& a = H.history.empty() ? std::vector<double>() : H.history[0].grad;
                const auto& b = K < int(H.history.size()) ? H.history[size_t(K)].grad : std::vector<double>();
                g0.push_back(i < a.size() && std::isfinite(a[i]) ? a[i] : 0.0);
                gk.push_back(i < b.size() && std::isfinite(b[i]) ? b[i] : 0.0);
            }
            ui::hint("Adjoint gradient of the misfit with respect to the normalised parameters u in [0, 1]: large bars = the data\n"
                     "constrain the parameter strongly; the sign gives the direction that INCREASES the misfit.");
            if (!names.empty() && ImPlot::BeginPlot("##adjsens", ImVec2(-1, -1))) {
                std::vector<const char*> lab; for (auto& n : names) lab.push_back(n.c_str());
                std::vector<double> pos(names.size()); for (size_t i = 0; i < pos.size(); ++i) pos[i] = double(i);
                ImPlot::SetupAxes(nullptr, "d misfit / du", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisTicks(ImAxis_X1, pos.data(), int(pos.size()), lab.data());
                ImPlotSpec s0; s0.FillColor = ImVec4(0.6f, 0.6f, 0.65f, 1);
                ImPlot::PlotBars("prior (iteration 0)", g0.data(), int(g0.size()), 0.35, -0.19, s0);
                ImPlotSpec sk; sk.FillColor = col4(Pal::Accent);
                char l[48]; std::snprintf(l, sizeof l, "iteration %d", K);
                ImPlot::PlotBars(l, gk.data(), int(gk.size()), 0.35, 0.19, sk);
                ImPlot::EndPlot();
            }
        } else {
            ui::hint("Per-cell sensitivity maps are in 'Field maps' (Map: d misfit / d ln k, prior or final).\n"
                     "Below: the layer-wise sum of |d misfit / d ln k| - where in the reservoir the data are informative.");
            const int nz = std::max(1, F.nz), nxy = std::max(1, F.nx * F.ny);
            std::vector<double> s0(static_cast<size_t>(nz), 0.0), s1(static_cast<size_t>(nz), 0.0), kk(static_cast<size_t>(nz), 0.0);
            for (int k = 0; k < nz; ++k) {
                kk[size_t(k)] = k + 1;
                for (int c = 0; c < nxy; ++c) {
                    const size_t q = size_t(c) + size_t(nxy) * size_t(k);
                    if (q < F.gradK0.size() && std::isfinite(F.gradK0[q])) s0[size_t(k)] += std::fabs(F.gradK0[q]);
                    if (q < F.gradK.size() && std::isfinite(F.gradK[q])) s1[size_t(k)] += std::fabs(F.gradK[q]);
                }
            }
            if (ImPlot::BeginPlot("##adjlay", ImVec2(-1, -1))) {
                ImPlot::SetupAxes("sum |d misfit / d ln k|", "layer", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                ImPlotSpec a; a.LineColor = ImVec4(0.6f, 0.6f, 0.65f, 1); a.LineWeight = 2; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 3;
                ImPlot::PlotLine("prior", s0.data(), kk.data(), nz, a);
                ImPlotSpec b; b.LineColor = col4(Pal::Accent); b.LineWeight = 2.5f; b.Marker = ImPlotMarker_Circle; b.MarkerSize = 3;
                ImPlot::PlotLine("final", s1.data(), kk.data(), nz, b);
                ImPlot::EndPlot();
            }
        }
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

// ============================================================================ C27 inversion statistics (iTOUGH2-style)
// Linearised a-posteriori analysis of the history-match result: parameter standard deviations / 95 % confidence
// intervals, correlation matrix, eigen-analysis, composite scaled sensitivities, residual analysis (studentized
// residuals, outliers, objective share per data set). Called with study_.m held (hmParameters view).
void App::hmStatsPanel() {
    const float fs = ImGui::GetFontSize();
    if (!ImGui::CollapsingHeader("Inversion statistics (iTOUGH2-style linear uncertainty analysis)")) return;
    const auto& H = study_.hm; const auto& P = study_.hmParams;
    int act = 0; for (auto& p : P) if (p.active) ++act;
    ImGui::BeginDisabled(study_.running || !H.ok || H.best.size() != P.size());
    char b[200];
    std::snprintf(b, sizeof b, ICON_PLAY " Compute at the matched values  (%d runs)", act + 1);
    if (ui::accentButton(b)) {
        const Case base = case_; const auto hp = P; const auto v = H.best; HMStatsOptions o; o.fdStep = study_.statsFd;
        o.eco = study_.eco; o.predict.clear();
        for (int r = 0; r <= int(Response::MaxFailure); ++r) if ((study_.statsPredMask >> r) & 1u) o.predict.push_back(Response(r));
        study_.deferred = [this, base, hp, v, o] {
            startStudy("Inversion statistics", [base, hp, v, o](StudyState& S, const StudyCallbacks& cb) {
                HMStatistics st = hmStatistics(base, hp, v, o, cb);
                std::lock_guard<std::mutex> l(S.m); S.hmStats = std::move(st); S.hasStats = true;
            });
        };
    }
    ImGui::EndDisabled();
    ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6);
    ImGui::SliderFloat("FD step (fraction of range)", &study_.statsFd, 0.005f, 0.1f, "%.3f");
    // C28 FOSM predictions (evaluated on the same runs: no extra simulations)
    ImGui::AlignTextToFramePadding(); ImGui::TextUnformatted("FOSM predictions:");
    for (int r = 0; r <= int(Response::MaxFailure); ++r) {
        if (Response(r) == Response::Misfit) continue;
        bool on = (study_.statsPredMask >> r) & 1u;
        ImGui::SameLine();
        if (ImGui::Checkbox(toString(Response(r)), &on)) study_.statsPredMask = on ? (study_.statsPredMask | (1u << r)) : (study_.statsPredMask & ~(1u << r));
    }
    if (!study_.hasStats) {
        ui::hint("C_pp = s0^2 (J^T J)^-1 with J the sensitivities of the weighted residuals (sim - obs) / sigma in the search space "
                 "(log10 for log parameters) and s0^2 = SSR / (m - n). Shows how well the data determine each parameter, which "
                 "parameters trade off (correlation), and which observations dominate or look like outliers.");
        return;
    }
    const HMStatistics& S = study_.hmStats;
    if (!S.ok) { ImGui::TextColored(col4(Pal::Err), "%s", S.error.c_str()); return; }
    ImGui::Text("%d observations, %d parameters, objective %.5g, s0^2 = %.4g, condition number of J^T J = %.3g", S.nObs, S.nPar, S.objective, S.s02, S.condition);
    if (S.condition > 1e8) ImGui::TextColored(col4(Pal::Warn), "Nearly singular: at least one parameter combination is not determined by the data (see the eigenvectors).");
    auto fmtv = [](double v) { char t[48]; if (!std::isfinite(v)) std::snprintf(t, sizeof t, "inf"); else std::snprintf(t, sizeof t, "%.4g", v); return std::string(t); };
    // ---- parameters
    double cssMax = 0; for (double c : S.css) cssMax = std::max(cssMax, c);
    if (ui::beginTable("hmst", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ui::tableColumn("Parameter", 0, 2.2f); ui::tableColumn("Value", 0, 0.9f); ui::tableColumn("Std. dev. (search coord.)", 0, 1.3f);
        ui::tableColumn("95 % interval", 0, 1.6f); ui::tableColumn("Composite scaled sensitivity", 0, 1.6f); ui::tableColumn("Determined?", 0, 1.0f);
        ImGui::TableHeadersRow();
        for (int k = 0; k < S.nPar; ++k) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(S.name[size_t(k)].c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(fmtv(S.value[size_t(k)]).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::Text("%s%s", fmtv(S.stdDev[size_t(k)]).c_str(), S.logSpace[size_t(k)] ? " (log10)" : "");
            ImGui::TableSetColumnIndex(3); ImGui::Text("%s .. %s", fmtv(S.ci95lo[size_t(k)]).c_str(), fmtv(S.ci95hi[size_t(k)]).c_str());
            ImGui::TableSetColumnIndex(4);
            {
                ImDrawList* dl = ImGui::GetWindowDrawList();
                const ImVec2 a = ImGui::GetCursorScreenPos(); const float w = ImGui::GetContentRegionAvail().x * 0.6f, h = fs * 0.8f;
                const float f = cssMax > 0 ? float(S.css[size_t(k)] / cssMax) : 0.f;
                dl->AddRectFilled(ImVec2(a.x, a.y + 2), ImVec2(a.x + w * f, a.y + 2 + h), Pal::Accent, 2);
                ImGui::Dummy(ImVec2(w, h + 2)); ImGui::SameLine(); ImGui::TextUnformatted(fmtv(S.css[size_t(k)]).c_str());
            }
            ImGui::TableSetColumnIndex(5);
            // iTOUGH2 / Hill rule of thumb: CSS < 1 % of the largest -> not estimable; sigma (log10) > 0.5 -> poorly determined
            const bool weak = cssMax > 0 && S.css[size_t(k)] < 0.01 * cssMax;
            const bool poor = !std::isfinite(S.stdDev[size_t(k)]) || (S.logSpace[size_t(k)] ? S.stdDev[size_t(k)] > 0.5 : (S.value[size_t(k)] != 0 && S.stdDev[size_t(k)] > 0.5 * std::fabs(S.value[size_t(k)])));
            if (weak) ImGui::TextColored(col4(Pal::Err), "insensitive");
            else if (poor) ImGui::TextColored(col4(Pal::Warn), "poorly");
            else ImGui::TextColored(col4(Pal::Ok), "well");
        }
        ui::endTable();
    }
    // ---- correlation matrix (colour = |r|) + eigen-analysis
    if (S.nPar > 1 && ImGui::TreeNodeEx("Correlation matrix and eigen-analysis", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginTable("hmcorr", S.nPar + 1, ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit)) {
            ui::tableColumn("", 0, 0);
            for (int k = 0; k < S.nPar; ++k) { std::snprintf(b, sizeof b, "p%d", k + 1); ui::tableColumn(b, 0, 0); }
            ImGui::TableHeadersRow();
            for (int a = 0; a < S.nPar; ++a) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::Text("p%d %s", a + 1, S.name[size_t(a)].substr(0, 22).c_str());
                for (int c = 0; c < S.nPar; ++c) {
                    ImGui::TableSetColumnIndex(c + 1);
                    const double r = S.corr[size_t(a)][size_t(c)];
                    const float t = float(std::min(1.0, std::fabs(r)));
                    const ImU32 bg = r >= 0 ? IM_COL32(220, 60, 60, int(40 + 160 * t)) : IM_COL32(60, 110, 230, int(40 + 160 * t));
                    if (a != c) ImGui::TableSetBgColor(ImGuiTableBgTarget_CellBg, bg);
                    ImGui::Text("%+.2f", r);
                }
            }
            ui::endTable();
        }
        ImGui::TextDisabled("|r| > 0.9: the two parameters trade off - the data constrain only a combination of them.");
        for (size_t k = 0; k < S.eigVal.size() && k < 3; ++k) {
            std::string v;
            for (int a = 0; a < S.nPar; ++a) { std::snprintf(b, sizeof b, "%s%+.2f p%d", a ? "  " : "", S.eigVec[k][size_t(a)], a + 1); v += b; }
            ImGui::Text("eigenvalue %s (std %s):  %s", fmtv(S.eigVal[k]).c_str(), fmtv(std::sqrt(S.eigVal[k])).c_str(), v.c_str());
        }
        ImGui::TextDisabled("largest covariance eigenvalues first = least determined parameter combinations");
        ImGui::TreePop();
    }
    // ---- residual analysis
    if (ImGui::TreeNodeEx("Residual analysis", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("weighted residuals: mean %.3g, std %.3g, skewness %.3g, excess kurtosis %.3g;  %zu potential outliers (|studentized| > 1.96)",
                    S.resMean, S.resStd, S.resSkew, S.resKurt, S.outliers.size());
        if (std::fabs(S.resMean) > 0.5 * std::max(S.resStd, 1e-12)) ImGui::TextColored(col4(Pal::Warn), "The residuals are biased: a systematic model error (missing physics / wrong fixed parameter) is likely.");
        if (ui::beginTable("hmcof", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("Data set", 0, 2.0f); ui::tableColumn("Share of the objective", 0, 2.0f);
            ImGui::TableHeadersRow();
            for (const auto& c : S.contribution) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(c.first.c_str());
                ImGui::TableSetColumnIndex(1); ImGui::ProgressBar(float(c.second / 100.0), ImVec2(-1, fs * 0.9f), (fmtv(c.second) + " %").c_str());
            }
            ui::endTable();
        }
        if (!S.studentized.empty() && ImPlot::BeginPlot("Studentized residuals", ImVec2(-1, fs * 9))) {
            ImPlot::SetupAxes("observation", "w_i", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            std::vector<double> x(S.studentized.size()); for (size_t i = 0; i < x.size(); ++i) x[i] = double(i + 1);
            ImPlotSpec s; s.Marker = ImPlotMarker_Circle; s.MarkerSize = 3; s.LineWeight = 0; s.MarkerFillColor = col4(Pal::Accent);
            ImPlot::PlotScatter("w", x.data(), S.studentized.data(), int(x.size()), s);
            const double lim[2] = {-1.96, 1.96};
            ImPlotSpec sl; sl.LineColor = col4(Pal::Warn); sl.LineWeight = 1.5f; sl.Flags = ImPlotInfLinesFlags_Horizontal;
            ImPlot::PlotInfLines("+-1.96", lim, 2, sl);
            ImPlot::EndPlot();
        }
        ImGui::TreePop();
    }
    // ---- C28 FOSM prediction uncertainty + data worth
    if (!S.pred.empty() && ImGui::TreeNodeEx("Prediction uncertainty (FOSM) and data worth", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ui::beginTable("hmpred", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ui::tableColumn("Prediction", 0, 1.2f); ui::tableColumn("Value +- 1 sd", 0, 1.6f); ui::tableColumn("95 % interval (linearised)", 0, 2.0f);
            ImGui::TableHeadersRow();
            for (const auto& p : S.pred) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(p.name.c_str());
                ImGui::TableSetColumnIndex(1); ImGui::Text("%s +- %s", fmtv(p.value).c_str(), fmtv(p.stdDev).c_str());
                ImGui::TableSetColumnIndex(2); ImGui::Text("%s .. %s", fmtv(p.lo).c_str(), fmtv(p.hi).c_str());
            }
            ui::endTable();
        }
        if (!S.worth.empty() && ui::beginTable("hmworth", 3 + int(S.pred.size()), ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY, ImVec2(0, fs * 14))) {
            ui::tableColumn("Data set", 0, 1.6f); ui::tableColumn("Kind", 0, 1.3f); ui::tableColumn("Points", 0, 0.6f);
            for (const auto& p : S.pred) ui::tableColumn(("d sd " + p.name + " [%]").c_str(), 0, 1.0f);
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableHeadersRow();
            for (const auto& w : S.worth) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(w.set.c_str());
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(w.existing ? "observed (remove)" : "candidate (add)");
                ImGui::TableSetColumnIndex(2); ImGui::Text("%d", w.points);
                for (size_t j = 0; j < w.dStdPct.size(); ++j) {
                    ImGui::TableSetColumnIndex(3 + int(j));
                    const double d = w.dStdPct[j];
                    const bool big = std::isfinite(d) && std::fabs(d) > 10;
                    if (big) ImGui::TextColored(col4(w.existing ? Pal::Warn : Pal::Ok), "%+.1f", d); else ImGui::TextUnformatted(std::isfinite(d) ? fmtv(d).c_str() : "unbounded");
                }
            }
            ui::endTable();
        }
        ui::hint("sigma_z^2 = g^T C_pp g with g = dz/dx from the Jacobian runs. Observed sets: how much the prediction uncertainty would grow "
                 "without that data (its worth). Candidate sets: expected reduction if that vector were measured at the observation times "
                 "with the automatic error (5 % of the value, >= 3 % of the series maximum). Linear analysis - valid near the match.");
        ImGui::TreePop();
    }
}
