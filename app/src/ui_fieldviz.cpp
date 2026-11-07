// C30: field / section / fracture visualisation of every run type and post-processed geomechanics
//  * fieldsTab(): Results tab "Fields & sections" (black oil, compositional, thermal, MPI / CLI runs, loaded .rsres):
//      Section map       any result array (pressure, saturations, every engine array: geomechanics stress / strain /
//                        displacement / failure, fracture stress shadow HFDS*, phase-field damage PF*, tracers, salt...) on
//                        a plan layer or a vertical I / J section, 15 colour maps, fixed range over all report steps or a
//                        symmetric range around 0 for changes, fracture / well overlays, hover read-out
//      Fracture face     the hydraulic / EDFM fracture in its own plane (distance along strike x depth), coloured by any
//                        fracture-cell array (HFWIDTH, PFWIDTH, HFPNET, HFCOND, HFOPEN, pressure, saturation; thermal EDFM:
//                        aperture / pressure / temperature / conductivity)
//      Opening profiles  w(s) (or any fracture array) along the strike for every report step - how the opening develops
//      Phase field       the in-run phase-field model: damage 1 - v, stress components, displacement on the model mesh,
//                        crack-opening profile w(x) of every step
//      Well section      C31: curtain under the well trajectory coloured by any array + per-connection log tracks (TVD / MD)
//      Geomechanics      one-way FE post-processing (geoseq::postProcess) of the finished run: stress / strain /
//                        displacement / failure maps from the pressure (and temperature) changes, for runs without GEOMECH
//  * postGeoControls(): the post-processing form (also used for ECLIPSE / OPM results in the result viewer)
//  * pollPostGeo(): merges a finished post-processing job into its results
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "plotkit.h"
#include "implot.h"
#include "colormaps.h"
#include "ressim/geoseq.h"
#include "ressim/phasefield.h"
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <set>
#include <thread>

using namespace ressim;

struct App::PostGeoJob {
    std::thread th;
    std::atomic<int> done{0}, total{0};
    std::atomic<bool> finished{false}, cancel{false};
    std::vector<Snapshot> snaps;
    geoseq::PostResult res;
    int target = 0;                 // 0: main run results, 1: result-viewer case
    std::string caseKey;            // result-viewer case: name|path
    size_t nSnaps = 0;
    ~PostGeoJob() { cancel = true; if (th.joinable()) th.join(); }
};

// colour map of the 2-D field plots (sections, fracture faces, phase field, EGS / fracture panels in ui_egs.cpp)
Colormap g_fvMap = Colormap::Turbo;
bool g_fvRev = false;
namespace {

ImU32 cmU32(float t, float a = 1.0f) {
    float rgb[3]; colormapRGB(g_fvMap, t, rgb, g_fvRev);
    return ImGui::ColorConvertFloat4ToU32(ImVec4(rgb[0], rgb[1], rgb[2], a));
}
void cmapCombo(const char* id) {
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextItemWidth(fs * 8);
    if (ImGui::BeginCombo(id, colormapName(g_fvMap))) {
        for (int i = 0; i < int(Colormap::Count); ++i)
            if (ImGui::Selectable(colormapName(Colormap(i)), int(g_fvMap) == i)) g_fvMap = Colormap(i);
        ImGui::EndCombo();
    }
    ImGui::SameLine(); ImGui::Checkbox("Reverse##fvrev", &g_fvRev);
}
void colorBarH(const char* label, double lo, double hi, float width, const char* unit) {
    const float fs = ImGui::GetFontSize();
    ImGui::TextDisabled("%s", label); ImGui::SameLine();
    char b[48]; std::snprintf(b, sizeof(b), "%.4g", lo); ImGui::TextUnformatted(b); ImGui::SameLine();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < 48; ++i)
        dl->AddRectFilled(ImVec2(p.x + width * i / 48.0f, p.y + fs * 0.15f), ImVec2(p.x + width * (i + 1) / 48.0f + 0.5f, p.y + fs * 0.95f), cmU32((i + 0.5f) / 48.0f));
    ImGui::Dummy(ImVec2(width, fs)); ImGui::SameLine();
    std::snprintf(b, sizeof(b), "%.4g %s", hi, unit ? unit : ""); ImGui::TextUnformatted(b);
}
void exportScale(const char* label, double lo, double hi) {
    if (!pk::recording()) return;
    ImU32 lut[64];
    for (int k = 0; k < 64; ++k) lut[k] = cmU32(k / 63.0f);
    pk::colorScale(label, lo, hi, lut, 64);
}
// unit of a result array (case units)
std::string unitOfArray(const std::string& n, UnitSystem us) {
    const bool fld = us == UnitSystem::Field;
    const std::string p = fld ? "psi" : "bar";
    if (n == "PRESSURE" || n == "HFPNET" || n.rfind("HFDS", 0) == 0 || n == "HFDTAU" || n.rfind("PFDS", 0) == 0 || n.rfind("SIGEFF", 0) == 0 ||
        n == "SIGMIN" || n == "SQRTJ2") return p;
    if (n == "HFWIDTH" || n == "PFWIDTH" || n.rfind("HFDISP", 0) == 0 || n == "PFDISP") return "mm";
    if (n == "UPLIFT" || n == "DISPX" || n == "DISPY" || n == "DISPMAG") return "m";
    if (n == "HFSTRV") return "microstrain";
    if (n == "TEMP" || n == "DTEMP") return fld ? "degF" : "degC";
    return "";
}
bool signedArray(const std::string& n) {
    return n.rfind("HFDS", 0) == 0 || n.rfind("PFDS", 0) == 0 || n.rfind("HFDISP", 0) == 0 || n == "UPLIFT" || n == "DISPX" || n == "DISPY" ||
           n == "HFSTRV" || n == "STRAINV" || n == "DTEMP";
}
std::vector<std::string> arrayList(const RunResults& R) {
    std::vector<std::string> v = {"PRESSURE", "SWAT", "SOIL", "SGAS", "RS"};
    if (!R.snaps.empty()) for (const auto& e : R.snaps.back().extra) v.push_back(e.first);
    return v;
}
}  // namespace

// ---------------------------------------------------------------------------- post-processed geomechanics
void App::postGeoControls(int target, std::shared_ptr<const ressim::Grid> g, bool fieldUnits, const std::vector<ressim::Snapshot>* snaps) {
    const float fs = ImGui::GetFontSize();
    if (!postGeoInit_) {
        postGeoInit_ = true;
        postGeoIn_ = geomech::Input{};
        postGeoIn_.E = {10e9}; postGeoIn_.nu = {0.25}; postGeoIn_.biot = {1.0};
        postGeoIn_.overLayers = 4; postGeoIn_.underLayers = 2; postGeoIn_.sidePad = 2;
        postGeoIn_.overThick = 1000; postGeoIn_.underThick = 500; postGeoIn_.sideThick = 1000;
    }
    geomech::Input& in = postGeoIn_;
    ui::hint("One-way finite-element geomechanics of the finished run: the pressure (and temperature) change of every report "
             "step from the first one loads the mechanics mesh of the grid (Cartesian bricks or corner-point hexahedra), and the "
             "stress / strain / displacement / failure maps are added to the run - in the 3-D viewer (Engine output arrays), "
             "the Section map and the cell history. Works for every engine and for ECLIPSE / OPM restart results; for a fully "
             "coupled run use GEOMECH in the deck instead.");
    static double E = 10, nu = 0.25, biot = 1.0, alphaT = 1e-5, coh = 2, fric = 30, k0 = 0.8, sv = 22.6;
    static bool plastic = false;
    if (ui::beginProps("##pg", 0.5f)) {
        ui::propD("Young's modulus (GPa)", &E, "%.4g");
        ui::propD("Poisson's ratio", &nu, "%.3g");
        ui::propD("Biot coefficient", &biot, "%.3g");
        ui::propD("Thermal expansion (1/K)", &alphaT, "%.3g", "linear expansion coefficient; used when the run has TEMP");
        ui::propI("Overburden layers", &in.overLayers, 0, 40);
        ui::propD("Overburden thickness (m)", &in.overThick, "%.4g");
        ui::propI("Underburden layers", &in.underLayers, 0, 40);
        ui::propD("Underburden thickness (m)", &in.underThick, "%.4g");
        ui::propI("Side-burden elements", &in.sidePad, 0, 20);
        ui::propD("Side-burden width (m)", &in.sideThick, "%.4g");
        ui::propD("Overburden gradient (kPa/m)", &sv, "%.4g");
        ui::propD("K0 (sigma'h / sigma'v)", &k0, "%.3g");
        ui::propD("Cohesion (MPa)", &coh, "%.3g");
        ui::propD("Friction angle (deg)", &fric, "%.3g");
        ui::propBool("Drucker-Prager plasticity", &plastic);
        ui::endProps();
    }
    in.E = {std::max(1e-3, E) * 1e9}; in.nu = {std::clamp(nu, 0.0, 0.49)}; in.biot = {std::clamp(biot, 0.0, 1.0)}; in.alphaT = {alphaT};
    in.padE = in.E[0]; in.padNu = in.nu[0]; in.padAlphaT = alphaT;
    in.cohesion = coh * 1e6; in.frictionDeg = fric; in.k0 = k0; in.sigmaVGrad = sv * 1e3; in.plasticity = plastic ? 1 : 0;
    const bool busy = postGeo_ && !postGeo_->finished;
    if (busy) {
        const int d = postGeo_->done, t = std::max(1, int(postGeo_->total));
        ImGui::ProgressBar(float(d) / float(t), ImVec2(fs * 16, 0));
        ImGui::SameLine();
        if (ImGui::Button("Cancel##pg")) postGeo_->cancel = true;
        return;
    }
    if (!g) { ui::hint("No grid for this run."); return; }
    if (ImGui::Button(ICON_PLAY "  Compute geomechanics for this run")) {
        auto job = std::make_shared<PostGeoJob>();
        job->target = target;
        if (target == 0) {   // the caller holds the runner lock and passes the snapshots
            if (snaps) job->snaps = *snaps;
        } else if (rv_.c3d >= 0 && rv_.c3d < int(rv_.cases.size())) {
            RVCase& c = rv_.cases[size_t(rv_.c3d)];
            job->caseKey = c.name + "|" + c.path;
            if (c.rc) { try { job->snaps = ecl::snapshotsFromRestart(*c.rc, -1, 0); } catch (...) {} }
            else job->snaps = c.snaps;
        }
        job->nSnaps = job->snaps.size();
        if (job->snaps.empty()) { postGeoMsg_ = "no report steps with PRESSURE"; return; }
        const geomech::Input inCopy = in;
        PostGeoJob* J = job.get();
        job->th = std::thread([J, g, inCopy, fieldUnits]() {
            J->res = geoseq::postProcess(inCopy, *g, fieldUnits, J->snaps, [J](int i, int n) { J->done = i; J->total = n; return !J->cancel.load(); });
            J->finished = true;
        });
        postGeo_ = job;
    }
    if (!postGeoMsg_.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", postGeoMsg_.c_str()); }
}

void App::pollPostGeo() {
    if (!postGeo_ || !postGeo_->finished) return;
    std::shared_ptr<PostGeoJob> J = postGeo_;
    postGeo_.reset();
    if (J->th.joinable()) J->th.join();
    char b[256];
    if (!J->res.error.empty() && J->res.steps == 0) { postGeoMsg_ = "Geomechanics post-processing failed: " + J->res.error; notify(postGeoMsg_); return; }
    std::snprintf(b, sizeof(b), "Geomechanics: %d report steps in %.1f s, reservoir uplift max %.3g mm / subsidence max %.3g mm, max failure index %.3f%s",
                  J->res.steps, J->res.seconds, J->res.maxUplift * 1e3, J->res.maxSubsidence * 1e3, J->res.maxFailure,
                  J->res.error.empty() ? "" : (" (" + J->res.error + ")").c_str());
    postGeoMsg_ = b;
    static const std::set<std::string> mech = {"UPLIFT", "DISPX", "DISPY", "DISPMAG", "STRAINV", "SIGEFFM", "SIGEFFV", "SIGEFFH", "FAILURE",
                                               "POROM", "SIGMIN", "SQRTJ2", "EPSPLAS", "EPSPV", "KDAMAGE"};
    auto merge = [&](std::vector<Snapshot>& dst) {
        for (size_t s = 0; s < dst.size() && s < J->snaps.size(); ++s)
            for (const auto& e : J->snaps[s].extra) {
                if (!mech.count(e.first)) continue;
                bool rep = false;
                for (auto& d : dst[s].extra) if (d.first == e.first) { d.second = e.second; rep = true; break; }
                if (!rep) dst[s].extra.push_back(e);
            }
    };
    if (J->target == 0) {
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            RunResults& R = runner_.results();
            if (R.snaps.size() == J->nSnaps) merge(R.snaps);
        }
        runner_.touch();
    } else {
        for (auto& c : rv_.cases)
            if (c.name + "|" + c.path == J->caseKey) {
                if (c.rc) { c.post.resize(J->snaps.size()); merge(c.post); }
                else merge(c.snaps);
                rv_.key3d.clear(); rv_.grKey.clear();
            }
    }
    notify(postGeoMsg_);
}

// ---------------------------------------------------------------------------- Results tab
void App::fieldsTab(const RunResults& R) {
    const float fs = ImGui::GetFontSize();
    if (!R.grid || R.snaps.empty()) { ui::hint("This run has no grid / report steps."); return; }
    const Grid& g = *R.grid;
    const UnitSystem us = R.input.units;
    const bool fld = us == UnitSystem::Field;
    const char* lu = fld ? "ft" : "m";
    const int ns = int(R.snaps.size());
    static int step = -1;
    if (step < 0 || step >= ns) step = std::clamp(snap_, 0, ns - 1);
    const bool hasHf = !R.input.grid.hfPieces.empty() && g.numHfCells() == int(R.input.grid.hfPieces.size());
    const bool hasSeg = !R.fracSegs.empty() && R.fracSnap.size() == R.snaps.size();
    const bool hasPf = R.hfPf && std::any_of(R.hfPf->begin(), R.hfPf->end(), [](const Simulator::PfRecord& r) { return bool(r.field); });

    if (!ImGui::BeginTabBar("##fields")) return;
    // ======================================================== section map
    if (ui::tabItem("Section map")) {
        static std::string arr = "PRESSURE";
        static int orient = 0, idx = -1, rangeMode = 0;
        static bool showFrac = true, showWells = true;
        const std::vector<std::string> names = arrayList(R);
        if (std::find(names.begin(), names.end(), arr) == names.end()) arr = "PRESSURE";
        ImGui::SetNextItemWidth(fs * 10);
        if (ImGui::BeginCombo("Array##sm", arr.c_str(), ImGuiComboFlags_HeightLarge)) {
            for (const auto& n : names) {
                const std::string u = unitOfArray(n, us);
                if (ImGui::Selectable((n + (u.empty() ? "" : "  (" + u + ")")).c_str(), n == arr)) { arr = n; if (signedArray(n)) rangeMode = 2; }
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        static const char* ors[] = {"Plan (layer K)", "Section along X (row J)", "Section along Y (column I)"};
        ImGui::SetNextItemWidth(fs * 12); ImGui::Combo("##orient", &orient, ors, 3);
        const int nIdx = orient == 0 ? g.nz : orient == 1 ? g.ny : g.nx;
        if (idx < 0 || idx >= nIdx) idx = nIdx / 2;
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 9); ImGui::SliderInt(orient == 0 ? "K##idx" : orient == 1 ? "J##idx" : "I##idx", &idx, 0, std::max(0, nIdx - 1), "%d");
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 12); ImGui::SliderInt("Step##sm", &step, 0, ns - 1, "%d");
        ImGui::SameLine(); ImGui::TextDisabled("t = %.4g d", R.snaps[size_t(step)].time);
        cmapCombo("##cmsm");
        ImGui::SameLine();
        static const char* rms[] = {"Range: all steps", "Range: this step", "Range: symmetric around 0"};
        ImGui::SetNextItemWidth(fs * 13); ImGui::Combo("##rm", &rangeMode, rms, 3);
        ImGui::SameLine(); ImGui::Checkbox("Fractures", &showFrac);
        ImGui::SameLine(); ImGui::Checkbox("Wells", &showWells);
        const std::vector<float>& v = R.snaps[size_t(step)].field(arr);
        auto cellOf = [&](int a, int b) {   // section coordinates -> global cell
            return orient == 0 ? a + g.nx * (b + g.ny * idx) : orient == 1 ? a + g.nx * (idx + g.ny * b) : idx + g.nx * (a + g.ny * b);
        };
        const int na_ = orient == 0 ? g.nx : orient == 1 ? g.nx : g.ny, nb_ = orient == 0 ? g.ny : g.nz;
        double lo = 1e300, hi = -1e300;
        auto acc = [&](const std::vector<float>& a, bool onlySection) {
            if (onlySection) {
                for (int i = 0; i < na_; ++i) for (int j = 0; j < nb_; ++j) {
                    const int q = cellOf(i, j);
                    const int ac = size_t(q) < g.activeOf.size() ? g.activeOf[size_t(q)] : -1;
                    if (ac >= 0 && size_t(ac) < a.size() && std::isfinite(a[size_t(ac)])) { lo = std::min(lo, double(a[size_t(ac)])); hi = std::max(hi, double(a[size_t(ac)])); }
                }
            } else for (float x : a) if (std::isfinite(x)) { lo = std::min(lo, double(x)); hi = std::max(hi, double(x)); }
        };
        if (rangeMode == 1) acc(v, true);
        else {   // all steps: cached (Norne-size runs: millions of values per frame otherwise)
            static std::string rKey; static double rLo = 0, rHi = 0;
            char kb[112]; std::snprintf(kb, sizeof kb, "|%d|%.9g|%zu|%p|%u", ns, R.snaps.back().time, R.snaps.back().extra.size(), static_cast<const void*>(&R.snaps[0]), unsigned(runner_.version()));
            const std::string key = arr + kb;
            if (key != rKey) { for (int s = 0; s < ns; ++s) acc(R.snaps[size_t(s)].field(arr), false); rKey = key; rLo = lo; rHi = hi; }
            else { lo = rLo; hi = rHi; }
        }
        if (rangeMode == 2) { const double m = std::max(std::fabs(lo), std::fabs(hi)); lo = -m; hi = m; }
        if (!(hi >= lo)) { ui::hint("No values of this array."); ImGui::EndTabItem(); ImGui::EndTabBar(); return; }
        if (hi - lo < 1e-12) { hi = lo + 1; }
        const std::string unit = unitOfArray(arr, us);
        colorBarH(arr.c_str(), lo, hi, fs * 14, unit.c_str());
        fig::recorder().beginFigure("run.fields", 1, 1);
        if (ImPlot::BeginPlot("##secmap", ImVec2(-1, -1), (orient == 0 ? ImPlotFlags_Equal : 0) | ImPlotFlags_NoMenus)) {
            const std::string xl = std::string(orient == 2 ? "y (" : "x (") + lu + ")";
            const std::string yl = std::string(orient == 0 ? "y (" : "depth (") + lu + ")";
            recPanel((arr + (orient == 0 ? " - plan" : " - section")).c_str(), xl.c_str(), yl.c_str(), "");
            ImPlot::SetupAxes(xl.c_str(), yl.c_str(), 0, orient == 0 ? 0 : ImPlotAxisFlags_Invert);   // sections: depth down
            if (orient == 0) ImPlot::SetupAxesLimits(g.bbMin.x, g.bbMax.x, g.bbMin.y, g.bbMax.y, ImPlotCond_Once);
            else ImPlot::SetupAxesLimits(orient == 1 ? g.bbMin.x : g.bbMin.y, orient == 1 ? g.bbMax.x : g.bbMax.y, g.bbMin.z, g.bbMax.z, ImPlotCond_Once);
            const double bx[2] = {orient == 2 ? g.bbMin.y : g.bbMin.x, orient == 2 ? g.bbMax.y : g.bbMax.x};
            const double by[2] = {orient == 0 ? g.bbMin.y : g.bbMin.z, orient == 0 ? g.bbMax.y : g.bbMax.z};
            ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
            ImPlot::PlotScatter("##bb", bx, by, 2, sp);
            ImPlot::PushPlotClipRect();
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            int hovered = -1;
            for (int i = 0; i < na_; ++i)
                for (int j = 0; j < nb_; ++j) {
                    const int q = cellOf(i, j);
                    if (size_t(q) >= g.corners.size()) continue;
                    const auto& C = g.corners[size_t(q)];
                    double qx[4], qy[4];
                    if (orient == 0) { qx[0] = C[0].x; qy[0] = C[0].y; qx[1] = C[1].x; qy[1] = C[1].y; qx[2] = C[3].x; qy[2] = C[3].y; qx[3] = C[2].x; qy[3] = C[2].y; }
                    else if (orient == 1) { qx[0] = C[0].x; qy[0] = C[0].z; qx[1] = C[1].x; qy[1] = C[1].z; qx[2] = C[5].x; qy[2] = C[5].z; qx[3] = C[4].x; qy[3] = C[4].z; }
                    else { qx[0] = C[0].y; qy[0] = C[0].z; qx[1] = C[2].y; qy[1] = C[2].z; qx[2] = C[6].y; qy[2] = C[6].z; qx[3] = C[4].y; qy[3] = C[4].z; }
                    const int a = size_t(q) < g.activeOf.size() ? g.activeOf[size_t(q)] : -1;
                    ImU32 col = IM_COL32(70, 70, 76, 255);
                    if (a >= 0 && size_t(a) < v.size() && std::isfinite(v[size_t(a)])) col = cmU32(float((v[size_t(a)] - lo) / (hi - lo)));
                    else if (a < 0) continue;
                    pk::quad(qx, qy, col);
                    const double x0 = std::min({qx[0], qx[1], qx[2], qx[3]}), x1 = std::max({qx[0], qx[1], qx[2], qx[3]});
                    const double y0 = std::min({qy[0], qy[1], qy[2], qy[3]}), y1 = std::max({qy[0], qy[1], qy[2], qy[3]});
                    if (mp.x >= x0 && mp.x <= x1 && mp.y >= y0 && mp.y <= y1) hovered = q;
                }
            exportScale(arr.c_str(), lo, hi);
            if (showFrac) {
                const ImU32 fc = IM_COL32(255, 255, 255, 230);
                auto drawSeg = [&](double xa, double ya, double xb, double yb, double zt, double zb, int host) {
                    if (orient == 0) {
                        const int k = host >= 0 ? host / (g.nx * g.ny) : -1;
                        if (k == idx) pk::line(xa, ya, xb, yb, fc, 2.5f, "fracture");
                    } else {
                        const int j = host >= 0 ? (host / g.nx) % g.ny : -1, i = host >= 0 ? host % g.nx : -1;
                        if ((orient == 1 && j == idx) || (orient == 2 && i == idx)) {
                            const double s0 = orient == 1 ? std::min(xa, xb) : std::min(ya, yb), s1 = orient == 1 ? std::max(xa, xb) : std::max(ya, yb);
                            if (s1 - s0 < 1e-6 * (1 + std::fabs(s0))) pk::line(s0, zt, s0, zb, fc, 2.5f, "fracture");
                            else pk::rect(s0, zt, s1, zb, IM_COL32(255, 255, 255, 40), fc, 1.5f, "fracture");
                        }
                    }
                };
                if (hasHf) for (const auto& P : R.input.grid.hfPieces) drawSeg(P.xa, P.ya, P.xb, P.yb, P.zt, P.zb, P.host);
                else if (!R.fracSegs.empty()) for (const auto& s : R.fracSegs) drawSeg(s.xa, s.ya, s.xb, s.yb, s.zt, s.zb, s.host);
            }
            if (showWells)
                for (const WellInput& w : R.input.wells)
                    for (const Completion& c : w.completions) {
                        if (c.i < 0 || c.j < 0 || c.i >= g.nx || c.j >= g.ny) continue;
                        const ImU32 wc = w.type == WellType::Producer ? IM_COL32(60, 200, 110, 255) : IM_COL32(70, 140, 255, 255);
                        if (orient == 0) {
                            if (idx < c.k1 || idx > c.k2) continue;
                            const int q = c.i + g.nx * (c.j + g.ny * idx);
                            if (size_t(q) < g.center.size()) pk::marker(g.center[size_t(q)].x, g.center[size_t(q)].y, fs * 0.3f, wc, 0, 0, 1.0f, "well");
                        } else if ((orient == 1 && c.j == idx) || (orient == 2 && c.i == idx)) {
                            for (int k = std::max(0, c.k1); k <= std::min(g.nz - 1, c.k2); ++k) {
                                const int q = c.i + g.nx * (c.j + g.ny * k);
                                if (size_t(q) < g.center.size()) pk::marker(orient == 1 ? g.center[size_t(q)].x : g.center[size_t(q)].y, g.center[size_t(q)].z, fs * 0.28f, wc, 0, 0, 1.0f, "well");
                            }
                        }
                    }
            ImPlot::PopPlotClipRect();
            if (ImPlot::IsPlotHovered() && hovered >= 0) {
                const int a = g.activeOf[size_t(hovered)];
                const int I = hovered % g.nx, J = (hovered / g.nx) % g.ny, K = hovered / (g.nx * g.ny);
                if (a >= 0 && size_t(a) < v.size()) ImGui::SetTooltip("cell (%d, %d, %d)\n%s = %.5g %s", I + 1, J + 1, K + 1, arr.c_str(), v[size_t(a)], unit.c_str());
            }
            recEnd();
            ImPlot::EndPlot();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    // ======================================================== fracture face / opening profiles
    if ((hasHf || hasSeg) && ui::tabItem("Fracture face")) {
        static int frac = 0, qsel = 0;
        // fracture list and per-piece geometry (strike distance from the trace start, depth range)
        std::vector<std::string> fn;
        struct Pc { int frac; double s0, s1, zt, zb; int act; size_t seg; };
        std::vector<Pc> pcs;
        if (hasHf) {
            for (const auto& F : R.input.grid.hfracs) fn.push_back(F.name);
            for (size_t k = 0; k < R.input.grid.hfPieces.size(); ++k) {
                const auto& P = R.input.grid.hfPieces[k];
                const auto& F = R.input.grid.hfracs[size_t(P.frac)];
                const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0), tx = L > 0 ? (F.x1 - F.x0) / L : 1, ty = L > 0 ? (F.y1 - F.y0) / L : 0;
                const double sa = (P.xa - F.x0) * tx + (P.ya - F.y0) * ty, sb = (P.xb - F.x0) * tx + (P.yb - F.y0) * ty;
                const int a = g.activeOf[size_t(g.hfCells[k].global)];
                pcs.push_back({P.frac, std::min(sa, sb), std::max(sa, sb), P.zt, P.zb, a, k});
            }
        } else {
            fn = R.fracNames;
            std::vector<double> x0(fn.size() + 1, NAN), y0(fn.size() + 1, NAN);
            for (size_t i = 0; i < R.fracSegs.size(); ++i) {
                const auto& s = R.fracSegs[i];
                const size_t f = size_t(std::max(0, s.frac));
                if (f < x0.size() && !std::isfinite(x0[f])) { x0[f] = s.xa; y0[f] = s.ya; }
            }
            for (size_t i = 0; i < R.fracSegs.size(); ++i) {
                const auto& s = R.fracSegs[i];
                const size_t f = size_t(std::max(0, s.frac));
                const double sa = std::hypot(s.xa - x0[f], s.ya - y0[f]), sb = std::hypot(s.xb - x0[f], s.yb - y0[f]);
                pcs.push_back({s.frac, std::min(sa, sb), std::max(sa, sb), s.zt, s.zb, -1, i});
            }
        }
        if (fn.empty()) fn.push_back("fracture");
        frac = std::clamp(frac, 0, int(fn.size()) - 1);
        ImGui::SetNextItemWidth(fs * 9);
        if (ImGui::BeginCombo("Fracture##ff", fn[size_t(frac)].c_str())) {
            for (int i = 0; i < int(fn.size()); ++i) if (ImGui::Selectable(fn[size_t(i)].c_str(), i == frac)) frac = i;
            ImGui::EndCombo();
        }
        // quantities with values on the fracture cells
        std::vector<std::string> qn;
        if (hasHf) {
            qn = {"PRESSURE", "SWAT", "SGAS"};
            for (const auto& e : R.snaps.back().extra) {
                bool any = false;
                for (const auto& p : pcs) if (p.act >= 0 && size_t(p.act) < e.second.size() && e.second[size_t(p.act)] != 0) { any = true; break; }
                if (any) qn.push_back(e.first);
            }
        } else qn = {"Aperture (mm)", "Pressure", "Temperature", "Conductivity (mD m)"};
        if (qsel >= int(qn.size())) qsel = 0;
        if (hasHf && qsel == 0) for (int i = 0; i < int(qn.size()); ++i) if (qn[size_t(i)] == "HFWIDTH" || qn[size_t(i)] == "PFWIDTH") { qsel = i; break; }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 10);
        if (ImGui::BeginCombo("Quantity##ff", qn[size_t(qsel)].c_str())) {
            for (int i = 0; i < int(qn.size()); ++i) if (ImGui::Selectable(qn[size_t(i)].c_str(), i == qsel)) qsel = i;
            ImGui::EndCombo();
        }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 12); ImGui::SliderInt("Step##ff", &step, 0, ns - 1, "%d");
        ImGui::SameLine(); ImGui::TextDisabled("t = %.4g d", R.snaps[size_t(step)].time);
        cmapCombo("##cmff");
        auto valueAt = [&](const Pc& p, int s) -> double {
            if (hasHf) {
                const std::vector<float>& a = R.snaps[size_t(s)].field(qn[size_t(qsel)]);
                return p.act >= 0 && size_t(p.act) < a.size() ? double(a[size_t(p.act)]) : NAN;
            }
            const auto& st = R.fracSnap[size_t(s)];
            const std::vector<float>& a = qsel == 0 ? st.ap : qsel == 1 ? st.p : qsel == 2 ? st.T : st.cond;
            return p.seg < a.size() ? double(a[p.seg]) : NAN;
        };
        double lo = 1e300, hi = -1e300;
        for (int s = 0; s < ns; ++s)
            for (const auto& p : pcs) if (p.frac == frac) { const double x = valueAt(p, s); if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); } }
        if (!(hi >= lo)) { lo = 0; hi = 1; }
        if (hi - lo < 1e-12) hi = lo + 1;
        const std::string unit = hasHf ? unitOfArray(qn[size_t(qsel)], us) : std::string();
        colorBarH(qn[size_t(qsel)].c_str(), lo, hi, fs * 14, unit.c_str());
        const ImVec2 av = ImGui::GetContentRegionAvail();
        fig::recorder().beginFigure("run.fracface", 2, 1);
        if (ImPlot::BeginSubplots("##ffsub", 2, 1, ImVec2(av.x, av.y), 0)) {
            if (ImPlot::BeginPlot("Fracture face")) {
                const std::string xl = std::string("distance along strike (") + lu + ")", yl = std::string("depth (") + lu + ")";
                recPanel(("Fracture face: " + qn[size_t(qsel)]).c_str(), xl.c_str(), yl.c_str(), "");
                ImPlot::SetupAxes(xl.c_str(), yl.c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
                double s0 = 1e300, s1 = -1e300, z0 = 1e300, z1 = -1e300;
                for (const auto& p : pcs) if (p.frac == frac) { s0 = std::min(s0, p.s0); s1 = std::max(s1, p.s1); z0 = std::min(z0, p.zt); z1 = std::max(z1, p.zb); }
                if (s1 > s0) {
                    const double bx[2] = {s0, s1}, by[2] = {z0, z1};
                    ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
                    ImPlot::PlotScatter("##bbf", bx, by, 2, sp);
                }
                ImPlot::PushPlotClipRect();
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                const Pc* hov = nullptr;
                for (const auto& p : pcs) {
                    if (p.frac != frac) continue;
                    const double x = valueAt(p, step);
                    const ImU32 col = std::isfinite(x) ? cmU32(float((x - lo) / (hi - lo))) : IM_COL32(70, 70, 76, 255);
                    pk::rect(p.s0, p.zt, p.s1, p.zb, col, IM_COL32(0, 0, 0, 70), 1.0f);
                    if (mp.x >= p.s0 && mp.x <= p.s1 && mp.y >= p.zt && mp.y <= p.zb) hov = &p;
                }
                exportScale(qn[size_t(qsel)].c_str(), lo, hi);
                ImPlot::PopPlotClipRect();
                if (hov && ImPlot::IsPlotHovered()) ImGui::SetTooltip("s %.4g - %.4g %s, depth %.4g - %.4g %s\n%s = %.5g %s", hov->s0, hov->s1, lu, hov->zt, hov->zb, lu,
                                                                      qn[size_t(qsel)].c_str(), valueAt(*hov, step), unit.c_str());
                recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Profiles along the strike (all report steps)")) {
                const std::string xl = std::string("distance along strike (") + lu + ")";
                const std::string yl = qn[size_t(qsel)] + (unit.empty() ? "" : " (" + unit + ")");
                recPanel(("Profiles: " + qn[size_t(qsel)]).c_str(), xl.c_str(), yl.c_str(), "");
                ImPlot::SetupAxes(xl.c_str(), yl.c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                // profile = largest value over the depth of every strike position (the open layer of a contained fracture)
                std::vector<double> sKeys;
                for (const auto& p : pcs) if (p.frac == frac) sKeys.push_back(0.5 * (p.s0 + p.s1));
                std::sort(sKeys.begin(), sKeys.end());
                sKeys.erase(std::unique(sKeys.begin(), sKeys.end(), [](double a, double b) { return std::fabs(a - b) < 1e-6 * (1 + std::fabs(a)); }), sKeys.end());
                const int every = std::max(1, ns / 24);
                std::vector<int> drawn;
                for (int s = 0; s < ns; s += every) drawn.push_back(s);
                if (drawn.back() != ns - 1) drawn.push_back(ns - 1);
                if (std::find(drawn.begin(), drawn.end(), step) == drawn.end()) drawn.push_back(step);   // the selected step always
                for (int s : drawn) {
                    std::vector<double> y(sKeys.size(), -1e300);
                    for (const auto& p : pcs) {
                        if (p.frac != frac) continue;
                        const double m = 0.5 * (p.s0 + p.s1);
                        const size_t i = size_t(std::lower_bound(sKeys.begin(), sKeys.end(), m - 1e-6 * (1 + std::fabs(m))) - sKeys.begin());
                        const double x = valueAt(p, s);
                        if (i < y.size() && std::isfinite(x)) y[i] = std::max(y[i], x);
                    }
                    std::vector<double> xx, yy;
                    for (size_t i = 0; i < y.size(); ++i) if (y[i] > -1e299) { xx.push_back(sKeys[i]); yy.push_back(y[i]); }
                    if (xx.empty()) continue;
                    const ImU32 col = cmU32(ns > 1 ? float(s) / float(ns - 1) : 1.0f);
                    char lab[48]; std::snprintf(lab, sizeof(lab), "t = %.4g d", R.snaps[size_t(s)].time);
                    if (fig::recorder().recording()) recSeries(lab, xx.data(), yy.data(), int(xx.size()), col, s == step ? 2.5f : 1.2f);
                    ImPlotSpec sp; sp.LineColor = ImGui::ColorConvertU32ToFloat4(col); sp.LineWeight = s == step ? 2.8f : 1.2f;
                    if (s != step) sp.Flags = ImPlotItemFlags_NoLegend;
                    ImPlot::PlotLine(lab, xx.data(), yy.data(), int(xx.size()), sp);
                }
                recEnd(); ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    // ======================================================== phase field
    if (hasPf && ui::tabItem("Phase field")) {
        static int rec = -1, what = 0;
        std::vector<int> ids;
        for (int i = 0; i < int(R.hfPf->size()); ++i) if ((*R.hfPf)[size_t(i)].field) ids.push_back(i);
        if (rec < 0 || rec >= int(ids.size())) rec = int(ids.size()) - 1;
        ImGui::SetNextItemWidth(fs * 14); ImGui::SliderInt("Phase-field step", &rec, 0, int(ids.size()) - 1, "%d");
        const Simulator::PfRecord& PR = (*R.hfPf)[size_t(ids[size_t(rec)])];
        const phasefield::Field& F = *PR.field;
        ImGui::SameLine(); ImGui::TextDisabled("t = %.4g d, half-length %.4g m, crack volume %.4g m2/m, pressure %.4g MPa", PR.t, PR.halfLength, PR.volume, PR.pressure * 1e-6);
        static const char* ws[] = {"Damage 1 - v", "sigma_xx change (MPa)", "sigma_yy change (MPa, normal to the crack)", "sigma_xy change (MPa)", "|u| (mm)", "u_y opening displacement (mm)"};
        ImGui::SetNextItemWidth(fs * 16); ImGui::Combo("##pfw", &what, ws, 6);
        ImGui::SameLine(); cmapCombo("##cmpf");
        static bool core = true;
        ImGui::SameLine(); ImGui::Checkbox("Zoom to the crack", &core);
        const int nx1 = int(F.x.size()), ny1 = int(F.y.size()), nxe = nx1 - 1, nye = ny1 - 1;
        auto elemVal = [&](int i, int j) -> double {
            const size_t e = size_t(i) + size_t(nxe) * size_t(j);
            auto nodeAvg = [&](const std::vector<float>& a) {
                const size_t n0 = size_t(i) + size_t(nx1) * size_t(j);
                return 0.25 * (double(a[n0]) + a[n0 + 1] + a[n0 + size_t(nx1)] + a[n0 + size_t(nx1) + 1]);
            };
            switch (what) {
                case 0: return 1.0 - nodeAvg(F.v);
                case 1: return e < F.sxx.size() ? -F.sxx[e] * 1e-6 : NAN;   // compression +
                case 2: return e < F.syy.size() ? -F.syy[e] * 1e-6 : NAN;
                case 3: return e < F.sxy.size() ? -F.sxy[e] * 1e-6 : NAN;
                case 4: {
                    const size_t n0 = size_t(i) + size_t(nx1) * size_t(j);
                    if (F.v[n0] < 0.5f) return NAN;
                    return std::hypot(double(F.ux[n0]), double(F.uy[n0])) * 1e3;
                }
                default: { const size_t n0 = size_t(i) + size_t(nx1) * size_t(j); return F.v[n0] < 0.5f ? NAN : double(F.uy[n0]) * 1e3; }
            }
        };
        const double half = std::max(PR.halfLength, 1e-3), xr = core ? 1.6 * half : 1e300, yr = core ? 0.8 * half : 1e300;
        double lo = 1e300, hi = -1e300;
        for (int j = 0; j < nye; ++j) for (int i = 0; i < nxe; ++i) {
            const double xc = 0.5 * (F.x[size_t(i)] + F.x[size_t(i) + 1]), yc = 0.5 * (F.y[size_t(j)] + F.y[size_t(j) + 1]);
            if (std::fabs(xc) > xr || std::fabs(yc) > yr) continue;
            const double x = elemVal(i, j);
            if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); }
        }
        if (what == 0) { lo = 0; hi = 1; }
        else if (what != 4) { const double m = std::max(std::fabs(lo), std::fabs(hi)); lo = -m; hi = m; }
        if (!(hi > lo)) { lo = 0; hi = 1; }
        colorBarH(ws[what], lo, hi, fs * 14, "");
        const ImVec2 av = ImGui::GetContentRegionAvail();
        fig::recorder().beginFigure("run.phasefield", 1, 2);
        if (ImPlot::BeginSubplots("##pfsub", 1, 2, ImVec2(av.x, av.y), 0, nullptr, nullptr)) {
            if (ImPlot::BeginPlot("Phase-field model")) {
                recPanel(ws[what], "x along the strike (m)", "y normal (m)", "");
                ImPlot::SetupAxes("x along the strike (m)", "y normal (m)", 0, 0);
                if (core) ImPlot::SetupAxesLimits(-xr, xr, -yr, yr, ImPlotCond_Always);
                ImPlot::PushPlotClipRect();
                for (int j = 0; j < nye; ++j) for (int i = 0; i < nxe; ++i) {
                    if (std::fabs(F.x[size_t(i)]) > xr * 1.05 + 1e-9 && std::fabs(F.x[size_t(i) + 1]) > xr * 1.05 + 1e-9) continue;
                    if (std::fabs(F.y[size_t(j)]) > yr * 1.05 + 1e-9 && std::fabs(F.y[size_t(j) + 1]) > yr * 1.05 + 1e-9) continue;
                    const double x = elemVal(i, j);
                    const ImU32 col = std::isfinite(x) ? cmU32(float((x - lo) / (hi - lo))) : IM_COL32(40, 40, 44, 255);
                    pk::rect(F.x[size_t(i)], F.y[size_t(j)], F.x[size_t(i) + 1], F.y[size_t(j) + 1], col);
                }
                exportScale(ws[what], lo, hi);
                ImPlot::PopPlotClipRect();
                recEnd(); ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Crack opening w(x)")) {
                recPanel("Crack opening w(x)", "x along the strike (m)", "opening (mm)", "");
                ImPlot::SetupAxes("x along the strike (m)", "opening (mm)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                for (int r = 0; r <= rec; ++r) {
                    const phasefield::Field& Fr = *(*R.hfPf)[size_t(ids[size_t(r)])].field;
                    if (Fr.cod.empty()) continue;
                    std::vector<double> xx, yy;
                    for (size_t i = 0; i < Fr.cod.size(); ++i) if (std::fabs(Fr.codX[i]) <= xr) { xx.push_back(Fr.codX[i]); yy.push_back(Fr.cod[i] * 1e3); }
                    const ImU32 col = cmU32(rec > 0 ? float(r) / float(rec) : 1.0f);
                    char lab[48]; std::snprintf(lab, sizeof(lab), "t = %.4g d", (*R.hfPf)[size_t(ids[size_t(r)])].t);
                    if (fig::recorder().recording()) recSeries(lab, xx.data(), yy.data(), int(xx.size()), col, r == rec ? 2.5f : 1.0f);
                    ImPlotSpec sp; sp.LineColor = ImGui::ColorConvertU32ToFloat4(col); sp.LineWeight = r == rec ? 2.6f : 1.0f;
                    if (r != rec) sp.Flags = ImPlotItemFlags_NoLegend;
                    ImPlot::PlotLine(lab, xx.data(), yy.data(), int(xx.size()), sp);
                }
                recEnd(); ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    // ======================================================== C31 well section (curtain along the trajectory) + log tracks
    if (!R.input.wells.empty() && ui::tabItem("Well section")) {
        static int wsel = 0, rangeMode = 0, axisMode = 0, nTracks = 3;
        static float marginCells = 4.0f;
        static std::string arr = "SWAT";
        static std::string trk[4] = {"PRESSURE", "SWAT", "SOIL", "SGAS"};
        static bool showInit = true;
        const std::vector<std::string> names = arrayList(R);
        if (std::find(names.begin(), names.end(), arr) == names.end()) arr = "PRESSURE";
        if (wsel < 0 || wsel >= int(R.input.wells.size())) wsel = 0;
        ImGui::SetNextItemWidth(fs * 8);
        if (ImGui::BeginCombo("Well##ws", R.input.wells[size_t(wsel)].name.c_str(), ImGuiComboFlags_HeightLarge)) {
            for (int w = 0; w < int(R.input.wells.size()); ++w) if (ImGui::Selectable(R.input.wells[size_t(w)].name.c_str(), w == wsel)) wsel = w;
            ImGui::EndCombo();
        }
        auto arrCombo = [&](const char* id, std::string& a, float wfs) {
            ImGui::SetNextItemWidth(fs * wfs);
            if (ImGui::BeginCombo(id, a.c_str(), ImGuiComboFlags_HeightLarge)) {
                for (const auto& n : names) {
                    const std::string u = unitOfArray(n, us);
                    if (ImGui::Selectable((n + (u.empty() ? "" : "  (" + u + ")")).c_str(), n == a)) a = n;
                }
                ImGui::EndCombo();
            }
        };
        ImGui::SameLine(); arrCombo("Array##ws", arr, 9);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 12); ImGui::SliderInt("Step##ws", &step, 0, ns - 1, "%d");
        ImGui::SameLine(); ImGui::TextDisabled("t = %.4g d", R.snaps[size_t(step)].time);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 7); ImGui::SliderFloat("Margin (cells)##ws", &marginCells, 0.0f, 20.0f, "%.0f");
        cmapCombo("##cmws");
        ImGui::SameLine();
        static const char* rms[] = {"Range: all steps", "Range: this step", "Range: symmetric around 0"};
        ImGui::SetNextItemWidth(fs * 13); ImGui::Combo("##wsrm", &rangeMode, rms, 3);
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 5); ImGui::SliderInt("Tracks##ws", &nTracks, 0, 4);
        for (int t = 0; t < nTracks; ++t) { ImGui::SameLine(); char id[24]; std::snprintf(id, sizeof id, "##trk%d", t); arrCombo(id, trk[t], 7); }
        ImGui::SameLine(); ImGui::SetNextItemWidth(fs * 6);
        static const char* axs[] = {"TVD", "MD"};
        ImGui::Combo("Track axis##ws", &axisMode, axs, 2);
        ImGui::SameLine(); ImGui::Checkbox("Initial##ws", &showInit);
        ui::hint("Curtain: the cells under the well's plan trajectory (extended by the margin at both ends; a vertical well gets an "
                 "east-west curtain) coloured by the array; the trajectory and its connections (open / shut at this step) on top. "
                 "Tracks: the arrays in the connection cells against true vertical depth or measured depth (dashed: first report step).");
        const std::vector<WellGlyph> glyphs = wellGlyphs(g, R.input, step);
        const WellGlyph* G = nullptr;
        for (const auto& gl : glyphs) if (gl.name == R.input.wells[size_t(wsel)].name) G = &gl;
        if (!G || G->path.size() < 2) { ui::hint("This well has no connection in the grid."); ImGui::EndTabItem(); ImGui::EndTabBar(); return; }
        // ---- curtain polyline in plan (from the trajectory), extended at both ends
        const double cdx = (g.bbMax.x - g.bbMin.x) / std::max(1, g.nx), cdy = (g.bbMax.y - g.bbMin.y) / std::max(1, g.ny);
        const double cell = std::max(1e-9, std::min(cdx, cdy));
        std::vector<std::array<double, 2>> pl;
        std::vector<double> pathS(G->path.size(), 0.0);
        for (size_t k = 0; k < G->path.size(); ++k) {
            const std::array<double, 2> q{G->path[k][0], G->path[k][1]};
            if (pl.empty() || std::hypot(q[0] - pl.back()[0], q[1] - pl.back()[1]) > 1e-3 * cell) pl.push_back(q);
        }
        std::array<double, 2> d0{1, 0}, d1{1, 0};
        if (pl.size() >= 2) {
            const double l0 = std::hypot(pl[1][0] - pl[0][0], pl[1][1] - pl[0][1]), l1 = std::hypot(pl.back()[0] - pl[pl.size() - 2][0], pl.back()[1] - pl[pl.size() - 2][1]);
            d0 = {(pl[1][0] - pl[0][0]) / l0, (pl[1][1] - pl[0][1]) / l0};
            d1 = {(pl.back()[0] - pl[pl.size() - 2][0]) / l1, (pl.back()[1] - pl[pl.size() - 2][1]) / l1};
        }
        const double margin = std::max(double(marginCells), 0.0) * cell + (pl.size() < 2 ? 0.5 * cell : 0.0);
        pl.insert(pl.begin(), {pl.front()[0] - d0[0] * margin, pl.front()[1] - d0[1] * margin});
        pl.push_back({pl.back()[0] + d1[0] * margin, pl.back()[1] + d1[1] * margin});
        std::vector<double> cum(pl.size(), 0.0);
        for (size_t k = 1; k < pl.size(); ++k) cum[k] = cum[k - 1] + std::hypot(pl[k][0] - pl[k - 1][0], pl[k][1] - pl[k - 1][1]);
        const double Ltot = cum.back();
        auto sOf = [&](double x, double y) {   // curtain coordinate of a plan point (nearest point of the polyline)
            double best = 1e300, sb = 0;
            for (size_t k = 1; k < pl.size(); ++k) {
                const double ax = pl[k - 1][0], ay = pl[k - 1][1], bx = pl[k][0] - ax, by = pl[k][1] - ay, l2 = bx * bx + by * by;
                const double t = l2 > 0 ? std::clamp(((x - ax) * bx + (y - ay) * by) / l2, 0.0, 1.0) : 0.0;
                const double d = std::hypot(x - ax - t * bx, y - ay - t * by);
                if (d < best) { best = d; sb = cum[k - 1] + t * std::sqrt(l2); }
            }
            return sb;
        };
        for (size_t k = 0; k < G->path.size(); ++k) pathS[k] = sOf(G->path[k][0], G->path[k][1]);
        // ---- columns along the curtain (cached: nearest column centre of every sample)
        struct Seg { double s0, s1; int i, j; };
        static std::vector<Seg> segs; static std::string segKey;
        {
            char kb[160]; std::snprintf(kb, sizeof kb, "%s|%.3f|%p|%d|%d|%d|%zu", R.input.wells[size_t(wsel)].name.c_str(), marginCells, static_cast<const void*>(&g), g.nx, g.ny, g.nz, G->path.size());
            if (segKey != kb) {
                segKey = kb; segs.clear();
                std::vector<double> cx(size_t(g.nx) * g.ny, NAN), cy(size_t(g.nx) * g.ny, NAN);
                for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
                    for (int k = 0; k < g.nz; ++k) {   // first layer with corner data
                        const size_t q = size_t(i + g.nx * (j + g.ny * k));
                        if (q >= g.corners.size()) break;
                        const auto& C = g.corners[q];
                        cx[size_t(i + g.nx * j)] = 0.25 * (C[0].x + C[1].x + C[2].x + C[3].x);
                        cy[size_t(i + g.nx * j)] = 0.25 * (C[0].y + C[1].y + C[2].y + C[3].y);
                        if (q < g.activeOf.size() && g.activeOf[q] >= 0) break;
                    }
                }
                const double ds = cell / 4;
                int ci = -1, cj = -1; double s0 = 0;
                for (double sv = 0; sv <= Ltot + 1e-9; sv += ds) {
                    size_t k = 1; while (k + 1 < pl.size() && cum[k] < sv) ++k;
                    const double t = cum[k] > cum[k - 1] ? (sv - cum[k - 1]) / (cum[k] - cum[k - 1]) : 0.0;
                    const double x = pl[k - 1][0] + t * (pl[k][0] - pl[k - 1][0]), y = pl[k - 1][1] + t * (pl[k][1] - pl[k - 1][1]);
                    int bi = -1, bj = -1; double bd = 1e300;
                    for (int jj = 0; jj < g.ny; ++jj) for (int ii = 0; ii < g.nx; ++ii) {
                        const size_t c = size_t(ii + g.nx * jj);
                        if (!std::isfinite(cx[c])) continue;
                        const double d = (cx[c] - x) * (cx[c] - x) + (cy[c] - y) * (cy[c] - y);
                        if (d < bd) { bd = d; bi = ii; bj = jj; }
                    }
                    if (std::sqrt(bd) > 1.5 * std::max(cdx, cdy)) { bi = -1; bj = -1; }
                    if (bi != ci || bj != cj) {
                        if (ci >= 0) segs.push_back({s0, sv, ci, cj});
                        ci = bi; cj = bj; s0 = sv;
                    }
                }
                if (ci >= 0) segs.push_back({s0, Ltot, ci, cj});
            }
        }
        // ---- colour range
        const std::vector<float>& v = R.snaps[size_t(step)].field(arr);
        double lo = 1e300, hi = -1e300;
        auto accCurtain = [&](const std::vector<float>& a) {
            for (const Seg& sg : segs) for (int k = 0; k < g.nz; ++k) {
                const size_t q = size_t(sg.i + g.nx * (sg.j + g.ny * k));
                const int ac = q < g.activeOf.size() ? g.activeOf[q] : -1;
                if (ac >= 0 && size_t(ac) < a.size() && std::isfinite(a[size_t(ac)])) { lo = std::min(lo, double(a[size_t(ac)])); hi = std::max(hi, double(a[size_t(ac)])); }
            }
        };
        if (rangeMode == 1) accCurtain(v);
        else {
            static std::string rKey; static double rLo = 0, rHi = 0;
            char kb[200]; std::snprintf(kb, sizeof kb, "%s|%s|%d|%.9g|%p|%u", arr.c_str(), segKey.c_str(), ns, R.snaps.back().time, static_cast<const void*>(&R.snaps[0]), unsigned(runner_.version()));
            if (rKey != kb) { for (int s2 = 0; s2 < ns; ++s2) accCurtain(R.snaps[size_t(s2)].field(arr)); rKey = kb; rLo = lo; rHi = hi; }
            else { lo = rLo; hi = rHi; }
        }
        if (rangeMode == 2 || (signedArray(arr) && rangeMode != 1)) { const double m = std::max(std::fabs(lo), std::fabs(hi)); lo = -m; hi = m; }
        if (!(hi >= lo)) { ui::hint("No values of this array along the well."); ImGui::EndTabItem(); ImGui::EndTabBar(); return; }
        if (hi - lo < 1e-12) hi = lo + 1;
        const std::string unit = unitOfArray(arr, us);
        colorBarH(arr.c_str(), lo, hi, fs * 14, unit.c_str());
        // ---- connection cells (completion order) for the tracks
        struct Conn { int q = -1; double tvd = 0, md = 0; bool open = true; };
        std::vector<Conn> conns;
        {
            const WellInput& W = R.input.wells[size_t(wsel)];
            for (const auto& cp : W.completions) {
                if (cp.lgr >= 0) continue;
                const bool on = cp.open && cp.fromStep <= step && (cp.untilStep < 0 || step < cp.untilStep);
                for (int k = cp.k1; k <= cp.k2; ++k) {
                    if (cp.i < 0 || cp.i >= g.nx || cp.j < 0 || cp.j >= g.ny || k < 0 || k >= g.nz) continue;
                    Conn c; c.q = cp.i + g.nx * (cp.j + g.ny * k); c.open = on;
                    if (size_t(c.q) < g.center.size()) c.tvd = g.center[size_t(c.q)].z;
                    conns.push_back(c);
                }
            }
            double md = 0; std::array<double, 3> prev{G->path[0][0], G->path[0][1], G->path[0][2]};
            for (auto& c : conns) {
                if (size_t(c.q) >= g.center.size()) continue;
                const auto& P = g.center[size_t(c.q)];
                md += std::sqrt((P.x - prev[0]) * (P.x - prev[0]) + (P.y - prev[1]) * (P.y - prev[1]) + (P.z - prev[2]) * (P.z - prev[2]));
                prev = {double(P.x), double(P.y), double(P.z)};
                c.md = md;
            }
        }
        const int nt = conns.empty() ? 0 : nTracks;
        const ImVec2 av = ImGui::GetContentRegionAvail();
        fig::recorder().beginFigure("run.wellsection", 1, 1 + nt);
        std::vector<float> ratios(size_t(1 + nt), 1.0f); ratios[0] = 3.5f;
        if (ImPlot::BeginSubplots("##wssub", 1, 1 + nt, ImVec2(av.x, av.y), axisMode == 0 ? ImPlotSubplotFlags_LinkRows : 0, nullptr, ratios.data())) {
            double zMin = 1e300, zMax = -1e300;
            for (const Seg& sg : segs) for (int k = 0; k < g.nz; ++k) {
                const size_t q = size_t(sg.i + g.nx * (sg.j + g.ny * k));
                if (q >= g.corners.size() || q >= g.activeOf.size() || g.activeOf[q] < 0) continue;
                for (int c = 0; c < 8; ++c) { zMin = std::min(zMin, double(g.corners[q][c].z)); zMax = std::max(zMax, double(g.corners[q][c].z)); }
            }
            if (!(zMax > zMin)) { zMin = g.bbMin.z; zMax = g.bbMax.z; }
            const double zTop = std::min(zMin, double(G->path[0][2]));
            if (ImPlot::BeginPlot(("Curtain along " + R.input.wells[size_t(wsel)].name).c_str())) {
                const std::string xl = std::string("distance along the trajectory (") + lu + ")", yl = std::string("depth (") + lu + ")";
                recPanel((arr + " along " + R.input.wells[size_t(wsel)].name).c_str(), xl.c_str(), yl.c_str(), "");
                ImPlot::SetupAxes(xl.c_str(), yl.c_str(), 0, ImPlotAxisFlags_Invert);
                ImPlot::SetupAxesLimits(0, Ltot, zTop - 0.02 * (zMax - zTop), zMax + 0.02 * (zMax - zTop), ImPlotCond_Once);
                const double bx[2] = {0, Ltot}, by[2] = {zTop, zMax};
                ImPlotSpec sp; sp.Flags = ImPlotItemFlags_NoLegend; sp.MarkerSize = 0; sp.LineColor = ImVec4(0, 0, 0, 0);
                ImPlot::PlotScatter("##wsbb", bx, by, 2, sp);
                ImPlot::PushPlotClipRect();
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                int hovered = -1;
                for (const Seg& sg : segs)
                    for (int k = 0; k < g.nz; ++k) {
                        const size_t q = size_t(sg.i + g.nx * (sg.j + g.ny * k));
                        if (q >= g.corners.size()) continue;
                        const int a = q < g.activeOf.size() ? g.activeOf[q] : -1;
                        if (a < 0) continue;
                        const auto& C = g.corners[q];
                        const double zt = 0.25 * (C[0].z + C[1].z + C[2].z + C[3].z), zb = 0.25 * (C[4].z + C[5].z + C[6].z + C[7].z);
                        ImU32 col = IM_COL32(70, 70, 76, 255);
                        if (size_t(a) < v.size() && std::isfinite(v[size_t(a)])) col = cmU32(float((v[size_t(a)] - lo) / (hi - lo)));
                        const double qx[4] = {sg.s0, sg.s1, sg.s1, sg.s0}, qy[4] = {zt, zt, zb, zb};
                        pk::quad(qx, qy, col);
                        if (mp.x >= sg.s0 && mp.x <= sg.s1 && mp.y >= std::min(zt, zb) && mp.y <= std::max(zt, zb)) hovered = int(q);
                    }
                exportScale(arr.c_str(), lo, hi);
                // trajectory and connections
                const ImU32 tc = G->injector ? IM_COL32(70, 140, 255, 255) : IM_COL32(40, 40, 40, 255);
                for (size_t k = 1; k < G->path.size(); ++k) pk::line(pathS[k - 1], G->path[k - 1][2], pathS[k], G->path[k][2], tc, 2.5f, k == 1 ? "trajectory" : nullptr);
                for (const auto& pf : G->perfs) {
                    const ImU32 pc = pf.open ? (G->injector ? IM_COL32(90, 170, 255, 255) : IM_COL32(60, 220, 110, 255)) : IM_COL32(150, 150, 150, 255);
                    pk::line(sOf(pf.a[0], pf.a[1]), pf.a[2], sOf(pf.b[0], pf.b[1]), pf.b[2], pc, 6.0f, pf.open ? "open connection" : "shut connection");
                }
                pk::text(pathS[0], G->path[0][2], R.input.wells[size_t(wsel)].name.c_str(), tc, ImVec2(4, -2), 0);
                ImPlot::PopPlotClipRect();
                if (ImPlot::IsPlotHovered() && hovered >= 0) {
                    const int a = g.activeOf[size_t(hovered)];
                    const int I = hovered % g.nx, J = (hovered / g.nx) % g.ny, K = hovered / (g.nx * g.ny);
                    if (a >= 0 && size_t(a) < v.size()) ImGui::SetTooltip("cell (%d, %d, %d)\n%s = %.5g %s", I + 1, J + 1, K + 1, arr.c_str(), v[size_t(a)], unit.c_str());
                }
                recEnd();
                ImPlot::EndPlot();
            }
            for (int t = 0; t < nt; ++t) {
                const std::string& tn = trk[t];
                const std::string tu = unitOfArray(tn, us);
                if (ImPlot::BeginPlot(("##track" + std::to_string(t)).c_str())) {
                    const std::string xl = tn + (tu.empty() ? "" : " (" + tu + ")");
                    const std::string yl = std::string(axisMode == 0 ? "TVD (" : "MD (") + lu + ")";
                    recPanel(("Track " + tn).c_str(), xl.c_str(), yl.c_str(), "");
                    ImPlot::SetupAxes(xl.c_str(), yl.c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_Invert | (axisMode == 1 ? ImPlotAxisFlags_AutoFit : 0));
                    auto series = [&](int si, const char* lab, ImU32 col, float wgt, bool dash) {
                        const std::vector<float>& a = R.snaps[size_t(si)].field(tn);
                        std::vector<double> xx, yy;
                        for (const auto& c : conns) {
                            const int ac = size_t(c.q) < g.activeOf.size() ? g.activeOf[size_t(c.q)] : -1;
                            if (ac < 0 || size_t(ac) >= a.size() || !std::isfinite(a[size_t(ac)])) continue;
                            xx.push_back(a[size_t(ac)]); yy.push_back(axisMode == 0 ? c.tvd : c.md);
                        }
                        if (xx.empty()) return;
                        if (fig::recorder().recording()) recSeries(lab, xx.data(), yy.data(), int(xx.size()), col, wgt, false, true);
                        ImPlotSpec sp; sp.LineColor = ImGui::ColorConvertU32ToFloat4(col); sp.LineWeight = wgt;
                        sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = dash ? 0 : 3.5f; sp.MarkerFillColor = sp.LineColor;
                        if (dash) sp.Flags = ImPlotItemFlags_NoLegend;
                        ImPlot::PlotLine(lab, xx.data(), yy.data(), int(xx.size()), sp);
                    };
                    if (showInit && step > 0) series(0, "first step", IM_COL32(150, 150, 150, 255), 1.2f, true);
                    series(step, "this step", IM_COL32(220, 90, 40, 255), 2.2f, false);
                    for (const auto& c : conns)   // open / shut flags at the right edge
                        if (!c.open) pk::marker(ImPlot::GetPlotLimits().X.Max, axisMode == 0 ? c.tvd : c.md, fs * 0.25f, IM_COL32(150, 150, 150, 255), 0, 1, 1.0f, nullptr);
                    recEnd();
                    ImPlot::EndPlot();
                }
            }
            ImPlot::EndSubplots();
        }
        fig::recorder().endFigure();
        ImGui::EndTabItem();
    }
    // ======================================================== geomechanics post-processing
    if (ui::tabItem("Geomechanics (post)")) {
        bool have = false;
        for (const auto& e : R.snaps.back().extra) if (e.first == "UPLIFT") have = true;
        if (have) ui::hint("This run has geomechanics maps (GEOMECH in the deck or an earlier post-processing). Computing again "
                           "replaces them with the moduli below.");
        postGeoControls(0, R.grid, fld, &R.snaps);
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}
