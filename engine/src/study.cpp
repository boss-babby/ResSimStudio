// Multi-run studies: batch runner, sensitivity analysis and history matching.
#include "ressim/study.h"
#include "ressim/geostat.h"
#include "ressim/ecl.h"
#include "ressim/coupled_bo.h"
#include "ressim/coupled_deck.h"
#include <thread>
#include <chrono>
#include <memory>
#include <mutex>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <cstdio>
#include <cstring>

namespace ressim {

static const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ======================================================================== parameters
namespace {
// apply f to the SWOF (gas = false) or SGOF (gas = true) table of every saturation region
template <class F> void forEachTable(Case& c, bool gas, F f) {
    f(gas ? c.props.sgof : c.props.swof);
    for (auto& R : c.props.satExtra) f(gas ? R.sgof : R.swof);
}
}
const char* toString(ParamKind k) {
    switch (k) {
        case ParamKind::PermXY: return "Perm X/Y multiplier";
        case ParamKind::PermZ: return "Perm Z multiplier";
        case ParamKind::Poro: return "Porosity multiplier";
        case ParamKind::FaultMult: return "Fault transmissibility (MULTFLT)";
        case ParamKind::AquiferStrength: return "Aquifer strength multiplier";
        case ParamKind::KrwMax: return "krw end-point multiplier";
        case ParamKind::KrwExp: return "krw curvature exponent";
        case ParamKind::KroExp: return "krow curvature exponent";
        case ParamKind::WaterVisc: return "Water viscosity multiplier";
        case ParamKind::RockComp: return "Rock compressibility multiplier";
        case ParamKind::OWC: return "OWC depth";
        case ParamKind::DatumPressure: return "Initial (datum) pressure";
        case ParamKind::WellSkin: return "Well skin";
        case ParamKind::WellPI: return "Well PI multiplier (WPIMULT)";
        case ParamKind::OilVisc: return "Oil viscosity multiplier";
        case ParamKind::GOC: return "GOC depth";
        case ParamKind::Ntg: return "Net-to-gross multiplier";
        case ParamKind::MultZ: return "Vertical transmissibility (MULTZ)";
        case ParamKind::KrgExp: return "krg curvature exponent";
        case ParamKind::KrogExp: return "krog curvature exponent";
        case ParamKind::KrowMax: return "krow end-point multiplier";
        case ParamKind::PcMult: return "Capillary pressure multiplier";
        case ParamKind::InjRate: return "Injection rate multiplier";
        case ParamKind::ProdRate: return "Production rate multiplier";
        case ParamKind::GeoRange: return "Geomodel variogram range multiplier";
        case ParamKind::GeoSeed: return "Geomodel realisation (seed)";
    }
    return "?";
}
bool paramIsMultiplier(ParamKind k) {
    return !(k == ParamKind::OWC || k == ParamKind::DatumPressure || k == ParamKind::WellSkin || k == ParamKind::FaultMult ||
             k == ParamKind::KrwExp || k == ParamKind::KroExp || k == ParamKind::GOC || k == ParamKind::KrgExp || k == ParamKind::KrogExp ||
             k == ParamKind::GeoSeed);
}
bool paramUsesBox(ParamKind k) {
    return k == ParamKind::PermXY || k == ParamKind::PermZ || k == ParamKind::Poro || k == ParamKind::Ntg || k == ParamKind::MultZ;
}

std::string paramLabel(const StudyParam& p, const Case& c) {
    if (!p.label.empty()) return p.label;
    std::ostringstream o;
    switch (p.kind) {
        case ParamKind::PermXY: o << "KXY"; break;
        case ParamKind::PermZ: o << "KZ"; break;
        case ParamKind::Poro: o << "PORO"; break;
        case ParamKind::FaultMult: o << "MULTFLT " << (p.target < int(c.grid.faults.size()) ? c.grid.faults[p.target].name : "?"); break;
        case ParamKind::AquiferStrength: o << "AQ" << (p.target < int(c.aquifers.size()) ? c.aquifers[p.target].id : 0) << " size"; break;
        case ParamKind::KrwMax: o << "krw max"; break;
        case ParamKind::KrwExp: o << "krw exp"; break;
        case ParamKind::KroExp: o << "krow exp"; break;
        case ParamKind::WaterVisc: o << "mu_w"; break;
        case ParamKind::RockComp: o << "c_rock"; break;
        case ParamKind::OWC: o << "OWC"; break;
        case ParamKind::DatumPressure: o << "P_init"; break;
        case ParamKind::WellSkin: o << "skin " << (p.target < int(c.wells.size()) ? c.wells[p.target].name : "?"); break;
        case ParamKind::WellPI: o << "PI " << (p.target < int(c.wells.size()) ? c.wells[p.target].name : "?"); break;
        case ParamKind::OilVisc: o << "mu_o"; break;
        case ParamKind::GOC: o << "GOC"; break;
        case ParamKind::Ntg: o << "NTG"; break;
        case ParamKind::MultZ: o << "MULTZ"; break;
        case ParamKind::KrgExp: o << "krg exp"; break;
        case ParamKind::KrogExp: o << "krog exp"; break;
        case ParamKind::KrowMax: o << "krow max"; break;
        case ParamKind::PcMult: o << "Pcow"; break;
        case ParamKind::InjRate: o << "q_inj " << (p.target >= 0 && p.target < int(c.wells.size()) ? c.wells[p.target].name : "all"); break;
        case ParamKind::ProdRate: o << "q_prod " << (p.target >= 0 && p.target < int(c.wells.size()) ? c.wells[p.target].name : "all"); break;
        case ParamKind::GeoRange: case ParamKind::GeoSeed: {
            static const char* kP[] = {"PORO", "PERMX", "NTG", "FACIES"};
            o << (p.kind == ParamKind::GeoRange ? "range " : "seed ");
            if (p.target >= 0 && p.target < int(c.propModels.size())) o << kP[std::clamp(c.propModels[size_t(p.target)].property, 0, 3)] << " #" << p.target + 1;
            else o << "all models";
            break;
        }
    }
    if (paramUsesBox(p.kind)) {
        const bool full = p.i1 < 0 && p.j1 < 0 && p.i2 < 0 && p.j2 < 0;
        if (full && p.k1 >= 0 && p.k1 == p.k2) o << " L" << p.k1 + 1;
        else if (full && p.k1 >= 0) o << " L" << p.k1 + 1 << "-" << (p.k2 < 0 ? c.grid.nz : p.k2 + 1);
        else if (!full) o << " box";
    }
    return o.str();
}

double paramToUnit(const StudyParam& p, double v) {
    if (p.high == p.low) return 0.5;
    double u;
    if (p.log && p.low > 0 && p.high > 0 && v > 0) u = (std::log(v) - std::log(p.low)) / (std::log(p.high) - std::log(p.low));
    else u = (v - p.low) / (p.high - p.low);
    return std::max(0.0, std::min(1.0, u));
}
double paramFromUnit(const StudyParam& p, double u) {
    u = std::max(0.0, std::min(1.0, u));
    if (p.log && p.low > 0 && p.high > 0) return std::exp(std::log(p.low) + u * (std::log(p.high) - std::log(p.low)));
    return p.low + u * (p.high - p.low);
}

Case applyParams(const Case& base, const std::vector<StudyParam>& params, const std::vector<double>& values) {
    Case c = base;
    {   // C40 geomodel parameters first: re-simulate the property models, then the multipliers act on the new fields
        std::vector<PropModelInput> pm = c.propModels;
        std::vector<char> redo(pm.size(), 0);
        for (size_t i = 0; i < params.size(); ++i) {
            const StudyParam& p = params[i];
            if (p.kind != ParamKind::GeoRange && p.kind != ParamKind::GeoSeed) continue;
            const double v = (p.active && i < values.size()) ? values[i] : p.base;
            for (size_t m = 0; m < pm.size(); ++m) {
                if (p.target >= 0 && size_t(p.target) != m) continue;
                if (p.kind == ParamKind::GeoRange) {
                    if (v == 1.0) continue;
                    pm[m].rangeMax *= v; pm[m].rangeMin *= v; pm[m].rangeVert *= v;
                } else {
                    const unsigned sd = unsigned(std::max(1L, std::lround(v)));
                    if (sd == pm[m].seed) continue;
                    pm[m].seed = sd;
                }
                redo[m] = 1;
            }
        }
        for (size_t m = 0; m < pm.size(); ++m) if (redo[m]) geostat::applyPropModel(c, pm[m]);
    }
    GridInput& g = c.grid;
    bool anyPerm = false;
    for (auto& p : params) if (p.kind == ParamKind::PermXY || p.kind == ParamKind::PermZ) anyPerm = true;
    if (anyPerm && g.permz.empty()) g.permz = g.permx;            // keep kz independent of kx multipliers
    for (size_t i = 0; i < params.size(); ++i) {
        const StudyParam& p = params[i];
        const double v = (p.active && i < values.size()) ? values[i] : p.base;
        if (p.kind == ParamKind::GeoRange || p.kind == ParamKind::GeoSeed) continue;   // done above
        if (paramIsMultiplier(p.kind) && v == 1.0) continue;
        switch (p.kind) {
            case ParamKind::GeoRange: case ParamKind::GeoSeed: break;   // C40: applied before the loop
            case ParamKind::PermXY: case ParamKind::PermZ: case ParamKind::Poro: case ParamKind::Ntg: case ParamKind::MultZ: {
                const size_t ncell = size_t(g.nx) * g.ny * g.nz;
                if (p.kind == ParamKind::Ntg && g.ntg.size() != ncell) g.ntg.assign(ncell, 1.0);
                if (p.kind == ParamKind::MultZ && g.multz.size() != ncell) g.multz.assign(ncell, 1.0);
                const int i1 = p.i1 < 0 ? 0 : p.i1, i2 = p.i2 < 0 ? g.nx - 1 : std::min(p.i2, g.nx - 1);
                const int j1 = p.j1 < 0 ? 0 : p.j1, j2 = p.j2 < 0 ? g.ny - 1 : std::min(p.j2, g.ny - 1);
                const int k1 = p.k1 < 0 ? 0 : p.k1, k2 = p.k2 < 0 ? g.nz - 1 : std::min(p.k2, g.nz - 1);
                for (int k = k1; k <= k2; ++k)
                    for (int j = j1; j <= j2; ++j)
                        for (int ii = i1; ii <= i2; ++ii) {
                            const int q = g.index(ii, j, k);
                            if (p.kind == ParamKind::PermXY) { g.permx[q] *= v; if (!g.permy.empty()) g.permy[q] *= v; }
                            else if (p.kind == ParamKind::PermZ) g.permz[q] *= v;
                            else if (p.kind == ParamKind::Ntg) g.ntg[q] = std::min(1.0, g.ntg[q] * v);
                            else if (p.kind == ParamKind::MultZ) g.multz[q] *= v;
                            else g.poro[q] = std::min(1.0, g.poro[q] * v);
                        }
            } break;
            case ParamKind::FaultMult: if (p.target >= 0 && p.target < int(g.faults.size())) g.faults[p.target].mult = v; break;
            case ParamKind::AquiferStrength:
                if (p.target >= 0 && p.target < int(c.aquifers.size())) {
                    AquiferInput& a = c.aquifers[p.target];
                    if (a.type == AquiferType::Fetkovich) { a.volume *= v; a.pi *= v; } else { a.perm *= v; a.thickness *= v; }
                }
                break;
            case ParamKind::KrwMax: forEachTable(c, false, [&](std::vector<SatRow>& t) { for (auto& r : t) r.krw = std::min(1.0, r.krw * v); }); break;
            case ParamKind::KrowMax: forEachTable(c, false, [&](std::vector<SatRow>& t) { for (auto& r : t) r.kro = std::min(1.0, r.kro * v); }); break;
            case ParamKind::PcMult: forEachTable(c, false, [&](std::vector<SatRow>& t) { for (auto& r : t) r.pc *= v; }); break;
            case ParamKind::KrwExp: case ParamKind::KroExp: case ParamKind::KrgExp: case ParamKind::KrogExp: {
                const bool gas = p.kind == ParamKind::KrgExp || p.kind == ParamKind::KrogExp;
                const bool first = p.kind == ParamKind::KrwExp || p.kind == ParamKind::KrgExp;   // column 2 (krw / krg)
                forEachTable(c, gas, [&](std::vector<SatRow>& t) {
                    double mx = 0;
                    for (auto& r : t) mx = std::max(mx, first ? r.krw : r.kro);
                    if (mx > 0)
                        for (auto& r : t) { double& k = first ? r.krw : r.kro; k = mx * std::pow(std::max(0.0, k / mx), v); }
                });
            } break;
            case ParamKind::OilVisc: {
                auto scale = [&](std::vector<PVDRow>& pvdo, std::vector<PVTORecord>& pvto) {
                    for (auto& r : pvdo) r.mu *= v;
                    for (auto& r : pvto) { r.mu *= v; for (auto& u : r.under) u.mu *= v; }
                };
                scale(c.props.pvdo, c.props.pvto);
                for (auto& R : c.props.pvtExtra) scale(R.pvdo, R.pvto);
            } break;
            case ParamKind::GOC: c.equil.goc = v; break;
            case ParamKind::InjRate: case ParamKind::ProdRate: {
                const bool inj = p.kind == ParamKind::InjRate;
                const std::string wn = p.target >= 0 && p.target < int(c.wells.size()) ? c.wells[p.target].name : std::string();
                auto isInj = [&](const std::string& n) { for (auto& w : c.wells) if (w.name == n) return w.type != WellType::Producer; return false; };
                for (auto& st : c.schedule)
                    for (auto& wc : st.controls) {
                        if (wc.history) continue;
                        if (!wn.empty() ? wc.well != wn : isInj(wc.well) != inj) continue;
                        if (inj) { if (wc.rate > 0) wc.rate *= v; }
                        else for (double* q : {&wc.orat, &wc.wrat, &wc.grat, &wc.lrat, &wc.resv}) if (*q > 0) *q *= v;
                    }
            } break;
            case ParamKind::WaterVisc: c.props.pvtw.visc *= v; for (auto& R : c.props.pvtExtra) R.pvtw.visc *= v; break;
            case ParamKind::RockComp: c.props.rock.cr *= v; break;
            case ParamKind::OWC: c.equil.owc = v; break;
            case ParamKind::DatumPressure: c.equil.datumPressure = v; break;
            case ParamKind::WellSkin:
                if (p.target >= 0 && p.target < int(c.wells.size())) for (auto& cp : c.wells[p.target].completions) cp.skin = v;
                break;
            case ParamKind::WellPI: if (p.target >= 0 && p.target < int(c.wells.size())) c.wells[p.target].piMult *= v; break;
        }
    }
    return c;
}

// ======================================================================== runs
// C28: compositional / CO2STORE and thermal cases in studies (field summary + coupled KPIs; well vectors are not mapped)
static RunResult runSingleDeckEngine(const Case& c, std::atomic<bool>* cancel) {
    RunResult r;
    std::string err;
    coupledbo::Output co;
    auto keepGoing = [cancel]() { return !(cancel && cancel->load()); };
    bool ok = false;
    std::vector<SummaryPoint>& pts = r.summary;
    if (c.thermal.enabled) {
        thermal::DeckModel m;
        if (!thermal::deckModel(c, m, &err)) { r.error = "thermal model: " + err; return r; }
        r.wellNames = m.wellNames;
        std::vector<thermal::Summary> o; thermal::Stats st;
        ok = coupleddeck::runThermal(c, m, o, co, {}, [&](const thermal::Summary&) { return keepGoing(); }, &err, -1, &st);
        for (const auto& S : o) {
            SummaryPoint p; p.time = S.time; p.dt = pts.empty() ? S.time : S.time - pts.back().time;
            p.fopr = S.fopr; p.fwpr = S.fwpr; p.fgpr = S.fgpr; p.fwir = S.fwir; p.fgir = S.fgir;
            p.fopt = S.fopt; p.fwpt = S.fwpt; p.fgpt = S.fgpt; p.fwit = S.fwit; p.fgit = S.fgit; p.fpr = S.fpr; p.fprp = S.fpr; p.oip = S.foip;
            p.fwct = S.fopr + S.fwpr > 0 ? S.fwpr / (S.fopr + S.fwpr) : 0.0; p.fgor = S.fopr > 0 ? S.fgpr / S.fopr : 0.0;
            for (const auto& W : S.wells) {   // C49: per-well rates and BHP (were dropped: WBHP = 0 in regressions / studies)
                WellRates w; w.oil = W.oil; w.water = W.water; w.gas = W.gas; w.waterInj = W.winj; w.bhp = W.bhp; w.open = W.open;
                p.wells.push_back(w);
            }
            pts.push_back(p);
        }
        r.stats.timesteps = st.steps; r.stats.newton = st.newton; r.stats.linear = st.linear;
    } else {
        comp::DeckModel m;
        if (!comp::deckModel(c, m, &err)) { r.error = "compositional model: " + err; return r; }
        r.wellNames = m.wellNames;
        std::vector<comp::CompSummary> o; comp::CompStats st;
        ok = coupleddeck::runComp(c, m, o, co, {}, [&](const comp::CompSummary&) { return keepGoing(); }, &err, -1, &st);
        for (const auto& S : o) {
            SummaryPoint p; p.time = S.time; p.dt = pts.empty() ? S.time : S.time - pts.back().time;
            p.fopr = S.fopr; p.fwpr = S.fwpr; p.fgpr = S.fgpr; p.fwir = S.fwir; p.fgir = S.fgir;
            p.fopt = S.fopt; p.fwpt = S.fwpt; p.fgpt = S.fgpt; p.fwit = S.fwit; p.fgit = S.fgit; p.fpr = S.fpr; p.fprp = S.fpr;
            p.oip = S.foip; p.gip = S.fgip; p.wip = S.fwip;
            p.fwct = S.fopr + S.fwpr > 0 ? S.fwpr / (S.fopr + S.fwpr) : 0.0;   // C49 (were 0)
            p.fgor = S.fopr > 0 ? S.fgpr / S.fopr : 0.0;
            for (const auto& W : S.wells) {   // C49: per-well rates and BHP
                WellRates w; w.oil = W.oil; w.water = W.water; w.gas = W.gas; w.waterInj = W.winj; w.gasInj = W.ginj; w.bhp = W.bhp; w.open = W.open;
                p.wells.push_back(w);
            }
            pts.push_back(p);
        }
        r.stats.timesteps = st.steps; r.stats.newton = st.newton; r.stats.linear = st.linear;
    }
    r.kpi = co.kpi;
    r.oip0 = pts.empty() ? 0 : pts.front().oip;
    r.ok = ok && !pts.empty() && keepGoing();
    if (!r.ok) r.error = !keepGoing() ? "cancelled" : err.empty() ? "run failed" : err;
    return r;
}

RunResult runSingle(const Case& c, std::atomic<bool>* cancel) {
    if (c.thermal.enabled || c.comp.enabled) {
        try { return runSingleDeckEngine(c, cancel); }
        catch (const std::exception& e) { RunResult r; r.error = e.what(); return r; }
    }
    RunResult r;
    try {
        Simulator s;
        s.cancel = cancel;
        s.setup(c);
        s.initialize();
        if (coupledbo::needed(c)) {   // C27: studies evaluate the coupled mechanics too (seismicity / uplift responses)
            coupledbo::Output co;
            coupledbo::run(s, c, co, {}, cancel);
            r.kpi = co.kpi;
        } else
            while (!s.finished()) {
                if (!s.runReportStep()) break;
            }
        r.summary = s.summary();
        r.stats = s.stats();
        r.oip0 = r.summary.empty() ? 0 : r.summary.front().oip;
        r.ok = s.finished();
        if (!r.ok) r.error = "cancelled";
    } catch (const std::exception& e) { r.ok = false; r.error = e.what(); }
    return r;
}

namespace { BatchRunner g_batchRunner; std::mutex g_batchMx; }
void setBatchRunner(BatchRunner r) { std::lock_guard<std::mutex> lk(g_batchMx); g_batchRunner = std::move(r); }
bool hasBatchRunner() { std::lock_guard<std::mutex> lk(g_batchMx); return bool(g_batchRunner); }

std::vector<RunResult> runBatch(const std::vector<Case>& cases, int workers, std::atomic<bool>* cancel,
                                const std::function<void(int, int)>& progress) {
    const int n = int(cases.size());
    std::vector<RunResult> out(n);
    if (n == 0) return out;
    {
        BatchRunner ext;
        { std::lock_guard<std::mutex> lk(g_batchMx); ext = g_batchRunner; }
        if (ext && n > 1) {                           // C16: MPI cluster / external scheduler
            std::vector<RunResult> r = ext(cases, cancel, progress);
            if (int(r.size()) == n) return r;
            // runner failed as a whole (e.g. mpiexec not found): fall through to the local thread pool
        }
    }
    int hw = int(std::thread::hardware_concurrency()); if (hw < 1) hw = 1;
    int nw = workers > 0 ? workers : hw;
    nw = std::max(1, std::min(nw, n));
    std::atomic<int> next{0}, done{0};
    std::mutex mx;
    auto work = [&]() {
        for (;;) {
            const int i = next++;
            if (i >= n) break;
            if (cancel && cancel->load()) { out[i].error = "cancelled"; continue; }
            Case c = cases[i];
            c.tuning.threads = (nw >= hw) ? 1 : std::max(1, hw / nw);
            out[i] = runSingle(c, cancel);
            const int d = ++done;
            if (progress) { std::lock_guard<std::mutex> lk(mx); progress(d, n); }
        }
    };
    if (nw == 1) work();
    else {
        std::vector<std::thread> th;
        for (int t = 0; t < nw; ++t) th.emplace_back(work);
        for (auto& t : th) t.join();
    }
    return out;
}

// ======================================================================== observations
static int wellIndex(const Case& c, const std::string& name) {
    for (int i = 0; i < int(c.wells.size()); ++i) if (c.wells[i].name == name) return i;
    return -1;
}
static double pointValue(const SummaryPoint& s, ObsQty q, int w) {
    if (obsIsWell(q)) {
        if (w < 0 || w >= int(s.wells.size())) return kNaN;
        const WellRates& r = s.wells[w];
        switch (q) {
            case ObsQty::WOPR: return r.oil;
            case ObsQty::WWPR: return r.water;
            case ObsQty::WGPR: return r.gas;
            case ObsQty::WWIR: return r.waterInj;
            case ObsQty::WBHP: return r.bhp;
            case ObsQty::WWCT: return (r.oil + r.water) > 0 ? r.water / (r.oil + r.water) : 0.0;
            case ObsQty::WGOR: return r.oil > 0 ? r.gas / r.oil : 0.0;
            default: return kNaN;
        }
    }
    switch (q) {
        case ObsQty::FOPR: return s.fopr;
        case ObsQty::FWPR: return s.fwpr;
        case ObsQty::FGPR: return s.fgpr;
        case ObsQty::FPR: return s.fpr;
        case ObsQty::FWCT: return s.fwct;
        case ObsQty::FGOR: return s.fgor;
        default: return kNaN;
    }
}

double simulatedValue(const RunResult& r, const Case& c, ObsQty q, const std::string& well, double t) {
    if (r.summary.empty()) return kNaN;
    const int w = obsIsWell(q) ? wellIndex(c, well) : -1;
    if (obsIsWell(q) && w < 0) return kNaN;
    const auto& S = r.summary;
    // first point with time >= t
    size_t k = 0;
    while (k < S.size() && S[k].time < t - 1e-6) ++k;
    if (k >= S.size()) return pointValue(S.back(), q, w);
    const bool smooth = q == ObsQty::FPR || q == ObsQty::WBHP;
    if (smooth && k > 0 && S[k].time > S[k - 1].time) {
        const double a = (t - S[k - 1].time) / (S[k].time - S[k - 1].time);
        return (1 - a) * pointValue(S[k - 1], q, w) + a * pointValue(S[k], q, w);
    }
    if (k == 0 && S.size() > 1 && !smooth) return pointValue(S[1], q, w);   // rates at t=0: first step
    return pointValue(S[k], q, w);
}

std::vector<double> autoSigma(const std::vector<Observation>& obs) {
    std::map<std::string, double> mx;
    auto key = [](const Observation& o) { return std::string(toString(o.qty)) + ":" + o.well; };
    for (auto& o : obs) { double& m = mx[key(o)]; m = std::max(m, std::fabs(o.value)); }
    std::vector<double> s(obs.size());
    for (size_t i = 0; i < obs.size(); ++i) {
        const Observation& o = obs[i];
        double v;
        if (o.sigma > 0) v = o.sigma;
        else {
            v = std::max(0.05 * std::fabs(o.value), 0.03 * mx[key(o)]);
            if (o.qty == ObsQty::WWCT || o.qty == ObsQty::FWCT) v = std::max(v, 0.02);
            v = std::max(v, 1e-6);
        }
        // relative weight: w -> sigma / sqrt(w); w = 0 removes the observation from the objective
        s[i] = o.weight > 0 ? v / std::sqrt(o.weight) : 1e30;
    }
    return s;
}

std::vector<double> residuals(const RunResult& r, const Case& c, const std::vector<Observation>& obs, const std::vector<double>& sigma) {
    std::vector<double> res(obs.size(), 0.0);
    for (size_t i = 0; i < obs.size(); ++i) {
        const double s = r.ok ? simulatedValue(r, c, obs[i].qty, obs[i].well, obs[i].time) : kNaN;
        res[i] = std::isfinite(s) ? (s - obs[i].value) / sigma[i] : 30.0;
    }
    return res;
}

double misfit(const RunResult& r, const Case& c, const std::vector<Observation>& obs, const std::vector<double>& sigma) {
    if (obs.empty()) return 0;
    const auto res = residuals(r, c, obs, sigma);
    double s = 0; for (double v : res) s += v * v;
    return s / double(obs.size());
}

std::vector<Observation> observationsFromHistory(const Case& c) {
    std::vector<Observation> out;
    std::map<std::string, WellControl> cur;
    double t = 0;
    for (const ReportStep& s : c.schedule) {
        for (const WellControl& wc : s.controls)
            for (const WellInput& w : c.wells)
                if (w.name == wc.well || (!wc.well.empty() && wc.well.back() == '*' && w.name.compare(0, wc.well.size() - 1, wc.well, 0, wc.well.size() - 1) == 0))
                    cur[w.name] = wc;
        t += s.length;
        for (const WellInput& w : c.wells) {
            auto it = cur.find(w.name);
            if (it == cur.end() || !it->second.history || w.type != WellType::Producer || it->second.status != WellStatus::Open) continue;
            const WellControl& wc = it->second;
            if (wc.orat >= 0) out.push_back({ObsQty::WOPR, w.name, t, wc.orat, -1});
            if (wc.wrat >= 0 && c.phases.water) out.push_back({ObsQty::WWPR, w.name, t, wc.wrat, -1});
            if (wc.grat >= 0 && c.phases.gas) out.push_back({ObsQty::WGPR, w.name, t, wc.grat, -1});
        }
    }
    return out;
}

namespace {
long civilDays(int y, int m, int d) {
    y -= m <= 2;
    const long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = unsigned(y - era * 400);
    const unsigned doy = (153 * unsigned(m + (m > 2 ? -3 : 9)) + 2) / 5 + unsigned(d) - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + long(doe) - 719468;
}
}  // namespace

std::vector<Observation> observationsFromSummary(const ecl::Summary& s, const Case& c, const SummaryObsOptions& opt, std::string* note) {
    static const char* kw[] = {"WOPR", "WWPR", "WGPR", "WWIR", "WBHP", "WWCT", "WGOR", "FOPR", "FWPR", "FGPR", "FPR", "FWCT", "FGOR"};
    std::vector<Observation> out;
    const double shift = double(civilDays(s.startYear, s.startMonth, s.startDay) - civilDays(c.startYear, c.startMonth, c.startDay));
    std::vector<int> steps = s.reportIndex;
    if (steps.empty()) for (int i = 0; i < int(s.time.size()); ++i) steps.push_back(i);
    std::map<std::string, int> perQty;
    std::set<std::string> seen;
    int unknownWells = 0;
    auto known = [&](const std::string& w) {
        if (c.wells.empty()) return true;
        for (auto& x : c.wells) if (x.name == w) return true;
        return false;
    };
    for (int q = 0; q < 13; ++q) {
        if (!((opt.qtyMask >> q) & 1u)) continue;
        const ObsQty qty = ObsQty(q);
        const bool well = kw[q][0] == 'W';
        std::string name = kw[q];
        if (opt.history) { if (qty == ObsQty::FPR) continue; name += 'H'; }
        for (int v = 0; v < int(s.vars.size()); ++v) {
            const ecl::SummaryVar& sv = s.vars[size_t(v)];
            if (sv.keyword != name) continue;
            if (well && (sv.wgname.empty() || sv.wgname == ":+:+:+:+")) continue;
            if (well && !known(sv.wgname)) { ++unknownWells; continue; }
            if (!seen.insert(std::to_string(q) + ":" + sv.wgname).second) continue;   // duplicated vector in the SMSPEC
            // rate of the well / field used to skip shut periods (rates, ratios and BHP of a shut well are meaningless)
            int rateVar = -1;
            if (opt.skipShut) {
                std::string rk = well ? (qty == ObsQty::WWIR ? "WWIR" : "WLPR") : "FLPR";
                if (opt.history) rk += 'H';
                rateVar = s.find(well ? rk + ":" + sv.wgname : rk);
                if (rateVar < 0 && well && qty != ObsQty::WWIR) {      // no liquid-rate vector: oil + water
                    rateVar = s.find(std::string(opt.history ? "WOPRH:" : "WOPR:") + sv.wgname);
                }
            }
            double lastT = -1e300;
            for (int st : steps) {
                if (st < 0 || st >= int(s.time.size())) continue;
                const double t = s.time[size_t(st)] + shift;
                if (t <= 1e-9 || t - lastT < opt.minSpacingDays) continue;
                const double val = s.values[size_t(v)][size_t(st)];
                if (!std::isfinite(val)) continue;
                if (rateVar >= 0) {
                    double r = s.values[size_t(rateVar)][size_t(st)];
                    if (rateVar >= 0 && s.vars[size_t(rateVar)].keyword.compare(0, 4, "WOPR") == 0) {
                        const int wv = s.find(std::string(opt.history ? "WWPRH:" : "WWPR:") + sv.wgname);
                        if (wv >= 0) r += s.values[size_t(wv)][size_t(st)];
                    }
                    if (!(r > 0)) continue;
                }
                if (qty == ObsQty::WBHP && !(val > 0)) continue;
                Observation o; o.qty = qty; o.well = well ? sv.wgname : ""; o.time = t; o.value = val;
                out.push_back(o);
                lastT = t;
                perQty[kw[q]]++;
            }
        }
    }
    if (note) {
        std::string n = std::to_string(out.size()) + " observations from " + (opt.history ? "history" : "simulated") + " vectors";
        if (!perQty.empty()) { n += " ("; bool f = true; for (auto& [k, m] : perQty) { n += (f ? "" : ", ") + k + " " + std::to_string(m); f = false; } n += ")"; }
        if (unknownWells) n += "; " + std::to_string(unknownWells) + " vectors of wells not in the model skipped";
        if (shift != 0) n += "; times shifted by " + std::to_string(long(shift)) + " days to the model start date";
        const char* cu = c.units == UnitSystem::Metric ? "METRIC" : "FIELD";
        if (!s.unitSystem.empty() && s.unitSystem != cu) n += "; WARNING: summary units " + s.unitSystem + " differ from the model (" + cu + ")";
        *note = n;
    }
    return out;
}

std::vector<Observation> syntheticObservations(const RunResult& truth, const Case& c, const std::vector<ObsQty>& qtys,
                                               double relNoise, unsigned seed) {
    std::vector<Observation> out;
    if (!truth.ok) return out;
    std::mt19937 rng(seed); std::normal_distribution<double> N(0.0, 1.0);
    double t = 0;
    for (const ReportStep& s : c.schedule) {
        t += s.length;
        for (ObsQty q : qtys) {
            if (obsIsWell(q)) {
                for (const WellInput& w : c.wells) {
                    const bool prod = w.type == WellType::Producer;
                    if (q == ObsQty::WWIR ? (w.type != WellType::WaterInjector) : (q != ObsQty::WBHP && !prod)) continue;
                    const double v = simulatedValue(truth, c, q, w.name, t);
                    if (!std::isfinite(v)) continue;
                    out.push_back({q, w.name, t, v * (1.0 + relNoise * N(rng)), -1});
                }
            } else {
                const double v = simulatedValue(truth, c, q, "", t);
                if (std::isfinite(v)) out.push_back({q, "", t, v * (1.0 + relNoise * N(rng)), -1});
            }
        }
    }
    return out;
}

// ======================================================================== responses
const char* toString(Response r) {
    switch (r) {
        case Response::FOPT: return "Cumulative oil (FOPT)";
        case Response::RF: return "Recovery factor";
        case Response::FWPT: return "Cumulative water (FWPT)";
        case Response::FGPT: return "Cumulative gas (FGPT)";
        case Response::FWIT: return "Cumulative water injected";
        case Response::FPR: return "Final average pressure";
        case Response::FWCT: return "Final water cut";
        case Response::NPV: return "NPV";
        case Response::Misfit: return "History-match misfit";
        case Response::SeisEvents: return "Expected induced events (faults)";
        case Response::MaxMw: return "Largest magnitude Mw (faults)";
        case Response::MaxDCFS: return "Max Coulomb stress change (faults)";
        case Response::Uplift: return "Surface uplift [m] (geomechanics)";
        case Response::MaxFailure: return "Max failure index (geomechanics)";
    }
    return "?";
}

double responseValue(const RunResult& r, const Case& c, Response resp, const Economics& eco,
                     const std::vector<Observation>& obs, const std::vector<double>& sigma) {
    if (!r.ok || r.summary.empty()) return kNaN;
    const SummaryPoint& e = r.summary.back();
    switch (resp) {
        case Response::FOPT: return e.fopt;
        case Response::RF: return r.oip0 > 0 ? e.fopt / r.oip0 : 0;
        case Response::FWPT: return e.fwpt;
        case Response::FGPT: return e.fgpt;
        case Response::FWIT: return e.fwit;
        case Response::FPR: return e.fpr;
        case Response::FWCT: return e.fwct;
        case Response::NPV: {
            double npv = 0;
            for (const SummaryPoint& s : r.summary) {
                if (s.dt <= 0) continue;
                const double cash = eco.oilPrice * s.fopr + eco.gasPrice * s.fgpr - eco.waterCost * s.fwpr - eco.injCost * s.fwir;
                npv += cash * s.dt / std::pow(1.0 + eco.discount, (s.time - 0.5 * s.dt) / 365.25);
            }
            return npv;
        }
        case Response::Misfit: return misfit(r, c, obs, sigma);
        case Response::SeisEvents: return coupledbo::kpiSum(r.kpi, "_expected_events");
        case Response::MaxMw: return coupledbo::kpiMax(r.kpi, "_max_mw");
        case Response::MaxDCFS: return coupledbo::kpiMax(r.kpi, "_max_dcfs");
        case Response::Uplift: return coupledbo::kpiMax(r.kpi, "surface_uplift_m");
        case Response::MaxFailure: return coupledbo::kpiMax(r.kpi, "max_failure_index");
    }
    return kNaN;
}

// ======================================================================== helpers
namespace {
void say(const StudyCallbacks& cb, const std::string& s) { if (cb.log) cb.log(s); }
std::string f2s(const char* f, double a, double b = 0, double c = 0) { char buf[256]; std::snprintf(buf, sizeof(buf), f, a, b, c); return buf; }
bool cancelled(const StudyCallbacks& cb) { return cb.cancel && cb.cancel->load(); }
std::vector<double> baseValues(const std::vector<StudyParam>& P) { std::vector<double> v; for (auto& p : P) v.push_back(p.base); return v; }
void archive(const StudyCallbacks& cb, const char* tag, const std::vector<std::vector<double>>& vals, const std::vector<RunResult>& runs) {
    if (!cb.onRun) return;
    for (size_t i = 0; i < runs.size() && i < vals.size(); ++i)
        if (runs[i].error != "cancelled") cb.onRun(tag, vals[i], runs[i]);
}
void archive1(const StudyCallbacks& cb, const char* tag, const std::vector<double>& vals, const RunResult& run) {
    if (cb.onRun && run.error != "cancelled") cb.onRun(tag, vals, run);
}

// solve small dense system A x = b (n x n, row-major) with partial pivoting; returns false if singular
bool solveDense(std::vector<double> A, std::vector<double>& b, int n) {
    for (int k = 0; k < n; ++k) {
        int pr = k;
        for (int i = k + 1; i < n; ++i) if (std::fabs(A[i * n + k]) > std::fabs(A[pr * n + k])) pr = i;
        if (std::fabs(A[pr * n + k]) < 1e-300) return false;
        if (pr != k) { for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[pr * n + j]); std::swap(b[k], b[pr]); }
        for (int i = k + 1; i < n; ++i) {
            const double f = A[i * n + k] / A[k * n + k];
            if (f == 0) continue;
            for (int j = k; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) { double s = b[i]; for (int j = i + 1; j < n; ++j) s -= A[i * n + j] * b[j]; b[i] = s / A[i * n + i]; }
    return true;
}
double percentile(std::vector<double> v, double p) {
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
    if (v.empty()) return kNaN;
    std::sort(v.begin(), v.end());
    const double pos = p * (v.size() - 1);
    const size_t i = size_t(pos);
    const double f = pos - i;
    return i + 1 < v.size() ? v[i] * (1 - f) + v[i + 1] * f : v.back();
}
std::vector<double> ranks(const std::vector<double>& v) {
    std::vector<int> idx(v.size()); std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return v[a] < v[b]; });
    std::vector<double> r(v.size());
    for (size_t i = 0; i < idx.size();) {
        size_t j = i; while (j + 1 < idx.size() && v[idx[j + 1]] == v[idx[i]]) ++j;
        for (size_t k = i; k <= j; ++k) r[idx[k]] = 0.5 * (i + j);
        i = j + 1;
    }
    return r;
}
double pearson(const std::vector<double>& a, const std::vector<double>& b) {
    const size_t n = a.size(); if (n < 3) return 0;
    double ma = 0, mb = 0; for (size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n;
    double sab = 0, saa = 0, sbb = 0;
    for (size_t i = 0; i < n; ++i) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
    return (saa > 0 && sbb > 0) ? sab / std::sqrt(saa * sbb) : 0;
}
double fopAt(const RunResult& r, double t) {
    const auto& S = r.summary;
    if (S.empty()) return kNaN;
    if (t <= S.front().time) return S.front().fopt;
    for (size_t k = 1; k < S.size(); ++k)
        if (S[k].time >= t - 1e-9) {
            const double a = (t - S[k - 1].time) / std::max(1e-12, S[k].time - S[k - 1].time);
            return (1 - a) * S[k - 1].fopt + a * S[k].fopt;
        }
    return S.back().fopt;
}
}  // namespace

// ======================================================================== tornado
TornadoResult runTornado(const Case& base, const std::vector<StudyParam>& params, Response resp, const Economics& eco,
                         const StudyCallbacks& cb) {
    TornadoResult tr;
    const std::vector<Observation>& obs = base.observations;
    const std::vector<double> sig = autoSigma(obs);
    std::vector<Case> cases;
    std::vector<int> pidx;
    const std::vector<double> bv = baseValues(params);
    std::vector<std::vector<double>> vals{bv};
    cases.push_back(applyParams(base, params, bv));
    for (int i = 0; i < int(params.size()); ++i) {
        if (!params[i].active) continue;
        std::vector<double> lo = bv, hi = bv;
        lo[i] = params[i].low; hi[i] = params[i].high;
        cases.push_back(applyParams(base, params, lo));
        cases.push_back(applyParams(base, params, hi));
        vals.push_back(lo); vals.push_back(hi);
        pidx.push_back(i);
    }
    say(cb, f2s("Tornado: %g runs (%g parameters, low/high)", double(cases.size()), double(pidx.size())));
    auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
    archive(cb, "tornado", vals, runs);
    tr.baseRun = runs[0];
    tr.ok = runs[0].ok;
    tr.base = responseValue(runs[0], cases[0], resp, eco, obs, sig);
    for (size_t k = 0; k < pidx.size(); ++k) {
        TornadoBar b; b.param = pidx[k];
        const RunResult& rl = runs[1 + 2 * k]; const RunResult& rh = runs[2 + 2 * k];
        b.okLow = rl.ok; b.okHigh = rh.ok;
        b.low = responseValue(rl, cases[1 + 2 * k], resp, eco, obs, sig);
        b.high = responseValue(rh, cases[2 + 2 * k], resp, eco, obs, sig);
        if (!rl.ok) say(cb, "  low run of " + paramLabel(params[b.param], base) + " failed: " + rl.error);
        if (!rh.ok) say(cb, "  high run of " + paramLabel(params[b.param], base) + " failed: " + rh.error);
        // C28: a side identical to the base usually means an active constraint (e.g. a rate multiplier on a well that runs on
        // its BHP / THP limit - FAULT_SLIP: injector at the 280 bar limit, so x1.5 changes nothing) - say so instead of hiding it
        {
            const ParamKind pk = params[b.param].kind;
            const bool rate = pk == ParamKind::InjRate || pk == ParamKind::ProdRate;
            auto same = [&](double v) { return std::isfinite(v) && std::isfinite(tr.base) && std::fabs(v - tr.base) <= 1e-9 * std::max(1.0, std::fabs(tr.base)); };
            if (rl.ok && rh.ok && (same(b.low) || same(b.high)))
                say(cb, "  NOTE: " + paramLabel(params[b.param], base) + (same(b.high) && !same(b.low) ? ": the high side" : same(b.low) && !same(b.high) ? ": the low side" : ": both sides") +
                        " equal the base" + (rate ? " - the wells are probably on their BHP / THP limit (the rate target is not the active control)" : " - the response is insensitive in that range"));
        }
        tr.bars.push_back(b);
    }
    std::sort(tr.bars.begin(), tr.bars.end(), [&](const TornadoBar& a, const TornadoBar& b) {
        auto sw = [&](const TornadoBar& x) { const double l = std::isfinite(x.low) ? x.low : tr.base, h = std::isfinite(x.high) ? x.high : tr.base; return std::fabs(h - l); };
        return sw(a) > sw(b);
    });
    return tr;
}

// ======================================================================== Monte Carlo (LHS)
MonteCarloResult runMonteCarlo(const Case& base, const std::vector<StudyParam>& params, int ns, unsigned seed,
                               Response resp, const Economics& eco, const StudyCallbacks& cb) {
    MonteCarloResult mc;
    ns = std::max(2, ns);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const int np = int(params.size());
    mc.samples.assign(ns, baseValues(params));
    for (int i = 0; i < np; ++i) {
        if (!params[i].active) continue;
        std::vector<int> perm(ns); std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int s = 0; s < ns; ++s) mc.samples[s][i] = paramFromUnit(params[i], (perm[s] + U(rng)) / ns);
    }
    std::vector<Case> cases;
    for (int s = 0; s < ns; ++s) cases.push_back(applyParams(base, params, mc.samples[s]));
    say(cb, f2s("Monte Carlo: %g Latin-hypercube samples", double(ns)));
    auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
    archive(cb, "MC", mc.samples, runs);
    const std::vector<double> sig = autoSigma(base.observations);
    mc.response.resize(ns);
    for (int s = 0; s < ns; ++s) {
        mc.response[s] = responseValue(runs[s], cases[s], resp, eco, base.observations, sig);
        if (!runs[s].ok) ++mc.failed;
    }
    // industry convention: P90 = low (exceeded with 90% probability), P10 = high
    mc.p90 = percentile(mc.response, 0.10); mc.p50 = percentile(mc.response, 0.50); mc.p10 = percentile(mc.response, 0.90);
    double sum = 0; int cnt = 0; for (double v : mc.response) if (std::isfinite(v)) { sum += v; ++cnt; }
    mc.mean = cnt ? sum / cnt : kNaN;
    // Spearman rank correlation per parameter
    mc.spearman.assign(np, 0.0);
    std::vector<int> okIdx; for (int s = 0; s < ns; ++s) if (std::isfinite(mc.response[s])) okIdx.push_back(s);
    std::vector<double> rv; for (int s : okIdx) rv.push_back(mc.response[s]);
    const auto rr = ranks(rv);
    for (int i = 0; i < np; ++i) {
        if (!params[i].active) continue;
        std::vector<double> pv; for (int s : okIdx) pv.push_back(mc.samples[s][i]);
        mc.spearman[i] = pearson(ranks(pv), rr);
    }
    // FOPT fan chart at report times
    double t = 0; mc.times.push_back(0);
    for (auto& st : base.schedule) { t += st.length; mc.times.push_back(t); }
    for (double tt : mc.times) {
        std::vector<double> v; for (int s : okIdx) v.push_back(fopAt(runs[s], tt));
        mc.fanP90.push_back(percentile(v, 0.10)); mc.fanP50.push_back(percentile(v, 0.5)); mc.fanP10.push_back(percentile(v, 0.90));
    }
    say(cb, f2s("  P90 %.6g  P50 %.6g  P10 %.6g", mc.p90, mc.p50, mc.p10) + f2s("  (%g failed)", double(mc.failed)));
    return mc;
}

// ======================================================================== analytic gradients (discrete adjoint / forward sensitivities)
const char* toString(GradientMethod g) {
    switch (g) {
        case GradientMethod::FiniteDifference: return "Finite differences";
        case GradientMethod::Analytic: return "Analytic (forward sensitivities)";
    }
    return "?";
}

namespace {
// a forward run with adjoint recording at normalised active coordinates u
struct RecordedRun {
    std::unique_ptr<Simulator> sim;
    Case c;
    RunResult run;
    std::vector<double> values;
};
RecordedRun recordedRun(const Case& base, const std::vector<StudyParam>& params, const std::vector<int>& act,
                        const std::vector<double>& u, std::atomic<bool>* cancel) {
    RecordedRun R;
    R.values = baseValues(params);
    for (size_t k = 0; k < act.size(); ++k) R.values[size_t(act[k])] = paramFromUnit(params[size_t(act[k])], u[k]);
    R.c = applyParams(base, params, R.values);
    try {
        R.sim.reset(new Simulator);
        R.sim->cancel = cancel;
        R.sim->setAdjointRecording(true);
        R.sim->setup(R.c);
        R.sim->initialize();
        while (!R.sim->finished()) if (!R.sim->runReportStep()) break;
        R.run.summary = R.sim->summary();
        R.run.stats = R.sim->stats();
        R.run.oip0 = R.run.summary.empty() ? 0 : R.run.summary.front().oip;
        R.run.ok = R.sim->finished();
        if (!R.run.ok) R.run.error = "cancelled";
    } catch (const std::exception& e) { R.run.ok = false; R.run.error = e.what(); }
    return R;
}
// adjoint (forward = false) or forward-sensitivity evaluation of a recorded run; gradient in normalised coordinates.
// Observations with zero weight or of unknown wells are constant in the misfit and are left out.
struct AnalyticOut { bool ok = false; std::string error; std::vector<double> grad; std::vector<std::vector<double>> jac; double replayDev = 0, maxLinRes = 0; };
AnalyticOut analyticGradient(RecordedRun& R, const Case& base, const std::vector<StudyParam>& params, const std::vector<int>& act,
                             const std::vector<double>& u, const std::vector<double>& sig, bool forward, double hres) {
    AnalyticOut out;
    const auto& obs = base.observations;
    const int nd = int(obs.size()), n = int(act.size());
    std::vector<Simulator::AdjObs> ao; std::vector<int> aIdx;
    for (int i = 0; i < nd; ++i) {
        if (!(sig[size_t(i)] < 1e29)) continue;
        Simulator::AdjObs o; o.qty = obs[size_t(i)].qty; o.time = obs[size_t(i)].time; o.value = obs[size_t(i)].value; o.sigma = sig[size_t(i)];
        if (obsIsWell(o.qty)) { o.well = wellIndex(R.c, obs[size_t(i)].well); if (o.well < 0) continue; }
        ao.push_back(o); aIdx.push_back(i);
    }
    // C11: permeability / porosity / MULTZ box multipliers and well PI multipliers get analytic dR/dtheta from ONE reverse
    // pass (field gradients); the others use residual differences of perturbed models (cost per parameter)
    const Grid& grid = R.sim->grid();
    auto isField = [&](ParamKind k) {
        if (forward) return false;
        if (k == ParamKind::WellPI) return true;
        if (k == ParamKind::Poro) return grid.lgr.empty();
        return (k == ParamKind::PermXY || k == ParamKind::PermZ || k == ParamKind::MultZ) && grid.fieldDerivOk;
    };
    Simulator::AdjFieldReq fr;
    std::vector<int> repl;   // positions k (in act) using perturbed models
    for (int k = 0; k < n; ++k) {
        const ParamKind pk = params[size_t(act[size_t(k)])].kind;
        if (!isField(pk)) { repl.push_back(k); continue; }
        if (pk == ParamKind::PermXY || pk == ParamKind::PermZ) fr.perm = true;
        else if (pk == ParamKind::Poro) fr.poro = true;
        else if (pk == ParamKind::MultZ) fr.trans = true;
        else fr.wi = true;
    }
    const bool anyField = fr.perm || fr.poro || fr.trans || fr.wi;
    const int nr = int(repl.size());
    std::vector<Case> plus(static_cast<size_t>(nr)), minus(static_cast<size_t>(nr));
    std::vector<Simulator::AdjParam> ap(static_cast<size_t>(nr));
    for (int kr = 0; kr < nr; ++kr) {
        const int k = repl[size_t(kr)];
        std::vector<double> up = u, um = u;
        up[size_t(k)] = std::min(1.0, u[size_t(k)] + hres); um[size_t(k)] = std::max(0.0, u[size_t(k)] - hres);
        std::vector<double> vp = R.values, vm = R.values;
        vp[size_t(act[size_t(k)])] = paramFromUnit(params[size_t(act[size_t(k)])], up[size_t(k)]);
        vm[size_t(act[size_t(k)])] = paramFromUnit(params[size_t(act[size_t(k)])], um[size_t(k)]);
        plus[size_t(kr)] = applyParams(base, params, vp); minus[size_t(kr)] = applyParams(base, params, vm);
        ap[size_t(kr)].plus = &plus[size_t(kr)]; ap[size_t(kr)].minus = &minus[size_t(kr)]; ap[size_t(kr)].h = 0.5 * (up[size_t(k)] - um[size_t(k)]);
    }
    const Simulator::AdjResult a = R.sim->adjointGradient(ao, ap, forward, anyField ? &fr : nullptr);
    if (!a.ok) { out.error = a.error; return out; }
    out.ok = true; out.replayDev = a.replayDev; out.maxLinRes = a.maxLinRes;
    out.grad.assign(size_t(n), 0.0);
    for (int kr = 0; kr < nr; ++kr) out.grad[size_t(repl[size_t(kr)])] = a.grad[size_t(kr)] / std::max(1, nd);   // misfit = F / nd
    if (anyField)
        for (int k = 0; k < n; ++k) {
            const StudyParam& P = params[size_t(act[size_t(k)])];
            if (!isField(P.kind)) continue;
            const double v = R.values[size_t(act[size_t(k)])];
            if (!(v > 0)) continue;
            // dF / d ln v summed over the parameter's cells / faces / connections
            double gl = 0;
            const GridInput& g = R.c.grid;
            const int i1 = P.i1 < 0 ? 0 : P.i1, i2 = P.i2 < 0 ? g.nx - 1 : std::min(P.i2, g.nx - 1);
            const int j1 = P.j1 < 0 ? 0 : P.j1, j2 = P.j2 < 0 ? g.ny - 1 : std::min(P.j2, g.ny - 1);
            const int k1 = P.k1 < 0 ? 0 : P.k1, k2 = P.k2 < 0 ? g.nz - 1 : std::min(P.k2, g.nz - 1);
            auto inBox = [&](int q) {
                if (q < 0 || q >= g.nx * g.ny * g.nz) return false;
                const int ii = q % g.nx, jj = (q / g.nx) % g.ny, kk = q / (g.nx * g.ny);
                return ii >= i1 && ii <= i2 && jj >= j1 && jj <= j2 && kk >= k1 && kk <= k2;
            };
            if (P.kind == ParamKind::WellPI) {
                if (P.target >= 0 && P.target < int(a.gLogWI.size())) for (double x : a.gLogWI[size_t(P.target)]) gl += x;
            } else if (P.kind == ParamKind::MultZ) {
                for (size_t f = 0; f < grid.conns.size() && f < a.gLogTrans.size(); ++f)
                    if (grid.conns[f].dir == 'Z' && inBox(grid.conns[f].mq)) gl += a.gLogTrans[f];
            } else {
                for (int q = 0; q < g.nx * g.ny * g.nz; ++q) {
                    if (!inBox(q)) continue;
                    const int ac = grid.activeOf[size_t(q)];
                    if (ac < 0) continue;
                    if (P.kind == ParamKind::PermXY) gl += a.gLogKx[size_t(ac)] + a.gLogKy[size_t(ac)];
                    else if (P.kind == ParamKind::PermZ) gl += a.gLogKz[size_t(ac)];
                    else if (g.poro[size_t(q)] < 1.0) gl += a.gLogPoro[size_t(ac)];   // (poro * v clamped at 1)
                }
            }
            // chain rule to the normalised coordinate: d ln v / du
            const double uk = u[size_t(k)], hu = 1e-6;
            const double ua = std::max(0.0, uk - hu), ub = std::min(1.0, uk + hu);
            const double dvdu = (paramFromUnit(P, ub) - paramFromUnit(P, ua)) / (ub - ua);
            out.grad[size_t(k)] = gl * dvdu / v / std::max(1, nd);
        }
    if (forward) {   // Jacobian of the normalised residuals r_i = (y_i - d_i) / sigma_i
        out.jac.assign(size_t(nd), std::vector<double>(size_t(n), 0.0));
        for (size_t q = 0; q < aIdx.size(); ++q)
            for (int k = 0; k < n; ++k) out.jac[size_t(aIdx[q])][size_t(k)] = a.dydp[q][size_t(k)] / ao[q].sigma;
    }
    return out;
}
} // namespace

GradientCheck misfitGradient(const Case& base, const std::vector<StudyParam>& params, int methods, double fdStep,
                             const StudyCallbacks& cb, double resStep) {
    GradientCheck gc;
    const auto& obs = base.observations;
    if (obs.empty()) { gc.error = "no observed data"; return gc; }
    const std::vector<double> sig = autoSigma(obs);
    for (int i = 0; i < int(params.size()); ++i) if (params[size_t(i)].active) gc.act.push_back(i);
    const int n = int(gc.act.size());
    if (n == 0) { gc.error = "no active parameters"; return gc; }
    std::vector<double> u(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) u[size_t(k)] = paramToUnit(params[size_t(gc.act[size_t(k)])], params[size_t(gc.act[size_t(k)])].base);
    auto clock = [] { return std::chrono::steady_clock::now(); };
    auto secs = [](std::chrono::steady_clock::time_point a) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count(); };
    if (methods & 3) {
        for (int m = 0; m < 2; ++m) {
            if (!(methods & (1 << m))) continue;
            const auto t0 = clock();
            RecordedRun R = recordedRun(base, params, gc.act, u, cb.cancel);
            if (!R.run.ok) { gc.error = "base run failed: " + R.run.error; return gc; }
            gc.misfit = misfit(R.run, R.c, obs, sig);
            AnalyticOut a = analyticGradient(R, base, params, gc.act, u, sig, m == 1, resStep);
            if (!a.ok) { gc.error = a.error; return gc; }
            (m == 0 ? gc.adjoint : gc.forward) = a.grad;
            (m == 0 ? gc.secAdjoint : gc.secForward) = secs(t0);
            gc.replayDev = std::max(gc.replayDev, a.replayDev); gc.maxLinRes = std::max(gc.maxLinRes, a.maxLinRes);
            say(cb, f2s(m == 0 ? "  adjoint gradient: %.2f s (replay deviation %.1e, linear residual %.1e)"
                               : "  forward sensitivities: %.2f s (replay deviation %.1e, linear residual %.1e)", secs(t0), a.replayDev, a.maxLinRes));
        }
    }
    if (methods & 4) {
        const auto t0 = clock();
        std::vector<Case> cases; std::vector<double> du(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k)
            for (int s = 0; s < 2; ++s) {
                std::vector<double> v = baseValues(params);
                for (int q = 0; q < n; ++q) v[size_t(gc.act[size_t(q)])] = paramFromUnit(params[size_t(gc.act[size_t(q)])], u[size_t(q)]);
                const double uk = s == 0 ? std::min(1.0, u[size_t(k)] + fdStep) : std::max(0.0, u[size_t(k)] - fdStep);
                if (s == 0) du[size_t(k)] = uk; else du[size_t(k)] -= uk;
                v[size_t(gc.act[size_t(k)])] = paramFromUnit(params[size_t(gc.act[size_t(k)])], uk);
                cases.push_back(applyParams(base, params, v));
            }
        auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
        gc.fdRuns = int(runs.size());
        gc.fd.resize(size_t(n));
        for (int k = 0; k < n; ++k) {
            const RunResult& a = runs[size_t(2 * k)]; const RunResult& b = runs[size_t(2 * k + 1)];
            gc.fd[size_t(k)] = (a.ok && b.ok) ? (misfit(a, cases[size_t(2 * k)], obs, sig) - misfit(b, cases[size_t(2 * k + 1)], obs, sig)) / du[size_t(k)] : kNaN;
        }
        gc.secFd = secs(t0);
        if (!(methods & 3)) { RunResult r0 = runSingle(applyParams(base, params, baseValues(params)), cb.cancel); gc.misfit = misfit(r0, base, obs, sig); }
    }
    gc.ok = true;
    return gc;
}

// projected L-BFGS on [0,1]^n driven by the adjoint gradient (Byrd-Lu-Nocedal-Zhu style projection, Armijo backtracking)
HMResult runAdjointLBFGS(const Case& base, const std::vector<StudyParam>& params, const AdjointHMOptions& opt,
                         const StudyCallbacks& cb) {
    HMResult hm;
    const auto tStart = std::chrono::steady_clock::now();
    const auto& obs = base.observations;
    if (obs.empty()) { say(cb, "History matching needs observed data"); return hm; }
    const std::vector<double> sig = autoSigma(obs);
    std::vector<int> act;
    for (int i = 0; i < int(params.size()); ++i) if (params[size_t(i)].active) act.push_back(i);
    const int n = int(act.size());
    if (n == 0) { say(cb, "No active parameters"); return hm; }
    std::vector<double> u(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) u[size_t(k)] = paramToUnit(params[size_t(act[size_t(k)])], params[size_t(act[size_t(k)])].base);
    auto evalRun = [&](const std::vector<double>& uu, double& f) {
        RecordedRun R = recordedRun(base, params, act, uu, cb.cancel);
        hm.runs++;
        archive1(cb, "adjoint", R.values, R.run);
        f = R.run.ok ? misfit(R.run, R.c, obs, sig) : std::numeric_limits<double>::infinity();
        return R;
    };
    auto gradOf = [&](RecordedRun& R, const std::vector<double>& uu, std::vector<double>& g) {
        AnalyticOut a = analyticGradient(R, base, params, act, uu, sig, false, opt.resStep);
        if (!a.ok) { say(cb, "  adjoint failed: " + a.error); return false; }
        g = a.grad;
        return true;
    };
    double f = 0;
    RecordedRun R0 = evalRun(u, f);
    if (!R0.run.ok) { say(cb, "Base run failed: " + R0.run.error); return hm; }
    std::vector<double> g;
    if (!gradOf(R0, u, g)) return hm;
    hm.initialRun = R0.run; hm.initialMisfit = f; hm.bestMisfit = f; hm.best = R0.values; hm.bestRun = R0.run;
    // C13: full iteration record (projected-gradient norm, cumulative runs / time, gradient per parameter, the run)
    auto record = [&](int it, double fv, double alpha, const std::vector<double>& gu, const RunResult& run) {
        HMIteration h; h.iter = it; h.misfit = fv; h.lambda = alpha; h.values = hm.best; h.runs = hm.runs;
        h.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
        h.grad.assign(params.size(), std::numeric_limits<double>::quiet_NaN());
        double gn = 0;
        for (int k = 0; k < n && size_t(k) < gu.size(); ++k) { h.grad[size_t(act[size_t(k)])] = gu[size_t(k)]; gn += gu[size_t(k)] * gu[size_t(k)]; }
        h.gradNorm = std::sqrt(gn); h.run = run;
        hm.history.push_back(std::move(h));
    };
    record(0, f, 0, g, R0.run);
    say(cb, f2s("Adjoint L-BFGS: %g parameters, %g observations, initial misfit %.5g", n, double(obs.size()), f));
    R0.sim.reset();
    std::vector<std::vector<double>> S, Y;
    auto dot = [](const std::vector<double>& a, const std::vector<double>& b) { double s = 0; for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i]; return s; };
    double pg0 = -1;
    for (int it = 1; it <= opt.maxIter && !cancelled(cb); ++it) {
        std::vector<char> fr(size_t(n), 1);
        double pgn = 0;
        for (int k = 0; k < n; ++k) {
            if ((u[size_t(k)] <= 0 && g[size_t(k)] > 0) || (u[size_t(k)] >= 1 && g[size_t(k)] < 0)) fr[size_t(k)] = 0;
            else pgn += g[size_t(k)] * g[size_t(k)];
        }
        pgn = std::sqrt(pgn);
        if (pg0 < 0) pg0 = pgn;
        if (pgn <= 1e-10 || pgn < opt.tol * std::max(pg0, 1e-300) * 1e-2) { say(cb, "  converged (projected gradient)"); break; }
        // two-loop recursion on the free variables
        std::vector<double> q(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k) q[size_t(k)] = fr[size_t(k)] ? g[size_t(k)] : 0.0;
        const int m = int(S.size());
        std::vector<double> al(static_cast<size_t>(m));
        for (int i = m - 1; i >= 0; --i) { al[size_t(i)] = dot(S[size_t(i)], q) / dot(Y[size_t(i)], S[size_t(i)]); for (int k = 0; k < n; ++k) q[size_t(k)] -= al[size_t(i)] * Y[size_t(i)][size_t(k)]; }
        const double gam = m > 0 ? dot(S.back(), Y.back()) / dot(Y.back(), Y.back()) : 1.0;
        for (double& v : q) v *= gam;
        for (int i = 0; i < m; ++i) { const double be = dot(Y[size_t(i)], q) / dot(Y[size_t(i)], S[size_t(i)]); for (int k = 0; k < n; ++k) q[size_t(k)] += S[size_t(i)][size_t(k)] * (al[size_t(i)] - be); }
        std::vector<double> d(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k) d[size_t(k)] = fr[size_t(k)] ? -q[size_t(k)] : 0.0;
        if (dot(d, g) >= 0) { for (int k = 0; k < n; ++k) d[size_t(k)] = fr[size_t(k)] ? -g[size_t(k)] : 0.0; S.clear(); Y.clear(); }
        double dmax = 0; for (double v : d) dmax = std::max(dmax, std::fabs(v));
        const double cap = S.empty() ? opt.maxStep : 2.0 * opt.maxStep;
        if (dmax > cap) for (double& v : d) v *= cap / dmax;
        // backtracking line search on the projected path
        double alpha = 1.0; bool acc = false;
        for (int ls = 0; ls < 6 && !cancelled(cb); ++ls) {
            std::vector<double> ut(static_cast<size_t>(n));
            for (int k = 0; k < n; ++k) ut[size_t(k)] = std::max(0.0, std::min(1.0, u[size_t(k)] + alpha * d[size_t(k)]));
            std::vector<double> st(static_cast<size_t>(n)); for (int k = 0; k < n; ++k) st[size_t(k)] = ut[size_t(k)] - u[size_t(k)];
            const double slope = dot(g, st);
            double ft = 0;
            RecordedRun Rt = evalRun(ut, ft);
            if (Rt.run.ok && ft <= f + 1e-4 * slope) {
                std::vector<double> gt;
                if (!gradOf(Rt, ut, gt)) break;
                std::vector<double> yv(static_cast<size_t>(n)); for (int k = 0; k < n; ++k) yv[size_t(k)] = gt[size_t(k)] - g[size_t(k)];
                if (dot(st, yv) > 1e-12 * std::sqrt(dot(st, st) * dot(yv, yv))) {
                    S.push_back(st); Y.push_back(yv);
                    if (int(S.size()) > std::max(1, opt.memory)) { S.erase(S.begin()); Y.erase(Y.begin()); }
                }
                const double rel = (f - ft) / std::max(f, 1e-300);
                u = ut; f = ft; g = gt;
                hm.bestMisfit = f; hm.best = Rt.values; hm.bestRun = Rt.run;
                std::string vs;
                for (int k = 0; k < n; ++k) vs += " " + paramLabel(params[size_t(act[size_t(k)])], base) + "=" + f2s("%.4g", hm.best[size_t(act[size_t(k)])]);
                say(cb, f2s("  iter %g: misfit %.5g  step %.3g ", it, f, alpha) + vs);
                record(it, f, alpha, g, Rt.run);
                acc = true;
                if (rel < opt.tol) { say(cb, "  converged (relative improvement below tolerance)"); it = opt.maxIter + 1; }
                break;
            }
            // quadratic interpolation of the step, safeguarded to [0.1, 0.5] alpha
            double an = 0.5 * alpha;
            if (Rt.run.ok && std::isfinite(ft) && slope < 0) {
                const double den = 2.0 * (ft - f - slope);
                if (den > 0) an = std::max(0.1 * alpha, std::min(0.5 * alpha, -slope * alpha / den));
            }
            alpha = an;
        }
        if (!acc) {
            if (!S.empty()) { S.clear(); Y.clear(); say(cb, "  line search failed - restarting from steepest descent"); continue; }
            say(cb, "  no further improvement - stopping");
            break;
        }
    }
    hm.ok = true;
    say(cb, f2s("Adjoint L-BFGS done: misfit %.5g -> %.5g in %g runs", hm.initialMisfit, hm.bestMisfit, hm.runs));
    return hm;
}

// ======================================================================== Levenberg-Marquardt
HMResult runLevenbergMarquardt(const Case& base, const std::vector<StudyParam>& params, const LMOptions& opt,
                               const StudyCallbacks& cb) {
    HMResult hm;
    const auto tStart = std::chrono::steady_clock::now();
    auto stamp = [&](HMIteration& h, const RunResult& run) {
        h.runs = hm.runs; h.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count(); h.run = run;
    };
    const auto& obs = base.observations;
    if (obs.empty()) { say(cb, "History matching needs observed data"); return hm; }
    const std::vector<double> sig = autoSigma(obs);
    std::vector<int> act;
    for (int i = 0; i < int(params.size()); ++i) if (params[i].active) act.push_back(i);
    const int n = int(act.size());
    if (n == 0) { say(cb, "No active parameters"); return hm; }
    const int nd = int(obs.size());
    auto valuesOf = [&](const std::vector<double>& u) {
        std::vector<double> v = baseValues(params);
        for (int k = 0; k < n; ++k) v[act[k]] = paramFromUnit(params[act[k]], u[k]);
        return v;
    };
    std::vector<double> u(n);
    for (int k = 0; k < n; ++k) u[k] = paramToUnit(params[act[k]], params[act[k]].base);
    // base run
    Case c0 = applyParams(base, params, valuesOf(u));
    hm.initialRun = runSingle(c0, cb.cancel); hm.runs++;
    archive1(cb, "base", valuesOf(u), hm.initialRun);
    if (!hm.initialRun.ok) { say(cb, "Base run failed: " + hm.initialRun.error); return hm; }
    std::vector<double> r = residuals(hm.initialRun, c0, obs, sig);
    auto sq = [&](const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return s / nd; };
    double f = sq(r);
    hm.initialMisfit = f; hm.bestMisfit = f; hm.best = valuesOf(u); hm.bestRun = hm.initialRun;
    hm.history.push_back({0, f, opt.lambda0, hm.best, 0});
    stamp(hm.history.back(), hm.initialRun);
    say(cb, f2s("LM: %g parameters, %g observations, initial misfit %.5g", n, nd, f));
    double lambda = opt.lambda0;
    for (int it = 1; it <= opt.maxIter && !cancelled(cb); ++it) {
        std::vector<double> J(size_t(nd) * n);
        bool haveJ = false;
        if (opt.gradient == GradientMethod::Analytic) {
            // exact Jacobian of the discrete model: one recorded run + forward sensitivities (no extra simulations)
            RecordedRun R = recordedRun(base, params, act, u, cb.cancel);
            hm.runs++;
            if (R.run.ok) {
                AnalyticOut a = analyticGradient(R, base, params, act, u, sig, true, opt.resStep);
                if (a.ok) {
                    for (int i = 0; i < nd; ++i) for (int k = 0; k < n; ++k) J[size_t(i) * n + k] = a.jac[size_t(i)][size_t(k)];
                    haveJ = true;
                } else say(cb, "  analytic sensitivities failed (" + a.error + ") - finite differences used");
            }
        }
        if (!haveJ) {
            // finite-difference Jacobian (parallel)
            std::vector<Case> cases; std::vector<double> hstep(n);
            std::vector<std::vector<double>> jv;
            for (int k = 0; k < n; ++k) {
                std::vector<double> up = u;
                hstep[k] = (u[k] + opt.fdStep <= 1.0) ? opt.fdStep : -opt.fdStep;
                up[k] += hstep[k];
                jv.push_back(valuesOf(up));
                cases.push_back(applyParams(base, params, jv.back()));
            }
            auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
            archive(cb, "LM", jv, runs);
            hm.runs += n;
            for (int k = 0; k < n; ++k) {
                const auto rk = residuals(runs[k], cases[k], obs, sig);
                for (int i = 0; i < nd; ++i) J[size_t(i) * n + k] = runs[k].ok ? (rk[i] - r[i]) / hstep[k] : 0.0;
            }
        }
        std::vector<double> JtJ(size_t(n) * n, 0.0), g(n, 0.0);
        for (int i = 0; i < nd; ++i)
            for (int a = 0; a < n; ++a) {
                const double ja = J[size_t(i) * n + a];
                g[a] += ja * r[i];
                for (int b = 0; b < n; ++b) JtJ[size_t(a) * n + b] += ja * J[size_t(i) * n + b];
            }
        {   // gradient of the misfit (sum r^2 / nd) at the current iterate: 2 J^T r / nd
            HMIteration& h = hm.history.back();
            h.grad.assign(params.size(), std::numeric_limits<double>::quiet_NaN());
            double gn = 0;
            for (int k = 0; k < n; ++k) { const double gk = 2.0 * g[k] / nd; h.grad[size_t(act[k])] = gk; gn += gk * gk; }
            h.gradNorm = std::sqrt(gn);
        }
        // three damping values tried in parallel
        bool accepted = false;
        for (int attempt = 0; attempt < 4 && !accepted && !cancelled(cb); ++attempt) {
            const double lams[3] = {lambda * 0.1, lambda, lambda * 10};
            std::vector<std::vector<double>> trialU;
            std::vector<Case> tc;
            for (double lam : lams) {
                std::vector<double> A = JtJ, d(n);
                for (int a = 0; a < n; ++a) { A[size_t(a) * n + a] += lam * std::max(JtJ[size_t(a) * n + a], 1e-6) + 1e-12; d[a] = -g[a]; }
                if (!solveDense(A, d, n)) continue;
                std::vector<double> un = u;
                for (int a = 0; a < n; ++a) un[a] = std::max(0.0, std::min(1.0, u[a] + d[a]));
                trialU.push_back(un);
                tc.push_back(applyParams(base, params, valuesOf(un)));
            }
            auto tr = runBatch(tc, cb.workers, cb.cancel, cb.progress);
            { std::vector<std::vector<double>> tv; for (auto& x : trialU) tv.push_back(valuesOf(x)); archive(cb, "LM", tv, tr); }
            hm.runs += int(tc.size());
            int bi = -1; double bf = f; std::vector<double> br;
            for (size_t q = 0; q < tc.size(); ++q) {
                if (!tr[q].ok) continue;
                auto rq = residuals(tr[q], tc[q], obs, sig);
                const double fq = sq(rq);
                if (fq < bf) { bf = fq; bi = int(q); br = rq; }
            }
            if (bi >= 0) {
                const double rel = (f - bf) / std::max(f, 1e-300);
                u = trialU[bi]; r = br; f = bf;
                lambda = std::max(1e-6, lams[bi] * 0.5);
                hm.bestMisfit = f; hm.best = valuesOf(u); hm.bestRun = tr[bi];
                accepted = true;
                std::string vs;
                for (int k = 0; k < n; ++k) vs += " " + paramLabel(params[act[k]], base) + "=" + f2s("%.4g", hm.best[act[k]]);
                say(cb, f2s("  iter %g: misfit %.5g  lambda %.3g ", it, f, lambda) + vs);
                hm.history.push_back({it, f, lambda, hm.best, 0});
                stamp(hm.history.back(), tr[bi]);
                if (rel < opt.tol) { say(cb, "  converged (relative improvement below tolerance)"); it = opt.maxIter + 1; }
            } else lambda *= 100;
        }
        if (!accepted) { say(cb, "  no further improvement - stopping"); break; }
    }
    hm.ok = true;
    say(cb, f2s("LM done: misfit %.5g -> %.5g in %g runs", hm.initialMisfit, hm.bestMisfit, hm.runs));
    return hm;
}

// ======================================================================== ES-MDA
HMResult runESMDA(const Case& base, const std::vector<StudyParam>& params, const ESMDAOptions& opt, const StudyCallbacks& cb) {
    HMResult hm;
    const auto& obs = base.observations;
    if (obs.empty()) { say(cb, "History matching needs observed data"); return hm; }
    const std::vector<double> sig = autoSigma(obs);
    std::vector<int> act;
    for (int i = 0; i < int(params.size()); ++i) if (params[i].active) act.push_back(i);
    const int n = int(act.size()), nd = int(obs.size()), ne = std::max(4, opt.ensemble), na = std::max(1, opt.assimilations);
    if (n == 0) { say(cb, "No active parameters"); return hm; }
    auto logit = [](double u) { u = std::max(1e-4, std::min(1 - 1e-4, u)); return std::log(u / (1 - u)); };
    auto sigm = [](double z) { return 1.0 / (1.0 + std::exp(-z)); };
    auto valuesOf = [&](const std::vector<double>& z) {
        std::vector<double> v = baseValues(params);
        for (int k = 0; k < n; ++k) v[act[k]] = paramFromUnit(params[act[k]], sigm(z[k]));
        return v;
    };
    std::mt19937 rng(opt.seed); std::normal_distribution<double> N(0.0, 1.0);
    std::vector<std::vector<double>> Z(ne, std::vector<double>(n));
    for (int k = 0; k < n; ++k) {
        const double z0 = logit(std::max(0.02, std::min(0.98, paramToUnit(params[act[k]], params[act[k]].base))));
        for (int j = 0; j < ne; ++j) Z[j][k] = z0 + opt.priorStd * N(rng);
    }
    // base run for reference
    Case c0 = applyParams(base, params, baseValues(params));
    hm.initialRun = runSingle(c0, cb.cancel); hm.runs++;
    archive1(cb, "base", baseValues(params), hm.initialRun);
    hm.initialMisfit = misfit(hm.initialRun, c0, obs, sig);
    say(cb, f2s("ES-MDA: ensemble %g, %g assimilations, %g observations", ne, na, nd) + f2s(", base misfit %.5g", hm.initialMisfit));
    const double alpha = double(na);
    std::vector<double> dobs(nd);
    for (int i = 0; i < nd; ++i) dobs[i] = obs[i].value / sig[i];
    std::vector<RunResult> runs;
    std::vector<Case> cases;
    auto runEnsemble = [&]() {
        cases.clear();
        std::vector<std::vector<double>> ev;
        for (int j = 0; j < ne; ++j) { ev.push_back(valuesOf(Z[j])); cases.push_back(applyParams(base, params, ev.back())); }
        runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
        archive(cb, "ES-MDA", ev, runs);
        hm.runs += ne;
    };
    for (int a = 0; a <= na && !cancelled(cb); ++a) {
        runEnsemble();
        if (cancelled(cb)) break;
        std::vector<int> ok;
        std::vector<std::vector<double>> D(ne);
        std::vector<double> mf(ne, kNaN);
        for (int j = 0; j < ne; ++j) {
            if (!runs[j].ok) continue;
            ok.push_back(j);
            D[j].resize(nd);
            for (int i = 0; i < nd; ++i) {
                const double s = simulatedValue(runs[j], cases[j], obs[i].qty, obs[i].well, obs[i].time);
                D[j][i] = std::isfinite(s) ? s / sig[i] : dobs[i] + 30.0;
            }
            mf[j] = misfit(runs[j], cases[j], obs, sig);
        }
        // statistics
        double mean = 0, best = 1e300; int bj = -1;
        for (int j : ok) { mean += mf[j]; if (mf[j] < best) { best = mf[j]; bj = j; } }
        mean = ok.empty() ? kNaN : mean / ok.size();
        double spread = 0;
        for (int k = 0; k < n; ++k) {
            double m = 0, s2 = 0;
            for (int j = 0; j < ne; ++j) m += sigm(Z[j][k]);
            m /= ne;
            for (int j = 0; j < ne; ++j) s2 += (sigm(Z[j][k]) - m) * (sigm(Z[j][k]) - m);
            spread += std::sqrt(s2 / (ne - 1));
        }
        spread /= n;
        HMIteration hi; hi.iter = a; hi.misfit = mean; hi.spread = spread;
        if (bj >= 0) hi.values = valuesOf(Z[bj]);
        hm.history.push_back(hi);
        say(cb, std::string(a == 0 ? "  prior" : a == na ? "  posterior" : "  update ") + (a > 0 && a < na ? f2s("%g", a) : std::string()) +
                f2s(": mean misfit %.5g  best %.5g  spread(u) %.3f", mean, best, spread) + f2s("  failed %g", double(ne - int(ok.size()))));
        if (a == na) {   // final (posterior) ensemble
            hm.ensemble.clear(); hm.ensembleMisfit = mf; hm.ensembleRuns = runs;
            for (int j = 0; j < ne; ++j) hm.ensemble.push_back(valuesOf(Z[j]));
            if (bj >= 0) { hm.best = valuesOf(Z[bj]); hm.bestMisfit = best; hm.bestRun = runs[bj]; hm.ok = true; }
            break;
        }
        const int ns = int(ok.size());
        if (ns < 3) { say(cb, "  too many failed runs - stopping"); break; }
        // ensemble anomalies
        std::vector<double> zm(n, 0.0), dm(nd, 0.0);
        for (int j : ok) { for (int k = 0; k < n; ++k) zm[k] += Z[j][k]; for (int i = 0; i < nd; ++i) dm[i] += D[j][i]; }
        for (double& v : zm) v /= ns;
        for (double& v : dm) v /= ns;
        const double sc = 1.0 / std::sqrt(double(ns - 1));
        std::vector<double> dZ(size_t(n) * ns), dD(size_t(nd) * ns);   // column-major per member
        for (int c = 0; c < ns; ++c) {
            const int j = ok[c];
            for (int k = 0; k < n; ++k) dZ[size_t(c) * n + k] = (Z[j][k] - zm[k]) * sc;
            for (int i = 0; i < nd; ++i) dD[size_t(c) * nd + i] = (D[j][i] - dm[i]) * sc;
        }
        // G = dD^T dD + alpha I  (ns x ns)
        std::vector<double> G(size_t(ns) * ns, 0.0);
        for (int p = 0; p < ns; ++p)
            for (int q = p; q < ns; ++q) {
                double s = 0;
                for (int i = 0; i < nd; ++i) s += dD[size_t(p) * nd + i] * dD[size_t(q) * nd + i];
                G[size_t(p) * ns + q] = G[size_t(q) * ns + p] = s;
            }
        for (int p = 0; p < ns; ++p) G[size_t(p) * ns + p] += alpha;
        for (int c = 0; c < ns; ++c) {
            const int j = ok[c];
            std::vector<double> innov(nd), y(ns, 0.0);
            for (int i = 0; i < nd; ++i) innov[i] = dobs[i] + std::sqrt(alpha) * N(rng) - D[j][i];
            for (int p = 0; p < ns; ++p) { double s = 0; for (int i = 0; i < nd; ++i) s += dD[size_t(p) * nd + i] * innov[i]; y[p] = s; }
            if (!solveDense(G, y, ns)) continue;
            for (int k = 0; k < n; ++k) {
                double s = 0; for (int p = 0; p < ns; ++p) s += dZ[size_t(p) * n + k] * y[p];
                Z[j][k] = std::max(-7.0, std::min(7.0, Z[j][k] + s));
            }
        }
    }
    if (hm.ok) say(cb, f2s("ES-MDA done: best misfit %.5g (base %.5g), %g runs", hm.bestMisfit, hm.initialMisfit, hm.runs));
    return hm;
}


// ======================================================================== multi-objective
std::vector<ObjectiveTerm> objectiveTerms(const std::vector<Observation>& obs, ObjectiveGrouping g) {
    std::vector<ObjectiveTerm> out;
    std::map<std::string, size_t> at;
    auto cls = [](ObsQty q) -> const char* {
        switch (q) {
            case ObsQty::WOPR: case ObsQty::FOPR: return "Oil rate";
            case ObsQty::WWPR: case ObsQty::FWPR: case ObsQty::WWIR: return "Water rate";
            case ObsQty::WGPR: case ObsQty::FGPR: return "Gas rate";
            case ObsQty::WBHP: case ObsQty::FPR: return "Pressure";
            case ObsQty::WWCT: case ObsQty::FWCT: return "Water cut";
            case ObsQty::WGOR: case ObsQty::FGOR: return "GOR";
        }
        return "Other";
    };
    for (size_t i = 0; i < obs.size(); ++i) {
        const Observation& o = obs[i];
        std::string k;
        switch (g) {
            case ObjectiveGrouping::Quantity: k = toString(o.qty); break;
            case ObjectiveGrouping::Well: k = o.well.empty() ? std::string("FIELD") : o.well; break;
            case ObjectiveGrouping::QuantityWell: k = std::string(toString(o.qty)) + (o.well.empty() ? "" : " " + o.well); break;
            case ObjectiveGrouping::PhaseClass: k = cls(o.qty); break;
        }
        auto it = at.find(k);
        if (it == at.end()) { at[k] = out.size(); out.push_back({k, {}}); it = at.find(k); }
        out[it->second].idx.push_back(int(i));
    }
    return out;
}

std::vector<double> misfitTerms(const RunResult& r, const Case& c, const std::vector<Observation>& obs,
                                const std::vector<double>& sigma, const std::vector<ObjectiveTerm>& terms) {
    const auto res = residuals(r, c, obs, sigma);
    std::vector<double> out;
    for (const auto& t : terms) {
        double s = 0; int n = 0;
        for (int i : t.idx) if (i >= 0 && i < int(res.size()) && sigma[size_t(i)] < 1e29) { s += res[size_t(i)] * res[size_t(i)]; ++n; }
        out.push_back(n ? s / n : 0.0);
    }
    return out;
}

std::vector<int> paretoFront(const std::vector<std::vector<double>>& f) {
    std::vector<int> out;
    auto finite = [](const std::vector<double>& v) { for (double x : v) if (!std::isfinite(x)) return false; return !v.empty(); };
    for (size_t i = 0; i < f.size(); ++i) {
        if (!finite(f[i])) continue;
        bool dominated = false;
        for (size_t j = 0; j < f.size() && !dominated; ++j) {
            if (i == j || !finite(f[j]) || f[j].size() != f[i].size()) continue;
            bool le = true, lt = false;
            for (size_t k = 0; k < f[i].size(); ++k) { if (f[j][k] > f[i][k]) { le = false; break; } if (f[j][k] < f[i][k]) lt = true; }
            dominated = le && lt;
        }
        if (!dominated) out.push_back(int(i));
    }
    return out;
}

// ======================================================================== experimental design
const char* toString(DesignKind k) {
    switch (k) {
        case DesignKind::LatinHypercube: return "Latin hypercube";
        case DesignKind::FullFactorial2: return "Full factorial (2-level + centre)";
        case DesignKind::CentralComposite: return "Central composite (face-centred)";
        case DesignKind::PlackettBurman: return "Plackett-Burman (screening)";
        case DesignKind::Random: return "Uniform random";
    }
    return "?";
}

std::vector<std::vector<double>> designPoints(const std::vector<StudyParam>& params, DesignKind kind, int samples, unsigned seed) {
    std::vector<int> act;
    for (int i = 0; i < int(params.size()); ++i) if (params[i].active) act.push_back(i);
    const int n = int(act.size());
    std::vector<std::vector<double>> U;   // unit coordinates of the active parameters
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> R(0.0, 1.0);
    if (n == 0) return {};
    switch (kind) {
        case DesignKind::LatinHypercube: case DesignKind::Random: {
            const int ns = std::max(2, samples);
            U.assign(ns, std::vector<double>(n));
            for (int k = 0; k < n; ++k) {
                std::vector<int> perm(ns); std::iota(perm.begin(), perm.end(), 0);
                std::shuffle(perm.begin(), perm.end(), rng);
                for (int s = 0; s < ns; ++s) U[s][k] = kind == DesignKind::Random ? R(rng) : (perm[s] + R(rng)) / ns;
            }
        } break;
        case DesignKind::FullFactorial2: case DesignKind::CentralComposite: {
            const int nf = std::min(n, 10);
            for (int m = 0; m < (1 << nf); ++m) {
                std::vector<double> u(n, 0.5);
                for (int k = 0; k < nf; ++k) u[k] = (m >> k) & 1 ? 1.0 : 0.0;
                U.push_back(u);
            }
            if (kind == DesignKind::CentralComposite)
                for (int k = 0; k < n; ++k) for (double a : {0.0, 1.0}) { std::vector<double> u(n, 0.5); u[k] = a; U.push_back(u); }
            U.push_back(std::vector<double>(n, 0.5));
        } break;
        case DesignKind::PlackettBurman: {
            static const char* g12 = "++-+++---+-";
            static const char* g20 = "++--++++-+-+----++-";
            static const char* g24 = "+++++-+-++--++--+-+----";
            std::vector<std::vector<int>> M;
            auto cyclic = [&](const char* g) {
                const int L = int(std::strlen(g));
                for (int r = 0; r < L; ++r) { std::vector<int> row(L); for (int c2 = 0; c2 < L; ++c2) row[c2] = g[(c2 - r + L) % L] == '+' ? 1 : -1; M.push_back(row); }
                M.push_back(std::vector<int>(L, -1));
            };
            auto sylvester = [&](int N) {
                std::vector<std::vector<int>> H{{1}};
                while (int(H.size()) < N) {
                    const int m = int(H.size());
                    std::vector<std::vector<int>> H2(2 * m, std::vector<int>(2 * m));
                    for (int a = 0; a < m; ++a) for (int b = 0; b < m; ++b) { H2[a][b] = H2[a][b + m] = H2[a + m][b] = H[a][b]; H2[a + m][b + m] = -H[a][b]; }
                    H = H2;
                }
                for (auto& row : H) M.push_back(std::vector<int>(row.begin() + 1, row.end()));   // drop the constant column
            };
            if (n <= 7) sylvester(8);
            else if (n <= 11) cyclic(g12);
            else if (n <= 15) sylvester(16);
            else if (n <= 19) cyclic(g20);
            else if (n <= 23) cyclic(g24);
            else { int N = 32; while (N - 1 < n) N *= 2; sylvester(N); }
            for (auto& row : M) { std::vector<double> u(n); for (int k = 0; k < n; ++k) u[k] = row[size_t(k)] > 0 ? 1.0 : 0.0; U.push_back(u); }
            U.push_back(std::vector<double>(n, 0.5));
        } break;
    }
    std::vector<std::vector<double>> out;
    for (auto& u : U) {
        std::vector<double> v = baseValues(params);
        for (int k = 0; k < n; ++k) v[act[k]] = paramFromUnit(params[act[k]], u[k]);
        out.push_back(v);
    }
    return out;
}

std::vector<RunResult> runDesign(const Case& base, const std::vector<StudyParam>& params,
                                 const std::vector<std::vector<double>>& values, const char* tag, const StudyCallbacks& cb) {
    std::vector<Case> cases;
    for (auto& v : values) cases.push_back(applyParams(base, params, v));
    say(cb, f2s("%g runs", double(cases.size())) + " (" + tag + ")");
    auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
    archive(cb, tag, values, runs);
    int failed = 0; for (auto& r : runs) if (!r.ok) ++failed;
    if (failed) say(cb, f2s("  %g runs failed", double(failed)));
    return runs;
}

// ======================================================================== response-surface proxy
namespace {
std::vector<double> proxyFeatures(const std::vector<double>& u, bool quad) {
    std::vector<double> x{1.0};
    for (double a : u) x.push_back(a - 0.5);
    if (quad)
        for (size_t k = 0; k < u.size(); ++k)
            for (size_t l = k; l < u.size(); ++l) x.push_back((u[k] - 0.5) * (u[l] - 0.5));
    return x;
}
std::vector<double> unitOf(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values) {
    std::vector<double> u;
    for (int k : m.act) u.push_back(paramToUnit(params[size_t(k)], k < int(values.size()) ? values[size_t(k)] : params[size_t(k)].base));
    return u;
}
// C28 Gaussian process: k(u, v) = sf2 exp(-1/2 sum ((u_k - v_k) / ell_k)^2)
double gpKernel(const ProxyModel& m, const std::vector<double>& u, const std::vector<double>& v) {
    double q = 0;
    for (size_t k = 0; k < u.size() && k < m.gpEll.size(); ++k) { const double d = (u[k] - v[k]) / m.gpEll[k]; q += d * d; }
    return m.gpSf2 * std::exp(-0.5 * q);
}
double gpPredictUnit(const ProxyModel& m, const std::vector<double>& u, double* sd) {
    const size_t ns = m.gpX.size();
    std::vector<double> ks(ns);
    double f = 0;
    for (size_t i = 0; i < ns; ++i) { ks[i] = gpKernel(m, u, m.gpX[i]); f += ks[i] * m.gpAlpha[i]; }
    if (sd) {   // v = L^-1 k*, var = sf2 - v.v
        std::vector<double> v(ns);
        for (size_t i = 0; i < ns; ++i) { double t = ks[i]; for (size_t j = 0; j < i; ++j) t -= m.gpL[i * ns + j] * v[j]; v[i] = t / m.gpL[i * ns + i]; }
        double vv = 0; for (double x : v) vv += x * x;
        *sd = m.gpScale * std::sqrt(std::max(0.0, m.gpSf2 - vv));
    }
    return m.gpMean + m.gpScale * f;
}
double predictUnit(const ProxyModel& m, const std::vector<double>& u) {
    if (m.kind == 1) return gpPredictUnit(m, u, nullptr);
    const auto x = proxyFeatures(u, m.quadratic);
    double s = 0; for (size_t i = 0; i < x.size() && i < m.coef.size(); ++i) s += x[i] * m.coef[i];
    return s;
}
} // namespace

ProxyModel fitProxy(const std::vector<StudyParam>& params, const std::vector<std::vector<double>>& values,
                    const std::vector<double>& y, bool quadratic) {
    ProxyModel m;
    for (int i = 0; i < int(params.size()); ++i) if (params[i].active) m.act.push_back(i);
    const int n = int(m.act.size());
    if (n == 0) return m;
    std::vector<std::vector<double>> X; std::vector<double> Y;
    for (size_t s = 0; s < values.size() && s < y.size(); ++s) {
        if (!std::isfinite(y[s])) continue;
        X.push_back(proxyFeatures(unitOf(m, params, values[s]), quadratic));
        Y.push_back(y[s]);
    }
    const int pq = 1 + n + n * (n + 1) / 2;
    if (quadratic && int(Y.size()) < pq + 2) {   // not enough runs for a quadratic surface
        quadratic = false;
        for (auto& x : X) x.resize(size_t(1 + n));
    }
    m.quadratic = quadratic;
    const int p = int(X.empty() ? 0 : X[0].size()), ns = int(Y.size());
    m.samples = ns;
    if (ns < p + 1) return m;
    std::vector<double> A(size_t(p) * p, 0.0), b(size_t(p), 0.0);
    for (int s = 0; s < ns; ++s)
        for (int a = 0; a < p; ++a) {
            b[size_t(a)] += X[size_t(s)][size_t(a)] * Y[size_t(s)];
            for (int c2 = 0; c2 < p; ++c2) A[size_t(a) * p + c2] += X[size_t(s)][size_t(a)] * X[size_t(s)][size_t(c2)];
        }
    double tr = 0; for (int a = 0; a < p; ++a) tr += A[size_t(a) * p + a];
    for (int a = 1; a < p; ++a) A[size_t(a) * p + a] += 1e-9 * tr / p;   // light ridge (not on the intercept)
    std::vector<double> coef = b;
    if (!solveDense(A, coef, p)) return m;
    m.coef = coef;
    // inverse of the normal matrix for leave-one-out residuals
    std::vector<double> inv(size_t(p) * p);
    for (int c2 = 0; c2 < p; ++c2) {
        std::vector<double> e(size_t(p), 0.0); e[size_t(c2)] = 1;
        if (!solveDense(A, e, p)) return m;
        for (int a = 0; a < p; ++a) inv[size_t(a) * p + c2] = e[size_t(a)];
    }
    double mean = 0; for (double v : Y) mean += v; mean /= ns;
    double sst = 0, sse = 0, press = 0;
    for (int s = 0; s < ns; ++s) {
        double f = 0; for (int a = 0; a < p; ++a) f += X[size_t(s)][size_t(a)] * coef[size_t(a)];
        const double e = Y[size_t(s)] - f;
        double h = 0;
        for (int a = 0; a < p; ++a) { double t = 0; for (int c2 = 0; c2 < p; ++c2) t += inv[size_t(a) * p + c2] * X[size_t(s)][size_t(c2)]; h += X[size_t(s)][size_t(a)] * t; }
        const double el = e / std::max(1e-6, 1 - h);
        sst += (Y[size_t(s)] - mean) * (Y[size_t(s)] - mean); sse += e * e; press += el * el;
    }
    m.r2 = sst > 0 ? 1 - sse / sst : 1.0;
    m.q2 = sst > 0 ? 1 - press / sst : 1.0;
    m.rmse = std::sqrt(sse / ns);
    m.ok = true;
    for (int k = 0; k < n; ++k) {
        std::vector<double> u(size_t(n), 0.5), v = u;
        u[size_t(k)] = 0; v[size_t(k)] = 1;
        double lo = predictUnit(m, u), hi = predictUnit(m, v), mid = predictUnit(m, std::vector<double>(size_t(n), 0.5));
        m.importance.push_back(std::max({std::fabs(hi - lo), std::fabs(hi - mid), std::fabs(lo - mid)}));
    }
    return m;
}

// ---- C28 Gaussian-process proxy
namespace {
// Cholesky of the dense SPD matrix A (row-major, lower factor in place); false if not positive definite
bool cholesky(std::vector<double>& A, size_t n) {
    for (size_t j = 0; j < n; ++j) {
        double d = A[j * n + j];
        for (size_t k = 0; k < j; ++k) d -= A[j * n + k] * A[j * n + k];
        if (!(d > 0)) return false;
        d = std::sqrt(d); A[j * n + j] = d;
        for (size_t i = j + 1; i < n; ++i) {
            double t = A[i * n + j];
            for (size_t k = 0; k < j; ++k) t -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = t / d;
        }
        for (size_t k = j + 1; k < n; ++k) A[j * n + k] = 0;
    }
    return true;
}
// negative log marginal likelihood for the hyper-parameters th = (log ell_1..n, log sf2, log sn2); fills L / alpha
double gpNll(ProxyModel& m, const std::vector<double>& th, const std::vector<double>& ys, bool keep) {
    const size_t n = m.act.size(), ns = m.gpX.size();
    ProxyModel t = m;
    t.gpEll.assign(n, 1.0);
    for (size_t k = 0; k < n; ++k) t.gpEll[k] = std::exp(std::clamp(th[k], std::log(0.02), std::log(20.0)));
    t.gpSf2 = std::exp(std::clamp(th[n], std::log(0.05), std::log(20.0)));
    t.gpSn2 = std::exp(std::clamp(th[n + 1], std::log(1e-8), std::log(1.0)));
    std::vector<double> K(ns * ns);
    for (size_t i = 0; i < ns; ++i)
        for (size_t j = 0; j <= i; ++j) { const double k = gpKernel(t, t.gpX[i], t.gpX[j]); K[i * ns + j] = K[j * ns + i] = k; }
    for (size_t i = 0; i < ns; ++i) K[i * ns + i] += t.gpSn2 + 1e-10 * t.gpSf2;
    if (!cholesky(K, ns)) return 1e300;
    std::vector<double> z(ns), a(ns);
    for (size_t i = 0; i < ns; ++i) { double v = ys[i]; for (size_t j = 0; j < i; ++j) v -= K[i * ns + j] * z[j]; z[i] = v / K[i * ns + i]; }
    for (size_t ii = ns; ii-- > 0;) { double v = z[ii]; for (size_t j = ii + 1; j < ns; ++j) v -= K[j * ns + ii] * a[j]; a[ii] = v / K[ii * ns + ii]; }
    double nll = 0;
    for (size_t i = 0; i < ns; ++i) nll += 0.5 * ys[i] * a[i] + std::log(K[i * ns + i]);
    nll += 0.5 * double(ns) * std::log(2 * 3.141592653589793);
    if (keep) { t.gpL = std::move(K); t.gpAlpha = std::move(a); t.gpLogLik = -nll; m = std::move(t); }
    return nll;
}
// Nelder-Mead minimisation (unconstrained; the bounds are applied by clamping inside the objective)
std::vector<double> nelderMead(const std::function<double(const std::vector<double>&)>& f, std::vector<double> x0, double step, int maxEval) {
    const size_t n = x0.size();
    std::vector<std::vector<double>> P(n + 1, x0); std::vector<double> F(n + 1);
    for (size_t i = 0; i < n; ++i) P[i + 1][i] += step;
    for (size_t i = 0; i <= n; ++i) F[i] = f(P[i]);
    int ev = int(n) + 1;
    while (ev < maxEval) {
        std::vector<size_t> o(n + 1); std::iota(o.begin(), o.end(), size_t(0));
        std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return F[a] < F[b]; });
        std::vector<std::vector<double>> Q; std::vector<double> G; for (size_t i : o) { Q.push_back(P[i]); G.push_back(F[i]); }
        P = Q; F = G;
        if (std::fabs(F[n] - F[0]) < 1e-9 * (1 + std::fabs(F[0]))) break;
        std::vector<double> c(n, 0.0);
        for (size_t i = 0; i < n; ++i) for (size_t k = 0; k < n; ++k) c[k] += P[i][k] / double(n);
        auto lin = [&](double t) { std::vector<double> x(n); for (size_t k = 0; k < n; ++k) x[k] = c[k] + t * (P[n][k] - c[k]); return x; };
        const auto xr = lin(-1.0); const double fr = f(xr); ++ev;
        if (fr < F[0]) { const auto xe = lin(-2.0); const double fe = f(xe); ++ev; if (fe < fr) { P[n] = xe; F[n] = fe; } else { P[n] = xr; F[n] = fr; } }
        else if (fr < F[n - 1]) { P[n] = xr; F[n] = fr; }
        else {
            const auto xc = fr < F[n] ? lin(-0.5) : lin(0.5); const double fc = f(xc); ++ev;
            if (fc < std::min(fr, F[n])) { P[n] = xc; F[n] = fc; }
            else { for (size_t i = 1; i <= n; ++i) { for (size_t k = 0; k < n; ++k) P[i][k] = P[0][k] + 0.5 * (P[i][k] - P[0][k]); F[i] = f(P[i]); ++ev; } }
        }
    }
    size_t b = size_t(std::min_element(F.begin(), F.end()) - F.begin());
    return P[b];
}
}  // namespace

ProxyModel fitProxyGP(const std::vector<StudyParam>& params, const std::vector<std::vector<double>>& values, const std::vector<double>& y) {
    ProxyModel m; m.kind = 1; m.quadratic = false;
    for (int i = 0; i < int(params.size()); ++i) if (params[i].active) m.act.push_back(i);
    const size_t n = m.act.size();
    if (n == 0) return m;
    std::vector<double> Y;
    for (size_t s = 0; s < values.size() && s < y.size(); ++s) {
        if (!std::isfinite(y[s])) continue;
        m.gpX.push_back(unitOf(m, params, values[s])); Y.push_back(y[s]);
    }
    const size_t ns = Y.size();
    m.samples = int(ns);
    if (ns < 3) return m;
    double mean = 0; for (double v : Y) mean += v; mean /= double(ns);
    double var = 0; for (double v : Y) var += (v - mean) * (v - mean); var /= double(ns);
    m.gpMean = mean; m.gpScale = var > 0 ? std::sqrt(var) : 1.0;
    std::vector<double> ys(ns); for (size_t i = 0; i < ns; ++i) ys[i] = (Y[i] - mean) / m.gpScale;
    auto f = [&](const std::vector<double>& th) { return gpNll(m, th, ys, false); };
    std::vector<double> best; double fb = 1e300;
    for (double ell0 : {0.3, 1.0}) {   // two starts: short / long correlation
        std::vector<double> th0(n + 2, std::log(ell0)); th0[n] = 0.0; th0[n + 1] = std::log(1e-4);
        const auto th = nelderMead(f, th0, 0.7, 250 * int(n + 2));
        const double v = f(th);
        if (v < fb) { fb = v; best = th; }
    }
    if (fb >= 1e299) return m;
    gpNll(m, best, ys, true);
    // fit + exact leave-one-out residuals: e_i = alpha_i / [K^-1]_ii
    std::vector<double> Kinv(ns * ns, 0.0);
    for (size_t c2 = 0; c2 < ns; ++c2) {
        std::vector<double> z(ns), a(ns);
        for (size_t i = 0; i < ns; ++i) { double v = i == c2 ? 1.0 : 0.0; for (size_t j = 0; j < i; ++j) v -= m.gpL[i * ns + j] * z[j]; z[i] = v / m.gpL[i * ns + i]; }
        for (size_t ii = ns; ii-- > 0;) { double v = z[ii]; for (size_t j = ii + 1; j < ns; ++j) v -= m.gpL[j * ns + ii] * a[j]; a[ii] = v / m.gpL[ii * ns + ii]; }
        for (size_t i = 0; i < ns; ++i) Kinv[i * ns + c2] = a[i];
    }
    double sst = 0, sse = 0, press = 0;
    for (size_t i = 0; i < ns; ++i) {
        const double fit = ys[i] - (m.gpSn2 + 1e-10 * m.gpSf2) * m.gpAlpha[i];
        const double e = ys[i] - fit, el = m.gpAlpha[i] / Kinv[i * ns + i];
        sst += ys[i] * ys[i]; sse += e * e; press += el * el;
    }
    m.r2 = sst > 0 ? 1 - sse / sst : 1.0;
    m.q2 = sst > 0 ? 1 - press / sst : 1.0;
    m.rmse = m.gpScale * std::sqrt(sse / double(ns));
    m.coef.clear();
    m.ok = true;
    for (size_t k = 0; k < n; ++k) {
        std::vector<double> u(n, 0.5), v = u; u[k] = 0; v[k] = 1;
        const double lo = predictUnit(m, u), hi = predictUnit(m, v), mid = predictUnit(m, std::vector<double>(n, 0.5));
        m.importance.push_back(std::max({std::fabs(hi - lo), std::fabs(hi - mid), std::fabs(lo - mid)}));
    }
    return m;
}

double proxyStd(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values) {
    if (!m.ok || m.kind != 1) return kNaN;
    double sd = 0; gpPredictUnit(m, unitOf(m, params, values), &sd);
    return sd;
}

double proxyPredict(const ProxyModel& m, const std::vector<StudyParam>& params, const std::vector<double>& values) {
    if (!m.ok) return kNaN;
    return predictUnit(m, unitOf(m, params, values));
}

std::vector<double> proxyOptimum(const ProxyModel& m, const std::vector<StudyParam>& params, bool maximise, unsigned seed, double* predicted) {
    std::vector<double> v = baseValues(params);
    if (!m.ok) return v;
    const int n = int(m.act.size());
    const double sg = maximise ? -1 : 1;
    std::mt19937 rng(seed); std::uniform_real_distribution<double> R(0.0, 1.0);
    std::vector<double> best(size_t(n), 0.5); double bf = sg * predictUnit(m, best);
    for (int it = 0; it < 20000; ++it) {
        std::vector<double> u(static_cast<size_t>(n), 0.0); for (auto& a : u) a = R(rng);
        const double f = sg * predictUnit(m, u);
        if (f < bf) { bf = f; best = u; }
    }
    double step = 0.1;
    for (int pass = 0; pass < 60 && step > 1e-4; ++pass) {
        bool imp = false;
        for (int k = 0; k < n; ++k)
            for (double d : {step, -step}) {
                std::vector<double> u = best; u[size_t(k)] = std::max(0.0, std::min(1.0, u[size_t(k)] + d));
                const double f = sg * predictUnit(m, u);
                if (f < bf) { bf = f; best = u; imp = true; }
            }
        if (!imp) step *= 0.5;
    }
    for (int k = 0; k < n; ++k) v[size_t(m.act[size_t(k)])] = paramFromUnit(params[size_t(m.act[size_t(k)])], best[size_t(k)]);
    if (predicted) *predicted = sg * bf;
    return v;
}

} // namespace ressim

// ======================================================================== C11: per-cell field history matching (adjoint)
namespace ressim {

Case applyFieldMultipliers(const Case& base, const std::vector<double>& multK, const std::vector<double>& multPhi) {
    Case c = base;
    GridInput& g = c.grid;
    const size_t n = size_t(g.nx) * g.ny * g.nz;
    if (multK.size() == n && g.permx.size() >= n) {
        if (g.permz.empty()) g.permz = g.permx;
        for (size_t q = 0; q < n; ++q) {
            g.permx[q] *= multK[q];
            if (q < g.permy.size()) g.permy[q] *= multK[q];
            if (q < g.permz.size()) g.permz[q] *= multK[q];
        }
    }
    if (multPhi.size() == n && g.poro.size() >= n)
        for (size_t q = 0; q < n; ++q) g.poro[q] = std::min(1.0, g.poro[q] * multPhi[q]);
    return c;
}

FieldHMResult runAdjointFieldHM(const Case& base, const FieldHMOptions& opt, const StudyCallbacks& cb) {
    FieldHMResult out;
    const auto tStart = std::chrono::steady_clock::now();
    try {
        const auto& obs = base.observations;
        if (obs.empty()) throw std::runtime_error("no observed data");
        if (!opt.perm && !opt.poro) throw std::runtime_error("nothing to match (enable permeability and / or porosity)");
        const std::vector<double> sig = autoSigma(obs);
        const int nd = std::max(1, int(obs.size()));
        Grid G; { Case gc = base; expandLgrArrays(gc); G.build(gc); }
        if (opt.perm && !G.fieldDerivOk) throw std::runtime_error("per-cell permeability matching needs a grid without LGRs / dual porosity");
        const int na = G.numActive(), ncell = base.grid.nx * base.grid.ny * base.grid.nz;
        if (G.numLgrCells() > 0) throw std::runtime_error("per-cell matching does not support local grid refinements yet");
        const int nb = (opt.perm ? 1 : 0) + (opt.poro ? 1 : 0), nvar = na * nb;
        const double xlo = std::log(opt.lo), xhi = std::log(opt.hi);
        // smoothness pairs (grid neighbours)
        std::vector<std::pair<int, int>> pairs;
        for (const Connection& cn : G.conns) if (cn.dir == 'X' || cn.dir == 'Y' || cn.dir == 'Z') pairs.push_back({cn.a, cn.b});
        const double nf = std::max<size_t>(1, pairs.size());
        auto mults = [&](const std::vector<double>& x, std::vector<double>& mK, std::vector<double>& mP) {
            mK.clear(); mP.clear();
            if (opt.perm) mK.assign(size_t(ncell), 1.0);
            if (opt.poro) mP.assign(size_t(ncell), 1.0);
            for (int a = 0; a < na; ++a) {
                const int q = G.globalOf[size_t(a)];
                int o = 0;
                if (opt.perm) mK[size_t(q)] = std::exp(x[size_t(a) * nb + o++]);
                if (opt.poro) mP[size_t(q)] = std::exp(x[size_t(a) * nb + o]);
            }
        };
        struct Eval { bool ok = false; double J = 0, mis = 0, reg = 0; std::vector<double> g, gd; RunResult run; };
        auto evaluate = [&](const std::vector<double>& x, bool wantGrad) {
            Eval E;
            std::vector<double> mK, mP; mults(x, mK, mP);
            const Case c = applyFieldMultipliers(base, mK, mP);
            Simulator S; S.cancel = cb.cancel; S.setAdjointRecording(wantGrad);
            const auto t0 = std::chrono::steady_clock::now();
            try {
                S.setup(c); S.initialize();
                while (!S.finished()) if (!S.runReportStep()) break;
            } catch (const std::exception& e) { E.run.error = e.what(); return E; }
            out.secForward += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            ++out.runs;
            E.run.summary = S.summary(); E.run.stats = S.stats(); E.run.ok = S.finished();
            if (!E.run.ok) return E;
            std::vector<Simulator::AdjObs> ao;
            for (int i = 0; i < int(obs.size()); ++i) {
                if (!(sig[size_t(i)] < 1e29)) continue;
                Simulator::AdjObs o; o.qty = obs[size_t(i)].qty; o.time = obs[size_t(i)].time; o.value = obs[size_t(i)].value; o.sigma = sig[size_t(i)];
                if (obsIsWell(o.qty)) { o.well = wellIndex(c, obs[size_t(i)].well); if (o.well < 0) continue; }
                ao.push_back(o);
            }
            double reg = 0;
            for (auto& pr : pairs)
                for (int o = 0; o < nb; ++o) { const double d = x[size_t(pr.first) * nb + o] - x[size_t(pr.second) * nb + o]; reg += d * d; }
            reg *= opt.regWeight / nf;
            E.reg = reg;
            if (!wantGrad) {   // misfit only
                E.mis = misfit(E.run, c, obs, sig);
                E.J = E.mis + reg; E.ok = std::isfinite(E.J); return E;
            }
            Simulator::AdjFieldReq fr; fr.perm = opt.perm; fr.poro = opt.poro;
            const auto t1 = std::chrono::steady_clock::now();
            const Simulator::AdjResult A = S.adjointGradient(ao, {}, false, &fr);
            out.secAdjoint += std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            if (!A.ok) throw std::runtime_error("adjoint: " + A.error);
            E.mis = A.F / nd; E.J = E.mis + reg;
            E.g.assign(size_t(nvar), 0.0);
            for (int a = 0; a < na; ++a) {
                int o = 0;
                if (opt.perm) E.g[size_t(a) * nb + o++] = (A.gLogKx[size_t(a)] + A.gLogKy[size_t(a)] + A.gLogKz[size_t(a)]) / nd;
                if (opt.poro) E.g[size_t(a) * nb + o] = A.gLogPoro[size_t(a)] / nd;
            }
            E.gd = E.g;   // data-misfit part (sensitivity maps)
            for (auto& pr : pairs)
                for (int o = 0; o < nb; ++o) {
                    const double d = 2 * opt.regWeight / nf * (x[size_t(pr.first) * nb + o] - x[size_t(pr.second) * nb + o]);
                    E.g[size_t(pr.first) * nb + o] += d; E.g[size_t(pr.second) * nb + o] -= d;
                }
            E.ok = std::isfinite(E.J);
            return E;
        };
        std::vector<double> x(static_cast<size_t>(nvar), 0.0);
        Eval cur = evaluate(x, true);
        if (!cur.ok) throw std::runtime_error("the base run failed: " + cur.run.error);
        out.initialMisfit = cur.mis; out.history.push_back(cur.mis); out.initialRun = cur.run;
        out.nx = base.grid.nx; out.ny = base.grid.ny; out.nz = base.grid.nz;
        // C13: per-iteration record (statistics, ln-multiplier snapshots as float, the run)
        auto toGlobal = [&](const std::vector<double>& v, int o, std::vector<float>& dst, float fill) {
            dst.assign(size_t(ncell), fill);
            for (int a = 0; a < na; ++a) dst[size_t(G.globalOf[size_t(a)])] = float(v[size_t(a) * nb + o]);
        };
        auto record = [&](int it, const Eval& E, const std::vector<double>& xv, double step) {
            FieldHMResult::Iter r; r.iter = it; r.runs = out.runs; r.misfit = E.mis; r.objective = E.J; r.reg = E.reg; r.step = step;
            r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
            double gn = 0; for (double v : E.g) gn += v * v; r.gradNorm = std::sqrt(gn);
            int o = 0;
            auto stats = [&](int oo, double& mn, double& mx, double& me) {
                mn = 1e300; mx = -1e300; me = 0;
                for (int a = 0; a < na; ++a) { const double m = std::exp(xv[size_t(a) * nb + oo]); mn = std::min(mn, m); mx = std::max(mx, m); me += m; }
                me /= std::max(1, na);
            };
            if (opt.perm) { stats(o, r.kMin, r.kMax, r.kMean); toGlobal(xv, o, r.lnK, 0.0f); ++o; }
            if (opt.poro) { stats(o, r.phiMin, r.phiMax, r.phiMean); toGlobal(xv, o, r.lnPhi, 0.0f); }
            r.run = E.run;
            out.iters.push_back(std::move(r));
        };
        auto gradMaps = [&](const Eval& E, std::vector<float>& gk, std::vector<float>& gp) {
            const float nan = std::numeric_limits<float>::quiet_NaN();
            int o = 0;
            if (opt.perm) { toGlobal(E.gd, o, gk, nan); ++o; } else gk.clear();
            if (opt.poro) toGlobal(E.gd, o, gp, nan); else gp.clear();
        };
        record(0, cur, x, 0.0);
        gradMaps(cur, out.gradK0, out.gradPhi0);
        say(cb, f2s("Field HM (adjoint): %g cells x %g fields = %g unknowns, initial misfit %.5g", na, nb, cur.mis));
        say(cb, f2s("  observations %g", double(obs.size())));
        std::vector<std::vector<double>> Sm, Ym;
        auto clampv = [&](double v) { return std::min(xhi, std::max(xlo, v)); };
        for (int it = 0; it < opt.maxIter; ++it) {
            if (cancelled(cb)) break;
            // two-loop recursion
            std::vector<double> q = cur.g;
            const int m = int(Sm.size());
            std::vector<double> al(static_cast<size_t>(m)), rho(static_cast<size_t>(m));
            auto dot = [](const std::vector<double>& a, const std::vector<double>& b) { double s = 0; for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i]; return s; };
            for (int i = m - 1; i >= 0; --i) { rho[size_t(i)] = 1.0 / dot(Ym[size_t(i)], Sm[size_t(i)]); al[size_t(i)] = rho[size_t(i)] * dot(Sm[size_t(i)], q); for (int k = 0; k < nvar; ++k) q[size_t(k)] -= al[size_t(i)] * Ym[size_t(i)][size_t(k)]; }
            double gam = m > 0 ? dot(Sm[size_t(m) - 1], Ym[size_t(m) - 1]) / dot(Ym[size_t(m) - 1], Ym[size_t(m) - 1]) : 0.0;
            if (m == 0) { double gmx = 0; for (double v : cur.g) gmx = std::max(gmx, std::fabs(v)); gam = gmx > 0 ? 0.25 * opt.maxLogStep / gmx : 1.0; }
            for (double& v : q) v *= gam;
            for (int i = 0; i < m; ++i) { const double be = rho[size_t(i)] * dot(Ym[size_t(i)], q); for (int k = 0; k < nvar; ++k) q[size_t(k)] += Sm[size_t(i)][size_t(k)] * (al[size_t(i)] - be); }
            // direction -q, limited step, projection on the bounds; fall back to steepest descent if not a descent direction
            std::vector<double> d(static_cast<size_t>(nvar));
            for (int k = 0; k < nvar; ++k) d[size_t(k)] = -q[size_t(k)];
            if (dot(d, cur.g) >= 0) { Sm.clear(); Ym.clear(); double gmx = 0; for (double v : cur.g) gmx = std::max(gmx, std::fabs(v)); for (int k = 0; k < nvar; ++k) d[size_t(k)] = -cur.g[size_t(k)] * (gmx > 0 ? 0.25 * opt.maxLogStep / gmx : 0.0); }
            double dmx = 0; for (double v : d) dmx = std::max(dmx, std::fabs(v));
            double alpha = dmx > opt.maxLogStep ? opt.maxLogStep / dmx : 1.0;
            bool accepted = false;
            double alphaAcc = 0;
            std::vector<double> xn(static_cast<size_t>(nvar));
            Eval nx;
            for (int ls = 0; ls < 8 && !accepted; ++ls, alpha *= 0.5) {
                for (int k = 0; k < nvar; ++k) xn[size_t(k)] = clampv(x[size_t(k)] + alpha * d[size_t(k)]);
                double dec = 0; for (int k = 0; k < nvar; ++k) dec += cur.g[size_t(k)] * (xn[size_t(k)] - x[size_t(k)]);
                if (dec >= 0) continue;
                nx = evaluate(xn, true);
                if (nx.ok && nx.J <= cur.J + 1e-4 * dec) { accepted = true; alphaAcc = alpha * dmx; }
            }
            if (!accepted) { say(cb, "  line search failed - stopping"); break; }
            std::vector<double> sv(static_cast<size_t>(nvar)), yv(static_cast<size_t>(nvar));
            for (int k = 0; k < nvar; ++k) { sv[size_t(k)] = xn[size_t(k)] - x[size_t(k)]; yv[size_t(k)] = nx.g[size_t(k)] - cur.g[size_t(k)]; }
            if (dot(sv, yv) > 1e-12 * std::sqrt(dot(sv, sv) * dot(yv, yv))) { Sm.push_back(sv); Ym.push_back(yv); if (int(Sm.size()) > opt.memory) { Sm.erase(Sm.begin()); Ym.erase(Ym.begin()); } }
            const double rel = (cur.J - nx.J) / std::max(std::fabs(cur.J), 1e-300);
            x = xn; cur = std::move(nx); ++out.iterations;
            out.history.push_back(cur.mis);
            record(out.iterations, cur, x, alphaAcc);
            say(cb, f2s("  iteration %g: misfit %.5g (objective %.5g)", it + 1, cur.mis, cur.J));
            if (rel < opt.tol) break;
        }
        gradMaps(cur, out.gradK, out.gradPhi);
        mults(x, out.multK, out.multPhi);
        out.best = applyFieldMultipliers(base, out.multK, out.multPhi);
        out.bestMisfit = cur.mis; out.bestObjective = cur.J; out.bestRun = cur.run;
        out.ok = true;
        say(cb, f2s("Field HM done: misfit %.5g -> %.5g in %g runs", out.initialMisfit, out.bestMisfit, out.runs));
    } catch (const std::exception& e) { out.ok = false; out.error = e.what(); }
    return out;
}


// ======================================================================== C27 inversion statistics (iTOUGH2-style)
namespace {
// symmetric eigen-decomposition (cyclic Jacobi): A = V diag(w) V^T, columns of V in vec[k][i] = V(i, k)
void symEigen(std::vector<double> A, int n, std::vector<double>& w, std::vector<std::vector<double>>& vec) {
    std::vector<double> V(size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) V[size_t(i) * n + i] = 1;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) off += A[size_t(i) * n + j] * A[size_t(i) * n + j];
        if (off < 1e-30) break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                const double apq = A[size_t(p) * n + q];
                if (std::fabs(apq) < 1e-300) continue;
                const double app = A[size_t(p) * n + p], aqq = A[size_t(q) * n + q];
                const double th = 0.5 * (aqq - app) / apq;
                const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), sn = t * c;
                for (int k = 0; k < n; ++k) {
                    const double akp = A[size_t(k) * n + p], akq = A[size_t(k) * n + q];
                    A[size_t(k) * n + p] = c * akp - sn * akq; A[size_t(k) * n + q] = sn * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const double apk = A[size_t(p) * n + k], aqk = A[size_t(q) * n + k];
                    A[size_t(p) * n + k] = c * apk - sn * aqk; A[size_t(q) * n + k] = sn * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    const double vkp = V[size_t(k) * n + p], vkq = V[size_t(k) * n + q];
                    V[size_t(k) * n + p] = c * vkp - sn * vkq; V[size_t(k) * n + q] = sn * vkp + c * vkq;
                }
            }
    }
    w.assign(size_t(n), 0.0); vec.assign(size_t(n), std::vector<double>(size_t(n), 0.0));
    for (int k = 0; k < n; ++k) { w[size_t(k)] = A[size_t(k) * n + k]; for (int i = 0; i < n; ++i) vec[size_t(k)][size_t(i)] = V[size_t(i) * n + k]; }
}
}  // namespace

HMStatistics hmStatistics(const Case& base, const std::vector<StudyParam>& params, const std::vector<double>& values,
                          const HMStatsOptions& opt, const StudyCallbacks& cb) {
    HMStatistics S;
    const auto& obs = base.observations;
    if (obs.empty()) { S.error = "inversion statistics need observed data"; return S; }
    if (values.size() != params.size()) { S.error = "one value per parameter expected"; return S; }
    for (int i = 0; i < int(params.size()); ++i) if (params[size_t(i)].active) S.param.push_back(i);
    const int n = int(S.param.size()), m = int(obs.size());
    if (n == 0) { S.error = "no active parameters"; return S; }
    if (m <= n) { S.error = f2s("%g observations for %g parameters: at least n + 1 observations are needed", m, n); return S; }
    S.nObs = m; S.nPar = n;
    const std::vector<double> sig = autoSigma(obs);
    // search coordinate: log10 for log parameters
    auto toX = [&](const StudyParam& P, double v) { return P.log ? std::log10(std::max(v, 1e-300)) : v; };
    auto fromX = [&](const StudyParam& P, double x) { return P.log ? std::pow(10.0, x) : x; };
    std::vector<double> x(static_cast<size_t>(n)), h(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) {
        const StudyParam& P = params[size_t(S.param[size_t(k)])];
        S.name.push_back(paramLabel(P, base)); S.value.push_back(values[size_t(S.param[size_t(k)])]); S.logSpace.push_back(P.log);
        x[size_t(k)] = toX(P, S.value.back());
        const double range = std::fabs(toX(P, P.high) - toX(P, P.low));
        h[size_t(k)] = std::max(1e-6, opt.fdStep * (range > 0 ? range : std::max(1.0, std::fabs(x[size_t(k)]))));
        if (toX(P, P.high) < x[size_t(k)] + h[size_t(k)] && toX(P, P.low) <= x[size_t(k)] - h[size_t(k)]) h[size_t(k)] = -h[size_t(k)];   // stay inside the box
    }
    std::vector<Case> cases; std::vector<std::vector<double>> vals;
    vals.push_back(values); cases.push_back(applyParams(base, params, values));
    for (int k = 0; k < n; ++k) {
        std::vector<double> v = values;
        const StudyParam& P = params[size_t(S.param[size_t(k)])];
        v[size_t(S.param[size_t(k)])] = fromX(P, x[size_t(k)] + h[size_t(k)]);
        vals.push_back(v); cases.push_back(applyParams(base, params, v));
    }
    say(cb, f2s("Inversion statistics: %g parameters, %g observations, %g runs", n, m, n + 1));
    const auto runs = runBatch(cases, cb.workers, cb.cancel, cb.progress);
    archive(cb, "stats", vals, runs);
    S.runs = int(runs.size());
    if (cancelled(cb)) { S.error = "cancelled"; return S; }
    if (!runs[0].ok) { S.error = "run at the analysis point failed: " + runs[0].error; return S; }
    const std::vector<double> r = residuals(runs[0], cases[0], obs, sig);
    std::vector<double> J(size_t(m) * n, 0.0);
    for (int k = 0; k < n; ++k) {
        if (!runs[size_t(k) + 1].ok) { S.error = "perturbed run of " + S.name[size_t(k)] + " failed: " + runs[size_t(k) + 1].error; return S; }
        const auto rk = residuals(runs[size_t(k) + 1], cases[size_t(k) + 1], obs, sig);
        for (int i = 0; i < m; ++i) J[size_t(i) * n + k] = (rk[size_t(i)] - r[size_t(i)]) / h[size_t(k)];
    }
    S.residual = r;
    double ssr = 0; for (double v : r) ssr += v * v;
    S.objective = ssr; S.s02 = ssr / double(m - n);
    // Fisher matrix J^T J, eigen-analysis, pseudo-inverse (directions with no sensitivity are reported, not inverted)
    std::vector<double> F(size_t(n) * n, 0.0);
    for (int i = 0; i < m; ++i) for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) F[size_t(a) * n + b] += J[size_t(i) * n + a] * J[size_t(i) * n + b];
    std::vector<double> w; std::vector<std::vector<double>> V;
    symEigen(F, n, w, V);
    double wmax = 0, wmin = 1e300; for (double e : w) { wmax = std::max(wmax, e); wmin = std::min(wmin, e); }
    S.condition = wmin > 0 ? wmax / wmin : std::numeric_limits<double>::infinity();
    std::vector<double> Finv(size_t(n) * n, 0.0);
    for (int k = 0; k < n; ++k) {
        if (w[size_t(k)] <= 1e-12 * std::max(wmax, 1e-300)) continue;
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) Finv[size_t(a) * n + b] += V[size_t(k)][size_t(a)] * V[size_t(k)][size_t(b)] / w[size_t(k)];
    }
    S.cov.assign(size_t(n), std::vector<double>(size_t(n), 0.0)); S.corr = S.cov;
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) S.cov[size_t(a)][size_t(b)] = S.s02 * Finv[size_t(a) * n + b];
    for (int a = 0; a < n; ++a) {
        const double sd = std::sqrt(std::max(0.0, S.cov[size_t(a)][size_t(a)]));
        S.stdDev.push_back(w.empty() || sd == 0 ? std::numeric_limits<double>::infinity() : sd);
        const StudyParam& P = params[size_t(S.param[size_t(a)])];
        const double z = 1.96 * S.stdDev.back();
        S.ci95lo.push_back(std::isfinite(z) ? fromX(P, x[size_t(a)] - z) : (P.log ? 0.0 : -std::numeric_limits<double>::infinity()));
        S.ci95hi.push_back(std::isfinite(z) ? fromX(P, x[size_t(a)] + z) : std::numeric_limits<double>::infinity());
    }
    for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) {
        const double d = std::sqrt(std::max(0.0, S.cov[size_t(a)][size_t(a)] * S.cov[size_t(b)][size_t(b)]));
        S.corr[size_t(a)][size_t(b)] = d > 0 ? S.cov[size_t(a)][size_t(b)] / d : (a == b ? 1.0 : 0.0);
    }
    // covariance eigenvalues = s0^2 / Fisher eigenvalues; worst-determined direction first
    std::vector<int> ord(static_cast<size_t>(n)); std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return w[size_t(a)] < w[size_t(b)]; });
    for (int k : ord) { S.eigVal.push_back(w[size_t(k)] > 0 ? S.s02 / w[size_t(k)] : std::numeric_limits<double>::infinity()); S.eigVec.push_back(V[size_t(k)]); }
    // composite scaled sensitivity (derivative with respect to ln p for log parameters, p * d/dp for linear ones)
    for (int k = 0; k < n; ++k) {
        const double scale = S.logSpace[size_t(k)] ? 1.0 / std::log(10.0) : std::fabs(S.value[size_t(k)]);
        double q = 0; for (int i = 0; i < m; ++i) { const double d = J[size_t(i) * n + k] * scale; q += d * d; }
        S.css.push_back(std::sqrt(q / m));
    }
    // leverage, studentized residuals, prediction standard deviation
    const double s0 = std::sqrt(S.s02);
    for (int i = 0; i < m; ++i) {
        double hii = 0;
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) hii += J[size_t(i) * n + a] * Finv[size_t(a) * n + b] * J[size_t(i) * n + b];
        hii = std::clamp(hii, 0.0, 0.999999);
        S.leverage.push_back(hii);
        const double st = s0 > 0 ? r[size_t(i)] / (s0 * std::sqrt(1.0 - hii)) : 0.0;
        S.studentized.push_back(st);
        if (std::fabs(st) > 1.96) S.outliers.push_back(i);
        S.predStd.push_back(sig[size_t(i)] * s0 * std::sqrt(hii));
    }
    // residual moments
    double mean = 0; for (double v : r) mean += v; mean /= m;
    double m2 = 0, m3 = 0, m4 = 0;
    for (double v : r) { const double d = v - mean; m2 += d * d; m3 += d * d * d; m4 += d * d * d * d; }
    m2 /= m; m3 /= m; m4 /= m;
    S.resMean = mean; S.resStd = std::sqrt(m2);
    S.resSkew = m2 > 0 ? m3 / std::pow(m2, 1.5) : 0; S.resKurt = m2 > 0 ? m4 / (m2 * m2) - 3.0 : 0;
    // objective share per data set
    std::map<std::string, double> share;
    for (int i = 0; i < m; ++i) {
        std::string key = toString(obs[size_t(i)].qty); if (!obs[size_t(i)].well.empty()) key += " " + obs[size_t(i)].well;
        share[key] += r[size_t(i)] * r[size_t(i)];
    }
    for (const auto& e : share) S.contribution.push_back({e.first, ssr > 0 ? 100.0 * e.second / ssr : 0.0});
    std::sort(S.contribution.begin(), S.contribution.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    // ---- C28 FOSM prediction uncertainty and data worth (Finsterle, iTOUGH2: sigma_z^2 = g^T C_pp g; worth of a data set =
    // change of sigma_z when its rows are removed from / added to the Fisher matrix). Directions without sensitivity make
    // a prediction that depends on them unbounded (reported as infinity).
    auto quadStd = [&](const std::vector<double>& Fm, const std::vector<double>& g) {
        std::vector<double> we; std::vector<std::vector<double>> Ve;
        symEigen(Fm, n, we, Ve);
        double wm = 0; for (double e : we) wm = std::max(wm, e);
        double gn = 0; for (double v : g) gn += v * v; gn = std::sqrt(gn);
        double q = 0;
        for (int k = 0; k < n; ++k) {
            double pr = 0; for (int a = 0; a < n; ++a) pr += Ve[size_t(k)][size_t(a)] * g[size_t(a)];
            if (we[size_t(k)] <= 1e-12 * std::max(wm, 1e-300)) { if (std::fabs(pr) > 1e-9 * std::max(gn, 1e-300)) return std::numeric_limits<double>::infinity(); continue; }
            q += pr * pr / we[size_t(k)];
        }
        return std::sqrt(std::max(0.0, S.s02 * q));
    };
    std::vector<std::vector<double>> G;   // [prediction][param]
    for (Response rp : opt.predict) {
        HMStatistics::Prediction P; P.name = toString(rp);
        P.value = responseValue(runs[0], cases[0], rp, opt.eco, obs, sig);
        std::vector<double> g(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k) g[size_t(k)] = (responseValue(runs[size_t(k) + 1], cases[size_t(k) + 1], rp, opt.eco, obs, sig) - P.value) / h[size_t(k)];
        P.stdDev = quadStd(F, g); P.lo = P.value - 1.96 * P.stdDev; P.hi = P.value + 1.96 * P.stdDev;
        S.pred.push_back(P); G.push_back(g);
    }
    if (opt.dataWorth && !G.empty()) {
        auto addRows = [&](std::vector<double>& Fm, const std::vector<double>& rows, int nr, double sign) {
            for (int i = 0; i < nr; ++i) for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b)
                Fm[size_t(a) * n + b] += sign * rows[size_t(i) * n + a] * rows[size_t(i) * n + b];
        };
        auto pct = [&](const std::vector<double>& Fm, HMStatistics::DataWorth& D) {
            for (size_t j = 0; j < G.size(); ++j) {
                const double s1 = quadStd(Fm, G[j]), s0p = S.pred[j].stdDev;
                D.dStdPct.push_back(!std::isfinite(s1) ? std::numeric_limits<double>::infinity() : s0p > 0 && std::isfinite(s0p) ? 100.0 * (s1 / s0p - 1.0) : 0.0);
            }
        };
        // existing data sets: remove their rows
        std::map<std::string, std::vector<int>> sets;
        for (int i = 0; i < m; ++i) { std::string key = toString(obs[size_t(i)].qty); if (!obs[size_t(i)].well.empty()) key += " " + obs[size_t(i)].well; sets[key].push_back(i); }
        for (const auto& e : sets) {
            std::vector<double> rows; for (int i : e.second) for (int a = 0; a < n; ++a) rows.push_back(J[size_t(i) * n + a]);
            std::vector<double> Fm = F; addRows(Fm, rows, int(e.second.size()), -1.0);
            HMStatistics::DataWorth D; D.set = e.first; D.existing = true; D.points = int(e.second.size()); pct(Fm, D); S.worth.push_back(D);
        }
        // candidate data sets: vectors not observed yet, at the observation times (sigma by the automatic rule)
        std::vector<double> times; for (const auto& o : obs) times.push_back(o.time);
        std::sort(times.begin(), times.end()); times.erase(std::unique(times.begin(), times.end(), [](double a, double b) { return std::fabs(a - b) < 1e-6; }), times.end());
        std::vector<std::pair<ObsQty, std::string>> cand;
        auto observed = [&](ObsQty q, const std::string& w) { for (const auto& o : obs) if (o.qty == q && o.well == w) return true; return false; };
        for (const auto& w : base.wells) {
            const bool prod = w.type == WellType::Producer;
            const std::vector<ObsQty> qs = prod ? std::vector<ObsQty>{ObsQty::WBHP, ObsQty::WWCT, ObsQty::WGOR, ObsQty::WOPR} : std::vector<ObsQty>{ObsQty::WBHP, ObsQty::WWIR};
            for (ObsQty q : qs) if (!observed(q, w.name)) cand.push_back({q, w.name});
        }
        if (!observed(ObsQty::FPR, "")) cand.push_back({ObsQty::FPR, ""});
        for (const auto& [q, wn] : cand) {
            std::vector<Observation> pseudo;
            for (double t : times) { Observation o; o.qty = q; o.well = wn; o.time = t; o.value = simulatedValue(runs[0], cases[0], q, wn, t); if (std::isfinite(o.value)) pseudo.push_back(o); }
            if (pseudo.empty()) continue;
            const std::vector<double> sg = autoSigma(pseudo);
            std::vector<double> rows; int nr = 0; double tot = 0;
            for (size_t i = 0; i < pseudo.size(); ++i) {
                std::vector<double> row(static_cast<size_t>(n));
                bool fin = true;
                for (int k = 0; k < n; ++k) {
                    const double yk = simulatedValue(runs[size_t(k) + 1], cases[size_t(k) + 1], q, wn, pseudo[i].time);
                    row[size_t(k)] = (yk - pseudo[i].value) / h[size_t(k)] / sg[i];
                    if (!std::isfinite(row[size_t(k)])) fin = false;
                    tot += std::fabs(row[size_t(k)]);
                }
                if (fin) { rows.insert(rows.end(), row.begin(), row.end()); ++nr; }
            }
            if (nr == 0 || tot == 0) continue;   // no sensitivity at all (e.g. GOR of a dead-oil case): no worth
            std::vector<double> Fm = F; addRows(Fm, rows, nr, 1.0);
            HMStatistics::DataWorth D; D.set = std::string(toString(q)) + (wn.empty() ? "" : " " + wn); D.existing = false; D.points = nr; pct(Fm, D); S.worth.push_back(D);
        }
        // most valuable first (by the first prediction): existing sets by increase, candidates by reduction
        std::stable_sort(S.worth.begin(), S.worth.end(), [](const HMStatistics::DataWorth& a, const HMStatistics::DataWorth& b) {
            if (a.existing != b.existing) return a.existing;
            const double va = a.dStdPct.empty() ? 0 : a.dStdPct[0], vb = b.dStdPct.empty() ? 0 : b.dStdPct[0];
            return a.existing ? va > vb : va < vb;
        });
    }
    S.ok = true;
    say(cb, f2s("  objective %.5g, s0^2 %.4g, condition number %.3g", S.objective, S.s02, S.condition));
    for (const auto& P : S.pred) say(cb, "  prediction " + P.name + f2s(" = %.6g +- %.4g (1 sd, FOSM)", P.value, P.stdDev));
    return S;
}

} // namespace ressim
