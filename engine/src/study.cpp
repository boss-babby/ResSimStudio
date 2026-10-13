// Multi-run studies: batch runner, sensitivity analysis and history matching.
#include "ressim/study.h"
#include "ressim/ecl.h"
#include <thread>
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
    }
    return "?";
}
bool paramIsMultiplier(ParamKind k) {
    return !(k == ParamKind::OWC || k == ParamKind::DatumPressure || k == ParamKind::WellSkin || k == ParamKind::FaultMult ||
             k == ParamKind::KrwExp || k == ParamKind::KroExp || k == ParamKind::GOC || k == ParamKind::KrgExp || k == ParamKind::KrogExp);
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
    GridInput& g = c.grid;
    bool anyPerm = false;
    for (auto& p : params) if (p.kind == ParamKind::PermXY || p.kind == ParamKind::PermZ) anyPerm = true;
    if (anyPerm && g.permz.empty()) g.permz = g.permx;            // keep kz independent of kx multipliers
    for (size_t i = 0; i < params.size(); ++i) {
        const StudyParam& p = params[i];
        const double v = (p.active && i < values.size()) ? values[i] : p.base;
        if (paramIsMultiplier(p.kind) && v == 1.0) continue;
        switch (p.kind) {
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
RunResult runSingle(const Case& c, std::atomic<bool>* cancel) {
    RunResult r;
    try {
        Simulator s;
        s.cancel = cancel;
        s.setup(c);
        s.initialize();
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

std::vector<RunResult> runBatch(const std::vector<Case>& cases, int workers, std::atomic<bool>* cancel,
                                const std::function<void(int, int)>& progress) {
    const int n = int(cases.size());
    std::vector<RunResult> out(n);
    if (n == 0) return out;
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

// ======================================================================== Levenberg-Marquardt
HMResult runLevenbergMarquardt(const Case& base, const std::vector<StudyParam>& params, const LMOptions& opt,
                               const StudyCallbacks& cb) {
    HMResult hm;
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
    say(cb, f2s("LM: %g parameters, %g observations, initial misfit %.5g", n, nd, f));
    double lambda = opt.lambda0;
    for (int it = 1; it <= opt.maxIter && !cancelled(cb); ++it) {
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
        std::vector<double> J(size_t(nd) * n);
        for (int k = 0; k < n; ++k) {
            const auto rk = residuals(runs[k], cases[k], obs, sig);
            for (int i = 0; i < nd; ++i) J[size_t(i) * n + k] = runs[k].ok ? (rk[i] - r[i]) / hstep[k] : 0.0;
        }
        std::vector<double> JtJ(size_t(n) * n, 0.0), g(n, 0.0);
        for (int i = 0; i < nd; ++i)
            for (int a = 0; a < n; ++a) {
                const double ja = J[size_t(i) * n + a];
                g[a] += ja * r[i];
                for (int b = 0; b < n; ++b) JtJ[size_t(a) * n + b] += ja * J[size_t(i) * n + b];
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
double predictUnit(const ProxyModel& m, const std::vector<double>& u) {
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
