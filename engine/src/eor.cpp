// CO2 / gas EOR laboratory tools: swelling test, P-X diagram / FCMP, mixing-cell MMP, MMP correlations, WAG generator
#include "ressim/eor.h"
#include "ressim/case.h"
#include "ressim/compsim.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace ressim {
namespace eor {

namespace {
std::vector<double> normalized(std::vector<double> z, size_t nc) {
    z.resize(nc, 0.0);
    double s = 0; for (double& v : z) { v = std::max(0.0, v); s += v; }
    if (s > 0) for (double& v : z) v /= s;
    return z;
}
std::vector<double> mix(const std::vector<double>& a, const std::vector<double>& b, double wa) {
    std::vector<double> z(a.size());
    for (size_t i = 0; i < a.size(); ++i) z[i] = wa * a[i] + (1.0 - wa) * b[i];
    return z;
}
// saturation pressure of z at T: bubble or dew point according to the mixture label, the other one as a fallback
double psatOf(const eos::Model& m, double T, const std::vector<double>& z, double guess, bool& dew) {
    bool vap = false;
    try { vap = eos::isVapourLike(m, z, T); } catch (...) {}
    for (int pass = 0; pass < 2; ++pass) {
        const bool bubble = pass == 0 ? !vap : vap;
        double p = -1;
        try { p = eos::saturationPressure(m, T, z, bubble, guess, 1e5, 1.5e8); } catch (...) { p = -1; }
        if (!(p > 0) && guess > 0) { try { p = eos::saturationPressure(m, T, z, bubble, 0, 1e5, 1.5e8); } catch (...) { p = -1; } }
        if (p > 1.0e5) { dew = !bubble; return p; }
    }
    // near-critical mixtures: the saturation-pressure Newton can fall into the trivial solution - locate the upper
    // boundary of the two-phase region directly (stability-tested flash, scan down then bisection)
    const eos::TParams tp = eos::atTemperature(m, T);
    auto two = [&](double P) { return eos::flash(m, tp, P, z).phases == 2; };
    double pOne = guess > 0 ? 1.6 * guess : 1.0e8, pTwo = -1;
    if (two(pOne)) return -1;
    for (double P = pOne / 1.04; P > 1e5; P /= 1.04) { if (two(P)) { pTwo = P; break; } pOne = P; }
    if (pTwo < 0) return -1;
    for (int it = 0; it < 50 && pOne - pTwo > 1e-7 * pOne; ++it) { const double pm = 0.5 * (pOne + pTwo); (two(pm) ? pTwo : pOne) = pm; }
    const double ps = 0.5 * (pOne + pTwo);
    const eos::FlashResult f = eos::flash(m, tp, pTwo, z);
    dew = f.phases == 2 ? f.beta > 0.5 : vap;   // incipient phase is the liquid when the mixture is mostly vapour
    return ps;
}
}  // namespace

// ============================================================================================ swelling test / P-X
SwellingResult swellingTest(const eos::Model& m, double T, const std::vector<double>& oilIn, const std::vector<double>& gasIn,
                            const std::vector<double>& alphas, double stdT, double stdP) {
    SwellingResult R;
    const size_t nc = size_t(m.nc());
    if (nc == 0) { R.message = "no components"; return R; }
    const std::vector<double> oil = normalized(oilIn, nc), gas = normalized(gasIn, nc);
    const eos::TParams tp = eos::atTemperature(m, T);
    bool dew0 = false;
    R.psatOil = psatOf(m, T, oil, 0, dew0);
    if (!(R.psatOil > 0)) { R.message = "no saturation pressure for the original oil at this temperature"; return R; }
    eos::EvalOptions eo; eo.dn = false; eo.dP = false; eo.root = dew0 ? eos::Root::Vapour : eos::Root::Liquid;
    const double vmOil = eos::evalPhase(m, tp, R.psatOil, oil, eo).vm;
    const double vStd = eos::Rgas * stdT / stdP;   // ideal-gas standard molar volume
    double guess = R.psatOil;
    for (double a : alphas) {
        if (!(a >= 0 && a < 1)) continue;
        SwellRow row; row.alpha = a; row.gasOilRatio = a / (1 - a);
        row.gasAddedStd = row.gasOilRatio * vStd / vmOil;
        const std::vector<double> z = mix(gas, oil, a);
        bool dew = false;
        row.psat = a == 0 ? R.psatOil : psatOf(m, T, z, guess, dew);
        row.dew = a == 0 ? dew0 : dew;
        if (row.psat > 0) {
            guess = row.psat;
            eos::EvalOptions e2 = eo; e2.root = row.dew ? eos::Root::Vapour : eos::Root::Liquid;
            const eos::PhaseEval pe = eos::evalPhase(m, tp, row.psat, z, e2);
            row.swelling = pe.vm / (1 - a) / vmOil;
            row.rhoSat = pe.rho;
            if (row.psat > R.fcmp) { R.fcmp = row.psat; R.fcmpAlpha = a; }
        }
        R.rows.push_back(row);
    }
    // refine the cricondenbar of the pseudo-binary (golden section around the largest tabulated saturation pressure)
    if (R.fcmp > 0 && R.rows.size() >= 3) {
        size_t k = 0;
        for (size_t i = 0; i < R.rows.size(); ++i) if (R.rows[i].psat == R.fcmp) k = i;
        if (k > 0 && k + 1 < R.rows.size() && R.rows[k + 1].psat > 0) {
            double lo = R.rows[k - 1].alpha, hi = R.rows[k + 1].alpha;
            const double g = 0.5 * (std::sqrt(5.0) - 1);
            auto f = [&](double a) { bool d; return psatOf(m, T, mix(gas, oil, a), R.fcmp, d); };
            double a1 = hi - g * (hi - lo), a2 = lo + g * (hi - lo), f1 = f(a1), f2 = f(a2);
            for (int it = 0; it < 25 && hi - lo > 1e-5; ++it) {
                if (f1 > f2) { hi = a2; a2 = a1; f2 = f1; a1 = hi - g * (hi - lo); f1 = f(a1); }
                else { lo = a1; a1 = a2; f1 = f2; a2 = lo + g * (hi - lo); f2 = f(a2); }
            }
            const double best = std::max(f1, f2);
            if (best > R.fcmp) { R.fcmp = best; R.fcmpAlpha = f1 > f2 ? a1 : a2; }
        }
    }
    R.ok = true;
    char b[200]; std::snprintf(b, sizeof b, "Psat(oil) %.4g Pa, FCMP %.4g Pa at %.1f mol%% gas", R.psatOil, R.fcmp, 100 * R.fcmpAlpha);
    R.message = b;
    return R;
}

// ============================================================================================ mixing-cell MMP
// Ahmadi & Johns (2011): fresh gas is mixed with the liquid of the first cell, the vapour of each cell moves on and
// mixes with the liquid of the next cell, the vapour of the last cell mixes with fresh oil - one more cell per contact.
// The key tie lines (plateaus of the tie-line length profile) develop with the contacts; at the MMP the shortest key
// tie line becomes critical (zero length).
double mixingCellMinTieLine(const eos::Model& m, const eos::TParams& tp, double P, const std::vector<double>& oilIn,
                            const std::vector<double>& gasIn, int contacts, double ratio, std::vector<double>* tls, int* flashes) {
    const size_t nc = size_t(m.nc());
    const std::vector<double> oil = normalized(oilIn, nc), gas = normalized(gasIn, nc);
    struct Cell { bool hasL = false, hasV = false, two = false; std::vector<double> x, y, K; double tl = 0; };
    std::vector<Cell> cells;
    auto flashCell = [&](const std::vector<double>& z, const Cell* prev) {
        Cell c;
        eos::FlashOptions fo;
        if (prev && prev->two && prev->K.size() == nc) fo.Kinit = &prev->K;
        const eos::FlashResult f = eos::flash(m, tp, P, z, fo);
        if (flashes) ++*flashes;
        if (f.phases == 2) {
            c.two = c.hasL = c.hasV = true; c.x = f.x; c.y = f.y; c.K = f.K;
            double s = 0; for (size_t i = 0; i < nc; ++i) s += (f.x[i] - f.y[i]) * (f.x[i] - f.y[i]);
            c.tl = std::sqrt(s);
        } else if (f.vapourLabel) { c.hasV = true; c.y = z; }
        else { c.hasL = true; c.x = z; }
        return c;
    };
    cells.push_back(flashCell(mix(gas, oil, ratio), nullptr));
    for (int k = 1; k < contacts; ++k) {
        std::vector<Cell> nw; nw.reserve(cells.size() + 1);
        const size_t n = cells.size();
        for (size_t i = 0; i <= n; ++i) {
            const std::vector<double>* up = i == 0 ? &gas : (cells[i - 1].hasV ? &cells[i - 1].y : nullptr);
            const std::vector<double>* dn = i == n ? &oil : (cells[i].hasL ? &cells[i].x : nullptr);
            if (!up && !dn) { nw.push_back(Cell()); continue; }
            const std::vector<double> z = !up ? *dn : !dn ? *up : mix(*up, *dn, ratio);
            nw.push_back(flashCell(z, i < n ? &cells[i] : (n ? &cells[n - 1] : nullptr)));
        }
        cells.swap(nw);
    }
    // minimum tie line over the two-phase cells; a single-phase cell inside the two-phase profile = the path crossed
    // the critical region (zero-length key tie line)
    double tmin = 1e300;
    int first = -1, last = -1;
    for (size_t i = 0; i < cells.size(); ++i) if (cells[i].two) { if (first < 0) first = int(i); last = int(i); tmin = std::min(tmin, cells[i].tl); }
    if (first < 0) tmin = 0;
    for (int i = first + 1; i < last; ++i) if (!cells[size_t(i)].two && (cells[size_t(i)].hasL || cells[size_t(i)].hasV)) tmin = 0;
    if (tls) { tls->clear(); for (auto& c : cells) tls->push_back(c.two ? c.tl : 0.0); }
    return tmin;
}

MmpResult mixingCellMmp(const eos::Model& m, double T, const std::vector<double>& oil, const std::vector<double>& gas, const MmpOptions& o) {
    MmpResult R;
    if (m.nc() == 0) { R.message = "no components"; return R; }
    const eos::TParams tp = eos::atTemperature(m, T);
    const int nP = std::max(4, o.nP);
    auto eval = [&](double P, std::vector<double>* prof) { return mixingCellMinTieLine(m, tp, P, oil, gas, std::max(5, o.contacts), o.ratio, prof, &R.flashes); };
    double tlMax = 0;
    for (int k = 0; k < nP; ++k) {
        const double P = o.pLo + (o.pHi - o.pLo) * k / (nP - 1);
        R.p.push_back(P); R.tlMin.push_back(eval(P, nullptr)); tlMax = std::max(tlMax, R.tlMin.back());
    }
    if (!(tlMax > 0)) { R.message = "no immiscible pressure in the interval: MMP is below the lower bound"; R.mmp = o.pLo; return R; }
    // MMP: zero intercept of the minimum tie-line length. TL_min falls almost linearly with pressure; close to the MMP
    // the finite number of contacts leaves a slowly converging tail (and key tie lines may switch), so the straight
    // line is fitted on the monotonic part between 50 % and 10 % of the low-pressure length (Ahmadi & Johns 2011:
    // extrapolation of the minimum tie line to zero length)
    const double tl0 = R.tlMin[0];
    const double zeroTol = 0.02 * tlMax;
    size_t seg = 0;   // end of the monotonic decreasing segment
    while (seg + 1 < R.p.size() && R.tlMin[seg + 1] < R.tlMin[seg] && R.tlMin[seg + 1] > zeroTol) ++seg;
    for (size_t k = 1; k < R.p.size(); ++k) if (R.tlMin[k] <= zeroTol) { R.mmpBracket = R.p[k]; break; }
    auto fitWindow = [&](std::vector<std::pair<double, double>>& pts) {
        pts.clear();
        for (size_t k = 0; k <= seg; ++k) if (R.tlMin[k] <= 0.5 * tl0 && R.tlMin[k] >= 0.1 * tl0) pts.push_back({R.p[k], R.tlMin[k]});
    };
    std::vector<std::pair<double, double>> pts;
    fitWindow(pts);
    if (pts.size() < 3 && seg >= 1) {   // densify the window
        double a = R.p[0], b = R.p[seg];
        for (size_t k = 0; k <= seg; ++k) { if (R.tlMin[k] > 0.5 * tl0) a = R.p[k]; }
        for (size_t k = seg + 1; k-- > 0;) { if (R.tlMin[k] < 0.1 * tl0) b = R.p[k]; }
        if (seg + 1 < R.p.size()) b = std::max(b, R.p[seg + 1]);
        for (int j = 1; j <= 4; ++j) {
            const double P = a + (b - a) * j / 5.0;
            const double t = eval(P, nullptr);
            const size_t at = size_t(std::lower_bound(R.p.begin(), R.p.end(), P) - R.p.begin());
            R.p.insert(R.p.begin() + std::ptrdiff_t(at), P); R.tlMin.insert(R.tlMin.begin() + std::ptrdiff_t(at), t);
        }
        seg = 0;
        while (seg + 1 < R.p.size() && R.tlMin[seg + 1] < R.tlMin[seg] && R.tlMin[seg + 1] > zeroTol) ++seg;
        fitWindow(pts);
    }
    if (pts.size() < 2) { pts.clear(); for (size_t k = seg >= 2 ? seg - 2 : 0; k <= seg; ++k) if (R.tlMin[k] > zeroTol) pts.push_back({R.p[k], R.tlMin[k]}); }
    if (pts.size() >= 2) {
        const double n = double(pts.size());
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (auto& q : pts) { sx += q.first; sy += q.second; sxx += q.first * q.first; sxy += q.first * q.second; }
        const double den = n * sxx - sx * sx;
        const double slope = den != 0 ? (n * sxy - sx * sy) / den : 0, icpt = (sy - slope * sx) / n;
        if (slope < 0) R.mmp = std::max(-icpt / slope, pts.back().first);
    }
    if (R.mmp < 0 && R.mmpBracket > 0) R.mmp = R.mmpBracket;
    // diagnostic tie-line profile at the highest immiscible pressure
    if (seg < R.p.size()) { R.profileP = R.p[seg]; eval(R.profileP, &R.profileTL); }
    R.ok = R.mmp > 0;
    char b[240];
    if (R.ok) std::snprintf(b, sizeof b, "mixing-cell MMP %.4g Pa (%s), %d flashes", R.mmp, R.mmp <= o.pHi ? "minimum tie line extrapolated to zero" : "extrapolated above the scanned interval", R.flashes);
    else std::snprintf(b, sizeof b, "MMP not found in [%.3g, %.3g] Pa - raise the upper bound", o.pLo, o.pHi);
    R.message = b;
    return R;
}

// ============================================================================================ slim tube
SlimTubeResult slimTubeMmp(const eos::Model& m, double T, const std::vector<double>& oilIn, const std::vector<double>& gasIn,
                           const std::vector<double>& pressures, const SlimTubeOptions& o) {
    SlimTubeResult R;
    const size_t nc = size_t(m.nc());
    const std::vector<double> oil = normalized(oilIn, nc), gas = normalized(gasIn, nc);
    const eos::TParams tp = eos::atTemperature(m, T);
    const int n = std::max(10, o.cells);
    for (double P : pressures) {
        comp::CompCase c;
        c.nx = n; c.ny = 1; c.nz = 1;
        c.dx = {o.length / n}; c.dy = {0.05}; c.dz = {0.05};
        c.poro = {0.35}; c.permx = {5e-12}; c.permy = {5e-12}; c.permz = {5e-12};
        c.fluid = m; c.T = T; c.gravity = false; c.water = false; c.swc = 0;
        c.sorg = 0; c.sgc = 0; c.no = o.corey; c.ng = o.corey; c.kro0 = 1; c.krg0 = 1;
        c.p0 = {P}; c.z0 = {oil};
        const double pv = o.length * 0.05 * 0.05 * 0.35;
        eos::EvalOptions eo; eo.dn = false; eo.dP = false;
        const eos::FlashResult fg = eos::flash(m, tp, P, gas);
        const double vmg = fg.phases == 2 ? (1 - fg.beta) * fg.L.vm + fg.beta * fg.V.vm : (fg.vapourLabel ? fg.V.vm : fg.L.vm);
        const double dur = o.days * 86400.0;
        comp::CompWell inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {1e-11};
        inj.rate = o.pvInjected * pv / vmg / dur; inj.bhpMax = 3 * P + 5e7; inj.zinj = gas;
        comp::CompWell prd; prd.name = "PRD"; prd.cells = {n - 1}; prd.wi = {1e-11}; prd.control = 0; prd.bhp = P;
        c.wells = {inj, prd};
        std::vector<double> prof;
        double rf = -1;
        try {
            comp::CompSimulator sim; sim.setup(c);
            double o0, g0, w0; sim.inPlaceSurface(o0, g0, w0);
            bool good = true;
            for (int k = 1; k <= 10 && good; ++k) {
                good = sim.advanceTo(dur * k / 10.0, dur / 40.0);
                double o1, g1, w1; sim.inPlaceSurface(o1, g1, w1);
                prof.push_back(o0 > 0 ? 1 - o1 / o0 : 0.0);
            }
            if (good && !prof.empty()) rf = prof.back();
        } catch (...) { rf = -1; }
        R.p.push_back(P); R.rf.push_back(rf); R.profiles.push_back(prof);
    }
    // MMP: break-over (two straight lines, the break point minimising the total squared error) and the RF threshold
    std::vector<std::pair<double, double>> pts;
    for (size_t k = 0; k < R.p.size(); ++k) if (R.rf[k] >= 0) pts.push_back({R.p[k], R.rf[k]});
    std::sort(pts.begin(), pts.end());
    for (size_t k = 1; k < pts.size(); ++k)
        if (pts[k - 1].second < o.rfMiscible && pts[k].second >= o.rfMiscible) {
            R.mmpRf = pts[k - 1].first + (o.rfMiscible - pts[k - 1].second) * (pts[k].first - pts[k - 1].first) / (pts[k].second - pts[k - 1].second);
            break;
        }
    auto fit = [&](size_t a, size_t b, double& sl, double& ic) {   // [a, b)
        double sx = 0, sy = 0, sxx = 0, sxy = 0; const double nn = double(b - a);
        for (size_t i = a; i < b; ++i) { sx += pts[i].first; sy += pts[i].second; sxx += pts[i].first * pts[i].first; sxy += pts[i].first * pts[i].second; }
        const double den = nn * sxx - sx * sx;
        sl = den != 0 ? (nn * sxy - sx * sy) / den : 0; ic = (sy - sl * sx) / nn;
        double e = 0; for (size_t i = a; i < b; ++i) { const double r = pts[i].second - (sl * pts[i].first + ic); e += r * r; }
        return e;
    };
    double best = 1e300;
    for (size_t k = 2; k + 2 <= pts.size(); ++k) {
        double s1, i1, s2, i2;
        const double e = fit(0, k, s1, i1) + fit(k, pts.size(), s2, i2);
        if (e < best && s1 > s2 && std::fabs(s1 - s2) > 1e-30) {
            const double x = (i2 - i1) / (s1 - s2);
            if (x >= pts.front().first && x <= pts.back().first) { best = e; R.mmpBreak = x; }
        }
    }
    R.ok = R.mmpBreak > 0 || R.mmpRf > 0;
    char b[200];
    std::snprintf(b, sizeof b, "slim tube (%d cells, %.2f HCPV): break-over %.4g Pa, RF %.0f %% at %.4g Pa", n, o.pvInjected, R.mmpBreak, 100 * o.rfMiscible, R.mmpRf);
    R.message = b;
    return R;
}

// ============================================================================================ correlations
double mmpYelligMetcalfe(double T) {
    const double tf = T * 9.0 / 5.0 - 459.67;
    const double psia = 1833.7217 + 2.2518055 * tf + 0.01800674 * tf * tf - 103949.93 / tf;
    return eos::psiaToPa(psia);
}
double mmpCronquist(double T, double mwC5plus, double yLight) {
    const double tf = T * 9.0 / 5.0 - 459.67;
    const double psia = 15.988 * std::pow(tf, 0.744206 + 0.0011038 * mwC5plus + 0.0015279 * 100.0 * yLight);
    return eos::psiaToPa(psia);
}

// ============================================================================================ WAG schedule
int applyWagSchedule(Case& c, const WagSpec& w, std::string* message) {
    auto fail = [&](const std::string& s) { if (message) *message = s; return 0; };
    const WellInput* wi = nullptr;
    for (const auto& q : c.wells) if (q.name == w.well) wi = &q;
    if (!wi) return fail("WAG: well '" + w.well + "' not found");
    if (wi->type == WellType::Producer) return fail("WAG: '" + w.well + "' is a producer - choose an injector");
    if (w.cycles <= 0 || !(w.gasDays > 0) || !(w.waterDays > 0) || w.start < 0) return fail("WAG: cycles and half-cycle lengths must be positive");
    const double eps = 1e-6;
    auto starts = [&]() { std::vector<double> t; double s = 0; for (const auto& r : c.schedule) { t.push_back(s); s += r.length; } t.push_back(s); return t; };
    // index of the report step starting at time t (splitting / appending steps as needed)
    auto stepAt = [&](double t) -> size_t {
        std::vector<double> st = starts();
        const double end = st.back();
        if (t >= end - eps) {
            if (t > end + eps) { ReportStep r; r.length = t - end; c.schedule.push_back(r); }
            ReportStep r; r.length = std::max(1.0, c.schedule.empty() ? 30.0 : c.schedule.back().length);
            c.schedule.push_back(r);   // the step that starts at t (gets trimmed by the next boundary when one follows)
            return c.schedule.size() - 1;
        }
        for (size_t s = 0; s + 1 < st.size(); ++s) {
            if (std::fabs(st[s] - t) <= eps) return s;
            if (t > st[s] && t < st[s + 1]) {
                ReportStep tail; tail.length = st[s + 1] - t;
                c.schedule[s].length = t - st[s];
                c.schedule.insert(c.schedule.begin() + std::ptrdiff_t(s) + 1, tail);
                return s + 1;
            }
        }
        return c.schedule.size() - 1;
    };
    // boundaries
    std::vector<std::pair<double, int>> marks;   // time, 2 = gas, 1 = water, 0 = shut
    double t = w.start, gl = w.gasDays;
    for (int k = 0; k < w.cycles; ++k) {
        for (int h = 0; h < 2; ++h) {
            const bool gas = (h == 0) == w.gasFirst;
            marks.push_back({t, gas ? 2 : 1});
            t += gas ? gl : w.waterDays;
        }
        gl *= w.taper > 0 ? w.taper : 1.0;
    }
    if (w.chaseWaterDays > 0) { marks.push_back({t, 1}); t += w.chaseWaterDays; }
    const double tEnd = t;
    // the schedule must reach the end of the last half-cycle: split / append boundaries first, then trim a trailing
    // helper step that was only appended to host a boundary
    const double oldEnd = starts().back();
    if (tEnd > oldEnd + eps) {   // extend with steps of the last report length (the boundaries split them below)
        const double len = std::max(1.0, c.schedule.empty() ? 30.0 : c.schedule.back().length);
        double e = oldEnd;
        while (e < tEnd - eps) { ReportStep r; r.length = std::min(len, tEnd - e); if (tEnd - e - r.length < 0.25 * len) r.length = tEnd - e; c.schedule.push_back(r); e += r.length; }
    }
    for (const auto& mk : marks) stepAt(mk.first);
    if (w.shutAfter && tEnd >= oldEnd - eps) stepAt(tEnd);   // one step after the WAG period to host the shut-in
    // remove the well's controls / WINJGAS from the start on
    {
        const std::vector<double> st = starts();
        for (size_t s = 0; s < c.schedule.size(); ++s) {
            if (st[s] < w.start - eps) continue;
            auto& r = c.schedule[s];
            r.controls.erase(std::remove_if(r.controls.begin(), r.controls.end(), [&](const WellControl& q) { return q.well == w.well; }), r.controls.end());
            r.injStreams.erase(std::remove_if(r.injStreams.begin(), r.injStreams.end(), [&](const std::pair<std::string, std::string>& q) { return q.first == w.well; }), r.injStreams.end());
        }
    }
    int written = 0;
    for (const auto& mk : marks) {
        const size_t s = stepAt(mk.first);
        WellControl wc; wc.well = w.well; wc.mode = Control::RATE;
        wc.rate = mk.second == 2 ? w.gasRate : w.waterRate; wc.bhp = w.bhpMax > 0 ? w.bhpMax : -1; wc.phase = mk.second;
        c.schedule[s].controls.push_back(wc);
        if (mk.second == 2 && !w.stream.empty()) c.schedule[s].injStreams.push_back({w.well, w.stream});
        ++written;
    }
    if (w.shutAfter && tEnd < starts().back() - eps) {
        const size_t s = stepAt(tEnd);
        WellControl wc; wc.well = w.well; wc.status = WellStatus::Shut; wc.mode = Control::RATE; wc.rate = 0;
        c.schedule[s].controls.push_back(wc);
    }
    if (message) {
        char b[200]; std::snprintf(b, sizeof b, "WAG on %s: %d half-cycles from day %.1f to day %.1f (%d report steps)", w.well.c_str(), written, w.start, tEnd, int(c.schedule.size()));
        *message = b;
    }
    return written;
}

}  // namespace eor
}  // namespace ressim
