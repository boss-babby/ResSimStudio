// C20: rate-and-state seismicity forecasting + convex minimum-seismicity injection designer (see seismicity.h)
#include "ressim/seismicity.h"
#include "ressim/simulator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <cstdio>
#include <limits>

namespace ressim {
namespace seis {

namespace {
constexpr double kYear = 365.25;
inline double logaddexp0(double x) { return x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)); }
inline double sigmoid(double x) { return x >= 0 ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x)); }
// exact integral of exp(linear) over a segment of length ds with end exponents x0, x1 (already shifted), and its
// partial derivatives with respect to x0 and x1
inline void segInt(double ds, double x0, double x1, double& S, double& d0, double& d1) {
    const double d = x1 - x0, e0 = std::exp(x0);
    if (std::fabs(d) < 1e-3) {
        S = ds * e0 * (1 + d / 2 + d * d / 6);
        d0 = ds * e0 * (0.5 + d / 6 + d * d / 24);
        d1 = ds * e0 * (0.5 + d / 3 + d * d / 8);
    } else {
        const double e1 = std::exp(x1), h = (e1 - e0) / d;
        S = ds * h; d0 = ds * (h - e0) / d; d1 = ds * (e1 - h) / d;
    }
}
struct PatchEval { double E = 0; double lnI = 0; };
// excess events of one patch for dS (MPa) at times t (days); optional dE/d dS_i
PatchEval patchExcess(const std::vector<double>& t, const double* dS, const Patch& p, double* dEddS) {
    const size_t n = t.size();
    const double a = std::max(1e-9, p.aSigma), ta = std::max(1e-9, p.taYears * kYear), r = p.rate / kYear;
    std::vector<double> e(n);
    double mx = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) { e[i] = (t[i] - t[0]) / ta + dS[i] / a; mx = std::max(mx, e[i]); }
    double sum = 0;
    std::vector<double> de(dEddS ? n : 0, 0.0);
    for (size_t k = 0; k + 1 < n; ++k) {
        double S, d0, d1; segInt(t[k + 1] - t[k], e[k] - mx, e[k + 1] - mx, S, d0, d1);
        sum += S;
        if (dEddS) { de[k] += d0; de[k + 1] += d1; }
    }
    PatchEval out;
    const double T = t[n - 1] - t[0];
    if (sum <= 0) { out.lnI = -1e300; out.E = -r * T; if (dEddS) std::fill(dEddS, dEddS + n, 0.0); return out; }
    out.lnI = mx + std::log(sum);
    const double z = out.lnI - std::log(ta);
    out.E = r * ta * (logaddexp0(z) - T / ta);
    if (dEddS) {
        const double coef = r * ta * sigmoid(z) / sum / a;   // dE/de_i = r ta sigma(z) (dI/de_i)/I ; de_i/ddS_i = 1/a
        for (size_t i = 0; i < n; ++i) dEddS[i] = coef * de[i];
    }
    return out;
}
} // namespace

const char* toString(Objective o) {
    switch (o) { case Objective::Events: return "Expected events"; case Objective::Peak: return "Peak rate"; default: return "Blend"; }
}

Forecast forecast(const std::vector<double>& t, const std::vector<double>& dS, const Patch& p, double tShutIn) {
    Forecast f; f.t = t; f.dS = dS;
    const size_t n = t.size();
    if (n == 0 || dS.size() != n) return f;
    const double a = std::max(1e-9, p.aSigma), ta = std::max(1e-9, p.taYears * kYear), r = p.rate / kYear;
    std::vector<double> e(n);
    double mx = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) { e[i] = (t[i] - t[0]) / ta + dS[i] / a; mx = std::max(mx, e[i]); f.maxDS = i ? std::max(f.maxDS, dS[i]) : dS[i]; }
    f.R.resize(n); f.N.resize(n);
    double sum = 0;
    for (size_t i = 0; i < n; ++i) {
        if (i) { double S, d0, d1; segInt(t[i] - t[i - 1], e[i - 1] - mx, e[i] - mx, S, d0, d1); sum += S; }
        const double z = sum > 0 ? mx + std::log(sum) - std::log(ta) : -1e300;
        const double L = sum > 0 ? logaddexp0(z) : 0.0;
        f.R[i] = p.rate * std::exp(e[i] - L);                 // events / year
        f.N[i] = r * ta * L;                                  // events
        if (f.R[i] > f.peakR) { f.peakR = f.R[i]; f.tPeak = t[i]; }
    }
    f.excess = f.N[n - 1] - r * (t[n - 1] - t[0]);
    f.excessAtShutIn = f.excess;
    if (tShutIn > t[0] && tShutIn < t[n - 1]) {
        size_t k = size_t(std::upper_bound(t.begin(), t.end(), tShutIn) - t.begin());
        k = std::min(std::max<size_t>(k, 1), n - 1);
        const double w = (tShutIn - t[k - 1]) / std::max(1e-300, t[k] - t[k - 1]);
        const double N = f.N[k - 1] + w * (f.N[k] - f.N[k - 1]);
        f.excessAtShutIn = N - r * (tShutIn - t[0]);
    }
    return f;
}

std::vector<double> modelStress(const LinearModel& m, const std::vector<double>& x, int p) {
    const size_t n = m.t.size(), nv = size_t(m.nv);
    std::vector<double> s(n);
    const auto& G = m.G[size_t(p)];
    for (size_t i = 0; i < n; ++i) {
        double v = m.dS0.empty() ? 0.0 : m.dS0[size_t(p)][i];
        const double* g = &G[i * nv];
        for (size_t j = 0; j < nv; ++j) v += g[j] * x[j];
        s[i] = v;
    }
    return s;
}

double modelEvents(const LinearModel& m, const std::vector<double>& x, std::vector<double>* grad) {
    const size_t n = m.t.size(), nv = size_t(m.nv);
    if (grad) grad->assign(nv, 0.0);
    double E = 0;
    std::vector<double> dE(n);
    for (size_t p = 0; p < m.patches.size(); ++p) {
        const std::vector<double> s = modelStress(m, x, int(p));
        const double w = m.patches[p].weight;
        const PatchEval pe = patchExcess(m.t, s.data(), m.patches[p], grad ? dE.data() : nullptr);
        E += w * pe.E;
        if (grad) {
            const auto& G = m.G[p];
            for (size_t i = 0; i < n; ++i) {
                const double c = w * dE[i];
                if (c == 0) continue;
                const double* g = &G[i * nv];
                for (size_t j = 0; j < nv; ++j) (*grad)[j] += c * g[j];
            }
        }
    }
    return E;
}

double modelPeak(const LinearModel& m, const std::vector<double>& x, double beta, std::vector<double>* grad) {
    const size_t n = m.t.size(), nv = size_t(m.nv);
    std::vector<std::vector<double>> z(m.patches.size());
    double mx = -std::numeric_limits<double>::infinity();
    for (size_t p = 0; p < m.patches.size(); ++p) {
        z[p] = modelStress(m, x, int(p));
        for (double& v : z[p]) { v /= std::max(1e-9, m.patches[p].aSigma); mx = std::max(mx, v); }
    }
    double sum = 0;
    for (const auto& zp : z) for (double v : zp) sum += std::exp(beta * (v - mx));
    const double P = mx + std::log(sum) / beta;
    if (grad) {
        grad->assign(nv, 0.0);
        for (size_t p = 0; p < m.patches.size(); ++p) {
            const double ia = 1.0 / std::max(1e-9, m.patches[p].aSigma);
            for (size_t i = 0; i < n; ++i) {
                const double wgt = std::exp(beta * (z[p][i] - mx)) / sum * ia;
                if (wgt < 1e-300) continue;
                const double* g = &m.G[p][i * nv];
                for (size_t j = 0; j < nv; ++j) (*grad)[j] += wgt * g[j];
            }
        }
    }
    return P;
}

void project(const Constraints& c, std::vector<double>& x) {
    const size_t nv = x.size();
    int ng = int(c.target.size());
    for (size_t v = 0; v < nv; ++v) if (c.group.empty() || c.group[v] < 0 || c.group[v] >= ng) x[v] = std::min(std::max(x[v], c.lb[v]), c.ub[v]);
    for (int g = 0; g < ng; ++g) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        bool any = false;
        for (size_t v = 0; v < nv; ++v) {
            if (c.group[v] != g) continue;
            any = true;
            const double w = c.weight[v];
            lo = std::min(lo, (x[v] - c.ub[v]) / w); hi = std::max(hi, (x[v] - c.lb[v]) / w);
        }
        if (!any) continue;
        auto sumAt = [&](double lam) {
            double s = 0;
            for (size_t v = 0; v < nv; ++v) if (c.group[v] == g) s += c.weight[v] * std::min(std::max(x[v] - lam * c.weight[v], c.lb[v]), c.ub[v]);
            return s;
        };
        lo -= 1; hi += 1;
        for (int it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (sumAt(mid) > c.target[size_t(g)]) lo = mid; else hi = mid;
            if (hi - lo < 1e-15 * (1 + std::fabs(lo))) break;
        }
        const double lam = 0.5 * (lo + hi);
        for (size_t v = 0; v < nv; ++v) if (c.group[v] == g) x[v] = std::min(std::max(x[v] - lam * c.weight[v], c.lb[v]), c.ub[v]);
    }
}

OptResult optimize(const LinearModel& m, const Constraints& c, const std::vector<double>& x0, const OptOptions& o) {
    OptResult res;
    const size_t nv = size_t(m.nv);
    std::vector<double> x = x0; project(c, x);
    const double E0 = modelEvents(m, x), P0 = modelPeak(m, x, o.peakSharpness);
    const double sE = std::max(1e-12, std::fabs(E0)), sP = std::max(1e-3, std::fabs(P0));
    auto F = [&](const std::vector<double>& v, std::vector<double>* g) {
        std::vector<double> ge, gp;
        double f = 0;
        if (o.obj != Objective::Peak) { const double e = modelEvents(m, v, g ? &ge : nullptr); f += (o.obj == Objective::Blend ? 1 - o.blend : 1.0) * e / sE; }
        if (o.obj != Objective::Events) { const double p = modelPeak(m, v, o.peakSharpness, g ? &gp : nullptr); f += (o.obj == Objective::Blend ? o.blend : 1.0) * p / sP; }
        if (g) {
            g->assign(nv, 0.0);
            const double we = o.obj == Objective::Blend ? (1 - o.blend) / sE : 1.0 / sE, wp = o.obj == Objective::Blend ? o.blend / sP : 1.0 / sP;
            for (size_t j = 0; j < nv; ++j) (*g)[j] = (ge.empty() ? 0.0 : we * ge[j]) + (gp.empty() ? 0.0 : wp * gp[j]);
        }
        return f;
    };
    std::vector<double> g, gy, y = x, xn(nv), xold = x;
    double fx = F(x, &g);
    res.F0 = fx; res.events0 = E0; res.peak0 = P0;
    double xs = 0, gn = 0;
    for (size_t j = 0; j < nv; ++j) { xs = std::max(xs, std::fabs(x[j])); gn = std::max(gn, std::fabs(g[j])); }
    for (size_t j = 0; j < nv; ++j) xs = std::max(xs, c.ub[j] < 1e300 ? c.ub[j] : 0.0);
    double L = gn > 0 ? gn / (0.05 * std::max(xs, 1e-12)) : 1.0;
    double tk = 1, fy = fx; gy = g;
    int stall = 0;
    for (int it = 0; it < o.maxIter; ++it) {
        res.iters = it + 1;
        double fn = 0;
        for (int bt = 0; bt < 60; ++bt) {
            for (size_t j = 0; j < nv; ++j) xn[j] = y[j] - gy[j] / L;
            project(c, xn);
            fn = F(xn, nullptr);
            double q = fy;
            for (size_t j = 0; j < nv; ++j) { const double d = xn[j] - y[j]; q += gy[j] * d + 0.5 * L * d * d; }
            if (fn <= q + 1e-15 * std::fabs(q)) break;
            L *= 2;
        }
        if (fn > fx) {                   // adaptive restart (O'Donoghue & Candes)
            y = x; fy = fx; gy = g; tk = 1;
            if (++stall > 3) break;
            continue;
        }
        double dx = 0, xm = 1e-300;
        for (size_t j = 0; j < nv; ++j) { dx = std::max(dx, std::fabs(xn[j] - x[j])); xm = std::max(xm, std::fabs(xn[j])); }
        const double rel = (fx - fn) / std::max(1e-300, std::fabs(fx));
        const double tk1 = 0.5 * (1 + std::sqrt(1 + 4 * tk * tk));
        for (size_t j = 0; j < nv; ++j) y[j] = xn[j] + (tk - 1) / tk1 * (xn[j] - x[j]);
        project(c, y);
        x = xn; fx = fn; tk = tk1;
        fy = F(y, &gy);
        L *= 0.95;
        if (dx / xm < o.tol && rel < o.tol) { stall++; if (stall > 5) { res.converged = true; break; } } else stall = 0;
    }
    fx = F(x, &g);
    // KKT residual: projected-gradient step relative to the step itself
    std::vector<double> pg(nv);
    for (size_t j = 0; j < nv; ++j) pg[j] = x[j] - g[j] / L;
    project(c, pg);
    double num = 0, den = 1e-300;
    for (size_t j = 0; j < nv; ++j) { num = std::max(num, std::fabs(pg[j] - x[j])); den = std::max(den, std::fabs(x[j])); }
    res.kkt = num / den;
    if (!res.converged && res.kkt < 1e-6) res.converged = true;
    res.x = x; res.F = fx; res.events = modelEvents(m, x); res.peak = modelPeak(m, x, o.peakSharpness);
    return res;
}

// ------------------------------------------------------------------------------------------- engine-driven designer
namespace {
std::vector<double> stepLengths(const DesignSetup& s) {
    std::vector<double> L;
    const int ni = std::max(1, s.nCtrl) * std::max(1, s.subSteps);
    const double dtI = s.injDays / ni;
    for (int i = 0; i < ni; ++i) L.push_back(dtI);
    const int np = std::max(0, s.postSteps);
    if (np > 0 && s.postDays > 0) {
        double rho = 1;
        if (s.postDays > np * dtI * 1.0001) {
            double lo = 1, hi = 50;
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                const double sum = dtI * (std::pow(mid, np) - 1) / (mid - 1);
                if (sum > s.postDays) hi = mid; else lo = mid;
            }
            rho = 0.5 * (lo + hi);
            double l = dtI, tot = 0;
            for (int i = 0; i < np; ++i) { L.push_back(l); tot += l; l *= rho; }
            L.back() += s.postDays - tot;
        } else {
            for (int i = 0; i < np; ++i) L.push_back(s.postDays / np);
        }
    }
    return L;
}
const WellControl* firstControl(const Case& c, const std::string& w) {
    for (const auto& st : c.schedule) for (const auto& wc : st.controls) if (wc.well == w) return &wc;
    return nullptr;
}
int wellIndex(const Case& c, const std::string& w) {
    for (size_t i = 0; i < c.wells.size(); ++i) if (c.wells[i].name == w) return int(i);
    return -1;
}
Case buildCase(const Case& base, const DesignSetup& s, const std::vector<std::vector<double>>& q, bool continueAfter) {
    Case c = base;
    const std::vector<double> L = stepLengths(s);
    c.schedule.assign(L.size(), ReportStep{});
    for (size_t i = 0; i < L.size(); ++i) c.schedule[i].length = L[i];
    if (!base.schedule.empty()) {
        c.schedule[0].gprod = base.schedule[0].gprod; c.schedule[0].ginj = base.schedule[0].ginj;
        c.schedule[0].guiderat = base.schedule[0].guiderat;
    }
    for (const auto& w : base.wells) {
        if (std::find(s.wells.begin(), s.wells.end(), w.name) != s.wells.end()) continue;
        if (const WellControl* wc = firstControl(base, w.name)) c.schedule[0].controls.push_back(*wc);
    }
    const int sub = std::max(1, s.subSteps);
    double qmax = 0;
    for (const auto& qw : q) for (double v : qw) qmax = std::max(qmax, v);
    for (size_t wi = 0; wi < s.wells.size(); ++wi) {
        const int idx = wellIndex(base, s.wells[wi]);
        if (idx < 0) throw std::runtime_error("seismicity design: unknown well " + s.wells[wi]);
        if (base.wells[size_t(idx)].type == WellType::Producer) throw std::runtime_error("seismicity design: well " + s.wells[wi] + " is not an injector");
        const WellControl* b = firstControl(base, s.wells[wi]);
        for (int j = 0; j < s.nCtrl; ++j) {
            WellControl wc;
            if (b) { wc.bhp = b->bhp; wc.phase = b->phase; }
            wc.well = s.wells[wi]; wc.mode = Control::RATE; wc.history = false; wc.deckMode.clear();
            wc.rate = q[wi][size_t(j)];
            wc.status = wc.rate > 1e-12 * std::max(qmax, 1e-300) ? WellStatus::Open : WellStatus::Shut;
            if (wc.status == WellStatus::Shut) wc.rate = 0;
            c.schedule[size_t(j * sub)].controls.push_back(wc);
        }
        const size_t post = size_t(s.nCtrl * sub);
        if (post < c.schedule.size()) {
            WellControl wc;
            if (b) { wc.bhp = b->bhp; wc.phase = b->phase; }
            wc.well = s.wells[wi]; wc.mode = Control::RATE;
            if (continueAfter) { wc.rate = q[wi].back(); wc.status = wc.rate > 0 ? WellStatus::Open : WellStatus::Shut; }
            else { wc.rate = 0; wc.status = WellStatus::Shut; }
            c.schedule[post].controls.push_back(wc);
        }
    }
    return c;
}
double toMPa(const Case& c) { return c.units == UnitSystem::Field ? 0.00689475729 : 0.1; }
struct Capture { std::vector<double> t; std::vector<std::vector<double>> p; };   // p[patch][i] (MPa)
Capture runCapture(const Case& c, const std::vector<Patch>& patches, std::atomic<bool>* cancel) {
    Simulator sim;
    sim.cancel = cancel;
    sim.setup(c);
    sim.initialize();
    while (!sim.finished()) if (!sim.runReportStep()) break;
    if (!sim.finished()) throw std::runtime_error("seismicity design: simulation cancelled or failed");
    const Grid& g = sim.grid();
    std::vector<int> act;
    for (const auto& p : patches) {
        if (p.i < 0 || p.j < 0 || p.k < 0 || p.i >= g.nx || p.j >= g.ny || p.k >= g.nz) throw std::runtime_error("seismicity design: patch " + p.name + " outside the grid");
        const int a = g.activeOf[size_t(p.i + g.nx * (p.j + g.ny * p.k))];
        if (a < 0) throw std::runtime_error("seismicity design: patch " + p.name + " is in an inactive cell");
        act.push_back(a);
    }
    Capture cap; cap.p.resize(patches.size());
    const double s = toMPa(c);
    for (const auto& sn : sim.snapshots()) {
        cap.t.push_back(sn.time);
        for (size_t k = 0; k < patches.size(); ++k) cap.p[k].push_back(double(sn.pressure[size_t(act[k])]) * s);
    }
    return cap;
}
double interp(const std::vector<double>& t, const std::vector<double>& v, double x) {
    if (x <= t.front()) return x < t.front() ? 0.0 : v.front();
    if (x >= t.back()) return v.back();
    const size_t k = size_t(std::upper_bound(t.begin(), t.end(), x) - t.begin());
    const double w = (x - t[k - 1]) / (t[k] - t[k - 1]);
    return v[k - 1] + w * (v[k] - v[k - 1]);
}
} // namespace

Case designCase(const Case& base, const DesignSetup& s, const std::vector<std::vector<double>>& q) { return buildCase(base, s, q, false); }

std::string designScheduleDeck(const Case& base, const DesignSetup& s, const std::vector<std::vector<double>>& q) {
    std::string out;
    char b[512];
    std::snprintf(b, sizeof b, "-- ResSim C20 seismicity-aware injection design: %d intervals of %.6g d, then %.6g d shut-in horizon\n",
                  s.nCtrl, s.injDays / std::max(1, s.nCtrl), s.postDays);
    out += b;
    const double dt = s.injDays / std::max(1, s.nCtrl);
    auto phaseOf = [&](const std::string& w) {
        const int i = wellIndex(base, w);
        return i >= 0 && base.wells[size_t(i)].type == WellType::GasInjector ? "GAS" : "WATER";
    };
    for (int j = 0; j < s.nCtrl; ++j) {
        out += "WCONINJE\n";
        for (size_t w = 0; w < s.wells.size(); ++w) {
            const WellControl* bc = firstControl(base, s.wells[w]);
            const double qq = q[w][size_t(j)];
            if (qq > 0) {
                if (bc && bc->bhp > 0) std::snprintf(b, sizeof b, " '%s' %s OPEN RATE %.6g 1* %.6g /\n", s.wells[w].c_str(), phaseOf(s.wells[w]), qq, bc->bhp);
                else std::snprintf(b, sizeof b, " '%s' %s OPEN RATE %.6g /\n", s.wells[w].c_str(), phaseOf(s.wells[w]), qq);
            } else std::snprintf(b, sizeof b, " '%s' %s SHUT RATE 0 /\n", s.wells[w].c_str(), phaseOf(s.wells[w]));
            out += b;
        }
        std::snprintf(b, sizeof b, "/\nTSTEP\n %d*%.6g /\n\n", std::max(1, s.subSteps), dt / std::max(1, s.subSteps));
        out += b;
    }
    out += "WCONINJE\n";
    for (const auto& w : s.wells) { std::snprintf(b, sizeof b, " '%s' %s SHUT RATE 0 /\n", w.c_str(), phaseOf(w)); out += b; }
    out += "/\n";
    if (s.postDays > 0) { std::snprintf(b, sizeof b, "TSTEP\n %.6g /\n", s.postDays); out += b; }
    return out;
}

// ---------------------------------------------------------------- C22 (H1f) closed-form schedule gain
namespace {
template <class F> double simpson(F f, double a, double b, int n) {
    if (!(b > a)) return 0.0;
    n = std::max(2, n + (n & 1));
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += f(a + i * h) * (i & 1 ? 4.0 : 2.0);
    return s * h / 3.0;
}
// int_0^inf (exp(P e^{-x}) - 1) dx
double tailP(double P) {
    if (P <= 0) return 0.0;
    auto f = [P](double x) { return std::expm1(P * std::exp(-x)); };
    return simpson(f, 0.0, 4.0, 800) + simpson(f, 4.0, 60.0, 600);
}
struct ArcCost { double cost, sp, se; };
// ramp at u to s*, plateau s*, burst at u of length tb, decay; costs of excess exp(Pi s/sc) - 1 (Ti = 1)
ArcCost threeArc(double Pi, double T, double u, double tb) {
    const double sc = -std::expm1(-1.0 / T);
    auto volErr = [&](double ss) { const double tr = T * std::log(u / (u - ss)); return u * tr + ss * (1.0 - tr - tb) + u * tb - 1.0; };
    double lo = 1e-12, hi = std::min(u * (1 - 1e-12), u * -std::expm1(-(1.0 - tb) / T));
    if (volErr(lo) > 0 || volErr(hi) < 0) return {std::numeric_limits<double>::infinity(), 0, 0};
    for (int it = 0; it < 200; ++it) { const double mid = 0.5 * (lo + hi); (volErr(mid) > 0 ? hi : lo) = mid; }
    const double ss = 0.5 * (lo + hi), tr = T * std::log(u / (u - ss));
    if (tr + tb > 1.0) return {std::numeric_limits<double>::infinity(), 0, 0};
    const double ramp = simpson([&](double t) { return std::expm1(Pi * u * -std::expm1(-t / T) / sc); }, 0.0, tr, 400);
    const double plat = (1.0 - tr - tb) * std::expm1(Pi * ss / sc);
    const double burst = simpson([&](double t) { return std::expm1(Pi * (u - (u - ss) * std::exp(-t / T)) / sc); }, 0.0, tb, 400);
    const double se = u - (u - ss) * std::exp(-tb / T);
    return {ramp + plat + burst + T * tailP(Pi * se / sc), ss / sc, se / sc};
}
} // namespace

GainEstimate scheduleGain(double Pi, double theta, double uMax) {
    GainEstimate g;
    if (!(Pi > 0) || !(theta > 0) || !(uMax > 1.0)) return g;
    const double T = theta, u = uMax, sc = -std::expm1(-1.0 / T);
    // constant rate: s = 1 - e^{-t/T}
    const double t1 = std::min(1.0, 12 * T);
    auto fc = [&](double t) { return std::expm1(Pi * -std::expm1(-t / T) / sc); };
    const double Jc = simpson(fc, 0.0, t1, 1200) + simpson(fc, t1, 1.0, 400) + T * tailP(Pi);
    const ArcCost c0 = threeArc(Pi, T, u, 0.0);
    // golden section over the burst length
    double a = 0.0, b = std::min(0.5, 3 * T);
    const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = b - gr * (b - a), x2 = a + gr * (b - a);
    double f1 = threeArc(Pi, T, u, x1).cost, f2 = threeArc(Pi, T, u, x2).cost;
    for (int it = 0; it < 80; ++it) {
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = threeArc(Pi, T, u, x1).cost; }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = threeArc(Pi, T, u, x2).cost; }
    }
    const double tb = 0.5 * (a + b);
    ArcCost c3 = threeArc(Pi, T, u, tb);
    double tbBest = tb;
    if (!(c3.cost <= c0.cost)) { c3 = c0; tbBest = 0; }
    if (std::isfinite(c3.cost) && Jc > 0) {
        g.gain = std::max(0.0, 1.0 - c3.cost / Jc);        // constant rate is always feasible: a negative value means the 3-arc family is wrong here
        g.gainNoBurst = std::isfinite(c0.cost) ? std::max(0.0, 1.0 - c0.cost / Jc) : 0.0;
        g.exact = theta <= 0.1 || (uMax <= 3.0 + 1e-9 && theta <= 0.4);   // CYCLE4_REPORT F1-F3 validation range
        g.burstFraction = u * tbBest / (1.0 + 0.0);   // volume fraction injected in the burst (total = 1 per unit time window)
        g.sPlateau = c3.sp; g.sEnd = c3.se;
    }
    return g;
}

DesignOutput designInjection(const Case& base, const DesignSetup& s, std::atomic<bool>* cancel,
                             const std::function<void(const std::string&)>& log) {
    DesignOutput out;
    const auto t0 = std::chrono::steady_clock::now();
    auto say = [&](const std::string& m) { out.logLines.push_back(m); if (log) log(m); };
    try {
        const int nw = int(s.wells.size()), nc = std::max(1, s.nCtrl);
        if (nw == 0) throw std::runtime_error("seismicity design: no injector selected");
        if (s.patches.empty()) throw std::runtime_error("seismicity design: no fault patch defined");
        if (!s.wellVolume.empty() && int(s.wellVolume.size()) != nw) throw std::runtime_error("seismicity design: one volume per well expected");
        if (s.wellVolume.empty() && !(s.totalVolume > 0)) throw std::runtime_error("seismicity design: total volume must be > 0");
        const double dtc = s.injDays / nc;
        out.ctrlStart.resize(size_t(nc + 1));
        for (int j = 0; j <= nc; ++j) out.ctrlStart[size_t(j)] = j * dtc;
        std::vector<double> meanRate(static_cast<size_t>(nw));
        for (int w = 0; w < nw; ++w) meanRate[size_t(w)] = (s.wellVolume.empty() ? s.totalVolume / nw : s.wellVolume[size_t(w)]) / s.injDays;
        // 1. base run (designed wells shut) and one step-response run per well
        std::vector<std::vector<double>> qz(size_t(nw), std::vector<double>(size_t(nc), 0.0));
        say("Base run (designed injectors shut) ...");
        const Capture base0 = runCapture(buildCase(base, s, qz, false), s.patches, cancel); out.runs++;
        const size_t nt = base0.t.size();
        std::vector<std::vector<std::vector<double>>> step(static_cast<size_t>(nw));   // [well][patch][i] MPa per unit rate
        for (int w = 0; w < nw; ++w) {
            auto q = qz; for (double& v : q[size_t(w)]) v = meanRate[size_t(w)];
            say("Step-response run for " + s.wells[size_t(w)] + " ...");
            const Capture cw = runCapture(buildCase(base, s, q, true), s.patches, cancel); out.runs++;
            if (cw.t.size() != nt) throw std::runtime_error("seismicity design: report times differ between runs");
            step[size_t(w)].resize(s.patches.size());
            for (size_t p = 0; p < s.patches.size(); ++p) {
                auto& sp = step[size_t(w)][p]; sp.resize(nt);
                const double k = s.patches[p].friction * s.patches[p].cfsFactor / meanRate[size_t(w)];
                for (size_t i = 0; i < nt; ++i) sp[i] = k * (cw.p[p][i] - base0.p[p][i]);
            }
        }
        // 2. linear model (time-invariant superposition of step responses)
        LinearModel m;
        m.t = base0.t; m.patches = s.patches; m.nv = nw * nc;
        m.dS0.resize(s.patches.size()); m.G.resize(s.patches.size());
        for (size_t p = 0; p < s.patches.size(); ++p) {
            const double k = s.patches[p].friction * s.patches[p].cfsFactor;
            m.dS0[p].resize(nt);
            for (size_t i = 0; i < nt; ++i) m.dS0[p][i] = k * (base0.p[p][i] - base0.p[p][0]);
            m.G[p].assign(nt * size_t(m.nv), 0.0);
            for (size_t i = 0; i < nt; ++i)
                for (int w = 0; w < nw; ++w)
                    for (int j = 0; j < nc; ++j) {
                        const double ti = m.t[i] - m.t[0];
                        const double v = interp(m.t, step[size_t(w)][p], ti - j * dtc + m.t[0]) - interp(m.t, step[size_t(w)][p], ti - (j + 1) * dtc + m.t[0]);
                        m.G[p][i * size_t(m.nv) + size_t(w * nc + j)] = v;
                    }
        }
        // 3. constraints and the convex optimisation
        Constraints c;
        c.lb.assign(size_t(m.nv), 0.0); c.ub.resize(size_t(m.nv)); c.group.resize(size_t(m.nv)); c.weight.assign(size_t(m.nv), dtc);
        for (int w = 0; w < nw; ++w)
            for (int j = 0; j < nc; ++j) {
                const size_t v = size_t(w * nc + j);
                c.ub[v] = std::max(1.0, s.maxRateFactor) * meanRate[size_t(w)];
                if (j == nc - 1 && s.lastRateFactor > 0) c.ub[v] = std::min(c.ub[v], s.lastRateFactor * meanRate[size_t(w)]);
                c.group[v] = s.wellVolume.empty() ? 0 : w;
            }
        if (s.wellVolume.empty()) c.target = {s.totalVolume};
        else c.target = s.wellVolume;
        std::vector<double> x0(size_t(m.nv));
        for (int w = 0; w < nw; ++w) for (int j = 0; j < nc; ++j) x0[size_t(w * nc + j)] = meanRate[size_t(w)];
        say("Convex optimisation (" + std::string(toString(s.opt.obj)) + ", " + std::to_string(m.nv) + " controls) ...");
        const OptResult r = optimize(m, c, x0, s.opt);
        out.iters = r.iters;
        out.qConst.assign(size_t(nw), {}); out.qOpt.assign(size_t(nw), {});
        for (int w = 0; w < nw; ++w)
            for (int j = 0; j < nc; ++j) { out.qConst[size_t(w)].push_back(x0[size_t(w * nc + j)]); out.qOpt[size_t(w)].push_back(r.x[size_t(w * nc + j)]); }
        for (size_t p = 0; p < s.patches.size(); ++p) {
            out.predConst.push_back(forecast(m.t, modelStress(m, x0, int(p)), s.patches[p], s.injDays + m.t[0]));
            out.predOpt.push_back(forecast(m.t, modelStress(m, r.x, int(p)), s.patches[p], s.injDays + m.t[0]));
            out.eventsConst += s.patches[p].weight * out.predConst.back().excess;
            out.eventsOpt += s.patches[p].weight * out.predOpt.back().excess;
            out.peakConst = std::max(out.peakConst, out.predConst.back().peakR / std::max(1e-300, s.patches[p].rate));
            out.peakOpt = std::max(out.peakOpt, out.predOpt.back().peakR / std::max(1e-300, s.patches[p].rate));
        }
        {   // C22 (H1f): design potential from the closed form on the governing patch (largest weighted constant-rate excess)
            size_t pg = 0; double best = -1;
            for (size_t p = 0; p < s.patches.size(); ++p) {
                const double e = s.patches[p].weight * out.predConst[p].excess;
                if (e > best) { best = e; pg = p; }
            }
            if (!out.predConst.empty()) {
                const Forecast& F = out.predConst[pg];
                const double tShut = s.injDays + m.t[0];
                double sShut = 0;
                for (size_t i = 0; i < F.t.size(); ++i) if (F.t[i] <= tShut + 1e-9) sShut = F.dS[i];
                double t63 = s.injDays;
                for (size_t i = 0; i < F.t.size(); ++i)
                    if (sShut > 0 && F.dS[i] >= (1.0 - std::exp(-1.0)) * sShut) { t63 = F.t[i] - m.t[0]; break; }
                out.Pi = F.maxDS / std::max(1e-12, s.patches[pg].aSigma);
                out.theta = std::max(1e-4, t63 / std::max(1e-12, s.injDays));
                // The closed form assumes a near-field patch whose stress tracks the injection (first-order response, peak at shut-in).
                // A far-field patch (diffusive delay: stress keeps rising after shut-in) is a different regime - there the optimum
                // front-loads the injection (C22 finding) and the closed-form number would be misleading.
                out.farField = F.maxDS > 1.05 * sShut;
                if (!out.farField) {
                    const GainEstimate ge = scheduleGain(out.Pi, out.theta, std::max(1.0001, s.maxRateFactor));
                    out.gainTheory = ge.gain; out.gainTheoryNoBurst = ge.gainNoBurst; out.gainTheoryIsBound = !ge.exact;
                }
            }
            for (int w = 0; w < nw; ++w) {
                const size_t v = size_t(w * nc + nc - 1);
                if (r.x[v] >= 0.95 * c.ub[v] && r.x[v] > 1.2 * meanRate[size_t(w)]) out.endsWithBurst = true;
            }
        }
        char b[400];
        if (out.farField)
            std::snprintf(b, sizeof b, "Design potential: Pi %.3g, far-field patch (stress peaks after shut-in) - closed form not applicable; "
                          "expect a front-loaded optimum", out.Pi);
        else
            std::snprintf(b, sizeof b, "Design potential (first-order reference, H1f): Pi %.3g, theta %.3g -> %s%.1f %% fewer events than constant rate "
                          "(%.1f %% without a shut-in burst; diffusive responses reach less, H1g)%s", out.Pi, out.theta, out.gainTheoryIsBound ? "at least " : "up to ", 100 * out.gainTheory, 100 * out.gainTheoryNoBurst,
                          out.endsWithBurst ? "; the optimum ends with a max-rate burst at shut-in" : "");
        say(b);
        std::snprintf(b, sizeof b, "Optimum: %d iterations, KKT %.2e; predicted excess events %.4g -> %.4g (%+.1f%%), peak R/r %.3g -> %.3g",
                      r.iters, r.kkt, out.eventsConst, out.eventsOpt, 100 * (out.eventsOpt / (std::fabs(out.eventsConst) > 0 ? out.eventsConst : 1) - 1),
                      out.peakConst, out.peakOpt);
        say(b);
        // 4. nonlinear verification: forward runs of the constant and the optimal schedule
        bool optRunFailed = false;
        for (int which = 0; which < 2; ++which) {
            say(which == 0 ? "Verification run: constant schedule ..." : "Verification run: optimal schedule ...");
            Capture cv;
            try {
                cv = runCapture(designCase(base, s, which == 0 ? out.qConst : out.qOpt), s.patches, cancel); out.runs++;
            } catch (const std::exception& e) {
                if (which == 0) throw;
                // the optimal schedule cannot be simulated (e.g. rates that trigger fracture growth / time-step collapse):
                // fall back to the verified constant schedule instead of failing the whole design
                say(std::string("WARNING: the optimal schedule could not be simulated (") + e.what() +
                    ") - returning the constant schedule; reduce 'Max rate / mean rate'");
                out.qOpt = out.qConst; out.simOpt = out.simConst; out.predOpt = out.predConst;
                out.simEventsOpt = out.simEventsConst; out.eventsOpt = out.eventsConst; out.peakOpt = out.peakConst;
                optRunFailed = true;
                break;
            }
            double ev = 0;
            for (size_t p = 0; p < s.patches.size(); ++p) {
                const double k = s.patches[p].friction * s.patches[p].cfsFactor;
                std::vector<double> d(cv.t.size());
                for (size_t i = 0; i < d.size(); ++i) d[i] = k * (cv.p[p][i] - cv.p[p][0]);
                Forecast f = forecast(cv.t, d, s.patches[p], s.injDays + cv.t[0]);
                ev += s.patches[p].weight * f.excess;
                (which == 0 ? out.simConst : out.simOpt).push_back(std::move(f));
            }
            (which == 0 ? out.simEventsConst : out.simEventsOpt) = ev;
        }
        if (optRunFailed) {
            out.ok = true;
            out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            return out;
        }
        auto relErr = [](double a, double b) { return std::fabs(a - b) / std::max(1e-12, std::fabs(b)); };
        out.linearityError = std::max(relErr(out.simEventsConst, out.eventsConst), relErr(out.simEventsOpt, out.eventsOpt));
        std::snprintf(b, sizeof b, "Simulated (nonlinear) excess events: constant %.4g, optimal %.4g (%+.1f%%); linear-model error %.2f%%",
                      out.simEventsConst, out.simEventsOpt, 100 * (out.simEventsOpt / (std::fabs(out.simEventsConst) > 0 ? out.simEventsConst : 1) - 1),
                      100 * out.linearityError);
        say(b);
        if (out.linearityError > 0.2)
            say("WARNING: the reservoir response is strongly nonlinear for these rates (e.g. growing hydraulic fractures, BHP limits "
                "taking over, two-phase effects). The correction passes and the verification runs decide; a smaller 'Max rate / mean "
                "rate' keeps the design in the near-linear range");
        // 5. sequential convex programming with model correction: dS0 += (simulated - predicted) stress of the last verified
        //    schedule, re-optimise, re-verify. The best VERIFIED schedule is returned (never worse than the first optimum).
        std::vector<double> x = r.x;
        double bestSim = out.simEventsOpt, lastPred = out.eventsOpt;
        out.corrPred.push_back(out.eventsOpt); out.corrSim.push_back(out.simEventsOpt);
        std::vector<Forecast> lastSim = out.simOpt;
        out.finalModelError = relErr(out.simEventsOpt, out.eventsOpt);
        for (int it = 0; it < s.correctionIters && relErr(out.corrSim.back(), lastPred) > s.correctionTol; ++it) {
            if (cancel && cancel->load()) break;
            for (size_t p = 0; p < s.patches.size(); ++p) {
                const std::vector<double> pred = modelStress(m, x, int(p));
                const Forecast& fs = lastSim[p];
                for (size_t i = 0; i < nt; ++i) m.dS0[p][i] += interp(fs.t, fs.dS, m.t[i]) - pred[i];
            }
            const OptResult rc = optimize(m, c, x, s.opt);
            x = rc.x; out.iters += rc.iters;
            std::vector<std::vector<double>> q(static_cast<size_t>(nw));
            for (int w = 0; w < nw; ++w) for (int j = 0; j < nc; ++j) q[size_t(w)].push_back(x[size_t(w * nc + j)]);
            Capture cv;
            try {
                cv = runCapture(designCase(base, s, q), s.patches, cancel); out.runs++; out.corrections++;
            } catch (const std::exception& e) {
                std::snprintf(b, sizeof b, "Correction pass %d: forward run failed (%s) - keeping the best verified schedule", it + 1, e.what());
                say(b);
                break;
            }
            double ev = 0, evp = 0;
            std::vector<Forecast> fsim, fpred;
            for (size_t p = 0; p < s.patches.size(); ++p) {
                const double k = s.patches[p].friction * s.patches[p].cfsFactor;
                std::vector<double> d(cv.t.size());
                for (size_t i = 0; i < d.size(); ++i) d[i] = k * (cv.p[p][i] - cv.p[p][0]);
                fsim.push_back(forecast(cv.t, d, s.patches[p], s.injDays + cv.t[0]));
                fpred.push_back(forecast(m.t, modelStress(m, x, int(p)), s.patches[p], s.injDays + m.t[0]));
                ev += s.patches[p].weight * fsim.back().excess; evp += s.patches[p].weight * fpred.back().excess;
            }
            lastPred = evp; out.corrPred.push_back(evp); out.corrSim.push_back(ev); lastSim = fsim;
            std::snprintf(b, sizeof b, "Correction pass %d: predicted %.4g, simulated %.4g (model error %.2f%%)", it + 1, evp, ev, 100 * relErr(ev, evp));
            say(b);
            if (ev < bestSim) {
                bestSim = ev; out.qOpt = q; out.simOpt = fsim; out.predOpt = fpred; out.simEventsOpt = ev; out.eventsOpt = evp;
                out.finalModelError = relErr(ev, evp);
                out.peakOpt = 0;
                for (size_t p = 0; p < s.patches.size(); ++p) out.peakOpt = std::max(out.peakOpt, fpred[p].peakR / std::max(1e-300, s.patches[p].rate));
            }
        }
        if (out.simEventsOpt > out.simEventsConst) {   // never return a verified schedule that is worse than constant rate
            say("No designed schedule beat constant rate in the verification runs - returning the constant schedule");
            out.qOpt = out.qConst; out.simOpt = out.simConst; out.predOpt = out.predConst;
            out.simEventsOpt = out.simEventsConst; out.eventsOpt = out.eventsConst; out.peakOpt = out.peakConst;
            out.finalModelError = relErr(out.simEventsConst, out.eventsConst);
        }
        std::snprintf(b, sizeof b, "Returned schedule: simulated excess events %.4g vs constant %.4g (%+.1f%%), model error %.2f%%",
                      out.simEventsOpt, out.simEventsConst, 100 * (out.simEventsOpt / (std::fabs(out.simEventsConst) > 0 ? out.simEventsConst : 1) - 1),
                      100 * out.finalModelError);
        say(b);
        out.ok = true;
    } catch (const std::exception& e) { out.ok = false; out.error = e.what(); say(std::string("ERROR: ") + e.what()); }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

// ---------------------------------------------------------------------------------------------------------------------
// C23 research Cycle 5 (F14): stage-count kernel fits
namespace {
struct StageFitter {
    const std::vector<StageObs>& S;
    std::vector<int> tests;
    explicit StageFitter(const std::vector<StageObs>& s) : S(s) {
        for (const auto& o : S) if (std::find(tests.begin(), tests.end(), o.test) == tests.end()) tests.push_back(o.test);
    }
    // Poisson log-likelihood with the per-test MLE factor for the given unnormalised expectations
    double ll(const std::vector<double>& m) const {
        double L = 0;
        for (int t : tests) {
            double sn = 0, sm = 0;
            for (size_t i = 0; i < S.size(); ++i) if (S[i].test == t) { sn += S[i].count; sm += m[i]; }
            const double a = sm > 0 ? sn / sm : 0;
            for (size_t i = 0; i < S.size(); ++i) if (S[i].test == t) {
                const double mu = std::max(1e-300, a * m[i]);
                L += S[i].count * std::log(mu) - mu - std::lgamma(S[i].count + 1.0);
            }
        }
        return L;
    }
};
}  // namespace

KernelFit fitStageKernels(const std::vector<StageObs>& stagesIn) {
    KernelFit F;
    std::vector<StageObs> S = stagesIn;
    if (S.empty()) { F.verdict = "no stages"; return F; }
    {   // cumulative volume in the given order per test
        std::vector<std::pair<int, double>> acc;
        for (auto& o : S) {
            auto it = std::find_if(acc.begin(), acc.end(), [&](const std::pair<int, double>& q) { return q.first == o.test; });
            if (it == acc.end()) { acc.push_back({o.test, 0.0}); it = acc.end() - 1; }
            if (o.cumVolume < 0) o.cumVolume = it->second;
            it->second += std::max(0.0, o.volume);
        }
    }
    StageFitter fit(S);
    const int nt = int(fit.tests.size());
    const size_t n = S.size();
    auto add = [&](const std::string& nm, double L, int extra, double e1 = 0, double e2 = 0) {
        KernelModel m; m.name = nm; m.logL = L; m.params = nt + extra; m.aic = 2 * m.params - 2 * L; m.exponent = e1; m.exponent2 = e2; F.models.push_back(m);
    };
    bool hasE = true, hasQ = true, hasD = true;
    for (const auto& o : S) { hasE = hasE && o.energy > 0; hasQ = hasQ && o.rate > 0; hasD = hasD && o.duration > 0; }
    std::vector<double> V(n), Q(n), C(n), m(n);
    for (size_t i = 0; i < n; ++i) { V[i] = std::max(1e-12, S[i].volume); Q[i] = std::max(1e-12, S[i].rate); C[i] = 1 + std::max(0.0, S[i].cumVolume); }
    for (size_t i = 0; i < n; ++i) m[i] = 1;
    add("constant per stage", fit.ll(m), 0);
    if (hasD) { for (size_t i = 0; i < n; ++i) m[i] = S[i].duration; add("duration", fit.ll(m), 0); }
    add("volume", fit.ll(V), 0);
    if (hasE) { for (size_t i = 0; i < n; ++i) m[i] = S[i].energy; add("hydraulic energy", fit.ll(m), 0); }
    auto llBG = [&](double b, double g) { for (size_t i = 0; i < n; ++i) m[i] = V[i] * std::pow(Q[i], b) * std::pow(C[i], g); return fit.ll(m); };
    // grids (profile likelihood); exponents in [-2, 3]
    const int NB = 101, NG = 101;
    auto bAt = [](int k) { return -2.0 + 5.0 * k / 100.0; };
    double bestG = 0, LG = -1e300;
    for (int k = 0; k < NG; ++k) { const double L = llBG(0, bAt(k)); if (L > LG) { LG = L; bestG = bAt(k); } }
    add("volume x (1 + V_cum)^gamma", LG, 1, bestG);
    if (hasQ) {
        double bestB = 0, LB = -1e300;
        for (int k = 0; k < NB; ++k) { const double L = llBG(bAt(k), 0); if (L > LB) { LB = L; bestB = bAt(k); } }
        add("volume x rate^beta", LB, 1, bestB);
        F.rateExponentAlone = bestB;
        std::vector<double> prof(NB, -1e300);
        double LJ = -1e300, bJ = 0, gJ = 0;
        for (int kb = 0; kb < NB; ++kb)
            for (int kg = 0; kg < NG; ++kg) {
                const double L = llBG(bAt(kb), bAt(kg));
                prof[size_t(kb)] = std::max(prof[size_t(kb)], L);
                if (L > LJ) { LJ = L; bJ = bAt(kb); gJ = bAt(kg); }
            }
        add("volume x rate^beta x (1 + V_cum)^gamma", LJ, 2, bJ, gJ);
        F.rateExponent = bJ; F.cumExponent = gJ;
        F.rateLo = 1e300; F.rateHi = -1e300;
        for (int kb = 0; kb < NB; ++kb) if (prof[size_t(kb)] >= LJ - 1.92) { F.rateLo = std::min(F.rateLo, bAt(kb)); F.rateHi = std::max(F.rateHi, bAt(kb)); }
        double mx = 0, my = 0, sxy = 0, sxx = 0, syy = 0;
        for (size_t i = 0; i < n; ++i) { mx += std::log(Q[i]); my += std::log(C[i]); }
        mx /= n; my /= n;
        for (size_t i = 0; i < n; ++i) { const double a = std::log(Q[i]) - mx, b = std::log(C[i]) - my; sxy += a * b; sxx += a * a; syy += b * b; }
        F.rateCumCorrelation = (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0;
    }
    std::sort(F.models.begin(), F.models.end(), [](const KernelModel& a, const KernelModel& b) { return a.aic < b.aic; });
    F.best = F.models.front().name;
    char buf[400];
    if (!hasQ) F.verdict = "no rates given: kernel ranking only";
    else if (F.rateLo > 0)
        std::snprintf(buf, sizeof buf, "lower rates reduce events at fixed volume (N ~ V q^%.2f, 95%% [%.2f, %.2f]): rate scheduling can help", F.rateExponent, F.rateLo, F.rateHi), F.verdict = buf;
    else if (F.rateHi < 0)
        std::snprintf(buf, sizeof buf, "higher rates give FEWER events per volume (beta %.2f, 95%% [%.2f, %.2f])", F.rateExponent, F.rateLo, F.rateHi), F.verdict = buf;
    else
        std::snprintf(buf, sizeof buf, "no resolvable rate effect at fixed volume (beta %.2f, 95%% [%.2f, %.2f]; rate vs earlier-volume correlation %.2f): "
                      "scheduling cannot be shown to reduce events - volume and site factors dominate", F.rateExponent, F.rateLo, F.rateHi, F.rateCumCorrelation), F.verdict = buf;
    return F;
}

} // namespace seis
} // namespace ressim
