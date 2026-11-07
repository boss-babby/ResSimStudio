// C23 geomechanics G3: quasi-dynamic rate-and-state friction on an embedded fault (see rsf.h for the model).
#include "ressim/rsf.h"
#include "ressim/efem.h"
#include "ressim/geocouple.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>

namespace ressim {
namespace rsf {

namespace {

double expint1(double x) {   // E1(x), x > 0
    if (x <= 0) return 1e300;
    if (x < 1.0) {
        double sum = 0, term = 1;
        for (int k = 1; k < 60; ++k) {
            term *= -x / k;
            sum += -term / k;
            if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
        }
        return -0.5772156649015329 - std::log(x) + sum;
    }
    // continued fraction (Lentz)
    double b = x + 1, c = 1e300, d = 1 / b, h = d;
    for (int i = 1; i < 200; ++i) {
        const double an = -double(i) * i;
        b += 2;
        d = 1 / (an * d + b);
        c = b + an / c;
        const double del = c * d;
        h *= del;
        if (std::fabs(del - 1) < 1e-15) break;
    }
    return h * std::exp(-x);
}

inline double sexp(double x) { return std::exp(std::min(x, 700.0)); }

struct Model {
    const Input& in;
    int n;
    double L, dx, eta, cs;
    std::vector<double> xc, a, sn, tau0, K;   // K row-major n x n
    std::vector<double> sWarm;                // warm start of the frictional stress per element
    std::vector<double> dp, dtau, tau, sig, V;
    explicit Model(const Input& i) : in(i), n(std::max(1, i.n)) {}

    void forcingAt(double t, std::vector<double>& p, std::vector<double>& dt) const {
        p.assign(size_t(n), 0.0); dt.assign(size_t(n), 0.0);
        if (in.inj.on) {
            const double xw = in.inj.at * L;
            auto src = [&](double tt, double sgn) {
                if (tt <= 0) return;
                for (int i = 0; i < n; ++i) {
                    const double dxw = xc[size_t(i)] - xw, r2 = dxw * dxw + in.inj.offset * in.inj.offset;
                    p[size_t(i)] += sgn * in.inj.pStar * expint1(std::max(r2, 1e-6) / (4 * in.inj.diffusivity * tt));
                }
            };
            src(t - in.inj.tStart, 1.0);
            if (t > in.inj.tStop) src(t - in.inj.tStop, -1.0);
        }
        if (in.forceT.size() >= 1 && (in.forceP.size() == in.forceT.size() || in.forceTau.size() == in.forceT.size())) {
            const auto& T = in.forceT;
            size_t k = size_t(std::upper_bound(T.begin(), T.end(), t) - T.begin());
            size_t k0, k1; double w;
            if (k == 0) { k0 = k1 = 0; w = 0; }
            else if (k >= T.size()) { k0 = k1 = T.size() - 1; w = 0; }
            else { k0 = k - 1; k1 = k; w = (t - T[k0]) / std::max(1e-30, T[k1] - T[k0]); }
            auto add = [&](const std::vector<std::vector<double>>& F, std::vector<double>& o) {
                if (F.size() != T.size()) return;
                const auto& A = F[k0]; const auto& B = F[k1];
                for (int i = 0; i < n && size_t(i) < A.size() && size_t(i) < B.size(); ++i) o[size_t(i)] += (1 - w) * A[size_t(i)] + w * B[size_t(i)];
            };
            add(in.forceP, p); add(in.forceTau, dt);
        }
        if (in.forcing) {
            std::vector<double> p2(size_t(n), 0.0), t2(size_t(n), 0.0);
            in.forcing(t, p2, t2);
            for (int i = 0; i < n; ++i) { p[size_t(i)] += i < int(p2.size()) ? p2[size_t(i)] : 0; dt[size_t(i)] += i < int(t2.size()) ? t2[size_t(i)] : 0; }
        }
    }

    // frictional stress s from tau - eta V(s) = s, V(s) = 2 V0 exp(-psi/a) sinh(s / (sig a))
    double solveV(int i, double ta, double sg, double psi) {
        const double ai = a[size_t(i)], V0 = in.V0;
        auto Vof = [&](double s, double& dV) {
            const double x = s / (sg * ai), e = -psi / ai;
            const double ep = sexp(x + e), em = sexp(-x + e);
            dV = V0 * (ep + em) / (sg * ai);
            return V0 * (ep - em);
        };
        double lo = std::min(0.0, ta), hi = std::max(0.0, ta);
        double s = std::clamp(sWarm[size_t(i)], lo, hi);
        for (int it = 0; it < 80; ++it) {
            double dV; const double v = Vof(s, dV);
            const double h = s + eta * v - ta;
            if (h > 0) hi = s; else lo = s;
            const double dh = 1 + eta * dV;
            double sn = s - h / dh;
            if (!(sn > lo && sn < hi)) sn = 0.5 * (lo + hi);
            if (std::fabs(sn - s) <= 1e-13 * std::max(1.0, std::fabs(ta)) + 1e-12) { s = sn; break; }
            s = sn;
            if (hi - lo <= 1e-14 * std::max(1.0, std::fabs(ta))) break;
        }
        sWarm[size_t(i)] = s;
        double dV; return Vof(s, dV);
    }

    // y = [delta(n), psi(n)] -> dy
    void rhs(double t, const std::vector<double>& y, std::vector<double>& f, bool alphaTerm) {
        forcingAt(t, dp, dtau);
        f.assign(size_t(2 * n), 0.0);
        tau.assign(size_t(n), 0.0); sig.assign(size_t(n), 0.0); V.assign(size_t(n), 0.0);
        std::vector<double> dsig;
        if (alphaTerm && in.alpha > 0) {   // d sigma' / dt = -dp/dt (central difference)
            const double h = std::max(1.0, 1e-6 * std::fabs(t));
            std::vector<double> pa, pb, ta_, tb_;
            forcingAt(t + h, pa, ta_); forcingAt(std::max(0.0, t - h), pb, tb_);
            dsig.resize(size_t(n));
            for (int i = 0; i < n; ++i) dsig[size_t(i)] = -(pa[size_t(i)] - pb[size_t(i)]) / (t + h - std::max(0.0, t - h));
        }
        for (int i = 0; i < n; ++i) {
            double ti = tau0[size_t(i)] + in.tauRate * t + dtau[size_t(i)];
            const double* Ki = &K[size_t(i) * size_t(n)];
            for (int j = 0; j < n; ++j) ti += Ki[j] * (y[size_t(j)] - in.Vpl * t);
            const double sg = std::max(1e5, sn[size_t(i)] - in.p0 - dp[size_t(i)]);
            const double psi = y[size_t(n + i)];
            const double v = solveV(i, ti, sg, psi);
            tau[size_t(i)] = ti; sig[size_t(i)] = sg; V[size_t(i)] = v;
            f[size_t(i)] = v;
            double dpsi = in.b * in.V0 / in.Dc * (sexp((in.f0 - psi) / in.b) - std::fabs(v) / in.V0);
            if (!dsig.empty()) dpsi -= in.alpha * dsig[size_t(i)] / sg;
            f[size_t(n + i)] = dpsi;
        }
    }
};

}  // namespace

double nucleationLength(double G, double nu, double a, double b, double Dc, double sig) {
    if (!(b > a) || sig <= 0) return 1e300;
    const double Gp = G / (1 - nu);
    return 2 * Gp * Dc * b / (M_PI * sig * (b - a) * (b - a));
}

bool run(const Input& in, Output& out, std::atomic<bool>* cancel) {
    const auto t0c = std::chrono::steady_clock::now();
    out = Output{};
    Model M(in);
    const int n = M.n;
    M.L = std::hypot(in.x1 - in.x0, in.y1 - in.y0);
    if (!(M.L > 0)) { out.error = "zero fault length"; return false; }
    if (!(in.Dc > 0) || !(in.V0 > 0) || !(in.b >= 0) || !(in.G > 0)) { out.error = "invalid friction / elastic parameters"; return false; }
    M.dx = M.L / n;
    M.cs = std::sqrt(in.G / in.rho);
    M.eta = in.G / (2 * M.cs);
    M.xc.resize(size_t(n)); M.a.resize(size_t(n)); M.sn.resize(size_t(n)); M.tau0.resize(size_t(n)); M.sWarm.assign(size_t(n), 0.0);
    for (int i = 0; i < n; ++i) {
        M.xc[size_t(i)] = (i + 0.5) * M.dx;
        const double fr = M.xc[size_t(i)] / M.L;
        M.a[size_t(i)] = (int(in.aProf.size()) == n) ? in.aProf[size_t(i)] : ((fr >= in.vwFrom && fr <= in.vwTo) ? in.aVW : in.aVS);
        M.sn[size_t(i)] = (int(in.sigmaNProf.size()) == n) ? in.sigmaNProf[size_t(i)] : in.sigmaN;
        if (!(M.a[size_t(i)] > 0)) { out.error = "a must be > 0"; return false; }
    }
    // kernel
    const double c = in.G / (2 * M_PI * (1 - in.nu));
    M.K.resize(size_t(n) * size_t(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double xp = (j + 1) * M.dx, xm = j * M.dx, x = M.xc[size_t(i)];
            M.K[size_t(i) * size_t(n) + size_t(j)] = c * (1 / (x - xp) - 1 / (x - xm));
        }
    // initial state: steady state at vInit; tau0 from the ratio, the profile, or the steady state
    std::vector<double> y(size_t(2 * n));
    std::vector<double> p0f, t0f;
    M.forcingAt(0.0, p0f, t0f);
    for (int i = 0; i < n; ++i) {
        const double sg = std::max(1e5, M.sn[size_t(i)] - in.p0 - p0f[size_t(i)]);
        const double ai = M.a[size_t(i)];
        const double psi = in.f0 - in.b * std::log(in.vInit / in.V0);
        y[size_t(i)] = 0; y[size_t(n + i)] = psi;
        double t0v;
        if (int(in.tau0Prof.size()) == n) t0v = in.tau0Prof[size_t(i)];
        else if (in.tauRatio >= 0) t0v = in.tauRatio * sg;
        else t0v = sg * ai * std::asinh(in.vInit / (2 * in.V0) * std::exp(psi / ai)) + M.eta * in.vInit;   // exact balance at vInit
        M.tau0[size_t(i)] = t0v - t0f[size_t(i)];   // so that tau(0) = t0v including a non-zero initial stress forcing
        M.sWarm[size_t(i)] = t0v;
    }
    // diagnostics
    {
        double sgm = 0, bmA = 0; int nv = 0;
        for (int i = 0; i < n; ++i) { sgm += std::max(1e5, M.sn[size_t(i)] - in.p0); if (in.b > M.a[size_t(i)]) { bmA += in.b - M.a[size_t(i)]; ++nv; } }
        sgm /= n;
        const double aw = nv ? in.b - bmA / nv : in.aVW;
        out.hStar = nucleationLength(in.G, in.nu, aw, in.b, in.Dc, sgm);
        out.Lb = in.G / (1 - in.nu) * in.Dc / (std::max(1e-12, in.b) * sgm);
        out.kc = std::max(0.0, in.b - aw) * sgm / in.Dc;
        out.kSelf = -M.K[0];
        out.dx = M.dx;
    }
    out.s = M.xc; out.a = M.a;

    // Dormand-Prince 5(4)
    static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static const double a21 = 1. / 5;
    static const double a31 = 3. / 40, a32 = 9. / 40;
    static const double a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9;
    static const double a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729;
    static const double a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656;
    static const double b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84;
    static const double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
    const size_t N2 = size_t(2 * n);
    std::vector<double> k1, k2, k3, k4, k5, k6, k7, yt(N2), yn(N2);
    const bool al = in.alpha > 0;
    double t = 0;
    M.rhs(t, y, k1, al);
    double dt = std::min(in.dtMax, std::max(1e-3, 0.01 * in.Dc / std::max(1e-15, *std::max_element(M.V.begin(), M.V.end()))));
    bool inEv = false;
    Event ev; std::vector<double> dStart, tauStart;
    double lastSnapT = -1e300; long lastSnapStep = -1000000;
    const double snapDt = in.tEnd / std::max(1, in.snapshots);
    auto snap = [&](double tt, const std::vector<double>& yy) {
        if (out.snapT.size() >= 6000) return;
        std::vector<double> lv(static_cast<size_t>(n)), sl(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) { lv[size_t(i)] = std::log10(std::max(1e-30, std::fabs(M.V[size_t(i)]))); sl[size_t(i)] = yy[size_t(i)]; }
        out.snapT.push_back(tt); out.snapLogV.push_back(std::move(lv)); out.snapSlip.push_back(std::move(sl));
    };
    snap(0, y); lastSnapT = 0;
    long recEvery = 1;
    auto vmaxOf = [&](int& im) { double vm = 0; im = 0; for (int i = 0; i < n; ++i) if (std::fabs(M.V[size_t(i)]) > vm) { vm = std::fabs(M.V[size_t(i)]); im = i; } return vm; };
    while (t < in.tEnd) {
        if (cancel && cancel->load()) { out.error = "cancelled"; break; }
        if (out.steps >= in.maxSteps) { out.error = "step limit reached"; break; }
        dt = std::min(dt, in.tEnd - t);
        auto stage = [&](std::vector<double>& kk, double cc, std::initializer_list<std::pair<double, const std::vector<double>*>> terms) {
            for (size_t q = 0; q < N2; ++q) { double s = y[q]; for (const auto& tm : terms) s += dt * tm.first * (*tm.second)[q]; yt[q] = s; }
            M.rhs(t + cc * dt, yt, kk, al);
        };
        stage(k2, c2, {{a21, &k1}});
        stage(k3, c3, {{a31, &k1}, {a32, &k2}});
        stage(k4, c4, {{a41, &k1}, {a42, &k2}, {a43, &k3}});
        stage(k5, c5, {{a51, &k1}, {a52, &k2}, {a53, &k3}, {a54, &k4}});
        stage(k6, 1.0, {{a61, &k1}, {a62, &k2}, {a63, &k3}, {a64, &k4}, {a65, &k5}});
        for (size_t q = 0; q < N2; ++q) yn[q] = y[q] + dt * (b1 * k1[q] + b3 * k3[q] + b4 * k4[q] + b5 * k5[q] + b6 * k6[q]);
        M.rhs(t + dt, yn, k7, al);
        double err = 0;
        bool bad = false;
        for (int i = 0; i < 2 * n; ++i) {
            const size_t q = size_t(i);
            const double e = dt * (e1 * k1[q] + e3 * k3[q] + e4 * k4[q] + e5 * k5[q] + e6 * k6[q] + e7 * k7[q]);
            const double sc = i < n ? in.rtol * std::max(std::fabs(yn[q]), in.Dc) : in.rtol * std::max(std::fabs(yn[q]), in.b > 0 ? in.b : 0.01);
            const double r = std::fabs(e) / sc;
            if (!std::isfinite(r)) bad = true;
            err = std::max(err, r);
        }
        if (bad || err > 1.0) {
            ++out.rejected;
            dt *= bad ? 0.1 : std::max(0.1, 0.9 * std::pow(err, -0.2));
            if (dt < 1e-9) { out.error = "time step underflow"; break; }
            M.rhs(t, y, k1, al);   // restore M.V / tau for the current state
            continue;
        }
        // accept
        t += dt; y.swap(yn); k1.swap(k7);   // FSAL: k7 evaluated at the new state, M.V/tau/sig are from it
        ++out.steps;
        int im; const double vm = vmaxOf(im);
        double pm = 0; for (double v : M.dp) pm = std::max(pm, v);
        if (out.steps % recEvery == 0) {
            out.t.push_back(t); out.vmax.push_back(vm); out.pMax.push_back(pm);
            if (out.t.size() > 20000) {   // decimate
                std::vector<double> a2, b2, c2v;
                for (size_t q = 0; q < out.t.size(); q += 2) { a2.push_back(out.t[q]); b2.push_back(out.vmax[q]); c2v.push_back(out.pMax[q]); }
                out.t.swap(a2); out.vmax.swap(b2); out.pMax.swap(c2v); recEvery *= 2;
            }
        }
        // events
        if (!inEv && vm > in.vEvent) {
            inEv = true; ev = Event{}; ev.tStart = t; ev.hypo = im; ev.sHypo = M.xc[size_t(im)]; ev.vMax = vm;
            dStart.assign(y.begin(), y.begin() + n); tauStart = M.tau;
        } else if (inEv) {
            ev.vMax = std::max(ev.vMax, vm);
            if (vm < in.vEvent) {
                inEv = false; ev.tEnd = t;
                double smax = 0; for (int i = 0; i < n; ++i) smax = std::max(smax, y[size_t(i)] - dStart[size_t(i)]);
                double sum = 0, dsum = 0; int nr = 0;
                for (int i = 0; i < n; ++i) {
                    const double sl = y[size_t(i)] - dStart[size_t(i)];
                    if (sl > 0.01 * smax) { sum += sl * M.dx; dsum += tauStart[size_t(i)] - M.tau[size_t(i)]; ++nr; }
                }
                ev.slipMax = smax; ev.length = nr * M.dx; ev.slipMean = nr ? sum / ev.length : 0; ev.stressDrop = nr ? dsum / nr : 0;
                ev.moment = in.G * sum * ev.length;
                ev.Mw = ev.moment > 0 ? (2.0 / 3.0) * (std::log10(ev.moment) - 9.1) : -10;
                out.events.push_back(ev);
            }
        }
        if (t - lastSnapT >= snapDt || (inEv && out.steps - lastSnapStep >= 4)) { snap(t, y); lastSnapT = t; lastSnapStep = out.steps; }
        // next step
        const double fac = err > 0 ? std::min(5.0, 0.9 * std::pow(err, -0.2)) : 5.0;
        dt = std::min(in.dtMax, dt * fac);
    }
    snap(t, y);
    out.finalSlip.assign(y.begin(), y.begin() + n); out.finalV = M.V;
    out.ok = out.error.empty();
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0c).count();
    return out.ok;
}

bool fromEfem(const efem::Input& efin, const efem::Result& r, int frac, Input& in, std::string& err, double capFrac) {
    if (!r.ok) { err = "the contact solve has no result"; return false; }
    if (frac < 0 || frac >= int(efin.fracs.size())) { err = "no such fracture"; return false; }
    const auto& F = efin.fracs[size_t(frac)];
    std::vector<std::pair<double, std::pair<double, double>>> pts;   // s -> (sigma'_n, |tau|)
    for (const auto& c : r.cut) if (c.frac == frac) pts.push_back({c.s, {-c.tn, std::fabs(c.tt)}});
    if (pts.empty()) { err = "the fracture is not cut by the mesh"; return false; }
    std::sort(pts.begin(), pts.end());
    in.x0 = F.x0; in.y0 = F.y0; in.x1 = F.x1; in.y1 = F.y1;
    const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
    const int n = std::max(1, in.n);
    in.sigmaNProf.assign(size_t(n), 0); in.tau0Prof.assign(size_t(n), 0);
    in.p0 = F.pressure;
    for (int i = 0; i < n; ++i) {
        const double s = (i + 0.5) * L / n - 0.5 * L;
        double sg, ta;
        if (s <= pts.front().first) { sg = pts.front().second.first; ta = pts.front().second.second; }
        else if (s >= pts.back().first) { sg = pts.back().second.first; ta = pts.back().second.second; }
        else {
            size_t k = 1; while (k < pts.size() && pts[k].first < s) ++k;
            const double w = (s - pts[k - 1].first) / std::max(1e-30, pts[k].first - pts[k - 1].first);
            sg = (1 - w) * pts[k - 1].second.first + w * pts[k].second.first;
            ta = (1 - w) * pts[k - 1].second.second + w * pts[k].second.second;
        }
        sg = std::max(sg, 1e5);
        // cap at a fraction of the steady-state strength at vInit (rsf a at this position)
        const double fr = (i + 0.5) / n;
        const double ai = (int(in.aProf.size()) == n) ? in.aProf[size_t(i)] : ((fr >= in.vwFrom && fr <= in.vwTo) ? in.aVW : in.aVS);
        const double fss = in.f0 + (ai - in.b) * std::log(in.vInit / in.V0);
        ta = std::min(ta, capFrac * fss * sg);
        in.sigmaNProf[size_t(i)] = sg + in.p0;   // total = effective + fracture pressure
        in.tau0Prof[size_t(i)] = ta;
    }
    return true;
}

bool fromCoupled(const geocouple::Output& co, int frac, Input& in, std::string& err) {
    if (!co.ok || co.steps.empty()) { err = "the coupled run has no steps"; return false; }
    std::vector<size_t> idx;
    for (size_t k = 0; k < co.frac.size(); ++k) if (co.frac[k] == frac) idx.push_back(k);
    if (idx.empty()) { err = "no pieces on that fracture"; return false; }
    // distinct s positions (layers averaged)
    std::map<double, std::vector<size_t>> byS;
    for (size_t k : idx) byS[std::round(co.s[k] * 1e6) / 1e6].push_back(k);
    std::vector<double> S; for (auto& kv : byS) S.push_back(kv.first);
    const double L = std::hypot(in.x1 - in.x0, in.y1 - in.y0);
    const int n = std::max(1, in.n);
    in.forceT.clear(); in.forceP.clear(); in.forceTau.clear();
    std::vector<double> base;
    for (size_t st = 0; st < co.steps.size(); ++st) {
        const auto& P = co.steps[st].p;
        std::vector<double> ps;
        for (auto& kv : byS) { double s = 0; int c = 0; for (size_t k : kv.second) if (k < P.size()) { s += P[k]; ++c; } ps.push_back(c ? s / c : 0); }
        if (st == 0) base = ps;
        std::vector<double> row(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            const double s = (i + 0.5) * L / n - 0.5 * L;
            double v;
            if (S.size() == 1 || s <= S.front()) v = ps.front() - base.front();
            else if (s >= S.back()) v = ps.back() - base.back();
            else {
                size_t k = 1; while (k < S.size() && S[k] < s) ++k;
                const double w = (s - S[k - 1]) / std::max(1e-30, S[k] - S[k - 1]);
                v = (1 - w) * (ps[k - 1] - base[k - 1]) + w * (ps[k] - base[k]);
            }
            row[size_t(i)] = v;
        }
        in.forceT.push_back(co.steps[st].t * 86400.0);
        in.forceP.push_back(std::move(row));
    }
    return true;
}

}  // namespace rsf
}  // namespace ressim
