// C15: fault / fracture stability, Barton-Bandis aperture and rate-and-state slip (see faultstab.h)
#include "ressim/faultstab.h"
#include <algorithm>
#include <cmath>

namespace ressim {
namespace stab {

namespace {
const double kPi = 3.14159265358979323846;
double rad(double d) { return d * kPi / 180.0; }
void stressTensor(const StressState& s, double S[3][3]) {
    const double th = rad(s.azimuthDeg);
    const double a[3] = {std::sin(th), std::cos(th), 0}, b[3] = {std::cos(th), -std::sin(th), 0}, z[3] = {0, 0, 1};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) S[i][j] = s.SH * a[i] * a[j] + s.Sh * b[i] * b[j] + s.Sv * z[i] * z[j];
}
void planeVectors(const Plane& pl, double n[3], double sv[3], double dd[3]) {
    const double ph = rad(pl.strikeDeg), de = rad(pl.dipDeg);
    sv[0] = std::sin(ph); sv[1] = std::cos(ph); sv[2] = 0;                        // strike
    const double dh[3] = {std::cos(ph), -std::sin(ph), 0};                         // dip direction (strike + 90)
    for (int i = 0; i < 3; ++i) dd[i] = std::cos(de) * dh[i] + std::sin(de) * (i == 2 ? 1.0 : 0.0);   // down-dip
    n[0] = sv[1] * dd[2] - sv[2] * dd[1]; n[1] = sv[2] * dd[0] - sv[0] * dd[2]; n[2] = sv[0] * dd[1] - sv[1] * dd[0];
    const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (int i = 0; i < 3; ++i) n[i] /= l;
}
}  // namespace

StressState fromGradients(double depth, double svGrad, double sHGrad, double shGrad, double azimuthDeg, double p, double biot) {
    StressState s;
    s.Sv = svGrad * depth; s.SH = sHGrad * depth; s.Sh = shGrad * depth;
    s.azimuthDeg = azimuthDeg; s.p = p; s.biot = biot;
    return s;
}

void stressPath(StressState& s, double dp, double nu, double dT, double E, double alphaT) {
    const double dS = stressPathCoeff(s.biot, nu) * dp + (E > 0 && alphaT != 0 ? E * alphaT * dT / (1 - nu) : 0.0);
    s.SH += dS; s.Sh += dS; s.p += dp;
}

Resolved resolve(const StressState& s, const Plane& pl, const Strength& st) {
    Resolved r;
    double S[3][3], n[3], sv[3], dd[3];
    stressTensor(s, S);
    planeVectors(pl, n, sv, dd);
    double t[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) t[i] += S[i][j] * n[j];
    r.sn = t[0] * n[0] + t[1] * n[1] + t[2] * n[2];
    double ts[3];
    for (int i = 0; i < 3; ++i) ts[i] = t[i] - r.sn * n[i];
    r.tau = std::sqrt(ts[0] * ts[0] + ts[1] * ts[1] + ts[2] * ts[2]);
    r.snEff = r.sn - s.biot * s.p;
    r.slipTend = r.snEff > 0 ? r.tau / r.snEff : 1e9;
    const double s1 = std::max({s.Sv, s.SH, s.Sh}), s3 = std::min({s.Sv, s.SH, s.Sh});
    r.dilTend = s1 - s3 > 0 ? (s1 - r.sn) / (s1 - s3) : 0.0;
    r.cff = r.tau - st.mu * r.snEff - st.cohesion;
    r.pCrit = st.mu > 0 && s.biot > 0 ? (r.sn - (r.tau - st.cohesion) / st.mu) / s.biot : 1e30;
    r.dpCrit = r.pCrit - s.p;
    // shear traction acts on the hanging wall: rake from the strike direction (positive = up-dip component reversed sign)
    if (r.tau > 0) {
        const double cs = (ts[0] * sv[0] + ts[1] * sv[1] + ts[2] * sv[2]) / r.tau, cd = (ts[0] * dd[0] + ts[1] * dd[1] + ts[2] * dd[2]) / r.tau;
        r.rakeDeg = std::atan2(-cd, cs) * 180.0 / kPi;
    }
    for (int i = 0; i < 3; ++i) r.n[i] = n[i];
    return r;
}

std::vector<double> slipTendencyGrid(const StressState& s, double dStrike, double dDip, int& nStrike, int& nDip) {
    nStrike = std::max(2, int(std::lround(360.0 / dStrike)) + 1);
    nDip = std::max(2, int(std::lround(90.0 / dDip)) + 1);
    std::vector<double> out(size_t(nStrike) * size_t(nDip), 0.0);
    Strength st;
    for (int j = 0; j < nDip; ++j)
        for (int i = 0; i < nStrike; ++i) {
            Plane pl; pl.strikeDeg = i * 360.0 / (nStrike - 1); pl.dipDeg = j * 90.0 / (nDip - 1);
            out[size_t(j) * size_t(nStrike) + size_t(i)] = std::min(5.0, resolve(s, pl, st).slipTend);
        }
    return out;
}

double bbAperture(const BBParams& b, double snEff, double slip) {
    const double s = std::max(0.0, snEff);
    const double closure = (b.Kn0 > 0 && b.Vm > 0) ? s * b.Vm / (b.Kn0 * b.Vm + s) : 0.0;   // -> Vm as s -> inf
    const double dil = std::min(b.dilMax, std::max(0.0, slip) * std::tan(rad(b.psiDeg)));
    return std::max(1e-7, b.w0 + b.Vm - closure + dil);   // w0: residual aperture at full closure
}

namespace {
double interp(const std::vector<double>& t, const std::vector<double>& v, double x) {
    if (t.empty() || v.size() != t.size()) return 0.0;
    if (x <= t.front()) return v.front();
    if (x >= t.back()) return v.back();
    const size_t k = size_t(std::upper_bound(t.begin(), t.end(), x) - t.begin());
    const double f = (x - t[k - 1]) / std::max(1e-300, t[k] - t[k - 1]);
    return v[k - 1] + f * (v[k] - v[k - 1]);
}
}  // namespace

bool simulateRsf(const RsfInput& in, RsfResult& out) {
    out = RsfResult();
    if (!(in.a > 0) || !(in.b > 0) || !(in.Dc > 0) || !(in.V0 > 0) || !(in.G > 0) || !(in.radius > 0) || !(in.sn0 > 0)) {
        out.error = "rate-and-state: a, b, Dc, V0, G, radius and the initial effective normal stress must be > 0";
        return false;
    }
    const double eta = in.G / (2 * std::max(1.0, in.vs));
    const double k = 7 * kPi * in.G / (16 * in.radius);
    out.kPatch = k; out.kCrit = in.sn0 * (in.b - in.a) / in.Dc;
    const double tEnd = in.tEnd > 0 ? in.tEnd : (in.t.empty() ? 0.0 : in.t.back());
    if (!(tEnd > 0)) { out.error = "rate-and-state: empty loading history"; return false; }
    // initial state: locked (V = Vinit), friction equal to the applied stress -> state variable
    const double Vi = std::max(1e-25, in.Vinit);
    const double fInit = in.tau0 / in.sn0;
    const double lnTh0 = std::log(in.Dc / in.V0) + (fInit - in.f0 - in.a * std::log(Vi / in.V0)) / in.b;
    // state y = (slip, ln theta); V from the force balance (monotone in ln V -> safeguarded Newton)
    double lnVprev = std::log(Vi);
    auto solveV = [&](double t, double slip, double lnTh, double& lnV) -> bool {
        const double sn = std::max(1e3, in.sn0 + interp(in.t, in.dSn, t));
        const double tauA = in.tau0 + interp(in.t, in.dTau, t) - k * slip;
        const double c0 = in.f0 + in.b * (std::log(in.V0 / in.Dc) + lnTh);
        double x = lnV, lo = -80, hi = 20;   // ln V bounds (V up to ~5e8 m/s is never reached: damping caps it)
        for (int it = 0; it < 100; ++it) {
            const double V = std::exp(x);
            const double g = sn * (c0 + in.a * (x - std::log(in.V0))) + eta * V - tauA;
            const double dg = sn * in.a + eta * V;
            if (g > 0) hi = x; else lo = x;
            double xn = x - g / dg;
            if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
            if (std::fabs(xn - x) < 1e-12) { x = xn; lnV = x; return true; }
            x = xn;
        }
        lnV = x;
        return std::isfinite(x);
    };
    auto rhs = [&](double t, const double y[2], double dy[2], double& lnV) -> bool {
        lnV = lnVprev;
        if (!solveV(t, y[0], y[1], lnV)) return false;
        const double V = std::exp(lnV);
        dy[0] = V;
        dy[1] = (1.0 - V * std::exp(y[1]) / in.Dc) / std::exp(y[1]);   // d ln theta / dt = (1 - V theta / Dc) / theta
        return true;
    };
    double t = 0, y[2] = {0.0, lnTh0};
    double h = std::min(tEnd * 1e-6, 1.0);
    auto record = [&](double tt, const double yy[2], double lnV) {
        const double sn = std::max(1e3, in.sn0 + interp(in.t, in.dSn, tt));
        const double V = std::exp(lnV);
        const double tau = in.tau0 + interp(in.t, in.dTau, tt) - k * yy[0] - eta * V;
        RsfSample s{tt, yy[0], V, tau, sn, std::exp(yy[1]), tau / sn, bbAperture(in.bb, sn, yy[0])};
        out.samples.push_back(s);
    };
    double lnV0 = std::log(Vi);
    record(0, y, lnV0);
    bool inEvent = false;
    RsfEvent ev;
    double lastRec = 0;
    // Dormand-Prince 5(4) with error control on (slip / Dc, ln theta)
    static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9;
    static const double a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729;
    static const double a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656;
    static const double b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84;
    static const double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
    while (t < tEnd && out.steps < in.maxSteps) {
        if (t + h > tEnd) h = tEnd - t;
        double k1[2], k2[2], k3[2], k4[2], k5[2], k6[2], k7[2], yt[2], y5[2], lv;
        bool ok = rhs(t, y, k1, lv);
        for (int i = 0; i < 2 && ok; ++i) yt[i] = y[i] + h * a21 * k1[i];
        ok = ok && rhs(t + c2 * h, yt, k2, lv);
        for (int i = 0; i < 2 && ok; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        ok = ok && rhs(t + c3 * h, yt, k3, lv);
        for (int i = 0; i < 2 && ok; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        ok = ok && rhs(t + c4 * h, yt, k4, lv);
        for (int i = 0; i < 2 && ok; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        ok = ok && rhs(t + c5 * h, yt, k5, lv);
        for (int i = 0; i < 2 && ok; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        ok = ok && rhs(t + h, yt, k6, lv);
        for (int i = 0; i < 2 && ok; ++i) y5[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        double lnVnew = lv;
        ok = ok && rhs(t + h, y5, k7, lnVnew);
        double err = 1e30;
        if (ok) {
            const double sc[2] = {in.Dc, 1.0};
            err = 0;
            for (int i = 0; i < 2; ++i) {
                const double e = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
                err = std::max(err, std::fabs(e) / (1e-6 * sc[i] + 1e-6 * std::fabs(y5[i]) * (i == 1 ? 0.0 : 1.0)));
            }
        }
        if (ok && err <= 1.0) {
            t += h; y[0] = y5[0]; y[1] = y5[1]; lnVprev = lnVnew; ++out.steps;
            const double V = std::exp(lnVnew);
            if (!inEvent && V > in.Vseis) { inEvent = true; ev = RsfEvent(); ev.tStart = t; ev.slip = -y[0]; }
            if (inEvent) {
                ev.Vmax = std::max(ev.Vmax, V);
                if (V < 0.1 * in.Vseis) {
                    inEvent = false; ev.tEnd = t; ev.slip += y[0];
                    ev.M0 = in.G * kPi * in.radius * in.radius * ev.slip;
                    ev.Mw = ev.M0 > 0 ? (2.0 / 3.0) * (std::log10(ev.M0) - 9.1) : -99;
                    out.events.push_back(ev);
                }
            }
            // keep the output compact: log-spaced in time plus every step at seismic rates
            if (V > 1e-6 || t - lastRec > 1e-3 * tEnd || t >= tEnd) { record(t, y, lnVnew); lastRec = t; }
            const double fac = err > 0 ? 0.9 * std::pow(err, -0.2) : 5.0;
            h *= std::min(5.0, std::max(0.2, fac));
        } else {
            h *= ok ? std::max(0.1, 0.9 * std::pow(err, -0.25)) : 0.25;
            if (h < 1e-12) { out.error = "rate-and-state: step size underflow at t = " + std::to_string(t) + " s"; return false; }
        }
        h = std::min(h, std::max(1e-9, 0.02 * tEnd));
    }
    if (inEvent) { ev.tEnd = t; ev.slip += y[0]; ev.M0 = in.G * kPi * in.radius * in.radius * ev.slip; ev.Mw = ev.M0 > 0 ? (2.0 / 3.0) * (std::log10(ev.M0) - 9.1) : -99; out.events.push_back(ev); }
    if (out.steps >= in.maxSteps) out.error = "rate-and-state: step limit reached at t = " + std::to_string(t) + " s (results up to there)";
    return true;
}

} // namespace stab
} // namespace ressim
