// C14: hydraulic fracture propagation (PKN / KGD), see hydrofrac.h
#include "ressim/hydrofrac.h"
#include "ressim/grid.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace ressim {
namespace hf {

namespace {
// Gaussian elimination with partial pivoting (A row-major n x n, overwritten); false if singular
bool solveDense(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int k = 0; k < n; ++k) {
        int piv = k; double best = std::fabs(A[size_t(k) * n + k]);
        for (int i = k + 1; i < n; ++i) { const double v = std::fabs(A[size_t(i) * n + k]); if (v > best) { best = v; piv = i; } }
        if (!(best > 0)) return false;
        if (piv != k) { for (int j = 0; j < n; ++j) std::swap(A[size_t(k) * n + j], A[size_t(piv) * n + j]); std::swap(b[size_t(k)], b[size_t(piv)]); }
        const double d = A[size_t(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double f = A[size_t(i) * n + k] / d;
            if (f == 0) continue;
            for (int j = k; j < n; ++j) A[size_t(i) * n + j] -= f * A[size_t(k) * n + j];
            b[size_t(i)] -= f * b[size_t(k)];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[size_t(i)];
        for (int j = i + 1; j < n; ++j) s -= A[size_t(i) * n + j] * b[size_t(j)];
        b[size_t(i)] = s / A[size_t(i) * n + i];
    }
    return true;
}
// tridiagonal (lower l, diagonal d, upper u; rhs r -> solution)
bool solveTri(std::vector<double> l, std::vector<double> d, std::vector<double> u, std::vector<double>& r) {
    const size_t n = d.size();
    for (size_t i = 1; i < n; ++i) {
        if (d[i - 1] == 0) return false;
        const double f = l[i] / d[i - 1];
        d[i] -= f * u[i - 1]; r[i] -= f * r[i - 1];
    }
    if (d[n - 1] == 0) return false;
    r[n - 1] /= d[n - 1];
    for (size_t i = n - 1; i-- > 0;) r[i] = (r[i] - u[i] * r[i + 1]) / d[i];
    return true;
}
}  // namespace

double Result::meanWidth() const {
    return snaps.empty() || snaps.back().L <= 0 || height <= 0 ? 0.0 : volFracture / (2 * snaps.back().L * height);
}

double nordgrenLength(const Input& in, double t) {
    const double Ep = in.E / (1 - in.nu * in.nu), q = 0.5 * in.rate;   // per wing
    return 0.68 * std::pow(Ep * q * q * q / (2 * in.mu * std::pow(in.height, 4)), 0.2) * std::pow(t, 0.8);
}
double nordgrenWidth(const Input& in, double t) {
    const double Ep = in.E / (1 - in.nu * in.nu), q = 0.5 * in.rate;
    return 2.5 * std::pow(2 * in.mu * q * q / (Ep * in.height), 0.2) * std::pow(t, 0.2);
}
double kgdLength(const Input& in, double t) {
    const double Ep = in.E / (1 - in.nu * in.nu), Q = in.rate / in.height;   // total, per unit height
    return 0.6152 * std::pow(Ep * Q * Q * Q / (12 * in.mu), 1.0 / 6.0) * std::pow(t, 2.0 / 3.0);
}
double kgdWidth(const Input& in, double t) {
    const double Ep = in.E / (1 - in.nu * in.nu), Q = in.rate / in.height;
    return std::sqrt(2.1 * std::sqrt(12 * in.mu * Q / Ep) * kgdLength(in, t));
}

// ------------------------------------------------------------------ C20 P3D: equilibrium height (symmetric barriers)
// A(H) = (pi / (2 E')) [K sqrt(2/pi) H^1.5 + (2 dsigma / pi) h sqrt(H^2 - h^2)]  (no cancellation near H = h)
double p3dArea(double H, double h, double ds, double kic, double Ep) {
    if (H <= h) H = h;
    return M_PI / (2 * Ep) * (kic * std::sqrt(2 / M_PI) * std::pow(H, 1.5) + (2 * ds / M_PI) * h * std::sqrt(std::max(0.0, H * H - h * h)));
}
namespace {
double p3dAreaDH(double H, double h, double ds, double kic, double Ep) {
    const double s = std::sqrt(std::max(H * H - h * h, 1e-300));
    return M_PI / (2 * Ep) * (kic * std::sqrt(2 / M_PI) * 1.5 * std::sqrt(H) + (2 * ds / M_PI) * h * H / s);
}
// phi - sin(2 phi) / 2 with phi = acos(h / H) (barrier term of the area relation), stable for phi -> 0
double barrierTerm(double H, double h) {
    const double phi = std::atan2(std::sqrt(std::max(0.0, H * H - h * h)), h);
    if (phi < 1e-3) { const double p2 = phi * phi; return phi * p2 * (2.0 / 3.0 - p2 * 2.0 / 15.0); }
    return phi - 0.5 * std::sin(2 * phi);
}
}
P3DSection p3dSection(double A, double h, double ds, double kic, double Ep, double Hmax) {
    P3DSection S;
    A = std::max(A, 0.0);
    const double pc = kic / std::sqrt(0.5 * M_PI * h);              // contained until the tip K_I reaches K_IC
    const double Ac = M_PI * h * h / (2 * Ep) * pc;
    const double kC = 2 * Ep / (M_PI * h * h);
    if (A <= Ac || (ds <= 0 && kic <= 0)) { S.H = h; S.pnet = kC * A; S.dpdA = kC; return S; }
    if (A >= p3dArea(Hmax, h, ds, kic, Ep)) {                       // height cap: fixed H, pressure from the area relation
        S.H = Hmax; S.capped = true;
        S.dpdA = 2 * Ep / (M_PI * Hmax * Hmax);
        S.pnet = S.dpdA * A + (2 * ds / M_PI) * barrierTerm(Hmax, h);
        return S;
    }
    // A(H) increases monotonically: bracketed Newton with the analytic slope
    double lo = h, hi = Hmax;
    double H = h * std::pow(std::max(A / std::max(Ac, 1e-300), 1.0), 2.0 / 3.0);
    if (!(H > lo && H < hi)) H = 0.5 * (lo + hi);
    for (int it = 0; it < 200; ++it) {
        const double f = p3dArea(H, h, ds, kic, Ep) - A;
        if (f > 0) hi = H; else lo = H;
        if (std::fabs(f) <= 1e-14 * A || hi - lo <= 1e-15 * H) break;
        double Hn = H - f / p3dAreaDH(H, h, ds, kic, Ep);
        if (!(Hn > lo && Hn < hi)) Hn = 0.5 * (lo + hi);
        H = Hn;
    }
    S.H = H;
    // p_net from the area relation (insensitive to round-off in H, unlike K / sqrt(pi H / 2) + (2 dsigma / pi) acos(h / H))
    S.pnet = 2 * Ep * A / (M_PI * H * H) + (2 * ds / M_PI) * barrierTerm(H, h);
    const double sq = std::sqrt(std::max(H * H - h * h, 1e-300));
    const double dpdH = -0.5 * kic * std::sqrt(2 / M_PI) * std::pow(H, -1.5) + (2 * ds / M_PI) * h / (H * sq);
    const double dAdH = p3dAreaDH(H, h, ds, kic, Ep);
    S.dHdA = 1 / dAdH;
    S.dpdA = dpdH / dAdH;
    return S;
}

static bool simulateP3D(const Input& in, Result& out) {
    auto fail = [&](const std::string& m) { out.error = m; return false; };
    if (!(in.stressContrast >= 0)) return fail("stress contrast must be >= 0");
    const int N = in.cells;
    const double dx = in.maxLength / N, Ep = in.E / (1 - in.nu * in.nu), h = in.height, ds = in.stressContrast, K = in.kic;
    const double Hmax = in.maxHeight > h ? in.maxHeight : 20 * h;
    const double q = 0.5 * in.rate;                                 // per wing, m3/s
    out.height = h;
    std::vector<double> A(size_t(N), 0.0), An, H(size_t(N), h), P(size_t(N), 0.0), dP(size_t(N), 0.0), kA(size_t(N), 0.0), dkA(size_t(N), 0.0);
    std::vector<double> tau(size_t(N), -1.0), Lv(size_t(N), 0.0), R(size_t(N), 0.0), Al, Ad, Au, dA;
    std::vector<char> closed0(size_t(N), 1), canOpen(size_t(N), 0);
    bool capped = false;
    // section state: H, p_net, dp/dA, conductance A^3/(pi^2 H^2 mu) and its derivative
    auto section = [&](int i) {
        const double a = A[size_t(i)];
        if (a <= 0) { H[size_t(i)] = h; P[size_t(i)] = 0; dP[size_t(i)] = 2 * Ep / (M_PI * h * h); kA[size_t(i)] = 0; dkA[size_t(i)] = 0; return false; }
        const P3DSection S = p3dSection(a, h, ds, K, Ep, Hmax);
        H[size_t(i)] = S.H; P[size_t(i)] = S.pnet; dP[size_t(i)] = S.dpdA;
        const double c = 1 / (M_PI * M_PI * S.H * S.H * in.mu);
        kA[size_t(i)] = a * a * a * c;
        dkA[size_t(i)] = 3 * a * a * c - 2 * a * a * a * c / S.H * S.dHdA;
        return S.capped;
    };
    std::vector<double> tout;
    const int no = std::max(1, in.outputs);
    for (int k = 1; k <= no; ++k) tout.push_back(in.tEnd * std::pow(10.0, -3.0 * double(no - k) / std::max(1, no - 1)));
    if (no == 1) tout = {in.tEnd};
    size_t nextOut = 0;
    double t = 0, dt = std::min(in.tEnd * 1e-6, 1e-3), volInj = 0, volLeak = 0;
    auto snapshot = [&]() {
        Snapshot S; S.t = t;
        int tip = -1;
        double amax = 0; for (double v : A) amax = std::max(amax, v);
        for (int i = 0; i < N; ++i) if (A[size_t(i)] > 1e-9 * amax && amax > 0) tip = i;
        std::vector<double> w(size_t(N), 0.0);
        for (int i = 0; i < N; ++i) { section(i); w[size_t(i)] = 4 * A[size_t(i)] / (M_PI * H[size_t(i)]); }
        if (tip >= 0) {
            double L = (tip + 0.5) * dx;
            if (tip > 1) {   // PKN near-tip asymptote w^3 linear in the distance
                const double a = std::pow(w[size_t(tip - 1)], 3.0), b = std::pow(w[size_t(tip - 2)], 3.0);
                L = (tip - 0.5) * dx + std::clamp(b > a ? a / (b - a) * dx : dx, 0.0, 1.5 * dx);
            } else L = (tip + 1) * dx;
            S.L = L;
        }
        S.w0 = N > 1 ? std::max(0.0, 1.5 * w[0] - 0.5 * w[1]) : w[0];
        S.pnet0 = N > 1 ? 1.5 * P[0] - 0.5 * P[1] : P[0];
        double vf = 0, hs = 0; int no2 = 0;
        for (int i = 0; i < N; ++i) vf += dx * A[size_t(i)];
        vf *= 2;
        S.eff = volInj > 0 ? vf / volInj : 1.0;
        for (int i = 0; i <= tip; ++i) {
            S.x.push_back((i + 0.5) * dx); S.w.push_back(w[size_t(i)]); S.p.push_back(P[size_t(i)]); S.H.push_back(H[size_t(i)]);
            hs += H[size_t(i)]; ++no2;
        }
        out.volFracture = vf;
        out.height = no2 ? hs / no2 : h;
        out.heightWell = H[0];
        out.snaps.push_back(std::move(S));
    };
    int lastTip = 0;
    while (t < in.tEnd * (1 - 1e-12)) {
        double dtS = std::min(dt, in.tEnd - t);
        if (nextOut < tout.size() && t + dtS > tout[nextOut]) dtS = std::max(1e-12, tout[nextOut] - t);
        An = A;
        int tip = -1;
        for (int i = 0; i < N; ++i) { section(i); closed0[size_t(i)] = An[size_t(i)] <= 0; if (!closed0[size_t(i)]) tip = i; }
        for (int i = 0; i < N; ++i) canOpen[size_t(i)] = 1;
        if (K > 0 && tip >= 0 && tip + 1 < N) {   // lateral tip: PKN criterion with the local height
            if (P[size_t(tip)] < K / std::sqrt(0.5 * M_PI * H[size_t(tip)])) canOpen[size_t(tip + 1)] = 0;
            for (int j = tip + 2; j < N; ++j) canOpen[size_t(j)] = 0;
        }
        for (int i = 0; i < N; ++i) {   // Carter leak-off over both faces of the current height (explicit)
            Lv[size_t(i)] = 0;
            if (in.cl > 0 && tau[size_t(i)] >= 0)
                Lv[size_t(i)] = std::min(2 * H[size_t(i)] * dx * in.cl * 2 * (std::sqrt(t + dtS - tau[size_t(i)]) - std::sqrt(std::max(0.0, t - tau[size_t(i)]))),
                                         dx * An[size_t(i)] + dtS * q);
        }
        bool ok = false;
        for (int it = 0; it < 50; ++it) {
            bool cap = false;
            for (int i = 0; i < N; ++i) cap |= section(i);
            for (int i = 0; i < N; ++i) R[size_t(i)] = dx * (A[size_t(i)] - An[size_t(i)]) + Lv[size_t(i)] - (i == 0 ? dtS * q : 0.0);
            Ad.assign(size_t(N), dx); Al.assign(size_t(N), 0.0); Au.assign(size_t(N), 0.0);
            int last = 0;
            for (int i = 0; i < N; ++i) if (A[size_t(i)] > 0) last = i;
            const int iEnd = std::min(N - 1, last + 1);
            for (int i = 0; i < iEnd; ++i) {
                if (A[size_t(i)] <= 0 && A[size_t(i + 1)] <= 0) continue;
                if (!canOpen[size_t(i + 1)] && closed0[size_t(i + 1)]) continue;
                const double kf = 0.5 * (kA[size_t(i)] + kA[size_t(i + 1)]), dpd = (P[size_t(i)] - P[size_t(i + 1)]) / dx;
                const double F = kf * dpd;
                R[size_t(i)] += dtS * F; R[size_t(i + 1)] -= dtS * F;
                const double ga = 0.5 * dkA[size_t(i)] * dpd + kf * dP[size_t(i)] / dx;
                const double gb = 0.5 * dkA[size_t(i + 1)] * dpd - kf * dP[size_t(i + 1)] / dx;
                Ad[size_t(i)] += dtS * ga; Au[size_t(i)] += dtS * gb;
                Al[size_t(i + 1)] -= dtS * ga; Ad[size_t(i + 1)] -= dtS * gb;
            }
            double rmax = 0; for (double v : R) rmax = std::max(rmax, std::fabs(v));
            dA = R; for (double& v : dA) v = -v;
            if (!solveTri(Al, Ad, Au, dA)) break;
            double dmax = 0, amax = 0;
            for (int i = 0; i < N; ++i) {
                const double na = std::max(0.0, A[size_t(i)] + dA[size_t(i)]);
                dmax = std::max(dmax, std::fabs(na - A[size_t(i)]));
                A[size_t(i)] = na; amax = std::max(amax, na);
            }
            ++out.newton;
            if (dmax <= 1e-10 * std::max(amax, 1e-14) || (rmax <= 1e-13 * std::max(dtS * q, 1e-30) && it > 0)) { ok = true; capped |= cap; break; }
        }
        if (!ok) {
            A = An; dt *= 0.25;
            if (dt < in.tEnd * 1e-13) return fail("time step too small (Newton failed)");
            continue;
        }
        int newTip = 0;
        double stored = 0;
        for (int i = 0; i < N; ++i) {
            if (A[size_t(i)] > 0) { newTip = i; if (tau[size_t(i)] < 0) tau[size_t(i)] = t + 0.5 * dtS; }
            stored += dx * (A[size_t(i)] - An[size_t(i)]);
        }
        if (in.cl > 0) volLeak += std::max(0.0, dtS * q - stored) * 2;
        if (newTip >= N - 2) return fail("fracture reached the end of the computational domain - increase the domain length");
        volInj += in.rate * dtS;
        t += dtS; ++out.steps;
        const int adv = newTip - lastTip;
        lastTip = newTip;
        bool limited = false;
        if (K > 0 && adv >= 1 && newTip + 1 < N) { section(newTip); limited = P[size_t(newTip)] >= K / std::sqrt(0.5 * M_PI * H[size_t(newTip)]); }
        if (adv > 1 || limited) dt = dtS * 0.6; else if (adv == 1) dt = dtS * 1.1; else dt = dtS * 1.4;
        dt = std::min(dt, in.tEnd / 50);
        if (nextOut < tout.size() && t >= tout[nextOut] * (1 - 1e-9)) { snapshot(); ++nextOut; }
    }
    if (out.snaps.empty() || out.snaps.back().t < t * (1 - 1e-9)) snapshot();
    out.volInjected = volInj; out.volLeak = volLeak; out.heightRunaway = capped;
    return true;
}

bool simulate(const Input& in, Result& out) {
    out = Result();
    if (in.model == 2) {
        if (!(in.E > 0) || !(in.nu >= 0 && in.nu < 0.5)) { out.error = "elastic constants out of range"; return false; }
        if (!(in.height > 0) || !(in.mu > 0) || !(in.rate > 0) || !(in.tEnd > 0) || !(in.maxLength > 0)) { out.error = "height, viscosity, rate, time and domain must be positive"; return false; }
        if (in.cells < 10 || in.cells > 5000) { out.error = "cells per wing out of range"; return false; }
        if (in.cl < 0 || in.kic < 0) { out.error = "leak-off coefficient and toughness must be >= 0"; return false; }
        return simulateP3D(in, out);
    }
    auto fail = [&](const std::string& m) { out.error = m; return false; };
    if (!(in.E > 0) || !(in.nu >= 0 && in.nu < 0.5)) return fail("elastic constants out of range");
    if (!(in.height > 0) || !(in.mu > 0) || !(in.rate > 0) || !(in.tEnd > 0) || !(in.maxLength > 0)) return fail("height, viscosity, rate, time and domain must be positive");
    if (in.cells < 10 || in.cells > (in.model == 1 ? 600 : 5000)) return fail("cells per wing out of range");
    if (in.cl < 0 || in.kic < 0) return fail("leak-off coefficient and toughness must be >= 0");
    const int N = in.cells;
    out.height = in.height;
    const double dx = in.maxLength / N, Ep = in.E / (1 - in.nu * in.nu), H = in.height;
    const bool pkn = in.model == 0;
    const double q = pkn ? 0.5 * in.rate : 0.5 * in.rate / H;     // per wing (PKN m3/s, KGD m2/s)
    const double sc = pkn ? 0.25 * M_PI * H * dx : dx;            // storage per unit width
    const double kc = pkn ? M_PI * H / (64 * in.mu) : 1.0 / (12 * in.mu);
    const double la = pkn ? 2 * H * dx : 2 * dx;                  // leak-off area (both faces)
    const double eP = Ep / (2 * H);                               // PKN p_net / w
    const double vScale = pkn ? 2.0 : 2.0 * H;                    // both wings (KGD: times height) -> m3
    std::vector<double> C;
    if (!pkn) {
        C.assign(size_t(N) * N, 0.0);
        const double a = 0.5 * dx, f = Ep / (4 * M_PI);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                const double xi = (i + 0.5) * dx, xj = (j + 0.5) * dx;
                C[size_t(i) * N + j] = f * (1 / (xi - xj + a) - 1 / (xi - xj - a) + 1 / (xi + xj + a) - 1 / (xi + xj - a));
            }
    }
    auto pressure = [&](const std::vector<double>& w, std::vector<double>& p) {
        if (pkn) { for (int i = 0; i < N; ++i) p[size_t(i)] = eP * w[size_t(i)]; return; }
        for (int i = 0; i < N; ++i) { double s = 0; const double* row = &C[size_t(i) * N]; for (int j = 0; j < N; ++j) if (w[size_t(j)] > 0) s += row[j] * w[size_t(j)]; p[size_t(i)] = s; }
    };
    const double pCritPkn = in.kic / std::sqrt(0.5 * M_PI * H);
    std::vector<double> w(size_t(N), 0.0), wn, p(size_t(N), 0.0), tau(size_t(N), -1.0), Lv(size_t(N), 0.0), R(size_t(N), 0.0), dw;
    std::vector<char> closed0(size_t(N), 1), canOpen(size_t(N), 0);
    // output times: log-spaced over three decades, the last at tEnd
    std::vector<double> tout;
    const int no = std::max(1, in.outputs);
    for (int k = 1; k <= no; ++k) tout.push_back(in.tEnd * std::pow(10.0, -3.0 * double(no - k) / std::max(1, no - 1)));
    if (no == 1) tout = {in.tEnd};
    size_t nextOut = 0;
    double t = 0, dt = std::min(in.tEnd * 1e-6, 1e-3), volInj = 0, volLeak = 0;
    auto snapshot = [&]() {
        Snapshot S; S.t = t;
        int tip = -1;
        double wmax = 0; for (double v : w) wmax = std::max(wmax, v);
        for (int i = 0; i < N; ++i) if (w[size_t(i)] > 1e-9 * wmax && wmax > 0) tip = i;
        if (tip >= 0) {
            const double ex = pkn ? 3.0 : (in.kic > 0 ? 2.0 : 1.5);   // near-tip asymptote w^ex linear in the distance (KGD: viscous 2/3, LEFM 1/2)
            double L = (tip + 0.5) * dx;
            if (tip > 1) {
                // extrapolate w^ex linearly from the two cells behind the (partially filled) tip cell
                const double a = std::pow(w[size_t(tip - 1)], ex), b = std::pow(w[size_t(tip - 2)], ex);
                L = (tip - 0.5) * dx + std::clamp(b > a ? a / (b - a) * dx : dx, 0.0, 1.5 * dx);
            } else L = (tip + 1) * dx;
            S.L = L;
        }
        pressure(w, p);
        S.w0 = N > 1 ? std::max(0.0, 1.5 * w[0] - 0.5 * w[1]) : w[0];
        S.pnet0 = N > 1 ? 1.5 * p[0] - 0.5 * p[1] : p[0];
        double vf = 0; for (int i = 0; i < N; ++i) vf += sc * w[size_t(i)];
        vf *= vScale;
        S.eff = volInj > 0 ? vf / volInj : 1.0;
        for (int i = 0; i <= tip; ++i) { S.x.push_back((i + 0.5) * dx); S.w.push_back(w[size_t(i)]); S.p.push_back(p[size_t(i)]); }
        out.snaps.push_back(std::move(S));
        out.volFracture = vf;
    };
    int lastTip = 0;
    while (t < in.tEnd * (1 - 1e-12)) {
        double dtS = std::min(dt, in.tEnd - t);
        if (nextOut < tout.size() && t + dtS > tout[nextOut]) dtS = std::max(1e-12, tout[nextOut] - t);
        wn = w;
        pressure(wn, p);
        // a closed cell next to the tip can only receive fluid when the tip stress intensity reaches K_IC
        int tip = -1;
        for (int i = 0; i < N; ++i) { closed0[size_t(i)] = wn[size_t(i)] <= 0; if (!closed0[size_t(i)]) tip = i; }
        for (int i = 0; i < N; ++i) canOpen[size_t(i)] = 1;
        if (in.kic > 0 && tip >= 0 && tip + 1 < N) {
            bool grow;
            if (pkn) grow = p[size_t(tip)] >= pCritPkn;
            // displacement correlation on the tip element, calibrated for constant DD elements: a uniformly loaded crack
            // gives K_est / K_exact = 1.253 (n = 5..200 elements: 1.238..1.253)
            else grow = Ep * wn[size_t(tip)] * std::sqrt(2 * M_PI) / (8 * std::sqrt(0.5 * dx)) / 1.253 >= in.kic;
            if (!grow) canOpen[size_t(tip + 1)] = 0;
            // with toughness the criterion is re-evaluated after every new cell: at most one cell opens per step
            for (int j = tip + 2; j < N; ++j) canOpen[size_t(j)] = 0;
        }
        // Carter leak-off volume of the step for exposed cells (explicit in the exposure times)
        for (int i = 0; i < N; ++i) {
            Lv[size_t(i)] = 0;
            if (in.cl > 0 && tau[size_t(i)] >= 0)
                Lv[size_t(i)] = std::min(la * in.cl * 2 * (std::sqrt(t + dtS - tau[size_t(i)]) - std::sqrt(std::max(0.0, t - tau[size_t(i)]))), sc * wn[size_t(i)] + dtS * q);
        }
        bool ok = false;
        for (int it = 0; it < 40; ++it) {
            pressure(w, p);
            // pressure driving the fluid: KGD closed cells ahead of the tip carry the elastic (suction) pressure of the
            // opening behind them - conservative and converging from above (zero net pressure there under-predicts L by
            // ~5 % and loses mass through the w >= 0 clamp)
            const std::vector<double>& pf = p;
            for (int i = 0; i < N; ++i) R[size_t(i)] = sc * (w[size_t(i)] - wn[size_t(i)]) + Lv[size_t(i)] - (i == 0 ? dtS * q : 0.0);
            std::vector<double> Ad, Al, Au;
            std::vector<double> J;
            if (pkn) { Ad.assign(size_t(N), sc); Al.assign(size_t(N), 0.0); Au.assign(size_t(N), 0.0); }
            else { J.assign(size_t(N) * N, 0.0); for (int i = 0; i < N; ++i) J[size_t(i) * N + i] = sc; }
            int last = 0;
            for (int i = 0; i < N; ++i) if (w[size_t(i)] > 0) last = i;
            const int iEnd = std::min(N - 1, last + 1);
            for (int i = 0; i < iEnd; ++i) {
                const double wa = w[size_t(i)], wb = w[size_t(i + 1)];
                if (wa <= 0 && wb <= 0) continue;
                if (!canOpen[size_t(i + 1)] && closed0[size_t(i + 1)]) continue;
                const double wf3 = 0.5 * (wa * wa * wa + wb * wb * wb), dpd = (pf[size_t(i)] - pf[size_t(i + 1)]) / dx;
                const double F = kc * wf3 * dpd;
                R[size_t(i)] += dtS * F; R[size_t(i + 1)] -= dtS * F;
                const double dFa = kc * 1.5 * wa * wa * dpd, dFb = kc * 1.5 * wb * wb * dpd;
                if (pkn) {
                    const double ga = dFa + kc * wf3 * eP / dx, gb = dFb - kc * wf3 * eP / dx;
                    Ad[size_t(i)] += dtS * ga; Au[size_t(i)] += dtS * gb;
                    Al[size_t(i + 1)] -= dtS * ga; Ad[size_t(i + 1)] -= dtS * gb;
                } else {
                    double* Ja = &J[size_t(i) * N]; double* Jb = &J[size_t(i + 1) * N];
                    const double kf = dtS * kc * wf3 / dx;
                    for (int k = 0; k < N; ++k) {
                        const double ca = C[size_t(i) * N + k], cb = C[size_t(i + 1) * N + k];
                        const double g = kf * (ca - cb);
                        Ja[k] += g; Jb[k] -= g;
                    }
                    Ja[i] += dtS * dFa; Ja[i + 1] += dtS * dFb;
                    Jb[i] -= dtS * dFa; Jb[i + 1] -= dtS * dFb;
                }
            }
            double rmax = 0; for (double v : R) rmax = std::max(rmax, std::fabs(v));
            dw = R;
            for (double& v : dw) v = -v;
            const bool solved = pkn ? solveTri(Al, Ad, Au, dw) : solveDense(J, dw, N);
            if (!solved) break;
            double dmax = 0, wmax = 0;
            for (int i = 0; i < N; ++i) {
                const double nw = std::max(0.0, w[size_t(i)] + dw[size_t(i)]);
                dmax = std::max(dmax, std::fabs(nw - w[size_t(i)]));
                w[size_t(i)] = nw; wmax = std::max(wmax, nw);
            }
            ++out.newton;
            if (dmax <= 1e-10 * std::max(wmax, 1e-12) || (rmax <= 1e-13 * std::max(dtS * q, 1e-30) && it > 0)) { ok = true; break; }
        }
        if (!ok) {
            w = wn; dt *= 0.25;
            if (dt < in.tEnd * 1e-13) return fail("time step too small (Newton failed)");
            continue;
        }
        // exposure times of the newly opened cells; leak-off actually realised = injected - stored (a cell that runs dry
        // cannot leak its full Carter volume: w >= 0 is a complementarity condition)
        int newTip = 0;
        double stored = 0;
        for (int i = 0; i < N; ++i) {
            if (w[size_t(i)] > 0) { newTip = i; if (tau[size_t(i)] < 0) tau[size_t(i)] = t + 0.5 * dtS; }
            stored += sc * (w[size_t(i)] - wn[size_t(i)]);
        }
        if (in.cl > 0) volLeak += std::max(0.0, dtS * q - stored) * vScale;
        if (newTip >= N - 2) return fail("fracture reached the end of the computational domain - increase the domain length");
        volInj += in.rate * dtS;
        t += dtS; ++out.steps;
        // step size: about one cell of front advance per step
        const int adv = newTip - lastTip;
        lastTip = newTip;
        bool limited = false;   // toughness: the front advanced and could already advance again - the step limited it
        if (in.kic > 0 && adv >= 1 && newTip + 1 < N) {
            if (pkn) limited = eP * w[size_t(newTip)] >= pCritPkn;
            else { pressure(w, p); limited = Ep * w[size_t(newTip)] * std::sqrt(2 * M_PI) / (8 * std::sqrt(0.5 * dx)) / 1.253 >= in.kic; }
        }
        if (adv > 1 || limited) dt = dtS * 0.6; else if (adv == 1) dt = dtS * 1.1; else dt = dtS * 1.4;
        dt = std::min(dt, in.tEnd / 50);
        if (nextOut < tout.size() && t >= tout[nextOut] * (1 - 1e-9)) { snapshot(); ++nextOut; }
    }
    if (out.snaps.empty() || out.snaps.back().t < t * (1 - 1e-9)) snapshot();
    out.volInjected = volInj; out.volLeak = volLeak;
    return true;
}

// ---------------------------------------------------------------- C17 multi-cluster growth with stress shadowing
bool simulateClusters(const ClusterInput& ci, ClusterResult& out) {
    out = ClusterResult();
    const Input& in = ci.base;
    auto fail = [&](const std::string& m) { out.error = m; return false; };
    if (in.model != 0) return fail("multi-cluster growth uses the PKN model");
    if (!(in.E > 0) || !(in.nu >= 0 && in.nu < 0.5)) return fail("elastic constants out of range");
    if (!(in.height > 0) || !(in.mu > 0) || !(in.rate > 0) || !(in.tEnd > 0) || !(in.maxLength > 0)) return fail("height, viscosity, rate, time and domain must be positive");
    if (ci.clusters < 1 || ci.clusters > 50) return fail("1 - 50 clusters");
    if (in.cells < 10 || in.cells > 5000) return fail("cells per wing out of range");
    if (in.cl < 0 || !(ci.perfFriction >= 0) || !(ci.spacing > 0)) return fail("leak-off, perforation friction and spacing must be valid");
    const int N = in.cells, M = ci.clusters;
    const double dx = in.maxLength / N, Ep = in.E / (1 - in.nu * in.nu), H = in.height, Q = in.rate;
    const double sc = 0.25 * M_PI * H * dx, kc = M_PI * H / (64 * in.mu), la = 2 * H * dx, eP = Ep / (2 * H), kp = ci.perfFriction;
    std::vector<double> pos(static_cast<size_t>(M), 0.0), sOff(static_cast<size_t>(M), 0.0);
    for (int k = 0; k < M; ++k) { pos[size_t(k)] = k * ci.spacing; if (size_t(k) < ci.stressOffset.size()) sOff[size_t(k)] = ci.stressOffset[size_t(k)]; }
    std::vector<std::vector<double>> w(static_cast<size_t>(M), std::vector<double>(static_cast<size_t>(N), 0.0)), wn = w,
        tau(static_cast<size_t>(M), std::vector<double>(static_cast<size_t>(N), -1.0)), Lv = w;
    std::vector<double> q(static_cast<size_t>(M), Q / M), qn, shadow(static_cast<size_t>(M), 0.0), Lk(static_cast<size_t>(M), 0.0), pbar(static_cast<size_t>(M), 0.0);
    std::vector<int> lastTip(static_cast<size_t>(M), 0);
    double pw = 0, pwn = 0;
    for (int k = 0; k < M; ++k) pw = std::max(pw, sOff[size_t(k)]);
    std::vector<double> tout;
    const int no = std::max(1, in.outputs);
    for (int k = 1; k <= no; ++k) tout.push_back(in.tEnd * std::pow(10.0, -3.0 * double(no - k) / std::max(1, no - 1)));
    if (no == 1) tout = {in.tEnd};
    size_t nextOut = 0;
    double t = 0, dt = std::min(in.tEnd * 1e-6, 1e-3), volInj = 0, volLeak = 0;
    auto geometry = [&]() {   // half-lengths and mean net pressures (tip: last open cell, partial fill by the w^3 extrapolation)
        for (int k = 0; k < M; ++k) {
            const auto& wk = w[size_t(k)];
            int tip = -1; double sw = 0;
            for (int i = 0; i < N; ++i) if (wk[size_t(i)] > 0) { tip = i; sw += wk[size_t(i)]; }
            double L = 0;
            if (tip > 1) {
                const double a = std::pow(wk[size_t(tip - 1)], 3.0), b = std::pow(wk[size_t(tip - 2)], 3.0);
                L = (tip - 0.5) * dx + std::clamp(b > a ? a / (b - a) * dx : dx, 0.0, 1.5 * dx);
            } else if (tip >= 0) L = (tip + 1) * dx;
            Lk[size_t(k)] = L;
            pbar[size_t(k)] = tip >= 0 ? eP * sw / (tip + 1) : 0.0;
        }
    };
    auto record = [&]() {
        geometry();
        out.t.push_back(t); out.pw.push_back(pw); out.q.push_back(q); out.L.push_back(Lk);
    };
    std::vector<double> R(static_cast<size_t>(N), 0.0), Ad, Al, Au, u, v;
    std::vector<std::vector<double>> U(static_cast<size_t>(M)), V(static_cast<size_t>(M));
    while (t < in.tEnd * (1 - 1e-12)) {
        double dtS = std::min(dt, in.tEnd - t);
        if (nextOut < tout.size() && t + dtS > tout[nextOut]) dtS = std::max(1e-12, tout[nextOut] - t);
        wn = w; qn = q; pwn = pw;
        // stress shadow of the other clusters (explicit)
        geometry();
        for (int k = 0; k < M; ++k) {
            double s = 0;
            if (ci.shadow)
                for (int j = 0; j < M; ++j) {
                    if (j == k || !(Lk[size_t(j)] > 0)) continue;
                    const double c = 0.5 * std::min(H, 2 * Lk[size_t(j)]), xi = std::fabs(pos[size_t(k)] - pos[size_t(j)]) / c;
                    s += pbar[size_t(j)] * (1 - xi * xi * xi / std::pow(1 + xi * xi, 1.5));
                }
            shadow[size_t(k)] = s;
        }
        for (int k = 0; k < M; ++k)
            for (int i = 0; i < N; ++i) {
                double& lv = Lv[size_t(k)][size_t(i)]; lv = 0;
                const double tk = tau[size_t(k)][size_t(i)];
                if (in.cl > 0 && tk >= 0) lv = std::min(la * in.cl * 2 * (std::sqrt(t + dtS - tk) - std::sqrt(std::max(0.0, t - tk))), sc * wn[size_t(k)][size_t(i)] + dtS * 0.5 * Q);
            }
        bool ok = false;
        for (int it = 0; it < 50; ++it) {
            double sAlpha = 0, sBeta = 0, sq = 0;
            std::vector<double> alpha(static_cast<size_t>(M), 0.0), beta(static_cast<size_t>(M), 0.0);
            for (int k = 0; k < M; ++k) {
                auto& wk = w[size_t(k)]; const auto& wnk = wn[size_t(k)];
                for (int i = 0; i < N; ++i) R[size_t(i)] = sc * (wk[size_t(i)] - wnk[size_t(i)]) + Lv[size_t(k)][size_t(i)] - (i == 0 ? dtS * 0.5 * q[size_t(k)] : 0.0);
                Ad.assign(size_t(N), sc); Al.assign(size_t(N), 0.0); Au.assign(size_t(N), 0.0);
                int last = 0;
                for (int i = 0; i < N; ++i) if (wk[size_t(i)] > 0) last = i;
                const int iEnd = std::min(N - 1, last + 1);
                for (int i = 0; i < iEnd; ++i) {
                    const double wa = wk[size_t(i)], wb = wk[size_t(i + 1)];
                    if (wa <= 0 && wb <= 0) continue;
                    const double wf3 = 0.5 * (wa * wa * wa + wb * wb * wb), dpd = eP * (wa - wb) / dx, F = kc * wf3 * dpd;
                    R[size_t(i)] += dtS * F; R[size_t(i + 1)] -= dtS * F;
                    const double ga = kc * 1.5 * wa * wa * dpd + kc * wf3 * eP / dx, gb = kc * 1.5 * wb * wb * dpd - kc * wf3 * eP / dx;
                    Ad[size_t(i)] += dtS * ga; Au[size_t(i)] += dtS * gb;
                    Al[size_t(i + 1)] -= dtS * ga; Ad[size_t(i + 1)] -= dtS * gb;
                }
                u = R; for (double& x : u) x = -x;
                v.assign(size_t(N), 0.0); v[0] = 0.5 * dtS;   // dR_0 / dq_k = -dt / 2  ->  A v = dt / 2 e_0
                if (!solveTri(Al, Ad, Au, u) || !solveTri(Al, Ad, Au, v)) { ok = false; it = 99; break; }
                U[size_t(k)] = u; V[size_t(k)] = v;
                const double qk = q[size_t(k)];
                const double g = sOff[size_t(k)] + shadow[size_t(k)] + eP * wk[0] + kp * qk * qk - pw;
                if (qk <= 0 && g > 0) { alpha[size_t(k)] = -qk; beta[size_t(k)] = 0; }   // closed cluster: needs more than p_w
                else {
                    const double den = eP * v[0] + 2 * kp * qk;
                    alpha[size_t(k)] = (-g - eP * u[0]) / den; beta[size_t(k)] = 1 / den;
                }
                sAlpha += alpha[size_t(k)]; sBeta += beta[size_t(k)]; sq += qk;
            }
            if (it >= 99) break;
            if (!(sBeta > 0)) break;
            const double dpw = (Q - sq - sAlpha) / sBeta;
            pw += dpw;
            double dmax = 0, wmax = 0;
            for (int k = 0; k < M; ++k) {
                const double qNew = std::max(0.0, q[size_t(k)] + alpha[size_t(k)] + beta[size_t(k)] * dpw), dq = qNew - q[size_t(k)];
                q[size_t(k)] = qNew;
                auto& wk = w[size_t(k)];
                for (int i = 0; i < N; ++i) {
                    const double nw = std::max(0.0, wk[size_t(i)] + U[size_t(k)][size_t(i)] + V[size_t(k)][size_t(i)] * dq);
                    dmax = std::max(dmax, std::fabs(nw - wk[size_t(i)]));
                    wk[size_t(i)] = nw; wmax = std::max(wmax, nw);
                }
            }
            ++out.newton;
            if (it > 0 && dmax <= 1e-10 * std::max(wmax, 1e-12) && std::fabs(dpw) <= 1e-7 * std::max(1.0, std::fabs(pw))) { ok = true; break; }
        }
        if (!ok) {
            w = wn; q = qn; pw = pwn; dt *= 0.25;
            if (dt < in.tEnd * 1e-13) return fail("time step too small (Newton failed)");
            continue;
        }
        int adv = 0;
        for (int k = 0; k < M; ++k) {
            int newTip = 0; double stored = 0;
            for (int i = 0; i < N; ++i) {
                if (w[size_t(k)][size_t(i)] > 0) { newTip = i; if (tau[size_t(k)][size_t(i)] < 0) tau[size_t(k)][size_t(i)] = t + 0.5 * dtS; }
                stored += sc * (w[size_t(k)][size_t(i)] - wn[size_t(k)][size_t(i)]);
            }
            if (in.cl > 0) volLeak += std::max(0.0, dtS * 0.5 * q[size_t(k)] - stored) * 2;
            if (newTip >= N - 2) return fail("a cluster fracture reached the end of the computational domain - increase the domain length");
            adv = std::max(adv, newTip - lastTip[size_t(k)]);
            lastTip[size_t(k)] = newTip;
        }
        volInj += Q * dtS;
        t += dtS; ++out.steps;
        if (adv > 1) dt = dtS * 0.6; else if (adv == 1) dt = dtS * 1.1; else dt = dtS * 1.4;
        dt = std::min(dt, in.tEnd / 50);
        if (nextOut < tout.size() && t >= tout[nextOut] * (1 - 1e-9)) { record(); ++nextOut; }
    }
    if (out.t.empty() || out.t.back() < t * (1 - 1e-9)) record();
    geometry();
    out.Lf = Lk; out.qf = q; out.shadowf = shadow;
    out.w0f.assign(size_t(M), 0.0); out.volf.assign(size_t(M), 0.0);
    for (int k = 0; k < M; ++k) {
        const auto& wk = w[size_t(k)];
        out.w0f[size_t(k)] = N > 1 ? std::max(0.0, 1.5 * wk[0] - 0.5 * wk[1]) : wk[0];
        for (int i = 0; i < N; ++i) out.volf[size_t(k)] += 2 * sc * wk[size_t(i)];
        out.volFracture += out.volf[size_t(k)];
    }
    out.volInjected = volInj; out.volLeak = volLeak;
    return true;
}

bool simulateAdaptive(const Input& base, Result& out) {
    // domain from the no-leak-off similarity length (an upper bound for the half-length), then refit so the final
    // fracture spans most of the grid (the Carter tip exposure is resolved to first order in the cell size)
    Input in = base;
    const double L0 = in.model == 1 ? kgdLength(in, in.tEnd) : nordgrenLength(in, in.tEnd);   // P3D: PKN length is an upper bound
    in.maxLength = std::max(1.6 * L0, 1.0);
    if (in.cells <= 0) in.cells = in.model == 1 ? 150 : 400;
    for (int pass = 0; pass < 4; ++pass) {
        if (!simulate(in, out)) {
            if (out.error.find("end of the computational domain") == std::string::npos) return false;
            in.maxLength *= 2;
            continue;
        }
        const double L = out.snaps.back().L;
        if (L > 0.45 * in.maxLength || pass == 3) return true;
        in.maxLength = std::max(1.5 * L, 4 * in.maxLength / in.cells);
    }
    return !out.snaps.empty() && out.error.empty();
}


// ------------------------------------------------------------------------------------------------ C15 stage layout
std::vector<std::array<double, 3>> wellPath(const Case& c, const std::string& well) {
    std::vector<std::array<double, 3>> pts;
    const WellInput* W = nullptr;
    for (const WellInput& w : c.wells) if (w.name == well) W = &w;
    if (!W) return pts;
    Case c0 = c;
    c0.grid.hfracs.clear(); c0.grid.hfExpanded = 0; c0.grid.hfPieces.clear();
    Grid G;
    G.build(c0);
    const GridInput& g = c0.grid;
    for (const Completion& cp : W->completions) {
        if (cp.lgr >= 0) continue;
        for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
            if (cp.i < 0 || cp.i >= g.nx || cp.j < 0 || cp.j >= g.ny || k < 0 || k >= g.nz) continue;
            const int kk = (g.dual == 1 && k < g.nzMatrix()) ? k + g.nzMatrix() : k;
            const Vec3 v = G.center[size_t(g.index(cp.i, cp.j, kk))];
            const std::array<double, 3> p{double(v.x), double(v.y), double(v.z)};
            if (pts.empty() || std::fabs(pts.back()[0] - p[0]) + std::fabs(pts.back()[1] - p[1]) + std::fabs(pts.back()[2] - p[2]) > 1e-9) pts.push_back(p);
        }
    }
    return pts;
}

std::vector<GridInput::HFrac> layoutStages(const Case& c, const StageLayout& L) {
    const auto P = wellPath(c, L.well);
    if (P.empty()) throw std::runtime_error("well " + L.well + " has no completion in the grid");
    if (L.stages < 1) throw std::runtime_error("at least one stage");
    if (!(L.halfLength > 0) || !(L.height > 0) || !(L.width > 0)) throw std::runtime_error("half-length, height and width must be > 0");
    std::vector<double> s(P.size(), 0.0);
    for (size_t i = 1; i < P.size(); ++i) s[i] = s[i - 1] + std::sqrt(std::pow(P[i][0] - P[i - 1][0], 2) + std::pow(P[i][1] - P[i - 1][1], 2) + std::pow(P[i][2] - P[i - 1][2], 2));
    const double Ltot = s.back();
    std::vector<GridInput::HFrac> out;
    for (int k = 0; k < L.stages; ++k) {
        double sk;
        if (L.spacing > 0) sk = Ltot - L.toeOffset - k * L.spacing;
        else sk = Ltot * (k + 0.5) / L.stages;
        if (P.size() > 1 && (sk < -1e-9 || sk > Ltot + 1e-9)) throw std::runtime_error("stage " + std::to_string(k + 1) + " falls outside the completed length (" + std::to_string(Ltot) + ")");
        sk = std::min(std::max(sk, 0.0), Ltot);
        size_t seg = 0;
        while (seg + 2 < P.size() && s[seg + 1] < sk) ++seg;
        double x = P[0][0], y = P[0][1], z = P[0][2], dx = 1, dy = 0;
        if (P.size() > 1) {
            const double f = (sk - s[seg]) / std::max(1e-12, s[seg + 1] - s[seg]);
            x = P[seg][0] + f * (P[seg + 1][0] - P[seg][0]); y = P[seg][1] + f * (P[seg + 1][1] - P[seg][1]); z = P[seg][2] + f * (P[seg + 1][2] - P[seg][2]);
            dx = P[seg + 1][0] - P[seg][0]; dy = P[seg + 1][1] - P[seg][1];
        }
        double ux, uy;   // fracture trace direction
        if (L.strikeDeg >= 0) { const double a = L.strikeDeg * M_PI / 180.0; ux = std::sin(a); uy = std::cos(a); }
        else {
            const double l = std::hypot(dx, dy);
            if (l < 1e-9) { ux = 1; uy = 0; } else { ux = -dy / l; uy = dx / l; }
        }
        GridInput::HFrac F;
        F.name = L.prefix + std::to_string(k + 1);
        F.x0 = x - L.halfLength * ux; F.y0 = y - L.halfLength * uy; F.x1 = x + L.halfLength * ux; F.y1 = y + L.halfLength * uy;
        F.zTop = std::max(0.0, z - 0.5 * L.height); F.zBot = z + 0.5 * L.height;
        F.width = L.width; F.perm = L.perm; F.poro = L.poro; F.gamma = L.gamma;
        out.push_back(F);
    }
    return out;
}

} // namespace hf
} // namespace ressim
