// H7 (Research Cycle 7): is the induced-event count convex in the injection-rate history when the seismicity comes from
// a quasi-dynamic rate-and-state fault (elastic interaction, nucleation, stress drops) instead of the Dieterich rate model?
//
// Pre-registered design (CYCLE7_PLAN.md):
//   fault    2-D rate-and-state (rsf.h), 1500 m, 64 elements, velocity-weakening centre (25-75 %), no tectonic loading
//            (Vpl = 0): every event is injection-induced; initial steady state at vInit = 1e-9 m/s (tau0 = f_ss sigma')
//   forcing  line-source well 20 m from the fault centre, 2-D Theis superposition of a piecewise-constant rate history
//            q(t) = u_k q_ref over K = 4 intervals of 100 d, mean(u) = 1 (fixed injected volume), u_k in [0, 2.5]
//            dp(x, t) = P* sum_k (u_k - u_{k-1}) E1(r^2 / (4 c (t - t_k))),  P* = 1.5 MPa, c = 0.05 m2/s; window 400 d
//   test     M = 30 random pairs (A, B), midpoint C = (A + B) / 2 (same volume): violation if N(C) > (N(A) + N(B)) / 2
//            (strict; event count N over the window), also for the cumulative seismic moment M0 and for the
//            Dieterich expected count with the same pressure at the fault centre (P1 model, convex by F1)
//   verdict  convexity SUPPORTED if strict count violations <= 10 % of the pairs, REFUTED if > 25 %, else inconclusive
// Output: one line per run + summary; JSON-like lines "PAIR ..." for the report.
#include "ressim/rsf.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
using namespace ressim;

static const double day = 86400.0;
static const int K = 4;
static const double Tk = 100 * day, Pstar = 1.5e6, cdiff = 0.05, offset = 20.0, L = 1500.0;

static double E1(double x) { return x > 50 ? 0.0 : -std::expint(-x); }

static double dpAt(const std::vector<double>& u, double r2, double t) {
    double p = 0, prev = 0;
    for (int k = 0; k < K; ++k) {
        const double tk = k * Tk;
        if (t <= tk) break;
        p += (u[size_t(k)] - prev) * E1(r2 / (4 * cdiff * (t - tk)));
        prev = u[size_t(k)];
    }
    return Pstar * p;
}

struct Res { int n = 0; double m0 = 0, dieterich = 0, sec = 0; bool ok = false; };

static Res runOne(const std::vector<double>& u) {
    rsf::Input in; in.n = 64; in.x1 = L; in.Vpl = 0; in.tEnd = K * Tk; in.vwFrom = 0.25; in.vwTo = 0.75;
    in.dtMax = 2 * day;
    const double dx = L / in.n;
    std::vector<double> r2(size_t(in.n));
    for (int i = 0; i < in.n; ++i) { const double s = (i + 0.5) * dx - 0.5 * L; r2[size_t(i)] = s * s + offset * offset; }
    in.forcing = [&](double t, std::vector<double>& dp, std::vector<double>& dtau) {
        (void)dtau;
        for (int i = 0; i < in.n; ++i) dp[size_t(i)] += dpAt(u, r2[size_t(i)], t);
    };
    rsf::Output o; Res R;
    const auto t0 = std::chrono::steady_clock::now();
    R.ok = rsf::run(in, o);
    R.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    R.n = int(o.events.size());
    for (const auto& e : o.events) R.m0 += e.moment;
    // Dieterich (1994) expected count at the fault centre with the same pressure: S(t) = -mu dp (Coulomb, f 0.6),
    // A sigma = 0.01 * 50 MPa, background rate r = 1 / yr, t_a = A sigma / tau_dot (tau_dot = 1 kPa/yr) -> R/r = 1 / gamma
    const double Asig = 0.01 * 50e6, taudot = 1e3 / (365.25 * day), r = 1.0 / (365.25 * day);
    double gam = 1.0 / taudot, prevS = 0, cnt = 0;
    const int nt = 4000; const double dt = in.tEnd / nt;
    for (int j = 1; j <= nt; ++j) {
        const double t = j * dt, S = 0.6 * dpAt(u, r2[size_t(in.n / 2)], t), dS = (S - prevS) / dt + taudot; prevS = S;
        // d gamma / dt = (1 - gamma dS/dt) / A sigma  (exact update for constant dS over the step)
        const double a = dS / Asig;
        if (std::fabs(a) > 1e-30) gam = (gam - 1 / (Asig * a)) * std::exp(-a * dt) + 1 / (Asig * a); else gam += dt / Asig;
        cnt += r / (gam * taudot) * dt * 1e3;   // per 1000 background events
    }
    R.dieterich = cnt;
    return R;
}

int main(int argc, char** argv) {
    const int M = argc > 1 ? std::atoi(argv[1]) : 30;
    const unsigned seed = argc > 2 ? unsigned(std::atoi(argv[2])) : 7u;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    auto randomSchedule = [&]() {
        std::vector<double> u(K);
        for (;;) {
            double s = 0;
            for (double& v : u) { v = U(rng); s += v; }
            for (double& v : u) v *= K / s;
            if (*std::max_element(u.begin(), u.end()) <= 2.5) return u;
        }
    };
    const std::vector<double> ones(K, 1.0);
    const Res c0 = runOne(ones);
    std::printf("CONST n %d M0 %.4g Dieterich %.4g (%.1f s)\n", c0.n, c0.m0, c0.dieterich, c0.sec);
    int vN = 0, vNtol = 0, vM = 0, vD = 0, okPairs = 0, eqN = 0;
    double sec = c0.sec;
    for (int p = 0; p < M; ++p) {
        const auto A = randomSchedule(), B = randomSchedule();
        std::vector<double> C(K);
        for (int k = 0; k < K; ++k) C[size_t(k)] = 0.5 * (A[size_t(k)] + B[size_t(k)]);
        const Res ra = runOne(A), rb = runOne(B), rc = runOne(C);
        sec += ra.sec + rb.sec + rc.sec;
        if (!(ra.ok && rb.ok && rc.ok)) { std::printf("PAIR %d failed\n", p); continue; }
        ++okPairs;
        const double avgN = 0.5 * (ra.n + rb.n), avgM = 0.5 * (ra.m0 + rb.m0), avgD = 0.5 * (ra.dieterich + rb.dieterich);
        const bool vn = rc.n > avgN + 1e-9, vnt = rc.n > avgN + 1.0, vm = rc.m0 > avgM * (1 + 1e-6), vd = rc.dieterich > avgD * (1 + 1e-9);
        vN += vn; vNtol += vnt; vM += vm; vD += vd; eqN += (ra.n == rb.n && rb.n == rc.n);
        std::printf("PAIR %d A [%.2f %.2f %.2f %.2f] B [%.2f %.2f %.2f %.2f] N %d %d mid %d | M0 %.3g %.3g mid %.3g | D %.4g %.4g mid %.4g | %s%s%s\n", p,
                    A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], ra.n, rb.n, rc.n, ra.m0, rb.m0, rc.m0, ra.dieterich, rb.dieterich, rc.dieterich,
                    vn ? "N-VIOL " : "", vm ? "M0-VIOL " : "", vd ? "D-VIOL" : "");
        std::fflush(stdout);
    }
    std::printf("SUMMARY pairs %d | count violations strict %d (%.0f %%), > 1 event %d | moment violations %d (%.0f %%) | Dieterich violations %d | all-equal count pairs %d | %.0f s\n",
                okPairs, vN, 100.0 * vN / std::max(1, okPairs), vNtol, vM, 100.0 * vM / std::max(1, okPairs), vD, eqN, sec);
    const double f = double(vN) / std::max(1, okPairs);
    std::printf("VERDICT %s\n", f <= 0.10 ? "SUPPORTED (<= 10 %)" : f > 0.25 ? "REFUTED (> 25 %)" : "INCONCLUSIVE (10-25 %)");
}
