// H9 (Research Cycle 9, pre-registered in CYCLE9_PLAN.md): does a Dieterich surrogate, calibrated on a few rate-and-state
// (RSF) runs, find injection schedules as good as a stochastic search on the RSF fault itself?
//
// Objective (same as H7b / F22): total event count of an E = 16 ensemble of heterogeneous RSF faults (rsf.h, 1500 m,
// 64 elements, +-20 % normal-stress heterogeneity, common random numbers), 2-D Theis forcing of K = 4 rate intervals of
// 100 d, mean(u) = 1 (fixed volume), u_k in [0, 2.5], window 400 d.
//
// Surrogate: Dieterich (1994) expected count at the fault-centre pressure, S = 0.6 dp, parameters A sigma, background
// stressing rate tau_dot, scale c: N_D = c * integral r / (gamma tau_dot) dt. Fitted by Poisson deviance to the RSF
// totals of T = 12 random training schedules (Nelder-Mead in log space, 8 starts). Optimised on a 0.1 simplex grid +
// coordinate refinement (the Dieterich count is convex in u, F1). Cost: 12 + 1 = 13 RSF ensemble evaluations.
//
// Competitor: 20 random schedules + 40 steps of a (1+1)-ES (Gaussian step, 1/5 success rule, projection onto the feasible
// set) = 60 RSF ensemble evaluations; 3 seeds; median best reported.
//
// H9b (argv[4] = m0): the same with the ensemble-mean seismic moment as objective (log least-squares fit).
// Verdict: SUPPORTED if N_RSF(surrogate optimum) <= 1.05 x median ES best; REFUTED if > 1.15 x; otherwise INCONCLUSIVE.
#include "ressim/rsf.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>
#include <vector>
using namespace ressim;

static const double day = 86400.0;
static const int K = 4;
static const double Tk = 100 * day, Pstar = 1.5e6, cdiff = 0.05, offset = 20.0, L = 1500.0, Umax = 2.5;
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

static int E = 16;
static std::vector<std::vector<double>> realisations;
static void makeRealisations(unsigned seed) {   // identical to H7b (seed 1007)
    std::mt19937 g(seed); std::uniform_real_distribution<double> U(0, 1);
    realisations.clear();
    for (int e = 0; e < E; ++e) {
        std::vector<double> sn(64);
        double ph[4], am[4];
        for (int m = 0; m < 4; ++m) { ph[m] = 2 * M_PI * U(g); am[m] = (U(g) - 0.5) * 0.1; }
        for (int i = 0; i < 64; ++i) {
            double d = 0; const double x = (i + 0.5) / 64;
            for (int m = 0; m < 4; ++m) d += am[m] * std::sin(2 * M_PI * (m + 1) * x + ph[m]);
            sn[size_t(i)] = 50e6 * (1 + d);
        }
        realisations.push_back(sn);
    }
}
static int rsfCalls = 0;
static double rsfSec = 0;
static bool g_m0 = false;   // H9b: objective = ensemble-mean seismic moment (N m) instead of the total count
static double rsfCount(const std::vector<double>& u) {   // total count over the ensemble (H9) or mean moment (H9b)
    ++rsfCalls;
    double n = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto& sn : realisations) {
        rsf::Input in; in.sigmaNProf = sn; in.n = 64; in.x1 = L; in.Vpl = 0; in.tEnd = K * Tk; in.vwFrom = 0.25; in.vwTo = 0.75; in.dtMax = 2 * day;
        const double dx = L / in.n;
        std::vector<double> r2(size_t(in.n));
        for (int i = 0; i < in.n; ++i) { const double s = (i + 0.5) * dx - 0.5 * L; r2[size_t(i)] = s * s + offset * offset; }
        in.forcing = [&](double t, std::vector<double>& dp, std::vector<double>&) { for (int i = 0; i < in.n; ++i) dp[size_t(i)] += dpAt(u, r2[size_t(i)], t); };
        rsf::Output o;
        if (!rsf::run(in, o)) std::printf("  (rsf run failed)\n");
        if (g_m0) { for (const auto& e : o.events) n += e.moment / E; } else n += double(o.events.size());
    }
    rsfSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return n;
}

// Dieterich expected count (unit scale) at the fault centre; pressure history precomputed per schedule
struct DCache { std::vector<double> S; };
static const int NT = 2000;
static DCache stressOf(const std::vector<double>& u) {
    DCache d; d.S.resize(NT + 1);
    const double r2 = offset * offset;   // fault centre
    for (int j = 0; j <= NT; ++j) d.S[size_t(j)] = 0.6 * dpAt(u, r2, j * (K * Tk) / NT);
    return d;
}
static double dieterich(const DCache& d, double Asig, double taudot) {
    const double dt = (K * Tk) / NT, r = 1.0 / (365.25 * day);
    double gam = 1.0 / taudot, cnt = 0;
    for (int j = 1; j <= NT; ++j) {
        const double dS = (d.S[size_t(j)] - d.S[size_t(j - 1)]) / dt + taudot, a = dS / Asig;
        if (std::fabs(a) > 1e-30) gam = (gam - 1 / (Asig * a)) * std::exp(-a * dt) + 1 / (Asig * a); else gam += dt / Asig;
        cnt += r / (gam * taudot) * dt;
    }
    return cnt;
}

static std::vector<double> project(std::vector<double> u) {   // clip to [0, Umax], mean 1 (iterated)
    for (int it = 0; it < 50; ++it) {
        double s = 0; for (double& v : u) { v = std::clamp(v, 0.0, Umax); s += v; }
        if (s <= 0) { u.assign(K, 1.0); return u; }
        for (double& v : u) v *= K / s;
        if (*std::max_element(u.begin(), u.end()) <= Umax + 1e-12) break;
    }
    for (double& v : u) v = std::clamp(v, 0.0, Umax);
    return u;
}

static std::vector<double> nelderMead(std::function<double(const std::vector<double>&)> f, std::vector<double> x0, double step, int iters) {
    const int n = int(x0.size());
    std::vector<std::vector<double>> P(size_t(n + 1), x0); std::vector<double> F(size_t(n + 1));
    for (int i = 0; i < n; ++i) P[size_t(i + 1)][size_t(i)] += step;
    for (int i = 0; i <= n; ++i) F[size_t(i)] = f(P[size_t(i)]);
    for (int it = 0; it < iters; ++it) {
        std::vector<int> o(size_t(n + 1)); for (int i = 0; i <= n; ++i) o[size_t(i)] = i;
        std::sort(o.begin(), o.end(), [&](int a, int b) { return F[size_t(a)] < F[size_t(b)]; });
        std::vector<std::vector<double>> Q; std::vector<double> G; for (int i : o) { Q.push_back(P[size_t(i)]); G.push_back(F[size_t(i)]); } P = Q; F = G;
        std::vector<double> c(size_t(n), 0.0); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) c[size_t(j)] += P[size_t(i)][size_t(j)] / n;
        auto at = [&](double a) { std::vector<double> x(static_cast<size_t>(n)); for (int j = 0; j < n; ++j) x[size_t(j)] = c[size_t(j)] + a * (P[size_t(n)][size_t(j)] - c[size_t(j)]); return x; };
        auto xr = at(-1); double fr = f(xr);
        if (fr < F[0]) { auto xe = at(-2); double fe = f(xe); if (fe < fr) { P[size_t(n)] = xe; F[size_t(n)] = fe; } else { P[size_t(n)] = xr; F[size_t(n)] = fr; } }
        else if (fr < F[size_t(n - 1)]) { P[size_t(n)] = xr; F[size_t(n)] = fr; }
        else { auto xc = at(0.5); double fc = f(xc);
            if (fc < F[size_t(n)]) { P[size_t(n)] = xc; F[size_t(n)] = fc; }
            else for (int i = 1; i <= n; ++i) { for (int j = 0; j < n; ++j) P[size_t(i)][size_t(j)] = P[0][size_t(j)] + 0.5 * (P[size_t(i)][size_t(j)] - P[0][size_t(j)]); F[size_t(i)] = f(P[size_t(i)]); } }
    }
    int b = int(std::min_element(F.begin(), F.end()) - F.begin());
    return P[size_t(b)];
}

// ---- C48 H11: Gaussian-process (squared-exponential, isotropic; ML length scale and noise, profiled signal variance)
// surrogate of log M0 with expected-improvement acquisition, total budget = training + EI steps RSF evaluations.
struct GP {
    std::vector<std::vector<double>> X; std::vector<double> y; double mu = 0, sf2 = 1, ell = 1, sn2 = 1e-4;
    std::vector<double> L, alpha; int n = 0;
    double k(const std::vector<double>& a, const std::vector<double>& b) const {
        double d2 = 0; for (size_t i = 0; i < a.size(); ++i) d2 += (a[i] - b[i]) * (a[i] - b[i]);
        return std::exp(-0.5 * d2 / (ell * ell));
    }
    // factor K/sf2 + (sn2/sf2) I; returns the negative profiled log likelihood (sf2 in closed form)
    double factor() {
        n = int(X.size()); L.assign(size_t(n * n), 0.0);
        for (int i = 0; i < n; ++i) for (int j = 0; j <= i; ++j) L[size_t(i * n + j)] = k(X[size_t(i)], X[size_t(j)]) + (i == j ? sn2 : 0.0);
        for (int j = 0; j < n; ++j) {
            double d = L[size_t(j * n + j)]; for (int q = 0; q < j; ++q) d -= L[size_t(j * n + q)] * L[size_t(j * n + q)];
            if (!(d > 0)) return 1e300;
            d = std::sqrt(d); L[size_t(j * n + j)] = d;
            for (int i = j + 1; i < n; ++i) { double v = L[size_t(i * n + j)]; for (int q = 0; q < j; ++q) v -= L[size_t(i * n + q)] * L[size_t(j * n + q)]; L[size_t(i * n + j)] = v / d; }
        }
        // mean = GLS constant, alpha = R^-1 (y - mu)
        std::vector<double> one(static_cast<size_t>(n), 1.0); auto a1 = solve(one), ay = solve(y);
        double s1 = 0, sy = 0; for (int i = 0; i < n; ++i) { s1 += a1[size_t(i)]; sy += ay[size_t(i)]; }
        mu = sy / s1;
        std::vector<double> r(static_cast<size_t>(n)); for (int i = 0; i < n; ++i) r[size_t(i)] = y[size_t(i)] - mu;
        alpha = solve(r);
        double q = 0; for (int i = 0; i < n; ++i) q += r[size_t(i)] * alpha[size_t(i)];
        sf2 = std::max(1e-12, q / n);
        double ld = 0; for (int i = 0; i < n; ++i) ld += std::log(L[size_t(i * n + i)]);
        return 0.5 * n * std::log(sf2) + ld;
    }
    std::vector<double> solve(const std::vector<double>& b) const {
        std::vector<double> z(b);
        for (int i = 0; i < n; ++i) { for (int q = 0; q < i; ++q) z[size_t(i)] -= L[size_t(i * n + q)] * z[size_t(q)]; z[size_t(i)] /= L[size_t(i * n + i)]; }
        for (int i = n - 1; i >= 0; --i) { for (int q = i + 1; q < n; ++q) z[size_t(i)] -= L[size_t(q * n + i)] * z[size_t(q)]; z[size_t(i)] /= L[size_t(i * n + i)]; }
        return z;
    }
    void fit() {   // ML over a log grid of the length scale and the relative noise
        double bestNll = 1e300, bl = 1, bs = 1e-4;
        for (double le = -1.5; le <= 1.01; le += 0.1) for (double ls = -6; ls <= -0.99; ls += 0.5) {
            ell = std::pow(10.0, le); sn2 = std::pow(10.0, ls);
            const double v = factor(); if (v < bestNll) { bestNll = v; bl = ell; bs = sn2; }
        }
        ell = bl; sn2 = bs; factor();
    }
    void predict(const std::vector<double>& x, double& m, double& sd) const {
        std::vector<double> kx(static_cast<size_t>(n)); for (int i = 0; i < n; ++i) kx[size_t(i)] = k(x, X[size_t(i)]);
        m = mu; for (int i = 0; i < n; ++i) m += kx[size_t(i)] * alpha[size_t(i)];
        auto v = solve(kx); double s = 1.0; for (int i = 0; i < n; ++i) s -= kx[size_t(i)] * v[size_t(i)];
        sd = std::sqrt(std::max(1e-14, sf2 * s));
    }
};

template <class RS>
static void gpDesign(std::vector<std::vector<double>> TU, std::vector<double> TN, int budget, RS& randomSchedule, std::mt19937& rng, double nConst) {
    GP gp;
    if (const char* e = std::getenv("H11_SEED")) rng.seed(unsigned(std::atoi(e)));   // exploratory: EI candidate sampler seed
    auto bestIdx = [&]() { return int(std::min_element(TN.begin(), TN.end()) - TN.begin()); };
    while (int(TU.size()) < budget) {
        gp.X = TU; gp.y.clear(); for (double v : TN) gp.y.push_back(std::log(v));
        gp.fit();
        const double fmin = *std::min_element(gp.y.begin(), gp.y.end());
        auto ei = [&](const std::vector<double>& u) {
            double m, sd; gp.predict(u, m, sd);
            const double z = (fmin - m) / sd;
            return sd * (z * 0.5 * std::erfc(-z / std::sqrt(2.0)) + std::exp(-0.5 * z * z) / std::sqrt(2 * 3.14159265358979323846));
        };
        std::vector<double> xb; double eb = -1;
        for (int c = 0; c < 20000; ++c) { auto u = randomSchedule(rng); const double e = ei(u); if (e > eb) { eb = e; xb = u; } }
        for (double h : {0.1, 0.05, 0.02}) for (int it = 0; it < 30; ++it) {   // pairwise moves on the simplex
            bool imp = false;
            for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) if (i != j) {
                auto u = xb; u[size_t(i)] += h; u[size_t(j)] -= h;
                if (u[size_t(i)] > Umax || u[size_t(j)] < 0) continue;
                const double e = ei(u); if (e > eb + 1e-15) { eb = e; xb = u; imp = true; }
            }
            if (!imp) break;
        }
        double m, sd; gp.predict(xb, m, sd);
        const double f = rsfCount(xb);
        TU.push_back(xb); TN.push_back(f);
        std::printf("GP-EI %2d ell %.3g noise %.2g | x [%.2f %.2f %.2f %.2f] EI %.3g pred %.4g (+- %.2f in log) RSF %.6g | best %.6g\n",
                    int(TU.size()), gp.ell, gp.sn2, xb[0], xb[1], xb[2], xb[3], eb, std::exp(m), sd, f, TN[size_t(bestIdx())]);
        std::fflush(stdout);
    }
    const int b = bestIdx(); const double esMed = 2.21503e15, ratio = TN[size_t(b)] / esMed;
    std::printf("SUMMARY H11 const %.6g | GP-EI best %.6g at [%.2f %.2f %.2f %.2f] (%d RSF evaluations incl. 12 training) | ES median best (H9b) %.6g | ratio %.3f | vs const %.3f | RSF %.0f s\n",
                nConst, TN[size_t(b)], TU[size_t(b)][0], TU[size_t(b)][1], TU[size_t(b)][2], TU[size_t(b)][3], int(TU.size()), esMed, ratio, TN[size_t(b)] / nConst, rsfSec);
    std::printf("VERDICT H11 %s\n", ratio <= 1.05 ? "SUPPORTED (<= 1.05)" : ratio > 1.15 ? "REFUTED (> 1.15)" : "INCONCLUSIVE (1.05-1.15)");
}

int main(int argc, char** argv) {
    const int T = argc > 1 ? std::atoi(argv[1]) : 12;
    const int budget = argc > 2 ? std::atoi(argv[2]) : 60;
    const int seeds = argc > 3 ? std::atoi(argv[3]) : 3;
    const bool gpMode = argc > 4 && std::string(argv[4]) == "m0gp";   // C48 H11: GP + expected improvement on log M0
    g_m0 = argc > 4 && (std::string(argv[4]) == "m0" || gpMode);
    makeRealisations(7 + 1000);
    std::mt19937 rng(9);
    if (const char* e = std::getenv("H11_TRAIN_SEED")) rng.seed(unsigned(std::atoi(e)));   // C49 H11b: training-set resampling
    std::uniform_real_distribution<double> U(0, 1);
    auto randomSchedule = [&](std::mt19937& g) {
        std::vector<double> u(K);
        for (;;) { double s = 0; for (double& v : u) { v = U(g); s += v; } for (double& v : u) v *= K / s; if (*std::max_element(u.begin(), u.end()) <= Umax) return u; }
    };
    const std::vector<double> ones(K, 1.0);
    const double nConst = rsfCount(ones);
    std::printf("CONST u=1: RSF objective %.6g (%s)\n", nConst, g_m0 ? "mean moment N m" : "total count");

    // ---- surrogate: training
    std::vector<std::vector<double>> TU; std::vector<double> TN; std::vector<DCache> TD;
    for (int t = 0; t < T; ++t) {
        auto u = randomSchedule(rng); const double n = rsfCount(u);
        TU.push_back(u); TN.push_back(n); TD.push_back(stressOf(u));
        std::printf("TRAIN %2d [%.2f %.2f %.2f %.2f] RSF %.6g\n", t, u[0], u[1], u[2], u[3], n); std::fflush(stdout);
    }
    if (gpMode) { gpDesign(TU, TN, budget, randomSchedule, rng, nConst); return 0; }
    // Poisson deviance in log parameters (log10 A sigma [Pa], log10 tau_dot [Pa/s]); the scale c is optimal in closed form
    auto fitObj = [&](const std::vector<double>& p, double* cOut) {
        const double Asig = std::pow(10.0, p[0]), td = std::pow(10.0, p[1]);
        std::vector<double> m(static_cast<size_t>(T)); double sm = 0, sn = 0;
        for (int t = 0; t < T; ++t) { m[size_t(t)] = dieterich(TD[size_t(t)], Asig, td); sm += m[size_t(t)]; sn += TN[size_t(t)]; }
        if (!(sm > 0) || !std::isfinite(sm)) return 1e300;
        double c = sn / sm;
        if (g_m0) {   // H9b: least squares on log M0, closed-form scale
            double lc = 0; int k = 0;
            for (int t = 0; t < T; ++t) if (TN[size_t(t)] > 0 && m[size_t(t)] > 0) { lc += std::log(TN[size_t(t)] / m[size_t(t)]); ++k; }
            c = std::exp(lc / std::max(1, k)); if (cOut) *cOut = c;
            double e2 = 0;
            for (int t = 0; t < T; ++t) if (TN[size_t(t)] > 0 && m[size_t(t)] > 0) { const double r = std::log(c * m[size_t(t)] / TN[size_t(t)]); e2 += r * r; }
            return e2;
        }
        if (cOut) *cOut = c;
        double dev = 0;
        for (int t = 0; t < T; ++t) { const double mu = std::max(1e-12, c * m[size_t(t)]); const double y = TN[size_t(t)]; dev += 2 * ((y > 0 ? y * std::log(y / mu) : 0) - (y - mu)); }
        return dev;
    };
    std::vector<double> best = {5.7, -4.5}; double bestDev = 1e300;
    for (double a0 : {4.5, 5.0, 5.5, 6.0, 6.5, 7.0}) for (double t0 : {-6.5, -5.0, -3.5}) {
        auto p = nelderMead([&](const std::vector<double>& x) { return fitObj(x, nullptr); }, {a0, t0}, 0.4, 150);
        const double d = fitObj(p, nullptr); if (d < bestDev) { bestDev = d; best = p; }
    }
    double cFit = 1; fitObj(best, &cFit);
    const double AsigF = std::pow(10.0, best[0]), tdF = std::pow(10.0, best[1]);
    std::printf("FIT A_sigma %.4g Pa  tau_dot %.4g Pa/s (t_a %.3g d)  scale %.4g  deviance %.3f (%d schedules)\n", AsigF, tdF, AsigF / tdF / day, cFit, bestDev, T);
    double ssr = 0, sst = 0, mean = 0; for (double n : TN) mean += n; mean /= T;
    for (int t = 0; t < T; ++t) { const double pr = cFit * dieterich(TD[size_t(t)], AsigF, tdF); ssr += (pr - TN[size_t(t)]) * (pr - TN[size_t(t)]); sst += (TN[size_t(t)] - mean) * (TN[size_t(t)] - mean); }
    std::printf("FIT R2 on training %.3f\n", sst > 0 ? 1 - ssr / sst : 0.0);
    // optimise the surrogate: 0.1 grid on the simplex + coordinate refinement
    auto surr = [&](const std::vector<double>& u) { return cFit * dieterich(stressOf(u), AsigF, tdF); };
    std::vector<double> uS = ones; double fS = surr(ones);
    for (int a = 0; a <= 25; ++a) for (int b = 0; b <= 25; ++b) for (int c = 0; c <= 25; ++c) {
        const double d = 4.0 - 0.1 * (a + b + c); if (d < -1e-9 || d > Umax + 1e-9) continue;
        std::vector<double> u = {0.1 * a, 0.1 * b, 0.1 * c, std::max(0.0, d)};
        const double f = surr(u); if (f < fS) { fS = f; uS = u; }
    }
    for (double h : {0.05, 0.02, 0.01}) for (int it = 0; it < 30; ++it) {
        bool imp = false;
        for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) if (i != j) {
            auto u = uS; u[size_t(i)] += h; u[size_t(j)] -= h;
            if (u[size_t(i)] > Umax || u[size_t(j)] < 0) continue;
            const double f = surr(u); if (f < fS - 1e-12) { fS = f; uS = u; imp = true; }
        }
        if (!imp) break;
    }
    const double nS = rsfCount(uS);
    std::printf("SURROGATE optimum [%.2f %.2f %.2f %.2f] predicted %.4g, RSF %.6g (%d RSF evaluations incl. training)\n", uS[0], uS[1], uS[2], uS[3], fS, nS, T + 1);
    std::fflush(stdout);

    // ---- stochastic competitor
    std::vector<double> bestES;
    for (int s = 0; s < seeds; ++s) {
        std::mt19937 g(100 + s); std::normal_distribution<double> N01(0, 1);
        std::vector<double> xb; double fb = 1e300; int used = 0;
        const int nRand = std::min(20, budget);
        for (; used < nRand; ++used) { auto u = randomSchedule(g); const double f = rsfCount(u); if (f < fb) { fb = f; xb = u; } }
        double sigma = 0.4; int succ = 0, tries = 0;
        for (; used < budget; ++used) {
            auto u = xb; for (double& v : u) v += sigma * N01(g); u = project(u);
            const double f = rsfCount(u); ++tries;
            if (f <= fb) { if (f < fb) ++succ; fb = f; xb = u; }
            if (tries == 5) { sigma *= (succ >= 1) ? 1.5 : 0.6; sigma = std::clamp(sigma, 0.02, 1.0); tries = succ = 0; }
        }
        std::printf("ES seed %d best RSF %.6g at [%.2f %.2f %.2f %.2f] (%d evaluations)\n", s, fb, xb[0], xb[1], xb[2], xb[3], budget); std::fflush(stdout);
        bestES.push_back(fb);
    }
    std::sort(bestES.begin(), bestES.end());
    const double med = bestES[bestES.size() / 2];
    const double ratio = nS / med;
    std::printf("SUMMARY const %.6g | surrogate %.6g (13 evals) | ES median best %.6g (%d evals x %d seeds) | ratio %.3f | RSF calls %d, %.0f s\n",
                nConst, nS, med, budget, seeds, ratio, rsfCalls, rsfSec);
    std::printf("VERDICT %s\n", ratio <= 1.05 ? "SUPPORTED (<= 1.05)" : ratio > 1.15 ? "REFUTED (> 1.15)" : "INCONCLUSIVE (1.05-1.15)");
}
