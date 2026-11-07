#pragma once
// Nonlinear-solver toolkit shared by the reservoir engine and future physics modules (hysteresis, compositional
// flash, thermal, geomechanics, coupled multiphysics).
//
//  * Dual<N>        forward-mode dual numbers (alias of AD<N>), plus the elementwise functions physics code needs
//                   (abs, min/max, smooth max, pow with dual exponent, comparisons on values).
//  * jacobian()     ForwardDiff-style: evaluate a vector function written once for generic scalars and get its value
//                   and exact Jacobian (Julia's ForwardDiff.jacobian idiom); derivative() for scalar functions.
//  * localNewton()  damped Newton with exact AD Jacobian, pivoted LU and Armijo backtracking line search for small
//                   dense systems (phase equilibrium / flash, well equations, thermal equilibrium, local constitutive
//                   updates). Box bounds keep iterates physical.
//  * OscillationDetector  OPM-style detection of oscillating / stagnating Newton residuals with update relaxation.
//  * PidStepController    PID time-step control on a relative change measure (Valli / Gustafsson type).
#include "ressim/ad.h"
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <algorithm>

namespace ressim {

template <int N> using Dual = AD<N>;

// ---------------------------------------------------------------- elementwise functions
template <int N> inline AD<N> abs(const AD<N>& x) { return x.v >= 0 ? x : -x; }
template <int N> inline AD<N> max(const AD<N>& a, const AD<N>& b) { return a.v >= b.v ? a : b; }
template <int N> inline AD<N> min(const AD<N>& a, const AD<N>& b) { return a.v <= b.v ? a : b; }
template <int N> inline AD<N> max(const AD<N>& a, double b) { return a.v >= b ? a : AD<N>(b); }
template <int N> inline AD<N> min(const AD<N>& a, double b) { return a.v <= b ? a : AD<N>(b); }
// C1 approximation of max(a, b) with transition width eps (keeps Newton quadratic near kinks)
template <int N> inline AD<N> smoothMax(const AD<N>& a, const AD<N>& b, double eps) {
    const AD<N> d = a - b;
    return 0.5 * (a + b + sqrt(d * d + eps * eps));
}
template <int N> inline AD<N> pow(const AD<N>& x, const AD<N>& y) { return exp(y * log(x)); }
template <int N> inline AD<N> tanh(const AD<N>& x) { const double t = std::tanh(x.v); return chain(x, t, 1.0 - t * t); }
template <int N> inline bool operator<(const AD<N>& a, const AD<N>& b) { return a.v < b.v; }
template <int N> inline bool operator>(const AD<N>& a, const AD<N>& b) { return a.v > b.v; }
template <int N> inline bool operator<(const AD<N>& a, double b) { return a.v < b; }
template <int N> inline bool operator>(const AD<N>& a, double b) { return a.v > b; }
inline double value(double x) { return x; }
template <int N> inline double value(const AD<N>& x) { return x.v; }

// ---------------------------------------------------------------- ForwardDiff-style drivers
// F: callable  std::array<AD<N>, M> F(const std::array<AD<N>, N>&)  (write it as a template / generic lambda)
template <int N, int M, class F>
void jacobian(F&& f, const std::array<double, N>& x, std::array<double, M>& fx, std::array<std::array<double, N>, M>& J) {
    std::array<AD<N>, N> xa;
    for (int i = 0; i < N; ++i) xa[i] = AD<N>::variable(x[i], i);
    const std::array<AD<N>, M> r = f(xa);
    for (int k = 0; k < M; ++k) { fx[k] = r[k].v; for (int i = 0; i < N; ++i) J[k][i] = r[k].d[i]; }
}
template <class F> inline double derivative(F&& f, double x, double* fx = nullptr) {
    const AD<1> r = f(AD<1>::variable(x, 0));
    if (fx) *fx = r.v;
    return r.d[0];
}

// dense pivoted LU solve of J dx = b (in place; false if singular)
template <int N> inline bool denseSolve(std::array<std::array<double, N>, N> J, std::array<double, N>& b) {
    double sc = 0; for (auto& r : J) for (double v : r) sc = std::max(sc, std::fabs(v));
    if (!(sc > 0) || !std::isfinite(sc)) return false;
    for (int k = 0; k < N; ++k) {
        int p = k;
        for (int i = k + 1; i < N; ++i) if (std::fabs(J[i][k]) > std::fabs(J[p][k])) p = i;
        if (std::fabs(J[p][k]) <= 1e-14 * sc) return false;
        std::swap(J[p], J[k]); std::swap(b[p], b[k]);
        for (int i = k + 1; i < N; ++i) {
            const double f = J[i][k] / J[k][k];
            if (f == 0) continue;
            for (int j = k; j < N; ++j) J[i][j] -= f * J[k][j];
            b[i] -= f * b[k];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < N; ++j) s -= J[i][j] * b[j];
        b[i] = s / J[i][i];
    }
    return true;
}

struct LocalNewtonOptions {
    int maxIt = 50;
    double tol = 1e-10;            // on max |F|
    double stepTol = 1e-14;        // on max |dx| (relative to max(1,|x|))
    int maxBacktrack = 12;
    double armijo = 1e-4;
};
struct LocalNewtonResult { bool converged = false; int iterations = 0; double residual = 0; };

// Solve F(x) = 0, lo <= x <= hi.  F as for jacobian() with M == N.
template <int N, class F>
LocalNewtonResult localNewton(F&& f, std::array<double, N>& x, const LocalNewtonOptions& o = {},
                              const std::array<double, N>* lo = nullptr, const std::array<double, N>* hi = nullptr) {
    LocalNewtonResult res;
    auto clampX = [&](std::array<double, N>& y) {
        for (int i = 0; i < N; ++i) {
            if (lo) y[i] = std::max(y[i], (*lo)[i]);
            if (hi) y[i] = std::min(y[i], (*hi)[i]);
        }
    };
    auto merit = [](const std::array<double, N>& r) { double s = 0; for (double v : r) s += v * v; return 0.5 * s; };
    auto inf = [](const std::array<double, N>& r) { double s = 0; for (double v : r) s = std::max(s, std::fabs(v)); return s; };
    clampX(x);
    std::array<double, N> fx; std::array<std::array<double, N>, N> J;
    jacobian<N, N>(f, x, fx, J);
    for (res.iterations = 0; res.iterations < o.maxIt; ++res.iterations) {
        res.residual = inf(fx);
        if (res.residual < o.tol) { res.converged = true; return res; }
        std::array<double, N> dx = fx;
        if (!denseSolve<N>(J, dx)) {                 // singular: steepest descent direction J^T F
            for (int i = 0; i < N; ++i) { double s = 0; for (int k = 0; k < N; ++k) s += J[k][i] * fx[k]; dx[i] = s; }
        }
        const double m0 = merit(fx);
        double lam = 1.0;
        std::array<double, N> xn, fn; std::array<std::array<double, N>, N> Jn;
        bool accepted = false;
        for (int bt = 0; bt <= o.maxBacktrack; ++bt, lam *= 0.5) {
            for (int i = 0; i < N; ++i) xn[i] = x[i] - lam * dx[i];
            clampX(xn);
            jacobian<N, N>(f, xn, fn, Jn);
            bool fin = true; for (double v : fn) fin = fin && std::isfinite(v);
            if (fin && merit(fn) <= (1.0 - 2.0 * o.armijo * lam) * m0) { accepted = true; break; }
        }
        double step = 0;
        for (int i = 0; i < N; ++i) step = std::max(step, std::fabs(xn[i] - x[i]) / std::max(1.0, std::fabs(x[i])));
        if (!accepted) {                              // take the smallest step anyway if it is finite
            bool fin = true; for (double v : fn) fin = fin && std::isfinite(v);
            if (!fin) return res;
        }
        x = xn; fx = fn; J = Jn;
        if (step < o.stepTol) { res.residual = inf(fx); res.converged = res.residual < o.tol; return res; }
    }
    res.residual = inf(fx);
    res.converged = res.residual < o.tol;
    return res;
}

// ---------------------------------------------------------------- Newton stabilisation
// Tracks per-equation residual norms over Newton iterations. Oscillation: the residual of at least `minEqs`
// equations returned close to its value two iterations ago while the last change was large (OPM criterion).
// The caller multiplies its Newton update by omega().
class OscillationDetector {
public:
    explicit OscillationDetector(int nEq = 3, double relTol = 0.2, double omegaMin = 0.5, double dOmega = 0.1, int minEqs = 2)
        : nEq_(nEq), relTol_(relTol), omegaMin_(omegaMin), dOmega_(dOmega), minEqs_(std::min(minEqs, nEq)) {}
    void reset() { hist_.clear(); omega_ = 1.0; oscillating_ = stagnating_ = false; }
    // r: nEq residual norms of the current iteration
    void push(const double* r) {
        std::array<double, 8> v{}; for (int i = 0; i < nEq_ && i < 8; ++i) v[i] = r[i];
        hist_.push_back(v); if (hist_.size() > 3) hist_.pop_front();
        oscillating_ = stagnating_ = false;
        if (hist_.size() < 3) return;
        int osc = 0, stag = 0, active = 0;
        for (int i = 0; i < nEq_ && i < 8; ++i) {
            const double r0 = hist_[2][i], r1 = hist_[1][i], r2 = hist_[0][i];
            if (!(r0 > 0)) continue;
            ++active;
            const double d1 = std::fabs((r0 - r2) / r0), d2 = std::fabs((r0 - r1) / r0);
            if (d1 < relTol_ && relTol_ < d2) ++osc;
            if (d2 < 1e-3) ++stag;
        }
        oscillating_ = osc >= minEqs_;
        stagnating_ = active > 0 && stag == active;
        if (oscillating_ || stagnating_) omega_ = std::max(omega_ - dOmega_, omegaMin_);
    }
    double omega() const { return omega_; }
    bool oscillating() const { return oscillating_; }
    bool stagnating() const { return stagnating_; }
private:
    int nEq_; double relTol_, omegaMin_, dOmega_; int minEqs_;
    std::deque<std::array<double, 8>> hist_;
    double omega_ = 1.0; bool oscillating_ = false, stagnating_ = false;
};

// PID time-step control: e_n = max relative change of the primary variables over the step divided by its target.
// dt_{n+1} = dt_n * (e_{n-1}/e_n)^kP * (1/e_n)^kI * (e_{n-1}^2 / (e_n e_{n-2}))^kD, limited to [minFactor, maxFactor].
class PidStepController {
public:
    PidStepController(double kP = 0.075, double kI = 0.175, double kD = 0.01, double minFactor = 0.2, double maxFactor = 3.0)
        : kP_(kP), kI_(kI), kD_(kD), minF_(minFactor), maxF_(maxFactor) {}
    void reset() { e1_ = e2_ = 1.0; }
    double factor(double e) {
        e = std::max(e, 1e-10);
        double f = std::pow(e1_ / e, kP_) * std::pow(1.0 / e, kI_) * std::pow(e1_ * e1_ / (e * e2_), kD_);
        if (!std::isfinite(f)) f = 1.0;
        e2_ = e1_; e1_ = e;
        return std::min(maxF_, std::max(minF_, f));
    }
private:
    double kP_, kI_, kD_, minF_, maxF_;
    double e1_ = 1.0, e2_ = 1.0;
};

} // namespace ressim
