// Forward-mode automatic differentiation with a fixed number of derivatives.
// Used to build the exact Jacobian of the fully-implicit black-oil equations.
#pragma once
#include <cmath>
#include <algorithm>

namespace ressim {

template <int N>
struct AD {
    double v = 0.0;
    double d[N] = {};

    AD() = default;
    AD(double value) : v(value) {}                       // constant (zero derivatives)
    static AD variable(double value, int slot) { AD a(value); a.d[slot] = 1.0; return a; }

    AD& operator+=(const AD& b) { v += b.v; for (int i = 0; i < N; ++i) d[i] += b.d[i]; return *this; }
    AD& operator-=(const AD& b) { v -= b.v; for (int i = 0; i < N; ++i) d[i] -= b.d[i]; return *this; }
    AD& operator*=(const AD& b) { for (int i = 0; i < N; ++i) d[i] = d[i] * b.v + v * b.d[i]; v *= b.v; return *this; }
    AD& operator/=(const AD& b) {
        const double inv = 1.0 / b.v;
        for (int i = 0; i < N; ++i) d[i] = (d[i] * b.v - v * b.d[i]) * inv * inv;
        v *= inv; return *this;
    }
    AD& operator+=(double b) { v += b; return *this; }
    AD& operator-=(double b) { v -= b; return *this; }
    AD& operator*=(double b) { v *= b; for (int i = 0; i < N; ++i) d[i] *= b; return *this; }
    AD& operator/=(double b) { return (*this) *= (1.0 / b); }
};

template <int N> inline AD<N> operator-(AD<N> a) { a.v = -a.v; for (int i = 0; i < N; ++i) a.d[i] = -a.d[i]; return a; }
template <int N> inline AD<N> operator+(AD<N> a, const AD<N>& b) { return a += b; }
template <int N> inline AD<N> operator-(AD<N> a, const AD<N>& b) { return a -= b; }
template <int N> inline AD<N> operator*(AD<N> a, const AD<N>& b) { return a *= b; }
template <int N> inline AD<N> operator/(AD<N> a, const AD<N>& b) { return a /= b; }
template <int N> inline AD<N> operator+(AD<N> a, double b) { return a += b; }
template <int N> inline AD<N> operator-(AD<N> a, double b) { return a -= b; }
template <int N> inline AD<N> operator*(AD<N> a, double b) { return a *= b; }
template <int N> inline AD<N> operator/(AD<N> a, double b) { return a /= b; }
template <int N> inline AD<N> operator+(double a, AD<N> b) { return b += a; }
template <int N> inline AD<N> operator-(double a, const AD<N>& b) { return (-b) += a; }
template <int N> inline AD<N> operator*(double a, AD<N> b) { return b *= a; }
template <int N> inline AD<N> operator/(double a, const AD<N>& b) {
    AD<N> r(a / b.v); const double f = -a / (b.v * b.v);
    for (int i = 0; i < N; ++i) r.d[i] = f * b.d[i];
    return r;
}

// Chain rule helper: value f(x) with derivative df/dx.
template <int N> inline AD<N> chain(const AD<N>& x, double f, double dfdx) {
    AD<N> r(f); for (int i = 0; i < N; ++i) r.d[i] = dfdx * x.d[i]; return r;
}
// Two-argument chain rule: f(x,y).
template <int N> inline AD<N> chain2(const AD<N>& x, const AD<N>& y, double f, double fx, double fy) {
    AD<N> r(f); for (int i = 0; i < N; ++i) r.d[i] = fx * x.d[i] + fy * y.d[i]; return r;
}

template <int N> inline AD<N> sqrt(const AD<N>& x) { const double s = std::sqrt(x.v); return chain(x, s, s > 0 ? 0.5 / s : 0.0); }
template <int N> inline AD<N> exp(const AD<N>& x) { const double e = std::exp(x.v); return chain(x, e, e); }
template <int N> inline AD<N> log(const AD<N>& x) { return chain(x, std::log(x.v), 1.0 / x.v); }
template <int N> inline AD<N> pow(const AD<N>& x, double p) { const double f = std::pow(x.v, p); return chain(x, f, x.v != 0 ? p * f / x.v : 0.0); }

// Re-map derivatives of a small AD into a larger one (e.g. cell vars -> slot offset in a 2-cell flux).
template <int M, int N> inline AD<M> lift(const AD<N>& a, int offset) {
    AD<M> r(a.v); for (int i = 0; i < N; ++i) r.d[offset + i] = a.d[i]; return r;
}

using AD3 = AD<3>;
using AD6 = AD<6>;

} // namespace ressim
