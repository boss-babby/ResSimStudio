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

template <int N> inline AD<N> pow(const AD<N>& x, const AD<N>& y) {
    // x^y = exp(y ln x)
    const double f = std::pow(x.v, y.v);
    return chain2(x, y, f, x.v != 0 ? y.v * f / x.v : 0.0, x.v > 0 ? f * std::log(x.v) : 0.0);
}
template <int N> inline AD<N> abs(const AD<N>& x) { return x.v < 0 ? -x : x; }
template <int N> inline AD<N> tanh(const AD<N>& x) { const double t = std::tanh(x.v); return chain(x, t, 1.0 - t * t); }
template <int N> inline AD<N> sin(const AD<N>& x) { return chain(x, std::sin(x.v), std::cos(x.v)); }
template <int N> inline AD<N> cos(const AD<N>& x) { return chain(x, std::cos(x.v), -std::sin(x.v)); }
template <int N> inline AD<N> atan(const AD<N>& x) { return chain(x, std::atan(x.v), 1.0 / (1.0 + x.v * x.v)); }
template <int N> inline AD<N> cbrt(const AD<N>& x) { const double c = std::cbrt(x.v); return chain(x, c, c != 0 ? c / (3.0 * x.v) : 0.0); }

// branch selection keeps the derivative of the chosen argument (sub-gradient at ties: first argument)
template <int N> inline AD<N> max(const AD<N>& a, const AD<N>& b) { return a.v >= b.v ? a : b; }
template <int N> inline AD<N> min(const AD<N>& a, const AD<N>& b) { return a.v <= b.v ? a : b; }
template <int N> inline AD<N> max(const AD<N>& a, double b) { return a.v >= b ? a : AD<N>(b); }
template <int N> inline AD<N> min(const AD<N>& a, double b) { return a.v <= b ? a : AD<N>(b); }

// smooth max (C-infinity, |error| <= eps * ln 2): keeps Newton quadratic near switching points
template <int N> inline AD<N> smoothMax(const AD<N>& a, const AD<N>& b, double eps) {
    const AD<N> d = a - b;
    return 0.5 * (a + b + sqrt(d * d + eps * eps));
}

template <int N> inline bool operator<(const AD<N>& a, const AD<N>& b) { return a.v < b.v; }
template <int N> inline bool operator>(const AD<N>& a, const AD<N>& b) { return a.v > b.v; }
template <int N> inline bool operator<(const AD<N>& a, double b) { return a.v < b; }
template <int N> inline bool operator>(const AD<N>& a, double b) { return a.v > b; }
template <int N> inline double value(const AD<N>& a) { return a.v; }
inline double value(double a) { return a; }

// ---------------------------------------------------------------- runtime-sized AD
// For models whose number of unknowns per cell is only known at run time (compositional: nc components, thermal:
// + temperature, ...). Fixed inline capacity of kADVInline derivatives (no heap allocation in the assembly loops).
constexpr int kADVInline = 12;
struct ADV {
    double v = 0.0;
    int n = 0;
    double d[kADVInline] = {};

    ADV() = default;
    ADV(double value, int nder = 0) : v(value), n(nder) {}
    static ADV variable(double value, int nder, int slot) {
        ADV a(value, nder < kADVInline ? nder : kADVInline);
        if (slot < kADVInline) a.d[slot] = 1.0;
        return a;
    }
    int size() const { return n; }

    ADV& operator+=(const ADV& b) {
        v += b.v;
        grow(b.n);
        for (int i = 0; i < b.n; ++i) d[i] += b.d[i];
        return *this;
    }
    ADV& operator-=(const ADV& b) {
        v -= b.v;
        grow(b.n);
        for (int i = 0; i < b.n; ++i) d[i] -= b.d[i];
        return *this;
    }
    ADV& operator*=(const ADV& b) {
        grow(b.n);
        for (int i = 0; i < n; ++i) d[i] = d[i] * b.v + v * (i < b.n ? b.d[i] : 0.0);
        v *= b.v;
        return *this;
    }
    ADV& operator/=(const ADV& b) {
        grow(b.n);
        const double inv = 1.0 / b.v;
        for (int i = 0; i < n; ++i) d[i] = (d[i] * b.v - v * (i < b.n ? b.d[i] : 0.0)) * inv * inv;
        v *= inv;
        return *this;
    }
    ADV& operator*=(double b) { v *= b; for (int i = 0; i < n; ++i) d[i] *= b; return *this; }
    ADV& operator+=(double b) { v += b; return *this; }
    ADV& operator-=(double b) { v -= b; return *this; }

    void grow(int m) {
        if (m > kADVInline) m = kADVInline;
        if (m > n) {
            for (int i = n; i < m; ++i) d[i] = 0.0;
            n = m;
        }
    }
};

inline ADV operator+(ADV a, const ADV& b) { return a += b; }
inline ADV operator-(ADV a, const ADV& b) { return a -= b; }
inline ADV operator*(ADV a, const ADV& b) { return a *= b; }
inline ADV operator/(ADV a, const ADV& b) { return a /= b; }
inline ADV operator*(ADV a, double b) { return a *= b; }
inline ADV operator*(double b, ADV a) { return a *= b; }
inline ADV operator+(ADV a, double b) { return a += b; }
inline ADV operator-(ADV a, double b) { return a -= b; }
inline ADV operator-(ADV a) { a.v = -a.v; for (int i = 0; i < a.n; ++i) a.d[i] = -a.d[i]; return a; }

inline ADV chainV(const ADV& x, double f, double df) {
    ADV r(f, x.n);
    for (int i = 0; i < x.n; ++i) r.d[i] = df * x.d[i];
    return r;
}
inline ADV exp(const ADV& x) { const double e = std::exp(x.v); return chainV(x, e, e); }
inline ADV log(const ADV& x) { return chainV(x, std::log(x.v), 1.0 / x.v); }
inline ADV sqrt(const ADV& x) { const double s = std::sqrt(x.v); return chainV(x, s, s > 0 ? 0.5 / s : 0.0); }
inline ADV pow(const ADV& x, double p) { const double f = std::pow(x.v, p); return chainV(x, f, x.v != 0 ? p * f / x.v : 0.0); }

// Re-map derivatives of a small AD into a larger one (e.g. cell vars -> slot offset in a 2-cell flux).
template <int M, int N> inline AD<M> lift(const AD<N>& a, int offset) {
    AD<M> r(a.v); for (int i = 0; i < N; ++i) r.d[offset + i] = a.d[i]; return r;
}

using AD3 = AD<3>;
using AD6 = AD<6>;

} // namespace ressim
