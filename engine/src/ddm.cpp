// C26: 3-D displacement discontinuity method on triangular dislocation elements - see ddm.h
#include "ressim/ddm.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <omp.h>

namespace ressim {
namespace ddm {
namespace tde {
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#define WITHIN_KERNEL static inline
#include "tde_kernels.inc"
#undef WITHIN_KERNEL
#undef Real
#undef EPS
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
}  // namespace tde

namespace {
V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 cross(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double len(const V3& a) { return std::sqrt(dot(a, a)); }
V3 scale(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
bool finite6(const double* x, int n) { for (int i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false; return true; }
void pack(const Tri& t, double* v) { for (int d = 0; d < 3; ++d) { v[d] = t.a[size_t(d)]; v[3 + d] = t.b[size_t(d)]; v[6 + d] = t.c[size_t(d)]; } }
// traction of a symmetric tensor (xx yy zz xy xz yz) on n
V3 traction(const double S[6], const V3& n) {
    return {S[0] * n[0] + S[3] * n[1] + S[4] * n[2], S[3] * n[0] + S[1] * n[1] + S[5] * n[2], S[4] * n[0] + S[5] * n[1] + S[2] * n[2]};
}
}  // namespace

void triFrame(const Tri& t, V3& n, V3& s, V3& d) {
    n = cross(sub(t.b, t.a), sub(t.c, t.a));
    const double ln = len(n);
    n = ln > 0 ? scale(n, 1 / ln) : V3{0, 0, 1};
    s = cross(V3{0, 0, 1}, n);
    if (len(s) == 0) s = scale(V3{0, 1, 0}, n[2]);
    s = scale(s, 1 / len(s));
    d = cross(n, s);
}

bool triStrain(const V3& obs, const Tri& t, const V3& slip, double nu, bool hs, double e[6]) {
    double o[3] = {obs[0], obs[1], obs[2]}, v[9], s[3] = {slip[0], slip[1], slip[2]};
    pack(t, v);
    if (hs) tde::tde_strain_hs(o, v, s, nu, e); else tde::tde_strain_fs(o, v, s, nu, e);
    return finite6(e, 6);
}

bool triDisp(const V3& obs, const Tri& t, const V3& slip, double nu, bool hs, double u[3]) {
    double o[3] = {obs[0], obs[1], obs[2]}, v[9], s[3] = {slip[0], slip[1], slip[2]};
    pack(t, v);
    if (hs) tde::tde_disp_hs(o, v, s, nu, u); else tde::tde_disp_fs(o, v, s, nu, u);
    return finite6(u, 3);
}

bool triStress(const V3& obs, const Tri& t, const V3& slip, double G, double nu, bool hs, double S[6]) {
    double e[6];
    if (!triStrain(obs, t, slip, nu, hs, e)) { for (int i = 0; i < 6; ++i) S[i] = 0; return false; }
    const double lam = 2 * G * nu / (1 - 2 * nu), tr = e[0] + e[1] + e[2];
    S[0] = lam * tr + 2 * G * e[0]; S[1] = lam * tr + 2 * G * e[1]; S[2] = lam * tr + 2 * G * e[2];
    S[3] = 2 * G * e[3]; S[4] = 2 * G * e[4]; S[5] = 2 * G * e[5];
    return true;
}

Element quadElement(const V3& p0, const V3& p1, const V3& p2, const V3& p3) {
    Element e;
    const V3 nv = cross(sub(p2, p0), sub(p3, p1));
    const double ln = len(nv);
    e.area = 0.5 * ln;
    e.n = ln > 0 ? scale(nv, 1 / ln) : V3{0, 0, 1};
    for (int d = 0; d < 3; ++d) e.c[size_t(d)] = 0.25 * (p0[size_t(d)] + p1[size_t(d)] + p2[size_t(d)] + p3[size_t(d)]);
    auto add = [&](const V3& a, const V3& b, const V3& c) {
        Tri t{a, b, c};
        const V3 tn = cross(sub(b, a), sub(c, a));
        if (len(tn) <= 1e-12 * std::max(1e-30, e.area)) return;   // degenerate (collapsed corner-point face)
        if (dot(tn, e.n) < 0) std::swap(t.b, t.c);
        e.tris.push_back(t);
    };
    add(p0, p1, p2); add(p0, p2, p3);
    e.size = std::sqrt(std::max(0.0, e.area));
    return e;
}

bool elementStress(const V3& obs, const Element& e, const V3& b, double G, double nu, bool hs, double S[6]) {
    // The receiver is evaluated at +-delta along the SOURCE normal and averaged (delta = 1e-3 of the element size): the
    // tractions on planes parallel to the source are continuous across it, the in-plane components are averaged (mean of
    // the jump), and points on (or within rounding of) a triangle edge line - the centroid of a quad lies on its diagonal,
    // grid neighbours lie on the extension of each other's diagonals - never reach the singular TDE configuration, where
    // the 1/r edge terms of the two triangles would otherwise cancel catastrophically. Error O((delta / r)^2).
    for (int i = 0; i < 6; ++i) S[i] = 0;
    V3 t1 = sub(e.tris.empty() ? V3{1, 0, 0} : e.tris[0].b, e.tris.empty() ? V3{0, 0, 0} : e.tris[0].a);
    t1 = cross(e.n, cross(t1, e.n));
    const double l1 = len(t1); t1 = l1 > 0 ? scale(t1, 1 / l1) : V3{1, 0, 0};
    const double dl = 1e-3 * std::max(1e-9, e.size);
    // keep the half-space constraint z <= 0 for the shifted points
    const V3 nOff = scale(e.n, dl);
    for (int trial = 0; trial < 3; ++trial) {
        const V3 inPlane = trial == 0 ? V3{0, 0, 0} : scale(t1, (trial == 1 ? 1 : -1) * 7.3e-4 * std::max(1e-9, e.size));
        double acc[6] = {0, 0, 0, 0, 0, 0};
        bool ok = true;
        for (int side = -1; side <= 1 && ok; side += 2) {
            V3 o = {obs[0] + side * nOff[0] + inPlane[0], obs[1] + side * nOff[1] + inPlane[1], obs[2] + side * nOff[2] + inPlane[2]};
            if (hs && o[2] > 0) o[2] = 0;
            for (const Tri& t : e.tris) {
                V3 n, s, d; triFrame(t, n, s, d);
                const V3 sl = {dot(b, s), dot(b, d), dot(b, n)};
                double St[6];
                if (!triStress(o, t, sl, G, nu, hs, St)) { ok = false; break; }
                for (int i = 0; i < 6; ++i) acc[i] += 0.5 * St[i];
            }
        }
        if (ok) { for (int i = 0; i < 6; ++i) S[i] = acc[i]; return true; }
    }
    return false;
}

ShearMatrix buildShear(const std::vector<Element>& el, const std::vector<V3>& dirIn, double G, double nu, bool hs, double cutoff, int threads) {
    const auto t0 = std::chrono::steady_clock::now();
    ShearMatrix K;
    const int n = int(el.size());
    K.n = n; K.dir = dirIn; K.dir.resize(size_t(n), V3{1, 0, 0}); K.selfKt.assign(size_t(n), 0.0);
    auto ktkn = [&](int i, int j, const V3& dj, double& kt, double& kn) {
        double S[6];
        kt = kn = 0;
        if (!elementStress(el[size_t(i)].c, el[size_t(j)], dj, G, nu, hs, S)) return false;
        const V3 tr = traction(S, el[size_t(i)].n);           // tension +
        kn = -dot(tr, el[size_t(i)].n);                        // compression +
        kt = -dot(tr, K.dir[size_t(i)]);                       // comp-positive shear traction along d_i
        return true;
    };
    for (int j = 0; j < n; ++j) {   // project the rake into the plane and normalise
        V3& d = K.dir[size_t(j)];
        d = sub(d, scale(el[size_t(j)].n, dot(d, el[size_t(j)].n)));
        const double l = len(d);
        if (l > 0) d = scale(d, 1 / l);
        else { d = cross(el[size_t(j)].n, V3{0, 0, 1}); if (len(d) == 0) d = V3{1, 0, 0}; d = scale(d, 1 / len(d)); }
    }
    // global sense of the TDE discontinuity: the self term K_tau(j,j) is invariant to the sign of d_j and of the element
    // normal (both enter twice), so one check fixes the convention for all elements: slip along d relieves the shear along d
    double convention = 1.0;
    if (n > 0) {
        double kt, kn; ktkn(0, 0, K.dir[0], kt, kn);
        if (kt > 0) convention = -1.0;
    }
    K.sense = convention;
    std::vector<std::vector<int>> cols(static_cast<size_t>(n));
    std::vector<std::vector<double>> vt(static_cast<size_t>(n)), vn(static_cast<size_t>(n));
    long pairs = 0, dropped = 0;
    const int nt = threads > 0 ? threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 4) num_threads(nt) reduction(+ : pairs, dropped)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double r = len(sub(el[size_t(i)].c, el[size_t(j)].c));
            if (cutoff > 0 && i != j && r > cutoff * 0.5 * (el[size_t(i)].size + el[size_t(j)].size)) { ++dropped; continue; }
            double kt, kn;
            if (!ktkn(i, j, scale(K.dir[size_t(j)], convention), kt, kn)) { ++dropped; continue; }
            cols[size_t(i)].push_back(j); vt[size_t(i)].push_back(kt); vn[size_t(i)].push_back(kn);
            ++pairs;
        }
    }
    K.start.assign(size_t(n) + 1, 0);
    for (int i = 0; i < n; ++i) K.start[size_t(i) + 1] = K.start[size_t(i)] + int(cols[size_t(i)].size());
    K.col.reserve(size_t(K.start.back())); K.kt.reserve(size_t(K.start.back())); K.kn.reserve(size_t(K.start.back()));
    for (int i = 0; i < n; ++i) {
        K.col.insert(K.col.end(), cols[size_t(i)].begin(), cols[size_t(i)].end());
        K.kt.insert(K.kt.end(), vt[size_t(i)].begin(), vt[size_t(i)].end());
        K.kn.insert(K.kn.end(), vn[size_t(i)].begin(), vn[size_t(i)].end());
        for (int q = K.start[size_t(i)]; q < K.start[size_t(i) + 1]; ++q) if (K.col[size_t(q)] == i) K.selfKt[size_t(i)] = K.kt[size_t(q)];
    }
    K.pairs = pairs; K.dropped = dropped;
    K.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return K;
}

std::vector<double> buildOpening(const std::vector<Element>& el, double G, double nu, bool hs) {
    const int n = int(el.size());
    std::vector<double> K(size_t(n) * size_t(n), 0.0);
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double S[6];
            if (!elementStress(el[size_t(i)].c, el[size_t(j)], el[size_t(j)].n, G, nu, hs, S)) continue;
            K[size_t(i) * size_t(n) + size_t(j)] = dot(traction(S, el[size_t(i)].n), el[size_t(i)].n);
        }
    // opening sense: the discontinuity along +n must compress its own element (K_ii < 0); otherwise flip globally
    if (n > 0 && K[0] > 0) for (double& v : K) v = -v;
    return K;
}

SlipResult solveSlip(const ShearMatrix& K, const std::vector<double>& excess, const std::vector<double>& mu, std::vector<double>& ds,
                     double tol, int maxIter) {
    SlipResult R;
    const int n = K.n;
    ds.resize(size_t(n), 0.0);
    std::vector<double> diag(size_t(n), -1.0);
    for (int i = 0; i < n; ++i)
        for (int q = K.start[size_t(i)]; q < K.start[size_t(i) + 1]; ++q)
            if (K.col[size_t(q)] == i) diag[size_t(i)] = K.kt[size_t(q)] - mu[size_t(i)] * K.kn[size_t(q)];
    auto f = [&](int i) {
        double v = excess[size_t(i)];
        for (int q = K.start[size_t(i)]; q < K.start[size_t(i) + 1]; ++q)
            v += (K.kt[size_t(q)] - mu[size_t(i)] * K.kn[size_t(q)]) * ds[size_t(K.col[size_t(q)])];
        return v;
    };
    // only elements that can slip matter: start from the over-stressed ones, grow the set when neighbours fail
    for (R.iterations = 0; R.iterations < maxIter; ++R.iterations) {
        double res = 0;
        for (int i = 0; i < n; ++i) {
            const double fi = f(i);
            if (ds[size_t(i)] <= 0 && fi <= 0) continue;
            const double d = diag[size_t(i)] < 0 ? diag[size_t(i)] : -1.0;
            const double nv = std::max(0.0, ds[size_t(i)] - fi / d);
            res = std::max(res, ds[size_t(i)] > 0 || nv > 0 ? std::fabs(fi) : 0.0);
            ds[size_t(i)] = nv;
        }
        R.residual = res;
        if (res < tol) { R.converged = true; break; }
    }
    R.slipping = 0;
    for (int i = 0; i < n; ++i) if (ds[size_t(i)] > 0) ++R.slipping;
    return R;
}

namespace {
// dense LU solve (partial pivoting) of the sub-system on idx
bool denseSolve(std::vector<double> A, int m, std::vector<double>& b) {
    for (int k = 0; k < m; ++k) {
        int p = k; double best = std::fabs(A[size_t(k) * size_t(m) + size_t(k)]);
        for (int r = k + 1; r < m; ++r) { const double v = std::fabs(A[size_t(r) * size_t(m) + size_t(k)]); if (v > best) { best = v; p = r; } }
        if (!(best > 0)) return false;
        if (p != k) { for (int c = 0; c < m; ++c) std::swap(A[size_t(k) * size_t(m) + size_t(c)], A[size_t(p) * size_t(m) + size_t(c)]); std::swap(b[size_t(k)], b[size_t(p)]); }
        const double piv = A[size_t(k) * size_t(m) + size_t(k)];
        for (int r = k + 1; r < m; ++r) {
            const double l = A[size_t(r) * size_t(m) + size_t(k)] / piv;
            if (l == 0) continue;
            for (int c = k; c < m; ++c) A[size_t(r) * size_t(m) + size_t(c)] -= l * A[size_t(k) * size_t(m) + size_t(c)];
            b[size_t(r)] -= l * b[size_t(k)];
        }
    }
    for (int k = m - 1; k >= 0; --k) {
        double s = b[size_t(k)];
        for (int c = k + 1; c < m; ++c) s -= A[size_t(k) * size_t(m) + size_t(c)] * b[size_t(c)];
        b[size_t(k)] = s / A[size_t(k) * size_t(m) + size_t(k)];
    }
    return true;
}
}  // namespace

SlipResult solveOpening(const std::vector<double>& Ko, int n, const std::vector<double>& pNet, std::vector<double>& w, double tol, int maxIter) {
    // active set: open = elements with w > 0; solve K_oo w_o = -p_o; close negative openings; reopen closed elements whose
    // traction would be tensile (p + K w > 0); a handful of iterations for the convex problem
    SlipResult R;
    w.assign(size_t(n), 0.0);
    std::vector<char> open(size_t(n), 0);
    for (int i = 0; i < n; ++i) open[size_t(i)] = pNet[size_t(i)] > 0;
    for (R.iterations = 0; R.iterations < std::min(maxIter, 200); ++R.iterations) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i) if (open[size_t(i)]) idx.push_back(i);
        const int m = int(idx.size());
        std::fill(w.begin(), w.end(), 0.0);
        if (m > 0) {
            std::vector<double> A(static_cast<size_t>(m) * static_cast<size_t>(m)), b(static_cast<size_t>(m));
            for (int r = 0; r < m; ++r) {
                b[size_t(r)] = -pNet[size_t(idx[size_t(r)])];
                for (int c = 0; c < m; ++c) A[size_t(r) * size_t(m) + size_t(c)] = Ko[size_t(idx[size_t(r)]) * size_t(n) + size_t(idx[size_t(c)])];
            }
            if (!denseSolve(A, m, b)) return R;
            for (int r = 0; r < m; ++r) w[size_t(idx[size_t(r)])] = b[size_t(r)];
        }
        bool change = false;
        double res = 0;
        for (int i = 0; i < n; ++i) {
            double g = pNet[size_t(i)];
            for (int j = 0; j < n; ++j) g += Ko[size_t(i) * size_t(n) + size_t(j)] * w[size_t(j)];
            if (open[size_t(i)] && w[size_t(i)] < 0) { open[size_t(i)] = 0; change = true; }
            else if (!open[size_t(i)] && g > tol) { open[size_t(i)] = 1; change = true; }
            if (open[size_t(i)]) res = std::max(res, std::fabs(g));
        }
        R.residual = res;
        if (!change) { R.converged = true; break; }
    }
    for (double& v : w) v = std::max(0.0, v);
    R.slipping = 0;
    for (double v : w) if (v > 0) ++R.slipping;
    return R;
}

double pennyVolume(double p, double a, double E, double nu) { return 16.0 * (1 - nu * nu) * p * a * a * a / (3.0 * E); }
double pennyOpening(double p, double a, double r, double E, double nu) {
    return r >= a ? 0.0 : 8.0 * (1 - nu * nu) * p * std::sqrt(a * a - r * r) / (M_PI * E);
}

}  // namespace ddm
}  // namespace ressim
