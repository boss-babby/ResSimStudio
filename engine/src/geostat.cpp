// C28 geostatistics: variograms, kriging, sequential Gaussian simulation (see geostat.h)
#include "ressim/geostat.h"
#include "ressim/grid.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <unordered_map>

namespace ressim {
namespace geostat {
namespace {
constexpr double kPi = 3.141592653589793;

// anisotropic, range-normalised separation: |h| = 1 at the range in every direction
double scaledDist(const Variogram& v, double dx, double dy, double dz) {
    const double a = v.azimuthDeg * kPi / 180.0;
    const double ux = std::sin(a), uy = std::cos(a);                 // major axis (azimuth clockwise from north)
    const double hMaj = dx * ux + dy * uy, hMin = -dx * uy + dy * ux;
    const double rM = std::max(v.rangeMax, 1e-12), rm = std::max(v.rangeMin, 1e-12), rv = std::max(v.rangeVert, 1e-12);
    return std::sqrt((hMaj / rM) * (hMaj / rM) + (hMin / rm) * (hMin / rm) + (dz / rv) * (dz / rv));
}
double unitGamma(VarioModel m, double h) {
    switch (m) {
        case VarioModel::Spherical: return h >= 1 ? 1.0 : 1.5 * h - 0.5 * h * h * h;
        case VarioModel::Exponential: return 1.0 - std::exp(-3.0 * h);   // practical range
        case VarioModel::Gaussian: return 1.0 - std::exp(-3.0 * h * h);
    }
    return 1.0;
}
// transformed coordinates (range-normalised along the anisotropy axes) for the neighbour search
Xyz toScaled(const Variogram& v, const Xyz& p) {
    const double a = v.azimuthDeg * kPi / 180.0;
    const double ux = std::sin(a), uy = std::cos(a);
    return {(p[0] * ux + p[1] * uy) / std::max(v.rangeMax, 1e-12), (-p[0] * uy + p[1] * ux) / std::max(v.rangeMin, 1e-12), p[2] / std::max(v.rangeVert, 1e-12)};
}

// uniform bucket grid in scaled space (bucket = 0.5 range) holding conditioning points
struct Buckets {
    double h = 0.5;
    std::unordered_map<long long, std::vector<int>> map;
    static long long key(long long i, long long j, long long k) { return ((i + 1048576LL) << 42) ^ ((j + 1048576LL) << 21) ^ (k + 1048576LL); }
    void add(const Xyz& s, int id) { map[key(llround(std::floor(s[0] / h)), llround(std::floor(s[1] / h)), llround(std::floor(s[2] / h)))].push_back(id); }
    template <class F> void around(const Xyz& s, double r, F&& f) const {
        const long long i0 = llround(std::floor((s[0] - r) / h)), i1 = llround(std::floor((s[0] + r) / h));
        const long long j0 = llround(std::floor((s[1] - r) / h)), j1 = llround(std::floor((s[1] + r) / h));
        const long long k0 = llround(std::floor((s[2] - r) / h)), k1 = llround(std::floor((s[2] + r) / h));
        if ((i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1) > 200000) { for (const auto& e : map) for (int id : e.second) f(id); return; }
        for (long long i = i0; i <= i1; ++i) for (long long j = j0; j <= j1; ++j) for (long long k = k0; k <= k1; ++k) {
            auto it = map.find(key(i, j, k));
            if (it != map.end()) for (int id : it->second) f(id);
        }
    }
};

// dense solve (partial pivoting); false if singular
bool solve(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int c = 0; c < n; ++c) {
        int p = c; for (int r = c + 1; r < n; ++r) if (std::fabs(A[size_t(r) * n + c]) > std::fabs(A[size_t(p) * n + c])) p = r;
        if (std::fabs(A[size_t(p) * n + c]) < 1e-300) return false;
        if (p != c) { for (int k = 0; k < n; ++k) std::swap(A[size_t(c) * n + k], A[size_t(p) * n + k]); std::swap(b[size_t(c)], b[size_t(p)]); }
        for (int r = c + 1; r < n; ++r) {
            const double f = A[size_t(r) * n + c] / A[size_t(c) * n + c];
            if (f == 0) continue;
            for (int k = c; k < n; ++k) A[size_t(r) * n + k] -= f * A[size_t(c) * n + k];
            b[size_t(r)] -= f * b[size_t(c)];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = b[size_t(r)];
        for (int k = r + 1; k < n; ++k) s -= A[size_t(r) * n + k] * b[size_t(k)];
        b[size_t(r)] = s / A[size_t(r) * n + r];
    }
    return true;
}

// kriging at x from the conditioning points P (indices nb): ordinary (Lagrange) or simple (known mean)
bool krigeAt(const Variogram& vg, const Xyz& x, const std::vector<Point>& P, const std::vector<int>& nb, bool simple, double mean,
             double& est, double& var) {
    const int k = int(nb.size());
    const double c0 = vg.nugget + vg.sill;
    if (k == 0) { est = mean; var = c0; return false; }
    const int n = simple ? k : k + 1;
    std::vector<double> A(size_t(n) * n, 0.0), b(size_t(n), 0.0);
    for (int i = 0; i < k; ++i) {
        const Point& a = P[size_t(nb[size_t(i)])];
        for (int j = 0; j <= i; ++j) {
            const Point& c = P[size_t(nb[size_t(j)])];
            const double cv = i == j ? c0 : covariance(vg, a.x - c.x, a.y - c.y, a.z - c.z);
            A[size_t(i) * n + j] = A[size_t(j) * n + i] = cv;
        }
        // a data point exactly at x: the nugget is a measurement-scale effect - kriging is exact there (GSLIB convention)
        const double dd = std::fabs(a.x - x[0]) + std::fabs(a.y - x[1]) + std::fabs(a.z - x[2]);
        b[size_t(i)] = dd < 1e-9 ? c0 : covariance(vg, a.x - x[0], a.y - x[1], a.z - x[2]);
    }
    if (!simple) { for (int i = 0; i < k; ++i) { A[size_t(i) * n + k] = 1; A[size_t(k) * n + i] = 1; } b[size_t(k)] = 1; }
    for (int i = 0; i < k; ++i) A[size_t(i) * n + i] *= 1.0 + 1e-10;   // tiny regularisation (duplicate points)
    std::vector<double> rhs = b, w = b;
    if (!solve(A, w, n)) { est = mean; var = c0; return false; }
    est = simple ? mean : 0.0; var = c0;
    for (int i = 0; i < k; ++i) {
        const double v = P[size_t(nb[size_t(i)])].v;
        est += w[size_t(i)] * (simple ? v - mean : v);
        var -= w[size_t(i)] * rhs[size_t(i)];
    }
    if (!simple) var -= w[size_t(k)];
    var = std::max(0.0, var);
    return true;
}

std::vector<int> nearest(const Variogram& vg, const Buckets& B, const std::vector<Point>& P, const std::vector<Xyz>& S, const Xyz& x, int maxN, double radius) {
    const Xyz s = toScaled(vg, x);
    std::vector<std::pair<double, int>> c;
    B.around(s, radius, [&](int id) {
        const Xyz& q = S[size_t(id)];
        const double d2 = (q[0] - s[0]) * (q[0] - s[0]) + (q[1] - s[1]) * (q[1] - s[1]) + (q[2] - s[2]) * (q[2] - s[2]);
        if (d2 <= radius * radius) c.push_back({d2, id});
    });
    if (int(c.size()) > maxN) { std::nth_element(c.begin(), c.begin() + maxN, c.end()); c.resize(size_t(maxN)); }
    std::vector<int> out; for (const auto& e : c) out.push_back(e.second);
    (void)P;
    return out;
}

// standard normal quantile (Acklam's rational approximation, |error| < 1.2e-9)
double normInv(double p) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
    p = std::clamp(p, 1e-12, 1 - 1e-12);
    if (p < 0.02425) { const double q = std::sqrt(-2 * std::log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    if (p > 1 - 0.02425) { const double q = std::sqrt(-2 * std::log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    const double q = p - 0.5, r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
}  // namespace

double gamma(const Variogram& v, double dx, double dy, double dz) {
    const double h = scaledDist(v, dx, dy, dz);
    if (h <= 0) return 0.0;
    return v.nugget + v.sill * unitGamma(v.model, h);
}

KrigeOut krige(const std::vector<Xyz>& targets, const std::vector<Point>& data, const Options& opt) {
    KrigeOut o; o.est.assign(targets.size(), 0.0); o.var.assign(targets.size(), opt.vg.nugget + opt.vg.sill);
    double mean = opt.skMean;
    if (!opt.simpleKriging && !data.empty()) { mean = 0; for (const auto& p : data) mean += p.v; mean /= double(data.size()); }
    Buckets B; std::vector<Xyz> S(data.size());
    for (size_t i = 0; i < data.size(); ++i) { S[i] = toScaled(opt.vg, {data[i].x, data[i].y, data[i].z}); B.add(S[i], int(i)); }
    const double radius = std::max(0.5, opt.searchRanges);
    #pragma omp parallel for schedule(dynamic, 64)
    for (long t = 0; t < long(targets.size()); ++t) {
        const auto nb = nearest(opt.vg, B, data, S, targets[size_t(t)], std::max(1, opt.maxNbr), radius);
        double e = mean, v = 0;
        if (!krigeAt(opt.vg, targets[size_t(t)], data, nb, opt.simpleKriging, mean, e, v)) {
            #pragma omp atomic
            o.noData++;
        }
        o.est[size_t(t)] = e; o.var[size_t(t)] = v;
    }
    return o;
}

std::vector<double> sgs(const std::vector<Xyz>& targets, const std::vector<Point>& data, const Options& opt) {
    const size_t nt = targets.size();
    std::vector<double> out(nt, 0.0);
    // normal-score table (GSLIB nscore: plotting positions (i + 0.5) / n, ties get the mean score)
    std::vector<double> zs, ys;
    std::vector<Point> P;
    if (!data.empty()) {
        std::vector<size_t> o(data.size()); std::iota(o.begin(), o.end(), size_t(0));
        std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return data[a].v < data[b].v; });
        const double n = double(data.size());
        std::vector<double> score(data.size());
        for (size_t i = 0; i < o.size();) {
            size_t j = i; while (j + 1 < o.size() && data[o[j + 1]].v == data[o[i]].v) ++j;
            double s = 0; for (size_t q = i; q <= j; ++q) s += normInv((double(q) + 0.5) / n);
            s /= double(j - i + 1);
            for (size_t q = i; q <= j; ++q) score[o[q]] = s;
            i = j + 1;
        }
        for (size_t i = 0; i < data.size(); ++i) {
            Point p = data[i]; if (opt.normalScore) p.v = score[i];
            P.push_back(p);
        }
        for (size_t q : o) { zs.push_back(data[q].v); ys.push_back(score[q]); }
    }
    auto back = [&](double y) {
        if (!opt.normalScore || zs.empty()) return y;
        if (y <= ys.front()) return zs.front();     // tails: clamped to the data range (GSLIB ltail / utail = 1 with zmin / zmax = data limits)
        if (y >= ys.back()) return zs.back();
        const size_t k = size_t(std::upper_bound(ys.begin(), ys.end(), y) - ys.begin());
        const double t = (y - ys[k - 1]) / std::max(1e-300, ys[k] - ys[k - 1]);
        return zs[k - 1] + t * (zs[k] - zs[k - 1]);
    };
    // Gaussian-space variogram: total sill 1
    Variogram vg = opt.vg;
    if (opt.normalScore) { const double tot = std::max(1e-12, vg.nugget + vg.sill); vg.nugget /= tot; vg.sill /= tot; }
    const double mean = opt.normalScore ? 0.0 : (data.empty() ? opt.skMean : [&] { double m = 0; for (const auto& p : data) m += p.v; return m / double(data.size()); }());
    // conditioning set: data (snapped to the nearest target node: exact reproduction there) + simulated nodes
    Buckets B; std::vector<Xyz> S;
    std::vector<char> done(nt, 0);
    {
        Buckets TB; std::vector<Xyz> TS(nt);
        for (size_t t = 0; t < nt; ++t) { TS[t] = toScaled(vg, targets[t]); TB.add(TS[t], int(t)); }
        std::vector<Point> P2;
        for (const auto& p : P) {
            const Xyz s = toScaled(vg, {p.x, p.y, p.z});
            int best = -1; double bd = 1e300;
            for (double r : {0.25, 1.0, 4.0}) {
                TB.around(s, r, [&](int id) {
                    const Xyz& q = TS[size_t(id)];
                    const double dx = targets[size_t(id)][0] - p.x, dy = targets[size_t(id)][1] - p.y, dz = targets[size_t(id)][2] - p.z;
                    const double d2 = dx * dx + dy * dy + dz * dz; (void)q;
                    if (d2 < bd) { bd = d2; best = id; }
                });
                if (best >= 0) break;
            }
            Point q = p;
            if (best >= 0 && !done[size_t(best)]) { done[size_t(best)] = 1; out[size_t(best)] = p.v; q.x = targets[size_t(best)][0]; q.y = targets[size_t(best)][1]; q.z = targets[size_t(best)][2]; }
            P2.push_back(q);
        }
        P = std::move(P2);
    }
    for (size_t i = 0; i < P.size(); ++i) { S.push_back(toScaled(vg, {P[i].x, P[i].y, P[i].z})); B.add(S.back(), int(i)); }
    std::mt19937_64 rng(opt.seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    std::vector<size_t> path(nt); std::iota(path.begin(), path.end(), size_t(0));
    std::shuffle(path.begin(), path.end(), rng);
    const double radius = std::max(0.5, opt.searchRanges);
    for (size_t t : path) {
        if (done[t]) continue;
        const auto nb = nearest(vg, B, P, S, targets[t], std::max(1, opt.maxNbr), radius);
        double e = mean, v = vg.nugget + vg.sill;
        krigeAt(vg, targets[t], P, nb, true, mean, e, v);
        const double y = e + std::sqrt(std::max(0.0, v)) * N01(rng);
        out[t] = y; done[t] = 1;
        Point q; q.x = targets[t][0]; q.y = targets[t][1]; q.z = targets[t][2]; q.v = y;
        P.push_back(q); S.push_back(toScaled(vg, targets[t])); B.add(S.back(), int(P.size() - 1));
    }
    for (double& y : out) y = back(y);
    return out;
}

ExpVario experimental(const std::vector<Point>& pts, double lagWidth, int nLags, int direction, double azimuthDeg, long maxPairs) {
    ExpVario E;
    if (lagWidth <= 0 || nLags <= 0) return E;
    std::vector<double> s(size_t(nLags), 0.0); std::vector<long> c(size_t(nLags), 0);
    const size_t n = pts.size();
    const double a = azimuthDeg * kPi / 180.0, ux = std::sin(a), uy = std::cos(a), cosTol = std::cos(22.5 * kPi / 180.0);
    std::mt19937_64 rng(7);
    const long total = long(n) * long(n > 0 ? n - 1 : 0) / 2;
    const bool sample = total > maxPairs;
    std::uniform_int_distribution<size_t> U(0, n > 0 ? n - 1 : 0);
    auto addPair = [&](const Point& p, const Point& q) {
        const double dx = q.x - p.x, dy = q.y - p.y, dz = q.z - p.z;
        double h;
        if (direction == 1) { if (std::hypot(dx, dy) > 0.5 * lagWidth) return; h = std::fabs(dz); }
        else {
            if (std::fabs(dz) > 0.5 * lagWidth && std::fabs(dz) > 1e-9) return;
            h = std::hypot(dx, dy);
            if (direction == 2 && h > 0 && std::fabs(dx * ux + dy * uy) / h < cosTol) return;
        }
        const int k = int(std::floor(h / lagWidth + 0.5)) - 1;   // lag centres lagWidth, 2 lagWidth, ...
        if (k < 0 || k >= nLags) return;
        s[size_t(k)] += 0.5 * (p.v - q.v) * (p.v - q.v); ++c[size_t(k)];
    };
    if (!sample) { for (size_t i = 0; i < n; ++i) for (size_t j = i + 1; j < n; ++j) addPair(pts[i], pts[j]); }
    else for (long r = 0; r < maxPairs; ++r) { const size_t i = U(rng), j = U(rng); if (i != j) addPair(pts[i], pts[j]); }
    for (int k = 0; k < nLags; ++k) if (c[size_t(k)] > 0) { E.lag.push_back(lagWidth * (k + 1)); E.gamma.push_back(s[size_t(k)] / double(c[size_t(k)])); E.pairs.push_back(c[size_t(k)]); }
    return E;
}

ApplyReport applyToCase(Case& c, const ApplyInput& in) {
    ApplyReport R;
    GridInput& gi = c.grid;
    const size_t n = size_t(gi.nx) * size_t(gi.ny) * size_t(gi.nz);
    Grid g;
    try { g.build(c); } catch (const std::exception& e) { R.message = std::string("grid: ") + e.what(); return R; }
    if (g.center.size() < n) { R.message = "grid geometry not available"; return R; }
    std::vector<double>* prop = in.property == 0 ? &gi.poro : in.property == 1 ? &gi.permx : &gi.ntg;
    if (prop->size() < n) prop->resize(n, in.property == 0 ? 0.2 : in.property == 1 ? 100.0 : 1.0);
    const bool logk = in.property == 1 && in.logPerm;
    auto fwd = [&](double v) { return logk ? std::log10(std::max(v, 1e-6)) : v; };
    auto inv = [&](double v) { return logk ? std::pow(10.0, v) : v; };
    std::vector<Point> data;
    if (in.fromWells)
        for (const auto& w : c.wells)
            for (const auto& cp : w.completions)
                for (int k = cp.k1; k <= cp.k2; ++k) {
                    if (cp.i < 0 || cp.j < 0 || k < 0 || cp.i >= gi.nx || cp.j >= gi.ny || k >= gi.nz) continue;
                    const size_t q = size_t(cp.i) + size_t(gi.nx) * (size_t(cp.j) + size_t(gi.ny) * size_t(k));
                    if (g.activeOf[q] < 0) continue;
                    Point p; p.x = g.center[q].x; p.y = g.center[q].y; p.z = g.center[q].z; p.v = fwd((*prop)[q]);
                    data.push_back(p);
                }
    for (Point p : in.data) { p.v = fwd(p.v); data.push_back(p); }
    R.dataPoints = int(data.size());
    if (data.empty() && in.method == 0) { R.message = "kriging needs hard data (no completed cells / points)"; return R; }
    std::vector<Xyz> T; std::vector<size_t> gl;
    for (int a = 0; a < g.numActive(); ++a) {
        const size_t q = size_t(g.globalOf[size_t(a)]);
        if (q >= n) continue;   // local (LGR) cells keep their values
        T.push_back({g.center[q].x, g.center[q].y, g.center[q].z}); gl.push_back(q);
    }
    Options o = in.opt;
    if (logk && in.method == 1 && !o.normalScore && data.empty()) o.skMean = fwd(o.skMean);
    std::vector<double> v;
    if (in.method == 0) v = krige(T, data, o).est; else v = sgs(T, data, o);
    double s = 0, s2 = 0, mn = 1e300, mx = -1e300;
    for (size_t i = 0; i < T.size(); ++i) {
        double x = std::clamp(inv(v[i]), in.minV, in.maxV);
        if (in.property == 0) x = std::clamp(x, 0.0, 1.0);
        if (in.property == 2) x = std::clamp(x, 0.0, 1.0);
        if (in.property == 1) x = std::max(x, 0.0);
        (*prop)[gl[i]] = x;
        s += x; s2 += x * x; mn = std::min(mn, x); mx = std::max(mx, x);
    }
    if (in.property == 1) {
        if (gi.permy.size() < n) gi.permy.resize(n, 100.0);
        if (gi.permz.size() < n) gi.permz.resize(n, 10.0);
        for (size_t q : gl) { gi.permy[q] = gi.permx[q]; gi.permz[q] = in.kvkh * gi.permx[q]; }
    }
    R.cells = int(T.size());
    R.mean = T.empty() ? 0 : s / double(T.size());
    R.stdDev = T.empty() ? 0 : std::sqrt(std::max(0.0, s2 / double(T.size()) - R.mean * R.mean));
    R.minV = mn; R.maxV = mx; R.ok = true;
    char b[256];
    std::snprintf(b, sizeof b, "%s of %s on %d cells from %d hard data: mean %.4g, std %.4g, range %.4g .. %.4g", in.method == 0 ? "Kriging" : "SGS",
                  in.property == 0 ? "PORO" : in.property == 1 ? "PERMX" : "NTG", R.cells, R.dataPoints, R.mean, R.stdDev, mn, mx);
    R.message = b;
    return R;
}

// C41: facies by truncated Gaussian simulation (Matheron et al. 1987; Deutsch 2002 ch. 6) on the case grid
ApplyReport applyFaciesModel(Case& c, const PropModelInput& m) {
    ApplyReport R;
    GridInput& gi = c.grid;
    const int K = int(m.facProp.size());
    if (K < 2) { R.message = "facies model needs at least two facies proportions"; return R; }
    std::vector<double> p(m.facProp);
    double sp = 0;
    for (double& x : p) { x = std::max(0.0, x); sp += x; }
    if (!(sp > 0)) { R.message = "facies proportions are all zero"; return R; }
    std::vector<double> cum(size_t(K) + 1, 0.0), thr(size_t(K) - 1);
    for (int k = 0; k < K; ++k) { p[size_t(k)] /= sp; cum[size_t(k) + 1] = cum[size_t(k)] + p[size_t(k)]; }
    for (int k = 0; k + 1 < K; ++k) thr[size_t(k)] = normInv(cum[size_t(k) + 1]);
    auto gaussOf = [&](int f) { f = std::clamp(f, 0, K - 1); return normInv(std::clamp(0.5 * (cum[size_t(f)] + cum[size_t(f) + 1]), 1e-9, 1 - 1e-9)); };
    const size_t n = size_t(gi.nx) * size_t(gi.ny) * size_t(gi.nz);
    Grid g;
    try { g.build(c); } catch (const std::exception& e) { R.message = std::string("grid: ") + e.what(); return R; }
    if (g.center.size() < n) { R.message = "grid geometry not available"; return R; }
    std::vector<Point> data;
    if (m.fromWells && m.facSatnum && gi.satnum.size() >= n)
        for (const auto& w : c.wells)
            for (const auto& cp : w.completions)
                for (int k = cp.k1; k <= cp.k2; ++k) {
                    if (cp.i < 0 || cp.j < 0 || k < 0 || cp.i >= gi.nx || cp.j >= gi.ny || k >= gi.nz) continue;
                    const size_t q = size_t(cp.i) + size_t(gi.nx) * (size_t(cp.j) + size_t(gi.ny) * size_t(k));
                    if (g.activeOf[q] < 0) continue;
                    Point pt; pt.x = g.center[q].x; pt.y = g.center[q].y; pt.z = g.center[q].z; pt.v = gaussOf(gi.satnum[q] - 1);
                    data.push_back(pt);
                }
    for (const auto& d : m.data) { Point pt; pt.x = d[0]; pt.y = d[1]; pt.z = d[2]; pt.v = gaussOf(int(std::lround(d[3])) - 1); data.push_back(pt); }
    R.dataPoints = int(data.size());
    std::vector<Xyz> T; std::vector<size_t> gl;
    for (int a = 0; a < g.numActive(); ++a) {
        const size_t q = size_t(g.globalOf[size_t(a)]);
        if (q >= n) continue;
        T.push_back({g.center[q].x, g.center[q].y, g.center[q].z}); gl.push_back(q);
    }
    Options o;
    const double tot = std::max(1e-12, std::max(0.0, m.nugget) + std::max(0.0, m.sill));
    o.vg.model = VarioModel(std::clamp(m.model, 0, 2));
    o.vg.nugget = std::max(0.0, m.nugget) / tot; o.vg.sill = std::max(0.0, m.sill) / tot;   // unit-sill Gaussian field
    o.vg.rangeMax = std::max(1e-9, m.rangeMax); o.vg.rangeMin = std::max(1e-9, m.rangeMin); o.vg.rangeVert = std::max(1e-9, m.rangeVert);
    o.vg.azimuthDeg = m.azimuthDeg; o.seed = m.seed; o.maxNbr = std::max(1, m.maxNbr);
    o.normalScore = false; o.skMean = 0;
    const std::vector<double> y = sgs(T, data, o);
    if (gi.poro.size() < n) gi.poro.resize(n, 0.2);
    if (gi.permx.size() < n) gi.permx.resize(n, 100.0);
    if (gi.permy.size() < n) gi.permy.resize(n, 100.0);
    if (gi.permz.size() < n) gi.permz.resize(n, 10.0);
    if (m.facSatnum && gi.satnum.size() < n) gi.satnum.resize(n, 1);
    std::vector<int> count(size_t(K), 0);
    for (size_t i = 0; i < T.size(); ++i) {
        int f = 0;
        while (f + 1 < K && y[i] >= thr[size_t(f)]) ++f;
        ++count[size_t(f)];
        const size_t q = gl[i];
        if (size_t(f) < m.facPoro.size() && m.facPoro[size_t(f)] >= 0) gi.poro[q] = std::clamp(m.facPoro[size_t(f)], 0.0, 1.0);
        if (size_t(f) < m.facPerm.size() && m.facPerm[size_t(f)] >= 0) { gi.permx[q] = gi.permy[q] = m.facPerm[size_t(f)]; gi.permz[q] = m.kvkh * m.facPerm[size_t(f)]; }
        if (m.facSatnum) gi.satnum[q] = f + 1;
    }
    R.cells = int(T.size());
    R.ok = true;
    std::string s = "Truncated Gaussian facies on " + std::to_string(R.cells) + " cells from " + std::to_string(R.dataPoints) + " hard data; proportions";
    char b[64];
    for (int k = 0; k < K; ++k) { std::snprintf(b, sizeof b, " %d: %.3f (target %.3f)", k + 1, R.cells ? double(count[size_t(k)]) / R.cells : 0.0, p[size_t(k)]); s += b; }
    R.mean = R.cells ? double(count[0]) / R.cells : 0;
    R.message = s;
    return R;
}

ApplyReport applyPropModel(Case& c, const PropModelInput& m) {
    if (m.property == 3) return applyFaciesModel(c, m);
    ApplyInput in;
    in.property = std::clamp(m.property, 0, 2);
    in.method = m.method == 0 ? 0 : 1;
    in.opt.vg.model = VarioModel(std::clamp(m.model, 0, 2));
    in.opt.vg.nugget = std::max(0.0, m.nugget); in.opt.vg.sill = std::max(0.0, m.sill);
    in.opt.vg.rangeMax = std::max(1e-9, m.rangeMax); in.opt.vg.rangeMin = std::max(1e-9, m.rangeMin);
    in.opt.vg.rangeVert = std::max(1e-9, m.rangeVert); in.opt.vg.azimuthDeg = m.azimuthDeg;
    in.opt.seed = m.seed; in.opt.maxNbr = std::max(1, m.maxNbr);
    in.fromWells = m.fromWells;
    for (const auto& d : m.data) { Point p; p.x = d[0]; p.y = d[1]; p.z = d[2]; p.v = d[3]; in.data.push_back(p); }
    in.kvkh = m.kvkh; in.minV = m.minV; in.maxV = m.maxV;
    ApplyReport R = applyToCase(c, in);
    if (R.ok && in.property == 0 && m.permFromPoro) {   // permeability from the simulated porosity (cells the model touched)
        GridInput& g = c.grid;
        const size_t n = g.poro.size();
        g.permx.resize(n, 100.0); g.permy.resize(n, 100.0); g.permz.resize(n, 10.0);
        for (size_t q = 0; q < n; ++q) {
            const double k = std::clamp(std::pow(10.0, m.ppA + m.ppB * g.poro[q]), 1e-4, 1e5);
            g.permx[q] = g.permy[q] = k; g.permz[q] = m.kvkh * k;
        }
        R.message += "; PERM = 10^(a + b PORO)";
    }
    return R;
}

}  // namespace geostat
}  // namespace ressim
