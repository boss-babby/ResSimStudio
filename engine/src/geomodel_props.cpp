// C32 geomodelling: well-log blocking, facies (SIS / truncated Gaussian / object channels), petrophysics (SGS with
// optional collocated secondary), poro-perm transforms, volumetrics with contacts and saturation-height, percentiles
#include "ressim/geomodel.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
std::string upperS(std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; }
int defaultMethod(const std::string& m) {
    const std::string u = upperS(m);
    if (u.rfind("PERM", 0) == 0 || u == "K" || u == "KH" || u == "KV") return 1;
    if (u.rfind("FACIES", 0) == 0 || u.rfind("LITH", 0) == 0 || u == "ZONE") return 3;
    return 0;
}
// stratigraphic (sim-box) coordinates: plan centre and layer index x mean layer thickness - geostatistics is done in
// this space so faults and structural relief do not distort the spatial correlation
std::vector<geostat::Xyz> stratCoords(const GeoGrid& g) {
    const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
    double th = 0; int n = 0;
    for (size_t c = 0; c < g.bulk.size(); ++c) if (g.grid.actnum.empty() || g.grid.actnum[c]) { th += g.bot[c] - g.top[c]; ++n; }
    th = n ? th / n : 1.0;
    std::vector<geostat::Xyz> x(size_t(nx) * size_t(ny) * size_t(nz));
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const size_t c = size_t(i + nx * (j + ny * k));
        x[c] = {g.center[c][0], g.center[c][1], (k + 0.5) * th};
    }
    return x;
}
bool active(const GeoGrid& g, size_t c) { return g.grid.actnum.empty() || g.grid.actnum[c] != 0; }
double invNorm(double p) {   // Acklam's rational approximation of the standard normal quantile
    p = std::min(std::max(p, 1e-12), 1 - 1e-12);
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
    const double pl = 0.02425;
    if (p < pl) { const double q = std::sqrt(-2 * std::log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    if (p > 1 - pl) { const double q = std::sqrt(-2 * std::log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    const double q = p - 0.5, r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
// normal-score transform table of a sample (sorted values with their Gaussian scores); back-transform interpolates
struct NScore {
    std::vector<double> v, y;
    void build(std::vector<double> d) {
        std::sort(d.begin(), d.end());
        v = d; y.resize(d.size());
        for (size_t k = 0; k < d.size(); ++k) y[k] = invNorm((k + 0.5) / double(d.size()));
    }
    double fwd(double x) const {
        if (v.empty()) return 0;
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin()) return y.front();
        if (it == v.end()) return y.back();
        const size_t k = size_t(it - v.begin());
        const double f = v[k] > v[k - 1] ? (x - v[k - 1]) / (v[k] - v[k - 1]) : 0.5;
        return y[k - 1] + f * (y[k] - y[k - 1]);
    }
    double back(double g) const {
        if (v.empty()) return 0;
        if (v.size() == 1) return v[0];
        if (g <= y.front()) {   // lower tail: linear in the cdf down to the minimum minus half a class width
            const double lo = v.front() - 0.5 * (v[1] - v[0]);
            return lo + (v.front() - lo) * normCdf(g) / normCdf(y.front());
        }
        if (g >= y.back()) {
            const size_t n = v.size();
            const double hi = v.back() + 0.5 * (v[n - 1] - v[n - 2]);
            return v.back() + (hi - v.back()) * (normCdf(g) - normCdf(y.back())) / std::max(1e-12, 1 - normCdf(y.back()));
        }
        auto it = std::lower_bound(y.begin(), y.end(), g);
        const size_t k = size_t(it - y.begin());
        const double f = (g - y[k - 1]) / (y[k] - y[k - 1]);
        return v[k - 1] + f * (v[k] - v[k - 1]);
    }
};
// small dense solve (Gaussian elimination with partial pivoting); returns false if singular
bool solveDense(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int c = 0; c < n; ++c) {
        int p = c;
        for (int r = c + 1; r < n; ++r) if (std::fabs(A[size_t(r * n + c)]) > std::fabs(A[size_t(p * n + c)])) p = r;
        if (std::fabs(A[size_t(p * n + c)]) < 1e-14) return false;
        if (p != c) { for (int k = 0; k < n; ++k) std::swap(A[size_t(c * n + k)], A[size_t(p * n + k)]); std::swap(b[size_t(c)], b[size_t(p)]); }
        for (int r = c + 1; r < n; ++r) {
            const double f = A[size_t(r * n + c)] / A[size_t(c * n + c)];
            if (f == 0) continue;
            for (int k = c; k < n; ++k) A[size_t(r * n + k)] -= f * A[size_t(c * n + k)];
            b[size_t(r)] -= f * b[size_t(c)];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = b[size_t(r)];
        for (int k = r + 1; k < n; ++k) s -= A[size_t(r * n + k)] * b[size_t(k)];
        b[size_t(r)] = s / A[size_t(r * n + r)];
    }
    return true;
}
}  // namespace

// ================================================================ blocking
std::vector<Blocked> blockWells(const GeoGrid& g, const std::vector<GeoWell>& wells, const std::map<std::string, int>& methods) {
    std::vector<Blocked> out;
    if (g.grid.zcorn.empty()) return out;
    double th = 1e300;
    for (size_t c = 0; c < g.bulk.size(); ++c) if (active(g, c) && g.bot[c] > g.top[c]) th = std::min(th, g.bot[c] - g.top[c]);
    const double ds = std::max(0.05, std::min(th < 1e299 ? th / 4 : 1.0, std::min(g.spec.dx, g.spec.dy) / 8));
    for (const auto& w : wells) {
        const double L = w.totalMd();
        if (L <= 0 || w.logs.empty()) continue;
        struct Acc { double len = 0, mdSum = 0; std::map<std::string, std::vector<std::pair<double, double>>> vals; };
        std::map<int, Acc> acc;
        std::vector<int> order;
        for (double md = 0.5 * ds; md < L; md += ds) {
            const P3 p = w.at(md);
            const int c = g.cellOf(p[0], p[1], p[2]);
            if (c < 0 || !active(g, size_t(c))) continue;
            auto it = acc.find(c);
            if (it == acc.end()) { it = acc.emplace(c, Acc{}).first; order.push_back(c); }
            Acc& a = it->second;
            a.len += ds; a.mdSum += md * ds;
            for (const auto& lg : w.logs)
                for (const auto& cv : lg.curves) {
                    const double v = lg.valueAt(cv, md);
                    if (std::isfinite(v)) a.vals[upperS(cv.mnem)].push_back({v, ds});
                }
        }
        for (int c : order) {
            const Acc& a = acc[c];
            Blocked b; b.cell = c; b.well = w.name; b.md = a.len > 0 ? a.mdSum / a.len : 0;
            for (const auto& kv : a.vals) {
                int m = defaultMethod(kv.first);
                for (const auto& mm : methods) if (upperS(mm.first) == kv.first) m = mm.second;
                double num = 0, den = 0;
                if (m == 3) {
                    std::map<long, double> lenOf;
                    for (const auto& p : kv.second) lenOf[std::lround(p.first)] += p.second;
                    long best = 0; double bl = -1;
                    for (const auto& q : lenOf) if (q.second > bl) { bl = q.second; best = q.first; }
                    b.v[kv.first] = double(best);
                    continue;
                }
                for (const auto& p : kv.second) {
                    if (m == 1) { if (p.first > 0) { num += std::log(p.first) * p.second; den += p.second; } }
                    else if (m == 2) { if (p.first > 0) { num += p.second / p.first; den += p.second; } }
                    else { num += p.first * p.second; den += p.second; }
                }
                if (den <= 0) continue;
                b.v[kv.first] = m == 1 ? std::exp(num / den) : m == 2 ? den / num : num / den;
            }
            out.push_back(b);
        }
    }
    return out;
}

// ================================================================ facies
FaciesResult modelFacies(const GeoGrid& g, const std::vector<Blocked>& data, const FaciesSpec& spec) {
    FaciesResult R;
    const size_t n = g.bulk.size();
    const std::string cur = upperS(spec.curve);
    std::map<int, int> dataCode;   // cell -> code
    for (const auto& b : data) { auto it = b.v.find(cur); if (it != b.v.end() && b.cell >= 0) dataCode[b.cell] = int(std::lround(it->second)); }
    R.codes = spec.codes;
    if (R.codes.empty()) { std::set<int> s; for (auto& d : dataCode) s.insert(d.second); R.codes.assign(s.begin(), s.end()); }
    if (R.codes.empty()) R.codes = {0, 1};
    const size_t nc = R.codes.size();
    std::vector<double> prop = spec.proportions;
    if (prop.size() != nc) {
        prop.assign(nc, 0.0);
        for (auto& d : dataCode) for (size_t k = 0; k < nc; ++k) if (R.codes[k] == d.second) prop[k] += 1;
        double s = std::accumulate(prop.begin(), prop.end(), 0.0);
        if (s <= 0) { prop.assign(nc, 1.0 / double(nc)); s = 1; }
        for (double& p : prop) p /= s;
    }
    auto codeIdx = [&](int code) { for (size_t k = 0; k < nc; ++k) if (R.codes[k] == code) return int(k); return -1; };
    R.facies.assign(n, R.codes[0]);
    const auto X = stratCoords(g);
    std::mt19937 rng(spec.seed);
    std::uniform_real_distribution<double> U(0, 1);
    const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
    if (spec.method == 1) {   // truncated Gaussian: thresholds from the cumulative proportions
        std::vector<double> thr(nc + 1);
        thr[0] = -1e9; double cum = 0;
        for (size_t k = 0; k < nc; ++k) { cum += prop[k]; thr[k + 1] = k + 1 == nc ? 1e9 : invNorm(cum); }
        std::vector<geostat::Xyz> tg; std::vector<size_t> idx;
        for (size_t c = 0; c < n; ++c) if (active(g, c)) { tg.push_back(X[c]); idx.push_back(c); }
        std::vector<geostat::Point> d;
        for (auto& dc : dataCode) {
            const int k = codeIdx(dc.second);
            if (k < 0) continue;
            const double lo = std::max(thr[size_t(k)], -4.0), hi = std::min(thr[size_t(k) + 1], 4.0);
            const auto& p = X[size_t(dc.first)];
            d.push_back({p[0], p[1], p[2], 0.5 * (lo + hi)});
        }
        geostat::Options o; o.vg = spec.vg; o.vg.nugget = std::min(o.vg.nugget, 0.99); o.vg.sill = 1 - o.vg.nugget; o.maxNbr = spec.maxNbr; o.seed = spec.seed; o.normalScore = false;
        const auto y = geostat::sgs(tg, d, o);
        for (size_t q = 0; q < idx.size() && q < y.size(); ++q) {
            size_t k = 0;
            while (k + 1 < nc && y[q] >= thr[k + 1]) ++k;
            R.facies[idx[q]] = R.codes[k];
        }
    } else if (spec.method == 2 && nc >= 2) {   // object-based sinuous channels in the sim box
        double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300, zmax = 0;
        for (size_t c = 0; c < n; ++c) { xmin = std::min(xmin, X[c][0]); xmax = std::max(xmax, X[c][0]); ymin = std::min(ymin, X[c][1]); ymax = std::max(ymax, X[c][1]); zmax = std::max(zmax, X[c][2]); }
        const double target = prop[1];
        size_t nAct = 0; for (size_t c = 0; c < n; ++c) nAct += active(g, c);
        std::vector<char> ch(n, 0);
        size_t nch = 0;
        const double az = spec.chAzimuth * 3.14159265358979323846 / 180.0, ux = std::sin(az), uy = std::cos(az);   // flow direction (azimuth from north)
        const double cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax), half = 0.5 * std::hypot(xmax - xmin, ymax - ymin) + spec.chWidth;
        for (int t = 0; t < 400 && double(nch) < target * double(nAct); ++t) {
            const double off = (U(rng) * 2 - 1) * half, ph = U(rng) * 6.283185307, zTop = U(rng) * std::max(zmax - spec.chThick * 0.5, 0.0);
            for (size_t c = 0; c < n; ++c) {
                if (!active(g, c) || ch[c]) continue;
                const double rx = X[c][0] - cx, ry = X[c][1] - cy;
                const double s = rx * ux + ry * uy, l = -rx * uy + ry * ux;   // along / across the flow direction
                const double centre = off + spec.chAmplitude * std::sin(6.283185307 * s / std::max(1.0, spec.chWavelength) + ph);
                const double d = std::fabs(l - centre) / std::max(1e-9, 0.5 * spec.chWidth);
                if (d >= 1) continue;
                const double thick = spec.chThick * (1 - d * d);   // lens cross-section
                if (X[c][2] >= zTop && X[c][2] <= zTop + thick) { ch[c] = 1; ++nch; }
            }
        }
        for (size_t c = 0; c < n; ++c) R.facies[c] = ch[c] ? R.codes[1] : R.codes[0];
        for (auto& dc : dataCode) if (R.facies[size_t(dc.first)] != dc.second) { ++R.mismatches; R.facies[size_t(dc.first)] = dc.second; }
    } else {   // sequential indicator simulation (simple indicator kriging, one variogram for all indicators)
        std::vector<int> fac(n, -1);   // code index, -1 unknown
        for (auto& dc : dataCode) { const int k = codeIdx(dc.second); if (k >= 0 && active(g, size_t(dc.first))) fac[size_t(dc.first)] = k; }
        std::vector<size_t> path;
        for (size_t c = 0; c < n; ++c) if (active(g, c) && fac[c] < 0) path.push_back(c);
        std::shuffle(path.begin(), path.end(), rng);
        geostat::Variogram vg = spec.vg;
        const double ri = std::max(1, int(std::ceil(std::max(vg.rangeMax, vg.rangeMin) / g.spec.dx)));
        const double rj = std::max(1, int(std::ceil(std::max(vg.rangeMax, vg.rangeMin) / g.spec.dy)));
        double th = X.size() > size_t(nx * ny) ? X[size_t(nx * ny)][2] - X[0][2] : 1.0;
        const int wi = int(std::min(ri, 12.0)), wj = int(std::min(rj, 12.0)), wk = int(std::min(std::ceil(vg.rangeVert / std::max(th, 1e-9)), 6.0));
        const int maxN = std::max(4, spec.maxNbr);
        std::vector<double> A, b;
        for (size_t c : path) {
            const int i = int(c % size_t(nx)), j = int((c / size_t(nx)) % size_t(ny)), k = int(c / size_t(nx * ny));
            std::vector<std::pair<double, size_t>> nb;
            for (int kk = std::max(0, k - wk); kk <= std::min(nz - 1, k + wk); ++kk)
                for (int jj = std::max(0, j - wj); jj <= std::min(ny - 1, j + wj); ++jj)
                    for (int ii = std::max(0, i - wi); ii <= std::min(nx - 1, i + wi); ++ii) {
                        const size_t q = size_t(ii + nx * (jj + ny * kk));
                        if (fac[q] < 0) continue;
                        const double gm = geostat::gamma(vg, X[q][0] - X[c][0], X[q][1] - X[c][1], X[q][2] - X[c][2]);
                        if (gm >= vg.nugget + vg.sill - 1e-12) continue;   // beyond the range
                        nb.push_back({gm, q});
                    }
            std::vector<double> p = prop;
            if (!nb.empty()) {
                if (int(nb.size()) > maxN) { std::partial_sort(nb.begin(), nb.begin() + maxN, nb.end()); nb.resize(size_t(maxN)); }
                const int m = int(nb.size());
                A.assign(size_t(m * m), 0); b.assign(size_t(m), 0);
                for (int r = 0; r < m; ++r) {
                    const auto& pr = X[nb[size_t(r)].second];
                    b[size_t(r)] = geostat::covariance(vg, pr[0] - X[c][0], pr[1] - X[c][1], pr[2] - X[c][2]);
                    for (int s = 0; s < m; ++s) { const auto& ps = X[nb[size_t(s)].second]; A[size_t(r * m + s)] = geostat::covariance(vg, pr[0] - ps[0], pr[1] - ps[1], pr[2] - ps[2]); }
                }
                if (solveDense(A, b, m)) {
                    double s = 0;
                    for (size_t kf = 0; kf < nc; ++kf) {
                        double e = prop[kf];
                        for (int r = 0; r < m; ++r) e += b[size_t(r)] * ((fac[nb[size_t(r)].second] == int(kf) ? 1.0 : 0.0) - prop[kf]);
                        p[kf] = std::min(1.0, std::max(0.0, e)); s += p[kf];   // order-relation correction
                    }
                    if (s > 0) for (double& v : p) v /= s; else p = prop;
                }
            }
            double u = U(rng), cum = 0; int pick = int(nc) - 1;
            for (size_t kf = 0; kf < nc; ++kf) { cum += p[kf]; if (u <= cum) { pick = int(kf); break; } }
            fac[c] = pick;
        }
        for (size_t c = 0; c < n; ++c) if (fac[c] >= 0) R.facies[c] = R.codes[size_t(fac[c])];
    }
    // honoured data / proportions
    for (auto& dc : dataCode) if (R.facies[size_t(dc.first)] == dc.second) ++R.conditioned;
    R.proportions.assign(nc, 0.0);
    double na = 0;
    for (size_t c = 0; c < n; ++c) if (active(g, c)) { const int k = codeIdx(R.facies[c]); if (k >= 0) R.proportions[size_t(k)] += 1; na += 1; }
    if (na > 0) for (double& p : R.proportions) p /= na;
    return R;
}

// ================================================================ petrophysics
PetroResult modelProperty(const GeoGrid& g, const std::vector<Blocked>& data, const std::vector<int>* facies, const PetroSpec& spec) {
    PetroResult R;
    const size_t n = g.bulk.size();
    R.v.assign(n, kNaN);
    const std::string cur = upperS(spec.curve);
    const auto X = stratCoords(g);
    struct D { size_t cell; double v; };
    std::vector<D> all;
    for (const auto& b : data) {
        auto it = b.v.find(cur);
        if (it == b.v.end() || b.cell < 0 || !std::isfinite(it->second)) continue;
        double v = it->second;
        if (spec.logTransform) { if (v <= 0) continue; v = std::log10(v); }
        all.push_back({size_t(b.cell), v});
    }
    R.data = int(all.size());
    if (all.empty()) { R.note = "no blocked data for " + spec.curve; return R; }
    const bool sec = spec.secondary.size() == n && std::fabs(spec.corrSecondary) > 1e-6;
    std::vector<double> sNs;
    if (sec) {   // normal scores of the secondary over the active cells (rank transform)
        std::vector<std::pair<double, size_t>> r;
        for (size_t c = 0; c < n; ++c) if (active(g, c) && std::isfinite(spec.secondary[c])) r.push_back({spec.secondary[c], c});
        std::sort(r.begin(), r.end());
        sNs.assign(n, 0.0);
        for (size_t k = 0; k < r.size(); ++k) sNs[r[k].second] = invNorm((k + 0.5) / double(r.size()));
    }
    // groups: per facies code (when enough data) or one group
    std::map<int, std::vector<size_t>> cellsOf;   // group -> cells
    std::map<int, std::vector<D>> dataOf;
    const bool perF = spec.perFacies && facies && facies->size() == n;
    for (size_t c = 0; c < n; ++c) if (active(g, c)) cellsOf[perF ? (*facies)[c] : 0].push_back(c);
    for (const auto& d : all) dataOf[perF ? (*facies)[d.cell] : 0].push_back(d);
    int gi = 0;
    for (auto& grp : cellsOf) {
        std::vector<D> dd = dataOf[grp.first];
        std::string tag = perF ? "facies " + std::to_string(grp.first) : "all cells";
        if (dd.size() < 3) { dd = all; tag += " (too few data - all data used)"; }
        NScore ns; { std::vector<double> v; for (auto& d : dd) v.push_back(d.v); ns.build(v); }
        std::vector<geostat::Xyz> tg; for (size_t c : grp.second) tg.push_back(X[c]);
        const double rho = sec ? std::max(-0.99, std::min(0.99, spec.corrSecondary)) : 0.0;
        std::vector<geostat::Point> pts;
        for (auto& d : dd) { const auto& p = X[d.cell]; pts.push_back({p[0], p[1], p[2], ns.fwd(d.v) - (sec ? rho * sNs[d.cell] : 0.0)}); }
        geostat::Options o; o.vg = spec.vg; o.maxNbr = spec.maxNbr; o.seed = spec.seed + unsigned(97 * gi++); o.normalScore = false;
        const double tot = std::max(1e-9, o.vg.nugget + o.vg.sill), fr = (1 - rho * rho) / tot;   // residual variance 1 - rho^2
        o.vg.nugget *= fr; o.vg.sill *= fr;
        std::vector<double> y;
        if (spec.method == 0) { o.simpleKriging = true; o.skMean = 0; y = geostat::krige(tg, pts, o).est; }
        else y = geostat::sgs(tg, pts, o);
        for (size_t q = 0; q < grp.second.size() && q < y.size(); ++q) {
            const size_t c = grp.second[q];
            double gsc = y[q] + (sec ? rho * sNs[c] : 0.0);
            double v = ns.back(gsc);
            if (spec.logTransform) v = std::pow(10.0, v);
            R.v[c] = std::max(spec.minV, std::min(spec.maxV, v));
        }
        char b[160]; std::snprintf(b, sizeof b, "%s: %d data; ", tag.c_str(), int(dd.size()));
        R.note += b;
    }
    double s = 0, s2 = 0; int m = 0;
    for (size_t c = 0; c < n; ++c) if (std::isfinite(R.v[c])) { s += R.v[c]; s2 += R.v[c] * R.v[c]; ++m; }
    if (m) { R.mean = s / m; R.sd = std::sqrt(std::max(0.0, s2 / m - R.mean * R.mean)); }
    for (size_t c = 0; c < n; ++c) if (!std::isfinite(R.v[c])) R.v[c] = spec.logTransform ? std::max(spec.minV, 0.0) : std::max(spec.minV, std::min(spec.maxV, R.mean));
    return R;
}

std::vector<PoroPerm> fitPoroPerm(const std::vector<Blocked>& data, const std::string& phi, const std::string& perm, const std::vector<int>* facies,
                                  const std::vector<int>& codes) {
    const std::string P = upperS(phi), K = upperS(perm);
    auto fit = [&](int code, bool any) {
        PoroPerm f; double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0; int n = 0;
        for (const auto& b : data) {
            auto a = b.v.find(P), k = b.v.find(K);
            if (a == b.v.end() || k == b.v.end() || !(k->second > 0) || !std::isfinite(a->second)) continue;
            if (!any && facies && b.cell >= 0 && size_t(b.cell) < facies->size() && (*facies)[size_t(b.cell)] != code) continue;
            const double x = a->second, y = std::log10(k->second);
            sx += x; sy += y; sxx += x * x; sxy += x * y; syy += y * y; ++n;
        }
        f.n = n;
        if (n >= 2) {
            const double den = n * sxx - sx * sx;
            f.b = std::fabs(den) > 1e-30 ? (n * sxy - sx * sy) / den : 0.0;
            f.a = (sy - f.b * sx) / n;
            const double ssT = syy - sy * sy / n, ssR = ssT - f.b * (sxy - sx * sy / n);
            f.r2 = ssT > 0 ? 1 - ssR / ssT : 0.0;
        } else if (n == 1) f.a = sy;
        return f;
    };
    std::vector<PoroPerm> out;
    if (!facies || codes.empty()) { out.push_back(fit(0, true)); return out; }
    const PoroPerm global = fit(0, true);
    for (int c : codes) { PoroPerm f = fit(c, false); out.push_back(f.n >= 3 ? f : global); }
    return out;
}
std::vector<double> permFromPoro(const std::vector<double>& phi, const std::vector<PoroPerm>& fits, const std::vector<int>* facies,
                                 const std::vector<int>& codes, double scatterSd, unsigned seed) {
    std::vector<double> k(phi.size(), 0.0);
    std::mt19937 rng(seed);
    std::normal_distribution<double> N(0, 1);
    for (size_t c = 0; c < phi.size(); ++c) {
        size_t f = 0;
        if (facies && c < facies->size()) for (size_t q = 0; q < codes.size() && q < fits.size(); ++q) if (codes[q] == (*facies)[c]) f = q;
        if (fits.empty()) continue;
        const PoroPerm& P = fits[std::min(f, fits.size() - 1)];
        k[c] = std::pow(10.0, P.a + P.b * phi[c] + (scatterSd > 0 ? scatterSd * N(rng) : 0.0));
    }
    return k;
}

// ================================================================ volumetrics
VolResult volumetrics(const GeoGrid& g, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<int>* facies, const VolSpec& V) {
    VolResult R;
    const int nx = g.grid.nx, ny = g.grid.ny;
    std::map<std::string, VolRow> zr, fr, br;
    const double rbPerFt3 = 1.0 / 5.614583333;
    for (size_t c = 0; c < g.bulk.size(); ++c) {
        if (!active(g, c)) continue;
        if (!V.mask.empty() && (c >= V.mask.size() || !V.mask[c])) continue;
        const size_t colC = c % size_t(nx * ny);
        const int blkC = colC < g.block.size() ? g.block[colC] : 0;
        const Contacts& CT = blkC >= 0 && size_t(blkC) < V.blockContacts.size() ? V.blockContacts[size_t(blkC)] : V.contacts;
        const double phi = c < poro.size() ? poro[c] : 0.0, nt = c < ntg.size() ? ntg[c] : 1.0;
        VolRow r;
        r.bulk = g.bulk[c]; r.net = r.bulk * nt; r.pore = r.net * phi;
        const double z = g.center[c][2];
        double sw = 1.0;
        if (z < CT.owc) {
            if (V.satHeight) {
                const double pc = V.sh.drhoOW * (CT.owc - z);
                sw = pc > V.sh.pd ? V.sh.swirr + (1 - V.sh.swirr) * std::pow(V.sh.pd / pc, V.sh.lambda) : 1.0;
            } else sw = V.swConst;
        }
        const double hc = r.pore * (1 - sw);
        if (z < CT.goc) r.hcpvG = hc; else r.hcpvO = hc;
        if (V.field) { r.stoiip = r.hcpvO * rbPerFt3 / std::max(V.bo, 1e-12); r.giip = r.hcpvG / std::max(V.bg, 1e-12); }
        else { r.stoiip = r.hcpvO / std::max(V.bo, 1e-12); r.giip = r.hcpvG / std::max(V.bg, 1e-12); }
        auto add = [&](VolRow& a) { a.bulk += r.bulk; a.net += r.net; a.pore += r.pore; a.hcpvO += r.hcpvO; a.hcpvG += r.hcpvG; a.stoiip += r.stoiip; a.giip += r.giip; };
        const int k = int(c / size_t(nx * ny));
        const std::string zk = k < int(g.layerZone.size()) && size_t(g.layerZone[size_t(k)]) < g.spec.zones.size() ? g.spec.zones[size_t(g.layerZone[size_t(k)])].name : "zone";
        add(zr[zk]); add(R.total);
        if (facies && c < facies->size()) add(fr["facies " + std::to_string((*facies)[c])]);
        const size_t col = c % size_t(nx * ny);
        add(br["block " + std::to_string(col < g.block.size() ? g.block[col] + 1 : 1)]);
    }
    for (auto& kv : zr) { kv.second.key = kv.first; R.byZone.push_back(kv.second); }
    for (auto& kv : fr) { kv.second.key = kv.first; R.byFacies.push_back(kv.second); }
    for (auto& kv : br) { kv.second.key = kv.first; R.byBlock.push_back(kv.second); }
    R.total.key = "total";
    return R;
}
Percentiles percentiles(std::vector<double> v) {
    Percentiles P;
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
    if (v.empty()) return P;
    std::sort(v.begin(), v.end());
    auto q = [&](double f) { const double x = f * double(v.size() - 1); const size_t i = size_t(x); const double t = x - double(i); return i + 1 < v.size() ? v[i] * (1 - t) + v[i + 1] * t : v.back(); };
    P.p90 = q(0.10); P.p50 = q(0.50); P.p10 = q(0.90);
    P.mean = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
    return P;
}

}  // namespace geo
}  // namespace ressim
