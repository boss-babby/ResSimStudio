// C43 geomodel tools: mis-tie correction, polygons, layer-cake velocity model, seismic calculator / cropping
#include "ressim/geo_tools.h"
#include "ressim/geo_analysis.h"
#include "ressim/geostat.h"
#include "ressim/seismic3d.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <sstream>

namespace ressim {
namespace geo {
namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();

// dense symmetric solve (Gaussian elimination with partial pivoting); false if singular
bool solveDense(std::vector<double> A, std::vector<double> b, int n, std::vector<double>& x) {
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
    x.assign(size_t(n), 0.0);
    for (int r = n - 1; r >= 0; --r) {
        double s = b[size_t(r)];
        for (int k = r + 1; k < n; ++k) s -= A[size_t(r * n + k)] * x[size_t(k)];
        x[size_t(r)] = s / A[size_t(r * n + r)];
    }
    return true;
}
double pickAt(const LineHorizon& h, int line, int tr) {
    if (line < 0 || size_t(line) >= h.t.size() || tr < 0 || size_t(tr) >= h.t[size_t(line)].size()) return kNaN;
    return h.t[size_t(line)][size_t(tr)];
}
std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
geostat::Variogram idw() { geostat::Variogram v; v.rangeMax = v.rangeMin = 0; return v; }
}  // namespace

// ================================================================ P6 mis-ties
MistieFit fitMisties(const std::vector<SeismicLine>& lines, const std::vector<LineTie>& ties, const std::vector<LineHorizon>& hz, int refLine) {
    MistieFit F;
    const int n = int(lines.size());
    F.shift.assign(size_t(n), 0.0);
    F.component.assign(size_t(n), -1);
    struct Pair { int a, b; double d; };   // d = t_a - t_b
    std::vector<Pair> P;
    for (const auto& t : ties)
        for (const auto& h : hz) {
            const double ta = pickAt(h, t.a, t.ta), tb = pickAt(h, t.b, t.tb);
            if (std::isfinite(ta) && std::isfinite(tb) && t.a >= 0 && t.b >= 0 && t.a < n && t.b < n) P.push_back({t.a, t.b, ta - tb});
        }
    F.pairs = int(P.size());
    if (P.empty() || n == 0) return F;
    // connected components over the pairs
    std::vector<int> par(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) par[size_t(i)] = i;
    std::function<int(int)> root = [&](int i) { return par[size_t(i)] == i ? i : par[size_t(i)] = root(par[size_t(i)]); };
    std::vector<char> used(size_t(n), 0);
    for (const auto& p : P) { par[size_t(root(p.a))] = root(p.b); used[size_t(p.a)] = used[size_t(p.b)] = 1; }
    std::map<int, int> compId;
    for (int i = 0; i < n; ++i) if (used[size_t(i)]) { const int r = root(i); if (!compId.count(r)) { const int id = int(compId.size()); compId[r] = id; } F.component[size_t(i)] = compId[r]; }
    // reference line per component
    std::vector<int> ref(compId.size(), -1);
    for (int i = 0; i < n; ++i) if (F.component[size_t(i)] >= 0 && ref[size_t(F.component[size_t(i)])] < 0) ref[size_t(F.component[size_t(i)])] = i;
    if (refLine >= 0 && refLine < n && F.component[size_t(refLine)] >= 0) ref[size_t(F.component[size_t(refLine)])] = refLine;
    // normal equations over the unknown lines (used and not a reference)
    std::vector<int> idx(size_t(n), -1);
    int m = 0;
    for (int i = 0; i < n; ++i) {
        const int c = F.component[size_t(i)];
        if (c >= 0 && ref[size_t(c)] != i) idx[size_t(i)] = m++;
    }
    double sb = 0;
    for (const auto& p : P) sb += p.d * p.d;
    F.rmsBefore = std::sqrt(sb / double(P.size()));
    if (m > 0) {
        std::vector<double> A(size_t(m) * size_t(m), 0.0), b(size_t(m), 0.0), x;
        // residual r = d + s_a - s_b ; d/ds_a = +1, d/ds_b = -1
        for (const auto& p : P) {
            const int ia = idx[size_t(p.a)], ib = idx[size_t(p.b)];
            if (ia >= 0) { A[size_t(ia * m + ia)] += 1; b[size_t(ia)] -= p.d; }
            if (ib >= 0) { A[size_t(ib * m + ib)] += 1; b[size_t(ib)] += p.d; }
            if (ia >= 0 && ib >= 0) { A[size_t(ia * m + ib)] -= 1; A[size_t(ib * m + ia)] -= 1; }
        }
        if (solveDense(A, b, m, x))
            for (int i = 0; i < n; ++i) if (idx[size_t(i)] >= 0) F.shift[size_t(i)] = x[size_t(idx[size_t(i)])];
    }
    double sa = 0;
    for (const auto& p : P) { const double r = p.d + F.shift[size_t(p.a)] - F.shift[size_t(p.b)]; sa += r * r; }
    F.rmsAfter = std::sqrt(sa / double(P.size()));
    return F;
}

void applyLineShifts(std::vector<SeismicLine>& lines, std::vector<LineHorizon>& hz, const std::vector<double>& shift) {
    for (size_t l = 0; l < lines.size() && l < shift.size(); ++l) {
        if (shift[l] == 0) continue;
        lines[l].t0 += shift[l];
        for (auto& h : hz)
            if (l < h.t.size()) for (double& t : h.t[l]) if (std::isfinite(t)) t += shift[l];
    }
}

// ================================================================ P4 polygons
bool insidePolygon(const Polygon& p, double x, double y) {
    if (!p.closed || p.pts.size() < 3) return false;
    bool in = false;
    for (size_t i = 0, j = p.pts.size() - 1; i < p.pts.size(); j = i++) {
        const double xi = p.pts[i][0], yi = p.pts[i][1], xj = p.pts[j][0], yj = p.pts[j][1];
        if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) in = !in;
    }
    return in;
}
double polygonArea(const Polygon& p) {
    double a = 0;
    for (size_t i = 0, j = p.pts.size() ? p.pts.size() - 1 : 0; i < p.pts.size(); j = i++) a += p.pts[j][0] * p.pts[i][1] - p.pts[i][0] * p.pts[j][1];
    return 0.5 * std::fabs(a);
}
double polylineLength(const Polygon& p) {
    double L = 0;
    for (size_t i = 1; i < p.pts.size(); ++i) L += std::hypot(p.pts[i][0] - p.pts[i - 1][0], p.pts[i][1] - p.pts[i - 1][1]);
    if (p.closed && p.pts.size() > 2) L += std::hypot(p.pts.front()[0] - p.pts.back()[0], p.pts.front()[1] - p.pts.back()[1]);
    return L;
}
std::vector<Polygon> parsePolygons(const std::string& text) {
    std::vector<Polygon> out;
    Polygon cur;
    std::string pendingName, lastId;
    auto flush = [&]() {
        if (cur.pts.size() >= 2) {
            if (cur.name.empty()) cur.name = "Polygon " + std::to_string(out.size() + 1);
            // closed when the first point is repeated at the end (dropped) or when it has >= 3 points
            if (cur.pts.size() > 2 && std::fabs(cur.pts.front()[0] - cur.pts.back()[0]) < 1e-9 && std::fabs(cur.pts.front()[1] - cur.pts.back()[1]) < 1e-9) cur.pts.pop_back();
            cur.closed = cur.pts.size() >= 3;
            out.push_back(cur);
        }
        cur = Polygon{};
        cur.name = pendingName; pendingName.clear();
    };
    std::istringstream is(text);
    std::string line;
    while (std::getline(is, line)) {
        const std::string t = trim(line);
        if (t.empty() || t == ">" || t == "999" || t.rfind("999 ", 0) == 0 || t.rfind("999\t", 0) == 0) { flush(); continue; }
        if (t[0] == '#' || t[0] == '!') {
            std::string nm = trim(t.substr(1));
            if (cur.pts.empty()) cur.name = nm; else { flush(); cur.name = nm; }
            continue;
        }
        std::string s = t;
        for (char& c : s) if (c == ',' || c == ';') c = ' ';
        std::istringstream ls(s);
        double x, y;
        if (!(ls >> x >> y)) continue;
        if (std::fabs(x - 999) < 1e-9 && std::fabs(y - 999) < 1e-9) { flush(); continue; }
        double z;
        std::string id;
        if (ls >> z) { ls >> id; }
        if (!id.empty() && id != lastId && !lastId.empty()) flush();
        if (!id.empty()) { lastId = id; if (cur.name.empty()) cur.name = id; }
        cur.pts.push_back({x, y});
    }
    flush();
    return out;
}
std::vector<Polygon> readPolygons(const std::string& path, std::string& err) {
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return {}; }
    std::stringstream ss; ss << f.rdbuf();
    auto v = parsePolygons(ss.str());
    if (v.empty()) err = path + ": no polygon (need x y rows)";
    return v;
}
bool writePolygons(const std::string& path, const std::vector<Polygon>& polys) {
    std::ofstream f(path);
    if (!f) return false;
    for (const auto& p : polys) {
        f << "# " << p.name << "\n";
        char b[96];
        for (const auto& q : p.pts) { std::snprintf(b, sizeof b, "%.3f %.3f\n", q[0], q[1]); f << b; }
        if (p.closed && !p.pts.empty()) { std::snprintf(b, sizeof b, "%.3f %.3f\n", p.pts[0][0], p.pts[0][1]); f << b; }
        f << "999 999\n";
    }
    return bool(f);
}
Surface clipSurface(const Surface& s, const Polygon& p, bool keepInside) {
    Surface o = s;
    for (int j = 0; j < s.ny; ++j)
        for (int i = 0; i < s.nx; ++i) {
            double x, y; s.nodeXY(i, j, x, y);
            if (insidePolygon(p, x, y) != keepInside) o.z[size_t(i + s.nx * j)] = kNaN;
        }
    return o;
}
std::vector<char> polygonMask(const GeoGrid& g, const Polygon& p) {
    std::vector<char> m(g.center.size(), 0);
    for (size_t c = 0; c < g.center.size(); ++c) m[c] = insidePolygon(p, g.center[c][0], g.center[c][1]) ? 1 : 0;
    return m;
}

// ================================================================ P8 layer-cake velocity
Surface timeSurface(const SeismicCube& c, const TimeHorizon& h, const std::string& name) {
    // same lattice convention as horizonToDepth (seismic.cpp), z = TWT
    Surface s;
    s.name = name; s.nx = c.ni; s.ny = c.nj;
    s.x0 = c.x0; s.y0 = c.y0;
    s.dx = std::hypot(c.ix, c.iy); s.dy = std::hypot(c.jx, c.jy);
    s.rotDeg = std::atan2(c.iy, c.ix) * 180.0 / 3.14159265358979323846;
    s.z.assign(size_t(s.nx) * size_t(s.ny), kNaN);
    for (int i = 0; i < c.ni; ++i)
        for (int j = 0; j < c.nj; ++j) {
            const size_t k = size_t(i) * size_t(c.nj) + size_t(j);
            if (k < h.t.size()) s.z[size_t(i + s.nx * j)] = h.t[k];
        }
    if (c.ix * c.jy - c.iy * c.jx < 0) {   // left-handed crossline axis: mirror into the Surface frame
        Surface f = s; f.x0 = c.x0 + (c.nj - 1) * c.jx; f.y0 = c.y0 + (c.nj - 1) * c.jy;
        for (int i = 0; i < s.nx; ++i) for (int j = 0; j < s.ny; ++j) f.z[size_t(i + s.nx * j)] = s.z[size_t(i + s.nx * (s.ny - 1 - j))];
        return f;
    }
    return s;
}
std::vector<P3> lineHorizonTimePoints(const std::vector<SeismicLine>& lines, const LineHorizon& h, int step) {
    std::vector<P3> pts;
    step = std::max(1, step);
    for (size_t l = 0; l < lines.size() && l < h.t.size(); ++l)
        for (int t = 0; t < lines[l].ntr() && size_t(t) < h.t[l].size(); t += step)
            if (std::isfinite(h.t[l][size_t(t)])) pts.push_back({lines[l].x[size_t(t)], lines[l].y[size_t(t)], h.t[l][size_t(t)]});
    return pts;
}

LayerCakeResult layerCakeDepth(const LayerCakeInput& in, const std::vector<GeoWell>& wells) {
    LayerCakeResult R;
    const size_t nI = in.twt.size();
    if (nI == 0 || !in.twt[0].valid()) { R.warnings.push_back("no time surfaces"); return R; }
    const Surface& L = in.twt[0];   // output lattice
    const size_t nn = size_t(L.nx) * size_t(L.ny);
    std::vector<double> X(nn), Y(nn);
    for (int j = 0; j < L.ny; ++j) for (int i = 0; i < L.nx; ++i) L.nodeXY(i, j, X[size_t(i + L.nx * j)], Y[size_t(i + L.nx * j)]);
    // time of every interface on the output lattice
    std::vector<std::vector<double>> T(nI, std::vector<double>(nn, kNaN));
    for (size_t f = 0; f < nI; ++f) {
        if (!in.twt[f].valid()) { R.warnings.push_back("time surface " + std::to_string(f + 1) + " is not valid"); return R; }
        for (size_t k = 0; k < nn; ++k) T[f][k] = f == 0 ? L.z[k] : in.twt[f].at(X[k], Y[k]);
    }
    auto tAt = [&](size_t f, double x, double y) { return f == 0 ? L.at(x, y) : in.twt[f].at(x, y); };
    // well tops per interface: position (x, y, depth)
    struct WT { std::string well; double x, y, z; };
    std::vector<std::vector<WT>> tops(nI);
    auto upper = [](std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; };
    for (size_t f = 0; f < nI; ++f) {
        const std::string want = f < in.topName.size() ? upper(in.topName[f]) : std::string();
        if (want.empty()) continue;
        for (const auto& w : wells)
            for (const auto& t : w.tops)
                if (upper(t.surface) == want) { const P3 p = w.at(t.md); tops[f].push_back({w.name, p[0], p[1], p[2]}); }
    }
    // interval velocities per layer (layer f = between interface f-1 (datum for f = 0) and f)
    R.vint.resize(nI); R.vMean.assign(nI, 0.0); R.nWells.assign(nI, 0);
    std::vector<std::vector<double>> V(nI, std::vector<double>(nn, kNaN));
    for (size_t f = 0; f < nI; ++f) {
        std::vector<P3> vpts;
        // wells with the base top; the top of the layer from the same well (or the datum)
        for (const auto& b : tops[f]) {
            double zTop = in.datumZ, tTop = 0;
            if (f > 0) {
                bool found = false;
                for (const auto& a : tops[f - 1]) if (a.well == b.well) { zTop = a.z; found = true; }
                if (!found) continue;
                tTop = tAt(f - 1, b.x, b.y);
            }
            const double tBase = tAt(f, b.x, b.y);
            if (!std::isfinite(tTop) || !std::isfinite(tBase) || tBase - tTop < 1e-6) continue;
            const double v = 2000.0 * (b.z - zTop) / (tBase - tTop);
            if (v > 0 && std::isfinite(v)) vpts.push_back({b.x, b.y, v});
        }
        R.nWells[f] = int(vpts.size());
        double mean = 0;
        for (const auto& p : vpts) mean += p[2];
        if (!vpts.empty()) mean /= double(vpts.size());
        R.vMean[f] = mean;
        const int mode = f < in.mode.size() ? in.mode[f] : 0;
        const double vu = f < in.vUser.size() && in.vUser[f] > 0 ? in.vUser[f] : 0.0;
        Surface vs = L; vs.name = "Vint " + std::to_string(f + 1);
        if (mode == 2 || vpts.empty()) {
            const double v = vu > 0 ? vu : (mean > 0 ? mean : 2000.0);
            if (vpts.empty() && mode != 2) R.warnings.push_back("layer " + std::to_string(f + 1) + ": no well calibration - interval velocity " + std::to_string(int(v)));
            vs.z.assign(nn, v);
        } else if (mode == 1 && vpts.size() >= 2) {
            vs = gridPoints(vpts, L, idw(), 12);
            vs.name = "Vint " + std::to_string(f + 1);
            for (double& v : vs.z) if (!std::isfinite(v)) v = mean;
        } else vs.z.assign(nn, mean);
        V[f] = vs.z;
        R.vint[f] = vs;
    }
    // depth conversion
    R.depth.resize(nI);
    std::vector<double> zPrev(nn, in.datumZ), tPrev(nn, 0.0);
    for (size_t f = 0; f < nI; ++f) {
        Surface d = L;
        d.name = f < in.topName.size() && !in.topName[f].empty() ? in.topName[f] : (in.twt[f].name.empty() ? "Depth " + std::to_string(f + 1) : in.twt[f].name);
        for (size_t k = 0; k < nn; ++k) {
            const double dt = T[f][k] - tPrev[k];
            d.z[k] = std::isfinite(dt) && std::isfinite(zPrev[k]) ? zPrev[k] + V[f][k] * std::max(0.0, dt) / 2000.0 : kNaN;
        }
        zPrev = d.z; tPrev = T[f];
        R.depth[f] = d;
    }
    // residuals at the wells (before), residual correction (cumulative: interface f corrected by its own residual map)
    double sb = 0, sa = 0; int nb = 0;
    for (size_t f = 0; f < nI; ++f) {
        std::vector<P3> res;
        for (const auto& t : tops[f]) {
            const double zm = R.depth[f].at(t.x, t.y);
            if (!std::isfinite(zm)) continue;
            R.ties.push_back({t.well, int(f), t.x, t.y, t.z, zm, zm});
            res.push_back({t.x, t.y, t.z - zm});
            sb += (t.z - zm) * (t.z - zm); ++nb;
        }
        if (in.residual && !res.empty()) {
            Surface rs = gridPoints(res, R.depth[f], idw(), 12);
            for (size_t k = 0; k < nn; ++k) if (std::isfinite(R.depth[f].z[k]) && std::isfinite(rs.z[k])) R.depth[f].z[k] += rs.z[k];
        }
    }
    // keep the order of the interfaces after the correction (no negative isochores)
    for (size_t f = 1; f < nI; ++f)
        for (size_t k = 0; k < nn; ++k)
            if (std::isfinite(R.depth[f].z[k]) && std::isfinite(R.depth[f - 1].z[k])) R.depth[f].z[k] = std::max(R.depth[f].z[k], R.depth[f - 1].z[k]);
    for (auto& t : R.ties) { t.zAfter = R.depth[size_t(t.iface)].at(t.x, t.y); if (std::isfinite(t.zAfter)) sa += (t.zTop - t.zAfter) * (t.zTop - t.zAfter); }
    if (nb) { R.rmsBefore = std::sqrt(sb / nb); R.rmsAfter = std::sqrt(sa / nb); }
    return R;
}

VelocityModel fitCheckshots(const std::vector<double>& z, const std::vector<double>& twt, double datum, double* rmsMs) {
    VelocityModel best; best.datum = datum;
    const size_t n = std::min(z.size(), twt.size());
    if (n == 0) { if (rmsMs) *rmsMs = 0; return best; }
    // for a given k the TWT is monotone in v0: t(z) = 2000/k ln(1 + k (z - datum)/v0) (k > 0) - fit v0 by bisection
    // on the mean residual, then k by golden section on the rms
    auto timeOf = [&](double v0, double k, double zz) {
        const double h = std::max(0.0, zz - datum);
        return std::fabs(k) < 1e-9 ? 2000.0 * h / v0 : 2000.0 / k * std::log(1.0 + k * h / v0);
    };
    auto fitV0 = [&](double k) {
        double lo = 100, hi = 20000;
        for (int it = 0; it < 80; ++it) {
            const double mid = 0.5 * (lo + hi);
            double s = 0;
            for (size_t i = 0; i < n; ++i) s += timeOf(mid, k, z[i]) - twt[i];
            if (s > 0) lo = mid; else hi = mid;   // too slow (late) -> increase v0
        }
        return 0.5 * (lo + hi);
    };
    auto rms = [&](double k, double& v0) {
        v0 = fitV0(k);
        double s = 0;
        for (size_t i = 0; i < n; ++i) { const double r = timeOf(v0, k, z[i]) - twt[i]; s += r * r; }
        return std::sqrt(s / double(n));
    };
    double a = 0.0, b = 3.0;
    const double g = 0.5 * (std::sqrt(5.0) - 1);
    double c = b - g * (b - a), d = a + g * (b - a), vc, vd;
    double fc = rms(c, vc), fd = rms(d, vd);
    for (int it = 0; it < 60; ++it) {
        if (fc < fd) { b = d; d = c; fd = fc; vd = vc; c = b - g * (b - a); fc = rms(c, vc); }
        else { a = c; c = d; fc = fd; vc = vd; d = a + g * (b - a); fd = rms(d, vd); }
    }
    double v0;
    const double k = 0.5 * (a + b);
    const double r = rms(k, v0);
    double v00; const double r0 = rms(0.0, v00);   // the boundary k = 0 (constant velocity)
    if (r0 <= r) { best.v0 = v00; best.k = 0; if (rmsMs) *rmsMs = r0; }
    else { best.v0 = v0; best.k = k; if (rmsMs) *rmsMs = r; }
    return best;
}

// ================================================================ P9 seismic calculator / cropping
namespace {
// expression names of the attributes, in SeisAttr order (index 0 = amplitude)
const char* kAttrVar[] = {"AMP", "ENVELOPE", "PHASE", "FREQUENCY", "COSPHASE", "SWEETNESS", "RMS", "RELIMPEDANCE", "COHERENCE", "VARIANCE"};
static_assert(sizeof(kAttrVar) / sizeof(kAttrVar[0]) == size_t(SeisAttr::Count), "attribute variable names");
}  // namespace
std::vector<std::string> seismicCalculatorNames() {
    std::vector<std::string> v(std::begin(kAttrVar), std::end(kAttrVar));
    v.push_back("T"); v.push_back("I"); v.push_back("J");
    return v;
}

bool seismicCalculator(const SeismicCube& c, const std::string& expr, int winSamples, SeismicCube& out, std::string& err) {
    if (c.ni <= 0 || c.nj <= 0 || c.ns <= 0 || c.a.size() != size_t(c.ni) * size_t(c.nj) * size_t(c.ns)) { err = "no seismic volume"; return false; }
    // which attributes does the expression use? (whole-word, case-insensitive)
    std::string U = expr;
    for (char& ch : U) ch = char(std::toupper((unsigned char)ch));
    auto uses = [&](const std::string& w) {
        for (size_t p = U.find(w); p != std::string::npos; p = U.find(w, p + 1)) {
            const bool l = p == 0 || !(std::isalnum((unsigned char)U[p - 1]) || U[p - 1] == '_');
            const bool r = p + w.size() >= U.size() || !(std::isalnum((unsigned char)U[p + w.size()]) || U[p + w.size()] == '_');
            if (l && r) return true;
        }
        return false;
    };
    std::vector<std::pair<std::string, SeismicCube>> attr;
    for (int a = 1; a < int(SeisAttr::Count); ++a) {
        const std::string nm = kAttrVar[a];
        if (uses(nm)) attr.emplace_back(nm, cubeAttribute(c, SeisAttr(a), std::max(1, winSamples)));
    }
    out = c;
    // evaluate one inline at a time (bounded double memory)
    const size_t per = size_t(c.nj) * size_t(c.ns);
    std::vector<double> amp(per), T(per), I(per), J(per), r;
    std::vector<std::vector<double>> av(attr.size(), std::vector<double>(per));
    std::map<std::string, const std::vector<double>*> vars = {{"AMP", &amp}, {"T", &T}, {"I", &I}, {"J", &J}};
    for (size_t a = 0; a < attr.size(); ++a) vars[attr[a].first] = &av[a];
    for (int i = 0; i < c.ni; ++i) {
        const size_t off = size_t(i) * per;
        for (size_t k = 0; k < per; ++k) {
            amp[k] = c.a[off + k];
            T[k] = c.t0 + c.dt * double(k % size_t(c.ns));
            I[k] = i; J[k] = double(k / size_t(c.ns));
        }
        for (size_t a = 0; a < attr.size(); ++a) for (size_t k = 0; k < per; ++k) av[a][k] = attr[a].second.a[off + k];
        if (!evalExpression(expr, vars, per, r, err)) return false;
        for (size_t k = 0; k < per; ++k) out.a[off + k] = std::isfinite(r[k]) ? float(r[k]) : 0.0f;
    }
    out.name = c.name + " (calc)";
    return true;
}

SeismicCube cropCube(const SeismicCube& c, int i0, int i1, int j0, int j1, int s0, int s1, int di, int dj, int ds) {
    SeismicCube o;
    i0 = std::clamp(i0, 0, std::max(0, c.ni - 1)); i1 = std::clamp(i1, i0, std::max(0, c.ni - 1));
    j0 = std::clamp(j0, 0, std::max(0, c.nj - 1)); j1 = std::clamp(j1, j0, std::max(0, c.nj - 1));
    s0 = std::clamp(s0, 0, std::max(0, c.ns - 1)); s1 = std::clamp(s1, s0, std::max(0, c.ns - 1));
    di = std::max(1, di); dj = std::max(1, dj); ds = std::max(1, ds);
    o.name = c.name + " (crop)"; o.textHeader = c.textHeader; o.format = c.format;
    o.ni = (i1 - i0) / di + 1; o.nj = (j1 - j0) / dj + 1; o.ns = (s1 - s0) / ds + 1;
    o.il0 = c.il0 + i0 * c.dil; o.xl0 = c.xl0 + j0 * c.dxl; o.dil = c.dil * di; o.dxl = c.dxl * dj;
    o.t0 = c.t0 + s0 * c.dt; o.dt = c.dt * ds;
    c.xy(i0, j0, o.x0, o.y0);
    o.ix = c.ix * di; o.iy = c.iy * di; o.jx = c.jx * dj; o.jy = c.jy * dj;
    o.a.resize(size_t(o.ni) * size_t(o.nj) * size_t(o.ns));
    if (c.a.empty()) return o;
    for (int i = 0; i < o.ni; ++i)
        for (int j = 0; j < o.nj; ++j)
            for (int s = 0; s < o.ns; ++s)
                o.a[(size_t(i) * size_t(o.nj) + size_t(j)) * size_t(o.ns) + size_t(s)] = c.at(i0 + i * di, j0 + j * dj, s0 + s * ds);
    return o;
}

// C44: where a fault cuts an intersection path - at each depth level the fault's plan trace (stick interpolation) is
// intersected with the path segments; returns (distance along the path, depth) points, one polyline per crossing branch
std::vector<std::vector<std::array<double, 2>>> faultCutOnPath(const Fault& f, const Polygon& P, double z0, double z1, int nz) {
    std::vector<std::vector<std::array<double, 2>>> out;
    if (P.pts.size() < 2 || f.sticks.empty()) return out;
    std::vector<double> cum(P.pts.size(), 0.0);
    for (size_t i = 1; i < P.pts.size(); ++i) cum[i] = cum[i - 1] + std::hypot(P.pts[i][0] - P.pts[i - 1][0], P.pts[i][1] - P.pts[i - 1][1]);
    double fz0 = 1e300, fz1 = -1e300;   // depth extent of the sticks
    for (const auto& st : f.sticks) for (const auto& q : st) { fz0 = std::min(fz0, q[2]); fz1 = std::max(fz1, q[2]); }
    z0 = std::max(z0, fz0); z1 = std::min(z1, fz1);
    if (!(z1 > z0)) return out;
    std::vector<std::array<double, 2>> pts;
    for (int k = 0; k <= nz; ++k) {
        const double z = z0 + (z1 - z0) * k / nz;
        const auto tr = f.traceAt(z);
        for (size_t a = 1; a < tr.size(); ++a)
            for (size_t b = 1; b < P.pts.size(); ++b) {
                const double x1 = P.pts[b - 1][0], y1 = P.pts[b - 1][1], x2 = P.pts[b][0], y2 = P.pts[b][1];
                const double x3 = tr[a - 1][0], y3 = tr[a - 1][1], x4 = tr[a][0], y4 = tr[a][1];
                const double den = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
                if (std::fabs(den) < 1e-12) continue;
                const double t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / den;
                const double u = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / den;
                if (t < 0 || t > 1 || u < 0 || u > 1) continue;
                pts.push_back({cum[b - 1] + t * (cum[b] - cum[b - 1]), z});
            }
    }
    // chain the points into branches: next depth level, nearest distance (a fault crossing the path once gives one branch)
    std::sort(pts.begin(), pts.end(), [](const auto& A, const auto& B) { return A[1] < B[1] || (A[1] == B[1] && A[0] < B[0]); });
    // a trace through a path vertex hits both adjacent segments: drop the duplicate
    pts.erase(std::unique(pts.begin(), pts.end(), [&](const auto& A, const auto& B) { return A[1] == B[1] && std::fabs(A[0] - B[0]) < 1e-9 * (1 + cum.back()); }), pts.end());
    const double dz = (z1 - z0) / nz, L = cum.back();
    for (const auto& q : pts) {
        std::vector<std::array<double, 2>>* best = nullptr; double bd = 0.15 * L + 1;
        for (auto& br : out) {
            const auto& e = br.back();
            if (q[1] - e[1] > 1.5 * dz || q[1] <= e[1]) continue;
            const double d = std::fabs(q[0] - e[0]);
            if (d < bd) { bd = d; best = &br; }
        }
        if (best) best->push_back(q); else out.push_back({q});
    }
    return out;
}

}  // namespace geo
}  // namespace ressim
