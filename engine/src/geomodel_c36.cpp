// C36 geomodelling round 2 (geomodel.h): surface editing, multiple-point statistics (Direct Sampling) facies and the
// volumetric uncertainty loop (structural random field + contacts + petrophysics + stochastic realisations, tornado).
#include "ressim/geomodel.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

namespace ressim {
namespace geo {
namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kPi = 3.14159265358979323846;
}  // namespace

// ============================================================================ surface editing
bool surfaceIJ(const Surface& s, double x, double y, double& fi, double& fj) {
    const double c = std::cos(s.rotDeg * kPi / 180), sn = std::sin(s.rotDeg * kPi / 180);
    fi = (c * (x - s.x0) + sn * (y - s.y0)) / s.dx;
    fj = (-sn * (x - s.x0) + c * (y - s.y0)) / s.dy;
    return fi >= 0 && fj >= 0 && fi <= s.nx - 1 && fj <= s.ny - 1;
}

static double localMean(const Surface& s, int i, int j) {
    double sum = 0, w = 0;
    for (int b = -1; b <= 1; ++b) for (int a = -1; a <= 1; ++a) {
        const int ii = i + a, jj = j + b;
        if (ii < 0 || jj < 0 || ii >= s.nx || jj >= s.ny) continue;
        const double v = s.node(ii, jj);
        if (!std::isfinite(v)) continue;
        const double ww = (a == 0 ? 2.0 : 1.0) * (b == 0 ? 2.0 : 1.0);
        sum += ww * v; w += ww;
    }
    return w > 0 ? sum / w : kNaN;
}

int brushSurface(Surface& s, double x, double y, double radius, double amount, int mode) {
    if (!s.valid() || !(radius > 0)) return 0;
    double fi, fj; surfaceIJ(s, x, y, fi, fj);
    const int ri = int(std::ceil(radius / std::fabs(s.dx))) + 1, rj = int(std::ceil(radius / std::fabs(s.dy))) + 1;
    const int i0 = std::max(0, int(std::floor(fi)) - ri), i1 = std::min(s.nx - 1, int(std::ceil(fi)) + ri);
    const int j0 = std::max(0, int(std::floor(fj)) - rj), j1 = std::min(s.ny - 1, int(std::ceil(fj)) + rj);
    if (i0 > i1 || j0 > j1) return 0;
    const Surface src = s;   // smoothing reads the unmodified surface
    int n = 0;
    for (int j = j0; j <= j1; ++j)
        for (int i = i0; i <= i1; ++i) {
            double nx, ny; s.nodeXY(i, j, nx, ny);
            const double r2 = (nx - x) * (nx - x) + (ny - y) * (ny - y);
            if (r2 > radius * radius) continue;
            const double z = src.node(i, j);
            if (!std::isfinite(z)) continue;
            const double w = std::exp(-3.0 * r2 / (radius * radius));
            double zn = z;
            if (mode == 0) zn = z + amount * w;
            else if (mode == 1) { const double m = localMean(src, i, j); if (std::isfinite(m)) zn = z + std::clamp(amount, 0.0, 1.0) * w * (m - z); }
            else zn = z + w * (amount - z);
            if (zn != z) { s.z[size_t(i + s.nx * j)] = zn; ++n; }
        }
    return n;
}

void smoothSurface(Surface& s, int passes) {
    if (!s.valid()) return;
    for (int p = 0; p < passes; ++p) {
        const Surface src = s;
        for (int j = 0; j < s.ny; ++j)
            for (int i = 0; i < s.nx; ++i) {
                if (!std::isfinite(src.node(i, j))) continue;
                const double m = localMean(src, i, j);
                if (std::isfinite(m)) s.z[size_t(i + s.nx * j)] = m;
            }
    }
}

int fillSurfaceHoles(Surface& s) {
    if (!s.valid()) return 0;
    std::vector<size_t> holes;
    double mean = 0; int nd = 0;
    for (size_t k = 0; k < s.z.size(); ++k) { if (!std::isfinite(s.z[k])) holes.push_back(k); else { mean += s.z[k]; ++nd; } }
    if (holes.empty() || nd == 0) return 0;
    mean /= nd;
    // initial guess: nearest-ring average, then Gauss-Seidel on the Laplace equation (defined nodes fixed)
    for (size_t k : holes) s.z[k] = mean;
    for (int it = 0; it < 2000; ++it) {
        double dmax = 0;
        for (size_t k : holes) {
            const int i = int(k % size_t(s.nx)), j = int(k / size_t(s.nx));
            double sum = 0; int n = 0;
            if (i > 0) { sum += s.z[k - 1]; ++n; }
            if (i + 1 < s.nx) { sum += s.z[k + 1]; ++n; }
            if (j > 0) { sum += s.z[k - size_t(s.nx)]; ++n; }
            if (j + 1 < s.ny) { sum += s.z[k + size_t(s.nx)]; ++n; }
            const double v = sum / std::max(1, n);
            dmax = std::max(dmax, std::fabs(v - s.z[k]));
            s.z[k] = v;
        }
        if (dmax < 1e-6) break;
    }
    return int(holes.size());
}

int enforceThickness(std::vector<Surface>& hz, double minThick) {
    int n = 0;
    for (size_t h = 1; h < hz.size(); ++h) {
        Surface& s = hz[h];
        if (!s.valid()) continue;
        for (int j = 0; j < s.ny; ++j)
            for (int i = 0; i < s.nx; ++i) {
                double& z = s.z[size_t(i + s.nx * j)];
                if (!std::isfinite(z)) continue;
                double x, y; s.nodeXY(i, j, x, y);
                const double up = hz[h - 1].at(x, y);
                if (std::isfinite(up) && z < up + minThick) { z = up + minThick; ++n; }
            }
    }
    return n;
}

Surface isochoreSurface(const Surface& ref, double thickness, const std::string& name) {
    Surface s = ref; s.name = name;
    for (double& z : s.z) if (std::isfinite(z)) z += thickness;
    return s;
}

// ============================================================================ training images + Direct Sampling MPS
int TrainingImage::ncat() const { int m = 0; for (int x : v) m = std::max(m, x + 1); return m; }

TrainingImage channelTrainingImage(int nx, int ny, int nz, double width, double wavelength, double amplitude, double fraction, unsigned seed) {
    TrainingImage ti; ti.nx = std::max(1, nx); ti.ny = std::max(1, ny); ti.nz = std::max(1, nz);
    ti.v.assign(size_t(ti.nx) * size_t(ti.ny) * size_t(ti.nz), 0);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    size_t on = 0;
    const size_t target = size_t(std::clamp(fraction, 0.0, 0.95) * double(ti.v.size()));
    for (int c = 0; c < 400 && on < target; ++c) {
        const double y0 = U(rng) * ti.ny, ph = U(rng) * 2 * kPi, A = amplitude * (0.6 + 0.8 * U(rng)), L = std::max(2.0, wavelength * (0.7 + 0.6 * U(rng)));
        const double w = std::max(1.0, width * (0.8 + 0.4 * U(rng)));
        const double A2 = 0.25 * A * U(rng), ph2 = U(rng) * 2 * kPi;   // second harmonic: irregular meanders
        const int thick = ti.nz > 1 ? std::max(1, ti.nz / 3) : 1;
        const int zc = ti.nz > 1 ? int(U(rng) * ti.nz) : 0;
        for (int i = 0; i < ti.nx; ++i) {
            const double xx = i + 0.5;
            const double yc = y0 + A * std::sin(2 * kPi * xx / L + ph) + A2 * std::sin(4 * kPi * xx / L + ph2);
            const double slope = A * 2 * kPi / L * std::cos(2 * kPi * xx / L + ph) + A2 * 4 * kPi / L * std::cos(4 * kPi * xx / L + ph2);
            const double hw = 0.5 * w * std::sqrt(1 + slope * slope);   // constant width normal to the channel axis
            for (int j = std::max(0, int(std::floor(yc - hw))); j <= std::min(ti.ny - 1, int(std::ceil(yc + hw))); ++j) {
                if (std::fabs(j + 0.5 - yc) > hw) continue;
                // lens: thinner toward the banks
                const double f = 1 - std::fabs(j + 0.5 - yc) / hw;
                const int th = std::max(1, int(std::lround(thick * std::sqrt(f))));
                for (int k = std::max(0, zc - th / 2); k <= std::min(ti.nz - 1, zc - th / 2 + th - 1); ++k) {
                    int& v = ti.v[size_t(i) + size_t(ti.nx) * (size_t(j) + size_t(ti.ny) * size_t(k))];
                    if (v == 0) { v = 1; ++on; }
                }
            }
        }
    }
    return ti;
}

TrainingImage readTrainingImage(const std::string& path, std::string& err) {
    TrainingImage ti;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return ti; }
    std::string title; std::getline(f, title);
    std::vector<int> dims;
    { std::istringstream is(title); std::string tok; while (is >> tok) { char* e = nullptr; const long v = std::strtol(tok.c_str(), &e, 10); if (e && *e == 0 && v > 0) dims.push_back(int(v)); } }
    if (dims.size() < 2) { err = "training image: the title line must contain 'nx ny [nz]'"; return ti; }
    ti.nx = dims[0]; ti.ny = dims[1]; ti.nz = dims.size() >= 3 ? dims[2] : 1;
    std::string line; std::getline(f, line);
    const int nvar = std::max(1, std::atoi(line.c_str()));
    for (int v = 0; v < nvar; ++v) std::getline(f, line);
    const size_t n = size_t(ti.nx) * size_t(ti.ny) * size_t(ti.nz);
    ti.v.reserve(n);
    while (ti.v.size() < n && std::getline(f, line)) {
        std::istringstream is(line); double x;
        if (is >> x) ti.v.push_back(std::max(0, int(std::lround(x))));
    }
    if (ti.v.size() != n) { err = "training image: " + std::to_string(ti.v.size()) + " values, expected " + std::to_string(n); ti = TrainingImage(); }
    return ti;
}

bool writeTrainingImage(const std::string& path, const TrainingImage& ti) {
    std::ofstream f(path);
    if (!f) return false;
    f << ti.nx << " " << ti.ny << " " << ti.nz << "\n1\nfacies\n";
    for (int v : ti.v) f << v << "\n";
    return bool(f);
}

double indicatorContinuity(const std::vector<int>& v, int nx, int ny, int nz, int di, int dj, int cat) {
    long num = 0, den = 0;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int i2 = i + di, j2 = j + dj;
                if (i2 < 0 || j2 < 0 || i2 >= nx || j2 >= ny) continue;
                const size_t a = size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)), b = size_t(i2) + size_t(nx) * (size_t(j2) + size_t(ny) * size_t(k));
                if (v[a] < 0 || v[b] < 0) continue;
                if (v[a] == cat) { ++den; if (v[b] == cat) ++num; }
            }
    return den ? double(num) / double(den) : 0.0;
}

FaciesResult modelFaciesMps(const GeoGrid& g, const std::vector<Blocked>& data, const FaciesSpec& spec, const TrainingImage& ti, const MpsSpec& m) {
    FaciesResult R;
    const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
    const size_t n = size_t(nx) * size_t(ny) * size_t(nz);
    if (n == 0 || ti.v.empty()) return R;
    const int K = std::max(1, ti.ncat());
    // codes: spec / blocked data / 0..K-1
    std::vector<int> codes = spec.codes;
    if (codes.empty()) {
        std::set<int> cs;
        for (const auto& b : data) { auto it = b.v.find(spec.curve); if (it != b.v.end() && std::isfinite(it->second)) cs.insert(int(std::lround(it->second))); }
        codes.assign(cs.begin(), cs.end());
    }
    while (int(codes.size()) < K) { int c = codes.empty() ? 0 : codes.back() + 1; while (std::find(codes.begin(), codes.end(), c) != codes.end()) ++c; codes.push_back(c); }
    auto catOf = [&](int code) { for (size_t q = 0; q < codes.size(); ++q) if (codes[q] == code) return int(q); return -1; };
    std::vector<int> sim(n, -1);
    for (size_t c = 0; c < n; ++c) if (!g.grid.actnum.empty() && !g.grid.actnum[c]) sim[c] = -2;
    for (const auto& b : data) {
        if (b.cell < 0 || size_t(b.cell) >= n || sim[size_t(b.cell)] == -2) continue;
        auto it = b.v.find(spec.curve);
        if (it == b.v.end() || !std::isfinite(it->second)) continue;
        const int ct = catOf(int(std::lround(it->second)));
        if (ct >= 0 && ct < K && sim[size_t(b.cell)] < 0) { sim[size_t(b.cell)] = ct; ++R.conditioned; }
    }
    const std::vector<int> hard = sim;
    // grid index lag -> TI lag: cell axes in map coordinates, TI +x along the azimuth
    const double rot = g.spec.rotDeg * kPi / 180, az = m.azimuthDeg * kPi / 180;
    const double eix = std::cos(rot) * g.spec.dx, eiy = std::sin(rot) * g.spec.dx, ejx = -std::sin(rot) * g.spec.dy, ejy = std::cos(rot) * g.spec.dy;
    const double ax = std::sin(az), ay = std::cos(az), bx = -ay, by = ax;
    const double cell = m.tiCell > 0 ? m.tiCell : g.spec.dx;
    struct Lag { int di, dj, dk; int tx, ty, tz; double d2; };
    std::vector<Lag> tmpl;
    const int rk = ti.nz > 1 || nz > 1 ? std::max(0, m.radiusK) : 0;
    const int ri = int(std::ceil(m.radius * cell / std::min(std::fabs(g.spec.dx), std::fabs(g.spec.dy)))) + 1;
    for (int dk = -rk; dk <= rk; ++dk)
        for (int dj = -ri; dj <= ri; ++dj)
            for (int di = -ri; di <= ri; ++di) {
                if (!di && !dj && !dk) continue;
                const double px = di * eix + dj * ejx, py = di * eiy + dj * ejy;
                const double tx = (px * ax + py * ay) / cell, ty = (px * bx + py * by) / cell;
                const double pd2 = tx * tx + ty * ty;
                if (pd2 > m.radius * m.radius) continue;
                const double vz = rk > 0 ? double(dk) * m.radius / double(rk) : 0.0;
                tmpl.push_back({di, dj, dk, int(std::lround(tx)), int(std::lround(ty)), ti.nz > 1 ? int(std::lround(dk / std::max(1e-9, m.tiLayer))) : 0, pd2 + vz * vz});
            }
    std::sort(tmpl.begin(), tmpl.end(), [](const Lag& a, const Lag& b) { return a.d2 < b.d2; });
    // marginal of the TI (no-neighbour nodes)
    std::vector<double> marg(size_t(K), 0.0);
    for (int x : ti.v) if (x >= 0 && x < K) marg[size_t(x)] += 1;
    for (double& p : marg) p /= double(ti.v.size());
    std::mt19937 rng(m.seed);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<size_t> path;
    for (size_t c = 0; c < n; ++c) if (sim[c] == -1) path.push_back(c);
    std::shuffle(path.begin(), path.end(), rng);
    const size_t tiN = ti.v.size();
    const long nScan = std::max<long>(1, std::min<long>(m.maxScanAbs, long(m.maxScan * double(tiN))));
    std::vector<std::array<int, 4>> nb;   // tx, ty, tz, category
    nb.reserve(size_t(m.maxNbr));
    std::vector<int> bestK(static_cast<size_t>(K));
    std::vector<double> tgt = marg;
    if (int(m.target.size()) == K) { double t = 0; for (double x : m.target) t += std::max(0.0, x); if (t > 0) for (int q = 0; q < K; ++q) tgt[size_t(q)] = std::max(0.0, m.target[size_t(q)]) / t; }
    std::vector<long> cnt(size_t(K), 0); long nDone = 0;
    for (size_t c = 0; c < n; ++c) if (sim[c] >= 0) { ++cnt[size_t(sim[c])]; ++nDone; }
    for (size_t c : path) {
        const int i = int(c % size_t(nx)), j = int((c / size_t(nx)) % size_t(ny)), k = int(c / (size_t(nx) * size_t(ny)));
        nb.clear();
        for (const Lag& L : tmpl) {
            const int ii = i + L.di, jj = j + L.dj, kk = k + L.dk;
            if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz) continue;
            const int v = sim[size_t(ii) + size_t(nx) * (size_t(jj) + size_t(ny) * size_t(kk))];
            if (v < 0) continue;
            nb.push_back({L.tx, L.ty, L.tz, v});
            if (int(nb.size()) >= m.maxNbr) break;
        }
        int val = -1;
        if (nb.empty()) {
            double u = U(rng), acc = 0;
            for (int q = 0; q < K; ++q) { acc += marg[size_t(q)]; if (u <= acc) { val = q; break; } }
            if (val < 0) val = K - 1;
        } else {
            const int nn = int(nb.size());
            const int accept = int(std::floor(m.threshold * nn + 1e-9));
            // best candidate per centre category (servosystem chooses among them)
            std::fill(bestK.begin(), bestK.end(), nn + 1);
            int best = nn + 1;
            for (long s = 0; s < nScan; ++s) {
                const size_t p = size_t(rng() % tiN);
                const int cat = std::clamp(ti.v[p], 0, K - 1);
                const int ti_i = int(p % size_t(ti.nx)), ti_j = int((p / size_t(ti.nx)) % size_t(ti.ny)), ti_k = int(p / (size_t(ti.nx) * size_t(ti.ny)));
                int mis = 0;
                const int lim = bestK[size_t(cat)];
                for (const auto& q : nb) {
                    const int a = ti_i + q[0], b = ti_j + q[1], cc = ti_k + q[2];
                    if (a < 0 || b < 0 || cc < 0 || a >= ti.nx || b >= ti.ny || cc >= ti.nz || ti.at(a, b, cc) != q[3]) { if (++mis >= lim) break; }
                }
                if (mis < lim) { bestK[size_t(cat)] = mis; best = std::min(best, mis); if (best <= accept) break; }
            }
            double bs = 1e300;
            for (int q = 0; q < K; ++q) {
                if (bestK[size_t(q)] > nn) continue;
                const double pc = nDone > 0 ? double(cnt[size_t(q)]) / double(nDone) : tgt[size_t(q)];
                const double sc = bestK[size_t(q)] - m.servo * nn * (tgt[size_t(q)] - pc);
                if (sc < bs) { bs = sc; val = q; }
            }
            if (val < 0) val = 0;
        }
        sim[c] = std::clamp(val, 0, K - 1);
        ++cnt[size_t(sim[c])]; ++nDone;
    }
    R.codes = codes;
    R.facies.assign(n, codes[0]);
    R.proportions.assign(codes.size(), 0.0);
    long act = 0;
    for (size_t c = 0; c < n; ++c) {
        if (sim[c] < 0) continue;
        R.facies[c] = codes[size_t(sim[c])];
        R.proportions[size_t(sim[c])] += 1; ++act;
        if (hard[c] >= 0 && hard[c] != sim[c]) ++R.mismatches;
    }
    if (act) for (double& p : R.proportions) p /= double(act);
    return R;
}

// ============================================================================ uncertainty loop
double Tri::sample(double u) const {
    if (fixed()) return mode;
    const double a = lo, b = hi, c = std::clamp(mode, lo, hi);
    const double F = (c - a) / (b - a);
    u = std::clamp(u, 0.0, 1.0);
    return u < F ? a + std::sqrt(u * (b - a) * (c - a)) : b - std::sqrt((1 - u) * (b - a) * (b - c));
}

std::vector<Surface> perturbHorizons(const std::vector<Surface>& hz, double sd, double range, unsigned seed, const std::vector<std::array<double, 2>>& zeroAt, double shift) {
    std::vector<Surface> out = hz;
    if ((sd <= 0 && shift == 0) || hz.empty()) return out;
    const int M = 256;
    std::mt19937 rng(seed);
    std::normal_distribution<double> N(0, 1);
    std::uniform_real_distribution<double> U(0, 2 * kPi);
    const double kSd = std::sqrt(6.0) / std::max(1e-9, range);   // Gaussian covariance exp(-3 h^2 / r^2)
    std::vector<double> kx(M), ky(M), ph(M);
    for (int q = 0; q < M; ++q) { kx[size_t(q)] = kSd * N(rng); ky[size_t(q)] = kSd * N(rng); ph[size_t(q)] = U(rng); }
    const double xr = hz[0].x0, yr = hz[0].y0, amp = sd * std::sqrt(2.0 / M);
    auto field = [&](double x, double y) {
        double s = 0;
        for (int q = 0; q < M; ++q) s += std::cos(kx[size_t(q)] * (x - xr) + ky[size_t(q)] * (y - yr) + ph[size_t(q)]);
        return amp * s + shift;
    };
    // conditioning points (deduplicated)
    std::vector<geostat::Point> pts;
    for (const auto& p : zeroAt) {
        bool dup = false;
        for (const auto& q : pts) if (std::hypot(q.x - p[0], q.y - p[1]) < 1e-3 * range) { dup = true; break; }
        if (!dup) { geostat::Point g; g.x = p[0]; g.y = p[1]; g.z = 0; g.v = field(p[0], p[1]); pts.push_back(g); }
    }
    geostat::Options o;
    o.vg.model = geostat::VarioModel::Gaussian; o.vg.nugget = 0; o.vg.sill = 1; o.vg.rangeMax = o.vg.rangeMin = range; o.vg.rangeVert = 1e12;
    o.simpleKriging = true; o.skMean = 0; o.maxNbr = std::min<int>(16, int(pts.size())); o.searchRanges = 4.0;
    for (auto& s : out) {
        if (!s.valid()) continue;
        std::vector<geostat::Xyz> tg; tg.reserve(s.z.size());
        std::vector<double> f(s.z.size());
        for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) { double x, y; s.nodeXY(i, j, x, y); tg.push_back({x, y, 0.0}); f[size_t(i + s.nx * j)] = field(x, y); }
        if (!pts.empty()) {
            const auto K = geostat::krige(tg, pts, o);
            for (size_t q = 0; q < f.size() && q < K.est.size(); ++q) if (std::isfinite(K.est[q])) f[q] -= K.est[q];
        }
        for (size_t q = 0; q < s.z.size(); ++q) if (std::isfinite(s.z[q])) s.z[q] += f[q];
    }
    return out;
}

namespace {
std::vector<std::array<double, 2>> topPoints(const std::vector<GeoWell>& wells) {
    std::vector<std::array<double, 2>> p;
    for (const auto& w : wells) for (const auto& t : w.tops) { const P3 q = w.at(t.md); p.push_back({q[0], q[1]}); }
    return p;
}
bool needGrid(const UncInputs& in, const UncSample& s) { (void)in; return s.field || s.zShift != 0; }
// one realisation; cached grid / blocked data used when the structure is not perturbed
UncModel realiseWith(const UncInputs& in, const UncSample& s, const GeoGrid* cg, const std::vector<Blocked>* cb, std::string& err) {
    UncModel M;
    std::vector<Blocked> Bown;
    const std::vector<Blocked>* B = cb;
    if (needGrid(in, s) || !cg) {
        M.hz = (s.field || s.zShift != 0) ? perturbHorizons(in.hz, s.field ? s.zSd : 0.0, s.zRange, s.fieldSeed, topPoints(in.wells), s.zShift) : in.hz;
        M.g = buildGrid(in.gs, M.hz, in.faults, in.wells, err);
        if (!err.empty()) return M;
        Bown = blockWells(M.g, in.wells); B = &Bown;
    } else { M.g = *cg; M.hz = in.hz; }
    const GeoGrid& g = M.g;
    const bool fac = in.facies;
    if (fac) {
        FaciesSpec fs = in.fs; fs.seed = s.faciesSeed;
        if (fs.method == 3) {
            if (!in.ti || in.ti->v.empty()) { err = "MPS facies: no training image"; return M; }
            MpsSpec m = in.mps; m.seed = s.faciesSeed;
            M.fr = modelFaciesMps(g, *B, fs, *in.ti, m);
        } else M.fr = modelFacies(g, *B, fs);
    }
    PetroSpec ps = in.ps; ps.seed = s.petroSeed;
    M.poro = modelProperty(g, *B, fac ? &M.fr.facies : nullptr, ps).v;
    for (double& p : M.poro) p = std::clamp(p * s.poroMult, 0.0, 1.0);
    M.ntg = in.ntgConst.size() == M.poro.size() ? in.ntgConst : std::vector<double>(M.poro.size(), 1.0);
    for (double& v : M.ntg) v = std::clamp(v * s.ntgMult, 0.0, 1.0);
    if (fac && in.ntgFromFacies) for (size_t c = 0; c < M.ntg.size() && c < M.fr.facies.size(); ++c) if (M.fr.facies[c] == 0) M.ntg[c] = 0.0;
    VolSpec vs = in.vs;
    vs.contacts.owc = s.owc; vs.contacts.goc = s.goc;
    if (vs.satHeight) vs.sh.swirr = s.sw; else vs.swConst = s.sw;
    vs.bo = s.bo;
    M.vol = volumetrics(g, M.poro, M.ntg, fac ? &M.fr.facies : nullptr, vs);
    return M;
}
double pick(const Tri& t, double base, double u) {
    if (t.fixed()) return (t.mode == 0 && t.lo == 0) ? base : t.mode;
    return t.sample(u);
}
double modeOf(const Tri& t, double base) { return t.fixed() ? pick(t, base, 0.5) : std::clamp(t.mode, t.lo, t.hi); }
}  // namespace

UncModel realise(const UncInputs& in, const UncSample& s, std::string& err) {
    return realiseWith(in, s, nullptr, nullptr, err);
}

UncResult uncertaintyLoop(const UncInputs& in, const UncSpec& u, const std::function<bool(double)>& progress) {
    UncResult R;
    std::string err;
    const bool structural = u.zSd > 0;
    // shared grid + blocking when the structure is fixed
    GeoGrid G0; std::vector<Blocked> B0;
    G0 = buildGrid(in.gs, in.hz, in.faults, in.wells, err);
    if (!err.empty()) { R.err = err; return R; }
    B0 = blockWells(G0, in.wells);
    const VolSpec& vs = in.vs;
    const double swBase = vs.satHeight ? vs.sh.swirr : vs.swConst;
    auto baseSample = [&]() {
        UncSample s;
        s.owc = modeOf(u.owc, vs.contacts.owc); s.goc = modeOf(u.goc, vs.contacts.goc);
        s.poroMult = modeOf(u.poroMult, 1.0); s.ntgMult = modeOf(u.ntgMult, 1.0);
        s.sw = modeOf(u.sw, swBase); s.bo = modeOf(u.bo, vs.bo);
        s.zSd = u.zSd; s.zRange = u.zRange;
        if (s.poroMult == 0) s.poroMult = 1;
        if (s.ntgMult == 0) s.ntgMult = 1;
        s.faciesSeed = in.fs.seed; s.petroSeed = in.ps.seed;
        return s;
    };
    // ---- realisations (Latin hypercube or plain Monte Carlo, reproducible per seed)
    const int nr = std::max(1, u.nReal);
    R.real.resize(size_t(nr));
    std::vector<std::string> errs(static_cast<size_t>(nr));
    int done = 0; bool cancel = false;
    // uniform scores per input (6 Tri variables): plain Monte Carlo or Latin hypercube (one sample per stratum, independent permutations)
    std::vector<std::array<double, 6>> us(static_cast<size_t>(nr));
    {
        std::mt19937 rng(u.seed * 2654435761u + 17u);
        std::uniform_real_distribution<double> Ud(0, 1);
        for (int v = 0; v < 6; ++v) {
            std::vector<int> perm(static_cast<size_t>(nr));
            for (int r = 0; r < nr; ++r) perm[size_t(r)] = r;
            std::shuffle(perm.begin(), perm.end(), rng);
            for (int r = 0; r < nr; ++r) us[size_t(r)][size_t(v)] = u.latinHypercube ? (perm[size_t(r)] + Ud(rng)) / double(nr) : Ud(rng);
        }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int r = 0; r < nr; ++r) {
        if (cancel) continue;
        const auto& w = us[size_t(r)];
        UncSample s;
        s.owc = pick(u.owc, vs.contacts.owc, w[0]); s.goc = pick(u.goc, vs.contacts.goc, w[1]);
        s.poroMult = pick(u.poroMult, 1.0, w[2]); s.ntgMult = pick(u.ntgMult, 1.0, w[3]);
        s.sw = pick(u.sw, swBase, w[4]); s.bo = pick(u.bo, vs.bo, w[5]);
        if (s.poroMult == 0) s.poroMult = 1;
        if (s.ntgMult == 0) s.ntgMult = 1;
        s.field = structural; s.fieldSeed = u.seed * 7919u + unsigned(r); s.zSd = u.zSd; s.zRange = u.zRange;
        s.faciesSeed = u.stochastic ? in.fs.seed + 1000u + unsigned(r) : in.fs.seed;
        s.petroSeed = u.stochastic ? in.ps.seed + 5000u + unsigned(r) : in.ps.seed;
        std::string e;
        const UncModel M = realiseWith(in, s, &G0, &B0, e);
        UncReal& out = R.real[size_t(r)];
        out.s = s; out.stoiip = M.vol.total.stoiip; out.giip = M.vol.total.giip; out.pore = M.vol.total.pore; out.bulk = M.vol.total.bulk;
        errs[size_t(r)] = e;
#pragma omp critical(unc_progress)
        {
            ++done;
            if (progress && !progress(double(done) / double(nr + 2 * 7 + u.tornadoStochastic))) cancel = true;
        }
    }
    for (const auto& e : errs) if (!e.empty()) { R.err = e; return R; }
    if (cancel) { R.err = "cancelled"; return R; }
    std::vector<double> so, sg;
    for (const auto& x : R.real) { so.push_back(x.stoiip); sg.push_back(x.giip); }
    R.stoiip = percentiles(so); R.giip = percentiles(sg);
    auto closest = [&](double v) { int b = -1; double d = 1e300; for (size_t q = 0; q < R.real.size(); ++q) if (std::fabs(R.real[q].stoiip - v) < d) { d = std::fabs(R.real[q].stoiip - v); b = int(q); } return b; };
    R.iP90 = closest(R.stoiip.p90); R.iP50 = closest(R.stoiip.p50); R.iP10 = closest(R.stoiip.p10);
    // ---- tornado (one at a time, inputs at their P10 / P90 quantiles, everything else at the mode / base seed)
    const UncSample b0 = baseSample();
    std::string e;
    R.base = realiseWith(in, b0, &G0, &B0, e).vol.total.stoiip;
    struct Job { std::string name; UncSample lo, hi; };
    std::vector<Job> jobs;
    auto addTri = [&](const char* name, const Tri& t, double UncSample::*fld) {
        if (t.fixed()) return;
        Job j{name, b0, b0}; j.lo.*fld = t.sample(0.1); j.hi.*fld = t.sample(0.9); jobs.push_back(j);
    };
    addTri("OWC depth", u.owc, &UncSample::owc);
    addTri("GOC depth", u.goc, &UncSample::goc);
    addTri("Porosity multiplier", u.poroMult, &UncSample::poroMult);
    addTri("NTG multiplier", u.ntgMult, &UncSample::ntgMult);
    addTri(vs.satHeight ? "Swirr" : "Sw", u.sw, &UncSample::sw);
    addTri("Bo", u.bo, &UncSample::bo);
    if (structural) { Job j{"Structure (depth shift +-1.28 sd, tied at wells)", b0, b0}; j.lo.zShift = -1.2816 * u.zSd; j.hi.zShift = 1.2816 * u.zSd; jobs.push_back(j); }
    R.tornado.resize(jobs.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int q = 0; q < int(jobs.size()); ++q) {
        std::string e1, e2;
        const double a = realiseWith(in, jobs[size_t(q)].lo, &G0, &B0, e1).vol.total.stoiip;
        const double c = realiseWith(in, jobs[size_t(q)].hi, &G0, &B0, e2).vol.total.stoiip;
        R.tornado[size_t(q)] = TornadoBar{jobs[size_t(q)].name, std::min(a, c), std::max(a, c)};
    }
    if (u.stochastic && u.tornadoStochastic > 1) {
        std::vector<double> v(static_cast<size_t>(u.tornadoStochastic));
#pragma omp parallel for schedule(dynamic, 1)
        for (int q = 0; q < u.tornadoStochastic; ++q) {
            UncSample s = b0; s.faciesSeed = in.fs.seed + 90000u + unsigned(q); s.petroSeed = in.ps.seed + 95000u + unsigned(q);
            std::string e3;
            v[size_t(q)] = realiseWith(in, s, &G0, &B0, e3).vol.total.stoiip;
        }
        const Percentiles P = percentiles(v);
        R.tornado.push_back(TornadoBar{"Facies / porosity realisations (P90-P10)", P.p90, P.p10});
    }
    std::sort(R.tornado.begin(), R.tornado.end(), [](const TornadoBar& a, const TornadoBar& b) { return a.high - a.low > b.high - b.low; });
    if (progress) progress(1.0);
    // ---- report
    char b[300];
    const char* su = vs.field ? "stb" : "sm3";
    std::snprintf(b, sizeof b, "Uncertainty loop: %d realisations (%s)%s\n", nr, u.latinHypercube ? "Latin hypercube" : "Monte Carlo", structural ? ", structural random field" : ""); R.report += b;
    std::snprintf(b, sizeof b, "  STOIIP P90 %.5g  P50 %.5g  P10 %.5g  mean %.5g %s (base case %.5g)\n", R.stoiip.p90, R.stoiip.p50, R.stoiip.p10, R.stoiip.mean, su, R.base); R.report += b;
    if (R.giip.p10 > 0) { std::snprintf(b, sizeof b, "  GIIP   P90 %.5g  P50 %.5g  P10 %.5g\n", R.giip.p90, R.giip.p50, R.giip.p10); R.report += b; }
    for (const auto& t : R.tornado) { std::snprintf(b, sizeof b, "  tornado %-46s %.5g .. %.5g (swing %.3g %%)\n", t.name.c_str(), t.low, t.high, R.base > 0 ? 100 * (t.high - t.low) / R.base : 0.0); R.report += b; }
    return R;
}

}  // namespace geo
}  // namespace ressim
