// C32 geomodelling: structural pillar grid (fault blocks, stair-stepped faults, zones / layering, well tie) and its geometry
#include "ressim/geomodel.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <set>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kDeg = 3.14159265358979323846 / 180.0;
std::string upperS(std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; }

// fill undefined nodes by iterative neighbour averaging (inpainting) - keeps the defined values
void inpaint(Surface& s) {
    if (!s.valid()) return;
    std::vector<char> def(s.z.size());
    size_t nd = 0;
    for (size_t k = 0; k < s.z.size(); ++k) { def[k] = std::isfinite(s.z[k]); nd += def[k]; }
    if (nd == 0 || nd == s.z.size()) return;
    // initial: breadth-first nearest defined value, then smoothing of the filled nodes
    std::vector<int> q; q.reserve(s.z.size());
    std::vector<char> done = def;
    for (size_t k = 0; k < s.z.size(); ++k) if (def[k]) q.push_back(int(k));
    for (size_t h = 0; h < q.size(); ++h) {
        const int k = q[h], i = k % s.nx, j = k / s.nx;
        const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
        for (int d = 0; d < 4; ++d) {
            const int a = i + di[d], b = j + dj[d];
            if (a < 0 || b < 0 || a >= s.nx || b >= s.ny) continue;
            const size_t n = size_t(a + s.nx * b);
            if (done[n]) continue;
            done[n] = 1; s.z[n] = s.z[size_t(k)]; q.push_back(int(n));
        }
    }
    for (int it = 0; it < 60; ++it)
        for (int j = 0; j < s.ny; ++j)
            for (int i = 0; i < s.nx; ++i) {
                const size_t k = size_t(i + s.nx * j);
                if (def[k]) continue;
                double sum = 0; int n = 0;
                if (i > 0) { sum += s.z[k - 1]; ++n; }
                if (i + 1 < s.nx) { sum += s.z[k + 1]; ++n; }
                if (j > 0) { sum += s.z[k - size_t(s.nx)]; ++n; }
                if (j + 1 < s.ny) { sum += s.z[k + size_t(s.nx)]; ++n; }
                if (n) s.z[k] = sum / n;
            }
}
// bilinear with the lattice clamped (no NaN outside the surface extent once inpainted)
double atClamped(const Surface& s, double x, double y) {
    if (!s.valid()) return kNaN;
    const double c = std::cos(s.rotDeg * kDeg), sn = std::sin(s.rotDeg * kDeg);
    double u = (c * (x - s.x0) + sn * (y - s.y0)) / s.dx, v = (-sn * (x - s.x0) + c * (y - s.y0)) / s.dy;
    u = std::max(0.0, std::min(double(s.nx - 1), u)); v = std::max(0.0, std::min(double(s.ny - 1), v));
    const int i = std::min(int(u), s.nx - 2), j = std::min(int(v), s.ny - 2);
    const double fu = u - i, fv = v - j;
    return (1 - fu) * (1 - fv) * s.node(i, j) + fu * (1 - fv) * s.node(i + 1, j) + (1 - fu) * fv * s.node(i, j + 1) + fu * fv * s.node(i + 1, j + 1);
}
bool segX(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
    const double d1x = bx - ax, d1y = by - ay, d2x = dx - cx, d2y = dy - cy;
    const double den = d1x * d2y - d1y * d2x;
    if (std::fabs(den) < 1e-15) return false;
    const double t = ((cx - ax) * d2y - (cy - ay) * d2x) / den, u = ((cx - ax) * d1y - (cy - ay) * d1x) / den;
    return t >= 0 && t <= 1 && u >= 0 && u <= 1;
}
// signed distance to a polyline (+ left, - right), end segments extended
double signedDist(const std::vector<std::array<double, 2>>& tr, double x, double y) {
    double best = 1e300, sg = 1;
    for (size_t k = 0; k + 1 < tr.size(); ++k) {
        const double ax = tr[k][0], ay = tr[k][1], bx = tr[k + 1][0], by = tr[k + 1][1];
        const double ex = bx - ax, ey = by - ay, L2 = ex * ex + ey * ey;
        double t = L2 > 0 ? ((x - ax) * ex + (y - ay) * ey) / L2 : 0.0;
        if (k > 0) t = std::max(t, 0.0);
        if (k + 2 < tr.size()) t = std::min(t, 1.0);
        const double px = ax + t * ex, py = ay + t * ey, d = (x - px) * (x - px) + (y - py) * (y - py);
        if (d < best) { best = d; sg = (ex * (y - ay) - ey * (x - ax)) >= 0 ? 1 : -1; }
    }
    return sg * std::sqrt(best);
}
// side of a point relative to the nearest segment of a polyline (+1 left, -1 right)
int sideOf(const std::vector<std::array<double, 2>>& tr, double x, double y) {
    double best = 1e300; int sg = 1;
    for (size_t k = 0; k + 1 < tr.size(); ++k) {
        const double ax = tr[k][0], ay = tr[k][1], bx = tr[k + 1][0], by = tr[k + 1][1];
        const double ex = bx - ax, ey = by - ay, L2 = ex * ex + ey * ey;
        double t = L2 > 0 ? ((x - ax) * ex + (y - ay) * ey) / L2 : 0.0;
        // the end segments are extended (infinite lines beyond the tips) so the side is defined everywhere
        if (k > 0) t = std::max(t, 0.0);
        if (k + 2 < tr.size()) t = std::min(t, 1.0);
        const double px = ax + t * ex, py = ay + t * ey, d = (x - px) * (x - px) + (y - py) * (y - py);
        if (d < best) { best = d; sg = (ex * (y - ay) - ey * (x - ax)) >= 0 ? 1 : -1; }
    }
    return sg;
}
}  // namespace

void GeoGrid::cornerXY(int i, int j, double& x, double& y) const {
    const double c = std::cos(spec.rotDeg * kDeg), s = std::sin(spec.rotDeg * kDeg);
    const double u = i * spec.dx, v = j * spec.dy;
    x = spec.x0 + c * u - s * v; y = spec.y0 + s * u + c * v;
}
void GeoGrid::columnOf(double x, double y, double& fi, double& fj) const {
    const double c = std::cos(spec.rotDeg * kDeg), s = std::sin(spec.rotDeg * kDeg);
    fi = (c * (x - spec.x0) + s * (y - spec.y0)) / spec.dx;
    fj = (-s * (x - spec.x0) + c * (y - spec.y0)) / spec.dy;
}
void GeoGrid::pillarXY(int i, int j, double z, double& x, double& y) const {
    const int nx = grid.nx > 0 ? grid.nx : spec.nx, ny = grid.ny > 0 ? grid.ny : spec.ny;
    if (grid.coord.size() == size_t(6) * size_t(nx + 1) * size_t(ny + 1) && i >= 0 && j >= 0 && i <= nx && j <= ny) {
        const double* p = &grid.coord[size_t(6) * size_t(i + (nx + 1) * j)];
        const double w = std::fabs(p[5] - p[2]) > 1e-12 ? (z - p[2]) / (p[5] - p[2]) : 0.0;
        x = p[0] + w * (p[3] - p[0]); y = p[1] + w * (p[4] - p[1]);
        return;
    }
    cornerXY(i, j, x, y);
}
double GeoGrid::maxPillarSlope() const {
    double m = 0;
    for (size_t p = 0; p + 5 < grid.coord.size(); p += 6) {
        const double dz = grid.coord[p + 5] - grid.coord[p + 2];
        if (std::fabs(dz) > 1e-12) m = std::max(m, std::hypot(grid.coord[p + 3] - grid.coord[p], grid.coord[p + 4] - grid.coord[p + 1]) / std::fabs(dz));
    }
    return m;
}
int GeoGrid::cellOf(double x, double y, double z) const {
    const int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    double fi, fj; columnOf(x, y, fi, fj);
    // sloped / irregular pillars: invert the bilinear pillar offsets at this depth (fixed point, a few iterations)
    if (grid.coord.size() == size_t(6) * size_t(nx + 1) * size_t(ny + 1) && nx > 0 && ny > 0) {
        for (int it = 0; it < 6; ++it) {
            const int i = std::max(0, std::min(nx - 1, int(std::floor(fi)))), j = std::max(0, std::min(ny - 1, int(std::floor(fj))));
            const double u = std::max(0.0, std::min(1.0, fi - i)), v = std::max(0.0, std::min(1.0, fj - j));
            double sx = 0, sy = 0;
            for (int b = 0; b < 2; ++b) for (int a = 0; a < 2; ++a) {
                double px, py, cx, cy; pillarXY(i + a, j + b, z, px, py); cornerXY(i + a, j + b, cx, cy);
                const double w = (a ? u : 1 - u) * (b ? v : 1 - v);
                sx += w * (px - cx); sy += w * (py - cy);
            }
            double ni, nj; columnOf(x - sx, y - sy, ni, nj);
            const bool conv = std::fabs(ni - fi) < 1e-9 && std::fabs(nj - fj) < 1e-9;
            fi = ni; fj = nj;
            if (conv) break;
        }
    }
    if (fi < 0 || fj < 0 || fi >= nx || fj >= ny || grid.zcorn.empty()) return -1;
    const int i = std::min(int(fi), nx - 1), j = std::min(int(fj), ny - 1);
    const double u = fi - i, v = fj - j;
    auto zc = [&](int k, int f, int a, int b) { return grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j + b) * 2 * nx + 2 * i + a)]; };
    auto bil = [&](int k, int f) { return (1 - u) * (1 - v) * zc(k, f, 0, 0) + u * (1 - v) * zc(k, f, 1, 0) + (1 - u) * v * zc(k, f, 0, 1) + u * v * zc(k, f, 1, 1); };
    for (int k = 0; k < nz; ++k) {
        const double t = bil(k, 0), b = bil(k, 1);
        if (z >= t && z < b) return i + nx * (j + ny * k);
    }
    return -1;
}

void computeGeometry(GeoGrid& g) {
    const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
    const size_t n = size_t(nx) * size_t(ny) * size_t(nz);
    g.center.assign(n, P3{0, 0, 0}); g.bulk.assign(n, 0); g.top.assign(n, 0); g.bot.assign(n, 0);
    const double area = std::fabs(g.spec.dx * g.spec.dy);
    const bool useCoord = g.grid.coord.size() == size_t(6) * size_t(nx + 1) * size_t(ny + 1);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = size_t(i + nx * (j + ny * k));
                double t = 0, b = 0;
                for (int bb = 0; bb < 2; ++bb) for (int a = 0; a < 2; ++a) {
                    t += g.grid.zcorn[size_t(((2 * k) * 2 * ny + 2 * j + bb) * 2 * nx + 2 * i + a)];
                    b += g.grid.zcorn[size_t(((2 * k + 1) * 2 * ny + 2 * j + bb) * 2 * nx + 2 * i + a)];
                }
                t *= 0.25; b *= 0.25;
                g.top[c] = t; g.bot[c] = b;
                // vertical pillars on a parallelogram lattice: volume between two bilinear faces = area x mean corner thickness
                if (useCoord) {   // plan quadrilateral of the pillars at the cell mid-depth (exact for vertical pillars, any lattice; sloped: mid-depth area)
                    const double zc = 0.5 * (t + b);
                    double xs[4], ys[4];
                    g.pillarXY(i, j, zc, xs[0], ys[0]); g.pillarXY(i + 1, j, zc, xs[1], ys[1]); g.pillarXY(i + 1, j + 1, zc, xs[2], ys[2]); g.pillarXY(i, j + 1, zc, xs[3], ys[3]);
                    double ar = 0;
                    for (int q = 0; q < 4; ++q) ar += xs[q] * ys[(q + 1) % 4] - xs[(q + 1) % 4] * ys[q];
                    g.bulk[c] = 0.5 * std::fabs(ar) * std::max(0.0, b - t);
                    g.center[c] = P3{0.25 * (xs[0] + xs[1] + xs[2] + xs[3]), 0.25 * (ys[0] + ys[1] + ys[2] + ys[3]), 0.5 * (t + b)};
                } else {
                    g.bulk[c] = area * std::max(0.0, b - t);
                    double x0, y0, x1, y1;
                    g.cornerXY(i, j, x0, y0); g.cornerXY(i + 1, j + 1, x1, y1);
                    g.center[c] = P3{0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (t + b)};
                }
            }
}

GeoGrid buildGrid(const GridSpec& spec, std::vector<Surface> hz, const std::vector<Fault>& faults, const std::vector<GeoWell>& wells, std::string& err) {
    GeoGrid G; G.spec = spec;
    const int nx = spec.nx, ny = spec.ny;
    if (nx < 1 || ny < 1 || spec.dx <= 0 || spec.dy <= 0) { err = "grid lattice: nx, ny, dx, dy must be positive"; return G; }
    if (hz.size() < 2) { err = "at least two horizons (top and base) are needed"; return G; }
    for (auto& h : hz) { if (!h.valid()) { err = "horizon '" + h.name + "' is empty"; return G; } inpaint(h); }
    std::vector<Zone> zones = spec.zones;
    if (zones.empty()) for (size_t h = 0; h + 1 < hz.size(); ++h) { Zone z; z.name = hz[h].name + "-" + hz[h + 1].name; z.top = int(h); z.base = int(h + 1); z.nLayers = 5; zones.push_back(z); }
    for (const auto& z : zones)
        if (z.top < 0 || z.base < 0 || size_t(z.top) >= hz.size() || size_t(z.base) >= hz.size() || z.nLayers < 1) { err = "zone '" + z.name + "': bad horizon index / layer count"; return G; }
    int nz = 0;
    for (const auto& z : zones) nz += z.nLayers;
    G.spec.zones = zones;
    // ---- mid-depth fault traces, column sides, cut faces
    double zMid = 0; int nzm = 0;
    for (const auto& h : hz) for (double v : h.z) if (std::isfinite(v)) { zMid += v; ++nzm; }
    zMid = nzm ? zMid / nzm : 0;
    const size_t ncol = size_t(nx) * size_t(ny);
    std::vector<std::vector<std::array<double, 2>>> traces(faults.size());
    std::vector<std::vector<int>> side(faults.size(), std::vector<int>(ncol, 1));
    std::vector<int> hwSide(faults.size(), 1);
    auto colC = [&](int i, int j, double& x, double& y) { double x0, y0, x1, y1; G.cornerXY(i, j, x0, y0); G.cornerXY(i + 1, j + 1, x1, y1); x = 0.5 * (x0 + x1); y = 0.5 * (y0 + y1); };
    for (size_t f = 0; f < faults.size(); ++f) {
        traces[f] = faults[f].traceAt(zMid);
        if (traces[f].size() < 2) continue;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { double x, y; colC(i, j, x, y); side[f][size_t(i + nx * j)] = sideOf(traces[f], x, y); }
        // hanging wall: the side the fault moves toward with depth (normal-fault convention); vertical: left
        const auto deep = faults[f].traceAt(zMid + 100), shal = faults[f].traceAt(zMid - 100);
        double mx = 0, my = 0;
        for (size_t k = 0; k < deep.size() && k < shal.size(); ++k) { mx += deep[k][0] - shal[k][0]; my += deep[k][1] - shal[k][1]; }
        if (std::hypot(mx, my) > 1e-9) {
            double cx = 0, cy = 0;
            for (const auto& p : traces[f]) { cx += p[0]; cy += p[1]; }
            cx /= double(traces[f].size()); cy /= double(traces[f].size());
            hwSide[f] = sideOf(traces[f], cx + mx * 10, cy + my * 10);
        }
    }
    // ---- C36 fault-parallel pillars: per pillar dip vector (dx/dz, dy/dz) from the nearest fault planes, Gaussian taper
    std::vector<double> pSlope;   // 2 per pillar node; empty = vertical
    if (spec.pillars == 1 && !faults.empty()) {
        const int npx0 = nx + 1, npy0 = ny + 1;
        pSlope.assign(size_t(2) * size_t(npx0) * size_t(npy0), 0.0);
        std::vector<double> wsum(size_t(npx0) * size_t(npy0), 0.0);
        const double Lt = spec.pillarRange > 0 ? spec.pillarRange : 4.0 * std::max(spec.dx, spec.dy);
        for (size_t f = 0; f < faults.size(); ++f) {
            const auto& tr = traces[f];
            const auto deep = faults[f].traceAt(zMid + 50), shal = faults[f].traceAt(zMid - 50);
            if (tr.size() < 2 || deep.size() != tr.size() || shal.size() != tr.size()) continue;
            std::vector<std::array<double, 2>> sv(tr.size());
            for (size_t k = 0; k < tr.size(); ++k) sv[k] = {(deep[k][0] - shal[k][0]) / 100.0, (deep[k][1] - shal[k][1]) / 100.0};
            for (int jj = 0; jj < npy0; ++jj)
                for (int ii = 0; ii < npx0; ++ii) {
                    double px, py; G.cornerXY(ii, jj, px, py);
                    double best = 1e300, sx = 0, sy = 0;
                    for (size_t k = 0; k + 1 < tr.size(); ++k) {
                        const double ax = tr[k][0], ay = tr[k][1], bx = tr[k + 1][0] - ax, by = tr[k + 1][1] - ay;
                        const double L2 = bx * bx + by * by;
                        const double t = L2 > 0 ? std::max(0.0, std::min(1.0, ((px - ax) * bx + (py - ay) * by) / L2)) : 0.0;
                        const double d = std::hypot(px - ax - t * bx, py - ay - t * by);
                        if (d < best) { best = d; sx = (1 - t) * sv[k][0] + t * sv[k + 1][0]; sy = (1 - t) * sv[k][1] + t * sv[k + 1][1]; }
                    }
                    const double w = std::exp(-(best / Lt) * (best / Lt));
                    const size_t p = size_t(ii + npx0 * jj);
                    pSlope[2 * p] += w * sx; pSlope[2 * p + 1] += w * sy; wsum[p] += w;
                }
        }
        for (size_t p = 0; p < wsum.size(); ++p) if (wsum[p] > 1) { pSlope[2 * p] /= wsum[p]; pSlope[2 * p + 1] /= wsum[p]; }
    }
    auto pillarAt = [&](int ii, int jj, double z, double& x, double& y) {
        G.cornerXY(ii, jj, x, y);
        if (pSlope.empty()) return;
        const size_t p = size_t(ii + (nx + 1) * jj);
        x += pSlope[2 * p] * (z - zMid); y += pSlope[2 * p + 1] * (z - zMid);
    };
    // cut faces: X face (i,j)|(i+1,j), Y face (i,j)|(i,j+1); the centre-to-centre segment crosses the trace
    auto cuts = [&](size_t f, int i0, int j0, int i1, int j1) {
        const auto& tr = traces[f];
        if (tr.size() < 2) return false;
        double ax, ay, bx, by; colC(i0, j0, ax, ay); colC(i1, j1, bx, by);
        for (size_t k = 0; k + 1 < tr.size(); ++k) {
            double cx = tr[k][0], cy = tr[k][1], dx = tr[k + 1][0], dy = tr[k + 1][1];
            // extend the tips by half a cell so a fault ending on the grid edge closes
            if (k == 0) { const double L = std::hypot(dx - cx, dy - cy); if (L > 0) { cx -= (dx - cx) / L * 0.5 * spec.dx; cy -= (dy - cy) / L * 0.5 * spec.dx; } }
            if (k + 2 == tr.size()) { const double L = std::hypot(dx - cx, dy - cy); if (L > 0) { dx += (dx - cx) / L * 0.5 * spec.dx; dy += (dy - cy) / L * 0.5 * spec.dx; } }
            if (segX(ax, ay, bx, by, cx, cy, dx, dy)) return true;
        }
        return false;
    };
    std::vector<std::vector<char>> cutX(faults.size(), std::vector<char>(ncol, 0)), cutY(faults.size(), std::vector<char>(ncol, 0));
    std::vector<char> anyX(ncol, 0), anyY(ncol, 0);
    for (size_t f = 0; f < faults.size(); ++f)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = size_t(i + nx * j);
                if (i + 1 < nx && side[f][c] != side[f][c + 1] && cuts(f, i, j, i + 1, j)) { cutX[f][c] = 1; anyX[c] = 1; }
                if (j + 1 < ny && side[f][c] != side[f][c + size_t(nx)] && cuts(f, i, j, i, j + 1)) { cutY[f][c] = 1; anyY[c] = 1; }
            }
    // ---- fault blocks: connected components across uncut faces
    G.block.assign(ncol, -1);
    int nb = 0;
    for (size_t s = 0; s < ncol; ++s) {
        if (G.block[s] >= 0) continue;
        std::vector<size_t> st{s}; G.block[s] = nb;
        while (!st.empty()) {
            const size_t c = st.back(); st.pop_back();
            const int i = int(c % size_t(nx)), j = int(c / size_t(nx));
            auto go = [&](size_t n) { if (G.block[n] < 0) { G.block[n] = nb; st.push_back(n); } };
            if (i + 1 < nx && !anyX[c]) go(c + 1);
            if (i > 0 && !anyX[c - 1]) go(c - 1);
            if (j + 1 < ny && !anyY[c]) go(c + size_t(nx));
            if (j > 0 && !anyY[c - size_t(nx)]) go(c - size_t(nx));
        }
        ++nb;
    }
    G.nBlocks = nb;
    // throw hint: block entirely on the hanging-wall side of a fault (separated by it)
    std::vector<double> blockThrow(size_t(nb), 0.0);
    for (size_t f = 0; f < faults.size(); ++f) {
        if (faults[f].throwHint <= 0 || traces[f].size() < 2) continue;
        std::vector<int> sMin(size_t(nb), 2), sMax(size_t(nb), -2);
        for (size_t c = 0; c < ncol; ++c) { const int b = G.block[c]; sMin[size_t(b)] = std::min(sMin[size_t(b)], side[f][c]); sMax[size_t(b)] = std::max(sMax[size_t(b)], side[f][c]); }
        for (int b = 0; b < nb; ++b) if (sMin[size_t(b)] == sMax[size_t(b)] && sMin[size_t(b)] == hwSide[f]) blockThrow[size_t(b)] += faults[f].throwHint;
    }
    // ---- per pillar node and block: sample position (pulled into the block when several blocks meet)
    const int npx = nx + 1, npy = ny + 1;
    struct NodeB { int block; double x, y; };   // x, y: sample position, or for split nodes the direction into the block (split = true)
    std::vector<std::vector<NodeB>> nodeB(size_t(npx) * size_t(npy));
    for (int jj = 0; jj < npy; ++jj)
        for (int ii = 0; ii < npx; ++ii) {
            std::vector<std::pair<int, std::array<double, 2>>> adj;
            for (int b = -1; b <= 0; ++b) for (int a = -1; a <= 0; ++a) {
                const int i = ii + a, j = jj + b;
                if (i < 0 || j < 0 || i >= nx || j >= ny) continue;
                double x, y; colC(i, j, x, y);
                adj.push_back({G.block[size_t(i + nx * j)], {x, y}});
            }
            std::set<int> bl; for (auto& a : adj) bl.insert(a.first);
            double px, py; G.cornerXY(ii, jj, px, py);
            auto& out = nodeB[size_t(ii + npx * jj)];
            for (int b : bl) {
                if (bl.size() == 1) { out.push_back({b, px, py}); continue; }
                double cx = 0, cy = 0; int n = 0;
                for (auto& a : adj) if (a.first == b) { cx += a.second[0]; cy += a.second[1]; ++n; }
                cx /= n; cy /= n;
                out.push_back({b, cx - px, cy - py});   // direction into the block (resolved below)
            }
        }
    // ---- horizon depths per node / block, well tie (kriged residuals, simple kriging to zero away from wells)
    const size_t nh = hz.size();
    auto hDepth = [&](size_t h, double x, double y, int b) { return atClamped(hz[h], x, y) + blockThrow[size_t(b)]; };
    // split node (several blocks meet): the horizon is sampled at 0.75 and 1.5 cell sizes into the block and extrapolated
    // linearly back to the pillar, so a horizon that ramps across the fault (smooth / gridded through the fault gap) does
    // not smear the throw
    const double aStep = 0.75 * std::max(spec.dx, spec.dy);
    // side of each fault for every block (0: the block lies on both sides - fault terminates inside it)
    std::vector<std::vector<int>> blockSide(faults.size(), std::vector<int>(size_t(nb), 0));
    for (size_t f = 0; f < faults.size(); ++f) {
        if (traces[f].size() < 2) continue;
        std::vector<int> sMin(size_t(nb), 2), sMax(size_t(nb), -2);
        for (size_t c = 0; c < ncol; ++c) { const int b = G.block[c]; sMin[size_t(b)] = std::min(sMin[size_t(b)], side[f][c]); sMax[size_t(b)] = std::max(sMax[size_t(b)], side[f][c]); }
        for (int b = 0; b < nb; ++b) if (sMin[size_t(b)] == sMax[size_t(b)]) blockSide[f][size_t(b)] = sMin[size_t(b)];
    }
    double hzSpacing = 0;
    for (const auto& h : hz) hzSpacing = std::max(hzSpacing, std::max(h.dx, h.dy));
    auto nodeDepthXY = [&](size_t h, double px, double py, const NodeB& q, bool split) {
        if (!split) return hDepth(h, px, py, q.block);
        const double L = std::hypot(q.x, q.y);
        const double ux = L > 0 ? q.x / L : 0, uy = L > 0 ? q.y / L : 0;
        // the fault traces at this horizon's depth (dipping faults move with depth): the two samples must lie on the
        // block's side of every separating fault, at least one map cell away from it (outside the map's ramp)
        const double zEst = atClamped(hz[h], px, py);
        std::vector<std::pair<size_t, std::vector<std::array<double, 2>>>> trH;
        for (size_t f = 0; f < faults.size(); ++f) if (blockSide[f][size_t(q.block)] != 0) trH.push_back({f, faults[f].traceAt(zEst)});
        auto okAt = [&](double x, double y) {
            for (const auto& t : trH) {
                if (t.second.size() < 2) continue;
                const double d = signedDist(t.second, x, y) * blockSide[t.first][size_t(q.block)];
                if (d < 0.75 * hzSpacing) return false;
            }
            return true;
        };
        double s1 = -1;
        for (int k = 1; k <= 8 && s1 < 0; ++k) { const double s = 0.5 * k * aStep; if (okAt(px + ux * s, py + uy * s)) s1 = s; }
        if (s1 < 0) s1 = aStep;
        const double s2 = s1 + aStep;
        const double z1 = hDepth(h, px + ux * s1, py + uy * s1, q.block), z2 = hDepth(h, px + ux * s2, py + uy * s2, q.block);
        return z1 - (z2 - z1) * s1 / (s2 - s1);   // linear extrapolation back to the pillar
    };
    // sloped pillars: the horizon is met where the pillar line crosses it - fixed point z = h(x(z), y(z))
    auto nodeDepth = [&](size_t h, int ii, int jj, const NodeB& q, bool split) {
        double px, py; G.cornerXY(ii, jj, px, py);
        double z = nodeDepthXY(h, px, py, q, split);
        if (!pSlope.empty())
            for (int it = 0; it < 8; ++it) {
                pillarAt(ii, jj, z, px, py);
                const double zn = nodeDepthXY(h, px, py, q, split);
                const bool conv = std::fabs(zn - z) < 1e-6;
                z = zn;
                if (conv) break;
            }
        return z;
    };
    auto samplePos = [&](int ii, int jj, const NodeB& q, bool split, double& x, double& y) {
        G.cornerXY(ii, jj, x, y);
        if (!split) return;
        const double L = std::hypot(q.x, q.y);
        if (L > 0) { x += q.x / L * 0.25 * aStep; y += q.y / L * 0.25 * aStep; }
    };
    std::vector<std::vector<double>> tieCorr(nh);
    G.tieReport.clear();
    if (spec.tieToWells && !wells.empty()) {
        const double rng = spec.tieRange > 0 ? spec.tieRange : 10.0 * std::max(spec.dx, spec.dy);
        for (size_t h = 0; h < nh; ++h) {
            std::vector<geostat::Point> d;
            for (const auto& w : wells)
                for (const auto& t : w.tops) {
                    if (upperS(t.surface) != upperS(hz[h].name)) continue;
                    const P3 p = w.at(t.md);
                    double fi, fj; G.columnOf(p[0], p[1], fi, fj);
                    const int i = std::max(0, std::min(nx - 1, int(fi))), j = std::max(0, std::min(ny - 1, int(fj)));
                    const double zs = hDepth(h, p[0], p[1], G.block[size_t(i + nx * j)]);
                    d.push_back({p[0], p[1], 0.0, p[2] - zs});
                    char b[200]; std::snprintf(b, sizeof b, "%s @ %s: top %.2f, surface %.2f, residual %+.2f", hz[h].name.c_str(), w.name.c_str(), p[2], zs, p[2] - zs);
                    G.tieReport.push_back(b);
                }
            if (d.empty()) continue;
            std::vector<geostat::Xyz> tg;
            for (int jj = 0; jj < npy; ++jj) for (int ii = 0; ii < npx; ++ii) {
                const auto& nbv = nodeB[size_t(ii + npx * jj)];
                for (const auto& q : nbv) { double x, y; samplePos(ii, jj, q, nbv.size() > 1, x, y); tg.push_back({x, y, 0.0}); }
            }
            geostat::Options o; o.vg.model = geostat::VarioModel::Gaussian; o.vg.nugget = 0; o.vg.sill = 1; o.vg.rangeMax = o.vg.rangeMin = rng; o.vg.rangeVert = 1e9;
            o.simpleKriging = true; o.skMean = 0; o.maxNbr = 12; o.searchRanges = 3.0;
            if (d.size() == 1) { o.simpleKriging = true; }
            const auto K = geostat::krige(tg, d, o);
            tieCorr[h] = K.est;
            for (double& v : tieCorr[h]) if (!std::isfinite(v)) v = 0;
        }
    }
    // node / block depth table: [nodeSlot][h]
    std::vector<std::vector<std::vector<double>>> zN(nodeB.size());
    size_t slot = 0;
    for (size_t nd = 0; nd < nodeB.size(); ++nd) {
        zN[nd].resize(nodeB[nd].size());
        const int ii = int(nd % size_t(npx)), jj = int(nd / size_t(npx));
        for (size_t q = 0; q < nodeB[nd].size(); ++q, ++slot) {
            auto& v = zN[nd][q]; v.resize(nh);
            for (size_t h = 0; h < nh; ++h) {
                v[h] = nodeDepth(h, ii, jj, nodeB[nd][q], nodeB[nd].size() > 1);
                if (!tieCorr[h].empty() && slot < tieCorr[h].size()) v[h] += tieCorr[h][slot];
            }
            for (size_t h = 1; h < nh; ++h) v[h] = std::max(v[h], v[h - 1]);   // younger horizons truncate (no crossing)
        }
    }
    auto zAt = [&](int ii, int jj, int b, size_t h) {
        const size_t nd = size_t(ii + npx * jj);
        for (size_t q = 0; q < nodeB[nd].size(); ++q) if (nodeB[nd][q].block == b) return zN[nd][q][h];
        return zN[nd][0][h];
    };
    // ---- zones -> layers -> ZCORN
    GridInput& gi = G.grid;
    gi = GridInput();
    gi.nx = nx; gi.ny = ny; gi.nz = nz;
    gi.zcorn.assign(size_t(8) * ncol * size_t(nz), 0.0);
    G.layerZone.clear();
    for (size_t z = 0; z < zones.size(); ++z) for (int l = 0; l < zones[z].nLayers; ++l) G.layerZone.push_back(int(z));
    // mean zone thickness (for the parallel layering modes)
    std::vector<double> zMean(zones.size(), 0.0);
    for (size_t z = 0; z < zones.size(); ++z) {
        double s = 0; int n = 0;
        for (int jj = 0; jj < npy; ++jj) for (int ii = 0; ii < npx; ++ii) { s += zN[size_t(ii + npx * jj)][0][size_t(zones[z].base)] - zN[size_t(ii + npx * jj)][0][size_t(zones[z].top)]; ++n; }
        zMean[z] = n ? std::max(1e-6, s / n) : 1.0;
    }
    double zmin = 1e300, zmax = -1e300;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int b = G.block[size_t(i + nx * j)];
            for (int cb = 0; cb < 2; ++cb)
                for (int ca = 0; ca < 2; ++ca) {
                    int k = 0;
                    for (size_t z = 0; z < zones.size(); ++z) {
                        const Zone& Z = zones[z];
                        const double zt = zAt(i + ca, j + cb, b, size_t(Z.top)), zb = std::max(zt, zAt(i + ca, j + cb, b, size_t(Z.base)));
                        const double tl = zMean[z] / Z.nLayers;
                        auto lev = [&](int l) {
                            if (Z.mode == 1) return std::min(zt + l * tl, zb);
                            if (Z.mode == 2) return std::max(zb - (Z.nLayers - l) * tl, zt);
                            return zt + (zb - zt) * double(l) / Z.nLayers;
                        };
                        for (int l = 0; l < Z.nLayers; ++l, ++k) {
                            const double t = lev(l), bo = lev(l + 1);
                            gi.zcorn[size_t(((2 * k) * 2 * ny + 2 * j + cb) * 2 * nx + 2 * i + ca)] = t;
                            gi.zcorn[size_t(((2 * k + 1) * 2 * ny + 2 * j + cb) * 2 * nx + 2 * i + ca)] = bo;
                            zmin = std::min(zmin, t); zmax = std::max(zmax, bo);
                        }
                    }
                }
        }
    // ---- COORD (vertical pillars, or C36 fault-parallel sloped pillar lines)
    gi.coord.resize(size_t(6) * size_t(npx) * size_t(npy));
    for (int jj = 0; jj < npy; ++jj)
        for (int ii = 0; ii < npx; ++ii) {
            double xt, yt, xb, yb; pillarAt(ii, jj, zmin - 1, xt, yt); pillarAt(ii, jj, zmax + 1, xb, yb);
            double* p = &gi.coord[size_t(6) * size_t(ii + npx * jj)];
            p[0] = xt; p[1] = yt; p[2] = zmin - 1; p[3] = xb; p[4] = yb; p[5] = zmax + 1;
        }
    computeGeometry(G);
    gi.actnum.assign(size_t(ncol) * size_t(nz), 1);
    for (size_t c = 0; c < gi.actnum.size(); ++c)
        if (G.bot[c] - G.top[c] <= std::max(spec.minThickness, 1e-6)) gi.actnum[c] = 0;
    // ---- FAULTS (stair-step faces, all layers) + MULTFLT
    for (size_t f = 0; f < faults.size(); ++f) {
        FaultInput F; F.name = faults[f].name.empty() ? "F" + std::to_string(f + 1) : faults[f].name; F.mult = faults[f].mult;
        for (std::string::size_type p; (p = F.name.find(' ')) != std::string::npos;) F.name[p] = '_';
        if (F.name.size() > 8) F.name = F.name.substr(0, 8);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = size_t(i + nx * j);
                if (cutX[f][c]) { FaultFace ff; ff.i1 = ff.i2 = i; ff.j1 = ff.j2 = j; ff.k1 = 0; ff.k2 = nz - 1; ff.face = 'X'; F.faces.push_back(ff); }
                if (cutY[f][c]) { FaultFace ff; ff.i1 = ff.i2 = i; ff.j1 = ff.j2 = j; ff.k1 = 0; ff.k2 = nz - 1; ff.face = 'Y'; F.faces.push_back(ff); }
            }
        // merge runs of X faces along j / Y faces along i
        std::stable_sort(F.faces.begin(), F.faces.end(), [](const FaultFace& a, const FaultFace& b) {
            if (a.face != b.face) return a.face < b.face;
            return a.face == 'X' ? (a.i1 != b.i1 ? a.i1 < b.i1 : a.j1 < b.j1) : (a.j1 != b.j1 ? a.j1 < b.j1 : a.i1 < b.i1);
        });
        std::vector<FaultFace> merged;
        for (const auto& ff : F.faces) {
            if (!merged.empty()) {
                FaultFace& m = merged.back();
                if (ff.face == 'X' && m.face == 'X' && m.i1 == ff.i1 && m.j2 + 1 == ff.j1) { m.j2 = ff.j2; continue; }
                if (ff.face == 'Y' && m.face == 'Y' && m.j1 == ff.j1 && m.i2 + 1 == ff.i1) { m.i2 = ff.i2; continue; }
            }
            merged.push_back(ff);
        }
        F.faces = merged;
        if (!F.faces.empty()) gi.faults.push_back(F);
    }
    return G;
}

}  // namespace geo
}  // namespace ressim
