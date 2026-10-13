// Grid processing: Cartesian (DX/DY/DZ/TOPS) and corner-point (COORD/ZCORN) geometry,
// active cells, pore volumes, transmissibilities (incl. fault throw NNCs, MULTFLT, explicit NNC).
#include "ressim/grid.h"
#include "ressim/tools.h"
#include <stdexcept>
#include <string>
#include <cfloat>
#include <cmath>
#include <algorithm>

namespace ressim {

namespace {
struct V3 { double x, y, z; };
inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator*(V3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline V3 toV(const Vec3& v) { return {v.x, v.y, v.z}; }

// face corner cycles (consistent orientation, see corner numbering in grid.h)
const int kFaces[6][4] = {{0, 2, 6, 4}, {1, 5, 7, 3},   // -X, +X
                          {0, 4, 5, 1}, {2, 3, 7, 6},   // -Y, +Y
                          {0, 1, 3, 2}, {4, 6, 7, 5}};  // -Z (top), +Z (bottom)

double hexVolume(const std::array<V3, 8>& c) {
    V3 ctr{0, 0, 0};
    for (auto& p : c) ctr = ctr + p;
    ctr = ctr * 0.125;
    double v = 0;
    for (auto& f : kFaces) {
        V3 fc = (c[f[0]] + c[f[1]] + c[f[2]] + c[f[3]]) * 0.25;
        for (int e = 0; e < 4; ++e) {
            const V3 a = c[f[e]] - ctr, b = c[f[(e + 1) % 4]] - ctr, d = fc - ctr;
            v += dot(a, cross(b, d)) / 6.0;
        }
    }
    return std::fabs(v);
}
V3 faceCenter(const std::array<V3, 8>& c, int f) {
    return (c[kFaces[f][0]] + c[kFaces[f][1]] + c[kFaces[f][2]] + c[kFaces[f][3]]) * 0.25;
}
V3 faceAreaVec(const std::array<V3, 8>& c, int f) {
    const V3 d1 = c[kFaces[f][2]] - c[kFaces[f][0]], d2 = c[kFaces[f][3]] - c[kFaces[f][1]];
    return cross(d1, d2) * 0.5;
}

std::vector<double> fillOr(const std::vector<double>& v, size_t n, double def, const char* name, bool required) {
    if (v.empty()) {
        if (required) throw std::runtime_error(std::string("GRID: missing ") + name);
        return std::vector<double>(n, def);
    }
    if (v.size() != n) throw std::runtime_error(std::string("GRID: ") + name + " has " + std::to_string(v.size()) +
                                                " values, expected " + std::to_string(n));
    return v;
}

// Overlap of two vertical face strips spanning the same pair of pillars.
// Face L: depths (top,bottom) on pillar A = (ta1, ba1), pillar B = (tb1, bb1); face R likewise.
// Along the face, s in [0,1] interpolates linearly between the pillars. Returns the overlap
// "area" in the (s, z) plane (times width -> true area), and the centroid (s, z).
struct Overlap { double area = 0, s = 0.5, z = 0; };
Overlap stripOverlap(double ta1, double ba1, double tb1, double bb1, double ta2, double ba2, double tb2, double bb2) {
    auto lin = [](double a, double b, double s) { return a + (b - a) * s; };
    auto h = [&](double s) {
        const double top = std::max(lin(ta1, tb1, s), lin(ta2, tb2, s));
        const double bot = std::min(lin(ba1, bb1, s), lin(ba2, bb2, s));
        return bot - top;
    };
    // breakpoints where any pair of the four lines cross
    double pts[16]; int np = 0;
    pts[np++] = 0; pts[np++] = 1;
    const double L[4][2] = {{ta1, tb1}, {ta2, tb2}, {ba1, bb1}, {ba2, bb2}};
    for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b) {
            const double da = L[a][0] - L[b][0], db = L[a][1] - L[b][1];
            if ((da > 0) != (db > 0) && std::fabs(da - db) > 1e-14) {
                const double s = da / (da - db);
                if (s > 0 && s < 1) pts[np++] = s;
            }
        }
    std::sort(pts, pts + np);
    Overlap o; double sm = 0, zm = 0;
    for (int i = 0; i + 1 < np; ++i) {
        const double s0 = pts[i], s1 = pts[i + 1];
        if (s1 - s0 < 1e-14) continue;
        // Simpson on each piece (h is linear inside; s*h and z*h quadratic -> exact)
        const double sc = 0.5 * (s0 + s1);
        const double h0 = std::max(0.0, h(s0)), h1 = std::max(0.0, h(s1)), hc = std::max(0.0, h(sc));
        auto zmid = [&](double s) {
            return 0.5 * (std::max(lin(ta1, tb1, s), lin(ta2, tb2, s)) + std::min(lin(ba1, bb1, s), lin(ba2, bb2, s)));
        };
        const double w = (s1 - s0) / 6.0;
        o.area += w * (h0 + 4 * hc + h1);
        sm += w * (s0 * h0 + 4 * sc * hc + s1 * h1);
        zm += w * (zmid(s0) * h0 + 4 * zmid(sc) * hc + zmid(s1) * h1);
    }
    if (o.area > 0) { o.s = sm / o.area; o.z = zm / o.area; }
    return o;
}
}  // namespace

double Grid::faceArea(int g, char dir, bool minus) const {
    std::array<V3, 8> c;
    for (int v = 0; v < 8; ++v) c[v] = toV(corners[g][v]);
    const int f = (dir == 'X' ? 0 : dir == 'Y' ? 2 : 4) + (minus ? 0 : 1);
    const V3 a = faceAreaVec(c, f);
    return std::sqrt(dot(a, a));
}

void Grid::build(const Case& c) {
    if (c.grid.newtran && !c.grid.cornerPoint() && !c.grid.dx.empty()) {
        // NEWTRAN on a block-centred grid: build the equivalent corner-point geometry (each cell flat at its own
        // depth, as OPM Flow does) so displaced neighbour columns connect through their actual face overlaps
        Case cp = c;
        convertToCornerPoint(cp.grid);
        build(cp);
        return;
    }
    const GridInput& g = c.grid;
    nx = g.nx; ny = g.ny; nz = g.nz;
    if (nx < 1 || ny < 1 || nz < 1) throw std::runtime_error("GRID: dimensions must be >= 1");
    const size_t n = size_t(nx) * ny * nz;
    const size_t nl = size_t(nx) * ny;
    const Units U = Units::of(c.units);
    cornerPoint = g.cornerPoint();

    auto poro = fillOr(g.poro, n, 0, "PORO", true);
    auto kx = fillOr(g.permx, n, 0, "PERMX", true);
    auto ky = g.permy.empty() ? kx : fillOr(g.permy, n, 0, "PERMY", false);
    auto kz = g.permz.empty() ? kx : fillOr(g.permz, n, 0, "PERMZ", false);
    auto ntg = fillOr(g.ntg, n, 1.0, "NTG", false);
    auto mx = fillOr(g.multx, n, 1.0, "MULTX", false);
    auto my = fillOr(g.multy, n, 1.0, "MULTY", false);
    auto mz = fillOr(g.multz, n, 1.0, "MULTZ", false);

    // ---- dual porosity: effective fracture permeability (ECLIPSE multiplies it by the fracture porosity unless NODPPM)
    const bool dual = g.dual > 0;
    const int nzh = dual ? nz / 2 : nz;
    nzMatrix = dual ? nzh : 0; dualMode = g.dual; numMF = 0;
    permScale.clear();
    if (dual) {
        if (nz % 2 != 0) throw std::runtime_error("GRID: dual porosity needs an even number of layers (NZ = 2 x matrix layers)");
        if (!g.nodppm) {
            permScale.assign(n, 1.0);
            for (size_t q = size_t(nzh) * nl; q < n; ++q) { permScale[q] = poro[q]; kx[q] *= poro[q]; ky[q] *= poro[q]; kz[q] *= poro[q]; }
        }
    }

    // ---- fault multipliers (MULTFLT) folded into the + face multipliers of the low-side cell
    for (const FaultInput& f : g.faults) {
        for (const FaultFace& ff : f.faces)
            for (int k = std::max(0, ff.k1); k <= std::min(nz - 1, ff.k2); ++k)
                for (int j = std::max(0, ff.j1); j <= std::min(ny - 1, ff.j2); ++j)
                    for (int i = std::max(0, ff.i1); i <= std::min(nx - 1, ff.i2); ++i) {
                        int ii = i, jj = j, kk = k;
                        if (ff.minus) { if (ff.face == 'X') --ii; else if (ff.face == 'Y') --jj; else --kk; }
                        if (ii < 0 || jj < 0 || kk < 0) continue;
                        const int q = g.index(ii, jj, kk);
                        (ff.face == 'X' ? mx : ff.face == 'Y' ? my : mz)[q] *= f.mult;
                    }
    }

    // ---- geometry: 8 corners per cell
    corners.assign(n, {});
    std::vector<std::array<V3, 8>> C(n);
    std::vector<double> dxC, dyC, dzC;      // Cartesian sizes (kept for exact Cartesian transmissibility)
    if (cornerPoint) {
        if (g.coord.size() != size_t(6) * (nx + 1) * (ny + 1))
            throw std::runtime_error("GRID: COORD must have 6*(NX+1)*(NY+1) values");
        if (g.zcorn.size() != 8 * n) throw std::runtime_error("GRID: ZCORN must have 8*NX*NY*NZ values");
        auto pillarPoint = [&](int pi, int pj, double z) -> V3 {
            const double* p = &g.coord[6 * (size_t(pj) * (nx + 1) + pi)];
            const double dz = p[5] - p[2];
            if (std::fabs(dz) < 1e-12) return {p[0], p[1], z};
            const double t = (z - p[2]) / dz;
            return {p[0] + t * (p[3] - p[0]), p[1] + t * (p[4] - p[1]), z};
        };
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int q = g.index(i, j, k);
                    for (int kk = 0; kk < 2; ++kk)
                        for (int jj = 0; jj < 2; ++jj)
                            for (int ii = 0; ii < 2; ++ii) {
                                const size_t zi = (size_t(2 * k + kk) * 2 * ny + (2 * j + jj)) * 2 * nx + 2 * i + ii;
                                C[q][ii + 2 * jj + 4 * kk] = pillarPoint(i + ii, j + jj, g.zcorn[zi]);
                            }
                }
    } else {
        dxC = fillOr(g.dx, n, 0, "DX", true);
        dyC = fillOr(g.dy, n, 0, "DY", true);
        dzC = fillOr(g.dz, n, 0, "DZ", true);
        std::vector<double> top(n, 0.0);
        if (g.tops.size() == nl || g.tops.size() == n) {
            for (size_t q = 0; q < nl; ++q) top[q] = g.tops[q];
            for (int k = 1; k < nz; ++k)
                for (size_t q = 0; q < nl; ++q) {
                    const size_t gi = k * nl + q;
                    top[gi] = (g.tops.size() == n) ? g.tops[gi] : top[gi - nl] + dzC[gi - nl];
                }
        } else throw std::runtime_error("GRID: TOPS must have NX*NY or NX*NY*NZ values");
        std::vector<double> xs(n), ys(n);
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                double x = 0;
                for (int i = 0; i < nx; ++i) { const int q = g.index(i, j, k); xs[q] = x; x += dxC[q]; }
            }
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i) {
                double y = 0;
                for (int j = 0; j < ny; ++j) { const int q = g.index(i, j, k); ys[q] = y; y += dyC[q]; }
            }
        for (size_t q = 0; q < n; ++q) {
            const double x0 = xs[q], x1 = xs[q] + dxC[q], y0 = ys[q], y1 = ys[q] + dyC[q], z0 = top[q], z1 = top[q] + dzC[q];
            C[q] = {V3{x0, y0, z0}, V3{x1, y0, z0}, V3{x0, y1, z0}, V3{x1, y1, z0},
                    V3{x0, y0, z1}, V3{x1, y0, z1}, V3{x0, y1, z1}, V3{x1, y1, z1}};
        }
    }

    // ---- dual porosity: a fracture cell is the same block as its matrix cell -> identical geometry
    if (dual) {
        const size_t off = size_t(nzh) * nl;
        for (size_t q = off; q < n; ++q) {
            C[q] = C[q - off];
            if (!cornerPoint) { dxC[q] = dxC[q - off]; dyC[q] = dyC[q - off]; dzC[q] = dzC[q - off]; }
        }
    }

    // ---- per-cell geometry
    topDepth.assign(n, 0); cellDx.assign(n, 0); cellDy.assign(n, 0); cellDz.assign(n, 0); bulkVol.assign(n, 0); center.assign(n, {});
    std::vector<V3> ctr(n);
    bbMin = {FLT_MAX, FLT_MAX, FLT_MAX}; bbMax = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (size_t q = 0; q < n; ++q) {
        const auto& c8 = C[q];
        V3 m{0, 0, 0};
        for (int v = 0; v < 8; ++v) {
            m = m + c8[v];
            corners[q][v] = {float(c8[v].x), float(c8[v].y), float(c8[v].z)};
        }
        ctr[q] = m * 0.125;
        center[q] = {float(ctr[q].x), float(ctr[q].y), float(ctr[q].z)};
        if (cornerPoint) {
            auto len = [&](int f0, int f1) { const V3 d = faceCenter(c8, f1) - faceCenter(c8, f0); return std::sqrt(dot(d, d)); };
            cellDx[q] = len(0, 1); cellDy[q] = len(2, 3);
            cellDz[q] = 0.25 * ((c8[4].z - c8[0].z) + (c8[5].z - c8[1].z) + (c8[6].z - c8[2].z) + (c8[7].z - c8[3].z));
            bulkVol[q] = (cellDz[q] > 0) ? hexVolume(c8) : 0.0;
            topDepth[q] = 0.25 * (c8[0].z + c8[1].z + c8[2].z + c8[3].z);
        } else {
            cellDx[q] = dxC[q]; cellDy[q] = dyC[q]; cellDz[q] = dzC[q];
            bulkVol[q] = dxC[q] * dyC[q] * dzC[q];
            topDepth[q] = c8[0].z;
        }
    }

    // ---- active cells
    activeOf.assign(n, -1); globalOf.clear(); numPinchedCells = 0; numPinch = 0;
    for (size_t q = 0; q < n; ++q) {
        bool act = g.actnum.empty() || g.actnum[q] != 0;
        if (bulkVol[q] <= 0 || cellDz[q] <= 0) act = false;
        if (poro[q] * ntg[q] <= 0) act = false;
        if (g.minpv > 0 && bulkVol[q] * ntg[q] * poro[q] * U.poreVolume * (q < g.multpv.size() ? g.multpv[q] : 1.0) < g.minpv)
            act = false;   // MINPV
        if (act && g.pinch.on && cellDz[q] < g.pinch.threshold) { act = false; ++numPinchedCells; }
        if (act) { activeOf[q] = int(globalOf.size()); globalOf.push_back(int(q)); }
    }
    if (globalOf.empty()) throw std::runtime_error("GRID: no active cells");
    for (int a = 0; a < numActive(); ++a)
        for (auto& v : corners[globalOf[a]]) {
            bbMin.x = std::min(bbMin.x, v.x); bbMin.y = std::min(bbMin.y, v.y); bbMin.z = std::min(bbMin.z, v.z);
            bbMax.x = std::max(bbMax.x, v.x); bbMax.y = std::max(bbMax.y, v.y); bbMax.z = std::max(bbMax.z, v.z);
        }

    const int na = numActive();
    poreVolRef.resize(na); depth.resize(na); dzCell.resize(na);
    for (int a = 0; a < na; ++a) {
        const int q = globalOf[a];
        poreVolRef[a] = bulkVol[q] * ntg[q] * poro[q] * U.poreVolume;
        if (q < int(g.multpv.size())) poreVolRef[a] *= g.multpv[size_t(q)];
        depth[a] = cornerPoint ? ctr[q].z : topDepth[q] + 0.5 * dzC[q];
        dzCell[a] = cellDz[q];
    }

    // ---- connections
    conns.clear();
    numNNC = 0;
    auto add = [&](int q1, int q2, double t1, double t2, double mult, char dir, bool nnc) {
        const int a = activeOf[q1], b = activeOf[q2];
        if (a < 0 || b < 0 || t1 <= 0 || t2 <= 0 || mult <= 0) return;
        if (dual) {   // geometric neighbours never couple matrix to fracture; DUALPORO matrix blocks are isolated
            const bool m1 = q1 / int(nl) < nzh, m2 = q2 / int(nl) < nzh;
            if (m1 != m2 || (g.dual == 1 && m1)) return;
        }
        Connection cn;
        cn.a = a; cn.b = b; cn.trans = U.cDarcy * mult / (1.0 / t1 + 1.0 / t2);
        cn.dDepth = depth[b] - depth[a]; cn.dir = dir;
        conns.push_back(cn);
        if (nnc) ++numNNC;
    };
    if (!cornerPoint) {
        // OLDTRAN-style harmonic averaging of half-block transmissibilities (k-to-k neighbours)
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int q = g.index(i, j, k);
                    if (activeOf[q] < 0) continue;
                    if (i + 1 < nx) {
                        const int r = g.index(i + 1, j, k);
                        add(q, r, kx[q] * ntg[q] * dyC[q] * dzC[q] / (0.5 * dxC[q]),
                                  kx[r] * ntg[r] * dyC[r] * dzC[r] / (0.5 * dxC[r]), mx[q], 'X', false);
                    }
                    if (j + 1 < ny) {
                        const int r = g.index(i, j + 1, k);
                        add(q, r, ky[q] * ntg[q] * dxC[q] * dzC[q] / (0.5 * dyC[q]),
                                  ky[r] * ntg[r] * dxC[r] * dzC[r] / (0.5 * dyC[r]), my[q], 'Y', false);
                    }
                    if (k + 1 < nz) {
                        const int r = g.index(i, j, k + 1);
                        add(q, r, kz[q] * dxC[q] * dyC[q] / (0.5 * dzC[q]),
                                  kz[r] * dxC[r] * dyC[r] / (0.5 * dzC[r]), mz[q], 'Z', false);
                    }
                }
    } else {
        // NEWTRAN-style: T_half = k * |A . d| / |d|^2 with A the (overlap) face area vector and d centre->face
        auto halfT = [&](int q, double perm, double nt, const V3& A, const V3& fc) {
            const V3 d = fc - ctr[q];
            const double dd = dot(d, d);
            return dd > 0 ? perm * nt * std::fabs(dot(A, d)) / dd : 0.0;
        };
        // lateral connections between neighbouring columns, with throw (non k-to-k) overlaps
        for (int dir = 0; dir < 2; ++dir) {
            const int di = dir == 0 ? 1 : 0, dj = dir == 0 ? 0 : 1;
            // corner indices of the + face of the low cell and the - face of the high cell on pillars A and B
            const int lA_t = dir == 0 ? 1 : 2, lA_b = lA_t + 4, lB_t = 3, lB_b = 7;
            const int hA_t = 0, hA_b = 4, hB_t = dir == 0 ? 2 : 1, hB_b = hB_t + 4;
            for (int j = 0; j + dj < ny; ++j)
                for (int i = 0; i + di < nx; ++i) {
                    int lo = 0;
                    for (int k1 = 0; k1 < nz; ++k1) {
                        if (dual && k1 == nzh) lo = nzh;            // fracture half: restart the depth sweep
                        const int kEnd = dual && k1 < nzh ? nzh : nz;
                        const int q = g.index(i, j, k1);
                        if (activeOf[q] < 0) continue;
                        const auto& L = C[q];
                        const double lTop = std::min(L[lA_t].z, L[lB_t].z), lBot = std::max(L[lA_b].z, L[lB_b].z);
                        for (int k2 = lo; k2 < kEnd; ++k2) {
                            const int r = g.index(i + di, j + dj, k2);
                            const auto& R = C[r];
                            const double rTop = std::min(R[hA_t].z, R[hB_t].z), rBot = std::max(R[hA_b].z, R[hB_b].z);
                            if (rBot <= lTop) { if (k2 == lo) ++lo; continue; }
                            if (rTop >= lBot) break;
                            if (activeOf[r] < 0) continue;
                            const Overlap o = stripOverlap(L[lA_t].z, L[lA_b].z, L[lB_t].z, L[lB_b].z,
                                                           R[hA_t].z, R[hA_b].z, R[hB_t].z, R[hB_b].z);
                            if (o.area <= 1e-12) continue;
                            // geometry of the overlap patch: pillars A,B positions at its centroid depth
                            const V3 pa = L[lA_t] + (L[lA_b] - L[lA_t]) * ((L[lA_b].z - L[lA_t].z) > 1e-12 ? (o.z - L[lA_t].z) / (L[lA_b].z - L[lA_t].z) : 0.0);
                            const V3 pb = L[lB_t] + (L[lB_b] - L[lB_t]) * ((L[lB_b].z - L[lB_t].z) > 1e-12 ? (o.z - L[lB_t].z) / (L[lB_b].z - L[lB_t].z) : 0.0);
                            const V3 e = pb - pa;
                            const double w = std::sqrt(e.x * e.x + e.y * e.y);
                            if (w <= 0) continue;
                            V3 nrm = dir == 0 ? V3{e.y / w, -e.x / w, 0} : V3{-e.y / w, e.x / w, 0};
                            const V3 A = nrm * (o.area * w);
                            V3 fc = pa + e * o.s; fc.z = o.z;
                            const double pq = dir == 0 ? kx[q] : ky[q], pr = dir == 0 ? kx[r] : ky[r];
                            add(q, r, halfT(q, pq, ntg[q], A, fc), halfT(r, pr, ntg[r], A, fc),
                                (dir == 0 ? mx[q] : my[q]), dir == 0 ? 'X' : 'Y', k1 != k2);
                        }
                    }
                }
        }
        // vertical connections (k to k+1 in the same column)
        for (int k = 0; k + 1 < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int q = g.index(i, j, k), r = g.index(i, j, k + 1);
                    if (activeOf[q] < 0 || activeOf[r] < 0) continue;
                    const V3 A = faceAreaVec(C[q], 5);
                    const V3 fc = faceCenter(C[q], 5);
                    add(q, r, halfT(q, kz[q], 1.0, A, fc), halfT(r, kz[r], 1.0, A, fc), mz[q], 'Z', false);
                }
    }
    // PINCH: vertical NNCs across pinched-out / inactive layers. T = 1 / (1/T_top + 1/T_bottom) of the half-cell
    // transmissibilities of the two active cells (TOPBOT), times MULTZ of the upper cell (TOP) or the minimum over the gap (ALL)
    if (g.pinch.on) {
        auto halfZ = [&](int q, int face) {
            const V3 A = faceAreaVec(C[q], face), fc = faceCenter(C[q], face), d = fc - ctr[q];
            const double dd = dot(d, d);
            return dd > 0 ? kz[q] * std::fabs(dot(A, d)) / dd : 0.0;
        };
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                for (int k = 0; k + 1 < nz; ++k) {
                    const int q = g.index(i, j, k);
                    if (activeOf[q] < 0) continue;
                    const int kEnd = dual && k < nzh ? nzh : nz;   // never bridge matrix to fracture
                    double gap = 0, mmin = mz[q];
                    bool ok = true; int r = -1;
                    for (int k2 = k + 1; k2 < kEnd; ++k2) {
                        const int c2 = g.index(i, j, k2);
                        if (activeOf[c2] >= 0) { r = c2; break; }
                        const bool thin = cellDz[c2] < g.pinch.threshold;
                        if (!thin && !g.pinch.gap) { ok = false; break; }
                        gap += std::max(0.0, cellDz[c2]);
                        mmin = std::min(mmin, mz[c2]);
                    }
                    if (!ok || r < 0 || r == g.index(i, j, k + 1) || gap > g.pinch.maxGap) continue;
                    const double t1 = halfZ(q, 5), t2 = halfZ(r, 4);
                    const double mult = g.pinch.multzAll ? mmin : mz[q];
                    if (t1 <= 0 || t2 <= 0 || mult <= 0) continue;
                    Connection cn; cn.a = activeOf[q]; cn.b = activeOf[r];
                    cn.trans = U.cDarcy * mult / (1.0 / t1 + 1.0 / t2);
                    cn.dDepth = depth[cn.b] - depth[cn.a]; cn.dir = 'Z';
                    conns.push_back(cn); ++numNNC; ++numPinch;
                }
    }
    // dual porosity: matrix-fracture coupling T = CDARCY * V * sigma * K  (V bulk volume of the matrix cell);
    // LTOSIGMA: T = CDARCY * V * (fx Kx/lx^2 + fy Ky/ly^2 + fz Kz/lz^2)
    if (dual) {
        auto cellVal = [&](const std::vector<double>& v, size_t q, double def) {
            return q < v.size() && !std::isnan(v[q]) ? v[q] : def;
        };
        for (size_t qm = 0; qm < size_t(nzh) * nl; ++qm) {
            const size_t qf = qm + size_t(nzh) * nl;
            const int a = activeOf[qm], b = activeOf[qf];
            if (a < 0 || b < 0) continue;
            double sk = 0;
            const double lX = cellVal(g.lx, qm, 0), lY = cellVal(g.ly, qm, 0), lZ = cellVal(g.lz, qm, 0);
            if (g.ltosigma && (lX > 0 || lY > 0 || lZ > 0)) {
                if (lX > 0) sk += g.ltoF[0] * kx[qm] / (lX * lX);
                if (lY > 0) sk += g.ltoF[1] * ky[qm] / (lY * lY);
                if (lZ > 0) sk += g.ltoF[2] * kz[qm] / (lZ * lZ);
            } else sk = cellVal(g.sigmav, qm, g.sigma) * kx[qm];
            const double T = U.cDarcy * bulkVol[qm] * sk * cellVal(g.multmf, qm, 1.0);
            if (!(T > 0)) continue;
            Connection cn; cn.a = a; cn.b = b; cn.trans = T; cn.dDepth = depth[b] - depth[a]; cn.dir = 'M';
            conns.push_back(cn); ++numMF;
        }
    }
    // explicit NNCs (transmissibility given directly, already includes the Darcy constant)
    for (const NNCInput& e : g.nnc) {
        if (e.i1 < 0 || e.i1 >= nx || e.j1 < 0 || e.j1 >= ny || e.k1 < 0 || e.k1 >= nz ||
            e.i2 < 0 || e.i2 >= nx || e.j2 < 0 || e.j2 >= ny || e.k2 < 0 || e.k2 >= nz)
            throw std::runtime_error("GRID: NNC cell outside the grid");
        const int a = activeOf[g.index(e.i1, e.j1, e.k1)], b = activeOf[g.index(e.i2, e.j2, e.k2)];
        if (a < 0 || b < 0 || a == b || e.trans <= 0) continue;
        Connection cn; cn.a = a; cn.b = b; cn.trans = e.trans; cn.dDepth = depth[b] - depth[a]; cn.dir = 'N';
        conns.push_back(cn); ++numNNC;
    }
    // merge duplicate cell pairs (keeps the sparsity pattern clean)
    for (auto& cn : conns) if (cn.a > cn.b) { std::swap(cn.a, cn.b); cn.dDepth = -cn.dDepth; }
    std::sort(conns.begin(), conns.end(), [](const Connection& x, const Connection& y) { return x.a != y.a ? x.a < y.a : x.b < y.b; });
    size_t m = 0;
    for (size_t i = 0; i < conns.size(); ++i) {
        if (m > 0 && conns[m - 1].a == conns[i].a && conns[m - 1].b == conns[i].b) conns[m - 1].trans += conns[i].trans;
        else conns[m++] = conns[i];
    }
    conns.resize(m);

    // ---- adjacency
    adjStart.assign(na + 1, 0);
    for (auto& cn : conns) { adjStart[cn.a + 1]++; adjStart[cn.b + 1]++; }
    for (int a = 0; a < na; ++a) adjStart[a + 1] += adjStart[a];
    adjConn.assign(adjStart[na], 0); adjNbr.assign(adjStart[na], 0);
    std::vector<int> fill(adjStart.begin(), adjStart.end() - 1);
    for (int ci = 0; ci < int(conns.size()); ++ci) {
        const auto& cn = conns[ci];
        adjConn[fill[cn.a]] = ci; adjNbr[fill[cn.a]++] = cn.b;
        adjConn[fill[cn.b]] = ci; adjNbr[fill[cn.b]++] = cn.a;
    }
}

}  // namespace ressim
