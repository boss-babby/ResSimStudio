// Grid processing: Cartesian (DX/DY/DZ/TOPS) and corner-point (COORD/ZCORN) geometry,
// active cells, pore volumes, transmissibilities (incl. fault throw NNCs, MULTFLT, explicit NNC).
#include "ressim/grid.h"
#include "ressim/tools.h"
#include <stdexcept>
#include <string>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <algorithm>
#include <map>

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

std::vector<double> fillOr(const std::vector<double>& v, size_t n, double def, const char* name, bool required, size_t nAll = 0) {
    if (nAll < n) nAll = n;
    if (v.empty()) {
        if (required) throw std::runtime_error(std::string("GRID: missing ") + name);
        return std::vector<double>(nAll, def);
    }
    if (v.size() == nAll) return v;
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

// ---------------------------------------------------------------- local grid refinement: geometry and connections
// Local cells are trilinear subdivisions of their host hexahedron. Every host connection (neighbour, fault NNC, pinch-out)
// is redistributed over the local cells on the host face: the share of local pair (c, d) is its face overlap fraction
// times the ratio of the per-area two-point conductances, T_cd = T_AB * w_cd * [1/(1/tau_c + 1/tau_d)] / [1/(1/tau_A + 1/tau_B)],
// with tau the half-cell transmissibility per face area. All multipliers of the host connection (MULTX, MULTFLT,
// MULTREGT, ...) are thereby kept, also between two hosts of the same LGR. Local cells of one host are joined with their
// own half transmissibilities (no multiplier: MULTX etc. act on host faces).
static void buildLgrImpl(Grid& G, const Case& c, std::vector<std::array<V3, 8>>& C, std::vector<V3>& ctr,
                         const std::vector<double>& poro, const std::vector<double>& kx, const std::vector<double>& ky,
                         const std::vector<double>& kz, const std::vector<double>& ntg, const std::vector<double>& mx,
                         const std::vector<double>& my, const std::vector<double>& mz) {
    const GridInput& g = c.grid;
    const Units U = Units::of(c.units);
    const int n = G.nx * G.ny * G.nz;
    const bool cp = G.cornerPoint;
    int nL = 0;
    for (const LgrInput& L : g.lgrs) nL += L.size();
    const size_t nAll = size_t(n) + size_t(nL);
    G.hostLgr.assign(size_t(n), -1); G.lgrHost.assign(size_t(nL), -1); G.lgr.clear();
    C.resize(nAll); ctr.resize(nAll); G.corners.resize(nAll); G.center.resize(nAll);
    G.topDepth.resize(nAll, 0); G.cellDx.resize(nAll, 0); G.cellDy.resize(nAll, 0); G.cellDz.resize(nAll, 0); G.bulkVol.resize(nAll, 0);
    std::vector<std::array<double, 6>> fr(static_cast<size_t>(nL));          // fx0 fx1 fy0 fy1 fz0 fz1 inside the host
    std::vector<std::array<char, 6>> onFace(static_cast<size_t>(nL));        // local cell lies on host face f
    std::vector<std::vector<int>> kids(static_cast<size_t>(n));            // host -> local cells (only hosts filled)
    int off = 0;
    for (size_t l = 0; l < g.lgrs.size(); ++l) {
        const LgrInput& L = g.lgrs[l];
        std::string why;
        if (!L.valid(&why)) throw std::runtime_error("GRID: " + why);
        if (L.i2 >= G.nx || L.j2 >= G.ny || L.k2 >= G.nz) throw std::runtime_error("GRID: LGR " + L.name + " host box is outside the grid");
        for (int k = L.k1; k <= L.k2; ++k)
            for (int j = L.j1; j <= L.j2; ++j)
                for (int i = L.i1; i <= L.i2; ++i)
                    if (G.hostLgr[size_t(g.index(i, j, k))] >= 0)
                        throw std::runtime_error("GRID: LGR " + L.name + " overlaps LGR " + g.lgrs[size_t(G.hostLgr[size_t(g.index(i, j, k))])].name);
        std::vector<int> H[3]; std::vector<double> F0[3], F1[3];
        for (int d = 0; d < 3; ++d) L.split(d, H[d], F0[d], F1[d]);
        auto first = [&](int d, int m) { return m == 0 || H[d][size_t(m - 1)] != H[d][size_t(m)]; };
        auto last = [&](int d, int m) { const int N = d == 0 ? L.nx : d == 1 ? L.ny : L.nz; return m == N - 1 || H[d][size_t(m + 1)] != H[d][size_t(m)]; };
        Grid::LgrInfo I; I.name = L.name; I.offset = off; I.nx = L.nx; I.ny = L.ny; I.nz = L.nz;
        I.i1 = L.i1; I.i2 = L.i2; I.j1 = L.j1; I.j2 = L.j2; I.k1 = L.k1; I.k2 = L.k2;
        G.lgr.push_back(I);
        for (int k = 0; k < L.nz; ++k)
            for (int j = 0; j < L.ny; ++j)
                for (int i = 0; i < L.nx; ++i) {
                    const int o = off + L.index(i, j, k), gl = n + o;
                    const int h = g.index(H[0][size_t(i)], H[1][size_t(j)], H[2][size_t(k)]);
                    G.lgrHost[size_t(o)] = h; G.hostLgr[size_t(h)] = int(l); kids[size_t(h)].push_back(gl);
                    const std::array<double, 6> f = {F0[0][size_t(i)], F1[0][size_t(i)], F0[1][size_t(j)], F1[1][size_t(j)], F0[2][size_t(k)], F1[2][size_t(k)]};
                    fr[size_t(o)] = f;
                    onFace[size_t(o)] = {char(first(0, i)), char(last(0, i)), char(first(1, j)), char(last(1, j)), char(first(2, k)), char(last(2, k))};
                    const auto& Hc = C[size_t(h)];
                    for (int m = 0; m < 8; ++m) {
                        const double u = (m & 1) ? f[1] : f[0], v = (m & 2) ? f[3] : f[2], w = (m & 4) ? f[5] : f[4];
                        V3 r{0, 0, 0};
                        for (int q = 0; q < 8; ++q)
                            r = r + Hc[size_t(q)] * (((q & 1) ? u : 1 - u) * ((q & 2) ? v : 1 - v) * ((q & 4) ? w : 1 - w));
                        C[size_t(gl)][size_t(m)] = r;
                    }
                }
        off += L.size();
    }
    // per-cell geometry of the local cells
    for (size_t gl = size_t(n); gl < nAll; ++gl) {
        const auto& c8 = C[gl];
        V3 m{0, 0, 0};
        for (int v = 0; v < 8; ++v) { m = m + c8[size_t(v)]; G.corners[gl][size_t(v)] = {float(c8[size_t(v)].x), float(c8[size_t(v)].y), float(c8[size_t(v)].z)}; }
        ctr[gl] = m * 0.125;
        G.center[gl] = {float(ctr[gl].x), float(ctr[gl].y), float(ctr[gl].z)};
        const auto& f = fr[gl - size_t(n)];
        const int h = G.lgrHost[gl - size_t(n)];
        if (!cp) {
            G.cellDx[gl] = G.cellDx[size_t(h)] * (f[1] - f[0]); G.cellDy[gl] = G.cellDy[size_t(h)] * (f[3] - f[2]);
            G.cellDz[gl] = G.cellDz[size_t(h)] * (f[5] - f[4]);
            G.bulkVol[gl] = G.cellDx[gl] * G.cellDy[gl] * G.cellDz[gl];
            G.topDepth[gl] = c8[0].z;
        } else {
            auto len = [&](int f0, int f1) { const V3 d = faceCenter(c8, f1) - faceCenter(c8, f0); return std::sqrt(dot(d, d)); };
            G.cellDx[gl] = len(0, 1); G.cellDy[gl] = len(2, 3);
            G.cellDz[gl] = 0.25 * ((c8[4].z - c8[0].z) + (c8[5].z - c8[1].z) + (c8[6].z - c8[2].z) + (c8[7].z - c8[3].z));
            G.bulkVol[gl] = G.cellDz[gl] > 0 ? hexVolume(c8) : 0.0;
            G.topDepth[gl] = 0.25 * (c8[0].z + c8[1].z + c8[2].z + c8[3].z);
        }
    }
    // active cells: hosts leave, local cells of active hosts join (in local order after the global cells)
    const std::vector<int> oldGlobal = G.globalOf;
    const int naOld = int(oldGlobal.size());
    std::vector<int> remap(size_t(naOld), -1), glob;
    std::vector<double> pvN, depN, dzN;
    for (int a = 0; a < naOld; ++a) {
        const int q = oldGlobal[size_t(a)];
        if (G.hostLgr[size_t(q)] >= 0) continue;
        remap[size_t(a)] = int(glob.size()); glob.push_back(q);
        pvN.push_back(G.poreVolRef[size_t(a)]); depN.push_back(G.depth[size_t(a)]); dzN.push_back(G.dzCell[size_t(a)]);
    }
    std::vector<int> act(nAll, -1);
    for (int q = 0; q < n; ++q) if (G.activeOf[size_t(q)] >= 0) act[size_t(q)] = remap[size_t(G.activeOf[size_t(q)])];
    for (size_t gl = size_t(n); gl < nAll; ++gl) {
        const int h = G.lgrHost[gl - size_t(n)];
        bool on = G.activeOf[size_t(h)] >= 0;
        if (gl < g.actnum.size() && g.actnum[gl] == 0) on = false;
        if (G.bulkVol[gl] <= 0 || G.cellDz[gl] <= 0 || poro[gl] * ntg[gl] <= 0) on = false;
        if (!on) continue;
        double pv = G.bulkVol[gl] * ntg[gl] * poro[gl] * U.poreVolume;
        if (gl < g.multpv.size()) pv *= g.multpv[gl];
        act[gl] = int(glob.size()); glob.push_back(int(gl));
        pvN.push_back(pv); depN.push_back(ctr[gl].z); dzN.push_back(G.cellDz[gl]);
    }
    G.activeOf = act; G.globalOf = glob; G.poreVolRef = pvN; G.depth = depN; G.dzCell = dzN;

    // connections
    auto tau = [&](int q, int f) {
        const std::vector<double>& K = f < 2 ? kx : f < 4 ? ky : kz;
        const double nt = f < 4 ? ntg[size_t(q)] : 1.0;
        const V3 A = faceAreaVec(C[size_t(q)], f), fc = faceCenter(C[size_t(q)], f), d = fc - ctr[size_t(q)];
        const double dd = dot(d, d), aa = std::sqrt(dot(A, A));
        return (dd > 0 && aa > 0) ? K[size_t(q)] * nt * std::fabs(dot(A, d)) / dd / aa : 0.0;
    };
    struct El { int g; double l0, l1, s0, s1; };
    auto faceDepths = [&](int q, int f, double& top, double& bot) {
        top = bot = 0;
        for (int e = 0; e < 4; ++e) { const int v = kFaces[f][e]; (v < 4 ? top : bot) += 0.5 * C[size_t(q)][size_t(v)].z; }
    };
    auto elems = [&](int q, int f, std::vector<El>& out) {
        out.clear();
        if (G.hostLgr[size_t(q)] < 0) {
            El e{q, 0, 1, 0, 1};
            if (f < 4) faceDepths(q, f, e.s0, e.s1);
            out.push_back(e);
            return;
        }
        for (int gl : kids[size_t(q)]) {
            const size_t o = size_t(gl - n);
            if (!onFace[o][size_t(f)]) continue;
            const auto& r = fr[o];
            El e{gl, 0, 1, 0, 1};
            if (f < 2) { e.l0 = r[2]; e.l1 = r[3]; faceDepths(gl, f, e.s0, e.s1); }
            else if (f < 4) { e.l0 = r[0]; e.l1 = r[1]; faceDepths(gl, f, e.s0, e.s1); }
            else { e.l0 = r[0]; e.l1 = r[1]; e.s0 = r[2]; e.s1 = r[3]; }
            out.push_back(e);
        }
    };
    auto nearest = [&](int q, const V3& t) {
        int best = -1; double bd = 1e300;
        for (int gl : kids[size_t(q)]) {
            if (G.activeOf[size_t(gl)] < 0) continue;
            const V3 d = ctr[size_t(gl)] - t; const double dd = dot(d, d);
            if (dd < bd) { bd = dd; best = gl; }
        }
        return best;
    };
    std::vector<Connection> out;
    out.reserve(G.conns.size() + size_t(nL) * 3);
    auto push = [&](int ga, int gb, double T, char dir) {
        if (ga < 0 || gb < 0) return;
        const int a = G.activeOf[size_t(ga)], b = G.activeOf[size_t(gb)];
        if (a < 0 || b < 0 || a == b || !(T > 0)) return;
        Connection cn; cn.a = a; cn.b = b; cn.trans = T; cn.dDepth = G.depth[size_t(b)] - G.depth[size_t(a)]; cn.dir = dir;
        out.push_back(cn);
    };
    std::vector<El> SA, SB;
    for (const Connection& cn : G.conns) {
        const int ga = oldGlobal[size_t(cn.a)], gb = oldGlobal[size_t(cn.b)];
        const int hA = G.hostLgr[size_t(ga)], hB = G.hostLgr[size_t(gb)];
        if (hA < 0 && hB < 0) { push(ga, gb, cn.trans, cn.dir); continue; }
        int fA = -1;
        if (cn.dir == 'X' || cn.dir == 'Y' || cn.dir == 'Z') {
            int ia, ja, ka, ib, jb, kb;
            G.ijk(ga, ia, ja, ka); G.ijk(gb, ib, jb, kb);
            if (cn.dir == 'X' && ia != ib) fA = ib > ia ? 1 : 0;
            else if (cn.dir == 'Y' && ja != jb) fA = jb > ja ? 3 : 2;
            else if (cn.dir == 'Z' && ka != kb) fA = kb > ka ? 5 : 4;
        }
        if (fA < 0) {   // explicit NNC / unknown face: attach to the nearest local cell(s)
            const int xa = hA >= 0 ? nearest(ga, ctr[size_t(gb)]) : ga;
            const int xb = hB >= 0 ? nearest(gb, xa >= 0 ? ctr[size_t(xa)] : ctr[size_t(ga)]) : gb;
            push(xa, xb, cn.trans, cn.dir);
            continue;
        }
        const int fB = fA ^ 1;
        elems(ga, fA, SA); elems(gb, fB, SB);
        const double tA = tau(ga, fA), tB = tau(gb, fB);
        const double base = (tA > 0 && tB > 0) ? 1.0 / (1.0 / tA + 1.0 / tB) : 0.0;
        auto weight = [&](const El& x, const El& y, bool depthToo) {
            const double lo = std::max(0.0, std::min(x.l1, y.l1) - std::max(x.l0, y.l0));
            if (!depthToo) return lo;
            return lo * std::max(0.0, std::min(x.s1, y.s1) - std::max(x.s0, y.s0));
        };
        double tot = 0;
        for (const El& x : SA) for (const El& y : SB) tot += weight(x, y, true);
        const bool useDepth = tot > 0;
        if (!useDepth) for (const El& x : SA) for (const El& y : SB) tot += weight(x, y, false);
        if (!(tot > 0)) continue;
        for (const El& x : SA)
            for (const El& y : SB) {
                const double w = weight(x, y, useDepth);
                if (!(w > 0)) continue;
                const double tc = tau(x.g, fA), td = tau(y.g, fB);
                const double r = (tc > 0 && td > 0 && base > 0) ? (1.0 / (1.0 / tc + 1.0 / td)) / base : 0.0;
                push(x.g, y.g, cn.trans * (w / tot) * r, cn.dir);
            }
    }
    // inside every LGR: local neighbours
    for (size_t l = 0; l < G.lgr.size(); ++l) {
        const Grid::LgrInfo& L = G.lgr[l];
        for (int k = 0; k < L.nz; ++k)
            for (int j = 0; j < L.ny; ++j)
                for (int i = 0; i < L.nx; ++i) {
                    const int g1 = G.localCell(int(l), i, j, k);
                    if (G.activeOf[size_t(g1)] < 0) continue;
                    for (int d = 0; d < 3; ++d) {
                        const int i2 = i + (d == 0), j2 = j + (d == 1), k2 = k + (d == 2);
                        if (i2 >= L.nx || j2 >= L.ny || k2 >= L.nz) continue;
                        const int g2 = G.localCell(int(l), i2, j2, k2);
                        if (G.activeOf[size_t(g2)] < 0) continue;
                        // local cells of different hosts are joined through the redistributed host connection (keeps
                        // fault throws, MULTFLT / MULTX / MULTREGT inside the refinement)
                        if (G.lgrHost[size_t(g1 - n)] != G.lgrHost[size_t(g2 - n)]) continue;
                        const int f = 2 * d + 1;
                        const std::vector<double>& K = d == 0 ? kx : d == 1 ? ky : kz;
                        const V3 A = faceAreaVec(C[size_t(g1)], f), fc = faceCenter(C[size_t(g1)], f);
                        auto half = [&](int q) {
                            const V3 dv = fc - ctr[size_t(q)]; const double dd = dot(dv, dv);
                            return dd > 0 ? K[size_t(q)] * (d < 2 ? ntg[size_t(q)] : 1.0) * std::fabs(dot(A, dv)) / dd : 0.0;
                        };
                        const double t1 = half(g1), t2 = half(g2);
                        const double m = (d == 0 ? mx : d == 1 ? my : mz)[size_t(g1)];
                        if (t1 <= 0 || t2 <= 0 || m <= 0) continue;
                        push(g1, g2, U.cDarcy * m / (1.0 / t1 + 1.0 / t2), "XYZ"[d]);
                    }
                }
    }
    G.conns.swap(out);
}

// ---------------------------------------------------------------- local grid refinement input
void LgrInput::split(int dir, std::vector<int>& host, std::vector<double>& f0, std::vector<double>& f1) const {
    const int a = dir == 0 ? i1 : dir == 1 ? j1 : k1, b = dir == 0 ? i2 : dir == 1 ? j2 : k2;
    const int N = dir == 0 ? nx : dir == 1 ? ny : nz, H = b - a + 1;
    const std::vector<int>& fin = dir == 0 ? nxfin : dir == 1 ? nyfin : nzfin;
    const std::vector<double>& hw = dir == 0 ? hxfin : dir == 1 ? hyfin : hzfin;
    host.assign(size_t(N), a); f0.assign(size_t(N), 0.0); f1.assign(size_t(N), 1.0);
    std::vector<int> cnt(size_t(std::max(H, 0)), H > 0 ? N / H : 0);
    if (int(fin.size()) == H) cnt = fin;
    int l = 0;
    for (int h = 0; h < H && l < N; ++h) {
        const int l0 = l, l1 = std::min(N, l + std::max(cnt[size_t(h)], 0));
        double tot = 0;
        for (int m = l0; m < l1; ++m) tot += (int(hw.size()) == N && hw[size_t(m)] > 0) ? hw[size_t(m)] : 1.0;
        double acc = 0;
        for (int m = l0; m < l1; ++m) {
            const double w = (int(hw.size()) == N && hw[size_t(m)] > 0) ? hw[size_t(m)] : 1.0;
            host[size_t(m)] = a + h; f0[size_t(m)] = acc / tot; acc += w; f1[size_t(m)] = (m + 1 == l1) ? 1.0 : acc / tot;
        }
        l = l1;
    }
}

bool LgrInput::valid(std::string* why) const {
    auto bad = [&](const std::string& m) { if (why) *why = "LGR " + name + ": " + m; return false; };
    if (i1 < 0 || j1 < 0 || k1 < 0 || i2 < i1 || j2 < j1 || k2 < k1) return bad("invalid host box");
    if (nx < 1 || ny < 1 || nz < 1) return bad("NX, NY and NZ must be >= 1");
    const int H[3] = {i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1}, N[3] = {nx, ny, nz};
    const std::vector<int>* F[3] = {&nxfin, &nyfin, &nzfin};
    const std::vector<double>* W[3] = {&hxfin, &hyfin, &hzfin};
    const char* D = "XYZ";
    for (int d = 0; d < 3; ++d) {
        if (!F[d]->empty()) {
            if (int(F[d]->size()) != H[d]) return bad(std::string("N") + D[d] + "FIN needs one value per host cell");
            int sum = 0;
            for (int v : *F[d]) { if (v < 1) return bad(std::string("N") + D[d] + "FIN values must be >= 1"); sum += v; }
            if (sum != N[d]) return bad(std::string("N") + D[d] + "FIN must add up to N" + D[d]);
        } else if (N[d] % H[d] != 0)
            return bad(std::string("N") + D[d] + " = " + std::to_string(N[d]) + " is not a multiple of the " + std::to_string(H[d]) +
                       " host cells (give N" + D[d] + "FIN)");
        if (!W[d]->empty() && int(W[d]->size()) != N[d]) return bad(std::string("H") + D[d] + "FIN needs N" + D[d] + " values");
        for (double v : *W[d]) if (!(v > 0)) return bad(std::string("H") + D[d] + "FIN values must be positive");
    }
    if (parent != "GLOBAL" && !parent.empty()) return bad("nested refinements (parent " + parent + ") are not modelled");
    return true;
}

namespace {
struct CellArrayRef { const char* name; std::vector<double>* d; std::vector<int>* i; double def; bool inherit; };
std::vector<CellArrayRef> cellArrays(Case& c) {
    GridInput& g = c.grid;
    const double NaN = std::nan("");
    return {
        {"PORO", &g.poro, nullptr, 0, true}, {"PERMX", &g.permx, nullptr, 0, true}, {"PERMY", &g.permy, nullptr, NaN, true},
        {"PERMZ", &g.permz, nullptr, NaN, true}, {"NTG", &g.ntg, nullptr, 1, true}, {"MULTX", &g.multx, nullptr, 1, false},
        {"MULTY", &g.multy, nullptr, 1, false}, {"MULTZ", &g.multz, nullptr, 1, false}, {"MULTPV", &g.multpv, nullptr, 1, true},
        {"SWATINIT", &g.swatinit, nullptr, -1, true}, {"SWL", &g.swl, nullptr, NaN, true}, {"SWCR", &g.swcr, nullptr, NaN, true},
        {"SWU", &g.swu, nullptr, NaN, true}, {"SGL", &g.sgl, nullptr, NaN, true}, {"SGCR", &g.sgcr, nullptr, NaN, true},
        {"SGU", &g.sgu, nullptr, NaN, true}, {"SOWCR", &g.sowcr, nullptr, NaN, true}, {"SOGCR", &g.sogcr, nullptr, NaN, true},
        {"ISWL", &g.iswl, nullptr, NaN, true}, {"ISWCR", &g.iswcr, nullptr, NaN, true}, {"ISWU", &g.iswu, nullptr, NaN, true},
        {"ISGL", &g.isgl, nullptr, NaN, true}, {"ISGCR", &g.isgcr, nullptr, NaN, true}, {"ISGU", &g.isgu, nullptr, NaN, true},
        {"ISOWCR", &g.isowcr, nullptr, NaN, true}, {"ISOGCR", &g.isogcr, nullptr, NaN, true}, {"SPOLY", &g.spoly, nullptr, 0, true},
        {"SIGMAV", &g.sigmav, nullptr, NaN, true}, {"LX", &g.lx, nullptr, 0, true}, {"LY", &g.ly, nullptr, 0, true},
        {"LZ", &g.lz, nullptr, 0, true}, {"MULTMF", &g.multmf, nullptr, 1, true}, {"DZMTRXV", &g.dzmtrxv, nullptr, 0, true}, {"DISPERC", &g.disperc, nullptr, 0, true},
        {"PRESSURE", &c.init.pressure, nullptr, NaN, true}, {"SWAT", &c.init.swat, nullptr, NaN, true},
        {"SGAS", &c.init.sgas, nullptr, NaN, true}, {"RS", &c.init.rs, nullptr, NaN, true}, {"RV", &c.init.rv, nullptr, NaN, true},
        {"ACTNUM", nullptr, &g.actnum, 1, true}, {"SATNUM", nullptr, &g.satnum, 1, true}, {"PVTNUM", nullptr, &g.pvtnum, 1, true},
        {"EQLNUM", nullptr, &g.eqlnum, 1, true}, {"FIPNUM", nullptr, &g.fipnum, 1, true}, {"IMBNUM", nullptr, &g.imbnum, 1, true},
        {"PLMIXNUM", nullptr, &g.plmixnum, 1, true}, {"MULTNUM", nullptr, &g.multnum, 1, true}, {"FLUXNUM", nullptr, &g.fluxnum, 1, true},
        {"OPERNUM", nullptr, &g.opernum, 1, true}};
}
}  // namespace

void expandLgrArrays(Case& c) {
    GridInput& g = c.grid;
    if (g.lgrExpanded || g.lgrs.empty()) return;
    const size_t n = size_t(g.size());
    size_t nL = 0;
    std::vector<std::vector<int>> hostOf(g.lgrs.size());
    for (size_t l = 0; l < g.lgrs.size(); ++l) {
        const LgrInput& L = g.lgrs[l];
        std::string why;
        if (!L.valid(&why)) throw std::runtime_error("GRID: " + why);
        if (L.i2 >= g.nx || L.j2 >= g.ny || L.k2 >= g.nz) throw std::runtime_error("GRID: LGR " + L.name + " host box is outside the grid");
        std::vector<int> hi, hj, hk; std::vector<double> a, b;
        L.split(0, hi, a, b); L.split(1, hj, a, b); L.split(2, hk, a, b);
        auto& H = hostOf[l]; H.resize(size_t(L.size()));
        for (int k = 0; k < L.nz; ++k)
            for (int j = 0; j < L.ny; ++j)
                for (int i = 0; i < L.nx; ++i) H[size_t(L.index(i, j, k))] = g.index(hi[size_t(i)], hj[size_t(j)], hk[size_t(k)]);
        nL += size_t(L.size());
    }
    for (const CellArrayRef& r : cellArrays(c)) {
        bool any = false;
        for (const LgrInput& L : g.lgrs) if (L.props.count(r.name)) any = true;
        const size_t cur = r.d ? r.d->size() : r.i->size();
        if (cur != n) {
            if (!any || cur != 0) continue;   // absent (and not given locally) or malformed: left alone
            if (r.d) {   // only given inside a refinement: the global cells take the default (PERMY/PERMZ: PERMX)
                const bool kyz = (r.d == &g.permy || r.d == &g.permz) && g.permx.size() == n;
                if (kyz) *r.d = g.permx; else r.d->assign(n, r.def);
            } else r.i->assign(n, int(r.def));
        }
        for (size_t l = 0; l < g.lgrs.size(); ++l) {
            const LgrInput& L = g.lgrs[l];
            auto it = L.props.find(r.name);
            const std::vector<double>* loc = it != L.props.end() ? &it->second : nullptr;
            for (int m = 0; m < L.size(); ++m) {
                const int h = hostOf[l][size_t(m)];
                const bool own = loc && size_t(m) < loc->size() && !std::isnan((*loc)[size_t(m)]);
                if (r.d) {
                    double v = own ? (*loc)[size_t(m)] : r.inherit ? (*r.d)[size_t(h)] : r.def;
                    r.d->push_back(v);
                } else {
                    int v = own ? int(std::lround((*loc)[size_t(m)])) : r.inherit ? (*r.i)[size_t(h)] : int(r.def);
                    r.i->push_back(v);
                }
            }
        }
    }
    (void)nL;
    g.lgrExpanded = 1;
}

// MINC proximity function of a matrix block with fracture sets of spacing l[d] (fr[d]): V(x) = fraction of the block
// volume within distance x of a fracture, a_d(x) = dV/dx contribution of fracture set d (interface area per volume)
namespace {
struct MincBlock {
    double l[3] = {0, 0, 0}; bool fr[3] = {false, false, false};
    double xmax() const { double m = 1e300; for (int d = 0; d < 3; ++d) if (fr[d]) m = std::min(m, 0.5 * l[d]); return m; }
    double V(double x) const { double p = 1; for (int d = 0; d < 3; ++d) if (fr[d]) p *= std::max(0.0, 1 - 2 * x / l[d]); return 1 - p; }
    double a(int d, double x) const {
        if (!fr[d]) return 0;
        double p = 2 / l[d];
        for (int e = 0; e < 3; ++e) if (e != d && fr[e]) p *= std::max(0.0, 1 - 2 * x / l[e]);
        return p;
    }
    double depthOf(double F) const {   // x with V(x) = F (bisection, V monotone)
        double lo = 0, hi = xmax();
        for (int it = 0; it < 200; ++it) { const double m = 0.5 * (lo + hi); (V(m) < F ? lo : hi) = m; }
        return 0.5 * (lo + hi);
    }
};
}  // namespace

void Grid::buildMinc(const GridInput& g, double cDarcy, size_t n, size_t nl, int nzh, int nm, const std::vector<double>& kx,
                     const std::vector<double>& ky, const std::vector<double>& kz) {
    // volume fractions (outer -> inner)
    std::vector<double> f(size_t(nm), 0.0);
    if (g.mincVol.size() == size_t(nm)) f = g.mincVol;
    else {   // geometric, thin outer shells (transient near the fractures) without a coarse core: ratio 1 + 3 / n
        const double r = 1.0 + 3.0 / nm;
        for (int j = 0; j < nm; ++j) f[size_t(j)] = std::pow(r, j);
    }
    double sum = 0; for (double v : f) sum += v;
    for (double& v : f) v /= sum;
    mincN = nm; mincFrac = f;
    const size_t nh = size_t(nzh) * nl, ns = size_t(nm - 1), nM = nh * ns, base = n + numLgrCells();
    activeOf.resize(base + nM, -1);
    topDepth.resize(base + nM, 0); cellDx.resize(base + nM, 0); cellDy.resize(base + nM, 0); cellDz.resize(base + nM, 0);
    bulkVol.resize(base + nM, 0); center.resize(base + nM); corners.resize(base + nM);
    if (!permScale.empty()) permScale.resize(base + nM, 1.0);
    mincHost.assign(nM, -1); mincShell.assign(nM, 0);
    auto cellVal = [&](const std::vector<double>& v, size_t q, double def) { return q < v.size() && !std::isnan(v[q]) ? v[q] : def; };
    const bool dirOn[3][3] = {{false, false, true}, {true, true, false}, {true, true, true}};   // dims 1: Z, 2: X+Y, 3: all
    for (size_t qm = 0; qm < nh; ++qm) {
        const double V0 = bulkVol[qm];
        for (size_t s = 0; s < ns; ++s) {
            const size_t gq = base + qm * ns + s;
            mincHost[gq - base] = int(qm); mincShell[gq - base] = int(s) + 2;
            topDepth[gq] = topDepth[qm]; cellDx[gq] = cellDx[qm]; cellDy[gq] = cellDy[qm]; cellDz[gq] = cellDz[qm];
            bulkVol[gq] = V0 * f[s + 1];
            center[gq] = center[qm];
            for (auto& v : corners[gq]) v = center[qm];   // no own geometry: drawn as a point at the block centre
        }
        const int a = activeOf[qm];
        if (a < 0) continue;
        const int b = qm + nh < n ? activeOf[qm + nh] : -1;
        // block geometry: LX / LY / LZ, else a block from the shape factor (sigma = 4 dims / L^2)
        MincBlock B;
        const double lv[3] = {cellVal(g.lx, qm, 0), cellVal(g.ly, qm, 0), cellVal(g.lz, qm, 0)};
        bool any = false;
        for (int d = 0; d < 3; ++d) if (lv[d] > 0) { B.l[d] = lv[d]; B.fr[d] = true; any = true; }
        if (!any) {
            const double sg = cellVal(g.sigmav, qm, g.sigma);
            if (!(sg > 0)) throw std::runtime_error("GRID: MINC needs the matrix block size (LX / LY / LZ) or SIGMA > 0");
            const double L = std::sqrt(4.0 * g.mincDims / sg);
            for (int d = 0; d < 3; ++d) if (dirOn[g.mincDims - 1][d]) { B.l[d] = L; B.fr[d] = true; }
        }
        // shell boundaries x_1 .. x_nm (x_0 = 0 at the fracture), nodes at mid-depth of each shell
        std::vector<double> x(size_t(nm) + 1, 0.0), y(size_t(nm) + 1, 0.0);
        double F = 0;
        for (int j = 1; j < nm; ++j) { F += f[size_t(j - 1)]; x[size_t(j)] = B.depthOf(F); }
        x[size_t(nm)] = B.xmax();
        for (int j = 1; j <= nm; ++j) y[size_t(j)] = 0.5 * (x[size_t(j - 1)] + x[size_t(j)]);
        const double kd[3] = {kx[qm], ky[qm], kz[qm]}, mmf = cellVal(g.multmf, qm, 1.0);
        auto Tj = [&](int j, bool geomOnly = false) {   // between node j-1 (0 = fracture) and node j (geomOnly: A / d, no CDARCY K)
            double ak = 0;
            for (int d = 0; d < 3; ++d) ak += (geomOnly ? 1.0 : kd[d]) * B.a(d, x[size_t(j - 1)]);
            return (geomOnly ? 1.0 : cDarcy) * V0 * ak / std::max(1e-12, y[size_t(j)] - y[size_t(j - 1)]) * mmf;
        };
        // inner shells become active cells
        const double pv0 = poreVolRef[size_t(a)];
        std::vector<int> act(size_t(nm), a);
        for (size_t s = 0; s < ns; ++s) {
            const size_t gq = base + qm * ns + s;
            activeOf[gq] = int(globalOf.size()); globalOf.push_back(int(gq));
            poreVolRef.push_back(pv0 * f[s + 1]); depth.push_back(depth[size_t(a)]); dzCell.push_back(dzCell[size_t(a)]);
            act[s + 1] = activeOf[gq];
        }
        poreVolRef[size_t(a)] = pv0 * f[0];
        bulkVol[qm] = V0 * f[0];
        for (int j = 1; j <= nm; ++j) {
            const int p = j == 1 ? b : act[size_t(j - 2)], q = act[size_t(j - 1)];
            const double T = Tj(j);
            if (p < 0 || !(T > 0)) continue;
            Connection cn; cn.a = p; cn.b = q; cn.trans = T; cn.dDepth = 0; cn.dir = 'M';
            cn.geo = Tj(j, true);
            conns.push_back(cn); ++numMF;
        }
    }
}

int expandMincArrays(Case& c) {
    GridInput& g = c.grid;
    if (g.minc <= 1 || g.dual != 1 || g.mincExpanded || g.nz % 2) return 0;
    if (!g.lgrs.empty()) throw std::runtime_error("GRID: MINC with local grid refinements is not modelled");
    const size_t n = size_t(g.size()), half = n / 2, ns = size_t(g.minc - 1);
    for (const CellArrayRef& r : cellArrays(c)) {
        const size_t cur = r.d ? r.d->size() : r.i->size();
        if (cur != n) continue;
        if (r.d) { r.d->reserve(n + half * ns); for (size_t q = 0; q < half; ++q) for (size_t s = 0; s < ns; ++s) r.d->push_back(r.inherit ? (*r.d)[q] : r.def); }
        else { r.i->reserve(n + half * ns); for (size_t q = 0; q < half; ++q) for (size_t s = 0; s < ns; ++s) r.i->push_back(r.inherit ? (*r.i)[q] : int(r.def)); }
    }
    g.mincExpanded = 1;
    return int(half * ns);
}

int applyDpnum(Case& c) {
    GridInput& g = c.grid;
    if (!g.dual || g.dpnum.empty() || g.dpnumApplied) return 0;
    const size_t n = size_t(g.size()), half = n / 2;
    if (g.nz % 2 || g.lgrExpanded) return 0;
    std::vector<char> sp(half, 0);
    int count = 0;
    for (size_t q = 0; q < half && q < g.dpnum.size(); ++q) if (g.dpnum[q] == 0) { sp[q] = 1; ++count; }
    if (!count) { g.dpnumApplied = 1; return 0; }
    if (g.actnum.size() != n) g.actnum.assign(n, 1);
    // fracture permeability is multiplied by the fracture porosity by the grid unless NODPPM: undo it for the moved matrix values
    const bool scaleK = !g.nodppm;
    if (g.permy.size() != n && g.permx.size() == n) g.permy = g.permx;
    if (g.permz.size() != n && g.permx.size() == n) g.permz = g.permx;
    for (const CellArrayRef& r : cellArrays(c)) {
        if (r.d == &g.dzmtrxv || r.d == &g.sigmav || r.d == &g.multmf || r.d == &g.lx || r.d == &g.ly || r.d == &g.lz) continue;
        if (r.d && r.d->size() == n)
            for (size_t q = 0; q < half; ++q) if (sp[q]) (*r.d)[q + half] = (*r.d)[q];
        if (r.i && r.i->size() == n)
            for (size_t q = 0; q < half; ++q) if (sp[q]) (*r.i)[q + half] = (*r.i)[q];
    }
    for (size_t q = 0; q < half; ++q) {
        if (!sp[q]) continue;
        g.actnum[q] = 0;
        if (scaleK && q < g.poro.size()) {
            const double phi = g.poro[q + half];
            for (auto* k : {&g.permx, &g.permy, &g.permz})
                if (k->size() == n) (*k)[q + half] = phi > 0 ? (*k)[q + half] / phi : 0.0;
        }
    }
    // MULTX/Y/Z (and the other multipliers) were copied; NTG of fracture cells is normally 1, the moved cell takes the matrix NTG
    g.dpnumApplied = 1;
    return count;
}

// ------------------------------------------------------------------------------------------------ C15 hydraulic fractures
namespace {
// fracture pieces from the processed host geometry (works for Cartesian and corner-point cells: the trace is clipped
// against the cell's mid-depth quadrilateral, the depth interval against the cell's mean top / bottom)
void hfPiecesOf(const GridInput& g, const Grid& G, std::vector<GridInput::HfPiece>& out) {
    out.clear();
    const int n = G.numCells(), nl = G.nx * G.ny;
    for (int f = 0; f < int(g.hfracs.size()); ++f) {
        const GridInput::HFrac& F = g.hfracs[size_t(f)];
        const double Lx = F.x1 - F.x0, Ly = F.y1 - F.y0, L = std::hypot(Lx, Ly);
        if (!(L > 0)) throw std::runtime_error("HFRAC " + F.name + ": zero-length trace");
        if (!(F.width > 0)) throw std::runtime_error("HFRAC " + F.name + ": width must be > 0");
        const double tx = Lx / L, ty = Ly / L, nxv = -ty, nyv = tx;
        const double sx0 = std::min(F.x0, F.x1), sx1 = std::max(F.x0, F.x1), sy0 = std::min(F.y0, F.y1), sy1 = std::max(F.y0, F.y1);
        for (int q = 0; q < n; ++q) {
            if (G.activeOf[size_t(q)] < 0) continue;
            if (G.nzMatrix > 0 && q / nl < G.nzMatrix) continue;   // dual porosity: hosts are fracture-system cells
            const auto& C = G.corners[size_t(q)];
            static const int ord[4] = {0, 1, 3, 2};
            double px[4], py[4];
            double bx0 = 1e300, bx1 = -1e300, by0 = 1e300, by1 = -1e300;
            for (int m = 0; m < 4; ++m) {
                const int id = ord[m];
                px[m] = 0.5 * (double(C[size_t(id)].x) + double(C[size_t(id + 4)].x));
                py[m] = 0.5 * (double(C[size_t(id)].y) + double(C[size_t(id + 4)].y));
                bx0 = std::min(bx0, px[m]); bx1 = std::max(bx1, px[m]); by0 = std::min(by0, py[m]); by1 = std::max(by1, py[m]);
            }
            if (sx1 < bx0 || sx0 > bx1 || sy1 < by0 || sy0 > by1) continue;
            double a2 = 0;
            for (int m = 0; m < 4; ++m) a2 += px[m] * py[(m + 1) % 4] - px[(m + 1) % 4] * py[m];
            if (std::fabs(a2) < 1e-12) continue;
            const double orient = a2 > 0 ? 1.0 : -1.0;
            double t0 = 0, t1 = L;
            bool in = true;
            for (int m = 0; m < 4 && in; ++m) {   // Cyrus-Beck clip against the convex quadrilateral
                const int b = (m + 1) % 4;
                const double ex = px[b] - px[m], ey = py[b] - py[m];
                const double nEx = -ey * orient, nEy = ex * orient;
                const double num = (F.x0 - px[m]) * nEx + (F.y0 - py[m]) * nEy, den = tx * nEx + ty * nEy;
                if (std::fabs(den) < 1e-14) { if (num < 0) in = false; continue; }
                const double r = -num / den;
                if (den > 0) t0 = std::max(t0, r); else t1 = std::min(t1, r);
            }
            if (!in || t1 - t0 <= 1e-9 * L) continue;
            const double top = 0.25 * (double(C[0].z) + C[1].z + C[2].z + C[3].z), bot = 0.25 * (double(C[4].z) + C[5].z + C[6].z + C[7].z);
            const double zt = F.zTop < 0 ? top : std::max(top, F.zTop), zb = F.zBot < 0 ? bot : std::min(bot, F.zBot);
            if (zb - zt <= 1e-9 * std::max(1.0, bot - top)) continue;
            // mean normal distance of the host from the fracture plane (bilinear 12 x 12 midpoint rule, area weighted)
            double dsum = 0, wsum = 0;
            for (int ia = 0; ia < 12; ++ia)
                for (int ib = 0; ib < 12; ++ib) {
                    const double u = (ia + 0.5) / 12, v = (ib + 0.5) / 12;
                    const double x = (1 - u) * (1 - v) * px[0] + u * (1 - v) * px[1] + u * v * px[2] + (1 - u) * v * px[3];
                    const double y = (1 - u) * (1 - v) * py[0] + u * (1 - v) * py[1] + u * v * py[2] + (1 - u) * v * py[3];
                    const double xu = (1 - v) * (px[1] - px[0]) + v * (px[2] - px[3]), yu = (1 - v) * (py[1] - py[0]) + v * (py[2] - py[3]);
                    const double xv = (1 - u) * (px[3] - px[0]) + u * (px[2] - px[1]), yv = (1 - u) * (py[3] - py[0]) + u * (py[2] - py[1]);
                    const double J = std::fabs(xu * yv - xv * yu);
                    dsum += J * std::fabs((x - F.x0) * nxv + (y - F.y0) * nyv); wsum += J;
                }
            GridInput::HfPiece P;
            P.frac = f; P.host = q; P.ta = t0; P.tb = t1;
            P.xa = F.x0 + tx * t0; P.ya = F.y0 + ty * t0; P.xb = F.x0 + tx * t1; P.yb = F.y0 + ty * t1;
            P.zt = zt; P.zb = zb;
            P.dAvg = std::max(wsum > 0 ? dsum / wsum : 0.0, 1e-3 * std::sqrt(std::fabs(0.5 * a2)) + 1e-9);
            const double kx = q < int(g.permx.size()) ? G.effPerm(g.permx, q) : 0.0;
            const double ky = q < int(g.permy.size()) && !std::isnan(g.permy[size_t(q)]) ? G.effPerm(g.permy, q) : kx;
            const double ntg = q < int(g.ntg.size()) ? g.ntg[size_t(q)] : 1.0;
            P.kn = (kx * nxv * nxv + ky * nyv * nyv) * ntg;
            out.push_back(P);
        }
    }
}
}  // namespace

int prepareHydraulicFractures(Case& c) {
    GridInput& g = c.grid;
    if (g.hfracs.empty() || g.hfExpanded) return int(g.hfExpanded ? g.hfPieces.size() : 0);
    expandLgrArrays(c);
    expandMincArrays(c);
    applyDpnum(c);
    Grid G;
    {
        Case c0 = c;
        c0.grid.hfracs.clear();
        G.build(c0);
    }
    hfPiecesOf(g, G, g.hfPieces);
    if (g.hfPieces.empty()) throw std::runtime_error("HFRAC: no hydraulic fracture intersects an active cell");
    const size_t nAll0 = size_t(G.numAllCells()), nH = g.hfPieces.size();
    for (const CellArrayRef& r : cellArrays(c)) {
        const size_t cur = r.d ? r.d->size() : r.i->size();
        if (cur != nAll0) continue;
        if (r.d) { r.d->reserve(nAll0 + nH); for (const auto& P : g.hfPieces) r.d->push_back((*r.d)[size_t(P.host)]); }
        else { r.i->reserve(nAll0 + nH); for (const auto& P : g.hfPieces) r.i->push_back((*r.i)[size_t(P.host)]); }
    }
    // fracture cells are always active and take their own SATNUM when given
    if (g.actnum.size() == nAll0 + nH) for (size_t k = 0; k < nH; ++k) g.actnum[nAll0 + k] = 1;
    bool ownSat = false;
    for (const auto& F : g.hfracs) if (F.satnum > 0) ownSat = true;
    if (ownSat && g.satnum.empty()) g.satnum.assign(nAll0 + nH, 1);
    if (g.satnum.size() == nAll0 + nH)
        for (size_t k = 0; k < nH; ++k) { const int s = g.hfracs[size_t(g.hfPieces[k].frac)].satnum; if (s > 0) g.satnum[nAll0 + k] = s; }
    g.hfExpanded = 1;
    return int(nH);
}

void Grid::buildHydraulicFractures(const Case& c, double cDarcy) {
    const GridInput& g = c.grid;
    const Units U = Units::of(c.units);
    const size_t base = size_t(hfBase()), nH = g.hfPieces.size(), nT = base + nH;
    activeOf.resize(nT, -1);
    topDepth.resize(nT, 0); cellDx.resize(nT, 0); cellDy.resize(nT, 0); cellDz.resize(nT, 0); bulkVol.resize(nT, 0);
    center.resize(nT); corners.resize(nT);
    if (!permScale.empty()) permScale.resize(nT, 1.0);
    hfCells.clear();
    const double lenM = c.units == UnitSystem::Field ? 0.3048 : 1.0, mD = 9.869233e-16;
    std::vector<int> act(nH, -1);
    for (size_t i = 0; i < nH; ++i) {
        const GridInput::HfPiece& P = g.hfPieces[i];
        const GridInput::HFrac& F = g.hfracs[size_t(P.frac)];
        const int ha = P.host < int(activeOf.size()) ? activeOf[size_t(P.host)] : -1;
        const size_t gq = base + i;
        const double L = std::hypot(P.xb - P.xa, P.yb - P.ya), h = P.zb - P.zt, w = F.width;
        const double kf = F.perm > 0 ? F.perm : (w * lenM) * (w * lenM) / 12.0 / mD;
        HfCell hc; hc.frac = P.frac; hc.host = P.host; hc.global = int(gq); hc.kf = kf; hc.w = w; hc.len = L; hc.h = h; hc.gamma = F.gamma;
        hfCells.push_back(hc);
        // geometry: a thin slab around the fracture plane (drawn with a visible minimum thickness)
        const double tx = (P.xb - P.xa) / std::max(L, 1e-12), ty = (P.yb - P.ya) / std::max(L, 1e-12);
        const double half = 0.5 * std::max(w, 0.02 * std::min(L, h));
        const double nx_ = -ty * half, ny_ = tx * half;
        auto V = [](double x, double y, double z) { return Vec3{float(x), float(y), float(z)}; };
        corners[gq] = {V(P.xa - nx_, P.ya - ny_, P.zt), V(P.xb - nx_, P.yb - ny_, P.zt), V(P.xa + nx_, P.ya + ny_, P.zt), V(P.xb + nx_, P.yb + ny_, P.zt),
                       V(P.xa - nx_, P.ya - ny_, P.zb), V(P.xb - nx_, P.yb - ny_, P.zb), V(P.xa + nx_, P.ya + ny_, P.zb), V(P.xb + nx_, P.yb + ny_, P.zb)};
        center[gq] = V(0.5 * (P.xa + P.xb), 0.5 * (P.ya + P.yb), 0.5 * (P.zt + P.zb));
        topDepth[gq] = P.zt; cellDx[gq] = L; cellDy[gq] = w; cellDz[gq] = h; bulkVol[gq] = L * w * h;
        if (ha < 0) continue;
        activeOf[gq] = numActive(); globalOf.push_back(int(gq));
        act[i] = activeOf[gq];
        poreVolRef.push_back(std::max(1e-6, F.poro) * L * w * h * U.poreVolume);
        depth.push_back(0.5 * (P.zt + P.zb)); dzCell.push_back(h);
        // matrix - fracture
        const double T = cDarcy * 2.0 * L * h * P.kn / P.dAvg;
        if (T > 0) { Connection cn; cn.a = ha; cn.b = act[i]; cn.trans = T; cn.dDepth = depth.back() - depth[size_t(ha)]; cn.dir = 'F'; conns.push_back(cn); }
    }
    // fracture - fracture: along the trace (shared end, overlapping depth), vertically (overlapping trace, shared depth),
    // and at intersections of different fractures inside the same host (Moinfar et al. 2014)
    auto ff = [&](int a, int b, double ta, double tb) {
        if (a < 0 || b < 0 || !(ta > 0) || !(tb > 0)) return;
        Connection cn; cn.a = a; cn.b = b; cn.trans = cDarcy * ta * tb / (ta + tb); cn.dDepth = depth[size_t(b)] - depth[size_t(a)]; cn.dir = 'F';
        conns.push_back(cn);
    };
    for (size_t u = 0; u < nH; ++u)
        for (size_t v = u + 1; v < nH; ++v) {
            const GridInput::HfPiece& A = g.hfPieces[u]; const GridInput::HfPiece& B = g.hfPieces[v];
            const HfCell& HA = hfCells[u]; const HfCell& HB = hfCells[v];
            const double La = A.tb - A.ta, Lb = B.tb - B.ta, tol = 1e-6 * (La + Lb);
            if (A.frac == B.frac) {
                const double zov = std::min(A.zb, B.zb) - std::max(A.zt, B.zt);
                const double tov = std::min(A.tb, B.tb) - std::max(A.ta, B.ta);
                const double ztol = 1e-6 * ((A.zb - A.zt) + (B.zb - B.zt));
                if (zov > ztol && (std::fabs(A.tb - B.ta) < tol || std::fabs(B.tb - A.ta) < tol))
                    ff(act[u], act[v], HA.kf * HA.w * zov / (0.5 * La), HB.kf * HB.w * zov / (0.5 * Lb));
                else if (tov > tol && (std::fabs(A.zb - B.zt) < ztol || std::fabs(B.zb - A.zt) < ztol))
                    ff(act[u], act[v], HA.kf * HA.w * tov / (0.5 * HA.h), HB.kf * HB.w * tov / (0.5 * HB.h));
            } else if (A.host == B.host) {
                const GridInput::HFrac& Fa = g.hfracs[size_t(A.frac)]; const GridInput::HFrac& Fb = g.hfracs[size_t(B.frac)];
                const double ax = Fa.x1 - Fa.x0, ay = Fa.y1 - Fa.y0, bx = Fb.x1 - Fb.x0, by = Fb.y1 - Fb.y0;
                const double den = ax * by - ay * bx;
                if (std::fabs(den) < 1e-14) continue;
                const double sa = ((Fb.x0 - Fa.x0) * by - (Fb.y0 - Fa.y0) * bx) / den, sb = ((Fb.x0 - Fa.x0) * ay - (Fb.y0 - Fa.y0) * ax) / den;
                const double tA = sa * std::hypot(ax, ay), tB = sb * std::hypot(bx, by);
                if (tA < A.ta || tA > A.tb || tB < B.ta || tB > B.tb) continue;
                const double a1 = tA - A.ta, a2 = A.tb - tA, b1 = tB - B.ta, b2 = B.tb - tB;
                const double dA = std::max(1e-6, (a1 * a1 + a2 * a2) / (2 * (a1 + a2))), dB = std::max(1e-6, (b1 * b1 + b2 * b2) / (2 * (b1 + b2)));
                const double h = std::min(A.zb, B.zb) - std::max(A.zt, B.zt);
                if (h > 0) ff(act[u], act[v], HA.kf * HA.w * h / dA, HB.kf * HB.w * h / dB);
            }
        }
}

std::string Grid::cellLabel(int g) const {
    char b[160];
    int l, i, j, k;
    if (localIJK(g, l, i, j, k)) {
        int I, J, K; ijk(g, I, J, K);
        std::snprintf(b, sizeof b, "%s (%d, %d, %d) in (%d, %d, %d)", lgr[size_t(l)].name.c_str(), i + 1, j + 1, k + 1, I + 1, J + 1, K + 1);
    } else if (isHf(g)) {
        ijk(g, i, j, k);
        const HfCell& h = hfCells[size_t(g - hfBase())];
        std::snprintf(b, sizeof b, "fracture %d piece in (%d, %d, %d)", h.frac + 1, i + 1, j + 1, k + 1);
    } else if (isMincShell(g)) {
        ijk(g, i, j, k);
        std::snprintf(b, sizeof b, "(%d, %d, %d) MINC shell %d", i + 1, j + 1, k + 1, mincShell[size_t(g - numCells() - numLgrCells())]);
    } else { ijk(g, i, j, k); std::snprintf(b, sizeof b, "(%d, %d, %d)", i + 1, j + 1, k + 1); }
    return b;
}

std::vector<int> Grid::localCellsOnFace(int host, char dir, bool minus) const {
    std::vector<int> out;
    if (!isHost(host)) return out;
    const int l = hostLgr[size_t(host)];
    const LgrInfo& L = lgr[size_t(l)];
    for (int k = 0; k < L.nz; ++k)
        for (int j = 0; j < L.ny; ++j)
            for (int i = 0; i < L.nx; ++i) {
                const int gl = localCell(l, i, j, k);
                if (lgrHost[size_t(gl - numCells())] != host) continue;
                int ni = i, nj = j, nk = k;
                const int s = minus ? -1 : 1;
                if (dir == 'X') ni += s; else if (dir == 'Y') nj += s; else nk += s;
                const bool edge = ni < 0 || nj < 0 || nk < 0 || ni >= L.nx || nj >= L.ny || nk >= L.nz ||
                                  lgrHost[size_t(localCell(l, ni, nj, nk) - numCells())] != host;
                if (edge) out.push_back(gl);
            }
    return out;
}

std::vector<int> Grid::hostWellColumn(int host, char dir) const {
    std::vector<int> out;
    if (!isHost(host)) return out;
    const int l = hostLgr[size_t(host)];
    const LgrInfo& L = lgr[size_t(l)];
    const Vec3 hc = center[size_t(host)];
    // the local cell of this host nearest the host centre fixes the two transverse indices
    int bi = 0, bj = 0, bk = 0; double bd = 1e300;
    for (int k = 0; k < L.nz; ++k)
        for (int j = 0; j < L.ny; ++j)
            for (int i = 0; i < L.nx; ++i) {
                const int gl = localCell(l, i, j, k);
                if (lgrHost[size_t(gl - numCells())] != host) continue;
                const Vec3 p = center[size_t(gl)];
                const double d = double(p.x - hc.x) * (p.x - hc.x) + double(p.y - hc.y) * (p.y - hc.y) + double(p.z - hc.z) * (p.z - hc.z);
                if (d < bd) { bd = d; bi = i; bj = j; bk = k; }
            }
    const char d = char(std::toupper((unsigned char)dir));
    for (int k = 0; k < L.nz; ++k)
        for (int j = 0; j < L.ny; ++j)
            for (int i = 0; i < L.nx; ++i) {
                if (d == 'X' ? (j != bj || k != bk) : d == 'Y' ? (i != bi || k != bk) : (i != bi || j != bj)) continue;
                const int gl = localCell(l, i, j, k);
                if (lgrHost[size_t(gl - numCells())] == host) out.push_back(gl);
            }
    return out;
}

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
    if (c.grid.dual && !c.grid.dpnum.empty() && !c.grid.dpnumApplied) {
        Case cp = c;   // single porosity cells (DPNUM = 0) folded into their fracture-layer cell
        applyDpnum(cp);
        build(cp);
        return;
    }
    if (c.grid.minc > 1 && c.grid.dual == 1 && !c.grid.mincExpanded) {
        Case cp = c;   // MINC: per-cell arrays also cover the inner matrix shells
        expandMincArrays(cp);
        build(cp);
        return;
    }
    if (!c.grid.hfracs.empty() && !c.grid.hfExpanded) {
        Case cp = c;   // C15 hydraulic fractures: per-cell arrays also cover the fracture cells
        prepareHydraulicFractures(cp);
        build(cp);
        return;
    }
    if (!c.grid.lgrs.empty() && !c.grid.lgrExpanded) {
        // local grid refinements: build from a copy whose per-cell arrays carry the local cells
        Case cp = c;
        expandLgrArrays(cp);
        build(cp);
        return;
    }
    const GridInput& g = c.grid;
    nx = g.nx; ny = g.ny; nz = g.nz;
    if (nx < 1 || ny < 1 || nz < 1) throw std::runtime_error("GRID: dimensions must be >= 1");
    const size_t n = size_t(nx) * ny * nz;
    size_t nL = 0;
    for (const LgrInput& L : g.lgrs) nL += size_t(L.size());
    const int nmShells = (g.minc > 1 && g.dual == 1 && nz % 2 == 0) ? g.minc : 0;
    const size_t nM = nmShells ? n / 2 * size_t(nmShells - 1) : 0;
    const size_t nH = (g.hfExpanded && !g.hfracs.empty()) ? g.hfPieces.size() : 0;
    const size_t nAll = n + nL + nM + nH;
    mincHost.clear(); mincShell.clear(); mincFrac.clear(); mincN = 0; hfCells.clear();
    lgr.clear(); lgrHost.clear(); hostLgr.clear();
    if (nL && g.dual) throw std::runtime_error("GRID: local grid refinements in dual porosity models are not modelled by the ResSim engine");
    const size_t nl = size_t(nx) * ny;
    const Units U = Units::of(c.units);
    cornerPoint = g.cornerPoint();

    auto poro = fillOr(g.poro, n, 0, "PORO", true, nAll);
    auto kx = fillOr(g.permx, n, 0, "PERMX", true, nAll);
    auto ky = g.permy.empty() ? kx : fillOr(g.permy, n, 0, "PERMY", false, nAll);
    auto kz = g.permz.empty() ? kx : fillOr(g.permz, n, 0, "PERMZ", false, nAll);
    auto ntg = fillOr(g.ntg, n, 1.0, "NTG", false, nAll);
    auto mx = fillOr(g.multx, n, 1.0, "MULTX", false, nAll);
    auto my = fillOr(g.multy, n, 1.0, "MULTY", false, nAll);
    auto mz = fillOr(g.multz, n, 1.0, "MULTZ", false, nAll);
    // C38 EDIT-section TRANX / TRANY / TRANZ multipliers: Cartesian neighbour connections only (folded into the multiplier)
    auto tranEdit = [&](const std::vector<double>& v, const char* nm) {   // LGR cells (beyond n) keep 1
        auto r = fillOr(v, n, 1.0, nm, false, nAll);
        if (r.size() < size_t(nAll)) r.resize(size_t(nAll), 1.0);
        return r;
    };
    const auto tx = tranEdit(g.tranMultX, "TRANX"), ty = tranEdit(g.tranMultY, "TRANY"), tz = tranEdit(g.tranMultZ, "TRANZ");
    if (poro.size() != nAll || kx.size() != nAll || ky.size() != nAll || kz.size() != nAll)
        throw std::runtime_error("GRID: PORO / PERM arrays do not cover the local grid refinements");

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
    std::vector<char> isAqu(n, 0);   // C34 numerical aquifer cells are always active (OPM)
    for (const NumAquCell& a : g.numAquCells)
        if (a.i >= 0 && a.i < nx && a.j >= 0 && a.j < ny && a.k >= 0 && a.k < nz) isAqu[size_t(g.index(a.i, a.j, a.k))] = 1;
    for (size_t q = 0; q < n; ++q) {
        bool act = g.actnum.empty() || g.actnum[q] != 0;
        if (isAqu[q]) { activeOf[q] = int(globalOf.size()); globalOf.push_back(int(q)); continue; }
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

    int na = numActive();
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
    // g1 / g2: geometric half-cell values A / d without permeability / NTG (conduction geometry, C17)
    auto add = [&](int q1, int q2, double t1, double t2, double mult, char dir, bool nnc, double g1 = 0, double g2 = 0) {
        const int a = activeOf[q1], b = activeOf[q2];
        if (a < 0 || b < 0 || t1 <= 0 || t2 <= 0 || mult <= 0) return;
        if (dual) {   // geometric neighbours never couple matrix to fracture; DUALPORO matrix blocks are isolated
            const bool m1 = q1 / int(nl) < nzh, m2 = q2 / int(nl) < nzh;
            if (m1 != m2 || (g.dual == 1 && m1)) return;
        }
        Connection cn;
        cn.a = a; cn.b = b; cn.trans = U.cDarcy * mult / (1.0 / t1 + 1.0 / t2);
        cn.dDepth = depth[b] - depth[a]; cn.dir = dir;
        cn.dTa = cn.trans * (1.0 / t1) / (1.0 / t1 + 1.0 / t2); cn.dTb = cn.trans * (1.0 / t2) / (1.0 / t1 + 1.0 / t2);
        cn.mq = q1;
        if (g1 > 0 && g2 > 0) cn.geo = mult / (1.0 / g1 + 1.0 / g2);
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
                                  kx[r] * ntg[r] * dyC[r] * dzC[r] / (0.5 * dxC[r]), mx[q] * tx[q], 'X', false,
                            dyC[q] * dzC[q] / (0.5 * dxC[q]), dyC[r] * dzC[r] / (0.5 * dxC[r]));
                    }
                    if (j + 1 < ny) {
                        const int r = g.index(i, j + 1, k);
                        add(q, r, ky[q] * ntg[q] * dxC[q] * dzC[q] / (0.5 * dyC[q]),
                                  ky[r] * ntg[r] * dxC[r] * dzC[r] / (0.5 * dyC[r]), my[q] * ty[q], 'Y', false,
                            dxC[q] * dzC[q] / (0.5 * dyC[q]), dxC[r] * dzC[r] / (0.5 * dyC[r]));
                    }
                    if (k + 1 < nz) {
                        const int r = g.index(i, j, k + 1);
                        add(q, r, kz[q] * dxC[q] * dyC[q] / (0.5 * dzC[q]),
                                  kz[r] * dxC[r] * dyC[r] / (0.5 * dzC[r]), mz[q] * tz[q], 'Z', false,
                            dxC[q] * dyC[q] / (0.5 * dzC[q]), dxC[r] * dyC[r] / (0.5 * dzC[r]));
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
                                (dir == 0 ? mx[q] : my[q]) * (k1 == k2 ? (dir == 0 ? tx[q] : ty[q]) : 1.0), dir == 0 ? 'X' : 'Y', k1 != k2, halfT(q, 1.0, 1.0, A, fc), halfT(r, 1.0, 1.0, A, fc));
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
                    add(q, r, halfT(q, kz[q], 1.0, A, fc), halfT(r, kz[r], 1.0, A, fc), mz[q] * tz[q], 'Z', false, halfT(q, 1.0, 1.0, A, fc), halfT(r, 1.0, 1.0, A, fc));
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
                    double gap = 0, mmin = mz[q]; int qmin = q;
                    bool ok = true; int r = -1;
                    for (int k2 = k + 1; k2 < kEnd; ++k2) {
                        const int c2 = g.index(i, j, k2);
                        if (activeOf[c2] >= 0) { r = c2; break; }
                        const bool thin = cellDz[c2] < g.pinch.threshold;
                        // GAP bridges cells removed by MINPV / zero pore volume; a cell switched off with ACTNUM = 0 that is
                        // thicker than the threshold is a barrier (ECLIPSE writes no pinch-out NNC across it - Norne layer 4)
                        const bool userInactive = !g.actnum.empty() && g.actnum[size_t(c2)] == 0;
                        if (!thin && (!g.pinch.gap || userInactive)) { ok = false; break; }
                        gap += std::max(0.0, cellDz[c2]);
                        if (mz[c2] < mmin) { mmin = mz[c2]; qmin = c2; }
                    }
                    if (!ok || r < 0 || r == g.index(i, j, k + 1) || gap > g.pinch.maxGap) continue;
                    const double t1 = halfZ(q, 5), t2 = halfZ(r, 4);
                    const double mult = g.pinch.multzAll ? mmin : mz[q];
                    if (t1 <= 0 || t2 <= 0 || mult <= 0) continue;
                    Connection cn; cn.a = activeOf[q]; cn.b = activeOf[r];
                    cn.trans = U.cDarcy * mult / (1.0 / t1 + 1.0 / t2);
                    cn.dDepth = depth[cn.b] - depth[cn.a]; cn.dir = 'Z';
                    cn.dTa = cn.trans * (1.0 / t1) / (1.0 / t1 + 1.0 / t2); cn.dTb = cn.trans * (1.0 / t2) / (1.0 / t1 + 1.0 / t2);
                    cn.mq = g.pinch.multzAll ? qmin : q;
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
            double sk = 0, sg = 0;   // sigma K and the geometric sigma (conduction)
            const double lX = cellVal(g.lx, qm, 0), lY = cellVal(g.ly, qm, 0), lZ = cellVal(g.lz, qm, 0);
            if (g.ltosigma && (lX > 0 || lY > 0 || lZ > 0)) {
                if (lX > 0) { sk += g.ltoF[0] * kx[qm] / (lX * lX); sg += g.ltoF[0] / (lX * lX); }
                if (lY > 0) { sk += g.ltoF[1] * ky[qm] / (lY * lY); sg += g.ltoF[1] / (lY * lY); }
                if (lZ > 0) { sk += g.ltoF[2] * kz[qm] / (lZ * lZ); sg += g.ltoF[2] / (lZ * lZ); }
            } else { sg = cellVal(g.sigmav, qm, g.sigma); sk = sg * kx[qm]; }
            if (nmShells) continue;   // MINC: shell transmissibilities below
            const double T = U.cDarcy * bulkVol[qm] * sk * cellVal(g.multmf, qm, 1.0);
            if (!(T > 0)) continue;
            Connection cn; cn.a = a; cn.b = b; cn.trans = T; cn.dDepth = depth[b] - depth[a]; cn.dir = 'M';
            cn.geo = bulkVol[qm] * sg * cellVal(g.multmf, qm, 1.0);
            conns.push_back(cn); ++numMF;
        }
        if (nmShells) buildMinc(g, U.cDarcy, n, nl, nzh, nmShells, kx, ky, kz);
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
    // MULTREGT on every grid connection (neighbours, fault NNCs, pinch-outs); the last matching record wins
    numMultRegT = 0;
    if (!g.multregt.empty()) {
        const int nlay = nx * ny;
        for (auto& cn : conns) {
            if (cn.dir == 'M' || cn.dir == 'N') continue;
            int ga = globalOf[size_t(cn.a)], gb = globalOf[size_t(cn.b)];
            if (ga > gb) std::swap(ga, gb);
            const int d = gb - ga;
            const int ia = ga % nx, ja = (ga / nx) % ny;
            const bool neighbour = (d == 1 && ia + 1 < nx) || (d == nx && ja + 1 < ny) || d == nlay;
            for (size_t r = g.multregt.size(); r-- > 0;) {
                const GridInput::MultRegT& m = g.multregt[r];
                const std::vector<int>& R = m.set == 'F' ? g.fluxnum : m.set == 'O' ? g.opernum : g.multnum;
                if (R.size() < size_t(n)) continue;
                const int ra = R[size_t(ga)], rb = R[size_t(gb)];
                if (ra == rb) continue;
                const bool match = ((m.r1 < 0 || ra == m.r1) && (m.r2 < 0 || rb == m.r2)) || ((m.r1 < 0 || rb == m.r1) && (m.r2 < 0 || ra == m.r2));
                if (!match) continue;
                if ((m.nnc == 1 && neighbour) || (m.nnc == 2 && !neighbour)) continue;
                if (m.dirs.find(cn.dir) == std::string::npos) continue;
                cn.trans *= m.m; cn.dTa *= m.m; cn.dTb *= m.m; cn.geo *= m.m; ++numMultRegT;
                break;
            }
        }
        conns.erase(std::remove_if(conns.begin(), conns.end(), [](const Connection& c) { return c.trans <= 0; }), conns.end());
    }
    // C34 EDITNNC / EDITNNCR (EDIT section): multiply / replace the transmissibility of the non-neighbour connection between
    // two cells (fault juxtaposition, pinch-out or explicit NNC; Cartesian neighbours are skipped)
    numEditNnc = 0;
    if (!g.editnnc.empty() || !g.editnncr.empty()) {
        std::map<std::pair<int, int>, std::pair<double, double>> ed;   // (a, b) sorted -> (multiplier, replacement or -1)
        auto key = [&](const NNCInput& e, bool& ok) {
            ok = !(e.i1 < 0 || e.i1 >= nx || e.j1 < 0 || e.j1 >= ny || e.k1 < 0 || e.k1 >= nz || e.i2 < 0 || e.i2 >= nx || e.j2 < 0 ||
                   e.j2 >= ny || e.k2 < 0 || e.k2 >= nz);
            if (!ok) return std::make_pair(-1, -1);
            int a = activeOf[g.index(e.i1, e.j1, e.k1)], b = activeOf[g.index(e.i2, e.j2, e.k2)];
            if (a > b) std::swap(a, b);
            ok = a >= 0 && b >= 0 && a != b;
            return std::make_pair(a, b);
        };
        for (const NNCInput& e : g.editnnc) { bool ok; auto k = key(e, ok); if (!ok) continue; auto it = ed.emplace(k, std::make_pair(1.0, -1.0)).first; it->second.first *= e.trans; }
        for (const NNCInput& e : g.editnncr) { bool ok; auto k = key(e, ok); if (!ok) continue; auto it = ed.emplace(k, std::make_pair(1.0, -1.0)).first; it->second.second = e.trans; it->second.first = 1.0; }
        const int nlay = nx * ny;
        for (auto& cn : conns) {
            if (cn.dir == 'M') continue;
            auto it = ed.find(std::make_pair(std::min(cn.a, cn.b), std::max(cn.a, cn.b)));
            if (it == ed.end()) continue;
            int ga = globalOf[size_t(cn.a)], gb = globalOf[size_t(cn.b)];
            if (ga > gb) std::swap(ga, gb);
            const int d = gb - ga, ia = ga % nx, ja = (ga / nx) % ny;
            // ECLIPSE: EDITNNC addresses non-neighbour connections only (OPM would also scale a Cartesian neighbour; model2
            // 9_EDITNNC is unchanged either way - all its entries are fault juxtapositions)
            if (cn.dir != 'N' && ((d == 1 && ia + 1 < nx) || (d == nx && ja + 1 < ny) || d == nlay)) continue;   // Cartesian neighbours
            double f = it->second.first;
            if (it->second.second >= 0) f = cn.trans > 0 ? it->second.second / cn.trans : 0.0;
            cn.trans *= f; cn.dTa *= f; cn.dTb *= f; cn.geo *= f; ++numEditNnc;
        }
        conns.erase(std::remove_if(conns.begin(), conns.end(), [](const Connection& c) { return c.trans <= 0; }), conns.end());
    }
    // ---- C34 numerical aquifers (AQUNUM / AQUCON; OPM SingleNumericalAquifer): the aquifer cells lose their grid
    // connections and get pore volume A L phi and the aquifer depth; consecutive cells of one aquifer are joined with
    // T = 1 / (1/Tc1 + 1/Tc2), Tc = 2 k A / L, and the AQUCON faces connect to the first cell of the aquifer with
    // T = Tcon Tcell / (Tcon + Tcell) x mult, Tcon = 2 MULT k A_face NTG / d (TRANS_OPTION 1: Tcell = 2 k1 A_face / L1)
    numAquP0.clear(); numAquConns = 0;
    if (!g.numAquCells.empty()) {
        std::map<int, std::vector<const NumAquCell*>> aq;
        for (const NumAquCell& a : g.numAquCells) aq[a.id].push_back(&a);
        std::vector<char> aqA(size_t(numActive()), 0);
        for (const NumAquCell& a : g.numAquCells) { const int q = g.index(a.i, a.j, a.k); if (activeOf[size_t(q)] >= 0) aqA[size_t(activeOf[size_t(q)])] = 1; }
        conns.erase(std::remove_if(conns.begin(), conns.end(), [&](const Connection& cn) { return aqA[size_t(cn.a)] || aqA[size_t(cn.b)]; }), conns.end());
        auto tc = [](const NumAquCell* a) { return a->length > 0 ? 2.0 * a->perm * a->area / a->length : 0.0; };
        auto addC = [&](int a, int b, double T, char dir) {
            if (a < 0 || b < 0 || a == b || !(T > 0)) return;
            Connection cn; cn.a = a; cn.b = b; cn.trans = T; cn.dDepth = depth[size_t(b)] - depth[size_t(a)]; cn.dir = dir;
            cn.dTa = cn.dTb = 0.5 * T; conns.push_back(cn);
        };
        for (auto& [id, cells] : aq) {
            for (const NumAquCell* c : cells) {
                const int a = activeOf[size_t(g.index(c->i, c->j, c->k))];
                if (a < 0) continue;
                poreVolRef[size_t(a)] = c->area * c->length * c->poro * U.poreVolume;
                depth[size_t(a)] = c->depth;
                numAquP0.push_back({a, c->p0});
            }
            for (size_t m = 0; m + 1 < cells.size(); ++m) {
                const double t1 = tc(cells[m]), t2 = tc(cells[m + 1]);
                if (t1 > 0 && t2 > 0)
                    addC(activeOf[size_t(g.index(cells[m]->i, cells[m]->j, cells[m]->k))], activeOf[size_t(g.index(cells[m + 1]->i, cells[m + 1]->j, cells[m + 1]->k))],
                         U.cDarcy / (1.0 / t1 + 1.0 / t2), 'A');
            }
            const NumAquCell* c1 = cells.front();
            const int a1 = activeOf[size_t(g.index(c1->i, c1->j, c1->k))];
            if (a1 < 0) continue;
            for (const NumAquCon& cc : g.numAquCons) {
                if (cc.id != id) continue;
                for (int k = std::max(0, cc.k1); k <= std::min(nz - 1, cc.k2); ++k)
                    for (int j = std::max(0, cc.j1); j <= std::min(ny - 1, cc.j2); ++j)
                        for (int i = std::max(0, cc.i1); i <= std::min(nx - 1, cc.i2); ++i) {
                            const size_t q = size_t(g.index(i, j, k));
                            const int b = activeOf[q];
                            if (b < 0 || aqA[size_t(b)]) continue;
                            const int d = cc.face == 'X' ? 0 : cc.face == 'Y' ? 1 : 2;
                            int ni = i, nj = j, nk = k; (d == 0 ? ni : d == 1 ? nj : nk) += cc.minus ? -1 : 1;
                            if (!cc.allowInternal && ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
                                const int nb = activeOf[size_t(g.index(ni, nj, nk))];
                                if (nb >= 0 && !aqA[size_t(nb)]) continue;   // the face is internal to the reservoir
                            }
                            const double dx = cellDx[q], dy = cellDy[q], dz = cellDz[q];
                            const double A = d == 0 ? dy * dz : d == 1 ? dx * dz : dx * dy, L = d == 0 ? dx : d == 1 ? dy : dz;
                            const double k_ = (d == 0 ? kx : d == 1 ? ky : kz)[q];
                            const double mult = cc.minus ? 1.0 : (d == 0 ? mx : d == 1 ? my : mz)[q];
                            const double tcon = L > 0 ? 2.0 * mult * k_ * A * ntg[q] / L : 0.0;
                            const double tcell = cc.transOpt == 0 ? tc(c1) : (c1->length > 0 ? 2.0 * c1->perm * A / c1->length : 0.0);
                            if (tcon > 0 && tcell > 0) { addC(a1, b, U.cDarcy * tcon * tcell / (tcon + tcell) * cc.mult, 'A'); ++numAquConns; }
                        }
            }
        }
    }
    fieldDerivOk = !dual && g.numAquCells.empty();
    if (nL) { buildLgrImpl(*this, c, C, ctr, poro, kx, ky, kz, ntg, mx, my, mz); fieldDerivOk = false; }
    if (nH) { buildHydraulicFractures(c, U.cDarcy); fieldDerivOk = false; }
    // C38 NONNC: only Cartesian neighbours connect (X: i+-1 same j k, Y: j+-1 same i k); the cross-layer connections a fault
    // throw creates are dropped (ECLIPSE / OPM: no non-neighbour connections are generated)
    if (g.nonnc) {
        const size_t before = conns.size();
        const int na0 = int(globalOf.size());
        conns.erase(std::remove_if(conns.begin(), conns.end(), [&](const Connection& cn) {
            if (cn.a >= na0 || cn.b >= na0 || (cn.dir != 'X' && cn.dir != 'Y')) return false;
            const int qa = globalOf[size_t(cn.a)], qb = globalOf[size_t(cn.b)];
            const int ia = qa % nx, ja = (qa / nx) % ny, ka = qa / (nx * ny), ib = qb % nx, jb = (qb / nx) % ny, kb = qb / (nx * ny);
            if (ka != kb) return true;
            return cn.dir == 'X' ? !(ja == jb && std::abs(ia - ib) == 1) : !(ia == ib && std::abs(ja - jb) == 1);
        }), conns.end());
        numNonnc = int(before - conns.size());
    }
    // merge duplicate cell pairs (keeps the sparsity pattern clean)
    for (auto& cn : conns) if (cn.a > cn.b) { std::swap(cn.a, cn.b); std::swap(cn.dTa, cn.dTb); cn.dDepth = -cn.dDepth; }
    std::sort(conns.begin(), conns.end(), [](const Connection& x, const Connection& y) { return x.a != y.a ? x.a < y.a : x.b < y.b; });
    size_t m = 0;
    for (size_t i = 0; i < conns.size(); ++i) {
        if (m > 0 && conns[m - 1].a == conns[i].a && conns[m - 1].b == conns[i].b) {
            conns[m - 1].trans += conns[i].trans; conns[m - 1].dTa += conns[i].dTa; conns[m - 1].dTb += conns[i].dTb;
            conns[m - 1].geo = conns[m - 1].geo > 0 && conns[i].geo > 0 ? conns[m - 1].geo + conns[i].geo : 0.0;   // unknown part: fallback T / k
            if (conns[m - 1].mq != conns[i].mq) conns[m - 1].mq = -2;   // merged faces with different multiplier cells
        }
        else conns[m++] = conns[i];
    }
    conns.resize(m);

    conns.shrink_to_fit();   // C37 memory: drop the push_back slack (up to 2x of 64 B per connection)
    // ---- adjacency
    na = numActive();   // local grid refinements change the active cell count
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
