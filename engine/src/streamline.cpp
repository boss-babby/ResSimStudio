#include "ressim/streamline.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include "ressim/simulator.h"

namespace ressim {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

struct FaceLink { int nb; double out; };   // neighbour and outflux (>0 leaves this cell)

struct CellFaces {
    double U[3][2] = {{0, 0}, {0, 0}, {0, 0}};   // signed flux in +axis direction at face 0 / face 1, divided by PV
    std::vector<FaceLink> links[6];              // face = 2*axis + side
};

Vec3 mapPoint(const std::array<Vec3, 8>& c, double x, double y, double z) {
    auto L = [](const Vec3& a, const Vec3& b, double t) { return Vec3{float(a.x + (b.x - a.x) * t), float(a.y + (b.y - a.y) * t), float(a.z + (b.z - a.z) * t)}; };
    const Vec3 a = L(c[0], c[1], x), b = L(c[2], c[3], x), d = L(c[4], c[5], x), e = L(c[6], c[7], x);
    return L(L(a, b, y), L(d, e, y), z);
}

// time to leave [0,1] along one axis with linear velocity v(s) = u0 + (u1-u0) s, starting at s with v
double exitTime(double s, double u0, double u1, int& side) {
    const double v = u0 + (u1 - u0) * s, A = u1 - u0;
    side = -1;
    if (v > 1e-30) {
        if (u1 <= 0) return kInf;
        side = 1;
        return std::fabs(A) < 1e-12 * std::fabs(v) ? (1 - s) / v : std::log(u1 / v) / A;
    }
    if (v < -1e-30) {
        if (u0 >= 0) return kInf;
        side = 0;
        return std::fabs(A) < 1e-12 * std::fabs(v) ? s / -v : std::log(u0 / v) / A;
    }
    return kInf;
}

double advance(double s, double u0, double u1, double t) {
    const double v = u0 + (u1 - u0) * s, A = u1 - u0;
    if (std::fabs(A * t) < 1e-9) return s + v * t;
    return ((v * std::exp(A * t)) - u0) / A;
}

} // namespace

StreamlineResult traceStreamlines(const Simulator& sim, const StreamlineOptions& opt) {
    StreamlineResult R;
    const Grid& g = sim.grid();
    const int na = g.numActive();
    std::vector<double> qw, qo, qg;
    sim.connectionFluxes(qw, qo, qg);
    const auto& conns = g.conns;
    // ---- face fluxes (total reservoir volume rate) ----
    std::vector<CellFaces> F(na);
    auto faceOf = [&](int a, int b, char dir) {
        int ia, ja, ka, ib, jb, kb;
        const int ga = g.globalOf[a], gb = g.globalOf[b];
        if ((g.isLocal(ga) || g.isLocal(gb)) && (dir == 'X' || dir == 'Y' || dir == 'Z')) {
            // local grid refinement: side from the cell geometry (axis of cell a, centre of b)
            const auto& cr = g.corners[size_t(ga)];
            const int ax = dir == 'X' ? 1 : dir == 'Y' ? 2 : 4;
            float v[3] = {0, 0, 0};
            for (int m = 0; m < 8; ++m) { const float s = (m & ax) ? 0.25f : -0.25f; v[0] += s * cr[m].x; v[1] += s * cr[m].y; v[2] += s * cr[m].z; }
            const auto ca = g.center[size_t(ga)], cb = g.center[size_t(gb)];
            const float d = v[0] * (cb.x - ca.x) + v[1] * (cb.y - ca.y) + v[2] * (cb.z - ca.z);
            const int axis = dir == 'X' ? 0 : dir == 'Y' ? 1 : 2;
            return 2 * axis + (d > 0 ? 1 : 0);
        }
        g.ijk(ga, ia, ja, ka); g.ijk(gb, ib, jb, kb);
        int axis;
        if (dir == 'X') axis = 0; else if (dir == 'Y') axis = 1; else if (dir == 'Z') axis = 2;
        else {
            const int d[3] = {std::abs(ib - ia), std::abs(jb - ja), std::abs(kb - ka)};
            axis = d[0] >= d[1] && d[0] >= d[2] ? 0 : (d[1] >= d[2] ? 1 : 2);
        }
        const int da = axis == 0 ? ib - ia : axis == 1 ? jb - ja : kb - ka;
        return 2 * axis + (da > 0 ? 1 : 0);
    };
    std::vector<double> mfIn(na, 0.0);   // dual porosity: matrix -> fracture flux, a source (or sink) of the fracture cell
    for (size_t n = 0; n < conns.size(); ++n) {
        const auto& cn = conns[n];
        double q = qw[n] + qo[n] + qg[n];         // + from a to b
        if (opt.connFlux && opt.connFlux->size() == conns.size() && std::isfinite((*opt.connFlux)[n])) q = (*opt.connFlux)[n];
        if (cn.dir == 'M') {                       // streamlines are traced in the fracture system only (as FrontSim)
            const bool bFrac = g.isFracture(g.globalOf[cn.b]);
            mfIn[bFrac ? cn.b : cn.a] += bFrac ? q : -q;
            continue;
        }
        const int fa = faceOf(cn.a, cn.b, cn.dir), fb = faceOf(cn.b, cn.a, cn.dir);
        F[cn.a].links[fa].push_back({cn.b, q});
        F[cn.b].links[fb].push_back({cn.a, -q});
    }
    std::vector<double> pv(na);
    for (int a = 0; a < na; ++a) {
        pv[a] = std::max(1e-6, g.poreVolRef[a]);
        for (int f = 0; f < 6; ++f) {
            double out = 0;
            for (auto& l : F[a].links[f]) out += l.out;
            const int axis = f / 2, side = f % 2;
            F[a].U[axis][side] = (side == 1 ? out : -out) / pv[a];
        }
    }
    // ---- sources and sinks ----
    const auto wr = sim.wellCellRates();
    const auto& names = sim.wellNames();
    std::vector<double> cellSrc(na, 0.0);
    std::vector<int> cellSinkWell(na, -1), cellSrcId(na, -1);
    std::map<int, int> srcIndex, sinkIndex;   // key = well id (aquifers: -1-id)
    auto nameOf = [&](int w) { return w >= 0 ? (w < int(names.size()) ? names[w] : "W" + std::to_string(w)) : "AQUIFER " + std::to_string(-1 - w); };
    for (auto& c : wr) {
        if (c.rate < 0) {   // injection / aquifer influx
            if (!srcIndex.count(c.well)) { srcIndex[c.well] = int(R.sources.size()); R.sources.push_back(nameOf(c.well)); }
            cellSrc[c.cell] += -c.rate; cellSrcId[c.cell] = srcIndex[c.well];
        } else if (c.rate > 0) {
            if (!sinkIndex.count(c.well)) { sinkIndex[c.well] = int(R.sinks.size()); R.sinks.push_back(nameOf(c.well)); }
            cellSinkWell[c.cell] = sinkIndex[c.well];
        }
    }
    // matrix expulsion feeds the fracture streamlines (source "MATRIX")
    for (int a = 0; a < na; ++a) {
        if (mfIn[a] > 0 && cellSrc[a] <= 0) {
            if (!srcIndex.count(-1000000)) { srcIndex[-1000000] = int(R.sources.size()); R.sources.push_back("MATRIX"); }
            cellSrc[a] += mfIn[a]; cellSrcId[a] = srcIndex[-1000000];
        }   // imbibing fracture cells (flux into the matrix) are passed through, not sinks
    }
    R.pairFlux.assign(R.sources.size(), std::vector<double>(R.sinks.size() + 1, 0.0));   // last column: not reaching a producer
    double totalSrc = 0;
    for (int a = 0; a < na; ++a) totalSrc += cellSrc[a];
    R.cellTof.assign(na, -1.f); R.cellSink.assign(na, -1); R.cellSource.assign(na, -1);
    if (totalSrc <= 0) return R;
    std::vector<double> tofSum(na, 0.0), fluxSum(na, 0.0);
    std::vector<std::vector<double>> sinkFlux(na), srcFlux(na);
    const int maxCells = opt.maxCells > 0 ? opt.maxCells : 3 * na + 10;
    const int ppc = std::max(1, opt.pointsPerCell);

    // ---- seeding: per source cell, on its outflow faces, proportional to flux ----
    for (int a = 0; a < na; ++a) {
        if (cellSrc[a] <= 0) continue;
        double outTot = 0;
        for (int f = 0; f < 6; ++f) for (auto& l : F[a].links[f]) if (l.out > 0) outTot += l.out;
        if (outTot <= 0) continue;
        const int nCell = std::max(1, int(std::lround(opt.lines * cellSrc[a] / totalSrc)));
        double kept = 0;   // flux on faces that receive at least one seed (the rest is redistributed)
        for (int f = 0; f < 6; ++f) {
            double outF = 0;
            for (auto& l : F[a].links[f]) if (l.out > 0) outF += l.out;
            if (outF > 0 && std::lround(nCell * outF / outTot) > 0) kept += outF;
        }
        if (kept <= 0) continue;
        for (int f = 0; f < 6; ++f) {
            double outF = 0;
            for (auto& l : F[a].links[f]) if (l.out > 0) outF += l.out;
            if (outF <= 0) continue;
            const int nF = int(std::lround(nCell * outF / outTot));
            if (nF <= 0) continue;
            const int m = std::max(1, int(std::ceil(std::sqrt(double(nF)))));
            const double perLine = cellSrc[a] * (outF / kept) / (m * m);
            const int axis = f / 2, side = f % 2;
            for (int s1 = 0; s1 < m; ++s1)
                for (int s2 = 0; s2 < m; ++s2) {
                    double X[3];
                    const int ax1 = (axis + 1) % 3, ax2 = (axis + 2) % 3;
                    X[axis] = side; X[ax1] = (s1 + 0.5) / m; X[ax2] = (s2 + 0.5) / m;
                    Streamline L;
                    L.source = cellSrcId[a]; L.flux = perLine;
                    int cell = a;
                    double tof = 0;
                    L.pts.push_back(mapPoint(g.corners[g.globalOf[cell]], X[0], X[1], X[2])); L.tof.push_back(0.f); L.cells.push_back(cell);
                    // leave the source cell through face f: hop to the neighbour
                    int curFace = f;
                    bool first = true;
                    for (int it = 0; it < maxCells; ++it) {
                        if (!first) {
                            if (cellSinkWell[cell] >= 0) { L.sink = cellSinkWell[cell]; break; }
                            // Pollock step in this cell
                            const CellFaces& C = F[cell];
                            double tmin = kInf; int axMin = -1, sdMin = -1;
                            for (int d = 0; d < 3; ++d) {
                                int sd; const double t = exitTime(X[d], C.U[d][0], C.U[d][1], sd);
                                if (t < tmin) { tmin = t; axMin = d; sdMin = sd; }
                            }
                            if (!std::isfinite(tmin) || axMin < 0) break;   // stagnation
                            const auto& cr = g.corners[g.globalOf[cell]];
                            for (int sub = 1; sub <= ppc; ++sub) {
                                const double t = tmin * sub / ppc;
                                double Y[3];
                                for (int d = 0; d < 3; ++d) Y[d] = std::clamp(advance(X[d], C.U[d][0], C.U[d][1], t), 0.0, 1.0);
                                if (sub == ppc) { for (int d = 0; d < 3; ++d) X[d] = Y[d]; X[axMin] = sdMin; }
                                L.pts.push_back(mapPoint(cr, Y[0], Y[1], Y[2])); L.tof.push_back(float(tof + t));
                            }
                            const double t0 = tof; tof += tmin;
                            tofSum[cell] += L.flux * 0.5 * (t0 + tof); fluxSum[cell] += L.flux;
                            curFace = 2 * axMin + sdMin;
                        }
                        first = false;
                        // choose neighbour across curFace
                        const auto& links = F[cell].links[curFace];
                        int nb = -1; double best = 0;
                        if (links.size() == 1) { if (links[0].out > 0) nb = links[0].nb; }
                        else {
                            const Vec3 P = L.pts.back();
                            for (auto& l : links) {
                                if (l.out <= 0) continue;
                                // several neighbours across one face (fault throw, LGR host/child boundary): prefer the one
                                // whose bounding box contains the exit point, in depth first, then laterally
                                const auto& cb = g.corners[g.globalOf[l.nb]];
                                float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
                                for (auto& v : cb) {
                                    lo[0] = std::min(lo[0], v.x); hi[0] = std::max(hi[0], v.x); lo[1] = std::min(lo[1], v.y); hi[1] = std::max(hi[1], v.y);
                                    lo[2] = std::min(lo[2], v.z); hi[2] = std::max(hi[2], v.z);
                                }
                                auto in = [&](int d, double x) { const double tol = 1e-3 * (hi[d] - lo[d]) + 1e-6; return x >= lo[d] - tol && x <= hi[d] + tol; };
                                const double score = (in(2, P.z) ? 2e30 : 0) + (in(0, P.x) && in(1, P.y) ? 1e30 : 0) + l.out;
                                if (score > best) { best = score; nb = l.nb; }
                            }
                        }
                        if (nb < 0) break;
                        // entry coordinates in the neighbour
                        const int axis2 = curFace / 2, side2 = curFace % 2;
                        int ia, ja, ka, ib, jb, kb;
                        g.ijk(g.globalOf[cell], ia, ja, ka); g.ijk(g.globalOf[nb], ib, jb, kb);
                        X[axis2] = 1 - side2;
                        if (axis2 != 2 && kb != ka) {   // across a fault throw: map depth into the neighbour column
                            const auto& cb = g.corners[g.globalOf[nb]];
                            const Vec3 P = L.pts.back();
                            const Vec3 top = mapPoint(cb, X[0], X[1], 0), bot = mapPoint(cb, X[0], X[1], 1);
                            X[2] = std::clamp(double((P.z - top.z) / std::max(1e-6f, bot.z - top.z)), 0.0, 1.0);
                        }
                        cell = nb;
                        L.cells.push_back(cell);
                        if (cellSrc[cell] > 0 && it > 0 && cellSinkWell[cell] < 0) { /* passes another injector: continue */ }
                    }
                    if (L.sink < 0) ++R.stalled;
                    if (L.source >= 0) R.pairFlux[L.source][L.sink >= 0 ? L.sink : int(R.sinks.size())] += L.flux;
                    for (int c : L.cells) {
                        if (L.sink >= 0) { auto& v = sinkFlux[c]; if (v.size() < R.sinks.size()) v.resize(R.sinks.size(), 0); v[L.sink] += L.flux; }
                        if (L.source >= 0) { auto& v = srcFlux[c]; if (v.size() < R.sources.size()) v.resize(R.sources.size(), 0); v[L.source] += L.flux; }
                    }
                    R.maxTof = std::max(R.maxTof, double(L.tof.back()));
                    R.lines.push_back(std::move(L));
                }
        }
    }
    for (int a = 0; a < na; ++a) {
        if (fluxSum[a] > 0) R.cellTof[a] = float(tofSum[a] / fluxSum[a]);
        if (!sinkFlux[a].empty()) R.cellSink[a] = int(std::max_element(sinkFlux[a].begin(), sinkFlux[a].end()) - sinkFlux[a].begin());
        if (!srcFlux[a].empty()) R.cellSource[a] = int(std::max_element(srcFlux[a].begin(), srcFlux[a].end()) - srcFlux[a].begin());
    }
    return R;
}

} // namespace ressim
