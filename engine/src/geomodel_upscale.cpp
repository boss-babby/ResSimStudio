// C32 geomodelling: upscaling to a simulation grid (zone-respecting coarsening, PV / bulk weighting, arithmetic /
// harmonic / geometric / flow-based permeability), conversion to a simulation case and GRDECL export
#include "ressim/geomodel.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

namespace ressim {
namespace geo {

namespace {
// chunks of 'f' consecutive indices over [0, n)
std::vector<std::pair<int, int>> chunks(int n, int f) {
    std::vector<std::pair<int, int>> c;
    for (int a = 0; a < n; a += std::max(1, f)) c.push_back({a, std::min(n, a + std::max(1, f)) - 1});
    return c;
}
}  // namespace

// TPFA single-phase steady state on a box: p = 1 on the inlet face, 0 on the outlet face, no-flow elsewhere;
// k_eff = Q L / (A dp). Conjugate gradients with a Jacobi preconditioner (SPD).
double flowBasedPerm(int nx, int ny, int nz, const std::vector<double>& kIn, double dx, double dy, double dz, int dir) {
    const int n = nx * ny * nz;
    if (n <= 0 || int(kIn.size()) < n) return 0;
    std::vector<double> k(static_cast<size_t>(n));
    for (int c = 0; c < n; ++c) k[size_t(c)] = std::max(kIn[size_t(c)], 1e-12);
    const int N[3] = {nx, ny, nz};
    const double h[3] = {dx, dy, dz};
    auto id = [&](int i, int j, int l) { return i + nx * (j + ny * l); };
    auto harm = [](double a, double b) { return 2 * a * b / (a + b); };
    // face transmissibilities per direction d: area / distance
    auto area = [&](int d) { return d == 0 ? dy * dz : d == 1 ? dx * dz : dx * dy; };
    std::vector<double> diag(size_t(n), 0.0), rhs(size_t(n), 0.0);
    struct Off { int a, b; double t; };
    std::vector<Off> off;
    for (int l = 0; l < nz; ++l) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const int c = id(i, j, l), ijk[3] = {i, j, l};
        for (int d = 0; d < 3; ++d) {
            if (ijk[d] + 1 < N[d]) {
                const int o = d == 0 ? id(i + 1, j, l) : d == 1 ? id(i, j + 1, l) : id(i, j, l + 1);
                const double t = harm(k[size_t(c)], k[size_t(o)]) * area(d) / h[d];
                off.push_back({c, o, t}); diag[size_t(c)] += t; diag[size_t(o)] += t;
            }
        }
        // Dirichlet faces along 'dir' (half-cell transmissibility)
        if (ijk[dir] == 0) { const double t = 2 * k[size_t(c)] * area(dir) / h[dir]; diag[size_t(c)] += t; rhs[size_t(c)] += t * 1.0; }
        if (ijk[dir] == N[dir] - 1) { const double t = 2 * k[size_t(c)] * area(dir) / h[dir]; diag[size_t(c)] += t; }
    }
    auto mul = [&](const std::vector<double>& x, std::vector<double>& y) {
        for (int c = 0; c < n; ++c) y[size_t(c)] = diag[size_t(c)] * x[size_t(c)];
        for (const auto& o : off) { y[size_t(o.a)] -= o.t * x[size_t(o.b)]; y[size_t(o.b)] -= o.t * x[size_t(o.a)]; }
    };
    const size_t nn = static_cast<size_t>(n);
    std::vector<double> p(nn, 0.5), r(nn, 0.0), z(nn, 0.0), q(nn, 0.0), Ap(nn, 0.0);
    // initial guess: linear in the flow direction
    for (int l = 0; l < nz; ++l) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { const int ijk[3] = {i, j, l}; p[size_t(id(i, j, l))] = 1.0 - (ijk[dir] + 0.5) / N[dir]; }
    mul(p, Ap);
    double bn = 0;
    for (int c = 0; c < n; ++c) { r[size_t(c)] = rhs[size_t(c)] - Ap[size_t(c)]; bn += rhs[size_t(c)] * rhs[size_t(c)]; }
    bn = std::sqrt(std::max(bn, 1e-300));
    for (int c = 0; c < n; ++c) { z[size_t(c)] = r[size_t(c)] / diag[size_t(c)]; q[size_t(c)] = z[size_t(c)]; }
    double rz = 0; for (int c = 0; c < n; ++c) rz += r[size_t(c)] * z[size_t(c)];
    for (int it = 0; it < 20 * n + 100; ++it) {
        double rn = 0; for (int c = 0; c < n; ++c) rn += r[size_t(c)] * r[size_t(c)];
        if (std::sqrt(rn) < 1e-12 * bn) break;
        mul(q, Ap);
        double pAp = 0; for (int c = 0; c < n; ++c) pAp += q[size_t(c)] * Ap[size_t(c)];
        if (pAp <= 0) break;
        const double a = rz / pAp;
        for (int c = 0; c < n; ++c) { p[size_t(c)] += a * q[size_t(c)]; r[size_t(c)] -= a * Ap[size_t(c)]; }
        double rz2 = 0; for (int c = 0; c < n; ++c) { z[size_t(c)] = r[size_t(c)] / diag[size_t(c)]; rz2 += r[size_t(c)] * z[size_t(c)]; }
        const double b = rz2 / rz; rz = rz2;
        for (int c = 0; c < n; ++c) q[size_t(c)] = z[size_t(c)] + b * q[size_t(c)];
    }
    // outflow through the outlet face
    double Q = 0;
    for (int l = 0; l < nz; ++l) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const int ijk[3] = {i, j, l};
        if (ijk[dir] != N[dir] - 1) continue;
        const int c = id(i, j, l);
        Q += 2 * k[size_t(c)] * area(dir) / h[dir] * p[size_t(c)];
    }
    const double L = N[dir] * h[dir], A = (dir == 0 ? ny * dy * nz * dz : dir == 1 ? nx * dx * nz * dz : nx * dx * ny * dy);
    return Q * L / A;   // dp = 1
}

Upscaled upscale(const GeoGrid& F, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<double>& kx,
                 const std::vector<double>& ky, const std::vector<double>& kz, const UpscaleSpec& s) {
    Upscaled U;
    const int nx = F.grid.nx, ny = F.grid.ny, nz = F.grid.nz;
    const auto ci = chunks(nx, s.fx), cj = chunks(ny, s.fy);
    // layer groups inside each zone (coarse layers never straddle a zone boundary)
    std::vector<std::pair<int, int>> ck;
    std::vector<int> ckZone;
    for (int k = 0; k < nz;) {
        const int z = k < int(F.layerZone.size()) ? F.layerZone[size_t(k)] : 0;
        int e = k;
        while (e + 1 < nz && (e + 1 >= int(F.layerZone.size()) || F.layerZone[size_t(e + 1)] == z) && e + 1 - k < std::max(1, s.fz)) ++e;
        ck.push_back({k, e}); ckZone.push_back(z); k = e + 1;
    }
    const int NX = int(ci.size()), NY = int(cj.size()), NZ = int(ck.size());
    GeoGrid& G = U.g;
    G.spec = F.spec; G.spec.nx = NX; G.spec.ny = NY;
    // non-uniform last chunk: the coarse lattice keeps dx * fx (cornerXY); pillars of the last chunk use the fine edge
    G.spec.dx = F.spec.dx * std::max(1, s.fx); G.spec.dy = F.spec.dy * std::max(1, s.fy);
    G.layerZone = ckZone;
    GridInput& g = G.grid;
    g = GridInput(); g.nx = NX; g.ny = NY; g.nz = NZ;
    // COORD from the fine pillars at the chunk boundaries
    const int fpx = nx + 1;
    g.coord.resize(size_t(6) * size_t(NX + 1) * size_t(NY + 1));
    for (int J = 0; J <= NY; ++J) for (int I = 0; I <= NX; ++I) {
        const int fi = I < NX ? ci[size_t(I)].first : nx, fj = J < NY ? cj[size_t(J)].first : ny;
        for (int q = 0; q < 6; ++q) g.coord[size_t(6 * (I + (NX + 1) * J) + q)] = F.grid.coord[size_t(6 * (fi + fpx * fj) + q)];
    }
    g.zcorn.assign(size_t(8) * size_t(NX) * size_t(NY) * size_t(NZ), 0.0);
    auto fz = [&](int i, int j, int k, int f, int a, int b) { return F.grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j + b) * 2 * nx + 2 * i + a)]; };
    for (int K = 0; K < NZ; ++K) for (int J = 0; J < NY; ++J) for (int I = 0; I < NX; ++I)
        for (int f = 0; f < 2; ++f) for (int b = 0; b < 2; ++b) for (int a = 0; a < 2; ++a) {
            const int i = a ? ci[size_t(I)].second : ci[size_t(I)].first, j = b ? cj[size_t(J)].second : cj[size_t(J)].first, k = f ? ck[size_t(K)].second : ck[size_t(K)].first;
            g.zcorn[size_t(((2 * K + f) * 2 * NY + 2 * J + b) * 2 * NX + 2 * I + a)] = fz(i, j, k, f, a, b);
        }
    // geometry (lattice with the coarse spacing; last chunk volumes are approximate when nx % fx != 0)
    computeGeometry(G);
    const size_t Nc = size_t(NX) * size_t(NY) * size_t(NZ);
    U.poro.assign(Nc, 0); U.ntg.assign(Nc, 0); U.kx.assign(Nc, 0); U.ky.assign(Nc, 0); U.kz.assign(Nc, 0); U.facies.clear();
    g.actnum.assign(Nc, 0);
    // fault blocks of the coarse columns: block of the first fine column
    G.block.assign(size_t(NX) * size_t(NY), 0); G.nBlocks = F.nBlocks;
    for (int J = 0; J < NY; ++J) for (int I = 0; I < NX; ++I) G.block[size_t(I + NX * J)] = F.block.empty() ? 0 : F.block[size_t(ci[size_t(I)].first + nx * cj[size_t(J)].first)];
    for (int K = 0; K < NZ; ++K) for (int J = 0; J < NY; ++J) for (int I = 0; I < NX; ++I) {
        const size_t C = size_t(I + NX * (J + NY * K));
        double V = 0, Vn = 0, Vnp = 0, ka[3] = {0, 0, 0}, kh[3] = {0, 0, 0}, kg[3] = {0, 0, 0}, Va = 0;
        const int bx = ci[size_t(I)].second - ci[size_t(I)].first + 1, by = cj[size_t(J)].second - cj[size_t(J)].first + 1, bz = ck[size_t(K)].second - ck[size_t(K)].first + 1;
        std::vector<double> box[3];
        for (auto& v : box) v.assign(size_t(bx * by * bz), 0.0);
        double dzs = 0; int nact = 0;
        for (int k = ck[size_t(K)].first; k <= ck[size_t(K)].second; ++k)
            for (int j = cj[size_t(J)].first; j <= cj[size_t(J)].second; ++j)
                for (int i = ci[size_t(I)].first; i <= ci[size_t(I)].second; ++i) {
                    const size_t c = size_t(i + nx * (j + ny * k));
                    const bool act = F.grid.actnum.empty() || F.grid.actnum[c];
                    const size_t bi = size_t((i - ci[size_t(I)].first) + bx * ((j - cj[size_t(J)].first) + by * (k - ck[size_t(K)].first)));
                    if (!act) continue;
                    const double v = F.bulk[c], nt = c < ntg.size() ? ntg[c] : 1.0, ph = c < poro.size() ? poro[c] : 0.0;
                    V += v; Vn += v * nt; Vnp += v * nt * ph;
                    const double kk[3] = {c < kx.size() ? kx[c] : 0.0, c < ky.size() ? ky[c] : 0.0, c < kz.size() ? kz[c] : 0.0};
                    for (int d = 0; d < 3; ++d) {
                        box[d][bi] = kk[d];
                        ka[d] += v * kk[d]; kh[d] += kk[d] > 0 ? v / kk[d] : 1e30; kg[d] += v * std::log(std::max(kk[d], 1e-12));
                    }
                    Va += v; dzs += F.bot[c] - F.top[c]; ++nact;
                }
        if (V <= 0 || nact == 0) continue;
        g.actnum[C] = 1;
        // NTG bulk-weighted; porosity such that the coarse pore volume (coarse corner geometry) equals the fine sum exactly
        U.ntg[C] = Vn / V;
        const double Vc = G.bulk[C] > 0 ? G.bulk[C] : V;
        U.poro[C] = Vn > 0 ? Vnp / (U.ntg[C] * Vc) : 0.0;
        double* out[3] = {&U.kx[C], &U.ky[C], &U.kz[C]};
        for (int d = 0; d < 3; ++d) {
            if (s.permMethod == 0) *out[d] = ka[d] / Va;
            else if (s.permMethod == 1) *out[d] = Va / kh[d];
            else if (s.permMethod == 2) *out[d] = std::exp(kg[d] / Va);
            else *out[d] = flowBasedPerm(bx, by, bz, box[d], F.spec.dx, F.spec.dy, std::max(1e-6, dzs / nact), d);
        }
    }
    // faults: fine faces on coarse column boundaries
    for (const auto& ff : F.grid.faults) {
        FaultInput fc; fc.name = ff.name; fc.mult = ff.mult;
        std::set<std::tuple<int, int, char>> seen;
        for (const auto& f : ff.faces)
            for (int j = f.j1; j <= f.j2; ++j) for (int i = f.i1; i <= f.i2; ++i) {
                int I = -1, J = -1;
                for (int q = 0; q < NX; ++q) if (i >= ci[size_t(q)].first && i <= ci[size_t(q)].second) I = q;
                for (int q = 0; q < NY; ++q) if (j >= cj[size_t(q)].first && j <= cj[size_t(q)].second) J = q;
                if (I < 0 || J < 0) continue;
                const bool onB = f.face == 'X' ? i == ci[size_t(I)].second : f.face == 'Y' ? j == cj[size_t(J)].second : false;
                if (!onB || !seen.insert({I, J, f.face}).second) continue;
                FaultFace c; c.i1 = c.i2 = I; c.j1 = c.j2 = J; c.k1 = 0; c.k2 = NZ - 1; c.face = f.face; c.minus = f.minus;
                fc.faces.push_back(c);
            }
        if (!fc.faces.empty()) g.faults.push_back(fc);
    }
    return U;
}

void toCase(Case& c, const GeoGrid& g, const std::vector<double>& poro, const std::vector<double>& ntg, const std::vector<double>& kx,
            const std::vector<double>& ky, const std::vector<double>& kz, const std::vector<GeoWell>& wells) {
    const size_t n = g.bulk.size();
    GridInput gi = g.grid;
    gi.dx.clear(); gi.dy.clear(); gi.dz.clear(); gi.tops.clear();
    gi.poro = poro; gi.ntg = ntg; gi.permx = kx; gi.permy = ky; gi.permz = kz;
    gi.poro.resize(n, 0.0); gi.ntg.resize(n, 1.0); gi.permx.resize(n, 0.0); gi.permy.resize(n, 0.0); gi.permz.resize(n, 0.0);
    // carry the remaining case-level settings (regions etc. are re-sized by the caller when needed)
    gi.satnum.clear(); gi.pvtnum.clear(); gi.eqlnum.clear(); gi.fipnum.clear();
    c.grid = gi;
    const int nx = g.grid.nx, ny = g.grid.ny;
    for (const auto& w : wells) {
        const double L = w.totalMd();
        if (L <= 0) continue;
        std::vector<int> cells; std::vector<char> dirs;
        const double ds = std::max(0.1, std::min(g.spec.dx, g.spec.dy) / 10);
        for (double md = 0.5 * ds; md < L; md += ds) {
            const P3 p = w.at(md);
            const int cc = g.cellOf(p[0], p[1], p[2]);
            if (cc < 0 || (!g.grid.actnum.empty() && !g.grid.actnum[size_t(cc)])) continue;
            const P3 a = w.at(std::max(0.0, md - ds)), b = w.at(std::min(L, md + ds));
            double fa, fb, ga, gb; g.columnOf(a[0], a[1], fa, ga); g.columnOf(b[0], b[1], fb, gb);
            const double ex = std::fabs(fb - fa) * g.spec.dx, ey = std::fabs(gb - ga) * g.spec.dy, ez = std::fabs(b[2] - a[2]);
            const char d = ez >= std::max(ex, ey) ? 'Z' : ex >= ey ? 'X' : 'Y';
            if (cells.empty() || cells.back() != cc) { cells.push_back(cc); dirs.push_back(d); }
        }
        if (cells.empty()) continue;
        WellInput wi;
        for (auto& e : c.wells) if (e.name == w.name) wi = e;
        wi.name = w.name.size() > 8 ? w.name.substr(0, 8) : w.name;
        double fi, fj; g.columnOf(w.x, w.y, fi, fj);
        wi.i = std::max(0, std::min(nx - 1, int(fi))); wi.j = std::max(0, std::min(ny - 1, int(fj)));
        wi.completions.clear();
        std::set<int> seen;
        for (size_t q = 0; q < cells.size(); ++q) {
            if (!seen.insert(cells[q]).second) continue;
            Completion cp; cp.i = cells[q] % nx; cp.j = (cells[q] / nx) % ny; cp.k1 = cp.k2 = cells[q] / (nx * ny); cp.dir = dirs[q];
            wi.completions.push_back(cp);
        }
        bool replaced = false;
        for (auto& e : c.wells) if (e.name == wi.name) { e = wi; replaced = true; }
        if (!replaced) c.wells.push_back(wi);
    }
}

namespace {
void writeArr(FILE* f, const char* kw, const std::vector<double>& v, int prec) {
    std::fprintf(f, "%s\n", kw);
    int col = 0;
    for (size_t k = 0; k < v.size();) {
        size_t r = k + 1;
        while (r < v.size() && v[r] == v[k]) ++r;
        const size_t cnt = r - k;
        char b[64];
        if (cnt > 1) std::snprintf(b, sizeof b, "%zu*%.*g", cnt, prec, v[k]); else std::snprintf(b, sizeof b, "%.*g", prec, v[k]);
        std::fprintf(f, " %s", b);
        if (++col == 6) { std::fprintf(f, "\n"); col = 0; }
        k = r;
    }
    std::fprintf(f, "%s/\n\n", col ? "\n" : "");
}
}  // namespace

bool writeGrdecl(const std::string& path, const GeoGrid& g, const std::map<std::string, const std::vector<double>*>& props) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    const GridInput& gi = g.grid;
    std::fprintf(f, "-- ResGeo corner-point grid export (%d x %d x %d, %d fault blocks)\n\n", gi.nx, gi.ny, gi.nz, g.nBlocks);
    std::fprintf(f, "SPECGRID\n %d %d %d 1 F /\n\n", gi.nx, gi.ny, gi.nz);
    writeArr(f, "COORD", gi.coord, 10);
    writeArr(f, "ZCORN", gi.zcorn, 9);
    { std::vector<double> a(gi.actnum.begin(), gi.actnum.end()); writeArr(f, "ACTNUM", a, 3); }
    for (const auto& kv : props) if (kv.second && kv.second->size() == size_t(gi.nx) * size_t(gi.ny) * size_t(gi.nz)) writeArr(f, kv.first.c_str(), *kv.second, 6);
    if (!gi.faults.empty()) {
        std::fprintf(f, "FAULTS\n");
        for (const auto& F : gi.faults)
            for (const auto& x : F.faces)
                std::fprintf(f, " '%s' %d %d %d %d %d %d %c%s /\n", F.name.c_str(), x.i1 + 1, x.i2 + 1, x.j1 + 1, x.j2 + 1, x.k1 + 1, x.k2 + 1, x.face, x.minus ? "-" : "");
        std::fprintf(f, "/\n\nMULTFLT\n");
        for (const auto& F : gi.faults) std::fprintf(f, " '%s' %g /\n", F.name.c_str(), F.mult);
        std::fprintf(f, "/\n");
    }
    std::fclose(f);
    return true;
}

}  // namespace geo
}  // namespace ressim
