#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "ressim/tools.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>

namespace ressim {

double hallYarboroughZ(double ppr, double tpr) {
    const double t = 1.0 / tpr;
    const double A = 0.06125 * t * std::exp(-1.2 * (1 - t) * (1 - t));
    const double B = t * (14.76 - 9.76 * t + 4.58 * t * t);
    const double C = t * (90.7 - 242.2 * t + 42.4 * t * t);
    const double D = 2.18 + 2.82 * t;
    double y = 0.0125 * ppr * t * std::exp(-1.2 * (1 - t) * (1 - t));
    if (y <= 0) return 1.0;
    for (int it = 0; it < 100; ++it) {
        const double om = 1 - y;
        const double f = -A * ppr + (y + y * y + y * y * y - y * y * y * y) / (om * om * om) - B * y * y + C * std::pow(y, D);
        const double df = (1 + 4 * y + 4 * y * y - 4 * y * y * y + y * y * y * y) / (om * om * om * om) - 2 * B * y +
                          C * D * std::pow(y, D - 1);
        double yn = y - f / df;
        if (yn <= 0) yn = 0.5 * y;
        if (yn >= 1) yn = 0.5 * (y + 1);
        if (std::fabs(yn - y) < 1e-12) { y = yn; break; }
        y = yn;
    }
    return A * ppr / y;
}

void makeCorrelationPVT(const CorrelationInput& in, UnitSystem units, PropsInput& out, bool deadOil) {
    const double api = in.api, gg = in.gasGravity, T = in.temperatureF, TR = T + 459.67;
    if (api <= 5 || gg <= 0.5 || T <= 40 || in.pbMax <= 20 || in.pMax <= in.pbMax || in.nRows < 3)
        throw std::runtime_error("Correlation input out of range");
    const double go = 141.5 / (131.5 + api);
    auto rsStanding = [&](double p) {   // scf/stb
        return gg * std::pow((p / 18.2 + 1.4) * std::pow(10.0, 0.0125 * api - 0.00091 * T), 1.2048);
    };
    auto boStanding = [&](double rs) { return 0.9759 + 0.00012 * std::pow(rs * std::sqrt(gg / go) + 1.25 * T, 1.2); };
    const double muod = std::pow(10.0, std::pow(T, -1.163) * std::exp(6.9824 - 0.04658 * api)) - 1.0;
    auto muLive = [&](double rs) {
        const double A = 10.715 * std::pow(rs + 100, -0.515), B = 5.44 * std::pow(rs + 150, -0.338);
        return A * std::pow(muod, B);
    };
    // gas
    const double tpc = 169.2 + 349.5 * gg - 74.0 * gg * gg, ppc = 756.8 - 131.0 * gg - 3.6 * gg * gg;
    const double M = 28.97 * gg;
    auto zOf = [&](double p) { return hallYarboroughZ(p / ppc, TR / tpc); };
    auto bgOf = [&](double p) { return 5.035 * zOf(p) * TR / p; };       // rb/Mscf
    auto mugOf = [&](double p) {
        const double z = zOf(p);
        const double K = (9.4 + 0.02 * M) * std::pow(TR, 1.5) / (209 + 19 * M + TR);
        const double X = 3.5 + 986.0 / TR + 0.01 * M, Y = 2.4 - 0.2 * X;
        const double rho = 1.4935e-3 * p * M / (z * TR);
        return 1e-4 * K * std::exp(X * std::pow(rho, Y));
    };
    const bool metric = units == UnitSystem::Metric;
    const double fp = metric ? 1.0 / 14.5038 : 1.0;                   // psia -> bar
    const double frs = metric ? 0.178108 : 1e-3;                      // scf/stb -> sm3/sm3 or Mscf/stb
    const double fbg = metric ? 5.614583e-3 : 1.0;                    // rb/Mscf -> rm3/sm3
    const double fd = metric ? 16.0185 : 1.0;                         // lb/ft3 -> kg/m3

    out.density.oil = 62.37 * go * fd;
    out.density.gas = 0.0764 * gg * fd;
    out.pvto.clear(); out.pvdo.clear(); out.pvdg.clear();
    const int n = in.nRows;
    const double p0 = 14.7;
    if (deadOil) {
        const double rsb = rsStanding(in.pbMax);
        for (int i = 0; i < n; ++i) {
            const double p = p0 + (in.pMax - p0) * i / (n - 1);
            const double co = 1e-5;
            out.pvdo.push_back({p * fp, boStanding(rsb) * std::exp(-co * (p - in.pbMax)), muLive(rsb)});
        }
    } else {
        for (int i = 0; i < n; ++i) {
            const double pb = p0 + (in.pbMax - p0) * i / (n - 1);
            const double rs = std::max(rsStanding(pb), 0.5);
            PVTORecord r{rs * frs, pb * fp, boStanding(rs), muLive(rs), {}};
            if (i == n - 1 || i == n - 2) {        // undersaturated branches at the two highest Rs
                const double co = std::max(1e-6, (-1433 + 5 * rs + 17.2 * T - 1180 * gg + 12.61 * api) / (1e5 * 0.5 * (pb + in.pMax)));
                for (int k = 1; k <= 3; ++k) {
                    const double p = pb + (in.pMax - pb) * k / 3.0;
                    const double m = 2.6 * std::pow(p, 1.187) * std::exp(-11.513 - 8.98e-5 * p);
                    r.under.push_back({p * fp, r.bo * std::exp(-co * (p - pb)), r.mu * std::pow(p / pb, m)});
                }
            }
            out.pvto.push_back(r);
        }
    }
    for (int i = 0; i < n; ++i) {
        const double p = p0 + (in.pMax - p0) * std::pow(double(i) / (n - 1), 1.3);
        out.pvdg.push_back({p * fp, bgOf(p) * fbg, mugOf(p)});
    }
}

void makeCoreySWOF(const CoreyInput& c, std::vector<SatRow>& out) {
    out.clear();
    const int n = std::max(4, c.rows);
    const double mob = 1.0 - c.swc - c.sorw;
    if (mob <= 0) throw std::runtime_error("Corey: Swc + Sorw must be < 1");
    for (int i = 0; i <= n; ++i) {
        const double sw = c.swc + (1.0 - c.swc) * i / n;
        const double se = std::max(0.0, std::min(1.0, (sw - c.swc) / mob));
        const double so = std::max(0.0, std::min(1.0, (1.0 - sw - c.sorw) / mob));
        const double pc = c.pcEntry * std::pow(1.0 - std::min(1.0, (sw - c.swc) / (1.0 - c.swc)), c.pcExponent);
        out.push_back({sw, c.krwMax * std::pow(se, c.nw), c.kroMax * std::pow(so, c.now), pc});
    }
}

void makeCoreySGOF(const CoreyInput& c, std::vector<SatRow>& out) {
    out.clear();
    const int n = std::max(4, c.rows);
    const double sgMax = 1.0 - c.swc;
    const double mobG = 1.0 - c.swc - c.sgc - c.sorg, mobO = 1.0 - c.swc - c.sorg;
    if (mobG <= 0 || mobO <= 0) throw std::runtime_error("Corey: Swc + Sgc + Sorg must be < 1");
    for (int i = 0; i <= n; ++i) {
        const double sg = sgMax * i / n;
        const double sge = std::max(0.0, std::min(1.0, (sg - c.sgc) / mobG));
        const double so = std::max(0.0, std::min(1.0, (1.0 - c.swc - sg - c.sorg) / mobO));
        out.push_back({sg, c.krgMax * std::pow(sge, c.ng), c.kroMax * std::pow(so, c.nog), 0.0});
    }
}

void generateGeometry(const GeomInput& in, GridInput& g) {
    if (in.nx < 1 || in.ny < 1 || in.nz < 1) throw std::runtime_error("Grid dimensions must be >= 1");
    g.nx = in.nx; g.ny = in.ny; g.nz = in.nz;
    const int n = g.size();
    g.dx.assign(n, in.dx); g.dy.assign(n, in.dy); g.dz.resize(n);
    for (int k = 0; k < g.nz; ++k) {
        const double dz = in.dz.empty() ? 10.0 : in.dz[std::min<size_t>(k, in.dz.size() - 1)];
        for (int q = 0; q < g.nx * g.ny; ++q) g.dz[k * g.nx * g.ny + q] = dz;
    }
    g.tops.resize(size_t(g.nx) * g.ny);
    const double tx = std::tan(in.dipX * M_PI / 180.0), ty = std::tan(in.dipY * M_PI / 180.0);
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i)
            g.tops[j * g.nx + i] = in.top + tx * (i + 0.5) * in.dx + ty * (j + 0.5) * in.dy;
    if (g.actnum.size() != size_t(n)) g.actnum.clear();
}

void generateRock(const RockGenInput& in, GridInput& g) {
    const int nx = g.nx, ny = g.ny, nz = g.nz, n = g.size();
    auto layerVal = [&](const std::vector<double>& v, int k, double def) {
        if (v.empty()) return def;
        return v[std::min<size_t>(k, v.size() - 1)];
    };
    g.permx.assign(n, 0); g.permy.assign(n, 0); g.permz.assign(n, 0); g.poro.assign(n, 0); g.ntg.assign(n, 1);
    std::mt19937 rng(in.seed);
    std::normal_distribution<double> N01(0.0, 1.0);
    const int rx = std::max(0, int(std::ceil(2.0 * in.corrX))), ry = std::max(0, int(std::ceil(2.0 * in.corrY)));
    std::vector<double> white((nx + 2 * rx) * (ny + 2 * ry)), field(nx * ny), tmp((nx + 2 * rx) * ny);
    // Gaussian kernel weights
    auto kernel = [](int r, double L) {
        std::vector<double> w(2 * r + 1);
        double s = 0;
        for (int d = -r; d <= r; ++d) { w[d + r] = L > 0 ? std::exp(-0.5 * (d * d) / (0.25 * L * L + 1e-12)) : (d == 0); s += w[d + r] * w[d + r]; }
        for (double& x : w) x /= std::sqrt(s);      // unit-variance output for unit-variance white noise (separable)
        return w;
    };
    const auto wx = kernel(rx, in.corrX), wy = kernel(ry, in.corrY);
    const int W = nx + 2 * rx;
    for (int k = 0; k < nz; ++k) {
        for (double& v : white) v = N01(rng);
        // filter in y then x (separable)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < W; ++i) {
                double s = 0;
                for (int d = -ry; d <= ry; ++d) s += wy[d + ry] * white[(j + ry + d) * W + i];
                tmp[j * W + i] = s;
            }
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                double s = 0;
                for (int d = -rx; d <= rx; ++d) s += wx[d + rx] * tmp[j * W + i + rx + d];
                field[j * nx + i] = s;
            }
        const double kL = layerVal(in.layerPerm, k, 100), phiL = layerVal(in.layerPoro, k, 0.2);
        const double kvkh = layerVal(in.layerKvKh, k, 0.1), ntg = layerVal(in.layerNtg, k, 1.0);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int q = g.index(i, j, k);
                // mean-preserving lognormal: E[k] = kL
                const double kk = kL * std::exp(in.sigmaLn * field[j * nx + i] - 0.5 * in.sigmaLn * in.sigmaLn);
                g.permx[q] = kk; g.permy[q] = kk; g.permz[q] = kk * kvkh;
                g.poro[q] = std::min(0.45, std::max(0.01, phiL * std::pow(kk / kL, in.poroPermExp)));
                g.ntg[q] = ntg;
            }
    }
}

} // namespace ressim

// ---------------------------------------------------------------- corner-point helpers
namespace ressim {

void convertToCornerPoint(GridInput& g) {
    if (g.cornerPoint()) return;
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const size_t n = size_t(g.size()), nl = size_t(nx) * ny;
    if (g.dx.size() != n || g.dy.size() != n || g.dz.size() != n) throw std::runtime_error("convertToCornerPoint: DX/DY/DZ missing");
    std::vector<double> xs(nx + 1, 0.0), ys(ny + 1, 0.0);
    for (int i = 0; i < nx; ++i) xs[i + 1] = xs[i] + g.dx[g.index(i, 0, 0)];
    for (int j = 0; j < ny; ++j) ys[j + 1] = ys[j] + g.dy[g.index(0, j, 0)];
    std::vector<double> top(n);
    for (size_t q = 0; q < nl; ++q) top[q] = g.tops.empty() ? 0.0 : g.tops[q];
    for (int k = 1; k < nz; ++k)
        for (size_t q = 0; q < nl; ++q) {
            const size_t gi = k * nl + q;
            top[gi] = (g.tops.size() == n) ? g.tops[gi] : top[gi - nl] + g.dz[gi - nl];
        }
    double zmin = 1e300, zmax = -1e300;
    for (size_t q = 0; q < n; ++q) { zmin = std::min(zmin, top[q]); zmax = std::max(zmax, top[q] + g.dz[q]); }
    g.coord.resize(size_t(6) * (nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6 * (size_t(j) * (nx + 1) + i)];
            p[0] = xs[i]; p[1] = ys[j]; p[2] = zmin; p[3] = xs[i]; p[4] = ys[j]; p[5] = zmax;
        }
    g.zcorn.resize(8 * n);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int q = g.index(i, j, k);
                for (int kk = 0; kk < 2; ++kk)
                    for (int jj = 0; jj < 2; ++jj)
                        for (int ii = 0; ii < 2; ++ii)
                            g.zcorn[(size_t(2 * k + kk) * 2 * ny + (2 * j + jj)) * 2 * nx + 2 * i + ii] = top[q] + kk * g.dz[q];
            }
    g.dx.clear(); g.dy.clear(); g.dz.clear(); g.tops.clear();
}

void addFaultThrow(GridInput& g, int pos, double throwDepth, double mult, const std::string& name, bool alongJ) {
    if (!g.cornerPoint()) convertToCornerPoint(g);
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    if (pos <= 0 || pos >= (alongJ ? ny : nx)) throw std::runtime_error("addFaultThrow: fault position must be inside the grid");
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                if ((alongJ ? j : i) < pos) continue;
                for (int kk = 0; kk < 2; ++kk)
                    for (int jj = 0; jj < 2; ++jj)
                        for (int ii = 0; ii < 2; ++ii)
                            g.zcorn[(size_t(2 * k + kk) * 2 * ny + (2 * j + jj)) * 2 * nx + 2 * i + ii] += throwDepth;
            }
    double zmin = 1e300, zmax = -1e300;
    for (double z : g.zcorn) { zmin = std::min(zmin, z); zmax = std::max(zmax, z); }
    for (size_t p = 0; p * 6 < g.coord.size(); ++p) {
        g.coord[6 * p + 2] = std::min(g.coord[6 * p + 2], zmin);
        g.coord[6 * p + 5] = std::max(g.coord[6 * p + 5], zmax);
    }
    FaultInput f; f.name = name; f.mult = mult;
    FaultFace ff;
    if (!alongJ) { ff.i1 = ff.i2 = pos - 1; ff.j1 = 0; ff.j2 = ny - 1; ff.face = 'X'; }
    else { ff.j1 = ff.j2 = pos - 1; ff.i1 = 0; ff.i2 = nx - 1; ff.face = 'Y'; }
    ff.k1 = 0; ff.k2 = nz - 1;
    f.faces.push_back(ff);
    g.faults.push_back(f);
}


namespace {
template <class V> void forCellArrays(GridInput& g, V&& f) {
    for (auto* v : {&g.dx, &g.dy, &g.dz, &g.poro, &g.permx, &g.permy, &g.permz, &g.ntg, &g.multx, &g.multy, &g.multz,
                    &g.multpv, &g.swatinit, &g.swl, &g.swcr, &g.swu, &g.sgl, &g.sgcr, &g.sgu, &g.sowcr, &g.sogcr, &g.spoly,
                    &g.sigmav, &g.lx, &g.ly, &g.lz, &g.multmf, &g.dzmtrxv})
        f(*v);
    for (auto* v : {&g.actnum, &g.satnum, &g.pvtnum, &g.eqlnum, &g.fipnum, &g.plmixnum}) f(*v);
}
}  // namespace

void makeDualPorosity(Case& c, int mode, double fracPoro, double fracPerm, double blockL, double blockH, bool fractureRegion) {
    GridInput& g = c.grid;
    if (g.dual) { g.dual = mode; return; }
    const size_t n = size_t(g.size()), nl = size_t(g.nx) * g.ny;
    forCellArrays(g, [&](auto& v) { if (v.size() == n) { auto cp = v; v.insert(v.end(), cp.begin(), cp.end()); } });
    if (g.tops.size() == n) { auto cp = g.tops; g.tops.insert(g.tops.end(), cp.begin(), cp.end()); }
    if (!g.zcorn.empty()) { auto cp = g.zcorn; g.zcorn.insert(g.zcorn.end(), cp.begin(), cp.end()); }
    const int oldNz = g.nz;
    g.nz *= 2;
    const size_t n2 = size_t(g.size());
    if (g.permy.size() != n2) g.permy = g.permx;
    if (g.permz.size() != n2) g.permz = g.permx;
    for (size_t q = n; q < n2; ++q) { g.poro[q] = fracPoro; g.permx[q] = g.permy[q] = g.permz[q] = fracPerm; }
    if (!g.ntg.empty()) for (size_t q = n; q < n2; ++q) g.ntg[q] = 1.0;
    if (!g.swatinit.empty()) for (size_t q = n; q < n2; ++q) g.swatinit[q] = -1;
    for (auto& f : g.faults) for (auto& ff : f.faces) if (ff.k2 == oldNz - 1) ff.k2 = g.nz - 1;   // faults cut both media
    if (fractureRegion) {
        int maxReg = 1 + int(c.props.satExtra.size());
        if (g.satnum.size() != n2) g.satnum.assign(n2, 1);
        for (size_t q = n; q < n2; ++q) g.satnum[q] = maxReg + 1;
        SatRegion fr;
        fr.swof = {{0.0, 0.0, 1.0, 0.0}, {1.0, 1.0, 0.0, 0.0}};
        if (c.phases.gas) fr.sgof = {{0.0, 0.0, 1.0, 0.0}, {1.0, 1.0, 0.0, 0.0}};
        while (int(c.props.satExtra.size()) < maxReg - 1) c.props.satExtra.push_back(SatRegion{c.props.swof, c.props.sgof});
        c.props.satExtra.push_back(fr);
    }
    (void)nl;
    g.dual = mode < 1 ? 1 : mode;
    g.nodppm = true;
    g.ltosigma = false;
    g.sigma = 4.0 * (2.0 / (blockL * blockL) + 1.0 / (blockH * blockH));
}

void makeSinglePorosity(Case& c) {
    GridInput& g = c.grid;
    if (!g.dual) return;
    const size_t half = size_t(g.size()) / 2;
    forCellArrays(g, [&](auto& v) { if (v.size() == 2 * half) v.resize(half); });
    if (g.tops.size() == 2 * half) g.tops.resize(half);
    if (!g.zcorn.empty()) g.zcorn.resize(g.zcorn.size() / 2);
    g.nz /= 2;
    for (auto& f : g.faults) for (auto& ff : f.faces) { ff.k1 = std::min(ff.k1, g.nz - 1); ff.k2 = std::min(ff.k2, g.nz - 1); }
    for (auto& w : c.wells) for (auto& cp : w.completions) {
        if (cp.k1 >= g.nz) cp.k1 -= g.nz;
        if (cp.k2 >= g.nz) cp.k2 -= g.nz;
    }
    g.dual = 0; g.sigma = 0; g.sigmav.clear(); g.lx.clear(); g.ly.clear(); g.lz.clear(); g.multmf.clear(); g.dzmtrxv.clear();
    g.ltosigma = false; g.nodppm = false; g.gravdr = false;
}

}  // namespace ressim
