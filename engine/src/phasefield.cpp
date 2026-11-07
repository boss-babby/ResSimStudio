// C23 geomechanics G4: phase-field fluid-driven fracture (see phasefield.h).
#include "ressim/phasefield.h"
#include "ressim/linsolve.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace ressim {
namespace phasefield {

namespace {

std::vector<double> axisCoords(double core, double h, double half, double grade) {
    std::vector<double> pos{0.0};
    double x = 0;
    while (x < core - 1e-12) { x = std::min(core, x + h); pos.push_back(x); }
    double s = h;
    while (x < half - 1e-12) {
        s *= grade;
        double nx = x + s;
        if (half - nx < 0.5 * s * grade) nx = half;
        x = std::min(half, nx); pos.push_back(x);
    }
    std::vector<double> all;
    for (size_t k = pos.size(); k-- > 1;) all.push_back(-pos[k]);
    for (double p : pos) all.push_back(p);
    return all;
}

struct Mesh {
    std::vector<double> x, y;
    int nx = 0, ny = 0, nx1 = 0, ny1 = 0;
    int node(int i, int j) const { return i + nx1 * j; }
};

// node-stencil CSR with nd dofs per node (interleaved)
CSR pattern(const Mesh& M, int nd) {
    CSR A; A.n = M.nx1 * M.ny1 * nd;
    A.rowPtr.assign(size_t(A.n) + 1, 0);
    for (int j = 0; j < M.ny1; ++j)
        for (int i = 0; i < M.nx1; ++i)
            for (int c = 0; c < nd; ++c) {
                const int r = M.node(i, j) * nd + c;
                for (int dj = -1; dj <= 1; ++dj)
                    for (int di = -1; di <= 1; ++di) {
                        const int ii = i + di, jj = j + dj;
                        if (ii < 0 || jj < 0 || ii >= M.nx1 || jj >= M.ny1) continue;
                        for (int c2 = 0; c2 < nd; ++c2) A.col.push_back(M.node(ii, jj) * nd + c2);
                    }
                A.rowPtr[size_t(r) + 1] = int(A.col.size());
            }
    // rows are generated in order r = node*nd + c only if nd loops innermost per node: rebuild sorted per row
    for (int r = 0; r < A.n; ++r) std::sort(A.col.begin() + A.rowPtr[size_t(r)], A.col.begin() + A.rowPtr[size_t(r) + 1]);
    A.val.assign(A.col.size(), 0.0);
    A.finalizeDiag();
    return A;
}
inline int pos(const CSR& A, int r, int c) {
    const int* b = A.col.data() + A.rowPtr[size_t(r)];
    const int* e = A.col.data() + A.rowPtr[size_t(r) + 1];
    const int* p = std::lower_bound(b, e, c);
    return int(p - A.col.data());
}

const double gp[2] = {-0.5773502691896258, 0.5773502691896258};

// C26: preconditioned CG - both systems are symmetric positive definite once the Dirichlet rows AND columns are
// eliminated (symmetricDirichlet); one SpMV + one V-cycle per iteration instead of two of each for BiCGSTAB
LinResult pcg(const CSR& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt) {
    const int n = A.n;
    LinResult res;
    const size_t nz = static_cast<size_t>(n);
    std::vector<double> r(nz), z(nz), p(nz), q(nz);
    A.spmv(x.data(), q.data());
    double bn = 0;
    for (int i = 0; i < n; ++i) { r[size_t(i)] = b[size_t(i)] - q[size_t(i)]; bn += b[size_t(i)] * b[size_t(i)]; }
    bn = std::sqrt(bn);
    if (bn == 0) { std::fill(x.begin(), x.end(), 0.0); res.converged = true; return res; }
    M.apply(r.data(), z.data());
    p = z;
    double rz = 0; for (int i = 0; i < n; ++i) rz += r[size_t(i)] * z[size_t(i)];
    for (int it = 1; it <= maxIt; ++it) {
        A.spmv(p.data(), q.data());
        double pq = 0; for (int i = 0; i < n; ++i) pq += p[size_t(i)] * q[size_t(i)];
        if (!(pq > 0)) break;   // not SPD (or breakdown): caller falls back to BiCGSTAB
        const double al = rz / pq;
        double rn = 0;
        for (int i = 0; i < n; ++i) { x[size_t(i)] += al * p[size_t(i)]; r[size_t(i)] -= al * q[size_t(i)]; rn += r[size_t(i)] * r[size_t(i)]; }
        res.iterations = it; res.relResidual = std::sqrt(rn) / bn;
        if (res.relResidual < tol) { res.converged = true; break; }
        M.apply(r.data(), z.data());
        double rzn = 0; for (int i = 0; i < n; ++i) rzn += r[size_t(i)] * z[size_t(i)];
        const double be = rzn / rz; rz = rzn;
        for (int i = 0; i < n; ++i) p[size_t(i)] = z[size_t(i)] + be * p[size_t(i)];
    }
    return res;
}

// zero the columns of the Dirichlet unknowns (their value is 0, so no right-hand-side correction): symmetric system
void symmetricDirichlet(CSR& A, const std::vector<char>& fixedDof) {
    for (int r = 0; r < A.n; ++r)
        for (int q = A.rowPtr[size_t(r)]; q < A.rowPtr[size_t(r) + 1]; ++q)
            if (fixedDof[size_t(A.col[size_t(q)])] && A.col[size_t(q)] != r) A.val[size_t(q)] = 0;
}

bool solveSys(const CSR& A, const std::vector<double>& b, std::vector<double>& x, std::string& err, AMG* pre = nullptr) {
    AMG own; if (!pre) own.setup(A);
    AMG& amg = pre ? *pre : own;
    LinResult r = pcg(A, amg, b, x, 1e-10, 400);
    if (!r.converged) r = bicgstab(A, amg, b, x, 1e-10, 600);
    if (!r.converged && pre) {   // refreshed hierarchy no longer fits the damaged operator: rebuild once
        pre->setup(A);
        r = bicgstab(A, amg, b, x, 1e-10, 600);
    }
    if (!r.converged) {
        ILU0 ilu; ilu.factor(A);
        r = bicgstab(A, ilu, b, x, 1e-10, 4000);
    }
    if (std::getenv("PF_DEBUG")) std::fprintf(stderr, "    lin n %d it %d res %.2e conv %d\n", A.n, r.iterations, r.relResidual, int(r.converged));
    if (!r.converged && r.relResidual > 1e-6) { err = "linear solver did not converge"; return false; }
    return true;
}

}  // namespace

bool run(const Input& in, Output& out, std::atomic<bool>* cancel, const std::function<void(const std::string&)>& log) {
    const auto t0 = std::chrono::steady_clock::now();
    out = Output{};
    const double a = in.crackA, h = in.h;
    if (!(a > 0) || !(h > 0) || !(in.E > 0) || !(in.Gc > 0)) { out.error = "invalid input"; return false; }
    const double l = in.ell > 0 ? in.ell : 2 * h;
    // C26: crack direction (local s along, n normal) and the rotated core box that contains it
    const double th = in.angleDeg * M_PI / 180, cth = std::cos(th), sth = std::sin(th);
    const bool inclined = std::fabs(sth) > 1e-9;
    const double cxh = inclined ? in.coreX * a * std::fabs(cth) + in.coreY * a * std::fabs(sth) : in.coreX * a;
    const double cyh = inclined ? in.coreX * a * std::fabs(sth) + in.coreY * a * std::fabs(cth) : in.coreY * a;
    // load tensor T = p I - S (tension positive): T = p I1 + TS
    const double TSxx = -in.sxx, TSyy = -in.syy, TSxy = -in.sxy;
    const bool stressed = in.sxx != 0 || in.syy != 0 || in.sxy != 0;
    out.sigmaN = in.sxx * sth * sth + in.syy * cth * cth - 2 * in.sxy * sth * cth;
    out.tauS = (in.syy - in.sxx) * sth * cth + in.sxy * (cth * cth - sth * sth);
    out.controlUsed = in.control;
    Mesh M;
    M.x = axisCoords(cxh, h, std::max(in.domain * a, 2 * cxh), in.grade);
    M.y = axisCoords(cyh, h, std::max(in.domain * a, 2 * cyh), in.grade);
    M.nx1 = int(M.x.size()); M.ny1 = int(M.y.size()); M.nx = M.nx1 - 1; M.ny = M.ny1 - 1;
    const int nn = M.nx1 * M.ny1;
    if (nn > 400000) { out.error = "mesh too large"; return false; }
    out.x = M.x; out.y = M.y; out.nodes = nn; out.ell = l; out.h = h;
    const double Ep = in.E / (1 - in.nu * in.nu);
    out.Eprime = Ep; out.GcEff = in.Gc * (1 + 3 * h / (8 * l));
    // Griffith reference with the physical Gc: the measured p_c depends on l / a (0.85-0.87 at l/a = 0.2, 0.93-0.97 at 0.1)
    // and hardly on h at fixed l (1.5 %), so the Bourdin discretisation factor GcEff is reported but not applied
    const double pG = std::sqrt(Ep * in.Gc / (M_PI * a));
    out.pcTheory = out.sigmaN + std::sqrt(std::max(0.0, pG * pG - out.tauS * out.tauS));
    const double pMax = in.pMax > 0 ? in.pMax : out.sigmaN + 1.3 * (out.pcTheory - out.sigmaN);
    const double vMax = in.vMax > 0 ? in.vMax : 6 * 2 * M_PI * (out.pcTheory - out.sigmaN) * a * a / Ep;
    const double lam = in.E * in.nu / ((1 + in.nu) * (1 - 2 * in.nu)), mu = in.E / (2 * (1 + in.nu)), Kb = lam + mu;   // plane strain: K_2D = lam + mu
    int j0 = 0; for (int j = 0; j < M.ny1; ++j) if (std::fabs(M.y[size_t(j)]) < 1e-12) j0 = j;
    std::vector<char> crackNode(size_t(nn), 0), bnd(size_t(nn), 0);
    // initial crack: three node rows (|y| <= h) so that the elements on both sides of the crack line are fully broken
    // (one damaged row leaves the adjacent elements with g = O(0.01) at their Gauss points, which resists the opening the
    // more the finer the mesh - the Sneddon volume then decays like (l / a)^0.6)
    if (!inclined) {
        for (int i = 0; i < M.nx1; ++i)
            if (std::fabs(M.x[size_t(i)]) <= a + 1e-12)
                for (int dj = -1; dj <= 1; ++dj) if (j0 + dj >= 0 && j0 + dj < M.ny1) crackNode[size_t(M.node(i, j0 + dj))] = 1;
    } else {   // inclined: the band |n| <= h (|cos| + |sin|) keeps a connected, three-node-wide crack on the rectilinear mesh
        const double band = h * (std::fabs(cth) + std::fabs(sth)) + 1e-12;
        for (int j = 0; j < M.ny1; ++j) for (int i = 0; i < M.nx1; ++i) {
            const double X = M.x[size_t(i)], Y = M.y[size_t(j)];
            if (std::fabs(X * cth + Y * sth) <= a + 1e-12 && std::fabs(-X * sth + Y * cth) <= band) crackNode[size_t(M.node(i, j))] = 1;
        }
    }
    for (int j = 0; j < M.ny1; ++j) for (int i = 0; i < M.nx1; ++i) if (i == 0 || j == 0 || i == M.nx || j == M.ny) bnd[size_t(M.node(i, j))] = 1;
    CSR Ku = pattern(M, 2), Kv = pattern(M, 1);
    CSR Kdec = Ku;   // C26: component-decoupled copy for the AMG hierarchy
    std::vector<char> fixU(size_t(2 * nn), 0);
    for (int k = 0; k < nn; ++k) if (bnd[size_t(k)]) fixU[size_t(2 * k)] = fixU[size_t(2 * k) + 1] = 1;
    std::vector<double> u(size_t(2 * nn), 0.0), v(size_t(nn), 1.0), vn(size_t(nn), 1.0);
    for (int k = 0; k < nn; ++k) if (crackNode[size_t(k)]) v[size_t(k)] = vn[size_t(k)] = 0;
    const bool frozen = std::getenv("PF_FROZEN") != nullptr;   // debug: analytic AT2 profile, no phase-field solve
    if (frozen)
        for (int j = 0; j < M.ny1; ++j) for (int i = 0; i < M.nx1; ++i) {
            const double xx = std::max(0.0, std::fabs(M.x[size_t(i)]) - a), d = std::hypot(xx, std::max(0.0, std::fabs(M.y[size_t(j)]) - h));
            v[size_t(M.node(i, j))] = vn[size_t(M.node(i, j))] = 1 - std::exp(-d / (2 * l));
        }
    const size_t ng = size_t(M.nx) * M.ny * 4;
    std::vector<double> H(ng, 0.0), Hn(ng, 0.0);
    std::vector<double> Rinc(size_t(M.nx) * M.ny * 32, 0.0);   // per element Kaa^-1 Kau of the incompatible modes
    std::string err;

    auto elemLoop = [&](auto&& f) {
        for (int j = 0; j < M.ny; ++j)
            for (int i = 0; i < M.nx; ++i) {
                const double dx = M.x[size_t(i) + 1] - M.x[size_t(i)], dy = M.y[size_t(j) + 1] - M.y[size_t(j)];
                const int nd[4] = {M.node(i, j), M.node(i + 1, j), M.node(i + 1, j + 1), M.node(i, j + 1)};
                f(i, j, dx, dy, nd);
            }
    };
    static const double sx[4] = {-1, 1, 1, -1}, sy[4] = {-1, -1, 1, 1};
    auto shape = [](double xi, double et, double dx, double dy, double N[4], double Nx[4], double Ny[4]) {
        for (int q = 0; q < 4; ++q) {
            N[q] = 0.25 * (1 + sx[q] * xi) * (1 + sy[q] * et);
            Nx[q] = 0.25 * sx[q] * (1 + sy[q] * et) * 2 / dx;
            Ny[q] = 0.25 * sy[q] * (1 + sx[q] * xi) * 2 / dy;
        }
    };
    auto crackVolumeOf = [&](const std::vector<double>& u) {
        double V = 0;
        elemLoop([&](int, int, double dx, double dy, const int* nd) {
            for (int gy = 0; gy < 2; ++gy) for (int gx = 0; gx < 2; ++gx) {
                double N[4], Nx[4], Ny[4]; shape(gp[gx], gp[gy], dx, dy, N, Nx, Ny);
                double ux = 0, uy = 0, vx = 0, vy = 0;
                for (int q = 0; q < 4; ++q) { ux += N[q] * u[size_t(2 * nd[q])]; uy += N[q] * u[size_t(2 * nd[q] + 1)]; vx += Nx[q] * v[size_t(nd[q])]; vy += Ny[q] * v[size_t(nd[q])]; }
                V += (ux * vx + uy * vy) * 0.25 * dx * dy;
            }
        });
        return V;
    };
    auto crackVolume = [&]() { return crackVolumeOf(u); };
    // C26 general measure: largest distance from the centre of a v = 0.5 crossing (edge interpolation) or of a broken node;
    // also the direction of the furthest point on the + side (s >= 0) - follows a curving crack
    double tipAng = in.angleDeg;
    auto crackRadial = [&]() {
        double L = 0, Lp = -1;
        auto take = [&](double X, double Y) {
            const double r = std::hypot(X, Y);
            L = std::max(L, r);
            if (X * cth + Y * sth >= 0 && r > Lp) { Lp = r; tipAng = std::atan2(Y, X) * 180 / M_PI; }
        };
        for (int j = 0; j < M.ny1; ++j) for (int i = 0; i < M.nx1; ++i) {
            const double v1 = v[size_t(M.node(i, j))];
            if (v1 < 0.5) take(M.x[size_t(i)], M.y[size_t(j)]);
            if (i + 1 < M.nx1) { const double v2 = v[size_t(M.node(i + 1, j))];
                if ((v1 - 0.5) * (v2 - 0.5) < 0) take(M.x[size_t(i)] + (0.5 - v1) / (v2 - v1) * (M.x[size_t(i) + 1] - M.x[size_t(i)]), M.y[size_t(j)]); }
            if (j + 1 < M.ny1) { const double v2 = v[size_t(M.node(i, j + 1))];
                if ((v1 - 0.5) * (v2 - 0.5) < 0) take(M.x[size_t(i)], M.y[size_t(j)] + (0.5 - v1) / (v2 - v1) * (M.y[size_t(j) + 1] - M.y[size_t(j)])); }
        }
        return L;
    };
    // damage reached the edge of the uniform core (the graded mesh outside cannot carry the crack)
    auto leftCore = [&]() {
        for (int j = 0; j < M.ny1; ++j) for (int i = 0; i < M.nx1; ++i)
            if (v[size_t(M.node(i, j))] < 0.5 && (std::fabs(M.x[size_t(i)]) > 0.97 * cxh || std::fabs(M.y[size_t(j)]) > 0.97 * cyh)) return true;
        return false;
    };
    auto crackHalfRow = [&]() {
        double L = 0;
        for (int i = 0; i + 1 < M.nx1; ++i) {
            const double v1 = v[size_t(M.node(i, j0))], v2 = v[size_t(M.node(i + 1, j0))];
            if (v1 < 0.5) L = std::max(L, std::fabs(M.x[size_t(i)]));
            if (v2 < 0.5) L = std::max(L, std::fabs(M.x[size_t(i) + 1]));
            if ((v1 - 0.5) * (v2 - 0.5) < 0) {
                const double xc = M.x[size_t(i)] + (0.5 - v1) / (v2 - v1) * (M.x[size_t(i) + 1] - M.x[size_t(i)]);
                L = std::max(L, std::fabs(xc));
            }
        }
        return L;
    };

    auto crackHalf = [&]() { return (inclined || stressed) ? crackRadial() : crackHalfRow(); };
    const bool volCtl = in.control == 1;
    const bool general = inclined || stressed || volCtl;   // C23 path (horizontal, unstressed, pressure ramp) unchanged otherwise
    std::vector<double> u1(size_t(2 * nn), 0.0), uS(size_t(2 * nn), 0.0);   // unit-pressure and in-situ-stress responses
    AMG amgU, amgV;
    bool grown = false;
    double pPrev = 0, pPeak = -1e300;
    const bool external = volCtl && bool(in.nextVolume);   // C28 in-run coupling: open-ended, targets from the flow run
    for (int st = 1; st <= in.steps || external; ++st) {
        if (cancel && cancel->load()) { out.error = "cancelled"; break; }
        double p = volCtl ? 0.0 : pMax * st / in.steps;
        double Vt = vMax * st / std::max(1, in.steps);
        if (external && !in.nextVolume(st, Vt)) break;
        Step S; S.p = p;
        int it = 0; double dvMax = 1;
        for (; it < in.maxStagger; ++it) {
            // ---- momentum: Q4 + Wilson incompatible modes (1 - xi^2, 1 - eta^2; Taylor et al. 1976), condensed per element.
            // The incompatible modes remove the parasitic shear stiffness of the elongated elements of the graded mesh
            // (plain Q4 under-predicted the Sneddon volume by 40 % at aspect ratio ~40).
            std::fill(Ku.val.begin(), Ku.val.end(), 0.0);
            std::vector<double> fu(size_t(2 * nn), 0.0), fs(size_t(2 * nn), 0.0);   // unit pressure / in-situ stress loads
            {
                size_t ei = 0;
                elemLoop([&](int, int, double dx, double dy, const int* nd) {
                    double K[12][12] = {}; double fe[8] = {}, fse[8] = {};
                    for (int gy = 0; gy < 2; ++gy) for (int gx = 0; gx < 2; ++gx) {
                        double N[4], Nx[4], Ny[4]; shape(gp[gx], gp[gy], dx, dy, N, Nx, Ny);
                        double vg = 0, vx = 0, vy = 0;
                        for (int q = 0; q < 4; ++q) { vg += N[q] * v[size_t(nd[q])]; vx += Nx[q] * v[size_t(nd[q])]; vy += Ny[q] * v[size_t(nd[q])]; }
                        const double g = vg * vg + in.kres, w = 0.25 * dx * dy;
                        // shape-function derivatives of the 6 "nodes" (4 nodal + 2 incompatible), dofs: 2 per entry
                        double Gx[6], Gy[6];
                        for (int q = 0; q < 4; ++q) { Gx[q] = Nx[q]; Gy[q] = Ny[q]; }
                        Gx[4] = -2 * gp[gx] * 2 / dx; Gy[4] = 0;
                        Gx[5] = 0; Gy[5] = -2 * gp[gy] * 2 / dy;
                        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
                            K[2 * A][2 * B] += w * g * ((lam + 2 * mu) * Gx[A] * Gx[B] + mu * Gy[A] * Gy[B]);
                            K[2 * A][2 * B + 1] += w * g * (lam * Gx[A] * Gy[B] + mu * Gy[A] * Gx[B]);
                            K[2 * A + 1][2 * B] += w * g * (lam * Gy[A] * Gx[B] + mu * Gx[A] * Gy[B]);
                            K[2 * A + 1][2 * B + 1] += w * g * ((lam + 2 * mu) * Gy[A] * Gy[B] + mu * Gx[A] * Gx[B]);
                        }
                        // face traction (p I - S) grad v: unit pressure part and in-situ part
                        for (int A = 0; A < 4; ++A) {
                            fe[2 * A] += w * N[A] * vx; fe[2 * A + 1] += w * N[A] * vy;
                            fse[2 * A] += w * N[A] * (TSxx * vx + TSxy * vy); fse[2 * A + 1] += w * N[A] * (TSxy * vx + TSyy * vy);
                        }
                    }
                    // condense: Kaa^-1 Kau (4 x 8) by Gauss-Jordan on [Kaa | Kau]
                    double M4[4][12];
                    for (int r = 0; r < 4; ++r) { for (int c = 0; c < 4; ++c) M4[r][c] = K[8 + r][8 + c]; for (int c = 0; c < 8; ++c) M4[r][4 + c] = K[8 + r][c]; }
                    for (int c = 0; c < 4; ++c) {
                        int pv = c; for (int r = c + 1; r < 4; ++r) if (std::fabs(M4[r][c]) > std::fabs(M4[pv][c])) pv = r;
                        if (pv != c) for (int k = 0; k < 12; ++k) std::swap(M4[c][k], M4[pv][k]);
                        const double d = M4[c][c];
                        for (int k = 0; k < 12; ++k) M4[c][k] /= d;
                        for (int r = 0; r < 4; ++r) if (r != c) { const double f = M4[r][c]; for (int k = 0; k < 12; ++k) M4[r][k] -= f * M4[c][k]; }
                    }
                    double* R = &Rinc[ei * 32];
                    for (int r = 0; r < 4; ++r) for (int c = 0; c < 8; ++c) R[r * 8 + c] = M4[r][4 + c];
                    for (int A = 0; A < 8; ++A) {
                        const int r = 2 * nd[A / 2] + A % 2;
                        fu[size_t(r)] += fe[A]; fs[size_t(r)] += fse[A];
                        for (int B = 0; B < 8; ++B) {
                            double kab = K[A][B];
                            for (int q = 0; q < 4; ++q) kab -= K[A][8 + q] * R[q * 8 + B];
                            Ku.val[size_t(pos(Ku, r, 2 * nd[B / 2] + B % 2))] += kab;
                        }
                    }
                    ++ei;
                });
            }
            for (int k = 0; k < nn; ++k) if (bnd[size_t(k)])
                for (int c = 0; c < 2; ++c) {
                    const int r = 2 * k + c;
                    for (int q = Ku.rowPtr[size_t(r)]; q < Ku.rowPtr[size_t(r) + 1]; ++q) Ku.val[size_t(q)] = Ku.col[size_t(q)] == r ? 1.0 : 0.0;
                    fu[size_t(r)] = 0; fs[size_t(r)] = 0;
                }
            symmetricDirichlet(Ku, fixU);
            {   // hierarchy built at the first staggered iteration of a load step, values refreshed afterwards
                // unknown-based AMG: the hierarchy is built on the component-decoupled operator (x-y couplings dropped), so the
                // aggregates never mix the two displacement components (translations stay in the coarse space); C26: PCG
                // iterations 150 -> see ENGINES 6.20
                Kdec.val = Ku.val;
                for (int r = 0; r < Ku.n; ++r)
                    for (int q = Ku.rowPtr[size_t(r)]; q < Ku.rowPtr[size_t(r) + 1]; ++q) if ((Ku.col[size_t(q)] & 1) != (r & 1)) Kdec.val[size_t(q)] = 0;
                if (it == 0 || !amgU.refresh(Kdec)) amgU.setup(Kdec);
                if (!solveSys(Ku, fu, u1, err, &amgU)) { out.error = err; break; }
                if (stressed && !solveSys(Ku, fs, uS, err, &amgU)) { out.error = err; break; }
            }
            if (volCtl) {   // u = p u1 + uS with V(u) = V_target (V is linear in u)
                const double V1 = crackVolumeOf(u1), VS = stressed ? crackVolumeOf(uS) : 0.0;
                if (!(V1 > 0)) { out.error = "volume control: no crack opening under unit pressure"; break; }
                p = (Vt - VS) / V1;
            }
            for (size_t k = 0; k < u.size(); ++k) u[k] = p * u1[k] + (stressed ? uS[k] : 0.0);
            // ---- history
            {
                size_t gi = 0, ei = 0;
                elemLoop([&](int, int, double dx, double dy, const int* nd) {
                    double ue[8], al[4] = {0, 0, 0, 0};
                    for (int A = 0; A < 8; ++A) ue[A] = u[size_t(2 * nd[A / 2] + A % 2)];
                    const double* R = &Rinc[ei * 32];
                    for (int q = 0; q < 4; ++q) for (int B = 0; B < 8; ++B) al[q] -= R[q * 8 + B] * ue[B];   // al = (a1x, a1y, a2x, a2y)
                    ++ei;
                    for (int gy = 0; gy < 2; ++gy) for (int gx = 0; gx < 2; ++gx, ++gi) {
                        double N[4], Nx[4], Ny[4]; shape(gp[gx], gp[gy], dx, dy, N, Nx, Ny);
                        double exx = 0, eyy = 0, gxy = 0;
                        for (int q = 0; q < 4; ++q) { exx += Nx[q] * u[size_t(2 * nd[q])]; eyy += Ny[q] * u[size_t(2 * nd[q] + 1)]; gxy += Ny[q] * u[size_t(2 * nd[q])] + Nx[q] * u[size_t(2 * nd[q] + 1)]; }
                        const double p1x = -2 * gp[gx] * 2 / dx, p2y = -2 * gp[gy] * 2 / dy;
                        exx += p1x * al[0]; eyy += p2y * al[3]; gxy += p1x * al[1] + p2y * al[2];
                        const double tr = exx + eyy, tp = std::max(0.0, tr);
                        // plane strain: eps_zz = 0, dev eps in 3-D
                        const double m3 = tr / 3, dxx = exx - m3, dyy = eyy - m3, dzz = -m3;
                        const double Kv3 = lam + 2.0 * mu / 3.0;   // 3-D bulk modulus
                        const double psiP = 0.5 * Kv3 * tp * tp + mu * (dxx * dxx + dyy * dyy + dzz * dzz + 0.5 * gxy * gxy);
                        H[gi] = std::max(Hn[gi], psiP);
                    }
                });
                (void)Kb;
            }
            // ---- phase field
            std::fill(Kv.val.begin(), Kv.val.end(), 0.0);
            std::vector<double> fv(size_t(nn), 0.0);
            {
                size_t gi = 0;
                elemLoop([&](int, int, double dx, double dy, const int* nd) {
                    double Ke[4][4] = {}; double fe[4] = {};
                    for (int gy = 0; gy < 2; ++gy) for (int gx = 0; gx < 2; ++gx, ++gi) {
                        double N[4], Nx[4], Ny[4]; shape(gp[gx], gp[gy], dx, dy, N, Nx, Ny);
                        double ux = 0, uy = 0;
                        for (int q = 0; q < 4; ++q) { ux += N[q] * u[size_t(2 * nd[q])]; uy += N[q] * u[size_t(2 * nd[q] + 1)]; }
                        const double w = 0.25 * dx * dy, cm = 2 * H[gi] + in.Gc / (2 * l), cd = 2 * l * in.Gc;
                        for (int A = 0; A < 4; ++A) {
                            const double Txx = p + TSxx, Tyy = p + TSyy, Txy = TSxy;   // u . (p I - S) . grad N
                            fe[A] += w * (in.Gc / (2 * l) * N[A] + ux * (Txx * Nx[A] + Txy * Ny[A]) + uy * (Txy * Nx[A] + Tyy * Ny[A]));
                            for (int B = 0; B < 4; ++B) Ke[A][B] += w * (cm * N[A] * N[B] + cd * (Nx[A] * Nx[B] + Ny[A] * Ny[B]));
                        }
                    }
                    for (int A = 0; A < 4; ++A) {
                        fv[size_t(nd[A])] += fe[A];
                        for (int B = 0; B < 4; ++B) Kv.val[size_t(pos(Kv, nd[A], nd[B]))] += Ke[A][B];
                    }
                });
            }
            for (int k = 0; k < nn; ++k) if (crackNode[size_t(k)]) {
                for (int q = Kv.rowPtr[size_t(k)]; q < Kv.rowPtr[size_t(k) + 1]; ++q) Kv.val[size_t(q)] = Kv.col[size_t(q)] == k ? 1.0 : 0.0;
                fv[size_t(k)] = 0;
            }
            symmetricDirichlet(Kv, crackNode);
            std::vector<double> vNew = v;
            if (frozen) { ++it; break; }
            if (it == 0 || !amgV.refresh(Kv)) amgV.setup(Kv);
            if (!solveSys(Kv, fv, vNew, err, &amgV)) { out.error = err; break; }
            dvMax = 0;
            for (int k = 0; k < nn; ++k) {
                const double c = std::clamp(std::min(vNew[size_t(k)], vn[size_t(k)]), 0.0, 1.0);
                dvMax = std::max(dvMax, std::fabs(c - v[size_t(k)]));
                v[size_t(k)] = c;
            }
            if (dvMax < in.tol) { ++it; break; }
            if (crackHalf() > in.stopLength * in.coreX * a || (general && leftCore())) { ++it; break; }
        }
        if (!out.error.empty()) break;
        S.p = p;
        S.stagger = it; S.converged = dvMax < in.tol;
        Hn = H; vn = v;
        S.volume = crackVolume(); S.lengthHalf = crackHalf();
        crackRadial(); out.tipAngleDeg = tipAng;
        out.steps.push_back(S);
        if (in.onStep) in.onStep(S, tipAng);
        if (in.onField) {   // C30: nodal phase field / displacement, element stress and the opening profile
            Field F;
            F.x = M.x; F.y = M.y; F.p = p; F.volume = S.volume; F.half = S.lengthHalf; F.angleDeg = tipAng;
            F.v.resize(size_t(nn)); F.ux.resize(size_t(nn)); F.uy.resize(size_t(nn));
            for (int k = 0; k < nn; ++k) { F.v[size_t(k)] = float(v[size_t(k)]); F.ux[size_t(k)] = float(u[size_t(2 * k)]); F.uy[size_t(k)] = float(u[size_t(2 * k + 1)]); }
            const size_t ne = size_t(M.nx) * size_t(M.ny);
            F.sxx.assign(ne, 0.0f); F.syy.assign(ne, 0.0f); F.sxy.assign(ne, 0.0f);
            std::vector<double> col(size_t(M.nx), 0.0);
            elemLoop([&](int i, int j, double dx, double dy, const int* nd) {
                double N[4], Nx[4], Ny[4]; shape(0.0, 0.0, dx, dy, N, Nx, Ny);
                double exx = 0, eyy = 0, gxy = 0, vc = 0;
                for (int q = 0; q < 4; ++q) {
                    const double ux = u[size_t(2 * nd[q])], uy = u[size_t(2 * nd[q] + 1)];
                    exx += Nx[q] * ux; eyy += Ny[q] * uy; gxy += Ny[q] * ux + Nx[q] * uy; vc += 0.25 * v[size_t(nd[q])];
                }
                const double gdeg = vc * vc + in.kres, tr = exx + eyy;
                const size_t e = size_t(i) + size_t(M.nx) * size_t(j);
                F.sxx[e] = float(gdeg * (lam * tr + 2 * mu * exx));
                F.syy[e] = float(gdeg * (lam * tr + 2 * mu * eyy));
                F.sxy[e] = float(gdeg * mu * gxy);
                if (!inclined)
                    for (int gy = 0; gy < 2; ++gy) for (int gx = 0; gx < 2; ++gx) {
                        double N2[4], Nx2[4], Ny2[4]; shape(gp[gx], gp[gy], dx, dy, N2, Nx2, Ny2);
                        double ux = 0, uy = 0, vx = 0, vy = 0;
                        for (int q = 0; q < 4; ++q) { ux += N2[q] * u[size_t(2 * nd[q])]; uy += N2[q] * u[size_t(2 * nd[q] + 1)]; vx += Nx2[q] * v[size_t(nd[q])]; vy += Ny2[q] * v[size_t(nd[q])]; }
                        col[size_t(i)] += (ux * vx + uy * vy) * 0.25 * dy;   // per unit length along x
                    }
            });
            if (!inclined) {
                F.codX.resize(size_t(M.nx)); F.cod.resize(size_t(M.nx));
                for (int i = 0; i < M.nx; ++i) { F.codX[size_t(i)] = float(0.5 * (M.x[size_t(i)] + M.x[size_t(i) + 1])); F.cod[size_t(i)] = float(std::max(0.0, col[size_t(i)])); }
            }
            in.onField(F);
        }
        if (log) { char b[220]; std::snprintf(b, sizeof b, "%sp %.4g MPa: V %.4g m2, half-length %.4g m, %d staggered iterations", volCtl ? "volume step: " : "", p * 1e-6, S.volume, S.lengthHalf, it); log(b); }
        // linear-regime Sneddon check from the first (lowest-pressure) step: AT2 has no elastic threshold, so V/p creeps up
        // as the tips soften before the unstable jump
        if (out.steps.size() == 1 && S.p - out.sigmaN > 0) out.sneddonRatio = S.volume / (2 * M_PI * (S.p - out.sigmaN) * a * a / Ep);
        if (volCtl) {   // stable growth: the critical pressure is the peak of the pressure record
            if (p > pPeak) { pPeak = p; out.pcNumeric = p; out.pcBracketLo = pPrev; out.pcBracketHi = p; out.vAtPc = v; grown = true; }
            if (S.lengthHalf > in.stopLength * in.coreX * a || (general && leftCore())) break;
            pPrev = p;
            continue;
        }
        // unstable growth under pressure control = a jump of the half-length by more than max(4h, 0.1a) in one load step
        const double prevHalf = out.steps.size() >= 2 ? out.steps[out.steps.size() - 2].lengthHalf : S.lengthHalf;
        if (!grown && S.lengthHalf - prevHalf > std::max(4 * h, 0.1 * a)) {
            grown = true;
            out.pcNumeric = 0.5 * (pPrev + p);
            out.vAtPc = v;
            out.pcBracketLo = pPrev; out.pcBracketHi = p;
        }
        if (S.lengthHalf > in.stopLength * in.coreX * a || (general && leftCore())) break;
        pPrev = p;
    }
    out.v = v;
    out.ok = out.error.empty();
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out.ok;
}

}  // namespace phasefield
}  // namespace ressim
