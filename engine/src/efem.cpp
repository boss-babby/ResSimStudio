// C20 geomechanics G1: 2-D EFEM with frictional contact (see efem.h)
#include "ressim/efem.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace ressim {
namespace efem {
namespace {

// banded matrix (no pivoting: the condensed elasticity matrix is close to SPD), row-major band storage
struct Band {
    int n = 0, kl = 0, ku = 0, w = 0;
    std::vector<double> a;
    void init(int n_, int kl_, int ku_) { n = n_; kl = kl_; ku = ku_; w = kl + ku + 1; a.assign(size_t(n) * w, 0.0); }
    double& at(int i, int j) { return a[size_t(i) * w + (j - i + kl)]; }
    bool inBand(int i, int j) const { return j - i <= ku && i - j <= kl; }
    bool factor() {
        for (int k = 0; k < n; ++k) {
            const double piv = at(k, k);
            if (!(std::fabs(piv) > 0) || !std::isfinite(piv)) return false;
            const int iEnd = std::min(n - 1, k + kl), jEnd = std::min(n - 1, k + ku);
            for (int i = k + 1; i <= iEnd; ++i) {
                double& lik = at(i, k);
                if (lik == 0) continue;
                lik /= piv;
                const double l = lik;
                double* ri = &at(i, k + 1);
                const double* rk = &at(k, k + 1);
                for (int j = k + 1; j <= jEnd; ++j) ri[j - k - 1] -= l * rk[j - k - 1];
            }
        }
        return true;
    }
    void solve(std::vector<double>& b) {
        for (int i = 0; i < n; ++i) { const int j0 = std::max(0, i - kl); double s = b[size_t(i)]; for (int j = j0; j < i; ++j) s -= at(i, j) * b[size_t(j)]; b[size_t(i)] = s; }
        for (int i = n - 1; i >= 0; --i) { const int j1 = std::min(n - 1, i + ku); double s = b[size_t(i)]; for (int j = i + 1; j <= j1; ++j) s -= at(i, j) * b[size_t(j)]; b[size_t(i)] = s / at(i, i); }
    }
};

// C21 G1-c: fluid pressure of fracture F at (x, y) (uniform, or the along-fracture profile)
inline double pressureAt(const Fracture& F, double x, double y) {
    const size_t n = F.profS.size();
    if (n < 2 || F.profP.size() != n) return F.pressure;
    const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
    if (!(L > 0)) return F.profP[0];
    const double s = ((x - 0.5 * (F.x0 + F.x1)) * (F.x1 - F.x0) + (y - 0.5 * (F.y0 + F.y1)) * (F.y1 - F.y0)) / L;
    if (s <= F.profS[0]) return F.profP[0];
    if (s >= F.profS[n - 1]) return F.profP[n - 1];
    size_t k = 1; while (k < n - 1 && F.profS[k] < s) ++k;
    const double w = (s - F.profS[k - 1]) / std::max(F.profS[k] - F.profS[k - 1], 1e-300);
    return F.profP[k - 1] + w * (F.profP[k] - F.profP[k - 1]);
}
inline double pressureMax(const Fracture& F) {
    double m = std::fabs(F.pressure);
    if (F.profS.size() >= 2 && F.profP.size() == F.profS.size()) { m = 0; for (double v : F.profP) m = std::max(m, std::fabs(v)); }
    return m;
}

struct Cut {
    int e = 0, f = 0;
    double n[2] = {0, 0}, m[2] = {0, 0}, len = 0, area = 0, xc = 0, yc = 0, s = 0;
    double Kuw[8][2] = {}, Gwu[2][8] = {}, Gww[2][2] = {};
    double w[2] = {0, 0}, wtp = 0, t[2] = {0, 0}, T[2][2] = {};
    int state = 0;
};

inline bool inv2(const double A[2][2], double B[2][2]) {
    const double d = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (!(std::fabs(d) > 0) || !std::isfinite(d)) return false;
    B[0][0] = A[1][1] / d; B[0][1] = -A[0][1] / d; B[1][0] = -A[1][0] / d; B[1][1] = A[0][0] / d;
    return true;
}

// segment / axis-aligned rectangle clip (Liang-Barsky): parameter range [ta, tb] of p0 + t (p1 - p0) inside
bool clip(double x0, double y0, double x1, double y1, double xa, double ya, double xb, double yb, double& ta, double& tb) {
    ta = 0; tb = 1;
    const double dx = x1 - x0, dy = y1 - y0;
    const double p[4] = {-dx, dx, -dy, dy}, q[4] = {x0 - xa, xb - x0, y0 - ya, yb - y0};
    for (int k = 0; k < 4; ++k) {
        if (p[k] == 0) { if (q[k] < 0) return false; continue; }
        const double r = q[k] / p[k];
        if (p[k] < 0) ta = std::max(ta, r); else tb = std::min(tb, r);
    }
    return tb > ta;
}
}  // namespace

static bool solveConstantJump(const Input& in, Result& out, const Result* prev) {
    const auto t0 = std::chrono::steady_clock::now();
    out = Result();
    auto fail = [&](const std::string& m) { out.error = m; return false; };
    if (in.nx < 2 || in.ny < 2 || !(in.Lx > 0) || !(in.Ly > 0)) return fail("mesh / domain out of range");
    if (!(in.E > 0) || !(in.nu >= 0 && in.nu < 0.5)) return fail("elastic constants out of range");
    const int nx = in.nx, ny = in.ny, nnx = nx + 1, nn = nnx * (ny + 1), N = 2 * nn;
    const double hx = in.Lx / nx, hy = in.Ly / ny, A = hx * hy;
    const double c = in.E / ((1 + in.nu) * (1 - 2 * in.nu));
    const double D[3][3] = {{c * (1 - in.nu), c * in.nu, 0}, {c * in.nu, c * (1 - in.nu), 0}, {0, 0, c * (1 - 2 * in.nu) / 2}};
    const double Kpen = in.penalty * in.E / std::min(hx, hy);
    // reference element: B at the 2x2 Gauss points (rectangles: identical for all elements)
    const double g = 1 / std::sqrt(3.0), xiA[4] = {-1, 1, 1, -1}, etA[4] = {-1, -1, 1, 1};
    double Bg[4][3][8] = {}, dNx[4][4], dNy[4][4];
    for (int q = 0; q < 4; ++q) {
        const double xi = (q == 1 || q == 2) ? g : -g, et = q >= 2 ? g : -g;
        for (int a = 0; a < 4; ++a) {
            dNx[q][a] = xiA[a] * (1 + et * etA[a]) / 4 * 2 / hx;
            dNy[q][a] = etA[a] * (1 + xi * xiA[a]) / 4 * 2 / hy;
            Bg[q][0][2 * a] = dNx[q][a]; Bg[q][1][2 * a + 1] = dNy[q][a];
            Bg[q][2][2 * a] = dNy[q][a]; Bg[q][2][2 * a + 1] = dNx[q][a];
        }
    }
    const double wq = A / 4;
    double Ke[8][8] = {};
    for (int q = 0; q < 4; ++q) {
        double DB[3][8] = {};
        for (int r = 0; r < 3; ++r) for (int k = 0; k < 8; ++k) for (int s = 0; s < 3; ++s) DB[r][k] += D[r][s] * Bg[q][s][k];
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) { double v = 0; for (int r = 0; r < 3; ++r) v += Bg[q][r][i] * DB[r][j]; Ke[i][j] += wq * v; }
    }
    auto enodes = [&](int e, int nd[4]) { const int i = e % nx, j = e / nx; nd[0] = i + nnx * j; nd[1] = nd[0] + 1; nd[2] = nd[1] + nnx; nd[3] = nd[0] + nnx; };
    auto nodeXY = [&](int a, double& x, double& y) { x = (a % nnx) * hx; y = (a / nnx) * hy; };
    // ---- cut elements
    std::vector<Cut> cuts;
    std::vector<int> cutOf(size_t(nx) * ny, -1);
    for (size_t f = 0; f < in.fracs.size(); ++f) {
        const Fracture& F = in.fracs[f];
        const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
        if (!(L > 0)) return fail("fracture of zero length");
        const double mx = (F.x1 - F.x0) / L, my = (F.y1 - F.y0) / L;
        const int i0 = std::max(0, int(std::floor(std::min(F.x0, F.x1) / hx)) - 1), i1 = std::min(nx - 1, int(std::floor(std::max(F.x0, F.x1) / hx)) + 1);
        const int j0 = std::max(0, int(std::floor(std::min(F.y0, F.y1) / hy)) - 1), j1 = std::min(ny - 1, int(std::floor(std::max(F.y0, F.y1) / hy)) + 1);
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) {
                double ta, tb;
                if (!clip(F.x0, F.y0, F.x1, F.y1, i * hx, j * hy, (i + 1) * hx, (j + 1) * hy, ta, tb)) continue;
                const double len = (tb - ta) * L;
                if (len < 1e-6 * std::min(hx, hy)) continue;
                const int e = i + nx * j;
                if (cutOf[size_t(e)] >= 0) return fail("two fractures cut the same element (one embedded jump per element)");
                Cut C; C.e = e; C.f = int(f); C.len = len; C.area = A;
                C.m[0] = mx; C.m[1] = my; C.n[0] = -my; C.n[1] = mx;
                const double tc = 0.5 * (ta + tb);
                C.xc = F.x0 + tc * (F.x1 - F.x0); C.yc = F.y0 + tc * (F.y1 - F.y0); C.s = (tc - 0.5) * L;
                // grad phi at the Gauss points (phi = sum of the shape functions of the nodes on the + side)
                int nd[4]; enodes(e, nd);
                bool plus[4]; int np = 0;
                for (int a = 0; a < 4; ++a) { double x, y; nodeXY(nd[a], x, y); plus[a] = (x - C.xc) * C.n[0] + (y - C.yc) * C.n[1] > 0; np += plus[a]; }
                if (np == 0 || np == 4) continue;   // segment only touches the element (grazing): no jump here
                const double R[2][2] = {{C.n[0], C.m[0]}, {C.n[1], C.m[1]}};   // [[u]] = R w, w = (w_n, w_t)
                const double Nn[2][3] = {{C.n[0], 0, C.n[1]}, {0, C.n[1], C.n[0]}};   // sigma n = Nn sigma_voigt
                double RtNn[2][3] = {};
                for (int r = 0; r < 2; ++r) for (int k = 0; k < 3; ++k) for (int s = 0; s < 2; ++s) RtNn[r][k] += R[s][r] * Nn[s][k];
                for (int q = 0; q < 4; ++q) {
                    double gx = 0, gy = 0;
                    for (int a = 0; a < 4; ++a) if (plus[a]) { gx += dNx[q][a]; gy += dNy[q][a]; }
                    const double BwG[3][2] = {{-gx, 0}, {0, -gy}, {-gy, -gx}};
                    double Bw[3][2] = {};
                    for (int r = 0; r < 3; ++r) for (int k = 0; k < 2; ++k) for (int s = 0; s < 2; ++s) Bw[r][k] += BwG[r][s] * R[s][k];
                    double DBw[3][2] = {}, DB[3][8] = {};
                    for (int r = 0; r < 3; ++r) for (int s = 0; s < 3; ++s) {
                        for (int k = 0; k < 2; ++k) DBw[r][k] += D[r][s] * Bw[s][k];
                        for (int k = 0; k < 8; ++k) DB[r][k] += D[r][s] * Bg[q][s][k];
                    }
                    for (int i2 = 0; i2 < 8; ++i2) for (int k = 0; k < 2; ++k) { double v = 0; for (int r = 0; r < 3; ++r) v += Bg[q][r][i2] * DBw[r][k]; C.Kuw[i2][k] += wq * v; }
                    // mean traction operator (1 / A) int R^T Nn D (B u + Bw w)
                    for (int r = 0; r < 2; ++r) {
                        for (int k = 0; k < 8; ++k) { double v = 0; for (int s = 0; s < 3; ++s) v += RtNn[r][s] * DB[s][k]; C.Gwu[r][k] += v / 4; }
                        for (int k = 0; k < 2; ++k) { double v = 0; for (int s = 0; s < 3; ++s) v += RtNn[r][s] * DBw[s][k]; C.Gww[r][k] += v / 4; }
                    }
                }
                cutOf[size_t(e)] = int(cuts.size());
                cuts.push_back(C);
            }
    }
    if (prev && prev->cut.size() == cuts.size())
        for (size_t k = 0; k < cuts.size(); ++k) { cuts[k].wtp = prev->cut[k].wtp; cuts[k].w[0] = prev->cut[k].wn; cuts[k].w[1] = prev->cut[k].wt; }
    // ---- loads (consistent nodal forces of the far-field tractions) and constraints
    std::vector<double> fext(size_t(N), 0.0);
    for (int j = 0; j < ny; ++j) { const int a = nx + nnx * j, b = a + nnx; fext[size_t(2 * a)] += 0.5 * in.sxx * hy; fext[size_t(2 * b)] += 0.5 * in.sxx * hy; }
    for (int i = 0; i < nx; ++i) { const int a = i + nnx * ny, b = a + 1; fext[size_t(2 * a + 1)] += 0.5 * in.syy * hx; fext[size_t(2 * b + 1)] += 0.5 * in.syy * hx; }
    std::vector<char> fixedDof(size_t(N), 0);
    for (int j = 0; j <= ny; ++j) fixedDof[size_t(2 * (nnx * j))] = 1;   // left: u_x = 0
    for (int i = 0; i <= nx; ++i) fixedDof[size_t(2 * i + 1)] = 1;      // bottom: u_y = 0
    double fScale = 0;
    for (double v : fext) fScale += v * v;
    fScale = std::sqrt(fScale);
    for (const auto& C : cuts) fScale = std::max(fScale, pressureMax(in.fracs[size_t(C.f)]) * C.len);
    fScale = std::max(fScale, 1e-30);
    // contact law
    auto contact = [&](Cut& C) {
        const Fracture& F = in.fracs[size_t(C.f)];
        const double wn = C.w[0], wt = C.w[1];
        C.T[0][0] = C.T[0][1] = C.T[1][0] = C.T[1][1] = 0;
        if (wn > 0) { C.t[0] = 0; C.t[1] = 0; C.state = 0; return; }
        const double tn = Kpen * wn, ttr = Kpen * (wt - C.wtp), lim = F.cohesion + F.friction * std::fabs(tn);
        C.t[0] = tn; C.T[0][0] = Kpen;
        if (std::fabs(ttr) <= lim) { C.t[1] = ttr; C.T[1][1] = Kpen; C.state = 1; }
        else { const double sg = ttr > 0 ? 1.0 : -1.0; C.t[1] = sg * lim; C.T[1][0] = -sg * F.friction * Kpen; C.state = 2; }
    };
    const int kb = 2 * (nnx + 1) + 1;
    out.u.assign(size_t(N), 0.0);
    if (prev && prev->u.size() == size_t(N)) out.u = prev->u;
    std::vector<double>& u = out.u;
    Band K;
    std::vector<double> r(static_cast<size_t>(N)), du;
    bool converged = false;
    std::vector<int> lastState(cuts.size(), -1);
    for (int it = 0; it < in.maxNewton; ++it) {
        K.init(N, kb, kb);
        std::fill(r.begin(), r.end(), 0.0);
        for (size_t i = 0; i < r.size(); ++i) r[i] = -fext[i];
        double rw = 0;
        bool stateChanged = false;
        for (int e = 0; e < nx * ny; ++e) {
            int nd[4]; enodes(e, nd);
            int dof[8]; for (int a = 0; a < 4; ++a) { dof[2 * a] = 2 * nd[a]; dof[2 * a + 1] = 2 * nd[a] + 1; }
            double ue[8]; for (int k = 0; k < 8; ++k) ue[k] = u[size_t(dof[k])];
            double Kc[8][8], re[8];
            for (int i = 0; i < 8; ++i) { re[i] = 0; for (int j = 0; j < 8; ++j) { Kc[i][j] = Ke[i][j]; re[i] += Ke[i][j] * ue[j]; } }
            const int ci = cutOf[size_t(e)];
            if (ci >= 0) {
                Cut& C = cuts[size_t(ci)];
                contact(C);
                if (C.state != lastState[size_t(ci)]) stateChanged = true;
                lastState[size_t(ci)] = C.state;
                const double p = pressureAt(in.fracs[size_t(C.f)], C.xc, C.yc);
                double rwv[2];
                for (int k = 0; k < 2; ++k) {
                    double v = C.Gww[k][0] * C.w[0] + C.Gww[k][1] * C.w[1] - C.t[k] + (k == 0 ? p : 0.0);
                    for (int j = 0; j < 8; ++j) v += C.Gwu[k][j] * ue[j];
                    rwv[k] = v;
                }
                rw = std::max(rw, C.len * std::max(std::fabs(rwv[0]), std::fabs(rwv[1])));
                for (int i = 0; i < 8; ++i) re[i] += C.Kuw[i][0] * C.w[0] + C.Kuw[i][1] * C.w[1];
                double Am[2][2] = {{C.Gww[0][0] - C.T[0][0], C.Gww[0][1] - C.T[0][1]}, {C.Gww[1][0] - C.T[1][0], C.Gww[1][1] - C.T[1][1]}}, Ai[2][2] = {{0, 0}, {0, 0}};
                if (!inv2(Am, Ai)) return fail("singular embedded-jump equation in a cut element");
                double KA[8][2] = {};
                for (int i = 0; i < 8; ++i) for (int k = 0; k < 2; ++k) KA[i][k] = C.Kuw[i][0] * Ai[0][k] + C.Kuw[i][1] * Ai[1][k];
                for (int i = 0; i < 8; ++i) {
                    for (int j = 0; j < 8; ++j) Kc[i][j] -= KA[i][0] * C.Gwu[0][j] + KA[i][1] * C.Gwu[1][j];
                    re[i] -= KA[i][0] * rwv[0] + KA[i][1] * rwv[1];   // condensed residual
                }
            }
            for (int i = 0; i < 8; ++i) {
                r[size_t(dof[i])] += re[i];
                for (int j = 0; j < 8; ++j) K.at(dof[i], dof[j]) += Kc[i][j];
            }
        }
        // residual norm (true, uncondensed u-equation residual is r without the condensation terms: use both checks)
        double rn = 0;
        for (int i = 0; i < N; ++i) if (!fixedDof[size_t(i)]) rn += r[size_t(i)] * r[size_t(i)];
        rn = std::sqrt(rn);
        out.newton = it;
        if (it > 0 && !stateChanged && rn <= in.tol * fScale && rw <= in.tol * fScale) { converged = true; break; }
        for (int i = 0; i < N; ++i)
            if (fixedDof[size_t(i)]) {
                for (int j = std::max(0, i - kb); j <= std::min(N - 1, i + kb); ++j) { K.at(i, j) = 0; K.at(j, i) = 0; }
                K.at(i, i) = 1; r[size_t(i)] = 0;
            }
        if (!K.factor()) return fail("singular stiffness matrix (check the boundary conditions)");
        du = r; for (double& v : du) v = -v;
        K.solve(du);
        for (int i = 0; i < N; ++i) u[size_t(i)] += du[size_t(i)];
        // jump updates: dw = -A^{-1} (r_w + G_wu du), with r_w at the start of the iteration
        for (Cut& C : cuts) {
            int nd[4]; enodes(C.e, nd);
            double ue0[8], due[8];
            for (int a = 0; a < 4; ++a) for (int k = 0; k < 2; ++k) { due[2 * a + k] = du[size_t(2 * nd[a] + k)]; ue0[2 * a + k] = u[size_t(2 * nd[a] + k)] - due[2 * a + k]; }
            const double p = pressureAt(in.fracs[size_t(C.f)], C.xc, C.yc);
            double rwv[2];
            for (int k = 0; k < 2; ++k) {
                double v = C.Gww[k][0] * C.w[0] + C.Gww[k][1] * C.w[1] - C.t[k] + (k == 0 ? p : 0.0);
                for (int j = 0; j < 8; ++j) v += C.Gwu[k][j] * (ue0[j] + due[j]);
                rwv[k] = v;
            }
            double Am[2][2] = {{C.Gww[0][0] - C.T[0][0], C.Gww[0][1] - C.T[0][1]}, {C.Gww[1][0] - C.T[1][0], C.Gww[1][1] - C.T[1][1]}}, Ai[2][2] = {{0, 0}, {0, 0}};
            inv2(Am, Ai);
            C.w[0] -= Ai[0][0] * rwv[0] + Ai[0][1] * rwv[1];
            C.w[1] -= Ai[1][0] * rwv[0] + Ai[1][1] * rwv[1];
        }
    }
    if (!converged) return fail("EFEM contact iteration did not converge");
    // ---- output
    out.dofs = N;
    for (Cut& C : cuts) {
        contact(C);
        if (C.state == 2) C.wtp = C.w[1] - C.t[1] / Kpen;   // commit the plastic slip
        CutElement o; o.elem = C.e; o.frac = C.f; o.xc = C.xc; o.yc = C.yc; o.s = C.s; o.len = C.len;
        o.nx = C.n[0]; o.ny = C.n[1]; o.wn = C.w[0]; o.wt = C.w[1]; o.tn = C.t[0]; o.tt = C.t[1]; o.state = C.state; o.wtp = C.wtp;
        out.cut.push_back(o);
    }
    out.sxx.assign(size_t(nx) * ny, 0.0); out.syy = out.sxx; out.sxy = out.sxx;
    for (int e = 0; e < nx * ny; ++e) {
        int nd[4]; enodes(e, nd);
        double ue[8]; for (int a = 0; a < 4; ++a) { ue[2 * a] = u[size_t(2 * nd[a])]; ue[2 * a + 1] = u[size_t(2 * nd[a] + 1)]; }
        double sm[3] = {0, 0, 0};
        for (int q = 0; q < 4; ++q) {
            double eps[3] = {0, 0, 0};
            for (int rr = 0; rr < 3; ++rr) for (int k = 0; k < 8; ++k) eps[rr] += Bg[q][rr][k] * ue[k];
            const int ci = cutOf[size_t(e)];
            if (ci >= 0) {   // enhanced strain
                const Cut& C = cuts[size_t(ci)];
                bool plus[4];
                for (int a = 0; a < 4; ++a) { double x, y; nodeXY(nd[a], x, y); plus[a] = (x - C.xc) * C.n[0] + (y - C.yc) * C.n[1] > 0; }
                double gx = 0, gy = 0;
                for (int a = 0; a < 4; ++a) if (plus[a]) { gx += dNx[q][a]; gy += dNy[q][a]; }
                const double Jx = C.n[0] * C.w[0] + C.m[0] * C.w[1], Jy = C.n[1] * C.w[0] + C.m[1] * C.w[1];
                eps[0] -= gx * Jx; eps[1] -= gy * Jy; eps[2] -= gx * Jy + gy * Jx;
            }
            for (int rr = 0; rr < 3; ++rr) for (int s = 0; s < 3; ++s) sm[rr] += D[rr][s] * eps[s] / 4;
        }
        out.sxx[size_t(e)] = sm[0]; out.syy[size_t(e)] = sm[1]; out.sxy[size_t(e)] = sm[2];
    }
    out.ok = true;
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

}  // namespace efem
}  // namespace ressim

// ============================================================== continuous jumps (nodal Heaviside enrichment, default)
namespace ressim {
namespace efem {
namespace {
struct Poly { std::vector<double> x, y; };
// part of a convex polygon on the side s * ((x - px) nx + (y - py) ny) >= 0
Poly clipHalf(const Poly& P, double px, double py, double nx, double ny, double s) {
    Poly Q;
    const size_t n = P.x.size();
    for (size_t i = 0; i < n; ++i) {
        const size_t j = (i + 1) % n;
        const double di = s * ((P.x[i] - px) * nx + (P.y[i] - py) * ny), dj = s * ((P.x[j] - px) * nx + (P.y[j] - py) * ny);
        if (di >= 0) { Q.x.push_back(P.x[i]); Q.y.push_back(P.y[i]); }
        if ((di >= 0) != (dj >= 0)) {
            const double t = di / (di - dj);
            Q.x.push_back(P.x[i] + t * (P.x[j] - P.x[i])); Q.y.push_back(P.y[i] + t * (P.y[j] - P.y[i]));
        }
    }
    return Q;
}
}  // namespace

static bool solveContinuous(const Input& in, Result& out, const Result* prev) {
    const auto t0 = std::chrono::steady_clock::now();
    out = Result();
    auto fail = [&](const std::string& m) { out.error = m; return false; };
    if (in.nx < 2 || in.ny < 2 || !(in.Lx > 0) || !(in.Ly > 0)) return fail("mesh / domain out of range");
    if (!(in.E > 0) || !(in.nu >= 0 && in.nu < 0.5)) return fail("elastic constants out of range");
    const int nx = in.nx, ny = in.ny, nnx = nx + 1, nn = nnx * (ny + 1);
    const double hx = in.Lx / nx, hy = in.Ly / ny;
    const double c = in.E / ((1 + in.nu) * (1 - 2 * in.nu));
    const double D[3][3] = {{c * (1 - in.nu), c * in.nu, 0}, {c * in.nu, c * (1 - in.nu), 0}, {0, 0, c * (1 - 2 * in.nu) / 2}};
    const double Kpen = in.penalty * in.E / std::min(hx, hy);
    const double xiA[4] = {-1, 1, 1, -1}, etA[4] = {-1, -1, 1, 1};
    auto enodes = [&](int e, int nd[4]) { const int i = e % nx, j = e / nx; nd[0] = i + nnx * j; nd[1] = nd[0] + 1; nd[2] = nd[1] + nnx; nd[3] = nd[0] + nnx; };
    auto nodeX = [&](int a) { return (a % nnx) * hx; };
    auto nodeY = [&](int a) { return (a / nnx) * hy; };
    // shape functions / gradients of element e at (x, y)
    auto shape = [&](int e, double x, double y, double N[4], double gx[4], double gy[4]) {
        const double xc = (e % nx + 0.5) * hx, yc = (e / nx + 0.5) * hy, xi = 2 * (x - xc) / hx, et = 2 * (y - yc) / hy;
        for (int a = 0; a < 4; ++a) {
            N[a] = (1 + xi * xiA[a]) * (1 + et * etA[a]) / 4;
            gx[a] = xiA[a] * (1 + et * etA[a]) / 4 * 2 / hx;
            gy[a] = etA[a] * (1 + xi * xiA[a]) / 4 * 2 / hy;
        }
    };
    // ---- fractures: cut elements, enriched nodes
    struct CutC { int e = 0, f = 0; double xa = 0, ya = 0, xb = 0, yb = 0, len = 0, s = 0; double wtp[2] = {0, 0}; };
    std::vector<CutC> cuts;
    std::vector<int> cutOf(static_cast<size_t>(nx) * ny, -1), enrF(static_cast<size_t>(nn), -1);
    struct Frac { double nx = 0, ny = 0, mx = 0, my = 0, L = 0; };
    std::vector<Frac> FR(in.fracs.size());
    for (size_t f = 0; f < in.fracs.size(); ++f) {
        const Fracture& F = in.fracs[f];
        const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
        if (!(L > 0)) return fail("fracture of zero length");
        FR[f].L = L; FR[f].mx = (F.x1 - F.x0) / L; FR[f].my = (F.y1 - F.y0) / L; FR[f].nx = -FR[f].my; FR[f].ny = FR[f].mx;
        const int i0 = std::max(0, int(std::floor(std::min(F.x0, F.x1) / hx)) - 1), i1 = std::min(nx - 1, int(std::floor(std::max(F.x0, F.x1) / hx)) + 1);
        const int j0 = std::max(0, int(std::floor(std::min(F.y0, F.y1) / hy)) - 1), j1 = std::min(ny - 1, int(std::floor(std::max(F.y0, F.y1) / hy)) + 1);
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) {
                double ta, tb;
                if (!clip(F.x0, F.y0, F.x1, F.y1, i * hx, j * hy, (i + 1) * hx, (j + 1) * hy, ta, tb)) continue;
                const double len = (tb - ta) * L;
                if (len < 1e-6 * std::min(hx, hy)) continue;
                const int e = i + nx * j;
                int nd[4]; enodes(e, nd); int np = 0;
                for (int a = 0; a < 4; ++a) np += (nodeX(nd[a]) - F.x0) * FR[f].nx + (nodeY(nd[a]) - F.y0) * FR[f].ny > 0;
                if (np == 0 || np == 4) continue;
                if (cutOf[size_t(e)] >= 0) return fail("two fractures cut the same element");
                CutC C; C.e = e; C.f = int(f); C.len = len;
                C.xa = F.x0 + ta * (F.x1 - F.x0); C.ya = F.y0 + ta * (F.y1 - F.y0); C.xb = F.x0 + tb * (F.x1 - F.x0); C.yb = F.y0 + tb * (F.y1 - F.y0);
                C.s = (0.5 * (ta + tb) - 0.5) * L;
                cutOf[size_t(e)] = int(cuts.size());
                cuts.push_back(C);
            }
        // tips (end points inside the domain): nodes whose (open) support contains a tip are not enriched -> zero jump there
        auto tipInSupport = [&](int a) {
            for (int k = 0; k < 2; ++k) {
                const double tx = k ? F.x1 : F.x0, ty = k ? F.y1 : F.y0;
                if (tx <= 1e-12 || ty <= 1e-12 || tx >= in.Lx - 1e-12 || ty >= in.Ly - 1e-12) continue;   // edge crack end
                if (std::fabs(tx - nodeX(a)) < hx * (1 - 1e-9) && std::fabs(ty - nodeY(a)) < hy * (1 - 1e-9)) return true;
            }
            return false;
        };
        for (const CutC& C : cuts) {
            if (C.f != int(f)) continue;
            int nd[4]; enodes(C.e, nd);
            for (int a = 0; a < 4; ++a) {
                if (tipInSupport(nd[a])) continue;
                if (enrF[size_t(nd[a])] >= 0 && enrF[size_t(nd[a])] != int(f)) return fail("two fractures closer than one element (a node would carry two enrichments)");
                enrF[size_t(nd[a])] = int(f);
            }
        }
    }
    if (prev && prev->cut.size() == cuts.size())
        for (size_t k = 0; k < cuts.size(); ++k) cuts[k].wtp[0] = cuts[k].wtp[1] = prev->cut[k].wtp;
    // ---- dof numbering: per node 2 regular (+ 2 enriched)
    std::vector<int> dof0(static_cast<size_t>(nn)), dofE(static_cast<size_t>(nn), -1);
    int N = 0;
    for (int a = 0; a < nn; ++a) { dof0[size_t(a)] = N; N += 2; if (enrF[size_t(a)] >= 0) { dofE[size_t(a)] = N; N += 2; } }
    auto Hside = [&](int f, double x, double y) { return (x - in.fracs[size_t(f)].x0) * FR[size_t(f)].nx + (y - in.fracs[size_t(f)].y0) * FR[size_t(f)].ny > 0 ? 1.0 : 0.0; };
    // element dof list
    auto edofs = [&](int e, std::vector<int>& d, int nd[4]) {
        enodes(e, nd); d.clear();
        for (int a = 0; a < 4; ++a) { d.push_back(dof0[size_t(nd[a])]); d.push_back(dof0[size_t(nd[a])] + 1); }
        for (int a = 0; a < 4; ++a) if (dofE[size_t(nd[a])] >= 0) { d.push_back(dofE[size_t(nd[a])]); d.push_back(dofE[size_t(nd[a])] + 1); }
    };
    int kb = 1;
    {
        std::vector<int> d; int nd[4];
        for (int e = 0; e < nx * ny; ++e) { edofs(e, d, nd); const auto mm = std::minmax_element(d.begin(), d.end()); kb = std::max(kb, *mm.second - *mm.first); }
    }
    // quadrature points of an element: (x, y, weight, H per fracture-relevant evaluation done by the caller)
    struct QP { double x, y, w; };
    auto quad = [&](int e, std::vector<QP>& qp) {
        qp.clear();
        const int i = e % nx, j = e / nx;
        const int ci = cutOf[size_t(e)];
        if (ci < 0) {
            const double g = 1 / std::sqrt(3.0), xc = (i + 0.5) * hx, yc = (j + 0.5) * hy;
            for (int q = 0; q < 4; ++q) qp.push_back({xc + ((q == 1 || q == 2) ? g : -g) * hx / 2, yc + (q >= 2 ? g : -g) * hy / 2, hx * hy / 4});
            return;
        }
        const CutC& C = cuts[size_t(ci)];
        Poly R; R.x = {i * hx, (i + 1) * hx, (i + 1) * hx, i * hx}; R.y = {j * hy, j * hy, (j + 1) * hy, (j + 1) * hy};
        for (double sg : {1.0, -1.0}) {
            const Poly P = clipHalf(R, C.xa, C.ya, FR[size_t(C.f)].nx, FR[size_t(C.f)].ny, sg);
            for (size_t k = 1; k + 1 < P.x.size(); ++k) {   // fan triangulation, 3-point rule (degree 2)
                const double x0 = P.x[0], y0 = P.y[0], x1 = P.x[k], y1 = P.y[k], x2 = P.x[k + 1], y2 = P.y[k + 1];
                const double ar = 0.5 * std::fabs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
                if (ar <= 0) continue;
                const double L[3][3] = {{2. / 3, 1. / 6, 1. / 6}, {1. / 6, 2. / 3, 1. / 6}, {1. / 6, 1. / 6, 2. / 3}};
                for (int q = 0; q < 3; ++q) qp.push_back({L[q][0] * x0 + L[q][1] * x1 + L[q][2] * x2, L[q][0] * y0 + L[q][1] * y1 + L[q][2] * y2, ar / 3});
            }
        }
    };
    // extended B matrix at a point (3 x ndof of the element)
    auto bext = [&](int e, const int nd[4], double x, double y, std::vector<double>& B, int nde) {
        double Nn[4], gx[4], gy[4]; shape(e, x, y, Nn, gx, gy);
        B.assign(size_t(3) * nde, 0.0);
        auto put = [&](int col, double dx, double dy) { B[size_t(col)] = dx; B[size_t(nde + col + 1)] = dy; B[size_t(2 * nde + col)] = dy; B[size_t(2 * nde + col + 1)] = dx; };
        for (int a = 0; a < 4; ++a) put(2 * a, gx[a], gy[a]);
        int k = 8;
        for (int a = 0; a < 4; ++a) {
            const int f = enrF[size_t(nd[a])];
            if (f < 0) continue;
            const double sh = Hside(f, x, y) - Hside(f, nodeX(nd[a]), nodeY(nd[a]));
            put(k, sh * gx[a], sh * gy[a]);
            k += 2;
        }
    };
    // ---- loads, constraints
    std::vector<double> fext(static_cast<size_t>(N), 0.0);
    for (int j = 0; j < ny; ++j) { const int a = nx + nnx * j, b = a + nnx; fext[size_t(dof0[size_t(a)])] += 0.5 * in.sxx * hy; fext[size_t(dof0[size_t(b)])] += 0.5 * in.sxx * hy; }
    for (int i = 0; i < nx; ++i) { const int a = i + nnx * ny, b = a + 1; fext[size_t(dof0[size_t(a)] + 1)] += 0.5 * in.syy * hx; fext[size_t(dof0[size_t(b)] + 1)] += 0.5 * in.syy * hx; }
    std::vector<char> fixedDof(static_cast<size_t>(N), 0);
    for (int j = 0; j <= ny; ++j) fixedDof[size_t(dof0[size_t(nnx * j)])] = 1;
    for (int i = 0; i <= nx; ++i) fixedDof[size_t(dof0[size_t(i)] + 1)] = 1;
    double fScale = 0;
    for (double v : fext) fScale += v * v;
    fScale = std::sqrt(fScale);
    for (const auto& C : cuts) fScale = std::max(fScale, pressureMax(in.fracs[size_t(C.f)]) * C.len);
    fScale = std::max(fScale, 1e-30);
    // element stiffness matrices (constant): cache per element with enrichment, the regular one shared
    std::vector<std::vector<double>> Kel(static_cast<size_t>(nx) * ny);
    std::vector<double> Kreg;
    {
        std::vector<int> d; int nd[4]; std::vector<QP> qp; std::vector<double> B;
        for (int e = 0; e < nx * ny; ++e) {
            edofs(e, d, nd);
            const int nde = int(d.size());
            if (nde == 8 && !Kreg.empty()) continue;
            std::vector<double> Ke(size_t(nde) * nde, 0.0);
            quad(e, qp);
            for (const QP& q : qp) {
                bext(e, nd, q.x, q.y, B, nde);
                std::vector<double> DB(size_t(3) * nde, 0.0);
                for (int r = 0; r < 3; ++r) for (int s = 0; s < 3; ++s) { const double v = D[r][s]; if (v != 0) for (int k = 0; k < nde; ++k) DB[size_t(r) * nde + k] += v * B[size_t(s) * nde + k]; }
                for (int i = 0; i < nde; ++i) for (int j = 0; j < nde; ++j) { double v = 0; for (int r = 0; r < 3; ++r) v += B[size_t(r) * nde + i] * DB[size_t(r) * nde + j]; Ke[size_t(i) * nde + j] += q.w * v; }
            }
            if (nde == 8) Kreg = Ke; else Kel[size_t(e)] = Ke;
        }
    }
    // ---- Newton / active set on the contact states
    out.u.assign(static_cast<size_t>(N), 0.0);
    std::vector<double>& U = out.u;
    if (prev && prev->u.size() == size_t(N)) U = prev->u;
    const double gpS[2] = {0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0)};
    std::vector<int> lastState(cuts.size() * 2, -1), curState(cuts.size() * 2, 0);
    std::vector<double> gw(cuts.size() * 4, 0.0), gt(cuts.size() * 4, 0.0);   // per Gauss point (wn, wt) and (tn, tt)
    Band K;
    std::vector<double> r(static_cast<size_t>(N)), du;
    bool converged = false;
    for (int it = 0; it < in.maxNewton; ++it) {
        K.init(N, kb, kb);
        for (int i = 0; i < N; ++i) r[size_t(i)] = -fext[size_t(i)];
        std::vector<int> d; int nd[4];
        for (int e = 0; e < nx * ny; ++e) {
            edofs(e, d, nd);
            const int nde = int(d.size());
            const std::vector<double>& Ke = nde == 8 ? Kreg : Kel[size_t(e)];
            for (int i = 0; i < nde; ++i) {
                double s = 0;
                for (int j = 0; j < nde; ++j) { const double v = Ke[size_t(i) * nde + j]; s += v * U[size_t(d[size_t(j)])]; K.at(d[size_t(i)], d[size_t(j)]) += v; }
                r[size_t(d[size_t(i)])] += s;
            }
        }
        bool stateChanged = false;
        for (size_t ci = 0; ci < cuts.size(); ++ci) {
            CutC& C = cuts[ci];
            const Fracture& F = in.fracs[size_t(C.f)];
            const Frac& G = FR[size_t(C.f)];
            enodes(C.e, nd);
            for (int g = 0; g < 2; ++g) {
                const double x = C.xa + gpS[g] * (C.xb - C.xa), y = C.ya + gpS[g] * (C.yb - C.ya), wl = 0.5 * C.len;
                double Nn[4], gx[4], gy[4]; shape(C.e, x, y, Nn, gx, gy);
                double J[2] = {0, 0};
                for (int a = 0; a < 4; ++a) if (dofE[size_t(nd[a])] >= 0) { J[0] += Nn[a] * U[size_t(dofE[size_t(nd[a])])]; J[1] += Nn[a] * U[size_t(dofE[size_t(nd[a])] + 1)]; }
                const double wn = J[0] * G.nx + J[1] * G.ny, wt = J[0] * G.mx + J[1] * G.my;
                double t[2] = {0, 0}, T[2][2] = {{0, 0}, {0, 0}};
                int st = 0;
                if (wn <= 0) {   // closed (w_n = 0 included: the first iteration starts from the stick tangent)
                    const double tn = Kpen * wn, ttr = Kpen * (wt - C.wtp[g]), lim = F.cohesion + F.friction * std::fabs(tn);
                    t[0] = tn; T[0][0] = Kpen;
                    if (std::fabs(ttr) <= lim) { t[1] = ttr; T[1][1] = Kpen; st = 1; }
                    else { const double sg = ttr > 0 ? 1.0 : -1.0; t[1] = sg * lim; T[1][0] = -sg * F.friction * Kpen; st = 2; }
                }
                const size_t gi = ci * 2 + size_t(g);
                curState[gi] = st;
                if (st != lastState[gi]) stateChanged = true;
                gw[gi * 2] = wn; gw[gi * 2 + 1] = wt; gt[gi * 2] = t[0]; gt[gi * 2 + 1] = t[1];
                // global traction vector and tangent: t_glob = R (t - p e_n), Tg = R T R^T  (R = [n m])
                const double tl0 = t[0] - pressureAt(F, x, y), tl1 = t[1];
                const double tg[2] = {G.nx * tl0 + G.mx * tl1, G.ny * tl0 + G.my * tl1};
                const double Rm[2][2] = {{G.nx, G.mx}, {G.ny, G.my}};
                double Tg[2][2] = {{0, 0}, {0, 0}};
                for (int p = 0; p < 2; ++p) for (int q = 0; q < 2; ++q) for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) Tg[p][q] += Rm[p][k] * T[k][l] * Rm[q][l];
                for (int a = 0; a < 4; ++a) {
                    const int da = dofE[size_t(nd[a])];
                    if (da < 0) continue;
                    for (int p = 0; p < 2; ++p) r[size_t(da + p)] += wl * Nn[a] * tg[p];
                    for (int b = 0; b < 4; ++b) {
                        const int db = dofE[size_t(nd[b])];
                        if (db < 0) continue;
                        for (int p = 0; p < 2; ++p) for (int q = 0; q < 2; ++q) K.at(da + p, db + q) += wl * Nn[a] * Nn[b] * Tg[p][q];
                    }
                }
            }
        }
        lastState = curState;
        double rn = 0;
        for (int i = 0; i < N; ++i) if (!fixedDof[size_t(i)]) rn += r[size_t(i)] * r[size_t(i)];
        rn = std::sqrt(rn);
        out.newton = it;
        // the residual uses the return-mapped tractions of the current displacements: equilibrium is reached when it is
        // small, whatever the state flips of Gauss points with w_n ~ 0 at the tips (open <-> stick chatter)
        if (std::getenv("EFEM_DEBUG")) {
            int ns[3] = {0, 0, 0}; for (int v : curState) ns[v]++;
            std::fprintf(stderr, "  it %d  |r| %.3e (tol %.3e)  open %d stick %d slip %d  changed %d\n", it, rn, in.tol * fScale, ns[0], ns[1], ns[2], int(stateChanged));
        }
        if (it > 0 && rn <= in.tol * fScale) { converged = true; break; }
        for (int i = 0; i < N; ++i)
            if (fixedDof[size_t(i)]) {
                for (int j = std::max(0, i - kb); j <= std::min(N - 1, i + kb); ++j) { K.at(i, j) = 0; K.at(j, i) = 0; }
                K.at(i, i) = 1; r[size_t(i)] = 0;
            }
        if (!K.factor()) return fail("singular stiffness matrix (check the boundary conditions / enrichment)");
        du = r; for (double& v : du) v = -v;
        K.solve(du);
        for (int i = 0; i < N; ++i) U[size_t(i)] += du[size_t(i)];
    }
    if (!converged) return fail("contact iteration did not converge");
    out.dofs = N;
    for (size_t ci = 0; ci < cuts.size(); ++ci) {
        CutC& C = cuts[ci];
        CutElement o; o.elem = C.e; o.frac = C.f; o.len = C.len; o.s = C.s;
        o.xc = 0.5 * (C.xa + C.xb); o.yc = 0.5 * (C.ya + C.yb); o.nx = FR[size_t(C.f)].nx; o.ny = FR[size_t(C.f)].ny;
        o.wn = 0.5 * (gw[ci * 4] + gw[ci * 4 + 2]); o.wt = 0.5 * (gw[ci * 4 + 1] + gw[ci * 4 + 3]);
        o.tn = 0.5 * (gt[ci * 4] + gt[ci * 4 + 2]); o.tt = 0.5 * (gt[ci * 4 + 1] + gt[ci * 4 + 3]);
        o.state = std::max(curState[ci * 2], curState[ci * 2 + 1]);
        for (int g = 0; g < 2; ++g) if (curState[ci * 2 + size_t(g)] == 2) C.wtp[g] = gw[(ci * 2 + size_t(g)) * 2 + 1] - gt[(ci * 2 + size_t(g)) * 2 + 1] / Kpen;
        o.wtp = 0.5 * (C.wtp[0] + C.wtp[1]);
        out.cut.push_back(o);
    }
    out.sxx.assign(static_cast<size_t>(nx) * ny, 0.0); out.syy = out.sxx; out.sxy = out.sxx;
    {
        std::vector<int> d; int nd[4]; std::vector<QP> qp; std::vector<double> B;
        for (int e = 0; e < nx * ny; ++e) {
            edofs(e, d, nd); const int nde = int(d.size());
            quad(e, qp);
            double sm[3] = {0, 0, 0}, wsum = 0;
            for (const QP& q : qp) {
                bext(e, nd, q.x, q.y, B, nde);
                double eps[3] = {0, 0, 0};
                for (int rr = 0; rr < 3; ++rr) for (int k = 0; k < nde; ++k) eps[rr] += B[size_t(rr) * nde + k] * U[size_t(d[size_t(k)])];
                for (int rr = 0; rr < 3; ++rr) for (int s = 0; s < 3; ++s) sm[rr] += q.w * D[rr][s] * eps[s];
                wsum += q.w;
            }
            out.sxx[size_t(e)] = sm[0] / wsum; out.syy[size_t(e)] = sm[1] / wsum; out.sxy[size_t(e)] = sm[2] / wsum;
        }
    }
    out.ok = true;
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

bool solve(const Input& in, Result& out, const Result* prev) {
    return in.jumpModel == 1 ? solveConstantJump(in, out, prev) : solveContinuous(in, out, prev);
}

}  // namespace efem
}  // namespace ressim
