// Linear thermo-poroelastic finite elements (hex8) with fixed-stress coupling coefficients - see geomech.h.
#include "ressim/geomech.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace ressim {
namespace geomech {

namespace {
constexpr int kLoc[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// shape-function gradients of a brick (sizes a, b, c) at the natural point (xi, eta, zeta): G[node][dir]
void gradN(double a, double b, double c, double xi, double eta, double zeta, double G[8][3]) {
    for (int n = 0; n < 8; ++n) {
        const double sx = 2 * kLoc[n][0] - 1, sy = 2 * kLoc[n][1] - 1, sz = 2 * kLoc[n][2] - 1;
        G[n][0] = 0.125 * sx * (1 + sy * eta) * (1 + sz * zeta) * 2 / a;
        G[n][1] = 0.125 * sy * (1 + sx * xi) * (1 + sz * zeta) * 2 / b;
        G[n][2] = 0.125 * sz * (1 + sx * xi) * (1 + sy * eta) * 2 / c;
    }
}
// strain-displacement matrix (6 x 24; engineering shear strains xy, yz, xz)
void bMatrix(const double G[8][3], double B[6][24]) {
    std::memset(B, 0, sizeof(double) * 6 * 24);
    for (int n = 0; n < 8; ++n) {
        const int c = 3 * n;
        B[0][c] = G[n][0]; B[1][c + 1] = G[n][1]; B[2][c + 2] = G[n][2];
        B[3][c] = G[n][1]; B[3][c + 1] = G[n][0];
        B[4][c + 1] = G[n][2]; B[4][c + 2] = G[n][1];
        B[5][c] = G[n][2]; B[5][c + 2] = G[n][0];
    }
}
void dMatrix(double E, double nu, double D[6][6]) {
    std::memset(D, 0, sizeof(double) * 36);
    const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) D[i][j] = lam; D[i][i] = lam + 2 * mu; D[i + 3][i + 3] = mu; }
}
// eigenvalues of a symmetric 3x3 matrix (descending)
void eig3(const double s[6], double& l1, double& l2, double& l3) {
    // s = xx, yy, zz, xy, yz, xz
    const double p1 = s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    const double q = (s[0] + s[1] + s[2]) / 3;
    if (p1 < 1e-30 * (1 + q * q)) {
        double e[3] = {s[0], s[1], s[2]}; std::sort(e, e + 3);
        l1 = e[2]; l2 = e[1]; l3 = e[0]; return;
    }
    const double p2 = (s[0] - q) * (s[0] - q) + (s[1] - q) * (s[1] - q) + (s[2] - q) * (s[2] - q) + 2 * p1;
    const double p = std::sqrt(p2 / 6);
    const double b00 = (s[0] - q) / p, b11 = (s[1] - q) / p, b22 = (s[2] - q) / p, b01 = s[3] / p, b12 = s[4] / p, b02 = s[5] / p;
    const double detB = b00 * (b11 * b22 - b12 * b12) - b01 * (b01 * b22 - b12 * b02) + b02 * (b01 * b12 - b11 * b02);
    const double r = std::clamp(detB / 2, -1.0, 1.0);
    const double phi = std::acos(r) / 3;
    l1 = q + 2 * p * std::cos(phi);
    l3 = q + 2 * p * std::cos(phi + 2.0943951023931953);
    l2 = 3 * q - l1 - l3;
}
template <class T> T at(const std::vector<T>& v, size_t i, T d) { return v.empty() ? d : v.size() == 1 ? v[0] : (i < v.size() ? v[i] : d); }
} // namespace

void Solver::elemNodes(int i, int j, int k, int nd[8]) const {
    for (int n = 0; n < 8; ++n) nd[n] = node(i + kLoc[n][0], j + kLoc[n][1], k + kLoc[n][2]);
}

bool Solver::setup(const Input& in, int nx, int ny, int nz, const std::vector<double>& dx, const std::vector<double>& dy,
                   const std::vector<double>& dz, double top, const std::vector<double>& p0, const std::vector<double>& phi0,
                   std::string* error) {
    return setupImpl(in, nx, ny, nz, dx, dy, dz, top, p0, phi0, error, nullptr);
}

bool Solver::setupHex(const Input& in, int nx, int ny, int nz, const std::vector<double>& nodes, const std::vector<double>& p0,
                      const std::vector<double>& phi0, std::string* error) {
    const size_t nn = size_t(nx + 1) * size_t(ny + 1) * size_t(nz + 1);
    if (nx <= 0 || ny <= 0 || nz <= 0 || nodes.size() != 3 * nn) { if (error) *error = "geomechanics: corner-point node array size"; ready_ = false; return false; }
    auto nd = [&](int i, int j, int k) { return 3 * ((size_t(k) * size_t(ny + 1) + size_t(j)) * size_t(nx + 1) + size_t(i)); };
    // representative tensor sizes (used for padding thicknesses, face areas and loads): averages over the box
    std::vector<double> dx(size_t(nx), 0.0), dy(size_t(ny), 0.0), dz(size_t(nz), 0.0);
    for (int i = 0; i < nx; ++i) { double s = 0; int c = 0; for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) { s += std::hypot(nodes[nd(i + 1, j, k)] - nodes[nd(i, j, k)], nodes[nd(i + 1, j, k) + 1] - nodes[nd(i, j, k) + 1]); ++c; } dx[size_t(i)] = s / c; }
    for (int j = 0; j < ny; ++j) { double s = 0; int c = 0; for (int k = 0; k <= nz; ++k) for (int i = 0; i <= nx; ++i) { s += std::hypot(nodes[nd(i, j + 1, k)] - nodes[nd(i, j, k)], nodes[nd(i, j + 1, k) + 1] - nodes[nd(i, j, k) + 1]); ++c; } dy[size_t(j)] = s / c; }
    for (int k = 0; k < nz; ++k) { double s = 0; int c = 0; for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) { s += nodes[nd(i, j, k + 1) + 2] - nodes[nd(i, j, k) + 2]; ++c; } dz[size_t(k)] = std::max(1e-3, s / c); }
    // per-cell arrays expected by setupImpl (cell (i, 0, 0) / (0, j, 0) / (0, 0, k) carry the sizes)
    const size_t n = size_t(nx) * size_t(ny) * size_t(nz);
    std::vector<double> DX(n), DY(n), DZ(n);
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const size_t c = (size_t(k) * size_t(ny) + size_t(j)) * size_t(nx) + size_t(i);
        DX[c] = dx[size_t(i)]; DY[c] = dy[size_t(j)]; DZ[c] = dz[size_t(k)];
    }
    double top = 1e300; for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) top = std::min(top, nodes[nd(i, j, 0) + 2]);
    return setupImpl(in, nx, ny, nz, DX, DY, DZ, top, p0, phi0, error, &nodes);
}

double Solver::elemGrad(int i, int j, int k, double xi, double eta, double zeta, double G[8][3]) const {
    if (X_.empty()) {
        const double a = ex_[size_t(i)], b = ey_[size_t(j)], c = ez_[size_t(k)];
        gradN(a, b, c, xi, eta, zeta, G);
        return a * b * c / 8;
    }
    // isoparametric trilinear hexahedron: dN/dxi, J = sum x_n (dN_n/dxi)^T, G = dN/dxi J^-1
    double dN[8][3];
    for (int n = 0; n < 8; ++n) {
        const double sx = 2 * kLoc[n][0] - 1, sy = 2 * kLoc[n][1] - 1, sz = 2 * kLoc[n][2] - 1;
        dN[n][0] = 0.125 * sx * (1 + sy * eta) * (1 + sz * zeta);
        dN[n][1] = 0.125 * sy * (1 + sx * xi) * (1 + sz * zeta);
        dN[n][2] = 0.125 * sz * (1 + sx * xi) * (1 + sy * eta);
    }
    int nd[8]; elemNodes(i, j, k, nd);
    double J[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};   // J[d][r] = d x_d / d xi_r
    for (int n = 0; n < 8; ++n) for (int d = 0; d < 3; ++d) for (int r = 0; r < 3; ++r) J[d][r] += X_[3 * size_t(nd[n]) + size_t(d)] * dN[n][r];
    const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    const double id = std::fabs(det) > 1e-300 ? 1.0 / det : 0.0;
    double Ji[3][3];   // inverse: Ji[r][d] = d xi_r / d x_d
    Ji[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * id; Ji[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * id; Ji[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * id;
    Ji[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * id; Ji[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * id; Ji[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * id;
    Ji[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * id; Ji[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * id; Ji[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * id;
    for (int n = 0; n < 8; ++n) for (int d = 0; d < 3; ++d) G[n][d] = dN[n][0] * Ji[0][d] + dN[n][1] * Ji[1][d] + dN[n][2] * Ji[2][d];
    return std::fabs(det);
}

void Solver::totalStress(int cell, double S[6]) const {
    for (int r = 0; r < 6; ++r) S[r] = 0;
    if (cell < 0 || 6 * size_t(cell) + 5 >= sfull_.size()) return;
    for (int r = 0; r < 6; ++r) S[r] = sfull_[6 * size_t(cell) + size_t(r)];
    const double p = size_t(cell) < pAbs_.size() ? pAbs_[size_t(cell)] : 0.0;
    S[0] += p; S[1] += p; S[2] += p;
}

double Solver::minTotalStress(int cell) const {
    double S[6]; totalStress(cell, S);
    double l1, l2, l3; eig3(S, l1, l2, l3);
    return std::min(l1, std::min(l2, l3));
}

bool Solver::setupImpl(const Input& in, int nx, int ny, int nz, const std::vector<double>& dx, const std::vector<double>& dy,
                       const std::vector<double>& dz, double top, const std::vector<double>& p0, const std::vector<double>& phi0,
                       std::string* error, const std::vector<double>* resNodes) {
    auto fail = [&](const std::string& m) { if (error) *error = m; ready_ = false; return false; };
    in_ = in; ready_ = false;
    nx_ = nx; ny_ = ny; nz_ = nz;
    const int n = nx * ny * nz;
    if (n <= 0) return fail("geomechanics: empty grid");
    if (int(p0.size()) != n || int(phi0.size()) != n) return fail("geomechanics: initial pressure / porosity size");
    const int sp = std::max(0, in.sidePad);
    const int px = nx > 1 ? sp : 0, py = ny > 1 ? sp : 0;
    ox_ = px; oy_ = py; oz_ = std::max(0, in.overLayers);
    NX_ = nx + 2 * px; NY_ = ny + 2 * py; NZ_ = nz + oz_ + std::max(0, in.underLayers);
    auto cidx = [&](int i, int j, int k) { return size_t((k * ny + j) * nx + i); };
    // tensor-product element sizes (first row / column / pillar of the reservoir box)
    ex_.assign(size_t(NX_), 0); ey_.assign(size_t(NY_), 0); ez_.assign(size_t(NZ_), 0);
    for (int i = 0; i < NX_; ++i) ex_[size_t(i)] = (i < px || i >= px + nx) ? in.sideThick / std::max(1, sp) : at(dx, cidx(i - px, 0, 0), 1.0);
    for (int j = 0; j < NY_; ++j) ey_[size_t(j)] = (j < py || j >= py + ny) ? in.sideThick / std::max(1, sp) : at(dy, cidx(0, j - py, 0), 1.0);
    for (int k = 0; k < NZ_; ++k) {
        if (k < oz_) ez_[size_t(k)] = in.overThick / std::max(1, in.overLayers);
        else if (k >= oz_ + nz) ez_[size_t(k)] = in.underThick / std::max(1, in.underLayers);
        else ez_[size_t(k)] = at(dz, cidx(0, 0, k - oz_), 1.0);
    }
    // C25 general hexahedra: mechanics node coordinates (reservoir nodes + shifted padding nodes)
    X_.clear();
    if (resNodes) {
        const int nn = (NX_ + 1) * (NY_ + 1) * (NZ_ + 1);
        X_.assign(size_t(nn) * 3, 0.0);
        auto rn = [&](int i, int j, int k) { return 3 * ((size_t(k) * size_t(ny + 1) + size_t(j)) * size_t(nx + 1) + size_t(i)); };
        for (int K = 0; K <= NZ_; ++K)
            for (int J = 0; J <= NY_; ++J)
                for (int I = 0; I <= NX_; ++I) {
                    const int i = std::clamp(I - px, 0, nx), j = std::clamp(J - py, 0, ny), k = std::clamp(K - oz_, 0, nz);
                    double x = (*resNodes)[rn(i, j, k)], y = (*resNodes)[rn(i, j, k) + 1], z = (*resNodes)[rn(i, j, k) + 2];
                    for (int m = I; m < px; ++m) x -= ex_[size_t(m)];
                    for (int m = px + nx; m < I; ++m) x += ex_[size_t(m)];
                    for (int m = J; m < py; ++m) y -= ey_[size_t(m)];
                    for (int m = py + ny; m < J; ++m) y += ey_[size_t(m)];
                    for (int m = K; m < oz_; ++m) z -= ez_[size_t(m)];
                    for (int m = oz_ + nz; m < K; ++m) z += ez_[size_t(m)];
                    const size_t q = 3 * size_t(node(I, J, K));
                    X_[q] = x; X_[q + 1] = y; X_[q + 2] = z;
                }
    }
    const int ne = NX_ * NY_ * NZ_;
    eE_.assign(size_t(ne), in.padE); eNu_.assign(size_t(ne), in.padNu); eA_.assign(size_t(ne), 0.0); eB_.assign(size_t(ne), in.padAlphaT);
    resCell_.assign(size_t(ne), -1);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = cidx(i, j, k);
                const int e = elem(i + px, j + py, k + oz_);
                resCell_[size_t(e)] = int(c);
                eE_[size_t(e)] = at(in.E, c, 5e9); eNu_[size_t(e)] = at(in.nu, c, 0.25);
                eA_[size_t(e)] = at(in.biot, c, 1.0); eB_[size_t(e)] = at(in.alphaT, c, 1e-5);
                if (!(eE_[size_t(e)] > 0) || !(eNu_[size_t(e)] > -1 && eNu_[size_t(e)] < 0.5)) return fail("geomechanics: invalid Young's modulus / Poisson ratio");
            }
    // reservoir cell data: depths, initial stresses, coupling coefficients
    p0c_ = p0; phi0_ = phi0;
    depthC_.assign(size_t(n), 0); sx0_.assign(size_t(n) * 2, 0);
    cp_.assign(size_t(n), 0); cT_.assign(size_t(n), 0); invN_.assign(size_t(n), 0); betaPhi_.assign(size_t(n), 0);
    double zt = top; for (int k = 0; k < nz; ++k) {
        const double h = ez_[size_t(k + oz_)];
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) depthC_[cidx(i, j, k)] = zt + 0.5 * h;
        zt += h;
    }
    if (!X_.empty())   // general hexahedra: centre depth = mean node depth
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            int nd8[8]; elemNodes(i + px, j + py, k + oz_, nd8);
            double z = 0; for (int q = 0; q < 8; ++q) z += X_[3 * size_t(nd8[q]) + 2];
            depthC_[cidx(i, j, k)] = z / 8;
        }
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = cidx(i, j, k);
                const int e = elem(i + px, j + py, k + oz_);
                const double sv = in.sigmaVGrad * depthC_[c];
                sx0_[2 * c] = in.k0 * (sv - p0[c]) + p0[c];   // horizontal total
                sx0_[2 * c + 1] = sv;                          // vertical total
                const double E = eE_[size_t(e)], nu = eNu_[size_t(e)], a = eA_[size_t(e)], b = eB_[size_t(e)];
                const double Kdr = bulkModulus(E, nu);
                invN_[c] = in.Ks > 0 ? std::max(0.0, a - phi0[c]) / in.Ks : 0.0;
                betaPhi_[c] = 3 * b * (a - phi0[c]);
                cp_[c] = a * a / Kdr + invN_[c];
                cT_[c] = 3 * a * b - betaPhi_[c];
            }
    // ---- global stiffness (3 dofs per node), pattern: 27-node stencil
    const int nn = (NX_ + 1) * (NY_ + 1) * (NZ_ + 1);
    K_ = Csr(); K_.n = size_t(nn) * 3;
    K_.ptr.assign(K_.n + 1, 0);
    for (int k = 0; k <= NZ_; ++k)
        for (int j = 0; j <= NY_; ++j)
            for (int i = 0; i <= NX_; ++i) {
                std::vector<int> nb;
                for (int dk = -1; dk <= 1; ++dk) for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
                    const int a = i + di, b = j + dj, cc = k + dk;
                    if (a < 0 || b < 0 || cc < 0 || a > NX_ || b > NY_ || cc > NZ_) continue;
                    nb.push_back(node(a, b, cc));
                }
                std::sort(nb.begin(), nb.end());
                const int nd = node(i, j, k);
                for (int r = 0; r < 3; ++r) {
                    K_.ptr[size_t(3 * nd + r) + 1] = int(nb.size()) * 3;
                }
                (void)nd;
            }
    for (size_t r = 0; r < K_.n; ++r) K_.ptr[r + 1] += K_.ptr[r];
    K_.col.assign(size_t(K_.ptr[K_.n]), 0); K_.val.assign(K_.col.size(), 0.0); K_.diag.assign(K_.n, -1);
    for (int k = 0; k <= NZ_; ++k)
        for (int j = 0; j <= NY_; ++j)
            for (int i = 0; i <= NX_; ++i) {
                std::vector<int> nb;
                for (int dk = -1; dk <= 1; ++dk) for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
                    const int a = i + di, b = j + dj, cc = k + dk;
                    if (a < 0 || b < 0 || cc < 0 || a > NX_ || b > NY_ || cc > NZ_) continue;
                    nb.push_back(node(a, b, cc));
                }
                std::sort(nb.begin(), nb.end());
                const int nd = node(i, j, k);
                for (int r = 0; r < 3; ++r) {
                    int pos = K_.ptr[size_t(3 * nd + r)];
                    for (int m : nb) for (int s = 0; s < 3; ++s) { K_.col[size_t(pos)] = 3 * m + s; if (3 * m + s == 3 * nd + r) K_.diag[size_t(3 * nd + r)] = pos; ++pos; }
                }
            }
    auto addAt = [&](int r, int c, double v) {
        const int* b = &K_.col[size_t(K_.ptr[size_t(r)])];
        const int len = K_.ptr[size_t(r) + 1] - K_.ptr[size_t(r)];
        const int* p = std::lower_bound(b, b + len, c);
        K_.val[size_t(K_.ptr[size_t(r)] + (p - b))] += v;
    };
    const double gp = 1.0 / std::sqrt(3.0);
    for (int k = 0; k < NZ_; ++k)
        for (int j = 0; j < NY_; ++j)
            for (int i = 0; i < NX_; ++i) {
                const int e = elem(i, j, k);
                double D[6][6]; dMatrix(eE_[size_t(e)], eNu_[size_t(e)], D);
                double Ke[24][24]; std::memset(Ke, 0, sizeof(Ke));
                for (int g = 0; g < 8; ++g) {
                    double G[8][3], B[6][24], DB[6][24];
                    const double w = elemGrad(i, j, k, (2 * kLoc[g][0] - 1) * gp, (2 * kLoc[g][1] - 1) * gp, (2 * kLoc[g][2] - 1) * gp, G);
                    if (!(w > 0)) return fail("geomechanics: degenerate (inverted / zero-volume) element - check the corner-point geometry");
                    bMatrix(G, B);
                    for (int r = 0; r < 6; ++r) for (int q = 0; q < 24; ++q) { double s = 0; for (int t = 0; t < 6; ++t) s += D[r][t] * B[t][q]; DB[r][q] = s; }
                    for (int r = 0; r < 24; ++r) for (int q = 0; q < 24; ++q) { double s = 0; for (int t = 0; t < 6; ++t) s += B[t][r] * DB[t][q]; Ke[r][q] += s * w; }
                }
                int nd[8]; elemNodes(i, j, k, nd);
                for (int A = 0; A < 8; ++A) for (int r = 0; r < 3; ++r)
                    for (int Bn = 0; Bn < 8; ++Bn) for (int s = 0; s < 3; ++s) addAt(3 * nd[A] + r, 3 * nd[Bn] + s, Ke[3 * A + r][3 * Bn + s]);
            }
    // ---- boundary conditions per face: roller (normal component), free, fixed (all components)
    for (int f = 0; f < 6; ++f) bc_[f] = in.bcType[f] >= 0 && in.bcType[f] <= 3 ? in.bcType[f] : (f == 4 ? (in.topFree ? 1 : 0) : 0);
    {
        int nForce = 0;
        for (int f = 0; f < 6; ++f) if (bc_[f] == 3 && !in.plateDispControl[f]) ++nForce;
        if (nForce > 1) return fail("geomechanics: at most one force-controlled rigid plate (MECHBC type 3)");
    }
    {   // rigid-body guard: every direction needs a constraint somewhere (x: faces 0/1, y: 2/3, z: 4/5 or any fixed face)
        const bool anyFixed = std::find(bc_, bc_ + 6, 2) != bc_ + 6;
        if (!anyFixed && bc_[0] == 1 && bc_[1] == 1) bc_[0] = 0;
        if (!anyFixed && bc_[2] == 1 && bc_[3] == 1) bc_[2] = 0;
        if (!anyFixed && bc_[4] == 1 && bc_[5] == 1) bc_[5] = 0;
    }
    fixed_.assign(K_.n, 0);
    for (int k = 0; k <= NZ_; ++k)
        for (int j = 0; j <= NY_; ++j)
            for (int i = 0; i <= NX_; ++i) {
                const int nd = node(i, j, k);
                const bool on[6] = {i == 0, i == NX_, j == 0, j == NY_, k == 0, k == NZ_};
                for (int f = 0; f < 6; ++f) {
                    if (!on[f] || bc_[f] == 1) continue;
                    if (bc_[f] == 2) { for (int d = 0; d < 3; ++d) fixed_[size_t(3 * nd + d)] = 1; }
                    else fixed_[size_t(3 * nd + f / 2)] = 1;   // roller / rigid plate: normal component
                }
            }
    for (int f = 0; f < 6; ++f) bcStress_[f] = (bc_[f] == 1 || bc_[f] == 3) ? in.bcStress[f] : 0.0;
    // rigid plates: the normal dofs of the face (a node on two plates belongs to the first)
    for (int f = 0; f < 6; ++f) {
        plateDofs_[size_t(f)].clear(); plateU_[size_t(f)] = 0; plateF_[size_t(f)] = 0;
        if (bc_[f] != 3) continue;
        if (in.plateDispControl[f]) plateU_[size_t(f)] = in.bcDisp[f];
        const int d = f / 2;
        for (int k = 0; k <= NZ_; ++k)
            for (int j = 0; j <= NY_; ++j)
                for (int i = 0; i <= NX_; ++i) {
                    const int idx = d == 0 ? i : d == 1 ? j : k, lim = d == 0 ? NX_ : d == 1 ? NY_ : NZ_;
                    if (idx != (f % 2 == 0 ? 0 : lim)) continue;
                    const int dof = 3 * node(i, j, k) + d;
                    bool taken = false;
                    for (int g = 0; g < f; ++g) taken = taken || std::find(plateDofs_[size_t(g)].begin(), plateDofs_[size_t(g)].end(), dof) != plateDofs_[size_t(g)].end();
                    if (!taken) plateDofs_[size_t(f)].push_back(dof);
                }
    }
    // C16 plasticity: Drucker-Prager constants matched to Mohr-Coulomb, admissible initial state
    {
        const double pi = 3.14159265358979323846;
        const double fr = std::clamp(in.frictionDeg, 0.0, 89.0) * pi / 180, ps = std::clamp(std::min(in.dilationDeg, in.frictionDeg), 0.0, 89.0) * pi / 180;
        auto match = [&](double a, double& et, double& x) {
            if (in.dpMatch == 0) { et = 6 * std::sin(a) / (std::sqrt(3.0) * (3 - std::sin(a))); x = 6 * std::cos(a) / (std::sqrt(3.0) * (3 - std::sin(a))); }
            else if (in.dpMatch == 1) { et = 6 * std::sin(a) / (std::sqrt(3.0) * (3 + std::sin(a))); x = 6 * std::cos(a) / (std::sqrt(3.0) * (3 + std::sin(a))); }
            else { const double t = std::tan(a); et = 3 * t / std::sqrt(9 + 12 * t * t); x = 3 / std::sqrt(9 + 12 * t * t); }
        };
        double dummy = 0;
        match(fr, eta_, xi_); match(ps, etaB_, dummy);
        epsPc_.assign(size_t(ne) * 48, 0.0); epsPt_ = epsPc_;
        ebarC_.assign(size_t(ne) * 8, 0.0); ebarT_ = ebarC_;
        sgp_.assign(size_t(ne) * 48, 0.0);
        cAdj_.assign(size_t(ne), 0.0);
        if (in.plasticity == 1)
            for (int e = 0; e < ne; ++e) {
                const int c = resCell_[size_t(e)];
                if (c < 0) continue;
                const double sh = sx0_[2 * size_t(c)] - p0[size_t(c)], sv = sx0_[2 * size_t(c) + 1] - p0[size_t(c)];
                const double pm = (2 * sh + sv) / 3, J2 = ((sh - pm) * (sh - pm) * 2 + (sv - pm) * (sv - pm)) / 2;
                const double F0 = std::sqrt(J2) - eta_ * pm - xi_ * in.cohesion;
                if (F0 > 0 && xi_ > 0) cAdj_[size_t(e)] = F0 / xi_ * 1.001 + 1.0;
            }
        epc_.assign(size_t(n), 0.0); evp_.assign(size_t(n), 0.0); kmult_.assign(size_t(n), 1.0);
    }
    for (size_t r = 0; r < K_.n; ++r)
        for (int p = K_.ptr[r]; p < K_.ptr[r + 1]; ++p) {
            const size_t c = size_t(K_.col[size_t(p)]);
            if (fixed_[r] || fixed_[c]) K_.val[size_t(p)] = (r == c) ? 1.0 : 0.0;
        }
    u_.assign(K_.n, 0.0); f_.assign(K_.n, 0.0);
    ev_.assign(size_t(n), 0); phi_ = phi0; uxc_.assign(size_t(n), 0); uyc_.assign(size_t(n), 0); uzc_.assign(size_t(n), 0);
    sm_.assign(size_t(n), 0); szz_.assign(size_t(n), 0); sxx_.assign(size_t(n), 0); fail_.assign(size_t(n), 0);
    sfull_.assign(size_t(n) * 6, 0.0); pAbs_ = p0; sm0_.clear();
    ready_ = true;
    // initial stress state report (zero increments)
    std::vector<double> z0(static_cast<size_t>(n), 0.0);
    solve(z0, z0);
    return true;
}

int Solver::pcg(const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt) {
    const size_t N = K_.n;
    const Csr& A = K_;
    auto mv = [&](const std::vector<double>& v, std::vector<double>& y) {
#pragma omp parallel for schedule(static) if (N > 20000)
        for (long long r = 0; r < (long long)N; ++r) {
            double s = 0;
            for (int p = A.ptr[size_t(r)]; p < A.ptr[size_t(r) + 1]; ++p) s += A.val[size_t(p)] * v[size_t(A.col[size_t(p)])];
            y[size_t(r)] = s;
        }
    };
    // symmetric Gauss-Seidel preconditioner z = (D + U)^-1 D (D + L)^-1 r
    auto prec = [&](const std::vector<double>& r, std::vector<double>& z) {
        for (size_t i = 0; i < N; ++i) {
            double s = r[i];
            for (int p = A.ptr[i]; p < A.diag[i]; ++p) s -= A.val[size_t(p)] * z[size_t(A.col[size_t(p)])];
            z[i] = s / A.val[size_t(A.diag[i])];
        }
        for (size_t i = N; i-- > 0;) {
            double s = 0;
            for (int p = A.diag[i] + 1; p < A.ptr[i + 1]; ++p) s += A.val[size_t(p)] * z[size_t(A.col[size_t(p)])];
            z[i] -= s / A.val[size_t(A.diag[i])];
        }
    };
    std::vector<double> r(N), z(N), p(N), q(N);
    mv(x, q);
    double bn = 0; for (size_t i = 0; i < N; ++i) { r[i] = b[i] - q[i]; bn += b[i] * b[i]; }
    bn = std::sqrt(bn);
    if (bn < 1e-300) { std::fill(x.begin(), x.end(), 0.0); return 0; }
    double rn = 0; for (double v : r) rn += v * v;
    if (std::sqrt(rn) <= tol * bn) return 0;
    prec(r, z); p = z;
    double rz = 0; for (size_t i = 0; i < N; ++i) rz += r[i] * z[i];
    int it = 0;
    for (; it < maxIt; ++it) {
        mv(p, q);
        double pq = 0; for (size_t i = 0; i < N; ++i) pq += p[i] * q[i];
        if (!(pq > 0)) break;
        const double al = rz / pq;
        rn = 0;
        for (size_t i = 0; i < N; ++i) { x[i] += al * p[i]; r[i] -= al * q[i]; rn += r[i] * r[i]; }
        if (std::sqrt(rn) <= tol * bn) { ++it; break; }
        prec(r, z);
        double rz2 = 0; for (size_t i = 0; i < N; ++i) rz2 += r[i] * z[i];
        const double be = rz2 / rz; rz = rz2;
        for (size_t i = 0; i < N; ++i) p[i] = z[i] + be * p[i];
    }
    return it;
}

double Solver::faceArea(int f) const {
    double sx = 0, sy = 0, sz = 0;
    for (double v : ex_) sx += v;
    for (double v : ey_) sy += v;
    for (double v : ez_) sz += v;
    const int d = f / 2;
    return d == 0 ? sy * sz : d == 1 ? sx * sz : sx * sy;
}

void Solver::commit() {
    if (in_.plasticity == 1 && epsPt_.size() == epsPc_.size()) { epsPc_ = epsPt_; ebarC_ = ebarT_; }
}

namespace {
// Drucker-Prager return mapping (compression-positive effective stress s = xx yy zz xy yz xz); returns true if plastic.
// dEbar: equivalent plastic strain increment.
bool dpReturn(double s[6], double c, double eta, double xi, double etaB, double K, double G, double H, double& dEbar, double cRes = -1) {
    const double p = (s[0] + s[1] + s[2]) / 3;
    double d[6] = {s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]};
    const double J2 = 0.5 * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) + d[3] * d[3] + d[4] * d[4] + d[5] * d[5];
    const double q = std::sqrt(std::max(J2, 0.0));
    const double F = q - eta * p - xi * c;
    if (F <= 1e-12 * (std::fabs(p) + std::fabs(c) + 1.0)) return false;
    double den = G + K * eta * etaB + xi * xi * H;
    if (den < 1e-3 * G) den = 1e-3 * G;   // strong softening guard
    double dl = F / den;
    if (H < 0 && cRes >= 0 && c + H * xi * dl < cRes) {
        // softening would pass the residual cohesion within this increment: return to the residual cone instead
        // (perfectly plastic at c_res; dl from the residual yield function, exact for the final state)
        c = cRes;
        dl = std::max(0.0, (q - eta * p - xi * c) / (G + K * eta * etaB));
    }
    double pn;
    if (q - G * dl >= 0 || eta <= 0) {     // smooth cone
        const double sc = q > 0 ? 1 - G * dl / q : 0.0;
        for (double& v : d) v *= sc;
        pn = p + K * etaB * dl;
        dEbar = xi * dl;
    } else {                                 // apex (tensile vertex) p = -xi c / eta
        for (double& v : d) v = 0;
        pn = -xi * c / eta;
        dEbar = xi * q / G;
    }
    s[0] = d[0] + pn; s[1] = d[1] + pn; s[2] = d[2] + pn; s[3] = d[3]; s[4] = d[4]; s[5] = d[5];
    return true;
}
} // namespace

void Solver::internalForce(const std::vector<double>& dp, const std::vector<double>& dT, std::vector<double>& fint, bool storeGp) {
    fint.assign(K_.n, 0.0);
    const double gp = 1.0 / std::sqrt(3.0);
    const bool plastic = in_.plasticity == 1;
    int nPl = 0;
    for (int k = 0; k < NZ_; ++k)
        for (int j = 0; j < NY_; ++j)
            for (int i = 0; i < NX_; ++i) {
                const int e = elem(i, j, k), c = resCell_[size_t(e)];
                const double E = eE_[size_t(e)], nu = eNu_[size_t(e)];
                double D[6][6]; dMatrix(E, nu, D);
                const double Kb = bulkModulus(E, nu), G = shearModulus(E, nu);
                int nd[8]; elemNodes(i, j, k, nd);
                double ue[24];
                for (int nn = 0; nn < 8; ++nn) for (int d = 0; d < 3; ++d) ue[3 * nn + d] = u_[size_t(3 * nd[nn] + d)];
                const bool pl = plastic && c >= 0;
                double sLoad = 0, pAbs = 0, s0h = 0, s0v = 0;
                if (c >= 0) {
                    const size_t cc = size_t(c);
                    sLoad = eA_[size_t(e)] * dp[cc] + 3 * Kb * eB_[size_t(e)] * dT[cc];
                    pAbs = p0c_[cc] + dp[cc]; s0h = sx0_[2 * cc]; s0v = sx0_[2 * cc + 1];
                }
                const double coh = std::max(in_.residualCohesion, in_.cohesion + cAdj_[size_t(e)]);
                double fe[24] = {0};
                for (int g = 0; g < 8; ++g) {
                    double Gr[8][3], B[6][24];
                    const double w = elemGrad(i, j, k, (2 * kLoc[g][0] - 1) * gp, (2 * kLoc[g][1] - 1) * gp, (2 * kLoc[g][2] - 1) * gp, Gr);
                    bMatrix(Gr, B);
                    double eps[6] = {0};
                    for (int r = 0; r < 6; ++r) for (int q = 0; q < 24; ++q) eps[r] += B[r][q] * ue[q];
                    const size_t gi = (size_t(e) * 8 + size_t(g)) * 6, gb = size_t(e) * 8 + size_t(g);
                    double epl[6] = {0, 0, 0, 0, 0, 0};
                    if (pl) for (int r = 0; r < 6; ++r) epl[r] = epsPc_[gi + size_t(r)];
                    double ds[6];
                    for (int r = 0; r < 6; ++r) { double v = 0; for (int q = 0; q < 6; ++q) v += D[r][q] * (eps[q] - epl[q]); ds[r] = v; }
                    if (pl) {
                        // trial effective stress (compression +): sigma_0 - p - (D (eps - eps_p) - s m)
                        double se[6] = {s0h - pAbs - (ds[0] - sLoad), s0h - pAbs - (ds[1] - sLoad), s0v - pAbs - (ds[2] - sLoad), -ds[3], -ds[4], -ds[5]};
                        const double seTr[6] = {se[0], se[1], se[2], se[3], se[4], se[5]};
                        const double cNow = std::max(in_.residualCohesion, coh + in_.hardening * ebarC_[gb]);
                        double dEb = 0;
                        if (dpReturn(se, cNow, eta_, xi_, etaB_, Kb, G, in_.hardening, dEb, in_.residualCohesion)) {
                            ++nPl;
                            // d eps_p = C (s_new - s_trial)  (s = const - D (eps - eps_p))
                            const double dsn[6] = {se[0] - seTr[0], se[1] - seTr[1], se[2] - seTr[2], se[3] - seTr[3], se[4] - seTr[4], se[5] - seTr[5]};
                            double dep[6];
                            dep[0] = (dsn[0] - nu * (dsn[1] + dsn[2])) / E; dep[1] = (dsn[1] - nu * (dsn[0] + dsn[2])) / E; dep[2] = (dsn[2] - nu * (dsn[0] + dsn[1])) / E;
                            dep[3] = dsn[3] / G; dep[4] = dsn[4] / G; dep[5] = dsn[5] / G;
                            for (int r = 0; r < 6; ++r) epsPt_[gi + size_t(r)] = epl[r] + dep[r];
                            ebarT_[gb] = ebarC_[gb] + dEb;
                            for (int r = 0; r < 6; ++r) { double v = 0; for (int q = 0; q < 6; ++q) v += D[r][q] * (eps[q] - epsPt_[gi + size_t(q)]); ds[r] = v; }
                        } else {
                            for (int r = 0; r < 6; ++r) epsPt_[gi + size_t(r)] = epl[r];
                            ebarT_[gb] = ebarC_[gb];
                        }
                        if (storeGp) for (int r = 0; r < 6; ++r) sgp_[gi + size_t(r)] = se[r];
                    } else if (storeGp && c >= 0) {
                        const double se[6] = {s0h - pAbs - (ds[0] - sLoad), s0h - pAbs - (ds[1] - sLoad), s0v - pAbs - (ds[2] - sLoad), -ds[3], -ds[4], -ds[5]};
                        for (int r = 0; r < 6; ++r) sgp_[gi + size_t(r)] = se[r];
                    }
                    for (int q = 0; q < 24; ++q) { double v = 0; for (int r = 0; r < 6; ++r) v += B[r][q] * ds[r]; fe[q] += v * w; }
                }
                for (int nn = 0; nn < 8; ++nn) for (int d = 0; d < 3; ++d) fint[size_t(3 * nd[nn] + d)] += fe[3 * nn + d];
            }
    nPlast_ = nPl;
}

bool Solver::equilibrate(const std::vector<double>& dp, const std::vector<double>& dT, const std::vector<double>& fext) {
    for (int f = 0; f < 6; ++f) {
        if (bc_[f] != 3) continue;
        const double sg = (f % 2 == 0) ? 1.0 : -1.0;   // into the body
        for (int dof : plateDofs_[size_t(f)]) u_[size_t(dof)] = sg * plateU_[size_t(f)];
    }
    const size_t N = K_.n;
    std::vector<double> fint, r(N), du(N);
    double fextN = 0;
    for (size_t i = 0; i < N; ++i) if (!fixed_[i]) fextN += fext[i] * fext[i];
    fextN = std::sqrt(fextN);
    const int maxIt = in_.plasticity == 1 ? std::max(1, in_.plasticIter) : 4;
    lastConv_ = false; lastPIt_ = 0; lastCg_ = 0;
    std::vector<double> duOld;
    double omega = 1.0;
    for (int it = 0; it <= maxIt; ++it) {
        internalForce(dp, dT, fint, true);
        double rn = 0, fiN = 0;
        for (size_t i = 0; i < N; ++i) {
            fiN += fint[i] * fint[i];
            r[i] = fixed_[i] ? 0.0 : fext[i] - fint[i];
            rn += r[i] * r[i];
        }
        const double ref = std::max({fextN, std::sqrt(fiN), 1.0});
        if (std::sqrt(rn) <= std::max(in_.plasticTol, 1e-12) * ref) { lastConv_ = true; break; }
        if (it == maxIt) break;
        std::fill(du.begin(), du.end(), 0.0);
        lastCg_ += pcg(r, du, 1e-9, 4000);
        // C16: Aitken dynamic relaxation of the initial-stiffness (modified Newton) iteration - the elastic K is kept
        // (SPD, factor-free PCG) while the relaxation factor tracks the plastic tangent (Kuettler & Wall 2008)
        if (in_.plasticity == 1 && it > 0 && duOld.size() == N) {
            double num = 0, den = 0;
            for (size_t i = 0; i < N; ++i) { const double d = du[i] - duOld[i]; num += duOld[i] * d; den += d * d; }
            if (den > 0) omega = std::clamp(-omega * num / den, 0.05, 8.0);
        }
        duOld = du;
        for (size_t i = 0; i < N; ++i) u_[i] += omega * du[i];
        lastPIt_ = it + 1;
    }
    for (int f = 0; f < 6; ++f) {
        if (bc_[f] != 3) continue;
        const double sg = (f % 2 == 0) ? 1.0 : -1.0;
        double R = 0;
        for (int dof : plateDofs_[size_t(f)]) R += fint[size_t(dof)] - fext[size_t(dof)];
        plateF_[size_t(f)] = sg * R;
    }
    return lastConv_;
}

bool Solver::solve(const std::vector<double>& dp, const std::vector<double>& dT) {
    if (!ready_) return false;
    const int n = nx_ * ny_ * nz_;
    if (int(dp.size()) != n || int(dT.size()) != n) return false;
    std::fill(f_.begin(), f_.end(), 0.0);
    const double gp = 1.0 / std::sqrt(3.0);
    // loads f = sum_e s_e int B^T m dV,  s_e = alpha dp + 3 K_dr beta dT
    for (int k = 0; k < NZ_; ++k)
        for (int j = 0; j < NY_; ++j)
            for (int i = 0; i < NX_; ++i) {
                const int e = elem(i, j, k), c = resCell_[size_t(e)];
                if (c < 0) continue;
                const double s = eA_[size_t(e)] * dp[size_t(c)] + 3 * bulkModulus(eE_[size_t(e)], eNu_[size_t(e)]) * eB_[size_t(e)] * dT[size_t(c)];
                if (s == 0) continue;
                double g[24] = {0};
                for (int q = 0; q < 8; ++q) {
                    double G[8][3];
                    const double w = elemGrad(i, j, k, (2 * kLoc[q][0] - 1) * gp, (2 * kLoc[q][1] - 1) * gp, (2 * kLoc[q][2] - 1) * gp, G);
                    for (int nn = 0; nn < 8; ++nn) for (int d = 0; d < 3; ++d) g[3 * nn + d] += G[nn][d] * w;
                }
                int nd[8]; elemNodes(i, j, k, nd);
                for (int nn = 0; nn < 8; ++nn) for (int d = 0; d < 3; ++d) f_[size_t(3 * nd[nn] + d)] += s * g[3 * nn + d];
            }
    // boundary loads: normal stress increment s (compression +) on a free face = traction -s n, consistent nodal loads
    // of the bilinear face elements (area / 4 per node); rigid plates are handled by their net force below
    for (int f = 0; f < 6; ++f) {
        if (bcStress_[f] == 0 || bc_[f] != 1) continue;
        const int d = f / 2;                       // normal direction
        const double sgn = (f % 2 == 0) ? 1.0 : -1.0;   // outward normal -1 (min face) / +1 (max face): force = -s n
        const int fixedIdx = (f % 2 == 0) ? 0 : (d == 0 ? NX_ : d == 1 ? NY_ : NZ_);
        const int na = d == 0 ? NY_ : NX_, nb2 = d == 2 ? NY_ : NZ_;
        for (int b = 0; b < nb2; ++b)
            for (int a = 0; a < na; ++a) {
                int ii[4], jj[4], kk[4]; double area;
                const int ca[4] = {a, a + 1, a + 1, a}, cb[4] = {b, b, b + 1, b + 1};
                if (d == 0) { area = ey_[size_t(a)] * ez_[size_t(b)]; for (int q = 0; q < 4; ++q) { ii[q] = fixedIdx; jj[q] = ca[q]; kk[q] = cb[q]; } }
                else if (d == 1) { area = ex_[size_t(a)] * ez_[size_t(b)]; for (int q = 0; q < 4; ++q) { ii[q] = ca[q]; jj[q] = fixedIdx; kk[q] = cb[q]; } }
                else { area = ex_[size_t(a)] * ey_[size_t(b)]; for (int q = 0; q < 4; ++q) { ii[q] = ca[q]; jj[q] = cb[q]; kk[q] = fixedIdx; } }
                for (int q = 0; q < 4; ++q) f_[size_t(3 * node(ii[q], jj[q], kk[q]) + d)] += sgn * bcStress_[f] * area / 4;
            }
    }
    bool plates = false;
    for (int f = 0; f < 6; ++f) plates = plates || bc_[f] == 3;
    if (in_.plasticity != 1 && !plates) {      // linear elastic, ordinary supports: one PCG solve
        std::vector<double> rhs = f_;
        for (size_t r = 0; r < K_.n; ++r) if (fixed_[r]) rhs[r] = 0.0;
        lastCg_ = pcg(rhs, u_, 1e-8, 4000);
        lastPIt_ = 1; lastConv_ = true; nPlast_ = 0;
        cellResults(dp, dT, false);
        return true;
    }
    int fc = -1;
    for (int f = 0; f < 6; ++f) {
        if (bc_[f] != 3) continue;
        if (in_.plateDispControl[f]) plateU_[size_t(f)] = in_.bcDisp[f];
        else if (fc < 0) fc = f;
    }
    bool ok = true;
    if (fc < 0) ok = equilibrate(dp, dT, f_);
    else {
        // force-controlled plate: secant iteration on its displacement (exact in two solves when elastic)
        const double A = faceArea(fc), target = bcStress_[fc] * A;
        double Em = 0; int ne = 0;
        for (size_t e = 0; e < eE_.size(); ++e) if (resCell_[e] >= 0) { Em += eE_[e]; ++ne; }
        Em = ne ? Em / ne : in_.padE;
        const int d = fc / 2;
        double L = 0; for (double v : (d == 0 ? ex_ : d == 1 ? ey_ : ez_)) L += v;
        const double kEst = Em * A / std::max(L, 1e-9), tolF = 1e-9 * Em * A;
        double U0 = plateU_[size_t(fc)];
        ok = equilibrate(dp, dT, f_);
        double F0 = plateF_[size_t(fc)];
        if (std::fabs(F0 - target) > tolF) {
            double U1 = U0 + (target - F0) / kEst;
            for (int it = 0; it < 15; ++it) {
                plateU_[size_t(fc)] = U1;
                ok = equilibrate(dp, dT, f_);
                const double F1 = plateF_[size_t(fc)];
                if (std::fabs(F1 - target) <= tolF) break;
                const double slope = (F1 - F0) / (U1 - U0);
                const double U2 = std::fabs(slope) > 1e-6 * kEst ? U1 + (target - F1) / slope : U1 + (target - F1) / kEst;
                U0 = U1; F0 = F1; U1 = U2;
            }
        }
    }
    cellResults(dp, dT, true);
    return ok || in_.plasticity == 1;   // a non-converged plastic iteration keeps the best state (reported by lastConverged())
}

void Solver::cellResults(const std::vector<double>& dp, const std::vector<double>& dT, bool fromGp) {
    const double fr = in_.frictionDeg * 3.14159265358979323846 / 180;
    for (int k = 0; k < NZ_; ++k)
        for (int j = 0; j < NY_; ++j)
            for (int i = 0; i < NX_; ++i) {
                const int e = elem(i, j, k), c = resCell_[size_t(e)];
                if (c < 0) continue;
                double G[8][3], B[6][24], D[6][6];
                elemGrad(i, j, k, 0, 0, 0, G); bMatrix(G, B); dMatrix(eE_[size_t(e)], eNu_[size_t(e)], D);
                int nd[8]; elemNodes(i, j, k, nd);
                double ue[24];
                for (int nn = 0; nn < 8; ++nn) for (int d = 0; d < 3; ++d) ue[3 * nn + d] = u_[size_t(3 * nd[nn] + d)];
                double eps[6] = {0};
                for (int r = 0; r < 6; ++r) for (int q = 0; q < 24; ++q) eps[r] += B[r][q] * ue[q];
                const size_t cc = size_t(c);
                const double ev = eps[0] + eps[1] + eps[2];
                ev_[cc] = ev;
                phi_[cc] = std::clamp(phi0_[cc] + eA_[size_t(e)] * ev + invN_[cc] * dp[cc] - betaPhi_[cc] * dT[cc], 1e-4, 0.999);
                double ux = 0, uy = 0, uz = 0;
                for (int nn = 0; nn < 8; ++nn) { ux += ue[3 * nn]; uy += ue[3 * nn + 1]; uz += ue[3 * nn + 2]; }
                uxc_[cc] = ux / 8; uyc_[cc] = uy / 8; uzc_[cc] = uz / 8;
                double se[6];
                if (fromGp) {
                    for (int r = 0; r < 6; ++r) { double v = 0; for (int g = 0; g < 8; ++g) v += sgp_[(size_t(e) * 8 + size_t(g)) * 6 + size_t(r)]; se[r] = v / 8; }
                    double eb = 0, evp = 0;
                    for (int g = 0; g < 8; ++g) {
                        const size_t gi = (size_t(e) * 8 + size_t(g)) * 6;
                        eb += ebarT_[size_t(e) * 8 + size_t(g)];
                        evp += epsPt_[gi] + epsPt_[gi + 1] + epsPt_[gi + 2];
                    }
                    epc_[cc] = eb / 8; evp_[cc] = evp / 8;
                    kmult_[cc] = in_.permPlastic > 0 ? std::min(std::max(1.0, in_.permPlasticMax), std::exp(in_.permPlastic * epc_[cc])) : 1.0;
                } else {
                    const double s = eA_[size_t(e)] * dp[cc] + 3 * bulkModulus(eE_[size_t(e)], eNu_[size_t(e)]) * eB_[size_t(e)] * dT[cc];
                    double ds[6];
                    for (int r = 0; r < 6; ++r) { double v = 0; for (int q = 0; q < 6; ++q) v += D[r][q] * eps[q]; ds[r] = v - (r < 3 ? s : 0.0); }
                    // effective stress (compression positive): sigma0 - d sigma - p
                    const double p = p0c_[cc] + dp[cc];
                    const double t[6] = {sx0_[2 * cc] - ds[0] - p, sx0_[2 * cc] - ds[1] - p, sx0_[2 * cc + 1] - ds[2] - p, -ds[3], -ds[4], -ds[5]};
                    for (int r = 0; r < 6; ++r) se[r] = t[r];
                    epc_[cc] = 0; evp_[cc] = 0; kmult_[cc] = 1.0;
                }
                sxx_[cc] = se[0]; szz_[cc] = se[2]; sm_[cc] = (se[0] + se[1] + se[2]) / 3;
                for (int r = 0; r < 6; ++r) sfull_[6 * cc + size_t(r)] = se[r];
                pAbs_[cc] = p0c_[cc] + dp[cc];
                if (in_.permStress > 0) {   // C25 stress-dependent permeability (mean effective stress change)
                    if (sm0_.size() != sm_.size()) sm0_.assign(sm_.size(), std::numeric_limits<double>::quiet_NaN());
                    if (std::isnan(sm0_[cc])) sm0_[cc] = sm_[cc];
                    kmult_[cc] = std::clamp(kmult_[cc] * std::exp(-in_.permStress * (sm_[cc] - sm0_[cc])), 1e-4, std::max(1.0, in_.permPlasticMax));
                }
                double s1, s2, s3; eig3(se, s1, s2, s3);
                const double tau = 0.5 * (s1 - s3), smn = 0.5 * (s1 + s3);
                const double cAd = in_.cohesion + (cAdj_.empty() ? 0.0 : cAdj_[size_t(e)]);
                const double strength = cAd * std::cos(fr) + smn * std::sin(fr);
                double fi = strength > 0 ? tau / strength : 10.0;
                if (in_.tensile > 0 && s3 < 0) fi = std::max(fi, -s3 / in_.tensile);
                fail_[cc] = std::min(fi, 10.0);
            }
    uplift_ = 0; upliftRes_ = 0;
    for (int j = 0; j <= NY_; ++j)
        for (int i = 0; i <= NX_; ++i) {
            uplift_ = std::max(uplift_, -u_[size_t(3 * node(i, j, 0) + 2)]);
            upliftRes_ = std::max(upliftRes_, -u_[size_t(3 * node(i, j, oz_) + 2)]);
        }
}

} // namespace geomech
} // namespace ressim
