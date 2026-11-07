// Cubic EOS thermodynamics: see eos.h for the models and references.
#include "ressim/eos.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ressim {
namespace eos {

namespace {
constexpr double SQ2 = 1.4142135623730951;

// dense LU with partial pivoting (small systems: nc <= ~30)
bool luSolve(std::vector<double> A, std::vector<double>& b, int n) {
    std::vector<int> piv(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) {
        int p = k; double mx = std::fabs(A[size_t(k) * n + k]);
        for (int i = k + 1; i < n; ++i) { const double v = std::fabs(A[size_t(i) * n + k]); if (v > mx) { mx = v; p = i; } }
        if (!(mx > 1e-300)) return false;
        if (p != k) { for (int j = 0; j < n; ++j) std::swap(A[size_t(k) * n + j], A[size_t(p) * n + j]); std::swap(b[size_t(k)], b[size_t(p)]); }
        const double d = A[size_t(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double f = A[size_t(i) * n + k] / d;
            if (f == 0) continue;
            for (int j = k + 1; j < n; ++j) A[size_t(i) * n + j] -= f * A[size_t(k) * n + j];
            b[size_t(i)] -= f * b[size_t(k)];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[size_t(i)];
        for (int j = i + 1; j < n; ++j) s -= A[size_t(i) * n + j] * b[size_t(j)];
        b[size_t(i)] = s / A[size_t(i) * n + i];
    }
    return true;
}

// LU factorisation in place (partial pivoting) + solve for several right-hand sides
bool luFactor(std::vector<double>& A, std::vector<int>& piv, int n) {
    piv.resize(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) {
        int p = k; double mx = std::fabs(A[size_t(k) * n + k]);
        for (int i = k + 1; i < n; ++i) { const double v = std::fabs(A[size_t(i) * n + k]); if (v > mx) { mx = v; p = i; } }
        if (!(mx > 1e-300)) return false;
        piv[size_t(k)] = p;
        if (p != k) for (int j = 0; j < n; ++j) std::swap(A[size_t(k) * n + j], A[size_t(p) * n + j]);
        const double d = A[size_t(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            double& f = A[size_t(i) * n + k];
            f /= d;
            if (f == 0) continue;
            for (int j = k + 1; j < n; ++j) A[size_t(i) * n + j] -= f * A[size_t(k) * n + j];
        }
    }
    return true;
}
void luSolveF(const std::vector<double>& A, const std::vector<int>& piv, std::vector<double>& b, int n) {
    // full-row interchanges were applied during the factorisation (LAPACK getrf style): permute b first
    for (int k = 0; k < n; ++k) if (piv[size_t(k)] != k) std::swap(b[size_t(k)], b[size_t(piv[size_t(k)])]);
    for (int k = 0; k < n; ++k)
        for (int i = k + 1; i < n; ++i) b[size_t(i)] -= A[size_t(i) * n + k] * b[size_t(k)];
    for (int i = n - 1; i >= 0; --i) {
        double s = b[size_t(i)];
        for (int j = i + 1; j < n; ++j) s -= A[size_t(i) * n + j] * b[size_t(j)];
        b[size_t(i)] = s / A[size_t(i) * n + i];
    }
}

double omegaA(const Model& m, const Component& c) { return c.omegaA > 0 ? c.omegaA : m.type == EosType::PR ? 0.457235528921 : 0.427480233540; }
double omegaB(const Model& m, const Component& c) { return c.omegaB > 0 ? c.omegaB : m.type == EosType::PR ? 0.0777960739039 : 0.0866403499650; }
double mOf(const Model& m, double w) {
    if (m.type == EosType::SRK) return 0.480 + 1.574 * w - 0.176 * w * w;
    if (m.pr78 && w > 0.49) return 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w * w * w;
    return 0.37464 + 1.54226 * w - 0.26992 * w * w;
}
double vcOf(const Component& c) {
    if (c.vc > 0) return c.vc;
    const double zc = c.zc > 0 ? c.zc : 0.2905 - 0.085 * c.omega;
    return zc * Rgas * c.tc / c.pc;
}
}  // namespace

double Model::delta1() const { return type == EosType::PR ? 1.0 + SQ2 : 1.0; }
double Model::delta2() const { return type == EosType::PR ? 1.0 - SQ2 : 0.0; }

TParams atTemperature(const Model& m, double T) {
    const int nc = m.nc();
    TParams t; t.T = T;
    t.a.resize(static_cast<size_t>(nc)); t.dadT.resize(static_cast<size_t>(nc)); t.d2adT2.resize(static_cast<size_t>(nc)); t.b.resize(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) {
        const Component& c = m.comp[size_t(i)];
        const double mm = mOf(m, c.omega), s = std::sqrt(T / c.tc), u = 1.0 + mm * (1.0 - s);
        const double ac = omegaA(m, c) * Rgas * Rgas * c.tc * c.tc / c.pc;
        const double alpha = u * u, dal = -mm * u * s / T, d2al = mm * (mm * s * s + u * s) / (2.0 * T * T);
        t.a[size_t(i)] = ac * alpha; t.dadT[size_t(i)] = ac * dal; t.d2adT2[size_t(i)] = ac * d2al;
        t.b[size_t(i)] = omegaB(m, c) * Rgas * c.tc / c.pc;
    }
    t.aij.resize(size_t(nc) * nc); t.daijdT.resize(size_t(nc) * nc); t.d2aijdT2.resize(size_t(nc) * nc);
    for (int i = 0; i < nc; ++i)
        for (int j = 0; j < nc; ++j) {
            const double ai = t.a[size_t(i)], aj = t.a[size_t(j)], d1 = t.dadT[size_t(i)], d2 = t.dadT[size_t(j)];
            const double g = std::sqrt(ai * aj), kk = 1.0 - m.k(i, j);
            const double dg = (d1 * aj + ai * d2) / (2.0 * g);
            const double d2g = (t.d2adT2[size_t(i)] * aj + 2.0 * d1 * d2 + ai * t.d2adT2[size_t(j)]) / (2.0 * g) -
                               (d1 * aj + ai * d2) * (d1 * aj + ai * d2) / (4.0 * g * g * g);
            t.aij[size_t(i) * nc + j] = kk * g; t.daijdT[size_t(i) * nc + j] = kk * dg; t.d2aijdT2[size_t(i) * nc + j] = kk * d2g;
        }
    return t;
}

namespace {
// real roots of Z^3 + c2 Z^2 + c1 Z + c0, ascending, Newton-polished
std::vector<double> solveCubic(double c2, double c1, double c0) {
    std::vector<double> r;
    const double p = c1 - c2 * c2 / 3.0, q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    const double disc = 0.25 * q * q + p * p * p / 27.0;
    if (disc > 0) {
        const double sd = std::sqrt(disc);
        const double t = std::cbrt(-0.5 * q + sd) + std::cbrt(-0.5 * q - sd);
        r.push_back(t - c2 / 3.0);
    } else if (p < 0) {
        const double rr = 2.0 * std::sqrt(-p / 3.0);
        double arg = 3.0 * q / (p * rr);
        arg = std::max(-1.0, std::min(1.0, arg));
        const double ph = std::acos(arg) / 3.0;
        for (int k = 0; k < 3; ++k) r.push_back(rr * std::cos(ph - 2.0 * M_PI * k / 3.0) - c2 / 3.0);
    } else r.push_back(-c2 / 3.0);
    for (double& z : r)
        for (int it = 0; it < 3; ++it) {
            const double f = ((z + c2) * z + c1) * z + c0, d = (3.0 * z + 2.0 * c2) * z + c1;
            if (std::fabs(d) < 1e-300) break;
            z -= f / d;
        }
    std::sort(r.begin(), r.end());
    return r;
}

struct Mix { double nT = 0, B = 0, D = 0, DT = 0, DTT = 0; std::vector<double> x, Di, DiT; };
Mix mixture(const TParams& tp, const std::vector<double>& n, bool wantT) {
    const int nc = int(n.size());
    Mix M; M.x.resize(static_cast<size_t>(nc)); M.Di.assign(size_t(nc), 0.0); if (wantT) M.DiT.assign(size_t(nc), 0.0);
    for (double v : n) M.nT += v;
    for (int i = 0; i < nc; ++i) {
        M.x[size_t(i)] = n[size_t(i)] / M.nT;
        M.B += n[size_t(i)] * tp.b[size_t(i)];
        double s = 0, sT = 0, sTT = 0;
        for (int j = 0; j < nc; ++j) {
            s += n[size_t(j)] * tp.aij[size_t(i) * nc + j];
            if (wantT) { sT += n[size_t(j)] * tp.daijdT[size_t(i) * nc + j]; sTT += n[size_t(j)] * tp.d2aijdT2[size_t(i) * nc + j]; }
        }
        M.Di[size_t(i)] = 2.0 * s; M.D += n[size_t(i)] * s;
        if (wantT) { M.DiT[size_t(i)] = 2.0 * sT; M.DT += n[size_t(i)] * sT; M.DTT += n[size_t(i)] * sTT; }
    }
    return M;
}
}  // namespace

std::vector<double> cubicRoots(const Model& m, const TParams& tp, double P, const std::vector<double>& x) {
    const Mix M = mixture(tp, x, false);
    const double RT = Rgas * tp.T, A = M.D / (M.nT * M.nT) * P / (RT * RT), B = M.B / M.nT * P / RT;
    const double d1 = m.delta1(), d2 = m.delta2();
    return solveCubic((d1 + d2 - 1.0) * B - 1.0, A + d1 * d2 * B * B - (d1 + d2) * B * (B + 1.0), -(A * B + d1 * d2 * B * B * (B + 1.0)));
}

PhaseEval evalPhase(const Model& m, const TParams& tp, double P, const std::vector<double>& n, const EvalOptions& o) {
    const int nc = m.nc();
    const double T = tp.T, RT = Rgas * T;
    const bool wantT = o.dT || o.hres;
    const Mix M = mixture(tp, n, wantT);
    const double nT = M.nT;
    const double d1 = m.delta1(), d2 = m.delta2();
    const double Am = M.D / (nT * nT) * P / (RT * RT), Bm = M.B / nT * P / RT;
    std::vector<double> roots = solveCubic((d1 + d2 - 1.0) * Bm - 1.0, Am + d1 * d2 * Bm * Bm - (d1 + d2) * Bm * (Bm + 1.0),
                                           -(Am * Bm + d1 * d2 * Bm * Bm * (Bm + 1.0)));
    std::vector<double> ok;
    for (double z : roots) if (z > Bm * (1.0 + 1e-12)) ok.push_back(z);
    if (ok.empty()) ok.push_back(std::max(roots.back(), Bm * (1.0 + 1e-9)));
    double Z = ok.front();
    if (ok.size() > 1) {
        if (o.root == Root::Vapour) Z = ok.back();
        else if (o.root == Root::MinGibbs) {
            auto gres = [&](double z) {
                return z - 1.0 - std::log(z - Bm) - Am / (Bm * (d1 - d2)) * std::log((z + d1 * Bm) / (z + d2 * Bm));
            };
            if (gres(ok.back()) < gres(ok.front())) Z = ok.back();
        }
    }
    PhaseEval e;
    e.Z = Z;
    // reduced residual Helmholtz function F(T, V, n) and its derivatives (Michelsen & Mollerup)
    const double V = Z * nT * RT / P, B = M.B, D = M.D;
    const double VB = V - B, V1 = V + d1 * B, V2 = V + d2 * B;
    const double g = std::log(1.0 - B / V), gV = B / (V * VB), gB = -1.0 / VB;
    const double gVV = 1.0 / (V * V) - 1.0 / (VB * VB), gBV = 1.0 / (VB * VB), gBB = -1.0 / (VB * VB);
    const double f = std::log(V1 / V2) / (Rgas * B * (d1 - d2));
    const double fV = -1.0 / (Rgas * V1 * V2);
    const double fB = -(f + V * fV) / B;
    const double fVV = (2.0 * V + (d1 + d2) * B) / (Rgas * V1 * V1 * V2 * V2);
    const double fBV = -(2.0 * fV + V * fVV) / B;
    const double fBB = -(2.0 * fB + V * fBV) / B;
    const double FV = -nT * gV - D * fV / T, FD = -f / T;
    const double FnV = -gV, FnB = -gB, FVV = -nT * gVV - D * fVV / T, FBV = -nT * gBV - D * fBV / T, FBB = -nT * gBB - D * fBB / T;
    const double FDV = -fV / T, FBD = -fB / T;
    const double FB = -nT * gB - D * fB / T, Fn = -g;
    (void)FV;
    const double PV = -RT * FVV - nT * RT / (V * V);
    std::vector<double> Pi(static_cast<size_t>(nc)), FiV(static_cast<size_t>(nc));
    e.lnphi.resize(static_cast<size_t>(nc));
    const double lnZ = std::log(Z);
    double cmix = 0, mw = 0;
    for (int i = 0; i < nc; ++i) {
        const double bi = tp.b[size_t(i)], Di = M.Di[size_t(i)];
        const double Fi = Fn + FB * bi + FD * Di;
        FiV[size_t(i)] = FnV + FBV * bi + FDV * Di;
        Pi[size_t(i)] = -RT * FiV[size_t(i)] + RT / V;
        const double ci = m.comp[size_t(i)].sshift * bi;
        e.lnphi[size_t(i)] = Fi - lnZ - ci * P / RT;
        cmix += M.x[size_t(i)] * ci; mw += M.x[size_t(i)] * m.comp[size_t(i)].mw;
    }
    e.vm = Z * RT / P - cmix;
    e.mw = mw;
    e.rho = mw / e.vm;
    e.dvm_dP = 1.0 / (nT * PV);
    if (o.dn) {
        e.dlnphi_dn.resize(size_t(nc) * nc);
        e.dvm_dx.resize(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) {
            const double bi = tp.b[size_t(i)], Di = M.Di[size_t(i)];
            for (int j = 0; j < nc; ++j) {
                const double bj = tp.b[size_t(j)], Dj = M.Di[size_t(j)];
                const double Fij = FnB * (bi + bj) + FBB * bi * bj + FBD * (bi * Dj + bj * Di) + FD * 2.0 * tp.aij[size_t(i) * nc + j];
                e.dlnphi_dn[size_t(i) * nc + j] = nT * Fij + 1.0 + nT * Pi[size_t(i)] * Pi[size_t(j)] / (RT * PV);
            }
            const double Vi = -Pi[size_t(i)] / PV, ci = m.comp[size_t(i)].sshift * bi;
            e.dvm_dx[size_t(i)] = (Vi - ci) - e.vm;
        }
    }
    if (o.dP) {
        e.dlnphi_dP.resize(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) {
            const double Vi = -Pi[size_t(i)] / PV, ci = m.comp[size_t(i)].sshift * tp.b[size_t(i)];
            e.dlnphi_dP[size_t(i)] = Vi / RT - 1.0 / P - ci / RT;
        }
    }
    if (wantT) {
        // total temperature derivatives at constant V, n (a_ij depend on T)
        const double FVT = D * fV / (T * T) - fV * M.DT / T;
        const double PT = P / T - RT * FVT;
        e.dvm_dT = -PT / (nT * PV);
        e.dlnphi_dT.resize(static_cast<size_t>(nc));
        double hs = 0;
        for (int i = 0; i < nc; ++i) {
            const double bi = tp.b[size_t(i)], Di = M.Di[size_t(i)];
            const double FiT = D * fB / (T * T) * bi + f / (T * T) * Di - fB / T * bi * M.DT + FD * M.DiT[size_t(i)];
            const double Vi = -Pi[size_t(i)] / PV, ci = m.comp[size_t(i)].sshift * bi;
            e.dlnphi_dT[size_t(i)] = FiT + 1.0 / T - Vi * PT / RT + ci * P / (RT * T);
            hs += M.x[size_t(i)] * e.dlnphi_dT[size_t(i)];
        }
        e.hres = -RT * T * hs;
        if (!o.dT) e.dlnphi_dT.clear();
    }
    return e;
}

std::vector<double> wilsonK(const Model& m, double T, double P) {
    std::vector<double> K(size_t(m.nc()));
    for (int i = 0; i < m.nc(); ++i) {
        const Component& c = m.comp[size_t(i)];
        K[size_t(i)] = c.pc / P * std::exp(5.373 * (1.0 + c.omega) * (1.0 - c.tc / T));
    }
    return K;
}

// ---------------------------------------------------------------------------------------------- Rachford-Rice
double rachfordRice(const std::vector<double>& z, const std::vector<double>& K, bool clampToPhysical, int* iterations) {
    const int nc = int(z.size());
    double kmax = -1e300, kmin = 1e300;
    for (int i = 0; i < nc; ++i) if (z[size_t(i)] > 0) { kmax = std::max(kmax, K[size_t(i)]); kmin = std::min(kmin, K[size_t(i)]); }
    auto F = [&](double b, double* dF) {
        double s = 0, d = 0;
        for (int i = 0; i < nc; ++i) {
            const double t = K[size_t(i)] - 1.0, den = 1.0 + b * t;
            s += z[size_t(i)] * t / den; d -= z[size_t(i)] * t * t / (den * den);
        }
        if (dF) *dF = d;
        return s;
    };
    if (kmax <= 1.0) return clampToPhysical ? 0.0 : -1.0;   // all K <= 1: liquid
    if (kmin >= 1.0) return clampToPhysical ? 1.0 : 2.0;    // all K >= 1: vapour
    // window where all phase compositions stay positive (Whitson & Michelsen)
    double lo = 1.0 / (1.0 - kmax), hi = 1.0 / (1.0 - kmin);
    for (int i = 0; i < nc; ++i) {
        if (z[size_t(i)] <= 0) continue;
        const double Ki = K[size_t(i)];
        if (Ki > 1.0) lo = std::max(lo, (Ki * z[size_t(i)] - 1.0) / (Ki - 1.0));
        else if (Ki < 1.0) hi = std::min(hi, (1.0 - z[size_t(i)]) / (1.0 - Ki));
    }
    if (clampToPhysical) {
        if (F(0.0, nullptr) <= 0) return 0.0;
        if (F(1.0, nullptr) >= 0) return 1.0;
        lo = std::max(lo, 0.0); hi = std::min(hi, 1.0);
    }
    double b = 0.5 * (lo + hi);
    int it = 0;
    for (; it < 100; ++it) {
        double d; const double v = F(b, &d);
        if (v > 0) lo = b; else hi = b;
        double nb = b - v / d;
        if (!(nb > lo && nb < hi)) nb = 0.5 * (lo + hi);
        if (std::fabs(nb - b) < 1e-15 * std::max(1.0, std::fabs(b))) { b = nb; break; }
        b = nb;
        if (hi - lo < 1e-16) break;
    }
    if (iterations) *iterations = it;
    return b;
}

bool isVapourLike(const Model& m, const std::vector<double>& z, double T) {
    double num = 0, den = 0;
    for (int i = 0; i < m.nc(); ++i) { const double v = vcOf(m.comp[size_t(i)]); num += z[size_t(i)] * v * m.comp[size_t(i)].tc; den += z[size_t(i)] * v; }
    return den > 0 && T > num / den;
}

// ---------------------------------------------------------------------------------------------- stability
StabilityResult stability(const Model& m, const TParams& tp, double P, const std::vector<double>& z) {
    const int nc = m.nc();
    StabilityResult R;
    EvalOptions eo; eo.dn = false; eo.dP = false;
    const PhaseEval fz = evalPhase(m, tp, P, z, eo);
    std::vector<double> d(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) d[size_t(i)] = (z[size_t(i)] > 0 ? std::log(z[size_t(i)]) : -1e3) + fz.lnphi[size_t(i)];
    const std::vector<double> Kw = wilsonK(m, tp.T, P);
    double best = 0; std::vector<double> bestK;
    for (int trial = 0; trial < 2; ++trial) {
        std::vector<double> W(static_cast<size_t>(nc)), lnW(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) { W[size_t(i)] = trial == 0 ? z[size_t(i)] * Kw[size_t(i)] : z[size_t(i)] / Kw[size_t(i)]; }
        for (int i = 0; i < nc; ++i) lnW[size_t(i)] = std::log(std::max(W[size_t(i)], 1e-300));
        double tm = 0; bool trivial = false, conv = false;
        auto tmOf = [&](const std::vector<double>& w, const std::vector<double>& lphi) {
            double s = 1.0;
            for (int i = 0; i < nc; ++i) if (w[size_t(i)] > 0) s += w[size_t(i)] * (std::log(w[size_t(i)]) + lphi[size_t(i)] - d[size_t(i)] - 1.0);
            return s;
        };
        // successive substitution
        int it = 0;
        PhaseEval fw;
        for (; it < 40; ++it) {
            fw = evalPhase(m, tp, P, W, eo);
            double err = 0, triv = 0;
            for (int i = 0; i < nc; ++i) {
                const double nl = d[size_t(i)] - fw.lnphi[size_t(i)];
                err = std::max(err, std::fabs(nl - lnW[size_t(i)]));
                lnW[size_t(i)] = nl; W[size_t(i)] = std::exp(nl);
            }
            double sw = 0; for (double w : W) sw += w;
            for (int i = 0; i < nc; ++i) { const double r = W[size_t(i)] / sw - z[size_t(i)]; triv += r * r; }
            if (triv < 1e-12) { trivial = true; break; }
            if (err < 1e-10) { conv = true; break; }
            if (it >= 8 && err < 1e-3) break;   // second-order from here
        }
        R.iterations += it;
        if (!trivial && !conv) {
            // second-order (Michelsen): alpha_i = 2 sqrt(W_i), Levenberg-damped Newton on tm*(alpha)
            EvalOptions en; en.dn = true; en.dP = false;
            std::vector<double> alpha(static_cast<size_t>(nc));
            for (int i = 0; i < nc; ++i) alpha[size_t(i)] = 2.0 * std::sqrt(W[size_t(i)]);
            fw = evalPhase(m, tp, P, W, en);
            double tmCur = tmOf(W, fw.lnphi);
            double mu = 0;
            for (int k = 0; k < 40; ++k) {
                double sw = 0; for (double w : W) sw += w;
                std::vector<double> gvec(static_cast<size_t>(nc)), H(size_t(nc) * nc);
                double gn = 0;
                for (int i = 0; i < nc; ++i) {
                    const double si = std::sqrt(W[size_t(i)]);
                    gvec[size_t(i)] = si * (std::log(std::max(W[size_t(i)], 1e-300)) + fw.lnphi[size_t(i)] - d[size_t(i)]);
                    gn = std::max(gn, std::fabs(gvec[size_t(i)]));
                }
                if (gn < 1e-11) { conv = true; break; }
                for (int i = 0; i < nc; ++i)
                    for (int j = 0; j < nc; ++j)
                        H[size_t(i) * nc + j] = (i == j ? 1.0 : 0.0) +
                                                std::sqrt(W[size_t(i)] * W[size_t(j)]) * fw.dlnphi_dn[size_t(i) * nc + j] / sw;
                bool accepted = false;
                for (int t = 0; t < 12 && !accepted; ++t) {
                    std::vector<double> A = H, rhs(static_cast<size_t>(nc));
                    for (int i = 0; i < nc; ++i) { A[size_t(i) * nc + i] += mu; rhs[size_t(i)] = -gvec[size_t(i)]; }
                    if (!luSolve(A, rhs, nc)) { mu = std::max(1e-3, mu * 10); continue; }
                    std::vector<double> a2(static_cast<size_t>(nc)), W2(static_cast<size_t>(nc));
                    for (int i = 0; i < nc; ++i) { a2[size_t(i)] = std::max(1e-150, alpha[size_t(i)] + rhs[size_t(i)]); W2[size_t(i)] = 0.25 * a2[size_t(i)] * a2[size_t(i)]; }
                    const PhaseEval f2 = evalPhase(m, tp, P, W2, en);
                    const double tm2 = tmOf(W2, f2.lnphi);
                    if (tm2 <= tmCur + 1e-14 * std::fabs(tmCur)) {
                        alpha = a2; W = W2; fw = f2; tmCur = tm2; accepted = true; mu = mu * 0.1;
                        if (mu < 1e-10) mu = 0;
                    } else mu = std::max(1e-3, mu * 10);
                }
                ++R.iterations;
                if (!accepted) break;
                double triv = 0; double s2 = 0; for (double w : W) s2 += w;
                for (int i = 0; i < nc; ++i) { const double r = W[size_t(i)] / s2 - z[size_t(i)]; triv += r * r; }
                if (triv < 1e-12) { trivial = true; break; }
            }
        }
        if (trivial) continue;
        fw = evalPhase(m, tp, P, W, eo);
        tm = tmOf(W, fw.lnphi);
        if (tm < best) {
            best = tm;
            double sw = 0; for (double w : W) sw += w;
            bestK.assign(size_t(nc), 1.0);
            for (int i = 0; i < nc; ++i) {
                const double wi = W[size_t(i)] / sw;
                if (z[size_t(i)] > 0) bestK[size_t(i)] = trial == 0 ? wi / z[size_t(i)] : z[size_t(i)] / std::max(wi, 1e-300);
            }
        }
    }
    R.tpd = best;
    R.stable = !(best < -1e-9);
    if (!R.stable) R.K = bestK;
    return R;
}

// ---------------------------------------------------------------------------------------------- flash
FlashResult flash(const Model& m, const TParams& tp, double P, const std::vector<double>& zIn, const FlashOptions& o) {
    const int nc = m.nc();
    double zt = 0; for (double q : zIn) zt += q;
    std::vector<double> z = zIn;   // feed mole fractions (the flash works on a unit feed)
    if (zt > 0 && std::fabs(zt - 1.0) > 1e-15) for (double& q : z) q /= zt;
    FlashResult r;
    EvalOptions eo; eo.dn = false; eo.dP = false;
    EvalOptions full; full.dn = true; full.dP = true;
    auto single = [&](bool converged) {
        r.phases = 1; r.converged = converged;
        r.vapourLabel = isVapourLike(m, z, tp.T);
        const PhaseEval e = evalPhase(m, tp, P, z, full);
        r.x = z; r.y = z; r.K.assign(size_t(nc), 1.0);
        if (r.vapourLabel) { r.V = e; r.beta = 1.0; } else { r.L = e; r.beta = 0.0; }
        return r;
    };
    std::vector<double> K;
    bool tested = false;
    if (o.Kinit && int(o.Kinit->size()) == nc) {
        K = *o.Kinit;
        const double b = rachfordRice(z, K, false);
        if (!(b > 0 && b < 1) && !o.skipStability) K.clear();
    }
    if (K.empty()) {
        const StabilityResult s = stability(m, tp, P, z);
        tested = true;
        r.iterations += s.iterations;
        if (s.stable) return single(true);
        K = s.K;
    }
    std::vector<double> lnK(static_cast<size_t>(nc)), x(static_cast<size_t>(nc)), y(static_cast<size_t>(nc));
    double beta = 0.5, err = 1e300;
    PhaseEval L, V;
    auto comps = [&](double b) {
        double sx = 0, sy = 0;
        for (int i = 0; i < nc; ++i) {
            x[size_t(i)] = z[size_t(i)] / (1.0 + b * (K[size_t(i)] - 1.0)); y[size_t(i)] = K[size_t(i)] * x[size_t(i)];
            sx += x[size_t(i)]; sy += y[size_t(i)];
        }
        for (int i = 0; i < nc; ++i) { x[size_t(i)] /= sx; y[size_t(i)] /= sy; }
    };
    // Newton on the vapour mole numbers v (feed = 1 mol): g_i = ln f_iV - ln f_iL, full Jacobian (quadratic convergence).
    // On convergence L / V hold the phase evaluations (with derivatives) at the solution.
    std::vector<double> v(static_cast<size_t>(nc)), l(static_cast<size_t>(nc)), g(static_cast<size_t>(nc)), J(size_t(nc) * nc), dv(static_cast<size_t>(nc));
    std::vector<int> piv;
    auto newton = [&](int maxIt) -> bool {
        for (int k = 0; k < maxIt; ++k) {
            double Vt = 0, Lt = 0;
            for (int i = 0; i < nc; ++i) { Vt += v[size_t(i)]; Lt += l[size_t(i)]; }
            L = evalPhase(m, tp, P, l, full); V = evalPhase(m, tp, P, v, full);
            double gmax = 0;
            for (int i = 0; i < nc; ++i) {
                if (z[size_t(i)] <= 0) { g[size_t(i)] = 0; continue; }
                g[size_t(i)] = std::log(v[size_t(i)] / Vt) + V.lnphi[size_t(i)] - std::log(l[size_t(i)] / Lt) - L.lnphi[size_t(i)];
                gmax = std::max(gmax, std::fabs(g[size_t(i)]));
            }
            r.iterations++;
            if (gmax < o.tol) return true;
            for (int i = 0; i < nc; ++i)
                for (int j = 0; j < nc; ++j) {
                    // absent components (z = 0) are not unknowns: identity row and zero column
                    if (z[size_t(i)] <= 0 || z[size_t(j)] <= 0) { J[size_t(i) * nc + j] = i == j ? 1.0 : 0.0; continue; }
                    double a = -1.0 / Vt + V.dlnphi_dn[size_t(i) * nc + j] / Vt - 1.0 / Lt + L.dlnphi_dn[size_t(i) * nc + j] / Lt;
                    if (i == j) a += 1.0 / v[size_t(i)] + 1.0 / l[size_t(i)];
                    J[size_t(i) * nc + j] = a;
                }
            for (int i = 0; i < nc; ++i) dv[size_t(i)] = -g[size_t(i)];
            if (!luFactor(J, piv, nc)) return false;
            luSolveF(J, piv, dv, nc);
            // keep 0 < v_i < z_i
            double s = 1.0;
            for (int i = 0; i < nc; ++i) {
                if (z[size_t(i)] <= 0) continue;
                const double nv = v[size_t(i)] + dv[size_t(i)];
                if (nv <= 0) s = std::min(s, 0.9 * v[size_t(i)] / -dv[size_t(i)]);
                if (nv >= z[size_t(i)]) s = std::min(s, 0.9 * l[size_t(i)] / dv[size_t(i)]);
            }
            for (int i = 0; i < nc; ++i) if (z[size_t(i)] > 0) { v[size_t(i)] += s * dv[size_t(i)]; l[size_t(i)] = z[size_t(i)] - v[size_t(i)]; }
            double Vn = 0; for (double q : v) Vn += q;
            if (Vn < 1e-12 || Vn > 1.0 - 1e-12) return false;
        }
        return false;
    };
    auto twoPhase = [&](bool conv) -> FlashResult& {
        double Vt = 0; for (double q : v) Vt += q;
        r.phases = 2; r.converged = conv; r.beta = Vt;
        r.x.resize(static_cast<size_t>(nc)); r.y.resize(static_cast<size_t>(nc)); r.K.resize(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) {
            r.y[size_t(i)] = v[size_t(i)] / Vt; r.x[size_t(i)] = l[size_t(i)] / (1.0 - Vt);
            r.K[size_t(i)] = r.x[size_t(i)] > 0 ? r.y[size_t(i)] / r.x[size_t(i)] : 1.0;
        }
        if (conv) { r.L = std::move(L); r.V = std::move(V); }   // evaluated at the converged (scaled) phase moles
        else { r.L = evalPhase(m, tp, P, r.x, full); r.V = evalPhase(m, tp, P, r.y, full); }
        // absent components: equilibrium K = phi_L / phi_V (trace partition; warm start when the component arrives)
        for (int i = 0; i < nc; ++i)
            if (z[size_t(i)] <= 0 && r.L.lnphi.size() == size_t(nc) && r.V.lnphi.size() == size_t(nc))
                r.K[size_t(i)] = std::exp(r.L.lnphi[size_t(i)] - r.V.lnphi[size_t(i)]);
        // the lighter phase is the vapour
        if (r.L.rho < r.V.rho) { std::swap(r.L, r.V); std::swap(r.x, r.y); r.beta = 1.0 - r.beta; for (double& q : r.K) q = 1.0 / q; }
        return r;
    };
    // fast path (previous converged K, e.g. the last Newton iteration of the reservoir simulator): straight to Newton
    if (!tested && !K.empty()) {
        const double b0 = rachfordRice(z, K, false);
        if (b0 > 0 && b0 < 1) {
            comps(b0);
            for (int i = 0; i < nc; ++i) { v[size_t(i)] = b0 * y[size_t(i)]; l[size_t(i)] = z[size_t(i)] - v[size_t(i)]; }
            if (newton(8)) {
                double Vt = 0; for (double q : v) Vt += q;
                double triv = 0;   // reject the trivial solution (identical phases)
                if (Vt > 0 && Vt < 1)
                    for (int i = 0; i < nc; ++i)
                        if (z[size_t(i)] > 0) { const double lk = std::log((v[size_t(i)] / Vt) / (l[size_t(i)] / (1.0 - Vt))); triv += lk * lk; }
                if (Vt > 0 && Vt < 1 && triv > 1e-8) return twoPhase(true);
            }
        }
    }
    for (int i = 0; i < nc; ++i) lnK[size_t(i)] = std::log(std::max(K[size_t(i)], 1e-300));
    std::vector<double> d1v, d2v;   // GDEM history
    int it = 0;
    for (; it < o.maxSS; ++it) {
        for (int i = 0; i < nc; ++i) K[size_t(i)] = std::exp(lnK[size_t(i)]);
        beta = rachfordRice(z, K, false);
        comps(std::max(0.0, std::min(1.0, beta)) == beta ? beta : std::max(-0.5, std::min(1.5, beta)));
        L = evalPhase(m, tp, P, x, eo); V = evalPhase(m, tp, P, y, eo);
        std::vector<double> nl(static_cast<size_t>(nc));
        err = 0; double triv = 0;
        for (int i = 0; i < nc; ++i) {
            nl[size_t(i)] = L.lnphi[size_t(i)] - V.lnphi[size_t(i)];
            err = std::max(err, std::fabs(nl[size_t(i)] - lnK[size_t(i)]));
            triv += nl[size_t(i)] * nl[size_t(i)];
        }
        if (triv < 1e-8) {   // trivial solution: single phase
            r.iterations += it;
            if (!tested) { const StabilityResult s = stability(m, tp, P, z); if (s.stable) return single(true); }
            return single(false);
        }
        std::vector<double> dl(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) dl[size_t(i)] = nl[size_t(i)] - lnK[size_t(i)];
        lnK = nl;
        // GDEM (dominant eigenvalue) acceleration every 5th iteration
        if (!d1v.empty() && it % 5 == 4) {
            double a = 0, b = 0;
            for (int i = 0; i < nc; ++i) { a += dl[size_t(i)] * dl[size_t(i)]; b += d1v[size_t(i)] * dl[size_t(i)]; }
            const double lam = b != 0 ? a / b : 0;
            if (lam > 0 && lam < 0.97) for (int i = 0; i < nc; ++i) lnK[size_t(i)] += dl[size_t(i)] * lam / (1.0 - lam);
        }
        d1v = dl;
        if (err < 1e-4 && beta > 0 && beta < 1) break;
        if (err < o.tol) break;
    }
    r.iterations += it;
    for (int i = 0; i < nc; ++i) K[size_t(i)] = std::exp(lnK[size_t(i)]);
    beta = rachfordRice(z, K, false);
    if (!(beta > 0 && beta < 1)) {
        // negative-flash solution: the feed is single phase (or the K estimate is poor) -> decide by stability
        if (!tested) {
            const StabilityResult s = stability(m, tp, P, z);
            if (s.stable) return single(true);
            K = s.K; beta = rachfordRice(z, K, true);
        } else beta = std::max(1e-6, std::min(1.0 - 1e-6, beta));
        if (!(beta > 0 && beta < 1)) return single(false);
    }
    comps(beta);
    for (int i = 0; i < nc; ++i) { v[size_t(i)] = beta * y[size_t(i)]; l[size_t(i)] = z[size_t(i)] - v[size_t(i)]; }
    const bool conv = newton(o.maxNewton);
    double Vt = 0; for (double q : v) Vt += q;
    if (!(Vt > 0 && Vt < 1)) return single(false);
    return twoPhase(conv);
}


// ---------------------------------------------------------------------------------------------- flash sensitivities
bool flashSensitivities(const Model& m, const TParams& tp, double P, const std::vector<double>& z, const FlashResult& f, FlashSens& S, bool withT) {
    const int nc = m.nc();
    if (f.phases != 2) return false;
    const double V = f.beta, L = 1.0 - V;
    double zt = 0; for (double q : z) zt += q;
    std::vector<double> v(static_cast<size_t>(nc)), l(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) { v[size_t(i)] = V * f.y[size_t(i)] * zt; l[size_t(i)] = L * f.x[size_t(i)] * zt; }
    const double Vt = V * zt, Lt = L * zt;
    // the flash result already holds both phases with composition / pressure derivatives (scale invariant): reuse them
    // unless temperature derivatives are wanted
    const bool reuse = !withT && f.L.dlnphi_dn.size() == size_t(nc) * nc && f.V.dlnphi_dn.size() == size_t(nc) * nc &&
                       f.L.dlnphi_dP.size() == size_t(nc) && f.V.dlnphi_dP.size() == size_t(nc);
    PhaseEval ELo, EVo;
    if (!reuse) { EvalOptions o; o.dn = true; o.dP = true; o.dT = withT; ELo = evalPhase(m, tp, P, l, o); EVo = evalPhase(m, tp, P, v, o); }
    const PhaseEval& EL = reuse ? f.L : ELo;
    const PhaseEval& EV = reuse ? f.V : EVo;
    // J = d g / d v (g_i = ln f_iV(v) - ln f_iL(z - v)); A = d ln f_L / d l
    std::vector<double> J(size_t(nc) * nc), AL(size_t(nc) * nc);
    // absent components (z_i = 0): trace partition dv_i = c_i dz_i, c_i = K_i V / (L + K_i V), K_i = phi_iL / phi_iV
    std::vector<double> cTrace(static_cast<size_t>(nc), 0.0);
    std::vector<char> absent(static_cast<size_t>(nc), 0);
    for (int i = 0; i < nc; ++i)
        if (!(l[size_t(i)] > 0 && v[size_t(i)] > 0)) {
            absent[size_t(i)] = 1;
            const double Ki = std::exp(EL.lnphi[size_t(i)] - EV.lnphi[size_t(i)]);
            cTrace[size_t(i)] = Ki * Vt / (Lt + Ki * Vt);
        }
    for (int i = 0; i < nc; ++i)
        for (int j = 0; j < nc; ++j) {
            if (absent[size_t(i)]) { AL[size_t(i) * nc + j] = 0.0; J[size_t(i) * nc + j] = i == j ? 1.0 : 0.0; continue; }
            const double a = (i == j ? 1.0 / l[size_t(i)] : 0.0) - 1.0 / Lt + EL.dlnphi_dn[size_t(i) * nc + j] / Lt;
            const double b = (i == j ? 1.0 / v[size_t(i)] : 0.0) - 1.0 / Vt + EV.dlnphi_dn[size_t(i) * nc + j] / Vt;
            AL[size_t(i) * nc + j] = a; J[size_t(i) * nc + j] = a + b;
        }
    std::vector<int> piv;
    const bool fac = luFactor(J, piv, nc);
    auto solve = [&](std::vector<double> rhs) { if (!fac) rhs.assign(size_t(nc), 0.0); else luSolveF(J, piv, rhs, nc); return rhs; };
    auto finish = [&](const std::vector<double>& dv, const std::vector<double>& dzv, double& dbeta, double* dx, double* dy, size_t stride) {
        double dV = 0, dZ = 0; for (int i = 0; i < nc; ++i) { dV += dv[size_t(i)]; dZ += dzv[size_t(i)]; }
        const double dL = dZ - dV;
        dbeta = (dV * zt - Vt * dZ) / (zt * zt);
        for (int i = 0; i < nc; ++i) {
            dy[size_t(i) * stride] = (dv[size_t(i)] - f.y[size_t(i)] * dV) / Vt;
            dx[size_t(i) * stride] = (dzv[size_t(i)] - dv[size_t(i)] - f.x[size_t(i)] * dL) / Lt;
        }
    };
    const std::vector<double> zero(static_cast<size_t>(nc), 0.0);
    std::vector<double> r(static_cast<size_t>(nc));
    S.dx_dP.assign(size_t(nc), 0); S.dy_dP.assign(size_t(nc), 0); S.dx_dT.assign(size_t(nc), 0); S.dy_dT.assign(size_t(nc), 0);
    for (int i = 0; i < nc; ++i) r[size_t(i)] = absent[size_t(i)] ? 0.0 : -(EV.dlnphi_dP[size_t(i)] - EL.dlnphi_dP[size_t(i)]);
    finish(solve(r), zero, S.dbeta_dP, S.dx_dP.data(), S.dy_dP.data(), 1);
    if (withT && EV.dlnphi_dT.size() == size_t(nc) && EL.dlnphi_dT.size() == size_t(nc)) {
        for (int i = 0; i < nc; ++i) r[size_t(i)] = absent[size_t(i)] ? 0.0 : -(EV.dlnphi_dT[size_t(i)] - EL.dlnphi_dT[size_t(i)]);
        finish(solve(r), zero, S.dbeta_dT, S.dx_dT.data(), S.dy_dT.data(), 1);
    }
    S.dbeta_dz.assign(size_t(nc), 0); S.dx_dz.assign(size_t(nc) * nc, 0); S.dy_dz.assign(size_t(nc) * nc, 0);
    for (int j = 0; j < nc; ++j) {
        // d g_i / d z_j = - d ln f_iL / d l_j
        for (int i = 0; i < nc; ++i) r[size_t(i)] = absent[size_t(i)] ? (i == j ? cTrace[size_t(i)] : 0.0) : AL[size_t(i) * nc + j];
        std::vector<double> dzv(static_cast<size_t>(nc), 0.0); dzv[size_t(j)] = 1.0;
        finish(solve(r), dzv, S.dbeta_dz[size_t(j)], &S.dx_dz[size_t(j)], &S.dy_dz[size_t(j)], size_t(nc));
    }
    return true;
}

// ---------------------------------------------------------------------------------------------- saturation pressure
double saturationPressure(const Model& m, double T, const std::vector<double>& z, bool bubble, double pGuess, double pLo, double pHi,
                          std::vector<double>* incipient) {
    const int nc = m.nc();
    const TParams tp = atTemperature(m, T);
    double P = pGuess;
    std::vector<double> K;
    if (!(P > 0)) {   // Wilson estimate
        if (bubble) { P = 0; for (int i = 0; i < nc; ++i) { const auto& c = m.comp[size_t(i)]; P += z[size_t(i)] * c.pc * std::exp(5.373 * (1 + c.omega) * (1 - c.tc / T)); } }
        else { double s = 0; for (int i = 0; i < nc; ++i) { const auto& c = m.comp[size_t(i)]; s += z[size_t(i)] / (c.pc * std::exp(5.373 * (1 + c.omega) * (1 - c.tc / T))); } P = 1.0 / s; }
        P = std::max(pLo, std::min(pHi, P));
    }
    K = wilsonK(m, T, P);
    EvalOptions eo; eo.dn = false; eo.dP = true;
    std::vector<double> w(static_cast<size_t>(nc));
    for (int outer = 0; outer < 200; ++outer) {
        // incipient phase composition from the current K
        double s = 0;
        for (int i = 0; i < nc; ++i) { w[size_t(i)] = bubble ? z[size_t(i)] * K[size_t(i)] : z[size_t(i)] / K[size_t(i)]; s += w[size_t(i)]; }
        for (int i = 0; i < nc; ++i) w[size_t(i)] /= s;
        const PhaseEval fz = evalPhase(m, tp, P, z, eo), fw = evalPhase(m, tp, P, w, eo);
        // K = phiL / phiV; residual sum(z K) - 1 (bubble) or sum(z / K) - 1 (dew), Newton in ln P
        double S = 0, dS = 0, triv = 0;
        for (int i = 0; i < nc; ++i) {
            const double lnK = bubble ? fz.lnphi[size_t(i)] - fw.lnphi[size_t(i)] : fw.lnphi[size_t(i)] - fz.lnphi[size_t(i)];
            const double dlnK = bubble ? fz.dlnphi_dP[size_t(i)] - fw.dlnphi_dP[size_t(i)] : fw.dlnphi_dP[size_t(i)] - fz.dlnphi_dP[size_t(i)];
            K[size_t(i)] = std::exp(lnK);
            triv += lnK * lnK;
            if (bubble) { S += z[size_t(i)] * K[size_t(i)]; dS += z[size_t(i)] * K[size_t(i)] * dlnK; }
            else { S += z[size_t(i)] / K[size_t(i)]; dS -= z[size_t(i)] / K[size_t(i)] * dlnK; }
        }
        if (triv < 1e-6) break;   // collapsed to the trivial solution (near-critical / poor start): bracket instead
        const double res = std::log(S);
        if (std::fabs(res) < 1e-11 && outer > 2) { if (incipient) *incipient = w; return P; }
        // d ln S / d ln P = P dS / S (composition of the incipient phase lagged)
        const double der = P * dS / S;
        double dlnP = std::fabs(der) > 1e-12 ? -res / der : 0.0;
        dlnP = std::max(-0.5, std::min(0.5, dlnP));
        P = std::max(pLo, std::min(pHi, P * std::exp(dlnP)));
    }
    // Robust fallback: bracket the upper boundary of the two-phase region with the stability test (scan down from pHi),
    // then bisect. Converges for near-critical mixtures where the successive-substitution / Newton iteration collapses.
    double hi = pHi, lo = -1;
    StabilityResult sLo;
    if (!stability(m, tp, hi, z).stable) return -1.0;   // unstable even at pHi
    for (double p = pHi * 0.92; p >= pLo; p *= 0.92) {
        StabilityResult s = stability(m, tp, p, z);
        if (!s.stable) { lo = p; sLo = s; break; }
        hi = p;
    }
    if (lo < 0) return -1.0;   // single phase over the whole range
    for (int it = 0; it < 80 && hi - lo > 1e-10 * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        StabilityResult s = stability(m, tp, mid, z);
        if (s.stable) hi = mid; else { lo = mid; sLo = s; }
    }
    if (incipient && !sLo.K.empty()) {
        incipient->assign(size_t(nc), 0.0);
        double sum = 0;
        for (int i = 0; i < nc; ++i) { (*incipient)[size_t(i)] = bubble ? z[size_t(i)] * sLo.K[size_t(i)] : z[size_t(i)] / sLo.K[size_t(i)]; sum += (*incipient)[size_t(i)]; }
        for (double& q : *incipient) q /= sum;
    }
    return 0.5 * (lo + hi);
}

// ---------------------------------------------------------------------------------------------- LBC viscosity
double lbcViscosity(const Model& m, const std::vector<double>& x, double T, double molarDensity) {
    const int nc = m.nc();
    const double atm = 101325.0;
    double num = 0, den = 0, tcm = 0, pcm = 0, mwm = 0, vcm = 0;
    for (int i = 0; i < nc; ++i) {
        const Component& c = m.comp[size_t(i)];
        const double mwg = c.mw * 1e3, pca = c.pc / atm, tr = T / c.tc;
        const double xi = std::pow(c.tc, 1.0 / 6.0) / (std::sqrt(mwg) * std::pow(pca, 2.0 / 3.0));
        const double mui = (tr <= 1.5 ? 34e-5 * std::pow(tr, 0.94) : 17.78e-5 * std::pow(4.58 * tr - 1.67, 0.625)) / xi;   // cP
        num += x[size_t(i)] * mui * std::sqrt(mwg); den += x[size_t(i)] * std::sqrt(mwg);
        tcm += x[size_t(i)] * c.tc; pcm += x[size_t(i)] * pca; mwm += x[size_t(i)] * mwg;
        vcm += x[size_t(i)] * (c.vcVisc > 0 ? c.vcVisc : vcOf(c));
    }
    const double mu0 = num / den;
    const double xim = std::pow(tcm, 1.0 / 6.0) / (std::sqrt(mwm) * std::pow(pcm, 2.0 / 3.0));
    const double rr = molarDensity * vcm;
    const double poly = m.lbc[0] + rr * (m.lbc[1] + rr * (m.lbc[2] + rr * (m.lbc[3] + rr * m.lbc[4])));
    const double mu = mu0 + (std::pow(poly, 4) - 1e-4) / xim;   // cP
    return mu * 1e-3;
}


double lbcViscosity(const Model& m, const std::vector<double>& x, double T, double molarDensity, double* dmu_dx, double* dmu_dxi) {
    const int nc = m.nc();
    const double atm = 101325.0;
    std::vector<double> mui(static_cast<size_t>(nc)), sq(static_cast<size_t>(nc)), vci(static_cast<size_t>(nc));
    double num = 0, den = 0, tcm = 0, pcm = 0, mwm = 0, vcm = 0;
    for (int i = 0; i < nc; ++i) {
        const Component& c = m.comp[size_t(i)];
        const double mwg = c.mw * 1e3, pca = c.pc / atm, tr = T / c.tc;
        const double xi = std::pow(c.tc, 1.0 / 6.0) / (std::sqrt(mwg) * std::pow(pca, 2.0 / 3.0));
        mui[size_t(i)] = (tr <= 1.5 ? 34e-5 * std::pow(tr, 0.94) : 17.78e-5 * std::pow(4.58 * tr - 1.67, 0.625)) / xi;
        sq[size_t(i)] = std::sqrt(mwg);
        vci[size_t(i)] = c.vcVisc > 0 ? c.vcVisc : vcOf(c);
        num += x[size_t(i)] * mui[size_t(i)] * sq[size_t(i)]; den += x[size_t(i)] * sq[size_t(i)];
        tcm += x[size_t(i)] * c.tc; pcm += x[size_t(i)] * pca; mwm += x[size_t(i)] * mwg; vcm += x[size_t(i)] * vci[size_t(i)];
    }
    const double mu0 = num / den;
    const double xim = std::pow(tcm, 1.0 / 6.0) / (std::sqrt(mwm) * std::pow(pcm, 2.0 / 3.0));
    const double rr = molarDensity * vcm;
    const double* a = m.lbc;
    const double poly = a[0] + rr * (a[1] + rr * (a[2] + rr * (a[3] + rr * a[4])));
    const double dpoly = a[1] + rr * (2 * a[2] + rr * (3 * a[3] + rr * 4 * a[4]));
    const double term = (poly * poly * poly * poly - 1e-4) / xim;
    const double p3 = 4 * poly * poly * poly * dpoly / xim;
    if (dmu_dxi) *dmu_dxi = p3 * vcm * 1e-3;
    if (dmu_dx)
        for (int i = 0; i < nc; ++i) {
            const Component& c = m.comp[size_t(i)];
            const double dmu0 = (mui[size_t(i)] - mu0) * sq[size_t(i)] / den;
            const double dlnxim = c.tc / (6 * tcm) - c.mw * 1e3 / (2 * mwm) - 2 * (c.pc / atm) / (3 * pcm);
            dmu_dx[i] = (dmu0 + p3 * molarDensity * vci[size_t(i)] - term * dlnxim) * 1e-3;
        }
    return (mu0 + term) * 1e-3;
}

}  // namespace eos
}  // namespace ressim
