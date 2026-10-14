// Block CPR solver: see blocksolve.h.
#include "ressim/blocksolve.h"
#include <cmath>
#include <chrono>
#include <memory>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ressim {

static BlockAccelerator* g_accel = nullptr;
void setBlockAccelerator(BlockAccelerator* a) { g_accel = a; }
BlockAccelerator* blockAccelerator() { return g_accel; }

namespace {

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

constexpr int kPar = 6000;

// ---------------------------------------------------------------- small dense block helpers
template <int NB> inline void bmul(const double* A, const double* B, double* C) {         // C = A B
    for (int i = 0; i < NB; ++i)
        for (int j = 0; j < NB; ++j) {
            double s = 0;
            for (int k = 0; k < NB; ++k) s += A[i * NB + k] * B[k * NB + j];
            C[i * NB + j] = s;
        }
}
template <int NB> inline void bmulsub(const double* A, const double* B, double* C) {      // C -= A B
    for (int i = 0; i < NB; ++i)
        for (int j = 0; j < NB; ++j) {
            double s = 0;
            for (int k = 0; k < NB; ++k) s += A[i * NB + k] * B[k * NB + j];
            C[i * NB + j] -= s;
        }
}
template <int NB> inline void bmv(const double* A, const double* x, double* y) {          // y = A x
    for (int i = 0; i < NB; ++i) { double s = 0; for (int k = 0; k < NB; ++k) s += A[i * NB + k] * x[k]; y[i] = s; }
}
template <int NB> inline void bmvsub(const double* A, const double* x, double* y) {       // y -= A x
    for (int i = 0; i < NB; ++i) { double s = 0; for (int k = 0; k < NB; ++k) s += A[i * NB + k] * x[k]; y[i] -= s; }
}
// inverse by Gauss-Jordan with partial pivoting; false when singular
template <int NB> bool binv(const double* A, double* Ai) {
    double M[NB][2 * NB];
    for (int i = 0; i < NB; ++i)
        for (int j = 0; j < NB; ++j) { M[i][j] = A[i * NB + j]; M[i][NB + j] = i == j ? 1.0 : 0.0; }
    double scale = 0;
    for (int i = 0; i < NB * NB; ++i) scale = std::max(scale, std::fabs(A[i]));
    if (!(scale > 0) || !std::isfinite(scale)) return false;
    for (int k = 0; k < NB; ++k) {
        int pr = k;
        for (int i = k + 1; i < NB; ++i) if (std::fabs(M[i][k]) > std::fabs(M[pr][k])) pr = i;
        if (std::fabs(M[pr][k]) < 1e-14 * scale) return false;
        if (pr != k) for (int j = 0; j < 2 * NB; ++j) std::swap(M[k][j], M[pr][j]);
        const double d = 1.0 / M[k][k];
        for (int j = 0; j < 2 * NB; ++j) M[k][j] *= d;
        for (int i = 0; i < NB; ++i) {
            if (i == k || M[i][k] == 0) continue;
            const double f = M[i][k];
            for (int j = 0; j < 2 * NB; ++j) M[i][j] -= f * M[k][j];
        }
    }
    for (int i = 0; i < NB; ++i) for (int j = 0; j < NB; ++j) Ai[i * NB + j] = M[i][NB + j];
    return true;
}

// ---------------------------------------------------------------- vector kernels
static double dotv(const double* a, const double* b, int n) {
    double s = 0;
    if (useParallel(n > kPar)) {
        #pragma omp parallel for reduction(+:s) schedule(static)
        for (int i = 0; i < n; ++i) s += a[i] * b[i];
    } else for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}
template <class F> static void forN(int n, F f) {
    if (useParallel(n > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) f(i);
    } else for (int i = 0; i < n; ++i) f(i);
}

struct CoreBase {
    virtual ~CoreBase() = default;
    virtual void setup(const CSR& A, int na, int nw, const double* acc, const BlockSolverOptions& o) = 0;
    virtual void apply(const double* r, double* z) = 0;
    virtual int amgLevels() const = 0;
};

template <int NB> struct Core : CoreBase {
    static constexpr int BS = NB * NB;
    int na = 0, nw = 0, nc = 0, N = 0;
    // block-CSR of B' (cells)
    std::vector<int> rp, ci, dpos;
    std::vector<double> val;            // B' blocks (row = equation, column = unknown)
    std::vector<double> lu, invD;       // block ILU(0) factors
    // wells
    std::vector<int> cwPtr, cwW; std::vector<double> cwV;   // C: per cell, (well, nb values)
    std::vector<int> dwPtr, dwA; std::vector<double> dwV;   // D: per well, (cell, nb values)
    std::vector<double> ew;                                  // E diagonal
    // CPR
    bool useCpr = true;
    std::vector<double> wgt;            // na * NB decoupling weights
    CSR Ap;
    AMG amg;
    int reuseLeft = 0;
    std::vector<int> lastRowPtr, lastCol;
    std::vector<double> rp_, xp_, t_, rc_;
    std::vector<int> mark;
    int pIts = 8; double pTol = 1e-2; LinResult pres_;
    BlockAccelerator* accel = nullptr;

    int amgLevels() const override { return useCpr ? amg.levels() : 0; }

    void setup(const CSR& A, int na_, int nw_, const double* acc, const BlockSolverOptions& o) override {
        na = na_; nw = nw_; nc = na * NB; N = A.n;
        useCpr = o.cpr; pIts = o.pressureIts; pTol = o.pressureTol;
    amg.sweeps = 2; amg.overCorrect = true;
        buildBlocks(A);
        buildWells(A);
        foldWellDiagonal();
        accel = nullptr;
        if (o.accelerator && g_accel && g_accel->setup(na, NB, rp, ci, val)) {
            accel = g_accel;
            useCpr = false;
            return finishSetup();
        }
        factorILU();
        if (useCpr) {
            weights(acc);
            pressureMatrix(o.wellsInPressure);
            const bool same = Ap.rowPtr == lastRowPtr && Ap.col == lastCol;
            if (reuseLeft > 0 && same) --reuseLeft;
            else { amg.setup(Ap); reuseLeft = o.amgReuse; lastRowPtr = Ap.rowPtr; lastCol = Ap.col; }
        }
        finishSetup();
    }

    void finishSetup() {
        rp_.assign(size_t(na), 0.0);
        xp_.assign(size_t(na), 0.0);
        t_.assign(size_t(nc), 0.0);
        rc_.assign(size_t(nc), 0.0);
    }

    void buildBlocks(const CSR& A) {
        rp.assign(size_t(na) + 1, 0);
        mark.assign(size_t(na), -1);
        std::vector<std::vector<int>> rows; rows.resize(size_t(na));
        for (int a = 0; a < na; ++a) {
            auto& r = rows[size_t(a)];
            for (int e = 0; e < NB; ++e) {
                const int row = a * NB + e;
                for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                    const int c = A.col[p];
                    if (c >= nc) continue;
                    const int bc = c / NB;
                    if (mark[size_t(bc)] != a) { mark[size_t(bc)] = a; r.push_back(bc); }
                }
            }
            if (mark[size_t(a)] != a) { mark[size_t(a)] = a; r.push_back(a); }
            std::sort(r.begin(), r.end());
            rp[size_t(a) + 1] = rp[size_t(a)] + int(r.size());
        }
        ci.resize(size_t(rp[size_t(na)]));
        dpos.assign(size_t(na), -1);
        for (int a = 0; a < na; ++a) {
            int q = rp[size_t(a)];
            for (int bc : rows[size_t(a)]) { if (bc == a) dpos[size_t(a)] = q; ci[size_t(q++)] = bc; }
        }
        val.assign(ci.size() * BS, 0.0);
        std::fill(mark.begin(), mark.end(), -1);
        for (int a = 0; a < na; ++a) {
            for (int q = rp[size_t(a)]; q < rp[size_t(a) + 1]; ++q) mark[size_t(ci[size_t(q)])] = q;
            for (int e = 0; e < NB; ++e) {
                const int row = a * NB + e;
                for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                    const int c = A.col[p];
                    if (c >= nc) continue;
                    val[size_t(mark[size_t(c / NB)]) * BS + e * NB + c % NB] += A.val[p];
                }
            }
            for (int q = rp[size_t(a)]; q < rp[size_t(a) + 1]; ++q) mark[size_t(ci[size_t(q)])] = -1;
        }
    }

    void buildWells(const CSR& A) {
        cwPtr.assign(size_t(na) + 1, 0); cwW.clear(); cwV.clear();
        for (int a = 0; a < na; ++a) {
            const size_t start = cwW.size();
            for (int e = 0; e < NB; ++e) {
                const int row = a * NB + e;
                for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                    const int c = A.col[p];
                    if (c < nc) continue;
                    const int w = c - nc;
                    size_t k = start;
                    while (k < cwW.size() && cwW[k] != w) ++k;
                    if (k == cwW.size()) { cwW.push_back(w); cwV.resize(cwV.size() + NB, 0.0); }
                    cwV[k * NB + e] += A.val[p];
                }
            }
            cwPtr[size_t(a) + 1] = int(cwW.size());
        }
        dwPtr.assign(size_t(nw) + 1, 0); dwA.clear(); dwV.clear(); ew.assign(size_t(nw), 1.0);
        for (int w = 0; w < nw; ++w) {
            const int row = nc + w;
            const size_t start = dwA.size();
            for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                const int c = A.col[p];
                if (c >= nc) { if (c == row) ew[size_t(w)] = A.val[p]; continue; }
                const int a = c / NB;
                size_t k = start;
                while (k < dwA.size() && dwA[k] != a) ++k;
                if (k == dwA.size()) { dwA.push_back(a); dwV.resize(dwV.size() + NB, 0.0); }
                dwV[k * NB + c % NB] += A.val[p];
            }
            if (!(std::fabs(ew[size_t(w)]) > 1e-300) || !std::isfinite(ew[size_t(w)])) ew[size_t(w)] = 1.0;
            dwPtr[size_t(w) + 1] = int(dwA.size());
        }
    }

    // B' = B - (same-cell part of C E^-1 D)
    void foldWellDiagonal() {
        for (int w = 0; w < nw; ++w)
            for (int k = dwPtr[size_t(w)]; k < dwPtr[size_t(w) + 1]; ++k) {
                const int a = dwA[size_t(k)];
                const double* dv = &dwV[size_t(k) * NB];
                for (int m = cwPtr[size_t(a)]; m < cwPtr[size_t(a) + 1]; ++m) {
                    if (cwW[size_t(m)] != w) continue;
                    const double* cv = &cwV[size_t(m) * NB];
                    double* blk = &val[size_t(dpos[size_t(a)]) * BS];
                    for (int e = 0; e < NB; ++e) for (int v = 0; v < NB; ++v) blk[e * NB + v] -= cv[e] * dv[v] / ew[size_t(w)];
                }
            }
    }

    void factorILU() {
        lu = val;
        invD.assign(size_t(na) * BS, 0.0);
        std::vector<int> pos(size_t(na), -1);
        double tmp[BS];
        for (int i = 0; i < na; ++i) {
            for (int q = rp[size_t(i)]; q < rp[size_t(i) + 1]; ++q) pos[size_t(ci[size_t(q)])] = q;
            for (int q = rp[size_t(i)]; q < rp[size_t(i) + 1]; ++q) {
                const int k = ci[size_t(q)];
                if (k >= i) break;
                bmul<NB>(&lu[size_t(q) * BS], &invD[size_t(k) * BS], tmp);
                std::memcpy(&lu[size_t(q) * BS], tmp, sizeof(tmp));
                for (int qq = dpos[size_t(k)] + 1; qq < rp[size_t(k) + 1]; ++qq) {
                    const int j = ci[size_t(qq)];
                    const int pj = pos[size_t(j)];
                    if (pj >= 0) bmulsub<NB>(tmp, &lu[size_t(qq) * BS], &lu[size_t(pj) * BS]);
                }
            }
            double* Di = &invD[size_t(i) * BS];
            if (!binv<NB>(&lu[size_t(dpos[size_t(i)]) * BS], Di)) {
                // singular pivot block: diagonal fallback
                std::fill(Di, Di + BS, 0.0);
                for (int e = 0; e < NB; ++e) {
                    const double d = lu[size_t(dpos[size_t(i)]) * BS + e * NB + e];
                    Di[e * NB + e] = std::fabs(d) > 1e-300 ? 1.0 / d : 1.0;
                }
            }
            for (int q = rp[size_t(i)]; q < rp[size_t(i) + 1]; ++q) pos[size_t(ci[size_t(q)])] = -1;
        }
    }

    void iluApply(const double* r, double* z) const {
        // forward: y = L^-1 r (unit block diagonal)
        for (int i = 0; i < na; ++i) {
            double y[NB];
            for (int e = 0; e < NB; ++e) y[e] = r[i * NB + e];
            for (int q = rp[size_t(i)]; q < dpos[size_t(i)]; ++q) bmvsub<NB>(&lu[size_t(q) * BS], &z[ci[size_t(q)] * NB], y);
            for (int e = 0; e < NB; ++e) z[i * NB + e] = y[e];
        }
        // backward: z = U^-1 y
        for (int i = na - 1; i >= 0; --i) {
            double y[NB];
            for (int e = 0; e < NB; ++e) y[e] = z[i * NB + e];
            for (int q = dpos[size_t(i)] + 1; q < rp[size_t(i) + 1]; ++q) bmvsub<NB>(&lu[size_t(q) * BS], &z[ci[size_t(q)] * NB], y);
            bmv<NB>(&invD[size_t(i) * BS], y, &z[i * NB]);
        }
    }

    // true-IMPES (accumulation blocks) or quasi-IMPES (diagonal block) weights: D^T w = e_p, normalised so that
    // the pressure equation has the scale of the first (pressure) column
    void weights(const double* acc) {
        wgt.assign(size_t(na) * NB, 0.0);
        forN(na, [&](int a) {
            double D[BS], Di[BS];
            if (acc) std::memcpy(D, acc + size_t(a) * BS, sizeof(D));
            else std::memcpy(D, &val[size_t(dpos[size_t(a)]) * BS], sizeof(D));
            double* w = &wgt[size_t(a) * NB];
            if (binv<NB>(D, Di)) {
                // w = D^-T e0 -> w[e] = Di[0][e]
                for (int e = 0; e < NB; ++e) w[e] = Di[0 * NB + e];
                // scale: make sum |w_e D_e0| = |D_00|-like magnitude to keep AMG well scaled
                double s = 0; for (int e = 0; e < NB; ++e) s += std::fabs(w[e]);
                double ref = 0; for (int e = 0; e < NB; ++e) ref += std::fabs(D[e * NB]);
                if (s > 0 && ref > 0) { const double f = 1.0 / s; for (int e = 0; e < NB; ++e) w[e] *= f; }
            } else for (int e = 0; e < NB; ++e) w[e] = 1.0;
        });
    }

    void pressureMatrix(bool withWells) {
        // perforation-perforation couplings of each well: Ap(a,b) -= (w_a . C(a,w)) D(w,b)_p / E_w
        std::vector<std::vector<std::pair<int, double>>> extra;
        if (withWells && nw > 0) {
            extra.resize(size_t(na));
            for (int a = 0; a < na; ++a)
                for (int m = cwPtr[size_t(a)]; m < cwPtr[size_t(a) + 1]; ++m) {
                    const int w = cwW[size_t(m)];
                    double cp = 0;
                    for (int e = 0; e < NB; ++e) cp += wgt[size_t(a) * NB + e] * cwV[size_t(m) * NB + e];
                    if (cp == 0) continue;
                    for (int k = dwPtr[size_t(w)]; k < dwPtr[size_t(w) + 1]; ++k) {
                        const int b = dwA[size_t(k)];
                        if (b == a) continue;   // already in B'
                        extra[size_t(a)].push_back({b, -cp * dwV[size_t(k) * NB + 0] / ew[size_t(w)]});
                    }
                }
        }
        Ap = CSR{}; Ap.n = na; Ap.rowPtr.assign(size_t(na) + 1, 0);
        std::vector<int> cols; std::vector<double> vals;
        cols.reserve(ci.size() + 64); vals.reserve(ci.size() + 64);
        std::vector<int> pos(size_t(na), -1);
        for (int a = 0; a < na; ++a) {
            const size_t start = cols.size();
            for (int q = rp[size_t(a)]; q < rp[size_t(a) + 1]; ++q) {
                double s = 0;
                const double* blk = &val[size_t(q) * BS];
                for (int e = 0; e < NB; ++e) s += wgt[size_t(a) * NB + e] * blk[e * NB + 0];
                pos[size_t(ci[size_t(q)])] = int(cols.size());
                cols.push_back(ci[size_t(q)]); vals.push_back(s);
            }
            if (!extra.empty())
                for (auto& [b, v] : extra[size_t(a)]) {
                    if (pos[size_t(b)] < 0) { pos[size_t(b)] = int(cols.size()); cols.push_back(b); vals.push_back(0.0); }
                    vals[size_t(pos[size_t(b)])] += v;
                }
            // keep rows sorted (AMG/GS do not need it, but the diagonal lookup is cheaper)
            if (!extra.empty() && !extra[size_t(a)].empty()) {
                std::vector<std::pair<int, double>> row;
                for (size_t k = start; k < cols.size(); ++k) row.push_back({cols[k], vals[k]});
                std::sort(row.begin(), row.end());
                for (size_t k = 0; k < row.size(); ++k) { cols[start + k] = row[k].first; vals[start + k] = row[k].second; }
            }
            for (size_t k = start; k < cols.size(); ++k) pos[size_t(cols[k])] = -1;
            Ap.rowPtr[size_t(a) + 1] = int(cols.size());
        }
        Ap.col = std::move(cols); Ap.val = std::move(vals);
        Ap.finalizeDiag();
    }

    // z = M^-1 r for the cell system (two-stage CPR or ILU only)
    void cellPrec(const double* r, double* z) {
        if (accel) {
            if (accel->apply(r, z)) return;
            accel = nullptr;
            factorILU(); // device failure: fall back to the CPU for this solve
        }
        if (!useCpr) { iluApply(r, z); return; }
        forN(na, [&](int a) {
            double s = 0;
            for (int e = 0; e < NB; ++e) s += wgt[size_t(a) * NB + e] * r[a * NB + e];
            rp_[size_t(a)] = s;
        });
        amg.apply(rp_.data(), xp_.data());
        // t = r - B' (pressure correction)
        forN(na, [&](int a) {
            double y[NB];
            for (int e = 0; e < NB; ++e) y[e] = r[a * NB + e];
            for (int q = rp[size_t(a)]; q < rp[size_t(a) + 1]; ++q) {
                const double xp = xp_[size_t(ci[size_t(q)])];
                const double* blk = &val[size_t(q) * BS];
                for (int e = 0; e < NB; ++e) y[e] -= blk[e * NB + 0] * xp;
            }
            for (int e = 0; e < NB; ++e) t_[size_t(a) * NB + e] = y[e];
        });
        iluApply(t_.data(), z);
        for (int a = 0; a < na; ++a) z[a * NB] += xp_[size_t(a)];
    }

    // full system (cells + wells): block-triangular Schur preconditioner
    void apply(const double* r, double* z) override {
        // rc' = rc - C E^-1 rw
        std::memcpy(rc_.data(), r, sizeof(double) * size_t(nc));
        if (nw > 0)
            for (int a = 0; a < na; ++a)
                for (int m = cwPtr[size_t(a)]; m < cwPtr[size_t(a) + 1]; ++m) {
                    const int w = cwW[size_t(m)];
                    const double f = r[nc + w] / ew[size_t(w)];
                    for (int e = 0; e < NB; ++e) rc_[size_t(a) * NB + e] -= cwV[size_t(m) * NB + e] * f;
                }
        cellPrec(rc_.data(), z);
        for (int w = 0; w < nw; ++w) {
            double s = r[nc + w];
            for (int k = dwPtr[size_t(w)]; k < dwPtr[size_t(w) + 1]; ++k) {
                const int a = dwA[size_t(k)];
                for (int v = 0; v < NB; ++v) s -= dwV[size_t(k) * NB + v] * z[a * NB + v];
            }
            z[nc + w] = s / ew[size_t(w)];
        }
        for (int i = nc + nw; i < N; ++i) z[i] = r[i];
    }
};

} // namespace

struct BlockSolver::Impl {
    std::unique_ptr<CoreBase> core;
    int nb = 0;
};

BlockSolver::BlockSolver() : p_(new Impl) {}
BlockSolver::~BlockSolver() { delete p_; }

BlockSolverStats BlockSolver::solve(const CSR& A, int na, int nb, int nw, const double* acc, const std::vector<double>& b,
                                    std::vector<double>& x, const BlockSolverOptions& o) {
    BlockSolverStats st;
    const auto t0 = Clock::now();
    if (!p_->core || p_->nb != nb) {
        p_->nb = nb;
        switch (nb) {
            case 1: p_->core.reset(new Core<1>()); break;
            case 2: p_->core.reset(new Core<2>()); break;
            case 3: p_->core.reset(new Core<3>()); break;
            default: p_->core.reset(new Core<4>()); break;
        }
    }
    CoreBase& M = *p_->core;
    M.setup(A, na, nw, acc, o);
    st.amgLevels = M.amgLevels();
    st.setupSec = secs(t0);
    const auto t1 = Clock::now();
    const int n = A.n;
    x.assign(size_t(n), 0.0);
    const double bn = std::sqrt(dotv(b.data(), b.data(), n));
    if (!(bn > 0)) { st.converged = true; st.solveSec = secs(t1); return st; }

    // ---- right-preconditioned BiCGStab
    const size_t un = static_cast<size_t>(n);
    std::vector<double> r(b), rh(b), p(un, 0.0), v(un, 0.0), s(un, 0.0), t(un, 0.0), ph(un, 0.0), sh(un, 0.0);
    double rho = 1, alpha = 1, omega = 1;
    double best = 1.0;
    std::vector<double> xbest;
    bool breakdown = false;
    const bool useGmres = o.krylov == 2 || (o.krylov == 0 && (o.cpr || o.accelerator));
    for (int it = 1; it <= (useGmres ? 0 : o.maxIt); ++it) {
        const double rhoN = dotv(rh.data(), r.data(), n);
        if (!(std::fabs(rhoN) > 1e-300)) { breakdown = true; break; }
        const double beta = (rhoN / rho) * (alpha / omega);
        rho = rhoN;
        forN(n, [&](int i) { p[size_t(i)] = r[size_t(i)] + beta * (p[size_t(i)] - omega * v[size_t(i)]); });
        M.apply(p.data(), ph.data());
        A.spmv(ph.data(), v.data());
        const double rv = dotv(rh.data(), v.data(), n);
        if (!(std::fabs(rv) > 1e-300)) { breakdown = true; break; }
        alpha = rho / rv;
        forN(n, [&](int i) { s[size_t(i)] = r[size_t(i)] - alpha * v[size_t(i)]; });
        const double sn = std::sqrt(dotv(s.data(), s.data(), n));
        if (sn / bn < o.tol) {
            forN(n, [&](int i) { x[size_t(i)] += alpha * ph[size_t(i)]; });
            st.converged = true; st.iterations = it; st.relResidual = sn / bn; st.solveSec = secs(t1); return st;
        }
        M.apply(s.data(), sh.data());
        A.spmv(sh.data(), t.data());
        const double tt = dotv(t.data(), t.data(), n);
        omega = tt > 0 ? dotv(t.data(), s.data(), n) / tt : 0;
        forN(n, [&](int i) { x[size_t(i)] += alpha * ph[size_t(i)] + omega * sh[size_t(i)]; r[size_t(i)] = s[size_t(i)] - omega * t[size_t(i)]; });
        const double rn = std::sqrt(dotv(r.data(), r.data(), n));
        st.iterations = it; st.relResidual = rn / bn;
        if (!std::isfinite(rn)) { breakdown = true; break; }
        if (rn / bn < best) { best = rn / bn; if (best < 0.5) xbest = x; }
        if (rn / bn < o.tol) { st.converged = true; st.solveSec = secs(t1); return st; }
        if (omega == 0) { breakdown = true; break; }
    }
    // ---- restarted flexible GMRES from the best BiCGStab iterate
    if (!xbest.empty()) x = xbest; else x.assign(size_t(n), 0.0);
    (void)breakdown;
    const int m = std::max(5, o.restart);
    const size_t um = static_cast<size_t>(m);
    std::vector<std::vector<double>> V(um + 1, std::vector<double>(un, 0.0)), Z(um, std::vector<double>(un, 0.0));
    std::vector<double> H((um + 1) * um, 0.0), cs(um, 0.0), sn(um, 0.0), g(um + 1, 0.0), w(un, 0.0);
    auto Hat = [&](int i, int j) -> double& { return H[size_t(i) * m + j]; };
    int total = 0;
    const int budget = o.maxIt;
    while (total < budget) {
        A.spmv(x.data(), w.data());
        forN(n, [&](int i) { r[size_t(i)] = b[size_t(i)] - w[size_t(i)]; });
        const double beta = std::sqrt(dotv(r.data(), r.data(), n));
        st.relResidual = beta / bn;
        if (beta / bn < o.tol || !std::isfinite(beta)) break;
        forN(n, [&](int i) { V[0][size_t(i)] = r[size_t(i)] / beta; });
        std::fill(g.begin(), g.end(), 0.0); g[0] = beta;
        int j = 0;
        for (; j < m && total < budget; ++j, ++total) {
            M.apply(V[size_t(j)].data(), Z[size_t(j)].data());
            A.spmv(Z[size_t(j)].data(), w.data());
            for (int i = 0; i <= j; ++i) {           // modified Gram-Schmidt
                const double h = dotv(w.data(), V[size_t(i)].data(), n);
                Hat(i, j) = h;
                const double* vi = V[size_t(i)].data();
                forN(n, [&](int q) { w[size_t(q)] -= h * vi[q]; });
            }
            const double hn = std::sqrt(dotv(w.data(), w.data(), n));
            Hat(j + 1, j) = hn;
            if (hn > 0) { double* vj = V[size_t(j) + 1].data(); forN(n, [&](int q) { vj[q] = w[size_t(q)] / hn; }); }
            for (int i = 0; i < j; ++i) {
                const double tq = cs[size_t(i)] * Hat(i, j) + sn[size_t(i)] * Hat(i + 1, j);
                Hat(i + 1, j) = -sn[size_t(i)] * Hat(i, j) + cs[size_t(i)] * Hat(i + 1, j);
                Hat(i, j) = tq;
            }
            const double den = std::hypot(Hat(j, j), Hat(j + 1, j));
            cs[size_t(j)] = den > 0 ? Hat(j, j) / den : 1; sn[size_t(j)] = den > 0 ? Hat(j + 1, j) / den : 0;
            Hat(j, j) = den; Hat(j + 1, j) = 0;
            g[size_t(j) + 1] = -sn[size_t(j)] * g[size_t(j)]; g[size_t(j)] = cs[size_t(j)] * g[size_t(j)];
            st.relResidual = std::fabs(g[size_t(j) + 1]) / bn;
            if (st.relResidual < o.tol) { ++j; ++total; break; }
        }
        std::vector<double> y(static_cast<size_t>(j), 0.0);
        for (int i = j - 1; i >= 0; --i) {
            double sum = g[size_t(i)];
            for (int k = i + 1; k < j; ++k) sum -= Hat(i, k) * y[size_t(k)];
            y[size_t(i)] = Hat(i, i) != 0 ? sum / Hat(i, i) : 0;
        }
        for (int i = 0; i < j; ++i) { const double yi = y[size_t(i)]; const double* zi = Z[size_t(i)].data(); forN(n, [&](int q) { x[size_t(q)] += yi * zi[q]; }); }
        if (st.relResidual < o.tol) break;
    }
    // true residual
    A.spmv(x.data(), w.data());
    double rr = 0; for (int i = 0; i < n; ++i) { const double d = b[size_t(i)] - w[size_t(i)]; rr += d * d; }
    st.relResidual = std::sqrt(rr) / bn;
    st.converged = st.relResidual < o.tol * 1.0001;
    st.iterations += total;
    st.solveSec = secs(t1);
    return st;
}

} // namespace ressim
