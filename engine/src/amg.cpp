// Aggregation AMG (plain by default, optionally smoothed) and the two-stage CPR (Constrained Pressure Residual) preconditioner.
#include "ressim/linsolve.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

namespace ressim {

namespace {
const int kPar = 20000;          // parallel threshold for level kernels
const int kMaxCoarse = 400;      // dense LU below this size
const double kTheta = 0.08;      // strength threshold
// Plain (unsmoothed) aggregation by default, as in OPM Flow's CPR: the pressure matrix is non-symmetric and thin
// layers make it strongly anisotropic; prolongator smoothing then produces dense coarse operators (and sometimes
// non-positive coarse diagonals). Measured: SPE9 -30 %, thin-layer corner-point models 3-14x faster.
const double kRowTheta = 0.25;   // also strong: |a_ij| >= kRowTheta * max_k |a_ik| (every row keeps a strong neighbour)
const bool kSmooth = false;      // smoothed aggregation (with the safeguards below) when true
const double kFill = 2.0;        // smoothed level rejected when coarse nnz/row > kFill * max(8, fine nnz/row)

struct Sp { int nr = 0, nc = 0; std::vector<int> rp, ci; std::vector<double> v; };

Sp fromCSR(const CSR& A) { Sp s; s.nr = s.nc = A.n; s.rp = A.rowPtr; s.ci = A.col; s.v = A.val; return s; }

Sp transpose(const Sp& A) {
    Sp T; T.nr = A.nc; T.nc = A.nr; T.rp.assign(T.nr + 1, 0);
    for (int c : A.ci) T.rp[c + 1]++;
    for (int i = 0; i < T.nr; ++i) T.rp[i + 1] += T.rp[i];
    T.ci.resize(A.ci.size()); T.v.resize(A.v.size());
    std::vector<int> fill(T.rp.begin(), T.rp.end() - 1);
    for (int i = 0; i < A.nr; ++i)
        for (int p = A.rp[i]; p < A.rp[i + 1]; ++p) { const int q = fill[A.ci[p]]++; T.ci[q] = i; T.v[q] = A.v[p]; }
    return T;
}

Sp multiply(const Sp& A, const Sp& B) {
    Sp C; C.nr = A.nr; C.nc = B.nc; C.rp.assign(A.nr + 1, 0);
    std::vector<int> mark(B.nc, -1);
    // symbolic
    for (int i = 0; i < A.nr; ++i) {
        int cnt = 0;
        for (int p = A.rp[i]; p < A.rp[i + 1]; ++p)
            for (int q = B.rp[A.ci[p]]; q < B.rp[A.ci[p] + 1]; ++q)
                if (mark[B.ci[q]] != i) { mark[B.ci[q]] = i; ++cnt; }
        C.rp[i + 1] = C.rp[i] + cnt;
    }
    C.ci.resize(C.rp[A.nr]); C.v.assign(C.rp[A.nr], 0.0);
    std::fill(mark.begin(), mark.end(), -1);
    for (int i = 0; i < A.nr; ++i) {
        int pos = C.rp[i];
        const int start = pos;
        for (int p = A.rp[i]; p < A.rp[i + 1]; ++p) {
            const double a = A.v[p];
            for (int q = B.rp[A.ci[p]]; q < B.rp[A.ci[p] + 1]; ++q) {
                const int j = B.ci[q];
                if (mark[j] < start) { mark[j] = pos; C.ci[pos] = j; C.v[pos] = a * B.v[q]; ++pos; }
                else C.v[mark[j]] += a * B.v[q];
            }
        }
        // sort the row by column
        std::vector<std::pair<int, double>> row;
        row.reserve(pos - start);
        for (int p = start; p < pos; ++p) row.push_back({C.ci[p], C.v[p]});
        std::sort(row.begin(), row.end());
        for (int p = start; p < pos; ++p) { C.ci[p] = row[p - start].first; C.v[p] = row[p - start].second; }
    }
    return C;
}

CSR toCSR(const Sp& S) {
    CSR A; A.n = S.nr; A.rowPtr = S.rp; A.col = S.ci; A.val = S.v;
    // guarantee a diagonal entry in every row
    bool missing = false;
    for (int i = 0; i < A.n && !missing; ++i) {
        bool f = false;
        for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (A.col[p] == i) { f = true; break; }
        if (!f) missing = true;
    }
    if (missing) {
        CSR B; B.n = A.n; B.rowPtr.assign(A.n + 1, 0);
        for (int i = 0; i < A.n; ++i) {
            bool f = false;
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) {
                if (!f && A.col[p] > i) { B.col.push_back(i); B.val.push_back(0.0); f = true; }
                if (A.col[p] == i) f = true;
                B.col.push_back(A.col[p]); B.val.push_back(A.val[p]);
            }
            if (!f) { B.col.push_back(i); B.val.push_back(0.0); }
            B.rowPtr[i + 1] = int(B.col.size());
        }
        A = std::move(B);
    }
    A.finalizeDiag();
    for (int i = 0; i < A.n; ++i) if (std::fabs(A.val[A.diag[i]]) < 1e-300) A.val[A.diag[i]] = 1.0;
    return A;
}

void gaussSeidel(const CSR& A, const double* b, double* x, bool forward) {
    const int n = A.n;
    if (forward) {
        for (int i = 0; i < n; ++i) {
            double s = b[i];
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (A.col[p] != i) s -= A.val[p] * x[A.col[p]];
            x[i] = s / A.val[A.diag[i]];
        }
    } else {
        for (int i = n - 1; i >= 0; --i) {
            double s = b[i];
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (A.col[p] != i) s -= A.val[p] * x[A.col[p]];
            x[i] = s / A.val[A.diag[i]];
        }
    }
}

void spmvRect(const std::vector<int>& rp, const std::vector<int>& ci, const std::vector<double>& v, int nr,
              const double* x, double* y, bool add) {
    if (useParallel(nr > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nr; ++i) {
            double s = 0;
            for (int p = rp[i]; p < rp[i + 1]; ++p) s += v[p] * x[ci[p]];
            y[i] = add ? y[i] + s : s;
        }
    } else {
        for (int i = 0; i < nr; ++i) {
            double s = 0;
            for (int p = rp[i]; p < rp[i + 1]; ++p) s += v[p] * x[ci[p]];
            y[i] = add ? y[i] + s : s;
        }
    }
}
}  // namespace

// ======================================================================== AMG setup
void AMG::setup(const CSR& A0) { setup(CSR(A0)); }

void AMG::setup(CSR&& A0) {
    lv_.clear();
    lv_.emplace_back();
    lv_[0].A = std::move(A0);
    lv_[0].A.finalizeDiag();
    for (int i = 0; i < lv_[0].A.n; ++i) if (lv_[0].A.diag[i] < 0) { CSR t = toCSR(fromCSR(lv_[0].A)); lv_[0].A = std::move(t); break; }
    for (int i = 0; i < lv_[0].A.n; ++i) if (std::fabs(lv_[0].A.val[lv_[0].A.diag[i]]) < 1e-300) lv_[0].A.val[lv_[0].A.diag[i]] = 1.0;

    while (true) {
        Level& L = lv_.back();
        const CSR& A = L.A;
        const int n = A.n;
        if (n <= kMaxCoarse || lv_.size() >= 14) break;
        // ---- strength of connection
        std::vector<double> d(n);
        for (int i = 0; i < n; ++i) d[i] = A.val[A.diag[i]];
        std::vector<int> sp(n + 1, 0), sc;
        sc.reserve(A.col.size());
        for (int i = 0; i < n; ++i) {
            double rowMax = 0;
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (A.col[p] != i) rowMax = std::max(rowMax, std::fabs(A.val[p]));
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) {
                const int j = A.col[p];
                if (j == i) continue;
                if (std::fabs(A.val[p]) >= kTheta * std::sqrt(std::fabs(d[i] * d[j])) || std::fabs(A.val[p]) >= kRowTheta * rowMax) sc.push_back(j);
            }
            sp[i + 1] = int(sc.size());
        }
        // ---- greedy aggregation
        std::vector<int> agg(n, -1);
        int na = 0;
        for (int i = 0; i < n; ++i) {                       // pass 1: roots with free neighbourhoods
            if (agg[i] >= 0) continue;
            bool free = true;
            for (int p = sp[i]; p < sp[i + 1]; ++p) if (agg[sc[p]] >= 0) { free = false; break; }
            if (!free) continue;
            if (sp[i] == sp[i + 1]) continue;               // isolated: handled later
            agg[i] = na;
            for (int p = sp[i]; p < sp[i + 1]; ++p) agg[sc[p]] = na;
            ++na;
        }
        std::vector<int> agg1 = agg;
        for (int i = 0; i < n; ++i) {                       // pass 2: join a neighbouring aggregate
            if (agg[i] >= 0) continue;
            double best = 0; int ba = -1;
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) {
                const int j = A.col[p];
                if (j != i && agg1[j] >= 0 && std::fabs(A.val[p]) > best) { best = std::fabs(A.val[p]); ba = agg1[j]; }
            }
            if (ba >= 0 && sp[i] != sp[i + 1]) agg[i] = ba;
        }
        for (int i = 0; i < n; ++i) {                       // pass 3: leftovers form new aggregates
            if (agg[i] >= 0) continue;
            agg[i] = na;
            for (int p = sp[i]; p < sp[i + 1]; ++p) if (agg[sc[p]] < 0) agg[sc[p]] = na;
            ++na;
        }
        if (na >= int(0.85 * n) || na < 1) break;           // coarsening stalled
        // ---- spectral radius of D^-1 A (power iteration)
        std::vector<double> x(n), y(n);
        std::mt19937 rng(7); std::uniform_real_distribution<double> U(0.5, 1.5);
        for (double& v : x) v = U(rng);
        double rho = 1.0;
        for (int it = 0; kSmooth && it < 12; ++it) {
            A.spmv(x.data(), y.data());
            double nx = 0, ny = 0;
            for (int i = 0; i < n; ++i) { y[i] /= d[i]; nx += x[i] * x[i]; ny += y[i] * y[i]; }
            rho = std::sqrt(ny / std::max(nx, 1e-300));
            const double s = 1.0 / std::sqrt(std::max(ny, 1e-300));
            for (int i = 0; i < n; ++i) x[i] = y[i] * s;
        }
        // ---- (smoothed) prolongator  P = (I - omega D^-1 A) T
        // The pressure matrix is non-symmetric (true-IMPES weights + well rows); smoothing can then produce coarse
        // operators with non-positive diagonals on which Gauss-Seidel diverges. Such a level falls back to plain
        // (unsmoothed) aggregation, as used by OPM Flow's CPR.
        std::vector<double> dU(na, 0.0);                     // diagonal of the unsmoothed coarse operator
        for (int i = 0; i < n; ++i)
            for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (agg[A.col[p]] == agg[i]) dU[agg[i]] += A.val[p];
        Sp P, R, Ac;
        auto build = [&](double omega) {
            P = Sp{}; P.nr = n; P.nc = na; P.rp.assign(n + 1, 0);
            std::vector<double> acc(na, 0.0); std::vector<int> mark(na, -1), cols;
            for (int i = 0; i < n; ++i) {
                cols.clear();
                auto add = [&](int a, double v) { if (mark[a] != i) { mark[a] = i; acc[a] = 0; cols.push_back(a); } acc[a] += v; };
                add(agg[i], 1.0);
                if (omega != 0.0) {
                    const double f = omega / d[i];
                    for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) add(agg[A.col[p]], -f * A.val[p]);
                }
                std::sort(cols.begin(), cols.end());
                for (int a : cols) if (acc[a] != 0.0) { P.ci.push_back(a); P.v.push_back(acc[a]); }
                P.rp[i + 1] = int(P.ci.size());
            }
            R = transpose(P);
            Ac = multiply(R, multiply(fromCSR(A), P));
        };
        build(kSmooth ? 4.0 / (3.0 * std::max(rho, 1e-12)) : 0.0);
        if (kSmooth) {
            bool bad = double(Ac.ci.size()) / na > kFill * std::max(8.0, double(A.col.size()) / n);   // fill explosion
            for (int i = 0; i < na && !bad; ++i) {
                double dc = 0;
                for (int p = Ac.rp[i]; p < Ac.rp[i + 1]; ++p) if (Ac.ci[p] == i) { dc = Ac.v[p]; break; }
                if (dc * dU[i] <= 0 || std::fabs(dc) < 1e-3 * std::fabs(dU[i])) bad = true;
            }
            if (bad) build(0.0);
        }
        L.nc = na;
        L.pRow = P.rp; L.pCol = P.ci; L.pVal = P.v;
        L.rRow = R.rp; L.rCol = R.ci; L.rVal = R.v;
        Level next; next.A = toCSR(Ac);
        lv_.push_back(std::move(next));
    }
    for (auto& L : lv_) { L.x.assign(L.A.n, 0.0); L.b.assign(L.A.n, 0.0); L.r.assign(L.A.n, 0.0); }
    factorCoarsest();
}

// New values on the same sparsity pattern (next Newton iteration / time step): keep the aggregates and
// transfer operators, recompute the Galerkin coarse operators and the coarsest LU. Unsmoothed aggregation
// transfer operators depend only on the aggregates, so this is exact apart from the (older) aggregation.
bool AMG::refresh(const CSR& A0) {
    if (lv_.empty() || A0.n != lv_[0].A.n || A0.col.size() != lv_[0].A.col.size() || A0.col != lv_[0].A.col) return false;
    lv_[0].A.val = A0.val;
    for (int i = 0; i < lv_[0].A.n; ++i) if (std::fabs(lv_[0].A.val[lv_[0].A.diag[i]]) < 1e-300) lv_[0].A.val[lv_[0].A.diag[i]] = 1.0;
    for (size_t l = 0; l + 1 < lv_.size(); ++l) {
        const Level& L = lv_[l];
        Sp P; P.nr = L.A.n; P.nc = L.nc; P.rp = L.pRow; P.ci = L.pCol; P.v = L.pVal;
        Sp R; R.nr = L.nc; R.nc = L.A.n; R.rp = L.rRow; R.ci = L.rCol; R.v = L.rVal;
        CSR Ac = toCSR(multiply(R, multiply(fromCSR(L.A), P)));
        if (Ac.col != lv_[l + 1].A.col) return false;
        lv_[l + 1].A.val = Ac.val;
    }
    factorCoarsest();
    return true;
}

void AMG::factorCoarsest() {
    // ---- dense LU of the coarsest level (if small enough)
    const CSR& C = lv_.back().A;
    lu_.clear(); piv_.clear();
    if (C.n <= 2 * kMaxCoarse) {
        const int m = C.n;
        lu_.assign(size_t(m) * m, 0.0); piv_.resize(m);
        for (int i = 0; i < m; ++i) for (int p = C.rowPtr[i]; p < C.rowPtr[i + 1]; ++p) lu_[size_t(i) * m + C.col[p]] = C.val[p];
        for (int k = 0; k < m; ++k) {
            int pr = k; double mx = std::fabs(lu_[size_t(k) * m + k]);
            for (int i = k + 1; i < m; ++i) if (std::fabs(lu_[size_t(i) * m + k]) > mx) { mx = std::fabs(lu_[size_t(i) * m + k]); pr = i; }
            piv_[k] = pr;
            if (pr != k) for (int j = 0; j < m; ++j) std::swap(lu_[size_t(k) * m + j], lu_[size_t(pr) * m + j]);
            double& dk = lu_[size_t(k) * m + k];
            if (std::fabs(dk) < 1e-300) dk = 1e-300;
            for (int i = k + 1; i < m; ++i) {
                const double f = lu_[size_t(i) * m + k] / dk;
                lu_[size_t(i) * m + k] = f;
                if (f != 0) for (int j = k + 1; j < m; ++j) lu_[size_t(i) * m + j] -= f * lu_[size_t(k) * m + j];
            }
        }
    }
}

double AMG::operatorComplexity() const {
    if (lv_.empty()) return 0;
    double s = 0; for (auto& L : lv_) s += double(L.A.col.size());
    return s / double(lv_[0].A.col.size());
}

void AMG::cycle(int l) const {
    const Level& L = lv_[l];
    const int n = L.A.n;
    if (l + 1 == int(lv_.size())) {
        if (!lu_.empty()) {
            const int m = n;
            std::vector<double>& x = L.x;
            x = L.b;
            for (int k = 0; k < m; ++k) if (piv_[k] != k) std::swap(x[k], x[piv_[k]]);
            for (int i = 0; i < m; ++i) { double s = x[i]; for (int j = 0; j < i; ++j) s -= lu_[size_t(i) * m + j] * x[j]; x[i] = s; }
            for (int i = m - 1; i >= 0; --i) { double s = x[i]; for (int j = i + 1; j < m; ++j) s -= lu_[size_t(i) * m + j] * x[j]; x[i] = s / lu_[size_t(i) * m + i]; }
        } else {
            std::fill(L.x.begin(), L.x.end(), 0.0);
            for (int s = 0; s < 20; ++s) { gaussSeidel(L.A, L.b.data(), L.x.data(), true); gaussSeidel(L.A, L.b.data(), L.x.data(), false); }
        }
        return;
    }
    std::fill(L.x.begin(), L.x.end(), 0.0);
    gaussSeidel(L.A, L.b.data(), L.x.data(), true);
    L.A.spmv(L.x.data(), L.r.data());
    for (int i = 0; i < n; ++i) L.r[i] = L.b[i] - L.r[i];
    const Level& C = lv_[l + 1];
    spmvRect(L.rRow, L.rCol, L.rVal, L.nc, L.r.data(), C.b.data(), false);
    cycle(l + 1);
    spmvRect(L.pRow, L.pCol, L.pVal, n, C.x.data(), L.x.data(), true);
    gaussSeidel(L.A, L.b.data(), L.x.data(), false);
}

void AMG::apply(const double* r, double* z) const {
    if (lv_.empty()) return;
    const Level& L = lv_[0];
    std::copy(r, r + L.A.n, L.b.begin());
    cycle(0);
    std::copy(L.x.begin(), L.x.end(), z);
}

// ======================================================================== CPR
// condition estimate of a small dense block: ||D||_1 * ||D^-1||_1 (Gauss-Jordan; huge when singular)
static double blockCondition(const double D[3][3], int n) {
    double M[3][6] = {{0}};
    double n1 = 0;
    for (int j = 0; j < n; ++j) { double c = 0; for (int i = 0; i < n; ++i) c += std::fabs(D[i][j]); n1 = std::max(n1, c); }
    for (int i = 0; i < n; ++i) { for (int j = 0; j < n; ++j) M[i][j] = D[i][j]; M[i][n + i] = 1.0; }
    for (int k = 0; k < n; ++k) {
        int pr = k;
        for (int i = k + 1; i < n; ++i) if (std::fabs(M[i][k]) > std::fabs(M[pr][k])) pr = i;
        if (std::fabs(M[pr][k]) <= 1e-300 * std::max(1.0, n1)) return 1e300;
        if (pr != k) for (int j = 0; j < 2 * n; ++j) std::swap(M[k][j], M[pr][j]);
        const double piv = M[k][k];
        for (int j = 0; j < 2 * n; ++j) M[k][j] /= piv;
        for (int i = 0; i < n; ++i) if (i != k) { const double f = M[i][k]; for (int j = 0; j < 2 * n; ++j) M[i][j] -= f * M[k][j]; }
    }
    double ni = 0;
    for (int j = 0; j < n; ++j) { double c = 0; for (int i = 0; i < n; ++i) c += std::fabs(M[i][n + j]); ni = std::max(ni, c); }
    return n1 * ni;
}

void CPR::setup(const CSR& A, int nCells, int nv, int nWells, const double* accBlocks) {
    A_ = &A; na_ = nCells; nv_ = nv; nw_ = nWells;
    const int nc = na_ * nv_;
    // quasi-IMPES weights: D^T w = e_p for the diagonal block of each cell
    w_.assign(size_t(nc), 1.0);
    nFallback_ = 0;
    for (int a = 0; a < na_; ++a) {
        double D[3][3] = {{0}};
        auto jacBlock = [&] {
            for (int e = 0; e < nv_; ++e) {
                const int row = a * nv_ + e;
                const int p0 = A.diag[row] - e;              // position of column a*nv in this row
                for (int v = 0; v < nv_; ++v) D[e][v] = A.val[p0 + v];
            }
        };
        if (accBlocks) {
            const double* B = accBlocks + size_t(a) * nv_ * nv_;
            for (int e = 0; e < nv_; ++e) for (int v = 0; v < nv_; ++v) D[e][v] = B[e * nv_ + v];
            // true-IMPES needs a well-conditioned accumulation block. In cells where a phase is absent (gas cap
            // without oil, water-only cells) it is close to singular and its weights wreck the pressure system:
            // use the quasi-IMPES (Jacobian diagonal block) weights there.
            if (nv_ > 1 && blockCondition(D, nv_) > 1e10) { jacBlock(); ++nFallback_; }
        } else jacBlock();
        // solve D^T w = e0 (Gaussian elimination with pivoting on the transposed block)
        double M[3][4] = {{0}};
        for (int i = 0; i < nv_; ++i) { for (int j = 0; j < nv_; ++j) M[i][j] = D[j][i]; M[i][nv_] = (i == 0) ? 1.0 : 0.0; }
        bool ok = true;
        for (int k = 0; k < nv_ && ok; ++k) {
            int pr = k;
            for (int i = k + 1; i < nv_; ++i) if (std::fabs(M[i][k]) > std::fabs(M[pr][k])) pr = i;
            if (std::fabs(M[pr][k]) < 1e-14) { ok = false; break; }
            if (pr != k) for (int j = 0; j <= nv_; ++j) std::swap(M[k][j], M[pr][j]);
            for (int i = 0; i < nv_; ++i) {
                if (i == k) continue;
                const double f = M[i][k] / M[k][k];
                for (int j = k; j <= nv_; ++j) M[i][j] -= f * M[k][j];
            }
        }
        if (ok) for (int e = 0; e < nv_; ++e) w_[size_t(a) * nv_ + e] = M[e][nv_] / M[e][e];
        else {
            double s = 0; for (int e = 0; e < nv_; ++e) s += D[e][0];
            for (int e = 0; e < nv_; ++e) w_[size_t(a) * nv_ + e] = std::fabs(s) > 1e-300 ? 1.0 / s : 1.0;
        }
        // normalise: the pressure row keeps the scale of the largest weight (as OPM does); scaling a row of the
        // pressure system does not change its solution but keeps AMG strength measures comparable across cells
        double wm = 0; for (int e = 0; e < nv_; ++e) wm = std::max(wm, std::fabs(w_[size_t(a) * nv_ + e]));
        if (!(wm > 0) || !std::isfinite(wm)) { for (int e = 0; e < nv_; ++e) w_[size_t(a) * nv_ + e] = 1.0; }
    }
    // pressure matrix (cells + wells)
    const int np = na_ + nw_;
    Ap_ = CSR{}; Ap_.n = np; Ap_.rowPtr.assign(np + 1, 0);
    auto pcol = [&](int c) { return c < nc ? ((c % nv_) == 0 ? c / nv_ : -1) : na_ + (c - nc); };
    for (int a = 0; a < na_; ++a) {
        const int row = a * nv_;
        int cnt = 0;
        for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) if (pcol(A.col[p]) >= 0) ++cnt;
        Ap_.rowPtr[a + 1] = cnt;
    }
    for (int w = 0; w < nw_; ++w) {
        const int row = nc + w; int cnt = 0;
        for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) if (pcol(A.col[p]) >= 0) ++cnt;
        Ap_.rowPtr[na_ + w + 1] = cnt;
    }
    for (int i = 0; i < np; ++i) Ap_.rowPtr[i + 1] += Ap_.rowPtr[i];
    Ap_.col.resize(Ap_.rowPtr[np]); Ap_.val.assign(Ap_.rowPtr[np], 0.0);
    if (useParallel(na_ > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int a = 0; a < na_; ++a) {
            const int row = a * nv_;
            int q = Ap_.rowPtr[a];
            for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                const int c = pcol(A.col[p]);
                if (c < 0) continue;
                const int off = p - A.rowPtr[row];
                double s = 0;
                for (int e = 0; e < nv_; ++e) s += w_[size_t(a) * nv_ + e] * A.val[A.rowPtr[row + e] + off];
                Ap_.col[q] = c; Ap_.val[q] = s; ++q;
            }
        }
    } else {
        for (int a = 0; a < na_; ++a) {
            const int row = a * nv_;
            int q = Ap_.rowPtr[a];
            for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
                const int c = pcol(A.col[p]);
                if (c < 0) continue;
                const int off = p - A.rowPtr[row];
                double s = 0;
                for (int e = 0; e < nv_; ++e) s += w_[size_t(a) * nv_ + e] * A.val[A.rowPtr[row + e] + off];
                Ap_.col[q] = c; Ap_.val[q] = s; ++q;
            }
        }
    }
    for (int w = 0; w < nw_; ++w) {
        const int row = nc + w; int q = Ap_.rowPtr[na_ + w];
        for (int p = A.rowPtr[row]; p < A.rowPtr[row + 1]; ++p) {
            const int c = pcol(A.col[p]);
            if (c < 0) continue;
            Ap_.col[q] = c; Ap_.val[q] = A.val[p]; ++q;
        }
    }
    Ap_.finalizeDiag();
    amg_.setup(Ap_);
    ilu_.factor(A);
    rp_.assign(np, 0.0); xp_.assign(np, 0.0); z1_.assign(A.n, 0.0); t_.assign(A.n, 0.0);
}

void CPR::apply(const double* r, double* z) const {
    const int nc = na_ * nv_;
    const int N = A_->n;
    if (useParallel(na_ > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int a = 0; a < na_; ++a) {
            double s = 0;
            for (int e = 0; e < nv_; ++e) s += w_[size_t(a) * nv_ + e] * r[a * nv_ + e];
            rp_[a] = s;
        }
    } else {
        for (int a = 0; a < na_; ++a) {
            double s = 0;
            for (int e = 0; e < nv_; ++e) s += w_[size_t(a) * nv_ + e] * r[a * nv_ + e];
            rp_[a] = s;
        }
    }
    for (int w = 0; w < nw_; ++w) rp_[na_ + w] = r[nc + w];
    amg_.apply(rp_.data(), xp_.data());
    std::fill(z1_.begin(), z1_.end(), 0.0);
    for (int a = 0; a < na_; ++a) z1_[size_t(a) * nv_] = xp_[a];
    for (int w = 0; w < nw_; ++w) z1_[nc + w] = xp_[na_ + w];
    A_->spmv(z1_.data(), t_.data());
    if (useParallel(N > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; ++i) t_[i] = r[i] - t_[i];
    } else {
        for (int i = 0; i < N; ++i) t_[i] = r[i] - t_[i];
    }
    ilu_.apply(t_.data(), z);
    if (useParallel(N > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; ++i) z[i] += z1_[i];
    } else {
        for (int i = 0; i < N; ++i) z[i] += z1_[i];
    }
}

}  // namespace ressim
