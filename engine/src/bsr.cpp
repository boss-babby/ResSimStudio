#include <cstdlib>
#include "ressim/bsr.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace ressim {

namespace {
const int kPar = 6000;      // block rows above which kernels run in parallel

// y += a * x  for nb x nb block a (row-major) and nb-vectors
template <int NB> inline void gemvAdd(const double* a, const double* x, double* y, int) {
    for (int i = 0; i < NB; ++i) { double s = 0; for (int j = 0; j < NB; ++j) s += a[i * NB + j] * x[j]; y[i] += s; }
}
template <int NB> inline void gemvSub(const double* a, const double* x, double* y, int) {
    for (int i = 0; i < NB; ++i) { double s = 0; for (int j = 0; j < NB; ++j) s += a[i * NB + j] * x[j]; y[i] -= s; }
}
inline void gemvAddN(const double* a, const double* x, double* y, int nb) {
    for (int i = 0; i < nb; ++i) { double s = 0; for (int j = 0; j < nb; ++j) s += a[i * nb + j] * x[j]; y[i] += s; }
}
inline void gemvSubN(const double* a, const double* x, double* y, int nb) {
    for (int i = 0; i < nb; ++i) { double s = 0; for (int j = 0; j < nb; ++j) s += a[i * nb + j] * x[j]; y[i] -= s; }
}
// compile-time block sizes (compositional blocks are 4..10 wide): fully unrolled inner loops
template <int NB> inline void gemmT(const double* a, const double* b, double* c) {
    for (int i = 0; i < NB; ++i) {
        double r[NB] = {};
        for (int k = 0; k < NB; ++k) { const double f = a[i * NB + k]; for (int j = 0; j < NB; ++j) r[j] += f * b[k * NB + j]; }
        for (int j = 0; j < NB; ++j) c[i * NB + j] = r[j];
    }
}
template <int NB> inline void gemmSubT(const double* a, const double* b, double* c) {
    for (int i = 0; i < NB; ++i) {
        double r[NB] = {};
        for (int k = 0; k < NB; ++k) { const double f = a[i * NB + k]; for (int j = 0; j < NB; ++j) r[j] += f * b[k * NB + j]; }
        for (int j = 0; j < NB; ++j) c[i * NB + j] -= r[j];
    }
}
#define RESSIM_NB_SWITCH(NBV, F, ...) \
    switch (NBV) { case 2: F<2>(__VA_ARGS__); return; case 3: F<3>(__VA_ARGS__); return; case 4: F<4>(__VA_ARGS__); return; \
                   case 5: F<5>(__VA_ARGS__); return; case 6: F<6>(__VA_ARGS__); return; case 7: F<7>(__VA_ARGS__); return; \
                   case 8: F<8>(__VA_ARGS__); return; case 9: F<9>(__VA_ARGS__); return; case 10: F<10>(__VA_ARGS__); return; default: break; }
// c = a * b
inline void gemm(const double* a, const double* b, double* c, int nb) {
    RESSIM_NB_SWITCH(nb, gemmT, a, b, c)
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j) { double s = 0; for (int k = 0; k < nb; ++k) s += a[i * nb + k] * b[k * nb + j]; c[i * nb + j] = s; }
}
// c -= a * b
inline void gemmSub(const double* a, const double* b, double* c, int nb) {
    RESSIM_NB_SWITCH(nb, gemmSubT, a, b, c)
    for (int i = 0; i < nb; ++i)
        for (int j = 0; j < nb; ++j) { double s = 0; for (int k = 0; k < nb; ++k) s += a[i * nb + k] * b[k * nb + j]; c[i * nb + j] -= s; }
}

// dispatch a block-size templated kernel
template <int NB> inline void gemvAddT(const double* a, const double* x, double* y) { gemvAdd<NB>(a, x, y, NB); }
template <int NB> inline void gemvSubT(const double* a, const double* x, double* y) { gemvSub<NB>(a, x, y, NB); }
inline void gemvAddD(const double* a, const double* x, double* y, int nb) { RESSIM_NB_SWITCH(nb, gemvAddT, a, x, y) gemvAddN(a, x, y, nb); }
inline void gemvSubD(const double* a, const double* x, double* y, int nb) { RESSIM_NB_SWITCH(nb, gemvSubT, a, x, y) gemvSubN(a, x, y, nb); }
#define RESSIM_BLOCK_DISPATCH(NBV, CALL2, CALL3, CALLN) \
    switch (NBV) { case 2: CALL2; break; case 3: CALL3; break; default: CALLN; break; }
} // namespace

bool invertBlock(double* a, int nb) {
    double m[kMaxBlock][2 * kMaxBlock];
    for (int i = 0; i < nb; ++i) for (int j = 0; j < nb; ++j) { m[i][j] = a[i * nb + j]; m[i][nb + j] = (i == j) ? 1.0 : 0.0; }
    double scale = 0; for (int i = 0; i < nb * nb; ++i) scale = std::max(scale, std::fabs(a[i]));
    if (!(scale > 0) || !std::isfinite(scale)) return false;
    for (int k = 0; k < nb; ++k) {
        int pr = k;
        for (int i = k + 1; i < nb; ++i) if (std::fabs(m[i][k]) > std::fabs(m[pr][k])) pr = i;
        if (std::fabs(m[pr][k]) <= 1e-14 * scale) return false;
        if (pr != k) for (int j = 0; j < 2 * nb; ++j) std::swap(m[k][j], m[pr][j]);
        const double inv = 1.0 / m[k][k];
        for (int j = 0; j < 2 * nb; ++j) m[k][j] *= inv;
        for (int i = 0; i < nb; ++i) {
            if (i == k || m[i][k] == 0.0) continue;
            const double f = m[i][k];
            for (int j = 0; j < 2 * nb; ++j) m[i][j] -= f * m[k][j];
        }
    }
    for (int i = 0; i < nb; ++i) for (int j = 0; j < nb; ++j) a[i * nb + j] = m[i][nb + j];
    return true;
}

double blockConditionEstimate(const double* a, int nb) {
    double inv[kMaxBlock * kMaxBlock];
    std::memcpy(inv, a, sizeof(double) * nb * nb);
    if (!invertBlock(inv, nb)) return std::numeric_limits<double>::infinity();
    double n1 = 0, ni = 0;
    for (int j = 0; j < nb; ++j) {
        double c = 0, d = 0;
        for (int i = 0; i < nb; ++i) { c += std::fabs(a[i * nb + j]); d += std::fabs(inv[i * nb + j]); }
        n1 = std::max(n1, c); ni = std::max(ni, d);
    }
    return n1 * ni;
}

// ====================================================================== BSR
void BSR::finalizeDiag() {
    diag.assign(n, -1);
    for (int i = 0; i < n; ++i)
        for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p) if (col[p] == i) { diag[i] = p; break; }
}

void BSR::spmv(const double* x, double* y) const {
    const int b = nb, b2 = nb * nb;
    auto row = [&](int i) {
        double* yi = y + size_t(i) * b;
        for (int e = 0; e < b; ++e) yi[e] = 0;
        for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p) {
            const double* a = val.data() + size_t(p) * b2;
            const double* xj = x + size_t(col[p]) * b;
            RESSIM_BLOCK_DISPATCH(b, gemvAdd<2>(a, xj, yi, 2), gemvAdd<3>(a, xj, yi, 3), gemvAddD(a, xj, yi, b))
        }
    };
    if (useParallel(n > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) row(i);
    } else {
        for (int i = 0; i < n; ++i) row(i);
    }
}

CSR BSR::toCSR() const {
    CSR A; A.n = n * nb; A.rowPtr.assign(A.n + 1, 0);
    A.col.reserve(size_t(nnzb()) * nb * nb); A.val.reserve(size_t(nnzb()) * nb * nb);
    for (int i = 0; i < n; ++i)
        for (int e = 0; e < nb; ++e) {
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p)
                for (int v = 0; v < nb; ++v) { A.col.push_back(col[p] * nb + v); A.val.push_back(blk(p)[e * nb + v]); }
            A.rowPtr[i * nb + e + 1] = int(A.col.size());
        }
    A.finalizeDiag();
    return A;
}

BSR BSR::fromCSR(const CSR& A, int nBlocks, int b) {
    BSR S; S.n = nBlocks; S.nb = b; S.rowPtr.assign(nBlocks + 1, 0);
    const int nc = nBlocks * b;
    std::vector<int> mark(nBlocks, -1), cols;
    for (int i = 0; i < nBlocks; ++i) {
        cols.clear();
        for (int e = 0; e < b; ++e) {
            const int r = i * b + e;
            for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
                const int c = A.col[p];
                if (c >= nc) continue;
                const int j = c / b;
                if (mark[j] != i) { mark[j] = i; cols.push_back(j); }
            }
        }
        if (mark[i] != i) { mark[i] = i; cols.push_back(i); }
        std::sort(cols.begin(), cols.end());
        S.col.insert(S.col.end(), cols.begin(), cols.end());
        S.rowPtr[i + 1] = int(S.col.size());
    }
    S.val.assign(size_t(S.col.size()) * b * b, 0.0);
    std::vector<int> pos(nBlocks, -1);
    for (int i = 0; i < nBlocks; ++i) {
        for (int p = S.rowPtr[i]; p < S.rowPtr[i + 1]; ++p) pos[S.col[p]] = p;
        for (int e = 0; e < b; ++e) {
            const int r = i * b + e;
            for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
                const int c = A.col[p];
                if (c >= nc) continue;
                S.blk(pos[c / b])[e * b + c % b] += A.val[p];
            }
        }
    }
    S.finalizeDiag();
    return S;
}

// ====================================================================== BILU0
namespace {
// block ILU(0) factorisation in place: strictly-lower blocks hold L (unit diagonal implied), the diagonal blocks hold
// U_ii^-1, upper blocks U_ij (shared by BILU0 and the float32 BILU0f)
// one BILU0 row (rows k < i of the L pattern must be final); pos: per-thread scratch of size n filled with -1
inline int bilu0Row(BSR& LU, int i, std::vector<int>& pos) {
    const int b = LU.nb, b2 = b * b;
    int nReg = 0;
    double tmp[kMaxBlock * kMaxBlock];
    {
        const int r0 = LU.rowPtr[i], r1 = LU.rowPtr[i + 1];
        for (int p = r0; p < r1; ++p) pos[LU.col[p]] = p;
        for (int p = r0; p < r1; ++p) {
            const int k = LU.col[p];
            if (k >= i) break;
            // L_ik = A_ik * U_kk^-1 (the diagonal block of row k already holds its inverse)
            gemm(LU.blk(p), LU.blk(LU.diag[k]), tmp, b);
            std::memcpy(LU.blk(p), tmp, sizeof(double) * b2);
            for (int q = LU.diag[k] + 1; q < LU.rowPtr[k + 1]; ++q) {
                const int j = LU.col[q];
                if (pos[j] >= r0 && pos[j] < r1 && LU.col[pos[j]] == j) gemmSub(tmp, LU.blk(q), LU.blk(pos[j]), b);
            }
        }
        for (int p = r0; p < r1; ++p) pos[LU.col[p]] = -1;
        double* d = LU.blk(LU.diag[i]);
        std::memcpy(tmp, d, sizeof(double) * b2);
        if (!invertBlock(d, b)) {
            // singular pivot block (e.g. an isolated cell): shift its diagonal and retry, else use identity
            ++nReg;
            double sc = 0; for (int e = 0; e < b2; ++e) sc = std::max(sc, std::fabs(tmp[e]));
            if (!(sc > 0) || !std::isfinite(sc)) sc = 1.0;
            std::memcpy(d, tmp, sizeof(double) * b2);
            for (int e = 0; e < b; ++e) d[e * b + e] += (d[e * b + e] >= 0 ? 1e-8 : -1e-8) * sc + (d[e * b + e] == 0 ? sc * 1e-8 : 0);
            if (!invertBlock(d, b)) { for (int e = 0; e < b2; ++e) d[e] = 0; for (int e = 0; e < b; ++e) d[e * b + e] = 1.0 / sc; }
        }
    }
    return nReg;
}
// C45 lean BILU0 row: the factors are stored in float (F, pattern of 'pat'); row i is loaded into a double scratch row,
// eliminated against the (rounded) factors of the earlier rows exactly as bilu0Row does, its pivot block inverted in
// double, and stored back rounded
inline int bilu0RowF(const BSR& pat, float* F, int i, std::vector<int>& pos, std::vector<double>& buf) {
    const int b = pat.nb, b2 = b * b;
    int nReg = 0;
    double tmp[kMaxBlock * kMaxBlock], u[kMaxBlock * kMaxBlock];
    const int r0 = pat.rowPtr[i], r1 = pat.rowPtr[i + 1];
    buf.resize(size_t(r1 - r0) * b2);
    for (int p = r0; p < r1; ++p) {
        pos[pat.col[p]] = p;
        const float* f = F + size_t(p) * b2; double* d = buf.data() + size_t(p - r0) * b2;
        for (int e = 0; e < b2; ++e) d[e] = double(f[e]);
    }
    for (int p = r0; p < r1; ++p) {
        const int k = pat.col[p];
        if (k >= i) break;
        const float* fk = F + size_t(pat.diag[k]) * b2;
        for (int e = 0; e < b2; ++e) u[e] = double(fk[e]);
        double* lp = buf.data() + size_t(p - r0) * b2;
        gemm(lp, u, tmp, b);
        std::memcpy(lp, tmp, sizeof(double) * b2);
        for (int q = pat.diag[k] + 1; q < pat.rowPtr[k + 1]; ++q) {
            const int j = pat.col[q];
            if (pos[j] >= r0 && pos[j] < r1 && pat.col[pos[j]] == j) {
                const float* fq = F + size_t(q) * b2;
                for (int e = 0; e < b2; ++e) u[e] = double(fq[e]);
                gemmSub(tmp, u, buf.data() + size_t(pos[j] - r0) * b2, b);
            }
        }
    }
    for (int p = r0; p < r1; ++p) pos[pat.col[p]] = -1;
    double* d = buf.data() + size_t(pat.diag[i] - r0) * b2;
    std::memcpy(tmp, d, sizeof(double) * b2);
    if (!invertBlock(d, b)) {
        ++nReg;
        double sc = 0; for (int e = 0; e < b2; ++e) sc = std::max(sc, std::fabs(tmp[e]));
        if (!(sc > 0) || !std::isfinite(sc)) sc = 1.0;
        std::memcpy(d, tmp, sizeof(double) * b2);
        for (int e = 0; e < b; ++e) d[e * b + e] += (d[e * b + e] >= 0 ? 1e-8 : -1e-8) * sc + (d[e * b + e] == 0 ? sc * 1e-8 : 0);
        if (!invertBlock(d, b)) { for (int e = 0; e < b2; ++e) d[e] = 0; for (int e = 0; e < b; ++e) d[e * b + e] = 1.0 / sc; }
    }
    for (int p = r0; p < r1; ++p) {
        float* f = F + size_t(p) * b2; const double* x = buf.data() + size_t(p - r0) * b2;
        for (int e = 0; e < b2; ++e) f[e] = float(x[e]);
    }
    return nReg;
}
void bilu0InPlace(BSR& LU, int& nReg) {
    std::vector<int> pos(size_t(LU.n), -1);
    nReg = 0;
    for (int i = 0; i < LU.n; ++i) nReg += bilu0Row(LU, i, pos);
}
// level sets of a lower (forward) or upper (backward) dependency graph -> rows ordered by level, level starts
void levelSets(const BSR& A, bool forward, std::vector<int>& ord, std::vector<int>& start) {
    const int n = A.n;
    std::vector<int> lev(size_t(n), 0);
    int nl = 0;
    if (forward) {
        for (int i = 0; i < n; ++i) {
            int l = 0;
            for (int p = A.rowPtr[i]; p < A.diag[i]; ++p) l = std::max(l, lev[size_t(A.col[p])] + 1);
            lev[size_t(i)] = l; nl = std::max(nl, l + 1);
        }
    } else {
        for (int i = n - 1; i >= 0; --i) {
            int l = 0;
            for (int p = A.diag[i] + 1; p < A.rowPtr[i + 1]; ++p) l = std::max(l, lev[size_t(A.col[p])] + 1);
            lev[size_t(i)] = l; nl = std::max(nl, l + 1);
        }
    }
    start.assign(size_t(nl) + 1, 0);
    for (int i = 0; i < n; ++i) start[size_t(lev[size_t(i)]) + 1]++;
    for (int l = 0; l < nl; ++l) start[size_t(l) + 1] += start[size_t(l)];
    ord.assign(size_t(n), 0);
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < n; ++i) ord[size_t(fill[size_t(lev[size_t(i)])]++)] = i;   // ascending within a level
}
} // namespace

void BILU0::buildLevels() {
    levelSets(LU_, true, fOrd_, fStart_);
    levelSets(LU_, false, bOrd_, bStart_);
    // worth it only with enough rows per level to amortise one barrier per level
    const int nlf = int(fStart_.size()) - 1, nlb = int(bStart_.size()) - 1;
    // C21 H3: measured SLOWER than the serial sweep on 2 cores (54k cells: solve 22.2 s vs 19.3 s) - rows of one level are
    // scattered in memory. Opt-in until the level-permuted (contiguous) storage is done; see research/CYCLE3_REPORT.md
    const char* on = std::getenv("RESSIM_LEVEL_ILU");
    levelPar_ = (on && *on == '1') && LU_.n > 4000 && nlf > 0 && nlb > 0 && LU_.n / std::max(nlf, nlb) >= 96;
    // C28 contiguous variant: forced with RESSIM_LEVEL_ILU=2, off with 0/1; default ON for big matrices with enough rows per
    // level and more than one thread (-c28ilu, 54k cells / 116 rows per level, 2 threads: apply 2.09 ms vs 2.18 ms serial
    // vs 2.39 ms for the scattered C21 level scheduling; identical iterations). The gain grows with the core count and
    // the rows per level (3-D grids: nx*ny*nz / (nx+ny+nz) rows per level)
    contig_ = false;
    const bool force = on && *on == '2', off = on && (*on == '0' || *on == '1');
    const bool want = force || (!off && LU_.n >= 20000 && nlf > 0 && LU_.n / nlf >= 128 && useParallel(true));
    if (want && LU_.n > 1000 && nlf > 0 && LU_.n / nlf >= 32) buildContiguous();
}

void BILU0::buildContiguous() {
    const int n = LU_.n, b = LU_.nb;
    const size_t b2 = size_t(b) * b;
    // structural symmetry is required for the reversed forward levels to be valid backward levels
    for (int i = 0; i < n; ++i)
        for (int p = LU_.rowPtr[i]; p < LU_.rowPtr[i + 1]; ++p) {
            const int j = LU_.col[p];
            const int* c0 = LU_.col.data() + LU_.rowPtr[j], *c1 = LU_.col.data() + LU_.rowPtr[j + 1];
            if (!std::binary_search(c0, c1, i)) { contig_ = false; return; }
        }
    perm_ = fOrd_;
    std::vector<int> iperm(static_cast<size_t>(n));
    for (int k = 0; k < n; ++k) iperm[size_t(perm_[size_t(k)])] = k;
    P_ = BSR(); P_.n = n; P_.nb = b; P_.rowPtr.assign(size_t(n) + 1, 0);
    P_.col.reserve(LU_.col.size()); valMap_.clear(); valMap_.reserve(LU_.col.size());
    std::vector<std::pair<int, int>> row;
    for (int k = 0; k < n; ++k) {
        const int i = perm_[size_t(k)];
        row.clear();
        for (int p = LU_.rowPtr[i]; p < LU_.rowPtr[i + 1]; ++p) row.push_back({iperm[size_t(LU_.col[p])], p});
        std::sort(row.begin(), row.end());
        for (const auto& e : row) { P_.col.push_back(e.first); valMap_.push_back(e.second); }
        P_.rowPtr[size_t(k) + 1] = int(P_.col.size());
    }
    if (!lean_) P_.val.assign(P_.col.size() * b2, 0.0); else std::vector<double>().swap(P_.val);
    P_.finalizeDiag();
    // in the new numbering level l = rows [fStart_[l], fStart_[l+1]); every L entry of a row points to an earlier level
    for (int k = 0; k < n; ++k)
        for (int p = P_.rowPtr[k]; p < P_.diag[k]; ++p) if (P_.col[p] >= fStart_[size_t(std::upper_bound(fStart_.begin(), fStart_.end(), k) - fStart_.begin() - 1)]) { contig_ = false; return; }
    contig_ = true;
}

void BILU0::factor(const BSR& A) {
    if (LU_.n == A.n && LU_.nb == A.nb && LU_.col == A.col && LU_.rowPtr == A.rowPtr && LU_.diag.size() == size_t(A.n)) {
        if (!contig_ && !lean_) {                       // same pattern: values only (no reallocation)
            if (LU_.val.size() != A.val.size()) LU_.val = A.val; else std::copy(A.val.begin(), A.val.end(), LU_.val.begin());
        }
    } else if (lean_) {                                 // C45: pattern only - the values go straight into float storage
        LU_ = BSR{}; LU_.n = A.n; LU_.nb = A.nb; LU_.rowPtr = A.rowPtr; LU_.col = A.col; LU_.diag = A.diag;
        fOrd_.clear(); contig_ = false;
    } else { LU_ = A; fOrd_.clear(); contig_ = false; }
    if (LU_.diag.size() != size_t(LU_.n)) LU_.finalizeDiag();
    if (fOrd_.size() != size_t(LU_.n)) buildLevels();
    if (lean_) {
        if (!LU_.val.empty()) std::vector<double>().swap(LU_.val);
        if (!P_.val.empty()) std::vector<double>().swap(P_.val);
        const BSR& pat = contig_ ? P_ : LU_;
        const size_t b2 = size_t(pat.nb) * pat.nb, ne = pat.col.size();
        std::vector<float>& F = contig_ ? Pf_ : Lf_;
        if (contig_) std::vector<float>().swap(Lf_); else std::vector<float>().swap(Pf_);
        F.resize(ne * b2);
        #pragma omp parallel for schedule(static) if (useParallel(ne > 20000))
        for (long long e = 0; e < (long long)ne; ++e) {
            const double* src = &A.val[size_t(contig_ ? valMap_[size_t(e)] : int(e)) * b2];
            float* dst = &F[size_t(e) * b2];
            for (size_t q = 0; q < b2; ++q) dst[q] = float(src[q]);
        }
        int reg = 0;
        if (contig_ && useParallel(true)) {
            const int nl = int(fStart_.size()) - 1;
            #pragma omp parallel reduction(+ : reg)
            {
                std::vector<int> pos(size_t(pat.n), -1); std::vector<double> buf;
                for (int l = 0; l < nl; ++l) {
                    #pragma omp for schedule(static)
                    for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) reg += bilu0RowF(pat, F.data(), k, pos, buf);
                }
            }
        } else {
            std::vector<int> pos(size_t(pat.n), -1); std::vector<double> buf;
            for (int k = 0; k < pat.n; ++k) reg += bilu0RowF(pat, F.data(), k, pos, buf);
        }
        nReg_ = reg;
        return;
    }
    if (!Lf_.empty()) std::vector<float>().swap(Lf_);
    if (contig_) {   // C28: level-parallel factorisation in the level-ordered (contiguous) numbering
        // C37 memory: the factors live in P_ only - LU_ keeps the pattern (n, rowPtr, col, diag) for the reuse test
        if (!LU_.val.empty()) std::vector<double>().swap(LU_.val);
        const size_t b2 = size_t(P_.nb) * P_.nb;
        const int ne = int(valMap_.size());
        if (P_.val.size() != size_t(ne) * b2) P_.val.assign(size_t(ne) * b2, 0.0);   // C38: freed after a float factorisation
        if (!Pf_.empty()) std::vector<float>().swap(Pf_);
        #pragma omp parallel for schedule(static) if (useParallel(ne > 20000))
        for (int e = 0; e < ne; ++e) std::memcpy(&P_.val[size_t(e) * b2], &A.val[size_t(valMap_[size_t(e)]) * b2], sizeof(double) * b2);
        int reg = 0;
        const int nl = int(fStart_.size()) - 1;
        if (useParallel(true)) {
            #pragma omp parallel reduction(+ : reg)
            {
                std::vector<int> pos(size_t(P_.n), -1);
                for (int l = 0; l < nl; ++l) {
                    #pragma omp for schedule(static)
                    for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) reg += bilu0Row(P_, k, pos);
                }
            }
        } else { std::vector<int> pos(size_t(P_.n), -1); for (int k = 0; k < P_.n; ++k) reg += bilu0Row(P_, k, pos); }
        nReg_ = reg;
        if (floatStore_) {   // C38 float32 factors for the sweeps
            Pf_.resize(P_.val.size());
            const long long nv = (long long)P_.val.size();
            #pragma omp parallel for schedule(static) if (useParallel(nv > 200000))
            for (long long i = 0; i < nv; ++i) Pf_[size_t(i)] = float(P_.val[size_t(i)]);
            std::vector<double>().swap(P_.val);
        }
        return;
    }
    if (levelPar_ && useParallel(true)) {
        int reg = 0;
        const int nl = int(fStart_.size()) - 1;
        #pragma omp parallel reduction(+ : reg)
        {
            std::vector<int> pos(size_t(LU_.n), -1);
            for (int l = 0; l < nl; ++l) {
                #pragma omp for schedule(static)
                for (int q = fStart_[size_t(l)]; q < fStart_[size_t(l) + 1]; ++q) reg += bilu0Row(LU_, fOrd_[size_t(q)], pos);
            }
        }
        nReg_ = reg;
    } else bilu0InPlace(LU_, nReg_);
}

void BILU0::apply(const double* r, double* z) const {
    const int n = LU_.n, b = LU_.nb;
    if (contig_) {   // C28: gather -> level sweeps over contiguous rows -> scatter
        const int nl = int(fStart_.size()) - 1;
        rp_.resize(size_t(n) * b); zp_.resize(size_t(n) * b);
        // C33: one barrier per level - only worth it with many rows per level (thin corner-point models have hundreds
        // of levels with a few dozen rows; with a busy second core every barrier then costs milliseconds: OPM model2
        // ran 90 % of its wall time in gomp barriers)
        const bool par = useParallel(n > 20000 && n / std::max(1, nl) >= 128);
        if (floatFactors()) {   // C38: same sweeps with float32 factors (double vectors)
            const size_t b2 = size_t(b) * b;
            const float* F = Pf_.data();
            #pragma omp parallel if (par)
            {
                #pragma omp for schedule(static)
                for (int k = 0; k < n; ++k) for (int e = 0; e < b; ++e) rp_[size_t(k) * b + e] = r[size_t(perm_[size_t(k)]) * b + e];
                for (int l = 0; l < nl; ++l) {
                    #pragma omp for schedule(static)
                    for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) {
                        double* zi = zp_.data() + size_t(k) * b;
                        for (int e = 0; e < b; ++e) zi[e] = rp_[size_t(k) * b + e];
                        for (int p = P_.rowPtr[k]; p < P_.diag[k]; ++p) {
                            const float* a = F + size_t(p) * b2; const double* zk = zp_.data() + size_t(P_.col[p]) * b;
                            for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(a[e * b + f]) * zk[f]; zi[e] -= s; }
                        }
                    }
                }
                for (int l = nl - 1; l >= 0; --l) {
                    #pragma omp for schedule(static)
                    for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) {
                        double t[kMaxBlock];
                        double* zi = zp_.data() + size_t(k) * b;
                        for (int p = P_.diag[k] + 1; p < P_.rowPtr[k + 1]; ++p) {
                            const float* a = F + size_t(p) * b2; const double* zj = zp_.data() + size_t(P_.col[p]) * b;
                            for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(a[e * b + f]) * zj[f]; zi[e] -= s; }
                        }
                        const float* d = F + size_t(P_.diag[k]) * b2;
                        for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(d[e * b + f]) * zi[f]; t[e] = s; }
                        for (int e = 0; e < b; ++e) zi[e] = t[e];
                    }
                }
                #pragma omp for schedule(static)
                for (int k = 0; k < n; ++k) for (int e = 0; e < b; ++e) z[size_t(perm_[size_t(k)]) * b + e] = zp_[size_t(k) * b + e];
            }
            return;
        }
        #pragma omp parallel if (par)
        {
            #pragma omp for schedule(static)
            for (int k = 0; k < n; ++k) for (int e = 0; e < b; ++e) rp_[size_t(k) * b + e] = r[size_t(perm_[size_t(k)]) * b + e];
            for (int l = 0; l < nl; ++l) {
                #pragma omp for schedule(static)
                for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) {
                    double* zi = zp_.data() + size_t(k) * b;
                    for (int e = 0; e < b; ++e) zi[e] = rp_[size_t(k) * b + e];
                    for (int p = P_.rowPtr[k]; p < P_.diag[k]; ++p) {
                        const double* a = P_.blk(p); const double* zk = zp_.data() + size_t(P_.col[p]) * b;
                        RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zk, zi, 2), gemvSub<3>(a, zk, zi, 3), gemvSubD(a, zk, zi, b))
                    }
                }
            }
            for (int l = nl - 1; l >= 0; --l) {
                #pragma omp for schedule(static)
                for (int k = fStart_[size_t(l)]; k < fStart_[size_t(l) + 1]; ++k) {
                    double t[kMaxBlock];
                    double* zi = zp_.data() + size_t(k) * b;
                    for (int p = P_.diag[k] + 1; p < P_.rowPtr[k + 1]; ++p) {
                        const double* a = P_.blk(p); const double* zj = zp_.data() + size_t(P_.col[p]) * b;
                        RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zj, zi, 2), gemvSub<3>(a, zj, zi, 3), gemvSubD(a, zj, zi, b))
                    }
                    for (int e = 0; e < b; ++e) t[e] = 0;
                    gemvAddN(P_.blk(P_.diag[k]), zi, t, b);
                    for (int e = 0; e < b; ++e) zi[e] = t[e];
                }
            }
            #pragma omp for schedule(static)
            for (int k = 0; k < n; ++k) for (int e = 0; e < b; ++e) z[size_t(perm_[size_t(k)]) * b + e] = zp_[size_t(k) * b + e];
        }
        return;
    }
    if (!Lf_.empty()) {   // C45 lean float factors, original numbering (serial sweeps, double vectors)
        const size_t b2 = size_t(b) * b;
        const float* F = Lf_.data();
        for (int i = 0; i < n; ++i) {
            double* zi = z + size_t(i) * b;
            for (int e = 0; e < b; ++e) zi[e] = r[size_t(i) * b + e];
            for (int p = LU_.rowPtr[i]; p < LU_.diag[i]; ++p) {
                const float* a = F + size_t(p) * b2; const double* zk = z + size_t(LU_.col[p]) * b;
                for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(a[e * b + f]) * zk[f]; zi[e] -= s; }
            }
        }
        double t[kMaxBlock];
        for (int i = n - 1; i >= 0; --i) {
            double* zi = z + size_t(i) * b;
            for (int p = LU_.diag[i] + 1; p < LU_.rowPtr[i + 1]; ++p) {
                const float* a = F + size_t(p) * b2; const double* zj = z + size_t(LU_.col[p]) * b;
                for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(a[e * b + f]) * zj[f]; zi[e] -= s; }
            }
            const float* d = F + size_t(LU_.diag[i]) * b2;
            for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += double(d[e * b + f]) * zi[f]; t[e] = s; }
            for (int e = 0; e < b; ++e) zi[e] = t[e];
        }
        return;
    }
    if (levelPar_ && useParallel(n > 20000 && n / std::max(1, int(fStart_.size()) - 1) >= 128)) {   // C33 (see above)
        const int nlf = int(fStart_.size()) - 1, nlb = int(bStart_.size()) - 1;
        #pragma omp parallel
        {
            for (int l = 0; l < nlf; ++l) {
                #pragma omp for schedule(static)
                for (int q = fStart_[size_t(l)]; q < fStart_[size_t(l) + 1]; ++q) {
                    const int i = fOrd_[size_t(q)];
                    double* zi = z + size_t(i) * b;
                    for (int e = 0; e < b; ++e) zi[e] = r[size_t(i) * b + e];
                    for (int p = LU_.rowPtr[i]; p < LU_.diag[i]; ++p) {
                        const double* a = LU_.blk(p); const double* zk = z + size_t(LU_.col[p]) * b;
                        RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zk, zi, 2), gemvSub<3>(a, zk, zi, 3), gemvSubD(a, zk, zi, b))
                    }
                }
            }
            for (int l = 0; l < nlb; ++l) {
                #pragma omp for schedule(static)
                for (int q = bStart_[size_t(l)]; q < bStart_[size_t(l) + 1]; ++q) {
                    const int i = bOrd_[size_t(q)];
                    double t[kMaxBlock];
                    double* zi = z + size_t(i) * b;
                    for (int p = LU_.diag[i] + 1; p < LU_.rowPtr[i + 1]; ++p) {
                        const double* a = LU_.blk(p); const double* zj = z + size_t(LU_.col[p]) * b;
                        RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zj, zi, 2), gemvSub<3>(a, zj, zi, 3), gemvSubD(a, zj, zi, b))
                    }
                    for (int e = 0; e < b; ++e) t[e] = 0;
                    gemvAddN(LU_.blk(LU_.diag[i]), zi, t, b);
                    for (int e = 0; e < b; ++e) zi[e] = t[e];
                }
            }
        }
        return;
    }
    // forward: y_i = r_i - sum_{k<i} L_ik y_k   (y stored in z)
    for (int i = 0; i < n; ++i) {
        double* zi = z + size_t(i) * b;
        for (int e = 0; e < b; ++e) zi[e] = r[size_t(i) * b + e];
        for (int p = LU_.rowPtr[i]; p < LU_.diag[i]; ++p) {
            const double* a = LU_.blk(p); const double* zk = z + size_t(LU_.col[p]) * b;
            RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zk, zi, 2), gemvSub<3>(a, zk, zi, 3), gemvSubD(a, zk, zi, b))
        }
    }
    // backward: z_i = U_ii^-1 (y_i - sum_{j>i} U_ij z_j)
    double t[kMaxBlock];
    for (int i = n - 1; i >= 0; --i) {
        double* zi = z + size_t(i) * b;
        for (int p = LU_.diag[i] + 1; p < LU_.rowPtr[i + 1]; ++p) {
            const double* a = LU_.blk(p); const double* zj = z + size_t(LU_.col[p]) * b;
            RESSIM_BLOCK_DISPATCH(b, gemvSub<2>(a, zj, zi, 2), gemvSub<3>(a, zj, zi, 3), gemvSubD(a, zj, zi, b))
        }
        for (int e = 0; e < b; ++e) t[e] = 0;
        gemvAddN(LU_.blk(LU_.diag[i]), zi, t, b);
        for (int e = 0; e < b; ++e) zi[e] = t[e];
    }
}

// ====================================================================== C18 multicolour BILU0
Colouring colourBlocks(const BSR& A) {
    Colouring C;
    const int n = A.n;
    C.colour.assign(size_t(n), -1);
    // symmetric adjacency (the pattern of a reservoir Jacobian is structurally symmetric; be safe for one-sided entries)
    std::vector<std::vector<int>> adj(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        for (int p = A.rowPtr[size_t(i)]; p < A.rowPtr[size_t(i) + 1]; ++p) {
            const int j = A.col[size_t(p)];
            if (j != i) { adj[size_t(i)].push_back(j); adj[size_t(j)].push_back(i); }
        }
    for (auto& v : adj) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    // largest-degree-first greedy (fewest colours on irregular graphs: well couplings, NNCs, LGR)
    std::vector<int> order(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) order[size_t(i)] = i;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return adj[size_t(a)].size() > adj[size_t(b)].size(); });
    std::vector<int> mark;
    for (int i : order) {
        mark.assign(adj[size_t(i)].size() + 2, 0);
        for (int j : adj[size_t(i)]) { const int cj = C.colour[size_t(j)]; if (cj >= 0 && cj < int(mark.size())) mark[size_t(cj)] = 1; }
        int c = 0; while (mark[size_t(c)]) ++c;
        C.colour[size_t(i)] = c;
        C.colours = std::max(C.colours, c + 1);
    }
    C.start.assign(size_t(C.colours) + 1, 0);
    for (int i = 0; i < n; ++i) C.start[size_t(C.colour[size_t(i)]) + 1]++;
    for (int c = 0; c < C.colours; ++c) C.start[size_t(c) + 1] += C.start[size_t(c)];
    C.perm.assign(size_t(n), 0); C.iperm.assign(size_t(n), 0);
    std::vector<int> pos(C.start.begin(), C.start.end() - 1);
    for (int i = 0; i < n; ++i) { const int k = pos[size_t(C.colour[size_t(i)])]++; C.perm[size_t(k)] = i; C.iperm[size_t(i)] = k; }
    return C;
}

BSR permuteBsr(const BSR& A, const Colouring& C) {
    BSR P; P.n = A.n; P.nb = A.nb;
    const int n = A.n, b2 = A.nb * A.nb;
    P.rowPtr.assign(size_t(n) + 1, 0);
    for (int k = 0; k < n; ++k) { const int i = C.perm[size_t(k)]; P.rowPtr[size_t(k) + 1] = P.rowPtr[size_t(k)] + (A.rowPtr[size_t(i) + 1] - A.rowPtr[size_t(i)]); }
    P.col.resize(size_t(A.nnzb())); P.val.resize(A.val.size());
    std::vector<std::pair<int, int>> row;
    for (int k = 0; k < n; ++k) {
        const int i = C.perm[size_t(k)];
        row.clear();
        for (int p = A.rowPtr[size_t(i)]; p < A.rowPtr[size_t(i) + 1]; ++p) row.push_back({C.iperm[size_t(A.col[size_t(p)])], p});
        std::sort(row.begin(), row.end());
        int d = P.rowPtr[size_t(k)];
        for (const auto& e : row) {
            P.col[size_t(d)] = e.first;
            std::memcpy(P.val.data() + size_t(d) * b2, A.val.data() + size_t(e.second) * b2, sizeof(double) * size_t(b2));
            ++d;
        }
    }
    P.finalizeDiag();
    return P;
}

namespace {
// block ILU(0) of a colour-permuted matrix: rows of one colour are independent (their lower entries only reference earlier
// colours, whose rows are final) -> parallel over the rows of a colour
void bilu0Coloured(BSR& LU, const Colouring& C, int& nReg) {
    const int b = LU.nb, b2 = b * b;
    int reg = 0;
    for (int c = 0; c < C.colours; ++c) {
        const int r0c = C.start[size_t(c)], r1c = C.start[size_t(c) + 1];
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : reg) if (r1c - r0c > 2000)
        for (int i = r0c; i < r1c; ++i) {
            double tmp[kMaxBlock * kMaxBlock];
            const int r0 = LU.rowPtr[size_t(i)], r1 = LU.rowPtr[size_t(i) + 1];
            for (int p = r0; p < r1; ++p) {
                const int k = LU.col[size_t(p)];
                if (k >= i) break;
                gemm(LU.blk(p), LU.blk(LU.diag[size_t(k)]), tmp, b);
                std::memcpy(LU.blk(p), tmp, sizeof(double) * b2);
                // row i -= L_ik U_kj for j > k present in row i (both column lists sorted: merge)
                int q = LU.diag[size_t(k)] + 1, s = p + 1;
                const int qe = LU.rowPtr[size_t(k) + 1];
                while (q < qe && s < r1) {
                    const int jq = LU.col[size_t(q)], js = LU.col[size_t(s)];
                    if (jq == js) { gemmSub(tmp, LU.blk(q), LU.blk(s), b); ++q; ++s; }
                    else if (jq < js) ++q; else ++s;
                }
            }
            double* d = LU.blk(LU.diag[size_t(i)]);
            std::memcpy(tmp, d, sizeof(double) * b2);
            if (!invertBlock(d, b)) {
                ++reg;
                double sc = 0; for (int e = 0; e < b2; ++e) sc = std::max(sc, std::fabs(tmp[e]));
                if (!(sc > 0) || !std::isfinite(sc)) sc = 1.0;
                std::memcpy(d, tmp, sizeof(double) * b2);
                for (int e = 0; e < b; ++e) d[e * b + e] += (d[e * b + e] >= 0 ? 1e-8 : -1e-8) * sc + (d[e * b + e] == 0 ? sc * 1e-8 : 0);
                if (!invertBlock(d, b)) { for (int e = 0; e < b2; ++e) d[e] = 0; for (int e = 0; e < b; ++e) d[e * b + e] = 1.0 / sc; }
            }
        }
    }
    nReg = reg;
}
} // namespace

void MCBILU0::factor(const BSR& A) {
    if (patRow_ != A.rowPtr || patCol_ != A.col || int(C_.colour.size()) != A.n) { C_ = colourBlocks(A); patRow_ = A.rowPtr; patCol_ = A.col; }
    LU_ = permuteBsr(A, C_);
    if (keepPermuted) PA_ = LU_; else PA_ = BSR();
    bilu0Coloured(LU_, C_, nReg_);
}

void MCBILU0::apply(const double* r, double* z) const {
    const int n = LU_.n, b = LU_.nb;
    rp_.resize(size_t(n) * b); zp_.resize(size_t(n) * b);
    for (int k = 0; k < n; ++k) std::memcpy(&rp_[size_t(k) * b], r + size_t(C_.perm[size_t(k)]) * b, sizeof(double) * b);
    double* y = zp_.data();
    for (int c = 0; c < C_.colours; ++c) {   // forward: y_i = r_i - sum_{k < i} L_ik y_k  (k in earlier colours)
        const int r0c = C_.start[size_t(c)], r1c = C_.start[size_t(c) + 1];
#pragma omp parallel for schedule(static) if (r1c - r0c > 4000)
        for (int i = r0c; i < r1c; ++i) {
            double* yi = y + size_t(i) * b;
            for (int e = 0; e < b; ++e) yi[e] = rp_[size_t(i) * b + e];
            for (int p = LU_.rowPtr[size_t(i)]; p < LU_.diag[size_t(i)]; ++p) gemvSubD(LU_.blk(p), y + size_t(LU_.col[size_t(p)]) * b, yi, b);
        }
    }
    for (int c = C_.colours - 1; c >= 0; --c) {   // backward: z_i = U_ii^-1 (y_i - sum_{j > i} U_ij z_j)
        const int r0c = C_.start[size_t(c)], r1c = C_.start[size_t(c) + 1];
#pragma omp parallel for schedule(static) if (r1c - r0c > 4000)
        for (int i = r0c; i < r1c; ++i) {
            double* yi = y + size_t(i) * b;
            for (int p = LU_.diag[size_t(i)] + 1; p < LU_.rowPtr[size_t(i) + 1]; ++p) gemvSubD(LU_.blk(p), y + size_t(LU_.col[size_t(p)]) * b, yi, b);
            double t[kMaxBlock] = {};
            gemvAddN(LU_.blk(LU_.diag[size_t(i)]), yi, t, b);
            for (int e = 0; e < b; ++e) yi[e] = t[e];
        }
    }
    for (int k = 0; k < n; ++k) std::memcpy(z + size_t(C_.perm[size_t(k)]) * b, &zp_[size_t(k) * b], sizeof(double) * b);
}

// ====================================================================== WellSchur
bool WellSchur::build(const CSR& A, int na, int nv, int nw) {
    na_ = na; nv_ = nv; nw_ = nw; fill_ = 0;
    const int nc = na * nv;
    // ---- well rows: diagonal and cell couplings
    dInv_.assign(nw, 0.0);
    wcPtr_.assign(nw + 1, 0); wcCol_.clear(); wcVal_.clear();
    for (int w = 0; w < nw; ++w) {
        const int r = nc + w;
        double d = 0, rowMax = 0;
        for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
            const int c = A.col[p];
            if (c == r) d = A.val[p];
            else if (c < nc) { wcCol_.push_back(c); wcVal_.push_back(A.val[p]); rowMax = std::max(rowMax, std::fabs(A.val[p])); }
            else if (A.val[p] != 0.0) return false;                    // well-well coupling (groups): not supported
        }
        if (!(std::fabs(d) > 1e-12 * std::max(rowMax, 1e-300)) || !std::isfinite(d)) return false;
        dInv_[w] = 1.0 / d;
        wcPtr_[w + 1] = int(wcCol_.size());
    }
    // ---- cell rows: well columns
    cwPtr_.assign(nc + 1, 0); cwWell_.clear(); cwVal_.clear();
    for (int r = 0; r < nc; ++r) {
        for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p)
            if (A.col[p] >= nc && A.val[p] != 0.0) { cwWell_.push_back(A.col[p] - nc); cwVal_.push_back(A.val[p]); }
        cwPtr_[r + 1] = int(cwWell_.size());
    }
    // C37: the block pattern depends on the scalar pattern and on the (value-dependent) well couplings only - reuse it
    // while they are unchanged (every Newton iteration of a step, usually the whole run): no pattern rebuild, no
    // reallocation of S_ (~50 MB per 100k cells for three phases)
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](const std::vector<int>& v) {
        h ^= (unsigned long long)v.size() + 0x9e3779b97f4a7c15ULL; h *= 1099511628211ULL;
        for (int x : v) { h ^= (unsigned long long)(unsigned)x; h *= 1099511628211ULL; }
    };
    mix(A.rowPtr); mix(A.col); mix(wcPtr_); mix(wcCol_); mix(cwPtr_); mix(cwWell_);
    h ^= (unsigned long long)na * 31 + (unsigned long long)nv * 7 + (unsigned long long)nw; h *= 1099511628211ULL;
    const bool reuse = h == patHash_ && S_.n == na && S_.nb == nv && !S_.col.empty();
    if (reuse) {
        fill_ = int(S_.col.size()) - base_;
        std::fill(S_.val.begin(), S_.val.end(), 0.0);
    } else {
    patHash_ = h;
    // ---- block pattern: cell pattern + couplings through wells
    std::vector<std::vector<int>> wellCells(nw);
    for (int w = 0; w < nw; ++w) {
        auto& v = wellCells[w];
        for (int p = wcPtr_[w]; p < wcPtr_[w + 1]; ++p) v.push_back(wcCol_[p] / nv);
        std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    S_ = BSR{}; S_.n = na; S_.nb = nv; S_.rowPtr.assign(na + 1, 0);
    std::vector<int> mark(na, -1), cols, wl;
    int base = 0;
    for (int i = 0; i < na; ++i) {
        cols.clear(); wl.clear();
        for (int e = 0; e < nv; ++e) {
            const int r = i * nv + e;
            for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
                const int c = A.col[p];
                if (c < nc) { const int j = c / nv; if (mark[j] != i) { mark[j] = i; cols.push_back(j); } }
            }
            for (int p = cwPtr_[r]; p < cwPtr_[r + 1]; ++p) wl.push_back(cwWell_[p]);
        }
        if (mark[i] != i) { mark[i] = i; cols.push_back(i); }
        base += int(cols.size());
        std::sort(wl.begin(), wl.end()); wl.erase(std::unique(wl.begin(), wl.end()), wl.end());
        for (int w : wl) for (int j : wellCells[w]) if (mark[j] != i) { mark[j] = i; cols.push_back(j); }
        std::sort(cols.begin(), cols.end());
        S_.col.insert(S_.col.end(), cols.begin(), cols.end());
        S_.rowPtr[i + 1] = int(S_.col.size());
    }
    fill_ = int(S_.col.size()) - base;
    base_ = base;
    S_.col.shrink_to_fit();
    S_.val.assign(size_t(S_.col.size()) * nv * nv, 0.0);
    }
    // ---- values
    auto fillRow = [&](int i, std::vector<int>& pos) {
        for (int p = S_.rowPtr[i]; p < S_.rowPtr[i + 1]; ++p) pos[S_.col[p]] = p;
        for (int e = 0; e < nv; ++e) {
            const int r = i * nv + e;
            for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
                const int c = A.col[p];
                if (c < nc) S_.blk(pos[c / nv])[e * nv + c % nv] += A.val[p];
            }
            for (int p = cwPtr_[r]; p < cwPtr_[r + 1]; ++p) {
                const int w = cwWell_[p];
                const double f = cwVal_[p] * dInv_[w];
                for (int q = wcPtr_[w]; q < wcPtr_[w + 1]; ++q) {
                    const int c = wcCol_[q];
                    S_.blk(pos[c / nv])[e * nv + c % nv] -= f * wcVal_[q];
                }
            }
        }
    };
    if (useParallel(na > kPar)) {
        #pragma omp parallel
        {
            std::vector<int> pos(na, -1);
            #pragma omp for schedule(dynamic, 256)
            for (int i = 0; i < na; ++i) fillRow(i, pos);
        }
    } else {
        std::vector<int> pos(na, -1);
        for (int i = 0; i < na; ++i) fillRow(i, pos);
    }
    if (!reuse) S_.finalizeDiag();
    return true;
}

// C45 block assembly ------------------------------------------------------------------------------------------------
void WellSchur::bindCells(int na, int nv, const std::vector<int>& start, const std::vector<int>& cells) {
    na_ = na; nv_ = nv; nw_ = 0; fill_ = 0;
    S_ = BSR{}; S_.n = na; S_.nb = nv;
    S_.rowPtr.assign(start.begin(), start.begin() + na + 1);
    S_.col.assign(cells.begin(), cells.begin() + start[size_t(na)]);
    S_.val.assign(S_.col.size() * size_t(nv) * size_t(nv), 0.0);
    S_.finalizeDiag();
    map_.resize(S_.col.size());
    for (size_t k = 0; k < map_.size(); ++k) map_[k] = int(k);
    base_ = int(S_.col.size());
    patHash_ = 0;                                         // the first buildInPlace() lays out the fill
}

bool WellSchur::buildInPlace(const CSR& A, int nw, const std::vector<int>& start, const std::vector<int>& cells) {
    const int na = na_, nv = nv_, nc = na * nv;
    nw_ = nw;
    // ---- well rows (as build(); nothing in S_ is touched before every check has passed)
    std::vector<double> dInv(size_t(nw), 0.0);
    wcPtr_.assign(size_t(nw) + 1, 0); wcCol_.clear(); wcVal_.clear();
    for (int w = 0; w < nw; ++w) {
        const int r = nc + w;
        double d = 0, rowMax = 0;
        for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p) {
            const int c = A.col[p];
            if (c == r) d = A.val[p];
            else if (c < nc) { wcCol_.push_back(c); wcVal_.push_back(A.val[p]); rowMax = std::max(rowMax, std::fabs(A.val[p])); }
            else if (A.val[p] != 0.0) return false;
        }
        if (!(std::fabs(d) > 1e-12 * std::max(rowMax, 1e-300)) || !std::isfinite(d)) return false;
        dInv[size_t(w)] = 1.0 / d;
        wcPtr_[size_t(w) + 1] = int(wcCol_.size());
    }
    dInv_.swap(dInv);
    // ---- cell rows: well columns (the cell rows of A hold nothing else)
    cwPtr_.assign(size_t(nc) + 1, 0); cwWell_.clear(); cwVal_.clear();
    for (int r = 0; r < nc; ++r) {
        for (int p = A.rowPtr[r]; p < A.rowPtr[r + 1]; ++p)
            if (A.col[p] >= nc && A.val[p] != 0.0) { cwWell_.push_back(A.col[p] - nc); cwVal_.push_back(A.val[p]); }
        cwPtr_[size_t(r) + 1] = int(cwWell_.size());
    }
    // ---- fill pattern: reused while the well couplings keep their structure (the cell pattern is fixed by bindCells)
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](const std::vector<int>& v) {
        h ^= (unsigned long long)v.size() + 0x9e3779b97f4a7c15ULL; h *= 1099511628211ULL;
        for (int x : v) { h ^= (unsigned long long)(unsigned)x; h *= 1099511628211ULL; }
    };
    mix(wcPtr_); mix(wcCol_); mix(cwPtr_); mix(cwWell_);
    h ^= (unsigned long long)na * 31 + (unsigned long long)nv * 7 + (unsigned long long)nw; h *= 1099511628211ULL;
    if (h == 0) h = 1;
    if (h != patHash_) {
        patHash_ = h;
        std::vector<std::vector<int>> wellCells(static_cast<size_t>(nw));
        for (int w = 0; w < nw; ++w) {
            auto& v = wellCells[size_t(w)];
            for (int p = wcPtr_[size_t(w)]; p < wcPtr_[size_t(w) + 1]; ++p) v.push_back(wcCol_[size_t(p)] / nv);
            std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
        }
        BSR T; T.n = na; T.nb = nv; T.rowPtr.assign(size_t(na) + 1, 0);
        std::vector<int> mark(size_t(na), -1), cols, wl;
        for (int i = 0; i < na; ++i) {
            cols.clear(); wl.clear();
            for (int k = start[size_t(i)]; k < start[size_t(i) + 1]; ++k) { const int j = cells[size_t(k)]; mark[size_t(j)] = i; cols.push_back(j); }
            for (int e = 0; e < nv; ++e) {
                const int r = i * nv + e;
                for (int p = cwPtr_[size_t(r)]; p < cwPtr_[size_t(r) + 1]; ++p) wl.push_back(cwWell_[size_t(p)]);
            }
            std::sort(wl.begin(), wl.end()); wl.erase(std::unique(wl.begin(), wl.end()), wl.end());
            for (int w : wl) for (int j : wellCells[size_t(w)]) if (mark[size_t(j)] != i) { mark[size_t(j)] = i; cols.push_back(j); }
            std::sort(cols.begin(), cols.end());
            T.col.insert(T.col.end(), cols.begin(), cols.end());
            T.rowPtr[size_t(i) + 1] = int(T.col.size());
        }
        T.col.shrink_to_fit();
        if (T.rowPtr == S_.rowPtr && T.col == S_.col) {
            // same layout as the current one (e.g. back to the previous couplings): keep values and map
        } else {
            T.val.assign(T.col.size() * size_t(nv) * size_t(nv), 0.0);
            const size_t b2 = size_t(nv) * size_t(nv);
            for (int i = 0; i < na; ++i) {
                int q = T.rowPtr[size_t(i)];
                for (int k = start[size_t(i)]; k < start[size_t(i) + 1]; ++k) {
                    const int j = cells[size_t(k)];
                    while (T.col[size_t(q)] != j) ++q;      // both sorted, the cell columns are a subset
                    std::memcpy(T.val.data() + size_t(q) * b2, S_.val.data() + size_t(map_[size_t(k)]) * b2, sizeof(double) * b2);
                    map_[size_t(k)] = q;
                }
            }
            T.finalizeDiag();
            S_ = std::move(T);
        }
        fill_ = int(S_.col.size()) - base_;
    }
    // ---- S = Acc - Acw Dww^-1 Awc in place (rows with well columns only; same order of operations as build())
    std::vector<int> rows;
    for (int i = 0; i < na; ++i) if (cwPtr_[size_t(i + 1) * nv] > cwPtr_[size_t(i) * nv]) rows.push_back(i);
    auto subRow = [&](int i, std::vector<int>& pos) {
        for (int p = S_.rowPtr[size_t(i)]; p < S_.rowPtr[size_t(i) + 1]; ++p) pos[size_t(S_.col[size_t(p)])] = p;
        for (int e = 0; e < nv; ++e) {
            const int r = i * nv + e;
            for (int p = cwPtr_[size_t(r)]; p < cwPtr_[size_t(r) + 1]; ++p) {
                const int w = cwWell_[size_t(p)];
                const double f = cwVal_[size_t(p)] * dInv_[size_t(w)];
                for (int q = wcPtr_[size_t(w)]; q < wcPtr_[size_t(w) + 1]; ++q) {
                    const int c = wcCol_[size_t(q)];
                    S_.blk(pos[size_t(c / nv)])[e * nv + c % nv] -= f * wcVal_[size_t(q)];
                }
            }
        }
    };
    const int nr = int(rows.size());
    if (useParallel(nr > 256)) {
        #pragma omp parallel
        {
            std::vector<int> pos(size_t(na), -1);
            #pragma omp for schedule(dynamic, 16)
            for (int t = 0; t < nr; ++t) subRow(rows[size_t(t)], pos);
        }
    } else {
        std::vector<int> pos(size_t(na), -1);
        for (int t = 0; t < nr; ++t) subRow(rows[size_t(t)], pos);
    }
    return true;
}

size_t WellSchur::bytes() const {
    return S_.val.capacity() * sizeof(double) + (S_.col.capacity() + S_.rowPtr.capacity() + S_.diag.capacity() + map_.capacity() +
           wcPtr_.capacity() + wcCol_.capacity() + cwPtr_.capacity() + cwWell_.capacity()) * sizeof(int) +
           (wcVal_.capacity() + cwVal_.capacity() + dInv_.capacity()) * sizeof(double);
}

void WellSchur::reduceRhs(const double* b, double* bc) const {
    const int nc = na_ * nv_;
    for (int r = 0; r < nc; ++r) {
        double s = b[r];
        for (int p = cwPtr_[r]; p < cwPtr_[r + 1]; ++p) s -= cwVal_[p] * dInv_[cwWell_[p]] * b[nc + cwWell_[p]];
        bc[r] = s;
    }
}

void WellSchur::recover(const double* b, const double* xc, double* x) const {
    const int nc = na_ * nv_;
    if (x != xc) std::copy(xc, xc + nc, x);
    for (int w = 0; w < nw_; ++w) {
        double s = b[nc + w];
        for (int p = wcPtr_[w]; p < wcPtr_[w + 1]; ++p) s -= wcVal_[p] * xc[wcCol_[p]];
        x[nc + w] = s * dInv_[w];
    }
}

// ====================================================================== BlockCPR
void BlockCPR::setup(const BSR& S, const double* accBlocks, bool reuseAmg, std::vector<double>* releaseAcc) {
    S_ = &S;
    const int n = S.n, b = S.nb, b2 = b * b;
    w_.assign(size_t(n) * b, 1.0);
    int nFb = 0;
    // C21: independent per cell -> parallel (was ~30 % of the CPR set-up)
    #pragma omp parallel for schedule(static) reduction(+ : nFb) if (useParallel(n > kPar))
    for (int i = 0; i < n; ++i) {
        int& nFallback_ = nFb;
        double D[kMaxBlock * kMaxBlock];
        std::memcpy(D, S.blk(S.diag[i]), sizeof(double) * b2);
        bool inverted = false;
        if (accBlocks) {
            // condition estimate from the explicit inverse, which is then used directly for the weights
            const double* B = accBlocks + size_t(i) * b2;
            double inv[kMaxBlock * kMaxBlock];
            std::memcpy(inv, B, sizeof(double) * b2);
            if (invertBlock(inv, b)) {
                double n1 = 0, ni = 0;
                for (int j = 0; j < b; ++j) {
                    double c = 0, d = 0;
                    for (int r = 0; r < b; ++r) { c += std::fabs(B[r * b + j]); d += std::fabs(inv[r * b + j]); }
                    n1 = std::max(n1, c); ni = std::max(ni, d);
                }
                if (b == 1 || n1 * ni <= 1e10) { std::memcpy(D, inv, sizeof(double) * b2); inverted = true; }
                else ++nFallback_;
            } else if (b == 1) { std::memcpy(D, B, sizeof(double) * b2); }
            else ++nFallback_;
        }
        // D^T w = e0  ->  w = (D^-1)^T e0 = first row of D^-1
        double* w = w_.data() + size_t(i) * b;
        if (inverted || invertBlock(D, b)) {
            for (int e = 0; e < b; ++e) w[e] = D[0 * b + e];
        } else {
            double s = 0; for (int e = 0; e < b; ++e) s += S.blk(S.diag[i])[e * b + 0];
            for (int e = 0; e < b; ++e) w[e] = std::fabs(s) > 1e-300 ? 1.0 / s : 1.0;
        }
        double wm = 0; for (int e = 0; e < b; ++e) wm = std::max(wm, std::fabs(w[e]));
        if (!(wm > 0) || !std::isfinite(wm)) for (int e = 0; e < b; ++e) w[e] = 1.0;
    }
    nFallback_ = nFb;
    if (releaseAcc) std::vector<double>().swap(*releaseAcc);
    // pressure matrix: Ap_ij = sum_e w_ie S_ij[e][0]
    CSR Ap; Ap.n = n; Ap.rowPtr = S.rowPtr; Ap.col = S.col; Ap.val.assign(S.col.size(), 0.0);
    auto prow = [&](int i) {
        const double* w = w_.data() + size_t(i) * b;
        for (int p = S.rowPtr[i]; p < S.rowPtr[i + 1]; ++p) {
            const double* a = S.blk(p);
            double s = 0; for (int e = 0; e < b; ++e) s += w[e] * a[e * b];
            Ap.val[p] = s;
        }
    };
    if (useParallel(n > kPar)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) prow(i);
    } else for (int i = 0; i < n; ++i) prow(i);
    Ap.finalizeDiag();
    reused_ = reuseAmg && amg_.refresh(Ap);
    if (!reused_) amg_.setup(std::move(Ap));   // C45: the AMG owns the fine pressure matrix (was a second copy in Ap_)
    // C38: mixed precision uses the level-parallel BILU0 with float32 factors (the serial BILU0f only where the
    // contiguous level path is off)
    ilu_.setFloatFactors(mixed_);
    if (ilu_.lean()) { ilu_.factor(S); useF_ = false; }   // C45: float factors on every path, no double copy
    else if (mixed_) { ilu_.factor(S); useF_ = !ilu_.contiguousLevels(); if (useF_) iluf_.factor(S); }
    else { ilu_.factor(S); useF_ = false; }
    rp_.assign(n, 0); xp_.assign(n, 0); z1_.assign(size_t(n) * b, 0); t_.assign(size_t(n) * b, 0);
}

void BlockCPR::apply(const double* r, double* z) const {
    const BSR& S = *S_;
    const int n = S.n, b = S.nb;
    for (int i = 0; i < n; ++i) {
        double s = 0; for (int e = 0; e < b; ++e) s += w_[size_t(i) * b + e] * r[size_t(i) * b + e];
        rp_[i] = s;
    }
    amg_.apply(rp_.data(), xp_.data());
    // second stage on the residual r - S z1 with z1 = (xp, 0, ..): only the pressure column of each block contributes
    // (C12: 1 / nb of the full block SpMV)
#pragma omp parallel for schedule(static) if (n > 20000)
    for (int i = 0; i < n; ++i) {
        double acc[kMaxBlock];
        for (int e = 0; e < b; ++e) acc[e] = r[size_t(i) * b + e];
        for (int k = S.rowPtr[i]; k < S.rowPtr[i + 1]; ++k) {
            const double* B = S.blk(k);
            const double xk = xp_[size_t(S.col[k])];
            for (int e = 0; e < b; ++e) acc[e] -= B[e * b] * xk;
        }
        for (int e = 0; e < b; ++e) t_[size_t(i) * b + e] = acc[e];
    }
    if (useF_) iluf_.apply(t_.data(), z); else ilu_.apply(t_.data(), z);
    for (int i = 0; i < n; ++i) z[size_t(i) * b] += xp_[i];
}

// ============================================================================ C16 mixed precision / GPU interface
void BILU0f::factor(const BSR& A) {
    BSR LU = A;                                   // factor in double (stable pivots), store in float
    if (LU.diag.size() != size_t(LU.n)) LU.finalizeDiag();
    bilu0InPlace(LU, nReg_);
    n_ = LU.n; nb_ = LU.nb;
    rowPtr_ = LU.rowPtr; col_ = LU.col; diag_ = LU.diag;
    lu_.resize(LU.val.size());
    for (size_t e = 0; e < LU.val.size(); ++e) lu_[e] = float(LU.val[e]);
    w_.assign(size_t(n_) * nb_, 0.f);
}

void BILU0f::apply(const double* r, double* z) const {
    const int n = n_, b = nb_, b2 = b * b;
    float* y = w_.data();
    for (int i = 0; i < n; ++i) {
        float* yi = y + size_t(i) * b;
        for (int e = 0; e < b; ++e) yi[e] = float(r[size_t(i) * b + e]);
        for (int p = rowPtr_[i]; p < diag_[i]; ++p) {
            const float* a = lu_.data() + size_t(p) * b2; const float* yk = y + size_t(col_[p]) * b;
            for (int e = 0; e < b; ++e) { float s = 0; for (int f = 0; f < b; ++f) s += a[e * b + f] * yk[f]; yi[e] -= s; }
        }
    }
    float t[kMaxBlock];
    for (int i = n - 1; i >= 0; --i) {
        float* yi = y + size_t(i) * b;
        for (int p = diag_[i] + 1; p < rowPtr_[i + 1]; ++p) {
            const float* a = lu_.data() + size_t(p) * b2; const float* yj = y + size_t(col_[p]) * b;
            for (int e = 0; e < b; ++e) { float s = 0; for (int f = 0; f < b; ++f) s += a[e * b + f] * yj[f]; yi[e] -= s; }
        }
        const float* d = lu_.data() + size_t(diag_[i]) * b2;
        for (int e = 0; e < b; ++e) { float s = 0; for (int f = 0; f < b; ++f) s += d[e * b + f] * yi[f]; t[e] = s; }
        for (int e = 0; e < b; ++e) yi[e] = t[e];
    }
    for (size_t e = 0; e < size_t(n) * b; ++e) z[e] = double(y[e]);
}

namespace { GpuBlockSolver* g_gpu = nullptr; }
void setGpuBlockSolver(GpuBlockSolver* s) { g_gpu = s; }
GpuBlockSolver* gpuBlockSolver() { return g_gpu; }
namespace { GpuBlockSolver* g_dist = nullptr; }
void setDistributedSolver(GpuBlockSolver* s) { g_dist = s; }
GpuBlockSolver* distributedSolver() { return g_dist; }

bool CpuFloatBlockSolver::solve(const BSR& A, const std::vector<double>& bd, std::vector<double>& xd, double tol, int maxIt, int& iterations) {
    // exactly the device algorithm: float32 BSR SpMV, float32 BILU0 preconditioner, float32 BiCGSTAB vectors with
    // float dot products accumulated in double (the GPU reduction uses a two-level tree with the same effect)
    const int N = A.size(), b = A.nb, b2 = b * b;
    iterations = 0;
    std::vector<float> Av(A.val.size());
    for (size_t e = 0; e < Av.size(); ++e) Av[e] = float(A.val[e]);
    auto spmv = [&](const std::vector<float>& x, std::vector<float>& y) {
        for (int i = 0; i < A.n; ++i) {
            float acc[kMaxBlock] = {};
            for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k) {
                const float* B = Av.data() + size_t(k) * b2; const float* xk = x.data() + size_t(A.col[k]) * b;
                for (int e = 0; e < b; ++e) for (int f = 0; f < b; ++f) acc[e] += B[e * b + f] * xk[f];
            }
            for (int e = 0; e < b; ++e) y[size_t(i) * b + e] = acc[e];
        }
    };
    BILU0f M; MCBILU0 Mc; DeviceCprAmg Dc;
    const bool useCpr = cprAmg && Dc.build(A);
    if (useCpr) {}
    else if (multicolour) Mc.factor(A); else M.factor(A);
    // float copy of the multicolour factors (the device stores them in float32)
    std::vector<float> mcLU;
    if (multicolour && !useCpr) { mcLU.resize(Mc.factors().val.size()); for (size_t e = 0; e < mcLU.size(); ++e) mcLU[e] = float(Mc.factors().val[e]); }
    std::vector<double> tmpD(N), tmpZ(N);
    std::vector<float> yp(static_cast<size_t>(N));
    auto prec = [&](const std::vector<float>& r, std::vector<float>& z) {
        if (useCpr) { Dc.apply(A, Av, r.data(), z.data()); return; }
        if (multicolour) {   // exactly the device kernels: permute, per-colour forward / backward sweeps in float
            const BSR& F = Mc.factors(); const Colouring& C = Mc.colouring();
            for (int k = 0; k < A.n; ++k) for (int e = 0; e < b; ++e) yp[size_t(k) * b + e] = r[size_t(C.perm[size_t(k)]) * b + e];
            for (int c = 0; c < C.colours; ++c)
                for (int i = C.start[size_t(c)]; i < C.start[size_t(c) + 1]; ++i)
                    for (int p = F.rowPtr[size_t(i)]; p < F.diag[size_t(i)]; ++p) {
                        const float* a = mcLU.data() + size_t(p) * b2; const float* yk = yp.data() + size_t(F.col[size_t(p)]) * b;
                        for (int e = 0; e < b; ++e) { float s2 = 0; for (int f = 0; f < b; ++f) s2 += a[e * b + f] * yk[f]; yp[size_t(i) * b + e] -= s2; }
                    }
            for (int c = C.colours - 1; c >= 0; --c)
                for (int i = C.start[size_t(c)]; i < C.start[size_t(c) + 1]; ++i) {
                    float* yi = yp.data() + size_t(i) * b;
                    for (int p = F.diag[size_t(i)] + 1; p < F.rowPtr[size_t(i) + 1]; ++p) {
                        const float* a = mcLU.data() + size_t(p) * b2; const float* yj = yp.data() + size_t(F.col[size_t(p)]) * b;
                        for (int e = 0; e < b; ++e) { float s2 = 0; for (int f = 0; f < b; ++f) s2 += a[e * b + f] * yj[f]; yi[e] -= s2; }
                    }
                    const float* d = mcLU.data() + size_t(F.diag[size_t(i)]) * b2; float t2[kMaxBlock];
                    for (int e = 0; e < b; ++e) { float s2 = 0; for (int f = 0; f < b; ++f) s2 += d[e * b + f] * yi[f]; t2[e] = s2; }
                    for (int e = 0; e < b; ++e) yi[e] = t2[e];
                }
            for (int k = 0; k < A.n; ++k) for (int e = 0; e < b; ++e) z[size_t(C.perm[size_t(k)]) * b + e] = yp[size_t(k) * b + e];
            return;
        }
        for (int e = 0; e < N; ++e) tmpD[e] = r[e];
        M.apply(tmpD.data(), tmpZ.data());
        for (int e = 0; e < N; ++e) z[e] = float(tmpZ[e]);
    };
    auto dot = [&](const std::vector<float>& u, const std::vector<float>& v) { double s = 0; for (int e = 0; e < N; ++e) s += double(u[e]) * v[e]; return s; };
    std::vector<float> x(N, 0.f), r(N), r0(N), p(N, 0.f), v(N, 0.f), ph(N), s(N), sh(N), t(N);
    for (int e = 0; e < N; ++e) r[e] = float(bd[e]);
    r0 = r;
    const double bn = std::sqrt(dot(r, r));
    xd.assign(N, 0.0);
    if (!(bn > 0)) return true;
    double rho = 1, alpha = 1, omega = 1;
    for (int it = 1; it <= maxIt; ++it) {
        iterations = it;
        const double rhoN = dot(r0, r);
        if (rhoN == 0 || !std::isfinite(rhoN)) break;
        const double beta = (rhoN / rho) * (alpha / omega);
        for (int e = 0; e < N; ++e) p[e] = float(r[e] + beta * (p[e] - omega * v[e]));
        prec(p, ph); spmv(ph, v);
        const double r0v = dot(r0, v);
        if (r0v == 0 || !std::isfinite(r0v)) break;
        alpha = rhoN / r0v;
        for (int e = 0; e < N; ++e) s[e] = float(r[e] - alpha * v[e]);
        if (std::sqrt(dot(s, s)) <= tol * bn) { for (int e = 0; e < N; ++e) x[e] += float(alpha * ph[e]); break; }
        prec(s, sh); spmv(sh, t);
        const double tt = dot(t, t);
        omega = tt > 0 ? dot(t, s) / tt : 0;
        for (int e = 0; e < N; ++e) { x[e] += float(alpha * ph[e] + omega * sh[e]); r[e] = float(s[e] - omega * t[e]); }
        rho = rhoN;
        if (std::sqrt(dot(r, r)) <= tol * bn || omega == 0) break;
    }
    for (int e = 0; e < N; ++e) xd[e] = x[e];
    for (int e = 0; e < N; ++e) if (!std::isfinite(xd[e])) return false;
    return true;
}

// ============================================================================ C29 device CPR-AMG
void DevCsr::spmv(const float* x, float* y) const {
    for (int i = 0; i < rows; ++i) {
        float acc = 0;   // float accumulation as the device kernel (nb = 1 BSR SpMV)
        for (int k = rp[size_t(i)]; k < rp[size_t(i) + 1]; ++k) acc += v[size_t(k)] * x[size_t(ci[size_t(k)])];
        y[i] = acc;
    }
}

static DevCsr devFromCsr(const CSR& A) {
    DevCsr D; D.rows = D.cols = A.n; D.rp = A.rowPtr; D.ci = A.col; D.v.assign(A.val.begin(), A.val.end());
    return D;
}
static DevCsr devFromParts(int rows, int cols, const std::vector<int>& rp, const std::vector<int>& ci, const std::vector<double>& v) {
    DevCsr D; D.rows = rows; D.cols = cols; D.rp = rp; D.ci = ci; D.v.assign(v.begin(), v.end());
    return D;
}

bool DeviceCprAmg::build(const BSR& S) {
    n = S.n; nb = S.nb;
    lv.clear(); coarseInv = DevCsr{}; bjInv.clear();
    if (n <= 0 || nb < 1 || nb > kMaxBlock) return false;
    BlockCPR C;
    C.setup(S);
    const auto& w = C.weights();
    const auto& H = C.amg().hierarchy();
    if (H.empty() || w.size() != size_t(n) * nb) return false;
    // W: n x N, E: N x n (pressure = first unknown of each block)
    W.rows = n; W.cols = n * nb; W.rp.resize(size_t(n) + 1); W.ci.resize(size_t(n) * nb); W.v.resize(size_t(n) * nb);
    for (int i = 0; i <= n; ++i) W.rp[size_t(i)] = i * nb;
    for (int i = 0; i < n; ++i) for (int e = 0; e < nb; ++e) { W.ci[size_t(i) * nb + e] = i * nb + e; W.v[size_t(i) * nb + e] = float(w[size_t(i) * nb + e]); }
    E.rows = n * nb; E.cols = n; E.rp.assign(size_t(n) * nb + 1, 0); E.ci.clear(); E.v.clear();
    for (int r = 0; r < n * nb; ++r) {
        if (r % nb == 0) { E.ci.push_back(r / nb); E.v.push_back(1.0f); }
        E.rp[size_t(r) + 1] = int(E.ci.size());
    }
    const size_t L = H.size();
    lv.resize(L);
    for (size_t l = 0; l < L; ++l) {
        const auto& h = H[l];
        lv[l].A = devFromCsr(h.A);
        lv[l].dinv.resize(size_t(h.A.n));
        for (int i = 0; i < h.A.n; ++i) {
            double d = 0;
            for (int k = h.A.rowPtr[size_t(i)]; k < h.A.rowPtr[size_t(i) + 1]; ++k) if (h.A.col[size_t(k)] == i) d = h.A.val[size_t(k)];
            lv[l].dinv[size_t(i)] = std::fabs(d) > 1e-300 ? float(1.0 / d) : 0.0f;
        }
        if (l + 1 < L) {
            lv[l].P = devFromParts(h.A.n, h.nc, h.pRow, h.pCol, h.pVal);
            lv[l].R = devFromParts(h.nc, h.A.n, h.rRow, h.rCol, h.rVal);
        }
    }
    // explicit inverse of the coarsest operator (Gauss-Jordan, double) when small enough; else Jacobi sweeps
    const CSR& Cc = H.back().A;
    const int m = Cc.n;
    if (m > 0 && m <= 1500) {
        std::vector<double> M(size_t(m) * 2 * m, 0.0);
        for (int i = 0; i < m; ++i) {
            for (int k = Cc.rowPtr[size_t(i)]; k < Cc.rowPtr[size_t(i) + 1]; ++k) M[size_t(i) * 2 * m + size_t(Cc.col[size_t(k)])] += Cc.val[size_t(k)];
            M[size_t(i) * 2 * m + size_t(m + i)] = 1.0;
        }
        bool ok = true;
        for (int k = 0; k < m && ok; ++k) {
            int pr = k;
            for (int i = k + 1; i < m; ++i) if (std::fabs(M[size_t(i) * 2 * m + k]) > std::fabs(M[size_t(pr) * 2 * m + k])) pr = i;
            if (std::fabs(M[size_t(pr) * 2 * m + k]) < 1e-300) { ok = false; break; }
            if (pr != k) for (int j = 0; j < 2 * m; ++j) std::swap(M[size_t(k) * 2 * m + j], M[size_t(pr) * 2 * m + j]);
            const double piv = 1.0 / M[size_t(k) * 2 * m + k];
            for (int j = 0; j < 2 * m; ++j) M[size_t(k) * 2 * m + j] *= piv;
            for (int i = 0; i < m; ++i) {
                if (i == k) continue;
                const double f = M[size_t(i) * 2 * m + k];
                if (f == 0) continue;
                for (int j = 0; j < 2 * m; ++j) M[size_t(i) * 2 * m + j] -= f * M[size_t(k) * 2 * m + j];
            }
        }
        if (ok) {
            coarseInv.rows = coarseInv.cols = m; coarseInv.rp.resize(size_t(m) + 1);
            coarseInv.ci.resize(size_t(m) * m); coarseInv.v.resize(size_t(m) * m);
            for (int i = 0; i < m; ++i) {
                coarseInv.rp[size_t(i)] = i * m;
                for (int j = 0; j < m; ++j) { coarseInv.ci[size_t(i) * m + j] = j; coarseInv.v[size_t(i) * m + j] = float(M[size_t(i) * 2 * m + size_t(m + j)]); }
            }
            coarseInv.rp[size_t(m)] = m * m;
        }
    }
    // block-Jacobi second stage
    const int b2 = nb * nb;
    bjInv.resize(size_t(n) * b2);
    std::vector<double> blk(static_cast<size_t>(b2));
    for (int i = 0; i < n; ++i) {
        const double* D = S.blk(S.diag[size_t(i)]);
        std::copy(D, D + b2, blk.begin());
        if (!invertBlock(blk.data(), nb)) { std::fill(blk.begin(), blk.end(), 0.0); for (int e = 0; e < nb; ++e) blk[size_t(e) * nb + e] = 1.0 / std::max(1e-30, std::fabs(D[e * nb + e])); }
        for (int e = 0; e < b2; ++e) bjInv[size_t(i) * b2 + e] = float(blk[size_t(e)]);
    }
    x_.assign(L, {}); b_.assign(L, {}); t_.assign(L, {});
    for (size_t l = 0; l < L; ++l) { const size_t nl = size_t(lv[l].A.rows); x_[l].assign(nl, 0.f); b_[l].assign(nl, 0.f); t_[l].assign(nl, 0.f); }
    rp_.assign(size_t(n), 0.f); z1_.assign(size_t(n) * nb, 0.f); tt_.assign(size_t(n) * nb, 0.f);
    return true;
}

void DeviceCprAmg::apply(const BSR& S, const std::vector<float>& Av, const float* r, float* z) const {
    const int N = n * nb, b2 = nb * nb;
    const size_t L = lv.size();
    W.spmv(r, b_[0].data());
    // V-cycle, written as the device kernel sequence (spmv / diagonal scaling / axpy per level)
    for (size_t l = 0; l + 1 < L; ++l) {
        const Lvl& V = lv[l];
        std::vector<float>& x = x_[l]; std::vector<float>& b = b_[l]; std::vector<float>& t = t_[l];
        const size_t nl = x.size();
        for (size_t i = 0; i < nl; ++i) x[i] = omega * V.dinv[i] * b[i];                       // first sweep from x = 0
        for (int s = 1; s < pre; ++s) { V.A.spmv(x.data(), t.data()); for (size_t i = 0; i < nl; ++i) x[i] += omega * V.dinv[i] * (b[i] - t[i]); }
        V.A.spmv(x.data(), t.data());
        for (size_t i = 0; i < nl; ++i) t[i] = b[i] - t[i];
        V.R.spmv(t.data(), b_[l + 1].data());
    }
    {   // coarsest
        const Lvl& V = lv[L - 1];
        std::vector<float>& x = x_[L - 1]; std::vector<float>& b = b_[L - 1]; std::vector<float>& t = t_[L - 1];
        if (!coarseInv.rp.empty()) coarseInv.spmv(b.data(), x.data());
        else {
            const size_t nl = x.size();
            for (size_t i = 0; i < nl; ++i) x[i] = omega * V.dinv[i] * b[i];
            for (int s = 0; s < 30; ++s) { V.A.spmv(x.data(), t.data()); for (size_t i = 0; i < nl; ++i) x[i] += omega * V.dinv[i] * (b[i] - t[i]); }
        }
    }
    for (size_t l = L - 1; l-- > 0;) {
        const Lvl& V = lv[l];
        std::vector<float>& x = x_[l]; std::vector<float>& b = b_[l]; std::vector<float>& t = t_[l];
        const size_t nl = x.size();
        V.P.spmv(x_[l + 1].data(), t.data());
        for (size_t i = 0; i < nl; ++i) x[i] += t[i];
        for (int s = 0; s < post; ++s) { V.A.spmv(x.data(), t.data()); for (size_t i = 0; i < nl; ++i) x[i] += omega * V.dinv[i] * (b[i] - t[i]); }
    }
    E.spmv(x_[0].data(), z1_.data());
    // second stage: z = z1 + D^-1 (r - A z1)
    for (int i = 0; i < n; ++i) {
        float acc[kMaxBlock] = {};
        for (int k = S.rowPtr[size_t(i)]; k < S.rowPtr[size_t(i) + 1]; ++k) {
            const float* B = Av.data() + size_t(k) * b2; const float* xk = z1_.data() + size_t(S.col[size_t(k)]) * nb;
            for (int e = 0; e < nb; ++e) for (int f = 0; f < nb; ++f) acc[e] += B[e * nb + f] * xk[f];
        }
        for (int e = 0; e < nb; ++e) tt_[size_t(i) * nb + e] = r[size_t(i) * nb + e] - acc[e];
    }
    for (int i = 0; i < n; ++i) {
        const float* D = bjInv.data() + size_t(i) * b2;
        for (int e = 0; e < nb; ++e) {
            float s2 = 0;
            for (int f = 0; f < nb; ++f) s2 += D[e * nb + f] * tt_[size_t(i) * nb + f];
            z[size_t(i) * nb + e] = z1_[size_t(i) * nb + e] + s2;
        }
    }
    (void)N;
}

LinResult refineWithDevice(GpuBlockSolver& dev, const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxOuter, int innerIt) {
    LinResult res;
    const int N = A.size();
    x.assign(N, 0.0);
    std::vector<double> r = b, d, Ax(N);
    double bn = 0; for (double v : b) bn += v * v; bn = std::sqrt(bn);
    if (!(bn > 0)) { res.converged = true; return res; }
    // float32 resolves a residual reduction of ~1e-4 per device solve reliably; the double outer loop accumulates
    const double innerTol = std::max(1e-4, std::min(1e-2, tol * 10));
    for (int outer = 0; outer < maxOuter; ++outer) {
        double rn = 0; for (double v : r) rn += v * v; rn = std::sqrt(rn);
        res.relResidual = rn / bn;
        if (res.relResidual <= tol) { res.converged = true; return res; }
        int its = 0;
        if (!dev.solve(A, r, d, innerTol, innerIt, its)) return res;
        res.iterations += its;
        for (int e = 0; e < N; ++e) x[e] += d[e];
        A.spmv(x.data(), Ax.data());
        for (int e = 0; e < N; ++e) r[e] = b[e] - Ax[e];
    }
    double rn = 0; for (double v : r) rn += v * v; rn = std::sqrt(rn);
    res.relResidual = rn / bn; res.converged = res.relResidual <= tol;
    return res;
}

} // namespace ressim
