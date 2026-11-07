// C18 MPI domain-decomposed linear solver (see mpi_dd.h)
#include "ressim/mpi_dd.h"
#include "ressim/mpi_rt.h"
#include "ressim/linsolve.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>

namespace ressim {
namespace mpi {

std::vector<int> partitionRows(const BSR& A, int parts) {
    parts = std::max(1, std::min(parts, std::max(1, A.n)));
    std::vector<int> st(static_cast<size_t>(parts) + 1, A.n);
    st[0] = 0;
    const double total = double(A.nnzb()) + A.n;      // +1 per row: vector work
    int p = 1;
    for (int i = 0; i < A.n && p < parts; ++i) {
        const double done = double(A.rowPtr[size_t(i) + 1]) + (i + 1);
        if (done >= total * p / parts) st[size_t(p++)] = i + 1;
    }
    for (int q = 1; q <= parts; ++q) st[size_t(q)] = std::max(st[size_t(q)], st[size_t(q) - 1]);
    st[size_t(parts)] = A.n;
    return st;
}

namespace {

using clk = std::chrono::steady_clock;
double secs(clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); }

constexpr int kTagPat = 701, kTagVal = 702, kTagRhs = 703, kTagX = 704, kTagReq = 705, kTagHalo = 706, kTagPg = 708, kTagPs = 709, kTagWake = 710;

// quasi-IMPES weights of one block row: D^T w = e_0 (Gaussian elimination with partial pivoting); false if singular
bool qimpesWeights(const double* D, int nb, double* w) {
    double M[kMaxBlock][kMaxBlock + 1];
    for (int r = 0; r < nb; ++r) { for (int c = 0; c < nb; ++c) M[r][c] = D[c * nb + r]; M[r][nb] = r == 0 ? 1.0 : 0.0; }
    for (int c = 0; c < nb; ++c) {
        int pv = c;
        for (int r = c + 1; r < nb; ++r) if (std::fabs(M[r][c]) > std::fabs(M[pv][c])) pv = r;
        if (!(std::fabs(M[pv][c]) > 0)) return false;
        if (pv != c) for (int k = 0; k <= nb; ++k) std::swap(M[c][k], M[pv][k]);
        for (int r = 0; r < nb; ++r) if (r != c) {
            const double f = M[r][c] / M[c][c];
            if (f != 0) for (int k = c; k <= nb; ++k) M[r][k] -= f * M[c][k];
        }
    }
    for (int r = 0; r < nb; ++r) w[r] = M[r][nb] / M[r][r];
    return true;
}

struct Hdr { int cmd = 0, n = 0, nb = 0, nnzb = 0, maxIt = 0, patternNew = 0, cpr = 1, pad = 0; double tol = 0; };

// one rank's share of the distributed system
struct Local {
    int P = 1, me = 0, r0 = 0, r1 = 0, nb = 1, nl = 0;     // owned rows [r0, r1), nl = r1 - r0
    std::vector<int> starts;                               // P + 1 row offsets
    std::vector<int> rowPtr, colG, colX;                   // local CSR of blocks; colX: extended local index
    std::vector<double> val;
    std::vector<int> ghosts;                               // global rows of the ghost entries (sorted -> grouped by owner)
    std::vector<int> ghostStart;                           // P + 1 offsets of each owner's ghosts in 'ghosts'
    std::vector<std::vector<int>> sendIdx;                 // per rank: local rows it needs from me
    BSR diagBlk;                                           // local diagonal block (pattern built at setup)
    std::vector<int> diagMap;                              // local entry -> position in diagBlk (-1 = off-block)
    BILU0 prec;
    std::vector<double> sendBuf;
    DDStats* st = nullptr;
    // two-stage CPR: quasi-IMPES pressure rows; the global pressure matrix and its AMG live on rank 0 (which has A)
    bool cpr = true;
    std::vector<double> w;                                 // nl * nb weights
    AMG amg; CSR Ap; bool amgOk = false;
    std::vector<double> rpG, xpG, rp, xp, z1, az1, r2;

    void setupWeights() {
        w.assign(size_t(nl) * nb, 0.0);
        const size_t b2 = size_t(nb) * nb;
        for (int i = 0; i < nl; ++i) {
            int pd = -1;
            for (int p = rowPtr[size_t(i)]; p < rowPtr[size_t(i) + 1]; ++p) if (colX[size_t(p)] == i) { pd = p; break; }
            if (pd < 0 || !qimpesWeights(&val[size_t(pd) * b2], nb, &w[size_t(i) * nb])) { w[size_t(i) * nb] = 1.0; }
        }
    }
    // rank 0: global pressure matrix Ap_ij = w_i^T A_ij e_0 and its AMG hierarchy (refreshed when the pattern is kept)
    void setupPressure(const BSR& A) {
        const int n = A.n, b = A.nb; const size_t b2 = size_t(b) * b;
        const bool samePat = amgOk && Ap.n == n && Ap.rowPtr == A.rowPtr && Ap.col == A.col;
        Ap.n = n; Ap.rowPtr = A.rowPtr; Ap.col = A.col; Ap.val.assign(A.col.size(), 0.0);
        std::vector<double> wi(static_cast<size_t>(b));
        for (int i = 0; i < n; ++i) {
            const int pd = A.diag[size_t(i)];
            if (pd < 0 || !qimpesWeights(A.blk(pd), b, wi.data())) { std::fill(wi.begin(), wi.end(), 0.0); wi[0] = 1.0; }
            for (int p = A.rowPtr[size_t(i)]; p < A.rowPtr[size_t(i) + 1]; ++p) {
                double s = 0; for (int e = 0; e < b; ++e) s += wi[size_t(e)] * A.val[size_t(p) * b2 + size_t(e) * b];
                Ap.val[size_t(p)] = s;
            }
        }
        Ap.finalizeDiag();
        if (!samePat || !amg.refresh(Ap)) amg.setup(Ap);
        amgOk = true;
    }
    // z = M^-1 r, z sized NX (owned + ghosts); multiplicative two-stage: AMG on the pressure residual, then BILU0
    void precApply(const double* r, double* z) {
        if (!cpr) { prec.apply(r, z); return; }
        const size_t N = size_t(nl) * nb, NX = N + ghosts.size() * size_t(nb);
        rp.assign(size_t(std::max(nl, 1)), 0.0); xp.assign(size_t(std::max(nl, 1)), 0.0);
        for (int i = 0; i < nl; ++i) { double s = 0; for (int e = 0; e < nb; ++e) s += w[size_t(i) * nb + e] * r[size_t(i) * nb + e]; rp[size_t(i)] = s; }
        if (me == 0) {
            const int n = starts[size_t(P)];
            rpG.assign(size_t(n), 0.0); xpG.assign(size_t(n), 0.0);
            std::copy(rp.begin(), rp.begin() + nl, rpG.begin());
            for (int q = 1; q < P; ++q) recvDoubles(rpG.data() + starts[size_t(q)], starts[size_t(q) + 1] - starts[size_t(q)], q, kTagPg);
            amg.apply(rpG.data(), xpG.data());
            for (int q = 1; q < P; ++q) sendDoubles(xpG.data() + starts[size_t(q)], starts[size_t(q) + 1] - starts[size_t(q)], q, kTagPs);
            std::copy(xpG.begin(), xpG.begin() + nl, xp.begin());
        } else {
            sendDoubles(rp.data(), nl, 0, kTagPg);
            recvDoubles(xp.data(), nl, 0, kTagPs);
        }
        z1.assign(NX, 0.0); az1.assign(N, 0.0); r2.assign(N, 0.0);
        for (int i = 0; i < nl; ++i) z1[size_t(i) * nb] = xp[size_t(i)];
        spmv(z1, az1);
        for (size_t k = 0; k < N; ++k) r2[k] = r[k] - az1[k];
        prec.apply(r2.data(), z);
        for (size_t k = 0; k < N; ++k) z[k] += z1[k];
    }

    void setupPattern() {
        nl = r1 - r0;
        // ghosts
        ghosts.clear();
        for (int g : colG) if (g < r0 || g >= r1) ghosts.push_back(g);
        std::sort(ghosts.begin(), ghosts.end());
        ghosts.erase(std::unique(ghosts.begin(), ghosts.end()), ghosts.end());
        ghostStart.assign(size_t(P) + 1, 0);
        {
            size_t k = 0;
            for (int q = 0; q < P; ++q) {
                ghostStart[size_t(q)] = int(k);
                while (k < ghosts.size() && ghosts[k] < starts[size_t(q) + 1]) ++k;
            }
            ghostStart[size_t(P)] = int(ghosts.size());
        }
        colX.resize(colG.size());
        for (size_t p = 0; p < colG.size(); ++p) {
            const int g = colG[p];
            colX[p] = (g >= r0 && g < r1) ? g - r0 : nl + int(std::lower_bound(ghosts.begin(), ghosts.end(), g) - ghosts.begin());
        }
        // tell every owner which of its rows I need (pairwise rounds: send my request to dest, receive src's request)
        sendIdx.assign(size_t(P), {});
        for (int d = 1; d < P; ++d) {
            const int dest = (me + d) % P, src = (me - d + P) % P;
            int nReq = ghostStart[size_t(dest) + 1] - ghostStart[size_t(dest)], nIn = 0;
            sendrecvInts(&nReq, 1, dest, &nIn, 1, src, kTagReq);
            std::vector<int> in(static_cast<size_t>(std::max(nIn, 1)));
            sendrecvInts(ghosts.data() + ghostStart[size_t(dest)], nReq, dest, in.data(), nIn, src, kTagReq + 1);
            auto& si = sendIdx[size_t(src)];
            si.resize(size_t(nIn));
            for (int k = 0; k < nIn; ++k) si[size_t(k)] = in[size_t(k)] - r0;
        }
        // local diagonal block pattern
        diagBlk = BSR(); diagBlk.n = nl; diagBlk.nb = nb; diagBlk.rowPtr.assign(1, 0);
        diagMap.assign(colG.size(), -1);
        for (int i = 0; i < nl; ++i) {
            for (int p = rowPtr[size_t(i)]; p < rowPtr[size_t(i) + 1]; ++p)
                if (colX[size_t(p)] < nl) { diagMap[size_t(p)] = int(diagBlk.col.size()); diagBlk.col.push_back(colX[size_t(p)]); }
            diagBlk.rowPtr.push_back(int(diagBlk.col.size()));
        }
        diagBlk.val.assign(diagBlk.col.size() * size_t(nb) * nb, 0.0);
        diagBlk.finalizeDiag();
        if (st) st->patternSetups++;
    }
    void setupValues() {
        const size_t b2 = size_t(nb) * nb;
        for (size_t p = 0; p < colG.size(); ++p)
            if (diagMap[p] >= 0) std::memcpy(&diagBlk.val[size_t(diagMap[p]) * b2], &val[p * b2], sizeof(double) * b2);
        prec.factor(diagBlk);
        if (cpr) setupWeights();
    }
    // x: nl*nb owned entries followed by ghost entries (filled here)
    void halo(std::vector<double>& x) {
        const auto t0 = clk::now();
        for (int d = 1; d < P; ++d) {
            const int dest = (me + d) % P, src = (me - d + P) % P;
            const auto& si = sendIdx[size_t(dest)];
            sendBuf.resize(std::max<size_t>(1, si.size() * size_t(nb)));
            for (size_t k = 0; k < si.size(); ++k) std::memcpy(&sendBuf[k * size_t(nb)], &x[size_t(si[k]) * nb], sizeof(double) * size_t(nb));
            const int g0 = ghostStart[size_t(src)], gn = ghostStart[size_t(src) + 1] - g0;
            sendrecvDoubles(sendBuf.data(), int(si.size()) * nb, dest, x.data() + size_t(nl + g0) * nb, gn * nb, src, kTagHalo);
            if (st) st->haloDoubles += (long long)si.size() * nb;
        }
        if (st) st->haloSec += secs(t0);
    }
    void spmv(std::vector<double>& x, std::vector<double>& y) {   // y (owned) = A x  (x gets its halo)
        halo(x);
        const int b = nb; const size_t b2 = size_t(b) * b;
        for (int i = 0; i < nl; ++i) {
            double* yi = &y[size_t(i) * b];
            for (int e = 0; e < b; ++e) yi[e] = 0;
            for (int p = rowPtr[size_t(i)]; p < rowPtr[size_t(i) + 1]; ++p) {
                const double* a = &val[size_t(p) * b2]; const double* xj = &x[size_t(colX[size_t(p)]) * b];
                for (int e = 0; e < b; ++e) { double s = 0; for (int f = 0; f < b; ++f) s += a[e * b + f] * xj[f]; yi[e] += s; }
            }
        }
    }
    // distributed right-preconditioned BiCGSTAB; all ranks return the same (converged, iterations)
    bool bicgstab(const std::vector<double>& bl, std::vector<double>& xl, double tol, int maxIt, int& its) {
        const size_t N = size_t(nl) * nb, NX = N + ghosts.size() * size_t(nb);
        auto dot2 = [&](const std::vector<double>& a, const std::vector<double>& b1, const std::vector<double>& c, const std::vector<double>& d2, double& o1, double& o2) {
            double v[2] = {0, 0};
            for (size_t k = 0; k < N; ++k) { v[0] += a[k] * b1[k]; v[1] += c[k] * d2[k]; }
            allreduceSumN(v, 2);
            o1 = v[0]; o2 = v[1];
        };
        std::vector<double> r(bl), r0(bl), p(NX, 0.0), v(N, 0.0), s(N, 0.0), t(N, 0.0), ph(NX, 0.0), sh(NX, 0.0);
        xl.assign(N, 0.0);
        double bn = 0, dummy = 0;
        dot2(bl, bl, bl, bl, bn, dummy);
        bn = std::sqrt(bn);
        its = 0;
        if (!(bn > 0)) return true;
        double rho = 1, alpha = 1, omega = 1;
        for (int it = 0; it < maxIt; ++it) {
            double rhoN = 0, rr = 0;
            dot2(r0, r, r, r, rhoN, rr);
            if (std::sqrt(rr) <= tol * bn) return true;
            if (rhoN == 0 || !std::isfinite(rhoN)) return false;
            its = it + 1;
            const double beta = (rhoN / rho) * (alpha / omega);
            for (size_t k = 0; k < N; ++k) p[k] = r[k] + beta * (p[k] - omega * v[k]);
            precApply(p.data(), ph.data());
            spmv(ph, v);
            double r0v = 0;
            dot2(r0, v, r0, v, r0v, dummy);
            if (r0v == 0 || !std::isfinite(r0v)) return false;
            alpha = rhoN / r0v;
            for (size_t k = 0; k < N; ++k) s[k] = r[k] - alpha * v[k];
            precApply(s.data(), sh.data());
            spmv(sh, t);
            double tt = 0, ts = 0;
            dot2(t, t, t, s, tt, ts);
            omega = tt > 0 ? ts / tt : 0.0;
            for (size_t k = 0; k < N; ++k) { xl[k] += alpha * ph[k] + omega * sh[k]; r[k] = s[k] - omega * t[k]; }
            rho = rhoN;
            if (omega == 0) {
                double rn = 0; dot2(r, r, r, r, rn, dummy);
                return std::sqrt(rn) <= tol * bn;
            }
        }
        double rn = 0; dot2(r, r, r, r, rn, dummy);
        return std::sqrt(rn) <= tol * bn;
    }
};

Local g_local;
std::vector<int> g_patRow, g_patCol;   // rank 0: last scattered pattern
Local g_localC;                         // C28 collective (distributed-assembly) solves
std::vector<int> g_patRowC, g_patColC, g_startsC;

// C28: global pressure matrix on rank 0 from the ranks' own quasi-IMPES rows (S full-size on every rank, rows owned)
void setupPressureCollective(Local& L, const BSR& S) {
    const int b = S.nb; const size_t b2 = size_t(b) * b;
    const int e0 = S.rowPtr[size_t(L.r0)], e1 = S.rowPtr[size_t(L.r1)];
    std::vector<double> mine(size_t(std::max(1, e1 - e0)), 0.0);
    std::vector<double> wi(static_cast<size_t>(b));
    for (int i = L.r0; i < L.r1; ++i) {
        const int pd = S.diag[size_t(i)];
        if (pd < 0 || !qimpesWeights(S.blk(pd), b, wi.data())) { std::fill(wi.begin(), wi.end(), 0.0); wi[0] = 1.0; }
        for (int p = S.rowPtr[size_t(i)]; p < S.rowPtr[size_t(i) + 1]; ++p) {
            double s = 0; for (int e = 0; e < b; ++e) s += wi[size_t(e)] * S.val[size_t(p) * b2 + size_t(e) * b];
            mine[size_t(p - e0)] = s;
        }
    }
    if (L.me == 0) {
        const bool samePat = L.amgOk && L.Ap.n == S.n && L.Ap.rowPtr == S.rowPtr && L.Ap.col == S.col;
        L.Ap.n = S.n; L.Ap.rowPtr = S.rowPtr; L.Ap.col = S.col; L.Ap.val.assign(S.col.size(), 0.0);
        std::copy(mine.begin(), mine.begin() + (e1 - e0), L.Ap.val.begin() + e0);
        for (int q = 1; q < L.P; ++q) {
            const int a0 = S.rowPtr[size_t(L.starts[size_t(q)])], a1 = S.rowPtr[size_t(L.starts[size_t(q) + 1])];
            recvDoubles(L.Ap.val.data() + a0, a1 - a0, q, kTagPg + 20);
        }
        L.Ap.finalizeDiag();
        if (!samePat || !L.amg.refresh(L.Ap)) L.amg.setup(L.Ap);
        L.amgOk = true;
    } else sendDoubles(mine.data(), e1 - e0, 0, kTagPg + 20);
}

void receivePattern(Local& L, const Hdr& h) {
    L.P = size(); L.me = rank(); L.nb = h.nb;
    L.starts.assign(size_t(L.P) + 1, 0);
    bcast(L.starts.data(), int(sizeof(int)) * (L.P + 1), 0);
    L.r0 = L.starts[size_t(L.me)]; L.r1 = L.starts[size_t(L.me) + 1];
    const int nl = L.r1 - L.r0;
    L.rowPtr.assign(size_t(nl) + 1, 0);
    recvInts(L.rowPtr.data(), nl + 1, 0, kTagPat);
    L.colG.resize(size_t(std::max(1, L.rowPtr[size_t(nl)])));
    recvInts(L.colG.data(), L.rowPtr[size_t(nl)], 0, kTagPat + 1);
    L.colG.resize(size_t(L.rowPtr[size_t(nl)]));
    L.setupPattern();
}

} // namespace

bool DDSolver::solve(const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& iterations) {
    iterations = 0;
    if (rank() != 0) return false;
    if (A.nb > kMaxBlock || A.n <= 0) return false;   // before any collective: the serving ranks stay in their loop
    Local& L = g_local;
    L.st = &stats;
    const int P = size();
    stats.ranks = P;
    const auto t0 = clk::now();
    Hdr h; h.cmd = 1; h.n = A.n; h.nb = A.nb; h.nnzb = A.nnzb(); h.maxIt = maxIt; h.tol = tol; h.cpr = cpr ? 1 : 0;
    h.patternNew = (g_patRow != A.rowPtr || g_patCol != A.col || L.nb != A.nb) ? 1 : 0;
    for (int q = 1; q < size(); ++q) send(std::string(1, char(h.cmd)), q, kTagWake);   // C21: wakes the sleeping workers
    bcast(&h, int(sizeof(h)), 0);
    const size_t b2 = size_t(A.nb) * A.nb;
    if (h.patternNew) {
        g_patRow = A.rowPtr; g_patCol = A.col;
        L.P = P; L.me = 0; L.nb = A.nb;
        L.starts = partitionRows(A, P);
        bcast(L.starts.data(), int(sizeof(int)) * (P + 1), 0);
        std::vector<int> rp;
        for (int q = 1; q < P; ++q) {
            const int a = L.starts[size_t(q)], e = L.starts[size_t(q) + 1];
            rp.assign(size_t(e - a) + 1, 0);
            for (int i = a; i <= e; ++i) rp[size_t(i - a)] = A.rowPtr[size_t(i)] - A.rowPtr[size_t(a)];
            sendInts(rp.data(), e - a + 1, q, kTagPat);
            sendInts(A.col.data() + A.rowPtr[size_t(a)], A.rowPtr[size_t(e)] - A.rowPtr[size_t(a)], q, kTagPat + 1);
        }
        L.r0 = 0; L.r1 = L.starts[1];
        L.rowPtr.assign(A.rowPtr.begin(), A.rowPtr.begin() + L.r1 + 1);
        L.colG.assign(A.col.begin(), A.col.begin() + A.rowPtr[size_t(L.r1)]);
        L.setupPattern();
    }
    for (int q = 1; q < P; ++q) {
        const int a = L.starts[size_t(q)], e = L.starts[size_t(q) + 1];
        sendDoubles(A.val.data() + size_t(A.rowPtr[size_t(a)]) * b2, int(size_t(A.rowPtr[size_t(e)] - A.rowPtr[size_t(a)]) * b2), q, kTagVal);
        sendDoubles(b.data() + size_t(a) * A.nb, (e - a) * A.nb, q, kTagRhs);
    }
    L.val.assign(A.val.begin(), A.val.begin() + size_t(A.rowPtr[size_t(L.r1)]) * b2);
    std::vector<double> bl(b.begin(), b.begin() + size_t(L.r1) * A.nb), xl;
    L.cpr = cpr;
    L.setupValues();
    if (cpr) L.setupPressure(A);
    stats.scatterSec += secs(t0);
    const auto t1 = clk::now();
    int its = 0;
    const bool ok = L.bicgstab(bl, xl, tol, maxIt, its);
    stats.solveSec += secs(t1);
    x.assign(size_t(A.n) * A.nb, 0.0);
    std::copy(xl.begin(), xl.end(), x.begin());
    for (int q = 1; q < P; ++q) {
        const int a = L.starts[size_t(q)], e = L.starts[size_t(q) + 1];
        recvDoubles(x.data() + size_t(a) * A.nb, (e - a) * A.nb, q, kTagX);
    }
    iterations = its;
    stats.solves++; stats.iterations += its;
    return ok;
}

bool DDSolver::solveCollective(const BSR& S, const std::vector<double>& b, std::vector<double>& x, const std::vector<int>& starts,
                               double tol, int maxIt, int& iterations) {
    iterations = 0;
    const int P = size(), me = rank();
    // every rank sees the same S pattern and partition, so these early exits are taken on all ranks alike
    if (S.nb > kMaxBlock || S.n <= 0 || int(starts.size()) != P + 1 || starts.back() != S.n) return false;
    Local& L = g_localC;
    L.st = &stats;
    stats.ranks = P;
    const auto t0 = clk::now();
    const bool patNew = g_patRowC != S.rowPtr || g_patColC != S.col || g_startsC != starts || L.nb != S.nb;
    const size_t b2 = size_t(S.nb) * S.nb;
    if (patNew) {
        g_patRowC = S.rowPtr; g_patColC = S.col; g_startsC = starts;
        L.P = P; L.me = me; L.nb = S.nb; L.starts = starts;
        L.r0 = starts[size_t(me)]; L.r1 = starts[size_t(me) + 1];
        const int e0 = S.rowPtr[size_t(L.r0)];
        L.rowPtr.assign(size_t(L.r1 - L.r0) + 1, 0);
        for (int i = L.r0; i <= L.r1; ++i) L.rowPtr[size_t(i - L.r0)] = S.rowPtr[size_t(i)] - e0;
        L.colG.assign(S.col.begin() + e0, S.col.begin() + S.rowPtr[size_t(L.r1)]);
        L.amgOk = false;
        L.setupPattern();                                   // collective (ghost requests)
    }
    L.val.assign(S.val.begin() + size_t(S.rowPtr[size_t(L.r0)]) * b2, S.val.begin() + size_t(S.rowPtr[size_t(L.r1)]) * b2);
    std::vector<double> bl(b.begin() + size_t(L.r0) * S.nb, b.begin() + size_t(L.r1) * S.nb), xl;
    L.cpr = cpr;
    L.setupValues();
    if (cpr) setupPressureCollective(L, S);
    stats.scatterSec += secs(t0);
    const auto t1 = clk::now();
    int its = 0;
    const bool ok = L.bicgstab(bl, xl, tol, maxIt, its);
    stats.solveSec += secs(t1);
    x.assign(size_t(S.n) * S.nb, 0.0);
    std::copy(xl.begin(), xl.end(), x.begin() + size_t(L.r0) * S.nb);
    const auto t2 = clk::now();
    for (int q = 0; q < P; ++q) {                            // all-gather of the update (replicated Newton state)
        const int a = starts[size_t(q)], e = starts[size_t(q) + 1];
        if (e > a) bcast(x.data() + size_t(a) * S.nb, int(sizeof(double)) * (e - a) * S.nb, q);
    }
    stats.haloSec += secs(t2);
    iterations = its;
    stats.solves++; stats.iterations += its;
    return ok;
}

void DDSolver::shutdown() {
    if (rank() != 0) return;
    Hdr h; h.cmd = 0;
    for (int q = 1; q < size(); ++q) send(std::string(1, char(h.cmd)), q, kTagWake);   // C21: wakes the sleeping workers
    bcast(&h, int(sizeof(h)), 0);
}

int serveDD() {
    Local& L = g_local;
    DDStats local;
    L.st = &local;
    for (;;) {
        Hdr h;
        std::string wake;
        if (!recv(wake, 0, kTagWake, nullptr, nullptr)) return 3;   // sleeps between solves instead of spinning in MPI_Bcast
        if (!bcast(&h, int(sizeof(h)), 0)) return 3;
        if (h.cmd == 0) return 0;
        if (h.patternNew) receivePattern(L, h);
        const size_t b2 = size_t(L.nb) * L.nb;
        const int nl = L.r1 - L.r0;
        L.val.resize(std::max<size_t>(1, size_t(L.rowPtr[size_t(nl)]) * b2));
        recvDoubles(L.val.data(), int(size_t(L.rowPtr[size_t(nl)]) * b2), 0, kTagVal);
        std::vector<double> bl(std::max<size_t>(1, size_t(nl) * L.nb)), xl;
        recvDoubles(bl.data(), nl * L.nb, 0, kTagRhs);
        bl.resize(size_t(nl) * L.nb);
        L.cpr = h.cpr != 0;
        L.setupValues();
        int its = 0;
        L.bicgstab(bl, xl, h.tol, h.maxIt, its);
        sendDoubles(xl.data(), nl * L.nb, 0, kTagX);
    }
}

} // namespace mpi
} // namespace ressim
