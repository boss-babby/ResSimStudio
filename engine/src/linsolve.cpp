#include "ressim/linsolve.h"
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ressim {

static const int kParMin = 20000;  // below this size run serial (thread overhead dominates; C33: 8000 -> 20000, a barrier per
                                   // vector op made small runs 10x slower whenever another process shared the cores)

int physicalCores() {
    static const int cores = [] {
        int n = 0;
#ifdef _WIN32
        DWORD len = 0;
        GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &len);
        if (len > 0) {
            std::vector<char> buf(len);
            auto* info = reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buf.data());
            if (GetLogicalProcessorInformationEx(RelationProcessorCore, info, &len))
                for (DWORD off = 0; off < len;) {
                    auto* e = reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buf.data() + off);
                    if (e->Relationship == RelationProcessorCore) ++n;
                    off += e->Size;
                }
        }
#else
        std::ifstream f("/proc/cpuinfo");
        std::set<std::pair<int, int>> ids; std::string line; int phys = 0;
        while (std::getline(f, line)) {
            if (line.rfind("physical id", 0) == 0) phys = std::atoi(line.c_str() + line.find(':') + 1);
            else if (line.rfind("core id", 0) == 0) ids.insert({phys, std::atoi(line.c_str() + line.find(':') + 1)});
        }
        n = int(ids.size());
#endif
        const int logical = int(std::thread::hardware_concurrency());
        // VMs often report every vCPU as one core (core id 0): trust the topology only when it shows >= 2 cores
        if (n <= 1 || (logical > 0 && n > logical)) n = logical;
        return std::max(1, n);
    }();
    return cores;
}

void setNumThreads(int n) {
#ifdef _OPENMP
    // memory-bound sparse kernels gain little from hyper-threads while the extra threads spin in every serial section
    // C33: an explicit OMP_NUM_THREADS wins over the automatic choice (batch / regression runs share the machine)
    if (n <= 0) if (const char* e = std::getenv("OMP_NUM_THREADS"); e && std::atoi(e) > 0) n = std::atoi(e);
    omp_set_num_threads(n > 0 ? n : physicalCores());
#else
    (void)n;
#endif
}
bool useParallel(bool big) {
#ifdef _OPENMP
    return big && omp_get_max_threads() > 1;
#else
    (void)big;
    return false;
#endif
}
int numThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

void CSR::finalizeDiag() {
    diag.assign(n, -1);
    for (int i = 0; i < n; ++i)
        for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p)
            if (col[p] == i) { diag[i] = p; break; }
}

void CSR::spmv(const double* x, double* y) const {
    if (useParallel(n > kParMin)) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p) s += val[p] * x[col[p]];
            y[i] = s;
        }
    } else {
        for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p) s += val[p] * x[col[p]];
            y[i] = s;
        }
    }
}

void ILU0::factor(const CSR& A) {
    A_ = &A;
    lu_ = A.val;
    const int n = A.n;
    std::vector<int> pos(n, -1);
    for (int i = 0; i < n; ++i) {
        const int r0 = A.rowPtr[i], r1 = A.rowPtr[i + 1];
        for (int p = r0; p < r1; ++p) pos[A.col[p]] = p;
        for (int p = r0; p < r1; ++p) {
            const int k = A.col[p];
            if (k >= i) break;                      // columns are sorted
            double piv = lu_[A.diag[k]];
            const double lik = lu_[p] / piv;
            lu_[p] = lik;
            for (int q = A.diag[k] + 1; q < A.rowPtr[k + 1]; ++q) {
                const int j = A.col[q];
                const int t = pos[j];
                if (t >= 0) lu_[t] -= lik * lu_[q];
            }
        }
        double& d = lu_[A.diag[i]];
        if (std::fabs(d) < 1e-30) d = (d >= 0 ? 1e-30 : -1e-30);
        for (int p = r0; p < r1; ++p) pos[A.col[p]] = -1;
    }
}

void ILU0::apply(const double* r, double* z) const {
    const CSR& A = *A_;
    const int n = A.n;
    for (int i = 0; i < n; ++i) {                 // L y = r (unit lower)
        double s = r[i];
        for (int p = A.rowPtr[i]; p < A.diag[i]; ++p) s -= lu_[p] * z[A.col[p]];
        z[i] = s;
    }
    for (int i = n - 1; i >= 0; --i) {            // U z = y
        double s = z[i];
        for (int p = A.diag[i] + 1; p < A.rowPtr[i + 1]; ++p) s -= lu_[p] * z[A.col[p]];
        z[i] = s / lu_[A.diag[i]];
    }
}

// ---------------------------------------------------------------- vector kernels
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    const int n = int(a.size()); double s = 0;
    if (useParallel(n > kParMin)) {
        #pragma omp parallel for reduction(+:s) schedule(static)
        for (int i = 0; i < n; ++i) s += a[i] * b[i];
    } else {
        for (int i = 0; i < n; ++i) s += a[i] * b[i];
    }
    return s;
}
static double norm2(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

LinResult bicgstab(const LinearOperator& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x,
                   double tol, int maxIt) {
    const int n = A.size();
    LinResult res;
    x.assign(n, 0.0);
    std::vector<double> r(b), rh(b), p(n, 0), v(n, 0), s(n), t(n), ph(n), sh(n);
    const double bn = norm2(b);
    if (bn == 0) { res.converged = true; return res; }
    double rho = 1, alpha = 1, omega = 1;
    for (int it = 1; it <= maxIt; ++it) {
        const double rhoN = dot(rh, r);
        if (std::fabs(rhoN) < 1e-300) break;
        const double beta = (rhoN / rho) * (alpha / omega);
        rho = rhoN;
        if (useParallel(n > kParMin)) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
        } else {
            for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        M.apply(p.data(), ph.data());
        A.apply(ph.data(), v.data());
        const double rv = dot(rh, v);
        if (std::fabs(rv) < 1e-300) break;
        alpha = rho / rv;
        if (useParallel(n > kParMin)) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
        } else {
            for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
        }
        const double sn = norm2(s);
        if (sn / bn < tol) {
            if (useParallel(n > kParMin)) {
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
            } else {
                for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
            }
            res.converged = true; res.iterations = it; res.relResidual = sn / bn; return res;
        }
        M.apply(s.data(), sh.data());
        A.apply(sh.data(), t.data());
        const double tt = dot(t, t);
        omega = tt > 0 ? dot(t, s) / tt : 0;
        if (useParallel(n > kParMin)) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) { x[i] += alpha * ph[i] + omega * sh[i]; r[i] = s[i] - omega * t[i]; }
        } else {
            for (int i = 0; i < n; ++i) { x[i] += alpha * ph[i] + omega * sh[i]; r[i] = s[i] - omega * t[i]; }
        }
        const double rn = norm2(r);
        res.iterations = it; res.relResidual = rn / bn;
        if (rn / bn < tol) { res.converged = true; return res; }
        if (omega == 0) break;
    }
    return res;
}

LinResult gmres(const LinearOperator& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x,
                double tol, int maxIt, int m) {
    const int n = A.size();
    LinResult res;
    x.assign(n, 0.0);
    const double bn = norm2(b);
    if (bn == 0) { res.converged = true; return res; }
    std::vector<std::vector<double>> V(m + 1, std::vector<double>(n)), Z(m, std::vector<double>(n));
    std::vector<double> H((m + 1) * m), cs(m), sn(m), g(m + 1), r(n), w(n);
    auto Hat = [&](int i, int j) -> double& { return H[i * m + j]; };
    int total = 0;
    while (total < maxIt) {
        A.apply(x.data(), r.data());
        for (int i = 0; i < n; ++i) r[i] = b[i] - r[i];
        double beta = norm2(r);
        res.relResidual = beta / bn;
        if (beta / bn < tol) { res.converged = true; res.iterations = total; return res; }
        for (int i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        std::fill(g.begin(), g.end(), 0.0); g[0] = beta;
        int j = 0;
        for (; j < m && total < maxIt; ++j, ++total) {
            M.apply(V[j].data(), Z[j].data());
            A.apply(Z[j].data(), w.data());
            for (int i = 0; i <= j; ++i) {
                const double h = dot(w, V[i]);
                Hat(i, j) = h;
                for (int q = 0; q < n; ++q) w[q] -= h * V[i][q];
            }
            const double hn = norm2(w);
            Hat(j + 1, j) = hn;
            if (hn > 0) for (int q = 0; q < n; ++q) V[j + 1][q] = w[q] / hn;
            for (int i = 0; i < j; ++i) {
                const double t = cs[i] * Hat(i, j) + sn[i] * Hat(i + 1, j);
                Hat(i + 1, j) = -sn[i] * Hat(i, j) + cs[i] * Hat(i + 1, j);
                Hat(i, j) = t;
            }
            const double den = std::hypot(Hat(j, j), Hat(j + 1, j));
            cs[j] = den > 0 ? Hat(j, j) / den : 1; sn[j] = den > 0 ? Hat(j + 1, j) / den : 0;
            Hat(j, j) = den; Hat(j + 1, j) = 0;
            g[j + 1] = -sn[j] * g[j]; g[j] = cs[j] * g[j];
            res.relResidual = std::fabs(g[j + 1]) / bn;
            if (res.relResidual < tol) { ++j; ++total; break; }
        }
        // back-substitute and update x
        std::vector<double> y(j);
        for (int i = j - 1; i >= 0; --i) {
            double s = g[i];
            for (int k = i + 1; k < j; ++k) s -= Hat(i, k) * y[k];
            y[i] = s / Hat(i, i);
        }
        for (int i = 0; i < j; ++i) for (int q = 0; q < n; ++q) x[q] += y[i] * Z[i][q];
        if (res.relResidual < tol) { res.converged = true; res.iterations = total; return res; }
    }
    res.iterations = total;
    return res;
}

LinResult bicgstab(const CSR& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt) {
    return bicgstab(CsrOperator(A), M, b, x, tol, maxIt);
}
LinResult gmres(const CSR& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int m) {
    return gmres(CsrOperator(A), M, b, x, tol, maxIt, m);
}

} // namespace ressim
