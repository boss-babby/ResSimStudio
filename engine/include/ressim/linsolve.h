// Sparse linear algebra: CSR matrix, ILU(0), smoothed-aggregation AMG and CPR preconditioners,
// BiCGSTAB and restarted GMRES.
// Vector kernels and SpMV are OpenMP-parallel.
#pragma once
#include <vector>

namespace ressim {

struct CSR {
    int n = 0;
    std::vector<int> rowPtr, col, diag;   // diag[i] = position of (i,i)
    std::vector<double> val;
    void spmv(const double* x, double* y) const;
    void finalizeDiag();
};

class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    virtual void apply(const double* r, double* z) const = 0;   // z ~= A^-1 r
};

class ILU0 : public Preconditioner {
public:
    void factor(const CSR& A);
    void apply(const double* r, double* z) const override;   // z = (LU)^-1 r
private:
    const CSR* A_ = nullptr;
    std::vector<double> lu_;
};

// Smoothed-aggregation algebraic multigrid (one V-cycle per call): strength-based greedy aggregation,
// Jacobi-smoothed prolongator, Galerkin coarse operators, symmetric Gauss-Seidel smoothing, dense LU on the
// coarsest level.
class AMG : public Preconditioner {
public:
    void setup(const CSR& A);
    void apply(const double* r, double* z) const override;    // one V-cycle from z = 0
    int levels() const { return int(lv_.size()); }
    int coarsestSize() const { return lv_.empty() ? 0 : lv_.back().A.n; }
    double operatorComplexity() const;
    struct Level {
        CSR A;
        int nc = 0;                                 // columns of P (= size of next level)
        std::vector<int> pRow, pCol; std::vector<double> pVal;   // prolongator (n x nc), CSR
        std::vector<int> rRow, rCol; std::vector<double> rVal;   // restriction = P^T
        mutable std::vector<double> x, b, r;
    };
private:
    void cycle(int l) const;
    std::vector<Level> lv_;
    std::vector<double> lu_; std::vector<int> piv_;  // coarsest dense LU
};

// Two-stage Constrained Pressure Residual preconditioner for block systems with nv unknowns per cell
// (pressure first) followed by nw well (BHP) unknowns: quasi-IMPES decoupling -> AMG pressure solve ->
// ILU(0) on the full system.
class CPR : public Preconditioner {
public:
    // accBlocks (optional, nCells * nv * nv, row = equation, column = unknown): accumulation-term derivatives for
    // true-IMPES weights; without it the Jacobian diagonal block is used (quasi-IMPES)
    void setup(const CSR& A, int nCells, int nv, int nWells, const double* accBlocks = nullptr);
    void apply(const double* r, double* z) const override;
    const AMG& amg() const { return amg_; }
private:
    const CSR* A_ = nullptr;
    int na_ = 0, nv_ = 1, nw_ = 0;
    std::vector<double> w_;                          // na * nv decoupling weights
    CSR Ap_;
    AMG amg_;
    ILU0 ilu_;
    mutable std::vector<double> rp_, xp_, z1_, t_;
};

struct LinResult { bool converged = false; int iterations = 0; double relResidual = 0; };

LinResult bicgstab(const CSR& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x,
                   double tol, int maxIt);
LinResult gmres(const CSR& A, const Preconditioner& M, const std::vector<double>& b, std::vector<double>& x,
                double tol, int maxIt, int restart = 40);

void setNumThreads(int n);   // 0 = hardware default
int  numThreads();
// true when a loop should run in an OpenMP team. Small loops (or a 1-thread setting, e.g. batch runs) stay serial
// without entering the OpenMP runtime, whose per-region cost is high on Windows (libgomp / winpthreads).
bool useParallel(bool big);

} // namespace ressim
