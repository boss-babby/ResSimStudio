// Block-structured linear solver for the fully implicit black-oil system.
//
// Input: the scalar CSR Jacobian with cell-interleaved unknowns (nCells x nb, pressure first) followed by one
// BHP unknown per well. Internally:
//   * the cell-cell part B is stored block-CSR (nb x nb dense blocks);
//   * wells are folded in by a Schur complement: the preconditioner solves the reduced cell system
//     S = B - C E^-1 D approximately and back-substitutes the BHPs (block-triangular, exact for the wells);
//     the same-cell part of C E^-1 D is added to the diagonal blocks, the perforation-perforation couplings go
//     into the CPR pressure matrix (so AMG "sees" the wells);
//   * CPR: true-IMPES decoupling weights -> smoothed-aggregation AMG V-cycle on the pressure system ->
//     block-ILU(0) on the full system (multiplicative two-stage);
//   * Krylov: right-preconditioned BiCGStab with a restarted flexible GMRES fallback, on the full system.
// The Newton loop treats a non-converged but reduced residual as an inexact Newton step.
#pragma once
#include "ressim/linsolve.h"
#include <vector>

namespace ressim {

struct BlockSolverOptions {
    double tol = 1e-4;          // relative residual target
    int maxIt = 150;            // Krylov iterations (BiCGStab; GMRES gets the same budget)
    int restart = 40;           // GMRES restart length
    bool cpr = true;            // false: block-ILU(0) only
    bool wellsInPressure = true;
    int amgReuse = 0;
    int pressureIts = 8;    // inner AMG-BiCGStab iterations of the pressure stage (0 = one V-cycle)
    double pressureTol = 1e-2;
    bool accelerator = false; // use the registered BlockAccelerator (GPU) as the cell preconditioner
    int krylov = 0;         // 0 auto (CPR -> flexible GMRES, ILU -> BiCGStab), 1 BiCGStab, 2 FGMRES           // reuse the AMG hierarchy for this many further solves when the pattern is unchanged
};

struct BlockSolverStats {
    bool converged = false;
    int iterations = 0;
    double relResidual = 0;
    double setupSec = 0, solveSec = 0;
    int amgLevels = 0;
};

// Optional device accelerator (GPU) for the cell system. The engine stays portable: the application registers an
// implementation (D3D11 compute on Windows). It receives B' (block-CSR, wells folded in) once per linear solve and
// returns an approximate inverse application; the outer flexible GMRES runs in double precision on the CPU, so a
// single-precision device solve acts as mixed-precision iterative refinement.
class BlockAccelerator {
public:
    virtual ~BlockAccelerator() = default;
    virtual const char* name() const = 0;
    // false: device unavailable / failed -> the CPU preconditioner is used
    virtual bool setup(int nCells, int nb, const std::vector<int>& rowPtr, const std::vector<int>& col, const std::vector<double>& blocks) = 0;
    virtual bool apply(const double* r, double* z) = 0; // z ~= B'^-1 r (cells only, nCells * nb)
};

void setBlockAccelerator(BlockAccelerator* a); // nullptr = none (not owned)
BlockAccelerator* blockAccelerator();

class BlockSolver {
public:
    // accBlocks: nCells * nb * nb accumulation derivatives (row = equation, column = unknown) for true-IMPES
    // weights, or nullptr for quasi-IMPES (diagonal Jacobian block).
    BlockSolverStats solve(const CSR& A, int nCells, int nb, int nWells, const double* accBlocks,
                           const std::vector<double>& b, std::vector<double>& x, const BlockSolverOptions& o);
    struct Impl;
    BlockSolver();
    ~BlockSolver();
    BlockSolver(const BlockSolver&) = delete;
    BlockSolver& operator=(const BlockSolver&) = delete;
private:
    Impl* p_;
};

} // namespace ressim
