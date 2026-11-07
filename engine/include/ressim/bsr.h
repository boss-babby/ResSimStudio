#pragma once
// Block linear algebra for coupled multiphysics systems.
//
//  * BSR       - block compressed sparse row matrix (nb x nb dense blocks, row-major), the native layout of a
//                fully-implicit reservoir Jacobian (one block per cell connection, nb = unknowns per cell). The same
//                layout serves black-oil (nb = 2..3), solvent / polymer / thermal (nb = 4..5) and compositional
//                (nb = nc + 1) systems and maps one-to-one onto GPU block kernels.
//  * BILU0     - block incomplete LU(0): every diagonal block is inverted exactly (pivoted LU), so cells with a
//                vanishing phase (gas cap, water zone, state-switched cells) whose scalar diagonal is ~0 do not
//                break the factorisation as scalar ILU(0) does.
//  * WellSchur - exact Schur complement of the well (BHP) unknowns: the cell system becomes a pure block matrix
//                S = Acc - Acw Dww^-1 Awc (wells couple all their perforated cells), the BHPs are recovered after
//                the solve. Removes the badly scaled well rows from every preconditioner.
//  * BlockCPR  - two-stage constrained-pressure-residual preconditioner on a BSR system: quasi- or true-IMPES
//                decoupling -> AMG V-cycle on the pressure system -> BILU0 on the full block system.
#include "ressim/linsolve.h"
#include <vector>

namespace ressim {

constexpr int kMaxBlock = 16;         // largest supported block size (dense kernels use stack storage)

struct BSR {
    int n = 0, nb = 1;                // block rows (= block columns), block size
    std::vector<int> rowPtr, col, diag;
    std::vector<double> val;          // rowPtr[n] blocks of nb*nb, row-major (row = equation, column = unknown)
    int nnzb() const { return rowPtr.empty() ? 0 : rowPtr[n]; }
    int size() const { return n * nb; }
    double* blk(int k) { return val.data() + size_t(k) * nb * nb; }
    const double* blk(int k) const { return val.data() + size_t(k) * nb * nb; }
    void finalizeDiag();              // columns must be sorted within each row
    void spmv(const double* x, double* y) const;
    // scalar CSR with the same entries (debugging, AMG of nb = 1 systems)
    CSR toCSR() const;
    // block view of a scalar CSR whose first n*nb rows/columns are cell unknowns grouped in blocks of nb
    static BSR fromCSR(const CSR& A, int n, int nb);
};

struct BsrOperator : LinearOperator {
    const BSR& A;
    explicit BsrOperator(const BSR& a) : A(a) {}
    int size() const override { return A.size(); }
    void apply(const double* x, double* y) const override { A.spmv(x, y); }
};

// dense block helpers (row-major nb x nb)
bool invertBlock(double* a, int nb);                       // in place, partial pivoting; false if singular
double blockConditionEstimate(const double* a, int nb);    // 1-norm condition number (inf if singular)

class BILU0 : public Preconditioner {
public:
    void factor(const BSR& A);
    void apply(const double* r, double* z) const override;
    int regularisedBlocks() const { return nReg_; }        // singular diagonal blocks that were shifted
private:
    BSR LU_;                          // strictly lower blocks hold L (unit diagonal), diagonal blocks hold U_ii^-1
    int nReg_ = 0;
    // C21 level scheduling: rows grouped by dependency depth (forward = L part, backward = U part). Rows of one level are
    // independent, so factorisation and both sweeps run in parallel per level with the exact serial arithmetic
    // (bitwise identical result, no extra iterations - unlike the multicolour reordering).
    std::vector<int> fOrd_, fStart_, bOrd_, bStart_;
    bool levelPar_ = false;
    void buildLevels();
    // C28 (C23-b) level-ordered contiguous storage: the matrix is renumbered by forward level (a topological order of the
    // L dependencies, so every coupled pair keeps its relative order and ILU(0) is mathematically unchanged - same
    // iterations). Rows of a level are then contiguous in memory; with a structurally symmetric pattern the same levels
    // in reverse serve the backward sweep. Factorisation and both sweeps run level-parallel in the permuted space; apply
    // gathers / scatters the vectors once. RESSIM_LEVEL_ILU = 2 forces it on, 0 off; default: see buildLevels.
public:
    bool contiguousLevels() const { return contig_; }
    // C38: keep the level-ordered factors in float32 after the (double) factorisation - half the memory traffic of the
    // level-parallel sweeps and half the steady-state storage; vectors and arithmetic stay double (only the factors are
    // rounded, relative perturbation ~6e-8). Only the contiguous path uses it; the other paths stay double.
    void setFloatFactors(bool on) { floatStore_ = on; }
    bool floatFactors() const { return contig_ && (floatStore_ || lean_) && !Pf_.empty(); }
    size_t factorBytes() const { return contig_ ? (floatFactors() ? Pf_.size() * sizeof(float) : P_.val.size() * sizeof(double)) : (!Lf_.empty() ? Lf_.size() * sizeof(float) : LU_.val.size() * sizeof(double)); }
    // C45 lean mode (large models): the factorisation itself works in float32 storage - each row is eliminated in a
    // double scratch row against the already rounded factors of the earlier rows - so no double copy of the matrix
    // ever exists (the C38 float factors were rounded only after a double factorisation: same peak memory). Both
    // the contiguous level-parallel path and the serial path; vectors and the outer Krylov method stay double.
    void setLean(bool on) { lean_ = on; }
    bool lean() const { return lean_; }
    int forwardLevels() const { return int(fStart_.size()) - 1; }
private:
    bool contig_ = false;
    std::vector<int> perm_, valMap_;  // new -> old row; new block entry -> old block entry
    BSR P_;                           // permuted factors
    bool floatStore_ = false;
    std::vector<float> Pf_;           // C38 float32 copy of P_.val (P_.val freed between factorisations)
    bool lean_ = false;
    std::vector<float> Lf_;           // C45 lean float factors in the original numbering (non-contiguous path)
    mutable std::vector<double> rp_, zp_;
    void buildContiguous();
};

// C18 multicolour ordering: greedy distance-1 colouring of the block graph (largest-degree-first). Cells of one colour
// share no matrix entry, so in the colour-permuted matrix P A P^T every row of a colour depends only on earlier colours
// in the lower factor (later colours in the upper factor): ILU(0) factorisation and both triangular sweeps run in
// parallel over the rows of a colour (OpenMP on the CPU, one compute dispatch per colour on the GPU).
struct Colouring {
    int colours = 0;
    std::vector<int> colour;          // per original cell
    std::vector<int> perm;            // new index -> original cell (cells grouped by colour)
    std::vector<int> iperm;           // original cell -> new index
    std::vector<int> start;           // colour c occupies new indices [start[c], start[c+1])
};
Colouring colourBlocks(const BSR& A);
BSR permuteBsr(const BSR& A, const Colouring& C);   // P A P^T (columns sorted per row)
class MCBILU0 : public Preconditioner {
public:
    bool keepPermuted = false;        // also keep P A P^T (the GPU solver iterates in the permuted space)
    const BSR& permuted() const { return PA_; }
    void factor(const BSR& A);        // colouring (reused while the pattern is unchanged) + permuted block ILU(0)
    void apply(const double* r, double* z) const override;
    int colours() const { return C_.colours; }
    const Colouring& colouring() const { return C_; }
    const BSR& factors() const { return LU_; }    // permuted: L (unit diag), U_ii^-1 on the diagonal, U
    int regularisedBlocks() const { return nReg_; }
private:
    Colouring C_; BSR LU_, PA_; int nReg_ = 0;
    std::vector<int> patRow_, patCol_;
    mutable std::vector<double> rp_, zp_;
};

// C16 mixed precision: the same block ILU(0) stored and applied in float32 (half the memory traffic of the
// triangular sweeps, the layout of the GPU kernels). The outer Krylov method stays in double, so the converged
// solution is unchanged - only the preconditioner is approximated (relative perturbation ~1e-7).
class BILU0f : public Preconditioner {
public:
    void factor(const BSR& A);
    void apply(const double* r, double* z) const override;
    int regularisedBlocks() const { return nReg_; }
    size_t bytes() const { return lu_.size() * sizeof(float); }
private:
    int n_ = 0, nb_ = 1, nReg_ = 0;
    std::vector<int> rowPtr_, col_, diag_;
    std::vector<float> lu_;
    mutable std::vector<float> w_;
};
// C16 GPU offload: an external backend (the app registers a sokol_gfx / D3D11 compute implementation) may take over
// the Krylov solve of the Schur-reduced block system. The device works in float32; the engine wraps every device
// solve in double-precision iterative refinement (defect correction), so the Newton update keeps double accuracy.
class GpuBlockSolver {
public:
    virtual ~GpuBlockSolver() = default;
    virtual const char* name() const = 0;
    // approximately solve A x = b (relative residual tol in float32); false on device failure
    virtual bool solve(const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& iterations) = 0;
    // C18: running total of measured device busy time (GPU timestamp queries) in seconds; < 0 when not measurable
    virtual double deviceSeconds() const { return -1; }
};
void setGpuBlockSolver(GpuBlockSolver* s);    // nullptr = none (owned by the caller)
GpuBlockSolver* gpuBlockSolver();
// C18: double-precision distributed solver hook (MPI domain decomposition, see mpi_dd.h). When set, the black-oil
// Newton loop hands every Schur-reduced block system to it (CPU path as fallback when it does not converge).
void setDistributedSolver(GpuBlockSolver* s);
GpuBlockSolver* distributedSolver();
// C29 device CPR-AMG: a two-stage CPR preconditioner written only in device primitives - scalar / block SpMV, diagonal
// scaling and axpy - so that the whole V-cycle runs on the GPU (no triangular solves, no host round trips):
//   r_p = W r (quasi-IMPES weights, scalar n x N SpMV); x_p = V-cycle(A_p, r_p) with damped-Jacobi smoothing (omega 2/3,
//   pre / post sweeps), plain-aggregation prolongation P / restriction R = P^T of the host AMG hierarchy (amg.cpp) and the
//   explicit dense inverse of the coarsest operator (SpMV); z1 = E x_p (pressure slot); z = z1 + D_b^-1 (r - A z1)
//   (block-Jacobi second stage). Values float32 as on the device. Built on the host (double), uploaded by gpu_compute.cpp,
//   emulated exactly by CpuFloatBlockSolver::cprAmg (tests).
struct DevCsr { int rows = 0, cols = 0; std::vector<int> rp, ci; std::vector<float> v; void spmv(const float* x, float* y) const; };
struct DeviceCprAmg {
    struct Lvl { DevCsr A, P, R; std::vector<float> dinv; };
    DevCsr W, E;                       // pressure restriction (n x N) and prolongation (N x n)
    std::vector<Lvl> lv;               // lv[0].A = pressure matrix; lv.back() = coarsest (no P / R)
    DevCsr coarseInv;                  // dense inverse of the coarsest operator
    std::vector<float> bjInv;          // block-diagonal inverses of the full system (n * nb * nb)
    int n = 0, nb = 1, pre = 1, post = 1;
    float omega = 0.6667f;
    bool build(const BSR& S);          // false: no usable hierarchy (caller falls back)
    // the device preconditioner in float32 (r, z of size n * nb; Av = float values of S)
    void apply(const BSR& S, const std::vector<float>& Av, const float* r, float* z) const;
    mutable std::vector<std::vector<float>> x_, b_, t_;
    mutable std::vector<float> rp_, z1_, tt_;
};
// reference float32 implementation of the device algorithm (BiCGSTAB + BILU0f, all vectors float) - used by the tests
// and as the "GPU emulation" backend when no device is present
class CpuFloatBlockSolver : public GpuBlockSolver {
public:
    bool multicolour = false;          // C18: emulate the device MC-BILU0 (colour-permuted) instead of natural BILU0
    bool cprAmg = false;               // C29: emulate the device CPR-AMG preconditioner (DeviceCprAmg)
    const char* name() const override { return cprAmg ? "CPU float32 emulation (CPR-AMG)" : multicolour ? "CPU float32 emulation (MC-BILU0)" : "CPU float32 emulation"; }
    bool solve(const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& iterations) override;
};
// double-precision iterative refinement around a float device solve: x = sum of corrections until ||b - A x|| <= tol ||b||
LinResult refineWithDevice(GpuBlockSolver& dev, const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxOuter, int innerIt);
class WellSchur {
public:
    // A: scalar CSR, rows/cols 0..na*nv-1 are cells (blocks of nv), na*nv..na*nv+nw-1 one BHP unknown per well.
    // Returns false (caller keeps the full system) if a well equation has no usable BHP derivative.
    bool build(const CSR& A, int na, int nv, int nw);
    const BSR& reduced() const { return S_; }
    void reduceRhs(const double* b, double* bc) const;                 // bc = b_c - Acw Dww^-1 b_w
    void recover(const double* b, const double* xc, double* x) const;  // x = [xc ; Dww^-1 (b_w - Awc xc)]
    int fillBlocks() const { return fill_; }                           // blocks added by well coupling
    // C45 block assembly (no scalar copy of the Jacobian): the simulator assembles the cell-cell blocks straight into
    // S_ and keeps only the well columns / well rows in a small scalar CSR Aw (cell rows 0..na*nv-1 hold well columns
    // only, rows na*nv.. are the well rows as in build()). bindCells() sets S_ to the cell pattern (start / cells:
    // sorted block columns per cell, self included); cellRow(k, e) is row e of the block of cell-pattern entry k.
    // buildInPlace() adds the well fill blocks (pattern reused while the couplings keep their structure, else S_ is
    // re-laid out keeping the assembled values) and subtracts Acw Dww^-1 Awc in place - the arithmetic of build(), so
    // the reduced system is bitwise the same. On false, S_ still holds the plain cell blocks (see cellRow / cellBlocks).
    void bindCells(int na, int nv, const std::vector<int>& start, const std::vector<int>& cells);
    double* cellRow(size_t k, int e) { return S_.val.data() + (size_t(map_[k]) * nv_ + size_t(e)) * nv_; }
    const double* cellBlock(size_t k) const { return S_.val.data() + size_t(map_[k]) * nv_ * nv_; }
    void zeroValues() { std::fill(S_.val.begin(), S_.val.end(), 0.0); }
    bool buildInPlace(const CSR& Aw, int nw, const std::vector<int>& start, const std::vector<int>& cells);
    size_t bytes() const;                                              // storage of S_ and the well parts
private:
    std::vector<int> map_;                               // C45: cell-pattern entry -> block of S_
    int na_ = 0, nv_ = 1, nw_ = 0, fill_ = 0;
    BSR S_;
    std::vector<double> dInv_;                           // 1 / D_ww
    std::vector<int> wcPtr_, wcCol_; std::vector<double> wcVal_;   // Awc rows (scalar cell columns)
    std::vector<int> cwPtr_, cwWell_; std::vector<double> cwVal_;  // Acw by scalar cell row
    unsigned long long patHash_ = 0;                     // C37: hash of the scalar pattern + well couplings of S_
    int base_ = 0;                                       // cell blocks without the well fill
};

class BlockCPR : public Preconditioner {
public:
    // accBlocks (optional, n * nb * nb): accumulation derivatives for true-IMPES weights; cells whose block is
    // ill-conditioned (> 1e10) fall back to quasi-IMPES (Jacobian diagonal block) weights.
    // releaseAcc (optional): the vector holding accBlocks, freed right after the weights (C45 memory: before the AMG
    // and BILU0 set-up, which is the peak of the run)
    void setup(const BSR& S, const double* accBlocks = nullptr, bool reuseAmg = false, std::vector<double>* releaseAcc = nullptr);
    void apply(const double* r, double* z) const override;
    const AMG& amg() const { return amg_; }
    const CSR& pressureMatrix() const { return amg_.fineMatrix(); }
    int trueImpesFallbacks() const { return nFallback_; }
    int bilu0Regularised() const { return ilu_.regularisedBlocks(); }
    bool amgReused() const { return reused_; }
    const std::vector<double>& weights() const { return w_; }   // C29: quasi / true-IMPES weights (n * nb)
    void setMixed(bool m) { mixed_ = m; }     // C16: float32 BILU0 second stage
    void setLean(bool l) { ilu_.setLean(l); } // C45: BILU0 second stage factorised straight into float32 storage
    size_t factorBytes() const { return useF_ ? iluf_.bytes() : ilu_.factorBytes(); }
private:
    const BSR* S_ = nullptr;
    bool mixed_ = false, useF_ = false;
    BILU0f iluf_;
    int nFallback_ = 0; bool reused_ = false;
    std::vector<double> w_;
    AMG amg_;
    BILU0 ilu_;
    mutable std::vector<double> rp_, xp_, z1_, t_;
};

} // namespace ressim
