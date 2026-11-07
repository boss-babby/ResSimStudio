// C18: single-simulation MPI domain decomposition of the linear solve.
//
// Every rank owns a contiguous, nnz-balanced range of block rows (cells in natural ordering = slabs along the slowest
// grid axis, like ECLIPSE/OPM's default partition). Rank 0 runs the Newton loop and hands each Schur-reduced block
// system to DDSolver::solve; the matrix rows and right-hand side are scattered once per solve (the row pattern only
// when it changes), and all ranks run one distributed BiCGSTAB:
//   * SpMV with a halo exchange of exactly the ghost entries each rank needs (pairwise MPI_Sendrecv rounds),
//   * global dot products by MPI_Allreduce (two at a time),
//   * preconditioner = two-stage CPR: quasi-IMPES pressure residual gathered to rank 0, one AMG V-cycle on the global
//     pressure matrix (the coarse / global level that keeps iteration counts independent of the rank count), scattered
//     back; then block-Jacobi (non-overlapping additive Schwarz) BILU0 of each rank's diagonal block. cpr = false gives
//     the plain block-Jacobi BILU0.
// The solution is gathered on rank 0. Ranks > 0 sit in serveDD() until rank 0 calls shutdown().
//
//   mpiexec -n 4 ResSimCLI deck.DATA --mpi-dd [-o out.rsres]
//
// Scope: the black-oil engine's linear solve (the dominant cost on large models); assembly stays on rank 0 in solve(),
// and is distributed over the ranks with solveCollective() + Simulator::setDistributedAssembly (C28, ResSimCLI --mpi-da).
#pragma once
#include "ressim/bsr.h"
#include <vector>

namespace ressim {
namespace mpi {

struct DDStats {
    int ranks = 1, solves = 0, iterations = 0, patternSetups = 0;
    double scatterSec = 0, solveSec = 0, haloSec = 0;
    long long haloDoubles = 0;     // doubles exchanged in halos (all ranks of rank 0's view: own sends)
};

// contiguous block-row partition with balanced nonzeros: returns parts+1 row offsets
std::vector<int> partitionRows(const BSR& A, int parts);

class DDSolver : public GpuBlockSolver {
public:
    const char* name() const override { return cpr ? "MPI domain decomposition (BiCGSTAB + CPR: global AMG pressure + block-Jacobi BILU0)" : "MPI domain decomposition (BiCGSTAB + block-Jacobi BILU0)"; }
    // rank 0 only (collective with the ranks in serveDD)
    bool solve(const BSR& A, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& iterations) override;
    // C28 (C21-3) distributed assembly: called collectively on ALL ranks, each with a full-size S of which only its own
    // block rows [starts[rank], starts[rank+1]) are assembled (Simulator::setDistributedAssembly). No matrix scatter: the
    // ranks take their rows directly; for CPR only the scalar pressure rows go to rank 0 (1 / nb^2 of the matrix); the
    // solution is all-gathered (every rank gets x for the replicated Newton update).
    bool solveCollective(const BSR& S, const std::vector<double>& b, std::vector<double>& x, const std::vector<int>& starts,
                         double tol, int maxIt, int& iterations);
    void shutdown();               // rank 0: release the serving ranks
    bool cpr = true;               // two-stage CPR: global AMG pressure stage (rank 0) + block-Jacobi BILU0; false = BILU0 only
    DDStats stats;
};

// ranks > 0: serve solves until rank 0 shuts down; returns 0 on a clean shutdown
int serveDD();

} // namespace mpi
} // namespace ressim
