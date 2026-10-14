// Linear-solver benchmark on a dumped Jacobian (RESSIM_DUMP=file from the simulator).
#include "ressim/linsolve.h"
#include "ressim/blocksolve.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>

using namespace ressim;

#ifdef GPU_BENCH
void registerGpuSolver();
#endif
int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: lsbench dump.bin [old] [tol]\n");
        return 1;
    }
    FILE* f = std::fopen(argv[1], "rb");
    if (!f) return 2;
    int h[5];
    if (std::fread(h, sizeof(int), 5, f) != 5) return 3;
    CSR A; A.n = h[0];
    const int nnz = h[1], na = h[2], nv = h[3], nw = h[4];
    A.rowPtr.resize(A.n + 1);
    A.col.resize(nnz);
    A.val.resize(nnz);
    std::vector<double> b(A.n);
    size_t ok = std::fread(A.rowPtr.data(), sizeof(int), A.n + 1, f);
    ok += std::fread(A.col.data(), sizeof(int), nnz, f);
    ok += std::fread(A.val.data(), sizeof(double), nnz, f);
    ok += std::fread(b.data(), sizeof(double), A.n, f);
    int hasAcc = 0;
    ok += std::fread(&hasAcc, sizeof(int), 1, f);
    std::vector<double> acc;
    if (hasAcc) {
        acc.resize(size_t(na) * nv * nv);
        ok += std::fread(acc.data(), sizeof(double), acc.size(), f);
    }
    std::fclose(f);
    A.finalizeDiag();
    #ifdef GPU_BENCH
    registerGpuSolver();
#endif
    std::printf("n %d nnz %d cells %d nv %d wells %d acc %d\n", A.n, nnz, na, nv, nw, hasAcc);

    const double tol = argc > 3 ? std::atof(argv[3]) : 1e-4;
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };

    if (argc > 2 && !std::strcmp(argv[2], "old")) {
        auto t0 = now();
        CPR cpr;
        cpr.setup(A, na, nv, nw, hasAcc ? acc.data() : nullptr);
        auto t1 = now();
        std::vector<double> x;
        LinResult r = bicgstab(A, cpr, b, x, tol, 300);
        auto t2 = now();
        std::printf("OLD CPR+BiCGStab: setup %.2fs solve %.2fs its %d res %.2e conv %d levels %d\n",
                    sec(t0, t1), sec(t1, t2), r.iterations, r.relResidual, r.converged, cpr.amg().levels());
        return 0;
    }

    const int nVar = blockAccelerator() ? 5 : 4;
    for (int variant = 0; variant < nVar; ++variant) {
        BlockSolver S;
        BlockSolverOptions o;
        o.tol = tol;
        o.maxIt = 200;
        if (variant == 1) o.wellsInPressure = false;
        if (variant == 2) o.cpr = false;
        if (variant == 4) { o.cpr = false; o.accelerator = true; o.maxIt = 100; }
        const double* ac = hasAcc ? acc.data() : nullptr;
        if (variant == 3) ac = nullptr;
        std::vector<double> x;
        BlockSolverStats st = S.solve(A, na, nv, nw, ac, b, x, o);
        if (variant == 4) std::printf("GPU device: %s\n", blockAccelerator()->name());
        std::printf("NEW v%d (cpr %d wellsP %d): setup %.2fs solve %.2fs its %d res %.2e conv %d levels %d\n",
                    variant, o.cpr, o.wellsInPressure, st.setupSec, st.solveSec, st.iterations, st.relResidual, st.converged, st.amgLevels);
    }
}
