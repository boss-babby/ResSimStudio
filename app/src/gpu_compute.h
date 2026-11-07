// C16 GPU linear-solver backend: sokol_gfx compute shaders (D3D11 backend on Windows) implementing the engine's
// GpuBlockSolver interface - float32 BiCGSTAB with block-Jacobi preconditioning on the device, called by the engine
// inside double-precision iterative refinement. Uses its own D3D11 device (independent of the UI renderer), so a
// simulation thread can use it while the UI renders.
#pragma once
#include <string>

namespace ressim { class GpuBlockSolver; }

namespace gpu {
// C18: timing of the last device solve (one refinement inner solve). deviceBusy is measured with D3D11 timestamp
// queries (GPU clock, sum over compute passes); total / setup / iterate are wall-clock on the CPU.
struct SolveTiming {
    bool multicolour = false, converged = false;
    int colours = 0, iterations = 0, readbacks = 0;
    double setup = 0, iterate = 0, total = 0, deviceBusy = -1;   // seconds (deviceBusy -1 = timestamps unavailable)
};
SolveTiming lastTiming();
// preconditioner on the device: 1 = multicolour block ILU(0) (default), 0 = block Jacobi
void setPreconditioner(int mode);
int preconditioner();
// benchmark: CPU BILU0 / CPU MC-BILU0 (OpenMP) vs GPU block Jacobi / GPU MC-BILU0 on block systems of growing size
// (each solved to double accuracy; GPU inside the engine's iterative refinement). Returns a multi-line table.
std::string benchmark(int maxCells = 400000, int threads = 0);
// initialise on first use; false + reason when no compute-capable device / shader compiler is available
bool available(std::string* info = nullptr);
// adapter name, feature level, backend ("sokol_gfx D3D11 compute")
std::string deviceInfo();
// the solver object to register with ressim::setGpuBlockSolver (lives for the program lifetime)
ressim::GpuBlockSolver* solver();
// self test: random block system, device vs CPU double; returns a one-line report, ok = refined error < 1e-8
std::string selfTest(bool& ok);
} // namespace gpu
