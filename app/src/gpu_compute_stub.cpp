// C39: GPU block-solver interface on Linux / macOS. The C16 device solver (gpu_compute.cpp) runs sokol_gfx compute
// shaders on Direct3D 11; the portable renderer is OpenGL 3.3 core (no compute shaders, and macOS stops at GL 4.1), so
// off Windows the GPU solver reports itself unavailable and the engine keeps its CPU (OpenMP) solvers.
#include "gpu_compute.h"
#ifndef _WIN32
namespace gpu {
namespace { int g_prec = 1; }
SolveTiming lastTiming() { return SolveTiming(); }
void setPreconditioner(int mode) { g_prec = mode; }
int preconditioner() { return g_prec; }
std::string benchmark(int, int) { return "GPU solver: not available on this platform (Direct3D 11 compute only) - the CPU solvers are used.\n"; }
bool available(std::string* info) {
    if (info) *info = "GPU compute solver needs Direct3D 11 (Windows); this build uses the CPU solvers";
    return false;
}
std::string deviceInfo() { return "none (Direct3D 11 compute is Windows-only)"; }
ressim::GpuBlockSolver* solver() { return nullptr; }
std::string selfTest(bool& ok) { ok = false; return "GPU self test skipped: no Direct3D 11 compute device on this platform"; }
}  // namespace gpu
#endif
