// H2 driver: run a deck, print fracture growth events "t_days side length KI_MPa" and the wall time.
#include "ressim/simulator.h"
#include "ressim/io.h"
#include <chrono>
#include <cstdio>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 2) return 1;
    std::vector<std::string> warn;
    Case c = importEclipseDeck(argv[1], &warn);
    c.tuning.threads = 1;
    Simulator sim; sim.setup(c); sim.initialize();
    const auto t0 = std::chrono::steady_clock::now();
    while (sim.runReportStep()) {}
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    for (const auto& e : sim.hfGrowthEvents()) std::printf("EV %.6f %d %.4f %.4f\n", e.t, e.side, e.length, e.KI * 1e-6);
    const auto& st = sim.stats();
    std::printf("DONE steps %d newton %d cuts %d  FPR %.4f  %.2f s\n", st.timesteps, st.newton, st.cuts, sim.summary().back().fpr, sec);
}
