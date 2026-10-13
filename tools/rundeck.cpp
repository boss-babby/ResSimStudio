// rundeck deck.DATA [-export out.DATA] [-steps n] : import, report warnings, export, run with the ResSim engine
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/case.h"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 2) return 1;
    std::string exp; int maxSteps = 1 << 30; bool quiet = false;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-export") && i + 1 < argc) exp = argv[++i];
        else if (!std::strcmp(argv[i], "-steps") && i + 1 < argc) maxSteps = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-q")) quiet = true;
    }
    try {
        std::vector<std::string> warn;
        Case c = importEclipseDeck(argv[1], &warn);
        for (auto& w : warn) std::printf("WARN: %s\n", w.c_str());
        const std::string v = validateCase(c);
        std::printf("validate: %s\n", v.empty() ? "OK" : v.c_str());
        if (!exp.empty()) { std::ofstream o(exp); o << eclipseDeckToString(c); std::printf("exported %s\n", exp.c_str()); }
        Simulator s; s.log = [&](const std::string& m) { if (!quiet) std::printf("  | %s\n", m.c_str()); };
        auto t0 = std::chrono::steady_clock::now();
        s.setup(c); s.initialize();
        int n = 0;
        while (n++ < maxSteps && s.runReportStep()) {}
        const auto& e = s.summary().back();
        std::printf("DONE t=%.1f d  FOPT %.6g FGPT %.6g FWPT %.6g FPR %.6g  newton %d cuts %d  %.1f s\n", e.time, e.fopt, e.fgpt, e.fwpt, e.fpr,
                    s.stats().newton, s.stats().cuts, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    } catch (const std::exception& e) { std::printf("ERROR: %s\n", e.what()); return 2; }
}
