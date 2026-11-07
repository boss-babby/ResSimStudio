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
        std::FILE* csv = nullptr;
        for (int i = 1; i + 1 < argc; ++i) if (std::string(argv[i]) == "-csv") csv = std::fopen(argv[i + 1], "w");
        if (csv) std::fprintf(csv, "time,FOPT,FWPT,FGPT,FPR,FOPR,FWPR,FGPR,FOIP,FGIP,FWIP\n");
        std::FILE* wcsv = nullptr;   // -wcsv FILE: per-well rates / BHP (long format: time,well,mode,open,oil,water,gas,winj,ginj,bhp)
        for (int i = 1; i + 1 < argc; ++i) if (std::string(argv[i]) == "-wcsv") wcsv = std::fopen(argv[i + 1], "w");
        if (wcsv) std::fprintf(wcsv, "time,well,mode,open,oil,water,gas,winj,ginj,bhp\n");
        while (n++ < maxSteps && s.runReportStep()) {
            const auto& q = s.summary().back();
            std::printf("STEP %d t=%.1f d FOPT %.6g FPR %.5g newton %d cuts %d %.0f s\n", n, q.time, q.fopt, q.fpr, s.stats().newton, s.stats().cuts,
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            if (csv) { std::fprintf(csv, "%.4f,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", q.time, q.fopt, q.fwpt, q.fgpt, q.fpr, q.fopr, q.fwpr, q.fgpr, q.oip, q.gip, q.wip); std::fflush(csv); }
            if (wcsv) {
                for (size_t w = 0; w < q.wells.size() && w < s.wellNames().size(); ++w) {
                    const WellRates& r = q.wells[w];
                    std::fprintf(wcsv, "%.4f,%s,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", q.time, s.wellNames()[w].c_str(), int(r.mode), int(r.open), r.oil, r.water, r.gas,
                                 r.waterInj, r.gasInj, r.bhp);
                }
                std::fflush(wcsv);
            }
            std::fflush(stdout);
        }
        if (wcsv) std::fclose(wcsv);
        if (csv) std::fclose(csv);
        const auto& e = s.summary().back();
        std::printf("DONE t=%.1f d  FOPT %.6g FGPT %.6g FWPT %.6g FPR %.6g  newton %d cuts %d  %.1f s\n", e.time, e.fopt, e.fgpt, e.fwpt, e.fpr,
                    s.stats().newton, s.stats().cuts, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        const auto& st = s.stats();
        std::printf("PROFILE  steps %d  newton %d  linear %d  assemble %.2f s  precond %.2f s  krylov %.2f s  update %.2f s  (CPR %d, AMG levels %d)\n",
                    st.timesteps, st.newton, st.linear, st.secAssemble, st.secPrecond, st.secSolve, st.secUpdate, int(st.cpr), st.amgLevels);
    } catch (const std::exception& e) { std::printf("ERROR: %s\n", e.what()); return 2; }
}
