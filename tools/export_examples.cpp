// Builds every built-in example, optionally runs it, and writes it as an ECLIPSE deck:
//   export_examples <outdir> [-run] [-t seconds] [id ...]
#include "ressim/examples.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <cctype>
#include <algorithm>
#include <vector>

using namespace ressim;

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: export_examples <outdir> [-run] [-t sec] [id ...]\n"); return 1; }
    const std::string out = argv[1];
    bool run = false; double maxSec = 600; std::vector<std::string> ids;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-run")) run = true;
        else if (!std::strcmp(argv[i], "-t") && i + 1 < argc) maxSec = std::atof(argv[++i]);
        else ids.push_back(argv[i]);
    }
    int fails = 0;
    for (const ExampleInfo& e : listExamples()) {
        if (!ids.empty() && std::find(ids.begin(), ids.end(), e.id) == ids.end()) continue;
        try {
            const Case c = makeExample(e.id);
            const std::string err = validateCase(c);
            std::string name = e.id; for (char& ch : name) ch = char(std::toupper((unsigned char)ch));
            const std::string path = out + "/" + name + ".DATA";
            exportEclipseDeck(c, path);
            // round trip: the exported deck must import again with the same dimensions / wells
            const Case back = importEclipseDeck(path);
            const bool rt = back.grid.nx == c.grid.nx && back.grid.ny == c.grid.ny && back.grid.nz == c.grid.nz &&
                            back.wells.size() == c.wells.size() && back.schedule.size() == c.schedule.size();
            std::printf("%-10s %-13s %-46s %s%s\n", e.id.c_str(), toString(e.level), path.c_str(), err.empty() ? "valid" : err.c_str(),
                        rt ? "  round-trip OK" : "  ROUND-TRIP MISMATCH");
            if (!err.empty() || !rt) ++fails;
            if (!run) continue;
            Simulator s; s.setup(c); s.initialize();
            std::atomic<bool> cancel{false}; s.cancel = &cancel;
            const auto t0 = std::chrono::steady_clock::now();
            std::atomic<bool> done{false};
            std::thread wd([&] { while (!done) { std::this_thread::sleep_for(std::chrono::milliseconds(100));
                if (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > maxSec) cancel = true; } });
            int k = 0; while (s.runReportStep()) ++k;
            done = true; wd.join();
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            const SummaryPoint& S = s.summary().back();
            std::printf("    run: t=%.1f d steps %d/%zu  FPR=%.2f FOPR=%.4g FWCT=%.3f FGOR=%.4g FOPT=%.4g  ts=%d cuts=%d newton=%d  %.1fs\n",
                        S.time, k, c.schedule.size(), S.fpr, S.fopr, S.fwct, S.fgor, S.fopt, s.stats().timesteps, s.stats().cuts, s.stats().newton, sec);
            if (k < int(c.schedule.size())) ++fails;
        } catch (const std::exception& ex) { std::printf("%-10s ERROR %s\n", e.id.c_str(), ex.what()); ++fails; }
        std::fflush(stdout);
    }
    return fails ? 1 : 0;
}
