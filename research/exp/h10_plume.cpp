// H10 check: CO2 plume extent (column max of the historical maximum gas saturation) of an H10 deck at given times.
//   h10_plume deck.DATA NX NZ t1 [t2 ...]
#include "ressim/io.h"
#include "ressim/case.h"
#include "ressim/compsim.h"
#include "ressim/comp_deck.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace ressim;
int main(int argc, char** argv) {
    std::vector<std::string> w;
    Case c = importEclipseDeck(argv[1], &w);
    const int nx = std::atoi(argv[2]), nz = std::atoi(argv[3]);
    std::vector<double> ts; for (int a = 4; a < argc; ++a) ts.push_back(std::atof(argv[a]));
    comp::DeckModel m; std::string err;
    if (!comp::deckModel(c, m, &err)) { std::printf("ERROR %s\n", err.c_str()); return 1; }
    std::vector<comp::CompSummary> o;
    std::vector<double> hist;
    comp::runDeckModel(m, o, [&](const comp::CompSummary& S, comp::CompSimulator& sim) {
        { const auto sg = sim.gasSaturation(); if (hist.empty()) hist.assign(sg.size(), 0.0); for (size_t q = 0; q < sg.size(); ++q) hist[q] = std::max(hist[q], sg[q]); }
        for (double t : ts) if (std::fabs(S.time - t) < 1e-6) {
            const auto g = sim.gasSaturation();
            std::printf("t %.0f Sg column max:", S.time);
            int last = -1, lastH = -1;
            for (int i = 0; i < nx; ++i) for (int k = 0; k < nz; ++k) { size_t q = size_t(i + nx * k); if (q < hist.size() && hist[q] > 0.01) lastH = std::max(lastH, i); }
            for (int i = 0; i < nx; ++i) { double mx = 0; for (int k = 0; k < nz; ++k) { size_t q = size_t(i + nx * k); if (q < g.size()) mx = std::max(mx, g[q]); } std::printf(" %.2f", mx); if (mx > 0.01) last = i; }
            std::printf("  | front Sg > 0.01 now: column %d (x = %.0f m), ever: column %d\n", last + 1, (last + 1) * 10.0, lastH + 1);
        }
        return true; }, &err, -1, nullptr);
    return 0;
}
