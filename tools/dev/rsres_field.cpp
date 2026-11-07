// Developer tool (C47): print a per-cell field of a ResSim results file (.rsres) at the snapshots nearest to given times,
// as the column-wise maximum over layers (map view along x for 2-D sections) - e.g. the CO2 plume extent.
//   rsres_field FILE.rsres FIELD NX NY NZ t1 [t2 ...]      (FIELD: SGAS, PRESSURE, SWAT, XCO2, ...; all cells active)
// Build: g++ -std=c++17 -O1 -fopenmp -Iengine/include tools/dev/rsres_field.cpp build/linux/libressim_engine.a -o rsres_field -lpthread -ldl
#include "ressim/results_file.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 7) { std::fprintf(stderr, "usage: rsres_field FILE FIELD NX NY NZ t1 [t2 ...]\n"); return 1; }
    const ResultsFile R = loadResultsFile(argv[1]);
    const int nx = std::atoi(argv[3]), ny = std::atoi(argv[4]), nz = std::atoi(argv[5]);
    for (int a = 6; a < argc; ++a) {
        const double t = std::atof(argv[a]);
        size_t best = 0;
        for (size_t s = 0; s < R.snaps.size(); ++s) if (std::fabs(R.snaps[s].time - t) < std::fabs(R.snaps[best].time - t)) best = s;
        if (R.snaps.empty()) { std::printf("no snapshots\n"); return 1; }
        const auto& v = R.snaps[best].field(argv[2]);
        std::printf("t %.1f %s column max:", R.snaps[best].time, argv[2]);
        for (int i = 0; i < nx; ++i) {
            double m = -1e300;
            for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) { const size_t c = size_t(i + nx * (j + ny * k)); if (c < v.size()) m = std::max(m, double(v[c])); }
            std::printf(" %.3g", m);
        }
        std::printf("\n");
    }
    return 0;
}
