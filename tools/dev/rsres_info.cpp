// Developer tool (C48): summary of a ResSim results file - engine, summary rows, snapshots, and per array the size and
// min / max at the last snapshot.   rsres_info FILE.rsres
// Build: g++ -std=c++17 -O1 -fopenmp -Iengine/include tools/dev/rsres_info.cpp build/linux/libressim_engine.a -o rsres_info -lpthread -ldl
#include "ressim/results_file.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cstdio>
using namespace ressim;
static void line(const char* n, const std::vector<float>& v) {
    if (v.empty()) return;
    const auto mm = std::minmax_element(v.begin(), v.end());
    std::printf("  %-10s n %zu  min %.6g  max %.6g\n", n, v.size(), double(*mm.first), double(*mm.second));
}
int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: rsres_info FILE.rsres\n"); return 1; }
    const ResultsFile R = loadResultsFile(argv[1]);
    std::printf("engine %s | title %s | summary rows %zu | snapshots %zu\n", R.engine.c_str(), R.title.c_str(), R.summary.size(), R.snaps.size());
    if (R.snaps.empty()) return 0;
    const Snapshot& s = R.snaps.back();
    std::printf("last snapshot t = %.3f d, first t = %.3f d\n", s.time, R.snaps.front().time);
    line("PRESSURE", s.pressure); line("SWAT", s.sw); line("SGAS", s.sg); line("SOIL", s.so); line("XCO2", s.xco2);
    for (const auto& e : s.extra) line(e.first.c_str(), e.second);
    return 0;
}
