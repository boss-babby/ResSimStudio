// norne_regtran: sum of transmissibility between FIPNUM region pairs, ours vs ECLIPSE INIT/EGRID
// (TRANX/TRANY/TRANZ incl. pinch-outs + NNC TRANNNC). Usage: norne_regtran DECK INIT EGRID
#include "ressim/ecl.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include <cmath>
#include <cstdio>
#include <map>
using namespace ressim; using namespace ressim::ecl;
int main(int c, char** v) {
    ResultCase rc; std::string e; openResultCase(v[2], rc, &e);
    std::vector<std::string> warn; Case cs = importEclipseDeck(v[1], &warn);
    Simulator sim; sim.log = [](const std::string&) {}; sim.setup(cs);
    const Grid& G = sim.grid(); const int nx = G.nx, ny = G.ny, nz = G.nz;
    const auto& fip = cs.grid.fipnum;
    auto R = [&](int g) { return fip.empty() ? 1 : fip[size_t(g)]; };
    auto TX = rc.initArray("TRANX"), TY = rc.initArray("TRANY"), TZ = rc.initArray("TRANZ");
    std::map<std::pair<int, int>, double> se, so; std::map<std::pair<int, int>, int> ne, no;
    auto key = [&](int ga, int gb) { int a = R(ga), b = R(gb); return std::make_pair(std::min(a, b), std::max(a, b)); };
    const int N = nx * ny * nz;
    for (int g = 0; g < N; ++g) {
        const int a = G.activeOf[size_t(g)]; if (a < 0) continue;
        const int i = g % nx, j = (g / nx) % ny, k = g / (nx * ny);
        if (i + 1 < nx && G.activeOf[size_t(g + 1)] >= 0 && TX[size_t(a)] > 0) { auto kk = key(g, g + 1); if (kk.first != kk.second) { se[kk] += TX[size_t(a)]; ne[kk]++; } }
        if (j + 1 < ny && G.activeOf[size_t(g + nx)] >= 0 && TY[size_t(a)] > 0) { auto kk = key(g, g + nx); if (kk.first != kk.second) { se[kk] += TY[size_t(a)]; ne[kk]++; } }
        if (TZ[size_t(a)] > 0) { int gb = -1; for (int kb = k + 1; kb < nz; ++kb) { int q = g + (kb - k) * nx * ny; if (G.activeOf[size_t(q)] >= 0) { gb = q; break; } }
            if (gb >= 0) { auto kk = key(g, gb); if (kk.first != kk.second) { se[kk] += TZ[size_t(a)]; ne[kk]++; } } }
    }
    { File eg; eg.open(v[3]); auto n1 = eg.ints("NNC1"), n2 = eg.ints("NNC2"); auto tn = rc.init.numeric("TRANNNC");
      for (size_t q = 0; q < tn.size() && q < n1.size(); ++q) { auto kk = key(n1[q] - 1, n2[q] - 1); if (kk.first != kk.second) { se[kk] += tn[q]; ne[kk]++; } } }
    for (const auto& cn : G.conns) { int ga = G.globalOf[size_t(cn.a)], gb = G.globalOf[size_t(cn.b)]; auto kk = key(ga, gb); if (kk.first != kk.second) { so[kk] += cn.trans; no[kk]++; } }
    std::map<std::pair<int, int>, int> all; for (auto& [k, x] : se) all[k] = 1; for (auto& [k, x] : so) all[k] = 1;
    std::printf("FIPNUM pair   n ecl/ours   sum T ecl / ours   ratio\n");
    for (auto& [k, x] : all) {
        const double a = se[k], b = so[k];
        const bool flag = std::fabs(b - a) > 0.02 * std::max(a, b) && std::max(a, b) > 1e-3;
        std::printf("%3d-%-3d %6d/%-6d %12.5g / %-12.5g %7.4f %s\n", k.first, k.second, ne[k], no[k], a, b, a > 0 ? b / a : NAN, flag ? "<--" : "");
    }
}
