// dev tool: compare a ResSim snapshot with an ECLIPSE / OPM restart step, per cell (active cells mapped through PORV > 0
// of the reference INIT). usage: initcmp DECK.DATA REF.INIT REF.UNRST [step=0] [ARRAY=SWAT] [extra global arrays from
// the deck: SWL,SWATINIT...]  -> stats + the 15 worst cells
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/ecl.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <string>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: initcmp DECK INIT UNRST [step] [ARRAY]\n"); return 2; }
    const int step = argc > 4 ? std::atoi(argv[4]) : 0; const std::string A = argc > 5 ? argv[5] : "SWAT";
    Case c = importEclipseDeck(argv[1]);
    Simulator s; s.setup(c); s.initialize();
    for (int k = 0; k < step; ++k) s.runReportStep();
    const Snapshot& S = s.snapshots().back();
    ecl::File fi, fr; std::string err;
    if (!fi.open(argv[2], &err) || !fr.open(argv[3], &err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 1; }
    const std::vector<double> porv = fi.numeric("PORV");
    const std::vector<double> ref = fr.numeric(A, step);
    const auto& mine = S.field(A);
    std::vector<int> refAct(porv.size(), -1); int na = 0;
    for (size_t g = 0; g < porv.size(); ++g) if (porv[g] > 0) refAct[g] = na++;
    std::printf("snapshot t=%g, ours %zu active, ref %d active (%zu values)\n", S.time, mine.size(), na, ref.size());
    const Grid& G = s.grid();
    struct D { double d; int g; double a, b; };
    std::vector<D> ds; double sum = 0, sa = 0, sb = 0;
    for (int a = 0; a < int(mine.size()); ++a) {
        const int g = G.globalOf[size_t(a)];
        if (g < 0 || size_t(g) >= refAct.size() || refAct[size_t(g)] < 0) continue;
        const double r = ref[size_t(refAct[size_t(g)])];
        ds.push_back({mine[size_t(a)] - r, g, mine[size_t(a)], r}); sum += std::fabs(mine[size_t(a)] - r); sa += mine[size_t(a)]; sb += r;
    }
    std::printf("%s: mean |diff| %.4g, mean ours %.5g ref %.5g over %zu cells\n", A.c_str(), sum / std::max<size_t>(1, ds.size()), sa / ds.size(), sb / ds.size(), ds.size());
    if (const char* ce = std::getenv("CELLS")) {   // CELLS="i,j,k;i,j,k" (1-based): print these cells
        std::string L = ce; size_t p = 0;
        while (p < L.size()) {
            size_t q = L.find(';', p); if (q == std::string::npos) q = L.size();
            int i = 0, j = 0, k = 0; std::sscanf(L.substr(p, q - p).c_str(), "%d,%d,%d", &i, &j, &k);
            const int g = (i - 1) + c.grid.nx * ((j - 1) + c.grid.ny * (k - 1));
            for (const D& x : ds) if (x.g == g) std::printf("  CELL %d %d %d: ours %.5g ref %.5g\n", i, j, k, x.a, x.b);
            p = q + 1;
        }
    }
    std::sort(ds.begin(), ds.end(), [](const D& x, const D& y) { return std::fabs(x.d) > std::fabs(y.d); });
    const int nx = c.grid.nx, ny = c.grid.ny;
    for (size_t i = 0; i < ds.size() && i < 15; ++i) {
        const int g = ds[i].g;
        auto at = [&](const std::vector<double>& v) { return size_t(g) < v.size() ? v[size_t(g)] : -9.0; };
        std::printf("  cell %d %d %d: ours %.5g ref %.5g  SWL %.4g SWATINIT %.4g SGU %.4g\n", g % nx + 1, (g / nx) % ny + 1, g / (nx * ny) + 1,
                    ds[i].a, ds[i].b, at(c.grid.swl), at(c.grid.swatinit), at(c.grid.sgu));
    }
}
