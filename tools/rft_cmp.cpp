// rft_cmp: runs a deck and compares cell pressure / Sw / Sg with every RFT record (ECLIPSE .RFT) at its time.
// Usage: rft_cmp DECK FILE.RFT [tmax]
#include "ressim/ecl.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace ressim; using namespace ressim::ecl;
int main(int argc, char** argv) {
    File f; std::string err; if (!f.open(argv[2], &err)) { std::printf("%s\n", err.c_str()); return 2; }
    const double tmax = argc > 3 ? std::atof(argv[3]) : 1e9;
    std::vector<std::string> warn; Case cs = importEclipseDeck(argv[1], &warn);
    Simulator s; s.log = [](const std::string&) {}; s.setup(cs); s.initialize();
    const Grid& G = s.grid();
    for (int n = 0;; ++n) {
        const int it = f.find("TIME", n); if (it < 0) break;
        const double t = f.numeric(it)[0]; if (t > tmax) break;
        const int nxt = f.find("TIME", n + 1);
        auto ok = [&](int i) { return i >= 0 && (nxt < 0 || i < nxt); };
        const int iw = f.find("WELLETC", 0, it); auto we = f.strings(iw);
        const int ii = f.find("CONIPOS", 0, it), ij = f.find("CONJPOS", 0, it), ik = f.find("CONKPOS", 0, it);
        const int ip = f.find("PRESSURE", 0, it), is = f.find("SWAT", 0, it), ig = f.find("SGAS", 0, it);
        if (!ok(ip) || !ok(ii)) continue;
        while (s.summary().back().time < t - 1e-6) if (!s.runReportStep()) break;
        const Snapshot& sn = s.snapshots().back();
        auto I = f.ints(ii), J = f.ints(ij), K = f.ints(ik); auto P = f.numeric(ip);
        auto S = ok(is) ? f.numeric(is) : std::vector<double>(), Gs = ok(ig) ? f.numeric(ig) : std::vector<double>();
        double dp = 0, dsw = 0, dsg = 0; int m = 0;
        std::printf("\nRFT %s at %.0f d (ours %.0f d): I J K  P ours/ecl  SW ours/ecl  SG ours/ecl  FIPNUM\n", we.size() > 1 ? we[1].c_str() : "?", t, sn.time);
        for (size_t c = 0; c < I.size(); ++c) {
            const int g = (I[c] - 1) + G.nx * ((J[c] - 1) + G.ny * (K[c] - 1)); const int a = G.activeOf[size_t(g)];
            if (a < 0) { std::printf("  %3d %3d %3d inactive in ours\n", I[c], J[c], K[c]); continue; }
            const double po = sn.pressure[size_t(a)], so = sn.sw[size_t(a)], go = sn.sg[size_t(a)];
            std::printf("  %3d %3d %3d  %7.2f/%-7.2f  %.3f/%.3f  %.3f/%.3f  %d\n", I[c], J[c], K[c], po, P[c], so, c < S.size() ? S[c] : NAN, go, c < Gs.size() ? Gs[c] : NAN,
                        cs.grid.fipnum.empty() ? 0 : cs.grid.fipnum[size_t(g)]);
            dp += po - P[c]; if (c < S.size()) dsw += so - S[c]; if (c < Gs.size()) dsg += go - Gs[c]; ++m;
        }
        if (m) std::printf("  mean diff: P %+.2f  SW %+.4f  SG %+.4f\n", dp / m, dsw / m, dsg / m);
        std::fflush(stdout);
    }
}
