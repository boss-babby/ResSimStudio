// rft_dump: lists an ECLIPSE RFT file (time, well, connections with P / SWAT / SGAS). Usage: rft_dump FILE.RFT [WELL]
#include "ressim/ecl.h"
#include <cstdio>
#include <string>
using namespace ressim::ecl;
int main(int argc, char** argv) {
    File f; std::string err; if (!f.open(argv[1], &err)) { std::printf("%s\n", err.c_str()); return 2; }
    const std::string want = argc > 2 ? argv[2] : "";
    for (int n = 0;; ++n) {
        const int it = f.find("TIME", n); if (it < 0) break;
        const double t = f.numeric(it)[0];
        const int iw = f.find("WELLETC", 0, it); auto we = f.strings(iw);
        const std::string well = we.size() > 1 ? we[1] : "?", type = we.size() > 5 ? we[5] : "";
        const int ii = f.find("CONIPOS", 0, it), ij = f.find("CONJPOS", 0, it), ik = f.find("CONKPOS", 0, it);
        const int ip = f.find("PRESSURE", 0, it), is = f.find("SWAT", 0, it), ig = f.find("SGAS", 0, it), id = f.find("DEPTH", 0, it);
        const int nxt = f.find("TIME", n + 1);
        auto ok = [&](int i) { return i >= 0 && (nxt < 0 || i < nxt); };
        auto I = ok(ii) ? f.ints(ii) : std::vector<int>();
        std::printf("t %8.1f %-8s %s  %zu conns\n", t, well.c_str(), type.c_str(), I.size());
        if (!want.empty() && well == want && ok(ip)) {
            auto J = f.ints(ij), K = f.ints(ik); auto P = f.numeric(ip), D = f.numeric(id);
            auto S = ok(is) ? f.numeric(is) : std::vector<double>(), G = ok(ig) ? f.numeric(ig) : std::vector<double>();
            for (size_t c = 0; c < I.size(); ++c) std::printf("    %3d %3d %3d  D %.2f  P %.3f  SW %.4f  SG %.4f\n", I[c], J[c], K[c], D[c], P[c], c < S.size() ? S[c] : -1, c < G.size() ? G[c] : -1);
        }
    }
}
