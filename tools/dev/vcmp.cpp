// Developer tool: compare summary vectors of a ResSim run with an OPM / ECLIPSE reference at the reference times.
//   vcmp DECK REF.SMSPEC KEY[:WELL] ...      (EVERY=n prints every n-th reference time; columns: reference, ResSim)
// Build: g++ -std=c++17 -O1 -fopenmp -Iengine/include tools/dev/vcmp.cpp build/linux/libressim_engine.a -o vcmp -lpthread -ldl
#include "ressim/api.h"
#include "ressim/ecl.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace ressim;
static int find(const ecl::Summary& S, const std::string& k) {
    std::string kw = k, wg; auto p = k.find(':'); if (p != std::string::npos) { kw = k.substr(0, p); wg = k.substr(p + 1); }
    for (size_t v = 0; v < S.vars.size(); ++v) if (S.vars[v].keyword == kw && (wg.empty() || S.vars[v].wgname == wg)) return int(v);
    return -1;
}
static double at(const ecl::Summary& S, int v, double t) {   // rates are step-wise: value of the step containing t
    const auto& T = S.time; const auto& V = S.values[size_t(v)];
    if (T.empty()) return NAN;
    size_t i = 0; while (i < T.size() && T[i] < t) ++i;
    if (i == 0) return V[0]; if (i >= T.size()) return V.back();
    return V[i];
}
int main(int argc, char** argv) {
    if (argc < 4) { std::printf("usage: vcmp DECK REF.SMSPEC KEY[:WELL] ...\n"); return 1; }
    ecl::Summary R; std::string err; if (!ecl::readSummary(argv[2], R, &err)) { std::printf("ERR %s\n", err.c_str()); return 2; }
    api::Session S; if (!S.open(argv[1]) || !S.run()) { std::printf("ERR %s\n", S.error().c_str()); return 3; }
    const ecl::Summary& O = S.summary();
    const int every = std::getenv("EVERY") ? std::atoi(std::getenv("EVERY")) : 1;
    std::printf("%8s", "t");
    for (int a = 3; a < argc; ++a) std::printf(" | %22s", argv[a]);
    std::printf("\n");
    int n = 0;
    for (size_t i = 0; i < R.time.size(); ++i) {
        if (n++ % every) continue;
        const double t = R.time[i];
        std::printf("%8.1f", t);
        for (int a = 3; a < argc; ++a) {
            const int vr = find(R, argv[a]), vo = find(O, argv[a]);
            std::printf(" | %10.4g %10.4g", vr >= 0 ? double(R.values[size_t(vr)][i]) : NAN, vo >= 0 ? at(O, vo, t) : NAN);
        }
        std::printf("\n");
    }
}
