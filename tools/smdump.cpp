// smdump REF.SMSPEC DECK.DATA KEY... : OPM / ECLIPSE reference vs ResSim run side by side at the reference ministeps
// (STRIDE=n prints every n-th ministep, DUMPME=1 also lists the ResSim time steps). Diagnostic tool (C37 / C38).
#include "ressim/ecl.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 4) { std::printf("usage: smdump REF.SMSPEC DECK.DATA KEY...\n"); return 2; }
    ecl::Summary ref; std::string err;
    if (!ecl::readSummary(argv[1], ref, &err)) { std::printf("err %s\n", err.c_str()); return 1; }
    std::vector<DeckMessage> rep;
    Case c = importEclipseDeck(argv[2], nullptr, &rep);
    Simulator s;
    if (std::getenv("LOG")) s.log = [](const std::string& m) { std::printf("log: %s\n", m.c_str()); };
    s.setup(c); s.initialize(); while (s.runReportStep()) {}
    std::vector<std::string> wn; for (auto& w : c.wells) wn.push_back(w.name);
    ecl::Summary me = ecl::summaryFromRun(s.summary(), wn, c.units, ecl::componentNames(c));
    if (std::getenv("DUMPME")) for (auto& p : s.summary()) std::printf("me t=%.4f fopr=%.6g fpr=%.6g\n", p.time, p.fopr, p.fpr);
    const int stride = std::getenv("STRIDE") ? std::max(1, std::atoi(std::getenv("STRIDE"))) : 1;
    std::printf("%10s", "t");
    for (int k = 3; k < argc; ++k) std::printf(" | %16s ref %8s", argv[k], "ressim");
    std::printf("\n");
    for (size_t i = 0; i < ref.time.size(); i += size_t(stride)) {
        std::printf("%10.3f", ref.time[i]);
        for (int k = 3; k < argc; ++k) {
            const int a = ref.find(argv[k]), b = me.find(argv[k]);
            std::printf(" | %20.6g %8.6g", a >= 0 ? ref.values[size_t(a)][i] : NAN, b >= 0 ? me.at(b, ref.time[i]) : NAN);
        }
        std::printf("\n");
    }
}
