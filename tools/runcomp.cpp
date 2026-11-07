// Runs an ECLIPSE 300 compositional deck with the ResSim compositional core and writes field / well summaries.
// usage: runcomp DECK.DATA [-steps N] [-csv F.csv] [-wcsv W.csv]
#include "ressim/comp_deck.h"
#include "ressim/examples.h"
#include "ressim/io.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace ressim;
int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: runcomp DECK [-steps N] [-csv F.csv] [-wcsv W.csv]\n"); return 2; }
    int maxSteps = -1; const char* csvp = nullptr; const char* wcsvp = nullptr;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-steps") && i + 1 < argc) maxSteps = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-csv") && i + 1 < argc) csvp = argv[++i];
        else if (!std::strcmp(argv[i], "-wcsv") && i + 1 < argc) wcsvp = argv[++i];
    }
    comp::DeckModel m; std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    const std::string src = argv[1];
    if (src.rfind("ex:", 0) == 0 || src.rfind("rt:", 0) == 0) {   // built-in example (rt: = through an ECLIPSE deck round trip)
        Case c = makeExample(src.substr(3));
        if (src[0] == 'r') { const std::string txt = eclipseDeckToString(c); c = eclipseDeckFromString(txt); if (maxSteps == -2) std::printf("%s\n", txt.c_str()); }
        const std::string v = validateCase(c);
        if (!v.empty()) { std::fprintf(stderr, "invalid case: %s\n", v.c_str()); return 1; }
        if (!comp::deckModel(c, m, &err)) { std::fprintf(stderr, "error: %s\n", err.c_str()); return 1; }
    } else if (!comp::deckModelFromFile(argv[1], m, &err)) { std::fprintf(stderr, "error: %s\n", err.c_str()); return 1; }
    std::printf("compositional model: %d cells, %zu connections, %d components, %zu wells, %zu report steps, T %.2f K, %zu separator stages\n",
                m.base.size(), m.base.conns.size(), m.base.fluid.nc(), m.wellNames.size(), m.stepDays.size(), m.base.T, m.base.separator.size());
    for (const auto& n : m.notes) std::printf("  note: %s\n", n.c_str());
    std::FILE* csv = csvp ? std::fopen(csvp, "w") : nullptr;
    std::FILE* wcsv = wcsvp ? std::fopen(wcsvp, "w") : nullptr;
    if (csv) std::fprintf(csv, "time,FOPT,FWPT,FGPT,FPR,FOPR,FWPR,FGPR,FOIP,FGIP,FWIP,FGIT,FWIT,FGIR,FWIR\n");
    if (wcsv) std::fprintf(wcsv, "time,well,mode,open,oil,water,gas,winj,ginj,bhp\n");
    std::vector<comp::CompSummary> out; comp::CompStats st; double tA = 0, tS = 0, tL = 0;
    const bool ok = comp::runDeckModel(m, out, [&](const comp::CompSummary& S, const comp::CompSimulator& sim) {
        std::printf("t %8.1f d  FPR %8.2f  FOPR %9.2f  FGPR %9.2f  FWPR %8.2f  FGIR %9.2f  FWIR %8.2f  FOPT %.5g\n", S.time, S.fpr, S.fopr, S.fgpr, S.fwpr, S.fgir, S.fwir, S.fopt);
        if (m.base.co2store) std::printf("           CO2 in place: mobile %.4g t (structural %.4g), residual %.4g t, dissolved %.4g t, mineral %.4g t; krg=0 %.4g t; NaCl solid %.4g t\n",
                                         S.co2Mob / 1e3, S.co2Struct / 1e3, S.co2Immob / 1e3, S.co2Diss / 1e3, S.co2Min / 1e3, S.co2KrImmob / 1e3, S.saltSolid / 1e3);
        if (m.base.co2sol) std::printf("           CO2: in hydrocarbon phases %.5g t, dissolved in brine %.5g t\n", S.co2Mob / 1e3, S.co2Diss / 1e3);
        std::fflush(stdout);
        tA = sim.tAsm_; tS = sim.tSetup_; tL = sim.tSolve_;
        if (csv) { std::fprintf(csv, "%.4f,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", S.time, S.fopt, S.fwpt, S.fgpt, S.fpr, S.fopr, S.fwpr, S.fgpr,
                                S.foip, S.fgip, S.fwip, S.fgit, S.fwit, S.fgir, S.fwir); std::fflush(csv); }
        if (wcsv) for (size_t w = 0; w < S.wells.size(); ++w) {
            const auto& W = S.wells[w];
            std::fprintf(wcsv, "%.4f,%s,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", S.time, m.wellNames[w].c_str(), W.onBhp ? 1 : 0, W.open ? 1 : 0, W.oil, W.water, W.gas, W.winj, W.ginj, W.bhp);
        }
        return true;
    }, &err, maxSteps, &st);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%s: %d steps, %d Newton, %d linear, %d cuts, %d flashes, %d relaxed, %.1f s\n", ok ? "DONE" : "FAILED", st.steps, st.newton, st.linear, st.cuts, st.flashes, st.relaxed, sec);
    std::printf("time split: assembly+flash %.2f s, CPR setup %.2f s, GMRES %.2f s, other %.2f s\n", tA, tS, tL, sec - tA - tS - tL);
    if (!ok) std::fprintf(stderr, "error: %s\n", err.c_str());
    if (csv) std::fclose(csv);
    if (wcsv) std::fclose(wcsv);
    return ok ? 0 : 1;
}
