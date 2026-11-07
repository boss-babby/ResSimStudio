// Runs a built-in thermal example (or a deck) on the ResSim thermal core and prints per-report-step results:
//   runtherm <example-id | deck.DATA> [maxSteps]
// (thermal / THMC development driver; link with the engine objects)
#include "ressim/examples.h"
#include "ressim/therm_deck.h"
#include "ressim/io.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace ressim;

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: runtherm <example-id | deck.DATA> [maxSteps]\n"); return 2; }
    const std::string id = argv[1];
    const int maxSteps = argc > 2 ? std::atoi(argv[2]) : -1;
    thermal::DeckModel m;
    std::string err;
    try {
        bool ok;
        if (id.find('.') != std::string::npos) ok = thermal::deckModelFromFile(id, m, &err);
        else ok = thermal::deckModel(makeExample(id), m, &err);
        if (!ok) { std::fprintf(stderr, "model: %s\n", err.c_str()); return 1; }
    } catch (const std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 1; }
    for (const auto& n : m.notes) std::printf("note: %s\n", n.c_str());
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<thermal::Summary> out;
    thermal::Stats st;
    const bool ok = thermal::runDeckModel(m, out, [&](const thermal::Summary& s, const thermal::Simulator& sim) {
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double tmax = 0; for (double T : sim.temperature()) tmax = std::max(tmax, T);
        std::printf("t %7.1f d  FPR %7.2f  Tmax %7.1f C  FOPR %8.3f  FGPR %10.1f  FWIR %8.2f  Qrx %10.4g MW  cumRx %10.4g GJ  %5.1fs",
                    s.time, s.fpr, tmax - 273.15, s.fopr, s.fgpr, s.fwir, sim.reactionHeatRate() * 1e-6, sim.cumReactionHeat() * 1e-9, secs);
        if (const geomech::Solver* g = sim.mechanics()) {
            double fmax = 0; for (double f : g->failure()) fmax = std::max(fmax, f);
            std::printf("  uplift %.4f m (res top %.4f)  fail %.2f", g->surfaceUplift(), g->reservoirTopUplift(), fmax);
        }
        // per-component production (mol/s summed over producers)
        std::printf("\n   prod mol/s:");
        std::vector<double> tot(size_t(sim.nc()) + 1, 0.0);
        for (size_t w = 0; w < sim.wellReports().size(); ++w)
            if (!sim.model().wells[w].injector) for (size_t k = 0; k < tot.size(); ++k) tot[k] += sim.wellReports()[w].molar[k];
        for (size_t k = 0; k < tot.size(); ++k) std::printf(" %s %.3g", k < size_t(sim.nc()) ? sim.model().comps[k].name.c_str() : "H2O", tot[k]);
        std::printf("\n   wells:");
        for (size_t w = 0; w < s.wells.size(); ++w) std::printf(" %s T %.1f H %.4g bhp %.1f |", m.wellNames[w].c_str(), s.wells[w].tfluid, s.wells[w].hrate, s.wells[w].bhp);
        std::printf("\n");
        if (std::getenv("PROFILE")) {   // temperature / gas saturation / porosity profile along the first row
            const int stepc = std::max(1, sim.size() / 24);
            std::printf("   T:"); for (int q = 0; q < sim.size(); q += stepc) std::printf(" %.0f", sim.temperature()[size_t(q)] - 273.15);
            std::printf("\n   Sg:"); for (int q = 0; q < sim.size(); q += stepc) std::printf(" %.2f", sim.sg()[size_t(q)]);
            std::printf("\n   So:"); for (int q = 0; q < sim.size(); q += stepc) std::printf(" %.2f", sim.so()[size_t(q)]);
            const std::vector<double> pf = sim.fluidPorosity();
            std::printf("\n   phiF:"); for (int q = 0; q < sim.size(); q += stepc) std::printf(" %.3f", pf[size_t(q)]);
            std::printf("\n");
        }
        std::fflush(stdout);
        return true;
    }, &err, maxSteps, &st);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%s steps %d newton %d linear %d cuts %d mech %d/%d cg %d  %.1f s\n", ok ? "OK" : ("FAILED: " + err).c_str(), st.steps, st.newton,
                st.linear, st.cuts, st.mechSolves, st.mechIters, st.mechCg, secs);
    return ok ? 0 : 1;
}
