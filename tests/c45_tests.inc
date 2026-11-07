// C45 tests: low-memory linear algebra (lean BILU0 / block assembly), multisegment cross-flow option (MSWOPTS third value)
// with the restart state version 16
#include "ressim/io.h"
#include <memory>

static void c45SetEnv(const char* k, const char* v) {   // portable setenv / unsetenv (v = nullptr: remove)
#ifdef _WIN32
    _putenv_s(k, v ? v : "");
#else
    if (v) setenv(k, v, 1); else unsetenv(k);
#endif
}

static void testC45LeanSolver() {
    std::printf("== C45 low-memory linear solver (float32 BILU0 factorised in place) against the default\n");
    // the float factors live in the contiguous level-ordered BILU0 (default from 20k rows): forced here on the 20x20x4 demo
    c45SetEnv("RESSIM_LEVEL_ILU", "2");
    auto run = [](int lean, double& fopt, double& fwpt, int& lin) {
        Case c = makeDemoCase();
        if (c.schedule.size() > 8) c.schedule.resize(8);
        c.tuning.leanMemory = lean;
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        fopt = s.summary().empty() ? 0 : s.summary().back().fopt;
        fwpt = s.summary().empty() ? 0 : s.summary().back().fwpt;
        lin = s.stats().linear;
    };
    double f0, w0, f1, w1; int l0, l1;
    run(0, f0, w0, l0); run(1, f1, w1, l1);
    std::printf("    default FOPT %.10g FWPT %.6g (%d linear) | lean FOPT %.10g FWPT %.6g (%d linear), FOPT diff %.3g\n", f0, w0, l0, f1, w1, l1, f1 - f0);
    check(f0 > 0 && std::fabs(f1 - f0) <= 1e-5 * f0 && std::fabs(w1 - w0) <= 1e-4 * (1 + w0),
          "C45: lean BILU0 gives the default result to solver tolerance");
    check(f1 != f0 || l1 != l0, "C45: lean BILU0 is exercised (float32 factors change the round-off)");
    double fa, wa; int la; run(-1, fa, wa, la);
    c45SetEnv("RESSIM_LEVEL_ILU", nullptr);
    check(fa == f0 && wa == w0 && la == l0, "C45: leanMemory auto stays off below 500k cells (bit-identical)");
}

static void testC45MswCrossflow() {
    std::printf("== C45 multisegment cross-flow option (MSW-2D-VERT-02), MSWOPTS round trip, restart v16\n");
    const std::string deck = findDeck("decks/opm_msw/MSW-2D-VERT-02.DATA");
    if (deck.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    Case c = importEclipseDeck(deck);
    check(!c.mswOpt.xflow, "C45: multisegment cross-flow is off by default");
    // deck round trip of the third MSWOPTS value
    {
        Case x = c; x.mswOpt.xflow = true;
        const std::string out = "c45_msw_rt.DATA";
        exportEclipseDeck(x, out);
        std::ifstream f(out); std::stringstream ss; ss << f.rdbuf();
        const bool line = ss.str().find("-- RESSIM: MSWOPTS 1 1 1") != std::string::npos;
        const Case y = importEclipseDeck(out);
        exportEclipseDeck(c, out);
        std::ifstream g(out); std::stringstream s2; s2 << g.rdbuf();
        const bool noLine = s2.str().find("MSWOPTS") == std::string::npos;
        std::remove(out.c_str());
        check(line && y.mswOpt.xflow && noLine, "C45: MSWOPTS third value written / read; defaults not written");
    }
    auto runAll = [](const Case& cc, std::vector<std::pair<double, std::string>>* blobs) {
        auto s = std::make_unique<Simulator>(); s->log = [](const std::string&) {}; s->setup(cc); s->initialize();
        while (s->runReportStep()) if (blobs) blobs->push_back({s->summary().empty() ? 0.0 : s->summary().back().time, s->saveState()});
        return s;
    };
    const auto a = runAll(c, nullptr);
    Case cx = c; cx.mswOpt.xflow = true;
    std::vector<std::pair<double, std::string>> blobs;
    const auto b = runAll(cx, &blobs);
    const double ga = a->summary().back().fgpt, gb = b->summary().back().fgpt, ob = b->summary().back().fopt;
    std::printf("    FGPT default %.6g, with MSW cross-flow %.6g\n", ga, gb);
    check(ga > 0 && std::fabs(gb - ga) > 1e-4 * ga, "C45: MSW cross-flow option changes the run (connections re-inject segment fluid)");
    if (blobs.size() < 3) { check(false, "C45: report steps for the restart"); return; }
    const auto& mid = blobs[blobs.size() / 2];
    const bool v16 = mid.second.find("RESSIM-BO-STATE-16") != std::string::npos;
    Simulator r; r.log = [](const std::string&) {}; r.setup(cx); r.initialize(); r.loadState(mid.second);
    while (r.runReportStep()) {}
    const double gr = r.summary().back().fgpt, orr = r.summary().back().fopt;
    std::printf("    restart at %.4g d: FGPT %.8g / %.8g, FOPT %.8g / %.8g (state v16 %s)\n", mid.first, gb, gr, ob, orr, v16 ? "yes" : "no");
    check(v16 && std::fabs(gr - gb) <= 1e-6 * (1 + gb) && std::fabs(orr - ob) <= 1e-6 * (1 + ob),
          "C45: restart (state v16, segment mixtures) reproduces the continuous MSW cross-flow run");
}
