// C49 tests: converged (Picard-coupled) multisegment connection heads (MSWOPTS 4th token), the cured cross-flow
// flip-flop of MSW-2D-VERT-02, restart with converged heads, FAULTMECH item 10 (ramp width of the two-way factor)
#include <filesystem>
#include <fstream>

static std::string c49Read(const std::string& p) { std::ifstream f(p, std::ios::binary); std::stringstream ss; ss << f.rdbuf(); return ss.str(); }

static void testC49MswOptsRoundTrip() {
    std::printf("== C49 MSWOPTS fourth token (connection heads): parse, export, defaults not written\n");
    const std::string deck = findDeck("decks/opm_msw/MSW-2D-VERT-02.DATA");
    if (deck.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(deck);
    check(c.mswOpt.heads == -1, "C49: connection heads default (-1 = converged)");
    const std::string out = "c49_msw_rt.DATA";
    for (int h : {0, 1}) {
        Case x = c; x.mswOpt.heads = h;
        exportEclipseDeck(x, out);
        const std::string t = c49Read(out);
        const Case y = importEclipseDeck(out);
        const std::string want = std::string("-- RESSIM: MSWOPTS 1 1 0 ") + char('0' + h);
        check(t.find(want) != std::string::npos && y.mswOpt.heads == h, (std::string("C49: MSWOPTS heads ") + char('0' + h) + " written and read back").c_str());
    }
    {
        Case x = c; x.mswOpt.xflow = true;   // C45 three-token form unchanged when heads is defaulted
        exportEclipseDeck(x, out);
        const std::string t = c49Read(out);
        check(t.find("-- RESSIM: MSWOPTS 1 1 1\n") != std::string::npos, "C49: cross-flow alone still writes the three-token form");
    }
    exportEclipseDeck(c, out);
    check(c49Read(out).find("MSWOPTS") == std::string::npos, "C49: default MSW options are not written");
    std::remove(out.c_str());
}

static void testC49ConvergedHeads() {
    std::printf("== C49 converged multisegment heads: MSW-2D-VERT-02 with cross-flow no longer flip-flops; restart\n");
    const std::string deck = findDeck("decks/opm_msw/MSW-2D-VERT-02.DATA");
    if (deck.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(deck);
    auto runAll = [](const Case& cc, std::vector<std::pair<double, std::string>>* blobs) {
        auto s = std::make_unique<Simulator>(); s->log = [](const std::string&) {}; s->setup(cc); s->initialize();
        while (s->runReportStep()) if (blobs) blobs->push_back({s->summary().empty() ? 0.0 : s->summary().back().time, s->saveState()});
        return s;
    };
    Case base = c; base.mswOpt.heads = 0;                            // C44 behaviour: no cross-flow, explicit heads
    Case expl = c; expl.mswOpt.xflow = true; expl.mswOpt.heads = 0;  // C45: cross-flow with explicit heads (flip-flop)
    Case conv = c; conv.mswOpt.xflow = true;                         // C49 default heads = converged
    const auto a = runAll(base, nullptr);
    const auto e = runAll(expl, nullptr);
    std::vector<std::pair<double, std::string>> blobs;
    const auto v = runAll(conv, &blobs);
    const double ga = a->summary().back().fgpt, ge = e->summary().back().fgpt, gv = v->summary().back().fgpt;
    const double oa = a->summary().back().fopt, oe = e->summary().back().fopt, ov = v->summary().back().fopt;
    std::printf("    FGPT no cross-flow %.6g | cross-flow explicit %.6g (%+.2f %%) | cross-flow converged %.6g (%+.2f %%)\n", ga, ge,
                100 * (ge - ga) / ga, gv, 100 * (gv - ga) / ga);
    std::printf("    FOPT %.6g | %.6g | %.6g; Picard: %d steps, %d re-solves, %d not converged\n", oa, oe, ov, v->mswPicardSteps_, v->mswPicardIts_, v->mswPicardUnconv_);
    check(v->mswPicardSteps_ > 0 && v->mswPicardIts_ > 0, "C49: the converged-head iterations run (default)");
    check(e->mswPicardSteps_ == 0, "C49: heads 0 keeps the explicit start-of-step heads");
    check(std::fabs(gv - ga) < 0.5 * std::fabs(ge - ga), "C49: converged heads remove most of the cross-flow FGPT deviation of the explicit heads");
    check(std::fabs(ov - oa) <= 1e-3 * oa && std::fabs(oe - oa) > 1e-2 * oa, "C49: the flip-flop's extra oil (FOPT +4 %) is gone with converged heads");
    if (blobs.size() < 3) { check(false, "C49: report steps for the restart"); return; }
    const auto& mid = blobs[blobs.size() / 2];
    Simulator r; r.log = [](const std::string&) {}; r.setup(conv); r.initialize(); r.loadState(mid.second);
    while (r.runReportStep()) {}
    const double gr = r.summary().back().fgpt;
    std::printf("    restart at %.4g d: FGPT %.8g / %.8g\n", mid.first, gv, gr);
    check(std::fabs(gr - gv) <= 1e-6 * (1 + gv), "C49: restart reproduces the continuous run with converged heads");
}

static void testC49FaultRamp() {
    std::printf("== C49 FAULTMECH item 10 (ramp width of the two-way permeability factor): parse, export, default\n");
    const std::string path = findExample("CO2_GEOMECH.DATA");
    if (path.empty()) { std::printf("    (examples/CO2_GEOMECH.DATA not found - skipped)\n"); return; }
    const std::string text = c32ReadText(path);
    const std::string a = "-- RESSIM: FAULTMECH F1 0.6 0 0 100 5 1.0 -1.0 5\n";
    const size_t q = text.find(a);
    check(q != std::string::npos, "C49: CO2_GEOMECH layout (FAULTMECH line)");
    if (q == std::string::npos) return;
    try {
        const Case d = importEclipseDeck(path, nullptr, nullptr);
        check(!d.faultMech.faults.empty() && std::fabs(d.faultMech.faults[0].permRamp - 0.05) < 1e-12, "C49: ramp defaults to 0.05 sigma'n");
        const std::string out = "c49_ramp_out.DATA";
        exportEclipseDeck(d, out);
        check(c49Read(out).find("-- RESSIM: FAULTMECH F1 0.6 0 0 100 5 1 -1 5\n") != std::string::npos, "C49: default ramp is not written (C48 line unchanged)");
        std::string t2 = text; t2.replace(q, a.size(), "-- RESSIM: FAULTMECH F1 0.6 0 0 100 5 1.0 -1.0 5 0.5\n");
        { std::ofstream(out, std::ios::binary) << t2; }
        const Case e = importEclipseDeck(out, nullptr, nullptr);
        check(!e.faultMech.faults.empty() && std::fabs(e.faultMech.faults[0].permRamp - 0.5) < 1e-12, "C49: item 10 read (0.5)");
        exportEclipseDeck(e, out);
        const Case f = importEclipseDeck(out, nullptr, nullptr);
        check(!f.faultMech.faults.empty() && std::fabs(f.faultMech.faults[0].permRamp - 0.5) < 1e-12, "C49: ramp survives the round trip");
        std::remove(out.c_str());
    } catch (const std::exception& ex) { check(false, (std::string("C49: FAULTMECH ramp deck: ") + ex.what()).c_str()); }
}

// C49: STOP-status wells report their zero-inflow BHP as WBHP (OPM / ECLIPSE, opm-tests SPE1CASE1_SHUT_CONTROLS);
// SHUT wells report 0
static void testC49StoppedWellBhp() {
    std::printf("== C49 STOP-status wells: WBHP = zero-inflow BHP, SHUT wells 0\n");
    const std::string path = findDeck("examples/SPE1.DATA");
    std::string text = c49Read(path);
    const size_t ts = text.find("\nTSTEP");
    check(ts != std::string::npos, "C49: SPE1 example layout (TSTEP)");
    if (ts == std::string::npos) return;
    text = text.substr(0, ts) + "\nWELSPECS\n 'S1' 'PROD' 5 5 8400 OIL /\n 'X1' 'PROD' 7 7 8400 OIL /\n/\n"
           "COMPDAT\n 'S1' 5 5 3 3 OPEN 1* 1* 0.5 1* 0 1* Z /\n 'X1' 7 7 3 3 OPEN 1* 1* 0.5 1* 0 1* Z /\n/\n"
           "WCONPROD\n 'S1' STOP ORAT 1000 1* 1* 1* 1* 1000 /\n 'X1' SHUT ORAT 1000 1* 1* 1* 1* 1000 /\n/\n"
           "TSTEP\n 31 28 /\nEND\n";
    const std::string out = "c49_stop.DATA";
    { std::ofstream(out, std::ios::binary) << text; }
    try {
        const Case c = importEclipseDeck(out, nullptr, nullptr);
        const api::RunOutput r = api::runCase(c);
        check(r.ok, "C49: SPE1 + STOP / SHUT wells runs");
        if (r.ok) {
            const ecl::Summary s = api::toSummary(r, c);
            const int is = s.find("WBHP:S1"), ix = s.find("WBHP:X1"), io = s.find("WOPR:S1"), ip = s.find("FPR");
            check(is >= 0 && ix >= 0 && io >= 0 && ip >= 0, "C49: WBHP:S1 / WBHP:X1 / WOPR:S1 / FPR present");
            if (is >= 0 && ix >= 0 && io >= 0 && ip >= 0) {
                const double t = s.time.back();
                const double bs = s.at(is, t), bx = s.at(ix, t), qo = s.at(io, t), fpr = s.at(ip, t);
                std::printf("   WBHP:S1 %.2f  WBHP:X1 %.2f  WOPR:S1 %.3g  FPR %.2f\n", bs, bx, qo, fpr);
                check(bs > 0.9 * fpr && bs < 1.1 * fpr, "C49: STOP well reports a BHP near the reservoir pressure");
                check(bx == 0.0, "C49: SHUT well reports WBHP 0");
                check(qo == 0.0, "C49: STOP well does not flow");
            }
        }
    } catch (const std::exception& ex) { check(false, (std::string("C49: STOP deck: ") + ex.what()).c_str()); }
    std::remove(out.c_str());
}

// C49: compositional / CO2STORE runs through api::runCase (regression runner, studies) keep the per-well rates and BHP
// (they were dropped: WBHP:INJ = 0 in the OPM CO2STORE regression)
static void testC49CompWellSummary() {
    std::printf("== C49 CO2STORE via api::runCase: per-well WBHP / WGIR present\n");
    try {
        const Case c = importEclipseDeck(findDeck("examples/CO2STORE.DATA"), nullptr, nullptr);
        const api::RunOutput r = api::runCase(c);
        check(r.ok, "C49: CO2STORE example runs through api::runCase");
        if (!r.ok) return;
        check(!r.wellNames.empty() && !r.summary.empty() && r.summary.back().wells.size() == r.wellNames.size(), "C49: per-well data and names kept");
        const ecl::Summary s = api::toSummary(r, c);
        const int ib = s.find("WBHP:INJ"), ig = s.find("WGIR:INJ"), fg = s.find("FGIR");
        check(ib >= 0 && ig >= 0 && fg >= 0, "C49: WBHP:INJ / WGIR:INJ / FGIR present");
        if (ib >= 0 && ig >= 0 && fg >= 0) {
            // the example ends with a shut-in period: check the injection steps
            const auto B = s.series(ib), G = s.series(ig), F = s.series(fg);
            int nInj = 0, bhpOk = 0; double dmax = 0, fmax = 0;
            for (size_t k = 0; k < F.size() && k < B.size() && k < G.size(); ++k) {
                fmax = std::max(fmax, F[k]); dmax = std::max(dmax, std::fabs(G[k] - F[k]));
                if (F[k] > 0) { ++nInj; if (B[k] > 1.0) ++bhpOk; }
            }
            std::printf("   injection steps %d, with WBHP > 0: %d, max |WGIR-FGIR| %.3g of max FGIR %.4g\n", nInj, bhpOk, dmax, fmax);
            check(nInj > 0 && bhpOk == nInj, "C49: injector BHP reported (> 0) at every injection step");
            check(dmax <= 1e-3 * std::max(1.0, fmax), "C49: WGIR:INJ = FGIR (single injector)");
        }
    } catch (const std::exception& ex) { check(false, (std::string("C49: CO2STORE api run: ") + ex.what()).c_str()); }
}
