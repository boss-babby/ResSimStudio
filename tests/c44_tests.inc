// C44 tests: WPWE against OPM reference summaries, WELSPECS item 9 (STOP), connection-level WELOPEN on one layer of a
// multi-layer COMPDAT record
#include "ressim/api.h"

// WPWE events summed per reference report step: "WPWEn:WELL" -> {report index -> sum}
static std::map<std::string, std::map<int, double>> c44Events(const ecl::Summary& S, const std::vector<double>& rt) {
    std::map<std::string, std::map<int, double>> out;
    for (size_t v = 0; v < S.vars.size(); ++v) {
        const auto& V = S.vars[v];
        if (V.keyword.rfind("WPWE", 0) != 0 || V.keyword == "WPWEM") continue;
        for (size_t m = 0; m < S.time.size() && m < S.values[v].size(); ++m) {
            if (S.values[v][m] == 0) continue;
            int r = 0; while (r < int(rt.size()) && S.time[m] > rt[size_t(r)] + 1e-6) ++r;
            out[V.keyword + ":" + V.wgname][r] += S.values[v][m];
        }
    }
    return out;
}

static void testC44WpweReference() {
    std::printf("== C44 WPWE event timelines against OPM Flow (opm-tests wpwe/)\n");
    const char* decks[] = {"WPWE-WORKOVER-SHUT-STOP", "WPWE-COMPLUMP-CONNECTION-COUNT", "WPWE-DECK-CONNECTION-CHANGES",
                           "WPWE-DECK-STATUS-CONVERSION", "WPWE-MIXED-CLOSURE-CAUSES"};
    for (const char* d : decks) {
        const std::string deck = findDeck(std::string("decks/opm_wpwe/") + d + ".DATA");
        if (deck.empty()) { std::printf("  (deck %s not found - skipped)\n", d); continue; }
        const std::string ref = deck.substr(0, deck.find_last_of("/\\") + 1) + "opm-simulation-reference/flow/" + d + ".SMSPEC";
        ecl::Summary R; std::string err;
        api::Session S;
        const bool ok = ecl::readSummary(ref, R, &err) && S.open(deck) && S.run();
        int series = 0, diff = 0;
        if (ok) {
            std::vector<double> rt; for (int i : R.reportIndex) if (i >= 0 && size_t(i) < R.time.size()) rt.push_back(R.time[size_t(i)]);
            const auto er = c44Events(R, rt), eo = c44Events(S.summary(), rt);
            std::set<std::string> keys; for (auto& kv : er) keys.insert(kv.first); for (auto& kv : eo) keys.insert(kv.first);
            for (const auto& k : keys) {
                ++series;
                const auto a = er.find(k), b = eo.find(k);
                const bool same = a != er.end() && b != eo.end() && a->second == b->second;
                if (!same) { ++diff; std::printf("    differs: %s\n", k.c_str()); }
            }
        }
        std::printf("    %s: %d event series, %d differ\n", d, series, diff);
        check(ok && series > 0 && diff == 0, (std::string("C44: WPWE timeline identical to OPM Flow - ") + d).c_str());
    }
}

static void testC44ShutInStop() {
    std::printf("== C44 WELSPECS item 9 (automatic shut-in instruction)\n");
    auto run = [](const std::string& tail, const std::string& cecon) {
        Case c = smallCase();
        for (auto& w : c.wells) if (w.name == "PROD") w.specsTail = tail;
        ExtraKeyword k; k.name = "CECON"; k.section = DeckSection::Schedule; k.step = 0; k.body = cecon; c.extraKeywords.push_back(k);
        Simulator s; s.setup(c); s.initialize(); while (s.runReportStep()) {}
        int pw = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == "PROD") pw = int(w);
        std::array<int, 8> sum{}; for (const auto& p : s.summary()) for (int e = 0; e < 8; ++e) sum[size_t(e)] += p.wells[size_t(pw)].pwe[e];
        std::printf("    tail '%s', %s -> WPWE4 %d, WPWE7 %d\n", tail.c_str(), cecon.substr(0, cecon.size() - 1).c_str(), sum[4], sum[7]);
        return sum;
    };
    const auto a = run("1* 1* 'STOP'", " 'PROD' 0 0 1 1 1e-6 1* 1* 'WELL' 'NO' /\n");
    const auto b = run("1* 1* 'SHUT'", " 'PROD' 0 0 1 1 1e-6 1* 1* 'WELL' 'NO' /\n");
    const auto c = run("1* 1* 'STOP'", " 'PROD' 0 0 0 0 1e-6 1* 1* 'CON' 'NO' /\n");
    check(a[4] == 1 && a[7] == 0 && b[4] == 0 && b[7] == 1, "C44: WELL workover follows WELSPECS item 9 - STOP -> WPWE4, SHUT -> WPWE7");
    check(c[7] == 1 && c[4] == 0, "C44: a well left without open connections is shut (WPWE7) even with item 9 = STOP (OPM)");
}

static void testC44WelopenLayer() {
    std::printf("== C44 WELOPEN on one layer of a multi-layer COMPDAT record\n");
    const std::string deck = findDeck("decks/opm_wpwe/WPWE-MIXED-CLOSURE-CAUSES.DATA");
    if (deck.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    api::Session S;
    if (!S.open(deck)) { check(false, "C44: deck opens"); return; }
    int n = 0, open = 0;
    for (const auto& w : S.mutableCase().wells) if (w.name == "C-1H") for (const auto& cp : w.completions) { ++n; open += cp.open ? 1 : 0; }
    std::printf("    C-1H: %d connection records, %d open at the start (COMPDAT 1-5, WELOPEN SHUT K = 1, 2, 3)\n", n, open);
    check(n == 5 && open == 2, "C44: WELOPEN SHUT on single layers splits the COMPDAT 1-5 record; layers 4 and 5 stay open");
}

// C44 restart state v15: the economic connection closures (WTEST 'E' re-opening), the WPWE baseline and the injector
// recirculation totals are saved. A restart taken while the well is shut by CON on all connections must re-open them at the
// next retest exactly like the continuous run.
static void testC44RestartEcon() {
    std::printf("== C44 restart state v15 (economic connection closures, WPWE baseline)\n");
    Case c = smallCase();
    { ExtraKeyword k; k.name = "CECON"; k.section = DeckSection::Schedule; k.step = 0; k.body = " 'PROD' 0 0 0 0 1e-6 1* 1* 'CON' 'NO' /\n"; c.extraKeywords.push_back(k); }
    { ExtraKeyword t; t.name = "WTEST"; t.section = DeckSection::Schedule; t.step = 0; t.body = " 'PROD' 200 'E' /\n"; c.extraKeywords.push_back(t); }
    int pw = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == "PROD") pw = int(w);
    Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
    std::vector<std::pair<double, std::string>> blobs;
    while (a.runReportStep()) blobs.push_back({a.summary().empty() ? 0.0 : a.summary().back().time, a.saveState()});
    double tClose = -1, tOpen = -1;
    for (const SummaryPoint& p : a.summary()) {
        if (tClose < 0 && p.wells[size_t(pw)].pwe[2] > 0) tClose = p.time;
        if (tClose >= 0 && tOpen < 0 && p.wells[size_t(pw)].pwe[1] > 0) tOpen = p.time;
    }
    const std::string* blob = nullptr; double tB = -1;
    for (const auto& [t, b] : blobs) if (t >= tClose && t < tOpen) { blob = &b; tB = t; break; }
    if (!blob) {
        std::printf("    closed %.4g d, re-opened %.4g d, report times:", tClose, tOpen);
        for (const auto& bt : blobs) std::printf(" %.4g", bt.first);
        std::printf("\n");
        check(false, "C44: restart point between the closure and the retest"); return;
    }
    Simulator b; b.log = [](const std::string&) {}; b.setup(c); b.initialize(); b.loadState(*blob);
    while (b.runReportStep()) {}
    auto tail = [&](const Simulator& s, int& w1, int& w2, double& fopt) {
        w1 = w2 = 0; fopt = s.summary().empty() ? 0 : s.summary().back().fopt;
        for (const SummaryPoint& p : s.summary()) if (p.time > tB + 1e-9) { w1 += p.wells[size_t(pw)].pwe[1]; w2 += p.wells[size_t(pw)].pwe[2]; }
    };
    int a1, a2, b1, b2; double fa, fb; tail(a, a1, a2, fa); tail(b, b1, b2, fb);
    std::printf("    closed at %.4g d, restart at %.4g d, retest %.4g d: after restart WPWE1 %d / %d, WPWE2 %d / %d, FOPT %.6g / %.6g\n",
                tClose, tB, tOpen, a1, b1, a2, b2, fa, fb);
    check(b1 == a1 && b2 == a2 && a1 > 0 && std::fabs(fa - fb) <= 1e-6 * (1 + std::fabs(fa)),
          "C44: restart while shut by CON - WTEST 'E' re-opens the connections and events / FOPT match the continuous run");
}

// C44 EHYSTR models 1 / 3: the wetting phases (krw, krog) follow the imbibition curve (OPM EclHysteresisTwoPhaseLaw), models
// 0 / 2 keep them on drainage. With a larger imbibition critical water saturation (ISWCR) water is less mobile under model 1.
static void testC44HysteresisWetting() {
    std::printf("== C44 EHYSTR 1 / 3 wetting-phase imbibition\n");
    auto run = [](int model, double& fwpt, double& fopt, std::string& log) {
        Case c = makeExample("qfs");
        c.props.hyst.on = true; c.props.hyst.model = model; c.props.hyst.what = "KR";
        c.grid.endscale = true;
        const size_t n = size_t(c.grid.nx * c.grid.ny * c.grid.nz);
        c.grid.iswcr.assign(n, 0.4);
        Simulator s; s.log = [&](const std::string& m) { log += m + "\n"; }; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        fwpt = s.summary().back().fwpt; fopt = s.summary().back().fopt;
    };
    double w0, o0, w1, o1; std::string l0, l1;
    run(0, w0, o0, l0); run(1, w1, o1, l1);
    std::printf("    model 0: FWPT %.6g FOPT %.6g;  model 1: FWPT %.6g FOPT %.6g\n", w0, o0, w1, o1);
    check(l1.find("wetting phases (krw, krog) on the imbibition curve") != std::string::npos && l0.find("wetting phases on drainage") != std::string::npos,
          "C44: EHYSTR 1 reports imbibition wetting curves, model 0 drainage");
    check(w1 < w0 * 0.999, "C44: EHYSTR 1 - larger ISWCR on the imbibition krw curve reduces water production");
}

// C44 injector cross-flow: fluid back-produced at injector connections mixes in the wellbore and leaves through the injecting
// connections (no surface production). OPM opm-tests model2/GUIDERAT_POTENTIAL_FALLBACK: INJ2 back-flowed gas used to be
// reported as field gas production and starved the producers' group target.
static void testC44InjectorRecirculation() {
    std::printf("== C44 injector cross-flow recirculation (GUIDERAT_POTENTIAL_FALLBACK)\n");
    const std::string deck = findDeck("decks/opm_model2/GUIDERAT_POTENTIAL_FALLBACK.DATA");
    if (deck.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    const std::string ref = deck.substr(0, deck.rfind('/') + 1) + "opm-simulation-reference/flow/GUIDERAT_POTENTIAL_FALLBACK.SMSPEC";
    ecl::Summary R; std::string err;
    if (!ecl::readSummary(ref, R, &err)) { check(false, "C44: reference summary readable"); return; }
    api::Session S;
    if (!S.open(deck) || !S.run()) { check(false, "C44: GUIDERAT_POTENTIAL_FALLBACK runs"); return; }
    const ecl::Summary& O = S.summary();
    auto find = [](const ecl::Summary& X, const char* kw, const char* wg) {
        for (size_t v = 0; v < X.vars.size(); ++v) if (X.vars[v].keyword == kw && (!*wg || X.vars[v].wgname == wg)) return int(v);
        return -1;
    };
    auto at = [](const ecl::Summary& X, int v, double t) {
        size_t i = 0; while (i < X.time.size() && X.time[i] < t - 1e-9) ++i;
        return X.values[size_t(v)][std::min(i, X.time.size() - 1)];
    };
    const int rg = find(R, "FGPR", ""), og = find(O, "FGPR", ""), rb = find(R, "WBHP", "INJ2"), ob = find(O, "WBHP", "INJ2");
    const int og2 = find(O, "WGPR", "INJ2");
    if (rg < 0 || og < 0 || rb < 0 || ob < 0) { check(false, "C44: FGPR / WBHP:INJ2 present"); return; }
    double eg = 0, eb = 0, backGas = 0;
    for (size_t i = 1; i < R.time.size(); ++i) {
        const double t = R.time[i];
        const double rgv = R.values[size_t(rg)][i], rbv = R.values[size_t(rb)][i];
        eg = std::max(eg, std::fabs(double(at(O, og, t)) - rgv) / std::max(1.0, std::fabs(rgv)));
        eb = std::max(eb, std::fabs(double(at(O, ob, t)) - rbv) / std::max(1.0, std::fabs(rbv)));
        if (og2 >= 0) backGas = std::max(backGas, double(std::fabs(at(O, og2, t))));
    }
    std::printf("    max rel. difference vs OPM: FGPR %.3g %%, WBHP:INJ2 %.3g %%; INJ2 reported gas production %.4g\n", 100 * eg, 100 * eb, backGas);
    check(eg < 0.02 && eb < 0.01 && backGas < 1.0, "C44: injector back-flow recirculates - FGPR and WBHP:INJ2 match OPM Flow, no surface production from INJ2");
}

// C44 WPIMULT per connection and step (ECLIPSE / OPM): applied to the connections existing at its step, cumulative, not
// inherited by connections defined later; FBHPDEF default BHP limits. Both survive the deck round trip.
static void testC44WpimultFbhpdef() {
    std::printf("== C44 WPIMULT per connection / step, FBHPDEF\n");
    std::string d = eclipseDeckToString(smallCase());
    const size_t p0 = d.find("WCONPROD"), p1 = d.find("END", p0);
    if (p0 == std::string::npos || p1 == std::string::npos) { check(false, "C44: deck layout"); return; }
    d = d.substr(0, p0) +
        "FBHPDEF\n 50 500 /\n\nWPIMULT\n 'PROD' 0.5 /\n/\n\nWCONPROD\n 'PROD' OPEN ORAT 10000 /\n/\n\nWCONINJE\n 'INJ' WATER OPEN RATE 12000 /\n/\n\n"
        "TSTEP\n 2*91.3125 /\n\nWPIMULT\n 'PROD' 0.2 1* 1* 2 /\n/\n\nTSTEP\n 2*91.3125 /\n\n"
        "COMPDAT\n 'PROD' 9 9 1 1 OPEN 1* 1* 0.5 1* 0 1* Z /\n/\n\nTSTEP\n 2*91.3125 /\n\nEND\n";
    auto eff = [](const Case& c, int i, int j, int k, int s) {
        for (const auto& w : c.wells) if (w.name == "PROD")
            for (const auto& q : w.completions)
                if (q.open && q.i == i && q.j == j && k >= q.k1 && k <= q.k2 && s >= std::max(0, q.fromStep) && (q.untilStep < 0 || s < q.untilStep)) return q.mult * w.piMult;
        return -1.0;
    };
    auto bhpOf = [](const Case& c, const char* well) {
        for (const auto& wc : c.schedule[0].controls) if (wc.well == well) return wc.bhp;
        return -1.0;
    };
    try {
        const Case a = eclipseDeckFromString(d, "");
        const std::string out = eclipseDeckToString(a);
        const Case b = eclipseDeckFromString(out, "");
        for (const Case* c : {&a, &b}) {
            const double k1s0 = eff(*c, 9, 9, 0, 0), k2s1 = eff(*c, 9, 9, 1, 1), k2s3 = eff(*c, 9, 9, 1, 3), k1s3 = eff(*c, 9, 9, 0, 3), nw = eff(*c, 8, 8, 0, 5);
            std::printf("    %s: K1 %.3g (step 0) %.3g (step 3), K2 %.3g (step 1) %.3g (step 3), new connection %.3g; BHP limits PROD %.4g INJ %.4g\n",
                        c == &a ? "imported" : "round trip", k1s0, k1s3, k2s1, k2s3, nw, bhpOf(*c, "PROD"), bhpOf(*c, "INJ"));
            check(std::fabs(k1s0 - 0.5) < 1e-12 && std::fabs(k1s3 - 0.5) < 1e-12 && std::fabs(k2s1 - 0.5) < 1e-12 && std::fabs(k2s3 - 0.1) < 1e-12 && std::fabs(nw - 1) < 1e-12,
                  c == &a ? "C44: WPIMULT - cumulative on existing connections from its step, earlier periods and later connections unaffected"
                          : "C44: WPIMULT factors per connection and step survive the deck round trip");
            check(std::fabs(bhpOf(*c, "PROD") - 50) < 1e-9 && std::fabs(bhpOf(*c, "INJ") - 500) < 1e-9,
                  c == &a ? "C44: FBHPDEF sets the defaulted producer / injector BHP limits" : "C44: FBHPDEF defaults survive the round trip");
        }
        check(out.find("FBHPDEF") != std::string::npos, "C44: FBHPDEF kept in the exported deck");
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); check(false, "C44: WPIMULT / FBHPDEF deck parses"); }
}

// C44 geomodel: fault cut on an intersection / 2-D line path (stick interpolation per depth level)
#include "ressim/geo_tools.h"
static void testC44FaultCut() {
    std::printf("== C44 fault cut on a plan path\n");
    using namespace ressim::geo;
    Polygon path; path.closed = false; path.pts = {{0, 500}, {600, 500}, {1000, 500}};   // two segments along y = 500
    Fault vert; vert.name = "V";          // vertical plane x = 400
    vert.sticks = {{{400, 0, 1000}, {400, 0, 2000}}, {{400, 1000, 1000}, {400, 1000, 2000}}};
    auto br = faultCutOnPath(vert, path, 900, 2100, 20);   // clipped to the stick depths 1000..2000
    bool ok = br.size() == 1 && br[0].size() == 21;
    for (const auto& q : ok ? br[0] : std::vector<std::array<double, 2>>{}) ok = ok && std::fabs(q[0] - 400) < 1e-9 && q[1] >= 1000 - 1e-9 && q[1] <= 2000 + 1e-9;
    check(ok, (std::string("vertical fault: one branch of 21 points at s = 400 (got ") + std::to_string(br.size()) + " branches)").c_str());
    Fault dip; dip.name = "D";            // dipping 60 deg-ish: x = 300 + 0.5 (z - 1000), crossing the segment joint at x = 600
    dip.sticks = {{{300, 0, 1000}, {800, 0, 2000}}, {{300, 1000, 1000}, {800, 1000, 2000}}};
    br = faultCutOnPath(dip, path, 1000, 2000, 10);
    double err = 0; size_t n = 0;
    for (const auto& b : br) for (const auto& q : b) { err = std::max(err, std::fabs(q[0] - (300 + 0.5 * (q[1] - 1000)))); ++n; }
    std::printf("  dipping fault: %zu branch(es), %zu points, max |s - s_exact| %.2e\n", br.size(), n, err);
    check(br.size() == 1 && n >= 11 && err < 1e-6, "dipping fault: one continuous branch along s = 300 + 0.5 (z - 1000) across the path joint");
    Fault off; off.name = "O"; off.sticks = {{{2000, 0, 1000}, {2000, 0, 2000}}, {{2000, 1000, 1000}, {2000, 1000, 2000}}};
    check(faultCutOnPath(off, path, 1000, 2000, 10).empty(), "fault beyond the path end: no cut");
}
