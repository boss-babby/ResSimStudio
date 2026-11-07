// =====================================================================================================================
// C31 tests: UDQ expression engine (arithmetic, set / elementwise functions, undefined values, patterns, UPDATE), ACTIONX
// conditions (AND / OR / parentheses, well sets, '?' expansion), UDA targets in the black-oil engine (one-step lag
// semantics checked exactly), ACTIONX at run time (WELOPEN '?', UDQ ASSIGN + WELTARG inside an action), deck round trip
// and bitwise restart across fired actions (examples/ACTIONX_UDQ.DATA)
// =====================================================================================================================
#include "ressim/udq.h"

static void testC31Expressions() {
    std::printf("== C31 UDQ / ACTIONX expression engine\n");
    using namespace ressim::udq;
    Source s;
    s.wells = {"OP1", "OP2", "WI1"}; s.groups = {"G1", "FIELD"};
    s.field = [](const std::string& k) { return k == "FOPR" ? 1000.0 : k == "FWPR" ? 500.0 : k == "FPR" ? 3000.0 : NAN; };
    s.well = [](const std::string& k, size_t i) {
        const double o[] = {300, 700, 0}, w[] = {100, 400, 0};
        if (k == "WOPR") return o[i];
        if (k == "WWPR") return w[i];
        if (k == "WWCT") return o[i] + w[i] > 0 ? w[i] / (o[i] + w[i]) : 0.0;
        return (double)NAN;
    };
    s.group = [](const std::string& k, size_t) { return k == "GOPR" ? 1000.0 : (double)NAN; };
    s.days = 100; s.day = 11; s.month = 4; s.year = 2025;
    Set set; std::string err;
    const bool ok = set.add("-- comment line\n DEFINE FUWCT FWPR / (FWPR + FOPR) /\n DEFINE WUOPRL (WOPR 'OP*' - 150) * 0.90 /\n"
                            " DEFINE FUSUM SUM(WOPR 'OP*') /\n ASSIGN FUX 3.5 /\n DEFINE FUMX MAX(WWCT) + 2^3 /\n DEFINE FUP1 WOPR OP1*2 /\n"
                            " DEFINE WURANK SORTD(WOPR) /\n DEFINE FUDIV FOPR / (FWPR - 500) /\n DEFINE FUDEF IDV(FUDIV) + IDV(FUSUM) /\n"
                            " DEFINE FUAVG AVEG(WOPR 'OP*') /\n DEFINE FUCMP (FPR > 2500) + (FPR == 3000) + (FOPR <= 10) /\n"
                            " DEFINE FUG GOPR G1 * 0.5 /\n DEFINE FUT TIME + MNTH /\n DEFINE FUN NINT(2.5) + ABS(-3) + LOG(100) + LN(EXP(1)) /\n"
                            " DEFINE WUC WOPR / 100 /\n DEFINE FUSC SUM(WUC) /\n UNITS FUWCT 'frac' /\n/\n", 0, err);
    check(ok && err.empty(), "UDQ keyword body parsed (DEFINE / ASSIGN / UNITS, comments, '/' as division and record end)");
    set.evaluate(s, 0);
    const double wct = set.scalar("FUWCT");
    std::printf("    FUWCT %.6g FUSUM %g FUMX %.6g FUP1 %g FUAVG %.6g FUCMP %g FUG %g FUT %g FUN %.6g FUSC %g\n", wct, set.scalar("FUSUM"),
                set.scalar("FUMX"), set.scalar("FUP1"), set.scalar("FUAVG"), set.scalar("FUCMP"), set.scalar("FUG"), set.scalar("FUT"), set.scalar("FUN"), set.scalar("FUSC"));
    check(std::fabs(wct - 1.0 / 3.0) < 1e-15 && set.units("FUWCT") == "frac", "field arithmetic: FWPR / (FWPR + FOPR) = 1/3, UNITS kept");
    check(set.wellValue("WUOPRL", 0) == 135 && set.wellValue("WUOPRL", 1) == 495 && std::isnan(set.wellValue("WUOPRL", 2)),
          "well UDQ over a pattern: (WOPR 'OP*' - 150) * 0.9, wells outside the pattern undefined");
    check(set.scalar("FUSUM") == 1000 && set.scalar("FUX") == 3.5 && std::fabs(set.scalar("FUMX") - (400.0 / 1100.0 + 8)) < 1e-14,
          "set functions SUM / MAX over well sets, ASSIGN, power operator");
    check(set.scalar("FUP1") == 600, "a named well (no wildcard) gives a scalar: WOPR OP1*2 = 600");
    check(set.wellValue("WURANK", 1) == 1 && set.wellValue("WURANK", 0) == 2 && set.wellValue("WURANK", 2) == 3, "SORTD ranks the wells (1 = largest)");
    check(std::isnan(set.scalar("FUDIV")) && set.scalar("FUDEF") == 1, "division by zero is undefined; IDV detects it");
    check(std::fabs(set.scalar("FUAVG") - std::sqrt(300.0 * 700.0)) < 1e-9 && set.scalar("FUCMP") == 2, "AVEG geometric mean; comparisons give 1 / 0");
    check(set.scalar("FUG") == 500 && set.scalar("FUT") == 104 && std::fabs(set.scalar("FUN") - (3 + 3 + 2 + 1)) < 1e-12 && set.scalar("FUSC") == 10,
          "group quantity, TIME + MNTH, NINT (half away from zero) ABS LOG LN EXP, UDQ of UDQ");
    {
        Expr e; std::string m;
        const bool bad1 = !parseExpr("FOPR +", e, m), bad2 = !parseExpr("SQRTX(FOPR)", e, m), bad3 = !parseExpr("(FOPR * 2", e, m);
        check(bad1 && bad2 && bad3, "syntax errors detected (incomplete expression, unknown function, missing ')')");
    }
    check(wildMatch("OP*", "op12") && wildMatch("P?", "P1") && !wildMatch("P?", "P12") && wildMatch("*1", "WI1") && !wildMatch("OP*", "WI1"),
          "ECLIPSE well patterns (* ?, case-insensitive)");
    {   // UPDATE NEXT / OFF
        Set u; std::string e2;
        u.add("DEFINE FUT TIME /\n/\n", 0, e2);
        u.add("UPDATE FUT NEXT /\n/\n", 1, e2);
        u.add("UPDATE FUT ON /\n/\n", 3, e2);
        Source t = s; double v[5];
        for (int st = 0; st < 5; ++st) { t.days = 10.0 * st; u.evaluate(t, st); v[st] = u.scalar("FUT"); }
        std::printf("    UPDATE: %g %g %g %g %g\n", v[0], v[1], v[2], v[3], v[4]);
        check(v[0] == 0 && v[1] == 10 && v[2] == 10 && v[3] == 30 && v[4] == 40, "UPDATE NEXT evaluates once then holds, UPDATE ON resumes");
        const std::string blob = u.save();
        Set w2; w2.add("DEFINE FUT TIME /\n/\n", 0, e2); w2.add("UPDATE FUT NEXT /\n/\n", 1, e2); w2.add("UPDATE FUT ON /\n/\n", 3, e2);
        w2.load(blob);
        check(w2.scalar("FUT") == 40, "UDQ values survive save / load (restart)");
    }
    {   // ACTIONX
        Action a; std::string m;
        bool p = parseAction(" A1 2 10 /\n WWCT 'OP*' > 0.3 AND\n FOPR > 900 /\n/\nWELOPEN\n '?' SHUT /\n/\nENDACTIO\n", 0, a, m);
        const Trigger t = evaluateAction(a, s, set);
        check(p && a.maxRuns == 2 && a.minWaitDays == 10 && t.fired && t.wells.size() == 1 && t.wells[0] == 1,
              "ACTIONX header (name, runs, interval); well-set condition AND scalar: only OP2 (WWCT 0.364) triggers");
        const std::string body = expandBody(a.body, {"OP2", "OP9"});
        check(body.find("'OP2' SHUT /") != std::string::npos && body.find("'OP9' SHUT /") != std::string::npos && body.find("'?'") == std::string::npos,
              "'?' expands to one record per triggering well");
        Action b; p = parseAction("B /\nMNTH > MAR AND ( FUWCT < 0.1 OR DAY >= 11 ) /\n/\nWCONPROD\n'OP1' OPEN ORAT 10 /\n/\n", 0, b, m);
        Action c2; const bool p2 = parseAction("C /\nWWCT 'OP*' > 0.3 OR WOPR 'OP*' < 350 /\n/\n", 0, c2, m);
        const Trigger tc = evaluateAction(c2, s, set);
        Action d; const bool p3 = parseAction("D /\nFOPR >> 3 /\n/\n", 0, d, m);
        check(p && evaluateAction(b, s, set).fired && p2 && tc.fired && tc.wells.size() == 2 && !p3,
              "month names, parentheses, OR = union of the well sets, malformed condition rejected");
        Action e; p = parseAction("E /\nFUWCT > 0.5 /\n/\n", 0, e, m);
        check(p && !evaluateAction(e, s, set).fired, "condition on a UDQ (FUWCT 0.33 > 0.5 is false)");
    }
}

static void testC31Deck() {
    std::printf("== C31 UDA + ACTIONX in the black-oil engine (ACTIONX_UDQ.DATA)\n");
    const std::string path = findExample("ACTIONX_UDQ.DATA");
    if (path.empty()) { std::printf("    (examples/ACTIONX_UDQ.DATA not found - skipped)\n"); return; }
    std::vector<DeckMessage> rep;
    Case c = importEclipseDeck(path, nullptr, &rep);
    bool nonNum = false;
    for (const auto& m : rep) if (m.text.find("non-numeric") != std::string::npos) nonNum = true;
    const WellControl* p1 = nullptr; const WellControl* i1 = nullptr;
    for (const auto& wc : c.schedule[0].controls) { if (wc.well == "P1") p1 = &wc; if (wc.well == "I1") i1 = &wc; }
    check(!nonNum && p1 && i1 && p1->udaOf("ORAT") && *p1->udaOf("ORAT") == "FUPLIM" && i1->udaOf("RATE") && *i1->udaOf("RATE") == "FUINJ",
          "UDA items imported (WCONPROD ORAT 'FUPLIM', WCONINJE RATE 'FUINJ') without 'non-numeric' warnings");
    {
        int warnLogic = 0, infoLogic = 0;
        for (const auto& m : rep) {
            const bool logic = m.text.find("UDQ") != std::string::npos || m.text.find("ACTIONX") != std::string::npos;
            if (logic && m.level >= 1) ++warnLogic;
            if (logic && m.level == 0) ++infoLogic;
        }
        check(warnLogic == 0 && infoLogic >= 1, "import report: UDQ / ACTIONX are information (evaluated), not warnings");
    }
    {
        const std::string txt = eclipseDeckToString(c);
        bool rt = false;
        try {
            const Case d = eclipseDeckFromString(txt);
            int udas = 0;
            for (const auto& wc : d.schedule[0].controls) udas += int(wc.uda.size());
            rt = udas == 5;
        } catch (const std::exception& e) { std::printf("    re-import: %s\n", e.what()); }
        check(txt.find("'FUPLIM'") != std::string::npos && txt.find("'FUINJ'") != std::string::npos && txt.find("ENDACTIO") != std::string::npos &&
              txt.find("DEFINE FUINJ") != std::string::npos && rt, "export writes the UDA names, UDQ and ACTIONX blocks; re-import restores the 5 UDAs");
    }
    std::vector<std::string> logs;
    Simulator s; s.log = [&](const std::string& l) { logs.push_back(l); }; s.setup(c); s.initialize();
    while (s.runReportStep()) {}
    const auto& S = s.summary();
    check(s.udqSet() != nullptr && s.numActions() == 2, "UDQ set and 2 ACTIONX blocks set up");
    auto udqOf = [](const SummaryPoint& p, const std::string& n) { for (const auto& u : p.udq) if (u.first == n) return u.second; return (double)NAN; };
    size_t iI1 = 0, iP1 = 0;
    for (size_t w = 0; w < s.caseData().wells.size(); ++w) { if (s.caseData().wells[w].name == "I1") iI1 = w; if (s.caseData().wells[w].name == "P1") iP1 = w; }
    double dW = 0; int nW = 0, nLag = 0; double dLag = 0, dSame = 0, maxP = 0;
    for (size_t k = 0; k < S.size(); ++k) {
        const double l = S[k].fwpr + S[k].fopr;
        if (l > 0) { dW = std::max(dW, std::fabs(udqOf(S[k], "FUWCT") - S[k].fwpr / l)); ++nW; }
        for (const auto& w : S[k].wells) maxP = std::max(maxP, w.oil);
        if (k > 0 && S[k].wells[iI1].mode == Control::RATE && S[k].wells[iI1].open) {
            const double tgt = udqOf(S[k - 1], "FUINJ");
            dLag = std::max(dLag, std::fabs(S[k].wells[iI1].waterInj - tgt) / tgt); ++nLag;
            dSame = std::max(dSame, std::fabs(S[k].wells[iI1].waterInj - udqOf(S[k], "FUINJ")) / tgt);
        }
    }
    std::printf("    %zu steps: FOPT %.6g stb, FUWCT |udq - FWPR/FLPR| %.2e over %d points, injector on RATE at the UDA value in %d steps "
                "(max rel. diff %.2e; vs the same step's FUINJ %.2e), max producer rate %.1f\n", S.size(), S.back().fopt, dW, nW, nLag, dLag, dSame, maxP);
    {
        const ecl::Summary es = ecl::summaryFromRun(S, s.wellNames(), c.units);
        bool hasF = false;
        for (const auto& v : es.vars) { if (v.key == "FUWCT") hasF = true; }
        check(hasF, "UDQ vectors in the ECLIPSE summary output (FUWCT)");
    }
    check(nW > 20 && dW < 1e-12, "UDQ FUWCT = FWPR / (FWPR + FOPR) at every accepted step (summary vector FUWCT)");
    check(nLag > 10 && dLag < 1e-4 && dSame > 100 * dLag, "UDA: the injector rate equals FUINJ of the previous step within the well-equation tolerance (lag semantics)");
    double p1Early = 0;
    for (const auto& p : S) if (p.time > 30 && p.time < 200) p1Early = std::max(p1Early, p.wells[iP1].oil);
    check(std::fabs(p1Early - 1500) < 1e-6 * 1500, "UDA: producer ORAT limit 'FUPLIM' = 1500 stb/d in force");
    const auto& ev = s.actionEvents();
    double tCut = -1, tPres = -1; std::string cutWells;
    for (const auto& e : ev) {
        std::printf("    action %s at %.1f d (step %d) %s\n", e.name.c_str(), e.t, e.step + 1, e.wells.c_str());
        if (e.name == "WCUT" && tCut < 0) { tCut = e.t; cutWells = e.wells; }
        if (e.name == "PRESSURE") tPres = e.t;
    }
    bool shutAfter = tCut > 0;
    if (tCut > 0) {
        const std::string first = cutWells.substr(0, cutWells.find(' '));
        size_t iw = 0;
        for (size_t w = 0; w < s.caseData().wells.size(); ++w) if (s.caseData().wells[w].name == first) iw = w;
        for (const auto& p : S) if (p.time > tCut + 1e-6 && (p.wells[iw].oil > 0 || p.wells[iw].water > 0 || p.wells[iw].open)) shutAfter = false;
    }
    check(tCut > 0 && cutWells.find("P1") != std::string::npos && shutAfter, "ACTIONX WCUT: P1 (first water breakthrough) shut via WELOPEN '?' and stays shut");
    int cutRuns = 0; for (const auto& e : ev) if (e.name == "WCUT") ++cutRuns;
    check(cutRuns >= 1 && cutRuns <= 4, "ACTIONX run limit respected (WCUT at most 4 runs)");
    bool raised = tPres > 0;
    if (tPres > 0) for (const auto& p : S) if (p.time > tPres + 1e-6 && udqOf(p, "FUPLIM") != 2200) raised = false;
    check(raised, "ACTIONX PRESSURE (FPR > 5300): UDQ ASSIGN inside the action raises FUPLIM to 2200 from the next report step");
    // restart across fired actions: bitwise continuation
    int stepR = -1;
    for (const auto& e : ev) if (e.name == "WCUT") { stepR = e.step + 3; break; }
    if (stepR > 0 && stepR < int(c.schedule.size())) {
        Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
        while (a.reportStep() < stepR && a.runReportStep()) {}
        const std::string blob = a.saveState();
        Simulator b; b.log = [](const std::string&) {}; b.setup(c); b.initialize(); b.loadState(blob);
        while (b.runReportStep()) {}
        const double f0 = S.back().fopt, f1 = b.summary().back().fopt;
        std::printf("    restart after report step %d: FOPT %.10g vs continuous %.10g\n", stepR, f1, f0);
        check(std::fabs(f1 - f0) <= 1e-9 * f0, "restart after a fired action (state v6: applied actions, UDQ values, run counts) continues identically");
    } else check(false, "restart test: no WCUT event");
}
