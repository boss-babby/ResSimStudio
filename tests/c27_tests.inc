// =====================================================================================================================
// C27 tests (industry comparison with CMG / TOUGH3): target-change time-step control (CMG NORM), economic limits +
// well tests (WECON / WTEST), van Genuchten - Mualem tables (TOUGH), iTOUGH2-style inversion statistics, studies with
// the coupled mechanics (seismicity responses), example library listing
// =====================================================================================================================
#include "ressim/coupled_bo.h"
#include "ressim/tools.h"
#include "ressim/examples.h"
#include <set>

static void testC27TimestepTargets() {
    std::printf("== C27 target-change time-step control (CMG NORM / Aziz-Settari)\n");
    Case c = smallCase();
    c.tuning.maxDt = 60;
    const RunResult r0 = runSingle(c);
    Case t = c; t.tuning.dsTarget = 0.02; t.tuning.dpTarget = 0;
    const RunResult r1 = runSingle(t);
    Case u = c; u.tuning.dpTarget = 20; // case pressure units
    const RunResult r2 = runSingle(u);
    const double f0 = r0.summary.back().fopt, f1 = r1.summary.back().fopt, f2 = r2.summary.back().fopt;
    std::printf("    off: %d steps FOPT %.6g | dS 0.02: %d steps FOPT %.6g | dP 20: %d steps FOPT %.6g\n", r0.stats.timesteps, f0,
                r1.stats.timesteps, f1, r2.stats.timesteps, f2);
    check(r0.ok && r1.ok && r2.ok, "runs with target-change control complete");
    check(r1.stats.timesteps > r0.stats.timesteps, "a tight saturation target takes more (smaller) time steps");
    check(std::fabs(f1 / f0 - 1) < 0.02 && std::fabs(f2 / f0 - 1) < 0.02, "same production within 2 % (time-step control only)");
    check(r1.stats.lastDS > 0 && r1.stats.lastDS < 0.2, "saturation change of the last step is bounded");
    // deck round trip of the two new TUNING items
    const std::string txt = projectText(t);
    const Case back = projectFromText(txt, "");
    check(std::fabs(back.tuning.dsTarget - 0.02) < 1e-12 && back.tuning.dpTarget == 0, "RESSIM TUNING dsTarget / dpTarget round trip");
}

static Case c27EconCase(double maxWct, const std::string& action, double testInterval) {
    Case c = smallCase();
    ExtraKeyword k; k.name = "WECON"; k.section = DeckSection::Schedule; k.step = 0;
    char b[200]; std::snprintf(b, sizeof b, " 'PROD' 1* 1* %g 1* 1* '%s' 'NO' /\n", maxWct, action.c_str());
    k.body = b;
    c.extraKeywords.push_back(k);
    if (testInterval > 0) {
        ExtraKeyword w; w.name = "WTEST"; w.section = DeckSection::Schedule; w.step = 0;
        std::snprintf(b, sizeof b, " 'PROD' %g 'E' 3 /\n", testInterval);
        w.body = b;
        c.extraKeywords.push_back(w);
    }
    return c;
}

static void testC27Wecon() {
    std::printf("== C27 economic limits (WECON) and well tests (WTEST) in the ResSim engine\n");
    Case c0 = smallCase();
    const RunResult r0 = runSingle(c0);
    double wctMax = 0; for (const auto& s : r0.summary) wctMax = std::max(wctMax, s.fwct);
    const double lim = 0.5 * wctMax;
    std::printf("    reference: final FOPT %.6g, max water cut %.3f -> limit %.3f\n", r0.summary.back().fopt, wctMax, lim);
    auto runSim = [&](const Case& c, RunStats& st, double& fopt, double& lastQo) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        st = s.stats(); fopt = s.summary().back().fopt; lastQo = s.summary().back().fopr;
        return s.finished();
    };
    RunStats s1, s2, s3; double f1 = 0, f2 = 0, f3 = 0, q1 = 0, q2 = 0, q3 = 0;
    runSim(c27EconCase(lim, "WELL", 0), s1, f1, q1);
    std::printf("    WECON WCT %.3f WELL: shut-ins %d, FOPT %.6g, final oil rate %.4g\n", lim, s1.econShut, f1, q1);
    check(s1.econShut == 1 && q1 == 0 && f1 < r0.summary.back().fopt, "the producer is shut at the water-cut limit and stays shut");
    runSim(c27EconCase(lim, "NONE", 0), s2, f2, q2);
    check(s2.econShut == 0 && std::fabs(f2 / r0.summary.back().fopt - 1) < 1e-9, "workover NONE: the limit is reported only (identical run)");
    runSim(c27EconCase(lim, "WELL", 200), s3, f3, q3);
    std::printf("    with WTEST 200 d 'E' 3: shut-ins %d, re-openings %d, FOPT %.6g\n", s3.econShut, s3.econTests, f3);
    check(s3.econTests >= 1 && s3.econTests <= 3 && s3.econShut == s3.econTests + 1 && f3 > f1, "WTEST re-opens the well (at most 3 tests), it is shut again while above the limit");
    // restart from the middle of the run reproduces the continuous run bitwise (economic-limit state in the restart)
    {
        const Case c = c27EconCase(lim, "WELL", 200);
        Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
        std::string blob;
        for (int k = 0; a.runReportStep(); ++k) if (k == 11) blob = a.saveState();
        Simulator b; b.log = [](const std::string&) {}; b.setup(c); b.initialize(); b.loadState(blob);
        while (b.runReportStep()) {}
        const double da = a.summary().back().fopt, db = b.summary().back().fopt;
        std::printf("    restart at step 12: FOPT %.10g vs %.10g\n", da, db);
        check(da == db, "restart with an economically shut well is bitwise identical");
    }
    // deck round trip: WECON / WTEST stay verbatim keywords of the exported deck and drive the engine again
    {
        const Case c = c27EconCase(lim, "WELL", 200);
        const std::string path = "/tmp/c27_wecon.DATA";
        exportEclipseDeck(c, path);
        std::vector<std::string> warn;
        const Case d = importEclipseDeck(path, &warn);
        int n = 0; for (const auto& k : d.extraKeywords) if (k.name == "WECON" || k.name == "WTEST") ++n;
        RunStats s4; double f4 = 0, q4 = 0; runSim(d, s4, f4, q4);
        std::printf("    exported + re-imported deck: %d WECON/WTEST keywords, shut-ins %d, FOPT %.6g\n", n, s4.econShut, f4);
        check(n == 2 && s4.econShut == s3.econShut && std::fabs(f4 / f3 - 1) < 1e-6, "WECON / WTEST deck round trip");
    }
}

static void testC27VanGenuchten() {
    std::printf("== C27 van Genuchten - Mualem tables (TOUGH RPCAP / ICP 7)\n");
    VanGenuchtenInput v; v.swr = 0.2; v.snr = 0.05; v.m = 0.5; v.p0 = 0.2; v.pcMax = 5; v.rows = 24;
    std::vector<SatRow> t; makeVanGenuchtenSWOF(v, t);
    double worst = 0; bool mono = true;
    for (size_t i = 0; i < t.size(); ++i) {
        const double se = std::clamp((t[i].s - v.swr) / (1 - v.swr - v.snr), 0.0, 1.0);
        const double krw = se <= 0 ? 0 : std::sqrt(se) * std::pow(1 - std::pow(1 - std::pow(se, 1 / v.m), v.m), 2);
        worst = std::max(worst, std::fabs(t[i].krw - krw));
        if (i) mono = mono && t[i].krw >= t[i - 1].krw - 1e-15 && t[i].kro <= t[i - 1].kro + 1e-15 && t[i].pc <= t[i - 1].pc + 1e-15;
    }
    // hand value: Se = 0.5, m = 0.5 -> krw = sqrt(0.5) (1 - sqrt(0.75))^2 = 0.012692
    VanGenuchtenInput h = v; h.swr = 0; h.snr = 0; h.p0 = 0;
    std::vector<SatRow> th; makeVanGenuchtenSWOF(h, th);
    double kHalf = -1; for (size_t i = 1; i < th.size(); ++i) if (th[i - 1].s <= 0.5 && th[i].s >= 0.5) kHalf = th[i - 1].krw + (th[i].krw - th[i - 1].krw) * (0.5 - th[i - 1].s) / std::max(1e-12, th[i].s - th[i - 1].s);
    std::printf("    SWOF %zu rows: max |krw - formula| %.2g, Pc(Swr) %.3g (cap 5), krw(0.5) ~ %.4f (exact 0.01269, table interpolation)\n", t.size(), worst, t.front().pc, kHalf);
    check(worst < 1e-12 && mono && std::fabs(t.front().pc - 5) < 1e-12 && t.back().kro == 0 && t.back().pc == 0, "SWOF: Mualem krw, monotone curves, capped Pc, kro(1) = 0");
    check(std::fabs(kHalf - 0.012692) < 0.004, "krw(Se = 0.5, m = 0.5) matches the closed form");
    VanGenuchtenInput g = v; g.swr = 0.25; g.nonWetting = 1;
    std::vector<SatRow> sg; makeVanGenuchtenSGOF(g, 0.2, sg);
    check(sg.front().s == 0 && sg.front().krw == 0 && sg.front().pc == 0 && std::fabs(sg.back().s - 0.8) < 1e-12 && sg.back().kro == 0,
          "SGOF: krg = Pcgo = 0 at Sg = 0, krog = 0 at the liquid residual");
    Case c = smallCase(); makeVanGenuchtenSWOF(v, c.props.swof);
    const RunResult r = runSingle(c);
    std::printf("    small case with the VG SWOF: ok %d, FOPT %.6g\n", r.ok, r.ok ? r.summary.back().fopt : 0.0);
    check(r.ok, "a case runs with van Genuchten - Mualem tables");
}

static void testC27HmStatistics() {
    std::printf("== C27 iTOUGH2-style inversion statistics (covariance, correlation, CSS, residual analysis)\n");
    Case base = smallCase();
    std::vector<StudyParam> P(3);
    P[0].kind = ParamKind::PermXY; P[0].base = 1; P[0].low = 0.3; P[0].high = 3; P[0].log = true;
    P[1].kind = ParamKind::KrwMax; P[1].base = 1; P[1].low = 0.4; P[1].high = 1.6; P[1].log = true;
    P[2].kind = ParamKind::WaterVisc; P[2].base = 1; P[2].low = 0.5; P[2].high = 2; P[2].log = false;
    const std::vector<double> truthV = {1.8, 0.6, 1.0};
    const Case truthCase = applyParams(base, P, truthV);
    const RunResult truth = runSingle(truthCase);
    base.observations = syntheticObservations(truth, truthCase, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::FPR}, 0.01, 3);
    StudyCallbacks cb;
    HMStatsOptions o;
    const auto t0 = std::chrono::steady_clock::now();
    const HMStatistics S = hmStatistics(base, P, truthV, o, cb);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!S.ok) { std::printf("    error: %s\n", S.error.c_str()); check(false, "inversion statistics"); return; }
    std::printf("    m %d n %d runs %d (%.1fs): objective %.4g, s0^2 %.4g, cond %.3g\n", S.nObs, S.nPar, S.runs, sec, S.objective, S.s02, S.condition);
    for (int k = 0; k < S.nPar; ++k)
        std::printf("    %-28s value %.4g  std %.3g%s  95%% [%.4g, %.4g]  CSS %.3g\n", S.name[size_t(k)].c_str(), S.value[size_t(k)], S.stdDev[size_t(k)],
                    S.logSpace[size_t(k)] ? " (log10)" : "", S.ci95lo[size_t(k)], S.ci95hi[size_t(k)], S.css[size_t(k)]);
    std::printf("    corr(k, krw) %+.3f  corr(k, muw) %+.3f  corr(krw, muw) %+.3f;  outliers %zu;  residual mean %.3g std %.3g\n",
                S.corr[0][1], S.corr[0][2], S.corr[1][2], S.outliers.size(), S.resMean, S.resStd);
    bool sym = true, inCI = true; double share = 0;
    for (int a = 0; a < S.nPar; ++a) {
        sym = sym && std::fabs(S.corr[size_t(a)][size_t(a)] - 1) < 1e-9;
        for (int b2 = 0; b2 < S.nPar; ++b2) sym = sym && std::fabs(S.corr[size_t(a)][size_t(b2)] - S.corr[size_t(b2)][size_t(a)]) < 1e-9 && std::fabs(S.corr[size_t(a)][size_t(b2)]) <= 1 + 1e-9;
        inCI = inCI && S.ci95lo[size_t(a)] <= truthV[size_t(a)] && truthV[size_t(a)] <= S.ci95hi[size_t(a)] && std::isfinite(S.stdDev[size_t(a)]) && S.stdDev[size_t(a)] > 0;
    }
    for (const auto& c : S.contribution) share += c.second;
    check(S.nPar == 3 && S.runs == 4 && sym, "covariance / correlation matrix: symmetric, unit diagonal, |r| <= 1");
    check(inCI, "the 95 % intervals at the truth contain the true values (data with 1 % noise)");
    check(std::fabs(share - 100) < 1e-6 && S.contribution.size() >= 3, "objective shares per data set sum to 100 %");
    // s0^2 ~ 1 when the weights match the noise: the synthetic noise is 1 %, autoSigma is larger -> s0^2 < 1, finite
    check(S.s02 > 0 && std::isfinite(S.s02) && S.leverage.size() == size_t(S.nObs), "error variance and leverages");
    double hsum = 0; for (double v : S.leverage) hsum += v;
    std::printf("    sum of leverages %.4f (= n = 3 for a full-rank J)\n", hsum);
    check(std::fabs(hsum - 3.0) < 1e-6, "trace of the hat matrix equals the number of parameters");
}

static void testC27CoupledStudy() {
    std::printf("== C27 studies evaluate the coupled mechanics (FAULT_SLIP.DATA): seismicity responses\n");
    const std::string path = findExample("FAULT_SLIP.DATA");
    if (path.empty()) { std::printf("    (examples/FAULT_SLIP.DATA not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    check(coupledbo::needed(c), "FAULT_SLIP needs the coupled runner");
    const auto t0 = std::chrono::steady_clock::now();
    const RunResult r = runSingle(c);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Economics eco;
    const double ev = responseValue(r, c, Response::SeisEvents, eco), mw = responseValue(r, c, Response::MaxMw, eco),
                 dc = responseValue(r, c, Response::MaxDCFS, eco), up = responseValue(r, c, Response::Uplift, eco);
    std::printf("    runSingle: ok %d, %zu KPIs, expected events %.4g, Mw max %.3g, max dCFS %.4g, uplift %s (%.1fs)\n", r.ok, r.kpi.size(), ev, mw, dc,
                std::isfinite(up) ? "finite" : "n/a (no GEOMECH)", sec);
    check(r.ok && !r.kpi.empty() && std::isfinite(ev) && ev > 0 && std::isfinite(mw) && std::isfinite(dc) && dc > 0 && !std::isfinite(up),
          "runSingle returns the fault KPIs; responses without the physics are NaN");
    // same numbers as the CLI / C26 headless runner: max dCFS 5.296 bar
    check(std::fabs(dc - 5.296) < 0.01, "max dCFS identical to the C26 coupled CLI run (5.296)");
    // tornado on the injection rate: more injection -> more expected events
    std::vector<StudyParam> P(1);
    P[0].kind = ParamKind::InjRate; P[0].target = -1; P[0].base = 1; P[0].low = 0.5; P[0].high = 1.5; P[0].log = false;
    StudyCallbacks cb;
    const TornadoResult T = runTornado(c, P, Response::SeisEvents, eco, cb);
    std::printf("    tornado (injection rate x0.5 / x1.5): expected events %.4g / %.4g (base %.4g)\n", T.ok && !T.bars.empty() ? T.bars[0].low : NAN,
                T.ok && !T.bars.empty() ? T.bars[0].high : NAN, T.base);
    check(T.ok && T.bars.size() == 1 && T.bars[0].okLow && T.bars[0].okHigh && T.bars[0].high > T.bars[0].low, "seismicity tornado: more injection, more expected events");
}

static void testC27Examples() {
    std::printf("== C27 example library listing (CLI --list-examples / --example)\n");
    const auto L = listExamples();
    std::set<std::string> ids; for (const auto& e : L) ids.insert(e.id);
    std::printf("    %zu built-in examples\n", L.size());
    check(L.size() >= 10 && ids.size() == L.size(), "example ids are unique");
    const Case c = makeExample(L.front().id);
    exportEclipseDeck(c, "/tmp/c27_example.DATA");
    std::vector<std::string> warn;
    const Case d = importEclipseDeck("/tmp/c27_example.DATA", &warn);
    check(validateCase(d).empty(), "an exported example deck validates");
}

// ---------------------------------------------------------------- C27 passive tracers (CMG IMEX tracer module / TOUGH3 EOS1 + EOS7R decay)
static Case c27TracerCase(double halfLife, double kd, double init, bool inject) {
    Case c = smallCase();
    ExtraKeyword t; t.name = "TRACER"; t.section = DeckSection::Props; t.body = " 'TW' 'WAT' /\n/\n";
    c.extraKeywords.push_back(t);
    if (inject) {
        ExtraKeyword w; w.name = "WTRACER"; w.section = DeckSection::Schedule; w.step = 0; w.body = " 'INJ' 'TW' 1.0 /\n/\n";
        c.extraKeywords.push_back(w);
    }
    if (halfLife > 0 || kd > 0 || init > 0) { Case::TracerProp P; P.name = "TW"; P.halfLife = halfLife; P.kd = kd; P.init = init; c.tracerProps.push_back(P); }
    return c;
}

static void testC27Tracer() {
    std::printf("== C27 passive tracers in the black-oil engine (TRACER / WTRACER / RESSIM TRACERPROP)\n");
    struct Out { double inj = 0, prod = 0, dec = 0, inPlace0 = 0, inPlace = 0, cmax = 0, tBreak = -1, fopt = 0; std::vector<double> dts; };
    auto run = [](const Case& c, Out& o, std::string* blobAt11 = nullptr, const std::string* restartFrom = nullptr) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        if (restartFrom) s.loadState(*restartFrom);
        o.inPlace0 = s.tracerInPlace(0);
        for (int k = 0; s.runReportStep(); ++k) if (blobAt11 && k == 11) *blobAt11 = s.saveState();
        const auto& S = s.summary();
        o.inj = S.back().fcompIT.empty() ? 0 : S.back().fcompIT[0];
        o.prod = S.back().fcompPT.empty() ? 0 : S.back().fcompPT[0];
        o.dec = s.tracerDecayed(0); o.inPlace = s.tracerInPlace(0); o.fopt = S.back().fopt;
        for (double v : s.tracerConc()) o.cmax = std::max(o.cmax, v);
        for (const auto& p : S) { if (o.tBreak < 0 && !p.fcompP.empty() && p.fcompP[0] > 1e-9) o.tBreak = p.time; o.dts.push_back(p.dt); }
        return s.tracers().size();
    };
    // 1) conservative water tracer: mass balance, maximum principle, flow unchanged
    Out a; const size_t nt = run(c27TracerCase(0, 0, 0, true), a);
    const RunResult r0 = runSingle(smallCase());
    const double bal = a.inPlace + a.prod + a.dec - a.inj - a.inPlace0;
    std::printf("    TW: injected %.6g, produced %.6g, in place %.6g, balance error %.2e (rel), max c - 1 = %.2e, breakthrough %.4g d\n",
                a.inj, a.prod, a.inPlace, std::fabs(bal) / a.inj, a.cmax - 1.0, a.tBreak);
    check(nt == 1 && a.inj > 0 && a.prod > 0 && a.tBreak > 0, "a water tracer injected with the water reaches the producer");
    check(std::fabs(bal) <= 1e-8 * a.inj, "tracer mass balance: in place + produced + decayed = injected (1e-8)");
    check(a.cmax <= 1.0 + 1e-7, "maximum principle: no concentration above the injected one (upwind implicit, M-matrix; 1e-7 = solver tolerance)");
    check(std::fabs(a.fopt / r0.summary.back().fopt - 1) < 1e-12, "the passive tracer does not change the flow solution");
    // 2) adsorption retards the breakthrough, decay removes tracer - and both keep the balance
    Out b; run(c27TracerCase(0, 0.5, 0, true), b);
    Out d; run(c27TracerCase(200, 0, 0, true), d);
    const double balB = b.inPlace + b.prod + b.dec - b.inj - b.inPlace0, balD = d.inPlace + d.prod + d.dec - d.inj - d.inPlace0;
    std::printf("    kd 0.5: breakthrough %.4g d (vs %.4g), produced %.6g | half-life 200 d: decayed %.6g, produced %.6g, balance %.1e / %.1e\n",
                b.tBreak, a.tBreak, b.prod, d.dec, d.prod, std::fabs(balB) / b.inj, std::fabs(balD) / d.inj);
    check(b.tBreak > a.tBreak && b.prod < a.prod, "linear adsorption retards the breakthrough and lowers the production");
    check(d.dec > 0 && d.prod < a.prod && std::fabs(balB) <= 1e-8 * b.inj && std::fabs(balD) <= 1e-8 * d.inj, "decay: tracer is removed, balances hold");
    // 3) pure decay with the wells stopped: the integrating-factor decay gives M0 exp(-lambda t) exactly, for any step
    {
        Case c = c27TracerCase(365.25, 0.2, 2.0, false);
        for (auto& s : c.schedule) for (auto& wc : s.controls) wc.status = WellStatus::Shut;
        Out e; run(c, e);
        const double lam = std::log(2.0) / 365.25;
        double disc = e.inPlace0, t = 0;
        for (double dt : e.dts) if (dt > 0) { disc *= std::exp(-lam * dt); t += dt; }
        const double exact = e.inPlace0 * std::exp(-lam * t);
        std::printf("    decay only (%.0f d, 91-day steps = 0.25 half-lives): M/M0 %.10f, step product %.10f, exp(-lambda t) %.10f, max c %.8f\n",
                    t, e.inPlace / e.inPlace0, disc / e.inPlace0, exact / e.inPlace0, e.cmax);
        check(e.prod == 0 && std::fabs(e.inPlace / disc - 1) < 1e-9, "decay without flow: product of the per-step decay factors");
        check(std::fabs(e.inPlace / exact - 1) < 1e-9, "decay without flow: exact exp(-lambda t) for any step length (integrating factor)");
        check(std::fabs(e.cmax / (2.0 * std::exp(-lam * t)) - 1) < 1e-8, "uniform initial concentration (RESSIM TRACERPROP init) decays cell by cell");
    }
    // 4) restart (tracer state in the v4 restart) and deck round trip (TRACER / WTRACER verbatim + RESSIM TRACERPROP)
    {
        const Case c = c27TracerCase(200, 0.5, 0, true);
        std::string blob; Out f; run(c, f, &blob);
        Out g; run(c, g, nullptr, &blob);
        std::printf("    restart at step 12: in place %.12g vs %.12g, decayed %.12g vs %.12g\n", f.inPlace, g.inPlace, f.dec, g.dec);
        // the restarted run re-accumulates the summary from the restart step: compare the final in-place amount + decay
        check(f.inPlace == g.inPlace && f.dec == g.dec, "restart from the middle of the run: tracer state bitwise identical");
        const std::string path = "/tmp/c27_tracer.DATA";
        exportEclipseDeck(c, path);
        std::vector<std::string> warn;
        const Case dk = importEclipseDeck(path, &warn);
        int n = 0; for (const auto& k : dk.extraKeywords) if (k.name == "TRACER" || k.name == "WTRACER") ++n;
        Out h; run(dk, h);
        std::printf("    exported + re-imported deck: %d TRACER/WTRACER keywords, %zu TRACERPROP, produced %.8g vs %.8g\n", n, dk.tracerProps.size(), h.prod, f.prod);
        check(n == 2 && dk.tracerProps.size() == 1 && std::fabs(dk.tracerProps[0].kd - 0.5) < 1e-12 && std::fabs(h.prod / f.prod - 1) < 1e-6,
              "TRACER / WTRACER / RESSIM TRACERPROP deck round trip");
    }
}
