// Engine regression tests: material balance (3-phase demo) and Buckley-Leverett (1D waterflood).
#include "ressim/simulator.h"
#include "ressim/io.h"
#include "ressim/tools.h"
#include "ressim/study.h"
#include "ressim/streamline.h"
#include "ressim/examples.h"
#include "ressim/grid.h"
#include "ressim/ecl.h"
#include <cstdio>
#include <cmath>
#include <chrono>
#include <string>

using namespace ressim;

static int failures = 0;
static SummaryPoint runSteps(const Case& c, int steps, Simulator& sim);
static void check(bool ok, const char* what) { std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what); if (!ok) ++failures; }

static void testDemo(bool verbose) {
    std::printf("== Demo case: 3-phase live oil waterflood\n");
    Simulator sim;
    if (verbose) sim.log = [](const std::string& s) { std::printf("%s\n", s.c_str()); };
    sim.setup(makeDemoCase());
    sim.initialize();
    const SummaryPoint s0 = sim.summary().front();
    auto t0 = std::chrono::steady_clock::now();
    while (sim.runReportStep()) {
        const SummaryPoint& s = sim.summary().back();
        std::printf("  t=%7.1f d  FPR=%7.1f  FOPR=%8.1f  FWPR=%8.1f  FGOR=%6.3f  FWIR=%8.1f  BHP(P)=%7.1f %s  steps=%d cuts=%d\n",
                    s.time, s.fpr, s.fopr, s.fwpr, s.fgor, s.fwir, s.wells[0].bhp, toString(s.wells[0].mode),
                    sim.stats().timesteps, sim.stats().cuts);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const SummaryPoint& e = sim.summary().back();
    const double mbO = (s0.oip - e.oip - e.fopt) / s0.oip;
    const double mbW = (s0.wip - e.wip - e.fwpt + e.fwit) / s0.wip;
    const double mbG = (s0.gip - e.gip - e.fgpt + e.fgit) / s0.gip;
    std::printf("  material balance error: oil %.2e  water %.2e  gas %.2e\n", mbO, mbW, mbG);
    std::printf("  %d timesteps, %d cuts, %d Newton, %d linear its, %.2f s\n", sim.stats().timesteps, sim.stats().cuts,
                sim.stats().newton, sim.stats().linear, secs);
    check(sim.finished(), "run reached end of schedule");
    check(std::fabs(mbO) < 1e-4 && std::fabs(mbW) < 1e-4 && std::fabs(mbG) < 1e-4, "material balance < 1e-4");
}

static void testBL() {
    std::printf("== Buckley-Leverett 1D (200 cells)\n");
    const int n = 200;
    Simulator sim;
    sim.setup(makeBuckleyLeverettCase(n));
    sim.initialize();
    for (int s = 0; s < 5; ++s) sim.runReportStep();
    const Snapshot& sn = sim.snapshots().back();
    const double pv = 1000.0 * 100 * 10 * 0.2 / 5.614583;
    const double injected = sim.summary().back().fwit * 1.0;    // Bw ~ 1
    const double pvi = injected / pv;
    const double sf = 1.0 / std::sqrt(2.0);
    const double fwf = sf * sf / (sf * sf + (1 - sf) * (1 - sf));
    const double xf = pvi * (fwf / sf);                         // dimensionless front position
    int idx = n - 1;
    for (int i = 0; i < n; ++i) if (sn.sw[i] < 0.5 * sf) { idx = i; break; }
    const double xnum = (idx + 0.5) / n;
    double swMax = 0; for (int i = 0; i < n; ++i) swMax = std::max(swMax, double(sn.sw[i]));
    std::printf("  PVI=%.3f  analytic front x/L=%.3f  numeric=%.3f  Sw(inj cell)=%.3f\n", pvi, xf, xnum, swMax);
    check(std::fabs(xnum - xf) < 0.04, "front position within 4% of analytic");
}

static void testIO() {
    std::printf("== Project I/O round trip + ECLIPSE export\n");
    Case c = makeDemoCase();
    const std::string a = projectToString(c);
    const std::string b = projectToString(projectFromString(a));
    check(a == b, "save -> load -> save is identical");
    extern Case featureCase();
    const std::string fa = projectToString(featureCase());
    check(fa == projectToString(projectFromString(fa)) && fa.find("ZCORN") != std::string::npos && fa.find("AQUIFER") != std::string::npos,
          "corner-point/fault/aquifer project round trip is identical");
    const std::string deck = eclipseDeckToString(c);
    check(deck.find("PVTO") != std::string::npos && deck.find("WCONPROD") != std::string::npos, "ECLIPSE deck has PVTO/WCONPROD");
}

Case featureCase() {
    Case c = makeDemoCase();
    c.title = "Feature case: corner-point + fault + aquifers";
    c.startDay = 15; c.startMonth = 3; c.startYear = 2026;
    addFaultThrow(c.grid, 10, 12.5, 0.5, "F1");
    c.grid.nnc.push_back({0, 0, 0, 19, 19, 3, 0.25});
    AquiferInput a; a.id = 1; a.type = AquiferType::Fetkovich; a.volume = 2e9; a.ct = 1e-5; a.pi = 50;
    AquiferConnection ac; ac.i1 = 0; ac.i2 = 19; ac.j1 = 0; ac.j2 = 0; ac.k1 = 0; ac.k2 = 3; ac.face = 'Y'; ac.minus = true; a.conns = {ac};
    AquiferInput b = a; b.id = 2; b.type = AquiferType::CarterTracy; b.perm = 150; b.poro = 0.22; b.rInner = 3000; b.thickness = 100;
    b.angle = 90; b.visc = 0.4; b.p0 = 4900; ac.j1 = ac.j2 = 19; ac.minus = false; ac.mult = 0.5; b.conns = {ac};
    c.aquifers = {a, b};
    c.schedule[5].controls.push_back({"PROD", WellStatus::Open, Control::LRAT, 8000, 500, -1, 8500, -1, -1, 900, true});
    c.tuning.dpMaxRel = 0.15; c.tuning.threads = 2; c.tuning.linearSolver = 1;
    c.wells[0].piMult = 1.25;
    StudyParam sp; sp.kind = ParamKind::PermXY; sp.k1 = sp.k2 = 2; sp.low = 0.25; sp.high = 4; sp.label = "Layer 3 perm";
    StudyParam sq; sq.kind = ParamKind::OWC; sq.base = 9000; sq.low = 8950; sq.high = 9100; sq.log = false; sq.active = false;
    c.studyParams = {sp, sq};
    c.observations.push_back({ObsQty::WBHP, "PROD", 182.625, 3250.5, 25});
    c.observations.push_back({ObsQty::FPR, "", 365.25, 3900, -1});
    c.observations.push_back({ObsQty::WWCT, "PROD", 730.5, 0.125, 0.02});
    return c;
}

static void testCPR() {
    std::printf("== Linear solvers: ILU0-BiCGStab vs CPR-AMG\n");
    Case c = makeDemoCase();
    c.schedule.resize(6);
    double fpr[2]; int lin[2], nwt[2];
    for (int s = 0; s < 2; ++s) {
        c.tuning.linearSolver = s;
        Simulator sim; sim.setup(c); sim.initialize();
        while (sim.runReportStep()) {}
        fpr[s] = sim.summary().back().fpr; lin[s] = sim.stats().linear; nwt[s] = sim.stats().newton;
        std::printf("  %s: Newton %d, linear %d, FPR %.4f%s\n", s ? "CPR" : "ILU", nwt[s], lin[s], fpr[s], s && sim.stats().cpr ? " (CPR active)" : "");
    }
    check(std::fabs(fpr[0] - fpr[1]) < 0.01 && lin[1] < lin[0], "CPR gives the same solution with fewer linear iterations");
    const Case f = featureCase();
    const Case g = eclipseDeckFromString(eclipseDeckToString(f));
    check(g.tuning.linearSolver == 1 && std::fabs(g.wells[0].piMult - 1.25) < 1e-12 && g.observations.size() == 3 &&
              g.observations[0].qty == ObsQty::WBHP && g.observations[1].well.empty() && std::fabs(g.observations[2].sigma - 0.02) < 1e-12 &&
              g.studyParams.size() == 2 && g.studyParams[0].label == "Layer 3 perm" && g.studyParams[0].k1 == 2 && !g.studyParams[1].active &&
              g.studyParams[1].base == 9000,
          "deck round trip keeps SOLVER, WPIMULT, observations and study parameters");
}

static void testDeck() {
    std::printf("== ECLIPSE deck two-way (export -> import -> export)\n");
    for (int which = 0; which < 2; ++which) {
        const Case c = which == 0 ? makeDemoCase() : featureCase();
        const std::string a = eclipseDeckToString(c);
        std::vector<std::string> warn;
        const Case c2 = eclipseDeckFromString(a, "", &warn);
        const std::string b = eclipseDeckToString(c2);
        for (auto& w : warn) std::printf("    warning: %s\n", w.c_str());
        if (a != b) {   // show first difference
            size_t i = 0; while (i < a.size() && i < b.size() && a[i] == b[i]) ++i;
            std::printf("    first difference at %zu:\n----A\n%s\n----B\n%s\n", i, a.substr(i > 80 ? i - 80 : 0, 200).c_str(), b.substr(i > 80 ? i - 80 : 0, 200).c_str());
        }
        check(a == b && warn.empty(), which == 0 ? "demo deck round trip is identical" : "corner-point/fault/NNC/aquifer/history deck round trip is identical");
        Simulator s1, s2;
        const SummaryPoint r1 = runSteps(c, 6, s1), r2 = runSteps(c2, 6, s2);
        std::printf("    6 steps: FPR %.4f vs %.4f, FOPR %.4f vs %.4f, aquifer %.5g vs %.5g\n", r1.fpr, r2.fpr, r1.fopr, r2.fopr, r1.faqt, r2.faqt);
        check(std::fabs(r1.fpr - r2.fpr) < 1e-6 * r1.fpr && std::fabs(r1.fopr - r2.fopr) < 1e-6 * std::max(1.0, r1.fopr), "imported deck runs identically");
    }
    // a hand-written deck using features the exporter never writes
    {
        std::FILE* f = std::fopen("/tmp/rb_inc.inc", "wb");
        if (!f) f = std::fopen("rb_inc.inc", "wb");
        std::fputs("-- included file\nPORO\n 300*0.2 /\n", f); std::fclose(f);
    }
    const char* deck = R"(-- foreign deck
RUNSPEC
TITLE
   Foreign deck test
DIMENS
  10 10 3 /
OIL
WATER
FIELD
TABDIMS
 1 1 20 20 /
START
  1 'JAN' 2026 /
GRID
DXV
 10*100 /
DYV
 10*100 /
DZV
 20 30 50 /
TOPS
 100*5000 /
INCLUDE
  'rb_inc.inc' /
EQUALS
 'PERMX' 100 /
 'PERMX' 500 1 10 1 10 2 2 /   layer 2 is better
/
BOX
 1 5 1 10 3 3 /
MULTIPLY
 'PERMX' 0.1 /
/
ENDBOX
COPY
 'PERMX' 'PERMY' /
 'PERMX' 'PERMZ' /
/
MULTIPLY
 'PERMZ' 0.1 /
/
SATNUM
 300*1 /
PROPS
PVTW
 3000 1.01 3e-6 0.5 0 /
PVCDO
 3000 1.2 1e-5 2.0 0 /
ROCK
 3000 4e-6 /
DENSITY
 50 64 0.06 /
SWOF
 0.2 0 1 0
 0.5 0.1 0.3 0
 0.8 0.4 0 0 /
SOLUTION
EQUIL
 5050 3000 5200 0 1* 1* /
SUMMARY
FOPR
WBHP
 'P1' /
SCHEDULE
WELSPECS
 'P1' 'G1' 10 10 1* 'OIL' /
 'I1' 'G1' 1 1 1* 'WATER' /
/
COMPDAT
 'P1' 0 0 1 3 'OPEN' 2* 0.5 /
 'I1' 2* 1 3 'OPEN' 2* 0.5 /
/
WCONPROD
 'P1' 'OPEN' 'ORAT' 800 4* 1000 /
/
WCONINJE
 'I1' 'WATER' 'OPEN' 'RATE' 1000 1* 4500 /
/
DATES
 1 'FEB' 2026 /
 1 'MAR' 2026 /
/
WELOPEN
 'I1' 'SHUT' /
/
TSTEP
 2*30 /
END
)";
    std::vector<std::string> warn;
    Case c;
    try { c = eclipseDeckFromString(deck, "/tmp/", &warn); }
    catch (const std::exception& e) { std::printf("    import error: %s\n", e.what()); check(false, "foreign deck imports"); return; }
    for (auto& w : warn) std::printf("    warning: %s\n", w.c_str());
    const int q2 = c.grid.index(3, 3, 1), q3 = c.grid.index(2, 2, 2), q3b = c.grid.index(7, 2, 2);
    std::printf("    title '%s', steps %zu (%.0f %.0f %.0f %.0f), permx L2 %.0f L3(box) %.0f L3 %.0f, permz L2 %.0f, poro %.2f, pvdo %zu rows\n",
                c.title.c_str(), c.schedule.size(), c.schedule[0].length, c.schedule[1].length, c.schedule[2].length, c.schedule[3].length,
                c.grid.permx[q2], c.grid.permx[q3], c.grid.permx[q3b], c.grid.permz[q2], c.grid.poro[0], c.props.pvdo.size());
    check(c.title == "Foreign deck test" && c.schedule.size() == 4 && c.schedule[0].length == 31 && c.schedule[1].length == 28 &&
          c.grid.permx[q2] == 500 && std::fabs(c.grid.permx[q3] - 10) < 1e-9 && c.grid.permx[q3b] == 100 && std::fabs(c.grid.permz[q2] - 50) < 1e-9 &&
          c.grid.poro[0] == 0.2 && c.wells.size() == 2 && c.wells[0].completions[0].i == 9 && c.schedule[2].controls.size() == 1 &&
          c.schedule[2].controls[0].status == WellStatus::Shut && c.startYear == 2026,
          "foreign deck: DXV/EQUALS/BOX/MULTIPLY/COPY/INCLUDE/DATES/WELOPEN/PVCDO");
    Simulator sim; sim.setup(c); sim.initialize(); while (sim.runReportStep()) {}
    std::printf("    ran %d steps, FOPR %.1f, FPR %.1f\n", sim.stats().timesteps, sim.summary().back().fopr, sim.summary().back().fpr);
    check(sim.finished(), "foreign deck runs to the end");
}

static Case smallCase() {
    Case c = makeDemoCase();
    GeomInput gi; gi.nx = gi.ny = 10; gi.nz = 2; gi.dx = gi.dy = 300; gi.dz = {50, 50};
    generateGeometry(gi, c.grid);
    const size_t n = size_t(c.grid.size());
    c.grid.poro.assign(n, 0.2); c.grid.permx.assign(n, 200); c.grid.permy.clear(); c.grid.permz.assign(n, 20); c.grid.ntg.clear();
    for (auto& w : c.wells) {
        if (w.type == WellType::Producer) { w.i = w.j = 9; }
        for (auto& cp : w.completions) { cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = 1; }
    }
    for (auto& s : c.schedule) for (auto& wc : s.controls) { if (wc.well == "PROD") wc.orat = 6000; if (wc.well == "INJ") wc.rate = 8000; }
    c.schedule.resize(24);
    return c;
}

static void testStudy() {
    std::printf("== Studies: batch runner, tornado, Monte Carlo, LM and ES-MDA history matching\n");
    Case base = smallCase();
    std::vector<StudyParam> P(3);
    P[0].kind = ParamKind::PermXY; P[0].base = 1; P[0].low = 0.3; P[0].high = 3; P[0].log = true;
    P[1].kind = ParamKind::KrwMax; P[1].base = 1; P[1].low = 0.4; P[1].high = 1.6; P[1].log = true;
    P[2].kind = ParamKind::WaterVisc; P[2].base = 1; P[2].low = 0.5; P[2].high = 2; P[2].active = true;
    // truth
    const std::vector<double> truthV = {1.8, 0.6, 1.0};
    const Case truthCase = applyParams(base, P, truthV);
    const RunResult truth = runSingle(truthCase);
    std::printf("  truth run ok=%d steps=%d  FOPT %.4g  WCT %.3f\n", truth.ok, truth.stats.timesteps, truth.summary.back().fopt, truth.summary.back().fwct);
    base.observations = syntheticObservations(truth, truthCase, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::FPR}, 0.01, 3);
    std::printf("  %zu synthetic observations\n", base.observations.size());
    check(truth.ok && base.observations.size() > 40, "synthetic observed data generated");

    StudyCallbacks cb;
    Economics eco;
    auto t0 = std::chrono::steady_clock::now();
    TornadoResult tr = runTornado(base, P, Response::FOPT, eco, cb);
    std::printf("  tornado: base FOPT %.5g;", tr.base);
    for (auto& b : tr.bars) std::printf("  %s [%.5g, %.5g]", paramLabel(P[b.param], base).c_str(), b.low, b.high);
    std::printf("  (%.1fs)\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(tr.ok && tr.bars.size() == 3 && tr.bars[0].okLow && tr.bars[0].okHigh, "tornado runs (7 parallel runs)");

    t0 = std::chrono::steady_clock::now();
    MonteCarloResult mc = runMonteCarlo(base, P, 12, 5, Response::RF, eco, cb);
    std::printf("  Monte Carlo RF: P90 %.4f P50 %.4f P10 %.4f, spearman", mc.p90, mc.p50, mc.p10);
    for (double v : mc.spearman) std::printf(" %.2f", v);
    std::printf("  (%.1fs)\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(mc.failed == 0 && mc.p90 <= mc.p50 && mc.p50 <= mc.p10 && mc.fanP50.size() == base.schedule.size() + 1, "Latin-hypercube Monte Carlo + fan chart");

    P[2].active = false;
    t0 = std::chrono::steady_clock::now();
    LMOptions lo; lo.maxIter = 10;
    HMResult lm = runLevenbergMarquardt(base, P, lo, cb);
    std::printf("  LM: misfit %.4g -> %.4g, kxy %.3f (truth 1.8), krw %.3f (truth 0.6), %d runs, %.1fs\n", lm.initialMisfit, lm.bestMisfit,
                lm.best[0], lm.best[1], lm.runs, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(lm.ok && lm.bestMisfit < 0.02 * lm.initialMisfit && std::fabs(lm.best[0] / 1.8 - 1) < 0.1 && std::fabs(lm.best[1] / 0.6 - 1) < 0.1,
          "Levenberg-Marquardt recovers the truth parameters");

    t0 = std::chrono::steady_clock::now();
    ESMDAOptions eo; eo.ensemble = 16; eo.assimilations = 4;
    HMResult es = runESMDA(base, P, eo, cb);
    double m0 = es.history.front().misfit, m1 = es.history.back().misfit;
    double meanK = 0, meanW = 0; int cnt = 0;
    for (size_t j = 0; j < es.ensemble.size(); ++j) if (std::isfinite(es.ensembleMisfit[j])) { meanK += es.ensemble[j][0]; meanW += es.ensemble[j][1]; ++cnt; }
    meanK /= std::max(1, cnt); meanW /= std::max(1, cnt);
    std::printf("  ES-MDA: ensemble mean misfit %.4g -> %.4g, posterior mean kxy %.3f krw %.3f, best %.4g, %d runs, %.1fs\n", m0, m1, meanK, meanW,
                es.bestMisfit, es.runs, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(es.ok && m1 < 0.1 * m0 && std::fabs(meanK / 1.8 - 1) < 0.2 && std::fabs(meanW / 0.6 - 1) < 0.2, "ES-MDA ensemble converges to the truth");
}

static void testStreamlines() {
    std::printf("== Streamlines (Pollock) and time-of-flight\n");
    // 1) homogeneous 1D-like line drive: TOF at producer = PV / q exactly (incompressible limit)
    Case c = smallCase();
    c.schedule.resize(4);
    Simulator sim;
    sim.setup(c);
    sim.initialize();
    while (sim.runReportStep()) {}
    auto t0 = std::chrono::steady_clock::now();
    StreamlineResult S = traceStreamlines(sim);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    int reached = 0; double qsum = 0;
    for (auto& L : S.lines) { if (L.sink >= 0) ++reached; qsum += L.flux; }
    double inj = 0; for (auto& w : sim.wellCellRates()) if (w.rate < 0) inj -= w.rate;
    double swept = 0, pvTot = 0;
    for (int a = 0; a < sim.grid().numActive(); ++a) { pvTot += sim.grid().poreVolRef[a]; if (S.cellTof[a] >= 0) swept += sim.grid().poreVolRef[a]; }
    std::printf("  %zu lines, %d reach PROD, flux %.1f of %.1f, max TOF %.0f d, swept PV %.1f%%, %.1f ms\n", S.lines.size(), reached, qsum, inj,
                S.maxTof, 100 * swept / pvTot, ms);
    check(S.lines.size() > 100 && reached >= int(0.95 * S.lines.size()) && std::fabs(qsum / inj - 1) < 1e-6, "streamlines connect INJ to PROD, flux conserved");
    check(S.pairFlux.size() == 1 && std::fabs(S.pairFlux[0][0] / inj - 1) < 0.05, "well-pair allocation INJ->PROD");
    // flux-weighted mean TOF at the producer ~ total PV / q (all PV swept, steady incompressible approx)
    double mt = 0; for (auto& L : S.lines) if (L.sink >= 0) mt += L.flux * L.tof.back(); mt /= qsum;
    std::printf("  mean breakthrough TOF %.0f d vs PV/q %.0f d\n", mt, pvTot / inj);
    check(mt > 0.5 * pvTot / inj && mt < 1.5 * pvTot / inj, "mean TOF consistent with PV / q");
    // monotone TOF along each line
    bool mono = true; for (auto& L : S.lines) for (size_t k = 1; k < L.tof.size(); ++k) if (L.tof[k] < L.tof[k - 1] - 1e-3f) mono = false;
    check(mono, "TOF increases monotonically along every streamline");
    // restoring an earlier report gives a valid state
    sim.restoreReport(1);
    StreamlineResult S1 = traceStreamlines(sim);
    check(!S1.lines.empty(), "trace at restored report step 1");

    // 2) corner-point faulted grid: lines cross the fault through NNCs
    Case f = featureCase();
    Simulator fs;
    fs.setup(f);
    fs.initialize();
    for (int k = 0; k < 2 && fs.runReportStep(); ++k) {}
    StreamlineResult SF = traceStreamlines(fs);
    int sinkReach = 0; for (auto& L : SF.lines) if (L.sink >= 0) ++sinkReach;
    std::printf("  feature case (corner-point, NNC=%d): %zu lines, %d reach producers, %d stalled, sources %zu sinks %zu\n", fs.grid().numNNC,
                SF.lines.size(), sinkReach, SF.stalled, SF.sources.size(), SF.sinks.size());
    check(!SF.lines.empty() && sinkReach > 0, "streamlines on a faulted corner-point grid");
}

static void testCorrelations() {
    std::printf("== Correlation PVT + heterogeneous rock, 3 years\n");
    Case c = makeDemoCase();
    CorrelationInput ci; ci.api = 32; ci.gasGravity = 0.8; ci.temperatureF = 210; ci.pbMax = 5000; ci.pMax = 9000;
    makeCorrelationPVT(ci, c.units, c.props);
    std::printf("  Hall-Yarborough z(Ppr=2,Tpr=1.5)=%.4f (Standing-Katz chart ~0.82)\n", hallYarboroughZ(2.0, 1.5));
    std::printf("  PVTO last: Rs=%.3f Pb=%.1f Bo=%.3f mu=%.3f ; PVDG@%.0f Bg=%.3f mu=%.4f\n", c.props.pvto.back().rs,
                c.props.pvto.back().pb, c.props.pvto.back().bo, c.props.pvto.back().mu, c.props.pvdg[5].p, c.props.pvdg[5].B, c.props.pvdg[5].mu);
    RockGenInput rg; rg.layerPerm = {200, 50, 500, 100}; rg.layerPoro = {0.22, 0.18, 0.25, 0.2}; rg.sigmaLn = 1.0;
    generateRock(rg, c.grid);
    c.equil.initialRs = -1; c.equil.goc = 8300;    // GOC above top -> saturated at GOC pressure
    c.schedule.resize(12);
    Simulator sim; sim.setup(c); sim.initialize();
    const SummaryPoint s0 = sim.summary().front();
    while (sim.runReportStep()) {}
    const SummaryPoint& e = sim.summary().back();
    const double mbO = (s0.oip - e.oip - e.fopt) / s0.oip, mbG = (s0.gip - e.gip - e.fgpt + e.fgit) / s0.gip;
    std::printf("  end FPR=%.1f FOPR=%.1f FGOR=%.3f  steps=%d cuts=%d  mb oil %.1e gas %.1e  %.2fs\n", e.fpr, e.fopr, e.fgor,
                sim.stats().timesteps, sim.stats().cuts, mbO, mbG, sim.stats().seconds);
    check(sim.finished() && std::fabs(mbO) < 1e-4 && std::fabs(mbG) < 1e-4, "correlation case runs with material balance");
}

static SummaryPoint runSteps(const Case& c, int steps, Simulator& sim) {
    Case cc = c; cc.schedule.resize(steps);
    sim.setup(cc); sim.initialize();
    while (sim.runReportStep()) {}
    return sim.summary().back();
}

static void testCornerPoint() {
    std::printf("== Corner-point geometry, fault throw, MULTFLT\n");
    Case cart = makeDemoCase();
    Case cp = cart; convertToCornerPoint(cp.grid);
    Grid g1, g2; g1.build(cart); g2.build(cp);
    double pv1 = 0, pv2 = 0, t1 = 0, t2 = 0;
    for (double v : g1.poreVolRef) pv1 += v;
    for (double v : g2.poreVolRef) pv2 += v;
    for (auto& c : g1.conns) t1 += c.trans;
    for (auto& c : g2.conns) t2 += c.trans;
    std::printf("  PV cart %.8g cp %.8g | sum T cart %.8g cp %.8g | conns %zu / %zu\n", pv1, pv2, t1, t2, g1.conns.size(), g2.conns.size());
    check(std::fabs(pv1 - pv2) < 1e-9 * pv1 && std::fabs(t1 - t2) < 1e-9 * t1 && g1.conns.size() == g2.conns.size(),
          "COORD/ZCORN of a box grid reproduces Cartesian PV and transmissibility");
    Simulator s1, s2;
    const SummaryPoint a = runSteps(cart, 4, s1), b = runSteps(cp, 4, s2);
    check(std::fabs(a.fpr - b.fpr) < 1e-6 * a.fpr && std::fabs(a.fopr - b.fopr) < 1e-6 * a.fopr, "corner-point run matches Cartesian run");

    // throw of half a layer (12.5 ft) at I = 10 -> every lateral connection across the fault overlaps 2 cells
    Case f = cart; addFaultThrow(f.grid, 10, 12.5, 1.0, "F1");
    Grid gf; gf.build(f);
    std::printf("  fault throw 12.5: connections %zu, non-neighbour %d\n", gf.conns.size(), gf.numNNC);
    check(gf.numNNC == 20 * 3, "half-layer throw creates one NNC per layer pair along the fault");
    double tAcross = 0, tRef = 0;
    for (auto& c : gf.conns) { int i1, j1, k1, i2, j2, k2; gf.ijk(gf.globalOf[c.a], i1, j1, k1); gf.ijk(gf.globalOf[c.b], i2, j2, k2);
        if (std::min(i1, i2) == 9 && std::max(i1, i2) == 10) tAcross += c.trans; }
    for (auto& c : g1.conns) { int i1, j1, k1, i2, j2, k2; g1.ijk(g1.globalOf[c.a], i1, j1, k1); g1.ijk(g1.globalOf[c.b], i2, j2, k2);
        if (std::min(i1, i2) == 9 && std::max(i1, i2) == 10) tRef += c.trans; }
    std::printf("  T across fault %.6g (unfaulted %.6g)\n", tAcross, tRef);
    check(tAcross < tRef && tAcross > 0.5 * tRef, "juxtaposition reduces total cross-fault transmissibility");
    Case sealed = f; sealed.grid.faults[0].mult = 0.0;
    Grid gs; gs.build(sealed);
    double tS = 0; for (auto& c : gs.conns) { int i1, j1, k1, i2, j2, k2; gs.ijk(gs.globalOf[c.a], i1, j1, k1); gs.ijk(gs.globalOf[c.b], i2, j2, k2);
        if (std::min(i1, i2) == 9 && std::max(i1, i2) == 10) tS += c.trans; }
    check(tS == 0.0, "MULTFLT 0 seals the fault (incl. throw NNCs)");
    Simulator s3; Case f2 = f; f2.schedule.resize(6);
    s3.setup(f2); s3.initialize(); const SummaryPoint s0 = s3.summary().front();
    while (s3.runReportStep()) {}
    const SummaryPoint& e = s3.summary().back();
    const double mbO = (s0.oip - e.oip - e.fopt) / s0.oip, mbW = (s0.wip - e.wip - e.fwpt + e.fwit) / s0.wip;
    check(s3.finished() && std::fabs(mbO) < 1e-4 && std::fabs(mbW) < 1e-4, "faulted corner-point case runs with material balance");
}

static void testAquifers() {
    std::printf("== Analytic aquifers (Fetkovich, Carter-Tracy)\n");
    Case base = makeDemoCase();
    for (auto& s : base.schedule) for (auto& wc : s.controls) if (wc.well == "INJ") wc.status = WellStatus::Shut;
    AquiferConnection ac; ac.i1 = 0; ac.i2 = 19; ac.j1 = 0; ac.j2 = 19; ac.k1 = 3; ac.k2 = 3; ac.face = 'Z';
    Case fet = base; AquiferInput a; a.id = 1; a.type = AquiferType::Fetkovich; a.volume = 5e9; a.ct = 1e-5; a.pi = 200; a.conns = {ac};
    fet.aquifers = {a};
    Case ct = base; AquiferInput b; b.id = 2; b.type = AquiferType::CarterTracy; b.perm = 200; b.poro = 0.25; b.ct = 1e-5;
    b.rInner = 3400; b.thickness = 100; b.angle = 360; b.conns = {ac}; ct.aquifers = {b};
    Simulator s0, s1, s2;
    const SummaryPoint r0 = runSteps(base, 8, s0), r1 = runSteps(fet, 8, s1), r2 = runSteps(ct, 8, s2);
    std::printf("  FPR none %.1f  Fetkovich %.1f (We %.4g)  Carter-Tracy %.1f (We %.4g)\n", r0.fpr, r1.fpr, r1.faqt, r2.fpr, r2.faqt);
    check(r1.faqt > 0 && r1.fpr > r0.fpr + 10, "Fetkovich aquifer supports pressure");
    check(r2.faqt > 0 && r2.fpr > r0.fpr + 10, "Carter-Tracy aquifer supports pressure");
    for (Simulator* sp : {&s1, &s2}) {
        const SummaryPoint& i = sp->summary().front(); const SummaryPoint& e = sp->summary().back();
        const double mbW = (i.wip + e.faqtS - e.wip - e.fwpt + e.fwit) / i.wip;
        const double mbO = (i.oip - e.oip - e.fopt) / i.oip;
        std::printf("  material balance incl. influx: water %.2e oil %.2e\n", mbW, mbO);
        check(std::fabs(mbW) < 1e-4 && std::fabs(mbO) < 1e-4, "material balance with aquifer influx < 1e-4");
    }
    const auto st = s1.aquiferStates();
    check(!st.empty() && st[0].pressure < st[0].cumulative + 1e30 && st[0].cells == 400, "aquifer state: 400 bottom-face connections");
}

static void testExamples() {
    std::printf("== Example library (build, validate, deck round trip, quick runs)\n");
    const auto list = listExamples();
    int ok = 0;
    for (const ExampleInfo& e : list) {
        if (e.id == "hm") continue;                      // runs a truth model; covered below
        const Case c = makeExample(e.id);
        const std::string err = validateCase(c);
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        const bool same = err.empty() && b.grid.nx == c.grid.nx && b.grid.ny == c.grid.ny && b.grid.nz == c.grid.nz &&
                          b.wells.size() == c.wells.size() && b.schedule.size() == c.schedule.size() &&
                          b.grid.cornerPoint() == c.grid.cornerPoint() && b.aquifers.size() == c.aquifers.size() &&
                          b.equilExtra.size() == c.equilExtra.size() && b.props.satExtra.size() == c.props.satExtra.size();
        if (!same) std::printf("  %s: %s\n", e.id.c_str(), err.empty() ? "round-trip mismatch" : err.c_str());
        ok += same;
    }
    check(ok == int(list.size()) - 1, "every example validates and survives an ECLIPSE deck round trip");
    for (const char* id : {"spe1", "qfs", "regions", "faulted"}) {
        const Case c = makeExample(id);
        Simulator s; s.setup(c); s.initialize();
        int n = 0; while (s.runReportStep()) ++n;
        const SummaryPoint& i = s.summary().front(); const SummaryPoint& e = s.summary().back();
        const double mbO = (i.oip - e.oip - e.fopt) / i.oip;
        std::printf("  %-8s steps %d/%zu  FPR %.2f  FOPT %.4g  oil balance %.1e\n", id, n, c.schedule.size(), e.fpr, e.fopt, mbO);
        check(n == int(c.schedule.size()) && std::fabs(mbO) < 1e-4, (std::string(id) + ": runs to the end, oil balance < 1e-4").c_str());
        if (std::string(id) == "spe1") check(std::fabs(e.fpr - 3831.8) < 5, "SPE1 example: FPR after 10 years = 3831.8 psia (same as the imported SPE1CASE1 deck)");
    }
    // per-period injected phase (WAG-style): water first, gas from step 4
    Case w = makeExample("qfs");
    w.phases = Phases{true, true, true, true};
    {
        const Case d = makeDemoCase();
        w.props = d.props; w.equil = d.equil; w.equil.datumDepth = 5015; w.equil.datumPressure = 3000; w.equil.owc = 6000; w.equil.goc = 4000;
        w.equil.initialRs = 0.9;
    }
    WellControl g; g.well = "INJ"; g.mode = Control::RATE; g.rate = 300; g.bhp = 5000; g.phase = 2;
    w.schedule[4].controls.push_back(g);
    const Case wb = eclipseDeckFromString(eclipseDeckToString(w));
    int ph = 0; for (auto& wc : wb.schedule[4].controls) if (wc.well == "INJ") ph = wc.phase;
    check(ph == 2, "deck round trip keeps the injected phase of a WCONINJE period");
    Simulator sw; sw.setup(w); sw.initialize();
    for (int k = 0; k < 6; ++k) sw.runReportStep();
    const SummaryPoint& last = sw.summary().back();
    std::printf("  WAG: FWIT %.4g  FGIR %.4g\n", last.fwit, last.fgir);
    check(last.fwit > 0 && last.fgir > 0, "injector switches from water to gas injection");
}

static void testNewtran() {
    std::printf("== NEWTRAN (geometric transmissibility) + project format 2\n");
    // flat box grid: NEWTRAN == OLDTRAN
    Case flat = makeDemoCase(); Case flatN = flat; flatN.grid.newtran = true;
    Grid a, b; a.build(flat); b.build(flatN);
    double ta = 0, tb = 0; for (auto& c : a.conns) ta += c.trans; for (auto& c : b.conns) tb += c.trans;
    check(a.conns.size() == b.conns.size() && std::fabs(ta - tb) < 1e-9 * ta, "NEWTRAN on a flat grid equals OLDTRAN");
    // steep dip: column offset (53 ft) > layer thickness -> overlap NNCs, OLDTRAN keeps logical neighbours only
    Case dip = makeDemoCase();
    GridInput& g = dip.grid;
    for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) g.tops[size_t(j) * g.nx + i] = 8000 + 53.0 * i;
    Case dipN = dip; dipN.grid.newtran = true;
    Grid go, gn; go.build(dip); gn.build(dipN);
    std::printf("  dipping grid: OLDTRAN %zu conns (%d NNC), NEWTRAN %zu conns (%d NNC)\n", go.conns.size(), go.numNNC, gn.conns.size(), gn.numNNC);
    check(go.numNNC == 0 && gn.numNNC > 0, "NEWTRAN connects displaced columns through their overlaps");
    double pvo = 0, pvn = 0; for (double v : go.poreVolRef) pvo += v; for (double v : gn.poreVolRef) pvn += v;
    check(std::fabs(pvo - pvn) < 1e-9 * pvo, "NEWTRAN keeps the pore volume");
    // deck keyword round trip
    const std::string deck = eclipseDeckToString(dipN);
    check(deck.find("NEWTRAN") != std::string::npos && eclipseDeckFromString(deck).grid.newtran, "NEWTRAN keyword is exported and imported");
    check(eclipseDeckToString(dip).find("OLDTRAN") == std::string::npos, "OLDTRAN (default, rejected by OPM) is not written");
    // project format 2 keeps everything the deck has (regions, RSVD, MULTPV, NEWTRAN, study data)
    Case rc = makeExample("regions"); rc.grid.newtran = true; rc.grid.minpv = 1.5;
    rc.studyParams = featureCase().studyParams; rc.observations = featureCase().observations;
    const std::string path = "ressim_test_project.rsim";
    saveProject(rc, path);
    const Case back = loadProject(path);
    std::remove(path.c_str());
    check(eclipseDeckToString(back) == eclipseDeckToString(rc) && back.grid.newtran && back.studyParams.size() == 2 && back.observations.size() == 3,
          "project save -> load keeps regions / NEWTRAN / MINPV / study data");
    const std::string v1 = projectToString(makeDemoCase());
    { FILE* f = std::fopen(path.c_str(), "wb"); std::fwrite(v1.data(), 1, v1.size(), f); std::fclose(f); }
    const Case old = loadProject(path); std::remove(path.c_str());
    check(projectToString(old) == v1, "format 1 projects still load");
}

// Phase A: late COMPDAT timing, THPRES, SWATINIT, true-IMPES CPR weights
static void testPhaseA() {
    std::printf("Phase A physics (late COMPDAT, THPRES, SWATINIT, true-IMPES)\n");
    {   // late COMPDAT: injector perforated at report step 12, producer re-perforated at 24, original shut at 60
        Case c = makeExample("spe1");
        int inj = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].type != WellType::Producer) inj = int(w);
        const int prod = 1 - inj;
        c.wells[size_t(inj)].completions[0].fromStep = 12;
        Completion extra = c.wells[size_t(prod)].completions[0]; extra.k1 = extra.k2 = 1; extra.fromStep = 24;
        c.wells[size_t(prod)].completions.push_back(extra);
        c.wells[size_t(prod)].completions[0].untilStep = 60;
        const std::string d1 = eclipseDeckToString(c);
        // the first import turns the OPEN control of the not yet perforated injector into SHUT (ECLIPSE rule) and keeps
        // the control re-issued at step 12 -> same behaviour, and the deck is stable from then on
        const std::string d2 = eclipseDeckToString(eclipseDeckFromString(d1));
        check(eclipseDeckToString(eclipseDeckFromString(d2)) == d2, "late COMPDAT: deck import -> export is stable");
        const Case b = eclipseDeckFromString(d2);
        bool okc = b.wells[size_t(inj)].completions[0].fromStep == 12 && b.wells[size_t(prod)].completions.size() == 2;
        for (auto& q : b.wells[size_t(prod)].completions) okc = okc && ((q.k1 == 2 && q.untilStep == 60) || (q.k1 == 1 && q.fromStep == 24));
        check(okc, "late COMPDAT: connection open / shut steps survive the deck round trip");
        Simulator s; s.log = [](const std::string&) {}; s.setup(b); s.initialize();
        while (s.runReportStep()) {}
        double before = 0, after = 0;
        for (auto& p : s.summary()) { if (p.time < 360) before = std::max(before, p.fgir); if (p.time > 400) after = std::max(after, p.fgir); }
        std::printf("  max FGIR before step 12: %g, after: %g\n", before, after);
        check(before == 0 && after > 0, "late COMPDAT: injection starts when the connection opens");
    }
    auto regions = [](bool wells) {
        Case c = makeExample("regions");
        for (auto& f : c.grid.faults) f.mult = 1.0;
        if (!wells) for (auto& st : c.schedule) for (auto& w : st.controls) w.status = WellStatus::Shut;
        return c;
    };
    auto maxChange = [](Case c, double& dp, double& ds) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        for (int k = 0; k < 6 && s.runReportStep(); ++k) {}
        const auto& sn = s.snapshots(); dp = ds = 0;
        for (size_t a = 0; a < sn[0].pressure.size(); ++a) {
            dp = std::max(dp, double(std::fabs(sn.back().pressure[a] - sn[0].pressure[a])));
            ds = std::max(ds, double(std::fabs(sn.back().sw[a] - sn[0].sw[a])));
        }
        return sn[0].sw;
    };
    double dp0, ds0, dp1, ds1;
    maxChange(regions(false), dp0, ds0);
    Case t = regions(false); t.thpres.push_back({1, 2, -1});
    maxChange(t, dp1, ds1);
    std::printf("  compartments with an open fault: max dp %.3f -> %.3f with THPRES (defaulted)\n", dp0, dp1);
    check(dp0 > 2 && dp1 < 0.2 * dp0, "THPRES: defaulted threshold keeps the initial equilibrium of the EQLNUM regions");
    const Case tb = eclipseDeckFromString(eclipseDeckToString(t));
    check(tb.thpres.size() == 1 && tb.thpres[0].value < 0 && eclipseDeckToString(t).find("EQLOPTS") != std::string::npos, "THPRES: deck round trip (EQLOPTS + THPRES)");
    // SWATINIT with a non-zero Pc curve
    auto setPc = [](std::vector<SatRow>& tb_) { if (tb_.empty()) return; const double s0 = tb_.front().s, s1 = tb_.back().s;
        for (auto& r : tb_) { const double se = (r.s - s0) / std::max(1e-9, s1 - s0); r.pc = 3.0 * std::pow(1 - se, 3); } };
    Case w = t; setPc(w.props.swof); for (auto& r : w.props.satExtra) setPc(r.swof);
    w.grid.swatinit.assign(size_t(w.grid.size()), 0.35);
    double dp2, ds2;
    const std::vector<float> sw0 = maxChange(w, dp2, ds2);
    float smin = 1; for (float v : sw0) smin = std::min(smin, v);
    std::printf("  SWATINIT 0.35: initial min Sw %.4f, max dSw after 6 steps %.2e\n", smin, ds2);
    check(std::fabs(smin - 0.35) < 1e-4 && ds2 < 2e-3, "SWATINIT: initial water saturation honoured and in equilibrium (scaled Pc)");
    const Case wb = eclipseDeckFromString(eclipseDeckToString(w));
    check(wb.grid.swatinit.size() == w.grid.swatinit.size() && std::fabs(wb.grid.swatinit[0] - 0.35) < 1e-12, "SWATINIT: deck round trip");
    // true-IMPES vs quasi-IMPES CPR give the same answer
    double f[2];
    for (int m = 0; m < 2; ++m) {
        Case q = regions(true); q.tuning.linearSolver = 1; q.tuning.cprWeights = m;
        Simulator s; s.log = [](const std::string&) {}; s.setup(q); s.initialize();
        for (int k = 0; k < 8 && s.runReportStep(); ++k) {}
        f[m] = s.summary().back().fpr;
    }
    check(std::fabs(f[0] - f[1]) < 1e-3 * f[0], "CPR quasi-IMPES and true-IMPES weights give the same solution");
    Case q = makeExample("spe1"); q.tuning.cprWeights = 0;
    check(eclipseDeckFromString(eclipseDeckToString(q)).tuning.cprWeights == 0 && projectFromString(projectToString(q)).tuning.cprWeights == 0,
          "CPR weight option is kept in decks and projects");
}

static void testWetGas() {
    std::printf("Wet gas / condensate (VAPOIL, PVTG, RVVD, Rv state switching)\n");
    const Case c = makeExample("condensate");
    // deck round trip
    const std::string d1 = eclipseDeckToString(c);
    check(d1.find("VAPOIL") != std::string::npos && d1.find("PVTG") != std::string::npos && d1.find("RVVD") != std::string::npos,
          "VAPOIL deck export writes VAPOIL, PVTG and RVVD");
    const Case b = eclipseDeckFromString(d1);
    bool same = b.phases.vapoil && b.props.pvtg.size() == c.props.pvtg.size();
    for (size_t k = 0; same && k < c.props.pvtg.size(); ++k) {
        same = b.props.pvtg[k].rows.size() == c.props.pvtg[k].rows.size() && std::fabs(b.props.pvtg[k].pg - c.props.pvtg[k].pg) < 1e-9;
        for (size_t j = 0; same && j < c.props.pvtg[k].rows.size(); ++j)
            same = std::fabs(b.props.pvtg[k].rows[j].rv - c.props.pvtg[k].rows[j].rv) < 1e-9 && std::fabs(b.props.pvtg[k].rows[j].B - c.props.pvtg[k].rows[j].B) < 1e-9;
    }
    check(same, "PVTG survives the deck round trip");
    check(eclipseDeckToString(b) == eclipseDeckToString(eclipseDeckFromString(eclipseDeckToString(b))), "VAPOIL deck import -> export is stable");
    saveProject(c, "wetgas_test.rsim");
    const Case pj = loadProject("wetgas_test.rsim");
    std::remove("wetgas_test.rsim");
    check(pj.phases.vapoil && pj.props.pvtg.size() == c.props.pvtg.size() && pj.props.pvtg.back().rows.size() == c.props.pvtg.back().rows.size(),
          "VAPOIL / PVTG survive the project round trip");
    // initial state: single-phase gas above the dew point, all oil vaporised
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    const InPlace ip = s.inPlace();
    std::printf("  initial liquid oil %.4g, vaporised oil %.4g, free gas %.4g (Rv %.4f)\n", ip.oil, ip.oilVaporized, ip.gasFree, ip.oilVaporized / ip.gasFree);
    check(ip.oil < 1e-6 * ip.oilVaporized && std::fabs(ip.oilVaporized / ip.gasFree - 0.16705) < 2e-3, "condensate: initial oil is all vaporised at the saturated Rv");
    double fopt = 0, maxSo = 0, minRv = 1e9, totalEnd = 0;
    while (s.runReportStep()) {}
    fopt = s.summary().back().fopt;
    const InPlace ie = s.inPlace(); totalEnd = ie.oil + ie.oilVaporized;
    const auto& sn = s.snapshots();
    for (auto& z : sn) for (size_t a = 0; a < z.so.size(); ++a) maxSo = std::max(maxSo, double(z.so[a]));
    for (size_t a = 0; a < sn[120].rv.size(); ++a) minRv = std::min(minRv, double(sn[120].rv[a]));
    std::printf("  FOPT %.6g (OPM Flow 2.52567e6), max SOIL %.3f, min RV after cycling %.4f, oil balance error %.2e\n", fopt, maxSo, minRv,
                std::fabs(ip.oilVaporized + ip.oil - totalEnd - fopt) / (ip.oil + ip.oilVaporized));
    check(maxSo > 0.05, "condensate: liquid drops out below the dew point (Rv -> Sg/So state switch)");
    check(minRv < 0.05, "condensate: lean injection gas strips the vaporised oil (undersaturated Rv)");
    check(std::fabs(ip.oilVaporized + ip.oil - totalEnd - fopt) < 1e-4 * (ip.oil + ip.oilVaporized), "condensate: oil component mass balance");
    check(std::fabs(fopt - 2.52567e6) < 0.01 * 2.52567e6, "condensate: cumulative oil within 1 % of OPM Flow (SPE3)");
    // enumerated initial RV round trip
    Case e = c; const int n = e.grid.size();
    e.init.pressure.assign(size_t(n), 3500); e.init.swat.assign(size_t(n), 0.2); e.init.sgas.assign(size_t(n), 0.8);
    e.init.rs.assign(size_t(n), 0.0); e.init.rv.assign(size_t(n), 0.12);
    const Case eb = eclipseDeckFromString(eclipseDeckToString(e));
    check(eb.init.rv.size() == size_t(n) && std::fabs(eb.init.rv[0] - 0.12) < 1e-9, "enumerated RV array: deck round trip");
    Simulator s2; s2.log = [](const std::string&) {}; s2.setup(eb); s2.initialize();
    const InPlace i2 = s2.inPlace();
    check(std::fabs(i2.oilVaporized / i2.gasFree - 0.12) < 1e-6, "enumerated RV initialises the vaporised oil");
}

static void testPassThrough() {
    std::printf("Keyword pass-through for OPM Flow (groups, VFP, MSW, ACTIONX, UDQ, summary requests)\n");
    std::string d = eclipseDeckToString(makeExample("spe1"));
    auto ins = [&](const std::string& after, const std::string& text) {
        const size_t p = d.find(after); if (p == std::string::npos) { check(false, ("marker " + after).c_str()); return; }
        d.insert(p + after.size(), text);
    };
    ins("RUNSPEC\n\n", "VFPPDIMS\n 5 5 5 5 0 1 /\nWSEGDIMS\n 1 10 2 /\n");
    ins("SUMMARY\n\n", "GOPR\n 'FIELD' /\nFU_X\n");
    ins("SCHEDULE\n\n", "GRUPTREE\n 'PROD' 'FIELD' /\n 'INJ' 'FIELD' /\n/\n"
        "GCONPROD\n 'FIELD' ORAT 15000 /\n/\n"
        "VFPPROD\n 1 8400 LIQ WCT GOR THP ' ' FIELD BHP /\n 100 1000 /\n 100 500 /\n 0 /\n 1 /\n 0 /\n 1 1 1 1 1500 1600 /\n 2 1 1 1 2500 2600 /\n"
        "UDQ\n DEFINE FU_X FOPR * 2 /\n/\n"
        "WELSEGS\n 'PROD' 8400 8400 1* INC /\n 2 2 1 1 100 100 0.5 1e-4 /\n/\n"
        "COMPSEGS\n 'PROD' /\n 10 10 3 1 0 100 /\n/\n"
        "ACTIONX\n 'A1' 1 /\n FOPR < 100 /\n/\nWELOPEN\n 'PROD' SHUT /\n/\nENDACTIO\n");
    const size_t sc = d.find("SCHEDULE\n");
    const size_t tst = d.find("TSTEP", sc);
    d.insert(tst, "WCONPROD\n 'PROD' OPEN GRUP 1* 1* 1* 1* 1* 1000 50 1 /\n/\nWELOPEN\n 'PROD' SHUT 10 10 3 /\n/\n");
    std::vector<DeckMessage> rep;
    const Case c = eclipseDeckFromString(d, "", nullptr, &rep);
    auto has = [&](const char* n, DeckSection s) { for (auto& e : c.extraKeywords) if (e.name == n && e.section == s) return true; return false; };
    check(has("GRUPTREE", DeckSection::Schedule) && has("GCONPROD", DeckSection::Schedule) && has("VFPPROD", DeckSection::Schedule) &&
          has("UDQ", DeckSection::Schedule) && has("WELSEGS", DeckSection::Schedule) && has("COMPSEGS", DeckSection::Schedule) &&
          has("ACTIONX", DeckSection::Schedule) && has("VFPPDIMS", DeckSection::Runspec) && has("WSEGDIMS", DeckSection::Runspec) &&
          has("GOPR", DeckSection::Summary) && has("FU_X", DeckSection::Summary),
          "pass-through: unmodelled keywords are kept with their section");
    const std::string e1 = eclipseDeckToString(c);
    auto pos = [&](const char* s) { return e1.find(s); };
    check(pos("VFPPDIMS") < pos("\nGRID") && pos("GOPR") > pos("\nSUMMARY") && pos("GOPR") < pos("\nSCHEDULE") &&
          pos("GRUPTREE") > pos("\nSCHEDULE") && pos("VFPPROD") < pos("WCONPROD") && pos("ENDACTIO") != std::string::npos,
          "pass-through: written back in the right sections and order (VFP tables before the controls)");
    check(e1.find("GRUP 1* 1* 1* 1* 1* 1000 50 1") != std::string::npos, "pass-through: WCONPROD GRUP mode, THP limit and VFP table kept");
    bool shut = false;
    for (auto& w : c.wells) if (w.name == "PROD") for (auto& q : w.completions) if (q.k1 == 2 && q.k2 == 2 && (!q.open || q.untilStep >= 0)) shut = true;
    check(shut, "connection-level WELOPEN (I J K) shuts the connection natively");
    check(eclipseDeckToString(eclipseDeckFromString(e1)) == e1, "pass-through: deck import -> export is stable");
    const std::string ud = e1.substr(pos("UDQ"), 40);
    check(ud.find("DEFINE FU_X FOPR * 2") != std::string::npos, "pass-through: UDQ expressions keep their text");
}

static void testEndPointScaling() {
    std::printf("End-point scaling (ENDSCALE, SWL/SWCR/SWU/SGL/SGCR/SGU/SOWCR/SOGCR)\n");
    auto run = [](const Case& c, double* fwpt, double* wip0) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        if (wip0) *wip0 = s.inPlace().water;
        while (s.runReportStep()) {}
        if (fwpt) *fwpt = s.summary().back().fwpt;
        return s.summary().back().fopt;
    };
    const Case base = makeExample("qfs");
    double w0 = 0, wi0 = 0;
    const double f0 = run(base, &w0, &wi0);
    SatFunctions sf; sf.setup(base, 0);
    const SatFunctions::EndPoints T = sf.tableEndPoints();
    std::printf("  table end points: SWL %.3f SWCR %.3f SWU %.3f SOWCR %.3f\n", T.swl, T.swcr, T.swu, T.sowcr);
    const size_t n = size_t(base.grid.size());
    // identity scaling (cell end points = table end points) reproduces the unscaled run
    Case id = base; id.grid.endscale = true;
    id.grid.swl.assign(n, T.swl); id.grid.swcr.assign(n, T.swcr); id.grid.swu.assign(n, T.swu); id.grid.sowcr.assign(n, T.sowcr);
    double w1 = 0;
    const double f1 = run(id, &w1, nullptr);
    std::printf("  FOPT unscaled %.8g, identity scaling %.8g\n", f0, f1);
    check(std::fabs(f1 - f0) < 1e-7 * f0 && std::fabs(w1 - w0) <= 1e-7 * std::max(1.0, w0), "EPS: identity end points reproduce the unscaled run");
    // larger critical water saturation: later breakthrough, less water produced
    Case sc = base; sc.grid.endscale = true;
    sc.grid.swcr.assign(n, std::min(T.swu - 0.2, T.swcr + 0.12));
    double w2 = 0;
    const double f2 = run(sc, &w2, nullptr);
    std::printf("  SWCR +0.12: FOPT %.6g, FWPT %.6g (unscaled %.6g)\n", f2, w2, w0);
    check(w2 < w0, "EPS: a higher SWCR delays water and reduces produced water");
    // larger residual oil: less oil recovered
    Case so = base; so.grid.endscale = true;
    so.grid.sowcr.assign(n, T.sowcr + 0.1);
    const double f3 = run(so, nullptr, nullptr);
    std::printf("  SOWCR +0.10: FOPT %.6g\n", f3);
    check(f3 < f0, "EPS: a higher SOWCR leaves more oil behind");
    // SWL raises the connate water: more initial water in place (equilibration with the scaled Pc)
    Case wl = base; wl.grid.endscale = true;
    wl.grid.swl.assign(n, T.swl + 0.08);
    wl.grid.swcr.assign(n, std::max(T.swcr, T.swl + 0.08));
    Simulator s; s.log = [](const std::string&) {}; s.setup(wl); s.initialize();
    const double wiL = s.inPlace().water;
    std::printf("  initial water in place: unscaled %.6g, SWL +0.08 %.6g\n", wi0, wiL);
    check(wiL > wi0 * 1.05, "EPS: SWL raises the initial (connate) water in equilibrium");
    // ENDSCALE deck round trip; arrays without ENDSCALE are ignored
    Case rt = sc; rt.grid.swcr[0] = -1;   // a defaulted cell is written with the table value
    const std::string d = eclipseDeckToString(rt);
    check(d.find("ENDSCALE") != std::string::npos && d.find("\nSWCR\n") != std::string::npos, "EPS: deck export writes ENDSCALE and SWCR");
    const Case b = eclipseDeckFromString(d);
    check(b.grid.endscale && b.grid.swcr.size() == n && std::fabs(b.grid.swcr[1] - sc.grid.swcr[1]) < 1e-9 && std::fabs(b.grid.swcr[0] - T.swcr) < 1e-9,
          "EPS: deck round trip (defaulted cell -> table end point)");
    std::string d2 = d;
    const size_t e = d2.find("ENDSCALE");
    d2.erase(e, d2.find('/', e) + 2 - e);
    const Case b2 = eclipseDeckFromString(d2);
    check(!b2.grid.hasEndPointScaling(), "EPS: SWCR without ENDSCALE is ignored (as ECLIPSE / OPM Flow)");
}

static std::string findDeck(const std::string& rel) {
    for (const char* pre : {"", "../", "../../", "ressim/"}) {
        const std::string f = std::string(pre) + rel;
        if (FILE* fp = std::fopen(f.c_str(), "rb")) { std::fclose(fp); return f; }
    }
    return {};
}

static void testEquilAndWellbore() {
    std::printf("Equilibration phase pressures and OPM-style wellbore connection densities\n");
    // gas cap: the oil pressure of an oil-free cell follows the gas curve (Po = Pg - Pcgo), water zone: Po = Pw + Pcow
    {
        const Case c = makeExample("gascap");
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        const auto& sn = s.snapshots().back();
        const Grid& g = s.grid();
        double worstGas = 0, worstWat = 0; int nGas = 0, nWat = 0;
        const double grav = c.units == UnitSystem::Field ? 0.00694444 : 9.80665e-5;
        for (int a = 0; a < g.numActive(); ++a) {
            int i, j, k; g.ijk(g.globalOf[a], i, j, k);
            if (k + 1 >= g.nz) continue;
            const int gb = i + g.nx * (j + g.ny * (k + 1));
            if (g.activeOf[gb] < 0) continue;
            const int b = g.activeOf[gb];
            const double dz = g.depth[b] - g.depth[a];
            if (dz <= 0) continue;
            const double grad = (sn.pressure[b] - sn.pressure[a]) / (dz * grav);   // apparent density
            if (sn.so[a] < 1e-9 && sn.so[b] < 1e-9 && sn.sg[a] > 0.3 && sn.sg[b] > 0.3) { worstGas = std::max(worstGas, grad); ++nGas; }
            if (sn.sw[a] > 1 - 1e-9 && sn.sw[b] > 1 - 1e-9) { worstWat = std::max(worstWat, std::fabs(grad - 62.4) / 62.4); ++nWat; }
        }
        std::printf("  gas-cap pairs %d: max apparent density %.2f; water-zone pairs %d: max deviation from water gradient %.3f\n", nGas, worstGas, nWat, worstWat);
        check(nGas == 0 || worstGas < 25.0, "equilibrium: gas-cap pressure follows the gas gradient (not the oil curve)");
        check(nWat == 0 || worstWat < 0.1, "equilibrium: water-zone pressure follows the water gradient");
    }
    // POLYMER-03 (opm-tests): wet gas, gas cap producer, polymer flood; OPM Flow 2024.10 reference totals
    const std::string d3 = findDeck("decks/polymer/POLYMER-03.DATA");
    if (d3.empty()) { std::printf("  [SKIP] decks/polymer/POLYMER-03.DATA not found\n"); return; }
    Case c = importEclipseDeck(d3);
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    while (s.runReportStep()) {}
    const auto& e = s.summary().back();
    std::printf("  POLYMER-03: FOPT %.6g (OPM 1.65343e6), FGPT %.6g (OPM 4.87754e8), FWPT %.6g (OPM 7.24841e6)\n", e.fopt, e.fgpt, e.fwpt);
    check(std::fabs(e.fopt / 1.65343e6 - 1) < 0.01, "POLYMER-03: cumulative oil within 1 % of OPM Flow (connection densities, gas-cap equilibrium)");
    check(std::fabs(e.fgpt / 4.87754e8 - 1) < 0.01, "POLYMER-03: cumulative gas within 1 % of OPM Flow");
    check(std::fabs(e.fwpt / 7.24841e6 - 1) < 0.01, "POLYMER-03: cumulative water within 1 % of OPM Flow");
}

static void testStoneAndDrsdt() {
    std::printf("Three-phase oil kr (STONE1 / STONE2) and DRSDT\n");
    {
        Case c = makeExample("gascap");
        c.props.kroModel = 2;
        SatFunctions sf; sf.setup(c, 0);
        const double swc = sf.swco(), sg = 0.2;
        const double k2 = sf.kro(AD3(swc), AD3(sg)).v, krog = sf.kro(AD3(swc + 1e-12), AD3(sg)).v;
        c.props.kroModel = 0; SatFunctions d; d.setup(c, 0);
        check(std::fabs(k2 - d.kro(AD3(swc), AD3(sg)).v) < 1e-9 && std::fabs(krog - k2) < 1e-6, "STONE2 reduces to krog at connate water");
        c.props.kroModel = 1; c.props.stone1Exp = 1.5;
        const std::string deck = eclipseDeckToString(c);
        const Case r = eclipseDeckFromString(deck, "");
        check(r.props.kroModel == 1 && std::fabs(r.props.stone1Exp - 1.5) < 1e-12, "STONE1 / STONE1EX survive the deck round trip");
        c.props.kroModel = 2;
        check(eclipseDeckFromString(eclipseDeckToString(c), "").props.kroModel == 2, "STONE2 survives the deck round trip");
    }
    const double opmFopt[2] = {9.81402e6, 8.73187e6};   // OPM Flow 2024.10: GASCAP example with STONE1 / STONE2
    for (int m = 1; m <= 2; ++m) {
        Case c = makeExample("gascap"); c.props.kroModel = m;
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        const double f = s.summary().back().fopt;
        std::printf("  GASCAP STONE%d: FOPT %.6g (OPM %.6g)\n", m, f, opmFopt[m - 1]);
        check(std::fabs(f / opmFopt[m - 1] - 1) < 0.01, m == 1 ? "STONE1: cumulative oil within 1 % of OPM Flow" : "STONE2: cumulative oil within 1 % of OPM Flow");
    }
    {
        Case c = makeExample("spe1");
        c.schedule[0].drsdt = 0;
        const Case r = eclipseDeckFromString(eclipseDeckToString(c), "");
        check(!r.schedule.empty() && r.schedule[0].drsdt == 0 && !r.schedule[0].drsdtFree, "DRSDT survives the deck round trip");
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        const double fpr = s.summary().back().fpr;
        std::printf("  SPE1 + DRSDT 0: final FPR %.6g (OPM 3725.51, without DRSDT ~3832)\n", fpr);
        check(std::fabs(fpr / 3725.51 - 1) < 0.002, "SPE1 with DRSDT 0: final pressure within 0.2 % of OPM Flow");
    }
}

static void testStudyAdvanced() {
    std::printf("Studies: new parameters, multi-objective terms, Pareto front, experimental designs, proxy model, external state\n");
    Case c = smallCase();
    // parameters
    int inj = -1; for (int w = 0; w < int(c.wells.size()); ++w) if (c.wells[w].type != WellType::Producer) inj = w;
    double q0 = -1; for (auto& st : c.schedule) for (auto& wc : st.controls) if (inj >= 0 && wc.well == c.wells[inj].name && wc.rate > 0) { q0 = wc.rate; break; }
    std::vector<StudyParam> P(3);
    P[0].kind = ParamKind::InjRate; P[0].target = inj; P[0].base = 1;
    P[1].kind = ParamKind::Ntg; P[1].k1 = P[1].k2 = 0; P[1].base = 1;
    P[2].kind = ParamKind::OilVisc; P[2].base = 1;
    Case c2 = applyParams(c, P, {2.0, 0.5, 1.5});
    double q1 = -1; for (auto& st : c2.schedule) for (auto& wc : st.controls) if (inj >= 0 && wc.well == c2.wells[inj].name && wc.rate > 0) { q1 = wc.rate; break; }
    const double mu0 = !c.props.pvto.empty() ? c.props.pvto[0].mu : c.props.pvdo.empty() ? 0 : c.props.pvdo[0].mu;
    const double mu1 = !c2.props.pvto.empty() ? c2.props.pvto[0].mu : c2.props.pvdo.empty() ? 0 : c2.props.pvdo[0].mu;
    const size_t nc = size_t(c.grid.nx) * c.grid.ny * c.grid.nz;
    check(q0 > 0 && std::fabs(q1 / q0 - 2) < 1e-12 && c2.grid.ntg.size() == nc && std::fabs(c2.grid.ntg[0] - 0.5) < 1e-12 &&
          std::fabs(c2.grid.ntg[nc - 1] - 1.0) < 1e-12 && mu0 > 0 && std::fabs(mu1 / mu0 - 1.5) < 1e-12,
          "InjRate / NTG (layer box) / oil viscosity parameters are applied");
    // Pareto
    const auto pf = paretoFront({{1, 3}, {2, 2}, {3, 1}, {3, 3}, {NAN, 0}});
    check(pf == std::vector<int>({0, 1, 2}), "Pareto front of a 2-objective set");
    // designs
    std::vector<StudyParam> D(3);
    for (auto& d : D) { d.kind = ParamKind::PermXY; d.low = 0; d.high = 1; d.base = 0.5; d.log = false; }
    const size_t ff = designPoints(D, DesignKind::FullFactorial2, 0, 1).size(), cc = designPoints(D, DesignKind::CentralComposite, 0, 1).size(),
                 pb = designPoints(D, DesignKind::PlackettBurman, 0, 1).size(), lh = designPoints(D, DesignKind::LatinHypercube, 20, 1).size();
    std::printf("  designs: full factorial %zu, CCD %zu, Plackett-Burman %zu, LHS %zu\n", ff, cc, pb, lh);
    check(ff == 9 && cc == 15 && pb == 9 && lh == 20, "experimental design sizes (3 parameters)");
    // proxy on an analytic response
    const auto X = designPoints(D, DesignKind::LatinHypercube, 40, 3);
    std::vector<double> y;
    for (auto& v : X) y.push_back(2 + 3 * v[0] - v[1] + 0.5 * v[0] * v[0] + 0.2 * v[0] * v[2]);
    const ProxyModel m = fitProxy(D, X, y, true);
    double pred = 0; const auto opt = proxyOptimum(m, D, true, 1, &pred);
    std::printf("  proxy: R2 %.6f Q2 %.6f, optimum (%.3f %.3f %.3f) -> %.4f\n", m.r2, m.q2, opt[0], opt[1], opt[2], pred);
    check(m.ok && m.quadratic && m.r2 > 0.999999 && m.q2 > 0.9999 && opt[0] > 0.99 && opt[1] < 0.01 && std::fabs(pred - 5.7) < 1e-3,
          "quadratic proxy reproduces an analytic response and finds its maximum");
    // weights + terms
    Case h = c;
    h.observations = {{ObsQty::FPR, "", 100, 3000, -1, 1}, {ObsQty::FPR, "", 200, 2900, -1, 0}, {ObsQty::FOPR, "", 100, 500, -1, 1}};
    const auto sig = autoSigma(h.observations);
    const auto terms = objectiveTerms(h.observations, ObjectiveGrouping::PhaseClass);
    check(sig[1] > 1e29 && terms.size() == 2 && terms[0].idx.size() == 2, "observation weight 0 excludes it; objective terms by class");
    // external state -> same streamline fluxes as the simulator's own restored report
    Simulator A; A.setup(c); A.initialize(); while (A.runReportStep()) {}
    const int si = std::min(2, int(A.snapshots().size()) - 1);
    A.restoreReport(si);
    const StreamlineResult SA = traceStreamlines(A);
    std::vector<double> bhp;
    const SummaryPoint* sp = nullptr; for (auto& s2 : A.summary()) if (s2.time <= A.snapshots()[si].time + 1e-9) sp = &s2;
    for (size_t w = 0; w < A.wellNames().size(); ++w) bhp.push_back(sp && w < sp->wells.size() && sp->wells[w].open ? sp->wells[w].bhp : 0.0);
    Simulator B; B.setup(c); B.initialize();
    B.restoreExternal(A.snapshots()[si], si - 1, bhp);
    const StreamlineResult SB = traceStreamlines(B);
    double fa = 0, fb = 0; for (auto& L : SA.lines) fa += L.flux; for (auto& L : SB.lines) fb += L.flux;
    std::printf("  streamline flux: own state %.2f, external state %.2f\n", fa, fb);
    check(fa > 0 && std::fabs(fb / fa - 1) < 0.02, "restoreExternal (restart + BHP) reproduces the streamline flux");
}

static void testSPE10() {
    std::printf("SPE10 model 1 (section headers with trailing text, near-incompressible PVDG, injector at its BHP limit)\n");
    const std::string d = findDeck("decks/spe10/SPE10_MODEL1.DATA");
    if (d.empty()) { std::printf("  (deck not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(d, nullptr);
    bool propsClean = true;
    for (const auto& e : c.extraKeywords) if (e.section == DeckSection::Props) propsClean = false;
    check(propsClean, "'SUMMARY =====' is recognised as a section (no SUMMARY keywords in PROPS)");
    const Case r = eclipseDeckFromString(eclipseDeckToString(c), "");
    bool same = r.props.pvdg.size() == c.props.pvdg.size();
    for (size_t i = 0; same && i < c.props.pvdg.size(); ++i) same = r.props.pvdg[i].B == c.props.pvdg[i].B;
    bool decreasing = true;
    for (size_t i = 1; i < r.props.pvdg.size(); ++i) decreasing = decreasing && r.props.pvdg[i].B < r.props.pvdg[i - 1].B;
    check(same && decreasing, "PVDG Bg values differing in the 12th digit survive the export (OPM requires strictly decreasing Bg)");
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    bool ok = true;
    try { while (s.runReportStep()) {} } catch (const std::exception& e) { ok = false; std::printf("  %s\n", e.what()); }
    const auto& e = s.summary().back();
    std::printf("  FOPT %.6g (OPM 75660.9), FPR %.6g (OPM 1102.32), %d Newton\n", e.fopt, e.fpr, s.stats().newton);
    check(ok && std::fabs(e.fopt / 75660.9 - 1) < 0.01, "SPE10 model 1 runs to the end; cumulative oil within 1 % of OPM Flow");
    check(ok && std::fabs(e.fpr / 1102.32 - 1) < 0.01, "SPE10 model 1 final pressure within 1 % of OPM Flow");
}

static void testPolymer() {
    std::printf("Polymer flooding (POLYMER, Todd-Longstaff, adsorption, RRF, DPV, PLMIXNUM, PLYSHLOG)\n");
    {
        const Case e = makeExample("polymer");
        double fopt[2] = {0, 0};
        for (int pass = 0; pass < 2; ++pass) {
            Case d = e; if (pass) d.phases.polymer = false;
            Simulator s; s.log = [](const std::string&) {}; s.setup(d); s.initialize();
            while (s.runReportStep()) {}
            fopt[pass] = s.summary().back().fopt;
        }
        std::printf("  polymer example: FOPT %.6g (OPM Flow 545320), water flood %.6g (+%.1f %%)\n", fopt[0], fopt[1], 100 * (fopt[0] / fopt[1] - 1));
        check(std::fabs(fopt[0] / 545320 - 1) < 0.015, "polymer example: cumulative oil within 1.5 % of OPM Flow");
        check(fopt[0] > 1.03 * fopt[1], "polymer example: the polymer slug recovers more oil than the water flood");
    }
    const std::string d2 = findDeck("decks/polymer/POLYMER-02.DATA");
    if (d2.empty()) { std::printf("  [SKIP] decks/polymer/POLYMER-02.DATA not found\n"); return; }
    const Case c = importEclipseDeck(d2);
    check(c.phases.polymer && c.polymer.numMixRegions() == 2 && !c.polymer.plyshlog.empty() && !c.polymer.plyvisc.empty(),
          "POLYMER-02 import: POLYMER, 2 PLMIXNUM regions (PLYMAX/PLMIXPAR) and PLYSHLOG");
    const std::string t1 = eclipseDeckToString(c);
    const Case b = eclipseDeckFromString(t1);
    bool same = b.phases.polymer && b.polymer.numMixRegions() == 2 && b.polymer.plyshlog.size() == c.polymer.plyshlog.size() &&
                b.polymer.plyshlog[0].rows == c.polymer.plyshlog[0].rows && b.polymer.plyvisc == c.polymer.plyvisc &&
                b.polymer.plyads == c.polymer.plyads && b.polymer.mixing == c.polymer.mixing && b.polymer.cmax == c.polymer.cmax &&
                b.grid.plmixnum == c.grid.plmixnum;
    check(same, "polymer keywords survive the deck round trip (PLYVISC/PLYADS/PLYMAX/PLMIXPAR/PLYSHLOG/PLMIXNUM)");
    check(eclipseDeckToString(b) == t1, "polymer deck export is stable");
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    while (s.runReportStep()) {}
    const auto& e = s.summary().back();
    const double bal = e.fcit - e.fcpt - e.fcip - e.fcad;
    std::printf("  FOPT %.6g (OPM 20796.1), FCIT %.6g (OPM 290453), polymer balance error %.2e\n", e.fopt, e.fcit, std::fabs(bal) / std::max(1.0, e.fcit));
    check(std::fabs(bal) < 1e-4 * std::max(1.0, e.fcit), "polymer mass balance: injected - produced = in solution + adsorbed");
    check(std::fabs(e.fopt / 20796.1 - 1) < 0.015, "POLYMER-02: cumulative oil within 1.5 % of OPM Flow (shear thinning, 2 mixing regions)");
    check(std::fabs(e.fcit / 290453 - 1) < 0.015, "POLYMER-02: cumulative polymer injected within 1.5 % of OPM Flow (injectivity with shear)");
}

static void testPinch() {
    std::printf("== PINCH: pinch-out connections across thin / inactive layers\n");
    Case c = makeExample("xsec");
    GridInput& g = c.grid;
    const int nl = g.nx * g.ny;
    for (int q = 4 * nl; q < 5 * nl; ++q) g.dz[size_t(q)] = 1e-4;          // layer 5 pinches out
    g.pinch.on = true; g.pinch.threshold = 0.01;
    Grid G; G.build(c);
    int bridged = 0;
    for (auto& cn : G.conns) {
        const int ka = G.globalOf[cn.a] / nl, kb = G.globalOf[cn.b] / nl;
        if ((ka == 3 && kb == 5) || (ka == 5 && kb == 3)) ++bridged;
    }
    check(G.numPinchedCells == nl && G.numPinch == nl && bridged == nl, "thin layer made inactive and bridged by one NNC per column");
    Case d = c; d.grid.pinch.threshold = 0.0;                              // thin layer stays active, regular k-neighbours
    d.grid.actnum.assign(size_t(d.grid.size()), 1);
    for (int q = 6 * nl; q < 7 * nl; ++q) d.grid.actnum[size_t(q)] = 0;    // a thick inactive layer
    Grid D; D.build(d);
    d.grid.pinch.gap = false; Grid E; E.build(d);
    check(D.numPinch == nl && E.numPinch == 0, "GAP bridges thick inactive layers, NOGAP does not");
    const std::string s1 = eclipseDeckToString(c);
    const Case r = eclipseDeckFromString(s1);
    check(r.grid.pinch.on && std::fabs(r.grid.pinch.threshold - 0.01) < 1e-12 && eclipseDeckToString(r) == s1, "PINCH deck round trip");
    Simulator sim; const SummaryPoint e = runSteps(c, 3, sim);
    check(std::isfinite(e.fpr) && e.fopt > 0, "pinched model runs");
}

static void testObsFromSummary() {
    std::printf("== Observed data from a summary file (history / simulated vectors)\n");
    Case c = makeExample("qfs");
    c.schedule.resize(12);
    const RunResult r = runSingle(c);
    std::vector<std::string> names; for (auto& w : c.wells) names.push_back(w.name);
    ecl::Summary s = ecl::summaryFromRun(r.summary, names, c.units);
    s.startDay = c.startDay; s.startMonth = c.startMonth; s.startYear = c.startYear;
    SummaryObsOptions o; o.history = false; o.qtyMask = (1u << int(ObsQty::WOPR)) | (1u << int(ObsQty::WBHP)) | (1u << int(ObsQty::FPR));
    std::string note;
    const auto obs = observationsFromSummary(s, c, o, &note);
    double worst = 0; int n = 0;
    for (auto& ob : obs) {
        const double v = simulatedValue(r, c, ob.qty, ob.well, ob.time);
        if (std::isfinite(v)) { worst = std::max(worst, std::fabs(v - ob.value) / std::max(1.0, std::fabs(v))); ++n; }
    }
    std::printf("    %s; worst rel. diff %.2e\n", note.c_str(), worst);
    check(r.ok && n >= 24 && worst < 1e-5, "summary vectors -> observations reproduce the run");
    o.history = true;
    check(observationsFromSummary(s, c, o).empty(), "no history vectors in a ResSim summary -> no observations");
}

static void testDualPorosity() {
    std::printf("== Dual porosity / dual permeability (DUALPORO, DUALPERM, SIGMA, LTOSIGMA, NODPPM)\n");
    const Case c = makeExample("dualporo");
    check(validateCase(c).empty(), "dual porosity example validates");
    Grid G; G.build(c);
    const int nl = G.nx * G.ny;
    int mm = 0, mf = 0; double tmf = 0;
    for (auto& cn : G.conns) {
        const int ka = G.globalOf[cn.a] / nl, kb = G.globalOf[cn.b] / nl;
        if (cn.dir == 'M') { ++mf; if (G.globalOf[cn.a] == 0) tmf = cn.trans; }
        else if (ka < G.nzMatrix || kb < G.nzMatrix) ++mm;
    }
    check(mf == 675 && G.numMF == 675, "one matrix-fracture connection per block");
    check(mm == 0, "DUALPORO: no matrix-matrix or geometric matrix-fracture connections");
    const double cd = Units::of(c.units).cDarcy, V = 100.0 * 100.0 * 30.0;
    check(std::fabs(tmf - cd * V * c.grid.sigma * c.grid.permx[0]) < 1e-9 * tmf, "T_mf = CDARCY * V * sigma * K_matrix");
    check(std::fabs(G.effPerm(c.grid.permx, 3 * nl) - 4000.0 * 0.006) < 1e-9, "fracture permeability x fracture porosity (no NODPPM)");
    check(std::fabs(G.depth[G.activeOf[3 * nl]] - G.depth[G.activeOf[0]]) < 1e-9, "fracture cell takes the matrix geometry");
    {   // LTOSIGMA, anisotropic Kazemi form
        Case l = c; l.grid.sigma = 0; l.grid.ltosigma = true;
        l.grid.lx.assign(size_t(l.grid.size()), 12.0); l.grid.ly = l.grid.lx; l.grid.lz.assign(size_t(l.grid.size()), 8.0);
        Grid L; L.build(l);
        double t = 0; for (auto& cn : L.conns) if (cn.dir == 'M' && L.globalOf[cn.a] == 0) t = cn.trans;
        const double want = cd * V * (4 * 0.8 / 144 + 4 * 0.8 / 144 + 4 * 0.3 / 64);
        check(std::fabs(t - want) < 1e-9 * want, "LTOSIGMA: T = CDARCY V (fx Kx/lx^2 + fy Ky/ly^2 + fz Kz/lz^2)");
        Case p = c; p.grid.dual = 2; Grid P; P.build(p);
        int m2 = 0; for (auto& cn : P.conns) if (cn.dir != 'M' && P.globalOf[cn.b] / nl < P.nzMatrix) ++m2;
        check(m2 > 0, "DUALPERM: the matrix grid keeps its connections");
    }
    Simulator s1, s2, s3;
    Case lo = c; lo.grid.sigma *= 1e-4;
    const SummaryPoint base = runSteps(c, 16, s1), weak = runSteps(lo, 16, s2);
    Case dp = c; dp.grid.dual = 2; const SummaryPoint perm = runSteps(dp, 16, s3);
    std::printf("    FOPT base %.0f  sigma/1e4 %.0f  DUALPERM %.0f  (FWCT %.3f / %.3f)\n", base.fopt, weak.fopt, perm.fopt, base.fwct, weak.fwct);
    check(base.fopt > 3 * weak.fopt && weak.fwct > base.fwct, "imbibition from the matrix: weak coupling recovers far less, earlier water");
    check(perm.fopt >= base.fopt * 0.99, "DUALPERM recovers at least as much as DUALPORO");
    const std::string d = eclipseDeckToString(c);
    const Case r = eclipseDeckFromString(d);
    check(r.grid.dual == 1 && d.find("DUALPORO") != std::string::npos && d.find("SIGMA") != std::string::npos && eclipseDeckToString(r) == d,
          "DUALPORO / SIGMA deck round trip is exact");
}

static void testDeckRepair() {
    std::printf("== Deck repair: misplaced keywords, non-monotonic PVDG\n");
    Case c;
    for (auto& e : listExamples()) {
        Case x = makeExample(e.id);
        if (x.phases.gas && !x.phases.vapoil && x.props.pvdg.size() > 3) { c = x; break; }
    }
    check(c.props.pvdg.size() > 3, "found a dry-gas example");
    if (c.props.pvdg.size() <= 3) return;
    c.props.pvdg[2].B = c.props.pvdg[1].B;                       // flat Bg: OPM Flow rejects it
    std::string a = eclipseDeckToString(c);
    check(a.find("row(s) adjusted") != std::string::npos, "PVDG repair noted in the deck");
    // summary vectors placed in PROPS (no SUMMARY header in the source deck)
    const size_t pp = a.find("\nPROPS\n");
    check(pp != std::string::npos, "PROPS section found");
    a.insert(pp + 7, "\nRUNSUM\n\nWOPR\n/\n\nFOPT\n\n");
    std::vector<std::string> warn;
    const Case c2 = eclipseDeckFromString(a, "", &warn);
    bool dec = true;
    for (size_t i = 1; i < c2.props.pvdg.size(); ++i) dec = dec && c2.props.pvdg[i].B < c2.props.pvdg[i - 1].B;
    check(dec, "re-imported PVDG strictly decreasing");
    const std::string b = eclipseDeckToString(c2);
    const size_t ps = b.find("\nPROPS\n"), ss = b.find("\nSUMMARY\n");
    const size_t ru = b.find("\nRUNSUM\n");
    check(ps != std::string::npos && ss != std::string::npos && ru != std::string::npos && ru > ss,
          "RUNSUM written in SUMMARY, not PROPS");
    const std::string props = b.substr(ps, ss - ps);
    check(props.find("\nFOPT") == std::string::npos && props.find("\nWOPR") == std::string::npos, "no summary vectors left in PROPS");
}

// ECLIPSE decks with features ResSim cannot model must still open (lenient pass) and be marked external-only
static void testLenientImport() {
    std::printf("Lenient import / external-only decks\n");
    const std::string deck =
        "RUNSPEC\nDIMENS\n 3 3 1 /\nOIL\nWATER\nFIELD\nGRID\nDX\n 9*100 /\nDY\n 9*100 /\nDZ\n 9*20 /\nTOPS\n 9*5000 /\n"
        "PORO\n 9*0.2 /\nPERMX\n 9*100 /\nPERMY\n 9*100 /\nPERMZ\n 9*10 /\nPROPS\nDENSITY\n 50 62 0.06 /\nSOLUTION\n"
        "PRESSURE\n 9*3000 /\nSWAT\n 9*0.2 /\nSCHEDULE\nTSTEP\n 10 /\nEND\n";
    std::vector<std::string> warn;
    bool opened = false, refused = false;
    try {
        const Case c = eclipseDeckFromString(deck, "", &warn);
        opened = c.grid.nx == 3 && !c.externalOnly.empty();
        try {
            Simulator sim;
            sim.setup(c);
        } catch (const std::exception&) {
            refused = true;
        }
    } catch (const std::exception& e) {
        std::printf("  import threw: %s\n", e.what());
    }
    check(opened, "deck without PVT/SCAL tables opens leniently, marked external-only (OPM/ECLIPSE run it)");
    check(refused, "the ResSim engine refuses an external-only case with a clear message");
}

// Block solver (block ILU / CPR-AMG / FGMRES with Schur-folded wells) must give the same answer in every mode
static void testBlockSolver() {
    std::printf("Block linear solver\n");
    Case c = makeExample("spe1");
    c.schedule.resize(1);
    Simulator sim;
    sim.setup(c);
    sim.initialize();
    bool ok = true;
    for (int mode : {0, 1, 2}) {
        Case cm = c;
        cm.tuning.linearSolver = mode;
        Simulator s2;
        s2.setup(cm);
        s2.initialize();
        try {
            while (s2.runReportStep()) {}
        } catch (...) {
            ok = false;
        }
        if (s2.stats().newton <= 0) ok = false;
    }
    check(ok, "SPE1 report step runs with block ILU, CPR-AMG and automatic solver selection");

    // AD math: d/dx of pow(x, y) and smoothMax
    AD<2> x = AD<2>::variable(2.0, 0), y = AD<2>::variable(3.0, 1);
    const AD<2> f = pow(x, y);
    check(std::fabs(f.v - 8) < 1e-12 && std::fabs(f.d[0] - 12) < 1e-12 && std::fabs(f.d[1] - 8 * std::log(2.0)) < 1e-12,
          "AD pow(x, y): value and both partial derivatives exact");

    ADV a = ADV::variable(1.5, 4, 2), b = ADV::variable(0.5, 4, 3);
    const ADV g = exp(a * b) / (a + b);
    const double e = std::exp(0.75);
    check(std::fabs(g.d[2] - (0.5 * e * 2.0 - e) / 4.0) < 1e-12 &&
          std::fabs(g.d[3] - (1.5 * e * 2.0 - e) / 4.0) < 1e-12 &&
          g.d[0] == 0,
          "runtime-sized AD (compositional / thermal unknowns): exact derivatives");
}

static void testSummaryWriter() {
    std::printf("ECLIPSE binary summary writer\n");
    Case c = makeExample("qfs");
    Simulator sim;
    sim.setup(c);
    sim.initialize();
    while (sim.runReportStep()) {}
    std::vector<std::string> names;
    for (auto& w : c.wells) names.push_back(w.name);
    const ecl::Summary a = ecl::summaryFromRun(sim.summary(), names, c.units);
    std::string err;
    const bool wrote = ecl::writeSummary(a, "rt_smry", &err);
    ecl::Summary b;
    const bool read = wrote && ecl::readSummary("rt_smry.SMSPEC", b, &err);
    bool same = read && b.time.size() == a.time.size() && b.reportIndex.size() == a.reportIndex.size();
    for (const char* k : {"FOPT", "FPR", "FWCT"}) {
        const int ia = a.find(k), ib = b.find(k);
        if (ia < 0 || ib < 0) { same = false; continue; }
        const double va = a.values[size_t(ia)].back(), vb = b.values[size_t(ib)].back();
        if (std::fabs(va - vb) > 1e-6 * std::max(1.0, std::fabs(va))) same = false;
    }
    std::remove("rt_smry.SMSPEC");
    std::remove("rt_smry.UNSMRY");
    check(same, "SMSPEC/UNSMRY written by the own engine read back identically (vectors, times, report steps)");
}

int main(int argc, char** argv) {
    const bool verbose = argc > 1 && std::string(argv[1]) == "-v";
    try { testIO(); testDeck(); testCPR(); testLenientImport(); testBlockSolver(); testSummaryWriter(); testCornerPoint(); testNewtran(); testPhaseA(); testWetGas(); testEquilAndWellbore(); testPolymer(); testDualPorosity(); testPinch(); testDeckRepair(); testObsFromSummary(); testStoneAndDrsdt(); testSPE10(); testEndPointScaling(); testPassThrough(); testAquifers(); testExamples(); testStreamlines(); testStudy(); testStudyAdvanced(); testBL(); testCorrelations(); testDemo(verbose); }
    catch (const std::exception& ex) { std::printf("EXCEPTION: %s\n", ex.what()); return 2; }
    std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures);
    return failures ? 1 : 0;
}
