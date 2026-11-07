#include <type_traits>
#include <random>
#include <map>
#include "ressim/eor.h"
#include "ressim/gasdata.h"
#include "ressim/thermsim.h"
#include "ressim/therm_deck.h"
#include "ressim/steam.h"
#include "ressim/nonlinear.h"
#include "ressim/eos.h"
#include "ressim/compsim.h"
#include "ressim/comp_deck.h"
#include "ressim/deck_parser.h"
#include "ressim/bsr.h"
// Engine regression tests: material balance (3-phase demo) and Buckley-Leverett (1D waterflood).
#include "ressim/simulator.h"
#include "ressim/io.h"
#include "ressim/tools.h"
#include "ressim/study.h"
#include "ressim/components.h"
#include "ressim/streamline.h"
#include "ressim/examples.h"
#include "ressim/grid.h"
#include "ressim/ecl.h"
#include "ressim/results_file.h"
#include "ressim/binio.h"
#include "ressim/hydrofrac.h"
#include "ressim/faultstab.h"
#include "ressim/mpi_rt.h"
#include "ressim/mpi_dd.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <atomic>
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


// ---- multi-component product analysis (C9): species / heating values, stream analysis, component rate bookkeeping
static void testComponents() {
    std::printf("== Multi-component product analysis (C9)\n");
    using namespace ressim;
    check(recogniseSpecies("C1") == SP_CH4 && recogniseSpecies("nC4") == SP_C4H10 && recogniseSpecies(" h2 ") == SP_H2 &&
          recogniseSpecies("Water") == SP_H2O && recogniseSpecies("CO") == SP_CO && recogniseSpecies("C7+") == -1,
          "species recognition from CNAMES conventions");
    std::vector<ComponentInfo> C(7);
    const char* nm[] = {"H2", "CO", "CO2", "C1", "N2", "C20", "COAL"};
    for (int i = 0; i < 7; ++i) C[size_t(i)].name = nm[i];
    C[5].mw = 0.282; C[6].kind = 2; C[6].mw = 0.012;
    C.push_back(ComponentInfo{"H2O", 0.0180153, 3});
    assignHeatingValues(C);
    check(std::fabs(C[0].lhv - 241.82e3) < 1 && std::fabs(C[3].hhv - 890.36e3) < 1 && C[2].lhv == 0 && std::fabs(C[5].lhv - 44.0e6 * 0.282) < 1 &&
          C[6].lhv == 0 && C[7].lhv == 0 && std::fabs(C[0].mw - 2.01588e-3) < 1e-9, "standard heating values (H2, CH4, inert CO2, C5+ correlation, solid, water)");
    const std::vector<double> q = {1.0, 2.0, 0.5, 0.5, 1.0, 0.0, 0.0, 3.0}, inj = {0, 0, 1.0, 0, 0, 0, 0, 0};
    const StreamAnalysis A = analyseStream(C, q, &inj);
    const double lhv = 241.82e3 + 2 * 282.98e3 + 0.5 * 802.34e3;
    std::printf("  syngas %.3g mol/s, H2/CO %.3f, dry %.3g, y(H2) %.3f, LHV %.4g kW, gas LHV %.3f MJ/sm3, CO2 slip %.2f\n",
                A.syngas, A.h2co, A.dry, A.yDry[0], A.lhvRate * 1e-3, A.lhvGas, A.co2Slip);
    check(std::fabs(A.syngas - 3) < 1e-12 && std::fabs(A.h2co - 0.5) < 1e-12 && std::fabs(A.dry - 5) < 1e-12 && std::fabs(A.total - 8) < 1e-12 &&
          std::fabs(A.yDry[0] - 0.2) < 1e-12 && A.yDry[7] == 0 && std::fabs(A.lhvRate - lhv) < 1e-6 && std::fabs(A.co2Slip - 0.5) < 1e-12 &&
          std::fabs(A.lhvGas - lhv / 5 * kMolPerSm3 * 1e-6) < 1e-9, "stream analysis: syngas, H2/CO, dry composition, LHV power and MJ/sm3, CO2 slip");
    {   // incremental cumulatives
        SummaryPoint a, b; a.dt = 1; a.fcompP = {1, 2}; a.wells.resize(1); a.wells[0].comp = {1, -2};
        accumulateComponentPoint(nullptr, a);
        b = a; b.dt = 2; accumulateComponentPoint(&a, b);
        check(std::fabs(b.fcompPT[1] - 3 * 2 * 86400) < 1e-6 && std::fabs(b.wells[0].compCumP[0] - 3 * 86400) < 1e-6 &&
              std::fabs(b.wells[0].compCumI[1] - 6 * 86400) < 1e-6 && b.wells[0].compCumI[0] == 0, "component cumulatives (field and per well, production / injection split)");
    }
    {   // compositional example: report-step average well component rates integrate to the simulator's produced moles
        const Case c = makeExample("comp");
        comp::DeckModel m; std::string err; std::vector<comp::CompSummary> out;
        std::vector<double> sum, ref;
        double prevT = 0;
        const bool ok = comp::deckModel(c, m, &err) && comp::runDeckModel(m, out, [&](const comp::CompSummary& S, const comp::CompSimulator& sim) {
            const double dt = (S.time - prevT) * 86400; prevT = S.time;
            sum.resize(size_t(sim.ncw()), 0.0); ref.assign(size_t(sim.ncw()), 0.0);
            for (const auto& W : S.wells) for (size_t i = 0; i < W.molar.size() && i < sum.size(); ++i) sum[i] += std::max(0.0, W.molar[i]) * dt;
            for (int i = 0; i < sim.ncw(); ++i) ref[size_t(i)] = sim.produced(i);
            return true;
        }, &err, 4);
        double e = 0, mx = 0;
        for (size_t i = 0; i < ref.size(); ++i) { e = std::max(e, std::fabs(sum[i] - ref[i])); mx = std::max(mx, ref[i]); }
        std::printf("  compositional example: %zu components, produced %.4g kmol, |sum(rate dt) - produced| / max = %.1e\n", ref.size(), mx * 1e-3, mx > 0 ? e / mx : 1.0);
        check(ok && mx > 0 && e < 1e-9 * mx, "compositional: per-well component rates conserve the produced moles");
    }
    {   // UCG: per-component products from the thermal core (syngas, hydrogen) and the analysis of the producer stream
        const Case c = makeExample("ucg");
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out;
        std::vector<ComponentInfo> ci; std::vector<double> sum, ref; double prevT = 0;
        StreamAnalysis last;
        const bool ok = thermal::deckModel(c, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator& sim) {
            if (ci.empty()) {
                for (const auto& k : sim.model().comps) { ComponentInfo x; x.name = k.name; x.mw = k.mw; x.kind = k.type == 2 ? 2 : k.type == 1 ? 1 : 0; ci.push_back(x); }
                ci.push_back(ComponentInfo{"H2O", 0.0180153, 3});
                assignHeatingValues(ci);
            }
            const double dt = (S.time - prevT) * 86400; prevT = S.time;
            sum.resize(ci.size(), 0.0); ref.assign(ci.size(), 0.0);
            std::vector<double> fp(ci.size(), 0.0), fi(ci.size(), 0.0);
            for (const auto& W : S.wells) for (size_t i = 0; i < W.molar.size() && i < ci.size(); ++i) {
                sum[i] += std::max(0.0, W.molar[i]) * dt;
                if (W.molar[i] > 0) fp[i] += W.molar[i]; else fi[i] -= W.molar[i];
            }
            for (const auto& r : sim.wellReports()) for (size_t i = 0; i < r.cumMolP.size() && i < ref.size(); ++i) ref[i] += r.cumMolP[i];
            if (S.time > 0) last = analyseStream(ci, fp, &fi);
            return true;
        }, &err, 4);
        double e = 0, mx = 0;
        for (size_t i = 0; i < ref.size(); ++i) { e = std::max(e, std::fabs(sum[i] - ref[i])); mx = std::max(mx, ref[i]); }
        std::printf("  UCG products: H2 %.4g kmol/d (%.3f t/d), syngas %.4g kmol/d, H2/CO %.3f, gas LHV %.2f MJ/sm3, LHV power %.3g MW, O2 slip %.3f, mole balance %.1e\n",
                    last.h2 * 86.4, last.h2Mass * 86.4, last.syngas * 86.4, last.h2co, last.lhvGas, last.lhvRate * 1e-6, last.o2Slip, mx > 0 ? e / mx : 1.0);
        check(ok && last.h2 > 0.01 && last.co > 0.01 && last.lhvGas > 1 && last.lhvRate > 0 && e < 1e-9 * mx,
              "UCG: H2 / syngas identified in the produced stream, LHV of the product gas, component moles conserved");
    }
}

static void testAdjoint() {
    std::printf("== Adjoint / forward-sensitivity history-matching gradients (C9)\n");
    Case base = smallCase();
    base.tuning.cnvTol = 1e-8; base.tuning.mbTol = 1e-11; base.tuning.linTol = 1e-9;
    base.tuning.initialDt = 10; base.tuning.maxDt = 10; base.tuning.growth = 1.0;
    std::vector<StudyParam> P(5);
    P[0].kind = ParamKind::PermXY; P[0].base = 1; P[0].low = 0.3; P[0].high = 3; P[0].log = true;
    P[1].kind = ParamKind::KrwMax; P[1].base = 1; P[1].low = 0.4; P[1].high = 1.6; P[1].log = true;
    P[2].kind = ParamKind::WaterVisc; P[2].base = 1; P[2].low = 0.5; P[2].high = 2;
    P[3].kind = ParamKind::Poro; P[3].base = 1; P[3].low = 0.7; P[3].high = 1.3;
    P[4].kind = ParamKind::WellPI; P[4].base = 1; P[4].low = 0.5; P[4].high = 2; P[4].target = 0;
    for (auto& p : P) p.active = true;
    const std::vector<double> truthV = {1.8, 0.6, 1.2, 0.9, 1.3};
    const Case truthCase = applyParams(base, P, truthV);
    const RunResult truth = runSingle(truthCase);
    base.observations = syntheticObservations(truth, truthCase, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::FPR, ObsQty::WOPR, ObsQty::FGOR}, 0.01, 7);
    std::printf("  truth ok=%d, %zu observations, %d steps\n", truth.ok, base.observations.size(), truth.stats.timesteps);
    StudyCallbacks cb;
    cb.log = [](const std::string& s) { std::printf("%s\n", s.c_str()); };
    GradientCheck gc = misfitGradient(base, P, 7, 1e-3, cb);
    check(gc.ok, "gradient check runs (adjoint, forward sensitivities, finite differences)");
    if (!gc.ok) { std::printf("  error: %s\n", gc.error.c_str()); return; }
    std::printf("  misfit %.6g; adjoint %.2fs, forward %.2fs, FD %.2fs (%d runs)\n", gc.misfit, gc.secAdjoint, gc.secForward, gc.secFd, gc.fdRuns);
    double eAF = 0, eFD = 0, gmax = 0;
    for (size_t k = 0; k < gc.act.size(); ++k) gmax = std::max(gmax, std::fabs(gc.fd[k]));
    for (size_t k = 0; k < gc.act.size(); ++k) {
        std::printf("    %-12s adjoint %+.6e  forward %+.6e  FD %+.6e\n", paramLabel(P[size_t(gc.act[k])], base).c_str(), gc.adjoint[k], gc.forward[k], gc.fd[k]);
        eAF = std::max(eAF, std::fabs(gc.adjoint[k] - gc.forward[k]) / gmax);
        eFD = std::max(eFD, std::fabs(gc.adjoint[k] - gc.fd[k]) / gmax);
    }
    std::printf("  max |adjoint - forward| / |g| = %.2e, max |adjoint - FD| / |g| = %.2e, replay deviation %.1e, linear residual %.1e\n",
                eAF, eFD, gc.replayDev, gc.maxLinRes);
    check(gc.replayDev < 1e-6, "adjoint replay reproduces the forward run");
    check(eAF < 1e-6, "adjoint gradient == forward-sensitivity gradient (dual consistency)");
    check(eFD < 1e-2, "adjoint gradient matches central finite differences of full runs (incl. control-switch events)");
    {   // smooth horizon (no control / phase switching events): the discrete adjoint is exact to FD round-off
        Case sb = base; sb.schedule.resize(4);
        StudyCallbacks cs;
        GradientCheck g4 = misfitGradient(sb, P, 7, 1e-4, cs);
        double e4 = 0, ef4 = 0, gm = 0;
        for (size_t k = 0; k < g4.act.size() && g4.ok; ++k) gm = std::max(gm, std::fabs(g4.fd[k]));
        for (size_t k = 0; k < g4.act.size() && g4.ok; ++k) {
            e4 = std::max(e4, std::fabs(g4.adjoint[k] - g4.fd[k]) / gm);
            ef4 = std::max(ef4, std::fabs(g4.adjoint[k] - g4.forward[k]) / gm);
        }
        std::printf("  smooth horizon (40 steps): max |adjoint - FD| / |g| = %.2e, |adjoint - forward| = %.2e\n", e4, ef4);
        check(g4.ok && e4 < 1e-5 && ef4 < 1e-7, "adjoint gradient exact on a smooth horizon (well heads differentiated)");
    }

    // history matching with the analytic gradients (truth kxy 1.8, krw 0.6)
    for (int i = 2; i < 5; ++i) P[size_t(i)].active = false;
    Case hmBase = smallCase();
    const Case tc2 = applyParams(hmBase, P, {1.8, 0.6, 1.0, 1.0, 1.0});
    const RunResult tr2 = runSingle(tc2);
    hmBase.observations = syntheticObservations(tr2, tc2, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::FPR}, 0.01, 3);
    StudyCallbacks cq;
    auto t0 = std::chrono::steady_clock::now();
    AdjointHMOptions ao; ao.maxIter = 30;
    HMResult ad = runAdjointLBFGS(hmBase, P, ao, cq);
    std::printf("  adjoint L-BFGS: misfit %.4g -> %.4g, kxy %.3f krw %.3f, %d runs, %.1fs\n", ad.initialMisfit, ad.bestMisfit,
                ad.best.empty() ? 0.0 : ad.best[0], ad.best.empty() ? 0.0 : ad.best[1], ad.runs,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(ad.ok && ad.bestMisfit < 0.02 * ad.initialMisfit && std::fabs(ad.best[0] / 1.8 - 1) < 0.1 && std::fabs(ad.best[1] / 0.6 - 1) < 0.1,
          "adjoint L-BFGS recovers the truth parameters");
    {   // C13: iteration records for the dedicated adjoint HM view
        bool rec = ad.history.size() >= 2;
        for (size_t i = 0; i < ad.history.size() && rec; ++i) {
            const HMIteration& h = ad.history[i];
            rec = h.run.ok && !h.run.summary.empty() && h.grad.size() == P.size() && std::isfinite(h.gradNorm) && h.runs >= 1 &&
                  (i == 0 || (h.misfit <= ad.history[i - 1].misfit && h.runs >= ad.history[i - 1].runs && h.seconds >= ad.history[i - 1].seconds));
        }
        std::printf("  adjoint iteration records: %zu iterates, |g| %.3g -> %.3g\n", ad.history.size(), ad.history.empty() ? 0.0 : ad.history.front().gradNorm,
                    ad.history.empty() ? 0.0 : ad.history.back().gradNorm);
        check(rec && ad.history.back().gradNorm < ad.history.front().gradNorm, "adjoint L-BFGS records every iterate (run, gradient, runs, time); the gradient norm falls");
    }
    t0 = std::chrono::steady_clock::now();
    LMOptions lo; lo.maxIter = 10; lo.gradient = GradientMethod::Analytic;
    HMResult lm = runLevenbergMarquardt(hmBase, P, lo, cq);
    std::printf("  LM (analytic Jacobian): misfit %.4g -> %.4g, kxy %.3f krw %.3f, %d runs, %.1fs\n", lm.initialMisfit, lm.bestMisfit,
                lm.best[0], lm.best[1], lm.runs, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    check(lm.ok && lm.bestMisfit < 0.02 * lm.initialMisfit && std::fabs(lm.best[0] / 1.8 - 1) < 0.1 && std::fabs(lm.best[1] / 0.6 - 1) < 0.1,
          "Levenberg-Marquardt with the analytic Jacobian recovers the truth");
}

// C11: analytic field gradients (permeability / porosity / MULTZ / well PI): one reverse pass, no perturbed replays
static void testAdjointField() {
    std::printf("== C11: adjoint field gradients (analytic dR/dtheta for perm, poro, MULTZ, WI)\n");
    Case base = smallCase();
    base.tuning.cnvTol = 1e-8; base.tuning.mbTol = 1e-11; base.tuning.linTol = 1e-9;
    base.tuning.initialDt = 10; base.tuning.maxDt = 10; base.tuning.growth = 1.0;
    // heterogeneous so that the harmonic averages matter
    for (size_t q = 0; q < base.grid.permx.size(); ++q) base.grid.permx[q] = 60.0 + 25.0 * double((q * 7) % 13);
    // box parameters: PERM X/Y per quadrant, PERM Z, MULTZ, PORO per layer, PI of both wells
    std::vector<StudyParam> P;
    auto box = [&](ParamKind k, int i1, int i2, int j1, int j2, int k1, int k2, double lo, double hi) {
        StudyParam p; p.kind = k; p.i1 = i1; p.i2 = i2; p.j1 = j1; p.j2 = j2; p.k1 = k1; p.k2 = k2; p.base = 1; p.low = lo; p.high = hi; p.log = true; p.active = true; P.push_back(p);
    };
    box(ParamKind::PermXY, 0, 4, 0, 4, -1, -1, 0.3, 3); box(ParamKind::PermXY, 5, 9, 0, 4, -1, -1, 0.3, 3);
    box(ParamKind::PermXY, 0, 4, 5, 9, -1, -1, 0.3, 3); box(ParamKind::PermXY, 5, 9, 5, 9, -1, -1, 0.3, 3);
    box(ParamKind::PermZ, -1, -1, -1, -1, -1, -1, 0.1, 10); box(ParamKind::MultZ, 3, 9, 3, 9, 0, 0, 0.1, 10);
    box(ParamKind::Poro, -1, -1, -1, -1, 0, 0, 0.7, 1.3); box(ParamKind::Poro, -1, -1, -1, -1, 1, 1, 0.7, 1.3);
    for (int w = 0; w < 2; ++w) { StudyParam p; p.kind = ParamKind::WellPI; p.target = w; p.base = 1; p.low = 0.5; p.high = 2; p.active = true; P.push_back(p); }
    const std::vector<double> truthV = {1.6, 0.7, 1.2, 2.0, 0.5, 0.4, 0.9, 1.1, 1.4, 0.8};
    const Case truthCase = applyParams(base, P, truthV);
    const RunResult truth = runSingle(truthCase);
    base.observations = syntheticObservations(truth, truthCase, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::FPR, ObsQty::WOPR, ObsQty::FGOR, ObsQty::WWIR}, 0.01, 11);
    Case sb = base; sb.schedule.resize(4);   // smooth horizon: no control / phase switching between the perturbed runs
    StudyCallbacks cs;
    GradientCheck g = misfitGradient(sb, P, 7, 1e-4, cs);
    check(g.ok, "field-gradient check runs");
    if (!g.ok) { std::printf("  error: %s\n", g.error.c_str()); return; }
    double gm = 0, eAF = 0, eFD = 0;
    for (size_t k = 0; k < g.act.size(); ++k) gm = std::max(gm, std::fabs(g.fd[k]));
    for (size_t k = 0; k < g.act.size(); ++k) {
        std::printf("    %-22s adjoint(field) %+.6e  forward(replay) %+.6e  FD %+.6e\n", paramLabel(P[size_t(g.act[k])], sb).c_str(), g.adjoint[k], g.forward[k], g.fd[k]);
        eAF = std::max(eAF, std::fabs(g.adjoint[k] - g.forward[k]) / gm);
        eFD = std::max(eFD, std::fabs(g.adjoint[k] - g.fd[k]) / gm);
    }
    std::printf("  %zu parameters: adjoint %.2fs (analytic dR/dtheta), forward sensitivities %.2fs (perturbed replays), FD %.2fs\n",
                g.act.size(), g.secAdjoint, g.secForward, g.secFd);
    std::printf("  max |field adjoint - forward| / |g| = %.2e, max |field adjoint - FD| / |g| = %.2e\n", eAF, eFD);
    check(eAF < 1e-6, "field adjoint (analytic dR/dtheta) == forward sensitivities of perturbed models");
    check(eFD < 1e-5, "field adjoint matches central finite differences of full runs");

    // per-cell maps from one pass: the sum over all cells reproduces a global multiplier's gradient
    {
        Case c = sb;
        Simulator S; S.setAdjointRecording(true); S.setup(c); S.initialize();
        while (!S.finished()) if (!S.runReportStep()) break;
        const std::vector<double> sig = autoSigma(c.observations);
        std::vector<Simulator::AdjObs> ao;
        for (size_t i = 0; i < c.observations.size(); ++i) {
            const Observation& o = c.observations[i];
            Simulator::AdjObs a; a.qty = o.qty; a.time = o.time; a.value = o.value; a.sigma = sig[i];
            if (obsIsWell(o.qty)) {
                a.well = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == o.well) a.well = int(w);
                if (a.well < 0) continue;
            }
            ao.push_back(a);
        }
        // reference: global PERM X/Y and PORO multipliers by perturbed replays
        std::vector<StudyParam> G(2);
        G[0].kind = ParamKind::PermXY; G[1].kind = ParamKind::Poro;
        for (auto& p : G) { p.active = true; p.log = false; p.base = 1; p.low = 0.5; p.high = 1.5; }
        const double h = 1e-4;
        std::vector<Case> cp(2), cm(2);
        std::vector<Simulator::AdjParam> ap(2);
        for (int j = 0; j < 2; ++j) {
            std::vector<double> vp = {1.0, 1.0}, vm = {1.0, 1.0};
            vp[size_t(j)] += h; vm[size_t(j)] -= h;
            cp[size_t(j)] = applyParams(c, G, vp); cm[size_t(j)] = applyParams(c, G, vm);
            ap[size_t(j)].plus = &cp[size_t(j)]; ap[size_t(j)].minus = &cm[size_t(j)]; ap[size_t(j)].h = h;
        }
        Simulator::AdjFieldReq fr; fr.perm = fr.poro = fr.trans = fr.wi = true;
        auto t0 = std::chrono::steady_clock::now();
        const Simulator::AdjResult rf = S.adjointGradient(ao, {}, false, &fr);
        const double tf = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        t0 = std::chrono::steady_clock::now();
        const Simulator::AdjResult rr = S.adjointGradient(ao, ap, false);
        const double tr = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double sK = 0, sP = 0;
        for (size_t a = 0; a < rf.gLogKx.size(); ++a) { sK += rf.gLogKx[a] + rf.gLogKy[a]; sP += rf.gLogPoro[a]; }
        const double eK = std::fabs(sK - rr.grad[0]) / std::max(std::fabs(rr.grad[0]), 1e-30), eP = std::fabs(sP - rr.grad[1]) / std::max(std::fabs(rr.grad[1]), 1e-30);
        std::printf("  one pass: %zu cells x (kx, ky, kz, phi) + %zu faces + WI = %zu gradients in %.2fs (field contraction %.3fs); 2 replay parameters %.2fs\n",
                    rf.gLogKx.size(), rf.gLogTrans.size(), 4 * rf.gLogKx.size() + rf.gLogTrans.size() + 2 * 2, tf, rf.secField, tr);
        std::printf("  sum_cells dF/dln k = %.6e vs global multiplier %.6e (rel %.1e); sum_cells dF/dln phi = %.6e vs %.6e (rel %.1e)\n",
                    sK, rr.grad[0], eK, sP, rr.grad[1], eP);
        check(rf.ok && rr.ok && eK < 1e-6 && eP < 1e-6, "per-cell field gradients sum to the global multiplier gradients");
    }
    {   // per-cell field history matching: truth = a high-permeability channel along the diagonal, prior = uniform
        Case prior = smallCase();
        Case tc = prior;
        for (int k = 0; k < tc.grid.nz; ++k)
            for (int j = 0; j < tc.grid.ny; ++j)
                for (int i = 0; i < tc.grid.nx; ++i) {
                    const size_t q = size_t(tc.grid.index(i, j, k));
                    const double kx = std::abs(i - j) <= 1 ? 600.0 : 90.0;
                    tc.grid.permx[q] = kx; tc.grid.permz[q] = 0.1 * kx;
                }
        const RunResult tr = runSingle(tc);
        prior.observations = syntheticObservations(tr, tc, {ObsQty::WBHP, ObsQty::WWCT, ObsQty::WOPR, ObsQty::FPR}, 0.01, 5);
        FieldHMOptions fo; fo.maxIter = 15; fo.regWeight = 1e-3;
        StudyCallbacks cq;
        const auto t0 = std::chrono::steady_clock::now();
        const FieldHMResult fh = runAdjointFieldHM(prior, fo, cq);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double onC = 0, offC = 0; int nOn = 0, nOff = 0;
        for (int j = 0; j < prior.grid.ny && fh.ok; ++j)
            for (int i = 0; i < prior.grid.nx; ++i) {
                const double m = std::log(fh.multK[size_t(prior.grid.index(i, j, 0))]);
                if (std::abs(i - j) <= 1) { onC += m; ++nOn; } else { offC += m; ++nOff; }
            }
        if (nOn) onC /= nOn;
        if (nOff) offC /= nOff;
        std::printf("  field HM (%d unknowns): misfit %.4g -> %.4g in %d iterations / %d runs, %.1fs (forward %.1fs, adjoint %.1fs); mean ln-multiplier channel %+.2f, off-channel %+.2f\n",
                    fh.ok ? int(fh.multK.size()) : 0, fh.initialMisfit, fh.bestMisfit, fh.iterations, fh.runs, secs, fh.secForward, fh.secAdjoint, onC, offC);
        if (!fh.ok) std::printf("  error: %s\n", fh.error.c_str());
        check(fh.ok && fh.bestMisfit < 0.1 * fh.initialMisfit && onC > offC, "per-cell (field) adjoint history matching reduces the misfit > 10x and finds the channel");
        {   // C13: per-iteration records, snapshots and sensitivity maps of the dedicated field-HM view
            const size_t nc = size_t(prior.grid.size());
            bool rec = fh.ok && int(fh.iters.size()) == fh.iterations + 1 && fh.initialRun.ok && fh.gradK0.size() == nc && fh.gradK.size() == nc &&
                       fh.nx == prior.grid.nx && fh.nz == prior.grid.nz;
            for (size_t i = 0; i < fh.iters.size() && rec; ++i) {
                const auto& it = fh.iters[i];
                rec = it.lnK.size() == nc && it.run.ok && it.kMin <= it.kMean && it.kMean <= it.kMax && std::fabs(it.objective - it.misfit - it.reg) < 1e-9 * std::max(1.0, it.objective) &&
                      (i == 0 || it.objective <= fh.iters[i - 1].objective);
            }
            const auto& L = fh.iters.empty() ? FieldHMResult::Iter() : fh.iters.back();
            double e = 0; for (size_t q = 0; q < nc && q < L.lnK.size(); ++q) e = std::max(e, std::fabs(std::exp(double(L.lnK[q])) - fh.multK[q]) / fh.multK[q]);
            std::printf("  field HM records: %zu iterates, last snapshot vs final multipliers %.1e, |g| %.3g -> %.3g\n", fh.iters.size(), e,
                        fh.iters.empty() ? 0.0 : fh.iters.front().gradNorm, L.gradNorm);
            check(rec && e < 1e-5, "field HM records every iterate (objective = misfit + regularisation, multiplier snapshots, runs) and the sensitivity maps");
        }
    }
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
static void testMultRegtReopen() {
    std::printf("== MULTREGT on neighbour / fault-NNC connections, COMPDAT shut -> re-open intervals\n");
    Case cart = makeDemoCase();
    Case f = cart; addFaultThrow(f.grid, 10, 12.5, 1.0, "F1");
    const int nx = f.grid.nx, ny = f.grid.ny, nz = f.grid.nz;
    f.grid.fluxnum.assign(size_t(nx) * ny * nz, 1);
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 10; i < nx; ++i) f.grid.fluxnum[size_t(i + nx * (j + ny * k))] = 2;
    auto across = [](const Case& cc, double& tN, double& tNNC, size_t& n) {
        Grid g; g.build(cc); tN = tNNC = 0; n = g.conns.size();
        for (auto& c : g.conns) { int i1, j1, k1, i2, j2, k2; g.ijk(g.globalOf[c.a], i1, j1, k1); g.ijk(g.globalOf[c.b], i2, j2, k2);
            if (std::min(i1, i2) != 9 || std::max(i1, i2) != 10) continue;
            if (std::abs(g.globalOf[c.a] - g.globalOf[c.b]) == 1) tN += c.trans; else tNNC += c.trans; }
    };
    double n0, x0, n1, x1; size_t c0, c1;
    across(f, n0, x0, c0);
    Case a = f; a.grid.multregt.push_back({1, 2, 0.5, "XYZ", 0, 'F'});
    across(a, n1, x1, c1);
    check(n0 > 0 && x0 > 0 && std::fabs(n1 - 0.5 * n0) < 1e-9 * n0 && std::fabs(x1 - 0.5 * x0) < 1e-9 * x0, "MULTREGT ALL scales neighbour and fault-NNC transmissibility across the region boundary");
    Case b = f; b.grid.multregt.push_back({2, 1, 0.25, "XYZ", 1, 'F'});
    across(b, n1, x1, c1);
    check(std::fabs(n1 - n0) < 1e-9 * n0 && std::fabs(x1 - 0.25 * x0) < 1e-9 * x0, "MULTREGT NNC (reversed region order) scales only the non-neighbour connections");
    Case d = f; d.grid.multregt.push_back({1, 2, 0.0, "X", 2, 'F'});
    across(d, n1, x1, c1);
    check(n1 == 0 && std::fabs(x1 - x0) < 1e-9 * x0 && c1 < c0, "MULTREGT 0 NONNC removes the neighbour connections only");
    Case e = f; e.grid.multregt.push_back({1, 2, 0.5, "Y", 0, 'F'});
    across(e, n1, x1, c1);
    check(std::fabs(n1 - n0) < 1e-9 * n0 && std::fabs(x1 - x0) < 1e-9 * x0, "MULTREGT with direction Y leaves X connections untouched");
    {
        const std::string dk = eclipseDeckToString(b);
        const Case r = eclipseDeckFromString(dk);
        double n2, x2; size_t c2; across(r, n2, x2, c2);
        check(dk.find("MULTREGT") != std::string::npos && dk.find("FLUXNUM") != std::string::npos && r.grid.multregt.size() == 1 &&
              r.grid.multregt[0].nnc == 1 && std::fabs(x2 - 0.25 * x0) < 1e-9 * x0 && eclipseDeckToString(r) == dk, "MULTREGT + FLUXNUM: deck round trip");
    }
    // producer connection shut at step 10 and re-opened at step 20 (COMPDAT SHUT ... COMPDAT OPEN)
    Case c = makeExample("spe1");
    int prod = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].type == WellType::Producer) prod = int(w);
    auto& comps = c.wells[size_t(prod)].completions;
    Completion re = comps[0]; re.fromStep = 20; re.untilStep = -1;
    comps[0].untilStep = 10; comps.push_back(re);
    const std::string d1 = eclipseDeckToString(c);
    const Case rb = eclipseDeckFromString(d1);
    const auto& rc = rb.wells[size_t(prod)].completions;
    bool ok = rc.size() == 2 && eclipseDeckToString(rb) == d1;
    if (ok) ok = rc[0].fromStep <= 0 && rc[0].untilStep == 10 && rc[1].fromStep == 20 && rc[1].untilStep < 0;
    check(ok, "COMPDAT shut -> re-open: one open interval per period, deck round trip is stable");
    Simulator s; s.log = [](const std::string&) {}; s.setup(rb); s.initialize();
    double q15 = -1, q25 = -1;
    for (int st = 0; st < 26 && s.runReportStep(); ++st) { if (st == 14) q15 = s.summary().back().fopr; if (st == 24) q25 = s.summary().back().fopr; }
    std::printf("  producer FOPR at step 15 %.6g, step 25 %.6g\n", q15, q25);
    check(q15 == 0.0 && q25 > 0.0, "re-opened connection produces only in its open periods");
    {   // WCONHIST item 10 is the observed BHP (history data), not a limit; a limit comes from WELTARG
        Case h = makeExample("spe1");
        for (auto& st : h.schedule) for (auto& wc : st.controls) if (wc.mode != Control::RATE && wc.mode != Control::BHP) { wc.history = true; wc.mode = Control::RESV; wc.resv = -1; wc.bhp = -1; wc.bhpHist = 150; }
        std::string d = eclipseDeckToString(h);
        const Case h1 = eclipseDeckFromString(d);
        bool okh = false;
        for (auto& st : h1.schedule) for (auto& wc : st.controls) if (wc.history) okh = std::fabs(wc.bhpHist - 150) < 1e-9 && wc.bhp < 0;
        Case h2 = h1; for (auto& st : h2.schedule) for (auto& wc : st.controls) if (wc.history) wc.bhp = 1000;
        const std::string d2 = eclipseDeckToString(h2);
        const Case h3 = eclipseDeckFromString(d2);
        bool okl = d2.find("WELTARG") != std::string::npos && eclipseDeckToString(h3) == d2;
        for (auto& st : h3.schedule) for (auto& wc : st.controls) if (wc.history) okl = okl && std::fabs(wc.bhp - 1000) < 1e-9;
        check(okh && okl, "WCONHIST observed BHP is history data (limit stays 1 atm); a BHP limit round-trips as WELTARG");
    }
}

static void testThreePointDegenerate() {
    std::printf("== SCALECRS three-point krw scaling with a table residual oil of zero (SR_t = SWU_t, Norne SWOF)\n");
    Case c = makeExample("qfs");
    c.props.swof = {{0.0, 0.0, 1.0, 0.0}, {0.0001, 0.0, 0.999, 0.0}, {0.5, 0.0625, 0.05, 0.0}, {1.0, 0.5, 0.0, 0.0}};
    c.grid.scalecrs = true; c.grid.endscale = true;
    SatFunctions sf; sf.setup(c, 0);
    SatFunctions::EndPoints e; e.swl = 0.05; e.swcr = 0.13; e.swu = 1.0; e.sowcr = 0.13; e.sgl = 0.0;
    const double k1 = sf.krw(AD3(0.5), &e).v, k2 = sf.krw(AD3(0.9), &e).v, k0 = sf.krw(AD3(0.12), &e).v;
    const double st = 0.0001 + (0.5 - 0.13) * (1.0 - 0.0001) / (0.87 - 0.13), kx = 0.0625 + (st - 0.5) / 0.5 * 0.4375;
    std::printf("  krw(0.5) %.6g (expected %.6g), krw(0.9) %.6g, krw(0.12) %.3g\n", k1, kx, k2, k0);
    check(std::fabs(k1 - kx) < 1e-9 && std::fabs(k2 - 0.5) < 1e-12 && k0 == 0.0,
          "three-point krw: [SWCR, 1-SOWCR] -> [SWCR_t, SWU_t], maximum above SR (OPM / ECLIPSE, no two-point fallback)");
}

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
    // ENDSCALE deck round trip; arrays without ENDSCALE are applied as OPM Flow does (C34)
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
    check(b2.grid.hasEndPointScaling(), "EPS: SWCR without ENDSCALE is applied (as OPM Flow, C34)");
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
    {   // DRSDT followed by VAPPARS in the same step (Norne): VAPPARS switches DRSDT off (OPM OilVaporizationProperties)
        Case c = makeExample("spe1");
        Simulator s0; s0.log = [](const std::string&) {}; s0.setup(c); s0.initialize();
        while (s0.runReportStep()) {}
        c.schedule[0].drsdt = 0; c.schedule[0].vap1 = 0; c.schedule[0].vap2 = 0; c.schedule[0].vapAfterDrsdt = true;
        const std::string deck = eclipseDeckToString(c);
        const Case r = eclipseDeckFromString(deck, "");
        const size_t pd = deck.find("DRSDT"), pv = deck.find("VAPPARS");
        check(!r.schedule.empty() && r.schedule[0].vapAfterDrsdt && pd != std::string::npos && pv != std::string::npos && pd < pv,
              "DRSDT then VAPPARS: order survives the deck round trip");
        Simulator s; s.log = [](const std::string&) {}; s.setup(r); s.initialize();
        while (s.runReportStep()) {}
        std::printf("  SPE1 + DRSDT 0 + VAPPARS 0 0: final FPR %.6g (no DRSDT %.6g)\n", s.summary().back().fpr, s0.summary().back().fpr);
        check(std::fabs(s.summary().back().fpr - s0.summary().back().fpr) < 1e-6 * s0.summary().back().fpr,
              "VAPPARS after DRSDT switches the DRSDT limit off (last keyword wins, as OPM / ECLIPSE)");
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
        // C35: the adsorbing rock volume is now (1 - phi) V as in OPM (it was V - PV with PV in rb, V in ft3); 1.2 % -> 1.5 %
        check(std::fabs(fopt[0] / 545320 - 1) < 0.02, "polymer example: cumulative oil within 2 % of OPM Flow");
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
    Case m = d;
    d.grid.actnum.assign(size_t(d.grid.size()), 1);
    for (int q = 6 * nl; q < 7 * nl; ++q) d.grid.actnum[size_t(q)] = 0;    // a thick layer switched off with ACTNUM
    Grid A; A.build(d);
    check(A.numPinch == 0, "a thick ACTNUM = 0 layer is a barrier even with GAP (as ECLIPSE on Norne)");
    m.grid.minpv = 1e-3;                                                   // a thick layer removed by MINPV
    for (int q = 6 * nl; q < 7 * nl; ++q) m.grid.poro[size_t(q)] = 1e-12;
    Grid D; D.build(m);
    m.grid.pinch.gap = false; Grid E; E.build(m);
    check(D.numPinch == nl && E.numPinch == 0, "GAP bridges thick MINPV-inactive layers, NOGAP does not");
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


static void testNumerics() {
    std::printf("== Numerics framework: BSR / BILU0 / well Schur / block CPR / AMG reuse / AD / local Newton\n");
    // ---- coupled block system: 1-D chain of na cells x 3 unknowns + 2 wells (scalar CSR as assembled by the engine)
    const int na = 300, nv = 3, nw = 2, nc = na * nv, n = nc + nw;
    std::vector<std::map<int, double>> rows(n);
    std::mt19937 rng(11); std::uniform_real_distribution<double> U(-1, 1);
    for (int a = 0; a < na; ++a)
        for (int e = 0; e < nv; ++e) {
            const int r = a * nv + e;
            for (int v = 0; v < nv; ++v) rows[r][a * nv + v] = (e == v ? 8.0 : 0.0) + U(rng) * (e == 0 || v == 0 ? 2.0 : 0.3);
            if (e == 2) rows[r][a * nv + 2] = 1e-13;               // vanishing phase: scalar diagonal ~ 0
            for (int nb : {a - 1, a + 1}) if (nb >= 0 && nb < na) for (int v = 0; v < nv; ++v) rows[r][nb * nv + v] = -1.0 + 0.2 * U(rng);
        }
    const int perf[2][3] = {{10, 11, 12}, {200, 250, 290}};
    for (int w = 0; w < nw; ++w) {
        double d = 0;
        for (int k = 0; k < 3; ++k) {
            const int a = perf[w][k];
            for (int e = 0; e < nv; ++e) { rows[a * nv + e][nc + w] = -0.7; rows[a * nv + e][a * nv + e] += 0.7; }
            rows[nc + w][a * nv] = -0.7; d += 0.7;
        }
        rows[nc + w][nc + w] = d * 1.1;
    }
    CSR A; A.n = n; A.rowPtr.assign(n + 1, 0);
    for (int r = 0; r < n; ++r) { for (auto& [c, v] : rows[r]) { A.col.push_back(c); A.val.push_back(v); } A.rowPtr[r + 1] = int(A.col.size()); }
    A.finalizeDiag();
    std::vector<double> xt(n), b(n), r(n);
    for (int i = 0; i < n; ++i) xt[i] = std::sin(0.1 * i) + 2.0;
    A.spmv(xt.data(), b.data());
    // BSR view equals the CSR cell block
    {
        CSR Ac; Ac.n = nc; Ac.rowPtr.assign(nc + 1, 0);
        for (int i = 0; i < nc; ++i) { for (int p = A.rowPtr[i]; p < A.rowPtr[i + 1]; ++p) if (A.col[p] < nc) { Ac.col.push_back(A.col[p]); Ac.val.push_back(A.val[p]); } Ac.rowPtr[i + 1] = int(Ac.col.size()); }
        const BSR B = BSR::fromCSR(Ac, na, nv);
        std::vector<double> y1(nc), y2(nc); Ac.spmv(xt.data(), y1.data()); B.spmv(xt.data(), y2.data());
        double e = 0; for (int i = 0; i < nc; ++i) e = std::max(e, std::fabs(y1[i] - y2[i]));
        check(e < 1e-12, "BSR spmv equals scalar CSR spmv");
    }
    WellSchur W;
    check(W.build(A, na, nv, nw), "well Schur complement builds");
    const BSR& S = W.reduced();
    std::vector<double> bc(nc), xc, x(n);
    W.reduceRhs(b.data(), bc.data());
    auto fullErr = [&] { W.recover(b.data(), xc.data(), x.data()); double e = 0; for (int i = 0; i < n; ++i) e = std::max(e, std::fabs(x[i] - xt[i])); return e; };
    {
        BILU0 M; M.factor(S);
        const LinResult lr = gmres(BsrOperator(S), M, bc, xc, 1e-12, 200, 40);
        check(lr.converged && fullErr() < 1e-8, "Schur-reduced BILU0-GMRES recovers the exact full solution (cells + BHPs)");
        ILU0 I; I.factor(A); std::vector<double> xs;
        const LinResult ls = gmres(A, I, b, xs, 1e-12, 200, 40);
        std::printf("    block path %d its, scalar ILU0 path %d its (conv %d)\n", lr.iterations, ls.iterations, int(ls.converged));
        check(lr.iterations <= ls.iterations || !ls.converged, "block path needs no more iterations than scalar ILU0 with a vanishing diagonal");
    }
    {
        BlockCPR P; P.setup(S);
        const LinResult lr = bicgstab(BsrOperator(S), P, bc, xc, 1e-11, 200);
        check(lr.converged && fullErr() < 1e-7, "block CPR (AMG + BILU0) BiCGSTAB solves the reduced system");
        BSR S2 = S; for (double& v : S2.val) v *= 1.5;              // next Newton iteration: same pattern, new values
        BlockCPR P2; P2.setup(S2);
        P.setup(S2, nullptr, true);
        std::vector<double> z1(nc), z2(nc); P.apply(bc.data(), z1.data()); P2.apply(bc.data(), z2.data());
        double e = 0, m = 0; for (int i = 0; i < nc; ++i) { e = std::max(e, std::fabs(z1[i] - z2[i])); m = std::max(m, std::fabs(z2[i])); }
        check(P.amgReused() && e <= 1e-10 * m, "AMG hierarchy reuse (refresh) matches a fresh setup on the same pattern");
    }
    {   // block-tridiagonal: BILU0 is an exact factorisation
        BSR T; T.n = 50; T.nb = 3; T.rowPtr.assign(51, 0);
        for (int i = 0; i < 50; ++i) { for (int j = std::max(0, i - 1); j <= std::min(49, i + 1); ++j) T.col.push_back(j); T.rowPtr[i + 1] = int(T.col.size()); }
        T.val.resize(T.col.size() * 9);
        for (int i = 0; i < 50; ++i) for (int p = T.rowPtr[i]; p < T.rowPtr[i + 1]; ++p) for (int k = 0; k < 9; ++k)
            T.blk(p)[k] = (T.col[p] == i ? (k % 4 == 0 ? 6.0 : 1.0) : -0.5) + 0.1 * U(rng);
        T.finalizeDiag();
        BILU0 M; M.factor(T);
        std::vector<double> xx(150), bb(150), zz(150); for (int i = 0; i < 150; ++i) xx[i] = 1.0 + i % 7;
        T.spmv(xx.data(), bb.data()); M.apply(bb.data(), zz.data());
        double e = 0; for (int i = 0; i < 150; ++i) e = std::max(e, std::fabs(zz[i] - xx[i]));
        check(e < 1e-10, "BILU0 is exact on a block-tridiagonal matrix");
    }
    // ---- AD drivers
    {
        auto f = [](const auto& v) { using T = std::decay_t<decltype(v[0])>; return std::array<T, 2>{v[0] * v[0] * v[1], exp(v[0]) + log(v[1])}; };
        std::array<double, 2> fx; std::array<std::array<double, 2>, 2> J;
        jacobian<2, 2>(f, {1.5, 2.0}, fx, J);
        check(std::fabs(J[0][0] - 6.0) < 1e-12 && std::fabs(J[0][1] - 2.25) < 1e-12 && std::fabs(J[1][0] - std::exp(1.5)) < 1e-12 &&
                  std::fabs(J[1][1] - 0.5) < 1e-12, "ForwardDiff-style jacobian() is exact");
        const double d = derivative([](auto x) { return pow(x, 2.5) * smoothMax(x, decltype(x)(1.0), 1e-9); }, 2.0);
        check(std::fabs(d - 3.5 * std::pow(2.0, 2.5)) < 1e-6, "derivative() through pow / smoothMax");
    }
    {   // local Newton: circle-hyperbola intersection, and a Rachford-Rice flash with bounds
        auto g = [](const auto& v) { using T = std::decay_t<decltype(v[0])>; return std::array<T, 2>{v[0] * v[0] + v[1] * v[1] - 4.0, v[0] * v[1] - 1.0}; };
        std::array<double, 2> x0{3.0, 0.1};
        const LocalNewtonResult lr = localNewton<2>(g, x0);
        check(lr.converged && std::fabs(x0[0] * x0[1] - 1) < 1e-9 && std::fabs(x0[0] * x0[0] + x0[1] * x0[1] - 4) < 1e-9, "localNewton converges (line search)");
        const double z[3] = {0.5, 0.3, 0.2}, K[3] = {3.0, 0.8, 0.1};
        auto rr = [&](const auto& v) {
            using T = std::decay_t<decltype(v[0])>; T s(0.0);
            for (int i = 0; i < 3; ++i) s += z[i] * (K[i] - 1.0) / (1.0 + v[0] * (K[i] - 1.0));
            return std::array<T, 1>{s};
        };
        std::array<double, 1> V{0.9}, lo{1e-9}, hi{1 - 1e-9};
        const LocalNewtonResult lf = localNewton<1>(rr, V, {}, &lo, &hi);
        double s = 0; for (int i = 0; i < 3; ++i) s += z[i] * (K[i] - 1.0) / (1.0 + V[0] * (K[i] - 1.0));
        check(lf.converged && V[0] > 0 && V[0] < 1 && std::fabs(s) < 1e-10, "Rachford-Rice flash solved with bounded local Newton");
    }
    {
        OscillationDetector od(2);
        const double seq[5][2] = {{1, 1}, {0.5, 0.6}, {1.02, 0.98}, {0.51, 0.6}, {1.0, 1.0}};
        bool seen = false; for (auto& q : seq) { od.push(q); seen = seen || od.oscillating(); }
        check(seen && od.omega() < 1.0, "oscillating Newton residuals detected, update relaxed");
        PidStepController pid; pid.reset();
        check(pid.factor(0.5) > 1.0 && pid.factor(4.0) < 1.0, "PID step controller grows small-change steps, shrinks large-change ones");
    }
    // ---- whole simulation: block (Schur + BSR) path vs scalar path, ILU and CPR
    {
        Case c = makeDemoCase(); c.schedule.resize(6);
        double fpr[4];
        for (int k = 0; k < 4; ++k) {
            c.tuning.linearSolver = k & 1; c.tuning.blockSolver = (k >> 1) != 0;
            Simulator sim; sim.setup(c); sim.initialize();
            while (sim.runReportStep()) {}
            fpr[k] = sim.summary().back().fpr;
        }
        std::printf("    FPR scalar ILU %.5f  scalar CPR %.5f  block BILU0 %.5f  block CPR %.5f\n", fpr[0], fpr[1], fpr[2], fpr[3]);
        check(std::fabs(fpr[2] - fpr[0]) < 0.01 && std::fabs(fpr[3] - fpr[1]) < 0.01, "block solver path gives the same solution as the scalar path");
    }
}

static void testNornePhysics() {
    std::printf("Hysteresis (EHYSTR), SCALECRS, VAPPARS, VFP / THP and group control\n");
    {   // EQUALS records that start with an array name in column 1 (Norne): every record applies
        Case c = makeExample("spe1");
        std::string d = eclipseDeckToString(c);
        const size_t r = d.find("\nRUNSPEC"), p = d.find("\nPROPS");
        check(r != std::string::npos && p != std::string::npos, "spe1 deck has RUNSPEC / PROPS");
        d.insert(d.find('\n', p + 1) + 1, "EQUALS\nSWL 0.11 1 10 1 10 1 1 /\nSWL 0.22 1 10 1 10 2 2 / comment\nSWL 0.33 1 10 1 10 3 3 /\n/\nCOPY\n  SWL SWCR /\n/\n");
        d.insert(d.find('\n', r + 1) + 1, "ENDSCALE\n/\n");
        const Case e = eclipseDeckFromString(d, "");
        const size_t L = 100;
        check(e.grid.swl.size() == 3 * L && std::fabs(e.grid.swl[0] - 0.11) < 1e-12 && std::fabs(e.grid.swl[L] - 0.22) < 1e-12 && std::fabs(e.grid.swl[2 * L] - 0.33) < 1e-12
              && e.grid.swcr.size() == 3 * L && std::fabs(e.grid.swcr[2 * L] - 0.33) < 1e-12,
              "EQUALS: records 'SWL v box /' on consecutive lines are records, not keywords");
    }
    {   // VFP table: multilinear interpolation, exact at the nodes
        VfpTable t; t.inj = false; t.floV = {100, 1000}; t.thpV = {100, 200}; t.wfrV = {0, 1}; t.gfrV = {1, 2}; t.alqV = {0};
        t.bhp.resize(16);
        for (int it = 0; it < 2; ++it) for (int iw = 0; iw < 2; ++iw) for (int ig = 0; ig < 2; ++ig) for (int f = 0; f < 2; ++f)
            t.bhp[size_t(((it * 2 + iw) * 2 + ig) * 2 + f)] = t.thpV[size_t(it)] + 1000 * t.wfrV[size_t(iw)] + 10 * t.gfrV[size_t(ig)] + 2 * t.floV[size_t(f)];
        const double exact = 150 + 1000 * 0.25 + 10 * 1.5 + 2 * 550;   // linear function -> exact interpolation
        check(std::fabs(t.eval(550, 150, 0.25, 1.5, 0) - exact) < 1e-9 && std::fabs(t.eval(1000, 200, 1, 2, 0) - (200 + 1000 + 20 + 2000)) < 1e-9,
              "VFP: multilinear interpolation exact for a linear table");
        check(std::fabs(t.eval(2000, 100, 0, 1, 0) - (100 + 10 + 4000)) < 1e-9, "VFP: linear extrapolation in the flow axis");
    }
    {   // SCALECRS three-point scaling: the displacing-phase residual point maps onto the table's
        Case c = makeExample("spe1");
        c.grid.endscale = true; c.grid.scalecrs = true;
        SatFunctions sf; sf.setup(c, 0);
        SatFunctions::EndPoints e = sf.tableEndPoints();
        const SatFunctions::EndPoints t = e;
        e.sowcr = std::min(0.9, t.sowcr + 0.1);
        const double krAtR = sf.krw(AD3(1.0 - e.sowcr - e.sgl), &e).v, krTab = sf.krw(AD3(1.0 - t.sowcr - t.sgl), nullptr).v;
        check(sf.threePoint() && std::fabs(krAtR - krTab) < 1e-9, "SCALECRS: krw at 1 - SOWCR - SGL equals the table value at its residual point");
        c.grid.scalecrs = false; SatFunctions s2; s2.setup(c, 0);
        check(!s2.threePoint(), "two-point scaling without SCALECRS");
    }
    {   // hysteresis: imbibition curve with a larger critical gas saturation traps gas
        Case c = makeExample("gascap");
        c.props.hyst.on = true; c.props.hyst.model = 0; c.props.hyst.what = "KR";
        c.grid.endscale = true;
        const size_t n = size_t(c.grid.nx * c.grid.ny * c.grid.nz);
        c.grid.sgcr.assign(n, -1.0);
        c.grid.isgcr.assign(n, 0.25);
        c.props.vap1 = 0; c.props.vap2 = 0;
        const Case r = eclipseDeckFromString(eclipseDeckToString(c), "");
        check(r.props.hyst.on && r.grid.isgcr.size() == n && std::fabs(r.grid.isgcr[0] - 0.25) < 1e-12 && r.grid.endscale,
              "SATOPTS HYSTER / EHYSTR / ISGCR survive the deck round trip");
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        Case d = makeExample("gascap");
        Simulator s0; s0.log = [](const std::string&) {}; s0.setup(d); s0.initialize();
        while (s0.runReportStep()) {}
        const SummaryPoint& a = s.summary().back(); const SummaryPoint& b = s0.summary().back();
        std::printf("  GASCAP hysteresis: FGPT %.6g (no hysteresis %.6g), FOPT %.6g (%.6g)\n", a.fgpt, b.fgpt, a.fopt, b.fopt);
        check(a.time >= b.time - 1e-6 && a.fgpt <= b.fgpt * 1.0001, "hysteresis run completes; trapped gas does not increase gas production");
    }
    {   // VAPPARS round trip and effect on a condensate (wet gas) case
        Case c = makeExample("condensate");
        c.props.vap1 = 0.5; c.props.vap2 = 0;
        const Case r = eclipseDeckFromString(eclipseDeckToString(c), "");
        check(std::fabs(r.props.vap1 - 0.5) < 1e-12, "VAPPARS survives the deck round trip");
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        check(s.summary().back().time > 1000, "condensate case with VAPPARS runs to the end");
    }
    {   // group production limit: FIELD ORAT below the wells' total
        Case c = makeExample("spe9like");
        Simulator s0; s0.log = [](const std::string&) {}; s0.setup(c); s0.initialize();
        s0.runReportStep();
        const double free = s0.summary().back().fopr;
        GroupProdControl g; g.group = "FIELD"; g.mode = Control::ORAT; g.orat = 0.6 * free;
        c.schedule[0].gprod.push_back(g);
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        s.runReportStep(); s.runReportStep();
        const double q = s.summary().back().fopr;
        std::printf("  SPE9LIKE group ORAT %.6g: FOPR %.6g (uncontrolled %.6g)\n", g.orat, q, free);
        check(std::fabs(q / g.orat - 1) < 0.02, "GCONPROD FIELD ORAT target honoured within 2 %");
    }
    {   // group injection limit (GCONINJE FIELD WATER RATE) with a stopped injector (own RATE target 0, Norne 2006):
        // the zero-rate well must not take a share of the group target
        Case c = makeExample("qfs");
        int inj = -1;
        for (size_t i = 0; i < c.wells.size(); ++i) if (c.wells[i].type == WellType::WaterInjector) { inj = int(i); break; }
        double q0 = 0;
        for (auto& wc : c.schedule[0].controls) if (inj >= 0 && wc.well == c.wells[size_t(inj)].name) q0 = wc.rate;
        bool okSetup = inj >= 0 && q0 > 0;
        if (okSetup) {
            WellInput z = c.wells[size_t(inj)]; z.name = "IZERO";
            for (auto& cp : z.completions) { cp.i = std::max(0, c.grid.nx - 1 - cp.i); }
            c.wells.push_back(z);
            WellControl wz; wz.well = "IZERO"; wz.mode = Control::RATE; wz.rate = 0; wz.bhp = 1e9; wz.phase = 1;
            c.schedule[0].controls.push_back(wz);
            GroupInjControl g; g.group = "FIELD"; g.phase = 1; g.mode = "RATE"; g.rate = 3 * q0;
            c.schedule[0].ginj.push_back(g);
        }
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        s.runReportStep(); s.runReportStep();
        const SummaryPoint& p = s.summary().back();
        const double qz = okSetup && p.wells.size() == c.wells.size() ? p.wells.back().waterInj : -1;
        const double qi = okSetup && p.wells.size() == c.wells.size() ? p.wells[size_t(inj)].waterInj : -1;
        std::printf("  GCONINJE FIELD RATE %.6g: injector %.6g (target %.6g), zero-target injector %.6g\n", 3 * q0, qi, q0, qz);
        check(okSetup && qz == 0.0 && std::fabs(qi / q0 - 1) < 1e-6, "group injection limit: an injector with RATE target 0 stays stopped");
    }
    {   // a producer driven to its 1 atm BHP limit, stopped by a zero target (WCONHIST ORAT 0 0 0), then reopened with a
        // small target (Norne B-1H 2004): the rate target must be honoured, not the stale 1 atm BHP
        Case c = makeExample("spe1");
        WellControl base; int prod = -1;
        for (auto& wc : c.schedule[0].controls) for (size_t i = 0; i < c.wells.size(); ++i)
            if (c.wells[i].name == wc.well && c.wells[i].type == WellType::Producer) { base = wc; prod = int(i); }
        const double atm = c.units == UnitSystem::Field ? 14.7 : 1.01325;
        bool okSetup = prod >= 0 && c.schedule.size() > 3;
        if (okSetup) {
            WellControl a = base; a.mode = Control::ORAT; a.orat = 1e7; a.bhp = atm; a.history = true;
            WellControl z = a; z.orat = 0; z.wrat = 0; z.grat = 0;
            WellControl b = a; b.orat = 500;
            for (auto& wc : c.schedule[0].controls) if (wc.well == base.well) wc = a;
            c.schedule[1].controls.push_back(z);
            c.schedule[2].controls.push_back(b);
        }
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        for (int i = 0; i < 3 && s.runReportStep(); ++i) {}
        const auto& S = s.summary();
        const double q0 = okSetup ? S[1].wells[size_t(prod)].oil : -1, b0 = okSetup ? S[1].wells[size_t(prod)].bhp : -1;
        const double q = okSetup ? S.back().wells[size_t(prod)].oil : -1, bhp = okSetup ? S.back().wells[size_t(prod)].bhp : -1;
        std::printf("  stopped / reopened producer: at the BHP limit q %.6g (BHP %.4g); after ORAT 0 then ORAT 500: q %.6g, BHP %.6g\n", q0, b0, q, bhp);
        check(okSetup && std::fabs(b0 - atm) < 1e-6 && std::fabs(q / 500 - 1) < 1e-4 && bhp > 10 * atm,
              "a producer reopened after a zero target honours its new rate target (no stale 1 atm BHP)");
    }
    {   // WCONHIST: observed rates are data, not limits - a history well at its BHP limit must not switch to an observed
        // GRAT / WRAT (Norne D-1H 2001)
        Case c = makeExample("spe1");
        int prod = -1; const double atm = c.units == UnitSystem::Field ? 14.7 : 1.01325;
        for (auto& wc : c.schedule[0].controls) for (size_t i = 0; i < c.wells.size(); ++i)
            if (c.wells[i].name == wc.well && c.wells[i].type == WellType::Producer) {
                prod = int(i); wc.mode = Control::ORAT; wc.orat = 1e7; wc.grat = 1; wc.wrat = 1; wc.lrat = 1e7; wc.bhp = atm; wc.history = true;
            }
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        for (int i = 0; i < 2 && s.runReportStep(); ++i) {}
        const WellRates& r = s.summary().back().wells[size_t(prod >= 0 ? prod : 0)];
        std::printf("  WCONHIST ORAT 1e7 (observed GRAT 1): q %.6g, gas %.6g, BHP %.4g, mode %s\n", r.oil, r.gas, r.bhp, toString(r.mode));
        check(prod >= 0 && r.mode == Control::BHP && std::fabs(r.bhp - atm) < 1e-6 && r.oil > 1000,
              "WCONHIST: observed rates are not limits (the well stays on its BHP limit)");
    }
    {   // THP control: BHP = THP + c q (+ head) on a producer with a large rate target
        Case c = makeExample("spe1");
        VfpTable t; t.id = 1; t.inj = false; t.flo = "OIL"; t.floV = {0, 1e5}; t.thpV = {100, 2000}; t.wfrV = {0}; t.gfrV = {0}; t.alqV = {0};
        t.datum = 0;
        t.bhp = {100 + 0.0, 100 + 1e5 * 0.1, 2000 + 0.0, 2000 + 1e5 * 0.1};
        c.vfp.push_back(t);
        for (auto& st : c.schedule) for (auto& wc : st.controls)
            if (wc.mode != Control::RATE && wc.mode != Control::BHP) { wc.thp = 1500; wc.vfp = 1; wc.deckMode = "THP"; wc.orat = 1e6; wc.mode = Control::ORAT; }
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        for (int i = 0; i < 6 && s.runReportStep(); ++i) {}
        const SummaryPoint& p = s.summary().back();
        bool okAll = false;
        for (const WellRates& w : p.wells) if (w.oil > 0) {
            const double vfp = 1500 + 0.1 * w.oil;   // at the table datum (0); the well head term is added by the simulator
            std::printf("  THP producer: q %.6g, BHP %.6g, VFP(q) at datum %.6g\n", w.oil, w.bhp, vfp);
            okAll = w.bhp > vfp;    // reference depth below the datum: BHP = VFP + hydrostatic head
        }
        check(okAll, "THP control: the producer's BHP follows the lift curve");
    }
}



// ---------------------------------------------------------------------------------------------- dual porosity: GRAVDR, DPNUM
static void testDualGravityDpnum() {
    std::printf("== Dual porosity: gravity drainage (GRAVDR / DZMTRX) and single porosity blocks (DPNUM)\n");
    // one matrix block (Pc = 0) next to a water-filled fracture: without GRAVDR nothing drives water into the matrix; with it
    // water imbibes from below until the normalised water columns are equal (here Swm = Swf = total water / total PV)
    const std::string deck =
        "RUNSPEC\nDIMENS\n 1 1 2 /\nOIL\nWATER\nFIELD\nDUALPORO\nGRAVDR\nSTART\n 1 JAN 2025 /\n"
        "GRID\nDX\n 2*100 /\nDY\n 2*100 /\nDZ\n 2*50 /\nTOPS\n 2*5000 /\nPORO\n 0.2 0.05 /\nPERMX\n 1 100 /\nPERMY\n 1 100 /\nPERMZ\n 1 100 /\n"
        "SIGMA\n 0.04 /\nDZMTRX\n 50 /\nNODPPM\n"
        "PROPS\nPVTW\n 3000 1.0 3e-6 0.5 0 /\nROCK\n 3000 4e-6 /\nDENSITY\n 50 64 0.06 /\nPVDO\n 1000 1.0 2.0\n 5000 0.99 2.0 /\n"
        "SWOF\n 0.2 0 1 0\n 0.8 1 0 0 /\nSOLUTION\nPRESSURE\n 2*3000 /\nSWAT\n 0.2 0.8 /\nSCHEDULE\nTSTEP\n 40*50 /\nEND\n";
    auto runSw = [](const Case& c, double& swm, double& swf, double& mb) {
        Simulator s; s.setup(c); s.initialize();
        const InPlace a = s.inPlace();
        while (s.runReportStep()) {}
        const Snapshot& sn = s.snapshots().back();
        swm = sn.sw[0]; swf = sn.sw[1];
        mb = std::fabs(s.inPlace().water - a.water) / a.water;
    };
    Case g = eclipseDeckFromString(deck);
    const bool parsed = g.grid.dual == 1 && g.grid.gravdr && g.grid.dzmtrxv.size() == 2 && std::fabs(g.grid.dzmtrxv[0] - 50) < 1e-12;
    check(parsed, "GRAVDR / DZMTRX import");
    Case n = g; n.grid.gravdr = false;
    double m1, f1, b1, m0, f0, b0;
    runSw(g, m1, f1, b1); runSw(n, m0, f0, b0);
    const double eq = (0.2 * 0.2 + 0.05 * 0.8) / 0.25;
    std::printf("  GRAVDR: Sw matrix %.4f fracture %.4f (equilibrium %.4f), water balance %.1e;  without: Sw matrix %.4f\n", m1, f1, eq, b1, m0);
    check(std::fabs(m0 - 0.2) < 1e-3, "no gravity drainage term: Pc = 0 matrix block takes no water");
    check(std::fabs(m1 - eq) < 0.01 && std::fabs(f1 - eq) < 0.02 && b1 < 1e-6, "GRAVDR: water columns equalise (Swm = Swf = 0.32), water conserved");
    {   // streamline / flux output uses the same GRAVDR potential: the M-F water flux equals the matrix accumulation rate
        // (backward Euler: the flux at the end of a step is exactly the step's mean rate)
        Simulator s; s.setup(g); s.initialize();
        std::vector<double> qw, qo, qg;
        for (int k = 0; k < 3; ++k) s.runReportStep();
        s.connectionFluxes(qw, qo, qg);
        const auto& sn = s.snapshots();
        const double pvm = 100.0 * 100.0 * 50.0 * 0.2 / 5.614583, dsw = sn.back().sw[0] - sn[sn.size() - 2].sw[0];
        const double acc = pvm * dsw / 50.0, q = qw.empty() ? 0 : std::fabs(qw[0]);
        std::printf("  GRAVDR connection flux: water M-F %.4f rb/d vs matrix accumulation %.4f rb/d\n", q, acc);
        check(qw.size() == 1 && acc > 0 && std::fabs(q - acc) < 0.01 * acc, "connection fluxes (streamlines) include the GRAVDR term (flux = matrix accumulation)");
    }
    {   // half the block height: same end state, slower approach
        Case h = g; h.grid.dzmtrxv.assign(2, 25.0); h.schedule.resize(2);
        Case f = g; f.schedule.resize(2);
        double mh, fh, bh, mf, ff, bf;
        runSw(h, mh, fh, bh); runSw(f, mf, ff, bf);
        check(mh - 0.2 < mf - 0.2 && mh > 0.2, "GRAVDR rate scales with the matrix block height Z");
    }
    {   // export -> import keeps GRAVDR / DZMTRX
        const Case b = eclipseDeckFromString(eclipseDeckToString(g));
        check(b.grid.gravdr && b.grid.dzmtrxv.size() == 2 && std::fabs(b.grid.dzmtrxv[1] - 50) < 1e-12, "GRAVDR / DZMTRX deck round trip");
    }
    {   // gas-oil: gas in the fracture drains oil from the top of the block until the normalised gas columns are equal
        const std::string gd =
            "RUNSPEC\nDIMENS\n 1 1 2 /\nOIL\nWATER\nGAS\nFIELD\nDUALPORO\nGRAVDR\nSTART\n 1 JAN 2025 /\n"
            "GRID\nDX\n 2*100 /\nDY\n 2*100 /\nDZ\n 2*50 /\nTOPS\n 2*5000 /\nPORO\n 0.2 0.05 /\nPERMX\n 1 100 /\nPERMY\n 1 100 /\nPERMZ\n 1 100 /\n"
            "SIGMA\n 0.04 /\nDZMTRX\n 50 /\nNODPPM\n"
            "PROPS\nPVTW\n 3000 1.0 3e-6 0.5 0 /\nROCK\n 3000 4e-6 /\nDENSITY\n 50 64 0.06 /\nPVDO\n 1000 1.0 2.0\n 5000 0.99 2.0 /\n"
            "PVDG\n 1000 2.0 0.02\n 5000 0.6 0.03 /\nSWOF\n 0.2 0 1 0\n 0.8 1 0 0 /\nSGOF\n 0 0 1 0\n 0.8 1 0 0 /\n"
            "SOLUTION\nPRESSURE\n 2*3000 /\nSWAT\n 0.2 0.2 /\nSGAS\n 0 0.75 /\nSCHEDULE\nTSTEP\n 40*50 /\nEND\n";
        Case c = eclipseDeckFromString(gd);
        Simulator s; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        const Snapshot& sn = s.snapshots().back();
        Case c0 = c; c0.grid.gravdr = false;
        Simulator s0; s0.setup(c0); s0.initialize(); while (s0.runReportStep()) {}
        std::printf("  GRAVDR gas: Sg matrix %.4f fracture %.4f (without: %.4f)\n", sn.sg[0], sn.sg[1], s0.snapshots().back().sg[0]);
        check(sn.sg[0] > 0.1 && std::fabs(sn.sg[0] - sn.sg[1]) < 0.02 && s0.snapshots().back().sg[0] < 1e-3,
              "GRAVDR gas-oil: gas enters the top of the block, gas columns equalise");
    }
    // DPNUM: the dual porosity example with a single porosity strip (i = 1..5, every layer)
    Case d = makeExample("dualporo");
    const int nl = d.grid.nx * d.grid.ny, half = d.grid.size() / 2;
    d.grid.dpnum.assign(size_t(half), 1);
    int nsp = 0;
    for (int q = 0; q < half; ++q) if (q % d.grid.nx < 5) { d.grid.dpnum[size_t(q)] = 0; ++nsp; }
    check(validateCase(d).empty(), "DPNUM case validates");
    Grid G0, G1; G0.build(makeExample("dualporo")); G1.build(d);
    int mf1 = 0; for (auto& cn : G1.conns) if (cn.dir == 'M') ++mf1;
    double pv0 = 0, pv1 = 0; for (double v : G0.poreVolRef) pv0 += v; for (double v : G1.poreVolRef) pv1 += v;
    // the moved cell carries the matrix pore volume; the fracture pore volume of single porosity blocks disappears
    double pvFracSp = 0, pvMatSp = 0;
    {
        const Case& c0 = d;
        for (int q = 0; q < half; ++q) if (d.grid.dpnum[size_t(q)] == 0) {
            const int a = G0.activeOf[size_t(q)], b = G0.activeOf[size_t(q + half)];
            if (a >= 0) pvMatSp += G0.poreVolRef[size_t(a)];
            if (b >= 0) pvFracSp += G0.poreVolRef[size_t(b)];
        }
        (void)c0;
    }
    std::printf("  DPNUM: %d single porosity blocks, active %d -> %d, M connections %d -> %d, PV %.6g -> %.6g\n", nsp, G0.numActive(),
                G1.numActive(), G0.numMF, mf1, pv0, pv1);
    check(G1.numActive() == G0.numActive() - nsp && mf1 == G0.numMF - nsp && std::fabs(pv1 - (pv0 - pvFracSp)) < 1e-9 * pv0,
          "DPNUM = 0: one active cell per block with the matrix pore volume, no matrix-fracture coupling");
    {   // a single porosity cell connects laterally to the fracture cell of its dual neighbour, with the harmonic matrix/fracture T
        const int qsp = d.grid.index(4, 3, 0) + half, qdu = d.grid.index(5, 3, 0) + half;
        bool found = false;
        for (auto& cn : G1.conns) {
            const int ga = G1.globalOf[size_t(cn.a)], gb = G1.globalOf[size_t(cn.b)];
            if ((ga == qsp && gb == qdu) || (ga == qdu && gb == qsp)) found = true;
        }
        check(found && G1.activeOf[size_t(d.grid.index(4, 3, 0))] < 0, "single porosity block connects to the neighbouring fracture system");
        (void)nl;
    }
    {
        const Case b = eclipseDeckFromString(eclipseDeckToString(d));
        check(b.grid.dpnum == d.grid.dpnum, "DPNUM deck round trip");
        Simulator s; s.setup(d); s.initialize();
        const InPlace a = s.inPlace();
        int st = 0; while (st < 8 && s.runReportStep()) ++st;
        const SummaryPoint& e = s.summary().back();
        const double mbO = (a.oil - s.inPlace().oil - e.fopt) / a.oil;
        std::printf("  DPNUM run: FOPT %.6g  oil balance %.1e\n", e.fopt, mbO);
        check(st == 8 && std::fabs(mbO) < 1e-5, "DPNUM model runs, oil balance < 1e-5");
    }
}

// ---------------------------------------------------------------------------------------------- local grid refinement
static void testLgr() {
    std::printf("== Local grid refinement (CARFIN / COMPDATL)\n");
    auto tsum = [](const Case& c, double& pv, size_t& nc) {
        Grid g; g.build(c); double t = 0; for (auto& x : g.conns) t += x.trans;
        pv = 0; for (double v : g.poreVolRef) pv += v; nc = g.conns.size(); return t;
    };
    auto run = [](const Case& c, int steps, SummaryPoint& out, int* newton = nullptr) {
        Simulator s; s.setup(c); s.initialize();
        for (int i = 0; (steps < 0 || i < steps) && s.runReportStep(); ++i) {}
        out = s.summary().back(); if (newton) *newton = s.stats().newton;
    };
    auto rel = [](double a, double b) { return std::fabs(a - b) / std::max(1e-30, std::fabs(b)); };
    {   // 1. trivial 1x1x1 refinement reproduces the base model exactly
        Case a = makeDemoCase(), b = a;
        LgrInput L; L.name = "L1"; L.i1 = 7; L.i2 = 11; L.j1 = 7; L.j2 = 11; L.k1 = 0; L.k2 = 3; L.nx = 5; L.ny = 5; L.nz = 4;
        b.grid.lgrs.push_back(L);
        double pa, pb; size_t na, nb; const double ta = tsum(a, pa, na), tb = tsum(b, pb, nb);
        SummaryPoint ra, rb; run(a, 8, ra); run(b, 8, rb);
        std::printf("  trivial: T %.10g / %.10g  PV %.10g / %.10g  FOPT %.10g / %.10g\n", ta, tb, pa, pb, ra.fopt, rb.fopt);
        check(rel(tb, ta) < 1e-10 && rel(pb, pa) < 1e-12 && na == nb && rel(rb.fopt, ra.fopt) < 1e-7,
              "trivial LGR (1 child per host): same transmissibilities, pore volume and production as the base grid");
    }
    {   // 2. refining every cell 2x2x2 equals the native fine grid
        Case f = makeDemoCase();
        const double kl[4] = {200, 50, 500, 100};
        for (int q = 0; q < f.grid.size(); ++q) { const int k = q / 400; f.grid.permx[size_t(q)] = f.grid.permy[size_t(q)] = kl[k]; f.grid.permz[size_t(q)] = 0.1 * kl[k]; f.grid.poro[size_t(q)] = 0.2; }
        Case co = f; co.grid.nx = 10; co.grid.ny = 10; co.grid.nz = 2; const int nc = 200;
        co.grid.dx.assign(nc, 600); co.grid.dy.assign(nc, 600); co.grid.dz.assign(nc, 50); co.grid.tops.assign(100, 8325);
        co.grid.poro.assign(nc, 0.2); co.grid.permx.assign(nc, 100); co.grid.permy = co.grid.permx; co.grid.permz.assign(nc, 10);
        LgrInput F; F.name = "ALL"; F.i1 = 0; F.i2 = 9; F.j1 = 0; F.j2 = 9; F.k1 = 0; F.k2 = 1; F.nx = 20; F.ny = 20; F.nz = 4;
        std::vector<double> lp(1600), lz(1600);
        for (int q = 0; q < 1600; ++q) { lp[size_t(q)] = kl[q / 400]; lz[size_t(q)] = 0.1 * kl[q / 400]; }
        F.props["PERMX"] = lp; F.props["PERMY"] = lp; F.props["PERMZ"] = lz;
        co.grid.lgrs.push_back(F);
        for (auto& w : co.wells) { w.lgr = 0; for (auto& cp : w.completions) cp.lgr = 0; }
        double p1, p2; size_t n1, n2; const double t1 = tsum(f, p1, n1), t2 = tsum(co, p2, n2);
        SummaryPoint r1, r2; run(f, 8, r1); run(co, 8, r2);
        std::printf("  full 2x2x2: T %.10g / %.10g  FOPT fine %.10g  LGR %.10g\n", t1, t2, r1.fopt, r2.fopt);
        check(validateCase(co).empty() && rel(t2, t1) < 1e-9 && rel(p2, p1) < 1e-9 && rel(r2.fopt, r1.fopt) < 1e-7,
              "full 2x2x2 refinement with local PERMX/Y/Z and COMPDATL wells = the native fine grid");
    }
    {   // 3. corner-point: trivial refinement across a fault throw keeps every fault NNC and MULTFLT
        Case a = makeExample("faulted"), b = a;
        LgrInput L; L.name = "FLT"; L.i1 = 17; L.i2 = 22; L.j1 = 5; L.j2 = 14; L.k1 = 0; L.k2 = 4; L.nx = 6; L.ny = 10; L.nz = 5;
        b.grid.lgrs.push_back(L);
        double pa, pb; size_t na, nb; const double ta = tsum(a, pa, na), tb = tsum(b, pb, nb);
        Case r = a; LgrInput R = L; R.nx = 18; R.ny = 30; R.nz = 10; r.grid.lgrs.push_back(R);
        double pr; size_t nr; tsum(r, pr, nr);
        std::printf("  faulted CP: T %.10g / %.10g  conns %zu / %zu  refined 3x3x2: PV %.10g / %.10g\n", ta, tb, na, nb, pr, pa);
        check(rel(tb, ta) < 1e-10 && na == nb, "trivial LGR across a fault throw: identical transmissibilities and NNCs");
        check(rel(pr, pa) < 1e-9, "3x3x2 refinement across a fault conserves pore volume");
        SummaryPoint ra, rb; run(a, 10, ra); run(b, 10, rb);
        check(rel(rb.fopt, ra.fopt) < 1e-6 && rel(rb.fwpt + 1, ra.fwpt + 1) < 1e-4, "faulted model with a trivial LGR: same production");
    }
    {   // 4. tools: LGR around a well, WELSPECS <-> WELSPECL round trip
        Case c = makeExample("qfs");
        const LgrInput L = makeLgrAroundWell(c, 1, 2, 3, 2, true);
        c.grid.lgrs.push_back(L);
        const auto before = c.wells[1].completions;
        moveWellToLgr(c, 1, 0);
        const bool inside = c.wells[1].lgr == 0 && c.wells[1].completions.size() == 1 && c.wells[1].completions[0].lgr == 0 &&
                            c.wells[1].i == 7 && c.wells[1].completions[0].k1 == 0 && c.wells[1].completions[0].k2 == 1;
        SummaryPoint r, base; run(c, -1, r); run(makeExample("qfs"), -1, base);
        Case back = c; moveWellToLgr(back, 1, -1);
        const auto& cb = back.wells[1].completions;
        const bool same = back.wells[1].lgr < 0 && cb.size() == before.size() && cb[0].i == before[0].i && cb[0].j == before[0].j &&
                          cb[0].k1 == before[0].k1 && cb[0].k2 == before[0].k2;
        std::printf("  qfs + LGR %s (%dx%dx%d): FOPT %.6g (coarse %.6g)\n", L.name.c_str(), L.nx, L.ny, L.nz, r.fopt, base.fopt);
        check(L.i1 == 22 && L.i2 == 24 && L.nx == 9 && L.nz == 2 && inside && same, "LGR around a well: box, WELSPECL/COMPDATL conversion and back");
        check(rel(r.fopt, base.fopt) < 0.03, "near-well LGR changes QFS cumulative oil by < 3 %");
        bool threw = false;
        try { makeLgrAroundWell(c, 0, 30, 2, 1, true); } catch (const std::exception&) { threw = true; }
        check(threw, "overlapping LGR boxes are rejected");
    }
    {   // 5. deck: CARFIN / NXFIN / HZFIN / EQUALS + MULTIPLY in CARFIN, REFINE in REGIONS, WELSPECL / COMPDATL; exact round trip
        Case q = makeExample("qfs");
        std::string d = eclipseDeckToString(q);
        auto ins = [&](const std::string& after, const std::string& text) { const size_t p = d.find(after); if (p != std::string::npos) d.insert(p, text); };
        ins("\nGRID", "\nLGR\n 2 2000 /\n");
        ins("\nPROPS", "\nCARFIN\n 'NEARP' 21 25 21 25 1 1 15 15 3 /\nNXFIN\n 5 3 3 2 2 /\nNYFIN\n 5 3 3 2 2 /\nHZFIN\n 1 2 1 /\n"
                       "EQUALS\n PORO 0.25 1 15 1 15 2 2 /\n/\nENDFIN\n\nCARFIN\n 'NEARI' 1 2 1 2 1 1 6 6 1 /\nMULTIPLY\n PERMX 2 /\n PERMY 2 /\n/\nENDFIN\n");
        ins("\nSOLUTION", "\nREGIONS\n\nFIPNUM\n 625*1 /\nREFINE\n 'NEARP' /\nFIPNUM\n 675*2 /\nENDFIN\n");
        // the producer moves into NEARP: WELSPECL + COMPDATL in local (15, 15, 1-3)
        const size_t ws = d.find("WELSPECS"); const size_t cd = d.find("COMPDAT");
        bool ok = ws != std::string::npos && cd != std::string::npos;
        if (ok) {
            const size_t pl = d.find("'PROD'", ws); const size_t pe = d.find('\n', pl);
            d.replace(pl, pe - pl, "/\nWELSPECL\n 'PROD' 'PROD' 'NEARP' 15 15 1* 'OIL' /");
            const size_t c0 = d.find("COMPDAT", d.find("WELSPECL")); const size_t cl = d.find("'PROD'", c0); const size_t ce = d.find('\n', cl);
            d.replace(cl, ce - cl, "/\nCOMPDATL\n 'PROD' 'NEARP' 15 15 1 3 'OPEN' 1* 1* 0.5 /");
        }
        Case a; std::string err;
        try { a = eclipseDeckFromString(d); err = validateCase(a); } catch (const std::exception& e) { err = e.what(); }
        if (!err.empty()) std::printf("  deck error: %s\n", err.c_str());
        ok = ok && err.empty() && a.grid.lgrs.size() == 2 && a.wells.size() == 2;
        if (ok) {
            const LgrInput& P = a.grid.lgrs[0];
            ok = P.nxfin == std::vector<int>{5, 3, 3, 2, 2} && P.hzfin.size() == 3 && P.props.count("PORO") && P.props.count("FIPNUM") &&
                 a.grid.lgrs[1].props.count("PERMX") && std::fabs(a.grid.lgrs[1].props.at("PERMX")[0] - 300) < 1e-9;
            int pw = a.wells[0].name == "PROD" ? 0 : 1;
            ok = ok && a.wells[size_t(pw)].lgr == 0 && a.wells[size_t(pw)].completions.size() == 1 && a.wells[size_t(pw)].completions[0].k2 == 2;
        }
        check(ok, "CARFIN/NXFIN/HZFIN/EQUALS/MULTIPLY, REFINE FIPNUM and WELSPECL/COMPDATL import");
        if (ok) {
            const std::string e1 = eclipseDeckToString(a);
            const Case b = eclipseDeckFromString(e1);
            const bool textSame = eclipseDeckToString(b) == e1;
            SummaryPoint ra, rb; int nwa, nwb; run(a, 20, ra, &nwa); run(b, 20, rb, &nwb);
            std::printf("  deck LGR: FOPT %.10g / %.10g  newton %d / %d  export has REFINE %d COMPDATL %d\n", ra.fopt, rb.fopt, nwa, nwb,
                        int(e1.find("REFINE") != std::string::npos), int(e1.find("COMPDATL") != std::string::npos));
            check(textSame && ra.fopt == rb.fopt && nwa == nwb && e1.find("REFINE") != std::string::npos && e1.find("WELSPECL") != std::string::npos,
                  "LGR deck export -> import is exact (text fixed point, identical run)");
        }
    }
    {   // 5b. WELOPENL: local-grid connection status (local indices of the named LGR)
        Case q = makeExample("qfs");
        std::string d = eclipseDeckToString(q);
        auto ins = [&](const std::string& after, const std::string& text) { const size_t p = d.find(after); if (p != std::string::npos) d.insert(p, text); };
        ins("\nGRID", "\nLGR\n 1 2000 /\n");
        ins("\nPROPS", "\nCARFIN\n 'NEARP' 21 25 21 25 1 1 15 15 3 /\nENDFIN\n");
        const size_t ws = d.find("WELSPECS");
        bool ok = ws != std::string::npos;
        if (ok) {
            const size_t pl = d.find("'PROD'", ws); const size_t pe = d.find('\n', pl);
            d.replace(pl, pe - pl, "/\nWELSPECL\n 'PROD' 'PROD' 'NEARP' 15 15 1* 'OIL' /");
            const size_t c0 = d.find("COMPDAT", d.find("WELSPECL")); const size_t cl = d.find("'PROD'", c0); const size_t ce = d.find('\n', cl);
            d.replace(cl, ce - cl, "/\nCOMPDATL\n 'PROD' 'NEARP' 15 15 1 1 'OPEN' 1* 1* 0.5 /\n 'PROD' 'NEARP' 15 15 2 2 'OPEN' 1* 1* 0.5 /\n 'PROD' 'NEARP' 15 15 3 3 'OPEN' 1* 1* 0.5 /");
            const size_t t0 = d.find("TSTEP", d.find("COMPDATL")); const size_t te = d.find("/\n", t0);
            ok = t0 != std::string::npos && te != std::string::npos;
            if (ok) d.insert(te + 2, "\nWELOPENL\n 'PROD' 'NEARP' 'SHUT' 15 15 2 /\n 'PROD' 'NEARP' 'SHUT' 3 3 1 /\n/\n");
        }
        Case a; std::string err;
        try { a = eclipseDeckFromString(d); err = validateCase(a); } catch (const std::exception& e) { err = e.what(); }
        if (!err.empty()) std::printf("  deck error: %s\n", err.c_str());
        int shut = 0, still = 0;
        if (ok && err.empty())
            for (const auto& w : a.wells) if (w.name == "PROD") for (const auto& cp : w.completions) { if (cp.k1 == 1 && cp.untilStep == 60) ++shut; if (cp.k1 != 1 && cp.untilStep < 0) ++still; }
        std::printf("  WELOPENL: shut %d (layer 2 after the 60 steps), untouched %d\n", shut, still);
        if (shut != 1) for (const auto& w : a.wells) for (const auto& cp : w.completions) std::printf("    %s lgr %d ijk %d %d %d-%d from %d until %d open %d\n", w.name.c_str(), cp.lgr, cp.i, cp.j, cp.k1, cp.k2, cp.fromStep, cp.untilStep, int(cp.open));
        if (shut != 1) { const size_t p = d.find("WELOPENL"); std::printf("%s\n", d.substr(p > 300 ? p - 300 : 0, 500).c_str()); }
        check(ok && err.empty() && shut == 1 && still == 2, "WELOPENL shuts the local connection (15,15,2) of LGR NEARP only");
    }
    {   // 6. LGR example: well-posed, runs, local cells present
        const Case c = makeExample("lgr");
        Grid g; g.build(c);
        SummaryPoint r; run(c, -1, r);
        std::printf("  example lgr: active %d (local %zu)  FOPT %.6g  FWCT %.3f\n", g.numActive(), size_t(g.numActive()) - 1764 + 81, r.fopt, r.fwct);
        check(validateCase(c).empty() && g.numActive() == 3861 && r.fopt > 0, "LGR example builds (3861 cells) and runs to the end");
    }
}

// ------------------------------------------------------------------------------------------- compositional thermodynamics (C1)
static void testEos() {
    using namespace ressim::eos;
    std::printf("EOS thermodynamics (PR/SRK, stability, flash, saturation, LBC)\n");
    // SPE5 fluid entered as E300 PROPS keywords (FIELD units, PRCORR = PR78)
    const std::string deckText = R"(RUNSPEC
FIELD
COMPS
 6 /
PROPS
EOS
 PR /
PRCORR
CNAMES
 C1 C3 C6 C10 C15 C20 /
TCRIT
 343.0 665.7 913.4 1111.8 1270.0 1380.0 /
PCRIT
 667.8 616.3 436.9 304.0 200.0 162.0 /
ZCRIT
 0.290 0.277 0.264 0.257 0.245 0.235 /
ACF
 0.0130 0.1524 0.3007 0.4885 0.6500 0.8500 /
MW
 16.04 44.10 86.18 142.29 206.0 282.0 /
BIC
 0.0
 0.0 0.0
 0.0 0.0 0.0
 0.05 0.005 0.0 0.0
 0.05 0.005 0.0 0.0 0.0 /
RTEMP
 160 /
ZI
 0.5 0.03 0.07 0.2 0.15 0.05 /
SOLUTION
ZMFVD
 8000 0.5 0.03 0.07 0.2 0.15 0.05
 9000 0.5 0.03 0.07 0.2 0.15 0.05 /
)";
    deck::RawDeck rd = deck::parseDeckText(deckText, "");
    CompFluid fl; std::vector<std::string> msgs;
    const bool ok = fluidFromDeck(rd, true, fl, &msgs);
    for (const auto& m : msgs) std::printf("    note: %s\n", m.c_str());
    check(ok && fl.eosRegions.size() == 1 && fl.eosRegions[0].nc() == 6 && fl.pr78 && fl.names[3] == "C10", "E300 PROPS keywords build a 6-component PR78 fluid");
    if (!ok) return;
    const Model& m = fl.eosRegions[0];
    check(std::fabs(m.comp[0].tc - 343.0 * 5 / 9) < 1e-9 && std::fabs(m.comp[0].pc - 667.8 * 6894.757293168) < 1e-6 &&
          std::fabs(m.k(4, 0) - 0.05) < 1e-15 && std::fabs(m.k(1, 5) - 0.005) < 1e-15 && std::fabs(fl.resT - fahrenheitToK(160)) < 1e-9 &&
          fl.zmfvd.size() == 1 && fl.zmfvd[0].depth.size() == 2 && std::fabs(fl.zmfvd[0].depth[1] - 9000 * 0.3048) < 1e-9 && fl.zi.size() == 6,
          "FIELD units, BIC lower triangle, RTEMP, ZMFVD and ZI converted to SI");
    const std::vector<double> z = fl.zi;
    const double T = fl.resT;
    const double pb = saturationPressure(m, T, z, true, 0);
    std::printf("  SPE5 oil bubble point at 160 F: %.2f psia (published 2302.9)\n", paToPsia(pb));
    check(std::fabs(paToPsia(pb) - 2302.9) < 1.0, "SPE5 bubble point matches the published 2302.9 psia (PR78)");
    const TParams tp = atTemperature(m, T);
    // phase stability and flash consistency around the bubble point
    const StabilityResult s1 = stability(m, tp, psiaToPa(4000), z), s2 = stability(m, tp, psiaToPa(2000), z);
    check(s1.stable && !s2.stable && s2.tpd < 0, "Michelsen stability: stable at 4000 psia, unstable at 2000 psia");
    const FlashResult fa = flash(m, tp, pb * 1.002, z), fb = flash(m, tp, pb * 0.998, z);
    std::printf("  flash at Pb+0.2%%: %d phase(s); Pb-0.2%%: %d phases, vapour fraction %.2e\n", fa.phases, fb.phases, fb.beta);
    check(fa.phases == 1 && !fa.vapourLabel && fb.phases == 2 && fb.beta > 0 && fb.beta < 2e-3, "flash agrees with the bubble point (single phase above, incipient vapour below)");
    const FlashResult f = flash(m, tp, psiaToPa(1500), z);
    double fug = 0, mb = 0, sx = 0, sy = 0;
    for (int i = 0; i < 6; ++i) {
        fug = std::max(fug, std::fabs(std::log(f.y[i]) + f.V.lnphi[i] - std::log(f.x[i]) - f.L.lnphi[i]));
        mb = std::max(mb, std::fabs(f.beta * f.y[i] + (1 - f.beta) * f.x[i] - z[i]));
        sx += f.x[i]; sy += f.y[i];
    }
    std::printf("  flash 1500 psia: V %.5f  rhoL %.2f rhoV %.2f kg/m3  |dlnf| %.1e  mass balance %.1e  iterations %d\n", f.beta, f.L.rho, f.V.rho, fug, mb, f.iterations);
    check(f.phases == 2 && f.converged && fug < 1e-10 && mb < 1e-12 && std::fabs(sx - 1) < 1e-12 && std::fabs(sy - 1) < 1e-12 && f.L.rho > f.V.rho,
          "two-phase flash: equal fugacities, component balance, normalised phases");
    // warm start with the converged K (simulator path) reproduces the solution in a couple of Newton iterations
    FlashOptions wo; wo.Kinit = &f.K;
    const FlashResult fw = flash(m, tp, psiaToPa(1500), z, wo);
    check(fw.phases == 2 && std::fabs(fw.beta - f.beta) < 1e-10 && fw.iterations <= 4, "K-value warm start converges immediately");
    // analytic derivatives vs central differences, with a Peneloux shift active
    Model ms = m;
    for (int i = 0; i < 6; ++i) ms.comp[i].sshift = -0.15 + 0.05 * i;
    const TParams tps = atTemperature(ms, T);
    const std::vector<double> n = {0.31, 0.07, 0.12, 0.25, 0.15, 0.10};
    double eN = 0, eP = 0, eT = 0, gd = 0, eH = 0, eV = 0;
    for (int root = 0; root < 2; ++root) {
        EvalOptions o; o.dT = true; o.hres = true; o.root = root ? Root::Vapour : Root::Liquid;
        const double P = psiaToPa(root ? 300 : 3000);
        const PhaseEval e = evalPhase(ms, tps, P, n, o);
        for (int j = 0; j < 6; ++j) {
            auto np = n, nm = n; np[j] += 1e-6; nm[j] -= 1e-6;
            const PhaseEval a = evalPhase(ms, tps, P, np, o), b = evalPhase(ms, tps, P, nm, o);
            double tot = 0; for (double q : n) tot += q;
            for (int i = 0; i < 6; ++i) eN = std::max(eN, std::fabs((a.lnphi[i] - b.lnphi[i]) / 2e-6 * tot - e.dlnphi_dn[i * 6 + j]));
            double s = 0; for (int i = 0; i < 6; ++i) s += n[i] * e.dlnphi_dn[i * 6 + j];
            gd = std::max(gd, std::fabs(s));
        }
        const PhaseEval pp = evalPhase(ms, tps, P * (1 + 1e-5), n, o), pm = evalPhase(ms, tps, P * (1 - 1e-5), n, o);
        const TParams t1 = atTemperature(ms, T + 1e-3), t2 = atTemperature(ms, T - 1e-3);
        const PhaseEval tpE = evalPhase(ms, t1, P, n, o), tmE = evalPhase(ms, t2, P, n, o);
        double gp = 0, gm = 0, tot = 0;
        for (double q : n) tot += q;
        for (int i = 0; i < 6; ++i) {
            eP = std::max(eP, std::fabs((pp.lnphi[i] - pm.lnphi[i]) / (2e-5 * P) - e.dlnphi_dP[i]) / std::fabs(e.dlnphi_dP[i]));
            eT = std::max(eT, std::fabs((tpE.lnphi[i] - tmE.lnphi[i]) / 2e-3 - e.dlnphi_dT[i]) / std::fabs(e.dlnphi_dT[i]));
            gp += n[i] / tot * tpE.lnphi[i]; gm += n[i] / tot * tmE.lnphi[i];
        }
        const double hfd = -Rgas * T * T * (gp - gm) / 2e-3;
        eH = std::max(eH, std::fabs(hfd - e.hres) / std::fabs(e.hres));
        eV = std::max(eV, std::fabs((pp.vm - pm.vm) / (2e-5 * P) - e.dvm_dP) / std::fabs(e.dvm_dP));
    }
    std::printf("  derivative errors: n %.1e  P %.1e  T %.1e  Gibbs-Duhem %.1e  Hres %.1e  dV/dP %.1e\n", eN, eP, eT, gd, eH, eV);
    check(eN < 1e-6 && eP < 1e-5 && eT < 1e-5 && gd < 1e-10 && eH < 1e-5 && eV < 1e-5, "analytic lnphi / volume / residual-enthalpy derivatives match finite differences (Gibbs-Duhem exact)");
    // flash sensitivities (implicit-function theorem) vs finite differences of complete flashes
    {
        const double P = psiaToPa(1500);
        FlashSens S; const bool okS = flashSensitivities(m, tp, P, z, f, S);
        const double h = P * 1e-6;
        const FlashResult a = flash(m, tp, P + h, z), b = flash(m, tp, P - h, z);
        double eP2 = std::fabs((a.beta - b.beta) / (2 * h) - S.dbeta_dP) / std::fabs(S.dbeta_dP);
        for (int i = 0; i < 6; ++i) eP2 = std::max(eP2, std::fabs((a.x[i] - b.x[i]) / (2 * h) - S.dx_dP[i]) / (std::fabs(S.dx_dP[i]) + 1e-30));
        const TParams t1 = atTemperature(m, T + 1e-3), t2 = atTemperature(m, T - 1e-3);
        const FlashResult c1 = flash(m, t1, P, z), c2 = flash(m, t2, P, z);
        const double eT2 = std::fabs((c1.beta - c2.beta) / 2e-3 - S.dbeta_dT) / std::fabs(S.dbeta_dT);
        double eZ = 0;
        for (int j = 0; j < 6; ++j) {
            auto zp = z, zm = z; zp[j] += 1e-6; zm[j] -= 1e-6;   // mole numbers (unnormalised feed)
            const FlashResult A = flash(m, tp, P, zp), B = flash(m, tp, P, zm);
            eZ = std::max(eZ, std::fabs((A.beta - B.beta) / 2e-6 - S.dbeta_dz[j]));
            for (int i = 0; i < 6; ++i) eZ = std::max(eZ, std::max(std::fabs((A.y[i] - B.y[i]) / 2e-6 - S.dy_dz[i * 6 + j]), std::fabs((A.x[i] - B.x[i]) / 2e-6 - S.dx_dz[i * 6 + j])));
        }
        std::printf("  flash sensitivities: dP rel %.1e  dT rel %.1e  dz abs %.1e\n", eP2, eT2, eZ);
        check(okS && eP2 < 1e-5 && eT2 < 1e-5 && eZ < 1e-6, "flash sensitivities d(beta, x, y)/d(P, T, z) match finite differences");
    }
    // dew point of a gas condensate and flash consistency
    const std::vector<double> zg = {0.90, 0.05, 0.02, 0.02, 0.007, 0.003};
    const double pd = saturationPressure(m, T, zg, false, 0);
    const FlashResult ga = flash(m, tp, pd * 1.003, zg), gb = flash(m, tp, pd * 0.997, zg);
    std::printf("  condensate dew point %.1f psia; above: %d phase (vapour %d), below: V = %.5f\n", paToPsia(pd), ga.phases, int(ga.vapourLabel), gb.beta);
    check(pd > 0 && ga.phases == 1 && ga.vapourLabel && gb.phases == 2 && gb.beta > 0.98, "dew point and flash agree (retrograde liquid appears below the dew point)");
    // pure-component vapour pressure (equal-fugacity roots) and dilute-gas LBC viscosity
    auto psat = [&](const char* name, double Tk) {
        Model q; q.comp.push_back(standardComponent(name));
        const TParams t = atTemperature(q, Tk);
        const std::vector<double> one = {1.0};
        double lo = 1e3, hi = q.comp[0].pc * 0.9999;
        EvalOptions a; a.dn = false; a.root = Root::Liquid; EvalOptions b = a; b.root = Root::Vapour;
        for (int k = 0; k < 200; ++k) {
            const double P = std::sqrt(lo * hi);
            const PhaseEval l = evalPhase(q, t, P, one, a), v = evalPhase(q, t, P, one, b);
            double d = l.lnphi[0] - v.lnphi[0];
            if (std::fabs(l.Z - v.Z) < 1e-9) d = l.Z > 0.3 ? 1 : -1;
            (d > 0 ? lo : hi) = P;
        }
        return lo;
    };
    const double pc3 = psat("C3", 300.0), pco2 = psat("CO2", 280.0);
    Model c1; c1.comp.push_back(standardComponent("C1"));
    EvalOptions gv; gv.root = Root::Vapour;
    const PhaseEval g1 = evalPhase(c1, atTemperature(c1, 300.0), 1e5, {1.0}, gv);
    const double mu1 = lbcViscosity(c1, {1.0}, 300.0, 1.0 / g1.vm) * 1e3;
    std::printf("  PR Psat C3(300 K) %.3f bar (NIST 9.97), CO2(280 K) %.2f bar (NIST 41.6); LBC CH4 1 bar 300 K %.4f cP (0.0112)\n", pc3 / 1e5, pco2 / 1e5, mu1);
    check(std::fabs(pc3 / 1e5 - 9.97) < 0.15 && std::fabs(pco2 / 1e5 - 41.6) < 0.8 && std::fabs(mu1 - 0.0112) < 0.0006, "pure-component vapour pressures and dilute-gas viscosity");
    // Rachford-Rice: analytic two-component case beta = (z1 (K1 - 1) + z2 (K2 - 1)) ... -> closed form for 2 components
    const std::vector<double> zz = {0.4, 0.6}, KK = {3.0, 0.2};
    const double br = rachfordRice(zz, KK, false);
    const double exact = -(zz[0] / (KK[1] - 1) + zz[1] / (KK[0] - 1));   // 2-component closed form
    check(std::fabs(br - exact) < 1e-13, "Rachford-Rice solution (two-component closed form)");
    // the deck importer reports the compositional fluid and its saturation pressure (deck itself runs in E300 / OPM)
    {
        std::vector<DeckMessage> rep; bool found = false;
        try { (void)eclipseDeckFromString(deckText, "", nullptr, &rep); } catch (...) {}
        for (const auto& r : rep) if (r.keyword == "EOS" && r.text.find("saturation pressure 2302.9 psia") != std::string::npos) { found = true; std::printf("  import report: %s\n", r.text.c_str()); }
        check(found, "deck import report shows the EOS fluid and its bubble point");
    }
    // SRK model: bubble point sanity and the same flash identities
    Model srk = m; srk.type = EosType::SRK;
    const double pbs = saturationPressure(srk, T, z, true, 0);
    const FlashResult fs = flash(srk, atTemperature(srk, T), psiaToPa(1500), z);
    double fugs = 0; for (int i = 0; i < 6; ++i) fugs = std::max(fugs, std::fabs(std::log(fs.y[i]) + fs.V.lnphi[i] - std::log(fs.x[i]) - fs.L.lnphi[i]));
    std::printf("  SRK bubble point %.1f psia, flash 1500 psia V %.4f\n", paToPsia(pbs), fs.beta);
    check(pbs > psiaToPa(1800) && pbs < psiaToPa(2900) && fs.phases == 2 && fugs < 1e-10, "SRK equation of state");
}

// ------------------------------------------------------------------------------------------- compositional simulator core (C2)
static ressim::eos::Model spe5Fluid() {
    using namespace ressim::eos;
    Model m; m.pr78 = true;
    auto add = [&](const char* n, double tcR, double pcpsi, double w, double mw, double zc) {
        Component c; c.name = n; c.tc = rankineToK(tcR); c.pc = psiaToPa(pcpsi); c.omega = w; c.mw = mw * 1e-3; c.zc = zc; m.comp.push_back(c);
    };
    add("C1", 343.0, 667.8, 0.013, 16.04, 0.290); add("C3", 665.7, 616.3, 0.1524, 44.10, 0.277); add("C6", 913.4, 436.9, 0.3007, 86.18, 0.264);
    add("C10", 1111.8, 304.0, 0.4885, 142.29, 0.257); add("C15", 1270.0, 200.0, 0.65, 206.0, 0.245); add("C20", 1380.0, 162.0, 0.85, 282.0, 0.235);
    m.setK(0, 4, 0.05); m.setK(0, 5, 0.05); m.setK(1, 4, 0.005); m.setK(1, 5, 0.005);
    return m;
}
static void testCompositional() {
    using namespace ressim::eos;
    using namespace ressim::comp;
    std::printf("Compositional FIM core (p, N_i + volume balance; analytic Jacobian; BlockCPR)\n");
    const std::vector<double> zoil = {0.5, 0.03, 0.07, 0.2, 0.15, 0.05}, zgas = {0.77, 0.20, 0.03, 0, 0, 0};
    const double T = fahrenheitToK(160);
    // 1. analytic Jacobian vs finite differences: vertical 3-cell column, single-phase oil above two two-phase cells,
    //    gravity, a gas injector and a BHP producer; with immobile and with mobile water (+ water injector)
    for (int wat = 0; wat < 3; ++wat) {   // 0 immobile water, 1 mobile water, 2 + rate-controlled two-connection producer
        CompCase c; c.nx = 1; c.ny = 1; c.nz = 3; c.dx = {30}; c.dy = {30}; c.dz = {5}; c.poro = {0.2};
        c.permx = c.permy = c.permz = {1e-13}; c.top = 2500; c.cr = 5e-10; c.pref = 2e7; c.swc = 0.2; c.sorg = 0.05; c.sgc = 0.03;
        c.fluid = spe5Fluid(); c.T = T;
        c.p0 = {psiaToPa(2700), psiaToPa(1900), psiaToPa(1700)}; c.z0 = {zoil};
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {1e-12}; inj.rate = 0.5; inj.zinj = zgas;
        CompWell prd; prd.cells = {2}; prd.wi = {2e-12}; prd.bhp = psiaToPa(1500);
        c.wells = {inj, prd};
        if (wat) {
            c.water = true; c.sorw = 0.15; c.sw0 = {0.3, 0.45, 0.6}; c.nw = 2.5; c.now = 1.8; c.krw0 = 0.5;
            CompWell wi; wi.injector = true; wi.water = true; wi.cells = {1}; wi.wi = {5e-12}; wi.rate = 2.0; c.wells.push_back(wi);
        }
        if (wat == 2) { c.wells[1].cells = {1, 2}; c.wells[1].wi = {1e-12, 2e-12}; c.wells[1].control = 1; c.wells[1].target = 30; c.wells[1].bhp = psiaToPa(500); }
        CompSimulator sim; sim.setup(c);
        const int ncw = sim.ncw(), nv = ncw + 1;
        std::vector<double> p = sim.pressure(), N = sim.moles();
        sim.storeOld();
        for (size_t k = 0; k < N.size(); ++k) N[k] *= 1.0 + 0.003 * double(k % 5);   // off-equilibrium iterate
        sim.setState(p, N);
        const double dt = 86400.0;
        std::vector<double> R0, J;
        sim.residualAndJacobian(dt, R0, &J);
        const int m = int(R0.size());
        double err = 0, scale = 0;
        for (int col = 0; col < m; ++col) {
            const int q = col / nv, k = col % nv;
            double Nt = 0; for (int i = 0; i < ncw; ++i) Nt += N[size_t(q) * ncw + i];
            const double h = k == 0 ? p[size_t(q)] * 1e-7 : Nt * 1e-7;
            std::vector<double> pp = p, NN = N, pm = p, Nm = N;
            if (k == 0) { pp[size_t(q)] += h; pm[size_t(q)] -= h; } else { NN[size_t(q) * ncw + k - 1] += h; Nm[size_t(q) * ncw + k - 1] -= h; }
            std::vector<double> Rp, Rm;
            sim.setState(pp, NN); sim.residualAndJacobian(dt, Rp, nullptr);
            sim.setState(pm, Nm); sim.residualAndJacobian(dt, Rm, nullptr);
            for (int r = 0; r < m; ++r) {
                const double fd = (Rp[size_t(r)] - Rm[size_t(r)]) / (2 * h), an = J[size_t(r) * m + col];
                err = std::max(err, std::fabs(fd - an) * (k == 0 ? p[size_t(q)] : Nt));
                scale = std::max(scale, std::fabs(an) * (k == 0 ? p[size_t(q)] : Nt));
            }
        }
        sim.setState(p, N);
        std::printf("  %s: phases per cell %d %d %d; Jacobian vs central differences: max scaled error %.2e\n",
                    wat == 0 ? "immobile water" : wat == 1 ? "mobile water" : "rate-controlled producer", sim.phasesIn(0), sim.phasesIn(1), sim.phasesIn(2), err / scale);
        if (wat == 2) std::printf("    producer BHP %.1f psia (rate control: %s)\n", sim.wellBhp(1) / 6894.757, sim.wellOnBhp(1) ? "no" : "yes");
        check(sim.phasesIn(0) == 1 && sim.phasesIn(1) == 2 && err / scale < 1e-5 && (wat < 2 || !sim.wellOnBhp(1)),
              wat == 0 ? "analytic compositional Jacobian (flash sensitivities, LBC, gravity, wells) matches finite differences"
              : wat == 1 ? "three-phase compositional Jacobian (mobile water, ECLIPSE default kro) matches finite differences"
                         : "rate-controlled producer: BHP eliminated implicitly (coupled connections), Jacobian matches finite differences");
    }
    // 1b. general geometry (cell list + connections), SWFN/SGFN/SOF3 tables with Pcow / Pcgo, ECLIPSE default and Stone II
    //     kro, implicit two-connection gas injector (rate / BHP limited), two-connection surface-rate producer (ORAT / LRAT,
    //     separator), lagged wellbore heads
    auto tabCase = [&](int variant) {
        CompCase c; c.fluid = spe5Fluid(); c.T = T; c.water = true; c.swc = 0.2; c.cr = 5e-10; c.pref = 2e7;
        for (int k = 0; k < 4; ++k) { c.cellPV.push_back(900.0 * (1 + 0.1 * k)); c.cellDepth.push_back(2502.5 + 5 * k); }
        c.conns = {{0, 1, 1.8e-11}, {1, 2, 1.5e-11}, {2, 3, 1.8e-11}, {0, 3, 2e-13}};   // last: a non-neighbour connection
        c.sat.sw = {0.2, 0.4, 0.6, 0.8, 1.0}; c.sat.krw = {0, 0.05, 0.2, 0.45, 1}; c.sat.pcow = {3e5, 1e5, 4e4, 1e4, 0};
        c.sat.sg = {0, 0.05, 0.3, 0.6, 0.8}; c.sat.krg = {0, 0, 0.15, 0.5, 0.9}; c.sat.pcog = {0, 1e4, 3e4, 6e4, 9e4};
        c.sat.so = {0, 0.15, 0.4, 0.6, 0.8}; c.sat.krow = {0, 0, 0.2, 0.5, 0.9}; c.sat.krog = {0, 0, 0.1, 0.45, 0.9};
        c.kroModel = variant == 1 ? 1 : 0;
        c.separator = {{300.0, 5e5}, {288.71, 101325.0}};
        c.p0 = {psiaToPa(2700), psiaToPa(1900), psiaToPa(1800), psiaToPa(1700)}; c.z0 = {zoil}; c.sw0 = {0.25, 0.35, 0.5, 0.6};
        CompWell inj; inj.injector = true; inj.cells = {0, 1}; inj.wi = {2e-12, 1e-12}; inj.zinj = zgas;
        inj.rate = variant == 1 ? 2000.0 : 0.5; inj.bhpMax = variant == 1 ? psiaToPa(2700) + 5e5 : 0.0;
        CompWell prd; prd.cells = {2, 3}; prd.wi = {1e-12, 2e-12}; prd.bhp = psiaToPa(500);
        prd.control = variant == 1 ? 6 : 3; prd.target = 1e-6;
        c.wells = {inj, prd};
        return c;
    };
    for (int variant = 0; variant < 2; ++variant) {
        CompSimulator sim; sim.setup(tabCase(variant));
        const int ncw = sim.ncw(), nv = ncw + 1;
        std::vector<double> p = sim.pressure(), N = sim.moles();
        sim.storeOld();
        for (size_t k = 0; k < N.size(); ++k) N[k] *= 1.0 + 0.003 * double(k % 5);
        sim.setState(p, N);
        const double dt = 86400.0;
        std::vector<double> R0, J;
        sim.residualAndJacobian(dt, R0, &J);
        const bool injBhp = sim.wellOnBhp(0), prdBhp = sim.wellOnBhp(1);
        const double bhpI = sim.wellBhp(0);
        const int m = int(R0.size());
        double err = 0, scale = 0;
        for (int col = 0; col < m; ++col) {
            const int q = col / nv, k = col % nv;
            double Nt = 0; for (int i = 0; i < ncw; ++i) Nt += N[size_t(q) * ncw + i];
            const double h = k == 0 ? p[size_t(q)] * 1e-7 : Nt * 1e-7;
            std::vector<double> pp = p, NN = N, pm = p, Nm = N;
            if (k == 0) { pp[size_t(q)] += h; pm[size_t(q)] -= h; } else { NN[size_t(q) * ncw + k - 1] += h; Nm[size_t(q) * ncw + k - 1] -= h; }
            std::vector<double> Rp, Rm;
            sim.setState(pp, NN); sim.residualAndJacobian(dt, Rp, nullptr);
            sim.setState(pm, Nm); sim.residualAndJacobian(dt, Rm, nullptr);
            for (int r = 0; r < m; ++r) {
                const double fd = (Rp[size_t(r)] - Rm[size_t(r)]) / (2 * h), an = J[size_t(r) * m + col];
                err = std::max(err, std::fabs(fd - an) * (k == 0 ? p[size_t(q)] : Nt));
                scale = std::max(scale, std::fabs(an) * (k == 0 ? p[size_t(q)] : Nt));
            }
        }
        std::printf("  tables + Pc + %s, injector %s (BHP %.1f psia), producer %s on %s: Jacobian max scaled error %.2e\n", variant ? "Stone II" : "ECLIPSE kro",
                    injBhp ? "BHP-limited" : "rate", bhpI / 6894.757, variant ? "LRAT" : "ORAT", prdBhp ? "BHP" : "rate", err / scale);
        check(err / scale < 1e-5 && injBhp == (variant == 1) && !prdBhp,
              variant ? "Stone II + Pc tables, BHP-limited two-connection injector, LRAT producer: Jacobian matches finite differences"
                      : "SWFN/SGFN/SOF3 + Pc, general connections (NNC), implicit injector, ORAT producer with heads: Jacobian matches finite differences");
    }
    {   // run: injector switches to its BHP limit; ORAT producer meets its stock-tank target through the separator
        CompCase c; c.nx = 12; c.ny = 1; c.nz = 1; c.dx = {20}; c.dy = {20}; c.dz = {5}; c.poro = {0.2};
        c.permx = c.permy = c.permz = {2e-13}; c.top = 2500; c.gravity = false; c.fluid = spe5Fluid(); c.T = T;
        c.p0 = {psiaToPa(3000)}; c.z0 = {zoil}; c.separator = {{300.0, 7e5}, {288.71, 101325.0}};
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {2e-12}; inj.zinj = zgas; inj.rate = 200.0; inj.bhpMax = psiaToPa(3300);
        CompWell prd; prd.cells = {c.nx - 1}; prd.wi = {2e-12}; prd.bhp = psiaToPa(1500); prd.control = 3; prd.target = 2.0 / 86400.0;   // 2 sm3/d
        c.wells = {inj, prd};
        CompSimulator sim; sim.setup(c);
        const bool ok = sim.advanceTo(40 * 86400.0, 5 * 86400.0);
        const CompWellReport& I = sim.wellReports()[0]; const CompWellReport& P = sim.wellReports()[1];
        double vo, vg; sim.surfaceVolumes(zoil, vo, vg);
        std::printf("  run: injector BHP %.2f psia (limit 3300, on BHP %s), gas %.0f sm3/d; producer ORAT %.4f sm3/d (target 2), GOR %.1f sm3/sm3 (oil %.1f), FPR %.1f psia\n",
                    I.bhp / 6894.757, I.onBhp ? "yes" : "no", I.qgi * 86400, P.qo * 86400, P.qg / std::max(1e-30, P.qo), vg / vo, sim.hcpvPressure() / 6894.757);
        check(ok && I.onBhp && std::fabs(I.bhp - psiaToPa(3300)) < 1.0 && I.qgi > 0 && !P.onBhp && std::fabs(P.qo * 86400 / 2.0 - 1) < 0.02,
              "injector switches to its BHP limit; ORAT producer holds the stock-tank rate (separator factors)");
    }
    {   // waterflood of the undersaturated SPE5 oil vs the Buckley-Leverett solution (breakthrough time)
        CompCase c; c.nx = 100; c.ny = 1; c.nz = 1; c.dx = {1}; c.dy = {1}; c.dz = {1}; c.poro = {0.2};
        c.permx = c.permy = c.permz = {1e-11}; c.gravity = false; c.fluid = spe5Fluid(); c.T = T;
        c.p0 = {psiaToPa(3500)}; c.z0 = {zoil}; c.water = true; c.swc = 0.2; c.sorw = 0.2; c.nw = 2; c.now = 2; c.muw = 5e-4; c.cw = 0;
        double vw, mo;
        {
            const TParams tp = atTemperature(c.fluid, T);
            const PhaseEval e = evalPhase(c.fluid, tp, psiaToPa(3500), zoil);
            mo = lbcViscosity(c.fluid, zoil, T, 1.0 / e.vm);
            vw = 18.01528e-3 * c.bwRef / c.rhoWsc;
        }
        const double pv = 100 * 0.2, tPV = 200 * 86400.0;   // one pore volume in 200 days
        CompWell inj; inj.injector = true; inj.water = true; inj.cells = {0}; inj.wi = {1}; inj.rate = pv / vw / tPV;
        CompWell prd; prd.cells = {c.nx - 1}; prd.wi = {1e-10}; prd.bhp = psiaToPa(1000);
        prd.control = 2; prd.target = pv / tPV;   // voidage replacement: reservoir-volume rate = injection rate
        c.wells = {inj, prd};
        CompSimulator sim; sim.setup(c);
        // Buckley-Leverett: fw = 1 / (1 + kro mu_w / (krw mu_o)); front by the Welge tangent
        auto fw = [&](double S) { const double sn = std::clamp((S - 0.2) / 0.6, 0.0, 1.0); const double krw = sn * sn, kro = (1 - sn) * (1 - sn); return krw <= 0 ? 0.0 : 1.0 / (1.0 + kro * c.muw / (krw * mo)); };
        double best = 0, sf = 0;
        for (double S = 0.2005; S < 0.8; S += 1e-4) { const double r = fw(S) / (S - 0.2); if (r > best) { best = r; sf = S; } }
        const double tbtAn = (sf - 0.2) / fw(sf);
        double tbt = -1, lastW = 0, lastI = 0;
        const int nchk = 120;
        for (int k = 1; k <= nchk && tbt < 0; ++k) {
            sim.advanceTo(tPV * 1.2 * k / nchk, 2 * 86400.0);
            const double dW = sim.produced(6) - lastW, dI = sim.injected(6) - lastI;
            lastW = sim.produced(6); lastI = sim.injected(6);
            if (dW > 0.01 * dI) tbt = 1.2 * (k - 0.5) / nchk;
        }
        double pav = 0; for (double v : sim.pressure()) pav += v / c.nx;
        std::printf("  waterflood: mu_o %.3f cP, BL front Sw %.3f, breakthrough %.3f PV analytic vs %.3f PV simulated (%d steps, %d Newton)\n",
                    mo * 1e3, sf, tbtAn, tbt, sim.stats().steps, sim.stats().newton);
        std::printf("    RESV producer: BHP %.1f psia, mean reservoir pressure %.1f psia (initial 3500), on BHP: %s\n", sim.wellBhp(1) / 6894.757, pav / 6894.757, sim.wellOnBhp(1) ? "yes" : "no");
        check(tbt > 0 && std::fabs(tbt / tbtAn - 1) < 0.12, "compositional waterflood breakthrough agrees with Buckley-Leverett");
        check(!sim.wellOnBhp(1) && std::fabs(pav / 6894.757 - 3500) < 20, "RESV-controlled producer keeps voidage replacement (pressure maintained)");
    }
    // 2. slim-tube gas injection (SPE5 oil, SPE5 injection gas) above and below the first-contact miscibility pressure
    auto slim = [&](double ppsia, double& rfHeavy, double& mbErr, CompStats& st, double& sgFront) {
        CompCase c; c.nx = 40; c.ny = 1; c.nz = 1; c.dx = {2.5}; c.dy = {1}; c.dz = {1}; c.poro = {0.2};
        c.permx = c.permy = c.permz = {1e-13}; c.top = 2500; c.gravity = false; c.fluid = spe5Fluid(); c.T = T;
        c.p0 = {psiaToPa(ppsia)}; c.z0 = {zoil};
        const TParams tp = atTemperature(c.fluid, T);
        EvalOptions vo; vo.dn = false; vo.dP = false;
        const double vg = evalPhase(c.fluid, tp, psiaToPa(ppsia), zgas, vo).vm;
        const double hcpv = 100.0 * 0.2, tEnd = 365.0 * 86400.0;
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {1}; inj.rate = 1.2 * hcpv / vg / tEnd; inj.zinj = zgas;
        CompWell prd; prd.cells = {c.nx - 1}; prd.wi = {1e-12}; prd.bhp = psiaToPa(ppsia) - 2e5;
        c.wells = {inj, prd};
        CompSimulator sim; sim.setup(c);
        std::vector<double> in0(6); for (int i = 0; i < 6; ++i) in0[i] = sim.componentInPlace(i);
        const bool ok = sim.advanceTo(tEnd, 5 * 86400.0);
        mbErr = 0;
        for (int i = 0; i < 6; ++i) {
            const double bal = sim.componentInPlace(i) - in0[i] - sim.injected(i) + sim.produced(i);
            mbErr = std::max(mbErr, std::fabs(bal) / (in0[i] + sim.injected(i)));
        }
        rfHeavy = (sim.produced(3) + sim.produced(4) + sim.produced(5)) / (in0[3] + in0[4] + in0[5]);
        st = sim.stats();
        sgFront = sim.gasSaturation()[size_t(c.nx / 2)];
        return ok;
    };
    double rfH, mbH, sgH, rfL, mbL, sgL; CompStats sH, sL;
    const bool okH = slim(4500, rfH, mbH, sH, sgH), okL = slim(2600, rfL, mbL, sL, sgL);
    std::printf("  slim tube 4500 psia: C10+ recovery %.3f  balance %.1e  Sg(mid) %.3f  steps %d newton %d cuts %d lin %d flashes %d\n", rfH, mbH, sgH, sH.steps, sH.newton, sH.cuts, sH.linear, sH.flashes);
    std::printf("  slim tube 2600 psia: C10+ recovery %.3f  balance %.1e  Sg(mid) %.3f  steps %d newton %d cuts %d lin %d\n", rfL, mbL, sgL, sL.steps, sL.newton, sL.cuts, sL.linear);
    check(okH && okL && mbH < 1e-9 && mbL < 1e-9, "slim tube runs to 1.2 HCPV with exact component balance");
    check(rfH > rfL + 0.05 && rfH > 0.8, "miscible (4500 psia) displacement recovers more C10+ than the immiscible one (2600 psia)");
}

// ------------------------------------------------------------------------------------------- compositional decks (C3)

// ---------------------------------------------------------------------------------------------- CO2STORE (milestone C4)
static ressim::comp::CompCase co2Case(int variant) {
    using namespace ressim; using namespace ressim::comp;
    CompCase c;
    c.fluid.comp = {eos::standardComponent("CO2"), eos::standardComponent("H2O")};
    c.co2store = true; c.salinity = 1.0; c.T = 323.15; c.cr = 4.5e-10; c.pref = 2e7;
    c.sat.sw = {0.2, 0.4, 0.6, 0.8, 1.0}; c.sat.krw = {0, 0.05, 0.2, 0.5, 1}; c.sat.pcow = {2e4, 8e3, 4e3, 1e3, 0};
    c.sat.sg = {0, 0.05, 0.3, 0.6, 0.8}; c.sat.krg = {0, 0, 0.12, 0.45, 0.8}; c.sat.pcog = {0, 0, 0, 0, 0};
    if (variant == 0) {   // 4-cell column: CO2-rich top cell, two-phase cells, brine bottom cell
        for (int k = 0; k < 4; ++k) { c.cellPV.push_back(800.0 * (1 + 0.1 * k)); c.cellDepth.push_back(1502.5 + 5 * k); }
        c.conns = {{0, 1, 2e-11}, {1, 2, 1.6e-11}, {2, 3, 2e-11}, {0, 3, 3e-13}};
        c.p0 = {2.0e7, 2.005e7, 2.01e7, 2.015e7};
        c.z0 = {{0.97, 0.03}, {0.5, 0.5}, {0.2, 0.8}, {0.005, 0.995}};
        CompWell inj; inj.injector = true; inj.cells = {0, 1}; inj.wi = {2e-12, 1e-12}; inj.zinj = {1.0, 0.0}; inj.rate = 3.0; inj.bhpMax = 2.2e7;
        CompWell prd; prd.cells = {2, 3}; prd.wi = {1e-12, 2e-12}; prd.bhp = 1.5e7; prd.control = 5; prd.target = 2e-4;
        c.wells = {inj, prd};
    } else {             // 1-D aquifer, CO2 injection at one end, brine production at the other
        const int n = 40;
        for (int k = 0; k < n; ++k) { c.cellPV.push_back(0.25 * 25 * 50 * 10); c.cellDepth.push_back(1500); }
        for (int k = 0; k + 1 < n; ++k) c.conns.push_back({k, k + 1, 1e-12 * 50 * 10 / 25});
        c.p0 = {2e7}; c.z0 = {{1e-8, 1 - 1e-8}};
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {5e-12}; inj.zinj = {1.0, 0.0}; inj.rate = 2.0; inj.bhpMax = 3e7;
        CompWell prd; prd.cells = {n - 1}; prd.wi = {5e-12}; prd.bhp = 2e7;
        c.wells = {inj, prd};
    }
    return c;
}

static void testCo2Store() {
    using namespace ressim; using namespace ressim::comp;
    std::printf("CO2STORE: CO2-brine thermodynamics and the CO2 storage simulator (milestone C4)\n");
    {
        const co2::Fluid w = co2::makeFluid(323.15, 0.0), b = co2::makeFluid(323.15, 1.0);
        const double m0 = co2::co2Molality(w, 2e7), m1 = co2::co2Molality(b, 2e7);
        double x, y; co2::solubility(w, 2e7, x, y);
        const double rho = co2::co2Density(w, 2e7), mu = co2::co2Viscosity(323.15, 784.3), rb = co2::brineDensity(b, 2e7);
        std::printf("  50 C, 200 bar: CO2 solubility %.3f mol/kg (Duan-Sun 1.30), y_H2O %.4f; 1 m NaCl %.3f mol/kg (ratio %.3f)\n", m0, y, m1, m1 / m0);
        std::printf("  CO2 density %.1f kg/m3 (NIST 784.3), viscosity at 784.3 kg/m3 %.4f cP (NIST 0.069); brine (1 m) %.1f kg/m3, water mu %.3f cP\n",
                    rho, mu * 1e3, rb, w.muBrine * 1e3);
        check(std::fabs(m0 - 1.30) < 0.04 && m1 / m0 > 0.78 && m1 / m0 < 0.85 && y > 0.004 && y < 0.009, "Spycher-Pruess / Duan-Sun CO2-brine mutual solubility");
        check(std::fabs(rho / 784.3 - 1) < 0.04 && std::fabs(mu * 1e3 / 0.069 - 1) < 0.05 && std::fabs(rb - 1034) < 6 && std::fabs(w.muBrine * 1e3 - 0.547) < 0.01,
              "CO2 density (Redlich-Kwong), Fenghour viscosity, Batzle-Wang brine");
    }
    {   // analytic Jacobian vs central differences
        CompSimulator sim; sim.setup(co2Case(0));
        const int ncw = sim.ncw(), nv = ncw + 1;
        std::vector<double> p = sim.pressure(), N = sim.moles();
        sim.storeOld();
        for (size_t k = 0; k < N.size(); ++k) N[k] *= 1.0 + 0.002 * double(k % 3);
        sim.setState(p, N);
        const double dt = 86400.0;
        std::vector<double> R0, J;
        sim.residualAndJacobian(dt, R0, &J);
        const int m = int(R0.size());
        double err = 0, scale = 0;
        for (int col = 0; col < m; ++col) {
            const int q = col / nv, k = col % nv;
            double Nt = 0; for (int i = 0; i < ncw; ++i) Nt += N[size_t(q) * ncw + i];
            const double h = k == 0 ? p[size_t(q)] * 1e-7 : Nt * 1e-7;
            std::vector<double> pp = p, NN = N, pm = p, Nm = N;
            if (k == 0) { pp[size_t(q)] += h; pm[size_t(q)] -= h; } else { NN[size_t(q) * ncw + k - 1] += h; Nm[size_t(q) * ncw + k - 1] -= h; }
            std::vector<double> Rp, Rm;
            sim.setState(pp, NN); sim.residualAndJacobian(dt, Rp, nullptr);
            sim.setState(pm, Nm); sim.residualAndJacobian(dt, Rm, nullptr);
            for (int r = 0; r < m; ++r) {
                const double fd = (Rp[size_t(r)] - Rm[size_t(r)]) / (2 * h), an = J[size_t(r) * m + col];
                err = std::max(err, std::fabs(fd - an) * (k == 0 ? p[size_t(q)] : Nt));
                scale = std::max(scale, std::fabs(an) * (k == 0 ? p[size_t(q)] : Nt));
            }
        }
        std::printf("  4-cell column (phases %d %d %d %d), CO2 injector + brine WRAT producer: Jacobian vs FD max scaled error %.2e\n",
                    sim.phasesIn(0), sim.phasesIn(1), sim.phasesIn(2), sim.phasesIn(3), err / scale);
        check(sim.phasesIn(1) == 2 && sim.phasesIn(3) == 1 && err / scale < 1e-5, "CO2STORE analytic Jacobian (solubility, densities, Fenghour, Pc, wells) matches finite differences");
    }
    {   // 1-D injection: CO2 mass balance, dissolution, trapping inventory
        CompSimulator sim; sim.setup(co2Case(1));
        const double co2In0 = sim.componentInPlace(0);
        bool ok = true;
        for (int m = 1; m <= 12 && ok; ++m) ok = sim.advanceTo(m * 30 * 86400.0);
        double mob, imm, dis; sim.co2Inventory(mob, imm, dis);
        const double inj = sim.injected(0), prod = sim.produced(0), inPlace = sim.componentInPlace(0);
        const double bal = (inPlace - co2In0 - inj + prod) / inj;
        const auto sg = sim.gasSaturation();
        int front = 0; for (size_t q = 0; q < sg.size(); ++q) if (sg[q] > 1e-3) front = int(q);
        const double tot = (mob + imm + dis), injKg = inj * co2::MW_CO2;
        std::printf("  1-D aquifer, 360 d: injected %.4g t CO2, mobile %.4g t, immobile %.4g t, dissolved %.4g t (%.1f %%), produced %.3g t; gas front at cell %d; balance %.1e\n",
                    injKg / 1e3, mob / 1e3, imm / 1e3, dis / 1e3, 100 * dis / tot, prod * co2::MW_CO2 / 1e3, front, bal);
        std::printf("        %d steps, %d Newton, %d cuts; brine production %.4g sm3\n", sim.stats().steps, sim.stats().newton, sim.stats().cuts, sim.wellReports()[1].cumW);
        check(ok && std::fabs(bal) < 1e-8 && dis > 0.02 * tot && front > 3 && front < 39 && std::fabs(tot - (inPlace - co2In0) * co2::MW_CO2 - co2In0 * co2::MW_CO2) < 1e-6 * tot,
              "CO2STORE injection: CO2 mass balance, dissolution into brine, plume advance, trapping inventory");
    }
    // OPM CO2STORE deck (OIL = brine + GAS): import, run, export both ways, GAS + WATER variant
    const std::string d = findDeck("examples/CO2STORE.DATA");
    if (d.empty()) { std::printf("  (examples/CO2STORE.DATA not found - skipped)\n"); return; }
    const Case c0 = importEclipseDeck(d);
    check(c0.comp.enabled && c0.comp.co2store && std::fabs(c0.comp.fluid.salinity - 0.7) < 1e-12 && std::fabs(c0.comp.fluid.resT - 323.15) < 1e-9 &&
          c0.engineUnsupported.empty(), "OPM CO2STORE deck (OIL + GAS + CO2STORE, SGOF, SALINITY, RTEMPVD) imports as a CO2 storage case");
    auto run = [&](const Case& c, int steps, std::vector<comp::CompSummary>& out, DeckModel* mm = nullptr) {
        DeckModel m; std::string err;
        if (!comp::deckModel(c, m, &err)) { std::printf("  error: %s\n", err.c_str()); return false; }
        const bool ok = comp::runDeckModel(m, out, {}, &err, steps);
        if (!ok) std::printf("  run error: %s\n", err.c_str());
        if (mm) *mm = m;
        return ok;
    };
    std::vector<comp::CompSummary> o1;
    DeckModel m1;
    const bool ok1 = run(c0, 14, o1, &m1);
    if (ok1) {
        const auto& a = o1[12]; const auto& b = o1.back();
        const double gi = a.fgit, gip = a.fgip - o1[0].fgip;
        const double mass = a.co2Mob + a.co2Immob + a.co2Diss, injMass = gi * co2::gasSurfaceVolume(m1.base.stdT, m1.base.stdP) > 0 ? gi / co2::gasSurfaceVolume(m1.base.stdT, m1.base.stdP) * co2::MW_CO2 : 0;
        std::printf("  CO2STORE.DATA day %.0f: FGIT %.4g sm3, FGIP change %.4g sm3, FPR %.2f bar; CO2 mobile %.4g t, trapped %.4g t, dissolved %.4g t (injected %.4g t); day %.0f dissolved %.4g t\n",
                    a.time, gi, gip, a.fpr, a.co2Mob / 1e3, a.co2Immob / 1e3, a.co2Diss / 1e3, injMass / 1e3, b.time, b.co2Diss / 1e3);
        check(std::fabs(gip / gi - 1) < 1e-6 && std::fabs(mass / injMass - 1) < 1e-4 && a.co2Diss > 0.1 * mass && b.co2Diss > a.co2Diss && a.fgir > 0 && b.fgir == 0,
              "CO2STORE deck run: FGIP balance, trapping inventory = injected mass, continued dissolution after shut-in");
    }
    {   // export verbatim and structured; GAS + WATER form of the same model (SWFN / SGFN, GWC from EQUIL item 3)
        Case cs = c0; cs.comp.structured = true;
        const std::string tv = eclipseDeckToString(c0), ts = eclipseDeckToString(cs);
        const Case cv = eclipseDeckFromString(tv, ""), cS = eclipseDeckFromString(ts, "");
        const bool txtOk = ts.find("CO2STORE") != std::string::npos && ts.find("SALINITY") != std::string::npos && ts.find("RTEMPVD") != std::string::npos &&
                           ts.find("PVTW") == std::string::npos && ts.find("PVDO") == std::string::npos && tv.find("PVDO") == std::string::npos;
        Case gw = c0; gw.phases.oil = false; gw.phases.water = true; gw.comp.structured = true;
        for (auto& r : gw.props.swof) (void)r;
        {   // SGOF (oil = brine) -> SWOF (Sw = 1 - Sg, krw = krog) + SGOF (krg)
            std::vector<SatRow> sw;
            for (size_t i = gw.props.sgof.size(); i-- > 0;) sw.push_back({1 - gw.props.sgof[i].s, gw.props.sgof[i].kro, 0, 0});
            gw.props.swof = sw;
            gw.equil.owc = gw.equil.goc;
        }
        const std::string tg = eclipseDeckToString(gw);
        const Case cg = eclipseDeckFromString(tg, "");
        std::vector<comp::CompSummary> ov, oS, og;
        const bool r = cv.comp.co2store && cS.comp.co2store && cg.comp.co2store && run(cv, 3, ov) && run(cS, 3, oS) && run(cg, 3, og);
        bool same = r && txtOk;
        if (same) for (size_t i = 0; i < 4; ++i) same = same && std::fabs(ov[i].fpr - o1[i].fpr) < 1e-9 * o1[i].fpr && std::fabs(oS[i].fpr - o1[i].fpr) < 1e-6 * o1[i].fpr &&
                                                    std::fabs(og[i].fpr - o1[i].fpr) < 1e-6 * o1[i].fpr && std::fabs(og[i].co2Diss - o1[i].co2Diss) < 1e-5 * (1 + o1[i].co2Diss);
        std::printf("  export: verbatim / structured / GAS + WATER (SWFN + SGFN) decks re-run, FPR day 91 %.6f / %.6f / %.6f (original %.6f)\n",
                    r ? ov[3].fpr : 0.0, r ? oS[3].fpr : 0.0, r ? og[3].fpr : 0.0, o1.size() > 3 ? o1[3].fpr : 0.0);
        check(same, "CO2STORE deck two-way: verbatim and structured export, GAS + WATER form (SWFN / SGFN) give the same results");
    }
    {   // built-in examples: CO2 storage (10 y injection, 4 steps) and SPE5-fluid lean-gas injection
        for (const char* id : {"co2", "comp"}) {
            const Case c = makeExample(id);
            comp::DeckModel m; std::string err; std::vector<comp::CompSummary> out;
            const bool ok = validateCase(c).empty() && comp::deckModel(c, m, &err) && comp::runDeckModel(m, out, {}, &err, std::string(id) == "co2" ? 4 : 8);
            if (std::string(id) == "co2") {
                const comp::CompSummary& e = out.back();
                const double tot = e.co2Mob + e.co2Immob + e.co2Diss;                                                   // kg
                const double injKg = e.fgit * co2::MW_CO2 / co2::gasSurfaceVolume(c.comp.fluid.stdT, c.comp.fluid.stdP);   // sm3 -> kg
                std::printf("  example co2: %zu steps, injected %.4g t, in place %.4g t (mobile %.3g, trapped %.3g, dissolved %.3g)\n",
                            out.size(), injKg / 1e3, tot / 1e3, e.co2Mob / 1e3, e.co2Immob / 1e3, e.co2Diss / 1e3);
                check(ok && std::fabs(tot - injKg) < 0.02 * injKg && e.co2Diss > 0.1 * tot, "CO2 storage example: in-place CO2 = injected CO2, > 10% dissolved after 4 years");
            } else {
                std::printf("  example comp: %zu steps, FOPT %.4g, FPR %.1f\n", out.size(), ok ? out.back().fopt : 0.0, ok ? out.back().fpr : 0.0);
                check(ok && out.size() == 9 && std::fabs(out.back().fopt - 3000 * 8 * 91.3125) < 1, "compositional example (SPE5 fluid): ORAT 3000 held for 2 years");
            }
        }
    }
}

static void testCompDeck() {
    using namespace ressim::comp;
    std::printf("Compositional E300 decks (EQUIL + ZMFVD initialisation, schedule, separator, WELLSTRE / WINJGAS)\n");
    // 1. hydrostatic equilibrium: oil column with a capillary transition zone across the OWC, no wells, 2 years
    const std::string eq = R"(RUNSPEC
DIMENS
 2 1 12 /
OIL
WATER
GAS
COMPS
 6 /
METRIC
START
 1 'JAN' 2020 /
GRID
DX
 24*50 /
DY
 24*50 /
DZ
 24*4 /
TOPS
 2*2500 /
PORO
 24*0.25 /
PERMX
 24*200 /
PERMY
 24*200 /
PERMZ
 24*20 /
PROPS
EOS
 PR /
PRCORR
CNAMES
 C1 C3 C6 C10 C15 C20 /
TCRIT
 190.56 369.83 507.44 617.67 705.56 766.67 /
PCRIT
 46.04 42.49 30.12 20.96 13.79 11.17 /
ZCRIT
 0.290 0.277 0.264 0.257 0.245 0.235 /
ACF
 0.0130 0.1524 0.3007 0.4885 0.6500 0.8500 /
MW
 16.04 44.10 86.18 142.29 206.0 282.0 /
BIC
 0.0
 0.0 0.0
 0.0 0.0 0.0
 0.05 0.005 0.0 0.0
 0.05 0.005 0.0 0.0 0.0 /
RTEMP
 71 /
ZMFVD
 2500 0.5 0.03 0.07 0.2 0.15 0.05
 2550 0.5 0.03 0.07 0.2 0.15 0.05 /
PVTW
 250 1.02 4.5E-5 0.5 0 /
ROCK
 250 5E-5 /
DENSITY
 800 1020 0.9 /
SWFN
 0.2 0 2.5
 0.4 0.05 0.9
 0.6 0.2 0.3
 0.8 0.5 0.05
 1.0 1.0 0 /
SGFN
 0 0 0
 0.8 1 0 /
SOF3
 0 0 0
 0.2 0 0
 0.8 1 1 /
SOLUTION
EQUIL
 2520 280 2530 0 2400 0 /
SCHEDULE
TSTEP
 365 365 /
END
)";
    {
        std::vector<std::string> warn;
        const ressim::deck::RawDeck D = ressim::deck::parseDeckText(eq, "");
        const Case c = eclipseDeckFromString(eq, "", &warn);
        DeckModel m; std::string err;
        const bool ok = deckModel(D, c, m, &err);
        if (!ok) std::printf("  error: %s\n", err.c_str());
        check(ok && m.base.size() == 24 && m.base.water && m.stepDays.size() == 2, "E300 deck with EQUIL / ZMFVD / SWFN-SGFN-SOF3 builds a compositional model");
        if (ok) {
            CompSimulator sim; sim.setup(m.base);
            const std::vector<double> p0 = sim.pressure(), sw0 = sim.waterSaturation();
            double swTop = 1, swBot = 0; for (double v : sw0) { swTop = std::min(swTop, v); swBot = std::max(swBot, v); }
            sim.advanceTo(730 * 86400.0, 30 * 86400.0);
            const std::vector<double> p1 = sim.pressure(), sw1 = sim.waterSaturation();
            double dp = 0, ds = 0; for (size_t q = 0; q < p0.size(); ++q) { dp = std::max(dp, std::fabs(p1[q] - p0[q])); ds = std::max(ds, std::fabs(sw1[q] - sw0[q])); }
            // datum pressure check: interpolate the initial oil pressure at the datum depth (2520 m: between layers 5 and 6)
            const double pd = 0.5 * (p0[2 * 4] + p0[2 * 5]);
            std::printf("  equilibrium: Sw %.3f (top) .. %.3f (bottom), p(datum) %.4f bar (280), after 2 years max |dp| %.2e bar, max |dSw| %.2e (%d steps)\n",
                        swTop, swBot, pd / 1e5, dp / 1e5, ds, sim.stats().steps);
            check(std::fabs(pd / 1e5 - 280) < 0.05 && swTop > 0.3 && swTop < 0.45 && swBot > 0.99 && dp < 50.0 && ds < 2e-4,
                  "EQUIL + ZMFVD hydrostatic initialisation is in equilibrium (capillary transition zone, OWC)");
        }
    }
    // 2. SPE5 compositional WAG deck (Killough & Kossack 1987): depletion and the first water / gas cycles
    const std::string d = findDeck("decks/spe5/SPE5_E300.DATA");
    if (d.empty()) { std::printf("  (decks/spe5/SPE5_E300.DATA not found - skipped)\n"); return; }
    DeckModel m; std::string err;
    const bool ok = deckModelFromFile(d, m, &err);
    if (!ok) std::printf("  error: %s\n", err.c_str());
    check(ok && m.base.size() == 147 && m.wellNames.size() == 3 && m.stepDays.size() == 240 && m.base.separator.size() == 1 &&
          !m.stepWells[0].empty() && m.stepWells[36].size() == 3 && m.stepWells[36][2].open && m.stepWells[36][2].zinj.size() == 6 &&
          std::fabs(m.stepWells[36][2].zinj[0] - 0.77) < 1e-12, "SPE5 E300 deck: grid, wells, WAG schedule, WELLSTRE / WINJGAS stream, FIELDSEP");
    if (!ok) return;
    {   // two-way deck: export the imported case and rebuild the compositional model from the exported text
        const Case c0 = importEclipseDeck(d);
        const std::string txt = eclipseDeckToString(c0);
        DeckModel m2; std::string e2;
        const bool ok2 = deckModel(ressim::deck::parseDeckText(txt, ""), eclipseDeckFromString(txt, ""), m2, &e2);
        if (!ok2) std::printf("  re-import error: %s\n", e2.c_str());
        bool same = ok2 && m2.base.size() == m.base.size() && m2.stepDays.size() == m.stepDays.size() && m2.base.separator.size() == 1 &&
                    m2.stepWells[36][2].zinj == m.stepWells[36][2].zinj && std::fabs(m2.base.T - m.base.T) < 1e-9;
        if (same) for (int q = 0; q < m.base.size(); ++q) same = same && std::fabs(m2.base.p0[size_t(q)] - m.base.p0[size_t(q)]) < 1e-6 * m.base.p0[size_t(q)];
        check(same, "SPE5 E300 deck survives export -> re-import (EOS, ZMFVD, FIELDSEP, WELLSTRE / WINJGAS, equilibration)");
        // structured export (as after UI edits): the compositional keywords are written from Case::comp
        Case cs = c0; cs.comp.structured = true;
        const std::string ts = eclipseDeckToString(cs);
        const Case c3 = eclipseDeckFromString(ts, "");
        DeckModel m3; std::string e3;
        const bool ok3 = c3.comp.enabled && deckModel(c3, m3, &e3);
        if (!ok3) std::printf("  structured re-import error: %s\n", e3.c_str());
        int nKw = 0; for (const char* k : {"COMPS", "EOS", "TCRIT", "BIC", "ZMFVD", "FIELDSEP", "WELLSTRE", "WINJGAS"}) nKw += ts.find(std::string("\n") + k + "\n") != std::string::npos;
        bool same3 = ok3 && nKw == 8 && ts.find("PVTO") == std::string::npos && m3.base.fluid.nc() == 6 && c3.comp.streams.size() == c0.comp.streams.size();
        if (same3) {
            for (int i = 0; i < 6; ++i) {
                const auto &x = m3.base.fluid.comp[size_t(i)], &y = m.base.fluid.comp[size_t(i)];
                same3 = same3 && std::fabs(x.tc / y.tc - 1) < 1e-9 && std::fabs(x.pc / y.pc - 1) < 1e-9 && std::fabs(x.omega - y.omega) < 1e-12 && std::fabs(x.mw / y.mw - 1) < 1e-9;
                for (int j = 0; j < 6; ++j) same3 = same3 && std::fabs(m3.base.fluid.k(i, j) - m.base.fluid.k(i, j)) < 1e-12;
            }
            same3 = same3 && m3.stepWells[36][2].zinj == m.stepWells[36][2].zinj && m3.base.separator.size() == 1;
            for (int q = 0; q < m.base.size(); ++q) same3 = same3 && std::fabs(m3.base.p0[size_t(q)] - m.base.p0[size_t(q)]) < 1e-6 * m.base.p0[size_t(q)];
        }
        check(same3, "structured compositional export (COMPS / EOS / TCRIT / BIC / ZMFVD / FIELDSEP / WELLSTRE / WINJGAS from the case) re-imports identically");
    }
    std::vector<CompSummary> out; CompStats st;
    const bool run = runDeckModel(m, out, {}, &err, 48, &st);   // 4 years: 2 years depletion, 1 water, 1 gas
    const CompSummary& a = out[3]; const CompSummary& b = out[30]; const CompSummary& e = out.back();
    // component balance check via in-place volumes is implicit in the core; here: controls and phase behaviour
    std::printf("  SPE5: day %.0f FOPR %.1f GOR %.3f Mscf/stb (above Pb), FPR %.1f; day %.0f FWIR %.0f; day %.0f FGIR %.0f FPR %.1f FOPT %.4g FGPT %.4g FWPT %.4g\n",
                a.time, a.fopr, a.fgpr / std::max(1e-9, a.fopr), a.fpr, b.time, b.fwir, e.time, e.fgir, e.fpr, e.fopt, e.fgpt, e.fwpt);
    std::printf("        %d steps, %d Newton, %d linear, %d cuts\n", st.steps, st.newton, st.linear, st.cuts);
    check(run && std::fabs(a.fopr / 12000 - 1) < 1e-3 && std::fabs(a.fgpr / a.fopr - 0.556) < 0.03 && std::fabs(b.fwir / 12000 - 1) < 1e-3 &&
          std::fabs(e.fgir / 12000 - 1) < 1e-3 && e.fpr > 1000 && st.cuts < 10,
          "SPE5 compositional WAG runs: ORAT target, separator GOR, water and gas injection rates");
    {   // water vaporisation (VAPWAT / '-- RESSIM: WATERVAP'): water vapour in the gas phase, water material balance, round trip
        using namespace ressim;
        Case c = makeExample("comp");
        c.comp.fluid.resT = 400.0; c.comp.fluid.tempvd.clear();   // hot reservoir: y_w = Psat / p ~ 1 %
        c.comp.waterVap = true;
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        double prodW[2] = {0, 0}, vap[2] = {0, 0}, bal = 1; bool ran[2] = {false, false};
        for (int v = 0; v < 2; ++v) {
            Case cc = c; cc.comp.waterVap = v == 1;
            comp::DeckModel m; std::string err;
            if (!comp::deckModel(cc, m, &err)) { std::printf("  WATERVAP: %s\n", err.c_str()); continue; }
            const int nc = m.base.fluid.nc();
            double w0 = -1;
            std::vector<comp::CompSummary> out; comp::CompStats st;
            ran[v] = comp::runDeckModel(m, out, [&](const comp::CompSummary&, const comp::CompSimulator& sim) {
                if (w0 < 0) w0 = sim.componentInPlace(nc);
                prodW[v] = sim.produced(nc); vap[v] = sim.vaporisedWater();
                if (v == 1) bal = (sim.componentInPlace(nc) - w0 + sim.produced(nc) - sim.injected(nc)) / std::max(1.0, w0);
                return true;
            }, &err, 8, &st);
        }
        std::printf("  water vaporisation (400 K): water in the gas phase %.4g kmol (off %.4g), produced water %.4g kmol (off %.4g), water balance %.1e, round trip %d\n",
                    vap[1] / 1e3, vap[0] / 1e3, prodW[1] / 1e3, prodW[0] / 1e3, bal, int(b.comp.waterVap));
        check(ran[0] && ran[1] && std::fabs(bal) < 1e-8 && vap[1] > 0 && vap[0] == 0 && b.comp.waterVap,
              "compositional water vaporisation: water carried in the gas phase, water balance closes, deck round trip");
    }
}

// ---------------------------------------------------------------------------------------------- C5: trapping mechanisms
static void testTrapping() {
    using namespace ressim; using namespace ressim::comp;
    std::printf("Trapping mechanisms and interactions (round C5)\n");
    {   // Land / Carlson / Killough scanning curves
        bool ok = true; double worst = 0;
        for (int model = 0; model < 2; ++model) {
            CompCase c = co2Case(1); c.hyst.on = true; c.hyst.model = model; c.trap.sgtMax = 0.3;
            double t = 0;
            gasRelPerm(c, 0, -1, 0.8, 0.8, &t);
            ok = ok && std::fabs(t - 0.3) < 1e-12;                                    // Land at the end point
            const double Cl = 1 / (0.3 - 0.05) - 1 / (0.8 - 0.05);
            for (double sh : {0.2, 0.4, 0.6}) {
                double st = 0; gasRelPerm(c, 0, -1, sh, sh, &st);
                const double land = 0.05 + (sh - 0.05) / (1 + Cl * (sh - 0.05));
                worst = std::max(worst, std::fabs(st - land));
                const double kd = gasRelPerm(c, 0, -1, sh, 0.0), ks = gasRelPerm(c, 0, -1, sh - 1e-7, sh);
                ok = ok && std::fabs(ks - kd) < 1e-5 && gasRelPerm(c, 0, -1, st, sh) == 0.0;
                for (double f : {0.3, 0.6, 0.9}) { const double s = st + f * (sh - st); ok = ok && gasRelPerm(c, 0, -1, s, sh) <= gasRelPerm(c, 0, -1, s, 0.0) + 1e-12; }
            }
        }
        std::printf("  Land trapped saturation: max deviation %.2e; scanning curves continuous at Sg,max, zero at Sgt, below drainage (Carlson + Killough)\n", worst);
        check(ok && worst < 1e-12, "hysteresis: Land trapping relation and Carlson / Killough scanning curves");
    }
    auto chase = [&](bool hyst, double& resid, double& mob, double& bal, double& sgt0) {
        CompCase c = co2Case(1);
        c.hyst.on = hyst; c.hyst.model = 0; c.trap.sgtMax = 0.3;
        CompWell brine; brine.injector = true; brine.cells = {0}; brine.wi = {5e-12}; brine.zinj = {0.0, 1.0}; brine.rate = 20.0; brine.bhpMax = 3e7; brine.open = false;
        c.wells.push_back(brine);
        CompSimulator sim; sim.setup(c);
        const double n0 = sim.componentInPlace(0);
        bool ok = true;
        for (int m = 1; m <= 12 && ok; ++m) ok = sim.advanceTo(m * 30 * 86400.0);
        sim.wells()[0].open = false; sim.wells()[2].open = true;
        for (int m = 13; m <= 16 && ok; ++m) ok = sim.advanceTo(m * 30 * 86400.0);
        const Co2Inventory I = sim.co2Trapping();
        resid = I.residual; mob = I.mobile + I.structural;
        const double inPlace = sim.componentInPlace(0);
        bal = std::fabs(I.total() / co2::MW_CO2 - inPlace) / inPlace + std::fabs(inPlace - n0 - sim.injected(0) + sim.produced(0)) / sim.injected(0);
        sgt0 = sim.trappedGasSaturation()[1];
        return ok;
    };
    {
        double r0, m0, b0, s0, r1, m1, b1, s1;
        const bool ok = chase(false, r0, m0, b0, s0) && chase(true, r1, m1, b1, s1);
        std::printf("  CO2 injection 360 d + brine chase 120 d: residual %.4g t / mobile %.4g t without hysteresis, %.4g t / %.4g t with Land (Sgt,max 0.3); "
                    "Sgt cell 2 %.3f; balance %.1e / %.1e\n", r0 / 1e3, m0 / 1e3, r1 / 1e3, m1 / 1e3, s1, b0, b1);
        check(ok && r1 > 1.3 * r0 && m1 < m0 && b0 < 1e-7 && b1 < 1e-7 && s1 > 0.05, "residual trapping: imbibition traps CO2 behind the plume (hysteresis), mass balance of the inventory");
    }
    {   // molecular diffusion against the semi-infinite analytic solution
        CompCase c = co2Case(1);
        c.wells.clear(); c.cellPV.clear(); c.cellDepth.clear(); c.conns.clear(); c.z0.clear();
        const int n = 40; const double dx = 0.05, A = 1.0, phi = 0.25;
        for (int k = 0; k < n; ++k) { c.cellPV.push_back(phi * A * dx); c.cellDepth.push_back(1500); c.z0.push_back(k < n / 2 ? std::vector<double>{0.01, 0.99} : std::vector<double>{1e-9, 1 - 1e-9}); }
        for (int k = 0; k + 1 < n; ++k) c.conns.push_back({k, k + 1, 1e-13 * A / dx, phi * A / dx});
        c.trap.diffusion = true; c.trap.dCo2Water = 2e-9; c.trap.dH2oWater = 2e-9;
        CompSimulator sim; sim.setup(c);
        const double n0 = sim.componentInPlace(0);
        double left0 = 0; for (int k = 0; k < n / 2; ++k) left0 += sim.moles()[size_t(k) * sim.ncw()];
        const double tEnd = 365.25 * 86400.0;
        bool ok = sim.advanceTo(tEnd, 86400.0 * 5);
        double right = 0; for (int k = n / 2; k < n; ++k) right += sim.moles()[size_t(k) * sim.ncw()];
        const double analytic = left0 / (n / 2 * dx) * std::sqrt(c.trap.dCo2Water * tEnd / M_PI);
        CompCase c2 = c; c2.trap.diffusion = false;
        CompSimulator s2; s2.setup(c2); ok = ok && s2.advanceTo(tEnd, 86400.0 * 30);
        double right2 = 0; for (int k = n / 2; k < n; ++k) right2 += s2.moles()[size_t(k) * s2.ncw()];
        const double cons = std::fabs(sim.componentInPlace(0) - n0) / n0;
        std::printf("  diffusion, 1 year: CO2 crossed %.4g mol (analytic sqrt(Dt/pi) %.4g mol, %+.1f %%), without diffusion %.2e mol; conservation %.1e\n",
                    right, analytic, 100 * (right / analytic - 1), right2, cons);
        check(ok && std::fabs(right / analytic - 1) < 0.06 && right2 < 1e-3 * right && cons < 1e-9, "molecular diffusion: Fickian spreading matches the analytic solution, conservative");
    }
    {   // salt precipitation (dry-out) near a CO2 injector; closed system -> NaCl conserved
        CompCase c = co2Case(1);
        c.wells.clear(); c.cellPV.clear(); c.cellDepth.clear(); c.conns.clear();
        const int n = 20;
        for (int k = 0; k < n; ++k) { c.cellPV.push_back(k + 1 < n ? 0.5 * std::pow(1.6, k) : 1e7); c.cellDepth.push_back(1500); }
        for (int k = 0; k + 1 < n; ++k) c.conns.push_back({k, k + 1, 1e-11, 0.0});
        c.salinity = 4.0; c.trap.saltPrecip = true;
        for (double& v : c.sat.pcow) v = 0;   // no capillary back-flow of brine: classical dry-out
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {2e-11}; inj.zinj = {1.0, 0.0}; inj.rate = 0.5; inj.bhpMax = 4e7;
        c.wells = {inj};
        CompSimulator sim; sim.setup(c);
        double d0, s0; sim.saltInPlace(d0, s0);
        bool ok = true;
        for (int m = 1; m <= 12 && ok; ++m) ok = sim.advanceTo(m * 30 * 86400.0);
        double d1, s1; sim.saltInPlace(d1, s1);
        const auto& pf = sim.porosityFactor(); const auto& kf = sim.permFactorCells();
        const Co2Inventory I = sim.co2Trapping();
        std::printf("  salt precipitation, 360 d dry CO2 at 0.5 mol/s: solid NaCl %.4g kg (%.3g %% of the salt), phi/phi0 cell 1 %.4f, k/k0 %.4f; NaCl balance %.1e\n",
                    I.saltSolid, 100 * s1 / (d1 + s1), pf[0], kf[0], std::fabs(d1 + s1 - d0 - s0) / (d0 + s0));
        check(ok && s1 > 0 && pf[0] < 1 && kf[0] < pf[0] && std::fabs(d1 + s1 - d0 - s0) < 1e-9 * (d0 + s0) && s0 == 0,
              "salt precipitation: dry-out precipitates NaCl, porosity / permeability reduction, NaCl conserved");
    }
    {   // mineral trapping: CO2 balance including the mineralised CO2, porosity change
        CompCase c = co2Case(1);
        co2::MineralReaction mr = co2::mineralPresets().front(); mr.area = 2e4;
        c.trap.minerals = {mr};
        CompSimulator sim; sim.setup(c);
        const double n0 = sim.componentInPlace(0);
        bool ok = true;
        for (int m = 1; m <= 12 && ok; ++m) ok = sim.advanceTo(m * 30 * 86400.0);
        const Co2Inventory I = sim.co2Trapping();
        const double minMol = I.mineral / co2::MW_CO2;
        const double bal = (sim.componentInPlace(0) + minMol - n0 - sim.injected(0) + sim.produced(0)) / sim.injected(0);
        const double tot = (I.total() - I.saltSolid) / co2::MW_CO2;
        double pmin = 1; for (double v : sim.porosityFactor()) pmin = std::min(pmin, v);
        std::printf("  mineral trapping (%s, A = %.0f m2/m3): %.4g t CO2 mineralised (%.2f %% of injected), min phi/phi0 %.5f; balance %.1e, inventory %.1e\n",
                    mr.name.c_str(), mr.area, I.mineral / 1e3, 100 * minMol / sim.injected(0), pmin, bal,
                    std::fabs(tot - sim.componentInPlace(0) - minMol) / tot);
        check(ok && I.mineral > 0 && std::fabs(bal) < 1e-8 && pmin < 1 && std::fabs(tot - sim.componentInPlace(0) - minMol) < 1e-8 * tot,
              "mineral trapping: kinetic carbonation, CO2 balance incl. minerals, porosity change");
    }
    {   // secondary producer limits and separator trains (compositional example)
        Case c = makeExample("comp");
        auto runC = [&](const Case& cc, std::vector<CompSummary>& out) { DeckModel m; std::string e; return comp::deckModel(cc, m, &e) && comp::runDeckModel(m, out, {}, &e, 3); };
        std::vector<CompSummary> o0, o1, o2;
        bool ok = runC(c, o0);
        Case cl = c;
        const double glim = ok ? 0.6 * o0[1].fgpr : 1.0;
        for (auto& st : cl.schedule) for (auto& wc : st.controls) if (wc.mode == Control::ORAT) wc.grat = glim;
        ok = ok && runC(cl, o1);
        Case cs = c;
        CompositionalInput::SepTrain tr; tr.name = "LP"; tr.stages = {{288.71, 101325.0}};
        cs.comp.sepTrains = {tr};
        for (auto& w : cs.wells) if (w.type == WellType::Producer) cs.comp.wellSep[w.name] = "LP";
        ok = ok && runC(cs, o2);
        // export / import of SEPCOND + WSEPCOND + limits
        Case ce = cs; ce.comp.structured = true;
        const std::string txt = eclipseDeckToString(ce);
        const Case ci = eclipseDeckFromString(txt, "");
        const bool rt = ci.comp.sepTrains.size() == 1 && ci.comp.sepTrains[0].stages.size() == 1 && !ci.comp.wellSep.empty() &&
                        std::fabs(ci.comp.sepTrains[0].stages[0].second - 101325.0) < 1 && txt.find("WSEPCOND") != std::string::npos;
        std::printf("  secondary limit: FGPR %.5g -> %.5g (limit %.5g), FOPR %.5g -> %.5g; separator train LP: FGPR %.5g (field separator %.5g); FOPR %.5g; SEPCOND round trip %d\n",
                    ok ? o0[1].fgpr : 0.0, ok ? o1[1].fgpr : 0.0, glim, ok ? o0[1].fopr : 0.0, ok ? o1[1].fopr : 0.0, ok ? o2[1].fgpr : 0.0, ok ? o0[1].fgpr : 0.0, ok ? o2[1].fopr : 0.0, int(rt));
        check(ok && o1[1].fgpr <= glim * 1.0001 && o1[1].fopr < o0[1].fopr - 1 && std::fabs(o2[1].fgpr / o0[1].fgpr - 1) > 1e-3 && std::fabs(o2[1].fopr - o0[1].fopr) < 1e-3 * o0[1].fopr && rt,
              "producer secondary rate limits (WCONPROD) and SEPCOND / WSEPCOND separator trains, two-way deck");
    }
    {   // group control in the compositional core: FIELD ORAT target shared by two GRUP producers, GCONINJE gas RATE
        Case c = makeExample("comp");
        WellInput p2 = c.wells[0]; p2.name = "PROD2";
        for (auto& cp : p2.completions) { cp.i = 6; cp.j = 0; }
        c.wells.push_back(p2);
        c.schedule[0].controls.clear();
        for (const char* nm : {"PROD", "PROD2"}) { WellControl wc; wc.well = nm; wc.mode = Control::BHP; wc.bhp = 1000; wc.deckMode = "GRUP"; c.schedule[0].controls.push_back(wc); }
        GroupProdControl gp; gp.group = "FIELD"; gp.mode = Control::ORAT; gp.orat = 4000; c.schedule[0].gprod = {gp};
        for (auto& wc : c.schedule[4].controls) if (wc.well == "INJ") wc.deckMode = "GRUP";
        GroupInjControl gi; gi.group = "FIELD"; gi.phase = 2; gi.mode = "RATE"; gi.rate = 5000; c.schedule[4].ginj = {gi};
        comp::DeckModel m; std::string err; std::vector<CompSummary> out;
        const bool ok = comp::deckModel(c, m, &err) && comp::runDeckModel(m, out, {}, &err, 6);
        bool hold = ok && out.size() == 7;
        double w1 = 0, w2 = 0;
        if (hold) {
            for (size_t k = 1; k < out.size(); ++k) hold = hold && std::fabs(out[k].fopr - 4000) < 0.005 * 4000;
            w1 = out[2].wells[0].oil; w2 = out[2].wells[2].oil;
            hold = hold && w1 > 100 && w2 > 100 && std::fabs(out[6].fgir - 5000) < 0.005 * 5000 && std::fabs(out[5].fgir - 5000) < 0.005 * 5000;
        }
        if (!ok) std::printf("  error: %s\n", err.c_str());
        std::printf("  group control: FIELD ORAT 4000 -> FOPR %.1f / %.1f (PROD %.1f + PROD2 %.1f stb/d); GCONINJE gas RATE 5000 -> FGIR %.1f Mscf/d\n",
                    ok ? out[1].fopr : 0.0, ok ? out.back().fopr : 0.0, w1, w2, ok ? out.back().fgir : 0.0);
        check(hold, "compositional group control: GCONPROD / GCONINJE targets shared by guide rates among GRUP wells");
    }
    {   // GUIDERAT: the guide-rate formula replaces the potential-based split; A = 0 -> equal guide rates; deck round trip
        Case c = makeExample("comp");
        WellInput p2 = c.wells[0]; p2.name = "PROD2";
        for (auto& cp : p2.completions) { cp.i = 3; cp.j = 0; cp.skin = 20; }   // damaged completion: lower potential
        c.wells.push_back(p2);
        c.schedule[0].controls.clear();
        for (const char* nm : {"PROD", "PROD2"}) { WellControl wc; wc.well = nm; wc.mode = Control::BHP; wc.bhp = 1000; wc.deckMode = "GRUP"; c.schedule[0].controls.push_back(wc); }
        GroupProdControl gp; gp.group = "FIELD"; gp.mode = Control::ORAT; gp.orat = 4000; c.schedule[0].gprod = {gp};
        auto run = [&](const Case& cc, double& a, double& b, double& f) {
            comp::DeckModel m; std::string err; std::vector<CompSummary> out;
            const bool ok = comp::deckModel(cc, m, &err) && comp::runDeckModel(m, out, {}, &err, 4);
            if (!ok || out.size() < 5) return false;
            a = out[3].wells[0].oil; b = out[3].wells[2].oil; f = out[3].fopr; return true;
        };
        double a0 = 0, b0 = 0, f0 = 0, a1 = 0, b1 = 0, f1 = 0, a2 = 0, b2 = 0, f2 = 0;
        const bool r0 = run(c, a0, b0, f0);
        Case c1 = c; c1.schedule[0].guiderat.set = true; c1.schedule[0].guiderat.phase = "OIL"; c1.schedule[0].guiderat.A = 0; c1.schedule[0].guiderat.B = 1;
        const bool r1 = run(c1, a1, b1, f1);
        // penalise the gas-oil ratio strongly: GR = POT / (1 + 100 GOR^2), GOR in Mscf/stb
        Case c2 = c; c2.schedule[0].guiderat = c1.schedule[0].guiderat; c2.schedule[0].guiderat.A = 1; c2.schedule[0].guiderat.E = 100; c2.schedule[0].guiderat.F = 2;
        const bool r2 = run(c2, a2, b2, f2);
        const std::string path = "/tmp/ressim_guiderat.DATA";
        exportEclipseDeck(c2, path);
        const Case b = importEclipseDeck(path);
        const GuideRateSpec& g = b.schedule[0].guiderat;
        double a3 = 0, b3 = 0, f3 = 0;
        const bool rt = g.set && g.phase == "OIL" && g.A == 1 && g.B == 1 && g.E == 100 && g.F == 2 && run(b, a3, b3, f3) && std::fabs(a3 - a2) < 0.005 * 4000;
        std::printf("  GUIDERAT: potentials PROD %.1f / PROD2 %.1f; A=0 (equal) %.1f / %.1f; POT/(1+100 GOR^2) %.1f / %.1f (FOPR %.1f); round trip %s\n",
                    a0, b0, a1, b1, a2, b2, f2, rt ? "OK" : "FAILED");
        check(r0 && r1 && r2 && std::fabs(a1 - b1) < 0.01 * 4000 && std::fabs(a0 - b0) > 100 && std::fabs(f1 - 4000) < 20 && std::fabs(f2 - 4000) < 20 && rt,
              "GUIDERAT formula in the compositional group allocation (equal guide rates for A = 0, target held) + deck round trip");
        std::remove(path.c_str());
    }
    {   // deck round trip of all trapping keywords; full deck run with every mechanism on
        Case c = makeExample("co2");
        co2::TrappingOptions& T = c.comp.trap;
        T.diffusion = true; T.dCo2Water = 3e-9; T.saltPrecip = true; T.saltSolubility = 5.5; T.permFact = {{0.5, 0.1}, {0.8, 0.5}, {1.0, 1.0}};
        T.minerals = {co2::mineralPresets().front()}; T.sgtMax = 0.27; T.structTol = 0.05;
        c.props.hyst.on = true; c.props.hyst.model = 2;
        c.comp.structured = true;
        const std::string txt = eclipseDeckToString(c);
        const Case ci = eclipseDeckFromString(txt, "");
        const co2::TrappingOptions& U = ci.comp.trap;
        const bool rt = U.diffusion && std::fabs(U.dCo2Water - 3e-9) < 1e-15 && U.saltPrecip && std::fabs(U.saltSolubility - 5.5) < 1e-9 && U.permFact.size() == 3 &&
                        U.minerals.size() == 1 && std::fabs(U.minerals[0].logk25 - T.minerals[0].logk25) < 1e-9 && U.minerals[0].name == T.minerals[0].name &&
                        std::fabs(U.sgtMax - 0.27) < 1e-12 && std::fabs(U.structTol - 0.05) < 1e-12 && ci.props.hyst.on && ci.props.hyst.model == 2 &&
                        txt.find("PRECSALT") != std::string::npos && txt.find("DIFFCWAT") != std::string::npos && txt.find("PERMFACT") != std::string::npos;
        comp::DeckModel m; std::string err; std::vector<CompSummary> out;
        const bool run = comp::deckModel(ci, m, &err) && comp::runDeckModel(m, out, {}, &err, 3);
        const CompSummary& e = run ? out.back() : CompSummary();
        const double injKg = e.fgit * co2::MW_CO2 / co2::gasSurfaceVolume(c.comp.fluid.stdT, c.comp.fluid.stdP);
        const double tot = e.co2Mob + e.co2Immob + e.co2Diss + e.co2Min;
        std::printf("  all mechanisms (deck round trip %d, %s): day %.0f mobile %.4g t (structural %.3g), residual %.4g t, dissolved %.4g t, mineral %.4g t, NaCl solid %.4g t; in place / injected %.6f\n",
                    int(rt), m.base.hyst.on ? "Killough" : "no hysteresis", e.time, (e.co2Mob - e.co2Struct) / 1e3, e.co2Struct / 1e3, e.co2Immob / 1e3, e.co2Diss / 1e3, e.co2Min / 1e3, e.saltSolid / 1e3, tot / injKg);
        if (!run) std::printf("  error: %s\n", err.c_str());
        check(rt && run && m.base.hyst.on && m.base.hyst.model == 1 && std::fabs(tot / injKg - 1) < 0.02 && e.co2Min > 0,
              "trapping keywords two-way (DIFFUSE / DIFFCWAT / PRECSALT / SALTSOL / PERMFACT / EHYSTR / RESSIM MINERAL...), all mechanisms in one run");
    }
}

// ---- C6: CO2 as an ordinary EOS component that dissolves in the brine (E300 CO2SOL), speed-path regressions
static void testCo2Sol() {
    using namespace ressim; using namespace ressim::comp;
    std::printf("CO2SOL (CO2 as an EOS component dissolving in brine) and the fast flash path\n");
    const eos::CompFluid F = spe5Co2CompFluid(1.0);
    const eos::Model M = F.eosRegions[0];
    const double T = 323.15;
    {   // 1. equilibrium: brine under a CO2-rich gas cap -> x_aq = (f_CO2 / f_pure) x_sat (Spycher-Pruess, 1 m NaCl)
        CompCase c; c.nx = c.ny = c.nz = 1; c.dx = {10}; c.dy = {10}; c.dz = {10}; c.poro = {0.2}; c.permx = c.permy = c.permz = {1e-13};
        c.fluid = M; c.T = T; c.water = true; c.swc = 0.2; c.gravity = false;
        c.p0 = {150e5}; c.z0 = {{0.005, 0, 0, 0, 0, 0, 0.995}}; c.sw0 = {0.5};
        c.co2sol = true; c.solSalinity = 1.0;
        CompSimulator sim; sim.setup(c);
        const double xa = sim.dissolvedCO2()[0];
        double xs = 0, ys = 0; co2::solubility(co2::makeFluid(T, 1.0), 150e5, xs, ys);
        double hc = 0, aq = 0; sim.co2SolSplit(hc, aq);
        std::printf("  150 bar, 50 C, 1 m NaCl: x_aq %.5f vs Spycher-Pruess x_sat %.5f (gas 99.5 %% CO2); CO2 in brine %.4g kg, in gas %.4g kg\n", xa, xs, aq, hc);
        check(xa > 0.97 * xs && xa < 1.0001 * xs && aq > 0, "CO2SOL equilibrium: aqueous CO2 from fugacity equality matches the Spycher-Pruess solubility");
    }
    {   // 2. analytic Jacobian (implicit-function derivatives of the aqueous split) vs finite differences
        CompCase c; c.fluid = M; c.T = T; c.water = true; c.swc = 0.2; c.cr = 5e-10; c.pref = 2e7; c.gravity = true;
        for (int k = 0; k < 3; ++k) { c.cellPV.push_back(800.0 * (1 + 0.2 * k)); c.cellDepth.push_back(2500.0 + 4 * k); }
        c.conns = {{0, 1, 1.5e-11}, {1, 2, 1.2e-11}};
        c.sat.sw = {0.2, 0.5, 0.8, 1.0}; c.sat.krw = {0, 0.1, 0.5, 1}; c.sat.pcow = {2e5, 5e4, 1e4, 0};
        c.sat.sg = {0, 0.05, 0.4, 0.8}; c.sat.krg = {0, 0, 0.3, 0.9}; c.sat.pcog = {0, 1e4, 3e4, 6e4};
        c.sat.so = {0, 0.15, 0.5, 0.8}; c.sat.krow = {0, 0, 0.35, 0.9}; c.sat.krog = {0, 0, 0.3, 0.9};
        c.p0 = {130e5, 125e5, 120e5};
        c.z0 = {{0.30, 0.03, 0.07, 0.15, 0.10, 0.05, 0.30}, {0.45, 0.03, 0.07, 0.20, 0.15, 0.05, 0.05}, {0.40, 0.03, 0.07, 0.20, 0.15, 0.05, 0.10}};
        c.sw0 = {0.4, 0.5, 0.6};
        c.co2sol = true; c.solSalinity = 0.5;
        CompWell prd; prd.cells = {2}; prd.wi = {1e-12}; prd.bhp = 100e5; prd.control = 0;
        c.wells = {prd};
        CompSimulator sim; sim.setup(c);
        const int ncw = sim.ncw(), nv = ncw + 1;
        std::vector<double> p = sim.pressure(), N = sim.moles();
        sim.storeOld();
        for (size_t k = 0; k < N.size(); ++k) N[k] *= 1.0 + 0.004 * double(k % 3);
        sim.setState(p, N);
        const double dt = 86400.0;
        std::vector<double> R0, J;
        sim.residualAndJacobian(dt, R0, &J);
        const int m = int(R0.size());
        double err = 0, scale = 0;
        for (int col = 0; col < m; ++col) {
            const int q = col / nv, k = col % nv;
            double Nt = 0; for (int i = 0; i < ncw; ++i) Nt += N[size_t(q) * ncw + i];
            const double h = k == 0 ? p[size_t(q)] * 1e-7 : Nt * 1e-7;
            std::vector<double> pp = p, NN = N, pm = p, Nm = N;
            if (k == 0) { pp[size_t(q)] += h; pm[size_t(q)] -= h; } else { NN[size_t(q) * ncw + k - 1] += h; Nm[size_t(q) * ncw + k - 1] -= h; }
            std::vector<double> Rp, Rm;
            sim.setState(pp, NN); sim.residualAndJacobian(dt, Rp, nullptr);
            sim.setState(pm, Nm); sim.residualAndJacobian(dt, Rm, nullptr);
            for (int r = 0; r < m; ++r) {
                const double fd = (Rp[size_t(r)] - Rm[size_t(r)]) / (2 * h), an = J[size_t(r) * m + col];
                err = std::max(err, std::fabs(fd - an) * (k == 0 ? p[size_t(q)] : Nt));
                scale = std::max(scale, std::fabs(an) * (k == 0 ? p[size_t(q)] : Nt));
            }
        }
        sim.setState(p, N);
        std::printf("  3 cells, 7 components + water, CO2 in oil / gas / brine, BHP producer: Jacobian max scaled error %.2e\n", err / scale);
        check(err / scale < 1e-5, "CO2SOL: aqueous CO2 split derivatives (implicit function theorem) - Jacobian matches finite differences");
    }
    {   // 3. CO2 EOR example (WAG, CO2SOL): CO2 material balance through injection, dissolution and production
        Case cs = makeExample("co2eor");
        DeckModel dm; std::string err;
        const bool ok = deckModel(cs, dm, &err);
        std::vector<CompSummary> out; CompStats st;
        double inPlace0 = -1, inPlace = 0, inj = 0, prod = 0, diss = 0; int ic = -1;
        const bool ran = ok && runDeckModel(dm, out, [&](const CompSummary& S, const CompSimulator& sim) {
            ic = dm.base.co2Comp >= 0 ? dm.base.co2Comp : 6;
            if (inPlace0 < 0) inPlace0 = sim.componentInPlace(ic);
            inPlace = sim.componentInPlace(ic); inj = sim.injected(ic); prod = sim.produced(ic); diss = S.co2Diss;
            return true;
        }, &err, 12, &st);
        const double bal = inj > 0 ? (inPlace - inPlace0 - (inj - prod)) / inj : 1.0;
        std::printf("  co2eor 3 y (1 y CO2 + 1 y water): injected %.4g kmol CO2, produced %.4g kmol, dissolved in brine %.4g t; balance error %.2e; %d steps, %d Newton\n",
                    inj / 1e3, prod / 1e3, diss / 1e3, bal, st.steps, st.newton);
        check(ran && std::fabs(bal) < 1e-6 && diss > 0, "CO2 EOR WAG example (CO2SOL): CO2 material balance and dissolution in the brine");
    }
    {   // 4. fast flash path (warm-started Newton) reproduces the full successive-substitution + stability flash
        const eos::TParams tp = eos::atTemperature(M, 344.26);
        std::vector<double> z = {0.5, 0.03, 0.07, 0.2, 0.15, 0.05, 0.0};
        double maxd = 0; int two = 0;
        for (double P = 80e5; P <= 160e5; P += 10e5) {
            const eos::FlashResult a = eos::flash(M, tp, P, z);
            if (a.phases != 2) continue;
            ++two;
            std::vector<double> K = a.K; for (double& k : K) k *= 1.03;
            eos::FlashOptions fo; fo.Kinit = &K;
            const eos::FlashResult b = eos::flash(M, tp, P, z, fo);
            maxd = std::max(maxd, std::fabs(a.beta - b.beta));
            for (size_t i = 0; i < z.size(); ++i) maxd = std::max(maxd, std::fabs(a.x[i] - b.x[i]) + std::fabs(a.y[i] - b.y[i]));
        }
        std::printf("  warm-started flash vs cold flash: %d two-phase pressures, max difference %.2e\n", two, maxd);
        check(two >= 3 && maxd < 1e-8, "fast flash path (direct Newton from previous K) agrees with the full flash");
    }
    {   // 5. absent component (CO2 not yet arrived): Newton converges and the trace-partition sensitivities match FD
        const eos::TParams tp = eos::atTemperature(M, 344.26);
        std::vector<double> z = {0.5, 0.03, 0.07, 0.2, 0.15, 0.05, 0.0};
        const double P = 120e5;
        const eos::FlashResult f = eos::flash(M, tp, P, z);
        eos::FlashSens S; const bool ok = eos::flashSensitivities(M, tp, P, z, f, S, false);
        const size_t nc = z.size(), j = 6;
        const double h = 1e-7;
        std::vector<double> zp = z; zp[j] += h;
        std::vector<double> K = f.K; eos::FlashOptions fo; fo.Kinit = &K;
        const eos::FlashResult fp = eos::flash(M, tp, P, zp, fo);
        double err = std::fabs((fp.beta - f.beta) / h - S.dbeta_dz[j]), sc = std::fabs(S.dbeta_dz[j]);
        for (size_t i = 0; i < nc; ++i) {
            err = std::max(err, std::fabs((fp.x[i] - f.x[i]) / h - S.dx_dz[i * nc + j]));
            err = std::max(err, std::fabs((fp.y[i] - f.y[i]) / h - S.dy_dz[i * nc + j]));
            sc = std::max(sc, std::max(std::fabs(S.dx_dz[i * nc + j]), std::fabs(S.dy_dz[i * nc + j])));
        }
        const double Kco2 = f.K[j];
        std::printf("  absent CO2: flash converged %d (%d it), K_CO2 %.4f, d/dz_CO2 FD error %.2e (scale %.2e)\n", f.converged ? 1 : 0, f.iterations, Kco2, err, sc);
        check(f.converged && ok && std::isfinite(err) && err < 1e-4 * std::max(1.0, sc) && Kco2 > 0.5 && Kco2 < 10,
              "absent component: converged flash, phi-ratio K and trace-partition sensitivities agree with FD");
    }
}

static void testEorTools() {
    std::printf("CO2 EOR tools: swelling test / FCMP, mixing-cell MMP, WAG generator, CO2SOL deck round trip\n");
    const eos::CompFluid F = spe5Co2CompFluid(1.0);
    const eos::Model& M = F.eosRegions[0];
    const double T = eos::fahrenheitToK(160);
    std::vector<double> oil = F.zi, co2(7, 0.0); co2[6] = 1.0;
    {   // 1. swelling test
        std::vector<double> al; for (int i = 0; i <= 19; ++i) al.push_back(0.05 * i);
        const eor::SwellingResult S = eor::swellingTest(M, T, oil, co2, al);
        bool mono = true; double last = 0; int found = 0;
        for (auto& r : S.rows) { if (r.psat > 0) ++found; if (r.psat > 0 && !r.dew) { mono = mono && r.swelling >= last; last = r.swelling; } }
        std::printf("  swelling: Psat(oil) %.2f psia, FCMP %.1f psia at %.1f mol%% CO2, SF(50 mol%%) %.4f, %d / %zu saturation points\n",
                    eos::paToPsia(S.psatOil), eos::paToPsia(S.fcmp), 100 * S.fcmpAlpha, S.rows[10].swelling, found, S.rows.size());
        check(S.ok && std::fabs(eos::paToPsia(S.psatOil) - 2302.89) < 0.5 && mono && found == int(S.rows.size()) && S.fcmp > S.psatOil && S.rows[10].swelling > 1.2,
              "swelling test: bubble point of SPE5 oil, monotonic swelling factor, complete P-X diagram, FCMP above Pb");
        // 2. MMP: CO2 multi-contact MMP below the FCMP, close to the Cronquist correlation; SPE5 solvent MMP
        eor::MmpOptions o; o.pLo = eos::psiaToPa(1000); o.pHi = eos::psiaToPa(6000);
        const eor::MmpResult R = eor::mixingCellMmp(M, T, oil, co2, o);
        const double cq = eor::mmpCronquist(T, 169.0, 0.5), ym = eor::mmpYelligMetcalfe(T);
        std::vector<double> gas = {0.77, 0.20, 0.03, 0, 0, 0, 0};
        const eor::MmpResult R2 = eor::mixingCellMmp(M, T, oil, gas, o);
        std::printf("  MMP (160 F): CO2 %.0f psia (Cronquist %.0f, Yellig-Metcalfe %.0f), SPE5 solvent %.0f psia; %d + %d flashes\n",
                    eos::paToPsia(R.mmp), eos::paToPsia(cq), eos::paToPsia(ym), eos::paToPsia(R2.mmp), R.flashes, R2.flashes);
        check(R.ok && R.mmp < S.fcmp && R.mmp > S.psatOil && std::fabs(R.mmp / cq - 1) < 0.15 && R2.ok && R2.mmp > R.mmp && eos::paToPsia(R2.mmp) < 4000,
              "mixing-cell MMP: CO2 MMP between Pb and FCMP, within 15 % of Cronquist; lean hydrocarbon solvent needs a higher MMP");
        // slim-tube simulation (1-D EOS displacements, 1.2 HCPV): the break-over pressure agrees with the mixing-cell MMP
        std::vector<double> ps; for (double pp = 2000; pp <= 4000; pp += 400) ps.push_back(eos::psiaToPa(pp));
        const eor::SlimTubeResult ST = eor::slimTubeMmp(M, T, oil, co2, ps);
        std::printf("  slim tube: RF");
        for (size_t k = 0; k < ST.rf.size(); ++k) std::printf(" %.0f:%.3f", eos::paToPsia(ST.p[k]), ST.rf[k]);
        std::printf("; break-over MMP %.0f psia (mixing cell %.0f)\n", eos::paToPsia(ST.mmpBreak), eos::paToPsia(R.mmp));
        bool rising = true; for (size_t k = 1; k < ST.rf.size(); ++k) rising = rising && ST.rf[k] > ST.rf[k - 1] - 0.005;
        check(ST.ok && rising && ST.rf.back() > 0.97 && std::fabs(ST.mmpBreak / R.mmp - 1) < 0.10,
              "slim-tube MMP (1-D compositional displacements) within 10 % of the mixing-cell MMP; miscible recovery > 97 %");
    }
    {   // 3. WAG generator on the co2eor example + deck round trip of CO2SOL / streams / WAG
        Case c = makeExample("co2eor");
        const double total0 = [&] { double t = 0; for (auto& r : c.schedule) t += r.length; return t; }();
        eor::WagSpec w; w.well = "INJ"; w.stream = "CO2"; w.start = 365.25; w.gasDays = 120; w.waterDays = 60; w.cycles = 5;
        w.gasRate = 6000; w.waterRate = 5000; w.bhpMax = 6000; w.taper = 0.9; w.chaseWaterDays = 400;
        std::string msg;
        const int n = eor::applyWagSchedule(c, w, &msg);
        double t = 0; int gasSteps = 0, waterSteps = 0, streams = 0; bool aligned = true;
        double gl = 120, expect = 365.25; std::vector<double> marks;
        for (int k = 0; k < 5; ++k) { marks.push_back(expect); expect += gl; marks.push_back(expect); expect += 60; gl *= 0.9; }
        marks.push_back(expect);
        size_t mi = 0;
        for (auto& r : c.schedule) {
            for (auto& wc : r.controls) if (wc.well == "INJ") {
                (wc.phase == 2 ? gasSteps : waterSteps)++;
                aligned = aligned && mi < marks.size() && std::fabs(t - marks[mi]) < 1e-6; ++mi;
            }
            for (auto& q : r.injStreams) if (q.first == "INJ" && q.second == "CO2") ++streams;
            t += r.length;
        }
        std::printf("  %s; gas %d / water %d controls, %d WINJGAS, schedule %.1f d (was %.1f)\n", msg.c_str(), gasSteps, waterSteps, streams, t, total0);
        check(n == 11 && gasSteps == 5 && waterSteps == 6 && streams == 5 && aligned && std::fabs(t - total0) < 1e-6,
              "WAG generator: tapered 5-cycle WAG + chase, controls on the half-cycle boundaries, WINJGAS CO2 in the gas half-cycles");
        const std::string path = "/tmp/ressim_co2sol_rt.DATA";
        exportEclipseDeck(c, path);
        const Case d = importEclipseDeck(path);
        int g2 = 0, s2 = 0; double t2 = 0;
        for (auto& r : d.schedule) { for (auto& wc : r.controls) if (wc.well == "INJ" && wc.phase == 2) ++g2; for (auto& q : r.injStreams) if (q.second == "CO2") ++s2; t2 += r.length; }
        const bool fl = d.comp.enabled && d.comp.fluid.co2sol && std::fabs(d.comp.fluid.salinity - 1.0) < 1e-6 && d.comp.fluid.eosRegions.size() == 1 && d.comp.fluid.eosRegions[0].nc() == 7;
        std::printf("  deck round trip: CO2SOL %d, salinity %.3f, nc %d, gas controls %d, WINJGAS CO2 %d, %.1f d\n", d.comp.fluid.co2sol ? 1 : 0, d.comp.fluid.salinity,
                    d.comp.fluid.eosRegions.empty() ? 0 : d.comp.fluid.eosRegions[0].nc(), g2, s2, t2);
        check(fl && g2 == 5 && s2 >= 5 && std::fabs(t2 - t) < 1e-3, "CO2 EOR deck (CO2SOL + SALINITY + WELLSTRE / WINJGAS WAG) exports and re-imports");
        std::remove(path.c_str());
    }
    {   // 4. SOLUBILI table (E300 CO2SOL): half of the Spycher-Pruess solubility -> about half the dissolved CO2; deck round trip
        using namespace comp;
        auto runDiss = [&](const Case& cs) {
            DeckModel dm; std::string err; double diss = 0;
            if (!deckModel(cs, dm, &err)) return -1.0;
            std::vector<CompSummary> out; CompStats st;
            runDeckModel(dm, out, [&](const CompSummary& S, const CompSimulator&) { diss = S.co2Diss; return true; }, &err, 8, &st);
            return diss;
        };
        Case c0 = makeExample("co2eor");
        Case c1 = c0;
        const double T = c1.comp.fluid.resT;
        const co2::Fluid cf = co2::makeFluid(T, c1.comp.fluid.salinity);
        const double molGas = c1.comp.fluid.stdP / (eos::Rgas * c1.comp.fluid.stdT), molW = 1000.0 / 0.018015;
        for (double pb : {1.01325, 20.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0}) {
            double x = 0, y = 0; co2::solubility(cf, pb * 1e5, x, y);
            eos::CompFluid::SolRow r; r.p = pb * 1e5; r.rsw = 0.5 * x / (1 - x) * molW / molGas; r.bw = 1.02; r.cw = 4e-10; r.muw = 5.5e-4;
            c1.comp.fluid.solubili.push_back(r);
        }
        const double d0 = runDiss(c0), d1 = runDiss(c1);
        const std::string path = "/tmp/ressim_solubili_rt.DATA";
        exportEclipseDeck(c1, path);
        const Case b = importEclipseDeck(path);
        const auto& S1 = c1.comp.fluid.solubili; const auto& S2 = b.comp.fluid.solubili;
        bool same = S1.size() == S2.size();
        for (size_t i = 0; same && i < S1.size(); ++i) same = std::fabs(S2[i].p / S1[i].p - 1) < 1e-9 && std::fabs(S2[i].rsw - S1[i].rsw) < 1e-9 * (1 + S1[i].rsw) && std::fabs(S2[i].muw / S1[i].muw - 1) < 1e-9;
        std::printf("  SOLUBILI: dissolved CO2 after 2 y %.4g t (default) vs %.4g t (table at 0.5 x solubility), ratio %.3f; round trip %d rows %s\n",
                    d0 / 1e3, d1 / 1e3, d0 > 0 ? d1 / d0 : 0.0, int(S2.size()), same ? "identical" : "DIFFERENT");
        check(d0 > 0 && d1 > 0 && d1 / d0 > 0.35 && d1 / d0 < 0.65 && same, "CO2SOL SOLUBILI table drives the brine solubility and survives export / import");
        std::remove(path.c_str());
    }
}

// ------------------------------------------------------------------------------------------------ thermal core (C7)
static ressim::thermal::Component thermHeavy(double mu300) {
    ressim::thermal::Component q; q.name = "HEAVY"; q.mw = 0.5; q.rhoRef = 980; q.cpL = 7e-10; q.ct1 = 6e-4;
    q.bvisc = 5180; q.avisc = mu300 / std::exp(q.bvisc / 300.0); q.cpl1 = 1000; q.dead = true;
    return q;
}
static ressim::thermal::Case thermLine(int n, bool steam, double days) {
    (void)days;
    ressim::thermal::Case c; c.nx = n; c.dx = {5}; c.dy = {10}; c.dz = {5}; c.poro = {0.3}; c.permx = {5e-12}; c.gravity = false; c.cr = 5e-10;
    c.comps = {thermHeavy(1.0)}; c.p0 = {3e6}; c.T0 = {300}; c.sw0 = {0.25}; c.x0 = {{1.0}};
    ressim::thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {5e-12}; inj.zinj = {0, 1}; inj.bhpMax = 6e6;
    inj.rate = 10.0 / 86400 * 1000 / ressim::steam::kMw;   // 10 m3/d cold-water equivalent
    if (steam) { inj.quality = 0.7; inj.tinj = 565; } else inj.tinj = 300;
    ressim::thermal::Well pr; pr.name = "PROD"; pr.cells = {n - 1}; pr.wi = {5e-12}; pr.bhp = 2e6;
    c.wells = {inj, pr};
    return c;
}
static void testThermal() {
    std::printf("Thermal core (steam tables, PU flash, energy balance, heat loss, steamflood)\n");
    {
        const ressim::steam::Props a = ressim::steam::liquid(300, 3e6), b = ressim::steam::vapour(300, 3500);
        const bool ok = std::fabs(a.v - 0.00100215168) < 1e-12 && std::fabs(a.h - 115331.273) < 1e-2 && std::fabs(a.cp - 4173.01218) < 1e-3 &&
                        std::fabs(b.v - 39.4913866) < 1e-6 && std::fabs(b.h - 2549911.45) < 1 &&
                        std::fabs(ressim::steam::tsat(1e5) - 372.755919) < 1e-5 && std::fabs(ressim::steam::tsat(1e7) - 584.149488) < 1e-5 &&
                        std::fabs(ressim::steam::viscosity(298.15, 998) - 8.897351e-4) < 1e-9;
        std::printf("  IF97: v1 %.9g h1 %.6f v2 %.7g tsat(0.1 MPa) %.6f mu %.6e\n", a.v, a.h, b.v, ressim::steam::tsat(1e5), ressim::steam::viscosity(298.15, 998));
        check(ok, "IAPWS-IF97 regions 1/2/4 and IAPWS-2008 viscosity match the reference values");
    }
    {   // univariant pure water: steam injected into a small water-filled cell produced at 1 MPa -> T = Tsat(p), two-phase
        ressim::thermal::Case c; c.nx = 1; c.dx = {1}; c.dy = {1}; c.dz = {1}; c.poro = {0.3}; c.permx = {1e-12}; c.gravity = false;
        c.p0 = {1e6}; c.T0 = {420}; c.sw0 = {1.0};
        ressim::thermal::Well inj; inj.name = "I"; inj.injector = true; inj.cells = {0}; inj.wi = {1e-11}; inj.zinj = {1}; inj.quality = 1.0; inj.pinj = 1.5e6;
        inj.rate = 0.05; inj.bhpMax = 5e6;
        ressim::thermal::Well pr; pr.name = "P"; pr.cells = {0}; pr.wi = {1e-11}; pr.bhp = 1e6;
        c.wells = {inj, pr};
        ressim::thermal::Simulator s; std::string err;
        bool ok = s.setup(c, &err) && s.advanceTo(5 * 86400.0, 0.25 * 86400.0);
        const double Ts = ressim::steam::tsat(s.pressure()[0]);
        std::printf("  pure water: p %.4f MPa T %.4f K Tsat %.4f K Sg %.3f\n", s.pressure()[0] / 1e6, s.temperature()[0], Ts, s.sg()[0]);
        check(ok && std::fabs(s.temperature()[0] - Ts) < 1e-6 && s.sg()[0] > 0.05 && s.sg()[0] < 0.999, "pure-water (p, U) flash on the saturation line: T = Tsat(p), two-phase");
    }
    {   // closed two-cell conduction box
        ressim::thermal::Case c; c.nx = 2; c.dx = {10}; c.dy = {10}; c.dz = {10}; c.poro = {0.3}; c.permx = {1e-13}; c.gravity = false; c.cr = 5e-10;
        c.comps = {thermHeavy(5.0)}; c.p0 = {5e6}; c.T0 = {300, 400}; c.sw0 = {0.4}; c.x0 = {{1.0}};
        ressim::thermal::Simulator s; std::string err;
        bool ok = s.setup(c, &err);
        const double E0 = s.energyInPlace(), Nw = s.fluidMoles(1), No = s.fluidMoles(0);
        ok = ok && s.advanceTo(3e4 * 86400.0, 300 * 86400.0);
        const double dE = std::fabs(s.energyInPlace() - E0) / std::fabs(E0);
        std::printf("  conduction box: T %.4f / %.4f K, |dE|/E %.2e, mass drift %.2e\n", s.temperature()[0], s.temperature()[1], dE,
                    std::fabs(s.fluidMoles(1) - Nw) / Nw + std::fabs(s.fluidMoles(0) - No) / No);
        check(ok && dE < 1e-10 && std::fabs(s.temperature()[0] - s.temperature()[1]) < 0.05 && s.temperature()[0] > 300 && s.temperature()[0] < 400,
              "closed conduction box conserves energy exactly and equilibrates");
    }
    {   // Vinsome-Westerveld heat loss vs the analytic conduction into a half space (constant face temperature)
        ressim::thermal::Case c; c.nx = 1; c.dx = {100}; c.dy = {100}; c.dz = {10}; c.poro = {0.3}; c.permx = {1e-13}; c.gravity = false; c.cr = 5e-10;
        ressim::thermal::Component q; q.name = "H"; q.mw = 0.5; q.rhoRef = 980; q.cpl1 = 1000; q.dead = true; c.comps = {q};
        c.p0 = {5e6}; c.T0 = {450}; c.sw0 = {0.5}; c.x0 = {{1.0}}; c.rockCv = 1e13; c.heatLoss = true; c.obT = 300;
        ressim::thermal::Simulator s; std::string err;
        bool ok = s.setup(c, &err) && s.advanceTo(3650 * 86400.0, 30 * 86400.0);
        const double t = 3650 * 86400.0, kap = c.obK / c.obCv;
        const double an = 2 * (2 * c.obK * 1e4 * 150 * std::sqrt(t / (M_PI * kap)));
        std::printf("  heat loss after 10 y: %.4e J vs analytic %.4e J (ratio %.4f)\n", s.cumHeatLoss(), an, s.cumHeatLoss() / an);
        check(ok && std::fabs(s.cumHeatLoss() / an - 1) < 0.05, "Vinsome-Westerveld heat loss within 5% of the analytic half-space solution");
    }
    {   // 1D hot water / steam vs cold water flood of a 1 Pa s heavy oil
        double rf[2] = {0, 0}, ebal[2] = {1, 1}, tmax[2] = {0, 0}, dts[2] = {0, 0};
        bool ok = true, steamOk = false;
        for (int m = 0; m < 2; ++m) {
            ressim::thermal::Case c = thermLine(30, m == 1, 400);
            ressim::thermal::Simulator s; std::string err;
            if (!s.setup(c, &err)) { ok = false; break; }
            const double E0 = s.energyInPlace(), st = s.stockTankOilInPlace();
            if (!s.advanceTo(1000 * 86400.0, 20 * 86400.0)) { ok = false; break; }
            rf[m] = s.wellReports()[1].cumO / st;
            double hScale = 1.0; for (const auto& r : s.wellReports()) hScale += std::fabs(r.cumH);
            ebal[m] = std::fabs((s.energyInPlace() - E0) - (s.cumEnthalpyIn() - s.cumHeatLoss())) / hScale;
            for (double T : s.temperature()) tmax[m] = std::max(tmax[m], T);
            if (m == 1) {
                int nst = 0; bool atTs = true;
                for (int q = 0; q < s.size(); ++q) if (s.sg()[size_t(q)] > 1e-3) { ++nst; atTs = atTs && std::fabs(s.temperature()[size_t(q)] - ressim::steam::tsat(s.pressure()[size_t(q)])) < 1e-3; }
                steamOk = nst > 0 && atTs;
                std::printf("  steam zone: %d cell(s), all at Tsat(p): %s\n", nst, atTs ? "yes" : "no");
            }
            dts[m] = s.stats().steps;
        }
        std::printf("  1D flood 1000 d: RF cold %.4f / steam %.4f, energy balance %.1e / %.1e, Tmax %.1f / %.1f K, steps %.0f / %.0f\n",
                    rf[0], rf[1], ebal[0], ebal[1], tmax[0], tmax[1], dts[0], dts[1]);
        check(ok && ebal[0] < 1e-5 && ebal[1] < 1e-5, "thermal flood energy balance (in place = injected - produced - lost) < 1e-5");
        check(ok && rf[1] > 1.3 * rf[0] && steamOk, "steamflood forms a saturated steam zone and out-recovers the cold waterflood");
    }
    {   // 5. E300 THERMAL examples: deck round trip of the thermal data + short runs on the thermal core
        using namespace ressim;
        int rtOk = 0, runOk = 0;
        for (const char* id : {"steamflood", "css", "sagd"}) {
            const Case c = makeExample(id);
            const Case b = eclipseDeckFromString(eclipseDeckToString(c));
            const ThermalInput& A = c.thermal; const ThermalInput& B = b.thermal;
            bool same = B.enabled && B.nc() == A.nc() && B.winj.size() == A.winj.size() && B.heatLoss == A.heatLoss && B.traps.size() == A.traps.size() &&
                        B.oilViscT.size() == A.oilViscT.size() && b.engineUnsupported.empty();
            for (int i = 0; same && i < A.nc(); ++i) {
                const ThermalComponent& x = A.comps[size_t(i)]; const ThermalComponent& y = B.comps[size_t(i)];
                same = x.name == y.name && x.dead == y.dead && std::fabs(x.mw - y.mw) < 1e-9 * x.mw && std::fabs(x.dref - y.dref) < 1e-9 * x.dref &&
                       std::fabs(x.kvcr[1] - y.kvcr[1]) <= 1e-9 * std::fabs(x.kvcr[1]) && std::fabs(x.kvcr[3] - y.kvcr[3]) <= 1e-9 * std::fabs(x.kvcr[3]);
            }
            for (size_t r = 0; same && r < A.oilViscT.size(); ++r)
                for (size_t k = 0; k < A.oilViscT[r].size(); ++k) same = same && std::fabs(A.oilViscT[r][k] - B.oilViscT[r][k]) <= 1e-9 * std::fabs(A.oilViscT[r][k]);
            if (!same) std::printf("  %s: thermal deck round trip mismatch\n", id);
            rtOk += same;
            // short run: first two report steps, energy balance of fluids + rock against net injected enthalpy and heat loss
            thermal::DeckModel m; std::string err;
            std::vector<thermal::Summary> out; thermal::Stats st;
            double e0 = 0, e1 = 0, hin = 0, loss = 0, gross = 0; bool first = true;
            const auto t0 = std::chrono::steady_clock::now();
            const bool ran = thermal::deckModel(c, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
                if (first) { e0 = s.energyInPlace(); first = false; }   // U may be negative (T below the reference)
                e1 = s.energyInPlace(); hin = s.cumEnthalpyIn(); loss = s.cumHeatLoss();
                gross = 0; for (const auto& r : s.wellReports()) gross += std::fabs(r.cumH);
                return true;
            }, &err, std::string(id) == "sagd" ? 1 : 2, &st);
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            const double ebal = std::fabs((e1 - e0) - (hin - loss)) / std::max(1.0, gross + std::fabs(loss));   // vs gross well enthalpy
            const thermal::Summary& L = out.empty() ? thermal::Summary() : out.back();
            std::printf("  %-10s %.0f d: FOPT %.4g  FWIT(CWE) %.4g  Tmax %.1f  steam PV %.4f  energy balance %.1e  %d steps %d Newton %.1f s %s\n",
                        id, L.time, L.fopt, L.fwit, L.fmaxT, L.fsteamPV, ebal, st.steps, st.newton, sec, ran ? "" : err.c_str());
            const bool good = ran && ebal < 1e-5 && L.fwit > 0 && L.fmaxT > 100;
            runOk += good;
        }
        check(rtOk == 3, "thermal examples: E300 THERMAL deck export -> import keeps components, K-values, OILVISCT, WINJTEMP, heat loss");
        check(runOk == 3, "thermal examples (steamflood, CSS, SAGD): first report steps run, steam injected (T > 100 C), energy balance < 1e-5");
    }
    {   // 6. steam trap (C9): the SAGD producer with a 75 C subcool target - the dynamic BHP floor Psat(T_hot + 75) becomes
        //    active, the produced fluid is exactly 75 C subcooled, and the implicit d bhp / d X coupling is exact
        using namespace ressim;
        Case c = makeExample("sagd");
        c.thermal.traps = {{"PROD", 0, 75.0}};
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        double sub = 0, jac = -1; bool trap = false;
        const bool ran = thermal::deckModel(c, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator& s) {
            if (S.time <= 0) return true;
            const auto& R = s.wellReports()[0];
            double Th = 0; for (int q : s.model().wells[0].cells) Th = std::max(Th, s.temperature()[size_t(q)]);
            sub = steam::tsat(R.bhp) - Th; trap = R.trapActive;
            jac = const_cast<thermal::Simulator&>(s).jacobianCheck(86400.0, nullptr, nullptr, s.model().wells[0].cells[0]);
            return true;
        }, &err, 1, &st);
        std::printf("  steam trap 75 C: active %d, produced-fluid subcool %.5f C, Jacobian (producer cell columns) %.1e, %d steps %s\n", int(trap), sub, jac, st.steps, ran ? "" : err.c_str());
        check(ran && trap && std::fabs(sub - 75) < 1e-3 && jac >= 0 && jac < 1e-2, "steam-trap (subcool) producer constraint: implicit BHP floor holds the subcool exactly, Jacobian verified");
    }
}

// ------------------------------------------------------------------------------------------------ C11: near-critical water,
// injector connections at the BHP limit, small analytic verifications (conduction, Theis drawdown)
static double expint1(double u) {   // E1(u), series (u < 1) / continued fraction
    if (u < 1) {
        double s = -0.57721566490153286 - std::log(u), t = 1;
        for (int k = 1; k < 60; ++k) { t *= -u / k; s -= t / k; if (std::fabs(t / k) < 1e-17) break; }
        return s;
    }
    double b = u + 1, c = 1e300, d = 1 / b, h = d;
    for (int i = 1; i < 200; ++i) { const double a = -double(i) * i; b += 2; d = 1 / (a * d + b); c = b + a / c; const double e = c * d; h *= e; if (std::fabs(e - 1) < 1e-16) break; }
    return h * std::exp(-u);
}
// C11 verification: Lauwerier (hot water in a layer, conductive heat loss into the cap rock)
static void testLauwerier() {
    using namespace ressim;
    std::printf("== C11 verification: Lauwerier (hot-water injection with conductive heat loss to the cap rock)\n");
    const int nx = 60, nzc = 16, nz = nzc + 1;
    // thin layer with a higher conductivity: Lauwerier assumes a vertically uniform layer temperature (no conduction
    // resistance inside the layer), which a finite cell with conductivity lamA approximates to (H / 2 / lamA) / (sqrt(pi kappa t) / lam)
    const double dx = 2.0, dy = 10.0, H = 0.4, lam = 2.0, lamA = 6.0;
    thermal::Case c; c.nx = nx; c.ny = 1; c.nz = nz; c.gravity = false; c.cr = 0;
    c.dx = {dx}; c.dy = {dy};
    const size_t n = size_t(nx) * nz;
    std::vector<double> dzl(size_t(nz), 0.0);
    double d = 0.2;
    for (int k = nzc - 1; k >= 0; --k) { dzl[size_t(k)] = d; d *= 1.35; }   // cap rock: fine next to the layer
    dzl[size_t(nzc)] = H;                                                    // the aquifer (bottom layer)
    c.dz.assign(n, 0.0);
    for (size_t q = 0; q < n; ++q) c.dz[q] = dzl[q / size_t(nx)];
    c.poro.assign(n, 0.01); c.permx.assign(n, 1e-20); c.permz.assign(n, 1e-20);
    for (int i = 0; i < nx; ++i) { const size_t q = size_t(i + nx * nzc); c.poro[q] = 0.25; c.permx[q] = 1e-12; c.permz[q] = 1e-12; }
    c.p0 = {5e6}; c.T0 = {300}; c.sw0 = {1.0}; c.thconsf = 0;
    c.thconr.assign(n, lam);
    for (int i = 0; i < nx; ++i) c.thconr[size_t(i + nx * nzc)] = lamA;
    const double rate = 3.7;   // mol/s of water
    thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {nx * nzc}; inj.wi = {1e-10}; inj.zinj = {1.0};
    inj.rate = rate; inj.tinj = 320; inj.quality = -1;
    thermal::Well prod; prod.name = "PROD"; prod.cells = {nx * nzc + nx - 1}; prod.wi = {1e-10}; prod.bhp = 5e6;
    c.wells = {inj, prod};
    thermal::Simulator s; std::string err;
    const double t = 365.0 * 86400.0;
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = s.setup(c, &err) && s.advanceTo(t, 1.0 * 86400.0);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!ok) std::printf("  error: %s (t = %.3f d, steps %d, cuts %d)\n", err.c_str(), s.time() / 86400, s.stats().steps, s.stats().cuts);
    // analytic (one-sided loss, layer thickness H): T_D = erfc(lam x dy / (2 mdot cp sqrt(kappa_r t'))), t' = t - x / v_T
    const steam::Props W = steam::liquid(310, 5e6);
    const double mcp = rate * steam::kMw * W.cp;                                 // W/K advected
    const double rcA = 0.25 * W.cp / W.v + 0.75 * c.rockCv;                      // layer
    const double rcR = 0.01 * W.cp / W.v + 0.99 * c.rockCv, kap = lam / rcR;     // cap rock
    double worst = 0, front = 0;
    for (int i = 0; ok && i < nx - 10; ++i) {
        const double x = (i + 0.5) * dx;
        const double tau = rcA * H * dy * x / mcp;
        double TD = 0;
        if (t > tau) TD = std::erfc(lam * x * dy / (2 * mcp * std::sqrt(kap * (t - tau))));
        const double Ts = (s.temperature()[size_t(i + nx * nzc)] - 300.0) / 20.0;
        worst = std::max(worst, std::fabs(Ts - TD));
        if (TD > 0.5) front = x;
        if (std::getenv("LAUW_DEBUG") && i % 4 == 0) std::printf("    x %6.1f  sim %.4f  Lauwerier %.4f\n", x, Ts, TD);
    }
    std::printf("  1 year, %d x %d cells, %.1fs: max |T_D - T_D(Lauwerier)| = %.3f (T_D = 0.5 at x = %.0f m; without heat loss the front would be at %.0f m)\n",
                nx, nz, secs, worst, front, mcp * t / (rcA * H * dy));
    check(ok && worst < 0.05, "hot-water front with conductive cap-rock heat loss matches Lauwerier's solution (|dT_D| < 0.05)");
}
// C11 verification: Terzaghi 1D consolidation (iterative fixed-stress coupling of the thermal core with mechanics)
static void testTerzaghi() {
    using namespace ressim;
    std::printf("== C11 verification: Terzaghi consolidation (surface load on a drained-top column)\n");
    const int nz = 40; const double dz = 0.5, H = nz * dz, k = 1e-16, load = 1e6, p0 = 5e6;
    thermal::Case c; c.nx = 1; c.ny = 1; c.nz = nz; c.gravity = false; c.cr = 0;
    c.dx = {10}; c.dy = {10}; c.dz = {dz}; c.poro = {0.2}; c.permx = {k}; c.permz = {k};
    c.p0 = {p0}; c.T0 = {300}; c.sw0 = {1.0}; c.thconr = {2.0}; c.thconsf = 0;
    c.krw0 = 1.0; c.swc = 0.0; c.sorw = 0.0;   // single-phase water: kr = 1
    c.geo.enabled = true; c.geo.E = {5e9}; c.geo.nu = {0.25}; c.geo.biot = {1.0}; c.geo.alphaT = {0.0};
    c.geo.coupling = 1; c.geo.maxIter = 30; c.geo.tol = 1e-9;
    c.geo.bcType[4] = 1; c.geo.bcStress[4] = load;
    thermal::Well dr; dr.name = "DRAIN"; dr.cells = {0}; dr.wi = {1e-9}; dr.bhp = p0;   // drained top (cell 0)
    c.wells = {dr};
    // analytic: dp0 = alpha s / (alpha^2 + M S), c_v = k / (mu (S + alpha^2 / M)); drained at the centre of cell 0
    const double M = 6e9;
    const steam::Props a = steam::liquid(300, p0 - 1e5), b = steam::liquid(300, p0 + 1e5);
    const double cw = -(b.v - a.v) / (0.5 * (a.v + b.v)) / 2e5, S = 0.2 * cw;
    const double mu = steam::viscosity(300, 1.0 / steam::liquid(300, p0).v);   // Pa s
    const double dp0 = load / (1.0 + M * S), cv = k / (mu * (S + 1.0 / M)), He = H - 0.5 * dz;
    auto terz = [&](double zc, double t) {   // zc: distance below the drained point
        double u = 0;
        for (int m = 0; m < 200; ++m) {
            const double Mm = M_PI * (2 * m + 1) / 2;
            u += 2.0 / Mm * std::sin(Mm * zc / He) * std::exp(-Mm * Mm * cv * t / (He * He));
        }
        return dp0 * u;
    };
    if (std::getenv("TERZ_NOGEO")) {   // debug: the same diffusion problem without mechanics (rock compressibility alpha^2 / M)
        c.geo.enabled = false; c.cr = 1.0 / M; c.p0.assign(size_t(nz), p0 + dp0); c.p0[0] = p0;
    }
    thermal::Simulator sim; std::string err;
    bool ok = sim.setup(c, &err);
    if (!ok) std::printf("  setup: %s\n", err.c_str());
    double worst = 0, tv[3] = {0.05, 0.2, 0.5};
    for (int r = 0; r < 3 && ok; ++r) {
        const double t = tv[r] * He * He / cv;
        ok = sim.advanceTo(t, std::max(60.0, t / 200));
        if (!ok) break;
        double e = 0;
        for (int i = 1; i < nz; ++i) {
            const double zc = (i + 0.5) * dz - 0.5 * dz;
            e = std::max(e, std::fabs((sim.pressure()[size_t(i)] - p0) - terz(zc, t)) / dp0);
        }
        std::printf("  T_v = %.2f (t = %.2f d): max |dp - dp_Terzaghi| / dp0 = %.3f (bottom: sim %.4f MPa, analytic %.4f MPa)\n", tv[r], t / 86400, e,
                    (sim.pressure()[size_t(nz - 1)] - p0) * 1e-6, terz((nz - 0.5) * dz - 0.5 * dz, t) * 1e-6);
        worst = std::max(worst, e);
    }
    std::printf("  undrained response dp0 = %.4f MPa (Skempton-type, fluid compressibility %.2e 1/Pa), c_v = %.3e m2/s, steps %d, mechanics solves %d\n",
                dp0 * 1e-6, cw, cv, sim.stats().steps, sim.stats().mechSolves);
    check(ok && worst < 0.03, "Terzaghi consolidation: excess pore pressure within 3% of the series solution (T_v 0.05 / 0.2 / 0.5)");
}
// C16 verification: elastoplastic geomechanics (Drucker-Prager matched to Mohr-Coulomb) - plane-strain platen test
static void testGeomechPlastic() {
    using namespace ressim;
    std::printf("== C16 geomechanics plasticity: plane-strain platen test (Drucker-Prager, plane-strain match)\n");
    const double c0 = 2e6, phiDeg = 30, top = 1000, h = 10, p0 = 1e7;
    auto run = [&](double psiDeg, double H, double cRes, double permB, std::vector<double>& s1, std::vector<double>& ep, std::vector<double>& km, double& s3) {
        geomech::Input in; in.enabled = true; in.E = {10e9}; in.nu = {0.25}; in.biot = {1.0}; in.alphaT = {0.0};
        in.k0 = 1.0; in.cohesion = c0; in.frictionDeg = phiDeg; in.tensile = 0;
        in.plasticity = 1; in.dpMatch = 2; in.dilationDeg = psiDeg; in.hardening = H; in.residualCohesion = cRes;
        in.permPlastic = permB; in.permPlasticMax = 1e3; in.plasticIter = 400; in.plasticTol = 1e-9;
        in.bcType[1] = 1;                                   // x+ free (constant confining stress)
        in.bcType[4] = 3; in.plateDispControl[4] = true;    // rigid platen on top, displacement controlled
        geomech::Solver g; std::string e;
        if (!g.setup(in, 1, 1, 1, {h}, {h}, {h}, top, {p0}, {0.2}, &e)) { std::printf("  setup: %s\n", e.c_str()); return false; }
        const std::vector<double> z(1, 0.0);
        s3 = g.sigEffH()[0];
        const double A = h * h;
        for (int i = 1; i <= 150; ++i) {
            const double U = 0.15 * i / 150;  // up to 1.5 % axial strain (yield at ~0.3 %, the cone then evolves to the plane-strain limit)
            g.setFaceLoad(4, 0, U);
            if (!g.solve(z, z)) return false;
            g.commit();
            s1.push_back(s3 + g.plateForce(4) / A); ep.push_back(g.plasticStrain()[0]); km.push_back(g.permMultiplier()[0]);
        }
        return true;
    };
    const double sp = std::sin(phiDeg * M_PI / 180), Np = (1 + sp) / (1 - sp);
    std::vector<double> s1, ep, km; double s3 = 0;
    bool ok = run(phiDeg, 0, 0, 0, s1, ep, km, s3);
    const double want = Np * s3 + 2 * c0 * std::sqrt(Np);
    const double e1 = ok ? std::fabs(s1.back() - want) / want : 1;
    // elastic part: plane strain, free lateral face: d sigma_z / d eps_z = E / (1 - nu^2)
    const double eEl = ok ? std::fabs((s1[0] - s3) / (0.15 / 150 / h) - 10e9 / (1 - 0.0625)) / (10e9 / 0.9375) : 1;
    std::printf("  associated, perfect plasticity: sigma'_3 %.3f MPa, limit sigma'_1 %.4f MPa vs Mohr-Coulomb %.4f MPa (err %.2e); elastic slope err %.1e\n",
                s3 * 1e-6, ok ? s1.back() * 1e-6 : 0.0, want * 1e-6, e1, eEl);
    check(ok && e1 < 2e-3 && eEl < 1e-6, "plasticity: plane-strain DP collapse stress equals Mohr-Coulomb (associated flow); exact elastic branch");
    // softening to the residual cohesion with damage permeability
    std::vector<double> t1, tp, tk; double t3 = 0;
    const bool ok2 = run(phiDeg, -2e9, 0.5e6, 200, t1, tp, tk, t3);
    const double wantR = Np * t3 + 2 * 0.5e6 * std::sqrt(Np);
    double peak = 0; for (double v : t1) peak = std::max(peak, v);
    const double e2 = ok2 ? std::fabs(t1.back() - wantR) / wantR : 1;
    std::printf("  softening H = -2 GPa: peak %.3f MPa, final %.4f MPa vs residual Mohr-Coulomb %.4f MPa (err %.2e); eps_p %.4f, k/k0 %.2f\n",
                peak * 1e-6, ok2 ? t1.back() * 1e-6 : 0.0, wantR * 1e-6, e2, ok2 ? tp.back() : 0.0, ok2 ? tk.back() : 0.0);
    // softening starts where the stress path first meets the cone (before the out-of-plane stress has adjusted), so
    // the drop after first yield and the final residual Mohr-Coulomb strength are the checks
    bool drop = false; for (size_t i = 1; ok2 && i < t1.size(); ++i) drop = drop || (t1[i] < t1[i - 1] - 1e5 && tp[i - 1] > 0);
    check(ok2 && drop && e2 < 5e-3 && tk.back() > 1.5, "plasticity: cohesion softening (stress drop) to the residual Mohr-Coulomb strength; damage permeability grows");
    // deck round trip of GEOPLAST / rigid-plate MECHBC
    Case c = makeExample("sagd_thm");
    c.thermal.geo.plasticity = 1; c.thermal.geo.dpMatch = 0; c.thermal.geo.dilationDeg = 10; c.thermal.geo.hardening = -1e7;
    c.thermal.geo.residualCohesion = 3e5; c.thermal.geo.permPlastic = 50; c.thermal.geo.bcType[4] = 3; c.thermal.geo.bcDisp[4] = 0.01; c.thermal.geo.plateDispControl[4] = true;
    const Case b = eclipseDeckFromString(eclipseDeckToString(c));
    const auto& G = b.thermal.geo;
    const bool rt = G.plasticity == 1 && G.dpMatch == 0 && G.dilationDeg == 10 && G.hardening == -1e7 && G.residualCohesion == 3e5 && G.permPlastic == 50 &&
                    G.bcType[4] == 3 && G.bcDisp[4] == 0.01 && G.plateDispControl[4];
    check(rt, "plasticity: GEOPLAST and rigid-plate MECHBC deck round trip");
    {   // coupled: SAGD THM with weak caprock/reservoir (thermal stress drives shear yielding near the chamber)
        Case w = makeExample("sagd_thm");
        w.thermal.geo.plasticity = 1; w.thermal.geo.cohesion = 0.3e6; w.thermal.geo.frictionDeg = 25; w.thermal.geo.dilationDeg = 5;
        w.thermal.geo.hardening = -5e8; w.thermal.geo.residualCohesion = 0.1e6; w.thermal.geo.permPlastic = 300;
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        int npl = 0; double kmax = 1, epmax = 0, e0 = 0, e1 = 0, hin = 0, loss = 0, gross = 0; bool first = true;
        const bool ran = thermal::deckModel(w, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
            if (first) { e0 = s.energyInPlace(); first = false; }
            e1 = s.energyInPlace(); hin = s.cumEnthalpyIn(); loss = s.cumHeatLoss();
            gross = 0; for (const auto& r : s.wellReports()) gross += std::fabs(r.cumH);
            if (const geomech::Solver* g = s.mechanics()) {
                npl = g->plasticPoints();
                for (double v : g->permMultiplier()) kmax = std::max(kmax, v);
                for (double v : g->plasticStrain()) epmax = std::max(epmax, v);
            }
            return true;
        }, &err, 1, &st);
        const double ebal = std::fabs((e1 - e0) - (hin - loss)) / std::max(1.0, gross + std::fabs(loss));
        std::printf("  SAGD THM + plasticity 91 d: %d points yielding now, max eps_p %.2e, max k/k0 %.3f, energy balance %.1e, %d steps, %d mechanics solves %s\n",
                    npl, epmax, kmax, ebal, st.steps, st.mechSolves, ran ? "" : err.c_str());
        check(ran && epmax > 1e-5 && kmax > 1.01 && ebal < 1e-5, "plasticity in the coupled THM run: shear yielding, damage permeability, energy balance kept");
    }
}

// C16 verification: Mandel's problem (2-D consolidation under a rigid frictionless plate; Mandel-Cryer effect)
static void testMandel() {
    using namespace ressim;
    std::printf("== C16 verification: Mandel's problem (rigid plate, fixed-stress coupled thermal core + mechanics)\n");
    const int nx = 41, nz = 2, n = nx * nz;
    const double a = 10, k = 1e-15, E = 1e9, nu = 0.2, phi = 0.2, pb = 5e6, pu = 1e6;
    thermal::Case c; c.nx = nx; c.ny = 1; c.nz = nz; c.gravity = false; c.cr = 0;
    c.dx.assign(size_t(n), (a - 0.01) / (nx - 1)); for (int kk = 0; kk < nz; ++kk) c.dx[size_t(kk * nx + nx - 1)] = 0.01;
    c.dy = {1}; c.dz = {1}; c.poro = {phi}; c.permx = {k}; c.permy = {k}; c.permz = {k};
    c.p0 = {pb + pu}; c.T0 = {300}; c.sw0 = {1.0}; c.thconr = {2.0}; c.thconsf = 0;
    c.krw0 = 1.0; c.swc = 0.0; c.sorw = 0.0;
    c.geo.enabled = true; c.geo.E = {E}; c.geo.nu = {nu}; c.geo.biot = {1.0}; c.geo.alphaT = {0.0};
    c.geo.coupling = 1; c.geo.maxIter = 40; c.geo.tol = 1e-10;
    c.geo.bcType[1] = 1;                 // x = a: free, drained
    c.geo.bcType[4] = 3;                 // rigid plate, net force increment 0 (load applied at t = 0+)
    thermal::Well dr; dr.name = "DRAIN"; dr.bhp = pb;
    for (int kk = 0; kk < nz; ++kk) { dr.cells.push_back(kk * nx + nx - 1); dr.wi.push_back(1e-9); }
    c.wells = {dr};
    const steam::Props A = steam::liquid(300, pb + pu - 1e5), B = steam::liquid(300, pb + pu + 1e5);
    const double cw = -(B.v - A.v) / (0.5 * (A.v + B.v)) / 2e5;
    const double mu = steam::viscosity(300, 1.0 / steam::liquid(300, pb + pu).v);
    const double M = 1 / (phi * cw), K = E / (3 * (1 - 2 * nu)), G = E / (2 * (1 + nu)), Ku = K + M;
    const double Bs = M / Ku, nuU = (3 * Ku - 2 * G) / (2 * (3 * Ku + G));
    const double cc = 2 * (k / mu) * Bs * Bs * G * (1 - nu) * (1 + nuU) * (1 + nuU) / (9 * (1 - nuU) * (nuU - nu));
    const double Ar = (1 - nu) / (nuU - nu);
    std::vector<double> al;
    for (int m = 1; m <= 120; ++m) {
        double lo = (m - 1) * M_PI + 1e-9, hi = (m - 1) * M_PI + M_PI / 2 - 1e-9;
        auto f = [&](double x) { return std::sin(x) - Ar * x * std::cos(x); };
        double flo = f(lo);
        for (int it = 0; it < 200; ++it) { const double md = 0.5 * (lo + hi), fm = f(md); if ((fm < 0) == (flo < 0)) { lo = md; flo = fm; } else hi = md; }
        al.push_back(0.5 * (lo + hi));
    }
    auto mandel = [&](double x, double t) {
        double s = 0;
        for (double q : al) s += std::sin(q) / (q - std::sin(q) * std::cos(q)) * (std::cos(q * x / a) - std::cos(q)) * std::exp(-q * q * cc * t / (a * a));
        return 2 * s;
    };
    thermal::Simulator sim; std::string err;
    bool ok = sim.setup(c, &err);
    if (!ok) std::printf("  setup: %s\n", err.c_str());
    const double tds[4] = {0.02, 0.1, 0.4, 1.0};
    double worst = 0, peakSim = 0, peakAna = 0;
    for (int r = 0; r < 4 && ok; ++r) {
        const double t = tds[r] * a * a / cc;
        ok = sim.advanceTo(t, std::max(1.0, t / 100));
        if (!ok) break;
        double e = 0;
        for (int i = 0; i < nx - 1; ++i) {
            double xc = 0; for (int q = 0; q < i; ++q) xc += c.dx[size_t(q)];
            xc += 0.5 * c.dx[size_t(i)];
            const double sv = (sim.pressure()[size_t(i)] - pb) / pu, an = mandel(xc, t);
            e = std::max(e, std::fabs(sv - an));
        }
        const double s0 = (sim.pressure()[0] - pb) / pu, a0 = mandel(0.5 * c.dx[0], t);
        peakSim = std::max(peakSim, s0); peakAna = std::max(peakAna, a0);
        std::printf("  t_D = %.2f: centre p/p0 sim %.4f analytic %.4f; max |error| over x %.4f\n", tds[r], s0, a0, e);
        worst = std::max(worst, e);
    }
    const geomech::Solver* g = sim.mechanics();
    std::printf("  B = %.3f, nu_u = %.3f, c = %.3e m2/s; Mandel-Cryer rise: sim %.4f analytic %.4f; plate settlement %.3e m; steps %d, mechanics solves %d\n",
                Bs, nuU, cc, peakSim, peakAna, g ? g->plateDisp(4) : 0.0, sim.stats().steps, sim.stats().mechSolves);
    check(ok && worst < 0.03 && peakSim > 1.0 && std::fabs(peakSim - peakAna) < 0.02, "Mandel: pore pressure within 3% of the analytic series, Mandel-Cryer rise above the initial pressure reproduced");
}
// C16: mixed-precision BILU0, the GPU block-solver interface (float32 device algorithm emulated on the CPU) with
// double iterative refinement, and their use inside the black-oil Newton loop
static BSR makeBlockLaplace(int nx, int ny, int nb, unsigned seed) {
    std::mt19937 rng(seed); std::uniform_real_distribution<double> U(0.5, 2.0);
    BSR A; A.n = nx * ny; A.nb = nb; A.rowPtr.assign(1, 0);
    std::vector<double> kx(size_t(nx) * ny), ky(size_t(nx) * ny);
    for (auto& v : kx) v = std::pow(10.0, 2 * (U(rng) - 1.25)); for (auto& v : ky) v = std::pow(10.0, 2 * (U(rng) - 1.25));   // 2 decades
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const int c = j * nx + i;
        std::vector<std::pair<int, double>> nbr;
        if (j > 0) nbr.push_back({c - nx, std::sqrt(ky[c] * ky[c - nx])});
        if (i > 0) nbr.push_back({c - 1, std::sqrt(kx[c] * kx[c - 1])});
        if (i < nx - 1) nbr.push_back({c + 1, std::sqrt(kx[c] * kx[c + 1])});
        if (j < ny - 1) nbr.push_back({c + nx, std::sqrt(ky[c] * ky[c + nx])});
        double diag = 0;
        for (auto& q : nbr) diag += q.second;
        diag *= 1.02;                                  // accumulation term (compressibility * PV / dt) as in the scaled Jacobian
        std::vector<std::pair<int, double>> row = nbr; row.push_back({c, diag});
        std::sort(row.begin(), row.end());
        for (auto& q : row) {
            A.col.push_back(q.first);
            for (int e = 0; e < nb; ++e) for (int f = 0; f < nb; ++f) {
                double v = 0;
                if (q.first == c) v = (e == f) ? diag * (1.0 + 0.3 * e) : 0.05 * diag * (e == 0 ? 1.0 : 0.2) * (f + 1);   // pressure-coupled block
                else v = (e == f) ? -q.second * (e == 0 ? 1.0 : 0.6) : (e == 0 ? -0.02 * q.second : 0.0);
                A.val.push_back(v);
            }
        }
        A.rowPtr.push_back(int(A.col.size()));
    }
    A.finalizeDiag();
    return A;
}

static void testMixedPrecision() {
    std::printf("== C16 mixed-precision BILU0 + GPU block-solver interface (float32 device emulation)\n");
    const BSR A = makeBlockLaplace(60, 50, 3, 7);
    const int N = A.size();
    std::vector<double> xt(N), b(N);
    for (int e = 0; e < N; ++e) xt[e] = std::sin(0.01 * e) + 0.1 * (e % 3);
    A.spmv(xt.data(), b.data());
    BsrOperator op(A);
    BILU0 M; M.factor(A);
    BILU0f Mf; Mf.factor(A);
    std::vector<double> x1, x2;
    const LinResult r1 = bicgstab(op, M, b, x1, 1e-9, 2000), r2 = bicgstab(op, Mf, b, x2, 1e-9, 2000);
    double e1 = 0, e2 = 0, xn = 0;
    for (int e = 0; e < N; ++e) { e1 = std::max(e1, std::fabs(x1[e] - xt[e])); e2 = std::max(e2, std::fabs(x2[e] - xt[e])); xn = std::max(xn, std::fabs(xt[e])); }
    std::printf("  BILU0 double: %d its, err %.2e | BILU0f: %d its, err %.2e | factor memory %.0f%%\n", r1.iterations, e1 / xn, r2.iterations,
                e2 / xn, 100.0 * double(Mf.bytes()) / double(A.val.size() * sizeof(double)));
    check(r1.converged && r2.converged && e2 / xn < std::max(1e-7, 2.0 * e1 / xn) && r2.iterations <= r1.iterations + std::max(2, r1.iterations / 10),
          "float32 BILU0 preconditioner: same accuracy as the double BILU0 (solver tolerance), iteration count within 10%");
    CpuFloatBlockSolver dev;
    std::vector<double> xd; int its = 0;
    dev.solve(A, b, xd, 1e-12, 3000, its);
    double ed = 0; for (int e = 0; e < N; ++e) ed = std::max(ed, std::fabs(xd[e] - xt[e]));
    std::vector<double> xr;
    const LinResult rr = refineWithDevice(dev, A, b, xr, 1e-10, 8, 2000);
    double er = 0; for (int e = 0; e < N; ++e) er = std::max(er, std::fabs(xr[e] - xt[e]));
    std::printf("  pure float32 device solve: err %.2e (float limit) | with double refinement: res %.2e, err %.2e, %d inner its\n", ed / xn, rr.relResidual,
                er / xn, rr.iterations);
    check(ed / xn > 1e-8 && rr.converged && er / xn < 1e-8, "iterative refinement recovers double accuracy from the float32 device solver");
    // inside the Newton loop
    Case c = makeDemoCase();
    c.schedule.resize(6);
    c.tuning.blockSolver = true;
    double fopt[3], fpr[3]; int nw[3]; RunStats st[3];
    for (int k = 0; k < 3; ++k) {
        Case ck = c;
        ck.tuning.mixedPrecision = (k == 1);
        ck.tuning.gpu = (k == 2);
        setGpuBlockSolver(k == 2 ? &dev : nullptr);
        Simulator sim; sim.setup(ck); sim.initialize();
        while (sim.runReportStep()) {}
        fopt[k] = sim.summary().back().fopt; fpr[k] = sim.summary().back().fpr; nw[k] = sim.stats().newton; st[k] = sim.stats();
    }
    setGpuBlockSolver(nullptr);
    std::printf("  demo 6 steps: FOPT %.6g / %.6g / %.6g, Newton %d / %d / %d, device solves %d (fallbacks %d)\n", fopt[0], fopt[1], fopt[2], nw[0], nw[1],
                nw[2], st[2].gpuSolves, st[2].gpuFallbacks);
    check(std::fabs(fopt[1] - fopt[0]) < 1e-5 * fopt[0] && std::fabs(fpr[1] - fpr[0]) < 1e-5 * fpr[0], "mixed precision run = double run");
    check(st[2].gpuSolves > 0 && st[2].gpuFallbacks <= st[2].gpuSolves / 5 && std::fabs(fopt[2] - fopt[0]) < 1e-5 * fopt[0],
          "GPU-interface run (device emulation + refinement) = CPU run");
    Case r = c; r.tuning.mixedPrecision = true; r.tuning.gpu = 1;
    const Case r2d = eclipseDeckFromString(eclipseDeckToString(r));
    const Case r2p = projectFromText(projectText(r));
    check(r2d.tuning.mixedPrecision && r2d.tuning.gpu == 1 && r2p.tuning.mixedPrecision && r2p.tuning.gpu == 1, "TUNING / project round trip keeps the precision and GPU options");
}

// C16: run-time MPI loader, distributed batch protocol (single rank in-process; multi-rank via mpiexec in the release
// check), and the study batch-runner hook
static void testMpiRuntime() {
    std::printf("== C16 MPI runtime + distributed batch + study batch-runner hook\n");
    // study hook
    int calls = 0;
    setBatchRunner([&](const std::vector<Case>& cs, std::atomic<bool>*, const std::function<void(int, int)>& prog) {
        ++calls; std::vector<RunResult> out;
        for (size_t i = 0; i < cs.size(); ++i) { out.push_back(runSingle(cs[i], nullptr)); if (prog) prog(int(i) + 1, int(cs.size())); }
        return out;
    });
    Case c = makeDemoCase(); c.schedule.resize(3);
    std::vector<Case> cs = {c, c};
    cs[1].wells[0].piMult = 2.0;
    auto rs = runBatch(cs, 2, nullptr, nullptr);
    setBatchRunner(nullptr);
    check(calls == 1 && rs.size() == 2 && rs[0].ok && rs[1].ok && !hasBatchRunner(), "runBatch delegates to the registered (MPI) batch runner");
    std::string err;
    const bool have = mpi::load(&err);
    std::printf("  MPI library: %s\n", have ? mpi::libraryName().c_str() : err.c_str());
    if (!have) { check(!err.empty(), "no MPI installed: loader reports a clear message (skipping MPI run)"); return; }
    if (std::getenv("RESSIM_TEST_NO_MPI_INIT")) { check(true, "MPI library binds (init skipped)"); return; }
    check(mpi::init() && mpi::size() >= 1, "MPI_Init through the run-time loaded library");
    namespace fs = std::filesystem;
    const fs::path jd = fs::temp_directory_path() / "ressim_mpi_job";
    fs::remove_all(jd); fs::create_directories(jd);
    for (int i = 0; i < 2; ++i) { std::ofstream f(jd / ("case_000" + std::to_string(i) + ".rsim")); f << projectText(cs[size_t(i)]); }
    const int rc = mpi::runBatchJob(jd.string(), false);
    bool ok = rc == 0;
    for (int i = 0; i < 2 && ok; ++i) {
        const ResultsFile rf = loadResultsFile((jd / ("case_000" + std::to_string(i) + ".rsres")).string(), false);
        const double ref = rs[size_t(i)].summary.back().fopt, got = rf.summary.empty() ? -1 : rf.summary.back().fopt;
        std::printf("  case %d: %s FOPT %.6g (direct %.6g)\n", i, rf.title.c_str(), got, ref);
        ok = rf.title == "OK" && std::fabs(got - ref) < 1e-9 * std::max(1.0, ref);
    }
    check(ok && fs::exists(jd / "done.txt"), "MPI batch job (rank 0 alone) writes identical results");
    fs::remove_all(jd);
    mpi::finalize();
}

// C12: Langmuir sorption (competitive gas adsorption, nanoparticle retention, deposition permeability loss)
static void testSorption() {
    using namespace ressim;
    std::printf("== C12 Langmuir sorption (mass-action site model in the thermal core)\n");
    {   // closed coal cell: CH4 / CO2 competitive (extended Langmuir) adsorption with heat of adsorption
        thermal::Case c; c.nx = 1; c.dx = {10}; c.dy = {10}; c.dz = {10}; c.poro = {0.1}; c.permx = {1e-14}; c.gravity = false; c.cr = 1e-10;
        thermal::Component m; m.name = "CH4"; m.mw = 0.016; m.type = 1; m.cpg1 = 36;
        thermal::Component d; d.name = "CO2"; d.mw = 0.044; d.type = 1; d.cpg1 = 37;
        c.comps = {m, d};
        c.p0 = {5e6}; c.T0 = {320}; c.sw0 = {0.2}; c.sg0 = {0.8}; c.x0 = {{0.7, 0.3}}; c.y0 = {{0.7, 0.3}};
        thermal::Sorption a; a.comp = "CH4"; a.capacity = 100; a.bRef = 1 / 6690.0; a.tRef = 320; a.dH = 15e3; a.tau = 3600; a.initEquilibrium = false;
        thermal::Sorption b = a; b.comp = "CO2"; b.bRef = 1 / 2000.0; b.dH = 25e3;
        std::string err;
        bool ok = thermal::addSorption(c, a, &err) && thermal::addSorption(c, b, &err);
        thermal::Simulator s;
        ok = ok && s.setup(c, &err);
        const int nc = ok ? s.nc() : 0, nb = nc + 3;
        auto N = [&](int k) { return s.moles()[size_t(1 + k)]; };
        // components: CH4 0, CO2 1, SITES 2, CH4_ADS 3, CO2_ADS 4
        const double m0 = ok ? N(0) + N(3) : 0, d0 = ok ? N(1) + N(4) : 0, E0 = ok ? s.energyInPlace() : 0;
        const double cap = ok ? (N(2) + N(3) + N(4)) / s.bulkVolume(0) : 0;   // sites are set per initial pore volume (cr p0)
        const double jac = ok ? s.jacobianCheck(3600.0) : 1;
        ok = ok && s.advanceTo(20 * 86400.0, 86400.0);
        double eIso = 1, eM = 1, eE = 1, eS = 1;
        if (ok) {
            std::vector<double> x, y, aq; s.phaseCompositions(x, y, aq);
            const double p = s.pressure()[0] * 1e-3, T = s.temperature()[0], Vb = s.bulkVolume(0);
            const double bm = a.bRef * std::exp(a.dH / 8.314462618 * (1 / T - 1 / a.tRef)), bd = b.bRef * std::exp(b.dH / 8.314462618 * (1 / T - 1 / b.tRef));
            const double pm = y[0] * p, pd = y[1] * p, den = 1 + bm * pm + bd * pd;
            const double qm = cap * bm * pm / den, qd = cap * bd * pd / den;
            eIso = std::max(std::fabs(N(3) / Vb - qm) / qm, std::fabs(N(4) / Vb - qd) / qd);
            eM = std::max(std::fabs(N(0) + N(3) - m0) / m0, std::fabs(N(1) + N(4) - d0) / d0);
            eS = std::fabs(N(2) + N(3) + N(4) - cap * Vb) / (cap * Vb);
            eE = std::fabs(s.energyInPlace() - E0 - s.cumReactionHeat()) / std::max(1.0, std::fabs(s.cumReactionHeat()));
            std::printf("  capacity %.4f mol/m3 bulk (100 at the reference PV)\n", cap);
            std::printf("  coal cell: p 5 -> %.3f MPa, T %.2f K, sorbed CH4 %.2f / CO2 %.2f mol/m3 (extended Langmuir %.2f / %.2f), selectivity CO2/CH4 %.2f\n",
                        p * 1e-3, T, N(3) / Vb, N(4) / Vb, qm, qd, (qd / qm) / (pd / pm));
        }
        std::printf("  isotherm error %.1e, component mass drift %.1e, site balance %.1e, energy vs adsorption heat %.1e, Jacobian %.1e %s\n",
                    eIso, eM, eS, eE, jac, ok ? "" : err.c_str());
        (void)nb;
        check(ok && eIso < 1e-5 && eM < 1e-10 && eS < 1e-10 && eE < 1e-4 && jac < 1e-3,
              "competitive (extended) Langmuir adsorption: isotherm at the final p, T; exact mass and site balances, energy = adsorption heat (Newton tol.); Jacobian");
    }
    {   // nanoparticles: 1D column, Langmuir retention -> self-sharpening front at u C0 / (phi C0 + q(C0)); deposition -> Kozeny-Carman
        double tbt[2] = {0, 0}, tan = 0, kmErr = 1; bool ok = true;
        for (int mode = 0; mode < 2; ++mode) {
            thermal::Case c; c.nx = 100; c.dx = {1}; c.dy = {1}; c.dz = {1}; c.poro = {0.2}; c.permx = {1e-12}; c.gravity = false; c.cr = 0;
            c.krw0 = 1; c.swc = 0; c.sorw = 0;
            thermal::Component np; np.name = "NP"; np.mw = mode == 1 ? 100.0 : 0.1; np.type = 1; np.aqModel = 2; np.kh298 = 1e3; np.vInf = 0;
            c.comps = {np}; c.p0 = {1e7}; c.T0 = {300}; c.sw0 = {1.0}; c.x0 = {{0.0}};
            if (mode == 1) { c.permLaw = 1; c.permMultMax = 1e5; }
            const double z = 1e-5, rate = 1.0;
            thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {1e-11}; inj.zinj = {z, 1 - z}; inj.rate = rate; inj.bhpMax = 5e7;
            inj.quality = -1; inj.tinj = 300;
            thermal::Well pr; pr.name = "PROD"; pr.cells = {99}; pr.wi = {1e-11}; pr.bhp = 1e7;
            c.wells = {inj, pr};
            thermal::Sorption so; so.comp = "NP"; so.capacity = 0.2; so.bRef = 20; so.tRef = 300; so.dH = 0; so.tau = 300; so.phase = 0;
            so.adsRho = mode == 1 ? 2000 : 0;
            std::string err;
            thermal::Simulator s;
            bool r = thermal::addSorption(c, so, &err) && s.setup(c, &err);
            const double vw = 1.0 / steam::liquid(300, 1e7).v / 0.018015;   // mol / m3 of water
            const double C0 = z * vw, cb = 0.2 * C0, q0 = 0.2 * 20 * cb / (1 + 20 * cb), u = rate / vw;
            tan = 99.0 * (0.2 * C0 + q0) / (u * C0);
            double tPrev = 0, cPrev = 0;
            for (int k = 1; r && k <= 60; ++k) {
                const double t = k * 86400.0;
                r = s.advanceTo(t, 0.25 * 86400.0);
                if (!r) break;
                std::vector<double> x, y, aq; s.phaseCompositions(x, y, aq);
                const double cn = aq[size_t(99) * 4 + 0] * vw / C0;
                if (cn >= 0.5 && cPrev < 0.5 && tbt[mode] == 0) tbt[mode] = tPrev + (0.5 - cPrev) / (cn - cPrev) * (t - tPrev);
                tPrev = t; cPrev = cn;
            }
            if (r && mode == 1) {
                const double f = s.fluidPorosity()[0], f0 = 0.2;
                const double want = std::pow(f / f0, 3) * std::pow((1 - f0) / (1 - f), 2);
                kmErr = std::fabs(s.permMultiplier()[0] - want) / want;
                std::printf("  deposition (100 kg/mol clusters, 2000 kg/m3): inlet phi_f %.4f, k/k0 %.4f (Kozeny-Carman %.4f)\n", f, s.permMultiplier()[0], want);
            }
            if (!r) std::printf("  mode %d failed: %s\n", mode, err.c_str());
            ok = ok && r;
        }
        const double e = std::fabs(tbt[0] - tan) / tan;
        std::printf("  nanoparticle breakthrough (C/C0 = 0.5 at x = 99 m): %.2f d vs Langmuir shock %.2f d (error %.1f %%; unretarded %.2f d)\n",
                    tbt[0] / 86400, tan / 86400, 100 * e, tan / 86400 / 2.25);
        check(ok && e < 0.03 && kmErr < 1e-6, "nanoparticle Langmuir retention: shock arrival within 3 %; deposition permeability loss follows Kozeny-Carman");
    }
    {   // SORB deck keyword: round trip and the CO2-ECBM example (CO2 replaces sorbed CH4)
        const Case c = makeExample("ecbm");
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        const bool rt = b.thermal.sorb.size() == 2 && b.thermal.sorb[1].comp == "CO2" && std::fabs(b.thermal.sorb[1].b - 1 / 2900.0) < 1e-12 &&
                        b.thermal.sorb[0].tau == 2 * 86400.0 && b.thermal.sorb[0].site == "SITES";
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out;
        double ch4Ads0 = -1, ch4Ads = 0, co2Ads = 0, ch4Prod = 0, mb = 1; int nc = 0;
        const bool ok = thermal::deckModel(b, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator& sim) {
            (void)S; nc = sim.nc();
            if (nc != 5) return true;
            const double a3 = sim.solidMoles(3), a4 = sim.solidMoles(4);
            if (ch4Ads0 < 0) ch4Ads0 = a3;
            ch4Ads = a3; co2Ads = a4;
            ch4Prod = 0; for (const auto& r : sim.wellReports()) if (r.cumMolP.size() > 0) ch4Prod += r.cumMolP[0];
            const double site = sim.solidMoles(2) + a3 + a4;
            static double site0 = -1; if (site0 < 0) site0 = site;
            mb = std::fabs(site - site0) / site0;
            return true;
        }, &err, 6);
        std::printf("  CO2-ECBM: %d components, sorbed CH4 %.3g -> %.3g kmol, sorbed CO2 %.3g kmol, CH4 produced %.3g kmol, site balance %.1e, SORB round trip %d %s\n",
                    nc, ch4Ads0 * 1e-3, ch4Ads * 1e-3, co2Ads * 1e-3, ch4Prod * 1e-3, mb, int(rt), ok ? "" : err.c_str());
        check(ok && rt && nc == 5 && ch4Ads < 0.8 * ch4Ads0 && co2Ads > 0 && ch4Prod > 0 && mb < 1e-9,
              "SORB deck keyword round trip; CO2-ECBM example: CO2 sorbs, sorbed CH4 is released and produced, sites conserved");
    }
}
// C12: embedded discrete fractures (EDFM) with Barton-Bandis closure and hydro-shear in the thermal core
// ------------------------------------------------------------------------------------------------ SPE4 problem 1A (C12)
// Fourth SPE comparative solution project (Aziz, Ramesh, Woo, JPT 1987), problem 1A: three cycles of cyclic steam
// (10 d injection at 1000 STB/d CWE, 0.7 quality, 450 F; 7 d soak; 348 d production) in a 13 x 4 r-z grid. Reference
// values quoted in the paper text: OOIP 79,227-80,572 sm3, OWIP 65,693-66,986 sm3, 3,594 GJ injected per cycle, cumulative
// heat loss at the end of cycle 1 491-703 GJ (the plotted rates / cumulatives are not tabulated in the paper).
static ressim::thermal::Case spe4Problem1(std::vector<int>& wellCells) {
    using namespace ressim;
    const double ft = 0.3048, md = 9.869233e-16, psi = 6894.757;
    const double rb[14] = {0.3, 3, 13, 23, 33, 43, 53, 63, 73, 83, 93, 103, 143, 263};
    const double dzF[4] = {10, 20, 25, 25}, khMd[4] = {2000, 500, 1000, 2000};
    const int nr = 13, nzl = 4;
    thermal::Case c;
    std::vector<double> rn(static_cast<size_t>(nr));
    for (int i = 0; i < nr; ++i) { const double a = rb[i] * ft, b = rb[i + 1] * ft; rn[size_t(i)] = (b - a) / std::log(b / a); }   // log-mean node
    const double top = 1500 * ft, rw = 0.3 * ft;
    std::vector<double> zc(nzl), dz(nzl);
    double z = 0;
    for (int k = 0; k < nzl; ++k) { dz[size_t(k)] = dzF[k] * ft; zc[size_t(k)] = z + 0.5 * dz[size_t(k)]; z += dz[size_t(k)]; }
    const double phi = 0.3, pTop = 75 * psi, rhoMix = 972, g = 9.80665;
    auto id = [&](int i, int k) { return k * nr + i; };
    c.cellPV.resize(size_t(nr * nzl)); c.cellBulk.resize(size_t(nr * nzl)); c.cellDepth.resize(size_t(nr * nzl)); c.p0.resize(size_t(nr * nzl));
    for (int k = 0; k < nzl; ++k)
        for (int i = 0; i < nr; ++i) {
            const double a = rb[i] * ft, b = rb[i + 1] * ft, V = M_PI * (b * b - a * a) * dz[size_t(k)];
            const int q = id(i, k);
            c.cellBulk[size_t(q)] = V; c.cellPV[size_t(q)] = phi * V; c.cellDepth[size_t(q)] = top + zc[size_t(k)];
            c.p0[size_t(q)] = pTop + rhoMix * g * (zc[size_t(k)] - zc[0]);
        }
    for (int k = 0; k < nzl; ++k)
        for (int i = 0; i < nr; ++i) {
            const double kh = khMd[k] * md, h = dz[size_t(k)];
            if (i + 1 < nr) {
                const double L = std::log(rn[size_t(i + 1)] / rn[size_t(i)]);
                c.conns.push_back({id(i, k), id(i + 1, k), 2 * M_PI * kh * h / L, 2 * M_PI * h / L});
            }
            if (k + 1 < nzl) {
                const double a = rb[i] * ft, b = rb[i + 1] * ft, A = M_PI * (b * b - a * a);
                const double kv0 = 0.5 * kh, kv1 = 0.5 * khMd[k + 1] * md, h1 = dz[size_t(k + 1)];
                c.conns.push_back({id(i, k), id(i, k + 1), A / (0.5 * h / kv0 + 0.5 * h1 / kv1), A / (0.5 * h + 0.5 * h1)});
            }
        }
    c.heatLoss = true;
    for (int i = 0; i < nr; ++i) {
        const double a = rb[i] * ft, b = rb[i + 1] * ft, A = M_PI * (b * b - a * a);
        c.lossFaces.push_back({id(i, 0), A}); c.lossFaces.push_back({id(i, nzl - 1), A});
    }
    // rock: 24 Btu/(ft d F) = 1.730 W/(m K); 35 Btu/(ft3 F) = 2.347e6 J/(m3 K); 5e-4 1/psi
    const double lam = 24 * 1055.056 / (ft * 86400 * (5.0 / 9.0)), cv = 35 * 1055.056 / (ft * ft * ft * (5.0 / 9.0));
    c.thconr = {lam}; c.thconsf = 0; c.rockCv = cv; c.obCv = cv; c.obK = lam;
    c.cr = 5e-4 / psi; c.pref = pTop;
    // oil: 972 kg/m3 at 60 F, 5e-6 1/psi, 3.8e-4 1/R, 0.5 Btu/(lbm R), MW 600, Table A-1 viscosity
    thermal::Component o; o.name = "OIL"; o.mw = 0.6; o.dead = true; o.rhoRef = 972; o.cpL = 5e-6 / psi; o.ct1 = 3.8e-4 * 1.8;
    o.cpl1 = 0.5 * 4186.8 * 0.6;
    const double tF[8] = {75, 100, 150, 200, 250, 300, 350, 500}, muC[8] = {5780, 1380, 187, 47, 17.4, 8.5, 5.2, 2.5};
    for (int j = 0; j < 8; ++j) { o.viscT.push_back((tF[j] - 32) / 1.8 + 273.15); o.viscV.push_back(muC[j] * 1e-3); }
    c.comps = {o};
    c.tRef = (60 - 32) / 1.8 + 273.15; c.stdT = c.tRef; c.pRef = 101325; c.stdP = 101325;
    c.T0 = {(125 - 32) / 1.8 + 273.15}; c.sw0 = {0.45}; c.sg0 = {0}; c.x0 = {{1.0}};
    // relative permeability (A-1 .. A-4)
    const double swir = 0.45, sorw = 0.15, sorg = 0.10, sgc = 0.06;
    c.swc = swir; c.sorw = sorw; c.sorg = sorg; c.sgc = sgc;
    for (int j = 0; j <= 20; ++j) {
        const double sw = swir + (1 - sorw - swir) * j / 20.0;
        c.swT.push_back(sw);
        c.krwT.push_back(0.1 * std::pow((sw - swir) / (1 - sorw - swir), 2.5));
        c.krowT.push_back(0.4 * std::pow((1 - sorw - sw) / (1 - sorw - swir), 2.0));
    }
    c.swT.push_back(1.0); c.krwT.push_back(0.1 * std::pow((1 - swir) / (1 - sorw - swir), 2.5)); c.krowT.push_back(0.0);
    c.slT.push_back(swir); c.krgT.push_back(0.2 * std::pow((1 - swir - sgc) / (1 - swir - sgc), 1.5)); c.krogT.push_back(0.0);
    for (int j = 0; j <= 20; ++j) {   // SLT: liquid saturation from swir + sorg to 1
        const double sl = swir + sorg + (1 - swir - sorg) * j / 20.0, sg = 1 - sl;
        c.slT.push_back(sl);
        c.krgT.push_back(sg > sgc ? 0.2 * std::pow((sg - sgc) / (1 - swir - sgc), 1.5) : 0.0);
        c.krogT.push_back(0.4 * std::pow(std::max(0.0, (1 - swir - sorg - sg) / (1 - swir - sorg)), 2.0));
    }
    // wells on the axis: J = 2 pi k h / ln(r1 / rw), all layers, BHP datum at the centre of the top layer
    thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.zinj = {0, 1}; inj.quality = 0.7; inj.tinj = (450 - 32) / 1.8 + 273.15;
    inj.rate = 1000 * 0.158987294928 / 86400 * 1000 / steam::kMw; inj.bhpMax = 1000 * psi; inj.refDepth = top + zc[0];
    thermal::Well pr; pr.name = "PROD"; pr.control = 3; pr.target = 1000 * 0.158987294928 / 86400; pr.bhp = 17 * psi; pr.refDepth = top + zc[0];
    wellCells.clear();
    for (int k = 0; k < nzl; ++k) {
        const double wi = 2 * M_PI * khMd[k] * md * dz[size_t(k)] / std::log(rn[0] / rw);
        inj.cells.push_back(id(0, k)); inj.wi.push_back(wi); pr.cells.push_back(id(0, k)); pr.wi.push_back(wi);
        wellCells.push_back(id(0, k));
    }
    c.wells = {inj, pr};
    return c;
}
static void testSpe4() {
    std::printf("SPE4 problem 1A: cyclic steam, r-z grid, 3 cycles (Aziz et al. 1987)\n");
    std::vector<int> wc;
    ressim::thermal::Case c = spe4Problem1(wc);
    ressim::thermal::Simulator s; std::string err;
    bool ok = s.setup(c, &err);
    if (!ok) { std::printf("  setup failed: %s\n", err.c_str()); check(false, "SPE4 1A setup"); return; }
    const double stb = 0.158987294928, day = 86400.0;
    const double ooip = s.stockTankOilInPlace(), owip = s.fluidMoles(1) * ressim::steam::kMw / 999.0;
    std::printf("  OOIP %.0f sm3 (paper 79,227-80,572)  OWIP %.0f sm3 (paper 65,693-66,986)\n", ooip, owip);
    check(ooip > 79227 * 0.99 && ooip < 80572 * 1.01 && owip > 65693 * 0.99 && owip < 66986 * 1.01, "SPE4 1A initial oil / water in place within the participants' range (+-1%)");
    double hinj1 = 0, loss1 = 0, maxBhp = 0;
    double cumO[3] = {0, 0, 0}, cumW[3] = {0, 0, 0};
    auto& W = s.wells();
    auto t0 = std::chrono::steady_clock::now();
    for (int cyc = 0; cyc < 3 && ok; ++cyc) {
        const double tc = cyc * 365 * day;
        W[0].open = true; W[1].open = false; s.limitNextDt(0.05 * day);
        const double hBefore = s.cumEnthalpyIn();
        for (double t = tc + 1 * day; t <= tc + 10 * day + 1 && ok; t += day) {
            ok = s.advanceTo(t, 0.5 * day);
            if (ok) maxBhp = std::max(maxBhp, s.wellReports()[0].bhp);
        }
        if (cyc == 0) { hinj1 = s.cumEnthalpyIn() - hBefore; std::printf("  cycle 1 injection: CWE %.1f m3 (target 1590), %.1f kJ/kg\n", s.wellReports()[0].cumCwe, hinj1 / std::max(1.0, s.wellReports()[0].cumCwe * 1000) / 1e3); }
        W[0].open = false; s.limitNextDt(0.05 * day);
        ok = ok && s.advanceTo(tc + 17 * day, 1 * day);
        W[1].open = true; s.limitNextDt(0.05 * day);
        for (double t = tc + 27 * day; t <= tc + 365 * day + 1 && ok; t += 10 * day) ok = s.advanceTo(std::min(t, tc + 365 * day), 5 * day);
        ok = ok && s.advanceTo(tc + 365 * day, 5 * day);
        if (cyc == 0) loss1 = s.cumHeatLoss();
        cumO[cyc] = s.wellReports()[1].cumO; cumW[cyc] = s.wellReports()[1].cumW;
        double Tmax = 0; for (double T : s.temperature()) Tmax = std::max(Tmax, T);
        std::printf("  cycle %d: cum oil %.0f STB  cum water %.0f STB  heat loss %.0f GJ  Tmax %.1f C\n", cyc + 1, cumO[cyc] / stb, cumW[cyc] / stb,
                    s.cumHeatLoss() / 1e9, Tmax - 273.15);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const ressim::thermal::Stats& st = s.stats();
    std::printf("  heat injected cycle 1 %.0f GJ (paper 3,594), max injection BHP %.0f psia, %d steps %d Newton %d cuts, %.1f s\n", hinj1 / 1e9,
                maxBhp / 6894.757, st.steps, st.newton, st.cuts, secs);
    check(ok, "SPE4 1A runs three cycles");
    // ECLIPSE / CMG injectivity (total block mobility x injected molar density): the cold, viscous well block holds the
    // injector at its BHP limit for the first ~0.5 d (Chevron's kr = 1 / mu_inj convention injects at full rate at once)
    check(ok && hinj1 < 3.594e12 * 1.005 && hinj1 > 3.594e12 * 0.95, "SPE4 1A heat injected in cycle 1 within 5% of 3,594 GJ (IF97 datum; BHP-limited start-up)");
    check(ok && loss1 > 491e9 * 0.9 && loss1 < 703e9 * 1.1, "SPE4 1A cumulative heat loss after cycle 1 within the participants' 491-703 GJ (+-10%)");
    check(ok && cumO[0] > 0 && cumO[1] > cumO[0] && cumO[2] > cumO[1], "SPE4 1A oil produced in every cycle");
}

// ------------------------------------------------------------------------------------------------ SPE11A-like (C12)
// Reduced SPE11A (Nordbotten et al., SPE J 2024): 2.8 m x 1.2 m x 1 cm vertical slice at 20 C, 1.1 bar at the top,
// SPE11A facies 1 (seal), 3, 4 and 5 (k, phi, s_w,imm, entry pressure from Table 1; Brooks-Corey kr exponent 2, Pc exponent
// 2 capped at 9.5e4 Pa, s_n,imm = 0.1, D_w = 1e-9 m2/s), the two SPE11A wells and the SPE11A injection protocol
// (1.7e-7 kg/s; well 1 0-5 h, well 2 2.5-5 h), monitored to 5 days. The real facies geometry is replaced by a layered
// one (lower reservoir F5 / seal F1 / upper reservoir F4 / seal F1 / top F3), so this is a physics check, not a
// reproduction of the benchmark maps.
static void testSpe11a() {
    using namespace ressim; using namespace ressim::comp;
    std::printf("SPE11A-like CO2 storage slice (layered facies, SPE11A fluids / wells / schedule)\n");
    const int nx = std::getenv("SPE11_NX") ? std::atoi(std::getenv("SPE11_NX")) : 28, nz = nx * 3 / 7; const double Lx = 2.8, Lz = 1.2, W = 0.01, dx = Lx / nx, dz = Lz / nz;
    struct Fac { double k, phi, swi, pe; };
    const Fac F1{4e-11, 0.44, 0.32, 1500}, F3{1e-9, 0.44, 0.12, 100}, F4{2e-9, 0.45, 0.12, 25}, F5{4e-9, 0.43, 0.12, 10};
    auto facOf = [&](double zUp) -> int { return zUp < 0.5 ? 5 : zUp < 0.6 ? 1 : zUp < 1.0 ? 4 : zUp < 1.1 ? 1 : 3; };
    auto fac = [&](int f) -> const Fac& { return f == 1 ? F1 : f == 3 ? F3 : f == 4 ? F4 : F5; };
    CompCase c;
    c.fluid.comp = {eos::standardComponent("CO2"), eos::standardComponent("H2O")};
    c.co2store = true; c.salinity = 0; c.T = 293.15; c.cr = 0; c.pref = 1.1e5;
    c.trap.diffusion = true; c.trap.dCo2Water = 1e-9; c.trap.dH2oWater = 1e-9; c.trap.dCo2Gas = 1.6e-5; c.trap.dH2oGas = 1.6e-5;
    const int facList[4] = {1, 3, 4, 5};
    auto table = [&](const Fac& f) {
        SatTables t; const double sgi = 0.1;
        for (int j = 0; j <= 30; ++j) {
            const double sw = f.swi + (1 - f.swi) * j / 30.0, swn = (sw - f.swi) / (1 - f.swi);
            t.sw.push_back(sw); t.krw.push_back(swn * swn);
            t.pcow.push_back(std::min(9.5e4, swn > 0 ? f.pe / std::sqrt(swn) : 9.5e4));
        }
        for (int j = 0; j <= 30; ++j) {
            const double sg = (1 - f.swi) * j / 30.0, sgn = std::max(0.0, (sg - sgi) / (1 - sgi - f.swi));
            t.sg.push_back(sg); t.krg.push_back(sgn * sgn); t.pcog.push_back(0);
        }
        return t;
    };
    c.sat = table(F1);
    for (int r = 1; r < 4; ++r) c.satRegions.push_back(table(fac(facList[r])));
    auto region = [&](int f) { return f == 1 ? 0 : f == 3 ? 1 : f == 4 ? 2 : 3; };
    const double g = 9.80665, pTop = 1.1e5;
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i) {
            const double zUp = Lz - (k + 0.5) * dz;   // k = 0: top row
            const Fac& f = fac(facOf(zUp));
            c.cellPV.push_back(f.phi * dx * dz * W); c.cellDepth.push_back((k + 0.5) * dz);
            c.satnum.push_back(region(facOf(zUp)));
            c.p0.push_back(pTop + 998.2 * g * (k + 0.5) * dz);
        }
    auto id = [&](int i, int k) { return k * nx + i; };
    auto zUpOf = [&](int k) { return Lz - (k + 0.5) * dz; };
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i) {
            const Fac& a = fac(facOf(zUpOf(k)));
            if (i + 1 < nx) c.conns.push_back({id(i, k), id(i + 1, k), a.k * dz * W / dx, a.phi * dz * W / dx});
            if (k + 1 < nz) {
                const Fac& b = fac(facOf(zUpOf(k + 1)));
                const double T = dx * W / (0.5 * dz / a.k + 0.5 * dz / b.k), phiH = 2 * a.phi * b.phi / (a.phi + b.phi);
                c.conns.push_back({id(i, k), id(i, k + 1), T, phiH * dx * W / dz});
            }
        }
    c.z0 = {{1e-9, 1 - 1e-9}};
    auto cellAt = [&](double x, double zUp) { return id(std::min(nx - 1, int(x / dx)), std::min(nz - 1, int((Lz - zUp) / dz))); };
    const double rate = 1.7e-7 / co2::MW_CO2;   // mol/s
    auto mkInj = [&](const char* nm, int cell, const Fac& f) {
        CompWell w; w.name = nm; w.injector = true; w.cells = {cell}; w.wi = {2 * M_PI * f.k * W / std::log(0.2 * dx / 1e-3)};
        w.zinj = {1.0, 0.0}; w.rate = rate; w.bhpMax = 5e5; return w;
    };
    CompWell w1 = mkInj("W1", cellAt(0.9, 0.3), F5), w2 = mkInj("W2", cellAt(1.7, 0.7), F4);
    w2.open = false;
    CompWell top; top.name = "TOP"; top.bhp = pTop; top.refDepth = 0.5 * dz;   // open top boundary (constant pressure)
    for (int i = 0; i < nx; ++i) { top.cells.push_back(id(i, 0)); top.wi.push_back(F3.k * dx * W / (0.5 * dz)); }
    // the producer cannot inject: a BHP-controlled water injector on the same cells completes the Dirichlet boundary
    CompWell topIn = top; topIn.name = "TOPIN"; topIn.injector = true; topIn.injBhpControl = true; topIn.bhpMax = pTop; topIn.zinj = {0.0, 1.0}; topIn.rate = 1.0;
    c.wells = {w1, w2, top, topIn};
    CompSimulator sim; std::string err;
    sim.setup(c);
    // saturated aqueous CO2 at the top conditions vs Carroll et al. (1991): 0.0384 mol/kg at 20 C and 1 atm CO2 partial pressure
    const co2::Fluid fl = co2::makeFluid(293.15, 0.0);
    const double mSat = co2::co2Molality(fl, pTop), mRef = 0.0384 * (pTop - steam::psat(293.15)) / 101325.0;
    std::printf("  CO2 solubility at 1.1 bar, 20 C: %.4f mol/kg (Henry / Carroll %.4f)\n", mSat, mRef);
    check(std::fabs(mSat / mRef - 1) < 0.06, "SPE11A conditions: CO2 solubility at 1.1 bar / 20 C within 6% of the Henry-law value");
    const double co2In0 = sim.componentInPlace(0);
    bool ok = true;
    auto t0 = std::chrono::steady_clock::now();
    auto& Wl = sim.wells();
    auto column = [&](const char* tag) { if (!std::getenv("SPE11_DEBUG")) return; std::printf("    p column (%s):", tag); for (int k = 0; k < nz; ++k) std::printf(" %.0f", sim.pressure()[size_t(id(nx / 2, k))]); std::printf("\n"); };
    column("t=0");
    ok = sim.advanceTo(60, 60); column("t=60 s");
    ok = sim.advanceTo(9000, 300); column("t=2.5 h");
    Wl[1].open = true;
    ok = ok && sim.advanceTo(18000, 300);
    double mob5, imm5, dis5; sim.co2Inventory(mob5, imm5, dis5);
    // gas location at the end of injection: below the lower seal top (z < 0.6 m) / upper reservoir / above
    auto gasSplit = [&](double& low, double& mid, double& high) {
        const std::vector<double> sg = sim.gasSaturation(); low = mid = high = 0;
        for (int k = 0; k < nz; ++k) for (int i = 0; i < nx; ++i) {
            const double v = sg[size_t(id(i, k))] * c.cellPV[size_t(id(i, k))], z = zUpOf(k);
            (z < 0.6 ? low : z < 1.1 ? mid : high) += v;
        }
    };
    double gl, gm, gh; gasSplit(gl, gm, gh); column("t=5 h");
    Wl[0].open = Wl[1].open = false;
    for (double t = 18000 + 3600; t <= 432000 + 1 && ok; t += 3600 * 6) { ok = sim.advanceTo(std::min(t, 432000.0), 1800); if (std::fmod(t - 21600, 86400) < 1) column("day"); }
    ok = ok && sim.advanceTo(432000, 1800);
    column("t=5 d");
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double mob, imm, dis; sim.co2Inventory(mob, imm, dis);
    if (std::getenv("SPE11_DEBUG")) {
        const std::vector<double> sg = sim.gasSaturation(), xd = sim.dissolvedCO2(); std::vector<double> X, Y, Z; sim.compositions(X, Y, Z);
        std::vector<int> ord(sg.size()); for (size_t q = 0; q < ord.size(); ++q) ord[q] = int(q);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return sg[size_t(a)] > sg[size_t(b)]; });
        for (int j = 0; j < 8; ++j) { const int q = ord[size_t(j)]; const double* N = &sim.moles()[size_t(q) * sim.ncw()];
            std::printf("    cell %d (i %d k %d) sg %.4f p %.1f Pa N %.3e %.3e xCO2aq %.3e y %.4f z %.3e ph %d\n", q, q % nx, q / nx, sg[size_t(q)], sim.pressure()[size_t(q)], N[0], N[1], xd[size_t(q)], Y[size_t(q) * 2], Z[size_t(q) * 2], sim.phasesIn(q)); }
    }
    const double inj = sim.injected(0), prod = sim.produced(0), inPlace = sim.componentInPlace(0);
    const double bal = (inPlace - co2In0 - inj + prod) / std::max(inj, 1e-30);
    double gl2, gm2, gh2; gasSplit(gl2, gm2, gh2);
    std::printf("  injected %.3f g (protocol %.3f g); 5 h: mobile %.3f immobile %.3f dissolved %.3f g; gas volume below / between / above seals %.1f / %.1f / %.1f cm3\n",
                inj * co2::MW_CO2 * 1e3, 1.7e-7 * 27000 * 1e3, mob5 * 1e3, imm5 * 1e3, dis5 * 1e3, gl * 1e6, gm * 1e6, gh * 1e6);
    std::printf("  5 d: mobile %.3f immobile %.3f dissolved %.3f g (%.0f %% dissolved), out of the top %.4f g; gas %.1f / %.1f / %.1f cm3; balance %.1e\n",
                mob * 1e3, imm * 1e3, dis * 1e3, 100 * dis / std::max(1e-30, mob + imm + dis), prod * co2::MW_CO2 * 1e3, gl2 * 1e6, gm2 * 1e6, gh2 * 1e6, bal);
    std::printf("        %d steps, %d Newton, %d cuts, %.1f s\n", sim.stats().steps, sim.stats().newton, sim.stats().cuts, secs);
    check(ok && std::fabs(inj * co2::MW_CO2 / (1.7e-7 * 27000) - 1) < 1e-3, "SPE11A-like: injection protocol delivered (4.59 g CO2)");
    check(ok && std::fabs(bal) < 1e-6, "SPE11A-like: CO2 mass balance closes");
    check(ok && gh < 1e-3 * (gl + gm + gh) && gl > 0 && gm > 0, "SPE11A-like: free CO2 held below the seals (capillary entry barrier)");
    check(ok && dis > dis5 && dis > 0.2 * (mob + imm + dis), "SPE11A-like: dissolution continues after shut-in (>20% dissolved at 5 d)");
}

// ------------------------------------------------------------------------------------------------ C16 verifications
// Point-distributed (gridpoint) control volumes of a triangular symmetry element: the square of side d around every
// gridpoint clipped to the element, face lengths clipped likewise (exact for straight symmetry lines).
namespace c16v {
struct Tri { double ax, ay, bx, by, cx, cy; };   // counter-clockwise
inline bool inTri(const Tri& t, double x, double y, double eps = 1e-9) {
    auto cr = [](double x1, double y1, double x2, double y2, double x, double y) { return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1); };
    return cr(t.ax, t.ay, t.bx, t.by, x, y) >= -eps && cr(t.bx, t.by, t.cx, t.cy, x, y) >= -eps && cr(t.cx, t.cy, t.ax, t.ay, x, y) >= -eps;
}
// area of the square [x0, x1] x [y0, y1] inside the triangle (Sutherland-Hodgman)
inline double clipArea(const Tri& t, double x0, double y0, double x1, double y1) {
    std::vector<std::pair<double, double>> P = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    const double V[4][2] = {{t.ax, t.ay}, {t.bx, t.by}, {t.cx, t.cy}, {t.ax, t.ay}};
    for (int e = 0; e < 3 && !P.empty(); ++e) {
        const double ex = V[e + 1][0] - V[e][0], ey = V[e + 1][1] - V[e][1];
        auto side = [&](const std::pair<double, double>& p) { return ex * (p.second - V[e][1]) - ey * (p.first - V[e][0]); };
        std::vector<std::pair<double, double>> Q;
        for (size_t i = 0; i < P.size(); ++i) {
            const auto& a = P[i]; const auto& b = P[(i + 1) % P.size()];
            const double sa = side(a), sb = side(b);
            if (sa >= 0) Q.push_back(a);
            if ((sa >= 0) != (sb >= 0)) { const double u = sa / (sa - sb); Q.push_back({a.first + u * (b.first - a.first), a.second + u * (b.second - a.second)}); }
        }
        P = Q;
    }
    double A = 0;
    for (size_t i = 0; i < P.size(); ++i) { const auto& a = P[i]; const auto& b = P[(i + 1) % P.size()]; A += a.first * b.second - b.first * a.second; }
    return 0.5 * std::fabs(A);
}
// length of the segment (x0,y0)-(x1,y1) inside the triangle (sampled; straight cuts -> 1e-4 relative)
inline double clipLen(const Tri& t, double x0, double y0, double x1, double y1) {
    const int n = 20000; int in = 0;
    for (int k = 0; k < n; ++k) { const double u = (k + 0.5) / n; if (inTri(t, x0 + u * (x1 - x0), y0 + u * (y1 - y0))) ++in; }
    return std::hypot(x1 - x0, y1 - y0) * in / n;
}
struct Element {
    std::vector<std::pair<int, int>> pts;   // gridpoint (i, j)
    std::vector<double> area;               // m2
    struct F { int a, b; double len; };     // len: face length (m)
    std::vector<F> faces;
    int index(int i, int j) const { for (size_t q = 0; q < pts.size(); ++q) if (pts[q].first == i && pts[q].second == j) return int(q); return -1; }
};
// ninePoint: Yanosik-McCracken (1979) nine-point scheme against the grid-orientation effect of adverse-mobility (steam)
// displacements - axial transmissibilities x 2/3, diagonal ones T/6 (consistent: T_a + 2 T_d = T); a diagonal that runs
// along a symmetry edge is shared with the mirror element (x 1/2). Face 'len' is the effective length (T = k h len / d).
inline Element buildElement(const Tri& t, int ni, int nj, double d, bool ninePoint = false) {
    Element E;
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni; ++i)
            if (inTri(t, i * d, j * d, 1e-6 * d)) { E.pts.push_back({i, j}); E.area.push_back(clipArea(t, (i - 0.5) * d, (j - 0.5) * d, (i + 0.5) * d, (j + 0.5) * d)); }
    for (size_t q = 0; q < E.pts.size(); ++q) {
        const int i = E.pts[q].first, j = E.pts[q].second;
        const int r = E.index(i + 1, j), u = E.index(i, j + 1);
        if (r >= 0) { const double L = clipLen(t, (i + 0.5) * d, (j - 0.5) * d, (i + 0.5) * d, (j + 0.5) * d); if (L > 1e-9) E.faces.push_back({int(q), r, L}); }
        if (u >= 0) { const double L = clipLen(t, (i - 0.5) * d, (j + 0.5) * d, (i + 0.5) * d, (j + 0.5) * d); if (L > 1e-9) E.faces.push_back({int(q), u, L}); }
    }
    if (ninePoint) {
        for (auto& f : E.faces) f.len *= 2.0 / 3.0;
        auto onEdge = [&](double x0, double y0, double x1, double y1) {   // both points on one triangle edge
            const double V[4][2] = {{t.ax, t.ay}, {t.bx, t.by}, {t.cx, t.cy}, {t.ax, t.ay}};
            for (int e = 0; e < 3; ++e) {
                const double ex = V[e + 1][0] - V[e][0], ey = V[e + 1][1] - V[e][1], L = std::hypot(ex, ey);
                if (std::fabs(ex * (y0 - V[e][1]) - ey * (x0 - V[e][0])) < 1e-6 * L * d && std::fabs(ex * (y1 - V[e][1]) - ey * (x1 - V[e][0])) < 1e-6 * L * d) return true;
            }
            return false;
        };
        for (size_t q = 0; q < E.pts.size(); ++q) {
            const int i = E.pts[q].first, j = E.pts[q].second;
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const int r = E.index(i + 1, j + sgn);
                if (r < 0) continue;
                const double xm = (i + 0.5) * d, ym = (j + 0.5 * sgn) * d;
                double fac = 0;
                if (onEdge(i * d, j * d, (i + 1) * d, (j + sgn) * d)) fac = 0.5;
                else if (inTri(t, xm, ym, -1e-6 * d)) fac = 1.0;
                if (fac > 0) E.faces.push_back({int(q), r, d * fac / 6.0});
            }
        }
    }
    return E;
}
}  // namespace c16v

// SPE4 problem 2A (Aziz, Ramesh, Woo 1987): steamflood of the problem-1 heavy oil in a 2.5-acre inverted nine-spot, 1/8
// symmetry element on the specified 9 x 5 x 4 gridpoint grid (diagonal: injector -> far producer = 8 d, d = 29.17 ft; near
// producer at (4, 4)). Injection 300 STB/d CWE (full well) into the bottom layer, 0.7 quality at 450 F, BHP <= 1000 psia;
// producers 17 psia, <= 1000 STB/d liquid, <= 10 STB/d steam (CWE, lagged BHP back-pressure control), 10 years. Reference:
// Table 1 of the paper, time for each producer to reach 10 STB/d of steam: near 257-539 d, far 1,398-2,245 d.
static void testSpe4Problem2A() {
    using namespace ressim;
    std::printf("SPE4 problem 2A: steamflood, inverted nine-spot 1/8 element, 9 x 5 x 4 (Aziz et al. 1987)\n");
    const double ft = 0.3048, md = 9.869233e-16, psi = 6894.757, stb = 0.158987294928, day = 86400.0;
    const double d = 330.0 * std::sqrt(0.5) / 8 * ft;   // 29.17 ft
    const c16v::Tri tri{0, 0, 8 * d, 0, 4 * d, 4 * d};
    const bool nine = !std::getenv("SPE4_5PT");
    const c16v::Element E = c16v::buildElement(tri, 9, 5, d, nine);
    std::printf("  %s-point flux stencil\n", nine ? "nine (Yanosik-McCracken)" : "five");
    const int np = int(E.pts.size()), nzl = 4;
    double aSum = 0; for (double a : E.area) aSum += a;
    std::printf("  %d gridpoints per layer, element area %.1f m2 (pattern / 8 = %.1f m2)\n", np, aSum, 2.5 * 4046.8564224 / 8);
    check(np == 25 && std::fabs(aSum / (2.5 * 4046.8564224 / 8) - 1) < 1e-3, "SPE4 2A: 1/8 element of 25 gridpoints, area = pattern / 8");
    const double dzF[4] = {10, 20, 25, 25}, khMd[4] = {2000, 500, 1000, 2000};
    const double top = 1500 * ft, phi = 0.3, pTop = 75 * psi, rhoMix = 972, g = 9.80665;
    std::vector<double> dz(nzl), zc(nzl); double z = 0;
    for (int k = 0; k < nzl; ++k) { dz[size_t(k)] = dzF[k] * ft; zc[size_t(k)] = z + 0.5 * dz[size_t(k)]; z += dz[size_t(k)]; }
    // properties of problem 1
    std::vector<int> dummy; thermal::Case c = spe4Problem1(dummy);
    c.conns.clear(); c.lossFaces.clear(); c.wells.clear();
    auto id = [&](int q, int k) { return k * np + q; };
    const int n = np * nzl;
    c.cellPV.assign(size_t(n), 0); c.cellBulk.assign(size_t(n), 0); c.cellDepth.assign(size_t(n), 0); c.p0.assign(size_t(n), 0);
    for (int k = 0; k < nzl; ++k)
        for (int q = 0; q < np; ++q) {
            const double V = E.area[size_t(q)] * dz[size_t(k)];
            c.cellBulk[size_t(id(q, k))] = V; c.cellPV[size_t(id(q, k))] = phi * V; c.cellDepth[size_t(id(q, k))] = top + zc[size_t(k)];
            c.p0[size_t(id(q, k))] = pTop + rhoMix * g * (zc[size_t(k)] - zc[0]);
        }
    for (int k = 0; k < nzl; ++k) {
        for (const auto& f : E.faces) {
            const double h = dz[size_t(k)], kh = khMd[k] * md;
            c.conns.push_back({id(f.a, k), id(f.b, k), kh * h * f.len / d, h * f.len / d});
        }
        if (k + 1 < nzl)
            for (int q = 0; q < np; ++q) {
                const double A = E.area[size_t(q)], kv0 = 0.5 * khMd[k] * md, kv1 = 0.5 * khMd[k + 1] * md;
                c.conns.push_back({id(q, k), id(q, k + 1), A / (0.5 * dz[size_t(k)] / kv0 + 0.5 * dz[size_t(k + 1)] / kv1), A / (0.5 * dz[size_t(k)] + 0.5 * dz[size_t(k + 1)])});
            }
    }
    for (int q = 0; q < np; ++q) { c.lossFaces.push_back({id(q, 0), E.area[size_t(q)]}); c.lossFaces.push_back({id(q, nzl - 1), E.area[size_t(q)]}); }
    // wells: Peaceman r_e = 0.2 d (Elf's choice in the paper), well fraction theta = interior angle / 360
    const double rw = 0.3 * ft, thI = 1.0 / 8, thN = 1.0 / 4, thF = 1.0 / 8;
    auto wi = [&](double th, int k) { return th * 2 * M_PI * khMd[k] * md * dz[size_t(k)] / std::log(0.2 * d / rw); };
    const int qI = E.index(0, 0), qN = E.index(4, 4), qF = E.index(8, 0);
    thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.zinj = {0, 1}; inj.quality = 0.7; inj.tinj = (450 - 32) / 1.8 + 273.15;
    inj.rate = thI * 300 * stb / day * 1000 / steam::kMw; inj.bhpMax = 1000 * psi; inj.refDepth = top + zc[3];
    inj.cells = {id(qI, 3)}; inj.wi = {wi(thI, 3)};
    auto prod = [&](const char* nm, int qq, double th) {
        thermal::Well p; p.name = nm; p.control = 3; p.target = th * 1000 * stb / day; p.bhp = 17 * psi; p.refDepth = top + zc[0];
        for (int k = 0; k < nzl; ++k) { p.cells.push_back(id(qq, k)); p.wi.push_back(wi(th, k)); }
        return p;
    };
    c.wells = {inj, prod("NEAR", qN, thN), prod("FAR", qF, thF)};
    thermal::Simulator s; std::string err;
    if (!s.setup(c, &err)) { std::printf("  setup failed: %s\n", err.c_str()); check(false, "SPE4 2A setup"); return; }
    const double ooip = s.stockTankOilInPlace();
    std::printf("  OOIP %.0f STB (element), %.0f STB full pattern\n", ooip / stb, 8 * ooip / stb);
    auto& W = s.wells();
    const double steamMax[3] = {0, thN * 10 * stb / day, thF * 10 * stb / day}, th[3] = {thI, thN, thF};
    double tHit[3] = {-1, -1, -1};
    int nCtl[3] = {0, 0, 0}, nIn[3] = {0, 0, 0}, nLim[3] = {0, 0, 0};
    const bool lagged = std::getenv("SPE4_LAGGED") != nullptr;   // C17 default: implicit steam-rate limit in the well model
    if (!lagged) { W[1].steamMax = steamMax[1]; W[2].steamMax = steamMax[2]; }   // reports after breakthrough / with steam below 1.5 x the limit
    bool ok = true;
    auto t0 = std::chrono::steady_clock::now();
    const double dtRep = 5 * day;
    for (double t = dtRep; t <= 3650 * day + 1 && ok; t += dtRep) {
        ok = s.advanceTo(t, dtRep);
        if (!ok) break;
        const auto& R = s.wellReports();
        for (int w = 1; w < 3; ++w) {
            const double qs = R[size_t(w)].qsteam;
            if (tHit[w] < 0 && qs >= steamMax[w] * 0.999) tHit[w] = t / day;
            if (lagged) {   // old C16 variant: lagged back-pressure control on the 5-day reports
                if (qs > steamMax[w]) W[size_t(w)].bhp = std::min(W[size_t(w)].bhp * std::min(1.25, 1 + 0.5 * (qs / steamMax[w] - 1)), 900 * psi);
                else if (qs < 0.7 * steamMax[w]) W[size_t(w)].bhp = std::max(17 * psi, W[size_t(w)].bhp / (qs < 0.2 * steamMax[w] ? 1.2 : 1.05));
            }
            if (tHit[w] > 0) { ++nCtl[w]; if (qs < 1.02 * steamMax[w]) ++nIn[w]; if (R[size_t(w)].steamLimited) ++nLim[w]; }
        }
        if (std::fmod(t / day, 500) < 1)
            std::printf("  t %5.0f d: oil near %.1f far %.1f STB/d, steam near %.2f far %.2f STB/d (full well), BHP near %.0f psia, Np %.1f MSTB pattern\n",
                        t / day, R[1].qo / thN / stb * day, R[2].qo / thF / stb * day, R[1].qsteam / thN / stb * day, R[2].qsteam / thF / stb * day,
                        R[1].bhp / psi, 8 * (R[1].cumO + R[2].cumO) / stb / 1000);
    }
    (void)th;
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto& R = s.wellReports();
    std::printf("  time to 10 STB/d steam: near %.0f d (paper 257-539), far %.0f d (paper 1,398-2,245)\n", tHit[1], tHit[2]);
    std::printf("  %s steam limit: steam <= 1.02 x limit in %d / %d (near), %d / %d (far) reports after breakthrough; limit binding in %d / %d reports\n",
                lagged ? "lagged" : "implicit", nIn[1], nCtl[1], nIn[2], nCtl[2], nLim[1], nLim[2]);
    std::printf("  10 y: Np %.1f MSTB full pattern, injected %.1f MSTB CWE, cum OSR %.3f, heat loss %.0f MMBtu; %d steps %d Newton %d cuts, %.1f s\n",
                8 * (R[1].cumO + R[2].cumO) / stb / 1000, 8 * R[0].cumCwe / stb / 1000, (R[1].cumO + R[2].cumO) / std::max(1e-9, R[0].cumCwe),
                8 * s.cumHeatLoss() / 1.055056e9, s.stats().steps, s.stats().newton, s.stats().cuts, secs);
    check(ok, "SPE4 2A runs 10 years");
    check(ok && tHit[1] > 257 * 0.85 && tHit[1] < 539 * 1.15, "SPE4 2A: near producer reaches 10 STB/d steam within the participants' 257-539 d (+-15%)");
    check(ok && tHit[2] > 1398 * 0.85 && tHit[2] < 2245 * 1.15, "SPE4 2A: far producer reaches 10 STB/d steam within the participants' 1,398-2,245 d (+-15%)");
    if (!lagged) check(ok && nIn[1] >= nCtl[1] * 95 / 100 && nIn[2] >= nCtl[2] * 95 / 100, "SPE4 2A: implicit steam-rate limit holds (steam <= 1.02 x 10 STB/d in >= 95% of the reports after breakthrough)");
}

// SPE11A on the official geometry (spe11a.geo of the 11thSPE-CSP repository, rasterised by tools/spe11_raster.py to
// decks/spe11/spe11a_facies_NXxNZ.txt): all seven facies of Table 1 (facies 7 inactive), Brooks-Corey kr / Pc of Table 2,
// both wells and the injection protocol, top boundary at 1.1e5 Pa in contact with water; molecular diffusion (no
// dispersion). Sparse SPE11 quantities: Box A / B mobile, immobile, dissolved and seal CO2 at 72 h, the time of maximum
// mobile CO2 in Box A, seal CO2. Indicative comparison: the FluidFlower experiment and two forecasts (Landa-Marban et al.
// 2025, Table 4) - same geometry and protocol, but the physical rig differs from the SPE11A specification (thickness,
// sand properties), so only order-of-magnitude checks are made.
static void testSpe11aGeometry() {
    using namespace ressim; using namespace ressim::comp;
    const int nx = std::getenv("SPE11_NX") ? std::atoi(std::getenv("SPE11_NX")) : 56, nz = nx * 3 / 7;
    std::printf("SPE11A on the official geometry (%d x %d cells), dispersivity %.3g m\n", nx, nz, std::getenv("SPE11_DISP") ? std::atof(std::getenv("SPE11_DISP")) : 1e-2);
    char path[160]; std::snprintf(path, sizeof path, "decks/spe11/spe11a_facies_%dx%d.txt", nx, nz);
    std::ifstream in(path);
    if (!in) { std::printf("  %s not found (run from the ressim folder; tools/spe11_raster.py makes other resolutions)\n", path); check(false, "SPE11A geometry: facies map present"); return; }
    std::vector<int> fac; { std::string line; while (std::getline(in, line)) { if (line.empty() || line[0] == '#') continue; std::istringstream ss(line); int v; while (ss >> v) fac.push_back(v); } }
    if (int(fac.size()) != nx * nz) { check(false, "SPE11A geometry: facies map size"); return; }
    const double Lx = 2.8, Lz = 1.2, W = 0.01, dx = Lx / nx, dz = Lz / nz, g = 9.80665, pTop = 1.1e5;
    struct Fac { double k, phi, swi, pe; };
    const Fac F[7] = {{4e-11, 0.44, 0.32, 1500}, {5e-10, 0.43, 0.14, 300}, {1e-9, 0.44, 0.12, 100}, {2e-9, 0.45, 0.12, 25},
                      {4e-9, 0.43, 0.12, 10}, {1e-8, 0.46, 0.10, 1}, {0, 0, 0, 0}};
    CompCase c;
    c.fluid.comp = {eos::standardComponent("CO2"), eos::standardComponent("H2O")};
    c.co2store = true; c.salinity = 0; c.T = 293.15; c.cr = 0; c.pref = pTop;
    c.trap.diffusion = true; c.trap.dCo2Water = 1e-9; c.trap.dH2oWater = 1e-9; c.trap.dCo2Gas = 1.6e-5; c.trap.dH2oGas = 1.6e-5;
    const double alphaD = std::getenv("SPE11_DISP") ? std::atof(std::getenv("SPE11_DISP")) : 1e-2;   // m (C17: dispersion on)
    auto table = [&](const Fac& f) {
        SatTables t; const double sgi = 0.1;
        for (int j = 0; j <= 30; ++j) {
            const double sw = f.swi + (1 - f.swi) * j / 30.0, swn = (sw - f.swi) / (1 - f.swi);
            t.sw.push_back(sw); t.krw.push_back(swn * swn);
            t.pcow.push_back(std::min(9.5e4, swn > 0 ? f.pe / std::sqrt(swn) : 9.5e4));
        }
        for (int j = 0; j <= 30; ++j) {
            const double sg = (1 - f.swi) * j / 30.0, sgn = std::max(0.0, (sg - sgi) / (1 - sgi - f.swi));
            t.sg.push_back(sg); t.krg.push_back(sgn * sgn); t.pcog.push_back(0);
        }
        return t;
    };
    c.sat = table(F[0]);
    for (int r = 1; r < 6; ++r) c.satRegions.push_back(table(F[r]));
    std::vector<int> cellOf(size_t(nx * nz), -1), facOfCell;
    std::vector<double> xc, zc;   // cell centres (x, z up)
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i) {
            const int f = fac[size_t(k * nx + i)];
            if (f == 7) continue;
            const Fac& P = F[f - 1];
            cellOf[size_t(k * nx + i)] = int(c.cellPV.size());
            c.cellPV.push_back(P.phi * dx * dz * W); c.cellDepth.push_back((k + 0.5) * dz);
            c.satnum.push_back(f - 1); c.p0.push_back(pTop + 998.2 * g * (k + 0.5) * dz);
            facOfCell.push_back(f); xc.push_back((i + 0.5) * dx); zc.push_back(Lz - (k + 0.5) * dz);
        }
    const int n = int(c.cellPV.size());
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i) {
            const int a = cellOf[size_t(k * nx + i)];
            if (a < 0) continue;
            const Fac& A = F[fac[size_t(k * nx + i)] - 1];
            auto link = [&](int b, int fb, double area, double len) {
                const Fac& B = F[fb - 1];
                const double T = area / (0.5 * len / A.k + 0.5 * len / B.k), phiH = 2 * A.phi * B.phi / (A.phi + B.phi);
                c.conns.push_back({a, b, T, phiH * area / len, alphaD / len});   // C17: SPE11A dispersivity 1e-2 m
            };
            if (i + 1 < nx && cellOf[size_t(k * nx + i + 1)] >= 0) link(cellOf[size_t(k * nx + i + 1)], fac[size_t(k * nx + i + 1)], dz * W, dx);
            if (k + 1 < nz && cellOf[size_t((k + 1) * nx + i)] >= 0) link(cellOf[size_t((k + 1) * nx + i)], fac[size_t((k + 1) * nx + i)], dx * W, dz);
        }
    c.z0 = {{1e-9, 1 - 1e-9}};
    auto cellAt = [&](double x, double zUp) { return cellOf[size_t(std::min(nz - 1, int((Lz - zUp) / dz)) * nx + std::min(nx - 1, int(x / dx)))]; };
    const double rate = 1.7e-7 / co2::MW_CO2;
    auto mkInj = [&](const char* nm, int cell) {
        CompWell w; w.name = nm; w.injector = true; w.cells = {cell};
        w.wi = {2 * M_PI * F[facOfCell[size_t(cell)] - 1].k * W / std::log(0.2 * dx / 9e-4)};
        w.zinj = {1.0, 0.0}; w.rate = rate; w.bhpMax = 5e5; return w;
    };
    const int c1 = cellAt(0.9, 0.3), c2 = cellAt(1.7, 0.7);
    if (c1 < 0 || c2 < 0) { check(false, "SPE11A geometry: wells in active cells"); return; }
    std::printf("  %d active cells; well 1 in facies %d, well 2 in facies %d\n", n, facOfCell[size_t(c1)], facOfCell[size_t(c2)]);
    CompWell w1 = mkInj("W1", c1), w2 = mkInj("W2", c2); w2.open = false;
    CompWell topW; topW.name = "TOP"; topW.bhp = pTop; topW.refDepth = 0.5 * dz;
    for (int i = 0; i < nx; ++i) { const int q = cellOf[size_t(i)]; if (q < 0) continue; topW.cells.push_back(q); topW.wi.push_back(F[facOfCell[size_t(q)] - 1].k * dx * W / (0.5 * dz)); }
    CompWell topIn = topW; topIn.name = "TOPIN"; topIn.injector = true; topIn.injBhpControl = true; topIn.bhpMax = pTop; topIn.zinj = {0.0, 1.0}; topIn.rate = 1.0;
    c.wells = {w1, w2, topW, topIn};
    CompSimulator sim; sim.setup(c);
    const double co2In0 = sim.componentInPlace(0);
    // box inventories (g): mobile (s_g > s_n,imm = 0.1), immobile, dissolved, seal (all CO2 in facies 1)
    struct Box { double x0, z0, x1, z1; };
    const Box bA{1.1, 0.0, 2.8, 0.6}, bB{0.0, 0.6, 1.1, 1.2};
    auto inv = [&](const Box& b, double& mob, double& imm, double& dis, double& seal) {
        mob = imm = dis = seal = 0;
        const std::vector<double> sg = sim.gasSaturation(), xd = sim.dissolvedCO2();
        const int ncw = sim.ncw();
        for (int q = 0; q < n; ++q) {
            if (xc[size_t(q)] < b.x0 || xc[size_t(q)] > b.x1 || zc[size_t(q)] < b.z0 || zc[size_t(q)] > b.z1) continue;
            const double* N = &sim.moles()[size_t(q) * size_t(ncw)];
            const double x = std::clamp(xd[size_t(q)], 0.0, 0.5), nd = std::min(N[0], x / (1 - x) * N[1]), ngas = std::max(0.0, N[0] - nd);
            const double gm = co2::MW_CO2 * 1e3;
            dis += nd * gm;
            (sg[size_t(q)] > 0.1 ? mob : imm) += ngas * gm;
            if (facOfCell[size_t(q)] == 1) seal += N[0] * gm;
        }
    };
    auto pAt = [&](double x, double zUp) { const int q = cellAt(x, zUp); return q >= 0 ? sim.pressure()[size_t(q)] : 0.0; };
    auto& Wl = sim.wells();
    bool ok = true;
    double mobAmax = -1, tMobAmax = 0, p1max = 0, p2max = 0;
    auto t0 = std::chrono::steady_clock::now();
    double A72[4] = {0, 0, 0, 0}, B72[4] = {0, 0, 0, 0};
    double t = 0;
    while (t < 432000 - 1 && ok) {
        const double dtOut = t < 86400 ? 600 : 3600;
        double tn = std::min(432000.0, t + dtOut);
        if (t < 9000 && tn > 9000) tn = 9000;
        if (t < 18000 && tn > 18000) tn = 18000;
        if (t < 259200 && tn > 259200) tn = 259200;
        ok = sim.advanceTo(tn, std::min(dtOut, 600.0));
        t = tn;
        if (std::fabs(t - 9000) < 1) Wl[1].open = true;
        if (std::fabs(t - 18000) < 1) { Wl[0].open = Wl[1].open = false; }
        double m, i_, d_, s_; inv(bA, m, i_, d_, s_);
        if (m > mobAmax) { mobAmax = m; tMobAmax = t; }
        p1max = std::max(p1max, pAt(1.5, 0.5)); p2max = std::max(p2max, pAt(1.7, 1.1));
        if (std::fabs(t - 259200) < 1) { inv(bA, A72[0], A72[1], A72[2], A72[3]); inv(bB, B72[0], B72[1], B72[2], B72[3]); }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double inj = sim.injected(0), prod = sim.produced(0), inPlace = sim.componentInPlace(0);
    const double bal = (inPlace - co2In0 - inj + prod) / std::max(inj, 1e-30);
    double sealAll = 0; { const int ncw = sim.ncw(); for (int q = 0; q < n; ++q) if (facOfCell[size_t(q)] == 1) sealAll += sim.moles()[size_t(q) * size_t(ncw)] * co2::MW_CO2 * 1e3; }
    double mob, imm, dis; sim.co2Inventory(mob, imm, dis);
    std::printf("  injected %.3f g (protocol 4.590), mass balance %.1e, out of the top %.4f g; max POP1 %.0f Pa, POP2 %.0f Pa\n", inj * co2::MW_CO2 * 1e3, bal, prod * co2::MW_CO2 * 1e3, p1max, p2max);
    std::printf("  max mobile CO2 in Box A at %.0f s (FluidFlower 14,800; forecasts 17,700-18,000)\n", tMobAmax);
    std::printf("  72 h Box A: mobile %.3f (exp 0.23, fc 0.35-0.87) immobile %.3f dissolved %.3f (exp 3.10, fc 2.11-3.60) seal %.3f g (exp 0.38, fc 0.43-0.57)\n", A72[0], A72[1], A72[2], A72[3]);
    std::printf("  72 h Box B: mobile %.3f immobile %.3f dissolved %.3f (exp 0.78, fc 0.30-0.52) seal %.3f g\n", B72[0], B72[1], B72[2], B72[3]);
    std::printf("  5 d: CO2 in seal facies %.3f g; total mobile %.3f immobile %.3f dissolved %.3f g; %d steps %d Newton %d cuts, %.1f s\n",
                sealAll, mob * 1e3, imm * 1e3, dis * 1e3, sim.stats().steps, sim.stats().newton, sim.stats().cuts, secs);
    check(ok && std::fabs(inj * co2::MW_CO2 / (1.7e-7 * 27000) - 1) < 1e-3, "SPE11A geometry: injection protocol delivered (4.59 g CO2)");
    check(ok && std::fabs(bal) < 1e-6, "SPE11A geometry: CO2 mass balance closes");
    check(ok && tMobAmax > 1.0e4 && tMobAmax < 2.6e4, "SPE11A geometry: mobile CO2 in Box A peaks at the end of injection (1.0e4-2.6e4 s; FluidFlower 1.48e4)");
    check(ok && A72[2] > 1.0 && A72[2] < 6.2, "SPE11A geometry: dissolved CO2 in Box A at 72 h of the experiment's order (1.0-6.2 g; FluidFlower 3.10 g)");
    check(ok && A72[3] < 1.5 && B72[0] + B72[1] + B72[2] > 0, "SPE11A geometry: seal holds (Box A seal < 1.5 g) and CO2 reaches Box B through the fault");
}

// EGS paper layout: Pruess (2006, Geothermics 35; LBNL-60397) five-spot fractured reservoir - 1 km2 pattern, 305 m
// thick, 1/8 symmetry element; MINC (fracture spacing 50 m, three orthogonal sets, 2 % permeable volume with 50 %
// porosity, 50e-15 m2), rock 2650 kg/m3, 1000 J/kg/K, 2.1 W/m/K; 200 C / 500 bar; injection 20 C at 510 bar downhole,
// production at 490 bar. Water and (optionally) supercritical CO2 as the heat-transmission fluid. Paper: CO2 gives ~ 75 MW
// and ~ 280 kg/s (full well), ~ 50 % more heat than water and ~ 3.7 x the initial mass flow -> water ~ 50 MW.
// The MINC continua are built here with general connections (proximity function of cubic blocks, five matrix shells).
static bool pruessFiveSpot(bool co2, double years, std::vector<double>& tYr, std::vector<double>& heatMW, std::vector<double>& massKg, double& secs) {
    using namespace ressim;
    const double L = 500, H = 305, D = 50, kf = 50e-15, yr = 365.25 * 86400;
    const double rw = std::getenv("EGS_RW") ? std::atof(std::getenv("EGS_RW")) : 0.1;   // well radius (not given in the paper)
    const int ng = std::getenv("EGS_N") ? std::atoi(std::getenv("EGS_N")) : 20;
    const double d = L / ng;
    const c16v::Tri tri{0, 0, L, 0, L, L};
    const c16v::Element E = c16v::buildElement(tri, ng + 1, ng + 1, d);
    const int np = int(E.pts.size());
    // MINC: cumulative volume fractions measured from the fracture faces, V(x) / V0 = 1 - (1 - 2x/D)^3
    const std::vector<double> cum = {0.02, 0.06, 0.14, 0.30, 0.60, 1.0};
    const int ns = int(cum.size());   // continuum 0 = fracture, 1.. = matrix shells
    auto xOf = [&](double v) { return 0.5 * D * (1 - std::cbrt(1 - std::min(v, 1.0))); };
    auto areaPerVol = [&](double x) { const double s = 1 - 2 * x / D; return 6.0 / D * s * s; };
    std::vector<double> xNode(static_cast<size_t>(ns)), xIf(static_cast<size_t>(ns));
    for (int m = 0; m < ns; ++m) { const double lo = m ? cum[size_t(m - 1)] : 0.0; xIf[size_t(m)] = xOf(cum[size_t(m)]); xNode[size_t(m)] = m ? xOf(0.5 * (lo + cum[size_t(m)])) : 0.0; }
    thermal::Case c; c.gravity = false; c.heatLoss = false;
    const int n = np * ns;
    auto id = [&](int q, int m) { return m * np + q; };
    c.cellPV.assign(size_t(n), 0); c.cellBulk.assign(size_t(n), 0); c.cellDepth.assign(size_t(n), 3000); c.p0.assign(size_t(n), 500e5);
    const double phiF = 0.5, phiM = 1e-3, km = 1e-18;
    for (int q = 0; q < np; ++q)
        for (int m = 0; m < ns; ++m) {
            const double V = E.area[size_t(q)] * H * (cum[size_t(m)] - (m ? cum[size_t(m - 1)] : 0.0));
            c.cellBulk[size_t(id(q, m))] = V; c.cellPV[size_t(id(q, m))] = (m ? phiM : phiF) * V;
        }
    for (const auto& f : E.faces) c.conns.push_back({id(f.a, 0), id(f.b, 0), kf * H * f.len / d, 0.02 * H * f.len / d});   // global flow in the fracture continuum
    for (int q = 0; q < np; ++q)
        for (int m = 0; m + 1 < ns; ++m) {
            const double A = areaPerVol(xIf[size_t(m)]) * E.area[size_t(q)] * H, dist = std::max(1e-3, xNode[size_t(m + 1)] - (m ? xNode[size_t(m)] : xIf[0]));
            c.conns.push_back({id(q, m), id(q, m + 1), km * A / dist, A / dist});
        }
    c.thconr = {2.1}; c.thconsf = 0; c.rockCv = 2650 * 1000; c.cr = 1e-10; c.pref = 500e5;
    c.T0 = {473.15}; c.tRef = 293.15; c.stdT = 288.71;
    const double thI = 1.0 / 8, thP = 1.0 / 8;
    const int qI = E.index(0, 0), qP = E.index(ng, ng);
    auto wi = [&](double th) { return th * 2 * M_PI * kf * H / std::log(0.2 * d / rw); };
    thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.bhpControl = true; inj.bhpMax = 510e5; inj.quality = -1; inj.tinj = 293.15;
    inj.cells = {id(qI, 0)}; inj.wi = {wi(thI)};
    thermal::Well pr; pr.name = "PROD"; pr.control = 0; pr.bhp = 490e5; pr.cells = {id(qP, 0)}; pr.wi = {wi(thP)};
    if (!co2) {
        c.sw0 = {1.0}; c.sg0 = {0.0}; c.swc = 0; c.sorw = 0; c.krw0 = 1; c.nw = 1;
        thermal::Component o; o.name = "OIL"; o.dead = true; c.comps = {o}; c.x0 = {{1.0}};
        inj.zinj = {0, 1};
    } else {
        thermal::Component g; g.name = "CO2"; g.type = 1; g.mw = 0.04401; g.tcG = 304.13; g.pcG = 7.3773e6; g.omegaG = 0.22394;
        g.cpg1 = 37.1; g.cpg2 = 0.045; g.avg = 2e-5; c.gasPR = true;
        c.comps = {g}; c.sw0 = {0.0}; c.sg0 = {1.0}; c.x0 = {{1.0}}; c.y0 = {{1.0}};
        c.swc = 0; c.sgc = 0; c.sorg = 0; c.krg0 = 1; c.ng = 1;
        inj.zinj = {1, 0};
    }
    c.wells = {inj, pr};
    thermal::Simulator s; std::string err;
    if (!s.setup(c, &err)) { std::printf("  %s setup failed: %s\n", co2 ? "CO2" : "water", err.c_str()); return false; }
    auto t0 = std::chrono::steady_clock::now();
    bool ok = true;
    for (double t = 0.25 * yr; t <= years * yr + 1 && ok; t += 0.25 * yr) {
        ok = s.advanceTo(t, 0.05 * yr);
        if (!ok) break;
        const auto& R = s.wellReports();
        tYr.push_back(t / yr);
        heatMW.push_back((R[0].hRate + R[1].hRate) / thP * 1e-6);   // net heat extraction, full well
        double m = 0; for (size_t k = 0; k < R[1].molar.size(); ++k) m += std::max(0.0, R[1].molar[k]) * (k + 1 < R[1].molar.size() ? c.comps[k].mw : steam::kMw);
        massKg.push_back(m / thP);
    }
    secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return ok;
}
static void testPruessEgs() {
    std::printf("EGS five-spot (Pruess 2006): MINC fractured reservoir, water vs CO2 heat transmission fluid\n");
    std::vector<double> tw, hw, mw, tc, hc, mc; double sw = 0, sc = 0;
    const bool okW = pruessFiveSpot(false, 25, tw, hw, mw, sw);
    std::printf("  water: %s, %.1f s\n", okW ? "ok" : "FAILED", sw);
    for (size_t k = 0; k < tw.size(); k += 8) std::printf("    t %5.2f yr: net heat %.1f MW, mass %.1f kg/s\n", tw[k], hw[k], mw[k]);
    check(okW, "Pruess EGS water case runs 25 years");
    double avgW = 0; for (double h : hw) avgW += h; avgW /= std::max<size_t>(1, hw.size());
    std::printf("  water: initial %.1f MW / %.1f kg/s, 25-yr mean %.1f MW (paper: CO2 ~ 75 MW = 1.5 x water -> ~ 50 MW)\n", hw.empty() ? 0 : hw[0], mw.empty() ? 0 : mw[0], avgW);
    // absolute rates hinge on the well representation, which the paper does not give (Peaceman, rw = 0.1 m here): order only
    check(okW && !hw.empty() && hw[0] > 10 && hw[0] < 100, "Pruess EGS: water net heat extraction of the paper's order (10-100 MW, full well)");
    check(okW && hw.size() > 4 && hw.back() < hw[0], "Pruess EGS: heat extraction declines as the cooled region reaches the producer");
    if (!std::getenv("EGS_WATER_ONLY")) {
        const bool okC = pruessFiveSpot(true, 25, tc, hc, mc, sc);
        std::printf("  CO2: %s, %.1f s\n", okC ? "ok" : "FAILED", sc);
        for (size_t k = 0; k < tc.size(); k += 8) std::printf("    t %5.2f yr: net heat %.1f MW, mass %.1f kg/s\n", tc[k], hc[k], mc[k]);
        double avgC = 0; for (double h : hc) avgC += h; avgC /= std::max<size_t>(1, hc.size());
        if (okC && !hc.empty() && !hw.empty())
            std::printf("  CO2 / water: initial heat %.2f, 25-yr mean heat %.2f (paper ~1.5), initial mass flow %.2f (paper ~3.7); CO2 initial %.1f MW / %.1f kg/s (paper ~75 MW / ~280 kg/s)\n",
                        hc[0] / hw[0], avgC / std::max(1e-9, avgW), mc[0] / std::max(1e-9, mw[0]), hc[0], mc[0]);
        check(okC && !hc.empty() && hc[0] / hw[0] > 1.2 && hc[0] / hw[0] < 1.9, "Pruess EGS: CO2 extracts ~1.5 x the heat of water initially (1.2-1.9)");
        check(okC && !mc.empty() && mc[0] / mw[0] > 2.8 && mc[0] / mw[0] < 4.6, "Pruess EGS: CO2 mass flow ~3.7 x water initially (2.8-4.6)");
    }
}

// C13 verification: Gringarten single-fracture heat extraction (EDFM + conduction)
static void testGringarten() {
    using namespace ressim;
    std::printf("== C13 verification: Gringarten single-fracture heat extraction (EDFM)\n");
    {   // C13 verification: Gringarten (1975) single-fracture heat extraction - water through one vertical fracture between
        // impermeable conductive half-spaces: T_D = (T_r - T_out) / (T_r - T_inj) = erfc(lam A / (mdot c_w sqrt(kappa t))),
        // A = L H (one side). Fine (geometric) matrix cells next to the plane; a uniform coarse grid is reported for contrast.
        auto run = [&](bool fine, std::vector<double>& tdSim, std::vector<double>& tdAn, double& secs) {
            const int nx = std::getenv("GR_NX") ? std::atoi(std::getenv("GR_NX")) : 20; const double dx = 200.0 / nx, H = 20, lam = 2.5;
            std::vector<double> w;
            if (fine) {   // geometric matrix cells away from the plane (ratio GR_R, first GR_C m), to ~95 m each side
                const double r = std::getenv("GR_R") ? std::atof(std::getenv("GR_R")) : 1.5, c0 = std::getenv("GR_C") ? std::atof(std::getenv("GR_C")) : 0.25;
                std::vector<double> side; double d = c0, tot = 0;
                while (tot < 60) { side.push_back(d); tot += d; d *= r; }
                for (size_t k = side.size(); k-- > 0;) w.push_back(side[k]);
                w.push_back(c0);
                for (double v : side) w.push_back(v);
            }
            else w.assign(7, 25.0);
            const int ny = int(w.size());
            thermal::Case c; c.nx = nx; c.ny = ny; c.nz = 1; c.gravity = false; c.cr = 1e-10;
            const size_t n = size_t(nx) * size_t(ny);
            c.dx = {dx}; c.dz = {H}; c.dy.assign(n, 0.0);
            for (size_t q = 0; q < n; ++q) c.dy[q] = w[q / size_t(nx)];
            c.poro = {0.01}; c.permx = {1e-18}; c.permz = {1e-18};   // 1 uD: no advective heat, but cooled-water contraction drains
            c.krw0 = 1; c.swc = 0; c.sorw = 0; c.p0 = {30e6}; c.T0 = {470}; c.sw0 = {1.0}; c.thconr = {lam}; c.thconsf = 0;
            double yc = 0; for (int j = 0; j < ny / 2; ++j) yc += w[size_t(j)];
            yc += 0.5 * w[size_t(ny / 2)] + 1e-4;
            thermal::Fracture f; f.x0 = 0; f.x1 = nx * dx; f.y0 = f.y1 = yc; f.aperture = 1e-3;
            std::string err;
            bool ok = thermal::addFractures(c, {f}, thermal::FractureStress(), &err);
            const double rate = 27.75;   // mol/s (0.5 kg/s)
            thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {ok ? c.fracCells.front().cell : 0}; inj.wi = {1e-10}; inj.zinj = {1.0};
            inj.rate = rate; inj.bhpMax = 40e6; inj.quality = -1; inj.tinj = 320;
            thermal::Well pr; pr.name = "PROD"; pr.cells = {ok ? c.fracCells.back().cell : 0}; pr.wi = {1e-10}; pr.bhp = 29.9e6;
            c.wells = {inj, pr};
            thermal::Simulator s;
            ok = ok && s.setup(c, &err);
            const steam::Props Wl = steam::liquid(395, 30e6);
            const double mcp = rate * steam::kMw * (steam::liquid(470, 30e6).h - steam::liquid(320, 30e6).h) / 150.0;
            const double rc = 0.01 * Wl.cp / Wl.v + 0.99 * s.model().rockCv, kap = lam / rc, A = nx * dx * H;
            const auto t0 = std::chrono::steady_clock::now();
            for (double yr : {0.25, 0.5, 1.0, 2.0, 5.0}) {
                ok = ok && s.advanceTo(yr * 365.25 * 86400.0, (std::getenv("GR_DT") ? std::atof(std::getenv("GR_DT")) : 10.0) * 86400.0);
                if (!ok) break;
                const double tout = s.wellReports()[1].tFluid;
                tdSim.push_back((470 - tout) / 150.0);
                tdAn.push_back(std::erfc(lam * A / (mcp * std::sqrt(kap * yr * 365.25 * 86400.0))));
            }
            secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (std::getenv("GR_DEBUG")) std::printf("    %zu cells: steps %d newton %d linear %d cuts %d\n", s.temperature().size(), s.stats().steps, s.stats().newton, s.stats().linear, s.stats().cuts);
            if (!ok) std::printf("  error: %s\n", err.c_str());
            return ok;
        };
        std::vector<double> sF, aF, sC, aC; double tF = 0, tC = 0;
        const bool okF = run(true, sF, aF, tF), okC = run(false, sC, aC, tC);
        double eF = 0, eC = 0;
        std::printf("  Gringarten single fracture, T_D at 0.25 / 0.5 / 1 / 2 / 5 years:\n    analytic  ");
        for (double v : aF) std::printf(" %.3f", v);
        std::printf("\n    EDFM fine ");
        for (size_t i = 0; i < sF.size(); ++i) { std::printf(" %.3f", sF[i]); eF = std::max(eF, std::fabs(sF[i] - aF[i])); }
        std::printf("   (max |dT_D| %.3f, %.1fs)\n    EDFM 25 m ", eF, tF);
        for (size_t i = 0; i < sC.size(); ++i) { std::printf(" %.3f", sC[i]); eC = std::max(eC, std::fabs(sC[i] - aC[i])); }
        std::printf("   (max |dT_D| %.3f, uniform 25 m matrix cells: early cooling over-predicted)\n", eC);
        // (20 segments along the fracture: max |dT_D| 0.031; the remaining error is first-order upwinding in the plane and falls with NX)
        check(okF && okC && sF.size() == 5 && eF < 0.045 && eC > 2 * eF,
              "EDFM heat extraction matches Gringarten's single-fracture solution (|dT_D| < 0.045 with near-fracture refinement; coarse matrix much worse)");
    }
}
static void testEdfm() {
    using namespace ressim;
    std::printf("== C12 EDFM: embedded fractures, Barton-Bandis cubic law, hydro-shear\n");
    const double kmat = 1e-20, w0 = 2e-4, H = 10, Lc = 195.0, p1 = 12e6, p2 = 10e6;
    auto base = [&](bool mech, double angleDeg) {
        thermal::Case c; c.nx = 21; c.ny = 5; c.nz = 1; c.dx = {10}; c.dy = {10}; c.dz = {H}; c.top = 1000;
        c.poro = {0.05}; c.permx = {kmat}; c.permz = {kmat}; c.gravity = false; c.cr = 1e-10;
        c.krw0 = 1; c.swc = 0; c.sorw = 0;
        c.p0 = {11e6}; c.T0 = {350}; c.sw0 = {1.0}; c.thconr = {2.5}; c.thconsf = 0;
        thermal::Fracture f; f.name = "F1"; f.x0 = 5; f.y0 = 25; f.x1 = 205; f.y1 = 25; f.aperture = w0; f.mech = mech;
        f.kn0 = 2e10; f.vm = 3e-4; f.mu = 0.6; f.ks = 0; (void)angleDeg;
        thermal::FractureStress st; st.sHmaxGrad = 0.0; st.sHminGrad = 0.0;
        if (mech) { st.sHmaxGrad = 18e3; st.sHminGrad = 18e3; }   // isotropic 18 MPa at 1005 m: s'_n = 18 - p
        std::string err;
        if (!thermal::addFractures(c, {f}, st, &err)) std::printf("  addFractures: %s\n", err.c_str());
        return c;
    };
    // analytic steady flow through the fracture (series of half-cell resistances between the end pieces' centres)
    const double mu = steam::viscosity(350, 1.0 / steam::liquid(350, 11e6).v);
    double eRigid = 1, eBB = 1, qR = 0, qB = 0, qRa = 0, qBa = 0;
    for (int mode = 0; mode < 2; ++mode) {
        thermal::Case c = base(mode == 1, 0);
        const int nf = int(c.fracCells.size());
        if (nf != 21) { std::printf("  unexpected fracture cell count %d\n", nf); break; }
        const int fa = c.fracCells.front().cell, fb = c.fracCells.back().cell;
        thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {fa}; inj.wi = {1e-6}; inj.zinj = {1.0}; inj.bhpControl = true;
        inj.bhpMax = p1; inj.rate = 1e3; inj.quality = -1; inj.tinj = 350;
        thermal::Well pr; pr.name = "PROD"; pr.cells = {fb}; pr.wi = {1e-6}; pr.bhp = p2;
        c.wells = {inj, pr};
        thermal::Simulator s; std::string err;
        bool ok = s.setup(c, &err) && s.advanceTo(5 * 86400.0, 86400.0);
        if (!ok) { std::printf("  mode %d: %s\n", mode, err.c_str()); continue; }
        const double rho = 1.0 / steam::liquid(350, 11e6).v;
        const double q = s.wellReports()[1].molar[0] * 0.018015 / rho;   // produced m3/s
        const double pa = s.pressure()[size_t(fa)], pb = s.pressure()[size_t(fb)];
        double qa;
        if (mode == 0) qa = w0 * w0 / 12 * w0 * H * (pa - pb) / (mu * Lc);
        else {   // q = H / (12 mu Lc) int_pb^pa w(p)^3 dp,  w = wOpen - Vm s' / (Kn0 Vm + s'),  s' = 18.09 MPa - p
            const auto& F = c.fracCells.front();
            auto w = [&](double p) { const double se = std::max(0.0, F.sn0 - p); return F.wOpen - F.vm * se / (F.kn0 * F.vm + se); };
            double I = 0; const int M = 2000;
            for (int k = 0; k < M; ++k) { const double pm = pb + (k + 0.5) * (pa - pb) / M; I += std::pow(w(pm), 3) * (pa - pb) / M; }
            qa = H / (12 * mu * Lc) * I;
            const std::vector<double> wa = s.fractureApertures();
            std::printf("  Barton-Bandis apertures: inlet %.1f um, outlet %.1f um (initial %.0f um at 11 MPa)\n", wa.front() * 1e6, wa.back() * 1e6, w0 * 1e6);
        }
        const double e = std::fabs(q - qa) / qa;
        (mode == 0 ? eRigid : eBB) = e; (mode == 0 ? qR : qB) = q; (mode == 0 ? qRa : qBa) = qa;
    }
    std::printf("  steady fracture flow: rigid %.4e m3/s (cubic law %.4e, error %.2f %%), Barton-Bandis %.4e (integral %.4e, error %.2f %%)\n",
                qR, qRa, 100 * eRigid, qB, qBa, 100 * eBB);
    check(eRigid < 0.01 && eBB < 0.02, "EDFM: steady fracture flow matches the cubic law (rigid) and the Barton-Bandis stress-dependent integral");
    {   // fracture-matrix transfer: pressure diffusion into the matrix from a fracture held at constant pressure (1D linear
        // flow away from the plane, early-time sqrt(t) law): cumulative influx Q = 2 A dp phi c_t sqrt(eta t / pi) (two faces)
        const int refine = std::getenv("EDFM_REF") ? std::atoi(std::getenv("EDFM_REF")) : 1;
        thermal::Case c; c.nx = 1; c.ny = 81 * refine; c.nz = 1; c.dx = {10}; c.dy = {0.5 / refine}; c.dz = {10}; c.top = 1000;
        const double km = 1e-17, phi = 0.1, ct = 1e-9;
        c.poro = {phi}; c.permx = {km}; c.permy = {km}; c.permz = {km}; c.gravity = false; c.cr = ct - 4.4e-10;   // + water ~4.4e-10
        c.krw0 = 1; c.swc = 0; c.sorw = 0; c.p0 = {10e6}; c.T0 = {350}; c.sw0 = {1.0};
        thermal::Fracture f; f.x0 = -1; f.y0 = 20.25 + (refine % 2 == 0 ? 0.25 / refine : 0.0); f.x1 = 11; f.y1 = f.y0; f.aperture = 1e-3;
        const bool nofrac = std::getenv("EDFM_NOFRAC") != nullptr;
        std::string err; bool ok = nofrac ? true : thermal::addFractures(c, {f}, thermal::FractureStress(), &err);
        const int fc = nofrac ? 0 : ok ? c.fracCells.front().cell : 0;
        thermal::Well w; w.name = "W"; w.injector = true; w.cells = {fc}; w.wi = {1e-6}; w.zinj = {1.0}; w.bhpControl = true; w.bhpMax = 11e6; w.rate = 1e3;
        w.quality = -1; w.tinj = 350;
        c.wells = {w};
        thermal::Simulator s;
        ok = ok && s.setup(c, &err);
        const double cw = [&] { const steam::Props a = steam::liquid(350, 10.4e6), b = steam::liquid(350, 10.6e6); return -(b.v - a.v) / (0.5 * (a.v + b.v)) / 2e5; }();
        const double ctot = c.cr + cw, mu2 = steam::viscosity(350, 1.0 / steam::liquid(350, 10.5e6).v);
        const double eta = km / (mu2 * phi * ctot), tEnd = 2 * 86400.0;
        double M0 = ok ? s.fluidMoles(1) : 0;   // water (nc = 1: index 0 is ... the single water component) - use total moles below
        (void)M0;
        double N0 = 0; if (ok) for (int q = 0; q < s.size(); ++q) N0 += s.moles()[size_t(q) * 3 + 1];
        ok = ok && s.advanceTo(tEnd, std::getenv("EDFM_DT") ? std::atof(std::getenv("EDFM_DT")) : 600.0);
        double N1 = 0; if (ok) for (int q = 0; q < s.size(); ++q) N1 += s.moles()[size_t(q) * 3 + 1];
        const double rho = 1.0 / steam::liquid(350, 10.5e6).v / 0.018015;
        const double pFar = ok ? s.pressure()[size_t(c.ny - 1)] : 10e6;
        if (std::getenv("EDFM_PROF") && ok) { for (int q = 0; q < c.ny; q += std::max(1, c.ny / 20)) std::printf("    y %6.2f  dp %.4f\n", (q + 0.5) * 0.5 / refine, (s.pressure()[size_t(q)] - pFar) * 1e-6); }
        const double dpF = ok ? s.pressure()[size_t(fc)] - 10e6 : 0;   // step from the initial pressure
        std::printf("  (fracture %.4f MPa, far field %.4f MPa, t %.1f s, steps %d)\n", ok ? s.pressure()[size_t(fc)] * 1e-6 : 0.0, pFar * 1e-6, s.time(), s.stats().steps);
        const double Qa = (nofrac ? 0.5 : 1.0) * 2 * 10 * 10 * 2 * dpF * phi * ctot * std::sqrt(eta * tEnd / M_PI) * rho;   // two faces x 2 dp phi c_t sqrt(eta t / pi)
        const double e = ok ? std::fabs((N1 - N0) - Qa) / Qa : 1;
        std::printf("  fracture-matrix transfer: influx %.4g mol vs sqrt(t) solution %.4g mol (error %.1f %%; eta t / L^2 = %.2f) %s\n",
                    N1 - N0, Qa, 100 * e, eta * tEnd / 400.0, ok ? "" : err.c_str());
        check(ok && e < 0.05, "EDFM fracture-matrix transmissibility 2 A k / <d>: early-time linear-flow influx within 5 %");
    }
    {   // hydro-shear: cold water into a critically stressed fracture (SHmax at 30 deg, Coulomb) - cooling unloads the plane,
        // slip and irreversible dilation raise the conductivity, which stays after the injection stops
        thermal::Case c; c.nx = 21; c.ny = 5; c.nz = 1; c.dx = {10}; c.dy = {10}; c.dz = {H}; c.top = 3000;
        c.poro = {0.01}; c.permx = {1e-19}; c.permz = {1e-19}; c.gravity = false; c.cr = 1e-10;
        c.krw0 = 1; c.swc = 0; c.sorw = 0; c.p0 = {30e6}; c.T0 = {470}; c.sw0 = {1.0}; c.thconr = {3.0}; c.thconsf = 0;
        thermal::Fracture f; f.x0 = 5; f.y0 = 25; f.x1 = 205; f.y1 = 25; f.aperture = 1e-4; f.mech = true;
        f.kn0 = 5e10; f.vm = 2e-4; f.mu = 0.6; f.ks = 5e10; f.dilationDeg = 5;
        thermal::FractureStress st; st.sHmaxGrad = 26e3; st.sHminGrad = 14e3; st.azimuthDeg = 60; st.E = 50e9; st.nu = 0.25; st.alphaT = 8e-6;
        std::string err;
        bool ok = thermal::addFractures(c, {f}, st, &err);
        const auto F0 = ok ? c.fracCells.front() : thermal::Case::FracCell();
        const int fa = ok ? c.fracCells.front().cell : 0, fb = ok ? c.fracCells.back().cell : 0;
        thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {fa}; inj.wi = {1e-7}; inj.zinj = {1.0};
        inj.rate = 2.0; inj.bhpMax = 36e6; inj.quality = -1; inj.tinj = 320;
        thermal::Well pr; pr.name = "PROD"; pr.cells = {fb}; pr.wi = {1e-7}; pr.bhp = 28e6;
        c.wells = {inj, pr};
        thermal::Simulator s;
        ok = ok && s.setup(c, &err);
        const double margin0 = F0.tau0 - F0.mu * (F0.sn0 - F0.p0);
        ok = ok && s.advanceTo(60 * 86400.0, 2 * 86400.0);
        std::vector<double> wHot = ok ? s.fractureApertures() : std::vector<double>();
        double tmin = 1e9; if (ok) for (int k = 0; k < 21; ++k) tmin = std::min(tmin, s.temperature()[size_t(c.fracCells[size_t(k)].cell)]);
        double wsh = 0; if (ok) for (const auto& fc : s.model().fracCells) wsh = std::max(wsh, fc.wShear);
        if (ok) { s.wells()[0].open = false; ok = s.advanceTo(200 * 86400.0, 5 * 86400.0); }
        const std::vector<double> wEnd = ok ? s.fractureApertures() : std::vector<double>();
        std::printf("  hydro-shear: initial Coulomb margin %.2f MPa (stable), coldest fracture cell %.1f K, max shear dilation %.1f um, inlet aperture %.0f -> %.0f um (after shut-in %.0f um) %s\n",
                    margin0 * 1e-6, tmin, wsh * 1e6, F0.w0 * 1e6, wHot.empty() ? 0.0 : wHot.front() * 1e6, wEnd.empty() ? 0.0 : wEnd.front() * 1e6, ok ? "" : err.c_str());
        check(ok && margin0 < 0 && wsh > 1e-6 && !wEnd.empty() && wEnd.front() > F0.w0 * 1.05,
              "EDFM hydro-shear: thermal unloading triggers Coulomb slip, irreversible dilation persists after shut-in");
    }
    {   // FRACTURE / FRACSTRESS deck keywords: round trip and the EGS doublet example (well - fracture connections)
        const Case c = makeExample("egs_edfm");
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        const auto& F = b.thermal.fractures;
        const bool rt = F.size() == 3 && F[1].name == "OBL1" && std::fabs(F[1].y1 - 235) < 1e-9 && std::fabs(F[0].aperture - 0.4e-3) < 1e-15 &&
                        F[2].mech && std::fabs(F[2].ks - 2e9) < 1 && std::fabs(b.thermal.fracStress.sHminGrad - 17e3) < 1e-9 &&
                        std::fabs(b.thermal.fracStress.azimuthDeg - 30) < 1e-12;
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out;
        const bool built = thermal::deckModel(b, m, &err);
        int nfc = built ? int(m.base.fracCells.size()) : 0, nwc = 0;
        if (built) for (const auto& w : m.base.wells) for (int q : w.cells) if (q >= m.base.nx * m.base.ny * m.base.nz) ++nwc;
        double tMinF = 1e9, wMax = 0, w0 = 0, cumI = 0, cumP = 0;
        const bool ok = built && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& sim) {
            const auto& fc = sim.model().fracCells;
            const std::vector<double> wf = sim.fractureApertures();
            for (size_t f = 0; f < fc.size(); ++f) { tMinF = std::min(tMinF, sim.temperature()[size_t(fc[f].cell)]); wMax = std::max(wMax, wf[f]); w0 = std::max(w0, fc[f].w0); }
            cumI = cumP = 0;
            for (const auto& r : sim.wellReports()) {
                for (double v : r.cumMolI) cumI += v;
                for (double v : r.cumMolP) cumP += v;
            }
            return true;
        }, &err, 4);
        std::printf("  EGS + EDFM deck: round trip %d, %d fracture cells, %d well-fracture connections, coldest fracture cell %.1f C, max aperture %.0f um (initial %.0f), "
                    "injected %.3g kmol, produced %.3g kmol %s\n", int(rt), nfc, nwc, tMinF - 273.15, wMax * 1e6, w0 * 1e6, cumI * 1e-3, cumP * 1e-3, ok ? "" : err.c_str());
        check(ok && rt && nfc > 100 && nwc == 6 && tMinF < 473.15 - 50 && wMax > w0 && cumP > 0.5 * cumI,
              "FRACTURE / FRACSTRESS deck round trip; EGS doublet: wells connect to the fractures, cold-water circulation cools and opens them");
    }
    {   // C13: parametric EGS designer - horizontal doublet + transverse stages; well-stream temperatures and thermal power
        EgsDesign d; d.stages = 3; d.years = 1; d.cell = 30;
        const Case c = makeEgsCase(d);
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out;
        const bool built = thermal::deckModel(c, m, &err);
        int nwc = 0; bool geom = built;
        if (built) {
            for (const auto& w : m.base.wells) for (int q : w.cells) if (q >= m.base.nx * m.base.ny * m.base.nz) ++nwc;
            double Lx = 0; for (int i = 0; i < c.grid.nx; ++i) Lx += c.grid.dx[size_t(i)];
            for (const auto& fc : m.base.fracCells)
                geom = geom && fc.xa >= 0 && fc.xb <= Lx + 1e-6 && fc.zb > fc.zt && std::fabs(std::hypot(fc.xb - fc.xa, fc.yb - fc.ya) - fc.segLen) < 1e-6 &&
                       std::fabs(fc.zb - fc.zt - fc.segH) < 1e-6;
        }
        std::vector<double> tP, pw;
        double tI = 0;
        const bool ok = built && thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator&) {
            if (S.time <= 0) return true;
            tP.push_back(S.wells[0].tfluid); tI = S.wells[1].tfluid;
            pw.push_back((S.wells[0].hrate + S.wells[1].hrate) * 1e3 / 86400.0 * 1e-6);   // kJ/d -> MW
            return true;
        }, &err);
        bool mono = tP.size() >= 3;
        for (size_t i = 1; i < tP.size(); ++i) mono = mono && tP[i] <= tP[i - 1] + 1e-6;
        std::printf("  EGS designer (3 stages, 1 year): %zu fracture cells, %d well-fracture connections, production T %.1f -> %.1f C, injection %.1f C, \n"
                    "    net thermal power %.2f -> %.2f MW %s\n", built ? m.base.fracCells.size() : size_t(0), nwc, tP.empty() ? 0.0 : tP.front(),
                    tP.empty() ? 0.0 : tP.back(), tI, pw.empty() ? 0.0 : pw.front(), pw.empty() ? 0.0 : pw.back(), ok ? "" : err.c_str());
        check(ok && geom && nwc == 6 && mono && tP.front() < d.tRes + 1e-6 && tP.back() > d.tInj && std::fabs(tI - d.tInj) < 0.01 && pw.back() > 0 && pw.back() < pw.front(),
              "EGS designer: transverse fractures connect both horizontal wells; produced temperature falls monotonically between T_inj and T_res; "
              "positive, declining net thermal power");
    }
}
static void testNearCritical() {
    using namespace ressim;
    std::printf("== C11: IF97 region 3, continuous supercritical water, injector BHP limit, analytic verifications\n");
    {   // 1. IF97 region 3 (IAPWS-IF97 table 33 verification values) through the (T, p) property functions
        struct R { double T, p, rho, h, cp; } ref[] = {{650, 25.5837018e6, 500, 1863.43019e3, 0}, {650, 22.2930643e6, 200, 2375.12401e3, 44.6579342e3},
                                                         {750, 78.3095639e6, 500, 2258.68845e3, 0}};
        double worst = 0;
        for (const R& r : ref) {
            const steam::Props L = steam::liquid(r.T, r.p), V = steam::vapour(r.T, r.p);
            worst = std::max({worst, std::fabs(1 / L.v - r.rho) / r.rho, std::fabs(1 / V.v - r.rho) / r.rho, std::fabs(L.h - r.h) / r.h});
            if (r.cp > 0) worst = std::max(worst, std::fabs(L.cp - r.cp) / r.cp);
        }
        std::printf("  region 3 vs IF97 table 33: worst relative deviation %.2e\n", worst);
        check(worst < 1e-7, "IF97 region 3 (rho, h, cp) reproduces the IAPWS verification values");
        // continuity / monotonicity of h(T) on isobars through the critical region (the old region-1 / region-2
        // extrapolations gave negative volumes and enthalpies there)
        int bad = 0; double maxCp = 0, worstSame = 0;
        for (double p : {20e6, 22.1e6, 25e6, 35e6, 60e6}) {
            for (int br = 0; br < 2; ++br) {
                const double Tlo = br ? 560 : (p < steam::kPc ? steam::tsat(p) + 0.01 : steam::kTc + 1e-3);
                const double Thi = br ? (p < steam::kPc ? steam::tsat(p) - 0.01 : 900) : 1100;
                double hp = br ? steam::liquid(Tlo, p).h : steam::vapour(Tlo, p).h;
                for (double T = Tlo + 0.05; T < Thi; T += 0.05) {
                    const steam::Props P = br ? steam::liquid(T, p) : steam::vapour(T, p);
                    if (!(P.v > 0) || !(P.h > hp) || !std::isfinite(P.h)) ++bad;
                    maxCp = std::max(maxCp, (P.h - hp) / 0.05);
                    hp = P.h;
                    if (T > steam::kTc) worstSame = std::max(worstSame, std::fabs(steam::liquid(T, p).h - steam::vapour(T, p).h) / P.h);
                }
            }
        }
        std::printf("  isobars 20-60 MPa, 560-1100 K: %d non-monotone / invalid points, peak dh/dT %.3g J/kg/K, |h_L - h_V| above Tc %.1e\n", bad, maxCp, worstSame);
        check(bad == 0 && worstSame < 1e-12, "water properties continuous and monotone through the critical region; one fluid above Tc");
    }
    {   // 2. the flash across Tc: K_w = p_v(T) / p is C1 (the old K_w = psat / p -> 1e3 jump made near-critical cells
        //    oscillate between two flash branches); both flash paths (K-value / aqueous solutes)
        double jump[2] = {1, 1}, step[2] = {1, 1};
        for (int m = 0; m < 2; ++m) {
            thermal::Case c; thermal::Component g; g.name = "G"; g.type = 1; g.mw = 0.028;
            if (m == 1) { g.aqModel = 2; g.kh298 = 8e9; g.vInf = 3.5e-5; }
            c.comps = {g};
            const std::vector<double> z = {0.1, 0.9};
            const double p = 25e6;
            jump[m] = std::fabs(thermal::flashPT(c, p, steam::kTc + 1e-7, z).nW - thermal::flashPT(c, p, steam::kTc - 1e-7, z).nW);
            double prev = thermal::flashPT(c, p, 640, z).nW, mx = 0;
            for (double T = 640.001; T < 680; T += 0.001) { const double w = thermal::flashPT(c, p, T, z).nW; mx = std::max(mx, std::fabs(w - prev)); prev = w; }
            step[m] = mx;
        }
        std::printf("  flash at 25 MPa across Tc: |dnW| at Tc %.1e / %.1e, largest change per mK (640-680 K) %.1e / %.1e\n", jump[0], jump[1], step[0], step[1]);
        check(jump[0] < 1e-6 && jump[1] < 1e-6 && step[0] < 1e-2 && step[1] < 1e-2, "water / gas split continuous across the critical temperature (K-value and aqueous-solute flash)");
    }
    {   // 3. injector at its BHP limit in a closed box. The connection must never produce the injection stream: with the
        //    box above the limit nothing flows (the old hysteresis band let the connection back-produce a component the
        //    cell did not hold -> Newton stagnated); from below, the box fills to the limit and injection tapers to 0
        auto box = [](double p0, double rate = 2.0) {
            thermal::Case c; c.nx = 3; c.dx = {10}; c.dy = {10}; c.dz = {5}; c.poro = {0.3}; c.permx = {1e-12}; c.gravity = false; c.cr = 5e-10;
            thermal::Component g; g.name = "G"; g.type = 1; g.mw = 0.044;
            c.comps = {thermHeavy(1.0), g}; c.p0 = {p0}; c.T0 = {350}; c.sw0 = {0.3}; c.x0 = {{1.0, 0.0}};
            thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {1e-11}; inj.zinj = {0, 1, 0};
            inj.rate = rate; inj.bhpMax = 5e6; inj.tinj = 350; inj.quality = -1;
            c.wells = {inj};
            return c;
        };
        thermal::Simulator a, b; std::string err;
        const bool okA = a.setup(box(5e6 + 50), &err) && a.advanceTo(10 * 86400.0, 86400.0);
        const double injA = okA ? a.wellReports()[0].cumMolI[1] : -1;
        thermal::Simulator r;
        const bool okR = r.setup(box(3e6, 1e-3), &err) && r.advanceTo(0.01 * 86400.0, 0.01 * 86400.0);
        const double jacRate = okR ? r.jacobianCheck(600.0, nullptr, nullptr, 0) : 1;
        const bool rateMode = okR && !r.wellReports()[0].onBhp;
        bool okB = b.setup(box(4.9e6), &err);
        // BHP-mode Jacobian at day 10 (p ~10 Pa below the limit): later the converged state sits within the +-5 Pa
        // difference step of the ramp's g'' jump at x = 0, which central differences cannot resolve
        okB = okB && b.advanceTo(10 * 86400.0, 86400.0);
        double jacBhp = 1;
        if (okB) { thermal::Simulator bj = b; jacBhp = bj.jacobianCheck(86400.0, nullptr, nullptr, 0); }
        okB = okB && b.advanceTo(30 * 86400.0, 86400.0);
        const double pB = okB ? b.pressure()[0] : 0, injB = okB ? b.wellReports()[0].cumMolI[1] : 0;
        std::printf("  above the limit: ok %d, injected %.3e mol, cuts %d | from below: ok %d, p_inj %.2f Pa over the limit, injected %.1f mol, cuts %d, Jacobian rate %.1e (rate mode %d) / BHP %.1e\n",
                    int(okA), injA, okA ? a.stats().cuts : -1, int(okB), pB - 5e6, injB, okB ? b.stats().cuts : -1, jacRate, int(rateMode), jacBhp);
        check(okA && std::fabs(injA) < 1e-9 && a.stats().cuts == 0, "injector above its BHP limit: no back-flow of the injection stream, no time-step cuts");
        check(okB && pB - 5e6 < 1.0 && pB - 5e6 > -2e3 && injB > 0 && b.stats().cuts <= 2, "closed box fills to the injector BHP limit and injection stops smoothly");
        check(rateMode && jacRate < 1e-2 && jacBhp < 1e-2, "injector drawdown ramp: implicit rate-control BHP and BHP-mode Jacobians verified by differences");
    }
    {   // 4. verification: transient 1D conduction (step of 20 K in a rod) vs the erf solution
        thermal::Case c; c.nx = 100; c.dx = {0.1}; c.dy = {1}; c.dz = {1}; c.poro = {0.01}; c.permx = {1e-18}; c.gravity = false; c.cr = 0;
        c.p0 = {5e6}; c.sw0 = {1.0}; c.thconr = {2.0}; c.thconsf = 0;
        c.T0.assign(100, 300.0); for (int i = 0; i < 50; ++i) c.T0[size_t(i)] = 320.0;
        thermal::Simulator s; std::string err;
        const double t = 5 * 86400.0;
        bool ok = s.setup(c, &err) && s.advanceTo(t, 0.02 * 86400.0);
        const steam::Props W = steam::liquid(310, 5e6);
        const double rc = 0.01 * W.cp / W.v + 0.99 * c.rockCv, kap = 2.0 / rc;
        double worst = 0;
        for (int i = 0; ok && i < 100; ++i) {
            const double x = (i + 0.5) * 0.1 - 5.0;
            const double Ta = 310.0 - 10.0 * std::erf(x / (2 * std::sqrt(kap * t)));
            worst = std::max(worst, std::fabs(s.temperature()[size_t(i)] - Ta));
        }
        std::printf("  1D conduction after 5 d (diffusion length %.2f m): max |T - T_erf| %.4f K of a 20 K step\n", 2 * std::sqrt(kap * t), worst);
        check(ok && worst < 0.1, "transient heat conduction matches the analytic erf solution (< 0.5% of the step)");
    }
    {   // 5. verification: Theis line-source drawdown around a constant-rate producer (black-oil core, 61 x 61 grid)
        const char* deck =
            "RUNSPEC\nTITLE\nTHEIS\nDIMENS\n 61 61 1 /\nOIL\nWATER\nMETRIC\nTABDIMS\n/\nWELLDIMS\n 1 1 1 1 /\nSTART\n 1 JAN 2020 /\n"
            "GRID\nDX\n 3721*20 /\nDY\n 3721*20 /\nDZ\n 3721*10 /\nTOPS\n 3721*1000 /\nPERMX\n 3721*100 /\nPERMY\n 3721*100 /\nPERMZ\n 3721*100 /\nPORO\n 3721*0.2 /\n"
            "PROPS\nPVTW\n 200 1.0 1e-5 0.5 0 /\nPVCDO\n 200 1.0 1e-4 1.0 0 /\nROCK\n 200 1e-5 /\nDENSITY\n 800 1000 1 /\nSWOF\n 0 0 1 0\n 1 1 0 0 /\n"
            "SOLUTION\nEQUIL\n 1005 200 2000 0 0 0 /\nSCHEDULE\nWELSPECS\n P G 31 31 1005 OIL /\n/\nCOMPDAT\n P 31 31 1 1 OPEN 1* 1* 0.2 /\n/\n"
            "WCONPROD\n P OPEN ORAT 100 4* 1 /\n/\nTSTEP\n 40*0.025 /\nEND\n";
        bool ok = true; double worst = 1, dpRef = 0;
        try {
            const Case c = eclipseDeckFromString(deck);
            Simulator sim; sim.setup(c); sim.initialize();
            while (sim.runReportStep()) {}
            const auto& sn = sim.snapshots();
            const Grid& G = sim.grid();
            const double k = 100 * 9.869233e-16, h = 10, mu = 1e-3, phi = 0.2, ct = 1.1e-9, q = 100.0 / 86400;
            worst = 0;
            for (const Snapshot& S : sn) {
                if (!(std::fabs(S.time - 0.5) < 1e-6 || std::fabs(S.time - 1.0) < 1e-6)) continue;
                for (int off : {3, 5, 8}) {
                    const int a = G.activeOf[size_t(30 * 61 + 30 + off)], a0 = G.activeOf[size_t(30 * 61 + 30 + off)];
                    (void)a0;
                    const double r = off * 20.0, tt = S.time * 86400;
                    const double dpA = q * mu / (4 * M_PI * k * h) * expint1(phi * mu * ct * r * r / (4 * k * tt)) / 1e5;
                    const double dpN = sn.front().pressure[size_t(a)] - S.pressure[size_t(a)];
                    worst = std::max(worst, std::fabs(dpN - dpA) / dpA);
                    if (std::getenv("THEIS_DEBUG")) std::printf("    t %.2f r %.0f num %.5f theis %.5f\n", S.time, r, dpN, dpA);
                    if (off == 5 && std::fabs(S.time - 1.0) < 1e-6) dpRef = dpA;
                }
            }
        } catch (const std::exception& e) { ok = false; std::printf("  Theis deck: %s\n", e.what()); }
        std::printf("  Theis drawdown at r = 60 / 100 / 160 m, t = 0.5 / 1 d: worst deviation %.2f%% (dp(100 m, 1 d) = %.3f bar)\n", 100 * worst, dpRef);
        check(ok && worst < 0.03, "radial drawdown matches the Theis line-source solution within 3%");
    }
}

// ------------------------------------------------------------------------------------------------ THM / THMC (C8)
static void testThmc() {
    std::printf("[THM / THMC: geomechanics, reactions, combustion, UCG]\n");
    using namespace ressim;
    {   // 1. uniaxial-strain column (rollers on sides / base, free top): exact FEM answers for uniform dp and dT
        //    u_top = H (alpha dp + 3 K beta dT) / M,  M = K + 4G/3; Lagrangian porosity phi0 + alpha eps_v
        double err = 0, phiErr = 0, cpErr = 0;
        for (int mode = 0; mode < 2; ++mode) {
            geomech::Input in; in.enabled = true; in.E = {5e9}; in.nu = {0.25}; in.biot = {1.0}; in.alphaT = {1e-5};
            const int nx = 3, ny = 3, nz = 10;
            const int n = nx * ny * nz;
            geomech::Solver g; std::string e;
            if (!g.setup(in, nx, ny, nz, {50}, {50}, {10}, 1000, std::vector<double>(size_t(n), 1e7), std::vector<double>(size_t(n), 0.2), &e)) { err = 1; break; }
            std::vector<double> dp(size_t(n), mode == 0 ? 1e6 : 0.0), dT(size_t(n), mode == 1 ? 10.0 : 0.0);
            g.solve(dp, dT);
            const double M = 6e9, load = 1e6, H = 100;
            err = std::max(err, std::fabs(g.surfaceUplift() - load * H / M) / (load * H / M));
            const double ev = load / M;
            for (int q = 0; q < n; ++q) {
                const double want = mode == 0 ? 0.2 + ev : 0.2 + ev - 3e-5 * 10 * (1.0 - 0.2);   // beta_phi = 3 beta (alpha - phi0)
                phiErr = std::max(phiErr, std::fabs(g.porosity()[size_t(q)] - want));
            }
            cpErr = std::fabs(g.cp(0) - 1.0 / geomech::bulkModulus(5e9, 0.25)) * geomech::bulkModulus(5e9, 0.25);
        }
        std::printf("  uniaxial column: uplift error %.1e, porosity error %.1e, c_p error %.1e\n", err, phiErr, cpErr);
        check(err < 1e-6 && phiErr < 1e-9 && cpErr < 1e-12, "geomechanics: uniaxial-strain uplift (pressure and thermal loads), Lagrangian porosity, fixed-stress c_p exact");
    }
    {   // 1b. boundary load (BCMECH-like): 1 MPa surface load on the free top of the uniaxial column - exact settlement s H / M,
        //     vertical effective stress +s, horizontal +s nu / (1 - nu); fixed faces and the deck round trip of MECHBC
        geomech::Input in; in.enabled = true; in.E = {5e9}; in.nu = {0.25}; in.biot = {1.0}; in.alphaT = {1e-5};
        const int nx = 3, ny = 3, nz = 10, n = nx * ny * nz;
        geomech::Solver g0, g1; std::string e;
        geomech::Input in1 = in; in1.bcType[4] = 1; in1.bcStress[4] = 1e6;
        const std::vector<double> p0(size_t(n), 1e7), ph(size_t(n), 0.2), z(size_t(n), 0.0);
        const bool ok = g0.setup(in, nx, ny, nz, {50}, {50}, {10}, 1000, p0, ph, &e) && g1.setup(in1, nx, ny, nz, {50}, {50}, {10}, 1000, p0, ph, &e);
        if (ok) { g0.solve(z, z); g1.solve(z, z); }
        const double M = 6e9, s = 1e6, H = 100;
        double eU = 1, eV = 1, eH = 1;
        if (ok) {
            eU = std::fabs((g1.dispZ()[0] - g0.dispZ()[0]) - s * (H - 5) / M) / (s * H / M);   // top-layer centre (5 m below the surface), + down
            eV = 0; eH = 0;
            for (int q = 0; q < n; ++q) {
                eV = std::max(eV, std::fabs(g1.sigEffV()[size_t(q)] - g0.sigEffV()[size_t(q)] - s) / s);
                eH = std::max(eH, std::fabs(g1.sigEffH()[size_t(q)] - g0.sigEffH()[size_t(q)] - s / 3) / s);
            }
        }
        Case c = makeExample("sagd_thm");
        c.thermal.geo.bcType[1] = 1; c.thermal.geo.bcStress[1] = 2.5e6; c.thermal.geo.bcType[5] = 2;
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        const bool rt = b.thermal.geo.bcType[1] == 1 && b.thermal.geo.bcStress[1] == 2.5e6 && b.thermal.geo.bcType[5] == 2 && b.thermal.geo.bcType[0] == -1;
        std::printf("  surface load 1 MPa: settlement error %.1e, d sigma'v error %.1e, d sigma'h error %.1e, MECHBC round trip %d\n", eU, eV, eH, int(rt));
        check(ok && eU < 1e-6 && eV < 1e-6 && eH < 1e-6 && rt, "geomechanics boundary loads (BCMECH-like): exact oedometric settlement and stress change, deck round trip");
    }
    {   // 2. RESSIM deck round trip of reactions / solids / streams / permeability law / geomechanics
        int ok = 0;
        for (const char* id : {"combustion", "ucg", "sagd_thm"}) {
            const Case c = makeExample(id);
            const Case b = eclipseDeckFromString(eclipseDeckToString(c));
            const ThermalInput& A = c.thermal; const ThermalInput& B = b.thermal;
            bool same = B.enabled && B.nc() == A.nc() && B.reactions.size() == A.reactions.size() && B.permLaw == A.permLaw &&
                        B.streams.size() == A.streams.size() && B.geo.enabled == A.geo.enabled && B.solidConc.size() == A.solidConc.size();
            for (int i = 0; same && i < A.nc(); ++i) same = A.comps[size_t(i)].type == B.comps[size_t(i)].type;
            for (size_t r = 0; same && r < A.reactions.size(); ++r) {
                const ThermalReaction& x = A.reactions[r]; const ThermalReaction& y = B.reactions[r];
                same = x.name == y.name && std::fabs(x.A - y.A) <= 1e-9 * std::fabs(x.A) && std::fabs(x.Ea - y.Ea) <= 1e-9 * std::fabs(x.Ea) &&
                       std::fabs(x.dH - y.dH) <= 1e-9 * std::fabs(x.dH) && x.o2pp == y.o2pp && x.reversible == y.reversible && x.sto == y.sto && x.prod == y.prod;
            }
            for (size_t k = 0; same && k < A.solidConc.size(); ++k) same = std::fabs(A.solidConc[k] - B.solidConc[k]) <= 1e-9 * std::fabs(A.solidConc[k]);
            if (same && A.geo.enabled) same = std::fabs(A.geo.E[0] - B.geo.E[0]) < 1e-6 * A.geo.E[0] && A.geo.overLayers == B.geo.overLayers && std::fabs(A.geo.k0 - B.geo.k0) < 1e-12;
            size_t nInj = 0; for (const auto& s : b.schedule) nInj += s.injStreams.size();
            if (!A.streams.empty()) same = same && nInj > 0;
            if (!same) std::printf("  %s: RESSIM thermal deck round trip mismatch\n", id);
            ok += same;
        }
        check(ok == 3, "THMC deck: export -> import keeps reactions, solid types / concentrations, WELLSTRE / WINJGAS, PERMLAW, GEOMECH");
    }
    {   // 3. combustion tube (6 h): a combustion front propagates, oxygen is consumed, energy closes with reaction heat
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        double e0 = 0, e1 = 0, hin = 0, loss = 0, rx = 0, gross = 0, front = 0, o2in = 0, o2out = 0; bool first = true;
        const bool ran = thermal::deckModel(makeExample("combustion"), m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
            if (first) { e0 = s.energyInPlace(); first = false; }
            e1 = s.energyInPlace(); hin = s.cumEnthalpyIn(); loss = s.cumHeatLoss(); rx = s.cumReactionHeat();
            gross = 0; for (const auto& r : s.wellReports()) gross += std::fabs(r.cumH);
            front = 0; for (int q = 0; q < s.size(); ++q) if (s.temperature()[size_t(q)] > 573.15) front = q;
            o2in = o2out = 0;
            for (size_t w = 0; w < s.wellReports().size(); ++w) (s.model().wells[w].injector ? o2in : o2out) += std::fabs(s.wellReports()[w].molar[3]);
            return true;
        }, &err, 6, &st);
        const double ebal = std::fabs((e1 - e0) - (hin - loss + rx)) / std::max(1.0, gross + std::fabs(rx));
        std::printf("  combustion tube 6 h: hottest-front cell %.0f, O2 in %.3g / out %.3g mol/s, reaction heat %.3g MJ, energy balance %.1e, %d steps %s\n",
                    front, o2in, o2out, rx * 1e-6, ebal, st.steps, ran ? "" : err.c_str());
        check(ran && front >= 6 && o2out < 0.8 * o2in && rx > 0 && ebal < 1e-5, "in-situ combustion: front propagates (> 300 C beyond the igniter), O2 consumed, energy balance with reaction heat < 1e-5");
    }
    {   // 4. UCG (8 days): coal pyrolysis / char gasification produce syngas, the cavity opens (porosity, permeability)
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        double co = 0, h2 = 0, phiMax = 0, kMax = 0, tMax = 0;
        const bool ran = thermal::deckModel(makeExample("ucg"), m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
            co = h2 = 0;
            for (size_t w = 0; w < s.wellReports().size(); ++w) if (!s.model().wells[w].injector) { co += s.wellReports()[w].molar[5]; h2 += s.wellReports()[w].molar[6]; }
            for (double v : s.fluidPorosity()) phiMax = std::max(phiMax, v);
            for (double v : s.permMultiplier()) kMax = std::max(kMax, v);
            for (double T : s.temperature()) tMax = std::max(tMax, T);
            return true;
        }, &err, 4, &st);
        std::printf("  UCG 8 d: CO %.3g H2 %.3g mol/s, cavity porosity %.3f, k multiplier %.3g, Tmax %.0f C, %d steps %s\n", co, h2, phiMax, kMax, tMax - 273.15, st.steps, ran ? "" : err.c_str());
        check(ran && co > 0.01 && h2 > 0.01 && phiMax > 0.3 && kMax > 10 && tMax > 873, "UCG: ignition, syngas (CO, H2) production, coal consumption opens the cavity (porosity, Kozeny-Carman permeability)");
    }
    {   // 5. SAGD THM (1 report step): steam chamber heave and thermal stress change from the coupled mechanics
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        double up = 0, upRes = 0, e0 = 0, e1 = 0, hin = 0, loss = 0, gross = 0; bool first = true;
        const bool ran = thermal::deckModel(makeExample("sagd_thm"), m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
            if (first) { e0 = s.energyInPlace(); first = false; }
            e1 = s.energyInPlace(); hin = s.cumEnthalpyIn(); loss = s.cumHeatLoss();
            gross = 0; for (const auto& r : s.wellReports()) gross += std::fabs(r.cumH);
            if (const geomech::Solver* g = s.mechanics()) { up = g->surfaceUplift(); upRes = g->reservoirTopUplift(); }
            return true;
        }, &err, 1, &st);
        const double ebal = std::fabs((e1 - e0) - (hin - loss)) / std::max(1.0, gross + std::fabs(loss));
        std::printf("  SAGD THM 91 d: surface uplift %.2f mm, reservoir top %.2f mm, energy balance %.1e, %d steps, %d mechanics solves %s\n",
                    up * 1e3, upRes * 1e3, ebal, st.steps, st.mechSolves, ran ? "" : err.c_str());
        check(ran && up > 1e-4 && upRes > up && st.mechSolves >= st.steps && ebal < 1e-5, "SAGD THM: fixed-stress coupled heave (surface < reservoir top), energy balance < 1e-5");
    }
}

// ------------------------------------------------------------------------------------------------ C10: unified solubility / trapping
static ressim::thermal::Component thermCO2(bool soluble) {
    ressim::thermal::Component q; q.name = "CO2"; q.mw = 0.04401; q.type = 1; q.avg = 1.5e-5 / std::pow(300.0, 0.8); q.bvg = 0.8; q.cpg1 = 37.1;
    const ressim::gasdata::GasInfo* gi = ressim::gasdata::find("CO2");
    if (gi) {
        q.tcG = gi->tc; q.pcG = gi->pc; q.omegaG = gi->omega; q.formula = gi->formula;
        if (soluble) { q.aqModel = 1; q.hA = gi->hA; q.hB = gi->hB; q.hC = gi->hC; q.vInf = gi->vInf; q.ksalt = gi->ksalt; }
    }
    return q;
}
static void testSolubility() {
    using namespace ressim;
    std::printf("Unified solubility / PR gas / trapping (C10: IAPWS G7-04 Henry, three-phase convex flash, inventory)\n");
    {   // IAPWS G7-04 check values ln(kH / 1 GPa)
        struct R { const char* n; double T, v; };
        const R ref[] = {{"CO2", 300, -1.7508}, {"CO2", 400, -0.5450}, {"CO2", 500, -0.6524}, {"CO2", 600, -1.3489}, {"H2", 300, 1.9702}, {"H2", 400, 1.8464},
                         {"H2", 500, 1.0513}, {"H2", 600, -0.1848}, {"CH4", 300, 1.4034}, {"CH4", 500, 1.0342}, {"N2", 400, 2.3509}, {"N2", 600, 0.1647}};
        double worst = 0;
        for (const R& r : ref) {
            const gasdata::GasInfo* g = gasdata::find(r.n);
            const double v = g ? std::log(gasdata::henryIapws(g->hA, g->hB, g->hC, r.T) / 1e9) : 1e9;
            worst = std::max(worst, std::fabs(v - r.v));
        }
        // analytic d ln kH / dT vs central differences
        const gasdata::GasInfo* g = gasdata::find("CO2");
        const double T = 350, h = 1e-3;
        const double fd = (std::log(gasdata::henryIapws(g->hA, g->hB, g->hC, T + h)) - std::log(gasdata::henryIapws(g->hA, g->hB, g->hC, T - h))) / (2 * h);
        const double an = gasdata::henryIapwsDlnDT(g->hA, g->hB, g->hC, T);
        std::printf("  G7-04: worst |ln kH - ref| %.2e, dlnkH/dT %.6e vs FD %.6e, psat(373.124 K) %.1f Pa\n", worst, an, fd, gasdata::psatWagnerPruss(373.124));
        check(worst < 2e-4 && std::fabs(an - fd) < 1e-6 * std::fabs(fd) + 1e-9 && std::fabs(gasdata::psatWagnerPruss(373.124) - 101325) < 50,
              "IAPWS G7-04 Henry constants match the published check values; analytic temperature derivative");
    }
    {   // element formulas / balances
        std::map<std::string, double> el;
        const bool p1 = gasdata::parseFormula("Ca(OH)2", el) && el["Ca"] == 1 && el["O"] == 2 && el["H"] == 2;
        const double i1 = gasdata::elementImbalance({"CH4", "O2", "CO2", "H2O"}, {1, 2, 0, 0}, {0, 0, 1, 2});
        const double i2 = gasdata::elementImbalance({"C", "H2O", "CO", "H2"}, {1, 1, 0, 0}, {0, 0, 1, 0.5});   // H unbalanced
        std::printf("  formulas: Ca(OH)2 %s, CH4 combustion imbalance %.1e, broken WGS-type imbalance %.2f, CO3 in dolomite %.0f\n",
                    p1 ? "ok" : "bad", i1, i2, gasdata::carbonateCO2("CaMg(CO3)2"));
        check(p1 && i1 < 1e-12 && std::fabs(i2 - 1.0) < 1e-12 && gasdata::carbonateCO2("CaCO3") == 1 && gasdata::carbonateCO2("CaMg(CO3)2") == 2,
              "element formula parser, reaction element balance and carbonate CO2 content");
    }
    {   // the general three-phase flash reproduces the legacy free-water flash when nothing dissolves
        thermal::Case c;
        thermal::Component h = thermHeavy(1.0);
        thermal::Component l; l.name = "LIGHT"; l.mw = 0.044; l.kv2 = 2e9; l.kv4 = -2500; l.rhoRef = 600;
        thermal::Component n2; n2.name = "N2"; n2.mw = 0.028; n2.type = 1;
        c.comps = {h, l, n2};
        struct S { double p, T; std::vector<double> z; };
        const S st[] = {{3e6, 450, {0.3, 0.2, 0.05, 0.45}}, {3e6, 300, {0.3, 0.2, 0.05, 0.45}}, {1e6, 460, {0.01, 0.01, 0.0, 0.98}},
                        {5e6, 520, {0.2, 0.3, 0.1, 0.4}}, {2e5, 480, {0.0, 0.5, 0.2, 0.3}}, {3e6, 700, {0.2, 0.1, 0.1, 0.6}}};
        double worst = 0;
        for (const S& x : st) {
            const thermal::FlashState a = thermal::flashPT(c, x.p, x.T, x.z), b = thermal::flashGeneral(c, x.p, x.T, x.z);
            worst = std::max({worst, std::fabs(a.nW - b.nW), std::fabs(a.nO - b.nO), std::fabs(a.nG - b.nG)});
            for (size_t i = 0; i < a.x.size(); ++i) if (a.nO > 0) worst = std::max(worst, std::fabs(a.x[i] - b.x[i]));
            for (size_t i = 0; i < a.y.size(); ++i) if (a.nG > 0) worst = std::max(worst, std::fabs(a.y[i] - b.y[i]));
        }
        std::printf("  general vs free-water flash (6 states, insoluble components): max difference %.2e\n", worst);
        check(worst < 1e-9, "convex three-phase flash reproduces the legacy free-water Rachford-Rice flash");
    }
    {   // CO2 solubility in water with the PR gas vs Spycher & Pruess (2005) at 50 C
        thermal::Case c; c.comps = {thermCO2(true)}; c.gasPR = true;
        const double pr[3] = {8.35e6, 15.43e6, 24.4e6}, sp[3] = {0.01835, 0.02225, 0.02389};
        double worst = 0, xs[3] = {0, 0, 0};
        for (int k = 0; k < 3; ++k) {
            const thermal::FlashState f = thermal::flashPT(c, pr[k], 323.15, {0.2, 0.8});
            xs[k] = f.a.empty() ? 0.0 : f.a[0];
            worst = std::max(worst, std::fabs(xs[k] / sp[k] - 1));
        }
        c.salinity = 1.0;
        const thermal::FlashState fb = thermal::flashPT(c, pr[0], 323.15, {0.2, 0.8});
        const double ratio = fb.a[0] / xs[0];
        std::printf("  x_CO2(aq) 50 C: %.5f / %.5f / %.5f (Spycher-Pruess %.5f / %.5f / %.5f), worst %.1f %%, 1 m NaCl ratio %.3f\n",
                    xs[0], xs[1], xs[2], sp[0], sp[1], sp[2], 100 * worst, ratio);
        check(worst < 0.08 && ratio < 0.85 && ratio > 0.7, "CO2 solubility (Henry + Krichevsky-Kasarnovsky + PR fugacity) within 8 % of Spycher-Pruess; salting-out");
    }
    {   // closed box: CO2 gas over water dissolves; mass / energy conservation, inventory closes, Jacobian
        thermal::Case c; c.nx = 1; c.dx = {10}; c.dy = {10}; c.dz = {10}; c.poro = {0.3}; c.permx = {1e-13}; c.gravity = false; c.cr = 5e-10;
        c.comps = {thermCO2(true)}; c.gasPR = true; c.track = 0;
        c.p0 = {10e6}; c.T0 = {323.15}; c.sw0 = {0.7}; c.sg0 = {0.3}; c.x0 = {{1.0}}; c.y0 = {{1.0}};
        thermal::Simulator s; std::string err;
        bool ok = s.setup(c, &err);
        const double N0 = s.fluidMoles(0), E0 = s.energyInPlace(), sg0 = ok ? s.sg()[0] : 0;
        double jac = ok ? s.jacobianCheck(86400.0) : 1;
        ok = ok && s.advanceTo(100 * 86400.0, 10 * 86400.0);
        const thermal::Inventory I = s.trackedInventory();
        const double mass = s.fluidMoles(0) * c.comps[0].mw;
        const double dN = std::fabs(s.fluidMoles(0) - N0) / N0, dE = std::fabs(s.energyInPlace() - E0) / std::max(1.0, std::fabs(E0));
        std::vector<double> x, y, a; s.phaseCompositions(x, y, a);
        const thermal::FlashState f = thermal::flashPT(s.model(), s.pressure()[0], s.temperature()[0], {s.moles()[1], s.moles()[2]});
        std::printf("  closed box: Sg %.4f -> %.4f, p %.3f MPa, x_aq %.5f (flash %.5f), dissolved %.1f%%, mass drift %.1e, |dE|/E %.1e, inventory closure %.1e, Jacobian %.1e %s\n",
                    sg0, ok ? s.sg()[0] : 0.0, s.pressure()[0] / 1e6, a.empty() ? 0.0 : a[0], f.a.empty() ? 0.0 : f.a[0], 100 * I.dissolvedAq / mass, dN, dE,
                    std::fabs(I.total() - mass) / mass, jac, ok ? "" : err.c_str());
        check(ok && dN < 1e-10 && dE < 1e-9 && I.dissolvedAq > 0 && std::fabs(I.total() - mass) < 1e-8 * mass && std::fabs(a[0] - f.a[0]) < 1e-9 &&
              s.sg()[0] < 0.3 && jac < 1e-3, "CO2 dissolution in a closed thermal cell: exact mass / energy balance, inventory closes, Jacobian consistent");
    }
    {   // 1D CO2 injection then water chase: Land / Carlson residual trapping in the thermal core (hysteresis on vs off)
        double resid[2] = {0, 0}, inPlace[2] = {0, 0}, diss[2] = {0, 0}; bool ok = true, trapOk = true;
        for (int m = 0; m < 2; ++m) {
            thermal::Case c; c.nx = 20; c.dx = {10}; c.dy = {10}; c.dz = {5}; c.poro = {0.3}; c.permx = {2e-13}; c.gravity = false; c.cr = 5e-10;
            c.comps = {thermCO2(true)}; c.gasPR = true; c.track = 0; c.hyst.on = m == 1; c.hyst.sgtMax = 0.3; c.sgc = 0.05; c.swc = 0.2;
            c.p0 = {10e6}; c.T0 = {323.15}; c.sw0 = {1.0}; c.x0 = {{1.0}};
            thermal::Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {2e-12}; inj.zinj = {1, 0}; inj.rate = 1.0; inj.bhpMax = 25e6;
            inj.quality = -1; inj.tinj = 323.15;
            thermal::Well pr; pr.name = "PROD"; pr.cells = {19}; pr.wi = {2e-12}; pr.bhp = 10e6;
            c.wells = {inj, pr};
            thermal::Simulator s; std::string err;
            bool r = s.setup(c, &err) && s.advanceTo(60 * 86400.0, 5 * 86400.0);
            if (r) { s.wells()[0].zinj = {0, 1}; s.wells()[0].rate = 10.0; s.limitNextDt(3600.0); r = s.advanceTo(120 * 86400.0, 5 * 86400.0); }
            const thermal::Inventory I = r ? s.trackedInventory() : thermal::Inventory();
            const std::vector<double> st = r ? s.trappedGasSaturation() : std::vector<double>();
            for (size_t q = 0; q < st.size(); ++q) trapOk = trapOk && st[q] <= s.sg()[q] + 1e-12;
            resid[m] = I.residual; inPlace[m] = I.total(); diss[m] = I.dissolvedAq;
            if (!r) { ok = false; std::printf("    run %d failed: %s\n", m, err.c_str()); }
        }
        std::printf("  1D CO2 + water chase (60 + 60 d): residual %.3g / %.3g t, dissolved %.3g / %.3g t, CO2 in place %.3g / %.3g t (no hysteresis / Land)\n",
                    resid[0] / 1e3, resid[1] / 1e3, diss[0] / 1e3, diss[1] / 1e3, inPlace[0] / 1e3, inPlace[1] / 1e3);
        check(ok && trapOk && resid[1] > resid[0] && resid[1] > 0 && diss[0] > 0 && inPlace[1] > inPlace[0],
              "thermal core gas hysteresis: Land-trapped residual CO2 exceeds the critical-saturation case, more CO2 retained");
    }
    {   // UCG + CO2 example: deck round trip of the C10 lines, then gasification + quench + CO2 injection (20 reports)
        const Case c = makeExample("ucg_co2");
        const Case b = eclipseDeckFromString(eclipseDeckToString(c));
        const ThermalInput& A = c.thermal; const ThermalInput& B = b.thermal;
        bool same = B.nc() == A.nc() && B.gasPR == A.gasPR && B.hystGas == A.hystGas && std::fabs(B.sgtMax - A.sgtMax) < 1e-12 && B.track == A.track &&
                    std::fabs(B.salinity - A.salinity) < 1e-12 && B.reactions.size() == A.reactions.size();
        for (int i = 0; same && i < A.nc(); ++i) {
            const ThermalComponent& x = A.comps[size_t(i)]; const ThermalComponent& y = B.comps[size_t(i)];
            same = x.aqModel == y.aqModel && std::fabs(x.hA - y.hA) < 1e-9 && std::fabs(x.hC - y.hC) < 1e-9 && std::fabs(x.vInf - y.vInf) < 1e-15 &&
                   std::fabs(x.ksalt - y.ksalt) < 1e-12 && x.formula == y.formula && std::fabs(x.trackEq - y.trackEq) < 1e-12;
        }
        thermal::DeckModel m; std::string err; std::vector<thermal::Summary> out; thermal::Stats st;
        double coMax = 0, h2Max = 0; thermal::Inventory I; double ebal = 1, e0 = 0, e1 = 0, hin = 0, loss = 0, gross = 0, rx = 0; bool first = true;
        const bool ran = thermal::deckModel(b, m, &err) && thermal::runDeckModel(m, out, [&](const thermal::Summary&, const thermal::Simulator& s) {
            for (size_t w = 0; w < s.wellReports().size(); ++w)
                if (!s.model().wells[w].injector && s.wellReports()[w].molar.size() > 6) { coMax = std::max(coMax, s.wellReports()[w].molar[5]); h2Max = std::max(h2Max, s.wellReports()[w].molar[6]); }
            if (first) { e0 = s.energyInPlace(); first = false; }
            e1 = s.energyInPlace(); hin = s.cumEnthalpyIn(); loss = s.cumHeatLoss(); rx = s.cumReactionHeat();
            gross = 0; for (const auto& r : s.wellReports()) gross += std::fabs(r.cumH);
            I = s.trackedInventory();
            return true;
        }, &err, 20, &st);
        ebal = std::fabs((e1 - e0) - (hin - loss + rx)) / std::max(1.0, gross + std::fabs(loss) + std::fabs(rx));
        std::printf("  UCG + CO2 (round trip %s): syngas peak CO %.3g H2 %.3g mol/s; CO2 after injection: dissolved %.3g, residual %.3g, mobile %.3g, mineral %.3g t; "
                    "energy balance %.1e, %d steps, gas injected %.4g sm3 %s\n", same ? "ok" : "MISMATCH", coMax, h2Max, I.dissolvedAq / 1e3, I.residual / 1e3, (I.mobile + I.structural) / 1e3,
                    I.mineral / 1e3, ebal, st.steps, out.empty() ? 0.0 : out.back().fgit, ran ? "" : err.c_str());
        check(same, "C10 deck lines (AQSOL / GASEOS / SALINITY / HYSTGAS / TRACK / FORMULA) survive export -> import");
        check(ran && coMax > 0.1 && h2Max > 0.05 && I.dissolvedAq > 0 && I.residual > 0 && I.mineral > 0 && ebal < 1e-5 && !out.empty() && out.back().fgit > 0,
              "UCG syngas production and CO2 storage in one thermal run: dissolved, residual and mineral trapping, energy balance < 1e-5");
    }
}


// ---------------------------------------------------------------- C14: restart (exact continuation) + result files
static double relDiff(double a, double b) { return std::fabs(a - b) / std::max({1e-12, std::fabs(a), std::fabs(b)}); }
static double maxRelArr(const std::vector<float>& a, const std::vector<float>& b) {
    if (a.size() != b.size()) return 1.0;
    double m = 0; for (size_t i = 0; i < a.size(); ++i) m = std::max(m, relDiff(a[i], b[i]));
    return m;
}
static void testRestart() {
    std::printf("== C14 restart: continuing from a saved state reproduces the continuous run\n");
    const std::string tmp = (std::filesystem::temp_directory_path() / "ressim_rst_test").string();
    std::filesystem::create_directories(tmp);
    auto boCase = [&](const char* label, Case c) {
      try {
        Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
        while (a.runReportStep()) {}
        const int ns = int(c.schedule.size()), k = std::max(1, ns / 2);
        Simulator b; b.log = [](const std::string&) {}; b.setup(c); b.initialize();
        for (int i = 0; i < k; ++i) b.runReportStep();
        const std::string blob = b.saveState();
        // through the .rsres file
        ResultsFile rf; rf.caseText = projectText(c); rf.summary = b.summary(); rf.snaps = b.snapshots(); rf.wellNames = b.wellNames();
        ResultsFile::Restart R; R.report = k; R.time = b.time(); R.kind = "BO"; R.state = blob; R.summaryCount = b.summary().size();
        rf.restarts.push_back(R);
        const std::string path = tmp + "/" + label + ".rsres";
        saveResultsFile(rf, path);
        const ResultsFile back = loadResultsFile(path);
        const ResultsFile::Restart* rb = back.restartAt(k);
        Simulator r; r.log = [](const std::string&) {};
        const Case c2 = projectFromText(back.caseText);
        r.setup(c2); r.initialize();
        bool okLoad = rb != nullptr;
        if (okLoad) r.loadState(rb->state);
        while (okLoad && r.runReportStep()) {}
        const SummaryPoint& A = a.summary().back();
        const SummaryPoint B = okLoad && !r.summary().empty() ? r.summary().back() : SummaryPoint();
        const double d = std::max({relDiff(A.fopt, B.fopt), relDiff(A.fwpt, B.fwpt), relDiff(A.fgpt, B.fgpt), relDiff(A.fpr, B.fpr), relDiff(A.time, B.time)});
        const double dp = okLoad && !r.snapshots().empty() ? maxRelArr(a.snapshots().back().pressure, r.snapshots().back().pressure) : 1.0;
        const bool sumOk = back.summary.size() == rf.summary.size() && back.snaps.size() == rf.snaps.size() && back.caseText == rf.caseText &&
                           (back.summary.empty() || (back.summary.back().fopt == rf.summary.back().fopt && back.summary.back().wells.size() == rf.summary.back().wells.size()));
        std::printf("  %-10s restart at step %d / %d: state %zu kB, |d summary| %.1e, |d pressure| %.1e, steps %d+%d vs %d\n", label, k, ns,
                    blob.size() / 1024, d, dp, b.stats().timesteps, r.stats().timesteps, a.stats().timesteps);
        char msg[160]; std::snprintf(msg, sizeof msg, "black oil %s: restart from the .rsres state == continuous run (rel. diff <= 1e-10), file round trip", label);
        check(okLoad && sumOk && d <= 1e-10 && dp <= 1e-10, msg);
      } catch (const std::exception& e) { std::printf("  exception: %s\n", e.what()); check(false, "boCase threw"); }
    };
    boCase("spe1", makeExample("spe1"));
    boCase("polymer", makeExample("polymer"));
    boCase("coning", makeExample("coning"));
    { Case c = makeExample("gascap"); c.props.hyst.on = true; c.props.hyst.model = 0; c.props.hyst.what = "KR"; boCase("gascap-hys", c); }
    boCase("condensate", makeExample("condensate"));
    boCase("spe9like", makeExample("spe9like"));
    boCase("dualporo", makeExample("dualporo"));
    {   // mismatched model is rejected
        Simulator b; b.log = [](const std::string&) {}; b.setup(makeExample("spe1")); b.initialize(); b.runReportStep();
        const std::string blob = b.saveState();
        Simulator r; r.log = [](const std::string&) {}; r.setup(makeExample("qfs")); r.initialize();
        bool threw = false; try { r.loadState(blob); } catch (const std::exception&) { threw = true; }
        check(threw, "restart state of another model is rejected with an error");
    }
    // compositional
    auto compCase = [&](const char* id, int nmax) {
      try {
        const Case c = makeExample(id);
        comp::DeckModel m; std::string err;
        if (!comp::deckModel(c, m, &err)) { check(false, (std::string("comp deck model ") + id + ": " + err).c_str()); return; }
        const int ns = std::min<int>(nmax, int(m.stepDays.size())), k = std::max(1, ns / 2);
        std::vector<comp::CompSummary> A, B, Rr;
        std::string blob;
        DeckRestart save; save.saveSteps = {k}; save.onSave = [&](int, double, const std::string& b) { blob = b; };
        const bool okA = comp::runDeckModel(m, A, {}, &err, ns);
        const bool okB = comp::runDeckModel(m, B, {}, &err, k, nullptr, &save);
        DeckRestart from; from.startStep = k; from.state = blob;
        const bool okR = !blob.empty() && comp::runDeckModel(m, Rr, {}, &err, ns, nullptr, &from);
        double d = 1;
        if (okA && okR && !A.empty() && !Rr.empty()) {
            const auto& a = A.back(); const auto& r = Rr.back();
            d = std::max({relDiff(a.fopt, r.fopt), relDiff(a.fgpt, r.fgpt), relDiff(a.fwpt, r.fwpt), relDiff(a.fgit, r.fgit), relDiff(a.fpr, r.fpr), relDiff(a.time, r.time),
                          relDiff(a.co2Diss, r.co2Diss), relDiff(a.co2Mob, r.co2Mob)});
            for (size_t w = 0; w < a.wells.size() && w < r.wells.size(); ++w) d = std::max({d, relDiff(a.wells[w].oil, r.wells[w].oil), relDiff(a.wells[w].bhp, r.wells[w].bhp), relDiff(a.wells[w].ginj, r.wells[w].ginj)});
        }
        std::printf("  comp %-8s restart at step %d / %d: blob %zu kB, |d| %.1e %s\n", id, k, ns, blob.size() / 1024, d, err.c_str());
        char msg[160]; std::snprintf(msg, sizeof msg, "compositional %s: restart == continuous run (rel. diff <= 1e-10)", id);
        check(okA && okB && okR && d <= 1e-10, msg);
      } catch (const std::exception& e) { std::printf("  exception: %s\n", e.what()); check(false, "compCase threw"); }
    };
    compCase("comp", 6);
    compCase("co2trap", 6);
    // thermal (incl. EDFM fractures with geomechanics)
    auto thermCase = [&](const char* id, int nmax) {
      try {
        const Case c = makeExample(id);
        thermal::DeckModel m; std::string err;
        if (!thermal::deckModel(c, m, &err)) { check(false, (std::string("thermal deck model ") + id + ": " + err).c_str()); return; }
        const int ns = std::min<int>(nmax, int(m.stepDays.size())), k = std::max(1, ns / 2);
        std::vector<thermal::Summary> A, B, Rr;
        std::string blob;
        DeckRestart save; save.saveSteps = {k}; save.onSave = [&](int, double, const std::string& b) { blob = b; };
        const bool okA = thermal::runDeckModel(m, A, {}, &err, ns);
        const bool okB = thermal::runDeckModel(m, B, {}, &err, k, nullptr, &save);
        DeckRestart from; from.startStep = k; from.state = blob;
        const bool okR = !blob.empty() && thermal::runDeckModel(m, Rr, {}, &err, ns, nullptr, &from);
        double d = 1;
        if (okA && okR && !A.empty() && !Rr.empty()) {
            const auto& a = A.back(); const auto& r = Rr.back();
            d = std::max({relDiff(a.fopt, r.fopt), relDiff(a.fwpt, r.fwpt), relDiff(a.fwit, r.fwit), relDiff(a.fpr, r.fpr), relDiff(a.favgT, r.favgT), relDiff(a.fhlt, r.fhlt), relDiff(a.time, r.time)});
            for (size_t w = 0; w < a.wells.size() && w < r.wells.size(); ++w) d = std::max({d, relDiff(a.wells[w].tfluid, r.wells[w].tfluid), relDiff(a.wells[w].bhp, r.wells[w].bhp), relDiff(a.wells[w].hrate, r.wells[w].hrate)});
        }
        std::printf("  thermal %-14s restart at step %d / %d: blob %zu kB, |d| %.1e %s\n", id, k, ns, blob.size() / 1024, d, err.c_str());
        char msg[160]; std::snprintf(msg, sizeof msg, "thermal %s: restart == continuous run (rel. diff <= 1e-10)", id);
        check(okA && okB && okR && d <= 1e-10, msg);
      } catch (const std::exception& e) { std::printf("  exception: %s\n", e.what()); check(false, "thermCase threw"); }
    };
    thermCase("css", 8);
    thermCase("egs_edfm", 6);
    thermCase("sagd_thm", 4);

    std::printf("== C14 ECLIPSE result writer (EGRID / INIT / SMSPEC / UNSMRY / UNRST) read back by the ECLIPSE reader\n");
    auto eclCase = [&](const char* id) {
      try {
        const Case c = makeExample(id);
        Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
        while (a.runReportStep()) {}
        const std::string base = tmp + "/" + id;
        std::string err; ecl::WriteReport wr;
        const bool ok = ecl::writeResultCase(base, c, a.grid(), a.summary(), a.snapshots(), a.wellNames(), &err, &wr);
        ecl::ResultCase rc;
        const bool okR = ok && ecl::openResultCase(base + ".EGRID", rc, &err);
        bool geo = okR && rc.egrid.nx == a.grid().nx && rc.egrid.ny == a.grid().ny && rc.egrid.nz == a.grid().nz && rc.egrid.numActive() == a.grid().numActive();
        double dz = 0;
        if (geo && c.grid.cornerPoint()) for (size_t i = 0; i < c.grid.zcorn.size(); ++i) dz = std::max(dz, std::fabs(rc.egrid.zcorn[i] - c.grid.zcorn[i]) / std::max(1.0, std::fabs(c.grid.zcorn[i])));
        geo = geo && dz < 1e-6;
        // geometry of the display grid vs the simulation grid (cell centres)
        double dc = 0;
        if (geo && rc.grid) for (int g = 0; g < a.grid().numCells(); ++g) {
            const Vec3 p = a.grid().center[size_t(g)], q = rc.grid->center[size_t(g)];
            dc = std::max({dc, double(std::fabs(p.x - q.x)), double(std::fabs(p.y - q.y)), double(std::fabs(p.z - q.z))});
        }
        const int iF = rc.summary.find("FOPT"), iW = a.wellNames().empty() ? -1 : rc.summary.find("WBHP:" + a.wellNames()[0]);
        const double fopt = iF >= 0 ? rc.summary.series(iF).back() : -1;
        const bool sum = iF >= 0 && relDiff(fopt, a.summary().back().fopt) < 1e-6 && rc.summary.time.size() == a.summary().size() && iW >= 0 &&
                         rc.summary.reportIndex.size() + 1 >= a.snapshots().size() - 1;
        double dp = 1; size_t nst = rc.restart.steps().size();
        if (nst == a.snapshots().size()) {
            const std::vector<float> p = rc.restart.load(int(nst) - 1, "PRESSURE");
            dp = maxRelArr(p, a.snapshots().back().pressure);
        }
        const std::vector<float> poro = rc.initArray("PORO"), pv = rc.initArray("PORV");
        const bool init = rc.hasInit && poro.size() == size_t(a.grid().numActive()) && pv.size() == poro.size() && relDiff(pv[0], a.grid().poreVolRef[0]) < 1e-6;
        std::printf("  %-8s %zu files, %zu summary vectors x %zu, %zu restart steps, |dzcorn| %.1e, |dcentre| %.2g, |dP| %.1e %s\n", id, wr.files.size(),
                    rc.summary.vars.size(), rc.summary.time.size(), nst, dz, dc, dp, err.c_str());
        char msg[160]; std::snprintf(msg, sizeof msg, "ECLIPSE writer %s: grid, INIT, summary (FOPT, WBHP, report steps) and restart PRESSURE read back", id);
        check(ok && okR && geo && dc < 0.05 && sum && dp < 1e-6 && init, msg);
        rc = ecl::ResultCase();   // close the files first (Windows cannot delete open files)
        std::error_code ec;
        for (const auto& f : wr.files) std::filesystem::remove(f, ec);
      } catch (const std::exception& e) { std::printf("  exception: %s\n", e.what()); check(false, "eclCase threw"); }
    };
    eclCase("spe1");
    eclCase("faulted");
    eclCase("dualporo");
    std::error_code ec; std::filesystem::remove_all(tmp, ec);
}

// C14 MINC: a 10 m matrix slab between two horizontal fractures (1D) drained at constant fracture pressure. The
// average matrix pressure follows the slab diffusion solution F(t) = 1 - sum 8/((2m+1)^2 pi^2) exp(-(2m+1)^2 pi^2 D t / L^2);
// MINC (nested shells) must reproduce the transient, Warren-Root (pseudo-steady, Kazemi sigma) only the late decay.
static std::string mincDeck(int shells) {
    std::string d = R"(RUNSPEC
TITLE
MINC slab
DIMENS
 1 1 2 /
OIL
WATER
METRIC
DUALPORO
TABDIMS
 2 1 20 20 /
WELLDIMS
 1 1 1 1 /
START
 1 JAN 2020 /
GRID
DX
 2*10 /
DY
 2*10 /
DZ
 2*10 /
TOPS
 2000 /
PORO
 0.2 0.0005 /
PERMX
 0.001 5000 /
PERMY
 0.001 5000 /
PERMZ
 0.001 5000 /
NODPPM
SIGMA
 0.04 /
LZ
 2*10 /
)";
    if (shells > 1) d += "-- RESSIM: MINC " + std::to_string(shells) + " 1\n";
    d += R"(PROPS
PVCDO
 200 1.0 1e-4 1.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
DENSITY
 800 1000 1 /
ROCK
 200 0 /
SWOF
 0.2 0 1 0
 1.0 1 0 0 /
 0.0 0 1 0
 1.0 1 0 0 /
REGIONS
SATNUM
 1 2 /
SOLUTION
EQUIL
 2005 200 3000 0 /
SCHEDULE
WELSPECS
 P G 1 1 2005 OIL /
/
COMPDAT
 P 1 1 1 1 OPEN 1* 1* 0.1 /
/
WCONPROD
 P OPEN BHP 5* 150 /
/
TSTEP
 20*0.25 20*1 20*4 20*10 /
END
)";
    return d;
}

static void testHydrofrac() {
    using namespace ressim;
    std::printf("== C14: hydraulic fracture propagation (PKN / KGD) vs similarity solutions, EGS stimulation sizing\n");
    auto rel = [](double a, double b) { return a / b - 1; };
    {   // PKN, no leak-off, zero toughness: Nordgren (1972)
        hf::Input in; in.model = 0; in.E = 20e9; in.nu = 0.25; in.mu = 0.01; in.height = 20; in.rate = 0.05; in.tEnd = 600; in.maxLength = 400; in.cells = 400; in.outputs = 4;
        hf::Result r; const bool ok = hf::simulate(in, r);
        const double eL = ok ? rel(r.last().L, hf::nordgrenLength(in, in.tEnd)) : 1, eW = ok ? rel(r.last().w0, hf::nordgrenWidth(in, in.tEnd)) : 1;
        const double bal = ok ? std::fabs(r.volInjected - r.volFracture) / r.volInjected : 1;
        std::printf("  PKN 600 s: L %.2f m (Nordgren %.2f, %+.2f%%), w0 %.3f mm (%.3f, %+.2f%%), mass balance %.1e, %d steps %s\n", ok ? r.last().L : 0.0,
                    hf::nordgrenLength(in, in.tEnd), 100 * eL, ok ? r.last().w0 * 1e3 : 0.0, hf::nordgrenWidth(in, in.tEnd) * 1e3, 100 * eW, bal, r.steps, r.error.c_str());
        check(ok && std::fabs(eL) < 0.03 && std::fabs(eW) < 0.04 && bal < 1e-6, "PKN viscosity-dominated growth matches Nordgren's L ~ t^4/5 and w0 ~ t^1/5 (< 3 % / 4 %), mass conserved");
    }
    {   // KGD, viscosity vertex (Adachi & Detournay 2002), grid convergence
        hf::Input in; in.model = 1; in.E = 20e9; in.nu = 0.25; in.mu = 0.01; in.height = 50; in.rate = 0.05; in.tEnd = 300; in.maxLength = 100; in.outputs = 3;
        double e[2] = {1, 1}, ew = 1; bool ok = true;
        for (int g = 0; g < 2; ++g) {
            in.cells = g == 0 ? 100 : 200;
            hf::Result r; ok = ok && hf::simulate(in, r);
            if (ok) { e[g] = rel(r.last().L, hf::kgdLength(in, in.tEnd)); if (g == 1) ew = rel(r.last().w0, hf::kgdWidth(in, in.tEnd)); }
        }
        std::printf("  KGD M-vertex 300 s (L_ref %.2f m): L error %+.2f%% (100 cells), %+.2f%% (200 cells); w0 %+.2f%%\n", hf::kgdLength(in, in.tEnd), 100 * e[0], 100 * e[1], 100 * ew);
        check(ok && std::fabs(e[1]) < 0.06 && std::fabs(e[1]) < std::fabs(e[0]) && std::fabs(ew) < 0.05,
              "KGD (displacement-discontinuity elasticity) viscosity-dominated half-length within 6 % of the M-vertex, converging with the grid");
        // toughness vertex: L = 0.9324 (E' Q t / K')^(2/3), K' = 4 sqrt(2/pi) K_IC
        in.mu = 1e-4; in.kic = 3e6; in.maxLength = 150; in.cells = 150;
        hf::Result r; const bool okK = hf::simulate(in, r);
        const double Ep = in.E / (1 - in.nu * in.nu), Kp = 4 * std::sqrt(2 / M_PI) * in.kic;
        const double Lk = 0.9324 * std::pow(Ep * in.rate / in.height * in.tEnd / Kp, 2.0 / 3.0), eK = okK ? rel(r.last().L, Lk) : 1;
        std::printf("  KGD K-vertex (K_IC 3 MPa m^0.5, 0.1 cP): L %.2f m vs %.2f m (%+.2f%%)\n", okK ? r.last().L : 0.0, Lk, 100 * eK);
        check(okK && std::fabs(eK) < 0.03, "KGD toughness-dominated growth (calibrated tip-element K_I) within 3 % of the K-vertex");
    }
    {   // PKN with dominant Carter leak-off: L -> q sqrt(t) / (pi H C_L)
        hf::Input in; in.model = 0; in.E = 20e9; in.nu = 0.25; in.mu = 0.01; in.height = 20; in.rate = 0.05; in.tEnd = 600; in.cl = 5e-4; in.maxLength = 60; in.cells = 2400; in.outputs = 3;
        hf::Result r; const bool ok = hf::simulate(in, r);
        const double Lc = 0.5 * in.rate * std::sqrt(in.tEnd) / (M_PI * in.height * in.cl), e = ok ? rel(r.last().L, Lc) : 1;
        const double bal = ok ? std::fabs(r.volInjected - r.volFracture - r.volLeak) / r.volInjected : 1;
        std::printf("  PKN Carter leak-off (C_L 5e-4 m/s^0.5): L %.2f m vs high-leak-off limit %.2f m (%+.2f%%), efficiency %.3f, balance %.1e\n", ok ? r.last().L : 0.0, Lc, 100 * e,
                    ok ? r.last().eff : 0.0, bal);
        check(ok && e > -0.02 && e < 0.05 && r.last().eff < 0.05 && bal < 1e-9, "PKN with dominant leak-off approaches the Carter limit L = q t^1/2 / (pi H C_L); injected = stored + leaked");
    }
    {   // C17: simultaneous multi-cluster growth with stress shadowing
        hf::ClusterInput ci; ci.base.model = 0; ci.base.E = 20e9; ci.base.nu = 0.25; ci.base.mu = 0.01; ci.base.height = 20; ci.base.rate = 0.1;
        ci.base.tEnd = 600; ci.base.maxLength = 300; ci.base.cells = 300; ci.base.outputs = 4;
        auto spread = [](const std::vector<double>& v) { double a = 1e300, b = 0; for (double x : v) { a = std::min(a, x); b = std::max(b, x); } return b > 0 ? a / b : 0.0; };
        auto balance = [](const hf::ClusterResult& r) { return std::fabs(r.volInjected - r.volFracture - r.volLeak) / r.volInjected; };
        hf::ClusterResult r1, r4, r5, r5p, r5n;
        ci.clusters = 1; ci.base.maxLength = 600; ci.base.cells = 600;
        const bool ok1 = hf::simulateClusters(ci, r1);
        hf::Input one = ci.base;
        ci.base.maxLength = 300; ci.base.cells = 300;
        const double e1 = ok1 ? r1.Lf[0] / hf::nordgrenLength(one, one.tEnd) - 1 : 1;
        ci.clusters = 4; ci.spacing = 2000; const bool ok4 = hf::simulateClusters(ci, r4);
        hf::Input quarter = ci.base; quarter.rate = ci.base.rate / 4;
        double e4 = ok4 ? 0 : 1, qd = ok4 ? 0 : 1;
        if (ok4) for (int k = 0; k < 4; ++k) { e4 = std::max(e4, std::fabs(r4.Lf[size_t(k)] / hf::nordgrenLength(quarter, quarter.tEnd) - 1)); qd = std::max(qd, std::fabs(r4.qf[size_t(k)] / (ci.base.rate / 4) - 1)); }
        ci.clusters = 5; ci.spacing = 10; const bool ok5 = hf::simulateClusters(ci, r5);
        ci.shadow = false; const bool ok5n = hf::simulateClusters(ci, r5n); ci.shadow = true;
        ci.perfFriction = hf::perfFrictionCoef(1000, 4, 0.010, 0.6); const bool ok5p = hf::simulateClusters(ci, r5p);
        std::printf("  multi-cluster (PKN, Q 0.1 m3/s, H 20 m, 600 s): 1 cluster L %.1f m (Nordgren %+.2f%%); 4 far-apart clusters max |L/L_N - 1| %.3f, max |q/(Q/4) - 1| %.4f\n",
                    ok1 ? r1.Lf[0] : 0.0, 100 * e1, e4, qd);
        auto show = [](const char* n, const hf::ClusterResult& r) {
            std::printf("    %-34s L", n); for (double x : r.Lf) std::printf(" %6.1f", x);
            std::printf(" m | q share"); double Q = 0; for (double x : r.qf) Q += x; for (double x : r.qf) std::printf(" %5.3f", Q > 0 ? x / Q : 0.0);
            std::printf(" | shadow"); for (double x : r.shadowf) std::printf(" %5.2f", x * 1e-6);
            std::printf(" MPa | p_w %.2f MPa, balance %.1e %s\n", r.pw.empty() ? 0.0 : r.pw.back() * 1e-6, r.volInjected > 0 ? std::fabs(r.volInjected - r.volFracture - r.volLeak) / r.volInjected : 1.0, r.error.c_str());
        };
        if (ok5n) show("5 clusters, 10 m, no shadow", r5n);
        if (ok5) show("5 clusters, 10 m, stress shadow", r5);
        if (ok5p) show("  + limited entry (4 x 10 mm perfs)", r5p);
        check(ok1 && std::fabs(e1) < 0.03 && balance(r1) < 1e-6, "multi-cluster solver, one cluster: Nordgren half-length within 3 %, mass conserved");
        check(ok4 && e4 < 0.03 && qd < 0.01, "multi-cluster: far-apart identical clusters share the rate equally and each follows Nordgren (Q / 4)");
        check(ok5 && ok5n && spread(r5n.Lf) > 0.99 && r5.Lf[2] < r5.Lf[0] && r5.Lf[2] < r5n.Lf[2] && r5.qf[0] > r5.qf[2] && balance(r5) < 1e-6,
              "stress shadowing: interior clusters take less fluid and stay shorter than the outer ones (no shadow: even growth), mass conserved");
        check(ok5p && spread(r5p.Lf) > spread(r5.Lf) + 0.05 && balance(r5p) < 1e-6, "limited-entry perforation friction evens out the cluster growth under stress shadowing");
    }
    {   // EGS designer with stimulation sizing
        EgsDesign d; d.stages = 3; d.cell = 30; d.stimulate = true; d.stimRate = 0.05; d.stimTime = 3600; d.stimCl = 1e-5; d.stimRetain = 0.1;
        hf::Result r; const bool okS = egsStimulation(d, r);
        bool ok = okS; Case c; std::string err;
        try { c = makeEgsCase(d); } catch (const std::exception& ex) { ok = false; err = ex.what(); }
        thermal::DeckModel m; ok = ok && thermal::deckModel(c, m, &err);
        int nwc = 0; if (ok) for (const auto& w : m.base.wells) for (int q : w.cells) if (q >= m.base.nx * m.base.ny * m.base.nz) ++nwc;
        const double L = okS ? r.last().L : 0, wm = okS ? r.meanWidth() : 0;
        bool geom = ok && c.thermal.fractures.size() == 3;
        if (geom) for (const auto& f : c.thermal.fractures) geom = geom && std::fabs((f.y1 - f.y0) - 2 * L) < 1e-6 && std::fabs(f.aperture - d.stimRetain * wm) < 1e-12;
        std::printf("  EGS stimulation (PKN, 0.05 m3/s x 3600 s, H %.0f m): half-length %.1f m, mean width %.2f mm, net pressure %.2f MPa, efficiency %.3f -> aperture %.3f mm, "
                    "%d well-fracture connections %s\n", d.fracHeight, L, wm * 1e3, okS ? r.last().pnet0 * 1e-6 : 0.0, okS ? r.last().eff : 0.0, d.stimRetain * wm * 1e3, nwc, err.c_str());
        check(ok && geom && L > d.wellSep && nwc == 6, "EGS designer sizes every stage from the propagation simulation (plane = injector +- L, aperture = retained mean width) and connects both wells");
    }
}

static void testMinc() {
    std::printf("== C14 MINC (nested matrix shells) vs the slab diffusion solution\n");
    auto run = [&](int shells, double& errMax, double& errLate, int& ncell, double& tLate) -> bool {
        try {
            std::vector<std::string> warn;
            Case c = eclipseDeckFromString(mincDeck(shells), "", &warn);
            if (c.grid.minc != (shells > 1 ? shells : 0)) { std::printf("    MINC keyword not imported (%d)\n", c.grid.minc); return false; }
            Simulator sim; sim.log = [](const std::string&) {};
            sim.setup(c); sim.initialize();
            while (sim.runReportStep()) {}
            const Grid& G = sim.grid();
            ncell = G.numActive();
            // diffusivity of the oil-filled matrix (Swc water immobile): ct = So co + Sw cw
            const double ct = 0.8 * 1e-4 + 0.2 * 4e-5, D = 0.00852702 * 0.001 / (0.2 * 1.0 * ct), L = 10.0;
            errMax = 0; errLate = 0; tLate = 0;
            const auto& S = sim.snapshots();
            for (size_t k = 1; k < S.size(); ++k) {
                double pv = 0, pp = 0, pf = 0;
                for (int a = 0; a < G.numActive(); ++a) {
                    const int g = G.globalOf[size_t(a)];
                    if (G.isFracture(g)) { pf = S[k].pressure[size_t(a)]; continue; }
                    pv += G.poreVolRef[size_t(a)]; pp += G.poreVolRef[size_t(a)] * S[k].pressure[size_t(a)];
                }
                const double pm = pp / pv, t = S[k].time;
                const double Fs = (200.0 - pm) / (200.0 - pf);
                double Fa = 1;
                for (int m = 0; m < 200; ++m) { const double o = 2 * m + 1; Fa -= 8 / (o * o * M_PI * M_PI) * std::exp(-o * o * M_PI * M_PI * D * t / (L * L)); }
                errMax = std::max(errMax, std::fabs(Fs - Fa));
                if (Fa > 0.5) { errLate = std::max(errLate, std::fabs(Fs - Fa)); if (tLate == 0) tLate = t; }
                if (shells > 1 && (k == 1 || k == 5 || k == 20 || k == 40 || k == 60 || k == S.size() - 1))
                    std::printf("    t %6.2f d  F MINC %.4f  analytic %.4f  (fracture %.2f bar)\n", t, Fs, Fa, pf);
            }
            return true;
        } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); return false; }
    };
    double eW = 1, eWl = 1, eM = 1, eMl = 1, eM4 = 1, eM4l = 1, t0 = 0; int nW = 0, nM = 0, nM4 = 0;
    const bool okW = run(0, eW, eWl, nW, t0), okM4 = run(4, eM4, eM4l, nM4, t0), okM = run(10, eM, eMl, nM, t0);
    std::printf("    max |F - F_exact|: Warren-Root %.4f (late %.4f, %d cells) | MINC 4 shells %.4f (late %.4f, %d cells) | MINC 10 shells %.4f (late %.4f, %d cells)\n",
                eW, eWl, nW, eM4, eM4l, nM4, eM, eMl, nM);
    check(okW && okM && okM4 && nM == 11 && nM4 == 5, "MINC: shells appended as active cells (fracture + 10 / 4 nested shells)");
    check(okM && eM < 0.02 && eM < 0.5 * eW, "MINC 10 shells: slab drainage transient within 0.02 of the analytic solution (better than Warren-Root)");
    check(okM4 && eM4 < eW && eM <= eM4 + 1e-3, "MINC converges with the number of shells (4 -> 10)");
    // deck round trip of the extension line
    try {
        Case c = eclipseDeckFromString(mincDeck(6), "");
        c.grid.mincVol = {0.02, 0.05, 0.1, 0.18, 0.25, 0.4};
        const Case c2 = eclipseDeckFromString(eclipseDeckToString(c), "");
        check(c2.grid.minc == 6 && c2.grid.mincDims == 1 && c2.grid.mincVol.size() == 6 && std::fabs(c2.grid.mincVol[5] - 0.4) < 1e-12, "MINC: deck export / import round trip (-- RESSIM: MINC)");
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "MINC deck round trip"); }
}

// C17 MINC in the thermal engine: heat release of a 10 m matrix slab (fractures on both faces, held at 30 C by a huge
// fracture heat capacity) - MINC shells vs Warren-Root vs an explicit fine 1D model of the half slab (the reference:
// same engine and constitutive laws, 20 conduction cells), plus the slab conduction series with the fine model's
// diffusivity fitted from its late decay (the series' shape is then a parameter-free check of the transient).
// mode 0 Warren-Root, 1 MINC (shells), 2 fine single-porosity half slab
static Case thermSlabCase(int mode, int shells) {
    Case c = makeExample("steamflood");
    c.wells.clear();
    c.thermal.winj.clear(); c.thermal.traps.clear();
    std::vector<ReportStep> st;
    for (int i = 0; i < 60; ++i) { ReportStep r; r.length = i < 20 ? 1.0 : i < 40 ? 5.0 : 20.0; st.push_back(r); }
    c.schedule = st;
    GridInput g;
    const double hcR = 2347, lam = 150, Tm = 90, Tf = 30, kM = 1e-3;
    std::vector<double> hc, th, t0;
    if (mode < 2) {
        g.nx = g.ny = 1; g.nz = 2; g.dual = 1;
        g.dx = {10, 10}; g.dy = {10, 10}; g.dz = {10, 10}; g.tops = {300, 300};
        g.poro = {0.3, 0.001}; g.permx = {kM, 5000}; g.permy = g.permx; g.permz = g.permx; g.nodppm = true;
        g.lz = {10, 10}; g.ltosigma = true;   // one fracture set (Z), spacing 10 m
        g.ltoF[0] = g.ltoF[1] = 0; g.ltoF[2] = M_PI * M_PI;   // pseudo-steady slab shape factor sigma = pi^2 / L^2 (Lim & Aziz)
        if (mode == 1) { g.minc = shells; g.mincDims = 1; }
        hc = {hcR, 1e9}; th = {lam, lam}; t0 = {Tm, Tf};
        g.multpv = {1.0, 1e6};   // the fracture system is a constant-pressure / constant-temperature boundary
    } else {
        const int nzS = 20;   // half slab 5 m: 20 x 0.25 m below a thin fracture cell
        g.nx = g.ny = 1; g.nz = nzS + 1;
        g.dx.assign(size_t(g.nz), 10); g.dy.assign(size_t(g.nz), 10); g.dz.assign(size_t(g.nz), 0.25); g.dz[0] = 0.002;
        g.tops.assign(1, 300); double z = 300; for (int k = 1; k < g.nz; ++k) { z += g.dz[size_t(k - 1)]; g.tops.push_back(z); }
        g.poro.assign(size_t(g.nz), 0.3); g.poro[0] = 0.3; g.permx.assign(size_t(g.nz), kM); g.permx[0] = 5000; g.permy = g.permx; g.permz = g.permx;
        g.multpv.assign(size_t(g.nz), 1.0); g.multpv[0] = 1e9;
        hc.assign(size_t(g.nz), hcR); hc[0] = 1e12; th.assign(size_t(g.nz), lam); t0.assign(size_t(g.nz), Tm); t0[0] = Tf;
    }
    c.grid = g;
    c.thermal.heatcr = hc; c.thermal.heatcrt.clear(); c.thermal.thconr = th; c.thermal.thconsf.clear(); c.thermal.tempi = t0; c.thermal.tempvd.clear();
    c.thermal.xmf.clear(); c.thermal.ymf.clear();
    c.thermal.heatLoss = false;
    return c;
}

static void testThermMinc() {
    std::printf("== C17 MINC in the thermal engine: matrix-slab heat release (MINC vs Warren-Root vs fine 1D model)\n");
    const double Tm = 90, Tf = 30;
    // F(t) = heat-release fraction of the matrix, bulk-heat-capacity weighted (T0 - T_avg) / (T0 - T_f)
    auto run = [&](int mode, int shells, std::vector<double>& tt, std::vector<double>& F, int& ncell, std::string& err) -> bool {
        tt.clear(); F.clear();
        thermal::DeckModel m; std::vector<thermal::Summary> out;
        const Case c = thermSlabCase(mode, shells);
        if (!thermal::deckModel(c, m, &err)) return false;
        ncell = int(m.base.cellPV.size());
        const int nGrid = c.grid.nx * c.grid.ny * c.grid.nz;
        return thermal::runDeckModel(m, out, [&](const thermal::Summary& S, const thermal::Simulator& sim) {
            const auto& T = sim.temperature();
            double w = 0, wt = 0;
            for (size_t a = 0; a < T.size(); ++a) {
                const int gq = a < m.cellGlobal.size() ? m.cellGlobal[a] : -1;
                const bool frac = mode < 2 ? (gq >= nGrid / 2 && gq < nGrid) : gq == 0;
                if (frac) continue;
                const double V = m.base.cellBulk[a];
                w += V; wt += V * (T[a] - 273.15);
            }
            if (S.time > 0 && w > 0) { tt.push_back(S.time); F.push_back((Tm - wt / w) / (Tm - Tf)); }
            return true;
        }, &err);
    };
    std::vector<double> tR, FR, tW, FW, t4, F4, t8, F8;
    int nR = 0, nW = 0, n4 = 0, n8 = 0; std::string eR, eW, e4, e8;
    const bool okR = run(2, 0, tR, FR, nR, eR), okW = run(0, 0, tW, FW, nW, eW), ok4 = run(1, 4, t4, F4, n4, e4), ok8 = run(1, 8, t8, F8, n8, e8);
    auto err = [&](const std::vector<double>& F, const std::vector<double>& t, double& late) {
        double e = 0; late = 0;
        for (size_t k = 0; k < F.size() && k < FR.size(); ++k) {
            if (std::fabs(t[k] - tR[k]) > 1e-9) return 1.0;
            e = std::max(e, std::fabs(F[k] - FR[k]));
            if (FR[k] > 0.5) late = std::max(late, std::fabs(F[k] - FR[k]));
        }
        return F.size() == FR.size() && !F.empty() ? e : 1.0;
    };
    double lW, l4, l8;
    const double dW = err(FW, tW, lW), d4 = err(F4, t4, l4), d8 = err(F8, t8, l8);
    // diffusivity of the fine model from its late exponential decay: ln(1 - F) slope = -pi^2 alpha / L^2 (L = 10 m slab)
    double alpha = 0, eS = 1;
    if (FR.size() > 50) {
        const size_t k1 = FR.size() - 12, k2 = FR.size() - 2;
        const double s = (std::log(1 - FR[k2]) - std::log(1 - FR[k1])) / (tR[k2] - tR[k1]);
        alpha = -s * 100 / (M_PI * M_PI);
        eS = 0;
        for (size_t k = 0; k < FR.size(); ++k) {
            double Fa = 1;
            for (int j = 0; j < 400; ++j) { const double o = 2 * j + 1; Fa -= 8 / (o * o * M_PI * M_PI) * std::exp(-o * o * M_PI * M_PI * alpha * tR[k] / 100); }
            eS = std::max(eS, std::fabs(FR[k] - Fa));
        }
    }
    for (size_t k : {size_t(0), size_t(4), size_t(19), size_t(39), size_t(59)})
        if (k < FR.size() && k < F8.size() && k < FW.size() && k < F4.size())
            std::printf("    t %6.1f d  F fine %.4f | Warren-Root %.4f | MINC 4 %.4f | MINC 8 %.4f\n", tR[k], FR[k], FW[k], F4[k], F8[k]);
    std::printf("    max |F - F_fine|: Warren-Root %.4f (late %.4f, %d cells) | MINC 4 %.4f (late %.4f, %d cells) | MINC 8 %.4f (late %.4f, %d cells)\n",
                dW, lW, nW, d4, l4, n4, d8, l8, n8);
    std::printf("    fine model: alpha %.4f m2/d from the late decay, max |F - slab series| %.4f %s%s%s%s\n", alpha, eS, eR.c_str(), eW.c_str(), e4.c_str(), e8.c_str());
    check(okR && okW && ok4 && ok8 && n4 == 5 && n8 == 9, "thermal MINC: deck model builds the nested shells (fracture + 4 / 8 shells) and runs");
    check(okR && eS < 0.02, "thermal fine 1D slab follows the slab conduction series (diffusivity fitted from the late decay) within 0.02");
    check(ok8 && d8 < 0.015 && d8 < 0.25 * dW, "thermal MINC 8 shells: matrix heat release within 0.015 of the fine model (Warren-Root misses the transient)");
    check(ok4 && ok8 && d8 <= d4 + 1e-3, "thermal MINC converges with the number of shells (4 -> 8)");
    {   // conduction geometry of layers with different DZ and permeability: G = A / (dz1 / 2 + dz2 / 2) (C17 review fix; the old
        // T / k_harmonic gave A / 9.9 instead of A / 5.5 here)
        Case c = thermSlabCase(2, 0);
        GridInput& g = c.grid;
        g.nz = 2; g.dx = {10, 10}; g.dy = {10, 10}; g.dz = {1, 10}; g.tops = {300, 301}; g.poro = {0.3, 0.3};
        g.permx = {100, 1}; g.permy = g.permx; g.permz = g.permx; g.ntg = {0.5, 0.5}; g.multpv.clear();
        c.thermal.heatcr = {2347}; c.thermal.thconr = {150}; c.thermal.tempi = {90};
        thermal::DeckModel m; std::string e;
        const bool ok = thermal::deckModel(c, m, &e) && m.base.conns.size() == 1;
        const double G = ok ? m.base.conns[0].G : 0;
        std::printf("    layered conduction geometry: G %.4f m (exact A / d = %.4f) %s\n", G, 100 / 5.5, e.c_str());
        check(ok && std::fabs(G - 100 / 5.5) < 1e-9 * G, "thermal conduction geometry A / d is exact for layers of different thickness and permeability (NTG excluded)");
    }
}

// C17 MINC in the compositional engine: pressure drainage of a 10 m (32.8 ft) matrix slab, undersaturated SPE5 oil,
// fracture held at 3500 psia (huge pore volume) - MINC vs Warren-Root vs a fine 1D model of the half slab (X direction)
static Case compSlabCase(int mode, int shells) {
    Case c = makeExample("comp");
    c.wells.clear();
    std::vector<ReportStep> st;
    for (int i = 0; i < 60; ++i) { ReportStep r; r.length = i < 20 ? 0.1 : i < 40 ? 0.5 : 2.0; st.push_back(r); }
    c.schedule = st;
    GridInput g;
    const double L = 32.8084, kM = 1e-3, p0 = 4000, pf = 3500;
    std::vector<double> pi, sw;
    if (mode < 2) {
        g.nx = g.ny = 1; g.nz = 2; g.dual = 1;
        g.dx = {L, L}; g.dy = {L, L}; g.dz = {L, L}; g.tops = {8400, 8400};
        g.poro = {0.3, 0.001}; g.permx = {kM, 5000}; g.permy = g.permx; g.permz = g.permx; g.nodppm = true;
        g.lz = {L, L}; g.ltosigma = true; g.ltoF[0] = g.ltoF[1] = 0; g.ltoF[2] = M_PI * M_PI;
        if (mode == 1) { g.minc = shells; g.mincDims = 1; }
        g.multpv = {1.0, 1e6};
        pi = {p0, pf};
    } else {
        const int nS = 20;
        g.nx = nS + 1; g.ny = g.nz = 1;
        g.dx.assign(size_t(g.nx), 0.5 * L / nS); g.dx[0] = 0.0066; g.dy.assign(size_t(g.nx), L); g.dz.assign(size_t(g.nx), L); g.tops.assign(size_t(g.nx), 8400);
        g.poro.assign(size_t(g.nx), 0.3); g.permx.assign(size_t(g.nx), kM); g.permx[0] = 5000; g.permy = g.permx; g.permz = g.permx;
        g.multpv.assign(size_t(g.nx), 1.0); g.multpv[0] = 1e9;
        pi.assign(size_t(g.nx), p0); pi[0] = pf;
    }
    c.grid = g;
    c.init = InitialState();
    c.init.pressure = pi; c.init.swat.assign(pi.size(), 0.2);
    return c;
}

static void testCompMinc() {
    std::printf("== C17 MINC in the compositional engine: matrix-slab pressure drainage (MINC vs Warren-Root vs fine 1D model)\n");
    auto run = [&](int mode, int shells, std::vector<double>& tt, std::vector<double>& F, int& ncell, std::string& err) -> bool {
        tt.clear(); F.clear();
        comp::DeckModel m; std::vector<comp::CompSummary> out;
        const Case c = compSlabCase(mode, shells);
        if (!comp::deckModel(c, m, &err)) return false;
        ncell = int(m.base.cellPV.size());
        const int nGrid = c.grid.nx * c.grid.ny * c.grid.nz;
        const double p0 = 4000 * 6894.757, pf = 3500 * 6894.757;
        return comp::runDeckModel(m, out, [&](const comp::CompSummary& S, const comp::CompSimulator& sim) {
            const auto& P = sim.pressure();
            double v = 0, vp = 0;
            for (size_t a = 0; a < P.size(); ++a) {
                const int gq = a < m.cellGlobal.size() ? m.cellGlobal[a] : -1;
                if (mode < 2 ? (gq >= nGrid / 2 && gq < nGrid) : gq == 0) continue;
                v += m.base.cellPV[a]; vp += m.base.cellPV[a] * P[a];
            }
            if (S.time > 0 && v > 0) { tt.push_back(S.time); F.push_back((p0 - vp / v) / (p0 - pf)); }
            return true;
        }, &err);
    };
    std::vector<double> tR, FR, tW, FW, t4, F4, t8, F8;
    int nR = 0, nW = 0, n4 = 0, n8 = 0; std::string eR, eW, e4, e8;
    const bool okR = run(2, 0, tR, FR, nR, eR), okW = run(0, 0, tW, FW, nW, eW), ok4 = run(1, 4, t4, F4, n4, e4), ok8 = run(1, 8, t8, F8, n8, e8);
    auto err = [&](const std::vector<double>& F, const std::vector<double>& t) {
        if (F.size() != FR.size() || F.empty()) return 1.0;
        double e = 0;
        for (size_t k = 0; k < F.size(); ++k) { if (std::fabs(t[k] - tR[k]) > 1e-9) return 1.0; e = std::max(e, std::fabs(F[k] - FR[k])); }
        return e;
    };
    const double dW = err(FW, tW), d4 = err(F4, t4), d8 = err(F8, t8);
    for (size_t k : {size_t(0), size_t(9), size_t(19), size_t(39), size_t(59)})
        if (k < FR.size() && k < F8.size() && k < FW.size() && k < F4.size())
            std::printf("    t %6.2f d  F fine %.4f | Warren-Root %.4f | MINC 4 %.4f | MINC 8 %.4f\n", tR[k], FR[k], FW[k], F4[k], F8[k]);
    std::printf("    max |F - F_fine|: Warren-Root %.4f (%d cells) | MINC 4 %.4f (%d cells) | MINC 8 %.4f (%d cells) %s%s%s%s\n", dW, nW, d4, n4, d8, n8,
                eR.c_str(), eW.c_str(), e4.c_str(), e8.c_str());
    check(okR && okW && ok4 && ok8 && n4 == 5 && n8 == 9 && !FR.empty() && FR.back() > 0.9, "compositional MINC: shells built, all four models drain the matrix");
    check(ok8 && d8 < 0.015 && d8 < 0.25 * dW, "compositional MINC 8 shells: matrix drainage within 0.015 of the fine model (Warren-Root misses the transient)");
    check(ok4 && ok8 && d8 <= d4 + 1e-3, "compositional MINC converges with the number of shells (4 -> 8)");
}

// C17 mechanical dispersion (DISPERC): 1D miscible gas displacement (C1-C3, supercritical single phase). The front variance
// grows by sigma^2 = 2 alpha x_f on top of the upwind numerical dispersion (convection-dispersion, D = alpha v)
// C17 review: LGR around a producer in the compositional and thermal deck paths (well indices of LGR completions used to read
// PERMX past its end; LGR cells took default porosity / SATNUM)
static void testDeckLgr() {
    std::printf("== C17 LGR cells in the compositional and thermal deck builders\n");
    {
        Case base = makeExample("comp");
        base.schedule.resize(8);
        Case c = base;
        c.grid.lgrs.push_back(makeLgrAroundWell(c, 0, 1, 3, 1, true));
        moveWellToLgr(c, 0, 0);
        auto run = [&](const Case& cs, double& fopt, double& wi, std::string& err) {
            comp::DeckModel m; std::vector<comp::CompSummary> out;
            if (!comp::deckModel(cs, m, &err)) return false;
            wi = 0;
            for (const auto& sw : m.stepWells) for (const auto& w : sw) if (w.name == "PROD") for (double x : w.wi) wi += x;
            const bool ok = comp::runDeckModel(m, out, {}, &err);
            fopt = out.empty() ? 0 : out.back().fopt;
            return ok;
        };
        double f0 = 0, f1 = 0, w0 = 0, w1 = 0; std::string e0, e1;
        const bool ok0 = run(base, f0, w0, e0), ok1 = run(c, f1, w1, e1);
        std::printf("    compositional: FOPT %.6g coarse vs %.6g with a 3x3 LGR around PROD; sum WI %.4g vs %.4g (per step) %s%s\n", f0, f1, w0, w1, e0.c_str(), e1.c_str());
        check(ok0 && ok1 && std::isfinite(w1) && w1 > 0 && std::fabs(f1 / f0 - 1) < 0.05, "compositional deck: producer completed in an LGR has finite well indices, FOPT within 5 % of the coarse run");
    }
    {
        Case base = makeExample("steamflood");
        base.schedule.resize(std::min<size_t>(base.schedule.size(), 6));
        Case c = base;
        c.grid.lgrs.push_back(makeLgrAroundWell(c, 0, 1, 3, 1, true));
        moveWellToLgr(c, 0, 0);
        auto run = [&](const Case& cs, double& cum, std::string& err) {
            thermal::DeckModel m; std::vector<thermal::Summary> out;
            if (!thermal::deckModel(cs, m, &err)) return false;
            const bool ok = thermal::runDeckModel(m, out, {}, &err);
            cum = out.empty() ? 0 : out.back().fopt + out.back().fwpt;
            return ok;
        };
        double f0 = 0, f1 = 0; std::string e0, e1;
        const bool ok0 = run(base, f0, e0), ok1 = run(c, f1, e1);
        std::printf("    thermal: liquid produced %.6g coarse vs %.6g with a 3x3 LGR around well 1 %s%s\n", f0, f1, e0.c_str(), e1.c_str());
        check(ok0 && ok1 && f0 > 0 && std::fabs(f1 / f0 - 1) < 0.10, "thermal deck: well completed in an LGR, cumulative liquid within 10 % of the coarse run");
    }
}

static void testCompDispersion() {
    using namespace ressim::eos;
    using namespace ressim::comp;
    std::printf("== C17 mechanical dispersion (DISPERC) in the compositional engine\n");
    const double T = fahrenheitToK(160);
    const std::vector<double> z0 = {0.99, 0.01, 0, 0, 0, 0}, zi = {0.80, 0.20, 0, 0, 0, 0};
    const int nx = 200; const double dx = 0.5, phi = 0.2, tEnd = 50 * 86400.0;
    auto run = [&](double alpha, double& xf, double& var, double& mb) {
        CompCase c; c.nx = nx; c.ny = 1; c.nz = 1; c.dx = {dx}; c.dy = {1}; c.dz = {1}; c.poro = {phi};
        c.permx = c.permy = c.permz = {1e-11}; c.gravity = false; c.fluid = spe5Fluid(); c.T = T;
        c.p0 = {psiaToPa(4000)}; c.z0 = {z0}; c.dispersivity = alpha;
        const TParams tp = atTemperature(c.fluid, T);
        EvalOptions vo; vo.dn = false; vo.dP = false;
        const double vg = evalPhase(c.fluid, tp, psiaToPa(4000), zi, vo).vm;
        const double pv = nx * dx * phi, tPV = 100 * 86400.0;
        CompWell inj; inj.injector = true; inj.cells = {0}; inj.wi = {1}; inj.rate = pv / vg / tPV; inj.zinj = zi;
        CompWell prd; prd.cells = {nx - 1}; prd.wi = {1e-11}; prd.bhp = psiaToPa(4000) - 1e5;
        prd.control = 2; prd.target = pv / tPV;
        c.wells = {inj, prd};
        CompSimulator sim; sim.setup(c);
        std::vector<double> in0(6); for (int i = 0; i < 6; ++i) in0[i] = sim.componentInPlace(i);
        const bool ok = sim.advanceTo(tEnd, 0.25 * 86400.0);
        const auto& N = sim.moles(); const int ncw = sim.ncw();
        std::vector<double> Cn(static_cast<size_t>(nx));
        for (int q = 0; q < nx; ++q) { double t = 0; for (int i = 0; i < 6; ++i) t += N[size_t(q) * ncw + i]; Cn[size_t(q)] = (N[size_t(q) * ncw + 1] / t - 0.01) / 0.19; }
        // moments of -dC/dx (front distribution); cell faces at x = (q+1) dx
        double m0 = 0, m1 = 0, m2 = 0;
        for (int q = 0; q + 1 < nx; ++q) { const double w = Cn[size_t(q)] - Cn[size_t(q + 1)], x = (q + 1) * dx; m0 += w; m1 += w * x; m2 += w * x * x; }
        xf = m1 / m0; var = m2 / m0 - xf * xf;
        mb = 0;
        for (int i = 0; i < 2; ++i) mb = std::max(mb, std::fabs(sim.componentInPlace(i) - in0[i] - sim.injected(i) + sim.produced(i)) / (in0[i] + sim.injected(i)));
        return ok;
    };
    double x0, v0, mb0, x1, v1, mb1, x2, v2, mb2;
    const bool ok0 = run(0, x0, v0, mb0), ok1 = run(1.0, x1, v1, mb1), ok2 = run(2.0, x2, v2, mb2);
    const double a1 = (v1 - v0) / (2 * x1), a2 = (v2 - v0) / (2 * x2);
    std::printf("    front at %.1f / %.1f / %.1f m; variance %.1f (numerical) / %.1f (alpha 1 m) / %.1f m2 (alpha 2 m); "
                "recovered alpha %.3f / %.3f m; balance %.1e\n", x0, x1, x2, v0, v1, v2, a1, a2, std::max(mb0, std::max(mb1, mb2)));
    check(ok0 && ok1 && ok2 && std::max(mb0, std::max(mb1, mb2)) < 1e-9, "dispersion runs: component balance exact");
    check(std::fabs(a1 / 1.0 - 1) < 0.15 && std::fabs(a2 / 2.0 - 1) < 0.15, "front variance grows by 2 alpha x_f (convection-dispersion), alpha recovered within 15 %");
    // deck path: DISPERC round trip and alpha / d on the connections
    Case dc = makeExample("comp");
    dc.grid.disperc.assign(size_t(dc.grid.size()), 3.28084);
    std::string err; Case rc;
    try { rc = eclipseDeckFromString(eclipseDeckToString(dc)); } catch (const std::exception& e) { err = e.what(); }
    comp::DeckModel m; double tmax = 1e300, tvert = 0;
    const bool okd = err.empty() && comp::deckModel(rc, m, &err);
    if (okd) for (const auto& cn : m.base.conns) if (cn.Tdisp > 0) { tmax = std::min(tmax, cn.Tdisp); tvert = std::max(tvert, cn.Tdisp); }
    std::printf("    deck round trip: %zu DISPERC values, alpha/d lateral %.6f (expected 3.28084 / 500 = %.6f), vertical max %.6f (3.28084 / 25 = %.6f) %s\n", rc.grid.disperc.size(), tmax, 3.28084 / 500, tvert, 3.28084 / 25, err.c_str());
    check(okd && rc.grid.disperc.size() == dc.grid.disperc.size() && std::fabs(tmax - 3.28084 / 500) < 1e-6 && std::fabs(tvert - 3.28084 / 25) < 1e-5,
          "DISPERC: export -> import keeps the per-cell dispersivity; face alpha / d on the connections");
}

// ------------------------------------------------------------------------------------------------ C15 hydraulic fractures
static std::string hfDeck(bool frac, bool horizontal = false, const std::string& extra = "", const std::string& ctrl = "ORAT 5 4* 10") {
    // tight oil (Sw = Swc immobile, slightly compressible oil): single-phase drawdown of a fractured vertical well
    std::vector<double> half;
    double s = 0;
    for (int k = 1; s < 600; ++k) { half.push_back(2.0 * std::pow(1.35, k)); s += half.back(); }
    std::vector<double> rows;
    for (size_t k = half.size(); k-- > 0;) rows.push_back(half[k]);
    rows.push_back(2);
    for (double v : half) rows.push_back(v);
    std::string dy;
    for (double v : rows) dy += " 81*" + std::to_string(v);
    const int ny = int(2 * half.size() + 1), jc = int(half.size()) + 1;
    const double yc = s + 1.0;
    std::string d = "RUNSPEC\nTITLE\nEDFM drawdown\nDIMENS\n 81 " + std::to_string(ny) + " 1 /\nOIL\nWATER\nMETRIC\nTABDIMS\n 1 1 20 20 /\n"
        "WELLDIMS\n 1 100 1 1 /\nSTART\n 1 JAN 2020 /\nGRID\nDX\n " + std::to_string(81 * ny) + "*10 /\nDY\n" + dy + " /\nDZ\n " + std::to_string(81 * ny) + "*10 /\nTOPS\n " +
        std::to_string(81 * ny) + "*2000 /\nPORO\n " + std::to_string(81 * ny) + "*0.1 /\nPERMX\n " + std::to_string(81 * ny) + "*0.1 /\nPERMY\n " +
        std::to_string(81 * ny) + "*0.1 /\nPERMZ\n " + std::to_string(81 * ny) + "*0.1 /\n";
    char b[256];
    if (frac && !horizontal) { std::snprintf(b, sizeof b, "-- RESSIM: HFRAC F1 305 %.6f 505 %.6f -1 -1 0.003 0 0.5 0 0\n", yc, yc); d += b; }
    d += extra;
    d += "PROPS\nPVCDO\n 300 1.0 1e-4 1.0 0 /\nPVTW\n 300 1.0 4e-5 0.5 0 /\nDENSITY\n 800 1000 1 /\nROCK\n 300 0 /\nSWOF\n 0.2 0 1 0\n 1.0 1 0 0 /\n"
         "SOLUTION\nEQUIL\n 2005 300 3000 0 /\nSCHEDULE\nWELSPECS\n P G 41 " + std::to_string(jc) + " 2005 OIL /\n/\nCOMPDAT\n";
    if (horizontal) d += " P 11 " + std::to_string(jc) + " 1 1 OPEN 1* 1* 0.1 3* X /\n P 12 " + std::to_string(jc) + " 1 1 OPEN 1* 1* 0.1 3* X /\n";   // replaced below
    else d += " P 41 " + std::to_string(jc) + " 1 1 OPEN 1* 1* 0.1 /\n";
    d += "/\nWCONPROD\n P OPEN " + ctrl + " /\n/\nTSTEP\n 10*0.1 10*0.5 10*1 10*2 10*5 /\nEND\n";
    return d;
}
// Gringarten, Ramey & Raghavan (1974) uniform-flux vertical fracture, dimensionless pressure at x_D (0.732 ~ infinite conductivity)
static double gringartenUF(double tD, double xD) {
    auto Ei = [](double x) { return std::expint(x); };
    const double a = 1 - xD, c = 1 + xD, st = std::sqrt(tD);
    return 0.5 * std::sqrt(M_PI * tD) * (std::erf(a / (2 * st)) + std::erf(c / (2 * st))) - a / 4 * Ei(-a * a / (4 * tD)) - c / 4 * Ei(-c * c / (4 * tD));
}
static void testHydraulicFractureFlow() {
    std::printf("== C15 hydraulic fractures in the black-oil grid (EDFM) vs Gringarten et al. (1974)\n");
    const double ct = 0.8 * 1e-4 + 0.2 * 4e-5, eta = 0.00852702 * 0.1 / (0.1 * 1.0 * ct), xf = 100.0;
    auto run = [&](const std::string& deck, bool cp, std::vector<double>& tD, std::vector<double>& pD, int& nHf) -> bool {
        try {
            Case c = eclipseDeckFromString(deck, "");
            if (cp) convertToCornerPoint(c.grid);
            Simulator sim; sim.log = [](const std::string&) {};
            sim.setup(c); sim.initialize();
            while (sim.runReportStep()) {}
            nHf = sim.numHydraulicFractureCells();
            for (const SummaryPoint& s : sim.summary()) {
                if (s.time <= 0 || s.wells.empty()) continue;
                const double q = s.fopr;   // B ~ 1 (pref = initial pressure)
                if (q <= 0) continue;
                tD.push_back(eta * s.time / (xf * xf));
                pD.push_back(2 * M_PI * 0.00852702 * 0.1 * 10 * (300.0 - s.wells[0].bhp) / (q * 1.0));
            }
            return !tD.empty();
        } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); return false; }
    };
    std::vector<double> tF, pF, tN, pN, tC, pC;
    int nF = 0, nN = 0, nC = 0;
    const bool okF = run(hfDeck(true), false, tF, pF, nF), okN = run(hfDeck(false), false, tN, pN, nN), okC = run(hfDeck(true), true, tC, pC, nC);
    double eMax = 0, rNo = 1e9, eCp = 0;
    for (size_t k = 0; k < tF.size(); ++k) {
        if (tF[k] < 0.02 || tF[k] > 1.0) continue;
        const double pa = gringartenUF(tF[k], 0.732);
        eMax = std::max(eMax, std::fabs(pF[k] - pa) / pa);
        if (k < pN.size()) rNo = std::min(rNo, pN[k] / pF[k]);
        if (k < pC.size()) eCp = std::max(eCp, std::fabs(pC[k] - pF[k]) / pF[k]);
        if (k % 8 == 0) std::printf("    tD %.4f  pD EDFM %.4f  Gringarten (inf. cond.) %.4f  unfractured %.4f\n", tF[k], pF[k], pa, k < pN.size() ? pN[k] : 0.0);
    }
    std::printf("    fracture cells %d (corner-point %d); max |pD - pD_exact| / pD = %.3f for 0.02 <= tD <= 1; unfractured / fractured >= %.2f; corner-point vs Cartesian %.2e\n",
                nF, nC, eMax, rNo, eCp);
    check(okF && nF == 21 && nN == 0, "EDFM: 200 m fracture (x 305 - 505) -> 21 fracture cells appended to the black-oil grid");
    check(okF && eMax < 0.06, "EDFM: fractured-well drawdown within 6% of the infinite-conductivity solution (tD 0.02 - 1)");
    check(okN && rNo > 1.5, "EDFM: the fracture lowers the drawdown by > 1.5x (vs the unfractured well)");
    check(okC && nC == nF && eCp < 0.01, "EDFM on the corner-point form of the grid matches the Cartesian result (< 1%)");
}
static void testUnconventional() {
    std::printf("== C15 unconventional: multistage layout, Langmuir adsorbed gas, stress-dependent permeability, deck round trip\n");
    // ---- multistage layout along a horizontal well
    try {
        Case c = eclipseDeckFromString(hfDeck(false, false, "", "BHP 5* 100"), "");
        WellInput& w = c.wells[0];
        w.completions.clear();
        for (int i = 10; i <= 70; ++i) { Completion cp; cp.i = i; cp.j = c.wells[0].completions.empty() ? 0 : 0; cp.k1 = cp.k2 = 0; cp.dir = 'X'; cp.diameter = 0.1; w.completions.push_back(cp); }
        const int jc = c.grid.ny / 2;   // the 2 m centre row
        for (auto& cp : w.completions) cp.j = jc;
        Case base = c;
        hf::StageLayout L; L.well = w.name; L.stages = 6; L.halfLength = 80; L.height = 10; L.width = 0.003; L.poro = 0.35;
        const auto F = hf::layoutStages(c, L);
        bool perp = F.size() == 6;
        for (const auto& f : F) perp = perp && std::fabs(f.x0 - f.x1) < 1e-6 && std::fabs(std::fabs(f.y1 - f.y0) - 160) < 1e-6;
        c.grid.hfracs = F;
        const double spacing = F.size() > 1 ? std::fabs(F[1].x0 - F[0].x0) : 0;
        std::printf("    6 stages perpendicular to the lateral, spacing %.1f m\n", spacing);
        check(perp && std::fabs(spacing - 100.0) < 1e-6, "multistage layout: 6 transverse fractures evenly spaced along a 600 m lateral");
        auto cum = [](Case cc, int& nhf) { Simulator s; s.log = [](const std::string&) {}; s.setup(cc); s.initialize(); while (s.runReportStep()) {} nhf = s.numHydraulicFractureCells(); return s.summary().back().fopt; };
        int n0 = 0, n1 = 0;
        const double q0 = cum(base, n0), q1 = cum(c, n1);
        std::printf("    cumulative oil: unstimulated lateral %.1f m3, 6-stage %.1f m3 (%d fracture cells)\n", q0, q1, n1);
        check(n1 > 6 && q1 > 1.3 * q0, "multistage hydraulic fractures raise the horizontal-well recovery (> 1.3x)");
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "multistage layout"); }
    // ---- Langmuir adsorbed gas: closed shale-gas box, depletion; material balance with the adsorbed term
    try {
        const std::string d = R"(RUNSPEC
DIMENS
 5 5 1 /
OIL
WATER
GAS
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 1 1 1 1 /
START
 1 JAN 2020 /
GRID
DX
 25*50 /
DY
 25*50 /
DZ
 25*20 /
TOPS
 25*2000 /
PORO
 25*0.06 /
PERMX
 25*0.5 /
PERMY
 25*0.5 /
PERMZ
 25*0.5 /
-- RESSIM: LANGMUIR 3.0 60 2.5
PROPS
PVDO
 50 1.05 1.0
 400 1.00 1.0 /
PVDG
 20 0.06 0.015
 50 0.024 0.016
 100 0.012 0.018
 200 0.006 0.022
 300 0.0042 0.026 /
PVTW
 200 1.0 4e-5 0.5 0 /
DENSITY
 800 1000 0.8 /
ROCK
 200 0 /
SWOF
 0.25 0 1 0
 1.0 1 0 0 /
SGOF
 0 0 1 0
 0.75 0.9 0 0 /
SOLUTION
EQUIL
 2010 200 3000 0 2500 0 /
SCHEDULE
WELSPECS
 P G 3 3 2010 GAS /
/
COMPDAT
 P 3 3 1 1 OPEN 1* 1* 0.1 /
/
WCONPROD
 P OPEN GRAT 2* 20000 2* 30 /
/
TSTEP
 20*30 /
END
)";
        Case c = eclipseDeckFromString(d, "");
        check(std::fabs(c.grid.langVL - 3.0) < 1e-12 && std::fabs(c.grid.langPL - 60) < 1e-12, "LANGMUIR deck line imported");
        Simulator sim; sim.log = [](const std::string&) {};
        sim.setup(c); sim.initialize();
        const InPlace i0 = sim.inPlace();
        while (sim.runReportStep()) {}
        const InPlace i1 = sim.inPlace();
        const double gp = sim.summary().back().fgpt;
        const double ads0 = 25 * 50 * 50 * 20 * (1 - 0.06) * 2.5 * 3.0 * 200 / (60 + 200);
        const double mb = (i0.gasFree + i0.gasAdsorbed) - (i1.gasFree + i1.gasAdsorbed);
        std::printf("    adsorbed gas initial %.4g sm3 (analytic %.4g), free %.4g; produced %.4g sm3, in-place change %.4g (desorbed %.4g)\n",
                    i0.gasAdsorbed, ads0, i0.gasFree, gp, mb, i0.gasAdsorbed - i1.gasAdsorbed);
        check(std::fabs(i0.gasAdsorbed - ads0) < 1e-6 * ads0, "Langmuir: initial adsorbed gas = V_b (1 - phi) rho V_L p / (P_L + p)");
        check(std::fabs(mb - gp) < 2e-3 * gp && i0.gasAdsorbed - i1.gasAdsorbed > 0.05 * gp, "Langmuir: production = free + desorbed gas (material balance 0.2%), desorption contributes");
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "Langmuir"); }
    // ---- stress-dependent permeability lowers the rate of a BHP well; deck round trip of the extension lines
    try {
        const std::string ex = "-- RESSIM: PERMSTRESS 0.004\n";
        auto runQ = [](const std::string& deck) { Case c = eclipseDeckFromString(deck, ""); Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize(); while (s.runReportStep()) {} return s.summary().back().fopt; };
        const double q0 = runQ(hfDeck(true, false, "", "BHP 5* 100")), q1 = runQ(hfDeck(true, false, ex, "BHP 5* 100"));
        std::printf("    cumulative oil without / with stress-dependent permeability (gamma 0.004 1/bar): %.2f / %.2f m3\n", q0, q1);
        check(q1 < q0 && q1 > 0.3 * q0, "PERMSTRESS: pressure-dependent permeability reduces the cumulative production");
        Case c = eclipseDeckFromString(hfDeck(true, false, ex + "-- RESSIM: LANGMUIR 2 50 2.6\n"), "");
        c.grid.hfracs[0].gamma = 0.01; c.grid.hfracs[0].satnum = 0;
        const Case c2 = eclipseDeckFromString(eclipseDeckToString(c), "");
        check(c2.grid.hfracs.size() == 1 && std::fabs(c2.grid.hfracs[0].x1 - 505) < 1e-9 && std::fabs(c2.grid.hfracs[0].gamma - 0.01) < 1e-12 &&
              std::fabs(c2.grid.permGamma - 0.004) < 1e-12 && std::fabs(c2.grid.langPL - 50) < 1e-12, "HFRAC / PERMSTRESS / LANGMUIR: deck export / import round trip");
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "PERMSTRESS"); }
}
static void testFaultStability() {
    std::printf("== C15 fault / fracture stability and rate-and-state slip\n");
    using namespace stab;
    // Andersonian normal-faulting regime at 2 km: Sv 46 MPa, SH 36, Sh 30 (SH along North), hydrostatic p
    StressState s = fromGradients(2000, 23e3, 18e3, 15e3, 0, 20e6, 1.0);
    Strength st; st.mu = 0.6;
    // optimally oriented normal fault strikes parallel to SH (N-S) and dips 0.5 atan(1/mu) from sigma_1 -> 60.5 deg
    const double dipOpt = 90 - 0.5 * std::atan(1 / 0.6) * 180 / M_PI;
    Plane opt{0, dipOpt}, bad{90, 30};
    const Resolved r = resolve(s, opt, st), rb = resolve(s, bad, st);
    // 2-D Mohr check: sigma_n = (s1 + s3)/2 + (s1 - s3)/2 cos 2 beta, tau = (s1 - s3)/2 sin 2 beta, beta = angle normal - sigma_1 (vertical) = dip
    const double s1 = 46e6, s3 = 30e6, beta = dipOpt * M_PI / 180;   // angle between the normal and sigma_1
    const double snA = 0.5 * (s1 + s3) + 0.5 * (s1 - s3) * std::cos(2 * beta), tauA = 0.5 * (s1 - s3) * std::sin(2 * beta);
    std::printf("    optimal plane (dip %.1f): sigma_n %.3f MPa (Mohr %.3f), tau %.3f MPa (Mohr %.3f), Ts %.3f, dP_crit %.2f MPa; misoriented Ts %.3f\n",
                dipOpt, r.sn / 1e6, snA / 1e6, r.tau / 1e6, tauA / 1e6, r.slipTend, r.dpCrit / 1e6, rb.slipTend);
    check(std::fabs(r.sn - snA) < 1e-3 && std::fabs(r.tau - tauA) < 1e-3, "stability: resolved tractions equal the Mohr-circle values");
    StressState s2 = s; s2.p += r.dpCrit;
    const Resolved r2 = resolve(s2, opt, st);
    check(std::fabs(r2.cff) < 1e-3 * r.tau && rb.slipTend < r.slipTend, "stability: the critical pore-pressure increase brings the plane exactly onto the Coulomb line");
    // poroelastic stress path: depletion raises the horizontal stress less than the pore pressure drop
    StressState sp = s; stressPath(sp, -5e6, 0.25);
    check(std::fabs((sp.Sh - s.Sh) - (-5e6 * stressPathCoeff(1.0, 0.25))) < 1e-6, "stability: poroelastic stress path dS_h = alpha (1 - 2 nu) / (1 - nu) dp");
    // Barton-Bandis: zero stress -> w_res + V_m, high stress -> w_res
    BBParams bb;
    check(std::fabs(bbAperture(bb, 0, 0) - (bb.w0 + bb.Vm)) < 1e-12 && bbAperture(bb, 1e12, 0) < bb.w0 + 1e-6 && bbAperture(bb, 20e6, 1e-3) > bbAperture(bb, 20e6, 0),
          "Barton-Bandis closure and shear dilation of the aperture");
    // rate-and-state: injection raises p on a critically stressed velocity-weakening patch -> a seismic event; velocity
    // strengthening (a > b) creeps aseismically
    RsfInput in;
    in.sn0 = r.snEff; in.tau0 = 0.55 * r.snEff;
    in.t = {0, 86400.0 * 200}; in.dTau = {0, 0}; in.dSn = {0, -0.3 * r.snEff};
    in.radius = 200;
    RsfResult o;
    const bool ok = simulateRsf(in, o);
    RsfInput in2 = in; in2.a = 0.02; in2.b = 0.015;
    RsfResult o2;
    const bool ok2 = simulateRsf(in2, o2);
    double vmax2 = 0; for (const auto& sm : o2.samples) vmax2 = std::max(vmax2, sm.V);
    std::printf("    velocity weakening: k %.3g < k_c %.3g Pa/m, %zu event(s)", o.kPatch, o.kCrit, o.events.size());
    if (!o.events.empty()) std::printf(", first at %.1f d, slip %.3f m, Vmax %.2f m/s, Mw %.2f", o.events[0].tStart / 86400, o.events[0].slip, o.events[0].Vmax, o.events[0].Mw);
    std::printf(" (%d steps)\n    velocity strengthening: %zu events, max slip rate %.2e m/s\n", o.steps, o2.events.size(), vmax2);
    check(ok && o.kPatch < o.kCrit && !o.events.empty() && o.events[0].Vmax > 0.01 && o.events[0].slip > 0, "rate-and-state: pressurised velocity-weakening patch slips seismically");
    check(ok2 && o2.events.empty() && vmax2 < 1e-3, "rate-and-state: velocity-strengthening patch creeps aseismically");
}

#include "c18_tests.inc"
#include "c20_tests.inc"
#include "c21_tests.inc"
#include "c22_tests.inc"
#include "c23_tests.inc"
#include "c24_tests.inc"
#include "c25_tests.inc"
#include "c26_tests.inc"
#include "c27_tests.inc"
#include "c28_tests.inc"
#include "c29_tests.inc"
#include "c30_tests.inc"
#include "c31_tests.inc"
#include "c32_tests.inc"
#include "c33_tests.inc"
#include "c34_tests.inc"
#include "c35_tests.inc"
#include "c36_tests.inc"
#include "c37_tests.inc"
#include "c38_tests.inc"
#include "c40_tests.inc"
#include "c41_tests.inc"
#include "c42_tests.inc"
#include "c43_tests.inc"
#include "c44_tests.inc"
#include "c45_tests.inc"
#include "c46_tests.inc"
#include "c47_tests.inc"
#include "c48_tests.inc"
#include "c49_tests.inc"

int main(int argc, char** argv) {
    const bool verbose = argc > 1 && std::string(argv[1]) == "-v";
    if (argc > 1 && std::string(argv[1]) == "-mr") { testMultRegtReopen(); testThreePointDegenerate(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-dp") { testDualPorosity(); testDualGravityDpnum(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-mc") { testComponents(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-hf") { testHydrofrac(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-hfbo") { testHydraulicFractureFlow(); testUnconventional(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-stab") { testFaultStability(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-dlgr") { testDeckLgr(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-cadj") { testCompAdjoint(); testThermAdjoint(); testMulticolourIlu(); testMpiDD(false); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-seis") { testSeismicity(); testSeismicityRobust(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29salt") { testC29Salt(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29wb") { testC29Wellbore(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29net") { testC29Network(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29pf") { testC29PfTwoWay(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29gpu") { testC29DeviceCprAmg(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30net") { testC30Network(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30hff") { testC30HfField(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30pf") { testC30PhaseFieldOutput(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30post") { testC30PostGeomech(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30cmap") { testC30Colormaps(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c31expr") { testC31Expressions(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c32geo") { testC32GeoFormats(); testC32GeoWorkflow(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c32actx") { testC32ActionConnections(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c16mp") { testMixedPrecision(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-poly") { testPolymer(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c36geo") { testC36Geo(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c36seis") { testC36Seismic(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-facil") { testC35Facilities(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-glift") { testC35GasLift(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-solv") { testC35Solvent(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c35") { testC35Chem(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c34") { testC34Items(); testEndPointScaling(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c33") { testC33Report(); testC33Api(); testC33Regress(); testC33ScheduleItems(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c36net") { testC36Network(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c49") { testC49MswOptsRoundTrip(); testC49ConvergedHeads(); testC49FaultRamp(); testC49StoppedWellBhp(); testC49CompWellSummary(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c48") { testC48Leakc(); testC48RefDepthShut(); testC48GroupUda(); testC48CellSnapshots(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c47") { testC47GasPadding(); testC47Co2storeImport(); testC47FaultCsv(); testC47ActionWellsJly(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c46") { testC46GeoProject(); testC46CompletionNumbers(); testC46SchedMultflt(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c45") { testC45LeanSolver(); testC45MswCrossflow(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c44r") { testDeck(); testCPR(); testC44WpimultFbhpdef(); testC44RestartEcon(); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c44") { testC44FaultCut(); testC44ShutInStop(); testC44WelopenLayer(); testC44WpweReference(); testC44RestartEcon(); testC44HysteresisWetting(); testC44InjectorRecirculation(); testC44WpimultFbhpdef(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c43") { testC43Mistie(); testC43Polygons(); testC43LayerCake(); testC43Seismic(); testC43GeoCapi(); testC43Wpwe(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c42") { testC42Import(); testC42Wells(); testC42Seismic2D(); testC42Calc(); testC42Cecon(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c41") { testC41Workflow(); testC41Facies(); testC41FaultNetwork(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c40geo") { testC40GeoAnalysis(); testC40StudyGeoParams(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c38prod") { testC38ProdData(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c38ilu") { testC38FloatIlu(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c38tran") { testC38TranEdit(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c38udt") { testC38Udt(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c38nonnc") { testC38Nonnc(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c37opm") { testC37Opm(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c37py") { testC37Py(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c37msw") { testC37Msw(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c36msw") { testC36Msw(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c32msw") { testC32MswUnit(); testC32MswDeck(); testC32MswAicdItems(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c31deck") { testC31Deck(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c30") { testC30Network(); testC30HfField(); testC30PhaseFieldOutput(); testC30PostGeomech(); testC30Colormaps(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c29") { testC29Salt(); testC29Wellbore(); testC29Network(); testC29DeviceCprAmg(); testC29PfTwoWay(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c28ilu") { testC28LevelIlu(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c28pf") { testC28PhaseFieldInRun(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c28da") { testC28DistAssembly(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c28") { testC28Fosm(); testC28Tracers(); testC28Geostat(); testC28GpProxy(); testC28ProdAnalysis(); testC28FaultSlipEngines(); testC28DistAssembly(); testC28LevelIlu(); testC28PhaseFieldInRun(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c27") { testC27TimestepTargets(); testC27Wecon(); testC27VanGenuchten(); testC27HmStatistics(); testC27CoupledStudy(); testC27Examples(); testC27Tracer(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c26") { testC26Kernel(); testC26Penny(); testC26Shear(); testC26RoundTrip(); testC26FaultRun(); testC26HfWidth(); testC26IterCoupling(); testC26PhaseField(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c25") { testC25HexPatch(); testC25CornerPointMesh(); testC25RoundTrip(); testC25BlackOilCoupled(); testC25Comp(); testC25Thermal(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c24") { testFaultMechStress(); testFaultMechRun(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c23") { testRsfSpringSlider(); testRsfFaultCycle(); testRsfInjection(); testRsfEmbedded(); testKernelFitAspo(); testPhaseFieldSneddonGriffith(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c22") { testScheduleGain(); testSeismicityLastInterval(); testSeismicityNearField(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-c21") { testParallelPerf(); testHfSubgrid(); testGeoCouple(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-hfall") { testEdfm(); testHydrofrac(); testHydraulicFractureFlow(); testUnconventional(); testHfGrowth(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-geoc") { testGeoCouple(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-efem") { testEfem(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-p3d") { testP3D(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-hfadj") { testHfGrowthAdjoint(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-hfg") { testHfGrowth(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-cdisp") { testCompDispersion(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-cminc") { testCompMinc(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-tminc") { testThermMinc(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-minc") { testMinc(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-rst") { testRestart(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-spe4b") { testSpe4Problem2A(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-spe11g") { testSpe11aGeometry(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-egs") { testPruessEgs(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-spe11") { testSpe11a(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-spe4") { testSpe4(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-grin") { testGringarten(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-edfm") { testEdfm(); testGringarten(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-sorb") { testSorption(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-mpida") {   // C28: mpiexec -n P run_tests -mpida (all ranks run the Newton loop)
        std::string err;
        const bool have = mpi::load(&err) && mpi::init();
        if (!have) std::printf("MPI not available: %s\n", err.c_str());
        else testC28DistAssemblyMpi();
        const int myRank = have ? mpi::rank() : 0;
        if (have) mpi::finalize();
        if (myRank > 0) return 0;
        std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0;
    }
    if (argc > 1 && std::string(argv[1]) == "-mpidd") {   // C18: mpiexec -n P run_tests -mpidd
        std::string err;
        const bool have = mpi::load(&err) && mpi::init();
        if (have && mpi::rank() > 0) { const int rc = mpi::serveDD(); mpi::finalize(); return rc; }
        testMpiDD(have && mpi::size() > 1);
        if (have) mpi::finalize();
        std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0;
    }
    if (argc > 1 && std::string(argv[1]) == "-mp") { testMixedPrecision(); testMulticolourIlu(); testMpiRuntime(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-gmp") { testGeomechPlastic(); testMandel(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-ver") { testLauwerier(); testTerzaghi(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-adjf") { testAdjointField(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-adj") { testAdjoint(); testAdjointField(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-num") { testNumerics(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-drsdt") { testStoneAndDrsdt(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-comp") { testCompositional(); testCompDeck(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-eos") { testEos(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-lgr") { testLgr(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-trap") { testTrapping(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-thm") { testThmc(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-scw") { testNearCritical(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-sol") { testSolubility(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-thermal") { testThermal(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-eor") { testEorTools(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-co2sol") { testCo2Sol(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-co2") { testCo2Store(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    if (argc > 1 && std::string(argv[1]) == "-norne") { testNornePhysics(); std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures); return failures ? 1 : 0; }
    try { testIO(); testDeck(); testEos(); testCompositional(); testCompDeck(); testCo2Store(); testTrapping(); testCo2Sol(); testEorTools(); testThermal(); testThmc(); testSolubility(); testNearCritical(); testLauwerier(); testTerzaghi(); testGeomechPlastic(); testMandel(); testSorption(); testEdfm(); testGringarten(); testSpe4(); testSpe11a(); testRestart(); testCPR(); testMixedPrecision(); testMpiRuntime(); testNumerics(); testCornerPoint(); testLgr(); testMultRegtReopen(); testThreePointDegenerate(); testNewtran(); testPhaseA(); testWetGas(); testEquilAndWellbore(); testPolymer(); testDualPorosity(); testDualGravityDpnum(); testMinc(); testThermMinc(); testCompMinc(); testCompDispersion(); testCompAdjoint(); testThermAdjoint(); testMulticolourIlu(); testMpiDD(false); testDeckLgr(); testHydrofrac(); testHydraulicFractureFlow(); testUnconventional(); testHfGrowth(); testHfGrowthAdjoint(); testP3D(); testEfem(); testParallelPerf(); testHfSubgrid(); testGeoCouple(); testScheduleGain(); testSeismicityLastInterval(); testSeismicityNearField(); testRsfSpringSlider(); testRsfFaultCycle(); testRsfInjection(); testRsfEmbedded(); testKernelFitAspo(); testPhaseFieldSneddonGriffith(); testFaultMechStress(); testFaultMechRun(); testC25HexPatch(); testC25CornerPointMesh(); testC25RoundTrip(); testC25BlackOilCoupled(); testC25Comp(); testC25Thermal(); testC26Kernel(); testC26Penny(); testC26Shear(); testC26RoundTrip(); testC26FaultRun(); testC26HfWidth(); testC26IterCoupling(); testC26PhaseField(); testC27TimestepTargets(); testC27Wecon(); testC27VanGenuchten(); testC27HmStatistics(); testC27CoupledStudy(); testC27Examples(); testC27Tracer(); testC28Fosm(); testC28Tracers(); testC28Geostat(); testC28GpProxy(); testC28ProdAnalysis(); testC28FaultSlipEngines(); testC28DistAssembly(); testC28LevelIlu(); testC28PhaseFieldInRun(); testC29Salt(); testC29Wellbore(); testC29Network(); testC29DeviceCprAmg(); testC29PfTwoWay(); testC30Network(); testC30HfField(); testC30PhaseFieldOutput(); testC30PostGeomech(); testC30Colormaps(); testC31Expressions(); testC31Deck(); testC32MswUnit(); testC32MswDeck(); testC32MswAicdItems(); testC32GeoFormats(); testC32GeoWorkflow(); testC32ActionConnections(); testC33Report(); testC33Api(); testC33Regress(); testC33ScheduleItems(); testC34Items(); testC35Chem(); testC35Solvent(); testC35GasLift(); testC35Facilities(); testC36Seismic(); testC36Geo(); testC36Msw(); testC36Network(); testC37Msw(); testC37Py(); testC37Opm(); testC38ProdData(); testC38FloatIlu(); testC38Nonnc(); testC38Udt(); testC38TranEdit(); testC40GeoAnalysis(); testC40StudyGeoParams(); testC41Workflow(); testC41Facies(); testC41FaultNetwork(); testC42Import(); testC42Wells(); testC42Seismic2D(); testC42Calc(); testC42Cecon(); testC43Mistie(); testC43Polygons(); testC43LayerCake(); testC43Seismic(); testC43GeoCapi(); testC43Wpwe(); testC44FaultCut(); testC44ShutInStop(); testC44WelopenLayer(); testC44WpweReference(); testC44RestartEcon(); testC44HysteresisWetting(); testC44InjectorRecirculation(); testC44WpimultFbhpdef(); testC45LeanSolver(); testC45MswCrossflow(); testC46GeoProject(); testC46CompletionNumbers(); testC46SchedMultflt(); testC47GasPadding(); testC47Co2storeImport(); testC47FaultCsv(); testC47ActionWellsJly(); testC48Leakc(); testC48RefDepthShut(); testC48GroupUda(); testC48CellSnapshots(); testC49MswOptsRoundTrip(); testC49ConvergedHeads(); testC49FaultRamp(); testC49StoppedWellBhp(); testC49CompWellSummary(); testSeismicity(); testSeismicityRobust(); testFaultStability(); testPinch(); testDeckRepair(); testObsFromSummary(); testStoneAndDrsdt(); testNornePhysics(); testSPE10(); testEndPointScaling(); testPassThrough(); testAquifers(); testExamples(); testStreamlines(); testStudy(); testAdjoint(); testAdjointField(); testComponents(); testStudyAdvanced(); testBL(); testCorrelations(); testDemo(verbose); }
    catch (const std::exception& ex) { std::printf("EXCEPTION: %s\n", ex.what()); return 2; }
    std::printf(failures ? "%d test(s) FAILED\n" : "All tests passed\n", failures);
    return failures ? 1 : 0;
}
