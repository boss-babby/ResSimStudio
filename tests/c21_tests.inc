// ===================================================================================== C21 tests (included by run_tests.cpp)
#include "ressim/linsolve.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static void setEnvVar(const char* k, const char* v) {
#ifdef _WIN32
    _putenv_s(k, v ? v : "");
#else
    if (v) setenv(k, v, 1); else unsetenv(k);
#endif
}

// H3: level-scheduled BILU0 must reproduce the serial factorisation / sweeps exactly; thread default = physical cores
static void testParallelPerf() {
    std::printf("== C21 parallel performance fixes (level-scheduled BILU0, thread default)\n");
    const BSR A = makeBlockLaplace(120, 100, 3, 11);
    const int N = A.size();
    std::vector<double> r(static_cast<size_t>(N)), z0(static_cast<size_t>(N)), z1(static_cast<size_t>(N));
    for (int e = 0; e < N; ++e) r[size_t(e)] = std::cos(0.017 * e) + 0.05 * (e % 7);
    setNumThreads(2);
    setEnvVar("RESSIM_LEVEL_ILU", nullptr);
    BILU0 serial; serial.factor(A); serial.apply(r.data(), z0.data());
    setEnvVar("RESSIM_LEVEL_ILU", "1");
    BILU0 level; level.factor(A); level.apply(r.data(), z1.data());
    level.factor(A); level.apply(r.data(), z1.data());          // second factor: same pattern, values only
    setEnvVar("RESSIM_LEVEL_ILU", nullptr);
    double dmax = 0; for (int e = 0; e < N; ++e) dmax = std::max(dmax, std::fabs(z0[size_t(e)] - z1[size_t(e)]));
    std::printf("    %d unknowns: max |z_level - z_serial| = %.3g\n", N, dmax);
    check(dmax == 0.0, "level-scheduled BILU0 is bitwise identical to the serial BILU0");
    const int pc = physicalCores();
    setNumThreads(1); setNumThreads(0);
#ifdef _OPENMP
    std::printf("    physical cores %d, threads after setNumThreads(0) %d\n", pc, omp_get_max_threads());
    {   // C33: an explicit OMP_NUM_THREADS wins over the automatic choice
        const char* env = std::getenv("OMP_NUM_THREADS");
        const int want = env && std::atoi(env) > 0 ? std::atoi(env) : pc;
        check(pc >= 1 && omp_get_max_threads() == want, "setNumThreads(0) resets to OMP_NUM_THREADS or one thread per physical core (after a 1-thread run)");
    }
#else
    check(pc >= 1, "physicalCores() >= 1");
#endif
}

// H2: sub-grid fracture tip (continuous wing length, implicit K_I = K_IC inside the step)
static void testHfSubgrid() {
    std::printf("== C21 sub-grid hydraulic-fracture tip (HFGROW ... SUBGRID)\n");
    auto deck = [](const std::string& grow) {
        std::string d = hfDeck(true, false, grow);
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), "WCONINJE\n P WATER OPEN RATE 2 1* 700 /");
        const size_t q = d.find(" 0.2 0 1 0\n");
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        return d;
    };
    const std::string sg = "-- RESSIM: HFGROW F1 1.0 340 0 20 1 SUBGRID\n", dq = "-- RESSIM: HFGROW F1 1.0 340 0 20\n";
    try {
        Case c = eclipseDeckFromString(deck(sg), "");
        check(!c.grid.hfracs.empty() && c.grid.hfracs[0].subgrid == 1, "HFGROW ... SUBGRID is parsed");
        const std::string out = eclipseDeckToString(c);
        check(out.find("SUBGRID") != std::string::npos, "the deck writer keeps the SUBGRID flag");
        Simulator s; s.log = [](const std::string&) {};
        s.setup(c); s.initialize();
        const double L0 = s.hfWingLength(0, 1);
        std::printf("    initial wing length %.4f m (input 20)\n", L0);
        check(std::fabs(L0 - 20.0) < 1e-9, "sub-grid tip starts exactly at the initial half-length (no cell rounding)");
        std::string mid; int stepMid = -1;
        while (s.runReportStep()) if (s.reportStep() == 25) { mid = s.saveState(); stepMid = s.reportStep(); }
        const double L1 = s.hfWingLength(0, 1), Lm = s.hfWingLength(0, -1);
        double worst = 0, worstEdge = 0; int ev = 0, edgeEv = 0; bool mono = true; double last = 0;
        for (const auto& e : s.hfGrowthEvents()) {
            if (e.side > 0) { if (e.length < last - 1e-9) mono = false; last = e.length; }
            if (e.length < 99.9) {
                const double dv = std::fabs(e.KI / 1e6 - 1.0);
                // K_I jumps where the tip enters a new 10 m host cell (boundaries at 5 + 10k m from the well): there the
                // root finder stops on its 0.05 m length bracket (C21 tolerance), so K_I is not a meaningful test
                const double ph = std::fmod(e.length - 5.0 + 1000.0, 10.0);
                const bool atEdge = ph < 0.06 || ph > 10.0 - 0.06;
                if (dv > 0.08) std::printf("      event t %.4f d side %+d L %.3f m K_I/K_IC %.4f%s\n", e.t, e.side, e.length, e.KI / 1e6, atEdge ? " (host-cell edge)" : "");
                if (atEdge) { ++edgeEv; worstEdge = std::max(worstEdge, dv); } else worst = std::max(worst, dv);
            }
            ++ev;
        }
        std::printf("    final wings %.2f / %.2f m, %d growth steps, max |K_I/K_IC - 1| at the accepted lengths %.3f\n", Lm, L1, ev, worst);
        check(L1 > 30 && Lm > 30, "the sub-grid wings grow under injection");
        check(mono, "the continuous tip never recedes");
        std::printf("    %d accepted lengths on a host-cell edge (bracket-converged, max |K_I/K_IC - 1| %.3f)\n", edgeEv, worstEdge);
        check(worst < 0.10 && worstEdge < 0.25, "K_I = K_IC (within 10 %) at every accepted sub-grid length off the host-cell edges (where K_I jumps)");
        Case cq = eclipseDeckFromString(deck(dq), "");
        Simulator q; q.log = [](const std::string&) {};
        q.setup(cq); q.initialize(); while (q.runReportStep()) {}
        const double Lq = q.hfWingLength(0, 1);
        std::printf("    piece-quantised scheme: final wing %.2f m\n", Lq);
        check(std::fabs(L1 - Lq) <= 15.0, "sub-grid and piece-quantised final lengths agree within one 10 m piece + 5 m");
        if (stepMid > 0) {
            Simulator r; r.log = [](const std::string&) {};
            r.setup(c); r.loadState(mid); while (r.runReportStep()) {}
            std::printf("    restart from step %d: final wing %.6f vs %.6f m\n", stepMid, r.hfWingLength(0, 1), L1);
            check(std::fabs(r.hfWingLength(0, 1) - L1) < 1e-6, "restart (state v3) continues the sub-grid tip exactly");
        }
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "sub-grid tip run"); }
}

// G1-c: Barton-Bandis aperture law, EFEM pressure profiles, staggered EDFM <-> EFEM coupling
#include "ressim/geocouple.h"
static void testGeoCouple() {
    std::printf("== C21 G1-c flow <-> contact-mechanics coupling\n");
    {
        geocouple::BartonBandis bb;   // aRes 2e-5, Vm 3e-4, Kn0 2e10
        const double a0 = geocouple::aperture(bb, 0, 0, 0), ah = geocouple::aperture(bb, bb.Kn0 * bb.Vm, 0, 0), ainf = geocouple::aperture(bb, 1e15, 0, 0);
        const double ao = geocouple::aperture(bb, 5e6, 1e-3, 0);
        geocouple::BartonBandis bd = bb; bd.dilationDeg = 5;
        const double ad = geocouple::aperture(bd, 0, 0, 1e-3), adc = geocouple::aperture(bd, 0, 0, 1.0);
        std::printf("    BB: a(0) %.4g  a(Kn0 Vm) %.4g  a(inf) %.4g  open %.4g  dil %.4g  dil cap %.4g (m)\n", a0, ah, ainf, ao, ad, adc);
        check(std::fabs(a0 - 3.2e-4) < 1e-12 && std::fabs(ah - (2e-5 + 1.5e-4)) < 1e-12 && std::fabs(ainf - 2e-5) < 1e-8,
              "Barton-Bandis closure: a(0) = a_res + V_m, half closure at s'_n = K_n0 V_m, a -> a_res");
        check(std::fabs(ao - (3.2e-4 + 1e-3)) < 1e-12, "open fracture: a = a_res + V_m + w_n");
        check(std::fabs(ad - (3.2e-4 + 1e-3 * std::tan(5 * M_PI / 180))) < 1e-12 && std::fabs(adc - (3.2e-4 + 5e-4)) < 1e-12, "shear dilation |w_t| tan(psi), capped");
    }
    {   // EFEM: closed horizontal crack, t_n = s_yy + p; a constant profile reproduces the uniform pressure exactly
        efem::Input in; in.Lx = 10; in.Ly = 10; in.nx = 40; in.ny = 40; in.E = 10e9; in.nu = 0.25; in.sxx = -40e6; in.syy = -34e6;
        efem::Fracture F; F.x0 = 3; F.y0 = 5.125; F.x1 = 7; F.y1 = 5.125; F.pressure = 30e6; F.friction = 0.6;
        in.fracs = {F};
        efem::Result r0, r1, r2;
        const bool ok0 = efem::solve(in, r0);
        in.fracs[0].profS = {-2, 2}; in.fracs[0].profP = {30e6, 30e6};
        const bool ok1 = efem::solve(in, r1);
        double tnm = 0, d = 0; int nc = 0;
        for (size_t i = 0; i < r0.cut.size(); ++i) {
            if (std::fabs(r0.cut[i].s) < 1.0) { tnm += r0.cut[i].tn; ++nc; }   // central part (away from the crack-tip elements)
            if (i < r1.cut.size()) d = std::max(d, std::fabs(r0.cut[i].tn - r1.cut[i].tn));
        }
        tnm /= std::max(1, nc);
        std::printf("    closed crack: mean t_n %.4f MPa (expected -4), |profile - uniform| %.3g Pa\n", tnm / 1e6, d);
        check(ok0 && ok1 && std::fabs(tnm + 4e6) < 0.02 * 4e6, "EFEM closed crack: t_n = s_yy + p (-4 MPa within 2 %)");
        check(ok1 && d < 1e-3, "a constant pressure profile reproduces the uniform-pressure solution");
        in.fracs[0].profP = {37e6, 31e6};   // 37 MPa at the (3,5) end, 31 MPa at (7,5): only the left part exceeds 34 MPa
        const bool ok2 = efem::solve(in, r2);
        double wl = 0, wr = 0;
        for (const auto& C : r2.cut) { if (C.s < -1.0) wl = std::max(wl, C.wn); if (C.s > 1.0) wr = std::max(wr, C.wn); }
        std::printf("    linear profile 37 -> 31 MPa: max opening left %.3g mm, right %.3g mm\n", wl * 1e3, wr * 1e3);
        check(ok2 && wl > 0 && wl > 10 * wr, "linear pressure profile opens the high-pressure end of the crack only");
    }
    // coupled EDFM injection (hfDeck fracture 200 m along x, p0 ~ 300 bar)
    auto deck = [](const std::string& rate) {
        std::string d = hfDeck(true, false, "");
        const size_t hw = d.find(" 0.003 0 0.5 0 0");   // 0.1 mm HFRAC width: finite conductivity (F_CD ~ 10), the aperture matters
        if (hw != std::string::npos) d.replace(hw, 6, " 0.0001");
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), "WCONINJE\n P WATER OPEN RATE " + rate + " 1* 900 /");
        const size_t q = d.find(" 0.2 0 1 0\n");
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        return d;
    };
    auto run = [&](const std::string& rate, const geocouple::Input* gc, geocouple::Output* go, double& bhp) {
        try {
            Case c = eclipseDeckFromString(deck(rate), "");
            Simulator s; s.setup(c); s.initialize();
            if (gc) { if (!geocouple::runCoupled(s, *gc, *go)) { std::printf("    coupling error: %s\n", go->error.c_str()); return false; } }
            else while (s.runReportStep()) {}
            bhp = s.summary().back().wells[0].bhp;
            return true;
        } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); return false; }
    };
    double b0 = 0, bc = 0, bo = 0;
    geocouple::Input g; g.E = 20e9; g.nu = 0.25; g.sxx = -45e6; g.syy = -40e6; g.nx = 60; g.ny = 60;
    g.bb.aRes = 1e-5; g.bb.Vm = 1e-4;   // closed Barton-Bandis aperture ~ 3.5e-5 m < the 1e-4 m HFRAC width
    geocouple::Output oc, oo;
    const bool r0 = run("2", nullptr, nullptr, b0);
    const bool rc = run("2", &g, &oc, bc);
    if (rc && std::getenv("RESSIM_GCDBG")) for (size_t k = 0; k < oc.s.size(); ++k) std::printf("      piece %zu s %.1f p %.2f MPa tn %.3f MPa a %.3g m mult %.3g\n", k, oc.s[k], oc.steps.back().p[k] / 1e6, oc.steps.back().tn[k] / 1e6, oc.steps.back().a[k], oc.steps.back().mult[k]);
    std::printf("    rate 2: uncoupled BHP %.2f bar, coupled (closed, s'_n ~ 6 MPa) %.2f bar, conductance x %.3g .. %.3g, open pieces %d\n",
                b0, bc, rc ? oc.steps.back().minMult : 0, rc ? oc.steps.back().maxMult : 0, rc ? oc.steps.back().openPieces : -1);
    check(r0 && rc && oc.steps.size() > 10 && oc.steps.back().openPieces == 0 && oc.steps.back().maxMult < 1 && bc > b0 + 0.1,
          "closed coupled fracture: Barton-Bandis aperture < HFRAC width, BHP above the uncoupled one");
    geocouple::Input go; go.E = 20e9; go.nu = 0.25; go.sxx = -40e6; go.syy = -31.5e6; go.nx = 60; go.ny = 60;   // default BB: closed a >= 1e-4 m
    const bool ro = run("5", &go, &oo, bo);
    double bu = 0; const bool ru = run("5", nullptr, nullptr, bu);
    double aw = 0, at = 0; int nw = 0, nt = 0;
    if (ro) {
        const auto& st = oo.steps.back();
        double smax = 0; for (double s : oo.s) smax = std::max(smax, std::fabs(s));
        for (size_t k = 0; k < oo.s.size(); ++k) { if (std::fabs(oo.s[k]) < 15) { aw += st.a[k]; ++nw; } if (std::fabs(oo.s[k]) > smax - 15) { at += st.a[k]; ++nt; } }
        aw /= std::max(1, nw); at /= std::max(1, nt);
        std::printf("    rate 5, s_yy -31.5 MPa: uncoupled BHP %.2f, coupled %.2f bar; open pieces %d / %zu, aperture at well %.3g mm, at tips %.3g mm, max opening %.3g mm\n",
                    bu, bo, st.openPieces, oo.s.size(), aw * 1e3, at * 1e3, st.maxOpening * 1e3);
    }
    check(ro && ru && oo.steps.back().openPieces > 0 && aw > at && bo < bu, "pressurised coupled fracture opens near the injector (larger aperture there) and lowers the BHP");
}
