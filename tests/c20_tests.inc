// ===================================================================================== C20 tests (included by run_tests.cpp)
#include "ressim/seismicity.h"

// rate-and-state seismicity: exact integral forecast vs RK4 of the Heimisson (2019) ODE, gradients vs finite differences,
// convexity, the analytic 3D-diffusion design problem of research/CYCLE1 (window flip, global optimum) and the
// engine-driven designer (Green's functions + nonlinear verification) on a water-injection case
static void testSeismicity() {
    using namespace ressim::seis;
    std::printf("Rate-and-state seismicity + convex injection designer (C20)\n");
    // ---- 1. forecast vs ODE
    {
        Patch p; p.aSigma = 0.05; p.taYears = 0.5; p.rate = 10;
        std::vector<double> t, s;
        for (int i = 0; i <= 2000; ++i) {
            const double ti = 20.0 * i / 2000;
            t.push_back(ti);
            s.push_back(ti < 5 ? 0.2 * (1 - std::exp(-ti / 2)) : 0.2 * (1 - std::exp(-2.5)) * std::exp(-(ti - 5) / 3));
        }
        const Forecast f = forecast(t, s, p, 5.0);
        const double ta = p.taYears * 365.25, a = p.aSigma;
        auto dS = [&](double x) {
            const size_t k = std::min<size_t>(size_t(x / 0.01), 1999);
            const double w = (x - t[k]) / 0.01;
            return s[k] + w * (s[k + 1] - s[k]);
        };
        auto rhs = [&](double x, double u) { return std::exp((x - u) / ta + dS(x) / a); };
        double u = 0; const double h = 1e-3;
        for (int i = 0; i < 20000; ++i) {
            const double x = i * h;
            const double k1 = rhs(x, u), k2 = rhs(x + h / 2, u + h / 2 * k1), k3 = rhs(x + h / 2, u + h / 2 * k2), k4 = rhs(x + h, u + h * k3);
            u += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
        const double Node = p.rate / 365.25 * u, excessOde = Node - p.rate / 365.25 * 20;
        const double rel = std::fabs(f.excess - excessOde) / std::fabs(excessOde);
        std::printf("    excess events: integral %.6f  ODE %.6f  rel %.2e; peak R %.3g /yr at %.2f d\n", f.excess, excessOde, rel, f.peakR, f.tPeak);
        check(rel < 1e-4, "Heimisson-Segall exact integral matches the rate-and-state ODE (rel < 1e-4)");
        check(f.excessAtShutIn > 0 && f.excessAtShutIn < f.excess, "excess at shut-in is positive and below the full-window excess");
    }
    // ---- 2. analytic 3D diffusion kernel (research/exp/h1_seis_opt.py baseline)
    const double DAY = 86400.0, kperm = 1e-14, mu = 1e-3, cdiff = 0.1, fr = 0.6, rr = 300.0;
    auto pstep = [&](double tsec) { return tsec <= 0 ? 0.0 : mu / (4 * M_PI * kperm * rr) * std::erfc(rr / (2 * std::sqrt(cdiff * tsec))); };
    auto build = [&](double thFac, int nc) {
        LinearModel m; Patch p; p.aSigma = 0.05; p.taYears = 100; p.rate = 365.25; p.friction = fr;
        m.patches = {p}; m.nv = nc;
        const double Ti = 10.0, Th = thFac * Ti;
        for (int i = 0; i <= 800; ++i) m.t.push_back(Ti * i / 800);
        if (Th > Ti) for (int i = 1; i <= 800; ++i) { const double x = double(i) / 800; m.t.push_back(Ti + (Th - Ti) * x * x); }
        m.dS0 = {std::vector<double>(m.t.size(), 0.0)};
        m.G = {std::vector<double>(m.t.size() * size_t(nc))};
        for (size_t i = 0; i < m.t.size(); ++i)
            for (int j = 0; j < nc; ++j) {
                const double t0 = Ti * j / nc, t1 = Ti * (j + 1) / nc;
                // per unit rate in m3/day: m3/s = q/86400 ; stress in MPa
                m.G[0][i * size_t(nc) + size_t(j)] = fr * (pstep((m.t[i] - t0) * DAY) - pstep((m.t[i] - t1) * DAY)) / DAY / 1e6;
            }
        return m;
    };
    const int nc = 50; const double V = 8640.0, q0 = V / 10.0;           // m3, m3/day
    auto sched = [&](int kind) {
        std::vector<double> q(static_cast<size_t>(nc));
        for (int j = 0; j < nc; ++j) { const double x = (j + 0.5) / nc; q[size_t(j)] = kind == 0 ? 1.0 : std::floor(x * 4) + 1.0; }
        double s = 0; for (double v : q) s += v * 10.0 / nc;
        for (double& v : q) v *= V / s;
        return q;
    };
    Constraints c; c.lb.assign(size_t(nc), 0.0); c.ub.assign(size_t(nc), 3 * q0); c.group.assign(size_t(nc), 0); c.weight.assign(size_t(nc), 10.0 / nc); c.target = {V};
    {
        const LinearModel m1 = build(1.0, nc), m5 = build(5.0, nc);
        const std::vector<double> qc = sched(0), qs = sched(1);
        const double c1 = modelEvents(m1, qc), s1 = modelEvents(m1, qs), c5 = modelEvents(m5, qc), s5 = modelEvents(m5, qs);
        std::printf("    step-up vs constant: co-injection window %+.1f%%, with 4x Ti tail %+.1f%%\n", 100 * (s1 / c1 - 1), 100 * (s5 / c5 - 1));
        check(s1 < c1 && s5 > c5, "counting-window artefact: step-up wins during injection, loses once the tail is counted");
        // gradients
        std::vector<double> g; modelEvents(m5, qs, &g);
        std::vector<double> gp; const double beta = 40; modelPeak(m5, qs, beta, &gp);
        double worst = 0, worstP = 0;
        for (int j : {3, 17, 33, 48}) {
            std::vector<double> a = qs, b = qs; const double h = 1e-4 * q0; a[size_t(j)] += h; b[size_t(j)] -= h;
            const double fd = (modelEvents(m5, a) - modelEvents(m5, b)) / (2 * h);
            const double fdp = (modelPeak(m5, a, beta) - modelPeak(m5, b, beta)) / (2 * h);
            worst = std::max(worst, std::fabs(fd - g[size_t(j)]) / std::max(1e-30, std::fabs(fd)));
            worstP = std::max(worstP, std::fabs(fdp - gp[size_t(j)]) / std::max(1e-12, std::fabs(fdp) + 1e-9));
        }
        std::printf("    gradient vs FD: events rel %.2e, smooth peak rel %.2e\n", worst, worstP);
        check(worst < 1e-5 && worstP < 1e-4, "analytic gradients of the event count and the smooth peak match finite differences");
        // convexity on random pairs
        unsigned st = 12345; auto rnd = [&]() { st = st * 1664525u + 1013904223u; return (st >> 8) / double(1u << 24); };
        int viol = 0;
        for (int k = 0; k < 200; ++k) {
            std::vector<double> a(static_cast<size_t>(nc)), b(static_cast<size_t>(nc)), mid(static_cast<size_t>(nc));
            const double lam = rnd();
            for (int j = 0; j < nc; ++j) { a[size_t(j)] = 3 * q0 * std::pow(rnd(), 4); b[size_t(j)] = 3 * q0 * std::pow(rnd(), 4); mid[size_t(j)] = lam * a[size_t(j)] + (1 - lam) * b[size_t(j)]; }
            if (modelEvents(m5, mid) > lam * modelEvents(m5, a) + (1 - lam) * modelEvents(m5, b) + 1e-12 * std::fabs(modelEvents(m5, a))) ++viol;
        }
        check(viol == 0, "expected event count is convex in the rate history (200 random midpoint tests)");
        // global optimum from two starts
        OptOptions o; o.maxIter = 6000; o.tol = 1e-11;
        const OptResult r1 = optimize(m5, c, qc, o), r2 = optimize(m5, c, qs, o);
        double vol = 0; for (double v : r1.x) vol += v * 10.0 / nc;
        const double gain = r1.events / c5 - 1;
        std::printf("    optimum: %+.2f%% vs constant (Python reference -4.7%%), starts agree to %.1e, KKT %.1e, %d its, volume err %.1e\n",
                    100 * gain, std::fabs(r1.events - r2.events) / std::fabs(r1.events), r1.kkt, r1.iters, std::fabs(vol - V) / V);
        check(gain < -0.035 && gain > -0.065, "convex optimum reduces expected events by 3.5-6.5% (research reference -4.7%)");
        check(std::fabs(r1.events - r2.events) < 2e-4 * std::fabs(r1.events), "two different starts reach the same global optimum");
        check(std::fabs(vol - V) < 1e-8 * V, "injected volume constraint holds exactly");
        // peak objective (LP surrogate as a smooth maximum) reduces the maximal stress
        OptOptions op; op.obj = Objective::Peak; op.maxIter = 6000;
        const OptResult rp = optimize(m5, c, qc, op);
        double mc = 0, mp = 0;
        for (double v : modelStress(m5, qc, 0)) mc = std::max(mc, v);
        for (double v : modelStress(m5, rp.x, 0)) mp = std::max(mp, v);
        std::printf("    peak objective: max dS %.4f -> %.4f MPa (research LP: -20%%)\n", mc, mp);
        check(mp < 0.88 * mc, "peak-rate design lowers the maximal Coulomb stress by > 12%");
    }
    // ---- 3. engine-driven designer on a water-filled reservoir (single phase -> nearly linear)
    {
        const char* deckText = R"(RUNSPEC
DIMENS
 41 41 1 /
OIL
WATER
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 1 1 1 1 /
START
 1 JAN 2025 /
GRID
DX
 1681*20 /
DY
 1681*20 /
DZ
 1681*10 /
TOPS
 1681*2000 /
PORO
 1681*0.1 /
PERMX
 1681*10 /
PERMY
 1681*10 /
PERMZ
 1681*1 /
PROPS
PVCDO
 200 1.0 1e-4 1.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
DENSITY
 800 1000 1 /
ROCK
 200 4e-5 /
SWOF
 0.2 0 1 0
 0.3 0.02 0.7 0
 1.0 1 0 0 /
SOLUTION
EQUIL
 2005 200 1000 0 /
SCHEDULE
WELSPECS
 I G 21 21 2005 WATER /
/
COMPDAT
 I 21 21 1 1 OPEN 1* 1* 0.2 /
/
WCONINJE
 I WATER OPEN RATE 100 1* 900 /
/
TSTEP
 10*1 /
END
)";
        ressim::Case cs = ressim::eclipseDeckFromString(deckText, "");
        DesignSetup s;
        s.wells = {"I"};
        Patch p; p.name = "F1"; p.i = 26; p.j = 20; p.k = 0; p.aSigma = 0.1; p.taYears = 50; p.rate = 2;
        s.patches = {p};
        s.nCtrl = 6; s.injDays = 6; s.postDays = 12; s.subSteps = 2; s.postSteps = 8; s.totalVolume = 600;
        const DesignOutput d = designInjection(cs, s);
        for (const auto& l : d.logLines) std::printf("      %s\n", l.c_str());
        check(d.ok, "engine-driven seismicity design runs (base + step-response + 2 verification runs)");
        if (d.ok) {
            double vol = 0; for (double v : d.qOpt[0]) vol += v * s.injDays / s.nCtrl;
            check(std::fabs(vol - 600) < 1e-6 * 600, "designed schedule injects the requested volume");
            check(d.eventsOpt <= d.eventsConst * (1 + 1e-9), "optimal schedule predicted no worse than constant");
            check(d.simEventsOpt < d.simEventsConst, "nonlinear forward run confirms fewer events with the optimal schedule");
            std::printf("    first-pass model error %.1f%%, after %d correction pass(es) %.2f%%\n", 100 * d.linearityError, d.corrections, 100 * d.finalModelError);
            check(d.finalModelError < 0.03, "after model correction (SCP) the design model matches the nonlinear simulator within 3%");
            check(d.simEventsOpt <= d.corrSim.front() * (1 + 1e-12), "returned schedule is the best verified one");
            const std::string deck = designScheduleDeck(cs, s, d.qOpt);
            check(deck.find("WCONINJE") != std::string::npos && deck.find("TSTEP") != std::string::npos, "design exports as WCONINJE / TSTEP schedule text");
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// C20: discrete adjoint with in-run hydraulic-fracture growth (frozen topology per recorded step)
static void testHfGrowthAdjoint() {
    std::printf("== C20 HFGROW adjoint: gradients through a growing fracture (open set frozen per recorded step) vs full-run finite differences\n");
    std::string d = hfDeck(true, false, "-- RESSIM: HFGROW F1 1.0 320 0 20\n");
    {
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), "WCONINJE\n P WATER OPEN RATE 2 1* 700 /");
        const size_t q = d.find(" 0.2 0 1 0\n");
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        const size_t ts = d.find("TSTEP\n 10*0.1 10*0.5 10*1 10*2 10*5 /");   // shorter schedule (test time)
        if (ts != std::string::npos) d.replace(ts, 37, "TSTEP\n 10*0.1 10*0.5 10*1 5*2 /");
        // fixed time-step sequence (0.1 d, then 0.5 d): the finite-difference reference runs take the same steps
        const size_t sc = d.find("SCHEDULE\n");
        if (sc != std::string::npos) d.insert(sc + 9, "TUNING\n 0.1 0.5 0.1 /\n/\n/\n");
    }
    Case base;
    try { base = eclipseDeckFromString(d, ""); } catch (const std::exception& e) { std::printf("    deck: %s\n", e.what()); check(false, "HFGROW adjoint: deck"); return; }
    std::vector<StudyParam> G(2);
    G[0].kind = ParamKind::PermXY; G[1].kind = ParamKind::Poro;
    for (auto& p : G) { p.active = true; p.log = false; p.base = 1; p.low = 0.5; p.high = 1.5; }
    struct Run { std::vector<SummaryPoint> s; int ev = 0, steps = 0; double wl = 0; };
    auto run = [&](const std::vector<double>& v, Simulator* keep) {
        Run r;
        Case c = applyParams(base, G, v);
        Simulator local; Simulator& s = keep ? *keep : local;
        s.log = [](const std::string&) {};
        if (keep) s.setAdjointRecording(true);
        s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        r.s = s.summary(); r.ev = int(s.hfGrowthEvents().size()); r.steps = s.stats().timesteps; r.wl = s.hfWingLength(0, 1);
        return r;
    };
    Simulator S;
    Run b;
    try { b = run({1.0, 1.0}, &S); } catch (const std::exception& e) { std::printf("    base run: %s\n", e.what()); check(false, "HFGROW adjoint: base run"); return; }
    // observations: injector BHP and field pressure at several times (data = 1 % off the base response)
    std::vector<Simulator::AdjObs> obs;
    for (double t : {1.0, 3.0, 6.0, 12.0, 24.0}) {
        for (int q = 0; q < 2; ++q) {
            Simulator::AdjObs o; o.qty = q == 0 ? ObsQty::WBHP : ObsQty::FPR; o.well = q == 0 ? 0 : -1; o.time = t;
            double y = 0; for (const auto& p : b.s) if (std::fabs(p.time - t) < 1e-6) y = q == 0 ? (p.wells.empty() ? 0 : p.wells[0].bhp) : p.fpr;
            o.value = 1.01 * y; o.sigma = 1.0; obs.push_back(o);
        }
    }
    auto misfit = [&](const Run& r) {
        double F = 0;
        for (const auto& o : obs) {
            double y = 0; for (const auto& p : r.s) if (std::fabs(p.time - o.time) < 1e-6) y = o.qty == ObsQty::WBHP ? (p.wells.empty() ? 0 : p.wells[0].bhp) : p.fpr;
            F += (y - o.value) * (y - o.value) / (o.sigma * o.sigma);
        }
        return F;
    };
    const double h = 1e-4;
    std::vector<Case> cp(2), cm(2);
    std::vector<Simulator::AdjParam> ap(2);
    for (int j = 0; j < 2; ++j) {
        std::vector<double> vp = {1.0, 1.0}, vm = {1.0, 1.0};
        vp[size_t(j)] += h; vm[size_t(j)] -= h;
        cp[size_t(j)] = applyParams(base, G, vp); cm[size_t(j)] = applyParams(base, G, vm);
        ap[size_t(j)].plus = &cp[size_t(j)]; ap[size_t(j)].minus = &cm[size_t(j)]; ap[size_t(j)].h = h;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const Simulator::AdjResult g = S.adjointGradient(obs, ap, false);
    const double ta = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!g.ok) std::printf("    adjoint error: %s\n", g.error.c_str());
    check(g.ok, "HFGROW adjoint: gradient available with in-run fracture growth");
    if (!g.ok) return;
    // full-run central differences. The reference must stay on the recorded growth branch (same opening steps / step
    // sequence): try decreasing perturbations and use the largest one that keeps the branch
    double fd[2] = {0, 0}, e0 = 1, e1 = 1, hUsed = 0; bool same = false;
    for (double hf : {2e-3, 2e-4, 2e-5}) {
        bool sm = true; double f[2];
        for (int j = 0; j < 2; ++j) {
            std::vector<double> vp = {1.0, 1.0}, vm = {1.0, 1.0};
            vp[size_t(j)] += hf; vm[size_t(j)] -= hf;
            const Run rp = run(vp, nullptr), rm = run(vm, nullptr);
            sm = sm && rp.ev == b.ev && rm.ev == b.ev && rp.steps == b.steps && rm.steps == b.steps;
            f[j] = (misfit(rp) - misfit(rm)) / (2 * hf);
        }
        const double a0 = std::fabs(g.grad[0] - f[0]) / std::max(std::fabs(f[0]), 1e-12), a1 = std::fabs(g.grad[1] - f[1]) / std::max(std::fabs(f[1]), 1e-12);
        std::printf("    FD h %.0e: dF/dk %.6g, dF/dphi %.6g (rel %.2e / %.2e), same branch: %s\n", hf, f[0], f[1], a0, a1, sm ? "yes" : "no");
        if (sm) { fd[0] = f[0]; fd[1] = f[1]; e0 = a0; e1 = a1; hUsed = hf; same = true; break; }
    }
    std::printf("    %d growth events (wing %.1f m), %d steps; F %.5g (forward %.5g); dF/dk: adjoint %.6g FD %.6g (rel %.2e); dF/dphi: adjoint %.6g FD %.6g (rel %.2e);"
                " FD step %.0e; adjoint %.2f s\n",
                b.ev, b.wl, b.steps, g.F, misfit(b), g.grad[0], fd[0], e0, g.grad[1], fd[1], e1, hUsed, ta);
    check(b.ev > 0, "HFGROW adjoint: the fracture grows during the recorded run");
    check(std::fabs(g.F - misfit(b)) <= 1e-6 * std::max(1.0, misfit(b)), "HFGROW adjoint: replayed misfit = forward misfit");
    check(same && e0 < 1e-2 && e1 < 1e-2, "HFGROW adjoint: gradient matches full-run finite differences (1 %) on the same growth branch");
}

// ---------------------------------------------------------------------------------------------------------------------
// C20: P3D hydraulic fracture (equilibrium height growth through symmetric stress barriers)
static void testP3D() {
    std::printf("== C20 P3D hydraulic fracture: equilibrium height (K_I = K_IC at both vertical tips) + PKN lateral flow\n");
    // section relations: inverse consistency and the zero-toughness containment curve H = h / cos(pi p_net / (2 dsigma))
    {
        const double Ep = 32e9, h = 30, ds = 3e6;
        double eInv = 0, eCurve = 0;
        for (double H : {30.5, 35.0, 45.0, 80.0, 200.0}) {
            const double A = hf::p3dArea(H, h, ds, 1e6, Ep);
            const hf::P3DSection S = hf::p3dSection(A, h, ds, 1e6, Ep, 1e4);
            eInv = std::max(eInv, std::fabs(S.H - H) / H);
            const hf::P3DSection S0 = hf::p3dSection(hf::p3dArea(H, h, ds, 0.0, Ep), h, ds, 0.0, Ep, 1e4);
            eCurve = std::max(eCurve, std::fabs(h / std::cos(M_PI * S0.pnet / (2 * ds)) - H) / H);
        }
        std::printf("    section: A(H) inverse rel %.1e; K_IC = 0 containment curve rel %.1e\n", eInv, eCurve);
        check(eInv < 1e-8 && eCurve < 1e-8, "P3D: equilibrium-height section relations (inverse, containment curve)");
    }
    hf::Input in;
    in.E = 30e9; in.nu = 0.25; in.height = 30; in.mu = 0.05; in.rate = 8.0 / 60; in.cl = 2e-5; in.kic = 1e6; in.sigmaMin = 40e6;
    in.tEnd = 5400; in.cells = 0; in.outputs = 10;
    hf::Result pkn, big, mid, low;
    hf::Input a = in; a.model = 0;
    const bool okP = hf::simulateAdaptive(a, pkn);
    a.model = 2; a.stressContrast = 1e9; const bool okB = hf::simulateAdaptive(a, big);
    a.stressContrast = 10e6; const bool okM = hf::simulateAdaptive(a, mid);
    a.stressContrast = 6e6; const bool okL = hf::simulateAdaptive(a, low);
    if (!okP || !okB || !okM || !okL) std::printf("    errors: %s | %s | %s | %s\n", pkn.error.c_str(), big.error.c_str(), mid.error.c_str(), low.error.c_str());
    check(okP && okB && okM && okL, "P3D: simulations complete");
    if (!(okP && okB && okM && okL)) return;
    const double eL = std::fabs(big.last().L / pkn.last().L - 1), eW = std::fabs(big.last().w0 / pkn.last().w0 - 1), eE = std::fabs(big.last().eff - pkn.last().eff);
    std::printf("    PKN: L %.1f m, w0 %.2f mm, eff %.3f, p_net %.2f MPa | P3D strong barriers: L %.1f, w0 %.2f, eff %.3f, H_w %.2f m (rel L %.1e, w0 %.1e)\n",
                pkn.last().L, pkn.last().w0 * 1e3, pkn.last().eff, pkn.last().pnet0 / 1e6, big.last().L, big.last().w0 * 1e3, big.last().eff, big.heightWell, eL, eW);
    check(eL < 0.01 && eW < 0.01 && eE < 0.005 && std::fabs(big.heightWell - 30) < 0.01, "P3D: strong barriers recover PKN (length, width, efficiency within 1 %)");
    auto monotone = [](const hf::Snapshot& S) { for (size_t i = 1; i < S.H.size(); ++i) if (S.H[i] > S.H[i - 1] * (1 + 1e-9)) return false; return true; };
    std::printf("    dsigma 10 MPa: L %.1f m, H_w %.1f m, mean H %.1f, p_net %.2f MPa, eff %.3f | dsigma 6 MPa: L %.1f, H_w %.1f, p_net %.2f, runaway %d\n",
                mid.last().L, mid.heightWell, mid.height, mid.last().pnet0 / 1e6, mid.last().eff, low.last().L, low.heightWell, low.last().pnet0 / 1e6, int(low.heightRunaway));
    check(mid.heightWell > 30.5 && low.heightWell > mid.heightWell && mid.last().L < pkn.last().L && low.last().L < mid.last().L,
          "P3D: weaker barriers -> taller, shorter fractures");
    check(monotone(mid.last()) && monotone(low.last()), "P3D: height decreases from the wellbore to the tip");
    // mass balance without leak-off
    hf::Input m = in; m.model = 2; m.stressContrast = 6e6; m.cl = 0; m.cells = 0;
    hf::Result rm; const bool okm = hf::simulateAdaptive(m, rm);
    const double eV = okm ? std::fabs(rm.volFracture - rm.volInjected) / rm.volInjected : 1;
    std::printf("    no leak-off: fracture volume %.4g m3 vs injected %.4g m3 (rel %.1e)\n", rm.volFracture, rm.volInjected, eV);
    check(okm && eV < 1e-6, "P3D: volume balance without leak-off");
}

// ---------------------------------------------------------------------------------------------------------------------
// C20 geomechanics G1: embedded finite elements with frictional contact (2-D plane strain)
#include "ressim/efem.h"
static void testEfem() {
    std::printf("== C20 G1 EFEM: embedded fractures with Coulomb contact vs Sneddon (1946) and the inclined frictional crack\n");
    const double E = 10e9, nu = 0.25, Ep = E / (1 - nu * nu);
    // (1) Sneddon: pressurised crack, w_n(s) = 4 p sqrt(a^2 - s^2) / E', no far-field stress
    for (double ang : {0.0, 30.0}) {
        efem::Input in; in.Lx = 10; in.Ly = 10; in.nx = 100; in.ny = 100; in.E = E; in.nu = nu;
        const double h = in.Lx / in.nx, th = ang * M_PI / 180, xl = 5 - 12 * h, xr = 5 + 12 * h;   // end points on vertical grid lines
        const double yc = 5 + 0.37 * h;
        efem::Fracture F; F.x0 = xl; F.x1 = xr; F.y0 = yc - (xr - xl) / 2 * std::tan(th); F.y1 = yc + (xr - xl) / 2 * std::tan(th); F.pressure = 1e6;
        in.fracs.push_back(F);
        efem::Result r;
        const bool ok = efem::solve(in, r);
        if (!ok) { std::printf("    error: %s\n", r.error.c_str()); check(false, "EFEM Sneddon: solve"); continue; }
        const double a = std::hypot(F.x1 - F.x0, F.y1 - F.y0) / 2, p = F.pressure;
        double V = 0, wmax = 0, wtmax = 0, eProf = 0;
        for (const auto& c : r.cut) {
            V += c.wn * c.len; wmax = std::max(wmax, c.wn); wtmax = std::max(wtmax, std::fabs(c.wt));
            // element mean of the analytic profile over the element's part of the crack
            const double s0 = std::max(-a, c.s - c.len / 2), s1 = std::min(a, c.s + c.len / 2);
            auto I = [&](double s) { return 0.5 * (s * std::sqrt(std::max(0.0, a * a - s * s)) + a * a * std::asin(std::clamp(s / a, -1.0, 1.0))); };
            const double wa = 4 * p / Ep * (I(s1) - I(s0)) / std::max(s1 - s0, 1e-12);
            if (std::fabs(c.s) < 0.7 * a) eProf = std::max(eProf, std::fabs(c.wn - wa) / (4 * p * a / Ep));
        }
        const double Va = 2 * M_PI * p * a * a / Ep, w0 = 4 * p * a / Ep;
        std::printf("    Sneddon %2.0f deg: %zu cut elements, volume %.4e vs %.4e m2 (rel %.3f), max opening %.4e vs %.4e (rel %.3f), interior profile error %.3f of w0, |w_t| max %.1e, Newton %d, %.2f s\n",
                    ang, r.cut.size(), V, Va, V / Va - 1, wmax, w0, wmax / w0 - 1, eProf, wtmax, r.newton, r.seconds);
        check(std::fabs(V / Va - 1) < 0.06 && std::fabs(wmax / w0 - 1) < 0.06 && eProf < 0.06 && wtmax < 0.05 * w0,
              ang == 0 ? "EFEM: Sneddon pressurised crack (volume, opening profile within 6 %)" : "EFEM: Sneddon crack inclined to the mesh (30 deg) within 6 %");
    }
    // (2) inclined crack under uniaxial compression with Coulomb friction: slip b(s) = 4 (|tau| - mu |sigma_n|) sqrt(a^2 - s^2) / E'
    for (double mu : {0.3, 0.8}) {
        efem::Input in; in.Lx = 10; in.Ly = 10; in.nx = 100; in.ny = 100; in.E = E; in.nu = nu; in.syy = -10e6;
        const double h = in.Lx / in.nx, th = 30 * M_PI / 180, xl = 5 - 12 * h, xr = 5 + 12 * h, yc = 5 + 0.37 * h;
        efem::Fracture F; F.x0 = xl; F.x1 = xr; F.y0 = yc - (xr - xl) / 2 * std::tan(th); F.y1 = yc + (xr - xl) / 2 * std::tan(th); F.friction = mu;
        in.fracs.push_back(F);
        efem::Result r;
        const bool ok = efem::solve(in, r);
        if (!ok) { std::printf("    error: %s\n", r.error.c_str()); check(false, "EFEM inclined crack: solve"); continue; }
        const double a = std::hypot(F.x1 - F.x0, F.y1 - F.y0) / 2, sig = 10e6;
        const double tau = sig * std::sin(th) * std::cos(th), sn = sig * std::cos(th) * std::cos(th), drive = std::max(0.0, tau - mu * sn);
        double S = 0, smax = 0, pen = 0, eProf = 0; int nslip = 0;
        for (const auto& c : r.cut) {
            S += std::fabs(c.wt) * c.len; smax = std::max(smax, std::fabs(c.wt)); pen = std::min(pen, c.wn); nslip += c.state == 2;
            const double s0 = std::max(-a, c.s - c.len / 2), s1 = std::min(a, c.s + c.len / 2);
            auto I = [&](double s) { return 0.5 * (s * std::sqrt(std::max(0.0, a * a - s * s)) + a * a * std::asin(std::clamp(s / a, -1.0, 1.0))); };
            const double ba = 4 * drive / Ep * (I(s1) - I(s0)) / std::max(s1 - s0, 1e-12);
            if (drive > 0 && std::fabs(c.s) < 0.7 * a) eProf = std::max(eProf, std::fabs(std::fabs(c.wt) - ba) / (4 * drive * a / Ep));
        }
        const double Sa = 2 * M_PI * drive * a * a / Ep, b0 = 4 * drive * a / Ep;
        if (drive > 0) {
            std::printf("    inclined 30 deg, mu %.1f: slip integral %.4e vs %.4e (rel %.3f), max slip %.4e vs %.4e (rel %.3f), interior profile error %.3f, %d / %zu slipping, max penetration %.1e m, Newton %d, %.2f s\n",
                        mu, S, Sa, S / Sa - 1, smax, b0, smax / b0 - 1, eProf, nslip, r.cut.size(), pen, r.newton, r.seconds);
            check(std::fabs(S / Sa - 1) < 0.08 && std::fabs(smax / b0 - 1) < 0.08 && eProf < 0.08 && nslip >= int(r.cut.size()) * 2 / 3,
                  "EFEM: frictional slip of an inclined crack in compression (Coulomb, within 8 %)");
        } else {
            const double bRef = 4 * tau * a / Ep;   // frictionless slip scale
            std::printf("    inclined 30 deg, mu %.1f (tan 30 < mu: locked): max |slip| %.2e m (frictionless %.2e), %d slipping, Newton %d\n", mu, smax, bRef, nslip, r.newton);
            check(nslip == 0 && smax < 1e-3 * bRef, "EFEM: crack stays locked when |tau| < mu sigma_n");
        }
    }
    // F9: the classic element-constant EFEM (jumpModel 1) on the same inclined problems - report only
    {
        efem::Input in; in.Lx = 10; in.Ly = 10; in.nx = 100; in.ny = 100; in.E = E; in.nu = nu; in.syy = -10e6; in.jumpModel = 1;
        const double h = in.Lx / in.nx, th = 30 * M_PI / 180, xl = 5 - 12 * h, xr = 5 + 12 * h, yc = 5 + 0.37 * h;
        efem::Fracture F; F.x0 = xl; F.x1 = xr; F.y0 = yc - (xr - xl) / 2 * std::tan(th); F.y1 = yc + (xr - xl) / 2 * std::tan(th); F.friction = 0.3;
        in.fracs.push_back(F);
        efem::Result r;
        if (efem::solve(in, r)) {
            const double a = std::hypot(F.x1 - F.x0, F.y1 - F.y0) / 2, sig = 10e6;
            const double drive = sig * std::sin(th) * std::cos(th) - 0.3 * sig * std::cos(th) * std::cos(th);
            double S = 0; for (const auto& c : r.cut) S += std::fabs(c.wt) * c.len;
            std::printf("    F9 comparison, element-constant EFEM (SKON) on the same mesh: slip integral / analytic = %.3f (continuous jumps: see above)\n", S / (2 * M_PI * drive * a * a / Ep));
        }
    }
}


// C20 GUI smoke finding: on the HFGROW example the correction pass hit rates that trigger fracture growth and the forward run
// stopped (time step below minimum), which aborted the whole design. The designer must keep the best verified schedule.
static void testSeismicityRobust() {
    std::printf("== C20 seismicity designer robustness: strongly nonlinear case (growing hydraulic fracture)\n");
    std::string path;
    for (const char* d : {"examples/HFGROW_INJECTION.DATA", "../examples/HFGROW_INJECTION.DATA", "../../examples/HFGROW_INJECTION.DATA"})
        if (std::ifstream(d).good()) { path = d; break; }
    if (path.empty()) { std::printf("    (examples/HFGROW_INJECTION.DATA not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(path);
    seis::DesignSetup s; s.wells = {"P"}; seis::Patch p; p.i = 45; p.j = 15; p.k = 0; s.patches = {p}; s.totalVolume = 60;
    const seis::DesignOutput o = seis::designInjection(c, s);
    std::printf("    ok %d, linear-model error %.1f %%, simulated events constant %.4g, returned %.4g, %d runs\n", int(o.ok), 100 * o.linearityError,
                o.simEventsConst, o.simEventsOpt, o.runs);
    check(o.ok, "design completes although a correction run fails (keeps the best verified schedule)");
    check(o.ok && o.simEventsOpt <= o.simEventsConst * (1 + 1e-12), "returned schedule is never worse than constant rate in the verification runs");
}
