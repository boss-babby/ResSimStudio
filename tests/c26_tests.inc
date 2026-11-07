// =====================================================================================================================
// C26 tests: 3-D DDM (triangular dislocation elements), fault stress transfer + quasi-static slip in the run (FAULTSLIP),
// DDM apertures of in-run hydraulic fractures (HFGROW ... DDMWIDTH)
// =====================================================================================================================
#include "ressim/ddm.h"
#include "ressim/hf_ddm.h"

static void testC26Kernel() {
    std::printf("== C26 TDE kernel (Nikkhoo & Walter 2015): strain vs finite differences of the displacement, free surface\n");
    const ddm::Tri t{{-1.0, -0.6, -3.0}, {1.2, -0.2, -2.6}, {0.1, 0.9, -3.4}};
    const ddm::V3 slip = {0.7, -0.4, 0.25};
    const double nu = 0.25;
    for (int hs = 0; hs < 2; ++hs) {
        double worst = 0, scale = 0;
        const ddm::V3 pts[4] = {{0.3, 0.2, -1.5}, {2.5, -1.0, -3.2}, {-0.4, 0.1, -3.05}, {0.0, 0.0, -0.4}};
        for (const auto& o : pts) {
            double e[6]; ddm::triStrain(o, t, slip, nu, hs != 0, e);
            const double h = 1e-5;
            double J[3][3];
            for (int j = 0; j < 3; ++j) {
                ddm::V3 a = o, b = o; a[size_t(j)] += h; b[size_t(j)] -= h;
                double ua[3], ub[3]; ddm::triDisp(a, t, slip, nu, hs != 0, ua); ddm::triDisp(b, t, slip, nu, hs != 0, ub);
                for (int i = 0; i < 3; ++i) J[i][j] = (ua[i] - ub[i]) / (2 * h);
            }
            const double fd[6] = {J[0][0], J[1][1], J[2][2], 0.5 * (J[0][1] + J[1][0]), 0.5 * (J[0][2] + J[2][0]), 0.5 * (J[1][2] + J[2][1])};
            for (int r = 0; r < 6; ++r) { worst = std::max(worst, std::fabs(e[r] - fd[r])); scale = std::max(scale, std::fabs(fd[r])); }
        }
        std::printf("    %s: max |strain - FD| %.3g (scale %.3g)\n", hs ? "half space" : "full space", worst, scale);
        check(worst < 1e-6 * std::max(1.0, scale) + 1e-8, hs ? "TDE strain = symmetric gradient of the TDE displacement (half space)" : "TDE strain = symmetric gradient of the TDE displacement (full space)");
    }
    // free surface: traction sigma.e_z = 0 at z = 0 for the half-space solution
    double S[6]; ddm::triStress({0.4, -0.3, 0.0}, t, slip, 1.0, nu, true, S);
    double Sf[6]; ddm::triStress({0.4, -0.3, 0.0}, t, slip, 1.0, nu, false, Sf);
    std::printf("    traction at the free surface: half space (%.2g, %.2g, %.2g), full space (%.2g, %.2g, %.2g) per unit G\n", S[4], S[5], S[2], Sf[4], Sf[5], Sf[2]);
    check(std::fabs(S[2]) + std::fabs(S[4]) + std::fabs(S[5]) < 1e-9 * (std::fabs(Sf[2]) + std::fabs(Sf[4]) + std::fabs(Sf[5])) + 1e-12,
          "half-space TDE: traction-free surface z = 0");
}

// square elements on the plane z = z0 (normal +z), centres inside radius a (or all of an nx x ny patch)
static std::vector<ddm::Element> c26Plane(double h, int nx, int ny, double z0, double a = -1) {
    std::vector<ddm::Element> el;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const double x0 = (i - nx / 2.0) * h, y0 = (j - ny / 2.0) * h;
        const double cx = x0 + 0.5 * h, cy = y0 + 0.5 * h;
        if (a > 0 && cx * cx + cy * cy > a * a) continue;
        el.push_back(ddm::quadElement({x0, y0, z0}, {x0 + h, y0, z0}, {x0 + h, y0 + h, z0}, {x0, y0 + h, z0}));
    }
    return el;
}

static void testC26Penny() {
    std::printf("== C26 DDM opening: penny-shaped crack (Sneddon 1946) and a long rectangular crack (plane strain)\n");
    const double E = 10e9, nu = 0.25, G = E / (2 * (1 + nu)), p = 1e6, a = 1.0;
    const int nd = 24;
    const double h = 2 * a / nd;
    auto el = c26Plane(h, nd, nd, -1000.0, a);
    const int n = int(el.size());
    const auto t0 = std::chrono::steady_clock::now();
    const auto K = ddm::buildOpening(el, G, nu, false);
    std::vector<double> w; const auto R = ddm::solveOpening(K, n, std::vector<double>(size_t(n), p), w);
    double V = 0, A = 0, wc = 0, rc = 1e9;
    for (int i = 0; i < n; ++i) {
        V += w[size_t(i)] * el[size_t(i)].area; A += el[size_t(i)].area;
        const double r = std::hypot(el[size_t(i)].c[0], el[size_t(i)].c[1]);
        if (r < rc) { rc = r; wc = w[size_t(i)]; }
    }
    const double aEff = std::sqrt(A / M_PI);
    const double Vref = ddm::pennyVolume(p, aEff, E, nu), wref = ddm::pennyOpening(p, aEff, rc, E, nu);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("    %d elements (a_eff %.4f m), %d active-set its: volume %.5g vs %.5g m3 (%.2f %%), centre opening %.5g vs %.5g m (%.2f %%), %.1f s\n",
                n, aEff, R.iterations, V, Vref, 100 * (V / Vref - 1), wc, wref, 100 * (wc / wref - 1), sec);
    check(R.converged && R.slipping == n, "pressurised penny crack: all elements open");
    check(std::fabs(V / Vref - 1) < 0.06, "penny crack volume within 6 % of Sneddon (constant-DD elements, staircase edge)");
    check(std::fabs(wc / wref - 1) < 0.05, "penny crack centre opening within 5 % of Sneddon");
    // long rectangle L = 8 H, height H in z: centre opening -> plane-strain 2 p H (1 - nu^2) / E (PKN cross-section)
    const double H = 10.0, L = 80.0; const int nh = 8, nl = 32;
    std::vector<ddm::Element> R2;
    for (int j = 0; j < nh; ++j) for (int i = 0; i < nl; ++i) {
        const double x0 = -L / 2 + i * L / nl, z0 = -2000.0 - H / 2 + j * H / nh;
        R2.push_back(ddm::quadElement({x0, 0, z0}, {x0 + L / nl, 0, z0}, {x0 + L / nl, 0, z0 + H / nh}, {x0, 0, z0 + H / nh}));
    }
    const auto K2 = ddm::buildOpening(R2, G, nu, true);
    std::vector<double> w2; ddm::solveOpening(K2, int(R2.size()), std::vector<double>(R2.size(), p), w2);
    double wmid = 0; for (int j = nh / 2 - 1; j <= nh / 2; ++j) wmid += 0.5 * w2[size_t(j * nl + nl / 2)];
    const double wps = 2 * p * H * (1 - nu * nu) / E;
    std::printf("    rectangle %g x %g m (%d x %d), centre opening %.5g m vs plane strain %.5g m (%.2f %%)\n", L, H, nl, nh, wmid, wps, 100 * (wmid / wps - 1));
    check(std::fabs(wmid / wps - 1) < 0.08, "long rectangular crack: centre opening within 8 % of the plane-strain (PKN) value");
}

static void testC26Shear() {
    std::printf("== C26 DDM shear interaction: reciprocity, loading of the neighbours, frictional equilibrium after slip\n");
    const double G = 10e9, nu = 0.25, h = 10.0;
    const int nx = 12, ny = 12;
    auto el = c26Plane(h, nx, ny, -1500.0);
    const int n = int(el.size());
    std::vector<ddm::V3> dir(size_t(n), ddm::V3{1, 0, 0});
    const auto K = ddm::buildShear(el, dir, G, nu, true, 0.0);
    auto kt = [&](int i, int j) { for (int q = K.start[size_t(i)]; q < K.start[size_t(i) + 1]; ++q) if (K.col[size_t(q)] == j) return K.kt[size_t(q)]; return 0.0; };
    const int c = (ny / 2) * nx + nx / 2, east = c + 1, north = c + nx;
    const double kself = kt(c, c), kE = kt(east, c), kN = kt(north, c);
    double asym = 0; for (int i = 0; i < n; i += 7) for (int j = 0; j < n; j += 5) asym = std::max(asym, std::fabs(kt(i, j) - kt(j, i)) / std::fabs(kself));
    // single square element in full space, uniform slip: centre stress drop ~ G s / h x 0.66 (numerical constant of a square
    // dislocation loop; reported for information)
    std::printf("    self K_tau %.4g Pa/m (x h / G = %.4f), mode-II neighbour %+.4g, mode-III neighbour %+.4g, max |K_ij - K_ji| / |K_ii| %.2g, %ld pairs, %.2f s\n",
                kself, kself * h / G, kE, kN, asym, K.pairs, K.seconds);
    check(kself < 0, "own slip relieves own shear stress (K_tau(i,i) < 0)");
    check(kE > 0 && kN > 0, "slip loads the in-plane neighbours (mode II and mode III directions)");
    check(asym < 1e-6, "Betti reciprocity: K_tau symmetric for equal coplanar elements");
    // LCP: central 4 x 4 over-stressed by 1 MPa, friction 0.6, no normal coupling (coplanar) -> equilibrium
    std::vector<double> ex(size_t(n), -0.3e6), mu(size_t(n), 0.6), ds;
    for (int j = ny / 2 - 2; j < ny / 2 + 2; ++j) for (int i = nx / 2 - 2; i < nx / 2 + 2; ++i) ex[size_t(j * nx + i)] = 1e6;
    const auto R = ddm::solveSlip(K, ex, mu, ds, 10.0, 20000);
    double fmaxSlip = 0, fmaxLock = -1e300, smax = 0;
    for (int i = 0; i < n; ++i) {
        double f = ex[size_t(i)];
        for (int q = K.start[size_t(i)]; q < K.start[size_t(i) + 1]; ++q) f += (K.kt[size_t(q)] - mu[size_t(i)] * K.kn[size_t(q)]) * ds[size_t(K.col[size_t(q)])];
        if (ds[size_t(i)] > 0) fmaxSlip = std::max(fmaxSlip, std::fabs(f)); else fmaxLock = std::max(fmaxLock, f);
        smax = std::max(smax, ds[size_t(i)]);
    }
    std::printf("    patch 4 x 4 at +1 MPa: %d slipping, %d its, max slip %.4g mm, |excess| on slipping %.3g Pa, max excess on locked %.3g Pa\n",
                R.slipping, R.iterations, smax * 1e3, fmaxSlip, fmaxLock);
    check(R.converged && R.slipping >= 16, "quasi-static slip converges (projected Gauss-Seidel)");
    check(fmaxSlip < 20 && fmaxLock < 20, "frictional equilibrium: slipping faces at strength, locked faces below");
}

static void testC26FaultRun() {
    std::printf("== C26 fault stress transfer in the run (FAULT_INJECTION.DATA + FAULTSLIP): slip, moment, loading of neighbours\n");
    const std::string path = findExample("FAULT_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/FAULT_INJECTION.DATA not found - skipped)\n"); return; }
    struct R { double maxSlip = 0, moment = 0, momentChk = 0, expected = 0, maxDcfsLocked = -1e300, eqRes = 0, units = 1, maxLoadUnslipped = 0; int slipped = 0, ev2 = 0, crit0 = 0; long pairs = 0; double sec = 0; };
    auto run = [&](int mode) {   // 0 no slip, 1 slip (no weakening), 2 slip with weakening dMu 0.1
        R r;
        Case c = importEclipseDeck(path);
        c.faultMech.slipOn = mode > 0; c.faultMech.slipE = 20; c.faultMech.slipNu = 0.25; c.faultMech.slipDMu = mode == 2 ? 0.1 : 0.0;
        for (auto& f : c.faultMech.faults) f.rsfLayer = 0;   // RSF trace not needed here
        const auto t0 = std::chrono::steady_clock::now();
        Simulator sim; sim.setup(c); sim.initialize();
        faultmech::Monitor M; std::string err;
        if (!M.setup(c, sim.grid(), err)) { std::printf("    %s\n", err.c_str()); return r; }
        M.step(sim.time(), sim.snapshots().back().pressure);
        while (sim.runReportStep()) M.step(sim.time(), sim.snapshots().back().pressure);
        const auto& O = M.out();
        const auto& S = O.steps.back();
        std::vector<char> ever(O.faces.size(), 0);
        for (const auto& st : O.steps) for (size_t i = 0; i < st.slip.size(); ++i) if (st.slip[i] > 0) ever[i] = 1;
        for (size_t i = 0; i < O.faces.size(); ++i) {
            if (mode > 0 && S.slip[i] > 0) { ++r.slipped; r.maxSlip = std::max(r.maxSlip, double(S.slip[i])); r.momentChk += O.shearModulus * O.faces[i].area * S.slip[i]; }
            if (!ever[i]) r.maxDcfsLocked = std::max(r.maxDcfsLocked, double(S.dcfs[i]));
        }
        // equilibrium of the faces that slipped in the last slipping step: CFF - max(0, CFF0) ~ 0 (bar)
        for (size_t s = 1; s < O.steps.size(); ++s) {
            const auto& A = O.steps[s]; const auto& B = O.steps[s - 1];
            if (A.slip.empty()) continue;
            for (size_t i = 0; i < O.faces.size(); ++i)
                if (A.slip[i] > B.slip[i] && mode == 1) r.eqRes = std::max(r.eqRes, std::fabs(double(A.cff[i]) - std::max(0.0, O.faces[i].cff0)));
        }
        // stress transfer: faces not (yet) slipped in a step that carry a positive shear-stress change from slip elsewhere
        for (const auto& st : O.steps) for (size_t i = 0; i < st.slip.size() && i < st.dtauSlip.size(); ++i)
            if (st.slip[i] <= 0) r.maxLoadUnslipped = std::max(r.maxLoadUnslipped, double(st.dtauSlip[i]));
        r.moment = S.moment; r.expected = O.series[0].expected.back(); r.pairs = O.ddmPairs; r.crit0 = O.criticalAtStart; r.units = O.mPerUnit * O.mPerUnit;
        for (const auto& e : O.events) if (e.kind == 2) ++r.ev2;
        r.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return r;
    };
    const R a = run(0), b = run(1), w = run(2);
    std::printf("    faces critical at t0: %d\n", a.crit0);
    std::printf("    no slip   : expected events %.4g, max dCFS of never-slipping faces %.3f bar (%.1f s)\n", a.expected, a.maxDcfsLocked, a.sec);
    std::printf("    slip      : %d faces slipped, max slip %.4g mm, M0 %.4g N m (check %.4g), Mw-equivalent %.2f, expected events %.4g,\n"
                "                max dCFS of never-slipping faces %.3f bar, equilibrium residual %.3g bar, %ld DDM pairs (%.1f s)\n",
                b.slipped, b.maxSlip * 1e3, b.moment, b.momentChk * b.units, b.moment > 0 ? 2.0 / 3.0 * (std::log10(b.moment) - 9.1) : 0.0,
                b.expected, b.maxDcfsLocked > -1e30 ? b.maxDcfsLocked : 0.0 /* none: every face slipped */, b.eqRes, b.pairs, b.sec);
    std::printf("    weakening : dMu 0.1 -> %d faces, max slip %.4g mm, M0 %.4g N m, %d seismic slip events (%.1f s)\n", w.slipped, w.maxSlip * 1e3, w.moment, w.ev2, w.sec);
    check(b.slipped > 0 && b.maxSlip > 0, "faces at the Coulomb limit slip (quasi-static DDM)");
    check(std::fabs(b.moment - b.momentChk * b.units) <= 1e-6 * std::max(1.0, b.moment), "seismic moment = G sum(A slip)");
    check(b.eqRes < 0.5, "slipping faces end at their strength (CFF = max(0, CFF0) within 0.5 bar)");
    std::printf("    largest shear-stress gain from slip on a face that had not slipped yet: %.4g bar\n", b.maxLoadUnslipped);
    check(b.maxLoadUnslipped > 1e-3, "slip loads the faces that have not slipped (stress transfer, positive d tau)");
    check(w.moment > b.moment && w.ev2 > 0, "dynamic weakening: larger moment and seismic slip events (cascade)");
}

static void testC26HfWidth() {
    std::printf("== C26 DDM apertures of an in-run hydraulic fracture (HFGROW_INJECTION.DATA + DDMWIDTH)\n");
    const std::string path = findExample("HFGROW_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/HFGROW_INJECTION.DATA not found - skipped)\n"); return; }
    struct R { double fpr = 0, pmax = 0, fwit = 0, maxW = 0, vol = 0, maxMult = 0; int open = 0, el = 0; size_t ev = 0; double sec = 0; };
    auto run = [&](bool ddmw) {
        R r;
        Case c = importEclipseDeck(path);
        c.grid.hfracs[0].ddmWidth = ddmw ? 1 : 0; c.grid.hfracs[0].ddmE = 20; c.grid.hfracs[0].ddmNu = 0.25;
        const auto t0 = std::chrono::steady_clock::now();
        Simulator sim; sim.setup(c); sim.initialize();
        hfddm::Aperture A; std::string err;
        const bool ok = ddmw && A.setup(sim, err);
        if (ddmw && !ok) std::printf("    %s\n", err.c_str());
        if (ok) A.update(sim);
        while (sim.runReportStep()) if (ok) A.update(sim);
        r.fwit = sim.summary().back().fwit; r.fpr = sim.summary().back().fpr; for (double v : sim.snapshots().back().pressure) r.pmax = std::max(r.pmax, v); r.ev = sim.hfGrowthEvents().size();
        if (ok) { r.maxW = A.out().maxWidth; r.vol = A.out().volume; r.open = A.out().open; r.el = A.out().elements; for (double m : A.out().mult) r.maxMult = std::max(r.maxMult, m); }
        r.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return r;
    };
    const R a = run(false), b = run(true);
    std::printf("    fixed width : FWIT %.6g, FPR %.6g, max p %.6g, %zu growth events (%.1f s)\n", a.fwit, a.fpr, a.pmax, a.ev, a.sec);
    std::printf("    DDM width   : FWIT %.6g, FPR %.6g, max p %.6g, %zu growth events, %d / %d pieces open, max aperture %.3f mm, opening volume %.4g m3, max conductance x%.3g (%.1f s)\n",
                b.fwit, b.fpr, b.pmax, b.ev, b.open, b.el, b.maxW * 1e3, b.vol, b.maxMult, b.sec);
    check(b.el > 0 && b.open > 0 && b.maxW > 0, "DDM apertures computed for the open fracture pieces");
    check(b.maxW < 0.05, "apertures physically bounded (< 5 cm)");
    check(std::fabs(b.pmax - a.pmax) > 1e-6 * std::max(1.0, a.pmax) || b.ev != a.ev, "mechanical aperture changes the fracture conductance and the run");
}

static void testC26RoundTrip() {
    std::printf("== C26 deck round trip: FAULTSLIP, HFGROW ... DDMWIDTH\n");
    const std::string path = findExample("FAULT_INJECTION.DATA");
    const std::string p2 = findExample("HFGROW_INJECTION.DATA");
    if (path.empty() || p2.empty()) { std::printf("    (examples not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    c.faultMech.slipOn = true; c.faultMech.slipE = 15; c.faultMech.slipNu = 0.22; c.faultMech.slipDMu = 0.05; c.faultMech.slipRefMm = 3; c.faultMech.slipCutoff = 30; c.faultMech.slipHalfSpace = 0;
    const Case d = projectFromText(projectText(c));
    check(d.faultMech.slipOn && std::fabs(d.faultMech.slipE - 15) < 1e-9 && std::fabs(d.faultMech.slipNu - 0.22) < 1e-9 && std::fabs(d.faultMech.slipDMu - 0.05) < 1e-9 &&
          std::fabs(d.faultMech.slipRefMm - 3) < 1e-9 && std::fabs(d.faultMech.slipCutoff - 30) < 1e-9 && d.faultMech.slipHalfSpace == 0, "FAULTSLIP written and read back");
    Case h = importEclipseDeck(p2);
    h.grid.hfracs[0].ddmWidth = 1; h.grid.hfracs[0].ddmE = 25; h.grid.hfracs[0].ddmNu = 0.2;
    const Case h2 = projectFromText(projectText(h));
    check(!h2.grid.hfracs.empty() && h2.grid.hfracs[0].ddmWidth == 1 && std::fabs(h2.grid.hfracs[0].ddmE - 25) < 1e-9 && std::fabs(h2.grid.hfracs[0].ddmNu - 0.2) < 1e-9,
          "HFGROW ... DDMWIDTH E nu written and read back");
}

// ---------------------------------------------------------------- C26 (C25-4a): iterative two-way FE coupling
// GEOMECH coupling 2 repeats every report step until the pore-volume multipliers agree with the mechanics. Checks: the
// fixed point converges (passes > 1 on the loading steps, final change below tol), the history of the rewound simulator
// stays consistent (one summary point per time step, increasing times), and the coarse-step iterative answer lies closer
// to a reference with 5 d report steps (same flow time steps, <= 5 d, in every run) than the coarse explicit answer.
static void testC26IterCoupling() {
    std::printf("== C26 iterative two-way FE geomechanics (fixed point per report step) vs explicit\n");
    const std::string path = findExample("GEOMECH_CP.DATA");
    if (path.empty()) { std::printf("    (examples/GEOMECH_CP.DATA not found - skipped)\n"); return; }
    struct R { double pAvg = 0, upl = 0, fwit = 0; int passes = 0, steps = 0, notConv = 0; double lastDiff = 0; bool mono = true; size_t nSum = 0; double sec = 0; };
    auto run = [&](int coupling, int split) {
        R r;
        Case c = importEclipseDeck(path);
        c.grid.hfracs.clear();                      // no fracture growth: isolate the flow <-> mechanics coupling
        c.faultMech.svGrad = 0;
        c.geoSeq.E = {1.5e9};                       // soft rock: strong pore-volume coupling
        c.geoSeq.coupling = coupling; c.geoSeq.maxIter = 10; c.geoSeq.tol = 1e-5;
        {   // uniform report steps of 'split' days over the 450 d of the deck; flow time steps <= 5 d in every run, so the runs
            // differ only in how often (explicit) / how consistently (iterative) the mechanics is coupled
            const ReportStep first = c.schedule.front();
            double total = 0; for (const auto& q : c.schedule) total += q.length;
            std::vector<ReportStep> s2;
            for (double t = 0; t < total - 1e-9; t += split) { ReportStep q = s2.empty() ? first : ReportStep{}; q.length = std::min(double(split), total - t); s2.push_back(q); }
            c.schedule = s2;
            c.tuning.maxDt = 5; c.tuning.initialDt = 1;
        }
        const auto t0 = std::chrono::steady_clock::now();
        Simulator sim; sim.log = [](const std::string&) {}; sim.setup(c); sim.initialize();
        const Grid& G = sim.grid();
        const double paU = 1e5;
        geoseq::Coupler C; std::string err;
        std::vector<double> p(size_t(G.numActive())), phi(size_t(G.numActive()), 0.18);
        const auto& P0 = sim.snapshots().back().pressure; for (size_t a = 0; a < p.size(); ++a) p[a] = double(P0[a]) * paU;
        if (!C.setup(c.geoSeq, G, G.globalOf, 1.0, p, phi, err)) { std::printf("    setup: %s\n", err.c_str()); return r; }
        C.step(p);
        std::vector<double> pv = C.pvFactor(p, c.props.rock.cr / paU); sim.setPoreVolumeMultipliers(pv);
        for (;;) {
            if (coupling == 2) {
                const geoseq::IterStep S = geoseq::iterativeReportStep(sim, C, pv, paU, c.props.rock.cr / paU, c.geoSeq.maxIter, c.geoSeq.tol);
                if (!S.ran) break;
                r.passes += S.passes; ++r.steps; if (!S.converged) ++r.notConv; r.lastDiff = S.diff;
            } else {
                if (!sim.runReportStep()) break;
                const auto& P = sim.snapshots().back().pressure; for (size_t a = 0; a < p.size(); ++a) p[a] = double(P[a]) * paU;
                C.step(p); pv = C.pvFactor(p, c.props.rock.cr / paU); sim.setPoreVolumeMultipliers(pv);
                ++r.steps; ++r.passes;
            }
        }
        const auto& SS = sim.summary();
        for (size_t i = 1; i < SS.size(); ++i) if (!(SS[i].time > SS[i - 1].time)) r.mono = false;
        r.nSum = SS.size(); r.fwit = SS.back().fwit;
        double v = 0, pp = 0; for (int a = 0; a < G.numActive(); ++a) { v += G.poreVolRef[size_t(a)]; pp += sim.snapshots().back().pressure[size_t(a)] * G.poreVolRef[size_t(a)]; }
        r.pAvg = pp / v; r.upl = C.solver().surfaceUplift();
        r.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return r;
    };
    const R E = run(1, 150), I = run(2, 150), F = run(2, 5), EF = run(1, 5);
    std::printf("    explicit, 150 d report steps: p_avg %.4f bar  uplift %.4f mm  FWIT %.6g  %d steps  %.1f s\n", E.pAvg, E.upl * 1e3, E.fwit, E.steps, E.sec);
    std::printf("    iterative, 150 d report steps: p_avg %.4f bar  uplift %.4f mm  FWIT %.6g  %d steps, %.2f passes/step, %d at maxIter, last change %.2e  %.1f s\n",
                I.pAvg, I.upl * 1e3, I.fwit, I.steps, double(I.passes) / std::max(1, I.steps), I.notConv, I.lastDiff, I.sec);
    std::printf("    iterative, 5 d report steps (reference): p_avg %.4f bar  uplift %.4f mm  FWIT %.6g  %d steps  %.1f s\n", F.pAvg, F.upl * 1e3, F.fwit, F.steps, F.sec);
    std::printf("    explicit, 5 d report steps: p_avg %.4f bar  uplift %.4f mm  (|p_avg - reference| %.4f bar)\n", EF.pAvg, EF.upl * 1e3, std::fabs(EF.pAvg - F.pAvg));
    const double eE = std::fabs(E.pAvg - F.pAvg), eI = std::fabs(I.pAvg - F.pAvg);
    std::printf("    |p_avg - fine|: explicit %.4f, iterative %.4f bar\n", eE, eI);
    check(I.steps == E.steps && I.steps > 0, "C26 iterative coupling: same report steps as explicit");
    check(I.passes > I.steps, "C26 iterative coupling: report steps repeated (fixed-point passes > 1)");
    check(I.notConv == 0, "C26 iterative coupling: fixed point converged on every step (tol 1e-5, 10 passes)");
    check(I.mono && I.nSum > size_t(I.steps), "C26 iterative coupling: rewound history consistent (one point per accepted time step, increasing times)");
    check(eI <= eE + 1e-9, "C26 iterative coupling: closer to the 5 d reference than explicit at 150 d report steps");
}

// ---------------------------------------------------------------- C26 (C24-b): phase field - volume control, orientation, in-situ stress
static void testC26PhaseField() {
    std::printf("== C26 phase field: volume-controlled (stable) growth, inclined crack, anisotropic in-situ stress\n");
    auto go = [](double ang, double sxx, double syy, double sxy, int steps, double vFac, const char* tag) {
        phasefield::Input in; in.h = 0.01; in.control = 1; in.steps = steps; in.angleDeg = ang; in.sxx = sxx; in.syy = syy; in.sxy = sxy;
        {   // injected volume: vFac x the Sneddon volume at the Griffith net pressure
            const double Ep = in.E / (1 - in.nu * in.nu), pg = std::sqrt(Ep * in.Gc / (M_PI * in.crackA));
            in.vMax = vFac * 2 * M_PI * pg * in.crackA * in.crackA / Ep;
        }
        phasefield::Output o; phasefield::run(in, o);
        const double K = std::sqrt(in.Gc * o.Eprime), pG = std::sqrt(o.Eprime * in.Gc / (M_PI * in.crackA));
        std::printf("    %-34s %d steps, %d nodes, %.0f s: sigma_n %.2f tau %.2f MPa, peak p %.3f MPa (theory %.3f, net ratio %.3f), tip direction %.1f deg\n",
                    tag, int(o.steps.size()), o.nodes, o.seconds, o.sigmaN * 1e-6, o.tauS * 1e-6, o.pcNumeric * 1e-6, o.pcTheory * 1e-6,
                    (o.pcNumeric - o.sigmaN) / std::max(1e-9, o.pcTheory - o.sigmaN), o.tipAngleDeg);
        (void)K; (void)pG;
        return o;
    };
    const phasefield::Output A = go(0, 0, 0, 0, 14, 3.5, "horizontal, no stress:");
    // stable growth: after the peak the pressure falls while the crack grows; toughness regime p sqrt(pi a) -> K_IC
    bool ok = A.ok && A.steps.size() > 5, grows = false, falls = true; size_t iPk = 0;
    for (size_t i = 0; i < A.steps.size(); ++i) if (A.steps[i].p > A.steps[iPk].p) iPk = i;
    for (size_t i = iPk + 1; i < A.steps.size(); ++i) { if (A.steps[i].p > A.steps[i - 1].p * 1.02) falls = false; }
    const auto& last = A.steps.empty() ? phasefield::Step{} : A.steps.back();
    grows = last.lengthHalf > 1.5 * 0.1;
    const double K = std::sqrt(100.0 * A.Eprime), kr = last.p * std::sqrt(M_PI * last.lengthHalf) / K;
    for (size_t i = 0; i < A.steps.size(); i += std::max<size_t>(1, A.steps.size() / 8))
        std::printf("      V %.4g m2  p %.3f MPa  a %.4f m\n", A.steps[i].volume, A.steps[i].p * 1e-6, A.steps[i].lengthHalf);
    std::printf("      last step: a %.4f m, p sqrt(pi a) / K_IC = %.3f\n", last.lengthHalf, kr);
    const double rA = A.pcNumeric / A.pcTheory;
    check(ok && grows, "C26 phase field volume control: crack grows past 1.5 a (growth after the peak is followed)");
    check(falls && iPk + 1 < A.steps.size(), "C26 phase field volume control: pressure peaks and then falls (stable growth)");
    check(rA > 0.75 && rA < 1.05, "C26 phase field volume control: peak pressure within the known AT2 band of Griffith (l/a = 0.2)");
    check(kr > 0.7 && kr < 1.15, "C26 phase field volume control: toughness regime p sqrt(pi a) ~ K_IC during growth");
    const phasefield::Output B = go(30, 5e6, 5e6, 0, 12, 3.0, "30 deg, isotropic 5 MPa:");
    const double rB = (B.pcNumeric - B.sigmaN) / (B.pcTheory - B.sigmaN);
    check(B.ok && std::fabs(B.sigmaN - 5e6) < 1 && std::fabs(rB - rA) < 0.08 * rA, "C26 phase field: inclined crack under isotropic stress - net peak pressure as the horizontal case (within 8 %)");
    const phasefield::Output C = go(45, 6e6, 5e6, 0, 12, 3.0, "45 deg, Sxx 6 / Syy 5 MPa:");
    check(C.ok && std::fabs(C.sigmaN - 5.5e6) < 1 && std::fabs(std::fabs(C.tauS) - 0.5e6) < 1, "C26 phase field: normal / shear stress on the inclined crack (5.5 / 0.5 MPa)");
    check(C.ok && C.steps.size() > 5 && C.tipAngleDeg < 45 - 5, "C26 phase field: the crack turns towards the maximum compressive stress (x)");
}
